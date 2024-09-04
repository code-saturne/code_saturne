/*============================================================================
 * Methods for particle coupling
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

/*============================================================================
 * Functions dealing with Lagrangian coupling
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <limits.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_array.h"
#include "cs_base.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

#include "bft_mem.h"
#include "bft_error.h"

#include "cs_physical_constants.h"
#include "cs_time_step.h"
#include "cs_turbulence_model.h"

#include "cs_lagr.h"
#include "cs_lagr_particle.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_coupling.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Constants */

static const cs_real_t  _c_stephan = 5.6703e-8;

/* Tensor to vector (t2v) and vector to tensor (v2t) mask arrays */

static const cs_lnum_t _iv2t[6] = {0, 1, 2, 0, 1, 0};
static const cs_lnum_t _jv2t[6] = {0, 1, 2, 1, 2, 2};

static const cs_lnum_t _t2v[3][3] = {{0, 3, 5},
                                     {3, 1, 4},
                                     {5, 4, 2}};

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute source terms for Lagrangian 2-way coupling.
 *
 * \remark  Source terms are computed for the starting cell of a particle
 *          during a given iteration. Even if particle exits the domain,
 *          it s necessary to compute a source term matching the exchange
 *          between the carrier fluid and the particle at the beginning
 *          of the time step. If cs_glob_lagr_time_step->nor == 2 and the
 *          particle interacts with a boundary, then the source terms
 *          are computed as if nor == 1.
 *
 * \param[in]   taup    dynamic characteristic time
 * \param[in]   tempct  thermal characteristic time
 * \param[out]  tsfext  external forces
 * \param[in]   force_p forces per mass unit on particles (m/s^2)
 * \param[in]   cpgd1   devolatization term 1 for heterogeneous coal
 * \param[in]   cpgd2   devolatization term 2 for heterogeneous coal
 * \param[in]   cpght   combustion term for heterogeneous coal
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_coupling(const cs_real_t    taup[],
                 const cs_real_t    tempct[],
                 cs_real_t          tsfext[],
                 const cs_real_3_t *force_p,
                 const cs_real_t    cpgd1[],
                 const cs_real_t    cpgd2[],
                 const cs_real_t    cpght[])
{
  /*Note: t_* stands for temporary array, used in case of time moments */
  cs_real_t *st_p = NULL, *t_st_p = NULL;
  cs_real_3_t *st_vel = NULL, *t_st_vel = NULL;
  cs_real_t *st_imp_vel = NULL, *t_st_imp_vel = NULL;
  cs_real_6_t *st_rij = NULL, *t_st_rij = NULL;
  cs_real_t *st_k = NULL, *t_st_k = NULL;
  cs_real_t *st_t_e = NULL, *t_st_t_e = NULL;
  cs_real_t *st_t_i = NULL, *t_st_t_i = NULL;

  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  cs_real_t *cell_f_vol = mq->cell_f_vol;

  const int *restrict c_disable_flag = mq->c_disable_flag;
  cs_lnum_t has_dc = mq->has_disable_flag; /* Has cells disabled? */

  /* Initialization
     ============== */

  {
    cs_field_t *f = cs_field_by_name_try("lagr_st_pressure");
    if (f != NULL)
      st_p = f->val;
  }

  {
    cs_field_t *f = cs_field_by_name_try("lagr_st_velocity");
    if (f != NULL)
      st_vel = (cs_real_3_t *)(f->val);
  }

  {
    cs_field_t *f = cs_field_by_name_try("lagr_st_imp_velocity");
    if (f != NULL)
      st_imp_vel = f->val;
  }

  {
    cs_field_t *f = cs_field_by_name_try("lagr_st_rij");
    if (f != NULL)
      st_rij = (cs_real_6_t *)(f->val);
  }

  {
    cs_field_t *f = cs_field_by_name_try("lagr_st_k");
    if (f != NULL)
      st_k = f->val;
  }

  {
    cs_field_t *f = cs_field_by_name_try("lagr_st_temperature");
    if (f != NULL)
      st_t_e = f->val;
  }

  {
    cs_field_t *f = cs_field_by_name_try("lagr_st_imp_temperature");
    if (f != NULL)
      st_t_i = f->val;
  }

  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;
  cs_lagr_source_terms_t *lag_st = cs_glob_lagr_source_terms;

  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;

  cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  cs_lnum_t ncel = cs_glob_mesh->n_cells;
  cs_lnum_t nbpart = p_set->n_particles;

  cs_real_t dtp = cs_glob_lagr_time_step->dtp;

  cs_real_3_t *auxl;
  BFT_MALLOC(auxl, nbpart, cs_real_3_t);

  /* Number of passes for steady source terms */
  if (   cs_glob_lagr_time_scheme->isttio == 1
      && cs_glob_time_step->nt_cur >= lag_st->nstits)
    lag_st->npts += 1;

  bool is_time_averaged = (   cs_glob_lagr_time_scheme->isttio == 1
                           && lag_st->npts > 0);

  lag_st->ntxerr = 0;
  lag_st->vmax = 0.0;
  lag_st->tmamax = 0.0;

  cs_real_t *volp = NULL, *volm = NULL;
  BFT_MALLOC(volp, ncel, cs_real_t);
  BFT_MALLOC(volm, ncel, cs_real_t);
  for (cs_lnum_t c_id = 0; c_id < ncel; c_id++) {
    volp[c_id] = 0.0;
    volm[c_id] = 0.0;
  }

  /* Preliminary computations
     ======================== */

  /* Finalization of external forces (if the particle interacts with a
     domain boundary, revert to order 1). */

  for (cs_lnum_t p_id = 0; p_id < nbpart; p_id++) {

    cs_real_t aux1 = dtp / taup[p_id];
    cs_real_t p_mass = cs_lagr_particles_get_real(p_set, p_id, CS_LAGR_MASS);

    if (   cs_glob_lagr_time_scheme->t_order == 1
        || cs_lagr_particles_get_lnum(p_set, p_id, CS_LAGR_REBOUND_ID) == 0)
      tsfext[p_id] = (1.0 - exp(-aux1)) * p_mass * taup[p_id];

    else
      tsfext[p_id] +=  (1.0 - (1.0 - exp (-aux1)) / aux1) * taup[p_id]
                    * p_mass;

  }

  for (cs_lnum_t p_id = 0; p_id < nbpart; p_id++) {

    cs_real_t  p_stat_w = cs_lagr_particles_get_real(p_set, p_id,
                                                     CS_LAGR_STAT_WEIGHT);
    cs_real_t  p_mass   = cs_lagr_particles_get_real(p_set, p_id, CS_LAGR_MASS);
    cs_real_t *p_vel    =
      cs_lagr_particles_attr_get_ptr<cs_real_t>(p_set, p_id, CS_LAGR_VELOCITY);

    cs_real_t  prev_p_mass = cs_lagr_particles_get_real_n(p_set, p_id, 1,
                                                          CS_LAGR_MASS);
    cs_real_t *prev_p_vel  =
      cs_lagr_particles_attr_n_get_ptr<cs_real_t>(p_set, p_id, 1,
                                                  CS_LAGR_VELOCITY);

    //TODO tsfext should be computed elsewhere (in sde) and the mass of particle
    // may be the previous mass.
    for (cs_lnum_t i = 0; i < 3; i++)
      auxl[p_id][i] = p_stat_w * (p_mass * p_vel[i] - prev_p_mass * prev_p_vel[i]
                                 - force_p[p_id][i] * tsfext[p_id]) / dtp;

  }

  /* Momentum source terms
     ===================== */

  if (lag_st->ltsdyn == 1) {

    if (is_time_averaged) {
      BFT_MALLOC(t_st_vel, n_cells_ext, cs_real_3_t);
      BFT_MALLOC(t_st_imp_vel, n_cells_ext, cs_real_t);
    }
    else {
      t_st_vel = st_vel;
      t_st_imp_vel = st_imp_vel;
    }

    cs_array_real_fill_zero(3 * n_cells_ext, (cs_real_t *) t_st_vel);
    cs_array_real_fill_zero(n_cells_ext,  t_st_imp_vel);

    for (cs_lnum_t p_id = 0; p_id < nbpart; p_id++) {

      unsigned char *particle = p_set->p_buffer + p_am->extents * p_id;

      cs_real_t  p_stat_w = cs_lagr_particle_get_real(particle, p_am,
                                                      CS_LAGR_STAT_WEIGHT);

      cs_real_t  prev_p_diam = cs_lagr_particle_get_real_n(particle, p_am, 1,
                                                           CS_LAGR_DIAMETER);
      cs_real_t  prev_p_mass = cs_lagr_particle_get_real_n(particle, p_am, 1,
                                                           CS_LAGR_MASS);
      cs_real_t  p_mass = cs_lagr_particle_get_real(particle, p_am,
                                                    CS_LAGR_MASS);

      cs_lnum_t c_id = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_CELL_ID);

      /* Volume and mass of particles in cell */
      volp[c_id] += p_stat_w * cs_math_pi * pow(prev_p_diam, 3) / 6.0;
      volm[c_id] += p_stat_w * prev_p_mass;

      /* Momentum source term */
      cs_real_t dvol = 0.;
      if (has_dc * c_disable_flag[has_dc * c_id] == 0)
        dvol = 1. / cell_f_vol[c_id];

      for (cs_lnum_t i = 0; i < 3; i++)
        t_st_vel[c_id][i] -= dvol * auxl[p_id][i];

      t_st_imp_vel[c_id] -= 2.0 * dvol * p_stat_w * p_mass / taup[p_id];

    }

  /* Turbulence source terms
     ======================= */

    cs_real_3_t *vel = (cs_real_3_t *)extra->vel->val;

    if (extra->itytur == 2 || extra->itytur == 4 ||
        extra->itytur == 5 || extra->iturb == CS_TURB_K_OMEGA) {

      /* In v2f the Lagrangian STs only influence k and epsilon
         (difficult to write something for v2, which loses its meaning as
         "Rij component") */

      if (is_time_averaged)
        BFT_MALLOC(t_st_k, n_cells_ext, cs_real_t);
      else
        t_st_k = st_k;

      cs_array_real_fill_zero(n_cells_ext, t_st_k);

      for (cs_lnum_t p_id = 0; p_id < nbpart; p_id++) {

        unsigned char *particle = p_set->p_buffer + p_am->extents * p_id;

        cs_lnum_t  c_id         = cs_lagr_particle_get_lnum(particle, p_am,
                                                           CS_LAGR_CELL_ID);
        cs_real_t *prev_f_vel  =
          (cs_real_t *)cs_lagr_particle_attr_n(particle, p_am, 1,
                                               CS_LAGR_VELOCITY_SEEN);
        cs_real_t *f_vel       =
          (cs_real_t *)cs_lagr_particle_attr(particle, p_am,
                                             CS_LAGR_VELOCITY_SEEN);

        cs_real_3_t vel_s =
        { 0.5 * (prev_f_vel[0] + f_vel[0]),
          0.5 * (prev_f_vel[1] + f_vel[1]),
          0.5 * (prev_f_vel[2] + f_vel[2])};

        cs_real_t dvol = 0.;
        if (has_dc * c_disable_flag[has_dc * c_id] == 0)
          dvol = 1. / cell_f_vol[c_id];
        t_st_k[c_id] -= dvol * cs_math_3_dot_product(vel_s, auxl[p_id]);

      }

      for (cs_lnum_t c_id = 0; c_id < ncel; c_id++)
        t_st_k[c_id] -= cs_math_3_dot_product(vel[c_id], t_st_vel[c_id]);

    }
    else if (extra->itytur == 3) {

      if (is_time_averaged)
        BFT_MALLOC(t_st_rij, n_cells_ext, cs_real_6_t);
      else
        t_st_rij = st_rij;

      cs_array_real_fill_zero(n_cells_ext * 6, (cs_real_t *)t_st_rij);

      for (cs_lnum_t p_id = 0; p_id < nbpart; p_id++) {

        unsigned char *particle = p_set->p_buffer + p_am->extents * p_id;

        cs_lnum_t  c_id         = cs_lagr_particle_get_lnum(particle, p_am,
                                                           CS_LAGR_CELL_ID);

        cs_real_t *prev_f_vel  =
          (cs_real_t *)cs_lagr_particle_attr_n(particle, p_am, 1,
                                               CS_LAGR_VELOCITY_SEEN);
        cs_real_t *f_vel       =
          (cs_real_t *)cs_lagr_particle_attr(particle, p_am,
                                             CS_LAGR_VELOCITY_SEEN);

        cs_real_3_t vel_s =
        { 0.5 * (prev_f_vel[0] + f_vel[0]),
          0.5 * (prev_f_vel[1] + f_vel[1]),
          0.5 * (prev_f_vel[2] + f_vel[2])};

        cs_real_t dvol = 0.;
        if (has_dc * c_disable_flag[has_dc * c_id] == 0)
          dvol = 1. / cell_f_vol[c_id];

        for (cs_lnum_t ij = 0; ij < 6; ij++) {
          cs_lnum_t i = _iv2t[ij];
          cs_lnum_t j = _jv2t[ij];

          t_st_rij[c_id][ij] -= ( vel_s[i] * auxl[p_id][j]
                                + vel_s[j] * auxl[p_id][i])*dvol;
        }


      }
      for (cs_lnum_t c_id = 0; c_id < ncel; c_id++) {
        for (cs_lnum_t ij = 0; ij < 6; ij++) {
          cs_lnum_t i = _iv2t[ij];
          cs_lnum_t j = _jv2t[ij];

          t_st_rij[c_id][ij] -= ( vel[c_id][i] * t_st_vel[c_id][j]
                                + vel[c_id][j] * t_st_vel[c_id][i]);

        }
      }

    }

  }

  /* Mass source terms
     ================= */

  if (    lag_st->ltsmas == 1
      && (   cs_glob_lagr_specific_physics->impvar == 1
          || cs_glob_lagr_specific_physics->idpvar == 1
          || cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_CTWR )) {

    if (is_time_averaged)
      BFT_MALLOC(t_st_p, n_cells_ext, cs_real_t);
    else
      t_st_p = st_p;

    cs_array_real_fill_zero(n_cells_ext, t_st_p);

    for (cs_lnum_t p_id = 0; p_id < nbpart; p_id++) {

      unsigned char *particle = p_set->p_buffer + p_am->extents * p_id;

      cs_real_t  p_stat_w
        = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_STAT_WEIGHT);
      cs_real_t  prev_p_mass
        = cs_lagr_particle_get_real_n(particle, p_am, 1, CS_LAGR_MASS);
      cs_real_t  p_mass
        = cs_lagr_particle_get_real_n(particle, p_am, 0, CS_LAGR_MASS);

      /* Fluid mass source term > 0 -> add mass to fluid */
      cs_lnum_t c_id = cs_lagr_particle_get_lnum(particle, p_am,
                                                    CS_LAGR_CELL_ID);

      cs_real_t dvol = 0.;
      if (has_dc * c_disable_flag[has_dc * c_id] == 0)
        dvol = 1. / cell_f_vol[c_id];

      t_st_p[c_id] += - p_stat_w * (p_mass - prev_p_mass) / dtp * dvol;

    }

  }

  /* Thermal source terms
     ==================== */

  if (lag_st->ltsthe == 1) {

    if (   cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_HEAT
        && cs_glob_lagr_specific_physics->itpvar == 1) {

      if (is_time_averaged) {
        BFT_MALLOC(t_st_t_e, n_cells_ext, cs_real_t);
        BFT_MALLOC(t_st_t_i, n_cells_ext, cs_real_t);
      }
      else {
        t_st_t_e = st_t_e;
        t_st_t_i = st_t_i;
      }

      cs_array_real_fill_zero(n_cells_ext, t_st_t_e);
      cs_array_real_fill_zero(n_cells_ext, t_st_t_i);

      for (cs_lnum_t p_id = 0; p_id < nbpart; p_id++) {

        unsigned char *particle = p_set->p_buffer + p_am->extents * p_id;
        cs_lnum_t  c_id = cs_lagr_particle_get_lnum(particle, p_am,
                                                   CS_LAGR_CELL_ID);
        cs_real_t  p_mass = cs_lagr_particle_get_real_n(particle, p_am, 0,
                                                        CS_LAGR_MASS);
        cs_real_t  prev_p_mass = cs_lagr_particle_get_real_n(particle, p_am, 1,
                                                             CS_LAGR_MASS);
        cs_real_t  p_cp = cs_lagr_particle_get_real_n(particle, p_am, 0,
                                                      CS_LAGR_CP);
        cs_real_t  prev_p_cp = cs_lagr_particle_get_real_n(particle, p_am, 1,
                                                           CS_LAGR_CP);
        cs_real_t  p_tmp = cs_lagr_particle_get_real_n(particle, p_am, 0,
                                                       CS_LAGR_TEMPERATURE);
        cs_real_t  prev_p_tmp = cs_lagr_particle_get_real_n(particle, p_am, 1,
                                                            CS_LAGR_TEMPERATURE);
        cs_real_t  p_stat_w = cs_lagr_particle_get_real(particle, p_am,
                                                        CS_LAGR_STAT_WEIGHT);

        cs_real_t dvol = 0.;
        if (has_dc * c_disable_flag[has_dc * c_id] == 0)
          dvol = 1. / cell_f_vol[c_id];

        t_st_t_e[c_id] += - (p_mass * p_tmp * p_cp
                            - prev_p_mass * prev_p_tmp * prev_p_cp
                            ) / dtp * p_stat_w * dvol;
        t_st_t_i[c_id] += tempct[nbpart + p_id] * p_stat_w; //FIXME not homogeneous

      }
      if (extra->radiative_model > 0) {

        for (cs_lnum_t p_id = 0; p_id < nbpart; p_id++) {

          unsigned char *particle = p_set->p_buffer + p_am->extents * p_id;
          cs_lnum_t  c_id = cs_lagr_particle_get_lnum(particle, p_am,
                                                     CS_LAGR_CELL_ID);
          cs_real_t  p_diam = cs_lagr_particle_get_real_n(particle, p_am, 0,
                                                          CS_LAGR_DIAMETER);
          cs_real_t  p_eps = cs_lagr_particle_get_real_n(particle, p_am, 0,
                                                         CS_LAGR_EMISSIVITY);
          cs_real_t  p_tmp = cs_lagr_particle_get_real_n(particle, p_am, 0,
                                                         CS_LAGR_TEMPERATURE);
          cs_real_t  p_stat_w = cs_lagr_particle_get_real(particle, p_am,
                                                          CS_LAGR_STAT_WEIGHT);

          cs_real_t dvol = 0.;
          if (has_dc * c_disable_flag[has_dc * c_id] == 0)
            dvol = 1. / cell_f_vol[c_id];

          cs_real_t aux1 = cs_math_pi * p_diam * p_diam * p_eps * dvol
                          * (extra->rad_energy->val[c_id]
                             - 4.0 * _c_stephan * cs_math_pow4(p_tmp));

          t_st_t_e[c_id] += aux1 * p_stat_w;

        }

      }

    }
    else if (cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_COAL) {
      if (cs_glob_lagr_const_dim->nlayer > 1)
        bft_error(__FILE__, __LINE__, 0,
                  _("Thermal coupling not implemented in multi-layer case"));

      else {

        if (is_time_averaged) {
          BFT_MALLOC(t_st_t_e, n_cells_ext, cs_real_t);
          BFT_MALLOC(t_st_t_i, n_cells_ext, cs_real_t);
        }
        else {
          t_st_t_e = st_t_e;
          t_st_t_i = st_t_i;
        }

        cs_array_real_fill_zero(n_cells_ext, t_st_t_e);
        cs_array_real_fill_zero(n_cells_ext, t_st_t_i);

        for (cs_lnum_t p_id = 0; p_id < nbpart; p_id++) {

          unsigned char *particle = p_set->p_buffer + p_am->extents * p_id;

          cs_lnum_t  c_id = cs_lagr_particle_get_lnum(particle, p_am,
                                                     CS_LAGR_CELL_ID);

          cs_real_t  p_mass = cs_lagr_particle_get_real_n(particle, p_am, 0,
                                                          CS_LAGR_MASS);
          cs_real_t  p_tmp = cs_lagr_particle_get_real(particle, p_am,
                                                       CS_LAGR_TEMPERATURE);
          cs_real_t  p_cp = cs_lagr_particle_get_real_n(particle, p_am, 0,
                                                        CS_LAGR_CP);

          cs_real_t  prev_p_mass = cs_lagr_particle_get_real_n
                                     (particle, p_am, 1, CS_LAGR_MASS);
          cs_real_t  prev_p_tmp  = cs_lagr_particle_get_real_n
                                     (particle, p_am, 1, CS_LAGR_TEMPERATURE);
          cs_real_t  prev_p_cp   = cs_lagr_particle_get_real_n
                                     (particle, p_am, 1, CS_LAGR_CP);

          cs_real_t  p_stat_w = cs_lagr_particle_get_real
                                  (particle, p_am, CS_LAGR_STAT_WEIGHT);

          cs_real_t dvol = 0.;
          if (has_dc * c_disable_flag[has_dc * c_id] == 0)
            dvol = 1. / cell_f_vol[c_id];

          t_st_t_e[c_id] += - (  p_mass * p_tmp * p_cp
                              - prev_p_mass * prev_p_tmp * prev_p_cp
                              ) / dtp * p_stat_w * dvol;
          t_st_t_i[c_id] += p_stat_w * p_mass * p_cp * dvol
                          / tempct[nbpart + p_id];

        }

      }

    }
    else if (cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_CTWR) {

      if (is_time_averaged) {
        BFT_MALLOC(t_st_t_e, n_cells_ext, cs_real_t);
        BFT_MALLOC(t_st_t_i, n_cells_ext, cs_real_t);
      }
      else {
        t_st_t_e = st_t_e;
        t_st_t_i = st_t_i;
      }

      cs_array_real_fill_zero(n_cells_ext, t_st_t_e);
      cs_array_real_fill_zero(n_cells_ext, t_st_t_i);

      for (cs_lnum_t p_id = 0; p_id < nbpart; p_id++) {

        unsigned char *particle = p_set->p_buffer + p_am->extents * p_id;

        cs_lnum_t  c_id = cs_lagr_particle_get_lnum(particle, p_am,
                                                   CS_LAGR_CELL_ID);

        cs_real_t  p_mass = cs_lagr_particle_get_real_n(particle, p_am, 0,
                                                        CS_LAGR_MASS);
        cs_real_t  p_tmp = cs_lagr_particle_get_real(particle, p_am,
                                                     CS_LAGR_TEMPERATURE);
        cs_real_t  p_cp = cs_lagr_particle_get_real_n(particle, p_am, 0,
                                                      CS_LAGR_CP);

        cs_real_t  prev_p_mass = cs_lagr_particle_get_real_n
                                   (particle, p_am, 1, CS_LAGR_MASS);
        cs_real_t  prev_p_tmp  = cs_lagr_particle_get_real_n
                                   (particle, p_am, 1, CS_LAGR_TEMPERATURE);
        cs_real_t  prev_p_cp   = cs_lagr_particle_get_real_n
                                   (particle, p_am, 1, CS_LAGR_CP);

        cs_real_t  p_stat_w = cs_lagr_particle_get_real
                                (particle, p_am, CS_LAGR_STAT_WEIGHT);

        cs_real_t dvol = 0.;
        if (has_dc * c_disable_flag[has_dc * c_id] == 0)
          dvol = 1. / cell_f_vol[c_id];

        t_st_t_e[c_id] += - (  p_mass * p_tmp * p_cp
                             - prev_p_mass * prev_p_tmp * prev_p_cp
                             ) / dtp * p_stat_w * dvol;
        t_st_t_i[c_id] += p_stat_w * p_mass * p_cp * dvol
                          / tempct[nbpart + p_id];

      }

    }

  }

  /* Check that the maximum admissible volume fraction of particles
     is not exceeded in some cells.
     ============================================================== */


  const cs_real_t tvmax = 0.8;
  const cs_real_t *cell_vol = cs_glob_mesh_quantities->cell_vol;

  for (cs_lnum_t c_id = 0; c_id < ncel; c_id++) {

    cs_real_t mf   = cell_vol[c_id] * extra->cromf->val[c_id];
    cs_real_t tauv = volp[c_id] / cell_vol[c_id];
    cs_real_t taum = volm[c_id] / mf;

    if (tauv > tvmax) {
      lag_st->ntxerr++;

      /* Note: it was not temporary array but directly st_val.
       * By consistency with st_vel, we put temporary array. */
      if (t_st_p != NULL)
        t_st_p[c_id] = 0.0;

      if (t_st_vel != NULL) {
        for (cs_lnum_t j = 0; j < 3; j++)
          t_st_vel[c_id][j] = 0.0;
      }

      if (t_st_imp_vel != NULL)
        t_st_imp_vel[c_id] = 0.0;

      if (t_st_k != NULL)
        t_st_k[c_id] = 0.0;

      if (t_st_rij != NULL) {
        for (cs_lnum_t j = 0; j < 6; j++)
          t_st_rij[c_id][j] = 0.0;
      }

      if (t_st_t_e != NULL)
        t_st_t_e[c_id] = 0.0;

      if (t_st_t_i != NULL)
        t_st_t_i[c_id] = 0.0;

    }

    lag_st->vmax = CS_MAX(tauv, lag_st->vmax);
    lag_st->tmamax = CS_MAX(taum, lag_st->tmamax);

  }

  /* Time average of source terms
     ============================ */

  if (is_time_averaged) {

    if (st_p != NULL) {
      for (cs_lnum_t c_id = 0; c_id < ncel; c_id++) {
        st_p[c_id]
          = (t_st_p[c_id] + (lag_st->npts - 1.0) * st_p[c_id])
          / lag_st->npts;
      }
    }

    if (st_vel != NULL) {
      for (cs_lnum_t c_id = 0; c_id < ncel; c_id++) {
        for (cs_lnum_t j = 0; j < 3; j++) {
          st_vel[c_id][j]
            =    (t_st_vel[c_id][j] + (lag_st->npts - 1.0) * st_vel[c_id][j])
               / lag_st->npts;
        }
      }
    }

    if (st_imp_vel != NULL) {
      for (cs_lnum_t c_id = 0; c_id < ncel; c_id++) {
        st_imp_vel[c_id]
          = (t_st_imp_vel[c_id] + (lag_st->npts - 1.0) * st_imp_vel[c_id])
          / lag_st->npts;
      }
    }

    if (st_k != NULL) {
      for (cs_lnum_t c_id = 0; c_id < ncel; c_id++) {
        st_k[c_id]
          = (t_st_k[c_id] + (lag_st->npts - 1.0) * st_k[c_id])
          / lag_st->npts;
      }
    }

    if (st_rij != NULL) {
      for (cs_lnum_t c_id = 0; c_id < ncel; c_id++) {
        for (cs_lnum_t j = 0; j < 6; j++) {
          st_rij[c_id][j]
            =    (t_st_rij[c_id][j] + (lag_st->npts - 1.0) * st_rij[c_id][j])
               / lag_st->npts;
        }
      }
    }

    if (st_t_e != NULL) {
      for (cs_lnum_t c_id = 0; c_id < ncel; c_id++) {
        st_t_e[c_id]
          = (t_st_t_e[c_id] + (lag_st->npts - 1.0) * st_t_e[c_id])
          / lag_st->npts;
      }
    }

    if (st_t_i != NULL) {
      for (cs_lnum_t c_id = 0; c_id < ncel; c_id++) {
        st_t_i[c_id]
          = (t_st_t_i[c_id] + (lag_st->npts - 1.0) * st_t_i[c_id])
          / lag_st->npts;
      }
    }

  }

  if (t_st_p != st_p)
    BFT_FREE(t_st_p);

  if (t_st_vel != st_vel)
    BFT_FREE(t_st_vel);

  if (t_st_imp_vel != st_imp_vel)
    BFT_FREE(t_st_imp_vel);

  if (t_st_k != st_k)
    BFT_FREE(t_st_k);

  if (t_st_rij != st_rij)
    BFT_FREE(t_st_rij);

  if (t_st_t_e != st_t_e)
    BFT_FREE(t_st_t_e);

  if (t_st_t_i != st_t_i)
    BFT_FREE(t_st_t_i);

  BFT_FREE(volp);
  BFT_FREE(volm);

  BFT_FREE(auxl);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
