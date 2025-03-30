/*============================================================================
 * Methods for particle coupling
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

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

#include "base/cs_array.h"
#include "base/cs_base.h"
#include "base/cs_math.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"

#include "bft/bft_mem.h"
#include "bft/bft_error.h"

#include "base/cs_physical_constants.h"
#include "base/cs_time_step.h"
#include "turb/cs_turbulence_model.h"

#include "lagr/cs_lagr.h"
#include "lagr/cs_lagr_particle.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "lagr/cs_lagr_coupling.h"

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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Prepare source terms for Lagrangian 2-way coupling.
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_coupling_initialize(void)
{
  cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;
  cs_lagr_source_terms_t *lag_st = cs_glob_lagr_source_terms;

  /*Note: t_* stands for temporary array, used in case of time moments */
  cs_real_t    *t_st_p = nullptr;
  cs_real_3_t  *t_st_vel= nullptr;
  cs_real_t    *t_st_imp_vel = nullptr;
  cs_real_6_t  *t_st_rij = nullptr;
  cs_real_t    *t_st_k = nullptr;
  cs_real_t    *t_st_t_e = nullptr;
  cs_real_t    *t_st_t_i = nullptr;
  cs_real_t    *volp = nullptr;
  cs_real_t    *volm = nullptr;

  cs_field_t *f_st_p = cs_field_by_name_try("lagr_st_pressure");
  cs_field_t *f_st_vel = cs_field_by_name_try("lagr_st_velocity");
  cs_field_t *f_st_imp_vel = cs_field_by_name_try("lagr_st_imp_velocity");
  cs_field_t *f_st_rij = cs_field_by_name_try("lagr_st_rij");
  cs_field_t *f_st_k = cs_field_by_name_try("lagr_st_k");
  cs_field_t *f_st_t_e = cs_field_by_name_try("lagr_st_temperature");
  cs_field_t *f_st_t_i = cs_field_by_name_try("lagr_st_imp_temperature");

  bool is_time_averaged = (   cs_glob_lagr_time_scheme->isttio == 1
                           && lag_st->npts > 0);

  /* Init particle volume and mass in cell */
  CS_MALLOC(volp, n_cells_ext, cs_real_t);
  CS_MALLOC(volm, n_cells_ext, cs_real_t);
  cs_array_real_fill_zero(n_cells_ext, volp);
  cs_array_real_fill_zero(n_cells_ext, volm);

  if (lag_st->ltsdyn == 1) {
    /* Momentum source terms
       ===================== */
    if (is_time_averaged) {
      CS_MALLOC(t_st_vel, n_cells_ext, cs_real_3_t);
      CS_MALLOC(t_st_imp_vel, n_cells_ext, cs_real_t);
    }
    else {
      if (f_st_vel != nullptr)
        t_st_vel = (cs_real_3_t *)(f_st_vel->val);
      if (f_st_imp_vel != nullptr)
        t_st_imp_vel = f_st_imp_vel->val;
    }

    cs_array_real_fill_zero(3 * n_cells_ext, (cs_real_t *) t_st_vel);
    cs_array_real_fill_zero(n_cells_ext, t_st_imp_vel);

  /* Turbulence source terms
     ======================= */

    if (extra->itytur == 2 || extra->itytur == 4 ||
        extra->itytur == 5 || extra->iturb == CS_TURB_K_OMEGA) {

      /* In v2f the Lagrangian STs only influence k and epsilon
         (difficult to write something for v2, which loses its meaning as
         "Rij component") */
      if (is_time_averaged)
        CS_MALLOC(t_st_k, n_cells_ext, cs_real_t);
      else if (f_st_k != nullptr)
        t_st_k = f_st_k->val;
      cs_array_real_fill_zero(n_cells_ext, t_st_k);
    }
    else if (extra->itytur == 3) {
      if (is_time_averaged)
        CS_MALLOC(t_st_rij, n_cells_ext, cs_real_6_t);
      else if (f_st_rij != nullptr)
        t_st_rij = (cs_real_6_t*)(f_st_rij->val);
      cs_array_real_fill_zero(n_cells_ext * 6, (cs_real_t *)t_st_rij);
    }
  }

  /* Mass source terms
     ================= */
  if (    lag_st->ltsmas == 1
      && (   cs_glob_lagr_specific_physics->solve_mass == 1
          || cs_glob_lagr_specific_physics->solve_diameter == 1
          || cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_CTWR )) {
    if (t_st_p == nullptr) {
      if (is_time_averaged)
        CS_MALLOC(t_st_p, n_cells_ext, cs_real_t);
      else if (f_st_p != nullptr)
        t_st_p = f_st_p->val;
    }
    cs_array_real_fill_zero(n_cells_ext, t_st_p);
  }

  /* Thermal source terms
     ==================== */

  if (lag_st->ltsthe == 1) {
    if (   cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_COAL
        && cs_glob_lagr_const_dim->nlayer > 1)
      bft_error(__FILE__, __LINE__, 0,
                _("Thermal coupling not implemented in multi-layer case"));

    if (  (   cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_HEAT
           && cs_glob_lagr_specific_physics->solve_temperature == 1)
        || cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_COAL
        || cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_CTWR) {

      if (is_time_averaged) {
        CS_MALLOC(t_st_t_e, n_cells_ext, cs_real_t);
        CS_MALLOC(t_st_t_i, n_cells_ext, cs_real_t);
      }
      else {
        t_st_t_e = f_st_t_e->val;
        t_st_t_i = f_st_t_i->val;
      }

      cs_array_real_fill_zero(n_cells_ext, t_st_t_e);
      cs_array_real_fill_zero(n_cells_ext, t_st_t_i);
    }
  }
  lag_st->t_st_p = t_st_p;
  lag_st->t_st_vel= t_st_vel;
  lag_st->t_st_imp_vel = t_st_imp_vel;
  lag_st->t_st_rij = t_st_rij;
  lag_st->t_st_k = t_st_k;
  lag_st->t_st_t_e = t_st_t_e;
  lag_st->t_st_t_i = t_st_t_i;
  lag_st->volp = volp;
  lag_st->volm = volm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Increment the source terms for Lagrangian 2-way coupling with
 *         quantites attached to a given particle.
 *
 * \remark  Source terms are computed for the starting cell of a particle
 *          during a given iteration. Even if particle exits the domain,
 *          it s necessary to compute a source term matching the exchange
 *          between the carrier fluid and the particle at the beginning
 *          of the time step. If nor == 2 and the particle interacts with a
 *          boundary, then the source terms are computed as if nor == 1.
 *
 * \param[in]   p_set   pointer to particle set
 * \param[in]   p_id    particle id
 * \param[in]   dt_part remaining time step associated to the particle
 * \param[in]   rebound true if a rebound occured over last trajectory step
 * \param[in]   taup    dynamic characteristic time
 * \param[in]   force_p forces per mass unit on particles (m/s^2)
 * \param[in]   tempct  thermal characteristic time
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_coupling_increment_part_contrib(cs_lagr_particle_set_t       *p_set,
                                        const cs_lnum_t               p_id,
                                        const cs_real_t               dt_part,
                                        const bool                    rebound,
                                        const cs_real_t               taup,
                                        const cs_real_3_t             force_p,
                                        const cs_real_2_t             tempct)
{
  /* WARNING : Only based on the first continuous phase */

  /* The cell_id incremented is the cell_id at the begining of the time step */
  cs_lnum_t c_id = cs_lagr_particles_get_lnum_n(p_set, p_id, 1, CS_LAGR_CELL_ID);
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  cs_real_t *cell_f_vol = mq->cell_vol;
  const int *restrict c_disable_flag = mq->c_disable_flag;
  cs_lnum_t has_dc = mq->has_disable_flag; /* Has cells disabled? */

  cs_real_t dvol = 0.;
  if (has_dc * c_disable_flag[has_dc * c_id] != 0)
    return;
  else
    dvol = 1. / cell_f_vol[c_id];

  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;
  cs_lagr_source_terms_t *lag_st = cs_glob_lagr_source_terms;

  /*Note: t_* stands for temporary array, used in case of time moments */
  cs_real_t    *t_st_p = lag_st->t_st_p;
  cs_real_3_t  *t_st_vel= lag_st->t_st_vel;
  cs_real_t    *t_st_imp_vel = lag_st->t_st_imp_vel;
  cs_real_6_t  *t_st_rij = lag_st->t_st_rij;
  cs_real_t    *t_st_k = lag_st->t_st_k;
  cs_real_t    *t_st_t_e = lag_st->t_st_t_e;
  cs_real_t    *t_st_t_i = lag_st->t_st_t_i;
  cs_real_t    *volp = lag_st->volp;
  cs_real_t    *volm = lag_st->volm;

  cs_real_t  p_stat_w = cs_lagr_particles_get_real(p_set, p_id,
                                                   CS_LAGR_STAT_WEIGHT);
  cs_real_t *p_vel    =
    cs_lagr_particles_attr_get_ptr<cs_real_t>(p_set, p_id, CS_LAGR_VELOCITY);
  cs_real_t *prev_p_vel  =
    cs_lagr_particles_attr_n_get_ptr<cs_real_t>(p_set, p_id, 1,
                                                CS_LAGR_VELOCITY);

  cs_real_t  p_mass = cs_lagr_particles_get_real(p_set, p_id, CS_LAGR_MASS);
  cs_real_t  prev_p_mass = cs_lagr_particles_get_real_n(p_set, p_id, 1,
                                                        CS_LAGR_MASS);

  cs_real_t dtp = cs_glob_lagr_time_step->dtp;

  cs_real_t rel_dt = dt_part / dtp;

  /* Preliminary computations
     ======================== */

  /* Finalization of external forces (if the particle interacts with a
     domain boundary, revert to order 1). based on the first carrier phase */

  cs_real_t aux1 = dt_part / taup;

  cs_real_t tsfext;
  if (   cs_glob_lagr_time_scheme->t_order == 1 || rebound)
    tsfext = (1.0 - exp(-aux1)) * p_mass * taup;

  else
    tsfext =  (1.0 - (1.0 - exp (-aux1)) / aux1) * taup * p_mass;

    //TODO tsfext should be computed elsewhere (in sde) and the mass of particle
    // may be the previous mass.
  cs_real_3_t auxl;
  for (cs_lnum_t i = 0; i < 3; i++)
    auxl[i] = p_stat_w * ((p_mass * p_vel[i] - prev_p_mass * prev_p_vel[i])
                           - force_p[i] * tsfext) / dtp;

  /* Momentum source terms
     ===================== */

  if (lag_st->ltsdyn == 1) {

    cs_real_t  prev_p_diam = cs_lagr_particles_get_real_n(p_set, p_id, 1,
                                                          CS_LAGR_DIAMETER);
    volp[c_id] += p_stat_w * cs_math_pi * pow(prev_p_diam, 3) / 6.0 * rel_dt;
    volm[c_id] += p_stat_w * prev_p_mass * rel_dt;

    /* Momentum source term */

    for (cs_lnum_t i = 0; i < 3; i++)
      t_st_vel[c_id][i] -= dvol * auxl[i];

    t_st_imp_vel[c_id] -= 2.0 * dvol * p_stat_w * p_mass / taup * rel_dt;

    /* Turbulence source terms
       ======================= */

    cs_real_t *prev_vel_s  =
      cs_lagr_particles_attr_n_get_ptr<cs_real_t>(p_set, p_id, 1,
                                                  CS_LAGR_VELOCITY_SEEN);
    cs_real_t *new_vel_s   =
      cs_lagr_particles_attr_get_ptr<cs_real_t>(p_set, p_id,
                                                CS_LAGR_VELOCITY_SEEN);

    cs_real_3_t vel_s= {0.5 * (prev_vel_s[0] + new_vel_s[0]),
                        0.5 * (prev_vel_s[1] + new_vel_s[1]),
                        0.5 * (prev_vel_s[2] + new_vel_s[2])};

    if (extra->itytur == 2 || extra->itytur == 4 ||
        extra->itytur == 5 || extra->iturb == CS_TURB_K_OMEGA)
      /* In v2f the Lagrangian STs only influence k and epsilon
         (difficult to write something for v2, which loses its meaning as
         "Rij component") */
      t_st_k[c_id] -= dvol * cs_math_3_dot_product(vel_s, auxl);

    else if (extra->itytur == 3) {
      for (cs_lnum_t ij = 0; ij < 6; ij++) {
        cs_lnum_t i = _iv2t[ij];
        cs_lnum_t j = _jv2t[ij];

        t_st_rij[c_id][ij] -= ( vel_s[i] * auxl[j]
                              + vel_s[j] * auxl[i])*dvol;
      }
    }
  }

  /* Mass source terms
     ================= */

  if (    lag_st->ltsmas == 1
      && (   cs_glob_lagr_specific_physics->solve_mass == 1
          || cs_glob_lagr_specific_physics->solve_diameter == 1
          || cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_CTWR )) {

      t_st_p[c_id] += - p_stat_w * (p_mass - prev_p_mass) / dtp * dvol;
  }

  /* Thermal source terms
     ==================== */

  if (   lag_st->ltsthe == 1
      && (  (    cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_HEAT
            && cs_glob_lagr_specific_physics->solve_temperature == 1)
          || cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_COAL
          || cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_CTWR)) {

    cs_real_t  p_cp = cs_lagr_particles_get_real_n(p_set, p_id, 0,
                                                  CS_LAGR_CP);
    cs_real_t  prev_p_cp = cs_lagr_particles_get_real_n(p_set, p_id, 1,
                                                       CS_LAGR_CP);
    cs_real_t  p_tmp = cs_lagr_particles_get_real_n(p_set, p_id, 0,
                                                   CS_LAGR_TEMPERATURE);
    cs_real_t  prev_p_tmp = cs_lagr_particles_get_real_n(p_set, p_id, 1,
                                                        CS_LAGR_TEMPERATURE);


      t_st_t_e[c_id] += - (p_mass * p_tmp * p_cp
                          - prev_p_mass * prev_p_tmp * prev_p_cp
                          ) / dtp * p_stat_w * dvol;
      t_st_t_i[c_id] += p_stat_w * p_mass * p_cp * dvol
                      / tempct[1];

    if (   cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_HEAT
        && cs_glob_lagr_specific_physics->solve_temperature == 1
        && extra->radiative_model > 0) {

      cs_real_t  p_diam = cs_lagr_particles_get_real_n(p_set, p_id, 0,
                                                      CS_LAGR_DIAMETER);
      cs_real_t  p_eps = cs_lagr_particles_get_real_n(p_set, p_id, 0,
                                                     CS_LAGR_EMISSIVITY);
      cs_real_t aux2 = cs_math_pi * p_diam * p_diam * p_eps * dvol
                      * (extra->rad_energy->val[c_id]
                         - 4.0 * _c_stephan * cs_math_pow4(p_tmp)) * rel_dt;
      t_st_t_e[c_id] += aux2 * p_stat_w * rel_dt;

    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize source terms for Lagrangian 2-way coupling by treating
 *         cell-attached properties.
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_coupling_finalize(void)
{
  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;
  cs_lnum_t ncel = cs_glob_mesh->n_cells;

  cs_lagr_source_terms_t *lag_st = cs_glob_lagr_source_terms;
  bool is_time_averaged = (   cs_glob_lagr_time_scheme->isttio == 1
                           && lag_st->npts > 0);

  /*Note: t_* stands for temporary array, used in case of time moments */
  cs_real_t    *t_st_p = lag_st->t_st_p;
  cs_real_3_t  *t_st_vel = lag_st->t_st_vel;
  cs_real_t    *t_st_imp_vel = lag_st->t_st_imp_vel;
  cs_real_6_t  *t_st_rij = lag_st->t_st_rij;
  cs_real_t    *t_st_k = lag_st->t_st_k;
  cs_real_t    *t_st_t_e = lag_st->t_st_t_e;
  cs_real_t    *t_st_t_i = lag_st->t_st_t_i;
  cs_real_t    *volp = lag_st->volp;
  cs_real_t    *volm = lag_st->volm;

  /* Number of passes for steady source terms */
  if (   cs_glob_lagr_time_scheme->isttio == 1
      && cs_glob_time_step->nt_cur >= cs_glob_lagr_source_terms->nstits)
    lag_st->npts += 1;

  lag_st->ntxerr = 0;
  lag_st->vmax = 0.0;
  lag_st->tmamax = 0.0;

  /* Check that the maximum admissible volume fraction of particles
     is not exceeded in some cells.
     ============================================================== */

  const cs_real_t tvmax = 0.8;
  const cs_real_t *cell_vol = cs_glob_mesh_quantities->cell_vol;

  if (cs_glob_lagr_source_terms->ltsdyn == 1) {
    cs_real_3_t *vel = (cs_real_3_t *)extra->vel->val;
    if (extra->itytur == 2 || extra->itytur == 4 ||
        extra->itytur == 5 || extra->iturb == CS_TURB_K_OMEGA) {

      /* In v2f the Lagrangian STs only influence k and epsilon
         (difficult to write something for v2, which loses its meaning as
       "Rij component") */

      for (cs_lnum_t c_id = 0; c_id < ncel; c_id++)
        t_st_k[c_id] -= cs_math_3_dot_product(vel[c_id], t_st_vel[c_id]);
    }
    else if (extra->itytur == 3) {
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
  for (cs_lnum_t c_id = 0; c_id < ncel; c_id++) {

    cs_real_t mf   = cell_vol[c_id] * extra->cromf->val[c_id];
    cs_real_t tauv = volp[c_id] / cell_vol[c_id];
    cs_real_t taum = volm[c_id] / mf;

    if (tauv > tvmax) {
      lag_st->ntxerr++;
      cs_glob_lagr_source_terms->ntxerr++;;

      if (t_st_vel != nullptr) {
        for (cs_lnum_t j = 0; j < 3; j++)
          t_st_vel[c_id][j] = 0.0;
      }

      if (t_st_imp_vel != nullptr)
        t_st_imp_vel[c_id] = 0.0;

      if (t_st_k != nullptr)
        t_st_k[c_id] = 0.0;

      if (t_st_rij != nullptr) {
        for (cs_lnum_t j = 0; j < 6; j++)
          t_st_rij[c_id][j] = 0.0;
      }

      if (t_st_t_e != nullptr)
        t_st_t_e[c_id] = 0.0;

      if (t_st_t_i != nullptr)
        t_st_t_i[c_id] = 0.0;
    }
    lag_st->vmax = cs::max(tauv, lag_st->vmax);
    lag_st->tmamax = cs::max(taum, lag_st->tmamax);
  }

  /* Time average of source terms
     ============================ */

  cs_real_t   *st_p       = nullptr;
  cs_real_3_t *st_vel     = nullptr;
  cs_real_t   *st_imp_vel = nullptr;
  cs_real_6_t *st_rij     = nullptr;
  cs_real_t   *st_k       =  nullptr;
  cs_real_t   *st_t_e     =  nullptr;
  cs_real_t   *st_t_i     =  nullptr;

  if (is_time_averaged) {
    cs_field_t *f_st_p = cs_field_by_name_try("lagr_st_pressure");
    cs_field_t *f_st_vel = cs_field_by_name_try("lagr_st_velocity");
    cs_field_t *f_st_imp_vel = cs_field_by_name_try("lagr_st_imp_velocity");
    cs_field_t *f_st_rij = cs_field_by_name_try("lagr_st_rij");
    cs_field_t *f_st_k = cs_field_by_name_try("lagr_st_k");
    cs_field_t *f_st_t_e = cs_field_by_name_try("lagr_st_temperature");
    cs_field_t *f_st_t_i = cs_field_by_name_try("lagr_st_imp_temperature");

    if (f_st_p != nullptr) {
      st_p = f_st_p->val;
      for (cs_lnum_t c_id = 0; c_id < ncel; c_id++) {
        st_p[c_id]
          = (t_st_p[c_id] + (lag_st->npts - 1.0) * st_p[c_id])
          / lag_st->npts;
      }
    }

    if (f_st_vel != nullptr) {
      st_vel = (cs_real_3_t*)(f_st_vel->val);
      for (cs_lnum_t c_id = 0; c_id < ncel; c_id++) {
        for (cs_lnum_t j = 0; j < 3; j++) {
          st_vel[c_id][j]
            =    (t_st_vel[c_id][j] + (lag_st->npts - 1.0) * st_vel[c_id][j])
               / lag_st->npts;
        }
      }
    }

    if (f_st_imp_vel != nullptr) {
      st_imp_vel = f_st_imp_vel->val;
      for (cs_lnum_t c_id = 0; c_id < ncel; c_id++) {
        st_imp_vel[c_id]
          = (t_st_imp_vel[c_id] + (lag_st->npts - 1.0) * st_imp_vel[c_id])
          / lag_st->npts;
      }
    }

    if (f_st_k != nullptr) {
      st_k = f_st_k->val;
      for (cs_lnum_t c_id = 0; c_id < ncel; c_id++) {
        st_k[c_id]
          = (t_st_k[c_id] + (lag_st->npts - 1.0) * st_k[c_id])
          / lag_st->npts;
      }
    }

    if (f_st_rij != nullptr) {
      st_rij = (cs_real_6_t*)f_st_rij->val;
      for (cs_lnum_t c_id = 0; c_id < ncel; c_id++) {
        for (cs_lnum_t j = 0; j < 6; j++) {
          st_rij[c_id][j]
            =    (t_st_rij[c_id][j] + (lag_st->npts - 1.0) * st_rij[c_id][j])
               / lag_st->npts;
        }
      }
    }

    if (f_st_t_e != nullptr) {
      st_t_e = f_st_t_e->val;
      for (cs_lnum_t c_id = 0; c_id < ncel; c_id++) {
        st_t_e[c_id]
          = (t_st_t_e[c_id] + (lag_st->npts - 1.0) * st_t_e[c_id])
          / lag_st->npts;
      }
    }

    if (f_st_t_i != nullptr) {
      st_t_i = f_st_t_i->val;
      for (cs_lnum_t c_id = 0; c_id < ncel; c_id++) {
        st_t_i[c_id]
          = (t_st_t_i[c_id] + (lag_st->npts - 1.0) * st_t_i[c_id])
          / lag_st->npts;
      }
    }

    if (t_st_p != st_p)
      CS_FREE(t_st_p);

    if (t_st_vel != st_vel)
      CS_FREE(t_st_vel);

    if (t_st_imp_vel != st_imp_vel)
      CS_FREE(t_st_imp_vel);

    if (t_st_k != st_k)
      CS_FREE(t_st_k);

    if (t_st_rij != st_rij)
      CS_FREE(t_st_rij);

    if (t_st_t_e != st_t_e)
      CS_FREE(t_st_t_e);

    if (t_st_t_i != st_t_i)
      CS_FREE(t_st_t_i);
  }

  CS_FREE(volp);
  CS_FREE(volm);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
