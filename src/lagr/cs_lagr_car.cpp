/*============================================================================
 * Compute particle characteristics: Tp, TL and PI
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
 * Functions dealing with particle tracking
 *============================================================================*/

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

#include "bft/bft_printf.h"
#include "bft/bft_error.h"
#include "bft/bft_mem.h"

#include "fvm/fvm_periodicity.h"

#include "atmo/cs_atmo.h"
#include "base/cs_base.h"
#include "base/cs_defs.h"
#include "base/cs_math.h"
#include "base/cs_field_operator.h"
#include "base/cs_halo.h"
#include "base/cs_interface.h"
#include "base/cs_math.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_order.h"
#include "base/cs_parall.h"
#include "base/cs_search.h"
#include "base/cs_timer_stats.h"
#include "base/cs_thermal_model.h"
#include "turb/cs_turbulence_model.h"
#include "base/cs_velocity_pressure.h"

#include "base/cs_field.h"

#include "base/cs_physical_constants.h"
#include "pprt/cs_physical_model.h"

#include "lagr/cs_lagr.h"
#include "lagr/cs_lagr_new.h"
#include "lagr/cs_lagr_particle.h"
#include "lagr/cs_lagr_stat.h"
#include "lagr/cs_lagr_precipitation_model.h"
#include "lagr/cs_lagr_prototypes.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "lagr/cs_lagr_car.h"

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute particle characteristics: Tp, TL and PI as well as
 * covariance and variance tensors for the stochasit model
 *
 * \param[in] iprev             time step indicator for fields
 *                                0: use fields at current time step
 *                                1: use fields at previous time step
 * \param[in]  phase_id         carrier phase id
 * \param[in]  dt               time step (per cell)
 * \param[out] taup             dynamic characteristic time
 * \param[out] tlag             fluid characteristic time
 * \param[out] piil             term in integration of up sde
 * \param[out] bx               turbulence characteristics
 * \param[out] tempct           thermal characteristic time
 * \param[out] beta             for the extended scheme
 * \param[in]  gradvf           fluid velocity gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_car(int              iprev,
            int              phase_id,
            const cs_real_t  dt[],
            cs_real_t        *taup,
            cs_real_3_t      *tlag,
            cs_real_3_t      *piil,
            cs_real_33_t     *bx,
            cs_real_t        tempct[],
            cs_real_3_t      *beta,
            cs_real_33_t     gradvf[])
{
  /* Particles management */

  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;

  /* Mesh */

  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_lnum_t n_cells = mesh->n_cells;
  cs_lnum_t ncelet = mesh->n_cells_with_ghosts;

  cs_lagr_extra_module_t *extra_i = cs_get_lagr_extra_module();
  cs_lagr_extra_module_t *extra = extra_i;
  iprev = CS_MIN(iprev, extra->vel->n_time_vals -1);

  /* Initialization
     ---------------*/

  bool turb_disp_model = false;
  if (   cs_glob_lagr_model->modcpl > 0
      && cs_glob_time_step->nt_cur > cs_glob_lagr_model->modcpl
      && cs_glob_time_step->nt_cur > cs_glob_lagr_stat_options->idstnt)
    turb_disp_model = true;

  cs_lnum_t nor = cs_glob_lagr_time_step->nor;

  cs_real_t rec    = 1000.0;
  cs_real_t c0     = cs_turb_crij_c0;

  cs_real_t cl     = 1.0 / (0.5 + 0.75 * c0);
  cs_real_t cb     = 0.8;
  cs_real_t cbcb   = 0.64;
  cs_real_t d6spi  = 6.0 / cs_math_pi;
  cs_real_t d1s3   = 1.0 / 3.0;

  const cs_real_t *grav = cs_glob_physical_constants->gravity;

  cs_real_t diftl0 = -1;
  if (   cs_glob_physical_model_flag[CS_COMBUSTION_EBU] == 0
      || cs_glob_physical_model_flag[CS_COMBUSTION_EBU] == 2)
    diftl0 = cs_field_get_key_double(cs_field_by_name("enthalpy"),
                                     cs_field_key_id("diffusivity_ref"));

  /* Compute Tp and Tc in case of thermal model
     -------------------------------------------*/

  for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

    unsigned char *particle = p_set->p_buffer + p_am->extents * ip;

    cs_real_t      p_diam   = cs_lagr_particle_get_real(particle, p_am,
                                                        CS_LAGR_DIAMETER);
    cs_lnum_t      cell_id  = cs_lagr_particle_get_lnum(particle, p_am,
                                                        CS_LAGR_CELL_ID);

    /* FIXME we may still need to do computations here */
    if (cs_lagr_particles_get_flag(p_set, ip, CS_LAGR_PART_FIXED))
      continue;

    cs_real_t p_mass = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS);
    cs_real_t p_rom  = p_mass * d6spi / pow(p_diam, 3.0);

    cs_real_t rom, xnul, rel_vel_norm = 0, rep, fdr;
    rom = extra_i[phase_id].cromf->val[cell_id];
    xnul = extra_i[phase_id].viscl->val[cell_id] / rom;

    cs_real_t *part_vel_seen
      = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                 CS_LAGR_VELOCITY_SEEN);
    cs_real_t *part_vel
      = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                 CS_LAGR_VELOCITY);

    for (int idim_ = 0; idim_ < 3; idim_++) {
      rel_vel_norm += (part_vel_seen[3*phase_id + idim_] - part_vel[idim_])
                    * (part_vel_seen[3*phase_id + idim_] - part_vel[idim_]);
    }

    rel_vel_norm = sqrt(rel_vel_norm);

    /* Compute the local Reynolds number */
    rep = rel_vel_norm * p_diam / xnul;

    cs_real_t d2 = cs_math_sq(p_diam); /* drag coefficient */

    if (rep <= rec) {
      fdr = 18.0 * xnul * (1.0 + 0.15 * pow (rep, 0.687)) / d2;
    }
    else
      fdr = 0.44 * 3.0 / 4.0 * rel_vel_norm / p_diam;

    /* Tp computation */
    taup[ip] = p_rom / rom / fdr;

    /* Added-mass term? */
    if (cs_glob_lagr_time_scheme->iadded_mass == 1){
      taup[ip] *= (  1.0
                   + 0.5 * cs_glob_lagr_time_scheme->added_mass_const
                         * rom / p_rom);
    }

    /* Tp user computation */

    cs_user_lagr_rt(ip, rep, rel_vel_norm, rom, p_rom, xnul, taup, dt);

    /* Tc computation : based on the first carrier flow in ncfd */

    if (   (   cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_HEAT
            && cs_glob_lagr_specific_physics->itpvar == 1)
        || cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_COAL
        || cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_CTWR ) {

      /* Fluid Cp */
      cs_real_t xcp;
      if (extra->icp > 0)
        xcp = extra->cpro_cp->val[0];
      else
        xcp = cs_glob_fluid_properties->cp0;

      /* Local Nu computation */
      /* a priori in gas or pulverized coal combustion,
         diffusivity is always constant */
      cs_real_t xrkl;
      if (diftl0 >= 0)
        xrkl = diftl0 / rom;
      else if (extra->cpro_viscls != nullptr)
        xrkl = extra->cpro_viscls->val[cell_id] / (rom * xcp);
      else
        xrkl = extra->visls0 / (rom * xcp);

      cs_real_t prt  = xnul / xrkl;
      cs_real_t fnus = 2.0 + 0.55 * pow (rep, 0.5) * pow (prt, (d1s3));

      cs_real_t p_cp = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_CP);

      /* Thermal characteristic time Tc computation */
      tempct[ip] = d2 * p_rom * p_cp / (fnus * 6.0 * rom * xcp * xrkl);

      /* User computation for Tc */
      cs_user_lagr_rt_t(ip, rep, rel_vel_norm, rom, p_rom,
                        xnul, xcp, xrkl, tempct, dt);

      /* Implicit source term for return thermal coupling */
      tempct[p_set->n_particles + ip]
        = fnus * cs_math_pi * p_diam * xrkl * rom;
    }
  }

  /* Compute TL
     ---------- */

  /* Compute turbulent kinetic energy and dissipation
     based on turbulence model */

  if (cs_glob_lagr_model->idistu == 1) {
    /* Compute properties associated to the Eulerian mesh
     * used in the complete model for discrete particles
       ---------- */

    cs_real_3_t *cell_tlag_et = nullptr;
    cs_real_3_t *cell_grad_lagr_time_r = nullptr;
    cs_real_3_t *cell_bbi = nullptr;
    cs_real_3_t *cell_bxi_tl = nullptr;
    cs_real_t mean_uvwdif;
    cs_real_3_t dir;

    cs_real_t *energi = extra_i[phase_id].cvar_k->val;
    cs_real_t *lagr_time = extra_i[phase_id].lagr_time->val;
    cs_real_6_t *rij_pre = nullptr;
    if (extra_i[phase_id].itytur == 3)
      rij_pre = (cs_real_6_t *)extra_i[phase_id].cvar_rij->vals[iprev];

    /* Determine the quantities associated to the cell in the local ref. */
    if (turb_disp_model) {
      BFT_MALLOC(cell_bbi, n_cells, cs_real_3_t);
      BFT_MALLOC(cell_bxi_tl, n_cells, cs_real_3_t);
      BFT_MALLOC(cell_tlag_et, n_cells, cs_real_3_t);

      /* compute grad_lagr_time in the local referential associated to the cell
       else we set bbi = {1.,1.,1.} and stay in the absolute referential */
      if (beta != nullptr) {
        BFT_MALLOC(cell_grad_lagr_time_r, n_cells, cs_real_3_t);
      }

      int stat_type = cs_lagr_stat_type_from_attr_id(CS_LAGR_VELOCITY);
      cs_field_t *stat_vel
        = cs_lagr_stat_get_moment(stat_type,
                                  CS_LAGR_STAT_GROUP_PARTICLE,
                                  CS_LAGR_MOMENT_MEAN,
                                  0,
                                  -1);

      cs_field_t *stat_w = cs_lagr_stat_get_stat_weight(0);

      for (cs_lnum_t cell_id = 0; cell_id < cs_glob_mesh->n_cells; cell_id++) {
        /* Initialize cell_tlag_et : */
        for (int i = 0; i < 3; i++)
          cell_tlag_et[cell_id][i] = 0.;

        if (stat_w->val[cell_id] > cs_glob_lagr_stat_options->threshold) {
          /* compute mean relative velocity <Up> - <Uf>*/
          for (int i = 0; i < 3; i++)
            dir[i] = stat_vel->val[cell_id * 3 + i]
                   - extra_i[phase_id].vel->vals[iprev][cell_id * 3 + i];

          /* Compute and store the mean relative velocity |<U_r>| = |<Up>-Uf|*/
          mean_uvwdif = cs_math_3_square_norm(dir);
          mean_uvwdif = (3.0 * mean_uvwdif) / (2.0 * energi[cell_id]);

          /* FIXME add proper isotropic behavior */

          /* Crossing trajectory in the direction of "<u_f>-<u_p>"
           * and in the span-wise direction */
          cs_real_t an, at;
          an = (1.0 + cbcb * mean_uvwdif);
          at = (1.0 + 4.0 * cbcb * mean_uvwdif);

          cell_bbi[cell_id][0] = sqrt(an); /* First direction, n,
                                    in the local reference frame */
          cell_bbi[cell_id][1] = sqrt(at); /* Second and third direction,
                                    orthogonal to n */
          cell_bbi[cell_id][2] = sqrt(at);

          /* Compute the timescale in parallel and transverse directions */
          for (int id = 0; id < 3; id++)
            cell_tlag_et[cell_id][id] = lagr_time[cell_id]
                                      / cell_bbi[cell_id][id];

          /* Compute the main direction in the global reference
           * frame */
           cs_math_3_normalize(dir, dir);

          /* Compute k_tilde in each cell */
          cs_real_t ktil;
          if (extra_i[phase_id].itytur == 3) {
            cs_real_t *rij = rij_pre[cell_id];
            /* Note that n.R.n = R : n(x)n */
            cs_real_t rnn = cs_math_3_sym_33_3_dot_product(dir, rij, dir);
            cs_real_t tr_r = cs_math_6_trace(rij);
            // bbn * R : n(x)n + bbt * R : (1 - n(x)n)
            ktil = 3.0 * (rnn * cell_bbi[cell_id][0]
                         + (tr_r -rnn) * cell_bbi[cell_id][1])
                 / (2.0 * (cell_bbi[cell_id][0]
                          + cell_bbi[cell_id][1] + cell_bbi[cell_id][2]));
            /* cell_bbi[1] == cell_bbi[2] is used */
          }
          else if (   extra_i[phase_id].itytur == 2
              || extra_i[phase_id].itytur == 4
              || extra_i[phase_id].itytur == 5
              || extra_i[phase_id].iturb == CS_TURB_K_OMEGA) {
            ktil  = energi[cell_id];
          }
          for (int i = 0; i <3; i++) {
            cell_bxi_tl[cell_id][i] = cl * (  (c0  * ktil )
                + ((ktil- energi[cell_id] / cell_bbi[cell_id][i])
                  * 2.0 / 3.0));
            cell_bxi_tl[cell_id][i] =
                CS_MAX(cell_bxi_tl[cell_id][i], cs_math_epzero);
          }

          if (cell_grad_lagr_time_r != nullptr) {
            cs_real_33_t trans_m;
            /* Rotate the frame of reference with respect to the
             * mean relative velocity direction.
             * This referential differs the referential associated directly
             * to the particle for non spheric particle
             * TODO extend extended scheme to non spheric particle*/

            // The rotation axis is the result of the cross product between
            // the new direction vector and the main axis.
            cs_real_t n_rot[3];
            /* the direction in the local reference frame "_r" is (1, 0, 0)
             * by convention */
            const cs_real_t dir_r[3] = {1.0, 0.0, 0.0};

            // Use quaternion (cos(theta/2), sin(theta/2) n_rot)
            // where n_rot = dir ^ dir_r normalised
            // so also       dir ^ (dir + dir_r)
            //
            // cos(theta/2) = || dir + dir_r|| / 2
            cs_real_t dir_p_dir_r[3] = {dir[0] + dir_r[0],
                                        dir[1] + dir_r[1],
                                        dir[2] + dir_r[2]};
            cs_real_t dir_p_dir_r_normed[3];
            cs_math_3_normalize(dir_p_dir_r, dir_p_dir_r_normed);

            /* dir ^(dir + dir_r) / || dir + dir_r|| = sin(theta/2) n_rot
             * for the quaternion */
            cs_math_3_cross_product(dir, dir_p_dir_r_normed, n_rot);

            /* quaternion, could be normalized afterwards
             *
             * Note that the division seems stupid but is not
             * in case of degenerated case where dir is null
             * */
            const cs_real_t euler[4]
              = {  cs_math_3_norm(dir_p_dir_r)
                 / (cs_math_3_norm(dir) + cs_math_3_norm(dir_r)),
                   n_rot[0],
                   n_rot[1],
                   n_rot[2]};

            trans_m[0][0] = 2.*(euler[0]*euler[0]+euler[1]*euler[1]-0.5);
            trans_m[0][1] = 2.*(euler[1]*euler[2]+euler[0]*euler[3]);
            trans_m[0][2] = 2.*(euler[1]*euler[3]-euler[0]*euler[2]);
            trans_m[1][0] = 2.*(euler[1]*euler[2]-euler[0]*euler[3]);
            trans_m[1][1] = 2.*(euler[0]*euler[0]+euler[2]*euler[2]-0.5);
            trans_m[1][2] = 2.*(euler[2]*euler[3]+euler[0]*euler[1]);
            trans_m[2][0] = 2.*(euler[1]*euler[3]+euler[0]*euler[2]);
            trans_m[2][1] = 2.*(euler[2]*euler[3]-euler[0]*euler[1]);
            trans_m[2][2] = 2.*(euler[0]*euler[0]+euler[3]*euler[3]-0.5);

            /* transform grad(Tl) in the local
             * reference frame */

          cs_math_33_3_product
            (trans_m, extra_i[phase_id].grad_lagr_time[cell_id],
             cell_grad_lagr_time_r[cell_id]);
          }
        }
        else {
          for (int i = 0; i < 3; i++) {
            cell_bbi[cell_id][i] = 1.;
            cell_bxi_tl[cell_id][i] = cl * c0 * energi[cell_id];
          }
          if (cell_tlag_et != nullptr) {
            for (int i = 0; i < 3; i++) {
              cell_bbi[cell_id][i] = 1.;
              cell_bxi_tl[cell_id][i] = cl * c0 * energi[cell_id];
            }
            if (cell_tlag_et != nullptr) {
              for (int i = 0; i < 3; i++) {
                cell_tlag_et[cell_id][i] = lagr_time[cell_id];
              }
            }
            if (cell_grad_lagr_time_r != nullptr) {
              for (int i = 0; i < 3; i++) {
                cell_grad_lagr_time_r[cell_id][i]
                  = extra_i[phase_id].grad_lagr_time[cell_id][i];
              }
            }
          }
        }
      }//end loop on cell
    }//endif init array associated to mean quantities for turb_disp_model

    /* -> Calculation of TL BX and potentially beta     */

    for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

      if (cs_lagr_particles_get_flag(p_set, ip, CS_LAGR_PART_FIXED)) {
        for (int id = 0; id < 3; id++) {
          tlag[ip][id]      = cs_math_epzero;
          bx[ip][id][nor-1] = 0.0;
        }
      }

      unsigned char *particle = p_set->p_buffer + p_am->extents * ip;

      cs_lnum_t cell_id = cs_lagr_particle_get_lnum(particle, p_am,
                                                    CS_LAGR_CELL_ID);

      if (   lagr_time[cell_id] > cs_math_epzero
          && energi[cell_id] > cs_math_epzero) {
        if (turb_disp_model) {

          for (cs_lnum_t id = 0; id < 3; id++)
            tlag[ip][id] = cell_tlag_et[cell_id][id];

          for (cs_lnum_t id = 0; id < 3; id++) {
            tlag[ip][id] = CS_MAX(tlag[ip][id], cs_math_epzero);
            cs_real_t bxi = cell_bxi_tl[cell_id][id] / tlag[ip][id];
            if (bxi > 0.0)
              bx[ip][id][nor-1] = sqrt(bxi);
            else
              bx[ip][id][nor-1] = 0.0;
          }

          /* Compute beta_i in the local referential
           *
           * Note: the extended time scheme should
           * compute grad(Tl*_i) = grad(Tl)/b_i - Tl grad(b_i)/b_i^2
           *
           * In the following, only the first term is kept
           * to avoid computing grad(b_i) where no particles are present.
           *
           * The other possibility would be to compute b_i everywhere
           * (taking <u_pi> from the statistic "correctly" initialized)
           *
           */
          if (beta != nullptr) {

            for (cs_lnum_t id = 0; id < 3; id++) {
              if (CS_ABS(tlag[ip][id] - taup[ip])
                  < cs_math_epzero * taup[ip])
                beta[ip][id] = 0.;
              else
                beta[ip][id] = cell_grad_lagr_time_r[cell_id][id]
                               / cell_bbi[cell_id][id]
                               * cs_math_pow2(bx[ip][id][nor-1])
                               / (tlag[ip][id]-taup[ip]);
            }
          }
        }
        else { // fluid particles
          tlag[ip][0] = lagr_time[cell_id];
          tlag[ip][0] = CS_MAX(tlag[ip][0], cs_math_epzero);

          if (cs_glob_lagr_model->idiffl == 0) {
            cs_real_t *part_vel_seen
              = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                         CS_LAGR_VELOCITY_SEEN);
            cs_real_t *part_vel
              = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                         CS_LAGR_VELOCITY);

            cs_real_t uvwdif = 0.;
            for (int i =0; i < 3; i++)
              uvwdif += cs_math_sq(  part_vel_seen[3 * phase_id + i]
                                   - part_vel[i]);
            uvwdif = sqrt((3.0 * uvwdif) / (2.0 * energi[cell_id]));
            tlag[ip][0] = lagr_time[cell_id]
                          / (1.0 + cb * uvwdif);
          }

          tlag[ip][1] = tlag[ip][0];
          tlag[ip][2] = tlag[ip][0];

          bx[ip][0][nor-1]
            = sqrt (c0 * cl * energi[cell_id] / tlag[ip][0]);
          bx[ip][1][nor-1] = bx[ip][0][nor-1];
          bx[ip][2][nor-1] = bx[ip][0][nor-1];

          /* compute beta without in the global referential
           * under the assumption that bbi={1.,1.,1.}*/
          if ( beta != nullptr ) {
            for (cs_lnum_t id = 0; id < 3; id++) {
              if (CS_ABS(tlag[ip][id] - taup[ip])
                  < cs_math_epzero * taup[ip])
                beta[ip][id] = 0.;
              else
                beta[ip][id] = extra_i[phase_id].grad_lagr_time[cell_id][id]
                             * cs_math_pow2(bx[ip][id][nor-1])
                             / (tlag[ip][id]-taup[ip]);
            }
          }
        }
      }
      else {

        /* FIXME we may still need to do computations here */
        if (cs_lagr_particles_get_flag(p_set, ip, CS_LAGR_PART_FIXED))
          continue;

        for (int id = 0; id < 3; id++) {
          tlag[ip][id] = cs_math_epzero;
          bx[ip][id][nor-1] = 0.0;
        }
        if (beta != nullptr) {
          for (int id = 0; id < 3; id++)
            beta[ip][id] = 0.;
        }
      }
    }

    if ( cell_tlag_et != nullptr)
      BFT_FREE(cell_tlag_et);
    if (turb_disp_model){
      BFT_FREE(cell_bbi);
      BFT_FREE(cell_bxi_tl);
    }
    BFT_FREE(cell_grad_lagr_time_r);
  } // end if idistu == 1

  else {
    for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {
      /* FIXME we may still need to do computations here */
      if (cs_lagr_particles_get_flag(p_set, ip, CS_LAGR_PART_FIXED))
        continue;

      for (cs_lnum_t id = 0; id < 3; id++ ) {
        tlag[ip][id] = cs_math_epzero;
        bx[ip][id][nor-1] = 0.0;
      }
      if (beta != nullptr) {
        for (int id = 0; id < 3; id++)
          beta[ip][id] = 0.;
      }
    }
  }

  /* Before computing Pii: compute gradients of covariance */
  /* Various initializations  */
  if (cs_glob_lagr_model->cs_used == 0) {
    cs_var_cal_opt_t var_cal_opt;
    int key_cal_opt_id = cs_field_key_id("var_cal_opt");
    cs_field_get_key_struct(extra->pressure, key_cal_opt_id, &var_cal_opt);

    /* Now compute the gradients */
    cs_real_t *f_inter_cov;
    BFT_MALLOC(f_inter_cov, ncelet, cs_real_t);

    /* Get the variable we want to compute the gradients
     * from (covariance velocity seen/velocity) */
    int stat_type = cs_lagr_stat_type_from_attr_id
                      (CS_LAGR_VELOCITY_SEEN_VELOCITY_COV);

    cs_field_t *stat_cov_skp
      = cs_lagr_stat_get_moment(stat_type,
                                CS_LAGR_STAT_GROUP_PARTICLE,
                                CS_LAGR_MOMENT_MEAN,
                                0,
                                -phase_id-1);

    for (int i = 0; i < 9; i++){
      for (int iel_ = 0; iel_ < n_cells; iel_++){
        f_inter_cov[iel_] = stat_cov_skp->val[9 * iel_ + i];
      }
      if (mesh->halo != nullptr)
        cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, f_inter_cov);

      cs_gradient_scalar("intermediate cov skp gradient [Lagrangian module]",
                          CS_GRADIENT_GREEN_ITER,
                          CS_HALO_STANDARD,
                          1,
                          0,
                          0,
                          1,
                          var_cal_opt.verbosity,
                          static_cast<cs_gradient_limit_t>(var_cal_opt.imligr),
                          var_cal_opt.epsrgr,
                          var_cal_opt.climgr,
                          0,
                          nullptr,
                          f_inter_cov,
                          nullptr,
                          nullptr,
                          extra_i[phase_id].grad_cov_skp[i]);
    }
    /* Get the variable we want to compute the gradients from (velocity seen) */
    stat_type = cs_lagr_stat_type_from_attr_id(CS_LAGR_VELOCITY_SEEN);

    cs_field_t *stat_cov_sk
      = cs_lagr_stat_get_moment(stat_type,
                                CS_LAGR_STAT_GROUP_PARTICLE,
                                CS_LAGR_MOMENT_VARIANCE,
                                0,
                                -phase_id-1);

    for (int i = 0; i < 6; i++) {
      for (int iel_ = 0; iel_ < n_cells; iel_++){
        f_inter_cov[iel_] = stat_cov_sk->val[6 * iel_ + i];
      }
      if (mesh->halo != nullptr)
        cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, f_inter_cov);

      cs_gradient_scalar("intermediate cov sk gradient [Lagrangian module]",
                          CS_GRADIENT_GREEN_ITER,
                          CS_HALO_STANDARD,
                          1,
                          0, // n_r_sweeps (number of reconstruction sweeps)
                          0,
                          1,
                          var_cal_opt.verbosity,
                          static_cast<cs_gradient_limit_t>(var_cal_opt.imligr),
                          var_cal_opt.epsrgr,
                          var_cal_opt.climgr,
                          0,
                          nullptr,
                          f_inter_cov,
                          nullptr,
                          nullptr,
                          extra_i[phase_id].grad_cov_sk[i]);
    }

    BFT_FREE(f_inter_cov);
  }

  /* Compute Pii
     ----------- */

  if (turb_disp_model && cs_glob_lagr_model->cs_used == 0) {
    /* add grad(<Vf>)*(<Up>-<Us>) if there are enough particles*/
    for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

      unsigned char *particle = p_set->p_buffer + p_am->extents * ip;

      cs_lnum_t cell_id = cs_lagr_particle_get_lnum(particle, p_am,
                                                    CS_LAGR_CELL_ID);
      cs_real_t *part_vel_seen
        = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                   CS_LAGR_VELOCITY_SEEN);

      /* New piil for the fluctuation equation */

      /* Initialize piil to 0. */
      for (cs_lnum_t id = 0; id < 3; id++){
        piil[ip][id] = 0.;
      }

      for (cs_lnum_t id = 0; id < 3; id++) {
        for (int j = 0; j < 3; j++){
          if (extra_i[phase_id].iturb == CS_TURB_RIJ_EPSILON_SSG) {
            /* Arcen & Tanière model (2009) */
            piil[ip][id]
              += extra_i[phase_id].grad_cov_skp[3 * id + j][cell_id][j];
          }
          else if (extra_i[phase_id].iturb == CS_TURB_K_EPSILON) {
            /* Arcen & Tanière model  (2009) */
            piil[ip][id]
              += extra_i[phase_id].grad_cov_skp[3 * id + j][cell_id][j];
          }
          else {
            bft_error
              (__FILE__, __LINE__, 0,
              _("Lagrangian turbulent dispersion is not compatible with\n"
                "the selected turbulence model.\n"
                "\n"
                "Turbulent dispersion is taken into account with idistu = %d\n"
                " Activated turbulence model is %d, when only k-eps, LES, "
                "Rij-eps,\n"
                " V2f or k-omega are handled."),
              (int)cs_glob_lagr_model->idistu,
              (int)extra_i[phase_id].iturb);
          }
        }
        for (int j = 0; j < 3; j++){
          /* Instantaneous velocity model */
          piil[ip][id] -= (part_vel_seen[3 * phase_id + j]
              - extra_i[phase_id].vel->vals[iprev][3 * cell_id + j])
            * extra_i[phase_id].grad_vel[cell_id][id][j];
        }
      }
    }
  }
  else if (cs_glob_lagr_model->cs_used == 0) {
    for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

      unsigned char *particle = p_set->p_buffer + p_am->extents * ip;

      cs_lnum_t cell_id = cs_lagr_particle_get_lnum(particle, p_am,
                                                    CS_LAGR_CELL_ID);
      cs_real_t *part_vel_seen
        = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                   CS_LAGR_VELOCITY_SEEN);

      cs_real_t romf = extra->cromf->val[cell_id];
      for (cs_lnum_t id = 0; id < 3; id++){
        /*  Compute: II = ( -grad(P)/Rom(f)+g) */
        piil[ip][id] = - extra_i[phase_id].grad_pr[cell_id][id] / romf + grav[id];

        /* add grad(<Vf>)*(Us-<Vf>) */
        for (int j = 0; j < 3; j++){
          piil[ip][id]
            -= (part_vel_seen[3 * phase_id + j]
                - extra_i[phase_id].vel->vals[iprev][3 * cell_id + j])
            * extra_i[phase_id].grad_vel[cell_id][id][j];
        }
      }
    }
  }

  /* Compute Pii for code_saturne models
     ----------- */
  else {
    /* Compute: II = ( -grad(P)/Rom(f)+g) */
    for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {
      cs_lnum_t      cell_id  = cs_lagr_particles_get_lnum(p_set, ip ,
                                                           CS_LAGR_CELL_ID);

      cs_real_t romf = extra->cromf->val[cell_id];
      for (int id = 0; id < 3; id++)
        piil[ip][id] = -extra->grad_pr[cell_id][id] / romf + grav[id];
    }

    if (turb_disp_model) {
      /* add grad(<Vf>)*(<Up>-<Us>) if there is enough particles */

      int stat_type = cs_lagr_stat_type_from_attr_id(CS_LAGR_VELOCITY);

      cs_field_t *stat_vel
        = cs_lagr_stat_get_moment(stat_type,
                                  CS_LAGR_STAT_GROUP_PARTICLE,
                                  CS_LAGR_MOMENT_MEAN,
                                  0,
                                  -1);

      stat_type = cs_lagr_stat_type_from_attr_id(CS_LAGR_VELOCITY_SEEN);

      cs_field_t *stat_vel_s
        = cs_lagr_stat_get_moment(stat_type,
                                  CS_LAGR_STAT_GROUP_PARTICLE,
                                  CS_LAGR_MOMENT_MEAN,
                                  0,
                                  -1);

      cs_field_t *stat_w = cs_lagr_stat_get_stat_weight(0);

      for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {
        cs_lnum_t cell_id  = cs_lagr_particles_get_lnum(p_set, ip ,
                                                        CS_LAGR_CELL_ID);

        if (stat_w->val[cell_id] > cs_glob_lagr_stat_options->threshold) {
          for (int i = 0; i < 3; i++) {
            for (cs_lnum_t j = 0; j < 3; j++) {
              cs_real_t vpm   = stat_vel->val[cell_id*3 + j];
              cs_real_t vsm   = stat_vel_s->val[cell_id*3 + j];
              piil[ip][i] += gradvf[cell_id][i][j] * (vpm - vsm);
            }
          }
        }
      }
    }
  }

  /* Add buoyancy effects based on Boussinesq approximation */
  if (    (    cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_COAL
            || cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_CTWR
            || (   cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_HEAT
                && cs_glob_lagr_specific_physics->solve_temperature_seen == 1))
       && cs_field_by_name_try("thermal_expansion") != nullptr
       && cs_field_by_name("thermal_expansion") != nullptr
       && cs_glob_velocity_pressure_model->idilat == 0
       && cs_glob_lagr_model->cs_used == 1) {

    const cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();
    for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {
      cs_lnum_t      cell_id  = cs_lagr_particles_get_lnum(p_set, ip ,
                                                           CS_LAGR_CELL_ID);
      cs_real_t temp_ref = phys_pro->t0;
      cs_real_t temp_s   =
        cs_lagr_particles_get_real(p_set, ip, CS_LAGR_FLUID_TEMPERATURE);

      cs_real_t expansion_coef
        = cs_field_by_name("thermal_expansion")->val[cell_id];

      cs_real_t  _tkelvi = cs_physical_constants_celsius_to_kelvin;
      if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] > -1) {
        /* potential temp at ref */
        cs_real_t pref = cs_glob_atmo_constants->ps;
        cs_real_t rair = phys_pro->r_pg_cnst;
        cs_real_t cp0 = phys_pro->cp0;
        cs_real_t rscp = rair/cp0;
        temp_ref = cs_glob_atmo_option->meteo_t0 *
          pow(pref/ cs_glob_atmo_option->meteo_psea, rscp) - _tkelvi;
      }
      else if (cs_glob_thermal_model->itpscl == CS_TEMPERATURE_SCALE_KELVIN)
          temp_ref -= _tkelvi;

      cs_real_t buoyancy_fac = - expansion_coef * (temp_s - temp_ref);

      for (int id = 0; id < 3; id++)
        piil[ip][id] += buoyancy_fac * grav[id];

      if (   cs_glob_lagr_time_scheme->interpol_field > 0
          && extra->grad_tempf != nullptr) {
        /* Interpolate the local hydrostatic pressure gradient so its is in
         * equillibrium with the interpolated temperature at the position of the
         * particle and not in the center of the cell */
        cs_real_t *part_coord =
          cs_lagr_particles_attr_get_ptr<cs_real_t>(p_set, ip, CS_LAGR_COORDS);
        cs_real_t *cell_cen = cs_glob_mesh_quantities->cell_cen[cell_id];
        for (int id = 0; id < 3; id++) {
          for (int i = 0; i < 3; i++)
            piil[ip][id] += expansion_coef * grav[id]
              * extra->grad_tempf[cell_id][i] * (part_coord[i] - cell_cen[i]);
        }
      }
    }
  }

  /* Add particle back effect seen by mean fluid velocity
   * (two way coupling terms)
   * ==================================================== */
  if (   cs_glob_lagr_source_terms->ltsdyn == 1
      && cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING
      && cs_glob_lagr_model->cs_used == 1) {

    const cs_real_3_t *lagr_st_vel
      = (const cs_real_3_t *)cs_field_by_name("lagr_st_velocity")->val;

    for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {
      cs_lnum_t      cell_id  = cs_lagr_particles_get_lnum(p_set, ip ,
                                                           CS_LAGR_CELL_ID);
      /* Add: II =  -alpha_p rho_p <(U_s -U_p)/tau_p> / (alpha_f rho_f)  */

      cs_real_t romf = extra->cromf->val[cell_id];

      for (int i = 0; i < 3; i++)
        piil[ip][i] += lagr_st_vel[cell_id][i] / romf;

    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
