/*============================================================================
 * Compute particle characteristics: Tp, TL and PI
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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

#include "bft_printf.h"
#include "bft_error.h"
#include "bft_mem.h"

#include "fvm_periodicity.h"

#include "cs_base.h"
#include "cs_defs.h"
#include "cs_math.h"
#include "cs_halo.h"
#include "cs_interface.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_search.h"
#include "cs_timer_stats.h"
#include "cs_thermal_model.h"

#include "cs_field.h"
#include "cs_field_pointer.h"

#include "cs_physical_constants.h"
#include "cs_physical_model.h"

#include "cs_lagr.h"
#include "cs_lagr_new.h"
#include "cs_lagr_particle.h"
#include "cs_lagr_stat.h"
#include "cs_lagr_precipitation_model.h"
#include "cs_lagr_prototypes.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_car.h"

/*----------------------------------------------------------------------------*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute particle characteristics: Tp, TL and PI
 *
 * \param[in] iprev     time step indicator for fields
 *                        0: use fields at current time step
 *                        1: use fields at previous time step
 * \param[in]  dt       time step (per cell)
 * \param[out] taup     dynamic characteristic time
 * \param[out] tlag     fluid characteristic time
 * \param[out] piil     term in integration of up sde
 * \param[out] bx       turbulence characteristics
 * \param[out] tempct   thermal charactersitic time
 * \param[in]  gradpr   pressure gradient
 * \param[in]  gradvf   fluid velocity gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_car(int              iprev,
            const cs_real_t  dt[],
            cs_real_t        taup[],
            cs_real_3_t      tlag[],
            cs_real_3_t      piil[],
            cs_real_33_t     bx[],
            cs_real_t        tempct[],
            cs_real_3_t      gradpr[],
            cs_real_33_t     gradvf[])
{
  /* Particles management */

  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;

  /* Mesh */

  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_lnum_t ncel = mesh->n_cells;

  cs_lagr_extra_module_t *extra = cs_get_lagr_extra_module();
  iprev = CS_MIN(iprev, extra->vel->n_time_vals -1);

  /* Initialization
     ---------------*/

  bool turb_disp_model = false;
  if (   cs_glob_lagr_model->modcpl > 0
      && cs_glob_time_step->nt_cur > cs_glob_lagr_model->modcpl
      && cs_glob_time_step->nt_cur > cs_glob_lagr_stat_options->idstnt)
    turb_disp_model = true;

  cs_lnum_t nor = cs_glob_lagr_time_step->nor;

  cs_real_t bbi[3] = {0.0, 0.0, 0.0};
  cs_real_t ktil   = 0.0;

  cs_real_t rec    = 1000.0;
  cs_real_t c0     = 2.1;
  cs_real_t cl     = 1.0 / (0.5 + (3.0 / 4.0) * c0);
  cs_real_t cb     = 0.8;
  cs_real_t cbcb   = 0.64;
  cs_real_t d6spi  = 6.0 / cs_math_pi;
  cs_real_t d1s3   = 1.0 / 3.0;

  const cs_real_t *grav = cs_glob_physical_constants->gravity;

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

    cs_real_t  rom           = extra->cromf->val[cell_id];
    cs_real_t  xnul          = extra->viscl->val[cell_id] / rom;
    cs_real_t *part_vel_seen = cs_lagr_particle_attr(particle, p_am,
                                                     CS_LAGR_VELOCITY_SEEN);
    cs_real_t *part_vel      = cs_lagr_particle_attr(particle, p_am,
                                                     CS_LAGR_VELOCITY);

    cs_real_t rel_vel_norm = cs_math_3_distance(part_vel_seen, part_vel);

    /* Compute the local Reynolds number */
    cs_real_t rep = rel_vel_norm * p_diam / xnul; /* local Reynolds number */

    cs_real_t d2 = cs_math_sq(p_diam); /* drag coefficient */

    cs_real_t fdr;
    if (rep <= rec)
      fdr = 18.0 * xnul * (1.0 + 0.15 * pow (rep, 0.687)) / d2;
    else
      fdr = 0.44 * 3.0 / 4.0 * rel_vel_norm / p_diam;

    /* Tp computation */
    taup[ip] = p_rom / rom / fdr;

    /* Added-mass term? */
    if (cs_glob_lagr_time_scheme->iadded_mass == 1)
      taup[ip] *= (  1.0
                   + 0.5 * cs_glob_lagr_time_scheme->added_mass_const
                         * rom / p_rom);

    /* Tp user computation */

    cs_user_lagr_rt(ip, rep, rel_vel_norm, rom, p_rom, xnul, taup, dt);

    /* Tc computation */

    if (   (   cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_HEAT
            && cs_glob_lagr_specific_physics->itpvar == 1)
        || (cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_COAL) ) {

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
      if (   cs_glob_physical_model_flag[CS_COMBUSTION_EBU] == 0
          || cs_glob_physical_model_flag[CS_COMBUSTION_EBU] == 2)
        xrkl = extra->diftl0 / rom;
      else if (extra->cpro_viscls != NULL)
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

    cs_real_t  *energi = NULL, *dissip = NULL;
    BFT_MALLOC(energi, ncel, cs_real_t);
    BFT_MALLOC(dissip, ncel, cs_real_t);

    if (extra->itytur == 2 || extra->itytur == 4 || extra->iturb == 50) {

      for (cs_lnum_t cell_id = 0; cell_id < ncel; cell_id++) {
        energi[cell_id] = extra->cvar_k->vals[iprev][cell_id];
        dissip[cell_id] = extra->cvar_ep->vals[iprev][cell_id];
      }

    }
    else if (extra->itytur == 3) {

      if (extra->cvar_rij == NULL) {
        /* Deprecated irijco = 0 */
        for (cs_lnum_t cell_id = 0; cell_id < ncel; cell_id++) {

          energi[cell_id] = 0.5 * (  extra->cvar_r11->vals[iprev][cell_id]
                                   + extra->cvar_r22->vals[iprev][cell_id]
                                   + extra->cvar_r33->vals[iprev][cell_id]);
          dissip[cell_id] = extra->cvar_ep->vals[iprev][cell_id];

        }
      }
      else {
        /* irijco = 1 */
        for (cs_lnum_t cell_id = 0; cell_id < ncel; cell_id++) {

          energi[cell_id] = 0.5 * ( extra->cvar_rij->vals[iprev][6*cell_id]
                                  + extra->cvar_rij->vals[iprev][6*cell_id + 1]
                                  + extra->cvar_rij->vals[iprev][6*cell_id + 2]
                                  );
          dissip[cell_id] = extra->cvar_ep->vals[iprev][cell_id];
        }
      }
    }
    else if (extra->iturb == 60) {

      for (cs_lnum_t cell_id = 0; cell_id < ncel; cell_id++) {
        energi[cell_id] = extra->cvar_k->vals[iprev][cell_id];
        dissip[cell_id] = extra->cmu * energi[cell_id]
                                     * extra->cvar_omg->vals[iprev][cell_id];
      }

    }
    else {

      bft_error
        (__FILE__, __LINE__, 0,
         _("Lagrangian turbulent dispersion is not compatible with\n"
           "the selected turbulence model.\n"
           "\n"
           "Turbulent dispersion is taken into account with idistu = %d\n"
           " Activated turbulence model is %d, when only k-eps, LES, Rij-eps,\n"
           " V2f or k-omega are handled."),
         (int)cs_glob_lagr_model->idistu,
         (int)extra->iturb);

    }

    /* -> Calculation of TL and BX     */

    for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

      if (cs_lagr_particles_get_flag(p_set, ip, CS_LAGR_PART_FIXED)) {
        for (cs_lnum_t id = 0; id < 3; id++ ) {
          tlag[ip][id]      = cs_math_epzero;
          bx[ip][id][nor-1] = 0.0;
        }
      }

      unsigned char *particle = p_set->p_buffer + p_am->extents * ip;

      cs_lnum_t cell_id = cs_lagr_particle_get_lnum(particle, p_am,
                                                    CS_LAGR_CELL_ID);

      cs_real_t mean_part_vel[3], fluid_vel[3];

      if (dissip[cell_id] > 0.0 && energi[cell_id] > 0.0) {

        cs_real_t *part_vel_seen = cs_lagr_particle_attr(particle, p_am,
                                                         CS_LAGR_VELOCITY_SEEN);
        cs_real_t *part_vel      = cs_lagr_particle_attr(particle, p_am,
                                                         CS_LAGR_VELOCITY);

        cs_real_t tl  = cl * energi[cell_id] / dissip[cell_id];
        tl  = CS_MAX(tl, cs_math_epzero);

        for (cs_lnum_t i = 0; i < 3; i++) {
          mean_part_vel[i] = part_vel[i];
          fluid_vel[i] = part_vel_seen[i];
        }

        /* Compute the main direction in the global reference
         * frame */
        cs_real_t dir[3] = {0, 0, 0};

        if (turb_disp_model) {

          int stat_type = cs_lagr_stat_type_from_attr_id(CS_LAGR_VELOCITY);

          cs_field_t *stat_vel
            = cs_lagr_stat_get_moment(stat_type,
                                      CS_LAGR_STAT_GROUP_PARTICLE,
                                      CS_LAGR_MOMENT_MEAN,
                                      0,
                                      -1);

          cs_field_t *stat_w = cs_lagr_stat_get_stat_weight(0);

          if (stat_w->val[cell_id] > cs_glob_lagr_stat_options->threshold) {

            for (cs_lnum_t i = 0; i < 3; i++) {
              mean_part_vel[i] = stat_vel->val[cell_id * 3 + i];
              fluid_vel[i] = extra->vel->vals[iprev][cell_id * 3 + i];
            }

          }

          for (cs_lnum_t i = 0; i < 3; i++)
            dir[i] = mean_part_vel[i] - fluid_vel[i];
          cs_math_3_normalize(dir, dir);

        }

        cs_real_t uvwdif = 0.;
        for  (cs_lnum_t i = 0; i < 3; i++)
          uvwdif += cs_math_sq(fluid_vel[i] - mean_part_vel[i]);

        uvwdif = (3.0 * uvwdif) / (2.0 * energi[cell_id]);

        if (turb_disp_model) {

          /* Crossing trajectory in the direction of "<u_f>-<u_p>"
           * and in the span-wise direction */
          cs_real_t an, at;
          an = (1.0 + cbcb * uvwdif);
          at = (1.0 + 4.0 * cbcb * uvwdif);

          bbi[0] = sqrt(an); /* First direction, n, in the local reference
                                frame */
          bbi[1] = sqrt(at); /* Second and third direction, orthogonal to n */
          bbi[2] = sqrt(at);

          /* Compute the timescale of the fluid velocities seen by discrete
           * particles in parallel and transverse directions */
          for (cs_lnum_t id = 0; id < 3; id++)
            tlag[ip][id] = tl / bbi[id];

          if (extra->itytur == 3) {

            /* Deprecated irijco = 0 */
            if (extra->cvar_rij == NULL) {
              cs_real_t rij[6] = {
                extra->cvar_r11->vals[iprev][cell_id],
                extra->cvar_r22->vals[iprev][cell_id],
                extra->cvar_r33->vals[iprev][cell_id],
                0,//FIXME extra->cvar_r12->vals[iprev][cell_id],
                0,//extra->cvar_r23->vals[iprev][cell_id],
                0};//extra->cvar_r13->vals[iprev][cell_id]};
              /* Note that n.R.n = R : n(x)n */
              cs_real_t rnn = cs_math_3_sym_33_3_dot_product(dir, rij, dir);
              cs_real_t tr_r = cs_math_6_trace(rij);
              // bbn * R : n(x)n + bbt * R : (1 - n(x)n)
              ktil = 3.0 * (rnn * bbi[0] + (tr_r -rnn) * bbi[1]) /* bbi[1] == bbi[2] is used */
                         / (2.0 * (bbi[0] + bbi[1] + bbi[2]));
            } else {
              cs_real_t *rij = &(extra->cvar_rij->vals[iprev][6*cell_id]);
              /* Note that n.R.n = R : n(x)n */
              cs_real_t rnn = cs_math_3_sym_33_3_dot_product(dir, rij, dir);
              cs_real_t tr_r = cs_math_6_trace(rij);
              // bbn * R : n(x)n + bbt * R : (1 - n(x)n)
              ktil = 3.0 * (rnn * bbi[0] + (tr_r -rnn) * bbi[1]) /* bbi[1] == bbi[2] is used */
                         / (2.0 * (bbi[0] + bbi[1] + bbi[2]));
            }

          }
          else if (   extra->itytur == 2 || extra->itytur == 4
                   || extra->iturb == 50 || extra->iturb == 60) {
            ktil = energi[cell_id];
          }

          for (cs_lnum_t id = 0; id < 3; id++) {

            cs_real_t bxi
              = dissip[cell_id] * (  (c0  * bbi[id] * ktil / energi[cell_id])
                                   + (  (bbi[id] * ktil / energi[cell_id] - 1.0)
                                      * 2.0 / 3.0));
            if (bxi > 0.0)
              bx[ip][id][nor-1] = sqrt(bxi);
            else
              bx[ip][id][nor-1] = 0.0;

          }

        }
        else {

          for (cs_lnum_t id = 0; id < 3; id++)
            tlag[ip][id] = tl;

          if (cs_glob_lagr_model->idiffl == 0) {
            uvwdif      = sqrt(uvwdif);
            tlag[ip][0] = tl / (1.0 + cb * uvwdif);
            tlag[ip][1] = tlag[ip][0];
            tlag[ip][2] = tlag[ip][0];
          }

          bx[ip][0][nor-1] = sqrt (c0 * dissip[cell_id]);
          bx[ip][1][nor-1] = bx[ip][0][nor-1];
          bx[ip][2][nor-1] = bx[ip][0][nor-1];

        }

      }

    }

    BFT_FREE(energi);
    BFT_FREE(dissip);

  }
  else {

    for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

      /* FIXME we may still need to do computations here */
      if (cs_lagr_particles_get_flag(p_set, ip, CS_LAGR_PART_FIXED))
        continue;

      for (cs_lnum_t id = 0; id < 3; id++ ) {
        tlag[ip][id] = cs_math_epzero;
        bx[ip][id][nor-1] = 0.0;
      }

    }

  }

  /* Compute Pii
     ----------- */

  if (turb_disp_model) {

    int stat_type = cs_lagr_stat_type_from_attr_id(CS_LAGR_VELOCITY);

    cs_field_t *stat_vel
      = cs_lagr_stat_get_moment(stat_type,
                                CS_LAGR_STAT_GROUP_PARTICLE,
                                CS_LAGR_MOMENT_MEAN,
                                0,
                                -1);

    cs_field_t *stat_w = cs_lagr_stat_get_stat_weight(0);

    for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

      unsigned char *particle = p_set->p_buffer + p_am->extents * ip;

      cs_lnum_t cell_id = cs_lagr_particle_get_lnum(particle, p_am,
                                                    CS_LAGR_CELL_ID);

      /* Compute: II = ( -grad(P)/Rom(f)+grad(<Vf>)*(<Up>-<Uf>) + g ) */

      cs_real_t romf = extra->cromf->val[cell_id];

      for (cs_lnum_t id = 0; id < 3; id++) {

        piil[ip][id] = -gradpr[cell_id][id] / romf + grav[id];

        if (stat_w->val[cell_id] > cs_glob_lagr_stat_options->threshold) {

          for (cs_lnum_t i = 0; i < 3; i++) {
            cs_real_t vpm   = stat_vel->val[cell_id*3 + i];
            cs_real_t fluid_vel = extra->vel->vals[iprev][cell_id*3 + i];
            piil[ip][id] += gradvf[cell_id][id][i] * (vpm - fluid_vel);
          }

        }

      }
    }
  }
  else {

    for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

      unsigned char *particle = p_set->p_buffer + p_am->extents * ip;

      cs_lnum_t cell_id = cs_lagr_particle_get_lnum(particle, p_am,
                                                    CS_LAGR_CELL_ID);

      /* Compute: II = ( -grad(P)/Rom(f) + g) */

      cs_real_t romf = extra->cromf->val[cell_id];

      for (cs_lnum_t id = 0; id < 3; id++)
        piil[ip][id] = -gradpr[cell_id][id] / romf + grav[id];

    }

  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
