/*============================================================================
 * Compute particle characteristics: Tp, TL and PI
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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
#include "cs_prototypes.h"
#include "cs_search.h"
#include "cs_timer_stats.h"
#include "cs_thermal_model.h"

#include "cs_field.h"
#include "cs_field_pointer.h"

#include "cs_gui_particles.h"

#include "cs_physical_constants.h"
#include "cs_physical_model.h"

#include "cs_lagr.h"
#include "cs_lagr_geom.h"
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
 * \param[in]  energi   work array
 * \param[in]  dissip   work array
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
            cs_real_33_t     gradvf[],
            cs_real_t        energi[],
            cs_real_t        dissip[])
{
  /* Particles management */

  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;

  /* Mesh */

  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_lnum_t ncel = mesh->n_cells;

  cs_lagr_extra_module_t *extra = cs_get_lagr_extra_module();

  /* Initialization
     ---------------*/

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
    cs_lnum_t      cell_id  = cs_lagr_particle_get_cell_id(particle, p_am);

    if (cell_id >= 0) {

      cs_real_t p_mass = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS);
      cs_real_t p_rom  = p_mass * d6spi / pow(p_diam, 3.0);

      cs_real_t  rom           = extra->cromf->val[cell_id];
      cs_real_t  xnul          = extra->viscl->val[cell_id] / rom;
      cs_real_t *part_vel_seen = cs_lagr_particle_attr(particle, p_am,
                                                       CS_LAGR_VELOCITY_SEEN);
      cs_real_t *part_vel      = cs_lagr_particle_attr(particle, p_am, CS_LAGR_VELOCITY);

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

      /**************************
       *    Tc computation      *
       **************************/
      if (   (   cs_glob_lagr_model->physical_model == 1
              && cs_glob_lagr_specific_physics->itpvar == 1)
          || (cs_glob_lagr_model->physical_model == 2) ) {

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

  }

  /* ====================================================================
   * 3. CALCUL DE TL
   * ====================================================================   */

  /* Compute turbulent kinetic energy and dissipation
     based on turbulence model */

  if (cs_glob_lagr_time_scheme->idistu == 1) {

    if (extra->itytur == 2 || extra->iturb == 50) {

      for (cs_lnum_t cell_id = 0; cell_id < ncel; cell_id++) {

        energi[cell_id] = extra->cvar_k->vals[iprev][cell_id];
        dissip[cell_id] = extra->cvar_ep->vals[iprev][cell_id];

      }

    }

    else if (extra->itytur == 3) {

      if (extra->cvar_rij == NULL) {
        /* irijco = 0 */

        for (cs_lnum_t cell_id = 0; cell_id < ncel; cell_id++) {

          energi[cell_id] = 0.5 * (  extra->cvar_r11->vals[iprev][cell_id]
                                   + extra->cvar_r22->vals[iprev][cell_id]
                                   + extra->cvar_r33->vals[iprev][cell_id]);
          dissip[cell_id] = extra->cvar_ep->vals[iprev][cell_id];

        }
      } else {
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

      bft_printf(_("\n WARNING: STOP AT LAGRANGIAN MODULE EXECUTION\n"));
      bft_printf
        (_("The lagrangian module is not compatible with the selected turbulence model.\n"
           "\n"
           "Turbulent dispersion is taken into account with IDISTU = %d\n"
           " Activated turbulence model is %d, when only k-eps, Rij-eps,\n"
           " V2f or k-omega are compatible with turbulent dispersion and Lagrangian module.\n"
           "\n"),
         (int)cs_glob_lagr_time_scheme->idistu,
         (int)extra->iturb);
      cs_exit(1);

    }

    /* -> Calcul de TL et BX     */

    for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

      unsigned char *particle = p_set->p_buffer + p_am->extents * ip;

      cs_lnum_t cell_id = cs_lagr_particle_get_cell_id(particle, p_am);

      if (cell_id >= 0) {

        cs_real_t vpart[3], vflui[3];

        if (dissip[cell_id] > 0.0 && energi[cell_id] > 0.0) {

          cs_real_t *part_vel_seen = cs_lagr_particle_attr(particle, p_am,
                                                           CS_LAGR_VELOCITY_SEEN);
          cs_real_t *part_vel      = cs_lagr_particle_attr(particle, p_am,
                                                           CS_LAGR_VELOCITY);

          cs_real_t tl  = cl * energi[cell_id] / dissip[cell_id];
          tl  = CS_MAX(tl, cs_math_epzero);

          for (cs_lnum_t i = 0; i < 3; i++){

            vpart[i] = part_vel[i];
            vflui[i] = part_vel_seen[i];

          }

          if (   cs_glob_lagr_time_scheme->modcpl > 0
              && cs_glob_time_step->nt_cur > cs_glob_lagr_time_scheme->modcpl) {

            int stat_type = cs_lagr_stat_type_from_attr_id(CS_LAGR_VELOCITY);

            cs_field_t *stat_vel
              = cs_lagr_stat_get_moment(stat_type,
                                        CS_LAGR_MOMENT_MEAN,
                                        0,
                                        -1);

            cs_field_t *stat_w = cs_lagr_stat_get_stat_weight(0);

            if (stat_w->val[cell_id] > cs_glob_lagr_stat_options->threshold) {

              for (cs_lnum_t i = 0; i < 3; i++){
                vpart[i] = stat_vel->val[cell_id * 3 + i];
                vflui[i] = extra->vel->vals[iprev][cell_id * 3 + i];
              }

            }

          }

          cs_real_t uvwdif = 0.;
          for  (cs_lnum_t i = 0; i < 3; i++)
            uvwdif += cs_math_sq(vflui[i] - vpart[i]);

          uvwdif   = (3.0 * uvwdif) / (2.0 * energi[cell_id]);

          if (   cs_glob_lagr_time_scheme->modcpl > 0
              && cs_glob_time_step->nt_cur > cs_glob_lagr_time_scheme->modcpl) {

            if (cs_glob_lagr_time_scheme->idirla == 1) {

              bbi[0] = sqrt (1.0 + cbcb * uvwdif);
              bbi[1] = sqrt (1.0 + 4.0 * cbcb * uvwdif);
              bbi[2] = sqrt (1.0 + 4.0 * cbcb * uvwdif);
              for (cs_lnum_t id = 0; id < 3; id++)
                tlag[ip][id] = tl / bbi[id];

            }

            else if (cs_glob_lagr_time_scheme->idirla == 2) {

              bbi[0] = sqrt (1.0 + 4.0 * cbcb * uvwdif);
              bbi[1] = sqrt (1.0 + cbcb * uvwdif);
              bbi[2] = sqrt (1.0 + 4.0 * cbcb * uvwdif);
              for (cs_lnum_t id = 0; id < 3; id++)
                tlag[ip][id] = tl / bbi[id];

            }

            else if (cs_glob_lagr_time_scheme->idirla == 3) {

              bbi[0] = sqrt (1.0 + 4.0 * cbcb * uvwdif);
              bbi[1] = sqrt (1.0 + 4.0 * cbcb * uvwdif);
              bbi[2] = sqrt (1.0 + cbcb * uvwdif);
              for (cs_lnum_t id = 0; id < 3; id++)
                tlag[ip][id] = tl / bbi[id];

            }

            else {

              bft_printf(_("\n WARNING: STOP AT LAGRANGIAN MODULE EXECUTION\n"));
              bft_printf
                (_("Complete model direction choice value is not correct (in LAGCAR).\n"
                   " IDIRLA should be an integer ranging from 1 to 3:\n"
                   " =1 for a flow along x-axis\n"
                   " =2 for a flow along y-axis\n"
                   " =3 for a flow along z-axis\n"
                   " Here IDIRLA = %d\n"),
                 (int)cs_glob_lagr_time_scheme->idirla);
              cs_exit(1);

            }

            if (extra->itytur == 3) {

              if (extra->cvar_rij == NULL) {
                cs_real_t r11  = extra->cvar_r11->vals[iprev][cell_id];
                cs_real_t r22  = extra->cvar_r22->vals[iprev][cell_id];
                cs_real_t r33  = extra->cvar_r33->vals[iprev][cell_id];
                ktil = 3.0 * (r11 * bbi[0] + r22 * bbi[1] + r33 * bbi[2])
                           / (2.0 * (bbi[0] + bbi[1] + bbi[2]));
              } else {
                cs_real_t r11  = extra->cvar_rij->vals[iprev][6*cell_id    ];
                cs_real_t r22  = extra->cvar_rij->vals[iprev][6*cell_id + 1];
                cs_real_t r33  = extra->cvar_rij->vals[iprev][6*cell_id + 2];

                ktil = 3.0 * (r11 * bbi[0] + r22 * bbi[1] + r33 * bbi[2])
                           / (2.0 * (bbi[0] + bbi[1] + bbi[2]));
              }

            }

            else if (extra->itytur == 2 || extra->iturb == 50 || extra->iturb == 60)
              ktil = energi[cell_id];

            for (cs_lnum_t id = 0; id < 3; id++) {

              cs_real_t bxi = dissip[cell_id] * (  (c0  * bbi[id] * ktil / energi[cell_id])
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

            if (cs_glob_lagr_time_scheme->idiffl == 0) {

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

        else {

          for (cs_lnum_t id = 0; id < 3; id++ ){

            tlag[ip][id]    = cs_math_epzero;
            bx[ip][id][nor-1] = 0.0;

          }

        }

      }

    }

  }

  else {

    for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

      unsigned char *particle = p_set->p_buffer + p_am->extents * ip;

      if (cs_lagr_particle_get_cell_id(particle, p_am) >= 0) {

        for (cs_lnum_t id = 0; id < 3; id++ ){

          tlag[ip][id] = cs_math_epzero;
          bx[ip][id][nor-1] = 0.0;

        }

      }

    }

  }

  /* ====================================================================   */
  /* 4. CALCUL DE PII     */
  /* ====================================================================   */

  for (cs_lnum_t id = 0; id < 3; id++) {

    for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

      unsigned char *particle = p_set->p_buffer + p_am->extents * ip;

      cs_lnum_t cell_id = cs_lagr_particle_get_cell_id(particle, p_am);

      if (cell_id >= 0) {

        /* -->   Compute: II = ( -grad(P)/Rom(f)+grad(<Vf>)*(<Up>-<Uf>) + g )
         *       or
         *       Compute: II = ( -grad(P)/Rom(f) + g) */

        cs_real_t romf = extra->cromf->val[cell_id];

        piil[ip][id] = -gradpr[cell_id][id] / romf + grav[id];

        if (   cs_glob_lagr_time_scheme->modcpl > 0
            &&   cs_glob_time_step->nt_cur
               > cs_glob_lagr_time_scheme->modcpl) {

          int stat_type = cs_lagr_stat_type_from_attr_id(CS_LAGR_VELOCITY);

          cs_field_t *stat_vel
            = cs_lagr_stat_get_moment(stat_type,
                                      CS_LAGR_MOMENT_MEAN,
                                      0,
                                      -1);

          cs_field_t *stat_w = cs_lagr_stat_get_stat_weight(0);

          if (stat_w->val[cell_id] > cs_glob_lagr_stat_options->threshold) {

            for (cs_lnum_t i = 0; i < 3; i++) {

              cs_real_t vpm   = stat_vel->val[i + cell_id * 3];
              cs_real_t vflui = extra->vel->vals[iprev][i + cell_id * 3];
              piil[ip][id] += gradvf[cell_id][id][i] * (vpm - vflui);

            }

          }

        }

      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
