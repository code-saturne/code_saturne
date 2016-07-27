/*============================================================================
 * Methods for particle resuspension
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
 * Functions dealing with lagrangian resuspension
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

#include "cs_base.h"
#include "cs_math.h"
#include "cs_prototypes.h"

#include "bft_mem.h"
#include "bft_error.h"

#include "cs_physical_constants.h"
#include "cs_thermal_model.h"

#include "cs_lagr.h"
#include "cs_lagr_tracking.h"
#include "cs_lagr_roughness.h"
#include "cs_lagr_adh.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_resuspension.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_lagr_resuspension.c
        Lagrangian resuspension model.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local type definitions
 *============================================================================*/


/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the particle resuspension
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_resuspension(void)
{
  const cs_real_t tkelvi = 273.15;

  cs_lagr_particle_set_t *p_set = cs_lagr_get_particle_set();
  const cs_lagr_attribute_map_t *p_am = p_set->p_am;

  cs_lagr_boundary_interactions_t *lag_bi = cs_glob_lagr_boundary_interactions;
  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  cs_lnum_t n_faces = mesh->n_b_faces;

  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;

  const cs_real_3_t *restrict i_face_normal
        = (const cs_real_3_t *restrict)fvq->i_face_normal;
  const cs_real_3_t *restrict b_face_normal
        = (const cs_real_3_t *restrict)fvq->b_face_normal;

/* ================================================================  */
/* 1.    Resuspension sub model   */
/* ================================================================  */

  cs_lnum_t test_colli;
  cs_real_t adhesion_energ;

  for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

    unsigned char *part = p_set->p_buffer + p_am->extents * ip;
    cs_lnum_t face_id       = cs_lagr_particle_get_lnum(part, p_am, CS_LAGR_NEIGHBOR_FACE_ID);
    cs_real_t p_mass        = cs_lagr_particle_get_real(part, p_am, CS_LAGR_MASS);
    cs_real_t p_stat_weight = cs_lagr_particle_get_real(part, p_am, CS_LAGR_STAT_WEIGHT);
    cs_real_t p_diam        = cs_lagr_particle_get_real(part, p_am, CS_LAGR_DIAMETER);

    cs_real_t *part_vel     = cs_lagr_particle_attr(part, p_am, CS_LAGR_VELOCITY);
    cs_real_t *prev_part_vel= cs_lagr_particle_attr_n(part, p_am, 1, CS_LAGR_VELOCITY);

    test_colli = 0;

    cs_lnum_t iel = cs_lagr_particle_get_cell_id(part, p_am);

    cs_real_t temp;

    if (extra->scal_t != NULL) {

      if (   cs_glob_thermal_model->itherm == 1
          && cs_glob_thermal_model->itpscl == 2)
        temp = extra->scal_t->val[iel] + tkelvi;

      else if (   cs_glob_thermal_model->itherm == 1
               && cs_glob_thermal_model->itpscl == 1)
        temp = extra->scal_t->val[iel];

      else if (cs_glob_thermal_model->itherm == 2){

        cs_lnum_t mode = 1;
        CS_PROCF (usthht,USTHHT)(&mode, &(extra->scal_t->val[iel]), &temp);

      }

    }
    else
      temp = cs_glob_fluid_properties->t0;

    cs_lnum_t flag = cs_lagr_particle_get_lnum(part, p_am, CS_LAGR_DEPOSITION_FLAG);
    if (flag == CS_LAGR_PART_DEPOSITED)
      /* The particle has just deposited     */
      /* The adhesion force is calculated    */
      cs_lagr_adh(ip, temp, &adhesion_energ);

    else if (flag == CS_LAGR_PART_ROLLING) {

      /* The particle is rolling   */

      /* if the number of great asperities   */
      /* is null it is marked for a possible collision */
      if (cs_lagr_particle_get_lnum(part, p_am, CS_LAGR_N_LARGE_ASPERITIES) == 0)
        test_colli = 1;

      cs_real_t disp_norm = cs_lagr_particle_get_real(part, p_am, CS_LAGR_DISPLACEMENT_NORM);

      if (disp_norm > p_diam && disp_norm < 2.0 * p_diam) {

        /* If the particle has a displacement approximately   *
         * equal to a diameter, recalculation of the adhesion force     */

        cs_lagr_particle_set_real(part, p_am, CS_LAGR_DISPLACEMENT_NORM, 0.0);

        cs_lagr_adh(ip, temp, &adhesion_energ);

        if (test_colli == 1
            && cs_lagr_particle_get_lnum(part, p_am, CS_LAGR_N_LARGE_ASPERITIES) > 0) {

          cs_real_t kinetic_energy =  0.5 * p_mass
                                   * cs_math_3_dot_product(part_vel, part_vel);

          if (kinetic_energy > adhesion_energ) {

            /* The particle is resuspended and its kinetic energy is totally converted
             * along the wall-normal distance                                           */


            cs_lagr_particle_set_lnum(part, p_am, CS_LAGR_DEPOSITION_FLAG,
                                      CS_LAGR_PART_IN_FLOW);
            cs_lagr_particle_set_real(part, p_am, CS_LAGR_ADHESION_FORCE, 0.0);
            cs_lagr_particle_set_real(part, p_am, CS_LAGR_ADHESION_TORQUE, 0.0);
            cs_lagr_particle_set_lnum(part, p_am, CS_LAGR_N_LARGE_ASPERITIES, 0);
            cs_lagr_particle_set_lnum(part, p_am, CS_LAGR_N_SMALL_ASPERITIES, 0);
            cs_lagr_particle_set_real(part, p_am, CS_LAGR_DISPLACEMENT_NORM, 0.0);

            cs_real_t norm_face = fvq->b_face_surf[face_id];

            cs_real_t norm_velocity = sqrt(cs_math_3_dot_product(part_vel,
                                                                 part_vel));

            for (int id = 0; id < 3; id++)
              part_vel[id] = -norm_velocity / norm_face
                            * b_face_normal[face_id][id];

            /* Update of the number and weight of resuspended particles     */
            p_set->n_part_resusp += 1;
            p_set->weight_resusp += p_stat_weight;

            bound_stat[lag_bi->ires   * n_faces + face_id] +=   p_stat_weight;
            bound_stat[lag_bi->iflres * n_faces + face_id] +=   p_stat_weight * p_mass / norm_face;
            bound_stat[lag_bi->iflm   * n_faces + face_id] += - p_stat_weight * p_mass / norm_face;

          }

        }

      }
      else if (disp_norm >= 2.0 * p_diam) {

        cs_lnum_t ndiam = (cs_lnum_t)(disp_norm / p_diam);

        cs_lnum_t ii = 1;
        while (   ii <= ndiam
               && (cs_lagr_particle_get_lnum(part, p_am, CS_LAGR_DEPOSITION_FLAG)
                 == CS_LAGR_PART_ROLLING)) {

          cs_lagr_adh(ip, temp, &adhesion_energ);

          /* Reconstruct an estimate of the particle velocity   */
          /* at the current sub-time-step assuming linear variation  */
          /* (constant acceleration)   */

          cs_real_t v_part_t = sqrt(cs_math_3_dot_product(prev_part_vel,
                                                          prev_part_vel));
          cs_real_t v_part_t_dt = sqrt(cs_math_3_dot_product(part_vel,
                                                             part_vel));

          cs_real_t sub_dt = cs_glob_lagr_time_step->dtp / ndiam;

          cs_real_t v_part_inst = v_part_t + sub_dt * (v_part_t_dt + v_part_t) / cs_glob_lagr_time_step->dtp;

          /* Reconstruct an estimate of the angular velocity    */
          /* at the current sub-time-step   */

          cs_real_t omep = v_part_inst / (p_diam * 0.5);

          /* Variation of the angular velocity due to */
          /* the update of the adhesion torque   */

          cs_real_t domep =  cs_lagr_particle_get_real(part, p_am, CS_LAGR_ADHESION_TORQUE)
                           / (7.0 / 5.0 * p_mass * pow(p_diam * 0.5, 2));

          if ((domep * sub_dt) > omep) {

            cs_lagr_particle_set_lnum(part, p_am, CS_LAGR_DEPOSITION_FLAG,
                                      CS_LAGR_PART_NO_MOTION);
            part_vel[0] = 0.0;
            part_vel[1] = 0.0;
            part_vel[2] = 0.0;

          }

          if (   test_colli == 1
              && cs_lagr_particle_get_lnum(part, p_am, CS_LAGR_N_LARGE_ASPERITIES) > 0) {

            cs_real_t kinetic_energy =  0.5 * p_mass
                                      * (  pow(part_vel[0], 2)
                                         + pow(part_vel[1], 2)
                                         + pow(part_vel[2], 2));

            if (kinetic_energy > adhesion_energ) {

              /* The particle is resuspended    */
              /* and its kinetic energy is totally converted   */
              /* along the wall-normal distance */
              cs_lagr_particle_set_lnum(part, p_am, CS_LAGR_DEPOSITION_FLAG,
                                        CS_LAGR_PART_IN_FLOW);
              cs_lagr_particle_set_real(part, p_am, CS_LAGR_ADHESION_FORCE, 0.0);
              cs_lagr_particle_set_real(part, p_am, CS_LAGR_ADHESION_TORQUE, 0.0);
              cs_lagr_particle_set_lnum(part, p_am, CS_LAGR_N_LARGE_ASPERITIES, 0);
              cs_lagr_particle_set_lnum(part, p_am, CS_LAGR_N_SMALL_ASPERITIES, 0);
              cs_lagr_particle_set_real(part, p_am, CS_LAGR_DISPLACEMENT_NORM, 0.0);

              cs_real_t norm_face = fvq->b_face_surf[face_id];

              cs_real_t norm_velocity = sqrt(cs_math_3_dot_product(part_vel,
                                                                   part_vel));

              for (int id = 0; id < 3; id++)
                part_vel[id] = -norm_velocity / norm_face
                              * b_face_normal[face_id][id];

              /* Update of the number and weight of resuspended particles     */
              p_set->n_part_resusp += 1;
              p_set->weight_resusp += p_stat_weight;

              bound_stat[lag_bi->ires   * n_faces + face_id] +=   p_stat_weight;
              bound_stat[lag_bi->iflres * n_faces + face_id] +=   p_stat_weight * p_mass / norm_face;
              bound_stat[lag_bi->iflm   * n_faces + face_id] += - p_stat_weight * p_mass / norm_face;

            }

            if (cs_lagr_particle_get_lnum(part, p_am, CS_LAGR_N_LARGE_ASPERITIES) == 0)
              test_colli = 1;

          }

          ii++;

        }

      }

    }
    /* Treatment of internal deposition and user imposed motion */
    else if (flag == CS_LAGR_PART_IMPOSED_MOTION) {

      /* Reorient the face so that it is the outwarding normal */
      int reorient_face = 1;
      if (iel == mesh->i_face_cells[face_id][1])
        reorient_face = -1;

      /* Adhesion forces not implemented */

      /* Gravity forces */
      cs_real_3_t gravity = {cs_glob_physical_constants->gx,
                             cs_glob_physical_constants->gy,
                             cs_glob_physical_constants->gz};

      cs_real_t fgrav
        =   p_mass * reorient_face
          * cs_math_3_dot_product(gravity, i_face_normal[face_id]);

      /* Forces due to pressure difference */
      cs_lnum_t c_id1 = mesh->i_face_cells[face_id][0];
      cs_lnum_t c_id2 = mesh->i_face_cells[face_id][1];

      if (iel == c_id2) {
        c_id2 = c_id1;
        c_id1 = iel;
      }

      /* Warning: dynamic pressure, FIXME for iphydr = 1 */
      cs_real_t press_out = cs_glob_lagr_extra_module->pressure->val[c_id1];
      cs_real_t press_in = cs_glob_lagr_extra_module->pressure->val[c_id2];

      cs_real_t fpres = (press_out - press_in) * cs_math_pi * pow(p_diam, 2) * 0.25
        * cs_lagr_particle_get_real(part, p_am, CS_LAGR_FOULING_INDEX);

      /* Resuspension criterion: Fgrav + Fpres < 0 */
      if ((fgrav + fpres) < 0.) {
        cs_lagr_particle_set_lnum(part, p_am, CS_LAGR_DEPOSITION_FLAG,
                                  CS_LAGR_PART_IN_FLOW);
        /* TODO: impose particle velocity? Do some stats? */
      }
    }
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
