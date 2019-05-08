/*============================================================================
 * Methods for particle resuspension
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

#include "bft_mem.h"
#include "bft_error.h"

#include "cs_physical_constants.h"
#include "cs_prototypes.h"
#include "cs_random.h"
#include "cs_thermal_model.h"

#include "cs_lagr.h"
#include "cs_lagr_stat.h"
#include "cs_lagr_tracking.h"
#include "cs_lagr_roughness.h"
#include "cs_lagr_adh.h"
#include "cs_lagr_event.h"

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
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a resulspension event
 *
 * TODO add additional info to events.
 *
 * \param[in]  events             pointer to events set
 * \param[in]  particles          pointer to particle set
 * \param[in]  p_id               particle id
 * \param[in]  face_id            associated face id
 * \param[in]  particle_velocity  velocity after event
 */
/*----------------------------------------------------------------------------*/

static void
_add_resuspension_event(cs_lagr_event_set_t     *events,
                        cs_lagr_particle_set_t  *particles,
                        cs_lnum_t                p_id,
                        cs_lnum_t                face_id,
                        cs_real_t                particle_velocity[3])
{
  cs_lnum_t event_id = events->n_events;
  if (event_id >= events->n_events_max) {
    /* flush events */
    cs_lagr_stat_update_event(events,
                              CS_LAGR_STAT_GROUP_TRACKING_EVENT);
    events->n_events = 0;
    event_id = 0;
  }

  cs_lagr_event_init_from_particle(events, particles, event_id, p_id);

  cs_lagr_events_set_lnum(events,
                          event_id,
                          CS_LAGR_E_FACE_ID,
                          face_id);

  cs_lnum_t *e_flag = cs_lagr_events_attr(events,
                                          event_id,
                                          CS_LAGR_E_FLAG);

  cs_real_t *e_vel_post = cs_lagr_events_attr(events,
                                              event_id,
                                              CS_LAGR_E_VELOCITY);

  *e_flag = *e_flag | CS_EVENT_RESUSPENSION;

  for (int k = 0; k < 3; k++)
    e_vel_post[k] = particle_velocity[k];

  events->n_events += 1;
}

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
  const cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

  cs_lagr_particle_set_t *p_set = cs_lagr_get_particle_set();
  const cs_lagr_attribute_map_t *p_am = p_set->p_am;

  cs_lagr_boundary_interactions_t *lag_bi = cs_glob_lagr_boundary_interactions;
  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  cs_lnum_t n_faces = mesh->n_b_faces;

  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;

  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)fvq->i_face_normal;

  cs_lagr_event_set_t  *events = NULL;

  if (cs_lagr_stat_is_active(CS_LAGR_STAT_GROUP_TRACKING_EVENT)) {
    events = cs_lagr_event_set_boundary_interaction();
    /* Event set "expected" size: n boundary faces*2 */
    cs_lnum_t events_min_size = mesh->n_b_faces * 2;
    if (events->n_events_max < events_min_size)
      cs_lagr_event_set_resize(events, events_min_size);
  }

  /* ================================================================  */
  /* 1.    Resuspension sub model   */
  /* ================================================================  */

  cs_lnum_t test_colli;
  cs_real_t adhesion_energ;

  for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

    unsigned char *part = p_set->p_buffer + p_am->extents * ip;
    cs_lnum_t face_id
      = cs_lagr_particle_get_lnum(part, p_am, CS_LAGR_NEIGHBOR_FACE_ID);
    cs_real_t p_mass = cs_lagr_particle_get_real(part, p_am, CS_LAGR_MASS);
    cs_real_t p_stat_weight
      = cs_lagr_particle_get_real(part, p_am, CS_LAGR_STAT_WEIGHT);
    cs_real_t p_diam = cs_lagr_particle_get_real(part, p_am, CS_LAGR_DIAMETER);

    cs_real_t *part_vel = cs_lagr_particle_attr(part, p_am, CS_LAGR_VELOCITY);
    cs_real_t *prev_part_vel
      = cs_lagr_particle_attr_n(part, p_am, 1, CS_LAGR_VELOCITY);

    test_colli = 0;

    cs_lnum_t iel = cs_lagr_particle_get_lnum(part, p_am, CS_LAGR_CELL_ID);

    cs_real_t temp;

    if (extra->scal_t != NULL) {

      if (   cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_TEMPERATURE
          && cs_glob_thermal_model->itpscl == CS_TEMPERATURE_SCALE_CELSIUS)
        temp = extra->scal_t->val[iel] + tkelvi;

      else if (   cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_TEMPERATURE
               && cs_glob_thermal_model->itpscl ==
                          CS_TEMPERATURE_SCALE_KELVIN)
        temp = extra->scal_t->val[iel];

      else if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_ENTHALPY) {

        cs_lnum_t mode = 1;
        CS_PROCF (usthht,USTHHT)(&mode, &(extra->scal_t->val[iel]), &temp);

      }

    }
    else
      temp = cs_glob_fluid_properties->t0;

    cs_lnum_t flag = cs_lagr_particle_get_lnum(part, p_am, CS_LAGR_P_FLAG);
    cs_real_t diam_mean = cs_glob_lagr_clogging_model->diam_mean;

    /* Treatment of internal deposition and user imposed motion */
    if ((flag & CS_LAGR_PART_IMPOSED_MOTION) && face_id > -1) {

      /* Reorient the face so that it is the outwarding normal */
      int reorient_face = 1;
      if (iel == mesh->i_face_cells[face_id][1])
        reorient_face = -1;

      /* Adhesion forces not implemented */

      /* Gravity forces */
      cs_real_3_t gravity = {cs_glob_physical_constants->gravity[0],
        cs_glob_physical_constants->gravity[1],
        cs_glob_physical_constants->gravity[2]};

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
        cs_lagr_particles_unset_flag(p_set, ip, CS_LAGR_PART_DEPOSITION_FLAGS);
        /* TODO: impose particle velocity? Do some stats? */
      }
    }

    /* Monolayer resuspension model */
    if ((    cs_glob_lagr_model->clogging == 0 && face_id > 0)
         || (   cs_glob_lagr_model->clogging == 1 && face_id > 0
             && bound_stat[face_id + n_faces * lag_bi->ihdepm] < diam_mean)) {

      if (flag & CS_LAGR_PART_DEPOSITED)
        /* The particle has just deposited     */
        /* The adhesion force is calculated    */
        cs_lagr_adh(ip, temp, &adhesion_energ);

      else if (flag & CS_LAGR_PART_ROLLING) {

        /* The particle is rolling   */

        /* if the number of great asperities   */
        /* is null it is marked for a possible collision */
        if (cs_lagr_particle_get_lnum(part, p_am, CS_LAGR_N_LARGE_ASPERITIES) == 0)
          test_colli = 1;

        cs_real_t disp_norm
          = cs_lagr_particle_get_real(part, p_am, CS_LAGR_DISPLACEMENT_NORM);

        if (disp_norm > p_diam && disp_norm < 2.0 * p_diam) {

          /* If the particle has a displacement approximately   *
           * equal to a diameter, recalculation of the adhesion force     */

          cs_lagr_particle_set_real(part, p_am, CS_LAGR_DISPLACEMENT_NORM, 0.0);

          cs_lagr_adh(ip, temp, &adhesion_energ);

          if (   test_colli == 1
              && cs_lagr_particle_get_lnum(part, p_am,
                                           CS_LAGR_N_LARGE_ASPERITIES) > 0) {

          cs_real_t kinetic_energy =  0.5 * p_mass
                                   * cs_math_3_dot_product(part_vel, part_vel);

            if (kinetic_energy > adhesion_energ) {

              /* The particle is resuspended
               * with an angle (determined using the large-scale asperity radius) */

              cs_lagr_particles_unset_flag(p_set, ip, CS_LAGR_PART_DEPOSITION_FLAGS);
              cs_lagr_particle_set_real(part, p_am, CS_LAGR_ADHESION_FORCE, 0.0);
              cs_lagr_particle_set_real(part, p_am, CS_LAGR_ADHESION_TORQUE, 0.0);
              cs_lagr_particle_set_lnum(part, p_am, CS_LAGR_N_LARGE_ASPERITIES, 0);
              cs_lagr_particle_set_lnum(part, p_am, CS_LAGR_N_SMALL_ASPERITIES, 0);
              cs_lagr_particle_set_real(part, p_am, CS_LAGR_DISPLACEMENT_NORM, 0.0);

              cs_real_t norm_face = fvq->b_face_surf[face_id];

              cs_real_t norm_velocity = cs_math_3_norm(part_vel);

              cs_real_t norm_reent = sqrt((kinetic_energy-adhesion_energ)*2.0/p_mass);
              cs_real_t angle_reent =
                acos(p_diam * 0.5 /
                     (p_diam * 0.5 + cs_glob_lagr_reentrained_model->rayasg));

              for (cs_lnum_t id = 0; id < 3; id++) {
                part_vel[id] *= norm_reent / norm_velocity * cos(angle_reent);
                part_vel[id] -= norm_reent * sin(angle_reent) / norm_face
                  * fvq->b_face_normal[face_id * 3 +id];
              }

              /* Update of the number and weight of resuspended particles     */
              p_set->n_part_resusp += 1;
              p_set->weight_resusp += p_stat_weight;

              if (events != NULL)
                _add_resuspension_event(events,
                                        p_set,
                                        ip,
                                        face_id,
                                        part_vel);

            }

          }

        }
        else if (disp_norm >= 2.0 * p_diam) {

          cs_lnum_t ndiam = (cs_lnum_t)(disp_norm / p_diam);

          cs_lnum_t ii = 1;
          while (   ii <= ndiam
                 && (cs_lagr_particles_get_flag(p_set, ip,
                                                CS_LAGR_PART_ROLLING))) {

            cs_lagr_adh(ip, temp, &adhesion_energ);

            /* Reconstruct an estimate of the particle velocity   */
            /* at the current sub-time-step assuming linear variation  */
            /* (constant acceleration)   */

            if (   test_colli == 1
                && cs_lagr_particle_get_lnum(part, p_am,
                                             CS_LAGR_N_LARGE_ASPERITIES) > 0) {

              cs_real_t kinetic_energy =  0.5 * p_mass
                * cs_math_3_square_norm(part_vel);

              if (kinetic_energy > adhesion_energ) {

                /* The particle is resuspended    */
                /* and its kinetic energy is totally converted   */
                /* along the wall-normal distance */
                cs_lagr_particles_unset_flag(p_set, ip,
                                             CS_LAGR_PART_DEPOSITION_FLAGS);
                cs_lagr_particle_set_real(part, p_am, CS_LAGR_ADHESION_FORCE, 0.0);
                cs_lagr_particle_set_real(part, p_am, CS_LAGR_ADHESION_TORQUE, 0.0);
                cs_lagr_particle_set_lnum(part, p_am, CS_LAGR_N_LARGE_ASPERITIES, 0);
                cs_lagr_particle_set_lnum(part, p_am, CS_LAGR_N_SMALL_ASPERITIES, 0);
                cs_lagr_particle_set_real(part, p_am, CS_LAGR_DISPLACEMENT_NORM, 0.0);

                cs_real_t norm_face = cs_glob_mesh_quantities->b_face_surf[face_id];

                cs_real_t norm_velocity = cs_math_3_norm(part_vel);

                cs_real_t norm_reent = sqrt((kinetic_energy-adhesion_energ)*2.0/p_mass);
                cs_real_t angle_reent =
                  acos(p_diam * 0.5 /
                       (p_diam * 0.5 + cs_glob_lagr_reentrained_model->rayasg));

                for (cs_lnum_t id = 0; id < 3; id++) {
                  part_vel[id] *= norm_reent / norm_velocity * cos(angle_reent);
                  part_vel[id] -= norm_reent * sin(angle_reent) / norm_face
                    * cs_glob_mesh_quantities->b_face_normal[face_id * 3 +id];
                }

                /* Update of the number and weight of resuspended particles     */
                p_set->n_part_resusp += 1;
                p_set->weight_resusp += p_stat_weight;

                if (events != NULL)
                  _add_resuspension_event(events,
                                          p_set,
                                          ip,
                                          face_id,
                                          part_vel);

              }

              if (cs_lagr_particle_get_lnum(part, p_am, CS_LAGR_N_LARGE_ASPERITIES) == 0)
                test_colli = 1;

            }

            ii++;

          }

        }

      }
    } /* Enf of monolayer resuspension */
       else if (    cs_glob_lagr_model->clogging == 1
                 && bound_stat[face_id + n_faces * lag_bi->ihdepm] >= diam_mean) {

      /* Treatment of multilayer resuspension */

      if (flag == CS_LAGR_PART_DEPOSITED)
        /* The particle has just deposited     */
        /* The adhesion force is calculated    */
        cs_lagr_adh(ip, temp, &adhesion_energ);

      else if (flag == CS_LAGR_PART_ROLLING) {

        /* The cluster is rolling   */

        cs_real_t norm_face = cs_glob_mesh_quantities->b_face_surf[face_id];

        /*  Mean distance between protruding clusters */
        cs_real_t cluster_spacing;
        if (bound_stat[face_id + n_faces * lag_bi->inclg] <= 0)
          cluster_spacing = sqrt(norm_face);

        else
          cluster_spacing = sqrt( 2.0 * norm_face /
                                  bound_stat[face_id + n_faces * lag_bi->inclg] );

        cs_real_t disp_norm = cs_lagr_particle_get_real(part, p_am,
                                                        CS_LAGR_DISPLACEMENT_NORM);

        cs_lnum_t ndiam = (cs_lnum_t)(disp_norm / cluster_spacing);

        cs_lnum_t ii = 1;
        while (   ii <= ndiam
               && (cs_lagr_particles_get_flag(p_set, ip,
                                              CS_LAGR_PART_ROLLING))) {

          /* Reconstruct an estimate of the kinetic energy   */
          /* at the current sub-time-step assuming a linear variation  */

          cs_real_t kinetic_energy_prev = 0.5 * p_mass
            * cs_math_3_square_norm(prev_part_vel);
          cs_real_t kinetic_energy_cur = 0.5 * p_mass
            * cs_math_3_square_norm(part_vel);

          cs_real_t kinetic_energy = kinetic_energy_prev
            + ii / ndiam * (kinetic_energy_cur - kinetic_energy_prev);

          /* Adhesion energy upon rocking */

          cs_real_t adhes_energ, adhes_force, adhes_torque;
          cs_lagr_adh_pp(p_diam, temp, &adhes_energ, &adhes_force);
          /* Average number of contact in a cluster */
          cs_real_t ncont_pp = pow(p_diam/diam_mean, 2);

          cs_lnum_t ncont = 1;

          if (ncont_pp > 600.0) {
            cs_real_t rtmp;
            cs_random_normal(1, &rtmp);
            ncont = (int)ncont_pp + sqrt(ncont_pp) * rtmp;
          }
          else {
            cs_random_poisson(1, ncont_pp, &ncont);
          }
          ncont = CS_MAX(1, ncont);
          cs_lagr_particle_set_lnum(part, p_am, CS_LAGR_N_SMALL_ASPERITIES, ncont);

          adhes_energ *= ncont;
          adhes_force *= ncont ;
          cs_lagr_particle_set_real(part, p_am, CS_LAGR_ADHESION_FORCE, adhes_force);

          adhes_torque = adhes_force * p_diam * 0.5;
          cs_lagr_particle_set_real(part, p_am, CS_LAGR_ADHESION_TORQUE, adhes_torque);

          if (kinetic_energy > adhesion_energ) {

            /* The particle is resuspended    */
            /* and its kinetic energy is totally converted   */
            /* along the wall-normal distance */
            cs_lagr_particles_unset_flag(p_set, ip,
                                         CS_LAGR_PART_DEPOSITION_FLAGS);
            cs_lagr_particle_set_real(part, p_am, CS_LAGR_ADHESION_FORCE, 0.0);
            cs_lagr_particle_set_real(part, p_am, CS_LAGR_ADHESION_TORQUE, 0.0);
            cs_lagr_particle_set_lnum(part, p_am, CS_LAGR_N_LARGE_ASPERITIES, 0);
            cs_lagr_particle_set_lnum(part, p_am, CS_LAGR_N_SMALL_ASPERITIES, 0);
            cs_lagr_particle_set_real(part, p_am, CS_LAGR_DISPLACEMENT_NORM, 0.0);

            cs_real_t norm_velocity = cs_math_3_norm(part_vel);

            cs_real_t norm_reent = sqrt((kinetic_energy-adhesion_energ)*2.0/p_mass);
            cs_real_t angle_reent =
              acos(p_diam * 0.5 /
                   (p_diam * 0.5 + cs_glob_lagr_reentrained_model->rayasg));

            for (cs_lnum_t id = 0; id < 3; id++) {
              part_vel[id] *= norm_reent / norm_velocity * cos(angle_reent);
              part_vel[id] -= norm_reent * sin(angle_reent) / norm_face
                * cs_glob_mesh_quantities->b_face_normal[face_id * 3 +id];
            }

            /* Update of the number and weight of resuspended particles     */
            p_set->n_part_resusp += 1;
            p_set->weight_resusp += p_stat_weight;

            if (events != NULL)
              _add_resuspension_event(events,
                                      p_set,
                                      ip,
                                      face_id,
                                      part_vel);

          }

          ii++;

        }

      }


    } /* End of multilayer resuspension */

  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
