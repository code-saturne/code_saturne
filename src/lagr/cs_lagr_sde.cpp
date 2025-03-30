/*============================================================================
 * Methods for particle deposition
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
 * Functions dealing with particle deposition
 *============================================================================*/

#include "base/cs_defs.h"
#include "base/cs_math.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <math.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_mem.h"
#include "bft/bft_printf.h"

#include "base/cs_math.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_physical_constants.h"
#include "pprt/cs_physical_model.h"
#include "base/cs_prototypes.h"
#include "base/cs_random.h"
#include "base/cs_rotation.h"
#include "base/cs_thermal_model.h"
#include "turb/cs_turbulence_model.h"

#include "lagr/cs_lagr.h"
#include "lagr/cs_lagr_adh.h"
#include "lagr/cs_lagr_deposition_model.h"
#include "lagr/cs_lagr_event.h"
#include "lagr/cs_lagr_roughness.h"
#include "lagr/cs_lagr_tracking.h"
#include "lagr/cs_lagr_prototypes.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "lagr/cs_lagr_sde.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Boltzmann constant */
static const double _k_boltz = 1.38e-23;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Private function prototypes (definitions follow)
 *============================================================================*/

/*============================================================================
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
                        const cs_real_t          particle_velocity[3])
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

  auto *e_flag = cs_lagr_events_attr_get_ptr<cs_lnum_t>(events,
                                                        event_id,
                                                        CS_LAGR_E_FLAG);

  auto *e_vel_post = cs_lagr_events_attr_get_ptr<cs_real_t>(events,
                                                            event_id,
                                                            CS_LAGR_E_VELOCITY);

  *e_flag = *e_flag | CS_EVENT_RESUSPENSION;

  for (int k = 0; k < 3; k++)
    e_vel_post[k] = particle_velocity[k];

  events->n_events += 1;
}

/*----------------------------------------------------------------------------*/
/*! \brief Integration of SDEs by 1st order time scheme for one particle
 *
 * \param[in]  p_id      particle index in set
 * \param[in]  dt_part   remaining time step associated to the particle
 * \param[in]  nor       current step id (for 2nd order scheme)
 * \param[in]  taup      dynamic characteristic time
 * \param[in]  tlag      lagrangian fluid characteristic time
 * \param[in]  piil      term in integration of up sdes
 * \param[in]  bx        turbulence characteristics
 * \param[in]  vagaus    gaussian random variables
 * \param[in]  brgaus    gaussian random variables
 * \param[in]  force_p   forces per mass unit on particles (m/s^2)
 * \param[in]  beta      proportional to the gradient of T_lag
 */
/*----------------------------------------------------------------------------*/
void
cs_sde_vels_pos_1_st_order_time_integ(cs_lnum_t                       p_id,
                                      cs_real_t                       dt_part,
                                      int                             nor,
                                      const cs_real_t                *taup,
                                      const cs_real_3_t              *tlag,
                                      const cs_real_3_t              *piil,
                                      const cs_real_33_t             *bx,
                                      const cs_real_3_t              *vagaus,
                                      const cs_real_6_t               brgaus,
                                      const cs_real_3_t               force_p,
                                      const cs_real_3_t               beta)
{
  /* Particles management */
  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;
  unsigned char *particle = p_set->p_buffer + p_am->extents * p_id;

  /* use previous step for t_order == 1 or prediction step
   * and current one for correction step */
  cs_lnum_t cell_id = cs_lagr_particle_get_lnum_n(particle, p_set->p_am, 2-nor,
                                                  CS_LAGR_CELL_ID);


  if (cell_id <0)
    return;
  cs_lagr_extra_module_t *extra_i = cs_get_lagr_extra_module();
  cs_lagr_extra_module_t *extra = extra_i;
  int n_phases = extra->n_phases;

  /* Initialisations*/

  cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

  cs_real_t aux1, aux1m, aux2, aux3, aux3m, aux4,
            aux44, aux5, aux6, aux7, aux8,
            aux9, aux10, aux11;
  cs_real_t ter1f[n_phases], ter2f[n_phases], ter3f[n_phases];
  cs_real_t ter1p, ter2p[n_phases], ter3p[n_phases], ter4p, ter5p, ter6p[n_phases];
  cs_real_t ter1x, ter2x[n_phases], ter3x[n_phases], ter4x, ter5x, ter6x[n_phases];
  // Terms of the reduction matrix
  cs_real_t p11[n_phases], p21[n_phases], p22, p31[n_phases], p32, p33;
  // Terms of the covariance matrix for the Gauss vector
  cs_real_t gam2[n_phases];
  cs_real_t gam_gagam[n_phases], gagam2;
  cs_real_t gam_ome[n_phases], gagam_ome, ome2;
  cs_real_t tbrix1, tbrix2, tbriu;
  cs_real_t ter7x, ter7p, ter7f;
  ter7x = 0.;
  ter7p = 0.;
  ter7f = 0.;

  const int _prev_id = (extra->vel->n_time_vals > 1) ? 1 : 0;

  const cs_temperature_scale_t t_scl = cs_glob_thermal_model->temperature_scale;

  const cs_real_3_t **cvar_vel = nullptr;
  CS_MALLOC(cvar_vel, n_phases, const cs_real_3_t *);
  for (int phase_id = 0; phase_id < n_phases; phase_id++)
    cvar_vel[phase_id] =
      (const cs_real_3_t *)(extra_i[phase_id].vel->vals[_prev_id]);

  /* Fluid particles need no change of reference frame, others do */

  bool local_reference_frame = true;
  if (   cs_glob_lagr_model->shape == CS_LAGR_SHAPE_SPHERE_MODEL
      && cs_glob_lagr_model->modcpl == 0)
    local_reference_frame = false;

  /* Obtain the mean particle velocity for each cell, if present */

  cs_field_t *stat_vel = nullptr;
  cs_field_t *stat_vel_s = nullptr;

  if (   cs_glob_lagr_model->shape == CS_LAGR_SHAPE_SPHERE_MODEL
      && cs_glob_lagr_model->modcpl) {
    int stat_type = cs_lagr_stat_type_from_attr_id(CS_LAGR_VELOCITY);
    stat_vel = cs_lagr_stat_get_moment(stat_type,
                                       CS_LAGR_STAT_GROUP_PARTICLE,
                                       CS_LAGR_MOMENT_MEAN,
                                       0, -1);
    stat_type = cs_lagr_stat_type_from_attr_id(CS_LAGR_VELOCITY_SEEN);
    stat_vel_s = cs_lagr_stat_get_moment(stat_type,
                                       CS_LAGR_STAT_GROUP_PARTICLE,
                                       CS_LAGR_MOMENT_MEAN,
                                       0, -1);
  }

  cs_real_t *lambda;
  cs_real_3_t *loc_fluid_vel = nullptr;
  cs_real_3_t *part_vel_seen_r = nullptr;
  cs_real_3_t *old_part_vel_seen_r = nullptr;
  cs_real_3_t *fluid_vel_r = nullptr;
  cs_real_3_t *piil_r = nullptr;
  cs_real_3_t *tlag_r = nullptr;
  cs_real_3_t *taup_r = nullptr;
  CS_MALLOC(lambda, n_phases, cs_real_t);
  CS_MALLOC(loc_fluid_vel, n_phases, cs_real_3_t);
  CS_MALLOC(part_vel_seen_r, n_phases, cs_real_3_t);
  CS_MALLOC(old_part_vel_seen_r, n_phases, cs_real_3_t);
  CS_MALLOC(fluid_vel_r, n_phases, cs_real_3_t);
  CS_MALLOC(piil_r, n_phases, cs_real_3_t);
  CS_MALLOC(tlag_r, n_phases, cs_real_3_t);
  CS_MALLOC(taup_r, n_phases, cs_real_3_t);

  /* Integrate SDE's over particles
   * Note: new particles will be integrated at the next time step, otherwise
   * positions might be overwritten */

  const cs_real_t *cell_cen = cs_glob_mesh_quantities->cell_cen[cell_id];

  /* Get particle coordinates, velocity and velocity seen*/

  /* --> lambda represents the fraction of time spent in phase_id (normalized)
  In this variable we can take into account several different physical Properties
  such as: "the particle cannot exit phase x" for example */
  for (int phase_id = 0; phase_id < n_phases; phase_id++) {
    if (cs_glob_lagr_model->cs_used == 0)
      lambda[phase_id] = extra_i[phase_id].alpha->val[cell_id];
    else
      lambda[phase_id] = 1.;
  }

  /* Get particle coordinates, velocity and velocity seen*/

  auto *old_part_vel      =
    cs_lagr_particle_attr_n_get_ptr<cs_real_t>(particle, p_am, 1,
                                               CS_LAGR_VELOCITY);
  auto *old_part_vel_seen =
    cs_lagr_particle_attr_n_get_ptr<cs_real_t>(particle, p_am, 1,
                                               CS_LAGR_VELOCITY_SEEN);
  auto *old_part_coords   =
    cs_lagr_particle_attr_n_get_ptr<cs_real_t>(particle, p_am, 1,
                                               CS_LAGR_COORDS);
  auto *part_vel          =
    cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                             CS_LAGR_VELOCITY);
  auto *part_vel_seen     =
    cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                             CS_LAGR_VELOCITY_SEEN);
  auto *part_coords       =
    cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                             CS_LAGR_COORDS);

  int imposed_motion = cs_lagr_particles_get_flag(p_set, p_id,
                                                  CS_LAGR_PART_IMPOSED_MOTION);
  if (imposed_motion) {
    cs_real_t disp[3] = {0., 0., 0.};

    cs_user_lagr_imposed_motion(p_set,
                                p_id,
                                old_part_coords,
                                dt_part,
                                disp);

    for (cs_lnum_t id = 0; id < 3; id++) {

      part_coords[id] = old_part_coords[id] + disp[id];

      for (int phase_id = 0; phase_id < n_phases; phase_id++)
        part_vel_seen[3 * phase_id + id] =  0.0;

      part_vel[id] = disp[id] / dt_part;

    }
    return;
  } /* End IMPOSED_MOTION */

  /* resolve SDEs*/
  for (int phase_id = 0; phase_id < n_phases; phase_id++) {
    if (cs_glob_lagr_time_scheme->interpol_field == 1) {
      for (int i = 0; i < 3; i++) {
        loc_fluid_vel[phase_id][i] = cvar_vel[phase_id][cell_id][i];
        for (int j = 0; j < 3; j++)
          loc_fluid_vel[phase_id][i] += extra_i[phase_id].grad_vel[cell_id][i][j]
                            * (part_coords[j] - cell_cen[j]);
      }
    }
    else {
      for (int i = 0; i < 3; i++) {
        loc_fluid_vel[phase_id][i] = cvar_vel[phase_id][cell_id][i];
      }
    }
  }

  /* Initialize (without change of frame)*/

  cs_real_t part_vel_r[3] = {part_vel[0], part_vel[1], part_vel[2]};
  cs_real_t old_part_vel_r[3] = {old_part_vel[0],
                                 old_part_vel[1],
                                 old_part_vel[2]};

  for (int phase_id = 0; phase_id < n_phases; phase_id++) {
    for (int idim_ = 0; idim_ < 3; idim_++) {
      part_vel_seen_r[phase_id][idim_] = part_vel_seen[phase_id * 3 + idim_];
      old_part_vel_seen_r[phase_id][idim_] = old_part_vel_seen[phase_id * 3 + idim_];
      fluid_vel_r[phase_id][idim_] = loc_fluid_vel[phase_id][idim_];
    }
  }

  cs_real_t force_p_r[3] = {force_p[0], force_p[1], force_p[2]};
  cs_real_t taup_rm[3] = {0., 0., 0.};

  cs_real_t displ_r[3];
  cs_real_t trans_m[3][3];
  for (int phase_id = 0; phase_id < n_phases; phase_id ++) {
    for (int idim_ = 0; idim_ < 3; idim_ ++) {
      piil_r[phase_id][idim_] = piil[phase_id][idim_];
      tlag_r[phase_id][idim_] = tlag[phase_id][idim_];
      taup_r[phase_id][idim_] = taup[phase_id];
      taup_rm[idim_] += lambda[phase_id] / taup_r[phase_id][idim_];
    }
  }

  for (int idim_ = 0; idim_ < 3; idim_ ++)
    taup_rm[idim_] = 1./taup_rm[idim_];

  if (local_reference_frame) {

    /* Global reference frame --> local reference frame depending on model
       =================================================================== */

    switch (cs_glob_lagr_model->shape) {

    case CS_LAGR_SHAPE_SPHERE_MODEL:
      {
        /* Compute the main direction in the global reference
         * frame */
        cs_real_t dir[3];

        /* Obtain the mean particle velocity for each cell */
        const cs_real_t *mean_part_vel_p = stat_vel->val + (cell_id * 3);
        const cs_real_t *mean_part_vel_s = stat_vel_s->val + (cell_id * 3);

        /* Get Relative mean velocity */
        if (cs_glob_lagr_model->cs_used == 1) {
          /* Relative mean velocity
           * <u_p> - <u_s> = <u_r>
           * See Minier et al 2024.
           * */
          for (cs_lnum_t i = 0; i < 3; i++)
            dir[i] = mean_part_vel_p[i] - mean_part_vel_s[i];
        }
        else { /* neptune is used */
          /* Change the direction since mean part vel contains
           * mean velocity fluctuation */
          cs_real_t fluid_vel_mel[3] = {0., 0., 0.};
          for (int phase_id = 0; phase_id < n_phases; phase_id ++) {
            const cs_real_t *fluid_vel = cvar_vel[phase_id][cell_id];
            for (cs_lnum_t i = 0; i < 3; i++) {
              fluid_vel_mel[i] += lambda[phase_id] * fluid_vel[i];
            }
          }
          for (cs_lnum_t i = 0; i < 3; i++)
            dir[i] = mean_part_vel_p[i] - fluid_vel_mel[i];
        }

        cs_math_3_normalize(dir, dir);

        // Rotate the frame of reference with respect to the
        // mean relative velocity direction.

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
        const cs_real_t euler[4] =
        {  cs_math_3_norm(dir_p_dir_r)
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

      }
      break;

    case CS_LAGR_SHAPE_SPHEROID_JEFFERY_MODEL:
      {
        // Use Euler angles for spheroids (Jeffery)
        const cs_real_t *euler =
          cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                   CS_LAGR_EULER);

        trans_m[0][0] = 2.*(euler[0]*euler[0]+euler[1]*euler[1]-0.5);
        trans_m[0][1] = 2.*(euler[1]*euler[2]+euler[0]*euler[3]);
        trans_m[0][2] = 2.*(euler[1]*euler[3]-euler[0]*euler[2]);
        trans_m[1][0] = 2.*(euler[1]*euler[2]-euler[0]*euler[3]);
        trans_m[1][1] = 2.*(euler[0]*euler[0]+euler[2]*euler[2]-0.5);
        trans_m[1][2] = 2.*(euler[2]*euler[3]+euler[0]*euler[1]);
        trans_m[2][0] = 2.*(euler[1]*euler[3]+euler[0]*euler[2]);
        trans_m[2][1] = 2.*(euler[2]*euler[3]-euler[0]*euler[1]);
        trans_m[2][2] = 2.*(euler[0]*euler[0]+euler[3]*euler[3]-0.5);
      }
      break;

    case CS_LAGR_SHAPE_SPHEROID_STOC_MODEL:
      {
        // Use rotation matrix for stochastic model
          cs_real_t *orient_loc  =
            cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                     CS_LAGR_ORIENTATION);
        cs_real_t singularity_axis[3] = {1.0, 0.0, 0.0};
        // Get vector for rotation
        cs_real_t n_rot[3];
        cs_math_3_cross_product(orient_loc, singularity_axis, n_rot);
        cs_math_3_normalize(n_rot, n_rot);
        // Compute rotation angle
        cs_real_t cosa = cs_math_3_dot_product(orient_loc, singularity_axis);
        cs_real_t sina = sin(acos(cosa));

        // Compute the rotation matrix
        trans_m[0][0] = cosa + cs_math_pow2(n_rot[0])*(1.0 - cosa);
        trans_m[0][1] = n_rot[0]*n_rot[1]*(1.0 - cosa) - n_rot[2]*sina;
        trans_m[0][2] = n_rot[0]*n_rot[2]*(1.0 - cosa) + n_rot[1]*sina;
        trans_m[1][0] = n_rot[0]*n_rot[1]*(1.0 - cosa) + n_rot[2]*sina;
        trans_m[1][1] = cosa + cs_math_pow2(n_rot[1])*(1.0 - cosa);
        trans_m[1][2] = n_rot[1]*n_rot[2]*(1.0 - cosa) - n_rot[0]*sina;
        trans_m[2][0] = n_rot[0]*n_rot[2]*(1.0 - cosa) - n_rot[1]*sina;
        trans_m[2][1] = n_rot[1]*n_rot[2]*(1.0 - cosa) + n_rot[0]*sina;
        trans_m[2][2] = cosa + cs_math_pow2(n_rot[2])*(1.0 - cosa);
      }
        break;

      default:
        assert(0);
        bft_error
          (__FILE__, __LINE__, 0,
           _("Change of reference frame not implemented\n"));

      }

    /* 1.1 - particle velocity */

    cs_math_33_3_product(trans_m, old_part_vel, old_part_vel_r);

    /* 1.2 - flow-seen velocity  */

    for (int phase_id = 0; phase_id < n_phases; phase_id++)
      cs_math_33_3_product(trans_m,
          old_part_vel_seen + 3*phase_id,
          old_part_vel_seen_r[phase_id]);

    /* 1.3 Particle force: - pressure gradient/romp + external force + g */

    cs_math_33_3_product(trans_m, force_p, force_p_r);

    for (int phase_id = 0; phase_id < n_phases; phase_id ++) {
      /* 1.4 - flow velocity  (current and previous) */

      cs_math_33_3_product(trans_m,
                           loc_fluid_vel[phase_id],
                           fluid_vel_r[phase_id]);

      /* 1.5 - "piil" term    */

      cs_math_33_3_product(trans_m, piil[phase_id], piil_r[phase_id]);

    }

    /* 1.6 - taup  */

    if (   cs_glob_lagr_model->shape == CS_LAGR_SHAPE_SPHEROID_STOC_MODEL
          || cs_glob_lagr_model->shape == CS_LAGR_SHAPE_SPHEROID_JEFFERY_MODEL) {

        cs_real_t *radii =
          cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                   CS_LAGR_RADII);

      cs_real_t *s_p =
        cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                   CS_LAGR_SHAPE_PARAM);

      for (int phase_id = 0; phase_id < n_phases; phase_id ++) {
        taup_r[phase_id][0] = 3.0 / 8.0 * taup[phase_id]
          * (radii[0] * radii[0] * s_p[0] + s_p[3])
          / pow(radii[0] * radii[1] * radii[2], 2.0 / 3.0);
        taup_r[phase_id][1] = 3.0 / 8.0 * taup[phase_id]
          * (radii[1] * radii[1] * s_p[1] + s_p[3])
          / pow(radii[0] * radii[1] * radii[2], 2.0 / 3.0);
        taup_r[phase_id][2] = 3.0 / 8.0 * taup[phase_id]
          * (radii[2] * radii[2] * s_p[2] + s_p[3])
          / pow(radii[0] * radii[1] * radii[2], 2.0 / 3.0);
      }
    }
  }

  /* Integration of the SDE on the particles
   * ======================================= */

  for (cs_lnum_t id = 0; id < 3; id++) {
    gagam2 = 0.;
    gagam_ome = 0.;
    ome2 = 0.;

    /* Preliminary computation:
       ------------------------ */

    cs_real_t tci;
    /* velocity induced by a force by unit mass */
    cs_real_t v_lim = force_p_r[id] * taup_rm[id];

    /* Compute deterministic coefficients/terms
       ---------------------------------------- */

    aux1m = exp(-dt_part / taup_rm[id]);
    cs_real_t bb, cc, dd, ff;

    cs_real_t aa = taup_rm[id] * (1.0 - aux1m);
    cs_real_t ee = 1.0 - aux1m;

    /* --> first and second increment for particle velocity */
    ter1p = old_part_vel_r[id] * aux1m;
    ter4p = v_lim * ee;

    /* --> first and second increment for particle position */
    ter1x = aa * old_part_vel_r[id];
    ter4x = (dt_part - aa) * v_lim;

    /* Integral for the particles velocity */
    aux10 = 0.5 * taup_rm[id] * (1.0 - aux1m * aux1m);

    for (int phase_id = 0; phase_id < n_phases; phase_id ++) {
      tci = piil_r[phase_id][id] * tlag_r[phase_id][id];
      if (cs_glob_lagr_model->cs_used == 1)
        tci += fluid_vel_r[phase_id][id];

      cs_real_t multiphase_coef = lambda[phase_id] * taup_rm[id] / taup_r[phase_id][id];
      cs_real_t multiphase_coef2 = cs_math_pow2(multiphase_coef);

      /* Compute deterministic coefficients/terms
      ---------------------------------------- */

      aux1 = exp(-dt_part / taup_r[phase_id][id]);
      aux2 = exp(-dt_part / tlag_r[phase_id][id]);
      aux3 = tlag_r[phase_id][id]
        / (tlag_r[phase_id][id] - taup_r[phase_id][id]);
      aux3m = tlag_r[phase_id][id]
        / (tlag_r[phase_id][id] - taup_rm[id]);
      aux4 = tlag_r[phase_id][id]
        / (tlag_r[phase_id][id] + taup_r[phase_id][id]);
      aux44 = taup_rm[id] / (tlag_r[phase_id][id] + taup_rm[id]);

      aux5 = tlag_r[phase_id][id] * (1.0 - aux2);
      aux6 = cs_math_pow2(bx[phase_id][id][nor-1])
        * tlag_r[phase_id][id];
      aux7 = tlag_r[phase_id][id] - taup_r[phase_id][id];
      aux8 = cs_math_pow2(bx[phase_id][id][nor-1])
        * cs_math_pow2(aux3m);
      /* --> trajectory terms */
      bb = (aux5 - aa) * aux3m;
      cc = dt_part - aa - bb;
      ff = lambda[phase_id] / taup_r[phase_id][id] * taup_rm[id];

      ter2f[phase_id] = tci * (1.0 - aux2);

      /* Integral on flow velocity seen */
      gam2[phase_id]  = 0.5 * (1.0 - aux2 * aux2);
      p11[phase_id]   = sqrt(gam2[phase_id] * aux6);
      ter3f[phase_id] = p11[phase_id] * vagaus[phase_id][id];

      /* Terms for particle velocity */
      dd = aux3m * (aux2 - aux1m);

      ter3p[phase_id] = tci * ff * (ee - dd);

      /* Integral for the particles velocity */
      aux9 = 0.5 * tlag_r[phase_id][id] * (1.0 - aux2 * aux2);
      aux11 =   taup_rm[id] * tlag_r[phase_id][id]
              * (1.0 - aux1m * aux2)
              / (taup_rm[id] + tlag_r[phase_id][id]);

      ter6p[phase_id] = lambda[phase_id] / taup_r[phase_id][id]
        * aa * fluid_vel_r[phase_id][id];

      ter6x[phase_id] = ff * (dt_part - aa) * fluid_vel_r[phase_id][id];

      ter3x[phase_id] = ff * cc * tci;

      if (cs_glob_lagr_model->cs_used == 1) {
        /* Flow-seen velocity terms */
        ter1f[phase_id] = old_part_vel_seen_r[phase_id][id] * aux2;

        /* Terms for particle velocity */
        ter2p[phase_id] = old_part_vel_seen_r[phase_id][id] * dd;
        ter6p[phase_id] = 0.0;
        ter6x[phase_id] = 0.0;

        /* Terms for particle position */
        ter2x[phase_id] = bb * old_part_vel_seen_r[phase_id][id] ;

        /* Compute the other terms for the covariance matrix of the Gauss vector */
        gam_gagam[phase_id] = (aux9 - aux11) * (aux8 / aux3);

        gagam2 = (aux9 - 2.0 * aux11 + aux10) * aux8;

        gam_ome[phase_id] = ( (tlag_r[phase_id][id] - taup_r[phase_id][id])
            * (aux5 - aa)
              - tlag_r[phase_id][id] * aux9
              - taup_r[phase_id][id] * aux10
              + (tlag_r[phase_id][id] + taup_r[phase_id][id]) * aux11)
              * aux8;

        gagam_ome = aux3 * (  (tlag_r[phase_id][id] - taup_r[phase_id][id])
            * (1.0 - aux2)
                     - 0.5 * tlag_r[phase_id][id] * (1.0 - aux2 * aux2)
                     + cs_math_pow2(taup_r[phase_id][id])
                     / (tlag_r[phase_id][id] + taup_r[phase_id][id])
                     * (1.0 - aux1 * aux2)) * aux6;

        ome2 = aux7 * (aux7 * dt_part - 2.0
              * (tlag_r[phase_id][id] * aux5 - taup_r[phase_id][id] * aa))
             + 0.5 * tlag_r[phase_id][id] * tlag_r[phase_id][id] * aux5
             * (1.0 + aux2)
             + 0.5 * taup_r[phase_id][id] * taup_r[phase_id][id] * aa * (1.0 + aux1)
             - 2.0 * aux4 * tlag_r[phase_id][id]
             * taup_r[phase_id][id] * taup_r[phase_id][id]
                   * (1.0 - aux1 * aux2);
        ome2 = ome2 * aux8;
      } else {
        /* Flow-seen velocity terms */
        ter1f[phase_id] = (old_part_vel_seen_r[phase_id][id]
                          - fluid_vel_r[phase_id][id]) * aux2;

        /* Terms for particle velocity */
        ter2p[phase_id] = (old_part_vel_seen_r[phase_id][id]
            - fluid_vel_r[phase_id][id])
          * dd * ff;

        /* Terms for particle position */
        ter2x[phase_id] = ff * bb
          * (old_part_vel_seen_r[phase_id][id] - fluid_vel_r[phase_id][id]);

        /* Compute the other terms for the covariance matrix of the Gauss vector */
        gam_gagam[phase_id] = aux6 * multiphase_coef * aux3m
                      * (0.5 * (1. - aux2 * aux2)
                      - aux44 * (1 - aux1m * aux2));

        gagam2 += multiphase_coef2 * cs_math_pow2(aux3m) * aux6
                * (0.5 * (1. - aux2 * aux2)
                - 2 * aux44 * (1. - aux2 * aux1m)
                + 0.5 * taup_rm[id] / tlag_r[phase_id][id] * (1. - aux1m * aux1m));

        gam_ome[phase_id] = multiphase_coef * aux6 * aux3m *
              (- 0.5 * tlag_r[phase_id][id] * (1. - aux2 * aux2)
               + tlag_r[phase_id][id] * (1. - aux2)
               + taup_rm[id] * aux44 * (1. - aux2 * aux1m)
               - taup_rm[id] * (1. - aux2));

        gagam_ome += multiphase_coef2 * cs_math_pow2(aux3m)
          * aux6 / tlag_r[phase_id][id]
               *((tlag_r[phase_id][id] - taup_rm[id]) * (aux5 - aa)
               - 0.5*cs_math_pow2(tlag_r[phase_id][id]) * (1. - aux2 * aux2)
               - 0.5*cs_math_pow2(taup_rm[id]) * (1. - aux1m * aux1m)
               + tlag_r[phase_id][id] * taup_rm[id]*(1. - aux1m * aux2 ));

        ome2 += multiphase_coef2 * cs_math_pow2(aux3m)
          * aux6 / tlag_r[phase_id][id]
               *((tlag_r[phase_id][id] - taup_rm[id])
                   * ((tlag_r[phase_id][id] - taup_rm[id]) * dt_part
                     - 2.*(tlag_r[phase_id][id]*aux5 - taup_rm[id]*aa))
               +0.5*cs_math_pow2(tlag_r[phase_id][id]) * aux5 * (1 + aux2)
               +0.5*cs_math_pow2(taup_rm[id]) * aa * (1 + aux1m)
               -2.*tlag_r[phase_id][id] * tlag_r[phase_id][id] * taup_rm[id] * taup_rm[id]
                        /(tlag_r[phase_id][id] + taup_rm[id]) * (1 - aux1m * aux2));

      }
    }

    /* Additional terms when gradient of Tl is not negligible
     * */
    /* particle positions term */
    cs_real_t aux5b = taup_r[0][id] * (1.0 - aux1);
    ter7x = beta[id] * (
        aux4 * cs_math_pow2(taup_r[0][id]) * dt_part
      + cs_math_pow2(tlag_r[0][id]) * dt_part * aux2
      - cs_math_pow3(tlag_r[0][id])*(1. - aux2)
      + 0.5 * tlag_r[0][id] * aux7 * (dt_part - taup_r[0][id] * (1. - aux2))
      + aux3 * taup_r[0][id] * (2. * tlag_r[0][id] - taup_r[0][id]) * (aux5 - aux5b)
      - 0.5 * cs_math_pow3(tlag_r[0][id]) / (tlag_r[0][id] - 2. * taup_r[0][id])
        * (aux9 - aux5b)
      - cs_math_pow2(cs_math_pow2(taup_r[0][id]) * aux4) * (1. - aux1*aux2)
      );

    /* particle velocity term */
    ter7p = beta[id] * (
      cs_math_pow2(taup_r[0][id]) * aux4 * (1. - aux1 * aux2)
      - tlag_r[0][id] * dt_part * aux2
      + 0.5 * tlag_r[0][id] * (tlag_r[0][id] - 2. * taup_r[0][id]) * (1. - aux1)
      + taup_r[0][id] * aux3 * (2. * tlag_r[0][id] - taup_r[0][id]) * (aux2 - aux1)
      - 0.5 * cs_math_pow3(tlag_r[0][id]) / (tlag_r[0][id] - 2. * taup_r[0][id])
        * (aux2 * aux2 - aux1)
        );
    /* velocity seen by the particle */
    ter7f = beta[id] * (
        tlag_r[0][id] * aux7 * (1. - (1. + dt_part / tlag_r[0][id]) * aux2)
        - 0.5 *cs_math_pow2(aux5)
        + cs_math_pow2(taup_r[0][id]) / (tlag_r[0][id] + taup_r[0][id])
          * (aux5 - aux2 * aux5b)
        );

    /* Now we can compute the pij and then ter5p, ter5x
    pij are the components of the matrix from the Cholesky decomposition of Cij
    (covariance matrix of the Gauss vector) */

    p22 = gagam2;
    p32 = gagam_ome;
    p33 = ome2;
    for (int phase_id = 0; phase_id < n_phases; phase_id++){
      if (cs::abs(p11[phase_id]) > cs_math_epzero) {
        p21[phase_id] = gam_gagam[phase_id] / p11[phase_id];
        p31[phase_id] = gam_ome[phase_id] / p11[phase_id];
      } else {
        p21[phase_id] = 0.;
        p31[phase_id] = 0.;
      }
    }
    for (int phase_id = 0; phase_id < n_phases; phase_id++){
      p22 -= cs_math_pow2(p21[phase_id]);
      p32 -= p31[phase_id] * p21[phase_id];
      p33 -= cs_math_pow2(p31[phase_id]);
    }

    p22 = sqrt(cs::max(0.0, p22));

    if (p22 > cs_math_epzero)
      p32 /= p22;
    else
      p32 = 0.;

    p33 -= cs_math_pow2(p32);
    p33 = sqrt(cs::max(0.0, p33));
    ter5p = 0.;
    ter5x = 0.;
    for (int phase_id = 0; phase_id < n_phases; phase_id++){
      ter5p += p21[phase_id] * vagaus[phase_id][id];
      ter5x += p31[phase_id] * vagaus[phase_id][id];
    }
    ter5p += p22 * vagaus[n_phases][id];
    ter5x += p32 * vagaus[n_phases][id] + p33 * vagaus[n_phases + 1][id];

      /* (2.3) Compute terms in the Brownian movement case */
      /* TODO: Warning based on the first carrier phase with N fluids */
      if (cs_glob_lagr_brownian->lamvbr == 1) {

        /* Compute temperature seen */
        cs_real_t tempf;
        if (extra->temperature != nullptr) {
          if (t_scl == CS_TEMPERATURE_SCALE_KELVIN)
            tempf = extra->temperature->val[cell_id];
          else /* if (t_scl == CS_TEMPERATURE_SCALE_CELSIUS) */
            tempf = extra->temperature->val[cell_id] + tkelvi;
        }

        else {
          tempf = cs_glob_fluid_properties->t0;
          if (t_scl == CS_TEMPERATURE_SCALE_CELSIUS)
            tempf += tkelvi;
        }

      cs_real_t p_mass = cs_lagr_particle_get_real(particle, p_am,
                                                   CS_LAGR_MASS);

      cs_real_t ddbr = sqrt(2.0 * _k_boltz * tempf / (p_mass * taup_r[0][id]));

      cs_real_t tix2 =   cs_math_pow2((taup_r[0][id] * ddbr))
        * (dt_part - taup_r[0][id] * (1.0 - aux1) * (3.0 - aux1) / 2.0);
      cs_real_t tiu2 =   ddbr * ddbr * taup_r[0][id]
                       * (1.0 - exp(-2.0 * dt_part / taup_r[0][id])) / 2.0;

      cs_real_t tixiu  =
        cs_math_pow2((ddbr * taup_r[0][id] * (1.0 - aux1))) / 2.0;

      tbrix2 = tix2 - (tixiu * tixiu) / tiu2;

      if (tbrix2 > 0.0)
        tbrix2    = sqrt(tbrix2) * brgaus[id];
      else
        tbrix2    = 0.0;

      if (tiu2 > 0.0)
        tbrix1    = tixiu / sqrt(tiu2) * brgaus[3];
      else
        tbrix1    = 0.0;

      if (tiu2 > 0.0) {
        tbriu      = sqrt(tiu2) * brgaus[id + 3];
        if (cs_glob_lagr_time_scheme->t_order == 2)
          cs_lagr_particles_set_real(p_set, p_id, CS_LAGR_BROWN_STATE_1,
                                     sqrt(tiu2));

      }
      else {
        tbriu     = 0.0;
        if (cs_glob_lagr_time_scheme->t_order == 2)
          cs_lagr_particles_set_real(p_set, p_id, CS_LAGR_BROWN_STATE_1, 0.);
      }
    }
    else {
      tbrix1  = 0.0;
      tbrix2  = 0.0;
      tbriu   = 0.0;
    }

    /* Finalize increments */

    /* trajectory  */
    /* Initialized with the four terms not depending on continuous phases */
    displ_r[id] = ter1x + ter4x + ter5x + ter7x + tbrix1 + tbrix2;

    /* particles velocity */
    /* Initialized with the four terms not depending on phases */
    part_vel_r[id] = ter1p + ter4p + ter5p + ter7p + tbriu;

    for (int phase_id = 0; phase_id < n_phases; phase_id++){
      /* trajectory  */
      displ_r[id] += ter2x[phase_id] + ter3x[phase_id] + ter6x[phase_id];

      /* flow-seen velocity */
      part_vel_seen_r[phase_id][id]
        = ter1f[phase_id] + ter2f[phase_id] + ter3f[phase_id] + ter7f;
      if (cs_glob_lagr_model->cs_used == 0)
        part_vel_seen_r[phase_id][id] += fluid_vel_r[phase_id][id];

      /* particles velocity */
      part_vel_r[id] += ter2p[phase_id] + ter3p[phase_id] + ter6p[phase_id];
    }
  }

  /* Reference frame change: --> back to global reference frame
   * NB: Inverse transformation: transpose of trans_m
   * ================================================ */

  if (local_reference_frame) {

    /* Displacement */
    cs_real_t displ[3];
    cs_math_33t_3_product(trans_m, displ_r, displ);
    for (cs_lnum_t id = 0; id < 3; id++)
      part_coords[id] = old_part_coords[id] + displ[id];

    /* Particle velocity */
    cs_math_33t_3_product(trans_m, part_vel_r, part_vel);

    /* Flow-seen velocity */
    for (int phase_id = 0; phase_id < n_phases; phase_id++){
      cs_real_t part_vel_seen_p[3] = {0., 0., 0.};
      cs_math_33t_3_product(trans_m, part_vel_seen_r[phase_id], part_vel_seen_p);
      for (int id = 0; id < 3; id++) {
        if (lambda[phase_id] < 1.e-8) {
          part_vel_seen[phase_id * 3 + id] = cvar_vel[phase_id][cell_id][id];
        } else {
          part_vel_seen[phase_id * 3 + id] = part_vel_seen_p[id];
        }
      }
    }
  }

  else { /* local_reference_frame == false */
    for (cs_lnum_t id = 0; id < 3; id++) {
      part_coords[id] = old_part_coords[id] + displ_r[id];
      part_vel[id] = part_vel_r[id];
      for (int phase_id = 0; phase_id < n_phases; phase_id ++){
        if (lambda[phase_id] < 1.0e-8){
          part_vel_seen[phase_id * 3 + id] = cvar_vel[phase_id][cell_id][id];
        } else {
          part_vel_seen[phase_id * 3 + id] = part_vel_seen_r[phase_id][id];
        }
      }
    }
  }

  /* Save the covariance field for model determination in neptune_cfd */
  if (   cs_glob_time_step->nt_cur > cs_glob_lagr_stat_options->idstnt
      && cs_glob_lagr_model->cs_used == 0) {
    cs_real_t *vel_seen_vel_cov
      = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
          CS_LAGR_VELOCITY_SEEN_VELOCITY_COV);
    for (int phase_id = 0; phase_id < n_phases; phase_id++) {
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          vel_seen_vel_cov[9 * phase_id + (3 * i + j)]
            = (part_vel_seen[phase_id * 3 + i] - cvar_vel[phase_id][cell_id][i])
            * (part_vel[j] - stat_vel->val[cell_id * 3 + j]);
        }
      }
    }
  }
  CS_FREE(lambda);
  CS_FREE(loc_fluid_vel);
  CS_FREE(part_vel_seen_r);
  CS_FREE(old_part_vel_seen_r);
  CS_FREE(fluid_vel_r);
  CS_FREE(piil_r);
  CS_FREE(tlag_r);
  CS_FREE(taup_r);
  CS_FREE(cvar_vel);
}

/*----------------------------------------------------------------------------*/
/*! \brief Integration of SDEs by 2nd order scheme
 *
 * When there has beed interaction with a boundary face, the velocity and
 * velocity seen computations are forced to 1st order.
 *
 * \param[in]  p_id      particle index in set
 * \param[in]  dt_part   remaining time step associated to the particle
 * \param[in]  nor       current step id (for 2nd order scheme)
 * \param[in]  taup      dynamic characteristic time
 * \param[in]  tlag      lagrangian fluid characteristic time
 * \param[in]  piil      term in integration of up sdes
 * \param[in]  bx        turbulence characteristics
 * \param[out] tsfext    info for return coupling source terms
 * \param[in]  vagaus    gaussian random variables
 * \param[in]  brgaus    gaussian random variables
 * \param[in]  force_p   forces per mass unit on particles (m/s^2)
 * \param[in]  beta      proportional to the gradient of T_lag
 */
/*----------------------------------------------------------------------------*/

static void
cs_sde_vels_pos_2_nd_order_time_integ(cs_lnum_t                       p_id,
                                      cs_real_t                       dt_part,
                                      int                             nor,
                                      const cs_real_t                *taup,
                                      const cs_real_3_t              *tlag,
                                      const cs_real_3_t              *piil,
                                      const cs_real_33_t             *bx,
                                      cs_real_t                      *tsfext,
                                      const cs_real_3_t              *vagaus,
                                      const cs_real_6_t               brgaus,
                                      const cs_real_3_t               force_p,
                                      const cs_real_3_t               beta)
{
  cs_real_t  aux0, aux1, aux2, aux3, aux4, aux5, aux6, aux7, aux8, aux9;
  cs_real_t  aux10, aux11, aux12, aux17, aux18, aux19;
  cs_real_t  aux20;
  cs_real_t  ter1, ter2, ter3, ter4, ter5;
  cs_real_t  sige, tapn, gamma2;
  cs_real_t  grgam2, gagam;
  cs_real_t  p11, p21, p22;
  cs_real_t  tbriu = 0.;

  /* Particles management */
  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;
  unsigned char *particle = p_set->p_buffer + p_am->extents * p_id;

  cs_lagr_extra_module_t *extra_i = cs_get_lagr_extra_module();
  cs_lagr_extra_module_t *extra = extra_i;

  /* Initializations
     ===========================================================================*/

  /* Integration of the SDE on partilces
     ===========================================================================*/

  /* Computation each sub-time step
     ---------------------------------------------------------------------------*/
  /* Compute tau_p*A_p and II*TL+<u>:
   * ------------------------------- */
  cs_real_6_t auxl;

  /* use previous step for t_order == 1 or prediction step
   * and current one for correction step */
  cs_lnum_t cell_id = cs_lagr_particle_get_lnum_n(particle, p_set->p_am, 2-nor,
                                                  CS_LAGR_CELL_ID);

  for (cs_lnum_t id = 0; id < 3; id++) {
    auxl[id] = force_p[id];

    if (nor == 1)
      auxl[id + 3] =   piil[0][id] * tlag[0][id]
                     + extra->vel->vals[1][cell_id * 3 + id];
    else
      auxl[id + 3] =   piil[0][id] * tlag[0][id]
                     + extra->vel->vals[0][cell_id * 3 + id];

  }

  /* Prediction step
  ===========================================================================*/

  if (nor == 1) {

    /* Save tau_p^n */
    cs_lagr_particle_set_real(particle, p_am, CS_LAGR_TAUP_AUX, taup[0]);

    /* Save coupling */
    if (cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING) {

      aux0     = -dt_part / taup[0];
      aux1     =  exp(aux0);
      *tsfext  =   taup[0]
                   * cs_lagr_particles_get_real(p_set, p_id, CS_LAGR_MASS)
                   * (-aux1 + (aux1 - 1.0) / aux0);
    }

    /* Load terms at t = t_n : */
    auto *old_part_vel      =
      cs_lagr_particle_attr_n_get_ptr<cs_real_t>(particle, p_am,
                                                 1, CS_LAGR_VELOCITY);
    auto *old_part_vel_seen =
      cs_lagr_particle_attr_n_get_ptr<cs_real_t>(particle, p_am,
                                                 1, CS_LAGR_VELOCITY_SEEN);
    auto *pred_part_vel_seen =
      cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                               CS_LAGR_PRED_VELOCITY_SEEN);
    auto *pred_part_vel =
      cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                               CS_LAGR_PRED_VELOCITY);

    for (cs_lnum_t id = 0; id < 3; id++) {

      aux0    =  -dt_part / taup[0];
      aux1    =  -dt_part / tlag[0][id];
      aux2    = exp(aux0);
      aux3    = exp(aux1);
      aux4    = tlag[0][id] / (tlag[0][id] - taup[0]);
      aux5    = aux3 - aux2;

      pred_part_vel_seen[id] =   0.5 * old_part_vel_seen[id]
                               * aux3 + auxl[id + 3]
                               * (-aux3 + (aux3 - 1.0) / aux1);

      ter1    = 0.5 * old_part_vel[id] * aux2;
      ter2    = 0.5 * old_part_vel_seen[id] * aux4 * aux5;
      ter3    =   auxl[id + 3]
                * (  -aux2 + ((tlag[0][id] + taup[0]) / dt_part) * (1.0 - aux2)
                   - (1.0 + tlag[0][id] / dt_part) * aux4 * aux5);
      ter4    = auxl[id] * (-aux2 + (aux2 - 1.0) / aux0);
      pred_part_vel[id] = ter1 + ter2 + ter3 + ter4;

    }

    /* Euler scheme */

    cs_sde_vels_pos_1_st_order_time_integ(p_id,
                                          dt_part,
                                          nor,
                                          taup,
                                          tlag,
                                          piil,
                                          bx,
                                          vagaus,
                                          brgaus,
                                          force_p,
                                          beta);
  }
  else {

    /* Correction stage  (not called if a specific treatment at a face has
     *                    occured during previous trajectography step)
     * ----------------*/

    /* Compute Us */
    /* FIXME bugged? should not be CS_LAGR_REBOUND_ID != 0 */
    if (   !cs_lagr_particles_get_flag(p_set, p_id, CS_LAGR_PART_IMPOSED_MOTION)
        && cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_REBOUND_ID) == 0) {
      cs_real_t *part_vel
        = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                   CS_LAGR_VELOCITY);
      cs_real_t *part_vel_seen
        = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                   CS_LAGR_VELOCITY_SEEN);
      cs_real_t *old_part_vel
        = cs_lagr_particle_attr_n_get_ptr<cs_real_t>(particle, p_am,
                                                     1, CS_LAGR_VELOCITY);
      cs_real_t *old_part_vel_seen
        = cs_lagr_particle_attr_n_get_ptr<cs_real_t>(particle, p_am,
                                                     1, CS_LAGR_VELOCITY_SEEN);
      cs_real_t *pred_part_vel_seen
        = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                   CS_LAGR_PRED_VELOCITY_SEEN);
      cs_real_t *pred_part_vel
        = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                   CS_LAGR_PRED_VELOCITY);

      for (cs_lnum_t id = 0; id < 3; id++) {

        aux0    =  -dt_part / taup[0];
        aux1    =  -dt_part / tlag[0][id];
        aux2    = exp(aux0);
        aux3    = exp(aux1);
        aux4    = tlag[0][id] / (tlag[0][id] - taup[0]);
        aux5    = aux3 - aux2;
        aux6    = aux3 * aux3;

        ter1    = 0.5 * old_part_vel_seen[id] * aux3;
        ter2    = auxl[id + 3] * (1.0 - (aux3 - 1.0) / aux1);
        ter3    =  -aux6 + (aux6 - 1.0) / (2.0 * aux1);
        ter4    = 1.0 - (aux6 - 1.0) / (2.0 * aux1);

        sige    =   (  ter3 * bx[0][id][0]
                     + ter4 * bx[0][id][1] )
                  * (1.0 / (1.0 - aux6));

        ter5    = 0.5 * tlag[0][id] * (1.0 - aux6);

        part_vel_seen[id] =   pred_part_vel_seen[id] + ter1 + ter2
                            + sige * sqrt(ter5) * vagaus[0][id];

        /* compute Up */
        ter1    = 0.5 * old_part_vel[id] * aux2;
        ter2    = 0.5 * old_part_vel_seen[id] * aux4 * aux5;
        ter3    =   auxl[id + 3]
                * (   1.0 - ((tlag[0][id] + taup[0]) / dt_part) * (1.0 - aux2)
                    + (tlag[0][id] / dt_part) * aux4 * aux5)
          + auxl[id] * (1.0 - (aux2 - 1.0) / aux0);

        tapn    = cs_lagr_particles_get_real(p_set, p_id, CS_LAGR_TAUP_AUX);

        aux7    = exp(-dt_part / tapn);
        aux8    = 1.0 - aux3 * aux7;
        aux9    = 1.0 - aux6;
        aux10   = 1.0 - aux7 * aux7;
        aux11   = tapn / (tlag[0][id] + tapn);
        aux12   = tlag[0][id] / (tlag[0][id] - tapn);
        aux17   = sige * sige * aux12 * aux12;
        aux18   = 0.5 * tlag[0][id] * aux9;
        aux19   = 0.5 * tapn * aux10;
        aux20   = tlag[0][id] * aux11 * aux8;

        /* compute correlation matrix */
        gamma2  = sige * sige * aux18;
        grgam2  = aux17 * (aux18 - 2.0 * aux20 + aux19);
        gagam   = sige * sige * aux12 * (aux18 - aux20);

        /* Gaussian vector simulation */

        p11  = sqrt(cs::max(0.0, gamma2));
        if (p11 > cs_math_epzero) {
          p21  = gagam / p11;
          p22  = grgam2 - p21 * p21;
          p22  = sqrt(cs::max(0.0, p22));
        }
        else {
          p21  = 0.0;
          p22  = 0.0;
        }

        ter4    = p21 * vagaus[0][id] + p22 * vagaus[1][id];

        /* Compute terms in Brownian movement */
        if (cs_glob_lagr_brownian->lamvbr == 1)
          tbriu = brgaus[id + 3]
            * cs_lagr_particles_get_real(p_set, p_id, CS_LAGR_BROWN_STATE_1);
        else
          tbriu = 0.0;

        /* finalize writing */

        part_vel[id] = pred_part_vel[id] + ter1 + ter2 + ter3 + ter4 + tbriu;

      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*! \brief Deposition submodel
 *
 *  1/ Modification of the coordinate system (global ->local)
 *  2/ Call of subroutine lagcli
 *  3/ Integration of the stochastic differential equations
 *     in the 2 directions different from the normal to the boundary face
 *  4/ Modification of the coordinate system (local ->global)
 *  5/ Update of the particle position
 *
 * \param[in]  dt_part   remaining time step associated to the particle
 * \param[in]  p_id      particle index in set
 * \param[in]  nor       current step id (for 2nd order scheme)
 * \param[in]  taup      dynamic characteristic time
 * \param[in]  piil      term in integration of UP SDEs
 * \param[in]  vagaus    gaussian random variables
 * \param[in]  romp      particles associated density
 * \param[in]  force_p   forces per mass unit on particles (m/s^2)
 * \param[in]  tempf     temperature of the fluid (K)
 * \param[in]  vislen    FIXME
 * \param[in]  events    associated events set
 * \param[in]  depint    interface location near-wall/core-flow
 * \param[in]  nresnew
 *
 */
/*----------------------------------------------------------------------------*/

static void
_lagesd(cs_real_t                       dt_part,
        cs_lnum_t                       p_id,
        int                             nor,
        const cs_real_t                 taup,
        const cs_real_3_t               piil,
        const cs_real_3_t              *vagaus,
        const cs_real_t                 romp,
        const cs_real_3_t               force_p,
        cs_real_t                       tempf,
        const cs_real_t                 vislen[],
        cs_lagr_event_set_t            *events,
        cs_real_t                      *depint,
        cs_lnum_t                      *nresnew)
{
  /* mesh and mesh quantities */
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  /* Particles management */
  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;

  cs_lagr_extra_module_t *extra_i = cs_get_lagr_extra_module();
  cs_lagr_extra_module_t *extra = extra_i;

  cs_lagr_boundary_interactions_t  *bi = cs_glob_lagr_boundary_interactions;

  /* Hydrodynamic drag and torque on a deposited particle     */
  cs_real_t drag_force[3];
  cs_real_t drag_torque[3];

  /* Lift force and torque on a deposited particle   */
  cs_real_t lift_force[1];
  cs_real_t lift_torque[3];

  /* Gravity force and torque on a deposited particle   */
  cs_real_t grav_force[3];
  cs_real_t grav_torque[3];

  /* Adhesion force and torque on a deposited particle   */
  cs_real_t adhes_torque[3];

  /* Map field arrays */
  const int _prev_id = (extra->vel->n_time_vals > 1) ? 1 : 0;
  const cs_real_3_t *cvar_vel
    = (const cs_real_3_t *)(extra->vel->vals[_prev_id]);

  /* ===========================================================================
   * 1. INITIALISATIONS
   * ======================================================================== */

  const cs_real_t  *grav  = cs_glob_physical_constants->gravity;

  /* particle data */
  unsigned char *particle = p_set->p_buffer + p_am->extents * p_id;

  cs_real_t p_mass = cs_lagr_particle_get_real(particle, p_am,
                                                CS_LAGR_MASS);
  cs_real_t p_diam = cs_lagr_particle_get_real(particle, p_am,
                                                CS_LAGR_DIAMETER);
  cs_real_t p_stat_w = cs_lagr_particle_get_real(particle, p_am,
                                                  CS_LAGR_STAT_WEIGHT);

  /* use previous step for t_order == 1 or prediction step
   * and current one for correction step */
  cs_lnum_t cell_id = cs_lagr_particle_get_lnum_n(particle, p_set->p_am, 2-nor,
                                                  CS_LAGR_CELL_ID);

  cs_lnum_t face_id = cs_lagr_particle_get_lnum(particle, p_am,
                                                 CS_LAGR_NEIGHBOR_FACE_ID);

  assert(face_id > -1);

  cs_real_t ustar = extra->ustar->val[face_id];
  cs_real_t lvisq = vislen[face_id];

  cs_real_t tvisq;
  if (ustar > 0.0)
    tvisq = lvisq / ustar;
  else
    tvisq = cs_math_big_r;

  /* Constants for the calculation of bxp and tlp  */
  cs_real_t c0 = cs_turb_crij_c0;

  cs_real_t cl = 1.0 / (0.5 + 0.75 * c0);

  /* Pointer on the density w.r.t the flow    */
  cs_real_t romf = extra->cromf->val[cell_id];

  cs_real_t visccf = extra->viscl->val[cell_id] / romf;

  cs_real_t yplus = cs_lagr_particle_get_real(particle, p_am,
                                              CS_LAGR_YPLUS);

  /* Turbulent kinetic energy and dissipation w.r.t y+  */
  cs_real_t energi, dissip;
  if (yplus <= 5.0) {

    energi = 0.1 * cs_math_pow2(yplus) * cs_math_pow2(ustar);
    dissip = 0.2 * cs_math_pow4(ustar) / visccf;

  }

  else if (yplus <= 30.0) {

    energi = cs_math_pow2(ustar) / sqrt(0.09);
    dissip = 0.2 * cs_math_pow4(ustar) / visccf;

  }
  else {

    assert(yplus <= 100.0); /* should not arrive here otherwise */

    energi   = cs_math_pow2(ustar) / sqrt(0.09);
    dissip   = cs_math_pow4(ustar) / (0.41 * yplus * visccf);

  }

  /* ===========================================================================
   * 2. Reference frame change:
   * --------------------------
   * global reference frame --> local reference frame for the boundary face
   * ======================================================================== */

  const cs_real_3_t *rot_m
    = (const cs_real_3_t *)cs_glob_lagr_b_face_proj[face_id];

  /* 2.1 - particle velocity   */

  cs_real_t *old_part_vel =
    cs_lagr_particle_attr_n_get_ptr<cs_real_t>(particle, p_am, 1,
                                               CS_LAGR_VELOCITY);
  cs_real_t vpart[3];

  cs_math_33_3_product(rot_m, old_part_vel, vpart);

  cs_real_t vpart0[3] = {vpart[0], vpart[1], vpart[2]};

  /* 2.2 - flow-seen velocity  */

  cs_real_t *old_part_vel_seen =
    cs_lagr_particle_attr_n_get_ptr<cs_real_t>(particle, p_am, 1,
                                               CS_LAGR_VELOCITY_SEEN);
  cs_real_t vvue[3];

  cs_math_33_3_product(rot_m, old_part_vel_seen, vvue);

  /* 2.3 - Gravity vector */

  cs_real_t ggp[3];

  cs_math_33_3_product(rot_m, grav, ggp);

  /* 2.4 - flow velocity  */

  cs_real_t vflui[3];

  cs_math_33_3_product(rot_m, cvar_vel[cell_id], vflui);

  cs_real_t norm = sqrt(cs_math_pow2(vflui[1]) + cs_math_pow2(vflui[2]));
  cs_real_t inv_norm = ((norm > cs_math_zero_threshold) ?  1. / norm : 0);

  /* Velocity norm w.r.t y+ */
  cs_real_t norm_vit;

  if (yplus <= 5.0)
    norm_vit = yplus * ustar;

  else if (yplus <= 30.0)
    norm_vit = ( -3.05 + 5.0 * log(yplus)) * ustar;

  else {
    assert(yplus < 100.0);
    norm_vit = (2.5 * log (yplus) + 5.5) * ustar;
  }

  vflui[1] = norm_vit * vflui[1] * inv_norm;
  vflui[2] = norm_vit * vflui[2] * inv_norm;

  /* 2.5 Particle force: - pressure gradient/romp + external force + g   */

  cs_real_t force_p_r[3];

  cs_math_33_3_product(rot_m, force_p, force_p_r);

  /* 2.6 - "piil" term    */

  cs_real_t piilp[3];

  cs_math_33_3_product(rot_m, piil, piilp);

  /* 2.7 - tlag */

  cs_real_t tlp = cs_math_epzero;
  if (dissip > cs_math_zero_threshold) {

    tlp = cl * energi / dissip;
    tlp = cs::max(tlp, cs_math_epzero);

  }

  /* 2.8 - bx   */
  cs_real_t bxp = sqrt(c0 * dissip);

  /* =========================================================================
   * 3. Integration of the EDS on the particles
   * =========================================================================*/

  /* Retrieve of the turbulent kinetic energy */
  cs_real_t  enertur;
  if (extra->itytur == 2 || extra->itytur == 4 ||
      extra->iturb == 50 || extra->iturb == 60)
    enertur  = extra->cvar_k->vals[_prev_id][cell_id];

  else if (extra->itytur == 3) {
    enertur  = 0.5 * (  extra->cvar_rij->vals[_prev_id][6*cell_id    ]
                      + extra->cvar_rij->vals[_prev_id][6*cell_id + 1]
                      + extra->cvar_rij->vals[_prev_id][6*cell_id + 2]);
  }

  cs_lnum_t marko  = cs_lagr_particle_get_lnum(particle, p_am,
                                                CS_LAGR_MARKO_VALUE);
  cs_real_t interf = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_INTERF);
  cs_real_t depl[3];

  cs_lagr_deposition(dt_part,
                     &marko,
                     tempf,
                     lvisq,
                     tvisq,
                     vpart,
                     vvue,
                     depl,
                     &p_diam,
                     romp,
                     taup,
                     &yplus,
                     &interf,
                     &enertur,
                     ggp,
                     vflui,
                     force_p_r,
                     piilp,
                     depint);

  cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_MARKO_VALUE, marko);

  if (cs_lagr_particles_get_flag(p_set, p_id, CS_LAGR_PART_DEPOSITION_FLAGS)) {

    depl[0]  = 0.0;
    vpart[0] = 0.0;

    /* Integration in the 2 other directions */

    for (cs_lnum_t id = 1; id < 3; id++) { // FIXME seems strange

      cs_lnum_t i0 = id-1;

      cs_real_t tci   = piilp[id] * tlp + vflui[id];
      cs_real_t aux2  = exp(-dt_part / tlp);
      cs_real_t aux6  = bxp * bxp * tlp;

      /* --> Terms for the flow-seen velocity     */
      cs_real_t ter1f = vvue[id] * aux2;
      cs_real_t ter2f = tci * (1.0 - aux2);

      /* --> (2.3) Coefficients computation for the stochastic integrals:  */
      cs_real_t gama2 = 0.5 * (1.0 - aux2 * aux2);

      /* --> Integral for the flow-seen velocity  */
      cs_real_t p11   = sqrt (gama2 * aux6);
      cs_real_t ter3f = p11 * vagaus[0][i0];

      /* --> flow-seen velocity    */
      vvue[id] = ter1f + ter2f + ter3f;

    }

  }

  else { /* Particle in flow */

    for (cs_lnum_t id = 1; id < 3; id++) {

      cs_lnum_t i0 = id-1; //FIXME strange

      cs_real_t tci   = piilp[id] * tlp + vflui[id];
      cs_real_t v_lim = force_p_r[id] * taup;
      cs_real_t aux1  = exp(-dt_part / taup);
      cs_real_t aux2  = exp(-dt_part / tlp);
      cs_real_t aux3  = tlp / (tlp - taup);
      cs_real_t aux4  = tlp / (tlp + taup);
      cs_real_t aux5  = tlp * (1.0 - aux2);
      cs_real_t aux6  = bxp * bxp * tlp;
      cs_real_t aux7  = tlp - taup;
      cs_real_t aux8  = bxp * bxp * cs_math_pow2(aux3);

      /* --> Terms for the trajectory   */
      cs_real_t aa    = taup * (1.0 - aux1);
      cs_real_t bb    = (aux5 - aa) * aux3;
      cs_real_t cc    = dt_part - aa - bb;
      cs_real_t ter1x = aa * vpart[id];
      cs_real_t ter2x = bb * vvue[id];
      cs_real_t ter3x = cc * tci;
      cs_real_t ter4x = (dt_part - aa) * v_lim;

      /* --> Terms for the flow-seen velocity     */
      cs_real_t ter1f = vvue[id] * aux2;
      cs_real_t ter2f = tci * (1.0 - aux2);

      /* --> Terms for the particles velocity     */
      cs_real_t dd    = aux3 * (aux2 - aux1);
      cs_real_t ee    = 1.0 - aux1;
      cs_real_t ter1p = vpart[id] * aux1;
      cs_real_t ter2p = vvue[id] * dd;
      cs_real_t ter3p = tci * (ee - dd);
      cs_real_t ter4p = v_lim * ee;

      /* --> (2.3) Coefficients computation for the stochastic integrals:  */
      cs_real_t gama2  = 0.5 * (1.0 - aux2 * aux2);
      cs_real_t omegam = aux3 * ( (tlp - taup) * (1.0 - aux2)
                                  - 0.5 * tlp * (1.0 - aux2 * aux2)
                                  + cs_math_pow2(taup) / (tlp + taup)
                                  * (1.0 - aux1 * aux2)
                                  ) * aux6;
      cs_real_t omega2 =   aux7
                         * (aux7 * dt_part - 2.0 * (tlp * aux5 - taup * aa))
                         + 0.5 * tlp * tlp * aux5 * (1.0 + aux2)
                         + 0.5 * cs_math_pow2(taup) * aa * (1.0 + aux1)
                         - 2.0 * aux4 * tlp * cs_math_pow2(taup) * (1.0 - aux1 * aux2);
      omega2 *= aux8;

      cs_real_t  p11, p21, p22, p31, p32, p33;

      if (cs::abs(gama2) >cs_math_epzero) {

        p21    = omegam / sqrt (gama2);
        p22    = omega2 - cs_math_pow2(p21);
        p22    = sqrt(cs::max(0.0, p22));

      }
      else {

        p21    = 0.0;
        p22    = 0.0;

      }

      cs_real_t ter5x = p21 * vagaus[0][i0] + p22 * vagaus[1][i0];

      /* --> Integral for the flow-seen velocity  */

      p11   = sqrt(gama2 * aux6);
      cs_real_t ter3f = p11 * vagaus[0][i0];

      /* --> Integral for particles velocity */

      cs_real_t aux9  = 0.5 * tlp * (1.0 - aux2 * aux2);
      cs_real_t aux10 = 0.5 * taup * (1.0 - aux1 * aux1);
      cs_real_t aux11 = taup * tlp * (1.0 - aux1 * aux2) / (taup + tlp);
      cs_real_t grga2 = (aux9 - 2.0 * aux11 + aux10) * aux8;
      cs_real_t gagam = (aux9 - aux11) * (aux8 / aux3);
      cs_real_t gaome = (  (tlp - taup) * (aux5 - aa)
                         - tlp * aux9
                         - taup * aux10
                         + (tlp + taup) * aux11) * aux8;

      if (p11 > cs_math_epzero)
        p31 = gagam / p11;

      else
        p31 = 0.0;

      if (p22 > cs_math_epzero)
        p32 = (gaome - p31 * p21) / p22;

      else
        p32 = 0.0;

      p33 = grga2 - cs_math_pow2(p31) - cs_math_pow2(p32);
      p33 = sqrt (cs::max(0.0, p33));

      cs_real_t ter5p =   p31 * vagaus[0][i0]
                        + p32 * vagaus[1][i0]
                        + p33 * vagaus[2][i0];

      /* --> trajectory  */
      depl[id] = ter1x + ter2x + ter3x + ter4x + ter5x;

      /* --> flow-seen velocity    */
      vvue[id] = ter1f + ter2f + ter3f;

      /* --> particles velocity    */
      vpart[id] = ter1p + ter2p + ter3p + ter4p + ter5p;

    }

  }

  if (cs_glob_lagr_model->resuspension == 1) {

    cs_real_t p_height = cs_lagr_particle_get_real(particle, p_am,
                                                   CS_LAGR_HEIGHT);

    cs_lnum_t iresusp = 0;

    if (cs_lagr_particles_get_flag(p_set, p_id, CS_LAGR_PART_DEPOSITION_FLAGS)) {

      cs_lnum_t n_f_id
        = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_NEIGHBOR_FACE_ID);

      cs_lnum_t nfabor = cs_glob_mesh->n_b_faces;

      /* Resuspension model
       * differentiation between monolayer and multilayer resuspension
       * based on the mean deposit height
       * (another possibility is to use the jamming limit...) */
      cs_real_t diam_mean = cs_glob_lagr_clogging_model->diam_mean;

      if (   (cs_glob_lagr_model->clogging == 0 && iresusp == 0)
          || (   cs_glob_lagr_model->clogging == 1
              &&   bound_stat[n_f_id + nfabor * bi->ihdepm]
                 < diam_mean && iresusp == 0)) {

        /* Monolayer resuspension */

        /* Calculation of the hydrodynamic drag and torque
         * applied on the deposited particle   */

        drag_force[0] =  3.0 * cs_math_pi * p_diam
          * (vvue[0] - vpart[0]) * visccf * romf * 3.39;
        drag_torque[0] = 0.0;

        for (cs_lnum_t id = 1; id < 3; id++) {

          drag_force[id]   =  3.0 * cs_math_pi * p_diam
            * (vvue[id] - vpart[id]) * visccf * romf * 1.7;
          drag_torque[id] = 1.4 * drag_force[id] * p_diam * 0.5;

        }

        /* Calculation of lift force and torque */
        lift_force[0] = - 20.0 * cs_math_pow2(visccf) * romf *
          pow(ustar * p_diam * 0.5 / visccf, 2.31);

        for (cs_lnum_t id = 1; id < 3; id++) {

          lift_torque[id] = -lift_force[0] * p_diam * 0.5
            * vvue[id] / sqrt(cs_math_pow2(vvue[1]) + cs_math_pow2(vvue[2]) );

        }

        /* Calculation of gravity force and torque */
        for (cs_lnum_t id = 0; id < 3; id++) {
          grav_force[id] = 4./3. * cs_math_pi * cs_math_pow3(p_diam/2)
                                 * romp * ggp[id];
        }

        for (cs_lnum_t id = 1; id < 3; id++) {

          cs_real_t sign11, sign22;
          if (grav_force[1]*vvue[1] < 0 )
            sign11 = -1.0;
          else
            sign11 = 1.0;
          if (grav_force[2]*vvue[2] < 0 )
            sign22 = -1.0;
          else
            sign22 = 1.0;

          grav_torque[id] = (-grav_force[0] +
                              grav_force[1]*sign11 +
                              grav_force[2]*sign22 ) * p_diam * 0.5
            * vvue[id] / sqrt(cs_math_pow2(vvue[1]) + cs_math_pow2(vvue[2]));

        }

        /* Case with consolidation */
        if (    cs_glob_lagr_consolidation_model->iconsol > 0
             &&   cs_lagr_particle_get_real(particle, p_am, CS_LAGR_CONSOL_HEIGHT)
                > 0.01 * diam_mean) {

          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_ADHESION_FORCE,
                                    cs_glob_lagr_consolidation_model->force_consol);
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_ADHESION_TORQUE,
                                    cs_glob_lagr_consolidation_model->force_consol
                                    * p_diam * 0.5);

        }

        cs_real_t adhes_force =
          cs_lagr_particle_get_real(particle, p_am, CS_LAGR_ADHESION_FORCE);

        /* Is there direct wall-normal lift-off of the particle ? */
        if (   (adhes_force + grav_force[0] + lift_force[0] + drag_force[0]) < 0
            && iresusp == 0 ) {

          /* Update of the number and weight of resuspended particles     */
          p_set->n_part_resusp += 1;
          p_set->weight_resusp += p_stat_w;

          if (events != nullptr) {
            const cs_real_t *part_vel
              = cs_lagr_particles_attr_get_const_ptr<cs_real_t>(p_set, p_id,
                                                                CS_LAGR_VELOCITY);
            _add_resuspension_event(events,
                                    p_set,
                                    p_id,
                                    face_id,
                                    part_vel);
          }

          /* Update of surface covered and deposit height
           * (for non-rolling particles) */
          if (cs_lagr_particles_get_flag(p_set, p_id, CS_LAGR_PART_DEPOSITED)) {

            bound_stat[n_f_id + nfabor * bi->iscovc]
              -=   cs_math_pi * cs_math_pow2(p_diam) * p_stat_w
                 * 0.25 / mq->b_face_surf[n_f_id];

            bound_stat[n_f_id + nfabor * bi->ihdepm]
              -=   cs_math_pi * p_height * cs_math_pow2(p_diam) * p_stat_w
                 * 0.25 / mq->b_face_surf[n_f_id];

            bound_stat[n_f_id + nfabor * bi->ihdepv]
              -=   cs_math_pow2(cs_math_pi * p_height * cs_math_pow2(p_diam) * p_stat_w
                 * 0.25 / mq->b_face_surf[n_f_id]);

            bound_stat[n_f_id + nfabor * bi->inclg]
              -= p_stat_w;

          }

          /* The particle is resuspended  */
          vpart[0] = cs::min(-1.0 / p_mass * dt_part
                             * cs::abs(drag_force[0] - adhes_force),
                             0.001);
          if (cs_lagr_particles_get_flag(p_set, p_id, CS_LAGR_PART_DEPOSITED)) {
            vpart[1] = 0.0;
            vpart[2] = 0.0;
          }

          cs_lagr_particles_unset_flag(p_set, p_id, CS_LAGR_PART_DEPOSITION_FLAGS);
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_ADHESION_FORCE, 0.0);
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_ADHESION_TORQUE, 0.0);

          if (cs_glob_lagr_model->clogging == 1) {
            if (cs::abs(p_height-p_diam)/p_diam > 1.0e-6) {
              cs_real_t d_resusp = pow(0.75 * cs_math_pow2(p_diam) * p_height, 1.0/3.0);
              cs_lagr_particle_set_real(particle, p_am, CS_LAGR_DIAMETER, d_resusp);
              cs_lagr_particle_set_real(particle, p_am, CS_LAGR_HEIGHT, d_resusp);
            }
          }

          if (p_am->count[0][CS_LAGR_N_LARGE_ASPERITIES] > 0)
            cs_lagr_particle_set_lnum(particle, p_am,
                                      CS_LAGR_N_LARGE_ASPERITIES, 0);

          if (p_am->count[0][CS_LAGR_N_SMALL_ASPERITIES] > 0)
            cs_lagr_particle_set_lnum(particle, p_am,
                                      CS_LAGR_N_SMALL_ASPERITIES, 0);

          if (p_am->count[0][CS_LAGR_DISPLACEMENT_NORM] > 0)
            cs_lagr_particle_set_real(particle, p_am,
                                      CS_LAGR_DISPLACEMENT_NORM, 0.0);

          iresusp = 1;
        }
        /* No direct normal lift-off */
        else if (iresusp == 0) {

          /* Calculation of the norm of the hydrodynamic
           * torque and drag (tangential) */

          for (cs_lnum_t id = 1; id < 3; id++) {

            adhes_torque[id]  = - cs_lagr_particle_get_real(particle, p_am,
                                                            CS_LAGR_ADHESION_TORQUE)
              * vvue[id] / sqrt(cs_math_pow2(vvue[1]) + cs_math_pow2(vvue[2]));
          }

          cs_real_t iner_tor = (7.0 / 5.0) * p_mass * cs_math_pow2((p_diam * 0.5));
          cs_real_t cst_1, cst_4;
          cst_4 =   6 * cs_math_pi * visccf
            * romf * 1.7 * 1.4
            * cs_math_pow2(p_diam * 0.5);
          cst_1 = cst_4 * (p_diam * 0.5) / iner_tor;

          for (cs_lnum_t id = 1; id < 3; id++) {

            vpart0[id] = vpart[id];
            vpart[id]  =  vpart0[id] * exp(-cst_1 * dt_part)
                        + (vvue[id] + (adhes_torque[id] + lift_torque[id]
                                       + grav_torque[id]) / cst_4)
                        * (1.0 - exp(-cst_1 * dt_part) );

          }

          cs_real_t scalax  = vpart[1] * vvue[1] + vpart[2] * vvue[2];

          if (scalax > 0.0) {

            /* The calculated particle velocity is aligned   */
            /* with the flow seen   */
            /* --> The particle starts or keep on rolling    */

            /* If the particle stars rolling:
             * update of the surface covered and deposit height */

            if (cs_lagr_particles_get_flag(p_set, p_id, CS_LAGR_PART_DEPOSITED)) {

              bound_stat[n_f_id + nfabor * bi->iscovc]
                -=   cs_math_pi * cs_math_pow2(p_diam) * p_stat_w
                   * 0.25 / mq->b_face_surf[n_f_id];

              bound_stat[n_f_id + nfabor * bi->ihdepm]
                -=   cs_math_pi * p_height * cs_math_pow2(p_diam) * p_stat_w
                   * 0.25 / mq->b_face_surf[n_f_id];

              bound_stat[n_f_id + nfabor * bi->ihdepv]
                -=   cs_math_pow2(cs_math_pi * p_height
                   * cs_math_pow2(p_diam) * p_stat_w
                   * 0.25 / mq->b_face_surf[n_f_id]);

              bound_stat[n_f_id + nfabor * bi->inclg]
                -= p_stat_w;

            }

            cs_lagr_particles_unset_flag(p_set, p_id, CS_LAGR_PART_DEPOSITION_FLAGS);
            cs_lagr_particles_set_flag(p_set, p_id, CS_LAGR_PART_ROLLING);

            vpart[0] = 0.0;

            for (cs_lnum_t id = 1; id < 3; id++) {

              if (cs::abs (vpart[id]) > cs::abs (vvue[id]))
                /* The velocity of the rolling particle cannot   */
                /* exceed the surrounding fluid velocity    */
                vpart[id] = vvue[id];

              cs_real_t kkk = vvue[id] + adhes_torque[id] / cst_4;
              cs_real_t kk  = vpart0[id] - kkk;

              depl[id] = (kkk * dt_part + kk / cst_1 * (1. - exp(-cst_1 * dt_part)));

            }

            iresusp = 1;

          }
          /* if (scalax..) */
          else {

            /* The particle is not set into motion or stops
             * the flag is set to 1 and velocity and displacement are null */

            /* Simplified treatment:
             * no update of iscovc, ihdepm for particles that stop */
            if (cs_lagr_particles_get_flag(p_set, p_id, CS_LAGR_PART_ROLLING))
              iresusp = 1;

            for (cs_lnum_t id = 1; id < 3; id++) {
              depl[id]     = 0.0;
              vpart[id]    = 0.0;
            }

          } /* if (scalax..)   */

        }

      } /* End of monolayer resuspension*/

      else if (   cs_glob_lagr_model->clogging == 1
               &&    bound_stat[n_f_id + nfabor * bi->ihdepm]
                  >= diam_mean
               && iresusp == 0) {

        /* Multilayer resuspension model */

        cs_real_t mean_depo_height
          = bound_stat[n_f_id + nfabor * bi->ihdepm];
        /* Calculation of the hydrodynamic forces and torques
         * applied on the deposited cluster :
         * Principle: drag>0 for protruding clusters, 0 otherwise */

        if (   p_height > mean_depo_height
            && cs_lagr_particles_get_flag(p_set, p_id, CS_LAGR_PART_DEPOSITED)) {

          /* Calculation of drag force and torque*/
          drag_force[0] =   3.0 * cs_math_pi * (p_height - mean_depo_height)
                          * (vvue[0] - vpart[0]) * visccf * romf * 3.39;
          drag_torque[0] = 0.0;

          for (cs_lnum_t id = 1; id < 3; id++) {

            drag_force[id]  =   3.0 * cs_math_pi * (p_height - mean_depo_height)
                              * (vvue[id] - vpart[id]) * visccf * romf * 1.7;
            drag_torque[id] = 1.4 * drag_force[id] * p_diam * 0.5;

          }

          /* Calculation of lift force and torque */
          lift_force[0] = - 20.0 * cs_math_pow2(visccf) * romf *
            pow(ustar * (p_height - mean_depo_height) * 0.5 / visccf, 2.31);

          for (cs_lnum_t id = 1; id < 3; id++) {

            lift_torque[id] = -lift_force[0] * p_diam * 0.5
              * vvue[id] / sqrt(cs_math_pow2(vvue[1]) + cs_math_pow2(vvue[2]) );

          }

          /* Calculation of gravity force and torque */
          for (cs_lnum_t id = 0; id < 3; id++) {
            grav_force[id] = p_diam * ggp[id]
              *(p_height - mean_depo_height) / p_height;
          }

          for (cs_lnum_t id = 1; id < 3; id++) {

            cs_real_t sign11, sign22;
            if (grav_force[1]*vvue[1] < 0 )
              sign11 = -1.0;
            else
              sign11 = 1.0;
            if (grav_force[2]*vvue[2] < 0 )
              sign22 = -1.0;
            else
              sign22 = 1.0;

            grav_torque[id] = (-grav_force[0] +
                                grav_force[1]*sign11 +
                                grav_force[2]*sign22 ) * p_diam * 0.5
              * vvue[id] / sqrt(cs_math_pow2(vvue[1]) + cs_math_pow2(vvue[2]) );

          }

        }
        else if (   p_height <= mean_depo_height
                 && cs_lagr_particles_get_flag(p_set, p_id,
                                               CS_LAGR_PART_DEPOSITED)) {

          /* Calculation of drag force and torque*/
          for (cs_lnum_t id = 0; id < 3; id++) {
            drag_force[id]  = 0.0;
            drag_torque[id] = 0.0;
          }

          /* Calculation of lift force and torque */
          lift_force[0] = 0.0;

          for (cs_lnum_t id = 1; id < 3; id++) {
            lift_torque[id] = 0.0;
          }

          /* Calculation of gravity force and torque */
          for (cs_lnum_t id = 0; id < 3; id++) {
            grav_force[id] = p_diam * ggp[id];
          }

          for (cs_lnum_t id = 1; id < 3; id++) {

            cs_real_t sign11, sign22;
            if (grav_force[1]*vvue[1] < 0 )
              sign11 = -1.0;
            else
              sign11 = 1.0;
            if (grav_force[2]*vvue[2] < 0 )
              sign22 = -1.0;
            else
              sign22 = 1.0;

            grav_torque[id] = (-grav_force[0] +
                                grav_force[1]*sign11 +
                                grav_force[2]*sign22 ) * p_diam * 0.5
              * vvue[id] / sqrt(cs_math_pow2(vvue[1]) + cs_math_pow2(vvue[2]) );

          }

        }
        else if (cs_lagr_particles_get_flag(p_set, p_id, CS_LAGR_PART_ROLLING)) {

          /* Calculation of drag force and torque*/
          drag_force[0] =  3.0 * cs_math_pi * p_diam
            * (vvue[0] - vpart[0]) * visccf * romf * 3.39;
          drag_torque[0] = 0.0;

          for (cs_lnum_t id = 1; id < 3; id++) {

            drag_force[id]   =  3.0 * cs_math_pi * p_diam
              * (vvue[id] - vpart[id]) * visccf * romf * 1.7;
            drag_torque[id] = 1.4 * drag_force[id] * p_diam * 0.5;

          }

          /* Calculation of lift force and torque */
          lift_force[0] = - 20.0 * cs_math_pow2(visccf) * romf *
            pow(ustar * p_diam * 0.5 / visccf, 2.31);

          for (cs_lnum_t id = 1; id < 3; id++) {

            lift_torque[id] = -lift_force[0] * p_diam * 0.5
              * vvue[id] / sqrt(cs_math_pow2(vvue[1]) + cs_math_pow2(vvue[2]) );

          }

          /* Calculation of gravity force and torque */
          for (cs_lnum_t id = 0; id < 3; id++) {
            grav_force[id] = p_diam * ggp[id];
          }

          for (cs_lnum_t id = 1; id < 3; id++) {

            cs_real_t sign11, sign22;
            if (grav_force[1]*vvue[1] < 0 )
              sign11 = -1.0;
            else
              sign11 = 1.0;
            if (grav_force[2]*vvue[2] < 0 )
              sign22 = -1.0;
            else
              sign22 = 1.0;

            grav_torque[id] = (-grav_force[0] +
                                grav_force[1]*sign11 +
                                grav_force[2]*sign22 ) * p_diam * 0.5
              * vvue[id] / sqrt(cs_math_pow2(vvue[1]) + cs_math_pow2(vvue[2]) );

          }
        }

        /* Adhesion force and torque for multilayered structures:
         * equal to adhesion force between single particles
         * times the number of particle-particle contacts  */

        cs_real_t adhes_energ, adhes_force;
        cs_lagr_adh_pp(p_diam, tempf, &adhes_energ, &adhes_force);
        /* Average number of contact in a cluster */
        cs_real_t ncont_pp = cs_math_pow2(p_diam/diam_mean);

        /* Case with consolidation */
        if (   cs_glob_lagr_consolidation_model->iconsol > 0
            && cs_lagr_particles_get_flag(p_set, p_id,
                                          CS_LAGR_PART_DEPOSITED)) {

          const cs_real_t consol_height
            = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_CONSOL_HEIGHT);

          if (consol_height > 0.01 * diam_mean) {
            adhes_force = cs_glob_lagr_consolidation_model->force_consol +
              (adhes_force - cs_glob_lagr_consolidation_model->force_consol) * 0.5
            * (1.0 + tanh((mean_depo_height - consol_height)
                          /(0.1 * consol_height)));
          }
          else {
            adhes_force *= ncont_pp;
          }
          cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_N_SMALL_ASPERITIES,
                                    ncont_pp);
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_ADHESION_FORCE,
                                    adhes_force);
          cs_real_t adhes_tor = adhes_force * p_diam * 0.5;
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_ADHESION_TORQUE,
                                    adhes_tor);
        }
        else {
          /* Case without consolidation */

          int ncont = 1;

          if (ncont_pp > 600.0) {

            cs_real_t rtmp;

            cs_random_normal(1, &rtmp);

            ncont = (int)ncont_pp + sqrt(ncont_pp) * rtmp;

          }
          else {

            cs_random_poisson(1, ncont_pp, &ncont);

          }
          ncont = cs::max(1, ncont);
          cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_N_SMALL_ASPERITIES,
                                    ncont);

          adhes_energ *= ncont;
          adhes_force *= ncont;
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_ADHESION_FORCE,
                                    adhes_force);

          cs_real_t adhes_tor = adhes_force * p_diam * 0.5;
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_ADHESION_TORQUE,
                                    adhes_tor);

        }

        for (cs_lnum_t id = 1; id < 3; id++) {
          adhes_torque[id] =
            - cs_lagr_particle_get_real(particle, p_am, CS_LAGR_ADHESION_TORQUE)
            * vvue[id] / sqrt(cs_math_pow2(vvue[1]) + cs_math_pow2(vvue[2]) );
        }


        /* Is there direct wall-normal lift-off of the cluster ?  */
        if ((adhes_force + grav_force[0] + lift_force[0] + drag_force[0]) < 0
             && iresusp == 0 ) {

          /* Update of the number and weight of resuspended particles     */
          p_set->n_part_resusp += 1;
          p_set->weight_resusp += p_stat_w;

          if (events != nullptr) {
            const cs_real_t *part_vel
              = cs_lagr_particles_attr_get_const_ptr<cs_real_t>(p_set, p_id,
                                                                CS_LAGR_VELOCITY);
            _add_resuspension_event(events,
                                    p_set,
                                    p_id,
                                    face_id,
                                    part_vel);
          }

          /* Update of surface covered and deposit height
           * (for non-rolling particles) */
          if (cs_lagr_particles_get_flag(p_set, p_id, CS_LAGR_PART_DEPOSITED)) {

            bound_stat[n_f_id + nfabor * bi->iscovc] -=
              cs_math_pi * cs_math_pow2(p_diam) * p_stat_w
              * 0.25 / mq->b_face_surf[n_f_id];

            bound_stat[n_f_id + nfabor * bi->ihdepm] -=
              cs_math_pi * p_height * cs_math_pow2(p_diam) * p_stat_w
              * 0.25 / mq->b_face_surf[n_f_id];

            bound_stat[n_f_id + nfabor * bi->ihdepv] -=
              cs_math_pow2(cs_math_pi * p_height * cs_math_pow2(p_diam) * p_stat_w
                   * 0.25 / mq->b_face_surf[n_f_id]);

            bound_stat[n_f_id + nfabor * bi->inclg] -=
              p_stat_w;

          }

          /* The particle is resuspended    */
          vpart[0] = cs::min(-1.0 / p_mass * dt_part
                            * cs::abs(-lift_force[0] - drag_force[0]
                                      - adhes_force -grav_force[0] ),
                            0.001);
          if (cs_lagr_particles_get_flag(p_set, p_id, CS_LAGR_PART_DEPOSITED)) {
            vpart[1] = 0.0;
            vpart[2] = 0.0;
          }

          cs_lagr_particles_unset_flag(p_set, p_id, CS_LAGR_PART_DEPOSITION_FLAGS);
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_ADHESION_FORCE, 0.0);
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_ADHESION_TORQUE, 0.0);

          if (cs::abs(p_height-p_diam)/p_diam > 1.0e-6) {
            cs_real_t d_resusp = pow(0.75 * cs_math_pow2(p_diam) * p_height, 1.0/3.0);
            cs_lagr_particle_set_real(particle, p_am, CS_LAGR_DIAMETER, d_resusp);
            cs_lagr_particle_set_real(particle, p_am, CS_LAGR_HEIGHT, d_resusp);
          }

          if (p_am->count[0][CS_LAGR_N_LARGE_ASPERITIES] > 0)
            cs_lagr_particle_set_lnum(particle, p_am,
                                      CS_LAGR_N_LARGE_ASPERITIES, 0);

          if (p_am->count[0][CS_LAGR_N_SMALL_ASPERITIES] > 0)
            cs_lagr_particle_set_lnum(particle, p_am,
                                      CS_LAGR_N_SMALL_ASPERITIES, 0);

          if (p_am->count[0][CS_LAGR_DISPLACEMENT_NORM] > 0)
            cs_lagr_particle_set_real(particle, p_am,
                                      CS_LAGR_DISPLACEMENT_NORM, 0.0);

          iresusp = 1;
        }

        /* No direct normal lift-off */
        else if (iresusp == 0) {

          /* Calculation of the norm of the hydrodynamic
           * torque and drag (tangential) */

          cs_real_t drag_tor_norm  = sqrt(  cs_math_pow2(drag_torque[1])
                                          + cs_math_pow2(drag_torque[2]));
          cs_real_t lift_tor_norm  = sqrt(  cs_math_pow2(lift_torque[1])
                                          + cs_math_pow2(lift_torque[2]));
          cs_real_t grav_tor_norm  = sqrt(  cs_math_pow2(grav_torque[1])
                                          + cs_math_pow2(grav_torque[2]));

          /* Differentiation between two cases:
           *  a) Already rolling clusters
           *  b) Deposited clusters being broken */

          if (cs_lagr_particles_get_flag(p_set, p_id, CS_LAGR_PART_ROLLING)) {

            cs_real_t iner_tor = (7.0 / 5.0) * p_mass
                               * cs_math_pow2((p_diam * 0.5));
            cs_real_t cst_1, cst_4;
            cst_4 =   6 * cs_math_pi * visccf
              * romf * 1.7 * 1.4
              * cs_math_pow2(p_diam * 0.5);
            cst_1 = cst_4 * (p_diam * 0.5) / iner_tor;

            for (cs_lnum_t id = 1; id < 3; id++) {

              vpart0[id] = vpart[id];
              vpart[id]  =  vpart0[id] * exp(-cst_1 * dt_part)
                + (   vvue[id] + (adhes_torque[id] + lift_torque[id]
                    + grav_torque[id]) / cst_4)
                * (1.0 - exp(-cst_1 * dt_part) );

            }

            cs_real_t scalax  = vpart[1] * vvue[1] + vpart[2] * vvue[2];

            if (scalax > 0.0) {

              /* The cluster continues to roll (depo_flag = 2)  */

              vpart[0] = 0.0;

              for (cs_lnum_t id = 1; id < 3; id++) {

                if (cs::abs (vpart[id]) > cs::abs (vvue[id]))
                  /* The velocity of the rolling particle cannot   */
                  /* exceed the surrounding fluid velocity    */
                  vpart[id] = vvue[id];

                cs_real_t kkk = vvue[id] + adhes_torque[id] / cst_4;
                cs_real_t kk  = vpart0[id] - kkk;

                depl[id] = kkk * dt_part + kk / cst_1
                                              * (1. - exp(-cst_1 * dt_part));

              }

              iresusp = 1;

            }
            /* if (scalax..) */
            else {

              /* The cluster stops moving
               * the flag is set to 1 and velocity and displacement are null */

              /* Simplified treatment:
               * no update of iscovc, ihdepm for particles that stop */
              iresusp = 1;

              for (cs_lnum_t id = 1; id < 3; id++) {

                depl[id]     = 0.0;
                vpart[id]    = 0.0;

              }

            } /* if (scalax..)   */

          }

          else if (   cs_lagr_particles_get_flag(p_set, p_id,
                                                 CS_LAGR_PART_DEPOSITED)
                   && iresusp == 0 ) {
            /* Cluster being broken:
             * we first check if there is resuspension */
            cs_real_t cond_resusp[2];
            for (cs_lnum_t id = 0; id < 2; id++) {
              cond_resusp[id] = (adhes_torque[id] + grav_torque[id]
                              + lift_torque[id] + drag_torque[id] ) * vvue[id];
            }

            if (cond_resusp[0] > 0.0 || cond_resusp[1] > 0.0) {
              iresusp = 1;
              cs_real_t clust_consol_height;
              cs_real_t height_reent;
              cs_real_t random;
              /* Sample of a possible break line */
              if (  cs_lagr_particle_get_real(particle, p_am, CS_LAGR_CONSOL_HEIGHT)
                  < diam_mean) {
                cs_random_uniform(1, &random);
                clust_consol_height = 0.0;
              }
              else {
                cs_real_t param = 1.0 - ((drag_tor_norm + lift_tor_norm
                                          - grav_tor_norm) /p_diam
                                          - 2.0 * adhes_force ) /
                  (cs_glob_lagr_consolidation_model->force_consol - adhes_force);
                if (param >= 1.)
                  clust_consol_height = p_height; // Very high consolidation
                else if (param <= -1.)
                  clust_consol_height = 0.; // Very high hydrodynamic forces
                else
                  clust_consol_height =
                    cs_lagr_particle_get_real(particle, p_am, CS_LAGR_CONSOL_HEIGHT)
                    * (1 + cs_glob_lagr_consolidation_model->slope_consol
                       * 0.5 * log((1.0+param)/(1.0-param) ) );
              }
              cs_random_uniform(1, &random);
              height_reent = random * (p_height - clust_consol_height);
              height_reent = cs::min(height_reent, p_height);
              height_reent = cs::max(height_reent, 0.0);

              /* Treatment of the new rolling particle */
              cs_lnum_t itreated = 0;
              cs_lnum_t nb_part_reent = height_reent / p_height *
              cs_lagr_particle_get_real(particle, p_am, CS_LAGR_CLUSTER_NB_PART);

              if (nb_part_reent < 1.0 && itreated == 0) {
                /* No resuspension (cluster too small)*/
                itreated = 1;

                for (cs_lnum_t id = 0; id < 3; id++) {
                  vpart[id] = 0.0;
                  depl[id]  = 0.0;
                }
              }
              else if ((cs_lagr_particle_get_real(particle, p_am, CS_LAGR_CLUSTER_NB_PART)
                         -nb_part_reent) < 1.0 && itreated == 0) {
                /* The whole cluster starts rolling*/
                cs_lagr_particles_unset_flag(p_set, p_id, CS_LAGR_PART_DEPOSITION_FLAGS);
                cs_lagr_particles_set_flag(p_set, p_id, CS_LAGR_PART_ROLLING);

                /* Update of surface covered and deposit height */
                bound_stat[n_f_id + nfabor * bi->iscovc] -=
                  cs_math_pi * cs_math_pow2(p_diam) * p_stat_w
                  * 0.25 / mq->b_face_surf[n_f_id];

                bound_stat[n_f_id + nfabor * bi->ihdepm] -=
                  cs_math_pi * p_height * cs_math_pow2(p_diam) * p_stat_w
                  * 0.25 / mq->b_face_surf[n_f_id];

                bound_stat[n_f_id + nfabor * bi->ihdepv] -=
                  cs_math_pow2(cs_math_pi * p_height * cs_math_pow2(p_diam) * p_stat_w
                       * 0.25 / mq->b_face_surf[n_f_id]);

                bound_stat[n_f_id + nfabor * bi->inclg] -=
                  p_stat_w;

                itreated = 1;
                cs_real_t d_resusp = pow(0.75 * cs_math_pow2(p_diam) * p_height, 1.0/3.0);
                cs_lagr_particle_set_real(particle, p_am, CS_LAGR_DIAMETER, d_resusp);
                cs_lagr_particle_set_real(particle, p_am, CS_LAGR_HEIGHT, d_resusp);

                /* Treatment of cluster motion */
                cs_real_t iner_tor = (7.0 / 5.0) * p_mass * cs_math_pow2((p_diam * 0.5));
                cs_real_t cst_1, cst_4;
                cst_4 =   6 * cs_math_pi * visccf
                  * romf * 1.7 * 1.4
                  * cs_math_pow2(p_diam * 0.5);
                cst_1 = cst_4 * (p_diam * 0.5) / iner_tor;

                for (cs_lnum_t id = 1; id < 3; id++) {
                  vpart0[id] = vpart[id];
                  vpart[id]  =  vpart0[id] * exp(-cst_1 * dt_part)
                    + (vvue[id] + (  adhes_torque[id] + lift_torque[id]
                                   + grav_torque[id]) / cst_4)
                    * (1.0 - exp(-cst_1 * dt_part) );
                }

                vpart[0] = 0.0;
                for (cs_lnum_t id = 1; id < 3; id++) {
                  if (cs::abs (vpart[id]) > cs::abs (vvue[id]))
                    /* The velocity of the rolling particle cannot   */
                    /* exceed the surrounding fluid velocity    */
                    vpart[id] = vvue[id];
                  cs_real_t kkk = vvue[id] + adhes_torque[id] / cst_4;
                  cs_real_t kk  = vpart0[id] - kkk;
                  depl[id] = kkk * dt_part + kk / cst_1
                                                * (1. - exp(-cst_1 * dt_part));
                }

              }
              else if (itreated == 0) {
                /* Duplication of the particle */
                *nresnew = *nresnew + 1;
                cs_lagr_particle_set_resize(p_set->n_particles + *nresnew);
                cs_lagr_part_copy(p_set->n_particles + *nresnew, p_id);

                /* We split both particles:
                * Part p_id stays while the new one starts rolling */
                unsigned char *new_part = p_set->p_buffer
                  + p_am->extents * p_set->n_particles+*nresnew;
                cs_real_t nb_resusp =  height_reent / p_height
                  * cs_lagr_particle_get_real(particle, p_am, CS_LAGR_CLUSTER_NB_PART);
                cs_lagr_particle_set_real(new_part, p_am, CS_LAGR_CLUSTER_NB_PART,
                                          nb_resusp);
                cs_real_t m_resusp = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS)
                  * cs_lagr_particle_get_real(new_part, p_am, CS_LAGR_CLUSTER_NB_PART)
                  / cs_lagr_particle_get_real(particle, p_am, CS_LAGR_CLUSTER_NB_PART);
                cs_lagr_particle_set_real(new_part, p_am, CS_LAGR_MASS, m_resusp);
                cs_real_t d_resusp = pow(0.75 * cs_math_pow2(p_diam) * p_height, 1.0/3.0);
                cs_lagr_particle_set_real(new_part, p_am, CS_LAGR_DIAMETER, d_resusp);
                cs_lagr_particle_set_real(new_part, p_am, CS_LAGR_HEIGHT, d_resusp);

                cs_lagr_particles_unset_flag(p_set, p_id, CS_LAGR_PART_DEPOSITION_FLAGS);
                cs_lagr_particles_set_flag(p_set, p_id, CS_LAGR_PART_DEPOSITED);
                vpart[0] = 0.0;
                for (cs_lnum_t id = 1; id < 3; id++) {
                  vpart[id] = 0.0;
                  depl[id]  = 0.0;
                }

                /* Update of deposit height */
                cs_real_t d_stay = cs_lagr_particle_get_real(particle, p_am,
                                                              CS_LAGR_HEIGHT);
                cs_lagr_particle_set_real(particle, p_am, CS_LAGR_HEIGHT, d_stay);

                bound_stat[n_f_id + nfabor * bi->ihdepm]
                  -= cs_math_pi * height_reent * cs_math_pow2(p_diam) * p_stat_w
                     * 0.25 / mq->b_face_surf[n_f_id];

                bound_stat[n_f_id + nfabor * bi->ihdepv]
                  -= cs_math_pow2(cs_math_pi * height_reent * cs_math_pow2(p_diam) * p_stat_w
                    * 0.25 / mq->b_face_surf[n_f_id]);

                cs_real_t nb_stay
                  =   cs_lagr_particle_get_real(particle, p_am,
                                                 CS_LAGR_CLUSTER_NB_PART)
                    - cs_lagr_particle_get_real(new_part, p_am,
                                                CS_LAGR_CLUSTER_NB_PART);
                cs_lagr_particle_set_real(particle, p_am, CS_LAGR_CLUSTER_NB_PART, nb_stay);

                cs_real_t mp_stay
                  = p_mass - cs_lagr_particle_get_real(new_part, p_am,
                                                       CS_LAGR_MASS);
                cs_lagr_particle_set_real(particle, p_am, CS_LAGR_MASS, mp_stay);

                /* The new particle starts rolling */
                cs_lagr_particles_unset_flag(p_set, p_id,
                                             CS_LAGR_PART_DEPOSITION_FLAGS);
                cs_lagr_particles_set_flag(p_set, p_id,
                                           CS_LAGR_PART_ROLLING);
                cs_lagr_particle_set_real(new_part, p_am,
                                          CS_LAGR_ADHESION_FORCE, adhes_force);
                cs_lagr_particle_set_lnum(new_part, p_am,
                                          CS_LAGR_N_SMALL_ASPERITIES,
                                          cs::max(1,ncont_pp));
                cs_lagr_particle_set_real(new_part, p_am,
                                          CS_LAGR_ADHESION_TORQUE,
                                          adhes_force * d_resusp * 0.5);
                cs_lagr_particle_set_real(new_part, p_am, CS_LAGR_DEPO_TIME, 0.0);
                cs_lagr_particle_set_real(new_part, p_am, CS_LAGR_CONSOL_HEIGHT, 0.0);

                itreated = 1;
              }

            } /* End of condition for resuspension */

            else {
              /* No rolling occurs */
              iresusp = 1;
              for (cs_lnum_t id = 1; id < 3; id++) {
                vpart[id] = 0.0;
                depl[id]  = 0.0;
              }
            }

          }

        }

      } /* End of multilayer resuspension */

    }

  }
  /* if ireent.eq.0 --> Motionless deposited particle */
  else {

    if (cs_lagr_particles_get_flag(p_set, p_id,
                                   CS_LAGR_PART_DEPOSITION_FLAGS)) {
      for (cs_lnum_t id = 1; id < 3; id++) {
        vpart[id] = 0.0;
        vvue[id]  = 0.0;
        depl[id]  = 0.0;
      }
    }

  }

  /* Reference frame change:
   * -----------------------
   * local reference frame for the boundary face --> global reference frame
   * NB: Inverse transformation: transpose of rot_m
   * ============================================== */

  /* 3.1 - Displacement   */

  cs_real_t depg[3];

  cs_math_33t_3_product(rot_m, depl, depg);

  /* 3.2 - Particle velocity   */

  cs_real_t *part_vel =
    cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                             CS_LAGR_VELOCITY);

  cs_math_33t_3_product(rot_m, vpart, part_vel);

  /* 3.3 - flow-seen velocity  */

  cs_real_t *part_vel_seen =
    cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                             CS_LAGR_VELOCITY_SEEN);

  cs_math_33t_3_product(rot_m, vvue, part_vel_seen);

  /* Computation of the new particle position
   * ======================================== */

  cs_real_t *part_coords =
    cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                             CS_LAGR_COORDS);
  for (cs_lnum_t id = 0; id < 3; id++)
    part_coords[id] += depg[id];
}

/*----------------------------------------------------------------------------*/
/*! \brief Deposition submodel
 *
 *  Main subroutine of the submodel
 *   1/ Calculation of the normalized wall-normal distance of
 *           the boundary-cell particles
 *   2/ Sorting of the particles with respect to their normalized
 *           wall-normal distance
 *         * If y^+ > depint : the standard Langevin model is applied
 *         * If y^+ < depint : the deposition submodel is applied
 *
 * \param[in]  p_id      particle index in set
 * \param[in]  dt_part   remaining time step associated to the particle
 * \param[in]  nor       current step id (for 2nd order scheme)
 * \param[in]  taup      dynamic characteristic time
 * \param[in]  tlag      lagrangian fluid characteristic time
 * \param[in]  piil      term in integration of up sdes
 * \param[in]  bx        turbulence characteristics
 * \param[in]  vagaus    gaussian random variables
 * \param[in]  romp      particles associated density
 * \param[in]  force_p   forces per mass unit on particles (m/s^2)
 * \param[in]  vislen    nu/u* = y/y+
 * \param[out] nresnew
 *
 */
/*----------------------------------------------------------------------------*/

static void
cs_sde_vels_pos_time_integ_depot(cs_lnum_t                       p_id,
                                 cs_real_t                       dt_part,
                                 int                             nor,
                                 const cs_real_t                 taup,
                                 const cs_real_3_t               tlag,
                                 const cs_real_3_t               piil,
                                 const cs_real_33_t              bx,
                                 const cs_real_3_t              *vagaus,
                                 const cs_real_t                 romp,
                                 const cs_real_3_t               force_p,
                                 const cs_real_t                 vislen[],
                                 cs_lnum_t                      *nresnew)
{
  /* Particles management */
  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;

  cs_lagr_extra_module_t *extra = cs_get_lagr_extra_module();

  /* Initialisations*/

  cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

  cs_real_t vitf = 0.0;

  cs_real_t aux1, aux2, aux3, aux4, aux5, aux6, aux7, aux8, aux9, aux10, aux11;
  cs_real_t ter1f, ter2f, ter3f;
  cs_real_t ter1p, ter2p, ter3p, ter4p, ter5p;
  cs_real_t ter1x, ter2x, ter3x, ter4x, ter5x;
  cs_real_t p11, p21, p22, p31, p32, p33;
  cs_real_t omega2, gama2, omegam;
  cs_real_t grga2, gagam, gaome;

  const cs_temperature_scale_t t_scl = cs_glob_thermal_model->temperature_scale;

  const int _prev_id = (extra->vel->n_time_vals > 1) ? 1 : 0;
  const cs_real_3_t *cvar_vel
    = (const cs_real_3_t *)(extra->vel->vals[_prev_id]);

  /* Interface location between near-wall region */
  /* and core of the flow (normalized units) */

  cs_real_t depint      = 100.0;

  /* Tracking events if requested */

  cs_lagr_event_set_t  *events = nullptr;

  if (cs_lagr_stat_is_active(CS_LAGR_STAT_GROUP_TRACKING_EVENT))
    events = cs_lagr_event_set_boundary_interaction();

  /* Loop on the particles
   * Note: new particles will be integrated at the next time step, otherwise
   * positions might be overwritten */

    unsigned char *particle = p_set->p_buffer + p_am->extents * p_id;

  int imposed_motion = cs_lagr_particles_get_flag(p_set, p_id,
                                                  CS_LAGR_PART_IMPOSED_MOTION);

  if (! imposed_motion) {

    /* use previous step for t_order == 1 or prediction step
     * and current one for correction step */
    cs_lnum_t cell_id = cs_lagr_particle_get_lnum_n(particle, p_set->p_am, 2-nor,
                                                    CS_LAGR_CELL_ID);

    cs_real_t *old_part_vel      =
      cs_lagr_particle_attr_n_get_ptr<cs_real_t>(particle, p_am, 1,
                                                 CS_LAGR_VELOCITY);
    cs_real_t *old_part_vel_seen =
      cs_lagr_particle_attr_n_get_ptr<cs_real_t>(particle, p_am, 1,
                                                 CS_LAGR_VELOCITY_SEEN);
    cs_real_t *part_vel          =
      cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                               CS_LAGR_VELOCITY);
    cs_real_t *part_vel_seen     =
      cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                               CS_LAGR_VELOCITY_SEEN);
    cs_real_t *part_coords       =
      cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                               CS_LAGR_COORDS);
    cs_real_t *old_part_coords   =
      cs_lagr_particle_attr_n_get_ptr<cs_real_t>(particle, p_am, 1,
                                                 CS_LAGR_COORDS);

    /* Fluid temperature computation depending on the type of flow  */
    cs_real_t tempf;

    if (   extra->temperature != nullptr
        && t_scl == CS_TEMPERATURE_SCALE_CELSIUS)
      tempf = extra->temperature->val[cell_id] + tkelvi;

    else if (   extra->temperature != nullptr
             && t_scl == CS_TEMPERATURE_SCALE_KELVIN)
      tempf = extra->temperature->val[cell_id];
    else {
      tempf = cs_glob_fluid_properties->t0;
      if (t_scl == CS_TEMPERATURE_SCALE_CELSIUS)
        tempf += tkelvi;
    }

    /* If y^+ is greater than the interface location,
       the standard model is applied
       ============================================== */

    cs_lnum_t face_id = cs_lagr_particle_get_lnum(particle, p_am,
                                                  CS_LAGR_NEIGHBOR_FACE_ID);
    cs_real_t yplus = cs_lagr_particle_get_real(particle, p_am,
                                                CS_LAGR_YPLUS);

    int deposition_flags
      = cs_lagr_particles_get_flag(p_set, p_id,
                                   CS_LAGR_PART_DEPOSITION_FLAGS);

    if (face_id < 0 || (yplus > depint && deposition_flags == 0)) {

      cs_lagr_particle_set_lnum(particle,
                                p_am,
                                CS_LAGR_MARKO_VALUE,
                                CS_LAGR_COHERENCE_STRUCT_BULK);

      for (cs_lnum_t id = 0; id < 3; id++) {

        vitf = cvar_vel[cell_id][id];

        /* Preliminary computations
           ------------------------
           compute II*TL+<u> and [(grad<P>/rhop+g)*tau_p+<Uf>] ?  */

        cs_real_t tci = piil[id] * tlag[id] + vitf;
        cs_real_t v_lim = force_p[id] * taup;

        /* Compute deterministic coefficients/terms
           ---------------------------------------- */

        aux1 = exp(-dt_part / taup);
        aux2 = exp(-dt_part / tlag[id]);
        aux3 = tlag[id] / (tlag[id] - taup);
        aux4 = tlag[id] / (tlag[id] + taup);
        aux5 = tlag[id] * (1.0 - aux2);
        aux6 = cs_math_pow2(bx[id][nor-1]) * tlag[id];
        aux7 = tlag[id] - taup;
        aux8 = cs_math_pow2(bx[id][nor-1]) * cs_math_pow2(aux3);

        /* --> trajectory terms */
        cs_real_t aa = taup * (1.0 - aux1);
        cs_real_t bb = (aux5 - aa) * aux3;
        cs_real_t cc = dt_part - aa - bb;

        ter1x = aa * old_part_vel[id];
        ter2x = bb * old_part_vel_seen[id];
        ter3x = cc * tci;
        ter4x = (dt_part - aa) * v_lim;

        /* --> flow-seen velocity terms   */
        ter1f = old_part_vel_seen[id] * aux2;
        ter2f = tci * (1.0 - aux2);

        /* --> termes pour la vitesse des particules     */
        cs_real_t dd = aux3 * (aux2 - aux1);
        cs_real_t ee = 1.0 - aux1;

        ter1p = old_part_vel[id] * aux1;
        ter2p = old_part_vel_seen[id] * dd;
        ter3p = tci * (ee - dd);
        ter4p = v_lim * ee;

        /* Coefficients computation for the stochastic integral */
        /* Integral for particles position */
        gama2  = 0.5 * (1.0 - aux2 * aux2);
        omegam = aux3 * ( (tlag[id] - taup) * (1.0 - aux2)
                                - 0.5 * tlag[id] * (1.0 - aux2 * aux2)
                                + cs_math_pow2(taup) / (tlag[id] + taup)
                                * (1.0 - aux1 * aux2)
                                ) * aux6;
        omega2 =  aux7 * (aux7 * dt_part - 2.0 * (tlag[id] * aux5 - taup * aa))
                 + 0.5 * tlag[id] * tlag[id] * aux5 * (1.0 + aux2)
                 + 0.5 * taup * taup * aa * (1.0 + aux1)
                 - 2.0 * aux4 * tlag[id] * taup * taup
                       * (1.0 - aux1 * aux2);
        omega2 = aux8 * omega2;

        if (cs::abs(gama2) > cs_math_epzero) {

          p21 = omegam / sqrt(gama2);
          p22 = omega2 - cs_math_pow2(p21);
          p22 = sqrt(cs::max(0.0, p22));

        }
        else {

          p21 = 0.0;
          p22 = 0.0;

        }

        ter5x = p21 * vagaus[0][id] + p22 * vagaus[1][id];

        /* --> integral for the flow-seen velocity  */
        p11   = sqrt(gama2 * aux6);
        ter3f = p11 * vagaus[0][id];

        /* --> integral for the particles velocity  */
        aux9  = 0.5 * tlag[id] * (1.0 - aux2 * aux2);
        aux10 = 0.5 * taup * (1.0 - aux1 * aux1);
        aux11 =   taup * tlag[id]
                * (1.0 - aux1 * aux2)
                / (taup + tlag[id]);

        grga2 = (aux9 - 2.0 * aux11 + aux10) * aux8;
        gagam = (aux9 - aux11) * (aux8 / aux3);
        gaome = ( (tlag[id] - taup) * (aux5 - aa)
                  - tlag[id] * aux9
                  - taup * aux10
                  + (tlag[id] + taup) * aux11)
                * aux8;

        if (p11 > cs_math_epzero)
          p31 = gagam / p11;
        else
          p31 = 0.0;

        if (p22 > cs_math_epzero)
          p32 = (gaome - p31 * p21) / p22;
        else
          p32 = 0.0;

        p33 = grga2 - cs_math_pow2(p31) - cs_math_pow2(p32);
        p33 = sqrt(cs::max(0.0, p33));
        ter5p =   p31 * vagaus[0][id]
                + p32 * vagaus[1][id]
                + p33 * vagaus[2][id];

        /* Update of the particle state-vector */

        part_coords[id] =   old_part_coords[id]
                          + ter1x + ter2x + ter3x + ter4x + ter5x;

        part_vel_seen[id] =  ter1f + ter2f + ter3f;

        part_vel[id]      = ter1p + ter2p + ter3p + ter4p + ter5p;

      }

    }

    /* Otherwise, the deposition submodel is applied
     * ============================================= */

    else if (! (deposition_flags & CS_LAGR_PART_TO_DELETE)) {

      cs_lnum_t *marko_value
        = (cs_lnum_t *)cs_lagr_particle_attr(particle,
                                             p_am,
                                             CS_LAGR_MARKO_VALUE);

      if (yplus< cs_lagr_particle_get_real(particle, p_am,
                                            CS_LAGR_INTERF)) {

        if (*marko_value < 0)
          *marko_value = CS_LAGR_COHERENCE_STRUCT_DEGEN_INNER_ZONE_DIFF;
        else
          *marko_value = CS_LAGR_COHERENCE_STRUCT_INNER_ZONE_DIFF;

      }
      else {

        if (*marko_value < 0)
          *marko_value = CS_LAGR_COHERENCE_STRUCT_DEGEN_SWEEP;

        else if (*marko_value == CS_LAGR_COHERENCE_STRUCT_INNER_ZONE_DIFF)
          *marko_value = CS_LAGR_COHERENCE_STRUCT_DEGEN_EJECTION;

      }

      _lagesd(dt_part,
              p_id,
              nor,
              taup,
              piil,
              vagaus,
              romp,
              force_p,
              tempf,
              vislen,
              events,
              &depint,
              nresnew);

    }

  }

  /* Specific treatment for particles with imposed motion */

  else if (imposed_motion) {

      cs_real_t disp[3] = {0., 0., 0.};

      cs_real_t *old_part_coords =
        cs_lagr_particle_attr_n_get_ptr<cs_real_t>(particle, p_am, 1,
                                                   CS_LAGR_COORDS);
      cs_real_t *part_coords =
        cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                 CS_LAGR_COORDS);

      cs_real_t *part_vel_seen =
        cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                 CS_LAGR_VELOCITY_SEEN);

      cs_real_t *part_vel =
        cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                 CS_LAGR_VELOCITY);

    cs_user_lagr_imposed_motion(p_set,
                                p_id,
                                old_part_coords,
                                dt_part,
                                disp);

    for (cs_lnum_t id = 0; id < 3; id++) {

      part_coords[id] = old_part_coords[id] + disp[id];

      part_vel_seen[id] =  0.0;

      part_vel[id] = disp[id] / dt_part;

    }
  }

}


/*----------------------------------------------------------------------------*/
/*!
 * \brief Integration of particle equations of motion:
 *
 * - Standard Model : First or second order
 * - Deposition submodel (Guingo & Minier, 2007) if needed
 *
 * \param[in]  p_id      particle index in set
 * \param[in]  dt_part   remaining time step associated to the particle
 * \param[in]  nor       current step id (for 2nd order scheme)
 * \param[in]  taup      dynamic characteristic time
 * \param[in]  tlag      fluid characteristic time
 * \param[in]  piil      term in integration of U-P SDEs
 * \param[in]  bx        turbulence characteristics
 * \param[out] tsfext    info for return coupling source terms
 * \param[out] force_p   forces per mass unit on particles (m/s^2)
 * \param[in]  vislen    nu/u* = y/y+
 * \param[in]  beta      proportional to the gradient of T_lag
 * \param[out] vagaus    gaussian random variables
 * \param[out] brgaus    gaussian random variables
 * \param[in]  nresnew
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_sde(cs_lnum_t                        p_id,
            cs_real_t                        dt_part,
            int                              nor,
            const cs_real_t                 *taup,
            const cs_real_3_t               *tlag,
            const cs_real_3_t               *piil,
            const cs_real_33_t              *bx,
            cs_real_t                       *tsfext,
            const cs_real_3_t                force_p,
            const cs_real_t                  vislen[],
            const cs_real_3_t                beta,
            cs_real_3_t                     *vagaus,
            cs_real_6_t                      brgaus,
            cs_lnum_t                       *nresnew)
{
  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;

  cs_real_t romp;

  /* Computation of particle density */
  cs_real_t aa = 6.0 / cs_math_pi;
  cs_real_t d3 = cs_math_pow3(cs_lagr_particles_get_real(p_set, p_id,
                                                         CS_LAGR_DIAMETER));
  romp = aa * cs_lagr_particles_get_real(p_set, p_id, CS_LAGR_MASS) / d3;

  /* First order
     ----------- */

  if (cs_glob_lagr_time_scheme->t_order == 1) {

    /* If no deposition sub-model is activated, call of subroutine lages1
       for every particle */

    if (cs_glob_lagr_model->deposition <= 0)
      cs_sde_vels_pos_1_st_order_time_integ(p_id,
                                            dt_part,
                                            nor,
                                            taup,
                                            tlag,
                                            piil,
                                            bx,
                                            vagaus,
                                            brgaus,
                                            force_p,
                                            beta);

    /* Management of the deposition submodel */

    /* TODO extend to multiphase flow */
    else
      cs_sde_vels_pos_time_integ_depot(p_id,
                                       dt_part,
                                       nor,
                                       taup[0],
                                       tlag[0],
                                       piil[0],
                                       bx[0],
                                       vagaus,
                                       romp,
                                       force_p,
                                       vislen,
                                       nresnew);

  }

  /* Second order
     ------------ */

  else {

    /* TODO extend to multiphase flow */
    cs_sde_vels_pos_2_nd_order_time_integ(p_id,
                                          dt_part,
                                          nor,
                                          taup,
                                          tlag,
                                          piil,
                                          bx,
                                          tsfext,
                                          vagaus,
                                          brgaus,
                                          force_p,
                                          beta);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Integration of a stochastic differential equation (SDE) for
 *        a user particle variable (attribute).
 *
 * \f[
 *  \frac{dV}{dt} = \frac{V - PIP}{TCARAC}
 * \f]
 *
 * When there is interaction with a boundary face, the integration
 * degenerates to order 1 (even if the 2nd order scheme is active).
 *
 * \param[in]  attr    attribute/variable
 * \param[in]  p_id    particle id
 * \param[in]  nor     current step id (for 2nd order scheme)
 * \param[in]  dt_part remaining time step associated to the particle
 * \param[in]  tcarac  variable characteristic time
 * \param[in]  pip     right-hand side associated with SDE
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_sde_attr(cs_lagr_attribute_t   attr,
                 const cs_lnum_t       p_id,
                 int                   nor,
                 const cs_real_t       dt_part,
                 cs_real_t             tcarac,
                 cs_real_t             pip)
{
  /* Particles management */
  cs_lagr_particle_set_t         *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am  = p_set->p_am;
  unsigned char *particle = p_set->p_buffer + p_am->extents * p_id;

  int ltsvar = 0;

  if (p_set->p_am->source_term_displ != nullptr) {
    if (p_set->p_am->source_term_displ[attr] >= 0)
      ltsvar = 1;
  }

  assert(nor == 1 || nor == 2);

  if (nor == 1) {

    if (tcarac <= 0.0)
      bft_error
        (__FILE__, __LINE__, 0,
         _("The characteristic time for the stochastic differential equation\n"
           "of variable %d should be > 0.\n\n"
           "Here, for particle %ld, its value is %e11.4."),
         attr, (long)p_id, tcarac);

    cs_real_t aux1 = dt_part / tcarac;
    if (cs_glob_lagr_time_scheme->cell_wise_integ == 1)
      aux1 *= cs_lagr_particles_get_real(p_set, p_id,
                                         CS_LAGR_REMAINING_INTEG_TIME);
    cs_real_t aux2 = exp(-aux1);

    cs_real_t ter1 = cs_lagr_particle_get_real_n(particle, p_am, 1, attr)*aux2;
    cs_real_t ter2 = pip * (1.0 - aux2);

    /* Pour le cas NORDRE= 1 ou s'il y a rebond,     */
    /* le ETTP suivant est le resultat final    */
    cs_lagr_particle_set_real(particle, p_am, attr, ter1 + ter2);

    /* Pour le cas NORDRE= 2, on calcule en plus TSVAR pour NOR= 2  */
    if (ltsvar) {
      cs_real_t *part_ptsvar = cs_lagr_particles_source_terms(p_set, p_id,
                                                              attr);
      cs_real_t ter3 = (-aux2 + (1.0 - aux2) / aux1) * pip;
      *part_ptsvar = 0.5 * ter1 + ter3;

    }

  }
  else if (nor == 2) {

    if (cs_lagr_particles_get_lnum(p_set, p_id, CS_LAGR_REBOUND_ID) > 0)
    return;

    if (tcarac <= 0.0)
      bft_error
        (__FILE__, __LINE__, 0,
         _("The characteristic time for the stochastic differential equation\n"
           "of variable %d should be > 0.\n\n"
           "Here, for particle %ld, its value is %e11.4."),
         attr, (long)p_id, tcarac);

    cs_real_t aux1   = dt_part / tcarac;
    if (cs_glob_lagr_time_scheme->cell_wise_integ == 1)
      aux1 *= cs_lagr_particles_get_real(p_set, p_id,
                                         CS_LAGR_REMAINING_INTEG_TIME);
    cs_real_t aux2   = exp(-aux1);
    cs_real_t ter1   = 0.5 * cs_lagr_particle_get_real_n(particle, p_am, 1,
                                                         attr) * aux2;
    cs_real_t ter2   = pip  * (1.0 - (1.0 - aux2) / aux1);

    /* Pour le cas NORDRE= 2, le ETTP suivant est le resultat final */
    cs_real_t *part_ptsvar = cs_lagr_particles_source_terms(p_set, p_id, attr);
    cs_lagr_particle_set_real(particle, p_am, attr,
                              *part_ptsvar + ter1 + ter2);

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
