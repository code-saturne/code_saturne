/*============================================================================
 * Lagrangian boundary condition types and injection definitions.
 *============================================================================*/

/* VERS */

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local (user defined) function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Computation of particle injection profile.
 *
 * Note: if the input pointer is non-NULL, it must point to valid data
 * when the selection function is called, so that value or structure should
 * not be temporary (i.e. local);
 *
 * parameters:
 *   zone_id     <-- id of associated mesh zone
 *   location_id <-- id of associated mesh location
 *   input       <-- pointer to optional (untyped) value or structure.
 *   n_elts      <-- number of zone elements
 *   elt_ids     <-- ids of zone elements
 *   profile     <-- weight of a given zone element (size: n_elts)
 *----------------------------------------------------------------------------*/

/*! [lagr_bc_profile_func_2] */
static void
_injection_profile(int               zone_id,
                   int               location_id,
                   const void       *input,
                   cs_lnum_t         n_elts,
                   const cs_lnum_t   elt_ids[],
                   cs_real_t         profile[])
{
  const cs_real_3_t  *b_face_coords
    = (const cs_real_3_t *)cs_glob_mesh_quantities->b_face_cog;

  const int itmx = 8;

  /* Data initializations with experimental measurements
     --------------------------------------------------- */

  /* transverse coordinate */

  cs_real_t zi[] = {0.e-3, 1.e-3, 1.5e-3, 2.0e-3, 2.5e-3, 3.0e-3, 3.5e-3,
                    4.0e-3, 4.5e-3, 5.0e-3};

  /* particle volume fraction */

  cs_real_t lvf[] = {0.377e-4, 2.236e-4, 3.014e-4, 4.306e-4, 5.689e-4,
                     8.567e-4, 7.099e-4, 4.520e-4, 2.184e-4, 0.377e-4};

  /* vertical mean velocity of the particles */

  cs_real_t ui[] = {5.544, 8.827, 9.068, 9.169, 8.923, 8.295, 7.151, 6.048,
                    4.785, 5.544};

  /* Loop en elements
     ---------------- */

  for (cs_lnum_t ei = 0; ei < n_elts; ei++) {

    /* Face center */

    const cs_lnum_t face_id = elt_ids[ei];

    const cs_real_t z = b_face_coords[face_id][2];

    /* Interpolation */

    int i = 0;

    if (z > zi[0]) {
      for (i = 0; i < itmx; i++) {
        if (z >= zi[i] && z < zi[i+1])
          break;
      }
    }

    /* Compute volume fraction and statistical weight */

    cs_real_t up = ui[i] +(z-zi[i])*(ui[i+1]-ui[i])/(zi[i+1]-zi[i]);
    cs_real_t lvfp = lvf[i] + (z-zi[i])*(lvf[i+1]-lvf[i])/(zi[i+1]-zi[i]);

    /* number of particles in the cell */

    profile[ei] = lvfp * up;
  }
}
/*! [lagr_bc_profile_func_2] */

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define particle boundary conditions.
 *
 * This is used for the definition of inlet and other boundaries,
 * based on predefined boundary zones (\ref cs_zone_t).
 *
 * \param[in] bc_type    type of the boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_boundary_conditions(const int  bc_type[])
{
  cs_mesh_t *mesh = cs_glob_mesh;

  /*! [lagr_bc_variables] */
  cs_lagr_zone_data_t *lagr_bcs = cs_lagr_get_boundary_conditions();
  /*! [lagr_bc_variables] */

  /* Zone types
     ========== */

  /* For every boundary zone, we define the associated type:

       CS_LAGR_INLET     -> zone of particle inlet
       CS_LAGR_OUTLET    -> particle outlet
       CS_LAGR_REBOUND   -> rebound of the particles
       CS_LAGR_DEPO1     -> definitive deposition
       CS_LAGR_DEPO2     -> definitive deposition, but the particle remains
                            in memory (useful only if iensi2 = 1)
       CS_LAGR_DEPO_DLVO -> deposition of the particle with DLVO forces
       CS_LAGR_FOULING   -> fouling (coal only physical_model = 2)
       CS_LAGR_SYM       -> symmetry condition for the particles (zero flux)

 */

  /* define zone types */

  /*! [lagr_bc_define_type_1] */
  {
     const cs_zone_t  *z;
     int n_zones = cs_boundary_zone_n_zones();

     /* default: rebound for all types */
     for (int z_id = 0; z_id < n_zones; z_id++) {
       lagr_bcs->zone_type[z_id] = CS_LAGR_REBOUND;
     }

     /* inlet and outlet for specified zones */
     z = cs_boundary_zone_by_name("inlet");
     lagr_bcs->zone_type[z->id] = CS_LAGR_INLET;

     z = cs_boundary_zone_by_name("outlet");
     lagr_bcs->zone_type[z->id] = CS_LAGR_OUTLET;
  }
  /*! [lagr_bc_define_type_1] */

  /* Injection per particle set into the calculation domain
     ====================================================== */

  /* For every injection (usually inlet) zone,
     we provide the following information:

     *   n_inject: number of particles injected per set and per zone
     *   injection_frequency: injection frequency. If injection_frequency = 0,
                              then the injection occurs only at the first
                              absolute iteration.

     *   injection_profile_func: optional pointer to profile definition
     *   injection_profile_input: associated input, or NULL

     *   cluster: number of the group to which the particle belongs
                  (only if one wishes to calculate statistics per group)

     *   velocity_profile: type of condition on the velocity
              = -1 imposed flow velocity
              =  0 imposed velocity along the normal direction of the
                    boundary face, with norm equal to velocity[0] (m/s)
              =  1 imposed velocity: we prescribe velocity[0] (m/s)
                                                  velocity[1] (m/s)
                                                  velocity[2] (m/s)
              =  2 user-defined profile

     *   temperature_profile: type of temperature condition
              = 0 fluid temperature
              = 1 imposed temperature: we prescribe the temperature

     *   density

     *   coal_number: number of the coal of the particle
                      (only if physical_model = 2) */

  /* Access an injection set for selected zone
     (created automatically if not previously accessed) */

  /*! [lagr_bc_define_injection_1] */
  {
    const cs_zone_t  *z = cs_boundary_zone_by_name("inlet");
    int set_id = 0;
    cs_lagr_injection_set_t *zis
      = cs_lagr_get_injection_set(lagr_bcs, z->id, set_id);

    /* Now define parameters for this class and set */

    zis->n_inject = 100;
    zis->injection_frequency = 1;

    /* Assign other attributes (could be done through the GUI) */

    if (cs_glob_lagr_model->n_stat_classes > 0)
      zis->cluster = set_id + 1;

    zis->velocity_profile = 0;
    zis->velocity_magnitude = 1.1;

    zis->stat_weight = 1.0;
    zis->flow_rate = 0.0;

    /* Mean value and standard deviation of the diameter */
    zis->diameter = 5e-05;
    zis->diameter_variance = 0.0;

    /* Density */

    zis->density = 2500.0;

    zis->fouling_index = 100.0;

    /* Temperature and Cp */

    if (cs_glob_lagr_specific_physics->itpvar == 1) {
      zis->temperature_profile = 1;
      zis->temperature = 20.0;

      zis->cp = 1400.;
      zis->emissivity = 0.7;
    }

  }
  /*! [lagr_bc_define_injection_1] */

  /*! [lagr_bc_define_injection_2] */
  {
    const cs_zone_t  *z = cs_boundary_zone_by_name("inlet");
    int set_id = 1;
    cs_lagr_injection_set_t *zis
      = cs_lagr_get_injection_set(lagr_bcs, z->id, set_id);

    /* Assign injection profile function */

    zis->injection_profile_func = _injection_profile;
    zis->injection_profile_input = NULL; /* default */

  }
  /*! [lagr_bc_define_injection_2] */

  /* same procedure for additional injections at other zones or sets... */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handling of a particle interaction with a boundary of type
 *        \ref CS_LAGR_BC_USER.
 *
 * In this example, the particle is simply deposited and marked for
 * elimination.
 *
 * \param[in, out]  particles       pointer to particle set
 * \param[in]       p_id            particle id
 * \param[in]       face_id         boundary face id
 * \param[in]       face_norm       unit face (or face subdivision) normal
 * \param[in]       c_intersect     coordinates of intersection with the face
 * \param[in]       t_intersect     relative distance (in [0, 1]) of the
 *                                  intersection point with the face relative
 *                                  to the initial trajectory segment
 * \param[in]       b_zone_id       boundary zone id of the matching face
 * \param[in, out]  event_flag      event flag in case events are available
 * \param[in, out]  tracking_state  particle tracking state
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_user_boundary_interaction(cs_lagr_particle_set_t    *particles,
                                  cs_lnum_t                  p_id,
                                  cs_lnum_t                  face_id,
                                  const cs_real_t            face_norm[3],
                                  const cs_real_t            c_intersect[3],
                                  cs_real_t                  t_intersect,
                                  int                        b_zone_id,
                                  int                       *event_flag,
                                  cs_lagr_tracking_state_t  *tracking_state)
{
  /* Update deposition-related counter */

  /*! [update] */

# pragma omp atomic
  particles->n_part_dep += 1;

# pragma omp atomic
  particles->weight_dep += cs_lagr_particle_get_real(particles,
                                                     p_id,
                                                     CS_LAGR_STAT_WEIGHT);

  /* Mark particle as deposited and update its coordinates */

  cs_lagr_particles_set_flag(particles, p_id, CS_LAGR_PART_DEPOSITED);

  cs_real_t  *particle_coord
    = cs_lagr_particles_attr(particles, p_id, CS_LAGR_COORDS);
  for (int k = 0; k < 3; k++)
    particle_coord[k] = c_intersect[k];

  /* Update event and particle state */

  *event_flag = *event_flag | (CS_EVENT_OUTFLOW | CS_EVENT_DEPOSITION);
  *tracking_state = CS_LAGR_PART_OUT;

  /*! [update] */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
