/*============================================================================
 * Methods for lagrangian module
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

#include "bft_printf.h"
#include "bft_error.h"
#include "bft_mem.h"

#include "fvm_periodicity.h"

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_interface.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_prototypes.h"
#include "cs_random.h"
#include "cs_search.h"
#include "cs_time_step.h"
#include "cs_timer_stats.h"
#include "cs_thermal_model.h"

#include "cs_field.h"
#include "cs_field_pointer.h"

#include "cs_lagr.h"
#include "cs_lagr_particle.h"
#include "cs_lagr_stat.h"
#include "cs_lagr_geom.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Local (user defined) function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User modification of newly injected particles.
 *
 * This function is called after the initialization of the new particles in
 * order to modify them according to new particle profiles (injection
 * profiles, position of the injection point, statistical weights,
 * correction of the diameter if the standard-deviation option is activated).
 *
 * This function is called for each injection zone and class. Particles
 * with ids between \c pset->n_particles and \c n_elts are initialized
 * but may be modidied by this function.
 *
 * \param[in,out]  particles         particle set
 * \param[in]      zis               zone injection set data
 * \param[in]      particle_range    start and past-the-end ids of new particles
 *                                   for this zone and class
 * \param[in]      particle_face_id  face ids of new particles if zone is
 *                                   a boundary,  NULL otherwise
 * \param[in]      visc_length       viscous layer thickness
 *                                   (size: number of mesh boundary faces)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_in(cs_lagr_particle_set_t         *particles,
                const cs_lagr_injection_set_t  *zis,
                const cs_lnum_t                 particle_range[2],
                const cs_lnum_t                 particle_face_id[],
                const cs_real_t                 visc_length[])
{
  /* Initialization
     -------------- */

  const int itmx = 8;

  /* transverse coordinate */
  cs_real_t  zi[] = {0.e-3 , 1.e-3 , 1.5e-3, 2.0e-3, 2.5e-3,
                     3.0e-3, 3.5e-3, 4.0e-3, 4.5e-3, 5.0e-3};

  /* particle volume fraction */

  cs_real_t lvf[] = {0.377e-4, 2.236e-4, 3.014e-4, 4.306e-4, 5.689e-4,
                     8.567e-4, 7.099e-4, 4.520e-4, 2.184e-4, 0.377e-4};

  /* vertical mean velocity of the particles */
  cs_real_t  ui[] = {5.544e0, 8.827e0, 9.068e0, 9.169e0, 8.923e0,
                     8.295e0, 7.151e0, 6.048e0, 4.785e0, 5.544e0};

  /* transverse mean velocity of the particles */
  cs_real_t  wi[] = {0.e0   , 0.179e0, 0.206e0, 0.221e0, 0.220e0,
                     0.223e0, 0.206e0, 0.190e0, 0.195e0, 0.504e0};

  /* fluctuation of the vertical velocity of the particles */
  cs_real_t  uf[] = {0.352e0, 0.352e0, 0.275e0, 0.252e0, 0.367e0,
                     0.516e0, 0.657e0, 0.872e0, 1.080e0, 0.792e0};

  /* fluctuation of the transverse velocity of the particles */
  cs_real_t  wf[] = {0.058e0, 0.058e0, 0.056e0, 0.056e0, 0.060e0,
                     0.063e0, 0.058e0, 0.072e0, 0.091e0, 0.232e0};

#if 0
  /* shear-stress (currently not used) of the particle velocity */
  cs_real_t  uvi[] = {0.0017e0, 0.0017e0,  0.0016e0,  0.0027e0,  0.0077e0,
                      0.0146e0, 0.0206e0,  0.0447e0,  0.0752e0,  0.1145e0};
#endif

  const int ntcabs = cs_glob_time_step->nt_cur;
  const cs_real_t dtp = cs_glob_lagr_time_step->dtp;

  cs_lagr_zone_data_t  *lagr_bcs = cs_lagr_get_boundary_conditions();

  /* Loop on new particles
     --------------------- */

  for (cs_lnum_t p_id = particle_range[0]; p_id < particle_range[1]; p_id++) {

    // cs_lnum_t face_id = particle_face_id[p_id - particle_range[0]];

    /* Data initializations with experimental measurements */

    const cs_real_t *part_coords = cs_lagr_particles_attr_const(particles,
                                                                p_id,
                                                                CS_LAGR_COORDS);
    cs_real_t z = part_coords[2];

    /* Interpolation */

    int it = 0;

    if (z > zi[0]) {
      for (it = 0; it < itmx; it++) {
        if (z >= zi[it] && z < zi[it+1])
          break;
      }
    }

    /* Calculation of particle velocity */

    cs_real_t up  =   ui[it] + (z - zi[it]) * (ui[it+1] - ui[it])
                    / (zi[it+1] - zi[it]);

    /* The value of the mean transverse velocity is currently set to zero
     * due to uncertainties on this variable */

    cs_real_t wp;
    if (false)
      wp  = wi[it] + (z - zi[it]) * (wi[it+1] - wi[it]) / (zi[it+1] - zi[it]);
    else
      wp = 0.0;

    cs_real_t upp = uf[it] +   (z - zi[it]) * (uf[it+1] - uf[it])
                             / (zi[it+1] - zi[it]);
    cs_real_t wpp = wf[it] +   (z - zi[it]) * (wf[it+1] - wf[it])
                             / (zi[it+1] - zi[it]);

    /* Calculations of the instantaneous particle velocity */

    cs_real_t vgauss[3];

    cs_random_normal(3, vgauss);

    cs_real_t *part_vel
      = cs_lagr_particles_attr(particles, p_id, CS_LAGR_VELOCITY);
    part_vel[0] = up  + vgauss[0] * upp;
    part_vel[1] = 0.0;
    part_vel[2] = wp + vgauss[1] * wpp;

    /* Diameter */

    cs_real_t dp = 49.3e-6 + vgauss[2] * 4.85e-6;
    if (dp < 30.e-6) dp = 30.e-6;
    else if (dp > 70.e-6) dp = 70.e-6;

    cs_lagr_particles_set_real(particles, p_id, CS_LAGR_DIAMETER, dp);

  }

  /* Trick to average the statistics at iteration nstist
   * starting from an unsteady two-coupling calculation */

  if (cs_glob_time_step->nt_cur > cs_glob_lagr_stat_options->nstist) {

    cs_glob_lagr_source_terms->nstits = cs_glob_lagr_stat_options->nstist;
    cs_glob_lagr_time_scheme->isttio = 1;

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
