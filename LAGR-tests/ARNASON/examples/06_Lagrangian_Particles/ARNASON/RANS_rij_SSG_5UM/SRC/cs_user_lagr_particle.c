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
  /* Loop on new particles
     --------------------- */

  for (cs_lnum_t p_id = particle_range[0]; p_id < particle_range[1]; p_id++) {

    /* particles are injected along pipe axis (x=0,y=0) */

    cs_real_t *part_coords = cs_lagr_particles_attr(particles,
                                                    p_id,
                                                    CS_LAGR_COORDS);

    part_coords[0] = 0.; // x = 0
    part_coords[1] = 0.; // y = 0
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
