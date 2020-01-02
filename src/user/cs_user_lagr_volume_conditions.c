/*============================================================================
 * Lagrangian volume injection definitions.
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
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define particle volume conditions.
 *
 * This is used for the definition of volume injections,
 * based on predefined volume zones (\ref cs_zone_t).
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_volume_conditions(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handling of a particle interaction with a interior face of type
 *        \ref CS_LAGR_BC_USER.
 *
 * \param[in, out]  particles       pointer to particle set
 * \param[in]       p_id            particle id
 * \param[in]       face_id         interior face id
 * \param[in]       face_norm       unit face (or face subdivision) normal
 * \param[in]       c_intersect     coordinates of intersection with the face
 * \param[in]       t_intersect     relative distance (in [0, 1]) of the
 *                                  intersection point with the face relative
 *                                  to the initial trajectory segment
 * \param[in, out]  tracking_state  particle tracking state
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_user_internal_interaction(cs_lagr_particle_set_t    *particles,
                                  cs_lnum_t                  p_id,
                                  cs_lnum_t                  face_id,
                                  const cs_real_t            face_norm[3],
                                  const cs_real_t            c_intersect[3],
                                  cs_real_t                  t_intersect,
                                  cs_lagr_tracking_state_t  *tracking_state)
{

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
