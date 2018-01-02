#ifndef __CS_LAGR_TRACKING_H__
#define __CS_LAGR_TRACKING_H__

/*============================================================================
 * Functions and types for the Lagrangian module
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include "cs_interface.h"

#include "assert.h"

#include "cs_lagr_particle.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize particle tracking subsystem
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_tracking_initialize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply one particle movement step.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_tracking_particle_movement(const cs_real_t  visc_length[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize Lagrangian module.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_tracking_finalize(void);

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculates the distance to the wall y+.
 */
/*----------------------------------------------------------------------------*/

void
_test_wall_cell(const void                     *particle,
                const cs_lagr_attribute_map_t  *p_am,
                const cs_real_t                 visc_length[],
                cs_real_t                      *yplus,
                cs_lnum_t                      *face_id);
/*----------------------------------------------------------------------------*/


END_C_DECLS

#endif /* __CS_LAGR_TRACKING_H__ */
