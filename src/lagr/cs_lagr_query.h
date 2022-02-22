#ifndef __CS_LAGR_QUERY_H__
#define __CS_LAGR_QUERY_H__

/*============================================================================
 * Caller interaction with Lagrangian module.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
  Global variables
  ============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return Lagranian module status.
 *
 * \return 0 if module is not active, > 0 if active
 */
/*----------------------------------------------------------------------------*/

int
cs_lagr_model_type(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return Lagranian particle restart status.
 *
 * \return 1 if particles restart is available, 0 otherwise
 */
/*----------------------------------------------------------------------------*/

int
cs_lagr_particle_restart(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_QUERY_H__ */
