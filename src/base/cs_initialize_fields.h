#ifndef __CS_INITIALIZE_FIELDS_H__
#define __CS_INITIALIZE_FIELDS_H__

/*============================================================================
 * Initialize fields.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computed variable and property initialization.
 *
 * The time step, the indicator of wall distance computation are also
 * initialized just before reading a restart file or use the user
 * initializations.
 */
/*----------------------------------------------------------------------------*/

void
cs_initialize_fields_stage_0(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize variable, time step, and wall distance fields.
 */
/*----------------------------------------------------------------------------*/

void
cs_initialize_fields_stage_1(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_INITIALIZE_FIELDS_H__ */
