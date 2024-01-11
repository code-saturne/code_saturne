#ifndef __CS_LAGR_OPTIONS_H__
#define __CS_LAGR_OPTIONS_H__

/*============================================================================
 * Functions and types for the Lagrangian module
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

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Lagrangian module options definition.
 *
 * - default initialization
 * - read user settings
 * - check settings coherency
 * - initialize some structures relative to Lagrangian module
 *
 * \param[in]       isuite
 * \param[in]       have_thermal_model
 * \param[in]       dtref
 * \param[in, out]  iccvfg
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_options_definition(int         isuite,
                           int         have_thermal_model,
                           cs_real_t   dtref,
                           int        *iccvfg);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_OPTIONS_H__ */
