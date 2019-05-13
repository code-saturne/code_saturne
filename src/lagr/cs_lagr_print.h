#ifndef __CS_LAGR_PRINT_H__
#define __CS_LAGR_PRINT_H__

/*============================================================================
 * Functions and types for lagrangian specific prints
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
 * \brief Write Lagrangian particle info file
 *
 * This file logs, for each time step:
 *  - number of particles in the domain
 *  - number of entered particles
 *  - number of exited particles
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_print(cs_real_t ttcabs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Close Lagrangian particle info file
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_print_finalize(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_PRINT_H__ */
