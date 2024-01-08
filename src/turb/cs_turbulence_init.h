#ifndef __CS_TURBULENCE_INIT_H__
#define __CS_TURBULENCE_INIT_H__

/*============================================================================
 * Turbulence variables initialization for various models.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

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
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief First pass at initialization of turbulence variables.
 *
 * If the reference velocity is negative (incorrect or not initialized),
 * values of k, Rij, epsilon, and omega are set to -10*cs_math_big_r.
 * We will later test if those values were modified by user initialization
 * or reading a restart file.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_init_by_ref_quantities(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Clipping of turbulence variables initialization
 *
 * If the user has prescribed "correct" values (in the sense that k, eps,
 * Rii, ...  > 0), we clip these values if needed for consistency.
 *
 * If the prescribed values are manifestly incorrect (negative values),
 * we assume this is an error, print warnings, and return a positive
 * error count.
 *
 * The same logic is applied to computation restarts.
 */
/*----------------------------------------------------------------------------*/

int
cs_turbulence_init_clip_and_verify(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TURBULENCE_INIT_H__ */
