#ifndef __CS_CF_BOUNDARY_CONDITIONS_H__
#define __CS_CF_BOUNDARY_CONDITIONS_H__

/*============================================================================
 * Compressible flow boundary conditions.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Automatic boundary condition for compressible flows
 *
 * \param[in]  bc_type  type of boundary for each face
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_boundary_conditions(int  bc_type[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Allocate boundary flux indicator arrays.
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_boundary_conditions_init(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Provide access to boundary face indicator array of convection flux
 *        - 0 upwind scheme
 *        - 1 imposed flux
 */
/*----------------------------------------------------------------------------*/

int *
cs_cf_boundary_conditions_get_icvfli(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Provide access to imposed thermal flux indicator at the boundary
 *        (some boundary contributions of the total energy eq. have to be
 *         cancelled)
 */
/*----------------------------------------------------------------------------*/

int *
cs_cf_boundary_conditions_get_ifbet(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Prepare (reset) condition coefficients specific to compressible flows.
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_boundary_conditions_reset(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CF_BOUNDARY_CONDITIONS_H__ */
