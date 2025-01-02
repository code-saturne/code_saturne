#ifndef __CS_BOUNDARY_CONDITIONS_COUPLING_H__
#define __CS_BOUNDARY_CONDITIONS_COUPLING_H__

/*============================================================================
 * Update boundary conditions for thermal field.
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Input data for 1D wall thermal coupling
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_coupling_t_in(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Input data for 1D wall thermal coupling
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief output data for 1D wall thermal coupling
 *
 * \param[in, out]  hbord  exchange coefficients for boundary
 * \param[in, out]  tbord  boundary temperature
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_coupling_t_out(cs_real_t  hbord[],
                                      cs_real_t  tbord[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BOUNDARY_CONDITIONS_COUPLING_H__ */
