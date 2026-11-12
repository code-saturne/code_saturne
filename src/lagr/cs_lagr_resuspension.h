#ifndef CS_LAGR_RESUSPENSION_H
#define CS_LAGR_RESUSPENSION_H

/*============================================================================
 * Methods for particle resuspension
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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
 * \brief Calculation of the adhesion force of the particles lying on a
 *        monolayer deposit whose adhesion force is not known yet.
 *
 * To be called before the trajectory step with the cell-wise integration,
 * which integrates the particles injected in the current time step.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_resuspension_init_deposited(void);


/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the particle resuspension
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_resuspension(void);

/*----------------------------------------------------------------------------*/

#endif /* CS_LAGR_RESUSPENSION_H */
