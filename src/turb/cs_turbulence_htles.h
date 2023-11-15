#ifndef __CS_TURBULENCE_HTLES_H__
#define __CS_TURBULENCE_HTLES_H__

/*============================================================================
 * Hybrid Temporal LES turbulence model.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve the HTLES method.
 *
 * Solve the HTLES for one time step.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_htles(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialization of the fields for HTLES and Dmax calculation
 */
/*----------------------------------------------------------------------------*/

void
cs_htles_initialization(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TURBULENCE_HTLES_H__ */
