#ifndef __CS_SOLVE_ALL_H__
#define __CS_SOLVE_ALL_H__

/*============================================================================
 * Solve the Navier-Stokes equations, source term convection
 * diffusion equations for scalars ... .
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

#include "cs_defs.h"

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
 * \brief Resolution of incompressible Navier Stokes, scalar transport
 *        equations... for a time step.
 *
 * \param[in]     itrale        ALE iteration number
 */
/*----------------------------------------------------------------------------*/

void
cs_solve_all(const int   itrale);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SOLVE_ALL_H__ */
