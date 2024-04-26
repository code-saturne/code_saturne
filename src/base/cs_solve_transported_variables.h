#ifndef __CS_SOLVE_TRANSPORTED_VARIABLES__
#define __CS_SOLVE_TRANSPORTED_VARIABLES__

/*============================================================================
 * Resolution of source term convection diffusion equations
 * for scalars in a time step.
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

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \file cs_solve_transported_variables.c
 *
 * \brief Resolution of source term convection diffusion equations
 *        for scalars in a time step.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*
 * \param[in]     nscal         total number of scalars
 * \param[in]     iterns        Navier-Stokes iteration number
 */
/*----------------------------------------------------------------------------*/

void
cs_solve_transported_variables(int nscal,
                               int iterns);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SOLVE_TRANSPORTED_VARIABLES__ */
