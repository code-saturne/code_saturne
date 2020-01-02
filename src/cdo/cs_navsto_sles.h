#ifndef __CS_NAVSTO_SLES_H__
#define __CS_NAVSTO_SLES_H__

/*============================================================================
 * Routines to handle SLES structure and PETSc interfaces for solving the
 * Navier-Stokes system of equations
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#if defined(HAVE_PETSC)
#include <petscversion.h>
#include <petscksp.h>
#endif

#include "cs_equation_param.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

#if defined(HAVE_PETSC)
/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of multiplicative block preconditioner for a CG
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] a        pointer to PETSc Matrix context
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

void
cs_navsto_sles_amg_block_hook(void     *context,
                              Mat       a,
                              KSP       ksp);
#endif

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_NAVSTO_SLES_H__ */
