#ifndef __CS_NAVSTO_SLES_H__
#define __CS_NAVSTO_SLES_H__

/*============================================================================
 * Functions to handle SLES structures used during the resolution of the
 * Navier-Stokes system of equations
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
#include <petscksp.h>
#endif

#include "cs_iter_algo.h"
#include "cs_navsto_param.h"
#include "cs_saddle_itsol.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * User function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User-defined algorithm to solve a saddle point problem (the system is
 *        stored in a hybrid way). Please refer to cs_saddle_system_t structure
 *        and cs_saddle_block_precond_t structure definitions.
 *
 * \param[in]      nslesp  pointer to a cs_navsto_param_sles_t structure
 * \param[in]      ssys    pointer to a cs_saddle_system_t structure
 * \param[in, out] sbp     Block-preconditioner for the Saddle-point problem
 * \param[in, out] x1      array for the first part
 * \param[in, out] x2      array for the second part
 * \param[in, out] ia      pointer to a cs_iter_algo_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_navsto_sles_solve(const cs_navsto_param_sles_t    *nslesp,
                          cs_saddle_system_t              *ssys,
                          cs_saddle_block_precond_t       *sbp,
                          cs_real_t                       *x1,
                          cs_real_t                       *x2,
                          cs_iter_algo_t                  *ia);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

#if defined(HAVE_PETSC)
/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of multiplicative block preconditioner for a CG
 *
 * \param[in, out] context    pointer to optional (untyped) value or structure
 * \param[in, out] ksp_struct pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

void
cs_navsto_sles_amg_block_hook(void     *context,
                              void     *ksp_struct);
#endif

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_NAVSTO_SLES_H__ */
