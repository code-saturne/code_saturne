#ifndef __CS_SLES_PC_CUDA_H__
#define __CS_SLES_PC_CUDA_H__

/*============================================================================
 * Sparse Linear Equation Solver Preconditioners using CUDA
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_matrix.h"
#include "cs_sles.h"
#include "cs_sles_pc.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * User function prototypes
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Function for application of a null-preconditioner.
 *
 * In cases where it is desired that the preconditioner modify a vector
 * "in place", x_in should be set to NULL, and x_out contain the vector to
 * be modified (\f$x_{out} \leftarrow M^{-1}x_{out})\f$).
 *
 * parameters:
 *   context       <-> pointer to preconditioner context
 *   x_in          <-- input vector
 *   x_out         <-> input/output vector
 *
 * returns:
 *   preconditioner application status
 *----------------------------------------------------------------------------*/

cs_sles_pc_state_t
cs_sles_pc_cuda_apply_none(void                *context,
                           const cs_real_t     *x_in,
                           cs_real_t           *x_out);

/*----------------------------------------------------------------------------
 * Function for application of a Jacobi preconditioner.
 *
 * In cases where it is desired that the preconditioner modify a vector
 * "in place", x_in should be set to NULL, and x_out contain the vector to
 * be modified (\f$x_{out} \leftarrow M^{-1}x_{out})\f$).
 *
 * parameters:
 *   context       <-> pointer to preconditioner context
 *   x_in          <-- input vector
 *   x_out         <-> input/output vector
 *
 * returns:
 *   preconditioner application status
 *----------------------------------------------------------------------------*/

cs_sles_pc_state_t
cs_sles_pc_cuda_apply_jacobi(void                *context,
                             const cs_real_t     *x_in,
                             cs_real_t           *x_out);

/*----------------------------------------------------------------------------
 * Function for application of a polynomial preconditioner.
 *
 * In cases where it is desired that the preconditioner modify a vector
 * "in place", x_in should be set to NULL, and x_out contain the vector to
 * be modified (\f$x_{out} \leftarrow M^{-1}x_{out})\f$).
 *
 * parameters:
 *   context       <-> pointer to preconditioner context
 *   x_in          <-- input vector
 *   x_out         <-> input/output vector
 *
 * returns:
 *   preconditioner application status
 *----------------------------------------------------------------------------*/

cs_sles_pc_state_t
cs_sles_pc_cuda_apply_poly(void                *context,
                           const cs_real_t     *x_in,
                           cs_real_t           *x_out);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SLES_PC_CUDA_H__ */
