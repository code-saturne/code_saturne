#ifndef __CS_MULTIGRID_SMOOTHER_H__
#define __CS_MULTIGRID_SMOOTHER_H__

/*============================================================================
 * Sparse Linear Equation Solvers
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

#include "cs_base.h"
#include "cs_halo_perio.h"
#include "cs_matrix.h"
#include "cs_sles.h"
#include "cs_sles_it.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Iterative linear solver context (opaque) */

typedef struct _cs_multigrid_smoother_t  cs_multigrid_smoother_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create iterative sparse linear system solver info and context used as
 *        a smoother in a multigrid solver
 *
 * parameters:
 * \param[in]  smoother_type   type of smoother (CG, Jacobi, ...)
 * \param[in]  poly_degree     preconditioning polynomial degree
 *                             (0: diagonal; -1: non-preconditioned;
 *                             see \ref sles_it for details)
 * \param[in]  n_iter          number of iterations to perform
 *
 * \return a pointer to newly created smoother info object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_it_t *
cs_multigrid_smoother_create(cs_sles_it_type_t     smoother_type,
                             int                   poly_degree,
                             int                   n_iter);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy iterative sparse linear system solver info and context used
 *        as a smoother in a multigrid solver
 *
 * \param[in, out]  context   pointer to iterative sparse linear solver info
 *                            (actual type: cs_sles_it_t  **)
 */
/*----------------------------------------------------------------------------*/

void
cs_multigrid_smoother_destroy(void  **context);

/*----------------------------------------------------------------------------
 * Setup iterative sparse linear equation solver.
 *
 * parameters:
 *   context   <-> pointer to iterative sparse linear solver info
 *                 (actual type: cs_multigrid_smoother_t  *)
 *   name      <-- pointer to system name
 *   a         <-- associated matrix
 *   verbosity <-- verbosity level
 *----------------------------------------------------------------------------*/

void
cs_multigrid_smoother_setup(void               *context,
                            const char         *name,
                            const cs_matrix_t  *a,
                            int                 verbosity);

/*----------------------------------------------------------------------------
 * Call iterative sparse linear equation solver.
 *
 * parameters:
 *   context       <-> pointer to iterative sparse linear solver info
 *                     (actual type: cs_multigrid_smoother_t  *)
 *   name          <-- pointer to system name
 *   a             <-- matrix
 *   verbosity     <-- verbosity level
 *   precision     <-- solver precision
 *   r_norm        <-- residual normalization
 *   n_iter        --> number of iterations
 *   residual      --> residual
 *   rhs           <-- right hand side
 *   vx            <-> system solution
 *   aux_size      <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors   --- optional working area (internal allocation if NULL)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_multigrid_smoother_solve(void                *context,
                            const char          *name,
                            const cs_matrix_t   *a,
                            int                  verbosity,
                            double               precision,
                            double               r_norm,
                            int                 *n_iter,
                            double              *residual,
                            const cs_real_t     *rhs,
                            cs_real_t           *vx,
                            size_t               aux_size,
                            void                *aux_vectors);

/*----------------------------------------------------------------------------
 * Free iterative sparse linear equation solver setup context.
 *
 * This function frees resolution-related data, such as
 * buffers and preconditioning but does not free the whole context,
 * as info used for logging (especially performance data) is maintained.

 * parameters:
 *   context <-> pointer to iterative sparse linear solver info
 *               (actual type: cs_multigrid_smoother_t  *)
 *----------------------------------------------------------------------------*/

void
cs_multigrid_smoother_free(void  *context);

/*----------------------------------------------------------------------------
 * Log sparse linear equation solver info.
 *
 * parameters:
 *   context  <-> pointer to iterative sparse linear solver info
 *                (actual type: cs_multigrid_smoother_t  *)
 *   log_type <-- log type
 *----------------------------------------------------------------------------*/

void
cs_multigrid_smoother_log(const void  *context,
                          cs_log_t     log_type);

/*----------------------------------------------------------------------------
 * Return iterative solver type.
 *
 * parameters:
 *   context <-- pointer to iterative solver info and context
 *
 * returns:
 *   selected solver type
 *----------------------------------------------------------------------------*/

cs_sles_it_type_t
cs_multigrid_smoother_get_type(const cs_multigrid_smoother_t  *context);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MULTIGRID_SMOOTHER_H__ */
