#ifndef CS_SLES_HYPRE_H
#define CS_SLES_HYPRE_H

/*============================================================================
 * Sparse Linear Equation Solvers using HYPRE
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

/*----------------------------------------------------------------------------
 * HYPRE headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_base.h"
#include "alge/cs_matrix.h"
#include "alge/cs_sles.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Solver and preconditioner types
 *----------------------------------------------------------------------------*/

typedef enum {

  /* Solver or preconditioner */

  CS_SLES_HYPRE_BOOMERAMG,     /*!< BoomerAMG (algebraic multigrid) */
  CS_SLES_HYPRE_HYBRID,        /*!< Hybrid solver */
  CS_SLES_HYPRE_ILU,           /*!< hypre-ILU (incomplete LU) */

  CS_SLES_HYPRE_BICGSTAB,      /*!< BiCGSTAB */
  CS_SLES_HYPRE_GMRES,         /*!< GMRES */
  CS_SLES_HYPRE_FLEXGMRES,     /*!< Flexible GMRES */
  CS_SLES_HYPRE_LGMRES,        /*!< LGMRES */
  CS_SLES_HYPRE_PCG,           /*!< Preconditioned Congugate Gradient */

  /* Preconditioner only */

  CS_SLES_HYPRE_EUCLID,        /*!< hypre-ILU (incomplete LU) */
  CS_SLES_HYPRE_PARASAILS,     /*!< ParaSails (sparse approximate inverse)) */

  /* End  of allowable solvers */

  CS_SLES_HYPRE_NONE           /*!< No solver or preconditioner */

} cs_sles_hypre_type_t;

/*----------------------------------------------------------------------------
 * Function pointer for user settings of a HYPRE solver setup.
 *
 * This function is called during the setup stage for a HYPRE solver.
 *
 * When first called, the solver argument is NULL, and must be created
 * using HYPRE functions.
 *
 * Note: if the context pointer is non-null, it must point to valid data
 * when the selection function is called so that value or structure should
 * not be temporary (i.e. local);
 *
 * parameters:
 *   verbosity <-- verbosity level
 *   context   <-> pointer to optional (untyped) value or structure
 *   solver    <-> handle to HYPRE solver (to be cast as HYPRE_Solver)
 *----------------------------------------------------------------------------*/

typedef void
(cs_sles_hypre_setup_hook_t) (int   verbosity,
                              void  *context,
                              void  *solver);

/* HYPRE linear solver context (opaque) */

typedef struct _cs_sles_hypre_t  cs_sles_hypre_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * Define and associate a HYPRE linear system solver for a given field or
 * equation name.
 *
 * If this system did not previously exist, it is added to the list of "known"
 * systems. Otherwise, its definition is replaced by the one defined here.
 *
 * This is a utility function: if finer control is needed, see cs_sles_define
 * and cs_sles_petsc_create.
 *
 * The associated solver required that the matrix passed to it is a HYPRE
 * matrix (see cs_matrix_set_type_hypre).
 *
 * Note that this function returns a pointer directly to the iterative solver
 * management structure. This may be used to set further options. If needed,
 * cs_sles_find may be used to obtain a pointer to the matching cs_sles_t
 * container.
 */
/*----------------------------------------------------------------------------*/

extern "C" cs_sles_hypre_t *
cs_sles_hypre_define(int                          f_id,
                     const char                  *name,
                     cs_sles_hypre_type_t         solver_type,
                     cs_sles_hypre_type_t         precond_type,
                     cs_sles_hypre_setup_hook_t  *setup_hook,
                     void                        *context);

/*----------------------------------------------------------------------------*/
/*
 * Create HYPRE linear system solver info and context.
 *
 * In case of rotational periodicity for a block (non-scalar) matrix, the
 * matrix type will be forced to MATSHELL ("shell") regardless of the option
 * used.
 */
/*----------------------------------------------------------------------------*/

extern "C" cs_sles_hypre_t *
cs_sles_hypre_create(cs_sles_hypre_type_t         solver_type,
                     cs_sles_hypre_type_t         precond_type,
                     cs_sles_hypre_setup_hook_t  *setup_hook,
                     void                        *context);

/*----------------------------------------------------------------------------*/
/*
 * Destroy HYPRE linear system solver info and context.
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_sles_hypre_destroy(void  **context);

/*----------------------------------------------------------------------------*/
/*
 * Create HYPRE linear system solver info and context based on existing info
 * and context.
 */
/*----------------------------------------------------------------------------*/

extern "C" void *
cs_sles_hypre_copy(const void  *context);

/*----------------------------------------------------------------------------*/
/*
 * Setup HYPRE linear equation solver.
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_sles_hypre_setup(void               *context,
                    const char         *name,
                    const cs_matrix_t  *a,
                    int                 verbosity);

/*----------------------------------------------------------------------------*/
/*
 * Call HYPRE linear equation solver.
 */
/*----------------------------------------------------------------------------*/

extern "C" cs_sles_convergence_state_t
cs_sles_hypre_solve(void                *context,
                    const char          *name,
                    const cs_matrix_t   *a,
                    int                  verbosity,
                    double               precision,
                    double               r_norm,
                    int                 *n_iter,
                    double              *residual,
                    const cs_real_t     *rhs,
                    cs_real_t           *vx_ini,
                    cs_real_t           *vx,
                    size_t               aux_size,
                    void                *aux_vectors);

/*----------------------------------------------------------------------------*/
/*
 * Free HYPRE linear equation solver setup context.
 *
 * This function frees resolution-related data, such as buffers and
 * preconditioning but does not free the whole context, as info used for
 * logging (especially performance data) is maintained.
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_sles_hypre_free(void  *context);

/*----------------------------------------------------------------------------*/
/*
 * Error handler for HYPRE solver.
 *
 * In case of divergence or breakdown, this error handler outputs an error
 * message. It does nothing in case the maximum iteration count is reached.
 */
/*----------------------------------------------------------------------------*/

extern "C" bool
cs_sles_hypre_error_post_and_abort(cs_sles_t                    *sles,
                                   cs_sles_convergence_state_t   state,
                                   const cs_matrix_t            *a,
                                   const cs_real_t              *rhs,
                                   cs_real_t                    *vx);

/*----------------------------------------------------------------------------*/
/*
 * Log sparse linear equation solver info.
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_sles_hypre_log(const void  *context,
                  cs_log_t     log_type);

/*----------------------------------------------------------------------------*/
/*
 * Set the max. number of iterations associated to the given HYPRE context
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_sles_hypre_set_n_max_iter(cs_sles_hypre_t   *context,
                             int                n_max_iter);

/*----------------------------------------------------------------------------*/
/*
 * Define whether the solver should run on the host or accelerated device.
 *
 * If no device is available, this setting is ignored.
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_sles_hypre_set_host_device(cs_sles_hypre_t   *context,
                              int                use_device);

/*----------------------------------------------------------------------------*/
/*
 * Query whether the solver should run on the host or accelerated device.
 */
/*----------------------------------------------------------------------------*/

extern "C" int
cs_sles_hypre_get_host_device(const cs_sles_hypre_t   *context);

/*----------------------------------------------------------------------------*/
/*
 * Print information on hypre library.
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_sles_hypre_library_info(cs_log_t  log_type);

/*----------------------------------------------------------------------------*/

#endif /* CS_SLES_HYPRE_H */
