#ifndef __CS_SLES_HYPRE_H__
#define __CS_SLES_HYPRE_H__

/*============================================================================
 * Sparse Linear Equation Solvers using HYPRE
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
 * HYPRE headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_matrix.h"
#include "cs_sles.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

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
 * Note: if the context pointer is non-NULL, it must point to valid data
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
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define and associate a HYPRE linear system solver
 *        for a given field or equation name.
 *
 * If this system did not previously exist, it is added to the list of
 * "known" systems. Otherwise, its definition is replaced by the one
 * defined here.
 *
 * This is a utility function: if finer control is needed, see
 * \ref cs_sles_define and \ref cs_sles_petsc_create.
 *
 * The associated solver required that the matrix passed to it is a HYPRE
 * matrix (see cs_matrix_set_type_hypre).
 *
 * Note that this function returns a pointer directly to the iterative solver
 * management structure. This may be used to set further options.
 * If needed, \ref cs_sles_find may be used to obtain a pointer to the matching
 * \ref cs_sles_t container.
 *
 * \param[in]      f_id          associated field id, or < 0
 * \param[in]      name          associated name if f_id < 0, or NULL
 * \param[in]      solver_type   HYPRE solver type
 * \param[in]      precond_type  HYPRE preconditioner type
 * \param[in]      setup_hook    pointer to optional setup epilogue function
 * \param[in,out]  context       pointer to optional (untyped) value or
 *                               structure for setup_hook, or NULL
 *
 * \return  pointer to newly created iterative solver info object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_hypre_t *
cs_sles_hypre_define(int                          f_id,
                     const char                  *name,
                     cs_sles_hypre_type_t         solver_type,
                     cs_sles_hypre_type_t         precond_type,
                     cs_sles_hypre_setup_hook_t  *setup_hook,
                     void                        *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create HYPRE linear system solver info and context.
 *
 * In case of rotational periodicity for a block (non-scalar) matrix,
 * the matrix type will be forced to MATSHELL ("shell") regardless
 * of the option used.
 *
 * \param[in]  solver_type   HYPRE solver type
 * \param[in]  precond_type  HYPRE preconditioner type
 * \param[in]  setup_hook    pointer to optional setup epilogue function
 * \param[in]  context       pointer to optional (untyped) value or structure
 *                           for setup_hook, or NULL
 *
 * \return  pointer to newly created linear system object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_hypre_t *
cs_sles_hypre_create(cs_sles_hypre_type_t         solver_type,
                     cs_sles_hypre_type_t         precond_type,
                     cs_sles_hypre_setup_hook_t  *setup_hook,
                     void                        *context);

/*----------------------------------------------------------------------------
 * Destroy HYPRE linear system solver info and context.
 *
 * parameters:
 *   context  <-> pointer to HYPRE linear solver info
 *                (actual type: cs_sles_hypre_t  **)
 *----------------------------------------------------------------------------*/

void
cs_sles_hypre_destroy(void  **context);

/*----------------------------------------------------------------------------
 * Create HYPRE linear system solver info and context
 * based on existing info and context.
 *
 * parameters:
 *   context <-- pointer to reference info and context
 *               (actual type: cs_sles_hypre_t  *)
 *
 * returns:
 *   pointer to newly created solver info object
 *   (actual type: cs_sles_hypre_t  *)
 *----------------------------------------------------------------------------*/

void *
cs_sles_hypre_copy(const void  *context);

/*----------------------------------------------------------------------------
 * Setup HYPRE linear equation solver.
 *
 * parameters:
 *   context   <-> pointer to HYPRE linear solver info
 *                 (actual type: cs_sles_hypre_t  *)
 *   name      <-- pointer to system name
 *   a         <-- associated matrix
 *   verbosity <-- verbosity level
 *----------------------------------------------------------------------------*/

void
cs_sles_hypre_setup(void               *context,
                    const char         *name,
                    const cs_matrix_t  *a,
                    int                 verbosity);

/*----------------------------------------------------------------------------
 * Call HYPRE linear equation solver.
 *
 * parameters:
 *   context       <-> pointer to HYPRE linear solver info
 *                     (actual type: cs_sles_hypre_t  *)
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
cs_sles_hypre_solve(void                *context,
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
 * Free HYPRE linear equation solver setup context.
 *
 * This function frees resolution-related data, such as
 * buffers and preconditioning but does not free the whole context,
 * as info used for logging (especially performance data) is maintained.

 * parameters:
 *   context <-> pointer to HYPRE linear solver info
 *               (actual type: cs_sles_hypre_t  *)
 *----------------------------------------------------------------------------*/

void
cs_sles_hypre_free(void  *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Error handler for HYPRE solver.
 *
 * In case of divergence or breakdown, this error handler outputs an error
 * message
 * It does nothing in case the maximum iteration count is reached.

 * \param[in, out]  sles           pointer to solver object
 * \param[in]       state          convergence state
 * \param[in]       a              matrix
 * \param[in]       rhs            right hand side
 * \param[in, out]  vx             system solution
 *
 * \return  false (do not attempt new solve)
 */
/*----------------------------------------------------------------------------*/

bool
cs_sles_hypre_error_post_and_abort(cs_sles_t                    *sles,
                                   cs_sles_convergence_state_t   state,
                                   const cs_matrix_t            *a,
                                   const cs_real_t              *rhs,
                                   cs_real_t                    *vx);

/*----------------------------------------------------------------------------
 * Log sparse linear equation solver info.
 *
 * parameters:
 *   context  <-> pointer to HYPRE linear solver info
 *                (actual type: cs_sles_hypre_t  *)
 *   log_type <-- log type
 *----------------------------------------------------------------------------*/

void
cs_sles_hypre_log(const void  *context,
                  cs_log_t     log_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the max. number of iterations associated to the given HYPRE
 *        contrext
 *
 * \param[in,out]  context       pointer to HYPRE linear solver info
 * \param[in]      n_max_iter    max. number of iterations
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_hypre_set_n_max_iter(cs_sles_hypre_t   *context,
                             int                n_max_iter);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define whether the solver should run on the host or
 *        accelerated device.
 *
 * If no device is available, this setting is ignored.
 *
 * \param[in,out]  context       pointer to HYPRE linear solver info
 * \param[in]      use_device    0 for host, 1 for device (GPU)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_hypre_set_host_device(cs_sles_hypre_t   *context,
                              int                use_device);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Query whether the solver should run on the host or accelerated device.
 *
 * \param[in,out]  context       pointer to HYPRE linear solver info
 *
 * \return   0 for host, 1 for device (GPU)
 */
/*----------------------------------------------------------------------------*/

int
cs_sles_hypre_get_host_device(const cs_sles_hypre_t   *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print information on hypre library.
 *
 * \param[in]  log_type  log type
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_hypre_library_info(cs_log_t  log_type);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SLES_HYPRE_H__ */
