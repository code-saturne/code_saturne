#ifndef __CS_SLES_PETSC_H__
#define __CS_SLES_PETSC_H__

/*============================================================================
 * Sparse Linear Equation Solvers using PETSc
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
 * Local headers
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
 * Function pointer for user settings of a PETSc KSP solver setup.
 *
 * This function is called at the end of the setup stage for a KSP solver.
 *
 * Note that using the advanced KSPSetPostSolve and KSPSetPreSolve functions,
 * this also allows setting further function pointers for pre and post-solve
 * operations (see the PETSc documentation).
 *
 * Note: if the context pointer is non-NULL, it must point to valid data
 * when the selection function is called so that value or structure should
 * not be temporary (i.e. local);
 *
 * parameters:
 *   context <-> pointer to optional (untyped) value or structure
 *   ksp     <-> pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

typedef void
(cs_sles_petsc_setup_hook_t) (void  *context,
                              void  *ksp);

/* Iterative linear solver context (opaque) */

typedef struct _cs_sles_petsc_t  cs_sles_petsc_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * User function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Function pointer for user settings of a PETSc KSP solver setup.
 *
 * This function is called at the end of the setup stage for a KSP solver.
 *
 * Note that using the advanced KSPSetPostSolve and KSPSetPreSolve functions,
 * this also allows setting further function pointers for pre and post-solve
 * operations (see the PETSc documentation).
 *
 * Note: if the context pointer is non-NULL, it must point to valid data
 * when the selection function is called so that value or structure should
 * not be temporary (i.e. local);
 *
 * parameters:
 *   context <-> pointer to optional (untyped) value or structure
 *   ksp     <-> pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

void
cs_user_sles_petsc_hook(void   *context,
                        void   *ksp);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize PETSc if needed (calls cs_matrix_petsc_ensure_init).
 *----------------------------------------------------------------------------*/

void
cs_sles_petsc_init(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define and associate a PETSc linear system solver
 *        for a given field or equation name.
 *
 * If this system did not previously exist, it is added to the list of
 * "known" systems. Otherwise, its definition is replaced by the one
 * defined here.
 *
 * This is a utility function: if finer control is needed, see
 * \ref cs_sles_define and \ref cs_sles_petsc_create.
 *
 * In case of rotational periodicity for a block (non-scalar) matrix,
 * the matrix type will be forced to MATSHELL ("shell") regardless
 * of the option used.
 *
 * Note that this function returns a pointer directly to the iterative solver
 * management structure. This may be used to set further options.
 * If needed, \ref cs_sles_find may be used to obtain a pointer to the matching
 * \ref cs_sles_t container.
 *
 * \param[in]      f_id          associated field id, or < 0
 * \param[in]      name          associated name if f_id < 0, or NULL
 * \param[in]      matrix_type   PETSc matrix type
 * \param[in]      setup_hook    pointer to optional setup epilogue function
 * \param[in,out]  context       pointer to optional (untyped) value or
 *                               structure for setup_hook, or NULL
 *
 * \return  pointer to newly created iterative solver info object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_petsc_t *
cs_sles_petsc_define(int                          f_id,
                     const char                  *name,
                     const char                  *matrix_type,
                     cs_sles_petsc_setup_hook_t  *setup_hook,
                     void                        *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create PETSc linear system solver info and context.
 *
 * In case of rotational periodicity for a block (non-scalar) matrix,
 * the matrix type will be forced to MATSHELL ("shell") regardless
 * of the option used.
 *
 * \param[in]      matrix_type   PETSc matrix type
 * \param[in]      setup_hook    pointer to optional setup epilogue function
 * \param[in,out]  context       pointer to optional (untyped) value or
 *                               structure for setup_hook, or NULL
 *
 * \return  pointer to associated linear system object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_petsc_t *
cs_sles_petsc_create(const char                  *matrix_type,
                     cs_sles_petsc_setup_hook_t  *setup_hook,
                     void                        *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create PETSc linear system solver info and context
 *        based on existing info and context.
 *
 * \param[in]  context  pointer to reference info and context
 *                     (actual type: cs_sles_petsc_t  *)
 *
 * \return  pointer to newly created solver info object.
 *          (actual type: cs_sles_petsc_t  *)
 */
/*----------------------------------------------------------------------------*/

void *
cs_sles_petsc_copy(const void  *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy PETSc linear system solver info and context.
 *
 * \param[in, out]  context  pointer to iterative solver info and context
 *                           (actual type: cs_sles_petsc_t  **)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_petsc_destroy(void  **context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup PETSc linear equation solver.
 *
 * \param[in, out]  context    pointer to iterative solver info and context
 *                             (actual type: cs_sles_petsc_t  *)
 * \param[in]       name       pointer to system name
 * \param[in]       a          associated matrix
 * \param[in]       verbosity  associated verbosity
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_petsc_setup(void               *context,
                    const char         *name,
                    const cs_matrix_t  *a,
                    int                 verbosity);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Call PETSc linear equation solver.
 *
 * \param[in, out]  context        pointer to iterative solver info and context
 *                                 (actual type: cs_sles_petsc_t  *)
 * \param[in]       name           pointer to system name
 * \param[in]       a              matrix
 * \param[in]       verbosity      associated verbosity
 * \param[in]       precision      solver precision
 * \param[in]       r_norm         residual normalization
 * \param[out]      n_iter         number of "equivalent" iterations
 * \param[out]      residual       residual
 * \param[in]       rhs            right hand side
 * \param[in, out]  vx             system solution
 * \param[in]       aux_size       number of elements in aux_vectors (in bytes)
 * \param           aux_vectors    optional working area
 *                                 (internal allocation if NULL)
 *
 * \return  convergence state
 */
/*----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_sles_petsc_solve(void                *context,
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free PETSc linear equation solver setup context.
 *
 * This function frees resolution-related data, such as
 * buffers and preconditioning but does not free the whole context,
 * as info used for logging (especially performance data) is maintained.
 *
 * \param[in, out]  context  pointer to iterative solver info and context
 *                           (actual type: cs_sles_petsc_t  *)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_petsc_free(void  *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Error handler for PETSc solver.
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
cs_sles_petsc_error_post_and_abort(cs_sles_t                    *sles,
                                   cs_sles_convergence_state_t   state,
                                   const cs_matrix_t            *a,
                                   const cs_real_t              *rhs,
                                   cs_real_t                    *vx);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log sparse linear equation solver info.
 *
 * \param[in]  context   pointer to iterative solver info and context
 *                       (actual type: cs_sles_petsc_t  *)
 * \param[in]  log_type  log type
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_petsc_log(const void  *context,
                  cs_log_t     log_type);

/*----------------------------------------------------------------------------
 * \brief Output the settings of a KSP structure
 *
 * \param[in]  ksp     Krylov SubSpace structure
 *----------------------------------------------------------------------------*/

void
cs_sles_petsc_log_setup(void  *ksp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the parameters driving the termination of an iterative process
 *        associated to a KSP structure
 *
 * \param[in, out] context   pointer to iterative solver info and context
 *                           (actual type: cs_sles_petsc_t  *)
 * \param[in]      rtol      relative tolerance
 * \param[in]      atol      absolute tolerance
 * \param[in]      dtol      divergence tolerance
 * \param[in]      max_it    max. number of iterations
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_petsc_set_cvg_criteria(const void      *context,
                               double           rtol,
                               double           atol,
                               double           dtol,
                               int              max_it);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return matrix type associated with PETSc linear system solver
 *        info and context.
 *
 * \param[in, out]  context  pointer to iterative solver info and context
 *                           (actual type: cs_sles_petsc_t  **)
 *
 * \return  pointer to matrix type name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_sles_petsc_get_mat_type(void  *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print information on PETSc library.
 *
 * \param[in]  log_type  log type
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_petsc_library_info(cs_log_t  log_type);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SLES_PETSC_H__ */
