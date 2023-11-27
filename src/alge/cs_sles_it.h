#ifndef __CS_SLES_IT_H__
#define __CS_SLES_IT_H__

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

/*----------------------------------------------------------------------------
 * Solver types
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_SLES_PCG,                 /*!< Preconditioned conjugate gradient */
  CS_SLES_FCG,                 /*!< Preconditions flexible conjugate gradient,
                                    described in \cite Notay:2015 */
  CS_SLES_IPCG,                /*!< Preconditions inexact conjugate gradient */
  CS_SLES_JACOBI,              /*!< Jacobi */
  CS_SLES_BICGSTAB,            /*!< Preconditioned BiCGstab
                                    (biconjugate gradient stabilized) */
  CS_SLES_BICGSTAB2,           /*!< Preconditioned BiCGstab2 */
  CS_SLES_GCR,                 /*!< Generalized conjugate residual  */
  CS_SLES_GMRES,               /*!< Preconditioned GMRES
                                    (generalized minimal residual) */
  CS_SLES_P_GAUSS_SEIDEL,      /*!< Process-local Gauss-Seidel */
  CS_SLES_P_SYM_GAUSS_SEIDEL,  /*!< Process-local symmetric Gauss-Seidel */
  CS_SLES_PCR3,                /*!< 3-layer conjugate residual */
  CS_SLES_USER_DEFINED,        /*!< User-defined iterative solver */

  CS_SLES_N_IT_TYPES,          /*!< Number of resolution algorithms
                                    excluding smoother only */

  CS_SLES_TS_F_GAUSS_SEIDEL,   /*!< Truncated forward Gauss-Seidel smoother */
  CS_SLES_TS_B_GAUSS_SEIDEL,   /*!< Truncated backward Gauss-Seidel smoother */

  CS_SLES_N_SMOOTHER_TYPES     /*!< Number of resolution algorithms
                                    including smoother only */

} cs_sles_it_type_t;

/* Iterative linear solver context (opaque) */

typedef struct _cs_sles_it_t  cs_sles_it_t;

/* Forward type declarations */

typedef struct _cs_sles_it_convergence_t  cs_sles_it_convergence_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Names for iterative solver types */

extern const char *cs_sles_it_type_name[];

/*=============================================================================
 * User function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using a user-defined iterative solver
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- matrix
 *   diag_block_size <-- diagonal block size (unused here)
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_user_sles_it_solver(cs_sles_it_t              *c,
                       const cs_matrix_t         *a,
                       cs_lnum_t                  diag_block_size,
                       cs_sles_it_convergence_t  *convergence,
                       const cs_real_t           *rhs,
                       cs_real_t                 *restrict vx,
                       size_t                     aux_size,
                       void                      *aux_vectors);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define and associate an iterative sparse linear system solver
 * for a given field or equation name.
 *
 * If this system did not previously exist, it is added to the list of
 * "known" systems. Otherwise, its definition is replaced by the one
 * defined here.
 *
 * This is a utility function: if finer control is needed, see
 * cs_sles_define() and cs_sles_it_create().
 *
 * Note that this function returns a pointer directly to the iterative solver
 * management structure. This may be used to set further options,
 * for example using cs_sles_set_plot_options(). If needed, cs_sles_find()
 * may be used to obtain a pointer to the matching cs_sles_t container.
 *
 * parameters:
 *   f_id         <-- associated field id, or < 0
 *   name         <-- associated name if f_id < 0, or NULL
 *   solver_type  <-- type of solver (PCG, Jacobi, ...)
 *   poly_degree  <-- preconditioning polynomial degree
 *                    (0: diagonal; -1: non-preconditioned)
 *   n_max_iter   <-- maximum number of iterations
 *
 * returns:
 *   pointer to newly created iterative solver info object.
 *----------------------------------------------------------------------------*/

cs_sles_it_t *
cs_sles_it_define(int                 f_id,
                  const char         *name,
                  cs_sles_it_type_t   solver_type,
                  int                 poly_degree,
                  int                 n_max_iter);

/*----------------------------------------------------------------------------
 * Create iterative sparse linear system solver info and context.
 *
 * parameters:
 *   solver_type  <-- type of solver (PCG, Jacobi, ...)
 *   poly_degree  <-- preconditioning polynomial degree
 *                    (0: diagonal; -1: non-preconditioned)
 *   n_max_iter   <-- maximum number of iterations
 *   update_stats <-- automatic solver statistics indicator
 *
 * returns:
 *   pointer to newly created solver info object.
 *----------------------------------------------------------------------------*/

cs_sles_it_t *
cs_sles_it_create(cs_sles_it_type_t   solver_type,
                  int                 poly_degree,
                  int                 n_max_iter,
                  bool                update_stats);

/*----------------------------------------------------------------------------
 * Destroy iterative sparse linear system solver info and context.
 *
 * parameters:
 *   context  <-> pointer to iterative sparse linear solver info
 *                (actual type: cs_sles_it_t  **)
 *----------------------------------------------------------------------------*/

void
cs_sles_it_destroy(void  **context);

/*----------------------------------------------------------------------------
 * Create iterative sparse linear system solver info and context
 * based on existing info and context.
 *
 * parameters:
 *   context <-- pointer to reference info and context
 *               (actual type: cs_sles_it_t  *)
 *
 * returns:
 *   pointer to newly created solver info object
 *   (actual type: cs_sles_it_t  *)
 *----------------------------------------------------------------------------*/

void *
cs_sles_it_copy(const void  *context);

/*----------------------------------------------------------------------------
 * Setup iterative sparse linear equation solver.
 *
 * parameters:
 *   context   <-> pointer to iterative sparse linear solver info
 *                 (actual type: cs_sles_it_t  *)
 *   name      <-- pointer to system name
 *   a         <-- associated matrix
 *   verbosity <-- verbosity level
 *----------------------------------------------------------------------------*/

void
cs_sles_it_setup(void               *context,
                 const char         *name,
                 const cs_matrix_t  *a,
                 int                 verbosity);

/*----------------------------------------------------------------------------
 * Call iterative sparse linear equation solver.
 *
 * parameters:
 *   context       <-> pointer to iterative sparse linear solver info
 *                     (actual type: cs_sles_it_t  *)
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
cs_sles_it_solve(void                *context,
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
 *               (actual type: cs_sles_it_t  *)
 *----------------------------------------------------------------------------*/

void
cs_sles_it_free(void  *context);

/*----------------------------------------------------------------------------
 * Log sparse linear equation solver info.
 *
 * parameters:
 *   context  <-> pointer to iterative sparse linear solver info
 *                (actual type: cs_sles_it_t  *)
 *   log_type <-- log type
 *----------------------------------------------------------------------------*/

void
cs_sles_it_log(const void  *context,
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
cs_sles_it_get_type(const cs_sles_it_t  *context);

/*----------------------------------------------------------------------------
 * Return the initial residual for the previous solve operation with a solver.
 *
 * This is useful for convergence tests when this solver is used as
 * a preconditioning smoother.
 *
 * This operation is only valid between calls to cs_sles_it_setup()
 * (or cs_sles_it_solve()) and cs_sles_it_free().
 * It returns -1 otherwise.
 *
 * parameters:
 *   context <-- pointer to iterative solver info and context
 *
 * returns:
 *   initial residual from last call to \ref cs_sles_solve with this solver
 *----------------------------------------------------------------------------*/

double
cs_sles_it_get_last_initial_residual(const cs_sles_it_t  *context);

/*----------------------------------------------------------------------------
 * Return a preconditioner context for an iterative sparse linear
 * equation solver.
 *
 * This allows modifying parameters of a non default (Jacobi or polynomial)
 * preconditioner.
 *
 * parameters:
 *   context <-- pointer to iterative solver info and context
 *
 * returns:
 *   pointer to preconditoner context
 *----------------------------------------------------------------------------*/

cs_sles_pc_t  *
cs_sles_it_get_pc(cs_sles_it_t  *context);

/*----------------------------------------------------------------------------
 * Assign a preconditioner to an iterative sparse linear equation
 * solver, transfering its ownership to to solver context.
 *
 * This allows assigning a non default (Jacobi or polynomial) preconditioner.
 *
 * The input pointer is set to NULL to make it clear the caller does not
 * own the preconditioner anymore, though the context can be accessed using
 * cs_sles_it_get_cp().
 *
 * parameters:
 *   context <->  pointer to iterative solver info and context
 *   pc      <->  pointer to preconditoner context
 *----------------------------------------------------------------------------*/

void
cs_sles_it_transfer_pc(cs_sles_it_t     *context,
                       cs_sles_pc_t    **pc);

/*----------------------------------------------------------------------------
 * Copy options from one iterative sparse linear system solver info
 * and context to another.
 *
 * Optional plotting contexts are shared between the source and destination
 * contexts.
 *
 * Preconditioner settings are to be handled separately.
 *
 * parameters:
 *   src  <-- pointer to source info and context
 *   dest <-> pointer to destination info and context
 *----------------------------------------------------------------------------*/

void
cs_sles_it_transfer_parameters(const cs_sles_it_t  *src,
                               cs_sles_it_t        *dest);

/*----------------------------------------------------------------------------
 * Associate a similar info and context object with which some setup
 * data may be shared.
 *
 * This is especially useful for sharing preconditioning data between
 * similar solver contexts (for example ascending and descending multigrid
 * smoothers based on the same matrix).
 *
 * For preconditioning data to be effectively shared, cs_sles_it_setup()
 * (or cs_sles_it_solve()) must be called on "shareable" before being
 * called on "context" (without cs_sles_it_free() being called in between,
 * of course).
 *
 * It is the caller's responsibility to ensure the context is not used
 * for a cs_sles_it_setup() or cs_sles_it_solve() operation  after the
 * shareable object has been destroyed (normally by cs_sles_it_destroy()).
 *
 * parameters:
 *   context   <-> pointer to iterative sparse linear system solver info
 *   shareable <-- pointer to iterative solver info and context
 *                 whose context may be shared
 *----------------------------------------------------------------------------*/

void
cs_sles_it_set_shareable(cs_sles_it_t        *context,
                         const cs_sles_it_t  *shareable);

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set MPI communicator for global reductions.
 *
 * The system is solved only on ranks with a non-NULL communicator or
 * if the caller communicator has less than 2 ranks. convergence info
 * is broadcast across the caller communicator.
 *
 * \param[in, out]  context      pointer to iterative solver info and context
 * \param[in]       comm         MPI communicator for solving
 * \param[in]       caller_comm  MPI communicator of caller
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_set_mpi_reduce_comm(cs_sles_it_t  *context,
                               MPI_Comm       comm,
                               MPI_Comm       caller_comm);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Assign ordering to iterative solver.
 *
 * The solver context takes ownership of the order array (i.e. it will
 * handle its later deallocation).
 *
 * This is useful only for Block Gauss-Seidel.
 *
 * parameters:
 *   context <-> pointer to iterative solver info and context
 *   order   <-> pointer to ordering array
 *----------------------------------------------------------------------------*/

void
cs_sles_it_assign_order(cs_sles_it_t   *context,
                        cs_lnum_t     **order);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the threshold value under which a breakdown happens in
 *        solvers like BiCGStab or BiCGStab2
 *
 * \return the value of the threshold
 */
/*----------------------------------------------------------------------------*/

double
cs_sles_it_get_breakdown_threshold(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the threshold value under which a breakdown happens in
 *        solvers like BiCGStab or BiCGStab2
 *
 * \param[in]  threshold  value of the threshold
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_set_breakdown_threshold(double  threshold);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define convergence level under which the fallback to another
 *        solver may be used if applicable.
 *
 * Currently, this mechanism is used by default for solvers which may exhibit
 * breakdown, such as BiCGstab and 3-layer conjugate residual solvers, which
 * may fall back to a a more robust preconditioned GMRES solver.
 *
 * For those solvers, the default threshold is \ref CS_SLES_MAX_ITERATION,
 * meaning that reaching breakdown will lead to the use of the
 * fallback mechanism.
 *
 * \param[in, out]  context    pointer to iterative solver info and context
 * \param[in]       threshold  convergence level under which fallback is used
 * \param[in]       n_iter_max  maximum number of iterations fo fallback solver
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_set_fallback_threshold(cs_sles_it_t                 *context,
                                  cs_sles_convergence_state_t   threshold,
                                  int                           n_iter_max);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the number of iterations to be done before restarting the
 *        solver. Useful only for GCR or GMRES algorithms.
 *
 * \param[in, out]  context   pointer to iterative solver info and context
 * \param[in]       interval  convergence level under which fallback is used
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_set_restart_interval(cs_sles_it_t  *context,
                                int            interval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the max. number of iterations before stopping the algorithm
 *
 * \param[in, out]  context     pointer to iterative solver info and context
 * \param[in]       n_max_iter  max. number of iterations
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_set_n_max_iter(cs_sles_it_t  *context,
                          int            n_max_iter);

/*----------------------------------------------------------------------------
 * Query mean number of rows under which Conjugate Gradient algorithm
 * uses the single-reduction variant.
 *
 * The single-reduction variant requires only one parallel sum per
 * iteration (instead of 2), at the cost of additional vector operations,
 * so it tends to be more expensive when the number of matrix rows per
 * MPI rank is high, then becomes cheaper when the MPI latency cost becomes
 * more significant.
 *
 * This option is ignored for non-parallel runs, so 0 is returned.
 *
 * return:
 *   mean number of rows per active rank under which the
 *   single-reduction variant will be used
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_sles_it_get_pcg_single_reduction(void);

/*----------------------------------------------------------------------------
 * Set mean number of rows under which Conjugate Gradient algorithm
 * should use the single-reduction variant.
 *
 * The single-reduction variant requires only one parallel sum per
 * iteration (instead of 2), at the cost of additional vector operations,
 * so it tends to be more expensive when the number of matrix rows per
 * MPI rank is high, then becomes cheaper when the MPI latency cost becomes
 * more significant.
 *
 * This option is ignored for non-parallel runs.
 *
 * parameters:
 *   threshold <-- mean number of rows per active rank under which the
 *                 single-reduction variant will be used
 *----------------------------------------------------------------------------*/

void
cs_sles_it_set_pcg_single_reduction(cs_lnum_t  threshold);

/*----------------------------------------------------------------------------
 * Log the current global settings relative to parallelism.
 *----------------------------------------------------------------------------*/

void
cs_sles_it_log_parallel_options(void);

/*----------------------------------------------------------------------------
 * Error handler for iterative sparse linear equation solver.
 *
 * In case of divergence or breakdown, this error handler outputs
 * postprocessing data to assist debugging, then aborts the run.
 * It does nothing in case the maximum iteration count is reached.
 *
 * parameters:
 *   sles          <-> pointer to solver object
 *   state         <-- convergence state
 *   a             <-- matrix
 *   rhs           <-- right hand side
 *   vx            <-> system solution
 *
 * returns:
 *   false (do not attempt new solve)
 *----------------------------------------------------------------------------*/

bool
cs_sles_it_error_post_and_abort(cs_sles_t                    *sles,
                                cs_sles_convergence_state_t   state,
                                const cs_matrix_t            *a,
                                const cs_real_t              *rhs,
                                cs_real_t                    *vx);

/*----------------------------------------------------------------------------
 * Set plotting options for an iterative sparse linear equation solver.
 *
 * parameters:
 *   context       <-> pointer to iterative solver info and context
 *   base_name     <-- base plot name to activate, NULL otherwise
 *   use_iteration <-- if true, use iteration as time stamp
 *                     otherwise, use wall clock time
 *----------------------------------------------------------------------------*/

void
cs_sles_it_set_plot_options(cs_sles_it_t  *context,
                            const char    *base_name,
                            bool           use_iteration);

/*----------------------------------------------------------------------------
 * Convergence test.
 *
 * parameters:
 *   c           <-- pointer to solver context info
 *   n_iter      <-- Number of iterations done
 *   residual    <-- Non normalized residual
 *   convergence <-> Convergence information structure
 *
 * returns:
 *   convergence status.
 *----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_sles_it_convergence_test(cs_sles_it_t              *c,
                            unsigned                   n_iter,
                            double                     residual,
                            cs_sles_it_convergence_t  *convergence);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SLES_IT_H__ */
