#ifndef __CS_MULTIGRID_H__
#define __CS_MULTIGRID_H__

/*============================================================================
 * Multigrid solver.
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
#include "cs_grid.h"
#include "cs_sles.h"
#include "cs_sles_it.h"
#include "cs_sles_pc.h"
#include "cs_time_plot.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Multigrid types
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_MULTIGRID_V_CYCLE,          /*!< Use a V-cycle */
  CS_MULTIGRID_K_CYCLE,          /*!< Use a K-cycle */
  CS_MULTIGRID_K_CYCLE_HPC,      /*!< Use a K-cycle tuned for HPC systems */
  CS_MULTIGRID_N_TYPES           /*!< Number of multigrid types */

} cs_multigrid_type_t;

/* Multigrid linear solver context (opaque) */

typedef struct _cs_multigrid_t  cs_multigrid_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Names for multigrid types */

extern const char *cs_multigrid_type_name[];

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize multigrid solver API.
 *----------------------------------------------------------------------------*/

void
cs_multigrid_initialize(void);

/*----------------------------------------------------------------------------
 * Finalize multigrid solver API.
 *----------------------------------------------------------------------------*/

void
cs_multigrid_finalize(void);

/*----------------------------------------------------------------------------
 * Define and associate a multigrid sparse linear system solver
 * for a given field or equation name.
 *
 * If this system did not previously exist, it is added to the list of
 * "known" systems. Otherwise, its definition is replaced by the one
 * defined here.
 *
 * This is a utility function: if finer control is needed, see
 * cs_sles_define() and cs_multigrid_create().
 *
 * Note that this function returns a pointer directly to the multigrid solver
 * management structure. This may be used to set further options, for
 * example calling cs_multigrid_set_coarsening_options() and
 * cs_multigrid_set_solver_options().
 * If needed, cs_sles_find() may be used to obtain a pointer to the
 * matching cs_sles_t container.
 *
 * \param[in]  f_id     associated field id, or < 0
 * \param[in]  name     associated name if f_id < 0, or NULL
 * \param[in]  mg_type  type of multigrid algorithm to use
 *
 * \return  pointer to new multigrid info and context
 */
/*----------------------------------------------------------------------------*/

cs_multigrid_t *
cs_multigrid_define(int                   f_id,
                    const char           *name,
                    cs_multigrid_type_t   mg_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create multigrid linear system solver info and context.
 *
 * The multigrid variant is an ACM (Additive Corrective Multigrid) method.
 *
 * \param[in]  mg_type  type of multigrid algorithm to use
 *
 * \return  pointer to new multigrid info and context
 */
/*----------------------------------------------------------------------------*/

cs_multigrid_t *
cs_multigrid_create(cs_multigrid_type_t  mg_type);

/*----------------------------------------------------------------------------
 * Destroy multigrid linear system solver info and context.
 *
 * parameters:
 *   context  <-> pointer to multigrid linear solver info
 *                (actual type: cs_multigrid_t  **)
 *----------------------------------------------------------------------------*/

void
cs_multigrid_destroy(void  **context);

/*----------------------------------------------------------------------------
 * Create multigrid sparse linear system solver info and context
 * based on existing info and context.
 *
 * parameters:
 *   context <-- pointer to reference info and context
 *               (actual type: cs_multigrid_t  *)
 *
 * returns:
 *   pointer to newly created solver info object
 *   (actual type: cs_multigrid_t  *)
 *----------------------------------------------------------------------------*/

void *
cs_multigrid_copy(const void  *context);

/*----------------------------------------------------------------------------
 * Set multigrid coarsening parameters.
 *
 * parameters:
 *   mg                     <-> pointer to multigrid info and context
 *   aggregation_limit      <-- maximum allowed fine rows per coarse cell
 *   coarsening_type        <-- coarsening type; see cs_grid_coarsening_t
 *   n_max_levels           <-- maximum number of grid levels
 *   min_g_rows             <-- global number of rows on coarse grids
 *                              under which no coarsening occurs
 *   p0p1_relax             <-- p0/p1 relaxation_parameter
 *   postprocess_block_size <-- if > 0, postprocess coarsening
 *                              (using coarse cell numbers modulo this value)
 *----------------------------------------------------------------------------*/

void
cs_multigrid_set_coarsening_options(cs_multigrid_t  *mg,
                                    int              aggregation_limit,
                                    int              coarsening_type,
                                    int              n_max_levels,
                                    cs_gnum_t        min_g_rows,
                                    double           p0p1_relax,
                                    int              postprocess_block_size);

/*----------------------------------------------------------------------------
 * Set multigrid parameters for associated iterative solvers.
 *
 * parameters:
 *   mg                     <-> pointer to multigrid info and context
 *   descent_smoother_type  <-- type of smoother for descent
 *   ascent_smoother_type   <-- type of smoother for ascent
 *   coarse_solver_type     <-- type of solver
 *   n_max_cycles           <-- maximum number of cycles
 *   n_max_iter_descent     <-- maximum iterations per descent phase
 *   n_max_iter_ascent      <-- maximum iterations per descent phase
 *   n_max_iter_coarse      <-- maximum iterations per coarsest solution
 *   poly_degree_descent    <-- preconditioning polynomial degree
 *                              for descent phases (0: diagonal)
 *   poly_degree_ascent     <-- preconditioning polynomial degree
 *                              for ascent phases (0: diagonal)
 *   poly_degree_coarse     <-- preconditioning polynomial degree
 *                              for coarse solver  (0: diagonal)
 *   precision_mult_descent <-- precision multiplier for descent phases
 *                              (levels >= 1)
 *   precision_mult_ascent  <-- precision multiplier for ascent phases
 *   precision_mult_coarse  <-- precision multiplier for coarsest grid
 *----------------------------------------------------------------------------*/

void
cs_multigrid_set_solver_options(cs_multigrid_t     *mg,
                                cs_sles_it_type_t   descent_smoother_type,
                                cs_sles_it_type_t   ascent_smoother_type,
                                cs_sles_it_type_t   coarse_solver_type,
                                int                 n_max_cycles,
                                int                 n_max_iter_descent,
                                int                 n_max_iter_ascent,
                                int                 n_max_iter_coarse,
                                int                 poly_degree_descent,
                                int                 poly_degree_ascent,
                                int                 poly_degree_coarse,
                                double              precision_mult_descent,
                                double              precision_mult_ascent,
                                double              precision_mult_coarse);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the max. number of cycles for a multigrid
 *
 * \param[in, out]  mg              pointer to multigrid info and context
 * \param[in]       n_max_cycles    maximum number of cycles
 */
/*----------------------------------------------------------------------------*/

void
cs_multigrid_set_max_cycles(cs_multigrid_t     *mg,
                            int                 n_max_cycles);

/*----------------------------------------------------------------------------
 * Return solver type used on fine mesh.
 *
 * parameters:
 *   mg <-- pointer to multigrid info and context
 *
 * returns:
 *   type of smoother for descent (used for fine mesh)
 *----------------------------------------------------------------------------*/

cs_sles_it_type_t
cs_multigrid_get_fine_solver_type(const cs_multigrid_t  *mg);

/*----------------------------------------------------------------------------
 * Setup multigrid sparse linear equation solver.
 *
 * parameters:
 *   context   <-> pointer to multigrid info and context
 *                 (actual type: cs_multigrid_t  *)
 *   name      <-- pointer to name of linear system
 *   a         <-- associated matrix
 *   verbosity <-- associated verbosity
 *----------------------------------------------------------------------------*/

void
cs_multigrid_setup(void               *context,
                   const char         *name,
                   const cs_matrix_t  *a,
                   int                 verbosity);

/*----------------------------------------------------------------------------
 * Setup multigrid sparse linear equation solver with separate
 * convection-diffusion matrixes
 *
 * parameters:
 *   context   <-> pointer to multigrid info and context
 *                 (actual type: cs_multigrid_t  *)
 *   name      <-- pointer to name of linear system
 *   a         <-- associated matrix
 *   conv_diff <-- convection-diffusion mode
 *   verbosity <-- associated verbosity
 *----------------------------------------------------------------------------*/

void
cs_multigrid_setup_conv_diff(void               *context,
                             const char         *name,
                             const cs_matrix_t  *a,
                             bool                conv_diff,
                             int                 verbosity);

/*----------------------------------------------------------------------------
 * Call multigrid sparse linear equation solver.
 *
 * parameters:
 *   context       <-> pointer to iterative sparse linear solver info
 *                     (actual type: cs_multigrid_t  *)
 *   name          <-- pointer to name of linear system
 *   a             <-- matrix
 *   verbosity     <-- associated verbosity
 *   precision     <-- solver precision
 *   r_norm        <-- residual normalization
 *   n_iter        --> number of iterations
 *   residual      --> residual
 *   rhs           <-- right hand side
 *   vx            <-> system solution
 *   aux_size      <-- number of elements in aux_vectors
 *   aux_vectors   --- optional working area (internal allocation if NULL)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_multigrid_solve(void                *context,
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
 * Note that this function should free resolution-related data, such as
 * buffers and preconditioning but does not free the whole context,
 * as info used for logging (especially performance data) is maintained.

 * parameters:
 *   context <-> pointer to iterative sparse linear solver info
 *               (actual type: cs_multigrid_t  *)
 *----------------------------------------------------------------------------*/

void
cs_multigrid_free(void  *context);

/*----------------------------------------------------------------------------
 * Log sparse linear equation solver info.
 *
 * parameters:
 *   context  <-> pointer to iterative sparse linear solver info
 *                (actual type: cs_multigrid_t  *)
 *   log_type <-- log type
 *----------------------------------------------------------------------------*/

void
cs_multigrid_log(const void  *context,
                 cs_log_t     log_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a multigrid preconditioner.
 *
 * \param[in]  mg_type  type of multigrid algorithm to use
 *
 * \return  pointer to newly created preconditioner object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_pc_t *
cs_multigrid_pc_create(cs_multigrid_type_t  mg_type);

/*----------------------------------------------------------------------------
 * Error handler for multigrid sparse linear equation solver.
 *
 * In case of divergence or breakdown, this error handler outputs
 * postprocessing data to assist debugging, then aborts the run.
 * It does nothing in case the maximum iteration count is reached.
 *
 * parameters:
 *   sles          <-> pointer to solver object
 *   state         <-- convergence status
 *   a             <-- matrix
 *   rhs           <-- right hand side
 *   vx            <-> system solution
 *
 * returns:
 *   false (do not attempt new solve)
 *----------------------------------------------------------------------------*/

bool
cs_multigrid_error_post_and_abort(cs_sles_t                    *sles,
                                  cs_sles_convergence_state_t   state,
                                  const cs_matrix_t            *a,
                                  const cs_real_t              *rhs,
                                  cs_real_t                    *vx);

/*----------------------------------------------------------------------------
 * Set plotting options for multigrid.
 *
 * parameters:
 *   mg            <-> pointer to multigrid info and context
 *   base_name     <-- base plot name to activate, NULL otherwise
 *   use_iteration <-- if true, use iteration as time stamp
 *                     otherwise, use wall clock time
 *----------------------------------------------------------------------------*/

void
cs_multigrid_set_plot_options(cs_multigrid_t  *mg,
                              const char      *base_name,
                              bool             use_iteration);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Query the global multigrid parameters for parallel grid merging.
 *
 * \param[in]   mg                   pointer to multigrid info and context
 * \param[out]  rank_stride          number of ranks over which merging
 *                                   takes place, or NULL
 * \param[out]  rows_mean_threshold  mean number of rows under which merging
 *                                   should be applied, or NULL
 * \param[out]  rows_glob_threshold  global number of rows under which
 *                                   merging should be applied, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_multigrid_get_merge_options(const cs_multigrid_t  *mg,
                               int                   *rank_stride,
                               int                   *rows_mean_threshold,
                               cs_gnum_t             *rows_glob_threshold);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set global multigrid parameters for parallel grid merging behavior.
 *
 * \param[in, out]  mg                   pointer to multigrid info and context
 * \param[in]       rank_stride          number of ranks over which merging
 *                                       takes place
 * \param[in]       rows_mean_threshold  mean number of rows under which
 *                                       merging should be applied
 * \param[in]       rows_glob_threshold  global number of rows under which
 *                                       merging should be applied
 */
/*----------------------------------------------------------------------------*/

void
cs_multigrid_set_merge_options(cs_multigrid_t  *mg,
                               int              rank_stride,
                               int              rows_mean_threshold,
                               cs_gnum_t        rows_glob_threshold);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a grid associated with a given multigrid
 *        setup and level.
 *
 * If the multigrid hierarchy is not set up, or a level coarser than the
 * coarsest level is requested, NULL is returned.

 * \param[in]  mg     pointer to multigrid info and context
 * \param[in]  level  level of the requested grid (or -1 for coarsest)
 *
 * \return  pointer to grid of requested level (NULL id not present)
 */
/*----------------------------------------------------------------------------*/

const cs_grid_t *
cs_multigrid_get_grid(const cs_multigrid_t  *mg,
                      int                    level);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MULTIGRID_H__ */
