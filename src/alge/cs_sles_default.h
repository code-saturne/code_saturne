#ifndef __CS_SLES_DEFAULT_H__
#define __CS_SLES_DEFAULT_H__

/*============================================================================
 * Sparse Linear Equation Solvers
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_base.h"
#include "base/cs_halo_perio.h"
#include "alge/cs_matrix.h"
#include "alge/cs_sles.h"

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
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Default initializations for sparse linear equation solver API.
 *----------------------------------------------------------------------------*/

void
cs_sles_default_log_setup(void);

/*----------------------------------------------------------------------------
 * Default definition of a sparse linear equation solver
 *
 * parameters:
 *   f_id <-- associated field id, or < 0
 *   name <-- associated name if f_id < 0, or NULL
 *   a    <-- matrix
 *----------------------------------------------------------------------------*/

void
cs_sles_default(int                 f_id,
                const char         *name,
                const cs_matrix_t  *a);

/*----------------------------------------------------------------------------
 * Default setup setup for sparse linear equation solver API.
 *
 * This includes setup logging.
 *----------------------------------------------------------------------------*/

void
cs_sles_default_setup(void);

/*----------------------------------------------------------------------------
 * Return default verbosity associated to a field id, name couple.
 *
 * parameters:
 *   f_id <-- associated field id, or < 0
 *   name <-- associated name if f_id < 0, or NULL
 *
 * returns:
 *   verbosity associated with field or name
 *----------------------------------------------------------------------------*/

int
cs_sles_default_get_verbosity(int          f_id,
                              const char  *name);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return pointer to matrix structure matching equation solve.
 *
 * Some matrix properties (such as assembly options and geometric
 * association) are set immediately, but coefficients are not
 * assigned at this stage.
 *
 * \param[in]  f_id       associated field id, or < 0
 * \param[in]  name       associated name if f_id < 0, or nullptr
 * \param[in]  db_size    block sizes for diagonal
 * \param[in]  eb_size    block sizes for extra diagonal
 * \param[in]  symmetric  indicate if matrix is symmetric
 *
 * \return  pointer to matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_matrix_t *
cs_sles_default_get_matrix(int          f_id,
                           const char  *name,
                           cs_lnum_t    db_size,
                           cs_lnum_t    eb_size,
                           bool         symmetric);

/*----------------------------------------------------------------------------
 * Release of destroy matrix depending on whether is is cached or not.
 *
 * Matrices built by assembler are destroyed.
 *
 * parameters:
 *   matrix <-> pointer to matrix structure
 *----------------------------------------------------------------------------*/

void
cs_sles_default_release_matrix(cs_matrix_t  **m);

/*----------------------------------------------------------------------------
 * Default finalization for sparse linear equation solver API.
 *
 * This includes performance data logging output.
 *----------------------------------------------------------------------------*/

void
cs_sles_default_finalize(void);

/*----------------------------------------------------------------------------
 * Call sparse linear equation solver setup for convection-diffusion
 * systems
 *
 * parameters:
 *   f_id                   associated field id, or < 0
 *   name                   associated name if f_id < 0, or NULL
 *   diag_block_size        block sizes for diagonal
 *   extra_diag_block_size  block sizes for extra diagonal
 *   da                     diagonal values (NULL if zero)
 *   xa                     extradiagonal values (NULL if zero)
 *   conv_diff              convection-diffusion mode
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_setup_native_conv_diff(int                  f_id,
                               const char          *name,
                               const cs_lnum_t      diag_block_size,
                               const cs_lnum_t      extra_diag_block_size,
                               const cs_real_t     *da,
                               const cs_real_t     *xa,
                               bool                 conv_diff);

/*----------------------------------------------------------------------------*/
/*
 * \brief Call sparse linear equation solver for general colocated
 *        cell-centered finite volume scheme.
 *
 * The initial solution is assumed to be 0 (and does not need to
 * be initialized before calling this function).
 *
 * \param[in]       sc                     solver context
 * \param[in]       a                      matrix
 * \param[in]       precision              solver precision
 * \param[in]       r_norm                 residual normalization
 * \param[out]      n_iter                 number of "equivalent" iterations
 * \param[out]      residual               residual
 * \param[in]       rhs                    right hand side
 * \param[out]      vx                     system solution
 *
 * \return  convergence state
 */
/*----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_sles_solve_ccc_fv(cs_sles_t           *sc,
                     cs_matrix_t         *a,
                     double               precision,
                     double               r_norm,
                     int                 *n_iter,
                     double              *residual,
                     const cs_real_t     *rhs,
                     cs_real_t           *vx);

/*----------------------------------------------------------------------------
 * Call sparse linear equation solver using native matrix arrays.
 *
 * parameters:
 *   f_id                   <-- associated field id, or < 0
 *   name                   <-- associated name if f_id < 0, or NULL
 *   symmetric              <-- indicates if matrix coefficients are symmetric
 *   diag_block_size        <-- block sizes for diagonal
 *   extra_diag_block_size  <-- block sizes for extra diagonal
 *   da                     <-- diagonal values (NULL if zero)
 *   xa                     <-- extradiagonal values (NULL if zero)
 *   r_epsilon              <-- precision
 *   r_norm                 <-- residual normalization
 *   n_iter                 --> number of iterations
 *   residual               --> residual
 *   rhs                    <-- right hand side
 *   vx                     <-> system solution
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_sles_solve_native(int                  f_id,
                     const char          *name,
                     bool                 symmetric,
                     cs_lnum_t            diag_block_size,
                     cs_lnum_t            extra_diag_block_size,
                     const cs_real_t     *da,
                     const cs_real_t     *xa,
                     double               precision,
                     double               r_norm,
                     int                 *n_iter,
                     double              *residual,
                     const cs_real_t     *rhs,
                     cs_real_t           *vx);

/*----------------------------------------------------------------------------
 * Free sparse linear equation solver setup using native matrix arrays.
 *
 * parameters:
 *   f_id                   <-- associated field id, or < 0
 *   name                   <-- associated name if f_id < 0, or NULL
 *----------------------------------------------------------------------------*/

void
cs_sles_free_native(int          f_id,
                    const char  *name);

/*----------------------------------------------------------------------------
 * Error handler attempting fallback to alternative solution procedure for
 * sparse linear equation solver.
 *
 * In case of divergence with an iterative solver, this error handler
 * switches to a default preconditioner, then resets the solution vector.
 *
 * The default error for the solver type handler is then  set, in case
 * the solution fails again.
 *
 * Note that this error handler may rebuild solver contexts, so should not
 * be used in conjunction with shared contexts (such as multigrid
 * ascent/descent contexts), but only for "outer" solvers.
 *
 * parameters:
 *   sles          <-> pointer to solver object
 *   state         <-- convergence status
 *   a             <-- matrix
 *   rhs           <-- right hand side
 *   vx            <-> system solution
 *
 * returns:
 *   true if fallback solution is possible, false otherwise
 *----------------------------------------------------------------------------*/

bool
cs_sles_default_error(cs_sles_t                    *sles,
                      cs_sles_convergence_state_t   state,
                      const cs_matrix_t            *a,
                      const cs_real_t               rhs[],
                      cs_real_t                     vx[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SLES_DEFAULT_H__ */
