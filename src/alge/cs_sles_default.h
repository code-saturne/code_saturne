#ifndef __CS_SLES_DEFAULT_H__
#define __CS_SLES_DEFAULT_H__

/*============================================================================
 * Sparse Linear Equation Solvers
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_halo_perio.h"
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

/*----------------------------------------------------------------------------
 * Default finalization for sparse linear equation solver API.
 *
 * This includes performance data logging output.
 *----------------------------------------------------------------------------*/

void
cs_sles_default_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Call sparse linear equation solver setup for systems with internal
 *        coupling.
 *
 * \param[in]  f_id                   associated field id, or < 0
 * \param[in]  name                   associated name if f_id < 0, or NULL
 * \param[in]  symmetric              indicates if matrix coefficients
 *                                    are symmetric
 * \param[in]  diag_block_size        block sizes for diagonal, or NULL
 * \param[in]  extra_diag_block_size  block sizes for extra diagonal,
 *                                    or NULL
 * \param[in]  da                     diagonal values (NULL if zero)
 * \param[in]  xa                     extradiagonal values (NULL if zero)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_setup_native_coupling(int               f_id,
                              const char       *name,
                              bool              symmetric,
                              const int        *diag_block_size,
                              const int        *extra_diag_block_size,
                              const cs_real_t  *da,
                              const cs_real_t  *xa);

/*----------------------------------------------------------------------------
 * Call sparse linear equation solver setup for convection-diffusion
 * systems
 *
 * parameters:
 *   f_id                   associated field id, or < 0
 *   name                   associated name if f_id < 0, or NULL
 *   diag_block_size        block sizes for diagonal, or NULL
 *   extra_diag_block_size  block sizes for extra diagonal, or NULL
 *   da                     diagonal values (NULL if zero)
 *   xa                     extradiagonal values (NULL if zero)
 *   da_conv                diagonal values (NULL if zero)
 *   xa_conv                extradiagonal values (NULL if zero)
 *   da_diff                diagonal values (NULL if zero)
 *   xa_diff                extradiagonal values (NULL if zero)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_setup_native_conv_diff(int                  f_id,
                               const char          *name,
                               const int           *diag_block_size,
                               const int           *extra_diag_block_size,
                               const cs_real_t     *da,
                               const cs_real_t     *xa,
                               const cs_real_t     *da_conv,
                               const cs_real_t     *xa_conv,
                               const cs_real_t     *da_diff,
                               const cs_real_t     *xa_diff);

/*----------------------------------------------------------------------------
 * Call sparse linear equation solver using native matrix arrays.
 *
 * parameters:
 *   f_id                   <-- associated field id, or < 0
 *   name                   <-- associated name if f_id < 0, or NULL
 *   symmetric              <-- indicates if matrix coefficients are symmetric
 *   diag_block_size        <-- block sizes for diagonal, or NULL
 *   extra_diag_block_size  <-- block sizes for extra diagonal, or NULL
 *   da                     <-- diagonal values (NULL if zero)
 *   xa                     <-- extradiagonal values (NULL if zero)
 *   rotation_mode          <-- halo update option for rotational periodicity
 *   r_epsilon              <-- precision
 *   r_norm                 <-- residue normalization
 *   n_iter                 --> number of iterations
 *   residue                --> residue
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
                     const int           *diag_block_size,
                     const int           *extra_diag_block_size,
                     const cs_real_t     *da,
                     const cs_real_t     *xa,
                     cs_halo_rotation_t   rotation_mode,
                     double               precision,
                     double               r_norm,
                     int                 *n_iter,
                     double              *residue,
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
 *   rotation_mode <-- halo update option for rotational periodicity
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
                      cs_halo_rotation_t            rotation_mode,
                      const cs_real_t               rhs[],
                      cs_real_t                     vx[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SLES_DEFAULT_H__ */
