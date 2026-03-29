#ifndef __CS_SLES_IT_PRIV_H__
#define __CS_SLES_IT_PRIV_H__

/*============================================================================
 * Sparse Linear Equation Solvers: private elements.
 *
 * These elements are shared between iterative solvers and smoother
 * both for host and device implementations, but are not accessible to
 * calling code.
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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#if defined(HAVE_NCCL)
#include <nccl.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_base.h"
#include "alge/cs_blas.h"
#include "base/cs_file.h"
#include "base/cs_log.h"
#include "base/cs_halo.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "alge/cs_matrix.h"
#include "alge/cs_matrix_default.h"
#include "alge/cs_matrix_util.h"
#include "base/cs_post.h"
#include "base/cs_timer.h"
#include "base/cs_time_plot.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "alge/cs_sles.h"
#include "alge/cs_sles_it.h"
#include "alge/cs_sles_pc.h"

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#if !defined(HUGE_VAL)
#define HUGE_VAL  1.E+12
#endif

#define DB_SIZE_MAX 9

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Function pointer for actual resolution of a linear system.
 *
 * parameters:
 *   c             <-- pointer to solver context info
 *   a             <-- linear equation matrix
 *   convergence   <-- convergence information structure
 *   rhs           <-- right hand side
 *   vx_ini        <-- initial system solution
 *                     (vx if nonzero, nullptr if zero)
 *   vx            <-> system solution
 *   aux_size      <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors   --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence status
 *----------------------------------------------------------------------------*/

typedef cs_sles_convergence_state_t
(cs_sles_it_solve_t) (cs_sles_it_t              *c,
                      const cs_matrix_t         *a,
                      cs_lnum_t                  diag_block_size,
                      cs_sles_it_convergence_t  *convergence,
                      const cs_real_t           *rhs,
                      cs_real_t                 *vx_ini,
                      cs_real_t                 *vx,
                      size_t                     aux_size,
                      void                      *aux_vectors);

/* Solver setup data */
/*-------------------*/

typedef struct _cs_sles_it_setup_t {

  double               initial_residual; /* last initial residual value */

  cs_lnum_t            n_rows;           /* number of associated rows */

  const cs_real_t     *ad_inv;           /* pointer to diagonal inverse */
  cs_real_t           *_ad_inv;          /* private pointer to
                                            diagonal inverse */

  void                *pc_context;       /* preconditioner context */
  cs_sles_pc_apply_t  *pc_apply;         /* preconditioner apply */

} cs_sles_it_setup_t;

/* Solver additional data */
/*------------------------*/

typedef struct _cs_sles_it_add_t {

  cs_lnum_t           *order;            /* ordering */

} cs_sles_it_add_t;

/* Basic per linear system options and logging */
/*---------------------------------------------*/

struct _cs_sles_it_t {

  /* Base settings */

  cs_sles_it_type_t    type;               /* Solver type */

  bool                 on_device;          /* SpMV on device ? */

  bool                 update_stats;       /* do stats need to be updated ? */
  bool                 ignore_convergence; /* ignore convergence for some
                                              solvers used as preconditioners */

  int                  n_max_iter;         /* maximum number of iterations */
  int                  restart_interval;   /* maximum number of iterations
                                              before restarting the algorithm
                                              (only applicable for GMRES or GCR
                                              algorithm up to now) */

  cs_sles_it_solve_t  *solve;              /* pointer to solve function */

  cs_sles_pc_t        *pc;                 /* pointer to possibly shared
                                              preconditioner object */
  cs_sles_pc_t        *_pc;                /* pointer to owned
                                              preconditioner object */

  /* Performance data */

  unsigned             n_setups;           /* Number of times system setup */
  unsigned             n_solves;           /* Number of times system solved */

  unsigned             n_iterations_last;  /* Number of iterations for last
                                              system resolution */
  unsigned             n_iterations_min;   /* Minimum number ot iterations
                                              in system resolution history */
  unsigned             n_iterations_max;   /* Maximum number ot iterations
                                              in system resolution history */
  unsigned long long   n_iterations_tot;   /* Total accumulated number of
                                              iterations */

  cs_timer_counter_t   t_setup;            /* Total setup */
  cs_timer_counter_t   t_solve;            /* Total time used */

  /* Plot info */

  int                  plot_time_stamp;    /* Plot time stamp */
  cs_time_plot_t      *plot;               /* Pointer to plot structure,
                                              which may be owned or shared */
  cs_time_plot_t      *_plot;              /* Pointer to own plot structure */

  /* Communicator used for reduction operations
     (if left at NULL, main communicator will be used) */

# if defined(HAVE_MPI)
  MPI_Comm comm;
  MPI_Comm caller_comm;
  int      caller_n_ranks;
# endif

#if defined(HAVE_NCCL)
  ncclComm_t nccl_comm;
#endif

  const struct _cs_sles_it_t  *shared;     /* pointer to context sharing some
                                              setup and preconditioner data,
                                              or NULL */

  cs_sles_it_add_t            *add_data;   /* additional data */

  cs_sles_it_setup_t          *setup_data; /* setup data */

  /* Alternative solvers (fallback or heuristics) */

  cs_sles_convergence_state_t  fallback_cvg;  /* threshold for fallback
                                                 convergence */
  int                          fallback_n_max_iter; /* number of maximum iteration
                                                       for fallback solver */

  cs_sles_it_t                *fallback;      /* fallback solver */


};

/* Convergence testing and tracking */
/*----------------------------------*/

struct _cs_sles_it_convergence_t {

  const char          *name;               /* Pointer to name string */

  int                  verbosity;          /* Verbosity level */

  unsigned             n_iterations;       /* Current number of iterations */
  unsigned             n_iterations_max;   /* Maximum number of iterations */

  double               precision;          /* Precision limit */
  double               r_norm;             /* Residual normalization */
  double               residual;           /* Current residual */

};

/*============================================================================
 * Inline static function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Block Jacobi utilities.
 * Compute forward and backward to solve an LU 3*3 system.
 *
 * parameters:
 *   mat   <-- 3*3*dim matrix
 *   x     --> solution
 *   b     --> 1st part of RHS (c - b)
 *   c     --> 2nd part of RHS (c - b)
 *----------------------------------------------------------------------------*/

inline static void
_mat_c_m_b_33(const cs_real_t    mat[],
              cs_real_t        *restrict x,
              const cs_real_t  *restrict b,
              const cs_real_t  *restrict c)
{
  /* c - b */
  cs_real_t  aux[3];
  for (cs_lnum_t ii = 0; ii < 3; ii++) {
    aux[ii] = (c[ii] - b[ii]);
  }

  for (cs_lnum_t ii = 0; ii < 3; ii++) {
    x[ii] =   mat[3*ii + 0]*aux[0]
            + mat[3*ii + 1]*aux[1]
            + mat[3*ii + 2]*aux[2];
  }
}

/*----------------------------------------------------------------------------
 * Block Jacobi utilities.
 * Compute mat.(c-b) product.
 *
 * parameters:
 *   mat     <-- P*P*dim matrix
 *   db_size <-- matrix size
 *   x       --> solution
 *   b       --> 1st part of RHS (c - b)
 *   c       --> 2nd part of RHS (c - b)
 *----------------------------------------------------------------------------*/

inline static void
_mat_c_m_b(const cs_real_t   mat[],
           cs_lnum_t         db_size,
           cs_real_t        *restrict x,
           const cs_real_t  *restrict b,
           const cs_real_t  *restrict c)
{
  assert(db_size <= DB_SIZE_MAX);
  cs_real_t aux[DB_SIZE_MAX];

  /* c - b */
  for (cs_lnum_t ii = 0; ii < db_size; ii++) {
    aux[ii] = (c[ii] - b[ii]);
  }

  for (cs_lnum_t ii = 0; ii < db_size; ii++) {
    x[ii] = 0;
    for (cs_lnum_t jj = 0; jj < db_size; jj++) {
      x[ii] += aux[jj]*mat[ii*db_size + jj];
    }
  }
}

/*============================================================================
 * Semi-private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Set execution location (host or device) and stream if appplicable.
 *
 * \param[in, out]  ctx            reference to dispatch context
 * \param[in]       a              pointer to matrix
 * \param[out]      local_stream   do we force a local stream ?
 * \param[out]      stream         stream associated with current solve
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_set_exec_location
(
  [[maybe_unused]] cs_dispatch_context  &ctx,
  const cs_matrix_t                     *a,
  bool                                  &local_stream,
  [[maybe_unused]] cs_stream_t          &stream
);

/*----------------------------------------------------------------------------*/
/*
 * \brief Restore stream to previous settings if appplicable.
 *
 * \param[in]  local_stream   do we force a local stream ?
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_restore_exec_location([[maybe_unused]] bool  &local_stream);

/*----------------------------------------------------------------------------*/
/*
 * \brief Initialize or reset convergence info structure.
 *        At this stage, the initial residual is set to HUGE_VAL, as it is
 *        unknown.
 *
 * \param[in, out]  convergence   convergence info structure
 * \param[in]       solver_name   solver name
 * \param[in]       verbosity     verbosity level
 * \param[in]       n_iter_max    maximum number of iterations
 * \param[in]       precision     precision limit
 * \param[in]       r_norm        residual normalization
 * \param[in, out]  residual      initial residual
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_convergence_init(cs_sles_it_convergence_t  *convergence,
                            const char                *solver_name,
                            int                        verbosity,
                            unsigned                   n_iter_max,
                            double                     precision,
                            double                     r_norm,
                            double                    *residual);

/*----------------------------------------------------------------------------
 * Setup context for iterative linear solver.
 *
 * This function is common to most solvers
 *
 * parameters:
 *   c                <-> pointer to solver context info
 *   name             <-- pointer to system name
 *   a                <-- matrix
 *   verbosity        <-- verbosity level
 *   diag_block_size  <-- diagonal block size
 *   block_nn_inverse <-- if diagonal block size is 3 or 6, compute inverse of
 *                        block if true, inverse of block diagonal otherwise
 *   l1_inv           <-- if diagonal block size is 1, compute
 *                        inverse of L1 norm instead of diagonal
 *----------------------------------------------------------------------------*/

void
cs_sles_it_setup_priv(cs_sles_it_t       *c,
                      const char         *name,
                      const cs_matrix_t  *a,
                      int                 verbosity,
                      int                 diag_block_size,
                      bool                block_nn_inverse,
                      bool                l1_inv);

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*----------------------------------------------------------------------------*/

#endif /* __CS_SLES_IT_PRIV_H__ */
