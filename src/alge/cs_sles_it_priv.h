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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_blas.h"
#include "cs_file.h"
#include "cs_log.h"
#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_matrix.h"
#include "cs_matrix_default.h"
#include "cs_matrix_util.h"
#include "cs_post.h"
#include "cs_timer.h"
#include "cs_time_plot.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_sles.h"
#include "cs_sles_it.h"
#include "cs_sles_pc.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

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
 *   vx            --> system solution
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
                      cs_real_t                 *restrict vx,
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

  /* Solver setup */

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
 * Compute dot product, summing result over all ranks.
 *
 * parameters:
 *   c      <-- pointer to solver context info
 *   x      <-- first vector in s = x.y
 *   y      <-- second vector in s = x.y
 *
 * returns:
 *   result of s = x.y
 *----------------------------------------------------------------------------*/

inline static double
_dot_product(const cs_sles_it_t  *c,
             const cs_real_t     *x,
             const cs_real_t     *y)
{
  double s = cs_dot(c->setup_data->n_rows, x, y);

#if defined(HAVE_MPI)

  if (c->comm != MPI_COMM_NULL) {
    double _sum;
    MPI_Allreduce(&s, &_sum, 1, MPI_DOUBLE, MPI_SUM, c->comm);
    s = _sum;
  }

#endif /* defined(HAVE_MPI) */

  return s;
}

/*----------------------------------------------------------------------------
 * Compute dot product x.x, summing result over all ranks.
 *
 * parameters:
 *   c      <-- pointer to solver context info
 *   x      <-- vector in s = x.x
 *
 * returns:
 *   result of s = x.x
 *----------------------------------------------------------------------------*/

inline static double
_dot_product_xx(const cs_sles_it_t  *c,
                const cs_real_t     *x)
{
  double s;

  s = cs_dot_xx(c->setup_data->n_rows, x);

#if defined(HAVE_MPI)

  if (c->comm != MPI_COMM_NULL) {
    double _sum;
    MPI_Allreduce(&s, &_sum, 1, MPI_DOUBLE, MPI_SUM, c->comm);
    s = _sum;
  }

#endif /* defined(HAVE_MPI) */

  return s;
}

/*----------------------------------------------------------------------------
 * Compute 2 dot products x.x and x.y, summing result over all ranks.
 *
 * parameters:
 *   c      <-- pointer to solver context info
 *   x      <-- vector in s1 = x.x and s2 = x.y
 *   y      <-- vector in s2 = x.y
 *   s1     --> result of s1 = x.x
 *   s2     --> result of s2 = x.y
 *----------------------------------------------------------------------------*/

inline static void
_dot_products_xx_xy(const cs_sles_it_t  *c,
                    const cs_real_t     *x,
                    const cs_real_t     *y,
                    double              *s1,
                    double              *s2)
{
  double s[2];

  cs_dot_xx_xy(c->setup_data->n_rows, x, y, s, s+1);

#if defined(HAVE_MPI)

  if (c->comm != MPI_COMM_NULL) {
    double _sum[2];
    MPI_Allreduce(s, _sum, 2, MPI_DOUBLE, MPI_SUM, c->comm);
    s[0] = _sum[0];
    s[1] = _sum[1];
  }

#endif /* defined(HAVE_MPI) */

  *s1 = s[0];
  *s2 = s[1];
}

/*----------------------------------------------------------------------------
 * Compute 2 dot products x.x and x.y, summing result over all ranks.
 *
 * parameters:
 *   c      <-- pointer to solver context info
 *   x      <-- vector in s1 = x.y
 *   y      <-- vector in s1 = x.y and s2 = y.z
 *   z      <-- vector in s2 = y.z
 *   s1     --> result of s1 = x.y
 *   s2     --> result of s2 = y.z
 *----------------------------------------------------------------------------*/

inline static void
_dot_products_xy_yz(const cs_sles_it_t  *c,
                    const cs_real_t     *x,
                    const cs_real_t     *y,
                    const cs_real_t     *z,
                    double              *s1,
                    double              *s2)
{
  double s[2];

  cs_dot_xy_yz(c->setup_data->n_rows, x, y, z, s, s+1);

#if defined(HAVE_MPI)

  if (c->comm != MPI_COMM_NULL) {
    double _sum[2];
    MPI_Allreduce(s, _sum, 2, MPI_DOUBLE, MPI_SUM, c->comm);
    s[0] = _sum[0];
    s[1] = _sum[1];
  }

#endif /* defined(HAVE_MPI) */

  *s1 = s[0];
  *s2 = s[1];
}

/*----------------------------------------------------------------------------
 * Compute 3 dot products, summing result over all ranks.
 *
 * parameters:
 *   c      <-- pointer to solver context info
 *   x      <-- first vector
 *   y      <-- second vector
 *   z      <-- third vector
 *   s1     --> result of s1 = x.x
 *   s2     --> result of s2 = x.y
 *   s3     --> result of s3 = y.z
 *----------------------------------------------------------------------------*/

inline static void
_dot_products_xx_xy_yz(const cs_sles_it_t  *c,
                       const cs_real_t     *x,
                       const cs_real_t     *y,
                       const cs_real_t     *z,
                       double              *s1,
                       double              *s2,
                       double              *s3)
{
  double s[3];

  cs_dot_xx_xy_yz(c->setup_data->n_rows, x, y, z, s, s+1, s+2);

#if defined(HAVE_MPI)

  if (c->comm != MPI_COMM_NULL) {
    double _sum[3];

    MPI_Allreduce(s, _sum, 3, MPI_DOUBLE, MPI_SUM, c->comm);
    s[0] = _sum[0];
    s[1] = _sum[1];
    s[2] = _sum[2];
  }

#endif /* defined(HAVE_MPI) */

  *s1 = s[0];
  *s2 = s[1];
  *s3 = s[2];
}

/*----------------------------------------------------------------------------
 * Compute 5 dot products, summing result over all ranks.
 *
 * parameters:
 *   c      <-- pointer to solver context info
 *   x      <-- first vector
 *   y      <-- second vector
 *   z      <-- third vector
 *   xx     --> result of x.x
 *   yy     --> result of y.y
 *   xy     --> result of x.y
 *   xz     --> result of x.z
 *   yz     --> result of y.z
 *----------------------------------------------------------------------------*/

inline static void
_dot_products_xx_yy_xy_xz_yz(const cs_sles_it_t  *c,
                             const cs_real_t     *x,
                             const cs_real_t     *y,
                             const cs_real_t     *z,
                             double              *xx,
                             double              *yy,
                             double              *xy,
                             double              *xz,
                             double              *yz)
{
  double s[5];

  cs_dot_xx_yy_xy_xz_yz(c->setup_data->n_rows, x, y, z, s, s+1, s+2, s+3, s+4);

#if defined(HAVE_MPI)

  if (c->comm != MPI_COMM_NULL) {
    double _sum[5];
    MPI_Allreduce(s, _sum, 5, MPI_DOUBLE, MPI_SUM, c->comm);
    memcpy(s, _sum, 5*sizeof(double));
  }

#endif /* defined(HAVE_MPI) */

  *xx = s[0];
  *yy = s[1];
  *xy = s[2];
  *xz = s[3];
  *yz = s[4];
}

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
_fw_and_bw_lu33(const cs_real_t  mat[],
                cs_real_t        x[restrict],
                const cs_real_t  b[restrict],
                const cs_real_t  c[restrict])
{
  cs_real_t  aux[3];

  aux[0] = (c[0] - b[0]);
  aux[1] = (c[1] - b[1]) - aux[0]*mat[3];
  aux[2] = (c[2] - b[2]) - aux[0]*mat[6] - aux[1]*mat[7];

  x[2] = aux[2]/mat[8];
  x[1] = (aux[1] - mat[5]*x[2])/mat[4];
  x[0] = (aux[0] - mat[1]*x[1] - mat[2]*x[2])/mat[0];
}

/*----------------------------------------------------------------------------
 * Block Jacobi utilities.
 * Compute forward and backward to solve an LU P*P system.
 *
 * parameters:
 *   mat     <-- P*P*dim matrix
 *   db_size <-- matrix size
 *   x       --> solution
 *   b       --> 1st part of RHS (c - b)
 *   c       --> 2nd part of RHS (c - b)
 *----------------------------------------------------------------------------*/

inline static void
_fw_and_bw_lu(const cs_real_t  mat[],
              int              db_size,
              cs_real_t        x[restrict],
              const cs_real_t  b[restrict],
              const cs_real_t  c[restrict])
{
  assert(db_size <= DB_SIZE_MAX);
  cs_real_t aux[DB_SIZE_MAX];

  /* forward */
  for (int ii = 0; ii < db_size; ii++) {
    aux[ii] = (c[ii] - b[ii]);
    for (int jj = 0; jj < ii; jj++) {
      aux[ii] -= aux[jj]*mat[ii*db_size + jj];
    }
  }

  /* backward */
  for (int ii = db_size - 1; ii >= 0; ii-=1) {
    x[ii] = aux[ii];
    for (int jj = db_size - 1; jj > ii; jj-=1) {
      x[ii] -= x[jj]*mat[ii*db_size + jj];
    }
    x[ii] /= mat[ii*(db_size + 1)];
  }
}

/*----------------------------------------------------------------------------
 * Block Gauss-Seidel utilities.
 * Compute forward and backward to solve an LU P*P system.
 *
 * parameters:
 *   mat     <-- P*P*dim matrix
 *   db_size <-- matrix size
 *   x       --> solution
 *   b       <-> RHS in, work array
 *----------------------------------------------------------------------------*/

inline static void
_fw_and_bw_lu_gs(const cs_real_t  mat[],
                 int              db_size,
                 cs_real_t        x[restrict],
                 const cs_real_t  b[restrict])
{
  assert(db_size <= DB_SIZE_MAX);

  /* forward */
  for (int ii = 0; ii < db_size; ii++) {
    x[ii] = b[ii];
    for (int jj = 0; jj < ii; jj++)
      x[ii] -= x[jj]*mat[ii*db_size + jj];
  }

  /* backward */
  for (int ii = db_size - 1; ii >= 0; ii--) {
    for (int jj = db_size - 1; jj > ii; jj--)
      x[ii] -= x[jj]*mat[ii*db_size + jj];
    x[ii] /= mat[ii*(db_size + 1)];
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
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
 *----------------------------------------------------------------------------*/

void
cs_sles_it_setup_priv(cs_sles_it_t       *c,
                      const char         *name,
                      const cs_matrix_t  *a,
                      int                 verbosity,
                      int                 diag_block_size,
                      bool                block_nn_inverse);

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SLES_IT_PRIV_H__ */
