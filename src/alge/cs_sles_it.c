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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_sles_it.h"
#include "cs_sles_it_priv.h"

#if defined(HAVE_CUDA)
#include "cs_sles_it_cuda.h"
#endif

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_sles_it.c
        Iterative linear solvers

 \page sles_it Iterative linear solvers.

 For Krylov space solvers, default preconditioning is based
 on a Neumann polynomial of degree \a poly_degree, with a negative value
 meaning no preconditioning, and 0 diagonal preconditioning.

 For positive values of \a poly_degree, the preconditioning is explained here:
 \a D being the diagonal part of matrix \a A and \a X its extra-diagonal
 part, it can be written \f$A=D(Id+D^{-1}X)\f$. Therefore
 \f$A^{-1}=(Id+D^{-1}X)^{-1}D^{-1}\f$. A series development of
 \f$Id+D^{-1}X\f$ can then be used which yields, symbolically,
 \f[
 Id+\sum\limits_{I=1}^{poly\_degree}\left(-D^{-1}X\right)^{I}
 \f]

 The efficiency of the polynomial preconditioning will vary depending
 on the system type. In most cases, diagonal or degree 1 provide
 best results. Each polynomial preconditioning degree above 0 adds one
 matrix-vector product per inital matrix-vector product of the algorithm.
 Switching from diagonal to polynomial degree 1 often divides the number of
 required iterations by approximately 2, but each iteration then costs
 close to 2 times that of diagonal preconditioning (other vector operations
 and reductions are not doubled), so the net gain is often about 10%.
 Higher degree polynomials usually lead to diminishing returns.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define RINFIN  1.E+30

/* SIMD unit size to ensure SIMD alignement (2 to 8 required on most
 * current architectures, so 16 should be enough on most architectures) */

#define CS_SIMD_SIZE(s) (((s-1)/16+1)*16)

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

static int _thread_debug = 0;

/* Mean system rows threshold under which we use single-reduce version of PCG */

static cs_lnum_t _pcg_sr_threshold = 512;

/* Value of the threshold under which BiCGStab and BiCGStab2 break down */

static double  _epzero = 1.e-30; /* smaller than epzero */

/* Sparse linear equation solver type names */

const char *cs_sles_it_type_name[]
  = {N_("Conjugate Gradient"),
     N_("Flexible Conjugate Gradient"),
     N_("Inexact Preconditioned Conjugate Gradient"),
     N_("Jacobi"),
     N_("BiCGstab"),
     N_("BiCGstab2"),
     N_("GCR"),
     N_("GMRES"),
     N_("Gauss-Seidel"),
     N_("Symmetric Gauss-Seidel"),
     N_("3-layer conjugate residual"),
     N_("User-defined iterative solver"),
     N_("None"), /* Smoothers beyond this */
     N_("Truncated forward Gauss-Seidel"),
     N_("Truncated backwards Gauss-Seidel"),
};

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return a solver category based on possible preconditioning type.
 *
 * parameters:
 *   solver_type  <-- solver type
 *
 * returns:
 *   1 for Krylov-type solvers, 2 for Jacobi or Gauss-seidel, 0 otherwise
 *----------------------------------------------------------------------------*/

static int
_setup_category(cs_sles_it_type_t  solver_type)
{
  int retval = 1;

  if (   solver_type == CS_SLES_JACOBI
      || solver_type == CS_SLES_P_GAUSS_SEIDEL
      || solver_type == CS_SLES_P_SYM_GAUSS_SEIDEL)
    retval = 2;

  else if (solver_type >= CS_SLES_USER_DEFINED)
    retval = 0;

  return retval;
}

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

static cs_sles_convergence_state_t
_convergence_test(cs_sles_it_t              *c,
                  unsigned                   n_iter,
                  double                     residual,
                  cs_sles_it_convergence_t  *convergence)
{
  const int verbosity = convergence->verbosity;
  const cs_sles_it_setup_t  *s = c->setup_data;

  const char final_fmt[]
    = N_("  n_iter: %5d, res_abs: %11.4e, res_nor: %11.4e, norm: %11.4e,"
         " res_init: %11.4e\n");

  /* Update conversion info structure */

  convergence->n_iterations = n_iter;
  convergence->residual = residual;

  /* Plot convergence if requested */

  if (c->plot != NULL) {
    double vals = residual;
    double wall_time = cs_timer_wtime();
    c->plot_time_stamp += 1;
    cs_time_plot_vals_write(c->plot,             /* time plot structure */
                            c->plot_time_stamp,  /* current iteration number */
                            wall_time,           /* current time */
                            1,                   /* number of values */
                            &vals);              /* values */

  }

  /* If converged */

  if (residual < convergence->precision * convergence->r_norm) {
    if (verbosity > 1)
      bft_printf(_(final_fmt), n_iter, residual, residual/convergence->r_norm,
                 convergence->r_norm, s->initial_residual);
    return CS_SLES_CONVERGED;
  }

  /* If not converged */
  else if (n_iter >= convergence->n_iterations_max) {
    bool final_iteration = (c->fallback_cvg <= CS_SLES_MAX_ITERATION);
    if (verbosity > -1) {
      if (verbosity <= 1) { /* Already output if verbosity > 1 */
        if (final_iteration)
          bft_printf("%s [%s]:\n", cs_sles_it_type_name[c->type],
                     convergence->name);
      }
      if (verbosity > 1 || final_iteration) {
        if (convergence->r_norm > 0.)
          bft_printf(_(final_fmt),
                     n_iter, residual, residual/convergence->r_norm,
                     convergence->r_norm, s->initial_residual);
        else
          bft_printf(_("  n_iter : %5d, res_abs : %11.4e\n"),
                     n_iter, residual);
      }
      if (convergence->precision > 0. && final_iteration)
        bft_printf(_(" @@ Warning: non convergence\n"));
    }
    return CS_SLES_MAX_ITERATION;
  }

  /* If diverged */
  else {
    int diverges = 0;
    if (residual > s->initial_residual * 10000.0 && residual > 100.)
      diverges = 1;
#if (__STDC_VERSION__ >= 199901L)
    else if (isnan(residual) || isinf(residual)) {
      diverges = 1;
    }
#endif
    if (diverges == 1) {
      bft_printf(_("\n\n"
                   "%s [%s]: divergence after %u iterations:\n"
                   "  initial residual: %11.4e; current residual: %11.4e\n"),
                 cs_sles_it_type_name[c->type], convergence->name,
                 convergence->n_iterations,
                 s->initial_residual, convergence->residual);
      return CS_SLES_DIVERGED;
    }
  }

  return CS_SLES_ITERATING;
}

/*----------------------------------------------------------------------------
 * Compute 4 dot products, summing result over all ranks.
 *
 * parameters:
 *   c      <-- pointer to solver context info
 *   v      <-- first vector
 *   r      <-- second vector
 *   w      <-- third vector
 *   q      <-- fourth vector
 *   s1     --> result of s1 = v.r
 *   s2     --> result of s2 = v.w
 *   s3     --> result of s3 = v.q
 *   s4     --> result of s4 = r.r
 *----------------------------------------------------------------------------*/

inline static void
_dot_products_vr_vw_vq_rr(const cs_sles_it_t  *c,
                          const cs_real_t     *v,
                          const cs_real_t     *r,
                          const cs_real_t     *w,
                          const cs_real_t     *q,
                          double              *s1,
                          double              *s2,
                          double              *s3,
                          double              *s4)
{
  double s[4];

  /* Use two separate call as cs_blas.c does not yet hav matching call */

  cs_dot_xy_yz(c->setup_data->n_rows, w, v, q, s+1, s+2);
  cs_dot_xx_xy(c->setup_data->n_rows, r, v, s+3, s);

#if defined(HAVE_MPI)

  if (c->comm != MPI_COMM_NULL) {
    double _sum[4];
    MPI_Allreduce(s, _sum, 4, MPI_DOUBLE, MPI_SUM, c->comm);
    memcpy(s, _sum, 4*sizeof(double));
  }

#endif /* defined(HAVE_MPI) */

  *s1 = s[0];
  *s2 = s[1];
  *s3 = s[2];
  *s4 = s[3];
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using preconditioned conjugate gradient.
 *
 * Parallel-optimized version, groups dot products, at the cost of
 * computation of the preconditionning for n+1 iterations instead of n.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- matrix
 *   diag_block_size <-- diagonal block size
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_conjugate_gradient(cs_sles_it_t              *c,
                    const cs_matrix_t         *a,
                    cs_lnum_t                  diag_block_size,
                    cs_sles_it_convergence_t  *convergence,
                    const cs_real_t           *rhs,
                    cs_real_t                 *restrict vx,
                    size_t                     aux_size,
                    void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg;
  double  ro_0, ro_1, alpha, rk_gkm1, rk_gk, beta, residual;
  cs_real_t  *_aux_vectors;
  cs_real_t  *restrict rk, *restrict dk, *restrict gk;
  cs_real_t *restrict zk;

  unsigned n_iter = 0;

  /* Allocate or map work arrays */
  /*-----------------------------*/

  assert(c->setup_data != NULL);

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 4;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      BFT_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    rk = _aux_vectors;
    dk = _aux_vectors + wa_size;
    gk = _aux_vectors + wa_size*2;
    zk = _aux_vectors + wa_size*3;
  }

  /* Initialize iterative calculation */
  /*----------------------------------*/

  /* Residual and descent direction */

  cs_matrix_vector_multiply(a, vx, rk);  /* rk = A.x0 */

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_rows; ii++)
    rk[ii] -= rhs[ii];

  /* Preconditioning */

  c->setup_data->pc_apply(c->setup_data->pc_context,
                          rk,
                          gk);

  /* Descent direction */
  /*-------------------*/

#if defined(HAVE_OPENMP)

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_rows; ii++)
    dk[ii] = gk[ii];

#else

  memcpy(dk, gk, n_rows * sizeof(cs_real_t));

#endif

  _dot_products_xx_xy(c, rk, gk, &residual, &rk_gkm1);
  residual = sqrt(residual);

  /* If no solving required, finish here */

  c->setup_data->initial_residual = residual;
  cvg = _convergence_test(c, n_iter, residual, convergence);

  if (cvg == CS_SLES_ITERATING) {

    n_iter = 1;

    cs_matrix_vector_multiply(a, dk, zk);

    /* Descent parameter */

    _dot_products_xy_yz(c, rk, dk, zk, &ro_0, &ro_1);

    cs_real_t d_ro_1 = (CS_ABS(ro_1) > DBL_MIN) ? 1. / ro_1 : 0.;
    alpha =  - ro_0 * d_ro_1;

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        vx[ii] += (alpha * dk[ii]);

#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        rk[ii] += (alpha * zk[ii]);
    }

    /* Convergence test */

    residual = sqrt(_dot_product_xx(c, rk));
    cvg = _convergence_test(c, n_iter, residual, convergence);

    /* Current Iteration */
    /*-------------------*/

  }

  while (cvg == CS_SLES_ITERATING) {

    /* Preconditioning */

    c->setup_data->pc_apply(c->setup_data->pc_context, rk, gk);

    /* Compute residual and prepare descent parameter */

    _dot_products_xx_xy(c, rk, gk, &residual, &rk_gk);

    residual = sqrt(residual);

    /* Convergence test for end of previous iteration */

    if (n_iter > 1)
      cvg = _convergence_test(c, n_iter, residual, convergence);

    if (cvg != CS_SLES_ITERATING)
      break;

    n_iter += 1;

    /* Complete descent parameter computation and matrix.vector product */

    beta = (CS_ABS(rk_gkm1) > DBL_MIN) ? rk_gk / rk_gkm1 : 0.;
    rk_gkm1 = rk_gk;

#   pragma omp parallel for firstprivate(alpha) if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      dk[ii] = gk[ii] + (beta * dk[ii]);

    cs_matrix_vector_multiply(a, dk, zk);

    _dot_products_xy_yz(c, rk, dk, zk, &ro_0, &ro_1);

    cs_real_t d_ro_1 = (CS_ABS(ro_1) > DBL_MIN) ? 1. / ro_1 : 0.;
    alpha =  - ro_0 * d_ro_1;

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        vx[ii] += (alpha * dk[ii]);

#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        rk[ii] += (alpha * zk[ii]);
    }

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using flexible preconditioned conjugate gradient.
 *
 * Compared to standard PCG, FCG supports variable preconditioners.
 *
 * This variant, described in \cite Notay:2015, allows computing the
 * required inner products with a single global communication.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- matrix
 *   diag_block_size <-- diagonal block size
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_flexible_conjugate_gradient(cs_sles_it_t              *c,
                             const cs_matrix_t         *a,
                             cs_lnum_t                  diag_block_size,
                             cs_sles_it_convergence_t  *convergence,
                             const cs_real_t           *rhs,
                             cs_real_t                 *restrict vx,
                             size_t                     aux_size,
                             void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg = CS_SLES_ITERATING;
  cs_real_t  *_aux_vectors;
  cs_real_t  *restrict rk, *restrict vk, *restrict wk;
  cs_real_t  *restrict dk, *restrict qk;

  unsigned n_iter = 0;

  /* Allocate or map work arrays */
  /*-----------------------------*/

  assert(c->setup_data != NULL);

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 5;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      BFT_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    rk = _aux_vectors;
    vk = _aux_vectors + wa_size;
    wk = _aux_vectors + wa_size*2;
    dk = _aux_vectors + wa_size*3;
    qk = _aux_vectors + wa_size*4;
  }

  /* Initialize iterative calculation */
  /*----------------------------------*/

  /* Residual and descent direction */

  cs_matrix_vector_multiply(a, vx, rk);  /* rk = A.x0 */

# pragma omp parallel if(n_rows > CS_THR_MIN)
  {
#   pragma omp for nowait
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      rk[ii] = rhs[ii] - rk[ii];

#   pragma omp for nowait
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      qk[ii] = 0;
  }

  double rho_km1 = 0;

  while (cvg == CS_SLES_ITERATING) {

    /* Preconditioning */

    c->setup_data->pc_apply(c->setup_data->pc_context, rk, vk);

    cs_matrix_vector_multiply(a, vk, wk);

    /* Compute residual and prepare descent parameter */

    double alpha_k, beta_k, gamma_k, residual;

    _dot_products_vr_vw_vq_rr(c, vk, rk, wk, qk,
                              &alpha_k, &beta_k, &gamma_k, &residual);

    residual = sqrt(residual);

    /* Convergence test for end of previous iteration */

    if (n_iter > 0)
      cvg = _convergence_test(c, n_iter, residual, convergence);
    else
      c->setup_data->initial_residual = residual;

    if (cvg != CS_SLES_ITERATING)
      break;

    /* Complete descent parameter computation and matrix.vector product */

    if (n_iter > 0) {

      cs_real_t gk_rk1 = (CS_ABS(rho_km1) > DBL_MIN) ? gamma_k / rho_km1 : 0.;
      cs_real_t rho_k = beta_k - gamma_k * gk_rk1;
      cs_real_t ak_rk = (CS_ABS(rho_k) > DBL_MIN) ? alpha_k / rho_k : 0.;

#     pragma omp parallel if(n_rows > CS_THR_MIN)
      {
#       pragma omp for nowait
        for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
          dk[ii] = vk[ii] - gk_rk1 * dk[ii];
          vx[ii] = vx[ii] + ak_rk * dk[ii];
        }

#       pragma omp for nowait
        for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
          qk[ii] = wk[ii] - gk_rk1 * qk[ii];
          rk[ii] = rk[ii] - ak_rk * qk[ii];
        }
      }

      rho_km1 = rho_k;
    }
    else { /* n_iter == 0 */

      cs_real_t rho_k = beta_k;
      cs_real_t ak_rk = (CS_ABS(rho_k) > DBL_MIN) ? alpha_k / rho_k : 0.;

#     pragma omp parallel if(n_rows > CS_THR_MIN)
      {
#       pragma omp for nowait
        for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
          dk[ii] = vk[ii];
          vx[ii] = vx[ii] + ak_rk * vk[ii];
        }

#       pragma omp for nowait
        for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
          qk[ii] = wk[ii];
          rk[ii] = rk[ii] - ak_rk * wk[ii];
        }
      }

      rho_km1 = rho_k;
    }

    n_iter += 1;
  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using flexible preconditioned conjugate gradient.
 *
 * Compared to standard PCG, IPCG supports variable preconditioners, at
 * the expense of storing the residual at iterations k (rk) and k-1 (rkm1)
 * to compute the Beta coefficient. When the preconditioner is constant
 * across the iterations, IPCG is equivalent to PCG.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- matrix
 *   diag_block_size <-- diagonal block size
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_conjugate_gradient_ip(cs_sles_it_t              *c,
                       const cs_matrix_t         *a,
                       cs_lnum_t                  diag_block_size,
                       cs_sles_it_convergence_t  *convergence,
                       const cs_real_t           *rhs,
                       cs_real_t                 *restrict vx,
                       size_t                     aux_size,
                       void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg;
  double  ro_0, ro_1, alpha, rk_gk_m1, rkm1_gk, rk_gk, beta, residual;
  cs_real_t  *_aux_vectors;
  cs_real_t  *restrict rk, *restrict rkm1, *restrict dk, *restrict gk;
  cs_real_t  *restrict zk;

  unsigned n_iter = 0;

  /* Allocate or map work arrays */
  /*-----------------------------*/

  assert(c->setup_data != NULL);

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 5;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      BFT_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    rk    = _aux_vectors;
    rkm1  = _aux_vectors + wa_size;
    dk    = _aux_vectors + wa_size*2;
    gk    = _aux_vectors + wa_size*3;
    zk    = _aux_vectors + wa_size*4;
  }

  /* Initialize iterative calculation */
  /*----------------------------------*/

  /* Residual and descent direction */

  cs_matrix_vector_multiply(a, vx, rk);  /* rk = A.x0 */

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_rows; ii++)
    rk[ii] -= rhs[ii];

  /* Preconditioning */

  c->setup_data->pc_apply(c->setup_data->pc_context, rk, gk);

  /* Descent direction */
  /*-------------------*/

#if defined(HAVE_OPENMP)

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_rows; ii++)
    dk[ii] = gk[ii];

#else

  memcpy(dk, gk, n_rows * sizeof(cs_real_t));

#endif

  _dot_products_xx_xy(c, rk, gk, &residual, &rk_gk_m1);
  residual = sqrt(residual);

  /* If no solving required, finish here */

  c->setup_data->initial_residual = residual;
  cvg = _convergence_test(c, n_iter, residual, convergence);

  if (cvg == CS_SLES_ITERATING) {

    n_iter = 1;

    cs_matrix_vector_multiply(a, dk, zk);

    /* Descent parameter */

    _dot_products_xy_yz(c, rk, dk, zk, &ro_0, &ro_1);

    cs_real_t d_ro_1 = (CS_ABS(ro_1) > DBL_MIN) ? 1. / ro_1 : 0.;
    alpha =  - ro_0 * d_ro_1;

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        vx[ii] += (alpha * dk[ii]);

#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
        rkm1[ii] = rk[ii];
        rk[ii] += (alpha * zk[ii]);
      }
    }

    /* Convergence test */

    residual = sqrt(_dot_product_xx(c, rk));
    cvg = _convergence_test(c, n_iter, residual, convergence);

    /* Current Iteration */
    /*-------------------*/

  }

  while (cvg == CS_SLES_ITERATING) {

    /* Preconditioning */

    c->setup_data->pc_apply(c->setup_data->pc_context, rk, gk);

    /* Compute residual and prepare descent parameter */

    _dot_products_xx_xy_yz(c, rk, gk, rkm1, &residual, &rk_gk, &rkm1_gk);

    residual = sqrt(residual);

    /* Convergence test for end of previous iteration */

    if (n_iter > 1)
      cvg = _convergence_test(c, n_iter, residual, convergence);

    if (cvg != CS_SLES_ITERATING)
      break;

    n_iter += 1;

    /* Complete descent parameter computation and matrix.vector product */

    beta = (CS_ABS(rk_gk_m1) > DBL_MIN) ? (rk_gk - rkm1_gk) / rk_gk_m1 : 0.;
    rk_gk_m1 = rk_gk;

#   pragma omp parallel for firstprivate(alpha) if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      dk[ii] = gk[ii] + (beta * dk[ii]);

    cs_matrix_vector_multiply(a, dk, zk);

    _dot_products_xy_yz(c, rk, dk, zk, &ro_0, &ro_1);

    cs_real_t d_ro_1 = (CS_ABS(ro_1) > DBL_MIN) ? 1. / ro_1 : 0.;
    alpha =  - ro_0 * d_ro_1;

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        vx[ii] += (alpha * dk[ii]);

#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
        rkm1[ii] = rk[ii];
        rk[ii] += (alpha * zk[ii]);
      }
    }

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using preconditioned conjugate gradient
 * with single reduction.
 *
 * For more information, see Lapack Working note 56, at
 * http://www.netlib.org/lapack/lawnspdf/lawn56.pdf)
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- matrix
 *   diag_block_size <-- block size of element ii, ii
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_conjugate_gradient_sr(cs_sles_it_t              *c,
                       const cs_matrix_t         *a,
                       cs_lnum_t                  diag_block_size,
                       cs_sles_it_convergence_t  *convergence,
                       const cs_real_t           *rhs,
                       cs_real_t                 *restrict vx,
                       size_t                     aux_size,
                       void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg;
  double  ro_0, ro_1, alpha, rk_gkm1, rk_gk, gk_sk, beta, residual;
  cs_real_t *_aux_vectors;
  cs_real_t  *restrict rk, *restrict dk, *restrict gk, *restrict sk;
  cs_real_t *restrict zk;

  unsigned n_iter = 0;

  /* Allocate or map work arrays */
  /*-----------------------------*/

  assert(c->setup_data != NULL);

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 5;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      BFT_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    rk = _aux_vectors;
    dk = _aux_vectors + wa_size;
    gk = _aux_vectors + wa_size*2;
    zk = _aux_vectors + wa_size*3;
    sk = _aux_vectors + wa_size*4;
  }

  /* Initialize iterative calculation */
  /*----------------------------------*/

  /* Residual and descent direction */

  cs_matrix_vector_multiply(a, vx, rk);  /* rk = A.x0 */

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_rows; ii++)
    rk[ii] -= rhs[ii];

  /* Preconditionning */

  c->setup_data->pc_apply(c->setup_data->pc_context, rk, gk);

  /* Descent direction */
  /*-------------------*/

#if defined(HAVE_OPENMP)

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_rows; ii++)
    dk[ii] = gk[ii];

#else

  memcpy(dk, gk, n_rows * sizeof(cs_real_t));

#endif

  cs_matrix_vector_multiply(a, dk, zk); /* zk = A.dk */

  /* Descent parameter */

  _dot_products_xx_xy_yz(c, rk, dk, zk, &residual, &ro_0, &ro_1);
  residual = sqrt(residual);

  c->setup_data->initial_residual = residual;

  /* If no solving required, finish here */

  cvg = _convergence_test(c, n_iter, residual, convergence);

  if (cvg == CS_SLES_ITERATING) {

    n_iter = 1;

    cs_real_t d_ro_1 = (CS_ABS(ro_1) > DBL_MIN) ? 1. / ro_1 : 0.;
    alpha =  - ro_0 * d_ro_1;

    rk_gkm1 = ro_0;

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        vx[ii] += (alpha * dk[ii]);

#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        rk[ii] += (alpha * zk[ii]);
    }

    /* Convergence test */

    residual = sqrt(_dot_product_xx(c, rk));

    cvg = _convergence_test(c, n_iter, residual, convergence);

  }

  /* Current Iteration */
  /*-------------------*/

  while (cvg == CS_SLES_ITERATING) {

    /* Preconditionning */

    c->setup_data->pc_apply(c->setup_data->pc_context, rk, gk);

    cs_matrix_vector_multiply(a, gk, sk);  /* sk = A.gk */

    /* Compute residual and prepare descent parameter */

    _dot_products_xx_xy_yz(c, rk, gk, sk, &residual, &rk_gk, &gk_sk);

    residual = sqrt(residual);

    /* Convergence test for end of previous iteration */

    if (n_iter > 1)
      cvg = _convergence_test(c, n_iter, residual, convergence);

    if (cvg != CS_SLES_ITERATING)
      break;

    n_iter += 1;

    /* Complete descent parameter computation and matrix.vector product */

    beta = (CS_ABS(rk_gkm1) > DBL_MIN) ? rk_gk / rk_gkm1 : 0.;
    rk_gkm1 = rk_gk;

    ro_1 = gk_sk - beta*beta*ro_1;
    ro_0 = rk_gk;

    cs_real_t d_ro_1 = (CS_ABS(ro_1) > DBL_MIN) ? 1. / ro_1 : 0.;
    alpha =  - ro_0 * d_ro_1;

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
        dk[ii] = gk[ii] + (beta * dk[ii]);
        vx[ii] += alpha * dk[ii];
      }
#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
        zk[ii] = sk[ii] + (beta * zk[ii]);
        rk[ii] += alpha * zk[ii];
      }
    }

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using non-preconditioned conjugate gradient.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- matrix
 *   diag_block_size <-- block size of element ii, ii
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_conjugate_gradient_npc(cs_sles_it_t              *c,
                        const cs_matrix_t         *a,
                        cs_lnum_t                  diag_block_size,
                        cs_sles_it_convergence_t  *convergence,
                        const cs_real_t           *rhs,
                        cs_real_t                 *restrict vx,
                        size_t                     aux_size,
                        void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg;
  double  ro_0, ro_1, alpha, rk_rkm1, rk_rk, beta, residual;
  cs_real_t *_aux_vectors;
  cs_real_t  *restrict rk, *restrict dk, *restrict zk;

  unsigned n_iter = 0;

  /* Allocate or map work arrays */
  /*-----------------------------*/

  assert(c->setup_data != NULL);

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 3;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      BFT_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    rk = _aux_vectors;
    dk = _aux_vectors + wa_size;
    zk = _aux_vectors + wa_size*2;
  }

  /* Initialize iterative calculation */
  /*----------------------------------*/

  /* Residual and descent direction */

  cs_matrix_vector_multiply(a, vx, rk);  /* rk = A.x0 */

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_rows; ii++)
    rk[ii] -= rhs[ii];

  /* Descent direction */
  /*-------------------*/

#if defined(HAVE_OPENMP)

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_rows; ii++)
    dk[ii] = rk[ii];

#else

  memcpy(dk, rk, n_rows * sizeof(cs_real_t));

#endif

  rk_rkm1 = _dot_product_xx(c, rk);
  residual = sqrt(rk_rkm1);

  /* If no solving required, finish here */

  c->setup_data->initial_residual = residual;
  cvg = _convergence_test(c, n_iter, residual, convergence);

  if (cvg == CS_SLES_ITERATING) {

    n_iter = 1;

    cs_matrix_vector_multiply(a, dk, zk);

    /* Descent parameter */

    _dot_products_xy_yz(c, rk, dk, zk, &ro_0, &ro_1);

    cs_real_t d_ro_1 = (CS_ABS(ro_1) > DBL_MIN) ? 1. / ro_1 : 0.;
    alpha =  - ro_0 * d_ro_1;

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        vx[ii] += (alpha * dk[ii]);

#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        rk[ii] += (alpha * zk[ii]);
    }

    /* Convergence test */

    residual = sqrt(_dot_product(c, rk, rk));

    cvg = _convergence_test(c, n_iter, residual, convergence);

  }

  /* Current Iteration */
  /*-------------------*/

  while (cvg == CS_SLES_ITERATING) {

    /* Compute residual and prepare descent parameter */

    _dot_products_xx_xy(c, rk, rk, &residual, &rk_rk);

    residual = sqrt(residual);

    /* Convergence test for end of previous iteration */

    if (n_iter > 1)
      cvg = _convergence_test(c, n_iter, residual, convergence);

    if (cvg != CS_SLES_ITERATING)
      break;

    n_iter += 1;

    /* Complete descent parameter computation and matrix.vector product */

    beta = (CS_ABS(rk_rkm1) > DBL_MIN) ? rk_rk / rk_rkm1 : 0.;
    rk_rkm1 = rk_rk;

#   pragma omp parallel for firstprivate(alpha) if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      dk[ii] = rk[ii] + (beta * dk[ii]);

    cs_matrix_vector_multiply(a, dk, zk);

    _dot_products_xy_yz(c, rk, dk, zk, &ro_0, &ro_1);

    cs_real_t d_ro_1 = (CS_ABS(ro_1) > DBL_MIN) ? 1. / ro_1 : 0.;
    alpha =  - ro_0 * d_ro_1;

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        vx[ii] += (alpha * dk[ii]);

#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        rk[ii] += (alpha * zk[ii]);
    }

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using non-preconditioned conjugate gradient
 * with single reduction.
 *
 * For more information, see Lapack Working note 56, at
 * http://www.netlib.org/lapack/lawnspdf/lawn56.pdf)
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- matrix
 *   diag_block_size <-- block size of element ii, ii
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_conjugate_gradient_npc_sr(cs_sles_it_t              *c,
                           const cs_matrix_t         *a,
                           cs_lnum_t                  diag_block_size,
                           cs_sles_it_convergence_t  *convergence,
                           const cs_real_t           *rhs,
                           cs_real_t                 *restrict vx,
                           size_t                     aux_size,
                           void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg;
  double  ro_0, ro_1, alpha, rk_rkm1, rk_rk, rk_sk, beta, residual;
  cs_real_t *_aux_vectors;
  cs_real_t  *restrict rk, *restrict dk, *restrict sk;
  cs_real_t *restrict zk;

  unsigned n_iter = 0;

  /* Allocate or map work arrays */
  /*-----------------------------*/

  assert(c->setup_data != NULL);

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 4;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      BFT_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    rk = _aux_vectors;
    dk = _aux_vectors + wa_size;
    zk = _aux_vectors + wa_size*2;
    sk = _aux_vectors + wa_size*3;
  }

  /* Initialize iterative calculation */
  /*----------------------------------*/

  /* Residual and descent direction */

  cs_matrix_vector_multiply(a, vx, rk);  /* rk = A.x0 */

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_rows; ii++)
    rk[ii] = rk[ii] - rhs[ii];

  /* Descent direction */
  /*-------------------*/

#if defined(HAVE_OPENMP)

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_rows; ii++)
    dk[ii] = rk[ii];

#else

  memcpy(dk, rk, n_rows * sizeof(cs_real_t));

#endif

  cs_matrix_vector_multiply(a, dk, zk); /* zk = A.dk */

  /* Descent parameter */

  _dot_products_xx_xy_yz(c, rk, dk, zk, &residual, &ro_0, &ro_1);
  residual = sqrt(residual);

  /* If no solving required, finish here */

  c->setup_data->initial_residual = residual;

  cvg = _convergence_test(c, n_iter, residual, convergence);

  if (cvg == CS_SLES_ITERATING) {

    n_iter = 1;

    cs_real_t d_ro_1 = (CS_ABS(ro_1) > DBL_MIN) ? 1. / ro_1 : 0.;
    alpha =  - ro_0 * d_ro_1;

    rk_rkm1 = ro_0;

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        vx[ii] += (alpha * dk[ii]);

#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        rk[ii] += (alpha * zk[ii]);
    }

    /* Convergence test */

    residual = sqrt(_dot_product_xx(c, rk));

    cvg = _convergence_test(c, n_iter, residual, convergence);

  }

  /* Current Iteration */
  /*-------------------*/

  while (cvg == CS_SLES_ITERATING) {

    cs_matrix_vector_multiply(a, rk, sk);  /* sk = A.zk */

    /* Compute residual and prepare descent parameter */

    _dot_products_xx_xy(c, rk, sk, &residual, &rk_sk);

    rk_rk = residual;

    residual = sqrt(residual);

    /* Convergence test for end of previous iteration */

    if (n_iter > 1)
      cvg = _convergence_test(c, n_iter, residual, convergence);

    if (cvg != CS_SLES_ITERATING)
      break;

    n_iter += 1;

    /* Complete descent parameter computation and matrix.vector product */

    beta = (CS_ABS(rk_rkm1) > DBL_MIN) ? rk_rk / rk_rkm1 : 0.;
    rk_rkm1 = rk_rk;

    ro_1 = rk_sk - beta*beta*ro_1;
    ro_0 = rk_rk;

    cs_real_t d_ro_1 = (CS_ABS(ro_1) > DBL_MIN) ? 1. / ro_1 : 0.;
    alpha =  - ro_0 * d_ro_1;

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
        dk[ii] = rk[ii] + (beta * dk[ii]);
        vx[ii] += alpha * dk[ii];
      }
#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
        zk[ii] = sk[ii] + (beta * zk[ii]);
        rk[ii] += alpha * zk[ii];
      }
    }

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using preconditioned 3-layer conjugate residual.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- matrix
 *   diag_block_size <-- block size of element ii, ii
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_conjugate_residual_3(cs_sles_it_t              *c,
                      const cs_matrix_t         *a,
                      cs_lnum_t                  diag_block_size,
                      cs_sles_it_convergence_t  *convergence,
                      const cs_real_t           *rhs,
                      cs_real_t                 *restrict vx,
                      size_t                     aux_size,
                      void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg;
  cs_lnum_t  ii;
  double  residual;
  double  ak, bk, ck, dk, ek, denom, alpha, tau;
  cs_real_t *_aux_vectors;
  cs_real_t  *restrict vxm1;
  cs_real_t  *restrict rk, *restrict rkm1;
  cs_real_t  *restrict wk, *restrict zk;
  cs_real_t  *restrict tmp;

  unsigned n_iter = 0;

  /* Allocate or map work arrays */
  /*-----------------------------*/

  assert(c->setup_data != NULL);

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 6;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      BFT_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    vxm1 = _aux_vectors;
    rk = _aux_vectors + wa_size;
    rkm1 = _aux_vectors + wa_size*2;
    tmp = _aux_vectors + wa_size*3;
    wk = _aux_vectors + wa_size*4;
    zk = _aux_vectors + wa_size*5;

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (ii = 0; ii < n_rows; ii++) {
      vxm1[ii] = vx[ii];
      rkm1[ii] = 0.0;
    }
  }

  /* Initialize iterative calculation */
  /*----------------------------------*/

  /* Residual */

  cs_matrix_vector_multiply(a, vx, rk);  /* rk = A.x0 */

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (ii = 0; ii < n_rows; ii++)
    rk[ii] -= rhs[ii];

  residual = _dot_product(c, rk, rk);
  residual = sqrt(residual);

  /* If no solving required, finish here */

  c->setup_data->initial_residual = residual;
  cvg = _convergence_test(c, n_iter, residual, convergence);

  /* Current Iteration */
  /*-------------------*/

  while (cvg == CS_SLES_ITERATING) {

    /* Preconditionning */

    c->setup_data->pc_apply(c->setup_data->pc_context, rk, wk);

    cs_matrix_vector_multiply(a, wk, zk);

    _dot_products_xy_yz(c, rk, zk, rkm1, &ak, &bk);

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (ii = 0; ii < n_rows; ii++)
      tmp[ii] = rk[ii] - rkm1[ii];

    _dot_products_xy_yz(c, rk, tmp, rkm1, &ck, &dk);

    ek = _dot_product_xx(c, zk);

    denom = (ck-dk)*ek - ((ak-bk)*(ak-bk));

    if (fabs(denom) < 1.e-30)
      alpha = 1.0;
    else
      alpha = ((ak-bk)*bk - dk*ek)/denom;

    if (fabs(alpha) < 1.e-30 || fabs(alpha - 1.) < 1.e-30) {
      alpha = 1.0;
      tau = ak/ek;
    }
    else
      tau = ak/ek + ((1 - alpha)/alpha) * bk/ek;

    cs_real_t c0 = (1 - alpha);
    cs_real_t c1 = -alpha*tau;

#   pragma omp parallel firstprivate(alpha, tau, c0, c1) if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (ii = 0; ii < n_rows; ii++) {
        cs_real_t trk = rk[ii];
        rk[ii] = alpha*rk[ii] + c0*rkm1[ii] + c1*zk[ii];
        rkm1[ii] = trk;
      }
#     pragma omp for nowait
      for (ii = 0; ii < n_rows; ii++) {
        cs_real_t tvx = vx[ii];
        vx[ii] = alpha*vx[ii] + c0*vxm1[ii] + c1*wk[ii];
        vxm1[ii] = tvx;
      }
    }

    /* Compute residual */

    residual = sqrt(_dot_product(c, rk, rk));

    /* Convergence test for end of previous iteration */

    if (n_iter > 1)
      cvg = _convergence_test(c, n_iter, residual, convergence);

    if (cvg != CS_SLES_ITERATING)
      break;

    n_iter += 1;

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using Jacobi.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- diagonal block size
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_jacobi(cs_sles_it_t              *c,
        const cs_matrix_t         *a,
        cs_lnum_t                  diag_block_size,
        cs_sles_it_convergence_t  *convergence,
        const cs_real_t           *rhs,
        cs_real_t                 *restrict vx,
        size_t                     aux_size,
        void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg;
  cs_real_t *_aux_vectors;
  cs_real_t *restrict rk;

  double residual = -1;
  unsigned n_iter = 0;

  /* Allocate or map work arrays */
  /*-----------------------------*/

  assert(c->setup_data != NULL);

  const cs_real_t  *restrict ad_inv = c->setup_data->ad_inv;

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 1;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      BFT_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    rk = _aux_vectors;
  }

  const cs_real_t  *restrict ad = cs_matrix_get_diagonal(a);

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
    rk[ii] = vx[ii];
  }

  cvg = CS_SLES_ITERATING;

  /* Current iteration */
  /*-------------------*/

  while (cvg == CS_SLES_ITERATING) {

    n_iter += 1;

    /* Compute Vx <- Vx - (A-diag).Rk and residual. */

    cs_matrix_vector_multiply_partial(a, CS_MATRIX_SPMV_E, rk, vx);

    double  res2 = 0.0;

    if (convergence->precision > 0. || c->plot != NULL) {

#     pragma omp parallel for reduction(+:res2) if(n_rows > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
        vx[ii] = (rhs[ii]-vx[ii])*ad_inv[ii];
        double r = ad[ii] * (vx[ii]-rk[ii]);
        res2 += (r*r);
        rk[ii] = vx[ii];
      }

#if defined(HAVE_MPI)

      if (c->comm != MPI_COMM_NULL) {
        double _sum;
        MPI_Allreduce(&res2, &_sum, 1, MPI_DOUBLE, MPI_SUM, c->comm);
        res2 = _sum;
      }

#endif /* defined(HAVE_MPI) */

      residual = sqrt(res2); /* Actually, residual of previous iteration */

    }
    else {

#     pragma omp parallel for if(n_rows > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
        vx[ii] = (rhs[ii]-vx[ii])*ad_inv[ii];
        rk[ii] = vx[ii];
      }

    }

    /* Convergence test */

    if (n_iter == 1)
      c->setup_data->initial_residual = residual;

    cvg = _convergence_test(c, n_iter, residual, convergence);

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using block Jacobi.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- diagonal block size (unused here)
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              --> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_block_3_jacobi(cs_sles_it_t              *c,
                const cs_matrix_t         *a,
                cs_lnum_t                  diag_block_size,
                cs_sles_it_convergence_t  *convergence,
                const cs_real_t           *rhs,
                cs_real_t                 *restrict vx,
                size_t                     aux_size,
                void                      *aux_vectors)
{
  CS_UNUSED(diag_block_size);

  assert(diag_block_size == 3);

  cs_sles_convergence_state_t cvg;
  double  res2, residual;
  cs_real_t *_aux_vectors;
  cs_real_t  *restrict rk, *restrict vxx;

  unsigned n_iter = 0;

  /* Allocate or map work arrays */
  /*-----------------------------*/

  assert(c->setup_data != NULL);

  const cs_real_t  *restrict ad_inv = c->setup_data->ad_inv;

  const cs_lnum_t n_rows = c->setup_data->n_rows;
  const cs_lnum_t n_blocks = c->setup_data->n_rows / 3;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 2;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      BFT_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    rk  = _aux_vectors;
    vxx = _aux_vectors + wa_size;
  }

  const cs_real_t  *restrict ad = cs_matrix_get_diagonal(a);

  cvg = CS_SLES_ITERATING;

  /* Current iteration */
  /*-------------------*/

  while (cvg == CS_SLES_ITERATING) {

    n_iter += 1;
    memcpy(rk, vx, n_rows * sizeof(cs_real_t));  /* rk <- vx */

    /* Compute vxx <- vx - (a-diag).rk and residual. */

    cs_matrix_vector_multiply_partial(a, CS_MATRIX_SPMV_E, rk, vxx);

    res2 = 0.0;

    /* Compute vx <- diag^-1 . (vxx - rhs) and residual. */
#   pragma omp parallel for reduction(+:res2) if(n_blocks > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_blocks; ii++) {
      _fw_and_bw_lu33(ad_inv + 9*ii,
                      vx + 3*ii,
                      vxx + 3*ii,
                      rhs + 3*ii);
      for (cs_lnum_t jj = 0; jj < 3; jj++) {
        register double r = 0.0;
        for (cs_lnum_t kk = 0; kk < 3; kk++)
          r +=    ad[ii*9 + jj*3 + kk]
               * (vx[ii*3 + kk] - rk[ii*3 + kk]);
        res2 += (r*r);
      }
    }

#if defined(HAVE_MPI)

    if (c->comm != MPI_COMM_NULL) {
      double _sum;
      MPI_Allreduce(&res2, &_sum, 1, MPI_DOUBLE, MPI_SUM, c->comm);
      res2 = _sum;
    }

#endif /* defined(HAVE_MPI) */

    residual = sqrt(res2); /* Actually, residual of previous iteration */

    if (n_iter == 1)
      c->setup_data->initial_residual = residual;

    /* Convergence test */

    cvg = _convergence_test(c, n_iter, residual, convergence);

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);
  return cvg;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using block Jacobi.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- block size of diagonal elements
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              --> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_block_jacobi(cs_sles_it_t              *c,
              const cs_matrix_t         *a,
              cs_lnum_t                  diag_block_size,
              cs_sles_it_convergence_t  *convergence,
              const cs_real_t           *rhs,
              cs_real_t                 *restrict vx,
              size_t                     aux_size,
              void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg;
  double  res2, residual;
  cs_real_t *_aux_vectors;
  cs_real_t  *restrict rk, *restrict vxx;

  unsigned n_iter = 0;

  /* Call setup if not already done, allocate or map work arrays */
  /*-------------------------------------------------------------*/
  assert(c->setup_data != NULL);

  const cs_lnum_t db_size = cs_matrix_get_diag_block_size(a);
  const cs_lnum_t db_size_2 = db_size * db_size;

  const cs_real_t  *restrict ad_inv = c->setup_data->ad_inv;

  const cs_lnum_t n_rows = c->setup_data->n_rows;
  const cs_lnum_t n_blocks = c->setup_data->n_rows / diag_block_size;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 2;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      BFT_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    rk  = _aux_vectors;
    vxx = _aux_vectors + wa_size;
  }

  const cs_real_t  *restrict ad = cs_matrix_get_diagonal(a);

  cvg = CS_SLES_ITERATING;

  /* Current iteration */
  /*-------------------*/

  while (cvg == CS_SLES_ITERATING) {

    n_iter += 1;
    memcpy(rk, vx, n_rows * sizeof(cs_real_t));  /* rk <- vx */

    /* Compute Vx <- Vx - (A-diag).Rk and residual. */

    cs_matrix_vector_multiply_partial(a, CS_MATRIX_SPMV_E, rk, vxx);

    res2 = 0.0;

#   pragma omp parallel for reduction(+:res2) if(n_blocks > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_blocks; ii++) {
      _fw_and_bw_lu(ad_inv + db_size_2*ii,
                    db_size,
                    vx + db_size*ii,
                    vxx + db_size*ii,
                    rhs + db_size*ii);
      for (cs_lnum_t jj = 0; jj < db_size; jj++) {
        register double r = 0.0;
        for (cs_lnum_t kk = 0; kk < db_size; kk++)
          r +=    ad[ii*db_size_2 + jj*db_size + kk]
               * (vx[ii*db_size + kk] - rk[ii*db_size + kk]);
        res2 += (r*r);
      }
    }

#if defined(HAVE_MPI)

    if (c->comm != MPI_COMM_NULL) {
      double _sum;
      MPI_Allreduce(&res2, &_sum, 1, MPI_DOUBLE, MPI_SUM, c->comm);
      res2 = _sum;
    }

#endif /* defined(HAVE_MPI) */

    residual = sqrt(res2); /* Actually, residual of previous iteration */

    if (n_iter == 1)
      c->setup_data->initial_residual = residual;

    /* Convergence test */

    cvg = _convergence_test(c, n_iter, residual, convergence);

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  return cvg;
}


/*----------------------------------------------------------------------------
 * Test for (and eventually report) breakdown.
 *
 * parameters:
 *   c           <-- pointer to solver context info
 *   convergence <-- convergence information structure
 *   coeff       <-- coefficient name
 *   coeff       <-- coefficient to test
 *   epsilon     <-- value to test against
 *   n_iter      <-- current number of iterations
 *   cvg         <-> convergence status
 *
 * returns:
 *   true in case of breakdown, false otherwise
 *----------------------------------------------------------------------------*/

static inline bool
_breakdown(cs_sles_it_t                 *c,
           cs_sles_it_convergence_t     *convergence,
           const char                   *coeff_name,
           double                        coeff,
           double                        epsilon,
           double                        residual,
           int                           n_iter,
           cs_sles_convergence_state_t  *cvg)
{
  bool retval = false;

  if (CS_ABS(coeff) < epsilon) {

    bft_printf
      (_("\n\n"
         "%s [%s]:\n"
         " @@ Warning: non convergence\n"
         "\n"
         "    norm of coefficient \"%s\" is lower than %12.4e\n"
         "\n"
         "    The resolution does not progress anymore."),
       cs_sles_it_type_name[c->type], convergence->name, coeff_name, epsilon);
    bft_printf(_("  n_iter : %5u, res_abs : %11.4e, res_nor : %11.4e\n"),
               n_iter, residual, residual/convergence->r_norm);

    *cvg = CS_SLES_BREAKDOWN;
    retval = true;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using preconditioned Bi-CGSTAB.
 *
 * Parallel-optimized version, groups dot products, at the cost of
 * computation of the preconditionning for n+1 iterations instead of n.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   name            <-- pointer to system name
 *   a               <-- matrix
 *   diag_block_size <-- block size of diagonal elements
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_bi_cgstab(cs_sles_it_t              *c,
           const cs_matrix_t         *a,
           cs_lnum_t                  diag_block_size,
           cs_sles_it_convergence_t  *convergence,
           const cs_real_t           *rhs,
           cs_real_t                 *restrict vx,
           size_t                     aux_size,
           void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg;
  double  ro_0, ro_1, alpha, beta, betam1, gamma, omega, ukres0;
  double  residual;
  cs_real_t  *_aux_vectors;
  cs_real_t  *restrict res0, *restrict rk, *restrict pk, *restrict zk;
  cs_real_t  *restrict uk, *restrict vk;

  unsigned n_iter = 0;

  /* Allocate or map work arrays */
  /*-----------------------------*/

  assert(c->setup_data != NULL);

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 6;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      BFT_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    res0 = _aux_vectors;
    rk = _aux_vectors + wa_size;
    pk = _aux_vectors + wa_size*2;
    zk = _aux_vectors + wa_size*3;
    uk = _aux_vectors + wa_size*4;
    vk = _aux_vectors + wa_size*5;
  }

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
    pk[ii] = 0.0;
    uk[ii] = 0.0;
  }

  /* Initialize iterative calculation */
  /*----------------------------------*/

  cs_matrix_vector_multiply(a, vx, res0);

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
    res0[ii] = -res0[ii] + rhs[ii];
    rk[ii] = res0[ii];
  }

  alpha = 1.0;
  betam1 = 1.0;
  gamma = 1.0;

  cvg = CS_SLES_ITERATING;

  /* Current Iteration */
  /*-------------------*/

  while (cvg == CS_SLES_ITERATING) {

    /* Compute beta and omega;
       group dot products for new iteration's beta
       and previous iteration's residual to reduce total latency */

    if (n_iter == 0) {
      beta = _dot_product_xx(c, rk); /* rk == res0 here */
      residual = sqrt(beta);
      c->setup_data->initial_residual = residual;
    }
    else {
      _dot_products_xx_xy(c, rk, res0, &residual, &beta);
      residual = sqrt(residual);
    }

    /* Convergence test */
    cvg = _convergence_test(c, n_iter, residual, convergence);
    if (cvg != CS_SLES_ITERATING)
      break;

    n_iter += 1;

    if (_breakdown(c, convergence, "beta", beta, _epzero,
                   residual, n_iter, &cvg))
      break;

    if (_breakdown(c, convergence, "alpha", alpha, _epzero,
                   residual, n_iter, &cvg))
      break;

    omega = beta*gamma / (alpha*betam1);
    betam1 = beta;

    /* Compute pk */

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      pk[ii] = rk[ii] + omega*(pk[ii] - alpha*uk[ii]);

    /* Compute zk = c.pk */

    c->setup_data->pc_apply(c->setup_data->pc_context, pk, zk);

    /* Compute uk = A.zk */

    cs_matrix_vector_multiply(a, zk, uk);

    /* Compute uk.res0 and gamma */

    ukres0 = _dot_product(c, uk, res0);

    gamma = beta / ukres0;

    /* First update of vx and rk */

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        vx[ii] += (gamma * zk[ii]);

#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        rk[ii] -= (gamma * uk[ii]);
    }

    /* Compute zk = C.rk (zk is overwritten, vk is a working array */

    c->setup_data->pc_apply(c->setup_data->pc_context, rk, zk);

    /* Compute vk = A.zk and alpha */

    cs_matrix_vector_multiply(a, zk, vk);

    _dot_products_xx_xy(c, vk, rk, &ro_1, &ro_0);

    if (_breakdown(c, convergence, "rho1", ro_1, _epzero,
                   residual, n_iter, &cvg))
      break;

    cs_real_t d_ro_1 = (CS_ABS(ro_1) > DBL_MIN) ? 1. / ro_1 : 0.;
    alpha = ro_0 * d_ro_1;

    /* Final update of vx and rk */

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        vx[ii] += (alpha * zk[ii]);

#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        rk[ii] -= (alpha * vk[ii]);
    }

    /* Convergence test at beginning of next iteration so
       as to group dot products for better parallel performance */
  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Solution of (ad+ax).vx = Rhs using (not yet preconditioned) Bi-CGSTAB2.
 *
 * This Krylov method is based on an implemantation of J.P. Caltagirone
 * in his file bibpac6.f90 for Aquillon. He refers to it as Bi-CGSTAB2,
 * but a proper reference for the method has yet to be found.
 * (Gutknecht's BiCGstab2?)
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- matrix
 *   diag_block_size <-- block size of element ii, ii
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_bicgstab2(cs_sles_it_t              *c,
           const cs_matrix_t         *a,
           cs_lnum_t                  diag_block_size,
           cs_sles_it_convergence_t  *convergence,
           const cs_real_t           *rhs,
           cs_real_t                 *restrict vx,
           size_t                     aux_size,
           void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg;
  double  ro_0, ro_1, alpha, beta, gamma;
  double  omega_1, omega_2, mu, nu, tau;
  double  residual;
  cs_real_t  *_aux_vectors;
  cs_real_t  *restrict res0, *restrict qk, *restrict rk, *restrict sk;
  cs_real_t  *restrict tk, *restrict uk, *restrict vk, *restrict wk;
  cs_real_t  *restrict zk;

  unsigned n_iter = 0;

  /* Allocate or map work arrays */
  /*-----------------------------*/

  assert(c->setup_data != NULL);

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    size_t  n_wa = 9;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      BFT_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    res0 = _aux_vectors;
    zk = _aux_vectors + wa_size;
    qk = _aux_vectors + wa_size*2;
    rk = _aux_vectors + wa_size*3;
    sk = _aux_vectors + wa_size*4;
    tk = _aux_vectors + wa_size*5;
    uk = _aux_vectors + wa_size*6;
    vk = _aux_vectors + wa_size*7;
    wk = _aux_vectors + wa_size*8;
  }

  /* Initialize work arrays */

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
    uk[ii] = 0.0;
  }

  /* Initialize iterative calculation */
  /*----------------------------------*/

  cs_matrix_vector_multiply(a, vx, res0);

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
    res0[ii] = -res0[ii] + rhs[ii];
    rk[ii] = res0[ii];
    qk[ii] = rk[ii];
  }

  ro_0    = 1.0;
  alpha   = 0.0;
  omega_2 = 1.0;

  cvg = CS_SLES_ITERATING;

  /* Current Iteration */
  /*-------------------*/

  while (cvg == CS_SLES_ITERATING) {

    /* Compute beta and omega;
       group dot products for new iteration's beta
       and previous iteration's residual to reduce total latency */

    double mprec = 1.0e-60;

    if (n_iter == 0) {
      residual = sqrt(_dot_product_xx(c, rk)); /* rk == res0 here */
      c->setup_data->initial_residual = residual;
    }
    else
      residual = sqrt(_dot_product_xx(c, rk));

    /* Convergence test */
    cvg = _convergence_test(c, n_iter, residual, convergence);
    if (cvg != CS_SLES_ITERATING)
        break;

    n_iter += 1;

    ro_0 = -omega_2*ro_0;
    ro_1 = _dot_product(c, qk, rk);

    if (_breakdown(c, convergence, "rho0", ro_0, 1.e-60,
                   residual, n_iter, &cvg))
      break;

    if (_breakdown(c, convergence, "rho1", ro_1, _epzero,
                   residual, n_iter, &cvg))
      break;

    beta = alpha*ro_1/ro_0;
    ro_0 = ro_1;

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      uk[ii] = rk[ii] - beta* uk[ii];

    /* Compute vk =  A*uk */

    cs_matrix_vector_multiply(a, uk, vk);

    /* Preconditionning */

    c->setup_data->pc_apply(c->setup_data->pc_context, vk, zk);

    /* Compute gamma and alpha */

    gamma = _dot_product(c, qk, vk);

    if (_breakdown(c, convergence, "gamma", gamma, 1.e-60,
                   residual, n_iter, &cvg))
      break;

    alpha = ro_0/gamma;

    /* Update rk */

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
      rk[ii] -= alpha*vk[ii];
      vx[ii] += alpha*uk[ii];
    }

    /* p = A*r */

    cs_matrix_vector_multiply(a, rk, sk);

    c->setup_data->pc_apply(c->setup_data->pc_context, sk, zk);

    ro_1 = _dot_product(c, qk, sk);
    beta = alpha*ro_1/ro_0;
    ro_0 = ro_1;

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        vk[ii] = sk[ii] - beta*vk[ii];
#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        uk[ii] = rk[ii] - beta*uk[ii];
    }

    /* wk = A*vk */

    cs_matrix_vector_multiply(a, vk, wk);

    c->setup_data->pc_apply(c->setup_data->pc_context, wk, zk);

    gamma = _dot_product(c, qk, wk);
    alpha = (ro_0+mprec)/(gamma+mprec);

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
      rk[ii] -= alpha*vk[ii];
      sk[ii] -= alpha*wk[ii];
    }

    /* tk = A*sk */

    cs_matrix_vector_multiply(a, sk, tk);

    c->setup_data->pc_apply(c->setup_data->pc_context, tk, zk);

    _dot_products_xx_yy_xy_xz_yz(c, sk, tk, rk,
                                 &mu, &tau, &nu, &omega_1, &omega_2);

    tau = tau - (nu*nu/mu);
    omega_2 = (omega_2 - ((nu+mprec)*(omega_1+mprec)/(mu+mprec)))/(tau+mprec);

    omega_1 = (omega_1 - nu*omega_2)/mu;

    /* sol <- sol + omega_1*r + omega_2*s + alpha*u */
#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      vx[ii] += omega_1*rk[ii] + omega_2*sk[ii] + alpha*uk[ii];

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
      /* r <- r - omega_1*s - omega_2*t */
#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        rk[ii] += - omega_1*sk[ii] - omega_2*tk[ii];
      /* u <- u - omega_1*v - omega_2*w */
#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        uk[ii] += - omega_1*vk[ii] - omega_2*wk[ii];
    }

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Transform using Givens rotations a system Ax=b where A is an upper
 * triangular matrix of size n*(n-1) with a lower diagonal to an equivalent
 * system A'x=b' where A' is an upper triangular matrix.
 *
 * parameters:
 *   a              <-> dense matrix (size: size a_size*a_size):
 *                      a(i,j) = a[i + j*a_size]
 *                      input: A; output: A'
 *   a_size         <-- matrix dim
 *   b              <-> pre-allocated vector of a_size elements
 *                      input: b; output: b'
 *   givens_coeff   <-> pre-allocated vector of a_size elements
 *                      input: previous Givens coefficients
 *                      output: updated Givens coefficients
 *   update_rank    <-- rank of first non null coefficient on lower diagonal
 *   end_update     <-- rank of last non null coefficient on lower diagonal
 *----------------------------------------------------------------------------*/

static void
_givens_rot_update(cs_real_t    *restrict a,
                   cs_lnum_t              a_size,
                   cs_real_t    *restrict b,
                   cs_real_t    *restrict givens_coeff,
                   cs_lnum_t              update_rank,
                   cs_lnum_t              end_update)
{
  for (cs_lnum_t i = 0; i < update_rank; ++i) {
    for (cs_lnum_t j = update_rank; j < end_update; ++j) {

      cs_real_t _aux =   givens_coeff[i]*a[j*a_size + i]
                       + givens_coeff[i + a_size] * a[j*a_size + i+1];

      a[j*a_size + i+1] =   givens_coeff[i] * a[i+1 + j*a_size]
                          - givens_coeff[i + a_size] * a[j*a_size + i];

      a[j*a_size + i] = _aux;
    }
  }

  for (cs_lnum_t i = update_rank; i < end_update; ++i) {

    cs_real_t norm = sqrt(  a[i*a_size + i]   * a[i*a_size + i]
                          + a[i*a_size + i+1] * a[i*a_size + i+1]);

    givens_coeff[a_size + i] = a[i*a_size + i+1]/norm;
    givens_coeff[i] = a[i*a_size + i]/norm;

    b[i+1] = -b[i]*givens_coeff[a_size + i];
    b[i] = b[i]*givens_coeff[i];

    /* i == j */
    {
      cs_real_t _aux =   givens_coeff[i] * a[i*a_size + i]
                       + givens_coeff[a_size+i] * a[i*a_size + i+1];
      a[i*a_size+i+1] = 0;
      a[i*a_size + i] = _aux;
    }

    for (cs_lnum_t j = i+1; j < end_update; j++) {

      cs_real_t _aux =   givens_coeff[i] * a[j*a_size + i]
                       + givens_coeff[a_size+i] * a[j*a_size + i+1];

      a[i+1 + j*a_size] =   givens_coeff[i]*a[i+1 + j*a_size]
                          - givens_coeff[a_size + i]*a[i + j*a_size];

      a[j*a_size + i] = _aux;
    }
  }
}

/*----------------------------------------------------------------------------
 * Compute solution of Ax = b where A is an upper triangular matrix.
 *
 * As the system solved by GMRES will grow with iteration number, we
 * preallocate a allocated size, and solve only the useful part:
 *
 *   |       |       |   |x1|   |b1|
 *   |   A   |  pad  |   |x2|   |b2|
 *   |_______|       |   |x3|   |b3|
 *   |               | * |p | = |p |
 *   |     pad       |   |a |   |a |
 *   |               |   |d |   |d |
 *
 * parameters:
 *   a          <-- dense upper triangular matrix A (size: glob_size*glob_size)
 *                 a(i,j) = a[i + j*glob_size]
 *   a_size     <-- system size
 *   alloc_size <-- a_size + halo size
 *   b          <-- pre-allocated vector of a_size elements
 *   x          --> system solution, pre-allocated vector of a_size elements
 *
 * returns:
 *   0 in case of success, 1 in case of zero-pivot.
 *----------------------------------------------------------------------------*/

static int
_solve_diag_sup_halo(cs_real_t  *restrict a,
                     cs_lnum_t            a_size,
                     cs_lnum_t            alloc_size,
                     cs_real_t  *restrict b,
                     cs_real_t  *restrict x)
{
  for (cs_lnum_t i = a_size - 1; i > -1; i--) {

    x[i] = b[i];

    for (cs_lnum_t j = i + 1; j < a_size; j++)
      x[i] = x[i] - a[j*alloc_size + i]*x[j];

    x[i] /= a[i*alloc_size + i];
  }

  return 0; /* We should add a check for zero-pivot */
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using optimised preconditioned GCR.
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

#if defined(__INTEL_COMPILER)
#pragma optimization_level 2 /* Bug with O3 or above at least with icc 19
                                on Xeon(R) Platinum 8260 */
#endif

static cs_sles_convergence_state_t
_gcr(cs_sles_it_t              *c,
     const cs_matrix_t         *a,
     cs_lnum_t                  diag_block_size,
     cs_sles_it_convergence_t  *convergence,
     const cs_real_t           *rhs,
     cs_real_t                 *restrict vx,
     size_t                     aux_size,
     void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg = CS_SLES_ITERATING;

  cs_real_t *_aux_vectors = NULL, *alpha = NULL;
  cs_real_t *restrict rk, *restrict zk, *restrict ck;
  cs_real_t *restrict gkj, *restrict gkj_inv;

  /* In case of the standard GCR, n_k_per_restart --> Inf,
   * or stops until convergence*/
  const unsigned n_k_per_restart = c->restart_interval;

  size_t wa_size;

  unsigned n_restart = 0;

  /* Allocate or map work arrays */
  /*-----------------------------*/

  assert(c->setup_data != NULL);
  const cs_lnum_t n_rows = c->setup_data->n_rows;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 1 + n_k_per_restart * 2;
    wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      BFT_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    rk = _aux_vectors;                               /* store residuals  */
    zk = _aux_vectors + wa_size;                     /* store inv(M)*r   */
    ck = _aux_vectors + wa_size * (1 + n_k_per_restart);   /* store A*zk */
  }

  BFT_MALLOC(alpha, n_k_per_restart, cs_real_t);

  /* gkj stores the upper triangle matrix Gamma of crossed inner-products
   * Not necessary to allocate the full matrix size
   * gkj_inv stores the inverse of gkj */

  BFT_MALLOC(gkj, (n_k_per_restart + 1) * n_k_per_restart / 2, cs_real_t);
  BFT_MALLOC(gkj_inv, (n_k_per_restart + 1) * n_k_per_restart / 2, cs_real_t);

  /* Current Restart */
  while (cvg == CS_SLES_ITERATING) {

    unsigned n_c_iter = 0;

    /* Initialize iterative calculation */
    /*----------------------------------*/

    cs_matrix_vector_multiply(a, vx, rk);  /* rk = A.x0 */

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      rk[ii] -= rhs[ii];

    double residual = sqrt(_dot_product_xx(c, rk));

    if (n_restart == 0)
      c->setup_data->initial_residual = residual;

    /* Current Iteration on k */
    /* ---------------------- */

    while (cvg == CS_SLES_ITERATING && n_c_iter < n_k_per_restart) {

      /* Preconditionning */

      cs_real_t *zk_n = zk + n_c_iter * wa_size;
      cs_real_t *ck_n = ck + n_c_iter * wa_size;

      c->setup_data->pc_apply(c->setup_data->pc_context, rk, zk_n);

      cs_matrix_vector_multiply(a, zk_n, ck_n);

      for (cs_lnum_t jj = 0; jj < (int)n_c_iter; jj++) {
        cs_real_t *ck_j = ck + jj * wa_size;

        cs_lnum_t ii_jn = (n_c_iter + 1) * n_c_iter / 2 + jj;
        gkj[ii_jn] = _dot_product(c, ck_j, ck_n);
#       pragma omp parallel for if(n_rows > CS_THR_MIN)
        for (cs_lnum_t ii = 0; ii < n_rows; ii++)
          ck_n[ii] += - gkj[ii_jn] * ck_j[ii];
      }

      const int  iter_shift = (n_c_iter+1) * n_c_iter / 2 + n_c_iter;
      gkj[iter_shift] = sqrt(_dot_product(c, ck_n, ck_n));

      if (fabs(gkj[iter_shift]) > 0) {

        cs_real_t scale = 1. / gkj[iter_shift];
#       pragma omp parallel for if(n_rows > CS_THR_MIN)
        for (cs_lnum_t ii = 0; ii < n_rows; ii++)
          ck_n[ii] *= scale;

        alpha[n_c_iter] = _dot_product(c, ck_n, rk);

      }
      else
        alpha[n_c_iter] = 0.;

#     pragma omp parallel for if(n_rows > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        rk[ii] += - alpha[n_c_iter] * ck_n[ii];

      /* Compute residual */

      residual = sqrt(_dot_product_xx(c, rk));

      n_c_iter += 1;

      /* Convergence test of current iteration */

      cvg = _convergence_test(c, (n_restart * n_k_per_restart) + n_c_iter,
                              residual, convergence);

      if (cvg != CS_SLES_ITERATING)
        break;

    } /* Needs iterating or k < n_restart */

    /* Inversion of Gamma */

    if (n_c_iter == 1 && !(fabs(alpha[0]) > 0))
      gkj_inv[0] = 1.0;

    else {

      cs_lnum_t n_g_inv = (n_k_per_restart + 1) * n_k_per_restart / 2;
      for (cs_lnum_t jj = 0; jj < n_g_inv; jj++)
        gkj_inv[jj] = 0.0;

      for (cs_lnum_t kk = 0; kk < (int)n_c_iter; kk++) {
        for (cs_lnum_t ii = 0; ii < kk; ii++) {
          for (cs_lnum_t jj = 0; jj < kk; jj++)
            gkj_inv[(kk + 1) * kk / 2 + ii]
              +=   ((ii <= jj) ? gkj_inv[(jj + 1) * jj / 2 + ii] : 0.0)
                 * gkj[(kk + 1) * kk / 2  + jj];
        }

        for (cs_lnum_t jj = 0; jj < kk; jj++)
          gkj_inv[(kk + 1) * kk / 2 + jj] /= - gkj[(kk + 1) * kk / 2 + kk];

        gkj_inv[(kk + 1) * kk / 2 + kk] = 1.0 / gkj[(kk + 1) * kk / 2 + kk];
      }

    } /* n_c_iter > 1 */

    /* Compute the solution */

#   pragma omp parallel if (n_rows > CS_THR_MIN)
    {
      cs_lnum_t s_id, e_id;
      cs_parall_thread_range(n_rows, sizeof(cs_real_t), &s_id, &e_id);

      for (cs_lnum_t kk = 0; kk < (int)n_c_iter; kk++) {
        for(cs_lnum_t jj = 0; jj <= kk; jj++) {
          const cs_real_t *zk_j = zk + jj*wa_size;
          for (cs_lnum_t ii = s_id; ii < e_id; ii++)
            vx[ii] -=    alpha[kk] * zk_j[ii]
                      *  gkj_inv[(kk + 1) * kk / 2 + jj];
        }
      }
    }

    n_restart += 1;

  } /* Needs iterating */

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  BFT_FREE(alpha);
  BFT_FREE(gkj);
  BFT_FREE(gkj_inv);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using preconditioned GMRES.
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

static cs_sles_convergence_state_t
_gmres(cs_sles_it_t              *c,
       const cs_matrix_t         *a,
       cs_lnum_t                  diag_block_size,
       cs_sles_it_convergence_t  *convergence,
       const cs_real_t           *rhs,
       cs_real_t                 *restrict vx,
       size_t                     aux_size,
       void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg = CS_SLES_ITERATING;
  int l_iter, l_old_iter;
  double    beta, dot_prod, residual;
  cs_real_t  *_aux_vectors;
  cs_real_t *restrict _krylov_vectors, *restrict _h_matrix;
  cs_real_t *restrict _givens_coeff, *restrict _beta;
  cs_real_t *restrict dk, *restrict gk;
  cs_real_t *restrict bk, *restrict fk, *restrict krk;

  cs_lnum_t krylov_size_max = c->restart_interval;
  unsigned n_iter = 0;

  /* Allocate or map work arrays */
  /*-----------------------------*/

  assert(c->setup_data != NULL);

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  /* Allocate work arrays */

  int krylov_size = sqrt(n_rows*diag_block_size)*1.5 + 1;
  if (krylov_size > krylov_size_max)
    krylov_size = krylov_size_max;

#if defined(HAVE_MPI)
  if (c->comm != MPI_COMM_NULL) {
    int _krylov_size = krylov_size;
    MPI_Allreduce(&_krylov_size,
                  &krylov_size,
                  1,
                  MPI_INT,
                  MPI_MIN,
                  c->comm);
  }
#endif

  int     check_freq = (int)(krylov_size/10) + 1;
  double  epsi = 1.e-15;
  int scaltest = 0;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;

    size_t _aux_r_size;
    size_t  n_wa = 4;
    size_t  wa_size = n_cols < krylov_size? krylov_size : n_cols;

    wa_size = CS_SIMD_SIZE(wa_size);
    _aux_r_size =   wa_size*n_wa
                  + (krylov_size-1)*(n_rows + krylov_size) + 3*krylov_size;

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < _aux_r_size)
      BFT_MALLOC(_aux_vectors, _aux_r_size, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    dk = _aux_vectors;
    gk = _aux_vectors + wa_size;
    bk = _aux_vectors + 2*wa_size;
    fk = _aux_vectors + 3*wa_size;
    _krylov_vectors = _aux_vectors + n_wa*wa_size;
    _h_matrix = _aux_vectors + n_wa*wa_size + (krylov_size - 1)*n_rows;
    _givens_coeff =   _aux_vectors + n_wa*wa_size
                    + (krylov_size - 1)*(n_rows + krylov_size);
    _beta =   _aux_vectors + n_wa*wa_size
            + (krylov_size - 1)*(n_rows + krylov_size) + 2*krylov_size;
  }

  for (cs_lnum_t ii = 0; ii < krylov_size*(krylov_size - 1); ii++)
    _h_matrix[ii] = 0.;

  cvg = CS_SLES_ITERATING;

  while (cvg == CS_SLES_ITERATING) {

    /* Compute  rk <- a*vx (vx = x0) */

    cs_matrix_vector_multiply(a, vx, dk);

    /* Compute  rk <- rhs - rk (r0 = b-A*x0) */

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      dk[ii] = rhs[ii] - dk[ii];

    if (n_iter == 0) {
      residual = sqrt(_dot_product_xx(c, dk));
      c->setup_data->initial_residual = residual;
      cvg = _convergence_test(c, n_iter, residual, convergence);
      if (cvg != CS_SLES_ITERATING)
        break;
    }

    /* beta = ||r0|| */
    beta = sqrt(_dot_product(c, dk, dk));
    dot_prod = beta;

    _beta[0] = beta;
    for (cs_lnum_t ii = 1; ii < krylov_size; ii++)
      _beta[ii] = 0.;

    /* Lap */

    l_iter = 0;
    l_old_iter = 0;

    for (cs_lnum_t ii = 0; ii < krylov_size - 1; ii++) {

      /* krk = k-ieth col of _krylov_vector = vi */

      krk = _krylov_vectors + ii*n_rows;

#     pragma omp parallel for if(n_rows > CS_THR_MIN)
      for (cs_lnum_t jj = 0; jj < n_rows; jj++)
        krk[jj] = dk[jj]/dot_prod;

      c->setup_data->pc_apply(c->setup_data->pc_context, krk, gk);

      /* Compute w = dk <- A*vj */

      cs_matrix_vector_multiply(a, gk, dk);

      for (cs_lnum_t jj = 0; jj < ii + 1; jj++) {

        /* Compute h(k,i) = <w,vi> = <dk,vi> */

        _h_matrix[ii*krylov_size + jj]
          = _dot_product(c, dk, (_krylov_vectors + jj*n_rows));

        /* Compute w = dk <- w - h(i,k)*vi */

        cs_axpy(n_rows,
                -_h_matrix[ii*krylov_size+jj],
                (_krylov_vectors + jj*n_rows),
                dk);

      }

      /* Compute h(i+1,i) = sqrt<w,w> */

      dot_prod = sqrt(_dot_product(c, dk, dk));
      _h_matrix[ii*krylov_size + ii + 1] = dot_prod;

      if (dot_prod < epsi) scaltest = 1;

      if (   (l_iter + 1)%check_freq == 0
          || l_iter == krylov_size - 2
          || scaltest == 1) {

          /* H matrix to diagonal sup matrix */

        _givens_rot_update(_h_matrix,
                           krylov_size,
                           _beta,
                           _givens_coeff,
                           l_old_iter,
                           l_iter + 1);

        l_old_iter = l_iter + 1;

        /* Solve diag sup system */

        _solve_diag_sup_halo(_h_matrix, l_iter + 1, krylov_size, _beta, gk);

#       pragma omp parallel for if(n_rows > CS_THR_MIN)
        for (cs_lnum_t jj = 0; jj < n_rows; jj++) {
          fk[jj] = 0.0;
          for (cs_lnum_t kk = 0; kk <= l_iter; kk++)
            fk[jj] += _krylov_vectors[kk*n_rows + jj] * gk[kk];
        }

        c->setup_data->pc_apply(c->setup_data->pc_context, fk, gk);

#       pragma omp parallel for if(n_rows > CS_THR_MIN)
        for (cs_lnum_t jj = 0; jj < n_rows; jj++)
          fk[jj] = vx[jj] + gk[jj];

        cs_matrix_vector_multiply(a, fk, bk);

        /* Compute residual = | Ax - b |_1 */

#       pragma omp parallel for if(n_rows > CS_THR_MIN)
        for (cs_lnum_t jj = 0; jj < n_rows; jj++)
          bk[jj] -= rhs[jj];

        residual = sqrt(_dot_product_xx(c, bk));

        cvg = _convergence_test(c, n_iter, residual, convergence);

      }

      n_iter++;
      l_iter++;

      if (cvg == CS_SLES_CONVERGED || cvg == CS_SLES_MAX_ITERATION ||
          l_iter == krylov_size - 1   || scaltest == 1) {
#       pragma omp parallel for if (n_rows > CS_THR_MIN)
        for (cs_lnum_t jj = 0; jj < n_rows; jj++)
          vx[jj] = fk[jj];
        break;
      }
    }
  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using Process-local Gauss-Seidel.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- diagonal block size
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (unused here)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_p_ordered_gauss_seidel_msr(cs_sles_it_t              *c,
                            const cs_matrix_t         *a,
                            cs_lnum_t                  diag_block_size,
                            cs_sles_it_convergence_t  *convergence,
                            const cs_real_t           *rhs,
                            cs_real_t                 *restrict vx)
{
  cs_sles_convergence_state_t cvg;
  double  res2, residual;

  unsigned n_iter = 0;

  const cs_lnum_t n_rows = cs_matrix_get_n_rows(a);

  const cs_halo_t *halo = cs_matrix_get_halo(a);

  const cs_real_t  *restrict ad_inv = c->setup_data->ad_inv;

  const cs_real_t  *restrict ad = cs_matrix_get_diagonal(a);

  const cs_lnum_t  *a_row_index, *a_col_id;
  const cs_real_t  *a_d_val, *a_x_val;

  const cs_lnum_t db_size = cs_matrix_get_diag_block_size(a);
  const cs_lnum_t db_size_2 = db_size * db_size;

  cs_matrix_get_msr_arrays(a, &a_row_index, &a_col_id, &a_d_val, &a_x_val);

  const cs_lnum_t  *order = c->add_data->order;

  cvg = CS_SLES_ITERATING;

  /* Current iteration */
  /*-------------------*/

  while (cvg == CS_SLES_ITERATING) {

    n_iter += 1;

    /* Synchronize ghost cells first */

    if (halo != NULL)
      cs_matrix_pre_vector_multiply_sync(a, vx);

    /* Compute Vx <- Vx - (A-diag).Rk and residual. */

    res2 = 0.0;

    if (diag_block_size == 1) {

#     pragma omp parallel for reduction(+:res2)      \
                          if(n_rows > CS_THR_MIN && !_thread_debug)
      for (cs_lnum_t ll = 0; ll < n_rows; ll++) {

        cs_lnum_t ii = order[ll];

        const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
        const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
        const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

        cs_real_t vxm1 = vx[ii];
        cs_real_t vx0 = rhs[ii];

        for (cs_lnum_t jj = 0; jj < n_cols; jj++)
          vx0 -= (m_row[jj]*vx[col_id[jj]]);

        vx0 *= ad_inv[ii];

        register double r = ad[ii] * (vx0-vxm1);

        vx[ii] = vx0;

        res2 += (r*r);
      }

    }
    else {

#     pragma omp parallel for reduction(+:res2) \
                          if(n_rows > CS_THR_MIN  && !_thread_debug)
      for (cs_lnum_t ll = 0; ll < n_rows; ll++) {

        cs_lnum_t ii = order[ll];

        const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
        const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
        const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

        cs_real_t vx0[DB_SIZE_MAX], vxm1[DB_SIZE_MAX], _vx[DB_SIZE_MAX];

        for (cs_lnum_t kk = 0; kk < db_size; kk++) {
          vxm1[kk] = vx[ii*db_size + kk];
          vx0[kk] = rhs[ii*db_size + kk];
        }

        for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
          for (cs_lnum_t kk = 0; kk < db_size; kk++)
            vx0[kk] -= (m_row[jj]*vx[col_id[jj]*db_size + kk]);
        }

        _fw_and_bw_lu_gs(ad_inv + db_size_2*ii,
                         db_size,
                         _vx,
                         vx0);

        double rr = 0;
        for (cs_lnum_t kk = 0; kk < db_size; kk++) {
          register double r = ad[ii*db_size + kk] * (_vx[kk]-vxm1[kk]);
          rr += (r*r);
          vx[ii*db_size + kk] = _vx[kk];
        }
        res2 += rr;

      }

    }

#if defined(HAVE_MPI)

    if (c->comm != MPI_COMM_NULL) {
      double _sum;
      MPI_Allreduce(&res2, &_sum, 1, MPI_DOUBLE, MPI_SUM, c->comm);
      res2 = _sum;
    }

#endif /* defined(HAVE_MPI) */

    residual = sqrt(res2); /* Actually, residual of previous iteration */

    /* Convergence test */

    if (n_iter == 1)
      c->setup_data->initial_residual = residual;

    cvg = _convergence_test(c, n_iter, residual, convergence);

  }

  return cvg;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using Process-local Gauss-Seidel.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- diagonal block size
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (unused here)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_p_gauss_seidel_msr(cs_sles_it_t              *c,
                    const cs_matrix_t         *a,
                    cs_lnum_t                  diag_block_size,
                    cs_sles_it_convergence_t  *convergence,
                    const cs_real_t           *rhs,
                    cs_real_t                 *restrict vx)
{
  cs_sles_convergence_state_t cvg;
  double  res2, residual;

  unsigned n_iter = 0;

  const cs_lnum_t n_rows = cs_matrix_get_n_rows(a);

  const cs_halo_t *halo = cs_matrix_get_halo(a);

  const cs_real_t  *restrict ad_inv = c->setup_data->ad_inv;

  const cs_real_t  *restrict ad = cs_matrix_get_diagonal(a);

  const cs_lnum_t  *a_row_index, *a_col_id;
  const cs_real_t  *a_d_val, *a_x_val;

  const cs_lnum_t db_size = cs_matrix_get_diag_block_size(a);
  const cs_lnum_t db_size_2 = db_size * db_size;

  cs_matrix_get_msr_arrays(a, &a_row_index, &a_col_id, &a_d_val, &a_x_val);

  cvg = CS_SLES_ITERATING;

  /* Current iteration */
  /*-------------------*/

  while (cvg == CS_SLES_ITERATING) {

    n_iter += 1;

    /* Synchronize ghost cells first */

    if (halo != NULL)
      cs_matrix_pre_vector_multiply_sync(a, vx);

    /* Compute Vx <- Vx - (A-diag).Rk and residual. */

    res2 = 0.0;

    if (diag_block_size == 1) {

#     pragma omp parallel for reduction(+:res2)      \
                          if(n_rows > CS_THR_MIN && !_thread_debug)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

        const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
        const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
        const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

        cs_real_t vxm1 = vx[ii];
        cs_real_t vx0 = rhs[ii];

        for (cs_lnum_t jj = 0; jj < n_cols; jj++)
          vx0 -= (m_row[jj]*vx[col_id[jj]]);

        vx0 *= ad_inv[ii];

        register double r = ad[ii] * (vx0-vxm1);
        res2 += (r*r);

        vx[ii] = vx0;
      }

    }
    else {

#     pragma omp parallel for reduction(+:res2) \
                          if(n_rows > CS_THR_MIN && !_thread_debug)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

        const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
        const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
        const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

        cs_real_t vx0[DB_SIZE_MAX], vxm1[DB_SIZE_MAX], _vx[DB_SIZE_MAX];

        for (cs_lnum_t kk = 0; kk < db_size; kk++) {
          vxm1[kk] = vx[ii*db_size + kk];
          vx0[kk] = rhs[ii*db_size + kk];
        }

        for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
          for (cs_lnum_t kk = 0; kk < db_size; kk++)
            vx0[kk] -= (m_row[jj]*vx[col_id[jj]*db_size + kk]);
        }

        _fw_and_bw_lu_gs(ad_inv + db_size_2*ii,
                         db_size,
                         _vx,
                         vx0);

        double rr = 0;
        for (cs_lnum_t kk = 0; kk < db_size; kk++) {
          register double r = ad[ii*db_size + kk] * (_vx[kk]-vxm1[kk]);
          rr += (r*r);
          vx[ii*db_size + kk] = _vx[kk];
        }
        res2 += rr;

      }

    }

    if (convergence->precision > 0. || c->plot != NULL) {

#if defined(HAVE_MPI)

      if (c->comm != MPI_COMM_NULL) {
        double _sum;
        MPI_Allreduce(&res2, &_sum, 1, MPI_DOUBLE, MPI_SUM, c->comm);
        res2 = _sum;
      }

#endif /* defined(HAVE_MPI) */

      residual = sqrt(res2); /* Actually, residual of previous iteration */

      /* Convergence test */

      if (n_iter == 1)
        c->setup_data->initial_residual = residual;

      cvg = _convergence_test(c, n_iter, residual, convergence);

    }
    else if (n_iter >= convergence->n_iterations_max) {
      convergence->n_iterations = n_iter;
      cvg = CS_SLES_MAX_ITERATION;
    }

  }

  return cvg;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using Process-local symmetric Gauss-Seidel.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- diagonal block size
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (unused here)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_p_sym_gauss_seidel_msr(cs_sles_it_t              *c,
                        const cs_matrix_t         *a,
                        cs_lnum_t                  diag_block_size,
                        cs_sles_it_convergence_t  *convergence,
                        const cs_real_t           *rhs,
                        cs_real_t                 *restrict vx,
                        size_t                     aux_size,
                        void                      *aux_vectors)
{
  CS_UNUSED(aux_size);
  CS_UNUSED(aux_vectors);

  cs_sles_convergence_state_t cvg;
  double  res2, residual;

  /* Check matrix storage type */

  if (cs_matrix_get_type(a) != CS_MATRIX_MSR)
    bft_error
      (__FILE__, __LINE__, 0,
       _("Symmetric Gauss-Seidel Jacobi hybrid solver only supported with a\n"
         "matrix using %s storage."),
       "MSR");

  unsigned n_iter = 0;

  const cs_lnum_t n_rows = cs_matrix_get_n_rows(a);

  const cs_halo_t *halo = cs_matrix_get_halo(a);

  const cs_real_t  *restrict ad_inv = c->setup_data->ad_inv;

  const cs_real_t  *restrict ad = cs_matrix_get_diagonal(a);

  const cs_lnum_t  *a_row_index, *a_col_id;
  const cs_real_t  *a_d_val, *a_x_val;

  const cs_lnum_t db_size = cs_matrix_get_diag_block_size(a);
  const cs_lnum_t db_size_2 = db_size * db_size;

  cs_matrix_get_msr_arrays(a, &a_row_index, &a_col_id, &a_d_val, &a_x_val);

  cvg = CS_SLES_ITERATING;

  /* Current iteration */
  /*-------------------*/

  while (cvg == CS_SLES_ITERATING) {

    n_iter += 1;

    /* Synchronize ghost cells first */

    if (halo != NULL)
      cs_matrix_pre_vector_multiply_sync(a, vx);

    /* Compute Vx <- Vx - (A-diag).Rk and residual: forward step */

    if (diag_block_size == 1) {

#     pragma omp parallel for if(n_rows > CS_THR_MIN && !_thread_debug)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

        const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
        const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
        const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

        cs_real_t vx0 = rhs[ii];

        for (cs_lnum_t jj = 0; jj < n_cols; jj++)
          vx0 -= (m_row[jj]*vx[col_id[jj]]);

        vx[ii] = vx0 * ad_inv[ii];

      }

    }
    else {

#     pragma omp parallel for if(n_rows > CS_THR_MIN && !_thread_debug)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

        const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
        const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
        const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

        cs_real_t vx0[DB_SIZE_MAX], _vx[DB_SIZE_MAX];

        for (cs_lnum_t kk = 0; kk < diag_block_size; kk++)
          vx0[kk] = rhs[ii*db_size + kk];

        for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
          for (cs_lnum_t kk = 0; kk < diag_block_size; kk++)
            vx0[kk] -= (m_row[jj]*vx[col_id[jj]*db_size + kk]);
        }

        _fw_and_bw_lu_gs(ad_inv + db_size_2*ii,
                         db_size,
                         _vx,
                         vx0);

        for (cs_lnum_t kk = 0; kk < diag_block_size; kk++)
          vx[ii*db_size + kk] = _vx[kk];

      }

    }

    /* Synchronize ghost cells again */

    if (halo != NULL)
      cs_matrix_pre_vector_multiply_sync(a, vx);

    /* Compute Vx <- Vx - (A-diag).Rk and residual: backward step */

    res2 = 0.0;

    if (diag_block_size == 1) {

#     pragma omp parallel for reduction(+:res2)      \
                          if(n_rows > CS_THR_MIN && !_thread_debug)
      for (cs_lnum_t ii = n_rows - 1; ii > - 1; ii--) {

        const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
        const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
        const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

        cs_real_t vxm1 = vx[ii];
        cs_real_t vx0 = rhs[ii];

        for (cs_lnum_t jj = 0; jj < n_cols; jj++)
          vx0 -= (m_row[jj]*vx[col_id[jj]]);

        vx0 *= ad_inv[ii];

        register double r = ad[ii] * (vx0-vxm1);
        res2 += (r*r);

        vx[ii] = vx0;
      }

    }
    else {

#     pragma omp parallel for reduction(+:res2) \
                          if(n_rows > CS_THR_MIN && !_thread_debug)
      for (cs_lnum_t ii = n_rows - 1; ii > - 1; ii--) {

        const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
        const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
        const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

        cs_real_t vx0[DB_SIZE_MAX], vxm1[DB_SIZE_MAX], _vx[DB_SIZE_MAX];

        for (cs_lnum_t kk = 0; kk < db_size; kk++) {
          vxm1[kk] = vx[ii*db_size + kk];
          vx0[kk] = rhs[ii*db_size + kk];
        }

        for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
          for (cs_lnum_t kk = 0; kk < db_size; kk++)
            vx0[kk] -= (m_row[jj]*vx[col_id[jj]*db_size + kk]);
        }

        _fw_and_bw_lu_gs(ad_inv + db_size_2*ii,
                         db_size,
                         _vx,
                         vx0);

        double rr = 0;
        for (cs_lnum_t kk = 0; kk < db_size; kk++) {
          register double r = ad[ii*db_size + kk] * (_vx[kk]-vxm1[kk]);
          rr += (r*r);
          vx[ii*db_size + kk] = _vx[kk];
        }
        res2 += rr;

      }

    }

    if (convergence->precision > 0. || c->plot != NULL) {

#if defined(HAVE_MPI)

      if (c->comm != MPI_COMM_NULL) {
        double _sum;
        MPI_Allreduce(&res2, &_sum, 1, MPI_DOUBLE, MPI_SUM, c->comm);
        res2 = _sum;
      }

#endif /* defined(HAVE_MPI) */

      residual = sqrt(res2); /* Actually, residual of previous iteration */

      /* Convergence test */

      if (n_iter == 1)
        c->setup_data->initial_residual = residual;

      cvg = _convergence_test(c, n_iter, residual, convergence);

    }
    else if (n_iter >= convergence->n_iterations_max) {
      convergence->n_iterations = n_iter;
      cvg = CS_SLES_MAX_ITERATION;
    }

  }

  return cvg;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using Process-local symmetric Gauss-Seidel.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- diagonal block size
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (unused here)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_p_gauss_seidel(cs_sles_it_t              *c,
                const cs_matrix_t         *a,
                cs_lnum_t                  diag_block_size,
                cs_sles_it_convergence_t  *convergence,
                const cs_real_t           *rhs,
                cs_real_t                 *restrict vx,
                size_t                     aux_size,
                void                      *aux_vectors)
{
  CS_UNUSED(aux_size);
  CS_UNUSED(aux_vectors);

  cs_sles_convergence_state_t cvg;

  /* Check matrix storage type */

  if (cs_matrix_get_type(a) != CS_MATRIX_MSR)
    bft_error
      (__FILE__, __LINE__, 0,
       _("Gauss-Seidel Jacobi hybrid solver only supported with a\n"
         "matrix using %s storage."),
       "MSR");

  /* Allocate or map work arrays */
  /*-----------------------------*/

  assert(c->setup_data != NULL);

  /* Check for ordered variant */

  const cs_lnum_t  *order = NULL;

  if (c->add_data != NULL)
    order = c->add_data->order;

  if (order != NULL)
    cvg = _p_ordered_gauss_seidel_msr(c,
                                      a,
                                      diag_block_size,
                                      convergence,
                                      rhs,
                                      vx);

  else
    cvg = _p_gauss_seidel_msr(c,
                              a,
                              diag_block_size,
                              convergence,
                              rhs,
                              vx);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Switch to fallback solver if defined.
 *
 * vx is reset to zero in this case.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- matrix
 *   solver_type     <-- fallback solver type
 *   prev_state      <-- previous convergence state
 *   convergence     <-> convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area
 *                       (internal allocation if NULL)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_fallback(cs_sles_it_t                    *c,
          cs_sles_it_type_t                solver_type,
          const cs_matrix_t               *a,
          cs_sles_convergence_state_t      prev_state,
          const cs_sles_it_convergence_t  *convergence,
          int                             *n_iter,
          double                          *residual,
          const cs_real_t                 *rhs,
          cs_real_t                       *restrict vx,
          size_t                           aux_size,
          void                            *aux_vectors)
{
  cs_sles_convergence_state_t cvg = CS_SLES_ITERATING;

  /* Check if fallback was already defined for this case */

  if (c->fallback == NULL) {

    /* Share context if possible */

    if (_setup_category(c->type) == _setup_category(solver_type)) {
      c->fallback = cs_sles_it_create(solver_type,
                                      -1, /* poly_degree */
                                      c->fallback_n_max_iter,
                                      c->update_stats);
      cs_sles_it_set_shareable(c->fallback, c);
    }
    else
      c->fallback = cs_sles_it_create(solver_type,
                                      0, /* poly_degree */
                                      c->fallback_n_max_iter,
                                      c->update_stats);

    c->fallback->plot = c->plot;

  }

  c->fallback->plot_time_stamp = c->plot_time_stamp;

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  if (prev_state < CS_SLES_BREAKDOWN) {
#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_rows; i++)
      vx[i] = 0;
  }

  cvg = cs_sles_it_solve(c->fallback,
                         convergence->name,
                         a,
                         convergence->verbosity,
                         convergence->precision,
                         convergence->r_norm,
                         n_iter,
                         residual,
                         rhs,
                         vx,
                         aux_size,
                         aux_vectors);

  cs_sles_it_free(c->fallback);

  *n_iter += convergence->n_iterations;

  c->plot_time_stamp = c->fallback->plot_time_stamp;

  return cvg;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

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
                       void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg = CS_SLES_CONVERGED;

  CS_UNUSED(c);
  CS_UNUSED(a);
  CS_UNUSED(diag_block_size);
  CS_UNUSED(convergence);
  CS_UNUSED(rhs);
  CS_UNUSED(vx);
  CS_UNUSED(aux_size);
  CS_UNUSED(aux_vectors);

  return cvg;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define and associate an iterative sparse linear system solver
 *        for a given field or equation name.
 *
 * If this system did not previously exist, it is added to the list of
 * "known" systems. Otherwise, its definition is replaced by the one
 * defined here.
 *
 * This is a utility function: if finer control is needed, see
 * \ref cs_sles_define and \ref cs_sles_it_create.
 *
 * Note that this function returns a pointer directly to the iterative solver
 * management structure. This may be used to set further options,
 * for example using \ref cs_sles_it_set_plot_options. If needed,
 * \ref cs_sles_find may be used to obtain a pointer to the matching
 * \ref cs_sles_t container.
 *
 * \param[in]  f_id          associated field id, or < 0
 * \param[in]  name          associated name if f_id < 0, or NULL
 * \param[in]  solver_type   type of solver (PCG, Jacobi, ...)
 * \param[in]  poly_degree   preconditioning polynomial degree
 *                           (0: diagonal; -1: non-preconditioned)
 * \param[in]  n_max_iter    maximum number of iterations
 *
 * \return  pointer to newly created iterative solver info object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_it_t *
cs_sles_it_define(int                 f_id,
                  const char         *name,
                  cs_sles_it_type_t   solver_type,
                  int                 poly_degree,
                  int                 n_max_iter)
{
  /* Test for environment variables here */

  const char *s = getenv("CS_THREAD_DEBUG");
  if (s != NULL) {
    if (atoi(s) > 0)
      _thread_debug = true;
  }

  /* Now define solver */

  cs_sles_it_t *
    c = cs_sles_it_create(solver_type,
                          poly_degree,
                          n_max_iter,
                          true); /* update stats */

  cs_sles_t *sc = cs_sles_define(f_id,
                                 name,
                                 c,
                                 "cs_sles_it_t",
                                 cs_sles_it_setup,
                                 cs_sles_it_solve,
                                 cs_sles_it_free,
                                 cs_sles_it_log,
                                 cs_sles_it_copy,
                                 cs_sles_it_destroy);

  cs_sles_set_error_handler(sc,
                            cs_sles_it_error_post_and_abort);

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create iterative sparse linear system solver info and context.
 *
 * \param[in]  solver_type   type of solver (PCG, Jacobi, ...)
 * \param[in]  poly_degree   preconditioning polynomial degree
 *                           (0: diagonal; -1: non-preconditioned;
 *                           see \ref sles_it for details)
 * \param[in]  n_max_iter    maximum number of iterations
 * \param[in]  update_stats  automatic solver statistics indicator
 *
 * \return  pointer to newly created solver info object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_it_t *
cs_sles_it_create(cs_sles_it_type_t   solver_type,
                  int                 poly_degree,
                  int                 n_max_iter,
                  bool                update_stats)
{
  cs_sles_it_t *c;

  BFT_MALLOC(c, 1, cs_sles_it_t);

  c->type = solver_type;
  c->solve = NULL;

  switch(c->type) {
  case CS_SLES_JACOBI:
  case CS_SLES_P_GAUSS_SEIDEL:
  case CS_SLES_P_SYM_GAUSS_SEIDEL:
    c->_pc = NULL;
    break;
  default:
    if (poly_degree < 0) {
       /* specific implementation for non-preconditioned PCG */
      if (c->type == CS_SLES_PCG)
        c->_pc = NULL;
      else
        c->_pc = cs_sles_pc_none_create();
    }
    else if (poly_degree == 0)
      c->_pc = cs_sles_pc_jacobi_create();
    else if (poly_degree == 1)
      c->_pc =cs_sles_pc_poly_1_create();
    else
      c->_pc =cs_sles_pc_poly_2_create();
  }
  c->pc = c->_pc;

  c->on_device = false;
  c->update_stats = update_stats;
  c->ignore_convergence = false;

  c->n_max_iter = n_max_iter;
  c->restart_interval = 20; /* Default value commonly found in the literature */

  c->n_setups = 0;
  c->n_solves = 0;

  c->n_iterations_min = 0;
  c->n_iterations_max = 0;
  c->n_iterations_last = 0;
  c->n_iterations_tot = 0;

  CS_TIMER_COUNTER_INIT(c->t_setup);
  CS_TIMER_COUNTER_INIT(c->t_solve);

  c->plot_time_stamp = 0;
  c->plot = NULL;
  c->_plot = NULL;

#if defined(HAVE_MPI)
  c->comm = cs_glob_mpi_comm;
  c->caller_comm = cs_glob_mpi_comm;
  c->caller_n_ranks = cs_glob_n_ranks;
  if (c->caller_n_ranks < 2) {
    c->comm = MPI_COMM_NULL;
    c->caller_comm = cs_glob_mpi_comm;
  }
#endif

  c->setup_data = NULL;
  c->add_data = NULL;
  c->shared = NULL;

  /* Fallback mechanism */

  switch(c->type) {
  case CS_SLES_BICGSTAB:
  case CS_SLES_BICGSTAB2:
  case CS_SLES_PCR3:
    c->fallback_cvg = CS_SLES_MAX_ITERATION;
    break;
  default:
    c->fallback_cvg = CS_SLES_DIVERGED;
  }
  c->fallback_n_max_iter = n_max_iter;

  c->fallback = NULL;

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy iterative sparse linear system solver info and context.
 *
 * \param[in, out]  context  pointer to iterative solver info and context
 *                           (actual type: cs_sles_it_t  **)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_destroy(void **context)
{
  cs_sles_it_t *c = (cs_sles_it_t *)(*context);
  if (c != NULL) {
    if (c->fallback != NULL) {
      void *f = c->fallback;
      cs_sles_it_destroy(&f);
      c->fallback = f;
    }
    cs_sles_pc_destroy(&(c->_pc));
    cs_sles_it_free(c);
    if (c->_plot != NULL) {
      cs_time_plot_finalize(&(c->_plot));
      c->plot = NULL;
    }
    if (c->add_data != NULL) {
      BFT_FREE(c->add_data->order);
      BFT_FREE(c->add_data);
    }
    BFT_FREE(c);
    *context = c;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create iterative sparse linear system solver info and context
 *        based on existing info and context.
 *
 * \param[in]  context  pointer to reference info and context
 *                     (actual type: cs_sles_it_t  *)
 *
 * \return  pointer to newly created solver info object.
 *          (actual type: cs_sles_it_t  *)
 */
/*----------------------------------------------------------------------------*/

void *
cs_sles_it_copy(const void  *context)
{
  cs_sles_it_t *d = NULL;

  if (context != NULL) {
    const cs_sles_it_t *c = context;
    d = cs_sles_it_create(c->type,
                          -1,
                          c->n_max_iter,
                          c->update_stats);
    if (c->pc != NULL && c->_pc != NULL) {
      d->_pc = cs_sles_pc_clone(c->_pc);
      d->pc = d->_pc;
    }
    else {
      d->_pc = c->_pc;
      d->pc = c->pc;
    }

    /* If useful, copy the restart interval */
    if (c->type == CS_SLES_GMRES || c->type == CS_SLES_GCR)
      d->restart_interval = c->restart_interval;

#if defined(HAVE_MPI)
    d->comm = c->comm;
#endif
  }

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log sparse linear equation solver info.
 *
 * \param[in]  context   pointer to iterative solver info and context
 *                       (actual type: cs_sles_it_t  *)
 * \param[in]  log_type  log type
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_log(const void  *context,
               cs_log_t     log_type)
{
  const cs_sles_it_t  *c = context;

  if (log_type == CS_LOG_SETUP) {

    cs_log_printf(log_type,
                  _("  Solver type:                       %s\n"),
                  _(cs_sles_it_type_name[c->type]));
    if (c->pc != NULL)
      cs_log_printf(log_type,
                    _("  Preconditioning:                   %s\n"),
                    _(cs_sles_pc_get_type_name(c->pc)));
    if (c->type == CS_SLES_GMRES || c->type == CS_SLES_GCR)
      cs_log_printf(log_type,
                    "  Restart interval:                  %d\n",
                    c->restart_interval);
    cs_log_printf(log_type,
                  _("  Maximum number of iterations:      %d\n"),
                  c->n_max_iter);
  }

  else if (log_type == CS_LOG_PERFORMANCE) {

    int n_calls = c->n_solves;
    int n_it_min = c->n_iterations_min;
    int n_it_max = c->n_iterations_max;
    int n_it_mean = 0;

    if (n_it_min < 0)
      n_it_min = 0;

    if (n_calls > 0)
      n_it_mean = (int)(  c->n_iterations_tot
                         / ((unsigned long long)n_calls));

    if (n_it_mean == 0)
      cs_log_printf(log_type, _("\n  No resolution\n"));

    else {

      cs_log_printf(log_type,
                    _("\n"
                      "  Solver type:                   %s\n"),
                    _(cs_sles_it_type_name[c->type]));

      if (c->pc != NULL)
        cs_log_printf(log_type,
                      _("  Preconditioning:               %s\n"),
                      _(cs_sles_pc_get_type_name(c->pc)));
      cs_log_printf(log_type,
                    _("  Number of setups:              %12d\n"
                      "  Number of calls:               %12d\n"
                      "  Minimum number of iterations:  %12d\n"
                      "  Maximum number of iterations:  %12d\n"
                      "  Mean number of iterations:     %12d\n"
                      "  Total setup time:              %12.3f\n"
                      "  Total solution time:           %12.3f\n"),
                    c->n_setups, n_calls, n_it_min, n_it_max, n_it_mean,
                    c->t_setup.nsec*1e-9,
                    c->t_solve.nsec*1e-9);

      if (c->fallback != NULL) {

        n_calls = c->fallback->n_solves;
        n_it_min = c->fallback->n_iterations_min;
        n_it_max = c->fallback->n_iterations_max;
        n_it_mean = 0;

        if (n_it_min < 0)
          n_it_min = 0;

        if (n_calls > 0)
          n_it_mean = (int)(  c->fallback->n_iterations_tot
                              / ((unsigned long long)n_calls));

        cs_log_printf(log_type,
                      _("\n"
                        "  Backup solver type:            %s\n"),
                      _(cs_sles_it_type_name[c->fallback->type]));

        cs_log_printf(log_type,
                      _("  Number of calls:               %12d\n"
                        "  Minimum number of iterations:  %12d\n"
                        "  Maximum number of iterations:  %12d\n"
                        "  Mean number of iterations:     %12d\n"
                        "  Total solution time:           %12.3f\n"),
                      n_calls, n_it_min, n_it_max, n_it_mean,
                      c->fallback->t_solve.nsec*1e-9);

      } /* If used */

    }

  }

  if (c->pc != NULL)
    cs_sles_pc_log(c->pc, log_type);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup iterative sparse linear equation solver.
 *
 * \param[in, out]  context    pointer to iterative solver info and context
 *                             (actual type: cs_sles_it_t  *)
 * \param[in]       name       pointer to system name
 * \param[in]       a          associated matrix
 * \param[in]       verbosity  associated verbosity
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_setup(void               *context,
                 const char         *name,
                 const cs_matrix_t  *a,
                 int                 verbosity)
{
  cs_sles_it_t  *c = context;

  cs_timer_t t0;
  if (c->update_stats == true)
    t0 = cs_timer_time();

#if defined(HAVE_ACCEL)
  bool on_device = (cs_matrix_get_alloc_mode(a) > CS_ALLOC_HOST) ?
    true : false;
#endif

  const int diag_block_size = (cs_matrix_get_diag_block_size(a));

  if (verbosity > 1) {
    bft_printf(_("\n Setup of solver for linear system \"%s\"\n"),
               name);
    cs_matrix_log_info(a, verbosity);
  }

  bool block_nn_inverse = false;

  if (   c->type == CS_SLES_JACOBI
      || (   c->type >= CS_SLES_P_GAUSS_SEIDEL
          && c->type <= CS_SLES_P_SYM_GAUSS_SEIDEL)) {
    /* Force to Jacobi in case matrix type is not adapted */
    if (cs_matrix_get_type(a) != CS_MATRIX_MSR) {
      c->type = CS_SLES_JACOBI;
    }
    block_nn_inverse = true;
  }

  switch (c->type) {

  case CS_SLES_PCR3:
    c->solve = _conjugate_residual_3;
    break;

  case CS_SLES_PCG:
    /* Check for single-reduction */
    {
      bool single_reduce = false;
#if defined(HAVE_MPI)
      cs_gnum_t n_m_rows = cs_matrix_get_n_rows(a) * diag_block_size;
      if (c->comm != MPI_COMM_NULL) {
        int size;
        cs_gnum_t _n_m_rows;
        MPI_Allreduce(&n_m_rows, &_n_m_rows, 1, CS_MPI_GNUM, MPI_SUM, c->comm);
        MPI_Comm_size(c->comm, &size);
        n_m_rows = _n_m_rows / (cs_gnum_t)size;
      }
      if (c->comm != c->caller_comm)
        MPI_Bcast(&n_m_rows, 1, CS_MPI_GNUM, 0, cs_glob_mpi_comm);
      if (n_m_rows < (cs_gnum_t)_pcg_sr_threshold)
        single_reduce = true;
#endif
      if (!single_reduce) {
        if (c->pc != NULL)
          c->solve = _conjugate_gradient;
        else
          c->solve = _conjugate_gradient_npc;
        break;
      }
      else {
        if (c->pc != NULL)
          c->solve = _conjugate_gradient_sr;
        else
          c->solve = _conjugate_gradient_npc_sr;
      }
    }
    break;

  case CS_SLES_FCG:
    c->solve = _flexible_conjugate_gradient;
#if defined(HAVE_CUDA)
    if (on_device) {
      c->on_device = true;
      c->solve = cs_sles_it_cuda_fcg;
    }
#endif
    break;

  case CS_SLES_IPCG:
    c->solve = _conjugate_gradient_ip;
    break;

  case CS_SLES_JACOBI:
    if (diag_block_size == 1)
      c->solve = _jacobi;
    else if (diag_block_size == 3)
      c->solve = _block_3_jacobi;
    else
      c->solve = _block_jacobi;

#if defined(HAVE_CUDA)
    if (on_device) {
      c->on_device = true;
      if (diag_block_size == 1)
        c->solve = cs_sles_it_cuda_jacobi;
      else
        c->solve = cs_sles_it_cuda_block_jacobi;
    }
#endif

    break;

  case CS_SLES_BICGSTAB:
    c->solve = _bi_cgstab;
    break;

  case CS_SLES_BICGSTAB2:
    c->solve = _bicgstab2;
    break;

  case CS_SLES_GCR:
    assert(c->restart_interval > 1);
    c->solve = _gcr;
#if defined(HAVE_CUDA)
    if (on_device) {
      c->on_device = true;
      c->solve = cs_sles_it_cuda_gcr;
    }
#endif

    break;

  case CS_SLES_GMRES:
    assert(c->restart_interval > 1);
    c->solve = _gmres;
    break;

  case CS_SLES_P_GAUSS_SEIDEL:
    c->solve = _p_gauss_seidel;
    break;
  case CS_SLES_P_SYM_GAUSS_SEIDEL:
    c->solve = _p_sym_gauss_seidel_msr;
    break;

  case CS_SLES_USER_DEFINED:
    c->solve = cs_user_sles_it_solver;
    break;

  default:
    bft_error
      (__FILE__, __LINE__, 0,
       _("Setup of linear equation on \"%s\"\n"
         "with solver type %d, which is not defined)."),
       name, (int)c->type);
    break;
  }

  /* Setup preconditioner and/or auxiliary data */

  cs_sles_it_setup_priv(c, name, a, verbosity, diag_block_size,
                        block_nn_inverse);

  /* Now finish */

  if (c->update_stats == true) {
    cs_timer_t t1 = cs_timer_time();
    c->n_setups += 1;
    cs_timer_counter_add_diff(&(c->t_setup), &t0, &t1);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Call iterative sparse linear equation solver.
 *
 * \param[in, out]  context        pointer to iterative solver info and context
 *                                 (actual type: cs_sles_it_t  *)
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
                 void                *aux_vectors)
{
  cs_sles_it_t  *c = context;

  cs_sles_convergence_state_t cvg = CS_SLES_ITERATING;

  cs_timer_t t0 = {0, 0}, t1;

  unsigned _n_iter = 0;
  cs_sles_it_convergence_t  convergence;

  if (c->update_stats == true)
    t0 = cs_timer_time();

  const cs_lnum_t diag_block_size = cs_matrix_get_diag_block_size(a);

  /* Initialize number of iterations and residual,
     check for immediate return,
     and solve sparse linear system */

  *n_iter = 0;

  /* Setup if not already done */

  if (c->setup_data == NULL) {

    if (c->update_stats) { /* Stop solve timer to switch to setup timer */
      t1 = cs_timer_time();
      cs_timer_counter_add_diff(&(c->t_solve), &t0, &t1);
    }

    cs_sles_it_setup(c, name, a, verbosity);

    if (c->update_stats) /* Restart solve timer */
      t0 = cs_timer_time();

  }

  if (c->pc != NULL)
    cs_sles_pc_set_tolerance(c->pc, precision, r_norm);

  /* Solve sparse linear system */

  cs_sles_it_convergence_init(&convergence,
                              name,
                              verbosity,
                              c->n_max_iter,
                              precision,
                              r_norm,
                              residual);

  c->setup_data->initial_residual = -1;

  if (verbosity > 1)
    cs_log_printf(CS_LOG_DEFAULT,
                  _(" RHS norm:          %15.9e\n\n"), r_norm);

  /* Only call solver for "active" ranks */

  bool local_solve = true;
#if defined(HAVE_MPI)
  if (c->comm == MPI_COMM_NULL) {
    cs_lnum_t n_rows = cs_matrix_get_n_rows(a);
    if (n_rows == 0) {
      local_solve = false;
      cvg = CS_SLES_CONVERGED;
    }
  }
#endif

  if (local_solve) {

    const cs_real_t *_rhs = rhs;
    cs_real_t       *_vx = vx;

#if defined(HAVE_ACCEL)

    /* When using accelerators, ideally, the arrays passed here should
       already be usable on the device, but in case they are allocated with
       Fortran intrinsics, we copy them here. When using unified shared
       memory, we let the called solver handle prefetching rather than
       doing it here, as it is not strictly needed, and the solver may
       exploit additional asynchronicity */

    cs_lnum_t v_size = diag_block_size * cs_matrix_get_n_columns(a);
    cs_real_t *_rhs_w = NULL;
    cs_alloc_mode_t amode_vx = CS_ALLOC_HOST, amode_rhs = CS_ALLOC_HOST;

    if (c->on_device) {
      amode_vx = cs_check_device_ptr(vx);
      amode_rhs = cs_check_device_ptr(rhs);

      if (amode_vx == CS_ALLOC_HOST) {
        CS_MALLOC_HD(_vx, v_size, cs_real_t, CS_ALLOC_DEVICE);
        cs_copy_h2d(_vx, vx, v_size*sizeof(cs_real_t));
      }
      else if (amode_vx < CS_ALLOC_HOST_DEVICE_SHARED)
        cs_sync_h2d(_vx);

      if (amode_rhs == CS_ALLOC_HOST) {
        CS_MALLOC_HD(_rhs_w, v_size, cs_real_t, CS_ALLOC_DEVICE);
        cs_copy_h2d(_rhs_w, rhs, v_size*sizeof(cs_real_t));
        _rhs = _rhs_w;
      }
      else if (amode_rhs < CS_ALLOC_HOST_DEVICE_SHARED)
        cs_sync_h2d(_rhs);
    }

#endif

    /* Call selected solver function */

    cvg = c->solve(c,
                   a, diag_block_size, &convergence,
                   _rhs, _vx,
                   aux_size, aux_vectors);

#if defined(HAVE_ACCEL)

    /* Free device copies if needed */

    if (c->on_device) {
      if (amode_vx == CS_ALLOC_HOST) {
        cs_copy_d2h(vx, _vx, v_size*sizeof(cs_real_t));
        CS_FREE_HD(_vx);
      }
      else if (amode_vx < CS_ALLOC_HOST_DEVICE_SHARED)
        cs_sync_d2h(_vx);

      CS_FREE_HD(_rhs_w);
    }

#endif

  }

  /* Broadcast convergence info from "active" ranks to others*/

#if defined(HAVE_MPI)
  if (c->comm != c->caller_comm && c->ignore_convergence == false) {
    /* cvg is signed, so shift (with some margin) before copy to unsigned. */
    unsigned buf[2] = {(unsigned)cvg+10, convergence.n_iterations};
    MPI_Bcast(buf, 2, MPI_UNSIGNED, 0, c->caller_comm);
    MPI_Bcast(&convergence.residual, 1, MPI_DOUBLE, 0, c->caller_comm);
    cvg = (cs_sles_convergence_state_t)(buf[0] - 10);
    convergence.n_iterations = buf[1];
  }
#endif

  /* Update return values */

  _n_iter = convergence.n_iterations;

  *n_iter = convergence.n_iterations;
  *residual = convergence.residual;

  cs_sles_it_type_t fallback_type = CS_SLES_N_IT_TYPES;
  if (cvg < c->fallback_cvg)
    fallback_type = (c->on_device) ? CS_SLES_GCR : CS_SLES_GMRES;

  if (c->update_stats == true) {

    t1 = cs_timer_time();

    c->n_solves += 1;

    if (c->n_iterations_tot == 0)
      c->n_iterations_min = _n_iter;
    else if (c->n_iterations_min > _n_iter)
      c->n_iterations_min = _n_iter;
    if (c->n_iterations_max < _n_iter)
      c->n_iterations_max = _n_iter;

    c->n_iterations_last = _n_iter;
    c->n_iterations_tot += _n_iter;

    cs_timer_counter_add_diff(&(c->t_solve), &t0, &t1);

  }

  if (fallback_type != CS_SLES_N_IT_TYPES)
    cvg = _fallback(c,
                    fallback_type,
                    a,
                    cvg,
                    &convergence,
                    n_iter,
                    residual,
                    rhs,
                    vx,
                    aux_size,
                    aux_vectors);

  return cvg;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free iterative sparse linear equation solver setup context.
 *
 * This function frees resolution-related data, such as
 * buffers and preconditioning but does not free the whole context,
 * as info used for logging (especially performance data) is maintained.
 *
 * \param[in, out]  context  pointer to iterative solver info and context
 *                           (actual type: cs_sles_it_t  *)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_free(void  *context)
{
  cs_sles_it_t  *c = context;

  cs_timer_t t0;
  if (c->update_stats == true)
    t0 = cs_timer_time();

  if (c->fallback != NULL)
    cs_sles_it_free(c->fallback);

  if (c->_pc != NULL)
    cs_sles_pc_free(c->_pc);

  if (c->setup_data != NULL) {
    CS_FREE_HD(c->setup_data->_ad_inv);
    BFT_FREE(c->setup_data);
  }

  if (c->update_stats == true) {
    cs_timer_t t1 = cs_timer_time();
    cs_timer_counter_add_diff(&(c->t_setup), &t0, &t1);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return iterative solver type.
 *
 * \param[in]  context  pointer to iterative solver info and context
 *
 * \return  selected solver type
 */
/*----------------------------------------------------------------------------*/

cs_sles_it_type_t
cs_sles_it_get_type(const cs_sles_it_t  *context)
{
  return context->type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the initial residual for the previous solve operation
 *        with a solver.
 *
 * This is useful for convergence tests when this solver is used as
 * a preconditioning smoother.
 *
 * This operation is only valid between calls to \ref cs_sles_it_setup
 * (or \ref cs_sles_it_solve) and \ref cs_sles_it_free.
 * It returns -1 otherwise.
 *
 * \param[in]  context  pointer to iterative solver info and context
 *
 * \return initial residual from last call to \ref cs_sles_solve with this
 *         solver
 */
/*----------------------------------------------------------------------------*/

double
cs_sles_it_get_last_initial_residual(const cs_sles_it_t  *context)
{
  double retval = 1;
  if (context->setup_data != NULL)
    retval = context->setup_data->initial_residual;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a preconditioner context for an iterative sparse linear
 *        equation solver.
 *
 * This allows modifying parameters of a non default (Jacobi or polynomial)
 * preconditioner.
 *
 * \param[in]  context   pointer to iterative solver info and context
 *
 * \return  pointer to preconditoner context
 */
/*----------------------------------------------------------------------------*/

cs_sles_pc_t  *
cs_sles_it_get_pc(cs_sles_it_t  *context)
{
  cs_sles_pc_t  *pc = NULL;

  if (context != NULL) {
    cs_sles_it_t  *c = context;
    pc = c->pc;
  }

  return pc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a preconditioner to an iterative sparse linear equation
 *        solver, transfering its ownership to to solver context.
 *
 * This allows assigning a non default (Jacobi or polynomial) preconditioner.
 *
 * The input pointer is set to NULL to make it clear the caller does not
 * own the preconditioner anymore, though the context can be accessed using
 * \ref cs_sles_it_get_pc.
 *
 * \param[in, out]  context   pointer to iterative solver info and context
 * \param[in, out]  pc        pointer to preconditoner context
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_transfer_pc(cs_sles_it_t     *context,
                       cs_sles_pc_t    **pc)
{
  if (context != NULL) {
    cs_sles_it_t  *c = context;
    c->pc = NULL;
    cs_sles_pc_destroy(&(c->_pc));
    if (pc != NULL) {
      c->_pc = *pc;
      c->pc = *pc;
    }
  }
  else if (pc != NULL)
    cs_sles_pc_destroy(pc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy options from one iterative sparse linear system solver info
 *        and context to another.
 *
 * Optional plotting contexts are shared between the source and destination
 * contexts.
 *
 * Preconditioner settings are to be handled separately.
 *
 * \param[in]       src   pointer to source info and context
 * \param[in, out]  dest  pointer to destination info and context
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_transfer_parameters(const cs_sles_it_t  *src,
                               cs_sles_it_t        *dest)
{
  if (dest != NULL && src != NULL) {

    dest->update_stats = src->update_stats;
    dest->n_max_iter = src->n_max_iter;
    dest->restart_interval = src->restart_interval;

    dest->plot_time_stamp = src->plot_time_stamp;
    dest->plot = src->plot;
    if (dest->_plot != NULL)
      cs_time_plot_finalize(&(dest->_plot));

#if defined(HAVE_MPI)
    dest->comm = src->comm;
#endif

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Associate a similar info and context object with which some setup
 *        data may be shared.
 *
 * This is especially useful for sharing preconditioning data between
 * similar solver contexts (for example ascending and descending multigrid
 * smoothers based on the same matrix).
 *
 * For preconditioning data to be effectively shared, \ref cs_sles_it_setup
 * (or \ref cs_sles_it_solve) must be called on \p shareable before being
 * called on \p context (without \ref cs_sles_it_free being called in between,
 * of course).
 *
 * It is the caller's responsibility to ensure the context is not used
 * for a \ref cs_sles_it_setup or \ref cs_sles_it_solve operation after the
 * shareable object has been destroyed (normally by \ref cs_sles_it_destroy).
 *
 * \param[in, out]  context    pointer to iterative solver info and context
 * \param[in]       shareable  pointer to iterative solver info and context
 *                             whose context may be shared
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_set_shareable(cs_sles_it_t        *context,
                         const cs_sles_it_t  *shareable)
{
  cs_sles_it_t  *c = context;

  c->shared = shareable;

  c->pc = shareable->pc;

  if (c->pc != c->_pc && c->_pc != NULL)
    cs_sles_pc_destroy(&(c->_pc));
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set MPI communicator for global reductions.
 *
 * The system is solved only on ranks with a non-NULL communicator or
 * if the caller communicator has less than 2 ranks. convergence info
 * is the broadcast across the caller communicator.
 *
 * \param[in, out]  context      pointer to iterative solver info and context
 * \param[in]       comm         MPI communicator for solving
 * \param[in]       caller_comm  MPI communicator of caller
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_set_mpi_reduce_comm(cs_sles_it_t  *context,
                               MPI_Comm       comm,
                               MPI_Comm       caller_comm)
{
  cs_sles_it_t  *c = context;

  static int flag = -1;

  if (flag < 0)
    flag = cs_halo_get_use_barrier();

  c->comm = comm;
  c->caller_comm = caller_comm;

  if (c->caller_comm != MPI_COMM_NULL)
    MPI_Comm_size(c->caller_comm, &(c->caller_n_ranks));
  else
    c->caller_n_ranks = 1;

  if (comm != cs_glob_mpi_comm)
    cs_halo_set_use_barrier(0);
  else {
    cs_halo_set_use_barrier(flag);
    if (cs_glob_n_ranks < 2)
      c->comm = MPI_COMM_NULL;
  }
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign ordering to iterative solver.
 *
 * The solver context takes ownership of the order array (i.e. it will
 * handle its later deallocation).
 *
 * This is useful only for Process-local Gauss-Seidel.
 *
 * \param[in, out]  context  pointer to iterative solver info and context
 * \param[in, out]  order    pointer to ordering array
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_assign_order(cs_sles_it_t   *context,
                        cs_lnum_t     **order)
{
  if (context->type != CS_SLES_P_GAUSS_SEIDEL)
    BFT_FREE(*order);

  else {

    if (context->add_data == NULL) {
      BFT_MALLOC(context->add_data, 1, cs_sles_it_add_t);
      context->add_data->order = NULL;
    }

    BFT_FREE(context->add_data->order);

    context->add_data->order = *order;

    *order = NULL;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the threshold value under which a breakdown happens in
 *        solvers like BiCGStab or BiCGStab2
 *
 * \return the value of the threshold
 */
/*----------------------------------------------------------------------------*/

double
cs_sles_it_get_breakdown_threshold(void)
{
  return _epzero;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the threshold value under which a breakdown happens in
 *        solvers like BiCGStab or BiCGStab2
 *
 * \param[in]  threshold  value of the threshold
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_set_breakdown_threshold(double  threshold)
{
  _epzero = threshold;
}

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
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_set_fallback_threshold(cs_sles_it_t                 *context,
                                  cs_sles_convergence_state_t   threshold,
                                  int                           n_iter_max)
{
  context->fallback_cvg = threshold;
  context->fallback_n_max_iter = n_iter_max;
}

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
                                int            interval)
{
  if (context == NULL)
    return;

  context->restart_interval = interval;
}

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
                          int            n_max_iter)
{
  if (context == NULL)
    return;

  context->n_max_iter = n_max_iter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Query mean number of rows under which Conjugate Gradient algorithm
 *        uses the single-reduction variant.
 *
 * The single-reduction variant requires only one parallel sum per
 * iteration (instead of 2), at the cost of additional vector operations,
 * so it tends to be more expensive when the number of matrix rows per
 * MPI rank is high, then becomes cheaper when the MPI latency cost becomes
 * more significant.
 *
 * This option is ignored for non-parallel runs, so 0 is returned.
 *
 * \returns  mean number of rows per active rank under which the
 *           single-reduction variant will be used
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_sles_it_get_pcg_single_reduction(void)
{
#if defined(HAVE_MPI)
  return _pcg_sr_threshold;
#else
  return 0;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set mean number of rows under which Conjugate Gradient algorithm
 *        should use the single-reduction variant.
 *
 * The single-reduction variant requires only one parallel sum per
 * iteration (instead of 2), at the cost of additional vector operations,
 * so it tends to be more expensive when the number of matrix rows per
 * MPI rank is high, then becomes cheaper when the MPI latency cost becomes
 * more significant.
 *
 * This option is ignored for non-parallel runs.
 *
 * \param[in]  threshold  mean number of rows per active rank under which the
 *             single-reduction variant will be used
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_set_pcg_single_reduction(cs_lnum_t  threshold)
{
#if defined(HAVE_MPI)
  _pcg_sr_threshold = threshold;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the current global settings relative to parallelism.
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_log_parallel_options(void)
{
#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)
    cs_log_printf(CS_LOG_SETUP,
                  _("\n"
                    "Iterative linear solvers parallel parameters:\n"
                    "  PCG single-reduction threshold:     %ld\n"),
                  (long)_pcg_sr_threshold);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Error handler for iterative sparse linear equation solver.
 *
 * In case of divergence or breakdown, this error handler outputs
 * postprocessing data to assist debugging, then aborts the run.
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
cs_sles_it_error_post_and_abort(cs_sles_t                    *sles,
                                cs_sles_convergence_state_t   state,
                                const cs_matrix_t            *a,
                                const cs_real_t              *rhs,
                                cs_real_t                    *vx)
{
  if (state >= CS_SLES_BREAKDOWN)
    return false;

  const cs_sles_it_t  *c = cs_sles_get_context(sles);
  const char *name = cs_sles_get_name(sles);

  int mesh_id = cs_post_init_error_writer_cells();

  cs_sles_post_error_output_def(name,
                                mesh_id,
                                a,
                                rhs,
                                vx);

  cs_post_finalize();

  const char *error_type[] = {N_("divergence"), N_("breakdown")};
  int err_id = (state == CS_SLES_BREAKDOWN) ? 1 : 0;

  bft_error(__FILE__, __LINE__, 0,
            _("%s: error (%s) solving for %s"),
            _(cs_sles_it_type_name[c->type]),
            _(error_type[err_id]),
            name);

  return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set plotting options for an iterative sparse linear equation solver.
 *
 * \param[in, out]  context        pointer to iterative solver info and context
 * \param[in]       base_name      base plot name to activate, NULL otherwise
 * \param[in]       use_iteration  if true, use iteration as time stamp
 *                                 otherwise, use wall clock time
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_set_plot_options(cs_sles_it_t  *context,
                            const char    *base_name,
                            bool           use_iteration)
{
  if (context != NULL) {
    if (cs_glob_rank_id < 1 && base_name != NULL) {

      /* Destroy previous plot if options reset */
      if (context->_plot != NULL)
        cs_time_plot_finalize(&(context->_plot));

      /* Create new plot */
      cs_file_mkdir_default("monitoring");
      const char *probe_names[] = {base_name};
      context->_plot = cs_time_plot_init_probe(base_name,
                                               "monitoring/residual_",
                                               CS_TIME_PLOT_CSV,
                                               use_iteration,
                                               -1.0,  /* force flush */
                                               0,     /* no buffer */
                                               1,     /* n_probes */
                                               NULL,  /* probe_list */
                                               NULL,  /* probe_coords */
                                               probe_names);
      context->plot = context->_plot;
      context->plot_time_stamp = 0;
    }
  }
}

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
                            cs_sles_it_convergence_t  *convergence)
{
  return _convergence_test(c, n_iter, residual, convergence);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
