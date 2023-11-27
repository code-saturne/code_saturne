/*============================================================================
 * Sparse Linear Equation Solvers: Multigrid smoothers
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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_sles_it.h"
#include "cs_sles_it_priv.h"
#include "cs_multigrid_smoother.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_multigrid_smoother.c
        Iterative linear solvers used as multigrid smoothers only.

  These smoothers are based on iterative solvers, but are simplified so
  as to avoid the cost of some residual computation and convergence testing
  operations.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define RINFIN  1.E+30

#if !defined(HUGE_VAL)
#define HUGE_VAL  1.E+12
#endif

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

/*============================================================================
 * Private function definitions
 *============================================================================*/

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
  double  ro_0, ro_1, alpha, rk_gkm1, rk_gk, beta;
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

  {
    /* Preconditioning */

    c->setup_data->pc_apply(c->setup_data->pc_context, rk, gk);

    /* Descent direction */

    cs_matrix_vector_multiply(a, gk, zk);  /* dk = gk */

    /* Descent parameter */

    _dot_products_xy_yz(c, rk, gk, zk, &ro_0, &ro_1);

    rk_gkm1 = ro_0;
    alpha = -ro_0 / ro_1;

    if (convergence->n_iterations_max < 2) {
#     pragma omp parallel for if(n_rows > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        vx[ii] += (alpha * gk[ii]);
    }
    else {
#     pragma omp parallel if(n_rows > CS_THR_MIN)
      {
#       pragma omp for nowait
        for (cs_lnum_t ii = 0; ii < n_rows; ii++)
          dk[ii] = gk[ii];
#       pragma omp for nowait
        for (cs_lnum_t ii = 0; ii < n_rows; ii++)
          vx[ii] += (alpha * gk[ii]);
#       pragma omp for nowait
        for (cs_lnum_t ii = 0; ii < n_rows; ii++)
          rk[ii] += (alpha * zk[ii]);
      }
    }
  }

  for (n_iter = 1; n_iter < convergence->n_iterations_max; n_iter++) {

    /* Preconditioning */

    c->setup_data->pc_apply(c->setup_data->pc_context, rk, gk);

    /* Prepare descent parameter */

    rk_gk = _dot_product(c, rk, gk);

    /* Complete descent parameter computation and matrix.vector product */

    beta = rk_gk / rk_gkm1;
    rk_gkm1 = rk_gk;

#   pragma omp parallel for firstprivate(alpha) if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      dk[ii] = gk[ii] + (beta * dk[ii]);

    cs_matrix_vector_multiply(a, dk, zk);

    _dot_products_xy_yz(c, rk, dk, zk, &ro_0, &ro_1);

    alpha = -ro_0 / ro_1;

    if (n_iter + 1 < convergence->n_iterations_max) {
#     pragma omp parallel if(n_rows > CS_THR_MIN)
      {
#       pragma omp for nowait
        for (cs_lnum_t ii = 0; ii < n_rows; ii++)
          vx[ii] += (alpha * dk[ii]);
#       pragma omp for nowait
        for (cs_lnum_t ii = 0; ii < n_rows; ii++)
          rk[ii] += (alpha * zk[ii]);
      }
    }
    else { /* last iteration */
#     pragma omp parallel for if(n_rows > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        vx[ii] += alpha * dk[ii];
    }
  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  convergence->n_iterations = n_iter;

  return CS_SLES_MAX_ITERATION;
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
  double  ro_0, ro_1, alpha, rk_gkm1, rk_gk, gk_sk, beta;
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

  {

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      rk[ii] -= rhs[ii];

    /* Preconditionning */

    c->setup_data->pc_apply(c->setup_data->pc_context, rk, gk);

    /* Descent direction */

    cs_matrix_vector_multiply(a, gk, zk); /* zk = A.dk
                                                            dk == gk */

    /* Descent parameter */

    _dot_products_xy_yz(c, rk, gk, zk, &ro_0, &ro_1);

    alpha = -ro_0 / ro_1;
    rk_gkm1 = ro_0;

    if (convergence->n_iterations_max < 2) {
#     pragma omp parallel for if(n_rows > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        vx[ii] += (alpha * gk[ii]);
    }
    else {
#     pragma omp parallel if(n_rows > CS_THR_MIN)
      {
#       pragma omp for nowait
        for (cs_lnum_t ii = 0; ii < n_rows; ii++)
          dk[ii] = gk[ii];
#       pragma omp for nowait
        for (cs_lnum_t ii = 0; ii < n_rows; ii++)
          vx[ii] += (alpha * gk[ii]);
#       pragma omp for nowait
        for (cs_lnum_t ii = 0; ii < n_rows; ii++)
          rk[ii] += (alpha * zk[ii]);
      }
    }
  }

  /* Current Iteration */
  /*-------------------*/

  for (n_iter = 1; n_iter < convergence->n_iterations_max; n_iter++) {

    /* Preconditionning */

    c->setup_data->pc_apply(c->setup_data->pc_context, rk, gk);

    cs_matrix_vector_multiply(a, gk, sk);  /* sk = A.gk */

    /* compute residual and prepare descent parameter */

    _dot_products_xy_yz(c, rk, gk, sk, &rk_gk, &gk_sk);

    /* Complete descent parameter computation and matrix.vector product */

    beta = rk_gk / rk_gkm1;
    rk_gkm1 = rk_gk;

    ro_1 = gk_sk - beta*beta*ro_1;
    ro_0 = rk_gk;

    alpha = -ro_0 / ro_1;

    if (n_iter + 1 < convergence->n_iterations_max) {
#     pragma omp parallel if(n_rows > CS_THR_MIN)
      {
#       pragma omp for nowait
        for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
          dk[ii] = gk[ii] + (beta * dk[ii]);
          vx[ii] += alpha * dk[ii];
        }
#       pragma omp for nowait
        for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
          zk[ii] = sk[ii] + (beta * zk[ii]);
          rk[ii] += alpha * zk[ii];
        }
      }
    }
    else { /* last iteration */
#     pragma omp parallel for if(n_rows > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        vx[ii] += alpha * (gk[ii] + (beta * dk[ii]));
    }

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  convergence->n_iterations = n_iter;

  return CS_SLES_MAX_ITERATION;
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
  double  ro_0, ro_1, alpha, rk_rkm1, rk_rk, beta;
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

  {
    /* Descent direction */

    cs_matrix_vector_multiply(a, rk, zk); /* rk == dk */

    /* Descent parameter */

    _dot_products_xy_yz(c, rk, rk, zk, &ro_0, &ro_1);

    alpha = -ro_0 / ro_1;
    rk_rkm1 = ro_0;
  }

  if (convergence->n_iterations_max > 1) {
#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        dk[ii] = rk[ii];
#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        vx[ii] += (alpha * rk[ii]);
#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        rk[ii] += (alpha * zk[ii]);
    }
  }
  else {
#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      vx[ii] += (alpha * rk[ii]);
  }

  /* Current Iterations */

  for (n_iter = 1; n_iter < convergence->n_iterations_max; n_iter++) {

    /* Prepare descent parameter */

    rk_rk = _dot_product_xx(c, rk);

    /* Complete descent parameter computation and matrix.vector product */

    beta = rk_rk / rk_rkm1;
    rk_rkm1 = rk_rk;

#   pragma omp parallel for firstprivate(alpha) if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      dk[ii] = rk[ii] + (beta * dk[ii]);

    cs_matrix_vector_multiply(a, dk, zk);

    _dot_products_xy_yz(c, rk, dk, zk, &ro_0, &ro_1);

    alpha = -ro_0 / ro_1;

    if (n_iter + 1 < convergence->n_iterations_max) {
#     pragma omp parallel if(n_rows > CS_THR_MIN)
      {
#       pragma omp for nowait
        for (cs_lnum_t ii = 0; ii < n_rows; ii++)
          vx[ii] += (alpha * dk[ii]);

#       pragma omp for nowait
        for (cs_lnum_t ii = 0; ii < n_rows; ii++)
          rk[ii] += (alpha * zk[ii]);
      }
    }
    else { /* last iteration */
#     pragma omp parallel for if(n_rows > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        vx[ii] += (alpha * dk[ii]);
    }

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  convergence->n_iterations = n_iter;

  return CS_SLES_MAX_ITERATION;
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
  double  ro_0, ro_1, alpha, rk_rkm1, rk_rk, rk_sk, beta;
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

  {
    /* Descent direction */

    cs_matrix_vector_multiply(a, rk, zk); /* zk = A.dk
                                                            dk == rk*/

    /* Descent parameter */

    _dot_products_xy_yz(c, rk, rk, zk, &ro_0, &ro_1);

    alpha = -ro_0 / ro_1;

    rk_rkm1 = ro_0;

    if (convergence->n_iterations_max < 2) {
#     pragma omp parallel for if(n_rows > CS_THR_MIN)
        for (cs_lnum_t ii = 0; ii < n_rows; ii++)
          vx[ii] += (alpha * rk[ii]);
    }
    else {
#     pragma omp parallel for if(n_rows > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
        dk[ii] = rk[ii];
        vx[ii] += (alpha * rk[ii]);
        rk[ii] += (alpha * zk[ii]);
      }
    }
  }

  /* Current Iteration */

  for (n_iter = 1; n_iter < convergence->n_iterations_max; n_iter++) {

    cs_matrix_vector_multiply(a, rk, sk);  /* sk = A.zk */

    /* Descent parameter */

    _dot_products_xx_xy(c, rk, sk, &rk_rk, &rk_sk);

    beta = rk_rk / rk_rkm1;
    rk_rkm1 = rk_rk;

    ro_1 = rk_sk - beta*beta*ro_1;
    ro_0 = rk_rk;

    alpha = -ro_0 / ro_1;

    /* Complete descent parameter computation and matrix.vector product */

    if (n_iter + 1 < convergence->n_iterations_max) {
#     pragma omp parallel if(n_rows > CS_THR_MIN)
      {
#       pragma omp for nowait
        for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
          dk[ii] = beta*dk[ii] + rk[ii];
          vx[ii] += alpha * dk[ii];
        }
#       pragma omp for nowait
        for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
          zk[ii] = beta*zk[ii] + sk[ii];
          rk[ii] += alpha * zk[ii];
        }
      }
    }
    else { /* last iteration */
#     pragma omp parallel for if(n_rows > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        vx[ii] += alpha * (beta*dk[ii] + rk[ii]);
    }

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  convergence->n_iterations = n_iter;

  return CS_SLES_MAX_ITERATION;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using preconditioned 3-layer conjugate residual.
 *
 * Tuned version when used as smoother.
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
  cs_lnum_t  ii;
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

  /* Current Iteration */
  /*-------------------*/

  for (n_iter = 0; n_iter < convergence->n_iterations_max; n_iter++) {

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

    if (n_iter + 1 < convergence->n_iterations_max) {

#     pragma omp parallel firstprivate(alpha, tau, c0, c1) \
  if (n_rows > CS_THR_MIN)
      {
#       pragma omp for nowait
        for (ii = 0; ii < n_rows; ii++) {
          cs_real_t trk = rk[ii];
          rk[ii] = alpha*rk[ii] + c0*rkm1[ii] + c1*zk[ii];
          rkm1[ii] = trk;
        }
#       pragma omp for nowait
        for (ii = 0; ii < n_rows; ii++) {
          cs_real_t tvx = vx[ii];
          vx[ii] = alpha*vx[ii] + c0*vxm1[ii] + c1*wk[ii];
          vxm1[ii] = tvx;
        }
      }

    }
    else {

#     pragma omp parallel firstprivate(alpha, tau, c0, c1) \
  if (n_rows > CS_THR_MIN)
      {
#       pragma omp for nowait
        for (ii = 0; ii < n_rows; ii++) {
          cs_real_t tvx = vx[ii];
          vx[ii] = alpha*vx[ii] + c0*vxm1[ii] + c1*wk[ii];
          vxm1[ii] = tvx;
        }
      }

    }

  } /* Loop on iterations */

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  convergence->n_iterations = n_iter;

  return CS_SLES_MAX_ITERATION;
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
  cs_lnum_t  ii;
  cs_real_t *_aux_vectors;
  cs_real_t *restrict rk;

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

  /* Current iteration */
  /*-------------------*/

  for (n_iter = 0; n_iter < convergence->n_iterations_max; n_iter++) {

#if defined(HAVE_OPENMP)

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (ii = 0; ii < n_rows; ii++)
      rk[ii] = vx[ii];

#else

    memcpy(rk, vx, n_rows * sizeof(cs_real_t));   /* rk <- vx */

#endif

    /* Compute Vx <- Vx - (A-diag).Rk */

    cs_matrix_vector_multiply_partial(a, CS_MATRIX_SPMV_E, rk, vx);

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (ii = 0; ii < n_rows; ii++) {
      vx[ii] = (rhs[ii]-vx[ii])*ad_inv[ii];
    }

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  convergence->n_iterations = n_iter;

  return CS_SLES_MAX_ITERATION;
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

  /* Current iteration */
  /*-------------------*/

  for (n_iter = 0; n_iter < convergence->n_iterations_max; n_iter++) {

    memcpy(rk, vx, n_rows * sizeof(cs_real_t));  /* rk <- vx */

    /* Compute vxx <- vx - (a-diag).rk */

    cs_matrix_vector_multiply_partial(a, CS_MATRIX_SPMV_E, rk, vxx);

    /* Compute vx <- diag^-1 . (vxx - rhs) */
#   pragma omp parallel for if(n_blocks > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_blocks; ii++) {
      _fw_and_bw_lu33(ad_inv + 9*ii,
                      vx + 3*ii,
                      vxx + 3*ii,
                      rhs + 3*ii);
    }

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  convergence->n_iterations = n_iter;

  return CS_SLES_MAX_ITERATION;
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

  /* Current iteration */
  /*-------------------*/

  for (n_iter = 0; n_iter < convergence->n_iterations_max; n_iter++) {

    memcpy(rk, vx, n_rows * sizeof(cs_real_t));  /* rk <- vx */

    /* Compute Vx <- Vx - (A-diag).Rk */

    cs_matrix_vector_multiply_partial(a, CS_MATRIX_SPMV_E, rk, vxx);

#   pragma omp parallel for if(n_blocks > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_blocks; ii++) {
      _fw_and_bw_lu(ad_inv + db_size_2*ii,
                    db_size,
                    vx + db_size*ii,
                    vxx + db_size*ii,
                    rhs + db_size*ii);
    }

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  convergence->n_iterations = n_iter;

  return CS_SLES_MAX_ITERATION;
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
  unsigned n_iter = 0;

  const cs_lnum_t n_rows = cs_matrix_get_n_rows(a);
  const cs_halo_t *halo = cs_matrix_get_halo(a);
  const cs_real_t  *restrict ad_inv = c->setup_data->ad_inv;

  const cs_lnum_t  *a_row_index, *a_col_id;
  const cs_real_t  *a_d_val, *a_x_val;

  const cs_lnum_t db_size = cs_matrix_get_diag_block_size(a);

  cs_matrix_get_msr_arrays(a, &a_row_index, &a_col_id, &a_d_val, &a_x_val);

  const cs_lnum_t  *order = c->add_data->order;

  /* Current iteration */
  /*-------------------*/

  for (n_iter = 0; n_iter < convergence->n_iterations_max; n_iter++) {

    /* Synchronize ghost cells first */

    if (halo != NULL)
      cs_matrix_pre_vector_multiply_sync(a, vx);

    /* Compute Vx <- Vx - (A-diag).Rk */

    if (diag_block_size == 1) {

#     pragma omp parallel for if(n_rows > CS_THR_MIN && !_thread_debug)
      for (cs_lnum_t ll = 0; ll < n_rows; ll++) {

        cs_lnum_t ii = order[ll];

        const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
        const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
        const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

        cs_real_t vx0 = rhs[ii];

        for (cs_lnum_t jj = 0; jj < n_cols; jj++)
          vx0 -= (m_row[jj]*vx[col_id[jj]]);

        vx0 *= ad_inv[ii];
        vx[ii] = vx0;

      }

    }
    else {

      const cs_lnum_t db_size_2 = db_size * db_size;

#     pragma omp parallel for if(n_rows > CS_THR_MIN  && !_thread_debug)
      for (cs_lnum_t ll = 0; ll < n_rows; ll++) {

        cs_lnum_t ii = order[ll];

        const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
        const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
        const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

        cs_real_t vx0[DB_SIZE_MAX], _vx[DB_SIZE_MAX];

        for (cs_lnum_t kk = 0; kk < db_size; kk++)
          vx0[kk] = rhs[ii*db_size + kk];

        for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
          for (cs_lnum_t kk = 0; kk < db_size; kk++)
            vx0[kk] -= (m_row[jj]*vx[col_id[jj]*db_size + kk]);
        }

        _fw_and_bw_lu_gs(ad_inv + db_size_2*ii,
                         db_size,
                         _vx,
                         vx0);

      }

    }

  }

  convergence->n_iterations = n_iter;

  return CS_SLES_MAX_ITERATION;
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
  unsigned n_iter = 0;

  const cs_lnum_t n_rows = cs_matrix_get_n_rows(a);
  const cs_halo_t *halo = cs_matrix_get_halo(a);
  const cs_real_t  *restrict ad_inv = c->setup_data->ad_inv;

  const cs_lnum_t  *a_row_index, *a_col_id;
  const cs_real_t  *a_d_val, *a_x_val;

  const cs_lnum_t db_size = cs_matrix_get_diag_block_size(a);
  cs_matrix_get_msr_arrays(a, &a_row_index, &a_col_id, &a_d_val, &a_x_val);

  /* Current iteration */
  /*-------------------*/

  for (n_iter = 0; n_iter < convergence->n_iterations_max; n_iter++) {

    /* Synchronize ghost cells first */

    if (halo != NULL)
      cs_matrix_pre_vector_multiply_sync(a, vx);

    /* Compute Vx <- Vx - (A-diag).Rk */

    if (diag_block_size == 1) {

#     pragma omp parallel for if(n_rows > CS_THR_MIN && !_thread_debug)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

        const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
        const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
        const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

        cs_real_t vx0 = rhs[ii];

        for (cs_lnum_t jj = 0; jj < n_cols; jj++)
          vx0 -= (m_row[jj]*vx[col_id[jj]]);

        vx0 *= ad_inv[ii];
        vx[ii] = vx0;

      }

    }
    else {

      const cs_lnum_t db_size_2 = db_size * db_size;

#     pragma omp parallel for if(n_rows > CS_THR_MIN && !_thread_debug)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

        const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
        const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
        const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

        cs_real_t vx0[DB_SIZE_MAX], _vx[DB_SIZE_MAX];

        for (cs_lnum_t kk = 0; kk < db_size; kk++)
          vx0[kk] = rhs[ii*db_size + kk];

        for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
          for (cs_lnum_t kk = 0; kk < db_size; kk++)
            vx0[kk] -= (m_row[jj]*vx[col_id[jj]*db_size + kk]);
        }

        _fw_and_bw_lu_gs(ad_inv + db_size_2*ii,
                         db_size,
                         _vx,
                         vx0);

        for (cs_lnum_t kk = 0; kk < db_size; kk++)
          vx[ii*db_size + kk] = _vx[kk];

      }

    }

  }

  convergence->n_iterations = n_iter;

  return CS_SLES_MAX_ITERATION;
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

  /* Check matrix storage type */

  if (cs_matrix_get_type(a) != CS_MATRIX_MSR)
    bft_error
      (__FILE__, __LINE__, 0,
       _("Symmetric Gauss-Seidel Jacobi hybrid solver only supported with a\n"
         "matrix using %s storage."),
       _("MSR"));

  unsigned n_iter = 0;

  const cs_lnum_t n_rows = cs_matrix_get_n_rows(a);
  const cs_halo_t *halo = cs_matrix_get_halo(a);
  const cs_real_t  *restrict ad_inv = c->setup_data->ad_inv;

  const cs_lnum_t  *a_row_index, *a_col_id;
  const cs_real_t  *a_d_val, *a_x_val;

  const cs_lnum_t db_size = cs_matrix_get_diag_block_size(a);
  cs_matrix_get_msr_arrays(a, &a_row_index, &a_col_id, &a_d_val, &a_x_val);

  /* Current iteration */
  /*-------------------*/

  for (n_iter = 0; n_iter < convergence->n_iterations_max; n_iter++) {

    /* Synchronize ghost cells first */

    if (halo != NULL)
      cs_matrix_pre_vector_multiply_sync(a, vx);

    /* Compute Vx <- Vx - (A-diag).Rk: forward step */

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

      const cs_lnum_t db_size_2 = db_size*db_size;

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

    if (diag_block_size == 1) {

#     pragma omp parallel for if(n_rows > CS_THR_MIN && !_thread_debug)
      for (cs_lnum_t ii = n_rows - 1; ii > - 1; ii--) {

        const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
        const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
        const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

        cs_real_t vx0 = rhs[ii];

        for (cs_lnum_t jj = 0; jj < n_cols; jj++)
          vx0 -= (m_row[jj]*vx[col_id[jj]]);

        vx0 *= ad_inv[ii];
        vx[ii] = vx0;

      }

    }
    else {

      const cs_lnum_t db_size_2 = db_size*db_size;

#     pragma omp parallel for if(n_rows > CS_THR_MIN && !_thread_debug)
      for (cs_lnum_t ii = n_rows - 1; ii > - 1; ii--) {

        const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
        const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
        const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

        cs_real_t vx0[DB_SIZE_MAX], _vx[DB_SIZE_MAX];

        for (cs_lnum_t kk = 0; kk < db_size; kk++)
          vx0[kk] = rhs[ii*db_size + kk];

        for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
          for (cs_lnum_t kk = 0; kk < db_size; kk++)
            vx0[kk] -= (m_row[jj]*vx[col_id[jj]*db_size + kk]);
        }

        _fw_and_bw_lu_gs(ad_inv + db_size_2*ii,
                         db_size,
                         _vx,
                         vx0);

        for (cs_lnum_t kk = 0; kk < db_size; kk++)
          vx[ii*db_size + kk] = _vx[kk];

      }

    }

  }

  convergence->n_iterations = n_iter;

  return CS_SLES_MAX_ITERATION;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using Truncated forward Gauss-Seidel.
 *
 * This variant is intended for smoothing with a fixed number of
 * iterations, so does not compute a residual or run a convergence test.
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
_ts_f_gauss_seidel_msr(cs_sles_it_t              *c,
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

  const cs_lnum_t n_rows = cs_matrix_get_n_rows(a);
  const cs_lnum_t n_cols_ext = cs_matrix_get_n_columns(a);

  const cs_real_t  *restrict ad_inv = c->setup_data->ad_inv;

  const cs_lnum_t  *a_row_index, *a_col_id;
  const cs_real_t  *a_d_val, *a_x_val;

  const cs_lnum_t db_size = cs_matrix_get_diag_block_size(a);
  cs_matrix_get_msr_arrays(a, &a_row_index, &a_col_id, &a_d_val, &a_x_val);

  /* Single iteration */
  /*------------------*/

  /* Zeroe ghost cell values first */

  cs_lnum_t s_id = n_rows*diag_block_size;
  cs_lnum_t e_id = n_cols_ext*diag_block_size;

  for (cs_lnum_t ii = s_id; ii < e_id; ii++)
    vx[ii] = 0.0;

  /* Compute Vx <- Vx - (A-diag).Rk */

  if (diag_block_size == 1) {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN && !_thread_debug)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
      const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
      const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

      cs_real_t vx0 = rhs[ii];

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        if (col_id[jj] > ii) break;
        vx0 -= (m_row[jj]*vx[col_id[jj]]);
      }

      vx0 *= ad_inv[ii];

      vx[ii] = vx0;
    }

  }
  else {

    const cs_lnum_t db_size_2 = db_size*db_size;

#   pragma omp parallel for  if(n_rows > CS_THR_MIN && !_thread_debug)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
      const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
      const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

      cs_real_t vx0[DB_SIZE_MAX], _vx[DB_SIZE_MAX];

      for (cs_lnum_t kk = 0; kk < db_size; kk++) {
        vx0[kk] = rhs[ii*db_size + kk];
      }

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        if (col_id[jj] > ii) break;
        for (cs_lnum_t kk = 0; kk < db_size; kk++)
          vx0[kk] -= (m_row[jj]*vx[col_id[jj]*db_size + kk]);
      }

      _fw_and_bw_lu_gs(ad_inv + db_size_2*ii,
                       db_size,
                       _vx,
                       vx0);

      for (cs_lnum_t kk = 0; kk < db_size; kk++) {
        vx[ii*db_size + kk] = _vx[kk];
      }

    }

  }

  convergence->n_iterations = 1;

  return CS_SLES_MAX_ITERATION;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using Truncated backward Gauss-Seidel.
 *
 * This variant is intended for smoothing with a fixed number of
 * iterations, so does not compute a residual or run a convergence test.
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
_ts_b_gauss_seidel_msr(cs_sles_it_t              *c,
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

  const cs_lnum_t n_rows = cs_matrix_get_n_rows(a);
  const cs_lnum_t n_cols_ext = cs_matrix_get_n_columns(a);

  const cs_real_t  *restrict ad_inv = c->setup_data->ad_inv;

  const cs_lnum_t  *a_row_index, *a_col_id;
  const cs_real_t  *a_d_val, *a_x_val;

  const cs_lnum_t db_size = cs_matrix_get_diag_block_size(a);
  cs_matrix_get_msr_arrays(a, &a_row_index, &a_col_id, &a_d_val, &a_x_val);

  /* Single iteration */
  /*------------------*/

  /* Zeroe ghost cell values first */

  cs_lnum_t s_id = n_rows*diag_block_size;
  cs_lnum_t e_id = n_cols_ext*diag_block_size;

  for (cs_lnum_t ii = s_id; ii < e_id; ii++)
    vx[ii] = 0.0;

  /* Compute Vx <- Vx - (A-diag).Rk */

  if (diag_block_size == 1) {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN && !_thread_debug)
    for (cs_lnum_t ii = n_rows - 1; ii > - 1; ii--) {

      const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
      const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
      const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

      cs_real_t vx0 = rhs[ii];

      for (cs_lnum_t jj = n_cols-1; jj > -1; jj--) {
        if (col_id[jj] < ii) break;
        vx0 -= (m_row[jj]*vx[col_id[jj]]);
      }

      vx0 *= ad_inv[ii];

      vx[ii] = vx0;
    }

  }
  else {

    const cs_lnum_t db_size_2 = db_size*db_size;

#   pragma omp parallel for  if(n_rows > CS_THR_MIN && !_thread_debug)
    for (cs_lnum_t ii = n_rows - 1; ii > - 1; ii--) {

      const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
      const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
      const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

      cs_real_t vx0[DB_SIZE_MAX], _vx[DB_SIZE_MAX];

      for (cs_lnum_t kk = 0; kk < db_size; kk++) {
        vx0[kk] = rhs[ii*db_size + kk];
      }

      for (cs_lnum_t jj = n_cols-1; jj > -1; jj--) {
        if (col_id[jj] < ii) break;
        for (cs_lnum_t kk = 0; kk < db_size; kk++)
          vx0[kk] -= (m_row[jj]*vx[col_id[jj]*db_size + kk]);
      }

      _fw_and_bw_lu_gs(ad_inv + db_size_2*ii,
                       db_size,
                       _vx,
                       vx0);

      for (cs_lnum_t kk = 0; kk < db_size; kk++) {
        vx[ii*db_size + kk] = _vx[kk];
      }

    }

  }

  convergence->n_iterations = 1;

  return CS_SLES_MAX_ITERATION;
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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create iterative sparse linear system solver info and context used as
 *        a smoother in a multigrid solver
 *
 * parameters:
 * \param[in]  smoother_type   type of smoother (CG, Jacobi, ...)
 * \param[in]  poly_degree     preconditioning polynomial degree
 *                             (0: diagonal; -1: non-preconditioned;
 *                             see \ref sles_it for details)
 * \param[in]  n_iter          number of iterations to perform
 *
 * \return a pointer to newly created smoother info object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_it_t *
cs_multigrid_smoother_create(cs_sles_it_type_t    smoother_type,
                             int                  poly_degree,
                             int                  n_iter)
{
  cs_sles_it_t *c;

  BFT_MALLOC(c, 1, cs_sles_it_t);

  c->solve = NULL;
  c->_pc = NULL;

  /* Predefined settings */
  c->type = smoother_type;
  c->on_device = false;
  c->update_stats = false;
  c->ignore_convergence = true;

  c->fallback_cvg = CS_SLES_DIVERGED;
  c->fallback_n_max_iter = 0;
  c->fallback = NULL;

  switch (smoother_type) {      /* Valid choices */

  case CS_SLES_JACOBI:
  case CS_SLES_P_GAUSS_SEIDEL:
  case CS_SLES_P_SYM_GAUSS_SEIDEL:
  case CS_SLES_TS_F_GAUSS_SEIDEL:
  case CS_SLES_TS_B_GAUSS_SEIDEL:
    break;

  case CS_SLES_PCG:
    if (poly_degree < 0)
      c->_pc = cs_sles_pc_none_create();
    else if (poly_degree == 0)
      c->_pc = cs_sles_pc_jacobi_create();
    else if (poly_degree == 1)
      c->_pc = cs_sles_pc_poly_1_create();
    else
      c->_pc = cs_sles_pc_poly_2_create();
    break;

  case CS_SLES_PCR3:
    if (poly_degree < 0)
      c->_pc = cs_sles_pc_none_create();
    else if (poly_degree == 0)
      c->_pc = cs_sles_pc_jacobi_create();
    else if (poly_degree == 1)
      c->_pc = cs_sles_pc_poly_1_create();
    else
      c->_pc = cs_sles_pc_poly_2_create();
    break;

  default: /* Other iterative solvers are not tuned for smoothing */
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid smoother.", __func__);
    break;

  } /* Smoother type */

  c->pc = c->_pc;

  c->n_max_iter = n_iter;
  c->restart_interval = 20;

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

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup iterative sparse linear equation solver in case of used as a
 *        multigrid smoother
 *
 * \param[in, out]  context    pointer to iterative solver info and context
 *                             (actual type: cs_sles_it_t  *)
 * \param[in]       name       pointer to system name
 * \param[in]       a          associated matrix
 * \param[in]       verbosity  associated verbosity
 */
/*----------------------------------------------------------------------------*/

void
cs_multigrid_smoother_setup(void               *context,
                            const char         *name,
                            const cs_matrix_t  *a,
                            int                 verbosity)
{
  cs_sles_it_t  *c = context;

  const int diag_block_size = cs_matrix_get_diag_block_size(a);

  if (verbosity > 1) {
    bft_printf(_("\n Setup of solver for linear system \"%s\"\n"),
               name);
    cs_matrix_log_info(a, verbosity);
  }

  if (c->type == CS_SLES_JACOBI)
    cs_sles_it_setup_priv(c, name, a, verbosity, diag_block_size, true);

  else if (   c->type == CS_SLES_P_GAUSS_SEIDEL
           || c->type == CS_SLES_P_SYM_GAUSS_SEIDEL) {
    /* Force to Jacobi type in case matrix type is not adapted */
    if (cs_matrix_get_type(a) != CS_MATRIX_MSR)
      c->type = CS_SLES_JACOBI;
    cs_sles_it_setup_priv(c, name, a, verbosity, diag_block_size, true);
  }

  else if (   c->type == CS_SLES_TS_F_GAUSS_SEIDEL
           || c->type == CS_SLES_TS_B_GAUSS_SEIDEL) {
    /* Force to closest Jacobi type in case matrix type is not adapted */
    if (cs_matrix_get_type(a) != CS_MATRIX_MSR) {
      c->type = CS_SLES_JACOBI;
      c->n_max_iter = 2;
    }
    cs_sles_it_setup_priv(c, name, a, verbosity, diag_block_size, true);
  }

  else
    cs_sles_it_setup_priv(c, name, a, verbosity, diag_block_size, false);

  switch (c->type) {

  case CS_SLES_PCG:
    /* Check for single-reduction */
    {
      bool single_reduce = false;
#if defined(HAVE_MPI)
      cs_gnum_t n_m_rows = c->setup_data->n_rows;
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

  case CS_SLES_PCR3:
    c->solve = _conjugate_residual_3;
    break;

  case CS_SLES_JACOBI:
    if (diag_block_size == 1)
      c->solve = _jacobi;
    else if (diag_block_size == 3)
      c->solve = _block_3_jacobi;
    else
      c->solve = _block_jacobi;
    break;

  case CS_SLES_P_GAUSS_SEIDEL:
    c->solve = _p_gauss_seidel;
    break;
  case CS_SLES_P_SYM_GAUSS_SEIDEL:
    c->solve = _p_sym_gauss_seidel_msr;
    break;

  case CS_SLES_TS_F_GAUSS_SEIDEL:
    c->solve = _ts_f_gauss_seidel_msr;
    break;
  case CS_SLES_TS_B_GAUSS_SEIDEL:
    c->solve = _ts_b_gauss_seidel_msr;
    break;

  default:
    bft_error
      (__FILE__, __LINE__, 0,
       _(" %s: Setup of linear equation on \"%s\"\n"
         "with smoother type %d, which is not allowed or available)."),
       __func__, name, (int)c->type);
    break;
  }

  /* Now finish */
  assert(c->update_stats == false);
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
cs_multigrid_smoother_solve(void                *context,
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

  cs_sles_it_convergence_t  convergence;

  const cs_lnum_t diag_block_size = cs_matrix_get_diag_block_size(a);

  /* Initialize number of iterations and residual,
     and smooth linear system */

  *n_iter = 0;
  *residual = -1.0; /* Don't use this quantity when dealing with smoothers */

  /* Setup if not already done */

  if (c->setup_data == NULL)
    cs_sles_it_setup(c, name, a, verbosity);

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
                  _(" RHS norm:          %11.4e\n\n"), r_norm);

  /* Only call solver for "active" ranks */

#if defined(HAVE_MPI)
  if (c->caller_n_ranks < 2 || c->comm != MPI_COMM_NULL) {
#endif

    cvg = c->solve(c,
                   a, diag_block_size, &convergence,
                   rhs, vx,
                   aux_size, aux_vectors);

#if defined(HAVE_MPI)
  }
#endif

  /* Update return values */

  *n_iter = convergence.n_iterations;
  *residual = convergence.residual;

  return cvg;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
