/*============================================================================
 * Sparse Linear Equation Solvers: Multigrid smoothers
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

#include "base/cs_mem.h"
#include "base/cs_dispatch.h"

#if defined(__CUDACC__)
#include "alge/cs_blas_cuda.h"
#include "alge/cs_matrix_spmv.h"
#include "alge/cs_matrix_spmv_cuda.h"
#endif

#if defined(__HIPCC__)
#include "alge/cs_blas_hip.h"
#include "alge/cs_matrix_spmv.h"
#include "alge/cs_matrix_spmv_hip.h"
#endif

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "alge/cs_sles_it.h"
#include "alge/cs_sles_it_priv.h"
#include "alge/cs_multigrid_smoother.h"

#if (HAVE_CUDA)
#include "alge/cs_sles_it_cuda.h"
#endif

#if (HAVE_HIP)
#include "alge/cs_sles_it_hip.h"
#endif

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_multigrid_smoother.cpp
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

[[maybe_unused]] static const int _thread_debug = 0;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using Jacobi.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- diagonal block size
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx_ini          <-- initial system solution
 *                       (vx if nonzero, nullptr if zero)
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
        cs_real_t                 *restrict vx_ini,
        cs_real_t                 *restrict vx,
        size_t                     aux_size,
        void                      *aux_vectors)
{
  cs_real_t *_aux_vectors;
  cs_real_t *restrict rk;

  unsigned n_iter = 0;

  /* Allocate or map work arrays */
  /*-----------------------------*/

  assert(c->setup_data != nullptr);

  const cs_real_t  *restrict ad_inv = c->setup_data->ad_inv;

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 1;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == nullptr || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      CS_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = (cs_real_t *)aux_vectors;

    rk = _aux_vectors;
  }

  int iter_ini = 0;

  /* First iteration simplified if vx == 0
     ------------------------------------- */

  if (vx_ini != vx) {
    assert(vx_ini == nullptr);

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
      vx[ii] = rhs[ii]*ad_inv[ii];
    }

    iter_ini = 1;
  }

  /* Current iteration */
  /*-------------------*/

  for (n_iter = iter_ini; n_iter < convergence->n_iterations_max; n_iter++) {

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      rk[ii] = vx[ii];

    /* Compute Vx <- Vx - (A-diag).Rk */

    cs_matrix_vector_multiply_partial(a, CS_MATRIX_SPMV_E, rk, vx);

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
      vx[ii] = (rhs[ii]-vx[ii])*ad_inv[ii];
    }

  }

  if (_aux_vectors != aux_vectors)
    CS_FREE(_aux_vectors);

  convergence->n_iterations = n_iter;

  return CS_SLES_MAX_ITERATION;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using relaxed (weighed) Jacobi.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- diagonal block size
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx_ini          <-- initial system solution
 *                       (vx if nonzero, nullptr if zero)
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_relaxed_jacobi(cs_sles_it_t              *c,
                const cs_matrix_t         *a,
                cs_lnum_t                  diag_block_size,
                cs_sles_it_convergence_t  *convergence,
                const cs_real_t           *rhs,
                cs_real_t                 *vx_ini,
                cs_real_t                 *vx,
                size_t                     aux_size,
                void                      *aux_vectors)
{
  cs_real_t *_aux_vectors;
  cs_real_t *restrict vx_tmp;

  cs_real_t wk = 1.0;
  int wk_m = 1;

  cs_real_t w[3] = {2./3., 1.0, 1.0};

  // Scheduled relaxation variants:
  // https://dx.doi.org/10.1016/j.jcp.2016.12.010
  // https://doi.org/10.1007/s12046-023-02407-6

  if (c->type == CS_SLES_RJ2) {
    w[0] = 1.7319;
    w[1] = 0.5695;
    wk_m = 2;
  }
  if (c->type == CS_SLES_RJ3) {
    w[0] = 2.2473;
    w[1] = 0.8573;
    w[2] = 0.5296;
    wk_m = 3;
  }

  int n_iter = 0;

  cs_dispatch_context ctx;

#if defined(HAVE_ACCEL)
  bool local_stream;
  cs_stream_t stream;
  cs_sles_it_set_exec_location(ctx, a, local_stream, stream);
#endif

  /* Allocate or map work arrays */
  /*-----------------------------*/

  assert(c->setup_data != nullptr);

  const cs_real_t  *restrict ad_inv = c->setup_data->ad_inv;

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 1;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == nullptr || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      CS_MALLOC_HD(_aux_vectors, wa_size * n_wa, cs_real_t, ctx.alloc_mode());
    else
      _aux_vectors = (cs_real_t *)aux_vectors;

    vx_tmp = _aux_vectors;
  }

  int iter_ini = 0;

  /* First iteration simplified if vx == 0
     ------------------------------------- */

  cs_real_t *restrict vx_n = vx;
  cs_real_t *restrict vx_np = vx_tmp;

  if (vx_ini != vx) {
    assert(vx_ini == nullptr);
    wk = w[0];

    ctx.parallel_for(n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
      vx_n[ii] = wk * rhs[ii]*ad_inv[ii];
    });

    iter_ini = 1;
  }

  /* Current iteration */
  /*-------------------*/

  for (n_iter = iter_ini; n_iter < c->n_max_iter; n_iter++) {

    wk = w[n_iter % wk_m];
    cs_real_t o_m_wk = 1.0 - wk;

    /* Swap vx_np and vx_n */
    {cs_real_t *vx_t = vx_n; vx_n = vx_np; vx_np = vx_t;}

    /* Compute Vx <- (1-w).Vx + w.(A-diag).Rk */

    cs_matrix_vector_multiply_partial(ctx, a, CS_MATRIX_SPMV_E, vx_np, vx_n);

    ctx.parallel_for(n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
      vx_n[ii] = o_m_wk*vx_np[ii] + wk*(rhs[ii]-vx_n[ii])*ad_inv[ii];
    });

  }

  if (vx_n != vx) {
    ctx.parallel_for(n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
      vx[ii] = vx_n[ii];
    });
  }

  ctx.wait();

#if defined(HAVE_ACCEL)
  cs_sles_it_restore_exec_location(local_stream);
#endif

  if (_aux_vectors != aux_vectors)
    CS_FREE(_aux_vectors);

  convergence->n_iterations = n_iter;

  return CS_SLES_MAX_ITERATION;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using block Jacobi.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- diagonal block size (unused here)
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx_ini          <-- initial system solution
 *                       (vx if nonzero, nullptr if zero)
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
                cs_real_t                 *restrict vx_ini,
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

  assert(c->setup_data != nullptr);

  const cs_real_t  *restrict ad_inv = c->setup_data->ad_inv;

  const cs_lnum_t n_blocks = c->setup_data->n_rows / 3;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 2;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == nullptr || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      CS_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = static_cast<cs_real_t *>(aux_vectors);

    rk  = _aux_vectors;
    vxx = _aux_vectors + wa_size;
  }

  int iter_ini = 0;

  /* First iteration simplified if vx == 0
     ------------------------------------- */

  if (vx_ini != vx) {
    assert(vx_ini == nullptr);

    /* Compute vx <- diag^-1 . (vxx - rhs) */
#   pragma omp parallel for if(n_blocks > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_blocks; ii++) {
      for (cs_lnum_t jj = 0; jj < 3; jj++) {
        rk[3*ii + jj] = 0;
        vxx[3*ii + jj] = 0;
      }

      _mat_c_m_b_33(ad_inv + 9*ii,
                    vx + 3*ii,
                    vxx + 3*ii,
                    rhs + 3*ii);
    }

    iter_ini = 1;
  }

  /* Current iteration */
  /*-------------------*/

  for (n_iter = iter_ini; n_iter < convergence->n_iterations_max; n_iter++) {

    /* Compute vxx <- vx - (a-diag).rk */

    cs_matrix_vector_multiply_partial(a, CS_MATRIX_SPMV_E, vx, vxx);

    /* Compute vx <- diag^-1 . (vxx - rhs) */
#   pragma omp parallel for if(n_blocks > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_blocks; ii++) {
      for (cs_lnum_t kk = 0; kk < 3; kk++)
        rk[ii*3 + kk] = vx[ii*3 + kk];
      _mat_c_m_b_33(ad_inv + 9*ii,
                    vx + 3*ii,
                    vxx + 3*ii,
                    rhs + 3*ii);
    }

  }

  if (_aux_vectors != aux_vectors)
    CS_FREE(_aux_vectors);

  convergence->n_iterations = n_iter;

  return CS_SLES_MAX_ITERATION;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using block Jacobi.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- block size of diagonal elements
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx_ini          <-- initial system solution
 *                       (vx if nonzero, nullptr if zero)
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
              cs_real_t                 *restrict vx_ini,
              cs_real_t                 *restrict vx,
              size_t                     aux_size,
              void                      *aux_vectors)
{
  cs_real_t *_aux_vectors;
  cs_real_t  *restrict rk, *restrict vxx;

  unsigned n_iter = 0;

  /* Call setup if not already done, allocate or map work arrays */
  /*-------------------------------------------------------------*/
  assert(c->setup_data != nullptr);

  const cs_lnum_t db_size = cs_matrix_get_diag_block_size(a);
  const cs_lnum_t db_size_2 = db_size * db_size;

  const cs_real_t  *restrict ad_inv = c->setup_data->ad_inv;

  const cs_lnum_t n_blocks = c->setup_data->n_rows / diag_block_size;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 2;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == nullptr || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      CS_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = static_cast<cs_real_t *>(aux_vectors);

    rk  = _aux_vectors;
    vxx = _aux_vectors + wa_size;
  }

  int iter_ini = 0;

  /* First iteration simplified if vx == 0
     ------------------------------------- */

  if (vx_ini != vx) {
    assert(vx_ini == nullptr);

    /* Compute vx <- diag^-1 . (vxx - rhs) */
#   pragma omp parallel for if(n_blocks > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_blocks; ii++) {
      assert(db_size <= DB_SIZE_MAX);
      cs_real_t _vxx[DB_SIZE_MAX];

      for (cs_lnum_t jj = 0; jj < db_size; jj++) {
        rk[db_size*ii + jj] = 0;
        _vxx[jj] = 0;
      }

      _mat_c_m_b(ad_inv + db_size_2*ii,
                 db_size,
                 vx + db_size*ii,
                 _vxx,
                 rhs + db_size*ii);
    }

    iter_ini = 1;
  }

  /* Current iteration */
  /*-------------------*/

  for (n_iter = iter_ini; n_iter < convergence->n_iterations_max; n_iter++) {

    /* Compute Vx <- Vx - (A-diag).Rk */

    cs_matrix_vector_multiply_partial(a, CS_MATRIX_SPMV_E, vx, vxx);

#   pragma omp parallel for if(n_blocks > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_blocks; ii++) {
      for (cs_lnum_t kk = 0; kk < 3; kk++)
        rk[ii*db_size + kk] = vx[ii*db_size + kk];
      _mat_c_m_b(ad_inv + db_size_2*ii,
                 db_size,
                 vx + db_size*ii,
                 vxx + db_size*ii,
                 rhs + db_size*ii);
    }

  }

  if (_aux_vectors != aux_vectors)
    CS_FREE(_aux_vectors);

  convergence->n_iterations = n_iter;

  return CS_SLES_MAX_ITERATION;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using relaxed block Jacobi.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- block size of diagonal elements
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx_ini          <-- initial system solution
 *                       (vx if nonzero, nullptr if zero)
 *   vx              --> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_block_relaxed_jacobi(cs_sles_it_t              *c,
                      const cs_matrix_t         *a,
                      cs_lnum_t                  diag_block_size,
                      cs_sles_it_convergence_t  *convergence,
                      const cs_real_t           *rhs,
                      cs_real_t                 *restrict vx_ini,
                      cs_real_t                 *restrict vx,
                      size_t                     aux_size,
                      void                      *aux_vectors)
{
  cs_real_t *_aux_vectors = nullptr;
  cs_real_t *restrict vx_tmp = nullptr, *vxx = nullptr;

  cs_real_t wk = 1.0;
  int wk_m = 1;

  cs_real_t w[3] = {2./3., 1.0, 1.0};

  // Scheduled relaxation variants:
  // https://dx.doi.org/10.1016/j.jcp.2016.12.010
  // https://doi.org/10.1007/s12046-023-02407-6
  // Need to check if they are appropriate/relevant here.

  if (c->type == CS_SLES_RJ2) {
    w[0] = 1.7319;
    w[1] = 0.5695;
    wk_m = 2;
  }
  if (c->type == CS_SLES_RJ3) {
    w[0] = 2.2473;
    w[1] = 0.8573;
    w[2] = 0.5296;
    wk_m = 3;
  }

  unsigned n_iter = 0;

  cs_dispatch_context ctx;

#if defined(HAVE_ACCEL)
  bool local_stream;
  cs_stream_t stream;
  cs_sles_it_set_exec_location(ctx, a, local_stream, stream);
#endif

  /* Allocate or map work arrays */
  /*-----------------------------*/

  assert(c->setup_data != nullptr);
  assert(diag_block_size == cs_matrix_get_diag_block_size(a));

  const cs_lnum_t db_size = diag_block_size;
  const cs_lnum_t db_size_2 = db_size * db_size;

  const cs_real_t  *restrict ad_inv = c->setup_data->ad_inv;

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * db_size;
    const size_t n_wa = 2;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == nullptr || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      CS_MALLOC_HD(_aux_vectors, wa_size * n_wa, cs_real_t, ctx.alloc_mode());
    else
      _aux_vectors = static_cast<cs_real_t *>(aux_vectors);

    vx_tmp = _aux_vectors;
    vxx = _aux_vectors + wa_size;
  }

  int iter_ini = 0;

  /* First iteration simplified if vx == 0
     ------------------------------------- */

  cs_real_t *restrict vx_n = vx;
  cs_real_t *restrict vx_np = vx_tmp;

  if (vx_ini != vx) {
    assert(vx_ini == nullptr);
    wk = w[0];

    /* Compute vx <- diag^-1 . (vxx - rhs) */
    ctx.parallel_for(n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t r_ii) {
      const cs_lnum_t ii = r_ii / db_size;
      const cs_lnum_t jj = r_ii % db_size;

      const cs_real_t *_ad_inv = ad_inv + db_size_2*ii + db_size*jj;
      const cs_real_t *_rhs = rhs + db_size*ii;
      cs_real_t _vx = 0;

      for (cs_lnum_t kk = 0; kk < db_size; kk++) {
        _vx += _rhs[kk] * _ad_inv[kk];
      }

      vx_n[r_ii] = wk * _vx;
    });

    iter_ini = 1;
  }

  /* Current iteration */
  /*-------------------*/

  for (n_iter = iter_ini; n_iter < convergence->n_iterations_max; n_iter++) {

    wk = w[n_iter % wk_m];
    cs_real_t o_m_wk = 1.0 - wk;

    /* Swap vx_np and vx_n */
    {cs_real_t *vx_t = vx_n; vx_n = vx_np; vx_np = vx_t;}

    /* Compute Vx <- (1-w).Vx + w.(A-diag).Rk */

    cs_matrix_vector_multiply_partial(ctx, a, CS_MATRIX_SPMV_E, vx_np, vxx);

    ctx.parallel_for(n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t r_ii) {
      const cs_lnum_t ii = r_ii / db_size;
      const cs_lnum_t jj = r_ii % db_size;

      const cs_real_t *_ad_inv = ad_inv + db_size_2*ii + db_size*jj;
      const cs_real_t *_rhs = rhs + db_size*ii;
      const cs_real_t *_vxx = vxx + db_size*ii;

      assert(db_size <= DB_SIZE_MAX);

      cs_real_t _vx = 0;
      for (cs_lnum_t kk = 0; kk < db_size; kk++) {
        _vx += _ad_inv[kk] * (_rhs[kk] - _vxx[kk]);
      }

      vx_n[r_ii] = o_m_wk*vx_np[r_ii] + wk*_vx;
    });

  }

  if (vx_n != vx) {
    ctx.parallel_for(n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
      vx[ii] = vx_n[ii];
    });
  }

  ctx.wait();

#if defined(HAVE_ACCEL)
  cs_sles_it_restore_exec_location(local_stream);
#endif

  if (_aux_vectors != aux_vectors)
    CS_FREE(_aux_vectors);

  convergence->n_iterations = n_iter;

  return CS_SLES_MAX_ITERATION;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using Jacobi.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- diagonal block size
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx_ini          <-- initial system solution
 *                       (vx if nonzero, nullptr if zero)
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_l1_jacobi(cs_sles_it_t              *c,
           const cs_matrix_t         *a,
           cs_lnum_t                  diag_block_size,
           cs_sles_it_convergence_t  *convergence,
           const cs_real_t           *rhs,
           cs_real_t                 *vx_ini,
           cs_real_t                 *vx,
           size_t                     aux_size,
           void                      *aux_vectors)
{
  cs_real_t *_aux_vectors = nullptr;
  cs_real_t *vx_tmp = nullptr;

  unsigned n_iter = 0;

  cs_dispatch_context ctx;

#if defined(HAVE_ACCEL)
  bool local_stream;
  cs_stream_t stream;
  cs_sles_it_set_exec_location(ctx, a, local_stream, stream);
#endif

  /* Allocate or map work arrays */
  /*-----------------------------*/

  assert(c->setup_data != nullptr);

  const cs_real_t  *restrict l1_inv = c->setup_data->ad_inv;

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 1;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == nullptr || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      CS_MALLOC_HD(_aux_vectors, wa_size * n_wa, cs_real_t, ctx.alloc_mode());
    else
      _aux_vectors = (cs_real_t *)aux_vectors;

    vx_tmp = _aux_vectors;
  }

  int iter_ini = 0;

  /* First iteration simplified if vx == 0
     ------------------------------------- */

  cs_real_t *restrict vx_n = vx;
  cs_real_t *restrict vx_np = vx_tmp;

  if (vx_ini != vx) {
    assert(vx_ini == nullptr);

    ctx.parallel_for(n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
      vx_n[ii] = rhs[ii]*l1_inv[ii];
    });

    iter_ini = 1;
  }

  /* Current iteration */
  /*-------------------*/

  for (n_iter = iter_ini; n_iter < convergence->n_iterations_max; n_iter++) {

    /* Swap vx_np and vx_n */
    {cs_real_t *vx_t = vx_n; vx_n = vx_np; vx_np = vx_t;}

    /* Compute Vx <- Vx - (B - A.Rk)/L1 */

    cs_matrix_vector_multiply(a, vx_np, vx_n);

    ctx.parallel_for(n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
      vx_n[ii] = vx_np[ii] + (rhs[ii]-vx_n[ii]) * l1_inv[ii];
    });

  }

  if (vx_n != vx) {
    ctx.parallel_for(n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
      vx[ii] = vx_n[ii];
    });
  }

  ctx.wait();

  if (_aux_vectors != aux_vectors)
    CS_FREE(_aux_vectors);

  convergence->n_iterations = n_iter;

#if defined(HAVE_ACCEL)
  cs_sles_it_restore_exec_location(local_stream);
#endif

  return CS_SLES_MAX_ITERATION;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using Process-local Gauss-Seidel.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- diagonal block size
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx_ini          <-- initial system solution
 *                       (vx if nonzero, nullptr if zero)
 *   vx              <-> system solution
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
                            cs_real_t                 *restrict vx_ini,
                            cs_real_t                 *restrict vx)
#if defined(__has_feature)
#  if __has_feature(thread_sanitizer)
  __attribute__((no_sanitize("thread")))
#  endif
#endif
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

    if (n_iter == 0 && vx_ini != vx) {
      const cs_lnum_t _n_cols = cs_matrix_get_n_columns(a)*diag_block_size;
#     pragma omp parallel for if(_n_cols > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < _n_cols; ii++) {
        vx[ii] = 0;
      }
    }
    else if (halo != nullptr) {
      /* Synchronize ghost cells first */
      cs_matrix_pre_vector_multiply_sync(a, vx);
    }

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

        const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];
        const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
        const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
        const cs_real_t *restrict _ad_inv = ad_inv + db_size_2*ii;

        cs_real_t vx0[DB_SIZE_MAX], _vx[DB_SIZE_MAX];

        for (cs_lnum_t kk = 0; kk < db_size; kk++)
          vx0[kk] = rhs[ii*db_size + kk];

        for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
          for (cs_lnum_t kk = 0; kk < db_size; kk++)
            vx0[kk] -= (m_row[jj]*vx[col_id[jj]*db_size + kk]);
        }

        for (cs_lnum_t jj = 0; jj < db_size; jj++) {
          _vx[jj] = 0;
          for (cs_lnum_t kk = 0; kk < db_size; kk++)
            _vx[jj] += _ad_inv[jj*db_size + kk] * vx0[kk];
        }

      }

    }

  }

  convergence->n_iterations = n_iter;

  return CS_SLES_MAX_ITERATION;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using Process-local Gauss-Seidel.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- diagonal block size
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx_ini          <-- initial system solution
 *                       (vx if nonzero, nullptr if zero)
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
                    cs_real_t                 *restrict vx_ini,
                    cs_real_t                 *restrict vx)
#if defined(__has_feature)
#  if __has_feature(thread_sanitizer)
  __attribute__((no_sanitize("thread")))
#  endif
#endif
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

    if (n_iter == 0 && vx_ini != vx) {
      const cs_lnum_t _n_cols = cs_matrix_get_n_columns(a)*diag_block_size;
#     pragma omp parallel for if(_n_cols > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < _n_cols; ii++) {
        vx[ii] = 0;
      }
    }
    else if (halo != nullptr) {
      /* Synchronize ghost cells first */
      cs_matrix_pre_vector_multiply_sync(a, vx);
    }

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

        const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];
        const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
        const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
        const cs_real_t *restrict _ad_inv = ad_inv + db_size_2*ii;

        cs_real_t vx0[DB_SIZE_MAX], _vx[DB_SIZE_MAX];

        for (cs_lnum_t kk = 0; kk < db_size; kk++)
          vx0[kk] = rhs[ii*db_size + kk];

        for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
          for (cs_lnum_t kk = 0; kk < db_size; kk++)
            vx0[kk] -= (m_row[jj]*vx[col_id[jj]*db_size + kk]);
        }

        for (cs_lnum_t jj = 0; jj < db_size; jj++) {
          _vx[jj] = 0;
          for (cs_lnum_t kk = 0; kk < db_size; kk++)
            _vx[jj] += _ad_inv[jj*db_size + kk] * vx0[kk];
        }

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
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- diagonal block size
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx_ini          <-- initial system solution
 *                       (vx if nonzero, nullptr if zero)
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
                        cs_real_t                 *restrict vx_ini,
                        cs_real_t                 *restrict vx,
                        size_t                     aux_size,
                        void                      *aux_vectors)
#if defined(__has_feature)
#  if __has_feature(thread_sanitizer)
  __attribute__((no_sanitize("thread")))
#  endif
#endif
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

    if (n_iter == 0 && vx_ini != vx) {
      const cs_lnum_t _n_cols = cs_matrix_get_n_columns(a)*diag_block_size;
#     pragma omp parallel for if(_n_cols > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < _n_cols; ii++) {
        vx[ii] = 0;
      }
    }
    else if (halo != nullptr) {
      /* Synchronize ghost cells first */
      cs_matrix_pre_vector_multiply_sync(a, vx);
    }

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

        const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];
        const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
        const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
        const cs_real_t *restrict _ad_inv = ad_inv + db_size_2*ii;

        cs_real_t vx0[DB_SIZE_MAX], _vx[DB_SIZE_MAX];

        for (cs_lnum_t jj = 0; jj < diag_block_size; jj++)
          vx0[jj] = rhs[ii*db_size + jj];

        for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
          for (cs_lnum_t kk = 0; kk < diag_block_size; kk++)
            vx0[kk] -= (m_row[jj]*vx[col_id[jj]*db_size + kk]);
        }

        for (cs_lnum_t jj = 0; jj < db_size; jj++) {
          _vx[jj] = 0;
          for (cs_lnum_t kk = 0; kk < db_size; kk++) {
            _vx[jj] += _ad_inv[jj*db_size + kk] * vx0[kk];
            assert(fabs(_ad_inv[jj*db_size + kk] * vx0[kk]) < 1e30);
          }
        }

        for (cs_lnum_t jj = 0; jj < diag_block_size; jj++)
          vx[ii*db_size + jj] = _vx[jj];

      }

    }

    /* Synchronize ghost cells again */

    if (halo != nullptr)
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

        const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];
        const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
        const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
        const cs_real_t *restrict _ad_inv = ad_inv + db_size_2*ii;

        cs_real_t vx0[DB_SIZE_MAX], _vx[DB_SIZE_MAX];

        for (cs_lnum_t jj = 0; jj < db_size; jj++)
          vx0[jj] = rhs[ii*db_size + jj];

        for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
          for (cs_lnum_t kk = 0; kk < db_size; kk++)
            vx0[kk] -= (m_row[jj]*vx[col_id[jj]*db_size + kk]);
        }

        for (cs_lnum_t jj = 0; jj < db_size; jj++) {
          _vx[jj] = 0;
          for (cs_lnum_t kk = 0; kk < db_size; kk++) {
            _vx[jj] += _ad_inv[jj*db_size + kk] * vx0[kk];
            assert(fabs(_ad_inv[jj*db_size + kk] * vx0[kk]) < 1e30);
          }
        }

        for (cs_lnum_t jj = 0; jj < db_size; jj++) {
          vx[ii*db_size + jj] = _vx[jj];
          assert(fabs(vx[ii*db_size + jj]) < 1e30);
        }

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
 *   vx_ini          <-- initial system solution
 *                       (vx if nonzero, nullptr if zero)
 *   vx              --> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (unused here)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_ts_f_gauss_seidel_msr(cs_sles_it_t                *c,
                       const cs_matrix_t           *a,
                       cs_lnum_t                    diag_block_size,
                       cs_sles_it_convergence_t    *convergence,
                       const cs_real_t             *rhs,
                       cs_real_t                   *restrict vx_ini,
                       cs_real_t                   *restrict vx,
                       size_t                       aux_size,
                       void                        *aux_vectors)
#if defined(__has_feature)
#  if __has_feature(thread_sanitizer)
  __attribute__((no_sanitize("thread")))
#  endif
#endif
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

  cs_lnum_t s_id = (vx_ini == vx) ? n_rows*diag_block_size : 0;
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

      const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];
      const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
      const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
      const cs_real_t *restrict _ad_inv = ad_inv + db_size_2*ii;

      cs_real_t vx0[DB_SIZE_MAX], _vx[DB_SIZE_MAX];

      for (cs_lnum_t jj = 0; jj < db_size; jj++) {
        vx0[jj] = rhs[ii*db_size + jj];
      }

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        if (col_id[jj] > ii) break;
        for (cs_lnum_t kk = 0; kk < db_size; kk++)
          vx0[kk] -= (m_row[jj]*vx[col_id[jj]*db_size + kk]);
      }

      for (cs_lnum_t jj = 0; jj < db_size; jj++) {
        _vx[jj] = 0;
        for (cs_lnum_t kk = 0; kk < db_size; kk++)
          _vx[jj] += _ad_inv[jj*db_size + kk] * vx0[kk];
      }

      for (cs_lnum_t jj = 0; jj < db_size; jj++) {
        vx[ii*db_size + jj] = _vx[jj];
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
 *   vx_ini          <-- initial system solution
 *                       (vx if nonzero, nullptr if zero)
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
                       cs_real_t                 *restrict vx_ini,
                       cs_real_t                 *restrict vx,
                       size_t                     aux_size,
                       void                      *aux_vectors)
#if defined(__has_feature)
#  if __has_feature(thread_sanitizer)
  __attribute__((no_sanitize("thread")))
#  endif
#endif
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

  /* Zero ghost cell values first */

  cs_lnum_t s_id = (vx_ini == vx) ? n_rows*diag_block_size : 0;
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

      const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];
      const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
      const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
      const cs_real_t *restrict _ad_inv = ad_inv + db_size_2*ii;

      cs_real_t vx0[DB_SIZE_MAX], _vx[DB_SIZE_MAX];

      for (cs_lnum_t jj = 0; jj < db_size; jj++) {
        vx0[jj] = rhs[ii*db_size + jj];
      }

      for (cs_lnum_t jj = n_cols-1; jj > -1; jj--) {
        if (col_id[jj] < ii) break;
        for (cs_lnum_t kk = 0; kk < db_size; kk++)
          vx0[kk] -= (m_row[jj]*vx[col_id[jj]*db_size + kk]);
      }

      for (cs_lnum_t jj = 0; jj < db_size; jj++) {
        _vx[jj] = 0;
        for (cs_lnum_t kk = 0; kk < db_size; kk++)
          _vx[jj] += _ad_inv[jj*db_size + kk] * vx0[kk];
      }

      for (cs_lnum_t jj = 0; jj < db_size; jj++) {
        vx[ii*db_size + jj] = _vx[jj];
      }

    }

  }

  convergence->n_iterations = 1;

  return CS_SLES_MAX_ITERATION;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using Process-local symmetric Gauss-Seidel.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- diagonal block size
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx_ini          <-- initial system solution
 *                       (vx if nonzero, nullptr if zero)
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
                cs_real_t                 *restrict vx_ini,
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

  assert(c->setup_data != nullptr);

  /* Check for ordered variant */

  const cs_lnum_t  *order = nullptr;

  if (c->add_data != nullptr)
    order = c->add_data->order;

  if (order != nullptr)
    cvg = _p_ordered_gauss_seidel_msr(c,
                                      a,
                                      diag_block_size,
                                      convergence,
                                      rhs,
                                      vx_ini,
                                      vx);

  else
    cvg = _p_gauss_seidel_msr(c,
                              a,
                              diag_block_size,
                              convergence,
                              rhs,
                              vx_ini,
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
 * \param[in]  smoother_type   type of smoother (CG, Jacobi, ...)
 * \param[in]  poly_degree     preconditioning polynomial degree
 *                             (0: diagonal; -1: non-preconditioned;
 *                             see \ref sles_it for details)
 * \param[in]  n_iter          number of iterations to perform
 *
 * \return a pointer to newly created smoother info object, or nullptr
 *         if not available for this solver type.
 */
/*----------------------------------------------------------------------------*/

cs_sles_it_t *
cs_multigrid_smoother_create(cs_sles_it_type_t      smoother_type,
                             [[maybe_unused]] int   poly_degree,
                             int                    n_iter)
{
  switch (smoother_type) {      /* Valid choices */

  case CS_SLES_JACOBI:
    [[fallthrough]];
  case CS_SLES_P_GAUSS_SEIDEL:
    [[fallthrough]];
  case CS_SLES_P_SYM_GAUSS_SEIDEL:
    [[fallthrough]];
  case CS_SLES_L1_JACOBI:
    [[fallthrough]];
  case CS_SLES_R_JACOBI:
    [[fallthrough]];
  case CS_SLES_RJ2:
    [[fallthrough]];
  case CS_SLES_RJ3:
    [[fallthrough]];
  case CS_SLES_TS_F_GAUSS_SEIDEL:
    [[fallthrough]];
  case CS_SLES_TS_B_GAUSS_SEIDEL:
    break;

  default: /* Other iterative solvers are not tuned for smoothing */
    return nullptr;
    break;

  } /* Smoother type */

  cs_sles_it_t *c;

  CS_MALLOC(c, 1, cs_sles_it_t);

  c->solve = nullptr;
  c->_pc = nullptr;

  /* Predefined settings */
  c->type = smoother_type;
  c->on_device = false;
  c->update_stats = false;
  c->ignore_convergence = true;

  c->fallback_cvg = CS_SLES_DIVERGED;
  c->fallback_n_max_iter = 0;
  c->fallback = nullptr;

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
  c->plot = nullptr;
  c->_plot = nullptr;

#if defined(HAVE_MPI)
  c->comm = cs_glob_mpi_comm;
  c->caller_comm = cs_glob_mpi_comm;
  c->caller_n_ranks = cs_glob_n_ranks;
  if (c->caller_n_ranks < 2) {
    c->comm = MPI_COMM_NULL;
    c->caller_comm = cs_glob_mpi_comm;
  }
#endif

  c->setup_data = nullptr;
  c->add_data = nullptr;
  c->shared = nullptr;

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
  cs_sles_it_t  *c = (cs_sles_it_t *)context;

#if defined(HAVE_ACCEL)
  bool on_device = (cs_matrix_get_alloc_mode(a) > CS_ALLOC_HOST) ?
    true : false;
#endif

  const int diag_block_size = cs_matrix_get_diag_block_size(a);

  if (verbosity > 1) {
    bft_printf(_("\n Setup of solver for linear system \"%s\"\n"),
               name);
    cs_matrix_log_info(a, verbosity);
  }

  bool block_nn_inverse = false;
  bool l1_inverse = false;

  if (   c->type == CS_SLES_JACOBI
      || (   c->type >= CS_SLES_P_GAUSS_SEIDEL
          && c->type <= CS_SLES_P_SYM_GAUSS_SEIDEL)) {
    /* Force to Jacobi in case matrix type is not adapted */
    if (cs_matrix_get_type(a) != CS_MATRIX_MSR) {
      c->type = CS_SLES_JACOBI;
    }
    block_nn_inverse = true;
  }

  else if (   c->type >= CS_SLES_L1_JACOBI
           && c->type <= CS_SLES_TS_B_GAUSS_SEIDEL) {
    /* Force to closest Jacobi type in case matrix type is not adapted */
    if (   c->type == CS_SLES_TS_B_GAUSS_SEIDEL
        || c->type <= CS_SLES_TS_B_GAUSS_SEIDEL) {
      if (cs_matrix_get_type(a) != CS_MATRIX_MSR) {
        c->type = CS_SLES_JACOBI;
        c->n_max_iter = 2;
      }
    }
    block_nn_inverse = true;
  }

  switch (c->type) {

  case CS_SLES_L1_JACOBI:
    if (diag_block_size == 1) {
      c->solve = _l1_jacobi;
      l1_inverse = true;
    }
    else
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Setup of linear equation on \"%s\"\n"
                  "L1 jacobi not available for block diagonal size > 1."),
                __func__, name);
    break;

  case CS_SLES_R_JACOBI:
    [[fallthrough]];
  case CS_SLES_RJ2:
    [[fallthrough]];
  case CS_SLES_RJ3:
    {
      if (diag_block_size == 1)
        c->solve = _relaxed_jacobi;
      else
        c->solve = _block_relaxed_jacobi;
      unsigned wk_m = 1;
      if (c->type == CS_SLES_RJ2)
        wk_m = 2;
      else if (c->type == CS_SLES_RJ3)
        wk_m = 3;
      if (c->n_max_iter % wk_m)
        c->n_max_iter += wk_m - (c->n_max_iter%wk_m);
    }
#if defined(HAVE_ACCEL)
    c->on_device = on_device;
#endif
    break;

  case CS_SLES_JACOBI:
    if (diag_block_size == 1) {
      c->solve = _jacobi;
#if defined(HAVE_CUDA)
      if (on_device) {
        c->solve = cs_sles_it_cuda_jacobi;
      }
#elif defined(HAVE_HIP)
      if (on_device) {
        c->solve = cs_sles_it_hip_jacobi;
      }
#endif
    }
    else if (diag_block_size == 3) {
      c->solve = _block_3_jacobi;
#if defined(HAVE_CUDA)
      if (on_device) {
        c->solve = cs_sles_it_cuda_block_jacobi;
      }
#elif defined(HAVE_HIP)
      if (on_device) {
        c->solve = cs_sles_it_hip_block_jacobi;
      }
#endif
    }
    else {
      c->solve = _block_jacobi;
#if defined(HAVE_CUDA)
      if (on_device) {
        c->solve = cs_sles_it_cuda_block_jacobi;
      }
#elif defined(HAVE_HIP)
      if (on_device) {
        c->solve = cs_sles_it_hip_block_jacobi;
      }
#endif
    }
#if defined(HAVE_ACCEL)
    c->on_device = on_device;
#endif
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

  /* Setup preconditioner and/or auxiliary data */

  cs_sles_it_setup_priv(c, name, a, verbosity, diag_block_size,
                        block_nn_inverse, l1_inverse);

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
 * \param[in]       vx_ini         initial system solution
 *                                 (vx if nonzero, nullptr if zero)
 * \param[in, out]  vx             system solution
 * \param[in]       aux_size       number of elements in aux_vectors (in bytes)
 * \param           aux_vectors    optional working area
 *                                 (internal allocation if nullptr)
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
                            cs_real_t           *vx_ini,
                            cs_real_t           *vx,
                            size_t               aux_size,
                            void                *aux_vectors)
{
  cs_sles_it_t  *c = (cs_sles_it_t *)context;

  cs_sles_convergence_state_t cvg = CS_SLES_ITERATING;

  cs_sles_it_convergence_t  convergence;

  const cs_lnum_t diag_block_size = cs_matrix_get_diag_block_size(a);

  /* Initialize number of iterations and residual,
     and smooth linear system */

  *n_iter = 0;
  *residual = -1.0; /* Don't use this quantity when dealing with smoothers */

  /* Setup if not already done */

  if (c->setup_data == nullptr)
    cs_sles_it_setup(c, name, a, verbosity);

  if (c->pc != nullptr)
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

  bool local_solve = true;
#if defined(HAVE_MPI)
  if (c->comm == MPI_COMM_NULL) {
    cs_lnum_t n_rows = cs_matrix_get_n_rows(a);
    if (n_rows == 0) {
      local_solve = false;
    }
  }
#endif

  if (local_solve) {
    cvg = c->solve(c,
                   a, diag_block_size, &convergence,
                   rhs, vx_ini, vx,
                   aux_size, aux_vectors);
  }

  /* Update return values */

  *n_iter = convergence.n_iterations;
  *residual = convergence.residual;

  return cvg;
}

/*----------------------------------------------------------------------------*/
