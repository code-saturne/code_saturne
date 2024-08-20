/*============================================================================
 * Sparse Linear Equation Solvers using CUDA
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <algorithm>
#include <list>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

#if defined(HAVE_CUBLAS)
#include <cublas_v2.h>
#endif

#include <cooperative_groups.h>
#if (CUDART_VERSION >= 11000)
#include <cooperative_groups/reduce.h>
#endif
namespace cg = cooperative_groups;

/* Use graph capture ? */

#if (CUDART_VERSION > 9020)
#define HAVE_GRAPH_CAPTURE 1
#else
#define HAVE_GRAPH_CAPTURE 0
#endif

#include <cfloat>

// #include "mpi-ext.h"

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "cs_base.h"
#include "cs_base_accel.h"
#include "cs_blas_cuda.h"
#include "cs_log.h"
#include "cs_matrix.h"
#include "cs_matrix_priv.h"
#include "cs_matrix_spmv_cuda.h"
#include "cs_sles_it.h"
#include "cs_sles_pc.h"
#include "cs_sles_it_priv.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_sles_it_cuda.h"

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

 /* SIMD unit size to ensure SIMD alignement (based on warp size) */

#define CS_SIMD_SIZE(s) (((s-1)/32+1)*32)
#define CS_BLOCKSIZE 256

/*----------------------------------------------------------------------------
 * Compatibility macro for __ldg (load from generic memory) intrinsic,
 * forcing load from read-only texture cache.
 *
 * This was not available in (very old) CUDA architectures.
 *----------------------------------------------------------------------------*/

#if __CUDA_ARCH__ < 350
#define __ldg(ptr) *(ptr);
#endif

/*============================================================================
 *  Global variables
 *============================================================================*/

static bool _use_cublas = false;

/*============================================================================
 * Private function and kernel definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Kernel for sum reduction within a warp (for warp size 32).
 *
 * \param[in, out]  stmp  shared value to reduce
 * \param[in, out]  tid   thread id
 */
/*----------------------------------------------------------------------------*/

template <size_t block_size, typename T>
__device__ static void __forceinline__
_warp_reduce_sum(volatile T  *stmp,
                 size_t       tid)
{
  if (block_size >= 64) stmp[tid] += stmp[tid + 32];
  if (block_size >= 32) stmp[tid] += stmp[tid + 16];
  if (block_size >= 16) stmp[tid] += stmp[tid +  8];
  if (block_size >=  8) stmp[tid] += stmp[tid +  4];
  if (block_size >=  4) stmp[tid] += stmp[tid +  2];
  if (block_size >=  2) stmp[tid] += stmp[tid +  1];
}

/*----------------------------------------------------------------------------
 * Compute Vx <- Vx - (A-diag).Rk for Jacobi solver.
 *
 * parameters:
 *   n_rows    <-- number of rows
 *   ad_inv    <-- inverse of diagonal
 *   rhs       <-- right hand side
 *   vx        <-> solution at current iteration
 *   rk        <-> solution at previous iteraton
 *----------------------------------------------------------------------------*/

template <size_t block_size>
__global__ static void
_jacobi_compute_vx(cs_lnum_t                       n_rows,
                   const cs_real_t   *__restrict__ ad_inv,
                   const cs_real_t   *__restrict__ rhs,
                   cs_real_t         *__restrict__ vx,
                   cs_real_t         *__restrict__ rk)
{
  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;

  if (ii < n_rows) {
    vx[ii] = (rhs[ii]-vx[ii])*ad_inv[ii];
    rk[ii] = vx[ii];
  }
}

/*----------------------------------------------------------------------------
 * Compute Vx <- Vx - (A-diag).Rk for Jacobi solver when initial Vx is zero.
 *
 * parameters:
 *   n_rows    <-- number of rows
 *   ad_inv    <-- inverse of diagonal
 *   ad        <-- diagonal
 *   rhs       <-- right hand side
 *   vx        <-> solution at current iteration
 *   rk        <-> solution at previous iteraton
 *----------------------------------------------------------------------------*/

template <size_t block_size>
__global__ static void
_jacobi_compute_vx_ini0(cs_lnum_t                       n_rows,
                        const cs_real_t   *__restrict__ ad_inv,
                        const cs_real_t   *__restrict__ rhs,
                        cs_real_t         *__restrict__ vx,
                        cs_real_t         *__restrict__ rk)
{
  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;

  if (ii < n_rows) {
    vx[ii] = rhs[ii]*ad_inv[ii];
    rk[ii] = vx[ii];
  }
}

/*----------------------------------------------------------------------------
 * Compute Vx <- Vx - (A-diag).Rk and residual for Jacobi solver.
 *
 * This call must be followed by
 * cs_blas_cuda_reduce_single_block<block_size><<<1, block_size, 0>>>
 *  (grid_size, sum_block, _r_reduce);
 * Also, block_size must be a power of 2.
 *
 * parameters:
 *   n_rows    <-- number of rows
 *   ad_inv    <-- inverse of diagonal
 *   ad        <-- diagonal
 *   rhs       <-- right hand side
 *   vx        <-> solution at current iteration
 *   rk        <-> solution at previous iteraton
 *   sum_block --> contribution to residual
 *----------------------------------------------------------------------------*/

template <size_t block_size>
__global__ static void
_jacobi_compute_vx_and_residual(cs_lnum_t                       n_rows,
                                const cs_real_t   *__restrict__ ad_inv,
                                const cs_real_t   *__restrict__ ad,
                                const cs_real_t   *__restrict__ rhs,
                                cs_real_t         *__restrict__ vx,
                                cs_real_t         *__restrict__ rk,
                                double            *__restrict__ sum_block)
{
  __shared__ double sdata[block_size];

  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
  size_t tid = threadIdx.x;

  if (ii < n_rows) {
    vx[ii] = (rhs[ii]-vx[ii])*ad_inv[ii];
    double r = ad[ii] * (vx[ii]-rk[ii]);
    sdata[tid] = r*r;
    rk[ii] = vx[ii];
  }
  else
    sdata[tid] = 0.0;

  cs_blas_cuda_block_reduce_sum<block_size, 1>(sdata, tid, sum_block);
}

/*----------------------------------------------------------------------------
 * Compute Vx <- Vx - (A-diag).Rk and residual for Jacobi solver
 * when the initial value of Vx is zero.
 *
 * This call must be followed by
 * cs_blas_cuda_reduce_single_block<block_size><<<1, block_size, 0>>>
 *  (grid_size, sum_block, _r_reduce);
 * Also, block_size must be a power of 2.
 *
 * parameters:
 *   n_rows    <-- number of rows
 *   ad_inv    <-- inverse of diagonal
 *   rhs       <-- right hand side
 *   vx        <-> solution at current iteration
 *   rk        <-> solution at previous iteraton
 *   sum_block --> contribution to residual
 *----------------------------------------------------------------------------*/

template <size_t block_size>
__global__ static void
_jacobi_compute_vx_and_residual_ini0(cs_lnum_t                       n_rows,
                                     const cs_real_t   *__restrict__ ad_inv,
                                     const cs_real_t   *__restrict__ rhs,
                                     cs_real_t         *__restrict__ vx,
                                     cs_real_t         *__restrict__ rk,
                                     double            *__restrict__ sum_block)
{
  __shared__ double sdata[block_size];

  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
  size_t tid = threadIdx.x;

  if (ii < n_rows) {
    vx[ii] = rhs[ii]*ad_inv[ii];
    double r = rhs[ii];
    sdata[tid] = r*r;
    rk[ii] = vx[ii];
  }
  else
    sdata[tid] = 0.0;

  cs_blas_cuda_block_reduce_sum<block_size, 1>(sdata, tid, sum_block);
}

/*----------------------------------------------------------------------------
 * Block Jacobi utilities.
 * Compute forward and backward to solve an LU 3*3 system.
 *
 * parameters:
 *   mat   <-- 3*3*dim matrix
 *   x     <-> 1st part of RHS (c - b) in, solution out
 *   c     --> 2nd part of RHS (c - b)
 *----------------------------------------------------------------------------*/

__device__ static void  __forceinline__
_fw_and_bw_lu33_cuda(const cs_real_t  *__restrict__ mat,
                     cs_real_t        *__restrict__ x,
                     const cs_real_t  *__restrict__ c)
{
  cs_real_t  aux[3];

  aux[0] = (c[0] - x[0]);
  aux[1] = (c[1] - x[1]) - aux[0]*mat[3];
  aux[2] = (c[2] - x[2]) - aux[0]*mat[6] - aux[1]*mat[7];

  x[2] = aux[2]/mat[8];
  x[1] = (aux[1] - mat[5]*x[2])/mat[4];
  x[0] = (aux[0] - mat[1]*x[1] - mat[2]*x[2])/mat[0];
}

/*----------------------------------------------------------------------------
 * Compute Vx <- Vx - (A-diag).Rk and residual for Jacobi with
 * 3x3 block-diagonal matrix.
 *
 * This call must be followed by
 * cs_blas_cuda_reduce_single_block<block_size><<<1, block_size, 0>>>
 *  (grid_size, sum_block, _r_reduce);
 * Also, block_size must be a power of 2.
 *
 * parameters:
 *   n_b_rows  <-- number of block rows
 *   ad_inv    <-- inverse of diagonal
 *   ad        <-- diagonal
 *   rhs       <-- right hand side
 *   vx        <-> 1st part of RHS (c - b) in, solution at current iteration
 *   rk        <-> solution at previous iteraton
 *   sum_block --> contribution to residual
 *----------------------------------------------------------------------------*/

template <size_t block_size>
__global__ static void
_block_3_jacobi_compute_vx_and_residual
(
 cs_lnum_t                       n_b_rows,
 const cs_real_t   *__restrict__ ad_inv,
 const cs_real_t   *__restrict__ ad,
 const cs_real_t   *__restrict__ rhs,
 cs_real_t         *__restrict__ vx,
 cs_real_t         *__restrict__ rk,
 double            *__restrict__ sum_block
)
{
  __shared__ cs_real_t sdata[block_size];

  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
  size_t tid = threadIdx.x;
  sdata[tid] = 0.0;

  if (ii < n_b_rows) {
    _fw_and_bw_lu33_cuda(ad_inv + 9*ii, vx + 3*ii, rhs + 3*ii);

    #pragma unroll
    for (cs_lnum_t jj = 0; jj < 3; jj++) {
      double r = 0.0;
      #pragma unroll
      for (cs_lnum_t kk = 0; kk < 3; kk++)
        r +=  ad[ii*9 + jj*3 + kk] * (vx[ii*3 + kk] - rk[ii*3 + kk]);

      sdata[tid] += (r*r);
    }

    #pragma unroll
    for (cs_lnum_t kk = 0; kk < 3; kk++)
      rk[ii*3 + kk] = vx[ii*3 + kk];
  }

  cs_blas_cuda_block_reduce_sum<block_size, 1>(sdata, tid, sum_block);
}

/*----------------------------------------------------------------------------
 * Compute Vx <- Vx - (A-diag).Rk and residual for Jacobi with
 * 3x3 block-diagonal matrix when Vx is initially 0.
 *
 * This call must be followed by
 * cs_blas_cuda_reduce_single_block<block_size><<<1, block_size, 0>>>
 *  (grid_size, sum_block, _r_reduce);
 * Also, block_size must be a power of 2.
 *
 * parameters:
 *   n_b_rows  <-- number of block rows
 *   ad_inv    <-- inverse of diagonal
 *   ad        <-- diagonal
 *   rhs       <-- right hand side
 *   vx        <-> 1st part of RHS (c - b) in, solution at current iteration
 *   rk        <-> solution at previous iteraton
 *   sum_block --> contribution to residual
 *----------------------------------------------------------------------------*/

template <size_t block_size>
__global__ static void
_block_3_jacobi_compute_vx_and_residual_ini0
(
 cs_lnum_t                       n_b_rows,
 const cs_real_t   *__restrict__ ad_inv,
 const cs_real_t   *__restrict__ ad,
 const cs_real_t   *__restrict__ rhs,
 cs_real_t         *__restrict__ vx,
 cs_real_t         *__restrict__ rk,
 double            *__restrict__ sum_block
)
{
  __shared__ cs_real_t sdata[block_size];

  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
  size_t tid = threadIdx.x;
  sdata[tid] = 0.0;

  if (ii < n_b_rows) {
    #pragma unroll
    for (cs_lnum_t kk = 0; kk < 3; kk++) {
      vx[ii*3 + kk] = 0;
    }

    _fw_and_bw_lu33_cuda(ad_inv + 9*ii, vx + 3*ii, rhs + 3*ii);

    #pragma unroll
    for (cs_lnum_t jj = 0; jj < 3; jj++) {
      double r = 0.0;
      #pragma unroll
      for (cs_lnum_t kk = 0; kk < 3; kk++)
        r +=  ad[ii*9 + jj*3 + kk] * vx[ii*3 + kk];

      sdata[tid] += (r*r);
    }

    #pragma unroll
    for (cs_lnum_t kk = 0; kk < 3; kk++)
      rk[ii*3 + kk] = vx[ii*3 + kk];
  }

  cs_blas_cuda_block_reduce_sum<block_size, 1>(sdata, tid, sum_block);
}

/*----------------------------------------------------------------------------
 * Block Jacobi utilities.
 * Compute forward and backward to solve an LU system.
 *
 * parameters:
 *   mat     <-- p*p*dim matrix
 *   db_size <-- block (p) size
 *   x       <-> 1st part of RHS (c - b) in, solution out
 *   c       --> 2nd part of RHS (c - b)
 *----------------------------------------------------------------------------*/

__device__ static void  __forceinline__
_fw_and_bw_lu_cuda(const cs_real_t  *__restrict__ mat,
                   cs_lnum_t                      db_size,
                   cs_real_t        *__restrict__ x,
                   const cs_real_t  *__restrict__ c)
{
  assert(db_size <= 9);
  cs_real_t aux[9];

  /* forward */
  for (cs_lnum_t ii = 0; ii < db_size; ii++) {
    aux[ii] = (c[ii] - x[ii]);
    for (cs_lnum_t jj = 0; jj < ii; jj++) {
      aux[ii] -= aux[jj]*mat[ii*db_size + jj];
    }
  }

  /* backward */
  for (cs_lnum_t ii = db_size - 1; ii >= 0; ii-=1) {
    x[ii] = aux[ii];
    for (cs_lnum_t jj = db_size - 1; jj > ii; jj-=1) {
      x[ii] -= x[jj]*mat[ii*db_size + jj];
    }
    x[ii] /= mat[ii*(db_size + 1)];
  }
}

/*----------------------------------------------------------------------------
 * Compute Vx <- Vx - (A-diag).Rk and residual for Jacobi with
 * 3x3 block-diagonal matrix.
 *
 * parameters:
 *   n_b_rows  <-- number of block rows
 *   ad_inv    <-- inverse of diagonal
 *   ad        <-- diagonal
 *   rhs       <-- right hand side
 *   vx        <-> 1st part of RHS (c - b) in, solution at current iteration
 *   rk        <-> solution at previous iteraton
 *   sum_block --> contribution to residual
 *----------------------------------------------------------------------------*/

template <size_t block_size>
__global__ static void
_block_jacobi_compute_vx_and_residual(cs_lnum_t                       n_b_rows,
                                      cs_lnum_t                       db_size,
                                      const cs_real_t   *__restrict__ ad_inv,
                                      const cs_real_t   *__restrict__ ad,
                                      const cs_real_t   *__restrict__ rhs,
                                      cs_real_t         *__restrict__ vx,
                                      cs_real_t         *__restrict__ rk,
                                      double            *__restrict__ sum_block)
{
  __shared__ cs_real_t sdata[block_size];

  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
  size_t tid = threadIdx.x;
  sdata[tid] = 0.0;

  cs_lnum_t n = db_size;
  cs_lnum_t nn = db_size*db_size;

  if (ii < n_b_rows) {
    _fw_and_bw_lu_cuda(ad_inv + nn*ii, n, vx + n*ii, rhs + n*ii);

    #pragma unroll
    for (cs_lnum_t jj = 0; jj < n; jj++) {
      double r = 0.0;
      #pragma unroll
      for (cs_lnum_t kk = 0; kk < n; kk++)
        r +=  ad[ii*nn + jj*n + kk] * (vx[ii*n + kk] - rk[ii*n + kk]);

      sdata[tid] += (r*r);
    }

    #pragma unroll
    for (cs_lnum_t kk = 0; kk < n; kk++)
      rk[ii*n + kk] = vx[ii*n + kk];
  }

  cs_blas_cuda_block_reduce_sum<block_size, 1>(sdata, tid, sum_block);
}

/*----------------------------------------------------------------------------
 * Compute Vx <- Vx - (A-diag).Rk and residual for Jacobi with
 * block-diagonal matrix when Vx is initially 0.
 *
 * parameters:
 *   n_b_rows  <-- number of block rows
 *   ad_inv    <-- inverse of diagonal
 *   ad        <-- diagonal
 *   rhs       <-- right hand side
 *   vx        <-> 1st part of RHS (c - b) in, solution at current iteration
 *   rk        <-> solution at previous iteraton
 *   sum_block --> contribution to residual
 *----------------------------------------------------------------------------*/

template <size_t block_size>
__global__ static void
_block_jacobi_compute_vx_and_residual_ini0
(
 cs_lnum_t                       n_b_rows,
 cs_lnum_t                       db_size,
 const cs_real_t   *__restrict__ ad_inv,
 const cs_real_t   *__restrict__ ad,
 const cs_real_t   *__restrict__ rhs,
 cs_real_t         *__restrict__ vx,
 cs_real_t         *__restrict__ rk,
 double            *__restrict__ sum_block
)
{
  __shared__ cs_real_t sdata[block_size];

  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
  size_t tid = threadIdx.x;
  sdata[tid] = 0.0;

  cs_lnum_t n = db_size;
  cs_lnum_t nn = db_size*db_size;

  if (ii < n_b_rows) {
    #pragma unroll
    for (cs_lnum_t kk = 0; kk < n; kk++) {
      vx[ii*n + kk] = 0;
    }

    _fw_and_bw_lu_cuda(ad_inv + nn*ii, n, vx + n*ii, rhs + n*ii);

    #pragma unroll
    for (cs_lnum_t jj = 0; jj < n; jj++) {
      double r = 0.0;
      #pragma unroll
      for (cs_lnum_t kk = 0; kk < n; kk++)
        r +=  ad[ii*nn + jj*n + kk] * vx[ii*n + kk];

      sdata[tid] += (r*r);
    }

    #pragma unroll
    for (cs_lnum_t kk = 0; kk < n; kk++)
      rk[ii*n + kk] = vx[ii*n + kk];
  }

  cs_blas_cuda_block_reduce_sum<block_size, 1>(sdata, tid, sum_block);
}

/*----------------------------------------------------------------------------
 * FCG initialization: rk <- rhs -rk; qk = 0.
 *
 * parameters:
 *   n         <-- number of elements
 *   rhs       <-- vector of elements
 *   rk        <-> vector of elements
 *   qk        --> vector of elements
 *----------------------------------------------------------------------------*/

__global__ static void
_fcg_init(cs_lnum_t                      n,
          const cs_real_t  *__restrict__ rhs,
          cs_real_t        *__restrict__ rk,
          cs_real_t        *__restrict__ qk)
{
  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;

  if (ii < n) {
    rk[ii] = rhs[ii] - rk[ii];
    qk[ii] = 0.;
  }
}

/*----------------------------------------------------------------------------
 * FCG initialization when inital solution is zero: vx <- 0; rk <- rhs; qk = 0.
 *
 * parameters:
 *   n         <-- number of elements
 *   rhs       <-- vector of elements
 *   vx        --> vector of elements
 *   rk        --> vector of elements
 *   qk        --> vector of elements
 *----------------------------------------------------------------------------*/

__global__ static void
_fcg_init_vx0(cs_lnum_t                      n,
              const cs_real_t  *__restrict__ rhs,
              cs_real_t        *__restrict__ vx,
              cs_real_t        *__restrict__ rk,
              cs_real_t        *__restrict__ qk)
{
  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;

  if (ii < n) {
    vx[ii] = 0.;
    rk[ii] = rhs[ii];
    qk[ii] = 0.;
  }
}

/*----------------------------------------------------------------------------
 * FCG update of dk, vk, qk, rk at iterations > 0.
 *
 * parameters:
 *   n      <-- number of elements
 *   gk_rk1 <-- descent parameter
 *   ak_rk  <-- descent parameter
 *   wk     <-- vector of elements
 *   dk     <-> vector of elements
 *   qk     <-- vector of elements
 *   rk     <-- vector of elements
 *   vk     <-> vector of elements
 *   vx     <-> vector of elements
 *----------------------------------------------------------------------------*/

__global__ static void
_fcg_update_n(cs_lnum_t                      n,
              cs_real_t                      gk_rk1,
              cs_real_t                      ak_rk,
              const cs_real_t  *__restrict__ wk,
              cs_real_t        *__restrict__ dk,
              cs_real_t        *__restrict__ qk,
              cs_real_t        *__restrict__ rk,
              cs_real_t        *__restrict__ vk,
              cs_real_t        *__restrict__ vx)
{
  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;

  if (ii < n) {
    dk[ii] = vk[ii] - gk_rk1 * dk[ii];
    vx[ii] = vx[ii] + ak_rk * dk[ii];

    qk[ii] = wk[ii] - gk_rk1 * qk[ii];
    rk[ii] = rk[ii] - ak_rk * qk[ii];
  }
}

/*----------------------------------------------------------------------------
 * FCG update of dk, vk, qk, rk at iteration 0.
 *
 * parameters:
 *   n      <-- number of elements
 *   ak_rk  <-- descent parameter
 *   wk     <-- vector of elements
 *   dk     <-> vector of elements
 *   qk     <-- vector of elements
 *   rk     <-- vector of elements
 *   vk     <-> vector of elements
 *   vx     <-> vector of elements
 *----------------------------------------------------------------------------*/

__global__ static void
_fcg_update_0(cs_lnum_t                      n,
              cs_real_t                      ak_rk,
              const cs_real_t  *__restrict__ wk,
              cs_real_t        *__restrict__ dk,
              cs_real_t        *__restrict__ qk,
              cs_real_t        *__restrict__ rk,
              cs_real_t        *__restrict__ vk,
              cs_real_t        *__restrict__ vx)
{
  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;

  if (ii < n) {
    dk[ii] = vk[ii];
    vx[ii] = vx[ii] + ak_rk * vk[ii];

    qk[ii] = wk[ii];
    rk[ii] = rk[ii] - ak_rk * wk[ii];
  }
}

/*----------------------------------------------------------------------------
 * Compute y <- y - x and stage 1 of resulting y.y.
 *
 * This call must be followed by
 * cs_blas_cuda_reduce_single_block<block_size><<<1, block_size, 0>>>
 *  (grid_size, sum_block, _r_reduce);
 * Also, block_size must be a power of 2.
 *
 * parameters:
 *   n         <-- number of elements
 *   x         <-- vector of elements
 *   y         <-> vector of elements
 *   sum_block --> contribution to residual
 *----------------------------------------------------------------------------*/

template <size_t block_size>
__global__ static void
_ymx_dot_yy(cs_lnum_t                       n,
            const cs_real_t   *__restrict__ x,
            cs_real_t         *__restrict__ y,
            double            *__restrict__ sum_block)
{
  __shared__ double sdata[block_size];

  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
  size_t tid = threadIdx.x;
  size_t grid_size = blockDim.x*gridDim.x;

  sdata[tid] = 0.0;

  while (ii < n) {
    y[ii] = y[ii] - x[ii];
    double r = y[ii];
    sdata[tid] += r*r;

    ii += grid_size;
  }

  cs_blas_cuda_block_reduce_sum<block_size, 1>(sdata, tid, sum_block);
}

/*----------------------------------------------------------------------------
 * Set x to 0 and y to -r and compute stage 1 of resulting y.y.
 *
 * This call must be followed by
 * cs_blas_cuda_reduce_single_block<block_size><<<1, block_size, 0>>>
 *  (grid_size, sum_block, _r_reduce);
 * Also, block_size must be a power 2.
 *
 * parameters:
 *   n         <-- number of elements
 *   r         <-- vector of elements
 *   x         <-- vector of elements
 *   y         <-> vector of elements
 *   sum_block --> contribution to residual
 *----------------------------------------------------------------------------*/

template <size_t block_size>
__global__ static void
_x0_ymr_dot_yy(cs_lnum_t                       n,
               const cs_real_t   *__restrict__ r,
               cs_real_t         *__restrict__ x,
               cs_real_t         *__restrict__ y,
               double            *__restrict__ sum_block)
{
  __shared__ double sdata[block_size];

  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
  size_t tid = threadIdx.x;
  size_t grid_size = blockDim.x*gridDim.x;

  sdata[tid] = 0.0;

  while (ii < n) {
    x[ii] = 0;
    y[ii] = -r[ii];
    double rii = r[ii];
    sdata[tid] += rii*rii;

    ii += grid_size;
  }

  cs_blas_cuda_block_reduce_sum<block_size, 1>(sdata, tid, sum_block);
}

/*----------------------------------------------------------------------------
 * Compute y <- -alpha.x + y and stage 1 of resulting y.y.
 *
 * This call must be followed by
 * cs_blas_cuda_reduce_single_block<block_size><<<1, block_size, 0>>>
 *  (grid_size, sum_block, _r_reduce);
 *
 * parameters:
 *   n         <-- number of elements
 *   x         <-- vector of elements
 *   y         <-> vector of elements
 *   sum_block --> contribution to residual
 *----------------------------------------------------------------------------*/

template <size_t block_size>
__global__ static void
_smaxpy_dot_yy(cs_lnum_t                       n,
               const cs_real_t   *__restrict__ alpha,
               const cs_real_t   *__restrict__ x,
               cs_real_t         *__restrict__ y,
               double            *__restrict__ sum_block)
{
  __shared__ double sdata[block_size];

  cs_real_t _alpha = - *alpha;

  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
  size_t tid = threadIdx.x;
  size_t grid_size = blockDim.x*gridDim.x;

  sdata[tid] = 0.0;

  while (ii < n) {
    cs_real_t r = _alpha*x[ii] + y[ii];
    y[ii] = r;
    sdata[tid] += r*r;

    ii += grid_size;
  }

  cs_blas_cuda_block_reduce_sum<block_size, 1>(sdata, tid, sum_block);
}

/*----------------------------------------------------------------------------
 * Compute y <- alpha.y and stage 1 of following x.y dot product.
 *
 * This call must be followed by
 * cs_blas_cuda_reduce_single_block<block_size><<<1, block_size, 0>>>
 *  (grid_size, sum_block, _r_reduce);
 *
 * parameters:
 *   n         <-- number of elements
 *   alpha     <-- scaling factor
 *   x         <-- vector of elements
 *   y         <-> vector of elements
 *   sum_block --> contribution to residual
 *----------------------------------------------------------------------------*/

template <size_t block_size>
__global__ static void
_y_scale_dot_xy(cs_lnum_t                      n,
                const cs_real_t  *__restrict__ alpha,
                const cs_real_t  *__restrict__ x,
                cs_real_t        *__restrict__ y,
                double           *__restrict__ sum_block)
{
  __shared__ double sdata[block_size];

  cs_real_t _alpha = *alpha;

  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
  size_t tid = threadIdx.x;
  size_t grid_size = blockDim.x*gridDim.x;

  sdata[tid] = 0.0;

  while (ii < n) {
    cs_real_t r = _alpha * y[ii];
    y[ii] = r;
    sdata[tid] += x[ii]*r;

    ii += grid_size;
  }

  cs_blas_cuda_block_reduce_sum<block_size, 1>(sdata, tid, sum_block);
}

/*----------------------------------------------------------------------------
 * Invert Gamma for GCR algorithm, for multiple restart directions.
 *
 * This is a basically small, serial operation, done only on one thread.
 *
 * parameters:
 *   n_c_iter <-- number of restart iterations
 *   n_gkj    <-- number of gkj terms (diagonal matrix)
 *   mgkj     <-- - Gamma
 *   gkj_inv  <-> Inverse of Gamma
 *----------------------------------------------------------------------------*/

__global__ static void
_gcr_gkj_inv(int                            n_c_iter,
             int                            n_gkj,
             const cs_real_t  *__restrict__ mgkj,
             cs_real_t        *__restrict__ gkj_inv)
{
  int kk = blockIdx.x*blockDim.x + threadIdx.x;
  size_t grid_size = blockDim.x*gridDim.x;

  while (kk < n_gkj) {
    gkj_inv[kk] = 0;
    kk += grid_size;
  }

  kk = blockIdx.x*blockDim.x + threadIdx.x;
  if (kk != 0)
    return;

  for (int kk = 0; kk < (int)n_c_iter; kk++) {
    for (int ii = 0; ii < kk; ii++) {
      for (int jj = 0; jj < kk; jj++)
        gkj_inv[(kk + 1) * kk / 2 + ii]
          -=   ((ii <= jj) ? gkj_inv[(jj + 1) * jj / 2 + ii] : 0.0)
             * mgkj[(kk + 1) * kk / 2  + jj];
    }

    for (int jj = 0; jj < kk; jj++)
      gkj_inv[(kk + 1) * kk / 2 + jj] /= mgkj[(kk + 1) * kk / 2 + kk];

    gkj_inv[(kk + 1) * kk / 2 + kk] = 1.0 / -mgkj[(kk + 1) * kk / 2 + kk];
  }
}

/*----------------------------------------------------------------------------
 * Invert Gamma for GCR algorithm, for single restart direction
 *
 * This is a basically small, serial operation, done only on one thread.
 *
 * parameters:
 *   gkj      <-- Gamma
 *   gkj_inv  <-> Inverse of Gamma
 *----------------------------------------------------------------------------*/

__global__ static void
_gcr_gkj_inv_1(const cs_real_t  *__restrict__ gkj,
               cs_real_t        *__restrict__ gkj_inv)
{
  int ii = blockIdx.x*blockDim.x + threadIdx.x;
  if (ii != 0)
    return;

  if (abs(gkj[0]) <= 0)
    gkj_inv[0] = 1.0;
  else
    gkj_inv[0] = -1.0 / gkj[0];
}

/*----------------------------------------------------------------------------
 * Update solution vector of GCR
 *
 * parameters:
 *   n_rows   <-- number of rows
 *   n_c_iter <-- number of iterations in restart block
 *   wa_size  <-- work array vector size (may include padding, >= n_rows)
 *   alpha    <-- alpha coefficients
 *   gkj_inv  <-- inverse of Gamma
 *   zk       <-- zk work vectors (preconditioned rk -> zk)
 *   vx       <-> solution vector
 *----------------------------------------------------------------------------*/

__global__ static void
_gcr_update_vx(cs_lnum_t                      n_rows,
               int                            n_c_iter,
               size_t                         wa_size,
               const cs_real_t  *__restrict__ alpha,
               const cs_real_t  *__restrict__ gkj_inv,
               const cs_real_t  *__restrict__ zk,
               cs_real_t        *__restrict__ vx)
{
  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;

  if (ii < n_rows) {
    double sii = 0.0;
    #pragma unroll(2)
    for(cs_lnum_t kk = 0; kk < n_c_iter; kk++){
      for(cs_lnum_t jj = 0; jj <= kk; jj++){
        const cs_real_t *zk_j = zk + jj*wa_size;
        sii += alpha[kk] * zk_j[ii] * gkj_inv[(kk + 1) * kk / 2 + jj];
      }
    }
    vx[ii] -= sii;
  }
}

/*----------------------------------------------------------------------------
 * Syncronize a reduction sum globally.
 *
 * The associated stream will be synchronized first, then a global
 * MPI reduction will be used if needed.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c          <-- pointer to solver context info
 *   stream     <-- associated stream
 *   tuple_size <-- number of values in reduced tuple
 *   res        <-> local sum in, globally sum out (host)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static void
_sync_reduction_sum(const cs_sles_it_t  *c,
                    cudaStream_t         stream,
                    cs_lnum_t            tuple_size,
                    double               res[])
{
  cudaStreamSynchronize(stream);

#if defined(HAVE_MPI)

  if (c->comm != MPI_COMM_NULL)
    MPI_Allreduce(MPI_IN_PLACE, res, tuple_size, MPI_DOUBLE, MPI_SUM, c->comm);

#endif /* defined(HAVE_MPI) */
}

/*----------------------------------------------------------------------------
 * Compute dot product, summing result over all participating ranks.
 *
 * The associated stream will be synchronized first, then a global
 * MPI reduction will be used if needed.
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
             cudaStream_t         stream,
             const cs_real_t     *x,
             const cs_real_t     *y)
{
  double s;

  /* Alternatives */

  if (_use_cublas == false)
    s = cs_blas_cuda_dot(c->setup_data->n_rows, x, y);

#if defined(HAVE_CUBLAS)
  if (_use_cublas)
    s = cs_blas_cublas_dot(c->setup_data->n_rows, x, y);
#endif

#if defined(HAVE_MPI)

  if (c->comm != MPI_COMM_NULL)
    MPI_Allreduce(MPI_IN_PLACE, &s, 1, MPI_DOUBLE, MPI_SUM, c->comm);

#endif /* defined(HAVE_MPI) */

  return s;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute dot product x.y, summing result over all threads of a block.
 *
 * blockSize must be a power of 2.
 *
 * \param[in]   n      array size
 * \param[in]   v      v vector
 * \param[in]   w      w vector
 * \param[in]   q      q vector
 * \param[in]   r      r vector
 * \param[out]  b_res  result of s = x.x
 */
/*----------------------------------------------------------------------------*/

template <size_t blockSize, typename T>
__global__ static void
_dot_products_vr_vw_vq_rr_stage_1_of_2(cs_lnum_t    n,
                                       const T     *v,
                                       const T     *w,
                                       const T     *q,
                                       const T     *r,
                                       double      *b_res)
{
  __shared__ double stmp[blockSize*4];

  cs_lnum_t tid = threadIdx.x;
  size_t grid_size = blockDim.x*gridDim.x;

  stmp[tid*4] = 0.;
  stmp[tid*4 + 1] = 0.;
  stmp[tid*4 + 2] = 0.;
  stmp[tid*4 + 3] = 0.;

  for (cs_lnum_t i = blockIdx.x*(blockDim.x) + tid;
       i < n;
       i += grid_size) {
    stmp[tid*4]     += static_cast<double>(v[i] * r[i]);
    stmp[tid*4 + 1] += static_cast<double>(v[i] * w[i]);
    stmp[tid*4 + 2] += static_cast<double>(v[i] * q[i]);
    stmp[tid*4 + 3] += static_cast<double>(r[i] * r[i]);
  }

  // Output: b_res for this block

  cs_blas_cuda_block_reduce_sum<blockSize, 4>(stmp, tid, b_res);
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
                          cudaStream_t         stream,
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

  /* Alternatives (need to set option for this) */

  cudaStreamSynchronize(stream);

  cs_lnum_t n = c->setup_data->n_rows;

  const unsigned int block_size = 256;
  unsigned int grid_size = cs_cuda_grid_size(n, block_size);

  if (_use_cublas == false) {
    double *sum_block, *_s;
    cs_blas_cuda_get_2_stage_reduce_buffers(n, 4, grid_size, sum_block, _s);

    _dot_products_vr_vw_vq_rr_stage_1_of_2
      <block_size><<<grid_size, block_size, 0, stream>>>
      (n, v, w, q, r, sum_block);
    cs_blas_cuda_reduce_single_block<block_size, 4><<<1, block_size, 0, stream>>>
      (grid_size, sum_block, _s);

    _sync_reduction_sum(c, stream, 4, _s);

    for (int i = 0; i < 4; i++)
      s[i] = _s[i];
  }

#if defined(HAVE_CUBLAS)

  if (_use_cublas) {

    s[0] = cs_blas_cublas_dot(c->setup_data->n_rows, v, r);
    s[1] = cs_blas_cublas_dot(c->setup_data->n_rows, v, w);
    s[2] = cs_blas_cublas_dot(c->setup_data->n_rows, v, q);
    s[3] = cs_blas_cublas_dot(c->setup_data->n_rows, r, r);

#if defined(HAVE_MPI)

    if (c->comm != MPI_COMM_NULL) {
      MPI_Allreduce(MPI_IN_PLACE, s, 4, MPI_DOUBLE, MPI_SUM, c->comm);
    }

#endif /* defined(HAVE_MPI) */

  }

#endif /* defined(HAVE_CUBLAS) */

  *s1 = s[0];
  *s2 = s[1];
  *s3 = s[2];
  *s4 = s[3];
}

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using Jacobi (CUDA version).
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- diagonal block size
 *   rotation_mode   <-- halo update option for rotational periodicity
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

cs_sles_convergence_state_t
cs_sles_it_cuda_jacobi(cs_sles_it_t              *c,
                       const cs_matrix_t         *a,
                       cs_lnum_t                  diag_block_size,
                       cs_sles_it_convergence_t  *convergence,
                       const cs_real_t           *rhs,
                       cs_real_t                 *vx_ini,
                       cs_real_t                 *vx,
                       size_t                     aux_size,
                       void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg= CS_SLES_ITERATING;
  unsigned n_iter = 0;

  int device_id = cs_get_device_id();

  bool local_stream = false;
  cudaStream_t stream;
  stream = cs_matrix_spmv_cuda_get_stream();
  if (stream == 0) {
    local_stream = true;
    cudaStreamCreate(&stream);
  }

  const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;

  size_t vec_size = n_cols * sizeof(cs_real_t);

  /* Prefetch in case it is needed; actually, the reported allocation
     mode may be incorrect if the array are sub-arrays of a greater allocation,
     such as for multigrid, but in this case no prefetching should be needed,
     or it should have bed done by the caller. */
  {
    cs_alloc_mode_t amode_vx = cs_check_device_ptr(vx);
    cs_alloc_mode_t amode_rhs = cs_check_device_ptr(rhs);

    if (amode_vx == CS_ALLOC_HOST_DEVICE_SHARED && vx_ini == vx)
      cs_cuda_prefetch_h2d(vx, vec_size);

    if (amode_rhs == CS_ALLOC_HOST_DEVICE_SHARED)
      cs_cuda_prefetch_h2d(rhs, vec_size);
  }

  const cs_real_t  *__restrict__ ad
    =  cs_get_device_ptr_const(cs_matrix_get_diagonal(a));
  const cs_real_t *__restrict__ ad_inv
    = cs_get_device_ptr_const(c->setup_data->ad_inv);

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  double residual = -1.;

  if (convergence->precision < 0);
    residual = convergence->precision * convergence->r_norm * 2;

  /* Allocate or map work arrays
     --------------------------- */

  assert(c->setup_data != NULL);

  cs_real_t *_aux_vectors = NULL;
  if (n_cols > 0) {
    const size_t n_wa = 1;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (   aux_vectors == NULL
        || cs_cuda_is_device_ptr(aux_vectors) == false
        || aux_size/sizeof(cs_real_t) < (wa_size * n_wa)) {
#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
      cudaMalloc(&_aux_vectors, wa_size * n_wa *sizeof(cs_real_t));
#else
      cudaMallocManaged(&_aux_vectors, wa_size * n_wa *sizeof(cs_real_t));
#endif
    }
    else
      _aux_vectors = (cs_real_t *)aux_vectors;
  }

  cs_real_t *__restrict__ rk = _aux_vectors;

  const unsigned int blocksize = 256;
  unsigned int gridsize = cs_cuda_grid_size(n_rows, blocksize);

  double *sum_block, *res;
  cs_blas_cuda_get_2_stage_reduce_buffers(n_rows, 1, gridsize, sum_block, res);

  if (local_stream)
    cs_matrix_spmv_cuda_set_stream(stream);

  /* First iteration simplified if vx == 0
     ------------------------------------- */

  if (vx_ini != vx) {

    n_iter += 1;

    if (convergence->precision > 0. || c->plot != NULL) {
      _jacobi_compute_vx_and_residual_ini0
        <blocksize><<<gridsize, blocksize, 0, stream>>>
          (n_rows, ad_inv, rhs, vx, rk, sum_block);
      cs_blas_cuda_reduce_single_block<blocksize, 1><<<1, blocksize, 0, stream>>>
        (gridsize, sum_block, res);
      _sync_reduction_sum(c, stream, 1, res);
      residual = sqrt(*res); /* Actually, residual of previous iteration */
    }
    else {
      _jacobi_compute_vx_ini0<blocksize><<<gridsize, blocksize, 0, stream>>>
        (n_rows, ad_inv, rhs, vx, rk);
      cudaStreamSynchronize(stream);
    }

    /* Convergence test */
    if (n_iter == 1)
      c->setup_data->initial_residual = residual;
    cvg = cs_sles_it_convergence_test(c, n_iter, residual, convergence);

  }
  else
    cudaMemcpyAsync(rk, vx, n_rows*sizeof(cs_real_t),
                    cudaMemcpyDeviceToDevice, stream);

#if HAVE_GRAPH_CAPTURE > 0

  cudaGraph_t graph;
  cudaGraphExec_t graph_exec = NULL;

  /* Capture graph for a portion of kernels used here. */

  if (cs_glob_cuda_allow_graph) {
    cudaGraphNode_t graph_node;
    cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);

    /* Compute Vx <- Vx - (A-diag).Rk and residual. */

    if (convergence->precision > 0. || c->plot != NULL) {
      _jacobi_compute_vx_and_residual<blocksize><<<gridsize, blocksize, 0, stream>>>
        (n_rows, ad_inv, ad, rhs, vx, rk, sum_block);
      cs_blas_cuda_reduce_single_block<blocksize, 1><<<1, blocksize, 0, stream>>>
        (gridsize, sum_block, res);
    }
    else {
      _jacobi_compute_vx<blocksize><<<gridsize, blocksize, 0, stream>>>
        (n_rows, ad_inv, rhs, vx, rk);
    }

    cudaStreamEndCapture(stream, &graph);
    cudaError_t status = cudaGraphInstantiate(&graph_exec, graph, &graph_node,
                                              nullptr, 0);
    assert(status == cudaSuccess);
  }

#endif // HAVE_GRAPH_CAPTURE > 0

  /* Current iteration
     ----------------- */

  while (cvg == CS_SLES_ITERATING) {

    n_iter += 1;

    /* SpmV done outside of graph capture as halo synchronization
       currently synchonizes separate streams and may not be captured. */

    cs_matrix_vector_multiply_partial_d(a, CS_MATRIX_SPMV_E, rk, vx);

    /* Compute Vx <- Vx - (A-diag).Rk and residual. */

#if HAVE_GRAPH_CAPTURE > 0
    if (cs_glob_cuda_allow_graph)
      cudaGraphLaunch(graph_exec, stream);

    else
#endif // HAVE_GRAPH_CAPTURE > 0
    {
      if (convergence->precision > 0. || c->plot != NULL) {
        _jacobi_compute_vx_and_residual<blocksize><<<gridsize, blocksize, 0, stream>>>
          (n_rows, ad_inv, ad, rhs, vx, rk, sum_block);
        cs_blas_cuda_reduce_single_block<blocksize, 1><<<1, blocksize, 0, stream>>>
          (gridsize, sum_block, res);
      }
      else
        _jacobi_compute_vx<blocksize><<<gridsize, blocksize, 0, stream>>>
          (n_rows, ad_inv, rhs, vx, rk);
    }

    if (convergence->precision > 0. || c->plot != NULL) {
      _sync_reduction_sum(c, stream, 1, res);
      residual = sqrt(*res); /* Actually, residual of previous iteration */
    }
    else
      cudaStreamSynchronize(stream);

    /* Convergence test */
    if (n_iter == 1)
      c->setup_data->initial_residual = residual;
    cvg = cs_sles_it_convergence_test(c, n_iter, residual, convergence);

  }

#if HAVE_GRAPH_CAPTURE > 0
  if (cs_glob_cuda_allow_graph) {
    cudaGraphDestroy(graph);
    cudaGraphExecDestroy(graph_exec);
  }
#endif // HAVE_GRAPH_CAPTURE > 0

  if (_aux_vectors != (cs_real_t *)aux_vectors)
    cudaFree(_aux_vectors);

  if (local_stream) {
    cs_matrix_spmv_cuda_set_stream(0);
    cudaStreamDestroy(stream);
  }

  return cvg;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using block Jacobi (CUDA version).
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- diagonal block size
 *   rotation_mode   <-- halo update option for rotational periodicity
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

cs_sles_convergence_state_t
cs_sles_it_cuda_block_jacobi(cs_sles_it_t              *c,
                             const cs_matrix_t         *a,
                             cs_lnum_t                  diag_block_size,
                             cs_sles_it_convergence_t  *convergence,
                             const cs_real_t           *rhs,
                             cs_real_t                 *__restrict__ vx_ini,
                             cs_real_t                 *__restrict__ vx,
                             size_t                     aux_size,
                             void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg= CS_SLES_ITERATING;
  unsigned n_iter = 0;

  int device_id = cs_get_device_id();

  bool local_stream = false;
  cudaStream_t stream;
  stream = cs_matrix_spmv_cuda_get_stream();
  if (stream == 0) {
    local_stream = true;
    cudaStreamCreate(&stream);
  }

  const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;

  size_t vec_size = n_cols * sizeof(cs_real_t);

  /* Prefetch in case it is needed; actually, the reported allocation
     mode may be incorrect if the array are sub-arrays of a greater allocation,
     such as for multigrid, but in this case no prefetching should be needed,
     or it should have bed done by the caller. */
  {
    cs_alloc_mode_t amode_vx = cs_check_device_ptr(vx);
    cs_alloc_mode_t amode_rhs = cs_check_device_ptr(rhs);

    if (amode_vx == CS_ALLOC_HOST_DEVICE_SHARED && vx_ini == vx)
      cs_cuda_prefetch_h2d(vx, vec_size);

    if (amode_rhs == CS_ALLOC_HOST_DEVICE_SHARED)
      cs_cuda_prefetch_h2d(rhs, vec_size);
  }

  const cs_real_t  *__restrict__ ad
    =  cs_get_device_ptr_const(cs_matrix_get_diagonal(a));
  const cs_real_t *__restrict__ ad_inv
    = cs_get_device_ptr_const(c->setup_data->ad_inv);

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  double residual = -1;

  /* Allocate or map work arrays
     --------------------------- */

  assert(c->setup_data != NULL);

  cs_real_t *_aux_vectors = NULL;
  if (n_cols > 0) {
    const size_t n_wa = 1;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (   aux_vectors == NULL
        || cs_cuda_is_device_ptr(aux_vectors) == false
        || aux_size/sizeof(cs_real_t) < (wa_size * n_wa)) {
#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
      cudaMalloc(&_aux_vectors, wa_size * n_wa *sizeof(cs_real_t));
#else
      cudaMallocManaged(&_aux_vectors, wa_size * n_wa *sizeof(cs_real_t));
#endif
    }
    else
      _aux_vectors = (cs_real_t *)aux_vectors;
  }

  cs_real_t *__restrict__ rk = _aux_vectors;

  const unsigned int blocksize = 256;
  cs_lnum_t n_b_rows = n_rows / diag_block_size;

  unsigned int gridsize = cs_cuda_grid_size(n_b_rows, blocksize);

  double *sum_block, *res;
  cs_blas_cuda_get_2_stage_reduce_buffers(n_b_rows, 1, gridsize, sum_block, res);

  if (local_stream)
    cs_matrix_spmv_cuda_set_stream(stream);

  /* First iteration simplified if vx == 0
     ------------------------------------- */

  if (vx_ini != vx) {

    n_iter += 1;

    /* Compute Vx <- Vx - (A-diag).Rk and residual. */

    if (diag_block_size == 3)
      _block_3_jacobi_compute_vx_and_residual_ini0
        <blocksize><<<gridsize, blocksize, 0, stream>>>
        (n_b_rows, ad_inv, ad, rhs, vx, rk, sum_block);
    else
      _block_jacobi_compute_vx_and_residual_ini0
        <blocksize><<<gridsize, blocksize, 0, stream>>>
        (n_b_rows, diag_block_size, ad_inv, ad, rhs, vx, rk, sum_block);

    if (convergence->precision > 0. || c->plot != NULL) {
      cs_blas_cuda_reduce_single_block<blocksize, 1><<<1, blocksize, 0, stream>>>
        (gridsize, sum_block, res);
      _sync_reduction_sum(c, stream, 1, res);
      residual = sqrt(*res); /* Actually, residual of previous iteration */
    }
    else
      cudaStreamSynchronize(stream);

    /* Convergence test */
    if (n_iter == 1)
      c->setup_data->initial_residual = residual;
    cvg = cs_sles_it_convergence_test(c, n_iter, residual, convergence);

  }

  else
    cudaMemcpyAsync(rk, vx, n_rows*sizeof(cs_real_t),
                    cudaMemcpyDeviceToDevice, stream);

#if HAVE_GRAPH_CAPTURE > 0

  cudaGraph_t graph;
  cudaGraphExec_t graph_exec = NULL;

  /* Capture graph for a portion of kernels used here. */

  if (cs_glob_cuda_allow_graph) {
    cudaGraphNode_t graph_node;
    cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);

    /* Compute Vx <- Vx - (A-diag).Rk and residual. */

    if (diag_block_size == 3)
      _block_3_jacobi_compute_vx_and_residual
        <blocksize><<<gridsize, blocksize, 0, stream>>>
        (n_b_rows, ad_inv, ad, rhs, vx, rk, sum_block);
    else
      _block_jacobi_compute_vx_and_residual
        <blocksize><<<gridsize, blocksize, 0, stream>>>
        (n_b_rows, diag_block_size, ad_inv, ad, rhs, vx, rk, sum_block);

    if (convergence->precision > 0. || c->plot != NULL) {
      cs_blas_cuda_reduce_single_block<blocksize, 1><<<1, blocksize, 0, stream>>>
        (gridsize, sum_block, res);
    }

    cudaStreamEndCapture(stream, &graph);
    cudaError_t status = cudaGraphInstantiate(&graph_exec, graph, &graph_node,
                                              nullptr, 0);
    assert(status == cudaSuccess);
  }

#endif // HAVE_GRAPH_CAPTURE > 0

  /* Current iteration
     ----------------- */

  while (cvg == CS_SLES_ITERATING) {

    n_iter += 1;

    /* SpmV done outside of graph capture as halo synchronization
       currently synchonizes separate streams and may not be captured. */

    cs_matrix_vector_multiply_partial_d(a, CS_MATRIX_SPMV_E, rk, vx);

    /* Compute Vx <- Vx - (A-diag).Rk and residual. */

#if HAVE_GRAPH_CAPTURE > 0
    if (cs_glob_cuda_allow_graph)
      cudaGraphLaunch(graph_exec, stream);

    else
#endif // HAVE_GRAPH_CAPTURE > 0
    {
      if (diag_block_size == 3)
        _block_3_jacobi_compute_vx_and_residual
          <blocksize><<<gridsize, blocksize, 0, stream>>>
          (n_b_rows, ad_inv, ad, rhs, vx, rk, sum_block);
      else
        _block_jacobi_compute_vx_and_residual
          <blocksize><<<gridsize, blocksize, 0, stream>>>
          (n_b_rows, diag_block_size, ad_inv, ad, rhs, vx, rk, sum_block);

      if (convergence->precision > 0. || c->plot != NULL) {
        cs_blas_cuda_reduce_single_block<blocksize, 1><<<1, blocksize, 0, stream>>>
          (gridsize, sum_block, res);
      }
    }

    if (convergence->precision > 0. || c->plot != NULL) {
      _sync_reduction_sum(c, stream, 1, res);
      residual = sqrt(*res); /* Actually, residual of previous iteration */
    }
    else
      cudaStreamSynchronize(stream);

    /* Convergence test */
    if (n_iter == 1)
      c->setup_data->initial_residual = residual;
    cvg = cs_sles_it_convergence_test(c, n_iter, residual, convergence);

  }

#if HAVE_GRAPH_CAPTURE > 0
  if (cs_glob_cuda_allow_graph) {
    cudaGraphDestroy(graph);
    cudaGraphExecDestroy(graph_exec);
  }
#endif // HAVE_GRAPH_CAPTURE > 0

  if (_aux_vectors != (cs_real_t *)aux_vectors)
    cudaFree(_aux_vectors);

  if (local_stream) {
    cs_matrix_spmv_cuda_set_stream(0);
    cudaStreamDestroy(stream);
  }

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
 *   vx_ini          <-- initial system solution
 *                       (vx if nonzero, nullptr if zero)
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_sles_it_cuda_fcg(cs_sles_it_t              *c,
                    const cs_matrix_t         *a,
                    cs_lnum_t                  diag_block_size,
                    cs_sles_it_convergence_t  *convergence,
                    const cs_real_t           *rhs,
                    cs_real_t                 *restrict vx_ini,
                    cs_real_t                 *restrict vx,
                    size_t                     aux_size,
                    void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg = CS_SLES_ITERATING;

  int device_id = cs_get_device_id();

  bool local_stream = false;
  cudaStream_t stream;
  stream = cs_matrix_spmv_cuda_get_stream();
  if (stream == 0) {
    local_stream = true;
    cudaStreamCreate(&stream);
  }

  assert(c->setup_data != NULL);

  cs_real_t  *_aux_vectors;
  cs_real_t  *__restrict__ rk, *__restrict__ vk, *__restrict__ wk;
  cs_real_t  *__restrict__ dk, *__restrict__ qk;

  unsigned n_iter = 0;

  /* Allocate or map work arrays */
  /*-----------------------------*/

  assert(c->setup_data != NULL);

  const cs_lnum_t n_rows = c->setup_data->n_rows;
  const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;

  size_t vec_size = n_cols * sizeof(cs_real_t);

  /* Prefetch in case it is needed; actually, the reported allocation
     mode may be incorrect if the array are sub-arrays of a greater allocation,
     such as for multigrid, but in this case no prefetching should be needed,
     or it should have bed done by the caller. */
  {
    cs_alloc_mode_t amode_vx = cs_check_device_ptr(vx);
    cs_alloc_mode_t amode_rhs = cs_check_device_ptr(rhs);

    if (amode_vx == CS_ALLOC_HOST_DEVICE_SHARED && vx_ini == vx)
      cs_cuda_prefetch_h2d(vx, vec_size);

    if (amode_rhs == CS_ALLOC_HOST_DEVICE_SHARED)
      cs_cuda_prefetch_h2d(rhs, vec_size);
  }

  {
    const size_t n_wa = 5;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (   aux_vectors == nullptr
        || cs_cuda_is_device_ptr(aux_vectors) == false
        || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
       CS_MALLOC_HD(_aux_vectors, wa_size * n_wa, cs_real_t,
                    CS_ALLOC_HOST_DEVICE_SHARED);
    else
      _aux_vectors = (cs_real_t *)aux_vectors;

    rk = _aux_vectors;
    vk = _aux_vectors + wa_size;
    wk = _aux_vectors + wa_size*2;
    dk = _aux_vectors + wa_size*3;
    qk = _aux_vectors + wa_size*4;
  }

  const unsigned int blocksize = CS_BLOCKSIZE;

  cs_cuda_grid_size(n_rows, blocksize);

  unsigned int gridsize = cs_cuda_grid_size(n_rows, blocksize);
  unsigned int gridsize_blas1 = min(gridsize, 640);

  cs_blas_cuda_set_stream(stream);
  if (local_stream)
    cs_matrix_spmv_cuda_set_stream(stream);

  /* Initialize iterative calculation */
  /*----------------------------------*/

  /* Residual and descent direction */

  if (vx_ini == vx) {
    cs_matrix_vector_multiply_d(a, vx, rk);  /* rk = A.x0 */

    _fcg_init<<<gridsize, blocksize, 0, stream>>>
      (n_rows, rhs, rk, qk);
  }
  else
    _fcg_init_vx0<<<gridsize, blocksize, 0, stream>>>
      (n_rows, rhs, vx, rk, qk);

  double rho_km1 = 0;

  while (cvg == CS_SLES_ITERATING) {

    cudaStreamSynchronize(stream);

    /* Preconditioning */

    c->setup_data->pc_apply(c->setup_data->pc_context, rk, vk);

    cs_matrix_vector_multiply_d(a, vk, wk);

    /* Compute residual and prepare descent parameter */

    double alpha_k, beta_k, gamma_k, residual;

    _dot_products_vr_vw_vq_rr(c, stream, vk, rk, wk, qk,
                              &alpha_k, &beta_k, &gamma_k, &residual);

    residual = sqrt(residual);

    /* Convergence test for end of previous iteration */

    if (n_iter > 0)
      cvg = cs_sles_it_convergence_test(c, n_iter, residual, convergence);
    else
      c->setup_data->initial_residual = residual;

    if (cvg != CS_SLES_ITERATING)
      break;

    /* Complete descent parameter computation and matrix.vector product */

    if (n_iter > 0) {

      cs_real_t gk_rk1 = (abs(rho_km1) > DBL_MIN) ? gamma_k / rho_km1 : 0.;
      cs_real_t rho_k = beta_k - gamma_k * gk_rk1;
      cs_real_t ak_rk = (abs(rho_k) > DBL_MIN) ? alpha_k / rho_k : 0.;

      _fcg_update_n<<<gridsize, blocksize, 0, stream>>>
        (n_rows, gk_rk1, ak_rk, wk, dk, qk, rk, vk, vx);

      rho_km1 = rho_k;

    }
    else { /* n_iter == 0 */

      cs_real_t rho_k = beta_k;
      cs_real_t ak_rk = (abs(rho_k) > DBL_MIN) ? alpha_k / rho_k : 0.;

      _fcg_update_0<<<gridsize, blocksize, 0, stream>>>
        (n_rows, ak_rk, wk, dk, qk, rk, vk, vx);

      rho_km1 = rho_k;

    }

    n_iter += 1;

  } /* Needs iterating */

  if (_aux_vectors != aux_vectors)
    CS_FREE_HD(_aux_vectors);

  cs_blas_cuda_set_stream(0);
  if (local_stream) {
    cs_matrix_spmv_cuda_set_stream(0);
    cudaStreamDestroy(stream);
  }

  return cvg;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using optimised preconditioned GCR (CUDA version).
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- matrix
 *   diag_block_size <-- diagonal block size (unused here)
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

cs_sles_convergence_state_t
cs_sles_it_cuda_gcr(cs_sles_it_t              *c,
                    const cs_matrix_t         *a,
                    cs_lnum_t                  diag_block_size,
                    cs_sles_it_convergence_t  *convergence,
                    const cs_real_t           *rhs,
                    cs_real_t                 *vx_ini,
                    cs_real_t                 *vx,
                    size_t                     aux_size,
                    void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg= CS_SLES_ITERATING;

  int device_id = cs_get_device_id();

  bool local_stream = false;
  cudaStream_t stream;
  stream = cs_matrix_spmv_cuda_get_stream();
  if (stream == 0) {
    local_stream = true;
    cudaStreamCreate(&stream);
  }

  assert(c->setup_data != NULL);

  const cs_lnum_t n_rows = c->setup_data->n_rows;
  const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;

  size_t vec_size = n_cols * sizeof(cs_real_t);

  /* Prefetch in case it is needed; actually, the reported allocation
     mode may be incorrect if the array are sub-arrays of a greater allocation,
     such as for multigrid, but in this case no prefetching should be needed,
     or it should have bed done by the caller. */
  {
    cs_alloc_mode_t amode_vx = cs_check_device_ptr(vx);
    cs_alloc_mode_t amode_rhs = cs_check_device_ptr(rhs);

    if (amode_vx == CS_ALLOC_HOST_DEVICE_SHARED && vx_ini == vx)
      cs_cuda_prefetch_h2d(vx, vec_size);

    if (amode_rhs == CS_ALLOC_HOST_DEVICE_SHARED)
      cs_cuda_prefetch_h2d(rhs, vec_size);
  }

  double  residual = -1;

  /* In case of the standard GCR, n_k_per_restart --> Inf,
   * or stops until convergence*/

  const unsigned n_k_per_restart = c->restart_interval;
  size_t wa_size;
  unsigned n_restart = 0;

  /* Allocate or map work arrays */
  /*-----------------------------*/

  cs_real_t *_aux_vectors = NULL;
  {
    const size_t n_wa = 1 + n_k_per_restart * 2;
    wa_size = CS_SIMD_SIZE(n_cols);

    if (   aux_vectors == nullptr
        || cs_cuda_is_device_ptr(aux_vectors) == false
        || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      CS_MALLOC_HD(_aux_vectors, wa_size * n_wa, cs_real_t, CS_ALLOC_DEVICE);
    else
      _aux_vectors = (cs_real_t *)aux_vectors;

  }
  cs_real_t *__restrict__ rk = _aux_vectors;             /* store residuals  */
  cs_real_t *__restrict__ zk = _aux_vectors + wa_size;   /* store inv(M)*r   */
  cs_real_t *__restrict__ ck = _aux_vectors + wa_size * (1 + n_k_per_restart);
                                                         /* store A*zk */

  /* Use unified memory for arrays which may require an MPI reduction,
     device memory for others; also use double instead of generic
     cs_real_t (usually the same) for those arrays to avoid loss of precision
     through reduction if we ever switch cs_real_t to float. */

  /* gkj stores the upper triangle matrix Gamma of crossed inner-products
   * Not necessary to allocate the full matrix size
   * gkj_inv stores the inverse of gkj */

  cs_lnum_t n_gkj = (n_k_per_restart + 1) * n_k_per_restart / 2;
  cs_lnum_t aux_arrays_size = n_gkj + (n_k_per_restart+1) + 1 + 1;

  double *_aux_arrays;
  CS_MALLOC_HD(_aux_arrays, aux_arrays_size, double,
               CS_ALLOC_HOST_DEVICE_SHARED);

  double *__restrict__ mgkj = _aux_arrays;
  double *__restrict__ alpha = _aux_arrays + n_gkj;
  double *__restrict__ scale = alpha + n_k_per_restart+1;

  cs_real_t *gkj_inv;
  CS_MALLOC_HD(gkj_inv, n_gkj, cs_real_t, CS_ALLOC_DEVICE);

  const unsigned int blocksize = CS_BLOCKSIZE;
  const unsigned int blocksize_rsb = 512; /* cs_blas_cuda_reduce_single_block */

  cs_cuda_grid_size(n_rows, blocksize);

  unsigned int gridsize = cs_cuda_grid_size(n_rows, blocksize);
  unsigned int gridsize_blas1 = min(gridsize, 640);

  double *sum_block, *res;
  cs_blas_cuda_get_2_stage_reduce_buffers(n_rows, 1, gridsize, sum_block, res);

  cs_blas_cuda_set_stream(stream);
  if (local_stream)
    cs_matrix_spmv_cuda_set_stream(stream);

  /* Current Restart */

  cs_real_t *_vx = vx_ini;

  while (cvg == CS_SLES_ITERATING) {

    unsigned n_c_iter = 0;

    /* Initialize iterative calculation */
    /*----------------------------------*/

    if (_vx == vx) {
      cs_matrix_vector_multiply_d(a, vx, rk);  /* rk = A.x0 */

      _ymx_dot_yy<blocksize_rsb><<<gridsize_blas1, blocksize_rsb, 0, stream>>>
        (n_rows, rhs, rk, sum_block);
    }
    else {
      _x0_ymr_dot_yy<blocksize_rsb><<<gridsize_blas1, blocksize_rsb, 0, stream>>>
        (n_rows, rhs, vx, rk, sum_block);
    }

    cs_blas_cuda_reduce_single_block
      <blocksize_rsb, 1><<<1, blocksize_rsb, 0, stream>>>
      (gridsize_blas1, sum_block, res);

    _sync_reduction_sum(c, stream, 1, res);

    cudaMemcpy(&residual, res, sizeof(double), cudaMemcpyDeviceToHost);
    residual = sqrt(residual);

    if (n_restart == 0)
      c->setup_data->initial_residual = residual;

    /* Current Iteration on k */
    /* ---------------------- */

    while (cvg == CS_SLES_ITERATING && n_c_iter < n_k_per_restart) {

      /* Preconditionning */

      cs_real_t *zk_n = zk + n_c_iter * wa_size;
      cs_real_t *ck_n = ck + n_c_iter * wa_size;

      /* Preconditionning */

      c->setup_data->pc_apply(c->setup_data->pc_context, rk, zk_n);

      cs_matrix_vector_multiply_d(a, zk_n, ck_n);

      /* Compute the ck_n direction;

         Remarks:

         - Compared to the host-only code_saturne implementation, we change the
           sign of computed Gamma (gkj) versions to allow for axpy. To avoid
           extra host/device synchronizations, we apply the sign change to all
           subsequent operations involving those values.

         - Compared to the PETSc implementation, we use one less array
           and axpy operation, but we cannot group dot products, leading to lower
           computational cost but higher expected MPI overhead.
           We also need separate saxpy operations, while grouping them would
           almost double compute intensity (a same vector is used in n_c_iter
           operations), possibly compensating for the additional axpy series.
      */

      for (cs_lnum_t jj = 0; jj < (int)n_c_iter; jj++) {
        cs_real_t *ck_j = ck + jj * wa_size;

        int ii_jn = (n_c_iter + 1) * n_c_iter / 2 + jj;
        mgkj[ii_jn] = - _dot_product(c, stream, ck_j, ck_n);

        /* changed sign directly above to allow use of axpy */

        cs_blas_cuda_axpy(n_rows, &mgkj[ii_jn], ck_j, ck_n);
      }

      const int iter_shift = (n_c_iter+1) * n_c_iter / 2 + n_c_iter;
      mgkj[iter_shift] = - sqrt(_dot_product(c, stream, ck_n, ck_n));

      if (abs(mgkj[iter_shift]) > 0) {
        *scale = 1. / (- mgkj[iter_shift]);

        _y_scale_dot_xy
          <blocksize_rsb><<<gridsize_blas1, blocksize_rsb, 0, stream>>>
          (n_rows, scale, rk, ck_n, sum_block);
        cs_blas_cuda_reduce_single_block
          <blocksize_rsb, 1><<<1, blocksize_rsb, 0, stream>>>
          (gridsize, sum_block, &alpha[n_c_iter]);
        _sync_reduction_sum(c, stream, 1, &alpha[n_c_iter]);

      }
      else
        alpha[n_c_iter] = 0.;

      _smaxpy_dot_yy<blocksize_rsb><<<gridsize_blas1, blocksize_rsb, 0, stream>>>
        (n_rows, &alpha[n_c_iter], ck_n, rk, sum_block);
      cs_blas_cuda_reduce_single_block
        <blocksize_rsb, 1><<<1, blocksize_rsb, 0, stream>>>
        (gridsize, sum_block, res);
      _sync_reduction_sum(c, stream, 1, res);

      cudaMemcpy(&residual, res, sizeof(double), cudaMemcpyDeviceToHost);
      residual = sqrt(residual);

      n_c_iter += 1;

      /* Convergence test of current iteration */
      cvg = cs_sles_it_convergence_test(c,
                                        (n_restart*n_k_per_restart) + n_c_iter,
                                        residual,
                                        convergence);

      if (cvg != CS_SLES_ITERATING)
        break;

    } /* Needs iterating or k < n_restart */

    /* Inversion of Gamma;

     * This is a basically small, serial operation, done only on one thread.
     * we should check whether it is more efficient to run this on the
     * GPU (with slower code but avoiding memory transfer latency)
     * or on host CPU. */

    if (n_c_iter == 1)
      _gcr_gkj_inv_1<<<1, CS_CUDA_WARP_SIZE, 0, stream>>>
        (mgkj, gkj_inv);
    else
      _gcr_gkj_inv<<<1, CS_CUDA_WARP_SIZE, 0, stream>>>
        ((int)n_c_iter, n_gkj, mgkj, gkj_inv);

    /* Compute the solution */
    _gcr_update_vx<<<gridsize, blocksize, 0, stream>>>
      (n_rows, n_c_iter, wa_size, alpha, gkj_inv, zk, vx);

    n_restart += 1;

  } /* Needs iterating */

  if (_aux_vectors != aux_vectors)
    CS_FREE_HD(_aux_vectors);

  CS_FREE_HD(_aux_arrays);
  CS_FREE_HD(gkj_inv);

  cs_blas_cuda_set_stream(0);
  if (local_stream) {
    cs_matrix_spmv_cuda_set_stream(0);
    cudaStreamDestroy(stream);
  }

  return cvg;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
