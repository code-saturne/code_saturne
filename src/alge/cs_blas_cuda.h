#ifndef __CS_BLAS_CUDA_H__
#define __CS_BLAS_CUDA_H__

/*============================================================================
 * BLAS (Basic Linear Algebra Subroutine) functions
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
 * External library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_base_cuda.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

END_C_DECLS

/*============================================================================
 * Templated function definitions
 *============================================================================*/

#if defined(__CUDACC__)

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Kernel for sum reduction within a warp (for warp size 32).
 *
 * \param[in, out]  stmp  shared value to reduce
 * \param[in, out]  tid   thread id
 */
/*----------------------------------------------------------------------------*/

template <size_t blockSize, typename T>
__device__ static void __forceinline__
cs_blas_cuda_warp_reduce_sum(volatile T  *stmp,
                             size_t       tid)
{
  if (blockSize >= 64) stmp[tid] += stmp[tid + 32];
  if (blockSize >= 32) stmp[tid] += stmp[tid + 16];
  if (blockSize >= 16) stmp[tid] += stmp[tid +  8];
  if (blockSize >=  8) stmp[tid] += stmp[tid +  4];
  if (blockSize >=  4) stmp[tid] += stmp[tid +  2];
  if (blockSize >=  2) stmp[tid] += stmp[tid +  1];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute dot product x.y, summing result over all threads of a block.
 *
 * blockSize must be a power of 2.
 *
 * This kernel uses explicit loop unrolling, which seems to offer equal or
 * slightly better performance than the non-unrolled version on an Ampere
 * architecure with CUDA 11.
 *
 * \param[in, out]  stmp       shared value to reduce
 * \param[in, out]  tid        thread id
 * \param[out]      sum_block  contribution of this block to the sum
 */
/*----------------------------------------------------------------------------*/

template <size_t blockSize, typename T>
__device__ static void __forceinline__
cs_blas_cuda_block_reduce_sum(T       *stmp,
                              size_t   tid,
                              double  *sum_block)
{
  __syncthreads();

  /* Loop explicitely unrolled */

  if (blockSize >= 1024) {
    if (tid < 512) {
      stmp[tid] += stmp[tid + 512];
    }
    __syncthreads();
  }
  if (blockSize >=  512) {
    if (tid < 256) {
      stmp[tid] += stmp[tid + 256];
    }
    __syncthreads();
  }
  if (blockSize >=  256) {
    if (tid < 128) {
      stmp[tid] += stmp[tid + 128];
    } __syncthreads();
  }
  if (blockSize >=  128) {
    if (tid <  64) {
      stmp[tid] += stmp[tid +  64];
    } __syncthreads();
  }

  if (tid < 32)
    cs_blas_cuda_warp_reduce_sum<blockSize>(stmp, tid);

  // Output: b_res for this block

  if (tid == 0) sum_block[blockIdx.x] = stmp[0];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Kernel for sum reduction within a block.
 *
 * \param[in, out]  n        number of values to reduce
 * \param[in, out]  g_idata  input array (size n)
 * \param[in, out]  g_odata  onput array (size 1)
 */
/*----------------------------------------------------------------------------*/

template <size_t blockSize, typename T>
__global__ static void
cs_blas_cuda_reduce_single_block(size_t   n,
                                 T       *g_idata,
                                 T       *g_odata)
{
  __shared__ T sdata[blockSize];

  size_t tid = threadIdx.x;
  T r_s = 0;

  for (int i = threadIdx.x; i < n; i+= blockSize)
    r_s += g_idata[i];

  sdata[tid] = r_s;
  __syncthreads();

  for (int j = blockSize/2; j > CS_CUDA_WARP_SIZE; j /= 2) {
    if (tid < j) {
      sdata[tid] += sdata[tid + j];
    }
    __syncthreads();
  }

  if (tid < 32) cs_blas_cuda_warp_reduce_sum<blockSize>(sdata, tid);
  if (tid == 0) *g_odata = sdata[0];
}

#endif /* defined(__CUDACC__) */

BEGIN_C_DECLS

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize CUDA BLAS API.
 *
 * This frees resources such as the cuBLAS handle, if used.
 */
/*----------------------------------------------------------------------------*/

void
cs_blas_cuda_finalize(void);

#if defined(__CUDACC__)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign CUDA stream for next CUDA-based blas operations.
 *
 * If a stream other than the default stream (0) is used, it will not be
 * synchronized automatically after sparse matrix-vector products (so as to
 * avoid the corresponding overhead), so the caller will need to manage
 * stream syncronization manually.
 *
 * This function is callable only from CUDA code.
 */
/*----------------------------------------------------------------------------*/

void
cs_blas_cuda_set_stream(cudaStream_t  stream);

#endif /* defined(__CUDACC__) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to reduction buffer needed for 2-stage reductions.
 *
 * This buffer is used internally by all cs_blas_cuda 2-stage operations,
 * allocated and resized updon demand, and freed when calling
 * cs_blas_cuda_finalize, so it is assumed no two operations (in different
 * threads) use this simultaneously.
 *
 * Also check initialization of work arrays.
 *
 * \param[in]  n           size of arrays
 * \param[in]  tuple_size  number of values per tuple simultaneously reduced
 * \param[in]  grid_size   associated grid size
 *
 * \return  pointer to reduction bufffer.
 */
/*----------------------------------------------------------------------------*/

double *
cs_blas_cuda_get_2_stage_reduce_buffer(cs_lnum_t     n,
                                       cs_lnum_t     tuple_size,
                                       unsigned int  grid_size);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the absolute sum of vector values using CUDA.
 *
 * \param[in]  n  size of array x
 * \param[in]  x  array of floating-point values (on device)
 *
 * \return  sum of absolute array values
 */
/*----------------------------------------------------------------------------*/

double
cs_blas_cuda_asum(cs_lnum_t        n,
                 const cs_real_t  x[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the dot product of 2 vectors: x.y using CUDA.
 *
 * \param[in]  n  size of arrays x and y
 * \param[in]  x  array of floating-point values (on device)
 * \param[in]  y  array of floating-point values (on device)
 *
 * \return  dot product
 */
/*----------------------------------------------------------------------------*/

double
cs_blas_cuda_dot(cs_lnum_t        n,
                 const cs_real_t  x[],
                 const cs_real_t  y[]);

#if defined(HAVE_CUBLAS)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the absolute sum of vector values using cuBLAS.
 *
 * \param[in]  n  size of arrays x and y
 * \param[in]  x  array of floating-point values (on device)
 *
 * \return  sum of absolute array values
 */
/*----------------------------------------------------------------------------*/

double
cs_blas_cublas_asum(cs_lnum_t        n,
                    const cs_real_t  x[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the dot product of 2 vectors: x.y using cuBLAS.
 *
 * \param[in]  n  size of arrays x and y
 * \param[in]  x  array of floating-point values (on device)
 * \param[in]  y  array of floating-point values (on device)
 *
 * \return  dot product
 */
/*----------------------------------------------------------------------------*/

double
cs_blas_cublas_dot(cs_lnum_t        n,
                   const cs_real_t  x[],
                   const cs_real_t  y[]);

#endif  /* defined(HAVE_CUBLAS) */

/*----------------------------------------------------------------------------
 * Compute y <- alpha.x + y
 *
 * This function may be set to use either cuBLAS or a local kernel.
 *
 * parameters:
 *   n      <-- number of elements
 *   alpha  <-- constant value (on device)
 *   x      <-- vector of elements (on device)
 *   y      <-> vector of elements (on device)
 *----------------------------------------------------------------------------*/

void
cs_blas_cuda_axpy(cs_lnum_t         n,
                  const cs_real_t  *alpha,
                  const cs_real_t  *restrict x,
                  cs_real_t        *restrict y);

/*----------------------------------------------------------------------------
 * Compute x <- alpha.x
 *
 * This function may be set to use either cuBLAS or a local kernel.
 *
 * parameters:
 *   n      <-- number of elements
 *   alpha  <-- constant value (on device)
 *   x      <-> vector of elements (on device)
 *----------------------------------------------------------------------------*/

void
cs_blas_cuda_scal(cs_lnum_t         n,
                  const cs_real_t  *alpha,
                  cs_real_t        *restrict x);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BLAS_CUDA_H__ */
