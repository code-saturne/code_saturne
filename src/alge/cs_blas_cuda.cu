/*============================================================================
 * BLAS (Basic Linear Algebra Subroutine) functions using CUDA.
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
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <algorithm>
#include <assert.h>
#include <cuda.h>

#if defined(HAVE_CUBLAS)
#include <cublas_v2.h>
#endif

#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#if (CUDART_VERSION >= 11000)
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
namespace cg = cooperative_groups;
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_base_accel.h"
#include "cs_base_cuda.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_blas_cuda.h"

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define BLOCKSIZE 256

/*============================================================================
 *  Global variables
 *============================================================================*/

static double  *_r_reduce = NULL;
static double  *_r_grid = NULL;
static unsigned int  _r_grid_size = 0;

#if defined(HAVE_CUBLAS)

static cublasHandle_t  _handle = NULL;

#endif

/*============================================================================
 * Private kernel definitions
 *============================================================================*/

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
 * \param[in]   n      array size
 * \param[in]   x      x vector
 * \param[out]  b_res  result of s = x.x
 */
/*----------------------------------------------------------------------------*/

template <size_t blockSize, typename T>
__global__ static void
_asum_stage_1_of_2_u(cs_lnum_t    n,
                     const T     *x,
                     double      *b_res)
{
  __shared__ double stmp[blockSize];

  cs_lnum_t tid = threadIdx.x;
  size_t grid_size = blockDim.x*gridDim.x;

  stmp[tid] = 0;
  for (int i = blockIdx.x*(blockDim.x) + tid;
       i < n;
       i += grid_size) {
    stmp[tid] += fabs(static_cast<double>(x[i]));
  }

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

  if (tid == 0) b_res[blockIdx.x] = stmp[0];
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
 * \param[in]   n      array size
 * \param[in]   x      x vector
 * \param[in]   y      y vector
 * \param[out]  b_res  result of s = x.x
 */
/*----------------------------------------------------------------------------*/

template <size_t blockSize, typename T>
__global__ static void
_dot_xy_stage_1_of_2_u(cs_lnum_t    n,
                       const T     *x,
                       const T     *y,
                       double      *b_res)
{
  __shared__ double stmp[blockSize];

  cs_lnum_t tid = threadIdx.x;
  size_t grid_size = blockDim.x*gridDim.x;

  stmp[tid] = 0;
  for (int i = blockIdx.x*(blockDim.x) + tid;
       i < n;
       i += grid_size) {
    stmp[tid] += static_cast<double>(x[i] * y[i]);
  }

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

  if (tid == 0) b_res[blockIdx.x] = stmp[0];
}

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute grid size for given array and block sizes.
 *
 * Also check initialization of work arrays.
 *
 * \param[in]  n           size of arrays
 * \param[in]  block_size  block size for kernels
 *
 * \return  grid size for kernels
 */
/*----------------------------------------------------------------------------*/

unsigned int
_grid_size(cs_lnum_t     n,
           unsigned int  block_size)
{
  unsigned int grid_size = (n % block_size) ?  n/block_size : n/block_size + 1;

  if (_r_reduce == NULL) {
    // large enough for 3-way combined dot product.
    CS_REALLOC_HD(_r_reduce, 3, double, CS_ALLOC_HOST_DEVICE_SHARED);
  }

  if (_r_grid_size < grid_size) {
    CS_REALLOC_HD(_r_grid, grid_size, double, CS_ALLOC_DEVICE);
    _r_grid_size = grid_size;
  }

  return grid_size;
}

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize CUDA BLAS API.
 *
 * This frees resources such as the cuBLAS handle, if used.
 */
/*----------------------------------------------------------------------------*/

void
cs_blas_cuda_finalize(void)
{
  CS_FREE_HD(_r_grid);
  _r_grid_size = 0;

#if defined(HAVE_CUBLAS)

  if (_handle != NULL) {
    cublasDestroy(_handle);
    _handle = NULL;
  }

#endif
}

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
                 const cs_real_t  x[])
{
  const unsigned int block_size = 256;
  unsigned int grid_size = _grid_size(n, block_size);

  _asum_stage_1_of_2_u<block_size><<<grid_size, block_size, 0>>>
    (n, x, _r_grid);
  cs_blas_cuda_reduce_single_block<block_size><<<1, block_size, 0>>>
    (grid_size, _r_grid, _r_reduce);

  cudaStreamSynchronize(0);

  return _r_reduce[0];
}

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
                 const cs_real_t  y[])
{
  const unsigned int block_size = 256;
  unsigned int grid_size = _grid_size(n, block_size);

  _dot_xy_stage_1_of_2_u<block_size><<<grid_size, block_size, 0>>>
    (n, x, y, _r_grid);
  cs_blas_cuda_reduce_single_block<block_size><<<1, block_size, 0>>>
    (grid_size, _r_grid, _r_reduce);

  cudaStreamSynchronize(0);

  return _r_reduce[0];
}

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
                    const cs_real_t  x[])
{
  double dot = 0.;

  cublasStatus_t status = CUBLAS_STATUS_SUCCESS;

  if (_handle == NULL) {
    status = cublasCreate(&_handle);
    if (status != CUBLAS_STATUS_SUCCESS)
#if CUBLAS_VERSION >= 11600
      bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
                __func__, cublasGetStatusString(status));
#else
      bft_error(__FILE__, __LINE__, 0, _("%s: cuBLAS error %d."),
                __func__, (int)status);
#endif
  }

  if (sizeof(cs_real_t) == 8)
    cublasDasum(_handle, n, (double *)x, 1, &dot);
  else if (sizeof(cs_real_t) == 4) {
    float fdot;
    cublasSasum(_handle, n, (float *)x, 1, &fdot);
    dot = fdot;
  }

  if (status != CUBLAS_STATUS_SUCCESS)
#if CUBLAS_VERSION >= 11600
    bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
              __func__, cublasGetStatusString(status));
#else
    bft_error(__FILE__, __LINE__, 0, _("%s: cuBLAS error %d."),
              __func__, (int)status);
#endif

  return dot;
}

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
                   const cs_real_t  y[])
{
  double dot = 0.;

  cublasStatus_t status = CUBLAS_STATUS_SUCCESS;

  if (_handle == NULL) {
    status = cublasCreate(&_handle);
    if (status != CUBLAS_STATUS_SUCCESS)
#if CUBLAS_VERSION >= 11600
      bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
                __func__, cublasGetStatusString(status));
#else
      bft_error(__FILE__, __LINE__, 0, _("%s: cuBLAS error %d."),
                __func__, (int)status);
#endif
  }

  if (sizeof(cs_real_t) == 8)
    cublasDdot(_handle, n, (double *)x, 1, (double *)y, 1, &dot);
  else if (sizeof(cs_real_t) == 4) {
    float fdot;
    cublasSdot(_handle, n, (float *)x, 1, (float *)y, 1, &fdot);
    dot = fdot;
  }

  if (status != CUBLAS_STATUS_SUCCESS)
#if CUBLAS_VERSION >= 11600
    bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
              __func__, cublasGetStatusString(status));
#else
    bft_error(__FILE__, __LINE__, 0, _("%s: cuBLAS error %d."),
              __func__, (int)status);
#endif

  return dot;
}

#endif // defined(HAVE_CUBLAS)

/*----------------------------------------------------------------------------*/

END_C_DECLS
