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

static cudaStream_t _stream = 0;

static double  *_r_reduce = NULL;
static double  *_r_grid = NULL;
static unsigned int  _r_grid_size = 0;

#if defined(HAVE_CUBLAS)

static bool            _prefer_cublas = false;
static cublasHandle_t  _handle = NULL;

#endif

/*============================================================================
 * Private kernel definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the absolute sum of vector values,
 *        summing result over all threads of a block.
 *
 * blockSize must be a power of 2.
 *
 * \param[in]   n      array size
 * \param[in]   x      x vector
 * \param[out]  b_res  result of s = x.x
 */
/*----------------------------------------------------------------------------*/

template <size_t blockSize, typename T>
__global__ static void
_asum_stage_1_of_2(cs_lnum_t    n,
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

  // Output: b_res for this block

  cs_blas_cuda_block_reduce_sum<blockSize>(stmp, tid, b_res);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute dot product x.y, summing result over all threads of a block.
 *
 * blockSize must be a power of 2.
 *
 * \param[in]   n      array size
 * \param[in]   x      x vector
 * \param[in]   y      y vector
 * \param[out]  b_res  result of s = x.x
 */
/*----------------------------------------------------------------------------*/

template <size_t blockSize, typename T>
__global__ static void
_dot_xy_stage_1_of_2(cs_lnum_t    n,
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

  // Output: b_res for this block

  cs_blas_cuda_block_reduce_sum<blockSize>(stmp, tid, b_res);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute dot product x.x, summing result over all threads of a block.
 *
 * blockSize must be a power of 2.
 *
 * \param[in]   n      array size
 * \param[in]   x      x vector
 * \param[in]   y      y vector
 * \param[out]  b_res  result of s = x.x
 */
/*----------------------------------------------------------------------------*/

template <size_t blockSize, typename T>
__global__ static void
_dot_xx_stage_1_of_2(cs_lnum_t    n,
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
    stmp[tid] += static_cast<double>(x[i] * x[i]);
  }

  // Output: b_res for this block

  cs_blas_cuda_block_reduce_sum<blockSize>(stmp, tid, b_res);
}

/*----------------------------------------------------------------------------
 * Compute y <- alpha.x + y (device kernel)
 *
 * parameters:
 *   n     <-- number of elements
 *   alpha <-- constant value
 *   x     <-- vector of elements
 *   y     <-> vector of elements
 *----------------------------------------------------------------------------*/

__global__ static void
_axpy(cs_lnum_t         n,
      const cs_real_t  *alpha,
      const cs_real_t  *restrict x,
      cs_real_t        *restrict y)
 {
   cs_real_t _alpha = *alpha;
   cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;

   size_t grid_size = blockDim.x*gridDim.x;
   while (ii < n){
     y[ii] += _alpha * x[ii];
     ii += grid_size;
   }
 }

/*----------------------------------------------------------------------------
 * Compute x <- alpha.x (device kernel)
 *
 * parameters:
 *   n     <-- number of elements
 *   alpha <-- constant value
 *   x     <-> vector of elements
 *----------------------------------------------------------------------------*/

__global__ static void
_scal(cs_lnum_t         n,
      const cs_real_t  *alpha,
      cs_real_t        *restrict x)
 {
   cs_real_t _alpha = *alpha;
   cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;

   size_t grid_size = blockDim.x*gridDim.x;
   while (ii < n){
     x[ii] *= _alpha;
     ii += grid_size;
   }
 }

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize cuBLAS and return a handle.
 *
 * \return  pointer to global cuBLAS handle
 */
/*----------------------------------------------------------------------------*/

static cublasHandle_t
_init_cublas(void)
{
  cublasStatus_t status = cublasCreate(&_handle);

  if (status != CUBLAS_STATUS_SUCCESS)
#if CUBLAS_VERSION >= 11600
    bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
              __func__, cublasGetStatusString(status));
#else
  bft_error(__FILE__, __LINE__, 0, _("%s: cuBLAS error %d."),
            __func__, (int)status);
#endif

  return _handle;
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
  unsigned int grid_size = cs_cuda_grid_size(n, block_size);

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
cs_blas_cuda_set_stream(cudaStream_t  stream)
{
  _stream = stream;

#if defined(HAVE_CUBLAS)

  if (_handle != NULL) {
    cublasSetStream(_handle, stream);
  }

#endif
}

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
                                       unsigned int  grid_size)
{
  unsigned int t_grid_size = grid_size * tuple_size;

  if (_r_grid_size < t_grid_size) {
    CS_REALLOC_HD(_r_grid, t_grid_size, double, CS_ALLOC_DEVICE);
    _r_grid_size = t_grid_size;
  }

  return _r_grid;
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

  _asum_stage_1_of_2<block_size><<<grid_size, block_size, 0, _stream>>>
    (n, x, _r_grid);
  cs_blas_cuda_reduce_single_block<block_size><<<1, block_size, 0, _stream>>>
    (grid_size, _r_grid, _r_reduce);

  /* Need to synchronize stream in all cases so as to
     have up-to-date value in returned _r_reduce[0] */

  cudaStreamSynchronize(_stream);

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

  if (x != y)
    _dot_xy_stage_1_of_2<block_size><<<grid_size, block_size, 0, _stream>>>
      (n, x, y, _r_grid);
  else
    _dot_xx_stage_1_of_2<block_size><<<grid_size, block_size, 0, _stream>>>
      (n, x, _r_grid);
  cs_blas_cuda_reduce_single_block<block_size><<<1, block_size, 0, _stream>>>
    (grid_size, _r_grid, _r_reduce);

  /* Need to synchronize stream in all cases so as to
     have up-to-date value in returned _r_reduce[0] */

  cudaStreamSynchronize(_stream);

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

  if (_handle == NULL)
    _handle = _init_cublas();

  if (sizeof(cs_real_t) == 8)
    status = cublasDasum(_handle, n, (double *)x, 1, &dot);
  else if (sizeof(cs_real_t) == 4) {
    float fdot;
    status = cublasSasum(_handle, n, (float *)x, 1, &fdot);
    dot = fdot;
  }

  if (status != CUBLAS_STATUS_SUCCESS) {
#if CUBLAS_VERSION >= 11600
    bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
              __func__, cublasGetStatusString(status));
#else
    bft_error(__FILE__, __LINE__, 0, _("%s: cuBLAS error %d."),
              __func__, (int)status);
#endif
  }

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

  if (_handle == NULL)
    _handle = _init_cublas();

  if (sizeof(cs_real_t) == 8)
    cublasDdot(_handle, n, (double *)x, 1, (double *)y, 1, &dot);
  else if (sizeof(cs_real_t) == 4) {
    float fdot;
    cublasSdot(_handle, n, (float *)x, 1, (float *)y, 1, &fdot);
    dot = fdot;
  }

  if (status != CUBLAS_STATUS_SUCCESS) {
#if CUBLAS_VERSION >= 11600
    bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
              __func__, cublasGetStatusString(status));
#else
    bft_error(__FILE__, __LINE__, 0, _("%s: cuBLAS error %d."),
              __func__, (int)status);
#endif
  }

  return dot;
}

#endif // defined(HAVE_CUBLAS)

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
                  cs_real_t        *restrict y)
{
#if defined(HAVE_CUBLAS)

  if (_prefer_cublas && _handle != NULL) {
    _handle = _init_cublas();

    cublasStatus_t status = CUBLAS_STATUS_SUCCESS;

    if (sizeof(cs_real_t) == 8)
      status = cublasDaxpy(_handle, n,
                           (double *)alpha,
                           (double *)x, 1,
                           (double *)y, 1);
    else if (sizeof(cs_real_t) == 4)
      status = cublasSaxpy(_handle, n,
                           (float *)alpha,
                           (float *)x, 1,
                           (float *)y, 1);

    if (status != CUBLAS_STATUS_SUCCESS) {
#if CUBLAS_VERSION >= 11600
      bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
                __func__, cublasGetStatusString(status));
#else
      bft_error(__FILE__, __LINE__, 0, _("%s: cuBLAS error %d."),
                __func__, (int)status);
#endif
    }

    return;
  }

#endif /* defined(HAVE_CUBLAS) */

  /* If not using cuBLAS for this operation */

  const unsigned int blocksize = 256;  // try 640 also

  unsigned int gridsize = cs_cuda_grid_size(n, blocksize);
  // gridsize = min(gridsize, 640;

  _axpy<<<gridsize, blocksize, 0, _stream>>>(n, alpha, x, y);
}

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
                  cs_real_t        *restrict x)
{
#if defined(HAVE_CUBLAS)

  if (_prefer_cublas && _handle != NULL) {
    _handle = _init_cublas();

    cublasStatus_t status = CUBLAS_STATUS_SUCCESS;

    if (sizeof(cs_real_t) == 8)
      status = cublasDscal(_handle, n,
                           (double *)alpha,
                           (double *)x, 1);
    else if (sizeof(cs_real_t) == 4)
      status = cublasSscal(_handle, n,
                           (float *)alpha,
                           (float *)x, 1);

    if (status != CUBLAS_STATUS_SUCCESS) {
#if CUBLAS_VERSION >= 11600
      bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
                __func__, cublasGetStatusString(status));
#else
      bft_error(__FILE__, __LINE__, 0, _("%s: cuBLAS error %d."),
                __func__, (int)status);
#endif
    }

    return;
  }

#endif /* defined(HAVE_CUBLAS) */

  /* If not using cuBLAS for this operation */

  const unsigned int blocksize = 256;  // try 640 also

  unsigned int gridsize = cs_cuda_grid_size(n, blocksize);
  // gridsize = min(gridsize, 640;

  _scal<<<gridsize, blocksize, 0, _stream>>>(n, alpha, x);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
