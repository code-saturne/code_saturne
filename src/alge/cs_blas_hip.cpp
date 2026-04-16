/*============================================================================
 * BLAS (Basic Linear Algebra Subroutine) functions using HIP.
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
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <algorithm>
#include <assert.h>

#if defined(HAVE_HIPBLAS)
#include <hipblas.h>
#endif

#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"

#include "base/cs_base.h"
#include "base/cs_mem.h"
#include "base/cs_base_accel.h"
#include "base/cs_base_hip.h"
#include "base/cs_dispatch.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "alge/cs_blas_hip.h"

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define BLOCKSIZE 256

/*============================================================================
 *  Global variables
 *============================================================================*/

static double  *_r_reduce = nullptr;
static double  *_r_grid = nullptr;
static unsigned int  _r_grid_size = 0;
static unsigned int  _r_tuple_size = 0;

#if defined(HAVE_HIPBLAS)

static bool             _prefer_hipblas = false;
static hipblasHandle_t  _handle = nullptr;

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

  cs_blas_hip_block_reduce_sum<blockSize, 1>(stmp, tid, b_res);
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

  cs_blas_hip_block_reduce_sum<blockSize, 1>(stmp, tid, b_res);
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

  cs_blas_hip_block_reduce_sum<blockSize, 1>(stmp, tid, b_res);
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

/*============================================================================
 * Private function definitions
 *============================================================================*/

#if defined(HAVE_HIPBLAS)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check hipBLAS status.
 *
 * \param[in]  func_name  name of caller function
 * \param[in]  file_name  name of source file from which error handler called
 * \param[in]  line_num   line of source file from which error handler called
 * \param[in]  status     hipBLAS status code
 */
/*----------------------------------------------------------------------------*/

static void
_check_hipblas_status(const char       *func_name,
                      const char       *file_name,
                      int               line_num,
                      hipblasStatus_t   status)
{
  if (status != HIPBLAS_STATUS_SUCCESS)
    bft_error(file_name, line_num, 0, _("%s: hipBLAS error %d."),
              func_name, (int)status);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize hipBLAS and return a handle.
 *
 * \return  pointer to global hipBLAS handle
 */
/*----------------------------------------------------------------------------*/

static hipblasHandle_t
_init_hipblas(void)
{
  hipblasStatus_t status = hipblasCreate(&_handle);
  _check_hipblas_status(__func__, __FILE__, __LINE__, status);

  return _handle;
}

/*----------------------------------------------------------------------------
 * Compute x <- alpha.x
 *
 * This function may be set to use either hipBLAS or a local kernel.
 *
 * parameters:
 *   stream <-- associated CUDA stream
 *   n      <-- number of elements
 *   alpha  <-- constant value (on device)
 *   x      <-> vector of elements (on device)
 *----------------------------------------------------------------------------*/

static void
_set_hipblas_stream(hipStream_t  stream)
{
  if (_handle == nullptr)
    _handle = _init_hipblas();

  hipStream_t  hipblas_stream = 0;
  hipblasStatus_t status = hipblasGetStream(_handle, &hipblas_stream);
  _check_hipblas_status(__func__, __FILE__, __LINE__ -1, status);

  if (status == HIPBLAS_STATUS_SUCCESS && hipblas_stream != stream) {
    status = hipblasSetStream(_handle, stream);
    _check_hipblas_status(__func__, __FILE__, __LINE__ -1, status);
  }
}

#endif /* defined(HAVE_HIPBLAS) */

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
  unsigned int grid_size = cs_hip_grid_size(n, block_size);

  if (_r_reduce == nullptr) {
    // large enough for 4-way combined dot product.
    _r_tuple_size = 4;
    CS_REALLOC_HD(_r_reduce, _r_tuple_size, double, CS_ALLOC_HOST_DEVICE_SHARED);
  }

  if (_r_grid_size < grid_size) {
    CS_REALLOC_HD(_r_grid, grid_size, double, CS_ALLOC_DEVICE);
    _r_grid_size = grid_size;
  }

  return grid_size;
}

/*----------------------------------------------------------------------------*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize HIP BLAS API.
 *
 * This frees resources such as the hipBLAS handle, if used.
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_blas_hip_finalize(void)
{
  CS_FREE(_r_reduce);
  CS_FREE(_r_grid);
  _r_grid_size = 0;

#if defined(HAVE_HIPBLAS)

  if (_handle != nullptr) {
    hipblasDestroy(_handle);
    _handle = nullptr;
  }

#endif
}

#if defined(__HIPCC__)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the absolute sum of vector values using HIP.
 *
 * \param[in]  stream  associated HIP stream
 * \param[in]  n       size of array x
 * \param[in]  x       array of floating-point values (on device)
 *
 * \return  sum of absolute array values
 */
/*----------------------------------------------------------------------------*/

double
cs_blas_hip_asum(hipStream_t      stream,
                 cs_lnum_t        n,
                 const cs_real_t  x[])
{
  const unsigned int block_size = 256;
  unsigned int grid_size = _grid_size(n, block_size);

  _asum_stage_1_of_2<block_size><<<grid_size, block_size, 0, stream>>>
    (n, x, _r_grid);
  cs_blas_hip_reduce_single_block<block_size, 1><<<1, block_size, 0, stream>>>
    (grid_size, _r_grid, _r_reduce);

  /* Need to synchronize stream in all cases so as to
     have up-to-date value in returned _r_reduce[0] */

  CS_HIP_CHECK(hipStreamSynchronize(stream));
  CS_HIP_CHECK(hipGetLastError());

  return _r_reduce[0];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the dot product of 2 vectors: x.y using HIP.
 *
 * \param[in]  stream  associated HIP stream
 * \param[in]  n       size of arrays x and y
 * \param[in]  x       array of floating-point values (on device)
 * \param[in]  y       array of floating-point values (on device)
 *
 * \return  dot product
 */
/*----------------------------------------------------------------------------*/

double
cs_blas_hip_dot(hipStream_t      stream,
                cs_lnum_t        n,
                const cs_real_t  x[],
                const cs_real_t  y[])
{
  const unsigned int block_size = 256;
  unsigned int grid_size = _grid_size(n, block_size);

  if (x != y)
    _dot_xy_stage_1_of_2<block_size><<<grid_size, block_size, 0, stream>>>
      (n, x, y, _r_grid);
  else
    _dot_xx_stage_1_of_2<block_size><<<grid_size, block_size, 0, stream>>>
      (n, x, _r_grid);

  CS_HIP_CHECK(hipGetLastError());

  cs_blas_hip_reduce_single_block<block_size, 1><<<1, block_size, 0, stream>>>
    (grid_size, _r_grid, _r_reduce);

  CS_HIP_CHECK(hipGetLastError());

  /* Need to synchronize stream in all cases so as to
     have up-to-date value in returned _r_reduce[0] */

  CS_HIP_CHECK(hipStreamSynchronize(stream));

  return _r_reduce[0];
}

#if defined(HAVE_HIPBLAS)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the absolute sum of vector values using hipBLAS.
 *
 * \param[in]  stream  associated HIP stream
 * \param[in]  n       size of arrays x and y
 * \param[in]  x       array of floating-point values (on device)
 *
 * \return  sum of absolute array values
 */
/*----------------------------------------------------------------------------*/

double
cs_blas_hipblas_asum(hipStream_t      stream,
                     cs_lnum_t        n,
                     const cs_real_t  x[])
{
  double dot = 0.;

  hipblasStatus_t status = HIPBLAS_STATUS_SUCCESS;

  _set_hipblas_stream(stream);

  if (sizeof(cs_real_t) == 8)
    status = hipblasDasum(_handle, n, (const double *)x, 1, &dot);
  else if (sizeof(cs_real_t) == 4) {
    float fdot;
    status = hipblasSasum(_handle, n, (const float *)x, 1, &fdot);
    dot = fdot;
  }

  _check_hipblas_status(__func__, __FILE__, __LINE__, status);

  return dot;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the dot product of 2 vectors: x.y using hipBLAS.
 *
 * \param[in]  stream  associated HIP stream
 * \param[in]  n       size of arrays x and y
 * \param[in]  x       array of floating-point values (on device)
 * \param[in]  y       array of floating-point values (on device)
 *
 * \return  dot product
 */
/*----------------------------------------------------------------------------*/

double
cs_blas_hipblas_dot(hipStream_t      stream,
                    cs_lnum_t        n,
                    const cs_real_t  x[],
                    const cs_real_t  y[])
{
  double dot = 0.;

  _set_hipblas_stream(stream);

  hipblasStatus_t status = HIPBLAS_STATUS_SUCCESS;
  if (sizeof(cs_real_t) == 8)
    hipblasDdot(_handle, n, (const double *)x, 1, (const double *)y, 1, &dot);
  else if (sizeof(cs_real_t) == 4) {
    float fdot;
    hipblasSdot(_handle, n, (const float *)x, 1, (const float *)y, 1, &fdot);
    dot = fdot;
  }
  _check_hipblas_status(__func__, __FILE__, __LINE__, status);

  return dot;
}

#endif // defined(HAVE_HIPBLAS)

/*----------------------------------------------------------------------------
 * Compute x <- alpha.x
 *
 * This function may be set to use either hipBLAS or a local kernel.
 *
 * \param[in]  stream  associated HIP stream
 * \param[in]  n       number of elements
 * \param[in]  alpha   constant value (on device)
 * \param[in]  x       vector of elements (on device)
 *----------------------------------------------------------------------------*/

void
cs_blas_hip_scal(hipStream_t       stream,
                 cs_lnum_t         n,
                 const cs_real_t  *alpha,
                 cs_real_t        *x)
{
#if defined(HAVE_HIPBLAS)

  if (_prefer_hipblas && _handle != nullptr) {
    _set_hipblas_stream(stream);

    hipblasStatus_t status = HIPBLAS_STATUS_SUCCESS;
    if (sizeof(cs_real_t) == 8)
      status = hipblasDscal(_handle, n,
                           (const double *)alpha,
                           (double *)x, 1);
    else if (sizeof(cs_real_t) == 4)
      status = hipblasSscal(_handle, n,
                           (const float *)alpha,
                           (float *)x, 1);
    _check_hipblas_status(__func__, __FILE__, __LINE__, status);

    return;
  }

#endif /* defined(HAVE_HIPBLAS) */

  /* If not using hipBLAS for this operation */

  const unsigned int blocksize = 256;  // try 640 also

  unsigned int gridsize = cs_hip_grid_size(n, blocksize);
  // gridsize = min(gridsize, 640;

  _scal<<<gridsize, blocksize, 0, stream>>>(n, alpha, x);
}

#endif /* defined(__HIPCC__)

/*----------------------------------------------------------------------------*/
