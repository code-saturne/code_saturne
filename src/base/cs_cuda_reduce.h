#ifndef CS_CUDA_REDUCE_H
#define CS_CUDA_REDUCE_H

/*============================================================================
 * BLAS (Basic Linear Algebra Subroutine) functions
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
 * External library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_base.h"
#include "base/cs_base_cuda.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Templated function definitions
 *============================================================================*/

#if defined(__CUDACC__)

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Kernel for sum reduction within a warp (for warp size 32).
 *
 * \tparam  stride     vector stride
 * \tparam  T          data type
 *
 * \param[in, out]  stmp  shared value to reduce
 * \param[in, out]  tid   thread id
 */
/*----------------------------------------------------------------------------*/

template <size_t stride, typename T>
__device__ static void __forceinline__
cs_cuda_reduce_warp_sum(volatile T  *stmp,
                        size_t       tid)
{
  T v[stride];
  #pragma unroll
  for (size_t i = 0; i < stride; i++)
    v[i] = stmp[tid*stride + i];
  __syncwarp();

  {
    #pragma unroll
    for (size_t i = 0; i < stride; i++)
      v[i] += stmp[(tid + 16)*stride + i];
    __syncwarp();
    #pragma unroll
    for (size_t i = 0; i < stride; i++)
      stmp[tid*stride + i] = v[i];
    __syncwarp();
  }
  {
    #pragma unroll
    for (size_t i = 0; i < stride; i++)
      v[i] += stmp[(tid + 8)*stride + i];
    __syncwarp();
    #pragma unroll
    for (size_t i = 0; i < stride; i++)
      stmp[tid*stride + i] = v[i];
    __syncwarp();
  }
  {
    #pragma unroll
    for (size_t i = 0; i < stride; i++)
      v[i] += stmp[(tid + 4)*stride + i];
    __syncwarp();
    #pragma unroll
    for (size_t i = 0; i < stride; i++)
      stmp[tid*stride + i] = v[i];
    __syncwarp();
  }
  {
    #pragma unroll
    for (size_t i = 0; i < stride; i++)
      v[i] += stmp[(tid + 2)*stride + i];
    __syncwarp();
    #pragma unroll
    for (size_t i = 0; i < stride; i++)
      stmp[tid*stride + i] = v[i];
    __syncwarp();
  }
  {
    #pragma unroll
    for (size_t i = 0; i < stride; i++)
      v[i] += stmp[(tid + 1)*stride + i];
    __syncwarp();
    #pragma unroll
    for (size_t i = 0; i < stride; i++)
      stmp[tid*stride + i] = v[i];
    __syncwarp();
  }
}

template <typename T>
__device__ static void __forceinline__
cs_cuda_reduce_warp_sum(volatile T  *stmp,
                        size_t       tid)
{
  constexpr unsigned int full_mask{0xffffffff};
  T sum = stmp[tid];
  for (size_t offset = CS_CUDA_WARP_SIZE/2; offset > 0; offset >>= 1)
    sum += __shfl_down_sync(full_mask, sum, offset);
  stmp[tid] = sum;
}

#if __CUDA_ARCH__ >= 800
// Specialize warpReduceFunc for int inputs to use __reduce_add_sync intrinsic
// when on SM 8.0 or higher

__device__ static void __forceinline__
cs_cuda_reduce_warp_sum(volatile int  *stmp,
                        size_t         tid)
{
  constexpr unsigned int full_mask{0xffffffff};
  int sum = __reduce_add_sync(full_mask, stmp[tid]);
  stmp[tid] = sum;
}
#endif

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
 * \tparam  blockSize  size of CUDA block
 * \tparam  stride     vector stride
 * \tparam  T          data type
 *
 * \param[in, out]  stmp       shared value to reduce
 * \param[in, out]  tid        thread id
 * \param[out]      sum_block  contribution of this block to the sum
 */
/*----------------------------------------------------------------------------*/

template <size_t blockSize, size_t stride, typename T>
__device__ static void __forceinline__
cs_cuda_reduce_block_sum(T       *stmp,
                         size_t   tid,
                         T       *sum_block)
{
  __syncthreads();

  /* Loop explicitely unrolled */

  if (stride == 1) { /* Scalar data */

    if (blockSize >= 1024) {
      if (tid < 512) {
        stmp[tid] += stmp[tid + 512];
      }
      __syncthreads();
    }
    if (blockSize >= 512) {
      if (tid < 256) {
        stmp[tid] += stmp[tid + 256];
      }
      __syncthreads();
    }
    if (blockSize >= 256) {
      if (tid < 128) {
        stmp[tid] += stmp[tid + 128];
      } __syncthreads();
    }
    if (blockSize >= 128) {
      if (tid <  64) {
        stmp[tid] += stmp[tid +  64];
      } __syncthreads();
    }
    if (blockSize >= 64) {
      if (tid <  32) {
        stmp[tid] += stmp[tid +  32];
      } __syncthreads();
    }

    if (tid < 32) {
      cs_cuda_reduce_warp_sum<stride>(stmp, tid);
    }

    // Output: b_res for this block

    if (tid == 0) sum_block[blockIdx.x] = stmp[0];

  }

  else { /* Vector data */

    if (blockSize >= 1024) {
      if (tid < 512) {
        #pragma unroll
        for (size_t i = 0; i < stride; i++)
          stmp[tid*stride + i] += stmp[(tid + 512)*stride + i];
      }
      __syncthreads();
    }
    if (blockSize >= 512) {
      if (tid < 256) {
        #pragma unroll
        for (size_t i = 0; i < stride; i++)
          stmp[tid*stride + i] += stmp[(tid + 256)*stride + i];
      }
      __syncthreads();
    }
    if (blockSize >= 256) {
      if (tid < 128) {
        #pragma unroll
        for (size_t i = 0; i < stride; i++)
          stmp[tid*stride + i] += stmp[(tid + 128)*stride + i];
      } __syncthreads();
    }
    if (blockSize >= 128) {
      if (tid <  64) {
        #pragma unroll
        for (size_t i = 0; i < stride; i++)
          stmp[tid*stride + i] += stmp[(tid + 64)*stride + i];
      } __syncthreads();
    }
    if (blockSize >= 64) {
      if (tid <  32) {
        #pragma unroll
        for (size_t i = 0; i < stride; i++)
          stmp[tid*stride + i] += stmp[(tid + 32)*stride + i];
      } __syncthreads();
    }

    if (tid < 32)
      cs_cuda_reduce_warp_sum<stride>(stmp, tid);

    // Output: b_res for this block

    if (tid == 0) {
      #pragma unroll
      for (size_t i = 0; i < stride; i++)
        sum_block[(blockIdx.x)*stride + i] = stmp[i];
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Kernel for sum reduction of vector data within a block.
 *
 * \tparam  stride  number of vector components
 *
 * \param[in, out]  n        number of vector values to reduce
 * \param[in, out]  g_idata  input array (size n)
 * \param[in, out]  g_odata  onput array (size 1)
 */
/*----------------------------------------------------------------------------*/

template <size_t blockSize, size_t stride, typename T>
__global__ static void
cs_cuda_reduce_sum_single_block(size_t   n,
                                T       *g_idata,
                                T       *g_odata)
{
  __shared__ T sdata[blockSize * stride];

  size_t tid = threadIdx.x;
  T r_s[stride];

  if (stride == 1) {

    r_s[0] = 0;
    sdata[tid] = 0.;

    for (size_t i = threadIdx.x; i < n; i+= blockSize)
      r_s[0] += g_idata[i];

    sdata[tid] = r_s[0];
    __syncthreads();

    for (size_t j = blockSize/2; j >= CS_CUDA_WARP_SIZE; j /= 2) {
      if (tid < j) {
        sdata[tid] += sdata[tid + j];
      }
      __syncthreads();
    }

    if (tid < 32)
      cs_cuda_reduce_warp_sum<stride>(sdata, tid);
    if (tid == 0) *g_odata = sdata[0];

  }
  else {

    #pragma unroll
    for (size_t k = 0; k < stride; k++) {
      r_s[k] = 0;
      sdata[tid*stride + k] = 0.;
    }

    for (size_t i = threadIdx.x; i < n; i+= blockSize) {
      #pragma unroll
      for (size_t k = 0; k < stride; k++) {
        r_s[k] += g_idata[i*stride + k];
      }
    }

    #pragma unroll
    for (size_t k = 0; k < stride; k++)
      sdata[tid*stride + k] = r_s[k];
    __syncthreads();

    for (size_t j = blockSize/2; j >= CS_CUDA_WARP_SIZE; j /= 2) {
      if (tid < j) {
        #pragma unroll
        for (size_t k = 0; k < stride; k++)
            sdata[tid*stride + k] += sdata[(tid + j)*stride + k];
      }
      __syncthreads();
    }

    if (tid < 32)
      cs_cuda_reduce_warp_sum<stride>(sdata, tid);
    if (tid == 0) {
      #pragma unroll
      for (size_t k = 0; k < stride; k++)
        g_odata[k] = sdata[k];
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Kernel for general reduction within a warp (for warp size 32).
 *
 * \tparam  R          reducer class
 * \tparam  T          data type
 *
 * \param[in, out]  stmp     shared value to reduce
 * \param[in, out]  tid      thread id
 */
/*----------------------------------------------------------------------------*/

template <typename R, typename T>
__device__ static void __forceinline__
cs_cuda_reduce_warp_reduce(T           *stmp,
                           size_t       tid)
{
  R reducer;

  T v = stmp[tid];
  __syncwarp();

  {
    reducer.combine(v, stmp[tid + 16]); __syncwarp();
    stmp[tid] = v;                      __syncwarp();
  }
  {
    reducer.combine(v, stmp[tid +  8]); __syncwarp();
    stmp[tid] = v;                      __syncwarp();
  }
  {
    reducer.combine(v, stmp[tid +  4]); __syncwarp();
    stmp[tid] = v;                      __syncwarp();
  }
  {
    reducer.combine(v, stmp[tid +  2]); __syncwarp();
    stmp[tid] = v;                      __syncwarp();
  }
  {
    reducer.combine(v, stmp[tid +  1]); __syncwarp();
    stmp[tid] = v;                      __syncwarp();
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute general reduction, combining result over all
 *        threads of a block.
 *
 * blockSize must be a power of 2.
 *
 * \tparam  blockSize  size of CUDA block
 * \tparam  R          reducer class
 * \tparam  T          data type
 *
 * \param[in, out]  stmp      shared value to reduce
 * \param[in, out]  tid       thread id
 * \param[out]      rd_block  contribution of this block to the reduction
 */
/*----------------------------------------------------------------------------*/

template <size_t blockSize, typename R, typename T>
__device__ static void __forceinline__
cs_cuda_reduce_block_reduce(T       *stmp,
                            size_t   tid,
                            T       *rd_block)
{
  R reducer;
  __syncthreads();

  /* Loop explicitely unrolled */

  if (blockSize >= 1024) {
    if (tid < 512) {
      reducer.combine(stmp[tid], stmp[tid + 512]);
    }
    __syncthreads();
  }
  if (blockSize >= 512) {
    if (tid < 256) {
      reducer.combine(stmp[tid], stmp[tid + 256]);
    }
    __syncthreads();
  }
  if (blockSize >= 256) {
    if (tid < 128) {
      reducer.combine(stmp[tid], stmp[tid + 128]);
    } __syncthreads();
  }
  if (blockSize >= 128) {
    if (tid <  64) {
      reducer.combine(stmp[tid], stmp[tid +  64]);
    } __syncthreads();
  }
  if (blockSize >= 64) {
    if (tid <  32) {
      reducer.combine(stmp[tid], stmp[tid +  32]);
    } __syncthreads();
  }

  if (tid < 32) {
    cs_cuda_reduce_warp_reduce<R>(stmp, tid);
  }

  // Output: rd_block for this block

  if (tid == 0) rd_block[blockIdx.x] = stmp[0];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Kernel for reduction of vector data within a block.
 *
 * \tparam  stride  number of vector components
 * \tparam  stride  number of vector components
 * \tparam  stride  number of vector components
 *
 * \param[in, out]  n        number of vector values to reduce
 * \param[in, out]  g_idata  input array (size n)
 * \param[in, out]  g_odata  onput array (size 1)
 */
/*----------------------------------------------------------------------------*/

template <size_t blockSize, typename R, typename T>
__global__ static void
cs_cuda_reduce_single_block(size_t   n,
                            T       *g_idata,
                            T       *g_odata)
{
  extern __shared__  int p_stmp[];
  T *sdata = reinterpret_cast<T *>(p_stmp);
  R  reducer;

  size_t tid = threadIdx.x;
  T r_s[1];

  reducer.identity(r_s[0]);
  reducer.identity(sdata[tid]);

  for (size_t i = threadIdx.x; i < n; i+= blockSize)
    reducer.combine(r_s[0], g_idata[i]);

  sdata[tid] = r_s[0];
  __syncthreads();

  for (size_t j = blockSize/2; j >= CS_CUDA_WARP_SIZE; j /= 2) {
    if (tid < j) {
      reducer.combine(sdata[tid], sdata[tid + j]);
    }
    __syncthreads();
  }

  if (tid < 32) cs_cuda_reduce_warp_reduce<R>(sdata, tid);
  if (tid == 0) *g_odata = sdata[0];
}

#endif /* defined(__CUDACC__) */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

#endif /* CS_CUDA_REDUCE_H */
