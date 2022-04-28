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

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BLAS_CUDA_H__ */
