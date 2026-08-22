#ifndef CS_BLAS_CUDA_H
#define CS_BLAS_CUDA_H

/*============================================================================
 * BLAS (Basic Linear Algebra Subroutine) functions using CUDA.
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
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Finalize CUDA BLAS API.
 *
 * This frees resources such as the cuBLAS handle, if used.
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_blas_cuda_finalize(void);

#if defined(__CUDACC__)

/*----------------------------------------------------------------------------*/
/*
 * \brief Return the absolute sum of vector values using CUDA.
 *
 * \param[in]  stream  associated CUDA stream
 * \param[in]  n       size of array x
 * \param[in]  x       array of floating-point values (on device)
 *
 * \return  sum of absolute array values
 */
/*----------------------------------------------------------------------------*/

double
cs_blas_cuda_asum(cudaStream_t     stream,
                  cs_lnum_t        n,
                  const cs_real_t  x[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return the dot product of 2 vectors: x.y using CUDA.
 *
 * \param[in]  stream  associated CUDA stream
 * \param[in]  n       size of arrays x and y
 * \param[in]  x       array of floating-point values (on device)
 * \param[in]  y       array of floating-point values (on device)
 *
 * \return  dot product
 */
/*----------------------------------------------------------------------------*/

double
cs_blas_cuda_dot(cudaStream_t     stream,
                 cs_lnum_t        n,
                 const cs_real_t  x[],
                 const cs_real_t  y[]);

#if defined(HAVE_CUBLAS)

/*----------------------------------------------------------------------------*/
/*
 * \brief Return the absolute sum of vector values using cuBLAS.
 *
 * \param[in]  stream  associated CUDA stream
 * \param[in]  n       size of arrays x and y
 * \param[in]  x       array of floating-point values (on device)
 *
 * \return  sum of absolute array values
 */
/*----------------------------------------------------------------------------*/

double
cs_blas_cublas_asum(cudaStream_t     stream,
                    cs_lnum_t        n,
                    const cs_real_t  x[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return the dot product of 2 vectors: x.y using cuBLAS.
 *
 * \param[in]  stream  associated CUDA stream
 * \param[in]  n       size of arrays x and y
 * \param[in]  x       array of floating-point values (on device)
 * \param[in]  y       array of floating-point values (on device)
 *
 * \return  dot product
 */
/*----------------------------------------------------------------------------*/

double
cs_blas_cublas_dot(cudaStream_t     stream,
                   cs_lnum_t        n,
                   const cs_real_t  x[],
                   const cs_real_t  y[]);

#endif  /* defined(HAVE_CUBLAS) */

/*----------------------------------------------------------------------------
 * Compute x <- alpha.x
 *
 * This function may be set to use either cuBLAS or a local kernel.
 *
 * \param[in]  stream  associated CUDA stream
 * \param[in]  n       number of elements
 * \param[in]  alpha   constant value (on device)
 * \param[in]  x       vector of elements (on device)
 *----------------------------------------------------------------------------*/

void
cs_blas_cuda_scal(cudaStream_t      stream,
                  cs_lnum_t         n,
                  const cs_real_t  *alpha,
                  cs_real_t        *x);

#endif /* defined(__CUDACC__) */

/*----------------------------------------------------------------------------*/

#endif /* CS_BLAS_CUDA_H */
