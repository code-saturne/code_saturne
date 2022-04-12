#ifndef __CS_MATRIX_SPMV_CUDA_H__
#define __CS_MATRIX_SPMV_CUDA_H__

/*============================================================================
 * Sparse Matrix-vector multiplication kernels using CUDA.
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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*! \file cs_matrix_spmv_cuda.h
 *
 * \brief Sparse Matrix SpMV operations with CUDA.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * Public CUDA kernel prototypes
 *============================================================================*/

#if defined(__CUDACC__)

#endif /* defined(__CUDACC__) */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize CUDA matrix API.
 *
 * This frees resources such as the cuSPARSE handle, if used.
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_cuda_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with CSR matrix, scalar CUDA version.
 *
 * \param[in]   matrix        pointer to matrix structure
 * \param[in]   exclude_diag  exclude diagonal if true,
 * \param[in]   sync          synchronize ghost cells if true
 * \param[in]   x             multipliying vector values
 * \param[out]  y             resulting vector
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_cuda_p_l_csr(const cs_matrix_t  *matrix,
                            bool                exclude_diag,
                            bool                sync,
                            cs_real_t           x[restrict],
                            cs_real_t           y[restrict]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with CSR matrix, scalar cuSPARSE version.
 *
 * \param[in]   matrix        pointer to matrix structure
 * \param[in]   exclude_diag  exclude diagonal if true,
 * \param[in]   sync          synchronize ghost cells if true
 * \param[in]   x             multipliying vector values
 * \param[out]  y             resulting vector
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_cuda_p_l_csr_cusparse(cs_matrix_t  *matrix,
                                     bool          exclude_diag,
                                     bool          sync,
                                     cs_real_t     x[restrict],
                                     cs_real_t     y[restrict]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MATRIX_SPMV__CUDA_H__ */
