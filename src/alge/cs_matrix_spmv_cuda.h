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

#if defined(__CUDACC__)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign CUDA stream for next CUDA-based SpMV operations.
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
cs_matrix_spmv_cuda_set_stream(cudaStream_t  stream);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return stream used for CUDA-based SpMV operations.
 *
 * This function is callable only from CUDA code.
 */
/*----------------------------------------------------------------------------*/

cudaStream_t
cs_matrix_spmv_cuda_get_stream(void);

#endif /* defined(__CUDACC__) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with CSR matrix, scalar CUDA version.
 *
 * \param[in]   matrix        pointer to matrix structure
 * \param[in]   exclude_diag  exclude diagonal if true,
 * \param[in]   sync          synchronize ghost cells if true
 * \param[in]   d_x           multipliying vector values (on device)
 * \param[out]  d_y           resulting vector (on device)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_cuda_csr(const cs_matrix_t  *matrix,
                        bool                exclude_diag,
                        bool                sync,
                        cs_real_t           d_x[restrict],
                        cs_real_t           d_y[restrict]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with CSR matrix, scalar cuSPARSE version.
 *
 * \param[in]   matrix        pointer to matrix structure
 * \param[in]   exclude_diag  exclude diagonal if true,
 * \param[in]   sync          synchronize ghost cells if true
 * \param[in]   d_x           multipliying vector values (on device)
 * \param[out]  d_y           resulting vector (on device)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_cuda_csr_cusparse(cs_matrix_t  *matrix,
                                 bool          exclude_diag,
                                 bool          sync,
                                 cs_real_t     d_x[restrict],
                                 cs_real_t     d_y[restrict]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with MSR matrix, scalar CUDA version.
 *
 * \param[in]   matrix        pointer to matrix structure
 * \param[in]   exclude_diag  exclude diagonal if true,
 * \param[in]   sync          synchronize ghost cells if true
 * \param[in]   d_x           multipliying vector values (on device)
 * \param[out]  d_y           resulting vector (on device)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_cuda_msr(const cs_matrix_t  *matrix,
                        bool                exclude_diag,
                        bool                sync,
                        cs_real_t           d_x[restrict],
                        cs_real_t           d_y[restrict]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with MSR matrix, scalar cuSPARSE version.
 *
 * \param[in]   matrix        pointer to matrix structure
 * \param[in]   exclude_diag  exclude diagonal if true,
 * \param[in]   sync          synchronize ghost cells if true
 * \param[in]   d_x           multipliying vector values (on device)
 * \param[out]  d_y           resulting vector (on device)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_cuda_msr_cusparse(cs_matrix_t  *matrix,
                                 bool          exclude_diag,
                                 bool          sync,
                                 cs_real_t     d_x[restrict],
                                 cs_real_t     d_y[restrict]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with MSR matrix, block diagonal
 *        CUDA version.
 *
 * \param[in]   matrix        pointer to matrix structure
 * \param[in]   exclude_diag  exclude diagonal if true,
 * \param[in]   sync          synchronize ghost cells if true
 * \param[in]   d_x           multipliying vector values (on device)
 * \param[out]  d_y           resulting vector (on device)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_cuda_msr_b(cs_matrix_t  *matrix,
                          bool          exclude_diag,
                          bool          sync,
                          cs_real_t     d_x[restrict],
                          cs_real_t     d_y[restrict]);

#if defined(HAVE_CUSPARSE_GENERIC_API)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with MSR matrix, block diagonal
 *        cuSPARSE version.
 *
 * \param[in]   matrix        pointer to matrix structure
 * \param[in]   exclude_diag  exclude diagonal if true,
 * \param[in]   sync          synchronize ghost cells if true
 * \param[in]   d_x           multipliying vector values (on device)
 * \param[out]  d_y           resulting vector (on device)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_cuda_msr_b_cusparse(cs_matrix_t  *matrix,
                                   bool          exclude_diag,
                                   bool          sync,
                                   cs_real_t     d_x[restrict],
                                   cs_real_t     d_y[restrict]);

#endif /* defined(HAVE_CUSPARSE_GENERIC_API) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with MSR matrix, block
 *        cuSPARSE version.
 *
 * Remmark: this functions is available with older cuSPARSE versions not
 *          providing the generic API, because they
 *          assume dense matrixes are always in column-major order, while
 *          row-major is needed with interleaved blocks.
 *
 * \param[in]   matrix        pointer to matrix structure
 * \param[in]   exclude_diag  exclude diagonal if true,
 * \param[in]   sync          synchronize ghost cells if true
 * \param[in]   d_x           multipliying vector values (on device)
 * \param[out]  d_y           resulting vector (on device)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_cuda_msr_bb_cusparse(cs_matrix_t  *matrix,
                                    bool          exclude_diag,
                                    bool          sync,
                                    cs_real_t     d_x[restrict],
                                    cs_real_t     d_y[restrict]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MATRIX_SPMV__CUDA_H__ */
