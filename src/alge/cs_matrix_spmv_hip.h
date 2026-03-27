#ifndef __CS_MATRIX_SPMV_HIP_H__
#define __CS_MATRIX_SPMV_HIP_H__

/*============================================================================
 * Sparse Matrix-vector multiplication kernels using HIP.
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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*! \file
 *
 * \brief Sparse Matrix SpMV operations with HIP.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * Public HIP kernel prototypes
 *============================================================================*/

#if defined(__HIPCC__)

#endif /* defined(__HIPCC__) */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize HIP matrix API.
 *
 * This frees resources such as the rocSPARSE handle, if used.
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_hip_finalize(void);

#if defined(__HIPCC__)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign HIP stream for next HIP-based SpMV operations.
 *
 * If a stream other than the default stream (0) is used, it will not be
 * synchronized automatically after sparse matrix-vector products (so as to
 * avoid the corresponding overhead), so the caller will need to manage
 * stream syncronization manually.
 *
 * This function is callable only from HIP code.
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_hip_set_stream(hipStream_t  stream);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return stream used for HIP-based SpMV operations.
 *
 * This function is callable only from HIP code.
 */
/*----------------------------------------------------------------------------*/

hipStream_t
cs_matrix_spmv_hip_get_stream(void);

#endif /* defined(__HIPCC__) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with MSR matrix, scalar HIP version.
 *
 * \param[in]   matrix        pointer to matrix structure
 * \param[in]   exclude_diag  exclude diagonal if true,
 * \param[in]   sync          synchronize ghost cells if true
 * \param[in]   d_x           multipliying vector values (on device)
 * \param[out]  d_y           resulting vector (on device)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_hip_native(cs_matrix_t  *matrix,
                           bool          exclude_diag,
                           bool          sync,
                           cs_real_t     d_x[],
                           cs_real_t     d_y[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with CSR matrix, scalar HIP version.
 *
 * \param[in]   matrix        pointer to matrix structure
 * \param[in]   exclude_diag  exclude diagonal if true,
 * \param[in]   sync          synchronize ghost cells if true
 * \param[in]   d_x           multipliying vector values (on device)
 * \param[out]  d_y           resulting vector (on device)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_hip_csr(cs_matrix_t  *matrix,
                        bool          exclude_diag,
                        bool          sync,
                        cs_real_t     d_x[],
                        cs_real_t     d_y[]);

#if defined(HAVE_ROCSPARSE)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with CSR matrix, scalar rocSPARSE version.
 *
 * \param[in]   matrix        pointer to matrix structure
 * \param[in]   exclude_diag  exclude diagonal if true,
 * \param[in]   sync          synchronize ghost cells if true
 * \param[in]   d_x           multipliying vector values (on device)
 * \param[out]  d_y           resulting vector (on device)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_hip_csr_rocsparse(cs_matrix_t  *matrix,
                                 bool          exclude_diag,
                                 bool          sync,
                                 cs_real_t     d_x[],
                                 cs_real_t     d_y[]);

#endif /* defined(HAVE_ROCSPARSE) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with MSR matrix, scalar HIP version.
 *
 * \param[in]   matrix        pointer to matrix structure
 * \param[in]   exclude_diag  exclude diagonal if true,
 * \param[in]   sync          synchronize ghost cells if true
 * \param[in]   d_x           multipliying vector values (on device)
 * \param[out]  d_y           resulting vector (on device)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_hip_msr(cs_matrix_t  *matrix,
                        bool          exclude_diag,
                        bool          sync,
                        cs_real_t     d_x[],
                        cs_real_t     d_y[]);

#if defined(HAVE_ROCSPARSE)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with MSR matrix,
 *        scalar rocSPARSE version.
 *
 * \param[in]   matrix        pointer to matrix structure
 * \param[in]   exclude_diag  exclude diagonal if true,
 * \param[in]   sync          synchronize ghost cells if true
 * \param[in]   d_x           multipliying vector values (on device)
 * \param[out]  d_y           resulting vector (on device)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_hip_msr_rocsparse(cs_matrix_t  *matrix,
                                 bool          exclude_diag,
                                 bool          sync,
                                 cs_real_t     d_x[],
                                 cs_real_t     d_y[]);

#endif /* defined(HAVE_ROCSPARSE) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with MSR matrix, block diagonal
 *        HIP version.
 *
 * \param[in]   matrix        pointer to matrix structure
 * \param[in]   exclude_diag  exclude diagonal if true,
 * \param[in]   sync          synchronize ghost cells if true
 * \param[in]   d_x           multipliying vector values (on device)
 * \param[out]  d_y           resulting vector (on device)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_hip_msr_b(cs_matrix_t  *matrix,
                          bool          exclude_diag,
                          bool          sync,
                          cs_real_t     d_x[],
                          cs_real_t     d_y[]);

#if defined(HAVE_ROCSPARSE)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with MSR matrix, block diagonal
 *        rocSPARSE version.
 *
 * Remark: this functions is available with older rocSPARSE versions not
 *         providing the generic API, because they
 *         assume dense matrixes are always in column-major order, while
 *         row-major is needed with interleaved blocks.
 *
 * \param[in]   matrix        pointer to matrix structure
 * \param[in]   exclude_diag  exclude diagonal if true,
 * \param[in]   sync          synchronize ghost cells if true
 * \param[in]   d_x           multipliying vector values (on device)
 * \param[out]  d_y           resulting vector (on device)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_hip_msr_b_rocsparse(cs_matrix_t  *matrix,
                                   bool          exclude_diag,
                                   bool          sync,
                                   cs_real_t     d_x[],
                                   cs_real_t     d_y[]);

#endif /* defined(HAVE_ROCSPARSE) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with MSR matrix, block
 *        rocSPARSE version.
 *
 * Remmark: this functions is available with older rocSPARSE versions not
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
cs_matrix_spmv_hip_msr_bb_rocsparse(cs_matrix_t  *matrix,
                                    bool          exclude_diag,
                                    bool          sync,
                                    cs_real_t     d_x[],
                                    cs_real_t     d_y[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MATRIX_SPMV__HIP_H__ */
