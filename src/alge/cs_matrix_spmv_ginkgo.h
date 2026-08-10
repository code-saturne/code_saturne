#ifndef CS_MATRIX_SPMV_GINKGO_H
#define CS_MATRIX_SPMV_GINKGO_H

/*============================================================================
 * Sparse Matrix SpMV operations using Ginkgo.
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
#include "base/cs_log.h"
#include "alge/cs_matrix.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print Ginkgo library information.
 *
 * \param[in]  log_type  destination log type
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_ginkgo_library_info(cs_log_t  log_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix-vector product y = A.x with CSR matrix, scalar Ginkgo version.
 *
 * \param[in]   matrix        pointer to matrix structure
 * \param[in]   exclude_diag  exclude diagonal if true
 * \param[in]   sync          synchronize ghost cells if true
 * \param[in]   x             multiplying vector values
 * \param[out]  y             resulting vector
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_ginkgo_csr(cs_matrix_t  *matrix,
                          bool          exclude_diag,
                          bool          sync,
                          cs_real_t     x[],
                          cs_real_t     y[]);

#if defined(HAVE_CUDA) || defined(HAVE_HIP)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix-vector product y = A.x with CSR matrix, scalar Ginkgo device.
 *
 * \param[in]   matrix        pointer to matrix structure
 * \param[in]   exclude_diag  exclude diagonal if true
 * \param[in]   sync          synchronize ghost cells if true
 * \param[in]   d_x           multiplying vector values (on device)
 * \param[out]  d_y           resulting vector (on device)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_ginkgo_device_csr(cs_matrix_t  *matrix,
                                 bool          exclude_diag,
                                 bool          sync,
                                 cs_real_t     d_x[],
                                 cs_real_t     d_y[]);

#endif /* defined(HAVE_CUDA) || defined(HAVE_HIP) */

END_C_DECLS

/*----------------------------------------------------------------------------*/

#endif /* CS_MATRIX_SPMV_GINKGO_H */
