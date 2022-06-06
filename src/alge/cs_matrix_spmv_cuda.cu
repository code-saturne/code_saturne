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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#if defined(HAVE_CUSPARSE)
#include <cusparse_v2.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_base_cuda.h"
#include "cs_cuda_contrib.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_matrix.h"
#include "cs_matrix_priv.h"
#include "cs_matrix_spmv.h"

/*----------------------------------------------------------------------------*/
/*! \file cs_matrix_spmv_cuda.c
 *
 * \brief Sparse Matrix SpMV operations with CUDA.
 */
/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compatibility macro for __ldg (load from generic memory) intrinsic,
 * forcing load from read-only texture cache.
 *
 * This was not available in (very old) CUDA architectures.
 *----------------------------------------------------------------------------*/

#if __CUDA_ARCH__ < 350
#define __ldg(ptr) *(ptr);
#endif

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

#if defined(HAVE_CUSPARSE)

/* Mapping of matrix coefficients and structure to cuSPARSE */
/*----------------------------------------------------------*/

typedef struct _cs_matrix_cusparse_map_t {

#if defined(HAVE_CUSPARSE_GENERIC_API)

  bool  block_diag;             /* Use identity blocks diagonal structure ? */

  cusparseSpMatDescr_t  matA;   /* Handle to cusparse Matrix */

  cusparseDnMatDescr_t  matX;   /* Handle to cusparse Matrix (blocked vector) */
  cusparseDnMatDescr_t  matY;   /* Handle to cusparse Matrix (blocked vector) */

  cusparseDnVecDescr_t  vecX;   /* Handle to cusparse Vector */
  cusparseDnVecDescr_t  vecY;   /* Handle to cusparse output Vector */

  void  *vecXValues;            /* Pointer to vector values */
  void  *vecYValues;            /* Pointer to vector values */

  void  *dBuffer;               /* Associated buffer */

#endif

  /* When not using generic API */

  int  nnz;                   /* Number of nonzeroes */
  cusparseMatDescr  *descrA;  /* Handle to cusparse Matrix description */

  void  *d_row_index;         /* Pointer to row index */
  void  *d_col_id;            /* Pointer to column ids */
  void  *d_e_val;             /* Pointer to matrix extradiagonal values */

} cs_matrix_cusparse_map_t;

#endif // defined(HAVE_CUSPARSE)

/*============================================================================
 *  Global variables
 *============================================================================*/

static cudaStream_t _stream = 0;

#if defined(HAVE_CUSPARSE)

static cusparseHandle_t  _handle = NULL;

#endif

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/* \brief Local matrix.vector product y = A.x with CSR matrix arrays.
 *
 * \param[in]   n_rows     number of local rows
 * \param[in]   row_index  pointer to matrix rows index
 * \param[in]   col_id     pointer to matrix column id
 * \param[in]   val        pointer to matrix values
 * \param[in]   x          multipliying vector values
 * \param[out]  y          resulting vector
 */
/*----------------------------------------------------------------------------*/

__global__ static void
_mat_vect_p_l_csr(cs_lnum_t         n_rows,
                  const cs_lnum_t  *__restrict__ row_index,
                  const cs_lnum_t  *__restrict__ col_id,
                  const cs_real_t  *__restrict__ val,
                  const cs_real_t  *__restrict__ x,
                  cs_real_t        *__restrict__ y)
{
  cs_lnum_t ii = blockIdx.x * blockDim.x + threadIdx.x;
  cs_lnum_t jj;

  if (ii < n_rows) {
    cs_real_t sii = 0.0;
    const cs_lnum_t *__restrict__ _col_id = col_id + row_index[ii];
    const cs_real_t *__restrict__ m_row  = val + row_index[ii];
    cs_lnum_t n_cols = row_index[ii + 1] - row_index[ii];
#pragma unroll
    for (jj = 0; jj < n_cols; jj++) {
      sii += m_row[jj] * __ldg(x + _col_id[jj]);
    }
    y[ii] = sii;
  }
}

/*----------------------------------------------------------------------------*/
/* \brief Local matrix.vector product y = A.x with CSR matrix arrays,
 *        excluding diagonal part.
 *
 * \param[in]   n_rows     number of local rows
 * \param[in]   row_index  pointer to matrix rows index
 * \param[in]   col_id     pointer to matrix column id
 * \param[in]   val        pointer to matrix values
 * \param[in]   x          multipliying vector values
 * \param[out]  y          resulting vector
 */
/*----------------------------------------------------------------------------*/

__global__ static void
_mat_vect_p_l_csr_exdiag(cs_lnum_t         n_rows,
                         const cs_lnum_t  *__restrict__ row_index,
                         const cs_lnum_t  *__restrict__ col_id,
                         const cs_real_t  *__restrict__ val,
                         const cs_real_t  *__restrict__ x,
                         cs_real_t        *__restrict__ y)
{
  cs_lnum_t ii = blockIdx.x * blockDim.x + threadIdx.x;

  if (ii < n_rows) {
    cs_real_t        sii            = 0.0;
    const cs_lnum_t *__restrict__ _col_id = col_id + row_index[ii];
    const cs_real_t *__restrict__ m_row  = val + row_index[ii];
    cs_lnum_t n_cols = row_index[ii + 1] - row_index[ii];
#pragma unroll
    for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
      cs_lnum_t c_id = _col_id[jj];
      if (c_id != ii)
        sii += m_row[jj] * __ldg(x + c_id);
    }
    y[ii] = sii;
  }
}

/*----------------------------------------------------------------------------*/
/* \brief Substract local diagonal contribution with CSR matrix arrays.
 *
 * \param[in]   n_rows     number of local rows
 * \param[in]   row_index  pointer to matrix rows index
 * \param[in]   col_id     pointer to matrix column id
 * \param[in]   val        pointer to matrix values
 * \param[in]   x          multipliying vector values
 * \param[out]  y          resulting vector
 */
/*----------------------------------------------------------------------------*/

__global__ static void
_mat_vect_p_l_csr_substract_diag(cs_lnum_t         n_rows,
                                 const cs_lnum_t  *__restrict__ row_index,
                                 const cs_lnum_t  *__restrict__ col_id,
                                 const cs_real_t  *__restrict__ val,
                                 const cs_real_t  *__restrict__ x,
                                 cs_real_t        *__restrict__ y)
{
  cs_lnum_t ii = blockIdx.x * blockDim.x + threadIdx.x;

  if (ii < n_rows) {
    const cs_lnum_t *__restrict__ _col_id = col_id + row_index[ii];
    const cs_real_t *__restrict__ m_row  = val + row_index[ii];
    cs_lnum_t n_cols = row_index[ii + 1] - row_index[ii];
#pragma unroll
    for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
      cs_lnum_t c_id = _col_id[jj];
      if (c_id == ii) {
        y[ii] -= m_row[jj] * x[ii];
        break;
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/* \brief Local matrix.vector product y = A.x with MSR matrix arrays.
 *
 * \param[in]   n_rows     number of local rows
 * \param[in]   row_index  pointer to matrix rows index
 * \param[in]   col_id     pointer to matrix column id
 * \param[in]   d_val      pointer to diagonal matrix values
 * \param[in]   x_val      pointer to extradiagonal matrix values
 * \param[in]   x          multipliying vector values
 * \param[out]  y          resulting vector
 */
/*----------------------------------------------------------------------------*/

__global__ static void
_mat_vect_p_l_msr(cs_lnum_t         n_rows,
                  const cs_lnum_t  *__restrict__ row_index,
                  const cs_lnum_t  *__restrict__ col_id,
                  const cs_real_t  *__restrict__ d_val,
                  const cs_real_t  *__restrict__ x_val,
                  const cs_real_t  *__restrict__ x,
                  cs_real_t        *__restrict__ y)
{
  cs_lnum_t ii = blockIdx.x * blockDim.x + threadIdx.x;

  if (ii < n_rows) {
    const cs_lnum_t *__restrict__ _col_id = col_id + row_index[ii];
    const cs_real_t *__restrict__ m_row  = x_val + row_index[ii];

    cs_lnum_t n_cols = row_index[ii + 1] - row_index[ii];

    cs_real_t sii = 0.0;

    for (cs_lnum_t jj = 0; jj < n_cols; jj++)
      sii += m_row[jj] * __ldg(x + _col_id[jj]);

    y[ii] = sii + d_val[ii] * x[ii];
  }
}

/*----------------------------------------------------------------------------*/
/* \brief Local diagonal contribution y = Da.x  + y.
 *
 * This can be combined with a cuSPARSE CSR SpMV product with the
 * extra-diagonal portion of an MSR or distributed matrix.
 *
 * \param[in]       n_rows  number of local rows
 * \param[in]       d_val   pointer to diagonal matrix values
 * \param[in]       x       multipliying vector values
 * \param[in, out]  y       resulting vector
 */
/*----------------------------------------------------------------------------*/

__global__ static void
_mat_vect_p_l_msr_diag(cs_lnum_t         n_rows,
                       const cs_real_t  *__restrict__ d_val,
                       const cs_real_t  *__restrict__ x,
                       cs_real_t        *__restrict__ y)
{
  cs_lnum_t ii = blockIdx.x * blockDim.x + threadIdx.x;

  if (ii < n_rows)
    y[ii] = d_val[ii] * x[ii];
}

/*----------------------------------------------------------------------------*/
/* \brief Local matrix.vector product y = A.x with MSR matrix,
 *        3x3 blocked diagonal version.
 *
 * \param[in]   n_rows     number of local rows
 * \param[in]   d_val      pointer to diagonal matrix values
 * \param[in]   x          multipliying vector values
 * \param[out]  y          resulting vector
 */
/*----------------------------------------------------------------------------*/

__global__ static void
_b_3_3_spmv_diag(cs_lnum_t        n_rows,
                 const cs_real_t  *__restrict__ d_val,
                 const cs_real_t  *__restrict__ x,
                 cs_real_t        *__restrict__ y)
{
  cs_lnum_t ii = blockIdx.x * blockDim.x + threadIdx.x;
  if (ii < n_rows) {

#   pragma unroll
    for (cs_lnum_t kk = 0; kk < 3; kk++) {
      y[ii*3 + kk] =   d_val[ii * 9 + kk * 3]     * x[ii * 3]
                     + d_val[ii * 9 + kk * 3 + 1] * x[ii * 3 + 1]
                     + d_val[ii * 9 + kk * 3 + 2] * x[ii * 3 + 2];
    }

  }
}

/*----------------------------------------------------------------------------*/
/* \brief Local matrix.vector product y = A.x with MSR matrix,
 *        templated blocked diagonal version.
 *
 * \param[in]   n_rows     number of local rows
 * \param[in]   d_val      pointer to diagonal matrix values
 * \param[in]   x          multipliying vector values
 * \param[out]  y          resulting vector
 */
/*----------------------------------------------------------------------------*/

template <const int n>
__global__ static void
_b_spmv_diag(cs_lnum_t        n_rows,
             const cs_real_t  *__restrict__ d_val,
             const cs_real_t  *__restrict__ x,
             cs_real_t        *__restrict__ y)
{
  cs_lnum_t ii = blockIdx.x * blockDim.x + threadIdx.x;

  if (ii < n_rows) {
    const cs_lnum_t nn = n*n;

    cs_real_t sii[n];

    for (cs_lnum_t kk = 0; kk < n; kk++)
      sii[kk] = 0.;

    for (cs_lnum_t kk = 0; kk < n; kk++) {
      for (cs_lnum_t ll = 0; ll < n; ll++) {
        sii[kk] += d_val[ii*nn + kk*n + ll] * x[ii*n + ll];
      }
    }

    for (cs_lnum_t kk = 0; kk < n; kk++)
      y[ii*n + kk] = sii[kk];
  }
}

/*----------------------------------------------------------------------------*/
/* \brief Local matrix.vector product y = A.x with MSR matrix,
 *        3x3 blocked diagonal version.
 *
 * \param[in]   n_rows     number of local rows
 * \param[in]   row_index  pointer to matrix rows index
 * \param[in]   col_id     pointer to matrix column id
 * \param[in]   d_val      pointer to diagonal matrix values
 * \param[in]   x_val      pointer to extradiagonal matrix values
 * \param[in]   x          multipliying vector values
 * \param[out]  y          resulting vector
 */
/*----------------------------------------------------------------------------*/

__global__ static void
_b_3_3_mat_vect_p_l_msr(cs_lnum_t        n_rows,
                        const cs_lnum_t  *__restrict__ col_id,
                        const cs_lnum_t  *__restrict__ row_index,
                        const cs_real_t  *__restrict__ d_val,
                        const cs_real_t  *__restrict__ x_val,
                        const cs_real_t  *__restrict__ x,
                        cs_real_t        *__restrict__ y)
{
  cs_lnum_t ii = blockIdx.x * blockDim.x + threadIdx.x;

  if (ii < n_rows) {
    const cs_lnum_t *__restrict__ _col_id = col_id + row_index[ii];
    const cs_real_t *__restrict__ m_row  = x_val + row_index[ii];
    cs_lnum_t n_cols = row_index[ii + 1] - row_index[ii];
    cs_real_t sii[3];
    for (cs_lnum_t kk = 0; kk < 3; kk++) {
      sii[kk] =   d_val[ii * 9 + kk * 3]     * x[ii * 3]
                + d_val[ii * 9 + kk * 3 + 1] * x[ii * 3 + 1]
                + d_val[ii * 9 + kk * 3 + 2] * x[ii * 3 + 2];
    }

    for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
      for (cs_lnum_t kk = 0; kk < 3; kk++)
        sii[kk] += m_row[jj] * __ldg(x + (_col_id[jj]*3 + kk));
    }

    y[ii*3]     = sii[0];
    y[ii*3 + 1] = sii[1];
    y[ii*3 + 2] = sii[2];
  }
}

/*----------------------------------------------------------------------------*/
/* \brief Local matrix.vector product y = A.x with MSR matrix,
 *        excluding 3x3 blocked diagonal.
 *
 * \param[in]   n_rows     number of local rows
 * \param[in]   row_index  pointer to matrix rows index
 * \param[in]   col_id     pointer to matrix column id
 * \param[in]   d_val      pointer to diagonal matrix values
 * \param[in]   x_val      pointer to extradiagonal matrix values
 * \param[in]   x          multipliying vector values
 * \param[out]  y          resulting vector
 */
/*----------------------------------------------------------------------------*/

__global__ static void
_b_3_3_mat_vect_p_l_msr_exdiag(cs_lnum_t        n_rows,
                               const cs_lnum_t  *__restrict__ col_id,
                               const cs_lnum_t  *__restrict__ row_index,
                               const cs_real_t  *__restrict__ d_val,
                               const cs_real_t  *__restrict__ x_val,
                               const cs_real_t  *__restrict__ x,
                               cs_real_t        *__restrict__ y)
{
  cs_lnum_t ii = blockIdx.x * blockDim.x + threadIdx.x;

  if (ii < n_rows) {
    const cs_lnum_t *__restrict__ _col_id = col_id + row_index[ii];
    const cs_real_t *__restrict__ m_row  = x_val + row_index[ii];
    cs_lnum_t n_cols = row_index[ii + 1] - row_index[ii];
    cs_real_t sii[3];
    for (cs_lnum_t kk = 0; kk < 3; kk++)
      sii[kk] = 0.;

    for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
      for (cs_lnum_t kk = 0; kk < 3; kk++)
        sii[kk] += m_row[jj] * __ldg(x + (_col_id[jj]*3 + kk));
    }

    y[ii * 3]     = sii[0];
    y[ii * 3 + 1] = sii[1];
    y[ii * 3 + 2] = sii[2];
  }
}

/*----------------------------------------------------------------------------*/
/* \brief Local matrix.vector product y = A.x with MSR matrix,
 *        blocked diagonal version.
 *
 * \param[in]   n_rows     number of local rows
 * \param[in]   row_index  pointer to matrix rows index
 * \param[in]   col_id     pointer to matrix column id
 * \param[in]   d_val      pointer to diagonal matrix values
 * \param[in]   x_val      pointer to extradiagonal matrix values
 * \param[in]   x          multipliying vector values
 * \param[out]  y          resulting vector
 */
/*----------------------------------------------------------------------------*/

template <const int n>
__global__ static void
_b_mat_vect_p_l_msr(cs_lnum_t        n_rows,
                    const cs_lnum_t  *__restrict__ col_id,
                    const cs_lnum_t  *__restrict__ row_index,
                    const cs_real_t  *__restrict__ d_val,
                    const cs_real_t  *__restrict__ x_val,
                    const cs_real_t  *__restrict__ x,
                    cs_real_t        *__restrict__ y)
{
  cs_lnum_t ii = blockIdx.x * blockDim.x + threadIdx.x;

  if (ii < n_rows) {
    const cs_lnum_t nn = n*n;

    const cs_lnum_t *__restrict__ _col_id = col_id + row_index[ii];
    const cs_real_t *__restrict__ m_row  = x_val + row_index[ii];
    cs_lnum_t n_cols = row_index[ii + 1] - row_index[ii];
    cs_real_t sii[n];

    for (cs_lnum_t kk = 0; kk < n; kk++)
      sii[kk] = 0.;

    for (cs_lnum_t kk = 0; kk < n; kk++) {
      for (cs_lnum_t ll = 0; ll < n; ll++) {
        sii[kk] += d_val[ii*nn + kk*n + ll] * x[ii*n + ll];
      }
    }

    for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
      for (cs_lnum_t kk = 0; kk < n; kk++)
        sii[kk] += m_row[jj] * __ldg(x + (_col_id[jj]*n + kk));
    }

    for (cs_lnum_t kk = 0; kk < n; kk++)
      y[ii*n + kk] = sii[kk];
  }
}

/*----------------------------------------------------------------------------*/
/* \brief Local matrix.vector product y = A.x with MSR matrix,
 *        excluding blocked diagonal.
 *
 * \param[in]   n_rows     number of local rows
 * \param[in]   row_index  pointer to matrix rows index
 * \param[in]   col_id     pointer to matrix column id
 * \param[in]   d_val      pointer to diagonal matrix values
 * \param[in]   x_val      pointer to extradiagonal matrix values
 * \param[in]   x          multipliying vector values
 * \param[out]  y          resulting vector
 */
/*----------------------------------------------------------------------------*/

template <const int n>
__global__ static void
_b_mat_vect_p_l_msr_exdiag(cs_lnum_t        n_rows,
                           const cs_lnum_t  *__restrict__ col_id,
                           const cs_lnum_t  *__restrict__ row_index,
                           const cs_real_t  *__restrict__ d_val,
                           const cs_real_t  *__restrict__ x_val,
                           const cs_real_t  *__restrict__ x,
                           cs_real_t        *__restrict__ y)
{
  cs_lnum_t ii = blockIdx.x * blockDim.x + threadIdx.x;

  if (ii < n_rows) {
    const cs_lnum_t *__restrict__ _col_id = col_id + row_index[ii];
    const cs_real_t *__restrict__ m_row  = x_val + row_index[ii];
    cs_lnum_t n_cols = row_index[ii + 1] - row_index[ii];
    cs_real_t sii[n];
    for (cs_lnum_t kk = 0; kk < n; kk++)
      sii[kk] = 0.;

    for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
      for (cs_lnum_t kk = 0; kk < n; kk++)
        sii[kk] += m_row[jj] * __ldg(x + (_col_id[jj]*n + kk));
    }

    for (cs_lnum_t kk = 0; kk < n; kk++)
      y[ii*n + kk] = sii[kk];
  }
}

/*----------------------------------------------------------------------------
 * Start synchronization of ghost values prior to matrix.vector product.
 *
 * Values are packed on the device, so:
 * - If MPI is CUDA-aware, no values need to go through the host
 * - Otherwise, only halo values need to go through the host, not the
 *   whole array.
 *
 * parameters:
 *   matrix   <-- pointer to matrix structure
 *   d_x      <-> multipliying vector values (ghost values updated)
 *
 * returns:
 *   halo state to use for synchronisation finalisation.
 *----------------------------------------------------------------------------*/

static cs_halo_state_t *
_pre_vector_multiply_sync_x_start(const cs_matrix_t   *matrix,
                                  cs_real_t            d_x[restrict])
{
 cs_halo_state_t *hs = NULL;

  if (matrix->halo != NULL) {

    if (_stream != 0)
      cudaStreamSynchronize(_stream);

    hs = cs_halo_state_get_default();

    cs_halo_sync_pack_d(matrix->halo,
                        CS_HALO_STANDARD,
                        CS_REAL_TYPE,
                        matrix->db_size,
                        d_x,
                        NULL,
                        hs);

    cs_halo_sync_start(matrix->halo, d_x, hs);

  }

  return hs;
}

#if defined(HAVE_CUSPARSE)

/*----------------------------------------------------------------------------
 * Unset matrix cuSPARSE mapping.
 *
 * parameters:
 *   matrix    <-- pointer to matrix structure
 *----------------------------------------------------------------------------*/

static void
_unset_cusparse_map(cs_matrix_t   *matrix)
{
  cs_matrix_cusparse_map_t *csm
    = (cs_matrix_cusparse_map_t *)matrix->ext_lib_map;

  if (csm == NULL)
    return;

#if defined(HAVE_CUSPARSE_GENERIC_API)

  cusparseDestroySpMat(csm->matA);

  if (csm->block_diag == false) {
    if (csm->vecXValues != NULL)
      cusparseDestroyDnVec(csm->vecX);
    if (csm->vecYValues != NULL)
      cusparseDestroyDnVec(csm->vecY);
  }
  else {
    if (csm->vecXValues != NULL)
      cusparseDestroyDnMat(csm->matX);
    if (csm->vecYValues != NULL)
      cusparseDestroyDnMat(csm->matY);
  }

  if (csm->dBuffer != NULL) {
    CS_CUDA_CHECK(cudaFree(csm->dBuffer));
    csm->dBuffer = NULL;
  }

  csm->block_diag = false;

  csm->vecXValues = NULL;
  csm->vecYValues = NULL;

#else

  cusparseDestroyMatDescr(csm->descrA);

#endif

  csm->nnz = 0;
  csm->d_row_index = NULL;
  csm->d_col_id = NULL;
  csm->d_e_val = NULL;

  BFT_FREE(matrix->ext_lib_map);
  matrix->destroy_adaptor = NULL;
}

/*----------------------------------------------------------------------------
 * Set matrix cuSPARSE mapping.
 *
 * parameters:
 *   matrix    <-- pointer to matrix structure
 *----------------------------------------------------------------------------*/

static cs_matrix_cusparse_map_t *
_set_cusparse_map(cs_matrix_t   *matrix)
{
  cs_matrix_cusparse_map_t *csm
    = (cs_matrix_cusparse_map_t *)matrix->ext_lib_map;

  if (csm != NULL) {
    _unset_cusparse_map(matrix);
  }
  else {
    BFT_MALLOC(csm, 1, cs_matrix_cusparse_map_t);
    matrix->ext_lib_map = (void *)csm;
  }
  matrix->destroy_adaptor = _unset_cusparse_map;

  const void *row_index, *col_id;
  const void *e_val;
  cs_lnum_t nnz = 0;

  if (matrix->type == CS_MATRIX_CSR) {
    const cs_matrix_struct_csr_t *ms
      = (const cs_matrix_struct_csr_t  *)matrix->structure;
    const cs_matrix_coeff_csr_t *mc
      = (const cs_matrix_coeff_csr_t *)matrix->coeffs;
    nnz = ms->row_index[matrix->n_rows];
    row_index = cs_get_device_ptr_const_pf
                  (const_cast<cs_lnum_t *>(ms->row_index));
    col_id = cs_get_device_ptr_const_pf
               (const_cast<cs_lnum_t *>(ms->col_id));
    e_val = cs_get_device_ptr_const_pf
              (const_cast<cs_real_t *>(mc->val));
  }
  else {
    const cs_matrix_struct_dist_t *ms
      = (const cs_matrix_struct_dist_t *)matrix->structure;
    const cs_matrix_coeff_dist_t *mc
      = (const cs_matrix_coeff_dist_t *)matrix->coeffs;
    nnz = ms->e.row_index[matrix->n_rows];
    row_index = cs_get_device_ptr_const_pf
                  (const_cast<cs_lnum_t *>(ms->e.row_index));
    col_id = cs_get_device_ptr_const_pf
               (const_cast<cs_lnum_t *>(ms->e.col_id));
    e_val = cs_get_device_ptr_const_pf
              (const_cast<cs_real_t *>(mc->e_val));
  }

  cusparseStatus_t status = CUSPARSE_STATUS_SUCCESS;

  if (_handle == NULL)
    status = cusparseCreate(&_handle);

#if defined(HAVE_CUSPARSE_GENERIC_API)

  if (matrix->db_size > matrix->eb_size)
    csm->block_diag = true;
  else
    csm->block_diag = false;

  if (CUSPARSE_STATUS_SUCCESS != status)
    bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
              __func__, cusparseGetErrorString(status));

  if (matrix->eb_size == 1) {

    cusparseIndexType_t index_dtype
      = (sizeof(cs_lnum_t) == 4) ? CUSPARSE_INDEX_32I : CUSPARSE_INDEX_64I;
    cudaDataType_t val_dtype
      = (sizeof(cs_real_t) == 8) ? CUDA_R_64F : CUDA_R_32F;

    csm->vecXValues = NULL;  /* Pointer to vector values */
    csm->vecYValues = NULL;  /* Pointer to vector values */
    csm->dBuffer = NULL;

    status = cusparseCreateCsr(&(csm->matA),
                               matrix->n_rows,
                               matrix->n_cols_ext,
                               nnz,
                               const_cast<void *>(row_index),
                               const_cast<void *>(col_id),
                               const_cast<void *>(e_val),
                               index_dtype,
                               index_dtype,
                               CUSPARSE_INDEX_BASE_ZERO,
                               val_dtype);

  }
  else {

    csm->nnz = nnz;
    csm->d_e_val = const_cast<void *>(e_val);

    csm->d_row_index = const_cast<void *>(row_index);
    csm->d_col_id = const_cast<void *>(col_id);

    status = cusparseCreateMatDescr(&(csm->descrA));

  }

  if (CUSPARSE_STATUS_SUCCESS != status)
    bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
              __func__, cusparseGetErrorString(status));

#else

  if (CUSPARSE_STATUS_SUCCESS != status)
    bft_error(__FILE__, __LINE__, 0, _("%s: cuSPARSE error %d."),
              __func__, (int)status);

  csm->nnz = nnz;
  csm->d_e_val = const_cast<void *>(e_val);

  csm->d_row_index = const_cast<void *>(row_index);
  csm->d_col_id = const_cast<void *>(col_id);

  status = cusparseCreateMatDescr(&(csm->descrA));

  if (CUSPARSE_STATUS_SUCCESS != status)
    bft_error(__FILE__, __LINE__, 0, _("%s: cuSPARSE error %d."),
              __func__, (int)status);

  cusparseSetMatIndexBase(csm->descrA, CUSPARSE_INDEX_BASE_ZERO);
  cusparseSetMatType(csm->descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatDiagType(csm->descrA, CUSPARSE_DIAG_TYPE_NON_UNIT);

#endif

  return csm;
}

/*----------------------------------------------------------------------------
 * Update matrix cuSPARSE mapping.
 *
 * parameters:
 *   csm       <-> cuSPARSE matrix mapping
 *   matrix    <-- pointer to matrix structure
 *   d_x       <-- pointer to input vector (on device)
 *   d_y       <-- pointer to output vector (on device)
 *----------------------------------------------------------------------------*/

static void
_update_cusparse_map(cs_matrix_cusparse_map_t  *csm,
                     const cs_matrix_t         *matrix,
                     void                      *d_x,
                     void                      *d_y)
{
  assert(csm != NULL);

#if defined(HAVE_CUSPARSE_GENERIC_API)

  cusparseStatus_t status = CUSPARSE_STATUS_SUCCESS;
  cudaDataType_t val_dtype
    = (sizeof(cs_real_t) == 8) ? CUDA_R_64F : CUDA_R_32F;

  if (d_x != csm->vecXValues) {
    if (csm->vecXValues != NULL)
      cusparseDestroyDnVec(csm->vecX);

    status = cusparseCreateDnVec(&(csm->vecX),
                                 matrix->n_cols_ext,
                                 d_x,
                                 val_dtype);

    if (CUSPARSE_STATUS_SUCCESS != status)
      bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
                __func__, cusparseGetErrorString(status));

    csm->vecXValues = d_x;
  }

  if (d_y != csm->vecYValues) {
    if (csm->vecYValues != NULL)
      cusparseDestroyDnVec(csm->vecY);

    status = cusparseCreateDnVec(&(csm->vecY),
                                 matrix->n_rows,
                                 d_y,
                                 val_dtype);

    if (CUSPARSE_STATUS_SUCCESS != status)
      bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
                __func__, cusparseGetErrorString(status));

    csm->vecYValues = d_y;
  }

  if (csm->dBuffer == NULL) {
    size_t bufferSize = 0;
    cs_real_t alpha = 1.0;
    cs_real_t beta = 1.0;  /* 0 should be enough for SmPV, 1 needed for
                              y = A.x + b.y
                              which is useful when y is initialized by
                              a separate diagonal da.x product */

    status = cusparseSpMV_bufferSize(_handle,
                                     CUSPARSE_OPERATION_NON_TRANSPOSE,
                                     &alpha,
                                     csm->matA,
                                     csm->vecX,
                                     &beta,
                                     csm->vecY,
                                     val_dtype,
                                     CUSPARSE_MV_ALG_DEFAULT,
                                     &bufferSize);

    CS_CUDA_CHECK(cudaMalloc(&(csm->dBuffer), bufferSize));
  }

#endif
}

/*----------------------------------------------------------------------------
 * Update matrix cuSPARSE mapping in block diagonal case.
 *
 * parameters:
 *   csm       <-> cuSPARSE matrix mapping
 *   matrix    <-- pointer to matrix structure
 *   d_x       <-- pointer to input vector (on device)
 *   d_y       <-- pointer to output vector (on device)
 *----------------------------------------------------------------------------*/

static void
_update_cusparse_map_block_diag(cs_matrix_cusparse_map_t  *csm,
                                const cs_matrix_t         *matrix,
                                void                      *d_x,
                                void                      *d_y)
{
  assert(csm != NULL);

#if defined(HAVE_CUSPARSE_GENERIC_API)

  cusparseStatus_t status = CUSPARSE_STATUS_SUCCESS;
  cudaDataType_t val_dtype
    = (sizeof(cs_real_t) == 8) ? CUDA_R_64F : CUDA_R_32F;

  if (d_x != csm->vecXValues) {
    if (csm->vecXValues != NULL)
      cusparseDestroyDnMat(csm->matX);

    status = cusparseCreateDnMat(&(csm->matX),
                                 matrix->n_cols_ext,
                                 matrix->db_size,
                                 matrix->db_size, // n_cols_ext, /* leading dimension */
                                 d_x,
                                 val_dtype,
                                 CUSPARSE_ORDER_ROW);

    if (CUSPARSE_STATUS_SUCCESS != status)
      bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
                __func__, cusparseGetErrorString(status));

    csm->vecXValues = d_x;
  }

  if (d_y != csm->vecYValues) {
    if (csm->vecYValues != NULL)
      cusparseDestroyDnMat(csm->matY);

    status = cusparseCreateDnMat(&(csm->matY),
                                 matrix->n_rows,
                                 matrix->db_size,
                                 matrix->db_size, //n_rows, /* leading dimension */
                                 d_y,
                                 val_dtype,
                                 CUSPARSE_ORDER_ROW);

    if (CUSPARSE_STATUS_SUCCESS != status)
      bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
                __func__, cusparseGetErrorString(status));

    csm->vecYValues = d_y;
  }

  if (csm->dBuffer == NULL) {
    size_t bufferSize = 0;
    cs_real_t alpha = 1.0;
    cs_real_t beta = 1.0;  /* 0 should be enough for SmPV, 1 needed for
                              y = A.x + b.y
                              which is useful when y is initialized by
                              a separate diagonal da.x product */

    status = cusparseSpMM_bufferSize(_handle,
                                     CUSPARSE_OPERATION_NON_TRANSPOSE,
                                     CUSPARSE_OPERATION_NON_TRANSPOSE,
                                     &alpha,
                                     csm->matA,
                                     csm->matX,
                                     &beta,
                                     csm->matY,
                                     val_dtype,
                                     CUSPARSE_SPMM_ALG_DEFAULT,
                                     &bufferSize);

    if (CUSPARSE_STATUS_SUCCESS != status)
      bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
                __func__, cusparseGetErrorString(status));

    CS_CUDA_CHECK(cudaMalloc(&(csm->dBuffer), bufferSize));
  }

#endif
}

#endif // defined(HAVE_CUSPARSE)

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize CUDA matrix API.
 *
 * This frees resources such as the cuSPARSE handle, if used.
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_cuda_finalize(void)
{
  _stream = 0;

#if defined(HAVE_CUSPARSE)

  if (_handle != NULL) {
    cusparseDestroy(_handle);
    _handle = NULL;
  }

#endif
}

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
cs_matrix_spmv_cuda_set_stream(cudaStream_t  stream)
{
  _stream = stream;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return stream used for CUDA-based SpMV operations.
 *
 * This function is callable only from CUDA code.
 */
/*----------------------------------------------------------------------------*/

cudaStream_t
cs_matrix_spmv_cuda_get_stream(void)
{
  return _stream;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with CSR matrix, scalar CUDA version.
 *
 * \param[in]   matrix        pointer to matrix structure
 * \param[in]   exclude_diag  exclude diagonal if true,
 * \param[in]   sync          synchronize ghost cells if true
 * \param[in]   d_x           multipliying vector values (device)
 * \param[out]  d_y           resulting vector (device)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_cuda_csr(const cs_matrix_t  *matrix,
                        bool                exclude_diag,
                        bool                sync,
                        cs_real_t           d_x[restrict],
                        cs_real_t           d_y[restrict])
{
  const cs_matrix_struct_csr_t *ms
    = (const cs_matrix_struct_csr_t *)matrix->structure;
  const cs_matrix_coeff_csr_t *mc
    = (const cs_matrix_coeff_csr_t  *)matrix->coeffs;

  const cs_lnum_t *__restrict__ row_index
    = (const cs_lnum_t *)cs_get_device_ptr_const_pf
                           (const_cast<cs_lnum_t *>(ms->row_index));
  const cs_lnum_t *__restrict__ col_id
    = (const cs_lnum_t *)cs_get_device_ptr_const_pf
                           (const_cast<cs_lnum_t *>(ms->col_id));
  const cs_real_t *__restrict__ val
    = (const cs_real_t *)cs_get_device_ptr_const_pf
                           (const_cast<cs_real_t *>(mc->val));

  /* Ghost cell communication */

  if (sync) {
    cs_halo_state_t *hs = _pre_vector_multiply_sync_x_start(matrix, d_x);
    cs_halo_sync_wait(matrix->halo, d_x, hs);
  }

  /* Compute SpMV */

  unsigned int blocksize = 256;
  unsigned int gridsize
    = (unsigned int)ceil((double)ms->n_rows / blocksize);

  if (!exclude_diag)
    _mat_vect_p_l_csr<<<gridsize, blocksize, 0, _stream>>>
      (ms->n_rows, row_index, col_id, val, d_x, d_y);
  else
    _mat_vect_p_l_csr_exdiag<<<gridsize, blocksize, 0, _stream>>>
      (ms->n_rows, row_index, col_id, val, d_x, d_y);

  if (_stream == 0) {
    cudaStreamSynchronize(0);
    CS_CUDA_CHECK(cudaGetLastError());
  }
}

#if defined(HAVE_CUSPARSE)

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
                                 cs_real_t     d_y[restrict])
{
  cs_matrix_cusparse_map_t *csm
    = (cs_matrix_cusparse_map_t *)matrix->ext_lib_map;

  if (csm == NULL) {
    matrix->ext_lib_map = _set_cusparse_map(matrix);
    csm = (cs_matrix_cusparse_map_t *)matrix->ext_lib_map;
  }

  /* Ghost cell communication */

  if (sync) {
    cs_halo_state_t *hs = _pre_vector_multiply_sync_x_start(matrix,
                                                            (cs_real_t *)d_x);
    cs_halo_sync_wait(matrix->halo, (cs_real_t *)d_x, hs);
  }

  _update_cusparse_map(csm, matrix, d_x, d_y);

  cs_real_t alpha = 1.0;
  cs_real_t beta = 0.0;

  cusparseSetStream(_handle, _stream);

  cusparseStatus_t status = CUSPARSE_STATUS_SUCCESS;

#if defined(HAVE_CUSPARSE_GENERIC_API)

  cudaDataType_t val_dtype
    = (sizeof(cs_real_t) == 8) ? CUDA_R_64F : CUDA_R_32F;

  status = cusparseSpMV(_handle,
                        CUSPARSE_OPERATION_NON_TRANSPOSE,
                        &alpha,
                        csm->matA,
                        csm->vecX,
                        &beta,
                        csm->vecY,
                        val_dtype,
                        CUSPARSE_MV_ALG_DEFAULT,
                        csm->dBuffer);

  if (CUSPARSE_STATUS_SUCCESS != status)
    bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
              __func__, cusparseGetErrorString(status));

#else

  if (sizeof(cs_real_t) == 8) {
    double _alpha = alpha;
    double _beta = beta;

    status = cusparseDcsrmv(_handle,
                            CUSPARSE_OPERATION_NON_TRANSPOSE,
                            matrix->n_rows,
                            matrix->n_cols_ext,
                            csm->nnz,
                            &_alpha,
                            csm->descrA,
                            (const double *)csm->d_e_val,
                            (const int *)csm->d_row_index,
                            (const int *)csm->d_col_id,
                            (const double *)d_x,
                            &_beta,
                            (double *)d_y);
  }

  else if (sizeof(cs_real_t) == 4) {
    float _alpha = alpha;
    float _beta = beta;

    status = cusparseScsrmv(_handle,
                            CUSPARSE_OPERATION_NON_TRANSPOSE,
                            matrix->n_rows,
                            matrix->n_cols_ext,
                            csm->nnz,
                            &_alpha,
                            csm->descrA,
                            (const float *)csm->d_e_val,
                            (const int *)csm->d_row_index,
                            (const int *)csm->d_col_id,
                            (const float *)d_x,
                            &_beta,
                            (float *)d_y);

  if (CUSPARSE_STATUS_SUCCESS != status)
    bft_error(__FILE__, __LINE__, 0, _("%s: cuSPARSE error %d."),
              __func__, (int)status);

  }

#endif

  if (exclude_diag) {

    const cs_matrix_struct_csr_t *ms
      = (const cs_matrix_struct_csr_t *)matrix->structure;
    const cs_matrix_coeff_csr_t *mc
      = (const cs_matrix_coeff_csr_t  *)matrix->coeffs;
    const cs_lnum_t *__restrict__ d_row_index
      = (const cs_lnum_t *)cs_get_device_ptr_const_pf
                             (const_cast<cs_lnum_t *>(ms->row_index));
    const cs_lnum_t *__restrict__ d_col_id
      = (const cs_lnum_t *)cs_get_device_ptr_const_pf
                             (const_cast<cs_lnum_t *>(ms->col_id));
    const cs_real_t *__restrict__ d_val
      = (const cs_real_t *)cs_get_device_ptr_const_pf
                             (const_cast<cs_real_t *>(mc->val));

    unsigned int blocksize = 256;
    unsigned int gridsize
      = (unsigned int)ceil((double)ms->n_rows / blocksize);

    _mat_vect_p_l_csr_substract_diag<<<gridsize, blocksize, 0, _stream>>>
      (ms->n_rows, d_row_index, d_col_id, d_val,
       (const cs_real_t *)d_x, (cs_real_t *)d_y);

  }

  if (_stream == 0) {
    cudaStreamSynchronize(0);
    CS_CUDA_CHECK(cudaGetLastError());
  }
}

#endif /* defined(HAVE_CUSPARSE) */

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
                        cs_real_t           d_y[restrict])
{
  const cs_matrix_struct_dist_t *ms
    = (const cs_matrix_struct_dist_t *)matrix->structure;
  const cs_matrix_coeff_dist_t *mc
    = (const cs_matrix_coeff_dist_t *)matrix->coeffs;

  const cs_lnum_t *__restrict__ row_index
    = (const cs_lnum_t *)cs_get_device_ptr_const_pf
                           (const_cast<cs_lnum_t *>(ms->e.row_index));
  const cs_lnum_t *__restrict__ col_id
    = (const cs_lnum_t *)cs_get_device_ptr_const_pf
                          (const_cast<cs_lnum_t *>(ms->e.col_id));

  const cs_real_t *__restrict__ d_val
    = (const cs_real_t *)cs_get_device_ptr_const_pf
                           (const_cast<cs_real_t *>(mc->d_val));
  const cs_real_t *__restrict__ x_val
    = (const cs_real_t *)cs_get_device_ptr_const_pf
                           (const_cast<cs_real_t *>(mc->e_val));

  /* Ghost cell communication */

  if (sync) {
    cs_halo_state_t *hs = _pre_vector_multiply_sync_x_start(matrix, d_x);
    cs_halo_sync_wait(matrix->halo, d_x, hs);
  }

  /* Compute SpMV */

  unsigned int blocksize = 256;
  unsigned int gridsize
    = (unsigned int)ceil((double)ms->n_rows / blocksize);

  if (!exclude_diag)
    _mat_vect_p_l_msr<<<gridsize, blocksize, 0, _stream>>>
      (ms->n_rows, row_index, col_id, d_val, x_val, d_x, d_y);
  else
    _mat_vect_p_l_csr<<<gridsize, blocksize, 0, _stream>>>
      (ms->n_rows, row_index, col_id, x_val, d_x, d_y);

  if (_stream == 0) {
    cudaStreamSynchronize(0);
    CS_CUDA_CHECK(cudaGetLastError());
  }
}

#if defined(HAVE_CUSPARSE)

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
                                 cs_real_t     d_y[restrict])
{
  cs_matrix_cusparse_map_t *csm
    = (cs_matrix_cusparse_map_t *)matrix->ext_lib_map;

  if (csm == NULL) {
    matrix->ext_lib_map = _set_cusparse_map(matrix);
    csm = (cs_matrix_cusparse_map_t *)matrix->ext_lib_map;
  }

  /* Ghost cell communication */

  if (sync) {
    cs_halo_state_t *hs = _pre_vector_multiply_sync_x_start(matrix,
                                                            (cs_real_t *)d_x);
    cs_halo_sync_wait(matrix->halo, (cs_real_t *)d_x, hs);
  }

  _update_cusparse_map(csm, matrix, d_x, d_y);

  cs_real_t alpha = 1.;
  cs_real_t beta = 0.;

  if (!exclude_diag) {

    const cs_matrix_struct_dist_t *ms
      = (const cs_matrix_struct_dist_t *)matrix->structure;
    const cs_matrix_coeff_dist_t *mc
      = (const cs_matrix_coeff_dist_t *)matrix->coeffs;
    const cs_real_t *__restrict__ d_val
      = (const cs_real_t *)cs_get_device_ptr_const_pf
                             (const_cast<cs_real_t *>(mc->d_val));

    unsigned int blocksize = 256;
    unsigned int gridsize
      = (unsigned int)ceil((double)ms->n_rows / blocksize);

    _mat_vect_p_l_msr_diag<<<gridsize, blocksize, 0, _stream>>>
      (ms->n_rows, d_val, (const cs_real_t *)d_x, (cs_real_t *)d_y);

    beta = 1.;

  }

  cusparseSetStream(_handle, _stream);

  cusparseStatus_t status = CUSPARSE_STATUS_SUCCESS;

#if defined(HAVE_CUSPARSE_GENERIC_API)

  cudaDataType_t val_dtype
    = (sizeof(cs_real_t) == 8) ? CUDA_R_64F : CUDA_R_32F;

  status = cusparseSpMV(_handle,
                        CUSPARSE_OPERATION_NON_TRANSPOSE,
                        &alpha,
                        csm->matA,
                        csm->vecX,
                        &beta,
                        csm->vecY,
                        val_dtype,
                        CUSPARSE_MV_ALG_DEFAULT,
                        csm->dBuffer);

  if (CUSPARSE_STATUS_SUCCESS != status)
    bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
              __func__, cusparseGetErrorString(status));

#else

  if (sizeof(cs_real_t) == 8) {
    double _alpha = alpha;
    double _beta = beta;
    status = cusparseDcsrmv(_handle,
                            CUSPARSE_OPERATION_NON_TRANSPOSE,
                            matrix->n_rows,
                            matrix->n_cols_ext,
                            csm->nnz,
                            &_alpha,
                            csm->descrA,
                            (const double *)csm->d_e_val,
                            (const int *)csm->d_row_index,
                            (const int *)csm->d_col_id,
                            (const double *)d_x,
                            &_beta,
                            (double *)d_y);
  }

  else if (sizeof(cs_real_t) == 4) {
    float _alpha = alpha;
    float _beta = beta;

    status = cusparseScsrmv(_handle,
                            CUSPARSE_OPERATION_NON_TRANSPOSE,
                            matrix->n_rows,
                            matrix->n_cols_ext,
                            csm->nnz,
                            &_alpha,
                            csm->descrA,
                            (const float *)csm->d_e_val,
                            (const int *)csm->d_row_index,
                            (const int *)csm->d_col_id,
                            (const float *)d_x,
                            &_beta,
                            (float *)d_y);
  }

  if (CUSPARSE_STATUS_SUCCESS != status)
    bft_error(__FILE__, __LINE__, 0, _("%s: cuSPARSE error %d."),
              __func__, (int)status);

#endif

  if (_stream == 0) {
    cudaStreamSynchronize(0);
    CS_CUDA_CHECK(cudaGetLastError());
  }
}

#endif /* defined(HAVE_CUSPARSE) */

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
                          cs_real_t     d_y[restrict])
{
  const cs_matrix_struct_dist_t *ms
    = (const cs_matrix_struct_dist_t *)matrix->structure;
  const cs_matrix_coeff_dist_t *mc
    = (const cs_matrix_coeff_dist_t *)matrix->coeffs;

  const cs_lnum_t *__restrict__ row_index
    = (const cs_lnum_t *)cs_get_device_ptr_const_pf
                           (const_cast<cs_lnum_t *>(ms->e.row_index));
  const cs_lnum_t *__restrict__ col_id
    = (const cs_lnum_t *)cs_get_device_ptr_const_pf
                           (const_cast<cs_lnum_t *>(ms->e.col_id));

  const cs_real_t *__restrict__ d_val
    = (const cs_real_t *)cs_get_device_ptr_const_pf
                           (const_cast<cs_real_t *>(mc->d_val));
  const cs_real_t *__restrict__ x_val
    = (const cs_real_t *)cs_get_device_ptr_const_pf
                           (const_cast<cs_real_t *>(mc->e_val));

  /* Ghost cell communication */

  if (sync) {
    cs_halo_state_t *hs = _pre_vector_multiply_sync_x_start(matrix, d_x);
    cs_halo_sync_wait(matrix->halo, d_x, hs);
  }

  /* Compute SpMV */

  unsigned int blocksize = 128;
  unsigned int gridsize
    = (unsigned int)ceil((double)ms->n_rows / blocksize);

  if (!exclude_diag) {

    if (matrix->db_size == 3)
      _b_3_3_mat_vect_p_l_msr<<<gridsize, blocksize, 0, _stream>>>
        (ms->n_rows, col_id, row_index, d_val, x_val, d_x, d_y);
    else if (matrix->db_size == 6)
      _b_mat_vect_p_l_msr<6><<<gridsize, blocksize, 0, _stream>>>
        (ms->n_rows, col_id, row_index, d_val, x_val, d_x, d_y);
    else if (matrix->db_size == 9)
      _b_mat_vect_p_l_msr<9><<<gridsize, blocksize, 0, _stream>>>
        (ms->n_rows, col_id, row_index, d_val, x_val, d_x, d_y);
    else
      bft_error(__FILE__, __LINE__, 0, _("%s: block size %d not implemented."),
                __func__, (int)matrix->db_size);

  }
  else {

    if (matrix->db_size == 3)
      _b_3_3_mat_vect_p_l_msr_exdiag<<<gridsize, blocksize, 0, _stream>>>
        (ms->n_rows, col_id, row_index, d_val, x_val, d_x, d_y);
    else if (matrix->db_size == 6)
      _b_mat_vect_p_l_msr_exdiag<6><<<gridsize, blocksize, 0, _stream>>>
        (ms->n_rows, col_id, row_index, d_val, x_val, d_x, d_y);
    else if (matrix->db_size == 9)
      _b_mat_vect_p_l_msr_exdiag<9><<<gridsize, blocksize, 0, _stream>>>
        (ms->n_rows, col_id, row_index, d_val, x_val, d_x, d_y);
    else
      bft_error(__FILE__, __LINE__, 0, _("%s: block size %d not implemented."),
                __func__, (int)matrix->db_size);

  }

  if (_stream == 0) {
    CS_CUDA_CHECK(cudaStreamSynchronize(0));
    CS_CUDA_CHECK(cudaGetLastError());
  }
}

#if defined(HAVE_CUSPARSE_GENERIC_API)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with MSR matrix, block diagonal
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
cs_matrix_spmv_cuda_msr_b_cusparse(cs_matrix_t  *matrix,
                                   bool          exclude_diag,
                                   bool          sync,
                                   cs_real_t     d_x[restrict],
                                   cs_real_t     d_y[restrict])
{
  cs_matrix_cusparse_map_t *csm
    = (cs_matrix_cusparse_map_t *)matrix->ext_lib_map;

  if (csm == NULL) {
    matrix->ext_lib_map = _set_cusparse_map(matrix);
    csm = (cs_matrix_cusparse_map_t *)matrix->ext_lib_map;
  }

  /* Ghost cell communication */

  if (sync) {
    cs_halo_state_t *hs = _pre_vector_multiply_sync_x_start(matrix,
                                                            (cs_real_t *)d_x);
    cs_halo_sync_wait(matrix->halo, (cs_real_t *)d_x, hs);
  }

  _update_cusparse_map_block_diag(csm, matrix, d_x, d_y);

  cs_real_t alpha = 1.;
  cs_real_t beta = 0.;

  if (!exclude_diag) {

    const cs_matrix_struct_dist_t *ms
      = (const cs_matrix_struct_dist_t *)matrix->structure;
    const cs_matrix_coeff_dist_t *mc
      = (const cs_matrix_coeff_dist_t *)matrix->coeffs;
    const cs_real_t *__restrict__ d_val
      = (const cs_real_t *)cs_get_device_ptr_const_pf
                             (const_cast<cs_real_t *>(mc->d_val));

    unsigned int blocksize = 128;
    unsigned int gridsize
      = (unsigned int)ceil((double)ms->n_rows / blocksize);

    if (matrix->db_size == 3)
      _b_3_3_spmv_diag<<<gridsize, blocksize, 0, _stream>>>
        (ms->n_rows, d_val, (const cs_real_t *)d_x, (cs_real_t *)d_y);
    else if (matrix->db_size == 6)
      _b_spmv_diag<6><<<gridsize, blocksize, 0, _stream>>>
        (ms->n_rows, d_val, (const cs_real_t *)d_x, (cs_real_t *)d_y);
    else if (matrix->db_size == 9)
      _b_spmv_diag<9><<<gridsize, blocksize, 0, _stream>>>
        (ms->n_rows, d_val, (const cs_real_t *)d_x, (cs_real_t *)d_y);
    else
      bft_error(__FILE__, __LINE__, 0, _("%s: block size %d not implemented."),
                __func__, (int)matrix->db_size);

    beta = 1.;

  }

  cusparseSetStream(_handle, _stream);

  cusparseStatus_t status = CUSPARSE_STATUS_SUCCESS;

  cudaDataType_t val_dtype
    = (sizeof(cs_real_t) == 8) ? CUDA_R_64F : CUDA_R_32F;

  status = cusparseSpMM(_handle,
                        CUSPARSE_OPERATION_NON_TRANSPOSE,
                        CUSPARSE_OPERATION_NON_TRANSPOSE,
                        &alpha,
                        csm->matA,
                        csm->matX,
                        &beta,
                        csm->matY,
                        val_dtype,
                        CUSPARSE_SPMM_ALG_DEFAULT,
                        csm->dBuffer);

  if (CUSPARSE_STATUS_SUCCESS != status)
    bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
              __func__, cusparseGetErrorString(status));

  if (_stream == 0) {
    cudaStreamSynchronize(0);
    CS_CUDA_CHECK(cudaGetLastError());
  }
}

#endif /* defined(HAVE_CUSPARSE_GENERIC_API) */

#if defined(HAVE_CUSPARSE)

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
                                    cs_real_t     d_y[restrict])
{
  cs_matrix_cusparse_map_t *csm
    = (cs_matrix_cusparse_map_t *)matrix->ext_lib_map;

  if (csm == NULL) {
    matrix->ext_lib_map = _set_cusparse_map(matrix);
    csm = (cs_matrix_cusparse_map_t *)matrix->ext_lib_map;
  }

  /* Ghost cell communication */

  if (sync) {
    cs_halo_state_t *hs = _pre_vector_multiply_sync_x_start(matrix,
                                                            (cs_real_t *)d_x);
    cs_halo_sync_wait(matrix->halo, (cs_real_t *)d_x, hs);
  }

  /* no update_cusparse_map type function call here as only
     the non-generic API is available here */

  cs_real_t alpha = 1.;
  cs_real_t beta = 0.;

  if (!exclude_diag) {

    const cs_matrix_struct_dist_t *ms
      = (const cs_matrix_struct_dist_t *)matrix->structure;
    const cs_matrix_coeff_dist_t *mc
      = (const cs_matrix_coeff_dist_t *)matrix->coeffs;
    const cs_real_t *__restrict__ d_val
      = (const cs_real_t *)cs_get_device_ptr_const_pf
                             (const_cast<cs_real_t *>(mc->d_val));

    unsigned int blocksize = 128;
    unsigned int gridsize
      = (unsigned int)ceil((double)ms->n_rows / blocksize);

    if (matrix->db_size == 3)
      _b_3_3_spmv_diag<<<gridsize, blocksize, 0, _stream>>>
        (ms->n_rows, d_val, (const cs_real_t *)d_x, (cs_real_t *)d_y);
    else if (matrix->db_size == 6)
      _b_spmv_diag<6><<<gridsize, blocksize, 0, _stream>>>
        (ms->n_rows, d_val, (const cs_real_t *)d_x, (cs_real_t *)d_y);
    else if (matrix->db_size == 9)
      _b_spmv_diag<9><<<gridsize, blocksize, 0, _stream>>>
        (ms->n_rows, d_val, (const cs_real_t *)d_x, (cs_real_t *)d_y);
    else
      bft_error(__FILE__, __LINE__, 0, _("%s: block size %d not implemented."),
                __func__, (int)matrix->db_size);

    beta = 1.;

  }

  cusparseSetStream(_handle, _stream);

  cusparseStatus_t status = CUSPARSE_STATUS_SUCCESS;

  if (sizeof(cs_real_t) == 8) {
    double _alpha = alpha;
    double _beta = beta;
    status = cusparseDbsrmv(_handle,
                            CUSPARSE_DIRECTION_ROW,
                            CUSPARSE_OPERATION_NON_TRANSPOSE,
                            matrix->n_rows,
                            matrix->n_cols_ext,
                            csm->nnz,
                            &_alpha,
                            csm->descrA,
                            (const double *)csm->d_e_val,
                            (const int *)csm->d_row_index,
                            (const int *)csm->d_col_id,
                            matrix->eb_size,
                            (const double *)d_x,
                            &_beta,
                            (double *)d_y);
  }

  else if (sizeof(cs_real_t) == 4) {
    float _alpha = alpha;
    float _beta = beta;

    status = cusparseSbsrmv(_handle,
                            CUSPARSE_DIRECTION_ROW,
                            CUSPARSE_OPERATION_NON_TRANSPOSE,
                            matrix->n_rows,
                            matrix->n_cols_ext,
                            csm->nnz,
                            &_alpha,
                            csm->descrA,
                            (const float *)csm->d_e_val,
                            (const int *)csm->d_row_index,
                            (const int *)csm->d_col_id,
                            matrix->eb_size,
                            (const float *)d_x,
                            &_beta,
                            (float *)d_y);
  }

  if (CUSPARSE_STATUS_SUCCESS != status)
    bft_error(__FILE__, __LINE__, 0, _("%s: cuSPARSE error %d."),
              __func__, (int)status);

  if (_stream == 0) {
    cudaStreamSynchronize(0);
    CS_CUDA_CHECK(cudaGetLastError());
  }
}

#endif /* defined(HAVE_CUSPARSE) */

/*----------------------------------------------------------------------------*/

END_C_DECLS
