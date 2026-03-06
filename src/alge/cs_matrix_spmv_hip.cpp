/*============================================================================
 * Sparse Matrix-vector multiplication kernels using HIP.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#if defined(HAVE_HIPSPARSE)
#include <hipsparse.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_base.h"
#include "base/cs_base_hip.h"
#include "base/cs_halo.h"
#include "base/cs_halo_perio.h"
#include "base/cs_log.h"
#include "base/cs_mem.h"
#include "base/cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "alge/cs_matrix.h"
#include "alge/cs_matrix_priv.h"
#include "alge/cs_matrix_spmv.h"

/*----------------------------------------------------------------------------*/
/*! \file
 *
 * \brief Sparse Matrix SpMV operations with HIP.
 */
/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

#if defined(HAVE_HIPSPARSE)

/* Mapping of matrix coefficients and structure to cuSPARSE */
/*----------------------------------------------------------*/

typedef struct _cs_matrix_cusparse_map_t {

  bool  block_diag;             /* Use identity blocks diagonal structure ? */

  hipsparseSpMatDescr_t  matA;   /* Handle to cusparse Matrix */

  hipsparseDnMatDescr_t  matX;   /* Handle to cusparse Matrix (blocked vector) */
  hipsparseDnMatDescr_t  matY;   /* Handle to cusparse Matrix (blocked vector) */

  hipsparseDnVecDescr_t  vecX;   /* Handle to cusparse Vector */
  hipsparseDnVecDescr_t  vecY;   /* Handle to cusparse output Vector */

  void  *vecXValues;            /* Pointer to vector values */
  void  *vecYValues;            /* Pointer to vector values */

  void  *dBuffer;               /* Associated buffer */

  /* When not using generic API */

  int  nnz;                     /* Number of nonzeroes */
  hipsparseMatDescr_t *descrA;  /* Handle to cusparse Matrix description */

  void  *d_row_index;         /* Pointer to row index */
  void  *d_col_id;            /* Pointer to column ids */
  void  *d_e_val;             /* Pointer to matrix extradiagonal values */

} cs_matrix_cusparse_map_t;

#endif // defined(HAVE_HIPSPARSE)

/*============================================================================
 *  Global variables
 *============================================================================*/

static hipStream_t _stream = 0;

#if defined(HAVE_HIPSPARSE)

static hipsparseHandle_t  _handle = nullptr;

#endif

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/* \brief Zero range of elements.
 *
 * \param[in]   n   number of elements
 * \param[out]  x   resulting vector
 */
/*----------------------------------------------------------------------------*/

__global__ static void
_zero_range(cs_lnum_t    n,
            cs_real_t   *__restrict__ x)
{
  cs_lnum_t ii = blockIdx.x * blockDim.x + threadIdx.x;

  if (ii < n)
    x[ii] = 0;
}

/*----------------------------------------------------------------------------*/
/* \brief Local diagonal contribution y = Da.x  + y.
 *
 * \param[in]   n_rows      number of local rows
 * \param[in]   n_cols_ext  number of local columns (with ghosts)
 * \param[in]   d_val       pointer to diagonal matrix values
 * \param[in]   x           multipliying vector values
 * \param[out]  y           resulting vector
 */
/*----------------------------------------------------------------------------*/

__global__ static void
_mat_vect_p_l_native_diag(cs_lnum_t         n_rows,
                          cs_lnum_t         n_cols_ext,
                          const cs_real_t  *__restrict__ d_val,
                          const cs_real_t  *__restrict__ x,
                          cs_real_t        *__restrict__ y)
{
  cs_lnum_t ii = blockIdx.x * blockDim.x + threadIdx.x;

  if (ii < n_rows)
    y[ii] = d_val[ii] * x[ii];

  else if (ii < n_cols_ext)
    y[ii] = 0;
}

/*----------------------------------------------------------------------------
 * SpMV extradiagonal terms using native to face-based array and scatter
 * approach, handling conflicts through atomic add.
 *
 * Non-symmetric matrix case.
 *
 * parameters:
 *   n_edges  <-- local number of internal graph edges (mesh faces)
 *   edges    <-- edges (mesh face -> cells) connectivity
 *   xa       <-- extradiagonal values
 *   x        <-- vector
 *   y        <-> vector
 *----------------------------------------------------------------------------*/

__global__ static void
_mat_vect_p_l_native_exdiag(cs_lnum_t           n_edges,
                            const cs_lnum_2_t  * __restrict__ edges,
                            const cs_real_t    *__restrict__ xa,
                            const cs_real_t    *__restrict__ x,
                            cs_real_t          *__restrict__ y)
{
  cs_lnum_t edge_id = blockIdx.x * blockDim.x + threadIdx.x;

  if (edge_id < n_edges) {
    cs_lnum_t ii = edges[edge_id][0];
    cs_lnum_t jj = edges[edge_id][1];
    cs_real_t x_ii = __ldg(x + ii);
    cs_real_t x_jj = __ldg(x + jj);
#if 0
    atomicAdd(&y[ii], xa[edge_id*2]     * x_jj);
    atomicAdd(&y[jj], xa[edge_id*2 + 1] * x_ii);
#else
    cs_real_t aii_x = xa[edge_id*2]     * x_jj;
    cs_real_t ajj_x = xa[edge_id*2+1]   * x_ii;

    using sum_v = assembled_value<cs_real_t>;
    sum_v vii, vjj;

    vii.get() = aii_x;
    vjj.get() = ajj_x;

    sum_v::ref(y[ii]).conflict_free_add(-1u, vii);
    sum_v::ref(y[jj]).conflict_free_add(-1u, vjj);
#endif
  }
}

/*----------------------------------------------------------------------------
 * SpMV extradiagonal terms using native to face-based array and scatter
 * approach, handling conflicts through atomic add.
 *
 * Symmetric matrix case.
 *
 * parameters:
 *   n_edges  <-- local number of internal graph edges (mesh faces)
 *   edges    <-- edges (mesh face -> cells) connectivity
 *   xa       <-- extradiagonal values
 *   x        <-- vector
 *   y        <-> vector
 *----------------------------------------------------------------------------*/

__global__ static void
_mat_vect_p_l_native_exdiag_sym(cs_lnum_t           n_edges,
                                const cs_lnum_2_t  * __restrict__ edges,
                                const cs_real_t    *__restrict__ xa,
                                const cs_real_t    *__restrict__ x,
                                cs_real_t          *__restrict__ y)
{
  cs_lnum_t edge_id = blockIdx.x * blockDim.x + threadIdx.x;

  if (edge_id < n_edges) {
    cs_lnum_t ii = edges[edge_id][0];
    cs_lnum_t jj = edges[edge_id][1];
    cs_real_t x_ii = __ldg(x + ii);
    cs_real_t x_jj = __ldg(x + jj);
#if 0
    atomicAdd(&y[ii], xa[edge_id] * x_jj);
    atomicAdd(&y[jj], xa[edge_id] * x_ii);
#else
    cs_real_t aii_x = xa[edge_id] * x_jj;
    cs_real_t ajj_x = xa[edge_id] * x_ii;

    using sum_v = assembled_value<cs_real_t>;
    sum_v vii, vjj;

    vii.get() = aii_x;
    vjj.get() = ajj_x;

    sum_v::ref(y[ii]).conflict_free_add(-1u, vii);
    sum_v::ref(y[jj]).conflict_free_add(-1u, vjj);
#endif
  }
}

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
  cs_lnum_t r_ii = blockIdx.x * blockDim.x + threadIdx.x;
  if (r_ii < n_rows) {
    const cs_lnum_t ii = r_ii/3;
    const cs_lnum_t kk = r_ii%3;

    y[r_ii] =   d_val[ii * 9 + kk * 3]     * x[ii * 3]
              + d_val[ii * 9 + kk * 3 + 1] * x[ii * 3 + 1]
              + d_val[ii * 9 + kk * 3 + 2] * x[ii * 3 + 2];
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
_b_spmv_diag(cs_lnum_t         n_rows,
             const cs_real_t  *__restrict__ d_val,
             const cs_real_t  *__restrict__ x,
             cs_real_t        *__restrict__ y)
{
  cs_lnum_t r_ii = blockIdx.x * blockDim.x + threadIdx.x;
  if (r_ii < n_rows) {
    const cs_lnum_t ii = r_ii/n;
    const cs_lnum_t kk = r_ii%n;
    const cs_lnum_t nn = n*n;

    cs_real_t sii = 0.;
    for (cs_lnum_t ll = 0; ll < n; ll++) {
      sii += d_val[ii*nn + kk*n + ll] * x[ii*n + ll];
    }

    y[r_ii] = sii;
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
  cs_lnum_t r_ii = blockIdx.x * blockDim.x + threadIdx.x;
  if (r_ii < n_rows) {
    const cs_lnum_t ii = r_ii/3;
    const cs_lnum_t kk = r_ii%3;

    const cs_lnum_t *__restrict__ _col_id = col_id + row_index[ii];
    const cs_real_t *__restrict__ m_row  = x_val + row_index[ii];
    cs_lnum_t n_cols = row_index[ii + 1] - row_index[ii];
    cs_real_t sii =   d_val[ii * 9 + kk * 3]     * x[ii * 3]
                    + d_val[ii * 9 + kk * 3 + 1] * x[ii * 3 + 1]
                    + d_val[ii * 9 + kk * 3 + 2] * x[ii * 3 + 2];

    for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
      sii += m_row[jj] * __ldg(x + (_col_id[jj]*3 + kk));
    }

    y[r_ii] = sii;
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
  cs_lnum_t r_ii = blockIdx.x * blockDim.x + threadIdx.x;
  if (r_ii < n_rows) {
    const cs_lnum_t ii = r_ii/3;
    const cs_lnum_t kk = r_ii%3;

    const cs_lnum_t *__restrict__ _col_id = col_id + row_index[ii];
    const cs_real_t *__restrict__ m_row  = x_val + row_index[ii];
    cs_lnum_t n_cols = row_index[ii + 1] - row_index[ii];
    cs_real_t sii = 0.;

    for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
      sii += m_row[jj] * __ldg(x + (_col_id[jj]*3 + kk));
    }

    y[r_ii] = sii;
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
  cs_lnum_t r_ii = blockIdx.x * blockDim.x + threadIdx.x;
  if (r_ii < n_rows) {
    const cs_lnum_t ii = r_ii/n;
    const cs_lnum_t kk = r_ii%n;
    const cs_lnum_t nn = n*n;

    const cs_lnum_t *__restrict__ _col_id = col_id + row_index[ii];
    const cs_real_t *__restrict__ m_row  = x_val + row_index[ii];
    cs_lnum_t n_cols = row_index[ii + 1] - row_index[ii];
    cs_real_t sii = 0;

    for (cs_lnum_t ll = 0; ll < n; ll++) {
      sii += d_val[ii*nn + kk*n + ll] * x[ii*n + ll];
    }

    for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
      sii += m_row[jj] * __ldg(x + (_col_id[jj]*n + kk));
    }

    y[r_ii] = sii;
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
  cs_lnum_t r_ii = blockIdx.x * blockDim.x + threadIdx.x;
  if (r_ii < n_rows) {
    const cs_lnum_t ii = r_ii/n;
    const cs_lnum_t kk = r_ii%n;

    const cs_lnum_t *__restrict__ _col_id = col_id + row_index[ii];
    const cs_real_t *__restrict__ m_row  = x_val + row_index[ii];
    cs_lnum_t n_cols = row_index[ii + 1] - row_index[ii];
    cs_real_t sii = 0.;

    for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
      sii += m_row[jj] * __ldg(x + (_col_id[jj]*n + kk));
    }

    y[r_ii] = sii;
  }
}

/*----------------------------------------------------------------------------
 * Start synchronization of ghost values prior to matrix.vector product.
 *
 * Values are packed on the device, so:
 * - If MPI is HIP-aware, no values need to go through the host
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
                                  cs_real_t            d_x[])
{
  cs_halo_state_t *hs = nullptr;

  if (matrix->halo != nullptr) {

    if (_stream != 0) {
      CS_HIP_CHECK(hipStreamSynchronize(_stream));
      CS_HIP_CHECK(hipGetLastError());
    }

    hs = cs_halo_state_get_default();

    cs_halo_sync_pack_d(matrix->halo,
                        CS_HALO_STANDARD,
                        CS_REAL_TYPE,
                        matrix->db_size,
                        d_x,
                        nullptr,
                        hs);

    cs_halo_sync_start(matrix->halo, d_x, hs);

  }

  return hs;
}

/*----------------------------------------------------------------------------
 * Synchronize ghost values prior to matrix.vector product
 *
 * parameters:
 *   matrix        <-- pointer to matrix structure
 *   x             <-> multipliying vector values (ghost values updated)
 *----------------------------------------------------------------------------*/

static void
_pre_vector_multiply_sync_x_end(const cs_matrix_t   *matrix,
                                cs_halo_state_t     *hs,
                                cs_real_t           *restrict x)
{
  if (hs != nullptr) {

    assert(matrix->halo != nullptr);

    cs_halo_sync_wait(matrix->halo, x, hs);

    /* Synchronize periodic values */

#if !defined(_CS_UNIT_MATRIX_TEST) /* unit tests do not link with full library */

    // FIXME: ensure this is done on the GPU.

    if (matrix->halo->n_transforms > 0) {
      if (matrix->db_size == 3)
        cs_halo_perio_sync_var_vect(matrix->halo,
                                    CS_HALO_STANDARD,
                                    x,
                                    matrix->db_size);
      else if (matrix->db_size == 6)
        cs_halo_perio_sync_var_sym_tens(matrix->halo,
                                        CS_HALO_STANDARD,
                                        x);
    }

#endif
  }
}

#if defined(HAVE_HIPSPARSE)

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

  if (csm == nullptr)
    return;

  hipsparseDestroySpMat(csm->matA);

  if (csm->block_diag == false) {
    if (csm->vecXValues != nullptr)
      hipsparseDestroyDnVec(csm->vecX);
    if (csm->vecYValues != nullptr)
      hipsparseDestroyDnVec(csm->vecY);
  }
  else {
    if (csm->vecXValues != nullptr)
      hipsparseDestroyDnMat(csm->matX);
    if (csm->vecYValues != nullptr)
      hipsparseDestroyDnMat(csm->matY);
  }

  if (csm->dBuffer != nullptr) {
    CS_HIP_CHECK(hipFree(csm->dBuffer));
    csm->dBuffer = nullptr;
  }

  csm->block_diag = false;

  csm->vecXValues = nullptr;
  csm->vecYValues = nullptr;

  csm->nnz = 0;
  csm->d_row_index = nullptr;
  csm->d_col_id = nullptr;
  csm->d_e_val = nullptr;

  CS_FREE(matrix->ext_lib_map);
  matrix->destroy_adaptor = nullptr;
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

  if (csm != nullptr) {
    _unset_cusparse_map(matrix);
  }
  else {
    CS_MALLOC(csm, 1, cs_matrix_cusparse_map_t);
    matrix->ext_lib_map = (void *)csm;
  }
  matrix->destroy_adaptor = _unset_cusparse_map;

  const void *row_index, *col_id;
  const void *e_val;
  cs_lnum_t nnz = 0;

  if (matrix->type == CS_MATRIX_CSR) {
    const cs_matrix_struct_csr_t *ms
      = (const cs_matrix_struct_csr_t  *)matrix->structure;
    const cs_matrix_coeff_t *mc
      = (const cs_matrix_coeff_t *)matrix->coeffs;
    nnz = ms->row_index[matrix->n_rows];
    row_index = cs_get_device_ptr_const
                  (const_cast<cs_lnum_t *>(ms->row_index));
    col_id = cs_get_device_ptr_const
               (const_cast<cs_lnum_t *>(ms->col_id));
    e_val = cs_get_device_ptr_const
              (const_cast<cs_real_t *>(mc->val));
  }
  else {
    const cs_matrix_struct_dist_t *ms
      = (const cs_matrix_struct_dist_t *)matrix->structure;
    const cs_matrix_coeff_t *mc
      = (const cs_matrix_coeff_t *)matrix->coeffs;
    nnz = ms->e.row_index[matrix->n_rows];
    row_index = cs_get_device_ptr_const
                  (const_cast<cs_lnum_t *>(ms->e.row_index));
    col_id = cs_get_device_ptr_const
               (const_cast<cs_lnum_t *>(ms->e.col_id));
    e_val = cs_get_device_ptr_const
              (const_cast<cs_real_t *>(mc->e_val));
  }

  hipsparseStatus_t status = HIPSPARSE_STATUS_SUCCESS;

  if (_handle == nullptr)
    status = hipsparseCreate(&_handle);

  if (matrix->db_size > matrix->eb_size)
    csm->block_diag = true;
  else
    csm->block_diag = false;

  if (HIPSPARSE_STATUS_SUCCESS != status)
    bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
              __func__, hipsparseGetErrorString(status));

  csm->vecXValues = nullptr;  /* Pointer to vector values */
  csm->vecYValues = nullptr;  /* Pointer to vector values */
  csm->dBuffer = nullptr;

  csm->d_e_val = nullptr;
  csm->d_row_index = nullptr;
  csm->d_col_id = nullptr;

  if (matrix->eb_size == 1) {

    hipsparseIndexType_t index_dtype
      = (sizeof(cs_lnum_t) == 4) ? HIPSPARSE_INDEX_32I : HIPSPARSE_INDEX_64I;
    hipDataType val_dtype
      = (sizeof(cs_real_t) == 8) ? HIP_R_64F : HIP_R_32F;

    status = hipsparseCreateCsr(&(csm->matA),
                               matrix->n_rows,
                               matrix->n_cols_ext,
                               nnz,
                               const_cast<void *>(row_index),
                               const_cast<void *>(col_id),
                               const_cast<void *>(e_val),
                               index_dtype,
                               index_dtype,
                               HIPSPARSE_INDEX_BASE_ZERO,
                               val_dtype);

  }
  else {

    csm->nnz = nnz;
    csm->d_e_val = const_cast<void *>(e_val);

    csm->d_row_index = const_cast<void *>(row_index);
    csm->d_col_id = const_cast<void *>(col_id);

    status = hipsparseCreateMatDescr(&(csm->descrA));

  }

  if (HIPSPARSE_STATUS_SUCCESS != status)
    bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
              __func__, hipsparseGetErrorString(status));

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
  assert(csm != nullptr);

  hipsparseSpMVAlg_t spmv_alg_type = HIPSPARSE_SPMV_ALG_DEFAULT;

  hipsparseStatus_t status = HIPSPARSE_STATUS_SUCCESS;
  hipDataType val_dtype
    = (sizeof(cs_real_t) == 8) ? HIP_R_64F : HIP_R_32F;

  if (d_x != csm->vecXValues) {
    if (csm->vecXValues != nullptr)
      hipsparseDestroyDnVec(csm->vecX);

    status = hipsparseCreateDnVec(&(csm->vecX),
                                 matrix->n_cols_ext,
                                 d_x,
                                 val_dtype);

    if (HIPSPARSE_STATUS_SUCCESS != status)
      bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
                __func__, hipsparseGetErrorString(status));

    csm->vecXValues = d_x;
  }

  if (d_y != csm->vecYValues) {
    if (csm->vecYValues != nullptr)
      hipsparseDestroyDnVec(csm->vecY);

    status = hipsparseCreateDnVec(&(csm->vecY),
                                 matrix->n_rows,
                                 d_y,
                                 val_dtype);

    if (HIPSPARSE_STATUS_SUCCESS != status)
      bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
                __func__, hipsparseGetErrorString(status));

    csm->vecYValues = d_y;
  }

  if (csm->dBuffer == nullptr) {
    size_t bufferSize = 0;
    cs_real_t alpha = 1.0;
    cs_real_t beta = 1.0;  /* 0 should be enough for SmPV, 1 needed for
                              y = A.x + b.y
                              which is useful when y is initialized by
                              a separate diagonal da.x product */

    status = hipsparseSpMV_bufferSize(_handle,
                                     HIPSPARSE_OPERATION_NON_TRANSPOSE,
                                     &alpha,
                                     csm->matA,
                                     csm->vecX,
                                     &beta,
                                     csm->vecY,
                                     val_dtype,
                                     spmv_alg_type,
                                     &bufferSize);

    CS_HIP_CHECK(hipMalloc(&(csm->dBuffer), bufferSize));
  }

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
  assert(csm != nullptr);

  hipsparseStatus_t status = HIPSPARSE_STATUS_SUCCESS;
  hipDataType val_dtype
    = (sizeof(cs_real_t) == 8) ? HIP_R_64F : HIP_R_32F;

  if (d_x != csm->vecXValues) {
    if (csm->vecXValues != nullptr)
      hipsparseDestroyDnMat(csm->matX);

    status = hipsparseCreateDnMat(&(csm->matX),
                                 matrix->n_cols_ext,
                                 matrix->db_size,
                                 matrix->db_size,
                                 d_x,
                                 val_dtype,
                                 HIPSPARSE_ORDER_ROW);

    if (HIPSPARSE_STATUS_SUCCESS != status)
      bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
                __func__, hipsparseGetErrorString(status));

    csm->vecXValues = d_x;
  }

  if (d_y != csm->vecYValues) {
    if (csm->vecYValues != nullptr)
      hipsparseDestroyDnMat(csm->matY);

    status = hipsparseCreateDnMat(&(csm->matY),
                                 matrix->n_rows,
                                 matrix->db_size,
                                 matrix->db_size,
                                 d_y,
                                 val_dtype,
                                 HIPSPARSE_ORDER_ROW);

    if (HIPSPARSE_STATUS_SUCCESS != status)
      bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
                __func__, hipsparseGetErrorString(status));

    csm->vecYValues = d_y;
  }

  if (csm->dBuffer == nullptr) {
    size_t bufferSize = 0;
    cs_real_t alpha = 1.0;
    cs_real_t beta = 1.0;  /* 0 should be enough for SmPV, 1 needed for
                              y = A.x + b.y
                              which is useful when y is initialized by
                              a separate diagonal da.x product */

    status = hipsparseSpMM_bufferSize(_handle,
                                     HIPSPARSE_OPERATION_NON_TRANSPOSE,
                                     HIPSPARSE_OPERATION_NON_TRANSPOSE,
                                     &alpha,
                                     csm->matA,
                                     csm->matX,
                                     &beta,
                                     csm->matY,
                                     val_dtype,
                                     HIPSPARSE_SPMM_ALG_DEFAULT,
                                     &bufferSize);

    if (HIPSPARSE_STATUS_SUCCESS != status)
      bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
                __func__, hipsparseGetErrorString(status));

    CS_HIP_CHECK(hipMalloc(&(csm->dBuffer), bufferSize));
  }

}

#endif // defined(HAVE_HIPSPARSE)

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize HIP matrix API.
 *
 * This frees resources such as the cuSPARSE handle, if used.
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_hip_finalize(void)
{
  _stream = 0;

#if defined(HAVE_HIPSPARSE)

  if (_handle != nullptr) {
    hipsparseDestroy(_handle);
    _handle = nullptr;
  }

#endif
}

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
cs_matrix_spmv_hip_set_stream(hipStream_t  stream)
{
  _stream = stream;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return stream used for HIP-based SpMV operations.
 *
 * This function is callable only from HIP code.
 */
/*----------------------------------------------------------------------------*/

hipStream_t
cs_matrix_spmv_hip_get_stream(void)
{
  return _stream;
}

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
cs_matrix_spmv_hip_native(const cs_matrix_t  *matrix,
                           bool                exclude_diag,
                           bool                sync,
                           cs_real_t           d_x[],
                           cs_real_t           d_y[])
{
  const cs_matrix_struct_native_t  *ms
    = (const cs_matrix_struct_native_t *)matrix->structure;

  const cs_matrix_coeff_t  *mc
    = (const cs_matrix_coeff_t *)matrix->coeffs;

  const cs_real_t *__restrict__ da
    = (const cs_real_t *)cs_get_device_ptr_const
                           (const_cast<cs_real_t *>(mc->d_val));
  const cs_real_t *__restrict__ xa
    = (const cs_real_t *)cs_get_device_ptr_const
                           (const_cast<cs_real_t *>(mc->e_val));

  /* Ghost cell communication */

  cs_halo_state_t *hs = nullptr;
  if (sync)
    hs = _pre_vector_multiply_sync_x_start(matrix, d_x);

  /* Diagonal part of matrix.vector product */

  unsigned int blocksize = 256;
  unsigned int gridsize = cs_hip_grid_size(ms->n_cols_ext, blocksize);

  if (!exclude_diag)
    _mat_vect_p_l_native_diag<<<gridsize, blocksize, 0, _stream>>>
      (ms->n_rows, ms->n_cols_ext, da, d_x, d_y);
  else
    _zero_range<<<gridsize, blocksize, 0, _stream>>>
      (ms->n_cols_ext, d_y);

  /* Finalize ghost cell comunication if overlap used */

  if (hs != nullptr)
    cs_halo_sync_wait(matrix->halo, d_x, hs);

  CS_HIP_CHECK(hipStreamSynchronize(_stream));
  CS_HIP_CHECK(hipGetLastError());

  /* Non-diagonal terms */

  if (xa != nullptr) {
    gridsize = cs_hip_grid_size(ms->n_edges, blocksize);

    const cs_lnum_2_t *restrict edges
      = (const cs_lnum_2_t *)cs_get_device_ptr_const
                               (const_cast<cs_lnum_2_t *>(ms->edges));

#if 1
    if (mc->symmetric)
      _mat_vect_p_l_native_exdiag_sym<<<gridsize, blocksize, 0, _stream>>>
        (ms->n_edges, edges, xa, d_x, d_y);
    else
      _mat_vect_p_l_native_exdiag<<<gridsize, blocksize, 0, _stream>>>
        (ms->n_edges, edges, xa, d_x, d_y);

#else
    if (mc->symmetric) {
      for (cs_lnum_t e_id = 0; e_id < ms->n_edges; e_id++) {
        cs_lnum_t ii = edges[e_id][0];
        cs_lnum_t jj = edges[e_id][1];
        d_y[ii] += xa[e_id] * d_x[jj];
        d_y[jj] += xa[e_id] * d_x[ii];
      }

    }
    else {
      for (cs_lnum_t e_id = 0; e_id < ms->n_edges; e_id++) {
        cs_lnum_t ii = edges[e_id][0];
        cs_lnum_t jj = edges[e_id][1];
        d_y[ii] += xa[2*e_id] * d_x[jj];
        d_y[jj] += xa[2*e_id + 1] * d_x[ii];
      }
    }
#endif
  }

  CS_HIP_CHECK(hipStreamSynchronize(_stream));
  CS_HIP_CHECK(hipGetLastError());

  if (_stream == 0) {
    CS_HIP_CHECK(hipStreamSynchronize(_stream));
    CS_HIP_CHECK(hipGetLastError());
  }
}

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
                        cs_real_t     d_y[])
{
  const cs_matrix_struct_csr_t *ms
    = (const cs_matrix_struct_csr_t *)matrix->structure;
  const cs_matrix_coeff_t *mc
    = (const cs_matrix_coeff_t  *)matrix->coeffs;

  const cs_lnum_t *__restrict__ row_index
    = (const cs_lnum_t *)cs_get_device_ptr_const
                           (const_cast<cs_lnum_t *>(ms->row_index));
  const cs_lnum_t *__restrict__ col_id
    = (const cs_lnum_t *)cs_get_device_ptr_const
                           (const_cast<cs_lnum_t *>(ms->col_id));
  const cs_real_t *__restrict__ val
    = (const cs_real_t *)cs_get_device_ptr_const
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
    CS_HIP_CHECK(hipStreamSynchronize(_stream));
    CS_HIP_CHECK(hipGetLastError());
  }
}

#if defined(HAVE_HIPSPARSE)

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
cs_matrix_spmv_hip_csr_cusparse(cs_matrix_t  *matrix,
                                 bool          exclude_diag,
                                 bool          sync,
                                 cs_real_t     d_x[],
                                 cs_real_t     d_y[])
{
  cs_matrix_cusparse_map_t *csm
    = (cs_matrix_cusparse_map_t *)matrix->ext_lib_map;

  if (csm == nullptr) {
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

  hipsparseSetStream(_handle, _stream);

  hipsparseStatus_t status = HIPSPARSE_STATUS_SUCCESS;

  hipsparseSpMVAlg_t spmv_alg_type = HIPSPARSE_SPMV_ALG_DEFAULT;

  hipDataType val_dtype
    = (sizeof(cs_real_t) == 8) ? HIP_R_64F : HIP_R_32F;

  status = hipsparseSpMV(_handle,
                        HIPSPARSE_OPERATION_NON_TRANSPOSE,
                        &alpha,
                        csm->matA,
                        csm->vecX,
                        &beta,
                        csm->vecY,
                        val_dtype,
                        spmv_alg_type,
                        csm->dBuffer);

  if (HIPSPARSE_STATUS_SUCCESS != status)
    bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
              __func__, hipsparseGetErrorString(status));

  if (exclude_diag) {

    const cs_matrix_struct_csr_t *ms
      = (const cs_matrix_struct_csr_t *)matrix->structure;
    const cs_matrix_coeff_t *mc
      = (const cs_matrix_coeff_t *)matrix->coeffs;
    const cs_lnum_t *__restrict__ d_row_index
      = (const cs_lnum_t *)cs_get_device_ptr_const
                             (const_cast<cs_lnum_t *>(ms->row_index));
    const cs_lnum_t *__restrict__ d_col_id
      = (const cs_lnum_t *)cs_get_device_ptr_const
                             (const_cast<cs_lnum_t *>(ms->col_id));
    const cs_real_t *__restrict__ d_val
      = (const cs_real_t *)cs_get_device_ptr_const
                             (const_cast<cs_real_t *>(mc->val));

    unsigned int blocksize = 256;
    unsigned int gridsize
      = (unsigned int)ceil((double)ms->n_rows / blocksize);

    _mat_vect_p_l_csr_substract_diag<<<gridsize, blocksize, 0, _stream>>>
      (ms->n_rows, d_row_index, d_col_id, d_val,
       (const cs_real_t *)d_x, (cs_real_t *)d_y);

  }

  if (_stream == 0) {
    CS_HIP_CHECK(hipStreamSynchronize(_stream));
    CS_HIP_CHECK(hipGetLastError());
  }
}

#endif /* defined(HAVE_HIPSPARSE) */

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
                        cs_real_t     d_y[])
{
  const cs_matrix_struct_dist_t *ms
    = (const cs_matrix_struct_dist_t *)matrix->structure;
  const cs_matrix_coeff_t *mc
    = (const cs_matrix_coeff_t *)matrix->coeffs;

  const cs_lnum_t *__restrict__ row_index
    = (const cs_lnum_t *)cs_get_device_ptr_const
                           (const_cast<cs_lnum_t *>(ms->e.row_index));
  const cs_lnum_t *__restrict__ col_id
    = (const cs_lnum_t *)cs_get_device_ptr_const
                          (const_cast<cs_lnum_t *>(ms->e.col_id));

  const cs_real_t *__restrict__ d_val
    = (const cs_real_t *)cs_get_device_ptr_const
                           (const_cast<cs_real_t *>(mc->d_val));
  const cs_real_t *__restrict__ x_val
    = (const cs_real_t *)cs_get_device_ptr_const
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
    CS_HIP_CHECK(hipStreamSynchronize(_stream));
    CS_HIP_CHECK(hipGetLastError());
  }
}

#if defined(HAVE_HIPSPARSE)

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
cs_matrix_spmv_hip_msr_cusparse(cs_matrix_t  *matrix,
                                 bool          exclude_diag,
                                 bool          sync,
                                 cs_real_t     d_x[],
                                 cs_real_t     d_y[])
{
  cs_matrix_cusparse_map_t *csm
    = (cs_matrix_cusparse_map_t *)matrix->ext_lib_map;

  if (csm == nullptr) {
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
    const cs_matrix_coeff_t *mc
      = (const cs_matrix_coeff_t *)matrix->coeffs;
    const cs_real_t *__restrict__ d_val
      = (const cs_real_t *)cs_get_device_ptr_const
                             (const_cast<cs_real_t *>(mc->d_val));

    unsigned int blocksize = 256;
    unsigned int gridsize
      = (unsigned int)ceil((double)ms->n_rows / blocksize);

    _mat_vect_p_l_msr_diag<<<gridsize, blocksize, 0, _stream>>>
      (ms->n_rows, d_val, (const cs_real_t *)d_x, (cs_real_t *)d_y);

    beta = 1.;

  }

  hipsparseSetStream(_handle, _stream);

  hipsparseStatus_t status = HIPSPARSE_STATUS_SUCCESS;

  hipsparseSpMVAlg_t spmv_alg_type = HIPSPARSE_SPMV_ALG_DEFAULT;

  hipDataType val_dtype
    = (sizeof(cs_real_t) == 8) ? HIP_R_64F : HIP_R_32F;

  status = hipsparseSpMV(_handle,
                        HIPSPARSE_OPERATION_NON_TRANSPOSE,
                        &alpha,
                        csm->matA,
                        csm->vecX,
                        &beta,
                        csm->vecY,
                        val_dtype,
                        spmv_alg_type,
                        csm->dBuffer);

  if (HIPSPARSE_STATUS_SUCCESS != status)
    bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
              __func__, hipsparseGetErrorString(status));

  if (_stream == 0) {
    CS_HIP_CHECK(hipStreamSynchronize(_stream));
    CS_HIP_CHECK(hipGetLastError());
  }
}

#endif /* defined(HAVE_HIPSPARSE) */

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
                          cs_real_t     d_y[])
{
  const cs_matrix_struct_dist_t *ms
    = (const cs_matrix_struct_dist_t *)matrix->structure;
  const cs_matrix_coeff_t *mc
    = (const cs_matrix_coeff_t *)matrix->coeffs;

  const cs_lnum_t *__restrict__ row_index
    = (const cs_lnum_t *)cs_get_device_ptr_const
                           (const_cast<cs_lnum_t *>(ms->e.row_index));
  const cs_lnum_t *__restrict__ col_id
    = (const cs_lnum_t *)cs_get_device_ptr_const
                           (const_cast<cs_lnum_t *>(ms->e.col_id));

  const cs_real_t *__restrict__ d_val
    = (const cs_real_t *)cs_get_device_ptr_const
                           (const_cast<cs_real_t *>(mc->d_val));
  const cs_real_t *__restrict__ x_val
    = (const cs_real_t *)cs_get_device_ptr_const
                           (const_cast<cs_real_t *>(mc->e_val));

  /* Ghost cell communication */

  if (sync) {
    cs_halo_state_t *hs = _pre_vector_multiply_sync_x_start(matrix, d_x);
    _pre_vector_multiply_sync_x_end(matrix, hs, d_x);
  }

  /* Compute SpMV */

  cs_lnum_t n_r_rows = ms->n_rows*matrix->db_size;
  unsigned int blocksize = 128;
  unsigned int gridsize = cs_hip_grid_size(n_r_rows, blocksize);

  if (!exclude_diag) {

    if (matrix->db_size == 3)
      _b_3_3_mat_vect_p_l_msr<<<gridsize, blocksize, 0, _stream>>>
        (n_r_rows, col_id, row_index, d_val, x_val, d_x, d_y);
    else if (matrix->db_size == 6)
      _b_mat_vect_p_l_msr<6><<<gridsize, blocksize, 0, _stream>>>
        (n_r_rows, col_id, row_index, d_val, x_val, d_x, d_y);
    else if (matrix->db_size == 9)
      _b_mat_vect_p_l_msr<9><<<gridsize, blocksize, 0, _stream>>>
        (n_r_rows, col_id, row_index, d_val, x_val, d_x, d_y);
    else
      bft_error(__FILE__, __LINE__, 0, _("%s: block size %d not implemented."),
                __func__, (int)matrix->db_size);

  }
  else {

    if (matrix->db_size == 3)
      _b_3_3_mat_vect_p_l_msr_exdiag<<<gridsize, blocksize, 0, _stream>>>
        (n_r_rows, col_id, row_index, d_val, x_val, d_x, d_y);
    else if (matrix->db_size == 6)
      _b_mat_vect_p_l_msr_exdiag<6><<<gridsize, blocksize, 0, _stream>>>
        (n_r_rows, col_id, row_index, d_val, x_val, d_x, d_y);
    else if (matrix->db_size == 9)
      _b_mat_vect_p_l_msr_exdiag<9><<<gridsize, blocksize, 0, _stream>>>
        (n_r_rows, col_id, row_index, d_val, x_val, d_x, d_y);
    else
      bft_error(__FILE__, __LINE__, 0, _("%s: block size %d not implemented."),
                __func__, (int)matrix->db_size);

  }

  if (_stream == 0) {
    CS_HIP_CHECK(hipStreamSynchronize(_stream));
    CS_HIP_CHECK(hipGetLastError());
  }
}

#if defined(HAVE_HIPSPARSE)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with MSR matrix, block diagonal
 *        cuSPARSE version.
 *
 * Remark: this functions is available with older cuSPARSE versions not
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
cs_matrix_spmv_hip_msr_b_cusparse(cs_matrix_t  *matrix,
                                   bool          exclude_diag,
                                   bool          sync,
                                   cs_real_t     d_x[],
                                   cs_real_t     d_y[])
{
  cs_matrix_cusparse_map_t *csm
    = (cs_matrix_cusparse_map_t *)matrix->ext_lib_map;

  if (csm == nullptr) {
    matrix->ext_lib_map = _set_cusparse_map(matrix);
    csm = (cs_matrix_cusparse_map_t *)matrix->ext_lib_map;
  }

  /* Ghost cell communication */

  if (sync) {
    cs_halo_state_t *hs = _pre_vector_multiply_sync_x_start(matrix,
                                                            (cs_real_t *)d_x);
    _pre_vector_multiply_sync_x_end(matrix, hs, d_x);
  }

  _update_cusparse_map_block_diag(csm, matrix, d_x, d_y);

  cs_real_t alpha = 1.;
  cs_real_t beta = 0.;

  if (!exclude_diag) {

    const cs_matrix_struct_dist_t *ms
      = (const cs_matrix_struct_dist_t *)matrix->structure;
    const cs_matrix_coeff_t *mc
      = (const cs_matrix_coeff_t *)matrix->coeffs;
    const cs_real_t *__restrict__ d_val
      = (const cs_real_t *)cs_get_device_ptr_const
                             (const_cast<cs_real_t *>(mc->d_val));

    cs_lnum_t n_r_rows = ms->n_rows*matrix->db_size;
    unsigned int blocksize = 128;
    unsigned int gridsize = cs_hip_grid_size(n_r_rows, blocksize);

    if (matrix->db_size == 3)
      _b_3_3_spmv_diag<<<gridsize, blocksize, 0, _stream>>>
        (n_r_rows, d_val, (const cs_real_t *)d_x, (cs_real_t *)d_y);
    else if (matrix->db_size == 6)
      _b_spmv_diag<6><<<gridsize, blocksize, 0, _stream>>>
        (n_r_rows, d_val, (const cs_real_t *)d_x, (cs_real_t *)d_y);
    else if (matrix->db_size == 9)
      _b_spmv_diag<9><<<gridsize, blocksize, 0, _stream>>>
        (n_r_rows, d_val, (const cs_real_t *)d_x, (cs_real_t *)d_y);
    else
      bft_error(__FILE__, __LINE__, 0, _("%s: block size %d not implemented."),
                __func__, (int)matrix->db_size);

    beta = 1.;

  }

  hipsparseSetStream(_handle, _stream);

  hipsparseStatus_t status = HIPSPARSE_STATUS_SUCCESS;

  hipDataType val_dtype
    = (sizeof(cs_real_t) == 8) ? HIP_R_64F : HIP_R_32F;

  status = hipsparseSpMM(_handle,
                        HIPSPARSE_OPERATION_NON_TRANSPOSE,
                        HIPSPARSE_OPERATION_NON_TRANSPOSE,
                        &alpha,
                        csm->matA,
                        csm->matX,
                        &beta,
                        csm->matY,
                        val_dtype,
                        HIPSPARSE_SPMM_ALG_DEFAULT,
                        csm->dBuffer);

  if (HIPSPARSE_STATUS_SUCCESS != status)
    bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
              __func__, hipsparseGetErrorString(status));

  if (_stream == 0) {
    CS_HIP_CHECK(hipStreamSynchronize(_stream));
    CS_HIP_CHECK(hipGetLastError());
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with MSR matrix, block
 *        cuSPARSE version.
 *
 * Remark: this functions is available with older cuSPARSE versions not
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
cs_matrix_spmv_hip_msr_bb_cusparse(cs_matrix_t  *matrix,
                                    bool          exclude_diag,
                                    bool          sync,
                                    cs_real_t     d_x[],
                                    cs_real_t     d_y[])
{
  cs_matrix_cusparse_map_t *csm
    = (cs_matrix_cusparse_map_t *)matrix->ext_lib_map;

  if (csm == nullptr) {
    matrix->ext_lib_map = _set_cusparse_map(matrix);
    csm = (cs_matrix_cusparse_map_t *)matrix->ext_lib_map;
  }

  /* Ghost cell communication */

  if (sync) {
    cs_halo_state_t *hs = _pre_vector_multiply_sync_x_start(matrix,
                                                            (cs_real_t *)d_x);
    _pre_vector_multiply_sync_x_end(matrix, hs, d_x);
  }

  /* no update_cusparse_map type function call here as only
     the non-generic API is available here */

  cs_real_t alpha = 1.;
  cs_real_t beta = 0.;

  if (!exclude_diag) {

    const cs_matrix_struct_dist_t *ms
      = (const cs_matrix_struct_dist_t *)matrix->structure;
    const cs_matrix_coeff_t *mc
      = (const cs_matrix_coeff_t *)matrix->coeffs;
    const cs_real_t *__restrict__ d_val
      = (const cs_real_t *)cs_get_device_ptr_const
                             (const_cast<cs_real_t *>(mc->d_val));

    cs_lnum_t n_r_rows = ms->n_rows*matrix->db_size;
    unsigned int blocksize = 128;
    unsigned int gridsize = cs_hip_grid_size(n_r_rows, blocksize);

    if (matrix->db_size == 3)
      _b_3_3_spmv_diag<<<gridsize, blocksize, 0, _stream>>>
        (n_r_rows, d_val, (const cs_real_t *)d_x, (cs_real_t *)d_y);
    else if (matrix->db_size == 6)
      _b_spmv_diag<6><<<gridsize, blocksize, 0, _stream>>>
        (n_r_rows, d_val, (const cs_real_t *)d_x, (cs_real_t *)d_y);
    else if (matrix->db_size == 9)
      _b_spmv_diag<9><<<gridsize, blocksize, 0, _stream>>>
        (n_r_rows, d_val, (const cs_real_t *)d_x, (cs_real_t *)d_y);
    else
      bft_error(__FILE__, __LINE__, 0, _("%s: block size %d not implemented."),
                __func__, (int)matrix->db_size);

    beta = 1.;

  }

  hipsparseSetStream(_handle, _stream);

  hipsparseStatus_t status = HIPSPARSE_STATUS_SUCCESS;

  if (sizeof(cs_real_t) == 8) {
    double _alpha = alpha;
    double _beta = beta;
    status = hipsparseDbsrmv(_handle,
                            HIPSPARSE_DIRECTION_ROW,
                            HIPSPARSE_OPERATION_NON_TRANSPOSE,
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

    status = hipsparseSbsrmv(_handle,
                            HIPSPARSE_DIRECTION_ROW,
                            HIPSPARSE_OPERATION_NON_TRANSPOSE,
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

  if (HIPSPARSE_STATUS_SUCCESS != status)
    bft_error(__FILE__, __LINE__, 0, _("%s: cuSPARSE error %d."),
              __func__, (int)status);

  if (_stream == 0) {
    CS_HIP_CHECK(hipStreamSynchronize(_stream));
    CS_HIP_CHECK(hipGetLastError());
  }
}

#endif /* defined(HAVE_HIPSPARSE) */

/*----------------------------------------------------------------------------*/

END_C_DECLS
