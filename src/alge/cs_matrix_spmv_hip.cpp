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

#if defined(HAVE_ROCSPARSE)
#include <rocsparse/rocsparse.h>
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

#if defined(HAVE_ROCSPARSE)

/* Mapping of matrix coefficients and structure to rocSPARSE */
/*-----------------------------------------------------------*/

typedef struct _cs_matrix_rocsparse_map_t {

  bool  block_diag;              /* Use identity blocks diagonal structure ? */

  rocsparse_spmat_descr  mat_a;  /* Handle to rocSPARSE Matrix */

  rocsparse_dnmat_descr  mat_x;  /* Handle to rocSPARSE Matrix (blocked vector) */
  rocsparse_dnmat_descr  mat_y;  /* Handle to rocSPARSE Matrix (blocked vector) */

  rocsparse_dnvec_descr  vec_x;  /* Handle to rocSPARSE Vector */
  rocsparse_dnvec_descr  vec_y;  /* Handle to rocSPARSE output Vector */

  void  *vec_x_values;           /* Pointer to vector values */
  void  *vec_y_values;           /* Pointer to vector values */

  size_t  buffer_size;           /* SpMV buffer size */
  void   *buffer;                /* SpMV buffer */

  rocsparse_spmv_descr spmv_descr;  /* SpMV descriptor */

  int  nnz;                      /* Number of nonzeroes */

} cs_matrix_rocsparse_map_t;

#endif // defined(HAVE_ROCSPARSE)

/*============================================================================
 *  Global variables
 *============================================================================*/

static hipStream_t _stream = 0;

#if defined(HAVE_ROCSPARSE)

static rocsparse_handle  _handle = nullptr;

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

#pragma unroll 2
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
 * This can be combined with a rocSPARSE CSR SpMV product with the
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
 * - If MPI is ROC-aware, no values need to go through the host
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

#if defined(HAVE_ROCSPARSE)

/*----------------------------------------------------------------------------
 * Unset matrix rocSPARSE mapping.
 *
 * parameters:
 *   matrix    <-- pointer to matrix structure
 *----------------------------------------------------------------------------*/

static void
_unset_rocsparse_map(cs_matrix_t   *matrix)
{
  cs_matrix_rocsparse_map_t *csm
    = (cs_matrix_rocsparse_map_t *)matrix->ext_lib_map;

  if (csm == nullptr)
    return;

  rocsparse_destroy_spmv_descr(csm->spmv_descr);
  rocsparse_destroy_spmat_descr(csm->mat_a);

  if (csm->block_diag == false) {
    if (csm->vec_x_values != nullptr)
      rocsparse_destroy_dnvec_descr(csm->vec_x);
    if (csm->vec_y_values != nullptr)
      rocsparse_destroy_dnvec_descr(csm->vec_y);
  }
  else {
    if (csm->vec_x_values != nullptr)
      rocsparse_destroy_dnmat_descr(csm->mat_x);
    if (csm->vec_y_values != nullptr)
      rocsparse_destroy_dnmat_descr(csm->mat_y);
  }

  CS_FREE(csm->buffer);
  csm->buffer_size = 0;

  csm->block_diag = false;

  csm->vec_x_values = nullptr;
  csm->vec_y_values = nullptr;

  csm->nnz = 0;

  CS_FREE(matrix->ext_lib_map);
  matrix->destroy_adaptor = nullptr;
}

/*----------------------------------------------------------------------------
 * Update matrix rocSPARSE mapping.
 *
 * parameters:
 *   csm       <-> rocSPARSE matrix mapping
 *   matrix    <-- pointer to matrix structure
 *   d_x       <-- pointer to input vector (on device)
 *   d_y       <-- pointer to output vector (on device)
 *----------------------------------------------------------------------------*/

static void
_update_rocsparse_map(cs_matrix_rocsparse_map_t  *csm,
                      const cs_matrix_t          *matrix,
                      void                       *d_x,
                      void                       *d_y)
{
  assert(csm != nullptr);

  rocsparse_status status = rocsparse_status_success;
  rocsparse_datatype val_dtype = (sizeof(cs_real_t) == 8) ?
    rocsparse_datatype_f64_r : rocsparse_datatype_f32_r;

  if (matrix->eb_size == matrix->db_size) {

    if (d_x != csm->vec_x_values) {
      if (csm->vec_x_values != nullptr)
        status = rocsparse_destroy_dnvec_descr(csm->vec_x);

      if (rocsparse_status_success == status)
        status = rocsparse_create_dnvec_descr(&(csm->vec_x),
                                              matrix->n_cols_ext*matrix->db_size,
                                              d_x,
                                              val_dtype);
      if (rocsparse_status_success == status)
        csm->vec_x_values = d_x;
      else
        if (rocsparse_status_success != status)
          bft_error(__FILE__, __LINE__, 0, _("%s: rocSPARSE status %d."),
                    __func__, status);
    }

    if (d_y != csm->vec_y_values) {
      if (csm->vec_y_values != nullptr)
        status = rocsparse_destroy_dnvec_descr(csm->vec_y);

      if (rocsparse_status_success == status)
        status = rocsparse_create_dnvec_descr(&(csm->vec_y),
                                              matrix->n_rows*matrix->db_size,
                                              d_y,
                                              val_dtype);

      if (rocsparse_status_success == status)
        csm->vec_y_values = d_y;
      else
        bft_error(__FILE__, __LINE__, 0, _("%s: rocSPARSE status %d."),
                  __func__, status);
    }

  }

  else if (matrix->eb_size > 1) {

    if (d_x != csm->vec_x_values) {
      if (csm->vec_x_values != nullptr)
        status = rocsparse_destroy_dnmat_descr(csm->mat_x);

      if (rocsparse_status_success == status)
        status = rocsparse_create_dnmat_descr(&(csm->mat_x),
                                              matrix->n_cols_ext,
                                              matrix->db_size,
                                              matrix->db_size,
                                              d_x,
                                              val_dtype,
                                              rocsparse_order_row);

      if (rocsparse_status_success == status)
        csm->vec_x_values = d_x;
      else
        bft_error(__FILE__, __LINE__, 0, _("%s: rocSPARSE status %d."),
                  __func__, status);
    }

    if (d_y != csm->vec_y_values) {
      if (csm->vec_y_values != nullptr)
        status = rocsparse_destroy_dnmat_descr(csm->mat_y);

      if (rocsparse_status_success == status)
        status = rocsparse_create_dnmat_descr(&(csm->mat_y),
                                              matrix->n_rows,
                                              matrix->db_size,
                                              matrix->db_size,
                                              d_y,
                                              val_dtype,
                                              rocsparse_order_row);

      if (rocsparse_status_success == status)
        csm->vec_y_values = d_y;
      else
        bft_error(__FILE__, __LINE__, 0, _("%s: rocSPARSE status %d."),
                  __func__, status);
    }

  }
}

/*----------------------------------------------------------------------------
 * Set matrix rocSPARSE mapping.
 *
 * parameters:
 *   matrix    <-- pointer to matrix structure
 *----------------------------------------------------------------------------*/

static cs_matrix_rocsparse_map_t *
_set_rocsparse_map(cs_matrix_t   *matrix,
                   void          *d_x,
                   void          *d_y)
{
  cs_matrix_rocsparse_map_t *csm
    = (cs_matrix_rocsparse_map_t *)matrix->ext_lib_map;

  if (csm != nullptr) {
    _unset_rocsparse_map(matrix);
  }
  else {
    CS_MALLOC(csm, 1, cs_matrix_rocsparse_map_t);
    matrix->ext_lib_map = (void *)csm;
  }
  matrix->destroy_adaptor = _unset_rocsparse_map;

  if (matrix->eb_size == matrix->db_size)
    return csm;
  // Continued in _update_rocsparse_map_block_diag if the above is false.

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

  rocsparse_status status = rocsparse_status_success;

  /* For beta, 0 should be enough for SmPV, 1 is needed for
     y = A.x + b.y which is useful when y is initialized by a
     separate diagonal da.x product */

  const cs_real_t alpha = 1.0;
  const cs_real_t beta = (matrix->type == CS_MATRIX_CSR) ? 0. : 1.;

  if (_handle == nullptr)
    status = rocsparse_create_handle(&_handle);

  if (matrix->db_size > matrix->eb_size)
    csm->block_diag = true;
  else
    csm->block_diag = false;

  if (rocsparse_status_success != status)
    bft_error(__FILE__, __LINE__, 0, _("%s: rocSPARSE status %d."),
              __func__, status);

  csm->vec_x_values = nullptr;  /* Pointer to vector values */
  csm->vec_y_values = nullptr;  /* Pointer to vector values */
  csm->buffer = nullptr;

  rocsparse_indextype index_dtype = (sizeof(cs_lnum_t) == 4) ?
    rocsparse_indextype_i32 : rocsparse_indextype_i64;
  rocsparse_datatype val_dtype = (sizeof(cs_real_t) == 8) ?
    rocsparse_datatype_f64_r : rocsparse_datatype_f32_r;

  if (matrix->eb_size == 1) {
    status = rocsparse_create_csr_descr(&(csm->mat_a),
                                        matrix->n_rows,
                                        matrix->n_cols_ext,
                                        nnz,
                                        const_cast<void *>(row_index),
                                        const_cast<void *>(col_id),
                                        const_cast<void *>(e_val),
                                        index_dtype,
                                        index_dtype,
                                        rocsparse_index_base_zero,
                                        val_dtype);
  }
  else {
    status = rocsparse_create_bsr_descr(&(csm->mat_a),
                                        matrix->n_rows,
                                        matrix->n_cols_ext,
                                        nnz,
                                        rocsparse_direction_row,
                                        matrix->eb_size,
                                        const_cast<void *>(row_index),
                                        const_cast<void *>(col_id),
                                        const_cast<void *>(e_val),
                                        index_dtype,
                                        index_dtype,
                                        rocsparse_index_base_zero,
                                        val_dtype);
  }


  if (status == rocsparse_status_success)
    status = rocsparse_create_spmv_descr(&(csm->spmv_descr));

  if (rocsparse_status_success != status)
    bft_error(__FILE__, __LINE__, 0, _("%s: rocSPARSE status %d."),
              __func__, status);

  rocsparse_error p_error[1] = {};

  if (status == rocsparse_status_success) {
    const rocsparse_spmv_alg spmv_alg = rocsparse_spmv_alg_csr_adaptive;
    status = rocsparse_spmv_set_input(_handle,
                                      csm->spmv_descr,
                                      rocsparse_spmv_input_alg,
                                      &spmv_alg,
                                      sizeof(spmv_alg),
                                      p_error);
  }

  if (status == rocsparse_status_success) {
    const rocsparse_operation spmv_operation = rocsparse_operation_none;
    status = rocsparse_spmv_set_input(_handle,
                                      csm->spmv_descr,
                                      rocsparse_spmv_input_operation,
                                      &spmv_operation,
                                      sizeof(spmv_operation),
                                      p_error);
  }

  if (status == rocsparse_status_success) {
    const rocsparse_datatype val_dtype = (sizeof(cs_real_t) == 8) ?
      rocsparse_datatype_f64_r : rocsparse_datatype_f32_r;
    status = rocsparse_spmv_set_input(_handle,
                                      csm->spmv_descr,
                                      rocsparse_spmv_input_compute_datatype,
                                      &val_dtype,
                                      sizeof(val_dtype),
                                      p_error);
  }

  _update_rocsparse_map(csm, matrix, d_x, d_y);

  if (status == rocsparse_status_success) {
    size_t buffer_size;
    status = rocsparse_v2_spmv_buffer_size(_handle,
                                           csm->spmv_descr,
                                           csm->mat_a,
                                           csm->vec_x,
                                           csm->vec_y,
                                           rocsparse_v2_spmv_stage_analysis,
                                           &buffer_size,
                                           p_error);
  }

  char *buffer = nullptr;
  size_t buffer_size_a = 0, buffer_size_c = 0;

  if (status == rocsparse_status_success) {
    status = rocsparse_v2_spmv_buffer_size(_handle,
                                           csm->spmv_descr,
                                           csm->mat_a,
                                           csm->vec_x,
                                           csm->vec_y,
                                           rocsparse_v2_spmv_stage_analysis,
                                           &buffer_size_a,
                                           p_error);
  }

  if (status == rocsparse_status_success) {
    CS_MALLOC_HD(buffer, buffer_size_a, char, CS_ALLOC_DEVICE);

    status = rocsparse_v2_spmv(_handle,
                               csm->spmv_descr,
                               &alpha,
                               csm->mat_a,
                               csm->vec_x,
                               &beta,
                               csm->vec_y,
                               rocsparse_v2_spmv_stage_analysis,
                               buffer_size_a,
                               buffer,
                               p_error);
  }

  if (status == rocsparse_status_success) {
    status = rocsparse_v2_spmv_buffer_size(_handle,
                                           csm->spmv_descr,
                                           csm->mat_a,
                                           csm->vec_x,
                                           csm->vec_y,
                                           rocsparse_v2_spmv_stage_compute,
                                           &buffer_size_c,
                                           p_error);
  }

  if (status == rocsparse_status_success) {
    if (buffer_size_a != buffer_size_c) {
      CS_FREE(buffer);
      CS_MALLOC_HD(buffer, buffer_size_c, char, CS_ALLOC_DEVICE);
    }
  }
  else
    CS_FREE(buffer);

  csm->buffer_size = buffer_size_c;
  csm->buffer = buffer;

  if (rocsparse_status_success != status)
    bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
              __func__, rocsparse_error_get_message(p_error[0]));

  return csm;
}

/*----------------------------------------------------------------------------
 * Update matrix rocSPARSE mapping in block diagonal case.
 *
 * parameters:
 *   csm       <-> rocSPARSE matrix mapping
 *   matrix    <-- pointer to matrix structure
 *   d_x       <-- pointer to input vector (on device)
 *   d_y       <-- pointer to output vector (on device)
 *----------------------------------------------------------------------------*/

static void
_update_rocsparse_map_block_diag(cs_matrix_rocsparse_map_t  *csm,
                                 const cs_matrix_t          *matrix,
                                 void                       *d_x,
                                 void                       *d_y)
{
  assert(csm != nullptr);

  rocsparse_status status = rocsparse_status_success;
  rocsparse_datatype val_dtype = (sizeof(cs_real_t) == 8) ?
    rocsparse_datatype_f64_r : rocsparse_datatype_f32_r;

  if (d_x != csm->vec_x_values) {
    if (csm->vec_x_values != nullptr)
      rocsparse_destroy_dnmat_descr(csm->mat_x);

    status = rocsparse_create_dnmat_descr(&(csm->mat_x),
                                          matrix->n_cols_ext,
                                          matrix->db_size,
                                          matrix->db_size,
                                          d_x,
                                          val_dtype,
                                          rocsparse_order_row);

    if (rocsparse_status_success != status)
      bft_error(__FILE__, __LINE__, 0, _("%s: rocSPARSE status %d."),
                __func__, status);

    csm->vec_x_values = d_x;
  }

  if (d_y != csm->vec_y_values) {
    if (csm->vec_y_values != nullptr)
      rocsparse_destroy_dnmat_descr(csm->mat_y);

    status = rocsparse_create_dnmat_descr(&(csm->mat_y),
                                          matrix->n_rows,
                                          matrix->db_size,
                                          matrix->db_size,
                                          d_y,
                                          val_dtype,
                                          rocsparse_order_row);

    if (rocsparse_status_success != status)
      bft_error(__FILE__, __LINE__, 0, _("%s: rocSPARSE status %d."),
                __func__, status);

    csm->vec_y_values = d_y;
  }

  if (csm->buffer == nullptr) {
    size_t buffer_size = 0;
    cs_real_t alpha = 1.0;
    cs_real_t beta = 1.0;  /* 0 should be enough for SmPV, 1 needed for
                              y = A.x + b.y
                              which is useful when y is initialized by
                              a separate diagonal da.x product */

    status = rocsparse_spmm(_handle,
                            rocsparse_operation_none,
                            rocsparse_operation_none,
                            &alpha,
                            csm->mat_a,
                            csm->mat_x,
                            &beta,
                            csm->mat_y,
                            val_dtype,
                            rocsparse_spmm_alg_default,
                            rocsparse_spmm_stage_buffer_size,
                            &buffer_size,
                            nullptr);

    if (rocsparse_status_success != status)
      bft_error(__FILE__, __LINE__, 0, _("%s: rocSPARSE status %d."),
                __func__, status);

    char *buffer;
    CS_MALLOC_HD(buffer, buffer_size, char, CS_ALLOC_DEVICE);
    csm->buffer_size = buffer_size;
    csm->buffer = buffer;

    status = rocsparse_spmm(_handle,
                            rocsparse_operation_none,
                            rocsparse_operation_none,
                            &alpha,
                            csm->mat_a,
                            csm->mat_x,
                            &beta,
                            csm->mat_y,
                            val_dtype,
                            rocsparse_spmm_alg_default,
                            rocsparse_spmm_stage_preprocess,
                            &buffer_size,
                            csm->buffer);

    if (rocsparse_status_success != status)
      bft_error(__FILE__, __LINE__, 0, _("%s: rocSPARSE status %d."),
                __func__, status);
  }
}

#endif // defined(HAVE_ROCSPARSE)

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize ROC matrix API.
 *
 * This frees resources such as the rocSPARSE handle, if used.
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_hip_finalize(void)
{
  _stream = 0;

#if defined(HAVE_ROCSPARSE)

  if (_handle != nullptr) {
    rocsparse_destroy_handle(_handle);
    _handle = nullptr;
  }

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign ROC stream for next ROC-based SpMV operations.
 *
 * If a stream other than the default stream (0) is used, it will not be
 * synchronized automatically after sparse matrix-vector products (so as to
 * avoid the corresponding overhead), so the caller will need to manage
 * stream syncronization manually.
 *
 * This function is callable only from ROC code.
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_hip_set_stream(hipStream_t  stream)
{
  _stream = stream;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return stream used for ROC-based SpMV operations.
 *
 * This function is callable only from ROC code.
 */
/*----------------------------------------------------------------------------*/

hipStream_t
cs_matrix_spmv_hip_get_stream(void)
{
  return _stream;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with MSR matrix, scalar ROC version.
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
 * \brief Matrix.vector product y = A.x with CSR matrix, scalar ROC version.
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
    = (const cs_matrix_coeff_t *)matrix->coeffs;

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

#if defined(HAVE_ROCSPARSE)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with CSR matrix,
 *  scalar rocSPARSE version.
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
                                 cs_real_t     d_y[])
{
  cs_matrix_rocsparse_map_t *csm
    = (cs_matrix_rocsparse_map_t *)matrix->ext_lib_map;

  if (csm == nullptr) {
    matrix->ext_lib_map = _set_rocsparse_map(matrix, d_x, d_y);
    csm = (cs_matrix_rocsparse_map_t *)matrix->ext_lib_map;
  }

  /* Ghost cell communication */

  if (sync) {
    cs_halo_state_t *hs = _pre_vector_multiply_sync_x_start(matrix,
                                                            (cs_real_t *)d_x);
    cs_halo_sync_wait(matrix->halo, (cs_real_t *)d_x, hs);
  }

  _update_rocsparse_map(csm, matrix, d_x, d_y);

  cs_real_t alpha = 1.0;
  cs_real_t beta = 0.0;

  rocsparse_set_stream(_handle, _stream);

  rocsparse_status status = rocsparse_status_success;
  rocsparse_error p_error[1] = {};

  status = rocsparse_v2_spmv(_handle,
                             csm->spmv_descr,
                             &alpha,
                             csm->mat_a,
                             csm->vec_x,
                             &beta,
                             csm->vec_y,
                             rocsparse_v2_spmv_stage_compute,
                             csm->buffer_size,
                             csm->buffer,
                             p_error);

  if (rocsparse_status_success != status)
    bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
              __func__, rocsparse_error_get_message(p_error[0]));

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

#endif /* defined(HAVE_ROCSPARSE) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with MSR matrix, scalar ROC version.
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

#if defined(HAVE_ROCSPARSE)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with MSR matrix, scalar rocSPARSE version.
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
                                 cs_real_t     d_y[])
{
  cs_matrix_rocsparse_map_t *csm
    = (cs_matrix_rocsparse_map_t *)matrix->ext_lib_map;

  if (csm == nullptr) {
    matrix->ext_lib_map = _set_rocsparse_map(matrix, d_x, d_y);
    csm = (cs_matrix_rocsparse_map_t *)matrix->ext_lib_map;
  }

  /* Ghost cell communication */

  if (sync) {
    cs_halo_state_t *hs = _pre_vector_multiply_sync_x_start(matrix,
                                                            (cs_real_t *)d_x);
    cs_halo_sync_wait(matrix->halo, (cs_real_t *)d_x, hs);
  }

  _update_rocsparse_map(csm, matrix, d_x, d_y);

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

  rocsparse_set_stream(_handle, _stream);

  rocsparse_status status = rocsparse_status_success;
  rocsparse_error p_error[1] = {};

  status = rocsparse_v2_spmv(_handle,
                             csm->spmv_descr,
                             &alpha,
                             csm->mat_a,
                             csm->vec_x,
                             &beta,
                             csm->vec_y,
                             rocsparse_v2_spmv_stage_compute,
                             csm->buffer_size,
                             csm->buffer,
                             p_error);

  if (rocsparse_status_success != status)
    bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
              __func__, rocsparse_error_get_message(p_error[0]));

  if (_stream == 0) {
    CS_HIP_CHECK(hipStreamSynchronize(_stream));
    CS_HIP_CHECK(hipGetLastError());
  }
}

#endif /* defined(HAVE_ROCSPARSE) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with MSR matrix, block diagonal
 *        ROC version.
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
                                   cs_real_t     d_y[])
{
  cs_matrix_rocsparse_map_t *csm
    = (cs_matrix_rocsparse_map_t *)matrix->ext_lib_map;

  if (csm == nullptr) {
    matrix->ext_lib_map = _set_rocsparse_map(matrix, d_x, d_y);
    csm = (cs_matrix_rocsparse_map_t *)matrix->ext_lib_map;
  }

  /* Ghost cell communication */

  if (sync) {
    cs_halo_state_t *hs = _pre_vector_multiply_sync_x_start(matrix,
                                                            (cs_real_t *)d_x);
    _pre_vector_multiply_sync_x_end(matrix, hs, d_x);
  }

  _update_rocsparse_map_block_diag(csm, matrix, d_x, d_y);

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

  rocsparse_set_stream(_handle, _stream);

  rocsparse_status status = rocsparse_status_success;

  rocsparse_datatype val_dtype = (sizeof(cs_real_t) == 8) ?
    rocsparse_datatype_f64_r : rocsparse_datatype_f32_r;

  status = rocsparse_spmm(_handle,
                          rocsparse_operation_none,
                          rocsparse_operation_none,
                          &alpha,
                          csm->mat_a,
                          csm->mat_x,
                          &beta,
                          csm->mat_y,
                          val_dtype,
                          rocsparse_spmm_alg_default,
                          rocsparse_spmm_stage_compute,
                          &(csm->buffer_size),
                          csm->buffer);

  if (rocsparse_status_success != status)
    bft_error(__FILE__, __LINE__, 0, _("%s: rocSPARSE status %d."),
              __func__, status);

  if (_stream == 0) {
    CS_HIP_CHECK(hipStreamSynchronize(_stream));
    CS_HIP_CHECK(hipGetLastError());
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with MSR matrix, block
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
cs_matrix_spmv_hip_msr_bb_rocsparse(cs_matrix_t  *matrix,
                                    bool          exclude_diag,
                                    bool          sync,
                                    cs_real_t     d_x[],
                                    cs_real_t     d_y[])
{
  cs_matrix_rocsparse_map_t *csm
    = (cs_matrix_rocsparse_map_t *)matrix->ext_lib_map;

  if (csm == nullptr) {
    matrix->ext_lib_map = _set_rocsparse_map(matrix, d_x, d_y);
    csm = (cs_matrix_rocsparse_map_t *)matrix->ext_lib_map;
  }

  /* Ghost cell communication */

  if (sync) {
    cs_halo_state_t *hs = _pre_vector_multiply_sync_x_start(matrix,
                                                            (cs_real_t *)d_x);
    _pre_vector_multiply_sync_x_end(matrix, hs, d_x);
  }

  /* no update_rocsparse_map type function call here as only
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

  rocsparse_set_stream(_handle, _stream);

  rocsparse_status status = rocsparse_status_success;
  rocsparse_error p_error[1] = {};

  status = rocsparse_v2_spmv(_handle,
                             csm->spmv_descr,
                             &alpha,
                             csm->mat_a,
                             csm->vec_x,
                             &beta,
                             csm->vec_y,
                             rocsparse_v2_spmv_stage_compute,
                             csm->buffer_size,
                             csm->buffer,
                             p_error);

  if (rocsparse_status_success != status)
    bft_error(__FILE__, __LINE__, 0, _("%s: %s."),
              __func__, rocsparse_error_get_message(p_error[0]));

  if (_stream == 0) {
    CS_HIP_CHECK(hipStreamSynchronize(_stream));
    CS_HIP_CHECK(hipGetLastError());
  }
}

#endif /* defined(HAVE_ROCSPARSE) */

/*----------------------------------------------------------------------------*/

END_C_DECLS
