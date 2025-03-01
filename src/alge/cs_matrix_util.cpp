/*============================================================================
 * Utilitary functions for sparse matrixes.
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
 * Standard C library headers
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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "fvm/fvm_io_num.h"

#include "base/cs_base.h"
#include "alge/cs_blas.h"
#include "base/cs_file.h"
#include "base/cs_halo.h"
#include "base/cs_halo_perio.h"
#include "base/cs_log.h"
#include "base/cs_mem.h"
#include "base/cs_numbering.h"
#include "base/cs_order.h"
#include "base/cs_parall.h"
#include "base/cs_part_to_block.h"
#include "base/cs_prototypes.h"
#include "base/cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "alge/cs_matrix_util.h"

#include "alge/cs_matrix_priv.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * y[i] = abs(da[i]), with da possibly nullptr.
 *
 * parameters:
 *   da         <-- Pointer to coefficients array (usually matrix diagonal)
 *   dd         --> Resulting vector
 *   n_rows     <-- Number of rows
 *   n_cols_ext <-- Number of columns + ghost columns
*----------------------------------------------------------------------------*/

static void
_diag_dom_diag_contrib(const cs_real_t  *restrict da,
                       cs_real_t        *restrict dd,
                       cs_lnum_t         n_rows,
                       cs_lnum_t         n_cols_ext)
{
  cs_lnum_t  ii;

  if (da != nullptr) {
#   pragma omp parallel for
    for (ii = 0; ii < n_rows; ii++)
      dd[ii] = fabs(da[ii]);
    for (ii = n_rows; ii < n_cols_ext; ii++)
      dd[ii] = 0.0;
  }
  else {
#   pragma omp parallel for
    for (ii = 0; ii < n_cols_ext; ii++)
      dd[ii] = 0.0;
  }

}

/*----------------------------------------------------------------------------
 * Block diagonal contribution to diagonal dominance, with da possibly nullptr.
 *
 * parameters:
 *   da         <-- Pointer to coefficients array (usually matrix diagonal)
 *   dd         --> Resulting vector
 *   n_rows     <-- Number of rows
 *   n_cols_ext <-- Number of colmuns + ghost columns
 *   b_size     <-- block size
 *----------------------------------------------------------------------------*/

static void
_b_diag_dom_diag_contrib(const cs_real_t  *restrict da,
                         cs_real_t        *restrict dd,
                         cs_lnum_t         n_rows,
                         cs_lnum_t         n_cols_ext,
                         const cs_lnum_t   b_size)
{
  cs_lnum_t  ii, jj, kk;
  double  sign;
  const cs_lnum_t  b_size_2 = b_size * b_size;
  const cs_lnum_t  dd_size = n_cols_ext*b_size;

# pragma omp parallel for
  for (ii = 0; ii < dd_size; ii++)
    dd[ii] = 0.0;

  if (da != nullptr) {
#   pragma omp parallel for private(jj, kk, sign)
    for (ii = 0; ii < n_rows; ii++) {
      for (jj = 0; jj < b_size; jj++)
        dd[ii*b_size + jj] = 0.0;
      for (jj = 0; jj < b_size; jj++) {
        for (kk = 0; kk < b_size; kk++) {
          sign = (jj == kk) ? 1. : -1.;
          dd[ii*b_size + kk]
            += sign*fabs(da[ii*b_size_2 + jj*b_size + kk]);
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Normalize diagonal dominance, with da possibly nullptr.
 *
 * parameters:
 *   da         <-- Pointer to coefficients array (usually matrix diagonal)
 *   dd         --> Resulting vector
 *   n_rows     <-- Number of rows
*----------------------------------------------------------------------------*/

static void
_diag_dom_diag_normalize(const cs_real_t  *restrict da,
                         cs_real_t        *restrict dd,
                         cs_lnum_t         n_rows)
{
  cs_lnum_t  ii;

  if (da != nullptr) {
#   pragma omp parallel for
    for (ii = 0; ii < n_rows; ii++) {
      if (fabs(da[ii]) > 1.e-18)
        dd[ii] /= fabs(da[ii]);
      else if (dd[ii] > -1.e-18)
        dd[ii] = -1.e18;
      else
        dd[ii] = 0;
    }
  }
  else {
#   pragma omp parallel for
    for (ii = 0; ii < n_rows; ii++) {
      if (dd[ii] > -1.e-18)
        dd[ii] = -1.e18;
      else
        dd[ii] = 0;
    }
  }
}

/*----------------------------------------------------------------------------
 * Normalize block diagonal dominance, with da possibly nullptr.
 *
 * parameters:
 *   da      <-- Pointer to coefficients array (usually matrix diagonal)
 *   dd      --> Resulting vector
 *   n_rows  <-- Number of rows
 *   b_size  <-- block size
 *----------------------------------------------------------------------------*/

static void
_b_diag_dom_diag_normalize(const cs_real_t  *restrict da,
                           cs_real_t        *restrict dd,
                           cs_lnum_t         n_rows,
                           const cs_lnum_t   b_size)
{
  cs_lnum_t  ii, jj;
  double  d_val;

  if (da != nullptr) {
    const cs_lnum_t b_size_2 = b_size*b_size;
#   pragma omp parallel for private(jj, d_val)
    for (ii = 0; ii < n_rows; ii++) {
      for (jj = 0; jj < b_size; jj++) {
        d_val = fabs(da[ii*b_size_2 + jj*b_size + jj]);
        if (d_val > 1.e-18)
          dd[ii*b_size + jj] /= d_val;
        else if (dd[ii*b_size + jj] > -1.e-18)
          dd[ii*b_size + jj] = -1.e18;
        else
          dd[ii*b_size + jj] = 0;
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Measure Diagonal dominance of native matrix.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *   dd     --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_diag_dom_native(const cs_matrix_t  *matrix,
                 cs_real_t          *restrict dd)
{
  cs_lnum_t  ii, jj, edge_id;

  const auto ms =
    static_cast<const cs_matrix_struct_native_t *>(matrix->structure);
  const auto mc = static_cast<const cs_matrix_coeff_t *>(matrix->coeffs);

  const cs_real_t  *restrict xa = mc->e_val;

  /* diagonal contribution */

  _diag_dom_diag_contrib(mc->d_val, dd, ms->n_rows, ms->n_cols_ext);

  /* non-diagonal terms */

  if (mc->e_val != nullptr) {

    const cs_lnum_2_t *restrict face_cel_p = ms->edges;

    if (mc->symmetric) {

      for (edge_id = 0; edge_id < ms->n_edges; edge_id++) {
        ii = face_cel_p[edge_id][0];
        jj = face_cel_p[edge_id][1];
        dd[ii] -= fabs(xa[edge_id]);
        dd[jj] -= fabs(xa[edge_id]);
      }

    }
    else {

      for (edge_id = 0; edge_id < ms->n_edges; edge_id++) {
        ii = face_cel_p[edge_id][0];
        jj = face_cel_p[edge_id][1];
        dd[ii] -= fabs(xa[2*edge_id]);
        dd[jj] -= fabs(xa[2*edge_id + 1]);
      }

    }

  }

  _diag_dom_diag_normalize(mc->d_val, dd, ms->n_rows);
}

/*----------------------------------------------------------------------------
 * Measure Diagonal dominance of native matrix with block diagonal.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *   dd     --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_b_diag_dom_native(const cs_matrix_t  *matrix,
                   cs_real_t          *restrict dd)
{
  cs_lnum_t  ii, jj, kk, edge_id;

  const auto ms =
    static_cast<const cs_matrix_struct_native_t *>(matrix->structure);
  const auto mc = static_cast<const cs_matrix_coeff_t *>(matrix->coeffs);

  const cs_real_t  *restrict xa = mc->e_val;
  const cs_lnum_t db_size = matrix->db_size;

  /* block diagonal contribution */

  _b_diag_dom_diag_contrib(mc->d_val, dd, ms->n_rows, ms->n_cols_ext, db_size);

  /* non-diagonal terms */

  if (mc->e_val != nullptr) {

    const cs_lnum_2_t *restrict face_cel_p = ms->edges;

    if (mc->symmetric) {

      for (edge_id = 0; edge_id < ms->n_edges; edge_id++) {
        ii = face_cel_p[edge_id][0];
        jj = face_cel_p[edge_id][1];
        for (kk = 0; kk < db_size; kk++) {
          dd[ii*db_size + kk] -= fabs(xa[edge_id]);
          dd[jj*db_size + kk] -= fabs(xa[edge_id]);
        }
      }
    }
    else {

      for (edge_id = 0; edge_id < ms->n_edges; edge_id++) {
        ii = face_cel_p[edge_id][0];
        jj = face_cel_p[edge_id][1];
        for (kk = 0; kk < db_size; kk++) {
          dd[ii*db_size + kk] -= fabs(xa[2*edge_id]);
          dd[jj*db_size + kk] -= fabs(xa[2*edge_id + 1]);
        }
      }

    }

  }

  _b_diag_dom_diag_normalize(mc->d_val, dd, ms->n_rows, db_size);
}

/*----------------------------------------------------------------------------
 * Measure Diagonal dominance of native block matrix.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *   dd     --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_bb_diag_dom_native(const cs_matrix_t  *matrix,
                    cs_real_t          *restrict dd)
{
  cs_lnum_t  ii, jj, kk, ll, edge_id;

  const auto ms =
    static_cast<const cs_matrix_struct_native_t *>(matrix->structure);
  const auto mc = static_cast<const cs_matrix_coeff_t *>(matrix->coeffs);

  const cs_real_t  *restrict xa = mc->e_val;
  const cs_lnum_t db_size = matrix->db_size;
  const cs_lnum_t eb_size = matrix->eb_size;
  const cs_lnum_t eb_size_2 = eb_size*eb_size;

  /* block diagonal contribution */

  _b_diag_dom_diag_contrib(mc->d_val, dd, ms->n_rows, ms->n_cols_ext, db_size);

  /* non-diagonal terms */

  if (mc->e_val != nullptr) {

    const cs_lnum_2_t *restrict face_cel_p = ms->edges;

    if (mc->symmetric) {

      for (edge_id = 0; edge_id < ms->n_edges; edge_id++) {
        ii = face_cel_p[edge_id][0];
        jj = face_cel_p[edge_id][1];
        for (kk = 0; kk < eb_size; kk++) {
          for (ll = 0; ll < eb_size; ll++) {
            cs_lnum_t si = edge_id*eb_size_2 + kk*eb_size + ll;
            dd[ii*db_size + kk] -= fabs(xa[si]);
            dd[jj*db_size + kk] -= fabs(xa[si]);
          }
        }
      }
    }
    else {

      for (edge_id = 0; edge_id < ms->n_edges; edge_id++) {
        ii = face_cel_p[edge_id][0];
        jj = face_cel_p[edge_id][1];
        for (kk = 0; kk < db_size; kk++) {
          for (ll = 0; ll < eb_size; ll++) {
            cs_lnum_t si0 = 2*edge_id*eb_size_2 + kk*eb_size + ll;
            cs_lnum_t si1 = (2*edge_id+1)*eb_size_2 + kk*eb_size + ll;
            dd[ii*db_size + kk] -= fabs(xa[si0]);
            dd[jj*db_size + kk] -= fabs(xa[si1]);
          }
        }
      }

    }

  }

  _b_diag_dom_diag_normalize(mc->d_val, dd, ms->n_rows, db_size);
}

/*----------------------------------------------------------------------------
 * Measure Diagonal dominance of CSR matrix.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *   dd     --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_diag_dom_csr(const cs_matrix_t  *matrix,
              cs_real_t          *restrict dd)
{
  cs_lnum_t  ii, jj, n_cols;
  double  sii, sign, d_val;
  const cs_lnum_t  *restrict col_id;
  const cs_real_t  *restrict m_row;

  const auto ms =
    static_cast<const cs_matrix_struct_csr_t *>(matrix->structure);
  const auto mc = static_cast<const cs_matrix_coeff_t *>(matrix->coeffs);
  cs_lnum_t  n_rows = ms->n_rows;

  /* Standard case */

# pragma omp parallel for private(jj, col_id, m_row, n_cols, sii, sign, d_val)
  for (ii = 0; ii < n_rows; ii++) {
    col_id = ms->col_id + ms->row_index[ii];
    m_row = mc->val + ms->row_index[ii];
    n_cols = ms->row_index[ii+1] - ms->row_index[ii];
    sii = 0.0;
    d_val = 0.0;
    for (jj = 0; jj < n_cols; jj++) {
      if (col_id[jj] == ii) {
        sign = 1.;
        d_val = fabs(m_row[jj]);
      }
      else
        sign = -1.;
      sii += sign * fabs(m_row[jj]);
    }
    if (d_val > 1.e-18)
      sii /= d_val;
    else if (sii > -1.e-18)
      sii = -1.e18;
    else
      sii = 0;
    dd[ii] = sii;
  }
}

/*----------------------------------------------------------------------------
 * Measure Diagonal dominance of MSR or distributed matrix.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *   dd     --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_diag_dom_dist(const cs_matrix_t  *matrix,
              cs_real_t          *restrict dd)
{
  const auto ms =
    static_cast<const cs_matrix_struct_dist_t *>(matrix->structure);
  const auto mc = static_cast<const cs_matrix_coeff_t *>(matrix->coeffs);

  cs_lnum_t  n_rows = ms->n_rows;

  /* diagonal contribution */

  _diag_dom_diag_contrib(mc->d_val, dd, ms->n_rows, ms->n_cols_ext);

  /* extra-diagonal contribution */

  if (mc->e_val != nullptr) {

    const cs_matrix_struct_csr_t  *ms_e = &(ms->e);

#   pragma omp parallel for
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
      const cs_real_t *restrict m_row = mc->e_val + ms_e->row_index[ii];
      cs_lnum_t n_cols = ms_e->row_index[ii+1] - ms_e->row_index[ii];
      cs_real_t sii = 0.0;
      for (cs_lnum_t jj = 0; jj < n_cols; jj++)
        sii -= fabs(m_row[jj]);
      dd[ii] += sii;
    }

  }

  if (mc->h_val != nullptr) {

    const cs_matrix_struct_csr_t  *ms_h = &(ms->h);

#   pragma omp parallel for
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
      const cs_real_t *restrict m_row = mc->h_val + ms_h->row_index[ii];
      cs_lnum_t n_cols = ms_h->row_index[ii+1] - ms_h->row_index[ii];
      cs_real_t sii = 0.0;
      for (cs_lnum_t jj = 0; jj < n_cols; jj++)
        sii -= fabs(m_row[jj]);
      dd[ii] += sii;
    }

  }

  _diag_dom_diag_normalize(mc->d_val, dd, n_rows);
}

/*----------------------------------------------------------------------------
 * Measure Diagonal dominance of MSR or distributed matrix, blocked version.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *   dd     --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_b_diag_dom_dist(const cs_matrix_t  *matrix,
                 cs_real_t          *restrict dd)
{
  const auto ms =
    static_cast<const cs_matrix_struct_dist_t *>(matrix->structure);
  const auto mc = static_cast<const cs_matrix_coeff_t *>(matrix->coeffs);

  const cs_lnum_t db_size = matrix->db_size;
  const cs_lnum_t  n_rows = ms->n_rows;

  /* diagonal contribution */

  _b_diag_dom_diag_contrib(mc->d_val, dd, ms->n_rows, ms->n_cols_ext, db_size);

  /* extra-diagonal contribution */

  if (mc->e_val != nullptr) {

    const cs_matrix_struct_csr_t  *ms_e = &(ms->e);

#   pragma omp parallel for
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
      const cs_real_t *restrict m_row = mc->e_val + ms_e->row_index[ii];
      cs_lnum_t n_cols = ms_e->row_index[ii+1] - ms_e->row_index[ii];
      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        for (cs_lnum_t kk = 0; kk < db_size; kk++)
          dd[ii*db_size + kk] -= fabs(m_row[jj]);
      }
    }

  }

  if (mc->h_val != nullptr) {

    const cs_matrix_struct_csr_t  *ms_h = &(ms->h);

#   pragma omp parallel for
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
      const cs_real_t *restrict m_row = mc->h_val + ms_h->row_index[ii];
      cs_lnum_t n_cols = ms_h->row_index[ii+1] - ms_h->row_index[ii];
      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        for (cs_lnum_t kk = 0; kk < db_size; kk++)
          dd[ii*db_size + kk] -= fabs(m_row[jj]);
      }
    }

  }

  _b_diag_dom_diag_normalize(mc->d_val, dd, ms->n_rows, db_size);
}

/*----------------------------------------------------------------------------
 * Diagonal contribution to matrix dump.
 *
 * parameters:
 *   da         <-- Pointer to coefficients array (usually matrix diagonal)
 *   m_coo      <-- Matrix coefficient coordinates array
 *   m_val      <-- Matrix coefficient values array
 *   g_coo_num  <-- Global coordinate numbers
 *   n_rows     <-- Number of rows
*----------------------------------------------------------------------------*/

static void
_pre_dump_diag_contrib(const cs_real_t  *restrict da,
                       cs_gnum_t        *restrict m_coo,
                       cs_real_t        *restrict m_val,
                       const cs_gnum_t  *restrict g_coo_num,
                       cs_lnum_t         n_rows)
{
  cs_lnum_t  ii;

  if (da != nullptr) {
#   pragma omp parallel for
    for (ii = 0; ii < n_rows; ii++) {
      m_coo[ii*2] = g_coo_num[ii];
      m_coo[ii*2 + 1] = g_coo_num[ii];
      m_val[ii] = da[ii];
    }
  }
  else {
#   pragma omp parallel for
    for (ii = 0; ii < n_rows; ii++) {
      m_coo[ii*2] = g_coo_num[ii];
      m_coo[ii*2 + 1] = g_coo_num[ii];
      m_val[ii] = 0.0;
    }
  }
}

/*----------------------------------------------------------------------------
 * Block diagonal contribution to matrix dump.
 *
 * parameters:
 *   da        <-- Pointer to coefficients array (usually matrix diagonal)
 *   m_coo     <-- Matrix coefficient coordinates array
 *   m_val     <-- Matrix coefficient values array
 *   g_coo_num <-- Global coordinate numbers
 *   n_rows    <-- Number of rows
 *   b_size    <-- block size
 *----------------------------------------------------------------------------*/

static void
_b_pre_dump_diag_contrib(const cs_real_t  *restrict da,
                         cs_gnum_t        *restrict m_coo,
                         cs_real_t        *restrict m_val,
                         const cs_gnum_t  *restrict g_coo_num,
                         cs_lnum_t         n_rows,
                         cs_lnum_t         b_size)
{
  cs_lnum_t  ii, jj, kk;
  cs_lnum_t  b_size_2 = b_size * b_size;

  if (da != nullptr) {
#   pragma omp parallel for private(jj, kk)
    for (ii = 0; ii < n_rows; ii++) {
      for (jj = 0; jj < b_size; jj++) {
        for (kk = 0; kk < b_size; kk++) {
          m_coo[(ii*b_size_2 + jj*b_size + kk)*2]
            = g_coo_num[ii]*b_size + jj;
          m_coo[(ii*b_size_2 + jj*b_size + kk)*2 + 1]
            = g_coo_num[ii]*b_size + kk;
          m_val[ii*b_size_2 + jj*b_size + kk]
            = da[ii*b_size_2 + jj*b_size + kk];
        }
      }
    }
  }
  else {
#   pragma omp parallel for private(jj, kk)
    for (ii = 0; ii < n_rows; ii++) {
      for (jj = 0; jj < b_size; jj++) {
        for (kk = 0; kk < b_size; kk++) {
          m_coo[(ii*b_size_2 + jj*b_size + kk)*2]
            = g_coo_num[ii]*b_size + jj;
          m_coo[(ii*b_size_2 + jj*b_size + kk)*2 + 1]
            = g_coo_num[ii]*b_size + kk;
          m_val[ii*b_size_2 + jj*b_size + kk] = 0.0;
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Prepare dump of native matrix.
 *
 * parameters:
 *   matrix    <-- Pointer to matrix structure
 *   g_coo_num <-- Global coordinate numbers
 *   m_coo     --> Matrix coefficient coordinates array
 *   m_val     --> Matrix coefficient values array
 *
 * returns:
 *   number of matrix entries
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_pre_dump_native(const cs_matrix_t   *matrix,
                 const cs_gnum_t     *g_coo_num,
                 cs_gnum_t          **m_coo,
                 cs_real_t          **m_val)
{
  cs_lnum_t  ii, jj, edge_id, dump_id;

  cs_gnum_t   *restrict _m_coo;
  cs_real_t   *restrict _m_val;

  const auto ms =
    static_cast<const cs_matrix_struct_native_t *>(matrix->structure);
  const auto mc = static_cast<const cs_matrix_coeff_t *>(matrix->coeffs);

  const cs_real_t  *restrict xa = mc->e_val;

  cs_lnum_t  n_entries = ms->n_rows + ms->n_edges*2;

  /* Allocate arrays */

  CS_MALLOC(_m_coo, n_entries*2, cs_gnum_t);
  CS_MALLOC(_m_val, n_entries, double);

  *m_coo = _m_coo;
  *m_val = _m_val;

  /* diagonal contribution */

  _pre_dump_diag_contrib(mc->d_val, _m_coo, _m_val, g_coo_num, ms->n_rows);

  /* non-diagonal terms */

  dump_id = ms->n_rows;

  if (mc->e_val != nullptr) {

    const cs_lnum_2_t *restrict face_cel_p = ms->edges;

    if (mc->symmetric) {

      for (edge_id = 0; edge_id < ms->n_edges; edge_id++) {
        ii = face_cel_p[edge_id][0];
        jj = face_cel_p[edge_id][1];
        _m_coo[dump_id*2] = g_coo_num[ii];
        _m_coo[dump_id*2 + 1] = g_coo_num[jj];
        _m_val[dump_id] = xa[edge_id];
        _m_coo[dump_id*2 + 2] = g_coo_num[jj];
        _m_coo[dump_id*2 + 3] = g_coo_num[ii];
        _m_val[dump_id + 1] = xa[edge_id];
        dump_id += 2;
      }

    }
    else {

      for (edge_id = 0; edge_id < ms->n_edges; edge_id++) {
        ii = face_cel_p[edge_id][0];
        jj = face_cel_p[edge_id][1];
        _m_coo[dump_id*2] = g_coo_num[ii];
        _m_coo[dump_id*2 + 1] = g_coo_num[jj];
        _m_val[dump_id] = xa[2*edge_id];
        _m_coo[dump_id*2 + 2] = g_coo_num[jj];
        _m_coo[dump_id*2 + 3] = g_coo_num[ii];
        _m_val[dump_id + 1] = xa[2*edge_id + 1];
        dump_id += 2;
      }

    }

  }

  return n_entries;
}

/*----------------------------------------------------------------------------
 * Prepare dump of native matrix with block diagonal.
 *
 * parameters:
 *   matrix    <-- Pointer to matrix structure
 *   g_coo_num <-- Global coordinate numbers
 *   m_coo     --> Matrix coefficient coordinates array
 *   m_val     --> Matrix coefficient values array
 *
 * returns:
 *   number of matrix entries
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_b_pre_dump_native(const cs_matrix_t  *matrix,
                   const cs_gnum_t     *g_coo_num,
                   cs_gnum_t          **m_coo,
                   cs_real_t          **m_val)
{
  cs_lnum_t  ii, jj, kk, edge_id, dump_id;

  cs_gnum_t   *restrict _m_coo;
  cs_real_t   *restrict _m_val;

  const auto ms =
    static_cast<const cs_matrix_struct_native_t *>(matrix->structure);
  const auto mc = static_cast<const cs_matrix_coeff_t *>(matrix->coeffs);

  const cs_real_t  *restrict xa = mc->e_val;
  const cs_lnum_t db_size = matrix->db_size;

  cs_lnum_t  n_entries = (ms->n_rows*db_size + ms->n_edges*2) * db_size;

  /* Allocate arrays */

  CS_MALLOC(_m_coo, n_entries*2, cs_gnum_t);
  CS_MALLOC(_m_val, n_entries, double);

  *m_coo = _m_coo;
  *m_val = _m_val;

  /* block diagonal contribution */

  _b_pre_dump_diag_contrib(mc->d_val, _m_coo, _m_val,
                           g_coo_num, ms->n_rows, db_size);

  /* non-diagonal terms */

  dump_id = ms->n_rows*db_size*db_size;

  if (mc->e_val != nullptr) {

    const cs_lnum_2_t *restrict face_cel_p = ms->edges;

    if (mc->symmetric) {

      for (edge_id = 0; edge_id < ms->n_edges; edge_id++) {
        ii = face_cel_p[edge_id][0];
        jj = face_cel_p[edge_id][1];
        for (kk = 0; kk < db_size; kk++) {
          _m_coo[dump_id*2] = g_coo_num[ii]*db_size + kk;
          _m_coo[dump_id*2 + 1] = g_coo_num[jj]*db_size + kk;
          _m_val[dump_id] = xa[edge_id];
          _m_coo[dump_id*2 + 2] = g_coo_num[jj]*db_size + kk;
          _m_coo[dump_id*2 + 3] = g_coo_num[ii]*db_size + kk;
          _m_val[dump_id + 1] = xa[edge_id];
          dump_id += 2;
        }
      }
    }
    else {

      for (edge_id = 0; edge_id < ms->n_edges; edge_id++) {
        ii = face_cel_p[edge_id][0];
        jj = face_cel_p[edge_id][1];
        for (kk = 0; kk < db_size; kk++) {
          _m_coo[dump_id*2] = g_coo_num[ii]*db_size + kk;
          _m_coo[dump_id*2 + 1] = g_coo_num[jj]*db_size + kk;
          _m_val[dump_id] = xa[edge_id*2];
          _m_coo[dump_id*2 + 2] = g_coo_num[jj]*db_size + kk;
          _m_coo[dump_id*2 + 3] = g_coo_num[ii]*db_size + kk;
          _m_val[dump_id + 1] = xa[edge_id*2 + 1];
          dump_id += 2;
        }
      }

    }

  }

  return n_entries;
}

/*----------------------------------------------------------------------------
 * Prepare dump of CSR matrix.
 *
 * parameters:
 *   matrix    <-- Pointer to matrix structure
 *   g_coo_num <-- Global coordinate numbers
 *   m_coo     --> Matrix coefficient coordinates array
 *   m_val     --> Matrix coefficient values array
 *
 * returns:
 *   number of matrix entries
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_pre_dump_csr(const cs_matrix_t   *matrix,
              const cs_gnum_t     *g_coo_num,
              cs_gnum_t          **m_coo,
              cs_real_t          **m_val)
{
  cs_lnum_t  ii, jj, n_cols, dump_id;
  const cs_lnum_t  *restrict col_id;
  const cs_real_t  *restrict m_row;

  cs_gnum_t   *restrict _m_coo;
  cs_real_t   *restrict _m_val;

  const auto ms =
    static_cast<const cs_matrix_struct_csr_t *>(matrix->structure);
  const auto mc = static_cast<const cs_matrix_coeff_t *>(matrix->coeffs);
  cs_lnum_t  dump_id_shift = 0;
  cs_lnum_t  n_rows = ms->n_rows;

  cs_lnum_t  n_entries = ms->row_index[n_rows];

  /* Allocate arrays */

  CS_MALLOC(_m_coo, n_entries*2, cs_gnum_t);
  CS_MALLOC(_m_val, n_entries, double);

  *m_coo = _m_coo;
  *m_val = _m_val;

  /* Contribution */

  if (ms->have_diag == false) {
    assert(0);
#   pragma omp parallel for
    for (ii = 0; ii < n_rows; ii++) {
      _m_coo[ii*2] = g_coo_num[ii];
      _m_coo[ii*2+1] = g_coo_num[ii];
      _m_val[ii] = 0.0;
    }
    dump_id_shift = n_rows;
  }

# pragma omp parallel for private(jj, dump_id, col_id, m_row, n_cols)
  for (ii = 0; ii < n_rows; ii++) {
    col_id = ms->col_id + ms->row_index[ii];
    m_row = mc->val + ms->row_index[ii];
    n_cols = ms->row_index[ii+1] - ms->row_index[ii];
    for (jj = 0; jj < n_cols; jj++) {
      dump_id = ms->row_index[ii] + jj + dump_id_shift;
      _m_coo[dump_id*2] = g_coo_num[ii];
      _m_coo[dump_id*2+1] = g_coo_num[col_id[jj]];
      _m_val[dump_id] = m_row[jj];
    }
  }

  return n_entries;
}

/*----------------------------------------------------------------------------
 * Prepare dump of MSR matrix.
 *
 * parameters:
 *   matrix    <-- Pointer to matrix structure
 *   g_coo_num <-- Global coordinate numbers
 *   m_coo     --> Matrix coefficient coordinates array
 *   m_val     --> Matrix coefficient values array
 *
 * returns:
 *   number of matrix entries
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_pre_dump_msr(const cs_matrix_t   *matrix,
              const cs_gnum_t     *g_coo_num,
              cs_gnum_t          **m_coo,
              cs_real_t          **m_val)
{
  cs_lnum_t  ii, jj, n_cols, dump_id;
  const cs_lnum_t  *restrict col_id;
  const cs_real_t  *restrict m_row;

  cs_gnum_t   *restrict _m_coo;
  cs_real_t   *restrict _m_val;

  const auto ms =
    static_cast<const cs_matrix_struct_dist_t *>(matrix->structure);
  const auto ms_e = static_cast<const cs_matrix_struct_csr_t *>(&(ms->e));
  const auto mc = static_cast<const cs_matrix_coeff_t *>(matrix->coeffs);

  const cs_lnum_t  n_rows = ms->n_rows;
  cs_lnum_t  n_entries = ms_e->row_index[n_rows] + ms_e->n_rows;

  /* Allocate arrays */

  CS_MALLOC(_m_coo, n_entries*2, cs_gnum_t);
  CS_MALLOC(_m_val, n_entries, double);

  *m_coo = _m_coo;
  *m_val = _m_val;

  /* diagonal contribution */

  _pre_dump_diag_contrib(mc->d_val, _m_coo, _m_val, g_coo_num, ms->n_rows);

  /* extra-diagonal contribution */

  if (mc->e_val != nullptr) {
#   pragma omp parallel for private(jj, dump_id, col_id, m_row, n_cols)
    for (ii = 0; ii < n_rows; ii++) {
      col_id = ms_e->col_id + ms_e->row_index[ii];
      m_row = mc->e_val + ms_e->row_index[ii];
      n_cols = ms_e->row_index[ii+1] - ms_e->row_index[ii];
      for (jj = 0; jj < n_cols; jj++) {
        dump_id = ms_e->row_index[ii] + jj + ms_e->n_rows;
        _m_coo[dump_id*2] = g_coo_num[ii];
        _m_coo[dump_id*2+1] = g_coo_num[col_id[jj]];
        _m_val[dump_id] = m_row[jj];
      }
    }
  }
  else {
#   pragma omp parallel for private(jj, dump_id, col_id, n_cols)
    for (ii = 0; ii < n_rows; ii++) {
      col_id = ms_e->col_id + ms_e->row_index[ii];
      n_cols = ms_e->row_index[ii+1] - ms_e->row_index[ii];
      for (jj = 0; jj < n_cols; jj++) {
        dump_id = ms_e->row_index[ii] + jj + ms_e->n_rows;
        _m_coo[dump_id*2] = g_coo_num[ii];
        _m_coo[dump_id*2+1] = g_coo_num[col_id[jj]];
        _m_val[dump_id] = 0.0;
      }
    }
  }

  return n_entries;
}

/*----------------------------------------------------------------------------
 * Prepare dump of MSR matrix, blocked version.
 *
 * parameters:
 *   matrix    <-- Pointer to matrix structure
 *   g_coo_num <-- Global coordinate numbers
 *   m_coo     --> Matrix coefficient coordinates array
 *   m_val     --> Matrix coefficient values array
 *
 * returns:
 *   number of matrix entries
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_b_pre_dump_msr(const cs_matrix_t   *matrix,
                const cs_gnum_t     *g_coo_num,
                cs_gnum_t          **m_coo,
                cs_real_t          **m_val)
{
  cs_lnum_t  ii, jj, kk, n_cols, dump_id;
  const cs_lnum_t  *restrict col_id;
  const cs_real_t  *restrict m_row;

  cs_gnum_t   *restrict _m_coo;
  cs_real_t   *restrict _m_val;

  const auto ms =
    static_cast<const cs_matrix_struct_dist_t *>(matrix->structure);
  const auto       ms_e = static_cast<const cs_matrix_struct_csr_t *>(&(ms->e));
  const auto mc = static_cast<const cs_matrix_coeff_t *>(matrix->coeffs);
  const cs_lnum_t  db_size = matrix->db_size;
  const cs_lnum_t  n_rows = ms->n_rows;
  const cs_lnum_t  dump_id_shift = ms->n_rows*db_size*db_size;

  cs_lnum_t  n_entries =   ms_e->row_index[n_rows]*db_size
                         + ms_e->n_rows*db_size*db_size;

  /* Allocate arrays */

  CS_MALLOC(_m_coo, n_entries*2, cs_gnum_t);
  CS_MALLOC(_m_val, n_entries, double);

  *m_coo = _m_coo;
  *m_val = _m_val;

  /* diagonal contribution */

  _b_pre_dump_diag_contrib(mc->d_val, _m_coo, _m_val,
                           g_coo_num, ms->n_rows, db_size);


  /* extra-diagonal contribution */

  if (mc->e_val != nullptr) {
#   pragma omp parallel for private(jj, kk, dump_id, col_id, m_row, n_cols)
    for (ii = 0; ii < n_rows; ii++) {
      col_id = ms_e->col_id + ms_e->row_index[ii];
      m_row = mc->e_val + ms_e->row_index[ii];
      n_cols = ms_e->row_index[ii+1] - ms_e->row_index[ii];
      for (jj = 0; jj < n_cols; jj++) {
        for (kk = 0; kk < db_size; kk++) {
          dump_id = (ms_e->row_index[ii] + jj)*db_size + kk + dump_id_shift;
          _m_coo[dump_id*2] = g_coo_num[ii]*db_size + kk;
          _m_coo[dump_id*2+1] = g_coo_num[col_id[jj]]*db_size + kk;
          _m_val[dump_id] = m_row[jj];
        }
      }
    }
  }
  else {
#   pragma omp parallel for private(jj, kk, dump_id, col_id, n_cols)
    for (ii = 0; ii < n_rows; ii++) {
      col_id = ms_e->col_id + ms_e->row_index[ii];
      n_cols = ms_e->row_index[ii+1] - ms_e->row_index[ii];
      for (jj = 0; jj < n_cols; jj++) {
        for (kk = 0; kk < db_size; kk++) {
          dump_id = (ms_e->row_index[ii] + jj)*db_size + kk + dump_id_shift;
          _m_coo[dump_id*2] = g_coo_num[ii]*db_size + kk;
          _m_coo[dump_id*2+1] = g_coo_num[col_id[jj]]*db_size + kk;
          _m_val[dump_id] = 0.0;
        }
      }
    }
  }

  return n_entries;
}

/*----------------------------------------------------------------------------
 * Write header for dump of matrix to native file.
 *
 * 3 characters are used to define (in order):
 * - the size of cs_gnum_t (int, long, or long long)
 * - the size of floating point values (double, should be size _)
 * - 'l' for little-endian files, 'b' for big-endian files.
 *
 * parameters:
 *   f <-- Pointer to file structure
 *----------------------------------------------------------------------------*/

static void
_write_header_simple(cs_file_t  *f)
{
  unsigned char flags[3] = {(unsigned char) sizeof(cs_gnum_t),
                            (unsigned char) sizeof(double),
                            'b'};

  /* Check if system is "big-endian" or "little-endian" */
  {
    unsigned  int_endian = 0;
    *((char *)(&int_endian)) = '\1';
    if (int_endian == 1)
      flags[2] = 'l';
  }

  assert(sizeof(cs_gnum_t) < 255);
  assert(sizeof(double) < 255); /* should always be 8 */

  cs_file_write_global(f, flags, 1, 3);
}

/*----------------------------------------------------------------------------
 * Locally sort matrix dump data.
 *
 * parameters:
 *   n_entries <-- Number of local matrix entries
 *   m_coords  <-> Matrix coefficient coordinates (interlaced)
 *   m_vals    <-> Matrix coefficient values
 *----------------------------------------------------------------------------*/

static void
_sort_matrix_dump_data(cs_lnum_t   n_entries,
                       cs_gnum_t  *m_coords,
                       double     *m_vals)
{
  cs_lnum_t ii, jj;
  cs_gnum_t *_m_coords;
  double *_m_vals;

  cs_lnum_t *order = cs_order_gnum_s(nullptr, m_coords, 2, n_entries);

  CS_MALLOC(_m_coords, n_entries*2, cs_gnum_t);

  for (ii = 0; ii < n_entries; ii++) {
    jj = order[ii];
    _m_coords[ii*2] = m_coords[jj*2];
    _m_coords[ii*2+1] = m_coords[jj*2+1];
  }
  memcpy(m_coords, _m_coords, n_entries*2*sizeof(cs_gnum_t));

  CS_FREE(_m_coords);

  CS_MALLOC(_m_vals, n_entries, double);

  for (ii = 0; ii < n_entries; ii++)
    _m_vals[ii] = m_vals[order[ii]];

  memcpy(m_vals, _m_vals, n_entries*sizeof(double));

  CS_FREE(_m_vals);

  CS_FREE(order);
}

/*----------------------------------------------------------------------------
 * Prepare data for matrix dump.
 *
 * Build local arrays with global matrix coordinates and coefficients.
 *
 * parameters:
 *   m         <-- Pointer to matrix structure
 *   n_entries --> Number of local matrix entries
 *   m_coords  --> Matrix coefficient coordinates (interlaced)
 *   m_vals    --> Matrix coefficient values
 *----------------------------------------------------------------------------*/

static void
_prepare_matrix_dump_data(const cs_matrix_t   *m,
                          cs_lnum_t           *n_entries,
                          cs_gnum_t          **m_coords,
                          double             **m_vals)
{
  cs_lnum_t ii, jj;
  cs_lnum_t _n_entries = 0;
  cs_gnum_t coo_shift = 1, n_g_rows = 0;

  cs_gnum_t *g_coo_num = nullptr;
  cs_gnum_t *_m_coords = nullptr;
  double *_m_vals = nullptr;

  CS_MALLOC(g_coo_num, m->n_cols_ext, cs_gnum_t);

  n_g_rows = m->n_rows;

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    cs_gnum_t loc_shift = m->n_rows;
    MPI_Scan(&loc_shift, &coo_shift, 1, CS_MPI_GNUM, MPI_SUM, cs_glob_mpi_comm);
    coo_shift = coo_shift + 1 - loc_shift;
  }
#endif

  for (ii = 0; ii < m->n_rows; ii++)
    g_coo_num[ii] = ii + coo_shift;

  if (m->halo != nullptr)
    cs_halo_sync_untyped(m->halo,
                         CS_HALO_STANDARD,
                         sizeof(cs_gnum_t),
                         g_coo_num);

  /* Populate arrays based on matrix type */

  switch(m->type) {
  case CS_MATRIX_NATIVE:
    if (m->db_size == 1)
      _n_entries = _pre_dump_native(m, g_coo_num, &_m_coords, &_m_vals);
    else
      _n_entries = _b_pre_dump_native(m, g_coo_num, &_m_coords, &_m_vals);
    break;
  case CS_MATRIX_CSR:
    assert(m->db_size == 1);
    _n_entries = _pre_dump_csr(m, g_coo_num, &_m_coords, &_m_vals);
    break;
  case CS_MATRIX_MSR:
    if (m->db_size == 1)
      _n_entries = _pre_dump_msr(m, g_coo_num, &_m_coords, &_m_vals);
    else
      _n_entries = _b_pre_dump_msr(m, g_coo_num, &_m_coords, &_m_vals);
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Dump of matrixes in %s format\n"
                "is not operational yet."),
              _(cs_matrix_get_type_name(m)));
    break;
  }

  CS_FREE(g_coo_num);

  /* Sort values */

  _sort_matrix_dump_data(_n_entries, _m_coords, _m_vals);

  /* Remove references to extra periodicity ghost values if present */

  for (ii = 0, jj = 0; ii < _n_entries; ii++) {
    if (_m_coords[ii*2] <= n_g_rows && _m_coords[ii*2+1] <= n_g_rows) {
      _m_coords[jj*2] = _m_coords[ii*2];
      _m_coords[jj*2+1] = _m_coords[ii*2+1];
      _m_vals[jj] = _m_vals[ii];
      jj += 1;
    }
  }
  _n_entries = jj;

  /* Return arrays */

  *n_entries = _n_entries;
  *m_coords = _m_coords;
  *m_vals = _m_vals;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write matrix data to file in parallel mode.
 *
 * Arrays written are freed by this function.
 *
 * The full matrix is written, whether symmetric or not. Row and column
 * coordinates are 1-to-n based.
 *
 * File format:
 *   - number of matrix entries (1 integer of type cs_gnum_t).
 *   - array of row coordinates (n_entries integers of type cs_gnum_t)
 *   - array of column coordinates (n_entries integers of type cs_gnum_t)
 *   - array of matrix coefficients (n_entries values of type double)
 *
 * parameters:
 *   m <-- Pointer to matrix structure
 *   f <-> associated file pointer
 *----------------------------------------------------------------------------*/

static void
_write_matrix_g(const cs_matrix_t  *m,
                cs_file_t          *f)
{
  int  block_rank_step = 1, min_block_size = 0;
  cs_lnum_t  block_size = 0;
  cs_gnum_t  n_glob_ents = 0;

  cs_gnum_t  *b_coords = nullptr, *c_coords = nullptr, *r_coords = nullptr;
  double  *b_vals = nullptr;

  cs_block_dist_info_t bi;

  fvm_io_num_t *io_num = nullptr;
  cs_part_to_block_t *d = nullptr;

  cs_lnum_t  n_entries = 0;
  cs_gnum_t  *m_coords = nullptr;
  double  *m_vals = nullptr;

  const cs_datatype_t gnum_type
    = (sizeof(cs_gnum_t) == 8) ? CS_UINT64 : CS_UINT32;

  /* Initialization for matrix coefficients */

  _prepare_matrix_dump_data(m, &n_entries, &m_coords, &m_vals);

  /* Initialization for redistribution */

  io_num = fvm_io_num_create_from_adj_s(nullptr, m_coords, n_entries, 2);

  n_glob_ents = fvm_io_num_get_global_count(io_num);

  cs_file_get_default_comm(&block_rank_step, nullptr, nullptr);

  bi = cs_block_dist_compute_sizes(cs_glob_rank_id,
                                   cs_glob_n_ranks,
                                   block_rank_step,
                                   min_block_size/2,
                                   n_glob_ents);

  d = cs_part_to_block_create_by_gnum(cs_glob_mpi_comm,
                                      bi,
                                      n_entries,
                                      fvm_io_num_get_global_num(io_num));

  /* Write number of entries */

  cs_file_write_global(f, &n_glob_ents, sizeof(cs_gnum_t), 1);

  /* Distribute to blocks */

  block_size = (bi.gnum_range[1] - bi.gnum_range[0]);

  /* Distribute coordinate blocks on ranks */

  if (block_size > 0)
    CS_MALLOC(b_coords, block_size*2, cs_gnum_t);

  cs_part_to_block_copy_array(d, gnum_type, 2, m_coords, b_coords);

  CS_FREE(m_coords);

  /* De-interlace coordinates */

  if (block_size > 0) {
   cs_lnum_t  ii;
   CS_MALLOC(r_coords, block_size, cs_gnum_t);
   CS_MALLOC(c_coords, block_size, cs_gnum_t);
   for (ii = 0; ii < block_size; ii++) {
     r_coords[ii] = b_coords[ii*2];
     c_coords[ii] = b_coords[ii*2+1];
   }
   CS_FREE(b_coords);
  }

  /* Write coordinate blocks */

  cs_file_write_block_buffer(f,
                             r_coords,
                             sizeof(cs_gnum_t),
                             1,
                             bi.gnum_range[0],
                             bi.gnum_range[1]);

  cs_file_write_block_buffer(f,
                             c_coords,
                             sizeof(cs_gnum_t),
                             1,
                             bi.gnum_range[0],
                             bi.gnum_range[1]);

  CS_FREE(c_coords);
  CS_FREE(r_coords);

  /* Distribute value blocks on ranks */

  if (block_size > 0)
    CS_MALLOC(b_vals, block_size, double);

  cs_part_to_block_copy_array(d, CS_DOUBLE, 1, m_vals, b_vals);

  CS_FREE(m_vals);

  /* Write value blocks */

  cs_file_write_block_buffer(f,
                             b_vals,
                             sizeof(double),
                             1,
                             bi.gnum_range[0],
                             bi.gnum_range[1]);

  CS_FREE(b_vals);

  /* Free matrix coefficient distribution structures */

  cs_part_to_block_destroy(&d);
  io_num = fvm_io_num_destroy(io_num);
}

#endif /* #if defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Write matrix data to file in serial mode.
 *
 * Arrays written are freed by this function.
 *
 * The full matrix is written, whether symmetric or not. Row and column
 * coordinates are 1-to-n based.
 *
 * File format:
 *   - number of matrix entries (1 integer of type cs_gnum_t).
 *   - array of row coordinates (n_entries integers of type cs_gnum_t)
 *   - array of column coordinates (n_entries integers of type cs_gnum_t)
 *   - array of matrix coefficients (n_entries values of type double)
 *
 * parameters:
 *   m <-- Pointer to matrix structure
 *   f <-> associated file pointer
 *----------------------------------------------------------------------------*/

static void
_write_matrix_l(const cs_matrix_t  *m,
                cs_file_t          *f)
{
  cs_lnum_t  ii;
  cs_gnum_t  n_glob_ents = 0;

  cs_gnum_t  *c_coords = nullptr, *r_coords = nullptr;

  cs_lnum_t  n_entries = 0;
  cs_gnum_t  *m_coords = nullptr;
  double  *m_vals = nullptr;

  /* Initialization */

  _prepare_matrix_dump_data(m, &n_entries, &m_coords, &m_vals);

  n_glob_ents = n_entries;

  /* Write number of entries */

  cs_file_write_global(f, &n_glob_ents, sizeof(cs_gnum_t), 1);

  /* De-interlace coordinates */

  CS_MALLOC(r_coords, n_entries, cs_gnum_t);
  CS_MALLOC(c_coords, n_entries, cs_gnum_t);
  for (ii = 0; ii < n_entries; ii++) {
    r_coords[ii] = m_coords[ii*2];
    c_coords[ii] = m_coords[ii*2+1];
  }
  CS_FREE(m_coords);

  /* Write coordinate blocks */

  cs_file_write_global(f, r_coords, sizeof(cs_gnum_t), n_entries);
  cs_file_write_global(f, c_coords, sizeof(cs_gnum_t), n_entries);

  CS_FREE(r_coords);
  CS_FREE(c_coords);

  /* Write value blocks */

  cs_file_write_global(f, m_vals, sizeof(double), n_entries);

  CS_FREE(m_vals);
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Write vector data to file in parallel mode.
 *
 * File format:
 *   - vector size matrix entries (1 integer of type cs_gnum_t).
 *   - array of vector values (n_elts values of type double)
 *
 * parameters:
 *   n_elts <-- local number of matrix entries
 *   vals   <-- vector element values
 *   f      <-> associated file pointer
 *----------------------------------------------------------------------------*/

static void
_write_vector_g(cs_lnum_t         n_elts,
                const cs_real_t  *vals,
                cs_file_t        *f)
{
  cs_lnum_t  ii;

  int        block_rank_step = 1, min_block_size = 0;
  cs_lnum_t  block_size = 0;
  cs_gnum_t  coo_shift = 1;
  cs_gnum_t  local_max = 0, n_glob_ents = 0;

  cs_gnum_t  *g_elt_num = nullptr;
  double  *b_vals = nullptr;

  cs_block_dist_info_t bi;

  cs_part_to_block_t *d = nullptr;

  cs_gnum_t loc_shift = n_elts;
  MPI_Scan(&loc_shift, &coo_shift, 1, CS_MPI_GNUM, MPI_SUM, cs_glob_mpi_comm);
  coo_shift = coo_shift - loc_shift;

  /* Get maximum global number value */

  if (n_elts > 0)
    local_max = n_elts + coo_shift;

  MPI_Allreduce(&local_max, &n_glob_ents, 1, CS_MPI_GNUM, MPI_MAX,
                cs_glob_mpi_comm);

  CS_MALLOC(g_elt_num, n_elts, cs_gnum_t);

  for (ii = 0; ii < n_elts; ii++)
    g_elt_num[ii] = ii + coo_shift + 1;

  cs_file_get_default_comm(&block_rank_step, nullptr, nullptr);

  /* Redistribution structures */

  bi = cs_block_dist_compute_sizes(cs_glob_rank_id,
                                   cs_glob_n_ranks,
                                   block_rank_step,
                                   min_block_size/2,
                                   n_glob_ents);

  d = cs_part_to_block_create_by_gnum(cs_glob_mpi_comm,
                                      bi,
                                      n_elts,
                                      g_elt_num);

  /* Write number of entries */

  cs_file_write_global(f, &n_glob_ents, sizeof(cs_gnum_t), 1);

  /* Distribute to blocks */

  block_size = (bi.gnum_range[1] - bi.gnum_range[0]);

  if (block_size > 0)
    CS_MALLOC(b_vals, block_size, double);

  if (sizeof(cs_real_t) != sizeof(double)) {
    double  *p_vals = nullptr;
    CS_MALLOC(p_vals, n_elts, double);
    for (ii = 0; ii < n_elts; ii++)
      p_vals[ii] = vals[ii];
    cs_part_to_block_copy_array(d, CS_DOUBLE, 1, p_vals, b_vals);
    CS_FREE(p_vals);
  }
  else
    cs_part_to_block_copy_array(d, CS_DOUBLE, 1, vals, b_vals);

  /* Write value blocks */

  cs_file_write_block_buffer(f,
                             b_vals,
                             sizeof(double),
                             1,
                             bi.gnum_range[0],
                             bi.gnum_range[1]);

  CS_FREE(b_vals);

  /* Free distribution structures */

  cs_part_to_block_destroy(&d);

  CS_FREE(g_elt_num);
}

#endif /* #if defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Write vector data to file in serial mode.
 *
 * File format:
 *   - vector size matrix entries (1 integer of type cs_gnum_t).
 *   - array of vector values (n_elts values of type double)
 *
 * parameters:
 *   n_elts <-- local number of matrix entries
 *   vals   <-- vector element values
 *   f      <-> associated file pointer
 *----------------------------------------------------------------------------*/

static void
_write_vector_l(cs_lnum_t         n_elts,
                const cs_real_t  *vals,
                cs_file_t        *f)
{
  cs_gnum_t  n_glob_ents = n_elts;

  /* Write number of entries */

  cs_file_write_global(f, &n_glob_ents, sizeof(cs_gnum_t), 1);

  /* Write values */

  if (sizeof(cs_real_t) != sizeof(double)) {
    cs_lnum_t ii;
    double  *p_vals = nullptr;
    CS_MALLOC(p_vals, n_elts, double);
    for (ii = 0; ii < n_elts; ii++)
      p_vals[ii] = vals[ii];
    cs_file_write_global(f, p_vals, sizeof(double), n_elts);
    CS_FREE(p_vals);
  }
  else
    cs_file_write_global(f, vals, sizeof(double), n_elts);
}

/*----------------------------------------------------------------------------
 * Compute the Frobenius norm of a matrix.
 *
 * parameters:
 *   m <-- pointer to matrix structure
 *
 * returns:
 *    matrix Frobenius norm (or -1 if not computable)
 *----------------------------------------------------------------------------*/

static double
_frobenius_norm(const cs_matrix_t  *m)
{
  double retval = -1.;

  if (m == nullptr)
    return retval;

  cs_matrix_fill_type_t ft = m->fill_type;

  switch(m->type) {

  case CS_MATRIX_NATIVE:
    {
      cs_lnum_t  d_stride = m->db_size * m->db_size;
      cs_lnum_t  e_stride = m->eb_size * m->eb_size;
      const auto ms =
        static_cast<const cs_matrix_struct_native_t *>(m->structure);
      const auto mc = static_cast<const cs_matrix_coeff_t *>(m->coeffs);
      double e_mult = (m->eb_size == 1) ? m->db_size : 1;
      if (mc->symmetric)
        e_mult *= 2;
      else
        e_stride *= 2;

      retval = cs_dot_xx(d_stride*m->n_rows, mc->d_val);

      double ed_contrib = 0.;
      const cs_real_t  *restrict xa = mc->e_val;
#     pragma omp parallel reduction(+:ed_contrib) if (ms->n_edges > CS_THR_MIN)
      {
        double c = 0; /* Use Kahan compensated summation for
                         sum of block contributions (but not for local
                         block sums, for simplicity and performance) */
#       pragma omp for
        for (cs_lnum_t edge_id = 0; edge_id < ms->n_edges; edge_id++) {
          cs_lnum_t ii = ms->edges[edge_id][0];
          if (ii < ms->n_rows) {
            double bsum = 0;
            for (cs_lnum_t kk = 0; kk < e_stride; kk++) {
              bsum  += (  xa[edge_id*e_stride + kk]
                        * xa[edge_id*e_stride + kk]);
            }
            double z = bsum - c;
            double t = ed_contrib + z;
            c = (t - ed_contrib) - z;
            ed_contrib = t;
          }
        }
      }
      retval += ed_contrib*e_mult;

      cs_parall_sum(1, CS_DOUBLE, &retval);
    }
    break;

  case CS_MATRIX_CSR:
    {
      assert(   ft == CS_MATRIX_SCALAR
             || ft == CS_MATRIX_SCALAR_SYM
             || ft == CS_MATRIX_BLOCK);
      cs_lnum_t  stride = m->eb_size * m->eb_size;
      const auto ms = static_cast<const cs_matrix_struct_csr_t *>(m->structure);
      const auto *mc = static_cast<const cs_matrix_coeff_t *>(m->coeffs);
      cs_lnum_t n_vals = ms->row_index[m->n_rows];
      retval = cs_dot_xx(stride*n_vals, mc->val);
      cs_parall_sum(1, CS_DOUBLE, &retval);
    }
    break;

  case CS_MATRIX_MSR:
    {
      cs_lnum_t  d_stride = m->db_size * m->db_size;
      cs_lnum_t  e_stride = m->eb_size * m->eb_size;
      const auto ms =
        static_cast<const cs_matrix_struct_dist_t *>(m->structure);
      const auto ms_e   = static_cast<const cs_matrix_struct_csr_t *>(&(ms->e));
      const auto mc = static_cast<const cs_matrix_coeff_t  *>(m->coeffs);
      cs_lnum_t n_vals = ms_e->row_index[m->n_rows];
      double d_mult = (m->eb_size == 1) ? m->db_size : 1;
      retval = cs_dot_xx(d_stride*m->n_rows, mc->d_val);
      retval += d_mult * cs_dot_xx(e_stride*n_vals, mc->e_val);
      cs_parall_sum(1, CS_DOUBLE, &retval);
    }
    break;

  case CS_MATRIX_DIST:
    {
      cs_lnum_t  d_stride = m->db_size * m->db_size;
      cs_lnum_t  e_stride = m->eb_size * m->eb_size;
      const auto ms =
        static_cast<const cs_matrix_struct_dist_t *>(m->structure);
      const auto ms_e = static_cast<const cs_matrix_struct_csr_t *>(&(ms->e));
      const auto ms_h = static_cast<const cs_matrix_struct_csr_t *>(&(ms->h));
      const auto mc = static_cast<const cs_matrix_coeff_t  *>(m->coeffs);
      cs_lnum_t n_e_vals = ms_e->row_index[m->n_rows];
      cs_lnum_t n_h_vals = ms_h->row_index[m->n_rows];
      double d_mult = (m->eb_size == 1) ? m->db_size : 1;
      retval = cs_dot_xx(d_stride*m->n_rows, mc->d_val);
      retval += d_mult * cs_dot_xx(e_stride*n_e_vals, mc->e_val);
      retval += d_mult * cs_dot_xx(e_stride*n_h_vals, mc->h_val);
      cs_parall_sum(1, CS_DOUBLE, &retval);
    }
    break;

    default:
      retval = -1;
  }

  if (retval > 0)
    retval = sqrt(retval);

  return retval;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute diagonal dominance metric.
 *
 * parameters:
 *   matrix <-- pointer to matrix structure
 *   dd     --> diagonal dominance (normalized)
 *----------------------------------------------------------------------------*/

void
cs_matrix_diag_dominance(const cs_matrix_t  *matrix,
                         cs_real_t           dd[])
{
  const cs_halo_t *halo = matrix->halo;

  switch(matrix->type) {
  case CS_MATRIX_NATIVE:
    if (matrix->db_size == 1)
      _diag_dom_native(matrix, dd);
    else if (matrix->eb_size == 1)
      _b_diag_dom_native(matrix, dd);
    else
      _bb_diag_dom_native(matrix, dd);
    break;
  case CS_MATRIX_CSR:
    assert(matrix->db_size == 1);
    _diag_dom_csr(matrix, dd);
    break;
  case CS_MATRIX_MSR:
    if (matrix->db_size == 1)
      _diag_dom_dist(matrix, dd);
    else
      _b_diag_dom_dist(matrix, dd);
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Extraction of diagonal dominance of matrixes in %s format\n"
                "is not operational yet."),
              _(cs_matrix_get_type_name(matrix)));
    break;
  }

  /* Sync ghost rows as a precaution */

  if (halo != nullptr) {
    if (matrix->db_size == 1)
      cs_halo_sync_var(halo, CS_HALO_STANDARD, dd);
    else {
      cs_halo_sync_var_strided(halo, CS_HALO_STANDARD, dd, matrix->db_size);
      if (halo->n_transforms > 0 && matrix->db_size == 3)
        cs_halo_perio_sync_var_vect(halo,
                                    CS_HALO_STANDARD,
                                    dd,
                                    matrix->db_size);
    }
  }
}

/*----------------------------------------------------------------------------
 * Dump a linear system matrix and right-hand side to file.
 *
 * parameters:
 *   matrix <-- pointer to matrix structure
 *   rhs    <-- right hand side vector
 *   name   <-- identifier string used in file name
 *----------------------------------------------------------------------------*/

void
cs_matrix_dump_linear_system(const cs_matrix_t  *matrix,
                             const cs_real_t     rhs[],
                             const char         *name)
{
  char filename[64];
  cs_gnum_t n_g_rows = matrix->n_rows;

  cs_file_t  *f = nullptr;

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    cs_gnum_t n_l_rows = matrix->n_rows;
    MPI_Allreduce(&n_l_rows, &n_g_rows, 1, CS_MPI_GNUM, MPI_SUM,
                  cs_glob_mpi_comm);
  }
#endif

  snprintf(filename, 63, "%s_%010llu", name, (unsigned long long)n_g_rows);
  filename[63] = '\0';

  f = cs_file_open_default(filename, CS_FILE_MODE_WRITE);

  _write_header_simple(f);

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    _write_matrix_g(matrix, f);
    _write_vector_g(matrix->n_rows, rhs, f);
  }
#endif
  if (cs_glob_n_ranks == 1) {
    _write_matrix_l(matrix, f);
    _write_vector_l(matrix->n_rows, rhs, f);
  }

  f = cs_file_free(f);
}

/*----------------------------------------------------------------------------
 * Dump a vector to file.
 *
 * parameters:
 *   n_rows <-- number of local rows
 *   stride <-- associated stride
 *   x      <-- associated vector
 *   name   <-- identifier string used in file name
 *----------------------------------------------------------------------------*/

void
cs_matrix_dump_vector(const cs_lnum_t     n_rows,
                      const cs_lnum_t     stride,
                      const cs_real_t     x[],
                      const char         *name)
{
  char filename[64];
  cs_lnum_t n_l_rows = n_rows*stride;
  cs_gnum_t n_g_rows = n_l_rows;

  cs_file_t  *f = nullptr;

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    cs_gnum_t _n_g_rows = n_l_rows;
    MPI_Allreduce(&_n_g_rows, &n_g_rows, 1, CS_MPI_GNUM, MPI_SUM,
                  cs_glob_mpi_comm);
  }
#endif

  snprintf(filename, 63, "%s_%010llu", name, (unsigned long long)n_g_rows);
  filename[63] = '\0';

  f = cs_file_open_default(filename, CS_FILE_MODE_WRITE);

  _write_header_simple(f);

  /* Write the number of rows (may be useful for the treatment by an external
     tool) */

  cs_file_write_global(f, &n_g_rows, sizeof(cs_gnum_t), 1);

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    _write_vector_g(n_l_rows, x, f);
  }
#endif
  if (cs_glob_n_ranks == 1) {
    _write_vector_l(n_l_rows, x, f);
  }

  f = cs_file_free(f);
}

/*----------------------------------------------------------------------------
 * Dump a matrix to file.
 *
 * parameters:
 *   matrix <-- pointer to matrix structure
 *   name   <-- identifier string used in file name
 *----------------------------------------------------------------------------*/

void
cs_matrix_dump(const cs_matrix_t  *matrix,
               const char         *name)
{
  char filename[64];
  cs_gnum_t n_g_rows = matrix->n_rows;

  cs_file_t  *f = nullptr;

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    cs_gnum_t n_l_rows = matrix->n_rows;
    MPI_Allreduce(&n_l_rows, &n_g_rows, 1, CS_MPI_GNUM, MPI_SUM,
                  cs_glob_mpi_comm);
  }
#endif

  snprintf(filename, 63, "%s_%010llu", name, (unsigned long long)n_g_rows);
  filename[63] = '\0';

  f = cs_file_open_default(filename, CS_FILE_MODE_WRITE);

  _write_header_simple(f);

  /* Write the number of rows (may be useful for the treatment by an external
     tool) */

  cs_file_write_global(f, &n_g_rows, sizeof(cs_gnum_t), 1);


#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    _write_matrix_g(matrix, f);
  }
#endif
  if (cs_glob_n_ranks == 1) {
    _write_matrix_l(matrix, f);
  }

  f = cs_file_free(f);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log general info relative to matrix.
 *
 * \param[in]  matrix     pointer to matrix structure
 * \param[in]  verbosity  verbosity level
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_log_info(const cs_matrix_t  *matrix,
                   int                 verbosity)
{
  cs_log_t l = CS_LOG_DEFAULT;

  if (matrix == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              _("The matrix is not defined."));

  cs_log_printf(l,
                _("\n"
                  " Matrix info:\n"
                  "   type: %s\n"),
                cs_matrix_get_type_fullname(matrix));

  if (matrix->fill_type == CS_MATRIX_N_FILL_TYPES)
    return;

  cs_log_printf(l,
                _("   fill type: %s\n"),
                cs_matrix_fill_type_name[matrix->fill_type]);

  if (verbosity > 1) {
    double fnorm = _frobenius_norm(matrix);
    if (fnorm > -1)
      cs_log_printf(l,
                    _("   Frobenius norm: %15.9e\n"), fnorm);
  }

  cs_log_printf(l, "\n");
}

/*----------------------------------------------------------------------------
 * Test matrix dump operations.
 *
 * parameters:
 *   n_rows     <-- number of local rows
 *   n_cols_ext <-- number of columns including ghost columns (array size)
 *   n_edges    <-- local number of graph edges
 *   edges      <-- graph edges connectivity
 *   halo       <-- cell halo structure
 *   numbering  <-- vectorization or thread-related numbering info, or nullptr
 *----------------------------------------------------------------------------*/

void
cs_matrix_dump_test(cs_lnum_t              n_rows,
                    cs_lnum_t              n_cols_ext,
                    cs_lnum_t              n_edges,
                    const cs_lnum_2_t     *edges,
                    const cs_halo_t       *halo,
                    const cs_numbering_t  *numbering)
{
  cs_lnum_t  ii;
  int  test_id;

  cs_real_t  *da = nullptr, *xa = nullptr, *rhs = nullptr;
  cs_lnum_t  db_size = 3, db_size_2 = 3*3;
  cs_lnum_t  eb_size = 1;

  const int n_tests = 6;
  const char *name[] = {"matrix_native",
                        "matrix_native_sym",
                        "matrix_native_block",
                        "matrix_csr",
                        "matrix_msr",
                        "matrix_msr_block"};
  const cs_matrix_type_t type[] = {CS_MATRIX_NATIVE,
                                   CS_MATRIX_NATIVE,
                                   CS_MATRIX_NATIVE,
                                   CS_MATRIX_CSR,
                                   CS_MATRIX_MSR,
                                   CS_MATRIX_MSR};
  const bool sym_flag[] = {false, true, false, false, false, false};
  const int block_flag[] = {0, 0, 1, 0, 0, 1};

  /* Allocate and initialize  working arrays */
  /*-----------------------------------------*/

  CS_MALLOC(rhs, n_cols_ext*db_size, cs_real_t);

  CS_MALLOC(da, n_cols_ext*db_size_2, cs_real_t);
  CS_MALLOC(xa, n_edges*2, cs_real_t);

# pragma omp parallel for
  for (ii = 0; ii < n_cols_ext*db_size_2; ii++)
    da[ii] = 1.0 + ii*0.1/n_cols_ext;
# pragma omp parallel for
  for (ii = 0; ii < n_cols_ext*db_size; ii++)
    rhs[ii] = ii*0.1/n_cols_ext;

# pragma omp parallel for
  for (ii = 0; ii < n_edges; ii++) {
    xa[ii*2] = 0.5*(1.0 + ii*1.0/n_edges);
    xa[ii*2 + 1] = -0.5*(1.0 + ii*1.0/n_edges);
  }

  /* Loop on variant types */
  /*-----------------------*/

  for (test_id = 0; test_id < n_tests; test_id++) {

    cs_lnum_t _db_size = (block_flag[test_id]) ? db_size : 1;
    cs_lnum_t _eb_size = (block_flag[test_id]-1) ? eb_size : 1;

    cs_matrix_structure_t
      *ms = cs_matrix_structure_create(type[test_id],
                                       n_rows,
                                       n_cols_ext,
                                       n_edges,
                                       edges,
                                       halo,
                                       numbering);
    cs_matrix_t *m = cs_matrix_create(ms);

    cs_matrix_set_coefficients(m,
                               sym_flag[test_id],
                               _db_size,
                               _eb_size,
                               n_edges,
                               edges,
                               da,
                               xa);

    cs_matrix_dump_linear_system(m, rhs, name[test_id]);

    cs_matrix_release_coefficients(m);

    cs_matrix_destroy(&m);
    cs_matrix_structure_destroy(&ms);
  }

  CS_FREE(rhs);

  CS_FREE(da);
  CS_FREE(xa);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
