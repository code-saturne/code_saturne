/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2011 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

/*============================================================================
 * Utilitary functions for sparse matrixes.
 *============================================================================*/

#include "cs_defs.h"

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
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_blas.h"
#include "cs_halo.h"
#include "cs_log.h"
#include "cs_numbering.h"
#include "cs_prototypes.h"
#include "cs_perio.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_matrix_util.h"

#include "cs_matrix_priv.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

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
 * y[i] = abs(da[i]), with da possibly NULL.
 *
 * parameters:
 *   da          <-- Pointer to coefficients array (usually matrix diagonal)
 *   dd          --> Resulting vector
 *   n_cells     <-- Number of cells
 *   n_cells_ext <-- Number of cells + ghost cells
*----------------------------------------------------------------------------*/

static void
_diag_dom_diag_contrib(const cs_real_t  *restrict da,
                       cs_real_t        *restrict dd,
                       cs_lnum_t         n_cells,
                       cs_lnum_t         n_cells_ext)
{
  cs_lnum_t  ii;

  if (da != NULL) {
#   pragma omp parallel for
    for (ii = 0; ii < n_cells; ii++)
      dd[ii] = fabs(da[ii]);
    for (ii = n_cells; ii < n_cells_ext; ii++)
      dd[ii] = 0.0;
  }
  else {
#   pragma omp parallel for
    for (ii = 0; ii < n_cells_ext; ii++)
      dd[ii] = 0.0;
  }

}

/*----------------------------------------------------------------------------
 * Block diagonal contribution to diagonal dominance, with da possibly NULL.
 *
 * parameters:
 *   da          <-- Pointer to coefficients array (usually matrix diagonal)
 *   dd          --> Resulting vector
 *   n_cells     <-- Number of cells
 *   n_cells_ext <-- Number of cells + ghost cells
 *   b_size      <-- block size, including padding:
 *                   b_size[0]: useful block size
 *                   b_size[1]: vector block extents
 *                   b_size[2]: matrix line extents
 *                   b_size[3]: matrix line*column (block) extents
 *----------------------------------------------------------------------------*/

static void
_b_diag_dom_diag_contrib(const cs_real_t  *restrict da,
                         cs_real_t        *restrict dd,
                         cs_lnum_t         n_cells,
                         cs_lnum_t         n_cells_ext,
                         const int         b_size[4])
{
  cs_lnum_t  ii, jj, kk;
  double  sign;
  const cs_lnum_t  dd_size = n_cells_ext*b_size[1];

# pragma omp parallel for
  for (ii = 0; ii < dd_size; ii++)
    dd[ii] = 0.0;

  if (da != NULL) {
#   pragma omp parallel for private(jj, kk, sign)
    for (ii = 0; ii < n_cells; ii++) {
      for (jj = 0; jj < b_size[1]; jj++)
        dd[ii*b_size[1] + jj] = 0.0;
      for (jj = 0; jj < b_size[0]; ii++) {
        for (kk = 0; kk < b_size[0]; kk++) {
          sign = (jj == kk) ? 1. : -1.;
          dd[ii*b_size[1] + kk]
            += sign*fabs(da[ii*b_size[3] + jj*b_size[2] + kk]);
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Normalize diagonal dominance, with da possibly NULL.
 *
 * parameters:
 *   da          <-- Pointer to coefficients array (usually matrix diagonal)
 *   dd          --> Resulting vector
 *   n_cells     <-- Number of cells
*----------------------------------------------------------------------------*/

static void
_diag_dom_diag_normalize(const cs_real_t  *restrict da,
                         cs_real_t        *restrict dd,
                         cs_lnum_t         n_cells)
{
  cs_lnum_t  ii;

  if (da != NULL) {
#   pragma omp parallel for
    for (ii = 0; ii < n_cells; ii++) {
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
    for (ii = 0; ii < n_cells; ii++) {
      if (dd[ii] > -1.e-18)
        dd[ii] = -1.e18;
      else
        dd[ii] = 0;
    }
  }
}

/*----------------------------------------------------------------------------
 * Block diagonal contribution to diagonal dominance, with da possibly NULL.
 *
 * parameters:
 *   da      <-- Pointer to coefficients array (usually matrix diagonal)
 *   dd      --> Resulting vector
 *   n_cells <-- Number of cells
 *   b_size  <-- block size, including padding:
 *               b_size[0]: useful block size
 *               b_size[1]: vector block extents
 *               b_size[2]: matrix line extents
 *               b_size[3]: matrix line*column (block) extents
 *----------------------------------------------------------------------------*/

static void
_b_diag_dom_diag_normalize(const cs_real_t  *restrict da,
                           cs_real_t        *restrict dd,
                           cs_lnum_t         n_cells,
                           const int         b_size[4])
{
  cs_lnum_t  ii, jj;
  double  d_val;

  if (da != NULL) {
#   pragma omp parallel for private(jj, d_val)
    for (ii = 0; ii < n_cells; ii++) {
      for (jj = 0; jj < b_size[0]; ii++) {
        d_val = fabs(da[ii*b_size[3] + jj*b_size[2] + jj]);
        if (d_val > 1.e-18)
          dd[ii*b_size[1] + jj] /= d_val;
        else if (dd[ii*b_size[1] + jj] > -1.e-18)
          dd[ii*b_size[1] + jj] = -1.e18;
        else
          dd[ii*b_size[1] + jj] = 0;
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
  cs_lnum_t  ii, jj, face_id;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;

  const cs_real_t  *restrict xa = mc->xa;

  /* diagonal contribution */

  _diag_dom_diag_contrib(mc->da, dd, ms->n_cells, ms->n_cells_ext);

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    if (mc->symmetric) {

      const cs_lnum_t *restrict face_cel_p = ms->face_cell;

      for (face_id = 0; face_id < ms->n_faces; face_id++) {
        ii = face_cel_p[2*face_id] -1;
        jj = face_cel_p[2*face_id + 1] -1;
        dd[ii] -= fabs(xa[face_id]);
        dd[jj] -= fabs(xa[face_id]);
      }

    }
    else {

      const cs_lnum_t *restrict face_cel_p = ms->face_cell;

      for (face_id = 0; face_id < ms->n_faces; face_id++) {
        ii = face_cel_p[2*face_id] -1;
        jj = face_cel_p[2*face_id + 1] -1;
        dd[ii] -= fabs(xa[2*face_id]);
        dd[jj] -= fabs(xa[2*face_id + 1]);
      }

    }

  }

  _diag_dom_diag_normalize(mc->da, dd, ms->n_cells);
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
  cs_lnum_t  ii, jj, kk, face_id;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;

  const cs_real_t  *restrict xa = mc->xa;
  const int *b_size = matrix->b_size;

  /* block diagonal contribution */

  _b_diag_dom_diag_contrib(mc->da, dd, ms->n_cells, ms->n_cells_ext, b_size);

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    if (mc->symmetric) {

      const cs_lnum_t *restrict face_cel_p = ms->face_cell;

      for (face_id = 0; face_id < ms->n_faces; face_id++) {
        ii = face_cel_p[2*face_id] -1;
        jj = face_cel_p[2*face_id + 1] -1;
        for (kk = 0; kk < b_size[0]; kk++) {
          dd[ii*b_size[1] + kk] -= fabs(xa[face_id]);
          dd[jj*b_size[1] + kk] -= fabs(xa[face_id]);
        }
      }
    }
    else {

      const cs_lnum_t *restrict face_cel_p = ms->face_cell;

      for (face_id = 0; face_id < ms->n_faces; face_id++) {
        ii = face_cel_p[2*face_id] -1;
        jj = face_cel_p[2*face_id + 1] -1;
        for (kk = 0; kk < b_size[0]; kk++) {
          dd[ii*b_size[1] + kk] -= fabs(xa[2*face_id]);
          dd[jj*b_size[1] + kk] -= fabs(xa[2*face_id + 1]);
        }
      }

    }

  }

  _b_diag_dom_diag_normalize(mc->da, dd, ms->n_cells, b_size);
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
  cs_lnum_t  *restrict col_id;
  cs_real_t  *restrict m_row;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_csr_t  *mc = matrix->coeffs;
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
 * Measure Diagonal dominance of symmetric CSR matrix.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *   dd     --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_diag_dom_csr_sym(const cs_matrix_t  *matrix,
                  cs_real_t          *restrict dd)
{
  cs_lnum_t  ii, jj, n_cols;
  cs_lnum_t  *restrict col_id;
  cs_real_t  *restrict m_row;

  const cs_matrix_struct_csr_sym_t  *ms = matrix->structure;
  const cs_matrix_coeff_csr_sym_t  *mc = matrix->coeffs;
  cs_lnum_t  n_rows = ms->n_rows;

  cs_lnum_t sym_jj_start = 0;

  /* By construction, the matrix has either a full or an empty
     diagonal structure, so testing this on the first row is enough */

  if (ms->col_id[ms->row_index[0]] == 0)
    sym_jj_start = 1;

  /* Initialize dd */

  for (ii = 0; ii < ms->n_cols; ii++)
    dd[ii] = 0.0;

  /* Upper triangular + diagonal part in case of symmetric structure */

  for (ii = 0; ii < n_rows; ii++) {
    cs_real_t  sii = 0.0;
    col_id = ms->col_id + ms->row_index[ii];
    m_row = mc->val + ms->row_index[ii];
    n_cols = ms->row_index[ii+1] - ms->row_index[ii];
    for (jj = 0; jj < n_cols; jj++) {
      double sign = (col_id[jj] == ii) ? 1. : -1.;
      sii += sign * fabs(m_row[jj]);
    }
    dd[ii] += sii;
    for (jj = sym_jj_start; jj < n_cols; jj++)
      dd[col_id[jj]] -= fabs(m_row[jj]);
  }

  /* normalize */

  for (ii = 0; ii < n_rows; ii++) {
    cs_real_t sii = 0.0, d_val = 0.0;
    col_id = ms->col_id + ms->row_index[ii];
    m_row = mc->val + ms->row_index[ii];
    n_cols = ms->row_index[ii+1] - ms->row_index[ii];
    for (jj = 0; jj < n_cols; jj++) {
      if (col_id[jj] == ii)
        d_val = fabs(m_row[jj]);
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
 * Measure Diagonal dominance of MSR matrix.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *   dd     --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_diag_dom_msr(const cs_matrix_t  *matrix,
              cs_real_t          *restrict dd)
{
  cs_lnum_t  ii, jj, n_cols;
  cs_real_t  sii;
  cs_lnum_t  *restrict col_id;
  cs_real_t  *restrict m_row;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_msr_t  *mc = matrix->coeffs;
  cs_lnum_t  n_rows = ms->n_rows;

  /* diagonal contribution */

  _diag_dom_diag_contrib(mc->d_val, dd, ms->n_rows, ms->n_cols);

  /* extra-diagonal contribution */

  if (mc->x_val != NULL) {

#   pragma omp parallel for private(jj, col_id, m_row, n_cols, sii)
    for (ii = 0; ii < n_rows; ii++) {
      col_id = ms->col_id + ms->row_index[ii];
      m_row = mc->x_val + ms->row_index[ii];
      n_cols = ms->row_index[ii+1] - ms->row_index[ii];
      sii = 0.0;
      for (jj = 0; jj < n_cols; jj++)
        sii -= fabs(m_row[jj]);
      dd[ii] += sii;
    }

  }

  _diag_dom_diag_normalize(mc->d_val, dd, n_rows);
}

/*----------------------------------------------------------------------------
 * Measure Diagonal dominance of MSR matrix, blocked version.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *   dd     --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_b_diag_dom_msr(const cs_matrix_t  *matrix,
                cs_real_t          *restrict dd)
{
  cs_lnum_t  ii, jj, kk, n_cols;
  cs_lnum_t  *restrict col_id;
  cs_real_t  *restrict m_row;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_msr_t  *mc = matrix->coeffs;
  const int *b_size = matrix->b_size;
  const cs_lnum_t  n_rows = ms->n_rows;

  /* diagonal contribution */

  _b_diag_dom_diag_contrib(mc->d_val, dd, ms->n_rows, ms->n_cols, b_size);

  /* extra-diagonal contribution */

  if (mc->x_val != NULL) {

#   pragma omp parallel for private(jj, kk, col_id, m_row, n_cols, sii)
    for (ii = 0; ii < n_rows; ii++) {
      col_id = ms->col_id + ms->row_index[ii];
      m_row = mc->x_val + ms->row_index[ii];
      n_cols = ms->row_index[ii+1] - ms->row_index[ii];
      for (jj = 0; jj < n_cols; jj++) {
        for (kk = 0; kk < b_size[0]; kk++)
          dd[ii*b_size[1] + kk] -= fabs(m_row[jj]);
      }
    }

  }

  _b_diag_dom_diag_normalize(mc->d_val, dd, ms->n_rows, b_size);
}

/*----------------------------------------------------------------------------
 * Measure Diagonal dominance of symmetric MSR matrix.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *   dd     --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_diag_dom_msr_sym(const cs_matrix_t  *matrix,
                  cs_real_t          *restrict dd)
{
  cs_lnum_t  ii, jj, n_cols;
  cs_lnum_t  *restrict col_id;
  cs_real_t  *restrict m_row;

  const cs_matrix_struct_csr_sym_t  *ms = matrix->structure;
  const cs_matrix_coeff_msr_sym_t  *mc = matrix->coeffs;
  cs_lnum_t  n_rows = ms->n_rows;

  /* diagonal contribution */

  _diag_dom_diag_contrib(mc->d_val, dd, ms->n_rows, ms->n_cols);

  /* Initialize dd */

  for (ii = 0; ii < ms->n_cols; ii++)
    dd[ii] = 0.0;

  /* Upper triangular + diagonal part in case of symmetric structure */

  for (ii = 0; ii < n_rows; ii++) {
    cs_real_t  sii = 0.0;
    col_id = ms->col_id + ms->row_index[ii];
    m_row = mc->x_val + ms->row_index[ii];
    n_cols = ms->row_index[ii+1] - ms->row_index[ii];
    for (jj = 0; jj < n_cols; jj++) {
      double sign = (col_id[jj] == ii) ? 1. : -1.;
      sii += sign * fabs(m_row[jj]);
    }
    dd[ii] += sii;
    for (jj = 0; jj < n_cols; jj++)
      dd[col_id[jj]] -= fabs(m_row[jj]);
  }

}

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
  const cs_int_t n_cells = cs_glob_mesh->n_cells;
  const cs_halo_t *halo = cs_glob_mesh->halo;
  const cs_matrix_structure_t *ms = matrix->structure;

  switch(matrix->type) {
  case CS_MATRIX_NATIVE:
    if (matrix->b_size[3] == 1)
      _diag_dom_native(matrix, dd);
    else
      _b_diag_dom_native(matrix, dd);
    break;
  case CS_MATRIX_CSR:
    assert(matrix->b_size[3] == 1);
    _diag_dom_csr(matrix, dd);
    break;
  case CS_MATRIX_CSR_SYM:
    assert(matrix->b_size[3] == 1);
    _diag_dom_csr_sym(matrix, dd);
    break;
  case CS_MATRIX_MSR:
    if (matrix->b_size[3] == 1)
      _diag_dom_msr(matrix, dd);
    else
      _b_diag_dom_msr(matrix, dd);
    break;
    break;
  case CS_MATRIX_MSR_SYM:
    assert(matrix->b_size[3] == 1);
    _diag_dom_msr_sym(matrix, dd);
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Extraction of diagonal dominance of matrixes in %s format\n"
                "is not operational yet."),
              _(cs_matrix_type_name[matrix->type]));
    break;
  }

  /* Sync ghost cells as a precaution */

  if (halo != NULL) {
    if (matrix->b_size[3] == 1)
      cs_halo_sync_var(halo, CS_HALO_STANDARD, dd);
    else {
      cs_halo_sync_var_strided(halo, CS_HALO_STANDARD, dd, matrix->b_size[1]);
      if (halo->n_transforms > 0 && matrix->b_size[0] == 3)
        cs_perio_sync_var_vect(halo, CS_HALO_STANDARD, dd, matrix->b_size[1]);
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
