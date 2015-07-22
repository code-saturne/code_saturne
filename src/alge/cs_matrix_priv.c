/*============================================================================
 * Sparse Matrix Representation and Operations.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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

#if defined (HAVE_MKL)
#include <mkl_spblas.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_blas.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_numbering.h"
#include "cs_prototypes.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_matrix_priv.h"

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
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * \brief Destroy a CSR matrix structure.
 *
 * \param[inout]  matrix  pointer to a CSR matrix structure pointer
 *----------------------------------------------------------------------------*/

void
cs_matrix_destroy_struct_csr(cs_matrix_struct_csr_t  **matrix)
{
  if (matrix != NULL && *matrix !=NULL) {

    cs_matrix_struct_csr_t  *ms = *matrix;

    if (ms->row_index != NULL)
      BFT_FREE(ms->row_index);

    if (ms->col_id != NULL)
      BFT_FREE(ms->col_id);

    BFT_FREE(ms);

    *matrix = NULL;

  }
}

/*----------------------------------------------------------------------------
 * Create CSR matrix coefficients.
 *
 * returns:
 *   pointer to allocated CSR coefficients structure.
 *----------------------------------------------------------------------------*/

cs_matrix_coeff_csr_t *
cs_matrix_create_coeff_csr(void)
{
  cs_matrix_coeff_csr_t  *mc = NULL;

  /* Allocate */

  BFT_MALLOC(mc, 1, cs_matrix_coeff_csr_t);

  /* Initialize */

  mc->n_prefetch_rows = 0;

  mc->val = NULL;

  mc->x_prefetch = NULL;

  mc->d_val = NULL;
  mc->_d_val = NULL;

  return mc;
}

/*----------------------------------------------------------------------------
 * Release shared CSR matrix coefficients.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *----------------------------------------------------------------------------*/

void
cs_matrix_release_coeffs_csr(cs_matrix_t  *matrix)
{
  cs_matrix_coeff_csr_t  *mc = matrix->coeffs;
  if (mc != NULL)
    mc->d_val = NULL;
  return;
}

/*----------------------------------------------------------------------------
 * Copy diagonal of CSR matrix.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *   da     --> Diagonal (pre-allocated, size: n_rows)
 *----------------------------------------------------------------------------*/

void
cs_matrix_copy_diagonal_csr(const cs_matrix_t  *matrix,
                            cs_real_t          *restrict da)
{
  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_csr_t  *mc = matrix->coeffs;
  cs_lnum_t  n_rows = ms->n_rows;

  if (ms->have_diag == true) {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      const cs_lnum_t  *restrict col_id = ms->col_id + ms->row_index[ii];
      const cs_real_t  *restrict m_row = mc->val + ms->row_index[ii];
      cs_lnum_t  n_cols = ms->row_index[ii+1] - ms->row_index[ii];

      da[ii] = 0.0;
      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        if (col_id[jj] == ii) {
          da[ii] = m_row[jj];
          break;
        }
      }

    }

  }
  else { /* if (have_diag == false) */
#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      da[ii] = 0.0;
  }

}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with CSR matrix.
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- Pointer to matrix structure
 *   x            <-- Multipliying vector values
 *   y            --> Resulting vector
 *----------------------------------------------------------------------------*/

void
cs_matrix_vec_p_l_csr(bool                exclude_diag,
                      const cs_matrix_t  *matrix,
                      const cs_real_t    *restrict x,
                      cs_real_t          *restrict y)
{
  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_csr_t  *mc = matrix->coeffs;
  cs_lnum_t  n_rows = ms->n_rows;

  /* Standard case */

  if (!exclude_diag) {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      cs_lnum_t *restrict col_id = ms->col_id + ms->row_index[ii];
      cs_real_t *restrict m_row = mc->val + ms->row_index[ii];
      cs_lnum_t n_cols = ms->row_index[ii+1] - ms->row_index[ii];
      cs_real_t sii = 0.0;

      for (cs_lnum_t jj = 0; jj < n_cols; jj++)
        sii += (m_row[jj]*x[col_id[jj]]);

      y[ii] = sii;

    }

  }

  /* Exclude diagonal */

  else {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      cs_lnum_t *restrict col_id = ms->col_id + ms->row_index[ii];
      cs_real_t *restrict m_row = mc->val + ms->row_index[ii];
      cs_lnum_t n_cols = ms->row_index[ii+1] - ms->row_index[ii];
      cs_real_t sii = 0.0;

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        if (col_id[jj] != ii)
          sii += (m_row[jj]*x[col_id[jj]]);
      }

      y[ii] = sii;

    }
  }

}

#if defined (HAVE_MKL)

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with MSR matrix, using MKL
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- Pointer to matrix structure
 *   x            <-- Multipliying vector values
 *   y            --> Resulting vector
 *----------------------------------------------------------------------------*/

void
cs_matrix_vec_p_l_csr_mkl(bool                exclude_diag,
                          const cs_matrix_t  *matrix,
                          const cs_real_t    *restrict x,
                          cs_real_t          *restrict y)
{
  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_csr_t  *mc = matrix->coeffs;

  int n_rows = ms->n_rows;
  char transa[] = "n";

  if (exclude_diag)
    bft_error(__FILE__, __LINE__, 0,
              _(_no_exclude_diag_error_str), __func__);

  mkl_cspblas_dcsrgemv(transa,
                       &n_rows,
                       mc->val,
                       ms->row_index,
                       ms->col_id,
                       (double *)x,
                       y);
}

#endif /* defined (HAVE_MKL) */

/*----------------------------------------------------------------------------
 * Copy diagonal of native or MSR matrix.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *   da     --> Diagonal (pre-allocated, size: n_cells)
 *----------------------------------------------------------------------------*/

void
cs_matrix_copy_diagonal_separate(const cs_matrix_t  *matrix,
                                 cs_real_t          *restrict da)
{
  const cs_real_t *_da = NULL;
  if (matrix->type == CS_MATRIX_NATIVE) {
    const cs_matrix_coeff_native_t  *mc = matrix->coeffs;
    _da = mc->da;
  }
  else if (matrix->type == CS_MATRIX_MSR) {
    const cs_matrix_coeff_msr_t  *mc = matrix->coeffs;
    _da = mc->d_val;
  }
  const cs_lnum_t  n_cells = matrix->n_cells;

  /* Unblocked version */

  if (matrix->db_size[3] == 1) {

    if (_da != NULL) {
#     pragma omp parallel for  if(n_cells > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_cells; ii++)
        da[ii] = _da[ii];
    }
    else {
#     pragma omp parallel for  if(n_cells > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_cells; ii++)
        da[ii] = 0.0;
    }

  }

  /* Blocked version */

  else {

    const int *db_size = matrix->db_size;

    if (_da != NULL) {
#     pragma omp parallel for  if(n_cells*db_size[0] > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_cells; ii++) {
        for (cs_lnum_t jj = 0; jj < db_size[0]; jj++)
          da[ii*db_size[1] + jj] = _da[ii*db_size[3] + jj*db_size[2] + jj];
      }
    }
    else {
#     pragma omp parallel for  if(n_cells*db_size[1] > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_cells*db_size[1]; ii++)
        da[ii] = 0.0;
    }
  }
}

/*----------------------------------------------------------------------------
 * Create MSR matrix coefficients.
 *
 * returns:
 *   pointer to allocated MSR coefficients structure.
 *----------------------------------------------------------------------------*/

cs_matrix_coeff_msr_t *
cs_matrix_create_coeff_msr(void)
{
  cs_matrix_coeff_msr_t  *mc;

  /* Allocate */

  BFT_MALLOC(mc, 1, cs_matrix_coeff_msr_t);

  /* Initialize */

  mc->n_prefetch_rows = 0;
  mc->max_db_size = 0;
  mc->max_eb_size = 0;

  mc->d_val = NULL;

  mc->_d_val = NULL;
  mc->x_val = NULL;

  mc->x_prefetch = NULL;

  return mc;
}

/*----------------------------------------------------------------------------
 * Release shared MSR matrix coefficients.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *----------------------------------------------------------------------------*/

void
cs_matrix_release_coeffs_msr(cs_matrix_t  *matrix)
{
  cs_matrix_coeff_msr_t  *mc = matrix->coeffs;
  if (mc !=NULL) {
    /* Unmap shared values */
    mc->d_val = NULL;
  }
}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with MSR matrix.
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- Pointer to matrix structure
 *   x            <-- Multipliying vector values
 *   y            --> Resulting vector
 *----------------------------------------------------------------------------*/

void
cs_matrix_vec_p_l_msr(bool                exclude_diag,
                      const cs_matrix_t  *matrix,
                      const cs_real_t    *restrict x,
                      cs_real_t          *restrict y)
{
  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_msr_t  *mc = matrix->coeffs;
  cs_lnum_t  n_rows = ms->n_rows;

  /* Standard case */

  if (!exclude_diag && mc->d_val != NULL) {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      cs_lnum_t *restrict col_id = ms->col_id + ms->row_index[ii];
      cs_real_t *restrict m_row = mc->x_val + ms->row_index[ii];
      cs_lnum_t n_cols = ms->row_index[ii+1] - ms->row_index[ii];
      cs_real_t sii = 0.0;

      for (cs_lnum_t jj = 0; jj < n_cols; jj++)
        sii += (m_row[jj]*x[col_id[jj]]);

      y[ii] = sii + mc->d_val[ii]*x[ii];

    }

  }

  /* Exclude diagonal */

  else {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      cs_lnum_t *restrict col_id = ms->col_id + ms->row_index[ii];
      cs_real_t *restrict m_row = mc->x_val + ms->row_index[ii];
      cs_lnum_t n_cols = ms->row_index[ii+1] - ms->row_index[ii];
      cs_real_t sii = 0.0;

      for (cs_lnum_t jj = 0; jj < n_cols; jj++)
        sii += (m_row[jj]*x[col_id[jj]]);

      y[ii] = sii;

    }
  }

}

#if defined (HAVE_MKL)

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with MSR matrix, using MKL
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- Pointer to matrix structure
 *   x            <-- Multipliying vector values
 *   y            --> Resulting vector
 *----------------------------------------------------------------------------*/

void
cs_matrix_vec_p_l_msr_mkl(bool                exclude_diag,
                          const cs_matrix_t  *matrix,
                          const cs_real_t    *restrict x,
                          cs_real_t          *restrict y)
{
  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_msr_t  *mc = matrix->coeffs;

  int n_rows = ms->n_rows;
  char transa[] = "n";

  mkl_cspblas_dcsrgemv(transa,
                       &n_rows,
                       mc->x_val,
                       ms->row_index,
                       ms->col_id,
                       (double *)x,
                       y);

  /* Add diagonal contribution */

  if (!exclude_diag && mc->d_val != NULL) {
    cs_lnum_t ii;
    const double *restrict da = mc->d_val;
#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (ii = 0; ii < n_rows; ii++)
      y[ii] += da[ii] * x[ii];
  }
}

#endif /* defined (HAVE_MKL) */
