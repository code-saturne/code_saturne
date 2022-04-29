/*============================================================================
 * Sparse Matrix Representation and Operations using PETSc library.
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
 * PETSc headers
 *----------------------------------------------------------------------------*/

/* Avoid warnings due to previous values */
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_URL
#undef PACKAGE_VERSION

#include <petscmat.h>
#include <petscvec.h>
#include <petscversion.h>
#include <petscviewer.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_log.h"
#include "cs_numbering.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_matrix.h"
#include "cs_base_accel.h"
#include "cs_matrix_default.h"
#include "cs_matrix_petsc.h"
#include "cs_matrix_petsc_priv.h"
#include "cs_matrix_priv.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*! \file cs_matrix_petsc.c
 *
 * \brief Sparse Matrix Representation and Operations using PETSc.
 *
 * Setting a matrix type to PETc directly allows avoidind duplicating values,
 * and assigning them directly to PETSc.
 *
 * This should save memory, thoght performance might not be optimal,
 * as the current implementation requires many calls to MatSetValues.
 * MatSetValues can require at least one call per matrix block, as the
 * values arrays passed to it assume a dense rectangular block.
 * Block structures are not exploited yet.
 *
 * Performance could possibly be improved by assembling local arrays in a
 * manner similar to code_saturne's CS_MATRIX_DIST matrix, but using
 * PETSc types (to ensure sizes are the same), and assigning it to PETSc
 * using MatCreateMPIAIJWithSplitArrays. Use of this function is not
 * recommended by PETSc as it is considered cumbersome and inflexible.
 * This would not be too much of an issue here (since we already have that
 * capacity), but might not be compatible with all PETSc matrix formats.
 * An alternative would be to use MatCreateMPIAIJWithArrays with a temporary
 * copy, but would involve higher peak memory use. Finally, the
 * MatSetValuesBatch function mighr be used in the future, if it is extended
 * so as to allow multiple calls for a given matrix.
 */

/*----------------------------------------------------------------------------*/

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

static const char _petsc_ij_type_name[] = "PETSc, MATAIJ";
static const char _petsc_ij_type_fullname[] = "PETSc (MATAIJ series)";

static char _init_status = 0; /* 0 at start,  1 if initialized, 2 if finalized */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with PETSc matrix
 *
 * Note that since this function creates vectors, it may have significant
 * overhead. Since it is present for completednes, this should not be an issue.
 * If used more often, vectors could be maintained in the coeffs structure
 * along with the matrix.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_p_aij(const cs_matrix_t  *matrix,
               bool                exclude_diag,
               bool                sync,
               cs_real_t          *restrict x,
               cs_real_t          *restrict y)
{
  assert(exclude_diag == false);

  const PetscInt  n_rows = matrix->n_rows * matrix->db_size;
  const cs_matrix_coeffs_petsc_t  *coeffs = matrix->coeffs;

  /* Get pointers to structures through coefficients,
     and copy input values */

  PetscReal *_x = NULL;

  VecGetArray(coeffs->hx, &_x);
  for (PetscInt ii = 0; ii < n_rows; ii++) {
    _x[ii] = x[ii];;
  }
  VecRestoreArray(coeffs->hx, &_x);

  /* SpMv operation */

  MatMult(coeffs->hm, coeffs->hx, coeffs->hy);

  /* Copy data back */

  const PetscReal *_y = NULL;

  VecGetArrayRead(coeffs->hy, &_y);
  for (PetscInt ii = 0; ii < n_rows; ii++) {
    y[ii] = _y[ii];;
  }
  VecRestoreArrayRead(coeffs->hy, &_y);
}

/*----------------------------------------------------------------------------
 * Compute local and distant counts of matrix entries using an assembler
 *
 * The caller is responsible for freeing the returned diag_sizes and
 * offdiag_sizes arrays.
 *
 * parameters:
 *   ma            <-- pointer to matrix assembler structure
 *   diag_sizes    --> diagonal values
 *   offdiag_sizes --> off-diagonal values
 *----------------------------------------------------------------------------*/

static void
_compute_diag_sizes_assembler(const cs_matrix_assembler_t   *ma,
                              PetscInt                    **diag_sizes,
                              PetscInt                    **offdiag_sizes)
{
  const cs_lnum_t n_rows = cs_matrix_assembler_get_n_rows(ma);

  cs_lnum_t n_diag_add = (cs_matrix_assembler_get_separate_diag(ma)) ? 1 : 0;

  const cs_lnum_t *row_index = cs_matrix_assembler_get_row_index(ma);
  const cs_lnum_t *col_ids = cs_matrix_assembler_get_col_ids(ma);

  PetscInt *_diag_sizes, *_offdiag_sizes;
  BFT_MALLOC(_diag_sizes, n_rows, PetscInt);
  BFT_MALLOC(_offdiag_sizes, n_rows, PetscInt);

  /* Separate local and distant loops for better first touch logic */

# pragma omp parallel for  if(n_rows > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_rows; i++) {
    cs_lnum_t s_id = row_index[i];
    cs_lnum_t e_id = row_index[i+1];

    PetscInt n_r_diag = 0;
    PetscInt n_cols = e_id - s_id;

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      if (col_ids[j] < n_rows)
        n_r_diag++;
      else
        break;
    }

    _diag_sizes[i] = n_r_diag + n_diag_add;
    _offdiag_sizes[i] = n_cols - n_r_diag;
  }

  *diag_sizes = _diag_sizes;
  *offdiag_sizes = _offdiag_sizes;
}

/*----------------------------------------------------------------------------
 * Compute local and distant counts of matrix entries using an assembler
 * with full diagonal blocks and A.I extradiangonal blocks fill.
 *
 * The caller is responsible for freeing the returned diag_sizes and
 * offdiag_sizes arrays.
 *
 * parameters:
 *   ma            <-- pointer to matrix assembler structure
 *   db_size       <-- diagonal block size
 *   diag_sizes    --> diagonal values
 *   offdiag_sizes --> off-diagonal values
 *----------------------------------------------------------------------------*/

static void
_compute_diag_sizes_assembler_db(const cs_matrix_assembler_t   *ma,
                                 cs_lnum_t                     db_size,
                                 PetscInt                    **diag_sizes,
                                 PetscInt                    **offdiag_sizes)
{
  const cs_lnum_t n_rows = cs_matrix_assembler_get_n_rows(ma);

  cs_lnum_t n_diag_add
    = (cs_matrix_assembler_get_separate_diag(ma)) ? db_size : 0;

  const cs_lnum_t *row_index = cs_matrix_assembler_get_row_index(ma);
  const cs_lnum_t *col_ids = cs_matrix_assembler_get_col_ids(ma);

  PetscInt *_diag_sizes, *_offdiag_sizes;
  BFT_MALLOC(_diag_sizes, n_rows*db_size, PetscInt);
  BFT_MALLOC(_offdiag_sizes, n_rows*db_size, PetscInt);

  /* Separate local and distant loops for better first touch logic */

# pragma omp parallel for  if(n_rows > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_rows; i++) {
    cs_lnum_t s_id = row_index[i];
    cs_lnum_t e_id = row_index[i+1];

    PetscInt n_r_diag = 0;
    PetscInt n_cols = e_id - s_id;

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      if (col_ids[j] < n_rows) {
        n_r_diag++;
        if (col_ids[j] == i)
          n_r_diag += db_size -1;
      }
      else
        break;
    }

    for (cs_lnum_t j = 0; j < db_size; j++) {
      _diag_sizes[i*db_size + j] = n_r_diag + n_diag_add;
      _offdiag_sizes[i*db_size + j] = n_cols - n_r_diag;
    }
  }

  *diag_sizes = _diag_sizes;
  *offdiag_sizes = _offdiag_sizes;
}

/*----------------------------------------------------------------------------
 * Compute local and distant counts of matrix entries using an assembler
 * with full blocks fill.
 *
 * The caller is responsible for freeing the returned diag_sizes and
 * offdiag_sizes arrays.
 *
 * parameters:
 *   ma            <-- pointer to matrix assembler structure
 *   b_size        <-- diagonal block size
 *   diag_sizes    --> diagonal values
 *   offdiag_sizes --> off-diagonal values
 *----------------------------------------------------------------------------*/

static void
_compute_diag_sizes_assembler_b(const cs_matrix_assembler_t   *ma,
                                cs_lnum_t                      b_size,
                                PetscInt                    **diag_sizes,
                                PetscInt                    **offdiag_sizes)
{
  const cs_lnum_t n_rows = cs_matrix_assembler_get_n_rows(ma);

  cs_lnum_t n_diag_add
    = (cs_matrix_assembler_get_separate_diag(ma)) ? 1 : 0;

  const cs_lnum_t *row_index = cs_matrix_assembler_get_row_index(ma);
  const cs_lnum_t *col_ids = cs_matrix_assembler_get_col_ids(ma);

  PetscInt *_diag_sizes, *_offdiag_sizes;
  BFT_MALLOC(_diag_sizes, n_rows*b_size, PetscInt);
  BFT_MALLOC(_offdiag_sizes, n_rows*b_size, PetscInt);

  /* Separate local and distant loops for better first touch logic */

# pragma omp parallel for  if(n_rows > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_rows; i++) {
    cs_lnum_t s_id = row_index[i];
    cs_lnum_t e_id = row_index[i+1];

    PetscInt n_r_diag = 0;
    PetscInt n_cols = e_id - s_id;

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      if (col_ids[j] < n_rows) {
        n_r_diag++;
      }
      else
        break;
    }

    for (cs_lnum_t j = 0; j < b_size; j++) {
      _diag_sizes[i*b_size + j] = (n_r_diag + n_diag_add)*b_size;
      _offdiag_sizes[i*b_size + j] = (n_cols - n_r_diag)*b_size;
    }
  }

  *diag_sizes = _diag_sizes;
  *offdiag_sizes = _offdiag_sizes;
}

/*----------------------------------------------------------------------------
 * Compute local and distant counts of a native matrix's entries.
 *
 * The caller is responsible for freeing the returned diag_sizes and
 * offdiag_sizes arrays.
 *
 * parameters:
 *   matrix        <-- pointer to matrix structure
 *   have_diag     <-- does the matrix include a diagonal ?
 *   n_edges       <-- local number of graph edges
 *   edges         <-- edges (symmetric row <-> column) connectivity
 *   g_e_id        <-- global element ids
 *   da            <-- diagonal values
 *   diag_sizes    --> diagonal sizes
 *   offdiag_sizes --> off-diagonal sizes
 *----------------------------------------------------------------------------*/

static void
_compute_diag_sizes_native(cs_matrix_t        *matrix,
                           bool                have_diag,
                           cs_lnum_t           n_edges,
                           const cs_lnum_t     edges[restrict][2],
                           const cs_gnum_t     g_e_id[],
                           const cs_real_t     da[restrict],
                           PetscInt         **diag_sizes,
                           PetscInt         **offdiag_sizes)
{
  cs_lnum_t  n_rows = matrix->n_rows;
  cs_lnum_t  b_size = matrix->db_size;
  cs_lnum_t  e_size = matrix->eb_size;
  cs_lnum_t  b_stride = b_size * b_size;

  cs_lnum_t _n_rows = n_rows*b_size;

  PetscInt *_diag_sizes, *_offdiag_sizes;
  BFT_MALLOC(_diag_sizes, _n_rows, PetscInt);
  BFT_MALLOC(_offdiag_sizes, _n_rows, PetscInt);

  /* Case with b_size > e_size handled later */
  int n_diag = (have_diag && b_size == e_size) ? e_size : 0;

  cs_gnum_t g_id_lb = 0;
  cs_gnum_t g_id_ub = 0;

  if (_n_rows > 0) {
    g_id_lb = g_e_id[0];
    g_id_ub = g_e_id[n_rows-1] + 1;
  }

  /* Separate local and distant loops for better first touch logic */

# pragma omp parallel for  if(_n_rows > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < _n_rows; i++) {
    _diag_sizes[i] = n_diag;
    _offdiag_sizes[i] = 0;
  }

  const cs_lnum_t *group_index = NULL;
  if (matrix->numbering != NULL) {
    if (matrix->numbering->type == CS_NUMBERING_THREADS) {
      group_index = matrix->numbering->group_index;
    }
  }

  if (group_index != NULL) {

    const int n_threads = matrix->numbering->n_threads;
    const int n_groups = matrix->numbering->n_groups;

    for (int g_id = 0; g_id < n_groups; g_id++) {

#     pragma omp parallel for
      for (int t_id = 0; t_id < n_threads; t_id++) {

        for (cs_lnum_t edge_id = group_index[(t_id*n_groups + g_id)*2];
             edge_id < group_index[(t_id*n_groups + g_id)*2 + 1];
             edge_id++) {
          cs_lnum_t ii = edges[edge_id][0];
          cs_lnum_t jj = edges[edge_id][1];
          cs_gnum_t g_ii = g_e_id[ii];
          cs_gnum_t g_jj = g_e_id[jj];
          if (ii < n_rows) {
            if (g_jj >= g_id_lb && g_jj < g_id_ub)
              _diag_sizes[ii*b_size] += e_size;
            else
              _offdiag_sizes[ii*b_size] += e_size;
          }
          if (jj < n_rows) {
            if (g_ii >= g_id_lb && g_ii < g_id_ub)
              _diag_sizes[jj*b_size] += e_size;
            else
              _offdiag_sizes[jj*b_size] += e_size;
          }
        }
      }

    }

  }
  else {

    for (cs_lnum_t edge_id = 0; edge_id < n_edges; edge_id++) {
      cs_lnum_t ii = edges[edge_id][0];
      cs_lnum_t jj = edges[edge_id][1];
      cs_gnum_t g_ii = g_e_id[ii];
      cs_gnum_t g_jj = g_e_id[jj];
      if (ii < n_rows) {
        if (g_jj >= g_id_lb && g_jj < g_id_ub)
          _diag_sizes[ii*b_size] += e_size;
        else
          _offdiag_sizes[ii*b_size] += e_size;
      }
      if (jj < n_rows) {
        if (g_ii >= g_id_lb && g_ii < g_id_ub)
          _diag_sizes[jj*b_size] += e_size;
        else
          _offdiag_sizes[jj*b_size] += e_size;
      }
    }

  }

  /* Adjustment for "block" cases */

  if (b_size > 1) {

#   pragma omp parallel for  if(_n_rows > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_rows; i++) {
      for (cs_lnum_t j = 1; j < b_size; j++) {
        _diag_sizes[i*b_size + j] = _diag_sizes[i*b_size];
        _offdiag_sizes[i*b_size + j] = _offdiag_sizes[i*b_size];
      }
    }

    /* Delayed diagonal terms for block diagonal */
    if (have_diag && e_size == 1) {
#     pragma omp parallel for  if(_n_rows > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_rows; i++) {
        for (cs_lnum_t j = 0; j < b_size; j++) {
          for (cs_lnum_t k = 0; k < b_size; k++) {
            cs_real_t a = da[i*b_stride + j*b_size + k];
            if (a < -0.0 || a > 0.0)
              _diag_sizes[i*b_size + j] += 1;
          }
        }
      }
    }

  }

  *diag_sizes = _diag_sizes;
  *offdiag_sizes = _offdiag_sizes;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize coefficients context structure.
 *
 * \param[in, out]  matrix      pointer to matrix structure
 * \param[in]       type_name   string matching PETSc matrix type name,
 *                              defaults to "MATAIJ" if NULL
 */
/*----------------------------------------------------------------------------*/

static void
_setup_coeffs(cs_matrix_t  *matrix,
              const char   *type_name)
{
  if (matrix->coeffs == NULL) {
    cs_matrix_petsc_ensure_init();

    cs_matrix_coeffs_petsc_t  *coeffs;
    BFT_MALLOC(coeffs, 1, cs_matrix_coeffs_petsc_t);
    memset(coeffs, 0, sizeof(cs_matrix_coeffs_petsc_t));
    coeffs->matrix_state = 0;

    MatType matype = (type_name == NULL) ? MATAIJ : type_name;

    PetscStrallocpy(matype, &coeffs->matype_r);
    coeffs->matype = NULL;

    matrix->coeffs = coeffs;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function for initialization of PETSc matrix coefficients using
 *        local row ids and column indexes.
 *
 * \warning  The matrix pointer must point to valid data when the selection
 *           function is called, so the life cycle of the data pointed to
 *           should be at least as long as that of the assembler values
 *           structure.
 *
 * \param[in, out]  matrix_p  untyped pointer to matrix description structure
 * \param[in]       db_size   optional diagonal block sizes
 * \param[in]       eb_size   optional extra-diagonal block sizes
 */
/*----------------------------------------------------------------------------*/

static void
_assembler_values_init(void        *matrix_p,
                       cs_lnum_t    db_size,
                       cs_lnum_t    eb_size)
{
  cs_matrix_t  *matrix = (cs_matrix_t *)matrix_p;

  if (matrix->coeffs == NULL) {
    assert(0);  /* Prefer to setup coeff earlier, to pass type */
    _setup_coeffs(matrix, NULL);
  }

  /* Associated matrix assembler */

  const cs_matrix_assembler_t  *ma = matrix->assembler;

  cs_matrix_coeffs_petsc_t  *coeffs = matrix->coeffs;

  /* Create PETSc matrix */

  Mat hm = coeffs->hm;

  if (coeffs->matrix_state == 0) {

    const cs_gnum_t *l_range = cs_matrix_assembler_get_l_range(ma);
    const cs_gnum_t n_g_rows = cs_matrix_assembler_get_n_g_rows(ma);

    coeffs->l_range[0] = l_range[0];
    coeffs->l_range[1] = l_range[1];

    MPI_Comm comm = cs_glob_mpi_comm;
    if (comm == MPI_COMM_NULL)
      comm = MPI_COMM_SELF;

    MatCreate(comm, &hm);
    MatSetType(hm, coeffs->matype_r);

    PetscInt b_size = db_size;
    PetscInt n_rows = b_size * (l_range[1] - l_range[0]);
    PetscInt _n_g_rows = n_g_rows * b_size;

    MatSetSizes(hm,
                n_rows,      /* Number of local rows */
                n_rows,      /* Number of local columns */
                _n_g_rows,   /* Number of global rows */
                _n_g_rows);  /* Number of global columns */

    MatGetType(hm, &(coeffs->matype));

    coeffs->hm = hm;

    PetscInt *diag_sizes = NULL, *offdiag_sizes = NULL;

    if (db_size == 1)
      _compute_diag_sizes_assembler(ma,
                                    &diag_sizes,
                                    &offdiag_sizes);
    else if (eb_size == 1)
      _compute_diag_sizes_assembler_db(ma,
                                       db_size,
                                       &diag_sizes,
                                       &offdiag_sizes);
    else
      _compute_diag_sizes_assembler_b(ma,
                                      db_size,
                                      &diag_sizes,
                                      &offdiag_sizes);

    MatSeqAIJSetPreallocation(coeffs->hm, 0, diag_sizes);
    MatMPIAIJSetPreallocation(coeffs->hm, 0, diag_sizes, 0, offdiag_sizes);

    BFT_FREE(diag_sizes);
    BFT_FREE(offdiag_sizes);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add PETSc matrix coefficients using global row ids
 *        and column indexes, using intermediate copy for indexes and values.
 *
 * \param[in, out]  coeffs    PETSc Matrix coefficients handler
 * \param[in]       n         number of values to add
 * \param[in]       b_size    associated index block size
 * \param[in]       stride    associated data block size
 * \param[in]       row_g_id  associated global row ids
 * \param[in]       col_g_id  associated global column ids
 * \param[in]       vals      pointer to values (size: n*stride)
 */
/*----------------------------------------------------------------------------*/

static void
_assembler_values_add_block_cc(cs_matrix_coeffs_petsc_t  *coeffs,
                               PetscInt                   n,
                               PetscInt                   b_size,
                               PetscInt                   stride,
                               const cs_gnum_t            row_g_id[],
                               const cs_gnum_t            col_g_id[],
                               const cs_real_t            vals[])
{
  Mat hm = coeffs->hm;
  assert(hm != NULL);

  PetscInt h_b_size = b_size;
  PetscInt l_b = coeffs->l_range[0];
  PetscInt u_b = coeffs->l_range[1];

  for (PetscInt i = 0; i < n; i++) {
    PetscInt r_g_id = row_g_id[i];
    if (r_g_id >= l_b && r_g_id < u_b) {
      for (cs_lnum_t j = 0; j < b_size; j++) {
        for (cs_lnum_t k = 0; k < b_size; k++) {
          PetscInt idxm[] = {row_g_id[i]*h_b_size + (PetscInt)j};
          PetscInt idxn[] = {col_g_id[i]*h_b_size + (PetscInt)k};
          PetscScalar v[] = {vals[i*stride + j*b_size + k]};
          MatSetValues(hm, 1, idxm, 1, idxn, v, ADD_VALUES);
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add extradiagonal PETSc matrix coefficients using global row ids
 *        and column indexes, for fill type CS_MATRIX_BLOCK_D,
 *        CS_MATRIX_BLOCK_D_66, CS_MATRIX_BLOCK_D_SYM.
 *
 * \param[in, out]  coeffs    PETSc Matrix coefficients handler
 * \param[in]       n         number of values to add
 * \param[in]       b_size    associated data block size
 * \param[in]       stride    associated data block size
 * \param[in]       row_g_id  associated global row ids
 * \param[in]       col_g_id  associated global column ids
 * \param[in]       vals      pointer to values (size: n*stride)
 */
/*----------------------------------------------------------------------------*/

static void
_assembler_values_add_block_d_e(cs_matrix_coeffs_petsc_t  *coeffs,
                                PetscInt                  n,
                                PetscInt                  b_size,
                                const cs_gnum_t            row_g_id[],
                                const cs_gnum_t            col_g_id[],
                                const cs_real_t            vals[])
{
  Mat hm = coeffs->hm;
  assert(hm != NULL);

  PetscInt h_b_size = b_size;

  PetscInt l_b = coeffs->l_range[0];
  PetscInt u_b = coeffs->l_range[1];

  for (PetscInt i = 0; i < n; i++) {

    PetscInt r_g_id = row_g_id[ i];
    if (r_g_id >= l_b && r_g_id < u_b) {
      for (cs_lnum_t j = 0; j < b_size; j++) {
        PetscInt idxm[] = {row_g_id[i]*h_b_size + (PetscInt)j};
        PetscInt idxn[] = {col_g_id[i]*h_b_size + (PetscInt)j};
        PetscScalar v[] = {vals[i]};
        MatSetValues(hm, 1, idxm, 1, idxn, v, ADD_VALUES);
      }
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function for addition to PETSc matrix coefficients using
 *        local row ids and column indexes.
 *
 * This function can be used in all cases, including when
 *  sizeof(PetscInt) != sizeof(cs_gnum_t)
 *  sizeof(PetscReal) != sizeof(cs_real_t)
 *
 * Values whose associated row index is negative should be ignored;
 * Values whose column index is -1 are assumed to be assigned to a
 * separately stored diagonal. Other indexes should be valid.
 *
 * \warning  The matrix pointer must point to valid data when the selection
 *           function is called, so the life cycle of the data pointed to
 *           should be at least as long as that of the assembler values
 *           structure.
 *
 * \remark  Note that we pass column indexes (not ids) here; as the
 *          caller is already assumed to have identified the index
 *          matching a given column id.
 *
 * \param[in, out]  matrix_p  untyped pointer to matrix description structure
 * \param[in]       n         number of values to add
 * \param[in]       stride    associated data block size
 * \param[in]       row_g_id  associated global row ids
 * \param[in]       col_g_id  associated global column ids
 * \param[in]       vals      pointer to values (size: n*stride)
 */
/*----------------------------------------------------------------------------*/

static void
_assembler_values_add_g(void             *matrix_p,
                        cs_lnum_t         n,
                        cs_lnum_t         stride,
                        const cs_gnum_t   row_g_id[],
                        const cs_gnum_t   col_g_id[],
                        const cs_real_t   vals[])
{
  cs_matrix_t  *matrix = (cs_matrix_t *)matrix_p;
  cs_matrix_coeffs_petsc_t  *coeffs = matrix->coeffs;

  PetscInt nrows = n;
  PetscInt b_size = matrix->db_size;

  /* Scalar matrix
     ------------- */

  if (b_size == 1) {

    Mat hm = coeffs->hm;
    assert(hm != NULL);

    PetscInt l_b = coeffs->l_range[0];
    PetscInt u_b = coeffs->l_range[1];

    for (PetscInt i = 0; i < nrows; i++) {
      PetscInt r_g_id = row_g_id[i];
      if (r_g_id >= l_b && r_g_id < u_b) {
        PetscInt idxm[] = {r_g_id};
        PetscInt idxn[] = {col_g_id[i]};
        PetscScalar v[] = {vals[i]};
        MatSetValues(hm, 1, idxm, 1, idxn, v, ADD_VALUES);
      }
    }
  }

  /* Block matrix
     ------------ */

  else {

    /* Full blocks (including diagonal terms for diagonal fill) */

    if (   matrix->fill_type >= CS_MATRIX_BLOCK
        || row_g_id[0] == col_g_id[0])
      _assembler_values_add_block_cc(coeffs,
                                     nrows,
                                     b_size,
                                     stride,
                                     row_g_id,
                                     col_g_id,
                                     vals);

    /* Diagonal bloc extra-diagonal terms only */

    else if (matrix->fill_type >= CS_MATRIX_BLOCK_D)
      _assembler_values_add_block_d_e(coeffs,
                                      nrows,
                                      b_size,
                                      row_g_id,
                                      col_g_id,
                                      vals);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function to start the final assembly of matrix coefficients.
 *
 * \warning  The matrix pointer must point to valid data when the selection
 *           function is called, so the life cycle of the data pointed to
 *           should be at least as long as that of the assembler values
 *           structure.
 *
 * \param[in, out]  matrix_p  untyped pointer to matrix description structure
 */
/*----------------------------------------------------------------------------*/

static void
_assembler_values_begin(void  *matrix_p)
{
  cs_matrix_t  *matrix = (cs_matrix_t *)matrix_p;
  cs_matrix_coeffs_petsc_t  *coeffs = matrix->coeffs;

  MatAssemblyBegin(coeffs->hm, MAT_FINAL_ASSEMBLY);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function to complete the final assembly of matrix
 *        coefficients.
 *
 * \warning  The matrix pointer must point to valid data when the selection
 *           function is called, so the life cycle of the data pointed to
 *           should be at least as long as that of the assembler values
 *           structure.
 *
 * \param[in, out]  matrix_p  untyped pointer to matrix description structure
 */
/*----------------------------------------------------------------------------*/

static void
_assembler_values_end(void  *matrix_p)
{
  cs_matrix_t  *matrix = (cs_matrix_t *)matrix_p;
  cs_matrix_coeffs_petsc_t  *coeffs = matrix->coeffs;

  MatAssemblyEnd(coeffs->hm, MAT_FINAL_ASSEMBLY);

  if (coeffs->matrix_state == 0) {

    MPI_Comm comm = cs_glob_mpi_comm;
    if (comm == MPI_COMM_NULL) {
      comm = MPI_COMM_SELF;
    }

    /* Create associated vectors here also to avoid repeated creation
       (and possible overhead) where used */

    const PetscInt b_size = matrix->db_size;
    const PetscInt n_local = matrix->n_rows * b_size;

    Vec hx, hy;
    VecCreate(comm, &hx);
    VecSetSizes(hx, n_local, PETSC_DECIDE);
    VecSetBlockSize(hx, 1);

    if (   strcmp(coeffs->matype, MATAIJCUSPARSE) == 0
        || strcmp(coeffs->matype, MATSEQAIJCUSPARSE) == 0
        || strcmp(coeffs->matype, MATMPIAIJCUSPARSE) == 0)
      VecSetType(hx, VECCUDA);
    else
      VecSetType(hx, VECSTANDARD);

    /* VecSetFromOptions(Vec v); */

    VecDuplicate(hx, &hy);

    coeffs->hx = hx;
    coeffs->hy = hy;

  }

  /* Set stat flag */

  coeffs->matrix_state = 1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize a CSR matrix assembler values structure.
 *
 * The associated matrix's structure must have been created using
 * \ref cs_matrix_structure_create_from_assembler.
 *
 * \param[in, out]  matrix                 pointer to matrix structure
 * \param[in]       diag_block_size        block sizes for diagonal
 * \param[in]       extra_diag_block_size  block sizes for extra diagonal
 *
 * \return  pointer to initialized matrix assembler values structure;
 */
/*----------------------------------------------------------------------------*/

static cs_matrix_assembler_values_t *
_assembler_values_create_petsc(cs_matrix_t      *matrix,
                               const cs_lnum_t   diag_block_size,
                               const cs_lnum_t   extra_diag_block_size)
{
  cs_matrix_assembler_values_t *mav
    = cs_matrix_assembler_values_create(matrix->assembler,
                                        false,
                                        diag_block_size,
                                        extra_diag_block_size,
                                        (void *)matrix,
                                        _assembler_values_init,
                                        NULL,
                                        _assembler_values_add_g,
                                        _assembler_values_begin,
                                        _assembler_values_end);

  return mav;
}

/*----------------------------------------------------------------------------
 * Set PETSc matrix coefficients for block-diagonal cases.
 *
 * parameters:
 *   matrix    <-- pointer to matrix structure
 *   symmetric <-- indicates if extradiagonal values are symmetric
 *   n_edges   <-- local number of graph edges
 *   edges     <-- edges (symmetric row <-> column) connectivity
 *   da        <-- diagonal values
 *   xa        <-- extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_coeffs_ij_db(cs_matrix_t        *matrix,
                  bool                symmetric,
                  cs_lnum_t           n_edges,
                  const cs_lnum_t     edges[restrict][2],
                  const cs_real_t     da[restrict],
                  const cs_real_t     xa[restrict])
{
  const cs_lnum_t  n_rows = matrix->n_rows;

  const cs_gnum_t *g_id = cs_matrix_get_block_row_g_id(matrix);
  const bool have_diag = (xa != NULL) ? true : false;

  MPI_Comm comm = cs_glob_mpi_comm;
  if (comm == MPI_COMM_NULL)
    comm = MPI_COMM_SELF;

  assert(n_rows > 0);

  InsertMode insert_mode = ADD_VALUES; /* INSERT_VALUES possible
                                          with direct_assembly */

  cs_matrix_coeffs_petsc_t  *coeffs = matrix->coeffs;

  /* Sizes and buffers */

  PetscInt h_b_size = matrix->db_size;

  cs_lnum_t b_size = matrix->db_size;
  cs_lnum_t b_stride = b_size * b_size;

  /* Diagonal part
     ------------- */

  if (have_diag) {

    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      for (cs_lnum_t j = 0; j < b_size; j++) {
        for (cs_lnum_t k = 0; k < b_size; k++) {
          cs_real_t a = da[ii*b_stride + j*b_size + k];
          if (a < -0.0 || a > 0.0) {
            PetscInt idxm[] = {g_id[ii]*h_b_size + (PetscInt)j};
            PetscInt idxn[] = {g_id[ii]*h_b_size + (PetscInt)k};
            PetscScalar v[] = {a};
            MatSetValues(coeffs->hm, 1, idxm, 1, idxn, v, insert_mode);
          }
        }
      }

    }

  }

  /* Extradiagonal part */

  if (symmetric) {

    for (cs_lnum_t e_id = 0; e_id < n_edges; e_id++) {
      cs_lnum_t ii = edges[e_id][0];
      cs_lnum_t jj = edges[e_id][1];
      cs_gnum_t g_ii = g_id[ii];
      cs_gnum_t g_jj = g_id[jj];
      if (ii < n_rows) {
        for (cs_lnum_t k = 0; k < b_size; k++) {
          PetscInt idxm[] = {g_ii*h_b_size + (PetscInt)k};
          PetscInt idxn[] = {g_jj*h_b_size + (PetscInt)k};
          PetscScalar v[] = {xa[e_id]};
          MatSetValues(coeffs->hm, 1, idxm, 1, idxn, v, insert_mode);
        }
      }
      if (jj < n_rows) {
        for (cs_lnum_t k = 0; k < b_size; k++) {
          PetscInt idxm[] = {g_jj*h_b_size + (PetscInt)k};
          PetscInt idxn[] = {g_ii*h_b_size + (PetscInt)k};
          PetscScalar v[] = {xa[e_id]};
          MatSetValues(coeffs->hm, 1, idxm, 1, idxn, v, insert_mode);
        }
      }
    }

  }
  else { /* non-symmetric variant */

    for (cs_lnum_t e_id = 0; e_id < n_edges; e_id++) {
      cs_lnum_t ii = edges[e_id][0];
      cs_lnum_t jj = edges[e_id][1];
      cs_gnum_t g_ii = g_id[ii];
      cs_gnum_t g_jj = g_id[jj];
      if (ii < n_rows) {
        for (cs_lnum_t k = 0; k < b_size; k++) {
          PetscInt idxm[] = {g_ii*h_b_size + (PetscInt)k};
          PetscInt idxn[] = {g_jj*h_b_size + (PetscInt)k};
          PetscScalar v[] = {xa[e_id*2]};
          MatSetValues(coeffs->hm, 1, idxm, 1, idxn, v, insert_mode);
        }
      }
      if (jj < n_rows) {
        for (cs_lnum_t k = 0; k < b_size; k++) {
          PetscInt idxm[] = {g_jj*h_b_size + (PetscInt)k};
          PetscInt idxn[] = {g_ii*h_b_size + (PetscInt)k};
          PetscScalar v[] = {xa[e_id*2+1]};
          MatSetValues(coeffs->hm, 1, idxm, 1, idxn, v, insert_mode);
        }
      }
    }

  }
}

/*----------------------------------------------------------------------------
 * Set PETSc matrix coefficients for full block cases.
 *
 * parameters:
 *   matrix    <-- pointer to matrix structure
 *   symmetric <-- indicates if extradiagonal values are symmetric
 *   n_edges   <-- local number of graph edges
 *   edges     <-- edges (symmetric row <-> column) connectivity
 *   da        <-- diagonal values
 *   xa        <-- extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_coeffs_ij_b(cs_matrix_t        *matrix,
                 bool                symmetric,
                 cs_lnum_t           n_edges,
                 const cs_lnum_t     edges[restrict][2],
                 const cs_real_t     da[restrict],
                 const cs_real_t     xa[restrict])
{
  const cs_lnum_t  n_rows = matrix->n_rows;

  const cs_gnum_t *g_id = cs_matrix_get_block_row_g_id(matrix);
  const bool have_diag = (xa != NULL) ? true : false;

  MPI_Comm comm = cs_glob_mpi_comm;
  if (comm == MPI_COMM_NULL)
    comm = MPI_COMM_SELF;

  assert(n_rows > 0);

  InsertMode insert_mode = ADD_VALUES; /* INSERT_VALUES possible
                                          with direct_assembly */

  cs_matrix_coeffs_petsc_t  *coeffs = matrix->coeffs;

  /* Sizes and buffers */

  PetscInt h_b_size = matrix->db_size;

  cs_lnum_t b_size = matrix->db_size;
  cs_lnum_t b_stride = b_size * b_size;

  PetscInt idxm[9];
  PetscInt idxn[9];
  PetscScalar v[9*9];

  assert(b_size <= 9);

  /* Diagonal part
     ------------- */

  if (have_diag) {

    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
      for (cs_lnum_t j = 0; j < b_size; j++) {
        idxm[j] = g_id[ii]*h_b_size + (PetscInt)j;
        idxn[j] = g_id[ii]*h_b_size + (PetscInt)j;
        for (cs_lnum_t k = 0; k < b_size; k++) {
          v[j*b_size + k] = da[ii*b_stride + j*b_size + k];
        }
      }
      MatSetValues(coeffs->hm, h_b_size, idxm, h_b_size, idxn, v, insert_mode);
    }

  }  /* End of diagonal block addition */

  /* Extradiagonal part */

  if (symmetric) {

    for (cs_lnum_t e_id = 0; e_id < n_edges; e_id++) {

      cs_lnum_t ii = edges[e_id][0];
      cs_lnum_t jj = edges[e_id][1];
      cs_gnum_t g_ii = g_id[ii];
      cs_gnum_t g_jj = g_id[jj];
      if (ii < n_rows) {
        for (cs_lnum_t k = 0; k < b_size; k++) {
          idxm[k] = g_ii*h_b_size + (PetscInt)k;
          idxn[k] = g_jj*h_b_size + (PetscInt)k;
          for (cs_lnum_t l = 0; l < b_size; l++) {
            v[k*b_size + l] = xa[e_id*b_stride + k*b_size + l];
          }
        }
        MatSetValues(coeffs->hm, h_b_size, idxm, h_b_size, idxn, v, insert_mode);
      }
      if (jj < n_rows) {
        for (cs_lnum_t k = 0; k < b_size; k++) {
          idxm[k] = g_jj*h_b_size + (PetscInt)k;
          idxn[k] = g_ii*h_b_size + (PetscInt)k;
          for (cs_lnum_t l = 0; l < b_size; l++) {
            v[k*b_size + l] = xa[e_id*b_stride + k*b_size + l];
          }
        }
        MatSetValues(coeffs->hm, h_b_size, idxm, h_b_size, idxn, v, insert_mode);
      }

    }

  }
  else { /* non-symmetric variant */

    for (cs_lnum_t e_id = 0; e_id < n_edges; e_id++) {

      cs_lnum_t ii = edges[e_id][0];
      cs_lnum_t jj = edges[e_id][1];
      cs_gnum_t g_ii = g_id[ii];
      cs_gnum_t g_jj = g_id[jj];
      if (ii < n_rows) {
        for (cs_lnum_t k = 0; k < b_size; k++) {
          idxm[k] = g_ii*h_b_size + (PetscInt)k;
          idxn[k] = g_jj*h_b_size + (PetscInt)k;
          for (cs_lnum_t l = 0; l < b_size; l++) {
            v[k*b_size + l] = xa[e_id*2*b_stride + k*b_size + l];
          }
        }
        MatSetValues(coeffs->hm, h_b_size, idxm, h_b_size, idxn, v, insert_mode);
      }
      if (jj < n_rows) {
        for (cs_lnum_t k = 0; k < b_size; k++) {
          idxm[k] = g_jj*h_b_size + (PetscInt)k;
          idxn[k] = g_ii*h_b_size + (PetscInt)k;
          for (cs_lnum_t l = 0; l < b_size; l++) {
            v[k*b_size + l] = xa[(e_id*2+1)*b_stride + k*b_size + l];
          }
        }
        MatSetValues(coeffs->hm, h_b_size, idxm, h_b_size, idxn, v, insert_mode);
      }

    }

  }
}

/*----------------------------------------------------------------------------
 * Set PETSc matrix coefficients.
 *
 * parameters:
 *   matrix    <-- pointer to matrix structure
 *   symmetric <-- indicates if extradiagonal values are symmetric
 *   copy      <-- indicates if coefficients should be copied
 *   n_edges   <-- local number of graph edges
 *   edges     <-- edges (symmetric row <-> column) connectivity
 *   da        <-- diagonal values
 *   xa        <-- extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_coeffs(cs_matrix_t        *matrix,
            bool                symmetric,
            bool                copy,
            cs_lnum_t           n_edges,
            const cs_lnum_t     edges[restrict][2],
            const cs_real_t     da[restrict],
            const cs_real_t     xa[restrict])
{
  CS_UNUSED(copy);

  cs_lnum_t  n_rows = matrix->n_rows;
  cs_lnum_t  n_cols_ext = matrix->n_cols_ext;

  const cs_gnum_t *g_id = cs_matrix_get_block_row_g_id(matrix);
  const bool have_diag = (xa != NULL) ? true : false;

  MPI_Comm comm = cs_glob_mpi_comm;
  if (comm == MPI_COMM_NULL)
    comm = MPI_COMM_SELF;

  assert(n_rows > 0);

  InsertMode insert_mode = ADD_VALUES; /* INSERT_VALUES possible
                                          with direct_assembly */

  cs_matrix_coeffs_petsc_t  *coeffs = matrix->coeffs;

  /* Create PETSc matrix */

  Mat hm = coeffs->hm;

  PetscInt b_size = matrix->db_size;
  PetscInt e_size = matrix->eb_size;

  if (coeffs->matrix_state == 0) {

    cs_gnum_t  l_range[2] = {0, 0};
    if (n_rows > 0) {
      l_range[0] = g_id[0];
      l_range[1] = g_id[n_rows-1] + 1;
    }

    coeffs->l_range[0] = l_range[0];
    coeffs->l_range[1] = l_range[1];

    n_rows *= b_size;
    n_cols_ext *= b_size;

    MatCreate(comm, &hm);
    MatSetType(hm, coeffs->matype_r);

    MatSetSizes(hm,
                (PetscInt)n_rows,    /* Number of local rows */
                (PetscInt)n_rows,    /* Number of local columns */
                PETSC_DETERMINE,     /* Number of global rows */
                PETSC_DETERMINE);    /* Number of global columns */

    coeffs->hm = hm;

    MatGetType(hm, &(coeffs->matype));

    PetscInt *diag_sizes = NULL, *offdiag_sizes = NULL;

    _compute_diag_sizes_native(matrix,
                               have_diag,
                               n_edges,
                               edges,
                               g_id,
                               da,
                               &diag_sizes,
                               &offdiag_sizes);

    MatSeqAIJSetPreallocation(coeffs->hm, 0, diag_sizes);
    MatMPIAIJSetPreallocation(coeffs->hm, 0, diag_sizes, 0, offdiag_sizes);

    BFT_FREE(diag_sizes);
    BFT_FREE(offdiag_sizes);
  }

  /* Scalar case */

  if (b_size == 1) {

    if (have_diag) {

      for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
        PetscInt idxm[] = {g_id[ii]};
        PetscInt idxn[] = {g_id[ii]};
        PetscScalar v[] = {da[ii]};
        MatSetValues(coeffs->hm, 1, idxm, 1, idxn, v, insert_mode);
      }

    }

    if (symmetric) {

      for (cs_lnum_t e_id = 0; e_id < n_edges; e_id++) {
        cs_lnum_t ii = edges[e_id][0];
        cs_lnum_t jj = edges[e_id][1];
        cs_gnum_t g_ii = g_id[ii];
        cs_gnum_t g_jj = g_id[jj];
        if (ii < n_rows) {
          PetscInt idxm[] = {g_ii};
          PetscInt idxn[] = {g_jj};
          PetscScalar v[] = {xa[e_id]};
          MatSetValues(coeffs->hm, 1, idxm, 1, idxn, v, insert_mode);
        }
        if (jj < n_rows) {
          PetscInt idxm[] = {g_jj};
          PetscInt idxn[] = {g_ii};
          PetscScalar v[] = {xa[e_id]};
          MatSetValues(coeffs->hm, 1, idxm, 1, idxn, v, insert_mode);
        }
      }

    }
    else { /* non-symmetric variant */

      for (cs_lnum_t e_id = 0; e_id < n_edges; e_id++) {
        cs_lnum_t ii = edges[e_id][0];
        cs_lnum_t jj = edges[e_id][1];
        cs_gnum_t g_ii = g_id[ii];
        cs_gnum_t g_jj = g_id[jj];
        if (ii < n_rows) {
          PetscInt idxm[] = {g_ii};
          PetscInt idxn[] = {g_jj};
          PetscScalar v[] = {xa[e_id*2]};
          MatSetValues(coeffs->hm, 1, idxm, 1, idxn, v, insert_mode);
        }
        if (jj < n_rows) {
          PetscInt idxm[] = {g_jj};
          PetscInt idxn[] = {g_ii};
          PetscScalar v[] = {xa[e_id*2+1]};
          MatSetValues(coeffs->hm, 1, idxm, 1, idxn, v, insert_mode);
        }
      }

    }

  }

  /* Block diagonal only */

  else if (b_size > 1 && e_size == 1)
    _set_coeffs_ij_db(matrix,
                      symmetric,
                      n_edges,
                      edges,
                      da,
                      xa);

  /* Full block */

  else /* if (b_size > 1 && e_size > 1) */
    _set_coeffs_ij_b(matrix,
                     symmetric,
                     n_edges,
                     edges,
                     da,
                     xa);

  /* Finalize asembly and update state */

  _assembler_values_end(matrix);
}

/*----------------------------------------------------------------------------
 * Release PETSc matrix coefficients.
 *
 * Depending on current options and initialization, values will be copied
 * or simply mapped.
 *
 * parameters:
 *   matrix    <-- pointer to matrix structure
 *   symmetric <-- indicates if extradiagonal values are symmetric
 *   copy      <-- indicates if coefficients should be copied
 *   n_edges   <-- local number of graph edges
 *   edges     <-- edges (symmetric row <-> column) connectivity
 *   da        <-- diagonal values
 *   xa        <-- extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_release_coeffs(cs_matrix_t  *matrix)
{
  cs_matrix_coeffs_petsc_t  *coeffs = matrix->coeffs;

  if (matrix->coeffs != NULL) {

    if (coeffs->matrix_state > 0) {
      MatDestroy(&(coeffs->hm));

      VecDestroy(&(coeffs->hx));
      VecDestroy(&(coeffs->hy));

      coeffs->matrix_state = 0;
    }

    if (coeffs->matype_r != NULL)
      PetscFree(coeffs->matype_r);
  }
}

/*----------------------------------------------------------------------------
 * Release PETSc matrix coefficients.
 *
 * Depending on current options and initialization, values will be copied
 * or simply mapped.
 *
 * parameters:
 *   matrix    <-- pointer to matrix structure
 *   symmetric <-- indicates if extradiagonal values are symmetric
 *   copy      <-- indicates if coefficients should be copied
 *   n_edges   <-- local number of graph edges
 *   edges     <-- edges (symmetric row <-> column) connectivity
 *   da        <-- diagonal values
 *   xa        <-- extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_destroy_coeffs(cs_matrix_t  *matrix)
{
  if (matrix->coeffs != NULL) {
    _release_coeffs(matrix);
    BFT_FREE(matrix->coeffs);
  }
}

/*----------------------------------------------------------------------------
 * Copy diagonal of PETSc matrix.
 *
 * parameters:
 *   matrix <-- pointer to matrix structure
 *   da     --> diagonal (pre-allocated, size: n_rows)
 *----------------------------------------------------------------------------*/

static void
_copy_diagonal(const cs_matrix_t  *matrix,
               cs_real_t          *restrict da)
{
  const cs_matrix_coeffs_petsc_t  *coeffs = matrix->coeffs;

  const PetscInt b_size = matrix->db_size;

  const PetscInt n_rows = matrix->n_rows * b_size;

  Vec hd;
  VecDuplicate(coeffs->hx, &hd);

  MatGetDiagonal(coeffs->hm, hd);

  const PetscReal *_d = NULL;

  VecGetArrayRead(hd, &_d);

  for (PetscInt ii = 0; ii < n_rows; ii++) {
    da[ii] = _d[ii];
  }

  VecRestoreArrayRead(hd, &_d);
  VecDestroy(&hd);
}

/*============================================================================
 * Semi-private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize PETSc if needed
 *----------------------------------------------------------------------------*/

void
cs_matrix_petsc_ensure_init(void)
{
  if (_init_status > 0)
    return;

  /* Initialization must be called before setting options;
     it does not need to be called before calling
     cs_sles_petsc_define(), as this is handled automatically. */

  PetscBool is_initialized;
  PetscInitialized(&is_initialized);

  if (is_initialized == PETSC_FALSE) {
#if defined(HAVE_MPI)
    if (cs_glob_mpi_comm != MPI_COMM_NULL)
      PETSC_COMM_WORLD = cs_glob_mpi_comm;
    else {
      PETSC_COMM_WORLD = MPI_COMM_SELF;
      int flag = 0;
      MPI_Initialized(&flag);
      if (!flag)
        MPI_Init(NULL, NULL);
    }
#endif
    PetscInitializeNoArguments();
  }

  PetscPushErrorHandler(PetscAbortErrorHandler, NULL);
  _init_status = 1;
}

/*----------------------------------------------------------------------------
 * Finalize PETSc
 *----------------------------------------------------------------------------*/

void
cs_matrix_petsc_finalize(void)
{
  if (_init_status == 1) {
    PetscFinalize();
    _init_status = 2;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief return coefficients structure associated with PETSc matrix.
 *
 * \param[in]  matrix  pointer to matrix structure
 *
 * \return  pointer to matrix coefficients handler structure for PETSc matrix.
 */
/*----------------------------------------------------------------------------*/

cs_matrix_coeffs_petsc_t *
cs_matrix_petsc_get_coeffs(const cs_matrix_t  *matrix)
{
  return (cs_matrix_coeffs_petsc_t *)(matrix->coeffs);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Switch matrix type to PETSc.
 *
 * This releases previous coefficients if present, so should be called
 * just after matrix creation, before assigning coefficients.
 *
 * \param[in, out]  matrix      pointer to matrix structure
 * \param[in]       type_name   string matching PETSc matrix type name,
 *                              defaults to "MATAIJ" if NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_set_type_petsc(cs_matrix_t  *matrix,
                         const char   *type_name)
{
  matrix->type = CS_MATRIX_N_BUILTIN_TYPES;

  matrix->type_name = _petsc_ij_type_name;
  matrix->type_fname = _petsc_ij_type_fullname;

  /* Release previous coefficients if present */

  if (matrix->coeffs != NULL)
    matrix->destroy_coefficients(matrix);

  _setup_coeffs(matrix, type_name);

  cs_matrix_coeffs_petsc_t  *coeffs = matrix->coeffs;

  /* Set function pointers here */

  matrix->set_coefficients = _set_coeffs;
  matrix->release_coefficients = _release_coeffs;
  matrix->destroy_coefficients = _destroy_coeffs;
  matrix->assembler_values_create = _assembler_values_create_petsc;

  matrix->get_diagonal = NULL;

  /* Remark: block values are transformed into scalar values, so SpMv products
     should be possible, (and the function pointers updated). PETSc also has
     support for block matrixes, so we can exploit that in the future
     (we will need to make sure the correct type is used in this case). */

  for (int i = 0; i < CS_MATRIX_N_FILL_TYPES; i++) {
    matrix->vector_multiply[i][0] = _mat_vec_p_aij;
    matrix->vector_multiply[i][1] = NULL;
  }

  matrix->copy_diagonal = _copy_diagonal;

  /* Force MPI initialization if not already done.
   * The main code_saturne communicator is not modifed, as
   * this is purely for external libraries use. */

#if defined(HAVE_MPI)

  if (cs_glob_mpi_comm == MPI_COMM_NULL) {
    int mpi_flag;
    MPI_Initialized(&mpi_flag);
    if (mpi_flag == 0)
      MPI_Init(NULL, NULL);
  }

#endif

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
