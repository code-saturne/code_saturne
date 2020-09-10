/*============================================================================
 * Sparse Matrix Representation and Operations using HYPRE library.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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
 * HYPRE headers
 *----------------------------------------------------------------------------*/

#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE_utilities.h"

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
#include "cs_matrix_default.h"
#include "cs_matrix_hypre.h"
#include "cs_matrix_priv.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*! \file cs_matrix_hypre.c
 *
 * \brief Sparse Matrix Representation and Operations using HYPRE.
 */
/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* Note that most types are declared in cs_matrix_priv.h.
   only those only handled here are declared here. */

/* Adapter coefficients stucture for HYPRE */

typedef struct _cs_matrix_coeffs_hypre_t {

  HYPRE_IJMatrix hm;                       /* HYPRE matrix */
  HYPRE_IJVector hx;                       /* x (input) vector */
  HYPRE_IJVector hy;                       /* y (output) vector */

  int  matrix_state;                       /* Matrix state:
                                              0: not created
                                              1: created and assembled */

} cs_matrix_coeffs_hypre_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

static const char _hypre_ij_type_name[] = "HYPRE_PARCSR";
static const char _hypre_ij_type_fullname[] = "HYPRE IJ (HYPRE_ParCSR)";

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute local and distant counts of matrix entries.
 *
 * The caller is responsible for freeing the returned diag_sizes and
 * offdiag_sizes arrays.
 *
 * parameters:
 *   matrix        <-- pointer to matrix structure
 *   have_diag     <-- does the matrix include a diagonal ?
 *   n_edges       <-- local number of graph edges
 *   edges         <-- edges (symmetric row <-> column) connectivity
 *   diag_sizes    --> diagonal values
 *   offdiag_sizes --> off-diagonal values
 *----------------------------------------------------------------------------*/

static void
_compute_diag_sizes_native(cs_matrix_t        *matrix,
                           bool                have_diag,
                           cs_lnum_t           n_edges,
                           const cs_lnum_t     edges[restrict][2],
                           HYPRE_Int         **diag_sizes,
                           HYPRE_Int         **offdiag_sizes)
{
  cs_lnum_t  n_rows = matrix->n_rows;

  HYPRE_Int *_diag_sizes, *_offdiag_sizes;
  BFT_MALLOC(_diag_sizes, n_rows, HYPRE_Int);
  BFT_MALLOC(_offdiag_sizes, n_rows, HYPRE_Int);

  int n_diag = (have_diag) ? 1 : 0;

  /* Separate local and distant loops for better first touch logic */

# pragma omp parallel for  if(n_rows > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_rows; i++) {
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
          if (ii < n_rows) {
            if (jj < n_rows) {
              _diag_sizes[ii] += 1;
              _diag_sizes[jj] += 1;
            }
            else {
              _offdiag_sizes[ii] += 1;
            }
          }
          else {
            _offdiag_sizes[jj] += 1;
          }
        }
      }

    }

  }
  else {

    for (cs_lnum_t edge_id = 0; edge_id < n_edges; edge_id++) {
      cs_lnum_t ii = edges[edge_id][0];
      cs_lnum_t jj = edges[edge_id][1];
      if (ii < n_rows) {
        if (jj < n_rows) {
          _diag_sizes[ii] += 1;
          _diag_sizes[jj] += 1;
        }
        else {
          _offdiag_sizes[ii] += 1;
        }
      }
      else {
        _offdiag_sizes[jj] += 1;
      }
    }

  }

  *diag_sizes = _diag_sizes;
  *offdiag_sizes = _offdiag_sizes;
}

/*----------------------------------------------------------------------------
 * Set HYPRE ParCSR matrix coefficients.
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
_set_coeffs_ij(cs_matrix_t        *matrix,
               bool                symmetric,
               bool                copy,
               cs_lnum_t           n_edges,
               const cs_lnum_t     edges[restrict][2],
               const cs_real_t     da[restrict],
               const cs_real_t     xa[restrict])
{
  CS_UNUSED(copy);

  bool direct_assembly = false;

  const cs_lnum_t  n_rows = matrix->n_rows;

  const cs_gnum_t *g_id = cs_matrix_get_block_row_g_id(n_rows, matrix->halo);
  const bool have_diag = (xa != NULL) ? true : false;

  MPI_Comm comm = cs_glob_mpi_comm;
  if (comm == MPI_COMM_NULL)
    comm = MPI_COMM_WORLD;

  assert(n_rows > 0);

  cs_matrix_coeffs_hypre_t  *coeffs = matrix->coeffs;

  /* Create HYPRE matrix */

  HYPRE_IJMatrix hm = coeffs->hm;

  if (coeffs->matrix_state == 0) {

    HYPRE_BigInt ilower = g_id[0];
    HYPRE_BigInt iupper = g_id[n_rows-1];

    HYPRE_IJMatrixCreate(comm,
                         ilower,
                         iupper,
                         ilower,
                         iupper,
                         &hm);

    coeffs->hm = hm;

    HYPRE_IJMatrixSetObjectType(hm, HYPRE_PARCSR);

    HYPRE_Int *diag_sizes = NULL, *offdiag_sizes = NULL;

    _compute_diag_sizes_native(matrix,
                               have_diag,
                               n_edges,
                               edges,
                               &diag_sizes,
                               &offdiag_sizes);

    HYPRE_IJMatrixSetDiagOffdSizes(hm, diag_sizes, offdiag_sizes);
    HYPRE_IJMatrixSetMaxOffProcElmts(hm, 0);

    BFT_FREE(diag_sizes);
    BFT_FREE(offdiag_sizes);

    HYPRE_IJMatrixSetOMPFlag(hm, 0);

    HYPRE_MemoryLocation  memory_location = HYPRE_MEMORY_HOST;

    HYPRE_IJMatrixInitialize_v2(hm, memory_location);
  }

  HYPRE_Int max_chunk_size = 32768 - 1;

  HYPRE_BigInt *rows, *cols;
  HYPRE_Real *aij;
  BFT_MALLOC(rows, max_chunk_size+1, HYPRE_BigInt);
  BFT_MALLOC(cols, max_chunk_size+1, HYPRE_BigInt);
  BFT_MALLOC(aij, max_chunk_size+1, HYPRE_Real);

  if (have_diag) {
    cs_lnum_t s_id = 0;
    while (s_id < n_rows) {

      HYPRE_Int ic = 0;

      for (cs_lnum_t ii = s_id;
           ii < n_rows && ic < max_chunk_size;
           ii++) {
        rows[ic] = g_id[ii];
        cols[ic] = g_id[ii];
        aij[ic] = da[ii];
        ic++;
      }

      s_id = ic;

      if (direct_assembly)
        HYPRE_IJMatrixSetValues(hm, ic, NULL, rows, cols, aij);
      else
        HYPRE_IJMatrixAddToValues(hm, ic, NULL, rows, cols, aij);
    }
  }

  if (symmetric) {

    cs_lnum_t s_e_id = 0;

    while (s_e_id < n_edges) {

      HYPRE_Int ic = 0;

      for (cs_lnum_t e_id = s_e_id;
           e_id < n_edges && ic < max_chunk_size;
           e_id++) {
        cs_lnum_t ii = edges[e_id][0];
        cs_lnum_t jj = edges[e_id][1];

        if (ii < n_rows) {
          rows[ic] = g_id[ii];
          cols[ic] = g_id[jj];
          aij[ic] = xa[e_id];
          ic++;
        }
        if (jj < n_rows) {
          rows[ic] = g_id[jj];
          cols[ic] = g_id[ii];
          aij[ic] = xa[e_id];
          ic++;
        }
      }

      s_e_id = ic;

      if (direct_assembly)
        HYPRE_IJMatrixSetValues(hm, ic, NULL, rows, cols, aij);
      else
        HYPRE_IJMatrixAddToValues(hm, ic, NULL, rows, cols, aij);
    }

  }
  else { /* non-symmetric variant */

    cs_lnum_t s_e_id = 0;

    while (s_e_id < n_edges) {

      HYPRE_Int ic = 0;

      for (cs_lnum_t e_id = s_e_id;
           e_id < n_edges && ic < max_chunk_size;
           e_id++) {
        cs_lnum_t ii = edges[e_id][0];
        cs_lnum_t jj = edges[e_id][1];

        if (ii < n_rows) {
          rows[ic] = g_id[ii];
          cols[ic] = g_id[jj];
          aij[ic] = xa[e_id*2];
          ic++;
        }
        if (jj < n_rows) {
          rows[ic] = g_id[jj];
          cols[ic] = g_id[ii];
          aij[ic] = xa[e_id*2+1];
          ic++;
        }
      }

      s_e_id = ic;

      if (direct_assembly)
        HYPRE_IJMatrixSetValues(hm, ic, NULL, rows, cols, aij);
      else
        HYPRE_IJMatrixAddToValues(hm, ic, NULL, rows, cols, aij);

    }

  }

  BFT_FREE(rows);
  BFT_FREE(cols);
  BFT_FREE(aij);

  HYPRE_IJMatrixAssemble(hm);

  /* Create associated vectors here also to avoid repeated creation
     ² (and possible overhead) where used */

  if (coeffs->matrix_state == 0) {

    const HYPRE_Int  n_off_proc = matrix->n_cols_ext - matrix->n_rows;

    HYPRE_BigInt ilower, iupper, jlower, jupper;
    HYPRE_IJMatrixGetLocalRange(coeffs->hm,
                                &ilower,
                                &iupper,
                                &jlower,
                                &jupper);

    HYPRE_IJVectorCreate(comm, ilower, iupper - 1, &(coeffs->hx));
    HYPRE_IJVectorSetObjectType(coeffs->hx, HYPRE_PARCSR);
    HYPRE_IJVectorSetMaxOffProcElmts(coeffs->hx, n_off_proc);

    HYPRE_IJVectorCreate(comm, ilower, iupper - 1, &(coeffs->hy));
    HYPRE_IJVectorSetObjectType(coeffs->hy, HYPRE_PARCSR);
    HYPRE_IJVectorSetMaxOffProcElmts(coeffs->hy, n_off_proc);

    HYPRE_IJVectorInitialize(coeffs->hx);
    HYPRE_IJVectorInitialize(coeffs->hy);
  }

  /* Set stat flag */

  coeffs->matrix_state = 1;
}

/*----------------------------------------------------------------------------
 * Release HYPRE ParCSR matrix coefficients.
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
_release_coeffs_ij(cs_matrix_t  *matrix)
{
  cs_matrix_coeffs_hypre_t  *coeffs = matrix->coeffs;

  if (matrix->coeffs != NULL) {

    if (coeffs->matrix_state > 0) {
      HYPRE_IJMatrixDestroy(coeffs->hm);

      HYPRE_IJVectorDestroy(coeffs->hx);
      HYPRE_IJVectorDestroy(coeffs->hy);

      coeffs->matrix_state = 0;
    }
  }
}

/*----------------------------------------------------------------------------
 * Release HYPRE ParCSR matrix coefficients.
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
_destroy_coeffs_ij(cs_matrix_t  *matrix)
{
  cs_matrix_coeffs_hypre_t  *coeffs = matrix->coeffs;

  if (matrix->coeffs != NULL) {

    if (coeffs->matrix_state > 0) {
      HYPRE_IJMatrixDestroy(coeffs->hm);
      coeffs->matrix_state = 0;
    }

    BFT_FREE(coeffs);
  }
}

/*----------------------------------------------------------------------------
 * Copy diagonal of HYPRE ParCSR matrix.
 *
 * parameters:
 *   matrix <-- pointer to matrix structure
 *   da     --> diagonal (pre-allocated, size: n_rows)
 *----------------------------------------------------------------------------*/

static void
_copy_diagonal_ij(const cs_matrix_t  *matrix,
                  cs_real_t          *restrict da)
{
  const cs_matrix_coeffs_hypre_t  *coeffs = matrix->coeffs;

  const HYPRE_BigInt n_rows = matrix->n_rows;

  HYPRE_BigInt ilower, iupper, jlower, jupper;
  HYPRE_IJMatrixGetLocalRange(coeffs->hm,
                              &ilower,
                              &iupper,
                              &jlower,
                              &jupper);

  HYPRE_BigInt *rows, *cols;
  HYPRE_Int *n_rcols;
  HYPRE_Real *aij;

  BFT_MALLOC(n_rcols, n_rows, HYPRE_Int);
  BFT_MALLOC(rows, n_rows, HYPRE_BigInt);
  BFT_MALLOC(cols, n_rows, HYPRE_BigInt);

  if (sizeof(HYPRE_Real) == sizeof(cs_real_t))
    aij = da;
  else
    BFT_MALLOC(aij, n_rows, HYPRE_Real);

# pragma omp parallel for  if(n_rows > CS_THR_MIN)
  for (HYPRE_BigInt ii = 0; ii < n_rows; ii++) {
    n_rcols[ii] = 1;
    rows[ii] = ilower + ii;;
    cols[ii] = ilower + ii;;
  }

  HYPRE_IJMatrixGetValues(coeffs->hm,
                          matrix->n_rows,
                          n_rcols,
                          rows,
                          cols,
                          aij);

  if (aij != da) {
    for (HYPRE_BigInt ii = 0; ii < n_rows; ii++) {
      da[ii] = aij[ii];;
    }
    BFT_FREE(aij);
  }

  BFT_FREE(cols);
  BFT_FREE(rows);
  BFT_FREE(n_rcols);
}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with HYPRE matrix
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
_mat_vec_p_parcsr(const cs_matrix_t  *matrix,
                  bool                exclude_diag,
                  bool                sync,
                  cs_real_t          *restrict x,
                  cs_real_t          *restrict y)
{
  assert(exclude_diag == false);

  const cs_lnum_t  n_rows = matrix->n_rows;
  const cs_matrix_coeffs_hypre_t  *coeffs = matrix->coeffs;

  /* Get pointers to structures through coefficients,
     and copy input values */

  HYPRE_ParCSRMatrix p_a;
  HYPRE_IJMatrixGetObject(coeffs->hm, (void **)&p_a);

  HYPRE_Real *_t = NULL;

  if (sizeof(cs_real_t) == sizeof(HYPRE_Real)) {
    HYPRE_IJVectorSetValues(coeffs->hx, n_rows, NULL, x);
  }
  else {
    BFT_MALLOC(_t, n_rows, HYPRE_Real);
    for (HYPRE_BigInt ii = 0; ii < n_rows; ii++) {
      _t[ii] = x[ii];;
    }
    HYPRE_IJVectorSetValues(coeffs->hx, n_rows, NULL, _t);
  }
  if (sync)
    HYPRE_IJVectorAssemble(coeffs->hx);

  HYPRE_ParVector p_x, p_y;
  HYPRE_IJVectorGetObject(coeffs->hx, (void **)&p_x);
  HYPRE_IJVectorGetObject(coeffs->hy, (void **)&p_y);

  /* SpMv operation */

  HYPRE_ParCSRMatrixMatvec(1.0, p_a, p_x, 0, p_y);

  /* Copy data back */

  if (sizeof(cs_real_t) == sizeof(HYPRE_Real)) {
    HYPRE_IJVectorGetValues(coeffs->hy, n_rows, NULL, y);
  }
  else {
    HYPRE_IJVectorGetValues(coeffs->hy, n_rows, NULL, _t);
    for (HYPRE_BigInt ii = 0; ii < n_rows; ii++) {
      y[ii] = _t[ii];
    }
    BFT_FREE(_t);
  }

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Switch matrix type to hypre.
 *
 * This releases previous coefficients if present, so should be called
 * just after matrix creation, before assigning coefficients.
 *
 * \param[in, out]  matrix  pointer to matrix structure
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_set_type_hypre(cs_matrix_t  *matrix)
{
  matrix->type = CS_MATRIX_N_BUILTIN_TYPES;

  matrix->type_name = _hypre_ij_type_name;
  matrix->type_fname = _hypre_ij_type_fullname;

  /* Release previous coefficients if present */

  if (matrix->coeffs != NULL)
    matrix->destroy_coefficients(matrix);

  cs_matrix_coeffs_hypre_t  *coeffs;
  BFT_MALLOC(coeffs, 1, cs_matrix_coeffs_hypre_t);
  memset(coeffs, 0, sizeof(HYPRE_IJMatrix));
  coeffs->matrix_state = 0;

  matrix->coeffs = coeffs;

  /* Set function pointers here */

  matrix->set_coefficients = _set_coeffs_ij;
  matrix->release_coefficients = _release_coeffs_ij;
  matrix->destroy_coefficients = _destroy_coeffs_ij;
  matrix->get_diagonal = NULL;

  for (int i = 0; i < CS_MATRIX_N_FILL_TYPES; i++) {
    matrix->vector_multiply[i][0] = NULL;
    matrix->vector_multiply[i][1] = NULL;
  }

  cs_matrix_fill_type_t mft[] = {CS_MATRIX_SCALAR,
                                 CS_MATRIX_SCALAR_SYM};

  for (int i = 0; i < 2; i++)
    matrix->vector_multiply[mft[i]][0] = _mat_vec_p_parcsr;

  matrix->copy_diagonal = _copy_diagonal_ij;

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
