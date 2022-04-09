/*============================================================================
 * Sparse Matrix Representation and Operations.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

/*!
 * Notes:
 *
 * The aim of these structures and associated functions is multiple:
 *
 * - Provide an "opaque" matrix object for linear solvers, allowing possible
 *   choice of the matrix type based on run-time tuning at code initialization
 *   (depending on matrix size, architecture, and compiler, the most efficient
 *   structure for matrix.vector products may vary).
 *
 * - Provide at least a CSR matrix structure in addition to the "native"
 *   matrix structure, as this may allow us to leverage existing libraries.
 *
 * - Provide C interface functions required for interfacing with
 *   external libraries.
 *
 * Choice of a given matrix structure may be strongly related to the choice
 * of linear solver or smoother:
 *
 * Though many external libraries assume a classical CSR or block-CSR
 * structure, separate storage of diagonal elements (a.k.a. MSR format)
 * is very useful for Jacobi and Gauss-Seidel solvers and Jacobi or
 * polynomial preconditioners, in which the diagonal elements need to
 * be accessed and/or multiplied separately, so we select that approach.
 *
 * In distributed parallel mode, some matrix entries reference distant
 * rows (and ghost values). For a given local row i and column j,
 * if j > n_local_rows, a_ij references a \em distant value. For the
 * native (graph edge-based) format, symmetric storage is used, but for other
 * (matrix row-based) formats, only local rows are stored (i.e.*
 * a_ij is not stored if i is not in the local row range), so
 * the matrix storage is locally rectangular, even for a logically square
 * matrix.
 *
 * An additional optimization (also used in some external libraries) is to
 * store the matrix entries referencing distant rows (and ghost values)
 * separately, to simplify handling of computation-communication overlap.
 * The "distributed" format is thus a variation of the MSR format with
 * separate distant elements.
 *
 * The specific access requirements of Gauss-Seidel solvers and smoothers
 * lead us to only consider the MSR format for their implementation.
 * When requesting a Gauss-Seidel solver or smoother for another storage
 * format, a Jacobi solver or smoother may be substituted.
 */

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
#include "cs_sort.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_matrix.h"
#include "cs_matrix_priv.h"
#include "cs_matrix_spmv.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*! \file cs_matrix.c
 *
 * \brief Sparse Matrix Representation and Operations.
 *
 * Please refer to the
 * <a href="../../theory.pdf#matrix"><b>matrix</b></a> section of the
 * theory guide for more informations.
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

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Short names for matrix types */

static const char  *_matrix_type_name[] = {N_("native"),
                                           N_("CSR"),
                                           N_("MSR"),
                                           N_("distributed"),
                                           N_("external")};

/* Full names for matrix types */

static const char
*_matrix_type_fullname[] = {N_("diagonal + faces"),
                            N_("Compressed Sparse Row"),
                            N_("Modified Compressed Sparse Row"),
                            N_("Distributed (D+E+H)"),
                            N_("External")};

/* Fill type names for matrices */

const char  *cs_matrix_fill_type_name[] = {"CS_MATRIX_SCALAR",
                                           "CS_MATRIX_SCALAR_SYM",
                                           "CS_MATRIX_BLOCK_D",
                                           "CS_MATRIX_BLOCK_D_66",
                                           "CS_MATRIX_BLOCK_D_SYM",
                                           "CS_MATRIX_BLOCK"};

/*! Operation type type names for partial SpMV functions */

const char  *cs_matrix_spmv_type_name[] = {"y ← A.x",
                                           "y ← (A-D).x"};

#if defined (HAVE_MKL)

static char _no_exclude_diag_error_str[]
  = N_("Matrix product variant using function %s\n"
       "does not handle case with excluded diagonal.");

#endif

/* Tuning parameters */

cs_lnum_t _base_assembler_thr_min = 128;

/*============================================================================
 * Private function definitions
- *============================================================================*/

/*----------------------------------------------------------------------------
 * Set matrix fill metadata.
 *
 * parameters:
 *   matrix                <-> pointer to matrix structure
 *   symmetric             <-- indicates if matrix coefficients are symmetric
 *   diag_block_size       <-- block sizes for diagonal
 *   extra_diag_block_size <-- block sizes for extra diagonal
 *----------------------------------------------------------------------------*/

static void
_set_fill_info(cs_matrix_t   *matrix,
               bool           symmetric,
               cs_lnum_t      diag_block_size,
               cs_lnum_t      extra_diag_block_size)
{
  matrix->symmetric = symmetric;

  matrix->db_size = diag_block_size;
  matrix->eb_size = extra_diag_block_size;

  /* Set fill type */

  matrix->fill_type = cs_matrix_get_fill_type(symmetric,
                                              diag_block_size,
                                              extra_diag_block_size);
}

/*----------------------------------------------------------------------------
 * Clear matrix fill metadata.
 *
 * parameters:
 *   matrix <-> pointer to matrix structure
 *----------------------------------------------------------------------------*/

static void
_clear_fill_info(cs_matrix_t  *matrix)
{
  matrix->symmetric = false;

  matrix->db_size = 0;
  matrix->eb_size = 0;

  matrix->fill_type = CS_MATRIX_N_FILL_TYPES;
}

/*----------------------------------------------------------------------------
 * Start synchronization of ghost values prior to matrix.vector product
 *
 * parameters:
 *   matrix        <-- pointer to matrix structure
 *   x             <-> multiplying vector values (ghost values updated)
 *
 * returns:
 *   halo state to use for synchronisation finalisation.
 *----------------------------------------------------------------------------*/

static cs_halo_state_t *
_pre_vector_multiply_sync_x_start(const cs_matrix_t   *matrix,
                                  cs_real_t            x[restrict])
{
 cs_halo_state_t *hs = NULL;

  if (matrix->halo != NULL) {

    hs = cs_halo_state_get_default();

    /* Non-blocked version */

    cs_halo_sync_pack(matrix->halo,
                      CS_HALO_STANDARD,
                      CS_REAL_TYPE,
                      matrix->db_size,
                      x,
                      NULL,
                      hs);

    cs_halo_sync_start(matrix->halo, x, hs);

  }

  return hs;
}

/*----------------------------------------------------------------------------
 * Synchronize ghost values prior to matrix.vector product
 *
 * parameters:
 *   matrix        <-- pointer to matrix structure
 *   x             <-> multiplying vector values (ghost values updated)
 *----------------------------------------------------------------------------*/

static void
_pre_vector_multiply_sync_x_end(const cs_matrix_t   *matrix,
                                cs_halo_state_t     *hs,
                                cs_real_t            x[restrict])
{
  if (hs != NULL) {

    assert(matrix->halo != NULL);

    cs_halo_sync_wait(matrix->halo, x, hs);

    /* Synchronize periodic values */

#if !defined(_CS_UNIT_MATRIX_TEST) /* unit tests do not link with full library */

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

/*----------------------------------------------------------------------------
 * Synchronize ghost values prior to matrix.vector product
 *
 * parameters:
 *   matrix        <-- pointer to matrix structure
 *   x             <-> multiplying vector values (ghost values updated)
 *----------------------------------------------------------------------------*/

void
cs_matrix_pre_vector_multiply_sync(const cs_matrix_t   *matrix,
                                   cs_real_t           *x)
{
  if (matrix->halo != NULL) {
    cs_halo_state_t *hs = _pre_vector_multiply_sync_x_start(matrix, x);
    _pre_vector_multiply_sync_x_end(matrix, hs, x);
  }
}

/*----------------------------------------------------------------------------
 * Create native matrix structure.
 *
 * Note that the structure created maps to the given existing
 * face -> cell connectivity array, so it must be destroyed before this
 * array (usually the code's main face -> cell structure) is freed.
 *
 * parameters:
 *   n_rows      <-- number of local rows
 *   n_cols_ext  <-- number of local + ghost columns
 *   n_edges     <-- local number of graph edges
 *   edges       <-- edges (symmetric row <-> column) connectivity
 *
 * returns:
 *   pointer to allocated native matrix structure.
 *----------------------------------------------------------------------------*/

static cs_matrix_struct_native_t *
_create_struct_native(cs_lnum_t        n_rows,
                      cs_lnum_t        n_cols_ext,
                      cs_lnum_t        n_edges,
                      const cs_lnum_t  edges[][2])
{
  cs_matrix_struct_native_t  *ms;

  /* Allocate and map */

  BFT_MALLOC(ms, 1, cs_matrix_struct_native_t);

  /* Allocate and map */

  ms->n_rows = n_rows;
  ms->n_cols_ext = n_cols_ext;
  ms->n_edges = n_edges;

  ms->edges = edges;

  return ms;
}

/*----------------------------------------------------------------------------
 * Destroy native matrix structure.
 *
 * parameters:
 *   ms  <->  pointer to native matrix structure pointer
 *----------------------------------------------------------------------------*/

static void
_destroy_struct_native(void  **ms)
{
  if (ms != NULL && *ms !=NULL) {
    cs_matrix_struct_native_t  *_ms = *ms;

    BFT_FREE(_ms);

    *ms= NULL;
  }
}

/*----------------------------------------------------------------------------
 * Set Native matrix coefficients.
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
_set_coeffs_native(cs_matrix_t        *matrix,
                   bool                symmetric,
                   bool                copy,
                   cs_lnum_t           n_edges,
                   const cs_lnum_t     edges[restrict][2],
                   const cs_real_t     da[restrict],
                   const cs_real_t     xa[restrict])
{
  CS_UNUSED(n_edges);
  CS_UNUSED(edges);

  cs_matrix_coeff_dist_t  *mc = matrix->coeffs;
  const cs_matrix_struct_native_t  *ms = matrix->structure;

  mc->symmetric = symmetric;
  mc->db_size = matrix->db_size;
  mc->eb_size = matrix->eb_size;

  CS_FREE(mc->_d_val);
  CS_FREE(mc->_e_val);

  /* Map or copy values */

  if (da != NULL) {

    if (copy) {
      const cs_lnum_t b_size_2 = mc->db_size*mc->db_size;
      CS_MALLOC_HD(mc->_d_val, b_size_2*ms->n_rows, cs_real_t,
                   matrix->alloc_mode);
      memcpy(mc->_d_val, da, b_size_2 * sizeof(cs_real_t) * ms->n_rows);
      mc->d_val = mc->_d_val;
    }
    else
      mc->d_val = da;

  }
  else {
    mc->d_val = NULL;
  }

  if (xa != NULL) {

    size_t xa_n_vals = ms->n_edges;
    if (! symmetric)
      xa_n_vals *= 2;

    if (copy) {
      const cs_lnum_t b_size_2 = mc->eb_size*mc->eb_size;
      CS_MALLOC_HD(mc->_e_val, b_size_2*xa_n_vals, cs_real_t,
                   matrix->alloc_mode);
      memcpy(mc->_e_val, xa, b_size_2*xa_n_vals*sizeof(cs_real_t));
      mc->e_val = mc->_e_val;
    }
    else
      mc->e_val = xa;

  }
}

/*----------------------------------------------------------------------------
 * Copy diagonal of native, MSR, or distributed matrix.
 *
 * parameters:
 *   matrix <-- pointer to matrix structure
 *   da     --> diagonal (pre-allocated, size: n_rows)
 *----------------------------------------------------------------------------*/

static void
_copy_diagonal_separate(const cs_matrix_t  *matrix,
                        cs_real_t           da[restrict])
{
  const cs_real_t *_d_val = NULL;
  if (matrix->type == CS_MATRIX_NATIVE) {
    const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;
    _d_val = mc->d_val;
  }
  else if (   matrix->type == CS_MATRIX_MSR
           || matrix->type == CS_MATRIX_DIST) {
    const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;
    _d_val = mc->d_val;
  }
  const cs_lnum_t  n_rows = matrix->n_rows;

  /* Unblocked version */

  if (matrix->db_size == 1) {

    if (_d_val != NULL) {
#     pragma omp parallel for  if(n_rows > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        da[ii] = _d_val[ii];
    }
    else {
#     pragma omp parallel for  if(n_rows > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        da[ii] = 0.0;
    }

  }

  /* Blocked version */

  else {

    const cs_lnum_t db_size = matrix->db_size;
    const cs_lnum_t db_size_2 = db_size * db_size;

    if (_d_val != NULL) {
#     pragma omp parallel for  if(n_rows*db_size > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
        for (cs_lnum_t jj = 0; jj < db_size; jj++)
          da[ii*db_size + jj] = _d_val[ii*db_size_2 + jj*db_size + jj];
      }
    }
    else {
#     pragma omp parallel for  if(n_rows*db_size > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_rows*db_size; ii++)
        da[ii] = 0.0;
    }
  }
}

/*----------------------------------------------------------------------------
 * Destroy a CSR matrix structure.
 *
 * parameters:
 *   ms  <->  pointer to CSR matrix structure pointer
 *----------------------------------------------------------------------------*/

static void
_destroy_struct_csr(void  **ms)
{
  if (ms != NULL && *ms !=NULL) {
    cs_matrix_struct_csr_t  *_ms = *ms;

    BFT_FREE(_ms->_row_index);
    BFT_FREE(_ms->_col_id);
    BFT_FREE(_ms);

    *ms= NULL;
  }
}

/*----------------------------------------------------------------------------
 * Order and compact CSR matrix structure.
 *
 * parameters:
 *   ms         <-> CSR matrix structure
 *   alloc_mode <-- allocation mode
 *----------------------------------------------------------------------------*/

static void
_compact_struct_csr(cs_matrix_struct_csr_t  *ms,
                    cs_alloc_mode_t          alloc_mode)
{
  /* Sort line elements by column id (for better access patterns) */

  ms->direct_assembly = cs_sort_indexed(ms->n_rows,
                                        ms->_row_index,
                                        ms->_col_id);

  /* Compact elements if necessary */

  if (ms->direct_assembly == false) {

    cs_lnum_t *tmp_row_index = NULL;
    cs_lnum_t  kk = 0;

    BFT_MALLOC(tmp_row_index, ms->n_rows+1, cs_lnum_t);
    memcpy(tmp_row_index, ms->_row_index, (ms->n_rows+1)*sizeof(cs_lnum_t));

    kk = 0;

    for (cs_lnum_t ii = 0; ii < ms->n_rows; ii++) {
      cs_lnum_t *col_id = ms->_col_id + ms->_row_index[ii];
      cs_lnum_t n_cols = ms->_row_index[ii+1] - ms->_row_index[ii];
      cs_lnum_t col_id_prev = -1;
      ms->_row_index[ii] = kk;
      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        if (col_id_prev != col_id[jj]) {
          ms->_col_id[kk++] = col_id[jj];
          col_id_prev = col_id[jj];
        }
      }
    }
    ms->_row_index[ms->n_rows] = kk;

    assert(ms->_row_index[ms->n_rows] < tmp_row_index[ms->n_rows]);

    BFT_FREE(tmp_row_index);
    CS_REALLOC_HD(ms->_col_id,
                  (ms->_row_index[ms->n_rows]),
                  cs_lnum_t,
                  alloc_mode);

  }

  ms->row_index = ms->_row_index;
  ms->col_id = ms->_col_id;
}

/*----------------------------------------------------------------------------
 * Create a CSR matrix structure from a native matrix stucture.
 *
 * Note that the structure created maps global cell numbers to the given
 * existing face -> cell connectivity array, so it must be destroyed before
 * this array (usually the code's global cell numbering) is freed.
 *
 * parameters:
 *   alloc_mode  <-- allocation mode
 *   n_rows      <-- number of local rows
 *   n_cols_ext  <-- number of local + ghost columns
 *   n_edges     <-- local number of graph edges
 *   edges       <-- edges (symmetric row <-> column) connectivity
 *
 * returns:
 *   pointer to allocated CSR matrix structure.
 *----------------------------------------------------------------------------*/

static cs_matrix_struct_csr_t *
_create_struct_csr(cs_alloc_mode_t     alloc_mode,
                   cs_lnum_t           n_rows,
                   cs_lnum_t           n_cols_ext,
                   cs_lnum_t           n_edges,
                   const cs_lnum_2_t  *edges)
{
  cs_lnum_t ii, jj, face_id;
  const cs_lnum_t *restrict face_cel_p;

  cs_lnum_t  *ccount = NULL;

  cs_matrix_struct_csr_t  *ms;

  /* Allocate and map */

  BFT_MALLOC(ms, 1, cs_matrix_struct_csr_t);

  ms->n_rows = n_rows;
  ms->n_cols_ext = n_cols_ext;

  ms->direct_assembly = true;
  ms->have_diag = true;

  CS_MALLOC_HD(ms->_row_index, ms->n_rows + 1, cs_lnum_t, alloc_mode);
  ms->row_index = NULL;

  /* Count number of nonzero elements per row */

  BFT_MALLOC(ccount, ms->n_rows, cs_lnum_t);

  for (ii = 0; ii < ms->n_rows; ii++)  /* count starting with diagonal terms */
    ccount[ii] = 1;

  if (edges != NULL) {

    face_cel_p = (const cs_lnum_t *restrict)edges;

    for (face_id = 0; face_id < n_edges; face_id++) {
      ii = *face_cel_p++;
      jj = *face_cel_p++;
      if (ii < ms->n_rows)
        ccount[ii] += 1;
      if (jj < ms->n_rows)
        ccount[jj] += 1;
    }

  } /* if (edges != NULL) */

  ms->_row_index[0] = 0;
  for (ii = 0; ii < ms->n_rows; ii++) {
    ms->_row_index[ii+1] = ms->_row_index[ii] + ccount[ii];
    ccount[ii] = 1; /* pre-count for diagonal terms */
  }

  /* Build structure */

  CS_MALLOC_HD(ms->_col_id, (ms->_row_index[ms->n_rows]), cs_lnum_t, alloc_mode);
  ms->col_id = NULL;

  for (ii = 0; ii < ms->n_rows; ii++) {    /* diagonal terms */
    ms->_col_id[ms->_row_index[ii]] = ii;
  }

  if (edges != NULL) {                   /* non-diagonal terms */

    face_cel_p = (const cs_lnum_t *restrict)edges;

    for (face_id = 0; face_id < n_edges; face_id++) {
      ii = *face_cel_p++;
      jj = *face_cel_p++;
      if (ii < ms->n_rows) {
        ms->_col_id[ms->_row_index[ii] + ccount[ii]] = jj;
        ccount[ii] += 1;
      }
      if (jj < ms->n_rows) {
        ms->_col_id[ms->_row_index[jj] + ccount[jj]] = ii;
        ccount[jj] += 1;
      }
    }

  } /* if (edges != NULL) */

  BFT_FREE(ccount);

  /* Compact and finalize indexes */

  _compact_struct_csr(ms, alloc_mode);

  return ms;
}

/*----------------------------------------------------------------------------
 * Initialize a CSR structure.
 *
 * parameters:
 *   ms          <-> pointer to csr structure
 *   n_rows      <-- number of local rows
 *   n_cols_ext  <-- number of local + ghost columns
 *----------------------------------------------------------------------------*/

static void
_init_struct_csr(cs_matrix_struct_csr_t  *ms,
                 cs_lnum_t                n_rows,
                 cs_lnum_t                n_cols_ext)
{
  ms->n_rows = n_rows;
  ms->n_cols_ext = n_cols_ext;

  ms->direct_assembly = true;
  ms->have_diag = false;

  /* Count number of nonzero elements per row */

  ms->row_index = NULL;
  ms->col_id = NULL;

  ms->_row_index = NULL;
  ms->_col_id = NULL;
}

/*----------------------------------------------------------------------------
 * Initialize a CSR matrix structure from an index and an array related
 * to column id
 *
 * parameters:
 *   ms         <-- matrix structure to initialize
 *   have_diag  <-- indicates if the diagonal structure contains nonzeroes
 *   transfer   <-- transfer property of row_index and col_id
 *                  if true, map them otherwise
 *   ordered    <-- indicates if row entries are already ordered
 *   n_rows     <-- local number of rows
 *   n_cols_ext <-- local number of columns + ghosts
 *   row_index  <-- pointer to index on rows
 *   col_id     <-> pointer to array of column ids related to the row index
 *
 * returns:
 *    a pointer to a created CSR matrix structure
 *----------------------------------------------------------------------------*/

static void
_init_struct_csr_from_csr(cs_matrix_struct_csr_t  *ms,
                          bool                     have_diag,
                          bool                     transfer,
                          bool                     ordered,
                          cs_lnum_t                n_rows,
                          cs_lnum_t                n_cols_ext,
                          cs_lnum_t             **row_index,
                          cs_lnum_t             **col_id)
{
  cs_lnum_t  *_row_index = *row_index;
  cs_lnum_t  *_col_id = *col_id;

  /* Allocate and map */

  ms->n_rows = n_rows;
  ms->n_cols_ext = n_cols_ext;

  ms->direct_assembly = false; /* not relevant here */
  ms->have_diag = have_diag;

  ms->row_index = _row_index;
  ms->col_id = _col_id;

  ms->_row_index = NULL;
  ms->_col_id = NULL;

  if (transfer == true) {

    ms->_row_index = _row_index;
    ms->_col_id = _col_id;

    *row_index = NULL;
    *col_id = NULL;

    /* Sort line elements by column id (for better access patterns) */

    if (! ordered)
      cs_sort_indexed(ms->n_rows,
                      ms->_row_index,
                      ms->_col_id);

  }
}

/*----------------------------------------------------------------------------
 * Create a CSR matrix structure from an index and an array related
 * to column id
 *
 * parameters:
 *   have_diag  <-- indicates if the diagonal structure contains nonzeroes
 *   transfer   <-- transfer property of row_index and col_id
 *                  if true, map them otherwise
 *   ordered    <-- indicates if row entries are already ordered
 *   n_rows     <-- local number of rows
 *   n_cols_ext <-- local number of columns + ghosts
 *   row_index  <-- pointer to index on rows
 *   col_id     <-> pointer to array of column ids related to the row index
 *
 * returns:
 *    a pointer to a created CSR matrix structure
 *----------------------------------------------------------------------------*/

static cs_matrix_struct_csr_t *
_create_struct_csr_from_csr(bool         have_diag,
                            bool         transfer,
                            bool         ordered,
                            cs_lnum_t    n_rows,
                            cs_lnum_t    n_cols_ext,
                            cs_lnum_t  **row_index,
                            cs_lnum_t  **col_id)
{
  cs_matrix_struct_csr_t  *ms = NULL;

  BFT_MALLOC(ms, 1, cs_matrix_struct_csr_t);

  _init_struct_csr_from_csr(ms,
                            have_diag,
                            transfer,
                            ordered,
                            n_rows,
                            n_cols_ext,
                            row_index,
                            col_id);

  return ms;
}

/*----------------------------------------------------------------------------
 * Create a CSR matrix structure from an index and an array related
 * to column id
 *
 * parameters:
 *   have_diag       <-- indicates if diagonal structure contains nonzeroes
 *   direct_assembly <-- true if each value corresponds to a unique face
 *   n_rows          <- local number of rows
 *   n_cols_ext      <-- local number of columns + ghosts
 *   row_index       <-- index on rows
 *   col_id          <-> array of colum ids related to the row index
 *
 * returns:
 *    a pointer to a created CSR matrix structure
 *----------------------------------------------------------------------------*/

static cs_matrix_struct_csr_t *
_create_struct_csr_from_shared(bool              have_diag,
                               bool              direct_assembly,
                               cs_lnum_t         n_rows,
                               cs_lnum_t         n_cols_ext,
                               const cs_lnum_t  *row_index,
                               const cs_lnum_t  *col_id)
{
  cs_matrix_struct_csr_t  *ms = NULL;

  /* Allocate and map */

  BFT_MALLOC(ms, 1, cs_matrix_struct_csr_t);

  ms->n_rows = n_rows;
  ms->n_cols_ext = n_cols_ext;

  ms->direct_assembly = direct_assembly;
  ms->have_diag = have_diag;

  ms->row_index = row_index;
  ms->col_id = col_id;

  ms->_row_index = NULL;
  ms->_col_id = NULL;

  return ms;
}

/*----------------------------------------------------------------------------
 * Create a CSR matrix structure from the restriction to local rank of
 * another CSR matrix structure.
 *
 * parameters:
 *   src <-- base matrix structure
 *
 * returns:
 *    a pointer to a created CSR matrix structure
 *----------------------------------------------------------------------------*/

static cs_matrix_struct_csr_t *
_create_struct_csr_from_restrict_local(const cs_matrix_struct_csr_t  *src)
{
  cs_matrix_struct_csr_t  *ms = NULL;

  /* Allocate and map */

  BFT_MALLOC(ms, 1, cs_matrix_struct_csr_t);

  const cs_lnum_t n_rows = src->n_rows;

  ms->n_rows = n_rows;
  ms->n_cols_ext = n_rows;

  ms->direct_assembly = src->direct_assembly;
  ms->have_diag = src->have_diag;

  cs_alloc_mode_t amode = cs_check_device_ptr(src->row_index);

  CS_MALLOC_HD(ms->_row_index, ms->n_rows+1, cs_lnum_t, amode);
  CS_MALLOC_HD(ms->_col_id, src->row_index[ms->n_rows], cs_lnum_t, amode);

  ms->_row_index[0] = 0;

  cs_lnum_t k = 0;

  const cs_lnum_t *col_id_s = src->col_id;
  cs_lnum_t *col_id_d = ms->_col_id;

  for (cs_lnum_t i = 0; i < n_rows; i++) {
    const cs_lnum_t s_id = src->row_index[i];
    const cs_lnum_t e_id = src->row_index[i+1];
    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t c_id = col_id_s[j];
      if (c_id < n_rows) {
        col_id_d[k] = c_id;
        k += 1;
      }
    }
    ms->_row_index[i+1] = k;
  }

  CS_REALLOC_HD(ms->_col_id, ms->_row_index[n_rows], cs_lnum_t, amode);

  ms->row_index = ms->_row_index;
  ms->col_id = ms->_col_id;

  return ms;
}

/*----------------------------------------------------------------------------
 * Destroy CSR matrix coefficients.
 *
 * parameters:
 *   m  <->  pointer to matrix structure
 *----------------------------------------------------------------------------*/

static void
_destroy_coeff_csr(cs_matrix_t  *m)
{
  if (m->coeffs != NULL) {
    cs_matrix_coeff_csr_t  *mc = m->coeffs;

    BFT_FREE(mc->_val);
    BFT_FREE(mc->_d_val);

    BFT_FREE(m->coeffs);
  }
}

/*----------------------------------------------------------------------------
 * Create CSR matrix coefficients.
 *
 * returns:
 *   pointer to allocated CSR coefficients structure.
 *----------------------------------------------------------------------------*/

static cs_matrix_coeff_csr_t *
_create_coeff_csr(void)
{
  cs_matrix_coeff_csr_t  *mc = NULL;

  /* Allocate */

  BFT_MALLOC(mc, 1, cs_matrix_coeff_csr_t);

  /* Initialize */

  mc->val = NULL;
  mc->_val = NULL;

  mc->d_val = NULL;
  mc->_d_val = NULL;

  return mc;
}

/*----------------------------------------------------------------------------
 * Set CSR matrix coefficients to zero.
 *
 * The coefficients should already be allocated.
 *
 * Use of this function is preferrable to a simple loop, as its
 * threading behavior should be consistent with SpMW in NUMA cases.
 *
 * parameters:
 *   ms      <-> pointer to matrix structure
 *   eb_size <-> associated block size
 *   val     <-- value
 *----------------------------------------------------------------------------*/

static void
_zero_coeffs_csr(const cs_matrix_struct_csr_t  *ms,
                 const cs_lnum_t                eb_size,
                 cs_real_t                      val[])
{
  const cs_lnum_t  n_rows = ms->n_rows;

  if (eb_size == 1) {
#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
      const cs_lnum_t  n_cols = ms->row_index[ii+1] - ms->row_index[ii];
      cs_real_t  *m_row = val + ms->row_index[ii];
      for (cs_lnum_t jj = 0; jj < n_cols; jj++)
        m_row[jj] = 0.0;
    }
  }
  else {
    const cs_lnum_t eb_size_2 = eb_size * eb_size;
#   pragma omp parallel for  if(n_rows*eb_size > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
      const cs_lnum_t  n_cols = ms->row_index[ii+1] - ms->row_index[ii];
      cs_real_t  *m_row = val + ms->row_index[ii]*eb_size_2;
      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        for (cs_lnum_t kk = 0; kk < eb_size_2; kk++)
          m_row[jj*eb_size_2 + kk] = 0.0;
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Set CSR extradiagonal matrix coefficients for the case where direct
 * assignment is possible (i.e. when there are no multiple contributions
 * to a given coefficient).
 *
 * parameters:
 *   matrix      <-- pointer to matrix structure
 *   symmetric   <-- indicates if extradiagonal values are symmetric
 *   n_edges     <-- local number of graph edges
 *   edges       <-- edges (symmetric row <-> column) connectivity
 *   xa          <-- extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_xa_coeffs_csr_direct(cs_matrix_t        *matrix,
                          bool                symmetric,
                          cs_lnum_t           n_edges,
                          const cs_lnum_2_t  *edges,
                          const cs_real_t     xa[restrict])
{
  cs_lnum_t  ii, jj, face_id;
  cs_matrix_coeff_csr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  /* Copy extra-diagonal values */

  assert(edges != NULL);

  const cs_lnum_t *restrict edges_p
    = (const cs_lnum_t *restrict)(edges);

  if (symmetric == false) {

    for (face_id = 0; face_id < n_edges; face_id++) {
      cs_lnum_t kk, ll;
      ii = *edges_p++;
      jj = *edges_p++;
      if (ii < ms->n_rows) {
        for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
        mc->_val[kk] = xa[2*face_id];
      }
      if (jj < ms->n_rows) {
        for (ll = ms->row_index[jj]; ms->col_id[ll] != ii; ll++);
        mc->_val[ll] = xa[2*face_id + 1];
      }
    }

  }
  else { /* if symmetric == true */

    for (face_id = 0; face_id < n_edges; face_id++) {
      cs_lnum_t kk, ll;
      ii = *edges_p++;
      jj = *edges_p++;
      if (ii < ms->n_rows) {
        for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
        mc->_val[kk] = xa[face_id];
      }
      if (jj < ms->n_rows) {
        for (ll = ms->row_index[jj]; ms->col_id[ll] != ii; ll++);
        mc->_val[ll] = xa[face_id];
      }

    }

  } /* end of condition on coefficients symmetry */

}

/*----------------------------------------------------------------------------
 * Set CSR extradiagonal matrix coefficients for the case where there are
 * multiple contributions to a given coefficient).
 *
 * The matrix coefficients should have been initialized (i.e. set to 0)
 * some before using this function.
 *
 * parameters:
 *   matrix      <-- pointer to matrix structure
 *   symmetric   <-- indicates if extradiagonal values are symmetric
 *   n_edges     <-- local number of graph edges
 *   edges       <-- edges (symmetric row <-> column) connectivity
 *   xa          <-- extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_xa_coeffs_csr_increment(cs_matrix_t        *matrix,
                             bool                symmetric,
                             cs_lnum_t           n_edges,
                             const cs_lnum_2_t   edges[restrict],
                             const cs_real_t     xa[restrict])
{
  cs_lnum_t  ii, jj, face_id;
  cs_matrix_coeff_csr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  /* Copy extra-diagonal values */

  assert(edges != NULL);

  const cs_lnum_t *restrict edges_p
    = (const cs_lnum_t *restrict)(edges);

  if (symmetric == false) {

    for (face_id = 0; face_id < n_edges; face_id++) {
      cs_lnum_t kk, ll;
      ii = *edges_p++;
      jj = *edges_p++;
      if (ii < ms->n_rows) {
        for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
        mc->_val[kk] += xa[2*face_id];
      }
      if (jj < ms->n_rows) {
        for (ll = ms->row_index[jj]; ms->col_id[ll] != ii; ll++);
        mc->_val[ll] += xa[2*face_id + 1];
      }
    }

  }
  else { /* if symmetric == true */

    for (face_id = 0; face_id < n_edges; face_id++) {
      cs_lnum_t kk, ll;
      ii = *edges_p++;
      jj = *edges_p++;
      if (ii < ms->n_rows) {
        for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
        mc->_val[kk] += xa[face_id];
      }
      if (jj < ms->n_rows) {
        for (ll = ms->row_index[jj]; ms->col_id[ll] != ii; ll++);
        mc->_val[ll] += xa[face_id];
      }

    }

  } /* end of condition on coefficients symmetry */

}

/*----------------------------------------------------------------------------
 * Set CSR matrix coefficients.
 *
 * parameters:
 *   matrix      <-> pointer to matrix structure
 *   symmetric   <-- indicates if extradiagonal values are symmetric
 *   copy        <-- indicates if coefficients should be copied
 *   n_edges     <-- local number of graph edges
 *   edges       <-- edges (symmetric row <-> column) connectivity
 *   da          <-- diagonal values (NULL if all zero)
 *   xa          <-- extradiagonal values (NULL if all zero)
 *----------------------------------------------------------------------------*/

static void
_set_coeffs_csr(cs_matrix_t      *matrix,
                bool              symmetric,
                bool              copy,
                cs_lnum_t         n_edges,
                const cs_lnum_t   edges[restrict][2],
                const cs_real_t   da[restrict],
                const cs_real_t   xa[restrict])
{
  CS_UNUSED(copy);

  cs_matrix_coeff_csr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  if (mc->_val == NULL)
    CS_MALLOC_HD(mc->_val, ms->row_index[ms->n_rows], cs_real_t,
                 matrix->alloc_mode);
  mc->val = mc->_val;

  /* Initialize coefficients to zero if assembly is incremental */

  if (ms->direct_assembly == false || (n_edges > 0 && xa == NULL))
    _zero_coeffs_csr(ms, matrix->eb_size, mc->_val);

  /* Copy diagonal values */

  if (ms->have_diag == true) {

    if (da != NULL) {
      for (cs_lnum_t ii = 0; ii < ms->n_rows; ii++) {
        cs_lnum_t kk;
        for (kk = ms->row_index[ii]; ms->col_id[kk] != ii; kk++);
        mc->_val[kk] = da[ii];
      }
    }
    else {
      for (cs_lnum_t ii = 0; ii < ms->n_rows; ii++) {
        cs_lnum_t kk;
        for (kk = ms->row_index[ii]; ms->col_id[kk] != ii; kk++);
        mc->_val[kk] = 0.0;
      }
    }

  }

  /* Mark diagonal values as not queried (mc->_d_val not changed) */

  mc->d_val = NULL;

  /* Copy extra-diagonal values */

  if (edges != NULL && xa != NULL) {

    if (ms->direct_assembly == true)
      _set_xa_coeffs_csr_direct(matrix, symmetric, n_edges, edges, xa);
    else
      _set_xa_coeffs_csr_increment(matrix, symmetric, n_edges, edges, xa);

  }

}

/*----------------------------------------------------------------------------
 * Set CSR matrix coefficients provided in MSR form.
 *
 * If da and xa are equal to NULL, then initialize val with zeros.
 *
 * parameters:
 *   matrix           <-> pointer to matrix structure
 *   row_index        <-- MSR row index (0 to n-1)
 *   col_id           <-- MSR column id (0 to n-1)
 *   d_vals           <-- diagonal values (NULL if all zero)
 *   d_vals_transfer  <-- diagonal values whose ownership is transferred
 *                        (NULL or d_vals in, NULL out)
 *   x_vals           <-- extradiagonal values (NULL if all zero)
 *   x_vals_transfer  <-- extradiagonal values whose ownership is transferred
 *                        (NULL or x_vals in, NULL out)
 *----------------------------------------------------------------------------*/

static void
_set_coeffs_csr_from_msr(cs_matrix_t       *matrix,
                         const cs_lnum_t    row_index[],
                         const cs_lnum_t    col_id[],
                         const cs_real_t    d_vals[restrict],
                         cs_real_t        **d_vals_transfer,
                         const cs_real_t    x_vals[restrict],
                         cs_real_t        **x_vals_transfer)
{
  cs_matrix_coeff_csr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  const cs_lnum_t  n_rows = ms->n_rows;

  /* Sanity check */

  if (matrix->db_size > 1 || matrix->eb_size > 1)
    bft_error
      (__FILE__, __LINE__, 0,
       "%s:\n"
       "  case with diagonal block size %ld en extradiagonal block size %ld\n"
       "  not implemented.",
       __func__, (long)matrix->db_size, (long)matrix->eb_size);

  /* Special configuration where ownership is transferred directly */

  /* TODO: we should use metadata or check that the row_index and
     column id values are consistent, which should be true as long
     as columns are ordered in an identical manner */

  if (x_vals_transfer != NULL) {
    if (d_vals == NULL && *x_vals_transfer != NULL) {
      mc->_val = *x_vals_transfer;
      mc->val = mc->_val;
      *x_vals_transfer = NULL;
      return;
    }
  }

  /* Allocate local array */

  if (mc->_val == NULL)
    CS_MALLOC_HD(mc->_val, ms->row_index[ms->n_rows], cs_real_t,
                 matrix->alloc_mode);

  mc->val = mc->_val;

  /* Mark diagonal values as not queried (mc->_d_val not changed) */

  mc->d_val = NULL;

  /* Case with diagonal and extradiagonal values */

  if (d_vals != NULL && x_vals != NULL) {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      const cs_lnum_t  *restrict m_col_id = ms->col_id + ms->row_index[ii];
      cs_real_t  *restrict m_row = mc->_val + ms->row_index[ii];
      cs_lnum_t  n_cols = ms->row_index[ii+1] - ms->row_index[ii];

      const cs_lnum_t  *restrict s_col_id = col_id + row_index[ii];
      const cs_real_t  *restrict s_row = x_vals + ms->row_index[ii];
      cs_lnum_t  n_s_cols = row_index[ii+1] - row_index[ii];

      cs_lnum_t c_id_0 = 0;

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        if (m_col_id[jj] == ii)
          m_row[jj] = d_vals[ii];
        else {
          /* Optimize for ordered case */
          if (m_col_id[jj] == s_col_id[c_id_0]) {
            m_row[jj] = s_row[c_id_0];
            c_id_0++;
          }
          else {
            for (cs_lnum_t kk = c_id_0; kk < n_s_cols; kk++) {
              if (m_col_id[jj] == s_col_id[kk]) {
                m_row[jj] = s_row[kk];
                break;
              }
            }
          }
        }
      }

    }
  }

  /* Case with diagonal values only */

  else if (d_vals != NULL) {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      const cs_lnum_t  *restrict m_col_id = ms->col_id + ms->row_index[ii];
      cs_real_t  *restrict m_row = mc->_val + ms->row_index[ii];
      cs_lnum_t  n_cols = ms->row_index[ii+1] - ms->row_index[ii];

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        if (m_col_id[jj] == ii)
          m_row[jj] = d_vals[ii];
        else
          m_row[jj] = 0.;
      }

    }
  }

  /* Case with null-diagonal */

  else if (x_vals != NULL) {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      const cs_lnum_t  *restrict m_col_id = ms->col_id + ms->row_index[ii];
      cs_real_t  *restrict m_row = mc->_val + ms->row_index[ii];
      cs_lnum_t  n_cols = ms->row_index[ii+1] - ms->row_index[ii];

      const cs_lnum_t  *restrict s_col_id = col_id + row_index[ii];
      const cs_real_t  *restrict s_row = x_vals + ms->row_index[ii];
      cs_lnum_t  n_s_cols = row_index[ii+1] - row_index[ii];

      cs_lnum_t c_id_0 = 0;

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        if (m_col_id[jj] == ii)
          m_row[jj] = 0.;
        else {
          /* Optimize for ordered case */
          if (m_col_id[jj] == s_col_id[c_id_0]) {
            m_row[jj] = s_row[c_id_0];
            c_id_0++;
          }
          else {
            for (cs_lnum_t kk = c_id_0; kk < n_s_cols; kk++) {
              if (m_col_id[jj] == s_col_id[kk]) {
                m_row[jj] = s_row[kk];
                break;
              }
            }
          }
        }
      }

    }

  }

  else
    _zero_coeffs_csr(ms, matrix->eb_size, mc->_val);

  /* Now free transferred arrays */

  if (d_vals_transfer != NULL)
    BFT_FREE(*d_vals_transfer);
  if (x_vals_transfer != NULL)
    BFT_FREE(*x_vals_transfer);
}

/*----------------------------------------------------------------------------
 * Release shared CSR matrix coefficients.
 *
 * parameters:
 *   matrix <-- pointer to matrix structure
 *----------------------------------------------------------------------------*/

static void
_release_coeffs_csr(cs_matrix_t  *matrix)
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
 *   matrix <-- pointer to matrix structure
 *   da     --> diagonal (pre-allocated, size: n_rows)
 *----------------------------------------------------------------------------*/

static void
_copy_diagonal_csr(const cs_matrix_t  *matrix,
                   cs_real_t          *restrict da)
{
  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_csr_t  *mc = matrix->coeffs;
  cs_lnum_t  n_rows = ms->n_rows;

# pragma omp parallel for  if(n_rows > CS_THR_MIN)
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get matrix diagonal values for CSR matrix.
 *
 * In case of matrixes with block diagonal coefficients, a pointer to
 * the complete block diagonal is returned.
 *
 * \param[in]  matrix  pointer to matrix structure
 *
 * \return  pointer to matrix diagonal array
 */
/*----------------------------------------------------------------------------*/

static const cs_real_t *
_get_diagonal_csr(const cs_matrix_t  *matrix)
{
  const cs_real_t  *diag = NULL;

  cs_matrix_coeff_csr_t *mc = matrix->coeffs;
  assert(matrix->db_size == 1);
  if (mc->_d_val == NULL)
    CS_MALLOC_HD(mc->_d_val, matrix->n_rows, cs_real_t, matrix->alloc_mode);
  if (mc->d_val == NULL) {
    cs_matrix_copy_diagonal(matrix, mc->_d_val);
    mc->d_val = mc->_d_val;
  }
  diag = mc->d_val;

  return diag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function for initialization of CSR matrix coefficients using
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
_csr_assembler_values_init(void        *matrix_p,
                           cs_lnum_t    db_size,
                           cs_lnum_t    eb_size)
{
  CS_UNUSED(db_size);

  cs_matrix_t  *matrix = (cs_matrix_t *)matrix_p;

  cs_matrix_coeff_csr_t  *mc = matrix->coeffs;

  const cs_alloc_mode_t amode = matrix->alloc_mode;
  const cs_lnum_t n_rows = matrix->n_rows;
  cs_lnum_t e_size_2 = eb_size*eb_size;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  /* Initialize diagonal values */

  CS_FREE_HD(mc->_val);
  CS_MALLOC_HD(mc->_val, e_size_2*ms->row_index[ms->n_rows], cs_real_t, amode);
  mc->val = mc->_val;

# pragma omp parallel for  if(n_rows*db_size > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
    cs_lnum_t n_s_cols = (ms->row_index[ii+1] - ms->row_index[ii])*e_size_2;
    cs_lnum_t displ = ms->row_index[ii]*e_size_2;
    for (cs_lnum_t jj = 0; jj < n_s_cols; jj++)
      mc->_val[displ + jj] = 0;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function for addition to CSR matrix coefficients using
 *        local row ids and column indexes.
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
 * \param[in]       row_id    associated local row ids
 * \param[in]       col_idx   associated local column indexes
 * \param[in]       vals      pointer to values (size: n*stride)
 */
/*----------------------------------------------------------------------------*/

static void
_csr_assembler_values_add(void             *matrix_p,
                          cs_lnum_t         n,
                          cs_lnum_t         stride,
                          const cs_lnum_t   row_id[],
                          const cs_lnum_t   col_idx[],
                          const cs_real_t   vals[])
{
  cs_matrix_t  *matrix = (cs_matrix_t *)matrix_p;

  cs_matrix_coeff_csr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  bool sub_threads = (n*stride > _base_assembler_thr_min) ? true : false;
#if defined(_OPENMP)
  if (omp_in_parallel()) sub_threads = false;
#endif

  if (stride == 1) {

    /* Copy instead of test for OpenMP to avoid outlining for small sets */

    if (sub_threads == false) {
      for (cs_lnum_t ii = 0; ii < n; ii++) {
        if (row_id[ii] < 0)
          continue;
        else {
          cs_lnum_t r_id = row_id[ii];
          mc->_val[ms->row_index[r_id] + col_idx[ii]] += vals[ii];
        }
      }
    }

    else {
#     pragma omp parallel for
      for (cs_lnum_t ii = 0; ii < n; ii++) {
        if (row_id[ii] < 0)
          continue;
        else {
          cs_lnum_t r_id = row_id[ii];
#         pragma omp atomic
          mc->_val[ms->row_index[r_id] + col_idx[ii]] += vals[ii];
        }
      }
    }
  }

  else { /* if (stride > 1) */

    /* Copy instead of test for OpenMP to avoid outlining for small sets */

    if (sub_threads == false) {
      for (cs_lnum_t ii = 0; ii < n; ii++) {
        if (row_id[ii] < 0)
          continue;
        else {
          cs_lnum_t r_id = row_id[ii];
          cs_lnum_t displ = (ms->row_index[r_id] + col_idx[ii])*stride;
          for (cs_lnum_t jj = 0; jj < stride; jj++)
            mc->_val[displ + jj] += vals[ii*stride + jj];
        }
      }
    }

    else {
#     pragma omp parallel for
      for (cs_lnum_t ii = 0; ii < n; ii++) {
        if (row_id[ii] < 0)
          continue;
        else {
          cs_lnum_t r_id = row_id[ii];
          cs_lnum_t displ = (ms->row_index[r_id] + col_idx[ii])*stride;
          for (cs_lnum_t jj = 0; jj < stride; jj++)
#           pragma omp atomic
            mc->_val[displ + jj] += vals[ii*stride + jj];
        }
      }
    }
  }

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
_assembler_values_create_csr(cs_matrix_t  *matrix,
                             cs_lnum_t     diag_block_size,
                             cs_lnum_t     extra_diag_block_size)
{
  cs_matrix_assembler_values_t *mav
    = cs_matrix_assembler_values_create(matrix->assembler,
                                        false,
                                        diag_block_size,
                                        extra_diag_block_size,
                                        (void *)matrix,
                                        _csr_assembler_values_init,
                                        _csr_assembler_values_add,
                                        NULL,
                                        NULL,
                                        NULL);

  return mav;
}

/*----------------------------------------------------------------------------
 * Set MSR extradiagonal matrix coefficients for the case where direct
 * assignment is possible (i.e. when there are no multiple contributions
 * to a given coefficient).
 *
 * parameters:
 *   matrix      <-- pointer to matrix structure
 *   symmetric   <-- indicates if extradiagonal values are symmetric
 *   n_edges     <-- local number of graph edges
 *   edges       <-- edges (symmetric row <-> column) connectivity
 *   xa          <-- extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_e_coeffs_msr_direct(cs_matrix_t        *matrix,
                         bool                symmetric,
                         cs_lnum_t           n_edges,
                         const cs_lnum_2_t  *edges,
                         const cs_real_t    *restrict xa)
{
  cs_lnum_t  ii, jj, face_id;
  cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_matrix_struct_csr_t  *ms_e = &(ms->e);

  /* Copy extra-diagonal values */

  assert(edges != NULL || n_edges == 0);

  if (symmetric == false) {

    const cs_lnum_t *restrict edges_p
      = (const cs_lnum_t *restrict)(edges);

    for (face_id = 0; face_id < n_edges; face_id++) {
      cs_lnum_t kk, ll;
      ii = *edges_p++;
      jj = *edges_p++;
      if (ii < ms->n_rows) {
        for (kk = ms_e->row_index[ii]; ms_e->col_id[kk] != jj; kk++);
        mc->_e_val[kk] = xa[2*face_id];
      }
      if (jj < ms_e->n_rows) {
        for (ll = ms_e->row_index[jj]; ms_e->col_id[ll] != ii; ll++);
        mc->_e_val[ll] = xa[2*face_id + 1];
      }
    }

  }
  else { /* if symmetric == true */

    const cs_lnum_t *restrict edges_p
      = (const cs_lnum_t *restrict)(edges);

    for (face_id = 0; face_id < n_edges; face_id++) {
      cs_lnum_t kk, ll;
      ii = *edges_p++;
      jj = *edges_p++;
      if (ii < ms_e->n_rows) {
        for (kk = ms_e->row_index[ii]; ms_e->col_id[kk] != jj; kk++);
        mc->_e_val[kk] = xa[face_id];
      }
      if (jj < ms_e->n_rows) {
        for (ll = ms_e->row_index[jj]; ms_e->col_id[ll] != ii; ll++);
        mc->_e_val[ll] = xa[face_id];
      }

    }

  } /* end of condition on coefficients symmetry */

}

/*----------------------------------------------------------------------------
 * Set MSR extradiagonal block matrix coefficients for the case where direct
 * assignment is possible (i.e. when there are no multiple contributions
 * to a given coefficient).
 *
 * parameters:
 *   matrix      <-- pointer to matrix structure
 *   symmetric   <-- indicates if extradiagonal values are symmetric
 *   n_edges     <-- local number of graph edges
 *   edges       <-- edges (symmetric row <-> column) connectivity
 *   xa          <-- extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_e_coeffs_msr_direct_block(cs_matrix_t        *matrix,
                               bool                symmetric,
                               cs_lnum_t           n_edges,
                               const cs_lnum_2_t  *edges,
                               const cs_real_t    *restrict xa)
{
  cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_matrix_struct_csr_t  *ms_e = &(ms->e);

  const cs_lnum_t  b_size_2 = matrix->eb_size * matrix->eb_size;

  /* Copy extra-diagonal values */

  assert(edges != NULL || n_edges == 0);

  if (symmetric == false) {

    const cs_lnum_t *restrict edges_p
      = (const cs_lnum_t *restrict)(edges);

    for (cs_lnum_t face_id = 0; face_id < n_edges; face_id++) {
      cs_lnum_t ii = *edges_p++;
      cs_lnum_t jj = *edges_p++;
      cs_lnum_t kk, ll;
      if (ii < ms_e->n_rows) {
        for (kk = ms_e->row_index[ii]; ms_e->col_id[kk] != jj; kk++);
        for (cs_lnum_t pp = 0; pp < b_size_2; pp++)
          mc->_e_val[kk*b_size_2 + pp] = xa[2*face_id*b_size_2 + pp];
      }
      if (jj < ms_e->n_rows) {
        for (ll = ms_e->row_index[jj]; ms_e->col_id[ll] != ii; ll++);
        for (cs_lnum_t pp = 0; pp < b_size_2; pp++)
          mc->_e_val[ll*b_size_2 + pp] = xa[(2*face_id+1)*b_size_2 + pp];
      }
    }

  }
  else { /* if symmetric == true */

    const cs_lnum_t *restrict edges_p
      = (const cs_lnum_t *restrict)(edges);

    for (cs_lnum_t face_id = 0; face_id < n_edges; face_id++) {
      cs_lnum_t ii = *edges_p++;
      cs_lnum_t jj = *edges_p++;
      cs_lnum_t kk, ll;
      if (ii < ms_e->n_rows) {
        for (kk = ms_e->row_index[ii]; ms_e->col_id[kk] != jj; kk++);
        for (cs_lnum_t pp = 0; pp < b_size_2; pp++)
          mc->_e_val[kk*b_size_2 + pp] = xa[face_id*b_size_2 + pp];
      }
      if (jj < ms_e->n_rows) {
        for (ll = ms_e->row_index[jj]; ms_e->col_id[ll] != ii; ll++);
        for (cs_lnum_t pp = 0; pp < b_size_2; pp++)
          mc->_e_val[ll*b_size_2 + pp] = xa[face_id*b_size_2 + pp];
      }

    }

  } /* end of condition on coefficients symmetry */

}

/*----------------------------------------------------------------------------
 * Set MSR extradiagonal matrix coefficients for the case where there are
 * multiple contributions to a given coefficient.
 *
 * The matrix coefficients should have been initialized (i.e. set to 0)
 * some before using this function.
 *
 * parameters:
 *   matrix      <-- pointer to matrix structure
 *   symmetric   <-- indicates if extradiagonal values are symmetric
 *   n_edges     <-- local number of graph edges
 *   edges       <-- edges (symmetric row <-> column) connectivity
 *   xa          <-- extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_e_coeffs_msr_increment(cs_matrix_t        *matrix,
                            bool                symmetric,
                            cs_lnum_t           n_edges,
                            const cs_lnum_2_t  *edges,
                            const cs_real_t    *restrict xa)
{
  cs_lnum_t  ii, jj, face_id;
  cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_matrix_struct_csr_t  *ms_e = &(ms->e);

  /* Copy extra-diagonal values */

  assert(edges != NULL);

  if (symmetric == false) {

    const cs_lnum_t *restrict edges_p
      = (const cs_lnum_t *restrict)(edges);

    for (face_id = 0; face_id < n_edges; face_id++) {
      cs_lnum_t kk, ll;
      ii = *edges_p++;
      jj = *edges_p++;
      if (ii < ms_e->n_rows) {
        for (kk = ms_e->row_index[ii]; ms_e->col_id[kk] != jj; kk++);
        mc->_e_val[kk] += xa[2*face_id];
      }
      if (jj < ms_e->n_rows) {
        for (ll = ms_e->row_index[jj]; ms_e->col_id[ll] != ii; ll++);
        mc->_e_val[ll] += xa[2*face_id + 1];
      }
    }

  }
  else { /* if symmetric == true */

    const cs_lnum_t *restrict edges_p
      = (const cs_lnum_t *restrict)(edges);

    for (face_id = 0; face_id < n_edges; face_id++) {
      cs_lnum_t kk, ll;
      ii = *edges_p++;
      jj = *edges_p++;
      if (ii < ms_e->n_rows) {
        for (kk = ms_e->row_index[ii]; ms_e->col_id[kk] != jj; kk++);
        mc->_e_val[kk] += xa[face_id];
      }
      if (jj < ms_e->n_rows) {
        for (ll = ms_e->row_index[jj]; ms_e->col_id[ll] != ii; ll++);
        mc->_e_val[ll] += xa[face_id];
      }

    }

  } /* end of condition on coefficients symmetry */

}

/*----------------------------------------------------------------------------
 * Set MSR extradiagonal matrix coefficients for the case where there are
 * multiple contributions to a given coefficient.
 *
 * The matrix coefficients should have been initialized (i.e. set to 0)
 * some before using this function.
 *
 * parameters:
 *   matrix      <-- pointer to matrix structure
 *   symmetric   <-- indicates if extradiagonal values are symmetric
 *   n_edges     <-- local number of graph edges
 *   edges       <-- edges (symmetric row <-> column) connectivity
 *   xa          <-- extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_e_coeffs_msr_increment_block(cs_matrix_t        *matrix,
                                  bool                symmetric,
                                  cs_lnum_t           n_edges,
                                  const cs_lnum_2_t  *edges,
                                  const cs_real_t    *restrict xa)
{
  cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_matrix_struct_csr_t  *ms_e = &(ms->e);

  const cs_lnum_t  b_size_2 = matrix->eb_size * matrix->eb_size;

  /* Copy extra-diagonal values */

  assert(edges != NULL);

  if (symmetric == false) {

    const cs_lnum_t *restrict edges_p
      = (const cs_lnum_t *restrict)(edges);

    for (cs_lnum_t face_id = 0; face_id < n_edges; face_id++) {
      cs_lnum_t ii = *edges_p++;
      cs_lnum_t jj = *edges_p++;
      cs_lnum_t kk, ll;
      if (ii < ms_e->n_rows) {
        for (kk = ms_e->row_index[ii]; ms_e->col_id[kk] != jj; kk++);
        for (cs_lnum_t pp = 0; pp < b_size_2; pp++)
          mc->_e_val[kk*b_size_2 + pp] += xa[2*face_id*b_size_2 + pp];
      }
      if (jj < ms_e->n_rows) {
        for (ll = ms_e->row_index[jj]; ms_e->col_id[ll] != ii; ll++);
        for (cs_lnum_t pp = 0; pp < b_size_2; pp++)
          mc->_e_val[ll*b_size_2 + pp] += xa[(2*face_id+1)*b_size_2 + pp];
      }
    }

  }
  else { /* if symmetric == true */

    const cs_lnum_t *restrict edges_p
      = (const cs_lnum_t *restrict)(edges);

    for (cs_lnum_t face_id = 0; face_id < n_edges; face_id++) {
      cs_lnum_t ii = *edges_p++;
      cs_lnum_t jj = *edges_p++;
      cs_lnum_t kk, ll;
      if (ii < ms_e->n_rows) {
        for (kk = ms_e->row_index[ii]; ms_e->col_id[kk] != jj; kk++);
        for (cs_lnum_t pp = 0; pp < b_size_2; pp++)
          mc->_e_val[kk*b_size_2 + pp] += xa[face_id*b_size_2 + pp];
      }
      if (jj < ms_e->n_rows) {
        for (ll = ms_e->row_index[jj]; ms_e->col_id[ll] != ii; ll++);
        for (cs_lnum_t pp = 0; pp < b_size_2; pp++)
          mc->_e_val[ll*b_size_2 + pp] += xa[face_id*b_size_2 + pp];
      }

    }

  } /* end of condition on coefficients symmetry */

}

/*----------------------------------------------------------------------------
 * Map or copy MSR matrix diagonal coefficients.
 *
 * parameters:
 *   matrix           <-> pointer to matrix structure
 *   copy             <-- indicates if coefficients should be copied
 *   da               <-- diagonal values (NULL if all zero)
 *----------------------------------------------------------------------------*/

static void
_map_or_copy_d_coeffs_msr(cs_matrix_t      *matrix,
                          bool              copy,
                          const cs_real_t  *restrict da)
{
  cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_lnum_t n_rows = matrix->n_rows;

  /* Map or copy diagonal values */

  mc->db_size = matrix->db_size;
  CS_FREE(mc->_d_val);

  if (da != NULL) {

#if defined(HAVE_ACCEL)
    if (   cs_check_device_ptr(da) == CS_ALLOC_HOST
        && matrix->alloc_mode > CS_ALLOC_HOST)
      copy = true;
#endif

    if (copy) {
      const cs_lnum_t b_size = mc->db_size;
      const cs_lnum_t b_size_2 = b_size * b_size;
      CS_MALLOC_HD(mc->_d_val, b_size_2*n_rows, cs_real_t, matrix->alloc_mode);
#     pragma omp parallel for  if(n_rows*b_size > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
        for (cs_lnum_t jj = 0; jj < b_size_2; jj++)
          mc->_d_val[ii*b_size_2 + jj] = da[ii*b_size_2 + jj];
      }
      mc->d_val = mc->_d_val;
    }
    else
      mc->d_val = da;

  }
  else
    mc->d_val = NULL;
}

/*----------------------------------------------------------------------------
 * Map or copy MSR matrix extra diagonal coefficients.
 *
 * This assumes the xa values are already provided in MSR form.
 *
 * Setting xa = NULL and copy = true, this function also ensures allocation
 * of and zeroes extradiagonal coefficients.
 *
 * parameters:
 *   matrix           <-> pointer to matrix structure
 *   copy             <-- indicates if coefficients should be copied
 *   xa               <-- extradiagonal values (NULL if all zero)
 *----------------------------------------------------------------------------*/

static void
_map_or_copy_e_coeffs_msr(cs_matrix_t      *matrix,
                          bool              copy,
                          const cs_real_t  *restrict xa)
{
  cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_matrix_struct_csr_t  *ms_e = &(ms->e);

  const cs_lnum_t n_rows = matrix->n_rows;

  mc->eb_size = matrix->eb_size;
  CS_FREE(mc->_e_val);

  if (xa == NULL || copy) {

    const cs_lnum_t eb_size = mc->eb_size;
    const cs_lnum_t eb_size_2 = eb_size * eb_size;

    CS_MALLOC_HD(mc->_e_val,
                 eb_size_2*ms_e->row_index[ms_e->n_rows],
                 cs_real_t,
                 matrix->alloc_mode);
    mc->e_val = mc->_e_val;

    /* zero if required */
    if (xa == NULL)
      _zero_coeffs_csr(ms_e, mc->eb_size, mc->_e_val);

  }

  /* Map or copy extradiagonal values (we could use memcpy, but prefer
     to have a similar threading behavior to SpMV for NUMA performance) */

  if (xa != NULL) {

    if (copy) {
      if (mc->eb_size == 1) {
#       pragma omp parallel for  if(n_rows > CS_THR_MIN)
        for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
          const cs_lnum_t  n_cols = ms_e->row_index[ii+1] - ms_e->row_index[ii];
          const cs_real_t  *s_row = xa + ms_e->row_index[ii];
          cs_real_t  *m_row = mc->_e_val + ms_e->row_index[ii];
          for (cs_lnum_t jj = 0; jj < n_cols; jj++)
            m_row[jj] = s_row[jj];
        }
      }
      else {
        const cs_lnum_t b_size = mc->eb_size;
        const cs_lnum_t b_size_2 = b_size * b_size;
#       pragma omp parallel for  if(n_rows*b_size > CS_THR_MIN)
        for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
          const cs_lnum_t  n_cols = ms_e->row_index[ii+1] - ms_e->row_index[ii];
          const cs_real_t  *s_row = xa + ms_e->row_index[ii]*b_size_2;
          cs_real_t  *m_row = mc->_e_val + ms_e->row_index[ii]*b_size_2;
          for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
            for (cs_lnum_t kk = 0; kk < b_size_2; kk++)
              m_row[jj*b_size_2 + kk] = s_row[jj*b_size_2 + kk];
          }
        }
      }
    }

    else
      mc->e_val = xa;

  }
}

/*----------------------------------------------------------------------------
 * Set MSR matrix coefficients.
 *
 * parameters:
 *   matrix      <-> pointer to matrix structure
 *   symmetric   <-- indicates if extradiagonal values are symmetric
 *   copy        <-- indicates if coefficients should be copied
 *   n_edges     <-- local number of graph edges
 *   edges       <-- edges (symmetric row <-> column) connectivity
 *   da          <-- diagonal values (NULL if all zero)
 *   xa          <-- extradiagonal values (NULL if all zero)
 *----------------------------------------------------------------------------*/

static void
_set_coeffs_msr(cs_matrix_t         *matrix,
                bool                 symmetric,
                bool                 copy,
                cs_lnum_t            n_edges,
                const cs_lnum_2_t  *restrict edges,
                const cs_real_t    *restrict da,
                const cs_real_t    *restrict xa)
{
  cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_matrix_struct_csr_t  *ms_e = &(ms->e);

  /* Map or copy diagonal values */

  _map_or_copy_d_coeffs_msr(matrix, copy, da);

  /* Extradiagonal values */

  mc->eb_size = matrix->eb_size;

  const cs_lnum_t eb_size = mc->eb_size;
  const cs_lnum_t eb_size_2 = eb_size * eb_size;

  CS_FREE(mc->_e_val);
  CS_MALLOC_HD(mc->_e_val,
               eb_size_2*ms_e->row_index[ms_e->n_rows],
               cs_real_t,
               matrix->alloc_mode);
  mc->e_val = mc->_e_val;

  /* Copy extra-diagonal values if assembly is direct */

  if (ms_e->direct_assembly) {
    if (xa == NULL)
      _zero_coeffs_csr(ms_e, mc->eb_size, mc->_e_val);
    if (eb_size == 1)
      _set_e_coeffs_msr_direct(matrix, symmetric, n_edges, edges, xa);
    else
      _set_e_coeffs_msr_direct_block(matrix, symmetric, n_edges, edges, xa);
  }

  /* Initialize coefficients to zero if assembly is incremental */

  else {
    _zero_coeffs_csr(ms_e, mc->eb_size, mc->_e_val);
    if (eb_size == 1)
      _set_e_coeffs_msr_increment(matrix, symmetric, n_edges, edges, xa);
    else
      _set_e_coeffs_msr_increment_block(matrix, symmetric, n_edges, edges, xa);
  }
}

/*----------------------------------------------------------------------------
 * Set MSR matrix coefficients provided in the same form.
 *
 * If da and xa are equal to NULL, then initialize val with zeros.
 *
 * parameters:
 *   matrix           <-> pointer to matrix structure
 *   copy             <-- indicates if coefficients should be copied
 *                        when not transferred
 *   row_index        <-- MSR row index (0 to n-1)
 *   col_id           <-- MSR column id (0 to n-1)
 *   d_vals           <-- diagonal values (NULL if all zero)
 *   d_vals_transfer  <-- diagonal values whose ownership is transferred
 *                        (NULL or d_vals in, NULL out)
 *   x_vals           <-- extradiagonal values (NULL if all zero)
 *   x_vals_transfer  <-- extradiagonal values whose ownership is transferred
 *                        (NULL or x_vals in, NULL out)
 *----------------------------------------------------------------------------*/

static void
_set_coeffs_msr_from_msr(cs_matrix_t       *matrix,
                         bool               copy,
                         const cs_lnum_t    row_index[],
                         const cs_lnum_t    col_id[],
                         const cs_real_t   *d_vals,
                         cs_real_t        **d_vals_transfer,
                         const cs_real_t   *x_vals,
                         cs_real_t        **x_vals_transfer)
{
  CS_UNUSED(row_index);
  CS_UNUSED(col_id);

  cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  bool d_transferred = false, x_transferred = false;

  /* TODO: we should use metadata or check that the row_index and
     column id values are consistent, which should be true as long
     as columns are ordered in an identical manner */

  if (d_vals_transfer != NULL) {
    if (*d_vals_transfer != NULL) {
      mc->db_size = matrix->db_size;
      if (mc->_d_val != *d_vals_transfer) {
        CS_FREE(mc->_d_val);
        mc->_d_val = *d_vals_transfer;
      }
      mc->d_val = mc->_d_val;
      *d_vals_transfer = NULL;
      d_transferred = true;
    }
  }

  if (x_vals_transfer != NULL) {
    if (*x_vals_transfer != NULL) {
      mc->eb_size = matrix->eb_size;
      if (mc->_e_val != *x_vals_transfer) {
        BFT_FREE(mc->_e_val);
        mc->_e_val = *x_vals_transfer;
      }
      mc->e_val = mc->_e_val;
      *x_vals_transfer = NULL;
      x_transferred = true;
    }
  }

  if (d_transferred == false)
    _map_or_copy_d_coeffs_msr(matrix, copy, d_vals);

  if (x_transferred == false)
    _map_or_copy_e_coeffs_msr(matrix, copy, x_vals);

  /* Now free transferred arrays */

  if (d_vals_transfer != NULL)
    CS_FREE(*d_vals_transfer);
  if (x_vals_transfer != NULL)
    CS_FREE(*x_vals_transfer);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function for initialization of distributed matrix coefficients using
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
_msr_assembler_values_init(void              *matrix_p,
                           const cs_lnum_t    db_size,
                           const cs_lnum_t    eb_size)
{
  cs_matrix_t  *matrix = (cs_matrix_t *)matrix_p;

  cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_alloc_mode_t amode = matrix->alloc_mode;
  const cs_lnum_t n_rows = matrix->n_rows;

  mc->db_size = db_size;
  mc->eb_size = eb_size;

  const cs_lnum_t d_size_2 = mc->db_size * mc->db_size;
  const cs_lnum_t e_size_2 = mc->eb_size * mc->eb_size;

  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_matrix_struct_csr_t  *ms_e = &(ms->e);

  /* Allocate values and initialize to zero. */

  CS_FREE_HD(mc->d_idx);
  CS_FREE(mc->_d_val);
  CS_FREE(mc->_e_val);

  CS_FREE_HD(mc->_h_val);
  mc->h_val = NULL;

  CS_MALLOC_HD(mc->_d_val, d_size_2*n_rows, cs_real_t, amode);
  mc->d_val = mc->_d_val;

  cs_lnum_t nnz_e = e_size_2*ms_e->row_index[ms_e->n_rows];

  CS_MALLOC_HD(mc->_e_val, nnz_e, cs_real_t, amode);
  mc->e_val = mc->_e_val;

# pragma omp parallel for  if(n_rows*db_size > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
    for (cs_lnum_t jj = 0; jj < d_size_2; jj++)
      mc->_d_val[ii*d_size_2 + jj] = 0;
    cs_lnum_t n_e_cols = (ms_e->row_index[ii+1] - ms_e->row_index[ii])*e_size_2;
    cs_lnum_t displ_e = ms_e->row_index[ii]*e_size_2;
    for (cs_lnum_t jj = 0; jj < n_e_cols; jj++)
      mc->_e_val[displ_e + jj] = 0;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function for addition to MSR matrix coefficients using
 *        local row ids and column indexes.
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
 * \param[in]       row_id    associated local row ids
 * \param[in]       col_idx   associated local column indexes
 * \param[in]       vals      pointer to values (size: n*stride)
 */
/*----------------------------------------------------------------------------*/

static void
_msr_assembler_values_add(void             *matrix_p,
                          cs_lnum_t         n,
                          cs_lnum_t         stride,
                          const cs_lnum_t   row_id[],
                          const cs_lnum_t   col_idx[],
                          const cs_real_t   vals[])
{
  cs_matrix_t  *matrix = (cs_matrix_t *)matrix_p;

  cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_matrix_struct_csr_t  *ms_e = &(ms->e);

  bool sub_threads = (n*stride > _base_assembler_thr_min) ? true : false;
#if defined(_OPENMP)
  if (omp_in_parallel()) sub_threads = false;
#endif

  if (stride == 1) {

    /* Copy instead of test for OpenMP to avoid outlining for small sets */

    if (sub_threads == false) {
      for (cs_lnum_t ii = 0; ii < n; ii++) {
        cs_lnum_t r_id = row_id[ii];
        if (r_id < 0)
          continue;
        if (col_idx[ii] < 0) {
          mc->_d_val[r_id] += vals[ii];
        }
        else {
          mc->_e_val[ms_e->row_index[r_id] + col_idx[ii]] += vals[ii];
        }
      }
    }

    else {
#     pragma omp parallel for
      for (cs_lnum_t ii = 0; ii < n; ii++) {
        cs_lnum_t r_id = row_id[ii];
        if (r_id < 0)
          continue;
        if (col_idx[ii] < 0) {
#         pragma omp atomic
          mc->_d_val[r_id] += vals[ii];
        }
        else {
#         pragma omp atomic
          mc->_e_val[ms_e->row_index[r_id] + col_idx[ii]] += vals[ii];
        }
      }
    }
  }

  else { /* if (stride > 1) */

    /* Copy instead of test for OpenMP to avoid outlining for small sets */

    if (sub_threads == false) {
      for (cs_lnum_t ii = 0; ii < n; ii++) {
        cs_lnum_t r_id = row_id[ii];
        if (r_id < 0)
          continue;
        if (col_idx[ii] < 0) {
          for (cs_lnum_t jj = 0; jj < stride; jj++)
            mc->_d_val[r_id*stride + jj] += vals[ii*stride + jj];
        }
        else {
          cs_lnum_t displ = (ms_e->row_index[r_id] + col_idx[ii])*stride;
          for (cs_lnum_t jj = 0; jj < stride; jj++)
            mc->_e_val[displ + jj] += vals[ii*stride + jj];
        }
      }
    }

    else {
#     pragma omp parallel for
      for (cs_lnum_t ii = 0; ii < n; ii++) {
        cs_lnum_t r_id = row_id[ii];
        if (r_id < 0)
          continue;
        if (col_idx[ii] < 0) {
          for (cs_lnum_t jj = 0; jj < stride; jj++)
#           pragma omp atomic
            mc->_d_val[r_id*stride + jj] += vals[ii*stride + jj];
        }
        else {
          cs_lnum_t displ = (ms_e->row_index[r_id] + col_idx[ii])*stride;
          for (cs_lnum_t jj = 0; jj < stride; jj++)
#           pragma omp atomic
            mc->_e_val[displ + jj] += vals[ii*stride + jj];
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize an MSR matrix assembler values structure.
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
_assembler_values_create_msr(cs_matrix_t      *matrix,
                             const cs_lnum_t   diag_block_size,
                             const cs_lnum_t   extra_diag_block_size)
{
  cs_matrix_assembler_values_t *mav
    = cs_matrix_assembler_values_create(matrix->assembler,
                                        true,
                                        diag_block_size,
                                        extra_diag_block_size,
                                        (void *)matrix,
                                        _msr_assembler_values_init,
                                        _msr_assembler_values_add,
                                        NULL,
                                        NULL,
                                        NULL);

  return mav;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function for initialization of distributed matrix coefficients using
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
_dist_assembler_values_init(void        *matrix_p,
                            cs_lnum_t    db_size,
                            cs_lnum_t    eb_size)
{
  cs_matrix_t  *matrix = (cs_matrix_t *)matrix_p;

  cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_alloc_mode_t amode = matrix->alloc_mode;
  const cs_lnum_t n_rows = matrix->n_rows;

  mc->db_size = db_size;
  mc->eb_size = eb_size;

  const cs_lnum_t d_size_2 = mc->db_size * mc->db_size;
  const cs_lnum_t e_size_2 = mc->eb_size * mc->eb_size;

  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_matrix_struct_csr_t  *ms_e = &(ms->e);
  const cs_matrix_struct_csr_t  *ms_h = &(ms->h);

  /* Allocate values and initialize to zero. */

  CS_FREE_HD(mc->d_idx);
  CS_FREE(mc->_d_val);
  CS_FREE(mc->_e_val);
  CS_FREE_HD(mc->_h_val);

  CS_MALLOC_HD(mc->_d_val, d_size_2*n_rows, cs_real_t, amode);
  mc->d_val = mc->_d_val;

  cs_lnum_t nnz_e = e_size_2*ms_e->row_index[ms_e->n_rows];
  cs_lnum_t nnz_h = e_size_2*ms_h->row_index[ms_h->n_rows];

  CS_MALLOC_HD(mc->_e_val, nnz_e, cs_real_t, amode);
  CS_MALLOC_HD(mc->_h_val, nnz_h, cs_real_t, amode);
  mc->e_val = mc->_e_val;
  mc->h_val = mc->_h_val;

# pragma omp parallel for  if(n_rows*db_size > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
    for (cs_lnum_t jj = 0; jj < d_size_2; jj++)
      mc->_d_val[ii*d_size_2 + jj] = 0;
    cs_lnum_t n_e_cols = (ms_e->row_index[ii+1] - ms_e->row_index[ii])*e_size_2;
    cs_lnum_t displ_e = ms_e->row_index[ii]*e_size_2;
    for (cs_lnum_t jj = 0; jj < n_e_cols; jj++)
      mc->_e_val[displ_e + jj] = 0;
    cs_lnum_t n_h_cols = (ms_h->row_index[ii+1] - ms_h->row_index[ii])*e_size_2;
    cs_lnum_t displ_h = ms_h->row_index[ii]*e_size_2;
    for (cs_lnum_t jj = 0; jj < n_h_cols; jj++)
      mc->_h_val[displ_h + jj] = 0;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function for addition to distributed matrix coefficients using
 *        local row ids and column indexes.
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
 * \param[in]       row_id    associated local row ids
 * \param[in]       col_idx   associated local column indexes
 * \param[in]       vals      pointer to values (size: n*stride)
 */
/*----------------------------------------------------------------------------*/

static void
_dist_assembler_values_add(void             *matrix_p,
                           cs_lnum_t         n,
                           cs_lnum_t         stride,
                           const cs_lnum_t   row_id[],
                           const cs_lnum_t   col_idx[],
                           const cs_real_t   vals[])
{
  cs_matrix_t  *matrix = (cs_matrix_t *)matrix_p;

  cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_matrix_struct_csr_t  *ms_e = &(ms->e);
  const cs_matrix_struct_csr_t  *ms_h = &(ms->h);

  bool sub_threads = (n*stride > _base_assembler_thr_min) ? true : false;
#if defined(_OPENMP)
  if (omp_in_parallel()) sub_threads = false;
#endif

  if (stride == 1) {

    /* Copy instead of test for OpenMP to avoid outlining for small sets */

    if (sub_threads == false) {
      for (cs_lnum_t ii = 0; ii < n; ii++) {
        cs_lnum_t r_id = row_id[ii];
        if (r_id < 0)
          continue;
        if (col_idx[ii] < 0) {
          mc->_d_val[r_id] += vals[ii];
        }
        else {
          cs_lnum_t idx_s =   ms_e->row_index[r_id] + col_idx[ii]
                            - ms_e->row_index[r_id+1];
          if (idx_s < 0) {
            mc->_e_val[ms_e->row_index[r_id] + col_idx[ii]] += vals[ii];
          }
          else {
            mc->_h_val[ms_h->row_index[r_id] + idx_s] += vals[ii];
          }
        }
      }
    }

    else {
#     pragma omp parallel for
      for (cs_lnum_t ii = 0; ii < n; ii++) {
        cs_lnum_t r_id = row_id[ii];
        if (r_id < 0)
          continue;
        if (col_idx[ii] < 0) {
#         pragma omp atomic
          mc->_d_val[r_id] += vals[ii];
        }
        else {
          cs_lnum_t idx_s =   ms_e->row_index[r_id] + col_idx[ii]
                            - ms_e->row_index[r_id+1];
          if (idx_s < 0) {
#           pragma omp atomic
            mc->_e_val[ms_e->row_index[r_id] + col_idx[ii]] += vals[ii];
          }
          else {
#           pragma omp atomic
            mc->_h_val[ms_h->row_index[r_id] + idx_s] += vals[ii];
          }
        }
      }
    }
  }

  else { /* if (stride > 1) */

    /* Copy instead of test for OpenMP to avoid outlining for small sets */

    if (sub_threads == false) {
      for (cs_lnum_t ii = 0; ii < n; ii++) {
        cs_lnum_t r_id = row_id[ii];
        if (r_id < 0)
          continue;
        if (col_idx[ii] < 0) {
          for (cs_lnum_t jj = 0; jj < stride; jj++)
            mc->_d_val[r_id*stride + jj] += vals[ii*stride + jj];
        }
        else {
          cs_lnum_t idx_s =   ms_e->row_index[r_id] + col_idx[ii]
                            - ms_e->row_index[r_id+1];
          if (idx_s < 0) {
            cs_lnum_t displ = (ms_e->row_index[r_id] + col_idx[ii])*stride;
            for (cs_lnum_t jj = 0; jj < stride; jj++)
              mc->_e_val[displ + jj] += vals[ii*stride + jj];
          }
          else {
            cs_lnum_t displ = (ms_h->row_index[r_id] + idx_s)*stride;
            for (cs_lnum_t jj = 0; jj < stride; jj++)
              mc->_h_val[displ + jj] += vals[ii*stride + jj];
          }
        }
      }
    }

    else {
#     pragma omp parallel for
      for (cs_lnum_t ii = 0; ii < n; ii++) {
        cs_lnum_t r_id = row_id[ii];
        if (r_id < 0)
          continue;
        if (col_idx[ii] < 0) {
          for (cs_lnum_t jj = 0; jj < stride; jj++)
#           pragma omp atomic
            mc->_d_val[r_id*stride + jj] += vals[ii*stride + jj];
        }
        else {
          cs_lnum_t idx_s =   ms_e->row_index[r_id] + col_idx[ii]
                            - ms_e->row_index[r_id+1];
          if (idx_s < 0) {
            cs_lnum_t displ = (ms_e->row_index[r_id] + col_idx[ii])*stride;
            for (cs_lnum_t jj = 0; jj < stride; jj++)
#             pragma omp atomic
              mc->_e_val[displ + jj] += vals[ii*stride + jj];
          }
          else {
            cs_lnum_t displ = (ms_h->row_index[r_id] + idx_s)*stride;
            for (cs_lnum_t jj = 0; jj < stride; jj++)
#             pragma omp atomic
              mc->_h_val[displ + jj] += vals[ii*stride + jj];
          }
        }
      }
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize a distributed matrix assembler values structure.
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
_assembler_values_create_dist(cs_matrix_t  *matrix,
                              cs_lnum_t     diag_block_size,
                              cs_lnum_t     extra_diag_block_size)
{
  cs_matrix_assembler_values_t *mav
    = cs_matrix_assembler_values_create(matrix->assembler,
                                        true,
                                        diag_block_size,
                                        extra_diag_block_size,
                                        (void *)matrix,
                                        _dist_assembler_values_init,
                                        _dist_assembler_values_add,
                                        NULL,
                                        NULL,
                                        NULL);

  return mav;
}

/*----------------------------------------------------------------------------
 * Add variant
 *
 * parameters:
 *   type                 <-- matrix type
 *   mft                  <-- fill type tuned for
 *   func_names           <-- function names for SpMV operation types
 *   n_variants           <-> number of variants
 *   n_variants_max       <-> current maximum number of variants
 *   m_variant            <-> array of matrix variants
 *----------------------------------------------------------------------------*/

static void
_variant_add(const char                *name,
             cs_matrix_type_t           type,
             cs_matrix_fill_type_t      mft,
             const cs_numbering_t      *numbering,
             const char                *func_name,
             int                       *n_variants,
             int                       *n_variants_max,
             cs_matrix_variant_t      **m_variant)
{
  cs_matrix_variant_t  *v;
  int i = *n_variants;

  if (func_name == NULL)
    return;

  if (*n_variants_max == *n_variants) {
    if (*n_variants_max == 0)
      *n_variants_max = 8;
    else
      *n_variants_max *= 2;
    BFT_REALLOC(*m_variant, *n_variants_max, cs_matrix_variant_t);
  }

  v = (*m_variant) + i;

  for (int j = 0; j < CS_MATRIX_SPMV_N_TYPES; j++) {
    v->vector_multiply[j] = NULL;
    strncpy(v->name[j], name, 31);
    v->name[j][31] = '\0';
  }

  v->type = type;
  v->fill_type = mft;

  int retval = cs_matrix_spmv_set_func(type,
                                       mft,
                                       CS_MATRIX_SPMV_N_TYPES,
                                       numbering,
                                       func_name,
                                       v->vector_multiply,
                                       v->vector_multiply_xy_hd);

  if (retval == 0)
    *n_variants += 1;
}

/*----------------------------------------------------------------------------
 * Create a distributed matrix structure from a native matrix structure.
 *
 * Note that the structure created maps global cell numbers to the given
 * existing face -> cell connectivity array, so it must be destroyed before
 * this array (usually the code's global cell numbering) is freed.
 *
 * parameters:
 *   n     <-- number of elements
 *   index <-> count (shifted by 1) in, index out (size: n+1)
 *----------------------------------------------------------------------------*/

static void
_count_to_index(cs_lnum_t  n,
                cs_lnum_t  index[])
{
  for (cs_lnum_t i = 0; i < n; i++)
    index[i+1] += index[i];
}

/*----------------------------------------------------------------------------
 * Create a distributed matrix structure from a native matrix structure.
 *
 * Note that the structure created maps global cell numbers to the given
 * existing face -> cell connectivity array, so it must be destroyed before
 * this array (usually the code's global cell numbering) is freed.
 *
 * parameters:
 *   alloc_mode  <-- allocation mode
 *   n_rows      <-- number of local rows
 *   n_cols_ext  <-- number of local + ghost columns
 *   n_edges     <-- local number of graph edges
 *   edges       <-- edges (symmetric row <-> column) connectivity
 *
 * returns:
 *   pointer to allocated CSR matrix structure.
 *----------------------------------------------------------------------------*/

static cs_matrix_struct_dist_t *
_create_struct_dist(cs_alloc_mode_t     alloc_mode,
                    cs_lnum_t           n_rows,
                    cs_lnum_t           n_cols_ext,
                    cs_lnum_t           n_edges,
                    const cs_lnum_2_t  *edges)
{
  cs_matrix_struct_dist_t  *ms;

  /* Allocate and map */

  BFT_MALLOC(ms, 1, cs_matrix_struct_dist_t);

  ms->n_rows = n_rows;
  ms->n_cols_ext = n_cols_ext;

  if (n_edges == 0 || edges == NULL) {
    n_rows = 0;
    n_cols_ext = 0;
  }

  _init_struct_csr(&(ms->e), n_rows, n_rows);
  _init_struct_csr(&(ms->h), n_rows, n_cols_ext);

  ms->h_row_id = NULL;

  if (n_edges == 0 || edges == NULL)
    return ms;

  /* Count number of nonzero elements per row */

  const cs_lnum_t *restrict edges_p = (const cs_lnum_t *restrict)edges;

  cs_lnum_t index_size = n_rows+1;

  CS_MALLOC_HD(ms->h._row_index, n_rows+1, cs_lnum_t, alloc_mode);
  for (cs_lnum_t ii = 0; ii < index_size; ii++)
    ms->h._row_index[ii] = 0;

  CS_MALLOC_HD(ms->e._row_index, n_rows+1, cs_lnum_t, alloc_mode);

  for (cs_lnum_t ii = 0; ii < index_size; ii++)
    ms->e._row_index[ii] = 0;

  for (cs_lnum_t edge_id = 0; edge_id < n_edges; edge_id++) {
    cs_lnum_t ii = *edges_p++;
    cs_lnum_t jj = *edges_p++;
    if (ii < ms->n_rows) {
      if (jj < n_rows)
        ms->e._row_index[ii+1] += 1;
      else
        ms->h._row_index[ii+1] += 1;
    }
    if (jj < ms->n_rows) {
      if (ii < n_rows)
        ms->e._row_index[jj+1] += 1;
      else
        ms->h._row_index[jj+1] += 1;
    }
  }
  _count_to_index(n_rows, ms->e._row_index);
  _count_to_index(n_rows, ms->h._row_index);

  /* Now allocate column ids */

  CS_MALLOC_HD(ms->e._col_id, ms->e._row_index[n_rows], cs_lnum_t, alloc_mode);

  if (ms->h._row_index[n_rows] > 0)
    CS_MALLOC_HD(ms->h._col_id, ms->h._row_index[n_rows], cs_lnum_t, alloc_mode);
  else {
    CS_FREE_HD(ms->h._row_index);
    ms->h.n_rows = 0;
    ms->h.n_cols_ext = 0;
  }

  /* Prepare local count for shifts */

  cs_lnum_t count_size = 2*n_rows;
  cs_lnum_t *ccount;
  BFT_MALLOC(ccount, count_size, cs_lnum_t);

  for (cs_lnum_t ii = 0; ii < count_size; ii++)
    ccount[ii] = 0;

  /* Build structure */

  edges_p = (const cs_lnum_t *restrict)edges;

  for (cs_lnum_t edge_id = 0; edge_id < n_edges; edge_id++) {
    cs_lnum_t ii = *edges_p++;
    cs_lnum_t jj = *edges_p++;
    if (ii < ms->n_rows) {
      if (jj < n_rows) {
        ms->e._col_id[ms->e._row_index[ii] + ccount[ii*2]] = jj;
        ccount[ii*2] += 1;
      }
      else {
        ms->h._col_id[ms->h._row_index[ii] + ccount[ii*2+1]] = jj;
        ccount[ii*2+1] += 1;
      }
    }
    if (jj < ms->n_rows) {
      if (ii < n_rows) {
        ms->e._col_id[ms->e._row_index[jj] + ccount[jj*2]] = ii;
        ccount[jj*2] += 1;
      }
      else {
        ms->h._col_id[ms->h._row_index[jj] + ccount[jj*2+1]] = ii;
        ccount[jj*2+1] += 1;
      }
    }
  }

  BFT_FREE(ccount);

  /* Compact structures if needed */

  _compact_struct_csr(&(ms->e), alloc_mode);
  _compact_struct_csr(&(ms->h), alloc_mode);

  return ms;
}

/*----------------------------------------------------------------------------
 * Create an MSR matrix structure from a native matrix structure.
 *
 * Note that the structure created maps global cell numbers to the given
 * existing face -> cell connectivity array, so it must be destroyed before
 * this array (usually the code's global cell numbering) is freed.
 *
 * parameters:
 *   alloc_mode  <-- allocation mode
 *   n_rows      <-- number of local rows
 *   n_cols_ext  <-- number of local + ghost columns
 *   n_edges     <-- local number of graph edges
 *   edges       <-- edges (symmetric row <-> column) connectivity
 *
 * returns:
 *   pointer to a created distrubuted (MSR configured) matrix structure
 *----------------------------------------------------------------------------*/

static cs_matrix_struct_dist_t *
_create_struct_msr(cs_alloc_mode_t     alloc_mode,
                   cs_lnum_t           n_rows,
                   cs_lnum_t           n_cols_ext,
                   cs_lnum_t           n_edges,
                   const cs_lnum_2_t  *edges)
{
  cs_matrix_struct_dist_t  *ms;

  /* Allocate and map */

  BFT_MALLOC(ms, 1, cs_matrix_struct_dist_t);

  ms->n_rows = n_rows;
  ms->n_cols_ext = n_cols_ext;

  _init_struct_csr(&(ms->e), n_rows, n_cols_ext);
  _init_struct_csr(&(ms->h), 0, 0);

  cs_lnum_t index_size = n_rows+1;

  CS_MALLOC_HD(ms->e._row_index, index_size, cs_lnum_t, alloc_mode);

  for (cs_lnum_t ii = 0; ii < index_size; ii++)
    ms->e._row_index[ii] = 0;

  ms->e.row_index = ms->e._row_index;

  ms->h_row_id = NULL;

  if (n_edges == 0 || edges == NULL)
    return ms;

  /* Count number of nonzero elements per row */

  const cs_lnum_t *restrict edges_p = (const cs_lnum_t *restrict)edges;

  for (cs_lnum_t edge_id = 0; edge_id < n_edges; edge_id++) {
    cs_lnum_t ii = *edges_p++;
    cs_lnum_t jj = *edges_p++;
    if (ii < ms->n_rows)
      ms->e._row_index[ii+1] += 1;
    if (jj < ms->n_rows)
      ms->e._row_index[jj+1] += 1;
  }
  _count_to_index(n_rows, ms->e._row_index);

  /* Now allocate column ids */

  CS_MALLOC_HD(ms->e._col_id, ms->e._row_index[n_rows], cs_lnum_t, alloc_mode);

  /* Prepare local count for shifts */

  cs_lnum_t count_size = n_rows;
  cs_lnum_t *ccount;
  BFT_MALLOC(ccount, count_size, cs_lnum_t);

  for (cs_lnum_t ii = 0; ii < count_size; ii++)
    ccount[ii] = 0;

  /* Build structure */

  edges_p = (const cs_lnum_t *restrict)edges;

  for (cs_lnum_t edge_id = 0; edge_id < n_edges; edge_id++) {
    cs_lnum_t ii = *edges_p++;
    cs_lnum_t jj = *edges_p++;
    if (ii < ms->n_rows) {
      ms->e._col_id[ms->e._row_index[ii] + ccount[ii]] = jj;
      ccount[ii] += 1;
    }
    if (jj < ms->n_rows) {
      ms->e._col_id[ms->e._row_index[jj] + ccount[jj]] = ii;
      ccount[jj] += 1;
    }
  }

  BFT_FREE(ccount);

  /* Compact structures if needed */

  _compact_struct_csr(&(ms->e), alloc_mode);

  return ms;
}

/*----------------------------------------------------------------------------
 * Create an MSR matrix structure from an index and an array related
 * to column id
 *
 * parameters:
 *   separate_halo   <-- indicates if halo should be separated
 *   direct_assembly <-- true if each value corresponds to a unique face
 *   n_rows          <- local number of rows
 *   n_cols_ext      <-- local number of columns + ghosts
 *   row_index       <-- index on rows
 *   col_id          <-> array of column ids related to the row index
 *
 * returns:
 *    a pointer to a created distrubuted (MSR configured) matrix structure
 *----------------------------------------------------------------------------*/

static cs_matrix_struct_dist_t *
_create_struct_msr_from_shared(bool              direct_assembly,
                               cs_lnum_t         n_rows,
                               cs_lnum_t         n_cols_ext,
                               const cs_lnum_t  *row_index,
                               const cs_lnum_t  *col_id)
{
  cs_matrix_struct_dist_t  *ms = NULL;

  /* Allocate and map */

  BFT_MALLOC(ms, 1, cs_matrix_struct_dist_t);

  ms->n_rows = n_rows;
  ms->n_cols_ext = n_cols_ext;

  _init_struct_csr(&(ms->e), n_rows, n_cols_ext);
  _init_struct_csr(&(ms->h), 0, 0);

  ms->e.row_index = row_index;
  ms->e.col_id = col_id;

  ms->e._row_index = NULL;
  ms->e._col_id = NULL;

  ms->e.direct_assembly = direct_assembly;

  ms->h_row_id = NULL;

  return ms;
}

/*----------------------------------------------------------------------------
 * Create a CSR matrix structure from an index and an array related
 * to column id
 *
 * parameters:
 *   transfer   <-- transfer property of row_index and col_id
 *                  if true, map them otherwise
 *   ordered    <-- indicates if row entries are already ordered
 *   n_rows     <-- local number of rows
 *   n_cols_ext <-- local number of columns + ghosts
 *   row_index  <-- pointer to index on rows
 *   col_id     <-> pointer to array of column ids related to the row index
 *
 * returns:
 *    a pointer to a created CSR matrix structure
 *----------------------------------------------------------------------------*/

static cs_matrix_struct_dist_t *
_create_struct_msr_from_msr(bool         transfer,
                            bool         ordered,
                            cs_lnum_t    n_rows,
                            cs_lnum_t    n_cols_ext,
                            cs_lnum_t  **row_index,
                            cs_lnum_t  **col_id)
{
  cs_matrix_struct_dist_t  *ms = NULL;

  /* Allocate and map */

  BFT_MALLOC(ms, 1, cs_matrix_struct_dist_t);

  ms->n_rows = n_rows;
  ms->n_cols_ext = n_cols_ext;

  _init_struct_csr_from_csr(&(ms->e),
                            false,
                            transfer,
                            ordered,
                            n_rows,
                            n_cols_ext,
                            row_index,
                            col_id);

  _init_struct_csr(&(ms->h), 0, 0);

  ms->h_row_id = NULL;

  return ms;
}

/*----------------------------------------------------------------------------
 * Create a distributed matrix structure from a CSR representation.
 *
 * parameters:
 *   direct_assembly <-- true if each value corresponds to a unique face
 *   alloc_mode      <-- memory allocation mode
 *   n_rows          <- local number of rows
 *   n_cols_ext      <-- local number of columns + ghosts
 *   row_index       <-- index on rows
 *   col_id          <-> array of column ids related to the row index
 *
 * returns:
 *    a pointer to a created distrubuted (MSR configured) matrix structure
 *----------------------------------------------------------------------------*/

static cs_matrix_struct_dist_t *
_create_struct_msr_from_csr(bool              direct_assembly,
                            cs_alloc_mode_t   alloc_mode,
                            cs_lnum_t         n_rows,
                            cs_lnum_t         n_cols_ext,
                            const cs_lnum_t  *row_index,
                            const cs_lnum_t  *col_id)
{
  cs_matrix_struct_dist_t  *ms = NULL;

  /* Allocate and map */

  BFT_MALLOC(ms, 1, cs_matrix_struct_dist_t);

  ms->n_rows = n_rows;
  ms->n_cols_ext = n_cols_ext;

  ms->n_rows = n_rows;
  ms->n_cols_ext = n_cols_ext;

  _init_struct_csr(&(ms->e), n_rows, n_cols_ext);
  _init_struct_csr(&(ms->h), 0, 0);

  cs_lnum_t *_row_index, *_col_id;
  CS_MALLOC_HD(_row_index, n_rows + 1, cs_lnum_t, alloc_mode);
  CS_MALLOC_HD(_col_id, row_index[n_rows], cs_lnum_t, alloc_mode);
  _row_index[0] = 0;
  cs_lnum_t k = 0;
  for (cs_lnum_t i = 0; i < n_rows; i++) {
    cs_lnum_t n_cols = row_index[i+1] - row_index[i];
    const cs_lnum_t *s_c_id = col_id + row_index[i];
    for (cs_lnum_t j = 0; j < n_cols; j++) {
      if (s_c_id[j] != i)
        _col_id[k++] = s_c_id[j];
    }
    _row_index[i+1] = k;
  }
  CS_REALLOC_HD(_col_id, _row_index[n_rows], cs_lnum_t, alloc_mode);

  ms->e._row_index = _row_index;
  ms->e._col_id = _col_id;

  ms->e.row_index = ms->e._row_index;
  ms->e.col_id = ms->e._col_id;

  ms->e.direct_assembly = direct_assembly;

  ms->h_row_id = NULL;

  return ms;
}

/*----------------------------------------------------------------------------
 * Create a distributed matrix structure from a CSR representation.
 *
 * parameters:
 *   direct_assembly <-- true if each value corresponds to a unique face
 *   alloc_mode      <-- memory allocation mode
 *   n_rows          <- local number of rows
 *   n_cols_ext      <-- local number of columns + ghosts
 *   row_index       <-- index on rows
 *   col_id          <-> array of column ids related to the row index
 *
 * returns:
 *    a pointer to a created CSR matrix structure
 *----------------------------------------------------------------------------*/

static cs_matrix_struct_dist_t *
_create_struct_dist_from_csr(bool              direct_assembly,
                             cs_alloc_mode_t   alloc_mode,
                             cs_lnum_t         n_rows,
                             cs_lnum_t         n_cols_ext,
                             const cs_lnum_t  *row_index,
                             const cs_lnum_t  *col_id)
{
  cs_matrix_struct_dist_t  *ms = NULL;

  /* Allocate and map */

  BFT_MALLOC(ms, 1, cs_matrix_struct_dist_t);

  ms->n_rows = n_rows;
  ms->n_cols_ext = n_cols_ext;

  ms->n_rows = n_rows;
  ms->n_cols_ext = n_cols_ext;

  _init_struct_csr(&(ms->e), n_rows, n_rows);
  _init_struct_csr(&(ms->h), n_rows, n_cols_ext);

  cs_lnum_t *_e_row_index, *_e_col_id, *_h_row_index, *_h_col_id;
  CS_MALLOC_HD(_e_row_index, n_rows + 1, cs_lnum_t, alloc_mode);
  CS_MALLOC_HD(_h_row_index, n_rows + 1, cs_lnum_t, alloc_mode);

  _e_row_index[0] = 0;
  _h_row_index[0] = 0;

  cs_lnum_t k = 0, l = 0;
  for (cs_lnum_t i = 0; i < n_rows; i++) {
    cs_lnum_t n_cols = row_index[i+1] - row_index[i];
    const cs_lnum_t *s_c_id = col_id + row_index[i];
    for (cs_lnum_t j = 0; j < n_cols; j++) {
      if (s_c_id[j] < n_rows)
        k++;
      else if (s_c_id[j] != i)
        l++;
    }
    _e_row_index[i+1] = k;
    _h_row_index[i+1] = l;
  }

  CS_MALLOC_HD(_e_col_id, _e_row_index[n_rows], cs_lnum_t, alloc_mode);
  CS_MALLOC_HD(_h_col_id, _h_row_index[n_rows], cs_lnum_t, alloc_mode);

  k = 0; l = 0;
  for (cs_lnum_t i = 0; i < n_rows; i++) {
    cs_lnum_t n_cols = row_index[i+1] - row_index[i];
    const cs_lnum_t *s_c_id = col_id + row_index[i];
    for (cs_lnum_t j = 0; j < n_cols; j++) {
      if (s_c_id[j] < n_rows)
        _e_col_id[k++] = s_c_id[j];
      else if (s_c_id[j] != i)
        _h_col_id[l++] = s_c_id[j];
    }
  }

  ms->e._row_index = _e_row_index;
  ms->e._col_id = _e_col_id;
  ms->e.row_index = ms->e._row_index;
  ms->e._col_id = ms->e._col_id;
  ms->e.direct_assembly = direct_assembly;

  ms->h._row_index = _h_row_index;
  ms->h._col_id = _h_col_id;
  ms->h.row_index = ms->h._row_index;
  ms->h._col_id = ms->h._col_id;
  ms->h.direct_assembly = direct_assembly;

  ms->h_row_id = NULL;

  return ms;
}

/*----------------------------------------------------------------------------
 * Create a distributed matrix structure from a MSR representation.
 *
 * parameters:
 *   direct_assembly <-- true if each value corresponds to a unique face
 *   alloc_mode      <-- memory allocation mode
 *   n_rows          <- local number of rows
 *   n_cols_ext      <-- local number of columns + ghosts
 *   row_index       <-- index on rows
 *   col_id          <-> array of colum ids related to the row index
 *
 * returns:
 *    a pointer to a created CSR matrix structure
 *----------------------------------------------------------------------------*/

static cs_matrix_struct_dist_t *
_create_struct_dist_from_msr(bool              direct_assembly,
                             cs_alloc_mode_t   alloc_mode,
                             cs_lnum_t         n_rows,
                             cs_lnum_t         n_cols_ext,
                             const cs_lnum_t  *row_index,
                             const cs_lnum_t  *col_id)
{
  cs_matrix_struct_dist_t  *ms = NULL;

  /* Allocate and map */

  BFT_MALLOC(ms, 1, cs_matrix_struct_dist_t);

  ms->n_rows = n_rows;
  ms->n_cols_ext = n_cols_ext;

  ms->n_rows = n_rows;
  ms->n_cols_ext = n_cols_ext;

  _init_struct_csr(&(ms->e), n_rows, n_rows);
  _init_struct_csr(&(ms->h), n_rows, n_cols_ext);

  cs_lnum_t *_e_row_index, *_e_col_id, *_h_row_index, *_h_col_id;
  CS_MALLOC_HD(_e_row_index, n_rows + 1, cs_lnum_t, alloc_mode);
  CS_MALLOC_HD(_h_row_index, n_rows + 1, cs_lnum_t, alloc_mode);

  _e_row_index[0] = 0;
  _h_row_index[0] = 0;

  cs_lnum_t k = 0, l = 0;
  for (cs_lnum_t i = 0; i < n_rows; i++) {
    cs_lnum_t n_cols = row_index[i+1] - row_index[i];
    const cs_lnum_t *s_c_id = col_id + row_index[i];
    for (cs_lnum_t j = 0; j < n_cols; j++) {
      if (s_c_id[j] < n_rows)
        k++;
      else
        l++;
    }
    _e_row_index[i+1] = k;
    _h_row_index[i+1] = l;
  }

  CS_MALLOC_HD(_e_col_id, _e_row_index[n_rows], cs_lnum_t, alloc_mode);
  CS_MALLOC_HD(_h_col_id, _h_row_index[n_rows], cs_lnum_t, alloc_mode);

  k = 0; l = 0;
  for (cs_lnum_t i = 0; i < n_rows; i++) {
    cs_lnum_t n_cols = row_index[i+1] - row_index[i];
    const cs_lnum_t *s_c_id = col_id + row_index[i];
    for (cs_lnum_t j = 0; j < n_cols; j++) {
      if (s_c_id[j] < n_rows)
        _e_col_id[k++] = s_c_id[j];
      else
        _h_col_id[l++] = s_c_id[j];
    }
  }

  ms->e._row_index = _e_row_index;
  ms->e._col_id = _e_col_id;
  ms->e.row_index = ms->e._row_index;
  ms->e._col_id = ms->e._col_id;
  ms->e.direct_assembly = direct_assembly;

  ms->h._row_index = _h_row_index;
  ms->h._col_id = _h_col_id;
  ms->h.row_index = ms->h._row_index;
  ms->h._col_id = ms->h._col_id;
  ms->h.direct_assembly = direct_assembly;

  ms->h_row_id = NULL;

  return ms;
}

/*----------------------------------------------------------------------------
 * Destroy a distributed matrix structure.
 *
 * parameters:
 *   ms  <->  pointer to distributed matrix structure pointer
 *----------------------------------------------------------------------------*/

static void
_destroy_struct_dist(void  **ms)
{
  if (ms != NULL && *ms !=NULL) {
    cs_matrix_struct_dist_t  *_ms = *ms;

    CS_FREE_HD(_ms->e._row_index);
    CS_FREE_HD(_ms->e._col_id);
    CS_FREE_HD(_ms->h._row_index);
    CS_FREE_HD(_ms->h._col_id);

    CS_FREE_HD(_ms->h_row_id);

    BFT_FREE(_ms);

    *ms= NULL;
  }
}

/*----------------------------------------------------------------------------
 * Create distributed matrix coefficients.
 *
 * returns:
 *   pointer to allocated coefficients structure.
 *----------------------------------------------------------------------------*/

static cs_matrix_coeff_dist_t *
_create_coeff_dist(void)
{
  cs_matrix_coeff_dist_t  *mc;

  /* Allocate */

  BFT_MALLOC(mc, 1, cs_matrix_coeff_dist_t);

  /* Initialize */

  mc->symmetric = false;

  mc->db_size = 0;
  mc->eb_size = 0;

  mc->d_val = NULL;
  mc->e_val = NULL;
  mc->h_val = NULL;

  mc->_d_val = NULL;
  mc->_e_val = NULL;
  mc->_h_val = NULL;

  mc->d_idx = NULL;

  return mc;
}

/*----------------------------------------------------------------------------
 * Destroy distributed matrix coefficients.
 *
 * parameters:
 *   m  <->  pointer to matrix structure
 *----------------------------------------------------------------------------*/

static void
_destroy_coeff_dist(cs_matrix_t  *m)
{
  if (m->coeffs != NULL) {
    cs_matrix_coeff_dist_t  *mc = m->coeffs;

    CS_FREE_HD(mc->_h_val);
    CS_FREE(mc->_e_val);
    CS_FREE(mc->_d_val);
    CS_FREE_HD(mc->d_idx);

    BFT_FREE(m->coeffs);
  }
}

/*----------------------------------------------------------------------------
 * Set distributed extradiagonal matrix coefficients for the case where there are
 * multiple contributions to a given coefficient.
 *
 * The matrix coefficients should have been initialized (i.e. set to 0)
 * some before using this function.
 *
 * parameters:
 *   matrix      <-- pointer to matrix structure
 *   symmetric   <-- indicates if extradiagonal values are symmetric
 *   n_edges     <-- local number of graph edges
 *   edges       <-- edges (symmetric row <-> column) connectivity
 *   xa          <-- extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_e_coeffs_dist_increment(cs_matrix_t        *matrix,
                             bool                symmetric,
                             cs_lnum_t           n_edges,
                             const cs_lnum_2_t  *edges,
                             const cs_real_t    *restrict xa)
{
  cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_matrix_struct_dist_t  *ms = matrix->structure;

  const cs_matrix_struct_csr_t  *ms_e = &(ms->e);
  const cs_matrix_struct_csr_t  *ms_h = &(ms->h);

  /* Copy extra-diagonal values */

  assert(edges != NULL);

  const cs_lnum_t *restrict edges_p
    = (const cs_lnum_t *restrict)(edges);

  const cs_lnum_t xa_stride = (symmetric) ? 1 : 2;
  const cs_lnum_t xa_sj = xa_stride / 2;

  for (cs_lnum_t face_id = 0; face_id < n_edges; face_id++) {
    cs_lnum_t kk, ll;
    cs_lnum_t ii = *edges_p++;
    cs_lnum_t jj = *edges_p++;
    if (ii < ms->n_rows) {
      if (jj < ms->n_rows) {
        for (kk = ms_e->row_index[ii]; ms_e->col_id[kk] != jj; kk++);
        mc->_e_val[kk] += xa[xa_stride*face_id];
      }
      else {
        for (kk = ms_h->row_index[ii]; ms_h->col_id[kk] != jj; kk++);
        mc->_h_val[kk] += xa[xa_stride*face_id];
      }
    }
    if (jj < ms->n_rows) {
      if (ii < ms->n_rows) {
        for (ll = ms_e->row_index[jj]; ms_e->col_id[ll] != ii; ll++);
        mc->_e_val[ll] += xa[xa_stride*face_id + xa_sj];
      }
      else {
        for (ll = ms_h->row_index[jj]; ms_h->col_id[ll] != ii; ll++);
        mc->_h_val[ll] += xa[xa_stride*face_id + xa_sj];
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Set distributed extradiagonal matrix coefficients for the case where there are
 * multiple contributions to a given coefficient.
 *
 * The matrix coefficients should have been initialized (i.e. set to 0)
 * some before using this function.
 *
 * parameters:
 *   matrix      <-- pointer to matrix structure
 *   symmetric   <-- indicates if extradiagonal values are symmetric
 *   n_edges     <-- local number of graph edges
 *   edges       <-- edges (symmetric row <-> column) connectivity
 *   xa          <-- extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_e_coeffs_dist_increment_block(cs_matrix_t        *matrix,
                                   bool                symmetric,
                                   cs_lnum_t           n_edges,
                                   const cs_lnum_2_t  *edges,
                                   const cs_real_t    *restrict xa)
{
  cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_matrix_struct_dist_t  *ms = matrix->structure;

  const cs_lnum_t  b_size_2 = matrix->eb_size * matrix->eb_size;

  const cs_matrix_struct_csr_t  *ms_e = &(ms->e);
  const cs_matrix_struct_csr_t  *ms_h = &(ms->h);

  /* Copy extra-diagonal values */

  assert(edges != NULL);

  const cs_lnum_t *restrict edges_p
    = (const cs_lnum_t *restrict)(edges);

  const cs_lnum_t xa_stride = (symmetric) ? b_size_2 : 2*b_size_2;
  const cs_lnum_t xa_sj = (symmetric) ? 0 : b_size_2;

  for (cs_lnum_t face_id = 0; face_id < n_edges; face_id++) {
    cs_lnum_t ii = *edges_p++;
    cs_lnum_t jj = *edges_p++;
    cs_lnum_t kk, ll;
    if (ii < ms->n_rows) {
      if (jj < ms->n_rows) {
        for (kk = ms_e->row_index[ii]; ms_e->col_id[kk] != jj; kk++);
        for (cs_lnum_t pp = 0; pp < b_size_2; pp++)
          mc->_e_val[kk*b_size_2 + pp] += xa[xa_stride*face_id + pp];
      }
      else {
        for (kk = ms_h->row_index[ii]; ms_h->col_id[kk] != jj; kk++);
        for (cs_lnum_t pp = 0; pp < b_size_2; pp++)
          mc->_h_val[kk*b_size_2 + pp] += xa[xa_stride*face_id + pp];
      }
    }
    if (jj < ms->n_rows) {
      if (ii < ms->n_rows) {
        for (ll = ms_e->row_index[jj]; ms_e->col_id[ll] != ii; ll++);
        for (cs_lnum_t pp = 0; pp < b_size_2; pp++)
          mc->_e_val[ll*b_size_2 + pp] += xa[xa_stride*face_id + xa_sj + pp];
      }
      else {
        for (ll = ms_h->row_index[jj]; ms_h->col_id[ll] != ii; ll++);
        for (cs_lnum_t pp = 0; pp < b_size_2; pp++)
          mc->_h_val[ll*b_size_2 + pp] += xa[xa_stride*face_id + xa_sj + pp];
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Map or copy distributed matrix diagonal coefficients.
 *
 * parameters:
 *   matrix           <-> pointer to matrix structure
 *   copy             <-- indicates if coefficients should be copied
 *   da               <-- diagonal values (NULL if all zero)
 *----------------------------------------------------------------------------*/

static void
_map_or_copy_d_coeffs_dist(cs_matrix_t      *matrix,
                           bool              copy,
                           const cs_real_t  *restrict da)
{
  cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_lnum_t n_rows = matrix->n_rows;
  const cs_lnum_t db_size = matrix->db_size;
  const cs_lnum_t db_size_2 = matrix->db_size * matrix->db_size;

  /* Map or copy diagonal values */

  if (da != NULL) {

    if (copy) {
      CS_FREE_HD(mc->_d_val);
      CS_MALLOC_HD(mc->_d_val, db_size_2*n_rows, cs_real_t, matrix->alloc_mode);
#     pragma omp parallel for  if(n_rows*db_size > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
        for (cs_lnum_t jj = 0; jj < db_size_2; jj++)
          mc->_d_val[ii*db_size_2 + jj] = da[ii*db_size_2 + jj];
      }
      mc->d_val = mc->_d_val;
    }
    else
      mc->d_val = da;

  }
  else
    mc->d_val = NULL;

  mc->db_size = db_size;
}

/*----------------------------------------------------------------------------
 * Set distributed matrix coefficients.
 *
 * parameters:
 *   matrix      <-> pointer to matrix structure
 *   symmetric   <-- indicates if extradiagonal values are symmetric
 *   copy        <-- indicates if coefficients should be copied
 *   n_edges     <-- local number of graph edges
 *   edges       <-- edges (symmetric row <-> column) connectivity
 *   da          <-- diagonal values (NULL if all zero)
 *   xa          <-- extradiagonal values (NULL if all zero)
 *----------------------------------------------------------------------------*/

static void
_set_coeffs_dist(cs_matrix_t         *matrix,
                 bool                 symmetric,
                 bool                 copy,
                 cs_lnum_t            n_edges,
                 const cs_lnum_2_t  *restrict edges,
                 const cs_real_t    *restrict da,
                 const cs_real_t    *restrict xa)
{
  cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_lnum_t eb_size = matrix->eb_size;
  const cs_lnum_t eb_size_2 = matrix->eb_size * matrix->eb_size;

  /* Map or copy diagonal values */

  _map_or_copy_d_coeffs_dist(matrix, copy, da);

  /* Extradiagonal values */

  mc->eb_size = matrix->eb_size;

  CS_FREE_HD(mc->_e_val);
  CS_FREE_HD(mc->_h_val);
  CS_MALLOC_HD(mc->_e_val, eb_size_2*ms->e.row_index[ms->e.n_rows], cs_real_t,
               matrix->alloc_mode);
  if (ms->h.n_rows > 0)
    CS_MALLOC_HD(mc->_h_val, eb_size_2*ms->h.row_index[ms->h.n_rows], cs_real_t,
                 matrix->alloc_mode);
  else
    mc->_h_val = NULL;

  mc->e_val = mc->_e_val;
  mc->h_val = mc->_h_val;

  /* No direct assembly optimization (yet) */

  _zero_coeffs_csr(&(ms->e), matrix->eb_size, mc->_e_val);
  _zero_coeffs_csr(&(ms->h), matrix->eb_size, mc->_h_val);
  if (eb_size == 1)
    _set_e_coeffs_dist_increment(matrix, symmetric, n_edges, edges, xa);
  else
    _set_e_coeffs_dist_increment_block(matrix, symmetric, n_edges, edges, xa);
}

/*----------------------------------------------------------------------------
 * Release shared distributed matrix coefficients.
 *
 * parameters:
 *   matrix <-- pointer to matrix structure
 *----------------------------------------------------------------------------*/

static void
_release_coeffs_dist(cs_matrix_t  *matrix)
{
  CS_UNUSED(matrix);

  /* No shared values in distributed coefficients */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get matrix diagonal values for distributed matrix.
 *
 * In case of matrixes with block diagonal coefficients, a pointer to
 * the complete block diagonal is returned.
 *
 * \param[in]  matrix  pointer to matrix structure
 *
 * \return  pointer to matrix diagonal array
 */
/*----------------------------------------------------------------------------*/

static const cs_real_t *
_get_diagonal_dist(const cs_matrix_t  *matrix)
{
  const cs_real_t  *diag = NULL;

  cs_matrix_coeff_dist_t *mc = matrix->coeffs;
  if (mc->d_val == NULL) {
    cs_lnum_t n_rows = matrix->n_rows * matrix->db_size;
    if (mc->_d_val == NULL) {
      cs_lnum_t db_size_2 = matrix->db_size * matrix->db_size;
      CS_MALLOC_HD(mc->_d_val, db_size_2*matrix->n_rows, cs_real_t,
                   matrix->alloc_mode);
      mc->db_size = matrix->db_size;
    }
#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      mc->_d_val[ii] = 0.0;
    mc->d_val = mc->_d_val;
  }
  diag = mc->d_val;

  return diag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create matrix structure internals using a matrix assembler.
 *
 * Only CSR and MSR formats are handled.
 *
 * \param[in]  type  type of matrix considered
 * \param[in]  ma    pointer to matrix assembler structure
 *
 * \return  a pointer to created matrix structure internals
 */
/*----------------------------------------------------------------------------*/

static void *
_structure_from_assembler(cs_matrix_type_t        type,
                          cs_lnum_t               n_rows,
                          cs_lnum_t               n_cols_ext,
                          cs_matrix_assembler_t  *ma)
{
  void *structure = NULL;

  /* Get info on assembler structure */

  bool             ma_sep_diag = cs_matrix_assembler_get_separate_diag(ma);
  cs_alloc_mode_t  alloc_mode = cs_alloc_mode;
  const cs_lnum_t *row_index = cs_matrix_assembler_get_row_index(ma);
  const cs_lnum_t *col_id = cs_matrix_assembler_get_col_ids(ma);

  /* Define structure */

  switch(type) {

  case CS_MATRIX_CSR:
    /* Assume diagonal is present (should not be important
       for assembly using matrix assembler) */
    if (ma_sep_diag == false)
      structure = _create_struct_csr_from_shared(true, /* have_diag */
                                                 false, /* for safety */
                                                 n_rows,
                                                 n_cols_ext,
                                                 row_index,
                                                 col_id);
    else {
      cs_lnum_t *_row_index, *_col_id;
      CS_MALLOC_HD(_row_index, n_rows + 1, cs_lnum_t, alloc_mode);
      CS_MALLOC_HD(_col_id, row_index[n_rows] + n_rows, cs_lnum_t, alloc_mode);
      _row_index[0] = 0;
      for (cs_lnum_t i = 0; i < n_rows; i++) {
        cs_lnum_t n_cols = row_index[i+1] - row_index[i];
        cs_lnum_t j = 0, k = 0;
        const cs_lnum_t *s_c_id = col_id + row_index[i];
        cs_lnum_t *d_c_id = _col_id + row_index[i] + i;
        while (j < n_cols && s_c_id[j] < i)
          d_c_id[k++] = s_c_id[j++];
        d_c_id[k++] = i;
        while (j < n_cols)
          d_c_id[k++] = s_c_id[j++];
        _row_index[i+1] = row_index[i+1] + i + 1;
      }
      structure = _create_struct_csr_from_csr(true, /* have_idag */
                                              true,
                                              true,
                                              n_rows,
                                              n_cols_ext,
                                              &_row_index,
                                              &_col_id);
    }
    break;

  case CS_MATRIX_MSR:
    if (ma_sep_diag == true)
      structure = _create_struct_msr_from_shared(false, /* for safety */
                                                 n_rows,
                                                 n_cols_ext,
                                                 row_index,
                                                 col_id);
    else {
      structure = _create_struct_msr_from_csr(true,
                                              alloc_mode,
                                              n_rows,
                                              n_cols_ext,
                                              row_index,
                                              col_id);
    }
    break;

  case CS_MATRIX_DIST:
    if (ma_sep_diag == true)
      structure = _create_struct_dist_from_msr(false, /* for safety */
                                               alloc_mode,
                                               n_rows,
                                               n_cols_ext,
                                               row_index,
                                               col_id);
    else {
      structure = _create_struct_dist_from_csr(true,
                                               alloc_mode,
                                               n_rows,
                                               n_cols_ext,
                                               row_index,
                                               col_id);
    }
    break;

  default:
    if (type >= 0 && type < CS_MATRIX_N_BUILTIN_TYPES)
      bft_error(__FILE__, __LINE__, 0,
                _("%s: handling of matrices in %s format\n"
                  "is not operational yet."),
                __func__,
                _(_matrix_type_name[type]));
    else
      bft_error(__FILE__, __LINE__, 0,
                _("%s: handling of matrices in external format type %d\n"
                  "is not handled by this function."),
                __func__, (int)type);
    break;
  }

  return structure;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy matrix structure internals.
 *
 * \param[in]       type       matrix structure type
 * \param[in, out]  structure  pointer to matrix structure pointer
 */
/*----------------------------------------------------------------------------*/

static void
_structure_destroy(cs_matrix_type_t   type,
                   void             **structure)
{
  switch(type) {
  case CS_MATRIX_NATIVE:
    _destroy_struct_native(structure);
    break;
  case CS_MATRIX_CSR:
    _destroy_struct_csr(structure);
    break;
  case CS_MATRIX_MSR:
    _destroy_struct_dist(structure);
    break;
  case CS_MATRIX_DIST:
    _destroy_struct_dist(structure);
    break;
  default:
    assert(0);
    break;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a matrix container using a given type.
 *
 * \param[in]  type  chosen matrix type
 *
 * \return  pointer to created matrix structure;
 */
/*----------------------------------------------------------------------------*/

static cs_matrix_t *
_matrix_create(cs_matrix_type_t  type)
{
  cs_matrix_fill_type_t mft;
  cs_matrix_t *m;

  BFT_MALLOC(m, 1, cs_matrix_t);

  m->type = type;

  if (m->type >= 0 && m->type < CS_MATRIX_N_BUILTIN_TYPES) {
    m->type_name = _matrix_type_name[m->type];
    m->type_fname = _matrix_type_fullname[m->type];
  }
  else {
    m->type_name = _matrix_type_name[CS_MATRIX_N_BUILTIN_TYPES];
    m->type_fname = _matrix_type_fullname[CS_MATRIX_N_BUILTIN_TYPES];
  }

  /* Map shared structure */

  m->n_rows = 0;
  m->n_cols_ext = 0;

  m->symmetric = false;

  m->db_size = 0;
  m->eb_size = 0;

  m->alloc_mode = cs_alloc_mode;
  /* Native and distributed matrix types not on accelerator yet */
  if (m->type == CS_MATRIX_NATIVE || m->type == CS_MATRIX_DIST)
    m->alloc_mode = CS_ALLOC_HOST;

  m->fill_type = CS_MATRIX_N_FILL_TYPES;

  m->structure = NULL;
  m->_structure = NULL;

  m->halo = NULL;
  m->numbering = NULL;
  m->assembler = NULL;

  for (mft = 0; mft < CS_MATRIX_N_FILL_TYPES; mft++) {
    for (cs_matrix_spmv_type_t i = 0; i < CS_MATRIX_SPMV_N_TYPES; i++) {
      m->vector_multiply[mft][i] = NULL;
#if defined(HAVE_ACCEL)
      m->vector_multiply_h[mft][i] = NULL;
      m->vector_multiply_d[mft][i] = NULL;
#endif
    }
  }

  /* Define coefficients */

  switch(m->type) {
  case CS_MATRIX_NATIVE:
    m->coeffs = _create_coeff_dist();
    break;
  case CS_MATRIX_CSR:
    m->coeffs = _create_coeff_csr();
    break;
  case CS_MATRIX_MSR:
    m->coeffs = _create_coeff_dist();
    break;
  case CS_MATRIX_DIST:
    m->coeffs = _create_coeff_dist();
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Handling of matrixes in format type %d\n"
                "is not operational yet."),
              m->type);
    break;
  }

  m->xa = NULL;

  m->c2f_idx = NULL;
  m->c2f = NULL;
  m->c2f_sgn = NULL;

  m->cell_cen = NULL;
  m->cell_vol = NULL;
  m->face_normal = NULL;

  /* Mapping to external libraries */

  m->ext_lib_map = NULL;

  /* Set function pointers here */

  m->set_coefficients = NULL;
  m->destroy_adaptor = NULL;

  cs_matrix_spmv_set_defaults(m);

  switch(m->type) {

  case CS_MATRIX_NATIVE:

    m->set_coefficients = _set_coeffs_native;
    m->release_coefficients = _release_coeffs_dist;
    m->copy_diagonal = _copy_diagonal_separate;
    m->get_diagonal = _get_diagonal_dist;
    m->destroy_structure = _destroy_struct_native;
    m->destroy_coefficients = _destroy_coeff_dist;
    m->assembler_values_create = NULL;
 break;

  case CS_MATRIX_CSR:
    m->set_coefficients = _set_coeffs_csr;
    m->release_coefficients = _release_coeffs_csr;
    m->copy_diagonal = _copy_diagonal_csr;
    m->get_diagonal = _get_diagonal_csr;
    m->destroy_structure = _destroy_struct_csr;
    m->destroy_coefficients = _destroy_coeff_csr;
    m->assembler_values_create = _assembler_values_create_csr;
    break;

  case CS_MATRIX_MSR:
    m->set_coefficients = _set_coeffs_msr;
    m->release_coefficients = _release_coeffs_dist;
    m->copy_diagonal = _copy_diagonal_separate;
    m->get_diagonal = _get_diagonal_dist;
    m->destroy_structure = _destroy_struct_dist;
    m->destroy_coefficients = _destroy_coeff_dist;
    m->assembler_values_create = _assembler_values_create_msr;
    break;

  case CS_MATRIX_DIST:
    m->set_coefficients = _set_coeffs_dist;
    m->release_coefficients = _release_coeffs_dist;
    m->copy_diagonal = _copy_diagonal_separate;
    m->get_diagonal = _get_diagonal_dist;
    m->destroy_structure = _destroy_struct_dist;
    m->destroy_coefficients = _destroy_coeff_dist;
    m->assembler_values_create = _assembler_values_create_dist;
    break;

  default:
    assert(0);
    break;

  }

  for (int i = 0; i < CS_MATRIX_N_FILL_TYPES; i++) {
    if (m->vector_multiply[i][1] == NULL)
      m->vector_multiply[i][1] = m->vector_multiply[i][0];
#if defined(HAVE_ACCEL)
    if (m->vector_multiply_h[i][1] == NULL)
      m->vector_multiply_h[i][1] = m->vector_multiply_h[i][0];
    if (m->vector_multiply_d[i][1] == NULL)
      m->vector_multiply_d[i][1] = m->vector_multiply_d[i][0];
#endif
  }

  return m;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Semi private function definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a matrix structure.
 *
 * Note that the structure created usually maps to the given existing
 * cell global number, face -> cell connectivity arrays, and cell halo
 * structure, so it must be destroyed before they are freed
 * (usually along with the code's main face -> cell structure).
 *
 * Note that the resulting matrix structure will contain a full main
 * main diagonal, and that the extra-diagonal structure is always
 * symmetric (though the coefficients my not be, and we may choose a
 * matrix format that does not exploit this symmetry). If the edges
 * connectivity argument is NULL, the matrix will be purely diagonal.
 *
 * \param[in]  type        type of matrix considered
 * \param[in]  n_rows      local number of rows
 * \param[in]  n_cols_ext  number of local + ghost columns
 * \param[in]  n_edges     local number of (undirected) graph edges
 * \param[in]  edges       edges (symmetric row <-> column) connectivity
 * \param[in]  halo        halo structure associated with cells, or NULL
 * \param[in]  numbering   vectorization or thread-related numbering info,
 *                         or NULL
 *
 * \return  pointer to created matrix structure;
 */
/*----------------------------------------------------------------------------*/

cs_matrix_structure_t *
cs_matrix_structure_create(cs_matrix_type_t       type,
                           cs_lnum_t              n_rows,
                           cs_lnum_t              n_cols_ext,
                           cs_lnum_t              n_edges,
                           const cs_lnum_2_t     *edges,
                           const cs_halo_t       *halo,
                           const cs_numbering_t  *numbering)
{
  cs_matrix_structure_t *ms;

  BFT_MALLOC(ms, 1, cs_matrix_structure_t);

  ms->type = type;

  ms->n_rows = n_rows;
  ms->n_cols_ext = n_cols_ext;

  ms->alloc_mode = cs_alloc_mode;

  /* Define Structure */

  switch(ms->type) {
  case CS_MATRIX_NATIVE:
    ms->structure = _create_struct_native(n_rows,
                                          n_cols_ext,
                                          n_edges,
                                          edges);
    break;
  case CS_MATRIX_CSR:
    ms->structure = _create_struct_csr(ms->alloc_mode,
                                       n_rows,
                                       n_cols_ext,
                                       n_edges,
                                       edges);
    break;
  case CS_MATRIX_MSR:
    ms->structure = _create_struct_msr(ms->alloc_mode,
                                       n_rows,
                                       n_cols_ext,
                                       n_edges,
                                       edges);
    break;
  case CS_MATRIX_DIST:
    ms->structure = _create_struct_dist(ms->alloc_mode,
                                        n_rows,
                                        n_cols_ext,
                                        n_edges,
                                        edges);
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Handling of matrixes in format type %d\n"
                "is not operational yet."),
              type);
    break;
  }

  /* Set pointers to structures shared from mesh here */

  ms->halo = halo;
  ms->numbering = numbering;
  ms->assembler = NULL;

  return ms;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a matrix structure based on a MSR connectivity definition.
 *
 * Only CSR and MSR formats are handled.
 *
 * col_id is sorted row by row during the creation of this structure.
 *
 * In case the property of the row index and col_id arrays are transferred
 * to the structure, the arrays pointers passed as arguments are set to NULL,
 * to help ensure the caller does not use the original arrays directly after
 * this call.
 *
 * \param[in]       type        type of matrix considered
 * \param[in]       transfer    transfer property of row_index and col_id
 *                              if true, map them otherwise
 * \param[in]       have_diag   indicates if the structure includes the
 *                              diagonal (should be the same for all rows)
 * \param[in]       n_rows      local number of rows
 * \param[in]       n_cols_ext  local number of columns + ghosts
 * \param[in]       row_index   pointer to index on rows
 * \param[in, out]  col_id      pointer to array of column ids related to
 *                              the row index
 * \param[in]       halo        halo structure for synchronization, or NULL
 * \param[in]       numbering   vectorization or thread-related numbering info,
 *                              or NULL
 *
 * \return  a pointer to a created matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_matrix_structure_t *
cs_matrix_structure_create_msr(cs_matrix_type_t        type,
                               bool                    transfer,
                               bool                    have_diag,
                               cs_lnum_t               n_rows,
                               cs_lnum_t               n_cols_ext,
                               cs_lnum_t             **row_index,
                               cs_lnum_t             **col_id,
                               const cs_halo_t        *halo,
                               const cs_numbering_t   *numbering)
{
  cs_matrix_structure_t *ms = NULL;

  BFT_MALLOC(ms, 1, cs_matrix_structure_t);

  ms->type = type;

  ms->n_rows = n_rows;
  ms->n_cols_ext = n_cols_ext;

  /* Define Structure */

  switch(ms->type) {
  case CS_MATRIX_CSR:
    ms->structure = _create_struct_csr_from_csr(have_diag,
                                                transfer,
                                                false,
                                                n_rows,
                                                n_cols_ext,
                                                row_index,
                                                col_id);
    break;
  case CS_MATRIX_MSR:
    ms->structure = _create_struct_msr_from_msr(transfer,
                                                false,
                                                n_rows,
                                                n_cols_ext,
                                                row_index,
                                                col_id);
    break;
  default:
    if (type >= 0 && type < CS_MATRIX_N_BUILTIN_TYPES)
      bft_error(__FILE__, __LINE__, 0,
                _("%s: handling of matrices in %s format\n"
                  "is not operational yet."),
                __func__,
                _(_matrix_type_name[type]));
    else
      bft_error(__FILE__, __LINE__, 0,
                _("%s: handling of matrices in external format type %d\n"
                  "is not handled by this function."),
                __func__, (int)type);
    break;
  }

  /* Set pointers to structures shared from mesh here */

  ms->halo = halo;
  ms->numbering = numbering;
  ms->assembler = NULL;

  return ms;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create an MSR matrix structure sharing an existing connectivity
 * definition.
 *
 * Note that as the structure created maps to the given existing
 * cell global number, face -> cell connectivity arrays, and cell halo
 * structure, it must be destroyed before they are freed
 * (usually along with the code's main face -> cell structure).
 *
 * \param[in]  direct_assembly  true if each value corresponds to
                                a unique face
 * \param[in]  n_rows           local number of rows
 * \param[in]  n_cols_ext       local number of columns + ghosts
 * \param[in]  row_index        index on rows
 * \param[in]  col_id           array of column ids related to the row index
 * \param[in]  halo             halo structure for synchronization, or NULL
 * \param[in]  numbering        vectorization or thread-related numbering
 *                              info, or NULL
 *
 * \returns  a pointer to a created matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_matrix_structure_t *
cs_matrix_structure_create_msr_shared(bool                    direct_assembly,
                                      cs_lnum_t               n_rows,
                                      cs_lnum_t               n_cols_ext,
                                      const cs_lnum_t        *row_index,
                                      const cs_lnum_t        *col_id,
                                      const cs_halo_t        *halo,
                                      const cs_numbering_t   *numbering)
{
  cs_matrix_structure_t *ms = NULL;

  BFT_MALLOC(ms, 1, cs_matrix_structure_t);

  ms->type = CS_MATRIX_MSR;

  ms->n_rows = n_rows;
  ms->n_cols_ext = n_cols_ext;

  ms->alloc_mode = cs_check_device_ptr(col_id);

  /* Define Structure */

  ms->structure = _create_struct_msr_from_shared(direct_assembly,
                                                 n_rows,
                                                 n_cols_ext,
                                                 row_index,
                                                 col_id);

  /* Set pointers to structures shared from mesh here */

  ms->halo = halo;
  ms->numbering = numbering;
  ms->assembler = NULL;

  return ms;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a matrix structure using a matrix assembler.
 *
 * Only CSR and MSR formats are handled.
 *
 * \param[in]  type  type of matrix considered
 * \param[in]  ma    pointer to matrix assembler structure
 *
 * \return  a pointer to a created matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_matrix_structure_t *
cs_matrix_structure_create_from_assembler(cs_matrix_type_t        type,
                                          cs_matrix_assembler_t  *ma)
{
  cs_matrix_structure_t *ms = NULL;

  BFT_MALLOC(ms, 1, cs_matrix_structure_t);

  ms->type = type;

  ms->n_rows = cs_matrix_assembler_get_n_rows(ma);
  ms->n_cols_ext = cs_matrix_assembler_get_n_columns(ma);;

  ms->alloc_mode = cs_alloc_mode;

  /* Define internal structure */

  ms->structure = _structure_from_assembler(ms->type,
                                            ms->n_rows,
                                            ms->n_cols_ext,
                                            ma);

  /* Set pointers to structures shared from mesh here */

  ms->halo = cs_matrix_assembler_get_halo(ma);

  ms->numbering = NULL;

  ms->assembler = ma;

  return ms;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a matrix structure.
 *
 * \param[in, out]  ms  pointer to matrix structure pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_structure_destroy(cs_matrix_structure_t  **ms)
{
  if (ms != NULL && *ms != NULL) {

    cs_matrix_structure_t *_ms = *ms;

    _structure_destroy(_ms->type, &(_ms->structure));

    /* Now free main structure */

    BFT_FREE(*ms);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a matrix container using a given structure.
 *
 * Note that the matrix container maps to the assigned structure,
 * so it must be destroyed before that structure.
 *
 * \param[in]  ms  associated matrix structure
 *
 * \return  pointer to created matrix structure;
 */
/*----------------------------------------------------------------------------*/

cs_matrix_t *
cs_matrix_create(const cs_matrix_structure_t  *ms)
{
  assert(ms != NULL); /* Sanity check */

  cs_matrix_t *m = _matrix_create(ms->type);

  /* Map shared structure */

  m->n_rows = ms->n_rows;
  m->n_cols_ext = ms->n_cols_ext;

  m->structure = ms->structure;

  m->halo = ms->halo;
  m->numbering = ms->numbering;
  m->assembler = ms->assembler;

  return m;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a matrix directly from assembler.
 *
 * Only CSR and MSR formats are handled.
 *
 * \param[in]  type  type of matrix considered
 * \param[in]  ma    pointer to matrix assembler structure
 *
 * \return  a pointer to a created matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_matrix_t *
cs_matrix_create_from_assembler(cs_matrix_type_t        type,
                                cs_matrix_assembler_t  *ma)
{
  cs_matrix_t *m = _matrix_create(type);

  m->assembler = ma;

  m->type = type;

  m->n_rows = cs_matrix_assembler_get_n_rows(ma);
  m->n_cols_ext = cs_matrix_assembler_get_n_columns(ma);;

  /* Define internal structure */

  m->_structure = _structure_from_assembler(m->type,
                                            m->n_rows,
                                            m->n_cols_ext,
                                            ma);
  m->structure = m->_structure;

  /* Set pointers to structures shared from mesh here */

  m->halo = cs_matrix_assembler_get_halo(ma);

  m->numbering = NULL;

  m->assembler = ma;

  return m;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a matrix container by copying another
 *
 * Note that the matrix containers share the same assigned structure,
 * so they must be both destroyed before that structure.
 *
 * If assigned, coefficients are not copied.
 *
 * \param[in]  src  reference matrix structure
 *
 * \return  pointer to created matrix structure;
 */
/*----------------------------------------------------------------------------*/

cs_matrix_t *
cs_matrix_create_by_copy(cs_matrix_t   *src)
{
  cs_matrix_t *m;

  BFT_MALLOC(m, 1, cs_matrix_t);

  memcpy(m, src, sizeof(cs_matrix_t));

  /* Define coefficients */

  switch(m->type) {
  case CS_MATRIX_NATIVE:
    m->coeffs = _create_coeff_dist();
    break;
  case CS_MATRIX_CSR:
    m->coeffs = _create_coeff_csr();
    break;
  case CS_MATRIX_MSR:
    m->coeffs = _create_coeff_dist();
    break;
  case CS_MATRIX_DIST:
    m->coeffs = _create_coeff_dist();
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Handling of matrixes in format type %d\n"
                "is not operational yet."),
              m->type);
    break;
  }

  cs_matrix_release_coefficients(m);

  return m;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a matrix based on the local restriction of a base matrix.
 *
 * Coefficients are copied. Some coefficients may be shared with the
 * parent matrix, so the base matrix must not be destroyed before the
 * restriction matrix.
 *
 * \param[in]  src  reference matrix structure
 *
 * \return  pointer to created matrix structure;
 */
/*----------------------------------------------------------------------------*/

cs_matrix_t *
cs_matrix_create_by_local_restrict(const cs_matrix_t  *src)
{
  cs_matrix_t *m = NULL;

  const cs_lnum_t n_rows = src->n_rows;

  BFT_MALLOC(m, 1, cs_matrix_t);
  memcpy(m, src, sizeof(cs_matrix_t));
  m->n_cols_ext = m->n_rows;

  m->structure = NULL;
  m->_structure = NULL;

  m->halo = NULL;
  m->numbering = NULL;
  m->assembler = NULL;
  m->xa = NULL;
  m->coeffs = NULL;

  /* Define coefficients */

  cs_lnum_t eb_size_2 = src->eb_size * src->eb_size;

  switch(m->type) {
  case CS_MATRIX_MSR:
    {
      m->_structure = _create_struct_csr_from_restrict_local(src->structure);
      m->structure = m->_structure;
      m->coeffs = _create_coeff_dist();
      cs_matrix_coeff_dist_t  *mc = m->coeffs;
      cs_matrix_coeff_dist_t  *mc_src = src->coeffs;
      const cs_matrix_struct_csr_t *ms = m->structure;
      const cs_matrix_struct_csr_t *ms_src = src->structure;
      mc->d_val = mc_src->d_val;
      CS_MALLOC_HD(mc->_e_val, eb_size_2*ms->row_index[n_rows], cs_real_t,
                   src->alloc_mode);
      mc->e_val = mc->_e_val;
      for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
        const cs_lnum_t  n_cols = ms->row_index[ii+1] - ms->row_index[ii];
        const cs_real_t  *s_row =   mc_src->e_val
                                  + ms_src->row_index[ii]*eb_size_2;
        cs_real_t  *m_row = mc->_e_val + ms->row_index[ii]*eb_size_2;
        memcpy(m_row, s_row, sizeof(cs_real_t)*eb_size_2*n_cols);
      }
      mc->db_size = m->db_size;
      mc->eb_size = m->eb_size;
    }
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Handling of matrixes in %s format\n"
                "is not operational yet."),
              _(m->type_name));
    break;
  }

  return m;
}

/*----------------------------------------------------------------------------
 * Destroy a matrix structure.
 *
 * In the case of a compound matrix, sub-matrices are not destroyed.
 *
 * parameters:
 *   matrix <-> pointer to matrix structure pointer
 *----------------------------------------------------------------------------*/

void
cs_matrix_destroy(cs_matrix_t **matrix)
{
  if (matrix == NULL)
    return;

  cs_matrix_t *m = *matrix;

  if (m == NULL)
    return;

  m->destroy_coefficients(m);

  if (m->_structure != NULL) {
    m->destroy_structure(&(m->_structure));
    m->structure = NULL;
  }

  /* Now free main structure */

  BFT_FREE(*matrix);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return matrix type.
 *
 * \param[in]  matrix  pointer to matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_matrix_type_t
cs_matrix_get_type(const cs_matrix_t  *matrix)
{
  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0, _("The matrix is not defined."));
  return matrix->type;
}

/*----------------------------------------------------------------------------
 * Return matrix type name.
 *
 * parameters:
 *   matrix --> pointer to matrix structure
 *----------------------------------------------------------------------------*/

const char *
cs_matrix_get_type_name(const cs_matrix_t  *matrix)
{
  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0, _("%s: matrix not defined."), __func__);

  return matrix->type_name;
}

/*----------------------------------------------------------------------------
 * Return matrix type full name.
 *
 * parameters:
 *   matrix --> pointer to matrix structure
 *----------------------------------------------------------------------------*/

const char *
cs_matrix_get_type_fullname(const cs_matrix_t  *matrix)
{
  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0, _("%s: matrix not defined."), __func__);

  return matrix->type_fname;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of columns in a matrix.
 *
 * \param[in] matrix  pointer to matrix structure
 *
 * \return the number of columns
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_matrix_get_n_columns(const cs_matrix_t  *matrix)
{
  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0, _("The matrix is not defined."));

  return matrix->n_cols_ext;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the number of rows in matrix.
 *
 * \param[in] matrix  pointer to matrix structure
 *
 * \return the number of rows
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_matrix_get_n_rows(const cs_matrix_t  *matrix)
{
  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0, _("The matrix is not defined."));

  return matrix->n_rows;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the number of entries in matrix.
 *
 * When the block size is > 1, the number reported is the number of
 * entry blocks, not individual entries.
 *
 * \param[in] matrix  pointer to matrix structure
 *
 * \return the number of entries
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_matrix_get_n_entries(const cs_matrix_t  *matrix)
{
  cs_lnum_t retval = 0;

  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0, _("The matrix is not defined."));

  switch(matrix->type) {
  case CS_MATRIX_NATIVE:
    {
      const cs_matrix_struct_native_t  *ms = matrix->structure;
      retval = ms->n_edges*2 + ms->n_rows;
    }
    break;
  case CS_MATRIX_CSR:
    {
      const cs_matrix_struct_csr_t  *ms = matrix->structure;
      retval = ms->row_index[ms->n_rows];
    }
    break;
  case CS_MATRIX_MSR:
    {
      const cs_matrix_struct_dist_t  *ms = matrix->structure;
      retval = ms->e.row_index[ms->n_rows] + ms->n_rows;
    }
    break;
  case CS_MATRIX_DIST:
    {
      const cs_matrix_struct_dist_t  *ms = matrix->structure;
      retval += ms->e.row_index[ms->e.n_rows];
      if (ms->h.n_rows > 0)
        retval += ms->h.row_index[ms->h.n_rows];
      retval += ms->n_rows;
    }
    break;
  default:
    break;
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the size of the diagonal block for the given matrix
 *
 * \param[in] matrix  pointer to matrix structure
 *
 * \return the of the diagonal block
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_matrix_get_diag_block_size(const cs_matrix_t  *matrix)
{
  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0, _("The matrix is not defined."));

  return matrix->db_size;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the size of the extra-diagonal block for the given matrix
 *
 * \param[in] matrix  pointer to matrix structure
 *
 * \return the size of the extra-diagonal block
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_matrix_get_extra_diag_block_size(const cs_matrix_t  *matrix)
{
  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0, _("The matrix is not defined."));

  return matrix->eb_size;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the pointer to the halo structure for the given matrix
 *
 * \param[in] matrix  pointer to matrix structure
 *
 * \return pointer to the associated halo structure
 */
/*----------------------------------------------------------------------------*/

const cs_halo_t *
cs_matrix_get_halo(const cs_matrix_t  *matrix)
{
  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0, _("The matrix is not defined."));

  return matrix->halo;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to local global row range associated with a
 *        matrix, if available
 *
 * \param[in] matrix   pointer to matrix structure
 *
 * \return pointer to local range, or NULL
 */
/*----------------------------------------------------------------------------*/

const cs_gnum_t *
cs_matrix_get_l_range(const cs_matrix_t  *matrix)
{
  const cs_gnum_t *l_range = NULL;

  if (matrix != NULL) {
    if (matrix->assembler != NULL)
      l_range = cs_matrix_assembler_get_l_range(matrix->assembler);
  }

  return l_range;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Query matrix allocation mode.
 *
 * \param[in]  matrix  pointer to matrix structure
 *
 * \return  host/device allocation mode
 */
/*----------------------------------------------------------------------------*/

cs_alloc_mode_t
cs_matrix_get_alloc_mode(const cs_matrix_t  *matrix)
{
  return matrix->alloc_mode;
}

/*----------------------------------------------------------------------------*/
/*!
 *\brief Set matrix allocation mode.
 *
 * \param[in, out]  matrix      pointer to matrix structure
 * \param[in]       alloc_mode  host/device allocation mode
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_set_alloc_mode(cs_matrix_t       *matrix,
                         cs_alloc_mode_t   alloc_mode)
{
  matrix->alloc_mode = alloc_mode;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get matrix fill type, depending on block sizes.
 *
 * \param[in]  symmetric              indicates if matrix coefficients
 *                                    are symmetric
 * \param[in]  diag_block_size        block sizes for diagonal
 * \param[in]  extra_diag_block_size  block sizes for extra diagonal
 *
 * \return  matrix fill type
 */
/*----------------------------------------------------------------------------*/

cs_matrix_fill_type_t
cs_matrix_get_fill_type(bool       symmetric,
                        cs_lnum_t  diag_block_size,
                        cs_lnum_t  extra_diag_block_size)
{
  cs_matrix_fill_type_t fill_type = CS_MATRIX_N_FILL_TYPES;

  cs_lnum_t _db_size = diag_block_size;
  cs_lnum_t _eb_size = extra_diag_block_size;

  /* Set fill type */

  cs_base_check_bool(&symmetric);

  if (_db_size == 1) {
    if (symmetric)
      fill_type = CS_MATRIX_SCALAR_SYM;
    else
      fill_type = CS_MATRIX_SCALAR;
  }
  else if (_eb_size == 1) {
    if (symmetric)
      fill_type = CS_MATRIX_BLOCK_D_SYM;
    else if (_db_size == 6)
      fill_type = CS_MATRIX_BLOCK_D_66;
    else
      fill_type = CS_MATRIX_BLOCK_D;
  }
  else
    fill_type = CS_MATRIX_BLOCK;

  return fill_type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set matrix coefficients defined relative to a "native" edge graph,
 * sharing arrays with the caller when possible.
 *
 * With shared arrays, the matrix becomes unusable if the arrays passed as
 * arguments are not be modified (its coefficients should be unset first
 * to mark this).
 *
 * Depending on current options and initialization, values will be copied
 * or simply mapped.
 *
 * \param[in, out]  matrix                 pointer to matrix structure
 * \param[in]       symmetric              indicates if matrix coefficients
 *                                         are symmetric
 * \param[in]       diag_block_size        block sizes for diagonal
 * \param[in]       extra_diag_block_size  block sizes for extra diagonal
 * \param[in]       n_edges                local number of graph edges
 * \param[in]       edges                  edges (row <-> column) connectivity
 * \param[in]       da                     diagonal values (NULL if zero)
 * \param[in]       xa                     extradiagonal values (NULL if zero)
 *                                         casts as:
 *                                           xa[n_edges]    if symmetric,
 *                                           xa[n_edges][2] if non symmetric
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_set_coefficients(cs_matrix_t        *matrix,
                           bool                symmetric,
                           cs_lnum_t           diag_block_size,
                           cs_lnum_t           extra_diag_block_size,
                           const cs_lnum_t     n_edges,
                           const cs_lnum_2_t   edges[],
                           const cs_real_t    *da,
                           const cs_real_t    *xa)
{
  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0, _("The matrix is not defined."));

  cs_base_check_bool(&symmetric);

  /* Set the fill type */

  _set_fill_info(matrix,
                 symmetric,
                 diag_block_size,
                 extra_diag_block_size);

  /* Set the coefficients */

  if (matrix->set_coefficients != NULL) {
    matrix->xa = xa;
    matrix->set_coefficients(matrix, symmetric, false, n_edges, edges, da, xa);
  }
  else
    bft_error
      (__FILE__, __LINE__, 0,
       "Matrix format %s with fill type %s does not handle\n"
       "coefficient assignment from native (graph-edge) coefficients.",
       matrix->type_name,
       cs_matrix_fill_type_name[matrix->fill_type]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set matrix coefficients, copying values to private arrays.
 *
 * With private arrays, the matrix becomes independant from the
 * arrays passed as arguments.
 *
 * \param[in, out]  matrix                 pointer to matrix structure
 * \param[in]       symmetric              indicates if matrix coefficients
 *                                         are symmetric
 * \param[in]       diag_block_size        block sizes for diagonal
 * \param[in]       extra_diag_block_size  block sizes for extra diagonal
 * \param[in]       n_edges                local number of graph edges
 * \param[in]       edges                  edges (row <-> column) connectivity
 * \param[in]       da                     diagonal values (NULL if zero)
 * \param[in]       xa                     extradiagonal values (NULL if zero)
 *                                         casts as:
 *                                           xa[n_edges]    if symmetric,
 *                                           xa[n_edges][2] if non symmetric
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_copy_coefficients(cs_matrix_t        *matrix,
                            bool                symmetric,
                            cs_lnum_t           diag_block_size,
                            cs_lnum_t           extra_diag_block_size,
                            const cs_lnum_t     n_edges,
                            const cs_lnum_2_t   edges[],
                            const cs_real_t    *da,
                            const cs_real_t    *xa)
{
  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0, _("The matrix is not defined."));

  cs_base_check_bool(&symmetric);

  _set_fill_info(matrix,
                 symmetric,
                 diag_block_size,
                 extra_diag_block_size);

  if (matrix->set_coefficients != NULL)
    matrix->set_coefficients(matrix, symmetric, true, n_edges, edges, da, xa);
  else
    bft_error
      (__FILE__, __LINE__, 0,
       "Matrix format %s with fill type %s does not handle\n"
       "coefficient assignment from native (graph-edge) coefficients.",
       matrix->type_name,
       cs_matrix_fill_type_name[matrix->fill_type]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set matrix coefficients in an MSR format, transferring the
 * property of those arrays to the matrix.
 *
 * If the matrix is also in MSR format, this avoids an extra copy.
 * If it is in a different format, values are copied to the structure,
 * and the original arrays freed. In any case, the arrays pointers passed as
 * arguments are set to NULL, to help ensure the caller does not use the
 * original arrays directly after this call.
 *
 * \param[in, out]  matrix                 pointer to matrix structure
 * \param[in]       symmetric              indicates if matrix coefficients
 *                                         are symmetric
 * \param[in]       diag_block_size        block sizes for diagonal
 * \param[in]       extra_diag_block_size  block sizes for extra diagonal
 * \param[in]       row_index              MSR row index (0 to n-1)
 * \param[in]       col_id                 MSR column id (0 to n-1)
 * \param[in, out]  d_val                  diagonal values (NULL if zero)
 * \param[in, out]  x_val                  extradiagonal values (NULL if zero)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_transfer_coefficients_msr(cs_matrix_t         *matrix,
                                    bool                 symmetric,
                                    cs_lnum_t            diag_block_size,
                                    cs_lnum_t            extra_diag_block_size,
                                    const cs_lnum_t      row_index[],
                                    const cs_lnum_t      col_id[],
                                    cs_real_t          **d_val,
                                    cs_real_t          **x_val)
{
  const cs_real_t  *d_val_p = (d_val != NULL) ? *d_val : NULL;
  const cs_real_t  *x_val_p = (x_val != NULL) ? *x_val : NULL;

  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0, _("The matrix is not defined."));

  cs_base_check_bool(&symmetric);

  _set_fill_info(matrix,
                 symmetric,
                 diag_block_size,
                 extra_diag_block_size);

  switch(matrix->type) {

  case CS_MATRIX_CSR:
    _set_coeffs_csr_from_msr(matrix,
                             row_index,
                             col_id,
                             d_val_p,
                             d_val,
                             x_val_p,
                             x_val);
    break;

  case CS_MATRIX_MSR:
    _set_coeffs_msr_from_msr(matrix,
                             false, /* ignored in case of transfer */
                             row_index,
                             col_id,
                             d_val_p,
                             d_val,
                             x_val_p,
                             x_val);
    break;

  default:
    bft_error
      (__FILE__, __LINE__, 0,
       "Matrix format %s with fill type %s does not handle\n"
       "coefficient assignment from native (graph-edge) coefficients.",
       matrix->type_name,
       cs_matrix_fill_type_name[matrix->fill_type]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Release shared matrix coefficients.
 *
 * Pointers to mapped coefficients are set to NULL, while
 * coefficient copies owned by the matrix are not modified.
 *
 * This simply ensures the matrix does not maintain pointers
 * to nonexistant data.
 *
 * \param[in, out]  matrix  pointer to matrix structure
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_release_coefficients(cs_matrix_t  *matrix)
{
  /* Check API state */

  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0, _("The matrix is not defined."));

  if (matrix->destroy_adaptor != NULL) {
    matrix->destroy_adaptor(matrix);
  }

  if (matrix->release_coefficients != NULL) {
    matrix->xa = NULL;
    matrix->release_coefficients(matrix);
  }
  else {
    bft_error
      (__FILE__, __LINE__, 0,
       "Matrix format %s is missing a release_coefficients function.",
       matrix->type_name);
  }

  /* Set fill type to impossible value */

  _clear_fill_info(matrix);
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

cs_matrix_assembler_values_t *
cs_matrix_assembler_values_init(cs_matrix_t  *matrix,
                                cs_lnum_t     diag_block_size,
                                cs_lnum_t     extra_diag_block_size)
{
  cs_matrix_assembler_values_t *mav = NULL;

  /* Set fill type */

  _set_fill_info(matrix,
                 false, /* symmetric */
                 diag_block_size,
                 extra_diag_block_size);

  /* Create values assembler */

  if (matrix->assembler_values_create != NULL)
    mav = matrix->assembler_values_create(matrix,
                                          diag_block_size,
                                          extra_diag_block_size);

  else
    bft_error(__FILE__, __LINE__, 0,
              _("%s: direct assembly handling of matrices of type %s\n"
                "is not available."),
              __func__, _(matrix->type_name));

  return mav;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy matrix diagonal values.
 *
 * In case of matrixes with block diagonal coefficients, only the true
 * diagonal values are copied.
 *
 * \param[in]   matrix  pointer to matrix structure
 * \param[out]  da      diagonal (pre-allocated, size: n_rows*block_size)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_copy_diagonal(const cs_matrix_t  *matrix,
                        cs_real_t          *restrict da)
{
  /* Check API state */

  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0, _("The matrix is not defined."));

  if (matrix->copy_diagonal != NULL)
    matrix->copy_diagonal(matrix, da);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Query matrix coefficients symmetry
 *
 * \param[in]  matrix  pointer to matrix structure
 *
 * \return  true if coefficients are symmetric, false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_matrix_is_symmetric(const cs_matrix_t  *matrix)
{
  return matrix->symmetric;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Indicate whether coefficients were mapped from native face-based
 *        arrays.
 *
 * It is used in the current multgrid code, but should be removed as soon
 * as the dependency to the native format is removed.
 *
 * \param[in]  matrix  pointer to matrix structure
 *
 * \return  true if coefficients were mapped from native face-based arrays,
 *          false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_matrix_is_mapped_from_native(const cs_matrix_t  *matrix)
{
  bool retval = false;

  if (matrix->xa != NULL)
    retval = true;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get matrix diagonal values.
 *
 * In case of matrixes with block diagonal coefficients, a pointer to
 * the complete block diagonal is returned.
 *
 * \param[in]  matrix  pointer to matrix structure
 *
 * \return  pointer to matrix diagonal array
 */
/*----------------------------------------------------------------------------*/

const cs_real_t *
cs_matrix_get_diagonal(const cs_matrix_t  *matrix)
{
  const cs_real_t  *diag = NULL;

  if (matrix->get_diagonal != NULL)
    diag = matrix->get_diagonal(matrix);

  else
    bft_error(__FILE__, __LINE__, 0,
              _("%s: not available for matrix type: %s."),
              __func__, cs_matrix_get_type_name(matrix));

  return diag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get pointer to matrix extra-diagonal values in "native" format
 *
 * \deprecated
 *
 * This function only functions if the coefficients were mapped from native
 * coefficients using cs_matrix_set_coefficients(), in which case the pointer
 * returned is the same as the one passed to that function.
 *
 * It is used in the current multigrid code, but should be removed as soon
 * as the dependency to the native format is removed.
 *
 * \param[in]  matrix  pointer to matrix structure
 *
 * \return  pointer to matrix diagonal array
 */
/*----------------------------------------------------------------------------*/

const cs_real_t *
cs_matrix_get_extra_diagonal(const cs_matrix_t  *matrix)
{
  const cs_real_t  *exdiag = NULL;

  if (matrix->xa == NULL)
    bft_error
      (__FILE__, __LINE__, 0,
       _("Matrix coefficients were not mapped from native face-based arrays,\n"
         "so the extra-diagonal coefficients are not available in that form."));
  else
    exdiag = matrix->xa;

  return exdiag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize row info for a given matrix.
 *
 * \param[out]  r   row info structure
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_row_init(cs_matrix_row_info_t  *r)
{
  r->row_size = 0;
  r->buffer_size = 0;
  r->col_id = NULL;
  r->_col_id = NULL;
  r->vals = NULL;
  r->_vals = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize row info for a given matrix.
 *
 * \param[in, out]  r   row info structure
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_row_finalize(cs_matrix_row_info_t  *r)
{
  r->row_size = 0;
  r->buffer_size = 0;
  r->col_id = NULL;
  BFT_FREE(r->_col_id);
  r->vals = NULL;
  BFT_FREE(r->_vals);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get row values for a given matrix.
 *
 * This function may not work for all matrix types.
 *
 * In the case of blocked matrices, the true (non-blocked)
 * values are returned.
 *
 * The row information structure must have been previously initialized
 * using \ref cs_matrix_row_init, and should be finalized using
 * using \ref cs_matrix_row_finalize, so as to free buffers it may have
 * built for certain matrix formats.
 *
 * \param[in]       matrix     pointer to matrix structure
 * \param[in]       row_id     id of row to query
 * \param[in, out]  r          row info structure
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_get_row(const cs_matrix_t     *matrix,
                  const cs_lnum_t        row_id,
                  cs_matrix_row_info_t  *r)
{
  cs_lnum_t b_size = matrix->db_size;

  switch (matrix->type) {

  case CS_MATRIX_CSR:
    {
      const cs_matrix_struct_csr_t  *ms = matrix->structure;
      const cs_matrix_coeff_csr_t  *mc = matrix->coeffs;
      r->row_size = (ms->row_index[row_id+1] - ms->row_index[row_id])*b_size;
      r->col_id = ms->col_id + ms->row_index[row_id]*b_size;
      if (mc->val != NULL)
        r->vals = mc->val + ms->row_index[row_id]*b_size;
      else
        r->vals = NULL;
    }
    break;

  case CS_MATRIX_MSR:
    {
      const cs_lnum_t _row_id = row_id / b_size;
      const cs_matrix_struct_csr_t  *ms = matrix->structure;
      const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;
      const cs_lnum_t n_ed_cols =   ms->row_index[_row_id+1]
                                  - ms->row_index[_row_id];
      if (matrix->eb_size == 1)
        r->row_size = n_ed_cols + b_size;
      else
        r->row_size = (n_ed_cols+1)*b_size;
      if (r->buffer_size < r->row_size) {
        r->buffer_size = r->row_size*2;
        BFT_REALLOC(r->_col_id, r->buffer_size, cs_lnum_t);
        r->col_id = r->_col_id;
        BFT_REALLOC(r->_vals, r->buffer_size, cs_real_t);
        r->vals = r->_vals;
      }
      cs_lnum_t ii = 0, jj = 0;
      const cs_lnum_t *restrict c_id = ms->col_id + ms->row_index[_row_id];
      if (b_size == 1) {
        const cs_real_t *m_row = mc->e_val + ms->row_index[_row_id];
        for (jj = 0; jj < n_ed_cols && c_id[jj] < _row_id; jj++) {
          r->_col_id[ii] = c_id[jj];
          r->_vals[ii++] = m_row[jj];
        }
        r->_col_id[ii] = _row_id;
        r->_vals[ii++] = mc->d_val[_row_id];
        for (; jj < n_ed_cols; jj++) {
          r->_col_id[ii] = c_id[jj];
          r->_vals[ii++] = m_row[jj];
        }
      }
      else if (matrix->eb_size == 1) {
        const cs_lnum_t _sub_id = row_id % b_size;
        const cs_lnum_t db_size = matrix->db_size;
        const cs_lnum_t db_size_2 = matrix->db_size*matrix->db_size;
        const cs_real_t *m_row = mc->e_val + ms->row_index[_row_id];
        for (jj = 0; jj < n_ed_cols && c_id[jj] < _row_id; jj++) {
          r->_col_id[ii] = c_id[jj]*b_size + _sub_id;
          r->_vals[ii++] = m_row[jj];
        }
        for (cs_lnum_t kk = 0; kk < b_size; kk++) {
          r->_col_id[ii] = _row_id*b_size + kk;
          r->_vals[ii++] = mc->d_val[  _row_id*db_size_2
                                     + _sub_id*db_size + kk];
        }
        for (; jj < n_ed_cols; jj++) {
          r->_col_id[ii] = c_id[jj]*b_size + _sub_id;
          r->_vals[ii++] = m_row[jj];
        }
      }
      else {
        const cs_lnum_t _sub_id = row_id % b_size;
        const cs_lnum_t db_size = matrix->db_size;
        const cs_lnum_t eb_size = matrix->db_size;
        const cs_lnum_t db_size_2 = matrix->db_size*matrix->db_size;
        const cs_lnum_t eb_size_2 = matrix->eb_size*matrix->eb_size;
        const cs_real_t *m_row = mc->e_val + ms->row_index[_row_id]*eb_size_2;
        for (jj = 0; jj < n_ed_cols && c_id[jj] < _row_id; jj++) {
          for (cs_lnum_t kk = 0; kk < b_size; kk++) {
            r->_col_id[ii] = c_id[jj]*b_size + kk;
            r->_vals[ii++] = m_row[_sub_id*eb_size + kk];
          }
        }
        for (cs_lnum_t kk = 0; kk < b_size; kk++) {
          r->_col_id[ii] = _row_id*b_size + kk;
          r->_vals[ii++] = mc->d_val[  _row_id*db_size_2
                                     + _sub_id*db_size + kk];
        }
        for (; jj < n_ed_cols; jj++) {
          for (cs_lnum_t kk = 0; kk < b_size; kk++) {
            r->_col_id[ii] = c_id[jj]*b_size + kk;
            r->_vals[ii++] = m_row[_sub_id*eb_size + kk];
          }
        }
      }
    }
    break;

  case CS_MATRIX_DIST:
    {
      const cs_lnum_t _row_id = row_id / b_size;
      const cs_matrix_struct_dist_t  *ms = matrix->structure;
      const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;
      const cs_lnum_t n_ed_cols =   ms->e.row_index[_row_id+1]
                                  - ms->e.row_index[_row_id];
      const cs_lnum_t n_h_cols =    ms->h.row_index[_row_id+1]
                                  - ms->h.row_index[_row_id];
      if (matrix->eb_size == 1)
        r->row_size = n_ed_cols + n_h_cols + b_size;
      else
        r->row_size = (n_ed_cols + n_h_cols + 1)*b_size;
      if (r->buffer_size < r->row_size) {
        r->buffer_size = r->row_size*2;
        BFT_REALLOC(r->_col_id, r->buffer_size, cs_lnum_t);
        r->col_id = r->_col_id;
        BFT_REALLOC(r->_vals, r->buffer_size, cs_real_t);
        r->vals = r->_vals;
      }
      cs_lnum_t ii = 0, jj = 0;

      /* Local elements (including diagonal) */

      const cs_lnum_t *restrict c_id = ms->e.col_id + ms->e.row_index[_row_id];
      if (b_size == 1) {
        const cs_real_t *m_row = mc->e_val + ms->e.row_index[_row_id];
        for (jj = 0; jj < n_ed_cols && c_id[jj] < _row_id; jj++) {
          r->_col_id[ii] = c_id[jj];
          r->_vals[ii++] = m_row[jj];
        }
        r->_col_id[ii] = _row_id;
        r->_vals[ii++] = mc->d_val[_row_id];
        for (; jj < n_ed_cols; jj++) {
          r->_col_id[ii] = c_id[jj];
          r->_vals[ii++] = m_row[jj];
        }
      }
      else if (matrix->eb_size == 1) {
        const cs_lnum_t _sub_id = row_id % b_size;
        const cs_lnum_t db_size = matrix->db_size;
        const cs_lnum_t db_size_2 = matrix->db_size*matrix->db_size;
        const cs_real_t *m_row = mc->e_val + ms->e.row_index[_row_id];
        for (jj = 0; jj < n_ed_cols && c_id[jj] < _row_id; jj++) {
          r->_col_id[ii] = c_id[jj]*b_size + _sub_id;
          r->_vals[ii++] = m_row[jj];
        }
        for (cs_lnum_t kk = 0; kk < b_size; kk++) {
          r->_col_id[ii] = _row_id*b_size + kk;
          r->_vals[ii++] = mc->d_val[  _row_id*db_size_2
                                     + _sub_id*db_size + kk];
        }
        for (; jj < n_ed_cols; jj++) {
          r->_col_id[ii] = c_id[jj]*b_size + _sub_id;
          r->_vals[ii++] = m_row[jj];
        }
      }
      else {
        const cs_lnum_t _sub_id = row_id % b_size;
        const cs_lnum_t db_size = matrix->db_size;
        const cs_lnum_t eb_size = matrix->db_size;
        const cs_lnum_t db_size_2 = matrix->db_size*matrix->db_size;
        const cs_lnum_t eb_size_2 = matrix->eb_size*matrix->eb_size;
        const cs_real_t *m_row = mc->e_val + ms->e.row_index[_row_id]*eb_size_2;
        for (jj = 0; jj < n_ed_cols && c_id[jj] < _row_id; jj++) {
          for (cs_lnum_t kk = 0; kk < b_size; kk++) {
            r->_col_id[ii] = c_id[jj]*b_size + kk;
            r->_vals[ii++] = m_row[_sub_id*eb_size + kk];
          }
        }
        for (cs_lnum_t kk = 0; kk < b_size; kk++) {
          r->_col_id[ii] = _row_id*b_size + kk;
          r->_vals[ii++] = mc->d_val[  _row_id*db_size_2
                                     + _sub_id*db_size + kk];
        }
        for (; jj < n_ed_cols; jj++) {
          for (cs_lnum_t kk = 0; kk < b_size; kk++) {
            r->_col_id[ii] = c_id[jj]*b_size + kk;
            r->_vals[ii++] = m_row[_sub_id*eb_size + kk];
          }
        }
      }

      /* Halo elements */

      c_id = ms->h.col_id + ms->h.row_index[_row_id];
      if (b_size == 1) {
        const cs_real_t *m_row = mc->h_val + ms->h.row_index[_row_id];
        for (jj = 0; jj < n_h_cols; jj++) {
          r->_col_id[ii] = c_id[jj];
          r->_vals[ii++] = m_row[jj];
        }
      }
      else if (matrix->eb_size == 1) {
        const cs_lnum_t _sub_id = row_id % b_size;
        const cs_real_t *m_row = mc->h_val + ms->h.row_index[_row_id];
        for (jj = 0; jj < n_h_cols; jj++) {
          r->_col_id[ii] = c_id[jj]*b_size + _sub_id;
          r->_vals[ii++] = m_row[jj];
        }
      }
      else {
        const cs_lnum_t _sub_id = row_id % b_size;
        const cs_lnum_t b_size_2 = b_size * b_size;
        const cs_real_t *m_row = mc->h_val + ms->h.row_index[_row_id]*b_size_2;
        for (jj = 0; jj < n_h_cols; jj++) {
          for (cs_lnum_t kk = 0; kk < b_size; kk++) {
            r->_col_id[ii] = c_id[jj]*b_size + kk;
            r->_vals[ii++] = m_row[_sub_id*b_size + kk];
          }
        }
      }
    }

    break;

  default:
    bft_error
      (__FILE__, __LINE__, 0,
       _("Matrix format %s with fill type %s does not handle %s operation."),
       matrix->type_name,
       cs_matrix_fill_type_name[matrix->fill_type],
       __func__);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get arrays describing a matrix in native format.
 *
 * This function works for matrix in native format.
 *
 * Matrix block sizes can be obtained by cs_matrix_get_diag_block_size()
 * and cs_matrix_get_extra_diag_block_size().
 *
 * \param[in]   matrix     pointer to matrix structure
 * \param[out]  symmetric  true if symmetric
 * \param[out]  n_edges    number of associated faces
 * \param[out]  edges      edges (symmetric row <-> column) connectivity
 * \param[out]  d_val      diagonal values
 * \param[out]  x_val      extra-diagonal values
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_get_native_arrays(const cs_matrix_t   *matrix,
                            bool                *symmetric,
                            cs_lnum_t           *n_edges,
                            const cs_lnum_2_t  **edges,
                            const cs_real_t    **d_val,
                            const cs_real_t    **x_val)
{
  if (symmetric != NULL)
    *symmetric = false;
  if (n_edges != NULL)
    *n_edges = 0;
  if (edges != NULL)
    *edges = NULL;
  if (d_val != NULL)
    *d_val = NULL;
  if (x_val != NULL)
    *x_val = NULL;

  if (matrix->type == CS_MATRIX_NATIVE) {
    const cs_matrix_struct_native_t  *ms = matrix->structure;
    const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;
    if (n_edges != NULL)
      *n_edges = ms->n_edges;
    if (edges != NULL)
      *edges = ms->edges;
    if (mc != NULL) {
      if (symmetric != NULL)
        *symmetric = mc->symmetric;
      if (d_val != NULL)
        *d_val = mc->d_val;
      if (x_val != NULL)
        *x_val = mc->e_val;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get arrays describing a matrix in CSR format.
 *
 * This function only works for an CSR matrix (i.e. there is
 * no automatic conversion from another matrix type).
 *
 * Matrix block sizes can be obtained by cs_matrix_get_diag_block_size()
 * and cs_matrix_get_extra_diag_block_size().
 *
 * \param[in]   matrix     pointer to matrix structure
 * \param[out]  row_index  CSR row index
 * \param[out]  col_id     CSR column id
 * \param[out]  val        values
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_get_csr_arrays(const cs_matrix_t   *matrix,
                         const cs_lnum_t    **row_index,
                         const cs_lnum_t    **col_id,
                         const cs_real_t    **val)
{
  if (row_index != NULL)
    *row_index = NULL;
  if (col_id != NULL)
    *col_id = NULL;
  if (val != NULL)
    *val = NULL;

  if (matrix->type == CS_MATRIX_CSR) {
    const cs_matrix_struct_csr_t  *ms = matrix->structure;
    const cs_matrix_coeff_csr_t  *mc = matrix->coeffs;
    if (row_index != NULL)
      *row_index = ms->row_index;
    if (col_id != NULL)
      *col_id = ms->col_id;
    if (val != NULL && mc != NULL) {
      *val = mc->val;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get arrays describing a matrix in MSR format.
 *
 * This function only works for an MSR matrix (i.e. there is
 * no automatic conversion from another matrix type).
 *
 * Matrix block sizes can be obtained by cs_matrix_get_diag_block_size()
 * and cs_matrix_get_extra_diag_block_size().
 *
 * \param[in]   matrix     pointer to matrix structure
 * \param[out]  row_index  MSR row index
 * \param[out]  col_id     MSR column id
 * \param[out]  d_val      diagonal values
 * \param[out]  x_val      extra-diagonal values
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_get_msr_arrays(const cs_matrix_t   *matrix,
                         const cs_lnum_t    **row_index,
                         const cs_lnum_t    **col_id,
                         const cs_real_t    **d_val,
                         const cs_real_t    **x_val)
{
  if (row_index != NULL)
    *row_index = NULL;
  if (col_id != NULL)
    *col_id = NULL;
  if (d_val != NULL)
    *d_val = NULL;
  if (x_val != NULL)
    *x_val = NULL;

  if (matrix->type == CS_MATRIX_MSR) {
    const cs_matrix_struct_dist_t  *ms = matrix->structure;
    const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;
    if (row_index != NULL)
      *row_index = ms->e.row_index;
    if (col_id != NULL)
      *col_id = ms->e.col_id;
    if (mc != NULL) {
      if (d_val != NULL)
        *d_val = mc->d_val;
      if (x_val != NULL)
        *x_val = mc->e_val;
    }
  }
  else
    bft_error
      (__FILE__, __LINE__, 0,
       _("%s is not available for matrix using %s storage."),
       __func__, matrix->type_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Associate mesh information with a matrix.
 *
 * This may be useful for multigrid smoothing.
 *
 * At least cell centers and volumes are needed for relaxation, and face
 * adjacency and normals are needed for the "classical" option.
 *
 * Note that cells and faces here do not need to be primary mesh elements,
 * but could be dual mesh elements of some sort.
 *
 * The arrays passed to the matrix are shared, so should have a lifetime
 * at least as long as the matrix.
 *
 * \param[in, out]   matrix       pointer to matrix structure
 * \param[in]        c2f_idx      cell to faces index, or NULL
 * \param[in]        c2f          cell to faces adjacency, or NULL
 * \param[in]        c2f_sgn      cell to faces adjacency sign, or NULL
 * \param[in]        cell_cen     cell center coordinates
 * \param[in]        cell_vol     cell volumes
 * \param[in]        face_normal  face normal, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_set_mesh_association(cs_matrix_t         *matrix,
                               const cs_lnum_t     *c2f_idx,
                               const cs_lnum_t     *c2f,
                               const short int     *c2f_sgn,
                               const cs_real_3_t   *cell_cen,
                               const cs_real_t     *cell_vol,
                               const cs_real_3_t   *face_normal)
{
  matrix->c2f_idx = c2f_idx;
  matrix->c2f = c2f;
  matrix->c2f_sgn = c2f_sgn;

  matrix->cell_cen = cell_cen;
  matrix->cell_vol = cell_vol;
  matrix->face_normal = face_normal;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Query mesh information that me be associated with a matrix.
 *
 * This may be useful for multigrid smoothing.
 *
 * \param[in]   matrix       pointer to matrix structure
 * \param[out]  c2f_idx      cell to faces index, or NULL
 * \param[out]  c2f          cell to faces adjacency, or NULL
 * \param[out]  c2f_sgn      cell to faces adjacency sign, or NULL
 * \param[out]  cell_cen     cell center coordinates, or NULL
 * \param[out]  cell_vol     cell volumes, or NULL
 * \param[out]  face_normal  face normas, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_get_mesh_association(const cs_matrix_t   *matrix,
                               const cs_lnum_t    **c2f_idx,
                               const cs_lnum_t    **c2f,
                               const short int    **c2f_sgn,
                               const cs_real_3_t  **cell_cen,
                               const cs_real_t    **cell_vol,
                               const cs_real_3_t  **face_normal)
{
  if (c2f_idx != NULL)
    *c2f_idx = matrix->c2f_idx;
  if (c2f != NULL)
    *c2f = matrix->c2f;
  if (c2f_sgn != NULL)
    *c2f_sgn = matrix->c2f_sgn;

  if (cell_cen != NULL)
    *cell_cen = matrix->cell_cen;
  if (cell_vol != NULL)
    *cell_vol = matrix->cell_vol;
  if (face_normal != NULL)
    *face_normal = matrix->face_normal;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x
 *
 * This function includes a halo update of x prior to multiplication by A.
 *
 * \param[in]       matrix         pointer to matrix structure
 * \param[in, out]  x              multiplying vector values
 *                                 (ghost values updated)
 * \param[out]      y              resulting vector
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_vector_multiply(const cs_matrix_t   *matrix,
                          cs_real_t           *restrict x,
                          cs_real_t           *restrict y)
{
  assert(matrix != NULL);

  if (matrix->vector_multiply[matrix->fill_type][0] != NULL) {

#if defined(HAVE_ACCEL)
    if (   matrix->vector_multiply[matrix->fill_type][0]
        == matrix->vector_multiply_d[matrix->fill_type][0]) {
      cs_alloc_mode_t md_x = cs_check_device_ptr(x);
      cs_alloc_mode_t md_y = cs_check_device_ptr(y);
      if (md_x == CS_ALLOC_HOST || md_y == CS_ALLOC_HOST)
        matrix->vector_multiply_h[matrix->fill_type][0](matrix, false, true,
                                                        x, y);
      else {
        cs_real_t *d_x = (cs_real_t *)cs_get_device_ptr(x);
        cs_real_t *d_y = (cs_real_t *)cs_get_device_ptr(y);

        matrix->vector_multiply[matrix->fill_type][0](matrix, false, true,
                                                      d_x, d_y);
      }
    }
    else
      matrix->vector_multiply[matrix->fill_type][0](matrix, false, true, x, y);
#else
    matrix->vector_multiply[matrix->fill_type][0](matrix, false, true, x, y);
#endif

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              _("%s: Matrix of type: %s is missing a vector multiply\n"
                "function for fill type %s."),
              __func__, cs_matrix_get_type_name(matrix),
              cs_matrix_fill_type_name[matrix->fill_type]);
}

#if defined(HAVE_ACCEL)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x, on device
 *
 * This function includes a halo update of x prior to multiplication by A.
 *
 * \param[in]       matrix         pointer to matrix structure
 * \param[in, out]  x              multiplying vector values, on device
 *                                 (ghost values updated)
 * \param[out]      y              resulting vector (on device)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_vector_multiply_d(const cs_matrix_t   *matrix,
                            cs_real_t           *restrict x,
                            cs_real_t           *restrict y)
{
  assert(matrix != NULL);

  if (matrix->vector_multiply_d[matrix->fill_type][0] != NULL)
    matrix->vector_multiply_d[matrix->fill_type][0](matrix, false, true, x, y);

  else
    bft_error(__FILE__, __LINE__, 0,
              _("%s: Matrix of type: %s is missing a device vector multiply\n"
                "function for fill type %s."),
              __func__, cs_matrix_get_type_name(matrix),
              cs_matrix_fill_type_name[matrix->fill_type]);
}

#endif /* defined(HAVE_ACCEL) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x with no prior halo update of x.
 *
 * This function does not include a halo update of x prior to multiplication
 * by A, so it should be called only when the halo of x is known to already
 * be up to date (in which case we avoid the performance penalty of a
 * redundant update by using this variant of the matrix.vector product).
 *
 * \param[in]   matrix         pointer to matrix structure
 * \param[in]   x              multiplying vector values
 * \param[out]  y              resulting vector
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_vector_multiply_nosync(const cs_matrix_t  *matrix,
                                 cs_real_t          *restrict x,
                                 cs_real_t          *restrict y)
{
  assert(matrix != NULL);

  if (matrix->vector_multiply[matrix->fill_type][0] != NULL) {

#if defined(HAVE_ACCEL)
    if (   matrix->vector_multiply[matrix->fill_type][0]
        == matrix->vector_multiply_d[matrix->fill_type][0]) {
      cs_alloc_mode_t md_x = cs_check_device_ptr(x);
      cs_alloc_mode_t md_y = cs_check_device_ptr(y);
      if (md_x == CS_ALLOC_HOST || md_y == CS_ALLOC_HOST)
        matrix->vector_multiply_h[matrix->fill_type][0](matrix, false, false,
                                                        x, y);
      else {
        cs_real_t *d_x = (cs_real_t *)cs_get_device_ptr(x);
        cs_real_t *d_y = (cs_real_t *)cs_get_device_ptr(y);

        matrix->vector_multiply[matrix->fill_type][0](matrix, false, false,
                                                      d_x, d_y);
      }
    }
    else
      matrix->vector_multiply[matrix->fill_type][0](matrix, false, false, x, y);
#else
    matrix->vector_multiply[matrix->fill_type][0](matrix, false, false, x, y);
#endif

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              _("%s: Matrix of type: %s is missing a vector multiply\n"
                "function for fill type %s."),
              __func__, cs_matrix_get_type_name(matrix),
              cs_matrix_fill_type_name[matrix->fill_type]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Partial matrix.vector product.
 *
 * This function includes a halo update of x prior to multiplication,
 * except for the CS_MATRIX_SPMV_L operation type, which does not require it,
 * as halo adjacencies are only present and useful in the upper-diagonal part..
 *
 * \param[in]       matrix         pointer to matrix structure
 * \param[in]       op_type        SpMV operation type
 * \param[in, out]  x              multiplying vector values
 *                                 (ghost values updated)
 * \param[out]      y              resulting vector
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_vector_multiply_partial(const cs_matrix_t      *matrix,
                                  cs_matrix_spmv_type_t   op_type,
                                  cs_real_t              *restrict x,
                                  cs_real_t              *restrict y)
{
  assert(matrix != NULL);

  if (matrix->vector_multiply[matrix->fill_type][op_type] != NULL) {

#if defined(HAVE_ACCEL)
    if (   matrix->vector_multiply[matrix->fill_type][op_type]
        == matrix->vector_multiply_d[matrix->fill_type][op_type]) {
      cs_alloc_mode_t md_x = cs_check_device_ptr(x);
      cs_alloc_mode_t md_y = cs_check_device_ptr(y);
      if (md_x == CS_ALLOC_HOST || md_y == CS_ALLOC_HOST)
        matrix->vector_multiply_h[matrix->fill_type][op_type]
                  (matrix, true, true, x, y);
      else {
        cs_real_t *d_x = (cs_real_t *)cs_get_device_ptr(x);
        cs_real_t *d_y = (cs_real_t *)cs_get_device_ptr(y);

        matrix->vector_multiply[matrix->fill_type][op_type]
                  (matrix, true, true, d_x, d_y);
      }
    }
    else
      matrix->vector_multiply[matrix->fill_type][op_type]
                (matrix, true, true, x, y);
#else
    matrix->vector_multiply[matrix->fill_type][op_type]
              (matrix, true, true, x, y);
#endif
  }
  else
    bft_error(__FILE__, __LINE__, 0,
              _("%s: Matrix of type: %s is missing a partial SpMV\n"
                "(%s) function for fill type %s."),
              __func__, cs_matrix_get_type_name(matrix),
              cs_matrix_spmv_type_name[op_type],
              cs_matrix_fill_type_name[matrix->fill_type]);
}

#if defined(HAVE_ACCEL)

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Partial matrix.vector product, on device
 *
 * This function includes a halo update of x prior to multiplication,
 * except for the CS_MATRIX_SPMV_L operation type, which does not require it,
 * as halo adjacencies are only present and useful in the upper-diagonal part..
 *
 * \param[in]       matrix         pointer to matrix structure
 * \param[in]       op_type        SpMV operation type
 * \param[in, out]  x              multiplying vector values, on device
 *                                 (ghost values updated)
 * \param[out]      y              resulting vector, on device
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_vector_multiply_partial_d(const cs_matrix_t      *matrix,
                                    cs_matrix_spmv_type_t   op_type,
                                    cs_real_t              *restrict x,
                                    cs_real_t              *restrict y)
{
  assert(matrix != NULL);

  if (matrix->vector_multiply_d[matrix->fill_type][op_type] != NULL)
    matrix->vector_multiply_d[matrix->fill_type][op_type]
              (matrix, true, true, x, y);

  else
    bft_error(__FILE__, __LINE__, 0,
              _("%s: Matrix of type: %s is missing a device partial SpMV\n"
                "(%s) function for fill type %s."),
              __func__, cs_matrix_get_type_name(matrix),
              cs_matrix_spmv_type_name[op_type],
              cs_matrix_fill_type_name[matrix->fill_type]);
}

#endif /* defined(HAVE_ACCEL) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build matrix variant
 *
 * The variant will initially use default matrix-vector functions,
 * which can be later modified using cs_matrix_variant_set_func().
 *
 * \param[in]  m   pointer to matrix
 */
/*----------------------------------------------------------------------------*/

cs_matrix_variant_t *
cs_matrix_variant_create(cs_matrix_t  *m)
{
  cs_matrix_variant_t  *mv;

  BFT_MALLOC(mv, 1, cs_matrix_variant_t);

  mv->type = m->type;
  mv->fill_type = m->fill_type;

  for (int j = 0; j < 2; j++) {
    mv->vector_multiply[j] = NULL;
    strncpy(mv->name[j], "default", 31);
    mv->name[j][31] = '\0';
  }

  for (cs_matrix_fill_type_t mft = 0; mft < CS_MATRIX_N_FILL_TYPES; mft++) {
    for (cs_matrix_spmv_type_t spmv_type = 0;
         spmv_type < CS_MATRIX_SPMV_N_TYPES;
         spmv_type++)
      cs_matrix_spmv_set_func(m->type,
                              m->fill_type,
                              spmv_type,
                              m->numbering,
                              NULL, /* func_name */
                              mv->vector_multiply,
                              mv->vector_multiply_xy_hd);
  }

  return mv;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build list of variants for tuning or testing.
 *
 * The matrix coefficients should be assigned, so the fill type can
 * be determined.
 *
 * \param[in]   m             associated matrix
 * \param[out]  n_variants    number of variants
 * \param[out]  m_variant     array of matrix variants
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_variant_build_list(const cs_matrix_t       *m,
                             int                     *n_variants,
                             cs_matrix_variant_t    **m_variant)
{
  int  n_variants_max = 0;

  *n_variants = 0;
  *m_variant = NULL;

  if (m->type == CS_MATRIX_NATIVE) {

    _variant_add(_("native, baseline"),
                 m->type,
                 m->fill_type,
                 m->numbering,
                 "default",
                 n_variants,
                 &n_variants_max,
                 m_variant);

    if (m->numbering != NULL) {

#if defined(HAVE_OPENMP)

      if (m->numbering->type == CS_NUMBERING_THREADS) {
        _variant_add(_("native, OpenMP"),
                     m->type,
                     m->fill_type,
                     m->numbering,
                     "default",
                     n_variants,
                     &n_variants_max,
                     m_variant);
      }

      _variant_add(_("native, OpenMP atomic"),
                   m->type,
                   m->fill_type,
                   m->numbering,
                   "omp_atomic",
                   n_variants,
                   &n_variants_max,
                   m_variant);

#endif

      if (m->numbering->type == CS_NUMBERING_VECTORIZE) {
        _variant_add(_("native, vectorized"),
                     m->type,
                     m->fill_type,
                     m->numbering,
                     "default",
                     n_variants,
                     &n_variants_max,
                     m_variant);
      }

    }

  }

  if (m->type == CS_MATRIX_CSR) {

    _variant_add(_("CSR"),
                 m->type,
                 m->fill_type,
                 m->numbering,
                 "default",
                 n_variants,
                 &n_variants_max,
                 m_variant);

#if defined(HAVE_MKL)

    _variant_add(_("CSR, with MKL"),
                 m->type,
                 m->fill_type,
                 m->numbering,
                 "mkl",
                 n_variants,
                 &n_variants_max,
                 m_variant);

#endif /* defined(HAVE_MKL) */

#if defined(HAVE_CUDA)

    if (cs_get_device_id() > -1)
      _variant_add(_("CSR, CUDA"),
                   m->type,
                   m->fill_type,
                   m->numbering,
                   "cuda",
                   n_variants,
                   &n_variants_max,
                   m_variant);

    if (cs_get_device_id() > -1)
      _variant_add(_("CSR, with cuSPARSE"),
                   m->type,
                   m->fill_type,
                   m->numbering,
                   "cusparse",
                   n_variants,
                   &n_variants_max,
                   m_variant);

#endif /* defined(HAVE_CUDA) */

  }

  if (m->type == CS_MATRIX_MSR) {

    _variant_add(_("MSR"),
                 m->type,
                 m->fill_type,
                 m->numbering,
                 "default",
                 n_variants,
                 &n_variants_max,
                 m_variant);

#if defined(HAVE_MKL)

    _variant_add(_("MSR, with MKL"),
                 m->type,
                 m->fill_type,
                 m->numbering,
                 "mkl",
                 n_variants,
                 &n_variants_max,
                 m_variant);

#endif /* defined(HAVE_MKL) */

#if defined(HAVE_CUDA)

    if (cs_get_device_id() > -1)
      _variant_add(_("MSR, CUDA"),
                   m->type,
                   m->fill_type,
                   m->numbering,
                   "cuda",
                   n_variants,
                   &n_variants_max,
                   m_variant);

#endif /* defined(HAVE_CUDA) */

#if defined(HAVE_CUSPARSE)

    if (cs_get_device_id() > -1)
      _variant_add(_("MSR, with cuSPARSE"),
                   m->type,
                   m->fill_type,
                   m->numbering,
                   "cusparse",
                   n_variants,
                   &n_variants_max,
                   m_variant);

#endif /* defined(HAVE_CUSPARSE) */

#if defined(HAVE_OPENMP)

    if (omp_get_num_threads() > 1) {
      _variant_add(_("MSR, OpenMP scheduling"),
                   m->type,
                   m->fill_type,
                   m->numbering,
                   "omp_sched",
                   n_variants,
                   &n_variants_max,
                   m_variant);
    }

#endif /* defined(HAVE_OPENMP) */

  }

  n_variants_max = *n_variants;
  BFT_REALLOC(*m_variant, *n_variants, cs_matrix_variant_t);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a matrix variant structure.
 *
 * \param[in, out]  mv  pointer to matrix variant pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_variant_destroy(cs_matrix_variant_t  **mv)
{
  if (mv != NULL)
    BFT_FREE(*mv);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply a variant to a given matrix
 *
 * \param[in, out]  m   pointer to matrix
 * \param[in]       mv  pointer to matrix variant pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_variant_apply(cs_matrix_t          *m,
                        cs_matrix_variant_t  *mv)
{
  if (m == NULL || mv == NULL)
    return;

  if (m->destroy_adaptor != NULL)
    m->destroy_adaptor(m);

  if (   m->type < 0 || m->type > CS_MATRIX_N_BUILTIN_TYPES
      || m->fill_type < 0 || m->fill_type > CS_MATRIX_N_FILL_TYPES)
    return;

  for (int i = 0; i < 2; i++) {
    m->vector_multiply[m->fill_type][i] = mv->vector_multiply[i];

#if defined(HAVE_ACCEL)
    if (mv->vector_multiply_xy_hd[i] == 'h')
      m->vector_multiply_h[m->fill_type][i] = mv->vector_multiply[i];
    else if (mv->vector_multiply_xy_hd[i] == 'd')
      m->vector_multiply_d[m->fill_type][i] = mv->vector_multiply[i];
#endif
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply variants defined by tuning to a given matrix
 *
 * \param[in, out]  m   pointer to matrix
 * \param[in]       mv  pointer to matrix variant pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_variant_apply_tuned(cs_matrix_t          *m,
                              cs_matrix_variant_t  *mv)
{
  if (m == NULL || mv == NULL)
    return;

  if (   m->type < 0 || m->type > CS_MATRIX_N_BUILTIN_TYPES
      || m->fill_type < 0 || m->fill_type > CS_MATRIX_N_FILL_TYPES)
    return;

  if (m->destroy_adaptor != NULL)
    m->destroy_adaptor(m);

  for (int i = 0; i < 2; i++)
    m->vector_multiply[m->fill_type][i] = mv->vector_multiply[i];

#if defined(HAVE_ACCEL)
  if (cs_get_device_id() > -1) {
    for (int i = 0; i < 2; i++)
      m->vector_multiply_h[m->fill_type][i] = (mv+1)->vector_multiply[i];
    for (int i = 0; i < 2; i++)
      m->vector_multiply_d[m->fill_type][i] = (mv+2)->vector_multiply[i];
  }
  else {
    for (int i = 0; i < 2; i++)
      m->vector_multiply_h[m->fill_type][i] = mv->vector_multiply[i];
    for (int i = 0; i < 2; i++)
      m->vector_multiply_d[m->fill_type][i] = NULL;
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select the sparse matrix-vector product function to be used by a
 * matrix variant for a given fill type.
 *
 * Currently, possible variant functions are:
 *
 *   CS_MATRIX_NATIVE  (all fill types)
 *     default
 *     baseline
 *     omp             (for OpenMP with compatible numbering)
 *     omp_atomic      (for OpenMP with atomics)
 *     vector          (For vector machine with compatible numbering)
 *
 *   CS_MATRIX_CSR     (for CS_MATRIX_SCALAR or CS_MATRIX_SCALAR_SYM)
 *     default
 *     mkl             (with MKL)
 *
 *   CS_MATRIX_MSR     (all fill types)
 *     default
 *     mkl             (with MKL, for CS_MATRIX_SCALAR or CS_MATRIX_SCALAR_SYM)
 *     omp_sched       (For OpenMP with scheduling)
 *
 * parameters:
 *   mv         <->  pointer to matrix variant
 *   numbering  <--  mesh numbering info, or NULL
 *   fill type  <--  matrix fill type to merge from
 *   spmv_type  <--  SpMV operation type (full or sub-matrix)
 *                   (all types if CS_MATRIX_SPMV_N_TYPES)
 *   func_name  <--  function type name
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_variant_set_func(cs_matrix_variant_t     *mv,
                           cs_matrix_fill_type_t    fill_type,
                           cs_matrix_spmv_type_t    spmv_type,
                           const cs_numbering_t    *numbering,
                           const char              *func_name)
{
  int retcode = cs_matrix_spmv_set_func(mv->type,
                                        fill_type,
                                        spmv_type,
                                        numbering,
                                        func_name,
                                        mv->vector_multiply,
                                        mv->vector_multiply_xy_hd);

  if (retcode == 1) {
    if (spmv_type < CS_MATRIX_SPMV_N_TYPES)
      bft_error
        (__FILE__, __LINE__, 0,
         _("Assignment of matrix.vector product \"%s\" to variant \"%s\"\n"
           "of type \"%s\" for fill \"%s\" not allowed."),
         func_name, mv->name[spmv_type], _matrix_type_name[mv->type],
         cs_matrix_fill_type_name[fill_type]);
    else
      bft_error
        (__FILE__, __LINE__, 0,
         _("Assignment of matrix.vector product \"%s\" to variant \"%s\"\n"
           "of type \"%s\" not allowed."),
         func_name, mv->name[spmv_type], _matrix_type_name[mv->type]);
  }
  else if (retcode == 2)
    bft_error (__FILE__, __LINE__, 0,
               _("Matrix.vector product function type \"%s\"\n"
                 "is not available in this build."), func_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the type associated with a matrix variant.
 *
 * \param[in] mv  pointer to matrix variant structure
 *
 * \return the type of matrix
 */
/*----------------------------------------------------------------------------*/

cs_matrix_type_t
cs_matrix_variant_type(const cs_matrix_variant_t  *mv)
{
  return mv->type;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
