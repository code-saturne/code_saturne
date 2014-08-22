/*============================================================================
 * Sparse Matrix Representation and Operations.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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

/*
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
 *   matrix structure, as this may allow us to leverage existing librairies.
 *
 * - Provide a C interface, also so as to be able to interface more easily
 *   with external libraries.
 *
 * The structures used here could easily be extended to block matrixes,
 * using for example the same structure information with 3x3 blocks which
 * could arise from coupled velocity components. This would imply that the
 * corresponding vectors be interlaced (or an interlaced copy be used
 * for recurring operations such as sparse linear system resolution),
 * for better memory locality, and possible loop unrolling.
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
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_matrix.h"
#include "cs_matrix_priv.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* Variant default for Intel compiler and on Itanium (optimized by BULL)
   (Use compile flag -DNO_BULL_OPTIM to switch default to general code) */

#if (defined(__INTEL_COMPILER) && defined(__ia64__) && !defined(NO_BULL_OPTIM))
#define IA64_OPTIM
#endif

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* Note that most types are declared in cs_matrix_priv.h.
   only those only handled here are declared here. */

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Short names for matrix types */

const char  *cs_matrix_type_name[] = {N_("native"),
                                      N_("CSR"),
                                      N_("symmetric CSR"),
                                      N_("MSR")};

/* Full names for matrix types */

const char
*cs_matrix_type_fullname[] = {N_("diagonal + faces"),
                              N_("Compressed Sparse Row"),
                              N_("symmetric Compressed Sparse Row"),
                              N_("Modified Compressed Sparse Row")};

/* Fill type names for matrices */

const char  *cs_matrix_fill_type_name[] = {"CS_MATRIX_SCALAR",
                                           "CS_MATRIX_SCALAR_SYM",
                                           "CS_MATRIX_33_BLOCK_D",
                                           "CS_MATRIX_33_BLOCK_D_SYM",
                                           "CS_MATRIX_33_BLOCK"};

static char _no_exclude_diag_error_str[]
  = N_("Matrix product variant using function %s\n"
       "does not handle case with excluded diagonal.");

static const char *_matrix_operation_name[CS_MATRIX_N_FILL_TYPES * 2]
  = {N_("y <- A.x"),
     N_("y <- (A-D).x"),
     N_("Symmetric y <- A.x"),
     N_("Symmetric y <- (A-D).x"),
     N_("Block diagonal y <- A.x"),
     N_("Block diagonal y <- (A-D).x"),
     N_("Block diagonal symmetric y <- A.x"),
     N_("Block diagonal symmetric y <- (A-D).x"),
     N_("Block y <- A.x"),
     N_("Block y <- (A-D).x")};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute matrix-vector product for one dense block: y[i] = a[i].x[i]
 *
 * Vectors and blocks may be larger than their useful size, to
 * improve data alignment.
 *
 * parameters:
 *   b_id   <-- block id
 *   b_size <-- block size, including padding:
 *              b_size[0]: useful block size
 *              b_size[1]: vector block extents
 *              b_size[2]: matrix line extents
 *              b_size[3]: matrix line*column (block) extents
 *   a      <-- Pointer to block matrixes array (usually matrix diagonal)
 *   x      <-- Multipliying vector values
 *   y      --> Resulting vector
 *----------------------------------------------------------------------------*/

static inline void
_dense_b_ax(cs_lnum_t         b_id,
            const int         b_size[4],
            const cs_real_t  *restrict a,
            const cs_real_t  *restrict x,
            cs_real_t        *restrict y)
{
  cs_lnum_t   ii, jj;

# if defined(__xlc__) /* Tell IBM compiler not to alias */
# pragma disjoint(*x, *y, * a)
# endif

  for (ii = 0; ii < b_size[0]; ii++) {
    y[b_id*b_size[1] + ii] = 0.;
    for (jj = 0; jj < b_size[0]; jj++)
      y[b_id*b_size[1] + ii]
        +=   a[b_id*b_size[3] + ii*b_size[2] + jj]
           * x[b_id*b_size[1] + jj];
  }
}

/*----------------------------------------------------------------------------
 * Compute matrix-vector product for one dense block: y[i] = a[i].x[i]
 *
 * This variant uses a fixed 3x3 block, for better compiler optimization.
 *
 * parameters:
 *   b_id   <-- block id
 *   a      <-- Pointer to block matrixes array (usually matrix diagonal)
 *   x      <-- Multipliying vector values
 *   y      --> Resulting vector
 *----------------------------------------------------------------------------*/

static inline void
_dense_3_3_ax(cs_lnum_t         b_id,
              const cs_real_t  *restrict a,
              const cs_real_t  *restrict x,
              cs_real_t        *restrict y)
{
# if defined(__xlc__) /* Tell IBM compiler not to alias */
# pragma disjoint(*x, *y, * a)
# endif

  y[b_id*3]     =   a[b_id*9]         * x[b_id*3]
                  + a[b_id*9 + 1]     * x[b_id*3 + 1]
                  + a[b_id*9 + 2]     * x[b_id*3 + 2];

  y[b_id*3 + 1] =   a[b_id*9 + 3]     * x[b_id*3]
                  + a[b_id*9 + 3 + 1] * x[b_id*3 + 1]
                  + a[b_id*9 + 3 + 2] * x[b_id*3 + 2];

  y[b_id*3 + 2] =   a[b_id*9 + 6]     * x[b_id*3]
                  + a[b_id*9 + 6 + 1] * x[b_id*3 + 1]
                  + a[b_id*9 + 6 + 2] * x[b_id*3 + 2];
}

/*----------------------------------------------------------------------------
 * Compute matrix-vector product increment for one dense block:
 * y[i] += a[ij].x[j]
 *
 * Vectors and blocks may be larger than their useful size, to
 * improve data alignment.
 *
 * parameters:
 *   b_i    <-- block id for i
 *   b_j    <-- block id for j
 *   b_ij   <-- block id for matrix ij position
 *   b_size <-- block size, including padding:
 *              b_size[0]: useful block size
 *              b_size[1]: vector block extents
 *              b_size[2]: matrix line extents
 *              b_size[3]: matrix line*column (block) extents
 *   a      <-- Pointer to block matrixes array (usually matrix extra-diagonal)
 *   x      <-- Multipliying vector values
 *   y      --> Resulting vector
 *----------------------------------------------------------------------------*/

static inline void
_dense_eb_ax_add(cs_lnum_t         b_i,
                 cs_lnum_t         b_j,
                 cs_lnum_t         b_ij,
                 const int         b_size[4],
                 const cs_real_t  *restrict a,
                 const cs_real_t  *restrict x,
                 cs_real_t        *restrict y)
{
  cs_lnum_t   ii, jj;

# if defined(__xlc__) /* Tell IBM compiler not to alias */
# pragma disjoint(*x, *y, * a)
# endif

  for (ii = 0; ii < b_size[0]; ii++) {
    for (jj = 0; jj < b_size[0]; jj++)
      y[b_i*b_size[1] + ii]
        +=   a[b_ij*b_size[3] + ii*b_size[2] + jj]
           * x[b_j*b_size[1] + jj];
  }
}

/*----------------------------------------------------------------------------
 * y[i] = da[i].x[i], with da possibly NULL
 *
 * parameters:
 *   da     <-- Pointer to coefficients array (usually matrix diagonal)
 *   x      <-- Multipliying vector values
 *   y      --> Resulting vector
 *   n_elts <-- Array size
 *----------------------------------------------------------------------------*/

static inline void
_diag_vec_p_l(const cs_real_t  *restrict da,
              const cs_real_t  *restrict x,
              cs_real_t        *restrict y,
              cs_lnum_t         n_elts)
{
  cs_lnum_t  ii;

# if defined(__xlc__) /* Tell IBM compiler not to alias */
# pragma disjoint(*x, *y, *da)
# endif

  if (da != NULL) {
#   pragma omp parallel for  if(n_elts > CS_THR_MIN)
    for (ii = 0; ii < n_elts; ii++)
      y[ii] = da[ii] * x[ii];
  }
  else {
#   pragma omp parallel for  if(n_elts > CS_THR_MIN)
    for (ii = 0; ii < n_elts; ii++)
      y[ii] = 0.0;
  }

}

/*----------------------------------------------------------------------------
 * Block version of y[i] = da[i].x[i], with da possibly NULL
 *
 * parameters:
 *   da     <-- Pointer to coefficients array (usually matrix diagonal)
 *   x      <-- Multipliying vector values
 *   y      --> Resulting vector
 *   n_elts <-- Array size
 *   b_size <-- block size, including padding:
 *              b_size[0]: useful block size
 *              b_size[1]: vector block extents
 *              b_size[2]: matrix line extents
 *              b_size[3]: matrix line*column (block) extents
 *----------------------------------------------------------------------------*/

static inline void
_b_diag_vec_p_l(const cs_real_t  *restrict da,
                const cs_real_t  *restrict x,
                cs_real_t        *restrict y,
                cs_lnum_t         n_elts,
                const int         b_size[4])
{
  cs_lnum_t   ii;

  if (da != NULL) {
#   pragma omp parallel for  if(n_elts > CS_THR_MIN)
    for (ii = 0; ii < n_elts; ii++)
      _dense_b_ax(ii, b_size, da, x, y);
  }
  else {
#   pragma omp parallel for  if(n_elts*b_size[1] > CS_THR_MIN)
    for (ii = 0; ii < n_elts*b_size[1]; ii++)
      y[ii] = 0.0;
  }
}

/*----------------------------------------------------------------------------
 * Block version of y[i] = da[i].x[i], with da possibly NULL
 *
 * This variant uses a fixed 3x3 block, for better compiler optimization.
 *
 * parameters:
 *   da     <-- Pointer to coefficients array (usually matrix diagonal)
 *   x      <-- Multipliying vector values
 *   y      --> Resulting vector
 *   n_elts <-- Array size
 *----------------------------------------------------------------------------*/

static inline void
_3_3_diag_vec_p_l(const cs_real_t  *restrict da,
                  const cs_real_t  *restrict x,
                  cs_real_t        *restrict y,
                  cs_lnum_t         n_elts)
{
  cs_lnum_t   ii;

  if (da != NULL) {
#   pragma omp parallel for  if(n_elts*9 > CS_THR_MIN)
    for (ii = 0; ii < n_elts; ii++)
      _dense_3_3_ax(ii, da, x, y);
  }
  else {
#   pragma omp parallel for  if(n_elts*3 > CS_THR_MIN)
    for (ii = 0; ii < n_elts*3; ii++)
      y[ii] = 0.0;
  }
}

/*----------------------------------------------------------------------------
 * Set values from y[start_id] to y[end_id] to 0.
 *
 * parameters:
 *   y        --> Resulting vector
 *   start_id <-- start id in array
 *   end_id   <-- end id in array
 *----------------------------------------------------------------------------*/

static inline void
_zero_range(cs_real_t  *restrict y,
            cs_lnum_t   start_id,
            cs_lnum_t   end_id)
{
  cs_lnum_t   ii;

# pragma omp parallel for  if(end_id - start_id > CS_THR_MIN)
  for (ii = start_id; ii < end_id; ii++)
    y[ii] = 0.0;
}

/*----------------------------------------------------------------------------
 * Set values from y[start_id] to y[end_id] to 0, block version.
 *
 * parameters:
 *   y        --> resulting vector
 *   start_id <-- start id in array
 *   end_id   <-- end id in array
 *   b_size   <-- block size, including padding:
 *                b_size[0]: useful block size
 *                b_size[1]: vector block extents
 *----------------------------------------------------------------------------*/

static inline void
_b_zero_range(cs_real_t  *restrict y,
              cs_lnum_t   start_id,
              cs_lnum_t   end_id,
              const int   b_size[2])
{
  cs_lnum_t  ii;

# pragma omp parallel for  if((end_id-start_id)*b_size[1] > CS_THR_MIN)
  for (ii = start_id*b_size[1]; ii < end_id*b_size[1]; ii++)
    y[ii] = 0.0;
}

/*----------------------------------------------------------------------------
 * Set values from y[start_id] to y[end_id] to 0, block version.
 *
 * parameters:
 *   y        --> resulting vector
 *   start_id <-- start id in array
 *   end_id   <-- end id in array
 *----------------------------------------------------------------------------*/

static inline void
_3_3_zero_range(cs_real_t  *restrict y,
                cs_lnum_t   start_id,
                cs_lnum_t   end_id)
{
  cs_lnum_t  ii;

# pragma omp parallel for  if((end_id-start_id)*3 > CS_THR_MIN)
  for (ii = start_id*3; ii < end_id*3; ii++)
    y[ii] = 0.0;
}

/*----------------------------------------------------------------------------
 * Descend binary tree for the ordering of a cs_lnum_t (integer) array.
 *
 * parameters:
 *   number    <-> pointer to elements that should be ordered
 *   level     <-- level of the binary tree to descend
 *   n_elts    <-- number of elements in the binary tree to descend
 *----------------------------------------------------------------------------*/

inline static void
_sort_descend_tree(cs_lnum_t  number[],
                   size_t     level,
                   size_t     n_elts)
{
  size_t lv_cur;
  cs_lnum_t num_save;

  num_save = number[level];

  while (level <= (n_elts/2)) {

    lv_cur = (2*level) + 1;

    if (lv_cur < n_elts - 1)
      if (number[lv_cur+1] > number[lv_cur]) lv_cur++;

    if (lv_cur >= n_elts) break;

    if (num_save >= number[lv_cur]) break;

    number[level] = number[lv_cur];
    level = lv_cur;

  }

  number[level] = num_save;
}

/*----------------------------------------------------------------------------
 * Order an array of global numbers.
 *
 * parameters:
 *   number   <-> number of arrays to sort
 *   n_elts   <-- number of elements considered
 *----------------------------------------------------------------------------*/

static void
_sort_local(cs_lnum_t  number[],
            size_t     n_elts)
{
  size_t i, j, inc;
  cs_lnum_t num_save;

  if (n_elts < 2)
    return;

  /* Use shell sort for short arrays */

  if (n_elts < 20) {

    /* Compute increment */
    for (inc = 1; inc <= n_elts/9; inc = 3*inc+1);

    /* Sort array */
    while (inc > 0) {
      for (i = inc; i < n_elts; i++) {
        num_save = number[i];
        j = i;
        while (j >= inc && number[j-inc] > num_save) {
          number[j] = number[j-inc];
          j -= inc;
        }
        number[j] = num_save;
      }
      inc = inc / 3;
    }

  }

  else {

    /* Create binary tree */

    i = (n_elts / 2);
    do {
      i--;
      _sort_descend_tree(number, i, n_elts);
    } while (i > 0);

    /* Sort binary tree */

    for (i = n_elts - 1 ; i > 0 ; i--) {
      num_save   = number[0];
      number[0] = number[i];
      number[i] = num_save;
      _sort_descend_tree(number, 0, i);
    }
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
 *   n_cells     <-- Local number of participating cells
 *   n_cells_ext <-- Local number of cells + ghost cells sharing a face
 *   n_faces     <-- Local number of faces
 *   face_cell   <-- Face -> cells connectivity
 *
 * returns:
 *   pointer to allocated native matrix structure.
 *----------------------------------------------------------------------------*/

static cs_matrix_struct_native_t *
_create_struct_native(cs_lnum_t           n_cells,
                      cs_lnum_t           n_cells_ext,
                      cs_lnum_t           n_faces,
                      const cs_lnum_2_t  *face_cell)
{
  cs_matrix_struct_native_t  *ms;

  /* Allocate and map */

  BFT_MALLOC(ms, 1, cs_matrix_struct_native_t);

  /* Allocate and map */

  ms->n_cells = n_cells;
  ms->n_cells_ext = n_cells_ext;
  ms->n_faces = n_faces;

  ms->face_cell = face_cell;

  return ms;
}

/*----------------------------------------------------------------------------
 * Destroy native matrix structure.
 *
 * parameters:
 *   matrix  <->  Pointer to native matrix structure pointer
 *----------------------------------------------------------------------------*/

static void
_destroy_struct_native(cs_matrix_struct_native_t  **matrix)
{
  if (matrix != NULL && *matrix !=NULL) {

    BFT_FREE(*matrix);

  }
}

/*----------------------------------------------------------------------------
 * Create native matrix coefficients.
 *
 * returns:
 *   pointer to allocated native coefficients structure.
 *----------------------------------------------------------------------------*/

static cs_matrix_coeff_native_t *
_create_coeff_native(void)
{
  cs_matrix_coeff_native_t  *mc;

  /* Allocate */

  BFT_MALLOC(mc, 1, cs_matrix_coeff_native_t);

  /* Initialize */

  mc->symmetric = false;
  mc->max_db_size = 0;
  mc->max_eb_size = 0;

  mc->da = NULL;
  mc->xa = NULL;

  mc->_da = NULL;
  mc->_xa = NULL;

  return mc;
}

/*----------------------------------------------------------------------------
 * Destroy native matrix coefficients.
 *
 * parameters:
 *   coeff  <->  Pointer to native matrix coefficients pointer
 *----------------------------------------------------------------------------*/

static void
_destroy_coeff_native(cs_matrix_coeff_native_t **coeff)
{
  if (coeff != NULL && *coeff !=NULL) {

    cs_matrix_coeff_native_t  *mc = *coeff;

    if (mc->_xa != NULL)
      BFT_FREE(mc->_xa);

    if (mc->_da != NULL)
      BFT_FREE(mc->_da);

    BFT_FREE(*coeff);

  }
}

/*----------------------------------------------------------------------------
 * Set Native matrix coefficients.
 *
 * Depending on current options and initialization, values will be copied
 * or simply mapped.
 *
 * parameters:
 *   matrix           <-- Pointer to matrix structure
 *   symmetric        <-- Indicates if extradiagonal values are symmetric
 *   interleaved      <-- Indicates if matrix coefficients are interleaved
 *   copy             <-- Indicates if coefficients should be copied
 *   da               <-- Diagonal values
 *   xa               <-- Extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_coeffs_native(cs_matrix_t      *matrix,
                   bool              symmetric,
                   bool              interleaved,
                   bool              copy,
                   const cs_real_t  *da,
                   const cs_real_t  *xa)
{
  cs_matrix_coeff_native_t  *mc = matrix->coeffs;
  const cs_matrix_struct_native_t  *ms = matrix->structure;
  cs_lnum_t ii;
  mc->symmetric = symmetric;

  /* Map or copy values */

  if (da != NULL) {

    if (copy) {
      if (mc->_da == NULL || mc->max_db_size < matrix->db_size[3]) {
        BFT_REALLOC(mc->_da, matrix->db_size[3]*ms->n_cells, cs_real_t);
        mc->max_db_size = matrix->db_size[3];
      }
      memcpy(mc->_da, da, matrix->db_size[3]*sizeof(cs_real_t) * ms->n_cells);
      mc->da = mc->_da;
    }
    else
      mc->da = da;

  }
  else {
    mc->da = NULL;
  }

  if (xa != NULL) {

    if (interleaved || symmetric == true) {

      size_t xa_n_vals = ms->n_faces;
      if (! symmetric)
        xa_n_vals *= 2;

      if (copy) {
        if (mc->_xa == NULL || mc->max_eb_size < matrix->eb_size[3]) {
          BFT_MALLOC(mc->_xa, matrix->eb_size[3]*xa_n_vals, cs_real_t);
          mc->max_eb_size = matrix->eb_size[3];
        }
        memcpy(mc->_xa, xa, matrix->eb_size[3]*xa_n_vals*sizeof(cs_real_t));
        mc->xa = mc->_xa;
      }
      else
        mc->xa = xa;

    }
    else { /* !interleaved && symmetric == false */

      assert(matrix->db_size[3] == 1);
      assert(matrix->eb_size[3] == 1);

      if (mc->_xa == NULL)
        BFT_MALLOC(mc->_xa, 2*ms->n_faces, cs_real_t);

      for (ii = 0; ii < ms->n_faces; ++ii) {
        mc->_xa[2*ii] = xa[ii];
        mc->_xa[2*ii + 1] = xa[ms->n_faces + ii];
      }
      mc->xa = mc->_xa;

    }
  }
}

/*----------------------------------------------------------------------------
 * Release shared native matrix coefficients.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *----------------------------------------------------------------------------*/

static void
_release_coeffs_native(cs_matrix_t  *matrix)
{
  cs_matrix_coeff_native_t  *mc = matrix->coeffs;
  if (mc !=NULL) {
    mc->da = NULL;
    mc->xa = NULL;
  }
}

/*----------------------------------------------------------------------------
 * Copy diagonal of native or MSR matrix.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *   da     --> Diagonal (pre-allocated, size: n_cells)
 *----------------------------------------------------------------------------*/

static void
_copy_diagonal_separate(const cs_matrix_t  *matrix,
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
 * Local matrix.vector product y = A.x with native matrix.
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- Pointer to matrix structure
 *   x            <-- Multipliying vector values
 *   y            --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_p_l_native(bool                exclude_diag,
                    const cs_matrix_t  *matrix,
                    const cs_real_t    *restrict x,
                    cs_real_t          *restrict y)
{
  cs_lnum_t  ii, jj, face_id;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;

  const cs_real_t  *restrict xa = mc->xa;

  /* Tell IBM compiler not to alias */
# if defined(__xlc__)
# pragma disjoint(*x, *y, *xa)
# endif

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _diag_vec_p_l(mc->da, x, y, ms->n_cells);
    _zero_range(y, ms->n_cells, ms->n_cells_ext);
  }
  else
    _zero_range(y, 0, ms->n_cells_ext);

  /* Note: parallel and periodic synchronization could be delayed to here */

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    const cs_lnum_2_t *restrict face_cel_p = ms->face_cell;

    if (mc->symmetric) {

      for (face_id = 0; face_id < ms->n_faces; face_id++) {
        ii = face_cel_p[face_id][0];
        jj = face_cel_p[face_id][1];
        y[ii] += xa[face_id] * x[jj];
        y[jj] += xa[face_id] * x[ii];
      }

    }
    else {

      for (face_id = 0; face_id < ms->n_faces; face_id++) {
        ii = face_cel_p[face_id][0];
        jj = face_cel_p[face_id][1];
        y[ii] += xa[2*face_id] * x[jj];
        y[jj] += xa[2*face_id + 1] * x[ii];
      }

    }

  }
}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with native matrix.
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- Pointer to matrix structure
 *   x            <-- Multipliying vector values
 *   y            --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_b_mat_vec_p_l_native(bool                exclude_diag,
                      const cs_matrix_t  *matrix,
                      const cs_real_t    *restrict x,
                      cs_real_t          *restrict y)
{
  cs_lnum_t  ii, jj, kk, face_id;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;

  const cs_real_t  *restrict xa = mc->xa;
  const int *db_size = matrix->db_size;

  /* Tell IBM compiler not to alias */
# if defined(__xlc__)
# pragma disjoint(*x, *y, *xa)
# endif

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _b_diag_vec_p_l(mc->da, x, y, ms->n_cells, db_size);
    _b_zero_range(y, ms->n_cells, ms->n_cells_ext, db_size);
  }
  else
    _b_zero_range(y, 0, ms->n_cells_ext, db_size);

  /* Note: parallel and periodic synchronization could be delayed to here */

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    const cs_lnum_2_t *restrict face_cel_p = ms->face_cell;

    if (mc->symmetric) {

      for (face_id = 0; face_id < ms->n_faces; face_id++) {
        ii = face_cel_p[face_id][0];
        jj = face_cel_p[face_id][1];
        for (kk = 0; kk < db_size[0]; kk++) {
          y[ii*db_size[1] + kk] += xa[face_id] * x[jj*db_size[1] + kk];
          y[jj*db_size[1] + kk] += xa[face_id] * x[ii*db_size[1] + kk];
        }
      }
    }
    else {

      for (face_id = 0; face_id < ms->n_faces; face_id++) {
        ii = face_cel_p[face_id][0];
        jj = face_cel_p[face_id][1];
        for (kk = 0; kk < db_size[0]; kk++) {
          y[ii*db_size[1] + kk] += xa[2*face_id]     * x[jj*db_size[1] + kk];
          y[jj*db_size[1] + kk] += xa[2*face_id + 1] * x[ii*db_size[1] + kk];
        }
      }

    }

  }

}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with native matrix.
 *
 * This variant uses a fixed 3x3 block, for better compiler optimization.
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- Pointer to matrix structure
 *   x            <-- Multipliying vector values
 *   y            --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_3_3_mat_vec_p_l_native(bool                exclude_diag,
                        const cs_matrix_t  *matrix,
                        const cs_real_t    *restrict x,
                        cs_real_t          *restrict y)
{
  cs_lnum_t  ii, jj, kk, face_id;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;

  const cs_real_t  *restrict xa = mc->xa;

  assert(matrix->db_size[0] == 3 && matrix->db_size[3] == 9);

  /* Tell IBM compiler not to alias */
# if defined(__xlc__)
# pragma disjoint(*x, *y, *xa)
# endif

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _3_3_diag_vec_p_l(mc->da, x, y, ms->n_cells);
    _3_3_zero_range(y, ms->n_cells, ms->n_cells_ext);
  }
  else
    _3_3_zero_range(y, 0, ms->n_cells_ext);

  /* Note: parallel and periodic synchronization could be delayed to here */

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    const cs_lnum_2_t *restrict face_cel_p = ms->face_cell;

    if (mc->symmetric) {

      for (face_id = 0; face_id < ms->n_faces; face_id++) {
        ii = face_cel_p[face_id][0];
        jj = face_cel_p[face_id][1];
        for (kk = 0; kk < 3; kk++) {
          y[ii*3 + kk] += xa[face_id] * x[jj*3 + kk];
          y[jj*3 + kk] += xa[face_id] * x[ii*3 + kk];
        }
      }
    }
    else {

      for (face_id = 0; face_id < ms->n_faces; face_id++) {
        ii = face_cel_p[face_id][0];
        jj = face_cel_p[face_id][1];
        for (kk = 0; kk < 3; kk++) {
          y[ii*3 + kk] += xa[2*face_id]     * x[jj*3 + kk];
          y[jj*3 + kk] += xa[2*face_id + 1] * x[ii*3 + kk];
        }
      }

    }

  }

}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with native matrix.
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- Pointer to matrix structure
 *   x            <-- Multipliying vector values
 *   y            --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_bb_mat_vec_p_l_native(bool                exclude_diag,
                       const cs_matrix_t  *matrix,
                       const cs_real_t    *restrict x,
                       cs_real_t          *restrict y)
{
  cs_lnum_t  ii, jj, face_id;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;

  const cs_real_t  *restrict xa = mc->xa;
  const int *db_size = matrix->db_size;
  const int *eb_size = matrix->eb_size;

  /* Tell IBM compiler not to alias */
# if defined(__xlc__)
# pragma disjoint(*x, *y, *xa)
# endif

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _b_diag_vec_p_l(mc->da, x, y, ms->n_cells, db_size);
    _b_zero_range(y, ms->n_cells, ms->n_cells_ext, db_size);
  }
  else
    _b_zero_range(y, 0, ms->n_cells_ext, db_size);

  /* Note: parallel and periodic synchronization could be delayed to here */

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    const cs_lnum_2_t *restrict face_cel_p = ms->face_cell;

    if (mc->symmetric) {

      for (face_id = 0; face_id < ms->n_faces; face_id++) {
        ii = face_cel_p[face_id][0];
        jj = face_cel_p[face_id][1];
        _dense_eb_ax_add(ii, jj, face_id, eb_size, xa, x, y);
        _dense_eb_ax_add(jj, ii, face_id, eb_size, xa, x, y);
      }
    }
    else {

      for (face_id = 0; face_id < ms->n_faces; face_id++) {
        ii = face_cel_p[face_id][0];
        jj = face_cel_p[face_id][1];
        _dense_eb_ax_add(ii, jj, 2*face_id, eb_size, xa, x, y);
        _dense_eb_ax_add(jj, ii, 2*face_id + 1, eb_size, xa, x, y);
      }

    }

  }

}

#if defined(HAVE_OPENMP) /* OpenMP variants */

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with native matrix.
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- Pointer to matrix structure
 *   x            <-- Multipliying vector values
 *   y            --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_p_l_native_omp(bool                exclude_diag,
                        const cs_matrix_t  *matrix,
                        const cs_real_t    *restrict x,
                        cs_real_t          *restrict y)
{
  const int n_threads = matrix->numbering->n_threads;
  const int n_groups = matrix->numbering->n_groups;
  const cs_lnum_t *group_index = matrix->numbering->group_index;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;
  const cs_real_t  *restrict xa = mc->xa;

  assert(matrix->numbering->type == CS_NUMBERING_THREADS);

  /* Tell IBM compiler not to alias */

# if defined(__xlc__)
# pragma disjoint(*x, *y, *xa)
# endif

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _diag_vec_p_l(mc->da, x, y, ms->n_cells);
    _zero_range(y, ms->n_cells, ms->n_cells_ext);
  }
  else
    _zero_range(y, 0, ms->n_cells_ext);

  /* Note: parallel and periodic synchronization could be delayed to here */

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    const cs_lnum_2_t *restrict face_cel_p = ms->face_cell;

    if (mc->symmetric) {

      for (int g_id = 0; g_id < n_groups; g_id++) {

#       pragma omp parallel for
        for (int t_id = 0; t_id < n_threads; t_id++) {

          for (cs_lnum_t face_id = group_index[(t_id*n_groups + g_id)*2];
               face_id < group_index[(t_id*n_groups + g_id)*2 + 1];
               face_id++) {
            cs_lnum_t ii = face_cel_p[face_id][0];
            cs_lnum_t jj = face_cel_p[face_id][1];
            y[ii] += xa[face_id] * x[jj];
            y[jj] += xa[face_id] * x[ii];
          }
        }
      }
    }
    else {

      for (int g_id = 0; g_id < n_groups; g_id++) {

#       pragma omp parallel for
        for (int t_id = 0; t_id < n_threads; t_id++) {

          for (cs_lnum_t face_id = group_index[(t_id*n_groups + g_id)*2];
               face_id < group_index[(t_id*n_groups + g_id)*2 + 1];
               face_id++) {
            cs_lnum_t ii = face_cel_p[face_id][0];
            cs_lnum_t jj = face_cel_p[face_id][1];
            y[ii] += xa[2*face_id] * x[jj];
            y[jj] += xa[2*face_id + 1] * x[ii];
          }
        }
      }
    }

  }
}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with native matrix, blocked version
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- Pointer to matrix structure
 *   x            <-- Multipliying vector values
 *   y            --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_b_mat_vec_p_l_native_omp(bool                exclude_diag,
                          const cs_matrix_t  *matrix,
                          const cs_real_t    *restrict x,
                          cs_real_t          *restrict y)
{
  const int *db_size = matrix->db_size;

  const int n_threads = matrix->numbering->n_threads;
  const int n_groups = matrix->numbering->n_groups;
  const cs_lnum_t *group_index = matrix->numbering->group_index;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;
  const cs_real_t  *restrict xa = mc->xa;

  assert(matrix->numbering->type == CS_NUMBERING_THREADS);

  /* Tell IBM compiler not to alias */

# if defined(__xlc__)
# pragma disjoint(*x, *y, *xa)
# endif

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _b_diag_vec_p_l(mc->da, x, y, ms->n_cells, db_size);
    _b_zero_range(y, ms->n_cells, ms->n_cells_ext, db_size);
  }
  else
    _b_zero_range(y, 0, ms->n_cells_ext, db_size);

  /* Note: parallel and periodic synchronization could be delayed to here */

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    const cs_lnum_2_t *restrict face_cel_p = ms->face_cell;

    if (mc->symmetric) {

      for (int g_id = 0; g_id < n_groups; g_id++) {

#       pragma omp parallel for
        for (int t_id = 0; t_id < n_threads; t_id++) {

          for (cs_lnum_t face_id = group_index[(t_id*n_groups + g_id)*2];
               face_id < group_index[(t_id*n_groups + g_id)*2 + 1];
               face_id++) {
            cs_lnum_t ii = face_cel_p[face_id][0];
            cs_lnum_t jj = face_cel_p[face_id][1];
            for (cs_lnum_t kk = 0; kk < db_size[0]; kk++) {
              y[ii*db_size[1] + kk] += xa[face_id] * x[jj*db_size[1] + kk];
              y[jj*db_size[1] + kk] += xa[face_id] * x[ii*db_size[1] + kk];
            }
          }
        }
      }

    }
    else {

      for (int g_id = 0; g_id < n_groups; g_id++) {

#       pragma omp parallel for
        for (int t_id = 0; t_id < n_threads; t_id++) {

          for (cs_lnum_t face_id = group_index[(t_id*n_groups + g_id)*2];
               face_id < group_index[(t_id*n_groups + g_id)*2 + 1];
               face_id++) {
            cs_lnum_t ii = face_cel_p[face_id][0];
            cs_lnum_t jj = face_cel_p[face_id][1];
            for (cs_lnum_t kk = 0; kk < db_size[0]; kk++) {
              y[ii*db_size[1] + kk] += xa[2*face_id]     * x[jj*db_size[1] + kk];
              y[jj*db_size[1] + kk] += xa[2*face_id + 1] * x[ii*db_size[1] + kk];
            }
          }
        }
      }

    }

  }
}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with native matrix.
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- Pointer to matrix structure
 *   x            <-- Multipliying vector values
 *   y            --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_p_l_native_omp_atomic(bool                exclude_diag,
                               const cs_matrix_t  *matrix,
                               const cs_real_t    *restrict x,
                               cs_real_t          *restrict y)
{
  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;
  const cs_real_t  *restrict xa = mc->xa;

  /* Tell IBM compiler not to alias */

# if defined(__xlc__)
# pragma disjoint(*x, *y, *xa)
# endif

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _diag_vec_p_l(mc->da, x, y, ms->n_cells);
    _zero_range(y, ms->n_cells, ms->n_cells_ext);
  }
  else
    _zero_range(y, 0, ms->n_cells_ext);

  /* Note: parallel and periodic synchronization could be delayed to here */

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    const cs_lnum_2_t *restrict face_cel_p = ms->face_cell;

    if (mc->symmetric) {

#     pragma omp parallel for
      for (cs_lnum_t face_id = 0; face_id < ms->n_faces; face_id++) {
        cs_lnum_t ii = face_cel_p[face_id][0];
        cs_lnum_t jj = face_cel_p[face_id][1];
#       pragma omp atomic
        y[ii] += xa[face_id] * x[jj];
#       pragma omp atomic
        y[jj] += xa[face_id] * x[ii];
      }
    }
    else {

#     pragma omp parallel for
      for (cs_lnum_t face_id = 0; face_id < ms->n_faces; face_id++) {
        cs_lnum_t ii = face_cel_p[face_id][0];
        cs_lnum_t jj = face_cel_p[face_id][1];
#       pragma omp atomic
        y[ii] += xa[2*face_id] * x[jj];
#       pragma omp atomic
        y[jj] += xa[2*face_id + 1] * x[ii];
      }
    }

  }
}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with native matrix, blocked version
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- Pointer to matrix structure
 *   x            <-- Multipliying vector values
 *   y            --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_b_mat_vec_p_l_native_omp_atomic(bool                exclude_diag,
                                 const cs_matrix_t  *matrix,
                                 const cs_real_t    *restrict x,
                                 cs_real_t          *restrict y)
{
  const int *db_size = matrix->db_size;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;
  const cs_real_t  *restrict xa = mc->xa;

  assert(matrix->numbering->type == CS_NUMBERING_THREADS);

  /* Tell IBM compiler not to alias */

# if defined(__xlc__)
# pragma disjoint(*x, *y, *xa)
# endif

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _b_diag_vec_p_l(mc->da, x, y, ms->n_cells, db_size);
    _b_zero_range(y, ms->n_cells, ms->n_cells_ext, db_size);
  }
  else
    _b_zero_range(y, 0, ms->n_cells_ext, db_size);

  /* Note: parallel and periodic synchronization could be delayed to here */

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    const cs_lnum_2_t *restrict face_cel_p = ms->face_cell;

    if (mc->symmetric) {

#     pragma omp parallel for
      for (cs_lnum_t face_id = 0; face_id < ms->n_faces; face_id++) {
        cs_lnum_t ii = face_cel_p[face_id][0];
        cs_lnum_t jj = face_cel_p[face_id][1];
        for (cs_lnum_t kk = 0; kk < db_size[0]; kk++) {
#         pragma omp atomic
          y[ii*db_size[1] + kk] += xa[face_id] * x[jj*db_size[1] + kk];
#         pragma omp atomic
          y[jj*db_size[1] + kk] += xa[face_id] * x[ii*db_size[1] + kk];
        }
      }

    }
    else {

#     pragma omp parallel for
      for (cs_lnum_t face_id = 0; face_id < ms->n_faces; face_id++) {
        cs_lnum_t ii = face_cel_p[face_id][0];
        cs_lnum_t jj = face_cel_p[face_id][1];
        for (cs_lnum_t kk = 0; kk < db_size[0]; kk++) {
#         pragma omp atomic
          y[ii*db_size[1] + kk] += xa[2*face_id]   * x[jj*db_size[1] + kk];
#         pragma omp atomic
          y[jj*db_size[1] + kk] += xa[2*face_id+1] * x[ii*db_size[1] + kk];
        }
      }

    }

  }
}

#endif /* defined(HAVE_OPENMP) */

static void
_mat_vec_p_l_native_bull(bool                exclude_diag,
                         const cs_matrix_t  *matrix,
                         const cs_real_t    *restrict x,
                         cs_real_t          *restrict y)
{
  cs_lnum_t  ii, ii_prev, kk, face_id, kk_max;
  cs_real_t y_it, y_it_prev;
  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;
  const int ed_flag = (exclude_diag) ? 1 : 0;
  const cs_real_t  *restrict xa = mc->xa;
  const int l1_cache_size
    = (matrix->loop_length[matrix->fill_type][ed_flag] > 0) ?
        matrix->loop_length[matrix->fill_type][ed_flag] : 508;

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _diag_vec_p_l(mc->da, x, y, ms->n_cells);
    /* _zero_range(y, ms->n_cells, ms->n_cells_ext); */
  }
  else
    _zero_range(y, 0, ms->n_cells_ext);

  for (ii = ms->n_cells; ii < ms->n_cells_ext; y[ii++] = 0.0);

  /* Note: parallel and periodic synchronization could be delayed to here */

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    /*
     * 1/ Split y[ii] and y[jj] computation into 2 loops to remove compiler
     *    data dependency assertion between y[ii] and y[jj].
     * 2/ keep index (*face_cel_p) in L1 cache from y[ii] loop to y[jj] loop
     *    and xa in L2 cache.
     * 3/ break high frequency occurence of data dependency from one iteration
     *    to another in y[ii] loop (nonzero matrix value on the same line ii).
     */

    const cs_lnum_2_t *restrict face_cel_p = ms->face_cell;

    if (mc->symmetric) {

      for (face_id = 0;
           face_id < ms->n_faces;
           face_id += l1_cache_size) {

        kk_max = CS_MIN((ms->n_faces - face_id), l1_cache_size);

        /* sub-loop to compute y[ii] += xa[face_id] * x[jj] */

        ii = face_cel_p[0][0];
        ii_prev = ii;
        y_it_prev = y[ii_prev] + xa[face_id] * x[face_cel_p[0][1]];

        for (kk = 1; kk < kk_max; ++kk) {
          ii = face_cel_p[kk][0];
          /* y[ii] += xa[face_id+kk] * x[jj]; */
          if (ii == ii_prev) {
            y_it = y_it_prev;
          }
          else {
            y_it = y[ii];
            y[ii_prev] = y_it_prev;
          }
          ii_prev = ii;
          y_it_prev = y_it + xa[face_id+kk] * x[face_cel_p[kk][1]];
        }
        y[ii] = y_it_prev;

        /* sub-loop to compute y[ii] += xa[face_id] * x[jj] */

        for (kk = 0; kk < kk_max; ++kk) {
          y[face_cel_p[kk][1]] += xa[face_id+kk] * x[face_cel_p[kk][0]];
        }
        face_cel_p += l1_cache_size;
      }

    }
    else {

      for (face_id = 0;
           face_id < ms->n_faces;
           face_id+=l1_cache_size) {

        kk_max = CS_MIN((ms->n_faces - face_id),
                        l1_cache_size);

        /* sub-loop to compute y[ii] += xa[2*face_id] * x[jj] */

        ii = face_cel_p[0][0];
        ii_prev = ii;
        y_it_prev = y[ii_prev] + xa[2*face_id] * x[face_cel_p[0][1]];

        for (kk = 1; kk < kk_max; ++kk) {
          ii = face_cel_p[kk][0];
          /* y[ii] += xa[2*(face_id+i)] * x[jj]; */
          if (ii == ii_prev) {
            y_it = y_it_prev;
          }
          else {
            y_it = y[ii];
            y[ii_prev] = y_it_prev;
          }
          ii_prev = ii;
          y_it_prev = y_it + xa[2*(face_id+kk)] * x[face_cel_p[kk][1]];
        }
        y[ii] = y_it_prev;

        /* sub-loop to compute y[ii] += xa[2*face_id + 1] * x[jj] */

        for (kk = 0; kk < kk_max; ++kk) {
          y[face_cel_p[kk][1]] += xa[2*(face_id+kk) + 1] * x[face_cel_p[kk][0]];
        }
        face_cel_p += l1_cache_size;
      }

    }
  }
}

#if defined(SX) && defined(_SX) /* For vector machines */

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with native matrix.
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- Pointer to matrix structure
 *   x            <-- Multipliying vector values
 *   y            --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_p_l_native_vector(bool                exclude_diag,
                           const cs_matrix_t  *matrix,
                           const cs_real_t    *restrict x,
                           cs_real_t          *restrict y)
{
  cs_lnum_t  ii, jj, face_id;
  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;
  const cs_real_t  *restrict xa = mc->xa;

  assert(matrix->numbering->type == CS_NUMBERING_VECTORIZE);

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _diag_vec_p_l(mc->da, x, y, ms->n_cells);
    _zero_range(y, ms->n_cells, ms->n_cells_ext);
  }
  else
    _zero_range(y, 0, ms->n_cells_ext);

  /* Note: parallel and periodic synchronization could be delayed to here */

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    const cs_lnum_2_t *restrict face_cel_p = ms->face_cell;

    if (mc->symmetric) {

#     if defined(HAVE_OPENMP)
#       pragma omp simd safelen(CS_NUMBERING_SIMD_SIZE)
#     else
#       pragma dir nodep
#       pragma GCC ivdep
#     endif
      for (face_id = 0; face_id < ms->n_faces; face_id++) {
        ii = face_cel_p[face_id][0];
        jj = face_cel_p[face_id][1];
        y[ii] += xa[face_id] * x[jj];
        y[jj] += xa[face_id] * x[ii];
      }

    }
    else {

#     if defined(HAVE_OPENMP)
#       pragma omp simd safelen(CS_NUMBERING_SIMD_SIZE)
#     else
#       pragma dir nodep
#       pragma GCC ivdep
#     endif
      for (face_id = 0; face_id < ms->n_faces; face_id++) {
        ii = face_cel_p[face_id][0];
        jj = face_cel_p[face_id][1];
        y[ii] += xa[2*face_id] * x[jj];
        y[jj] += xa[2*face_id + 1] * x[ii];
      }

    }

  }
}

#endif /* Vector machine variant */

/*----------------------------------------------------------------------------
 * Create a CSR matrix structure from a native matrix stucture.
 *
 * Note that the structure created maps global cell numbers to the given
 * existing face -> cell connectivity array, so it must be destroyed before
 * this array (usually the code's global cell numbering) is freed.
 *
 * parameters:
 *   have_diag   <-- Indicates if the diagonal is nonzero
 *   n_cells     <-- Local number of participating cells
 *   n_cells_ext <-- Local number of cells + ghost cells sharing a face
 *   n_faces     <-- Local number of faces
 *   face_cell   <-- Face -> cells connectivity
 *
 * returns:
 *   pointer to allocated CSR matrix structure.
 *----------------------------------------------------------------------------*/

static cs_matrix_struct_csr_t *
_create_struct_csr(bool                have_diag,
                   cs_lnum_t           n_cells,
                   cs_lnum_t           n_cells_ext,
                   cs_lnum_t           n_faces,
                   const cs_lnum_2_t  *face_cell)
{
  int n_cols_max;
  cs_lnum_t ii, jj, face_id;
  const cs_lnum_t *restrict face_cel_p;

  cs_lnum_t  diag_elts = 1;
  cs_lnum_t  *ccount = NULL;

  cs_matrix_struct_csr_t  *ms;

  /* Allocate and map */

  BFT_MALLOC(ms, 1, cs_matrix_struct_csr_t);

  ms->n_rows = n_cells;
  ms->n_cols = n_cells_ext;

  ms->direct_assembly = true;
  ms->have_diag = have_diag;

  BFT_MALLOC(ms->row_index, ms->n_rows + 1, cs_lnum_t);

  /* Count number of nonzero elements per row */

  BFT_MALLOC(ccount, ms->n_cols, cs_lnum_t);

  if (have_diag == false)
    diag_elts = 0;

  for (ii = 0; ii < ms->n_rows; ii++)  /* count starting with diagonal terms */
    ccount[ii] = diag_elts;

  if (face_cell != NULL) {

    face_cel_p = (const cs_lnum_t *restrict)face_cell;

    for (face_id = 0; face_id < n_faces; face_id++) {
      ii = *face_cel_p++;
      jj = *face_cel_p++;
      ccount[ii] += 1;
      ccount[jj] += 1;
    }

  } /* if (face_cell != NULL) */

  n_cols_max = 0;

  ms->row_index[0] = 0;
  for (ii = 0; ii < ms->n_rows; ii++) {
    ms->row_index[ii+1] = ms->row_index[ii] + ccount[ii];
    if (ccount[ii] > n_cols_max)
      n_cols_max = ccount[ii];
    ccount[ii] = diag_elts; /* pre-count for diagonal terms */
  }

  ms->n_cols_max = n_cols_max;

  /* Build structure */

  BFT_MALLOC(ms->col_id, (ms->row_index[ms->n_rows]), cs_lnum_t);

  if (have_diag == true) {
    for (ii = 0; ii < ms->n_rows; ii++) {    /* diagonal terms */
      ms->col_id[ms->row_index[ii]] = ii;
    }
  }

  if (face_cell != NULL) {                   /* non-diagonal terms */

    face_cel_p = (const cs_lnum_t *restrict)face_cell;

    for (face_id = 0; face_id < n_faces; face_id++) {
      ii = *face_cel_p++;
      jj = *face_cel_p++;
      if (ii < ms->n_rows) {
        ms->col_id[ms->row_index[ii] + ccount[ii]] = jj;
        ccount[ii] += 1;
      }
      if (jj < ms->n_rows) {
        ms->col_id[ms->row_index[jj] + ccount[jj]] = ii;
        ccount[jj] += 1;
      }
    }

  } /* if (face_cell != NULL) */

  BFT_FREE(ccount);

  /* Sort line elements by column id (for better access patterns) */

  if (n_cols_max > 1) {

    for (ii = 0; ii < ms->n_rows; ii++) {
      cs_lnum_t *col_id = ms->col_id + ms->row_index[ii];
      cs_lnum_t n_cols = ms->row_index[ii+1] - ms->row_index[ii];
      cs_lnum_t col_id_prev = -1;
      _sort_local(col_id, ms->row_index[ii+1] - ms->row_index[ii]);
      for (jj = 0; jj < n_cols; jj++) {
        if (col_id[jj] == col_id_prev)
          ms->direct_assembly = false;
        col_id_prev = col_id[jj];
      }
    }

  }

  /* Compact elements if necessary */

  if (ms->direct_assembly == false) {

    cs_lnum_t *tmp_row_index = NULL;
    cs_lnum_t  kk = 0;

    BFT_MALLOC(tmp_row_index, ms->n_rows+1, cs_lnum_t);
    memcpy(tmp_row_index, ms->row_index, (ms->n_rows+1)*sizeof(cs_lnum_t));

    kk = 0;

    for (ii = 0; ii < ms->n_rows; ii++) {
      cs_lnum_t *col_id = ms->col_id + ms->row_index[ii];
      cs_lnum_t n_cols = ms->row_index[ii+1] - ms->row_index[ii];
      cs_lnum_t col_id_prev = -1;
      ms->row_index[ii] = kk;
      for (jj = 0; jj < n_cols; jj++) {
        if (col_id_prev != col_id[jj]) {
          ms->col_id[kk++] = col_id[jj];
          col_id_prev = col_id[jj];
        }
      }
    }
    ms->row_index[ms->n_rows] = kk;

    assert(ms->row_index[ms->n_rows] < tmp_row_index[ms->n_rows]);

    BFT_FREE(tmp_row_index);
    BFT_REALLOC(ms->col_id, (ms->row_index[ms->n_rows]), cs_lnum_t);

  }

  return ms;
}

/*----------------------------------------------------------------------------
 * Destroy CSR matrix structure.
 *
 * parameters:
 *   matrix  <->  Pointer to CSR matrix structure pointer
 *----------------------------------------------------------------------------*/

static void
_destroy_struct_csr(cs_matrix_struct_csr_t  **matrix)
{
  if (matrix != NULL && *matrix !=NULL) {

    cs_matrix_struct_csr_t  *ms = *matrix;

    if (ms->row_index != NULL)
      BFT_FREE(ms->row_index);

    if (ms->col_id != NULL)
      BFT_FREE(ms->col_id);

    BFT_FREE(ms);

    *matrix = ms;

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
  cs_matrix_coeff_csr_t  *mc;

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
 * Destroy CSR matrix coefficients.
 *
 * parameters:
 *   coeff  <->  Pointer to CSR matrix coefficients pointer
 *----------------------------------------------------------------------------*/

static void
_destroy_coeff_csr(cs_matrix_coeff_csr_t **coeff)
{
  if (coeff != NULL && *coeff !=NULL) {

    cs_matrix_coeff_csr_t  *mc = *coeff;

    if (mc->val != NULL)
      BFT_FREE(mc->val);

    if (mc->x_prefetch != NULL)
      BFT_FREE(mc->x_prefetch);

    if (mc->_d_val != NULL)
      BFT_FREE(mc->_d_val);

    BFT_FREE(*coeff);

  }
}

/*----------------------------------------------------------------------------
 * Set CSR extradiagonal matrix coefficients for the case where direct
 * assignment is possible (i.e. when there are no multiple contributions
 * to a given coefficient).
 *
 * parameters:
 *   matrix      <-- Pointer to matrix structure
 *   symmetric   <-- Indicates if extradiagonal values are symmetric
 *   interleaved <-- Indicates if matrix coefficients are interleaved
 *   xa          <-- Extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_xa_coeffs_csr_direct(cs_matrix_t      *matrix,
                          bool              symmetric,
                          bool              interleaved,
                          const cs_real_t  *restrict xa)
{
  cs_lnum_t  ii, jj, face_id;
  cs_matrix_coeff_csr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  /* Copy extra-diagonal values */

  assert(matrix->face_cell != NULL);

  const cs_lnum_t n_faces = matrix->n_faces;
  const cs_lnum_t *restrict face_cel_p
    = (const cs_lnum_t *restrict)(matrix->face_cell);

  if (symmetric == false) {

    if (interleaved == false) {
      const cs_real_t  *restrict xa1 = xa;
      const cs_real_t  *restrict xa2 = xa + matrix->n_faces;
      for (face_id = 0; face_id < n_faces; face_id++) {
        cs_lnum_t kk, ll;
        ii = *face_cel_p++;
        jj = *face_cel_p++;
        if (ii < ms->n_rows) {
          for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
          mc->val[kk] = xa1[face_id];
        }
        if (jj < ms->n_rows) {
          for (ll = ms->row_index[jj]; ms->col_id[ll] != ii; ll++);
          mc->val[ll] = xa2[face_id];
        }
      }
    }
    else { /* interleaved == true */
      for (face_id = 0; face_id < n_faces; face_id++) {
        cs_lnum_t kk, ll;
        ii = *face_cel_p++;
        jj = *face_cel_p++;
        if (ii < ms->n_rows) {
          for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
          mc->val[kk] = xa[2*face_id];
        }
        if (jj < ms->n_rows) {
          for (ll = ms->row_index[jj]; ms->col_id[ll] != ii; ll++);
          mc->val[ll] = xa[2*face_id + 1];
        }
      }
    }

  }
  else { /* if symmetric == true */

    for (face_id = 0; face_id < n_faces; face_id++) {
      cs_lnum_t kk, ll;
      ii = *face_cel_p++;
      jj = *face_cel_p++;
      if (ii < ms->n_rows) {
        for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
        mc->val[kk] = xa[face_id];
      }
      if (jj < ms->n_rows) {
        for (ll = ms->row_index[jj]; ms->col_id[ll] != ii; ll++);
        mc->val[ll] = xa[face_id];
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
 *   matrix      <-- Pointer to matrix structure
 *   symmetric   <-- Indicates if extradiagonal values are symmetric
 *   interleaved <-- Indicates if matrix coefficients are interleaved
 *   xa          <-- Extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_xa_coeffs_csr_increment(cs_matrix_t      *matrix,
                             bool              symmetric,
                             bool              interleaved,
                             const cs_real_t  *restrict xa)
{
  cs_lnum_t  ii, jj, face_id;
  cs_matrix_coeff_csr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  /* Copy extra-diagonal values */

  assert(matrix->face_cell != NULL);

  const cs_lnum_t n_faces = matrix->n_faces;
  const cs_lnum_t *restrict face_cel_p
    = (const cs_lnum_t *restrict)(matrix->face_cell);

  if (symmetric == false) {

    if (interleaved == false) {
      const cs_real_t  *restrict xa1 = xa;
      const cs_real_t  *restrict xa2 = xa + matrix->n_faces;
      for (face_id = 0; face_id < n_faces; face_id++) {
        cs_lnum_t kk, ll;
        ii = *face_cel_p++;
        jj = *face_cel_p++;
        if (ii < ms->n_rows) {
          for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
          mc->val[kk] += xa1[face_id];
        }
        if (jj < ms->n_rows) {
          for (ll = ms->row_index[jj]; ms->col_id[ll] != ii; ll++);
          mc->val[ll] += xa2[face_id];
        }
      }
    }
    else { /* interleaved == true */
      for (face_id = 0; face_id < n_faces; face_id++) {
        cs_lnum_t kk, ll;
        ii = *face_cel_p++;
        jj = *face_cel_p++;
        if (ii < ms->n_rows) {
          for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
          mc->val[kk] += xa[2*face_id];
        }
        if (jj < ms->n_rows) {
          for (ll = ms->row_index[jj]; ms->col_id[ll] != ii; ll++);
          mc->val[ll] += xa[2*face_id + 1];
        }
      }
    }
  }
  else { /* if symmetric == true */

    for (face_id = 0; face_id < n_faces; face_id++) {
      cs_lnum_t kk, ll;
      ii = *face_cel_p++;
      jj = *face_cel_p++;
      if (ii < ms->n_rows) {
        for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
        mc->val[kk] += xa[face_id];
      }
      if (jj < ms->n_rows) {
        for (ll = ms->row_index[jj]; ms->col_id[ll] != ii; ll++);
        mc->val[ll] += xa[face_id];
      }

    }

  } /* end of condition on coefficients symmetry */

}

/*----------------------------------------------------------------------------
 * Set CSR matrix coefficients.
 *
 * parameters:
 *   matrix           <-> Pointer to matrix structure
 *   symmetric        <-- Indicates if extradiagonal values are symmetric
 *   interleaved      <-- Indicates if matrix coefficients are interleaved
 *   copy             <-- Indicates if coefficients should be copied
 *   da               <-- Diagonal values (NULL if all zero)
 *   xa               <-- Extradiagonal values (NULL if all zero)
 *----------------------------------------------------------------------------*/

static void
_set_coeffs_csr(cs_matrix_t      *matrix,
                bool              symmetric,
                bool              interleaved,
                bool              copy,
                const cs_real_t  *restrict da,
                const cs_real_t  *restrict xa)
{
  cs_lnum_t  ii, jj;
  cs_matrix_coeff_csr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  if (mc->val == NULL)
    BFT_MALLOC(mc->val, ms->row_index[ms->n_rows], cs_real_t);

  /* Initialize coefficients to zero if assembly is incremental */

  if (ms->direct_assembly == false) {
    cs_lnum_t val_size = ms->row_index[ms->n_rows];
    for (ii = 0; ii < val_size; ii++)
      mc->val[ii] = 0.0;
  }

  /* Allocate prefetch buffer */

  mc->n_prefetch_rows
    = CS_MAX(matrix->loop_length[matrix->fill_type][0],
             matrix->loop_length[matrix->fill_type][1]);
  if (mc->n_prefetch_rows > 0 && mc->x_prefetch == NULL) {
    size_t prefetch_size = ms->n_cols_max * mc->n_prefetch_rows;
    size_t matrix_size = matrix->n_cells + (2 * matrix->n_faces);
    if (matrix_size > prefetch_size)
      prefetch_size = matrix_size;
    BFT_REALLOC(mc->x_prefetch, prefetch_size, cs_real_t);
  }

  /* Copy diagonal values */

  if (ms->have_diag == true) {

    if (da != NULL) {
      for (ii = 0; ii < ms->n_rows; ii++) {
        cs_lnum_t kk;
        for (kk = ms->row_index[ii]; ms->col_id[kk] != ii; kk++);
        mc->val[kk] = da[ii];
      }
    }
    else {
      for (ii = 0; ii < ms->n_rows; ii++) {
        cs_lnum_t kk;
        for (kk = ms->row_index[ii]; ms->col_id[kk] != ii; kk++);
        mc->val[kk] = 0.0;
      }
    }

  }

  /* Mark diagonal values as not queried (mc->_d_val not changed) */

  mc->d_val = NULL;

  /* Copy extra-diagonal values */

  if (matrix->face_cell != NULL) {

    if (xa != NULL) {

      if (ms->direct_assembly == true)
        _set_xa_coeffs_csr_direct(matrix, symmetric, interleaved, xa);
      else
        _set_xa_coeffs_csr_increment(matrix, symmetric, interleaved, xa);

    }
    else { /* if (xa == NULL) */

      for (ii = 0; ii < ms->n_rows; ii++) {
        const cs_lnum_t  *restrict col_id = ms->col_id + ms->row_index[ii];
        cs_real_t  *m_row = mc->val + ms->row_index[ii];
        cs_lnum_t  n_cols = ms->row_index[ii+1] - ms->row_index[ii];

        for (jj = 0; jj < n_cols; jj++) {
          if (col_id[jj] != ii)
            m_row[jj] = 0.0;
        }

      }

    }

  } /* (matrix->face_cell != NULL) */

}

/*----------------------------------------------------------------------------
 * Release shared CSR matrix coefficients.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
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
 *   matrix <-- Pointer to matrix structure
 *   da     --> Diagonal (pre-allocated, size: n_rows)
 *----------------------------------------------------------------------------*/

static void
_copy_diagonal_csr(const cs_matrix_t  *matrix,
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

static void
_mat_vec_p_l_csr(bool                exclude_diag,
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

static void
_mat_vec_p_l_csr_mkl(bool                exclude_diag,
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
 * Local matrix.vector product y = A.x with CSR matrix (prefetch).
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- Pointer to matrix structure
 *   x            <-- Multipliying vector values
 *   y            --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_p_l_csr_pf(bool                exclude_diag,
                    const cs_matrix_t  *matrix,
                    const cs_real_t    *restrict x,
                    cs_real_t          *restrict y)
{
  cs_lnum_t  start_row, ii, jj, n_cols;
  cs_lnum_t  *restrict col_id;
  cs_real_t  *restrict m_row;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_csr_t  *mc = matrix->coeffs;
  cs_lnum_t  n_rows = ms->n_rows;

  if (exclude_diag)
    bft_error(__FILE__, __LINE__, 0,
              _(_no_exclude_diag_error_str), __func__);

  /* Outer loop on prefetch lines */

  for (start_row = 0; start_row < n_rows; start_row += mc->n_prefetch_rows) {

    cs_lnum_t end_row = start_row + mc->n_prefetch_rows;

    cs_real_t  *restrict prefetch_p = mc->x_prefetch;

    /* Tell IBM compiler not to alias */
#   if defined(__xlc__)
#   pragma disjoint(*prefetch_p, *y, *m_row)
#   pragma disjoint(*prefetch_p, *x, *col_id)
#   endif

    if (end_row > n_rows)
      end_row = n_rows;

    /* Prefetch */

    for (ii = start_row; ii < end_row; ii++) {

      col_id = ms->col_id + ms->row_index[ii];
      n_cols = ms->row_index[ii+1] - ms->row_index[ii];

      for (jj = 0; jj < n_cols; jj++)
        *prefetch_p++ = x[col_id[jj]];

    }

    /* Compute */

    prefetch_p = mc->x_prefetch;

    for (ii = start_row; ii < end_row; ii++) {

      cs_real_t  sii = 0.0;

      m_row = mc->val + ms->row_index[ii];
      n_cols = ms->row_index[ii+1] - ms->row_index[ii];

      for (jj = 0; jj < n_cols; jj++)
        sii += *m_row++ * *prefetch_p++;

      y[ii] = sii;

    }

  }

}

/*----------------------------------------------------------------------------
 * Create a symmetric CSR matrix structure from a native matrix stucture.
 *
 * Note that the structure created maps global cell numbers to the given
 * existing face -> cell connectivity array, so it must be destroyed before
 * this array (usually the code's global cell numbering) is freed.
 *
 * parameters:
 *   have_diag   <-- Indicates if the diagonal is nonzero
 *                   (forced to true for symmetric variant)
 *   n_cells     <-- Local number of participating cells
 *   n_cells_ext <-- Local number of cells + ghost cells sharing a face
 *   n_faces     <-- Local number of faces
 *   face_cell   <-- Face -> cells connectivity
 *
 * returns:
 *   pointer to allocated CSR matrix structure.
 *----------------------------------------------------------------------------*/

static cs_matrix_struct_csr_sym_t *
_create_struct_csr_sym(bool                have_diag,
                       cs_lnum_t           n_cells,
                       cs_lnum_t           n_cells_ext,
                       cs_lnum_t           n_faces,
                       const cs_lnum_2_t  *face_cell)
{
  int n_cols_max;
  cs_lnum_t ii, jj, face_id;
  const cs_lnum_t *restrict face_cel_p;

  cs_lnum_t  diag_elts = 1;
  cs_lnum_t  *ccount = NULL;

  cs_matrix_struct_csr_sym_t  *ms;

  /* Allocate and map */

  BFT_MALLOC(ms, 1, cs_matrix_struct_csr_sym_t);

  ms->n_rows = n_cells;
  ms->n_cols = n_cells_ext;

  ms->have_diag = have_diag;
  ms->direct_assembly = true;

  BFT_MALLOC(ms->row_index, ms->n_rows + 1, cs_lnum_t);
  ms->row_index = ms->row_index;

  /* Count number of nonzero elements per row */

  BFT_MALLOC(ccount, ms->n_cols, cs_lnum_t);

  if (have_diag == false)
    diag_elts = 0;

  for (ii = 0; ii < ms->n_rows; ii++)  /* count starting with diagonal terms */
    ccount[ii] = diag_elts;

  if (face_cell != NULL) {

    face_cel_p = (const cs_lnum_t *restrict)face_cell;

    for (face_id = 0; face_id < n_faces; face_id++) {
      ii = *face_cel_p++;
      jj = *face_cel_p++;
      if (ii < jj)
        ccount[ii] += 1;
      else
        ccount[jj] += 1;
    }

  } /* if (face_cell != NULL) */

  n_cols_max = 0;

  ms->row_index[0] = 0;
  for (ii = 0; ii < ms->n_rows; ii++) {
    ms->row_index[ii+1] = ms->row_index[ii] + ccount[ii];
    if (ccount[ii] > n_cols_max)
      n_cols_max = ccount[ii];
    ccount[ii] = diag_elts; /* pre-count for diagonal terms */
  }

  ms->n_cols_max = n_cols_max;

  /* Build structure */

  BFT_MALLOC(ms->col_id, (ms->row_index[ms->n_rows]), cs_lnum_t);
  ms->col_id = ms->col_id;

  if (have_diag == true) {
    for (ii = 0; ii < ms->n_rows; ii++) {    /* diagonal terms */
      ms->col_id[ms->row_index[ii]] = ii;
    }
  }

  if (face_cell != NULL) {                   /* non-diagonal terms */

    face_cel_p = (const cs_lnum_t *restrict)face_cell;

    for (face_id = 0; face_id < n_faces; face_id++) {
      ii = *face_cel_p++;
      jj = *face_cel_p++;
      if (ii < jj && ii < ms->n_rows) {
        ms->col_id[ms->row_index[ii] + ccount[ii]] = jj;
        ccount[ii] += 1;
      }
      else if (ii > jj && jj < ms->n_rows) {
        ms->col_id[ms->row_index[jj] + ccount[jj]] = ii;
        ccount[jj] += 1;
      }
    }

  }

  BFT_FREE(ccount);

  /* Compact elements if necessary */

  if (ms->direct_assembly == false) {

    cs_lnum_t *tmp_row_index = NULL;
    cs_lnum_t  kk = 0;

    BFT_MALLOC(tmp_row_index, ms->n_rows+1, cs_lnum_t);
    memcpy(tmp_row_index, ms->row_index, (ms->n_rows+1)*sizeof(cs_lnum_t));

    kk = 0;

    for (ii = 0; ii < ms->n_rows; ii++) {
      cs_lnum_t *col_id = ms->col_id + ms->row_index[ii];
      cs_lnum_t n_cols = ms->row_index[ii+1] - ms->row_index[ii];
      cs_lnum_t col_id_prev = -1;
      ms->row_index[ii] = kk;
      for (jj = 0; jj < n_cols; jj++) {
        if (col_id_prev != col_id[jj]) {
          ms->col_id[kk++] = col_id[jj];
          col_id_prev = col_id[jj];
        }
      }
    }
    ms->row_index[ms->n_rows] = kk;

    assert(ms->row_index[ms->n_rows] < tmp_row_index[ms->n_rows]);

    BFT_FREE(tmp_row_index);
    BFT_REALLOC(ms->col_id, (ms->row_index[ms->n_rows]), cs_lnum_t);

  }

  return ms;
}

/*----------------------------------------------------------------------------
 * Destroy symmetric CSR matrix structure.
 *
 * parameters:
 *   matrix  <->  Pointer to CSR matrix structure pointer
 *----------------------------------------------------------------------------*/

static void
_destroy_struct_csr_sym(cs_matrix_struct_csr_sym_t  **matrix)
{
  if (matrix != NULL && *matrix !=NULL) {

    cs_matrix_struct_csr_sym_t  *ms = *matrix;

    if (ms->row_index != NULL)
      BFT_FREE(ms->row_index);

    if (ms->col_id != NULL)
      BFT_FREE(ms->col_id);

    BFT_FREE(ms);

    *matrix = ms;

  }
}

/*----------------------------------------------------------------------------
 * Create symmetric CSR matrix coefficients.
 *
 * returns:
 *   pointer to allocated CSR coefficients structure.
 *----------------------------------------------------------------------------*/

static cs_matrix_coeff_csr_sym_t *
_create_coeff_csr_sym(void)
{
  cs_matrix_coeff_csr_sym_t  *mc;

  /* Allocate */

  BFT_MALLOC(mc, 1, cs_matrix_coeff_csr_sym_t);

  /* Initialize */

  mc->val = NULL;

  mc->d_val = NULL;
  mc->_d_val = NULL;

  return mc;
}

/*----------------------------------------------------------------------------
 * Destroy symmetric CSR matrix coefficients.
 *
 * parameters:
 *   coeff  <->  Pointer to CSR matrix coefficients pointer
 *----------------------------------------------------------------------------*/

static void
_destroy_coeff_csr_sym(cs_matrix_coeff_csr_sym_t  **coeff)
{
  if (coeff != NULL && *coeff !=NULL) {

    cs_matrix_coeff_csr_sym_t  *mc = *coeff;

    if (mc->val != NULL)
      BFT_FREE(mc->val);

    if (mc->_d_val != NULL)
      BFT_FREE(mc->_d_val);

    BFT_FREE(*coeff);

  }
}

/*----------------------------------------------------------------------------
 * Set symmetric CSR extradiagonal matrix coefficients for the case where
 * direct assignment is possible (i.e. when there are no multiple
 * contributions to a given coefficient).
 *
 * parameters:
 *   matrix    <-- Pointer to matrix structure
 *   xa        <-- Extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_xa_coeffs_csr_sym_direct(cs_matrix_t      *matrix,
                              const cs_real_t  *restrict xa)
{
  cs_matrix_coeff_csr_sym_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_sym_t  *ms = matrix->structure;
  const cs_lnum_t n_faces = matrix->n_faces;
  const cs_lnum_2_t *restrict face_cel = matrix->face_cell;

  /* Copy extra-diagonal values */

  assert(matrix->face_cell != NULL);

# pragma omp parallel for  if(n_faces > CS_THR_MIN)
  for (cs_lnum_t face_id = 0; face_id < n_faces; face_id++) {
    cs_lnum_t ii = face_cel[face_id][0];
    cs_lnum_t jj = face_cel[face_id][1];
    if (ii < jj && ii < ms->n_rows) {
      cs_lnum_t kk;
      for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
      mc->val[kk] = xa[face_id];
    }
    else if (ii > jj && jj < ms->n_rows) {
      cs_lnum_t kk;
      for (kk = ms->row_index[jj]; ms->col_id[kk] != ii; kk++);
      mc->val[kk] = xa[face_id];
    }
  }
}

/*----------------------------------------------------------------------------
 * Set symmetric CSR extradiagonal matrix coefficients for the case where
 * there are multiple contributions to a given coefficient).
 *
 * The matrix coefficients should have been initialized (i.e. set to 0)
 * some before using this function.
 *
 * parameters:
 *   matrix    <-- Pointer to matrix structure
 *   symmetric <-- Indicates if extradiagonal values are symmetric
 *   xa        <-- Extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_xa_coeffs_csr_sym_increment(cs_matrix_t      *matrix,
                                 const cs_real_t  *restrict xa)
{
  cs_lnum_t  ii, jj, face_id;
  cs_matrix_coeff_csr_sym_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_sym_t  *ms = matrix->structure;
  const cs_lnum_t n_faces = matrix->n_faces;
  const cs_lnum_t *restrict face_cel_p
    = (const cs_lnum_t *restrict)(matrix->face_cell);

  /* Copy extra-diagonal values */

  assert(matrix->face_cell != NULL);

  for (face_id = 0; face_id < n_faces; face_id++) {
    cs_lnum_t kk;
    ii = *face_cel_p++;
    jj = *face_cel_p++;
    if (ii < jj && ii < ms->n_rows) {
      for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
      mc->val[kk] += xa[face_id];
    }
    else if (ii > jj && jj < ms->n_rows) {
      for (kk = ms->row_index[jj]; ms->col_id[kk] != ii; kk++);
      mc->val[kk] += xa[face_id];
    }
  }
}

/*----------------------------------------------------------------------------
 * Set symmetric CSR matrix coefficients.
 *
 * parameters:
 *   matrix           <-> Pointer to matrix structure
 *   symmetric        <-- Indicates if extradiagonal values are symmetric (true)
 *   interleaved      <-- Indicates if matrix coefficients are interleaved
 *   copy             <-- Indicates if coefficients should be copied
 *   da               <-- Diagonal values (NULL if all zero)
 *   xa               <-- Extradiagonal values (NULL if all zero)
 *----------------------------------------------------------------------------*/

static void
_set_coeffs_csr_sym(cs_matrix_t      *matrix,
                    bool              symmetric,
                    bool              interleaved,
                    bool              copy,
                    const cs_real_t  *restrict da,
                    const cs_real_t  *restrict xa)
{
  cs_matrix_coeff_csr_sym_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_sym_t  *ms = matrix->structure;

  if (mc->val == NULL)
    BFT_MALLOC(mc->val, ms->row_index[ms->n_rows], cs_real_t);

  /* Initialize coefficients to zero if assembly is incremental */

  if (ms->direct_assembly == false) {
    cs_lnum_t val_size = ms->row_index[ms->n_rows];
#   pragma omp parallel for  if(val_size > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < val_size; ii++)
      mc->val[ii] = 0.0;
  }

  /* Copy diagonal values */

  if (ms->have_diag == true) {

    const cs_lnum_t *_diag_index = ms->row_index;

    if (da != NULL) {
#     pragma omp parallel for  if(ms->n_rows > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < ms->n_rows; ii++)
        mc->val[_diag_index[ii]] = da[ii];
    }
    else {
#     pragma omp parallel for  if(ms->n_rows > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < ms->n_rows; ii++)
        mc->val[_diag_index[ii]] = 0.0;
    }

  }

  /* Copy extra-diagonal values */

  if (matrix->face_cell != NULL) {

    if (xa != NULL) {

      if (symmetric == false)
        bft_error(__FILE__, __LINE__, 0,
                  _("Assigning non-symmetric matrix coefficients to a matrix\n"
                    "in a symmetric CSR format."));

      if (ms->direct_assembly == true)
        _set_xa_coeffs_csr_sym_direct(matrix, xa);
      else
        _set_xa_coeffs_csr_sym_increment(matrix, xa);

    }
    else { /* if (xa == NULL) */

#     pragma omp parallel for  if(ms->n_rows > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < ms->n_rows; ii++) {
        const cs_lnum_t *restrict col_id = ms->col_id + ms->row_index[ii];
        cs_real_t *m_row = mc->val + ms->row_index[ii];
        cs_lnum_t  n_cols = ms->row_index[ii+1] - ms->row_index[ii];

        for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
          if (col_id[jj] != ii)
            m_row[jj] = 0.0;
        }

      }

    }

  } /* (matrix->face_cell != NULL) */

}

/*----------------------------------------------------------------------------
 * Release shared symmetric CSR matrix coefficients.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *----------------------------------------------------------------------------*/

static void
_release_coeffs_csr_sym(cs_matrix_t  *matrix)
{
  cs_matrix_coeff_csr_sym_t  *mc = matrix->coeffs;
  if (mc != NULL)
    mc->d_val = NULL;
  return;
}

/*----------------------------------------------------------------------------
 * Copy diagonal of symmetric CSR matrix.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *   da     --> Diagonal (pre-allocated, size: n_rows)
 *----------------------------------------------------------------------------*/

static void
_copy_diagonal_csr_sym(const cs_matrix_t  *matrix,
                       cs_real_t          *restrict da)
{
  cs_lnum_t  ii;
  const cs_matrix_struct_csr_sym_t  *ms = matrix->structure;
  const cs_matrix_coeff_csr_sym_t  *mc = matrix->coeffs;
  cs_lnum_t  n_rows = ms->n_rows;

  if (ms->have_diag == true) {

    /* As structure is symmetric, diagonal values appear first,
       so diag_index == row_index */

    const cs_lnum_t *diag_index = ms->row_index;

    for (ii = 0; ii < n_rows; ii++)
      da[ii] = mc->val[diag_index[ii]];

  }
  else { /* if (have_diag == false) */

    for (ii = 0; ii < n_rows; da[ii++] = 0.0);

  }

}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with symmetric CSR matrix.
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- Pointer to matrix structure
 *   x            <-- Multipliying vector values
 *   y            --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_p_l_csr_sym(bool                 exclude_diag,
                     const cs_matrix_t   *matrix,
                     const cs_real_t     *restrict x,
                     cs_real_t           *restrict y)
{
  cs_lnum_t  ii, jj, n_cols;
  cs_lnum_t  *restrict col_id;
  cs_real_t  *restrict m_row;

  const cs_matrix_struct_csr_sym_t  *ms = matrix->structure;
  const cs_matrix_coeff_csr_sym_t  *mc = matrix->coeffs;
  cs_lnum_t  n_rows = ms->n_rows;

  cs_lnum_t jj_start = 0;
  cs_lnum_t sym_jj_start = 0;

  /* Tell IBM compiler not to alias */
# if defined(__xlc__)
# pragma disjoint(*x, *y, *m_row, *col_id)
# endif

  /* By construction, the matrix has either a full or an empty
     diagonal structure, so testing this on the first row is enough */

  if (ms->col_id[ms->row_index[0]] == 0) {
    sym_jj_start = 1;
    if (exclude_diag)
      jj_start = 1;
  }

  /* Initialize y */

  for (ii = 0; ii < ms->n_cols; ii++)
    y[ii] = 0.0;

  /* Upper triangular + diagonal part in case of symmetric structure */

  for (ii = 0; ii < n_rows; ii++) {

    cs_real_t  sii = 0.0;

    col_id = ms->col_id + ms->row_index[ii];
    m_row = mc->val + ms->row_index[ii];
    n_cols = ms->row_index[ii+1] - ms->row_index[ii];

    for (jj = jj_start; jj < n_cols; jj++)
      sii += (m_row[jj]*x[col_id[jj]]);

    y[ii] += sii;

    for (jj = sym_jj_start; jj < n_cols; jj++)
      y[col_id[jj]] += (m_row[jj]*x[ii]);
  }

}

#if defined (HAVE_MKL)

static void
_mat_vec_p_l_csr_sym_mkl(bool                exclude_diag,
                         const cs_matrix_t  *matrix,
                         const cs_real_t    *restrict x,
                         cs_real_t          *restrict y)
{
  const cs_matrix_struct_csr_sym_t  *ms = matrix->structure;
  const cs_matrix_coeff_csr_sym_t  *mc = matrix->coeffs;

  int n_rows = ms->n_rows;
  char uplo[] = "u";

  if (exclude_diag)
    bft_error(__FILE__, __LINE__, 0,
              _(_no_exclude_diag_error_str), __func__);

  mkl_cspblas_dcsrsymv(uplo,
                       &n_rows,
                       mc->val,
                       ms->row_index,
                       ms->col_id,
                       (double *)x,
                       y);
}

#endif /* defined (HAVE_MKL) */

/*----------------------------------------------------------------------------
 * Create MSR matrix coefficients.
 *
 * returns:
 *   pointer to allocated MSR coefficients structure.
 *----------------------------------------------------------------------------*/

static cs_matrix_coeff_msr_t *
_create_coeff_msr(void)
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
 * Destroy MSR matrix coefficients.
 *
 * parameters:
 *   coeff  <->  Pointer to MSR matrix coefficients pointer
 *----------------------------------------------------------------------------*/

static void
_destroy_coeff_msr(cs_matrix_coeff_msr_t  **coeff)
{
  if (coeff != NULL && *coeff !=NULL) {

    cs_matrix_coeff_msr_t  *mc = *coeff;

    if (mc->x_prefetch != NULL)
      BFT_FREE(mc->x_prefetch);

    if (mc->x_val != NULL)
      BFT_FREE(mc->x_val);

    if (mc->_d_val != NULL)
      BFT_FREE(mc->_d_val);

    BFT_FREE(*coeff);

  }
}

/*----------------------------------------------------------------------------
 * Set MSR extradiagonal matrix coefficients for the case where direct
 * assignment is possible (i.e. when there are no multiple contributions
 * to a given coefficient).
 *
 * parameters:
 *   matrix      <-- Pointer to matrix structure
 *   symmetric   <-- Indicates if extradiagonal values are symmetric
 *   interleaved <-- Indicates if matrix coefficients are interleaved
 *   xa          <-- Extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_xa_coeffs_msr_direct(cs_matrix_t      *matrix,
                          bool              symmetric,
                          bool              interleaved,
                          const cs_real_t  *restrict xa)
{
  cs_lnum_t  ii, jj, face_id;
  cs_matrix_coeff_msr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  /* Copy extra-diagonal values */

  assert(matrix->face_cell != NULL);

  if (symmetric == false) {

    const cs_lnum_t n_faces = matrix->n_faces;
    const cs_lnum_t *restrict face_cel_p
      = (const cs_lnum_t *restrict)(matrix->face_cell);

    if (interleaved == false) {
      const cs_real_t  *restrict xa1 = xa;
      const cs_real_t  *restrict xa2 = xa + matrix->n_faces;
      for (face_id = 0; face_id < n_faces; face_id++) {
        cs_lnum_t kk, ll;
        ii = *face_cel_p++;
        jj = *face_cel_p++;
        if (ii < ms->n_rows) {
          for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
          mc->x_val[kk] = xa1[face_id];
        }
        if (jj < ms->n_rows) {
          for (ll = ms->row_index[jj]; ms->col_id[ll] != ii; ll++);
          mc->x_val[ll] = xa2[face_id];
        }
      }
    }
    else { /* interleaved == true */
      for (face_id = 0; face_id < n_faces; face_id++) {
        cs_lnum_t kk, ll;
        ii = *face_cel_p++;
        jj = *face_cel_p++;
        if (ii < ms->n_rows) {
          for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
          mc->x_val[kk] = xa[2*face_id];
        }
        if (jj < ms->n_rows) {
          for (ll = ms->row_index[jj]; ms->col_id[ll] != ii; ll++);
          mc->x_val[ll] = xa[2*face_id + 1];
        }
      }
    }

  }
  else { /* if symmetric == true */

    const cs_lnum_t n_faces = matrix->n_faces;
    const cs_lnum_t *restrict face_cel_p
      = (const cs_lnum_t *restrict)(matrix->face_cell);

    for (face_id = 0; face_id < n_faces; face_id++) {
      cs_lnum_t kk, ll;
      ii = *face_cel_p++;
      jj = *face_cel_p++;
      if (ii < ms->n_rows) {
        for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
        mc->x_val[kk] = xa[face_id];
      }
      if (jj < ms->n_rows) {
        for (ll = ms->row_index[jj]; ms->col_id[ll] != ii; ll++);
        mc->x_val[ll] = xa[face_id];
      }

    }

  } /* end of condition on coefficients symmetry */

}

/*----------------------------------------------------------------------------
 * Set MSR extradiagonal matrix coefficients for the case where there are
 * multiple contributions to a given coefficient).
 *
 * The matrix coefficients should have been initialized (i.e. set to 0)
 * some before using this function.
 *
 * parameters:
 *   matrix      <-- Pointer to matrix structure
 *   symmetric   <-- Indicates if extradiagonal values are symmetric
 *   interleaved <-- Indicates if matrix coefficients are interleaved
 *   xa          <-- Extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_xa_coeffs_msr_increment(cs_matrix_t      *matrix,
                             bool              symmetric,
                             bool              interleaved,
                             const cs_real_t  *restrict xa)
{
  cs_lnum_t  ii, jj, face_id;
  cs_matrix_coeff_msr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  /* Copy extra-diagonal values */

  assert(matrix->face_cell != NULL);

  if (symmetric == false) {

    const cs_lnum_t n_faces = matrix->n_faces;
    const cs_lnum_t *restrict face_cel_p
      = (const cs_lnum_t *restrict)(matrix->face_cell);

    const cs_real_t  *restrict xa1 = xa;
    const cs_real_t  *restrict xa2 = xa + matrix->n_faces;

    if (interleaved == false) {
      for (face_id = 0; face_id < n_faces; face_id++) {
        cs_lnum_t kk, ll;
        ii = *face_cel_p++;
        jj = *face_cel_p++;
        if (ii < ms->n_rows) {
          for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
          mc->x_val[kk] += xa1[face_id];
        }
        if (jj < ms->n_rows) {
          for (ll = ms->row_index[jj]; ms->col_id[ll] != ii; ll++);
          mc->x_val[ll] += xa2[face_id];
        }
      }
    }
    else { /* interleaved == true */
      for (face_id = 0; face_id < n_faces; face_id++) {
        cs_lnum_t kk, ll;
        ii = *face_cel_p++;
        jj = *face_cel_p++;
        if (ii < ms->n_rows) {
          for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
          mc->x_val[kk] += xa[2*face_id];
        }
        if (jj < ms->n_rows) {
          for (ll = ms->row_index[jj]; ms->col_id[ll] != ii; ll++);
          mc->x_val[ll] += xa[2*face_id + 1];
        }
      }
    }
  }
  else { /* if symmetric == true */

    const cs_lnum_t n_faces = matrix->n_faces;
    const cs_lnum_t *restrict face_cel_p
      = (const cs_lnum_t *restrict)(matrix->face_cell);

    for (face_id = 0; face_id < n_faces; face_id++) {
      cs_lnum_t kk, ll;
      ii = *face_cel_p++;
      jj = *face_cel_p++;
      if (ii < ms->n_rows) {
        for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
        mc->x_val[kk] += xa[face_id];
      }
      if (jj < ms->n_rows) {
        for (ll = ms->row_index[jj]; ms->col_id[ll] != ii; ll++);
        mc->x_val[ll] += xa[face_id];
      }

    }

  } /* end of condition on coefficients symmetry */

}

/*----------------------------------------------------------------------------
 * Set MSR matrix coefficients.
 *
 * parameters:
 *   matrix           <-> Pointer to matrix structure
 *   symmetric        <-- Indicates if extradiagonal values are symmetric
 *   interleaved      <-- Indicates if matrix coefficients are interleaved
 *   copy             <-- Indicates if coefficients should be copied
 *   da               <-- Diagonal values (NULL if all zero)
 *   xa               <-- Extradiagonal values (NULL if all zero)
 *----------------------------------------------------------------------------*/

static void
_set_coeffs_msr(cs_matrix_t      *matrix,
                bool              symmetric,
                bool              interleaved,
                bool              copy,
                const cs_real_t  *restrict da,
                const cs_real_t  *restrict xa)
{
  cs_lnum_t  ii, jj;
  cs_matrix_coeff_msr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  /* Allocate prefetch buffer if needed */

  mc->n_prefetch_rows
    =  CS_MAX(matrix->loop_length[matrix->fill_type][0],
              matrix->loop_length[matrix->fill_type][1]);
  if (mc->n_prefetch_rows > 0 && mc->x_prefetch == NULL) {
    size_t prefetch_size = ms->n_cols_max * mc->n_prefetch_rows;
    size_t matrix_size = matrix->n_cells + (2 * matrix->n_faces);
    if (matrix_size > prefetch_size)
      prefetch_size = matrix_size;
    BFT_REALLOC(mc->x_prefetch, prefetch_size, cs_real_t);
  }

  /* Map or copy diagonal values */

  if (da != NULL) {

    if (copy) {
      if (mc->_d_val == NULL || mc->max_db_size < matrix->db_size[3]) {
        BFT_REALLOC(mc->_d_val, matrix->db_size[3]*ms->n_rows, cs_real_t);
        mc->max_db_size = matrix->db_size[3];
      }
      memcpy(mc->_d_val, da, matrix->db_size[3]*sizeof(cs_real_t) * ms->n_rows);
      mc->d_val = mc->_d_val;
    }
    else
      mc->d_val = da;

  }
  else
    mc->d_val = NULL;

  /* Extradiagonal values */ //TODO with matrix->eb_size[3] > 1

  if (mc->x_val == NULL)
    BFT_MALLOC(mc->x_val, ms->row_index[ms->n_rows], cs_real_t);

  /* Initialize coefficients to zero if assembly is incremental */

  if (ms->direct_assembly == false) {
    cs_lnum_t val_size = ms->row_index[ms->n_rows];
    for (ii = 0; ii < val_size; ii++)

      mc->x_val[ii] = 0.0;
  }

  /* Copy extra-diagonal values */

  if (matrix->face_cell != NULL) {

    if (xa != NULL) {

      if (ms->direct_assembly == true)
        _set_xa_coeffs_msr_direct(matrix, symmetric, interleaved, xa);
      else
        _set_xa_coeffs_msr_increment(matrix, symmetric, interleaved, xa);

    }
    else { /* if (xa == NULL) */

      for (ii = 0; ii < ms->n_rows; ii++) {
        const cs_lnum_t  *restrict col_id = ms->col_id + ms->row_index[ii];
        cs_real_t  *m_row = mc->x_val + ms->row_index[ii];
        cs_lnum_t  n_cols = ms->row_index[ii+1] - ms->row_index[ii];

        for (jj = 0; jj < n_cols; jj++) {
          if (col_id[jj] != ii)
            m_row[jj] = 0.0;
        }

      }

    }

  } /* (matrix->face_cell != NULL) */

}

/*----------------------------------------------------------------------------
 * Release shared MSR matrix coefficients.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *----------------------------------------------------------------------------*/

static void
_release_coeffs_msr(cs_matrix_t  *matrix)
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

static void
_mat_vec_p_l_msr(bool                exclude_diag,
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

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        if (col_id[jj] != ii)
          sii += (m_row[jj]*x[col_id[jj]]);
      }

      y[ii] = sii;

    }
  }

}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with MSR matrix, blocked version.
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- Pointer to matrix structure
 *   x            <-- Multipliying vector values
 *   y            --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_b_mat_vec_p_l_msr(bool                exclude_diag,
                   const cs_matrix_t  *matrix,
                   const cs_real_t    *restrict x,
                   cs_real_t          *restrict y)
{
  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_msr_t  *mc = matrix->coeffs;
  const cs_lnum_t  n_rows = ms->n_rows;
  const int *db_size = matrix->db_size;

  /* Standard case */

  if (!exclude_diag && mc->d_val != NULL) {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      cs_lnum_t *restrict col_id = ms->col_id + ms->row_index[ii];
      cs_real_t *restrict m_row = mc->x_val + ms->row_index[ii];
      cs_lnum_t n_cols = ms->row_index[ii+1] - ms->row_index[ii];

      _dense_b_ax(ii, db_size, mc->d_val, x, y);

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        for (cs_lnum_t kk = 0; kk < db_size[0]; kk++) {
          y[ii*db_size[1] + kk]
            += (m_row[jj]*x[col_id[jj]*db_size[1] + kk]);
        }
      }

    }

  }

  /* Exclude diagonal */

  else {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      cs_lnum_t *restrict col_id = ms->col_id + ms->row_index[ii];
      cs_real_t *restrict m_row = mc->x_val + ms->row_index[ii];
      cs_lnum_t n_cols = ms->row_index[ii+1] - ms->row_index[ii];

      for (cs_lnum_t kk = 0; kk < db_size[0]; kk++)
        y[ii*db_size[1] + kk] = 0.;

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        for (cs_lnum_t kk = 0; kk < db_size[0]; kk++) {
          y[ii*db_size[1] + kk]
            += (m_row[jj]*x[col_id[jj]*db_size[1] + kk]);
        }
      }

    }
  }

}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with MSR matrix, using MKL
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- Pointer to matrix structure
 *   x            <-- Multipliying vector values
 *   y            --> Resulting vector
 *----------------------------------------------------------------------------*/

#if defined (HAVE_MKL)

static void
_mat_vec_p_l_msr_mkl(bool                exclude_diag,
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

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with MSR matrix (prefetch).
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- Pointer to matrix structure
 *   x            <-- Multipliying vector values
 *   y            --> Resulting vector
 *----------------------------------------------------------------------------*/

#if defined(__INTEL_COMPILER)
#pragma optimization_level 0 /* With icc 11.1.072, errors occur above this */
#endif

static void
_mat_vec_p_l_msr_pf(bool                exclude_diag,
                    const cs_matrix_t  *matrix,
                    const cs_real_t    *restrict x,
                    cs_real_t          *restrict y)
{
  cs_lnum_t  start_row, ii, jj, n_cols;
  cs_lnum_t  *restrict col_id;
  cs_real_t  *restrict m_row;
  cs_real_t  *restrict prefetch_p;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_msr_t  *mc = matrix->coeffs;
  cs_lnum_t  n_rows = ms->n_rows;

  /* Standard case */

  if (!exclude_diag && mc->d_val != NULL) {

    /* Outer loop on prefetch lines */

    for (start_row = 0; start_row < n_rows; start_row += mc->n_prefetch_rows) {

      cs_lnum_t end_row = start_row + mc->n_prefetch_rows;

      prefetch_p = mc->x_prefetch;

      if (end_row > n_rows)
        end_row = n_rows;

      /* Prefetch */

      for (ii = start_row; ii < end_row; ii++) {

        col_id = ms->col_id + ms->row_index[ii];
        n_cols = ms->row_index[ii+1] - ms->row_index[ii];

        for (jj = 0; jj < n_cols; jj++)
          *prefetch_p++ = x[col_id[jj]];

        *prefetch_p++ = x[ii];

      }

      /* Compute */

      prefetch_p = mc->x_prefetch;

      for (ii = start_row; ii < end_row; ii++) {

        cs_real_t  sii = 0.0;

        m_row = mc->x_val + ms->row_index[ii];
        n_cols = ms->row_index[ii+1] - ms->row_index[ii];

        for (jj = 0; jj < n_cols; jj++)
          sii += *m_row++ * *prefetch_p++;

        y[ii] = sii + (mc->d_val[ii] * *prefetch_p++);

      }

    }

  }

  /* Exclude diagonal */

  else {

    /* Outer loop on prefetch lines */

    for (start_row = 0; start_row < n_rows; start_row += mc->n_prefetch_rows) {

      cs_lnum_t end_row = start_row + mc->n_prefetch_rows;

      prefetch_p = mc->x_prefetch;

      if (end_row > n_rows)
        end_row = n_rows;

      /* Prefetch */

      for (ii = start_row; ii < end_row; ii++) {

        col_id = ms->col_id + ms->row_index[ii];
        n_cols = ms->row_index[ii+1] - ms->row_index[ii];

        for (jj = 0; jj < n_cols; jj++)
          *prefetch_p++ = x[col_id[jj]];

      }

      /* Compute */

      prefetch_p = mc->x_prefetch;

      for (ii = start_row; ii < end_row; ii++) {

        cs_real_t  sii = 0.0;

        m_row = mc->x_val + ms->row_index[ii];
        n_cols = ms->row_index[ii+1] - ms->row_index[ii];

        for (jj = 0; jj < n_cols; jj++)
          sii += *m_row++ * *prefetch_p++;

        y[ii] = sii;

      }
    }

  }
}

/*----------------------------------------------------------------------------
 * Synchronize ghost cells prior to matrix.vector product
 *
 * parameters:
 *   rotation_mode <-- Halo update option for rotational periodicity
 *   matrix        <-- Pointer to matrix structure
 *   x             <-> Multipliying vector values (ghost values updated)
 *   y             --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_pre_vector_multiply_sync(cs_halo_rotation_t   rotation_mode,
                          const cs_matrix_t   *matrix,
                          cs_real_t           *restrict x,
                          cs_real_t           *restrict y)
{
  size_t n_cells_ext = matrix->n_cells_ext;

  assert(matrix->halo != NULL);

  /* Non-blocked version */

  if (matrix->db_size[3] == 1) {

    /* Synchronize for parallelism and periodicity first */

    _zero_range(y, matrix->n_cells, n_cells_ext);

    /* Update distant ghost cells */

    if (matrix->halo != NULL)
      cs_halo_sync_component(matrix->halo,
                             CS_HALO_STANDARD,
                             rotation_mode,
                             x);

  }

  /* Blocked version */

  else { /* if (matrix->db_size[3] > 1) */

    const int *db_size = matrix->db_size;

    /* Synchronize for parallelism and periodicity first */

    _b_zero_range(y, matrix->n_cells, n_cells_ext, db_size);

    /* Update distant ghost cells */

    if (matrix->halo != NULL) {

      cs_halo_sync_var_strided(matrix->halo,
                               CS_HALO_STANDARD,
                               x,
                               db_size[1]);

      /* Synchronize periodic values */

      if (matrix->halo->n_transforms > 0 && db_size[0] == 3)
        cs_halo_perio_sync_var_vect(matrix->halo,
                                    CS_HALO_STANDARD,
                                    x,
                                    db_size[1]);

    }

  }
}

/*----------------------------------------------------------------------------
 * Copy array to reference for matrix computation check.
 *
 * parameters:
 *   n_elts      <-- number values to compare
 *   y           <-- array to copare or copy
 *   yr          <-- reference array
 *
 * returns:
 *   maximum difference between values
 *----------------------------------------------------------------------------*/

static double
_matrix_check_compare(cs_lnum_t        n_elts,
                      const cs_real_t  y[],
                      cs_real_t        yr[])
{
  cs_lnum_t  ii;

  double dmax = 0.0;

  for (ii = 0; ii < n_elts; ii++) {
    double d = CS_ABS(y[ii] - yr[ii]);
    if (d > dmax)
      dmax = d;
  }

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {
    double dmaxg;
    MPI_Allreduce(&dmax, &dmaxg, 1, MPI_DOUBLE, MPI_MAX, cs_glob_mpi_comm);
    dmax = dmaxg;
  }

#endif

  return dmax;
}

/*----------------------------------------------------------------------------
 * Check local matrix.vector product operations.
 *
 * parameters:
 *   t_measure   <-- minimum time for each measure
 *   n_variants  <-- number of variants in array
 *   n_cells     <-- number of local cells
 *   n_cells_ext <-- number of cells including ghost cells (array size)
 *   n_faces     <-- local number of internal faces
 *   cell_num    <-- Optional global cell numbers (1 to n), or NULL
 *   face_cell   <-- face -> cells connectivity (1 to n)
 *   halo        <-- cell halo structure
 *   numbering   <-- vectorization or thread-related numbering info, or NULL
 *   m_variant   <-> array of matrix variants
 *----------------------------------------------------------------------------*/

static void
_matrix_check(int                    n_variants,
              cs_lnum_t              n_cells,
              cs_lnum_t              n_cells_ext,
              cs_lnum_t              n_faces,
              const cs_gnum_t       *cell_num,
              const cs_lnum_2_t     *face_cell,
              const cs_halo_t       *halo,
              const cs_numbering_t  *numbering,
              cs_matrix_variant_t   *m_variant)
{
  cs_lnum_t  ii;
  int  v_id, f_id, ed_flag;

  bool print_subtitle = false;
  cs_real_t  *da = NULL, *xa = NULL, *x = NULL, *y = NULL;
  cs_real_t  *yr0 = NULL, *yr1 = NULL;
  cs_matrix_structure_t *ms = NULL;
  cs_matrix_t *m = NULL;
  int d_block_size[4] = {3, 3, 3, 9};
  int ed_block_size[4] = {3, 3, 3, 9};

  /* Allocate and initialize  working arrays */

  if (CS_MEM_ALIGN > 0) {
    BFT_MEMALIGN(x, CS_MEM_ALIGN, n_cells_ext*d_block_size[1], cs_real_t);
    BFT_MEMALIGN(y, CS_MEM_ALIGN, n_cells_ext*d_block_size[1], cs_real_t);
    BFT_MEMALIGN(yr0, CS_MEM_ALIGN, n_cells_ext*d_block_size[1], cs_real_t);
    BFT_MEMALIGN(yr1, CS_MEM_ALIGN, n_cells_ext*d_block_size[1], cs_real_t);
  }
  else {
    BFT_MALLOC(x, n_cells_ext*d_block_size[1], cs_real_t);
    BFT_MALLOC(y, n_cells_ext*d_block_size[1], cs_real_t);
    BFT_MALLOC(yr0, n_cells_ext*d_block_size[1], cs_real_t);
    BFT_MALLOC(yr1, n_cells_ext*d_block_size[1], cs_real_t);
  }

  BFT_MALLOC(da, n_cells_ext*d_block_size[3], cs_real_t);
  BFT_MALLOC(xa, n_faces*2*ed_block_size[3], cs_real_t);

  /* Initialize arrays */

# pragma omp parallel for
  for (ii = 0; ii < n_cells_ext*d_block_size[3]; ii++)
    da[ii] = 1.0 + cos(ii);

# pragma omp parallel for
  for (ii = 0; ii < n_faces*ed_block_size[3]; ii++) {
    xa[ii*2] = 0.5*(0.9 + cos(ii));
    xa[ii*2 + 1] = -0.5*(0.9 + cos(ii));
  }

# pragma omp parallel for
  for (ii = 0; ii < n_cells_ext*d_block_size[1]; ii++)
    x[ii] = sin(ii);

  /* Loop on fill options */

  for (f_id = 0; f_id < CS_MATRIX_N_FILL_TYPES; f_id++) {

    const int *_d_block_size
      = (f_id >= CS_MATRIX_33_BLOCK_D) ? d_block_size : NULL;
    const int *_ed_block_size
      = (f_id >= CS_MATRIX_33_BLOCK) ? ed_block_size : NULL;
    const cs_lnum_t _block_mult = (_d_block_size != NULL) ? d_block_size[1] : 1;
    const bool sym_coeffs = (   f_id == CS_MATRIX_SCALAR_SYM
                             || f_id == CS_MATRIX_33_BLOCK_D_SYM) ? true : false;

    /* Loop on diagonal exclusion options */

    for (ed_flag = 0; ed_flag < 2; ed_flag++) {

      print_subtitle = true;

      /* Loop on variant types */

      for (v_id = 0; v_id < n_variants; v_id++) {

        cs_matrix_variant_t *v = m_variant + v_id;

        cs_matrix_vector_product_t  *vector_multiply
          = v->vector_multiply[f_id][ed_flag];

        if (vector_multiply == NULL)
          continue;

        ms = cs_matrix_structure_create(v->type,
                                        true,
                                        n_cells,
                                        n_cells_ext,
                                        n_faces,
                                        cell_num,
                                        face_cell,
                                        halo,
                                        numbering);
        m = cs_matrix_create(ms);

        m->loop_length[f_id][ed_flag] = v->loop_length[f_id][ed_flag];

        cs_matrix_set_coefficients(m,
                                   sym_coeffs,
                                   _d_block_size,
                                   _ed_block_size,
                                   da,
                                   xa);

        /* Check multiplication */

        vector_multiply(ed_flag, m, x, y);
        if (v_id == 0)
          memcpy(yr0, y, n_cells*_block_mult*sizeof(cs_real_t));
        else {
          double dmax = _matrix_check_compare(n_cells*_block_mult, y, yr0);
          if (print_subtitle) {
            bft_printf("\n%s\n",
                       _matrix_operation_name[f_id*2 + ed_flag]);
            print_subtitle = false;
          }
          bft_printf("  %-32s : %12.5e\n",
                     v->name,
                     dmax);
          bft_printf_flush();
        }

        cs_matrix_release_coefficients(m);
        cs_matrix_destroy(&m);
        cs_matrix_structure_destroy(&ms);

      } /* end of loop on variants */

    } /* end of loop on ed_flag */

  } /* end of loop on fill types */

  BFT_FREE(yr1);
  BFT_FREE(yr0);

  BFT_FREE(y);
  BFT_FREE(x);

  BFT_FREE(xa);
  BFT_FREE(da);
}

/*----------------------------------------------------------------------------
 * Initialize local variant matrix.
 *
 * parameters:
 *   v  <-> pointer to matrix variant
 *----------------------------------------------------------------------------*/

static void
_variant_init(cs_matrix_variant_t  *v)
{
  v->matrix_create_cost = -1.;
  for (int i = 0; i < CS_MATRIX_N_FILL_TYPES; i++) {
    for (int j = 0; j < 2; j++) {
      v->vector_multiply[i][j] = NULL;
      v->loop_length[i][j] = 0;
      v->matrix_vector_cost[i][j] = -1.;
    }
    v->matrix_assign_cost[i] = -1.;
  }
}

/*----------------------------------------------------------------------------
 * Add variant
 *
 * parameters:
 *   name                 <-- matrix variant name
 *   type                 <-- matrix type
 *   n_fill_types         <-- number of fill types tuned for
 *   fill_types           <-- array of fill types tuned for
 *   ed_flag              <-- 0: with diagonal only, 1 exclude only; 2; both
 *   loop_length          <-- loop length option for some algorithms
 *   vector_multiply      <-- function pointer for A.x
 *   b_vector_multiply    <-- function pointer for block A.x
 *   bb_vector_multiply   <-- function pointer for block A.x
 *                             with block extra diag
 *   n_variants           <-> number of variants
 *   n_variants_max       <-> current maximum number of variants
 *   m_variant            <-> array of matrix variants
 *----------------------------------------------------------------------------*/

static void
_variant_add(const char                        *name,
             cs_matrix_type_t                   type,
             int                                n_fill_types,
             cs_matrix_fill_type_t              fill_types[],
             int                                ed_flag,
             int                                loop_length,
             cs_matrix_vector_product_t        *vector_multiply,
             cs_matrix_vector_product_t        *b_vector_multiply,
             cs_matrix_vector_product_t        *bb_vector_multiply,
             int                               *n_variants,
             int                               *n_variants_max,
             cs_matrix_variant_t              **m_variant)
{
  cs_matrix_variant_t  *v;
  int j;
  int i = *n_variants;

  if (*n_variants_max == *n_variants) {
    if (*n_variants_max == 0)
      *n_variants_max = 8;
    else
      *n_variants_max *= 2;
    BFT_REALLOC(*m_variant, *n_variants_max, cs_matrix_variant_t);
  }

  v = (*m_variant) + i;

  _variant_init(v);

  strcpy(v->name, name);
  v->type = type;

  for (j = 0; j < n_fill_types; j++) {

    cs_matrix_fill_type_t mft =  fill_types[j];

    v->loop_length[mft][0] = loop_length;
    v->loop_length[mft][1] = loop_length;

    switch(mft) {

    case CS_MATRIX_SCALAR:
    case  CS_MATRIX_SCALAR_SYM:
      if (ed_flag != 1)
        v->vector_multiply[mft][0] = vector_multiply;
      if (ed_flag != 0)
        v->vector_multiply[mft][1] = vector_multiply;
      break;

    case CS_MATRIX_33_BLOCK_D:
    case CS_MATRIX_33_BLOCK_D_SYM:
      if (ed_flag != 1)
        v->vector_multiply[mft][0] = b_vector_multiply;
      if (ed_flag != 0)
        v->vector_multiply[mft][1] = b_vector_multiply;
      break;

    case CS_MATRIX_33_BLOCK:
      if (ed_flag != 1)
        v->vector_multiply[mft][0] = bb_vector_multiply;
      if (ed_flag != 0)
        v->vector_multiply[mft][1] = bb_vector_multiply;
      break;

    default:
      assert(0);
      break;
    }

  }

  *n_variants += 1;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a matrix Structure.
 *
 * Note that the structure created maps to the given existing
 * cell global number, face -> cell connectivity arrays, and cell halo
 * structure, so it must be destroyed before they are freed
 * (usually along with the code's main face -> cell structure).
 *
 * Note that the resulting matrix structure will contain either a full or
 * an empty main diagonal, and that the extra-diagonal structure is always
 * symmetric (though the coefficients my not be, and we may choose a
 * matrix format that does not exploit ths symmetry). If the face_cell
 * connectivity argument is NULL, the matrix will be purely diagonal.
 *
 * parameters:
 *   type        <-- Type of matrix considered
 *   have_diag   <-- Indicates if the diagonal structure contains nonzeroes
 *   n_cells     <-- Local number of cells
 *   n_cells_ext <-- Local number of cells + ghost cells sharing a face
 *   n_faces     <-- Local number of internal faces
 *   cell_num    <-- Optional global cell numbers (1 to n), or NULL
 *   face_cell   <-- Face -> cells connectivity (1 to n)
 *   halo        <-- Halo structure associated with cells, or NULL
 *   numbering   <-- vectorization or thread-related numbering info, or NULL
 *
 * returns:
 *   pointer to created matrix structure;
 *----------------------------------------------------------------------------*/

cs_matrix_structure_t *
cs_matrix_structure_create(cs_matrix_type_t       type,
                           bool                   have_diag,
                           cs_lnum_t              n_cells,
                           cs_lnum_t              n_cells_ext,
                           cs_lnum_t              n_faces,
                           const cs_gnum_t       *cell_num,
                           const cs_lnum_2_t     *face_cell,
                           const cs_halo_t       *halo,
                           const cs_numbering_t  *numbering)
{
  cs_matrix_structure_t *ms;

  BFT_MALLOC(ms, 1, cs_matrix_structure_t);

  ms->type = type;

  ms->n_cells = n_cells;
  ms->n_cells_ext = n_cells_ext;
  ms->n_faces = n_faces;

  /* Define Structure */

  switch(ms->type) {
  case CS_MATRIX_NATIVE:
    ms->structure = _create_struct_native(n_cells,
                                          n_cells_ext,
                                          n_faces,
                                          face_cell);
    break;
  case CS_MATRIX_CSR:
    ms->structure = _create_struct_csr(have_diag,
                                       n_cells,
                                       n_cells_ext,
                                       n_faces,
                                       face_cell);
    break;
  case CS_MATRIX_CSR_SYM:
    ms->structure = _create_struct_csr_sym(have_diag,
                                           n_cells,
                                           n_cells_ext,
                                           n_faces,
                                           face_cell);
    break;
  case CS_MATRIX_MSR:
    ms->structure = _create_struct_csr(false,
                                       n_cells,
                                       n_cells_ext,
                                       n_faces,
                                       face_cell);
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Handling of matrixes in %s format\n"
                "is not operational yet."),
              _(cs_matrix_type_name[type]));
    break;
  }

  /* Set pointers to structures shared from mesh here */

  ms->face_cell = face_cell;
  ms->cell_num = cell_num;
  ms->halo = halo;
  ms->numbering = numbering;

  return ms;
}

/*----------------------------------------------------------------------------
 * Destroy a matrix structure.
 *
 * parameters:
 *   ms <-> Pointer to matrix structure pointer
 *----------------------------------------------------------------------------*/

void
cs_matrix_structure_destroy(cs_matrix_structure_t  **ms)
{
  if (ms != NULL && *ms != NULL) {

    cs_matrix_structure_t *_ms = *ms;

    switch(_ms->type) {
    case CS_MATRIX_NATIVE:
      {
        cs_matrix_struct_native_t *structure = _ms->structure;
        _destroy_struct_native(&structure);
      }
      break;
    case CS_MATRIX_CSR:
      {
        cs_matrix_struct_csr_t *structure = _ms->structure;
        _destroy_struct_csr(&structure);
      }
      break;
    case CS_MATRIX_CSR_SYM:
      {
        cs_matrix_struct_csr_sym_t *structure = _ms->structure;
        _destroy_struct_csr_sym(&structure);
      }
      break;
    case CS_MATRIX_MSR:
      {
        cs_matrix_struct_csr_t *structure = _ms->structure;
        _destroy_struct_csr(&structure);
      }
      break;
    default:
      assert(0);
      break;
    }
    _ms->structure = NULL;

    /* Now free main structure */

    BFT_FREE(*ms);
  }
}

/*----------------------------------------------------------------------------
 * Select the sparse matrix-vector product function to be used by a
 * matrix or variant for a given fill type.
 *
 * Currently, possible variant functions are:
 *
 *   CS_MATRIX_NATIVE  (all fill types)
 *     default
 *     standard
 *     3_3_diag        (for CS_MATRIX_33_BLOCK_D or CS_MATRIX_33_BLOCK_D_SYM)
 *     bull            (for CS_MATRIX_SCALAR or CS_MATRIX_SCALAR_SYM)
 *     omp             (for OpenMP with compatible numbering)
 *     vector          (For vector machine with compatible numbering)
 *
 *   CS_MATRIX_CSR     (for CS_MATRIX_SCALAR or CS_MATRIX_SCALAR_SYM)
 *     default
 *     standard
 *     prefetch
 *     mkl             (with MKL)
 *
 *   CS_MATRIX_CSR_SYM (for CS_MATRIX_SCALAR_SYM)
 *     default
 *     standard
 *     mkl             (with MKL)
 *
 *   CS_MATRIX_MSR     (all fill types except CS_MATRIX_33_BLOCK)
 *     standard
 *     prefetch
 *     mkl             (with MKL, for CS_MATRIX_SCALAR or CS_MATRIX_SCALAR_SYM)
 *
 * parameters:
 *   m_type          <-- Matrix type
 *   numbering       <-- mesh numbering type, or NULL
 *   fill type       <-- matrix fill type to merge from
 *   ed_flag         <-- 0: with diagonal only, 1 exclude only; 2; both
 *   func_name       <-- function type name, or NULL for default
 *   vector_multiply <-> multiplication function array
 *   loop_length     <-> loop length array
 *
 * returns:
 *   0 for success, 1 for incompatible function, 2 for compatible
 *   function not available in current build
 *----------------------------------------------------------------------------*/

static int
_set_spmv_func(cs_matrix_type_t             m_type,
               const cs_numbering_t        *numbering,
               cs_matrix_fill_type_t        fill_type,
               int                          ed_flag,
               const char                  *func_name,
               cs_matrix_vector_product_t  *vector_multiply[][2],
               int                          loop_length[][2])
{
  int retcode = 1;
  int standard = 0;

  cs_matrix_vector_product_t *spmv[2] = {NULL, NULL};
  int l_length[2] = {0, 0};

  if (func_name == NULL)
    standard = 2;
  else if (!strcmp(func_name, "default"))
    standard = 2;
  else if (!strcmp(func_name, "standard"))
    standard = 1;

  switch(m_type) {

  case CS_MATRIX_NATIVE:

    if (standard > 0) { /* standard or default */

      switch(fill_type) {
      case CS_MATRIX_SCALAR:
      case CS_MATRIX_SCALAR_SYM:
        spmv[0] = _mat_vec_p_l_native;
        spmv[1] = _mat_vec_p_l_native;
        break;
      case CS_MATRIX_33_BLOCK_D:
      case CS_MATRIX_33_BLOCK_D_SYM:
        spmv[0] = _b_mat_vec_p_l_native;
        spmv[1] = _b_mat_vec_p_l_native;
        break;
      case CS_MATRIX_33_BLOCK:
        spmv[0] = _bb_mat_vec_p_l_native;
        spmv[1] = _bb_mat_vec_p_l_native;
        break;
      default:
        break;
      }

      if (standard > 1) { /* default optimized variants */
        switch(fill_type) {
        case CS_MATRIX_SCALAR:
        case CS_MATRIX_SCALAR_SYM:
#if defined(IA64_OPTIM)
          spmv[0] = _mat_vec_p_l_native_bull;
          spmv[1] = _mat_vec_p_l_native_bull;
#endif
          if (numbering != NULL) {
#if defined(HAVE_OPENMP)
            if (numbering->type == CS_NUMBERING_THREADS) {
              spmv[0] = _mat_vec_p_l_native_omp;
              spmv[1] = _mat_vec_p_l_native_omp;
            }
#endif
#if defined(SX) && defined(_SX) /* For vector machines */
            if (numbering->type == CS_NUMBERING_VECTORIZE) {
              spmv[0] = _mat_vec_p_l_native_vector;
              spmv[1] = _mat_vec_p_l_native_vector;
            }
#endif
          }
          break;
        case CS_MATRIX_33_BLOCK_D:
        case CS_MATRIX_33_BLOCK_D_SYM:
          if (numbering != NULL) {
#if defined(HAVE_OPENMP)
            if (numbering->type == CS_NUMBERING_THREADS) {
              spmv[0] = _b_mat_vec_p_l_native_omp;
              spmv[1] = _b_mat_vec_p_l_native_omp;
            }
#endif
          }
          break;
        default:
          break;
        }
      }

    }

    else if (!strcmp(func_name, "3_3_diag")) {
      switch(fill_type) {
      case CS_MATRIX_33_BLOCK_D:
      case CS_MATRIX_33_BLOCK_D_SYM:
        spmv[0] = _3_3_mat_vec_p_l_native;
        spmv[1] = _3_3_mat_vec_p_l_native;
        break;
      default:
        break;
      }
    }

    else if (!strcmp(func_name, "bull")) {
      switch(fill_type) {
      case CS_MATRIX_SCALAR:
      case CS_MATRIX_SCALAR_SYM:
        l_length[0] = 508;
        l_length[1] = 508;
        spmv[0] = _mat_vec_p_l_native_bull;
        spmv[1] = _mat_vec_p_l_native_bull;
        break;
      default:
        break;
      }
    }

    else if (!strcmp(func_name, "omp")) {
#if defined(HAVE_OPENMP)
      if (numbering != NULL) {
        if (numbering->type == CS_NUMBERING_THREADS) {
          switch(fill_type) {
          case CS_MATRIX_SCALAR:
          case CS_MATRIX_SCALAR_SYM:
            spmv[0] = _mat_vec_p_l_native_omp;
            spmv[1] = _mat_vec_p_l_native_omp;
            break;
          case CS_MATRIX_33_BLOCK_D:
          case CS_MATRIX_33_BLOCK_D_SYM:
            spmv[0] = _b_mat_vec_p_l_native_omp;
            spmv[1] = _b_mat_vec_p_l_native_omp;
            break;
          default:
            break;
          }
        }
      }
#else
      retcode = 2;
#endif
    }

    else if (!strcmp(func_name, "vector")) {
#if defined(SX) && defined(_SX)
      switch(fill_type) {
      case CS_MATRIX_SCALAR:
      case CS_MATRIX_SCALAR_SYM:
        spmv[0] = _mat_vec_p_l_native_vector;
        spmv[1] = _mat_vec_p_l_native_vector;
        break;
      default:
        break;
      }
#else
      retcode = 2;
#endif
    }

    break;

  case CS_MATRIX_CSR:

    switch(fill_type) {
    case CS_MATRIX_SCALAR:
    case CS_MATRIX_SCALAR_SYM:
      if (standard > 0) {
        spmv[0] = _mat_vec_p_l_csr;
        spmv[1] = _mat_vec_p_l_csr;
      }
      else if (!strcmp(func_name, "prefetch")) {
        l_length[0] = 508;
        l_length[1] = 508;
        spmv[0] = _mat_vec_p_l_csr_pf;
        spmv[1] = _mat_vec_p_l_csr_pf;
      }
      else if (!strcmp(func_name, "mkl")) {
#if defined(HAVE_MKL)
        spmv[0] = _mat_vec_p_l_csr_mkl;
        spmv[1] = _mat_vec_p_l_csr_mkl;
#else
        retcode = 2;
#endif
      }
      break;
    default:
      break;
    }

    break;

  case CS_MATRIX_CSR_SYM:

    switch(fill_type) {
    case CS_MATRIX_SCALAR_SYM:
      if (standard > 0) {
        spmv[0] = _mat_vec_p_l_csr_sym;
        spmv[1] = _mat_vec_p_l_csr_sym;
      }
      else if (!strcmp(func_name, "mkl")) {
#if defined(HAVE_MKL)
        spmv[0] = _mat_vec_p_l_csr_sym_mkl;
        spmv[1] = _mat_vec_p_l_csr_sym_mkl;
#else
        retcode = 2;
#endif
      }
      break;
    default:
      break;
    }
    break;

 case CS_MATRIX_MSR:

    if (standard > 0) {
      switch(fill_type) {
      case CS_MATRIX_SCALAR:
      case CS_MATRIX_SCALAR_SYM:
        spmv[0] = _mat_vec_p_l_msr;
        spmv[1] = _mat_vec_p_l_msr;
        break;
      case CS_MATRIX_33_BLOCK_D:
      case CS_MATRIX_33_BLOCK_D_SYM:
        spmv[0] = _b_mat_vec_p_l_msr;
        spmv[1] = _b_mat_vec_p_l_msr;
        break;
      default:
        break;
      }
    }

    else if (!strcmp(func_name, "prefetch")) {
      switch(fill_type) {
      case CS_MATRIX_SCALAR:
      case CS_MATRIX_SCALAR_SYM:
        l_length[0] = 508;
        l_length[1] = 508;
        spmv[0] = _mat_vec_p_l_msr_pf;
        spmv[1] = _mat_vec_p_l_msr_pf;
        break;
      default:
        break;
      }
    }

    else if (!strcmp(func_name, "mkl")) {
#if defined(HAVE_MKL)
      switch(fill_type) {
      case CS_MATRIX_SCALAR:
      case CS_MATRIX_SCALAR_SYM:
        spmv[0] = _mat_vec_p_l_msr_mkl;
        spmv[1] = _mat_vec_p_l_msr_mkl;
        break;
      default:
        break;
      }
#else
      retcode = 2;
#endif
    }

    break;

  default:
    break;
  }

  if (ed_flag != 1 && spmv[0] != NULL) {
    vector_multiply[fill_type][0] = spmv[0];
    loop_length[fill_type][0] = l_length[0];
    retcode = 0;
  }
  if (ed_flag != 0 && spmv[0] != NULL) {
    vector_multiply[fill_type][1] = spmv[1];
    loop_length[fill_type][1] = l_length[1];
    retcode = 0;
  }

  return retcode;
}

/*----------------------------------------------------------------------------
 * Create a matrix container using a given structure.
 *
 * Note that the matrix container maps to the assigned structure,
 * so it must be destroyed before that structure.
 *
 * parameters:
 *   ms <-- Associated matrix structure
 *
 * returns:
 *   pointer to created matrix structure;
 *----------------------------------------------------------------------------*/

cs_matrix_t *
cs_matrix_create(const cs_matrix_structure_t  *ms)
{
  int i;
  cs_matrix_fill_type_t mft;
  cs_matrix_t *m;

  BFT_MALLOC(m, 1, cs_matrix_t);

  m->type = ms->type;

  /* Map shared structure */

  m->n_cells = ms->n_cells;
  m->n_cells_ext = ms->n_cells_ext;
  m->n_faces = ms->n_faces;

  for (i = 0; i < 4; i++) {
    m->db_size[i] = 1;
    m->eb_size[i] = 1;
  }
  m->fill_type = CS_MATRIX_N_FILL_TYPES;

  m->structure = ms->structure;

  m->face_cell = ms->face_cell;
  m->cell_num = ms->cell_num;
  m->halo = ms->halo;
  m->numbering = ms->numbering;

  for (mft = 0; mft < CS_MATRIX_N_FILL_TYPES; mft++) {
    for (i = 0; i < 2; i++) {
      m->loop_length[mft][i] = 0;
      m->vector_multiply[mft][i] = NULL;
    }
  }

  /* Define coefficients */

  switch(m->type) {
  case CS_MATRIX_NATIVE:
    m->coeffs = _create_coeff_native();
    break;
  case CS_MATRIX_CSR:
    m->coeffs = _create_coeff_csr();
    break;
  case CS_MATRIX_CSR_SYM:
    m->coeffs = _create_coeff_csr_sym();
    break;
  case CS_MATRIX_MSR:
    m->coeffs = _create_coeff_msr();
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Handling of matrixes in %s format\n"
                "is not operational yet."),
              _(cs_matrix_type_name[m->type]));
    break;
  }

  /* Set function pointers here */

  m->set_coefficients = NULL;

  for (mft = 0; mft < CS_MATRIX_N_FILL_TYPES; mft++)
    _set_spmv_func(m->type,
                   m->numbering,
                   mft,
                   2,    /* ed_flag */
                   NULL, /* func_name */
                   m->vector_multiply,
                   m->loop_length);

  switch(m->type) {

  case CS_MATRIX_NATIVE:

    m->set_coefficients = _set_coeffs_native;
    m->release_coefficients = _release_coeffs_native;
    m->copy_diagonal = _copy_diagonal_separate;
    break;

  case CS_MATRIX_CSR:
    m->set_coefficients = _set_coeffs_csr;
    m->release_coefficients = _release_coeffs_csr;
    m->copy_diagonal = _copy_diagonal_csr;
    break;

  case CS_MATRIX_CSR_SYM:
    m->set_coefficients = _set_coeffs_csr_sym;
    m->release_coefficients = _release_coeffs_csr_sym;
    m->copy_diagonal = _copy_diagonal_csr_sym;
    m->vector_multiply[CS_MATRIX_SCALAR_SYM][0] = _mat_vec_p_l_csr_sym;
    break;

  case CS_MATRIX_MSR:
    m->set_coefficients = _set_coeffs_msr;
    m->release_coefficients = _release_coeffs_msr;
    m->copy_diagonal = _copy_diagonal_separate;
    break;

  default:
    assert(0);
    break;

  }

  for (i = 0; i < CS_MATRIX_N_FILL_TYPES; i++) {
    if (m->vector_multiply[i][1] == NULL)
      m->vector_multiply[i][1] = m->vector_multiply[i][0];
  }

  return m;
}

/*----------------------------------------------------------------------------
 * Create a matrix container using a given variant.
 *
 * If the matrix variant is incompatible with the structure, it is ignored,
 * and defaults for that structure are used instead.
 *
 * parameters:
 *   ms <-- Associated matrix structure
 *   mv <-- Associated matrix variant
 *
 * returns:
 *   pointer to created matrix structure;
 *----------------------------------------------------------------------------*/

cs_matrix_t *
cs_matrix_create_by_variant(const cs_matrix_structure_t  *ms,
                            const cs_matrix_variant_t    *mv)
{
  cs_matrix_t *m = cs_matrix_create(ms);

  if (mv != NULL) {
    if (mv->type == ms->type) {
      for (int i = 0; i < CS_MATRIX_N_FILL_TYPES; i++) {
        for (int j = 0; j < 2; j++) {
          if (mv->vector_multiply[i][j] != NULL) {
            m->loop_length[i][j] = mv->loop_length[i][j];
            m->vector_multiply[i][j] = mv->vector_multiply[i][j];
          }
        }
      }
    }
  }

  return m;
}

/*----------------------------------------------------------------------------
 * Destroy a matrix structure.
 *
 * parameters:
 *   matrix <-> Pointer to matrix structure pointer
 *----------------------------------------------------------------------------*/

void
cs_matrix_destroy(cs_matrix_t **matrix)
{
  if (matrix != NULL && *matrix != NULL) {

    cs_matrix_t *m = *matrix;

    switch(m->type) {
    case CS_MATRIX_NATIVE:
      {
        cs_matrix_coeff_native_t *coeffs = m->coeffs;
        _destroy_coeff_native(&coeffs);
      }
      break;
    case CS_MATRIX_CSR:
      {
        cs_matrix_coeff_csr_t *coeffs = m->coeffs;
        _destroy_coeff_csr(&coeffs);
        m->coeffs = NULL;
      }
      break;
    case CS_MATRIX_CSR_SYM:
      {
        cs_matrix_coeff_csr_sym_t *coeffs = m->coeffs;
        _destroy_coeff_csr_sym(&coeffs);
        m->coeffs = NULL;
      }
      break;
    case CS_MATRIX_MSR:
      {
        cs_matrix_coeff_msr_t *coeffs = m->coeffs;
        _destroy_coeff_msr(&coeffs);
        m->coeffs = NULL;
      }
      break;
    default:
      assert(0);
      break;
    }

    m->coeffs = NULL;

    /* Now free main structure */

    BFT_FREE(*matrix);
  }
}

/*----------------------------------------------------------------------------
 * Return number of columns in matrix.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_matrix_get_n_columns(const cs_matrix_t  *matrix)
{
  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("The matrix is not defined."));
  return matrix->n_cells_ext;
}

/*----------------------------------------------------------------------------
 * Return number of rows in matrix.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_matrix_get_n_rows(const cs_matrix_t  *matrix)
{
  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("The matrix is not defined."));
  return matrix->n_cells;
}

/*----------------------------------------------------------------------------
 * Return matrix diagonal block sizes.
 *
 * Block sizes are defined by a array of 4 values:
 *   0: useful block size, 1: vector block extents,
 *   2: matrix line extents,  3: matrix line*column extents
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *
 * returns:
 *   pointer to block sizes
 *----------------------------------------------------------------------------*/

const int *
cs_matrix_get_diag_block_size(const cs_matrix_t  *matrix)
{
  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("The matrix is not defined."));

  return matrix->db_size;
}

/*----------------------------------------------------------------------------
 * Get matrix fill type, depending on block sizes.
 *
 * Block sizes are defined by an optional array of 4 values:
 *   0: useful block size, 1: vector block extents,
 *   2: matrix line extents,  3: matrix line*column extents
 *
 * parameters:
 *   symmetric              <-- Indicates if matrix coefficients are symmetric
 *   diag_block_size        <-- Block sizes for diagonal, or NULL
 *   extra_diag_block_size  <-- Block sizes for extra diagonal, or NULL
 *
 * returns:
 *   matrix fill type
 *----------------------------------------------------------------------------*/

cs_matrix_fill_type_t
cs_matrix_get_fill_type(bool        symmetric,
                        const int  *diag_block_size,
                        const int  *extra_diag_block_size)
{
  cs_matrix_fill_type_t fill_type = CS_MATRIX_N_FILL_TYPES;

  int _db_size = 1, _eb_size = 1;
  if (diag_block_size != NULL)
    _db_size = diag_block_size[0];

  if (extra_diag_block_size != NULL)
    _eb_size = extra_diag_block_size[0];

  /* Set fill type */

  if (_eb_size == 3)
    fill_type = CS_MATRIX_33_BLOCK;
  else if (_db_size == 3) {
    if (symmetric)
      fill_type = CS_MATRIX_33_BLOCK_D_SYM;
    else
      fill_type = CS_MATRIX_33_BLOCK_D;
  }
  else if (_db_size == 1) {
    if (symmetric)
      fill_type = CS_MATRIX_SCALAR_SYM;
    else
      fill_type = CS_MATRIX_SCALAR;
  }
  return fill_type;
}

/*----------------------------------------------------------------------------
 * Set matrix coefficients, sharing arrays with the caller when possible.
 *
 * With shared arrays, the matrix becomes unusable if the arrays passed as
 * arguments are not be modified (its coefficients should be unset first
 * to mark this).
 *
 * Depending on current options and initialization, values will be copied
 * or simply mapped.
 *
 * Block sizes are defined by an optional array of 4 values:
 *   0: useful block size, 1: vector block extents,
 *   2: matrix line extents,  3: matrix line*column extents
 *
 * parameters:
 *   matrix                 <-> Pointer to matrix structure
 *   symmetric              <-- Indicates if matrix coefficients are symmetric
 *   diag_block_size        <-- Block sizes for diagonal, or NULL
 *   extra_diag_block_size  <-- Block sizes for extra diagonal, or NULL
 *   da                     <-- Diagonal values (NULL if zero)
 *   xa                     <-- Extradiagonal values (NULL if zero)
 *----------------------------------------------------------------------------*/

void
cs_matrix_set_coefficients(cs_matrix_t      *matrix,
                           bool              symmetric,
                           const int        *diag_block_size,
                           const int        *extra_diag_block_size,
                           const cs_real_t  *da,
                           const cs_real_t  *xa)
{
  int i;

  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("The matrix is not defined."));

  if (diag_block_size == NULL) {
    for (i = 0; i < 4; i++)
      matrix->db_size[i] = 1;
  }
  else {
    for (i = 0; i < 4; i++)
      matrix->db_size[i] = diag_block_size[i];
  }

  if (extra_diag_block_size == NULL) {
    for (i = 0; i < 4; i++)
      matrix->eb_size[i] = 1;
  }
  else {
    for (i = 0; i < 4; i++)
      matrix->eb_size[i] = extra_diag_block_size[i];
  }

  /* Set fill type */

  if (matrix->eb_size[0] == 3)
    matrix->fill_type = CS_MATRIX_33_BLOCK;
  else if (matrix->db_size[0] == 3) {
    if (symmetric)
      matrix->fill_type = CS_MATRIX_33_BLOCK_D_SYM;
    else
      matrix->fill_type = CS_MATRIX_33_BLOCK_D;
  }
  else if (matrix->db_size[0] == 1) {
    if (symmetric)
      matrix->fill_type = CS_MATRIX_SCALAR_SYM;
    else
      matrix->fill_type = CS_MATRIX_SCALAR;
  }

  /* Set coefficients */

  if (matrix->set_coefficients != NULL)
    matrix->set_coefficients(matrix, symmetric, true, false, da, xa);
}

/*----------------------------------------------------------------------------
 * Set matrix coefficients in the non-interleaved case.
 *
 * In the symmetric case, there is no difference with the interleaved case.
 *
 * Depending on current options and initialization, values will be copied
 * or simply mapped (non-symmetric values will be copied).
 *
 * parameters:
 *   matrix    <-> Pointer to matrix structure
 *   symmetric <-- Indicates if matrix coefficients are symmetric
 *   da        <-- Diagonal values (NULL if zero)
 *   xa        <-- Extradiagonal values (NULL if zero)
 *----------------------------------------------------------------------------*/

void
cs_matrix_set_coefficients_ni(cs_matrix_t      *matrix,
                              bool              symmetric,
                              const cs_real_t  *da,
                              const cs_real_t  *xa)
{
  int i;

  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("The matrix is not defined."));

  for (i = 0; i < 4; i++)
    matrix->db_size[i] = 1;

  for (i = 0; i < 4; i++)
    matrix->eb_size[i] = 1;

  /* Set fill type */

  if (symmetric)
    matrix->fill_type = CS_MATRIX_SCALAR_SYM;
  else
    matrix->fill_type = CS_MATRIX_SCALAR;
  if (matrix->set_coefficients != NULL)
    matrix->set_coefficients(matrix, symmetric, false, false, da, xa);

  /* Set coefficients */
}

/*----------------------------------------------------------------------------
 * Set matrix coefficients, copying values to private arrays.
 *
 * With private arrays, the matrix becomes independant from the
 * arrays passed as arguments.
 *
 * Block sizes are defined by an optional array of 4 values:
 *   0: useful block size, 1: vector block extents,
 *   2: matrix line extents,  3: matrix line*column extents
 *
 * parameters:
 *   matrix                 <-> Pointer to matrix structure
 *   symmetric              <-- Indicates if matrix coefficients are symmetric
 *   diag_block_size        <-- Block sizes for diagonal, or NULL
 *   extra_diag_block_size  <-- Block sizes for extra diagonal, or NULL
 *   da                     <-- Diagonal values (NULL if zero)
 *   xa                     <-- Extradiagonal values (NULL if zero)
 *----------------------------------------------------------------------------*/

void
cs_matrix_copy_coefficients(cs_matrix_t      *matrix,
                            bool              symmetric,
                            const int        *diag_block_size,
                            const int        *extra_diag_block_size,
                            const cs_real_t  *da,
                            const cs_real_t  *xa)
{
  int i;

  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("The matrix is not defined."));

  if (diag_block_size == NULL) {
    for (i = 0; i < 4; i++)
      matrix->db_size[i] = 1;
  }
  else {
    for (i = 0; i < 4; i++)
      matrix->db_size[i] = diag_block_size[i];
  }

  if (extra_diag_block_size == NULL) {
    for (i = 0; i < 4; i++)
      matrix->eb_size[i] = 1;
  }
  else {
    for (i = 0; i < 4; i++)
      matrix->eb_size[i] = extra_diag_block_size[i];
  }

  if (matrix->set_coefficients != NULL)
    matrix->set_coefficients(matrix, symmetric, true, true, da, xa);

  /* Set fill type */

  if (matrix->eb_size[1] == 3)
    matrix->fill_type = CS_MATRIX_33_BLOCK;
  else if (matrix->db_size[1] == 3) {
    if (symmetric)
      matrix->fill_type = CS_MATRIX_33_BLOCK_D_SYM;
    else
      matrix->fill_type = CS_MATRIX_33_BLOCK_D;
  }
  else if (matrix->db_size[1] == 1) {
    if (symmetric)
      matrix->fill_type = CS_MATRIX_SCALAR_SYM;
    else
      matrix->fill_type = CS_MATRIX_SCALAR;
  }
}

/*----------------------------------------------------------------------------
 * Release shared matrix coefficients.
 *
 * Pointers to mapped coefficients are set to NULL, while
 * coefficient copies owned by the matrix are not modified.
 *
 * This simply ensures the matrix does not maintain pointers
 * to nonexistant data.
 *
 * parameters:
 *   matrix <-> Pointer to matrix structure
 *----------------------------------------------------------------------------*/

void
cs_matrix_release_coefficients(cs_matrix_t  *matrix)
{
  /* Check API state */

  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("The matrix is not defined."));

  if (matrix->release_coefficients != NULL)
    matrix->release_coefficients(matrix);

  /* Set fill type to impossible value */

  matrix->fill_type = CS_MATRIX_N_FILL_TYPES;
}

/*----------------------------------------------------------------------------
 * Copy matrix diagonal values.
 *
 * In case of matrixes with block diagonal coefficients, only the true
 * diagonal values are copied.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *   da     --> Diagonal (pre-allocated, size: n_cells)
 *----------------------------------------------------------------------------*/

void
cs_matrix_copy_diagonal(const cs_matrix_t  *matrix,
                        cs_real_t          *restrict da)
{
  /* Check API state */

  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("The matrix is not defined."));

  if (matrix->copy_diagonal != NULL)
    matrix->copy_diagonal(matrix, da);
}

/*----------------------------------------------------------------------------
 * Get matrix diagonal values.
 *
 * In case of matrixes with block diagonal coefficients, a pointer to
 * the complete block diagonal is returned.
 *
 * parameters:
 *   matrix --> Pointer to matrix structure
 *
 * returns:
 *   pointer to matrix diagonal array
 *----------------------------------------------------------------------------*/

const cs_real_t *
cs_matrix_get_diagonal(const cs_matrix_t  *matrix)
{
  cs_lnum_t ii;

  const cs_real_t  *diag = NULL;

  switch(matrix->type) {

  case CS_MATRIX_NATIVE:
    {
      cs_matrix_coeff_native_t *mc = matrix->coeffs;
      if (mc->da == NULL) {
        cs_lnum_t n_rows = matrix->n_cells * matrix->db_size[3];
        if (mc->_da == NULL || mc->max_db_size < matrix->db_size[3]) {
          BFT_REALLOC(mc->_da, matrix->db_size[3]*matrix->n_cells, cs_real_t);
          mc->max_db_size = matrix->db_size[3];
        }
#       pragma omp parallel for  if(n_rows > CS_THR_MIN)
        for (ii = 0; ii < n_rows; ii++)
          mc->_da[ii] = 0.0;
        mc->da = mc->_da;
      }
      diag = mc->da;
    }
    break;

  case CS_MATRIX_CSR:
    {
      cs_matrix_coeff_csr_t *mc = matrix->coeffs;
      assert(matrix->db_size[3] == 1);
      if (mc->_d_val == NULL)
        BFT_MALLOC(mc->_d_val, matrix->n_cells, cs_real_t);
      if (mc->d_val == NULL) {
        cs_matrix_copy_diagonal(matrix, mc->_d_val);
        mc->d_val = mc->_d_val;
      }
      diag = mc->d_val;
    }
    break;

  case CS_MATRIX_CSR_SYM:
    {
      cs_matrix_coeff_csr_sym_t *mc = matrix->coeffs;
      assert(matrix->db_size[3] == 1);
      if (mc->_d_val == NULL)
        BFT_MALLOC(mc->_d_val, matrix->n_cells, cs_real_t);
      if (mc->d_val == NULL) {
        cs_matrix_copy_diagonal(matrix, mc->_d_val);
        mc->d_val = mc->_d_val;
      }
      diag = mc->d_val;
    }
    break;

  case CS_MATRIX_MSR:
    {
      cs_matrix_coeff_msr_t *mc = matrix->coeffs;
      if (mc->d_val == NULL) {
        cs_lnum_t n_rows = matrix->n_cells * matrix->db_size[3];
        if (mc->_d_val == NULL || mc->max_db_size < matrix->db_size[3]) {
          BFT_REALLOC(mc->_d_val, matrix->db_size[3]*matrix->n_cells, cs_real_t);
          mc->max_db_size = matrix->db_size[3];
        }
#       pragma omp parallel for  if(n_rows > CS_THR_MIN)
        for (ii = 0; ii < n_rows; ii++)
          mc->_d_val[ii] = 0.0;
        mc->d_val = mc->_d_val;
      }
      diag = mc->d_val;
    }
    break;

  default:
    assert(0);
    break;
  }

  return diag;
}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x
 *
 * This function includes a halo update of x prior to multiplication by A.
 *
 * parameters:
 *   rotation_mode <-- Halo update option for rotational periodicity
 *   matrix        <-- Pointer to matrix structure
 *   x             <-> Multipliying vector values (ghost values updated)
 *   y             --> Resulting vector
 *----------------------------------------------------------------------------*/

void
cs_matrix_vector_multiply(cs_halo_rotation_t   rotation_mode,
                          const cs_matrix_t   *matrix,
                          cs_real_t           *restrict x,
                          cs_real_t           *restrict y)
{
  assert(matrix != NULL);

  if (matrix->halo != NULL)
    _pre_vector_multiply_sync(rotation_mode,
                              matrix,
                              x,
                              y);

  if (matrix->vector_multiply[matrix->fill_type][0] != NULL)
    matrix->vector_multiply[matrix->fill_type][0](false, matrix, x, y);
  else
    bft_error
      (__FILE__, __LINE__, 0,
       _("Matrix is missing a vector multiply function for fill type %s."),
       cs_matrix_fill_type_name[matrix->fill_type]);
}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with no prior halo update of x.
 *
 * This function does not include a halo update of x prior to multiplication
 * by A, so it should be called only when the halo of x is known to already
 * be up to date (in which case we avoid the performance penalty of a
 * redundant update by using this variant of the matrix.vector product).
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *   x      <-- Multipliying vector values
 *   y      --> Resulting vector
 *----------------------------------------------------------------------------*/

void
cs_matrix_vector_multiply_nosync(const cs_matrix_t  *matrix,
                                 const cs_real_t    *x,
                                 cs_real_t          *restrict y)
{
  assert(matrix != NULL);

  if (matrix->vector_multiply[matrix->fill_type][0] != NULL)
    matrix->vector_multiply[matrix->fill_type][0](false, matrix, x, y);
  else
    bft_error
      (__FILE__, __LINE__, 0,
       _("Matrix is missing a vector multiply function for fill type %s."),
       cs_matrix_fill_type_name[matrix->fill_type]);
}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = (A-D).x
 *
 * This function includes a halo update of x prior to multiplication by A.
 *
 * parameters:
 *   rotation_mode <-- Halo update option for rotational periodicity
 *   matrix        <-- Pointer to matrix structure
 *   x             <-> Multipliying vector values (ghost values updated)
 *   y             --> Resulting vector
 *----------------------------------------------------------------------------*/

void
cs_matrix_exdiag_vector_multiply(cs_halo_rotation_t   rotation_mode,
                                 const cs_matrix_t   *matrix,
                                 cs_real_t           *restrict x,
                                 cs_real_t           *restrict y)
{
  assert(matrix != NULL);

  if (matrix->halo != NULL)
    _pre_vector_multiply_sync(rotation_mode,
                              matrix,
                              x,
                              y);

  if (matrix->vector_multiply[matrix->fill_type][1] != NULL)
    matrix->vector_multiply[matrix->fill_type][1](true, matrix, x, y);
  else
    bft_error
      (__FILE__, __LINE__, 0,
       _("Matrix is missing a vector multiply function for fill type %s."),
       cs_matrix_fill_type_name[matrix->fill_type]);
}

/*----------------------------------------------------------------------------
 * Build matrix variant
 *
 * The variant will initially use default matrix-vector functions,
 * which can be later modified using cs_matrix_variant_set_func().
 *
 * parameters:
 *   type         <-- Type of matrix considered
 *   numbering    <-- vectorization or thread-related numbering info,
 *                    or NULL
 *----------------------------------------------------------------------------*/

cs_matrix_variant_t *
cs_matrix_variant_create(cs_matrix_type_t         type,
                         const cs_numbering_t    *numbering)

{
  cs_matrix_fill_type_t mft;
  cs_matrix_variant_t  *mv;

  BFT_MALLOC(mv, 1, cs_matrix_variant_t);

  _variant_init(mv);

  mv->type = type;

  strncpy(mv->name, cs_matrix_type_name[type], 31);
  mv->name[31] = '\0';

  for (mft = 0; mft < CS_MATRIX_N_FILL_TYPES; mft++) {
    (void) _set_spmv_func(type,
                          numbering,
                          mft,
                          2,
                          NULL, /* func_name */
                          mv->vector_multiply,
                          mv->loop_length);
  }

  return mv;
}

/*----------------------------------------------------------------------------
 * Build list of variants for tuning or testing.
 *
 * parameters:
 *   n_fill_types <-- number of fill types tuned for
 *   fill_types   <-- array of fill types tuned for
 *   type_filter  <-- true for matrix types tuned for, false for others
 *   numbering    <-- vectorization or thread-related numbering info,
 *                    or NULL
 *   n_variants   --> number of variants
 *   m_variant    --> array of matrix variants
 *----------------------------------------------------------------------------*/

void
cs_matrix_variant_build_list(int                      n_fill_types,
                             cs_matrix_fill_type_t    fill_types[],
                             bool                     type_filter[],
                             const cs_numbering_t    *numbering,
                             int                     *n_variants,
                             cs_matrix_variant_t    **m_variant)
{
  int  i;
  int  _n_fill_types;
  cs_matrix_fill_type_t  _fill_types[CS_MATRIX_N_FILL_TYPES];

  int  n_variants_max = 0;

  *n_variants = 0;
  *m_variant = NULL;

  if (type_filter[CS_MATRIX_NATIVE]) {

    _variant_add(_("Native, baseline"),
                 CS_MATRIX_NATIVE,
                 n_fill_types,
                 fill_types,
                 2, /* ed_flag */
                 0, /* loop_length */
                 _mat_vec_p_l_native,
                 _b_mat_vec_p_l_native,
                 _bb_mat_vec_p_l_native,
                 n_variants,
                 &n_variants_max,
                 m_variant);

    _variant_add(_("Native, 3x3 blocks"),
                 CS_MATRIX_NATIVE,
                 n_fill_types,
                 fill_types,
                 2, /* ed_flag */
                 0, /* loop_length */
                 NULL,
                 _3_3_mat_vec_p_l_native,
                 NULL,
                 n_variants,
                 &n_variants_max,
                 m_variant);

    _variant_add(_("Native, Bull algorithm"),
                 CS_MATRIX_NATIVE,
                 n_fill_types,
                 fill_types,
                 2, /* ed_flag */
                 508, /* loop_length */
                 _mat_vec_p_l_native_bull,
                 NULL,
                 NULL,
                 n_variants,
                 &n_variants_max,
                 m_variant);

    if (numbering != NULL) {

#if defined(HAVE_OPENMP)

      if (numbering->type == CS_NUMBERING_THREADS)
        _variant_add(_("Native, OpenMP"),
                     CS_MATRIX_NATIVE,
                     n_fill_types,
                     fill_types,
                     2, /* ed_flag */
                     0, /* loop_length */
                     _mat_vec_p_l_native_omp,
                     _b_mat_vec_p_l_native_omp,
                     NULL,
                     n_variants,
                     &n_variants_max,
                     m_variant);

      _variant_add(_("Native, OpenMP atomic"),
                   CS_MATRIX_NATIVE,
                   n_fill_types,
                   fill_types,
                   2, /* ed_flag */
                   0, /* loop_length */
                   _mat_vec_p_l_native_omp_atomic,
                   _b_mat_vec_p_l_native_omp_atomic,
                   NULL,
                   n_variants,
                   &n_variants_max,
                   m_variant);

#endif

#if defined(SX) && defined(_SX) /* For vector machines */
      if (numbering->type == CS_NUMBERING_VECTORIZE)
        _variant_add(_("Native, vectorized"),
                     CS_MATRIX_NATIVE,
                     n_fill_types,
                     fill_types,
                     2, /* ed_flag */
                     0, /* loop_length */
                     _mat_vec_p_l_native_vector,
                     NULL,
                     NULL,
                     n_variants,
                     &n_variants_max,
                     m_variant);
#endif

    }

  }

  if (type_filter[CS_MATRIX_CSR]) {

    _variant_add(_("CSR"),
                 CS_MATRIX_CSR,
                 n_fill_types,
                 fill_types,
                 2, /* ed_flag */
                 0, /* loop_length */
                 _mat_vec_p_l_csr,
                 NULL,
                 NULL,
                 n_variants,
                 &n_variants_max,
                 m_variant);

    _variant_add(_("CSR, with prefetch"),
                 CS_MATRIX_CSR,
                 n_fill_types,
                 fill_types,
                 0, /* ed_flag */
                 508, /* loop_length */
                 _mat_vec_p_l_csr_pf,
                 NULL,
                 NULL,
                 n_variants,
                 &n_variants_max,
                 m_variant);

#if defined(HAVE_MKL)

    _variant_add(_("CSR, with MKL"),
                 CS_MATRIX_CSR,
                 n_fill_types,
                 fill_types,
                 0, /* ed_flag */
                 0, /* loop_length */
                 _mat_vec_p_l_csr_mkl,
                 NULL,
                 NULL,
                 n_variants,
                 &n_variants_max,
                 m_variant);

#endif /* defined(HAVE_MKL) */

  }

  if (type_filter[CS_MATRIX_CSR_SYM]) {

    for (i = 0, _n_fill_types = 0; i < n_fill_types; i++) {
      if (fill_types[i] == CS_MATRIX_SCALAR_SYM) {
        _fill_types[_n_fill_types++] = fill_types[i];
      }
    }

    if (_n_fill_types > 0) {

      _variant_add(_("CSR_SYM"),
                   CS_MATRIX_CSR_SYM,
                   _n_fill_types,
                   _fill_types,
                   2, /* ed_flag */
                   0, /* loop_length */
                   _mat_vec_p_l_csr_sym,
                   NULL,
                   NULL,
                   n_variants,
                   &n_variants_max,
                   m_variant);

#if defined(HAVE_MKL)

      _variant_add(_("CSR_SYM, with MKL"),
                   CS_MATRIX_CSR_SYM,
                   _n_fill_types,
                   _fill_types,
                   0, /* ed_flag */
                   0, /* loop_length */
                   _mat_vec_p_l_csr_sym_mkl,
                   NULL,
                   NULL,
                   n_variants,
                   &n_variants_max,
                   m_variant);

#endif /* defined(HAVE_MKL) */

    }

  }

  if (type_filter[CS_MATRIX_MSR]) {

    _variant_add(_("MSR"),
                 CS_MATRIX_MSR,
                 n_fill_types,
                 fill_types,
                 2, /* ed_flag */
                 0, /* loop_length */
                 _mat_vec_p_l_msr,
                 _b_mat_vec_p_l_msr,
                 NULL,
                 n_variants,
                 &n_variants_max,
                 m_variant);

    _variant_add(_("MSR, with prefetch"),
                 CS_MATRIX_MSR,
                 n_fill_types,
                 fill_types,
                 2, /* ed_flag */
                 508, /* loop_length */
                 _mat_vec_p_l_msr_pf,
                 NULL,
                 NULL,
                 n_variants,
                 &n_variants_max,
                 m_variant);

#if defined(HAVE_MKL)

    _variant_add(_("MSR, with MKL"),
                 CS_MATRIX_MSR,
                 n_fill_types,
                 fill_types,
                 2, /* ed_flag */
                 0, /* loop_length */
                 _mat_vec_p_l_msr_mkl,
                 NULL,
                 NULL,
                 n_variants,
                 &n_variants_max,
                 m_variant);

#endif /* defined(HAVE_MKL) */

  }

  n_variants_max = *n_variants;
  BFT_REALLOC(*m_variant, *n_variants, cs_matrix_variant_t);
}

/*----------------------------------------------------------------------------
 * Destroy a matrix variant structure.
 *
 * parameters:
 *   mv <-> Pointer to matrix variant pointer
 *----------------------------------------------------------------------------*/

void
cs_matrix_variant_destroy(cs_matrix_variant_t  **mv)
{
  if (mv != NULL)
    BFT_FREE(*mv);
}

/*----------------------------------------------------------------------------
 * Select the sparse matrix-vector product function to be used by a
 * matrix variant for a given fill type.
 *
 * Currently, possible variant functions are:
 *
 *   CS_MATRIX_NATIVE  (all fill types)
 *     standard
 *     3_3_diag        (for CS_MATRIX_33_BLOCK_D or CS_MATRIX_33_BLOCK_D_SYM)
 *     bull            (for CS_MATRIX_SCALAR or CS_MATRIX_SCALAR_SYM)
 *     omp             (for OpenMP with compatible numbering)
 *     vector          (For vector machine with compatible numbering)
 *
 *   CS_MATRIX_CSR     (for CS_MATRIX_SCALAR or CS_MATRIX_SCALAR_SYM)
 *     standard
 *     prefetch
 *     mkl             (with MKL)
 *
 *   CS_MATRIX_CSR_SYM (for CS_MATRIX_SCALAR_SYM)
 *     standard
 *     mkl             (with MKL)
 *
 *   CS_MATRIX_MSR     (all fill types except CS_MATRIX_33_BLOCK)
 *     standard
 *     prefetch
 *     mkl             (with MKL, for CS_MATRIX_SCALAR or CS_MATRIX_SCALAR_SYM)
 *
 * parameters:
 *   mv        <-> Pointer to matrix variant
 *   numbering <-- mesh numbering info, or NULL
 *   fill type <-- matrix fill type to merge from
 *   ed_flag   <-- 0: with diagonal only, 1 exclude only; 2; both
 *   func_name <-- function type name
 *----------------------------------------------------------------------------*/

void
cs_matrix_variant_set_func(cs_matrix_variant_t     *mv,
                           const cs_numbering_t    *numbering,
                           cs_matrix_fill_type_t    fill_type,
                           int                      ed_flag,
                           const char              *func_name)
{
  int retcode = _set_spmv_func(mv->type,
                               numbering,
                               fill_type,
                               ed_flag,
                               func_name,
                               mv->vector_multiply,
                               mv->loop_length);

  if (retcode == 1)
    bft_error
      (__FILE__, __LINE__, 0,
       _("Assignment of matrix.vector product \"%s\" to matrix variant \"%s\"\n"
         "of type \"%s\" for fill \"%s\" not allowed."),
       func_name, mv->name, cs_matrix_type_name[mv->type],
       cs_matrix_fill_type_name[fill_type]);
  else if (retcode == 2)
    bft_error
      (__FILE__, __LINE__, 0,
       _("Matrix.vector product function type \"%s\"\n"
         "is not available in this build."),
       func_name);
}

/*----------------------------------------------------------------------------
 * Merge a functions to a matrix variant from another variant sharing
 * the same structure.
 *
 * Functions from the structure to merge for the selected fill type are
 * assigned to the main variant.
 *
 * This can be useful when tuning has been done separately for different fill
 * types, and the resulting selected structure is identical.
 *
 * parameters:
 *   mv        <-> Pointer to matrix variant
 *   mv_merge  <-- Pointer to matrix variant to merge
 *   fill type <-- matrix fill type to merge from
 *----------------------------------------------------------------------------*/

void
cs_matrix_variant_merge(cs_matrix_variant_t        *mv,
                        const cs_matrix_variant_t  *mv_merge,
                        cs_matrix_fill_type_t       fill_type)
{
  if (mv->type == mv_merge->type) {
    for (int i = 0; i < 2; i++) {
      mv->vector_multiply[fill_type][i]
        = mv_merge->vector_multiply[fill_type][i];
      mv->matrix_vector_cost[fill_type][i]
        = mv_merge->matrix_vector_cost[fill_type][i];
    }
    mv->matrix_assign_cost[fill_type] = mv_merge->matrix_assign_cost[fill_type];
  }
}

/*----------------------------------------------------------------------------
 * Get the type associated with a matrix variant.
 *
 * parameters:
 *   mv <-- Pointer to matrix variant structure
 *----------------------------------------------------------------------------*/

cs_matrix_type_t
cs_matrix_variant_type(const cs_matrix_variant_t  *mv)
{
  return mv->type;
}

/*----------------------------------------------------------------------------
 * Test local matrix.vector product operations.
 *
 * parameters:
 *   n_cells        <-- number of local cells
 *   n_cells_ext    <-- number of cells including ghost cells (array size)
 *   n_faces        <-- local number of internal faces
 *   cell_num       <-- Optional global cell numbers (1 to n), or NULL
 *   face_cell      <-- face -> cells connectivity (1 to n)
 *   halo           <-- cell halo structure
 *   numbering      <-- vectorization or thread-related numbering info, or NULL
 *
 * returns:
 *   pointer to tuning results structure
 *----------------------------------------------------------------------------*/

void
cs_matrix_variant_test(cs_lnum_t              n_cells,
                       cs_lnum_t              n_cells_ext,
                       cs_lnum_t              n_faces,
                       const cs_gnum_t       *cell_num,
                       const cs_lnum_2_t     *face_cell,
                       const cs_halo_t       *halo,
                       const cs_numbering_t  *numbering)
{
  int  n_variants = 0;
  bool type_filter[CS_MATRIX_N_TYPES] = {true, true, true, true};
  cs_matrix_fill_type_t  fill_types[] = {CS_MATRIX_SCALAR,
                                         CS_MATRIX_SCALAR_SYM,
                                         CS_MATRIX_33_BLOCK_D,
                                         CS_MATRIX_33_BLOCK_D_SYM,
                                         CS_MATRIX_33_BLOCK};
  cs_matrix_variant_t  *m_variant = NULL;

  /* Test basic flag combinations */

  bft_printf
    (_("\n"
       "Checking matrix structure and operation variants (diff/reference):\n"
       "------------------------------------------------\n"));

  /* Build variants array */

  cs_matrix_variant_build_list(CS_MATRIX_N_FILL_TYPES,
                               fill_types,
                               type_filter,
                               numbering,
                               &n_variants,
                               &m_variant);

  /* Run tests on variants */

  _matrix_check(n_variants,
                n_cells,
                n_cells_ext,
                n_faces,
                cell_num,
                face_cell,
                halo,
                numbering,
                m_variant);

  n_variants = 0;
  BFT_FREE(m_variant);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
