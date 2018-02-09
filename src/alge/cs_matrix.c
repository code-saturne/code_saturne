/*============================================================================
 * Sparse Matrix Representation and Operations.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
#include "cs_sort.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_matrix.h"
#include "cs_matrix_assembler.h"
#include "cs_matrix_priv.h"

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

/* Cache line multiple, in cs_real_t units */

#define CS_CL  (CS_CL_SIZE/8)

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
                                           "CS_MATRIX_BLOCK_D",
                                           "CS_MATRIX_BLOCK_D_66",
                                           "CS_MATRIX_BLOCK_D_SYM",
                                           "CS_MATRIX_BLOCK"};

#if defined (HAVE_MKL)

static char _no_exclude_diag_error_str[]
  = N_("Matrix product variant using function %s\n"
       "does not handle case with excluded diagonal.");

#endif

static const char *_matrix_operation_name[CS_MATRIX_N_FILL_TYPES][2]
  = {{N_("y <- A.x"),
      N_("y <- (A-D).x")},
     {N_("Symmetric y <- A.x"),
      N_("Symmetric y <- (A-D).x")},
     {N_("Block diagonal y <- A.x"),
      N_("Block diagonal y <- (A-D).x")},
     {N_("Block 6 diagonal y <- A.x"),
      N_("Block 6 diagonal y <- (A-D).x")},
     {N_("Block diagonal symmetric y <- A.x"),
      N_("Block diagonal symmetric y <- (A-D).x")},
     {N_("Block y <- A.x"),
      N_("Block y <- (A-D).x")}};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Set matrix fill metadata.
 *
 * Block sizes are defined by an optional array of 4 values:
 *   0: useful block size, 1: vector block extents,
 *   2: matrix line extents,  3: matrix line*column extents
 *
 * parameters:
 *   matrix               <-> pointer to matrix structure
 *   symmetric            <-- indicates if matrix coefficients are symmetric
 *   diag_block_size       <-- block sizes for diagonal, or NULL
 *   extra_diag_block_size <-- block sizes for extra diagonal, or NULL
 *----------------------------------------------------------------------------*/

static void
_set_fill_info(cs_matrix_t      *matrix,
               bool              symmetric,
               const cs_lnum_t   diag_block_size[4],
               const cs_lnum_t   extra_diag_block_size[4])
{
  matrix->symmetric = symmetric;

  if (diag_block_size == NULL) {
    for (int i = 0; i < 4; i++)
      matrix->db_size[i] = 1;
  }
  else {
    for (int i = 0; i < 4; i++)
      matrix->db_size[i] = diag_block_size[i];
  }

  if (extra_diag_block_size == NULL) {
    for (int i = 0; i < 4; i++)
      matrix->eb_size[i] = 1;
  }
  else {
    for (int i = 0; i < 4; i++)
      matrix->eb_size[i] = extra_diag_block_size[i];
  }

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
  if (matrix->type != CS_MATRIX_CSR_SYM)
    matrix->symmetric = false;

  for (int i = 0; i < 4; i++) {
    matrix->db_size[i] = 0;
    matrix->eb_size[i] = 0;
  }

  matrix->fill_type = CS_MATRIX_N_FILL_TYPES;
}

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
 *   a      <-- pointer to block matrixes array (usually matrix diagonal)
 *   x      <-- multipliying vector values
 *   y      --> resulting vector
 *----------------------------------------------------------------------------*/

static inline void
_dense_b_ax(cs_lnum_t         b_id,
            const cs_lnum_t   b_size[4],
            const cs_real_t   a[restrict],
            const cs_real_t   x[restrict],
            cs_real_t         y[restrict])
{
  cs_lnum_t   ii, jj;

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
 *   a      <-- pointer to block matrixes array (usually matrix diagonal)
 *   x      <-- multipliying vector values
 *   y      --> resulting vector
 *----------------------------------------------------------------------------*/

static inline void
_dense_3_3_ax(cs_lnum_t         b_id,
              const cs_real_t   a[restrict],
              const cs_real_t   x[restrict],
              cs_real_t         y[restrict])
{
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
 * Compute matrix-vector product for one dense block: y[i] = a[i].x[i]
 *
 * This variant uses a fixed 6x6 block, for better compiler optimization.
 *
 * parameters:
 *   b_id   <-- block id
 *   a      <-- pointer to block matrixes array (usually matrix diagonal)
 *   x      <-- multipliying vector values
 *   y      --> resulting vector
 *----------------------------------------------------------------------------*/

static inline void
_dense_6_6_ax(cs_lnum_t        b_id,
              const cs_real_t  a[restrict],
              const cs_real_t  x[restrict],
              cs_real_t        y[restrict])
{
  const cs_lnum_t b_id_6 = b_id*6, b_id_36 = b_id*36;

  y[b_id_6]     =   a[b_id_36]         * x[b_id_6]
                  + a[b_id_36 + 1]     * x[b_id_6 + 1]
                  + a[b_id_36 + 2]     * x[b_id_6 + 2]
                  + a[b_id_36 + 3]     * x[b_id_6 + 3]
                  + a[b_id_36 + 4]     * x[b_id_6 + 4]
                  + a[b_id_36 + 5]     * x[b_id_6 + 5];

  y[b_id_6 + 1] =   a[b_id_36 + 6]     * x[b_id_6]
                  + a[b_id_36 + 6 + 1] * x[b_id_6 + 1]
                  + a[b_id_36 + 6 + 2] * x[b_id_6 + 2]
                  + a[b_id_36 + 6 + 3] * x[b_id_6 + 3]
                  + a[b_id_36 + 6 + 4] * x[b_id_6 + 4]
                  + a[b_id_36 + 6 + 5] * x[b_id_6 + 5];

  y[b_id_6 + 2] =   a[b_id_36 + 12]     * x[b_id_6]
                  + a[b_id_36 + 12 + 1] * x[b_id_6 + 1]
                  + a[b_id_36 + 12 + 2] * x[b_id_6 + 2]
                  + a[b_id_36 + 12 + 3] * x[b_id_6 + 3]
                  + a[b_id_36 + 12 + 4] * x[b_id_6 + 4]
                  + a[b_id_36 + 12 + 5] * x[b_id_6 + 5];

  y[b_id_6 + 3] =   a[b_id_36 + 18]     * x[b_id_6]
                  + a[b_id_36 + 18 + 1] * x[b_id_6 + 1]
                  + a[b_id_36 + 18 + 2] * x[b_id_6 + 2]
                  + a[b_id_36 + 18 + 3] * x[b_id_6 + 3]
                  + a[b_id_36 + 18 + 4] * x[b_id_6 + 4]
                  + a[b_id_36 + 18 + 5] * x[b_id_6 + 5];

  y[b_id_6 + 4] =   a[b_id_36 + 24]     * x[b_id_6]
                  + a[b_id_36 + 24 + 1] * x[b_id_6 + 1]
                  + a[b_id_36 + 24 + 2] * x[b_id_6 + 2]
                  + a[b_id_36 + 24 + 3] * x[b_id_6 + 3]
                  + a[b_id_36 + 24 + 4] * x[b_id_6 + 4]
                  + a[b_id_36 + 24 + 5] * x[b_id_6 + 5];

  y[b_id_6 + 5] =   a[b_id_36 + 30]     * x[b_id_6]
                  + a[b_id_36 + 30 + 1] * x[b_id_6 + 1]
                  + a[b_id_36 + 30 + 2] * x[b_id_6 + 2]
                  + a[b_id_36 + 30 + 3] * x[b_id_6 + 3]
                  + a[b_id_36 + 30 + 4] * x[b_id_6 + 4]
                  + a[b_id_36 + 30 + 5] * x[b_id_6 + 5];

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
 *   a      <-- pointer to block matrixes array (usually matrix extra-diagonal)
 *   x      <-- multipliying vector values
 *   y      --> resulting vector
 *----------------------------------------------------------------------------*/

static inline void
_dense_eb_ax_add(cs_lnum_t        b_i,
                 cs_lnum_t        b_j,
                 cs_lnum_t        b_ij,
                 const cs_lnum_t  b_size[4],
                 const cs_real_t  a[restrict],
                 const cs_real_t  x[restrict],
                 cs_real_t        y[restrict])
{
  cs_lnum_t   ii, jj;

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
 *   da     <-- pointer to coefficients array (usually matrix diagonal)
 *   x      <-- multipliying vector values
 *   y      --> resulting vector
 *   n_elts <-- array size
 *----------------------------------------------------------------------------*/

static inline void
_diag_vec_p_l(const cs_real_t  da[restrict],
              const cs_real_t  x[restrict],
              cs_real_t        y[restrict],
              cs_lnum_t        n_elts)
{
  cs_lnum_t  ii;

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
 *   da     <-- pointer to coefficients array (usually matrix diagonal)
 *   x      <-- multipliying vector values
 *   y      --> resulting vector
 *   n_elts <-- array size
 *   b_size <-- block size, including padding:
 *              b_size[0]: useful block size
 *              b_size[1]: vector block extents
 *              b_size[2]: matrix line extents
 *              b_size[3]: matrix line*column (block) extents
 *----------------------------------------------------------------------------*/

static inline void
_b_diag_vec_p_l(const cs_real_t  da[restrict],
                const cs_real_t  x[restrict],
                cs_real_t        y[restrict],
                cs_lnum_t        n_elts,
                const cs_lnum_t  b_size[4])
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
 *   da     <-- pointer to coefficients array (usually matrix diagonal)
 *   x      <-- multipliying vector values
 *   y      --> resulting vector
 *   n_elts <-- array size
 *----------------------------------------------------------------------------*/

static inline void
_3_3_diag_vec_p_l(const cs_real_t  da[restrict],
                  const cs_real_t  x[restrict],
                  cs_real_t        y[restrict],
                  cs_lnum_t        n_elts)
{
  cs_lnum_t   ii;

  if (da != NULL) {
#   pragma omp parallel for  if(n_elts*3 > CS_THR_MIN)
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
 * Block version of y[i] = da[i].x[i], with da possibly NULL
 *
 * This variant uses a fixed 6x6 block, for better compiler optimization.
 *
 * parameters:
 *   da     <-- pointer to coefficients array (usually matrix diagonal)
 *   x      <-- multipliying vector values
 *   y      --> resulting vector
 *   n_elts <-- array size
 *----------------------------------------------------------------------------*/

static inline void
_6_6_diag_vec_p_l(const cs_real_t  da[restrict],
                  const cs_real_t  x[restrict],
                  cs_real_t        y[restrict],
                  cs_lnum_t        n_elts)
{
  cs_lnum_t   ii;

  if (da != NULL) {
#   pragma omp parallel for  if(n_elts*6 > CS_THR_MIN)
    for (ii = 0; ii < n_elts; ii++)
      _dense_6_6_ax(ii, da, x, y);
  }
  else {
#   pragma omp parallel for  if(n_elts*6 > CS_THR_MIN)
    for (ii = 0; ii < n_elts*6; ii++)
      y[ii] = 0.0;
  }
}

/*----------------------------------------------------------------------------
 * Set values from y[start_id] to y[end_id] to 0.
 *
 * parameters:
 *   y        --> resulting vector
 *   start_id <-- start id in array
 *   end_id   <-- end id in array
 *----------------------------------------------------------------------------*/

static inline void
_zero_range(cs_real_t   y[restrict],
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
_b_zero_range(cs_real_t        y[restrict],
              cs_lnum_t        start_id,
              cs_lnum_t        end_id,
              const cs_lnum_t  b_size[2])
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
_3_3_zero_range(cs_real_t  y[restrict],
                cs_lnum_t  start_id,
                cs_lnum_t  end_id)
{
  cs_lnum_t  ii;

# pragma omp parallel for  if((end_id-start_id)*3 > CS_THR_MIN)
  for (ii = start_id*3; ii < end_id*3; ii++)
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
_6_6_zero_range(cs_real_t  y[restrict],
                cs_lnum_t  start_id,
                cs_lnum_t  end_id)
{
  cs_lnum_t  ii;

# pragma omp parallel for  if((end_id-start_id)*6 > CS_THR_MIN)
  for (ii = start_id*6; ii < end_id*6; ii++)
    y[ii] = 0.0;
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
 *   matrix  <->  pointer to native matrix structure pointer
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
 *   coeff  <->  pointer to native matrix coefficients pointer
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

  cs_matrix_coeff_native_t  *mc = matrix->coeffs;
  const cs_matrix_struct_native_t  *ms = matrix->structure;
  mc->symmetric = symmetric;

  /* Map or copy values */

  if (da != NULL) {

    if (copy) {
      if (mc->_da == NULL || mc->max_db_size < matrix->db_size[3]) {
        BFT_REALLOC(mc->_da, matrix->db_size[3]*ms->n_rows, cs_real_t);
        mc->max_db_size = matrix->db_size[3];
      }
      memcpy(mc->_da, da, matrix->db_size[3]*sizeof(cs_real_t) * ms->n_rows);
      mc->da = mc->_da;
    }
    else
      mc->da = da;

  }
  else {
    mc->da = NULL;
  }

  if (xa != NULL) {

    size_t xa_n_vals = ms->n_edges;
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
}

/*----------------------------------------------------------------------------
 * Release shared native matrix coefficients.
 *
 * parameters:
 *   matrix <-- pointer to matrix structure
 *----------------------------------------------------------------------------*/

static void
_release_coeffs_native(cs_matrix_t  *matrix)
{
  cs_matrix_coeff_native_t  *mc = matrix->coeffs;
  if (mc != NULL) {
    mc->da = NULL;
    mc->xa = NULL;
  }
}

/*----------------------------------------------------------------------------
 * Copy diagonal of native or MSR matrix.
 *
 * parameters:
 *   matrix <-- pointer to matrix structure
 *   da     --> diagonal (pre-allocated, size: n_rows)
 *----------------------------------------------------------------------------*/

static void
_copy_diagonal_separate(const cs_matrix_t  *matrix,
                        cs_real_t           da[restrict])
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
  const cs_lnum_t  n_rows = matrix->n_rows;

  /* Unblocked version */

  if (matrix->db_size[3] == 1) {

    if (_da != NULL) {
#     pragma omp parallel for  if(n_rows > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        da[ii] = _da[ii];
    }
    else {
#     pragma omp parallel for  if(n_rows > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++)
        da[ii] = 0.0;
    }

  }

  /* Blocked version */

  else {

    const cs_lnum_t *db_size = matrix->db_size;

    if (_da != NULL) {
#     pragma omp parallel for  if(n_rows*db_size[0] > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
        for (cs_lnum_t jj = 0; jj < db_size[0]; jj++)
          da[ii*db_size[1] + jj] = _da[ii*db_size[3] + jj*db_size[2] + jj];
      }
    }
    else {
#     pragma omp parallel for  if(n_rows*db_size[1] > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_rows*db_size[1]; ii++)
        da[ii] = 0.0;
    }
  }
}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with native matrix.
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- pointer to matrix structure
 *   x            <-- multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_p_l_native(bool                exclude_diag,
                    const cs_matrix_t  *matrix,
                    const cs_real_t     x[restrict],
                    cs_real_t           y[restrict])
{
  cs_lnum_t  ii, jj, face_id;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;

  const cs_real_t  *restrict xa = mc->xa;

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _diag_vec_p_l(mc->da, x, y, ms->n_rows);
    _zero_range(y, ms->n_rows, ms->n_cols_ext);
  }
  else
    _zero_range(y, 0, ms->n_cols_ext);

  /* Note: parallel and periodic synchronization could be delayed to here */

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    const cs_lnum_2_t *restrict face_cel_p = ms->edges;

    if (mc->symmetric) {

      for (face_id = 0; face_id < ms->n_edges; face_id++) {
        ii = face_cel_p[face_id][0];
        jj = face_cel_p[face_id][1];
        y[ii] += xa[face_id] * x[jj];
        y[jj] += xa[face_id] * x[ii];
      }

    }
    else {

      for (face_id = 0; face_id < ms->n_edges; face_id++) {
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
 *   matrix       <-- pointer to matrix structure
 *   x            <-- multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_b_mat_vec_p_l_native(bool                exclude_diag,
                      const cs_matrix_t  *matrix,
                      const cs_real_t    x[restrict],
                      cs_real_t          y[restrict])
{
  cs_lnum_t  ii, jj, kk, face_id;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;

  const cs_real_t  *restrict xa = mc->xa;
  const cs_lnum_t *db_size = matrix->db_size;

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _b_diag_vec_p_l(mc->da, x, y, ms->n_rows, db_size);
    _b_zero_range(y, ms->n_rows, ms->n_cols_ext, db_size);
  }
  else
    _b_zero_range(y, 0, ms->n_cols_ext, db_size);

  /* Note: parallel and periodic synchronization could be delayed to here */

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    const cs_lnum_2_t *restrict face_cel_p = ms->edges;

    if (mc->symmetric) {

      for (face_id = 0; face_id < ms->n_edges; face_id++) {
        ii = face_cel_p[face_id][0];
        jj = face_cel_p[face_id][1];
        for (kk = 0; kk < db_size[0]; kk++) {
          y[ii*db_size[1] + kk] += xa[face_id] * x[jj*db_size[1] + kk];
          y[jj*db_size[1] + kk] += xa[face_id] * x[ii*db_size[1] + kk];
        }
      }
    }
    else {

      for (face_id = 0; face_id < ms->n_edges; face_id++) {
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
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- pointer to matrix structure
 *   x            <-- multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_bb_mat_vec_p_l_native(bool                exclude_diag,
                       const cs_matrix_t  *matrix,
                       const cs_real_t    x[restrict],
                       cs_real_t          y[restrict])
{
  cs_lnum_t  ii, jj, face_id;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;

  const cs_real_t  *restrict xa = mc->xa;
  const cs_lnum_t *db_size = matrix->db_size;
  const cs_lnum_t *eb_size = matrix->eb_size;

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _b_diag_vec_p_l(mc->da, x, y, ms->n_rows, db_size);
    _b_zero_range(y, ms->n_rows, ms->n_cols_ext, db_size);
  }
  else
    _b_zero_range(y, 0, ms->n_cols_ext, db_size);

  /* Note: parallel and periodic synchronization could be delayed to here */

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    const cs_lnum_2_t *restrict face_cel_p = ms->edges;

    if (mc->symmetric) {

      for (face_id = 0; face_id < ms->n_edges; face_id++) {
        ii = face_cel_p[face_id][0];
        jj = face_cel_p[face_id][1];
        _dense_eb_ax_add(ii, jj, face_id, eb_size, xa, x, y);
        _dense_eb_ax_add(jj, ii, face_id, eb_size, xa, x, y);
      }
    }
    else {

      for (face_id = 0; face_id < ms->n_edges; face_id++) {
        ii = face_cel_p[face_id][0];
        jj = face_cel_p[face_id][1];
        _dense_eb_ax_add(ii, jj, 2*face_id, eb_size, xa, x, y);
        _dense_eb_ax_add(jj, ii, 2*face_id + 1, eb_size, xa, x, y);
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
 *   matrix       <-- pointer to matrix structure
 *   x            <-- multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_3_3_mat_vec_p_l_native(bool                exclude_diag,
                        const cs_matrix_t  *matrix,
                        const cs_real_t     x[restrict],
                        cs_real_t           y[restrict])
{
  cs_lnum_t  ii, jj, kk, face_id;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;

  const cs_real_t  *restrict xa = mc->xa;

  assert(matrix->db_size[0] == 3 && matrix->db_size[3] == 9);

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _3_3_diag_vec_p_l(mc->da, x, y, ms->n_rows);
    _3_3_zero_range(y, ms->n_rows, ms->n_cols_ext);
  }
  else
    _3_3_zero_range(y, 0, ms->n_cols_ext);

  /* Note: parallel and periodic synchronization could be delayed to here */

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    const cs_lnum_2_t *restrict face_cel_p = ms->edges;

    if (mc->symmetric) {

      for (face_id = 0; face_id < ms->n_edges; face_id++) {
        ii = face_cel_p[face_id][0];
        jj = face_cel_p[face_id][1];
        for (kk = 0; kk < 3; kk++) {
          y[ii*3 + kk] += xa[face_id] * x[jj*3 + kk];
          y[jj*3 + kk] += xa[face_id] * x[ii*3 + kk];
        }
      }
    }
    else {

      for (face_id = 0; face_id < ms->n_edges; face_id++) {
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
 * This variant uses a fixed 6x6 block, for better compiler optimization.
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- pointer to matrix structure
 *   x            <-- multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_6_6_mat_vec_p_l_native(bool                exclude_diag,
                        const cs_matrix_t  *matrix,
                        const cs_real_t     x[restrict],
                        cs_real_t           y[restrict])
{
  cs_lnum_t  ii, jj, kk, face_id;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;

  const cs_real_t  *restrict xa = mc->xa;

  assert(matrix->db_size[0] == 3 && matrix->db_size[3] == 9);

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _6_6_diag_vec_p_l(mc->da, x, y, ms->n_rows);
    _6_6_zero_range(y, ms->n_rows, ms->n_cols_ext);
  }
  else
    _6_6_zero_range(y, 0, ms->n_cols_ext);

  /* Note: parallel and periodic synchronization could be delayed to here */

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    const cs_lnum_2_t *restrict face_cel_p = ms->edges;

    if (mc->symmetric) {

      for (face_id = 0; face_id < ms->n_edges; face_id++) {
        ii = face_cel_p[face_id][0];
        jj = face_cel_p[face_id][1];
        for (kk = 0; kk < 6; kk++) {
          y[ii*6 + kk] += xa[face_id] * x[jj*6 + kk];
          y[jj*6 + kk] += xa[face_id] * x[ii*6 + kk];
        }
      }
    }
    else {

      for (face_id = 0; face_id < ms->n_edges; face_id++) {
        ii = face_cel_p[face_id][0];
        jj = face_cel_p[face_id][1];
        for (kk = 0; kk < 6; kk++) {
          y[ii*6 + kk] += xa[2*face_id]     * x[jj*6 + kk];
          y[jj*6 + kk] += xa[2*face_id + 1] * x[ii*6 + kk];
        }
      }

    }

  }

}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with native matrix, blocked version.
 *
 * This variant uses fixed block size variants for common cases.
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- pointer to matrix structure
 *   x            <-- multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_b_mat_vec_p_l_native_fixed(bool                exclude_diag,
                            const cs_matrix_t  *matrix,
                            const cs_real_t     x[restrict],
                            cs_real_t           y[restrict])
{
  if (matrix->db_size[0] == 3 && matrix->db_size[3] == 9)
    _3_3_mat_vec_p_l_native(exclude_diag, matrix, x, y);

  else if (matrix->db_size[0] == 6 && matrix->db_size[3] == 36)
    _6_6_mat_vec_p_l_native(exclude_diag, matrix, x, y);

  else
    _b_mat_vec_p_l_native(exclude_diag, matrix, x, y);
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
                        const cs_real_t     x[restrict],
                        cs_real_t           y[restrict])
{
  const int n_threads = matrix->numbering->n_threads;
  const int n_groups = matrix->numbering->n_groups;
  const cs_lnum_t *group_index = matrix->numbering->group_index;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;
  const cs_real_t  *restrict xa = mc->xa;

  assert(matrix->numbering->type == CS_NUMBERING_THREADS);

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _diag_vec_p_l(mc->da, x, y, ms->n_rows);
    _zero_range(y, ms->n_rows, ms->n_cols_ext);
  }
  else
    _zero_range(y, 0, ms->n_cols_ext);

  /* Note: parallel and periodic synchronization could be delayed to here */

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    const cs_lnum_2_t *restrict face_cel_p = ms->edges;

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
                          const cs_real_t     x[restrict],
                          cs_real_t           y[restrict])
{
  const cs_lnum_t *db_size = matrix->db_size;

  const int n_threads = matrix->numbering->n_threads;
  const int n_groups = matrix->numbering->n_groups;
  const cs_lnum_t *group_index = matrix->numbering->group_index;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;
  const cs_real_t  *restrict xa = mc->xa;

  assert(matrix->numbering->type == CS_NUMBERING_THREADS);

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _b_diag_vec_p_l(mc->da, x, y, ms->n_rows, db_size);
    _b_zero_range(y, ms->n_rows, ms->n_cols_ext, db_size);
  }
  else
    _b_zero_range(y, 0, ms->n_cols_ext, db_size);

  /* Note: parallel and periodic synchronization could be delayed to here */

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    const cs_lnum_2_t *restrict face_cel_p = ms->edges;

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
                               const cs_real_t     x[restrict],
                               cs_real_t           y[restrict])
{
  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;
  const cs_real_t  *restrict xa = mc->xa;

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _diag_vec_p_l(mc->da, x, y, ms->n_rows);
    _zero_range(y, ms->n_rows, ms->n_cols_ext);
  }
  else
    _zero_range(y, 0, ms->n_cols_ext);

  /* Note: parallel and periodic synchronization could be delayed to here */

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    const cs_lnum_2_t *restrict face_cel_p = ms->edges;

    if (mc->symmetric) {

#     pragma omp parallel for
      for (cs_lnum_t face_id = 0; face_id < ms->n_edges; face_id++) {
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
      for (cs_lnum_t face_id = 0; face_id < ms->n_edges; face_id++) {
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
                                 const cs_real_t     x[restrict],
                                 cs_real_t           y[restrict])
{
  const cs_lnum_t *db_size = matrix->db_size;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;
  const cs_real_t  *restrict xa = mc->xa;

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _b_diag_vec_p_l(mc->da, x, y, ms->n_rows, db_size);
    _b_zero_range(y, ms->n_rows, ms->n_cols_ext, db_size);
  }
  else
    _b_zero_range(y, 0, ms->n_cols_ext, db_size);

  /* Note: parallel and periodic synchronization could be delayed to here */

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    const cs_lnum_2_t *restrict face_cel_p = ms->edges;

    if (mc->symmetric) {

#     pragma omp parallel for
      for (cs_lnum_t face_id = 0; face_id < ms->n_edges; face_id++) {
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
      for (cs_lnum_t face_id = 0; face_id < ms->n_edges; face_id++) {
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

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with native matrix.
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- pointer to matrix structure
 *   x            <-- multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_p_l_native_vector(bool                exclude_diag,
                           const cs_matrix_t  *matrix,
                           const cs_real_t     x[restrict],
                           cs_real_t           y[restrict])
{
  cs_lnum_t  ii, jj, face_id;
  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;
  const cs_real_t  *restrict xa = mc->xa;

  assert(matrix->numbering->type == CS_NUMBERING_VECTORIZE);

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _diag_vec_p_l(mc->da, x, y, ms->n_rows);
    _zero_range(y, ms->n_rows, ms->n_cols_ext);
  }
  else
    _zero_range(y, 0, ms->n_cols_ext);

  /* Note: parallel and periodic synchronization could be delayed to here */

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    const cs_lnum_2_t *restrict face_cel_p = ms->edges;

    if (mc->symmetric) {

#     if defined(HAVE_OPENMP_SIMD)
#       pragma omp simd safelen(CS_NUMBERING_SIMD_SIZE)
#     else
#       pragma dir nodep
#       pragma GCC ivdep
#       pragma _NEC ivdep
#     endif
      for (face_id = 0; face_id < ms->n_edges; face_id++) {
        ii = face_cel_p[face_id][0];
        jj = face_cel_p[face_id][1];
        y[ii] += xa[face_id] * x[jj];
        y[jj] += xa[face_id] * x[ii];
      }

    }
    else {

#     if defined(HAVE_OPENMP_SIMD)
#       pragma omp simd safelen(CS_NUMBERING_SIMD_SIZE)
#     else
#       pragma dir nodep
#       pragma GCC ivdep
#       pragma _NEC ivdep
#     endif
      for (face_id = 0; face_id < ms->n_edges; face_id++) {
        ii = face_cel_p[face_id][0];
        jj = face_cel_p[face_id][1];
        y[ii] += xa[2*face_id] * x[jj];
        y[jj] += xa[2*face_id + 1] * x[ii];
      }

    }

  }
}

/*----------------------------------------------------------------------------
 * Destroy a CSR matrix structure.
 *
 * parameters:
 *   matrix  <->  pointer to CSR matrix structure pointer
 *----------------------------------------------------------------------------*/

static void
_destroy_struct_csr(cs_matrix_struct_csr_t  **matrix)
{
  if (matrix != NULL && *matrix !=NULL) {

    cs_matrix_struct_csr_t  *ms = *matrix;

    BFT_FREE(ms->_row_index);

    BFT_FREE(ms->_col_id);

    BFT_FREE(ms);

    *matrix = NULL;

  }
}

/*----------------------------------------------------------------------------
 * Create a CSR matrix structure from a native matrix stucture.
 *
 * Note that the structure created maps global cell numbers to the given
 * existing face -> cell connectivity array, so it must be destroyed before
 * this array (usually the code's global cell numbering) is freed.
 *
 * parameters:
 *   have_diag   <-- indicates if the diagonal is nonzero
 *   n_rows      <-- number of local rows
 *   n_cols_ext  <-- number of local + ghost columns
 *   n_edges     <-- local number of graph edges
 *   edges       <-- edges (symmetric row <-> column) connectivity
 *
 * returns:
 *   pointer to allocated CSR matrix structure.
 *----------------------------------------------------------------------------*/

static cs_matrix_struct_csr_t *
_create_struct_csr(bool                have_diag,
                   cs_lnum_t           n_rows,
                   cs_lnum_t           n_cols_ext,
                   cs_lnum_t           n_edges,
                   const cs_lnum_2_t  *edges)
{
  cs_lnum_t ii, jj, face_id;
  const cs_lnum_t *restrict face_cel_p;

  cs_lnum_t  diag_elts = 1;
  cs_lnum_t  *ccount = NULL;

  cs_matrix_struct_csr_t  *ms;

  /* Allocate and map */

  BFT_MALLOC(ms, 1, cs_matrix_struct_csr_t);

  ms->n_rows = n_rows;
  ms->n_cols_ext = n_cols_ext;

  ms->direct_assembly = true;
  ms->have_diag = have_diag;

  BFT_MALLOC(ms->_row_index, ms->n_rows + 1, cs_lnum_t);
  ms->row_index = NULL;

  /* Count number of nonzero elements per row */

  BFT_MALLOC(ccount, ms->n_rows, cs_lnum_t);

  if (have_diag == false)
    diag_elts = 0;

  for (ii = 0; ii < ms->n_rows; ii++)  /* count starting with diagonal terms */
    ccount[ii] = diag_elts;

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
    ccount[ii] = diag_elts; /* pre-count for diagonal terms */
  }

  /* Build structure */

  BFT_MALLOC(ms->_col_id, (ms->_row_index[ms->n_rows]), cs_lnum_t);
  ms->col_id = NULL;

  if (have_diag == true) {
    for (ii = 0; ii < ms->n_rows; ii++) {    /* diagonal terms */
      ms->_col_id[ms->_row_index[ii]] = ii;
    }
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

    for (ii = 0; ii < ms->n_rows; ii++) {
      cs_lnum_t *col_id = ms->_col_id + ms->_row_index[ii];
      cs_lnum_t n_cols = ms->_row_index[ii+1] - ms->_row_index[ii];
      cs_lnum_t col_id_prev = -1;
      ms->_row_index[ii] = kk;
      for (jj = 0; jj < n_cols; jj++) {
        if (col_id_prev != col_id[jj]) {
          ms->_col_id[kk++] = col_id[jj];
          col_id_prev = col_id[jj];
        }
      }
    }
    ms->_row_index[ms->n_rows] = kk;

    assert(ms->_row_index[ms->n_rows] < tmp_row_index[ms->n_rows]);

    BFT_FREE(tmp_row_index);
    BFT_REALLOC(ms->_col_id, (ms->_row_index[ms->n_rows]), cs_lnum_t);

  }

  ms->row_index = ms->_row_index;
  ms->col_id = ms->_col_id;

  return ms;
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
 *   col_id     <-> pointer to array of colum ids related to the row index
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

  cs_lnum_t  *_row_index = *row_index;
  cs_lnum_t  *_col_id = *col_id;

  /* Allocate and map */

  BFT_MALLOC(ms, 1, cs_matrix_struct_csr_t);

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
 * Destroy CSR matrix coefficients.
 *
 * parameters:
 *   coeff  <->  pointer to CSR matrix coefficients pointer
 *----------------------------------------------------------------------------*/

static void
_destroy_coeff_csr(cs_matrix_coeff_csr_t **coeff)
{
  if (coeff != NULL && *coeff !=NULL) {

    cs_matrix_coeff_csr_t  *mc = *coeff;

    BFT_FREE(mc->_val);
    BFT_FREE(mc->_d_val);

    BFT_FREE(*coeff);

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
 *   matrix           <-> pointer to matrix structure
 *----------------------------------------------------------------------------*/

static void
_zero_coeffs_csr(cs_matrix_t  *matrix)
{
  cs_matrix_coeff_csr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  const cs_lnum_t  n_rows = ms->n_rows;
  const cs_lnum_t *eb_size = matrix->eb_size;

  if (eb_size[0] == 1) {
#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
      const cs_lnum_t  n_cols = ms->row_index[ii+1] - ms->row_index[ii];
      cs_real_t  *m_row = mc->_val + ms->row_index[ii];
      for (cs_lnum_t jj = 0; jj < n_cols; jj++)
        m_row[jj] = 0.0;
    }
  }
  else {
#   pragma omp parallel for  if(n_rows*eb_size[0] > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
      const cs_lnum_t  n_cols = ms->row_index[ii+1] - ms->row_index[ii];
      cs_real_t  *m_row = mc->_val + ms->row_index[ii]*eb_size[3];
      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        for (cs_lnum_t kk = 0; kk < eb_size[3]; kk++)
          m_row[jj*eb_size[3] + kk] = 0.0;
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

  cs_lnum_t  ii, jj;
  cs_matrix_coeff_csr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  if (mc->_val == NULL)
    BFT_MALLOC(mc->_val, ms->row_index[ms->n_rows], cs_real_t);
  mc->val = mc->_val;

  /* Initialize coefficients to zero if assembly is incremental */

  if (ms->direct_assembly == false)
    _zero_coeffs_csr(matrix);

  /* Copy diagonal values */

  if (ms->have_diag == true) {

    if (da != NULL) {
      for (ii = 0; ii < ms->n_rows; ii++) {
        cs_lnum_t kk;
        for (kk = ms->row_index[ii]; ms->col_id[kk] != ii; kk++);
        mc->_val[kk] = da[ii];
      }
    }
    else {
      for (ii = 0; ii < ms->n_rows; ii++) {
        cs_lnum_t kk;
        for (kk = ms->row_index[ii]; ms->col_id[kk] != ii; kk++);
        mc->_val[kk] = 0.0;
      }
    }

  }

  /* Mark diagonal values as not queried (mc->_d_val not changed) */

  mc->d_val = NULL;

  /* Copy extra-diagonal values */

  if (edges != NULL) {

    if (xa != NULL) {

      if (ms->direct_assembly == true)
        _set_xa_coeffs_csr_direct(matrix, symmetric, n_edges, edges, xa);
      else
        _set_xa_coeffs_csr_increment(matrix, symmetric, n_edges, edges, xa);

    }
    else { /* if (xa == NULL) */

      for (ii = 0; ii < ms->n_rows; ii++) {
        const cs_lnum_t  *restrict col_id = ms->col_id + ms->row_index[ii];
        cs_real_t  *m_row = mc->_val + ms->row_index[ii];
        cs_lnum_t  n_cols = ms->row_index[ii+1] - ms->row_index[ii];

        for (jj = 0; jj < n_cols; jj++) {
          if (col_id[jj] != ii)
            m_row[jj] = 0.0;
        }

      }

    }

  } /* (matrix->edges != NULL) */

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
 *   d_vals_transfer  <-- diagonal values whose ownership is trasferred
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

  if (matrix->db_size[0] > 1 || matrix->eb_size[0] > 1)
    bft_error
      (__FILE__, __LINE__, 0,
       "%s:\n"
       "  case with diagonal block size %d en extradiagonal block size %d\n"
       "  not implemented.\n",
       __func__, matrix->db_size[0], matrix->eb_size[0]);

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
    BFT_MALLOC(mc->_val, ms->row_index[ms->n_rows], cs_real_t);

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
    _zero_coeffs_csr(matrix);

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

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with CSR matrix.
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- pointer to matrix structure
 *   x            <-- multipliying vector values
 *   y            --> resulting vector
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

      const cs_lnum_t *restrict col_id = ms->col_id + ms->row_index[ii];
      const cs_real_t *restrict m_row = mc->val + ms->row_index[ii];
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

      const cs_lnum_t *restrict col_id = ms->col_id + ms->row_index[ii];
      const cs_real_t *restrict m_row = mc->val + ms->row_index[ii];
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
 * Create a symmetric CSR matrix structure from a native matrix stucture.
 *
 * Note that the structure created maps global cell numbers to the given
 * existing face -> cell connectivity array, so it must be destroyed before
 * this array (usually the code's global cell numbering) is freed.
 *
 * parameters:
 *   have_diag   <-- indicates if the diagonal is nonzero
 *                   (forced to true for symmetric variant)
 *   n_rows      <-- number of local rows
 *   n_cols_ext  <-- number of local + ghost columns
 *   n_edges     <-- local number of graph edges
 *   edges       <-- edges (symmetric row <-> column) connectivity
 *
 * returns:
 *   pointer to allocated CSR matrix structure.
 *----------------------------------------------------------------------------*/

static cs_matrix_struct_csr_sym_t *
_create_struct_csr_sym(bool                have_diag,
                       cs_lnum_t           n_rows,
                       cs_lnum_t           n_cols_ext,
                       cs_lnum_t           n_edges,
                       const cs_lnum_2_t  *edges)
{
  cs_lnum_t ii, jj, face_id;
  const cs_lnum_t *restrict edges_p;

  cs_lnum_t  diag_elts = 1;
  cs_lnum_t  *ccount = NULL;

  cs_matrix_struct_csr_sym_t  *ms;

  /* Allocate and map */

  BFT_MALLOC(ms, 1, cs_matrix_struct_csr_sym_t);

  ms->n_rows = n_rows;
  ms->n_cols = n_cols_ext;

  ms->have_diag = have_diag;
  ms->direct_assembly = true;

  /* Allocate index on size n_cols rather than n_rows, so
     as to prolong index, and make local matrix really "square".
     This is not important for the local code, but seems expected
     of external libraries such as MKL (i.e. not doing this
     causes crashes usign the MKL routines). */

  BFT_MALLOC(ms->row_index, ms->n_cols + 1, cs_lnum_t);
  ms->row_index = ms->row_index;

  /* Count number of nonzero elements per row */

  BFT_MALLOC(ccount, ms->n_cols, cs_lnum_t);

  if (have_diag == false)
    diag_elts = 0;

  for (ii = 0; ii < ms->n_rows; ii++)  /* count starting with diagonal terms */
    ccount[ii] = diag_elts;

  if (edges != NULL) {

    edges_p = (const cs_lnum_t *restrict)edges;

    for (face_id = 0; face_id < n_edges; face_id++) {
      ii = *edges_p++;
      jj = *edges_p++;
      if (ii < jj)
        ccount[ii] += 1;
      else
        ccount[jj] += 1;
    }

  } /* if (edges != NULL) */

  ms->row_index[0] = 0;
  for (ii = 0; ii < ms->n_rows; ii++) {
    ms->row_index[ii+1] = ms->row_index[ii] + ccount[ii];
    ccount[ii] = diag_elts; /* pre-count for diagonal terms */
  }

  /* Build structure */

  BFT_MALLOC(ms->col_id, (ms->row_index[ms->n_rows]), cs_lnum_t);
  ms->col_id = ms->col_id;

  if (have_diag == true) {
    for (ii = 0; ii < ms->n_rows; ii++) {    /* diagonal terms */
      ms->col_id[ms->row_index[ii]] = ii;
    }
  }

  if (edges != NULL) {                   /* non-diagonal terms */

    edges_p = (const cs_lnum_t *restrict)edges;

    for (face_id = 0; face_id < n_edges; face_id++) {
      ii = *edges_p++;
      jj = *edges_p++;
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

  /* Prolong index in case of use by external library, such as MKL. */

  for (ii = ms->n_rows; ii < ms->n_cols; ii++)
    ms->row_index[ii+1] = ms->row_index[ms->n_rows];

  return ms;
}

/*----------------------------------------------------------------------------
 * Destroy symmetric CSR matrix structure.
 *
 * parameters:
 *   matrix  <->  pointer to CSR matrix structure pointer
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
 *   coeff  <->  pointer to CSR matrix coefficients pointer
 *----------------------------------------------------------------------------*/

static void
_destroy_coeff_csr_sym(cs_matrix_coeff_csr_sym_t  **coeff)
{
  if (coeff != NULL && *coeff !=NULL) {

    cs_matrix_coeff_csr_sym_t  *mc = *coeff;

    BFT_FREE(mc->val);
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
 *   matrix    <-- pointer to matrix structure
 *   n_edges   <-- local number of graph edges
 *   edges     <-- edges (symmetric row <-> column) connectivity
 *   xa        <-- extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_xa_coeffs_csr_sym_direct(cs_matrix_t        *matrix,
                              cs_lnum_t           n_edges,
                              const cs_lnum_2_t  *restrict edges,
                              const cs_real_t    *restrict xa)
{
  cs_matrix_coeff_csr_sym_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_sym_t  *ms = matrix->structure;

  /* Copy extra-diagonal values */

  assert(edges != NULL);

# pragma omp parallel for  if(n_edges > CS_THR_MIN)
  for (cs_lnum_t face_id = 0; face_id < n_edges; face_id++) {
    cs_lnum_t ii = edges[face_id][0];
    cs_lnum_t jj = edges[face_id][1];
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
 *   matrix    <-- pointer to matrix structure
 *   n_edges   <-- local number of graph edges
 *   edges     <-- edges (symmetric row <-> column) connectivity
 *   xa        <-- extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_xa_coeffs_csr_sym_increment(cs_matrix_t        *matrix,
                                 cs_lnum_t           n_edges,
                                 const cs_lnum_2_t  *restrict edges,
                                 const cs_real_t    *restrict xa)
{
  cs_lnum_t  ii, jj, face_id;
  cs_matrix_coeff_csr_sym_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_sym_t  *ms = matrix->structure;
  const cs_lnum_t *restrict edges_p
    = (const cs_lnum_t *restrict)(edges);

  /* Copy extra-diagonal values */

  assert(edges != NULL);

  for (face_id = 0; face_id < n_edges; face_id++) {
    cs_lnum_t kk;
    ii = *edges_p++;
    jj = *edges_p++;
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
 *   matrix    <-> pointer to matrix structure
 *   symmetric <-- indicates if extradiagonal values are symmetric (true)
 *   copy      <-- indicates if coefficients should be copied
 *   n_edges   <-- local number of graph edges
 *   edges     <-- edges (symmetric row <-> column) connectivity
 *   da        <-- diagonal values (NULL if all zero)
 *   xa        <-- extradiagonal values (NULL if all zero)
 *----------------------------------------------------------------------------*/

static void
_set_coeffs_csr_sym(cs_matrix_t        *matrix,
                    bool                symmetric,
                    bool                copy,
                    cs_lnum_t           n_edges,
                    const cs_lnum_2_t  *restrict edges,
                    const cs_real_t    *restrict da,
                    const cs_real_t    *restrict xa)
{
  CS_UNUSED(copy);

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

  if (edges != NULL) {

    if (xa != NULL) {

      if (symmetric == false)
        bft_error(__FILE__, __LINE__, 0,
                  _("Assigning non-symmetric matrix coefficients to a matrix\n"
                    "in a symmetric CSR format."));

      if (ms->direct_assembly == true)
        _set_xa_coeffs_csr_sym_direct(matrix, n_edges, edges, xa);
      else
        _set_xa_coeffs_csr_sym_increment(matrix, n_edges, edges, xa);

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

  } /* (edges != NULL) */

}

/*----------------------------------------------------------------------------
 * Release shared symmetric CSR matrix coefficients.
 *
 * parameters:
 *   matrix <-- pointer to matrix structure
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
 *   matrix <-- pointer to matrix structure
 *   da     --> diagonal (pre-allocated, size: n_rows)
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
 *   matrix       <-- pointer to matrix structure
 *   x            <-- multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_p_l_csr_sym(bool                 exclude_diag,
                     const cs_matrix_t   *matrix,
                     const cs_real_t      x[restrict],
                     cs_real_t            y[restrict])
{
  cs_lnum_t  ii, jj, n_cols;
  cs_lnum_t  *restrict col_id;
  cs_real_t  *restrict m_row;

  const cs_matrix_struct_csr_sym_t  *ms = matrix->structure;
  const cs_matrix_coeff_csr_sym_t  *mc = matrix->coeffs;
  cs_lnum_t  n_rows = ms->n_rows;

  cs_lnum_t jj_start = 0;
  cs_lnum_t sym_jj_start = 0;

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

  mc->max_db_size = 0;
  mc->max_eb_size = 0;

  mc->d_val = NULL;
  mc->x_val = NULL;

  mc->_d_val = NULL;
  mc->_x_val = NULL;

  return mc;
}

/*----------------------------------------------------------------------------
 * Destroy MSR matrix coefficients.
 *
 * parameters:
 *   coeff  <->  pointer to MSR matrix coefficients pointer
 *----------------------------------------------------------------------------*/

static void
_destroy_coeff_msr(cs_matrix_coeff_msr_t  **coeff)
{
  if (coeff != NULL && *coeff !=NULL) {

    cs_matrix_coeff_msr_t  *mc = *coeff;

    BFT_FREE(mc->_x_val);

    BFT_FREE(mc->_d_val);

    BFT_FREE(*coeff);

  }
}

/*----------------------------------------------------------------------------
 * Set MSR matrix extradiagonal coefficients to zero.
 *
 * The coefficients should already be allocated.
 *
 * Use of this function is preferrable to a simple loop, as its
 * threading behavior should be consistent with SpMW in NUMA cases.
 *
 * parameters:
 *   matrix           <-> pointer to matrix structure
 *----------------------------------------------------------------------------*/

static void
_zero_x_coeffs_msr(cs_matrix_t  *matrix)
{
  cs_matrix_coeff_msr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  const cs_lnum_t  n_rows = ms->n_rows;
  const cs_lnum_t *eb_size = matrix->eb_size;

  if (eb_size[0] == 1) {
#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
      const cs_lnum_t  n_cols = ms->row_index[ii+1] - ms->row_index[ii];
      cs_real_t  *m_row = mc->_x_val + ms->row_index[ii];
      for (cs_lnum_t jj = 0; jj < n_cols; jj++)
        m_row[jj] = 0.0;
    }
  }
  else {
#   pragma omp parallel for  if(n_rows*eb_size[0] > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
      const cs_lnum_t  n_cols = ms->row_index[ii+1] - ms->row_index[ii];
      cs_real_t  *m_row = mc->_x_val + ms->row_index[ii]*eb_size[3];
      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        for (cs_lnum_t kk = 0; kk < eb_size[3]; kk++)
          m_row[jj*eb_size[3] + kk] = 0.0;
      }
    }
  }
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
_set_xa_coeffs_msr_direct(cs_matrix_t        *matrix,
                          bool                symmetric,
                          cs_lnum_t           n_edges,
                          const cs_lnum_2_t  *edges,
                          const cs_real_t    *restrict xa)
{
  cs_lnum_t  ii, jj, face_id;
  cs_matrix_coeff_msr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  /* Copy extra-diagonal values */

  assert(edges != NULL);

  if (symmetric == false) {

    const cs_lnum_t *restrict edges_p
      = (const cs_lnum_t *restrict)(edges);

    for (face_id = 0; face_id < n_edges; face_id++) {
      cs_lnum_t kk, ll;
      ii = *edges_p++;
      jj = *edges_p++;
      if (ii < ms->n_rows) {
        for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
        mc->_x_val[kk] = xa[2*face_id];
      }
      if (jj < ms->n_rows) {
        for (ll = ms->row_index[jj]; ms->col_id[ll] != ii; ll++);
        mc->_x_val[ll] = xa[2*face_id + 1];
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
      if (ii < ms->n_rows) {
        for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
        mc->_x_val[kk] = xa[face_id];
      }
      if (jj < ms->n_rows) {
        for (ll = ms->row_index[jj]; ms->col_id[ll] != ii; ll++);
        mc->_x_val[ll] = xa[face_id];
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
_set_xa_coeffs_msr_increment(cs_matrix_t        *matrix,
                             bool                symmetric,
                             cs_lnum_t           n_edges,
                             const cs_lnum_2_t  *edges,
                             const cs_real_t    *restrict xa)
{
  cs_lnum_t  ii, jj, face_id;
  cs_matrix_coeff_msr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  /* Copy extra-diagonal values */

  assert(edges != NULL);

  if (symmetric == false) {

    const cs_lnum_t *restrict edges_p
      = (const cs_lnum_t *restrict)(edges);

    for (face_id = 0; face_id < n_edges; face_id++) {
      cs_lnum_t kk, ll;
      ii = *edges_p++;
      jj = *edges_p++;
      if (ii < ms->n_rows) {
        for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
        mc->_x_val[kk] += xa[2*face_id];
      }
      if (jj < ms->n_rows) {
        for (ll = ms->row_index[jj]; ms->col_id[ll] != ii; ll++);
        mc->_x_val[ll] += xa[2*face_id + 1];
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
      if (ii < ms->n_rows) {
        for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
        mc->_x_val[kk] += xa[face_id];
      }
      if (jj < ms->n_rows) {
        for (ll = ms->row_index[jj]; ms->col_id[ll] != ii; ll++);
        mc->_x_val[ll] += xa[face_id];
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
_map_or_copy_da_coeffs_msr(cs_matrix_t      *matrix,
                           bool              copy,
                           const cs_real_t  *restrict da)
{
  cs_matrix_coeff_msr_t  *mc = matrix->coeffs;

  const cs_lnum_t n_rows = matrix->n_rows;
  const cs_lnum_t *db_size = matrix->db_size;

  /* Map or copy diagonal values */

  if (da != NULL) {

    if (copy) {
      if (mc->_d_val == NULL || mc->max_db_size < db_size[3]) {
        BFT_REALLOC(mc->_d_val, db_size[3]*n_rows, cs_real_t);
        mc->max_db_size = db_size[3];
      }
#     pragma omp parallel for  if(n_rows*db_size[0] > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
        for (cs_lnum_t jj = 0; jj < db_size[3]; jj++)
          mc->_d_val[ii*db_size[3] + jj] = da[ii*db_size[3] + jj];
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
_map_or_copy_xa_coeffs_msr(cs_matrix_t      *matrix,
                           bool              copy,
                           const cs_real_t  *restrict xa)
{
  cs_matrix_coeff_msr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_lnum_t n_rows = matrix->n_rows;
  const cs_lnum_t *eb_size = matrix->eb_size;

  if (xa == NULL || copy) {

    /* Ensure allocation */
    if (mc->_x_val == NULL || mc->max_eb_size < eb_size[3]) {
      BFT_REALLOC(mc->_d_val,
                  eb_size[3]*ms->row_index[ms->n_rows],
                  cs_real_t);
      mc->max_eb_size = eb_size[3];
    }
    mc->x_val = mc->_x_val;

    /* zero if required */
    if (xa == NULL)
      _zero_x_coeffs_msr(matrix);

  }

  /* Map or copy extradiagonal values (we could use memcpy, but prefer
     to have a similar threading behavior to SpMV for NUMA performance) */

  if (xa != NULL) {

    if (copy) {
      if (eb_size[0] == 1) {
#       pragma omp parallel for  if(n_rows > CS_THR_MIN)
        for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
          const cs_lnum_t  n_cols = ms->row_index[ii+1] - ms->row_index[ii];
          const cs_real_t  *s_row = xa + ms->row_index[ii];
          cs_real_t  *m_row = mc->_x_val + ms->row_index[ii];
          for (cs_lnum_t jj = 0; jj < n_cols; jj++)
            m_row[jj] = s_row[jj];
        }
      }
      else {
#       pragma omp parallel for  if(n_rows*eb_size[0] > CS_THR_MIN)
        for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
          const cs_lnum_t  n_cols = ms->row_index[ii+1] - ms->row_index[ii];
          const cs_real_t  *s_row = xa + ms->row_index[ii]*eb_size[3];
          cs_real_t  *m_row = mc->_x_val + ms->row_index[ii]*eb_size[3];
          for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
            for (cs_lnum_t kk = 0; kk < eb_size[3]; kk++)
              m_row[jj*eb_size[3] + kk] = s_row[jj*eb_size[3] + kk];
          }
        }
      }
    }

    else
      mc->x_val = xa;

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
  cs_matrix_coeff_msr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  /* Map or copy diagonal values */

  _map_or_copy_da_coeffs_msr(matrix, copy, da);

  /* Extradiagonal values */

  if (mc->_x_val == NULL)
    BFT_MALLOC(mc->_x_val, ms->row_index[ms->n_rows], cs_real_t);
  mc->x_val = mc->_x_val;

  /* Copy extra-diagonal values if assembly is direct */

  if (ms->direct_assembly)
    _set_xa_coeffs_msr_direct(matrix, symmetric, n_edges, edges, xa);

  /* Initialize coefficients to zero if assembly is incremental */

  else {
    _map_or_copy_xa_coeffs_msr(matrix, true, NULL);
    if (xa != NULL)
      _set_xa_coeffs_msr_increment(matrix, symmetric, n_edges, edges, xa);
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

  cs_matrix_coeff_msr_t  *mc = matrix->coeffs;

  bool d_transferred = false, x_transferred = false;

  /* TODO: we should use metadata or check that the row_index and
     column id values are consistent, which should be true as long
     as columns are ordered in an identical manner */

  if (d_vals_transfer != NULL) {
    if (*d_vals_transfer != NULL) {
      mc->max_db_size = matrix->db_size[0];
      if (mc->_d_val != *d_vals_transfer) {
        BFT_FREE(mc->_d_val);
        mc->_d_val = *d_vals_transfer;
      }
      mc->d_val = mc->_d_val;
      *d_vals_transfer = NULL;
      d_transferred = true;
    }
  }

  if (x_vals_transfer != NULL) {
    if (*x_vals_transfer != NULL) {
      mc->max_db_size = matrix->db_size[0];
      BFT_FREE(mc->_x_val);
      mc->_x_val = *x_vals_transfer;
      mc->x_val = mc->_x_val;
      *x_vals_transfer = NULL;
      x_transferred = true;
    }
  }

  if (d_transferred == false)
    _map_or_copy_da_coeffs_msr(matrix, copy, d_vals);

  if (x_transferred == false)
    _map_or_copy_xa_coeffs_msr(matrix, copy, x_vals);

  /* Now free transferred arrays */

  if (d_vals_transfer != NULL)
    BFT_FREE(*d_vals_transfer);
  if (x_vals_transfer != NULL)
    BFT_FREE(*x_vals_transfer);
}

/*----------------------------------------------------------------------------
 * Release shared MSR matrix coefficients.
 *
 * parameters:
 *   matrix <-- pointer to matrix structure
 *----------------------------------------------------------------------------*/

static void
_release_coeffs_msr(cs_matrix_t  *matrix)
{
  cs_matrix_coeff_msr_t  *mc = matrix->coeffs;
  if (mc !=NULL) {
    /* Unmap shared values */
    mc->d_val = NULL;
    mc->x_val = NULL;
  }
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
 * \param[in]       db size   optional diagonal block sizes
 * \param[in]       eb size   optional extra-diagonal block sizes
 */
/*----------------------------------------------------------------------------*/

static void
_csr_assembler_values_init(void             *matrix_p,
                           const cs_lnum_t   db_size[4],
                           const cs_lnum_t   eb_size[4])
{
  CS_UNUSED(db_size);

  cs_matrix_t  *matrix = (cs_matrix_t *)matrix_p;

  cs_matrix_coeff_csr_t  *mc = matrix->coeffs;

  const cs_lnum_t n_rows = matrix->n_rows;
  cs_lnum_t e_stride = 1;
  if (eb_size != NULL)
    e_stride = eb_size[3];

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  /* Initialize diagonal values */

  BFT_REALLOC(mc->_val, e_stride*ms->row_index[ms->n_rows], cs_real_t);
  mc->val = mc->_val;

# pragma omp parallel for  if(n_rows*db_size[0] > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
    cs_lnum_t n_s_cols = (ms->row_index[ii+1] - ms->row_index[ii])*e_stride;
    cs_lnum_t displ = ms->row_index[ii]*e_stride;
    for (cs_lnum_t jj = 0; jj < n_s_cols; jj++)
      mc->_val[displ + jj] = 0;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for addition to CSR matrix coefficients using
 *        local row ids and column indexes.
 *
 * Values whose associated row index is negative should be ignored;
 * Values whose column index is -1 are assumed to be assigned to a
 * separately stored diagonal. Other indexes shoudl be valid.
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
 * \param[in]       val       pointer to values (size: n*stride)
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

  const cs_lnum_t n_rows = matrix->n_rows;
  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  if (stride == 1) {

    /* Copy instead of test for OpenMP to avoid outlining for small sets */

    if (n_rows*stride <= CS_THR_MIN) {
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
#     pragma omp parallel for  if(n_rows*stride > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n; ii++) {
        if (row_id[ii] < 0)
          continue;
        else {
          cs_lnum_t r_id = row_id[ii];
          mc->_val[ms->row_index[r_id] + col_idx[ii]] += vals[ii];
        }
      }
    }
  }

  else { /* if (stride > 1) */

    /* Copy instead of test for OpenMP to avoid outlining for small sets */

    if (n_rows*stride <= CS_THR_MIN) {
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
#     pragma omp parallel for  if(n_rows*stride > CS_THR_MIN)
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
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function for initialization of MSR matrix coefficients using
 *        local row ids and column indexes.
 *
 * \warning  The matrix pointer must point to valid data when the selection
 *           function is called, so the life cycle of the data pointed to
 *           should be at least as long as that of the assembler values
 *           structure.
 *
 * \param[in, out]  matrix_p  untyped pointer to matrix description structure
 * \param[in]       db size   optional diagonal block sizes
 * \param[in]       eb size   optional extra-diagonal block sizes
 */
/*----------------------------------------------------------------------------*/

static void
_msr_assembler_values_init(void             *matrix_p,
                           const cs_lnum_t   db_size[4],
                           const cs_lnum_t   eb_size[4])
{
  cs_matrix_t  *matrix = (cs_matrix_t *)matrix_p;

  cs_matrix_coeff_msr_t  *mc = matrix->coeffs;

  const cs_lnum_t n_rows = matrix->n_rows;

  cs_lnum_t d_stride = 1;
  if (db_size != NULL)
    d_stride = db_size[3];
  cs_lnum_t e_stride = 1;
  if (eb_size != NULL)
    e_stride = eb_size[3];

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  /* Initialize diagonal values */

  BFT_REALLOC(mc->_d_val, d_stride*n_rows, cs_real_t);
  mc->d_val = mc->_d_val;
  mc->max_db_size = d_stride;

  BFT_REALLOC(mc->_x_val, e_stride*ms->row_index[ms->n_rows], cs_real_t);
  mc->x_val = mc->_x_val;
  mc->max_eb_size = e_stride;

# pragma omp parallel for  if(n_rows*db_size[0] > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_rows; ii++) {
    for (cs_lnum_t jj = 0; jj < d_stride; jj++)
      mc->_d_val[ii*d_stride + jj] = 0;
    cs_lnum_t n_s_cols = (ms->row_index[ii+1] - ms->row_index[ii])*e_stride;
    cs_lnum_t displ = ms->row_index[ii]*e_stride;
    for (cs_lnum_t jj = 0; jj < n_s_cols; jj++)
      mc->_x_val[displ + jj] = 0;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for addition to MSR matrix coefficients using
 *        local row ids and column indexes.
 *
 * Values whose associated row index is negative should be ignored;
 * Values whose column index is -1 are assumed to be assigned to a
 * separately stored diagonal. Other indexes shoudl be valid.
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
 * \param[in]       val       pointer to values (size: n*stride)
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

  cs_matrix_coeff_msr_t  *mc = matrix->coeffs;

  const cs_lnum_t n_rows = matrix->n_rows;
  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  if (stride == 1) {

    /* Copy instead of test for OpenMP to avoid outlining for small sets */

    if (n_rows*stride <= CS_THR_MIN) {
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
          mc->_x_val[ms->row_index[r_id] + col_idx[ii]] += vals[ii];
        }
      }
    }

    else {
#     pragma omp parallel for  if(n_rows*stride > CS_THR_MIN)
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
          mc->_x_val[ms->row_index[r_id] + col_idx[ii]] += vals[ii];
        }
      }
    }
  }

  else { /* if (stride > 1) */

    /* Copy instead of test for OpenMP to avoid outlining for small sets */

    if (n_rows*stride <= CS_THR_MIN) {
      for (cs_lnum_t ii = 0; ii < n; ii++) {
        cs_lnum_t r_id = row_id[ii];
        if (r_id < 0)
          continue;
        if (col_idx[ii] < 0) {
          for (cs_lnum_t jj = 0; jj < stride; jj++)
            mc->_d_val[r_id*stride + jj] += vals[ii*stride + jj];
        }
        else {
          cs_lnum_t displ = (ms->row_index[r_id] + col_idx[ii])*stride;
          for (cs_lnum_t jj = 0; jj < stride; jj++)
            mc->_x_val[displ + jj] += vals[ii*stride + jj];
        }
      }
    }

    else {
#     pragma omp parallel for  if(n_rows*stride > CS_THR_MIN)
      for (cs_lnum_t ii = 0; ii < n; ii++) {
        cs_lnum_t r_id = row_id[ii];
        if (r_id < 0)
          continue;
        if (col_idx[ii] < 0) {
          for (cs_lnum_t jj = 0; jj < stride; jj++)
            mc->_d_val[r_id*stride + jj] += vals[ii*stride + jj];
        }
        else {
          cs_lnum_t displ = (ms->row_index[r_id] + col_idx[ii])*stride;
          for (cs_lnum_t jj = 0; jj < stride; jj++)
            mc->_x_val[displ + jj] += vals[ii*stride + jj];
        }
      }
    }
  }

}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with MSR matrix.
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- pointer to matrix structure
 *   x            <-- multipliying vector values
 *   y            --> resulting vector
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

      const cs_lnum_t *restrict col_id = ms->col_id + ms->row_index[ii];
      const cs_real_t *restrict m_row = mc->x_val + ms->row_index[ii];
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

      const cs_lnum_t *restrict col_id = ms->col_id + ms->row_index[ii];
      const cs_real_t *restrict m_row = mc->x_val + ms->row_index[ii];
      cs_lnum_t n_cols = ms->row_index[ii+1] - ms->row_index[ii];
      cs_real_t sii = 0.0;

      for (cs_lnum_t jj = 0; jj < n_cols; jj++)
        sii += (m_row[jj]*x[col_id[jj]]);

      y[ii] = sii;

    }
  }

}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with MSR matrix.
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- pointer to matrix structure
 *   x            <-- multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_p_l_msr_omp_sched(bool                exclude_diag,
                           const cs_matrix_t  *matrix,
                           const cs_real_t    *restrict x,
                           cs_real_t          *restrict y)
{
  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_msr_t  *mc = matrix->coeffs;
  cs_lnum_t  n_rows = ms->n_rows;

  /* Standard case */

  if (!exclude_diag && mc->d_val != NULL) {

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {

      cs_lnum_t n_s_rows = cs_align(n_rows * 0.9, CS_CL);
      if (n_s_rows > n_rows)
        n_s_rows = n_rows;

#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_s_rows; ii++) {

        const cs_lnum_t *restrict col_id = ms->col_id + ms->row_index[ii];
        const cs_real_t *restrict m_row = mc->x_val + ms->row_index[ii];
        cs_lnum_t n_cols = ms->row_index[ii+1] - ms->row_index[ii];
        cs_real_t sii = 0.0;

        for (cs_lnum_t jj = 0; jj < n_cols; jj++)
          sii += (m_row[jj]*x[col_id[jj]]);

        y[ii] = sii + mc->d_val[ii]*x[ii];

      }

#     pragma omp for schedule(dynamic, CS_CL)
      for (cs_lnum_t ii = n_s_rows; ii < n_rows; ii++) {

        const cs_lnum_t *restrict col_id = ms->col_id + ms->row_index[ii];
        const cs_real_t *restrict m_row = mc->x_val + ms->row_index[ii];
        cs_lnum_t n_cols = ms->row_index[ii+1] - ms->row_index[ii];
        cs_real_t sii = 0.0;

        for (cs_lnum_t jj = 0; jj < n_cols; jj++)
          sii += (m_row[jj]*x[col_id[jj]]);

        y[ii] = sii + mc->d_val[ii]*x[ii];

      }

    }

  }

  /* Exclude diagonal */

  else {

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {

      cs_lnum_t n_s_rows = cs_align(n_rows * 0.9, CS_CL);
      if (n_s_rows > n_rows)
        n_s_rows = n_rows;

#     pragma omp for
      for (cs_lnum_t ii = 0; ii < n_s_rows; ii++) {

        const cs_lnum_t *restrict col_id = ms->col_id + ms->row_index[ii];
        const cs_real_t *restrict m_row = mc->x_val + ms->row_index[ii];
        cs_lnum_t n_cols = ms->row_index[ii+1] - ms->row_index[ii];
        cs_real_t sii = 0.0;

        for (cs_lnum_t jj = 0; jj < n_cols; jj++)
          sii += (m_row[jj]*x[col_id[jj]]);

        y[ii] = sii;

      }

#     pragma omp for schedule(dynamic, CS_CL)
      for (cs_lnum_t ii = n_s_rows; ii < n_rows; ii++) {

        const cs_lnum_t *restrict col_id = ms->col_id + ms->row_index[ii];
        const cs_real_t *restrict m_row = mc->x_val + ms->row_index[ii];
        cs_lnum_t n_cols = ms->row_index[ii+1] - ms->row_index[ii];
        cs_real_t sii = 0.0;

        for (cs_lnum_t jj = 0; jj < n_cols; jj++)
          sii += (m_row[jj]*x[col_id[jj]]);

        y[ii] = sii;

      }

    }

  }

}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with MSR matrix, blocked version.
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- pointer to matrix structure
 *   x            <-- multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_b_mat_vec_p_l_msr_generic(bool                exclude_diag,
                           const cs_matrix_t  *matrix,
                           const cs_real_t     x[restrict],
                           cs_real_t           y[restrict])
{
  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_msr_t  *mc = matrix->coeffs;
  const cs_lnum_t  n_rows = ms->n_rows;
  const cs_lnum_t *db_size = matrix->db_size;

  /* Standard case */

  if (!exclude_diag && mc->d_val != NULL) {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      const cs_lnum_t *restrict col_id = ms->col_id + ms->row_index[ii];
      const cs_real_t *restrict m_row = mc->x_val + ms->row_index[ii];
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

      const cs_lnum_t *restrict col_id = ms->col_id + ms->row_index[ii];
      const cs_real_t *restrict m_row = mc->x_val + ms->row_index[ii];
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
 * Local matrix.vector product y = A.x with MSR matrix, 3x3 blocked version.
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- pointer to matrix structure
 *   x            <-- multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_3_3_mat_vec_p_l_msr(bool                exclude_diag,
                     const cs_matrix_t  *matrix,
                     const cs_real_t    *restrict x,
                     cs_real_t          *restrict y)
{
  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_msr_t  *mc = matrix->coeffs;
  const cs_lnum_t  n_rows = ms->n_rows;

  assert(matrix->db_size[0] == 3 && matrix->db_size[3] == 9);

  /* Standard case */

  if (!exclude_diag && mc->d_val != NULL) {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      const cs_lnum_t *restrict col_id = ms->col_id + ms->row_index[ii];
      const cs_real_t *restrict m_row = mc->x_val + ms->row_index[ii];
      cs_lnum_t n_cols = ms->row_index[ii+1] - ms->row_index[ii];

      _dense_3_3_ax(ii, mc->d_val, x, y);

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        for (cs_lnum_t kk = 0; kk < 3; kk++)
          y[ii*3 + kk] += (m_row[jj]*x[col_id[jj]*3 + kk]);
      }

    }

  }

  /* Exclude diagonal */

  else {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      const cs_lnum_t *restrict col_id = ms->col_id + ms->row_index[ii];
      const cs_real_t *restrict m_row = mc->x_val + ms->row_index[ii];
      cs_lnum_t n_cols = ms->row_index[ii+1] - ms->row_index[ii];

      for (cs_lnum_t kk = 0; kk < 3; kk++)
        y[ii*3 + kk] = 0.;

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        for (cs_lnum_t kk = 0; kk < 3; kk++)
          y[ii*3 + kk] += (m_row[jj]*x[col_id[jj]*3 + kk]);
      }

    }
  }

}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with MSR matrix, 6x6 blocked version.
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- pointer to matrix structure
 *   x            <-- multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_6_6_mat_vec_p_l_msr(bool                exclude_diag,
                     const cs_matrix_t  *matrix,
                     const cs_real_t     x[restrict],
                     cs_real_t           y[restrict])
{
  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_msr_t  *mc = matrix->coeffs;
  const cs_lnum_t  n_rows = ms->n_rows;

  assert(matrix->db_size[0] == 6 && matrix->db_size[3] == 36);

  /* Standard case */

  if (!exclude_diag && mc->d_val != NULL) {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      const cs_lnum_t *restrict col_id = ms->col_id + ms->row_index[ii];
      const cs_real_t *restrict m_row = mc->x_val + ms->row_index[ii];
      cs_lnum_t n_cols = ms->row_index[ii+1] - ms->row_index[ii];

      _dense_6_6_ax(ii, mc->d_val, x, y);

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        for (cs_lnum_t kk = 0; kk < 6; kk++)
          y[ii*6 + kk] += (m_row[jj]*x[col_id[jj]*6 + kk]);
      }

    }

  }

  /* Exclude diagonal */

  else {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      const cs_lnum_t *restrict col_id = ms->col_id + ms->row_index[ii];
      const cs_real_t *restrict m_row = mc->x_val + ms->row_index[ii];
      cs_lnum_t n_cols = ms->row_index[ii+1] - ms->row_index[ii];

      for (cs_lnum_t kk = 0; kk < 6; kk++)
        y[ii*6 + kk] = 0.;

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        for (cs_lnum_t kk = 0; kk < 6; kk++)
          y[ii*6 + kk] += (m_row[jj]*x[col_id[jj]*6 + kk]);
      }

    }
  }

}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with MSR matrix, blocked version.
 *
 * This variant uses fixed block size variants for common cases.
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- pointer to matrix structure
 *   x            <-- multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_b_mat_vec_p_l_msr(bool                exclude_diag,
                   const cs_matrix_t  *matrix,
                   const cs_real_t     x[restrict],
                   cs_real_t           y[restrict])
{
  if (matrix->db_size[0] == 3 && matrix->db_size[3] == 9)
    _3_3_mat_vec_p_l_msr(exclude_diag, matrix, x, y);

  else if (matrix->db_size[0] == 6 && matrix->db_size[3] == 36)
    _6_6_mat_vec_p_l_msr(exclude_diag, matrix, x, y);

  else
    _b_mat_vec_p_l_msr_generic(exclude_diag, matrix, x, y);
}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with MSR matrix, using MKL
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- pointer to matrix structure
 *   x            <-- multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

#if defined (HAVE_MKL)

static void
_mat_vec_p_l_msr_mkl(bool                exclude_diag,
                     const cs_matrix_t  *matrix,
                     const cs_real_t     x[restrict],
                     cs_real_t           y[restrict])
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
 * Synchronize ghost values prior to matrix.vector product
 *
 * parameters:
 *   rotation_mode <-- halo update option for rotational periodicity
 *   matrix        <-- pointer to matrix structure
 *   x             <-> multipliying vector values (ghost values updated)
 *----------------------------------------------------------------------------*/

static void
_pre_vector_multiply_sync_x(cs_halo_rotation_t   rotation_mode,
                            const cs_matrix_t   *matrix,
                            cs_real_t            x[restrict])
{
  assert(matrix->halo != NULL);

  /* Non-blocked version */

  if (matrix->db_size[3] == 1) {

    if (matrix->halo != NULL)
      cs_halo_sync_component(matrix->halo,
                             CS_HALO_STANDARD,
                             rotation_mode,
                             x);

  }

  /* Blocked version */

  else { /* if (matrix->db_size[3] > 1) */

    const cs_lnum_t *db_size = matrix->db_size;

    /* Update distant ghost rows */

    if (matrix->halo != NULL) {

      cs_halo_sync_var_strided(matrix->halo,
                               CS_HALO_STANDARD,
                               x,
                               db_size[1]);

      /* Synchronize periodic values */

#if !defined(_CS_UNIT_MATRIX_TEST) /* unit tests do not link with full library */

      if (matrix->halo->n_transforms > 0) {
        if (db_size[0] == 3)
          cs_halo_perio_sync_var_vect(matrix->halo,
                                      CS_HALO_STANDARD,
                                      x,
                                      db_size[1]);
        else if (db_size[0] == 6)
          cs_halo_perio_sync_var_sym_tens(matrix->halo,
                                          CS_HALO_STANDARD,
                                          x);
      }

#endif

    }

  }
}

/*----------------------------------------------------------------------------
 * Zero ghost values prior to matrix.vector product
 *
 * parameters:
 *   matrix        <-- pointer to matrix structure
 *   y             --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_pre_vector_multiply_sync_y(const cs_matrix_t   *matrix,
                            cs_real_t            y[restrict])
{
  size_t n_cols_ext = matrix->n_cols_ext;

  if (matrix->db_size[3] == 1)
    _zero_range(y, matrix->n_rows, n_cols_ext);

  else
    _b_zero_range(y, matrix->n_rows, n_cols_ext, matrix->db_size);
}

/*----------------------------------------------------------------------------
 * Synchronize ghost values prior to matrix.vector product
 *
 * parameters:
 *   rotation_mode <-- halo update option for rotational periodicity
 *   matrix        <-- pointer to matrix structure
 *   x             <-> multipliying vector values (ghost values updated)
 *   y             --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_pre_vector_multiply_sync(cs_halo_rotation_t   rotation_mode,
                          const cs_matrix_t   *matrix,
                          cs_real_t           *restrict x,
                          cs_real_t           *restrict y)
{
  _pre_vector_multiply_sync_y(matrix, y);

  _pre_vector_multiply_sync_x(rotation_mode, matrix, x);
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
 *   n_rows      <-- local number of rows
 *   n_cols_ext  <-- number of local + ghost columns
 *   n_edges     <-- local number of (undirected) graph edges
 *   edges       <-- edges (symmetric row <-> column) connectivity
 *   halo        <-- cell halo structure
 *   numbering   <-- vectorization or thread-related numbering info, or NULL
 *   m_variant   <-> array of matrix variants
 *----------------------------------------------------------------------------*/

static void
_matrix_check(int                    n_variants,
              cs_lnum_t              n_rows,
              cs_lnum_t              n_cols_ext,
              cs_lnum_t              n_edges,
              const cs_lnum_2_t     *edges,
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
  cs_lnum_t d_block_size[4] = {3, 3, 3, 9};
  cs_lnum_t ed_block_size[4] = {3, 3, 3, 9};

  /* Allocate and initialize  working arrays */

  if (CS_MEM_ALIGN > 0) {
    BFT_MEMALIGN(x, CS_MEM_ALIGN, n_cols_ext*d_block_size[1], cs_real_t);
    BFT_MEMALIGN(y, CS_MEM_ALIGN, n_cols_ext*d_block_size[1], cs_real_t);
    BFT_MEMALIGN(yr0, CS_MEM_ALIGN, n_cols_ext*d_block_size[1], cs_real_t);
    BFT_MEMALIGN(yr1, CS_MEM_ALIGN, n_cols_ext*d_block_size[1], cs_real_t);
  }
  else {
    BFT_MALLOC(x, n_cols_ext*d_block_size[1], cs_real_t);
    BFT_MALLOC(y, n_cols_ext*d_block_size[1], cs_real_t);
    BFT_MALLOC(yr0, n_cols_ext*d_block_size[1], cs_real_t);
    BFT_MALLOC(yr1, n_cols_ext*d_block_size[1], cs_real_t);
  }

  BFT_MALLOC(da, n_cols_ext*d_block_size[3], cs_real_t);
  BFT_MALLOC(xa, n_edges*2*ed_block_size[3], cs_real_t);

  /* Initialize arrays */

# pragma omp parallel for
  for (ii = 0; ii < n_cols_ext*d_block_size[3]; ii++)
    da[ii] = 1.0 + cos(ii);

# pragma omp parallel for
  for (ii = 0; ii < n_edges*ed_block_size[3]; ii++) {
    xa[ii*2] = 0.5*(0.9 + cos(ii));
    xa[ii*2 + 1] = -0.5*(0.9 + cos(ii));
  }

# pragma omp parallel for
  for (ii = 0; ii < n_cols_ext*d_block_size[1]; ii++)
    x[ii] = sin(ii);

  /* Loop on fill options */

  for (f_id = 0; f_id < CS_MATRIX_N_FILL_TYPES; f_id++) {

    const cs_lnum_t *_d_block_size
      = (f_id >= CS_MATRIX_BLOCK_D) ? d_block_size : NULL;
    const cs_lnum_t *_ed_block_size
      = (f_id >= CS_MATRIX_BLOCK) ? ed_block_size : NULL;
    const cs_lnum_t _block_mult = (_d_block_size != NULL) ? d_block_size[1] : 1;
    const bool sym_coeffs = (   f_id == CS_MATRIX_SCALAR_SYM
                             || f_id == CS_MATRIX_BLOCK_D_SYM) ? true : false;

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
                                        n_rows,
                                        n_cols_ext,
                                        n_edges,
                                        edges,
                                        halo,
                                        numbering);
        m = cs_matrix_create(ms);

        cs_matrix_set_coefficients(m,
                                   sym_coeffs,
                                   _d_block_size,
                                   _ed_block_size,
                                   n_edges,
                                   edges,
                                   da,
                                   xa);

        /* Check multiplication */

        vector_multiply(ed_flag, m, x, y);
        if (v_id == 0)
          memcpy(yr0, y, n_rows*_block_mult*sizeof(cs_real_t));
        else {
          double dmax = _matrix_check_compare(n_rows*_block_mult, y, yr0);
          if (print_subtitle) {
            bft_printf("\n%s\n",
                       _matrix_operation_name[f_id][ed_flag]);
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
      v->matrix_vector_cost[i][j][0] = -1.;
      v->matrix_vector_cost[i][j][1] = -1.;
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

    switch(mft) {

    case CS_MATRIX_SCALAR:
    case  CS_MATRIX_SCALAR_SYM:
      if (ed_flag != 1)
        v->vector_multiply[mft][0] = vector_multiply;
      if (ed_flag != 0)
        v->vector_multiply[mft][1] = vector_multiply;
      break;

    case CS_MATRIX_BLOCK_D:
    case CS_MATRIX_BLOCK_D_66:
    case CS_MATRIX_BLOCK_D_SYM:
      if (ed_flag != 1)
        v->vector_multiply[mft][0] = b_vector_multiply;
      if (ed_flag != 0)
        v->vector_multiply[mft][1] = b_vector_multiply;
      break;

    case CS_MATRIX_BLOCK:
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
 *     omp             (for OpenMP with compatible numbering)
 *     vector          (For vector machine with compatible numbering)
 *
 *   CS_MATRIX_CSR     (for CS_MATRIX_SCALAR or CS_MATRIX_SCALAR_SYM)
 *     default
 *     standard
 *     mkl             (with MKL)
 *
 *   CS_MATRIX_CSR_SYM (for CS_MATRIX_SCALAR_SYM)
 *     default
 *     standard
 *     mkl             (with MKL)
 *
 *   CS_MATRIX_MSR     (all fill types except CS_MATRIX_33_BLOCK)
 *     standard
 *     omp_sched       (Improved scheduling for OpenMP)
 *     mkl             (with MKL, for CS_MATRIX_SCALAR or CS_MATRIX_SCALAR_SYM)
 *
 * parameters:
 *   m_type          <-- Matrix type
 *   numbering       <-- mesh numbering type, or NULL
 *   fill type       <-- matrix fill type to merge from
 *   ed_flag         <-- 0: with diagonal only, 1 exclude only; 2; both
 *   func_name       <-- function type name, or NULL for default
 *   vector_multiply <-> multiplication function array
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
               cs_matrix_vector_product_t  *vector_multiply[][2])
{
  int retcode = 1;
  int standard = 0;

  cs_matrix_vector_product_t *spmv[2] = {NULL, NULL};

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
      case CS_MATRIX_BLOCK_D:
      case CS_MATRIX_BLOCK_D_66:
      case CS_MATRIX_BLOCK_D_SYM:
        spmv[0] = _b_mat_vec_p_l_native;
        spmv[1] = _b_mat_vec_p_l_native;
        break;
      case CS_MATRIX_BLOCK:
        spmv[0] = _bb_mat_vec_p_l_native;
        spmv[1] = _bb_mat_vec_p_l_native;
        break;
        break;
      default:
        break;
      }

      if (standard > 1) { /* default optimized variants */
        switch(fill_type) {
        case CS_MATRIX_SCALAR:
        case CS_MATRIX_SCALAR_SYM:
          if (numbering != NULL) {
#if defined(HAVE_OPENMP)
            if (numbering->type == CS_NUMBERING_THREADS) {
              spmv[0] = _mat_vec_p_l_native_omp;
              spmv[1] = _mat_vec_p_l_native_omp;
            }
#endif
            if (numbering->type == CS_NUMBERING_VECTORIZE) {
              spmv[0] = _mat_vec_p_l_native_vector;
              spmv[1] = _mat_vec_p_l_native_vector;
            }
          }
          break;
        case CS_MATRIX_BLOCK_D:
        case CS_MATRIX_BLOCK_D_66:
        case CS_MATRIX_BLOCK_D_SYM:
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

    else if (!strcmp(func_name, "fixed")) {
      switch(fill_type) {
      case CS_MATRIX_BLOCK_D:
      case CS_MATRIX_BLOCK_D_66:
      case CS_MATRIX_BLOCK_D_SYM:
        spmv[0] = _b_mat_vec_p_l_native_fixed;
        spmv[1] = _b_mat_vec_p_l_native_fixed;
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
          case CS_MATRIX_BLOCK_D:
          case CS_MATRIX_BLOCK_D_66:
          case CS_MATRIX_BLOCK_D_SYM:
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
      switch(fill_type) {
      case CS_MATRIX_SCALAR:
      case CS_MATRIX_SCALAR_SYM:
        spmv[0] = _mat_vec_p_l_native_vector;
        spmv[1] = _mat_vec_p_l_native_vector;
        break;
      default:
        break;
      }
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
      case CS_MATRIX_BLOCK_D:
      case CS_MATRIX_BLOCK_D_66:
      case CS_MATRIX_BLOCK_D_SYM:
        spmv[0] = _b_mat_vec_p_l_msr;
        spmv[1] = _b_mat_vec_p_l_msr;
        break;
      default:
        break;
      }
    }

    else if (!strcmp(func_name, "generic")) {
      switch(fill_type) {
      case CS_MATRIX_BLOCK_D:
      case CS_MATRIX_BLOCK_D_66:
      case CS_MATRIX_BLOCK_D_SYM:
        spmv[0] = _b_mat_vec_p_l_msr_generic;
        spmv[1] = _b_mat_vec_p_l_msr_generic;
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

    else if (!strcmp(func_name, "omp_sched")) {
      switch(fill_type) {
      case CS_MATRIX_SCALAR:
      case CS_MATRIX_SCALAR_SYM:
        spmv[0] = _mat_vec_p_l_msr_omp_sched;
        spmv[1] = _mat_vec_p_l_msr_omp_sched;
        break;
      default:
        break;
      }
    }

    break;

  default:
    break;
  }

  if (ed_flag != 1 && spmv[0] != NULL) {
    vector_multiply[fill_type][0] = spmv[0];
    retcode = 0;
  }
  if (ed_flag != 0 && spmv[0] != NULL) {
    vector_multiply[fill_type][1] = spmv[1];
    retcode = 0;
  }

  return retcode;
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
      BFT_MALLOC(_row_index, n_rows + 1, cs_lnum_t);
      BFT_MALLOC(_col_id, row_index[n_rows] + n_rows, cs_lnum_t);
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
      structure = _create_struct_csr_from_shared(false,
                                                 false, /* for safety */
                                                 n_rows,
                                                 n_cols_ext,
                                                 row_index,
                                                 col_id);
    else {
      cs_lnum_t *_row_index, *_col_id;
      BFT_MALLOC(_row_index, n_rows + 1, cs_lnum_t);
      BFT_MALLOC(_col_id, row_index[n_rows], cs_lnum_t);
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
      BFT_REALLOC(_col_id, _row_index[n_rows], cs_lnum_t);
      structure = _create_struct_csr_from_csr(false,
                                              true,
                                              true,
                                              n_rows,
                                              n_cols_ext,
                                              &_row_index,
                                              &_col_id);
    }
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _("%s: handling of matrices in %s format\n"
                "is not operational yet."),
              __func__,
              _(cs_matrix_type_name[type]));
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
    {
      cs_matrix_struct_native_t *_structure = *structure;
      _destroy_struct_native(&_structure);
      *structure = _structure;
    }
    break;
  case CS_MATRIX_CSR:
    {
      cs_matrix_struct_csr_t *_structure = *structure;
      _destroy_struct_csr(&_structure);
      *structure = _structure;
    }
    break;
  case CS_MATRIX_CSR_SYM:
    {
      cs_matrix_struct_csr_sym_t *_structure = *structure;
      _destroy_struct_csr_sym(&_structure);
      *structure = _structure;
    }
    break;
  case CS_MATRIX_MSR:
    {
      cs_matrix_struct_csr_t *_structure = *structure;
      _destroy_struct_csr(&_structure);
      *structure = _structure;
    }
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
  int i;
  cs_matrix_fill_type_t mft;
  cs_matrix_t *m;

  BFT_MALLOC(m, 1, cs_matrix_t);

  m->type = type;

  /* Map shared structure */

  m->n_rows = 0;
  m->n_cols_ext = 0;

  if (m->type != CS_MATRIX_CSR_SYM)
    m->symmetric = false;
  else
    m->symmetric = true;

  for (i = 0; i < 4; i++) {
    m->db_size[i] = 0;
    m->eb_size[i] = 0;
  }
  m->fill_type = CS_MATRIX_N_FILL_TYPES;

  m->structure = NULL;
  m->_structure = NULL;

  m->halo = NULL;
  m->numbering = NULL;
  m->assembler = NULL;

  for (mft = 0; mft < CS_MATRIX_N_FILL_TYPES; mft++) {
    for (i = 0; i < 2; i++)
      m->vector_multiply[mft][i] = NULL;
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

  m->xa = NULL;

  /* Set function pointers here */

  m->set_coefficients = NULL;

  for (mft = 0; mft < CS_MATRIX_N_FILL_TYPES; mft++)
    _set_spmv_func(m->type,
                   m->numbering,
                   mft,
                   2,    /* ed_flag */
                   NULL, /* func_name */
                   m->vector_multiply);

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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

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
 * Note that the resulting matrix structure will contain either a full or
 * an empty main diagonal, and that the extra-diagonal structure is always
 * symmetric (though the coefficients my not be, and we may choose a
 * matrix format that does not exploit this symmetry). If the edges
 * connectivity argument is NULL, the matrix will be purely diagonal.
 *
 * \param[in]  type        type of matrix considered
 * \param[in]  have_diag   indicates if the diagonal structure
 *                         contains nonzeroes
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
                           bool                   have_diag,
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

  /* Define Structure */

  switch(ms->type) {
  case CS_MATRIX_NATIVE:
    ms->structure = _create_struct_native(n_rows,
                                          n_cols_ext,
                                          n_edges,
                                          edges);
    break;
  case CS_MATRIX_CSR:
    ms->structure = _create_struct_csr(have_diag,
                                       n_rows,
                                       n_cols_ext,
                                       n_edges,
                                       edges);
    break;
  case CS_MATRIX_CSR_SYM:
    ms->structure = _create_struct_csr_sym(have_diag,
                                           n_rows,
                                           n_cols_ext,
                                           n_edges,
                                           edges);
    break;
  case CS_MATRIX_MSR:
    ms->structure = _create_struct_csr(false,
                                       n_rows,
                                       n_cols_ext,
                                       n_edges,
                                       edges);
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Handling of matrixes in %s format\n"
                "is not operational yet."),
              _(cs_matrix_type_name[type]));
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
 * \param[in, out]  col_id      pointer to array of colum ids related to
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
    ms->structure = _create_struct_csr_from_csr(false,
                                                transfer,
                                                false,
                                                n_rows,
                                                n_cols_ext,
                                                row_index,
                                                col_id);
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _("%s: handling of matrices in %s format\n"
                "is not operational yet."),
              __func__,
              _(cs_matrix_type_name[type]));
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
 * definition as well as an optional edge-based definition.
 *
 * Note that as the structure created maps to the given existing
 * cell global number, face -> cell connectivity arrays, and cell halo
 * structure, it must be destroyed before they are freed
 * (usually along with the code's main face -> cell structure).
 *
 * \param[in]  have_diag        indicates if the structure includes the
 *                              diagonal (should be the same for all rows)
 * \param[in]  direct_assembly  true if each value corresponds to
                                a unique face
 * \param[in]  n_rows           local number of rows
 * \param[in]  n_cols_ext       local number of columns + ghosts
 * \param[in]  row_index        index on rows
 * \param[in]  col_id           array of colum ids related to the row index
 * \param[in]  halo             halo structure for synchronization, or NULL
 * \param[in]  numbering        vectorization or thread-related numbering
 *                              info, or NULL
 *
 * \returns  a pointer to a created CDO matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_matrix_structure_t *
cs_matrix_structure_create_msr_shared(bool                    have_diag,
                                      bool                    direct_assembly,
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

  /* Define Structure */

  ms->structure = _create_struct_csr_from_shared(have_diag,
                                                 direct_assembly,
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
 * \brief Create a matrix container using a given variant.
 *
 * If the matrix variant is incompatible with the structure, it is ignored,
 * and defaults for that structure are used instead.
 *
 * \param[in]  ms  associated matrix structure
 * \param[in]  mv  associated matrix variant
 *
 * \return  pointer to created matrix structure;
 */
/*----------------------------------------------------------------------------*/

cs_matrix_t *
cs_matrix_create_by_variant(const cs_matrix_structure_t  *ms,
                            const cs_matrix_variant_t    *mv)
{
  cs_matrix_t *m = cs_matrix_create(ms);

  m->assembler = ms->assembler;

  if (mv != NULL) {
    if (mv->type == ms->type) {
      for (int i = 0; i < CS_MATRIX_N_FILL_TYPES; i++) {
        for (int j = 0; j < 2; j++) {
          if (mv->vector_multiply[i][j] != NULL)
            m->vector_multiply[i][j] = mv->vector_multiply[i][j];
        }
      }
    }
  }

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

  cs_matrix_release_coefficients(m);

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

    if (m->_structure != NULL)
      _structure_destroy(m->type, &(m->_structure));

    /* Now free main structure */

    BFT_FREE(*matrix);
  }
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
    bft_error(__FILE__, __LINE__, 0,
              _("The matrix is not defined."));
  return matrix->type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of columns in a matrix.
 *
 * \param[in]  matrix  pointer to matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_matrix_get_n_columns(const cs_matrix_t  *matrix)
{
  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("The matrix is not defined."));
  return matrix->n_cols_ext;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of rows in matrix.
 *
 * \param[in]  matrix  pointer to matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_matrix_get_n_rows(const cs_matrix_t  *matrix)
{
  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("The matrix is not defined."));
  return matrix->n_rows;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of entries in matrix.
 *
 * When the block size is > 1, the number reported is the number of
 * entry blocks, not individual entries.
 *
 * \param[in]  matrix  pointer to matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_matrix_get_n_entries(const cs_matrix_t  *matrix)
{
  cs_lnum_t retval = 0;

  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("The matrix is not defined."));

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
  case CS_MATRIX_CSR_SYM:
    {
      const cs_matrix_struct_csr_t  *ms = matrix->structure;
      retval = ms->row_index[ms->n_rows]*2 - ms->n_rows;
    }
    break;
  case CS_MATRIX_MSR:
    {
      const cs_matrix_struct_csr_t  *ms = matrix->structure;
      retval = ms->row_index[ms->n_rows] + ms->n_rows;
    }
    break;
  default:
    break;
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return matrix diagonal block sizes.
 *
 * Block sizes are defined by a array of 4 values:
 *   0: useful block size, 1: vector block extents,
 *   2: matrix line extents,  3: matrix line*column extents
 *
 * \param[in]  matrix  pointer to matrix structure
 *
 * \return  pointer to block sizes
 */
/*----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_matrix_get_diag_block_size(const cs_matrix_t  *matrix)
{
  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("The matrix is not defined."));

  return matrix->db_size;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return matrix extra-diagonal block sizes.
 *
 * Block sizes are defined by a array of 4 values:
 *   0: useful block size, 1: vector block extents,
 *   2: matrix line extents,  3: matrix line*column extents
 *
 * \param[in]  matrix  pointer to matrix structure
 *
 * \return  pointer to block sizes
 */
/*----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_matrix_get_extra_diag_block_size(const cs_matrix_t  *matrix)
{
  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("The matrix is not defined."));

  return matrix->eb_size;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to matrix halo structure.
 *
 * \param[in]  matrix  pointer to matrix structure
 *
 * \return  pointer to halo strucuture
 */
/*----------------------------------------------------------------------------*/

const cs_halo_t *
cs_matrix_get_halo(const cs_matrix_t  *matrix)
{
  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("The matrix is not defined."));

  return matrix->halo;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get matrix fill type, depending on block sizes.
 *
 * Block sizes are defined by an optional array of 4 values:
 *   0: useful block size, 1: vector block extents,
 *   2: matrix line extents,  3: matrix line*column extents
 *
 * \param[in]  symmetric              indicates if matrix coefficients
 *                                    are symmetric
 * \param[in]  diag_block_size        block sizes for diagonal, or NULL
 * \param[in]  extra_diag_block_size  block sizes for extra diagonal, or NULL
 *
 * \return  matrix fill type
 */
/*----------------------------------------------------------------------------*/

cs_matrix_fill_type_t
cs_matrix_get_fill_type(bool              symmetric,
                        const cs_lnum_t  *diag_block_size,
                        const cs_lnum_t  *extra_diag_block_size)
{
  cs_matrix_fill_type_t fill_type = CS_MATRIX_N_FILL_TYPES;

  cs_lnum_t _db_size = 1, _eb_size = 1;
  if (diag_block_size != NULL)
    _db_size = diag_block_size[0];

  if (extra_diag_block_size != NULL)
    _eb_size = extra_diag_block_size[0];

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
 * Block sizes are defined by an optional array of 4 values:
 *   0: useful block size, 1: vector block extents,
 *   2: matrix line extents,  3: matrix line*column extents
 *
 * \param[in, out]  matrix                 pointer to matrix structure
 * \param[in]       symmetric              indicates if matrix coefficients
 *                                         are symmetric
 * \param[in]       diag_block_size        block sizes for diagonal, or NULL
 * \param[in]       extra_diag_block_size  block sizes for extra diagonal,
 *                                         or NULL
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
                           const cs_lnum_t    *diag_block_size,
                           const cs_lnum_t    *extra_diag_block_size,
                           const cs_lnum_t     n_edges,
                           const cs_lnum_2_t   edges[],
                           const cs_real_t    *da,
                           const cs_real_t    *xa)
{
  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("The matrix is not defined."));

  cs_base_check_bool(&symmetric);

  /* Set fill type */
  _set_fill_info(matrix,
                 symmetric,
                 diag_block_size,
                 extra_diag_block_size);

  /* Set coefficients */

  if (matrix->set_coefficients != NULL) {
    matrix->xa = xa;
    matrix->set_coefficients(matrix, symmetric, false, n_edges, edges, da, xa);
  }
  else
    bft_error
      (__FILE__, __LINE__, 0,
       "Matrix format %s with fill type %s does not handle\n"
       "coefficient assignment from native (graph-edge) coefficients.",
       cs_matrix_type_name[matrix->type],
       cs_matrix_fill_type_name[matrix->fill_type]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set matrix coefficients, copying values to private arrays.
 *
 * With private arrays, the matrix becomes independant from the
 * arrays passed as arguments.
 *
 * Block sizes are defined by an optional array of 4 values:
 *   0: useful block size, 1: vector block extents,
 *   2: matrix line extents,  3: matrix line*column extents
 *
 * \param[in, out]  matrix                 pointer to matrix structure
 * \param[in]       symmetric              indicates if matrix coefficients
 *                                         are symmetric
 * \param[in]       diag_block_size        block sizes for diagonal, or NULL
 * \param[in]       extra_diag_block_size  block sizes for extra diagonal,
 *                                         or NULL
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
                            const cs_lnum_t    *diag_block_size,
                            const cs_lnum_t    *extra_diag_block_size,
                            const cs_lnum_t     n_edges,
                            const cs_lnum_2_t   edges[],
                            const cs_real_t    *da,
                            const cs_real_t    *xa)
{
  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("The matrix is not defined."));

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
       cs_matrix_type_name[matrix->type],
       cs_matrix_fill_type_name[matrix->fill_type]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set matrix coefficients in an MSR format, transfering the
 * property of those arrays to the matrix.
 *
 * If the matrix is also in MSR format, this avoids an extra copy.
 * If it is in a different format, values are copied to the structure,
 * and the original arrays freed. In any case, the arrays pointers passed as
 * arguments are set to NULL, to help ensure the caller does not use the
 * original arrays directly after this call.
 *
 * Block sizes are defined by an optional array of 4 values:
 *   0: useful block size, 1: vector block extents,
 *   2: matrix line extents,  3: matrix line*column extents
 *
 * \param[in, out]  matrix                 pointer to matrix structure
 * \param[in]       symmetric              indicates if matrix coefficients
 *                                         are symmetric
 * \param[in]       diag_block_size        block sizes for diagonal, or NULL
 * \param[in]       extra_diag_block_size  block sizes for extra diagonal,
 *                                         or NULL
 * \param[in]       row_index              MSR row index (0 to n-1)
 * \param[in]       col_id                 MSR column id (0 to n-1)
 * \param[in, out]  d_val                  diagonal values (NULL if zero)
 * \param[in, out]  x_val                  extradiagonal values (NULL if zero)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_transfer_coefficients_msr(cs_matrix_t         *matrix,
                                    bool                 symmetric,
                                    const cs_lnum_t     *diag_block_size,
                                    const cs_lnum_t     *extra_diag_block_size,
                                    const cs_lnum_t      row_index[],
                                    const cs_lnum_t      col_id[],
                                    cs_real_t          **d_val,
                                    cs_real_t          **x_val)
{
  const cs_real_t  *d_val_p = (d_val != NULL) ? *d_val : NULL;
  const cs_real_t  *x_val_p = (x_val != NULL) ? *x_val : NULL;

  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("The matrix is not defined."));

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
       cs_matrix_type_name[matrix->type],
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
    bft_error(__FILE__, __LINE__, 0,
              _("The matrix is not defined."));

  if (matrix->release_coefficients != NULL) {
    matrix->xa = NULL;
    matrix->release_coefficients(matrix);
  }
  else
    bft_error
      (__FILE__, __LINE__, 0,
       "Matrix format %s is missing a release_coefficients function.",
       cs_matrix_type_name[matrix->type]);

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
 * Block sizes are defined by an optional array of 4 values:
 *   0: useful block size, 1: vector block extents,
 *   2: matrix line extents,  3: matrix line*column extents
 *
 * \param[in, out]  matrix                 pointer to matrix structure
 * \param[in]       diag_block_size        block sizes for diagonal, or NULL
 * \param[in]       extra_diag_block_size  block sizes for extra diagonal,
 *                                         or NULL
 *
 * \return  pointer to initialized matrix assembler values structure;
 */
/*----------------------------------------------------------------------------*/

cs_matrix_assembler_values_t *
cs_matrix_assembler_values_init(cs_matrix_t      *matrix,
                                const cs_lnum_t  *diag_block_size,
                                const cs_lnum_t  *extra_diag_block_size)
{
  cs_matrix_assembler_values_t *mav = NULL;

  /* Set fill type */

  _set_fill_info(matrix,
                 false, /* symmetric */
                 diag_block_size,
                 extra_diag_block_size);

  /* Create values assembler */

  switch(matrix->type) {
  case CS_MATRIX_CSR:
    mav = cs_matrix_assembler_values_create(matrix->assembler,
                                            false,
                                            diag_block_size,
                                            extra_diag_block_size,
                                            (void *)matrix,
                                            _csr_assembler_values_init,
                                            _csr_assembler_values_add,
                                            NULL,
                                            NULL,
                                            NULL);
    break;
  case CS_MATRIX_MSR:
    mav = cs_matrix_assembler_values_create(matrix->assembler,
                                            true,
                                            diag_block_size,
                                            extra_diag_block_size,
                                            (void *)matrix,
                                            _msr_assembler_values_init,
                                            _msr_assembler_values_add,
                                            NULL,
                                            NULL,
                                            NULL);
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _("%s: handling of matrices in %s format\n"
                "is not operational yet."),
              __func__,
              _(cs_matrix_type_name[matrix->type]));
    break;
  }

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
    bft_error(__FILE__, __LINE__, 0,
              _("The matrix is not defined."));

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
  cs_lnum_t ii;

  const cs_real_t  *diag = NULL;

  switch(matrix->type) {

  case CS_MATRIX_NATIVE:
    {
      cs_matrix_coeff_native_t *mc = matrix->coeffs;
      if (mc->da == NULL) {
        cs_lnum_t n_rows = matrix->n_rows * matrix->db_size[3];
        if (mc->_da == NULL || mc->max_db_size < matrix->db_size[3]) {
          BFT_REALLOC(mc->_da, matrix->db_size[3]*matrix->n_rows, cs_real_t);
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
        BFT_MALLOC(mc->_d_val, matrix->n_rows, cs_real_t);
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
        BFT_MALLOC(mc->_d_val, matrix->n_rows, cs_real_t);
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
        cs_lnum_t n_rows = matrix->n_rows * matrix->db_size[3];
        if (mc->_d_val == NULL || mc->max_db_size < matrix->db_size[3]) {
          BFT_REALLOC(mc->_d_val, matrix->db_size[3]*matrix->n_rows, cs_real_t);
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
 * It is used in the current multgrid code, but should be removed as soon
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
 * In the case of blocked matrixes, the true (non-blocked)
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
  cs_lnum_t b_size = matrix->db_size[0];

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
      const cs_matrix_coeff_csr_t  *mc = matrix->coeffs;
      const cs_lnum_t n_ed_cols =   ms->row_index[_row_id+1]
                                  - ms->row_index[_row_id];
      if (b_size == 1)
        r->row_size = n_ed_cols + 1;
      else if (matrix->eb_size[0] == 1)
        r->row_size = n_ed_cols*b_size;
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
        const cs_real_t *m_row = mc->val + ms->row_index[_row_id];
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
      else if (matrix->eb_size[0] == 1) {
        const cs_lnum_t _sub_id = row_id % b_size;
        const cs_lnum_t *db_size = matrix->db_size;
        const cs_real_t *m_row = mc->val + ms->row_index[_row_id];
        for (jj = 0; jj < n_ed_cols && c_id[jj] < _row_id; jj++) {
          r->_col_id[ii] = c_id[jj]*b_size + _sub_id;
          r->_vals[ii++] = m_row[jj];
        }
        for (cs_lnum_t kk = 0; kk < b_size; kk++) {
          r->_col_id[ii] = _row_id*b_size + kk;
          r->_vals[ii++] = mc->d_val[  _row_id*db_size[3]
                                     + _sub_id*db_size[2] + kk];
        }
        for (; jj < n_ed_cols; jj++) {
          r->_col_id[ii] = c_id[jj]*b_size + _sub_id;
          r->_vals[ii++] = m_row[jj];
        }
      }
      else {
        const cs_lnum_t _sub_id = row_id % b_size;
        const cs_lnum_t *db_size = matrix->db_size;
        const cs_lnum_t *eb_size = matrix->db_size;
        const cs_real_t *m_row = mc->val + ms->row_index[_row_id]*eb_size[3];
        for (jj = 0; jj < n_ed_cols && c_id[jj] < _row_id; jj++) {
          for (cs_lnum_t kk = 0; kk < b_size; kk++) {
            r->_col_id[ii] = c_id[jj]*b_size + kk;
            r->_vals[ii++] = m_row[_sub_id*eb_size[2] + kk];
          }
        }
        for (cs_lnum_t kk = 0; kk < b_size; kk++) {
          r->_col_id[ii] = _row_id*b_size + kk;
          r->_vals[ii++] = mc->d_val[  _row_id*db_size[3]
                                     + _sub_id*db_size[2] + kk];
        }
        for (; jj < n_ed_cols; jj++) {
          for (cs_lnum_t kk = 0; kk < b_size; kk++) {
            r->_col_id[ii] = c_id[jj]*b_size + kk;
            r->_vals[ii++] = m_row[_sub_id*eb_size[2] + kk];
          }
        }
      }
    }
    break;

  default:
    bft_error
      (__FILE__, __LINE__, 0,
       _("Matrix format %s with fill type %s does not handle %s operation."),
       cs_matrix_type_name[matrix->type],
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
    const cs_matrix_coeff_native_t  *mc = matrix->coeffs;
    if (n_edges != NULL)
      *n_edges = ms->n_edges;
    if (edges != NULL)
      *edges = ms->edges;
    if (mc != NULL) {
      if (symmetric != NULL)
        *symmetric = mc->symmetric;
      if (d_val != NULL)
        *d_val = mc->da;
      if (x_val != NULL)
        *x_val = mc->xa;
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
    const cs_matrix_struct_csr_t  *ms = matrix->structure;
    const cs_matrix_coeff_msr_t  *mc = matrix->coeffs;
    if (row_index != NULL)
      *row_index = ms->row_index;
    if (col_id != NULL)
      *col_id = ms->col_id;
    if (mc != NULL) {
      if (d_val != NULL)
        *d_val = mc->d_val;
      if (x_val != NULL)
        *x_val = mc->x_val;
    }
  }
  else
    bft_error
      (__FILE__, __LINE__, 0,
       _("%s is not available for matrix using %s storage."),
       __func__,
       cs_matrix_type_name[matrix->type]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = A.x
 *
 * This function includes a halo update of x prior to multiplication by A.
 *
 * \param[in]       rotation_mode  halo update option for
 *                                 rotational periodicity
 * \param[in]       matrix         pointer to matrix structure
 * \param[in, out]  x              multipliying vector values
 *                                 (ghost values updated)
 * \param[out]      y              resulting vector
 */
/*----------------------------------------------------------------------------*/

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
 * \param[in]   x              multipliying vector values
 * \param[out]  y              resulting vector
 */
/*----------------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix.vector product y = (A-D).x
 *
 * This function includes a halo update of x prior to multiplication by A.
 *
 * \param[in]       rotation_mode  halo update option for
 *                                 rotational periodicity
 * \param[in]       matrix         pointer to matrix structure
 * \param[in, out]  x              multipliying vector values
 *                                 (ghost values updated)
 * \param[out]      y              resulting vector
 */
/*----------------------------------------------------------------------------*/

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
 * Synchronize ghost values prior to matrix.vector product
 *
 * parameters:
 *   rotation_mode <-- halo update option for rotational periodicity
 *   matrix        <-- pointer to matrix structure
 *   x             <-> multipliying vector values (ghost values updated)
 *----------------------------------------------------------------------------*/

void
cs_matrix_pre_vector_multiply_sync(cs_halo_rotation_t   rotation_mode,
                                   const cs_matrix_t   *matrix,
                                   cs_real_t           *x)
{
  if (matrix->halo != NULL)
    _pre_vector_multiply_sync_x(rotation_mode, matrix, x);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build matrix variant
 *
 * The variant will initially use default matrix-vector functions,
 * which can be later modified using cs_matrix_variant_set_func().
 *
 * \param[in]  type       type of matrix considered
 * \param[in]  numbering  vectorization or thread-related numbering info,
 *                        or NULL
 */
/*----------------------------------------------------------------------------*/

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
                          mv->vector_multiply);
  }

  return mv;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build list of variants for tuning or testing.
 *
 * \param[in]   n_fill_types  number of fill types tuned for
 * \param[in]   fill_types    array of fill types tuned for
 * \param[in]   type_filter   true for matrix types tuned for, false for others
 * \param[in]   numbering     vectorization or thread-related numbering info,
 *                            or NULL
 * \param[out]  n_variants    number of variants
 * \param[out]  m_variant     array of matrix variants
 */
/*----------------------------------------------------------------------------*/

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
                 _mat_vec_p_l_native,
                 _b_mat_vec_p_l_native,
                 _bb_mat_vec_p_l_native,
                 n_variants,
                 &n_variants_max,
                 m_variant);

    _variant_add(_("Native, fixed blocks"),
                 CS_MATRIX_NATIVE,
                 n_fill_types,
                 fill_types,
                 2, /* ed_flag */
                 NULL,
                 _b_mat_vec_p_l_native_fixed,
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
                   _mat_vec_p_l_native_omp_atomic,
                   _b_mat_vec_p_l_native_omp_atomic,
                   NULL,
                   n_variants,
                   &n_variants_max,
                   m_variant);

#endif

      if (numbering->type == CS_NUMBERING_VECTORIZE)
        _variant_add(_("Native, vectorized"),
                     CS_MATRIX_NATIVE,
                     n_fill_types,
                     fill_types,
                     2, /* ed_flag */
                     _mat_vec_p_l_native_vector,
                     NULL,
                     NULL,
                     n_variants,
                     &n_variants_max,
                     m_variant);

    }

  }

  if (type_filter[CS_MATRIX_CSR]) {

    _variant_add(_("CSR"),
                 CS_MATRIX_CSR,
                 n_fill_types,
                 fill_types,
                 2, /* ed_flag */
                 _mat_vec_p_l_csr,
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
                 _mat_vec_p_l_msr,
                 _b_mat_vec_p_l_msr,
                 NULL,
                 n_variants,
                 &n_variants_max,
                 m_variant);

    _variant_add(_("MSR, generic"),
                 CS_MATRIX_MSR,
                 n_fill_types,
                 fill_types,
                 2, /* ed_flag */
                 NULL,
                 _b_mat_vec_p_l_msr_generic,
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
                 _mat_vec_p_l_msr_mkl,
                 NULL,
                 NULL,
                 n_variants,
                 &n_variants_max,
                 m_variant);

#endif /* defined(HAVE_MKL) */

    _variant_add(_("MSR, OpenMP scheduling"),
                 CS_MATRIX_MSR,
                 n_fill_types,
                 fill_types,
                 2, /* ed_flag */
                 _mat_vec_p_l_msr_omp_sched,
                 NULL,
                 NULL,
                 n_variants,
                 &n_variants_max,
                 m_variant);

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
 * \brief Select the sparse matrix-vector product function to be used by a
 * matrix variant for a given fill type.
 *
 * Currently, possible variant functions are:
 *
 *   CS_MATRIX_NATIVE  (all fill types)
 *     standard
 *     fixed           (for CS_MATRIX_??_BLOCK_D or CS_MATRIX_??_BLOCK_D_SYM)
 *     omp             (for OpenMP with compatible numbering)
 *     vector          (For vector machine with compatible numbering)
 *
 *   CS_MATRIX_CSR     (for CS_MATRIX_SCALAR or CS_MATRIX_SCALAR_SYM)
 *     standard
 *     mkl             (with MKL)
 *
 *   CS_MATRIX_CSR_SYM (for CS_MATRIX_SCALAR_SYM)
 *     standard
 *     mkl             (with MKL)
 *
 *   CS_MATRIX_MSR     (all fill types except CS_MATRIX_33_BLOCK)
 *     standard
 *     fixed           (for CS_MATRIX_??_BLOCK_D or CS_MATRIX_??_BLOCK_D_SYM)
 *     mkl             (with MKL, for CS_MATRIX_SCALAR or CS_MATRIX_SCALAR_SYM)
 *
 * parameters:
 *   mv        <-> Pointer to matrix variant
 *   numbering <-- mesh numbering info, or NULL
 *   fill type <-- matrix fill type to merge from
 *   ed_flag   <-- 0: with diagonal only, 1 exclude only; 2; both
 *   func_name <-- function type name
 */
/*----------------------------------------------------------------------------*/

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
                               mv->vector_multiply);

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Merge a functions to a matrix variant from another variant sharing
 * the same structure.
 *
 * Functions from the structure to merge for the selected fill type are
 * assigned to the main variant.
 *
 * This can be useful when tuning has been done separately for different fill
 * types, and the resulting selected structure is identical.
 *
 * \param[in, out]  mv         pointer to matrix variant
 * \param[in]       mv_merge   pointer to matrix variant to merge
 * \param[in]       fill_type  matrix fill type to merge from
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_variant_merge(cs_matrix_variant_t        *mv,
                        const cs_matrix_variant_t  *mv_merge,
                        cs_matrix_fill_type_t       fill_type)
{
  if (mv->type == mv_merge->type) {
    for (int i = 0; i < 2; i++) {
      mv->vector_multiply[fill_type][i]
        = mv_merge->vector_multiply[fill_type][i];
      mv->matrix_vector_cost[fill_type][i][0]
        = mv_merge->matrix_vector_cost[fill_type][i][0];
      mv->matrix_vector_cost[fill_type][i][1]
        = mv_merge->matrix_vector_cost[fill_type][i][1];
    }
    mv->matrix_assign_cost[fill_type] = mv_merge->matrix_assign_cost[fill_type];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the type associated with a matrix variant.
 *
 * \param[in]  mv  pointer to matrix variant structure
 */
/*----------------------------------------------------------------------------*/

cs_matrix_type_t
cs_matrix_variant_type(const cs_matrix_variant_t  *mv)
{
  return mv->type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Test local matrix.vector product operations.
 *
 * \param[in]  n_rows       local number of rows
 * \param[in]  n_cols_ext   number of local + ghost columns
 * \param[in]  n_edges      local number of (undirected) graph edges
 * \param[in]  edges        edges (symmetric row <-> column) connectivity
 * \param[in]  halo         cell halo structure
 * \param[in]  numbering    vectorization or thread-related numbering info,
 *                          or NULL
 *
 * \return  pointer to tuning results structure
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_variant_test(cs_lnum_t              n_rows,
                       cs_lnum_t              n_cols_ext,
                       cs_lnum_t              n_edges,
                       const cs_lnum_2_t     *edges,
                       const cs_halo_t       *halo,
                       const cs_numbering_t  *numbering)
{
  int  n_variants = 0;
  bool type_filter[CS_MATRIX_N_TYPES] = {true, true, true, true};
  cs_matrix_fill_type_t  fill_types[] = {CS_MATRIX_SCALAR,
                                         CS_MATRIX_SCALAR_SYM,
                                         CS_MATRIX_BLOCK_D,
                                         CS_MATRIX_BLOCK_D_66,
                                         CS_MATRIX_BLOCK_D_SYM,
                                         CS_MATRIX_BLOCK};
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
                n_rows,
                n_cols_ext,
                n_edges,
                edges,
                halo,
                numbering,
                m_variant);

  n_variants = 0;
  BFT_FREE(m_variant);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
