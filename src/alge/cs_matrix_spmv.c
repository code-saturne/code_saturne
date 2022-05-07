/*============================================================================
 * Sparse Matrix-vector multiplication kernels.
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
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_matrix.h"
#include "cs_matrix_priv.h"
#include "cs_matrix_spmv.h"

#if defined (HAVE_CUDA)
#include "cs_matrix_spmv_cuda.h"
#endif

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*! \file cs_matrix_spmv.c
 *
 * \brief Sparse Matrix SpMV kernels.
 */
/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* Cache line multiple, in cs_real_t units */

static const cs_lnum_t _cs_cl = (CS_CL_SIZE/8);

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* Note that most types are declared in cs_matrix_priv.h.
   only those only handled here are declared here. */

/*============================================================================
 *  Global variables
 *============================================================================*/

#if defined (HAVE_MKL)

static char _no_exclude_diag_error_str[]
  = N_("Matrix product variant using function %s\n"
       "does not handle case with excluded diagonal.");

#endif

/*============================================================================
 * Private function definitions
- *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute matrix-vector product for one dense block: y[i] = a[i].x[i]
 *
 * Vectors and blocks may be larger than their useful size, to
 * improve data alignment.
 *
 * parameters:
 *   b_id     <-- block id
 *   b_size   <-- block size
 *   b_size_2 <-- block size squared
 *   a        <-- pointer to block matrixes array (usually matrix diagonal)
 *   x        <-- multipliying vector values
 *   y        --> resulting vector
 *----------------------------------------------------------------------------*/

static inline void
_dense_b_ax(cs_lnum_t         b_id,
            cs_lnum_t         b_size,
            cs_lnum_t         b_size_2,
            const cs_real_t   a[restrict],
            const cs_real_t   x[restrict],
            cs_real_t         y[restrict])
{
  cs_lnum_t   ii, jj;

  for (ii = 0; ii < b_size; ii++) {
    y[b_id*b_size + ii] = 0.;
    for (jj = 0; jj < b_size; jj++)
      y[b_id*b_size + ii]
        +=   a[b_id*b_size_2 + ii*b_size + jj]
           * x[b_id*b_size + jj];
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
 *   b_i      <-- block id for i
 *   b_j      <-- block id for j
 *   b_ij     <-- block id for matrix ij position
 *   b_size   <-- block size
 *   b_size_2 <-- block size squared
 *   a        <-- pointer to block matrixes array (usually extra-diagonal)
 *   x        <-- multipliying vector values
 *   y        --> resulting vector
 *----------------------------------------------------------------------------*/

static inline void
_dense_eb_ax_add(cs_lnum_t        b_i,
                 cs_lnum_t        b_j,
                 cs_lnum_t        b_ij,
                 cs_lnum_t        b_size,
                 cs_lnum_t        b_size_2,
                 const cs_real_t  a[restrict],
                 const cs_real_t  x[restrict],
                 cs_real_t        y[restrict])
{
  cs_lnum_t   ii, jj;

  for (ii = 0; ii < b_size; ii++) {
    for (jj = 0; jj < b_size; jj++)
      y[b_i*b_size + ii]
        +=   a[b_ij*b_size_2 + ii*b_size + jj]
           * x[b_j*b_size + jj];
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
 *   da       <-- pointer to coefficients array (usually matrix diagonal)
 *   x        <-- multipliying vector values
 *   y        --> resulting vector
 *   n_elts   <-- array size
 *   b_size   <-- block size
 *----------------------------------------------------------------------------*/

static inline void
_b_diag_vec_p_l(const cs_real_t  da[restrict],
                const cs_real_t  x[restrict],
                cs_real_t        y[restrict],
                cs_lnum_t        n_elts,
                cs_lnum_t        b_size)
{
  if (da != NULL) {
    cs_lnum_t  b_size_2 = b_size*b_size;
#   pragma omp parallel for  if(n_elts > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_elts; ii++)
      _dense_b_ax(ii, b_size, b_size_2, da, x, y);
  }
  else {
#   pragma omp parallel for  if(n_elts*b_size > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_elts*b_size; ii++)
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
 *   b_size   <-- block size
 *   b_size_2 <-- block size squared
 *----------------------------------------------------------------------------*/

static inline void
_b_zero_range(cs_real_t  y[restrict],
              cs_lnum_t  start_id,
              cs_lnum_t  end_id,
              cs_lnum_t  b_size)
{
  cs_lnum_t  ii;

# pragma omp parallel for  if((end_id-start_id)*b_size > CS_THR_MIN)
  for (ii = start_id*b_size; ii < end_id*b_size; ii++)
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
 * Start synchronization of ghost values prior to matrix.vector product
 *
 * parameters:
 *   matrix        <-- pointer to matrix structure
 *   x             <-> multipliying vector values (ghost values updated)
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
 *   x             <-> multipliying vector values (ghost values updated)
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
 * Matrix.vector product y = A.x with native matrix.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_p_l_native(const cs_matrix_t  *matrix,
                    bool                exclude_diag,
                    bool                sync,
                    cs_real_t           x[restrict],
                    cs_real_t           y[restrict])
{
  cs_lnum_t  ii, jj, face_id;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_real_t  *restrict xa = mc->e_val;

  /* Initialize ghost cell communication */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _diag_vec_p_l(mc->d_val, x, y, ms->n_rows);
    _zero_range(y, ms->n_rows, ms->n_cols_ext);
  }
  else
    _zero_range(y, 0, ms->n_cols_ext);

  /* Finalize ghost cell comunication if overlap used */

  if (hs != NULL)
    cs_halo_sync_wait(matrix->halo, x, hs);

  /* non-diagonal terms */

  if (mc->e_val != NULL) {

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
 * Matrix.vector product y = A.x with native matrix.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_b_mat_vec_p_l_native(const cs_matrix_t  *matrix,
                      bool                exclude_diag,
                      bool                sync,
                      cs_real_t           x[restrict],
                      cs_real_t           y[restrict])
{
  cs_lnum_t  ii, jj, kk, face_id;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_real_t  *restrict xa = mc->e_val;
  const cs_lnum_t db_size = matrix->db_size;

  /* Initialize ghost cell communication */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _b_diag_vec_p_l(mc->d_val, x, y, ms->n_rows, db_size);
    _b_zero_range(y, ms->n_rows, ms->n_cols_ext, db_size);
  }
  else
    _b_zero_range(y, 0, ms->n_cols_ext, db_size);

  /* Finalize ghost cell comunication if overlap used */

  if (hs != NULL)
    _pre_vector_multiply_sync_x_end(matrix, hs, x);

  /* non-diagonal terms */

  if (mc->e_val != NULL) {

    const cs_lnum_2_t *restrict face_cel_p = ms->edges;

    if (mc->symmetric) {

      for (face_id = 0; face_id < ms->n_edges; face_id++) {
        ii = face_cel_p[face_id][0];
        jj = face_cel_p[face_id][1];
        for (kk = 0; kk < db_size; kk++) {
          y[ii*db_size + kk] += xa[face_id] * x[jj*db_size + kk];
          y[jj*db_size + kk] += xa[face_id] * x[ii*db_size + kk];
        }
      }
    }
    else {

      for (face_id = 0; face_id < ms->n_edges; face_id++) {
        ii = face_cel_p[face_id][0];
        jj = face_cel_p[face_id][1];
        for (kk = 0; kk < db_size; kk++) {
          y[ii*db_size + kk] += xa[2*face_id]     * x[jj*db_size + kk];
          y[jj*db_size + kk] += xa[2*face_id + 1] * x[ii*db_size + kk];
        }
      }

    }

  }

}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with native matrix.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_bb_mat_vec_p_l_native(const cs_matrix_t  *matrix,
                       bool                exclude_diag,
                       bool                sync,
                       cs_real_t           x[restrict],
                       cs_real_t           y[restrict])
{
  cs_lnum_t  ii, jj, face_id;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_real_t  *restrict xa = mc->e_val;
  const cs_lnum_t  db_size = matrix->db_size;
  const cs_lnum_t  eb_size = matrix->eb_size;
  const cs_lnum_t  eb_size_2 = eb_size*eb_size;

  /* Initialize ghost cell communication */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _b_diag_vec_p_l(mc->d_val, x, y, ms->n_rows, db_size);
    _b_zero_range(y, ms->n_rows, ms->n_cols_ext, db_size);
  }
  else
    _b_zero_range(y, 0, ms->n_cols_ext, db_size);

  /* Finalize ghost cell comunication if overlap used */

  if (hs != NULL)
    _pre_vector_multiply_sync_x_end(matrix, hs, x);

  /* non-diagonal terms */

  if (mc->e_val != NULL) {

    const cs_lnum_2_t *restrict face_cel_p = ms->edges;

    if (mc->symmetric) {

      for (face_id = 0; face_id < ms->n_edges; face_id++) {
        ii = face_cel_p[face_id][0];
        jj = face_cel_p[face_id][1];
        _dense_eb_ax_add(ii, jj, face_id, eb_size, eb_size_2, xa, x, y);
        _dense_eb_ax_add(jj, ii, face_id, eb_size, eb_size_2, xa, x, y);
      }
    }
    else {

      for (face_id = 0; face_id < ms->n_edges; face_id++) {
        ii = face_cel_p[face_id][0];
        jj = face_cel_p[face_id][1];
        _dense_eb_ax_add(ii, jj, 2*face_id, eb_size, eb_size_2, xa, x, y);
        _dense_eb_ax_add(jj, ii, 2*face_id + 1, eb_size, eb_size_2, xa, x, y);
      }

    }

  }

}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with native matrix.
 *
 * This variant uses a fixed 3x3 block, for better compiler optimization.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_3_3_mat_vec_p_l_native(const cs_matrix_t  *matrix,
                        bool                exclude_diag,
                        bool                sync,
                        cs_real_t           x[restrict],
                        cs_real_t           y[restrict])
{
  cs_lnum_t  ii, jj, kk, face_id;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_real_t  *restrict xa = mc->e_val;

  assert(matrix->db_size == 3);

  /* Initialize ghost cell communication */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _3_3_diag_vec_p_l(mc->d_val, x, y, ms->n_rows);
    _3_3_zero_range(y, ms->n_rows, ms->n_cols_ext);
  }
  else
    _3_3_zero_range(y, 0, ms->n_cols_ext);

  /* Finalize ghost cell comunication */

  if (hs != NULL)
    _pre_vector_multiply_sync_x_end(matrix, hs, x);

  /* non-diagonal terms */

  if (mc->e_val != NULL) {

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
 * Matrix.vector product y = A.x with native matrix.
 *
 * This variant uses a fixed 6x6 block, for better compiler optimization.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_6_6_mat_vec_p_l_native(const cs_matrix_t  *matrix,
                        bool                exclude_diag,
                        bool                sync,
                        cs_real_t           x[restrict],
                        cs_real_t           y[restrict])
{
  cs_lnum_t  ii, jj, kk, face_id;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_real_t  *restrict xa = mc->e_val;

  assert(matrix->db_size == 6);

  /* Initialize ghost cell communication */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _6_6_diag_vec_p_l(mc->d_val, x, y, ms->n_rows);
    _6_6_zero_range(y, ms->n_rows, ms->n_cols_ext);
  }
  else
    _6_6_zero_range(y, 0, ms->n_cols_ext);

  /* Finalize ghost cell comunication if overlap used */

  if (hs != NULL)
    _pre_vector_multiply_sync_x_end(matrix, hs, x);

  /* non-diagonal terms */

  if (mc->e_val != NULL) {

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
 * Matrix.vector product y = A.x with native matrix, blocked version.
 *
 * This variant uses fixed block size variants for common cases.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_b_mat_vec_p_l_native_fixed(const cs_matrix_t  *matrix,
                            bool                exclude_diag,
                            bool                sync,
                            cs_real_t           x[restrict],
                            cs_real_t           y[restrict])
{
  if (matrix->db_size == 3)
    _3_3_mat_vec_p_l_native(matrix, exclude_diag, sync, x, y);

  else if (matrix->db_size == 6)
    _6_6_mat_vec_p_l_native(matrix, exclude_diag, sync, x, y);

  else
    _b_mat_vec_p_l_native(matrix, exclude_diag, sync, x, y);
}

#if defined(HAVE_OPENMP) /* OpenMP variants */

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with native matrix.
 *
 * parameters:
 *   matrix       <-- Pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> Multipliying vector values
 *   y            --> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_p_l_native_omp(const cs_matrix_t  *matrix,
                        bool                exclude_diag,
                        bool                sync,
                        cs_real_t           x[restrict],
                        cs_real_t           y[restrict])
{
  const int n_threads = matrix->numbering->n_threads;
  const int n_groups = matrix->numbering->n_groups;
  const cs_lnum_t *group_index = matrix->numbering->group_index;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;
  const cs_real_t  *restrict xa = mc->e_val;

  assert(matrix->numbering->type == CS_NUMBERING_THREADS);

  /* Initialize ghost cell communication */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _diag_vec_p_l(mc->d_val, x, y, ms->n_rows);
    _zero_range(y, ms->n_rows, ms->n_cols_ext);
  }
  else
    _zero_range(y, 0, ms->n_cols_ext);

  /* Finalize ghost cell comunication if overlap used */

  if (hs != NULL)
    cs_halo_sync_wait(matrix->halo, x, hs);

  /* non-diagonal terms */

  if (mc->e_val != NULL) {

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
 * Matrix.vector product y = A.x with native matrix, blocked version
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_b_mat_vec_p_l_native_omp(const cs_matrix_t  *matrix,
                          bool                exclude_diag,
                          bool                sync,
                          cs_real_t           x[restrict],
                          cs_real_t           y[restrict])
{
  const cs_lnum_t db_size = matrix->db_size;

  const int n_threads = matrix->numbering->n_threads;
  const int n_groups = matrix->numbering->n_groups;
  const cs_lnum_t *group_index = matrix->numbering->group_index;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;
  const cs_real_t  *restrict xa = mc->e_val;

  assert(matrix->numbering->type == CS_NUMBERING_THREADS);

  /* Initialize ghost cell communication */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _b_diag_vec_p_l(mc->d_val, x, y, ms->n_rows, db_size);
    _b_zero_range(y, ms->n_rows, ms->n_cols_ext, db_size);
  }
  else
    _b_zero_range(y, 0, ms->n_cols_ext, db_size);

  /* Finalize ghost cell comunication if overlap used */

  if (hs != NULL)
    _pre_vector_multiply_sync_x_end(matrix, hs, x);

  /* non-diagonal terms */

  if (mc->e_val != NULL) {

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
            for (cs_lnum_t kk = 0; kk < db_size; kk++) {
              y[ii*db_size + kk] += xa[face_id] * x[jj*db_size + kk];
              y[jj*db_size + kk] += xa[face_id] * x[ii*db_size + kk];
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
            for (cs_lnum_t kk = 0; kk < db_size; kk++) {
              y[ii*db_size + kk] += xa[2*face_id]     * x[jj*db_size + kk];
              y[jj*db_size + kk] += xa[2*face_id + 1] * x[ii*db_size + kk];
            }
          }
        }
      }

    }

  }
}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with native matrix.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_p_l_native_omp_atomic(const cs_matrix_t  *matrix,
                               bool                exclude_diag,
                               bool                sync,
                               cs_real_t           x[restrict],
                               cs_real_t           y[restrict])
{
  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;
  const cs_real_t  *restrict xa = mc->e_val;

  /* Initialize ghost cell communication */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _diag_vec_p_l(mc->d_val, x, y, ms->n_rows);
    _zero_range(y, ms->n_rows, ms->n_cols_ext);
  }
  else
    _zero_range(y, 0, ms->n_cols_ext);

  /* Finalize ghost cell comunication if overlap used */

  if (hs != NULL)
    cs_halo_sync_wait(matrix->halo, x, hs);

  /* non-diagonal terms */

  if (mc->e_val != NULL) {

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
 * Matrix.vector product y = A.x with native matrix, blocked version
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_b_mat_vec_p_l_native_omp_atomic(const cs_matrix_t  *matrix,
                                 bool                exclude_diag,
                                 bool                sync,
                                 cs_real_t           x[restrict],
                                 cs_real_t           y[restrict])
{
  const cs_lnum_t db_size = matrix->db_size;

  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;
  const cs_real_t  *restrict xa = mc->e_val;

  /* Initialize ghost cell communication */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _b_diag_vec_p_l(mc->d_val, x, y, ms->n_rows, db_size);
    _b_zero_range(y, ms->n_rows, ms->n_cols_ext, db_size);
  }
  else
    _b_zero_range(y, 0, ms->n_cols_ext, db_size);

  /* Finalize ghost cell comunication if overlap used */

  if (hs != NULL)
    _pre_vector_multiply_sync_x_end(matrix, hs, x);

  /* non-diagonal terms */

  if (mc->e_val != NULL) {

    const cs_lnum_2_t *restrict face_cel_p = ms->edges;

    if (mc->symmetric) {

#     pragma omp parallel for
      for (cs_lnum_t face_id = 0; face_id < ms->n_edges; face_id++) {
        cs_lnum_t ii = face_cel_p[face_id][0];
        cs_lnum_t jj = face_cel_p[face_id][1];
        for (cs_lnum_t kk = 0; kk < db_size; kk++) {
#         pragma omp atomic
          y[ii*db_size + kk] += xa[face_id] * x[jj*db_size + kk];
#         pragma omp atomic
          y[jj*db_size + kk] += xa[face_id] * x[ii*db_size + kk];
        }
      }

    }
    else {

#     pragma omp parallel for
      for (cs_lnum_t face_id = 0; face_id < ms->n_edges; face_id++) {
        cs_lnum_t ii = face_cel_p[face_id][0];
        cs_lnum_t jj = face_cel_p[face_id][1];
        for (cs_lnum_t kk = 0; kk < db_size; kk++) {
#         pragma omp atomic
          y[ii*db_size + kk] += xa[2*face_id]   * x[jj*db_size + kk];
#         pragma omp atomic
          y[jj*db_size + kk] += xa[2*face_id+1] * x[ii*db_size + kk];
        }
      }

    }

  }
}

#endif /* defined(HAVE_OPENMP) */

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with native matrix.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_p_l_native_vector(const cs_matrix_t  *matrix,
                           bool                exclude_diag,
                           bool                sync,
                           cs_real_t           x[restrict],
                           cs_real_t           y[restrict])
{
  cs_lnum_t  ii, jj, face_id;
  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;
  const cs_real_t  *restrict xa = mc->e_val;

  assert(matrix->numbering->type == CS_NUMBERING_VECTORIZE);

  /* Initialize ghost cell communication */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;

  /* Diagonal part of matrix.vector product */

  if (! exclude_diag) {
    _diag_vec_p_l(mc->d_val, x, y, ms->n_rows);
    _zero_range(y, ms->n_rows, ms->n_cols_ext);
  }
  else
    _zero_range(y, 0, ms->n_cols_ext);

  /* Finalize ghost cell comunication if overlap used */

  if (hs != NULL)
    cs_halo_sync_wait(matrix->halo, x, hs);

  /* non-diagonal terms */

  if (mc->e_val != NULL) {

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
 * Matrix.vector product y = A.x with CSR matrix.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_p_l_csr(const cs_matrix_t  *matrix,
                 bool                exclude_diag,
                 bool                sync,
                 cs_real_t          *restrict x,
                 cs_real_t          *restrict y)
{
  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_csr_t  *mc = matrix->coeffs;
  cs_lnum_t  n_rows = ms->n_rows;

  /* Ghost cell communication */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;
  if (hs != NULL)
    cs_halo_sync_wait(matrix->halo, x, hs);

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
_mat_vec_p_l_csr_mkl(const cs_matrix_t  *matrix,
                     bool                exclude_diag,
                     bool                sync,
                     cs_real_t          *restrict x,
                     cs_real_t          *restrict y)
{
  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_csr_t  *mc = matrix->coeffs;

  /* Ghost cell communication */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;
  if (hs != NULL)
    cs_halo_sync_wait(matrix->halo, x, hs);

  /* MKL call */

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
 * Matrix.vector product y = A.x with MSR matrix.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_p_l_msr(const cs_matrix_t  *matrix,
                 bool                exclude_diag,
                 bool                sync,
                 cs_real_t          *restrict x,
                 cs_real_t          *restrict y)
{
  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_lnum_t  n_rows = ms->n_rows;

  const cs_lnum_t  *e_col_id = ms->e.col_id;
  const cs_lnum_t  *e_row_index = ms->e.row_index;

  /* Ghost cell communication */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;
  if (hs != NULL)
    cs_halo_sync_wait(matrix->halo, x, hs);

  /* Standard case */

  if (!exclude_diag && mc->d_val != NULL) {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      const cs_lnum_t *restrict col_id = e_col_id + e_row_index[ii];
      const cs_real_t *restrict m_row = mc->e_val + e_row_index[ii];
      cs_lnum_t n_cols = e_row_index[ii+1] - e_row_index[ii];
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

      const cs_lnum_t *restrict col_id = e_col_id + e_row_index[ii];
      const cs_real_t *restrict m_row = mc->e_val + e_row_index[ii];
      cs_lnum_t n_cols = e_row_index[ii+1] - e_row_index[ii];
      cs_real_t sii = 0.0;

      for (cs_lnum_t jj = 0; jj < n_cols; jj++)
        sii += (m_row[jj]*x[col_id[jj]]);

      y[ii] = sii;

    }
  }

}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with MSR matrix.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_p_l_msr_omp_sched(const cs_matrix_t  *matrix,
                           bool                exclude_diag,
                           bool                sync,
                           cs_real_t          *restrict x,
                           cs_real_t          *restrict y)
{
  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_lnum_t  n_rows = ms->n_rows;

  const cs_lnum_t  *e_col_id = ms->e.col_id;
  const cs_lnum_t  *e_row_index = ms->e.row_index;

  /* Ghost cell communication */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;
  if (hs != NULL)
    cs_halo_sync_wait(matrix->halo, x, hs);

  /* Standard case */

  if (!exclude_diag && mc->d_val != NULL) {

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {

      cs_lnum_t n_s_rows = cs_align(n_rows * 0.9, _cs_cl);
      if (n_s_rows > n_rows)
        n_s_rows = n_rows;

#     pragma omp for nowait
      for (cs_lnum_t ii = 0; ii < n_s_rows; ii++) {

        const cs_lnum_t *restrict col_id = e_col_id + e_row_index[ii];
        const cs_real_t *restrict m_row = mc->e_val + e_row_index[ii];
        cs_lnum_t n_cols = e_row_index[ii+1] - e_row_index[ii];
        cs_real_t sii = 0.0;

        for (cs_lnum_t jj = 0; jj < n_cols; jj++)
          sii += (m_row[jj]*x[col_id[jj]]);

        y[ii] = sii + mc->d_val[ii]*x[ii];

      }

#     pragma omp for schedule(dynamic, _cs_cl)
      for (cs_lnum_t ii = n_s_rows; ii < n_rows; ii++) {

        const cs_lnum_t *restrict col_id = e_col_id + e_row_index[ii];
        const cs_real_t *restrict m_row = mc->e_val + e_row_index[ii];
        cs_lnum_t n_cols = e_row_index[ii+1] - e_row_index[ii];
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

      cs_lnum_t n_s_rows = cs_align(n_rows * 0.9, _cs_cl);
      if (n_s_rows > n_rows)
        n_s_rows = n_rows;

#     pragma omp for
      for (cs_lnum_t ii = 0; ii < n_s_rows; ii++) {

        const cs_lnum_t *restrict col_id = e_col_id + e_row_index[ii];
        const cs_real_t *restrict m_row = mc->e_val + e_row_index[ii];
        cs_lnum_t n_cols = e_row_index[ii+1] - e_row_index[ii];
        cs_real_t sii = 0.0;

        for (cs_lnum_t jj = 0; jj < n_cols; jj++)
          sii += (m_row[jj]*x[col_id[jj]]);

        y[ii] = sii;

      }

#     pragma omp for schedule(dynamic, _cs_cl)
      for (cs_lnum_t ii = n_s_rows; ii < n_rows; ii++) {

        const cs_lnum_t *restrict col_id = e_col_id + e_row_index[ii];
        const cs_real_t *restrict m_row = mc->e_val + e_row_index[ii];
        cs_lnum_t n_cols = e_row_index[ii+1] - e_row_index[ii];
        cs_real_t sii = 0.0;

        for (cs_lnum_t jj = 0; jj < n_cols; jj++)
          sii += (m_row[jj]*x[col_id[jj]]);

        y[ii] = sii;

      }

    }

  }

}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with MSR matrix, blocked version.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_b_mat_vec_p_l_msr_generic(const cs_matrix_t  *matrix,
                           bool                exclude_diag,
                           bool                sync,
                           cs_real_t           x[restrict],
                           cs_real_t           y[restrict])
{
  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_lnum_t  n_rows = ms->n_rows;
  const cs_lnum_t  db_size = matrix->db_size;
  const cs_lnum_t  db_size_2 = db_size*db_size;

  const cs_lnum_t  *e_col_id = ms->e.col_id;
  const cs_lnum_t  *e_row_index = ms->e.row_index;

  /* Ghost cell communication */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;
  if (hs != NULL)
    _pre_vector_multiply_sync_x_end(matrix, hs, x);

  /* Standard case */

  if (!exclude_diag && mc->d_val != NULL) {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      const cs_lnum_t *restrict col_id = e_col_id + e_row_index[ii];
      const cs_real_t *restrict m_row = mc->e_val + e_row_index[ii];
      cs_lnum_t n_cols = e_row_index[ii+1] - e_row_index[ii];

      _dense_b_ax(ii, db_size, db_size_2, mc->d_val, x, y);

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        for (cs_lnum_t kk = 0; kk < db_size; kk++) {
          y[ii*db_size + kk]
            += (m_row[jj]*x[col_id[jj]*db_size + kk]);
        }
      }

    }

  }

  /* Exclude diagonal */

  else {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      const cs_lnum_t *restrict col_id = e_col_id + e_row_index[ii];
      const cs_real_t *restrict m_row = mc->e_val + e_row_index[ii];
      cs_lnum_t n_cols = e_row_index[ii+1] - e_row_index[ii];

      for (cs_lnum_t kk = 0; kk < db_size; kk++)
        y[ii*db_size + kk] = 0.;

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        for (cs_lnum_t kk = 0; kk < db_size; kk++) {
          y[ii*db_size + kk]
            += (m_row[jj]*x[col_id[jj]*db_size + kk]);
        }
      }

    }
  }

}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with MSR matrix, 3x3 blocked version.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_b_mat_vec_p_l_msr_3(const cs_matrix_t  *matrix,
                     bool                exclude_diag,
                     bool                sync,
                     cs_real_t          *restrict x,
                     cs_real_t          *restrict y)
{
  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_lnum_t  n_rows = ms->n_rows;

  const cs_lnum_t  *e_col_id = ms->e.col_id;
  const cs_lnum_t  *e_row_index = ms->e.row_index;

  assert(matrix->db_size == 3);

  /* Ghost cell communication */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;
  if (hs != NULL)
    _pre_vector_multiply_sync_x_end(matrix, hs, x);

  /* Standard case */

  if (!exclude_diag && mc->d_val != NULL) {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      const cs_lnum_t *restrict col_id = e_col_id + e_row_index[ii];
      const cs_real_t *restrict m_row = mc->e_val + e_row_index[ii];
      cs_lnum_t n_cols = e_row_index[ii+1] - e_row_index[ii];

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

      const cs_lnum_t *restrict col_id = e_col_id + e_row_index[ii];
      const cs_real_t *restrict m_row = mc->e_val + e_row_index[ii];
      cs_lnum_t n_cols = e_row_index[ii+1] - e_row_index[ii];

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
 * Matrix.vector product y = A.x with MSR matrix, 6x6 blocked version.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_b_mat_vec_p_l_msr_6(const cs_matrix_t  *matrix,
                     bool                exclude_diag,
                     bool                sync,
                     cs_real_t           x[restrict],
                     cs_real_t           y[restrict])
{
  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_lnum_t  n_rows = ms->n_rows;

  const cs_lnum_t  *e_col_id = ms->e.col_id;
  const cs_lnum_t  *e_row_index = ms->e.row_index;

  assert(matrix->db_size == 6);

  /* Ghost cell communication */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;
  if (hs != NULL)
    _pre_vector_multiply_sync_x_end(matrix, hs, x);

  /* Standard case */

  if (!exclude_diag && mc->d_val != NULL) {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      const cs_lnum_t *restrict col_id = e_col_id + e_row_index[ii];
      const cs_real_t *restrict m_row = mc->e_val + e_row_index[ii];
      cs_lnum_t n_cols = e_row_index[ii+1] - e_row_index[ii];

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

      const cs_lnum_t *restrict col_id = e_col_id + e_row_index[ii];
      const cs_real_t *restrict m_row = mc->e_val + e_row_index[ii];
      cs_lnum_t n_cols = e_row_index[ii+1] - e_row_index[ii];

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
 * Matrix.vector product y = A.x with MSR matrix, blocked version.
 *
 * This variant uses fixed block size variants for common cases.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_b_mat_vec_p_l_msr(const cs_matrix_t  *matrix,
                   bool                exclude_diag,
                   bool                sync,
                   cs_real_t           x[restrict],
                   cs_real_t           y[restrict])
{
  if (matrix->db_size == 3)
    _b_mat_vec_p_l_msr_3(matrix, exclude_diag, sync, x, y);

  else if (matrix->db_size == 6)
    _b_mat_vec_p_l_msr_6(matrix, exclude_diag, sync, x, y);

  else
    _b_mat_vec_p_l_msr_generic(matrix, exclude_diag, sync, x, y);
}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with MSR matrix, blocked version.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_bb_mat_vec_p_l_msr_3(const cs_matrix_t  *matrix,
                      bool                exclude_diag,
                      bool                sync,
                      cs_real_t           x[restrict],
                      cs_real_t           y[restrict])
{
  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_lnum_t  n_rows = ms->n_rows;

  const cs_lnum_t  *e_col_id = ms->e.col_id;
  const cs_lnum_t  *e_row_index = ms->e.row_index;

  /* Ghost cell communication */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;
  if (hs != NULL)
    _pre_vector_multiply_sync_x_end(matrix, hs, x);

  /* Standard case */

  if (!exclude_diag && mc->d_val != NULL) {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      const cs_lnum_t *restrict col_id = e_col_id + e_row_index[ii];
      const cs_real_t *restrict m_row =  mc->e_val + (e_row_index[ii]*9);
      cs_lnum_t n_cols = e_row_index[ii+1] - e_row_index[ii];

      _dense_3_3_ax(ii, mc->d_val, x, y);

      cs_real_t * _y = y + ii*3;

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        _y[0] += (  m_row[jj*9]         * x[col_id[jj]*3]
                  + m_row[jj*9 + 1]     * x[col_id[jj]*3 + 1]
                  + m_row[jj*9 + 2]     * x[col_id[jj]*3 + 2]);
        _y[1] += (  m_row[jj*9 + 3]     * x[col_id[jj]*3]
                  + m_row[jj*9 + 3 + 1] * x[col_id[jj]*3 + 1]
                  + m_row[jj*9 + 3 + 2] * x[col_id[jj]*3 + 2]);
        _y[2] += (  m_row[jj*9 + 6]     * x[col_id[jj]*3]
                  + m_row[jj*9 + 6 + 1] * x[col_id[jj]*3 + 1]
                  + m_row[jj*9 + 6 + 2] * x[col_id[jj]*3 + 2]);
      }

    }

  }

  /* Exclude diagonal */

  else {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      const cs_lnum_t *restrict col_id = e_col_id + e_row_index[ii];
      const cs_real_t *restrict m_row =  mc->e_val + (e_row_index[ii]*9);
      cs_lnum_t n_cols = e_row_index[ii+1] - e_row_index[ii];

      cs_real_t * _y = y + (ii*3);

      for (cs_lnum_t kk = 0; kk < 3; kk++)
        _y[kk] = 0.;

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        _y[0] += (  m_row[jj*9]         * x[col_id[jj]*3]
                  + m_row[jj*9 + 1]     * x[col_id[jj]*3 + 1]
                  + m_row[jj*9 + 2]     * x[col_id[jj]*3 + 2]);
        _y[1] += (  m_row[jj*9 + 3]     * x[col_id[jj]*3]
                  + m_row[jj*9 + 3 + 1] * x[col_id[jj]*3 + 1]
                  + m_row[jj*9 + 3 + 2] * x[col_id[jj]*3 + 2]);
        _y[2] += (  m_row[jj*9 + 6]     * x[col_id[jj]*3]
                  + m_row[jj*9 + 6 + 1] * x[col_id[jj]*3 + 1]
                  + m_row[jj*9 + 6 + 2] * x[col_id[jj]*3 + 2]);
      }

    }
  }

}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with MSR matrix, blocked version.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_bb_mat_vec_p_l_msr_generic(const cs_matrix_t  *matrix,
                            bool                exclude_diag,
                            bool                sync,
                            cs_real_t           x[restrict],
                            cs_real_t           y[restrict])
{
  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_lnum_t  n_rows = ms->n_rows;
  const cs_lnum_t  db_size = matrix->db_size;
  const cs_lnum_t  eb_size = matrix->eb_size;
  const cs_lnum_t  db_size_2 = db_size * db_size;
  const cs_lnum_t  eb_size_2 = eb_size * eb_size;

  const cs_lnum_t  *e_col_id = ms->e.col_id;
  const cs_lnum_t  *e_row_index = ms->e.row_index;

  /* Ghost cell communication */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;
  if (hs != NULL)
    _pre_vector_multiply_sync_x_end(matrix, hs, x);

  /* Standard case */

  if (!exclude_diag && mc->d_val != NULL) {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      const cs_lnum_t *restrict col_id = e_col_id + e_row_index[ii];
      const cs_real_t *restrict m_row =   mc->e_val
                                        + (e_row_index[ii]*eb_size_2);
      cs_lnum_t n_cols = e_row_index[ii+1] - e_row_index[ii];

      _dense_b_ax(ii, db_size, db_size_2, mc->d_val, x, y);

      cs_real_t * _y = y + (ii*db_size);

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        for (cs_lnum_t kk = 0; kk < db_size; kk++) {
          for (cs_lnum_t ll = 0; ll < db_size; ll++) {
            _y[kk] += (  m_row[jj*eb_size_2 + kk*eb_size + ll]
                       * x[col_id[jj]*db_size + ll]);
          }
        }
      }

    }

  }

  /* Exclude diagonal */

  else {

#   pragma omp parallel for  if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      const cs_lnum_t *restrict col_id = e_col_id + e_row_index[ii];
      const cs_real_t *restrict m_row =   mc->e_val
                                        + (e_row_index[ii]*eb_size_2);
      cs_lnum_t n_cols = e_row_index[ii+1] - e_row_index[ii];

      cs_real_t * _y = y + (ii*db_size);

      for (cs_lnum_t kk = 0; kk < db_size; kk++)
        _y[kk] = 0.;

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        for (cs_lnum_t kk = 0; kk < db_size; kk++) {
          for (cs_lnum_t ll = 0; ll < db_size; ll++) {
            _y[kk] += (  m_row[jj*eb_size_2 + kk*eb_size + ll]
                       * x[col_id[jj]*db_size + ll]);
          }
        }
      }

    }
  }

}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with MSR matrix, blocked version.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_bb_mat_vec_p_l_msr(const cs_matrix_t  *matrix,
                    bool                exclude_diag,
                    bool                sync,
                    cs_real_t           x[restrict],
                    cs_real_t           y[restrict])
{
  if (matrix->eb_size == 3)
    _bb_mat_vec_p_l_msr_3(matrix, exclude_diag, sync, x, y);

  else
    _bb_mat_vec_p_l_msr_generic(matrix, exclude_diag, sync, x, y);
}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with MSR matrix, using MKL
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- pointer to matrix structure
 *   hs           <-- halo state: if non-NULL, call cs_halo_sync_wait
 *                    locally (possibly allowing computation/communication
 *                    overlap)
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

#if defined (HAVE_MKL)

static void
_mat_vec_p_l_msr_mkl(const cs_matrix_t  *matrix,
                     bool                exclude_diag,
                     bool                sync,
                     cs_real_t           x[restrict],
                     cs_real_t           y[restrict])
{
  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  /* Ghost cell communication */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;
  if (hs != NULL)
    cs_halo_sync_wait(matrix->halo, x, hs);

  /* Call MKL function */

  int n_rows = ms->n_rows;
  char transa[] = "n";

  mkl_cspblas_dcsrgemv(transa,
                       &n_rows,
                       mc->e_val,
                       ms->e.row_index,
                       ms->e.col_id,
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
 * Matrix.vector product y = A.x with MSR matrix.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_p_l_dist(const cs_matrix_t  *matrix,
                  bool                exclude_diag,
                  bool                sync,
                  cs_real_t          *restrict x,
                  cs_real_t          *restrict y)
{
  /* Initialize halo synchronization */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;

  /* Compute local part */

  _mat_vec_p_l_msr(matrix,
                   exclude_diag,
                   false,
                   x,
                   y);

  /* Finalize halo communication */

  if (hs != NULL)
    cs_halo_sync_wait(matrix->halo, x, hs);

  /* Compute distant part */

  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_lnum_t  n_h_rows = ms->h.n_rows;

  /* Standard case */

  if (n_h_rows > 0) {

    const cs_lnum_t  *h_col_id = ms->h.col_id;
    const cs_lnum_t  *h_row_index = ms->h.row_index;
    const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

#   pragma omp parallel for  if(n_h_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_h_rows; ii++) {

      const cs_lnum_t *restrict col_id = h_col_id + h_row_index[ii];
      const cs_real_t *restrict m_row = mc->h_val + h_row_index[ii];
      cs_lnum_t n_cols = h_row_index[ii+1] - h_row_index[ii];
      cs_real_t sii = 0.0;

      for (cs_lnum_t jj = 0; jj < n_cols; jj++)
        sii += (m_row[jj]*x[col_id[jj]]);

      y[ii] += sii;

    }

  }
}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with MSR matrix.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_p_l_dist_omp_sched(const cs_matrix_t  *matrix,
                           bool                 exclude_diag,
                           bool                 sync,
                           cs_real_t           *restrict x,
                           cs_real_t           *restrict y)
{
  /* Initialize halo synchronization */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;

  /* Compute local part */

  _mat_vec_p_l_msr_omp_sched(matrix,
                             exclude_diag,
                             false,
                             x,
                             y);

  /* Finalize halo communication */

  if (hs != NULL)
    cs_halo_sync_wait(matrix->halo, x, hs);

  /* Compute distant part */

  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_lnum_t  n_h_rows = ms->h.n_rows;

  if (n_h_rows > 0) {

    const cs_lnum_t  *h_col_id = ms->h.col_id;
    const cs_lnum_t  *h_row_index = ms->h.row_index;
    const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

#   pragma omp parallel for  if(n_h_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_h_rows; ii++) {

      const cs_lnum_t *restrict col_id = h_col_id + h_row_index[ii];
      const cs_real_t *restrict m_row = mc->h_val + h_row_index[ii];
      cs_lnum_t n_cols = h_row_index[ii+1] - h_row_index[ii];
      cs_real_t sii = 0.0;

      for (cs_lnum_t jj = 0; jj < n_cols; jj++)
        sii += (m_row[jj]*x[col_id[jj]]);

      y[ii] += sii;

    }

  }
}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with MSR matrix, blocked version.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_b_mat_vec_p_l_dist_generic(const cs_matrix_t  *matrix,
                            bool                exclude_diag,
                            bool                sync,
                            cs_real_t           x[restrict],
                            cs_real_t           y[restrict])
{
  /* Initialize halo synchronization */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;

  /* Compute local part */

  _b_mat_vec_p_l_msr_generic(matrix,
                             exclude_diag,
                             false,
                             x,
                             y);

  /* Finalize halo communication */

  if (hs != NULL)
    _pre_vector_multiply_sync_x_end(matrix, hs, x);

  /* Compute distant part */

  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_lnum_t  n_h_rows = ms->h.n_rows;

  if (n_h_rows > 0) {

    const cs_lnum_t  db_size = matrix->db_size;

    const cs_lnum_t  *h_col_id = ms->h.col_id;
    const cs_lnum_t  *h_row_index = ms->h.row_index;

#   pragma omp parallel for  if(n_h_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_h_rows; ii++) {

      const cs_lnum_t *restrict col_id = h_col_id + h_row_index[ii];
      const cs_real_t *restrict m_row = mc->h_val + h_row_index[ii];
      cs_lnum_t n_cols = h_row_index[ii+1] - h_row_index[ii];

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        for (cs_lnum_t kk = 0; kk < db_size; kk++) {
          y[ii*db_size + kk]
            += (m_row[jj]*x[col_id[jj]*db_size + kk]);
        }
      }

    }

  }
}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with MSR matrix, 3x3 blocked version.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_b_mat_vec_p_l_dist_3(const cs_matrix_t  *matrix,
                      bool                exclude_diag,
                      bool                sync,
                      cs_real_t          *restrict x,
                      cs_real_t          *restrict y)
{
  /* Initialize halo synchronization */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;

  /* Compute local part */

  _b_mat_vec_p_l_msr_3(matrix,
                       exclude_diag,
                       false,
                       x,
                       y);

  /* Finalize halo communication */

  if (hs != NULL)
    _pre_vector_multiply_sync_x_end(matrix, hs, x);

  /* Compute distant part */

  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_lnum_t  n_h_rows = ms->h.n_rows;

  if (n_h_rows > 0) {

    const cs_lnum_t  *h_col_id = ms->h.col_id;
    const cs_lnum_t  *h_row_index = ms->h.row_index;

#   pragma omp parallel for  if(n_h_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_h_rows; ii++) {

      const cs_lnum_t *restrict col_id = h_col_id + h_row_index[ii];
      const cs_real_t *restrict m_row = mc->h_val + h_row_index[ii];
      cs_lnum_t n_cols = h_row_index[ii+1] - h_row_index[ii];

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        for (cs_lnum_t kk = 0; kk < 3; kk++)
          y[ii*3 + kk] += (m_row[jj]*x[col_id[jj]*3 + kk]);
      }

    }

  }
}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with MSR matrix, 6x6 blocked version.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_b_mat_vec_p_l_dist_6(const cs_matrix_t  *matrix,
                      bool                exclude_diag,
                      bool                sync,
                      cs_real_t           x[restrict],
                      cs_real_t           y[restrict])
{
  /* Initialize halo synchronization */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;

  /* Compute local part */

  _b_mat_vec_p_l_msr_6(matrix,
                       exclude_diag,
                       false,
                       x,
                       y);

  /* Finalize halo communication */

  if (hs != NULL)
    _pre_vector_multiply_sync_x_end(matrix, hs, x);

  /* Compute distant part */

  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_lnum_t  n_h_rows = ms->h.n_rows;

  if (n_h_rows > 0) {

    const cs_lnum_t  *h_col_id = ms->h.col_id;
    const cs_lnum_t  *h_row_index = ms->h.row_index;

#   pragma omp parallel for  if(n_h_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_h_rows; ii++) {

      const cs_lnum_t *restrict col_id = h_col_id + h_row_index[ii];
      const cs_real_t *restrict m_row = mc->h_val + h_row_index[ii];
      cs_lnum_t n_cols = h_row_index[ii+1] - h_row_index[ii];

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        for (cs_lnum_t kk = 0; kk < 6; kk++)
          y[ii*6 + kk] += (m_row[jj]*x[col_id[jj]*6 + kk]);
      }

    }

  }
}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with MSR matrix, blocked version.
 *
 * This variant uses fixed block size variants for common cases.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_b_mat_vec_p_l_dist(const cs_matrix_t  *matrix,
                    bool                exclude_diag,
                    bool                sync,
                    cs_real_t           x[restrict],
                    cs_real_t           y[restrict])
{
  if (matrix->db_size == 3)
    _b_mat_vec_p_l_dist_3(matrix, exclude_diag, sync, x, y);

  else if (matrix->db_size == 6)
    _b_mat_vec_p_l_dist_6(matrix, exclude_diag, sync, x, y);

  else
    _b_mat_vec_p_l_dist_generic(matrix, exclude_diag, sync, x, y);
}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with MSR matrix, blocked version.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_bb_mat_vec_p_l_dist_3(const cs_matrix_t  *matrix,
                       bool                exclude_diag,
                       bool                sync,
                       cs_real_t           x[restrict],
                       cs_real_t           y[restrict])
{
  /* Initialize halo synchronization */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;

  /* Compute local part */

  _bb_mat_vec_p_l_msr_3(matrix,
                        exclude_diag,
                        false,
                        x,
                        y);

  /* Finalize halo communication */

  if (hs != NULL)
    _pre_vector_multiply_sync_x_end(matrix, hs, x);

  /* Compute distant part */

  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_lnum_t  n_h_rows = ms->h.n_rows;

  if (n_h_rows > 0) {

    const cs_lnum_t  *h_col_id = ms->h.col_id;
    const cs_lnum_t  *h_row_index = ms->h.row_index;

#   pragma omp parallel for  if(n_h_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_h_rows; ii++) {

      const cs_lnum_t *restrict col_id = h_col_id + h_row_index[ii];
      const cs_real_t *restrict m_row =  mc->h_val + (h_row_index[ii]*9);
      cs_lnum_t n_cols = h_row_index[ii+1] - h_row_index[ii];

      cs_real_t * _y = y + ii*3;

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        _y[0] += (  m_row[jj*9]         * x[col_id[jj]*3]
                  + m_row[jj*9 + 1]     * x[col_id[jj]*3 + 1]
                  + m_row[jj*9 + 2]     * x[col_id[jj]*3 + 2]);
        _y[1] += (  m_row[jj*9 + 3]     * x[col_id[jj]*3]
                  + m_row[jj*9 + 3 + 1] * x[col_id[jj]*3 + 1]
                  + m_row[jj*9 + 3 + 2] * x[col_id[jj]*3 + 2]);
        _y[2] += (  m_row[jj*9 + 6]     * x[col_id[jj]*3]
                  + m_row[jj*9 + 6 + 1] * x[col_id[jj]*3 + 1]
                  + m_row[jj*9 + 6 + 2] * x[col_id[jj]*3 + 2]);
      }

    }

  }
}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with MSR matrix, blocked version.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_bb_mat_vec_p_l_dist_generic(const cs_matrix_t  *matrix,
                             bool                exclude_diag,
                             bool                sync,
                             cs_real_t           x[restrict],
                             cs_real_t           y[restrict])
{
  /* Initialize halo synchronization */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;

  /* Compute local part */

  _b_mat_vec_p_l_msr_generic(matrix,
                             exclude_diag,
                             false,
                             x,
                             y);

  /* Finalize halo communication */

  if (hs != NULL)
    _pre_vector_multiply_sync_x_end(matrix, hs, x);

  /* Compute distant part */

  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_lnum_t  n_h_rows = ms->h.n_rows;

  if (n_h_rows > 0) {

    const cs_lnum_t  *h_col_id = ms->h.col_id;
    const cs_lnum_t  *h_row_index = ms->h.row_index;

    const cs_lnum_t  b_size = matrix->eb_size;
    const cs_lnum_t  b_size_2 = b_size * b_size;

#   pragma omp parallel for  if(n_h_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_h_rows; ii++) {

      const cs_lnum_t *restrict col_id = h_col_id + h_row_index[ii];
      const cs_real_t *restrict m_row
        = mc->h_val + (h_row_index[ii]*b_size_2);
      cs_lnum_t n_cols = h_row_index[ii+1] - h_row_index[ii];

      cs_real_t * _y = y + (ii*b_size);

      for (cs_lnum_t jj = 0; jj < n_cols; jj++) {
        for (cs_lnum_t kk = 0; kk < b_size; kk++) {
          for (cs_lnum_t ll = 0; ll < b_size; ll++) {
            _y[kk] += (  m_row[jj*b_size_2 + kk*b_size + ll]
                       * x[col_id[jj]*b_size + ll]);
          }
        }
      }

    }

  }
}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with MSR matrix, blocked version.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true,
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_bb_mat_vec_p_l_dist(const cs_matrix_t  *matrix,
                    bool                exclude_diag,
                    bool                sync,
                    cs_real_t           x[restrict],
                    cs_real_t           y[restrict])
{
  if (matrix->eb_size == 3)
    _bb_mat_vec_p_l_dist_3(matrix, exclude_diag, sync, x, y);

  else
    _bb_mat_vec_p_l_dist_generic(matrix, exclude_diag, sync, x, y);
}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with MSR matrix, using MKL
 *
 * parameters:
 *   exclude_diag <-- exclude diagonal if true
 *   matrix       <-- pointer to matrix structure
 *   hs           <-- halo state: if non-NULL, call cs_halo_sync_wait
 *                    locally (possibly allowing computation/communication
 *                    overlap)
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

#if defined (HAVE_MKL)

static void
_mat_vec_p_l_dist_mkl(const cs_matrix_t  *matrix,
                      bool                exclude_diag,
                      bool                sync,
                      cs_real_t           x[restrict],
                      cs_real_t           y[restrict])
{
  /* Initialize halo synchronization */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : NULL;

  /* Compute local part */

  _mat_vec_p_l_msr_mkl(matrix,
                       exclude_diag,
                       sync,
                       x,
                       y);

  /* Finalize halo communication */

  if (hs != NULL)
    cs_halo_sync_wait(matrix->halo, x, hs);

  /* Compute distant part */

  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_lnum_t  n_h_rows = ms->h.n_rows;

  /* Standard case */

  if (n_h_rows > 0) {

    const cs_lnum_t  *h_col_id = ms->h.col_id;
    const cs_lnum_t  *h_row_index = ms->h.row_index;
    const cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

#   pragma omp parallel for  if(n_h_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_h_rows; ii++) {

      const cs_lnum_t *restrict col_id = h_col_id + h_row_index[ii];
      const cs_real_t *restrict m_row = mc->h_val + h_row_index[ii];
      cs_lnum_t n_cols = h_row_index[ii+1] - h_row_index[ii];
      cs_real_t sii = 0.0;

      for (cs_lnum_t jj = 0; jj < n_cols; jj++)
        sii += (m_row[jj]*x[col_id[jj]]);

      y[ii] += sii;

    }

  }
}

#endif /* defined (HAVE_MKL) */

#if defined(HAVE_ACCEL)

/*----------------------------------------------------------------------------
 * Select the sparse matrix-vector product function to be used on device by a
 * matrix or variant for a given fill type.
 *
 * Currently, possible variant functions are:
 *
 *   CS_MATRIX_CSR     (for CS_MATRIX_SCALAR or CS_MATRIX_SCALAR_SYM)
 *     default
 *     cuda            (CUDA-accelerated)
 *     cusparse        (with cuSPARSE)
 *
 *   CS_MATRIX_MSR
 *     default
 *     cuda            (CUDA-accelerated)
 *     cusparse        (with cuSPARSE)
 *
 * parameters:
 *   m_type      <--  Matrix type
 *   fill type   <--  matrix fill type to merge from
 *   spmv_type   <--  SpMV operation type (full or sub-matrix)
 *                    (all types if CS_MATRIX_SPMV_N_TYPES)
 *   numbering   <--  mesh numbering structure, or NULL
 *   func_name   <--  function type name, or NULL for default
 *   spmv        <->  multiplication function array
 *
 * returns:
 *   0 for success, 1 for incompatible function, 2 for compatible
 *   function not available in current build
 *----------------------------------------------------------------------------*/

static int
_matrix_spmv_set_func_d(cs_matrix_type_t             m_type,
                        cs_matrix_fill_type_t        fill_type,
                        cs_matrix_spmv_type_t        spmv_type,
                        const cs_numbering_t        *numbering,
                        const char                  *func_name,
                        cs_matrix_vector_product_t  *spmv[CS_MATRIX_SPMV_N_TYPES])
{
  const char *_func_name = func_name;

#if defined(HAVE_CUDA)

#if !defined(HAVE_CUSPARSE_GENERIC_API)
const char s_cuda[] = "cuda";
#endif

#if defined(HAVE_CUSPARSE)

const char s_cusparse[] = "cusparse";
const char *default_name = s_cusparse;

#if !defined(HAVE_CUSPARSE_GENERIC_API)
 if (fill_type >= CS_MATRIX_BLOCK_D && fill_type < CS_MATRIX_BLOCK)
   default_name = s_cuda;
#endif

#else

const char *default_name = s_cuda;

#endif /* defined(HAVE_CUSPARSE) */

#else

const char default_name[] = "not_implemented";

#endif /* defined(HAVE_CUDA) */

 if (_func_name == NULL)
   _func_name = default_name;
 else if (!strcmp(func_name, "default"))
   _func_name = default_name;

  char spmv_xy_hd[CS_MATRIX_SPMV_N_TYPES];

  int retcode = cs_matrix_spmv_set_func(m_type,
                                        fill_type,
                                        spmv_type,
                                        numbering,
                                        _func_name,
                                        spmv,
                                        spmv_xy_hd);

  return retcode;
}

#endif /* defined(HAVE_ACCEL) */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Assign default sparse matrix-vector product functions for a given matrix.
 *
 * parameters:
 *   m <-> pointer to matrix structure
 *----------------------------------------------------------------------------*/

void
cs_matrix_spmv_set_defaults(cs_matrix_t  *m)
{
  char spmv_xy_hd[CS_MATRIX_SPMV_N_TYPES];

  if (m->destroy_adaptor != NULL)
    m->destroy_adaptor(m);

  for (cs_matrix_fill_type_t mft = 0; mft < CS_MATRIX_N_FILL_TYPES; mft++) {
    for (cs_matrix_spmv_type_t spmv_type = 0;
         spmv_type < CS_MATRIX_SPMV_N_TYPES;
         spmv_type++) {
      cs_matrix_spmv_set_func(m->type,
                              mft,
                              spmv_type,
                              m->numbering,
                              NULL, /* func_name */
                              m->vector_multiply[mft],
                              spmv_xy_hd);
#if defined(HAVE_ACCEL)
      _matrix_spmv_set_func_d(m->type,
                              mft,
                              spmv_type,
                              m->numbering,
                              NULL, /* func_name */
                              m->vector_multiply_d[mft]);
#endif
    }
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
 *     baseline
 *     omp             (for OpenMP with compatible numbering)
 *     omp_atomic      (for OpenMP with atomic add)
 *     vector          (For vector machine with compatible numbering)
 *
 *   CS_MATRIX_CSR     (for CS_MATRIX_SCALAR or CS_MATRIX_SCALAR_SYM)
 *     default
 *     mkl             (with MKL)
 *     cuda            (CUDA-accelerated)
 *     cusparse        (with cuSPARSE)
 *
 *   CS_MATRIX_MSR
 *     default
 *     omp_sched       (Improved OpenMP scheduling, for CS_MATRIX_SCALAR*)
 *     mkl             (with MKL, for CS_MATRIX_SCALAR or CS_MATRIX_SCALAR_SYM)
 *     cuda            (CUDA-accelerated)
 *     cusparse        (with cuSPARSE)
 *
 *   CS_MATRIX_DIST
 *     default
 *     omp_sched       (Improved OpenMP scheduling, for CS_MATRIX_SCALAR*)
 *     mkl             (with MKL, for CS_MATRIX_SCALAR or CS_MATRIX_SCALAR_SYM)
 *
 * parameters:
 *   m_type      <--  Matrix type
 *   fill type   <--  matrix fill type to merge from
 *   spmv_type   <--  SpMV operation type (full or sub-matrix)
 *                    (all types if CS_MATRIX_SPMV_N_TYPES)
 *   numbering   <--  mesh numbering structure, or NULL
 *   func_name   <--  function type name, or NULL for default
 *   spmv        <->  multiplication function array
 *   spmv_xy_hd  <->  multiplication function x and y host/device location
 *
 * returns:
 *   0 for success, 1 for incompatible function, 2 for compatible
 *   function not available in current build
 *----------------------------------------------------------------------------*/

int
cs_matrix_spmv_set_func(cs_matrix_type_t             m_type,
                        cs_matrix_fill_type_t        fill_type,
                        cs_matrix_spmv_type_t        spmv_type,
                        const cs_numbering_t        *numbering,
                        const char                  *func_name,
                        cs_matrix_vector_product_t  *spmv[CS_MATRIX_SPMV_N_TYPES],
                        char                         spmv_xy_hd[CS_MATRIX_SPMV_N_TYPES])
{
  int retcode = 1;
  int standard = 0;

  cs_matrix_vector_product_t *_spmv[CS_MATRIX_SPMV_N_TYPES]
    = {NULL, NULL};

  char _spmv_xy_hd[CS_MATRIX_SPMV_N_TYPES] = {'h', 'h'};

  if (func_name == NULL)
    standard = 2;
  else if (!strcmp(func_name, "default"))
    standard = 2;
  else if (!strcmp(func_name, "baseline"))
    standard = 1;

  switch(m_type) {

  /* Native
     ------ */

  case CS_MATRIX_NATIVE:

    if (standard > 0) { /* standard or default */

      switch(fill_type) {
      case CS_MATRIX_SCALAR:
      case CS_MATRIX_SCALAR_SYM:
        _spmv[0] = _mat_vec_p_l_native;
        _spmv[1] = _mat_vec_p_l_native;
        break;
      case CS_MATRIX_BLOCK_D:
      case CS_MATRIX_BLOCK_D_66:
      case CS_MATRIX_BLOCK_D_SYM:
        _spmv[0] = _b_mat_vec_p_l_native_fixed;
        _spmv[1] = _b_mat_vec_p_l_native_fixed;
        break;
      case CS_MATRIX_BLOCK:
        _spmv[0] = _bb_mat_vec_p_l_native;
        _spmv[1] = _bb_mat_vec_p_l_native;
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
              _spmv[0] = _mat_vec_p_l_native_omp;
              _spmv[1] = _mat_vec_p_l_native_omp;
            }
#endif
            if (numbering->type == CS_NUMBERING_VECTORIZE) {
              _spmv[0] = _mat_vec_p_l_native_vector;
              _spmv[1] = _mat_vec_p_l_native_vector;
            }
          }
          break;
        case CS_MATRIX_BLOCK_D:
        case CS_MATRIX_BLOCK_D_66:
        case CS_MATRIX_BLOCK_D_SYM:
          if (numbering != NULL) {
#if defined(HAVE_OPENMP)
            if (numbering->type == CS_NUMBERING_THREADS) {
              _spmv[0] = _b_mat_vec_p_l_native_omp;
              _spmv[1] = _b_mat_vec_p_l_native_omp;
            }
#endif
          }
          break;
        default:
          break;
        }
      }

    }

    else if (!strcmp(func_name, "omp")) {
#if defined(HAVE_OPENMP)
      if (numbering != NULL) {
        if (numbering->type == CS_NUMBERING_THREADS) {
          switch(fill_type) {
          case CS_MATRIX_SCALAR:
          case CS_MATRIX_SCALAR_SYM:
            _spmv[0] = _mat_vec_p_l_native_omp;
            _spmv[1] = _mat_vec_p_l_native_omp;
            break;
          case CS_MATRIX_BLOCK_D:
          case CS_MATRIX_BLOCK_D_66:
          case CS_MATRIX_BLOCK_D_SYM:
            _spmv[0] = _b_mat_vec_p_l_native_omp;
            _spmv[1] = _b_mat_vec_p_l_native_omp;
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

    else if (!strcmp(func_name, "omp_atomic")) {
#if defined(HAVE_OPENMP)
      switch(fill_type) {
      case CS_MATRIX_SCALAR:
      case CS_MATRIX_SCALAR_SYM:
        _spmv[0] = _mat_vec_p_l_native_omp_atomic;
        _spmv[1] = _mat_vec_p_l_native_omp_atomic;
        break;
      case CS_MATRIX_BLOCK_D:
      case CS_MATRIX_BLOCK_D_66:
      case CS_MATRIX_BLOCK_D_SYM:
        _spmv[0] = _b_mat_vec_p_l_native_omp_atomic;
        _spmv[1] = _b_mat_vec_p_l_native_omp_atomic;
        break;
      default:
        break;
      }
#else
      retcode = 2;
#endif
    }

    else if (!strcmp(func_name, "vector")) {
      switch(fill_type) {
      case CS_MATRIX_SCALAR:
      case CS_MATRIX_SCALAR_SYM:
        _spmv[0] = _mat_vec_p_l_native_vector;
        _spmv[1] = _mat_vec_p_l_native_vector;
        break;
      default:
        break;
      }
    }

    break;

  /* CSR
     --- */

  case CS_MATRIX_CSR:

    switch(fill_type) {
    case CS_MATRIX_SCALAR:
    case CS_MATRIX_SCALAR_SYM:
      if (standard > 0) {
        _spmv[0] = _mat_vec_p_l_csr;
        _spmv[1] = _mat_vec_p_l_csr;
      }
      else if (!strcmp(func_name, "mkl")) {
#if defined(HAVE_MKL)
        _spmv[0] = _mat_vec_p_l_csr_mkl;
        _spmv[1] = _mat_vec_p_l_csr_mkl;
#else
        retcode = 2;
#endif
      }
      else if (!strcmp(func_name, "cuda")) {
#if defined(HAVE_CUDA)
        _spmv[0] = cs_matrix_spmv_cuda_csr;
        _spmv[1] = cs_matrix_spmv_cuda_csr;
        _spmv_xy_hd[0] = 'd';
        _spmv_xy_hd[1] = 'd';
#else
        retcode = 2;
#endif
      }
      else if (!strcmp(func_name, "cusparse")) {
#if defined(HAVE_CUSPARSE)
        _spmv[0] = (cs_matrix_vector_product_t *)
                     cs_matrix_spmv_cuda_csr_cusparse;
        _spmv[1] = (cs_matrix_vector_product_t *)
                     cs_matrix_spmv_cuda_csr_cusparse;
        _spmv_xy_hd[0] = 'd';
        _spmv_xy_hd[1] = 'd';
#else
        retcode = 2;
#endif
      }
      break;
    default:
      break;
    }

    break;

  /* MSR
     --- */

  case CS_MATRIX_MSR:

    if (standard > 0) {
      switch(fill_type) {
      case CS_MATRIX_SCALAR:
      case CS_MATRIX_SCALAR_SYM:
        _spmv[0] = _mat_vec_p_l_msr;
        _spmv[1] = _mat_vec_p_l_msr;
        break;
      case CS_MATRIX_BLOCK_D:
      case CS_MATRIX_BLOCK_D_66:
      case CS_MATRIX_BLOCK_D_SYM:
        _spmv[0] = _b_mat_vec_p_l_msr;
        _spmv[1] = _b_mat_vec_p_l_msr;
        break;
      case CS_MATRIX_BLOCK:
        _spmv[0] = _bb_mat_vec_p_l_msr;
        _spmv[1] = _bb_mat_vec_p_l_msr;
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
        _spmv[0] = _mat_vec_p_l_msr_mkl;
        _spmv[1] = _mat_vec_p_l_msr_mkl;
        break;
      default:
        break;
      }
#else
      retcode = 2;
#endif
    }

    else if (!strcmp(func_name, "cuda")) {
#if defined(HAVE_CUDA)
      switch(fill_type) {
      case CS_MATRIX_SCALAR:
      case CS_MATRIX_SCALAR_SYM:
        _spmv[0] = cs_matrix_spmv_cuda_msr;
        _spmv[1] = cs_matrix_spmv_cuda_msr;
        _spmv_xy_hd[0] = 'd';
        _spmv_xy_hd[1] = 'd';
        break;
      case CS_MATRIX_BLOCK_D:
      case CS_MATRIX_BLOCK_D_66:
      case CS_MATRIX_BLOCK_D_SYM:
        _spmv[0] = cs_matrix_spmv_cuda_msr_b;
        _spmv[1] = cs_matrix_spmv_cuda_msr_b;
        _spmv_xy_hd[0] = 'd';
        _spmv_xy_hd[1] = 'd';
        break;
      default:
        break;
      }
#else
      retcode = 2;
#endif
    }

    else if (!strcmp(func_name, "cusparse")) {
#if defined(HAVE_CUSPARSE)
      switch(fill_type) {
      case CS_MATRIX_SCALAR:
      case CS_MATRIX_SCALAR_SYM:
        _spmv[0] = cs_matrix_spmv_cuda_msr_cusparse;
        _spmv[1] = cs_matrix_spmv_cuda_msr_cusparse;
        _spmv_xy_hd[0] = 'd';
        _spmv_xy_hd[1] = 'd';
        break;
      case CS_MATRIX_BLOCK_D:
      case CS_MATRIX_BLOCK_D_66:
      case CS_MATRIX_BLOCK_D_SYM:
#if defined(HAVE_CUSPARSE_GENERIC_API)
        _spmv[0] = cs_matrix_spmv_cuda_msr_b_cusparse;
        _spmv[1] = cs_matrix_spmv_cuda_msr_b_cusparse;
        _spmv_xy_hd[0] = 'd';
        _spmv_xy_hd[1] = 'd';
#endif
        break;
      case CS_MATRIX_BLOCK:
        _spmv[0] = cs_matrix_spmv_cuda_msr_bb_cusparse;
        _spmv[1] = cs_matrix_spmv_cuda_msr_bb_cusparse;
        _spmv_xy_hd[0] = 'd';
        _spmv_xy_hd[1] = 'd';
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
        _spmv[0] = _mat_vec_p_l_msr_omp_sched;
        _spmv[1] = _mat_vec_p_l_msr_omp_sched;
        break;
      default:
        break;
      }
    }

    break;

  /* Distributed
     ----------- */

  case CS_MATRIX_DIST:

    if (standard > 0) {
      switch(fill_type) {
      case CS_MATRIX_SCALAR:
      case CS_MATRIX_SCALAR_SYM:
        _spmv[0] = _mat_vec_p_l_dist;
        _spmv[1] = _mat_vec_p_l_dist;
        break;
      case CS_MATRIX_BLOCK_D:
      case CS_MATRIX_BLOCK_D_66:
      case CS_MATRIX_BLOCK_D_SYM:
        _spmv[0] = _b_mat_vec_p_l_dist;
        _spmv[1] = _b_mat_vec_p_l_dist;
        break;
      case CS_MATRIX_BLOCK:
        _spmv[0] = _bb_mat_vec_p_l_dist;
        _spmv[1] = _bb_mat_vec_p_l_dist;
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
        _spmv[0] = _mat_vec_p_l_dist_mkl;
        _spmv[1] = _mat_vec_p_l_dist_mkl;
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
        _spmv[0] = _mat_vec_p_l_dist_omp_sched;
        _spmv[1] = _mat_vec_p_l_dist_omp_sched;
        break;
      default:
        break;
      }
    }

    break;

  default:
    break;
  }

  if (spmv_type < CS_MATRIX_SPMV_N_TYPES) {
    if (_spmv[spmv_type] != NULL) {
      spmv[spmv_type] = _spmv[spmv_type];
      spmv_xy_hd[spmv_type] = _spmv_xy_hd[spmv_type];
      retcode = 0;
    }
  }
  else {
    for (cs_matrix_spmv_type_t i = 0; i < CS_MATRIX_SPMV_N_TYPES; i++) {
      if (_spmv[i] != NULL) {
        spmv[i] = _spmv[i];
        spmv_xy_hd[i] = _spmv_xy_hd[i];
        retcode = 0;
      }
    }
  }

  return retcode;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
