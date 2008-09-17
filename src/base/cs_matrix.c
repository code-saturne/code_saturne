/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
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
 * Sparse Matrix Representation and Operations.
 *============================================================================*/

/*
 * Notes:
 *
 * The aim of these structures and associated functions is to multiple:
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(__STDC_VERSION__)      /* size_t */
#if (__STDC_VERSION__ == 199901L)
#    include <stddef.h>
#  else
#    include <stdlib.h>
#  endif
#else
#include <stdlib.h>
#endif

#if defined(_CS_HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_defs.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_prototypes.h"
#include "cs_perio.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_matrix.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* Variant for Intel compiler and on Itanium only (optimized by BULL)
   (Use compile flag -DNO_BULL_OPTIM to switch back to general code) */

#if (defined(__INTEL_COMPILER) && defined(__ia64__) && !defined(NO_BULL_OPTIM))

#define IA64_OPTIM
#define IA64_OPTIM_L1_CACHE_SIZE (508)
#define IA64_OPTIM_MIN(a, b) (((a) < (b)) ? (a) : (b))

#endif

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* Formats currently supported:
 *
 *  - Native
 *  - Compressed Sparse Row (CSR)
 */

/* Function pointer types */
/*------------------------*/

typedef void
(cs_matrix_set_coeffs_t) (cs_matrix_t      *matrix,
                          cs_bool_t         symmetric,
                          const cs_real_t  *restrict da,
                          const cs_real_t  *restrict xa);

typedef void
(cs_matrix_release_coeffs_t) (cs_matrix_t  *matrix);

typedef void
(cs_matrix_get_diagonal_t) (const cs_matrix_t  *matrix,
                            cs_real_t          *restrict da);

typedef void
(cs_matrix_vector_product_t) (const cs_matrix_t  *matrix,
                              const cs_real_t    *restrict x,
                              cs_real_t          *restrict y);

typedef void
(cs_matrix_alpha_a_x_p_beta_y_t) (cs_real_t           alpha,
                                  cs_real_t           beta,
                                  const cs_matrix_t  *matrix,
                                  const cs_real_t    *restrict x,
                                  cs_real_t          *restrict y);

/*----------------------------------------------------------------------------
 * Local Structure Definitions
 *----------------------------------------------------------------------------*/

/* Native matrix structure representation */
/*----------------------------------------*/

/* Note: the members of this structure are already available through the top
 *       matrix structure, but are replicated here in case of future removal
 *       from the top structure (which would require computation/assignment of
 *       matrix coefficients in another form) */

typedef struct _cs_matrix_struct_native_t {

  cs_int_t           n_cells;       /* Local number of cells */
  cs_int_t           n_cells_ext;   /* Local number of participating cells
                                       (cells + ghost cells sharing a face) */
  cs_int_t           n_faces;       /* Local number of faces
                                       (for extra-diagonal terms */

  /* Pointers to shared arrays */

  const cs_int_t    *face_cell;     /* Face -> cells connectivity (1 to n) */

} cs_matrix_struct_native_t;

/* Native matrix coefficients */
/*----------------------------*/

typedef struct _cs_matrix_coeff_native_t {

  cs_bool_t         symmetric;      /* Symmetry indicator */

  /* Pointers to possibly shared arrays */

  const cs_real_t   *da;            /* Diagonal terms */
  const cs_real_t   *xa;            /* Extra-diagonal terms */

  /* Pointers to private arrays (NULL if shared) */

  cs_real_t         *_da;           /* Diagonal terms */
  cs_real_t         *_xa;           /* Extra-diagonal terms */

} cs_matrix_coeff_native_t;

/* CSR (Compressed Sparse Row) matrix structure representation */
/*-------------------------------------------------------------*/

typedef struct _cs_matrix_struct_csr_t {

  cs_bool_t         symmetric;        /* Symmetry indicator */
  cs_int_t          n_rows;           /* Local number of rows */
  cs_int_t          n_cols;           /* Local number of columns
                                         (> n_rows in case of ghost cells) */
  cs_int_t          n_cols_max;       /* Maximum number of nonzero values
                                         on a given row */

  /* Pointers to structure arrays and info (row_index, col_id) */

  cs_bool_t         direct_assembly;  /* True if each value corresponds to
                                         a unique face ; false if multiple
                                         faces contribute to the same
                                         value (i.e. we have split faces) */

  cs_int_t         *row_index;        /* Row index (0 to n-1) */
  cs_int_t         *col_id;           /* Column id (0 to n-1) */

  /* Pointers to optional arrays (built if needed) */

  cs_int_t         *diag_index;       /* Diagonal index (0 to n-1) for
                                         direct access to diagonal terms */

} cs_matrix_struct_csr_t;

/* CSR matrix coefficients representation */
/*----------------------------------------*/

typedef struct _cs_matrix_coeff_csr_t {

  int               n_prefetch_rows;  /* Number of rows at a time for which
                                         the x values in y = Ax should be
                                         prefetched (0 if no prefetch) */

  cs_real_t        *val;              /* Matrix coefficients */

  cs_real_t        *x_prefetch;       /* Prefetch array for x in y = Ax */

} cs_matrix_coeff_csr_t;

/* Structure associated with Matrix (representation-independent part) */
/*--------------------------------------------------------------------*/

struct _cs_matrix_t {

  cs_matrix_type_t       type;         /* Matrix storage and definition type */

  cs_bool_t              periodic;     /* Periodicity indicator */
  cs_bool_t              have_diag;    /* Has non-zero diagonal */

  cs_int_t               n_cells;      /* Local number of cells */
  cs_int_t               n_cells_ext;  /* Local number of participating cells
                                          (cells + ghost cells sharing a face) */
  cs_int_t               n_faces;      /* Local Number of mesh faces
                                          (necessary to affect coefficients) */

  /* Pointer to possibly shared data */

  const void            *structure;    /* Matrix structure */

  /* Pointer to private data */

  void                  *_structure;   /* Matrix structure */
  void                  *coeffs;       /* Matrix coefficients */

  /* Pointers to shared arrays from mesh structure
     (face->cell connectivity for coefficient assignment,
     local->local cell numbering for future info or renumbering,
     and halo) */

  const cs_int_t        *face_cell;    /* Face -> cells connectivity (1 to n) */
  const fvm_gnum_t      *cell_num;     /* Global cell numbers */
  const cs_halo_t       *halo;         /* Parallel or periodic halo */

  /* Function pointers */

  cs_matrix_set_coeffs_t        *set_coefficients;
  cs_matrix_release_coeffs_t    *release_coefficients;
  cs_matrix_get_diagonal_t      *get_diagonal;

  cs_matrix_vector_product_t      *vector_multiply;
  cs_matrix_alpha_a_x_p_beta_y_t  *alpha_a_x_p_beta_y;
};

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Short names for matrix types */

const char  *cs_matrix_type_name[] = {N_("native"),
                                      N_("CSR")};

/* Full names for matrix types */

const char  *cs_matrix_type_fullname[] = {N_("diagonal + faces"),
                                          N_("Compressed Sparse Row")};

static int _cs_glob_matrix_prefetch_rows = 2048;

static char _cs_glob_perio_ignore_error_str[]
  = N_("Matrix product with CS_PERIO_IGNORE rotation mode not yet\n"
       "implemented in this case, use cs_matrix_vector_multiply_nosync\n"
       "with an external halo synchronization, prededed by a backup and\n"
       "followed by a restoration of the rotation halo.");

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Descend binary tree for the ordering of a fvm_gnum (integer) array.
 *
 * parameters:
 *   number    <-- pointer to numbers of elements that should be ordered.
 *                 (if NULL, a default 1 to n numbering is considered)
 *   level     <-- level of the binary tree to descend
 *   n_elts    <-- number of elements in the binary tree to descend
 *   order     <-> ordering array
 *----------------------------------------------------------------------------*/

inline static void
_order_descend_tree(const cs_int_t  number[],
                    size_t          level,
                    const size_t    n_elts,
                    fvm_lnum_t      order[])
{
  size_t i_save, i1, i2, lv_cur;

  i_save = (size_t)(order[level]);

  while (level <= (n_elts/2)) {

    lv_cur = (2*level) + 1;

    if (lv_cur < n_elts - 1) {

      i1 = (size_t)(order[lv_cur+1]);
      i2 = (size_t)(order[lv_cur]);

      if (number[i1] > number[i2]) lv_cur++;
    }

    if (lv_cur >= n_elts) break;

    i1 = i_save;
    i2 = (size_t)(order[lv_cur]);

    if (number[i1] >= number[i2]) break;

    order[level] = order[lv_cur];
    level = lv_cur;

  }

  order[level] = i_save;
}

/*----------------------------------------------------------------------------
 * Order an array of global numbers.
 *
 * parameters:
 *   number   <-- array of entity numbers (if NULL, a default 1 to n
 *                numbering is considered)
 *   order    <-- pre-allocated ordering table
 *   n_elts   <-- number of elements considered
 *----------------------------------------------------------------------------*/

static void
_order_local(const cs_int_t  number[],
             cs_int_t        order[],
             const size_t    n_elts)
{
  size_t i;
  cs_int_t o_save;

  /* Initialize ordering array */

  for (i = 0 ; i < n_elts ; i++)
    order[i] = i;

  if (n_elts < 2)
    return;

  /* Create binary tree */

  i = (n_elts / 2) ;
  do {
    i--;
    _order_descend_tree(number, i, n_elts, order);
  } while (i > 0);

  /* Sort binary tree */

  for (i = n_elts - 1 ; i > 0 ; i--) {
    o_save   = order[0];
    order[0] = order[i];
    order[i] = o_save;
    _order_descend_tree(number, 0, i, order);
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
 *   n_cells     --> Local number of participating cells
 *   n_cells_ext --> Local number of cells + ghost cells sharing a face
 *   n_faces     --> Local number of faces
 *   face_cell   --> Face -> cells connectivity (1 to n)
 *
 * returns:
 *   pointer to allocated native matrix structure.
 *----------------------------------------------------------------------------*/

static cs_matrix_struct_native_t *
_create_struct_native(int              n_cells,
                      int              n_cells_ext,
                      int              n_faces,
                      const cs_int_t  *face_cell)
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
_destroy_struct_native(cs_matrix_struct_native_t **matrix)
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
 *   matrix    --> Pointer to matrix structure
 *   symmetric --> Indicates if extradiagonal values are symmetric
 *   da        --> Diagonal values
 *   xa        --> Extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_coeffs_native(cs_matrix_t      *matrix,
                   cs_bool_t         symmetric,
                   const cs_real_t  *da,
                   const cs_real_t  *xa)
{
  cs_matrix_coeff_native_t  *mc = matrix->coeffs;
  const cs_matrix_struct_native_t  *ms = matrix->structure;

  mc->symmetric = symmetric;

  /* Map or copy values */

  if (da != NULL && matrix->have_diag == true) {

    if (mc->_da == NULL)
      mc->da = da;
    else {
      memcpy(mc->_da, da, sizeof(cs_real_t) * ms->n_cells);
      mc->da = mc->_da;
    }

  }
  else {
    mc->da = NULL;
  }

  if (xa != NULL) {

    if (mc->_xa == NULL)
      mc->xa = xa;
    else {
      size_t xa_n_bytes = sizeof(cs_real_t) * ms->n_faces;
      if (symmetric)
        xa_n_bytes *= 2;
      memcpy(mc->_xa, xa, xa_n_bytes);
      mc->xa = mc->_xa;
    }

  }

}

/*----------------------------------------------------------------------------
 * Release native matrix coefficients.
 *
 * parameters:
 *   matrix --> Pointer to matrix structure
 *----------------------------------------------------------------------------*/

static void
_release_coeffs_native(cs_matrix_t  *matrix)
{
  cs_matrix_coeff_native_t  *mc = matrix->coeffs;

  if (mc !=NULL) {

    /* Unmap values */

    mc->da = NULL;
    mc->xa = NULL;

    /* Possibly allocated mc->_da and mc->_xa arrays are not freed
       here, but simply unmapped from mc->da and mc->xa;
       they are remapped when _set_coeffs_native() is used. */

  }

}

/*----------------------------------------------------------------------------
 * Get diagonal of native matrix.
 *
 * parameters:
 *   matrix --> Pointer to matrix structure
 *   da    <-- Diagonal (pre-allocated, size: n_cells)
 *----------------------------------------------------------------------------*/

static void
_get_diagonal_native(const cs_matrix_t  *matrix,
                     cs_real_t          *restrict da)
{
  cs_int_t  ii;
  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;
  cs_int_t  n_cells = ms->n_cells;

  if (mc->da != NULL) {

    for (ii = 0; ii < n_cells; ii++)
      da[ii] = mc->da[ii];

  }
  else {

    for (ii = 0; ii < n_cells; ii++)
      da[ii] = 0.0;

  }

}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with native matrix.
 *
 * parameters:
 *   matrix --> Pointer to matrix structure
 *   x      --> Multipliying vector values
 *   y      <-- Resulting vector
 *----------------------------------------------------------------------------*/

#if !defined(IA64_OPTIM)  /* Standard variant */

static void
_mat_vec_p_l_native(const cs_matrix_t  *matrix,
                    const cs_real_t    *restrict x,
                    cs_real_t          *restrict y)
{
  cs_int_t  ii, jj, face_id;
  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;
  const cs_real_t  *restrict da = mc->da;
  const cs_real_t  *restrict xa1 = mc->xa;
  const cs_real_t  *restrict xa2 = mc->xa + ms->n_faces;

  /* Tell IBM compiler not to alias */
#if defined(__xlc__)
#pragma disjoint(*x, *y, *da, *xa1, *xa2)
#endif

  /* Diagonal part of matrix.vector product */

  /* Note: also try with BLAS: DNDOT(n_cells, 1, y, 1, 1, da, x, 1, 1) */

  if (mc->da != NULL) {
    for (ii = 0; ii < ms->n_cells; ii++)
      y[ii] = da[ii] * x[ii];
  }
  else {
    for (ii = 0; ii < ms->n_cells; y[ii++] = 0.0);
  }

  for (ii = ms->n_cells; ii < ms->n_cells_ext; y[ii++] = 0.0);

  /* Note: parallel and periodic synchronization could be delayed to here */

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    if (mc->symmetric) {

      const cs_int_t *restrict face_cel_p = ms->face_cell;

      for (face_id = 0; face_id < ms->n_faces; face_id++) {
        ii = *face_cel_p++ - 1;
        jj = *face_cel_p++ - 1;
        y[ii] += xa1[face_id] * x[jj];
        y[jj] += xa1[face_id] * x[ii];
      }

    }
    else {

      const cs_int_t *restrict face_cel_p = ms->face_cell;

      for (face_id = 0; face_id < ms->n_faces; face_id++) {
        ii = *face_cel_p++ - 1;
        jj = *face_cel_p++ - 1;
        y[ii] += xa1[face_id] * x[jj];
        y[jj] += xa2[face_id] * x[ii];
      }

    }

  }

}

#else /* defined(IA64_OPTIM), IA64 optimized variant (optimization by BULL) */

static void
_mat_vec_p_l_native(const cs_matrix_t  *matrix,
                    const cs_real_t    *restrict x,
                    cs_real_t          *restrict y)
{
  cs_int_t  ii, ii_prev, kk, face_id, kk_max;
  cs_real_t y_it, y_it_prev;
  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;
  const cs_real_t  *restrict da = mc->da;
  const cs_real_t  *restrict xa1 = mc->xa;
  const cs_real_t  *restrict xa2 = mc->xa + ms->n_faces;

  /* Diagonal part of matrix.vector product */

  /* Note: also try with BLAS: DNDOT(n_cells, 1, y, 1, 1, da, x, 1, 1) */

  if (mc->da != NULL) {
    for (ii = 0; ii < ms->n_cells; ii++)
      y[ii] = da[ii] * x[ii];
  }
  else {
    for (ii = 0; ii < ms->n_cells; y[ii++] = 0.0);
  }

  for (ii = ms->n_cells; ii < ms->n_cells_ext; y[ii++] = 0.0);

  /* Note: parallel and periodic synchronization could be delayed to here */

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    /*
     * 1/ Split y[ii] and y[jj] computation into 2 loops to remove compiler
     *    data dependency assertion between y[ii] and y[jj].
     * 2/ keep index (*face_cel_p) in L1 cache from y[ii] loop to y[jj] loop
     *    and xa1 in L2 cache.
     * 3/ break high frequency occurence of data dependency from one iteration
     *    to another in y[ii] loop (nonzero matrix value on the same line ii).
     */

    if (mc->symmetric) {

      const cs_int_t *restrict face_cel_p = ms->face_cell;

      for (face_id = 0;
           face_id < ms->n_faces;
           face_id += IA64_OPTIM_L1_CACHE_SIZE) {

        kk_max = IA64_OPTIM_MIN((ms->n_faces - face_id),
                                IA64_OPTIM_L1_CACHE_SIZE);

        /* sub-loop to compute y[ii] += xa1[face_id] * x[jj] */

        ii_prev = face_cel_p[0] - 1;
        y_it_prev = y[ii_prev] + xa1[face_id] * x[face_cel_p[1] - 1];

        for (kk = 1; kk < kk_max; ++kk) {
          ii = face_cel_p[2*kk] - 1;
          /* y[ii] += xa1[face_id+kk] * x[jj]; */
          if(ii == ii_prev) {
            y_it = y_it_prev;
          }
          else {
            y_it = y[ii];
            y[ii_prev] = y_it_prev;
          }
          ii_prev = ii;
          y_it_prev = y_it + xa1[face_id+kk] * x[face_cel_p[2*kk+1] - 1];
        }
        y[ii] = y_it_prev;

        /* sub-loop to compute y[ii] += xa1[face_id] * x[jj] */

        for (kk = 0; kk < kk_max; ++kk) {
          y[face_cel_p[2*kk+1] - 1]
            += xa1[face_id+kk] * x[face_cel_p[2*kk] - 1];
        }
        face_cel_p += 2 * IA64_OPTIM_L1_CACHE_SIZE;
      }

    }
    else {

      const cs_int_t *restrict face_cel_p = ms->face_cell;

      for (face_id = 0;
           face_id < ms->n_faces;
           face_id+=IA64_OPTIM_L1_CACHE_SIZE) {

        kk_max = IA64_OPTIM_MIN((ms->n_faces - face_id),
                                IA64_OPTIM_L1_CACHE_SIZE);

        /* sub-loop to compute y[ii] += xa1[face_id] * x[jj] */

        ii_prev = face_cel_p[0] - 1;
        y_it_prev = y[ii_prev] + xa1[face_id] * x[face_cel_p[1] - 1];

        for (kk = 1; kk < kk_max; ++kk) {
          ii = face_cel_p[2*kk] - 1;
          /* y[ii] += xa1[face_id+i] * x[jj]; */
          if(ii == ii_prev) {
            y_it = y_it_prev;
          }
          else {
            y_it = y[ii];
            y[ii_prev] = y_it_prev;
          }
          ii_prev = ii;
          y_it_prev = y_it + xa1[face_id+kk] * x[face_cel_p[2*kk+1] - 1];
        }
        y[ii] = y_it_prev;

        /* sub-loop to compute y[ii] += xa2[face_id] * x[jj] */

        for (kk = 0; kk < kk_max; ++kk) {
          y[face_cel_p[2*kk+1] - 1]
            += xa2[face_id+kk] * x[face_cel_p[2*kk] - 1];
        }
        face_cel_p += 2 * IA64_OPTIM_L1_CACHE_SIZE;
      }

    }

  }

}

#endif /* defined(IA64_OPTIM) */

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = alpha.A.x + beta.y with native matrix.
 *
 * parameters:
 *   alpha  --> Scalar, alpha in alpha.A.x + beta.y
 *   beta   --> Scalar, beta in alpha.A.x + beta.y
 *   matrix --> Pointer to matrix structure
 *   x      --> Multipliying vector values
 *   y      <-> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_alpha_a_x_p_beta_y_native(cs_real_t           alpha,
                           cs_real_t           beta,
                           const cs_matrix_t  *matrix,
                           const cs_real_t    *restrict x,
                           cs_real_t          *restrict y)
{
  cs_int_t  ii, jj, face_id;
  const cs_matrix_struct_native_t  *ms = matrix->structure;
  const cs_matrix_coeff_native_t  *mc = matrix->coeffs;
  const cs_real_t  *restrict da = mc->da;
  const cs_real_t  *restrict xa1 = mc->xa;
  const cs_real_t  *restrict xa2 = mc->xa + ms->n_faces;

  /* Tell IBM compiler not to alias */
#if defined(__xlc__)
#pragma disjoint(*x, *y, *da, *xa1, *xa2)
#endif

  /* Diagonal part of matrix.vector product */

  if (mc->da != NULL) {
    for (ii = 0; ii < ms->n_cells; ii++)
      y[ii] = (alpha * da[ii] * x[ii]) + (beta * y[ii]);
  }
  else {
    for (ii = 0; ii < ms->n_cells; ii++)
      y[ii] *= beta;
  }

  for (ii = ms->n_cells; ii < ms->n_cells_ext; y[ii++] = 0.0);

  /* Note: parallel and periodic synchronization could be delayed to here */

  /* non-diagonal terms */

  if (mc->xa != NULL) {

    if (mc->symmetric) {

      const cs_int_t *restrict face_cel_p = ms->face_cell;

      for (face_id = 0; face_id < ms->n_faces; face_id++) {
        ii = *face_cel_p++ - 1;
        jj = *face_cel_p++ - 1;
        y[ii] += alpha * xa1[face_id] * x[jj];
        y[jj] += alpha * xa1[face_id] * x[ii];
      }

    }
    else {

      const cs_int_t *restrict face_cel_p = ms->face_cell;

      for (face_id = 0; face_id < ms->n_faces; face_id++) {
        ii = *face_cel_p++ - 1;
        jj = *face_cel_p++ - 1;
        y[ii] += alpha * xa1[face_id] * x[jj];
        y[jj] += alpha * xa2[face_id] * x[ii];
      }

    }

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
 *   symmetric   --> Indicates if symmetric variant should be used
 *   have_diag   --> Indicates if the diagonal is nonzero
 *                   (forced to true for symmetric variant)
 *   n_cells     --> Local number of participating cells
 *   n_cells_ext --> Local number of cells + ghost cells sharing a face
 *   n_faces     --> Local number of faces
 *   cell_num    --> Global cell numbers (1 to n)
 *   face_cell   --> Face -> cells connectivity (1 to n)
 *
 * returns:
 *   pointer to allocated CSR matrix structure.
 *----------------------------------------------------------------------------*/

static cs_matrix_struct_csr_t *
_create_struct_csr(cs_bool_t         symmetric,
                   cs_bool_t         have_diag,
                   int               n_cells,
                   int               n_cells_ext,
                   int               n_faces,
                   const cs_int_t   *face_cell)
{
  int n_cols_max;
  cs_int_t ii, jj, face_id;
  const cs_int_t *restrict face_cel_p;

  cs_int_t  diag_elts = 1;
  cs_int_t  *ccount = NULL;

  cs_matrix_struct_csr_t  *ms;

  /* Allocate and map */

  BFT_MALLOC(ms, 1, cs_matrix_struct_csr_t);

  ms->symmetric = symmetric;

  ms->n_rows = n_cells;
  ms->n_cols = n_cells_ext;

  ms->direct_assembly = true;

  BFT_MALLOC(ms->row_index, ms->n_rows + 1, cs_int_t);
  ms->row_index = ms->row_index;

  ms->diag_index = NULL; /* Diagonal index only built if required */

  /* Count number of nonzero elements per row */

  BFT_MALLOC(ccount, ms->n_cols, cs_int_t);

  if (have_diag == false)
    diag_elts = 0;

  for (ii = 0; ii < ms->n_rows; ii++)  /* count starting with diagonal terms */
    ccount[ii] = diag_elts;

  if (face_cell != NULL) {

    face_cel_p = face_cell;

    if (symmetric == false) {

      for (face_id = 0; face_id < n_faces; face_id++) {
        ii = *face_cel_p++ - 1;
        jj = *face_cel_p++ - 1;
        ccount[ii] += 1;
        ccount[jj] += 1;
      }

    }
    else { /* if symmetric == true */

      for (face_id = 0; face_id < n_faces; face_id++) {
        ii = *face_cel_p++ - 1;
        jj = *face_cel_p++ - 1;
        if (ii < jj)
          ccount[ii] += 1;
        else
          ccount[jj] += 1;
      }

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

  BFT_MALLOC(ms->col_id, (ms->row_index[ms->n_rows]), cs_int_t);
  ms->col_id = ms->col_id;

  if (have_diag == true) {
    for (ii = 0; ii < ms->n_rows; ii++) {    /* diagonal terms */
      ms->col_id[ms->row_index[ii]] = ii;
    }
  }

  if (face_cell != NULL) {                   /* non-diagonal terms */

    face_cel_p = face_cell;

    if (symmetric == false) {

      for (face_id = 0; face_id < n_faces; face_id++) {
        ii = *face_cel_p++ - 1;
        jj = *face_cel_p++ - 1;
        if (ii < ms->n_rows) {
          ms->col_id[ms->row_index[ii] + ccount[ii]] = jj;
          ccount[ii] += 1;
        }
        if (jj < ms->n_rows) {
          ms->col_id[ms->row_index[jj] + ccount[jj]] = ii;
          ccount[jj] += 1;
        }
      }

    }
    else { /* if symmetric == true */

      for (face_id = 0; face_id < n_faces; face_id++) {
        ii = *face_cel_p++ - 1;
        jj = *face_cel_p++ - 1;
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

  } /* if (face_cell != NULL) */

  BFT_FREE(ccount);

  /* Sort line elements by column id (for better access patterns) */

  if (n_cols_max > 1) {

    cs_int_t  *order = NULL;
    cs_int_t  *new_col_id = NULL;

    BFT_MALLOC(order, n_cols_max, fvm_lnum_t);
    BFT_MALLOC(new_col_id, n_cols_max, cs_int_t);

    for (ii = 0; ii < ms->n_rows; ii++) {
      cs_int_t *col_id = ms->col_id + ms->row_index[ii];
      cs_int_t n_cols = ms->row_index[ii+1] - ms->row_index[ii];
      cs_int_t col_id_prev = -1;
      _order_local(col_id,
                   order,
                   ms->row_index[ii+1] - ms->row_index[ii]);
      for (jj = 0; jj < n_cols; jj++) {
        new_col_id[jj] = col_id[order[jj]];
        if (new_col_id[jj] == col_id_prev)
          ms->direct_assembly = false;
        col_id_prev = new_col_id[jj];
      }
      memcpy(col_id, new_col_id, n_cols*sizeof(cs_int_t));
    }

    BFT_FREE(new_col_id);
    BFT_FREE(order);

  }

  /* Compact elements if necessary */

  if (ms->direct_assembly == false) {

    cs_int_t *tmp_row_index = NULL;
    cs_int_t  kk = 0;

    BFT_MALLOC(tmp_row_index, ms->n_rows+1, cs_int_t);
    memcpy(tmp_row_index, ms->row_index, (ms->n_rows+1)*sizeof(cs_int_t));

    kk = 0;

    for (ii = 0; ii < ms->n_rows; ii++) {
      cs_int_t *col_id = ms->col_id + ms->row_index[ii];
      cs_int_t n_cols = ms->row_index[ii+1] - ms->row_index[ii];
      cs_int_t col_id_prev = -1;
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
    BFT_REALLOC(ms->col_id, (ms->row_index[ms->n_rows]), cs_int_t);

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
_destroy_struct_csr(cs_matrix_struct_csr_t **matrix)
{
  if (matrix != NULL && *matrix !=NULL) {

    cs_matrix_struct_csr_t  *ms = *matrix;

    if (ms->row_index != NULL)
      BFT_FREE(ms->row_index);

    if (ms->col_id != NULL)
      BFT_FREE(ms->col_id);

    if (ms->diag_index != NULL)
      BFT_FREE(ms->diag_index);

    BFT_FREE(ms);

    *matrix = ms;

  }
}

#if 0
/*----------------------------------------------------------------------------
 * Add a diagonal index to a CSR matrix structure, if applicable.
 *
 * parameters:
 *   ms      <->  Matrix structure
 *----------------------------------------------------------------------------*/

static void
_add_diag_index_struct_csr(cs_matrix_struct_csr_t  *ms)
{
  cs_int_t ii;

  /* No need to add a diagonal index for a symmetric matrix, as
     columns are sorted and thus diag_index[ii] == row_index[ii] */

  if (ms->symmetric == true || ms->diag_index != NULL)
    return;

  BFT_MALLOC(ms->diag_index, ms->n_rows, cs_int_t);

  for (ii = 0; ii < ms->n_rows; ii++) {
    cs_int_t kk;
    for (kk = ms->row_index[ii]; kk < ms->row_index[ii+1]; kk++);
    if (ms->col_id[kk] == ii) {
      ms->diag_index[ii] = kk;
      break;
    }
    if (kk == ms->row_index[ii+1]) { /* if some rows have no diagonal value */
      BFT_FREE(ms->diag_index);
      return;
    }
  }

}
#endif

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

  mc->n_prefetch_rows = _cs_glob_matrix_prefetch_rows;

  mc->val = NULL;

  mc->x_prefetch = NULL;

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

    BFT_FREE(*coeff);

  }
}

/*----------------------------------------------------------------------------
 * Set CSR extradiagonal matrix coefficients for the case where direct
 * assignment is possible (i.e. when there are no multiple contributions
 * to a given coefficient).
 *
 * parameters:
 *   matrix    --> Pointer to matrix structure
 *   symmetric --> Indicates if extradiagonal values are symmetric
 *   xa        --> Extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_xa_coeffs_csr_direct(cs_matrix_t      *matrix,
                          cs_bool_t         symmetric,
                          const cs_real_t  *restrict xa)
{
  cs_int_t  ii, jj, face_id;
  cs_matrix_coeff_csr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  /* Copy extra-diagonal values */

  assert(matrix->face_cell != NULL);

  if (symmetric == false) {

    const cs_int_t n_faces = matrix->n_faces;
    const cs_int_t *restrict face_cel_p = matrix->face_cell;

    const cs_real_t  *restrict xa1 = xa;
    const cs_real_t  *restrict xa2 = xa + matrix->n_faces;

    assert(ms->symmetric == false);

    for (face_id = 0; face_id < n_faces; face_id++) {
      cs_int_t kk, ll;
      ii = *face_cel_p++ - 1;
      jj = *face_cel_p++ - 1;
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
  else { /* if symmetric == true */

    if (ms->symmetric == true) {

      const cs_int_t n_faces = matrix->n_faces;
      const cs_int_t *restrict face_cel_p = matrix->face_cell;

      for (face_id = 0; face_id < n_faces; face_id++) {
        cs_int_t kk;
        ii = *face_cel_p++ - 1;
        jj = *face_cel_p++ - 1;
        if (ii < jj && ii < ms->n_rows) {
          for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
          mc->val[kk] = xa[face_id];
        }
        else if (ii > jj && jj < ms->n_rows) {
          for (kk = ms->row_index[jj]; ms->col_id[kk] != ii; kk++);
          mc->val[kk] = xa[face_id];
        }
      }

    }
    else { /* if (ms->symmetric == false) */

      const cs_int_t n_faces = matrix->n_faces;
      const cs_int_t *restrict face_cel_p = matrix->face_cell;

      for (face_id = 0; face_id < n_faces; face_id++) {
        cs_int_t kk, ll;
        ii = *face_cel_p++ - 1;
        jj = *face_cel_p++ - 1;
        if (ii < ms->n_rows) {
          for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
          mc->val[kk] = xa[face_id];
        }
        if (jj < ms->n_rows) {
          for (ll = ms->row_index[jj]; ms->col_id[ll] != ii; ll++);
          mc->val[ll] = xa[face_id];
        }

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
 *   matrix    --> Pointer to matrix structure
 *   symmetric --> Indicates if extradiagonal values are symmetric
 *   xa        --> Extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_xa_coeffs_csr_increment(cs_matrix_t      *matrix,
  cs_bool_t         symmetric,
  const cs_real_t  *restrict xa)
{
  cs_int_t  ii, jj, face_id;
  cs_matrix_coeff_csr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  /* Copy extra-diagonal values */

  assert(matrix->face_cell != NULL);

  if (symmetric == false) {

    const cs_int_t n_faces = matrix->n_faces;
    const cs_int_t *restrict face_cel_p = matrix->face_cell;

    const cs_real_t  *restrict xa1 = xa;
    const cs_real_t  *restrict xa2 = xa + matrix->n_faces;

    assert(ms->symmetric == false);

    for (face_id = 0; face_id < n_faces; face_id++) {
      cs_int_t kk, ll;
      ii = *face_cel_p++ - 1;
      jj = *face_cel_p++ - 1;
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
  else { /* if symmetric == true */

    if (ms->symmetric == true) {

      const cs_int_t n_faces = matrix->n_faces;
      const cs_int_t *restrict face_cel_p = matrix->face_cell;

      for (face_id = 0; face_id < n_faces; face_id++) {
        cs_int_t kk;
        ii = *face_cel_p++ - 1;
        jj = *face_cel_p++ - 1;
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
    else { /* if (ms->symmetric == false) */

      const cs_int_t n_faces = matrix->n_faces;
      const cs_int_t *restrict face_cel_p = matrix->face_cell;

      for (face_id = 0; face_id < n_faces; face_id++) {
        cs_int_t kk, ll;
        ii = *face_cel_p++ - 1;
        jj = *face_cel_p++ - 1;
        if (ii < ms->n_rows) {
          for (kk = ms->row_index[ii]; ms->col_id[kk] != jj; kk++);
          mc->val[kk] += xa[face_id];
        }
        if (jj < ms->n_rows) {
          for (ll = ms->row_index[jj]; ms->col_id[ll] != ii; ll++);
          mc->val[ll] += xa[face_id];
        }

      }

    }

  } /* end of condition on coefficients symmetry */

}

/*----------------------------------------------------------------------------
 * Set CSR matrix coefficients.
 *
 * parameters:
 *   matrix    --> Pointer to matrix structure
 *   symmetric --> Indicates if extradiagonal values are symmetric
 *   da        --> Diagonal values (NULL if all zero)
 *   xa        --> Extradiagonal values (NULL if all zero)
 *----------------------------------------------------------------------------*/

static void
_set_coeffs_csr(cs_matrix_t      *matrix,
                cs_bool_t         symmetric,
                const cs_real_t  *restrict da,
                const cs_real_t  *restrict xa)
{
  cs_int_t  ii, jj;
  cs_matrix_coeff_csr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  if (mc->val == NULL)
    BFT_MALLOC(mc->val, ms->row_index[ms->n_rows], cs_real_t);

  /* Initialize coefficients to zero if assembly is incremental */

  if (ms->direct_assembly == false) {
    cs_int_t val_size = ms->row_index[ms->n_rows];
    for (ii = 0; ii < val_size; ii++)
      mc->val[ii] = 0.0;
  }

  /* Allocate prefetch buffer */

  if (mc->n_prefetch_rows > 0 && mc->x_prefetch == NULL) {
    size_t prefetch_size = ms->n_cols_max * mc->n_prefetch_rows;
    size_t matrix_size = matrix->n_cells + (2 * matrix->n_faces);
    if (matrix_size > prefetch_size)
      prefetch_size = matrix_size;
    BFT_MALLOC(mc->x_prefetch, prefetch_size, cs_real_t);
  }

  /* Copy diagonal values */

  if (matrix->have_diag == true) {

    if (ms->symmetric == false && ms->diag_index == NULL) {

      if (da != NULL) {
        for (ii = 0; ii < ms->n_rows; ii++) {
          cs_int_t kk;
          for (kk = ms->row_index[ii]; ms->col_id[kk] != ii; kk++);
          mc->val[kk] = da[ii];
        }
      }
      else {
        for (ii = 0; ii < ms->n_rows; ii++) {
          cs_int_t kk;
          for (kk = ms->row_index[ii]; ms->col_id[kk] != ii; kk++);
          mc->val[kk] = 0.0;
        }
      }

    }
    else { /* If diagonal index is available, direct assignment */

      const cs_int_t *_diag_index
        = (ms->symmetric == true) ? ms->row_index : ms->diag_index;

      if (da != NULL) {
        for (ii = 0; ii < ms->n_rows; ii++)
          mc->val[_diag_index[ii]] = da[ii];
      }
      else {
        for (ii = 0; ii < ms->n_rows; ii++)
          mc->val[_diag_index[ii]] = 0.0;
      }

    }

  }

  /* Copy extra-diagonal values */

  if (matrix->face_cell != NULL) {

    if (xa != NULL) {

      if (symmetric == false && ms->symmetric == true)
        bft_error(__FILE__, __LINE__, 0,
                  _("Assigning non-symmetric matrix coefficients to a matrix\n"
                    "in a symmetric CSR format."));

      if (ms->direct_assembly == true)
        _set_xa_coeffs_csr_direct(matrix, symmetric, xa);
      else
        _set_xa_coeffs_csr_increment(matrix, symmetric, xa);

    }
    else { /* if (xa == NULL) */

      for (ii = 0; ii < ms->n_rows; ii++) {
        const cs_int_t  *restrict col_id = ms->col_id + ms->row_index[ii];
        cs_real_t  *m_row = mc->val + ms->row_index[ii];
        cs_int_t  n_cols = ms->row_index[ii+1] - ms->row_index[ii];

        for (jj = 0; jj < n_cols; jj++) {
          if (col_id[jj] != ii)
            m_row[jj] = 0.0;
        }

      }

    }

  } /* (matrix->face_cell != NULL) */

}

/*----------------------------------------------------------------------------
 * Release CSR matrix coefficients.
 *
 * parameters:
 *   matrix --> Pointer to matrix structure
 *----------------------------------------------------------------------------*/

static void
_release_coeffs_csr(cs_matrix_t  *matrix)
{
  cs_matrix_coeff_csr_t  *mc = matrix->coeffs;

  if (mc !=NULL) {
    if (mc->val != NULL)
      BFT_FREE(mc->val);
  }

}

/*----------------------------------------------------------------------------
 * Get diagonal of CSR matrix.
 *
 * parameters:
 *   matrix --> Pointer to matrix structure
 *   da     <-- Diagonal (pre-allocated, size: n_rows)
 *----------------------------------------------------------------------------*/

static void
_get_diagonal_csr(const cs_matrix_t  *matrix,
                  cs_real_t          *restrict da)
{
  cs_int_t  ii, jj;
  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_csr_t  *mc = matrix->coeffs;
  cs_int_t  n_rows = ms->n_rows;

  if (matrix->have_diag == true) {

    if (ms->symmetric == false && ms->diag_index == NULL) {

      /* Full rows for non-symmetric structure,
         upper triangular + diagonal part in case of symmetric structure */

      for (ii = 0; ii < n_rows; ii++) {

        const cs_int_t  *restrict col_id = ms->col_id + ms->row_index[ii];
        const cs_real_t  *restrict m_row = mc->val + ms->row_index[ii];
        cs_int_t  n_cols = ms->row_index[ii+1] - ms->row_index[ii];

        for (jj = 0; jj < n_cols; jj++) {
          if (col_id[jj] == ii) {
            da[ii] = m_row[jj];
            break;
          }
          da[ii] = 0.0;
        }

      }

    }

    else { /* if (ms->symmetric == true || ms->diag_index == true) */

      /* If structure is symmetric, diagonal values appear first,
         sor diag_index == row_index */

      const cs_int_t *diag_index
        = (ms->symmetric == true) ? ms->row_index : ms->diag_index;

      for (ii = 0; ii < n_rows; ii++)
        da[ii] = mc->val[diag_index[ii]];

    }

  }
  else { /* if (have_diag == false) */

    for (ii = 0; ii < n_rows; da[ii++] = 0.0);

  }

}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with CSR matrix.
 *
 * parameters:
 *   matrix --> Pointer to matrix structure
 *   x      --> Multipliying vector values
 *   y      <-- Resulting vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_p_l_csr(const cs_matrix_t  *matrix,
                 const cs_real_t    *restrict x,
                 cs_real_t          *restrict y)
{
  cs_int_t  ii, jj;
  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_csr_t  *mc = matrix->coeffs;
  cs_int_t  n_rows = ms->n_rows;

  assert(ms->symmetric == false);

  /* Full rows for non-symmetric structure */

  for (ii = 0; ii < n_rows; ii++) {

    const cs_int_t  *restrict col_id = ms->col_id + ms->row_index[ii];
    const cs_real_t  *restrict m_row = mc->val + ms->row_index[ii];
    cs_int_t  n_cols = ms->row_index[ii+1] - ms->row_index[ii];
    cs_real_t  sii = 0.0;

    /* Tell IBM compiler not to alias */
#if defined(__xlc__)
#pragma disjoint(*x, *y, *m_row, *col_id)
#endif

    for (jj = 0; jj < n_cols; jj++)
      sii += (m_row[jj]*x[col_id[jj]]);

    y[ii] = sii;

  }

}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with symmetric CSR matrix.
 *
 * parameters:
 *   matrix --> Pointer to matrix structure
 *   x      --> Multipliying vector values
 *   y      <-- Resulting vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_p_l_csr_sym(const cs_matrix_t   *matrix,
                     const cs_real_t     *restrict x,
                     cs_real_t           *restrict y)
{
  cs_int_t  ii, jj, n_cols;
  cs_int_t  *restrict col_id;
  cs_real_t  *restrict m_row;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_csr_t  *mc = matrix->coeffs;
  cs_int_t  n_rows = ms->n_rows;

  cs_int_t sym_jj_start = 0;

    /* Tell IBM compiler not to alias */
#if defined(__xlc__)
#pragma disjoint(*x, *y, *m_row, *col_id)
#endif

  assert(ms->symmetric == true);

  /* By construction, the matrix has either a full or an empty
     diagonal structure, so testing this on the first row is enough */

  for (ii = ms->row_index[0]; ii < ms->row_index[1]; ii++) {
    if (ms->col_id[ii] == 0)
      sym_jj_start = 1;
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

    for (jj = 0; jj < n_cols; jj++)
      sii += (m_row[jj]*x[col_id[jj]]);

    y[ii] += sii;

    for (jj = sym_jj_start; jj < n_cols; jj++)
      y[col_id[jj]] += (m_row[jj]*x[ii]);

  }

}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = A.x with CSR matrix (prefetch).
 *
 * parameters:
 *   matrix --> Pointer to matrix structure
 *   x      --> Multipliying vector values
 *   y      <-- Resulting vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_p_l_csr_pf(const cs_matrix_t  *matrix,
                    const cs_real_t    *restrict x,
                    cs_real_t          *restrict y)
{
  cs_int_t  start_row, ii, jj, n_cols;
  cs_int_t  *restrict col_id;
  cs_real_t  *restrict m_row;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_csr_t  *mc = matrix->coeffs;
  cs_int_t  n_rows = ms->n_rows;

  /* Outer loop on prefetch lines */

  for (start_row = 0; start_row < n_rows; start_row += mc->n_prefetch_rows) {

    cs_int_t end_row = start_row + mc->n_prefetch_rows;

    cs_real_t  *restrict prefetch_p = mc->x_prefetch;

      /* Tell IBM compiler not to alias */
#if defined(__xlc__)
#pragma disjoint(*prefetch_p, *y, *m_row)
#pragma disjoint(*prefetch_p, *x, *col_id)
#endif

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
 * Local matrix.vector product y = alpha.A.x + beta.y with CSR matrix.
 *
 * parameters:
 *   alpha  --> Scalar, alpha in alpha.A.x + beta.y
 *   beta   --> Scalar, beta in alpha.A.x + beta.y
 *   matrix --> Pointer to matrix structure
 *   x      --> Multipliying vector values
 *   y      <-> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_alpha_a_x_p_beta_y_csr(cs_real_t           alpha,
                        cs_real_t           beta,
                        const cs_matrix_t  *matrix,
                        const cs_real_t    *restrict x,
                        cs_real_t          *restrict y)
{
  cs_int_t  ii, jj;
  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_csr_t  *mc = matrix->coeffs;
  cs_int_t  n_rows = ms->n_rows;

  assert(ms->symmetric == false);

  /* Full rows for non-symmetric structure */

  for (ii = 0; ii < n_rows; ii++) {

    const cs_int_t  *restrict col_id = ms->col_id + ms->row_index[ii];
    const cs_real_t  *restrict m_row = mc->val + ms->row_index[ii];
    cs_int_t  n_cols = ms->row_index[ii+1] - ms->row_index[ii];
    cs_real_t  sii = 0.0;

    /* Tell IBM compiler not to alias */
#if defined(__xlc__)
#pragma disjoint(*x, *y, *m_row, *col_id)
#endif

    for (jj = 0; jj < n_cols; jj++)
      sii += (m_row[jj]*x[col_id[jj]]);

    y[ii] = (alpha * sii) + (beta * y[ii]);

  }

}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = alpha.A.x + beta.y
 * with symmetric CSR matrix.
 *
 * parameters:
 *   alpha  --> Scalar, alpha in alpha.A.x + beta.y
 *   beta   --> Scalar, beta in alpha.A.x + beta.y
 *   matrix --> Pointer to matrix structure
 *   x      --> Multipliying vector values
 *   y      <-> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_alpha_a_x_p_beta_y_csr_sym(cs_real_t           alpha,
                            cs_real_t           beta,
                            const cs_matrix_t  *matrix,
                            const cs_real_t    *restrict x,
                            cs_real_t          *restrict y)
{
  cs_int_t  ii, jj, n_cols;
  cs_int_t  *restrict col_id;
  cs_real_t  *restrict m_row;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_csr_t  *mc = matrix->coeffs;
  cs_int_t  n_rows = ms->n_rows;

    /* Tell IBM compiler not to alias */
#if defined(__xlc__)
#pragma disjoint(*x, *y, *m_row, *col_id)
#endif

  assert(ms->symmetric == true);

  for (ii = 0; ii < ms->n_rows; ii++)
    y[ii] *= beta;

  for (ii = ms->n_rows; ii < ms->n_cols; ii++)
    y[ii] = 0.0;

  /* Upper triangular + diagonal part in case of symmetric structure */

  for (ii = 0; ii < n_rows; ii++) {

    cs_real_t  sii = 0.0;

    col_id = ms->col_id + ms->row_index[ii];
    m_row = mc->val + ms->row_index[ii];
    n_cols = ms->row_index[ii+1] - ms->row_index[ii];

    for (jj = 0; jj < n_cols; jj++)
      sii += (m_row[jj]*x[col_id[jj]]);

    y[ii] += alpha * sii;

    for (jj = 1; jj < n_cols; jj++)
      y[col_id[jj]] += alpha * (m_row[jj]*x[ii]);

  }

}

/*----------------------------------------------------------------------------
 * Local matrix.vector product y = alpha.A.x + beta.y
 * with CSR matrix (prefetch).
 *
 * parameters:
 *   alpha  --> Scalar, alpha in alpha.A.x + beta.y
 *   beta   --> Scalar, beta in alpha.A.x + beta.y
 *   matrix --> Pointer to matrix structure
 *   x      --> Multipliying vector values
 *   y      <-> Resulting vector
 *----------------------------------------------------------------------------*/

static void
_alpha_a_x_p_beta_y_csr_pf(cs_real_t           alpha,
                           cs_real_t           beta,
                           const cs_matrix_t  *matrix,
                           const cs_real_t    *restrict x,
                           cs_real_t          *restrict y)
{
  cs_int_t  start_row, ii, jj, n_cols;
  cs_int_t  *restrict col_id;
  cs_real_t  *restrict m_row;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;
  const cs_matrix_coeff_csr_t  *mc = matrix->coeffs;
  cs_int_t  n_rows = ms->n_rows;

  /* Outer loop on prefetch lines */

  for (start_row = 0; start_row < n_rows; start_row += mc->n_prefetch_rows) {

    cs_int_t end_row = start_row + mc->n_prefetch_rows;
    cs_real_t  *restrict prefetch_p = mc->x_prefetch;

      /* Tell IBM compiler not to alias */
#if defined(__xlc__)
#pragma disjoint(*prefetch_p, *x, *col_id)
#pragma disjoint(*prefetch_p, *y, *m_row, *col_id)
#endif

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

      y[ii] = (alpha * sii) + (beta * y[ii]);

    }

  }

}

/*============================================================================
 *  Public function definitions for Fortran API
 *============================================================================*/

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
 *   type        --> Type of matrix considered
 *   symmetric   --> Indicates if a symmetric variant of the matrix type
 *                   should be used
 *   have_diag   --> Indicates if the diagonal structure contains nonzeroes
 *   periodic    --> Indicates if periodicity is present
 *   n_cells     --> Local number of cells
 *   n_cells_ext --> Local number of cells + ghost cells sharing a face
 *   n_faces     --> Local number of internal faces
 *   cell_num    --> Global cell numbers (1 to n)
 *   face_cell   --> Face -> cells connectivity (1 to n)
 *   halo        --> Halo structure associated with cells, or NULL
 *
 * returns:
 *   pointer to created matrix structure;
 *----------------------------------------------------------------------------*/

cs_matrix_t *
cs_matrix_create(cs_matrix_type_t   type,
                 cs_bool_t          symmetric,
                 cs_bool_t          have_diag,
                 cs_bool_t          periodic,
                 cs_int_t           n_cells,
                 cs_int_t           n_cells_ext,
                 cs_int_t           n_faces,
                 const fvm_gnum_t  *cell_num,
                 const cs_int_t    *face_cell,
                 const cs_halo_t   *halo)
{
  cs_matrix_t *m;

  BFT_MALLOC(m, 1, cs_matrix_t);

  m->type = type;
  m->periodic = periodic;
  m->have_diag = have_diag;

  m->n_cells = n_cells;
  m->n_cells_ext = n_cells_ext;
  m->n_faces = n_faces;

  /* Define Structure */

  switch(m->type) {
  case CS_MATRIX_NATIVE:
    m->_structure = _create_struct_native(n_cells,
                                          n_cells_ext,
                                          n_faces,
                                          face_cell);
    m->coeffs = _create_coeff_native();
    break;
  case CS_MATRIX_CSR:
    m->_structure = _create_struct_csr(symmetric,
                                       have_diag,
                                       n_cells,
                                       n_cells_ext,
                                       n_faces,
                                       face_cell);
    m->coeffs = _create_coeff_csr();
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Handling of matrixes in %s format\n"
                "is not operational yet."),
              _(cs_matrix_type_name[type]));
    break;
  }

  m->structure = m->_structure;

  /* Set pointers to structures shared from mesh here */

  m->face_cell = face_cell;
  m->cell_num = cell_num;
  m->halo = halo;

  /* Set function pointers here */

  switch(m->type) {

  case CS_MATRIX_NATIVE:
    m->set_coefficients = _set_coeffs_native;
    m->release_coefficients = _release_coeffs_native;
    m->get_diagonal = _get_diagonal_native;
    m->vector_multiply = _mat_vec_p_l_native;
    m->alpha_a_x_p_beta_y = _alpha_a_x_p_beta_y_native;
    break;

  case CS_MATRIX_CSR:
    m->set_coefficients = _set_coeffs_csr;
    m->release_coefficients = _release_coeffs_csr;
    m->get_diagonal = _get_diagonal_csr;
    if (symmetric == false) {
      if (_cs_glob_matrix_prefetch_rows > 0) {
        m->vector_multiply = _mat_vec_p_l_csr_pf;
        m->alpha_a_x_p_beta_y = _alpha_a_x_p_beta_y_csr_pf;
      }
      else {
        m->vector_multiply = _mat_vec_p_l_csr;
        m->alpha_a_x_p_beta_y = _alpha_a_x_p_beta_y_csr;
      }
    }
    else { /* if (symmetric == true) */
      m->vector_multiply = _mat_vec_p_l_csr_sym;
      m->alpha_a_x_p_beta_y = _alpha_a_x_p_beta_y_csr_sym;
    }
    break;

  default:
    assert(0);
    m->set_coefficients = NULL;
    m->vector_multiply = NULL;
    m->alpha_a_x_p_beta_y = NULL;

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
        cs_matrix_struct_native_t *_structure = m->_structure;
        cs_matrix_coeff_native_t *coeffs = m->coeffs;
        _destroy_struct_native(&_structure);
        _destroy_coeff_native(&coeffs);
        m->structure = NULL; m->coeffs = NULL;
      }
      break;
    case CS_MATRIX_CSR:
      {
        cs_matrix_struct_csr_t *_structure = m->_structure;
        cs_matrix_coeff_csr_t *coeffs = m->coeffs;
        _destroy_struct_csr(&_structure);
        _destroy_coeff_csr(&coeffs);
        m->structure = NULL; m->coeffs = NULL;
      }
      break;
    default:
      break;
    }

    /* Now free main structure */

    BFT_FREE(*matrix);
  }
}

/*----------------------------------------------------------------------------
 * Return number of columns in matrix.
 *
 * parameters:
 *   matrix --> Pointer to matrix structure
 *----------------------------------------------------------------------------*/

cs_int_t
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
 *   matrix --> Pointer to matrix structure
 *----------------------------------------------------------------------------*/

cs_int_t
cs_matrix_get_n_rows(const cs_matrix_t  *matrix)
{
  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("The matrix is not defined."));

  return matrix->n_cells;
}

/*----------------------------------------------------------------------------
 * Set matrix coefficients.
 *
 * Depending on current options and initialization, values will be copied
 * or simply mapped.
 *
 * parameters:
 *   matrix    <-> Pointer to matrix structure
 *   symmetric --> Indicates if matrix coefficients are symmetric
 *   da        --> Diagonal values (NULL if zero)
 *   xa        --> Extradiagonal values (NULL if zero)
 *----------------------------------------------------------------------------*/

void
cs_matrix_set_coefficients(cs_matrix_t      *matrix,
                           cs_bool_t         symmetric,
                           const cs_real_t  *da,
                           const cs_real_t  *xa)
{
  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("The matrix is not defined."));

  if (matrix->set_coefficients != NULL)
    matrix->set_coefficients(matrix, symmetric, da, xa);
}

/*----------------------------------------------------------------------------
 * Release matrix coefficients.
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
}

/*----------------------------------------------------------------------------
 * Get matrix diagonal values.
 *
 * parameters:
 *   matrix --> Pointer to matrix structure
 *   da     <-- Diagonal (pre-allocated, size: n_cells)
 *----------------------------------------------------------------------------*/

void
cs_matrix_get_diagonal(const cs_matrix_t  *matrix,
                       cs_real_t          *restrict da)
{
  /* Check API state */

  if (matrix == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("The matrix is not defined."));

  if (matrix->get_diagonal != NULL)
    matrix->get_diagonal(matrix, da);
}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x
 *
 * This function includes a halo update of x prior to multiplication by A.
 *
 * parameters:
 *   rotation_mode --> Halo update option for rotational periodicity
 *   matrix        --> Pointer to matrix structure
 *   x             <-> Multipliying vector values (ghost values updated)
 *   y             <-- Resulting vector
 *----------------------------------------------------------------------------*/

void
cs_matrix_vector_multiply(cs_perio_rota_t     rotation_mode,
                          const cs_matrix_t  *matrix,
                          cs_real_t          *restrict x,
                          cs_real_t          *restrict y)
{
  size_t ii;
  size_t n_cells_ext = matrix->n_cells_ext;

  /* Synchronize for parallelism and periodicity first */

  for (ii = matrix->n_cells; ii < n_cells_ext; y[ii++] = 0.);

  /* Update distant ghost cells */

  if (matrix->halo != NULL)
    cs_halo_sync_var(matrix->halo, CS_HALO_STANDARD, x);

  /* Synchronize periodic values */

  if (matrix->periodic) {
    if (rotation_mode == CS_PERIO_ROTA_IGNORE)
      bft_error(__FILE__, __LINE__, 0, _cs_glob_perio_ignore_error_str);
    cs_perio_sync_var_scal(matrix->halo, CS_HALO_STANDARD, rotation_mode, x);
  }

  /* Now call local matrix.vector product */

  if (matrix->vector_multiply != NULL)
    matrix->vector_multiply(matrix, x, y);

  else if (matrix->alpha_a_x_p_beta_y != NULL)
    matrix->alpha_a_x_p_beta_y(1.0, 0.0, matrix, x, y);
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
 *   matrix    --> Pointer to matrix structure
 *   x         --> Multipliying vector values
 *   y         <-- Resulting vector
 *----------------------------------------------------------------------------*/

void
cs_matrix_vector_multiply_nosync(const cs_matrix_t  *matrix,
                                 const cs_real_t    *x,
                                 cs_real_t          *restrict y)
{
  if (matrix != NULL) {

    if (matrix->vector_multiply != NULL)
      matrix->vector_multiply(matrix, x, y);

    else if (matrix->alpha_a_x_p_beta_y != NULL)
      matrix->alpha_a_x_p_beta_y(1.0, 0.0, matrix, x, y);

  }
}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = alpha.A.x + beta.y
 *
 * This function includes a halo update of x prior to multiplication by A.
 *
 * parameters:
 *   rotation_mode --> Halo update option for rotational periodicity
 *   alpha         --> Scalar, alpha in alpha.A.x + beta.y
 *   beta          --> Scalar, beta in alpha.A.x + beta.y
 *   matrix        --> Pointer to matrix structure
 *   x             <-> Multipliying vector values (ghost values updated)
 *   y             <-- Resulting vector
 *----------------------------------------------------------------------------*/

void
cs_matrix_alpha_a_x_p_beta_y(cs_perio_rota_t     rotation_mode,
                             cs_real_t           alpha,
                             cs_real_t           beta,
                             const cs_matrix_t  *matrix,
                             cs_real_t          *restrict x,
                             cs_real_t          *restrict y)
{
  /* Update distant ghost cells */

  if (matrix->halo != NULL)
    cs_halo_sync_var(matrix->halo, CS_HALO_STANDARD, x);

  /* Synchronize periodic values */

  if (matrix->periodic) {
    if (rotation_mode == CS_PERIO_ROTA_IGNORE)
      bft_error(__FILE__, __LINE__, 0, _cs_glob_perio_ignore_error_str);
    cs_perio_sync_var_scal(matrix->halo, CS_HALO_STANDARD, rotation_mode, x);
  }

  /* Now call local matrix.vector product */

  if (matrix->alpha_a_x_p_beta_y != NULL)
    matrix->alpha_a_x_p_beta_y(alpha, beta, matrix, x, y);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
