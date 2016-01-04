#ifndef __CS_MATRIX_PRIV_H__
#define __CS_MATRIX_PRIV_H__

/*============================================================================
 * Private types for sparse matrix representation and operations.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "cs_matrix.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Formats currently supported:
 *
 *  - Native
 *  - Compressed Sparse Row (CSR)
 *  - Symmetric Compressed Sparse Row (CSR_SYM)
 */

/*----------------------------------------------------------------------------
 * Function pointer types
 *----------------------------------------------------------------------------*/

typedef void
(cs_matrix_set_coeffs_t) (cs_matrix_t      *matrix,
                          bool              symmetric,
                          bool              copy,
                          const cs_real_t  *restrict da,
                          const cs_real_t  *restrict xa);

typedef void
(cs_matrix_release_coeffs_t) (cs_matrix_t  *matrix);

typedef void
(cs_matrix_copy_diagonal_t) (const cs_matrix_t  *matrix,
                             cs_real_t          *restrict da);

typedef void
(cs_matrix_vector_product_t) (bool                exclude_diag,
                              const cs_matrix_t  *matrix,
                              const cs_real_t    *restrict x,
                              cs_real_t          *restrict y);

/*----------------------------------------------------------------------------
 * Matrix types
 *----------------------------------------------------------------------------*/

/* Native matrix structure representation */
/*----------------------------------------*/

/* Note: the members of this structure are already available through the top
 *       matrix structure, but are replicated here in case of future removal
 *       from the top structure (which would require computation/assignment of
 *       matrix coefficients in another form) */

typedef struct _cs_matrix_struct_native_t {

  cs_lnum_t          n_cells;       /* Local number of cells */
  cs_lnum_t          n_cells_ext;   /* Local number of participating cells
                                       (cells + ghost cells sharing a face) */
  cs_lnum_t          n_faces;       /* Local number of faces
                                       (for extra-diagonal terms */

  /* Pointers to shared arrays */

  const cs_lnum_2_t  *face_cell;    /* Face -> cells connectivity (0 to n-1) */

} cs_matrix_struct_native_t;

/* Native matrix coefficients */
/*----------------------------*/

typedef struct _cs_matrix_coeff_native_t {

  bool              symmetric;       /* Symmetry indicator */
  int               max_db_size;     /* Current max allocated diag block size */
  int               max_eb_size;     /* Current max allocated extradiag block size */

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

  cs_lnum_t         n_rows;           /* Local number of rows */
  cs_lnum_t         n_cols;           /* Local number of columns
                                         (> n_rows in case of ghost cells) */
  cs_lnum_t         n_cols_max;       /* Maximum number of nonzero values
                                         on a given row */

  /* Pointers to structure arrays and info (row_index, col_id) */

  bool              have_diag;        /* Has non-zero diagonal */
  bool              direct_assembly;  /* True if each value corresponds to
                                         a unique face ; false if multiple
                                         faces contribute to the same
                                         value (i.e. we have split faces) */

  cs_lnum_t        *row_index;        /* Row index (0 to n-1) */
  cs_lnum_t        *col_id;           /* Column id (0 to n-1) */

} cs_matrix_struct_csr_t;

/* CSR matrix coefficients representation */
/*----------------------------------------*/

typedef struct _cs_matrix_coeff_csr_t {

  int               n_prefetch_rows;  /* Number of rows at a time for which
                                         the x values in y = Ax should be
                                         prefetched (0 if no prefetch) */

  cs_real_t        *val;              /* Matrix coefficients */

  cs_real_t        *x_prefetch;       /* Prefetch array for x in y = Ax */

  /* Pointers to auxiliary arrays used for queries */

  const cs_real_t  *d_val;            /* Pointer to diagonal matrix
                                         coefficients, if queried */
  cs_real_t        *_d_val;           /* Diagonal matrix coefficients,
                                         if queried */

} cs_matrix_coeff_csr_t;

/* CSR_SYM (Symmetric Compressed Sparse Row) matrix structure representation */
/*---------------------------------------------------------------------------*/

typedef struct _cs_matrix_struct_csr_sym_t {

  cs_lnum_t         n_rows;           /* Local number of rows */
  cs_lnum_t         n_cols;           /* Local number of columns
                                         (> n_rows in case of ghost cells) */
  cs_lnum_t         n_cols_max;       /* Maximum number of nonzero values
                                         on a given row */

  /* Pointers to structure arrays and info (row_index, col_id) */

  bool              have_diag;        /* Has non-zero diagonal */
  bool              direct_assembly;  /* True if each value corresponds to
                                         a unique face ; false if multiple
                                         faces contribute to the same
                                         value (i.e. we have split faces) */

  cs_lnum_t        *row_index;        /* Row index (0 to n-1) */
  cs_lnum_t        *col_id;           /* Column id (0 to n-1) */

} cs_matrix_struct_csr_sym_t;

/* symmetric CSR matrix coefficients representation */
/*--------------------------------------------------*/

typedef struct _cs_matrix_coeff_csr_sym_t {

  cs_real_t        *val;              /* Matrix coefficients */

  /* Pointers to auxiliary arrays used for queries */

  const cs_real_t  *d_val;            /* Pointer to diagonal matrix
                                         coefficients, if queried */
  cs_real_t        *_d_val;           /* Diagonal matrix coefficients,
                                         if queried */

} cs_matrix_coeff_csr_sym_t;

/* MSR matrix coefficients representation */
/*----------------------------------------*/

typedef struct _cs_matrix_coeff_msr_t {

  int              n_prefetch_rows;   /* Number of rows at a time for which
                                         the x values in y = Ax should be
                                         prefetched (0 if no prefetch) */
  int              max_db_size;       /* Current max allocated block size */
  int              max_eb_size;       /* Current max allocated extradiag block size */

  /* Pointers to possibly shared arrays */

  const cs_real_t  *d_val;            /* Diagonal matrix coefficients */

  /* Pointers to private arrays (NULL if shared) */

  cs_real_t        *_d_val;           /* Diagonal matrix coefficients */
  cs_real_t        *x_val;            /* Extra-diagonal matrix coefficients */

  cs_real_t        *x_prefetch;       /* Prefetch array for x in y = Ax */

} cs_matrix_coeff_msr_t;

/* Matrix structure (representation-independent part) */
/*----------------------------------------------------*/

struct _cs_matrix_structure_t {

  cs_matrix_type_t       type;         /* Matrix storage and definition type */

  cs_lnum_t              n_cells;      /* Local number of cells */
  cs_lnum_t              n_cells_ext;  /* Local number of participating cells
                                          (cells + ghost cells sharing a face) */
  cs_lnum_t              n_faces;      /* Local Number of mesh faces
                                          (necessary to affect coefficients) */

  void                  *structure;    /* Matrix structure */

  /* Pointers to shared arrays from mesh structure
     (face->cell connectivity for coefficient assignment,
     local->local cell numbering for future info or renumbering,
     and halo) */

  const cs_lnum_2_t     *face_cell;    /* Face -> cells connectivity */
  const cs_gnum_t       *cell_num;     /* Global cell numbers */
  const cs_halo_t       *halo;         /* Parallel or periodic halo */
  const cs_numbering_t  *numbering;    /* Vectorization or thread-related
                                          numbering information */
};

/* Structure associated with Matrix (representation-independent part) */
/*--------------------------------------------------------------------*/

struct _cs_matrix_t {

  cs_matrix_type_t       type;         /* Matrix storage and definition type */

  cs_lnum_t              n_cells;      /* Local number of cells */
  cs_lnum_t              n_cells_ext;  /* Local number of participating cells
                                          (cells + ghost cells sharing a face) */
  cs_lnum_t              n_faces;      /* Local Number of mesh faces
                                          (necessary to affect coefficients) */

  cs_matrix_fill_type_t  fill_type;    /* Matrix fill type */

  int                    db_size[4];   /* Diag Block size, including padding:
                                          0: useful block size
                                          1: vector block extents
                                          2: matrix line extents
                                          3: matrix line*column extents */

  int                    eb_size[4];   /* Extradiag block size, including padding:
                                          0: useful block size
                                          1: vector block extents
                                          2: matrix line extents
                                          3: matrix line*column extents */

  /* Pointer to shared structure */

  const void            *structure;    /* Matrix structure */

  /* Pointers to shared arrays from mesh structure
     (face->cell connectivity for coefficient assignment,
     local->local cell numbering for future info or renumbering,
     and halo) */

  const cs_lnum_2_t     *face_cell;    /* Face -> cells connectivity */
  const cs_gnum_t       *cell_num;     /* Global cell numbers */
  const cs_halo_t       *halo;         /* Parallel or periodic halo */
  const cs_numbering_t  *numbering;    /* Vectorization or thread-related
                                          numbering information */

  /* Pointers to shared arrays from coefficient assignment
     these may become NULL in the future if coefficients are passed
     directly from non "native" types, so precautions may need to be
     taken when accessing those */

  const cs_real_t      *xa;            /* Extra-diagonal terms */

  /* Pointer to private data */

  void                  *coeffs;       /* Matrix coefficients */

  /* Function pointers */

  cs_matrix_set_coeffs_t            *set_coefficients;
  cs_matrix_release_coeffs_t        *release_coefficients;
  cs_matrix_copy_diagonal_t         *copy_diagonal;

  /* Function pointer arrays, with CS_MATRIX_N_FILL_TYPES variants:
     fill_type*2 + exclude_diagonal_flag */

  cs_matrix_vector_product_t        *vector_multiply[CS_MATRIX_N_FILL_TYPES][2];

  /* Loop length parameter for some SpMv algorithms */

  int                                loop_length[CS_MATRIX_N_FILL_TYPES][2];

};

/* Structure used for tuning variants */
/*------------------------------------*/

struct _cs_matrix_variant_t {

  char                   name[32];     /* Variant name */

  cs_matrix_type_t       type;         /* Matrix storage and definition type */

  /* Loop length parameter for some SpMv algorithms */

  int                    loop_length[CS_MATRIX_N_FILL_TYPES][2];

  /* Function pointer arrays, with variants:
     fill_type + exclude_diagonal_flag */

  cs_matrix_vector_product_t   *vector_multiply[CS_MATRIX_N_FILL_TYPES][2];

  /* Measured structure creation cost, or -1 otherwise */

  double  matrix_create_cost;

  /* Measured assignment costs for each available fill type, or -1 otherwise */

  double  matrix_assign_cost[CS_MATRIX_N_FILL_TYPES];

  /* Measured operation costs for each available operation, or -1 otherwise
     fill_type*2 + exclude_diagonal_flag */

  double  matrix_vector_cost[CS_MATRIX_N_FILL_TYPES][2];

};

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Semi-private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MATRIX_PRIV_H__ */
