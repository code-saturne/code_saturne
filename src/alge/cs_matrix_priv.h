#ifndef __CS_MATRIX_PRIV_H__
#define __CS_MATRIX_PRIV_H__

/*============================================================================
 * Private types for sparse matrix representation and operations.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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
 *  - Modified Compressed Sparse Row (MSR), with separate diagonal
 *  - Symmetric Compressed Sparse Row (CSR_SYM)
 */

/*----------------------------------------------------------------------------
 * Function pointer types
 *----------------------------------------------------------------------------*/

typedef void
(cs_matrix_set_coeffs_t) (cs_matrix_t        *matrix,
                          bool                symmetric,
                          bool                copy,
                          cs_lnum_t           n_edges,
                          const cs_lnum_2_t  *restrict edges,
                          const cs_real_t    *restrict da,
                          const cs_real_t    *restrict xa);

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

  cs_lnum_t          n_rows;        /* Local number of rows */
  cs_lnum_t          n_cols_ext;    /* Local number of columns + ghosts */
  cs_lnum_t          n_edges;       /* Local number of graph edges
                                       (for extra-diagonal terms) */

  /* Pointers to shared arrays */

  const cs_lnum_2_t  *edges;        /* Edges (symmetric row <-> column)
                                       connectivity */

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
  cs_lnum_t         n_cols_ext;       /* Local number of columns + ghosts */

  /* Pointers to structure arrays and info (row_index, col_id) */

  bool              have_diag;        /* Has non-zero diagonal */
  bool              direct_assembly;  /* True if each value corresponds to
                                         a unique face ; false if multiple
                                         faces contribute to the same
                                         value (i.e. we have split faces) */

  const cs_lnum_t  *row_index;        /* Pointer to row index (0 to n-1) */
  const cs_lnum_t  *col_id;           /* Pointer to column id (0 to n-1) */

  cs_lnum_t        *_row_index;       /* Row index (0 to n-1), if owner */
  cs_lnum_t        *_col_id;          /* Column id (0 to n-1), if owner */

} cs_matrix_struct_csr_t;

/* CSR matrix coefficients representation */
/*----------------------------------------*/

typedef struct _cs_matrix_coeff_csr_t {

  /* Pointers to possibly shared arrays */

  const cs_real_t  *val;              /* Matrix coefficients */

  /* Pointers to private arrays (NULL if shared) */

  cs_real_t        *_val;             /* Diagonal matrix coefficients */

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
                                         (> n_rows in case of ghost columns) */

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

  int              max_db_size;       /* Current max allocated block size */
  int              max_eb_size;       /* Current max allocated extradiag
                                         block size */

  /* Pointers to possibly shared arrays */

  const cs_real_t  *d_val;            /* Diagonal matrix coefficients */
  const cs_real_t  *x_val;            /* Extra-diagonal matrix coefficients */

  /* Pointers to private arrays (NULL if shared) */

  cs_real_t        *_d_val;           /* Diagonal matrix coefficients */
  cs_real_t        *_x_val;           /* Extra-diagonal matrix coefficients */

} cs_matrix_coeff_msr_t;

/* Matrix structure (representation-independent part) */
/*----------------------------------------------------*/

struct _cs_matrix_structure_t {

  cs_matrix_type_t       type;         /* Matrix storage and definition type */

  cs_lnum_t              n_rows;       /* Local number of rows */
  cs_lnum_t              n_cols_ext;   /* Local number of columns + ghosts */

  void                  *structure;    /* Matrix structure */

  /* Pointers to shared arrays from mesh structure
     (face->cell connectivity for coefficient assignment,
     local->local cell numbering for future info or renumbering,
     and halo) */

  const cs_halo_t       *halo;         /* Parallel or periodic halo */
  const cs_numbering_t  *numbering;    /* Vectorization or thread-related
                                          numbering information */

  const cs_matrix_assembler_t  *assembler;   /* Associated matrix assembler */
};

/* Structure associated with Matrix (representation-independent part) */
/*--------------------------------------------------------------------*/

struct _cs_matrix_t {

  cs_matrix_type_t       type;         /* Matrix storage and definition type */

  cs_lnum_t              n_rows;       /* Local number of rows */
  cs_lnum_t              n_cols_ext;   /* Local number of columns + ghosts */

  cs_matrix_fill_type_t  fill_type;    /* Matrix fill type */

  bool                   symmetric;    /* true if coefficients are symmetric */

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

  const void            *structure;    /* Possibly shared matrix structure */
  void                  *_structure;   /* Private matrix structure */

  /* Pointers to arrays possibly shared from mesh structure
     (graph edges: face->cell connectivity for coefficient assignment,
     rows: local->local cell numbering for future info or renumbering,
     and halo) */

  const cs_halo_t       *halo;         /* Parallel or periodic halo */
  const cs_numbering_t  *numbering;    /* Vectorization or thread-related
                                          numbering information */

  const cs_matrix_assembler_t  *assembler;   /* Associated matrix assembler */

  /* Pointer to shared arrays from coefficient assignment from
     "native" type. This should be removed in the future, but requires
     removing the dependency to the native structure in the multigrid
     code first. */

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

};

/* Structure used for tuning variants */
/*------------------------------------*/

struct _cs_matrix_variant_t {

  char                   name[32];     /* Variant name */

  cs_matrix_type_t       type;         /* Matrix storage and definition type */

  /* Function pointer arrays, with variants:
     fill_type + exclude_diagonal_flag */

  cs_matrix_vector_product_t   *vector_multiply[CS_MATRIX_N_FILL_TYPES][2];

  /* Measured structure creation cost, or -1 otherwise */

  double  matrix_create_cost;

  /* Measured assignment costs for each available fill type, or -1 otherwise */

  double  matrix_assign_cost[CS_MATRIX_N_FILL_TYPES];

  /* Measured operation costs for each available operation, or -1 otherwise
     fill_type*2 + exclude_diagonal_flag */

  double  matrix_vector_cost[CS_MATRIX_N_FILL_TYPES][2][2];

};

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Semi-private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create CSR matrix coefficients.
 *
 * returns:
 *   pointer to allocated CSR coefficients structure.
 *----------------------------------------------------------------------------*/

cs_matrix_coeff_csr_t *
cs_matrix_create_coeff_csr(void);

/*----------------------------------------------------------------------------
 * Destroy CSR matrix coefficients.
 *
 * parameters:
 *   coeff  <->  Pointer to CSR matrix coefficients pointer
 *----------------------------------------------------------------------------*/

void
cs_matrix_destroy_coeff_csr(cs_matrix_coeff_csr_t **coeff);

/*----------------------------------------------------------------------------
 * Release shared CSR matrix coefficients.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *----------------------------------------------------------------------------*/

void
cs_matrix_release_coeffs_csr(cs_matrix_t  *matrix);

/*----------------------------------------------------------------------------
 * Copy diagonal of CSR matrix.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *   da     --> Diagonal (pre-allocated, size: n_rows)
 *----------------------------------------------------------------------------*/

void
cs_matrix_copy_diagonal_csr(const cs_matrix_t  *matrix,
                            cs_real_t          *restrict da);

/*----------------------------------------------------------------------------
 * Destroy CSR matrix structure.
 *
 * parameters:
 *   matrix  <->  Pointer to CSR matrix structure pointer
 *----------------------------------------------------------------------------*/

void
cs_matrix_destroy_struct_csr(cs_matrix_struct_csr_t  **matrix);

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
                      cs_real_t          *restrict y);

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
                          cs_real_t          *restrict y);

#endif /* defined (HAVE_MKL) */

/*----------------------------------------------------------------------------
 * Copy diagonal of native or MSR matrix.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *   da     --> Diagonal (pre-allocated, size: n_cols)
 *----------------------------------------------------------------------------*/

void
cs_matrix_copy_diagonal_separate(const cs_matrix_t  *matrix,
                                 cs_real_t          *restrict da);

/*----------------------------------------------------------------------------
 * Create MSR matrix coefficients.
 *
 * returns:
 *   pointer to allocated MSR coefficients structure.
 *----------------------------------------------------------------------------*/

cs_matrix_coeff_msr_t *
cs_matrix_create_coeff_msr(void);

/*----------------------------------------------------------------------------
 * Release shared MSR matrix coefficients.
 *
 * parameters:
 *   matrix <-- Pointer to matrix structure
 *----------------------------------------------------------------------------*/

void
cs_matrix_release_coeffs_msr(cs_matrix_t  *matrix);

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
                      cs_real_t          *restrict y);

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
                          cs_real_t          *restrict y);
#endif /* defined (HAVE_MKL) */

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MATRIX_PRIV_H__ */
