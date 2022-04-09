#ifndef __CS_MATRIX_PRIV_H__
#define __CS_MATRIX_PRIV_H__

/*============================================================================
 * Private types for sparse matrix representation and operations.
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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "cs_matrix.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
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
(cs_matrix_destroy_struct_t) (void  **ms);

typedef void
(cs_matrix_destroy_coeffs_t) (cs_matrix_t  *matrix);

typedef void
(cs_matrix_destroy_adaptor_t) (cs_matrix_t  *matrix);

typedef void
(cs_matrix_copy_diagonal_t) (const cs_matrix_t  *matrix,
                             cs_real_t          *restrict da);

typedef const cs_real_t *
(cs_matrix_get_diagonal_t)(const cs_matrix_t  *matrix);

typedef cs_matrix_assembler_values_t *
(cs_matrix_assembler_values_create_t) (cs_matrix_t  *matrix,
                                       cs_lnum_t     diag_block_size,
                                       cs_lnum_t     extra_diag_block_size);

/*----------------------------------------------------------------------------
 * Function pointer for matrix-veector product (y = A.x).
 *
 * parameters:
 *   matrix       <-> pointer to matrix structure
 *   exclude_diag <-- if true, compute (A-D).x instead of A.x
 *   sync         <-- if true, synchronize ghost values
 *   x            <-- x input vector (may be synchronized by this function)
 *   y            <-- y output vector
 *----------------------------------------------------------------------------*/

typedef void
(cs_matrix_vector_product_t) (cs_matrix_t  *matrix,
                              bool          exclude_diag,
                              bool          sync,
                              cs_real_t    *restrict x,
                              cs_real_t    *restrict y);

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

/* CSR (Compressed Sparse Row) matrix structure representation */
/*-------------------------------------------------------------*/

typedef struct _cs_matrix_struct_csr_t {

  cs_lnum_t         n_rows;           /* Local number of rows */
  cs_lnum_t         n_cols_ext;       /* Local number of columns + ghosts */

  /* Pointers to structure arrays and info (row_index, col_id) */

  bool              have_diag;        /* Has non-zero diagonal */
  bool              direct_assembly;  /* True if each value corresponds to
                                         a unique face; false if multiple
                                         faces contribute to the same
                                         value (i.e. we have split faces) */

  const cs_lnum_t  *row_index;        /* Pointer to row index (0 to n-1) */
  const cs_lnum_t  *col_id;           /* Pointer to column id (0 to n-1) */

  cs_lnum_t        *_row_index;       /* Row index (0 to n-1), if owner */
  cs_lnum_t        *_col_id;          /* Column id (0 to n-1), if owner */

} cs_matrix_struct_csr_t;

/* Distributed matrix structure representation */
/*---------------------------------------------*/

/* This storage assumes a representation in the following form:

   - D+E+H

   Where D is the diagonal, E the local extra-diagonal, and H is the part
   of the matrix referencing halo values (only the upper-part of the
   coefficients are needed on a given rank).

   In cases where a partial block structure is used for the diagonal
   (usually at boundary cells), the diagonal values may be indexed.

   Sub-matrices are stored in CSR format, to allow easy indexing
   or access to a given row's elements (which is useful for multigrid
   coarsening).

   Since the H portion of the matrix may be very sparse, with most rows
   empty, the matching row ids can be stored in addition to the column ids,
   defining a simple coordinate-type structure as a COO-type SpMv product
   can be expected to be cheaper than a CSR-based one in this case.
 */

typedef struct _cs_matrix_struct_dist_t {

  cs_lnum_t         n_rows;           /* Local number of rows */
  cs_lnum_t         n_cols_ext;       /* Local number of columns + ghosts */

  /* Pointers to structure arrays and info */

  cs_matrix_struct_csr_t  e;          /* E (local off-diagonal) CSR structure */
  cs_matrix_struct_csr_t  h;          /* H (halo) CSR structure */

  cs_lnum_t               *h_row_id;  /* Optional row id for coordinates
                                         format (col_id in h structure) */

} cs_matrix_struct_dist_t;

/* CSR matrix coefficients representation */
/*----------------------------------------*/

typedef struct _cs_matrix_coeff_csr_t {

  /* Pointers to possibly shared arrays */

  const cs_real_t  *val;              /* Matrix coefficients */

  /* Pointers to private arrays (NULL if shared) */

  cs_real_t        *_val;             /* Matrix coefficients */

  /* Pointers to auxiliary arrays used for queries */

  const cs_real_t  *d_val;            /* Pointer to diagonal matrix
                                         coefficients, if queried */
  cs_real_t        *_d_val;           /* Diagonal matrix coefficients,
                                         if queried */

} cs_matrix_coeff_csr_t;

/* Distributed matrix coefficients representation */
/*------------------------------------------------*/

/* Used for native, MSR, and distributed matrices */

typedef struct _cs_matrix_coeff_dist_t {

  bool             symmetric;         /* Symmetry indicator */

  int              db_size;           /* Diagonal block size */
  int              eb_size;           /* Extra-diagonal  block size */

  /* Pointers to possibly shared arrays */

  const cs_real_t  *d_val;            /* D (diagonal-only) coefficients */
  const cs_real_t  *e_val;            /* E (extra-diagonal) coefficients */
  const cs_real_t  *h_val;            /* H (halo-only) coefficients */

   /* Pointers to private arrays.
      NULL if shared:
      * If non-NULL, d_val, e_val, and h_val should point
        to matching private array.
      * In the case of CSR storage, where diagonal values can be stored in
        the e_val array, d_val will be NULL but _d_val may be used to store
        (cache) diagonal values. */

  cs_real_t        *_d_val;          /* D (diagonal) coefficients */
  cs_real_t        *_e_val;          /* E (local extra-diagonal) coefficients */
  cs_real_t        *_h_val;          /* H (halo) coefficients */

  /* Pointers to auxiliary matrix structure elements */

  cs_lnum_t        *d_idx;           /* Index for diagonal matrix coefficients
                                        in case of multiple block sizes */

} cs_matrix_coeff_dist_t;

/* Matrix structure (representation-independent part) */
/*----------------------------------------------------*/

struct _cs_matrix_structure_t {

  cs_matrix_type_t       type;         /* Matrix storage and definition type */

  cs_lnum_t              n_rows;       /* Local number of rows */
  cs_lnum_t              n_cols_ext;   /* Local number of columns + ghosts */

  cs_alloc_mode_t        alloc_mode;   /* Preferred allocation mode */

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

  const char            *type_name;    /* Pointer to matrix type name string */
  const char            *type_fname;   /* Pointer to matrix type
                                          full name string */

  cs_lnum_t              n_rows;       /* Local number of rows */
  cs_lnum_t              n_cols_ext;   /* Local number of columns + ghosts */

  cs_matrix_fill_type_t  fill_type;    /* Matrix fill type */

  bool                   symmetric;    /* true if coefficients are symmetric */

  cs_lnum_t              db_size;      /* Diagonal block size */

  cs_lnum_t              eb_size;      /* Extradiag block size */

  cs_alloc_mode_t        alloc_mode;   /* Preferred allocation mode */

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

  const cs_real_t       *xa;           /* Extra-diagonal terms */

  /* Pointer to associated connectivity and geometric data needed by
     multigrid smoothing. At least cell centers and volumes are needed for
     relaxation, and face normals are needed for the "classical" option.
     Cells and faces here do not need to be primary mesh elements,
     but could be dual mesh elements of some sort */

  const cs_lnum_t    *c2f_idx;         /* Cell to faces index (shared) */
  const cs_lnum_t    *c2f;             /* Cell to faces adjacency (shared) */
  const short int    *c2f_sgn;         /* Cell to faces orientation (shared) */

  const cs_real_3_t  *cell_cen;        /* Cell center (shared) */
  const cs_real_t    *cell_vol;        /* Cell volume (shared) */
  const cs_real_3_t  *face_normal;     /* Face normal (shared) */

  /* Pointer to private data */

  void                  *coeffs;       /* Matrix coefficients */

  void                  *ext_lib_map;  /* Mapping of structure and
                                          coefficients to external
                                          library, if needed */

  /* Function pointers */

  cs_matrix_set_coeffs_t               *set_coefficients;
  cs_matrix_release_coeffs_t           *release_coefficients;
  cs_matrix_copy_diagonal_t            *copy_diagonal;
  cs_matrix_get_diagonal_t             *get_diagonal;

  cs_matrix_destroy_struct_t           *destroy_structure;
  cs_matrix_destroy_coeffs_t           *destroy_coefficients;
  cs_matrix_destroy_adaptor_t          *destroy_adaptor;

  cs_matrix_assembler_values_create_t  *assembler_values_create;

  /* Function pointer arrays, with CS_MATRIX_N_FILL_TYPES variants:
     fill_type*4 + full, extra-diagonal, lower, upper parts */

  cs_matrix_vector_product_t  *vector_multiply[CS_MATRIX_N_FILL_TYPES]
                                              [CS_MATRIX_SPMV_N_TYPES];

#if defined(HAVE_ACCEL)

  /* Function pointer arrays, with CS_MATRIX_N_FILL_TYPES variants, for host
     or device only: fill_type*4 + full, extra-diagonal, lower, upper parts */

  cs_matrix_vector_product_t  *vector_multiply_h[CS_MATRIX_N_FILL_TYPES]
                                                [CS_MATRIX_SPMV_N_TYPES];
  cs_matrix_vector_product_t  *vector_multiply_d[CS_MATRIX_N_FILL_TYPES]
                                                [CS_MATRIX_SPMV_N_TYPES];

#endif /* defined(HAVE_ACCEL) */

};

/* Structure used for tuning variants */
/*------------------------------------*/

struct _cs_matrix_variant_t {

  char                   name[2][32]; /* Variant names */

  cs_matrix_type_t       type;        /* Matrix storage and definition type */
  cs_matrix_fill_type_t  fill_type;   /* Matrix storage and definition type */

  /* Function pointer arrays */

  cs_matrix_vector_product_t   *vector_multiply[CS_MATRIX_SPMV_N_TYPES];

  /* Associated vector host/device locations */

  char                          vector_multiply_xy_hd[CS_MATRIX_SPMV_N_TYPES];
};

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Semi-private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MATRIX_PRIV_H__ */
