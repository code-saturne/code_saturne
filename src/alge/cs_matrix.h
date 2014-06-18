#ifndef __CS_MATRIX_H__
#define __CS_MATRIX_H__

/*============================================================================
 * Sparse Matrix Representation and Operations
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "cs_halo.h"
#include "cs_numbering.h"
#include "cs_halo_perio.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Matrix structure representation types */

typedef enum {

  CS_MATRIX_NATIVE,     /* Native matrix format */
  CS_MATRIX_CSR,        /* Compressed Sparse Row storage format */
  CS_MATRIX_CSR_SYM,    /* Compressed Symmetric Sparse Row storage format */
  CS_MATRIX_MSR,        /* Modified Compressed Sparse Row storage format */
  CS_MATRIX_N_TYPES     /* Number of known matrix types */

} cs_matrix_type_t;

/* Matrix fill types (for tuning) */

typedef enum {

  CS_MATRIX_SCALAR,           /* Simple calar matrix */
  CS_MATRIX_SCALAR_SYM,       /* Simple scalar symmetric matrix */
  CS_MATRIX_33_BLOCK_D,       /* Matrix with 3x3 diagonal blocks
                                 (and 3.I extradiagonal blocks) */
  CS_MATRIX_33_BLOCK_D_SYM,   /* Symmetric matrix with 3x3 diagonal blocks
                                 (and 3.I extradiagonal blocks) */
  CS_MATRIX_33_BLOCK,         /* Matrix with 3x3 blocks
                                 (diagonal and extra-diagonal) */
  CS_MATRIX_N_FILL_TYPES      /* Number of possible matrix fill types */

} cs_matrix_fill_type_t;

/* Structure associated with opaque matrix structure object */

typedef struct _cs_matrix_structure_t cs_matrix_structure_t;

/* Structure associated with opaque matrix object */

typedef struct _cs_matrix_t cs_matrix_t;

/* Structure associated with opaque matrix tuning results object */

typedef struct _cs_matrix_variant_t cs_matrix_variant_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Short names for matrix types */

extern const char  *cs_matrix_type_name[];

/* Full names for matrix types */

extern const char  *cs_matrix_type_fullname[];

/* Fill type names for matrices */

extern const char  *cs_matrix_fill_type_name[];

/*=============================================================================
 * Public function prototypes
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
 *   face_cell   <-- Face -> cells connectivity
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
                           const cs_numbering_t  *numbering);

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
                            const cs_matrix_variant_t    *mv);

/*----------------------------------------------------------------------------
 * Destroy a matrix structure.
 *
 * parameters:
 *   ms <-> Pointer to matrix structure pointer
 *----------------------------------------------------------------------------*/

void
cs_matrix_structure_destroy(cs_matrix_structure_t  **ms);

/*----------------------------------------------------------------------------
 * Get the type associated with a matrix structure.
 *
 * parameters:
 *   ms <-- Associated matrix structure
 *----------------------------------------------------------------------------*/

cs_matrix_type_t
cs_matrix_structure_type(const cs_matrix_structure_t  *ms);

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
cs_matrix_create(const cs_matrix_structure_t  *ms);

/*----------------------------------------------------------------------------
 * Destroy a matrix structure.
 *
 * parameters:
 *   matrix <-> Pointer to matrix structure pointer
 *----------------------------------------------------------------------------*/

void
cs_matrix_destroy(cs_matrix_t **matrix);

/*----------------------------------------------------------------------------
 * Return number of columns in matrix.
 *
 * parameters:
 *   matrix --> Pointer to matrix structure
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_matrix_get_n_columns(const cs_matrix_t  *matrix);

/*----------------------------------------------------------------------------
 * Return number of rows in matrix.
 *
 * parameters:
 *   matrix --> Pointer to matrix structure
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_matrix_get_n_rows(const cs_matrix_t  *matrix);

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
cs_matrix_get_diag_block_size(const cs_matrix_t  *matrix);

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
                        const int  *extra_diag_block_size);

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
                           const cs_real_t  *xa);

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
                              const cs_real_t  *xa);

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
                            const cs_real_t  *xa);

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
cs_matrix_release_coefficients(cs_matrix_t  *matrix);

/*----------------------------------------------------------------------------
 * Copy matrix diagonal values.
 *
 * In case of matrixes with block diagonal coefficients, only the true
 * diagonal values are copied.
 *
 * parameters:
 *   matrix --> Pointer to matrix structure
 *   da     --> Diagonal (pre-allocated, size: n_cells)
 *----------------------------------------------------------------------------*/

void
cs_matrix_copy_diagonal(const cs_matrix_t  *matrix,
                        cs_real_t          *restrict da);

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
cs_matrix_get_diagonal(const cs_matrix_t  *matrix);

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x
 *
 * This function includes a halo update of x prior to multiplication by A.
 *
 * parameters:
 *   rotation_mode --> Halo update option for rotational periodicity
 *   matrix        --> Pointer to matrix structure
 *   x             <-> Multipliying vector values (ghost values updated)
 *   y             --> Resulting vector
 *----------------------------------------------------------------------------*/

void
cs_matrix_vector_multiply(cs_halo_rotation_t   rotation_mode,
                          const cs_matrix_t   *matrix,
                          cs_real_t           *restrict x,
                          cs_real_t           *restrict y);

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with no prior halo update of x.
 *
 * This function does not include a halo update of x prior to multiplication
 * by A, so it should be called only when the halo of x is known to already
 * be up to date (in which case we avoid the performance penalty of a
 * redundant update by using this variant of the matrix.vector product).
 *
 * parameters:
 *   matrix --> Pointer to matrix structure
 *   x      --> Multipliying vector values
 *   y      --> Resulting vector
 *----------------------------------------------------------------------------*/

void
cs_matrix_vector_multiply_nosync(const cs_matrix_t  *matrix,
                                 const cs_real_t    *x,
                                 cs_real_t          *restrict y);

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
                                 cs_real_t           *restrict y);

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
                             cs_matrix_variant_t    **m_variant);

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
                         const cs_numbering_t    *numbering);

/*----------------------------------------------------------------------------
 * Destroy a matrix variant structure.
 *
 * parameters:
 *   mv <-> Pointer to matrix variant pointer
 *----------------------------------------------------------------------------*/

void
cs_matrix_variant_destroy(cs_matrix_variant_t  **mv);

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
                           const char              *func_name);

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
                        cs_matrix_fill_type_t       fill_type);

/*----------------------------------------------------------------------------
 * Get the type associated with a matrix variant.
 *
 * parameters:
 *   mv <-- Pointer to matrix variant structure
 *----------------------------------------------------------------------------*/

cs_matrix_type_t
cs_matrix_variant_type(const cs_matrix_variant_t  *mv);

/*----------------------------------------------------------------------------
 * Test local matrix.vector product operations.
 *
 * parameters:
 *   n_cells        <-- number of local cells
 *   n_cells_ext    <-- number of cells including ghost cells (array size)
 *   n_faces        <-- local number of internal faces
 *   cell_num       <-- Optional global cell numbers (1 to n), or NULL
 *   face_cell      <-- face -> cells connectivity
 *   halo           <-- cell halo structure
 *   numbering      <-- vectorization or thread-related numbering info, or NULL
 *----------------------------------------------------------------------------*/

void
cs_matrix_variant_test(cs_lnum_t              n_cells,
                       cs_lnum_t              n_cells_ext,
                       cs_lnum_t              n_faces,
                       const cs_gnum_t       *cell_num,
                       const cs_lnum_2_t     *face_cell,
                       const cs_halo_t       *halo,
                       const cs_numbering_t  *numbering);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MATRIX_H__ */
