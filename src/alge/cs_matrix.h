#ifndef __CS_MATRIX_H__
#define __CS_MATRIX_H__

/*============================================================================
 * Sparse Matrix Representation and Operations
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "cs_halo.h"
#include "cs_numbering.h"
#include "cs_halo_perio.h"
#include "cs_matrix_assembler.h"

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

  CS_MATRIX_SCALAR,           /* Simple scalar matrix */
  CS_MATRIX_SCALAR_SYM,       /* Simple scalar symmetric matrix */
  CS_MATRIX_BLOCK_D,          /* Matrix with diagonal blocks
                                 (and m.I extradiagonal blocks) */
  CS_MATRIX_BLOCK_D_66,       /* Matrix with 6x6 diagonal blocks
                                 (and 6.I extradiagonal blocks;
                                 subcase of CS_MATRIX_BLOCK_D, allows
                                 separate tuning) */
  CS_MATRIX_BLOCK_D_SYM,      /* Symmetric matrix with diagonal blocks
                                 (and m.I extradiagonal blocks) */
  CS_MATRIX_BLOCK,            /* Block matrix */
  CS_MATRIX_N_FILL_TYPES      /* Number of possible matrix fill types */

} cs_matrix_fill_type_t;

/* Structure associated with opaque matrix structure object */

typedef struct _cs_matrix_structure_t cs_matrix_structure_t;

/* Structure associated with opaque matrix object */

typedef struct _cs_matrix_t cs_matrix_t;

/* Structure associated with opaque matrix tuning results object */

typedef struct _cs_matrix_variant_t cs_matrix_variant_t;

/* Information structure for extraction of matrix row */

typedef struct {

  cs_lnum_t          row_size;       /*< Row size from last call */
  cs_lnum_t          buffer_size;    /*< Allocated buffer size */
  const cs_lnum_t   *col_id;         /*< Pointer to local column ids */
  cs_lnum_t        *_col_id;         /*< Pointer to local column ids copy */
  const cs_real_t   *vals;           /*< Pointer to local row values */
  cs_real_t        *_vals;           /*< Pointer to local row values copy */

} cs_matrix_row_info_t;

/*! Function pointers */
/* -------------------*/

/*----------------------------------------------------------------------------
 * Matrix.vector product contribution y = A'.x with no prior halo update of x.
 *
 * This function does not include a halo update of x prior to multiplication
 * by A', so it should be called only when the halo of x is known to already
 * be up to date (in which case we avoid the performance penalty of a
 * redundant update by using this variant of the matrix.vector product).
 *
 * parameters:
 *   exclude_diag <-- if true, exclude diagonal from product
 *   input        <-- pointer to additional data
 *   x            <-- multipliying vector values
 *   y            <-> resulting vector
 *----------------------------------------------------------------------------*/

typedef void
(cs_matrix_vector_product_extend_t) (bool                exclude_diag,
                                     void               *input,
                                     const cs_real_t    *restrict x,
                                     cs_real_t          *restrict y);

/*----------------------------------------------------------------------------
 * Contribution to diagonal preconditionning associated with extended SpMV.
 *
 * parameters:
 *   input        <-- pointer to additional data
 *   ad           <-> diagonal part of linear equation matrix
 *----------------------------------------------------------------------------*/

typedef void
(cs_matrix_preconditioner_extend_t) (void               *input,
                                      cs_real_t          *restrict ad);

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
 * Create a matrix structure.
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
 * parameters:
 *   type        <-- type of matrix considered
 *   have_diag   <-- indicates if the diagonal structure contains nonzeroes
 *   n_rows      <-- local number of rows
 *   n_cols_ext  <-- number of columns + ghosts
 *   n_edges     <-- local number of (undirected) graph edges
 *   edges       <-- edges (symmetric row <-> column) connectivity
 *   halo        <-- halo structure associated with cells, or NULL
 *   numbering   <-- vectorization or thread-related numbering info, or NULL
 *
 * returns:
 *   pointer to created matrix structure;
 *----------------------------------------------------------------------------*/

cs_matrix_structure_t *
cs_matrix_structure_create(cs_matrix_type_t       type,
                           bool                   have_diag,
                           cs_lnum_t              n_rows,
                           cs_lnum_t              n_cols_ext,
                           cs_lnum_t              n_edges,
                           const cs_lnum_2_t     *edges,
                           const cs_halo_t       *halo,
                           const cs_numbering_t  *numbering);

/*----------------------------------------------------------------------------
 * Create a matrix structure based on a MSR connectivity definition.
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
 * parameters:
 *   type       <-- type of matrix considered
 *   transfer   <-- transfer property of row_index and col_id
 *                  if true, map them otherwise
 *   have_diag  <-- indicates if the structure includes the
 *                  diagonal (should be the same for all rows)
 *   n_rows     <-- local number of rows
 *   n_cols_ext <-- local number of columns + ghosts
 *   row_index  <-> pointer to index on rows
 *   col_id     <-> pointer to array of colum ids related to the row index
 *   halo       <-- halo structure for synchronization, or NULL
 *   numbering  <-- vectorization or thread-related numbering info, or NULL
 *
 * returns:
 *   a pointer to a created matrix structure
 *----------------------------------------------------------------------------*/

cs_matrix_structure_t *
cs_matrix_structure_create_msr(cs_matrix_type_t        type,
                               bool                    transfer,
                               bool                    have_diag,
                               cs_lnum_t               n_rows,
                               cs_lnum_t               n_cols_ext,
                               cs_lnum_t             **row_index,
                               cs_lnum_t             **col_id,
                               const cs_halo_t        *halo,
                               const cs_numbering_t   *numbering);

/*----------------------------------------------------------------------------
 * Create an MSR matrix structure sharing an existing connectivity definition
 * as well as an optional edge-based definition.
 *
 * Note that as the structure created maps to the given existing
 * cell global number, face -> cell connectivity arrays, and cell halo
 * structure, it must be destroyed before they are freed
 * (usually along with the code's main face -> cell structure).
 *
 * parameters:
 *   have_diag        <-- indicates if the structure includes the
 *                        diagonal (should be the same for all rows)
 *   direct_assembly  <-- true if each value corresponds to a unique face
 *   n_rows           <-- local number of rows
 *   n_cols_ext       <-- local number of columns + ghosts
 *   row_index        <-- pointer to index on rows
 *   col_id           <-- pointer to array of colum ids related to the row index
 *   halo             <-- halo structure for synchronization, or NULL
 *   numbering        <-- vectorization or thread-related numbering
 *                        info, or NULL
 *
 * returns:
 *   a pointer to a created CDO matrix structure
 *----------------------------------------------------------------------------*/

cs_matrix_structure_t *
cs_matrix_structure_create_msr_shared(bool                    have_diag,
                                      bool                    direct_assmbly,
                                      cs_lnum_t               n_rows,
                                      cs_lnum_t               n_cols_ext,
                                      const cs_lnum_t        *row_index,
                                      const cs_lnum_t        *col_id,
                                      const cs_halo_t        *halo,
                                      const cs_numbering_t   *numbering);

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
                                          cs_matrix_assembler_t  *ma);

/*----------------------------------------------------------------------------
 * Destroy a matrix structure.
 *
 * parameters:
 *   ms <-> pointer to matrix structure pointer
 *----------------------------------------------------------------------------*/

void
cs_matrix_structure_destroy(cs_matrix_structure_t  **ms);

/*----------------------------------------------------------------------------
 * Create a matrix container using a given structure.
 *
 * Note that the matrix container maps to the assigned structure,
 * so it must be destroyed before that structure.
 *
 * parameters:
 *   ms <-- associated matrix structure
 *
 * returns:
 *   pointer to created matrix structure;
 *----------------------------------------------------------------------------*/

cs_matrix_t *
cs_matrix_create(const cs_matrix_structure_t  *ms);

/*----------------------------------------------------------------------------
 * Create a matrix container using a given variant.
 *
 * If the matrix variant is incompatible with the structure, it is ignored,
 * and defaults for that structure are used instead.
 *
 * parameters:
 *   ms <-- associated matrix structure
 *   mv <-- associated matrix variant
 *
 * returns:
 *   pointer to created matrix structure;
 *----------------------------------------------------------------------------*/

cs_matrix_t *
cs_matrix_create_by_variant(const cs_matrix_structure_t  *ms,
                            const cs_matrix_variant_t    *mv);

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
cs_matrix_create_by_copy(cs_matrix_t   *src);

/*----------------------------------------------------------------------------
 * Destroy a matrix structure.
 *
 * In the case of a compoud matrix, sub-matrices are not destroyed.
 *
 * parameters:
 *   matrix <-> pointer to matrix structure pointer
 *----------------------------------------------------------------------------*/

void
cs_matrix_destroy(cs_matrix_t **matrix);

/*----------------------------------------------------------------------------
 * Return type of matrix.
 *
 * parameters:
 *   matrix --> pointer to matrix structure
 *----------------------------------------------------------------------------*/

cs_matrix_type_t
cs_matrix_get_type(const cs_matrix_t  *matrix);

/*----------------------------------------------------------------------------
 * Return number of columns in matrix.
 *
 * parameters:
 *   matrix --> pointer to matrix structure
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_matrix_get_n_columns(const cs_matrix_t  *matrix);

/*----------------------------------------------------------------------------
 * Return number of rows in matrix.
 *
 * parameters:
 *   matrix --> pointer to matrix structure
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
 *   matrix <-- pointer to matrix structure
 *
 * returns:
 *   pointer to block sizes
 *----------------------------------------------------------------------------*/

const int *
cs_matrix_get_diag_block_size(const cs_matrix_t  *matrix);

/*----------------------------------------------------------------------------
 * Return matrix extra-diagonal block sizes.
 *
 * Block sizes are defined by a array of 4 values:
 *   0: useful block size, 1: vector block extents,
 *   2: matrix line extents,  3: matrix line*column extents
 *
 * parameters:
 *   matrix <-- pointer to matrix structure
 *
 * returns:
 *   pointer to block sizes
 *----------------------------------------------------------------------------*/

const int *
cs_matrix_get_extra_diag_block_size(const cs_matrix_t  *matrix);

/*----------------------------------------------------------------------------
 * Return pointer to matrix halo structure.
 *
 * parameters:
 *   matrix <-- pointer to matrix structure
 *
 * returns:
 *   pointer to halo strucuture
 *----------------------------------------------------------------------------*/

const cs_halo_t *
cs_matrix_get_halo(const cs_matrix_t  *matrix);

/*----------------------------------------------------------------------------
 * Get matrix fill type, depending on block sizes.
 *
 * Block sizes are defined by an optional array of 4 values:
 *   0: useful block size, 1: vector block extents,
 *   2: matrix line extents,  3: matrix line*column extents
 *
 * parameters:
 *   symmetric              <-- indicates if matrix coefficients are symmetric
 *   diag_block_size        <-- block sizes for diagonal, or NULL
 *   extra_diag_block_size  <-- block sizes for extra diagonal, or NULL
 *
 * returns:
 *   matrix fill type
 *----------------------------------------------------------------------------*/

cs_matrix_fill_type_t
cs_matrix_get_fill_type(bool        symmetric,
                        const int  *diag_block_size,
                        const int  *extra_diag_block_size);

/*----------------------------------------------------------------------------
 * Set matrix coefficients defined relative to a "native" edge graph,
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
 * parameters:
 *   matrix                 <-> pointer to matrix structure
 *   symmetric              <-- indicates if matrix coefficients are symmetric
 *   diag_block_size        <-- block sizes for diagonal, or NULL
 *   extra_diag_block_size  <-- block sizes for extra diagonal, or NULL
 *   n_edges                <-- local number of graph edges
 *   edges                  <-- edges (row <-> column) connectivity
 *   da                     <-- diagonal values (NULL if zero)
 *   xa                     <-- extradiagonal values (NULL if zero)
 *                              casts as:
 *                                xa[n_edges]    if symmetric,
 *                                xa[n_edges][2] if non symmetric
 *----------------------------------------------------------------------------*/

void
cs_matrix_set_coefficients(cs_matrix_t        *matrix,
                           bool                symmetric,
                           const int          *diag_block_size,
                           const int          *extra_diag_block_size,
                           const cs_lnum_t     n_edges,
                           const cs_lnum_2_t   edges[],
                           const cs_real_t    *da,
                           const cs_real_t    *xa);

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
 *   matrix                 <-> pointer to matrix structure
 *   symmetric              <-- indicates if matrix coefficients are symmetric
 *   diag_block_size        <-- block sizes for diagonal, or NULL
 *   extra_diag_block_size  <-- block sizes for extra diagonal, or NULL
 *   n_edges                <-- local number of graph edges
 *   edges                  <-- edges (row <-> column) connectivity
 *   da                     <-- diagonal values (NULL if zero)
 *   xa                     <-- extradiagonal values (NULL if zero)
 *                              casts as:
 *                                xa[n_edges]    if symmetric,
 *                                xa[n_edges][2] if non symmetric
 *----------------------------------------------------------------------------*/

void
cs_matrix_copy_coefficients(cs_matrix_t        *matrix,
                            bool                symmetric,
                            const int          *diag_block_size,
                            const int          *extra_diag_block_size,
                            const cs_lnum_t     n_edges,
                            const cs_lnum_2_t   edges[],
                            const cs_real_t    *da,
                            const cs_real_t    *xa);

/*----------------------------------------------------------------------------
 * Set matrix coefficients in an MSR format, transferring the
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
 * parameters:
 *   matrix                 <-> pointer to matrix structure
 *   symmetric              <-- indicates if matrix coefficients are symmetric
 *   diag_block_size        <-- block sizes for diagonal, or NULL
 *   extra_diag_block_size  <-- block sizes for extra diagonal, or NULL
 *   row_index              <-- MSR row index (0 to n-1)
 *   col_id                 <-- MSR column id (0 to n-1)
 *   d_val                  <-> diagonal values (NULL if zero)
 *   x_val                  <-> extradiagonal values (NULL if zero)
 *----------------------------------------------------------------------------*/

void
cs_matrix_transfer_coefficients_msr(cs_matrix_t         *matrix,
                                    bool                 symmetric,
                                    const int           *diag_block_size,
                                    const int           *extra_diag_block_size,
                                    const cs_lnum_t      row_index[],
                                    const cs_lnum_t      col_id[],
                                    cs_real_t          **d_val,
                                    cs_real_t          **x_val);

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
cs_matrix_assembler_values_init(cs_matrix_t  *matrix,
                                const int    *diag_block_size,
                                const int    *extra_diag_block_size);

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
 *   matrix <-> pointer to matrix structure
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
 *   matrix --> pointer to matrix structure
 *   da     --> diagonal (pre-allocated, size: n_rows*block_size
 *----------------------------------------------------------------------------*/

void
cs_matrix_copy_diagonal(const cs_matrix_t  *matrix,
                        cs_real_t          *restrict da);

/*----------------------------------------------------------------------------
 * Query matrix coefficients symmetry
 *
 * parameters:
 *   matrix <-- pointer to matrix structure
 *
 * returns:
 *   true if coefficients are symmetric, false otherwise
 *----------------------------------------------------------------------------*/

bool
cs_matrix_is_symmetric(const cs_matrix_t  *matrix);

/*----------------------------------------------------------------------------
 * Get matrix diagonal values.
 *
 * In case of matrixes with block diagonal coefficients, a pointer to
 * the complete block diagonal is returned.
 *
 * parameters:
 *   matrix --> pointer to matrix structure
 *
 * returns:
 *   pointer to matrix diagonal array
 *----------------------------------------------------------------------------*/

const cs_real_t *
cs_matrix_get_diagonal(const cs_matrix_t  *matrix);

/*----------------------------------------------------------------------------
 * Get pointer to matrix extra-diagonal values in "native" format
 *
 * This function currently only functions if the matrix is in "native"
 * format or the coefficients were mapped from native coefficients using
 * cs_matrix_set_coefficients(), in which case the pointer returned is
 * the same as the one passed to that function.
 *
 * parameters:
 *   matrix --> pointer to matrix structure
 *
 * returns:
 *   pointer to matrix diagonal array
 *----------------------------------------------------------------------------*/

const cs_real_t *
cs_matrix_get_extra_diagonal(const cs_matrix_t  *matrix);

/*----------------------------------------------------------------------------
 * Initialize row info for a given matrix.
 *
 * parameters:
 *   r --> row info structure
 *----------------------------------------------------------------------------*/

void
cs_matrix_row_init(cs_matrix_row_info_t  *r);

/*----------------------------------------------------------------------------
 * Finalize row info for a given matrix.
 *
 * parameters:
 *   r <-> row info structure
 *----------------------------------------------------------------------------*/

void
cs_matrix_row_finalize(cs_matrix_row_info_t  *r);

/*----------------------------------------------------------------------------
 * Get row values for a given matrix.
 *
 * This function may not work for all matrix types.
 *
 * In the case of blocked matrixes, the true (non-blocked)
 * values are returned.
 *
 * The row information structure must have been previously initialized
 * using cs_matrix_row_init(), and should be finalized using
 * using cs_matrix_row_finalize(), so as to free buffers it may have
 * built for certain matrix formats.
 *
 * parameters:
 *   matrix    <-- pointer to matrix structure
 *   row_id    <-- id of row to query
 *   r         <-> row info structure
 *----------------------------------------------------------------------------*/

void
cs_matrix_get_row(const cs_matrix_t     *matrix,
                  const cs_lnum_t        row_id,
                  cs_matrix_row_info_t  *r);

/*----------------------------------------------------------------------------
 * Get arrays describing a matrix in native format.
 *
 * This function works for matrix in native format.
 *
 * Matrix block sizes can be obtained by cs_matrix_get_diag_block_size()
 * and cs_matrix_get_extra_diag_block_size().
 *
 * parameters:
 *   matrix    <-- pointer to matrix structure
 *   symmetric --> true if symmetric
 *   n_edges   --> number of associated faces
 *   edges     --> edges (symmetric row <-> column) connectivity
 *   d_val     --> diagonal values
 *   x_val     --> extra-diagonal values
 *----------------------------------------------------------------------------*/

void
cs_matrix_get_native_arrays(const cs_matrix_t    *matrix,
                            bool                *symmetric,
                            cs_lnum_t           *n_edges,
                            const cs_lnum_2_t  **edges,
                            const cs_real_t    **d_val,
                            const cs_real_t    **x_val);

/*----------------------------------------------------------------------------
 * Get arrays describing a matrix in CSR format.
 *
 * This function only works for an CSR matrix (i.e. there is
 * no automatic conversion from another matrix type).
 *
 * Matrix block sizes can be obtained by cs_matrix_get_diag_block_size()
 * and cs_matrix_get_extra_diag_block_size().
 *
 * parameters:
 *   matrix    <-- pointer to matrix structure
 *   row_index --> CSR row index
 *   col_id    --> CSR column id
 *   val       --> values
 *----------------------------------------------------------------------------*/

void
cs_matrix_get_csr_arrays(const cs_matrix_t   *matrix,
                         const cs_lnum_t    **row_index,
                         const cs_lnum_t    **col_id,
                         const cs_real_t    **val);

/*----------------------------------------------------------------------------
 * Get arrays describing a matrix in MSR format.
 *
 * This function only works for an MSR matrix (i.e. there is
 * no automatic conversion from another matrix type).
 *
 * Matrix block sizes can be obtained by cs_matrix_get_diag_block_size()
 * and cs_matrix_get_extra_diag_block_size().
 *
 * parameters:
 *   matrix    <-- pointer to matrix structure
 *   row_index --> MSR row index
 *   col_id    --> MSR column id
 *   d_val     --> diagonal values
 *   x_val     --> extra-diagonal values
 *----------------------------------------------------------------------------*/

void
cs_matrix_get_msr_arrays(const cs_matrix_t   *matrix,
                         const cs_lnum_t    **row_index,
                         const cs_lnum_t    **col_id,
                         const cs_real_t    **d_val,
                         const cs_real_t    **x_val);

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x
 *
 * This function includes a halo update of x prior to multiplication by A.
 *
 * parameters:
 *   rotation_mode --> halo update option for rotational periodicity
 *   matrix        --> pointer to matrix structure
 *   x             <-> multipliying vector values (ghost values updated)
 *   y             --> resulting vector
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
 *   matrix --> pointer to matrix structure
 *   x      --> multipliying vector values
 *   y      --> resulting vector
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
 *   rotation_mode <-- halo update option for rotational periodicity
 *   matrix        <-- pointer to matrix structure
 *   x             <-> multipliying vector values (ghost values updated)
 *   y             --> resulting vector
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
 *   type         <-- type of matrix considered
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
 *     generic         (for CS_MATRIX_??_BLOCK_D or CS_MATRIX_??_BLOCK_D_SYM)
 *     mkl             (with MKL, for CS_MATRIX_SCALAR or CS_MATRIX_SCALAR_SYM)
 *
 * parameters:
 *   mv        <-> pointer to matrix variant
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
 *   mv        <-> pointer to matrix variant
 *   mv_merge  <-- pointer to matrix variant to merge
 *   fill_type <-- matrix fill type to merge from
 *----------------------------------------------------------------------------*/

void
cs_matrix_variant_merge(cs_matrix_variant_t        *mv,
                        const cs_matrix_variant_t  *mv_merge,
                        cs_matrix_fill_type_t       fill_type);

/*----------------------------------------------------------------------------
 * Get the type associated with a matrix variant.
 *
 * parameters:
 *   mv <-- pointer to matrix variant structure
 *----------------------------------------------------------------------------*/

cs_matrix_type_t
cs_matrix_variant_type(const cs_matrix_variant_t  *mv);

/*----------------------------------------------------------------------------
 * Test local matrix.vector product operations.
 *
 * parameters:
 *   n_rows         <-- number of local rows
 *   n_cols_ext     <-- number of columns + ghosts
 *   n_edges        <-- local number of (undirected) graph edges
 *   edges          <-- edges (symmetric row <-> column) connectivity
 *   halo           <-- cell halo structure
 *   numbering      <-- vectorization or thread-related numbering info, or NULL
 *----------------------------------------------------------------------------*/

void
cs_matrix_variant_test(cs_lnum_t              n_rows,
                       cs_lnum_t              n_cols_ext,
                       cs_lnum_t              n_edges,
                       const cs_lnum_2_t     *edges,
                       const cs_halo_t       *halo,
                       const cs_numbering_t  *numbering);

/*----------------------------------------------------------------------------
 * Set coupling entity associated to matrix, or NULL if no coupling
 * entity has been set before.
 *
 * parameters:
 *   m                      <-- pointer to matrix structure
 *   vector_multiply_extend --> pointer to SpMV extension function
 *   preconditioner_extend  --> pointer to preconditioner extension function
 *   input_extend           --> pointer to associated data
 *----------------------------------------------------------------------------*/

void
cs_matrix_get_extend(const cs_matrix_t                   *m,
                     cs_matrix_vector_product_extend_t  **vector_multiply_extend,
                     cs_matrix_preconditioner_extend_t  **preconditioner_extend,
                     void                               **input_extend);

/*----------------------------------------------------------------------------
 * Set coupling entity associated to matrix, or NULL if no coupling
 * entity has been set before.
 *
 * parameters:
 *   m                      <-- pointer to matrix structure
 *   vector_multiply_extend <-- pointer to SpMV extension function
 *   preconditioner_extend  <-- pointer to preconditioner extension function
 *   input_extend           <-> pointer to associated data
 *----------------------------------------------------------------------------*/

void
cs_matrix_set_extend(cs_matrix_t                       *m,
                     cs_matrix_vector_product_extend_t *vector_multiply_extend,
                     cs_matrix_preconditioner_extend_t *preconditioner_extend,
                     void                              *input_extend);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MATRIX_H__ */
