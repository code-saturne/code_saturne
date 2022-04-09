#ifndef __CS_MATRIX_H__
#define __CS_MATRIX_H__

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "cs_halo.h"
#include "cs_numbering.h"
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

  CS_MATRIX_NATIVE,           /*!< Native (edge-based) matrix storage */
  CS_MATRIX_CSR,              /*!< Compressed Sparse Row storage */
  CS_MATRIX_MSR,              /*!< Modified Compressed Sparse Row storage
                                   (separate diagonal) */
  CS_MATRIX_DIST,             /*!< Distributed matrix storage
                                   (separate diagonal, off-diagonal, and
                                   distant coefficients) */

  CS_MATRIX_N_BUILTIN_TYPES,  /*!< Number of known and built-in matrix types */

  CS_MATRIX_N_TYPES           /*!< Number of known matrix types */

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

/*! SpMV operation types.
 *
 * Assuming a matrix is a sum of a diagonal (D), and extra-diagonal (E) part
 * - A = D + E
 */

typedef enum {

  CS_MATRIX_SPMV,             /*!< y ← A.x (full matrix)*/
  CS_MATRIX_SPMV_E,           /*!< y ← E.x (extra-diagonal only) */

  CS_MATRIX_SPMV_N_TYPES      /*!< Number of handled SpMV operations */

} cs_matrix_spmv_type_t;

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

/*============================================================================
 *  Global variables
 *============================================================================*/

/*! Fill type names for matrices */

extern const char  *cs_matrix_fill_type_name[];

/*! Operation type type names for partial SpMV functions */

extern const char  *cs_matrix_spmv_type_name[];

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
 * Note that the resulting matrix structure will contain a full main
 * main diagonal, and that the extra-diagonal structure is always
 * symmetric (though the coefficients my not be, and we may choose a
 * matrix format that does not exploit this symmetry). If the edges
 * connectivity argument is NULL, the matrix will be purely diagonal.
 *
 * parameters:
 *   type        <-- type of matrix considered
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
 * Create an MSR matrix structure sharing an existing connectivity definition.
 *
 * Note that as the structure created maps to the given existing
 * cell global number, face -> cell connectivity arrays, and cell halo
 * structure, it must be destroyed before they are freed
 * (usually along with the code's main face -> cell structure).
 *
 * parameters:
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
 *   a pointer to a created matrix structure
 *----------------------------------------------------------------------------*/

cs_matrix_structure_t *
cs_matrix_structure_create_msr_shared(bool                    direct_assmbly,
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
 * \return a pointer to a created matrix structure
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
                                cs_matrix_assembler_t  *ma);

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
cs_matrix_create_by_local_restrict(const cs_matrix_t  *src);

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
 * Return matrix type name.
 *
 * parameters:
 *   matrix --> pointer to matrix structure
 *----------------------------------------------------------------------------*/

const char *
cs_matrix_get_type_name(const cs_matrix_t  *matrix);

/*----------------------------------------------------------------------------
 * Return matrix type full name.
 *
 * parameters:
 *   matrix --> pointer to matrix structure
 *----------------------------------------------------------------------------*/

const char *
cs_matrix_get_type_fullname(const cs_matrix_t  *matrix);

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
cs_matrix_get_n_columns(const cs_matrix_t  *matrix);

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
cs_matrix_get_n_rows(const cs_matrix_t  *matrix);

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
cs_matrix_get_n_entries(const cs_matrix_t  *matrix);

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
cs_matrix_get_diag_block_size(const cs_matrix_t  *matrix);

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
cs_matrix_get_extra_diag_block_size(const cs_matrix_t  *matrix);

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
cs_matrix_get_halo(const cs_matrix_t  *matrix);

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
cs_matrix_get_l_range(const cs_matrix_t  *matrix);

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
cs_matrix_get_alloc_mode(const cs_matrix_t  *matrix);

/*----------------------------------------------------------------------------*/
/*!
 *\brief Set matrix allocation mode.
 *
 * \param[in, out]  matrix      pointer to matrix structure
 * \param[in]      alloc_mode  host/device allocation mode
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_set_alloc_mode(cs_matrix_t       *matrix,
                         cs_alloc_mode_t   alloc_mode);

/*----------------------------------------------------------------------------
 * Get matrix fill type, depending on block sizes.
 *
 * parameters:
 *   symmetric              <-- indicates if matrix coefficients are symmetric
 *   diag_block_size        <-- block sizes for diagonal
 *   extra_diag_block_size  <-- block sizes for extra diagonal
 *
 * returns:
 *   matrix fill type
 *----------------------------------------------------------------------------*/

cs_matrix_fill_type_t
cs_matrix_get_fill_type(bool       symmetric,
                        cs_lnum_t  diag_block_size,
                        cs_lnum_t  extra_diag_block_size);

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
 * parameters:
 *   matrix                 <-> pointer to matrix structure
 *   symmetric              <-- indicates if matrix coefficients are symmetric
 *   diag_block_size        <-- block sizes for diagonal
 *   extra_diag_block_size  <-- block sizes for extra diagonal
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
                           cs_lnum_t           diag_block_size,
                           cs_lnum_t           extra_diag_block_size,
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
 * parameters:
 *   matrix                 <-> pointer to matrix structure
 *   symmetric              <-- indicates if matrix coefficients are symmetric
 *   diag_block_size        <-- block sizes for diagonal
 *   extra_diag_block_size  <-- block sizes for extra diagonal
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
                            cs_lnum_t           diag_block_size,
                            cs_lnum_t           extra_diag_block_size,
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
 * parameters:
 *   matrix                 <-> pointer to matrix structure
 *   symmetric              <-- indicates if matrix coefficients are symmetric
 *   diag_block_size        <-- block sizes for diagonal
 *   extra_diag_block_size  <-- block sizes for extra diagonal
 *   row_index              <-- MSR row index (0 to n-1)
 *   col_id                 <-- MSR column id (0 to n-1)
 *   d_val                  <-> diagonal values (NULL if zero)
 *   x_val                  <-> extradiagonal values (NULL if zero)
 *----------------------------------------------------------------------------*/

void
cs_matrix_transfer_coefficients_msr(cs_matrix_t         *matrix,
                                    bool                 symmetric,
                                    cs_lnum_t            diag_block_size,
                                    cs_lnum_t            extra_diag_block_size,
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
                                cs_lnum_t     extra_diag_block_size);

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
cs_matrix_is_mapped_from_native(const cs_matrix_t  *matrix);

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
cs_matrix_get_native_arrays(const cs_matrix_t   *matrix,
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
                               const cs_real_3_t   *face_normal);

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
                               const cs_real_3_t  **face_normal);

/*----------------------------------------------------------------------------
 * Assign functions based on a variant to a given matrix.
 *
 * If the matrix variant is incompatible with the structure, it is ignored,
 * and defaults for that structure are used instead.
 *
 * parameters:
 *   m <-> associated matrix structure
 *   mv <-- associated matrix variant
 *
 * returns:
 *   pointer to created matrix structure;
 *----------------------------------------------------------------------------*/

void
cs_matrix_apply_variant(cs_matrix_t                *m,
                        const cs_matrix_variant_t  *mv);

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x
 *
 * This function includes a halo update of x prior to multiplication by A.
 *
 * parameters:
 *   matrix        --> pointer to matrix structure
 *   x             <-> multipliying vector values (ghost values updated)
 *   y             --> resulting vector
 *----------------------------------------------------------------------------*/

void
cs_matrix_vector_multiply(const cs_matrix_t   *matrix,
                          cs_real_t           *restrict x,
                          cs_real_t           *restrict y);

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
                            cs_real_t           *restrict y);

#endif /* defined(HAVE_ACCEL) */

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
 *   y      <-- resulting vector
 *----------------------------------------------------------------------------*/

void
cs_matrix_vector_multiply_nosync(const cs_matrix_t  *matrix,
                                 cs_real_t          *restrict x,
                                 cs_real_t          *restrict y);

/*----------------------------------------------------------------------------
 * Partial matrix.vector product.
 *
 * This function includes a halo update of x prior to multiplication,
 * except for the CS_MATRIX_SPMV_L operation type, which does not require it,
 * as halo adjacencies are only present and useful in the upper-diagonal part..
 *
 * parameters:
 *   matrix        <-- pointer to matrix structure
 *   op_type       <-- SpMV operation type
 *   x             <-> multipliying vector values (ghost values updated)
 *   y             --> resulting vector
 *----------------------------------------------------------------------------*/

void
cs_matrix_vector_multiply_partial(const cs_matrix_t      *matrix,
                                  cs_matrix_spmv_type_t   op_type,
                                  cs_real_t              *restrict x,
                                  cs_real_t              *restrict y);

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
                                    cs_real_t              *restrict y);

#endif /* defined(HAVE_ACCEL) */

/*----------------------------------------------------------------------------
 * Synchronize ghost values prior to matrix.vector product
 *
 * parameters:
 *   matrix        <-- pointer to matrix structure
 *   x             <-> multipliying vector values (ghost values updated)
 *----------------------------------------------------------------------------*/

void
cs_matrix_pre_vector_multiply_sync(const cs_matrix_t   *matrix,
                                   cs_real_t           *x);

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
                             cs_matrix_variant_t    **m_variant);

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
cs_matrix_variant_create(cs_matrix_t  *m);

/*----------------------------------------------------------------------------
 * Destroy a matrix variant structure.
 *
 * parameters:
 *   mv <-> Pointer to matrix variant pointer
 *----------------------------------------------------------------------------*/

void
cs_matrix_variant_destroy(cs_matrix_variant_t  **mv);

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
                        cs_matrix_variant_t  *mv);

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
                              cs_matrix_variant_t  *mv);

/*----------------------------------------------------------------------------
 * Select the sparse matrix-vector product function to be used by a
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
 *----------------------------------------------------------------------------*/

void
cs_matrix_variant_set_func(cs_matrix_variant_t     *mv,
                           cs_matrix_fill_type_t    fill_type,
                           cs_matrix_spmv_type_t    spmv_type,
                           const cs_numbering_t    *numbering,
                           const char              *func_name);

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
cs_matrix_variant_type(const cs_matrix_variant_t  *mv);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MATRIX_H__ */
