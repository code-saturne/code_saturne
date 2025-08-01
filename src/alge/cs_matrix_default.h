#ifndef __CS_MATRIX_DEFAULT_H__
#define __CS_MATRIX_DEFAULT_H__

/*============================================================================
 * Default Sparse Matrix structure and Tuning.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

#include "base/cs_field.h"
#include "base/cs_halo.h"
#include "alge/cs_matrix.h"
#include "base/cs_numbering.h"
#include "base/cs_halo_perio.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Matrix (native format) vector product
 *
 * parameters:
 *   symmetric     <-- Symmetry indicator:
 *   db_size       <-- block sizes for diagonal
 *   eb_size       <-- block sizes for extra diagonal
 *   f_id          <-- associated field id, or < 0
 *   dam           <-- Matrix diagonal
 *   xam           <-- Matrix extra-diagonal terms
 *   vx            <-- A*vx
 *   vy            <-> vy = A*vx
 *----------------------------------------------------------------------------*/

void
cs_matrix_vector_native_multiply(bool              symmetric,
                                 cs_lnum_t         db_size,
                                 cs_lnum_t         eb_size,
                                 int               f_id,
                                 const cs_real_t  *dam,
                                 const cs_real_t  *xam,
                                 cs_real_t        *vx,
                                 cs_real_t        *vy);

/*----------------------------------------------------------------------------
 * Initialize sparse matrix API.
 *----------------------------------------------------------------------------*/

void
cs_matrix_initialize(void);

/*----------------------------------------------------------------------------
 * Finalize sparse matrix API.
 *----------------------------------------------------------------------------*/

void
cs_matrix_finalize(void);

/*----------------------------------------------------------------------------
 * Update sparse matrix API in case of mesh modification.
 *----------------------------------------------------------------------------*/

void
cs_matrix_update_mesh(void);

/*----------------------------------------------------------------------------
 * Return default matrix for a given fill type
 *
 * parameters:
 *   symmetric              <-- Indicates if matrix coefficients are symmetric
 *   diag_block_size        <-- Block sizes for diagonal
 *   extra_diag_block_size  <-- Block sizes for extra diagonal
 *
 * returns:
 *   pointer to default matrix structure adapted to fill type
 *----------------------------------------------------------------------------*/

cs_matrix_t  *
cs_matrix_default(bool       symmetric,
                  cs_lnum_t  diag_block_size,
                  cs_lnum_t  extra_diag_block_size);

/*----------------------------------------------------------------------------
 * Return MSR matrix
 *
 * returns:
 *   pointer to MSR matrix
 *----------------------------------------------------------------------------*/

cs_matrix_t  *
cs_matrix_msr(void);

/*----------------------------------------------------------------------------
 * Return CSR matrix for a given fill type
 *
 * returns:
 *   pointer to CSR matrix
 *----------------------------------------------------------------------------*/

cs_matrix_t  *
cs_matrix_csr(void);

/*----------------------------------------------------------------------------
 * Return native matrix
 *
 * returns:
 *   pointer to native matrix
 *----------------------------------------------------------------------------*/

cs_matrix_t  *
cs_matrix_native(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return matrix wrapper for external library for a given fill type.
 *
 * \param[in]  type_name              Matrix type name
 * \param[in]  symmetric              Indicates if coefficients are symmetric
 * \param[in]  diag_block_size        Block sizes for diagonal
 * \param[in]  extra_diag_block_size  Block sizes for extra diagonal
 *
 * \return  Pointer to matrix matching requested type
 */
/*----------------------------------------------------------------------------*/

cs_matrix_t  *
cs_matrix_external(const char  *type_name,
                   bool         symmetric,
                   cs_lnum_t    diag_block_size,
                   cs_lnum_t    extra_diag_block_size);

/*----------------------------------------------------------------------------*/
/*
 * \brief Copy base matrix to external library matrix type for given fill type.
 *
 * Note that the matrix containers share the same assigned structure,
 * so they must be both destroyed before that structure.
 *
 * Coefficients and matching structures are not copied or created.
 *
 * This function is intended to allow sharing of a base structure or assembler
 * with an external library matrix wrapper, so as to allow efficient
 * coefficient assignment, but with external coefficient handling.
 *
 * The matrix shoud be converted to the desired external type after calling
 * this function, so that it can the be accessed using \ref cs_matrix_external.
 *
 * \param[in]  symmetric              Indicates if matrix coefficients are symmetric
 * \param[in]  diag_block_size        Block sizes for diagonal
 * \param[in]  extra_diag_block_size  Block sizes for extra diagonal
 *
 * \return  pointer to native matrix adapted to fill type
 */
/*----------------------------------------------------------------------------*/

cs_matrix_t  *
cs_matrix_copy_to_external(cs_matrix_t  *src,
                           bool          symmetric,
                           cs_lnum_t     diag_block_size,
                           cs_lnum_t     extra_diag_block_size);

/*----------------------------------------------------------------------------
 * Determine or apply default tuning for a given matrix type
 *
 * Information from the variant used fo this definition is copied,
 * so it may be freed after calling this function.
 *
 * parameters:
 *   fill type  <-- Fill type for which tuning behavior is set
 *   mv         <-- Matrix variant to use for this type
 *----------------------------------------------------------------------------*/

void
cs_matrix_default_set_tuned(cs_matrix_t  *m);

/*----------------------------------------------------------------------------
 * Set number of matrix computation runs for tuning.
 *
 * If this function is not called, defaults are:
 *  - minimum of 10 runs
 *
 * parameters:
 *   n_min_products <-- minimum number of SpM.V products for tuning.
 *----------------------------------------------------------------------------*/

void
cs_matrix_set_tuning_runs(int     n_min_products);

/*----------------------------------------------------------------------------
 * Get number of matrix computation runs for tuning.
 *
 * return:
 *   minimum number of SpM.V calls for tuning
 *----------------------------------------------------------------------------*/

int
cs_matrix_get_tuning_runs(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Set default matrix type for a given fill type.
 *
 * \param[in] fill type  Fill type for which tuning behavior is set
 * \param[in] type       Matrix type to use
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_default_set_type(cs_matrix_fill_type_t  fill_type,
                           cs_matrix_type_t       type);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return a (0-based) global block row numbering for a given matrix.
 *
 * The numbering is built or updated if not previously used, or if the
 * previous call considered a different matrix or halo, and is simply
 * returned otherwise.
 * In other words, this works as a matrix global numbering cache.
 *
 * The matrix's halo is used for the update.
 *
 * \param[in]  m  associated matrix
 *
 * \return  pointer to requested global numbering
 */
/*----------------------------------------------------------------------------*/

const cs_gnum_t *
cs_matrix_get_block_row_g_id(const cs_matrix_t  *m);

END_C_DECLS
#ifdef __cplusplus

/*----------------------------------------------------------------------------*/
/*
 * \brief Return a (0-based) global block row numbering for a given matrix.
 *
 * The numbering is built or updated if not previously used, or if the
 * previous call considered a different matrix or halo, and is simply
 * returned otherwise.
 * In other words, this works as a matrix global numbering cache.
 *
 * \param[in]  m     associated matrix
 * \param[in]  halo  associated halo
 *
 * \return  pointer to requested global numbering
 */
/*----------------------------------------------------------------------------*/

const cs_gnum_t *
cs_matrix_get_block_row_g_id(const cs_matrix_t  *m,
                             const cs_halo_t    *halo);

#endif /* cplusplus */
BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*
 * \brief Return matrix associated wiht a matrix assembler.
 *
 * Coefficients are not assigned at this stage.
 *
 * \param[in]  f                      pointer to associated field
 * \param[in]  type                   matrix type
 *
 * \return  pointer to associated matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_matrix_t *
cs_matrix_by_assembler(const cs_field_t  *f,
                       cs_matrix_type_t   type);

/*----------------------------------------------------------------------------*/
/*
 * \brief Assign coefficients to a matrix using a matrix assembler.
 *
 * \param[in]  f                      pointer to associated field
 * \param[in]  type                   matrix type
 * \param[in]  symmetric              is matrix symmetric ?
 * \param[in]  diag_block_size        block sizes for diagonal, or nullptr
 * \param[in]  extra_diag_block_size  block sizes for extra diagonal, or nullptr
 * \param[in]  da                     diagonal values (nullptr if zero)
 * \param[in]  xa                     extradiagonal values (nullptr if zero)
 *                                    casts as:
 *                                      xa[n_edges]    if symmetric,
 *                                      xa[n_edges][2] if non symmetric
 *
 * \return  pointer to associated matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_matrix_t *
cs_matrix_set_coefficients_by_assembler(const cs_field_t  *f,
                                        cs_matrix_type_t   type,
                                        bool               symmetric,
                                        cs_lnum_t          diag_block_size,
                                        cs_lnum_t          extra_diag_block_size,
                                        const cs_real_t   *da,
                                        const cs_real_t   *xa);

/*----------------------------------------------------------------------------
 * Release of destroy matrix depending on whether is is cached or not.
 *
 * Matrices built by assembler are destroyed.
 *
 * parameters:
 *   matrix <-> pointer to matrix structure pointer
 *----------------------------------------------------------------------------*/

void
cs_matrix_release(cs_matrix_t  **m);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MATRIX_DEFAULT_H__ */
