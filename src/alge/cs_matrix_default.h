#ifndef __CS_MATRIX_DEFAULT_H__
#define __CS_MATRIX_DEFAULT_H__

/*============================================================================
 * Default Sparse Matrix structure and Tuning.
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

#include "cs_field.h"
#include "cs_halo.h"
#include "cs_matrix.h"
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
 *   rotation_mode <-- halo update option for rotational periodicity
 *   f_id          <-- associated field id, or < 0
 *   dam           <-- Matrix diagonal
 *   xam           <-- Matrix extra-diagonal terms
 *   vx            <-- A*vx
 *   vy            <-> vy = A*vx
 *----------------------------------------------------------------------------*/

void
cs_matrix_vector_native_multiply(bool                symmetric,
                                 const int           db_size[4],
                                 const int           eb_size[4],
                                 cs_halo_rotation_t  rotation_mode,
                                 int                 f_id,
                                 const cs_real_t    *dam,
                                 const cs_real_t    *xam,
                                 cs_real_t          *vx,
                                 cs_real_t          *vy);

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
 *   diag_block_size        <-- Block sizes for diagonal, or NULL
 *   extra_diag_block_size  <-- Block sizes for extra diagonal, or NULL
 *
 * returns:
 *   pointer to default matrix structure adapted to fill type
 *----------------------------------------------------------------------------*/

cs_matrix_t  *
cs_matrix_default(bool        symmetric,
                  const int  *diag_block_size,
                  const int  *extra_diag_block_size);

/*----------------------------------------------------------------------------
 * Return MSR matrix for a given fill type
 *
 * parameters:
 *   symmetric              <-- Indicates if matrix coefficients are symmetric
 *   diag_block_size        <-- Block sizes for diagonal, or NULL
 *   extra_diag_block_size  <-- Block sizes for extra diagonal, or NULL
 *
 * returns:
 *   pointer to MSR matrix adapted to fill type
 *----------------------------------------------------------------------------*/

cs_matrix_t  *
cs_matrix_msr(bool        symmetric,
              const int  *diag_block_size,
              const int  *extra_diag_block_size);

/*----------------------------------------------------------------------------
 * Return native matrix for a given fill type
 *
 * parameters:
 *   symmetric              <-- Indicates if matrix coefficients are symmetric
 *   diag_block_size        <-- Block sizes for diagonal, or NULL
 *   extra_diag_block_size  <-- Block sizes for extra diagonal, or NULL
 *
 * returns:
 *   pointer to native matrix adapted to fill type
 *----------------------------------------------------------------------------*/

cs_matrix_t  *
cs_matrix_native(bool        symmetric,
                 const int  *diag_block_size,
                 const int  *extra_diag_block_size);

/*----------------------------------------------------------------------------
 * Force matrix variant for a given fill type
 *
 * Information from the variant used fo this definition is copied,
 * so it may be freed after calling this function.
 *
 * parameters:
 *   fill type  <-- Fill type for which tuning behavior is set
 *   mv         <-- Matrix variant to use for this type
 *----------------------------------------------------------------------------*/

void
cs_matrix_set_variant(cs_matrix_fill_type_t       fill_type,
                      const cs_matrix_variant_t  *mv);

/*----------------------------------------------------------------------------
 * Set matrix tuning behavior for a given fill type
 *
 * parameters:
 *   fill type  <-- Fill type for which tuning behavior is set
 *   tune       <-- 1 to activate tuning, 0 to deactivate
 *----------------------------------------------------------------------------*/

void
cs_matrix_set_tuning(cs_matrix_fill_type_t   fill_type,
                     int                     tune);

/*----------------------------------------------------------------------------
 * Return matrix tuning behavior for a given fill type.
 *
 * parameters:
 *   fill type  <-- Fill type for which tuning behavior is set
 *
 * returns:
 *   1 if tuning is active, 0 otherwise
 *----------------------------------------------------------------------------*/

int
cs_matrix_get_tuning(cs_matrix_fill_type_t   fill_type);

/*----------------------------------------------------------------------------
 * Set number of matrix computation runs for tuning.
 *
 * If this function is not called, defaults are:
 *  - minimum of 10 runs
 *  - minimum of 0.5 seconds of running
 *
 * parameters:
 *   n_min_products <-- minimum number of expected SpM.V products for
 *                      coefficients assign amortization.
 *   t_measure      <-- minimum running time per measure
 *----------------------------------------------------------------------------*/

void
cs_matrix_set_tuning_runs(int     n_min_products,
                          double  t_measure);

/*----------------------------------------------------------------------------
 * Get number of matrix computation runs for tuning.
 *
 * parameters:
 *   n_min_products --> minimum number of expected SpM.V products for
 *                      coefficients assign amortization.
 *   t_measure      --> minimum running time per measure, or NULL
 *----------------------------------------------------------------------------*/

void
cs_matrix_get_tuning_runs(int     *n_min_products,
                          double  *t_measure);

/*----------------------------------------------------------------------------
 * Return a (0-based) global block row numbering.
 *
 * The numbering is built if not previously present, and returned otherwise.
 *
 * Currently, the function only handles one n_rows/halo combination, and does
 * not check for consistency.
 *
 * parameters:
 *   n_rows <-- associated number of local rows
 *   halo   <-- associated halo, or NULL
 *
 * returns:
 *   pointer to requested global numbering
 *----------------------------------------------------------------------------*/

const cs_gnum_t *
cs_matrix_get_block_row_g_id(cs_lnum_t         n_rows,
                             const cs_halo_t  *halo);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign coefficients to a matrix using a matrix assembler.
 *
 * \param[in]  f                      pointer to associated field
 * \param[in]  type                   matrix type
 * \param[in]  symmetric              is matrix symmetric ?
 * \param[in]  diag_block_size        block sizes for diagonal, or NULL
 * \param[in]  extra_diag_block_size  block sizes for extra diagonal, or NULL
 * \param[in]  da                     diagonal values (NULL if zero)
 * \param[in]  xa                     extradiagonal values (NULL if zero)
 *                                    casts as:
 *                                      xa[n_edges]    if symmetric,
 *                                      xa[n_edges][2] if non symmetric
 *
 * \return  pointer to associated matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_matrix_t *
cs_matrix_set_coefficients_coupled(const cs_field_t  *f,
                                   cs_matrix_type_t   type,
                                   bool               symmetric,
                                   const int         *diag_block_size,
                                   const int         *extra_diag_block_size,
                                   const cs_real_t   *da,
                                   const cs_real_t   *xa);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MATRIX_DEFAULT_H__ */
