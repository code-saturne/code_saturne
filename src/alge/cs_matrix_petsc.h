#ifndef CS_MATRIX_PETSC_H
#define CS_MATRIX_PETSC_H

/*============================================================================
 * Sparse Matrix Representation and Operations using PETSc library.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include "base/cs_defs.h"

#include "alge/cs_matrix.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize PETSc if needed
 *----------------------------------------------------------------------------*/

extern "C" void
cs_matrix_petsc_ensure_init(void);

/*----------------------------------------------------------------------------
 * Finalize PETSc
 *----------------------------------------------------------------------------*/

extern "C" void
cs_matrix_petsc_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Switch matrix type to PETSc.
 *
 * This releases previous coefficients if present, so should be called
 * just after matrix creation, before assigning coefficients.
 *
 * \param[in, out]  matrix      pointer to matrix structure
 * \param[in]       type_name   string matching PETSc matrix type name,
 *                              defaults to "MATAIJ" if NULL
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_matrix_set_type_petsc(cs_matrix_t  *matrix,
                         const char   *type_name);

/*----------------------------------------------------------------------------*/
/*
 * \brief Function pointer for adding one row to a petsc matrix.
 *
 * \warning The matrix pointer must point to valid data when the selection
 *          function is called, so the life cycle of the data pointed to should
 *          be at least as long as that of the assembler values structure.
 *
 * \param[in]      n_cols      number of column values to add
 * \param[in]      row_g_id    row global id
 * \param[in]      col_g_ids   column global id list with n_cols elements
 * \param[in]      vals        value list with n_cols elements to add
 * \param[in, out] matrix_p    untyped pointer to matrix description structure
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_matrix_petsc_add_scal_row_values(const cs_lnum_t                n_cols,
                                    const cs_gnum_t                row_g_id,
                                    const cs_gnum_t               *col_g_ids,
                                    const cs_real_t               *vals,
                                    void                          *matrix_p);

/*----------------------------------------------------------------------------*/

#endif /* CS_MATRIX_PETSC_H */
