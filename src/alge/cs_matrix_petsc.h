#ifndef __CS_MATRIX_PETSC_H__
#define __CS_MATRIX_PETSC_H__

/*============================================================================
 * Sparse Matrix Representation and Operations using PETSc library.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

void
cs_matrix_petsc_ensure_init(void);

/*----------------------------------------------------------------------------
 * Finalize PETSc
 *----------------------------------------------------------------------------*/

void
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

void
cs_matrix_set_type_petsc(cs_matrix_t  *matrix,
                         const char   *type_name);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MATRIX_PETSC_H__ */
