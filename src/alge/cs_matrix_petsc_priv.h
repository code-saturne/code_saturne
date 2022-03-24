#ifndef __CS_MATRIX_PETSC_PRIV_H__
#define __CS_MATRIX_PETSC_PRIV_H__

/*============================================================================
 * Private types for sparse matrix representation and operations using PETSc.
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
 * PETSc headers
 *----------------------------------------------------------------------------*/

#include <petscvec.h>
#include <petscmat.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "cs_matrix_petsc.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Note that most types are declared in cs_matrix_priv.h.
   only those only handled here are declared here. */

/* Adapter coefficients stucture for PETSc */

typedef struct _cs_matrix_coeffs_petsc_t {

  PetscInt   l_range[2];                   /* Local rows range
                                              (block range if block size > 1) */

  Mat  hm;                                 /* PETSc matrix */
  Vec  hx;                                 /* x (input) vector */
  Vec  hy;                                 /* y (output) vector */

  int  matrix_state;                       /* Matrix state:
                                              0: not created
                                              1: created and assembled */

  char     *matype_r;                      /* requested PETSc matrix type */
  MatType   matype;                        /* actual PETSc matrix type */

} cs_matrix_coeffs_petsc_t;

/*=============================================================================
 * Semi-private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief return coefficients structure associated with PETSc matrix.
 *
 * \param[in]  matrix  pointer to matrix structure
 *
 * \return  pointer to matrix coefficients handler structure for PETSc matrix.
 */
/*----------------------------------------------------------------------------*/

cs_matrix_coeffs_petsc_t *
cs_matrix_petsc_get_coeffs(const cs_matrix_t  *matrix);

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MATRIX_PETSC_PRIV_H__ */
