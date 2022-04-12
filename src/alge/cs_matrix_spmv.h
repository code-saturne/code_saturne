#ifndef __CS_MATRIX_SPMV_H__
#define __CS_MATRIX_SPMV_H__

/*============================================================================
 * Sparse Matrix-vector multiplication functions and kernels.
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

#include "cs_halo.h"
#include "cs_numbering.h"
#include "cs_matrix_priv.h"
#include "cs_matrix_assembler.h"

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
 * Semi-private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Assign default sparse matrix-vector product functions for a given matrix.
 *
 * parameters:
 *   m <-> pointer to matrix structure
 *----------------------------------------------------------------------------*/

void
cs_matrix_spmv_set_defaults(cs_matrix_t  *m);

/*----------------------------------------------------------------------------
 * Select the sparse matrix-vector product function to be used by a
 * matrix or variant for a given fill type.
 *
 * Currently, possible variant functions are:
 *
 *   CS_MATRIX_NATIVE  (all fill types)
 *     default
 *     baseline
 *     omp             (for OpenMP with compatible numbering)
 *     omp_atomic      (for OpenMP with atomic add)
 *     vector          (For vector machine with compatible numbering)
 *
 *   CS_MATRIX_CSR     (for CS_MATRIX_SCALAR or CS_MATRIX_SCALAR_SYM)
 *     default
 *     mkl             (with MKL)
 *     cuda            (CUDA-accelerated)
 *     cusparse        (with cuSPARSE)
 *
 *   CS_MATRIX_MSR
 *     default
 *     omp_sched       (Improved OpenMP scheduling, for CS_MATRIX_SCALAR*)
 *     mkl             (with MKL, for CS_MATRIX_SCALAR or CS_MATRIX_SCALAR_SYM)
 *     cuda            (CUDA-accelerated)
 *     cusparse        (with cuSPARSE)
 *
 *   CS_MATRIX_DIST
 *     default
 *     omp_sched       (Improved OpenMP scheduling, for CS_MATRIX_SCALAR*)
 *     mkl             (with MKL, for CS_MATRIX_SCALAR or CS_MATRIX_SCALAR_SYM)
 *
 * parameters:
 *   m_type      <--  Matrix type
 *   fill type   <--  matrix fill type to merge from
 *   spmv_type   <--  SpMV operation type (full or sub-matrix)
 *                    (all types if CS_MATRIX_SPMV_N_TYPES)
 *   numbering   <--  mesh numbering structure, or NULL
 *   func_name   <--  function type name, or NULL for default
 *   spmv        <->  multiplication function array
 *   spmv_xy_hd  <->  multiplication function x and y host/device location
 *
 * returns:
 *   0 for success, 1 for incompatible function, 2 for compatible
 *   function not available in current build
 *----------------------------------------------------------------------------*/

int
cs_matrix_spmv_set_func(cs_matrix_type_t             m_type,
                        cs_matrix_fill_type_t        fill_type,
                        cs_matrix_spmv_type_t        spmv_type,
                        const cs_numbering_t        *numbering,
                        const char                  *func_name,
                        cs_matrix_vector_product_t  *spmv[CS_MATRIX_SPMV_N_TYPES],
                        char                   spmv_xy_hd[CS_MATRIX_SPMV_N_TYPES]);

/*======================================à=======================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MATRIX_SPMV_H__ */
