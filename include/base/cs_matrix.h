/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

#ifndef __CS_MATRIX_H__
#define __CS_MATRIX_H__

/*============================================================================
 * Sparse Matrix Representation and Operations
 *============================================================================*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_defs.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

#include "cs_halo.h"
#include "cs_perio.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Matrix types
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_MATRIX_NATIVE,     /* Native matrix format */
  CS_MATRIX_CSR,        /* Compressed Sparse Row storage format */
  CS_MATRIX_N_TYPES     /* Number of known matrix types */

} cs_matrix_type_t;

/* Structure associated with opaque matrix object */

typedef struct _cs_matrix_t cs_matrix_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Short names for matrix types */

extern const char  *cs_matrix_type_name[];

/* Full names for matrix types */

extern const char  *cs_matrix_type_fullname[];

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
 *   type        --> Type of matrix considered
 *   symmetric   --> Indicates if a symmetric variant of the matrix type
 *                   should be used
 *   have_diag   --> Indicates if the diagonal structure contains nonzeroes
 *   periodic    --> Indicates if periodicity is present
 *   n_cells     --> Local number of cells
 *   n_cells_ext --> Local number of cells + ghost cells sharing a face
 *   n_faces     --> Local number of internal faces
 *   cell_num    --> Global cell numbers (1 to n)
 *   face_cell   --> Face -> cells connectivity (1 to n)
 *   halo        --> Halo structure associated with cells, or NULL
 *
 * returns:
 *   pointer to created matrix structure;
 *----------------------------------------------------------------------------*/

cs_matrix_t *
cs_matrix_create(cs_matrix_type_t   type,
                 cs_bool_t          symmetric,
                 cs_bool_t          have_diag,
                 cs_bool_t          periodic,
                 cs_int_t           n_cells,
                 cs_int_t           n_cells_ext,
                 cs_int_t           n_faces,
                 const fvm_gnum_t  *cell_num,
                 const cs_int_t    *face_cell,
                 const cs_halo_t   *halo);

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

cs_int_t
cs_matrix_get_n_columns(const cs_matrix_t  *matrix);

/*----------------------------------------------------------------------------
 * Return number of rows in matrix.
 *
 * parameters:
 *   matrix --> Pointer to matrix structure
 *----------------------------------------------------------------------------*/

cs_int_t
cs_matrix_get_n_rows(const cs_matrix_t  *matrix);

/*----------------------------------------------------------------------------
 * Set matrix coefficients.
 *
 * Depending on current options and initialization, values will be copied
 * or simply mapped.
 *
 * parameters:
 *   matrix    <-> Pointer to matrix structure
 *   symmetric --> Indicates if matrix coefficients are symmetric
 *   da        --> Diagonal values (NULL if zero)
 *   xa        --> Extradiagonal values (NULL if zero)
 *----------------------------------------------------------------------------*/

void
cs_matrix_set_coefficients(cs_matrix_t      *matrix,
                           cs_bool_t         symmetric,
                           const cs_real_t  *da,
                           const cs_real_t  *xa);

/*----------------------------------------------------------------------------
 * Release matrix coefficients.
 *
 * parameters:
 *   matrix <-> Pointer to matrix structure
 *----------------------------------------------------------------------------*/

void
cs_matrix_release_coefficients(cs_matrix_t  *matrix);

/*----------------------------------------------------------------------------
 * Get matrix diagonal values.
 *
 * parameters:
 *   matrix --> Pointer to matrix structure
 *   da     <-- Diagonal (pre-allocated, size: n_cells)
 *----------------------------------------------------------------------------*/

void
cs_matrix_get_diagonal(const cs_matrix_t  *matrix,
                       cs_real_t          *restrict da);

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x
 *
 * This function includes a halo update of x prior to multiplication by A.
 *
 * parameters:
 *   rotation_mode --> Halo update option for rotational periodicity
 *   matrix        --> Pointer to matrix structure
 *   x             <-> Multipliying vector values (ghost values updated)
 *   y             <-- Resulting vector
 *----------------------------------------------------------------------------*/

void
cs_matrix_vector_multiply(cs_perio_rota_t     rotation_mode,
                          const cs_matrix_t  *matrix,
                          cs_real_t          *restrict x,
                          cs_real_t          *restrict y);

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
 *   y      <-- Resulting vector
 *----------------------------------------------------------------------------*/

void
cs_matrix_vector_multiply_nosync(const cs_matrix_t  *matrix,
                                 const cs_real_t    *x,
                                 cs_real_t          *restrict y);

/*----------------------------------------------------------------------------
 * Matrix.vector product y = alpha.A.x + beta.y
 *
 * This function includes a halo update of x prior to multiplication by A.
 *
 * parameters:
 *   rotation_mode --> Halo update option for rotational periodicity
 *   alpha         --> Scalar, alpha in alpha.A.x + beta.y
 *   beta          --> Scalar, beta in alpha.A.x + beta.y
 *   matrix        --> Pointer to matrix structure
 *   x             <-> Multipliying vector values (ghost values updated)
 *   y             <-- Resulting vector
 *----------------------------------------------------------------------------*/

void
cs_matrix_alpha_a_x_p_beta_y(cs_perio_rota_t     rotation_mode,
                             cs_real_t           alpha,
                             cs_real_t           beta,
                             const cs_matrix_t  *matrix,
                             cs_real_t          *restrict x,
                             cs_real_t          *restrict y);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_MATRIX_H__ */
