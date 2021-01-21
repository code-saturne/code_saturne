#ifndef __CS_SADDLE_ITSOL_H__
#define __CS_SADDLE_ITSOL_H__

/*============================================================================
 * In-house iterative solvers defined by blocks and associated to CDO
 * discretizations for saddle-point system
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_mesh_adjacencies.h"
#include "cs_iter_algo.h"
#include "cs_matrix.h"
#include "cs_range_set.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct {

  /*
   * A hybrid storage is used to represent the system described below
   *
   *      | M11 | M12 | | x1 |   |rhs1 |
   *  M = |-----------| |----| = |-----|
   *      | M21 |  0  | | x2 |   |rhs2 |
   *
   *  One assumes that M12 = M21^T
   */

  int              n_m11_matrices; /* 1, 3, 6 or 9 is expected */
  cs_matrix_t    **m11_matrices;
  cs_lnum_t        x1_size;        /* scatter view */
  cs_lnum_t        max_x1_size;    /* max(x1_size, n_m11_cols) */

  cs_real_t       *rhs1;

  cs_lnum_t        x2_size;        /* scatter view */
  cs_real_t       *rhs2;

  cs_real_t       *m21_unassembled;
  int              m21_stride;

  /* Indexed list used to scan the unassembled m21 operator */
  cs_adjacency_t  *m21_adjacency;

  /* Structure used for synchronisation (parallel or periodic). Enable to
     switch from a scatter view (the mesh view) to a gather view (the algebraic
     view) */
  cs_range_set_t  *rset;

} cs_saddle_system_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform a matrix-vector multiplication in case of scatter-view array
 *        as input parameter.  Thus, one performs a scatter --> gather (before
 *        the multiplication) and a gather --> scatter operation after the
 *        multiplication.  One assumes that matvec is allocated to the right
 *        size. No check is done.
 *
 *        The stride is equal to 1 for the matrix (db_size[3] = 1) and the
 *        vector
 *
 * \param[in]      rset      pointer to a cs_range_set_t structure
 * \param[in]      mat       matrix
 * \param[in, out] vec       vector
 * \param[in, out] matvec    resulting vector for the matrix-vector product
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_vector_multiply_gs_allocated(const cs_range_set_t      *rset,
                                       const cs_matrix_t         *mat,
                                       cs_real_t                 *vec,
                                       cs_real_t                 *matvec);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform a matrix-vector multiplication in case of scatter-view array
 *        as input parameter.  Thus, one performs a scatter --> gather (before
 *        the multiplication) and a gather --> scatter operation after the
 *        multiplication.  The output parameter matvec is not allocated. A
 *        check on the size is done for the input array.
 *
 *        The stride is equal to 1 for the matrix (db_size[3] = 1) and the
 *        vector
 *
 * \param[in]      rset      pointer to a cs_range_set_t structure
 * \param[in]      mat       matrix
 * \param[in]      vec_len   size of vec
 * \param[in, out] vec       vector of real numbers
 * \param[out]     matvec    resulting vector for the matrix-vector product
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_vector_multiply_gs(const cs_range_set_t      *rset,
                             const cs_matrix_t         *mat,
                             cs_lnum_t                  vec_len,
                             cs_real_t                 *vec,
                             cs_real_t                **p_matvec);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply the MINRES algorithm to a saddle point problem (the system is
 *        stored in a hybrid way). Please refer to cs_saddle_system_t structure
 *        definition.
 *        The stride is equal to 1 for the matrix (db_size[3] = 1) and the
 *        vector
 *
 * \param[in]      ssys      pointer to a cs_saddle_system_t structure
 * \param[in, out] x1        array for the first part
 * \param[in, out] x2        array for the second part
 * \param[in, out] info      pointer to a cs_iter_algo_info_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_minres(cs_saddle_system_t    *ssys,
                 cs_real_t             *x1,
                 cs_real_t             *x2,
                 cs_iter_algo_info_t   *info);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform elementary tests to assess this module
 *
 * \param[in]      ssys      pointer to a cs_saddle_system_t structure
 * \param[in, out] x1        array for the first part
 * \param[in, out] x2        array for the second part
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_test(cs_saddle_system_t   *ssys,
               cs_real_t            *x1,
               cs_real_t            *x2);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SADDLE_ITSOL_H__ */
