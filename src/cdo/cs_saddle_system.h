#ifndef __CS_SADDLE_SYSTEM_H__
#define __CS_SADDLE_SYSTEM_H__

/*============================================================================
 * Operations on cs_cdo_system_t structures: matrix/vector multiplications
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

#include "cs_cdo_system.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the inverse of the diagonal of the (1,1)-block matrix
 *        The storage of a matrix is in a gather view and the resulting array is
 *        in a scatter view.
 *
 * \param[in] b11_max_size  max size related to the (1,1) block
 * \param[in] sh            pointer to a system helper structure
 *
 * \return a pointer to the computed array (scatter view)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_saddle_system_b11_inv_diag(cs_lnum_t                b11_max_size,
                              cs_cdo_system_helper_t  *sh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform a matrix-vector multiplication for the (1,1) block when the
 *        input vector is in a scatter state.

 *        Thus, one performs a scatter --> gather (before the multiplication)
 *        and a gather --> scatter operation after the multiplication. One
 *        assumes that m11x1 is allocated to the right size. No check is done.
 *
 * \param[in]      sh        pointer to a system helper structure
 * \param[in, out] vec       vector
 * \param[in, out] matvec    resulting vector for the matrix-vector product
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_system_b11_matvec(const cs_cdo_system_helper_t  *sh,
                            cs_real_t                     *vec,
                            cs_real_t                     *matvec);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the resulting vector of the operation m12*x2
 *        The block (1,2) is stored in an unassembled way and corresponds to
 *        the (2,1) block. Thus, one considers a transposition of this block.
 *        x2 corresponds to a "scatter" view
 *
 * \param[in]      sh         pointer to a system helper structure
 * \param[in]      x2         array to be multiplied
 * \param[in, out] res        result array storing the matrix.vector operation
 * \param[in]      reset_res  false --> this is an update
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_system_b12_matvec(const cs_cdo_system_helper_t  *sh,
                            const cs_real_t               *x2,
                            cs_real_t                     *res,
                            bool                           reset_res);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the resulting vector of the operation m21*x1
 *        The block (2,1) is stored in an unassembled way.
 *        x1 corresponds to a "scatter" view
 *
 * \param[in]      sh   pointer to a system helper structure
 * \param[in]      x1   array to be multiplied
 * \param[in, out] res  result array storing the matrix.vector operation
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_system_b21_matvec(const cs_cdo_system_helper_t  *sh,
                            const cs_real_t               *x1,
                            cs_real_t                     *res);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the matrix-vector operation for a saddle-point system
 *        r1 = M11.x1 + M12.x2 (result for the first set)
 *        r2 = M21.x1          (result for the second set)
 *
 * \param[in]      sh  pointer to a system helper structure
 * \param[in, out] x1  array for the first part
 * \param[in, out] x2  array for the second part
 * \param[in, out] r1  result array for the first set of DoFs
 * \param[in, out] r2  result array for the second set of DoFs
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_system_matvec(const cs_cdo_system_helper_t  *sh,
                        cs_real_t                     *x1,
                        cs_real_t                     *x2,
                        cs_real_t                     *r1,
                        cs_real_t                     *r2);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the residual of the saddle-point system
 *        res1 = rhs1 - M11.x1 - M12.x2 (residual for the first set)
 *        res2 = rhs2 - M21.x1          (residual for the second set)
 *
 * \param[in]      sh    pointer to a system helper structure
 * \param[in, out] x1    array for the first part
 * \param[in, out] x2    array for the second part
 * \param[in, out] res1  resulting residual vector for the first set of DoFs
 * \param[in, out] res2  resulting residual vector for the second set of DoFs
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_system_residual(const cs_cdo_system_helper_t  *sh,
                          cs_real_t                     *x1,
                          cs_real_t                     *x2,
                          cs_real_t                     *res1,
                          cs_real_t                     *res2);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SADDLE_SYSTEM_H__ */
