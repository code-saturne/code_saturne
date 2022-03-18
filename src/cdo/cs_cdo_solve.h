#ifndef __CS_CDO_SOLVE_H__
#define __CS_CDO_SOLVE_H__

/*============================================================================
 * Set of helper functions to prepare or solve linear systems
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_matrix.h"
#include "cs_param_sles.h"
#include "cs_range_set.h"
#include "cs_sles.h"

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
 * \brief  Compute the value of the rhs norm used as a renormalization
 *         coefficient for the residual norm when solving the linear system. A
 *         pre-computation performed during the cellwise building of the
 *         algebraic has been done before this call according to the requested
 *         type of renormalization.
 *
 * \param[in]      type            type of renormalization
 * \param[in]      vol_tot         total volume of the computational domain
 * \param[in]      rhs_size        size of the rhs array
 * \param[in]      rhs             array related to the right-hand side
 * \param[in, out] normalization   value of the residual normalization
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_solve_sync_rhs_norm(cs_param_resnorm_type_t    type,
                           double                     vol_tot,
                           cs_lnum_t                  rhs_size,
                           const cs_real_t            rhs[],
                           double                    *normalization);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prepare a linear system and synchronize buffers in case of parallel
 *        or periodic computations. Transfer a mesh-based description of
 *        arrays x0 and rhs into an algebraic description for the linear system
 *        in x and b (scatter/gather views).
 *
 * \param[in]      stride     stride to apply to the range set operations
 * \param[in]      interlace  is data interlaced or not
 * \param[in]      x_size     size of the vector unknowns (scatter view)
 * \param[in]      rset       pointer to a range set structure
 * \param[in]      rhs_redux  do or not a parallel sum reduction on the RHS
 * \param[in, out] x          array of unknowns (in: initial guess)
 * \param[in, out] b          right-hand side
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_solve_prepare_system(int                     stride,
                            bool                    interlace,
                            cs_lnum_t               x_size,
                            const cs_range_set_t   *rset,
                            bool                    rhs_redux,
                            cs_real_t              *x,
                            cs_real_t              *b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a linear system arising with scalar-valued cell-based DoFs
 *         No rotation is taken into account when synchronizing halo.
 *
 * \param[in]  n_dofs         local number of DoFs
 * \param[in]  slesp          pointer to a cs_param_sles_t structure
 * \param[in]  matrix         pointer to a cs_matrix_t structure
 * \param[in]  normalization  value used for the residual normalization
 * \param[in, out] sles       pointer to a cs_sles_t structure
 * \param[in, out] x          solution of the linear system (in: initial guess)
 * \param[in, out] b          right-hand side (scatter/gather if needed)
 *
 * \return the number of iterations of the linear solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdo_solve_scalar_cell_system(cs_lnum_t                n_dofs,
                                const cs_param_sles_t   *slesp,
                                const cs_matrix_t       *matrix,
                                cs_real_t                normalization,
                                cs_sles_t               *sles,
                                cs_real_t               *x,
                                cs_real_t               *b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a linear system arising from CDO schemes with scalar-valued
 *         degrees of freedom
 *
 * \param[in]  n_scatter_dofs local number of DoFs (may be != n_gather_elts)
 * \param[in]  slesp          pointer to a cs_param_sles_t structure
 * \param[in]  matrix         pointer to a cs_matrix_t structure
 * \param[in]  rs             pointer to a cs_range_set_t structure
 * \param[in]  normalization  value used for the residual normalization
 * \param[in]  rhs_redux      do or not a parallel sum reduction on the RHS
 * \param[in, out] sles       pointer to a cs_sles_t structure
 * \param[in, out] x          solution of the linear system (in: initial guess)
 * \param[in, out] b          right-hand side (scatter/gather if needed)
 *
 * \return the number of iterations of the linear solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdo_solve_scalar_system(cs_lnum_t                     n_scatter_dofs,
                           const cs_param_sles_t        *slesp,
                           const cs_matrix_t            *matrix,
                           const cs_range_set_t         *rset,
                           cs_real_t                     normalization,
                           bool                          rhs_redux,
                           cs_sles_t                    *sles,
                           cs_real_t                    *x,
                           cs_real_t                    *b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a linear system arising from CDO schemes with vector-valued
 *         degrees of freedom (DoFs).
 *         Number of DoFs is equal to 3*n_scatter_elts
 *
 * \param[in]  n_scatter_elts local number of elements (may be != n_gather_elts)
 * \param[in]  interlace      way to arrange data (true/false)
 * \param[in]  slesp          pointer to a cs_param_sles_t structure
 * \param[in]  matrix         pointer to a cs_matrix_t structure
 * \param[in]  rset           pointer to a cs_range_set_t structure
 * \param[in]  normalization  value used for the residual normalization
 * \param[in]  rhs_redux      do or not a parallel sum reduction on the RHS
 * \param[in, out] sles       pointer to a cs_sles_t structure
 * \param[in, out] x          solution of the linear system (in: initial guess)
 * \param[in, out] b          right-hand side (scatter/gather if needed)
 *
 * \return the number of iterations of the linear solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdo_solve_vector_system(cs_lnum_t                     n_scatter_elts,
                           bool                          interlace,
                           const cs_param_sles_t        *slesp,
                           const cs_matrix_t            *matrix,
                           const cs_range_set_t         *rset,
                           cs_real_t                     normalization,
                           bool                          rhs_redux,
                           cs_sles_t                    *sles,
                           cs_real_t                    *x,
                           cs_real_t                    *b);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_SOLVE_H__ */
