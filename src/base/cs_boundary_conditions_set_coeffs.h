#ifndef __CS_BOUNDARY_CONDITIONS_SET_COEFFS_H__
#define __CS_BOUNDARY_CONDITIONS_SET_COEFFS_H__

/*============================================================================
 * Boundary condition management.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \file cs_boundary_conditions_set_coeffs.c
 *
 * \brief Translation of the boundary conditions given by the user in a form
 *        that fits to the solver.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Translation of the boundary conditions given by the user
 * in a form that fits to the solver.
 *
 * The values at a boundary face \f$ \fib \f$ stored in the face center
 * \f$ \centf \f$ of the variable \f$ P \f$ and its diffusive flux \f$ Q \f$
 * are written as:
 * \f[
 * P_{\face} = A_P^g + B_P^g P_{\centi}
 * \f]
 * and
 * \f[
 * Q_{\face} = A_P^f + B_P^f P_{\centi}
 * \f]
 * where \f$ P_\centi \f$ is the value of the variable \f$ P \f$ at the
 * neighboring cell.
 *
 * \warning
 * - If we consider an increment of a variable, the boundary conditions
 *   read:
 *   \f[
 *   \delta P_{\face} = B_P^g \delta P_{\centi}
 *   \f]
 *   and
 *   \f[
 *   \delta Q_{\face} = -B_P^f \delta P_{\centi}
 *   \f]
 *
 * - For a vector field such as the velocity \f$ \vect{u} \f$ the boundary
 *   conditions may read:
 *   \f[
 *   \vect{u}_{\face} = \vect{A}_u^g + \tens{B}_u^g \vect{u}_{\centi}
 *   \f]
 *   and
 *   \f[
 *   \vect{Q}_{\face} = \vect{A}_u^f + \tens{B}_u^f \vect{u}_{\centi}
 *   \f]
 *   where \f$ \tens{B}_u^g \f$ and \f$ \tens{B}_u^f \f$ are 3x3 tensor matrix
 *   which coupled velocity components next to a boundary.
 *
 * Please refer to the
 * <a href="../../theory.pdf#boundary"><b>boundary conditions</b></a> section
 * of the theory guide for more informations, as well as the
 * <a href="../../theory.pdf#condli"><b>condli</b></a> section.
 *
 * \param[in]     nvar          total number of variables
 * \param[in]     iterns        iteration number on Navier-Stokes equations
 * \param[in]     isvhb         indicator to save exchange coeffient
 *                               at the walls
 * \param[in]     itrale        ALE iteration number
 * \param[in]     italim        for ALE
 * \param[in]     itrfin        for ALE
 * \param[in]     ineefl        for ALE
 * \param[in]     itrfup        for ALE
 * \param[in,out] icodcl        face boundary condition code:
 *                               - 1 Dirichlet
 *                               - 2 Radiative outlet
 *                               - 3 Neumann
 *                               - 4 sliding and
 *                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
 *                               - 5 smooth wall and
 *                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
 *                               - 6 rough wall and
 *                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
 *                               - 9 free inlet/outlet
 *                                 (input mass flux blocked to 0)
 *                               - 10 Boundary value related to the next cell
 *                                 value by an affine function
 *                               - 11 Generalized Dirichlet for vectors
 *                               - 12 Dirichlet boundary value related to the
 *                                 next cell value by an affine function for
 *                                 the advection operator and Neumann for the
 *                                 diffusion operator
 *                               - 13 Dirichlet for the advection operator and
 *                                 Neumann for the diffusion operator
 *                               - 14 Generalized symmetry for vectors (used for
 *                                 Marangoni effects modeling)
 *                               - 15 Neumann for the advection operator and
 *                                 homogeneous Neumann for the diffusion
 *                                 operator (walls with hydro. pressure for
 *                                 the compressible module)
 * \param[in,out] isostd        indicator for standard outlet
 *                              and reference face index
 * \param[in]     dt            time step (per cell)
 * \param[in,out] rcodcl        boundary condition values:
 *                               - rcodcl(1) value of the Dirichlet
 *                               - rcodcl(2) value of the exterior exchange
 *                                 coefficient (infinite if no exchange)
 *                               - rcodcl(3) value flux density
 *                                 (negative if gain) in w/m2 or roughness
 *                                 in m if icodcl=6
 *                                 -# for the velocity \f$ (\mu+\mu_T)
 *                                    \gradv \vect{u} \cdot \vect{n}  \f$
 *                                 -# for the pressure \f$ \Delta t
 *                                    \grad P \cdot \vect{n}  \f$
 *                                 -# for a scalar \f$ cp \left( K +
 *                                     \dfrac{K_T}{\sigma_T} \right)
 *                                     \grad T \cdot \vect{n} \f$
 * \param[out]    visvdr        dynamic viscosity after V. Driest damping in
 *                              boundary cells
 * \param[out]    hbord         exchange coefficient at boundary
 * \param[out]    theipb        value of thermal scalar at \f$ \centip \f$
 *                              of boundary cells
 * \param[in]     nftcdt        Global indicator of condensation source terms
 *                              (ie. sum on the processors of nfbpcd) cells
 *                              associated to the face with condensation
 *                              phenomenon
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_coeffs(int        nvar,
                                  int        iterns,
                                  int        isvhb,
                                  int        itrale,
                                  int        italim,
                                  int        itrfin,
                                  int        ineefl,
                                  int        itrfup,
                                  int        isostd[],
                                  cs_real_t  visvdr[],
                                  cs_real_t  hbord[],
                                  cs_real_t  theipb[],
                                  int        nftcdt);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialization of boundary condition arrays.
 *
 * \param[in]     itrale        ALE iteration number
 * \param[in,out] isostd        indicator for standard outlet
 *                              and reference face index
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_coeffs_init(int  itrale,
                                       int  isostd[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BOUNDARY_CONDITIONS_SET_COEFFS_H__ */
