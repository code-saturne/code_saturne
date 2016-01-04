/*============================================================================
 * Set advanced numerical parameters for the current simulation
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

#include "cs_param_eq.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_user_cdo_numerics.c

  \brief Set advanced parameters about the numerical schemes for each
         equation to solve.
         Useful to change the default behaviour.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup advanced features concerning the way geometric quantities
 *         are built
 *
 * \return the type of computation to evaluate the cell center
 */
/*----------------------------------------------------------------------------*/

cs_cdo_cc_algo_t
cs_user_cdo_geometric_settings(void)
{
  /* Algorithm for computing cell centers */
  /* ==================================== */

  /* Choice between:
     CS_CDO_CC_MEANV: Cell center is computed as the mean of cell vertices
     CS_CDO_CC_BARYC: Cell center is computed as the real cell barycenter
     CS_CDO_CC_SATUR: Cell center is given by Code_Saturne
     CS_CDO_CC_ORTHO: Cell center is optimized to enforce orthogonality
                      between cell-face edge and face plane
   */

  return CS_CDO_CC_BARYC;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup advanced features concerning the numerical parameters
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_numeric_settings(void)
{
  /* Set the scheme used to discretize in space
     --> Choice between:
     CS_SPACE_SCHEME_CDOVB: CDO vertex-based scheme
     CS_SPACE_SCHEME_CDOFB: CDO cell-based scheme with hybridization
                            Degrees of freedom are located on faces
  */

  cs_param_eq_set_space_scheme("Laplace", // Equation name
                               CS_SPACE_SCHEME_CDOVB);

  /* Verbosity level (former IWARNI)
     --> Choice between:
     0  >> default
     1  >> detailed setup resume and coarse grain timer stats
     2  >> fine grain timer stats
     >2 >> more detailed (for debugging purpose)
   */
  cs_param_eq_set_verbosity_level("Laplace", 2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup advanced features concerning the discrete Hodge operators
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_hodge_settings(void)
{
  /* Set the algorithm used to build the discrete Hodge operator
     --> Choice between:
     CS_PARAM_HODGE_ALGO_VORONOI (only possible in specific cases: Cartesian
                                  meshes or Delaunay meshes for instance)
     CS_PARAM_HODGE_ALGO_COST    (default: splitting COnsistency/STabilization)
  */

  cs_param_eq_hodge_diffusion_set_algo("Laplace",  // Equation name
                                       CS_PARAM_HODGE_ALGO_COST); // Algo

  /* In the case of a COST algorithm, you may additionally specify the value
     of the stabilization coefficient (> 0)
     coef = 1./3       --> DGA scheme (Default)
     coef = 1./sqrt(3) --> SUSHI scheme
     coef = 1.         --> Generalized Crouzeix--Raviart scheme
     Other choices are possible leading to different schemes...
  */

  cs_param_eq_hodge_diffusion_set_coef("Laplace", 1./3.);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup advanced features concerning the resolution of linear systems
 *         within the CDO framework
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_itsol_settings(void)
{
  const char eqname[16] = "Laplace";

  /* Choice of the algorithm used to solve an equation among
     >> CS_PARAM_EQ_ALGO_CS_ITSOL     iter. solvers implemented in Code_Saturne
     >> CS_PARAM_EQ_ALGO_PETSC_ITSOL  iter. solvers implemented in PETSc
     >> ...

     Default is CS_PARAM_EQ_ALGO_CS_ITSOL
   */

  cs_param_eq_set_algo_type(eqname,
                            CS_PARAM_EQ_ALGO_CS_ITSOL);

  /* An iterative solver is defined thanks to
     >> a solver type among
         >> CS_PARAM_ITSOL_CG       conjuguate gradient
         >> CS_PARAM_ITSOL_BICG     bi-conjuguate gradient
         >> CS_PARAM_ITSOL_GMRES    GMRES algorithm

     >> a preconditioner type among
         >> CS_PARAM_PRECOND_DIAG   Jacobi (inverse of the diagonal)
         >> CS_PARAM_PRECOND_POLY   Neumann polynomial of order 1
         >> CS_PARAM_PRECOND_SSOR   Symmetric Successive Over-Relaxation
         >> CS_PARAM_PRECOND_ILU0   Incomplete LU factorisation (type 0)
         >> CS_PARAM_PRECOND_MG     Algebraic/Geometric Multigrid

     Default is CS_PARAM_ITSOL_CG / CS_PARAM_PRECOND_SSOR
  */

  cs_param_eq_set_itsol_type(eqname,
                             CS_PARAM_ITSOL_CG,
                             CS_PARAM_PRECOND_DIAG);

  /* Additionnal settings related to an iterative solver
     >> The required accuracy on the norm of the residual
     >> The max number of iterations of the iterative solver
     >> The normalization or not of the norm of residual (impact on the
         overall accuracy of the algorithm)
  */

  cs_param_eq_set_itsol_precision(eqname,  // equation name
                                  1e-12);  // precision

  cs_param_eq_set_itsol_max_iter(eqname,  // equation name
                                 2500);   // max. number of iterations

  cs_param_eq_set_itsol_normalization(eqname,  // equation name
                                      false);  // normalize residual ?
}
