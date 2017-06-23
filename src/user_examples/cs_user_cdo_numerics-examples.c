/*============================================================================
 * Set advanced numerical parameters for the current simulation when the CDO
 * kernel is used
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include "cs_equation.h"
#include "cs_domain.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_cdo_numerics-examples.c
 *
 * \brief Set advanced parameters about the numerical schemes for each
 *        equation to solve.
 *        Useful to change the default behaviour.
 */
/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

/*============================================================================
 * Private user function definitions
 *============================================================================*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup advanced features concerning the way geometric quantities
 *         are built
 *
 * \return the type of computation to evaluate the cell center
 */
/*----------------------------------------------------------------------------*/

cs_cdo_cell_center_algo_t
cs_user_cdo_geometric_settings(void)
{
  /* Algorithm for computing cell centers */
  /* ==================================== */

  /* Choice between:
     CS_CDO_CCENTER_MEANV: Cell center is computed as the mean of cell vertices
     CS_CDO_CCENTER_BARYC: Cell center is computed as the real cell barycenter
     CS_CDO_CCENTER_SATURNE: Cell center is given by Code_Saturne
   */

  return CS_CDO_CCENTER_BARYC;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup advanced features concerning the numerical parameters
 *         of the equation resolved during the computation
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_numeric_settings(void)
{
  /* Modify the setting of an equation using a generic process

     ***********  cs_equation_set_param(eq, key, "val")   ************

     CS_EQKEY_SPACE_SCHEME
     "cdo_vb"  for CDO vertex-based scheme
     "cdo_vcb" for CDO vertex+cell-based scheme
     "cdo_fb"  for CDO face-based scheme

     CS_EQKEY_VERBOSITY
     The higher the more detailed information is displayed
     "0" (default)
     "1" detailed setup resume and coarse grain timer stats
     "2" fine grain for timer stats

     CS_EQKEY_HODGE_DIFF_ALGO
     CS_EQKEY_HODGE_TIME_ALGO
     CS_EQKEY_HODGE_REAC_ALGO
     "voronoi" (default for time), leads to diagonal discrete Hodge operator
     but is not consistent for all meshes
     "cost" (default for diffusion) is more robust (i.e. it handles more
     general meshes but is is less efficient)
     "wbs" (default for reaction) is robust and accurate but is limited to
     the reconstruction of potential-like degrees of freedom and needs a correct
     computation of the cell barycenter

     CS_EQKEY_HODGE_DIFF_COEF
     CS_EQKEY_HODGE_TIME_COEF
     CS_EQKEY_HODGE_REAC_COEF
     This key is only useful if CS_EQKEY_HODGE_*_ALGO is set to "cost"
     val is either a name or a value: "dga", "sushi", "gcr" or "1.5", "9"..
     "dga" corresponds to the value 1./3.
     "sushi" corresponds to the value 1./sqrt(3.)
     "gcr" corresponds to the value 1.

     CS_EQKEY_SOLVER_FAMILY
     >> val: "cs" (default), "petsc"
     WARNING: For using "petsc" one needs to install Code_Saturne with PETSc

     CS_EQKEY_ITSOL
     >> val: "cg" is among the following choices:
     "cg" (default) is the standard conjuguate gradient algorithm
     "bicg" is Bi-CG algorithm (for non-symmetric linear systems)
     "bicgstab2" is BiCG-Stab2 algorithm (for non-symmetric linear systems)
     "cr3" is a 3-layer conjugate residual solver
     "gmres" is a robust iterative solver but not as efficient
     "amg" is an algebraic multigrid iterative solver

     CS_EQKEY_PRECOND
     >> val is among the following choices:
     "jacobi" diagonal preconditoner
     "block_jacobi"
     "poly1"  neumann polynomial of order 1
     "ssor"   symmetric successive over-relaxation (only with PETSC)
     "ilu0"   incomplete LU factorization
     "icc0"   incomplete Cholesky factorization (for symmetric matrices)
     "amg"    algebraic multigrid

     CS_EQKEY_ITSOL_EPS
     "1e-10" for instance

     CS_EQKEY_ITSOL_MAX_ITER
     "2000" for instance

     CS_EQKEY_ITSOL_RESNORM
     "true" or "false"

     CS_EQKEY_SLES_VERBOSITY
     "0", "1", "2" or higher

     CS_EQKEY_BC_ENFORCEMENT
     Set the type of enforcement of the boundary conditions
     "strong"       remove unknowns attached to a BC
     "penalization" weak enforcement using a huge penalization coefficient
     "weak"         weak enforcement using the Nitsche method
     "weak_sym"     weak enforcement keeping the symmetry of the system

     CS_EQKEY_BC_QUADRATURE
     Set the quadrature algorithm used for evaluating boundary conditions
     "bary"    used the barycenter approximation
     "higher"  used 4 Gauss points for approximating the integral
     "highest" used 5 Gauss points for approximating the integral
     Remark: "higher" and "highest" implies automatically a subdivision into
     tetrahedra of each cell

     CS_EQKEY_TIME_SCHEME
     "implicit": first-order in time (inconditionnally stable)
     "explicit":
     "crank_nicolson": second_order in time
     "theta_scheme": generic time scheme. One recovers "implicit" with theta
     equal to "1", "explicit" with "0", "crank_nicolson" with "0.5"

     CS_EQKEY_TIME_THETA
     Only useful if CS_EQKEY_TIME_SCHEME is set to "theta_scheme"
     >> val: "0.75" for instance (must be between 0 <=val<= 1)

     CS_EQKEY_ADV_FORMULATION
     "conservative"
     "non_conservative"

     CS_EQKEY_ADV_SCHEME
     "upwind"
     "centered"
     "samarskii" upwind/centered with a weight depending on the Peclet number
     "sg"        upwind/centered with a weight depending on the Peclet number
     "cip"       "continuous interior penalty" (only for VCB schemes

     CS_EQKEY_ADV_FLUX_QUADRA (see CS_EQKEY_BC_QUADRATURE)
     "bary" (default)
     "higher"
     "highest"

     CS_EQKEY_EXTRA_OP: Additional post-processing options:
     "peclet"  post-process an estimation of the Peclet number in each cell
     "upwind_coef" post-process an estimation of the upwinding coefficient

  */

  cs_equation_t  *eq = cs_equation_by_name("FVCA6.1");

  /* The modification of the space discretization should be apply first */
  cs_equation_set_param(eq, CS_EQKEY_SPACE_SCHEME, "cdo_vb");

  /* Modify other parameters than the space discretization */
  cs_equation_set_param(eq, CS_EQKEY_VERBOSITY, "2");
  cs_equation_set_param(eq, CS_EQKEY_HODGE_DIFF_ALGO, "cost");
  cs_equation_set_param(eq, CS_EQKEY_HODGE_DIFF_COEF, "dga");

  /* Linear algebra settings */
#if defined(HAVE_PETSC)
  cs_equation_set_param(eq, CS_EQKEY_SOLVER_FAMILY, "petsc");
  cs_equation_set_param(eq, CS_EQKEY_ITSOL, "cg");
  cs_equation_set_param(eq, CS_EQKEY_PRECOND, "amg");
#else
  cs_equation_set_param(eq, CS_EQKEY_SOLVER_FAMILY, "cs");
  cs_equation_set_param(eq, CS_EQKEY_PRECOND, "jacobi");
  cs_equation_set_param(eq, CS_EQKEY_ITSOL, "cg");
#endif
  cs_equation_set_param(eq, CS_EQKEY_ITSOL_MAX_ITER, "2500");
  cs_equation_set_param(eq, CS_EQKEY_ITSOL_EPS, "1e-12");
  cs_equation_set_param(eq, CS_EQKEY_ITSOL_RESNORM, "false");
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
