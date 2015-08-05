/*============================================================================
 * Set advanced numerical parameters for the current simulation
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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

cs_cdo_cc_algo_t
cs_user_cdo_geometric_settings(void)
{
  /* Algorithm for computing cell centers */
  /* ==================================== */

  /* Choice between:
     CS_CDO_CC_MEANV:   Cell center is computed as the mean of cell vertices
     CS_CDO_CC_BARYC:   Cell center is computed as the real cell barycenter
     CS_CDO_CC_SATURNE: Cell center is given by Code_Saturne
     CS_CDO_CC_ORTHO:   Cell center is optimized to enforce orthogonality
                        between cell-face edge and face plane
   */

  return CS_CDO_CC_BARYC;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup advanced features concerning the numerical parameters
 *         of the equation resolved during the computation
 *
 * \param[in, out]  domain  pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_numeric_settings(cs_domain_t   *domain)
{
  return; /* REMOVE_LINE_FOR_USE_OF_SUBROUTINE */

  /* Modify the setting of an equation using a generic process

                   cs_equation_set(eq, key, val)

     the couple (key,val) are strings among the following choices:
     >> key: "scheme_space"
         >> val: "cdo_vb" for CDO vertex-based scheme
         >> val: "cdo_fb" for CDO face-based scheme

     >> key: "verbosity"
        >> val: "0" (default), "1", "2", ...
        The higher the more detailed information is displayed
        "1" detailed setup resume and coarse grain timer stats
        "2" fine grain timer stats

     >> key: "hodge_diff_algo" or "hodge_time_algo"
       >> val: "voronoi" or "cost" (default)
       "voronoi" leads to diagonal discrete Hodge operator but is not
       consistent for all meshes
       "cost" is more robust (i.e. it handles more general meshes but is is
       less efficient)

     >> key: "hodge_diff_coef" or "hodge_time_coef"
        This key is only useful if "cost" is set as algorithm
        >> val: "dga", "sushi", "gcr" or "1.5", "9"..
        val is either a name or a value. Notice that
        "dga" corresponds to the value 1./3.
        "sushi" corresponds to the value 1./sqrt(3.)
        "gcr" corresponds to the value 1.

     >> key: "solver_family"
        >> val: "cs" (default), "petsc", "newton" (not implemented yet)
        For using "petsc" one needs to compile Code_Saturne with the PETSc
        library

     >> key: "itsol"
        >> val: "cg" (default), "bicg", "gmres", "amg"
        "cg" is the standard conjuguate gradient algorithm
        "bicg" is BiCG-Stab2 algorithm (for non-symmetric linear systems)
        "gmres" is a robust iterative solver but not as efficient
        "amg" is an algebraic multigrid iterative solver

     >> key: "precond"
        >> val: "jacobi", "poly1", "ssor", "ilu0", "icc0", "amg", "as"
        "jacobi" diagonal preconditoner
        "poly1"  neumann polynomial of order 1
        "ssor"   symmetric successive over-relaxation
        "ilu0"   incomplete LU factorization
        "icc0"   incomplete Cholesky factorization (for symmetric matrices)
        "amg"    algebraic multigrid
        "as"     additive schwarz method

     >> key: "itsol_max_iter"
     >> key: "itsol_eps"
     >> key: "itsol_resnorm"
        >> val: "true" or "false"
  */

  cs_equation_t  *eq = cs_domain_get_equation(domain, "FVCA6.1");

  if (eq != NULL) {
    cs_equation_set(eq, "space_scheme", "cdo_fb");
    cs_equation_set(eq, "verbosity", "2");
    cs_equation_set(eq, "hodge_diff_algo", "cost");
    cs_equation_set(eq, "hodge_diff_coef", "dga");
    cs_equation_set(eq, "solver_family", "petsc");
    cs_equation_set(eq, "itsol", "cg");
    cs_equation_set(eq, "precond", "amg");
    cs_equation_set(eq, "itsol_max_iter", "2500");
    cs_equation_set(eq, "itsol_eps", "1e-12");
    cs_equation_set(eq, "itsol_resnorm", "false");
  }

}
