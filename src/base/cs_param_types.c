/*============================================================================
 * Routines to handle basic parameter types
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

#include "cs_base.h"
#include "cs_log.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_param_types.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_PARAM_TYPES_DBG  0

/*============================================================================
 * Global variables
 *============================================================================*/

/* Separation lines: long, medium, short */
const char cs_sep_h1[80] =
  "=======================================================================\n";
const char cs_sep_h2[80] =
  "-----------------------------------------------------------------------\n";
const char cs_sepline[80] =
  "# =====================================================================\n";
const char cs_med_sepline[50] =
  "# ========================================\n";

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Global static variables
 *============================================================================*/

static const char
cs_param_space_scheme_name[CS_SPACE_N_SCHEMES][CS_BASE_STRING_LEN] =
  { N_("Legacy Finite Volume"),
    N_("CDO vertex-based"),
    N_("CDO vertex+cell-based"),
    N_("CDO edge-based"),
    N_("CDO face-based"),
    N_("CDO cell-based"),
    N_("HHO-P0"),
    N_("HHO-P1"),
    N_("HHO-P2")
  };

static const char
cs_param_time_scheme_name[CS_TIME_N_SCHEMES][CS_BASE_STRING_LEN] =
  { N_("Steady-state"),
    N_("1st order Forward Euler (Implicit)"),
    N_("1st order Backward Euler (Explicit)"),
    N_("Crank-Nicolson"),
    N_("Theta scheme"),
    N_("BDF2 (Implicit, 2nd order)")
  };

static const char
cs_param_adv_form_name[CS_PARAM_N_ADVECTION_FORMULATIONS][CS_BASE_STRING_LEN] =
  { N_("Conservative"),
    N_("Non-Conservative"),
    N_("Skew-symmetric"),
  };

static const char
cs_param_adv_scheme_name[CS_PARAM_N_ADVECTION_SCHEMES][CS_BASE_STRING_LEN] =
  { N_("Centered"),
    N_("Continuous interior penalty"),
    N_("Continuous interior penalty (cellwise)"),
    N_("Hybrid centered-upwind"),
    N_("Upwind with the Samarskii weight function "),
    N_("Upwind with the Scharfetter-Gummel weight function"),
    N_("Upwind"),
  };

static const char
cs_param_adv_strategy_name
[CS_PARAM_N_ADVECTION_STRATEGIES][CS_BASE_STRING_LEN] =
  { N_("Fully implicit"),
    N_("Linearized (implicit)"),
    N_("Explicit"),
  };

static const char
cs_param_adv_extrapol_name
[CS_PARAM_N_ADVECTION_EXTRAPOLATIONS][CS_BASE_STRING_LEN] =
  { N_("None"),
    N_("2nd order Taylor expansion"),
    N_("2nd order Adams-Bashforth technique"),
  };

static const char
cs_param_bc_type_name[CS_PARAM_N_BC_TYPES][CS_BASE_STRING_LEN] =
  { N_("Homogeneous Dirichlet"),
    N_("Dirichlet"),
    N_("Homogeneous Neumann"),
    N_("Neumann"),
    N_("Neumann (full)"),
    N_("Robin"),
    N_("Sliding"),
    N_("Wall prescribed")
  };

static const char
cs_param_bc_enforcement_name[CS_PARAM_N_BC_ENFORCEMENTS][CS_BASE_STRING_LEN] =
  { N_("weak using an algebraic manipulation"),
    N_("weak using a big penalization coefficient"),
    N_("weak using the Nitsche method"),
    N_("weak using the symmetrized Nitsche method") };

static const char
cs_param_nl_algo_name[CS_PARAM_N_NL_ALGOS][CS_BASE_STRING_LEN] =
  { N_("Linear algorithm"),
    N_("Picard (or fixed-point) algorithm"),
    N_("Anderson acceleration algorithm")
  };

static const char
cs_param_nl_algo_label[CS_PARAM_N_NL_ALGOS][CS_BASE_STRING_LEN] =
  { N_("None"),
    N_("Picard"),
    N_("Anderson")
  };

static const char
cs_param_precond_block_name[CS_PARAM_N_PCD_BLOCK_TYPES][CS_BASE_STRING_LEN] =
  { N_("No block preconditioner"),
    N_("Diagonal block preconditioner"),
    N_("Full diagonal block preconditioner"),
    N_("Full lower triangular block preconditioner"),
    N_("Full symmetric Gauss-Seidel block preconditioner"),
    N_("Full upper triangular block preconditioner"),
    N_("Lower triangular block preconditioner"),
    N_("Symmetric Gauss-Seidel block preconditioner"),
    N_("Upper triangular block preconditioner"),
    N_("Uzawa block preconditioner") };

static const char
cs_param_schur_approx_name[CS_PARAM_N_SCHUR_APPROX][CS_BASE_STRING_LEN] =
  { N_("None"),
    N_("Based on the diagonal"),
    N_("Identity matrix"),
    N_("Lumped inverse"),
    N_("Scaled mass matrix"),
    N_("Based on the diagonal + mass scaling"),
    N_("Lumped inverse + mass scaling") };

static const char
cs_param_dotprod_name[CS_PARAM_N_DOTPROD_TYPES][CS_BASE_STRING_LEN] =
  { N_("Classical Euclidean"),
    N_("Based on CDO quantities"),
  };

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Return true if the space scheme has degrees of freedom on faces,
 *          otherwise false
 *
 * \param[in] scheme      type of space scheme
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_param_space_scheme_is_face_based(cs_param_space_scheme_t    scheme)
{
  if (scheme == CS_SPACE_SCHEME_CDOFB  ||
      scheme == CS_SPACE_SCHEME_HHO_P0 ||
      scheme == CS_SPACE_SCHEME_HHO_P1 ||
      scheme == CS_SPACE_SCHEME_HHO_P2)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the space discretization scheme
 *
 * \param[in] scheme      type of space scheme
 *
 * \return the associated space scheme name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_space_scheme_name(cs_param_space_scheme_t    scheme)
{
  switch (scheme) {
  case CS_SPACE_SCHEME_LEGACY:
  case CS_SPACE_SCHEME_CDOVB:
  case CS_SPACE_SCHEME_CDOVCB:
  case CS_SPACE_SCHEME_CDOEB:
  case CS_SPACE_SCHEME_CDOFB:
  case CS_SPACE_SCHEME_CDOCB:
  case CS_SPACE_SCHEME_HHO_P0:
  case CS_SPACE_SCHEME_HHO_P1:
  case CS_SPACE_SCHEME_HHO_P2:
    return cs_param_space_scheme_name[scheme];

  default:
    return NULL;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the time discretization scheme
 *
 * \param[in] scheme      type of time scheme
 *
 * \return the associated time scheme name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_time_scheme_name(cs_param_time_scheme_t    scheme)
{
  switch (scheme) {
  case CS_TIME_SCHEME_STEADY:
  case CS_TIME_SCHEME_EULER_IMPLICIT:
  case CS_TIME_SCHEME_EULER_EXPLICIT:
  case CS_TIME_SCHEME_CRANKNICO:
  case CS_TIME_SCHEME_THETA:
  case CS_TIME_SCHEME_BDF2:
    return cs_param_time_scheme_name[scheme];

  default:
    return NULL;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the label associated to the advection formulation
 *
 * \param[in] adv_form      type of advection formulation
 *
 * \return the associated label
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_advection_form_name(cs_param_advection_form_t    adv_form)
{
  switch (adv_form) {
  case CS_PARAM_ADVECTION_FORM_CONSERV:
  case CS_PARAM_ADVECTION_FORM_NONCONS:
  case CS_PARAM_ADVECTION_FORM_SKEWSYM:
    return cs_param_adv_form_name[adv_form];

  default:
    return NULL;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the label of the advection scheme
 *
 * \param[in] scheme      type of advection scheme
 *
 * \return the associated advection scheme label
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_advection_scheme_name(cs_param_advection_scheme_t    scheme)
{
  switch (scheme) {
  case CS_PARAM_ADVECTION_SCHEME_CENTERED:
  case CS_PARAM_ADVECTION_SCHEME_CIP:
  case CS_PARAM_ADVECTION_SCHEME_CIP_CW:
  case CS_PARAM_ADVECTION_SCHEME_HYBRID_CENTERED_UPWIND:
  case CS_PARAM_ADVECTION_SCHEME_SAMARSKII:
  case CS_PARAM_ADVECTION_SCHEME_SG:
  case CS_PARAM_ADVECTION_SCHEME_UPWIND:
    return cs_param_adv_scheme_name[scheme];

  default:
    return NULL;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the label associated to the advection strategy
 *
 * \param[in] adv_stra      type of advection strategy
 *
 * \return the associated label
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_advection_strategy_name(cs_param_advection_strategy_t   adv_stra)
{
  switch (adv_stra) {
  case CS_PARAM_ADVECTION_IMPLICIT_FULL:
  case CS_PARAM_ADVECTION_IMPLICIT_LINEARIZED:
  case CS_PARAM_ADVECTION_EXPLICIT:
    return cs_param_adv_strategy_name[adv_stra];

  default:
    return NULL;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the label associated to the extrapolation used for the advection
 *         field
 *
 * \param[in] extrapol   type of extrapolation for the advection field
 *
 * \return the associated label
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_advection_extrapol_name(cs_param_advection_extrapol_t   extrapol)
{
  switch (extrapol) {
  case CS_PARAM_ADVECTION_EXTRAPOL_NONE:
  case CS_PARAM_ADVECTION_EXTRAPOL_TAYLOR_2:
  case CS_PARAM_ADVECTION_EXTRAPOL_ADAMS_BASHFORTH_2:
    return cs_param_adv_extrapol_name[extrapol];

  default:
    return NULL;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the type of boundary condition
 *
 * \param[in] type     type of boundary condition
 *
 * \return the associated bc name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_bc_name(cs_param_bc_type_t    type)
{
  switch (type) {
  case CS_PARAM_BC_HMG_DIRICHLET:
  case CS_PARAM_BC_DIRICHLET:
  case CS_PARAM_BC_HMG_NEUMANN:
  case CS_PARAM_BC_NEUMANN:
  case CS_PARAM_BC_NEUMANN_FULL:
  case CS_PARAM_BC_ROBIN:
  case CS_PARAM_BC_SLIDING:
  case CS_PARAM_BC_CIRCULATION:
  case CS_PARAM_BC_WALL_PRESCRIBED:
    return cs_param_bc_type_name[type];

  default:
    return NULL;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the type of enforcement of the boundary condition
 *
 * \param[in] type    type of enforcement of boundary conditions
 *
 * \return the associated name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_bc_enforcement_name(cs_param_bc_enforce_t  type)
{
  switch (type) {
  case CS_PARAM_BC_ENFORCE_ALGEBRAIC:
  case CS_PARAM_BC_ENFORCE_PENALIZED:
  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
    return cs_param_bc_enforcement_name[type];

  default:
    return NULL;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the non-linear algorithm
 *
 * \param[in] algo     type of algorithm
 *
 * \return the associated algorithm name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_nl_algo_name(cs_param_nl_algo_t   algo)
{
  switch (algo) {
  case CS_PARAM_NL_ALGO_NONE:
  case CS_PARAM_NL_ALGO_PICARD:
  case CS_PARAM_NL_ALGO_ANDERSON:
    return cs_param_nl_algo_name[algo];

  default:
    return NULL;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the label (short name) of the non-linear algorithm
 *
 * \param[in] algo     type of algorithm
 *
 * \return the associated algorithm label
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_nl_algo_label(cs_param_nl_algo_t   algo)
{
  switch (algo) {
  case CS_PARAM_NL_ALGO_NONE:
  case CS_PARAM_NL_ALGO_PICARD:
  case CS_PARAM_NL_ALGO_ANDERSON:
    return cs_param_nl_algo_label[algo];

  default:
    return NULL;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the type of dot product to apply
 *
 * \param[in] dp_type     type of dot product
 *
 * \return the associated type name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_dotprod_type_name(cs_param_dotprod_type_t   dp_type)
{
  switch (dp_type) {
  case CS_PARAM_DOTPROD_EUCLIDEAN:
  case CS_PARAM_DOTPROD_CDO:
    return cs_param_dotprod_name[dp_type];

  default:
    return NULL;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the solver
 *
 * \param[in] solver     type of iterative solver
 *
 * \return the associated solver name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_solver_name(cs_param_itsol_type_t  solver)
{
  switch (solver) {

  case CS_PARAM_ITSOL_NONE:
    return "None (PreOnly)";
    break;

  case CS_PARAM_ITSOL_AMG:
    return "Algebraic.Multigrid";
    break;
  case CS_PARAM_ITSOL_BICG:
    return "BiCG";
    break;
  case CS_PARAM_ITSOL_BICGSTAB2:
    return "BiCGstab2";
    break;
  case CS_PARAM_ITSOL_CG:
    return "CG";
    break;
  case CS_PARAM_ITSOL_CR3:
    return "Conjugate.Residual.3Layers";
    break;
  case CS_PARAM_ITSOL_FCG:
    return "Flexible.CG";
    break;
  case CS_PARAM_ITSOL_FGMRES:
    return "Flexible.GMRES";
    break;
  case CS_PARAM_ITSOL_GAUSS_SEIDEL:
    return "Gauss.Seidel";
    break;
  case CS_PARAM_ITSOL_GCR:
    return "Generalized Conjugate Residual";
    break;
  case CS_PARAM_ITSOL_GKB_CG:
    return "Golub-Kahan.BiOrthogonalization.with.CG.(inner.solver)";
    break;
  case CS_PARAM_ITSOL_GKB_GMRES:
    return "Golub-Kahan.BiOrthogonalization.with.GMRES.(inner.solver)";
    break;
  case CS_PARAM_ITSOL_GMRES:
    return "GMRES";
    break;
  case CS_PARAM_ITSOL_JACOBI:
    return "Jacobi";
    break;
  case CS_PARAM_ITSOL_MINRES:
    return "MinRes";
    break;
  case CS_PARAM_ITSOL_MUMPS:
    return "MUMPS";
    break;
  case CS_PARAM_ITSOL_SYM_GAUSS_SEIDEL:
    return "Symmetric.Gauss.Seidel";
    break;
  case CS_PARAM_ITSOL_USER_DEFINED:
    return "User-defined iterative solver";
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Invalid solver. Stop execution."), __func__);
  }

  return "";
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the preconditioner
 *
 * \param[in] precond     type of preconditioner
 *
 * \return the associated preconditioner name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_precond_name(cs_param_precond_type_t  precond)
{
  switch (precond) {

  case CS_PARAM_PRECOND_NONE:
    return  "None";
    break;

  case CS_PARAM_PRECOND_BJACOB_ILU0:
    return  "Block-Jacobi with ILU0 in each block";
    break;
  case CS_PARAM_PRECOND_BJACOB_SGS:
    return  "Block-Jacobi with symmetric Gauss-Seidel in each block";
    break;
  case CS_PARAM_PRECOND_AMG:
    return  "Algebraic.MultiGrid";
    break;
  case CS_PARAM_PRECOND_DIAG:
    return  "Diagonal";
    break;
  case CS_PARAM_PRECOND_GKB_CG:
    return "Golub-Kahan.BiOrthogonalization.with.CG.(inner.solver)";
    break;
  case CS_PARAM_PRECOND_GKB_GMRES:
    return "Golub-Kahan.BiOrthogonalization.with.GMRES.(inner.solver)";
    break;
  case CS_PARAM_PRECOND_LU:
    return  "LU";
    break;
  case CS_PARAM_PRECOND_ILU0:
    return  "ILU0";
    break;
  case CS_PARAM_PRECOND_ICC0:
    return  "ICC0";
    break;
  case CS_PARAM_PRECOND_MUMPS:
    return  "MUMPS";
    break;
  case CS_PARAM_PRECOND_POLY1:
    return  "Neumann.Poly.O1";
    break;
  case CS_PARAM_PRECOND_POLY2:
    return  "Neumann.Poly.O2";
    break;
  case CS_PARAM_PRECOND_SSOR:
    return  "SSOR";
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Invalid preconditioner. Stop execution."), __func__);
  }

  return "";
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the type of block preconditioning
 *
 * \param[in] type     type of block preconditioning
 *
 * \return the associated type name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_precond_block_name(cs_param_precond_block_t   type)
{
  switch (type) {
  case CS_PARAM_PRECOND_BLOCK_NONE:
  case CS_PARAM_PRECOND_BLOCK_DIAG:
  case CS_PARAM_PRECOND_BLOCK_FULL_DIAG:
  case CS_PARAM_PRECOND_BLOCK_FULL_LOWER_TRIANGULAR:
  case CS_PARAM_PRECOND_BLOCK_FULL_SYM_GAUSS_SEIDEL:
  case CS_PARAM_PRECOND_BLOCK_FULL_UPPER_TRIANGULAR:
  case CS_PARAM_PRECOND_BLOCK_LOWER_TRIANGULAR:
  case CS_PARAM_PRECOND_BLOCK_SYM_GAUSS_SEIDEL:
  case CS_PARAM_PRECOND_BLOCK_UPPER_TRIANGULAR:
  case CS_PARAM_PRECOND_BLOCK_UZAWA:
    return cs_param_precond_block_name[type];

  default:
    return NULL;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the type of Schur complement approximation
 *
 * \param[in] type     type of Schur complement approximation
 *
 * \return the associated type name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_schur_approx_name(cs_param_schur_approx_t   type)
{
  switch (type) {
  case CS_PARAM_SCHUR_NONE:
  case CS_PARAM_SCHUR_IDENTITY:
  case CS_PARAM_SCHUR_MASS_SCALED:
  case CS_PARAM_SCHUR_DIAG_INVERSE:
  case CS_PARAM_SCHUR_MASS_SCALED_DIAG_INVERSE:
  case CS_PARAM_SCHUR_LUMPED_INVERSE:
  case CS_PARAM_SCHUR_MASS_SCALED_LUMPED_INVERSE:
    return cs_param_schur_approx_name[type];

  default:
    return NULL;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the type of algebraic multigrid (AMG)
 *
 * \param[in] type     type of AMG
 *
 * \return the associated type name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_amg_type_name(cs_param_amg_type_t   type)
{
  switch (type) {

  case CS_PARAM_AMG_NONE:
    return  "None";
    break;
  case CS_PARAM_AMG_HYPRE_BOOMER_V:
    return  "Boomer V-cycle (Hypre)";
    break;
  case CS_PARAM_AMG_HYPRE_BOOMER_W:
    return  "Boomer W-cycle (Hypre)";
    break;
  case CS_PARAM_AMG_PETSC_GAMG_V:
    return  "GAMG V-cycle (PETSc)";
    break;
  case CS_PARAM_AMG_PETSC_GAMG_W:
    return  "GAMG W-cycle (PETSc)";
    break;
  case CS_PARAM_AMG_PETSC_PCMG:
    return  "PCMG (PETSc)";
    break;
  case CS_PARAM_AMG_HOUSE_V:
    return  "In-house (V-cycle)";
    break;
  case CS_PARAM_AMG_HOUSE_K:
    return  "In-house (K-cycle)";
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Invalid type of AMG. Stop execution."), __func__);
  }

  return "";
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
