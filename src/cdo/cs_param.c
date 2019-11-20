/*============================================================================
 * Routines to handle basic parameter settings
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_param.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_PARAM_DBG  0

/*============================================================================
 * Global variables
 *============================================================================*/

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
    N_("HHO-P0"),
    N_("HHO-P1"),
    N_("HHO-P2")
  };

static const char
cs_param_time_scheme_name[CS_TIME_N_SCHEMES][CS_BASE_STRING_LEN] =
  { N_("Steady-state"),
    N_("Implicit"),
    N_("Explicit"),
    N_("Crank-Nicolson"),
    N_("Theta scheme")
  };

static const char
cs_param_bc_type_name[CS_PARAM_N_BC_TYPES][CS_BASE_STRING_LEN] =
  { N_("Homogeneous Dirichlet"),
    N_("Dirichlet"),
    N_("Homogeneous Neumann"),
    N_("Neumann"),
    N_("Robin"),
    N_("Sliding")
  };

static const char
cs_param_bc_enforcement_name[CS_PARAM_N_BC_ENFORCEMENTS][CS_BASE_STRING_LEN] =
  { N_("weak using an algebraic manipulation"),
    N_("weak using a big penalization coefficient"),
    N_("weak using the Nitsche method"),
    N_("weak using the symmetrized Nitsche method") };

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

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
  if (scheme == CS_SPACE_N_SCHEMES)
    return NULL;
  else
    return cs_param_space_scheme_name[scheme];
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
  if (scheme == CS_TIME_N_SCHEMES)
    return NULL;
  else
    return cs_param_time_scheme_name[scheme];
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
  if (type == CS_PARAM_N_BC_TYPES)
    return NULL;
  else
    return cs_param_bc_type_name[type];
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
  if (type == CS_PARAM_N_BC_ENFORCEMENTS)
    return NULL;
  else
    return cs_param_bc_enforcement_name[type];
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
    return  "CG";
    break;
  case CS_PARAM_ITSOL_CR3:
    return "Conjugate.Residual.3Layers";
    break;
  case CS_PARAM_ITSOL_FCG:
    return  "Flexible.CG";
    break;
  case CS_PARAM_ITSOL_GAUSS_SEIDEL:
    return  "Gauss.Seidel";
    break;
  case CS_PARAM_ITSOL_GMRES:
    return "GMRES";
    break;
  case CS_PARAM_ITSOL_JACOBI:
    return  "Jacobi";
    break;
  case CS_PARAM_ITSOL_MINRES:
    return  "MinRes";
    break;
  case CS_PARAM_ITSOL_MUMPS:
    return  "MUMPS (LU factorization)";
    break;
  case CS_PARAM_ITSOL_MUMPS_LDLT:
    return  "MUMPS (LDLT factorization)";
    break;
  case CS_PARAM_ITSOL_SYM_GAUSS_SEIDEL:
    return  "Symmetric.Gauss.Seidel";
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

  case CS_PARAM_PRECOND_BJACOB:
    return  "Block-Jacobi";
    break;
  case CS_PARAM_PRECOND_AMG:
    return  "Algebraic.MultiGrid";
    break;
  case CS_PARAM_PRECOND_AMG_BLOCK:
    return  "Algebraic.MultiGrid.ByBlock";
    break;
  case CS_PARAM_PRECOND_AS:
    return  "Additive.Schwarz";
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
  case CS_PARAM_PRECOND_ILU0:
    return  "ILU0";
    break;
  case CS_PARAM_PRECOND_ICC0:
    return  "ICC0";
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
  case CS_PARAM_AMG_HYPRE_BOOMER:
    return  "Boomer (Hypre)";
    break;
  case CS_PARAM_AMG_PETSC_GAMG:
    return  "GAMG (PETSc)";
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
