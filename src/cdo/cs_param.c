/*============================================================================
 * Routines to handle the definition and usage of material properties
 *============================================================================*/

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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_mesh_location.h"
#include "cs_field.h"
#include "cs_cdo.h"

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
 * Global static variables
 *============================================================================*/

static const char
cs_param_def_type_name[CS_PARAM_N_DEF_TYPES][CS_BASE_STRING_LEN]=
  { N_("by analytic function"),
    N_("by array"),
    N_("by law (one argument)"),
    N_("by law (two arguments)"),
    N_("by law (two arguments: scalar+vector)"),
    N_("quantity over a volume"),
    N_("by time function"),
    N_("by user function"),
    N_("by value") };

static const char
cs_param_var_type_name[CS_PARAM_N_VAR_TYPES][CS_BASE_STRING_LEN]=
  { N_("scalar"),
    N_("vector"),
    N_("tensor") };

static const char
cs_param_bc_type_name[CS_PARAM_N_BC_TYPES][CS_BASE_STRING_LEN] =
  { N_("Homogeneous Dirichlet"),
    N_("Dirichlet"),
    N_("Homogeneous Neumann"),
    N_("Neumann"),
    N_("Robin") };

static const char
cs_param_boundary_type_name[CS_PARAM_N_BOUNDARY_TYPES][CS_BASE_STRING_LEN] =
  { N_("Wall"),
    N_("Inlet"),
    N_("Outlet"),
    N_("Symmetry") };

static const char
cs_param_hodge_type_desc[CS_PARAM_N_HODGE_TYPES][CS_BASE_STRING_LEN] =
  { "VpCd",
    "EpFd",
    "FpEd",
    "EdFp",
    "CpVd"  };

static const char
cs_param_hodge_algo_desc[CS_PARAM_N_HODGE_ALGOS][CS_BASE_STRING_LEN] =
  { "Voronoi",
    "Whitney on the Barycentric Subdivision (WBS)",
    "COnsistency-STabilization splitting (COST)" };

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name related to a type of variable
 *
 * \param[in] type     cs_param_var_type_t
 *
 * \return the name associated to this type
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_var_type_name(const cs_param_var_type_t   type)
{
  return cs_param_var_type_name[type];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name related to a type of definition
 *
 * \param[in] type     cs_param_def_type_t
 *
 * \return the name associated to this type
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_def_type_name(const cs_param_def_type_t   type)
{
  return cs_param_def_type_name[type];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a cs_def_t structure
 *
 * \param[in]      def_type   type of definition (by value, function...)
 * \param[in]      var_type   type of variables (scalar, vector, tensor...)
 * \param[in]      val        value to set
 * \param[in, out] def        pointer to cs_def_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_param_set_def(cs_param_def_type_t      def_type,
                 cs_param_var_type_t      var_type,
                 const void              *val,
                 cs_def_t                *def)
{
  assert(var_type != CS_PARAM_N_VAR_TYPES);

  switch (def_type) {

  case CS_PARAM_DEF_BY_VALUE:
  case CS_PARAM_DEF_BY_QOV:
    if (val == NULL) {
      if (var_type == CS_PARAM_VAR_SCAL)
        def->get.val = 0.0;
      else if (var_type == CS_PARAM_VAR_VECT)
        def->get.vect[0] = def->get.vect[1] = def->get.vect[2] = 0.0;
      else if (var_type == CS_PARAM_VAR_TENS) {
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
            def->get.tens[i][j] = 0.0;
      }
      else {
        assert(var_type == CS_PARAM_VAR_SYMTENS);
        for (int i = 0; i < 6; i++)
          def->get.twovects[i] = 0.0;
      }
    }
    else { // val != NULL
      if (var_type == CS_PARAM_VAR_SCAL)
        def->get.val = atof(val);
      else if (var_type == CS_PARAM_VAR_VECT) {
        char s[3][32];
        sscanf(val, "%s %s %s", s[0], s[1], s[2]);
        for (int i = 0; i < 3; i++)
          def->get.vect[i] = atof(s[i]);
      }
      else if (var_type == CS_PARAM_VAR_TENS) {
        char s[9][32];
        sscanf(val, "%s %s %s %s %s %s %s %s %s",
               s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8]);
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
            def->get.tens[i][j] = atof(s[3*i+j]);
      }
      else {
        assert(var_type == CS_PARAM_VAR_SYMTENS);
        char s[6][32];
        sscanf(val, "%s %s %s %s %s %s", s[0], s[1], s[2], s[3], s[4], s[5]);
        for (int i = 0; i < 6; i++)
          def->get.twovects[i] = atof(s[i]);
      }
    }
    break;

  case CS_PARAM_DEF_BY_ANALYTIC_FUNCTION:
    if (val == NULL)
      def->analytic = NULL;
    else
      def->analytic = (cs_analytic_func_t *)val;
    break;

  case CS_PARAM_DEF_BY_TIME_FUNCTION:
    if (val == NULL)
      def->time_func = NULL;
    else
      def->time_func = (cs_timestep_func_t *)val;
    break;

  case CS_PARAM_DEF_BY_LAW_ONESCA:
    if (val == NULL)
      def->law1_func = NULL;
    else
      def->law1_func = (cs_onevar_law_func_t *)val;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " This type of definition is not handled yet.\n"
              " Please modify your settings.");
    break;

  } /* end of switch on def_type */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a cs_get_t structure
 *
 * \param[in]      var_type   type of variables (scalar, vector, tensor...)
 * \param[in]      val        value to set
 * \param[in, out] get        pointer to cs_get_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_param_set_get(cs_param_var_type_t      var_type,
                 const char              *val,
                 cs_get_t                *get)
{
  assert(var_type != CS_PARAM_N_VAR_TYPES);

  if (val == NULL) {

    if (var_type == CS_PARAM_VAR_SCAL)
      get->val = 0.0;
    else if (var_type == CS_PARAM_VAR_VECT)
      get->vect[0] = get->vect[1] = get->vect[2] = 0.0;
    else if (var_type == CS_PARAM_VAR_TENS) {
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          get->tens[i][j] = 0.0;
    }
    else if (var_type == CS_PARAM_VAR_SYMTENS) {
      for (int i = 0; i < 6; i++)
        get->twovects[i] = 0.0;
    }
    else
      bft_error(__FILE__, __LINE__, 0, _(" Invalid type of variable."));

  }
  else { // val != NULL

    if (var_type == CS_PARAM_VAR_SCAL)
      get->val = atof(val);
    else if (var_type == CS_PARAM_VAR_VECT) {
      char s[3][32];
      sscanf(val, "%s %s %s", s[0], s[1], s[2]);
      for (int i = 0; i < 3; i++)
        get->vect[i] = atof(s[i]);
    }
    else if (var_type == CS_PARAM_VAR_TENS) {
      char s[9][32];
      sscanf(val, "%s %s %s %s %s %s %s %s %s",
             s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8]);
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          get->tens[i][j] = atof(s[3*i+j]);
    }
    else if (var_type == CS_PARAM_VAR_SYMTENS) {
      char s[6][32];
      sscanf(val, "%s %s %s %s %s %s", s[0], s[1], s[2], s[3], s[4], s[5]);
      for (int i = 0; i < 6; i++)
        get->twovects[i] = atof(s[i]);
    }
    else
      bft_error(__FILE__, __LINE__, 0, _(" Invalid type of variable."));

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a new cs_param_bc_t structure
 *
 * \param[in]  default_bc     default boundary condition
 *
 * \return a pointer to the new structure (free with cs_equation_param_t)
 */
/*----------------------------------------------------------------------------*/

cs_param_bc_t *
cs_param_bc_create(cs_param_bc_type_t  default_bc)
{
  cs_param_bc_t  *bc = NULL;

  BFT_MALLOC(bc, 1, cs_param_bc_t);

  bc->default_bc = default_bc;
  /* Initialization by default */
  bc->enforcement = CS_PARAM_BC_ENFORCE_STRONG;
  bc->quad_type = CS_QUADRATURE_BARY;
  bc->use_subdiv = false;

  bc->n_defs = 0;
  bc->defs = NULL;

  return bc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a cs_param_bc_def_t structure
 *
 * \param[in, out] bcpd       pointer to cs_param_bc_def_t struct. to set
 * \param[in]      loc_id     id related to a cs_mesh_location_t
 * \param[in]      bc_type    generic type of admissible boundary conditions
 * \param[in]      var_type   type of variables (scalar, vector, tensor...)
 * \param[in]      def_type   by value, function...
 * \param[in]      coef1      access to the value of the first coef
 * \param[in]      coef2      access to the value of the second coef (optional)
 */
/*----------------------------------------------------------------------------*/

void
cs_param_bc_def_set(cs_param_bc_def_t      *bcpd,
                    int                     loc_id,
                    cs_param_bc_type_t      bc_type,
                    cs_param_var_type_t     var_type,
                    cs_param_def_type_t     def_type,
                    const void             *coef1,
                    const void             *coef2)
{
  if (bcpd == NULL)
    return;

  /* Sanity checks */
  assert(def_type != CS_PARAM_N_DEF_TYPES);
  assert(bc_type != CS_PARAM_N_BC_TYPES);

  bcpd->loc_id = loc_id;
  bcpd->var_type = var_type;
  bcpd->bc_type = bc_type;
  bcpd->def_type = def_type;

  cs_param_set_def(def_type, var_type, coef1, &(bcpd->def_coef1));
  cs_param_set_def(def_type, var_type, coef2, &(bcpd->def_coef2));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the type of boundary condition
 *
 * \param[in] bc_type     type of boundary condition
 *
 * \return the associated bc name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_bc_name(cs_param_bc_type_t  bc)
{
  switch(bc) {

  case CS_PARAM_BC_HMG_DIRICHLET:
    return "Homogeneous Dirichlet";
    break;
  case CS_PARAM_BC_DIRICHLET:
    return "Dirichlet";
    break;
  case CS_PARAM_BC_HMG_NEUMANN:
    return "Homogeneous Neumann";
  case CS_PARAM_BC_NEUMANN:
    return "Neumann";
  case CS_PARAM_BC_ROBIN:
    return "Robin";
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid BC type. Stop execution."));
  }

  return "NULL"; // avoid a warning
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the type of enforcement of the boundary condition
 *
 * \param[in] bc_enforce    type of enforcement of boundary conditions
 *
 * \return the associated name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_get_bc_enforcement_name(cs_param_bc_enforce_t  type)
{
  switch(type) {

  case CS_PARAM_BC_ENFORCE_STRONG:
    return "strong";
    break;
  case CS_PARAM_BC_ENFORCE_WEAK_PENA:
    return "weak with a big penalization coefficient";
    break;
  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
    return "weak using the Nitsche method";
  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
    return "weak using the symmetrized Nitsche method";
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of enforcement. Stop execution."));
  }

  return "NULL"; // avoid a warning
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of the type of reaction term
 *
 * \param[in] r_type     type of reaction term
 *
 * \return the name associated with this type of reaction term
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_reaction_get_type_name(cs_param_reaction_type_t  r_type)
{
  switch (r_type) {
  case CS_PARAM_REACTION_TYPE_LINEAR:
    return "Linear";
    break;
  case CS_PARAM_N_REACTION_TYPES:
    return "Not set";
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of reaction term. Stop execution."));
  }

  return "NULL"; // avoid a warning
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of algorithm related to a discrete Hdoge operator
 *
 * \param[in] h_info     cs_param_hodge_t structure
 *
 * \return the name of the algorithm
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_hodge_get_algo_name(const cs_param_hodge_t   h_info)
{
  return cs_param_hodge_algo_desc[h_info.algo];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the type of discrete Hodge operator
 *
 * \param[in] h_info     cs_param_hodge_t structure
 *
 * \return the name of the type
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_hodge_get_type_name(const cs_param_hodge_t   h_info)
{
  return cs_param_hodge_type_desc[h_info.type];
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
  case CS_PARAM_ITSOL_CG:
    return  "CG";
    break;
  case CS_PARAM_ITSOL_BICG:
    return "BiCG";
    break;
  case CS_PARAM_ITSOL_BICGSTAB2:
    return "BiCGstab2";
    break;
  case CS_PARAM_ITSOL_CR3:
    return "Conjugate.Residual.3Layers";
    break;
  case CS_PARAM_ITSOL_GMRES:
    return "GMRES";
    break;
  case CS_PARAM_ITSOL_AMG:
    return "Algebraic.Multigrid";
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid solver. Stop execution."));
  }

  return "NULL";
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
  case CS_PARAM_PRECOND_DIAG:
    return  "Diagonal";
    break;
  case CS_PARAM_PRECOND_BJACOB:
    return  "Block-Jacobi";
    break;
  case CS_PARAM_PRECOND_POLY1:
    return  "Neumann.Poly.O1";
    break;
  case CS_PARAM_PRECOND_SSOR:
    return  "SSOR";
    break;
  case CS_PARAM_PRECOND_ILU0:
    return  "ILU0";
    break;
  case CS_PARAM_PRECOND_ICC0:
    return  "ICC0";
    break;
  case CS_PARAM_PRECOND_AMG:
    return  "Algebraic.MultiGrid";
    break;
  case CS_PARAM_PRECOND_AS:
    return  "Additive.Schwarz";
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid preconditioner. Stop execution."));
  }

  return "NULL";
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
