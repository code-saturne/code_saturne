/*============================================================================
 * Routines to handle cs_navsto_param_t structure
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_base.h"
#include "cs_equation.h"
#include "cs_log.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_navsto_param.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define CS_NAVSTO_PARAM_DBG  0

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private variables
 *============================================================================*/

static const char
cs_navsto_param_model_name[CS_NAVSTO_N_MODELS][CS_BASE_STRING_LEN] =
  { N_("Stokes velocity-pressure system"),
    N_("Oseen velocity-pressure system"),
    N_("Incompressible Navier-Stokes velocity-pressure system"),
    N_("Navier-Stokes with velocity-pressure unknowns and Boussinesq model")
  };

static const char
cs_navsto_param_time_state_name[CS_NAVSTO_N_TIME_STATES][CS_BASE_STRING_LEN] =
  { N_("Fully steady"),
    N_("Steaty-state as the limit of an unsteady computation"),
    N_("Fully unsteady")
  };

static const char
cs_navsto_param_coupling_name[CS_NAVSTO_N_COUPLINGS][CS_BASE_STRING_LEN] =
  { N_("Artificial compressibility algorithm"),
    N_("Artificial compressibility solved with the VPP_eps algorithm"),
    N_("Monolithic"),
    N_("Incremental projection algorithm"),
    N_("Uzawa-Augmented Lagrangian coupling")
  };

static const char _err_empty_nsp[] =
  N_(" %s: Stop setting an empty cs_navsto_param_t structure.\n"
     " Please check your settings.\n");

static const char
_space_scheme_key[CS_SPACE_N_SCHEMES][CS_BASE_STRING_LEN] =
  { "fv",
    "cdo_vb",
    "cdo_vcb",
    "cdo_eb",
    "cdo_fb",
    "hho_p0",
    "hho_p1",
    "hho_p2"
  };

static const char
_time_scheme_key[CS_TIME_N_SCHEMES][CS_BASE_STRING_LEN] =
  { "steady",
    "euler_implicit",
    "euler_explicit",
    "crank_nicolson",
    "theta_scheme"
  };

static const char
_dof_reduction_key[CS_PARAM_N_REDUCTIONS][CS_BASE_STRING_LEN] =
  { "derham",
    "average"
  };

static const char
_quad_type_key[CS_QUADRATURE_N_TYPES][CS_BASE_STRING_LEN] =
  { "none",
    "bary",
    "bary_subdiv",
    "higher",
    "highest"
  };

static const char
_adv_formulation_key[CS_PARAM_N_ADVECTION_FORMULATIONS][CS_BASE_STRING_LEN] =
  {
    "conservative",
    "non_conservative"
  };

static const char
_adv_scheme_key[CS_PARAM_N_ADVECTION_SCHEMES][CS_BASE_STRING_LEN] =
  {
    "centered",
    "cip",
    "cip_cw",
    "mix_centered_upwind",
    "samarskii",
    "sg",
    "upwind"
  };

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the \ref cs_equation_param_t structure related to the
 *         momentum equation according to the type of coupling
 *
 * \param[in]  nsp       pointer to a \ref cs_navsto_param_t structure
 *
 * \return a pointer to the corresponding \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

static inline cs_equation_param_t *
_get_momentum_param(cs_navsto_param_t    *nsp)
{
  switch (nsp->coupling) {

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY_VPP:
  case CS_NAVSTO_COUPLING_MONOLITHIC:
  case CS_NAVSTO_COUPLING_UZAWA:
    return cs_equation_param_by_name("momentum");

  case CS_NAVSTO_COUPLING_PROJECTION:
    return cs_equation_param_by_name("velocity_prediction");
    break;

  default:
    return NULL;

  }  /* Switch */
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a new structure to store all numerical parameters related
 *         to the resolution of the Navier-Stokes (NS) system
 *
 * \param[in]  boundaries     pointer to a cs_boundary_t structure
 * \param[in]  model          model related to the NS system to solve
 * \param[in]  time_state     state of the time for the NS equations
 * \param[in]  algo_coupling  algorithm used for solving the NS system
 *
 * \return a pointer to a new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_navsto_param_t *
cs_navsto_param_create(const cs_boundary_t              *boundaries,
                       cs_navsto_param_model_t           model,
                       cs_navsto_param_time_state_t      time_state,
                       cs_navsto_param_coupling_t        algo_coupling)
{
  cs_navsto_param_t  *param = NULL;
  BFT_MALLOC(param, 1, cs_navsto_param_t);

  param->boundaries = boundaries; /* shared structure */

  param->verbosity = 1;

  /* Default numerical settings */
  param->theta = 1.0;
  param->space_scheme = CS_SPACE_SCHEME_CDOFB;

  param->dof_reduction_mode = CS_PARAM_REDUCTION_AVERAGE;
  param->qtype = CS_QUADRATURE_BARY;

  /* Which equations are solved and which terms are needed */
  param->model = model;
  param->has_gravity = false;
  param->gravity[0] = param->gravity[1] = param->gravity[2] = 0.;
  param->time_state = time_state;
  if (time_state == CS_NAVSTO_TIME_STATE_FULL_STEADY)
    /* Forcing steady state in order to avoid inconsistencies */
    param->time_scheme = CS_TIME_SCHEME_STEADY;
  else
    param->time_scheme = CS_TIME_SCHEME_EULER_IMPLICIT;

  /* Resolution parameters */
  param->sles_strategy = CS_NAVSTO_SLES_EQ_WITHOUT_BLOCK;
  param->residual_tolerance = 1e-10;
  param->coupling = algo_coupling;
  param->gd_scale_coef = 1.0;    /* Default value if not set by the user */
  param->max_algo_iter = 20;

  param->adv_form   = CS_PARAM_N_ADVECTION_FORMULATIONS;
  param->adv_scheme = CS_PARAM_ADVECTION_SCHEME_UPWIND;

  /* Main set of properties */
  param->density = cs_property_add("density", CS_PROPERTY_ISO);
  param->lami_viscosity = cs_property_add("laminar_viscosity", CS_PROPERTY_ISO);

  /* Remark: As velocity and pressure may not be associated to an equation
     directly, one stores the definition of initial conditions and boundary
     conditions at this level */

  switch (algo_coupling) {

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY_VPP:
  case CS_NAVSTO_COUPLING_MONOLITHIC:
  case CS_NAVSTO_COUPLING_UZAWA:
    param->velocity_ic_is_owner = false;
    param->velocity_bc_is_owner = false;
    param->pressure_ic_is_owner = true;
    param->pressure_bc_is_owner = true;
    break;

  case CS_NAVSTO_COUPLING_PROJECTION:
    param->velocity_ic_is_owner = false;
    param->velocity_bc_is_owner = false;
    param->pressure_ic_is_owner = false;
    param->pressure_bc_is_owner = false;
    break;

  default:
    /* Nothing done */
    break;
  }

  /* Initial conditions for the pressure field */
  param->n_pressure_ic_defs = 0;
  param->pressure_ic_defs = NULL;

  /* Initial conditions for the velocity field */
  param->n_velocity_ic_defs = 0;
  param->velocity_ic_defs = NULL;

  /* Boundary conditions for the pressure field */
  param->n_pressure_bc_defs = 0;
  param->pressure_bc_defs = NULL;

  /* Boundary conditions for the velocity field */
  param->n_velocity_bc_defs = 0;
  param->velocity_bc_defs = NULL;

  /* Boundary conditions for the pressure field */
  return param;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a \ref cs_navsto_param_t structure
 *
 * \param[in, out]  param    pointer to a \ref cs_navsto_param_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_navsto_param_t *
cs_navsto_param_free(cs_navsto_param_t    *param)
{
  if (param == NULL)
    return param;

  /* Velocity initial conditions */

  if (param->n_velocity_ic_defs > 0) {

    /* Otherwise this is freed inside the related equation */
    if (param->velocity_ic_is_owner) {
      for (int i = 0; i < param->n_velocity_ic_defs; i++)
        param->velocity_ic_defs[i] = cs_xdef_free(param->velocity_ic_defs[i]);
    }
    BFT_FREE(param->velocity_ic_defs);
    param->velocity_ic_defs = NULL;

  }

  /* Pressure initial conditions */

  if (param->n_pressure_ic_defs > 0) {

    if (param->pressure_ic_is_owner) {
      for (int i = 0; i < param->n_pressure_ic_defs; i++)
        param->pressure_ic_defs[i] = cs_xdef_free(param->pressure_ic_defs[i]);
    }
    BFT_FREE(param->pressure_ic_defs);
    param->pressure_ic_defs = NULL;

  }

  /* Velocity boundary conditions */

  if (param->n_velocity_bc_defs > 0) {
    if (param->velocity_bc_is_owner) {
      for (int i = 0; i < param->n_velocity_bc_defs; i++)
        param->velocity_bc_defs[i] = cs_xdef_free(param->velocity_bc_defs[i]);
    }
    BFT_FREE(param->velocity_bc_defs);
    param->velocity_bc_defs = NULL;
  }

  /* Pressure boundary conditions */

  if (param->n_pressure_bc_defs > 0) {

    if (param->pressure_bc_is_owner) {
      for (int i = 0; i < param->n_pressure_bc_defs; i++)
        param->pressure_bc_defs[i] = cs_xdef_free(param->pressure_bc_defs[i]);
    }
    BFT_FREE(param->pressure_bc_defs);
    param->pressure_bc_defs = NULL;

  }

  /* Free the main structure */
  BFT_FREE(param);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a parameter attached to a keyname in a \ref cs_navsto_param_t
 *         structure
 *
 * \param[in, out] nsp      pointer to a \ref cs_navsto_param_t structure to set
 * \param[in]      key      key related to the member of eq to set
 * \param[in]      keyval   accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_param_set(cs_navsto_param_t    *nsp,
                    cs_navsto_key_t       key,
                    const char           *keyval)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  /* Conversion of the string to lower case */
  char val[CS_BASE_STRING_LEN];
  for (size_t i = 0; i < strlen(keyval); i++)
    val[i] = tolower(keyval[i]);
  val[strlen(keyval)] = '\0';

  switch(key) {

  case CS_NSKEY_ADVECTION_FORMULATION:
    if (strcmp(val, "conservative") == 0)
      nsp->adv_form = CS_PARAM_ADVECTION_FORM_CONSERV;
    else if (strcmp(val, "non_conservative") == 0)
      nsp->adv_form = CS_PARAM_ADVECTION_FORM_NONCONS;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Invalid val %s related to key"
                  " CS_NSKEY_ADVECTION_FORMULATION\n"
                  " Choice between conservative, non_conservative"),
                __func__, _val);
    }
    break;

  case CS_NSKEY_ADVECTION_SCHEME:
    if (strcmp(val, "upwind") == 0)
      nsp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_UPWIND;
    else if (strcmp(val, "samarskii") == 0)
      nsp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_SAMARSKII;
    else if (strcmp(val, "sg") == 0)
      nsp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_SG;
    else if (strcmp(val, "centered") == 0)
      nsp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_CENTERED;
    else if (strcmp(val, "mix_centered_upwind") == 0)
      nsp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_MIX_CENTERED_UPWIND;
    else if (strcmp(val, "cip") == 0) {
      nsp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_CIP;
      /* Automatically switch to a non-conservative formulation */
      nsp->adv_form = CS_PARAM_ADVECTION_FORM_NONCONS;
    }
    else if (strcmp(val, "cip_cw") == 0) {
      nsp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_CIP_CW;
      /* Automatically switch to a non-conservative formulation */
      nsp->adv_form = CS_PARAM_ADVECTION_FORM_NONCONS;
    }
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Invalid val %s related to key"
                  " CS_NSKEY_ADVECTION_SCHEME\n"
                  " Choice between upwind, samarskii, sg, centered, cip, cip_cw,"
                  " mix_centered_upwind"),
                __func__, _val);
    }
    break;

  case CS_NSKEY_DOF_REDUCTION:
    if (strcmp(val, "derham") == 0)
      nsp->dof_reduction_mode = CS_PARAM_REDUCTION_DERHAM;
    else if (strcmp(val, "average") == 0)
      nsp->dof_reduction_mode = CS_PARAM_REDUCTION_AVERAGE;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Invalid val %s related to key CS_NSKEY_DOF_REDUCTION\n"
                  " Choice between \"derham\" or \"average\"."),
                __func__, _val);
    }
    break;

  case CS_NSKEY_GD_SCALE_COEF:
    switch (nsp->coupling) {
    case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
    case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY_VPP:
    case CS_NAVSTO_COUPLING_UZAWA:
      nsp->gd_scale_coef = atof(val);
      break;

    case CS_NAVSTO_COUPLING_MONOLITHIC:
    case CS_NAVSTO_COUPLING_PROJECTION:
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(" %s: Trying to set the zeta parameter with the %s\n "
                 " although this will not have use in the algorithm.\n",
                 __func__, cs_navsto_param_coupling_name[nsp->coupling]);

    default:
      break;

    } /* End of switch */
    break;

  case CS_NSKEY_MAX_ALGO_ITER:
    nsp->max_algo_iter = atoi(val);
    break;

  case CS_NSKEY_QUADRATURE:
    {
      nsp->qtype = CS_QUADRATURE_NONE;

      if (strcmp(val, "bary") == 0)
        nsp->qtype = CS_QUADRATURE_BARY;
      else if (strcmp(val, "bary_subdiv") == 0)
        nsp->qtype = CS_QUADRATURE_BARY_SUBDIV;
      else if (strcmp(val, "higher") == 0)
        nsp->qtype = CS_QUADRATURE_HIGHER;
      else if (strcmp(val, "highest") == 0)
        nsp->qtype = CS_QUADRATURE_HIGHEST;
      else {
        const char *_val = val;
        bft_error(__FILE__, __LINE__, 0,
                  _(" %s: Invalid value \"%s\" for key CS_NSKEY_QUADRATURE\n"
                    " Valid choices are \"bary\", \"bary_subdiv\", \"higher\""
                    " and \"highest\"."), __func__, _val);
      }

    }
    break; /* Quadrature */

  case CS_NSKEY_RESIDUAL_TOLERANCE:
    nsp->residual_tolerance = atof(val);
    if (nsp->residual_tolerance < 0)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid value for the residual tolerance\n", __func__);
    break;

  case CS_NSKEY_SLES_STRATEGY:
    if (strcmp(val, "no_block") == 0) {
      nsp->sles_strategy = CS_NAVSTO_SLES_EQ_WITHOUT_BLOCK;
    }
    else if (strcmp(val, "block_amg_cg") == 0) {
      nsp->sles_strategy = CS_NAVSTO_SLES_BLOCK_MULTIGRID_CG;
    }
    else if (strcmp(val, "additive_gmres") == 0) {
      nsp->sles_strategy = CS_NAVSTO_SLES_ADDITIVE_GMRES_BY_BLOCK;
    }
    else if (strcmp(val, "diag_schur_gmres") == 0) {
      nsp->sles_strategy = CS_NAVSTO_SLES_DIAG_SCHUR_GMRES;
    }
    else if (strcmp(val, "upper_schur_gmres") == 0) {
      nsp->sles_strategy = CS_NAVSTO_SLES_UPPER_SCHUR_GMRES;
    }
    else if (strcmp(val, "gkb_gmres") == 0) {
      nsp->sles_strategy = CS_NAVSTO_SLES_GKB_GMRES;
    }
    else if (strcmp(val, "gkb") == 0) {
      nsp->sles_strategy = CS_NAVSTO_SLES_GKB;
    }
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Invalid val %s related to key CS_NSKEY_SLES_STRATEGY\n"
                  " Choice between \"no_block\", \"block_amg_cg\"..."),
                __func__, _val);
    }
    break;

  case CS_NSKEY_SPACE_SCHEME:
    if (strcmp(val, "cdo_fb") == 0) {
      nsp->space_scheme = CS_SPACE_SCHEME_CDOFB;
    }
    else if (strcmp(val, "hho_p0") == 0) {
      nsp->space_scheme = CS_SPACE_SCHEME_HHO_P0;
    }
    else if (strcmp(val, "hho_p1") == 0) {
      nsp->space_scheme = CS_SPACE_SCHEME_HHO_P1;
    }
    else if (strcmp(val, "hho_p2") == 0) {
      nsp->space_scheme = CS_SPACE_SCHEME_HHO_P2;
    }
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Invalid val %s related to key CS_NSKEY_SPACE_SCHEME\n"
                  " Choice between hho_{p0, p1, p2} or cdo_fb"),
                __func__, _val);
    }
    break;

  case CS_NSKEY_TIME_SCHEME:
    if (strcmp(val, "euler_implicit") == 0) {
      nsp->time_scheme = CS_TIME_SCHEME_EULER_IMPLICIT;
      nsp->theta = 1.;
    }
    else if (strcmp(val, "euler_explicit") == 0) {
      nsp->time_scheme = CS_TIME_SCHEME_EULER_EXPLICIT;
      nsp->theta = 0.;
    }
    else if (strcmp(val, "crank_nicolson") == 0) {
      nsp->time_scheme = CS_TIME_SCHEME_CRANKNICO;
      nsp->theta = 0.5;
    }
    else if (strcmp(val, "theta_scheme") == 0)
      nsp->time_scheme = CS_TIME_SCHEME_THETA;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Invalid value \"%s\" for CS_EQKEY_TIME_SCHEME\n"
                  " Valid choices are \"euler_implicit\","
                  " \"euler_explicit\"," " \"crank_nicolson\","
                  " and \"theta_scheme\"."), __func__, _val);
    }
    break;

  case CS_NSKEY_TIME_THETA:
    nsp->theta = atof(val);
    if (nsp->theta < 0. - cs_math_zero_threshold ||
        nsp->theta > 1.0 + cs_math_zero_threshold)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid value for theta\n", __func__);
    break;

  case CS_NSKEY_VERBOSITY:
    nsp->verbosity = atoi(val);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Invalid key for setting the Navier-Stokes system."),
              __func__);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Apply the numerical settings defined for the Navier-Stokes system
 *         to an equation related to this system.
 *
 * \param[in]       nsp    pointer to a \ref cs_navsto_param_t structure
 * \param[in, out]  eqp    pointer to a \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_param_transfer(const cs_navsto_param_t    *nsp,
                         cs_equation_param_t        *eqp)
{
  assert(nsp != NULL && eqp != NULL);

  /*  Set the space discretization scheme */
  const char  *ss_key = _space_scheme_key[nsp->space_scheme];

  cs_equation_set_param(eqp, CS_EQKEY_SPACE_SCHEME, ss_key);

  /*  Set the time discretization scheme */
  const char  *ts_key = _time_scheme_key[nsp->time_scheme];

  cs_equation_set_param(eqp, CS_EQKEY_TIME_SCHEME, ts_key);
  if (nsp->time_scheme == CS_TIME_SCHEME_THETA) {
    char  cvalue[36]; /* include '\0' */
    snprintf(cvalue, 35*sizeof(char), "%g", nsp->theta);
    cs_equation_set_param(eqp, CS_EQKEY_TIME_THETA, cvalue);
  }

  /*  Set the way DoFs are defined */
  const char  *dof_key = _dof_reduction_key[nsp->dof_reduction_mode];

  cs_equation_set_param(eqp, CS_EQKEY_DOF_REDUCTION, dof_key);

  /*  Set quadratures type */
  const char  *quad_key = _quad_type_key[nsp->qtype];

  /* If requested, add advection */
  if (nsp->adv_form != CS_PARAM_N_ADVECTION_FORMULATIONS) {
    /* If different from default value */
    const char *form_key = _adv_formulation_key[nsp->adv_form];
    cs_equation_set_param(eqp, CS_EQKEY_ADV_FORMULATION, form_key);

    assert(nsp->adv_scheme != CS_PARAM_N_ADVECTION_SCHEMES);
    const char *scheme_key = _adv_scheme_key[nsp->adv_scheme];
    cs_equation_set_param(eqp, CS_EQKEY_ADV_SCHEME, scheme_key);
  }

  cs_equation_set_param(eqp, CS_EQKEY_BC_QUADRATURE, quad_key);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of the main cs_navsto_param_t structure
 *
 * \param[in]  nsp    pointer to a cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_param_log(const cs_navsto_param_t    *nsp)
{
  if (nsp == NULL)
    return;

  /* Sanity checks */
  if (nsp->model == CS_NAVSTO_N_MODELS)
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid model for Navier-Stokes.\n",
              __func__);
  if (nsp->coupling == CS_NAVSTO_N_COUPLINGS)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid way of coupling the Navier-Stokes equations.\n",
              __func__);

  cs_log_printf(CS_LOG_SETUP, "  * NavSto | Verbosity: %d\n", nsp->verbosity);

  cs_log_printf(CS_LOG_SETUP, "  * NavSto | Model: %s\n",
                cs_navsto_param_model_name[nsp->model]);
  cs_log_printf(CS_LOG_SETUP, "  * NavSto | Time status: %s\n",
                cs_navsto_param_time_state_name[nsp->time_state]);

  cs_log_printf(CS_LOG_SETUP, "  * NavSto | Coupling: %s",
                cs_navsto_param_coupling_name[nsp->coupling]);
  switch (nsp->coupling) {

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY_VPP:
  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
  case CS_NAVSTO_COUPLING_UZAWA:
    cs_log_printf(CS_LOG_SETUP, " Tolerance: %5.3e\n", nsp->residual_tolerance);
    break;

  default:
    /* Nothing to print */
    cs_log_printf(CS_LOG_SETUP, "\n");
    break;
  }

  cs_log_printf(CS_LOG_SETUP, "  * NavSto | Gravity effect: %s\n",
                cs_base_strtf(nsp->has_gravity));
  if (nsp->has_gravity)
    cs_log_printf(CS_LOG_SETUP,
                  "  * NavSto | Gravity vector: [% 5.3e; % 5.3e; % 5.3e]\n",
                  nsp->gravity[0], nsp->gravity[1], nsp->gravity[2]);

  const char *space_scheme = cs_param_get_space_scheme_name(nsp->space_scheme);
  if (nsp->space_scheme < CS_SPACE_N_SCHEMES)
    cs_log_printf(CS_LOG_SETUP, "  * NavSto | Space scheme: %s\n",
                  space_scheme);
  else
    bft_error(__FILE__, __LINE__, 0,
              " %s: Undefined space scheme.", __func__);

  if (nsp->time_state != CS_NAVSTO_TIME_STATE_FULL_STEADY) {

    const char  *time_scheme = cs_param_get_time_scheme_name(nsp->time_scheme);
    if (time_scheme != NULL) {
      cs_log_printf(CS_LOG_SETUP, "  * NavSto | Time scheme: %s", time_scheme);
      if (nsp->time_scheme == CS_TIME_SCHEME_THETA)
        cs_log_printf(CS_LOG_SETUP, " with value %f\n", nsp->theta);
      else
        cs_log_printf(CS_LOG_SETUP, "\n");
    }
    else
      bft_error(__FILE__, __LINE__, 0, "%s: Invalid time scheme.", __func__);

  }

  /* Initial conditions for the velocity */
  char  prefix[256];

  cs_log_printf(CS_LOG_SETUP,
                "  * NavSto | Velocity.Init.Cond | Number of definitions %2d\n",
                nsp->n_velocity_ic_defs);

  for (int i = 0; i < nsp->n_velocity_ic_defs; i++) {
    sprintf(prefix, "  * NavSto | Velocity.Init.Cond | Definition %4d", i);
    cs_xdef_log(prefix, nsp->velocity_ic_defs[i]);
  }

  /* Initial conditions for the pressure */
  cs_log_printf(CS_LOG_SETUP,
                "  * NavSto | Pressure.Init.Cond | Number of definitions: %d\n",
                nsp->n_pressure_ic_defs);
  for (int i = 0; i < nsp->n_pressure_ic_defs; i++) {
    sprintf(prefix, "  * NavSto | Pressure.Init.Cond | Definition %4d", i);
    cs_xdef_log(prefix, nsp->pressure_ic_defs[i]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the name of the coupling algorithm
 *
 * \param[in]     coupling    a \ref cs_navsto_param_coupling_t
 *
 * \return the name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_navsto_param_get_coupling_name(cs_navsto_param_coupling_t  coupling)
{
  switch (coupling) {

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY_VPP:
  case CS_NAVSTO_COUPLING_MONOLITHIC:
  case CS_NAVSTO_COUPLING_PROJECTION:
  case CS_NAVSTO_COUPLING_UZAWA:
    return cs_navsto_param_coupling_name[coupling];

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid coupling.", __func__);
    break;

  }

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the velocity unknowns.
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set to a constant value
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      val       pointer to the value
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_add_velocity_ic_by_value(cs_navsto_param_t    *nsp,
                                   const char           *z_name,
                                   cs_real_t            *val)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  cs_xdef_t  *d = NULL;
  cs_equation_param_t *eqp = _get_momentum_param(nsp);

  if (eqp != NULL) { /* An equation related to the velocity is defined */

    d = cs_equation_add_ic_by_value(eqp, z_name, val);

  }
  else { /* No momentum equation available with the choice of velocity-pressure
            coupling */

    nsp->velocity_ic_is_owner = true;

    /* Add a new cs_xdef_t structure */
    int z_id = 0;
    if (z_name != NULL && strlen(z_name) > 0)
      z_id = (cs_volume_zone_by_name(z_name))->id;

    cs_flag_t  meta_flag = 0;
    if (z_id == 0)
      meta_flag |= CS_FLAG_FULL_LOC;

    d = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                              3,  /* dim */
                              z_id,
                              CS_FLAG_STATE_UNIFORM,
                              meta_flag,
                              val);
  }

  int  new_id = nsp->n_velocity_ic_defs;
  nsp->n_velocity_ic_defs += 1;
  BFT_REALLOC(nsp->velocity_ic_defs, nsp->n_velocity_ic_defs, cs_xdef_t *);
  nsp->velocity_ic_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the velocity unkowns.
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set according to an analytical function
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      analytic  pointer to an analytic function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_add_velocity_ic_by_analytic(cs_navsto_param_t      *nsp,
                                      const char             *z_name,
                                      cs_analytic_func_t     *analytic,
                                      void                   *input)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  cs_xdef_t  *d = NULL;
  cs_equation_param_t *eqp = _get_momentum_param(nsp);

  if (eqp != NULL) { /* An equation related to the velocity is defined */

    d = cs_equation_add_ic_by_analytic(eqp, z_name, analytic, input);

  }
  else { /* No momentum equation available with the choice of velocity-pressure
            coupling */

    nsp->velocity_ic_is_owner = true;

    /* Add a new cs_xdef_t structure */
    int z_id = 0;
    if (z_name != NULL && strlen(z_name) > 0)
      z_id = (cs_volume_zone_by_name(z_name))->id;

    cs_flag_t  meta_flag = 0;
    if (z_id == 0)
      meta_flag |= CS_FLAG_FULL_LOC;

    cs_xdef_analytic_input_t  anai = {.func = analytic,
                                      .input = input };

    d = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                              3,  /* dim */
                              z_id,
                              0,  /* state flag */
                              meta_flag,
                              &anai);
  }

  int  new_id = nsp->n_velocity_ic_defs;
  nsp->n_velocity_ic_defs += 1;
  BFT_REALLOC(nsp->velocity_ic_defs, nsp->n_velocity_ic_defs, cs_xdef_t *);
  nsp->velocity_ic_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the pressure unknowns.
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set to a constant value
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      val       pointer to the value
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_add_pressure_ic_by_value(cs_navsto_param_t    *nsp,
                                   const char           *z_name,
                                   cs_real_t            *val)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  /* Add a new cs_xdef_t structure */
  int z_id = 0;
  if (z_name != NULL && strlen(z_name) > 0)
    z_id = (cs_volume_zone_by_name(z_name))->id;

  cs_flag_t  meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                                        1,  /* dim */
                                        z_id,
                                        CS_FLAG_STATE_UNIFORM,
                                        meta_flag,
                                        val);

  int  new_id = nsp->n_pressure_ic_defs;
  nsp->n_pressure_ic_defs += 1;
  BFT_REALLOC(nsp->pressure_ic_defs, nsp->n_pressure_ic_defs, cs_xdef_t *);
  nsp->pressure_ic_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the pressure unkowns.
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set according to an analytical function
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      analytic  pointer to an analytic function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_add_pressure_ic_by_analytic(cs_navsto_param_t      *nsp,
                                      const char             *z_name,
                                      cs_analytic_func_t     *analytic,
                                      void                   *input)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  /* Add a new cs_xdef_t structure */
  int z_id = 0;
  if (z_name != NULL && strlen(z_name) > 0)
    z_id = (cs_volume_zone_by_name(z_name))->id;

  cs_flag_t  meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_analytic_input_t  anai = {.func = analytic,
                                    .input = input };

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                        1,  /* dim */
                                        z_id,
                                        0,  /* state flag */
                                        meta_flag,
                                        &anai);

  int  new_id = nsp->n_pressure_ic_defs;
  nsp->n_pressure_ic_defs += 1;
  BFT_REALLOC(nsp->pressure_ic_defs, nsp->n_pressure_ic_defs, cs_xdef_t *);
  nsp->pressure_ic_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add the definition of boundary conditions related to a fixed wall
 *         into the set of parameters for the management of the Navier-Stokes
 *         system of equations
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_set_fixed_walls(cs_navsto_param_t    *nsp)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);
  assert(nsp->boundaries != NULL);

  cs_equation_param_t *eqp = _get_momentum_param(nsp);
  cs_real_3_t  zero = {0, 0, 0};

  const cs_boundary_t  *bdy = nsp->boundaries;

  for (int i = 0; i < bdy->n_boundaries; i++) {
    if (    bdy->types[i] & CS_BOUNDARY_WALL
        && !(bdy->types[i] & CS_BOUNDARY_SLIDING_WALL)) {

      /* Homogeneous Dirichlet on the velocity field. Nothing to enforce on the
         pressure field */
      cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                                              3,    /* dim */
                                              bdy->zone_ids[i],
                                              CS_FLAG_STATE_UNIFORM, /* state */
                                              CS_CDO_BC_HMG_DIRICHLET,
                                              (void *)zero);
      int  new_id = nsp->n_velocity_bc_defs;

      nsp->n_velocity_bc_defs += 1;
      BFT_REALLOC(nsp->velocity_bc_defs, nsp->n_velocity_bc_defs, cs_xdef_t *);
      nsp->velocity_bc_defs[new_id] = d;

      cs_equation_add_xdef_bc(eqp, d);

      /* Homogeneous Neumann on the pressure field --> default BC */

    } /* Fixed wall */
  } /* Loop on domain boundaries */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add the definition of boundary conditions related to a symmetry
 *         into the set of parameters for the management of the Navier-Stokes
 *         system of equations
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_set_symmetries(cs_navsto_param_t    *nsp)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);
  assert(nsp->boundaries != NULL);

  cs_equation_param_t *eqp = _get_momentum_param(nsp);
  cs_real_t  zero = 0;

  const cs_boundary_t  *bdy = nsp->boundaries;

  for (int i = 0; i < bdy->n_boundaries; i++) {
    if (bdy->types[i] & CS_BOUNDARY_SYMMETRY) {

      /* Homogeneous Dirichlet on the normal component of the velocity field
         and homogenous Neumann on the normal stress (balance betwwen the
         pressure gradient and the viscous term) */
      cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                                              1,    /* dim */
                                              bdy->zone_ids[i],
                                              CS_FLAG_STATE_UNIFORM, /* state */
                                              CS_CDO_BC_SLIDING,
                                              (void *)&zero);

      cs_equation_add_xdef_bc(eqp, d);

      int  new_id = nsp->n_velocity_bc_defs;

      nsp->n_velocity_bc_defs += 1;
      BFT_REALLOC(nsp->velocity_bc_defs, nsp->n_velocity_bc_defs, cs_xdef_t *);
      nsp->velocity_bc_defs[new_id] = d;

      /* Homogeneous Neumann on the pressure field --> default BC (Nothing to
         do) */

    } /* Symmetry */
  } /* Loop on domain boundaries */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add the definition of boundary conditions related to outlets
 *         into the set of parameters for the management of the Navier-Stokes
 *         system of equations
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_set_outlets(cs_navsto_param_t    *nsp)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);
  assert(nsp->boundaries != NULL);

  cs_equation_param_t *eqp = _get_momentum_param(nsp);
  cs_real_33_t  zero = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

  const cs_boundary_t  *bdy = nsp->boundaries;

  int exclude_filter = CS_BOUNDARY_IMPOSED_P | CS_BOUNDARY_IMPOSED_VEL;

  for (int i = 0; i < bdy->n_boundaries; i++) {
    if (   bdy->types[i] & CS_BOUNDARY_OUTLET
        && ! (bdy->types[i] & exclude_filter)) {

      /* Add the homogeneous Neumann on the normal component */
      cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                                              9,    /* dim */
                                              bdy->zone_ids[i],
                                              CS_FLAG_STATE_UNIFORM, /* state */
                                              CS_CDO_BC_HMG_NEUMANN,
                                              (void *)&zero);
      cs_equation_add_xdef_bc(eqp, d);

      int  new_id = nsp->n_velocity_bc_defs;

      nsp->n_velocity_bc_defs += 1;
      BFT_REALLOC(nsp->velocity_bc_defs, nsp->n_velocity_bc_defs, cs_xdef_t *);
      nsp->velocity_bc_defs[new_id] = d;

      /* Homogeneous Neumann on the pressure field --> default BC.
         At the end of the day, we end up with a balance between the pressure
         gradient and the viscous term (and advection term in Navier-Stokes) */

    } /* Symmetry */
  } /* Loop on domain boundaries */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the pressure field on a boundary using a uniform value.
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" all
 *                           boundary faces are considered)
 * \param[in]      value     value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_set_pressure_bc_by_value(cs_navsto_param_t    *nsp,
                                   const char           *z_name,
                                   cs_real_t            *values)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  int  z_id = cs_get_bdy_zone_id(z_name);
  if (z_id < 0)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Zone \"%s\" does not exist.\n"
              " Please check your settings.", __func__, z_name);

  int  bdy_id = cs_boundary_id_by_zone_id(nsp->boundaries, z_id);
  if (bdy_id < 0)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Zone \"%s\" does not belong to an existing boundary.\n"
              " Please check your settings.", __func__, z_name);

  if (!(nsp->boundaries->types[bdy_id] & CS_BOUNDARY_IMPOSED_P))
    bft_error(__FILE__, __LINE__, 0,
              " %s: Zone \"%s\" is not related to a pressure boundary.\n"
              " Please check your settings.", __func__, z_name);

  /* Set the boundary condition for the pressure field */
  cs_xdef_t  *dp = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                                           1, /* dim */
                                           z_id,
                                           CS_FLAG_STATE_UNIFORM, /* state */
                                           CS_CDO_BC_DIRICHLET,
                                           (void *)values);

  int  pnew_id = nsp->n_pressure_bc_defs;

  nsp->n_pressure_bc_defs += 1;
  BFT_REALLOC(nsp->pressure_bc_defs, nsp->n_pressure_bc_defs, cs_xdef_t *);
  nsp->pressure_bc_defs[pnew_id] = dp;

  if (!nsp->pressure_bc_is_owner) {
    bft_error(__FILE__, __LINE__, 0, "%s: Not implemented yet", __func__);
#if 0 /* TODO */
    /* Retrieve the equation related to the pressure */
    cs_equation_param_t  *p_eqp = NULL;
    assert(p_eqp != NULL);
#endif
  }

  /* Add a new cs_xdef_t structure. For the momentum equation, this is a
   * homogeneous BC for the velocity */
  cs_real_33_t  zero = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

  cs_xdef_t  *du = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                                           9, /* dim */
                                           z_id,
                                           CS_FLAG_STATE_UNIFORM, /* state */
                                           CS_CDO_BC_HMG_NEUMANN,
                                           (void *)zero);

  int  unew_id = nsp->n_velocity_bc_defs;

  nsp->n_velocity_bc_defs += 1;
  BFT_REALLOC(nsp->velocity_bc_defs, nsp->n_velocity_bc_defs, cs_xdef_t *);
  nsp->velocity_bc_defs[unew_id] = du;

  cs_equation_param_t *u_eqp = _get_momentum_param(nsp);
  assert(u_eqp != NULL);
  cs_equation_add_xdef_bc(u_eqp, du);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the velocity field for a sliding wall boundary using a
 *         uniform value
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" all
 *                           boundary faces are considered)
 * \param[in]      values    array of three real values
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_set_velocity_wall_by_value(cs_navsto_param_t    *nsp,
                                     const char           *z_name,
                                     cs_real_t            *values)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  int  z_id = cs_get_bdy_zone_id(z_name);
  if (z_id < 0)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Zone \"%s\" does not exist.\n"
              " Please check your settings.", __func__, z_name);

  int  bdy_id = cs_boundary_id_by_zone_id(nsp->boundaries, z_id);
  if (bdy_id < 0)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Zone \"%s\" does not belong to an existing boundary.\n"
              " Please check your settings.", __func__, z_name);

  if (! (nsp->boundaries->types[bdy_id] & CS_BOUNDARY_SLIDING_WALL))
    bft_error(__FILE__, __LINE__, 0,
              " %s: Zone \"%s\" is not related to a sliding wall boundary.\n"
              " Please check your settings.", __func__, z_name);

  /* Add a new cs_xdef_t structure */
  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                                          3,    /* dim */
                                          z_id,
                                          CS_FLAG_STATE_UNIFORM, /* state */
                                          CS_CDO_BC_DIRICHLET,
                                          (void *)values);

  int  new_id = nsp->n_velocity_bc_defs;

  nsp->n_velocity_bc_defs += 1;
  BFT_REALLOC(nsp->velocity_bc_defs, nsp->n_velocity_bc_defs, cs_xdef_t *);
  nsp->velocity_bc_defs[new_id] = d;

  cs_equation_param_t *eqp = _get_momentum_param(nsp);
  assert(eqp != NULL);
  cs_equation_add_xdef_bc(eqp, d);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the velocity field for an inlet boundary using a uniform
 *         value
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" all
 *                           boundary faces are considered)
 * \param[in]      values    array of three real values
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_set_velocity_inlet_by_value(cs_navsto_param_t    *nsp,
                                      const char           *z_name,
                                      cs_real_t            *values)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  int  z_id = cs_get_bdy_zone_id(z_name);
  if (z_id < 0)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Zone \"%s\" does not exist.\n"
              " Please check your settings.", __func__, z_name);

  int  bdy_id = cs_boundary_id_by_zone_id(nsp->boundaries, z_id);
  if (bdy_id < 0)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Zone \"%s\" does not belong to an existing boundary.\n"
              " Please check your settings.", __func__, z_name);

  if (!(nsp->boundaries->types[bdy_id] & CS_BOUNDARY_IMPOSED_VEL))
    bft_error
      (__FILE__, __LINE__, 0,
       " %s: Zone \"%s\" is not related to an imposed velocity boundary.\n"
       " Please check your settings.", __func__, z_name);

  /* Add a new cs_xdef_t structure */
  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                                          3,    /* dim */
                                          z_id,
                                          CS_FLAG_STATE_UNIFORM, /* state */
                                          CS_CDO_BC_DIRICHLET,
                                          (void *)values);

  int  new_id = nsp->n_velocity_bc_defs;

  nsp->n_velocity_bc_defs += 1;
  BFT_REALLOC(nsp->velocity_bc_defs, nsp->n_velocity_bc_defs, cs_xdef_t *);
  nsp->velocity_bc_defs[new_id] = d;

  cs_equation_param_t *eqp = _get_momentum_param(nsp);
  assert(eqp != NULL);
  cs_equation_add_xdef_bc(eqp, d);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the velocity field for an inlet boundary using an analytical
 *         function
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" all
 *                           boundary faces are considered)
 * \param[in]      ana       pointer to an analytical function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_set_velocity_inlet_by_analytic(cs_navsto_param_t    *nsp,
                                         const char           *z_name,
                                         cs_analytic_func_t   *ana,
                                         void                 *input)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  int  z_id = cs_get_bdy_zone_id(z_name);
  if (z_id < 0)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Zone \"%s\" does not exist.\n"
              " Please check your settings.", __func__, z_name);

  int  bdy_id = cs_boundary_id_by_zone_id(nsp->boundaries, z_id);
  if (bdy_id < 0)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Zone \"%s\" does not belong to an existing boundary.\n"
              " Please check your settings.", __func__, z_name);

  if (!(nsp->boundaries->types[bdy_id] & CS_BOUNDARY_IMPOSED_VEL))
    bft_error
      (__FILE__, __LINE__, 0,
       " %s: Zone \"%s\" is not related to an imposed velocity boundary.\n"
       " Please check your settings.", __func__, z_name);

  /* Add a new cs_xdef_t structure */
  cs_xdef_analytic_input_t  anai = {.func = ana, .input = input };
  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                          3,    /* dim */
                                          z_id,
                                          0,    /* state */
                                          CS_CDO_BC_DIRICHLET,
                                          &anai);

  int  new_id = nsp->n_velocity_bc_defs;

  nsp->n_velocity_bc_defs += 1;
  BFT_REALLOC(nsp->velocity_bc_defs, nsp->n_velocity_bc_defs, cs_xdef_t *);
  nsp->velocity_bc_defs[new_id] = d;

  cs_equation_param_t *eqp = _get_momentum_param(nsp);
  cs_equation_add_xdef_bc(eqp, d);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new source term structure defined by an analytical function
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" all
 *                           cells are considered)
 * \param[in]      ana       pointer to an analytical function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_add_source_term_by_analytic(cs_navsto_param_t    *nsp,
                                      const char           *z_name,
                                      cs_analytic_func_t   *ana,
                                      void                 *input)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  cs_equation_param_t *eqp = _get_momentum_param(nsp);
  cs_xdef_t  *d = cs_equation_add_source_term_by_analytic(eqp,
                                                          z_name, ana, input);
  cs_xdef_set_quadrature(d, nsp->qtype);

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new source term structure defined by a constant value
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" all
 *                           cells are considered)
 * \param[in]      val       pointer to the value to set
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_add_source_term_by_val(cs_navsto_param_t    *nsp,
                                 const char           *z_name,
                                 cs_real_t            *val)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  cs_equation_param_t *eqp = _get_momentum_param(nsp);

  return cs_equation_add_source_term_by_val(eqp, z_name, val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new source term structure defined by an array
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" all
 *                           cells are considered)
 * \param[in]      loc       information to know where are located values
 * \param[in]      array     pointer to an array
 * \param[in]      is_owner  transfer the lifecycle to the cs_xdef_t structure
 *                           (true or false)
 * \param[in]      index     optional pointer to the array index
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_add_source_term_by_array(cs_navsto_param_t    *nsp,
                                   const char           *z_name,
                                   cs_flag_t             loc,
                                   cs_real_t            *array,
                                   bool                  is_owner,
                                   cs_lnum_t            *index)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  cs_equation_param_t *eqp = _get_momentum_param(nsp);

  return cs_equation_add_source_term_by_array(eqp, z_name, loc,
                                              array, is_owner,index);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a advection field for the Oseen problem
 *
 * \param[in, out]    nsp        pointer to a \ref cs_navsto_param_t
 * \param[in, out]    adv_fld    pointer to a \ref cs_adv_field_t
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_add_oseen_field(cs_navsto_param_t   *nsp,
                          cs_adv_field_t      *adv_fld)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  if (nsp->model != CS_NAVSTO_MODEL_OSEEN)
    bft_error(__FILE__, __LINE__, 0, " %s: Trying to set an external advection"
                                     " where there should not be one. Stopping",
                                     __func__);

  cs_equation_param_t *eqp = _get_momentum_param(nsp);

  cs_equation_add_advection(eqp, adv_fld);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
