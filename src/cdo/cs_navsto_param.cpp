/*============================================================================
 * Functions to handle cs_navsto_param_t structure
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

#include <assert.h>
#include <float.h>
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

#define CS_NAVSTO_PARAM_DBG 0

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private variables
 *============================================================================*/

static const char _err_empty_nsp[] =
  N_(" %s: Stop setting an empty cs_navsto_param_t structure.\n"
     " Please check your settings.\n");

static const char
cs_navsto_param_model_name[CS_NAVSTO_N_MODELS][CS_BASE_STRING_LEN] =
  { N_("Stokes equations"),
    N_("Oseen equations"),
    N_("Incompressible Navier-Stokes equations"),
  };

static const char
cs_navsto_param_coupling_name[CS_NAVSTO_N_COUPLINGS][CS_BASE_STRING_LEN] =
  { N_("Artificial compressibility algorithm"),
    N_("Monolithic"),
    N_("Incremental projection algo"),
  };

/* Keys to transfer settings from cs_param_navsto_t to cs_equation_param_t */

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
_dof_reduction_key[CS_PARAM_N_REDUCTIONS][CS_BASE_STRING_LEN] =
  { "derham",
    "average"
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
  case CS_NAVSTO_COUPLING_MONOLITHIC:
    return cs_equation_param_by_name("momentum");

  case CS_NAVSTO_COUPLING_PROJECTION:
    return cs_equation_param_by_name("velocity_prediction");
    break;

  default:
    return nullptr;

  }  /* Switch */
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a new structure to store all numerical parameters related
 *        to the resolution of the Navier-Stokes (NS) system
 *
 * \param[in] boundaries     pointer to a cs_boundary_t structure
 * \param[in] model          type of model related to the NS system
 * \param[in] model_flag     additional high-level model options
 * \param[in] algo_coupling  algorithm used for solving the NS system
 * \param[in] post_flag      predefined post-processings options
 *
 * \return a pointer to a new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_navsto_param_t *
cs_navsto_param_create(const cs_boundary_t          *boundaries,
                       cs_navsto_param_model_t       model,
                       cs_navsto_param_model_flag_t  model_flag,
                       cs_navsto_param_coupling_t    algo_coupling,
                       cs_navsto_param_post_flag_t   post_flag)
{
  cs_navsto_param_t *nsp = nullptr;
  BFT_MALLOC(nsp, 1, cs_navsto_param_t);

  /* Physical modelling */
  /* ------------------ */

  /* Which equations are solved and which terms are needed */

  nsp->model = model;
  nsp->model_flag = model_flag;

  /* Turbulence modelling (pointer to global structures) */

  nsp->turbulence = cs_turbulence_param_create();

  /* Main set of properties */
  /* ---------------------- */

  nsp->phys_constants = cs_get_glob_physical_constants();

  nsp->mass_density = cs_property_by_name(CS_PROPERTY_MASS_DENSITY);
  if (nsp->mass_density == nullptr)
    nsp->mass_density = cs_property_add(CS_PROPERTY_MASS_DENSITY,
                                        CS_PROPERTY_ISO);

  nsp->lam_viscosity = cs_property_add(CS_NAVSTO_LAM_VISCOSITY,
                                       CS_PROPERTY_ISO);

  if (nsp->turbulence->model->iturb == CS_TURB_NONE)
    nsp->tot_viscosity = nsp->lam_viscosity;
  else
    nsp->tot_viscosity = cs_property_add(CS_NAVSTO_TOTAL_VISCOSITY,
                                         CS_PROPERTY_ISO);

  /* Default numerical settings */
  /* -------------------------- */

  nsp->dof_reduction_mode = CS_PARAM_REDUCTION_AVERAGE;
  nsp->coupling = algo_coupling;
  nsp->space_scheme = CS_SPACE_SCHEME_CDOFB;

  /* Boussinesq term(s) */

  nsp->n_boussinesq_terms = 0;
  nsp->boussinesq_param   = nullptr;

  /* By default, one assumes a linearization of the non-linearities
   *
   * Set the default settings when a non-linear algorithm is used. This is only
   * useful if the advection term is implicit
   */

  nsp->nl_algo_type = CS_PARAM_NL_ALGO_NONE;

  nsp->nl_cvg_param.n_max_iter = 25;
  nsp->nl_cvg_param.rtol = 1e-5;
  nsp->nl_cvg_param.atol = 1e-5;
  nsp->nl_cvg_param.dtol = 1e3;

  nsp->anderson_param.n_max_dir = 6;
  nsp->anderson_param.starting_iter = 3;
  nsp->anderson_param.max_cond = -1; /* No test by default */
  nsp->anderson_param.beta = 1.0;    /* No damping by default */
  nsp->anderson_param.dp_type = CS_PARAM_DOTPROD_EUCLIDEAN;

  /* This first settings can be modified when initializing a coupling
     context (cf. cs_navsto_coupling.c) */

  /* Management of the outer resolution steps (i.e. the full system including
     the turbulence modelling or the the thermal system) */

  nsp->n_max_outer_iter = 5;
  nsp->delta_thermal_tolerance = 1e-2;

  /* Output indicators */

  nsp->verbosity = 1;
  nsp->post_flag = post_flag;

  /* Initial conditions
   * ------------------
   *
   * Remark: As velocity and pressure may not be associated to an equation
   * directly, one stores the definition of initial conditions and boundary
   * conditions at this level */

  switch (algo_coupling) {

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
    nsp->velocity_ic_is_owner = false;
    nsp->velocity_bc_is_owner = false;
    nsp->pressure_ic_is_owner = true;
    nsp->pressure_bc_is_owner = true;
    break;

  case CS_NAVSTO_COUPLING_MONOLITHIC:
    nsp->velocity_ic_is_owner = false;
    nsp->velocity_bc_is_owner = false;
    nsp->pressure_ic_is_owner = true;
    nsp->pressure_bc_is_owner = true;
    break;

  case CS_NAVSTO_COUPLING_PROJECTION:
    nsp->velocity_ic_is_owner = false;
    nsp->velocity_bc_is_owner = false;
    nsp->pressure_ic_is_owner = true;
    nsp->pressure_bc_is_owner = true;
    break;

  default:
    /* Nothing done */
    break;
  }

  /* Initial conditions for the pressure field */

  nsp->n_pressure_ic_defs = 0;
  nsp->pressure_ic_defs   = nullptr;

  /* Initial conditions for the velocity field */

  nsp->n_velocity_ic_defs = 0;
  nsp->velocity_ic_defs   = nullptr;

  /* Boundary conditions */
  /* ------------------- */

  /* Physical boundaries specific to the problem at stake */

  nsp->boundaries = boundaries; /* shared structure */

  /* Boundary conditions for the pressure field */

  nsp->n_pressure_bc_defs = 0;
  nsp->pressure_bc_defs   = nullptr;

  /* Boundary conditions for the velocity field */

  nsp->n_velocity_bc_defs = 0;
  nsp->velocity_bc_defs   = nullptr;

  /* Other conditions */
  /* ---------------- */

  /* Rescaling of the pressure */

  nsp->reference_pressure = 0.;

  return nsp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free a \ref cs_navsto_param_t structure
 *
 * \param[in, out] param  pointer to a \ref cs_navsto_param_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_navsto_param_t *
cs_navsto_param_free(cs_navsto_param_t *param)
{
  if (param == nullptr)
    return param;

  /* Turbulence modelling */

  BFT_FREE(param->turbulence);

  /* Boussinesq term(s) */

  if (param->n_boussinesq_terms > 0)
    BFT_FREE(param->boussinesq_param);

  /* Velocity initial conditions */

  if (param->n_velocity_ic_defs > 0) {

    if (param->velocity_ic_is_owner) /* Otherwise this is freed inside the
                                        related equation */
      for (int i = 0; i < param->n_velocity_ic_defs; i++)
        param->velocity_ic_defs[i] = cs_xdef_free(param->velocity_ic_defs[i]);

    BFT_FREE(param->velocity_ic_defs);
    param->velocity_ic_defs = nullptr;
  }

  /* Pressure initial conditions */

  if (param->n_pressure_ic_defs > 0) {

    if (param->pressure_ic_is_owner)
      for (int i = 0; i < param->n_pressure_ic_defs; i++)
        param->pressure_ic_defs[i] = cs_xdef_free(param->pressure_ic_defs[i]);

    BFT_FREE(param->pressure_ic_defs);
    param->pressure_ic_defs = nullptr;

  }

  /* Velocity boundary conditions */

  if (param->n_velocity_bc_defs > 0) {

    if (param->velocity_bc_is_owner)
      for (int i = 0; i < param->n_velocity_bc_defs; i++)
        param->velocity_bc_defs[i] = cs_xdef_free(param->velocity_bc_defs[i]);

    BFT_FREE(param->velocity_bc_defs);
    param->velocity_bc_defs = nullptr;

  }

  /* Pressure boundary conditions */

  if (param->n_pressure_bc_defs > 0) {

    if (param->pressure_bc_is_owner)
      for (int i = 0; i < param->n_pressure_bc_defs; i++)
        param->pressure_bc_defs[i] = cs_xdef_free(param->pressure_bc_defs[i]);

    BFT_FREE(param->pressure_bc_defs);
    param->pressure_bc_defs = nullptr;

  }

  BFT_FREE(param);

  return nullptr;
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
cs_navsto_param_set(cs_navsto_param_t *nsp,
                    cs_navsto_key_t    key,
                    const char        *keyval)
{
  if (nsp == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);
  if (keyval == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: Empty key value.\n", __func__);

  /* Conversion of the string to lower case */

  char val[CS_BASE_STRING_LEN];
  for (size_t i = 0; i < strlen(keyval); i++)
    val[i] = tolower(keyval[i]);
  val[strlen(keyval)] = '\0';

  switch (key) {

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
                __func__,
                _val);
    }
    break;

  case CS_NSKEY_NL_ALGO_MAX_ITER:
    nsp->nl_cvg_param.n_max_iter = atoi(val);
    break;

  case CS_NSKEY_NL_ALGO: {
    if (strcmp(val, "picard") == 0 || strcmp(val, "fixed-point") == 0)
      nsp->nl_algo_type = CS_PARAM_NL_ALGO_PICARD;
    else if (strcmp(val, "anderson") == 0)
      nsp->nl_algo_type = CS_PARAM_NL_ALGO_ANDERSON;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid value \"%s\" for key CS_NSKEY_NL_ALGO\n"
                " Valid choices are \"picard\", \"fixed-point\".",
                __func__,
                _val);
    }
  } break; /* Non-linear algorithm */

  case CS_NSKEY_NL_ALGO_ATOL:
    nsp->nl_cvg_param.atol = atof(val);
    if (nsp->nl_cvg_param.atol < 0)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid value for the absolute tolerance of the"
                " non-linear algorithm\n",
                __func__);
    break;

  case CS_NSKEY_NL_ALGO_DTOL:
    nsp->nl_cvg_param.dtol = atof(val);
    if (nsp->nl_cvg_param.dtol < 0)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid value for the divergence tolerance of the"
                " non-linear algorithm\n",
                __func__);
    break;

  case CS_NSKEY_NL_ALGO_RTOL:
    nsp->nl_cvg_param.rtol = atof(val);
    if (nsp->nl_cvg_param.rtol < 0)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid value for the relative tolerance of the"
                " non-linear algorithm\n",
                __func__);
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
    else if (strcmp(val, "mac_fb") == 0) {
      nsp->space_scheme = CS_SPACE_SCHEME_MACFB;
    }
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Invalid val %s related to key CS_NSKEY_SPACE_SCHEME\n"
                  " Choice between hho_{p0, p1, p2}, cdo_fb or mac_fb"),
                __func__,
                _val);
    }
    break;

  case CS_NSKEY_THERMAL_TOLERANCE:
    nsp->delta_thermal_tolerance = atof(val);
    /* If tolerance is set to a negative value then it stops the outer
       iteration process after the first iteration */
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
 * \brief Apply the numerical settings defined for the Navier-Stokes system to
 *        an equation related to this system. Be aware that the user-defined
 *        settings can be modified in this function when a different choice is
 *        set between the settings for the Navier-Stokes system and the
 *        settings for the momentum equation. The final choice is given by the
 *        settings for the Navier-Stokes system.
 *
 * \param[in]      nsp    pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] eqp    pointer to a \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_param_transfer(const cs_navsto_param_t *nsp, cs_equation_param_t *eqp)
{
  assert(nsp != nullptr && eqp != nullptr);

  /*  Set the space discretization scheme */

  const char *ss_key = _space_scheme_key[nsp->space_scheme];

  if (nsp->space_scheme != eqp->space_scheme)
    cs_equation_param_set(eqp, CS_EQKEY_SPACE_SCHEME, ss_key);

  /*  Set the way DoFs are defined */

  const char *dof_key = _dof_reduction_key[nsp->dof_reduction_mode];

  if (nsp->dof_reduction_mode != eqp->dof_reduction)
    cs_equation_param_set(eqp, CS_EQKEY_DOF_REDUCTION, dof_key);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of the main cs_navsto_param_t structure
 *
 * \param[in]  nsp    pointer to a cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_param_log(const cs_navsto_param_t *nsp)
{
  if (nsp == nullptr)
    return;

  const char navsto[16] = "  * NavSto |";

  if (nsp->model == CS_NAVSTO_N_MODELS)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid model for Navier-Stokes.\n",
              __func__);
  if (nsp->coupling == CS_NAVSTO_N_COUPLINGS)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid way of coupling the Navier-Stokes equations.\n",
              __func__);

  cs_log_printf(CS_LOG_SETUP, "%s Verbosity: %d\n", navsto, nsp->verbosity);

  /* Describe the physical modelling */

  cs_log_printf(CS_LOG_SETUP, "%s Model: %s\n",
                navsto, cs_navsto_param_get_model_name(nsp->model));

  if (nsp->model_flag & CS_NAVSTO_MODEL_GRAVITY_EFFECTS)
    cs_log_printf(CS_LOG_SETUP, "%s Model: Gravity effect activated\n", navsto);

  if (nsp->model_flag & CS_NAVSTO_MODEL_CORIOLIS_EFFECTS)
    cs_log_printf(CS_LOG_SETUP, "%s Model: Coriolis effect activated\n",
                  navsto);

  if (nsp->model_flag & CS_NAVSTO_MODEL_BOUSSINESQ) {

    cs_log_printf(CS_LOG_SETUP, "%s Model:"
                  " Boussinesq approximation activated (%d term(s))\n",
                  navsto, nsp->n_boussinesq_terms);

    for (int i = 0; i < nsp->n_boussinesq_terms; i++) {

      cs_navsto_param_boussinesq_t bp = nsp->boussinesq_param[i];

      cs_log_printf(CS_LOG_SETUP,
                    "%s Dilatation coef. = %5.3e; Reference value = %5.3e\n",
                    navsto, bp.beta, bp.var0);
    }

    if (nsp->n_boussinesq_terms == 0)
      bft_error(__FILE__, __LINE__, 0,
                "%s: Boussinesq model is activated but no parameter is given\n"
                " to define the Boussinesq term.\n",
                __func__);

  } /* Boussinesq term(s) added */

  /* Describe the space-time discretization */

  cs_log_printf(CS_LOG_SETUP, "%s Coupling: %s\n",
                navsto, cs_navsto_param_coupling_name[nsp->coupling]);

  if (cs_navsto_param_is_steady(nsp))
    cs_log_printf(CS_LOG_SETUP, "%s Time status: Steady\n", navsto);

  else
    cs_log_printf(CS_LOG_SETUP, "%s Time status: Unsteady\n", navsto);

  const char *space_scheme = cs_param_get_space_scheme_name(nsp->space_scheme);
  if (space_scheme != nullptr)
    cs_log_printf(CS_LOG_SETUP, "%s Space scheme: %s\n", navsto, space_scheme);
  else
    bft_error(__FILE__, __LINE__, 0, " %s: Undefined space scheme.", __func__);

  /* Way to define the degrees of freedom */

  if (nsp->dof_reduction_mode == CS_PARAM_REDUCTION_AVERAGE)
    cs_log_printf(CS_LOG_SETUP, "%s DoF reduction: %s\n", navsto,
                  "Average on the related entity");
  else if (nsp->dof_reduction_mode == CS_PARAM_REDUCTION_DERHAM)
    cs_log_printf(CS_LOG_SETUP, "%s DoF reduction: %s\n", navsto,
                  "potential (velociy and pressure)");
  else
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid way to define DoFs.",
              __func__);

  /* Describe if needed the SLES settings for the non-linear algorithm */

  if (nsp->nl_algo_type != CS_PARAM_NL_ALGO_NONE) {

    cs_log_printf(CS_LOG_SETUP, "%s Non-linear algo: %s\n",
                  navsto, cs_param_get_nl_algo_name(nsp->nl_algo_type));
    cs_log_printf(CS_LOG_SETUP, "%s Tolerances of non-linear algo:"
                  " rtol: %5.3e; atol: %5.3e; dtol: %5.3e; max_iter: %d\n",
                  navsto, nsp->nl_cvg_param.rtol, nsp->nl_cvg_param.atol,
                  nsp->nl_cvg_param.dtol, nsp->nl_cvg_param.n_max_iter);

    if (nsp->nl_algo_type == CS_PARAM_NL_ALGO_ANDERSON) {

      const cs_iter_algo_param_aac_t aap = nsp->anderson_param;

      cs_log_printf(CS_LOG_SETUP, "%s Anderson param: max. dir: %d; "
                    " start: %d; drop. tol: %5.3e; relax: %5.3e\n",
                    navsto, aap.n_max_dir, aap.starting_iter, aap.max_cond,
                    aap.beta);
      cs_log_printf(CS_LOG_SETUP, "%s Anderson param: Dot product type: %s\n",
                    navsto, cs_param_get_dotprod_type_name(aap.dp_type));
    }

  } /* A non-linear treatment is requested */

  /* Initial conditions for the velocity */

  cs_log_printf(CS_LOG_SETUP,
                "%s Velocity.Init.Cond | Number of definitions %2d\n",
                navsto, nsp->n_velocity_ic_defs);

  char prefix[256];
  for (int i = 0; i < nsp->n_velocity_ic_defs; i++) {
    sprintf(prefix, "%s Velocity.Init.Cond | Definition %2d", navsto, i);
    cs_xdef_log_setup(prefix, nsp->velocity_ic_defs[i]);
  }

  /* Initial conditions for the pressure */

  cs_log_printf(CS_LOG_SETUP,
                "%s Pressure.Init.Cond | Number of definitions: %d\n",
                navsto, nsp->n_pressure_ic_defs);
  for (int i = 0; i < nsp->n_pressure_ic_defs; i++) {
    sprintf(prefix, "%s Pressure.Init.Cond | Definition %2d", navsto, i);
    cs_xdef_log_setup(prefix, nsp->pressure_ic_defs[i]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a new Boussinesq term (source term for the momemtum equation)
 *
 * \param[in, out] nsp               pointer to a cs_navsto_param_t struct.
 * \param[in]      dilatation_coef   value of the dilatation coefficient
 * \param[in]      reference_value   reference value of the associated variable
 *
 * \return a pointer to the newly added structure
 */
/*----------------------------------------------------------------------------*/

cs_navsto_param_boussinesq_t *
cs_navsto_param_add_boussinesq_term(cs_navsto_param_t *nsp,
                                    cs_real_t          dilatation_coef,
                                    cs_real_t          reference_value)
{
  if (nsp == nullptr)
    return nullptr;

  nsp->model_flag |= CS_NAVSTO_MODEL_BOUSSINESQ;

  int b_id = nsp->n_boussinesq_terms;

  nsp->n_boussinesq_terms += 1;
  BFT_REALLOC(nsp->boussinesq_param,
              nsp->n_boussinesq_terms + 1,
              cs_navsto_param_boussinesq_t);

  cs_navsto_param_boussinesq_t *bp = nsp->boussinesq_param + b_id;

  bp->var0 = reference_value;
  bp->beta = dilatation_coef;

  /* bp->var (shared pointer) is set with another function */

  return bp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the array of values to consider in the Boussinesq term
 *
 * \param[in, out]  bp    pointer to a cs_navsto_param_boussinesq_t structure
 * \param[in]       var   shared pointer to the array of values to consider
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_param_set_boussinesq_array(cs_navsto_param_boussinesq_t *bp,
                                     const cs_real_t              *var)
{
  if (bp == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Boussinesq structure is empty\n", __func__);

  bp->var = var;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the \ref cs_equation_param_t structure related to the
 *         velocity equation (momentum equation in most of the cases)
 *
 * \param[in]  nsp    pointer to a cs_navsto_param_t structure
 *
 * \return a pointer to the set of parameters related to the momentum equation
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_navsto_param_get_velocity_param(const cs_navsto_param_t *nsp)
{
  cs_equation_param_t *eqp = nullptr;

  switch (nsp->coupling) {

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
  case CS_NAVSTO_COUPLING_MONOLITHIC:
    eqp = cs_equation_param_by_name("momentum");
    break;

  case CS_NAVSTO_COUPLING_PROJECTION:
          eqp = cs_equation_param_by_name("velocity_prediction");
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid coupling algorithm", __func__);
    break;
  }

  return eqp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the name of the model system of equations
 *
 * \param[in]   model    a \ref cs_navsto_param_model_t
 *
 * \return the corresponding name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_navsto_param_get_model_name(cs_navsto_param_model_t   model)
{
  switch (model) {

  case CS_NAVSTO_MODEL_STOKES:
  case CS_NAVSTO_MODEL_OSEEN:
  case CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES:
    return cs_navsto_param_model_name[model];

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid model.", __func__);
    break;
  }

  return nullptr;
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
  case CS_NAVSTO_COUPLING_MONOLITHIC:
  case CS_NAVSTO_COUPLING_PROJECTION:

    return cs_navsto_param_coupling_name[coupling];

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid coupling.", __func__);
    break;
  }

  return nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the value to consider for the reference pressure
 *
 * \param[in]  nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]  pref      value of the reference pressure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_set_reference_pressure(cs_navsto_param_t    *nsp,
                                 cs_real_t             pref)
{
  if (nsp == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  nsp->reference_pressure = pref;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply the given quadrature rule to all existing definitions under
 *        the cs_navsto_param_t structure
 *
 * \param[in, out] nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in]      qtype    type of quadrature to apply
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_param_set_quadrature_to_all(cs_navsto_param_t    *nsp,
                                      cs_quadrature_type_t  qtype)
{
  /* Loop on velocity ICs */

  for (int i = 0; i < nsp->n_velocity_ic_defs; i++) {

    cs_xdef_t *def = nsp->velocity_ic_defs[i];

    if (def->type == CS_XDEF_BY_ANALYTIC_FUNCTION) /* Otherwise not useful */
      cs_xdef_set_quadrature(def, qtype);

  }

  /* Loop on pressure ICs */

  for (int i = 0; i < nsp->n_pressure_ic_defs; i++) {

    cs_xdef_t *def = nsp->pressure_ic_defs[i];

    if (def->type == CS_XDEF_BY_ANALYTIC_FUNCTION) /* Otherwise not useful */
      cs_xdef_set_quadrature(def, qtype);

  }

  /* Loop on velocity BCs */

  for (int i = 0; i < nsp->n_velocity_bc_defs; i++) {

    cs_xdef_t *def = nsp->velocity_bc_defs[i];

    if (def->type == CS_XDEF_BY_ANALYTIC_FUNCTION) /* Otherwise not useful */
      cs_xdef_set_quadrature(def, qtype);

  }

  /* Loop on pressure BCs */

  for (int i = 0; i < nsp->n_pressure_bc_defs; i++) {

    cs_xdef_t *def = nsp->pressure_bc_defs[i];

    if (def->type == CS_XDEF_BY_ANALYTIC_FUNCTION) /* Otherwise not useful */
      cs_xdef_set_quadrature(def, qtype);

  }

  /* Set this level of quadrature to all terms of the momentum equation */

  cs_equation_param_t *eqp = _get_momentum_param(nsp);

  cs_equation_param_set_quadrature_to_all(eqp, qtype);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the velocity unknowns.
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set to a constant value
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if nullptr or "" if
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
  if (nsp == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  cs_xdef_t           *d   = nullptr;
  cs_equation_param_t *eqp = _get_momentum_param(nsp);

  if (eqp != nullptr) { /* An equation related to the velocity is defined */

    d = cs_equation_add_ic_by_value(eqp, z_name, val);
  }
  else { /* No momentum equation available with the choice of velocity-pressure
            coupling */

    nsp->velocity_ic_is_owner = true;

    /* Add a new cs_xdef_t structure */

    int z_id = 0;
    if (z_name != nullptr && strlen(z_name) > 0)
      z_id = (cs_volume_zone_by_name(z_name))->id;

    cs_flag_t meta_flag = 0;
    if (z_id == 0)
      meta_flag |= CS_FLAG_FULL_LOC;

    d = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                              3, /* dim */
                              z_id,
                              CS_FLAG_STATE_UNIFORM,
                              meta_flag,
                              val);
  }

  int new_id = nsp->n_velocity_ic_defs;
  nsp->n_velocity_ic_defs += 1;
  BFT_REALLOC(nsp->velocity_ic_defs, nsp->n_velocity_ic_defs, cs_xdef_t *);
  nsp->velocity_ic_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the velocity unknowns.
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set according to an analytical function
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if nullptr or "" if
 *                           all cells are considered)
 * \param[in]      analytic  pointer to an analytic function
 * \param[in]      input     nullptr or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_add_velocity_ic_by_analytic(cs_navsto_param_t  *nsp,
                                      const char         *z_name,
                                      cs_analytic_func_t *analytic,
                                      void               *input)
{
  if (nsp == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  cs_xdef_t           *d   = nullptr;
  cs_equation_param_t *eqp = _get_momentum_param(nsp);

  if (eqp != nullptr) { /* An equation related to the velocity is defined */

    d = cs_equation_add_ic_by_analytic(eqp, z_name, analytic, input);
  }
  else { /* No momentum equation available with the choice of velocity-pressure
            coupling */

    nsp->velocity_ic_is_owner = true;

    /* Add a new cs_xdef_t structure */

    int z_id = cs_volume_zone_id_by_name(z_name);

    cs_flag_t meta_flag = 0;
    if (z_id == 0)
      meta_flag |= CS_FLAG_FULL_LOC;

    cs_xdef_analytic_context_t anai = {
      .z_id = z_id, .func = analytic, .input = input, .free_input = nullptr
    };

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
 * \param[in]      z_name    name of the associated zone (if nullptr or "" if
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
  if (nsp == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  /* Add a new cs_xdef_t structure */

  int z_id = 0;
  if (z_name != nullptr && strlen(z_name) > 0)
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
 * \brief  Define the initial condition for the pressure unknowns.
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set according to an analytical function
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if nullptr or "" if
 *                           all cells are considered)
 * \param[in]      analytic  pointer to an analytic function
 * \param[in]      input     nullptr or pointer to a structure cast on-the-fly
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
  if (nsp == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  /* Add a new cs_xdef_t structure */

  int z_id = cs_volume_zone_id_by_name(z_name);

  cs_flag_t  meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_analytic_context_t ac
    = { .z_id = z_id, .func = analytic, .input = input, .free_input = nullptr };

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                        1,  /* dim */
                                        z_id,
                                        0,  /* state flag */
                                        meta_flag,
                                        &ac);

  int  new_id = nsp->n_pressure_ic_defs;
  nsp->n_pressure_ic_defs += 1;
  BFT_REALLOC(nsp->pressure_ic_defs, nsp->n_pressure_ic_defs, cs_xdef_t *);
  nsp->pressure_ic_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the definition of boundary conditions related to a fixed wall
 *        into the set of parameters for the management of the Navier-Stokes
 *        system of equations
 *
 * \param[in] nsp  pointer to a \ref cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_set_fixed_walls(cs_navsto_param_t *nsp)
{
  if (nsp == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);
  assert(nsp->boundaries != nullptr);

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
 * \brief Add the definition of boundary conditions related to a symmetry
 *        into the set of parameters for the management of the Navier-Stokes
 *        system of equations
 *
 * \param[in] nsp  pointer to a \ref cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_set_symmetries(cs_navsto_param_t *nsp)
{
  if (nsp == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);
  assert(nsp->boundaries != nullptr);

  cs_equation_param_t *eqp = _get_momentum_param(nsp);
  cs_real_t  zero = 0;

  const cs_boundary_t  *bdy = nsp->boundaries;

  for (int i = 0; i < bdy->n_boundaries; i++) {
    if (bdy->types[i] & CS_BOUNDARY_SYMMETRY) {

      /* Homogeneous Dirichlet on the normal component of the velocity field
         and homogeneous Neumann on the normal stress (balance between the
         pressure gradient and the viscous term) */

      cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                                              1,    /* dim */
                                              bdy->zone_ids[i],
                                              CS_FLAG_STATE_UNIFORM, /* state */
                                              CS_CDO_BC_SYMMETRY,
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
  if (nsp == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);
  assert(nsp->boundaries != nullptr);

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
                                              CS_CDO_BC_SYMMETRY,
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
 * \brief  Set the pressure field on a boundary using a uniform value.
 *
 * \param[in] nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in] z_name    name of the associated zone (if nullptr or "" all
 *                      boundary faces are considered)
 * \param[in] values    value to set
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_set_pressure_bc_by_value(cs_navsto_param_t    *nsp,
                                   const char           *z_name,
                                   cs_real_t            *values)
{
  if (nsp == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  int  z_id = cs_boundary_zone_id_by_name(z_name);
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

  /* Add a new cs_xdef_t structure. For the momentum equation, this is a
   * homogeneous Neumann BC for the velocity */

  cs_real_33_t  zero = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

  cs_xdef_t  *du = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                                           9, /* dim */
                                           z_id,
                                           CS_FLAG_STATE_UNIFORM, /* state */
                                           CS_CDO_BC_SYMMETRY,
                                           (void *)zero);

  int  unew_id = nsp->n_velocity_bc_defs;

  nsp->n_velocity_bc_defs += 1;
  BFT_REALLOC(nsp->velocity_bc_defs, nsp->n_velocity_bc_defs, cs_xdef_t *);
  nsp->velocity_bc_defs[unew_id] = du;

  cs_equation_param_t *u_eqp = _get_momentum_param(nsp);
  assert(u_eqp != nullptr);
  cs_equation_add_xdef_bc(u_eqp, du);

  return dp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the velocity field for a sliding wall boundary using a
 *         uniform value
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if nullptr or "" all
 *                           boundary faces are considered)
 * \param[in]      values    array of three real values
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_set_velocity_wall_by_value(cs_navsto_param_t    *nsp,
                                     const char           *z_name,
                                     cs_real_t            *values)
{
  if (nsp == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  int  z_id = cs_boundary_zone_id_by_name(z_name);
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
  assert(eqp != nullptr);
  cs_equation_add_xdef_bc(eqp, d);

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the velocity field for an inlet boundary using a uniform
 *         value
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if nullptr or "" all
 *                           boundary faces are considered)
 * \param[in]      values    array of three real values
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_set_velocity_inlet_by_value(cs_navsto_param_t    *nsp,
                                      const char           *z_name,
                                      cs_real_t            *values)
{
  if (nsp == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  int  z_id = cs_boundary_zone_id_by_name(z_name);
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
  assert(eqp != nullptr);
  cs_equation_add_xdef_bc(eqp, d);

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the velocity field for an inlet boundary using an analytical
 *         function
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if nullptr or "" all
 *                           boundary faces are considered)
 * \param[in]      ana       pointer to an analytical function
 * \param[in]      input     nullptr or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_set_velocity_inlet_by_analytic(cs_navsto_param_t    *nsp,
                                         const char           *z_name,
                                         cs_analytic_func_t   *ana,
                                         void                 *input)
{
  if (nsp == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  int  z_id = cs_boundary_zone_id_by_name(z_name);
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

  cs_xdef_analytic_context_t ac
    = { .z_id = z_id, .func = ana, .input = input, .free_input = nullptr };

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                          3,    /* dim */
                                          z_id,
                                          0,    /* state */
                                          CS_CDO_BC_DIRICHLET,
                                          &ac);

  int  new_id = nsp->n_velocity_bc_defs;
  nsp->n_velocity_bc_defs += 1;
  BFT_REALLOC(nsp->velocity_bc_defs, nsp->n_velocity_bc_defs, cs_xdef_t *);
  nsp->velocity_bc_defs[new_id] = d;

  cs_equation_param_t *eqp = _get_momentum_param(nsp);
  cs_equation_add_xdef_bc(eqp, d);

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the velocity field for an inlet boundary using an array
 *         of values
 *
 * \param[in]  nsp          pointer to a \ref cs_navsto_param_t structure
 * \param[in]  z_name       name of the associated zone (if nullptr or "" all
 *                          boundary faces are considered)
 * \param[in]  loc          information to know where are located values
 * \param[in]  array        pointer to an array
 * \param[in]  is_owner     transfer the lifecycle to the cs_xdef_t structure
 *                          (true or false)
 * \param[in]  full_length  if true, the size of "array" should be allocated
 *                          to the total numbers of entities related to the
 *                          given location. If false, a new list is allocated
 *                          and filled with the related subset indirection.
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_set_velocity_inlet_by_array(cs_navsto_param_t    *nsp,
                                      const char           *z_name,
                                      cs_flag_t             loc,
                                      cs_real_t            *array,
                                      bool                  is_owner,
                                      bool                  full_length)
{
  if (nsp == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  int  z_id = cs_boundary_zone_id_by_name(z_name);
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

  cs_xdef_array_context_t context = { .z_id           = z_id,
                                      .stride         = 3,
                                      .value_location = loc,
                                      .is_owner       = is_owner,
                                      .full_length    = full_length,
                                      .values         = array,
                                      /* Optional parameters */
                                      .full2subset = nullptr,
                                      .n_list_elts = 0,
                                      .elt_ids     = nullptr,
                                      .adjacency   = nullptr };

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_ARRAY,
                                          3,
                                          z_id,
                                          CS_FLAG_STATE_FACEWISE,
                                          CS_CDO_BC_DIRICHLET,
                                          (void *)&context);

  /* Build the indirection array if only a subset is used */

  if (!full_length)
    cs_xdef_array_build_full2subset(d);

  int  new_id = nsp->n_velocity_bc_defs;

  nsp->n_velocity_bc_defs += 1;
  BFT_REALLOC(nsp->velocity_bc_defs, nsp->n_velocity_bc_defs, cs_xdef_t *);
  nsp->velocity_bc_defs[new_id] = d;

  cs_equation_param_t *eqp = _get_momentum_param(nsp);
  assert(eqp != nullptr);
  cs_equation_add_xdef_bc(eqp, d);

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the velocity field for an inlet boundary using a DoF function
 *
 * \param[in]  nsp         pointer to a \ref cs_navsto_param_t structure
 * \param[in]  z_name      name of the associated zone (if nullptr or "" all
 *                         boundary faces are considered)
 * \param[in]  dof_loc     where are located DoFs
 * \param[in]  func        pointer to a cs_dof_function_t
 * \param[in]  func_input  nullptr or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_set_velocity_inlet_by_dof_func(cs_navsto_param_t    *nsp,
                                         const char           *z_name,
                                         cs_flag_t             dof_loc,
                                         cs_dof_func_t        *func,
                                         void                 *func_input)
{
  if (nsp == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  int  z_id = cs_boundary_zone_id_by_name(z_name);
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

  cs_xdef_dof_context_t dc = { .z_id         = z_id,
                               .dof_location = dof_loc,
                               .func         = func,
                               .input        = func_input,
                               .free_input   = nullptr };

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_DOF_FUNCTION,
                                          3,    /* dim */
                                          z_id,
                                          0,    /* state */
                                          CS_CDO_BC_DIRICHLET,
                                          &dc);

  int  new_id = nsp->n_velocity_bc_defs;
  nsp->n_velocity_bc_defs += 1;
  BFT_REALLOC(nsp->velocity_bc_defs, nsp->n_velocity_bc_defs, cs_xdef_t *);
  nsp->velocity_bc_defs[new_id] = d;

  cs_equation_param_t  *eqp = _get_momentum_param(nsp);
  cs_equation_add_xdef_bc(eqp, d);

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new source term structure defined by an analytical function
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if nullptr or "" all
 *                           cells are considered)
 * \param[in]      ana       pointer to an analytical function
 * \param[in]      input     nullptr or pointer to a structure cast on-the-fly
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
  if (nsp == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  cs_equation_param_t *eqp = _get_momentum_param(nsp);
  cs_xdef_t  *d = cs_equation_add_source_term_by_analytic(eqp,
                                                          z_name, ana, input);

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new source term structure defined by a constant value
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if nullptr or "" all
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
  if (nsp == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  cs_equation_param_t *eqp = _get_momentum_param(nsp);

  return cs_equation_add_source_term_by_val(eqp, z_name, val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new source term structure defined by an array
 *
 * \param[in] nsp           pointer to a \ref cs_navsto_param_t structure
 * \param[in] z_name        name of the associated zone (if nullptr or "" all
 *                          cells are considered)
 * \param[in] loc           information to know where are located values
 * \param[in] array         pointer to an array
 * \param[in] is_owner      transfer the lifecycle to the cs_xdef_t structure
 *                          (true or false)
 * \param[in] full_length   if true, array size is allocated and filled to
 *                          access the full-length array corresponding to
 *                          all locations where are defined the values
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
                                   bool                  full_length)
{
  if (nsp == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  cs_equation_param_t *eqp = _get_momentum_param(nsp);

  return cs_equation_add_source_term_by_array(eqp, z_name, loc,
                                              array, is_owner, full_length);
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
  if (nsp == nullptr)
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
