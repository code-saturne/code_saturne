/*============================================================================
 * Main functions dedicated to the "groundwater flow" module
 *============================================================================*/

/* VERS */

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
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_array.h"
#include "cs_equation_bc.h"
#include "cs_field.h"
#include "cs_gwf_priv.h"
#include "cs_gwf_sspf.h"
#include "cs_gwf_tpf.h"
#include "cs_gwf_uspf.h"
#include "cs_log.h"
#include "cs_parall.h"
#include "cs_param_types.h"
#include "cs_physical_model.h"
#include "cs_post.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gwf.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_gwf.c

  \brief Main high-level functions dedicated to groundwater flows when using
         CDO schemes
*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/* Redefined names of function from cs_math to get shorter names */

#define _dp3 cs_math_3_dot_product

#define CS_GWF_DBG 0

/*============================================================================
 * Local definitions
 *============================================================================*/

static const char
cs_gwf_model_name[CS_GWF_N_MODEL_TYPES][CS_BASE_STRING_LEN] =
  { N_("Saturated single-phase model"),
    N_("Unsaturated single-phase model"),
    N_("Miscible two-phase model"),
    N_("Immiscible two-phase model")
  };

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

static const char _err_empty_gw[] =
  " Stop execution. The structure related to the groundwater module is empty.\n"
  " Please check your settings.\n";

static cs_gwf_t *cs_gwf_main_structure = nullptr;

/* Pointer to shared structures (owned by a cs_domain_t structure) */

static const cs_cdo_quantities_t  *cs_cdo_quant;
static const cs_cdo_connect_t  *cs_cdo_connect;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the advection field related to the Darcy flux in the liquid
 *        phase
 *
 * \param[in]  gw     pointer to the main (high-level) GWF structure
 *
 * \return a pointer to a cs_adv_field_t structure or nullptr
 */
/*----------------------------------------------------------------------------*/

static cs_adv_field_t *
_get_l_adv_field(const cs_gwf_t   *gw)
{
  switch (gw->model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE: {
    cs_gwf_sspf_t *mc = (cs_gwf_sspf_t *)gw->model_context;

    return mc->darcy->adv_field;
  } break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE: {
    cs_gwf_uspf_t *mc = (cs_gwf_uspf_t *)gw->model_context;

    return mc->darcy->adv_field;
  } break;

  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE: {
    cs_gwf_tpf_t *mc = (cs_gwf_tpf_t *)gw->model_context;

    if (mc->l_darcy != nullptr)
      return mc->l_darcy->adv_field;
  } break;

  default:
    break;

  } /* Switch on the model */

  return nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a structure dedicated to manage groundwater flows
 *
 * \return a pointer to a new allocated cs_gwf_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_gwf_t *
_gwf_create(void)
{
  cs_gwf_t *gw = nullptr;

  BFT_MALLOC(gw, 1, cs_gwf_t);

  /* Default initialization */

  gw->verbosity = 0;
  gw->model = CS_GWF_N_MODEL_TYPES;
  gw->flag = 0;
  gw->post_flag = CS_GWF_POST_DARCY_FLUX_BALANCE;

  /* Property common to all model */

  gw->abs_permeability = nullptr;

  /* Modelling context */

  gw->model_context = nullptr;

  return gw;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the groundwater flow module has been activated
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_gwf_is_activated(void)
{
  if (cs_gwf_main_structure == nullptr)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the module dedicated to groundwater flows
 *
 * \param[in] model           type of physical modelling
 * \param[in] option_flag     optional flag to specify this module
 * \param[in] post_flag       optional automatic postprocessing
 *
 * \return a pointer to a new allocated groundwater flow structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_t *
cs_gwf_activate(cs_gwf_model_type_t      model,
                cs_flag_t                option_flag,
                cs_flag_t                post_flag)
{
  cs_gwf_t  *gw = _gwf_create();

  /* Set the physical model type */

  cs_glob_physical_model_flag[CS_GROUNDWATER] = 1;

  /* Store the pointer to the groundwater flow structure */

  cs_gwf_main_structure = gw;

  gw->model = model;
  gw->flag = option_flag;

  /* Add the porosity property */

  gw->soil_porosity = cs_property_add("soil_porosity", CS_PROPERTY_ISO);

  /* Allocate and initialize each model context (mc) */

  switch (model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    gw->model_context = cs_gwf_sspf_create();
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    gw->post_flag |= CS_GWF_POST_LIQUID_SATURATION;
    gw->model_context = cs_gwf_uspf_create();
    break;

  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE:
    gw->post_flag |= CS_GWF_POST_LIQUID_SATURATION;
    gw->model_context = cs_gwf_tpf_create(model);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for the GroundWater Flow module.\n",
              __func__);

  }

  /* Now one can set the post_flag (need the initializtion of the model
     context) */

  cs_gwf_set_post_options(post_flag, false);

  return gw;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all structures related to groundwater flows
 *
 * \return a nullptr pointer
 */
/*----------------------------------------------------------------------------*/

cs_gwf_t *
cs_gwf_destroy_all(void)
{
  if (cs_gwf_main_structure == nullptr)
    return nullptr;

  cs_gwf_t  *gw = cs_gwf_main_structure;

  switch (gw->model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE: {
    cs_gwf_sspf_t *mc = (cs_gwf_sspf_t *)gw->model_context;

    cs_gwf_sspf_free(&(mc));
  } break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE: {
    cs_gwf_uspf_t *mc = (cs_gwf_uspf_t *)gw->model_context;

    cs_gwf_uspf_free(&(mc));
  } break;

  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE: {
    cs_gwf_tpf_t *mc = (cs_gwf_tpf_t *)gw->model_context;

    cs_gwf_tpf_free(&(mc));
  } break;

  default:
    /* Nothing to do */
    break;
  }

  /* Free all soils */

  cs_gwf_soil_free_all();

  /* Free all tracers */

  cs_gwf_tracer_free_all();

  BFT_FREE(gw);

  /* Fields, equations, advection fields and properties are freed elsewhere */

  return nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of the main cs_gwf_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_log_setup(void)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == nullptr)
    return;
  assert(gw->model < CS_GWF_N_MODEL_TYPES);

  cs_log_printf(CS_LOG_SETUP, "\nSummary of the groundwater module\n");
  cs_log_printf(CS_LOG_SETUP, "%s", cs_sep_h1);

  /* Display information on the general options */

  if (gw->flag & CS_GWF_GRAVITATION)
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Gravitation: *True*\n");
  else
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Gravitation: *False*\n");

  if (gw->flag & CS_GWF_ENFORCE_DIVERGENCE_FREE)
    cs_log_printf(CS_LOG_SETUP,
                  "  * GWF | Enforce the divergence-free constraint"
                  " for the Darcy flux\n");
  if (gw->flag & CS_GWF_FORCE_RICHARDS_ITERATIONS)
    cs_log_printf(CS_LOG_SETUP,
                  "  * GWF | Force to solve Richards equation"
                  " at each time step\n");
  if (gw->flag & CS_GWF_RESCALE_HEAD_TO_ZERO_MEAN_VALUE)
    cs_log_printf(CS_LOG_SETUP,
                  "  * GWF | Rescale head w.r.t zero mean value\n");

  /* Display information on the post-processing options */

  bool  post_capacity =
    (gw->post_flag & CS_GWF_POST_SOIL_CAPACITY) ? true : false;
  bool  post_liquid_saturation =
    (gw->post_flag & CS_GWF_POST_LIQUID_SATURATION) ? true : false;
  bool  post_permeability =
    (gw->post_flag & CS_GWF_POST_PERMEABILITY) ? true : false;
  bool  post_gas_density =
    (gw->post_flag & CS_GWF_POST_GAS_MASS_DENSITY) ? true : false;
  bool  post_soil_state =
    (gw->post_flag & CS_GWF_POST_SOIL_STATE) ? true : false;

  cs_log_printf(CS_LOG_SETUP, "  * GWF | Post:"
                " Soil capacity %s Liquid saturation %s Permeability %s\n",
                cs_base_strtf(post_capacity),
                cs_base_strtf(post_liquid_saturation),
                cs_base_strtf(post_permeability));

  if (gw->model == CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE ||
      gw->model == CS_GWF_MODEL_MISCIBLE_TWO_PHASE)
    cs_log_printf(CS_LOG_SETUP, "  * GWF | Post: Gas mass density %s"
                  " Soil state %s\n",
                  cs_base_strtf(post_gas_density),
                  cs_base_strtf(post_soil_state));

  bool  do_balance =
    (gw->post_flag & CS_GWF_POST_DARCY_FLUX_BALANCE) ? true : false;
  bool  do_divergence =
    (gw->post_flag & CS_GWF_POST_DARCY_FLUX_DIVERGENCE) ? true : false;
  bool  post_boundary =
    (gw->post_flag & CS_GWF_POST_DARCY_FLUX_AT_BOUNDARY) ? true : false;

  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Darcy Flux:"
                " Balance %s Divergence %s At boundary faces: %s\n",
                cs_base_strtf(do_balance),
                cs_base_strtf(do_divergence),
                cs_base_strtf(post_boundary));

  /* Main options */

  cs_log_printf(CS_LOG_SETUP,
                "  * GWF | Model: **%s**\n", cs_gwf_model_name[gw->model]);

  switch(gw->model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    cs_gwf_sspf_log_setup((cs_gwf_sspf_t *)gw->model_context);
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    cs_gwf_uspf_log_setup((cs_gwf_uspf_t *)gw->model_context);
    break;

  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
    cs_gwf_tpf_log_setup((cs_gwf_tpf_t *)gw->model_context);
    break;

  default:
    break;

  } /* Type of model */

  /* Soils */

  cs_gwf_soil_log_setup();

  /* Tracers and decay chains */

  cs_gwf_tracer_log_all();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the main structure which manages a two-phase flow model
 *
 * \return a pointer to the structure cs_gwf_tpf_t
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tpf_t *
cs_gwf_get_two_phase_model(void)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  if (gw->model != CS_GWF_MODEL_MISCIBLE_TWO_PHASE &&
      gw->model != CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid model. One expects a two-phase flow model.\n",
              __func__);

  cs_gwf_tpf_t *mc = (cs_gwf_tpf_t *)gw->model_context;

  assert(mc != nullptr);
  return mc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the numerical options related to the two phase flow models
 *
 * \param[in] approx                          type of coefficient approximation
 * \param[in] solver                          type of solver
 * \param[in] use_incremental_solver          true/false
 * \param[in] use_diffusion_view_for_darcy    true/false
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_two_phase_numerical_options(cs_gwf_tpf_approx_type_t   approx,
                                       cs_gwf_tpf_solver_type_t   solver,
                                       bool       use_incremental_solver,
                                       bool       use_diffusion_view_for_darcy)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  cs_gwf_tpf_t *mc = (cs_gwf_tpf_t *)gw->model_context;
  assert(mc != nullptr);

  mc->approx_type = approx;
  mc->solver_type = solver;

  switch (solver) {

  case CS_GWF_TPF_SOLVER_PCPG_COUPLED:
    mc->use_coupled_solver = true;
    mc->use_diffusion_view_for_darcy = use_diffusion_view_for_darcy;
    mc->use_incremental_solver = use_incremental_solver;
    break;

  case CS_GWF_TPF_SOLVER_PLPC_COUPLED:
    mc->use_coupled_solver = true;
    mc->use_diffusion_view_for_darcy = true; /* No other choice */
    mc->use_incremental_solver = use_incremental_solver;

    if (!use_diffusion_view_for_darcy) {
      cs_base_warn(__FILE__, __LINE__);
      cs_log_printf(CS_LOG_WARNINGS,
                    "%s: Change an invalid user setting:\n"
                    "    Use a diffusion viewpoint for the Darcy term.\n",
                    __func__);
    }
    break;

  case CS_GWF_TPF_SOLVER_PLPG_SEGREGATED:
    mc->use_coupled_solver = false;
    mc->use_diffusion_view_for_darcy = true; /* No other choice */
    mc->use_incremental_solver = true;       /* No other choice */

    if (!use_diffusion_view_for_darcy) {
      cs_base_warn(__FILE__, __LINE__);
      cs_log_printf(CS_LOG_WARNINGS,
                    "%s: Change an invalid user setting:\n"
                    "    Use a diffusion viewpoint for the Darcy term.\n",
                    __func__);
    }

    if (!use_incremental_solver) {
      cs_base_warn(__FILE__, __LINE__);
      cs_log_printf(CS_LOG_WARNINGS,
                    "%s: Change an invalid user setting:\n"
                    "    Force an incremental resolution.\n", __func__);
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid setting", __func__);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the parameters defining the two-phase flow model in the miscible
 *        case. Use SI unit if not prescribed otherwise.
 *
 * \param[in] l_mass_density   mass density of the main liquid component
 * \param[in] l_viscosity      viscosity in the liquid phase (Pa.s)
 * \param[in] g_viscosity      viscosity in the gas phase (Pa.s)
 * \param[in] l_diffusivity_h  diffusivity of the main gas component in the
 *                             liquid phase
 * \param[in] h_molar_mass     molar mass of the main gas component
 * \param[in] ref_temperature  reference temperature in Kelvin
 * \param[in] henry_constant   constant in the Henry law
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_miscible_two_phase_model(cs_real_t       l_mass_density,
                                    cs_real_t       l_viscosity,
                                    cs_real_t       g_viscosity,
                                    cs_real_t       l_diffusivity_h,
                                    cs_real_t       h_molar_mass,
                                    cs_real_t       ref_temperature,
                                    cs_real_t       henry_constant)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));
  if (gw->model != CS_GWF_MODEL_MISCIBLE_TWO_PHASE)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid model.\n"
              "%s: One expects a miscible two-phase flow model.\n",
              __func__, __func__);

  cs_gwf_tpf_t *mc = (cs_gwf_tpf_t *)gw->model_context;

  assert(mc != nullptr);
  if (mc->is_miscible == false)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid model.\n"
              "%s: One expects a miscible two-phase flow model.\n",
              __func__, __func__);

  assert(l_mass_density > 0);
  assert(ref_temperature > 0);  /* In Kelvin */
  assert(h_molar_mass > 0);
  assert(l_viscosity > 0 && g_viscosity > 0);

  /* Set the parameters */

  mc->l_mass_density = l_mass_density;
  mc->l_viscosity = l_viscosity;
  mc->g_viscosity = g_viscosity;
  mc->l_diffusivity_h = l_diffusivity_h;
  mc->h_molar_mass = h_molar_mass;
  mc->ref_temperature = ref_temperature;
  mc->henry_constant = henry_constant;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the parameters defining the immiscible two-phase flow model.
 *        Use SI unit if not prescribed otherwise.
 *
 * \param[in] l_mass_density   mass density of the main liquid component
 * \param[in] l_viscosity      viscosity in the liquid phase (Pa.s)
 * \param[in] g_viscosity      viscosity in the gas phase (Pa.s)
 * \param[in] h_molar_mass     molar mass of the main gas component
 * \param[in] ref_temperature  reference temperature in Kelvin
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_immiscible_two_phase_model(cs_real_t       l_mass_density,
                                      cs_real_t       l_viscosity,
                                      cs_real_t       g_viscosity,
                                      cs_real_t       h_molar_mass,
                                      cs_real_t       ref_temperature)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));
  if (gw->model != CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid model.\n"
              "%s: One expects an immiscible two-phase flow model.\n",
              __func__, __func__);

  cs_gwf_tpf_t *mc = (cs_gwf_tpf_t *)gw->model_context;

  if (mc->is_miscible == true)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid model.\n"
              "%s: One expects an immiscible two-phase flow model.\n",
              __func__, __func__);

  assert(mc != nullptr);
  assert(l_mass_density > 0);
  assert(ref_temperature > 0);  /* In Kelvin */
  assert(h_molar_mass > 0);
  assert(l_viscosity > 0 && g_viscosity > 0);

  /* Set the parameters */

  mc->l_mass_density = l_mass_density;
  mc->l_viscosity = l_viscosity;
  mc->g_viscosity = g_viscosity;
  mc->l_diffusivity_h = 0;      /* immiscible case */
  mc->h_molar_mass = h_molar_mass;
  mc->ref_temperature = ref_temperature;
  mc->henry_constant = 1e-20;   /* immiscible case */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the flag dedicated to the post-processing of the GWF module
 *
 * \param[in] post_flag             flag to set
 * \param[in] reset                 reset post flag before
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_set_post_options(cs_flag_t       post_flag,
                        bool            reset)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  if (gw == nullptr)
    return;

  if (reset)
    gw->post_flag = post_flag;
  else
    gw->post_flag |= post_flag;

  if (gw->post_flag & CS_GWF_POST_DARCY_FLUX_AT_BOUNDARY) {

    cs_adv_field_t  *adv = _get_l_adv_field(gw);
    if (adv != nullptr)
      adv->status |= CS_ADVECTION_FIELD_DEFINE_AT_BOUNDARY_FACES;

    if (gw->model == CS_GWF_MODEL_MISCIBLE_TWO_PHASE ||
        gw->model == CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE) {

      cs_gwf_tpf_t  *mc = (cs_gwf_tpf_t *)gw->model_context;

      if (mc->g_darcy != nullptr) {
        adv = mc->g_darcy->adv_field;
        if (adv != nullptr)
          adv->status |= CS_ADVECTION_FIELD_DEFINE_AT_BOUNDARY_FACES;
      }
    }

  } /* Darcy flux at boundary */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the advection field related to the Darcy flux in the liquid
 *        phase
 *
 * \return a pointer to a cs_adv_field_t structure or nullptr
 */
/*----------------------------------------------------------------------------*/

cs_adv_field_t *
cs_gwf_get_adv_field(void)
{
  cs_gwf_t *gw = cs_gwf_main_structure;

  if (gw == nullptr)
    return nullptr;

  return _get_l_adv_field(gw);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and add a new cs_gwf_soil_t structure. An initialization by
 *        default of all members is performed.
 *        Case of a soil with an isotropic absolute permeability
 *
 * \param[in] z_name      name of the volume zone corresponding to the soil
 * \param[in] density     value of the bulk mass density
 * \param[in] k_abs       absolute (or intrisic) permeability (scalar-valued)
 * \param[in] porosity    value of the porosity (saturated moisture content)
 * \param[in] model       type of model for the soil behavior
 *
 * \return a pointer to the new allocated soil structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_soil_t *
cs_gwf_add_iso_soil(const char         *z_name,
                    double              density,
                    double              k_abs,
                    double              porosity,
                    cs_gwf_soil_model_t model)
{
  cs_gwf_t *gw = cs_gwf_main_structure;

  if (gw == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  const cs_zone_t  *zone = cs_volume_zone_by_name_try(z_name);

  if (zone == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Zone named \"%s\" is not defined.\n"
              " Stop adding a new soil.", __func__, z_name);

  assert(density > 0);
  assert(k_abs > 0);
  assert(porosity > 0);

  double  ktens_abs[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  ktens_abs[0][0] = ktens_abs[1][1] = ktens_abs[2][2] = k_abs;

  cs_gwf_soil_t  *soil = cs_gwf_soil_create(zone,
                                            gw->model, /* hydraulic model */
                                            model,     /* soil model */
                                            CS_PROPERTY_ISO,
                                            ktens_abs,
                                            porosity,
                                            density,
                                            gw->model_context);

  return soil;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and add a new cs_gwf_soil_t structure. An initialization by
 *        default of all members is performed.
 *
 * \param[in] z_name      name of the volume zone corresponding to the soil
 * \param[in] density     value of the bulk mass density
 * \param[in] k_abs       absolute (or intrisic) permeability (tensor-valued)
 * \param[in] porosity    value of the porosity (saturated moisture content)
 * \param[in] model       type of model for the soil behavior
 *
 * \return a pointer to the new allocated soil structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_soil_t *
cs_gwf_add_aniso_soil(const char                *z_name,
                      double                     density,
                      double                     k_abs[3][3],
                      double                     porosity,
                      cs_gwf_soil_model_t        model)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  const cs_zone_t  *zone = cs_volume_zone_by_name_try(z_name);

  if (zone == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Zone named \"%s\" is not defined.\n"
              " Stop adding a new soil.", __func__, z_name);

  assert(density > 0);
  assert(porosity > 0);
  assert(k_abs[0][0]*k_abs[0][0] +
         k_abs[1][1]*k_abs[1][1] +
         k_abs[2][2]*k_abs[2][2] > 0);

  cs_gwf_soil_t  *soil = cs_gwf_soil_create(zone,
                                            gw->model, /* hydraulic model */
                                            model,     /* soil model */
                                            CS_PROPERTY_ANISO,
                                            k_abs,
                                            porosity,
                                            density,
                                            gw->model_context);

  return soil;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new equation related to the groundwater flow module

 *         This equation is a particular type of unsteady advection-diffusion
 *         equation. The tracer is advected thanks to the darcian velocity and
 *         the diffusion property results from a physical modelling. Terms
 *         solved in this equation are activated according to predefined
 *         settings. The advection field corresponds to that of the liquid
 *         phase.
 *
 * \param[in]  tr_model   physical modelling to consider (0 = default settings)
 * \param[in]  eq_name    name of the tracer equation
 * \param[in]  var_name   name of the related variable
 *
 * \return a pointer to the new cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_t *
cs_gwf_add_tracer(cs_gwf_tracer_model_t     tr_model,
                  const char               *eq_name,
                  const char               *var_name)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  if (tr_model & CS_GWF_TRACER_USER)
    bft_error(__FILE__, __LINE__, 0,
              "%s: User-defined is not allowed in this context.\n"
              " Please consider cs_gwf_add_user_tracer() instead.", __func__);

  /* Set the advection field structure */

  cs_adv_field_t  *adv = _get_l_adv_field(gw);

  /* Set the function pointers */

  cs_gwf_tracer_init_setup_t  *init_setup = cs_gwf_tracer_default_init_setup;
  cs_gwf_tracer_finalize_setup_t *finalize_setup = nullptr;

  if (gw->model == CS_GWF_MODEL_SATURATED_SINGLE_PHASE)
    finalize_setup = cs_gwf_tracer_sat_finalize_setup;
  else
    finalize_setup = cs_gwf_tracer_unsat_finalize_setup;

  /* Call the main function to add a new tracer */

  assert(finalize_setup != nullptr);
  cs_gwf_tracer_t  *tracer = cs_gwf_tracer_add(tr_model,
                                               gw->model,
                                               eq_name,
                                               var_name,
                                               adv,
                                               0.,  /* lambda */
                                               -1,  /* chain position */
                                               -1,  /* chain id */
                                               init_setup,
                                               finalize_setup);

  return tracer;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new equation related to the groundwater flow module

 *         This equation is a particular type of unsteady advection-diffusion
 *         reaction equation. The tracer is advected thanks to the darcian
 *         velocity. The diffusion and reaction properties result from
 *         predefined physical modelling given by the parameter "tr_model".
 *         Other terms solved in this equation are activated according to
 *         predefined settings. The advection field corresponds to that of the
 *         liquid phase.
 *
 * \param[in]  tr_model   physical modelling to consider (0 = default settings)
 * \param[in]  eq_name    name of the tracer equation
 * \param[in]  var_name   name of the related variable
 * \param[in]  lambda     first order radioactive decay coefficient
 *
 * \return a pointer to the new cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_t *
cs_gwf_add_radioactive_tracer(cs_gwf_tracer_model_t     tr_model,
                              const char               *eq_name,
                              const char               *var_name,
                              double                    lambda)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  if (tr_model & CS_GWF_TRACER_USER)
    bft_error(__FILE__, __LINE__, 0,
              "%s: User-defined is not allowed in this context.\n"
              " Please consider cs_gwf_add_user_tracer() instead.", __func__);

  /* Set the advection field structure */

  cs_adv_field_t  *adv = _get_l_adv_field(gw);

  /* Set the function pointers */

  cs_gwf_tracer_init_setup_t  *init_setup = cs_gwf_tracer_default_init_setup;
  cs_gwf_tracer_finalize_setup_t *finalize_setup = nullptr;

  if (gw->model == CS_GWF_MODEL_SATURATED_SINGLE_PHASE)
    finalize_setup = cs_gwf_tracer_sat_finalize_setup;
  else
    finalize_setup = cs_gwf_tracer_unsat_finalize_setup;

  /* Call the main function to add a new tracer */

  cs_gwf_tracer_t  *tracer = cs_gwf_tracer_add(tr_model,
                                               gw->model,
                                               eq_name,
                                               var_name,
                                               adv,
                                               lambda,
                                               -1, /* chain position */
                                               -1, /* chain id */
                                               init_setup,
                                               finalize_setup);

  return tracer;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new equation related to the groundwater flow module
 *
 *         This equation is a particular type of unsteady advection-diffusion
 *         reaction equation.  Tracer is advected thanks to the darcian
 *         velocity and diffusion/reaction parameters result from a physical
 *         modelling. Terms are activated according to predefined settings.
 *         Modelling of the tracer parameters are left to the user
 *
 * \param[in] eq_name         name of the tracer equation
 * \param[in] var_name        name of the related variable
 * \param[in] init_setup      function pointer (predefined prototype)
 * \param[in] finalize_setup  function pointer (predefined prototype)
 *
 * \return a pointer to the new cs_gwf_tracer_t structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_t *
cs_gwf_add_user_tracer(const char                       *eq_name,
                       const char                       *var_name,
                       cs_gwf_tracer_init_setup_t       *init_setup,
                       cs_gwf_tracer_finalize_setup_t   *finalize_setup)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  /* Set the advection field structure */

  cs_adv_field_t  *adv = _get_l_adv_field(gw);

  /* Call the main function to add a new tracer */

  cs_gwf_tracer_t  *tracer = cs_gwf_tracer_add(CS_GWF_TRACER_USER,
                                               gw->model,
                                               eq_name,
                                               var_name,
                                               adv,
                                               0., /* not relevant here */
                                               -1, /* not relevant here */
                                               -1, /* not relevant here */
                                               init_setup,
                                               finalize_setup);

  return tracer;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a set of tracer equations corresponding to a radioactive decay
 *        chain in the groundwater flow module

 *        This equation is a particular type of unsteady advection-diffusion
 *        reaction equation. Tracer is advected thanks to the darcian velocity
 *        and diffusion/reaction parameters result from a physical modelling.
 *        Terms solved in this equation are activated according to predefined
 *        settings. The advection field corresponds to that of the liquid
 *        phase. A difference w.r.t. to standard tracer is the definition of
 *        specific source term taking into account the source/sink of the
 *        parent/current equation.
 *
 * \param[in] n_tracers    number of tracers equations
 * \param[in] unit         type of unit used in the tracer equations
 * \param[in] chain_name   name of the decay chain
 * \param[in] var_names    array of names of the related variable
 * \param[in] models       model associated to each tracer equation
 * \param[in] lambda_vals  set of first order radiactive decay coefficient
 *
 * \return a pointer to the new cs_gwf_tracer_decay_chain_t structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_tracer_decay_chain_t *
cs_gwf_add_decay_chain(int                       n_tracers,
                       cs_gwf_tracer_unit_t      unit,
                       const char               *chain_name,
                       const char               *var_names[],
                       cs_gwf_tracer_model_t     models[],
                       double                    lambda_vals[])
{
  cs_gwf_tracer_decay_chain_t *tdc = nullptr;

  if (n_tracers < 1)
    return tdc;

  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));
  assert(var_names != nullptr && models != nullptr && lambda_vals != nullptr);

  /* Create a structure for the decay chain */

  tdc = cs_gwf_tracer_create_decay_chain(n_tracers, chain_name, unit);

  /* Set the advection field structure */

  cs_adv_field_t  *adv = _get_l_adv_field(gw);

  /* Set the function pointers */

  cs_gwf_tracer_init_setup_t  *init_setup = cs_gwf_tracer_default_init_setup;
  cs_gwf_tracer_finalize_setup_t *finalize_setup = nullptr;

  if (gw->model == CS_GWF_MODEL_SATURATED_SINGLE_PHASE)
    finalize_setup = cs_gwf_tracer_sat_finalize_setup;
  else
    finalize_setup = cs_gwf_tracer_unsat_finalize_setup;

  size_t  max_len = strlen(var_names[0]);
  for (int i = 1; i < n_tracers; i++)
    if (strlen(var_names[i]) > max_len)
      max_len = strlen(var_names[i]);
  max_len += strlen("DecayChain%02d_") + 1;

  char *eqname = nullptr;
  BFT_MALLOC(eqname, max_len, char);

  for (int i = 0; i < n_tracers; i++) {

    /* Define automatically the equation name */

    const char  *varname = var_names[i];

    sprintf(eqname, "DecayChain%02d_%s", i, varname);

    /* Check the validity of the input data */

    if (lambda_vals[i] < 0) {

      cs_base_warn(__FILE__, __LINE__);
      cs_log_printf(CS_LOG_WARNINGS,
                    " %s: The decay coefficient for the tracer \"%s\" has a"
                    " negative value (%6.4e).\n",
                    __func__, varname, lambda_vals[i]);

    }

    /* Call the main function to add a new tracer */

    tdc->tracers[i] = cs_gwf_tracer_add(models[i],
                                        gw->model,
                                        eqname,
                                        varname,
                                        adv,
                                        lambda_vals[i],
                                        i,       /* position in the chain */
                                        tdc->id, /* id of the related chain */
                                        init_setup,
                                        finalize_setup);

  }

  BFT_FREE(eqname);

  return tdc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set shared pointers to main domain members
 *
 * \param[in] cdoq    pointer to additional mesh quantities for CDO schemes
 * \param[in] connect pointer to additional mesh connectivities for CDO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_init_sharing(const cs_cdo_quantities_t    *cdoq,
                    const cs_cdo_connect_t       *connect)
{
  /* Assign static const pointers */

  cs_cdo_quant = cdoq;
  cs_cdo_connect = connect;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the context of the model after the activation of the
 *        module and make first settings of the model parameters (physical and
 *        numerical). At this stage, cs_user_parameters() has not been called
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_init_model_context(void)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  if (cs_gwf_get_n_soils() < 1)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Groundwater module is activated but no soil is defined.",
              __func__);

  int dim = cs_gwf_soil_get_permeability_max_dim();
  cs_property_type_t  perm_type = CS_PROPERTY_ISO;
  if (dim == 9)
    perm_type = CS_PROPERTY_ANISO;
  else if (dim == 6)
    perm_type = CS_PROPERTY_ANISO_SYM;
  else if (dim == 3)
    perm_type = CS_PROPERTY_ORTHO;

  /* Add the absolute (or intrisic) permeability property (used in the
   * definition of the diffusion term in the conservation equations and in the
   * definition of the Darcy flux).
   */

  gw->abs_permeability = cs_property_add("absolute_permeability", perm_type);

  /* Continue the setup of the model and create new fields */

  switch (gw->model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    cs_gwf_sspf_init(
      (cs_gwf_sspf_t *)gw->model_context, gw->abs_permeability, gw->flag);
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    cs_gwf_uspf_init((cs_gwf_uspf_t *)gw->model_context, perm_type);
    break;

  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE:
    cs_gwf_tpf_init((cs_gwf_tpf_t *)gw->model_context, perm_type);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for the GroundWater Flow module.\n",
              __func__);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Predefined settings for the groundwater flow model and its related
 *        equations.
 *
 *        At this stage, all soils have been defined and equation parameters
 *        are set (cs_user_parameters() has been called and settings
 *        performed).
 *
 *        Create new cs_field_t structures according to the setting.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_init_setup(void)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  if (gw == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  int  permeability_dim = cs_property_get_dim(gw->abs_permeability);

  /* Continue the setup of the model and create new fields */

  switch (gw->model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    cs_gwf_sspf_init_setup(gw->flag, (cs_gwf_sspf_t *)gw->model_context);
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    cs_gwf_uspf_init_setup(gw->flag,
                           gw->post_flag,
                           permeability_dim,
                           (cs_gwf_uspf_t *)gw->model_context);
    break;

  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE:
    cs_gwf_tpf_init_setup(gw->post_flag, (cs_gwf_tpf_t *)gw->model_context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for the GroundWater Flow module.\n",
              __func__);
  }

  /* Add default post-processing related to groundwater flow module */

  cs_post_add_time_mesh_dep_output(cs_gwf_extra_post, gw);

  /* Same step for the tracer equations associated to the GWF module */

  cs_gwf_tracer_init_setup();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Last initialization step of the groundwater flow module. At this
 *        stage, the mesh quantities are defined.
 *
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_finalize_setup(const cs_cdo_connect_t     *connect,
                      const cs_cdo_quantities_t  *quant)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;
  if (gw == nullptr)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_gw));

  /* Define the soil porosity and the absolute permeability from the soil
     definition */

  cs_gwf_soil_define_shared_properties(gw->abs_permeability,
                                       gw->soil_porosity);

  switch (gw->model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    cs_gwf_sspf_finalize_setup(
      connect, quant, (cs_gwf_sspf_t *)gw->model_context);
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    cs_gwf_uspf_finalize_setup(
      connect, quant, gw->flag, (cs_gwf_uspf_t *)gw->model_context);
    break;

  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE:
    cs_gwf_tpf_finalize_setup(
      connect, quant, gw->flag, (cs_gwf_tpf_t *)gw->model_context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for the GroundWater Flow module.\n",
              __func__);
  }

  /* Finalize the soil setup */

  cs_gwf_soil_finalize_setup(gw->model, gw->post_flag, quant->n_cells);

  /* Finalize the tracer setup */

  cs_gwf_tracer_finalize_setup(connect, quant);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the groundwater system related to the hydraulic model:
 *        pressure head, head in law, moisture content, darcian velocity, soil
 *        capacity or permeability if needed.
 *        Quantities related to tracer model are updated elsewhere.
 *
 * \param[in] mesh          pointer to a cs_mesh_t structure
 * \param[in] connect       pointer to a cs_cdo_connect_t structure
 * \param[in] quant         pointer to a cs_cdo_quantities_t structure
 * \param[in] ts            pointer to a cs_time_step_t structure
 * \param[in] update_flag   metadata associated to the status of the update
 *                          step to perform
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_hydraulic_update(const cs_mesh_t             *mesh,
                        const cs_cdo_connect_t      *connect,
                        const cs_cdo_quantities_t   *quant,
                        const cs_time_step_t        *ts,
                        cs_flag_t                    update_flag)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Groundwater module is not allocated.", __func__);

  switch (gw->model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    cs_gwf_sspf_update(mesh,
                       connect,
                       quant,
                       ts,
                       update_flag,
                       gw->flag,
                       (cs_gwf_sspf_t *)gw->model_context);
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    cs_gwf_uspf_update(mesh,
                       connect,
                       quant,
                       ts,
                       update_flag,
                       gw->flag,
                       (cs_gwf_uspf_t *)gw->model_context);
    break;

  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE:
    cs_gwf_tpf_update(mesh,
                      connect,
                      quant,
                      ts,
                      update_flag,
                      gw->flag,
                      (cs_gwf_tpf_t *)gw->model_context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for the GroundWater Flow module.\n",
              __func__);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the GWF module (done after all the setup phase and after
 *        the initialization of all equations)
 *        One sets an initial value to all quantities related to this module.
 *
 * \param[in] mesh       pointer to a cs_mesh_t structure
 * \param[in] connect    pointer to a cs_cdo_connect_t structure
 * \param[in] quant      pointer to a cs_cdo_quantities_t structure
 * \param[in] ts         pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_init_values(const cs_mesh_t             *mesh,
                   const cs_cdo_connect_t      *connect,
                   const cs_cdo_quantities_t   *quant,
                   const cs_time_step_t        *ts)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  cs_gwf_hydraulic_update(mesh, connect, quant, ts, CS_FLAG_INITIALIZATION);

  /* Update the quantities related to tracer which depends on the "hydraulic"
   * state (the saturation, the hydraulic head, the Darcy velocity, etc.). Some
   * quantities are updated on-the-fly with functions associated to properties
   * but others needs to compute the new field or array to perform the update.
   *
   * For now, only the diffusion property associated to each tracer equation is
   * considered since the Darcy velocity may have changed.
   */

  cs_gwf_tracer_update_diff_pty(ts, mesh, connect, quant);

  /* Further steps dedicated to each model */

  switch (gw->model) {

  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE:
    cs_gwf_tpf_init_values(connect, quant, (cs_gwf_tpf_t *)gw->model_context);
    break;

  default:
    /* CS_GWF_MODEL_SATURATED_SINGLE_PHASE
       CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE */
    break; /* Nothing else to do */
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the steady-state of the groundwater flows module.
 *         Nothing is done if all equations are unsteady.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      time_step  pointer to a cs_time_step_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_compute_steady_state(const cs_mesh_t              *mesh,
                            const cs_time_step_t         *time_step,
                            const cs_cdo_connect_t       *connect,
                            const cs_cdo_quantities_t    *cdoq)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  /* Compute the new "hydraulic" state */
  /* --------------------------------- */

  switch (gw->model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    cs_gwf_sspf_compute_steady_state(mesh,
                                     connect,
                                     cdoq,
                                     time_step,
                                     gw->flag,
                                     (cs_gwf_sspf_t *)gw->model_context);
    break;

  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE:
    break; /* Nothing to do (the w_eq can be steady according to the numerical
              choices but this resolution is performed elsewhere) */

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for steady-state computations with the"
              " \"GroundWater Flow module\".\n",
              __func__);
  }

  /* Update the quantities related to tracer which depends on the "hydraulic"
   * state (the saturation, the hydraulic head, the Darcy velocity,
   * etc.). Some quantities are updated on-the-fly with functions associated
   * to properties but others needs to compute the new field or array to
   * perform the update.
   *
   * For now, only the diffusion property associated to each tracer equation
   * is considered since the Darcy velocity may have changed.
   */

  cs_gwf_tracer_update_diff_pty(time_step, mesh, connect, cdoq);

  /* Solve the tracer equations */
  /* -------------------------- */

  /* Once tracers are computed, one updates related quantities (the quantities
     associated to the tracer itsel) for instance the precipitation or the
     source term in the decay chain */

  cs_gwf_tracer_compute_steady_all(mesh, time_step, connect, cdoq);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the system related to groundwater flows module
 *
 * \param[in] mesh       pointer to a cs_mesh_t structure
 * \param[in] time_step  pointer to a cs_time_step_t structure
 * \param[in] connect    pointer to a cs_cdo_connect_t structure
 * \param[in] cdoq       pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_compute(const cs_mesh_t              *mesh,
               const cs_time_step_t         *time_step,
               const cs_cdo_connect_t       *connect,
               const cs_cdo_quantities_t    *cdoq)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  /* Compute the new "hydraulic" state */
  /* --------------------------------- */

  switch (gw->model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    cs_gwf_sspf_compute(mesh,
                        connect,
                        cdoq,
                        time_step,
                        gw->flag,
                        (cs_gwf_sspf_t *)gw->model_context);
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    cs_gwf_uspf_compute(mesh,
                        connect,
                        cdoq,
                        time_step,
                        gw->flag,
                        (cs_gwf_uspf_t *)gw->model_context);
    break;

  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE:
    cs_gwf_tpf_compute(mesh,
                       connect,
                       cdoq,
                       time_step,
                       gw->flag,
                       (cs_gwf_tpf_t *)gw->model_context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for the GroundWater Flow module.\n",
              __func__);
  }

  /* Solve the tracer equations */
  /* -------------------------- */

  /* Update the quantities related to tracer which depends on the "hydraulic"
   * state (the saturation, the hydraulic head, the Darcy velocity,
   * etc.). Some quantities are updated on-the-fly with functions associated
   * to properties but others needs to compute the new field or array to
   * perform the update.
   *
   * For now, only the diffusion property associated to each tracer equation
   * is considered since the Darcy velocity may have changed.
   */

  cs_gwf_tracer_update_diff_pty(time_step, mesh, connect, cdoq);

  /* Once tracers are computed, one updates related quantities (the quantities
     associated to the tracer itself) for instance the precipitation or the
     source term in the decay chain */

  cs_gwf_tracer_compute_all(mesh, time_step, connect, cdoq);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Predefined extra-operations for the groundwater flow module
 *
 * \param[in] connect   pointer to a cs_cdo_connect_t structure
 * \param[in] cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in] ts        pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_extra_op(const cs_cdo_connect_t      *connect,
                const cs_cdo_quantities_t   *cdoq,
                const cs_time_step_t        *ts)
{
  cs_gwf_t  *gw = cs_gwf_main_structure;

  if (gw == nullptr)
    return;

  switch (gw->model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    cs_gwf_sspf_extra_op(
      connect, cdoq, gw->post_flag, (cs_gwf_sspf_t *)gw->model_context);
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    cs_gwf_uspf_extra_op(
      connect, cdoq, gw->post_flag, (cs_gwf_uspf_t *)gw->model_context);
    break;

  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE:
    cs_gwf_tpf_extra_op(
      connect, cdoq, ts, gw->post_flag, (cs_gwf_tpf_t *)gw->model_context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for the GroundWater Flow module.\n",
              __func__);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Predefined post-processing output for the groundwater flow module.
 *        According to the model, additional postprocessing may be defined.
 *        Prototype of this function is given since it is a function pointer
 *        defined in cs_post.h (\ref cs_post_time_mesh_dep_output_t)
 *
 * \param[in, out] input        pointer to a optional structure (here a
 *                              cs_gwf_t structure)
 * \param[in]      mesh_id      id of the output mesh for the current call
 * \param[in]      cat_id       category id of the output mesh for this call
 * \param[in]      ent_flag     indicate global presence of cells (ent_flag[0]),
 *                              interior faces (ent_flag[1]), boundary faces
 *                              (ent_flag[2]), particles (ent_flag[3]) or probes
 *                              (ent_flag[4])
 * \param[in]      n_cells      local number of cells of post_mesh
 * \param[in]      n_i_faces    local number of interior faces of post_mesh
 * \param[in]      n_b_faces    local number of boundary faces of post_mesh
 * \param[in]      cell_ids     list of cells (0 to n-1)
 * \param[in]      i_face_ids   list of interior faces (0 to n-1)
 * \param[in]      b_face_ids   list of boundary faces (0 to n-1)
 * \param[in]      time_step    pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_extra_post(void                   *input,
                  int                     mesh_id,
                  int                     cat_id,
                  int                     ent_flag[5],
                  cs_lnum_t               n_cells,
                  cs_lnum_t               n_i_faces,
                  cs_lnum_t               n_b_faces,
                  const cs_lnum_t         cell_ids[],
                  const cs_lnum_t         i_face_ids[],
                  const cs_lnum_t         b_face_ids[],
                  const cs_time_step_t   *time_step)
{
  CS_UNUSED(cat_id);
  CS_UNUSED(ent_flag);
  CS_UNUSED(n_i_faces);
  CS_UNUSED(n_b_faces);
  CS_UNUSED(i_face_ids);
  CS_UNUSED(b_face_ids);

  if (input == nullptr)
    return;

  const cs_gwf_t  *gw = (const cs_gwf_t *)input;

  if (mesh_id == CS_POST_MESH_VOLUME) {

    if (gw->post_flag & CS_GWF_POST_DARCY_FLUX_DIVERGENCE) {

      cs_adv_field_t *adv = _get_l_adv_field(gw);

      /* Only case avalaible up to now */

      if (cs_advection_field_get_deftype(adv) == CS_XDEF_BY_ARRAY) {

        cs_real_t  *divergence =
          cs_advection_field_divergence_at_vertices(adv, time_step->t_cur);

        cs_post_write_vertex_var(mesh_id,
                                 CS_POST_WRITER_DEFAULT,
                                 "darcy_flux_divergence",
                                 1,
                                 false,
                                 false,
                                 CS_POST_TYPE_cs_real_t,
                                 divergence,
                                 time_step);

        BFT_FREE(divergence);

      }

    } /* Post-processing of the divergence is requested */

  }

  /* Additional post-processings */

  switch (gw->model) {

  case CS_GWF_MODEL_SATURATED_SINGLE_PHASE:
    cs_gwf_sspf_extra_post(mesh_id,
                           n_cells,
                           cell_ids,
                           gw->post_flag,
                           gw->abs_permeability,
                           (const cs_gwf_sspf_t *)gw->model_context,
                           time_step);
    break;

  case CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE:
    cs_gwf_uspf_extra_post(mesh_id,
                           n_cells,
                           cell_ids,
                           gw->post_flag,
                           (const cs_gwf_uspf_t *)gw->model_context,
                           time_step);
    break;

  case CS_GWF_MODEL_MISCIBLE_TWO_PHASE:
  case CS_GWF_MODEL_IMMISCIBLE_TWO_PHASE:
    cs_gwf_tpf_extra_post(mesh_id,
                          n_cells,
                          cell_ids,
                          gw->post_flag,
                          gw->abs_permeability,
                          (const cs_gwf_tpf_t *)gw->model_context,
                          cs_cdo_connect,
                          cs_cdo_quant,
                          time_step);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid model type for the GroundWater Flow module.\n",
              __func__);
  }
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
