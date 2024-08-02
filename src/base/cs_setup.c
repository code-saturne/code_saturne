/*============================================================================
 * Setup computation based on provided user data and functions.
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
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_ale.h"
#include "cs_atmo.h"
#include "cs_at_data_assim.h"
#include "cs_cf_thermo.h"
#include "cs_coal_read_data.h"
#include "cs_domain_setup.h"
#include "cs_elec_model.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_function_default.h"
#include "cs_gui.h"
#include "cs_gui_boundary_conditions.h"
#include "cs_gui_mobile_mesh.h"
#include "cs_gui_output.h"
#include "cs_gui_radiative_transfer.h"
#include "cs_gui_specific_physics.h"
#include "cs_gwf.h"
#include "cs_ibm.h"
#include "cs_internal_coupling.h"
#include "cs_lagr.h"
#include "cs_lagr_options.h"
#include "cs_mobile_structures.h"
#include "cs_parameters.h"
#include "cs_parameters_check.h"
#include "cs_porous_model.h"
#include "cs_prototypes.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_pressure_correction.h"
#include "cs_rad_transfer.h"
#include "cs_thermal_model.h"
#include "cs_turbulence_model.h"
#include "cs_rad_transfer_options.h"
#include "cs_restart.h"
#include "cs_runaway_check.h"
#include "cs_velocity_pressure.h"
#include "cs_wall_condensation.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_setup.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------
 * Local type definitions
 *----------------------------------------------------------------------------*/

/*============================================================================
 * External function prototypes
 *============================================================================*/

/* Bindings to Fortran routines */

void
cs_f_addfld(void);

void
cs_f_steady_laminar_flamelet_read_base(void);

void
cs_f_usipes(int *nmodpp);

void
cs_f_impini(void);

void
cs_f_indsui(void);

void
cs_f_colecd(void);

void
cs_f_fldini(void);

void
cs_f_usppmo(void);

void
cs_f_fldvar(int *nmodpp);

void
cs_f_atini1(void);

void
cs_f_solcat(int iappel);

void
cs_f_fldprp(void);

void
cs_f_usipsu(int *nmodpp);

void
cs_f_varpos(void);

void
cs_f_lagran_init_map(void);

void
cs_f_iniini(void);

void
cs_f_ppinii(void);

void
cs_f_ppini1(void);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief First initialization stages
 */
/*----------------------------------------------------------------------------*/

static void
_init_setup(void)
{
  cs_log_printf
    (CS_LOG_DEFAULT,
     _("\n\n"
       "===============================================================\n\n\n"
       "                   CALCULATION PREPARATION\n"
       "                   =======================\n\n\n"
       "===============================================================\n\n\n"));

  /* File for some specific physical models */

  cs_atmo_set_meteo_file_name("meteo");

  /* Handle some reference and physical values */

  cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();
  fp->pther = fp->p0;

  /* Other mappings */

  cs_f_iniini();
  cs_f_ppinii();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read specific physical model data file
 */
/*----------------------------------------------------------------------------*/

static void
_read_specific_physics_data(void)
{
  /* Diffusion flame - 3-point chemistry
   * premix flame    - EBU model
   * premix flame    - LWC model */
  if (   cs_glob_physical_model_flag[CS_COMBUSTION_3PT] != -1
      || cs_glob_physical_model_flag[CS_COMBUSTION_EBU] != -1
      || cs_glob_physical_model_flag[CS_COMBUSTION_LW] != -1)
    cs_f_colecd();

  /* Diffusion flame - steady laminar flamelet approach */
  if (cs_glob_physical_model_flag[CS_COMBUSTION_SLFM] != -1)
    cs_f_steady_laminar_flamelet_read_base();

  /* Pulverized coal combustion */
  if (cs_glob_physical_model_flag[CS_COMBUSTION_COAL] != -1) {
    cs_gui_coal_model();
    cs_coal_read_data();
  }

  /* Joule effect, electric arc or ionic conduction */
  if (   cs_glob_physical_model_flag[CS_JOULE_EFFECT] >= 1
      || cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] >= 1)
    cs_electrical_model_param();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function calling the user-defined functions for the definition
 *        of computation parameters: we run through here in any case.
 */
/*----------------------------------------------------------------------------*/

static void
_init_user
(
  int   *nmodpp
)
{
  cs_fluid_properties_t *fluid_props = cs_get_glob_fluid_properties();

  int icondb = cs_glob_wall_condensation->icondb;
  int condv = cs_glob_wall_condensation->icondv;

  /* Check for restart and read matching time steps and notebook values */
  cs_parameters_read_restart_info();

  /* Flow model selection through GUI */
  cs_gui_physical_model_select();

  /* Flow model selection through user Fortran subroutine */
  /* Warning : This is a user function maybe no need to translate it...
   * It is only called in the case 69_PARISFOG */
  cs_f_usppmo();
  cs_wall_condensation_set_onoff_state(icondb, condv);

  /* Other model selection through GUI */
  /* ALE parameters */
  cs_gui_ale_params();

  /* Thermal model */
  cs_gui_thermal_model();

  /* Turbulence model */
  cs_gui_turb_model();

  cs_gui_cp_params();
  cs_gui_dt();
  cs_gui_hydrostatic_pressure();

  /* Gravity and Coriolis
   * Presence or not of gravity may be needed to determine whether some fields
   * are created, so this is called before cs_user_model (to provide a
   * user call site allowing to modify GUI-defined value programatically
   * before property fields are created). */
  cs_gui_physical_constants();

  /* Activate radiative transfer model */
  cs_gui_radiative_transfer_parameters();
  cs_user_radiative_transfer_parameters();

  /* Flow and other model selection through user C routines */
  cs_user_model();

  /* Set type and order of the turbulence model */
  cs_set_type_order_turbulence_model();

  /* If CDO is active, initialize the context structures for models which
   * have been activated */
  if (cs_glob_param_cdo_mode != CS_PARAM_CDO_MODE_OFF) {
    /* Groundwater flow module */
    if (cs_gwf_is_activated())
      cs_gwf_init_model_context();
  }

  /* Activate CDO for ALE */
  if (cs_glob_ale == CS_ALE_CDO)
    cs_ale_activate();

  if (cs_glob_ale == CS_ALE_LEGACY)
    cs_gui_mobile_mesh_structures_add();

  /* Read thermomechanical data for specific physics */
  _read_specific_physics_data();

  /* Other model parameters, including user-defined scalars */
  cs_gui_user_variables();
  cs_gui_user_arrays();
  cs_gui_calculator_functions();

  /* Solid zones */
  cs_velocity_pressure_set_solid();

  /* Initialize parameters for specific physics */
  cs_rad_transfer_options();

  if (cs_glob_param_cdo_mode != CS_PARAM_CDO_MODE_ONLY) {
    cs_f_fldvar(nmodpp);

    /* Activate pressure correction model if CDO mode is not stand-alone */
    cs_pressure_correction_model_activate();
  }

  if (cs_glob_ale != CS_ALE_NONE)
    cs_gui_ale_diffusion_type();

  cs_gui_laminar_viscosity();

  /* Specific physics modules
   * Note: part of what is inside ppini1 could be moved here
   * so that usipsu / cs_user_parameters can be used by the user
   * to modify default settings */

  /* Atmospheric flows */
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] != -1)
    cs_f_atini1();

  /* Compressible flows */
  cs_gui_hydrostatic_equ_param();
  const cs_field_t *f_id = cs_field_by_name_try("velocity");
  if (f_id != NULL) {
    if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] != -1)
      cs_runaway_check_define_field_max(f_id->id, 1.e5);
    else
      cs_runaway_check_define_field_max(f_id->id, 1.e4);
  }

  /* Atmospheric module */
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] != -1) {
    /* Advanced init/allocation for the soil model */
    if (cs_glob_atmo_option->soil_cat >= 0)
      cs_f_solcat(1);
  }

  /* Initialization of global parameters */
  cs_gui_output_boundary();

  if (cs_glob_param_cdo_mode != CS_PARAM_CDO_MODE_ONLY)
    cs_f_fldprp();

  /* Initialization of additional user parameters */
  cs_gui_checkpoint_parameters();

  cs_gui_dt_param();

  /* Local numerical options */
  cs_gui_equation_parameters();

  if (cs_glob_param_cdo_mode != CS_PARAM_CDO_MODE_ONLY)
    cs_gui_numerical_options();

  /* Physical properties */
  cs_gui_physical_properties();

  /* Turbulence reference values (uref, almax) */
  cs_gui_turb_ref_values();

  /* Set turbulence constants according to model choices.
   * This can be overwritten by the user in cs_user_parameters() */
  cs_turb_compute_constants(-1);

  /* Scamin, scamax, turbulent flux model, diffusivities
   * (may change physical properties for some scalars, so called
   * after cs_gui_physical_properties). */
  cs_gui_scalar_model_settings();

  /* Porosity model */
  cs_gui_porous_model();

  /* Init fan */
  cs_gui_define_fans();

  /* Init error estimator */
  cs_gui_error_estimator();

  /* Initialize base evaluation functions */
  cs_function_default_define();

  /* User functions */
  cs_f_usipsu(nmodpp);
  cs_user_parameters(cs_glob_domain);

  /* If time step is local or variable, pass information to C layer, as it
   * may be needed for some field (or moment) definitions */
  if (cs_glob_time_step_options->idtvar != 0)
    cs_time_step_define_variable(1);

  if (cs_glob_time_step_options->idtvar == 2
    || cs_glob_time_step_options->idtvar == -1)
    cs_time_step_define_local(1);

  /* Initialize Fortran restarted computation flag (isuite) */
  cs_f_indsui();

  /* Default value of physical properties for the compressible model */
  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] != -1) {
    /* EOS has been set above with the GUI or in cs_user_model.
     * The variability of the thermal conductivity
     * (diffusivity_id for itempk) and the volume viscosity (iviscv) has
     * been set in fldprp.
     *
     * Compute cv0 according to chosen EOS */

    cs_real_t l_cp[1] = {fluid_props->cp0};
    cs_real_t l_xmasmr[1] = {fluid_props->xmasmr};
    cs_real_t l_cv[1] = {-1};
    cs_cf_thermo_cv(l_cp, l_xmasmr, l_cv, 1);
    fluid_props->cv0 = l_cv[0];
  }

  if (cs_glob_porosity_ibm_opt->porosity_mode > 0)
    cs_porous_model_set_model(3);

  /* Varpos
   * If CDO mode only, skip this stage */
  if (cs_glob_param_cdo_mode != CS_PARAM_CDO_MODE_ONLY)
    cs_f_varpos();

  /* Internal coupling */
  cs_gui_internal_coupling();
  cs_user_internal_coupling();

  cs_internal_coupling_setup();

  /* Mobile structures
   * After call to cs_gui_mobile_mesh_structures_add possible
   * call by user to cs_mobile_structures_add_n_structures */
  if (cs_glob_ale != CS_ALE_NONE)
    cs_mobile_structures_setup();
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup computation based on provided user data and functions.
 */
/*----------------------------------------------------------------------------*/

void
cs_setup(void)
{
  int is_restart = cs_restart_present();
  cs_real_t dtref = cs_glob_time_step->dt[0];
  cs_time_scheme_t *time_scheme = cs_get_glob_time_scheme();

  /* Initialize modules before user has access */
  _init_setup();

  int nmodpp = 0;
  for (int i = 1; i < CS_N_PHYSICAL_MODEL_TYPES; i++) {
    if (cs_glob_physical_model_flag[i] > -1)
      nmodpp ++;
  }

  /* User input, variable definitions */
  _init_user(&nmodpp);

  cs_f_ppini1();

  /* Map Fortran pointers to C global data */
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] != -1)
    cs_at_data_assim_initialize();

  /* Initialize lagr structures */
  cs_f_lagran_init_map();
  cs_lagr_map_specific_physics();

  int have_thermal_model = 0;
  if (cs_thermal_model_field() != NULL)
    have_thermal_model = 1;

  cs_lagr_options_definition(is_restart,
                             have_thermal_model,
                             dtref,
                             &time_scheme->iccvfg);
  cs_lagr_add_fields();

  if (cs_glob_param_cdo_mode != CS_PARAM_CDO_MODE_ONLY) {
    /* Additional fields if not in CDO mode only */
    cs_f_addfld();

    /* Changes after user initialization and additional fields dependent on
     * main fields options. */
    cs_parameters_global_complete();

    cs_f_fldini();
  }

  cs_parameters_eqp_complete();

  /* Time moments called after additional creation */
  cs_gui_time_moments();
  cs_user_time_moments();

  /* GUI based boundary condition definitions */
  /* TODO change NULL in nullptr */
  cs_gui_boundary_conditions_define(NULL);

  /* Some final settings */
  cs_gui_output();

  if (cs_glob_param_cdo_mode != CS_PARAM_CDO_MODE_ONLY) {
    /* Warning: Used in 0 validation cases ? */
    cs_f_usipes(&nmodpp);

    /* Avoid a second spurious call to this function
     * called in the C part if CDO is activated */
    if (cs_glob_param_cdo_mode == CS_PARAM_CDO_MODE_OFF) {
      cs_user_boundary_conditions_setup(cs_glob_domain);
      cs_user_finalize_setup(cs_glob_domain);
    }
  }

  cs_parameters_output_complete();

  /* Coherency checks */
  if (cs_glob_param_cdo_mode != CS_PARAM_CDO_MODE_ONLY)
    cs_parameters_check();

  cs_log_printf(CS_LOG_DEFAULT,
                "\n"
                "No error detected during the data verification.\n");

  /* Print output */
  cs_f_impini();
}

/*-----------------------------------------------------------------------------*/

END_C_DECLS
