/*============================================================================
 * Gas combustion model.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_mem.h"
#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_assert.h"
#include "base/cs_array.h"
#include "base/cs_base.h"
#include "base/cs_field_default.h"
#include "base/cs_field_pointer.h"
#include "base/cs_parameters_check.h"
#include "base/cs_physical_constants.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_physical_properties.h"
#include "base/cs_prototypes.h"
#include "base/cs_restart_default.h"
#include "base/cs_thermal_model.h"
#include "mesh/cs_mesh_location.h"
#include "turb/cs_turbulence_model.h"

#include "gui/cs_gui_specific_physics.h"
#include "pprt/cs_combustion_model.h"
#include "pprt/cs_physical_model.h"
#include "rayt/cs_rad_transfer.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cogz/cs_combustion_gas.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_combustion_gas.cpp
        Gas combustion model.
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

static const char *
_model_name(cs_combustion_gas_model_type_t  model_type)
{
  switch(model_type) {
  case 100:
    return "CS_COMBUSTION_3PT_ADIABATIC";
    break;
  case 101:
    return "CS_COMBUSTION_3PT_PERMEATIC";
    break;
  case 102:
    return "CS_COMBUSTION_BSH_ADIABATIC";
    break;
  case 103:
    return "CS_COMBUSTION_BSH_PERMEATIC";
    break;
  case 200:
    return "CS_COMBUSTION_SLFM_STEADY_ADIABATIC";
    break;
  case 201:
    return "CS_COMBUSTION_SLFM_STEADY_PERMEATIC";
    break;
  case 203:
    return "CS_COMBUSTION_SLFM_PROGRESS_PERMEATIC";
    break;
  case 300:
    return "CS_COMBUSTION_EBU_CONSTANT_ADIABATIC";
    break;
  case 301:
    return "CS_COMBUSTION_EBU_CONSTANT_PERMEATIC";
    break;
  case 302:
    return "CS_COMBUSTION_EBU_VARIABLE_ADIABATIC";
    break;
  case 303:
    return "CS_COMBUSTION_EBU_VARIABLE_PERMEATIC";
    break;
  case 400:
    return "CS_COMBUSTION_LW_2PEAK_ADIABATIC";
    break;
  case 401:
    return "CS_COMBUSTION_LW_2PEAK_PERMEATIC";
    break;
  case 402:
    return "CS_COMBUSTION_LW_3PEAK_ADIABATIC";
    break;
  case 403:
    return "CS_COMBUSTION_LW_3PEAK_PERMEATIC";
    break;
  case 404:
    return "CS_COMBUSTION_LW_4PEAK_ADIABATIC";
    break;
  case 405:
    return "CS_COMBUSTION_LW_4PEAK_PERMEATIC";
    break;
  default:
    break;
  }
  return "?";
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Global variables
 *============================================================================*/

/*! Combustion model parameters structure */

cs_combustion_gas_model_t  *cs_glob_combustion_gas_model = NULL;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize gas combustion model.
 */
/*----------------------------------------------------------------------------*/

static void
_combustion_gas_finalize(void)
{
  if (cs_glob_combustion_gas_model != NULL) {

    cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

    CS_FREE(cm->data_file_name);
    CS_FREE(cm->radiation_data_file_name);
    CS_FREE(cm->flamelet_library);
    CS_FREE(cm->radiation_library);
    CS_FREE(cm->rho_library);
    CS_FREE(cm);

    cs_glob_combustion_gas_model = cm;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a gas combustion model field.
 *
 * \param[in]  name   field base name
 * \param[in]  label  field base label
 *
 * \return  pointer to field;
 */
/*----------------------------------------------------------------------------*/

static cs_field_t *
_add_model_variable(const char  *name,
                    const char  *label)
{
  int f_id = cs_variable_field_create(name, label, CS_MESH_LOCATION_CELLS, 1);
  cs_field_t *f = cs_field_by_id(f_id);

  cs_add_model_field_indexes(f->id);

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a 1d property field at cells.
 *
 * \param[in]  name        field name
 * \param[in]  label       field label, or null
 *
 * \return  pointer to field
 */
/*----------------------------------------------------------------------------*/

static cs_field_t *
_add_property_1d(const char  *name,
                 const char  *label)
{
  const int keyvis = cs_field_key_id("post_vis");
  const int keylog = cs_field_key_id("log");
  const int keylbl = cs_field_key_id("label");

  if (cs_field_by_name_try(name) != NULL)
    cs_parameters_error(CS_ABORT_IMMEDIATE,
                        _("initial data setup"),
                        _("Field %s has already been assigned.\n"),
                        name);

  cs_physical_property_define_from_field(name,
                                         CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY,
                                         CS_MESH_LOCATION_CELLS,
                                         1,       // dim
                                         false);  // has previous

  int f_id = cs_physical_property_field_id_by_name(name);
  cs_field_t *f = cs_field_by_id(f_id);

  cs_field_set_key_int(f, keyvis, 0);
  cs_field_set_key_int(f, keylog, 1);
  if (label != nullptr)
    cs_field_set_key_str(f, keylbl, label);

  return f;
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Activate gas combustion model.
 *
 * \param[in]  type  gas combustion model type
 *
 * \return  pointer to gas combustion model structure.
 */
/*----------------------------------------------------------------------------*/

cs_combustion_gas_model_t  *
cs_combustion_gas_set_model(cs_combustion_gas_model_type_t  type)
{
  cs_glob_physical_model_flag[CS_COMBUSTION_3PT] = -1;
  cs_glob_physical_model_flag[CS_COMBUSTION_SLFM] = -1;
  cs_glob_physical_model_flag[CS_COMBUSTION_EBU] = -1;
  cs_glob_physical_model_flag[CS_COMBUSTION_LW] = -1;

  if (type == CS_COMBUSTION_GAS_NONE) {
    CS_FREE(cs_glob_combustion_gas_model);
    return NULL;
  }

  int macro_type = type / 100;
  int sub_type = type % 100;
  if (macro_type == 1)
    cs_glob_physical_model_flag[CS_COMBUSTION_3PT] = sub_type;
  else if (macro_type == 2)
    cs_glob_physical_model_flag[CS_COMBUSTION_SLFM] = sub_type;
  else if (macro_type == 3)
    cs_glob_physical_model_flag[CS_COMBUSTION_EBU] = sub_type;
  else if (macro_type == 4)
    cs_glob_physical_model_flag[CS_COMBUSTION_LW] = sub_type;

  if (cs_glob_combustion_gas_model != NULL) {
    cs_glob_combustion_gas_model->type = type;
    return cs_glob_combustion_gas_model;
  }

  /* Create and initialize model structure */

  cs_combustion_gas_model_t  *cm;

  CS_MALLOC(cm, 1, cs_combustion_gas_model_t);
  memset(cm, 0, sizeof(cs_combustion_gas_model_t));

  cs_glob_combustion_gas_model = cm;

  /* Members also present in coal combustion model */

  cm->n_gas_el_comp = 0;
  cm->n_gas_species = 0;
  cm->n_atomic_species = 0;
  cm->n_reactions = 0;
  cm->n_tab_points = 0;
  cm->pcigas = 0;
  cm->xco2 = 0;
  cm->xh2o = 0;

  for (int i = 0; i < CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS; i++) {
    cm->wmole[i] = 0;
  }

  for (int i = 0; i < CS_COMBUSTION_GAS_MAX_OXYDANTS; i++) {
    cm->oxyn2[i] = 0;
    cm->oxyh2o[i] = 0;
    cm->oxyco2[i] = 0;
  }

  for (int i = 0; i < CS_COMBUSTION_GAS_MAX_TABULATION_POINTS; i++) {
    cm->th[i] = 0;
  }

  /* Members specific to gas combustion model */

  cm->type = type;

  cm->data_file_name = nullptr;
  cm->radiation_data_file_name = nullptr;

  cm->use_janaf = true;

  cm->iic = 0;
  cm->iico2 = 0;
  cm->iio2 = 0;
  cm->isoot = -1;

  cm->hinfue = cs_math_big_r;
  cm->tinfue = - cs_math_big_r;
  cm->tinoxy = - cs_math_big_r;
  cm->hinoxy = cs_math_big_r;

  cm->xsoot = 0.;
  cm->rosoot = 0.;
  cm->lsp_fuel = 0.;

  cm->frmel = 0.;
  cm->tgf = 300.;
  cm->cebu = 2.5;

  for (int i = 0; i < CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES; i++) {
    cm->wmolg[i] = 0;
  }

  for (int i = 0; i < CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES; i++) {
    for (int j = 0; j < CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS; j++) {
      cm->coefeg[i][j] = 0;
    }
  }

  for (int i = 0; i < CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES; i++) {
    for (int j = 0; j < CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS; j++) {
      cm->compog[i][j] = 0;
    }
  }

  for (int i = 0; i < CS_COMBUSTION_GAS_MAX_GLOBAL_REACTIONS; i++) {
    cm->fs[i] = 0;
  }

  for (int i = 0; i < CS_COMBUSTION_GAS_MAX_TABULATION_POINTS; i++) {
    for (int j = 0; j < CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES; j++) {
      cm->cp_gas_g[i][j] = 0;
    }
  }

  /* Fields */

  cm->fm = nullptr;
  cm->fp2m = nullptr;
  cm->fsqm = nullptr;
  cm->pvm = nullptr;
  cm->ygfm = nullptr;
  cm->yfm = nullptr;
  cm->yfp2m = nullptr;
  cm->coyfp = nullptr;
  cm->npm = nullptr;
  cm->fsm = nullptr;

  for (int i = 0; i < CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES; i++) {
    cm->ym[i] = nullptr;
    cm->bym[i] = nullptr;
  }

  cm->t2m = nullptr;

  /* Numerical parameters */

  cm->srrom = 0.95;

  /* Legacy 3-pt model */
  const int nmaxh = 9, nmaxf = 9;

  cm->hstoea = 0.;
  for (int i = 0; i < nmaxh; i++)
    cm->hh[i] = 0.;
  for (int i = 0; i < nmaxf; i++)
    cm->ff[i] = 0.;
  for (int i = 0; i < nmaxh; i++) {
    for (int j = 0; j < nmaxf; j++)
      cm->tfh[i][j] = 0.;
  }

  /* Burke-Schumann model */
  for (int i = 0; i < 7; i++) {
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 5; k++)
        cm->coeff_therm[i][j][k] = 0.;
    }
  }

  /* Libby Williams model */
  cm->lw.vref = - cs_math_big_r;
  cm->lw.lref = - cs_math_big_r;
  cm->lw.ta = - cs_math_big_r;
  cm->lw.tstar = - cs_math_big_r;
  cm->lw.fmin = 0.;
  cm->lw.fmax = 1.;
  cm->lw.hmin = 0.;
  cm->lw.hmax = 0.;
  cm->lw.coeff1 = 0.;
  cm->lw.coeff2 = 0.;
  cm->lw.coeff3 = 0.;

  // Number of Diracs for LWC model

  if (cm->type / 100 == CS_COMBUSTION_LW) {
    int lw_model = cm->type % 100;
    if (lw_model == 0 || lw_model == 1)
      cm->lw.n_dirac = 2;
    else if (lw_model == 2 || lw_model == 3)
      cm->lw.n_dirac = 3;
    else if (lw_model == 4 || lw_model == 5)
      cm->lw.n_dirac = 4;
  }

  /*! Steady flamelet model */

  cm->n_gas_fl = -1;
  cm->nki = -1;
  cm->nxr = -1;
  cm->nzm = -1;
  cm->nzvar = -1;
  cm->nlibvar = -1;
  cm->ikimid = 1;
  cm->mode_fp2m = 1;

  cm->flamelet_library = nullptr;
  cm->radiation_library = nullptr;
  cm->rho_library = nullptr;

  cm->flamelet_zm = -1;
  cm->flamelet_zvar = -1;;
  cm->flamelet_xr = -1;
  cm->flamelet_ki = -1;
  cm->flamelet_temp = -1;
  cm->flamelet_rho = -1;
  cm->flamelet_vis = -1;
  cm->flamelet_dt = -1;
  cm->flamelet_temp2 = -1;
  cm->flamelet_hrr = -1;
  cm->flamelet_c = -1;
  cm->flamelet_omg_c = -1;

  cs_array_int_set_value(CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES,
                         -1,
                         cm->flamelet_species);

  /* Set finalization callback */

  cs_base_at_finalize(_combustion_gas_finalize);

  return cm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the thermochemical data file name.
 *
 * \param[in] file_name      name of the file.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_gas_set_thermochemical_data_file(const char  *file_name)
{
  cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  if (cm == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: gas combustion model not active."), __func__);

  CS_FREE(cm->data_file_name);

  if (file_name != NULL) {
    size_t l = strlen(file_name)+1;
    CS_REALLOC(cm->data_file_name, l, char);
    strncpy(cm->data_file_name, file_name, l);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the thermochemical data file name.
 *
 * \param[in] file_name      name of the file.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_gas_set_radiation_data_file(const char  *file_name)
{
  cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  if (cm == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: gas combustion model not active."), __func__);

  CS_FREE(cm->radiation_data_file_name);

  if (file_name != NULL) {
    size_t l = strlen(file_name)+1;
    CS_REALLOC(cm->radiation_data_file_name, l, char);
    strncpy(cm->radiation_data_file_name, file_name, l);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize steady flamelet library
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_slfm_init_library(void)
{
  if (cs_glob_combustion_gas_model != NULL) {

    cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

    size_t n = cm->nlibvar * cm->nxr * cm->nki * cm->nzvar * cm->nzm;

    CS_MALLOC(cm->flamelet_library, n, cs_real_t);

    cs_real_t *_flamelet_library = cm->flamelet_library;

    cs_array_real_fill_zero(n, _flamelet_library);

    size_t n_rho_lib = 5 * cm->nxr * cm->nki * cm->nzvar * cm->nzm;

    CS_MALLOC(cm->rho_library, n_rho_lib, cs_real_t);
    cs_array_real_fill_zero(n_rho_lib, cm->rho_library);

    if (cs_glob_rad_transfer_params->type > CS_RAD_TRANSFER_NONE) {
      size_t n_rad_lib = 2 * cs_glob_rad_transfer_params->nwsgg
                       * cm->nxr * cm->nki * cm->nzvar * cm->nzm;

      CS_MALLOC(cm->radiation_library, n_rad_lib, cs_real_t);

      cs_array_real_fill_zero(n_rad_lib, cm->radiation_library);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Specific setup operations for gas combustion models.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_gas_setup(void)
{
  if (cs_glob_combustion_gas_model == NULL)
    return;

  cs_combustion_gas_model_t  *cm = cs_glob_combustion_gas_model;
  int macro_type = cm->type / 100;

  const int kclvfl = cs_field_key_id("variance_clipping");
  const int keysca  = cs_field_key_id("scalar_id");
  const int kscavr = cs_field_key_id("first_moment_id");
  const int ksigmas = cs_field_key_id("turbulent_schmidt");
  const int kvisl0 = cs_field_key_id("diffusivity_ref");
  const cs_real_t viscl0 = cs_glob_fluid_properties->viscl0;

  /* Transported variables
     --------------------- */

  if (macro_type == CS_COMBUSTION_3PT) {
    // Clipping to scamax of the mixture fraction variance.
    cs_field_set_key_int(cm->fp2m, kclvfl, 2);
  }

  else if (macro_type == CS_COMBUSTION_SLFM) {
    if (cm->mode_fp2m == 0) {
      // Clipping by max of fm
      cs_field_set_key_int(cm->fp2m, kclvfl, 1);
    }
  }

  else if (macro_type == CS_COMBUSTION_LW) {
    cs_field_set_key_int(cm->fp2m, kclvfl, 0);
    cs_field_set_key_int(cm->yfp2m, kclvfl, 0);
  }

  // Physical or numerical values specific to combustion scalars

  const int n_fields = cs_field_n_fields();
  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (!(f->type & CS_FIELD_VARIABLE))
      continue;
    if (f->type & CS_FIELD_CDO || f->type & CS_FIELD_USER)
      continue;
    if (cs_field_get_key_int(f, keysca) <= 0)
      continue;

    int variance_id = cs_field_get_key_int(f, kscavr);
    if (f != CS_F_(h) && variance_id <= 0) {
      // We consider that turbulent viscosity domintes. We prohibit
      // the computation of laminar flames with  =/= 1
      cs_field_set_key_double(f, kvisl0, viscl0);
    }

    // Turbulent Schmidt or Prandtl
    cs_field_set_key_double(f, ksigmas, 0.7);
  }

  /* Additional information
     ---------------------- */

  // Compute ro0 based on t0 and p0.

  cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();

  if (macro_type != CS_COMBUSTION_SLFM)
    fp->ro0 = fp->pther * cm->wmolg[1] / (cs_physical_constants_r * fp->t0);
  else
    fp->ro0 = cm->flamelet_library[cm->flamelet_rho];

  fp->roref = fp->ro0;

  // Variable density
  fp->irovar = 1;

  /* User settings
     -------------*/

  cs_gui_combustion_gas_model();

  if (   macro_type == CS_COMBUSTION_3PT
      || macro_type == CS_COMBUSTION_SLFM)
    cs_gui_combustion_gas_model_temperatures();

  /* Parameter checks
     ---------------- */

  const char section_name[] = N_("Gas combustion model setup");

  if (macro_type == CS_COMBUSTION_3PT) {
    cs_parameters_is_greater_double(CS_ABORT_DELAYED, _(section_name),
                                    "tinfue", cm->tinfue, 0.);
    cs_parameters_is_greater_double(CS_ABORT_DELAYED, _(section_name),
                                    "tinoxy", cm->tinoxy, 0.);
  }

  else if (macro_type == CS_COMBUSTION_SLFM) {
    cs_parameters_is_greater_int(CS_ABORT_DELAYED,
                                 _(section_name),
                                 "flamelet_zm", cm->flamelet_zm, -1);
    cs_parameters_is_greater_int(CS_ABORT_DELAYED,
                                 _(section_name),
                                 "flamelet_zvar", cm->flamelet_zvar, -1);
    cs_parameters_is_greater_int(CS_ABORT_DELAYED,
                                 _(section_name),
                                 "flamelet_xr", cm->flamelet_xr, -1);
    cs_parameters_is_greater_int(CS_ABORT_DELAYED,
                                 _(section_name),
                                 "flamelet_temp", cm->flamelet_temp, -1);
    cs_parameters_is_greater_int(CS_ABORT_DELAYED,
                                 _(section_name),
                                 "flamelet_rho", cm->flamelet_rho, -1);
    cs_parameters_is_greater_int(CS_ABORT_DELAYED,
                                 _(section_name),
                                 "flamelet_vis", cm->flamelet_vis, -1);
    cs_parameters_is_greater_int(CS_ABORT_DELAYED,
                                 _(section_name),
                                 "flamelet_dt", cm->flamelet_dt, -1);

    if (cm->type%100 < 2)
      cs_parameters_is_greater_int(CS_ABORT_DELAYED,
                                     _(section_name),
                                     "flamelet_ki", cm->flamelet_ki, -1);
    else {
      cs_parameters_is_greater_int(CS_ABORT_DELAYED,
                                     _(section_name),
                                     "flamelet_c", cm->flamelet_c, -1);
      cs_parameters_is_greater_int(CS_ABORT_DELAYED,
                                     _(section_name),
                                     "flamelet_omg_c", cm->flamelet_omg_c, -1);
    }

    if (cm->hinfue < - cs_math_big_r)
      cs_parameters_error
        (CS_ABORT_DELAYED, _(section_name),
         _("hinfue must be set by the user, and not remain at %g.\n"),
         cm->hinfue);
    if (cm->hinoxy < - cs_math_big_r)
      cs_parameters_error
        (CS_ABORT_DELAYED, _(section_name),
         _("hinoxy must be set by the user, and not remain at %g.\n"),
         cm->hinoxy);
  }

  else if (macro_type == CS_COMBUSTION_EBU) {
    cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                    _(section_name),
                                    "cebu", cm->cebu, 0.);
  }

  else if (macro_type == CS_COMBUSTION_LW) {
    cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                    _(section_name),
                                    "vref", cm->lw.vref, 0.);
    cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                    _(section_name),
                                    "lref", cm->lw.lref, 0.);
    cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                    _(section_name),
                                    "ta", cm->lw.ta, 0.);
    cs_parameters_is_greater_double(CS_ABORT_DELAYED,
                                    _(section_name),
                                    "tstar", cm->lw.tstar, 0.);
  }

  cs_parameters_is_in_range_double(CS_ABORT_IMMEDIATE, _(section_name),
                                   "srrom", cm->srrom, 0., 1.);

  cs_rad_transfer_model_t rt_model_type = cs_glob_rad_transfer_params->type;

  if (cm->isoot > -1 && rt_model_type == CS_RAD_TRANSFER_NONE) {
    cs_parameters_error
      (CS_ABORT_DELAYED,
       _(section_name),
       _("Soot model (isoot = %d) used with no active radiation model.\n"
         "This has no use\n"),
       cm->isoot);
  }

  if (   (cs_glob_combustion_gas_model->type)%2 == 0
      && rt_model_type != CS_RAD_TRANSFER_NONE) {
    cs_parameters_error
      (CS_ABORT_DELAYED,
       _(section_name),
       _("Combustion models in adiabatic conditions are not compatible\n"
         "with a radiative model.\n"
         "  Active combustion model: %s\n"
         "  Active radiation model: %d\n\n"
         "You must use a permeatic conditions variant\n"
         "of the combustion model or deactive the radiation model.\n"),
       _model_name(cm->type), (int)rt_model_type);
  }

  cs_parameters_error_barrier();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the gas combustion module options to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_gas_log_setup(void)
{
  if (cs_glob_combustion_gas_model == NULL)
    return;

  cs_combustion_gas_model_t  *cm = cs_glob_combustion_gas_model;

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "Gas combustion module options\n"
                  "-----------------------------\n\n"));

  cs_log_printf(CS_LOG_SETUP,
                _("  Model type: %s\n\n"),
                _model_name(cm->type));

  const char *janaf_value_str[]
    = {N_("false (user enthalpy/temperature tabulation)"),
       N_("true (JANAF table)")};

  int janaf_msg_id = (cm->use_janaf) ? 1 : 0;
  const char _null_string[] = "<null>";
  const char *data_file = (cm->data_file_name != NULL) ?
    cm->data_file_name : _null_string;

  cs_log_printf(CS_LOG_SETUP,
                _("  Thermochemical data:\n"
                  "    data_file: %s\n"
                  "    use_janaf: %s\n\n"),
                data_file,
                janaf_value_str[janaf_msg_id]);

  switch(cm->isoot) {
  case -1:
    /* No Soot model */
    cs_log_printf(CS_LOG_SETUP,
                  _("    isoot:    -1 (No Soot model)\n\n"));
    break;
  case 0:
    /* constant fraction of product Xsoot */
    cs_log_printf(CS_LOG_SETUP,
                  _("    isoot:     0 (Constant soot yield)\n\n"));
    cs_log_printf(CS_LOG_SETUP,
                  _("  Parameters for the soot model:\n"
                    "    xsoot:  %14.5e (Fraction of product - Used only\n"
                    "            if the soot yield is not defined in the\n"
                    "            thermochemistry data file)\n"
                    "    rosoot: %14.5e (Soot density)\n\n"),
                  cm->xsoot,
                  cm->rosoot);
    break;
  case 1:
    /* 2 equations model of Moss et al. */
    cs_log_printf
      (CS_LOG_SETUP,
       _("    isoot:     1 (2 equations model of Moss et al.)\n\n"));
    cs_log_printf(CS_LOG_SETUP,
                  _("  Parameter for the soot model:\n"
                    "    rosoot: %14.5e (Soot density)\n\n"),
                  cm->rosoot);
    break;
  case 2:
    /* Smoke Point model */
    cs_log_printf(CS_LOG_SETUP, _("    isoot:     2 (Smoke Point model)\n\n"));
    cs_log_printf(CS_LOG_SETUP,
                  _("  Parameter for the soot model:\n"
                    "    rosoot: %14.5e (Soot density)\n"
                    "    lsp_fuel: %14.5e (Laminar smoke point)\n\n"),
                  cm->rosoot,
                  cm->lsp_fuel);
    break;
  default:
    break;
  }

  const char *ipthrm_value_str[] = {N_("0 (no mean pressure computation)"),
                                    N_("1 (mean pressure computation)")};
  cs_log_printf(CS_LOG_SETUP,
                _("    ipthrm:    %s\n\n"),
                _(ipthrm_value_str[cs_glob_fluid_properties->ipthrm]));

  const int *pm_flag = cs_glob_physical_model_flag;

  if (pm_flag[CS_COMBUSTION_3PT] != -1) {
  }
  else if (pm_flag[CS_COMBUSTION_SLFM] != -1) {
    cs_log_printf(CS_LOG_SETUP,
                  _("  Diffusion Flame: Steady laminar flamelet model\n"
                    "    option: %d\n\n"),
                  pm_flag[CS_COMBUSTION_SLFM]);
  }
  else if (pm_flag[CS_COMBUSTION_EBU] != -1) {
    cs_log_printf(CS_LOG_SETUP,
                  _("  Premixed flame: EBU model\n"
                    "    option: %d\n"
                    "    cebu: %14.5e\n\n"),
                  pm_flag[CS_COMBUSTION_EBU],
                  cm->cebu);
  }

  cs_log_printf(CS_LOG_SETUP,
                _("  Time stepping relaxation coefficient\n"
                  "    rho(n+1) = srrom*rho(n) + (1-srrom)*rho(n+1)\n"
                  "    srrom: %14.5e\n"),
                cm->srrom);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add variable fields for gas combustion models.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_gas_add_variable_fields(void)
{
  if (cs_glob_combustion_gas_model == NULL)
    return;

  cs_combustion_gas_model_t  *cm = cs_glob_combustion_gas_model;

  // Key ids for clipping
  const int kscmin = cs_field_key_id("min_scalar_clipping");
  const int kscmax = cs_field_key_id("max_scalar_clipping");

  // Key id for variance
  const int kscavr = cs_field_key_id("first_moment_id");

  /* Define variable fields
     ---------------------- */

  // Enthalpy

  if (cm->type%2 == 1) {
    int f_id = cs_variable_field_create("enthalpy", "Enthalpy",
                                        CS_MESH_LOCATION_CELLS, 1);
    cs_field_t *f = cs_field_by_id(f_id);
    cs_field_pointer_map(CS_ENUMF_(h), f);
    cs_add_model_field_indexes(f->id);

    cs_field_set_key_double(f, kscmin, -cs_math_big_r);
    cs_field_set_key_double(f, kscmax, cs_math_big_r);

    /* set thermal model */
    cs_thermal_model_t *thermal = cs_get_glob_thermal_model();
    thermal->thermal_variable = CS_THERMAL_MODEL_ENTHALPY;
  }

  switch (cm->type / 100) {

  case CS_COMBUSTION_3PT:
    /* Diffusion flame with 3-point chemistry
       -------------------------------------- */
    {
      // Mixture fraction and its variance
      cm->fm = _add_model_variable("mixture_fraction", "Fra_MEL");
      cs_field_set_key_double(cm->fm, kscmin, 0.);
      cs_field_set_key_double(cm->fm, kscmax, 1.);

      cm->fp2m = _add_model_variable("mixture_fraction_variance", "Var_FrMe");
      cs_field_set_key_int(cm->fp2m, kscavr, cm->fm->id);
    }
    break;

  case CS_COMBUSTION_SLFM:
    /* Diffusion flame with steady laminar flamelet approach
       ----------------------------------------------------- */
    {
      const int kivisl  = cs_field_key_id("diffusivity_id");
      const int key_coupled_with_vel_p = cs_field_key_id("coupled_with_vel_p");
      const int key_turb_diff = cs_field_key_id("turbulent_diffusivity_id");
      const int key_sgs_sca_coef = cs_field_key_id("sgs_scalar_flux_coef_id");

      // Mixture fraction and its variance
      cm->fm = _add_model_variable("mixture_fraction", "Fra_MEL");
      cs_field_set_key_double(cm->fm, kscmin, 0.);
      cs_field_set_key_double(cm->fm, kscmax, 1.);

      cs_field_set_key_int(cm->fm, key_coupled_with_vel_p, 1);
      cs_field_set_key_int(cm->fm, kivisl, 0);

      if (cs_glob_turb_model->model == CS_TURB_LES_SMAGO_DYN) {
        cs_field_set_key_int(cm->fm, key_turb_diff, 0);
        cs_field_set_key_int(cm->fm, key_sgs_sca_coef, 0);
      }

      if (cm->mode_fp2m == 0) {
        cm->fp2m = _add_model_variable("mixture_fraction_variance", "Var_FrMe");
        cs_field_set_key_int(cm->fp2m, kscavr, cm->fm->id);
        cs_field_set_key_int(cm->fp2m, key_coupled_with_vel_p, 1);
      }
      else if (cm->mode_fp2m == 1) {
        cm->fsqm = _add_model_variable("mixture_fraction_2nd_moment",
                                       "2nd_Moment_FrMe");
        cs_field_set_key_double(cm->fsqm, kscmin, 0.);
        cs_field_set_key_double(cm->fsqm, kscmax, 1.);
        cs_field_set_key_int(cm->fsqm, key_coupled_with_vel_p, 1);
        cs_field_set_key_int(cm->fsqm, kivisl, 0);

        if (cs_glob_turb_model->model == CS_TURB_LES_SMAGO_DYN) {
          cs_field_set_key_int(cm->fsqm, key_turb_diff, 0);
          cs_field_set_key_int(cm->fsqm, key_sgs_sca_coef, 0);
        }
      }

      if (CS_F_(h) != nullptr) {
        cs_field_set_key_int(CS_F_(h), key_coupled_with_vel_p, 1);
        cs_field_set_key_int(CS_F_(h), kivisl, 0);

        if (cs_glob_turb_model->model == CS_TURB_LES_SMAGO_DYN) {
          cs_field_set_key_int(CS_F_(h), key_turb_diff, 0);
          cs_field_set_key_int(CS_F_(h), key_sgs_sca_coef, 0);
        }
      }

      // Flamelet/Progress variable model
      if (cm->type%100 >= 2) {
        cm->pvm = _add_model_variable("progress_variable", "Prog_Var");
        cs_field_set_key_double(cm->pvm, kscmin, 0.);
        cs_field_set_key_double(cm->pvm, kscmax, cs_math_big_r);
        cs_field_set_key_int(cm->pvm, key_coupled_with_vel_p, 1);
        cs_field_set_key_int(cm->pvm, kivisl, 0);

        if (cs_glob_turb_model->model == CS_TURB_LES_SMAGO_DYN) {
          cs_field_set_key_int(cm->pvm, key_turb_diff, 0);
          cs_field_set_key_int(cm->pvm, key_sgs_sca_coef, 0);
        }
      }
    }
    break;

  case CS_COMBUSTION_EBU:
    /* Premixed flame: Eddy Break-Up model
       ----------------------------------- */
    {
      // Mass fraction of fresh gas
      cm->ygfm = _add_model_variable("fresh_gas_fraction", "Fra_GF");
      cs_field_set_key_double(cm->ygfm, kscmin, 0.);
      cs_field_set_key_double(cm->ygfm, kscmax, 1.);

      // Mixture fraction
      if (   cm->type == CS_COMBUSTION_EBU_VARIABLE_ADIABATIC
          || cm->type == CS_COMBUSTION_EBU_VARIABLE_PERMEATIC) {
        cm->fm = _add_model_variable("mixture_fraction", "Fra_MEL");
        cs_field_set_key_double(cm->fm, kscmin, 0.);
        cs_field_set_key_double(cm->fm, kscmax, 1.);
      }
    }
    break;

  case CS_COMBUSTION_LW:
    /* Premixed flame: Libby-Williams model
       ------------------------------------ */
    {
      cm->fm = _add_model_variable("mixture_fraction", "Fra_MEL");
      cs_field_set_key_double(cm->fm, kscmin, 0.);
      cs_field_set_key_double(cm->fm, kscmax, 1.);

      cm->fp2m = _add_model_variable("mixture_fraction_variance", "Var_FrMe");
      cs_field_set_key_int(cm->fp2m, kscavr, cm->fm->id);

      cm->yfm = _add_model_variable("mass_fraction", "Fra_Mas");
      cs_field_set_key_double(cm->yfm, kscmin, 0.);
      cs_field_set_key_double(cm->yfm, kscmax, 1.);

      cm->yfp2m = _add_model_variable("mass_fraction_variance", "Var_FMa");
      cs_field_set_key_int(cm->yfp2m, kscavr, cm->yfm->id);

      if (cm->type%100 >= 2) {
        cs_field_t *f = _add_model_variable("mass_fraction_covariance",
                                            "COYF_PP4");
        cs_field_set_key_double(f, kscmin, -0.25);
        cs_field_set_key_double(f, kscmax, 0.25);
      }
    }
    break;

  default:
    cs_assert(0);
    break;
  }

  // Map field pointers depending on model fields
  // (mapping nullptr not needed but OK (no-op) where no field is present)

  cs_field_pointer_map(CS_ENUMF_(fm), cm->fm);
  cs_field_pointer_map(CS_ENUMF_(fp2m), cm->fp2m);
  cs_field_pointer_map(CS_ENUMF_(ygfm), cm->ygfm);

  // Soot mass fraction and precursor number
  if (cm->isoot >= 1) {
    cm->fsm = _add_model_variable("soot_mass_fraction", "Fra_Soot");
    cs_field_set_key_double(cm->fsm, kscmin, 0.);
    cs_field_set_key_double(cm->fsm, kscmax, 1.);
    cs_field_pointer_map(CS_ENUMF_(fsm), cm->fsm);

    cm->npm = _add_model_variable("soot_precursor_number", "NPr_Soot");
    cs_field_set_key_double(cm->npm, kscmin, 0.);
    cs_field_set_key_double(cm->npm, kscmax, 1.);
    cs_field_pointer_map(CS_ENUMF_(npm), cm->npm);
  }

  /* Default numerical options
     ------------------------- */

  const int n_fields = cs_field_n_fields();
  const int keysca = cs_field_key_id("scalar_id");

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (!(f->type & CS_FIELD_VARIABLE))
      continue;
    if (f->type & CS_FIELD_CDO || f->type & CS_FIELD_USER)
      continue;
    if (cs_field_get_key_int(f, keysca) <= 0)
      continue;

    cs_equation_param_t *eqp = cs_field_get_equation_param(f);

    // Second order convective scheme
    eqp->blencv = 1.0;

    // Centered convective scheme
    eqp->ischcv = 1;

    // Automatic slope test
    eqp->isstpc = 0;

    // Reconstruct convection and diffusion fluxes at faces
    eqp->ircflu = 1;
  }

  if (cm->type%2 == 1) {
    if (   cm->type == CS_COMBUSTION_3PT_PERMEATIC
        || cm->type/100 == CS_COMBUSTION_EBU
        || cm->type/100 == CS_COMBUSTION_LW) {

      // Although we are in enthalpy formulation, we keep Cp constant.

      cs_fluid_properties_t *fluid_props = cs_get_glob_fluid_properties();
      fluid_props->icp = -1;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add property fields for gas combustion models.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_gas_add_property_fields(void)
{
  if (cs_glob_combustion_gas_model == NULL)
    return;

  cs_combustion_gas_model_t  *cm = cs_glob_combustion_gas_model;

  _add_property_1d("temperature", "Temperature");

  if (   cm->type == CS_COMBUSTION_BSH_ADIABATIC
      || cm->type == CS_COMBUSTION_BSH_PERMEATIC
      || cm->type/100 == CS_COMBUSTION_SLFM)
    cm->t2m = _add_property_1d("temperature_2", "Temperature_2");

  if (cm->type/100 != CS_COMBUSTION_SLFM) {
    cm->ym[0] = _add_property_1d("ym_fuel", "Ym_Fuel");
    cm->ym[1] = _add_property_1d("ym_oxyd", "Ym_Oxyd");
    cm->ym[2] = _add_property_1d("ym_prod", "Ym_Prod");
  }

  switch (cm->type / 100) {

  case CS_COMBUSTION_SLFM:
    /* Diffusion flame with steady laminar flamelet approach
       ----------------------------------------------------- */
    {
      /* In case of the classical steady laminar flamelet model,
         A progress variable is defined as propriety, measuring approximately
         the products (the progress as well)
         Otherwise, it should be transported. */

      if (cm->type%2 == 1)
        cm->xr = _add_property_1d("heat_loss", "Heat Loss");

      cm->hrr = _add_property_1d("heat_release_rate", "Heat Release Rate");

      if (cm->type%100 >= 2)
        cm->omgc = _add_property_1d("omega_c", "Omega C");

      if (cm->type%100 < 2)
        cm->totki = _add_property_1d("total_dissipation", "Total Dissip. Rate");

      if (cm->mode_fp2m == 1)
        cm->recvr = _add_property_1d("reconstructed_fp2m", "rec_fp2m");

      if (cm->n_gas_fl > 0) {
        int n = cs::min(cm->n_gas_fl, CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES - 1);
        for (int i = 0; i < n; i++) {
          char name[64], label[64];
          snprintf(name, 64, "fraction_%s", cm->flamelet_species_name[i]);
          name[63] = '\0';
          snprintf(label, 64, "Fraction_%s", cm->flamelet_species_name[i]);
          label[63] = '\0';
          cm->ym[i] = _add_property_1d(name, label);
        }
      }

      cs_field_t *f = _add_property_1d("ym_progress", "Ym_Progress");
      cm->ym[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES - 1] = f;

      // Add restart_file key info
      cs_field_set_key_int(f,
                           cs_field_key_id("restart_file"),
                           CS_RESTART_AUXILIARY);
    }
    break;

  case CS_COMBUSTION_LW:
    /* Premixed flame: Libby-Williams model
       ------------------------------------ */
    {
      cm->lw.mam = _add_property_1d("molar_mass", "Molar_Mass");
      cm->lw.tsc = _add_property_1d("source_term", "Source_Term");

      for (int i = 0; i < cm->lw.n_dirac; i++) {
        char name[64], label[64];

        snprintf(name, 64, "rho_local_%d", i+1); name[63] = '\0';
        snprintf(label, 64, "Rho_Local_%d", i+1); label[63] = '\0';
        cm->lw.rhol[i] = _add_property_1d(name, label);

        snprintf(name, 64, "temperature_local_%d", i+1); name[63] = '\0';
        snprintf(label, 64, "Temperature_Local_%d", i+1); label[63] = '\0';
        cm->lw.teml[i] = _add_property_1d(name, label);

        snprintf(name, 64, "ym_local_%d", i+1); name[63] = '\0';
        snprintf(label, 64, "Ym_Local_%d", i+1); label[63] = '\0';
        cm->lw.fmel[i] = _add_property_1d(name, label);

        snprintf(name, 64, "w_local_%d", i+1); name[63] = '\0';
        snprintf(label, 64, "w_Local_%d", i+1); label[63] = '\0';
        cm->lw.fmal[i] = _add_property_1d(name, label);

        snprintf(name, 64, "amplitude_local_%d", i+1); name[63] = '\0';
        snprintf(label, 64, "Amplitude_Local_%d", i+1); label[63] = '\0';
        cm->lw.ampl[i] = _add_property_1d(name, label);

        snprintf(name, 64, "chemical_st_local_%d", i+1); name[63] = '\0';
        snprintf(label, 64, "Chemical_ST_Local_%d", i+1); label[63] = '\0';
        cm->lw.tscl[i] = _add_property_1d(name, label);

        snprintf(name, 64, "molar_mass_local_%d", i+1); name[63] = '\0';
        snprintf(label, 64, "M_Local_%d", i+1); label[63] = '\0';
        cm->lw.maml[i] = _add_property_1d(name, label);
      }
    }
    break;

  default:
    break;
  }

  /* Boundary mass fractions
     ----------------------- */

  if (cm->type/100 != CS_COMBUSTION_SLFM) {
    int field_type = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY;
    cm->bym[0] = cs_field_create("boundary_ym_fuel", field_type,
                                 CS_MESH_LOCATION_BOUNDARY_FACES, 1, false);
    cm->bym[1] = cs_field_create("boundary_ym_oxydizer", field_type,
                                 CS_MESH_LOCATION_BOUNDARY_FACES, 1, false);
    cm->bym[2] = cs_field_create("boundary_ym_product", field_type,
                                 CS_MESH_LOCATION_BOUNDARY_FACES, 1, false);
  }

  /* Additional fields for radiation
     ------------------------------- */

  if (   cs_glob_rad_transfer_params->type != CS_RAD_TRANSFER_NONE
      && CS_F_(h) != nullptr) {
    cm->ckabs = _add_property_1d("kabs", "KABS");
    cm->t4m = _add_property_1d("temperature_4", "Temp4");
    cm->t3m = _add_property_1d("temperature_3", "Temp3");
  }
}

/*----------------------------------------------------------------------------*/
/*! \brief Compute molar and mass fractions of
 *         elementary species Ye, Xe (fuel, O2, CO2, H2O, N2) from
 *         global species Yg (fuel, oxidant, products)
 *
 * \param[in]     yg            global mass fraction
 * \param[out]    ye            elementary mass fraction
 * \param[out]    xe            elementary molar fraction
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_gas_yg2xye(const cs_real_t  yg[],
                         cs_real_t        ye[],
                         cs_real_t        xe[])
{
  const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;
  const int n_gas_e = cm->n_gas_el_comp;
  const int n_gas_g = cm->n_gas_species;

  /* Yg -> Ye conversion */

  for (int i = 0; i < n_gas_e; i++) {
    ye[i] = 0;
    for (int j = 0; j < n_gas_g; j++)
      ye[i] += cm->coefeg[j][i] * yg[j];
  }

  /* Verification */

  cs_real_t ytot = 0;
  for (int i = 0; i < n_gas_e; i++)
    ytot += ye[i];

  if (ytot < 0 || (1-ytot) < -cs_math_epzero)
    bft_printf(_(" Warning:\n"
                 " --------\n"
                 "   %s; mass fraction sum outside [0, 1] bounds\n"
                 "   sum_1=1,%d Yi = %g\n\n"),
               __func__, n_gas_e, ytot);

  /* Molar mass mixture */

  cs_real_t mm = 0;
  for (int i = 0; i < n_gas_e; i++) {
    mm += ye[i] / cm->wmole[i];
  }
  mm = 1 / mm;

  /* Ye -> Xe conversion */
  for (int i = 0; i < n_gas_e; i++) {
    xe[i] = ye[i] * mm / cm->wmole[i];
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
