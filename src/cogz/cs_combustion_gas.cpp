/*============================================================================
 * Gas combustion model.
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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------*/

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

#include "base/cs_base.h"
#include "pprt/cs_combustion_model.h"
#include "base/cs_physical_constants.h"
#include "pprt/cs_physical_model.h"
#include "base/cs_log.h"
#include "base/cs_math.h"

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
  \file cs_combustion_gas.c
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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Global variables
 *============================================================================*/

/*! Combustion model parameters structure */

cs_combustion_gas_model_t  *cs_glob_combustion_gas_model = NULL;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_co_models_init(void);

void
cs_f_ppincl_combustion_init(void);

void
cs_f_ppcpfu_models_init(void);

void
cs_f_thch_models_init(void);

void
cs_f_combustion_gas_get_data_file_name(int           name_max,
                                       const char  **name,
                                       int          *name_len);

void
cs_f_init_steady_laminar_flamelet_library(double  **flamelet_library_p);

/*----------------------------------------------------------------------------
 * Initialize steady flamelet library
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * TODO: when Fortran bindings are removed, this function can be simplified,
 * so as to keep only C/C++ part (instead of being removed).
 *----------------------------------------------------------------------------*/

void
cs_f_init_steady_laminar_flamelet_library(double  **flamelet_library_p)
{
  *flamelet_library_p = nullptr;

  if (cs_glob_combustion_gas_model != NULL) {

    cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

    size_t n = cm->nlibvar * cm->nxr * cm->nki * cm->nzvar * cm->nzm;

    CS_FREE(cm->flamelet_library_p);
    CS_MALLOC(cm->flamelet_library_p, n, cs_real_t);

    cs_real_t *_flamelet_library = cm->flamelet_library_p;
    for (size_t i = 0; i < n; i++)
      _flamelet_library[i] = 0;

    *flamelet_library_p = cm->flamelet_library_p;
  }
}

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
    CS_FREE(cm->flamelet_library_p);
    CS_FREE(cm);

    cs_glob_combustion_gas_model = cm;

  }
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return the name of the data file
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   name_max <-- maximum name length
 *   name     --> pointer to associated length
 *   name_len --> length of associated length
 *----------------------------------------------------------------------------*/

void
cs_f_combustion_gas_get_data_file_name(int           name_max,
                                       const char  **name,
                                       int          *name_len)
{
  assert(cs_glob_combustion_gas_model != NULL);

  *name = cs_glob_combustion_gas_model->data_file_name;
  *name_len = strlen(*name);
  if (*name_len > name_max) {
    bft_error
      (__FILE__, __LINE__, 0,
       _("Error retrieving thermochemistry data file  (\"%s\"):\n"
         "Fortran caller name length (%d) is too small for current length %d."),
       *name, name_max, *name_len);
  }
}

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
    BFT_FREE(cs_glob_combustion_gas_model);
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

  BFT_MALLOC(cm, 1, cs_combustion_gas_model_t);

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

  cm->data_file_name = NULL;

  cm->use_janaf = true;

  cm->iic = 0;
  cm->iico2 = 0;
  cm->iio2 = 0;
  cm->isoot = -1;

  cm->hinfue = cs_math_big_r;
  cm->tinfue = 0;
  cm->tinoxy = 0;
  cm->hinoxy = cs_math_big_r;

  cm->xsoot = 0.;
  cm->rosoot = 0.;
  cm->lsp_fuel = 0.;

  cm->frmel = 0.;
  cm->tgf = 300.;

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
      cm->cpgazg[i][j] = 0;
    }
  }

  cm->srrom = 0.95;

  /* Libby Williams model */
  cm->fmin_lwc = 0.;
  cm->fmax_lwc = 1.;
  cm->hmin_lwc = 0.;
  cm->hmax_lwc = 0.;

  /*! Steady flamelet model */

  cm->ngazfl = -1;
  cm->nki = -1;
  cm->nxr = -1;
  cm->nzm = -1;
  cm->nzvar = -1;
  cm->nlibvar = -1;
  cm->ikimid = 1;
  cm->mode_fp2m = 1;

  cm->flamelet_library_p = nullptr;

  /* Set finalization callback */

  cs_base_at_finalize(_combustion_gas_finalize);

  /* Set mappings with Fortran */

  cs_f_ppincl_combustion_init();
  cs_f_co_models_init();
  cs_f_ppcpfu_models_init();
  cs_f_thch_models_init();

  return cm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the thermochemical data file name.
 *
 * \param[in] file_name  name of the file.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_gas_set_thermochemical_data_file(const char  *file_name)
{
  cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  if (cm == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: gas combustion model not active."), __func__);

  BFT_FREE(cm->data_file_name);
  if (file_name != NULL) {
    size_t l = strlen(file_name)+1;
    BFT_REALLOC(cm->data_file_name, l, char);
    strncpy(cm->data_file_name, file_name, l);
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
                  "Combustion module options\n"
                    "-------------------------\n\n"));

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
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
