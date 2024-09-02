/*============================================================================
 * Management of the GUI parameters file: specific physics
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_selector.h"

#include "cs_atmo.h"
#include "cs_base.h"
#include "cs_cf_model.h"
#include "cs_coal.h"
#include "cs_combustion_gas.h"
#include "cs_gui_util.h"
#include "cs_gui.h"
#include "cs_mesh.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_parameters.h"
#include "cs_physical_model.h"
#include "cs_selector.h"
#include "cs_physical_properties.h"
#include "cs_physical_constants.h"
#include "cs_rad_transfer.h"
#include "cs_elec_model.h"
#include "cs_vof.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gui_specific_physics.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* debugging switch */
#define _XML_DEBUG_ 0

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Return an integer for class number.
 *
 * parameters:
 *   tn   <-- tree node associated with solid fuel
 *   type -->   type of diameter definition
 *
 * returns:
 *    number of class
 *----------------------------------------------------------------------------*/

static int
_get_nb_class(cs_tree_node_t  *tn,
              int              type)
{
  int value = 0;

  cs_tree_node_t  *tn_c = cs_tree_node_get_child(tn, "class");
  if (type == 1)
    value = cs_tree_get_node_count(tn_c, "diameter");
  else if (type == 2)
    value = cs_tree_get_node_count(tn_c, "mass_percent");

  return value;
}

/*-----------------------------------------------------------------------------
 * Get the kinetic model (CO2 or CO transport).
 *
 * parameters:
 *   tn_sf <-- tree node associated with solid fuels
 *   model <-> type of kinetic model to use
 *----------------------------------------------------------------------------*/

static void
_get_kinetic_model(cs_tree_node_t  *tn_sf,
                   int             *model)
{
  const char *s = cs_tree_node_get_child_value_str(tn_sf, "kinetic_model");

  if (s != NULL) {
    if (!strcmp(s, "unused"))
      *model = 0;
    else if (!strcmp(s, "co2_ym_transport"))
      *model = 1;
    else if (!strcmp(s, "co_ym_transport"))
      *model = 2;
    else
      bft_error(__FILE__, __LINE__, 0, _("Invalid value for node %s: %s"),
                "kinetic_model", s);
  }
}

/*-----------------------------------------------------------------------------
 * Return integer for diameter type.
 *
 * parameters:
 *   tn <-- parent node ("solid_fuel")
 *
 * returns:
 *    integer for diameter type
 *----------------------------------------------------------------------------*/

static int
_get_diameter_type(cs_tree_node_t  *tn)
{
  int ichoice = 0;

  const char *type_s = cs_tree_node_get_child_value_str(tn, "diameter_type");

  if (type_s == NULL) {
    ichoice = 1;
  }
  else {
    if (cs_gui_strcmp(type_s, "automatic"))
      ichoice = 1;
    else if (cs_gui_strcmp(type_s, "rosin-rammler_law"))
      ichoice = 2;
    else
      bft_error(__FILE__, __LINE__, 0, _("Invalid value for node %s: %s"),
                "diameter_type", type_s);
  }

  return ichoice;
}

/*-----------------------------------------------------------------------------
 * Return real value associated with coal species child node
 *
 * parameters:
 *    tn   <-- associated "solid_fuel" node
 *    type <-- type of element
 *
 * returns:
 *    composition on dry value
 *----------------------------------------------------------------------------*/

static double
_get_solid_fuel_child_real(cs_tree_node_t  *tn,
                           const char      *type)
{
  const cs_real_t *v_r = cs_tree_node_get_child_values_real(tn, type);

  if (v_r == NULL) {
    const char *id = cs_tree_node_get_tag(tn, "fuel_id");
    bft_error(__FILE__, __LINE__, 0,
              _("Missing %s/%s node or value for fuel %s"),
              tn->name, type, id);
  }

  return v_r[0];
}

/*-----------------------------------------------------------------------------
 * Return integer for PCI type and get associated value
 *
 * parameters:
 *    tn_coal   <-- "solid_fuel" node associated with coal
 *    pci_type  --> type
 *    pci_value <-> associated value if type < 6
 *----------------------------------------------------------------------------*/

static void
_get_pci_type_and_value(cs_tree_node_t  *tn_coal,
                        int             *pci_type,
                        cs_real_t       *pci_value)
{
  int ichoice = 0;

  const char path0[] = "Heating_model";
  cs_tree_node_t *tn = cs_tree_node_get_child(tn_coal, path0);

  if (tn == NULL)
    bft_error(__FILE__, __LINE__, 0, _("Missing %s child for node %s."),
              path0, tn_coal->name);

  const char *choice = cs_tree_node_get_tag(tn, "choice");

  if (choice != NULL) {

    if (cs_gui_strcmp(choice, "LHV")) {
      const char *type_s = cs_tree_node_get_child_value_str(tn, "type");
      if (type_s == NULL)
        ichoice = 1;
      else {
        if (cs_gui_strcmp(type_s, "dry_basis"))
          ichoice = 1;
        else if (cs_gui_strcmp(type_s, "dry_ash_free"))
          ichoice = 0;
        else if (cs_gui_strcmp(type_s, "as_received"))
          ichoice = 2;
        else
          bft_error(__FILE__, __LINE__, 0,
                    _("Invalid value for node %s/%s: %s"),
                    path0, "type", type_s);
      }
    }
    else if (cs_gui_strcmp(choice, "HHV")) {
      const char *type_s = cs_tree_node_get_child_value_str(tn, "type");
      if (type_s == NULL)
        ichoice = 4;
      else {
        if (cs_gui_strcmp(type_s, "dry_basis"))
          ichoice = 4;
        else if (cs_gui_strcmp(type_s, "dry_ash_free"))
          ichoice = 3;
        else if (cs_gui_strcmp(type_s, "as_received"))
          ichoice = 5;
        else
          bft_error(__FILE__, __LINE__, 0, _("Invalid value for node %s: %s"),
                    "Heating_model/type", type_s);
      }
    }
    else if (cs_gui_strcmp(choice, "IGT_correlation"))
      ichoice = 6;
    else
      bft_error(__FILE__, __LINE__, 0, _("Invalid value for node %s: %s"),
                "Heating_model/choice", choice);
  }

  *pci_type = ichoice;

  if (ichoice < 6)
    cs_gui_node_get_child_real(tn, "value", pci_value);
}

/*-----------------------------------------------------------------------------
 * Return integer for Y1/Y2 coeffficient type and get associated value
 *
 * parameters:
 *    tn_dv    <-- "devolatilisation_parameters" node associated with coal
 *    y1_type  --> type for y1
 *    y2_type  --> type for y2
 *    y1_value <-> associated y1 value if type > 0
 *    y2_value <-> associated y1 value if type > 0
 *----------------------------------------------------------------------------*/

static void
_get_y1y2_coefficient_values(cs_tree_node_t  *tn_dv,
                             int             *y1_type,
                             int             *y2_type,
                             cs_real_t       *y1_value,
                             cs_real_t       *y2_value)
{
  int ichoice = 0;

  const char path0[] = "stoichiometric_coefficient";
  cs_tree_node_t *tn = cs_tree_get_node(tn_dv, path0);

  if (tn == NULL)
    bft_error(__FILE__, __LINE__, 0, _("Missing %s descendant for node %s."),
              path0, tn_dv->name);

  const char *type_s = cs_tree_node_get_tag(tn, "type");

  if (type_s != NULL) {
    if (cs_gui_strcmp(type_s, "automatic_CHONS"))
      ichoice = 0;
    else if (cs_gui_strcmp(type_s, "user_define"))
      ichoice = 1;
    else if (cs_gui_strcmp(type_s, "automatic_formula"))
      ichoice = 2;
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Invalid value for node %s/%s: %s"),
                path0, "type", type_s);
  }

  *y1_type = ichoice;
  *y2_type = ichoice;

  if (ichoice > 0) {
    cs_gui_node_get_child_real(tn, "Y1", y1_value);
    cs_gui_node_get_child_real(tn, "Y2", y2_value);
  }
}

/*-----------------------------------------------------------------------------
 * Return double for char combustion value.
 *
 * parameters:
 *    tn_cc   <-- associated "char_combustion" node
 *    species <-- species choice
 *    key     <-- value name
 *
 * returns:
 *    queried value
 *----------------------------------------------------------------------------*/

static cs_real_t
_get_cc_specie_value(cs_tree_node_t  *tn_cc,
                     const char      *species,
                     const char      *key)
{
  cs_real_t retval = 0;

  cs_tree_node_t *tn = NULL;
  for (tn = cs_tree_node_get_child(tn_cc, "specie");
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn)) {
    if (cs_gui_strcmp(cs_tree_node_get_tag(tn, "nature"), species))
      break;
  }

  if (tn == NULL)
    bft_error(__FILE__, __LINE__, 0, _("Missing %s specie child for node %s."),
              species, tn_cc->name);

  retval = _get_solid_fuel_child_real(tn, key);

  return retval;
}

/*-----------------------------------------------------------------------------
 * Return integer for char combustion order of reaction.
 *
 * parameters:
 *    tn_cc   <-- associated "char_combustion" node
 *    species <-- species choice
 *    key     <-- value name
 *
 * returns:
 *    queried value
 *----------------------------------------------------------------------------*/

static int
_get_cc_reaction_order(cs_tree_node_t  *tn_cc,
                       const char      *species)
{
  int ichoice = 0;

  cs_tree_node_t *tn = NULL;
  for (tn = cs_tree_node_get_child(tn_cc, "specie");
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn)) {
    if (cs_gui_strcmp(cs_tree_node_get_tag(tn, "nature"), species))
      break;
  }

  if (tn == NULL)
    bft_error(__FILE__, __LINE__, 0, _("Missing %s specie child for node %s."),
              species, tn_cc->name);

  const char *choice = cs_tree_node_get_tag
                         (cs_tree_node_get_child(tn, "order_of_reaction"),
                          "choice");

  if (choice != NULL) {
    if (cs_gui_strcmp(choice, "0.5"))
      ichoice = 0;
    else if (cs_gui_strcmp(choice, "1"))
      ichoice = 1;
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Invalid value for node %s/%s: %s"),
                tn->name, "order_of_reaction/choice", choice);
  }

  return ichoice;
}

/*-----------------------------------------------------------------------------
 * Return for all oxidants (1 to 3) composition for a given species
 *
 * parameters:
 *    tn_oxi      <-- "oxidants" node
 *    species     <-- species key
 *    composition --> oxidants composition
 *----------------------------------------------------------------------------*/

static void
_get_oxidants_composition(cs_tree_node_t  *tn_oxi,
                          const char      *species_key,
                          cs_real_t        composition[3])
{
  int ioxy = 0;

  for (ioxy = 0; ioxy < 3; ioxy++)
    composition[ioxy] = 0;

  ioxy = 0;
  for (cs_tree_node_t *tn = cs_tree_node_get_child(tn_oxi, "oxidant");
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn), ioxy++) {

    int _ioxy = ioxy;
    const int *v_i = cs_tree_node_get_child_values_int(tn, "ox_id");
    if (v_i != NULL)
      _ioxy = v_i[0] - 1;
    if (_ioxy < 0 || _ioxy > 2)
      bft_error(__FILE__, __LINE__, 0,
                _("oxidant node id (%d) out of [1, 3] range."),
                _ioxy + 1);

    cs_gui_node_get_child_real(tn, species_key, &(composition[_ioxy]));
  }
}

/*-----------------------------------------------------------------------------
 * Return integer parameter for oxidant type.
 *
 * parameters:
 *   tn_oxi <-- "oxidants" tree node
 *
 * returns:
 *    parameter for oxidant type
 *----------------------------------------------------------------------------*/

static int
_get_oxidant_type(cs_tree_node_t  *tn_oxi)
{
  int   ichoice = 0;

  const char *o_type = cs_tree_node_get_child_value_str(tn_oxi, "oxidant_type");

  if (o_type != NULL) {
    if (! strcmp(o_type, "volumic_percent"))
      ichoice = 0;
    else if (! strcmp(o_type, "molar"))
      ichoice = 1;
    else
      bft_error(__FILE__, __LINE__, 0, _("Invalid value for node %s: %s"),
                "oxidant_type", o_type);
  }

  return ichoice;
}

/*-----------------------------------------------------------------------------
 * Return valur for absorption coefficient.
 *
 * returns:
 *    absorption coefficient
 *----------------------------------------------------------------------------*/

static cs_real_t
_get_absorption_coefficient(void)
{
  cs_real_t retval = 0;

  const char path_ac[] = "thermophysical_models/radiative_transfer/" \
                         "absorption_coefficient";

  cs_tree_node_t *tn = cs_tree_get_node(cs_glob_tree, path_ac);

  if (tn == NULL)
    bft_error(__FILE__, __LINE__, 0, _("Missing tree node %s."),
              path_ac);

  cs_gui_node_get_real(tn, &retval);

  return retval;
}

/*-----------------------------------------------------------------------------
 * Return integer for NOx feature status
 *
 *   parameters:
 *    tn_nox   <-- "nox_formation" tree node
 *    keyword  --> value of attribute node
 *----------------------------------------------------------------------------*/

static void
_get_nox_reburning(cs_tree_node_t  *tn_nox,
                   int             *keyword)
{
  const char *choice = cs_tree_node_get_child_value_str(tn_nox,
                                                        "reburning_model");

  if (cs_gui_strcmp(choice, "unused"))
    *keyword = 0;
  else if (cs_gui_strcmp(choice, "chen"))
    *keyword = 1;
  else
    *keyword = 2;
}

/*----------------------------------------------------------------------------
 * Atmospheric flows: read of meteorological file of data
 *----------------------------------------------------------------------------*/

static void
_gui_atmo_get_set_meteo_profile(void)
{
  const char path_af[] = "thermophysical_models/atmospheric_flows";

  cs_tree_node_t *tn = cs_tree_get_node(cs_glob_tree, path_af);

  if (tn == NULL)
    return;

  int is_meteo_file = 0;
  int is_large_scale_meteo = 0;

  cs_gui_node_get_child_status_int(tn, "read_meteo_data",
                                   &is_meteo_file);
  cs_gui_node_get_child_status_int(tn, "large_scale_meteo",
                                   &is_large_scale_meteo);

  if (is_meteo_file && !is_large_scale_meteo) {
    cs_glob_atmo_option->meteo_profile = 1;
    const char *cstr = cs_tree_node_get_child_value_str(tn, "meteo_data");

    /* Copy string */
    if (cstr != NULL)
      cs_atmo_set_meteo_file_name(cstr);

  }
  else if (is_large_scale_meteo && !is_meteo_file) {
    cs_glob_atmo_option->meteo_profile = 2;
    const char *str_latitude
      = cs_tree_node_get_child_value_str(tn, "latitude");
    const char *str_longitude
      = cs_tree_node_get_child_value_str(tn, "longitude");
    const char *str_domain_orient
      = cs_tree_node_get_child_value_str(tn, "domain_orientation");
    const char *str_wind_dir
      = cs_tree_node_get_child_value_str(tn, "wind_direction");

    const char *str_meteo_z0
      = cs_tree_node_get_child_value_str(tn, "meteo_z0");
    const char *str_meteo_uref
      = cs_tree_node_get_child_value_str(tn, "meteo_uref");
    const char *str_meteo_ustar
      = cs_tree_node_get_child_value_str(tn, "meteo_ustar");
    const char *str_meteo_dlmo
      = cs_tree_node_get_child_value_str(tn, "meteo_dlmo");
    const char *str_meteo_zref
      = cs_tree_node_get_child_value_str(tn, "meteo_zref");
    const char *str_meteo_psea
      = cs_tree_node_get_child_value_str(tn, "meteo_psea");
    const char *str_meteo_t0
      = cs_tree_node_get_child_value_str(tn, "meteo_t0");
    const char *str_meteo_qw0
      = cs_tree_node_get_child_value_str(tn, "meteo_qw0");
    const char *str_meteo_qwstar
      = cs_tree_node_get_child_value_str(tn, "meteo_qwstar");

    const char *str_syear = cs_tree_node_get_child_value_str(tn, "start_year");
    const char *str_sday = cs_tree_node_get_child_value_str(tn, "start_day");
    const char *str_shour = cs_tree_node_get_child_value_str(tn, "start_hour");
    const char *str_smin = cs_tree_node_get_child_value_str(tn, "start_min");
    const char *str_ssec = cs_tree_node_get_child_value_str(tn, "start_sec");

    if (str_latitude != NULL)
      cs_glob_atmo_option->latitude = atof(str_latitude);
    if (str_longitude != NULL)
      cs_glob_atmo_option->longitude = atof(str_longitude);
    if (str_domain_orient != NULL)
      cs_glob_atmo_option->domain_orientation = atof(str_domain_orient);
    if (str_wind_dir != NULL)
      cs_glob_atmo_option->meteo_angle = atof(str_wind_dir);

    if (str_meteo_z0 != NULL)
      cs_glob_atmo_option->meteo_z0  = atof(str_meteo_z0);
    if (str_meteo_uref != NULL)
      cs_glob_atmo_option->meteo_uref  = atof(str_meteo_uref);
    if (str_meteo_ustar != NULL)
      cs_glob_atmo_option->meteo_ustar0  = atof(str_meteo_ustar);
    if (str_meteo_dlmo != NULL)
      cs_glob_atmo_option->meteo_dlmo  = atof(str_meteo_dlmo);
    if (str_meteo_zref != NULL)
      cs_glob_atmo_option->meteo_zref  = atof(str_meteo_zref);
    if (str_meteo_psea != NULL)
      cs_glob_atmo_option->meteo_psea  = atof(str_meteo_psea);
    if (str_meteo_t0 != NULL)
      cs_glob_atmo_option->meteo_t0  = atof(str_meteo_t0);
    if (str_meteo_qw0 != NULL)
      cs_glob_atmo_option->meteo_qw0  = atof(str_meteo_qw0);
    if (str_meteo_qwstar != NULL)
      cs_glob_atmo_option->meteo_qwstar  = atof(str_meteo_qwstar);

    if(str_syear != NULL)
      cs_glob_atmo_option->syear = atoi(str_syear);
    if(str_sday != NULL)
      cs_glob_atmo_option->squant = atoi(str_sday);
    if(str_shour != NULL)
      cs_glob_atmo_option->shour = atoi(str_shour);
    if(str_smin != NULL)
      cs_glob_atmo_option->smin = atoi(str_smin);
    if(str_ssec != NULL)
      cs_glob_atmo_option->ssec = atof(str_ssec);
    /*TODO how to split which way profile definition we are using ? */

  }

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--meteo_profile  = %i\n", cs_glob_atmo_option->meteo_profile);
#endif
}

/*-----------------------------------------------------------------------------
 * Return 1 if a specific physics model is activated, 0 otherwise.
 *----------------------------------------------------------------------------*/

static int
_get_active_thermophysical_model(char  **model_name,
                                 char  **model_value)
{
  int isactiv = 0;

  if (*model_name != NULL && *model_value != NULL) {
    isactiv = 1;
    return isactiv;
  }
  else {
    BFT_FREE(*model_name);
    BFT_FREE(*model_value);
  }

  const char *model = NULL, *name = NULL;

  const char *name_m[] = {"solid_fuels",
                          "joule_effect",
                          "atmospheric_flows",
                          "compressible_model",
                          "groundwater_model",
                          "hgn_model"};
  const char *name_o[] = {"gas_combustion"};

  int n_names_m = sizeof(name_m) / sizeof(name_m[0]);
  int n_names_o = sizeof(name_o) / sizeof(name_o[0]);

  cs_tree_node_t *tn0 = cs_tree_get_node(cs_glob_tree, "thermophysical_models");

  /* Loop on model nodes to compare to the listed ones */

  if (tn0 != NULL) {
    for (cs_tree_node_t *tn = tn0->children;
         tn != NULL && name == NULL;
         tn = tn->next) {
      for (int j = 0; j < n_names_m && name == NULL; j++) {
        if (! strcmp(tn->name, name_m[j])) {
          model = cs_tree_node_get_tag(tn, "model");
          if (model != NULL && !cs_gui_strcmp(model, "off"))
            name = name_m[j];
        }
      }
      for (int j = 0; j < n_names_o && name == NULL; j++) {
        if (! strcmp(tn->name, name_o[j])) {
          model = cs_tree_node_get_tag(tn, "option");
          if (model != NULL && !cs_gui_strcmp(model, "off"))
            name = name_o[j];
        }
      }
    }
  }

  if (name != NULL) {
    BFT_MALLOC(*model_name, strlen(name)+1, char);
    strcpy(*model_name, name);

    BFT_MALLOC(*model_value, strlen(model)+1, char);
    strcpy(*model_value, model);

    isactiv = 1;
  }

  return isactiv;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Temperature for D3P Gas Combustion
 *
 * Fortran Interface:
 *
 * Toxy   <--   Oxidant temperature
 * Tfuel  <--   Fuel temperature
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicpi2, UICPI2) (double *const toxy,
                                double *const tfuel)
{
  cs_gui_fluid_properties_value("reference_oxydant_temperature", toxy);
  cs_gui_fluid_properties_value("reference_fuel_temperature", tfuel);

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--toxy  = %f\n", *toxy);
  bft_printf("--tfuel  = %f\n", *tfuel);
#endif
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Activate specific physical models based on XML settings.
 *----------------------------------------------------------------------------*/

void
cs_gui_physical_model_select(void)
{
  int isactiv = 0;

  /* Look for the active specific physics and give the value of the associated
     model attribute */

  char *model_name = NULL, *model_value = NULL;

  isactiv = _get_active_thermophysical_model(&model_name, &model_value);

  if (isactiv)  {

    if (cs_gui_strcmp(model_name, "solid_fuels")) {
      if (cs_gui_strcmp(model_value, "homogeneous_fuel"))
        cs_coal_model_set_model(CS_COMBUSTION_COAL_STANDARD);
      else if (cs_gui_strcmp(model_value,
                             "homogeneous_fuel_moisture"))
        cs_coal_model_set_model(CS_COMBUSTION_COAL_WITH_DRYING);
      else
        bft_error(__FILE__, __LINE__, 0,
                  _("Invalid coal model: %s."), model_value);
    }
    else if (cs_gui_strcmp(model_name, "gas_combustion")) {

      cs_tree_node_t *tn
        = cs_tree_get_node(cs_glob_tree, "thermophysical_models/gas_combustion");

      /* For gas combustion, use "option" instead of "model" tag,
         which has another other role */
      const char *model = cs_tree_node_get_tag(tn, "option");

      if (model != NULL && !cs_gui_strcmp(model, "off")) {

        if (cs_gui_strcmp(model_value, "adiabatic"))
          cs_combustion_gas_set_model(CS_COMBUSTION_3PT_ADIABATIC);
        else if (cs_gui_strcmp(model_value, "extended"))
          cs_combustion_gas_set_model(CS_COMBUSTION_3PT_PERMEATIC);
        else if (cs_gui_strcmp(model_value, "spalding"))
          cs_combustion_gas_set_model(CS_COMBUSTION_EBU_CONSTANT_ADIABATIC);
        else if (cs_gui_strcmp(model_value, "enthalpy_st"))
          cs_combustion_gas_set_model(CS_COMBUSTION_EBU_CONSTANT_PERMEATIC);
        else if (cs_gui_strcmp(model_value, "mixture_st"))
          cs_combustion_gas_set_model(CS_COMBUSTION_EBU_VARIABLE_ADIABATIC);
        else if (cs_gui_strcmp(model_value, "enthalpy_mixture_st"))
          cs_combustion_gas_set_model(CS_COMBUSTION_EBU_VARIABLE_PERMEATIC);
        else if (cs_gui_strcmp(model_value, "2-peak_adiabatic"))
          cs_combustion_gas_set_model(CS_COMBUSTION_LW_2PEAK_ADIABATIC);
        else if (cs_gui_strcmp(model_value, "2-peak_enthalpy"))
          cs_combustion_gas_set_model(CS_COMBUSTION_LW_2PEAK_PERMEATIC);
        else if (cs_gui_strcmp(model_value, "3-peak_adiabatic"))
          cs_combustion_gas_set_model(CS_COMBUSTION_LW_3PEAK_ADIABATIC);
        else if (cs_gui_strcmp(model_value, "3-peak_enthalpy"))
          cs_combustion_gas_set_model(CS_COMBUSTION_LW_3PEAK_PERMEATIC);
        else if (cs_gui_strcmp(model_value, "4-peak_adiabatic"))
          cs_combustion_gas_set_model(CS_COMBUSTION_LW_4PEAK_ADIABATIC);
        else if (cs_gui_strcmp(model_value, "4-peak_enthalpy"))
          cs_combustion_gas_set_model(CS_COMBUSTION_LW_4PEAK_PERMEATIC);
        else
          bft_error(__FILE__, __LINE__, 0,
                    _("Invalid gas combustion flow model: %s."),
                    model_value);

        cs_combustion_gas_model_t  *cm = cs_glob_combustion_gas_model;

        /* Read uniform variable thermodynamical pressure (ipthrm) */
        cs_fluid_properties_t *phys_pp = cs_get_glob_fluid_properties();
        cs_gui_node_get_child_status_int(tn, "thermodynamical_pressure",
                                         &(phys_pp->ipthrm));

        /* Read the soot model (isoot, rosoot, xsoot) */
        cs_tree_node_t *tn_soot  = cs_tree_get_node(tn, "soot_model");

        const char *model_soot = cs_tree_node_get_child_value_str(tn_soot,
                                                                  "model");

        if (model_soot != NULL && !cs_gui_strcmp(model_soot, "off")) {
          if (cs_gui_strcmp(model_soot, "constant_soot_yield")) {
            cm->isoot = 0;
            cs_gui_node_get_child_real(tn_soot, "soot_density",
                                       &(cm->rosoot));
            cs_gui_node_get_child_real(tn_soot, "soot_fraction",
                                       &(cm->xsoot));}
          else if (cs_gui_strcmp(model_soot, "moss")) {
            cm->isoot = 1;
            cs_gui_node_get_child_real(tn_soot, "soot_density",
                                       &(cm->rosoot));}
        }

        /* Read name of thermochemistry data file */

        const char *f_name = cs_tree_node_get_child_value_str(tn, "data_file");
        if (f_name != NULL)
          cs_combustion_gas_set_thermochemical_data_file(f_name);

      }
    }
    else if (cs_gui_strcmp(model_name, "atmospheric_flows")) {
      if (cs_gui_strcmp(model_value, "constant"))
        cs_glob_physical_model_flag[CS_ATMOSPHERIC] = CS_ATMO_CONSTANT_DENSITY;
      else if (cs_gui_strcmp(model_value, "dry"))
        cs_glob_physical_model_flag[CS_ATMOSPHERIC] = CS_ATMO_DRY;
      else if (cs_gui_strcmp(model_value, "humid"))
        cs_glob_physical_model_flag[CS_ATMOSPHERIC] = CS_ATMO_HUMID;

      /* Combination of the two modules */
      else if (cs_gui_strcmp(model_value, "humid_ctwr")) {
        cs_glob_physical_model_flag[CS_COOLING_TOWERS] = 1;
        cs_glob_physical_model_flag[CS_ATMOSPHERIC] = CS_ATMO_HUMID;
      }
      /* Cooling towers only */
      else if (cs_gui_strcmp(model_value, "ctwr"))
        cs_glob_physical_model_flag[CS_COOLING_TOWERS] = 1;
      else
        bft_error(__FILE__, __LINE__, 0,
                  _("Invalid atmospheric flow model: %s."),
                  model_value);

      /* Get and set meteo file if given */
      _gui_atmo_get_set_meteo_profile();

    }
    else if (cs_gui_strcmp(model_name, "joule_effect")) {
      if (cs_gui_strcmp(model_value, "joule")) {

        cs_tree_node_t *tn
          = cs_tree_get_node(cs_glob_tree,
                             "thermophysical_models/joule_effect/joule_model");
        const char *model = cs_tree_node_get_tag(tn, "model");

        if (cs_gui_strcmp(model, "AC/DC"))
          cs_glob_physical_model_flag[CS_JOULE_EFFECT] = 1;
        else if (cs_gui_strcmp(model, "three-phase"))
          cs_glob_physical_model_flag[CS_JOULE_EFFECT] = 2;
        else if (cs_gui_strcmp(model, "AC/DC+Transformer"))
          cs_glob_physical_model_flag[CS_JOULE_EFFECT] = 3;
        else if (cs_gui_strcmp(model, "three-phase+Transformer"))
          cs_glob_physical_model_flag[CS_JOULE_EFFECT] = 4;
        else
          bft_error(__FILE__, __LINE__, 0,
                    _("Invalid joule model: %s."),
                    model_value);

      }
      else if (cs_gui_strcmp(model_value, "arc"))
        cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] = 2;
      else
        bft_error(__FILE__, __LINE__, 0,
                  _("Invalid electrical model: %s."),
                  model_value);
    }
    else if (cs_gui_strcmp(model_name, "compressible_model")) {
      if (cs_gui_strcmp(model_value, "constant_gamma")) {
        cs_glob_physical_model_flag[CS_COMPRESSIBLE] = 0;
        cs_cf_model_t *cf_mdl = cs_get_glob_cf_model();
        cf_mdl->ieos = 1;
      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  _("Invalid compressible model: %s."),
                  model_value);
    }
    else if (cs_gui_strcmp(model_name, "groundwater_model")) {
      if (cs_gui_strcmp(model_value, "groundwater"))
        cs_glob_physical_model_flag[CS_GROUNDWATER] = 1;
    }
    else if (cs_gui_strcmp(model_name, "hgn_model")) {
      cs_vof_parameters_t *vof_param = cs_get_glob_vof_parameters();
      if (cs_gui_strcmp(model_value, "merkle_model")) {
        vof_param->vof_model = CS_VOF_ENABLED | CS_VOF_MERKLE_MASS_TRANSFER;
      } else {
        vof_param->vof_model = CS_VOF_ENABLED + CS_VOF_FREE_SURFACE;
      }
    }
  }

#if _XML_DEBUG_
  bft_printf("%s\n", __func__);
  if (isactiv) {
    bft_printf("--thermophysical model: %s\n", model_name);
    bft_printf("--thermophysical value: %s\n", model_value);
  }
#endif

  BFT_FREE(model_name);
  BFT_FREE(model_value);
}

/*----------------------------------------------------------------------------
 * Indirection between the solver numbering and the XML one
 * for physical properties of the activated specific physics
 * (pulverized solid fuels)
 *----------------------------------------------------------------------------*/

void
cs_gui_coal_model(void)
{
  cs_coal_model_t *cm = cs_glob_coal_model;
  assert(cm != NULL);

  /* Read gas mix absorption coefficient */
  if (cs_glob_rad_transfer_params->type > CS_RAD_TRANSFER_NONE)
    cm->ckabs0 = _get_absorption_coefficient();

  /* Solid fuel model node */

  const char path_sf[] = "thermophysical_models/solid_fuels";
  cs_tree_node_t *tn_sf = cs_tree_get_node(cs_glob_tree, path_sf);

  if (tn_sf == NULL)
    bft_error(__FILE__, __LINE__, 0, _("Missing tree node %s."),
              path_sf);

  /* Heterogeneous combustion options (shrinking sphere model) */
  cs_gui_node_get_child_status_int(tn_sf, "CO2_kinetics", &(cm->ihtco2));
  cs_gui_node_get_child_status_int(tn_sf, "H2O_kinetics", &(cm->ihth2o));

  /* Kinetic model (CO2 or CO transport) */
  _get_kinetic_model(tn_sf, &(cm->ieqco2));

  /* Solid fuel definitions
     ---------------------- */

  /* Number of coals */
  cm->n_coals = cs_tree_get_node_count(tn_sf, "solid_fuel");
  if (cm->n_coals > CS_COMBUSTION_MAX_COALS)
    bft_error(__FILE__, __LINE__, 0,
              _("Coal number is limited to %i\n"
                "In the parametric file it is %i.\n"
                "Calculation is interupted. Check the parametric file.\n"),
              cm->n_coals, CS_COMBUSTION_MAX_COALS);

  /* Loop on coal nodes */
  int iclag  = 0;
  int idecal = 0;

  int icha = 0;

  for (cs_tree_node_t *tn = cs_tree_node_get_child(tn_sf, "solid_fuel");
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn), icha++) {

    int itypdp = _get_diameter_type(tn);
    cm->n_classes_per_coal[icha] = _get_nb_class(tn, itypdp);
    if (cm->n_classes_per_coal[icha] > CS_COMBUSTION_MAX_CLASSES_PER_COAL)
      bft_error(__FILE__, __LINE__, 0,
                _("class number by coal is limited.\n"
                  "For coal %i it is %i \n in the parametric file \n"),
                  icha, cm->n_classes_per_coal[icha]);

    /* Compute number of classes and fill ICHCOR */
    cm->nclacp = cm->nclacp + cm->n_classes_per_coal[icha];

    for (int iclapc = 0; iclapc < cm->n_classes_per_coal[icha]; iclapc++) {
      int icla = iclapc + idecal;
      cm->ichcor[icla] = icha + 1;
    }
    idecal += cm->n_classes_per_coal[icha];

    /* Diameter type = 1 ---> given
                     = 2 ---> Rosin-Rammler law */

    if (itypdp == 1) {

      /* Loop on classes */

      int icla = 0;
      for (cs_tree_node_t *tn_d = cs_tree_get_node(tn, "class/diameter");
           tn_d != NULL;
           tn_d = cs_tree_node_get_next_of_name(tn_d), icla++)
        cs_gui_node_get_real(tn_d, &(cm->diam20[icla + iclag]));

    }
    else if (itypdp == 2) {

      int nbrf = cs_tree_get_node_count(tn, "refusal");

      cs_real_t *dprefus, *refus, *pourc;
      BFT_MALLOC(dprefus, nbrf, cs_real_t);
      BFT_MALLOC(refus, nbrf, cs_real_t);
      BFT_MALLOC(pourc, cm->n_classes_per_coal[icha], cs_real_t);

      for (int ii = 0; ii < nbrf; ii++) {
        dprefus[ii] = -1;
        refus[ii] = -1;
        pourc[ii] = 0;
      }

      int ii = 0;
      for (cs_tree_node_t *tn_r = cs_tree_get_node(tn, "refusal");
           tn_r != NULL;
           tn_r = cs_tree_node_get_next_of_name(tn_r), ii++) {

        cs_real_t _dprefus = -1;
        cs_gui_node_get_child_real(tn_r, "diameter", &_dprefus);
        if (_dprefus >= 0) dprefus[ii] = _dprefus * 1.e6; /* in microns */

        cs_gui_node_get_child_real(tn_r, "value", &(refus[ii]));

      }
      assert(ii == nbrf);

      ii = 0;
      for (cs_tree_node_t *tn_m = cs_tree_get_node(tn, "class/mass_percent");
           tn_m != NULL;
           tn_m = cs_tree_node_get_next_of_name(tn_m), ii++)
        cs_gui_node_get_real(tn_m, &(pourc[ii]));

      /* split classes */

      cs_real_t *rf;
      BFT_MALLOC(rf, cm->n_classes_per_coal[icha], cs_real_t);
      rf[0] = pourc[0] / 2.;

      for (int icla = 1; icla < cm->n_classes_per_coal[icha]; icla++)
        rf[icla] = rf[icla-1] + (pourc[icla] + pourc[icla-1]) / 2.;

      cs_real_t kk1 = 0., kk2 = 0., kk3 = 0., kk4 = 0.;
      for (ii = 0; ii < nbrf ; ii++) {
        kk1 = kk1 + log(dprefus[ii]);
        kk2 = kk2 + log(-log(refus[ii]));
        kk3 = kk3 + log(dprefus[ii])*log(dprefus[ii]);
        kk4 = kk4 + log(dprefus[ii])*log(-log(refus[ii]));
      }

      cs_real_t denom = (nbrf * kk3 - kk1 * kk1);

      if (fabs(denom) < 1e-16)
        bft_error(__FILE__, __LINE__, 0,
                  _("Rosin-Rammler refusal parameters for coal %d lead to "
                    "%g denominator."), icha+1, denom);

      cs_real_t qq  = (nbrf * kk4 - kk1 * kk2) / (nbrf * kk3 - kk1 * kk1);
      cs_real_t var = (kk2 * kk3 - kk1 * kk4) / (nbrf * kk3 - kk1 * kk1);
      cs_real_t xx  = exp(-var / qq);

      for (int icla = iclag; icla < iclag + cm->n_classes_per_coal[icha]; icla++)
        cm->diam20[icla]=  xx*pow((-log(1.-rf[icla-iclag])),(1./qq))*1.e-6; // metres

      bft_printf("** Rosin-Rammler results for the coal %i **\n"
                 "[ Checking of the Rosin-Rammler law ]\n"
                 "Diameter       refus given      refus computed\n\n", icha+1);

      for (int icla = 0; icla< nbrf; icla++)
        bft_printf("%f     %f     %f \n",
                   dprefus[icla], refus[icla],
                   exp(-pow((dprefus[icla]/xx),(qq))));

      bft_printf("\nRefus       diam. given      diam. computed\n");

      for (int icla = 0; icla< nbrf; icla++)
        bft_printf("%f     %f     %f \n",
                   refus[icla], dprefus[icla],
                   xx*pow((-log(refus[icla])),(1./qq)));

      bft_printf("\nDiameters computed by the Rosin-Rammler law\n");

      for (int icla = iclag; icla < iclag+cm->n_classes_per_coal[icha]; icla ++)
        bft_printf("%d     %f \n", icla-iclag, cm->diam20[icla]);

      BFT_FREE(pourc);
      BFT_FREE(refus);
      BFT_FREE(dprefus);
      BFT_FREE(rf);

    }
    else {
      bft_error(__FILE__, __LINE__, 0,
                _("type diameter value must be equal to 1 or 2."));
    }

    iclag = iclag + cm->n_classes_per_coal[icha];

    /* Elementary composition in C, H , O , N , S on dry (% mass) */
    cm->cch[icha] = _get_solid_fuel_child_real(tn, "C_composition_on_dry");
    cm->hch[icha] = _get_solid_fuel_child_real(tn, "H_composition_on_dry");
    cm->och[icha] = _get_solid_fuel_child_real(tn, "O_composition_on_dry");
    cm->nch[icha] = _get_solid_fuel_child_real(tn, "N_composition_on_dry");
    cm->sch[icha] = _get_solid_fuel_child_real(tn, "S_composition_on_dry");

    /* Elementary composition in C, H , O , N , S od coke (% mass) */
    cm->cck[icha] = _get_solid_fuel_child_real(tn, "C_coke_composition_on_dry");
    cm->hck[icha] = _get_solid_fuel_child_real(tn, "H_coke_composition_on_dry");
    cm->ock[icha] = _get_solid_fuel_child_real(tn, "O_coke_composition_on_dry");
    cm->nck[icha] = _get_solid_fuel_child_real(tn, "N_coke_composition_on_dry");
    cm->sck[icha] = _get_solid_fuel_child_real(tn, "S_coke_composition_on_dry");

    /* PCI on dry or pure coal based on IPCI value */
    _get_pci_type_and_value(tn, &(cm->ipci[icha]), &(cm->pcich[icha]));

    cm->h0ashc[icha] = _get_solid_fuel_child_real(tn, "ashes_enthalpy");
    cm->cpashc[icha] = _get_solid_fuel_child_real(tn, "ashes_thermal_capacity");

    /*  ---- CP moyen du charbon sec (J/kg/K) */
    cm->cp2ch[icha] = _get_solid_fuel_child_real(tn, "specific_heat_average");

    /* ---- Masse volumique initiale (kg/m3) */
    cm->rho0ch[icha] = _get_solid_fuel_child_real(tn, "density");

    /* ---- Thermal conductivity of the coal (W/m/K) */
    if (cm->type == CS_COMBUSTION_COAL_WITH_DRYING) {
      cm->thcdch[icha] = 1e-5; /* default */
      cs_gui_node_get_child_real(tn, "thermal_conductivity",
                                 &(cm->thcdch[icha]));
    }

    /* Ash characteristics */

    /* Ash fraction (kg/kg) in % */
    cm->xashsec[icha] = _get_solid_fuel_child_real(tn, "rate_of_ashes_on_mass");

    /* Transformation to kg/kg */
    cm->xashch[icha] = cm->xashsec[icha]/100.;

    /* Humidity fraction (kg/kg) in % */
    cm->xwatch[icha] = _get_solid_fuel_child_real(tn, "moisture");

    /* Transformation to kg/kg */
    cm->xwatch[icha] = cm->xwatch[icha]/100.;

    /* Transform ash ratio from dry to humid in kg/kg */
    cm->xashch[icha] = cm->xashch[icha]*(1.-cm->xwatch[icha]);

    /* Devolatilisation parameters (Kobayashi's model) */

    const char path_dv[] = "devolatilisation_parameters";
    cs_tree_node_t *tn_dv = cs_tree_get_node(tn, path_dv);

    if (tn_dv == NULL)
      bft_error(__FILE__, __LINE__, 0, _("Missing %s child for node %s."),
                path_dv, tn->name);

    _get_y1y2_coefficient_values(tn_dv,
                                 &(cm->iy1ch[icha]), &(cm->iy2ch[icha]),
                                 &(cm->y1ch[icha]), &(cm->y2ch[icha]));

#if _XML_DEBUG_
    /* volatile_matter used by GUI for y1y2 automatic formula (to compute
       (compute y1ch[icha] and y2ch[icha] when iy1ch[icha] = 2 and
       iy1ch[icha] = 2), but not used directly here. */
    cs_real_t volatile_matter
      = _get_solid_fuel_child_real(tn, "volatile_matter");
#endif

    cm->a1ch[icha]
      = _get_solid_fuel_child_real(tn_dv, "A1_pre-exponential_factor");
    cm->a2ch[icha]
      = _get_solid_fuel_child_real(tn_dv, "A2_pre-exponential_factor");
    cm->e1ch[icha]
      = _get_solid_fuel_child_real(tn_dv, "E1_energy_of_activation");
    cm->e2ch[icha]
      = _get_solid_fuel_child_real(tn_dv, "E2_energy_of_activation");

    /* Char combustion parameters */

    const char path_cc[] = "char_combustion";
    cs_tree_node_t *tn_cc = cs_tree_get_node(tn, path_cc);

    if (tn_cc == NULL)
      bft_error(__FILE__, __LINE__, 0, _("Missing %s child for node %s."),
                path_cc, tn->name);

    /* Heterogeneous combustion parameters for O2 (shrinking sphere model) */
    cm->ahetch[icha]
      = _get_cc_specie_value(tn_cc, "O2", "pre-exponential_constant");
    cm->ehetch[icha]
      = _get_cc_specie_value(tn_cc, "O2", "energy_of_activation");
    cm->iochet[icha] = _get_cc_reaction_order(tn_cc, "O2");

    /* Heterogeneous combustion parameters for CO2 (shrinking sphere model) */
    if (cm->ihtco2) {
      cm->ahetc2[icha]
        = _get_cc_specie_value(tn_cc, "CO2", "pre-exponential_constant");
      cm->ehetc2[icha]
        = _get_cc_specie_value(tn_cc, "CO2", "energy_of_activation");
      cm->ioetc2[icha] = _get_cc_reaction_order(tn_cc, "CO2");
    }

    /* Heterogeneous combustion parameters for H2O (shrinking sphere model) */
    if (cm->ihth2o) {
      cm->ahetwt[icha]
        = _get_cc_specie_value(tn_cc, "H2O", "pre-exponential_constant");
      cm->ehetwt[icha]
        = _get_cc_specie_value(tn_cc, "H2O", "energy_of_activation");
      cm->ioetwt[icha] = _get_cc_reaction_order(tn_cc, "H2O");
    }

    /* NOX model parameters
       QPR =  % of free nitrogen during devolatilization
            / % of density freed during devolatilization */

    cs_gui_node_get_child_status_int(tn_sf, "NOx_formation", &(cm->ieqnox));
    if (cm->ieqnox) {

      const char path_nox[] = "nox_formation";
      cs_tree_node_t *tn_nox = cs_tree_get_node(tn, path_nox);

      if (tn_nox == NULL)
        bft_error(__FILE__, __LINE__, 0, _("Missing %s child for node %s."),
                  path_nox, tn->name);

      cm->qpr[icha] = _get_solid_fuel_child_real(tn_nox, "nitrogen_fraction");
      cm->fn[icha] = _get_solid_fuel_child_real(tn_nox, "nitrogen_concentration");

      /* Distribution of nitrogen between HCN and NH3 */

      /* Under "devotilisation_parameters" node */

      cm->crepn1[icha][0] = _get_solid_fuel_child_real
                              (tn_dv, "HCN_NH3_partitionning_reaction_1");
      cm->crepn1[icha][1] = 1 - cm->crepn1[icha][0];
      cm->crepn2[icha][0] = _get_solid_fuel_child_real
                              (tn_dv, "HCN_NH3_partitionning_reaction_2");
      cm->crepn2[icha][1] = 1 - cm->crepn2[icha][0];

      /* Under "nox_formation" node */

      cm->repnck[icha] = _get_solid_fuel_child_real
                           (tn_nox, "percentage_HCN_char_combustion");
      cm->repnle[icha] = _get_solid_fuel_child_real
                           (tn_nox, "nitrogen_in_char_at_low_temperatures");
      cm->repnlo[icha] = _get_solid_fuel_child_real
                          (tn_nox, "nitrogen_in_char_at_high_temperatures");
      cs_gui_node_get_child_status_int
        (tn_nox, "improved_NOx_model", &(cm->imdnox));
      if (cm->imdnox)
        _get_nox_reburning(tn_nox, &(cm->irb));
    }
    else {
      cm->crepn1[icha][0] = 0.5;
      cm->crepn1[icha][1] = 0.5;
      cm->crepn2[icha][0] = 0.5;
      cm->crepn2[icha][1] = 0.5;
    }
  }

  /* Oxidant definitions
     ------------------- */

  const char path_oxy[] = "thermophysical_models/solid_fuels/oxidants";
  cs_tree_node_t *tn_oxi = cs_tree_get_node(cs_glob_tree, path_oxy);

  /* Numer of oxydants */
  cm->noxyd = cs_tree_get_node_count(tn_oxi, "oxidant");
  if (cm->noxyd < 1 || cm->noxyd > 3 ) {
    bft_error(__FILE__, __LINE__, 0,
              _("Oxidant count must be between 1 and 3.\n"
                "It is %d in the current setup."),
              cm->noxyd);
  }
  int itypoxy = _get_oxidant_type(tn_oxi);

  /* Composition in O2, N2, H2O, N2 */

  _get_oxidants_composition(tn_oxi, "O2_composition", cm->oxyo2);
  _get_oxidants_composition(tn_oxi, "N2_composition", cm->oxyn2);
  _get_oxidants_composition(tn_oxi, "H2O_composition", cm->oxyh2o);
  _get_oxidants_composition(tn_oxi, "CO2_composition", cm->oxyco2);

  if (itypoxy == 1) {
    /* transformation pourcentage volumique en nombre de mole */
    for (int ioxy = 0; ioxy < cm->noxyd; ioxy++) {
      cs_real_t coef = 100.;
      if (cm->oxyo2[ioxy] > 0.)
        coef = CS_MIN(coef, cm->oxyo2[ioxy]);
      if (cm->oxyn2[ioxy] > 0.)
        coef = CS_MIN(coef, cm->oxyn2[ioxy]);
      if (cm->oxyh2o[ioxy] > 0.)
        coef = CS_MIN(coef, cm->oxyh2o[ioxy]);
      if (cm->oxyco2[ioxy] > 0.)
        coef = CS_MIN(coef, cm->oxyco2[ioxy]);

      cm->oxyo2 [ioxy] /= coef;
      cm->oxyn2 [ioxy] /= coef;
      cm->oxyh2o[ioxy] /= coef;
      cm->oxyco2[ioxy] /= coef;
    }
  }

  /* Numerical parameters
     -------------------- */

  {
    cs_tree_node_t *tn
      = cs_tree_get_node(cs_glob_tree,
                         "numerical_parameters/density_relaxation");

    cs_gui_node_get_real(tn, &(cm->srrom)); // inactive line if tn == NULL
  }
}

/*----------------------------------------------------------------------------
 * Gas combustion model parameters
 *----------------------------------------------------------------------------*/

void
cs_gui_combustion_gas_model(void)
{
  cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;
  assert(cm != NULL);

  int model_type = cm->type / 100;

  /* 3-point, EBU, and LW */

  if (model_type == 1 || model_type == 3 || model_type == 4) {
    cs_tree_node_t *tn
      = cs_tree_get_node(cs_glob_tree,
                         "numerical_parameters/density_relaxation");

    cs_gui_node_get_real(tn, &(cm->srrom)); // inactive line if tn == NULL

#if _XML_DEBUG_
    bft_printf("==> %s\n", __func__);
    bft_printf("--srrom  = %f\n", cm->srrom);
#endif
  }
}

/*----------------------------------------------------------------------------
 * Electrical model: read parameters
 *----------------------------------------------------------------------------*/

void
cs_gui_elec_model(void)
{
  cs_elec_option_t *elec_opt = cs_get_glob_elec_option();

  /* Numerical parameters */

  cs_tree_node_t *tn0
    = cs_tree_get_node(cs_glob_tree, "numerical_parameters/density_relaxation");

  cs_gui_node_get_real(tn0, &(elec_opt->srrom));

  /* Model parameters */

  tn0 = cs_tree_get_node(cs_glob_tree, "thermophysical_models/joule_effect");
  if (tn0 == NULL)
    return;

  cs_gui_node_get_child_status_int(tn0, "variable_scaling",
                                   &(elec_opt->ielcor));

  int ieljou = cs_glob_physical_model_flag[CS_JOULE_EFFECT];
  int ielarc = cs_glob_physical_model_flag[CS_ELECTRIC_ARCS];

  if (ieljou > 0)
    cs_gui_node_get_child_real(tn0, "imposed_power",
                               &(elec_opt->puisim));

  if (ielarc > 0) {
    cs_gui_node_get_child_real(tn0, "imposed_current",
                               &(elec_opt->couimp));

    if (cs_glob_elec_option->ielcor > 0) {

      cs_tree_node_t *tn_rc = cs_tree_get_node(tn0, "recal_model");

      const char *model = cs_gui_node_get_tag(tn_rc, "model");

      if (! strcmp(model, "general_case"))
        elec_opt->modrec = 1;
      else if (! strcmp(model, "plane_define"))
        elec_opt->modrec = 2;
      else if (! strcmp(model, "user"))
        elec_opt->modrec = 3;
      else
        bft_error(__FILE__, __LINE__, 0, _("Invalid model: %s"), model);

      if (cs_glob_elec_option->modrec == 2) {

        const char *dir_s = cs_tree_node_get_child_value_str(tn_rc, "direction");
        if (cs_gui_strcmp(dir_s, "X"))
          elec_opt->idreca = 1;
        else if (cs_gui_strcmp(dir_s, "Y"))
          elec_opt->idreca = 2;
        else
          elec_opt->idreca = 3;

        cs_tree_node_t *tn_pl = cs_tree_node_get_child(tn_rc, "plane_definition");

        const char *key[] = {"A", "B", "C", "D", "epsilon"};
        for (int i = 0; i < 5; i++)
          cs_gui_node_get_child_real(tn_pl, key[i],
                                     &(elec_opt->crit_reca[i]));
      }
    }
  }

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--ielcor  = %i\n", cs_glob_elec_option->ielcor);
  bft_printf("--puisim  = %f\n", cs_glob_elec_option->puisim);
  bft_printf("--couimp  = %f\n", cs_glob_elec_option->couimp);
  bft_printf("--modrec  = %d\n", cs_glob_elec_option->modrec);
  bft_printf("--srrom   = %d\n", cs_glob_elec_option->srrom);
#endif
}

/*----------------------------------------------------------------------------
 * Electrical model: define plane for elreca
 *
 * Fortran Interface:
 *----------------------------------------------------------------------------*/

void
cs_gui_elec_model_rec(void)
{
  /* build list of interior faces */

  cs_lnum_t   n_selected_faces = 0;
  cs_lnum_t  *selected_faces = NULL;

  char crit[128] = "";

  cs_elec_option_t *elec_opt = cs_get_glob_elec_option();

  snprintf(crit, 127, "plane[%f, %f, %f, %f, epsilon=%6f]",
           elec_opt->crit_reca[0],
           elec_opt->crit_reca[1],
           elec_opt->crit_reca[2],
           elec_opt->crit_reca[3],
           elec_opt->crit_reca[4]);
  crit[127] = '\0';

  BFT_MALLOC(selected_faces, cs_glob_mesh->n_i_faces, cs_lnum_t);

  cs_selector_get_i_face_list(crit,
                              &n_selected_faces,
                              selected_faces);

  for (cs_lnum_t j = 0; j < n_selected_faces; j++)
    elec_opt->izreca[selected_faces[j]] = 1;

  BFT_FREE(selected_faces);
}

/*-----------------------------------------------------------------------------
 * Return the name of a thermophysical model.
 *
 * parameter:
 *   model_thermo  <--  thermophysical model to check
 *----------------------------------------------------------------------------*/

const char *
cs_gui_get_thermophysical_model(const char  *model_thermo)
{
  const char *retval = NULL;

  cs_tree_node_t *tn = cs_tree_get_node(cs_glob_tree, "thermophysical_models");
  tn = cs_tree_node_get_child(tn, model_thermo);

  if (tn != NULL) {
    if (! strcmp(model_thermo, "gas_combustion"))
      retval = cs_tree_node_get_tag(tn, "option");
    else
      retval = cs_tree_node_get_tag(tn, "model");
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * groundwater model : read parameters
 *
 * parameters:
 *   permeability    <--   permeability type
 *   unsteady        <--   steady flow
 *   unsaturated     <--   take into account unsaturated zone
 *----------------------------------------------------------------------------*/

void
cs_gui_gwf_model(int  *permeability,
                 int  *unsteady,
                 int  *unsaturated)
{
  cs_tree_node_t *tn0
    = cs_tree_get_node(cs_glob_tree, "thermophysical_models/groundwater_model");

  const char *mdl;

  /* Get flow type */
  mdl = cs_tree_node_get_tag(cs_tree_node_get_child(tn0, "flowType"),
                             "model");

  if (cs_gui_strcmp(mdl, "steady"))
    *unsteady = 0;
  else
    *unsteady = 1;

  /* Get permeability type */
  mdl = cs_tree_node_get_tag(cs_tree_node_get_child(tn0, "permeability"),
                             "model");

  if (cs_gui_strcmp(mdl, "anisotropic"))
    *permeability = 1;
  else
    *permeability = 0;

  /* Get the possible presence of unsaturated zone */
  mdl = cs_tree_node_get_tag(cs_tree_node_get_child(tn0, "unsaturatedZone"),
                             "model");

  if (cs_gui_strcmp(mdl, "true"))
    *unsaturated = 1;
  else
    *unsaturated = 0;

  /* Get first-order decay rate and chemistry model */

  const int key_decay = cs_field_key_id("fo_decay_rate");

  for (cs_tree_node_t *tn = cs_tree_get_node(tn0, "scalar");
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn)) {

    const char *name = cs_gui_node_get_tag(tn, "name");

    cs_field_t *f = cs_field_by_name_try(name);
    if (f == NULL) continue;

    if (    (f->type & CS_FIELD_VARIABLE)
         && (f->type & CS_FIELD_USER)) {

      /* get first-order decay rate */

      cs_real_t decay = cs_field_get_key_double(f, key_decay);
      cs_gui_node_get_child_real(tn, "fo_decay_rate", &decay);
      cs_field_set_key_double(f, key_decay, decay);

      /* get chemistry model */

      // const char *cmodel = cs_tree_node_get_tag(tn, "chemistry_model");

      /* TODO update for CDO-based GWF */
    }
  }

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--groundwater_anisotropic_permeability  = %d\n", *permeability);
  bft_printf("--groundwater_unsteady                  = %d\n", *unsteady);
  bft_printf("--groundwater_unsaturated               = %d\n", *unsaturated);
#endif
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
