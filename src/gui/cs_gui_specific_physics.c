/*============================================================================
 * Management of the GUI parameters file: specific physics
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include "mei_evaluate.h"

#include "cs_base.h"
#include "cs_gui_util.h"
#include "cs_gui.h"
#include "cs_gui_variables.h"
#include "cs_mesh.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_physical_model.h"
#include "cs_prototypes.h"
#include "cs_selector.h"
#include "cs_elec_model.h"
#include "cs_gwf_physical_properties.h"

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
#define LG_MAX 1000

/* debugging switch */
#define _XML_DEBUG_ 0

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Return the label or the name from a specific physics scalar.
 *
 * parameters:
 *   f        <-- pointer to field structure
 *   physics  <-- keyword: specific physic model required
 *   kw       <-- keyword: name of scalar
 *----------------------------------------------------------------------------*/

static void
_set_scalar_name_label(cs_field_t  *f,
                       const char  *physics,
                       const char  *kw)
{
  char *path = NULL;
  char *label  = NULL;

  const int klbl = cs_field_key_id("label");

  path = cs_xpath_short_path();
  cs_xpath_add_elements(&path, 3,
                        "thermophysical_models",
                         physics,
                        "variable");

  cs_xpath_add_test_attribute(&path, "name", kw);
  cs_xpath_add_attribute(&path, "label");

  label = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  cs_field_set_key_str(f, klbl, label);

  BFT_FREE(label);
}

/*-----------------------------------------------------------------------------
 * Return the name tor thermal scalar.
 *
 * parameters:
 *   f  <-- pointer to field structure
 *   kw <-- scalar name
 *----------------------------------------------------------------------------*/

static void
_set_thermal_scalar_name_label(cs_field_t  *f,
                               const char  *kw)
{
  char *path = NULL;
  char *label  = NULL;

  const int klbl = cs_field_key_id("label");

  path = cs_xpath_short_path();
  cs_xpath_add_elements(&path, 3,
                        "thermophysical_models",
                        "thermal_scalar",
                        "variable");
  cs_xpath_add_attribute(&path, kw);

  label = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  cs_field_set_key_str(f, klbl, label);

  BFT_FREE(label);
}

/*-----------------------------------------------------------------------------
 * Return an integer for class number.
 *
 * parameters:
 *    icha    -->   char number
 *    type    -->   type of diameter definition
 *
 * returns:
 *    number of class
 *----------------------------------------------------------------------------*/

static int
_cs_gui_get_nb_class(const int icha,
                     const int type)
{
  char *path;
  int value = 0;
  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  if (type == 1)
    cs_xpath_add_elements(&path, 2, "class", "diameter");
  else if (type == 2)
    cs_xpath_add_elements(&path, 2, "class", "mass_percent");
  value = cs_gui_get_nb_element(path);
  BFT_FREE(path);

  return value;
}

/*-----------------------------------------------------------------------------
 * Return an integer for the refusal number.
 *
 * parameters:
 *    icha    -->   char number
 *
 * returns:
 *    number of refusal
 *----------------------------------------------------------------------------*/

static int
_cs_gui_get_nb_refusal(const int icha)
{
  char *path;
  int value = 0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_element(&path, "refusal");
  value = cs_gui_get_nb_element(path);
  BFT_FREE(path);

  return value ;
}

/*-----------------------------------------------------------------------------
 * Return integer for diameter type.
 *
 * parameters:
 *    icha     -->  char number
 *
 * returns:
 *    integer for diameter type
 *----------------------------------------------------------------------------*/

static int
_get_diameter_type(const int icha)
{
  char *path;
  char *buff = NULL;
  int   ichoice = 0;
  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_element(&path, "diameter_type");
  cs_xpath_add_function_text(&path);
  buff = cs_gui_get_text_value(path);
  if (buff == NULL) {
    ichoice = 1;
  } else {
    if (cs_gui_strcmp(buff, "automatic"))
      ichoice = 1;
    else if (cs_gui_strcmp(buff, "rosin-rammler_law"))
      ichoice = 2;
    else
      bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);
  }
  BFT_FREE(path);
  BFT_FREE(buff);

  return ichoice;
}

/*-----------------------------------------------------------------------------
 * Return double for diameter value.
 *
 * parameters:
 *    icha     -->  char number
 *    icla     -->  class number
 *
 * returns:
 *    diameter value
 *----------------------------------------------------------------------------*/

static double
_get_solid_fuel_diameter(const int icha,
                         const int icla)
{
  char *path;
  double result;
  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_element(&path, "class");
  cs_xpath_add_element_num(&path, "diameter", icla);
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);
  BFT_FREE(path);

  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for refusal diameter.
 *
 * parameters:
 *    ii       -->  refusal number
 *    icha     -->  char number
 *
 * returns:
 *    refusal diameter
 *----------------------------------------------------------------------------*/

static double
_get_solid_fuel_refusal_diameter(const int ii,
                                 const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_element_num(&path, "refusal", ii);
  cs_xpath_add_element(&path, "diameter");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for refusal value.
 *
 * parameters:
 *    ii       -->  refusal number
 *    icha     -->  char number
 *
 * returns:
 *    refusal value
 *----------------------------------------------------------------------------*/

static double
_get_solid_fuel_refusal_value(const int ii,
                              const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_element_num(&path, "refusal", ii);
  cs_xpath_add_element(&path, "value");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for refusal percent.
 *
 * parameters:
 *    ii       -->  refusal number
 *    icha     -->  char number
 *
 * returns:
 *    refusal value
 *----------------------------------------------------------------------------*/

static double
_get_solid_fuel_pourc(const int ii,
                      const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_element(&path, "class");
  cs_xpath_add_element_num(&path, "mass_percent", ii);
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for the composition on dry.
 *
 * parameters:
 *    icha     -->   char number
 *    type     <--   type of elementaire
 * returns:
 *    composition on dry value
 *----------------------------------------------------------------------------*/

static double
_get_solid_fuel_composition_on_dry(const int icha,
                                   const char *const type)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_element(&path, type);
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for specific heat average.
 *
 * parameters:
 *    icha     -->   char number
 * returns:
 *    specific heat average
 *----------------------------------------------------------------------------*/

static double
_get_solid_fuel_specific_heat_average(const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_element(&path, "specific_heat_average");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for solid fuel density.
 *
 * parameters:
 *    icha     -->   char number
 * returns:
 *    density
 *----------------------------------------------------------------------------*/

static double
_get_solid_fuel_density(const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_element(&path, "density");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for solid fuel thermal conductivity
 *
 * parameters:
 *    icha     -->   char number
 * returns:
 *    thermal conductivity
 *----------------------------------------------------------------------------*/

static double
_get_solid_fuel_thermal_conductivity(const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_element(&path, "thermal_conductivity");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result)) {
    result = 1.e-5; /* default, used with Lagrangian model only */
  }

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for ashes rate
 *
 * parameters:
 *    icha     -->   char number
 * returns:
 *    ashes rate
 *----------------------------------------------------------------------------*/

static double
_get_ashes_rate(const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_element(&path, "rate_of_ashes_on_mass");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for humidity rate.
 *
 * parameters:
 *    icha     -->   char number
 * returns:
 *    humidity rate
 *----------------------------------------------------------------------------*/

static double
_get_humidity_rate(const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_element(&path, "moisture");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for volatile matter.
 *
 * parameters:
 *    icha     -->   char number
 * returns:
 *    volatile matter rate
 *----------------------------------------------------------------------------*/

static double
_get_volatile_matter(const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_element(&path, "volatile_matter");
  cs_xpath_add_function_text(&path);
  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for ashes enthalpy.
 *
 * parameters:
 *    icha     -->   char number
 * returns:
 *    ashes enthalpy
 *----------------------------------------------------------------------------*/

static double
_get_ashes_enthalpy(const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_element(&path, "ashes_enthalpy");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for ashes thermal capacity.
 *
 * parameters:
 *    icha     -->   char number
 * returns:
 *    ashes thermal capacity
 *----------------------------------------------------------------------------*/

static double
_get_ashes_thermal_capacity(const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_element(&path, "ashes_thermal_capacity");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return integer for PCI type.
 *
 * parameters:
 *    icha     -->   char number
 * returns:
 *    integer for PCI type
 *----------------------------------------------------------------------------*/

static int
_get_PCI_type(int icha)
{
  char *path;
  char *path2;
  char *buff = NULL;
  char *buff2 = NULL;
  int ichoice = 0.;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_element(&path, "Heating_model");
  cs_xpath_add_attribute(&path,"choice");
  buff = cs_gui_get_attribute_value(path);
  if (buff != NULL) {
    if (cs_gui_strcmp(buff, "LHV")) {
      path2 = cs_xpath_init_path();
      cs_xpath_add_elements(&path2, 2,"thermophysical_models", "solid_fuels");;
      cs_xpath_add_element_num(&path2, "solid_fuel", icha);
      cs_xpath_add_element(&path2, "Heating_model");
      cs_xpath_add_element(&path2, "type");
      cs_xpath_add_function_text(&path2);
      buff2 = cs_gui_get_text_value(path2);
      if (buff2 == NULL) {
        ichoice = 1;
      } else {
        if (cs_gui_strcmp(buff2, "dry_basis"))
          ichoice = 1;
        else if (cs_gui_strcmp(buff2, "dry_ash_free"))
          ichoice = 0;
        else if (cs_gui_strcmp(buff2, "as_received"))
          ichoice = 2;
        else
          bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path2);

        BFT_FREE(buff2);
      }
      BFT_FREE(path2);
    }
    else if (cs_gui_strcmp(buff, "HHV")) {
      path2 = cs_xpath_init_path();
      cs_xpath_add_elements(&path2, 2,"thermophysical_models", "solid_fuels");
      cs_xpath_add_element_num(&path2, "solid_fuel", icha);
      cs_xpath_add_element(&path2, "Heating_model");
      cs_xpath_add_element(&path2, "type");
      cs_xpath_add_function_text(&path2);
      buff2 = cs_gui_get_text_value(path2);
      if (buff2 == NULL) {
        ichoice = 4;
      } else {
        if (cs_gui_strcmp(buff2, "dry_basis"))
          ichoice = 4;
        else if (cs_gui_strcmp(buff2, "dry_ash_free"))
          ichoice = 3;
        else if (cs_gui_strcmp(buff2, "as_received"))
          ichoice = 5;
        else
          bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path2);

        BFT_FREE(buff2);
      }
      BFT_FREE(path2);
    }
    else if (cs_gui_strcmp(buff, "IGT_correlation"))
      ichoice = 6;
    else
      bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

    BFT_FREE(buff);
  }
  BFT_FREE(path);
  return ichoice;
}

/*-----------------------------------------------------------------------------
 * Return integer for Y1/Y2 coeffficient type.
 *
 * parameters:
 *    icha     -->   char number
 * returns:
 *    integer for Y1/Y2 coefficient type
 *----------------------------------------------------------------------------*/

static int
_get_Y1Y2_coefficient_type(const int icha)
{
  char *path;
  char *buff = NULL;
  int ichoice = 0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_elements(&path, 2,"devolatilisation_parameters",
                        "stoichiometric_coefficient");
  cs_xpath_add_attribute(&path,"type");
  buff = cs_gui_get_attribute_value(path);
  if (buff != NULL) {
    if (cs_gui_strcmp(buff, "automatic_CHONS"))
      ichoice = 0;
    else if (cs_gui_strcmp(buff, "user_define"))
      ichoice = 1;
    else if (cs_gui_strcmp(buff, "automatic_formula"))
      ichoice = 2;
    else
      bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);
  }
  BFT_FREE(path);
  BFT_FREE(buff);
  return ichoice;
}

/*-----------------------------------------------------------------------------
 * Return double for PCI value.
 *
 * parameters:
 *    icha     -->   char number
 * returns:
 *    PCI value
 *----------------------------------------------------------------------------*/

static double
_get_PCI_value(const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_element(&path,"Heating_model");
  cs_xpath_add_element(&path,"value");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for Y1 coefficient.
 *
 * parameters:
 *    icha     -->   char number
 * returns:
 *    Y1 coefficient
 *----------------------------------------------------------------------------*/

static double
_get_Y1_coefficient(const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_elements(&path, 3,"devolatilisation_parameters",
                        "stoichiometric_coefficient", "Y1");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for Y2 coefficient.
 *
 * parameters:
 *    icha     -->   char number
 * returns:
 *    Y2 coefficient
 *----------------------------------------------------------------------------*/

static double
_get_Y2_coefficient(const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_elements(&path, 3,"devolatilisation_parameters",
                        "stoichiometric_coefficient", "Y2");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for pre exponential factor.
 *
 * parameters:
 *    icha     -->   char number
 * returns:
 *    A1 pre exponentia factor
 *----------------------------------------------------------------------------*/

static double
_get_A1_pre_exponential_factor(const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_elements(&path, 2,"devolatilisation_parameters",
                        "A1_pre-exponential_factor");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for A2 pre exponential factor.
 *
 * parameters:
 *    icha     -->   char number
 * returns:
 *    A2 pre exponential factor
 *----------------------------------------------------------------------------*/

static double
_get_A2_pre_exponential_factor(const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_elements(&path, 2,"devolatilisation_parameters",
                        "A2_pre-exponential_factor");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for E1 energy of activation.
 *
 * parameters:
 *    icha     -->   char number
 * returns:
 *    E1 energy of activation
 *----------------------------------------------------------------------------*/

static double
_get_E1_energy_of_activation(const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_elements(&path, 2,"devolatilisation_parameters",
                        "E1_energy_of_activation");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for E2 energy of activation.
 *
 * parameters:
 *    icha     -->   char number
 * returns:
 *    E2 energy of activation
 *----------------------------------------------------------------------------*/

static double
_get_E2_energy_of_activation(const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_elements(&path, 2,"devolatilisation_parameters",
                        "E2_energy_of_activation");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for partition in HCN NO reaction number 1.
 *
 * parameters:
 *    icha     -->   char number
 * returns:
 *    partition in HCN/NO reaction number 1
 *----------------------------------------------------------------------------*/

static double
_get_nitrogen_partition_in_HCN_NH3_reaction_1(const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_elements(&path, 2, "devolatilisation_parameters",
                        "HCN_NH3_partitionning_reaction_1");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for partition in HCN/NO reaction number 2.
 *
 * parameters:
 *    icha     -->   char number
 * returns:
 *    partition in HCN/NO reaction number 2
 *----------------------------------------------------------------------------*/

static double
_get_nitrogen_partition_in_HCN_NH3_reaction_2(const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_elements(&path, 2, "devolatilisation_parameters",
                        "HCN_NH3_partitionning_reaction_2");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for nitrogen fraction.
 *
 * parameters:
 *    icha     -->   char number
 * returns:
 *    nitrogen fraction
 *----------------------------------------------------------------------------*/

static double
_get_nitrogen_fraction(const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_elements(&path, 2, "nox_formation", "nitrogen_fraction");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for nitrogen concentration.
 *
 * parameters:
 *    param    -->   name of parameter
 * returns:
 *    nitrogen concentration
 *----------------------------------------------------------------------------*/

static double
_get_nitrogen_concentration(const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_elements(&path, 2, "nox_formation", "nitrogen_concentration");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for nitrogen in char.
 *
 * parameters:
 *    param    -->   name of parameter
 * returns:
 *    nitrogen concentration
 *----------------------------------------------------------------------------*/

static double
_get_hcn_char_comb(const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_elements(&path, 2, "nox_formation", "percentage_HCN_char_combustion");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for nitrogen in char at low temperatures.
 *
 * parameters:
 *    param    -->   name of parameter
 * returns:
 *    nitrogen concentration
 *----------------------------------------------------------------------------*/

static double
_get_nitrogen_char_low_temp(const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_elements(&path, 2, "nox_formation", "nitrogen_in_char_at_low_temperatures");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for nitrogen in char at high temperatures.
 *
 * parameters:
 *    param    -->   name of parameter
 * returns:
 *    nitrogen concentration
 *----------------------------------------------------------------------------*/

static double
_get_nitrogen_char_high_temp(const int icha)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_elements(&path, 2, "nox_formation", "nitrogen_in_char_at_high_temperatures");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for pre exponential constant.
 *
 * parameters:
 *    icha     -->   char number
 *    species  -->   species choice
 * returns:
 *    pre exponential constant
 *----------------------------------------------------------------------------*/

static double
_get_pre_exponential_constant(const int icha,
                              const char *const species)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path,2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_elements(&path,2,"char_combustion","specie");
  cs_xpath_add_test_attribute(&path, "nature", species);
  cs_xpath_add_element(&path,"pre-exponential_constant");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return double for energy of activation.
 *
 * parameters:
 *    icha     -->   char number
 *    species  -->   species choice
 * returns:
 *    energy of activation
 *----------------------------------------------------------------------------*/

static double
_get_energy_of_activation(const int icha,
                          const char *const species)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path,2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_elements(&path,2,"char_combustion","specie");
  cs_xpath_add_test_attribute(&path, "nature", species);
  cs_xpath_add_element(&path,"energy_of_activation");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return integer for order of reaction.
 *
 * parameters:
 *    icha     -->   char number
 *    species  -->   species choice
 * returns:
 *    order of reaction
 *----------------------------------------------------------------------------*/

static int
_get_order_of_reaction(const int icha,
                       const char *const species)
{
  char *path;
  char *buff = NULL;
  int   ichoice = 0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path,2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_elements(&path,2,"char_combustion","specie");
  cs_xpath_add_test_attribute(&path, "nature", species);
  cs_xpath_add_element(&path,"order_of_reaction");
  cs_xpath_add_attribute(&path,"choice");
  buff = cs_gui_get_attribute_value(path);
  if (buff != NULL) {
    if (cs_gui_strcmp(buff, "0.5"))
      ichoice = 0;
    else if (cs_gui_strcmp(buff, "1"))
      ichoice = 1;
    else
      bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);
  }
  BFT_FREE(path);
  BFT_FREE(buff);
  return ichoice;
}

/*-----------------------------------------------------------------------------
 * Return double for oxydant composition.
 *
 *   parameters:
 *    icha     -->   char number
 *    species  -->   species choice
 * returns:
 *    oxydant composition
 *----------------------------------------------------------------------------*/

static double
_get_oxydant_composition(const int icha,
                         const char *const species)
{
  char *path;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path,3,
                        "thermophysical_models",
                        "solid_fuels",
                        "oxidants");
  cs_xpath_add_element_num(&path, "oxidant", icha);
  cs_xpath_add_element(&path, species);
  cs_xpath_add_function_text(&path);
  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return integer parameter for oxydant type.
 *
 * parameters:
 *    ichoice  <--   value of parameter
 * returns:
 *    parameter for oxydant type
 *----------------------------------------------------------------------------*/

static int
_get_oxydant_type(void) /*ici un char*/
{
  char *path;
  char *buff = NULL;
  int   ichoice = 0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path,3,
                        "thermophysical_models",
                        "solid_fuels",
                        "oxidants");
  cs_xpath_add_element(&path,"oxidant_type");
  cs_xpath_add_function_text(&path);
  buff = cs_gui_get_text_value(path);
  if (buff != NULL) {
    if (cs_gui_strcmp(buff, "volumic_percent"))
      ichoice = 0;
    else if (cs_gui_strcmp(buff, "molar"))
      ichoice = 1;
    else
      bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);
  }
  BFT_FREE(path);
  BFT_FREE(buff);
  return ichoice;
}

/*-----------------------------------------------------------------------------
 * Return double for absorption coefficient.
 *
 * parameters:
 *    result  <--   absorption coefficient
 * returns:
 *    absorption coefficient
 *----------------------------------------------------------------------------*/

static double
_get_absorption_coefficient(void) /*ici un char*/
{
  char *path;
  double   result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path,2,"thermophysical_models", "radiative_transfer");
  cs_xpath_add_element(&path,"absorption_coefficient");
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);
  BFT_FREE(path);
  return result;
}

/*-----------------------------------------------------------------------------
 * Return integer for CO2 kinetics status.
 *
 * parameters:
 *    icha     -->   char number
 *    keyword  <--  value of attribute node
 *----------------------------------------------------------------------------*/

static void
_getCO2KineticsStatus(int    *keyword)
{
  char *path   = NULL;
  int   status = 0;

  path = cs_xpath_init_path();

  cs_xpath_add_elements(&path, 3,
                        "thermophysical_models",
                        "solid_fuels",
                        "CO2_kinetics");

  cs_xpath_add_attribute(&path, "status");

  if (cs_gui_get_status(path, &status))
    *keyword = status;
  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Return integer for H2O kinetics status
 *
 * parameters:
 *    icha     -->   char number
 *    keyword  <--  value of attribute node
 *----------------------------------------------------------------------------*/

static void
_getH2OKineticsStatus(int   *keyword)
{
  char *path   = NULL;
  int   status = 0;

  path = cs_xpath_init_path();


  cs_xpath_add_elements(&path, 3,
                        "thermophysical_models",
                        "solid_fuels",
                        "H2O_kinetics");

  cs_xpath_add_attribute(&path, "status");

  if (cs_gui_get_status(path, &status))
    *keyword = status;
  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Return integer for NOx combustion status
 *
 *   parameters:
 *    icha     -->   char number
 *    keyword  <--  value of attribute node
 *----------------------------------------------------------------------------*/

static void
_getNOxStatus(int   *keyword)
{
  char *path   = NULL;
  int   status = 0;

  path = cs_xpath_init_path();

  cs_xpath_add_elements(&path,2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element(&path, "NOx_formation");

  cs_xpath_add_attribute(&path, "status");

  if (cs_gui_get_status(path, &status))
    *keyword = status;
  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Return integer for NOx feature status
 *
 *   parameters:
 *    icha     -->   char number
 *    keyword  <--  value of attribute node
 *----------------------------------------------------------------------------*/

static void
_getNOxFeatureStatus(const int icha,
                     int   *keyword)
{
  char *path   = NULL;
  int   status = 0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_elements(&path, 2, "nox_formation", "improved_NOx_model");

  cs_xpath_add_attribute(&path, "status");

  if (cs_gui_get_status(path, &status))
    *keyword = status;
  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Return integer for NOx feature status
 *
 *   parameters:
 *    icha     -->   char number
 *    keyword  <--  value of attribute node
 *----------------------------------------------------------------------------*/

static void
_getNOxReburning(const int icha,
                 int   *keyword)
{
  char *path   = NULL;
  char *choice = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_elements(&path, 2, "nox_formation", "reburning_model");

  cs_xpath_add_function_text(&path);
  choice = cs_gui_get_text_value(path);
  if (cs_gui_strcmp(choice, "unused"))
      *keyword = 0;
  else if (cs_gui_strcmp(choice, "chen"))
      *keyword = 1;
  else
      *keyword = 2;
  BFT_FREE(path);
  BFT_FREE(choice);
}

/*-----------------------------------------------------------------------------
 * Return the name of a joule effect model.
 *
 * parameter:
 *----------------------------------------------------------------------------*/

static char*
_get_joule_model(void)
{
  char *model = NULL;
  char *path = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "thermophysical_models",
                                  "joule_effect",
                                  "joule_model");
  cs_xpath_add_attribute(&path, "model");

  model = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return model;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Density under relaxation
 *
 * Fortran Interface:
 *
 * SUBROUTINE UICPI1 (SRROM, DIFTL0)
 * *****************
 * DOUBLE PRECISION SRROM   <--   density relaxation
 * DOUBLE PRECISION DIFTL0  <--   dynamic diffusion
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicpi1, UICPI1) (double *const srrom,
                                double *const diftl0)
{
  cs_var_t  *vars = cs_glob_var;

  cs_gui_numerical_double_parameters("density_relaxation", srrom);

  if (cs_gui_strcmp(vars->model, "gas_combustion") ||
      cs_gui_strcmp(vars->model, "solid_fuels")) {
    cs_gui_properties_value("dynamic_diffusion", diftl0);
  }

#if _XML_DEBUG_
  bft_printf("==>UICPI1\n");
  bft_printf("--srrom  = %f\n", *srrom);
  if (cs_gui_strcmp(vars->model, "gas_combustion")) {
    bft_printf("--diftl0  = %f\n", *diftl0);
  }
#endif
}

/*----------------------------------------------------------------------------
 * Temperature for D3P Gas Combustion
 *
 * Fortran Interface:
 *
 * SUBROUTINE UICPI2 (Toxy, Tfuel)
 * *****************
 * DOUBLE PRECISION Toxy   <--   Oxydant temperature
 * DOUBLE PRECISION Tfuel  <--   Fuel temperature
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicpi2, UICPI2) (double *const toxy,
                                double *const tfuel)
{
  cs_gui_reference_initialization("oxydant_temperature", toxy);
  cs_gui_reference_initialization("fuel_temperature", tfuel);
#if _XML_DEBUG_
  bft_printf("==>UICPI2\n");
  bft_printf("--toxy  = %f\n", *toxy);
  bft_printf("--tfuel  = %f\n", *tfuel);
#endif
}


/*----------------------------------------------------------------------------
 * Atmospheric flows: read of meteorological file of data
 *
 * Fortran Interface:
 *
 * subroutine uiati1
 * *****************
 * integer         imeteo   <--   on/off index
 * char(*)         fmeteo   <--   meteo file name
 * int             len      <--   meteo file name destination string length
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiati1, UIATI1) (int           *imeteo,
                                char          *fmeteo,
                                int           *len
                                CS_ARGF_SUPP_CHAINE)
{
  char *path   = NULL;
  int   status = 0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "thermophysical_models",
                                  "atmospheric_flows",
                                  "read_meteo_data");

  cs_xpath_add_attribute(&path, "status");
  if (cs_gui_get_status(path, &status))
    *imeteo = status;
  BFT_FREE(path);

  if (*imeteo) {

    int i, l;
    char *cstr = NULL;

    path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 3, "thermophysical_models",
                          "atmospheric_flows",
                          "meteo_data");

    cs_xpath_add_function_text(&path);
    cstr = cs_gui_get_text_value(path);

    BFT_FREE(path);

    /* Copy string */

    if (cstr != NULL) {

      /* Compute string length (removing start or end blanks) */

      l = strlen(cstr);
      if (l > *len)
        l = *len;

      for (i = 0; i < l; i++)
        fmeteo[i] = cstr[i];

      /* Pad with blanks if necessary */

      for (i = l; i < *len; i++)
        fmeteo[i] = ' ';

      BFT_FREE(cstr);

    }

  }

#if _XML_DEBUG_
  bft_printf("==>UIATI1\n");
  bft_printf("--imeteo  = %i\n", *imeteo);
#endif
}


/*----------------------------------------------------------------------------
 * Indirection between the solver numbering and the XML one
 * for physical properties of the activated specific physics
 * (pulverized solid fuels)
 *----------------------------------------------------------------------------*/

void CS_PROCF (uisofu, UISOFU) (const int    *const iirayo,
                                const int    *const iihmpr,
                                const int    *const ncharm,
                                int          *const ncharb,
                                int          *const nclpch,
                                int          *const nclacp,
                                const int    *const ncpcmx,
                                int          *const ichcor,
                                double       *const diam20,
                                double       *const cch,
                                double       *const hch,
                                double       *const och,
                                double       *const nch,
                                double       *const sch,
                                int          *const ipci,
                                double       *const pcich,
                                double       *const cp2ch,
                                double       *const rho0ch,
                                double       *const thcdch,
                                double       *const cck,
                                double       *const hck,
                                double       *const ock,
                                double       *const nck,
                                double       *const sck,
                                double       *const xashch,
                                double       *const xashsec,
                                double       *const xwatch,
                                double       *const h0ashc,
                                double       *const cpashc,
                                int          *const iy1ch,
                                double       *const y1ch,
                                int          *const iy2ch,
                                double       *const y2ch,
                                double       *const a1ch,
                                double       *const a2ch,
                                double       *const e1ch,
                                double       *const e2ch,
                                double       *const crepn1,
                                double       *const crepn2,
                                double       *const ahetch,
                                double       *const ehetch,
                                int          *const iochet,
                                double       *const ahetc2,
                                double       *const ehetc2,
                                int          *const ioetc2,
                                double       *const ahetwt,
                                double       *const ehetwt,
                                int          *const ioetwt,
                                int          *const ieqnox,
                                int          *const imdnox,
                                int          *const irb,
                                int          *const ihtco2,
                                int          *const ihth2o,
                                double       *const qpr,
                                double       *const fn,
                                double       *const ckabs1,
                                int          *const noxyd,
                                double       *const oxyo2,
                                double       *const oxyn2,
                                double       *const oxyh2o,
                                double       *const oxyco2,
                                double       *const repnck,
                                double       *const repnle,
                                double       *const repnlo)
{
  char *model = NULL;
  int iclag;
  int idecal;
  int itypdp;
  int itypoxy;
  int icha;
  int iclapc;
  int icla;
  int ii;
  double *pourc = NULL;
  double *refus= NULL;
  double *rf = NULL;
  double kk1, kk2, kk3, kk4;
  double qq, var, xx;
  int ioxy;
  double coef;
  int nbrf;
  double *dprefus= NULL;
  double volatile_matter = 0.;

  cs_var_t  *vars = cs_glob_var;

  if (*iihmpr != 1)
    cs_gui_load_file("dp_FCP.xml");

  /* ---- Nb de charbons */
  *ncharb = cs_gui_get_tag_count("/solid_fuels/solid_fuel", 1);
  if (*ncharb > *ncharm)
    bft_error(__FILE__, __LINE__, 0,
              _("Coal number is limited to %i\n"
                "In the parametric file it is %i.\n"
                "Calculation is interupted. Check the parametric file.\n"),
              *ncharb, *ncharm);

  /* --- BCLE SUR LES CHARBONS */
  iclag  = 0;
  idecal = 0;

  for (icha = 0; icha < *ncharb; icha++) {
    /* ---- Nb de classes */
    itypdp = _get_diameter_type(icha+1);
    nclpch[icha] = _cs_gui_get_nb_class(icha+1, itypdp);
    if (nclpch[icha] > *ncpcmx) {
      bft_error(__FILE__, __LINE__, 0,
                _("class number by coal is limited.\n"
                  "For coal %i it is %i \n in the parametric file \n"),
                  icha, nclpch[icha]);
    }

    /*---- Calcul du nb de classes et remplissage de ICHCOR */
    *nclacp = *nclacp + nclpch[icha];
    for (iclapc = 0; iclapc < nclpch[icha]; iclapc++) {
      icla = iclapc + idecal;
      ichcor[icla] = icha + 1;
    }
    idecal += nclpch[icha];

    /* Type de diametres  = 1 ---> diametre donnes
                          = 2 ---> loi de Rosin-Rammler */
    if (itypdp == 1) {
      for (icla = 0; icla < nclpch[icha]; icla++)
        diam20[icla + iclag] = _get_solid_fuel_diameter(icha+1,icla+1);
    } else if (itypdp == 2) {
      nbrf = _cs_gui_get_nb_refusal(icha+1);

      BFT_MALLOC(dprefus, nbrf,         double);
      BFT_MALLOC(refus,   nbrf,         double);
      BFT_MALLOC(pourc,   nclpch[icha], double);

      for (ii = 0; ii < nbrf; ii++) {
        dprefus[ii] = _get_solid_fuel_refusal_diameter(ii+1,icha+1)*1.e6;  //en microns
        refus[ii] = _get_solid_fuel_refusal_value(ii+1,icha+1);
      }
      for (ii = 0; ii<nclpch[icha]; ii++)
        pourc[ii] = _get_solid_fuel_pourc(ii+1,icha+1);

      /* decoupage des classes */
      BFT_MALLOC(rf, nclpch[icha], double);
      rf[0] = pourc[0] / 2.;

      for (icla = 1; icla < nclpch[icha]; icla++)
        rf[icla] = rf[icla-1] + (pourc[icla] + pourc[icla-1]) / 2.;

      kk1 = 0.;
      kk2 = 0.;
      kk3 = 0.;
      kk4 = 0.;
      for (ii = 0; ii < nbrf ; ii++) {
        kk1 = kk1 + log(dprefus[ii]);
        kk2 = kk2 + log(-log(refus[ii]));
        kk3 = kk3 + log(dprefus[ii])*log(dprefus[ii]);
        kk4 = kk4 + log(dprefus[ii])*log(-log(refus[ii]));
      }

      qq  = (nbrf * kk4 - kk1 * kk2) / (nbrf * kk3 - kk1 * kk1);
      var = (kk2 * kk3 - kk1 * kk4) / (nbrf * kk3 - kk1 * kk1);
      xx  = exp(-var / qq);

      for (icla = iclag; icla < iclag + nclpch[icha]; icla++)
        diam20[icla]=  xx*pow((-log(1.-rf[icla-iclag])),(1./qq))*1.e-6; // en metres

      bft_printf("** Rosin-Rammler results for the coal %i **\n"
                 "[ Checking of the Rosin-Rammler law ]\n"
                 "Diameter       refus given      refus computed\n\n", icha+1);

      for (icla = 0; icla< nbrf; icla++)
        bft_printf("%f     %f     %f \n", dprefus[icla], refus[icla],
            exp(-pow((dprefus[icla]/xx),(qq))));

      bft_printf("\nRefus       diam. given      diam. computed\n");

      for (icla = 0; icla< nbrf; icla++)
        bft_printf("%f     %f     %f \n", refus[icla], dprefus[icla],
            xx*pow((-log(refus[icla])),(1./qq)));

      bft_printf("\nDiameters computed by the Rosin-Rammler law\n");

      for (icla = iclag; icla <iclag+nclpch[icha]; icla ++)
        bft_printf("%d     %f \n", icla-iclag, diam20[icla]);

      BFT_FREE(pourc);
      BFT_FREE(refus);
      BFT_FREE(dprefus);
      BFT_FREE(rf);

    }
    else {
      bft_error(__FILE__, __LINE__, 0,
                _("type diameter value must be equal to 1 or 2.\n"
                  "Calculation is interupted \n"));
    }

    iclag = iclag + nclpch[icha];

    /* ---- Composition elementaire en C, H , O , N , S sur sec (% en masse) */
    cch[icha] = _get_solid_fuel_composition_on_dry(icha+1,"C_composition_on_dry");
    hch[icha] = _get_solid_fuel_composition_on_dry(icha+1,"H_composition_on_dry");
    och[icha] = _get_solid_fuel_composition_on_dry(icha+1,"O_composition_on_dry");
    nch[icha] = _get_solid_fuel_composition_on_dry(icha+1,"N_composition_on_dry");
    sch[icha] = _get_solid_fuel_composition_on_dry(icha+1,"S_composition_on_dry");


    /* ---- Composition elementaire en C, H , O , N , S du coke (% en masse) */
    cck[icha] = _get_solid_fuel_composition_on_dry(icha+1,"C_coke_composition_on_dry");
    hck[icha] = _get_solid_fuel_composition_on_dry(icha+1,"H_coke_composition_on_dry");
    ock[icha] = _get_solid_fuel_composition_on_dry(icha+1,"O_coke_composition_on_dry");
    nck[icha] = _get_solid_fuel_composition_on_dry(icha+1,"N_coke_composition_on_dry");
    sck[icha] = _get_solid_fuel_composition_on_dry(icha+1,"S_coke_composition_on_dry");

    /* ---- PCI sur charbon sec ou pur suivant la valeur de IPCI */
    ipci[icha] = _get_PCI_type(icha+1);
    if (ipci[icha] < 6)
      pcich[icha] = _get_PCI_value(icha+1);

    volatile_matter = _get_volatile_matter(icha+1);
    h0ashc[icha] = _get_ashes_enthalpy(icha+1);
    cpashc[icha] = _get_ashes_thermal_capacity(icha+1);

    /*  ---- CP moyen du charbon sec (J/kg/K) */
    cp2ch[icha] = _get_solid_fuel_specific_heat_average(icha+1);

    /* ---- Masse volumique initiale (kg/m3) */
    rho0ch[icha] = _get_solid_fuel_density(icha+1);

    /* ---- Thermal conductivity of the coal (W/m/K) */
    if (vars != NULL && thcdch != NULL)
      if (cs_gui_strcmp(vars->model_value, "homogeneous_fuel_moisture"))
        thcdch[icha] = _get_solid_fuel_thermal_conductivity(icha+1);

    /* ---- Caracteristiques cendres */

    /* ------ Taux de cendre (kg/kg) en % */
    xashsec[icha] = _get_ashes_rate(icha+1);

    /*      Transformation en kg/kg */
    xashch[icha] = xashsec[icha]/100.;

    /* ------ Taux d'humidite (kg/kg) en % */
    xwatch[icha] = _get_humidity_rate(icha+1);

    /*      Transformation en kg/kg */
    xwatch[icha] = xwatch[icha]/100.;

    /*      Transformation du taux de cendre de sec
            sur humide en kg/kg */
    xashch[icha] = xashch[icha]*(1.-xwatch[icha]);

    /* ---- Parametres de devolatilisation (modele de Kobayashi) */

    iy1ch[icha] = _get_Y1Y2_coefficient_type(icha+1);
    iy2ch[icha] = iy1ch[icha];
    if (iy1ch[icha] > 0) {
      y1ch[icha] = _get_Y1_coefficient(icha+1);
      y2ch[icha] = _get_Y2_coefficient(icha+1);
    }
    a1ch[icha] = _get_A1_pre_exponential_factor(icha+1);
    a2ch[icha] = _get_A2_pre_exponential_factor(icha+1);
    e1ch[icha] = _get_E1_energy_of_activation(icha+1);
    e2ch[icha] = _get_E2_energy_of_activation(icha+1);

    /*  ---- Parametres combustion heterogene pour O2
        (modele a sphere retrecissante) */
    ahetch[icha] = _get_pre_exponential_constant(icha+1,"O2");
    ehetch[icha] = _get_energy_of_activation(icha+1,"O2");
    iochet[icha] = _get_order_of_reaction(icha+1,"O2");

    /* ---- Parametres combustion heterogene pour CO2
       (modele a sphere retrecissante) */
    _getCO2KineticsStatus(ihtco2);
    if (*ihtco2) {
      ahetc2[icha] = _get_pre_exponential_constant(icha+1,"CO2");
      ehetc2[icha] = _get_energy_of_activation(icha+1,"CO2");
      ioetc2[icha] = _get_order_of_reaction(icha+1,"CO2");
    }

    /* ---- Parametres combustion heterogene pour H2O
       (modele a sphere retrecissante) */
    _getH2OKineticsStatus(ihth2o);
    if (*ihth2o) {
      ahetwt[icha] = _get_pre_exponential_constant(icha+1,"H2O");
      ehetwt[icha] = _get_energy_of_activation(icha+1,"H2O");
      ioetwt[icha] = _get_order_of_reaction(icha+1,"H2O");
    }

    /* ---- Parametres modele de NOX
       QPR =  %d'azote libere pendant la devol
       / %de MV libere pendant la devol */
    _getNOxStatus(ieqnox);
    if (*ieqnox) {
      qpr[icha] = _get_nitrogen_fraction(icha+1);
      fn[icha] = _get_nitrogen_concentration(icha+1);

      /* ---- Repartition de l'azote entre HCN et NH3 */
      crepn1[icha] = _get_nitrogen_partition_in_HCN_NH3_reaction_1(icha+1);
      crepn1[*ncharb+icha] = 1-crepn1[icha];
      crepn2[icha] = _get_nitrogen_partition_in_HCN_NH3_reaction_2(icha+1);
      crepn2[*ncharb+icha] = 1-crepn2[icha];

      repnck[icha] = _get_hcn_char_comb(icha+1);
      repnle[icha] = _get_nitrogen_char_low_temp(icha+1);
      repnlo[icha] = _get_nitrogen_char_high_temp(icha+1);
      _getNOxFeatureStatus(icha+1, imdnox);
      if (*imdnox)
        _getNOxReburning(icha+1, irb);
    }
    else {
          crepn1[icha] = 0.5;
      crepn1[*ncharb+icha] = 1-crepn1[icha];
      crepn2[icha] = 0.5;
      crepn2[*ncharb+icha] = 1-crepn2[icha];
    }
  }

  /*        ! --> Lecture rayonnement */
  /*  ! ---- Coefficient d'absorption du melange gazeux */
  if (*iirayo>0)
    *ckabs1 = _get_absorption_coefficient();

  /* --> Lecture caracteristiques Oxydants */

  /* ---- Nb d'oxydants */
  *noxyd = cs_gui_get_tag_count("/oxidants/oxidant", 1);
  if (*noxyd < 1 || *noxyd > 3 ) {
    bft_error(__FILE__, __LINE__, 0,
        _("Oxidant number must be between 1 and 3.\n"
          "It is  %i in the parametric file \n"
          "Calculation is interupted \n"),
        *noxyd);
  }
  itypoxy = _get_oxydant_type();

  /* ---- Composition en O2,N2,H2O,N2 */

  for (ioxy = 0; ioxy < 3; ioxy++) {
    oxyo2 [ioxy] = 0.;
    oxyn2 [ioxy] = 0.;
    oxyh2o[ioxy] = 0.;
    oxyco2[ioxy] = 0.;
  }

  for (ioxy = 0; ioxy < *noxyd; ioxy++) {
    oxyo2[ioxy] = _get_oxydant_composition(ioxy+1,"O2_composition");
    oxyn2[ioxy] = _get_oxydant_composition(ioxy+1,"N2_composition");
    oxyh2o[ioxy] = _get_oxydant_composition(ioxy+1,"H2O_composition");
    oxyco2[ioxy] = _get_oxydant_composition(ioxy+1,"CO2_composition");
  }

  if (itypoxy == 1) {
    /* transformation pourcentage volumique en nombre de mole */
    for (ioxy = 0; ioxy<*noxyd ; ioxy++) {
      coef = 100.;
      if (oxyo2[ioxy] > 0.)
        coef = CS_MIN(coef,oxyo2[ioxy]);
      if (oxyn2[ioxy] > 0.)
        coef = CS_MIN(coef,oxyn2[ioxy]);
      if (oxyh2o[ioxy] > 0.)
        coef = CS_MIN(coef,oxyh2o[ioxy]);
      if (oxyco2[ioxy] > 0.)
        coef = CS_MIN(coef,oxyco2[ioxy]);

      oxyo2 [ioxy] = oxyo2 [ioxy]/coef;
      oxyn2 [ioxy] = oxyn2 [ioxy]/coef;
      oxyh2o[ioxy] = oxyh2o[ioxy]/coef;
      oxyco2[ioxy] = oxyco2[ioxy]/coef;
    }
  }

  BFT_FREE(model);
}

/*----------------------------------------------------------------------------
 * Copy name of thermophysical data file from C to Fortran
 *----------------------------------------------------------------------------*/

void CS_PROCF(cfnmtd, CFNMTD) (char          *fstr,    /* --> Fortran string */
                               int           *len      /* --> String Length  */
                               CS_ARGF_SUPP_CHAINE)
{
  int i;
  int l = 0;
  char *cstr = NULL;
  char *path = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3,
                        "thermophysical_models",
                        "gas_combustion",
                        "data_file");

  cs_xpath_add_function_text(&path);
  cstr = cs_gui_get_text_value(path);

  BFT_FREE(path);

  /* Copy string */

  if (cstr != NULL) {

    /* Compute string length (removing start or end blanks) */

    l = strlen(cstr);
    if (l > *len)
      l = *len;

    for (i = 0; i < l; i++)
      fstr[i] = cstr[i];

    /* Pad with blanks if necessary */

    for (i = l; i < *len; i++)
      fstr[i] = ' ';

    BFT_FREE(cstr);

  }

}

/*----------------------------------------------------------------------------
 * groundwater model : read parameters
 *
 * Fortran Interface:
 *
 * subroutine uidai1
 * *****************
 * integer         permeability    <--   permeability type
 * integer         dispersion      <--   dispersion type
 * integer         unsteady        <--   steady flow
 * integer         gravity         <--   check if gravity is taken into account
 * integer         unsaturated     <--   take into account unsaturated zone
 *----------------------------------------------------------------------------*/

void CS_PROCF (uidai1, UIDAI1) (int    *const permeability,
                                int    *const dispersion,
                                int    *const unsteady,
                                int    *const gravity,
                                int    *const unsaturated)
{
  char *path   = NULL;
  char *mdl    = NULL;
  int   result;

  /* Get dispersion type */
  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "thermophysical_models",
                                  "groundwater_model",
                                  "dispersion");
  cs_xpath_add_attribute(&path, "model");
  mdl = cs_gui_get_attribute_value(path);
  BFT_FREE(path);

  if (cs_gui_strcmp(mdl, "anisotropic"))
    *dispersion = 1;
  else
    *dispersion = 0;
  BFT_FREE(mdl);

  /* Get flow type */
  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "thermophysical_models",
                                  "groundwater_model",
                                  "flowType");
  cs_xpath_add_attribute(&path, "model");
  mdl = cs_gui_get_attribute_value(path);
  BFT_FREE(path);

  if (cs_gui_strcmp(mdl, "steady"))
    *unsteady = 0;
  else
    *unsteady = 1;
  BFT_FREE(mdl);

  /* Get permeability type */
  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "thermophysical_models",
                                  "groundwater_model",
                                  "permeability");
  cs_xpath_add_attribute(&path, "model");
  mdl = cs_gui_get_attribute_value(path);
  BFT_FREE(path);

  if (cs_gui_strcmp(mdl, "anisotropic"))
    *permeability = 1;
  else
    *permeability = 0;
  BFT_FREE(mdl);

  /* Get gravity */
  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "thermophysical_models",
                                  "groundwater_model",
                                  "gravity");
  cs_xpath_add_attribute(&path, "status");

  if (cs_gui_get_status(path, &result))
    *gravity = result;
  BFT_FREE(path);

  /* Get the possible presence of unsaturated zone */
  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "thermophysical_models",
                                  "groundwater_model",
                                  "unsaturatedZone");
  cs_xpath_add_attribute(&path, "model");
  mdl = cs_gui_get_attribute_value(path);
  BFT_FREE(path);

  if (cs_gui_strcmp(mdl, "true"))
    *unsaturated = 1;
  else
    *unsaturated = 0;
  BFT_FREE(mdl);

  /* Get first-order decay rate and chemistry model */
  int n_fields = cs_field_n_fields();
  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (   (f->type & CS_FIELD_VARIABLE)
           && (f->type & CS_FIELD_USER)) {

      /* get first-order decay rate */
      cs_real_t decay = 0.;
      path = cs_xpath_init_path();
      cs_xpath_add_elements(&path, 2, "thermophysical_models",
                                      "groundwater_model");
      cs_xpath_add_element(&path, "scalar");
      cs_xpath_add_test_attribute(&path, "name", f->name);
      cs_xpath_add_element(&path, "fo_decay_rate");
      cs_xpath_add_function_text(&path);
      cs_gui_get_double(path, &decay);
      BFT_FREE(path);

      int key_decay = cs_field_key_id("fo_decay_rate");
      cs_field_set_key_double(f, key_decay, decay);

      /* get chemistry model */
      path = cs_xpath_init_path();
      cs_xpath_add_elements(&path, 2, "thermophysical_models",
                                      "groundwater_model");
      cs_xpath_add_element(&path, "scalar");
      cs_xpath_add_test_attribute(&path, "name", f->name);
      cs_xpath_add_attribute(&path, "chemistry_model");
      mdl = cs_gui_get_attribute_value(path);
      BFT_FREE(path);

      cs_gwf_soilwater_partition_t sorption_scal;
      int key_part = cs_field_key_id("gwf_soilwater_partition");
      cs_field_get_key_struct(f, key_part, &sorption_scal);
      if (cs_gui_strcmp(mdl, "EK"))
        sorption_scal.kinetic = 1;
      else
        sorption_scal.kinetic = 0;
      cs_field_set_key_struct(f, key_part, &sorption_scal);
      BFT_FREE(mdl);
    }
  }

#if _XML_DEBUG_
  bft_printf("==>UIDAI1\n");
  bft_printf("--groundwater_anisotropic_permeability  = %i\n", *permeability);
  bft_printf("--groundwater_anisotropic_dispersion    = %f\n", *dispersion);
  bft_printf("--groundwater_unsteady                  = %f\n", *unsteady);
  bft_printf("--groundwater_gravity                   = %f\n", *gravity);
  bft_printf("--groundwater_unsaturated               = %f\n", *unsaturated);
#endif
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Activate specific physical models based on XML settings.
 *
 * parameters:
 *   ieos    --> compressible
 *   ieqco2  --> CO2 massic fraction transport (for combustion only)
 *----------------------------------------------------------------------------*/

void
cs_gui_physical_model_select(cs_int_t  *ieos,
                             cs_int_t  *ieqco2)
{
  if (!cs_gui_file_is_loaded())
    return;

  int isactiv = 0;

  cs_var_t  *vars = cs_glob_var;

  *ieqco2 = 0;

  /* Look for the active specific physics and give the value of the associated
     model attribute */
  isactiv = cs_gui_get_activ_thermophysical_model();

  if (isactiv)  {

    if (cs_gui_strcmp(vars->model, "solid_fuels")) {
      if (cs_gui_strcmp(vars->model_value, "homogeneous_fuel"))
        cs_glob_physical_model_flag[CS_COMBUSTION_COAL] = 0;
      else if (cs_gui_strcmp(vars->model_value,
                             "homogeneous_fuel_moisture"))
        cs_glob_physical_model_flag[CS_COMBUSTION_COAL] = 1;
      else
        bft_error(__FILE__, __LINE__, 0,
                  _("Invalid coal model: %s.\n"), vars->model_value);
      /* ieqco2 fix to transport of CO2 mass fraction */
      *ieqco2 = 1;
    }
    else if  (cs_gui_strcmp(vars->model, "gas_combustion")) {
      char *path = NULL;
      path = cs_xpath_init_path();
      cs_xpath_add_elements(&path, 2, "thermophysical_models", "gas_combustion");
      cs_xpath_add_attribute(&path, "model");

      char *model = cs_gui_get_attribute_value(path);

      BFT_FREE(path);

      if (!cs_gui_strcmp(model, "off")) {

        if (cs_gui_strcmp(vars->model_value, "adiabatic"))
          cs_glob_physical_model_flag[CS_COMBUSTION_3PT] = 0;
        else if (cs_gui_strcmp(vars->model_value, "extended"))
          cs_glob_physical_model_flag[CS_COMBUSTION_3PT] = 1;
        else if (cs_gui_strcmp(vars->model_value, "spalding"))
          cs_glob_physical_model_flag[CS_COMBUSTION_EBU] = 0;
        else if (cs_gui_strcmp(vars->model_value, "enthalpy_st"))
          cs_glob_physical_model_flag[CS_COMBUSTION_EBU] = 1;
        else if (cs_gui_strcmp(vars->model_value, "mixture_st"))
          cs_glob_physical_model_flag[CS_COMBUSTION_EBU] = 2;
        else if (cs_gui_strcmp(vars->model_value, "enthalpy_mixture_st"))
          cs_glob_physical_model_flag[CS_COMBUSTION_EBU] = 3;
        else if (cs_gui_strcmp(vars->model_value, "2-peak_adiabatic"))
          cs_glob_physical_model_flag[CS_COMBUSTION_LW] = 0;
        else if (cs_gui_strcmp(vars->model_value, "2-peak_enthalpy"))
          cs_glob_physical_model_flag[CS_COMBUSTION_LW] = 1;
        else if (cs_gui_strcmp(vars->model_value, "3-peak_adiabatic"))
          cs_glob_physical_model_flag[CS_COMBUSTION_LW] = 2;
        else if (cs_gui_strcmp(vars->model_value, "3-peak_enthalpy"))
          cs_glob_physical_model_flag[CS_COMBUSTION_LW] = 3;
        else if (cs_gui_strcmp(vars->model_value, "4-peak_adiabatic"))
          cs_glob_physical_model_flag[CS_COMBUSTION_LW] = 4;
        else if (cs_gui_strcmp(vars->model_value, "4-peak_enthalpy"))
          cs_glob_physical_model_flag[CS_COMBUSTION_LW] = 5;
        else
          bft_error(__FILE__, __LINE__, 0,
                    _("Invalid gas combustion flow model: %s.\n"),
                    vars->model_value);
      }
    }
    else if (cs_gui_strcmp(vars->model, "atmospheric_flows")) {
      if (cs_gui_strcmp(vars->model_value, "constant"))
        cs_glob_physical_model_flag[CS_ATMOSPHERIC] = 0;
      else if (cs_gui_strcmp(vars->model_value, "dry"))
        cs_glob_physical_model_flag[CS_ATMOSPHERIC] = 1;
      else if (cs_gui_strcmp(vars->model_value, "humid"))
        cs_glob_physical_model_flag[CS_ATMOSPHERIC] = 2;
      else
        bft_error(__FILE__, __LINE__, 0,
            _("Invalid atmospheric flow model: %s.\n"),
            vars->model_value);
    }
    else if (cs_gui_strcmp(vars->model, "joule_effect")) {
      if (cs_gui_strcmp(vars->model_value, "joule")) {
        char *value = _get_joule_model();
        if (cs_gui_strcmp(value, "AC/DC"))
          cs_glob_physical_model_flag[CS_JOULE_EFFECT] = 1;
        else if (cs_gui_strcmp(value, "three-phase"))
          cs_glob_physical_model_flag[CS_JOULE_EFFECT] = 2;
        else if (cs_gui_strcmp(value, "AC/DC+Transformer"))
          cs_glob_physical_model_flag[CS_JOULE_EFFECT] = 3;
        else if (cs_gui_strcmp(value, "three-phase+Transformer"))
          cs_glob_physical_model_flag[CS_JOULE_EFFECT] = 4;
        else
          bft_error(__FILE__, __LINE__, 0,
              _("Invalid joule model: %s.\n"),
              vars->model_value);
        BFT_FREE(value);
      }
      else if (cs_gui_strcmp(vars->model_value, "arc"))
        cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] = 2;
      else
        bft_error(__FILE__, __LINE__, 0,
            _("Invalid electrical model: %s.\n"),
            vars->model_value);
    }
    else if (cs_gui_strcmp(vars->model, "compressible_model")) {
      if (cs_gui_strcmp(vars->model_value, "constant_gamma")) {
        cs_glob_physical_model_flag[CS_COMPRESSIBLE] = 0;
        *ieos = 1;
      }
      else
        bft_error(__FILE__, __LINE__, 0,
            _("Invalid compressible model: %s.\n"),
            vars->model_value);
    }
    else if (cs_gui_strcmp(vars->model, "groundwater_model")) {
      if (cs_gui_strcmp(vars->model_value, "groundwater"))
        cs_glob_physical_model_flag[CS_GROUNDWATER] = 1;
    }
  }

#if _XML_DEBUG_
  bft_printf("%s\n", __func__);
  if (isactiv) {
    bft_printf("--thermophysical model: %s\n", vars->model);
    bft_printf("--thermophysical value: %s\n", vars->model_value);
  }
#endif
}

/*----------------------------------------------------------------------------
 * Electrical model: read parameters
 *----------------------------------------------------------------------------*/

void
cs_gui_elec_model(void)
{
  char *path   = NULL;
  int   status = 0;
  double   result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "thermophysical_models",
                                  "joule_effect",
                                  "variable_scaling");

  cs_elec_option_t *elec_opt = cs_get_glob_elec_option();

  cs_xpath_add_attribute(&path, "status");
  if (cs_gui_get_status(path, &status))
    elec_opt->ielcor = status;
  BFT_FREE(path);

  int ieljou = cs_glob_physical_model_flag[CS_JOULE_EFFECT];
  int ielarc = cs_glob_physical_model_flag[CS_ELECTRIC_ARCS];

  if (ieljou > 0) {
    path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 3, "thermophysical_models",
                                    "joule_effect",
                                    "imposed_power");

    cs_xpath_add_function_text(&path);
    if (!cs_gui_get_double(path, &result))
      bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);
    elec_opt->puisim = result;
    BFT_FREE(path);
  }

  if (ielarc > 0) {
    path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 3, "thermophysical_models",
                                    "joule_effect",
                                    "imposed_current");
    cs_xpath_add_function_text(&path);
    if (!cs_gui_get_double(path, &result))
      bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);
    elec_opt->couimp = result;
    BFT_FREE(path);

    if (cs_glob_elec_option->ielcor > 0) {
      path = cs_xpath_init_path();
      char *choice;

      cs_xpath_add_elements(&path, 3, "thermophysical_models",
                                      "joule_effect",
                                      "recal_model");
      cs_xpath_add_attribute(&path, "model");
      choice = cs_gui_get_attribute_value(path);
      if (cs_gui_strcmp(choice, "general_case"))
        elec_opt->modrec = 1;
      else if (cs_gui_strcmp(choice, "plane_define"))
        elec_opt->modrec = 2;
      else if (cs_gui_strcmp(choice, "user"))
        elec_opt->modrec = 3;
      else
        bft_error(__FILE__, __LINE__, 0, _("Invalid model : %s\n"), choice);

      BFT_FREE(choice);

      if (cs_glob_elec_option->modrec == 2) {
        path = cs_xpath_init_path();
        cs_xpath_add_elements(&path, 4, "thermophysical_models",
                                        "joule_effect",
                                        "recal_model",
                                        "direction");
        cs_xpath_add_function_text(&path);
        choice = cs_gui_get_text_value(path);
        if (cs_gui_strcmp(choice, "X"))
          elec_opt->idreca = 1;
        else if (cs_gui_strcmp(choice, "Y"))
          elec_opt->idreca = 2;
        else
          elec_opt->idreca = 3;
        BFT_FREE(path);
        BFT_FREE(choice);

        double val;
        path = cs_xpath_init_path();
        cs_xpath_add_elements(&path, 5, "thermophysical_models",
                                        "joule_effect",
                                        "recal_model",
                                        "plane_definition",
                                        "A");
        cs_xpath_add_function_text(&path);
        if (!cs_gui_get_double(path, &val))
          bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);
        BFT_FREE(path);
        elec_opt->crit_reca[0] = val;

        path = cs_xpath_init_path();
        cs_xpath_add_elements(&path, 5, "thermophysical_models",
                                        "joule_effect",
                                        "recal_model",
                                        "plane_definition",
                                        "B");
        cs_xpath_add_function_text(&path);
        if (!cs_gui_get_double(path, &val))
          bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);
        BFT_FREE(path);
        elec_opt->crit_reca[1] = val;

        path = cs_xpath_init_path();
        cs_xpath_add_elements(&path, 5, "thermophysical_models",
                                        "joule_effect",
                                        "recal_model",
                                        "plane_definition",
                                        "C");
        cs_xpath_add_function_text(&path);
        if (!cs_gui_get_double(path, &val))
          bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);
        BFT_FREE(path);
        elec_opt->crit_reca[2] = val;

        path = cs_xpath_init_path();
        cs_xpath_add_elements(&path, 5, "thermophysical_models",
                                        "joule_effect",
                                        "recal_model",
                                        "plane_definition",
                                        "D");
        cs_xpath_add_function_text(&path);
        if (!cs_gui_get_double(path, &val))
          bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);
        BFT_FREE(path);
        elec_opt->crit_reca[3] = val;

        path = cs_xpath_init_path();
        cs_xpath_add_elements(&path, 5, "thermophysical_models",
                                        "joule_effect",
                                        "recal_model",
                                        "plane_definition",
                                        "epsilon");
        cs_xpath_add_function_text(&path);
        if (!cs_gui_get_double(path, &val))
          bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);
        BFT_FREE(path);
        elec_opt->crit_reca[4] = val;
      }
    }
    BFT_FREE(path);
  }

#if _XML_DEBUG_
  bft_printf("==>UIELI1\n");
  bft_printf("--ielcor  = %i\n", cs_glob_elec_option->ielcor);
  bft_printf("--puisim  = %f\n", cs_glob_elec_option->puisim);
  bft_printf("--couimp  = %f\n", cs_glob_elec_option->couimp);
  bft_printf("--modrec  = %f\n", cs_glob_elec_option->modrec);
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
  /* build list of cells */
  char *crit = NULL;

  cs_lnum_t   n_selected_faces = 0;
  cs_lnum_t  *selected_faces = NULL;

  BFT_MALLOC(crit, 66, char);

  char cVal[10];

  cs_elec_option_t *elec_opt = cs_get_glob_elec_option();

  strcpy(crit, "plane[");
  sprintf(cVal, "%f", elec_opt->crit_reca[0]);
  strcat(crit, cVal);
  strcat(crit, ",");
  sprintf(cVal, "%f", elec_opt->crit_reca[1]);
  strcat(crit, cVal);
  strcat(crit, ",");
  sprintf(cVal, "%f", elec_opt->crit_reca[2]);
  strcat(crit, cVal);
  strcat(crit, ",");
  sprintf(cVal, "%f", elec_opt->crit_reca[3]);
  strcat(crit, cVal);
  strcat(crit, ",epsilon=");
  sprintf(cVal, "%6f", elec_opt->crit_reca[4]);
  strcat(crit, cVal);
  strcat(crit, "]");

  BFT_MALLOC(selected_faces, cs_glob_mesh->n_i_faces, cs_lnum_t);

  cs_selector_get_i_face_list(crit,
                              &n_selected_faces,
                              selected_faces);

  for (int j=0; j < n_selected_faces; j++)
    elec_opt->izreca[selected_faces[j]] = 1;

  BFT_FREE(selected_faces);
  BFT_FREE(crit);

}

/*-----------------------------------------------------------------------------
 * Return the name of a thermophysical model.
 *
 * parameter:
 *   model_thermo          -->  thermophysical model
 *----------------------------------------------------------------------------*/

char*
cs_gui_get_thermophysical_model(const char *const model_thermo)
{
  char *model = NULL;
  char *path = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "thermophysical_models", model_thermo);
  if (cs_gui_strcmp(model_thermo, "gas_combustion"))
    cs_xpath_add_attribute(&path, "option");
  else
    cs_xpath_add_attribute(&path, "model");

  model = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return model;
}

/*-----------------------------------------------------------------------------
 * Return 1 if a specific physics model is activated. Store in the global
 * structure vars:
 *   vars->model         <= thermophysical model
 *   vars->model_value   <= related model name
 *----------------------------------------------------------------------------*/

int
cs_gui_get_activ_thermophysical_model(void)
{
  int i, isactiv = 0;

  if (cs_glob_var == NULL)
    cs_gui_init();

  cs_var_t  *vars = cs_glob_var;

  const char *name[] = { "solid_fuels",
                         "gas_combustion",
                         "joule_effect",
                         "atmospheric_flows",
                         "compressible_model",
                         "groundwater_model" };
  int name_nbr = sizeof(name) / sizeof(name[0]);

  if (vars->model != NULL && vars->model_value != NULL) {
    isactiv = 1;
    return isactiv;
  } else {
    BFT_FREE(vars->model);
    vars->model_value = NULL;
  }

  for (i = 0; i < name_nbr; i++) {
    char *value = cs_gui_get_thermophysical_model(name[i]);

    if (value && !cs_gui_strcmp(value, "off")) {
      BFT_MALLOC(vars->model, strlen(name[i])+1, char);
      strcpy(vars->model, name[i]);

      BFT_MALLOC(vars->model_value, strlen(value)+1, char);
      strcpy(vars->model_value, value);

      isactiv = 1;
      BFT_FREE(value);
      break;
    }

    BFT_FREE(value);
  }

  return isactiv;
}

/*------------------------------------------------------------------------------
 * Set GUI-defined labels for the atmospheric module
 *----------------------------------------------------------------------------*/

void
cs_gui_labels_atmospheric(void)
{
  cs_field_t *f = NULL;

  f = CS_F_(pot_t); /* Thermal scalar should be potential Temperature */
  if (f != NULL)
    _set_thermal_scalar_name_label(f, "label");

  f = CS_F_(totwt);
  if (f != NULL)
    _set_scalar_name_label(f, "atmospheric_flows", "total_water");

  f = CS_F_(ntdrp);
  if (f != NULL)
    _set_scalar_name_label(f, "atmospheric_flows", "number_of_droplets");
}

/*------------------------------------------------------------------------------
 * Set GUI-defined labels for the coal combustion module
 *
 * parameters:
 *   n_coals   <-- number of coals
 *   n_classes <-- number of coal classes
 *----------------------------------------------------------------------------*/

void
cs_gui_labels_coal_combustion(int  n_coals,
                              int  n_classes)
{
  int i;
  char name[64];
  cs_field_t *f = NULL;

  f = CS_F_(h); /* Thermal scalar should be Enthalpy */
  if (f != NULL)
    _set_thermal_scalar_name_label(f, "label");

  for (i = 0; i < n_classes; i++) {
    f = CS_FI_(h2, i);
    if (f != NULL) {
      snprintf(name, 63, "%s%2.2i", "x_p_h_", i+1); name[63] = '\0';
      _set_scalar_name_label(f, "solid_fuels", name);
    }
  }

  for (i = 0; i < n_classes; i++) {
    f = CS_FI_(np, i);
    if (f != NULL) {
      snprintf(name, 63, "%s%2.2i", "n_p_", i+1); name[63] = '\0';
      _set_scalar_name_label(f, "solid_fuels", name);
    }
  }

  for (i = 0; i < n_classes; i++) {
    f = CS_FI_(xch, i);
    if (f != NULL) {
      snprintf(name, 63, "%s%2.2i", "x_p_coal_", i+1); name[63] = '\0';
      _set_scalar_name_label(f, "solid_fuels", name);
    }
  }

  for (i = 0; i < n_classes; i++) {
    f = CS_FI_(xck, i);
    if (f != NULL) {
      snprintf(name, 63, "%s%2.2i", "x_p_char_", i+1); name[63] = '\0';
      _set_scalar_name_label(f, "solid_fuels", name);
    }
  }

  for (i = 0; i < n_classes; i++) {
    f = CS_FI_(xwt, i);
    if (f != NULL) {
      snprintf(name, 63, "%s%2.2i", "x_p_wt_", i+1); name[63] = '\0';
      _set_scalar_name_label(f, "solid_fuels", name);
    }
  }

  for (i = 0; i < n_coals; i++) {
    f = CS_FI_(f1m, i);
    if (f != NULL) {
      snprintf(name, 63, "%s%2.2i", "fr_mv1_", i+1); name[63] = '\0';
      _set_scalar_name_label(f, "solid_fuels", name);
    }
  }

  for (i = 0; i < n_coals; i++) {
    f = CS_FI_(f2m, i);
    if (f != NULL) {
      snprintf(name, 63, "%s%2.2i", "fr_mv2_", i+1); name[63] = '\0';
      _set_scalar_name_label(f, "solid_fuels", name);
    }
  }

  f = CS_F_(f4m);
  if (f != NULL)
    _set_scalar_name_label(f, "solid_fuels", "fr_oxyd2");

  f = CS_F_(f5m);
  if (f != NULL)
    _set_scalar_name_label(f, "solid_fuels", "fr_oxyd3");

  f = CS_F_(f6m);
  if (f != NULL)
    _set_scalar_name_label(f, "solid_fuels", "fr_h2o");

  f = CS_F_(f7m);
  if (f != NULL)
    _set_scalar_name_label(f, "solid_fuels", "fr_het_o2");

  f = CS_F_(f8m);
  if (f != NULL)
    _set_scalar_name_label(f, "solid_fuels", "fr_het_co2");

  f = CS_F_(f9m);
  if (f != NULL)
    _set_scalar_name_label(f, "solid_fuels", "fr_het_h2o");

  f = CS_F_(fvp2m);
  if (f != NULL)
    _set_scalar_name_label(f, "solid_fuels", "f1f2_variance");

  f = CS_F_(yco2);
  if (f != NULL)
    _set_scalar_name_label(f, "solid_fuels", "x_c_co2");

  f = CS_F_(yhcn);
  if (f != NULL)
    _set_scalar_name_label(f, "solid_fuels", "x_c_hcn");

  f = CS_F_(yno);
  if (f != NULL)
    _set_scalar_name_label(f, "solid_fuels", "x_c_no");

  f = CS_F_(ynh3);
  if (f != NULL)
    _set_scalar_name_label(f, "solid_fuels", "x_c_nh3");

  f = CS_F_(hox);
  if (f != NULL)
    _set_scalar_name_label(f, "solid_fuels", "x_c_h_ox");
}

/*------------------------------------------------------------------------------
 * Set GUI-defined labels for the compressible model variables
 *----------------------------------------------------------------------------*/

void
cs_gui_labels_compressible(void)
{
  cs_field_t *f = NULL;

  f = CS_F_(energy);
  if (f != NULL)
    _set_thermal_scalar_name_label(f, "label");

  f = CS_F_(t_kelvin);
  if (f != NULL)
    _set_scalar_name_label(f, "compressible_model", "temperature");

}

/*------------------------------------------------------------------------------
 * Set GUI-defined labels for the electric arcs variables
 *
 * parameters:
 *   n_gasses <-- number of constituent gasses
 *----------------------------------------------------------------------------*/

void
cs_gui_labels_electric_arcs(int  n_gasses)
{
  int i;
  char name[64];
  cs_field_t *f = NULL;

  f = CS_F_(h); /* Thermal scalar should be Enthalpy */
  if (f != NULL)
    _set_thermal_scalar_name_label(f, "label");

  f = CS_F_(potr);
  if (f != NULL)
    _set_scalar_name_label(f, "joule_effect", "elec_pot_r");

  f = CS_F_(poti);
  if (f != NULL)
    _set_scalar_name_label(f, "joule_effect", "elec_pot_i");

  f = CS_F_(potva);
  if (f != NULL)
    _set_scalar_name_label(f, "joule_effect", "vec_potential");

  for (i = 0; i < n_gasses - 1; i++) {
    f = CS_FI_(ycoel, i);
    if (f != NULL) {
      snprintf(name, 63, "%s%2.2i", "esl_fraction_", i+1); name[63] = '\0';
      _set_scalar_name_label(f, "joule_effect", name);
    }
  }
}

/*------------------------------------------------------------------------------
 * Set GUI-defined labels for the gas combustion variables
 *----------------------------------------------------------------------------*/

void
cs_gui_labels_gas_combustion(void)
{
  cs_field_t *f = NULL;

  f = CS_F_(h); /* Thermal scalar should be Enthalpy */
  if (f != NULL)
    _set_thermal_scalar_name_label(f, "label");

  f = CS_F_(fm);
  if (f != NULL)
    _set_scalar_name_label(f, "gas_combustion", "mixture_fraction");

  f = CS_F_(fp2m);
  if (f != NULL)
    _set_scalar_name_label(f, "gas_combustion", "mixture_fraction_variance");

  /* TODO: handle labels for "fsm" and "npm" key pointers */

  f = CS_F_(ygfm);
  if (f != NULL)
    _set_scalar_name_label(f, "gas_combustion", "fresh_gas_fraction");

  f = CS_F_(yfm);
  if (f != NULL)
    _set_scalar_name_label(f, "gas_combustion", "mass_fraction");

  f = CS_F_(yfp2m);
  if (f != NULL)
    _set_scalar_name_label(f, "gas_combustion", "mass_fraction_variance");

  f = CS_F_(coyfp);
  if (f != NULL)
    _set_scalar_name_label(f, "gas_combustion", "mass_fraction_covariance");
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
