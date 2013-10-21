/*============================================================================
 * Management of the GUI parameters file: specific physics
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2013 EDF S.A.

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
#include "cs_prototypes.h"
#include "cs_selector.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gui_specific_physics.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* debugging switch */
#define _XML_DEBUG_ 0

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Return the activated specific physics scalar number
 *----------------------------------------------------------------------------*/

static int
_scalar_number(const char* model)
{
  char *path = NULL;
  int   nb;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "thermophysical_models");
  cs_xpath_add_element(&path, model);
  cs_xpath_add_element(&path, "scalar");

  nb = cs_gui_get_nb_element(path);

  BFT_FREE(path);

  return nb;
}

/*-----------------------------------------------------------------------------
 * Return the label or the name from a specific physics scalar.
 *
 * parameters:
 *   physics              -->  keyword: specific physic model required
 *   kw                   -->  keyword: name of scalar
 *----------------------------------------------------------------------------*/

static char *_scalar_name_label(const char *physics, const char *kw)
{
  char *path = NULL;
  char *str  = NULL;

  path = cs_xpath_short_path();
  cs_xpath_add_elements(&path, 3,
                        "thermophysical_models",
                        physics,
                        "scalar");
  cs_xpath_add_test_attribute(&path, "name", kw);
  cs_xpath_add_attribute(&path, "label");

  str = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return str;
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
  int   ichoice;
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

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

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

static double
_get_PCI_type(const int icha)
{
  char *path;
  char *path2;
  char *buff = NULL;
  char *buff2 = NULL;
  double ichoice;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_element(&path, "Heating_model");
  cs_xpath_add_attribute(&path,"choice");
  buff = cs_gui_get_attribute_value(path);
  if (buff == NULL) {
    ichoice = 0.;
  } else {
    if (cs_gui_strcmp(buff, "LHV"))
      {
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
        }
      }
    else if (cs_gui_strcmp(buff, "HHV"))
      {
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
        }
      }
    else if (cs_gui_strcmp(buff, "IGT_correlation"))
      ichoice = 6;
    else
      bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);
  }
  BFT_FREE(path);
  BFT_FREE(path2);
  BFT_FREE(buff);
  BFT_FREE(buff2);
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
  int ichoice;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,"thermophysical_models", "solid_fuels");
  cs_xpath_add_element_num(&path, "solid_fuel", icha);
  cs_xpath_add_elements(&path, 2,"devolatilisation_parameters",
                        "stoichiometric_coefficient");
  cs_xpath_add_attribute(&path,"type");
  buff = cs_gui_get_attribute_value(path);
  if (buff == NULL) {
    ichoice = 0;
  } else {
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
  int   ichoice;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path,3,
                        "thermophysical_models",
                        "solid_fuels",
                        "oxidants");
  cs_xpath_add_element(&path,"oxidant_type");
  cs_xpath_add_function_text(&path);
  buff = cs_gui_get_text_value(path);
  if (buff == NULL) {
    ichoice = 0;
  } else {
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

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Predefined physics indicator.
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIPPMO
 * *****************
 *
 * INTEGER          IPPMOD <--  specific physics indicator array
 * INTEGER          ICOD3P  --> diffusion flame in fast complete chemistry
 * INTEGER          ICODEQ  --> diffusion flame in fast chemistry towards balance
 * INTEGER          ICOEBU  --> Eddy Break Up premixing flame
 * INTEGER          ICOBML  --> Bray - Moss - Libby premixing flame
 * INTEGER          ICOLWC  --> Libby Williams premixing flame
 * INTEGER          ICP3PL  --> Coal combustion. Combustible moyen local
 * INTEGER          ICPL3C  --> Coal combustion coupled with lagrangien approach
 * INTEGER          ICFUEL  --> Fuel combustion
 * INTEGER          IELJOU  --> Joule effect
 * INTEGER          IELARC  --> electrical arc
 * INTEGER          IELION  --> ionique mobility
 * INTEGER          ICOMPF  --> compressible without shock
 * INTEGER          IATMOS  --> atmospheric flows
 * INTEGER          IAEROS  --> cooling tower
 * INTEGER          INDJON  --> INDJON=1: a JANAF enthalpy-temperature
 *                              tabulation is used. INDJON=1: users tabulation
 * INTEGER          IEOS    --> compressible
 * INTEGER          IEQCO2  --> CO2 massic fraction transport
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uippmo, UIPPMO)(int *const ippmod,
                               int *const icod3p,
                               int *const icodeq,
                               int *const icoebu,
                               int *const icobml,
                               int *const icolwc,
                               int *const iccoal,
                               int *const icpl3c,
                               int *const icfuel,
                               int *const ieljou,
                               int *const ielarc,
                               int *const ielion,
                               int *const icompf,
                               int *const iatmos,
                               int *const iaeros,
                               int *const ieos,
                               int *const ieqco2)
{
  int isactiv = 0;
  int nscapp = 0;

  cs_var_t  *vars = cs_glob_var;

  ippmod[*icod3p - 1] = -1;
  ippmod[*icodeq - 1] = -1;
  ippmod[*icoebu - 1] = -1;
  ippmod[*icobml - 1] = -1;
  ippmod[*icolwc - 1] = -1;
  ippmod[*iccoal - 1] = -1;
  ippmod[*icpl3c - 1] = -1;
  ippmod[*icfuel - 1] = -1;
  ippmod[*ieljou - 1] = -1;
  ippmod[*ielarc - 1] = -1;
  ippmod[*ielion - 1] = -1;
  ippmod[*icompf - 1] = -1;
  ippmod[*iatmos - 1] = -1;
  ippmod[*iaeros - 1] = -1;

  *ieqco2 = 0;

  /* Look for the active specific physics and give the value of the associated
     model attribute */
  isactiv = cs_gui_get_activ_thermophysical_model();

  if (isactiv)
  {
    if (cs_gui_strcmp(vars->model, "solid_fuels"))
    {
      if (cs_gui_strcmp(vars->model_value, "homogeneous_fuel"))
        ippmod[*iccoal - 1] = 0;
      else if (cs_gui_strcmp(vars->model_value, "homogeneous_fuel_moisture") || cs_gui_strcmp(vars->model_value, "homogeneous_fuel_moisture_lagr"))
        ippmod[*iccoal - 1] = 1;
      else
        bft_error(__FILE__, __LINE__, 0,
                  _("Invalid coal model: %s.\n"), vars->model_value);
      // ieqco2 fix to transport of CO2 mass fraction
      *ieqco2 = 1;
    }
    else if  (cs_gui_strcmp(vars->model, "gas_combustion"))
    {
      if (cs_gui_strcmp(vars->model_value, "adiabatic"))
        ippmod[*icod3p - 1] = 0;
      else if (cs_gui_strcmp(vars->model_value, "extended"))
        ippmod[*icod3p - 1] = 1;
      else if (cs_gui_strcmp(vars->model_value, "spalding"))
        ippmod[*icoebu - 1] = 0;
      else if (cs_gui_strcmp(vars->model_value, "enthalpy_st"))
        ippmod[*icoebu - 1] = 1;
      else if (cs_gui_strcmp(vars->model_value, "mixture_st"))
        ippmod[*icoebu - 1] = 2;
      else if (cs_gui_strcmp(vars->model_value, "enthalpy_mixture_st"))
        ippmod[*icoebu - 1] = 3;
      else if (cs_gui_strcmp(vars->model_value, "2-peak_adiabatic"))
        ippmod[*icolwc - 1] = 0;
      else if (cs_gui_strcmp(vars->model_value, "2-peak_enthalpy"))
        ippmod[*icolwc - 1] = 1;
      else if (cs_gui_strcmp(vars->model_value, "3-peak_adiabatic"))
        ippmod[*icolwc - 1] = 2;
      else if (cs_gui_strcmp(vars->model_value, "3-peak_enthalpy"))
        ippmod[*icolwc - 1] = 3;
      else if (cs_gui_strcmp(vars->model_value, "4-peak_adiabatic"))
        ippmod[*icolwc - 1] = 4;
      else if (cs_gui_strcmp(vars->model_value, "4-peak_enthalpy"))
        ippmod[*icolwc - 1] = 5;
      else
        bft_error(__FILE__, __LINE__, 0,
                  _("Invalid gas combustion flow model: %s.\n"),
                  vars->model_value);
    }
    else if  (cs_gui_strcmp(vars->model, "atmospheric_flows"))
    {
      if (cs_gui_strcmp(vars->model_value, "constant"))
        ippmod[*iatmos - 1] = 0;
      else if (cs_gui_strcmp(vars->model_value, "dry"))
        ippmod[*iatmos - 1] = 1;
      else if (cs_gui_strcmp(vars->model_value, "humid"))
        ippmod[*iatmos - 1] = 2;
      else
        bft_error(__FILE__, __LINE__, 0,
                  _("Invalid atmospheric flow model: %s.\n"),
                  vars->model_value);
    }
    else if  (cs_gui_strcmp(vars->model, "joule_effect"))
    {
      if (cs_gui_strcmp(vars->model_value, "joule"))
      {
        char *value = _get_joule_model();
        if (cs_gui_strcmp(value, "AC/DC"))
          ippmod[*ieljou - 1] = 1;
        else if (cs_gui_strcmp(value, "three-phase"))
          ippmod[*ieljou - 1] = 2;
        else if (cs_gui_strcmp(value, "AC/DC+Transformer"))
          ippmod[*ieljou - 1] = 3;
        else if (cs_gui_strcmp(value, "three-phase+Transformer"))
          ippmod[*ieljou - 1] = 4;
        else
          bft_error(__FILE__, __LINE__, 0,
                    _("Invalid joule model: %s.\n"),
                    vars->model_value);
        BFT_FREE(value);
      }
      else if (cs_gui_strcmp(vars->model_value, "arc"))
        ippmod[*ielarc - 1] = 2;
      else
        bft_error(__FILE__, __LINE__, 0,
                  _("Invalid electrical model: %s.\n"),
                  vars->model_value);
    }
    else if  (cs_gui_strcmp(vars->model, "compressible_model"))
    {
      if (cs_gui_strcmp(vars->model_value, "constant_gamma")){
        ippmod[*icompf - 1] = 0;
        *ieos = 1;
      }
      else if (cs_gui_strcmp(vars->model_value, "variable_gamma'")){
        ippmod[*icompf - 1] = 0;
        *ieos = 2;
      }
      else if (cs_gui_strcmp(vars->model_value, "van_der_waals")){
        ippmod[*icompf - 1] = 0;
        *ieos = 3;
      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  _("Invalid compressible model: %s.\n"),
                  vars->model_value);
    }

    /* If the model is active, one only takes the specific physics scalars */
    nscapp = _scalar_number(vars->model);
  }

  vars->nscapp = nscapp;

#if _XML_DEBUG_
  bft_printf("==>UIPPMO\n");
  if (isactiv)
  {
    bft_printf("--thermophysical model: %s\n", vars->model);
    bft_printf("--thermophysical value: %s\n", vars->model_value);
    bft_printf("--model scalars number: %i\n", vars->nscapp);
  }
#endif
}

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

/*-----------------------------------------------------------------------------
 * Indirection between the solver numbering and the XML one
 * for physical properties of the activated specific physics (solid fuel)
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicppr, UICPPR) (const int *const nclass,
                                const int *const nsalpp,
                                const int *const ippmod,
                                const int *const iccoal,
                                const int *const ipppro,
                                const int *const ipproc,
                                const int *const ieqnox,
                                const int *const ihtco2,
                                const int *const ihth2o,
                                const int *const itemp1,
                                const int *const irom1,
                                const int *const iym1,
                                const int *const ighcn1,
                                const int *const ighcn2,
                                const int *const ignoth,
                                const int *const ignh31,
                                const int *const ignh32,
                                const int *const ifhcnd,
                                const int *const ifhcnc,
                                const int *const ifnh3d,
                                const int *const ifnh3c,
                                const int *const ifnohc,
                                const int *const ifnonh,
                                const int *const ifnoch,
                                const int *const ifnoth,
                                const int *const icnohc,
                                const int *const icnonh,
                                const int *const ifhcnr,
                                const int *const icnorb,
                                const int *const igrb,
                                const int *const immel,
                                const int *const itemp2,
                                const int *const ix2,
                                const int *const irom2,
                                const int *const idiam2,
                                const int *const igmdch,
                                const int *const igmdv1,
                                const int *const igmdv2,
                                const int *const igmhet,
                                const int *const ighco2,
                                const int *const ighh2o,
                                const int *const igmsec,
                                const int *const ibcarbone,
                                const int *const iboxygen,
                                const int *const ibhydrogen)
{
  int i = 0;
  int n;
  char *name = NULL;
  char *snumpp = NULL;

  cs_var_t  *vars = cs_glob_var;

  n = vars->nprop;
  vars->nprop  += *nsalpp;
  vars->nsalpp  = *nsalpp;

  BFT_REALLOC(vars->properties_ipp,  vars->nprop, int);
  BFT_REALLOC(vars->propce,          vars->nprop, int);
  BFT_REALLOC(vars->properties_name, vars->nprop, char*);

  /* ITEMP1 */
  vars->properties_ipp[n] = ipppro[ipproc[ *itemp1 -1] -1];
  vars->propce[n] = ipproc[ *itemp1 -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("Temp_GAZ") +1, char);
  strcpy(vars->properties_name[n++], "Temp_GAZ");

  /* IROM1 */
  vars->properties_ipp[n] = ipppro[ipproc[ *irom1 -1] -1];
  vars->propce[n] = ipproc[ *irom1 -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("ROM_GAZ") +1, char);
  strcpy(vars->properties_name[n++], "ROM_GAZ");

  /*  YM_CHX1M */
  vars->properties_ipp[n] = ipppro[ipproc[ iym1[0] -1] -1];
  vars->propce[n] = ipproc[iym1[0] -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("YM_CHx1m") +1, char);
  strcpy(vars->properties_name[n++], "YM_CHx1m");

  /*  YM_CHX2M */
  vars->properties_ipp[n] = ipppro[ipproc[iym1[1] -1] -1];
  vars->propce[n] = ipproc[iym1[1] -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("YM_CHx2m") +1, char);
  strcpy(vars->properties_name[n++], "YM_CHx2m");

  /*  YM_CO */
  vars->properties_ipp[n] = ipppro[ ipproc[iym1[2] -1] -1];
  vars->propce[n] = ipproc[iym1[2] -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("YM_CO") +1, char);
  strcpy(vars->properties_name[n++], "YM_CO");

  /*  YM_H2S */
  vars->properties_ipp[n] = ipppro[ipproc[iym1[3] -1] -1];
  vars->propce[n] = ipproc[iym1[3] -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("YM_H2S") +1, char);
  strcpy(vars->properties_name[n++], "YM_H2S");

  /*  YM_H2 */
  vars->properties_ipp[n] = ipppro[ipproc[iym1[4] -1] -1];
  vars->propce[n] = ipproc[iym1[4] -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("YM_H2") +1, char);
  strcpy(vars->properties_name[n++], "YM_H2");

  /*  YM_HCN */
  vars->properties_ipp[n] = ipppro[ipproc[iym1[5] -1] -1];
  vars->propce[n] = ipproc[iym1[5] -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("YM_HCN") +1, char);
  strcpy(vars->properties_name[n++], "YM_HCN");

  /*  YM_NH3 */
  vars->properties_ipp[n] = ipppro[ipproc[iym1[6] -1] -1];
  vars->propce[n] = ipproc[iym1[6] -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("YM_NH3") +1, char);
  strcpy(vars->properties_name[n++], "YM_NH3");

  /*  YM_O2 */
  vars->properties_ipp[n] = ipppro[ipproc[iym1[7] -1] -1];
  vars->propce[n] = ipproc[iym1[7] -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("YM_O2")+1, char);
  strcpy(vars->properties_name[n++], "YM_O2");

  /*  YM_CO2 */
  vars->properties_ipp[n] = ipppro[ipproc[iym1[8] -1] -1];
  vars->propce[n] = ipproc[iym1[8] -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("YM_CO2") +1, char);
  strcpy(vars->properties_name[n++], "YM_CO2");

  /*  YM_H2O */
  vars->properties_ipp[n] = ipppro[ ipproc[ iym1[9] -1 ] -1];
  vars->propce[n] = ipproc[iym1[9] -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("YM_H2O") +1, char);
  strcpy(vars->properties_name[n++], "YM_H2O");

  /*  YM_SO2 */
  vars->properties_ipp[n] = ipppro[ipproc[iym1[10] -1] -1];
  vars->propce[n] = ipproc[iym1[10] -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("YM_SO2") +1, char);
  strcpy(vars->properties_name[n++], "YM_SO2");

  /*  YM_N2 */
  vars->properties_ipp[n] = ipppro[ ipproc[iym1[11] -1] -1];
  vars->propce[n] = ipproc[iym1[11] -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("YM_N2") +1, char);
  strcpy(vars->properties_name[n++], "YM_N2");

  if (*ieqnox == 1) {
    /* IGHCN1 */
    vars->properties_ipp[n] = ipppro[ipproc[ *ighcn1 -1] -1];
    vars->propce[n] = ipproc[ *ighcn1 -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen("EXP1") +1, char);
    strcpy(vars->properties_name[n++], "EXP1");

    /* IGHCN2 */
    vars->properties_ipp[n] = ipppro[ipproc[ *ighcn2 -1] -1];
    vars->propce[n] = ipproc[ *ighcn2 -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen("EXP2") +1, char);
    strcpy(vars->properties_name[n++], "EXP2");

    /* ignoth */
    vars->properties_ipp[n] = ipppro[ipproc[ *ignoth -1 ] -1];
    vars->propce[n] = ipproc[*ignoth -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen("EXP3") +1, char);
    strcpy(vars->properties_name[n++], "EXP3");

    /* ignh31 */
    vars->properties_ipp[n] = ipppro[ipproc[ *ignh31 -1 ] -1];
    vars->propce[n] = ipproc[*ignh31 -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen("EXP4") +1, char);
    strcpy(vars->properties_name[n++], "EXP4");

    /* ignh32 */
    vars->properties_ipp[n] = ipppro[ipproc[ *ignh32 -1 ] -1];
    vars->propce[n] = ipproc[*ignh32 -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen("EXP5") +1, char);
    strcpy(vars->properties_name[n++], "EXP5");

    /* ifhcnd */
    vars->properties_ipp[n] = ipppro[ipproc[ *ifhcnd -1 ] -1];
    vars->propce[n] = ipproc[*ifhcnd -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen("F_HCN_DEV") +1, char);
    strcpy(vars->properties_name[n++], "F_HCN_DEV");

    /* ifhcnc */
    vars->properties_ipp[n] = ipppro[ipproc[ *ifhcnc -1 ] -1];
    vars->propce[n] = ipproc[*ifhcnc -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen("F_HCN_HET") +1, char);
    strcpy(vars->properties_name[n++], "F_HCN_HET");

    /* ifnh3d */
    vars->properties_ipp[n] = ipppro[ipproc[ *ifnh3d -1 ] -1];
    vars->propce[n] = ipproc[*ifnh3d -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen("F_NH3_DEV") +1, char);
    strcpy(vars->properties_name[n++], "F_NH3_DEV");

    /* ifnh3c */
    vars->properties_ipp[n] = ipppro[ipproc[ *ifnh3c -1 ] -1];
    vars->propce[n] = ipproc[*ifnh3c -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen("F_NH3_HET") +1, char);
    strcpy(vars->properties_name[n++], "F_NH3_HET");

    /* ifnohc */
    vars->properties_ipp[n] = ipppro[ipproc[ *ifnohc -1 ] -1];
    vars->propce[n] = ipproc[*ifnohc -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen("F_NO_HCN") +1, char);
    strcpy(vars->properties_name[n++], "F_NO_HCN");

    /* ifnonh */
    vars->properties_ipp[n] = ipppro[ipproc[ *ifnonh -1 ] -1];
    vars->propce[n] = ipproc[*ifnonh -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen("F_NO_NH3") +1, char);
    strcpy(vars->properties_name[n++], "F_NO_NH3");

    /* ifnoch */
    vars->properties_ipp[n] = ipppro[ipproc[ *ifnoch -1 ] -1];
    vars->propce[n] = ipproc[*ifnoch -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen("F_NO_HET") +1, char);
    strcpy(vars->properties_name[n++], "F_NO_HET");

    /* ifnoth */
    vars->properties_ipp[n] = ipppro[ipproc[ *ifnoth -1 ] -1];
    vars->propce[n] = ipproc[*ifnoth -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen("F_NO_THE") +1, char);
    strcpy(vars->properties_name[n++], "F_NO_THE");

    /* icnohc */
    vars->properties_ipp[n] = ipppro[ipproc[ *icnohc -1 ] -1];
    vars->propce[n] = ipproc[*icnohc -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen("C_NO_HCN") +1, char);
    strcpy(vars->properties_name[n++], "C_NO_HCN");

    /* icnonh */
    vars->properties_ipp[n] = ipppro[ipproc[ *icnonh -1 ] -1];
    vars->propce[n] = ipproc[*icnonh -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen("C_NO_NH3") +1, char);
    strcpy(vars->properties_name[n++], "C_NO_NH3");

    /* ifhcnr */
    vars->properties_ipp[n] = ipppro[ipproc[ *ifhcnr -1 ] -1];
    vars->propce[n] = ipproc[*ifhcnr -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen("F_HCN_RB") +1, char);
    strcpy(vars->properties_name[n++], "F_HCN_RB");

    /* icnorb */
    vars->properties_ipp[n] = ipppro[ipproc[ *icnorb -1 ] -1];
    vars->propce[n] = ipproc[*icnorb -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen("C_NO_RB") +1, char);
    strcpy(vars->properties_name[n++], "C_NO_RB");

    /* igrb */
    vars->properties_ipp[n] = ipppro[ipproc[ *igrb -1 ] -1];
    vars->propce[n] = ipproc[*igrb -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen("EXP_RB") +1, char);
    strcpy(vars->properties_name[n++], "EXP_RB");
  }

  /* IMEL */
  vars->properties_ipp[n] = ipppro[ipproc[ *immel -1] -1];
  vars->propce[n] = ipproc[*immel -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("XM") +1, char);
  strcpy(vars->properties_name[n++], "XM");

  /* ITEMP2 loop on classes */
  BFT_MALLOC(name, strlen("Temp_CP")+1 + 2, char);
  BFT_MALLOC(snumpp, 1 + 2, char);
  strcpy(name, "Temp_CP");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumpp, "%2.2i", i+1);
    strcat(name, snumpp);

    vars->properties_ipp[n] = ipppro[ipproc[itemp2[i] -1 ] -1];
    vars->propce[n] = ipproc[itemp2[i] -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen(name) +1, char);
    strcpy(vars->properties_name[n++], name);

    strcpy(name, "Temp_CP");
  }

  /* IX2 loop on classes */
  BFT_REALLOC(name, strlen("Frm_CP")+1 + 2, char);
  strcpy(name, "Frm_CP");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumpp, "%2.2i", i+1);
    strcat(name, snumpp);

    vars->properties_ipp[n] = ipppro[ipproc[ix2[i] -1] -1];
    vars->propce[n] = ipproc[ix2[i] -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen(name) +1, char);
    strcpy(vars->properties_name[n++], name);

    strcpy(name, "Frm_CP");
  }

  /* IROM2 loop on classes */
  BFT_REALLOC(name, strlen("Rho_CP")+1 + 2, char);
  strcpy(name, "Rho_CP");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumpp, "%2.2i", i+1);
    strcat(name, snumpp);

    vars->properties_ipp[n] = ipppro[ipproc[irom2[i] -1] -1];
    vars->propce[n] = ipproc[irom2[i] -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen(name) +1, char);
    strcpy(vars->properties_name[n++], name);

    strcpy(name, "Rho_CP");
  }

  /* IDIAM2 loop on classes */
  BFT_REALLOC(name, strlen("Dia_CK")+1 + 2, char);
  strcpy(name, "Dia_CK");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumpp, "%2.2i", i+1);
    strcat(name, snumpp);

    vars->properties_ipp[n] = ipppro[ipproc[idiam2[i] -1] -1];
    vars->propce[n] = ipproc[idiam2[i] -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen(name)+1, char);
    strcpy(vars->properties_name[n++], name);

    strcpy(name, "Dia_CK");
  }

  /* IGMDCH loop on classes */
  BFT_REALLOC(name, strlen("Ga_DCH")+1 + 2, char);
  strcpy(name, "Ga_DCH");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumpp, "%2.2i", i+1);
    strcat(name, snumpp);

    vars->properties_ipp[n] = ipppro[ipproc[igmdch[i] -1] -1];
    vars->propce[n] = ipproc[igmdch[i] -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen(name)+1, char);
    strcpy(vars->properties_name[n++], name);

    strcpy(name, "Ga_DCH");
  }

  /* IGMDV1 loop on classes */
  BFT_REALLOC(name, strlen("Ga_DV1")+1 + 2, char);
  strcpy(name, "Ga_DV1");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumpp, "%2.2i", i+1);
    strcat(name, snumpp);

    vars->properties_ipp[n] = ipppro[ipproc[igmdv1[i] -1] -1];
    vars->propce[n] = ipproc[igmdv1[i] -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen(name) +1, char);
    strcpy(vars->properties_name[n++], name);

    strcpy(name, "Ga_DV1");
  }

  /* IGMDV2 loop on classes */
  BFT_REALLOC(name, strlen("Ga_DV2")+1 + 2, char);
  strcpy(name, "Ga_DV2");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumpp, "%2.2i", i+1);
    strcat(name, snumpp);

    vars->properties_ipp[n] = ipppro[ipproc[igmdv2[i] -1 ] -1];
    vars->propce[n] = ipproc[igmdv2[i] -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen(name) +1, char);
    strcpy(vars->properties_name[n++], name);

    strcpy(name, "Ga_DV2");
  }

  /* IGMHET loop on classes */
  BFT_REALLOC(name, strlen("Ga_HET_O2")+1 + 2, char);
  strcpy(name, "Ga_HET_O2");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumpp, "%2.2i", i+1);
    strcat(name, snumpp);

    vars->properties_ipp[n] = ipppro[ipproc[igmhet[i] -1] -1];
    vars->propce[n] = ipproc[igmhet[i] -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen(name) +1, char);
    strcpy(vars->properties_name[n++], name);

    strcpy(name, "Ga_HET_O2");
  }

  if (*ihtco2 == 1) {
    /* IGHCO2 loop on classes */
    BFT_REALLOC(name, strlen("Ga_HET_CO2")+1 + 2, char);
    strcpy(name, "Ga_HET_CO2");
    for (i = 0; i < *nclass; i++)
    {
      sprintf(snumpp, "%2.2i", i+1);
      strcat(name, snumpp);

      vars->properties_ipp[n] = ipppro[ipproc[ighco2[i] -1] -1];
      vars->propce[n] = ipproc[ighco2[i] -1] -1;
      BFT_MALLOC(vars->properties_name[n], strlen(name) +1, char);
      strcpy(vars->properties_name[n++], name);

      strcpy(name, "Ga_HET_CO2");
    }
  }

  if (*ihth2o == 1) {
    /* IGNH2O loop on classes */
    BFT_REALLOC(name, strlen("Ga_HET_H2O")+1 + 2, char);
    strcpy(name, "Ga_HET_H2O");
    for (i = 0; i < *nclass; i++)
    {
      sprintf(snumpp, "%2.2i", i+1);
      strcat(name, snumpp);

      vars->properties_ipp[n] = ipppro[ipproc[ighh2o[i] -1] -1];
      vars->propce[n] = ipproc[ighh2o[i] -1] -1;
      BFT_MALLOC(vars->properties_name[n], strlen(name) +1, char);
      strcpy(vars->properties_name[n++], name);

      strcpy(name, "Ga_HET_H2O");
    }
  }

  if (ippmod[*iccoal -1] == 1) {
    /* IGMSEC loop on classes */
    BFT_REALLOC(name, strlen("Ga_SEC")+1 + 2, char);
    strcpy(name, "Ga_SEC");
    for (i = 0; i < *nclass; i++)
    {
      sprintf(snumpp, "%2.2i", i+1);
      strcat(name, snumpp);

      vars->properties_ipp[n] = ipppro[ipproc[igmsec[i] -1] -1];
      vars->propce[n] = ipproc[igmsec[i] -1] -1;
      BFT_MALLOC(vars->properties_name[n], strlen(name) +1, char);
      strcpy(vars->properties_name[n++], name);

      strcpy(name, "Ga_SEC");
    }
  }

  /* Bilan_C */
  vars->properties_ipp[n] = ipppro[ipproc[ *ibcarbone -1] -1];
  vars->propce[n] = ipproc[ *ibcarbone -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("Bilan_C")+1, char);
  strcpy(vars->properties_name[n++], "Bilan_C");

  /* Bilan_O */
  vars->properties_ipp[n] = ipppro[ipproc[*iboxygen -1] -1];
  vars->propce[n] = ipproc[ *iboxygen -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("Bilan_O")+1, char);
  strcpy(vars->properties_name[n++], "Bilan_O");

  /* Bilan_H */
  vars->properties_ipp[n] = ipppro[ipproc[ *ibhydrogen -1] -1];
  vars->propce[n] = ipproc[ *ibhydrogen -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("Bilan_H")+1, char);
  strcpy(vars->properties_name[n++], "Bilan_H");

  BFT_FREE(name);
  BFT_FREE(snumpp);

  if (n != vars->nsalpp)
    bft_error(__FILE__, __LINE__, 0,
              _("number of properties is not correct: %i instead of: %i\n"),
              n, vars->nsalpp);

#if _XML_DEBUG_
  bft_printf("==>UICPPR\n");
  bft_printf("-->nombre de proprietes = %i\n", vars->nprop);
  for (i=0 ; i<vars->nprop ; i++)
    bft_printf("-->properties_ipp[%i]: %i propce[%i]: %i "
               "properties_name[%i]: %s\n",
               i, vars->properties_ipp[i],
               i, vars->propce[i],
               i, vars->properties_name[i]);
#endif
}

/*-----------------------------------------------------------------------------
 * Indirection between the solver numbering and the XML one
 * for physical properties of the activated specific physics (gaz combustion)
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicopr, UICOPR) (const int *const nsalpp,
                                const int *const ippmod,
                                const int *const ipppro,
                                const int *const ipproc,
                                const int *const icolwc,
                                const int *const iirayo,
                                const int *const itemp,
                                const int *const imam,
                                const int *const iym,
                                const int *const ickabs,
                                const int *const it4m,
                                const int *const it3m,
                                const int *const itsc,
                                const int *const irhol,
                                const int *const iteml,
                                const int *const ifmel,
                                const int *const ifmal,
                                const int *const iampl,
                                const int *const itscl,
                                const int *const imaml)
{
  int n, idirac;
  int ndirac = 0;
  char *name = NULL;
  char *snumpp = NULL;

  cs_var_t  *vars = cs_glob_var;

  n = vars->nprop;
  vars->nprop  += *nsalpp;
  vars->nsalpp  = *nsalpp;

  BFT_REALLOC(vars->properties_ipp,  vars->nprop, int);
  BFT_REALLOC(vars->propce,          vars->nprop, int);
  BFT_REALLOC(vars->properties_name, vars->nprop, char*);
  BFT_MALLOC(snumpp, 1 + 2, char);
  /* Source Term */
  if (ippmod[*icolwc -1] >= 0) {
    vars->properties_ipp[n] = ipppro[ipproc[ *itsc -1] -1];
    vars->propce[n] = ipproc[ *itsc -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen("T.SOURCE") +1, char);
    strcpy(vars->properties_name[n++], "T.SOURCE");
  }

  /* Temperature */
  vars->properties_ipp[n] = ipppro[ipproc[ *itemp -1] -1];
  vars->propce[n] = ipproc[ *itemp -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("Temperature") +1, char);
  strcpy(vars->properties_name[n++], "Temperature");

  if (ippmod[*icolwc -1] >= 0) {
    /* Mass Molaire */
    vars->properties_ipp[n] = ipppro[ipproc[ *imam -1] -1];
    vars->propce[n] = ipproc[ *imam -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen("Mas_Mol") +1, char);
    strcpy(vars->properties_name[n++], "Mas_Mol");
  }

  /* Fuel mass fraction */
  vars->properties_ipp[n] = ipppro[ipproc[iym[0] -1] -1];
  vars->propce[n] = ipproc[iym[0] -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("YM_Fuel") +1, char);
  strcpy(vars->properties_name[n++], "YM_Fuel");

  /*  Oxydizer Mass fraction */
  vars->properties_ipp[n] = ipppro[ipproc[iym[1] -1] -1];
  vars->propce[n] = ipproc[iym[1] -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("YM_Oxyd") +1, char);
  strcpy(vars->properties_name[n++], "YM_Oxyd");

  /*  Product Mass fraction */
  vars->properties_ipp[n] = ipppro[ ipproc[iym[2] -1] -1];
  vars->propce[n] = ipproc[iym[2] -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("YM_Prod") +1, char);
  strcpy(vars->properties_name[n++], "YM_Prod");
  if (ippmod[*icolwc -1] > 0) {
    if (ippmod[*icolwc -1] == 0 || ippmod[*icolwc -1] == 1) {
      ndirac = 2;
    } else if (ippmod[*icolwc -1] == 2 || ippmod[*icolwc -1] == 3) {
      ndirac = 3;
    } else if (ippmod[*icolwc -1] == 4 || ippmod[*icolwc -1] == 5) {
      ndirac = 4;
    }
    BFT_MALLOC(name, strlen("RHOL")+1 + 2, char);
    strcpy(name, "RHOL");
    for (idirac = 0; idirac < ndirac; idirac++) {
      vars->properties_ipp[n] = ipppro[ipproc[irhol[idirac] -1] -1];
      vars->propce[n] = ipproc[irhol[idirac] -1] -1;
      BFT_MALLOC(vars->properties_name[n], strlen(name) +1 +2, char);
      sprintf(snumpp, "%2.2i", idirac+1);
      strcat(name, snumpp);
      strcpy(vars->properties_name[n++], name);
      strcpy(name, "TEML");


      vars->properties_ipp[n] = ipppro[ipproc[iteml[idirac] -1] -1];
      vars->propce[n] = ipproc[iteml[idirac] -1] -1;
      BFT_MALLOC(vars->properties_name[n], strlen(name) +1 +2, char);
      sprintf(snumpp, "%2.2i", idirac+1);
      strcat(name, snumpp);
      strcpy(vars->properties_name[n++], name);
      strcpy(name, "FMEL");

      vars->properties_ipp[n] = ipppro[ipproc[ifmel[idirac] -1] -1];
      vars->propce[n] = ipproc[ifmel[idirac] -1] -1;
      BFT_MALLOC(vars->properties_name[n], strlen(name) +1 +2, char);
      sprintf(snumpp, "%2.2i", idirac+1);
      strcat(name, snumpp);
      strcpy(vars->properties_name[n++], name);
      strcpy(name, "FMAL");

      vars->properties_ipp[n] = ipppro[ipproc[ifmal[idirac] -1] -1];
      vars->propce[n] = ipproc[ifmal[idirac] -1] -1;
      BFT_MALLOC(vars->properties_name[n], strlen(name) +1 +2, char);
      sprintf(snumpp, "%2.2i", idirac+1);
      strcat(name, snumpp);
      strcpy(vars->properties_name[n++], name);
      strcpy(name, "AMPL");

      vars->properties_ipp[n] = ipppro[ipproc[iampl[idirac] -1] -1];
      vars->propce[n] = ipproc[iampl[idirac] -1] -1;
      BFT_MALLOC(vars->properties_name[n], strlen(name)+1+2, char);
      sprintf(snumpp, "%2.2i", idirac+1);
      strcat(name, snumpp);
      strcpy(vars->properties_name[n++], name);
      strcpy(name, "TSCL");

      vars->properties_ipp[n] = ipppro[ipproc[itscl[idirac] -1] -1];
      vars->propce[n] = ipproc[itscl[idirac] -1] -1;
      BFT_MALLOC(vars->properties_name[n], strlen(name) +1 +2, char);
      sprintf(snumpp, "%2.2i", idirac+1);
      strcat(name, snumpp);
      strcpy(vars->properties_name[n++], name);
      strcpy(name, "MAML");

      vars->properties_ipp[n] = ipppro[ipproc[imaml[idirac] -1] -1];
      vars->propce[n] = ipproc[imaml[idirac] -1] -1;
      BFT_MALLOC(vars->properties_name[n], strlen(name) +1 +2, char);
      sprintf(snumpp, "%2.2i", idirac+1);
      strcat(name, snumpp);
      strcpy(vars->properties_name[n++], name);
      strcpy(name, "RHOL");
    }
  }
  if (*iirayo > 0) {
    /* Absoption coefficient */
    vars->properties_ipp[n] = ipppro[ipproc[ *ickabs -1] -1];
    vars->propce[n] = ipproc[ *ickabs -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen("KABS") +1, char);
    strcpy(vars->properties_name[n++], "KABS");

    /* Term T^4 */
    vars->properties_ipp[n] = ipppro[ipproc[ *it4m -1] -1];
    vars->propce[n] =ipproc[ *it4m -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen("TEMP4") +1, char);
    strcpy(vars->properties_name[n++], "TEMP4");

    /* Term T^3 */
    vars->properties_ipp[n] = ipppro[ipproc[ *it3m -1] -1];
    vars->propce[n] = ipproc[ *it3m -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen("TEMP3")+1, char);
    strcpy(vars->properties_name[n++], "TEMP3");
  }

  BFT_FREE(name);
  BFT_FREE(snumpp);

  if (n != vars->nsalpp)
    bft_error(__FILE__, __LINE__, 0,
              _("number of properties is not correct: %i instead of: %i\n"),
              n, vars->nsalpp);

#if _XML_DEBUG_
  bft_printf("==>UICPPR\n");
  bft_printf("-->nombre de proprietes = %i\n", vars->nprop);
  for (i=0 ; i<vars->nprop ; i++)
    bft_printf("-->properties_ipp[%i]: %i propce[%i]: %i "
               "properties_name[%i]: %s\n",
               i, vars->properties_ipp[i],
               i, vars->propce[i],
               i, vars->properties_name[i]);
#endif
}

/*------------------------------------------------------------------------------
 * Indirection between the solver numbering and the XML one
 * for the model scalar (gas combustion)
 *----------------------------------------------------------------------------*/
void CS_PROCF (uicosc, UICOSC) (const int *const ippmod,
                                const int *const icolwc,
                                const int *const icoebu,
                                const int *const icod3p,
                                const int *const ihm,
                                const int *const ifm,
                                const int *const ifp2m,
                                const int *const iygfm,
                                const int *const iyfm,
                                const int *const iyfp2m,
                                const int *const icoyfp)
{
  cs_var_t  *vars = cs_glob_var;
  char *label = NULL;

  if (vars->nscaus > 0) {
    BFT_REALLOC(vars->label, vars->nscapp + vars->nscaus, char*);
  } else {
    BFT_MALLOC(vars->label, vars->nscapp, char*);
  }
  //   model D3P
  if (ippmod[*icod3p-1] >=0) {
    label = _scalar_name_label("gas_combustion", "Fra_MEL");
    BFT_MALLOC(vars->label[*ifm -1], strlen(label)+1, char);
    strcpy(vars->label[*ifm -1], label);
    BFT_FREE(label);

    label = _scalar_name_label("gas_combustion", "Var_FMe");
    BFT_MALLOC(vars->label[*ifp2m -1], strlen(label)+1, char);
    strcpy(vars->label[*ifp2m -1], label);
    BFT_FREE(label);

    if (ippmod[*icod3p-1] == 1 ) {
      label = _scalar_name_label("gas_combustion", "Enthalpy");
      BFT_MALLOC(vars->label[*ihm -1], strlen(label)+1, char);
      strcpy(vars->label[*ihm -1], label);
      BFT_FREE(label);
    }
  }
  // model EBU
  if (ippmod[*icoebu-1] >= 0) {
    label = _scalar_name_label("gas_combustion", "Fra_GF");
    BFT_MALLOC(vars->label[*iygfm -1], strlen(label)+1, char);
    strcpy(vars->label[*iygfm -1], label);
    BFT_FREE(label);

    if (ippmod[*icoebu-1] >= 2) {
      label = _scalar_name_label("gas_combustion", "Fra_MEL");
      BFT_MALLOC(vars->label[*ifm -1], strlen(label)+1, char);
      strcpy(vars->label[*ifm -1], label);
      BFT_FREE(label);
    }

    if (ippmod[*icoebu-1] == 1 || ippmod[*icoebu-1] == 3) {
      label = _scalar_name_label("gas_combustion", "Enthalpy");
      BFT_MALLOC(vars->label[*ihm -1], strlen(label)+1, char);
      strcpy(vars->label[*ihm -1], label);
      BFT_FREE(label);
    }
  }
  // model LWC
  if (ippmod[*icolwc-1] >= 0) {
    label = _scalar_name_label("gas_combustion", "Fra_MEL");
    BFT_MALLOC(vars->label[*ifm -1], strlen(label)+1, char);
    strcpy(vars->label[*ifm -1], label);
    BFT_FREE(label);

    label = _scalar_name_label("gas_combustion", "Var_FMe");
    BFT_MALLOC(vars->label[*ifp2m -1], strlen(label)+1, char);
    strcpy(vars->label[*ifp2m -1], label);
    BFT_FREE(label);

    label = _scalar_name_label("gas_combustion", "Fra_Mas");
    BFT_MALLOC(vars->label[*iyfm -1], strlen(label)+1, char);
    strcpy(vars->label[*iyfm -1], label);
    BFT_FREE(label);

    label = _scalar_name_label("gas_combustion", "Var_FMa");
    BFT_MALLOC(vars->label[*iyfp2m -1], strlen(label)+1, char);
    strcpy(vars->label[*iyfp2m -1], label);
    BFT_FREE(label);
  }

  if (ippmod[*icolwc-1] >= 2) {
    label = _scalar_name_label("gas_combustion", "COYF_PP4");
    BFT_MALLOC(vars->label[*icoyfp -1], strlen(label)+1, char);
    strcpy(vars->label[*icoyfp -1], label);
    BFT_FREE(label);
  }
  if (ippmod[*icolwc-1] == 1 || ippmod[*icolwc-1] == 3 || ippmod[*icolwc-1] == 5) {
    label = _scalar_name_label("gas_combustion", "Enthalpy");
    BFT_MALLOC(vars->label[*ihm -1], strlen(label)+1, char);
    strcpy(vars->label[*ihm -1], label);
    BFT_FREE(label);
  }


#if _XML_DEBUG_
  bft_printf("==>UICPSC\n");
  for (i=0; i< vars->nscaus+vars->nscapp; i++)
    bft_printf("--label of scalar[%i]: %s\n", i, vars->label[i]);
#endif

}

/*----------------------------------------------------------------------------
 * Electrical model : read parameters
 *
 * Fortran Interface:
 *
 * subroutine uieli1
 * *****************
 * integer         ieljou    -->   joule model
 * integer         ielarc    -->   arc model
 * integer         ielcor    <--   scaling electrical variables
 * double          couimp    <--   imposed current intensity
 * double          puisim    <--   imposed power
 * integer         modrec    <--   scaling type for electric arc
 * integer         idreca    <--   current density component used to scaling
 *                                 (modrec ==2)
 * char            crit_reca <--   define criteria for plane used to scaling (modrec ==2)
 *----------------------------------------------------------------------------*/

void CS_PROCF (uieli1, UIELI1) (const int    *const ieljou,
                                const int    *const ielarc,
                                      int    *const ielcor,
                                      double *const couimp,
                                      double *const puisim,
                                      int    *const modrec,
                                      int    *const idreca,
                                      double *const crit_reca)
{
  char *path   = NULL;
  int   status = 0;
  double   result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "thermophysical_models",
                                  "joule_effect",
                                  "variable_scaling");

  cs_xpath_add_attribute(&path, "status");
  if (cs_gui_get_status(path, &status))
    *ielcor = status;
  BFT_FREE(path);

  if (*ieljou > 0) {
    path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 3, "thermophysical_models",
                                    "joule_effect",
                                    "imposed_power");

    cs_xpath_add_function_text(&path);
    if (!cs_gui_get_double(path, &result))
      bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);
    *puisim = result;
    BFT_FREE(path);
  }

  if (*ielarc > 0) {
    path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 3, "thermophysical_models",
                                    "joule_effect",
                                    "imposed_current");
    cs_xpath_add_function_text(&path);
    if (!cs_gui_get_double(path, &result))
      bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);
    *couimp = result;
    BFT_FREE(path);

    if (*ielcor > 0) {
      path = cs_xpath_init_path();
      char *choice;

      cs_xpath_add_elements(&path, 3, "thermophysical_models",
                                      "joule_effect",
                                      "recal_model");
      cs_xpath_add_attribute(&path, "model");
      choice = cs_gui_get_attribute_value(path);
      if (cs_gui_strcmp(choice, "general_case"))
        *modrec = 1;
      else if (cs_gui_strcmp(choice, "plane_define"))
        *modrec = 2;
      else if (cs_gui_strcmp(choice, "user"))
        *modrec = 3;
      else
        bft_error(__FILE__, __LINE__, 0, _("Invalid model : %s\n"), choice);

      BFT_FREE(choice);

      if (*modrec == 2) {
        path = cs_xpath_init_path();
        cs_xpath_add_elements(&path, 4, "thermophysical_models",
                                        "joule_effect",
                                        "recal_model",
                                        "direction");
        cs_xpath_add_function_text(&path);
        choice = cs_gui_get_text_value(path);
        if (cs_gui_strcmp(choice, "X"))
          *idreca = 1;
        else if (cs_gui_strcmp(choice, "Y"))
          *idreca = 2;
        else
          *idreca = 3;
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
        crit_reca[0] = val;

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
        crit_reca[1] = val;

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
        crit_reca[2] = val;

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
        crit_reca[3] = val;

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
        crit_reca[4] = val;
      }
    }
    BFT_FREE(path);
  }

#if _XML_DEBUG_
  bft_printf("==>UIELI1\n");
  bft_printf("--ielcor  = %i\n", *ielcor);
  bft_printf("--puisim  = %f\n", *puisim);
  bft_printf("--couimp  = %f\n", *couimp);
  bft_printf("--modrec  = %f\n", *modrec);
#endif
}

/*----------------------------------------------------------------------------
 * Electrical model : define plane for elreca
 *
 * Fortran Interface:
 *
 * subroutine uielrc
 * *****************
 * integer         izreca    <--   define plane used to scaling (modrec ==2)
 * char            crit_reca <--   define criteria for plane used to scaling (modrec ==2)
 *----------------------------------------------------------------------------*/

void CS_PROCF (uielrc, UIELRC) (int    *const izreca,
                                double *const crit_reca)
{
  /* build list of cells */
  char *crit = NULL;

  cs_lnum_t   n_selected_faces = 0;
  cs_lnum_t  *selected_faces = NULL;

  BFT_MALLOC(crit, 66, char);

  char cVal[10];

  strcpy(crit, "plane[");
  sprintf(cVal, "%f", crit_reca[0]);
  strcat(crit, cVal);
  strcat(crit, ",");
  sprintf(cVal, "%f", crit_reca[1]);
  strcat(crit, cVal);
  strcat(crit, ",");
  sprintf(cVal, "%f", crit_reca[2]);
  strcat(crit, cVal);
  strcat(crit, ",");
  sprintf(cVal, "%f", crit_reca[3]);
  strcat(crit, cVal);
  strcat(crit, ",epsilon=");
  sprintf(cVal, "%6f", crit_reca[4]);
  strcat(crit, cVal);
  strcat(crit, "]");

  BFT_MALLOC(selected_faces, cs_glob_mesh->n_i_faces, cs_lnum_t);

  cs_selector_get_i_face_list(crit,
                              &n_selected_faces,
                              selected_faces);

  for (int j=0; j < n_selected_faces; j++)
    izreca[selected_faces[j]] = 1;

  BFT_FREE(selected_faces);
  BFT_FREE(crit);

}

/*------------------------------------------------------------------------------
 * Indirection between the solver numbering and the XML one
 * for the model scalar
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicpsc, UICPSC) (const int *const ncharb,
                                const int *const nclass,
                                const int *const noxyd,
                                const int *const ippmod,
                                const int *const iccoal,
                                const int *const ieqnox,
                                const int *const ieqco2,
                                const int *const ihtco2,
                                const int *const ihth2o,
                                const int *const ihm,
                                const int *const inp,
                                const int *const ixch,
                                const int *const ixck,
                                const int *const ixwt,
                                const int *const ih2,
                                const int *const if1m,
                                const int *const if2m,
                                const int *const if4m,
                                const int *const if5m,
                                const int *const if6m,
                                const int *const if7m,
                                const int *const if8m,
                                const int *const ifvp2m,
                                const int *const iyco2,
                                const int *const if9m,
                                const int *const iyhcn,
                                const int *const iyno,
                                const int *const ihox,
                                const int *const iynh3)
{
  int i;
  char *name = NULL;
  char *snumsca = NULL;
  char *label = NULL;

  cs_var_t  *vars = cs_glob_var;

  if (vars->nscaus > 0) {
    BFT_REALLOC(vars->label, vars->nscapp + vars->nscaus, char*);
  } else {
    BFT_MALLOC(vars->label, vars->nscapp, char*);
  }

  /* IHM */
  label = _scalar_name_label("solid_fuels", "Enthalpy");
  BFT_MALLOC(vars->label[*ihm -1], strlen(label)+1, char);
  strcpy(vars->label[*ihm -1], label);
  BFT_FREE(label);

  /* Loop on classes IH2, INP, IXCH, IXCK */
  BFT_MALLOC(snumsca, 1 + 2, char);

  /* IH2 */
  BFT_MALLOC(name, strlen("ENT_CP")+1 + 2, char);
  for (i = 0; i < *nclass; i++) {
    strcpy(name, "ENT_CP");
    sprintf(snumsca,"%2.2i", i+1);
    strcat(name, snumsca);
    label = _scalar_name_label("solid_fuels", name);

    BFT_MALLOC(vars->label[ih2[i] -1], strlen(label)+1, char);
    strcpy(vars->label[ih2[i] -1], label);
    BFT_FREE(label);
  }

  /* INP */
  BFT_REALLOC(name, strlen("NP_CP") + 1 + 2, char);
  for (i = 0; i < *nclass; i++) {
    strcpy(name, "NP_CP");
    sprintf(snumsca,"%2.2i", i+1);
    strcat(name, snumsca);

    label = _scalar_name_label("solid_fuels", name);
    BFT_MALLOC(vars->label[inp[i] -1], strlen(label)+1, char);
    strcpy(vars->label[inp[i] -1], label);
    BFT_FREE(label);
  }

  /* IXCH */
  BFT_REALLOC(name, strlen("XCH_CP")+1 + 2, char);
  for (i = 0; i < *nclass; i++) {
    strcpy(name, "XCH_CP");
    sprintf(snumsca,"%2.2i", i+1);
    strcat(name, snumsca);

    label = _scalar_name_label("solid_fuels", name);
    BFT_MALLOC(vars->label[ixch[i] -1], strlen(label)+1, char);
    strcpy(vars->label[ixch[i] -1], label);
    BFT_FREE(label);
  }

  /* IXCK */
  BFT_REALLOC(name, strlen("XCK_CP")+1 + 2, char);
  for (i = 0; i < *nclass; i++) {
    strcpy(name, "XCK_CP");
    sprintf(snumsca,"%2.2i", i+1);
    strcat(name, snumsca);

    label = _scalar_name_label("solid_fuels", name);
    BFT_MALLOC(vars->label[ixck[i] -1], strlen(label)+1, char);
    strcpy(vars->label[ixck[i] -1], label);
    BFT_FREE(label);
  }

  /* Loop on coals IFM1 IFM2 */

  BFT_REALLOC(name, strlen("Fr_MV1")+1 + 2, char);
  for (i = 0; i < *ncharb; i++) {
    strcpy(name, "Fr_MV1");
    sprintf(snumsca,"%2.2i",i+1);
    strcat(name, snumsca);

    label = _scalar_name_label("solid_fuels", name);
    BFT_MALLOC(vars->label[if1m[i] -1], strlen(label)+1, char);
    strcpy(vars->label[if1m[i] -1], label);
    BFT_FREE(label);
  }

  BFT_REALLOC(name, strlen("Fr_MV2")+1 + 2, char);
  for (i = 0; i < *ncharb; i++) {
    strcpy(name, "Fr_MV2");
    sprintf(snumsca,"%2.2i",i+1);
    strcat(name, snumsca);

    label = _scalar_name_label("solid_fuels", name);
    BFT_MALLOC(vars->label[if2m[i] -1], strlen(label)+1, char);
    strcpy(vars->label[if2m[i] -1], label);
    BFT_FREE(label);
  }

  /* IF7M */
  label = _scalar_name_label("solid_fuels", "Fr_HET_O2");
  BFT_MALLOC(vars->label[*if7m -1], strlen(label)+1, char);
  strcpy(vars->label[*if7m -1], label);
  BFT_FREE(label);

  if (*ihtco2 == 1) {
    /* IF3MC2 */
    label = _scalar_name_label("solid_fuels", "Fr_HET_CO2");
    BFT_MALLOC(vars->label[*if8m -1], strlen(label)+1, char);
    strcpy(vars->label[*if8m -1], label);
    BFT_FREE(label);
  }

  if (*ihth2o == 1) {
    /* IF3MC2 */
    label = _scalar_name_label("solid_fuels", "Fr_HET_H2O");
    BFT_MALLOC(vars->label[*if9m -1], strlen(label)+1, char);
    strcpy(vars->label[*if9m -1],label );
    BFT_FREE(label);
  }

  /* IFVP2M */
  label = _scalar_name_label("solid_fuels", "Var_F1F2");
  BFT_MALLOC(vars->label[*ifvp2m -1], strlen(label)+1, char);
  strcpy(vars->label[*ifvp2m -1], label);
  BFT_FREE(label);

  if (ippmod[*iccoal -1] == 1)
  {
    /* IXWT */
    BFT_MALLOC(name, strlen("XWT_CP")+1 + 2, char);
    for (i = 0; i < *nclass; i++)
    {
      strcpy(name, "XWT_CP");
      sprintf(snumsca,"%2.2i", i+1);
      strcat(name, snumsca);

      label = _scalar_name_label("solid_fuels", name);
      BFT_MALLOC(vars->label[ixwt[i] -1], strlen(label)+1, char);
      strcpy(vars->label[ixwt[i] -1], label);
      BFT_FREE(label);
    }

    /* IF6M */
    label = _scalar_name_label("solid_fuels", "FR_H20");
    BFT_MALLOC(vars->label[*if6m -1], strlen(label)+1, char);
    strcpy(vars->label[*if6m -1], label);
    BFT_FREE(label);
  }

  if (*noxyd >= 2) {
    /* IF4M */
    label = _scalar_name_label("solid_fuels", "FR_OXYD2");
    BFT_MALLOC(vars->label[*if4m -1], strlen(label)+1, char);
    strcpy(vars->label[*if4m -1], label);
    BFT_FREE(label);
  }

  if (*noxyd == 3) {
    /* IF5M */
    label = _scalar_name_label("solid_fuels", "FR_OXYD3");
    BFT_MALLOC(vars->label[*if5m -1], strlen(label)+1, char);
    strcpy(vars->label[*if5m -1], label);
    BFT_FREE(label);
  }

  if (*ieqco2 == 1) {
    /* IYCO2 */
    label = _scalar_name_label("solid_fuels", "FR_CO2");
    BFT_MALLOC(vars->label[*iyco2 -1], strlen(label)+1, char);
    strcpy(vars->label[*iyco2 -1], label);
    BFT_FREE(label);
  }

  if (*ieqnox == 1) {
    /* FR_HCN */
    label = _scalar_name_label("solid_fuels", "FR_HCN");
    BFT_MALLOC(vars->label[*iyhcn -1], strlen(label)+1, char);
    strcpy(vars->label[*iyhcn -1], label);
    BFT_FREE(label);

    /* FR_NO */
    label = _scalar_name_label("solid_fuels", "FR_NO");
    BFT_MALLOC(vars->label[*iyno -1], strlen(label)+1, char);
    strcpy(vars->label[*iyno -1], label);
    BFT_FREE(label);

    /* Enth_Ox */
    label = _scalar_name_label("solid_fuels", "Enth_Ox");
    BFT_MALLOC(vars->label[*ihox -1], strlen(label)+1, char);
    strcpy(vars->label[*ihox -1], label);
    BFT_FREE(label);

    /* FR_NH3 */
    label = _scalar_name_label("solid_fuels", "FR_NH3");
    BFT_MALLOC(vars->label[*iynh3 -1], strlen(label)+1, char);
    strcpy(vars->label[*iynh3 -1], label);
    BFT_FREE(label);
  }

  BFT_FREE(name);
  BFT_FREE(snumsca);

#if _XML_DEBUG_
  bft_printf("==>UICPSC\n");
  for (i=0; i< vars->nscaus+vars->nscapp; i++)
    bft_printf("--label of scalar[%i]: %s\n", i, vars->label[i]);
#endif

}

/*------------------------------------------------------------------------------
 * Indirection between the solver numbering and the XML one
 * for the model scalar (compressible model)
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicfsc, UICFSC) (const int *const ienerg,
                                const int *const itempk)
{
  cs_var_t *vars = cs_glob_var;
  char *label = NULL;

  if (vars->nscaus > 0) {
    BFT_REALLOC(vars->label, vars->nscapp + vars->nscaus, char*);
  } else {
    BFT_MALLOC(vars->label, vars->nscapp, char*);
  }

  label = _scalar_name_label("compressible_model", "EnergieT");
  BFT_MALLOC(vars->label[*ienerg -1], strlen(label)+1, char);
  strcpy(vars->label[*ienerg -1], label);
  BFT_FREE(label);

  label = _scalar_name_label("compressible_model", "TempK");
  BFT_MALLOC(vars->label[*itempk -1], strlen(label)+1, char);
  strcpy(vars->label[*itempk -1], label);
  BFT_FREE(label);

#if _XML_DEBUG_
  bft_printf("==>UICPSC\n");
  for (i=0; i< vars->nscaus+vars->nscapp; i++)
    bft_printf("--label of scalar[%i]: %s\n", i, vars->label[i]);
#endif

}

/*-----------------------------------------------------------------------------
 * Indirection between the solver numbering and the XML one
 * for physical properties of the activated specific physics (electrical model)
 *----------------------------------------------------------------------------*/
void CS_PROCF (uielpr, UIELPR) (const int *const nsalpp,
                                const int *const ippmod,
                                const int *const ipppro,
                                const int *const ipproc,
                                const int *const ieljou,
                                const int *const ielarc,
                                const int *const itemp,
                                const int *const iefjou,
                                const int *const idjr,
                                const int *const idji,
                                const int *const ilapla,
                                const int *const idrad,
                                const int *const ixkabe)
{
  int n;
  char *name = NULL;
  char *snumpp = NULL;

  cs_var_t  *vars = cs_glob_var;

  n = vars->nprop;
  vars->nprop  += *nsalpp;
  vars->nsalpp  = *nsalpp;

  BFT_REALLOC(vars->properties_ipp,  vars->nprop, int);
  BFT_REALLOC(vars->propce,          vars->nprop, int);
  BFT_REALLOC(vars->properties_name, vars->nprop, char*);
  BFT_MALLOC(snumpp, 1 + 1, char);

  /* Temperature */
  vars->properties_ipp[n] = ipppro[ipproc[*itemp -1] -1];
  vars->propce[n] = ipproc[*itemp -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("Temperature") +1, char);
  strcpy(vars->properties_name[n++], "Temperature");

  /* Power */
  vars->properties_ipp[n] = ipppro[ipproc[*iefjou -1] -1];
  vars->propce[n] = ipproc[*iefjou -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("PuisJoul") +1, char);
  strcpy(vars->properties_name[n++], "PuisJoul");

  for (int idimve = 0; idimve < 3; idimve++) {
    vars->properties_ipp[n] = ipppro[ipproc[idjr[idimve] -1] -1];
    vars->propce[n] = ipproc[idjr[idimve] -1] -1;
    BFT_MALLOC(name, strlen("Cour_re") +1 + 1, char);
    BFT_MALLOC(vars->properties_name[n], strlen(name) +1, char);
    strcpy(name, "Cour_re");
    sprintf(snumpp,"%1.1i", idimve+1);
    strcat(name, snumpp);
    strcpy(vars->properties_name[n++], name);
  }

  if (ippmod[*ieljou - 1] == 2 || ippmod[*ieljou - 1] == 4)
    for (int idimve = 0; idimve < 3; idimve++) {
      vars->properties_ipp[n] = ipppro[ipproc[idji[idimve] -1] -1];
      vars->propce[n] = ipproc[idji[idimve] -1] -1;
      BFT_MALLOC(name, strlen("CouImag") +1 + 1, char);
      BFT_MALLOC(vars->properties_name[n], strlen(name) +1, char);
      strcpy(name, "CouImag");
      sprintf(snumpp,"%1.1i", idimve+1);
      strcat(name, snumpp);
      strcpy(vars->properties_name[n++], name);
    }

  if (ippmod[*ielarc - 1] >= 1) {
    for (int idimve = 0; idimve < 3; idimve++) {
      vars->properties_ipp[n] = ipppro[ipproc[ilapla[idimve] -1] -1];
      vars->propce[n] = ipproc[ilapla[idimve] -1] -1;
      BFT_MALLOC(name, strlen("For_Lap") +1 + 1, char);
      BFT_MALLOC(vars->properties_name[n], strlen(name) +1, char);
      strcpy(name, "For_Lap");
      sprintf(snumpp,"%1.1i", idimve+1);
      strcat(name, snumpp);
      strcpy(vars->properties_name[n++], name);
    }

    if (*ixkabe == 1) {
      vars->properties_ipp[n] = ipppro[ipproc[*idrad -1] -1];
      vars->propce[n] = ipproc[*idrad -1] -1;
      BFT_MALLOC(vars->properties_name[n], strlen("Coef_Abso") +1, char);
      strcpy(vars->properties_name[n++], "Coef_Abso");
    } else if (*ixkabe == 2) {
      vars->properties_ipp[n] = ipppro[ipproc[*idrad -1] -1];
      vars->propce[n] = ipproc[*idrad -1] -1;
      BFT_MALLOC(vars->properties_name[n], strlen("TS_radia") +1, char);
      strcpy(vars->properties_name[n++], "TS_radia");
    }
  }

  BFT_FREE(name);
  BFT_FREE(snumpp);

  if (n != vars->nsalpp)
    bft_error(__FILE__, __LINE__, 0,
              _("number of properties is not correct: %i instead of: %i\n"),
              n, vars->nsalpp);

#if _XML_DEBUG_
  bft_printf("==>UICPPR\n");
  bft_printf("-->nombre de proprietes = %i\n", vars->nprop);
  for (i=0 ; i<vars->nprop ; i++)
    bft_printf("-->properties_ipp[%i]: %i propce[%i]: %i "
               "properties_name[%i]: %s\n",
               i, vars->properties_ipp[i],
               i, vars->propce[i],
               i, vars->properties_name[i]);
#endif
}

/*------------------------------------------------------------------------------
 * Indirection between the solver numbering and the XML one
 * for the model scalar (electrical model)
 *----------------------------------------------------------------------------*/
void CS_PROCF (uielsc, UIELSC) (const int *const ippmod,
                                const int *const ieljou,
                                const int *const ielarc,
                                const int *const ngazg,
                                const int *const ihm,
                                const int *const ipotr,
                                const int *const iycoel,
                                const int *const ipoti,
                                const int *const ipotva)
{
  cs_var_t  *vars = cs_glob_var;
  char *snumsca = NULL;
  char *name = NULL;
  char *label = NULL;

  if (vars->nscaus > 0)
    BFT_REALLOC(vars->label, vars->nscapp + vars->nscaus, char*);
  else
    BFT_MALLOC(vars->label, vars->nscapp, char*);

  BFT_MALLOC(snumsca, 1 + 1, char);

  label = _scalar_name_label("joule_effect", "Enthalpy");
  BFT_MALLOC(vars->label[*ihm -1], strlen(label)+1, char);
  strcpy(vars->label[*ihm -1], label);
  BFT_FREE(label);

  label = _scalar_name_label("joule_effect", "PotElecReal");
  BFT_MALLOC(vars->label[*ipotr -1], strlen(label)+1, char);
  strcpy(vars->label[*ipotr -1], label);
  BFT_FREE(label);

  if (*ngazg > 1)
    for (int iesp = 0; iesp < *ngazg - 1; iesp++) {
      BFT_MALLOC(name, strlen("YM_ESL") +1 + 1, char);
      strcpy(name, "YM_ESL");
      sprintf(snumsca, "%2.2i", iesp +1);
      strcat(name, snumsca);
      label = _scalar_name_label("joule_effect", name);
      BFT_MALLOC(vars->label[iycoel[iesp] -1], strlen(label) +1, char);
      strcpy(vars->label[iycoel[iesp] -1], label);
      BFT_FREE(label);
    }

  if (ippmod[*ieljou - 1] == 2 || ippmod[*ieljou - 1] == 4) {
    label = _scalar_name_label("joule_effect", "POT_EL_I");
    BFT_MALLOC(vars->label[*ipoti -1], strlen(label) +1, char);
    strcpy(vars->label[*ipoti -1], label);
    BFT_FREE(label);
  }

  if (ippmod[*ielarc - 1] >= 2)
    for (int idimve = 0; idimve < 3; idimve++) {
      BFT_MALLOC(name, strlen("POT_VEC") +1 + 1, char);
      strcpy(name, "POT_VEC");
      sprintf(snumsca, "%2.2i", idimve +1);
      strcat(name, snumsca);
      label = _scalar_name_label("joule_effect", name);
      BFT_MALLOC(vars->label[ipotva[idimve] -1], strlen(label) +1, char);
      strcpy(vars->label[ipotva[idimve] -1], label);
      BFT_FREE(label);
    }

  BFT_FREE(snumsca);
  BFT_FREE(name);

#if _XML_DEBUG_
  bft_printf("==>UICPSC\n");
  for (i=0; i< vars->nscaus+vars->nscapp; i++)
    bft_printf("--label of scalar[%i]: %s\n", i, vars->label[i]);
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
 * Atmospheric flows: indirection between the solver numbering and the XML one
 * for physical properties.
 *
 * Fortran Interface:
 *
 * subroutine uiatpr
 * *****************
 * integer         nsalpp   -->
 * integer         ippmod   -->   specific physics indicator array
 * integer         iatmos   -->   index for atmospheric flow
 * integer         ipppro   -->
 * integer         ipproc   -->
 * integer         itempc   -->   index for real temperature
 * integer         iliqwt   -->   index for liquid water
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiatpr, UIATPR) (const int *const nsalpp,
                                const int *const ippmod,
                                const int *const iatmos,
                                const int *const ipppro,
                                const int *const ipproc,
                                const int *const itempc,
                                const int *const iliqwt)
{
  int n;
  cs_var_t *vars = cs_glob_var;

  n = vars->nprop;
  vars->nprop  += *nsalpp;
  vars->nsalpp  = *nsalpp;

  BFT_REALLOC(vars->properties_ipp,  vars->nprop, int);
  BFT_REALLOC(vars->propce,          vars->nprop, int);
  BFT_REALLOC(vars->properties_name, vars->nprop, char*);

  /* itempc */
  vars->properties_ipp[n] = ipppro[ipproc[ *itempc -1] -1];
  vars->propce[n] = ipproc[ *itempc -1] -1;
  BFT_MALLOC(vars->properties_name[n], strlen("real_temperature") +1, char);
  strcpy(vars->properties_name[n++], "real_temperature");

  if (ippmod[*iatmos -1] == 2) {
    /* iliqwt */
    vars->properties_ipp[n] = ipppro[ipproc[ *iliqwt -1] -1];
    vars->propce[n] = ipproc[ *iliqwt -1] -1;
    BFT_MALLOC(vars->properties_name[n], strlen("liquid_water") +1, char);
    strcpy(vars->properties_name[n++], "liquid_water");
  }
#if _XML_DEBUG_
  {
    int i;
    bft_printf("==>UIATPR\n");
    bft_printf("-->nombre de proprietes = %i\n", vars->nprop);
    for (i=0 ; i<vars->nprop ; i++)
      bft_printf("-->properties_ipp[%i]: %i propce[%i]: %i "
                 "properties_name[%i]: %s\n",
                 i, vars->properties_ipp[i],
                 i, vars->propce[i],
                 i, vars->properties_name[i]);
  }
#endif
}

/*----------------------------------------------------------------------------
 * Atmospheric flows: indirection between the solver numbering and the XML one
 * for models scalars.
 *
 * Fortran Interface:
 *
 * subroutine uiatsc
 * *****************
 * integer         ippmod   -->   specific physics indicator array
 * integer         iatmos   -->   index for atmospheric flow
 * integer         itempp   -->   index for potential temperature
 * integer         itempl   -->   index for liquid potential temperature
 * integer         itotwt   -->   index for total water content
 * integer         intdrp   -->   index for total number of droplets
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiatsc, UIATSC) (const int *const ippmod,
                                const int *const iatmos,
                                const int *const itempp,
                                const int *const itempl,
                                const int *const itotwt,
                                const int *const intdrp)
{
  cs_var_t  *vars = cs_glob_var;
  char *label = NULL;

  if (vars->nscaus > 0)
    BFT_REALLOC(vars->label, vars->nscapp + vars->nscaus, char*);
  else
    BFT_MALLOC(vars->label, vars->nscapp, char*);

  if (ippmod[*iatmos -1] == 1)
  {
    /* itempp */
    label = _scalar_name_label("atmospheric_flows", "potential_temperature");
    BFT_MALLOC(vars->label[*itempp -1], strlen(label)+1, char);
    strcpy(vars->label[*itempp -1], label);
    BFT_FREE(label);
  }
  else if (ippmod[*iatmos -1] == 2)
  {
    /* itempl */
    label = _scalar_name_label("atmospheric_flows", "liquid_potential_temperature");
    BFT_MALLOC(vars->label[*itempl -1], strlen(label)+1, char);
    strcpy(vars->label[*itempl -1], label);
    BFT_FREE(label);

    /* itotwt */
    label = _scalar_name_label("atmospheric_flows", "total_water");
    BFT_MALLOC(vars->label[*itotwt -1], strlen(label)+1, char);
    strcpy(vars->label[*itotwt -1], label);
    BFT_FREE(label);

    /* intdrp */
    label = _scalar_name_label("atmospheric_flows", "number_of_droplets");
    BFT_MALLOC(vars->label[*intdrp -1], strlen(label)+1, char);
    strcpy(vars->label[*intdrp -1], label);
    BFT_FREE(label);
  }
#if _XML_DEBUG_
  {
    int i;
    bft_printf("==>UIATSC\n");
    for (i=0; i< vars->nscaus+vars->nscapp; i++)
      bft_printf("--label of scalar[%i]: %s\n", i, vars->label[i]);
  }
#endif
}


/*----------------------------------------------------------------------------
 * Indirection between the solver numbering and the XML one
 * for physical properties of the activated specific physics
 * (pulverized solid fuels)
 *----------------------------------------------------------------------------*/

void CS_PROCF (uisofu, UISOFU) (const int    *const ippmod,
                                const int    *const iccoal,
                                const int    *const icpl3c,
                                const int    *const iirayo,
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
                                double       *const ipci,
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
  {
    cs_gui_load_file("dp_FCP.xml");

    if (ippmod[*iccoal - 1] == 0)
    {
      BFT_MALLOC(vars->model_value, strlen("homogeneous_fuel")+1, char);
      strcpy(vars->model, "homogeneous_fuel");
    } else if (ippmod[*iccoal - 1] == 1) {
      if (ippmod[*icpl3c - 1] > 0) {
        BFT_MALLOC(vars->model_value, strlen("homogeneous_fuel_moisture_lagr")+1, char);
        strcpy(vars->model, "homogeneous_fuel_moisture_lagr");
      } else {
        BFT_MALLOC(vars->model_value, strlen("homogeneous_fuel_moisture")+1, char);
        strcpy(vars->model, "homogeneous_fuel_moisture");
      }
    }
  }

  /* ---- Nb de charbons */
  *ncharb = cs_gui_get_tag_number("/solid_fuels/solid_fuel", 1);
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
    if (cs_gui_strcmp(vars->model_value, "homogeneous_fuel_moisture_lagr"))
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
  *noxyd = cs_gui_get_tag_number("/oxidants/oxidant", 1);
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

/*============================================================================
 * Public function definitions
 *============================================================================*/

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

  cs_var_t  *vars = cs_glob_var;

  const char *name[] = { "solid_fuels",
                         "gas_combustion",
                         "joule_effect",
                         "atmospheric_flows",
                         "compressible_model" };
  int name_nbr = sizeof(name) / sizeof(name[0]);

  if (vars->model != NULL && vars->model_value != NULL) {
    isactiv = 1;
    return isactiv;
  } else {
    vars->model = NULL;
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

/*----------------------------------------------------------------------------*/

END_C_DECLS
