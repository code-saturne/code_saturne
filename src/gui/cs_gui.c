/*============================================================================
 * Management of the GUI parameters file: main parameters
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
 * libxml2 library headers
 *----------------------------------------------------------------------------*/

#if defined(HAVE_LIBXML2)

#include <libxml/tree.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>

#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_selector.h"

#include "mei_evaluate.h"
#include "mei_math_util.h"

#include "cs_base.h"
#include "cs_file.h"
#include "cs_gui_util.h"
#include "cs_gui_variables.h"
#include "cs_gui_boundary_conditions.h"
#include "cs_gui_specific_physics.h"
#include "cs_gui_mobile_mesh.h"
#include "cs_mesh.h"
#include "cs_partition.h"
#include "cs_prototypes.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gui.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* debugging switch */
#define _XML_DEBUG_ 0

/*============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 * External global variables
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Management of the XML document
 *----------------------------------------------------------------------------*/

#if defined(HAVE_LIBXML2)
extern xmlXPathContextPtr xpathCtx;   /* Pointer on the Context       */
extern xmlNodePtr node;               /* Pointer on the root node     */
#endif

/*============================================================================
 * Private global variables
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer on the main variable structure */

cs_var_t    *cs_glob_var = NULL;
cs_label_t  *cs_glob_label = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Turbulence model parameters.
 *
 * parameters:
 *   param                -->  name of the parameters
 *   keyword             <--   turbulence model parameter
 *----------------------------------------------------------------------------*/

static void
cs_gui_advanced_options_turbulence(const char *const param,
                                         int  *const keyword)
{
  char *path = NULL;
  int  result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "thermophysical_models", "turbulence", param);

  if (cs_gui_strcmp("gravity_terms", param)) {
    cs_xpath_add_attribute(&path, "status");
    if (cs_gui_get_status(path, &result)) *keyword = result;

  } else if (cs_gui_strcmp("scale_model", param)) {
    cs_xpath_add_function_text(&path);
    if (cs_gui_get_int(path, &result)) *keyword = result;

  } else
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Return the name of the related scalar if the scalar "num_sca" is a variance
 *
 * parameter:
 *   num_sca           -->  scalar number
 *----------------------------------------------------------------------------*/

static char *
cs_gui_scalar_variance(const int num_sca)
{
  char *path = NULL;
  char *variance = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "additional_scalars");
  cs_xpath_add_element_num(&path, "scalar", num_sca);
  cs_xpath_add_element(&path, "variance");
  cs_xpath_add_function_text(&path);

  variance = cs_gui_get_text_value(path);

  BFT_FREE(path);

  return variance;
}

/*-----------------------------------------------------------------------------
 * Return the user thermal scalar indicator.
 *----------------------------------------------------------------------------*/

static int
cs_gui_thermal_scalar(void)
{
  char *model_name = NULL;
  int   test = 0;

  model_name = cs_gui_get_thermophysical_model("thermal_scalar");

  if (cs_gui_strcmp(model_name, "off"))
    test = 0;
  else {
    if (cs_gui_strcmp(model_name, "enthalpy"))
      test =  2 ;
    else if (cs_gui_strcmp(model_name, "temperature_kelvin"))
      test =  1 ;
    else if (cs_gui_strcmp(model_name, "temperature_celsius"))
      test = -1 ;
    else
      bft_error(__FILE__, __LINE__, 0,
          _("Invalid thermal model: %s\n"), model_name);
  }

  BFT_FREE(model_name);

  return test;
}

/*----------------------------------------------------------------------------
 * Get thermal user scalar number if it is exist.
 *
 * parameters:
 *   iscalt               <--  thermal scalar number order
 *   iscsth               <--  nature of the thermal scalar (C, K, J/kg)
 *----------------------------------------------------------------------------*/

static void
cs_gui_thermal_scalar_number(int *const iscalt,
                             int *const iscsth)
{
  int ind_thermal;
  int i, index, size;
  char *path = NULL;
  char **name = NULL;

  ind_thermal = cs_gui_thermal_scalar();

  if (ind_thermal) {
    path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 2, "additional_scalars", "/@type");
    name = cs_gui_get_attribute_values(path, &size);

    index = -1;
    for (i=0; i < size; i++)
      if (cs_gui_strcmp(name[i], "thermal"))
        index = i;

    *iscalt = index+1;
    iscsth[index] = ind_thermal;

    BFT_FREE(path);
    for (i=0; i < size; i++) BFT_FREE(name[i]);
    BFT_FREE(name);
  }
}

/*-----------------------------------------------------------------------------
 * Return the name of the diffusion_coefficient property for a scalar
 *
 * parameters:
 *   scalar_index   --> index of the scalar
 *----------------------------------------------------------------------------*/

static char *
_scalar_diffusion_coefficient_name(const int idx)
{
  int ncar = 0;
  char *name = NULL;
  char *suf = NULL;

  ncar = cs_gui_characters_number(idx+1);
  BFT_MALLOC(name, strlen("diffusion_coefficient") +2 +ncar, char);
  BFT_MALLOC(suf, 1 + ncar, char);
  sprintf(suf, "%i", idx+1);
  strcpy(name, "diffusion_coefficient");
  strcat(name, "_");
  strcat(name, suf);
  BFT_FREE(suf);
  return name;
}

/*-----------------------------------------------------------------------------
 * Return the value of choice for user scalar's property
 *
 * parameters:
 *   scalar_num     --> number of scalar
 *   choice         --> choice for property
 *----------------------------------------------------------------------------*/

static int
cs_gui_scalar_properties_choice(const int scalar_num, int *const choice)
{
  char *path = NULL;
  char *buff = NULL;
  int   ichoice;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "additional_scalars");
  cs_xpath_add_element_num(&path, "scalar", scalar_num);
  cs_xpath_add_element(&path, "property");
  cs_xpath_add_attribute(&path, "choice");

  buff = cs_gui_get_attribute_value(path);

  if (buff == NULL) {
    ichoice = 0;

  } else {
    ichoice = 1;

    if (cs_gui_strcmp(buff, "variable") || cs_gui_strcmp(buff, "user_law"))
      *choice = 1;
    else if (cs_gui_strcmp(buff, "constant"))
      *choice = 0;
    else
      bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);
  }

  BFT_FREE(path);
  BFT_FREE(buff);

  return ichoice;
}

/*-----------------------------------------------------------------------------
 * Return value of diffusion coefficient for user scalars
 *        return 1 if value exists
 *        return 0 if not
 *
 * parameters:
 *   num_sca  --> number of scalar
 *   value   <--  value of diffusion coefficient
 *----------------------------------------------------------------------------*/

static void
cs_gui_scalar_diffusion_value(const int           num_sca,
                                    double *const value)
{
  char  *path = NULL;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "additional_scalars");
  cs_xpath_add_element_num(&path, "scalar", num_sca);
  cs_xpath_add_element(&path, "property");
  cs_xpath_add_element(&path, "initial_value");
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &result))
    *value = result;

  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Get the status of steady management.
 *
 * parameter:
 *   keyword         <--  if 1 unsteady management else steady management
 *----------------------------------------------------------------------------*/

static void
cs_gui_get_steady_status(int *const keyword)
{
  char *path = NULL;
  int   result;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "steady_management");
  cs_xpath_add_attribute(&path, "status");

  if (cs_gui_get_status(path, &result))
    *keyword = result;
  else
    *keyword = 1;

  BFT_FREE(path);
}

/*----------------------------------------------------------------------------
 * Return the initialization choice of the turbulence variables.
 *----------------------------------------------------------------------------*/

static char *cs_gui_velocity_pressure_algo_choice(void)
{
  char *path = NULL;
  char *algo_choice;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,
                        "numerical_parameters",
                        "velocity_pressure_algo");
  cs_xpath_add_attribute(&path, "choice");

  algo_choice = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return algo_choice;
}

/*-----------------------------------------------------------------------------
 * Return  parameters for steady management.
 *
 * parameter:
 *   param           -->  steady parameter
 *   keyword         <--  new value for the steady parameter
 *----------------------------------------------------------------------------*/

static void
cs_gui_steady_parameters(const char   *const param,
                               double *const keyword)
{
  char   *path   = NULL;
  double  result = 0.0;
  int     status = 0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "analysis_control", "steady_management", param);

  if (cs_gui_strcmp(param,"zero_iteration")) {
    cs_xpath_add_attribute(&path, "status");
    if(cs_gui_get_status(path, &status))
      *keyword = status;
  } else {
    cs_xpath_add_function_text(&path);
    if (cs_gui_get_double(path, &result))
      *keyword = result;
  }
  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Modify time parameters.
 *
 * parameters:
 *   param              -->  time parameter
 *   keyword            <--  new value of the time parameter
 *----------------------------------------------------------------------------*/

static void
cs_gui_time_parameters(const char   *const param,
                             double *const keyword)
{
  char   *path   = NULL;
  double  result = 0.0;
  int     status = 0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "analysis_control", "time_parameters", param);

  if (cs_gui_strcmp(param,"zero_time_step") ||
      cs_gui_strcmp(param,"thermal_time_step")) {

    cs_xpath_add_attribute(&path, "status");
    if(cs_gui_get_status(path, &status))
      *keyword = status;

  } else {
    cs_xpath_add_function_text(&path);
    if (cs_gui_get_double(path, &result))
      *keyword = result;
  }
  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Modify restart parameters.
 *
 * parameters:
 *   param               -->  restart parameter
 *   keyword            <-->  new value of the restart parameter
 *----------------------------------------------------------------------------*/

static void
cs_gui_restart_parameters_status(const char *param, int *const keyword)
{
  int   result;
  char *path = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "calculation_management", "start_restart", param);

  if (cs_gui_strcmp(param, "restart_rescue")) {
    cs_xpath_add_function_text(&path);
    if (cs_gui_get_int(path, &result))
      *keyword = result;

  } else {
    cs_xpath_add_attribute(&path, "status");

    if (cs_gui_get_status(path, &result))
      *keyword = result;
  }

  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Return value of numerical parameter markup
 *
 * parameters:
 *   variable_type  --> name of variable
 *   value_type     --> name of numerical parameter parkup
 *   value          <-- value of numerical parameter
 *----------------------------------------------------------------------------*/

static void
cs_gui_variable_value(const char   *const variable_type,
                      const char   *const value_type,
                            double *const value)
{
  char  *path = NULL;
  double result;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "variable");
  cs_xpath_add_test_attribute(&path, "name", variable_type);
  cs_xpath_add_element(&path, value_type);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &result))
    *value = result;

  BFT_FREE(path);
}

/*----------------------------------------------------------------------------
 * Get the attribute value from the xpath query.
 *
 * parameters:
 *   path          --> path for xpath query
 *   child         --> child markup
 *   keyword      <--  value of attribute node
 *----------------------------------------------------------------------------*/

static void
_attribute_value(      char *      path,
                 const char *const child,
                       int  *const keyword)
{
  char *choice = NULL;
  int   result;

  assert(path != NULL);
  assert(child != NULL);

  if (cs_gui_strcmp(child, "order_scheme")) {

    /* *keyword = 1; */
    cs_xpath_add_attribute(&path, "choice");
    choice = cs_gui_get_attribute_value(path);

    if (cs_gui_strcmp(choice, "centered"))
      *keyword = 1;
    else if (cs_gui_strcmp(choice, "solu"))
      *keyword = 0;
    BFT_FREE(choice);

  } else {

    cs_xpath_add_attribute(&path, "status");

    if (cs_gui_get_status(path, &result)) {
      *keyword = result;

      if (cs_gui_strcmp(child, "slope_test")) {
        if (result == 1)
          *keyword = 0;
        if (result == 0)
          *keyword = 1;
      }
    }
  }
  BFT_FREE(path);
}

/*----------------------------------------------------------------------------
 * Get the attribute value associated to a child markup from a variable.
 *
 * parameters:
 *   name          -->  name of the variable markup
 *   child         -->  child markup
 *   keyword      <--   value of attribute node contained in the child markup
 *----------------------------------------------------------------------------*/

static void
cs_gui_variable_attribute(const char *const name,
                          const char *const child,
                                int  *const keyword)
{
  char *path = NULL;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "variable");
  cs_xpath_add_test_attribute(&path, "name", name);
  cs_xpath_add_element(&path, child);

  _attribute_value(path, child, keyword);
}

/*----------------------------------------------------------------------------
 * Return the attribute choice associated to a child markup from a variable.
 *
 * parameters:
 *   name          -->  name of the variable markup
 *   child         -->  child markup
 *----------------------------------------------------------------------------*/

static char *cs_gui_variable_choice(const char *const name,
                                    const char *const child)
{
  char *path = NULL;
  char *choice;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "variable");
  cs_xpath_add_test_attribute(&path, "name", name);
  cs_xpath_add_element(&path, child);
  cs_xpath_add_attribute(&path, "choice");

  choice = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return choice;
}

/*----------------------------------------------------------------------------
 * Return the attribute choice associated to a child markup from a variable.
 *
 * parameters:
 *   name          -->  name of the variable markup
 *   child         -->  child markup
 *----------------------------------------------------------------------------*/

static char *cs_gui_model_variable_choice(const char *const name,
                                          const char *const child)
{
  char *path = NULL;
  char *choice;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "scalar");
  cs_xpath_add_test_attribute(&path, "label", name);
  cs_xpath_add_element(&path, child);
  cs_xpath_add_attribute(&path, "choice");

  choice = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return choice;
}

/*----------------------------------------------------------------------------
 * Get the text value associated to a child markup from a scalar.
 *
 * parameters:
 *   label                -->  label of the scalar markup
 *   child                -->  name of the child markup
 *   value               <--   value of text node contained in the child markup
 *----------------------------------------------------------------------------*/

static void
cs_gui_scalar_value(const char   *const label,
                    const char   *const child,
                          double *const value)
{
  char   *path = NULL;
  double  result;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "scalar");
  cs_xpath_add_test_attribute(&path, "label", label);
  cs_xpath_add_element(&path, child);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &result))
    *value = result;

  BFT_FREE(path);
}

/*----------------------------------------------------------------------------
 * Get the attribute value associated to a child markup from a scalar.
 *
 * parameters:
 *   label         -->  name of the scalar markup
 *   child         -->  child markup
 *   keyword      <--   value of attribute node contained in the child markup
 *----------------------------------------------------------------------------*/

static void
cs_gui_scalar_attribute(const char *const label,
                        const char *const child,
                              int  *const keyword)
{
  char *path = NULL;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "scalar");
  cs_xpath_add_test_attribute(&path, "label", label);
  cs_xpath_add_element(&path, child);

  _attribute_value(path, child, keyword);
}

/*-----------------------------------------------------------------------------
 * Get values related the modelling scalar: min, max ...
 *
 * returns:
 *    1 if value exists
 *    0 otherwise
 *    the result is stored in "value"
 *----------------------------------------------------------------------------*/


static void
cs_gui_model_scalar_value(const   char *const model,
                          const   char *const name,
                          const   char *const keyword,
                                double *const value)
{
  char   *path = NULL;
  double  result;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "thermophysical_models");
  cs_xpath_add_element(&path, model);
  cs_xpath_add_element(&path, "scalar");
  cs_xpath_add_test_attribute(&path, "label", name);
  cs_xpath_add_element(&path, keyword);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path,&result))
    *value = result;

  BFT_FREE(path);
}

/*----------------------------------------------------------------------------
 * Get the attribute value associated to a child markup from a scalar.
 *
 * parameters:
 *   model         -->  model markup
 *   name          -->  name of the scalar markup
 *   child         -->  child markup
 *   keyword      <--   value of attribute node contained in the child markup
 *----------------------------------------------------------------------------*/

static void
cs_gui_model_scalar_output_status(const char *const model,
                                  const char *const name,
                                  const char *const child,
                                        int  *const keyword)
{
  char *path = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "thermophysical_models");
  cs_xpath_add_element(&path, model);
  cs_xpath_add_element(&path, "scalar");
  cs_xpath_add_test_attribute(&path, "label", name);
  cs_xpath_add_element(&path, child);

  _attribute_value(path, child, keyword);
}

/*-----------------------------------------------------------------------------
 * Modify double numerical parameters.
 *
 * parameters:
 *   param               -->  label of the numerical parameter
 *   keyword            <-->  value of the numerical parameter
 *----------------------------------------------------------------------------*/

void
cs_gui_numerical_double_parameters(const char   *const param,
                                         double *const keyword)
{
  char  *path = NULL;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "numerical_parameters");
  cs_xpath_add_element(&path, param);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &result))
    *keyword = result;

  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Modify integer numerical parameters.
 *
 * parameters:
 *   param               -->  label of the numerical parameter
 *   keyword            <-->  value of the numerical parameter
 *----------------------------------------------------------------------------*/

static void
cs_gui_numerical_int_parameters(const char *const param,
                                      int  *const keyword)
{
  char *path = NULL;
  char *choice = NULL;
  int   result;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "numerical_parameters");

  if (cs_gui_strcmp(param, "gradient_reconstruction")) {

    cs_xpath_add_element(&path, param);
    cs_xpath_add_attribute(&path, "choice");
    choice = cs_gui_get_attribute_value(path);
    if (choice)
      *keyword = atoi(choice);
    BFT_FREE(choice);

  } else if (cs_gui_strcmp(param,"piso_sweep_number")) {

    cs_xpath_add_element(&path, "velocity_pressure_algo");
    cs_xpath_add_element(&path, param);
    cs_xpath_add_function_text(&path);
    if (cs_gui_get_int(path, &result))
      *keyword = result;

  } else {

    cs_xpath_add_element(&path, param);
    cs_xpath_add_attribute(&path, "status");
    if (cs_gui_get_status(path, &result))
      *keyword = result;
  }
  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Modify gravity parameters.
 *
 * parameters:
 *   param               -->  gravity parameter (GX, GY, GZ)
 *   keyword            <-->  new value of the gravity parameter
 *----------------------------------------------------------------------------*/

static void
cs_gui_gravity_value(const char   *const param,
                           double *const value)
{
  char   *path = NULL;
  double  result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "physical_properties", "gravity", param);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &result))
    *value = result;

  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Modify coriolis source terms parameters.
 *
 * parameters:
 *   param               -->  coriolis parameter (OMEGAX, OMEGAY, OMEGAZ)
 *   keyword            <-->  new value of the coriolis parameter
 *----------------------------------------------------------------------------*/

static void
cs_gui_coriolis_value(const char   *const param,
                            double *const value)
{
  char   *path = NULL;
  double  result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "physical_properties", "omega", param);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &result))
    *value = result;

  BFT_FREE(path);
}

/*----------------------------------------------------------------------------
 * Return the value of the choice attribute from a property name.
 *
 * parameters:
 *   property_name        -->  name of the property
 *----------------------------------------------------------------------------*/

static char*
_properties_choice(const char *const property_name)
{
  char *path   = NULL;
  char *choice = NULL;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "property");
  cs_xpath_add_test_attribute(&path, "name", property_name);
  cs_xpath_add_attribute(&path, "choice");
  choice = cs_gui_get_attribute_value(path);
  BFT_FREE(path);
  return choice;
}

/*----------------------------------------------------------------------------
 * Get the value of the choice attribute from a property markup.
 * Return 1 if the xpath request has succeeded, 0 otherwise.
 *
 * parameters:
 *   property_name        -->  name of the property
 *   choice              <--   value of the attribute choice
 *----------------------------------------------------------------------------*/

static int
cs_gui_properties_choice(const char *const property_name, int *choice)
{
  char *buff = NULL;
  int   iok = 0;

  buff = _properties_choice(property_name);
  if (buff)
  {
    iok = 1;
    if (cs_gui_strcmp(buff, "variable") || cs_gui_strcmp(buff, "user_law"))
      *choice = 1;
    else if (cs_gui_strcmp(buff, "constant"))
      *choice = 0;
  }
  else
    iok = 0;
  BFT_FREE(buff);
  return iok;
}

/*----------------------------------------------------------------------------
 * Turbulence model parameters.
 *
 * parameters:
 *   param                -->  name of the parameters
 *   keyword             <--   turbulence model parameter
 *----------------------------------------------------------------------------*/

static void _option_turbulence_double(const char *const param,
                                          double *const keyword)
{
  char *path = NULL;
  double  result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "thermophysical_models", "turbulence", param);

  cs_xpath_add_function_text(&path);
  if (cs_gui_get_double(path, &result))
    *keyword = result;

  BFT_FREE(path);
}

/*----------------------------------------------------------------------------
 * Return the length choice for initialize turbulence
 *----------------------------------------------------------------------------*/

static char *cs_gui_reference_length_initialization_choice(void)
{
  char *path = NULL;
  char *initialization_choice;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3,
                        "thermophysical_models",
                        "reference_values",
                        "length");
  cs_xpath_add_attribute(&path, "choice");

  initialization_choice = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return initialization_choice;
}

/*----------------------------------------------------------------------------
 * Return the initialization choice of the turbulence variables.
 *
 * parameters:
 *   zone_id        -->  zone number
 *----------------------------------------------------------------------------*/

static char *cs_gui_turbulence_initialization_choice(const char* zone_id)
{
  char *path = NULL;
  char *initialization_choice;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3,
                        "thermophysical_models",
                        "turbulence",
                        "initialization");
  cs_xpath_add_test_attribute(&path, "zone_id", zone_id);
  cs_xpath_add_function_text(&path);
  cs_xpath_add_attribute(&path, "choice");

  initialization_choice = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return initialization_choice;
}

/*==================================
 * TREATMENTS FOR TIME AVERAGES
 *=================================*/

/*----------------------------------------------------------------------------
 * Return the number of variables and properties inside a given time average.
 *
 * parameters:
 *   id           -->  time average number (imom)
 *----------------------------------------------------------------------------*/

static int
_get_time_average_n_variables(const int id)
{
  char *path = NULL;
  int   number = 0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "analysis_control", "time_averages");
  cs_xpath_add_element_num(&path, "time_average", id);
  cs_xpath_add_element(&path, "var_prop");
  number = cs_gui_get_nb_element(path);

  BFT_FREE(path);

  return number;

}

/*----------------------------------------------------------------------------
 * Get value of a parameter for a given time avegare.
 *
 * parameters:
 *   id              -->  time average number (imom)
 *   param           -->  name of the parameter
 *   data           <--   value of the parameter
 *----------------------------------------------------------------------------*/

static void _get_time_average_data(const int         id,
                                   const char *const param,
                                         int  *const data)
{
  char *path = NULL;
  int   result = 0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "analysis_control", "time_averages");
  cs_xpath_add_element_num(&path, "time_average", id);
  cs_xpath_add_element(&path, param);

  cs_xpath_add_function_text(&path);
  if (cs_gui_get_int(path, &result))
    *data = result;

  BFT_FREE(path);
}

/*----------------------------------------------------------------------------
 * Return the name of a variable or a property for a given time average.
 *
 * parameters:
 *   id           -->  time average number (imom)
 *   nb           -->  variable or property number
 *----------------------------------------------------------------------------*/

static char *_get_time_average_variable_name(const int id, const int nb)
{
  char *path = NULL;
  char *name = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "analysis_control", "time_averages");
  cs_xpath_add_element_num(&path, "time_average", id);
  cs_xpath_add_element_num(&path, "var_prop", nb);
  cs_xpath_add_attribute(&path, "name");

  name = cs_gui_get_attribute_value(path);
  BFT_FREE(path);

  return name;
}

/*----------------------------------------------------------------------------
 * Return the label of a time average.
 *
 * parameters:
 *   id              -->  time average number (imom)
 *----------------------------------------------------------------------------*/

static char *_get_time_average_label(const int id)
{
  char *path = NULL;
  char *label = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "analysis_control", "time_averages");
  cs_xpath_add_element_num(&path, "time_average", id);
  cs_xpath_add_attribute(&path,"label");

  label = cs_gui_get_attribute_value(path);
  BFT_FREE(path);

  return label;
}

/*-----------------------------------------------------------------------------
 * Return label of variable
 *
 * parameters:
 *   variable   --> name of variable
 *----------------------------------------------------------------------------*/

static char *cs_gui_variable_label (const char *const variable)
{
  char *path = NULL;
  char *label = NULL;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "variable");
  cs_xpath_add_test_attribute(&path, "name", variable);
  cs_xpath_add_attribute(&path, "label");

  label = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return label;
}

/*-----------------------------------------------------------------------------
 * Return the label attribute of a property markup.
 *
 * parameters:
 *   property_name        -->  name of the property
 *----------------------------------------------------------------------------*/

static char *cs_gui_properties_label(const char *const property_name)
{
  char *path = NULL;
  char *label = NULL;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "property");
  cs_xpath_add_test_attribute(&path, "name", property_name);
  cs_xpath_add_attribute(&path, "label");

  label = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return label;
}

/*-----------------------------------------------------------------------------
 * Return the label or the name from a scalar.
 *
 * parameters:
 *   kw                   -->  keyword: 'label' or 'name'
 *   scalar_num          <--   number of the searching scalar
 *----------------------------------------------------------------------------*/

static char *_scalar_name_label(const char *kw, const int scalar_num)
{
  char *path = NULL;
  char *str  = NULL;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "additional_scalars");
  cs_xpath_add_element_num(&path, "scalar", scalar_num);
  cs_xpath_add_attribute(&path, kw);

  str = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return str;
}

/*-----------------------------------------------------------------------------
 * Return the name from a specific physic scalar.
 *
 * parameters:
 *   physics              -->  keyword: specific physic model required
 *   kw                   -->  scalar name
 *----------------------------------------------------------------------------*/

static char *_specific_physic_scalar_name_label(const char *physics, const char *kw)
{
  char *path = NULL;
  char *str  = NULL;

  path = cs_xpath_short_path();
  cs_xpath_add_elements(&path, 3,
                        "thermophysical_models",
                        physics,
                        "scalar");
  cs_xpath_add_test_attribute(&path, "label", kw);
  cs_xpath_add_attribute(&path, "name");

  str = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return str;
}

/*==========================
 * FOR VOLUMICS ZONES
 *==========================*/

/*-----------------------------------------------------------------------------
 * Return the name of the volumic zone
 *
 * parameters:
 *   ith_zone        -->  id of volumic zone
 *----------------------------------------------------------------------------*/

static char *cs_gui_volumic_zone_id(const int ith_zone)
{
  char *path = NULL;
  char *name = NULL;

  /* 1) get the name of the ith initialization zone */
  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "solution_domain", "volumic_conditions");
  cs_xpath_add_element_num(&path, "zone", ith_zone);
  cs_xpath_add_attribute(&path, "id");

  name = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return name;
}

/*-----------------------------------------------------------------------------
 * Return the localisation for the volumic zone with a given id
 *
 * parameters:
 *   zone_id      -->  volumic zone id
 *----------------------------------------------------------------------------*/

static char *cs_gui_volumic_zone_localization(const char *const zone_id)
{
  char *path = NULL;
  char *description = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "solution_domain",
                        "volumic_conditions",
                        "zone");
  cs_xpath_add_test_attribute(&path, "id", zone_id);
  cs_xpath_add_function_text(&path);

  description = cs_gui_get_text_value(path);

  BFT_FREE(path);

  return description;
}

/*-----------------------------------------------------------------------------
 * Return the initial value of variable for the volumic zone named name
 *
 * parameters:
 *   variable_name    -->  name of variable
 *   zone_id          -->  id of volumic zone
 *   initial_value    <--  initial value
 *----------------------------------------------------------------------------*/

#if (_XML_DEBUG_ > 0)

static void cs_gui_variable_initial_value(const char   *const variable_name,
                                          const char   *const zone_id,
                                                double *const initial_value)
{
  char *path = NULL;
  double result;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "variable");
  cs_xpath_add_test_attribute(&path, "name", variable_name);
  cs_xpath_add_element(&path, "initial_value");
  cs_xpath_add_test_attribute(&path, "zone_id", zone_id);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &result))
    *initial_value = result;
  else
    *initial_value = 0.0;

  BFT_FREE(path);
}

#endif /* (_XML_DEBUG_ > 0) */

/*----------------------------------------------------------------------------
 * Get label of 1D profile file name
 *
 * parameters:
 *   id           -->  number of order in list of 1D profile
 *----------------------------------------------------------------------------*/

static char *_get_profile(const char *kw, const int id)
{
  char *path = NULL;
  char *label = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "analysis_control", "profiles");
  cs_xpath_add_element_num(&path, "profile", id+1);
  cs_xpath_add_attribute(&path, kw);

  label = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return label;
}

/*----------------------------------------------------------------------------
 * Get number of variables or properties or scalar for 1D profile
 *
 * parameters:
 *   id           -->  number of 1D profile
 *----------------------------------------------------------------------------*/

static int _get_profile_names_number(const int id)
{
  char *path = NULL;
  int   number = 0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "analysis_control", "profiles");
  cs_xpath_add_element_num(&path, "profile", id+1);
  cs_xpath_add_element(&path, "var_prop");
  number = cs_gui_get_nb_element(path);

  BFT_FREE(path);

  return number;
}

/*----------------------------------------------------------------------------
 * Return the name of variables or properties or scalar for 1D profile
 *
 * parameters:
 *   id           -->  number of 1D profile
 *   nm           -->  number of the variable name of the idst 1D profile
 *----------------------------------------------------------------------------*/

static char *_get_profile_name(const int id, const int nm)
{
  char *path = NULL;
  char *name = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "analysis_control", "profiles");
  cs_xpath_add_element_num(&path, "profile", id+1);
  cs_xpath_add_element_num(&path, "var_prop", nm+1);
  cs_xpath_add_attribute(&path, "name");

  name = cs_gui_get_attribute_value(path);
  if (name == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Invalid xpath: %s\n name not found"), path);
  BFT_FREE(path);

  return name;
}

/*----------------------------------------------------------------------------
 * Return the label of variables or properties or scalar for 1D profile
 *
 * parameters:
 *   id           -->  number of 1D profile
 *   nm           -->  number of the variable name of the idst 1D profile
 *----------------------------------------------------------------------------*/

static char *_get_profile_label_name(const int id, const int nm)
{
  char *path = NULL;
  char *name = NULL;
  char *label = NULL;
  int j, shift;

  cs_var_t  *vars = cs_glob_var;
  name = _get_profile_name(id, nm);
  shift = vars->nvar - vars->nscapp - vars->nscaus;

  for (j=0 ; j < shift ; j++) {
    if (cs_gui_strcmp(name,  vars->name[j]))
      label = cs_gui_variable_label(name);
  }

  for (j=0 ; j < vars->nscaus + vars->nscapp; j++) {
    if (cs_gui_strcmp(name, vars->name[j+shift])) {
      BFT_MALLOC(label, strlen(vars->label[j])+1, char);
      strcpy(label, vars->label[j]);
    }
  }

  for (j=0; j < vars->nprop; j++) {
    if (cs_gui_strcmp(name, vars->properties_name[j])) {
      if (label != NULL)
        BFT_FREE(label);
      label = cs_gui_properties_label(vars->properties_name[j]);
    }
  }

  if (label == NULL)
    bft_error(__FILE__, __LINE__, 0,
        _("Invalid markup name: %s\n label not found"), name);

  BFT_FREE(path);
  BFT_FREE(name);

  return label;
}

/*----------------------------------------------------------------------------
 * Get coordinates or output frequency for 1D profile
 *
 * parameters:
 *   id           -->  number of 1D profile
 *    x          <--   name of the coordinate (x1, y1, z1, x2, y2, z2)
 *                     or the output frequency
 *----------------------------------------------------------------------------*/

static double _get_profile_coordinate(const int id, const char *const x)
{
  char *path = NULL;
  double coordinate = 0.0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "analysis_control", "profiles");
  cs_xpath_add_element_num(&path, "profile", id+1);
  cs_xpath_add_element(&path, x);
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &coordinate))
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);

  return coordinate;
}

/*----------------------------------------------------------------------------
 * Return the type of output frequency for 1D profile
 *
 * parameters:
 *   id           -->  number of average
 *----------------------------------------------------------------------------*/

static char *_get_profile_output_type(const int id)
{
    char *path = NULL;
    char *name = NULL;

    path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 2, "analysis_control", "profiles");
    cs_xpath_add_element_num(&path, "profile", id + 1);
    cs_xpath_add_element(&path, "output_type");
    cs_xpath_add_function_text(&path);

    name = cs_gui_get_text_value(path);
    if (name == NULL)
        bft_error(__FILE__, __LINE__, 0,
                  _("Invalid xpath: %s\n name not found"), path);
    BFT_FREE(path);

    return name;
}

/*----------------------------------------------------------------------------
 * Get output format for 1D profile
 *
 * parameters:
 *   id           -->  number of 1D profile
 *----------------------------------------------------------------------------*/

static int _get_profile_format(const int id)
{
  char *path = NULL, *format_s = NULL;
  int   format = 0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "analysis_control", "profiles");
  cs_xpath_add_element_num(&path, "profile", id+1);
  cs_xpath_add_element(&path, "format");
  cs_xpath_add_attribute(&path, "name");
  format_s = cs_gui_get_attribute_value(path);

  if (format_s != NULL) {
    if (cs_gui_strcmp(format_s, "CSV"))
      format = 1;
    else if (cs_gui_strcmp(format_s, "DAT"))
      format = 0;
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Invalid attribute value: %s \nXpath: %s\n"), format_s, path);
    BFT_FREE(format_s);
  }

  BFT_FREE(path);

  return format;
}

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialise the global 'vars' structure.
 *
 * Fortran Interface:
 *
 * subroutine uiinit
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiinit, UIINIT) (void)
{
  assert(cs_glob_var == NULL);

  BFT_MALLOC(cs_glob_var, 1, cs_var_t);

  cs_glob_var->model           = NULL;
  cs_glob_var->model_value     = NULL;
  cs_glob_var->head            = NULL;
  cs_glob_var->type            = NULL;
  cs_glob_var->name            = NULL;
  cs_glob_var->label           = NULL;
  cs_glob_var->rtp             = NULL;
  cs_glob_var->nvar            = 0;
  cs_glob_var->nscaus          = 0;
  cs_glob_var->nscapp          = 0;
  cs_glob_var->nprop           = 0;
  cs_glob_var->nsalpp          = 0;
  cs_glob_var->ntimaver        = 0;
  cs_glob_var->properties_name = NULL;
  cs_glob_var->properties_ipp  = NULL;
  cs_glob_var->propce          = NULL;

  BFT_MALLOC(cs_glob_label, 1, cs_label_t);

  cs_glob_label->_cs_gui_max_vars = 0;
  cs_glob_label->_cs_gui_last_var = 0;
  cs_glob_label->_cs_gui_var_name = NULL;
}

/*----------------------------------------------------------------------------
 * Turbulence model.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSTURB (ITURB, IDEUCH, IGRAKE, IGRAKI, XLOMLG)
 * *****************
 *
 * INTEGER          ITURB   <--   turbulence model
 * INTEGER          IDEUCH  <--   wall law treatment
 * INTEGER          IGRAKE  <--   k-eps gravity effects
 * INTEGER          IGRAKI  <--   Rij-eps gravity effects
 * DOUBLE PRECISION XLOMLG  <--   mixing_length_scale
 *----------------------------------------------------------------------------*/


void CS_PROCF (csturb, CSTURB) (int    *const iturb,
                                int    *const ideuch,
                                int    *const igrake,
                                int    *const igrari,
                                double *const xlomlg)
{
  char *model = NULL;
  char *flux_model = NULL;

  model = cs_gui_get_thermophysical_model("turbulence");
  if (model == NULL)
    return;

  if (cs_gui_strcmp(model, "off"))
    *iturb = 0;
  else if (cs_gui_strcmp(model, "mixing_length")) {
    *iturb = 10;
    _option_turbulence_double("mixing_length_scale", xlomlg);
  } else if (cs_gui_strcmp(model, "k-epsilon")) {
    *iturb = 20;
    cs_gui_advanced_options_turbulence("scale_model", ideuch);
    cs_gui_advanced_options_turbulence("gravity_terms", igrake);
  } else if (cs_gui_strcmp(model, "k-epsilon-PL")) {
    *iturb = 21;
    cs_gui_advanced_options_turbulence("scale_model", ideuch);
    cs_gui_advanced_options_turbulence("gravity_terms", igrake);
  } else if (cs_gui_strcmp(model, "Rij-epsilon")) {
    *iturb = 30;
    cs_gui_advanced_options_turbulence("scale_model", ideuch);
    cs_gui_advanced_options_turbulence("gravity_terms", igrari);
  } else if (cs_gui_strcmp(model, "Rij-SSG")) {
    *iturb = 31;
    cs_gui_advanced_options_turbulence("scale_model", ideuch);
    cs_gui_advanced_options_turbulence("gravity_terms", igrari);
  } else if (cs_gui_strcmp(model, "Rij-EBRSM")) {
    *iturb = 32;
    cs_gui_advanced_options_turbulence("scale_model", ideuch);
    cs_gui_advanced_options_turbulence("gravity_terms", igrari);
  } else if (cs_gui_strcmp(model, "LES_Smagorinsky")) {
    *iturb = 40;
  } else if (cs_gui_strcmp(model, "LES_dynamique")) {
    *iturb = 41;
  } else if (cs_gui_strcmp(model, "LES_WALE")) {
    *iturb = 42;
  } else if (cs_gui_strcmp(model, "v2f-phi")) {
    *iturb = 50;
    cs_gui_advanced_options_turbulence("scale_model", ideuch);
    cs_gui_advanced_options_turbulence("gravity_terms", igrake);
  } else if (cs_gui_strcmp(model, "v2f-BL-v2/k")) {
    *iturb = 51;
    cs_gui_advanced_options_turbulence("scale_model", ideuch);
    cs_gui_advanced_options_turbulence("gravity_terms", igrake);
  } else if (cs_gui_strcmp(model, "k-omega-SST")) {
    *iturb = 60;
    cs_gui_advanced_options_turbulence("scale_model", ideuch);
    cs_gui_advanced_options_turbulence("gravity_terms", igrake);
  } else if (cs_gui_strcmp(model, "Spalart-Allmaras")) {
    *iturb = 70;
  } else
    bft_error(__FILE__, __LINE__, 0,
        _("Invalid turbulence model: %s.\n"), model);

#if _XML_DEBUG_
  bft_printf("==>CSTURB\n");
  bft_printf("--model: %s\n", model);
  bft_printf("--iturb = %i\n", *iturb);
  bft_printf("--igrake = %i\n", *igrake);
  bft_printf("--igrari = %i\n", *igrari);
  bft_printf("--ideuch = %i\n", *ideuch);
  bft_printf("--xlomlg = %f\n", *xlomlg);
#endif

  BFT_FREE(model);
  BFT_FREE(flux_model);
}

/*----------------------------------------------------------------------------
 * Specific heat variable or constant indicator.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSCPVA (ICP)
 * *****************
 *
 * INTEGER          ICP     <--   specific heat variable or constant indicator
 *----------------------------------------------------------------------------*/

void CS_PROCF (cscpva, CSCPVA) (int *const icp)
{
  int choice;

  if (cs_gui_properties_choice("specific_heat", &choice)) *icp = choice;

#if _XML_DEBUG_
  bft_printf("==>CSCPVA\n");
  bft_printf("--icp = %i\n", *icp);
#endif
}

/*----------------------------------------------------------------------------
 * Volumic viscosity variable or constant indicator.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSCVVVA (ICP)
 * *****************
 *
 * INTEGER          IVISCV     <--   specific heat variable or constant indicator
 *----------------------------------------------------------------------------*/

void CS_PROCF (csvvva, CSVVVA) (int *const iviscv)
{
  int choice;

  if (cs_gui_properties_choice("volumic_viscosity", &choice))
    *iviscv = choice;

#if _XML_DEBUG_
  bft_printf("==>CSVVVA\n");
  bft_printf("--iviscv = %i\n", *iviscv);
#endif
}

/*----------------------------------------------------------------------------
 * User scalars number.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSNSCA (NSCAUS)
 * *****************
 *
 * INTEGER          NSCAUS     <--   user scalars number
 *----------------------------------------------------------------------------*/

void CS_PROCF (csnsca, CSNSCA) (int *const nscaus)
{
  int   i     = 0;
  char *label = NULL;

  cs_var_t *vars = cs_glob_var;

  *nscaus = cs_gui_get_tag_number("/additional_scalars/scalar", 1);

  cs_glob_var->nscaus = *nscaus;

  BFT_MALLOC(vars->label, *nscaus, char*);

  for (i=0; i < vars->nscaus; i++) {
    label = _scalar_name_label("label", i+1);
    BFT_MALLOC(cs_glob_var->label[i], strlen(label)+1, char);
    strcpy(cs_glob_var->label[i], label);
    BFT_FREE(label);
  }

#if _XML_DEBUG_
  bft_printf("==>CSNSCA\n");
  bft_printf("--user scalars number: %i\n", vars->nscaus);
  for (i=0; i<*nscaus; i++)
    bft_printf("--label of scalar[%i]: %s\n", i, vars->label[i]);
#endif
}

/*----------------------------------------------------------------------------
 * User scalars which are variance.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSISCA (ISCAVR)
 * *****************
 *
 * INTEGER          ISCAVR     <--   user scalars variance array
 *----------------------------------------------------------------------------*/

void CS_PROCF (csisca, CSISCA) (int *const iscavr)
{
  int i;
  int j;
  char *variance = NULL;

  cs_var_t  *vars = cs_glob_var;

  if (vars->nscaus > 0) {

    for (i=0; i < vars->nscaus; i++) {
      variance = cs_gui_scalar_variance(i+1);

      if (variance != NULL) {

        for (j=0 ; j < vars->nscaus ; j++) {

          if (cs_gui_strcmp(variance, vars->label[j])) {

            if (i == j)
              bft_error(__FILE__, __LINE__, 0,
                        _("Scalar: %i and its variance: %i are the same.\n"),
                        i, j);
            iscavr[i] = j + 1;
          }
        }
        BFT_FREE(variance);
      }
    }

#if _XML_DEBUG_
    bft_printf("==>CSISCA\n");
    for (i = 0 ; i < vars->nscaus ; i++)
      bft_printf("--iscavr[%i] = %i \n", i, iscavr[i]);
#endif
  }
  return;
}

/*----------------------------------------------------------------------------
 * Constant or variable indicator for the user scalar laminar viscosity.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSIVIS (ISCAVR, IVISLS, ISCALT, ISCSTH)
 * *****************
 *
 * INTEGER          ISCAVR  <-->  number of the related variance if any
 * INTEGER          IVISLS  <--   indicator for the user scalar viscosity
 * INTEGER          ISCALT  <-->  number of the user thermal scalar if any
 * INTEGER          ISCSTH  <-->  type of the user thermal scalar
 * INTEGER          ITEMPK   -->  rtp index for temperature (in K)
 *----------------------------------------------------------------------------*/

void CS_PROCF (csivis, CSIVIS) (int *const iscavr,
                                int *const ivisls,
                                int *const iscalt,
                                int *const iscsth,
                                int *const itempk)
{
  int i;
  int choice1, choice2;
  int test1, test2;

  cs_var_t  *vars = cs_glob_var;

  if (vars->nscaus > 0) {

    if (cs_gui_thermal_scalar()) {
      test1 = cs_gui_properties_choice("thermal_conductivity", &choice1);
      test2 = cs_gui_properties_choice("specific_heat", &choice2);

      if (test1 && test2) {
        cs_gui_thermal_scalar_number(iscalt, iscsth);

        if (choice1 || choice2)
          ivisls[*iscalt-1] = 1;
        else
          ivisls[*iscalt-1] = 0;
      }
    }

    for (i=0 ; i < vars->nscaus; i++) {
      if (iscavr[i] <= 0 ) {
        if (cs_gui_scalar_properties_choice(i+1, &choice1))
          if (*iscalt != i+1)
            ivisls[i] = choice1;
      }
    }

#if _XML_DEBUG_
    bft_printf("==>CSIVIS\n");
    for (i=0 ; i < vars->nscaus ; i++)
      bft_printf("--ivisls[%i] = %i\n", i, ivisls[i]);
#endif
  }

  if (cs_gui_strcmp(vars->model, "compressible_model"))
  {
    ivisls[*itempk -1] = 0;

    char *prop_choice = _properties_choice("thermal_conductivity");
    if (cs_gui_strcmp(prop_choice, "user_law"))
      ivisls[*itempk -1] = 1;
    BFT_FREE(prop_choice);
  }
}

/*----------------------------------------------------------------------------
 * Time passing parameter.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSIDTV (IDTVAR)
 * *****************
 *
 * INTEGER          IDTVAR  <--   fixed or variable time step
 *----------------------------------------------------------------------------*/

void CS_PROCF (csidtv, CSIDTV) (int *const idtvar)
{
  double param;
  int steady = 0;
  char* algo_choice = NULL;

  cs_gui_get_steady_status(&steady);
  if (steady) {
    algo_choice = cs_gui_velocity_pressure_algo_choice();
    if (cs_gui_strcmp(algo_choice, "simple"))
      *idtvar = -1;
    else
      *idtvar = 2;
  } else {
    param = (double) *idtvar;
    cs_gui_time_parameters("time_passing", &param);
    *idtvar = (int) param;
  }

#if _XML_DEBUG_
  bft_printf("==>CSIDTV\n");
  bft_printf("--idtvar = %i\n", *idtvar);
#endif
}

/*----------------------------------------------------------------------------
 * Hydrostatic pressure parameter.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSIPHY (IPHYDR)
 * *****************
 *
 * INTEGER          IPHYDR  <--   hydrostatic pressure
 *----------------------------------------------------------------------------*/

void CS_PROCF (csiphy, CSIPHY) (int *const iphydr)
{
  char *path = NULL;
  int   result;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "numerical_parameters");
  cs_xpath_add_element(&path, "hydrostatic_pressure");
  cs_xpath_add_attribute(&path, "status");

  if (cs_gui_get_status(path, &result))
    *iphydr = result;

  BFT_FREE(path);

#if _XML_DEBUG_
  bft_printf("==>CSIPHY\n");
  bft_printf("--iphydr = %i\n", *iphydr);
#endif
}

/*----------------------------------------------------------------------------
 * Hydrostatic equilibrium parameter.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSCFGP (icfgrp)
 * *****************
 *
 * INTEGER          icfgrp  <--   hydrostatic equilibrium
 *----------------------------------------------------------------------------*/

void CS_PROCF (cscfgp, CSCFGP) (int *const icfgrp)
{
  char *path = NULL;
  int   result;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "numerical_parameters");
  cs_xpath_add_element(&path, "hydrostatic_equilibrium");
  cs_xpath_add_attribute(&path, "status");

  if (cs_gui_get_status(path, &result)) *icfgrp = result;

  BFT_FREE(path);

#if _XML_DEBUG_
  bft_printf("==>CSCFGP\n");
  bft_printf("--icfgrp = %i\n", *icfgrp);
#endif
}

/*----------------------------------------------------------------------------
 * Constructs an indirection between an internal index vars->rtp and
 * the fortran array RTP.
 * This function is called after the third call to VARPOS routine.
 * Warning: in vars->rtp, variables are stored first and then
 * scalars (not like VARPOS because of ALE variables).
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSVNUM
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (csvnum, CSVNUM) (const int *const nvar,
                                const int *const iu,
                                const int *const iv,
                                const int *const iw,
                                const int *const ipr,
                                const int *const iturb,
                                const int *const ik,
                                const int *const iep,
                                const int *const ir11,
                                const int *const ir22,
                                const int *const ir33,
                                const int *const ir12,
                                const int *const ir13,
                                const int *const ir23,
                                const int *const iomg,
                                const int *const iphi,
                                const int *const ifb,
                                const int *const ial,
                                const int *const inusa,
                                const int *const iale,
                                const int *const iuma,
                                const int *const ivma,
                                const int *const iwma,
                                const int *const isca,
                                const int *const iscapp)
{
  int n = 0;
  int i, j, k;
  char *name = NULL;

  BFT_MALLOC(cs_glob_var->rtp,  *nvar, int);
  BFT_MALLOC(cs_glob_var->head, *nvar, char*);
  BFT_MALLOC(cs_glob_var->type, *nvar, char*);
  BFT_MALLOC(cs_glob_var->name, *nvar, char*);

  /* Warning!!  for scalars: vars->nscaus is already fill in CSNSCA */
  /*                         vars->nscapp is already fill in UIPPMO */

  cs_glob_var->nvar   = *nvar;

  /* 1) pressure and velocity variables */

  k = n;
  cs_glob_var->rtp[n] = *ipr -1;
  BFT_MALLOC(cs_glob_var->name[n], strlen("pressure")+1, char);
  strcpy(cs_glob_var->name[n++], "pressure");

  cs_glob_var->rtp[n] = *iu  -1;
  BFT_MALLOC(cs_glob_var->name[n], strlen("velocity_U")+1, char);
  strcpy(cs_glob_var->name[n++], "velocity_U");

  cs_glob_var->rtp[n] = *iv  -1;
  BFT_MALLOC(cs_glob_var->name[n], strlen("velocity_V")+1, char);
  strcpy(cs_glob_var->name[n++], "velocity_V");

  cs_glob_var->rtp[n] = *iw  -1;
  BFT_MALLOC(cs_glob_var->name[n], strlen("velocity_W")+1, char);
  strcpy(cs_glob_var->name[n++], "velocity_W");

  for (i=k; i < n; i++) {
    BFT_MALLOC(cs_glob_var->head[i], strlen("velocity_pressure")+1, char);
    strcpy(cs_glob_var->head[i], "velocity_pressure");
  }

  /* 2) turbulence variables */

  k = n;

  if (*iturb == 20 || *iturb == 21) {

    cs_glob_var->rtp[n] = *ik  -1;
    BFT_MALLOC(cs_glob_var->name[n], strlen("turb_k")+1, char);
    strcpy(cs_glob_var->name[n++], "turb_k");

    cs_glob_var->rtp[n] = *iep -1;
    BFT_MALLOC(cs_glob_var->name[n], strlen("turb_eps")+1, char);
    strcpy(cs_glob_var->name[n++], "turb_eps");

  } else if (*iturb == 30 || *iturb == 31 || *iturb == 32) {

    cs_glob_var->rtp[n] = *ir11 -1;
    BFT_MALLOC(cs_glob_var->name[n], strlen("component_R11")+1, char);
    strcpy(cs_glob_var->name[n++], "component_R11");

    cs_glob_var->rtp[n] = *ir22 -1;
    BFT_MALLOC(cs_glob_var->name[n], strlen("component_R22")+1, char);
    strcpy(cs_glob_var->name[n++], "component_R22");

    cs_glob_var->rtp[n] = *ir33 -1;
    BFT_MALLOC(cs_glob_var->name[n], strlen("component_R33")+1, char);
    strcpy(cs_glob_var->name[n++], "component_R33");

    cs_glob_var->rtp[n] = *ir12 -1;
    BFT_MALLOC(cs_glob_var->name[n], strlen("component_R12")+1, char);
    strcpy(cs_glob_var->name[n++], "component_R12");

    cs_glob_var->rtp[n] = *ir13 -1;
    BFT_MALLOC(cs_glob_var->name[n], strlen("component_R13")+1, char);
    strcpy(cs_glob_var->name[n++], "component_R13");

    cs_glob_var->rtp[n] = *ir23 -1;
    BFT_MALLOC(cs_glob_var->name[n], strlen("component_R23")+1, char);
    strcpy(cs_glob_var->name[n++], "component_R23");

    cs_glob_var->rtp[n] = *iep  -1;
    BFT_MALLOC(cs_glob_var->name[n], strlen("turb_eps")+1, char);
    strcpy(cs_glob_var->name[n++], "turb_eps");

    if (*iturb == 32) {
      cs_glob_var->rtp[n] = *ial  -1;
      BFT_MALLOC(cs_glob_var->name[n], strlen("turb_alpha")+1, char);
      strcpy(cs_glob_var->name[n++], "turb_alpha");
    }

  } else if (*iturb == 50) {

    cs_glob_var->rtp[n] = *ik   -1;
    BFT_MALLOC(cs_glob_var->name[n], strlen("turb_k")+1, char);
    strcpy(cs_glob_var->name[n++], "turb_k");

    cs_glob_var->rtp[n] = *iep  -1;
    BFT_MALLOC(cs_glob_var->name[n], strlen("turb_eps")+1, char);
    strcpy(cs_glob_var->name[n++], "turb_eps");

    cs_glob_var->rtp[n] = *iphi -1;
    BFT_MALLOC(cs_glob_var->name[n], strlen("turb_phi")+1, char);
    strcpy(cs_glob_var->name[n++], "turb_phi");

    cs_glob_var->rtp[n] = *ifb  -1;
    BFT_MALLOC(cs_glob_var->name[n], strlen("turb_fb")+1, char);
    strcpy(cs_glob_var->name[n++], "turb_fb");

  } else if (*iturb == 51) {

    cs_glob_var->rtp[n] = *ik   -1;
    BFT_MALLOC(cs_glob_var->name[n], strlen("turb_k")+1, char);
    strcpy(cs_glob_var->name[n++], "turb_k");

    cs_glob_var->rtp[n] = *iep  -1;
    BFT_MALLOC(cs_glob_var->name[n], strlen("turb_eps")+1, char);
    strcpy(cs_glob_var->name[n++], "turb_eps");

    cs_glob_var->rtp[n] = *iphi -1;
    BFT_MALLOC(cs_glob_var->name[n], strlen("turb_phi")+1, char);
    strcpy(cs_glob_var->name[n++], "turb_phi");

    cs_glob_var->rtp[n] = *ial  -1;
    BFT_MALLOC(cs_glob_var->name[n], strlen("turb_alpha")+1, char);
    strcpy(cs_glob_var->name[n++], "turb_alpha");

  } else if (*iturb == 60) {

    cs_glob_var->rtp[n] = *ik   -1;
    BFT_MALLOC(cs_glob_var->name[n], strlen("turb_k")+1, char);
    strcpy(cs_glob_var->name[n++], "turb_k");

    cs_glob_var->rtp[n] = *iomg -1;
    BFT_MALLOC(cs_glob_var->name[n], strlen("turb_omega")+1, char);
    strcpy(cs_glob_var->name[n++], "turb_omega");

  } else if (*iturb == 70) {

    cs_glob_var->rtp[n] = *inusa   -1;
    BFT_MALLOC(cs_glob_var->name[n], strlen("turb_nusa")+1, char);
    strcpy(cs_glob_var->name[n++], "turb_nusa");
  }

  for (i=k; i < n; i++) {
    BFT_MALLOC(cs_glob_var->head[i], strlen("turbulence")+1, char);
    strcpy(cs_glob_var->head[i], "turbulence");
  }

  /* 3) ALE variables */

  if (*iale) {
    k = n;

    cs_glob_var->rtp[n] = *iuma -1;
    BFT_MALLOC(cs_glob_var->name[n], strlen("mesh_velocity_U")+1, char);
    strcpy(cs_glob_var->name[n++], "mesh_velocity_U");

    cs_glob_var->rtp[n] = *ivma -1;
    BFT_MALLOC(cs_glob_var->name[n], strlen("mesh_velocity_V")+1, char);
    strcpy(cs_glob_var->name[n++], "mesh_velocity_V");

    cs_glob_var->rtp[n] = *iwma -1;
    BFT_MALLOC(cs_glob_var->name[n], strlen("mesh_velocity_W")+1, char);
    strcpy(cs_glob_var->name[n++], "mesh_velocity_W");

    for (i=k; i < n; i++) {
      BFT_MALLOC(cs_glob_var->head[i], strlen("ale_method")+1, char);
      strcpy(cs_glob_var->head[i], "ale_method");
    }
  }

  /* 4) update vars->type for variables */

  k = cs_glob_var->nvar -cs_glob_var->nscapp -cs_glob_var->nscaus;
  for (i=0; i < k; i++) {
    BFT_MALLOC(cs_glob_var->type[i], strlen("variable")+1, char);
    strcpy(cs_glob_var->type[i], "variable");
  }

  /* 5) user scalars */

  for (i=0; i < cs_glob_var->nscaus; i++) {
    cs_glob_var->rtp[n++] = isca[i] -1;

    name = _scalar_name_label("name", i+1);
    BFT_MALLOC(cs_glob_var->name[k+i], strlen(name) +1, char);
    strcpy(cs_glob_var->name[k+i], name);
    BFT_FREE(name);

    BFT_MALLOC(cs_glob_var->type[k+i], strlen("scalar")+1, char);
    strcpy(cs_glob_var->type[k+i], "scalar");

    BFT_MALLOC(cs_glob_var->head[k+i], strlen("additional_scalar")+1, char);
    strcpy(cs_glob_var->head[k+i], "additional_scalar");
  }

  /* 6) model scalars */

  k = cs_glob_var->nvar -cs_glob_var->nscaus - cs_glob_var->nscapp;
  for (i=0; i < cs_glob_var->nscapp; i++) {
    j = iscapp[i] -1;
    cs_glob_var->rtp[n++] = isca[j] -1;

    name = _specific_physic_scalar_name_label(cs_glob_var->model, cs_glob_var->label[j]);
    BFT_MALLOC(cs_glob_var->name[k+j], strlen(name) +1, char);
    strcpy(cs_glob_var->name[k+j], name);
    BFT_FREE(name);

    BFT_MALLOC(cs_glob_var->type[k+j], strlen("scalar")+1, char);
    strcpy(cs_glob_var->type[k+j], "scalar");

    BFT_MALLOC(cs_glob_var->head[k+j], strlen(cs_glob_var->model)+1, char);
    strcpy(cs_glob_var->head[k+j], cs_glob_var->model);
  }

  /* 7) check for errors */

  if (n != *nvar)
    bft_error(__FILE__, __LINE__, 0,
              _("The kernel variables number %i and the "
              "calculated one by the GUI %i are not the same.\n"),
              *nvar, n);

#if _XML_DEBUG_
  bft_printf("==>CSVNUM\n");
  bft_printf("--variables and scalars name: %i\n", cs_glob_var->nvar);
  for (i=0; i < cs_glob_var->nvar; i++)
    bft_printf("---name: %s\n", cs_glob_var->name[i]);
#endif

}

/*----------------------------------------------------------------------------
 * Restart parameters.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSISUI (NTSUIT, ILEAUX, ICCVFG)
 * *****************
 *
 * INTEGER          NTSUIT  <--   checkpoint frequency
 * INTEGER          ILEAUX  <--   restart with auxiliary
 * INTEGER          ICCFVG  <--   restart with frozen field
 *----------------------------------------------------------------------------*/

void CS_PROCF (csisui, CSISUI) (int *const ntsuit,
                                int *const ileaux,
                                int *const iccvfg)
{
  cs_gui_restart_parameters_status("restart_rescue",         ntsuit);
  cs_gui_restart_parameters_status("restart_with_auxiliary", ileaux);
  cs_gui_restart_parameters_status("frozen_field",           iccvfg);

#if _XML_DEBUG_
  bft_printf("==>CSISUI\n");
  bft_printf("--ntsuit = %i\n", *ntsuit);
  bft_printf("--ileaux = %i\n", *ileaux);
  bft_printf("--iccvfg = %i\n", *iccvfg);
#endif
}

/*----------------------------------------------------------------------------
 * Time passing parameters.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSTIME (INPDT0, IPTLTO, NTMABS, DTREF,
 * *****************  DTMIN,  DTMAX,  COUMAX, FOUMAX, VARRDT)
 *
 * INTEGER          INPDT0  <--   zero time step
 * INTEGER          IPTLTO  <--   thermal time step control
 * INTEGER          NTMABS  <--   iterations numbers
 * INTEGER          IDTVAR  <--   time steps'options
 * DOUBLE PRECISION DTREF   <--   time step
 * DOUBLE PRECISION DTMIN   <--   minimal time step
 * DOUBLE PRECISION DTMAX   <--   maximal time step
 * DOUBLE PRECISION COUMAX  <--   maximal courant number
 * DOUBLE PRECISION FOUMAX  <--   maximal fournier number
 * DOUBLE PRECISION VARRDT  <--   max time step variation between 2 iterations
 * DOUBLE PRECISION RELXST  <--   relaxation coefficient id idtvar = -1
 *----------------------------------------------------------------------------*/

void CS_PROCF (cstime, CSTIME) (int    *const inpdt0,
                                int    *const iptlro,
                                int    *const ntmabs,
                                int    *const idtvar,
                                double *const dtref,
                                double *const dtmin,
                                double *const dtmax,
                                double *const coumax,
                                double *const foumax,
                                double *const varrdt,
                                double *const relxst)
{
  double value;
  /* Default values for time step factor */
  double cdtmin = 0.1, cdtmax = 1000.;

  if (*idtvar == -1) {
    cs_gui_steady_parameters("relaxation_coefficient", relxst);

    value =(double) *inpdt0;
    cs_gui_steady_parameters("zero_iteration", &value);
    *inpdt0 = (int) value;

    value =(double) *ntmabs;
    cs_gui_steady_parameters("iterations", &value);
    *ntmabs = (int) value;
  } else {
    cs_gui_time_parameters("time_step_ref", dtref);
    cs_gui_time_parameters("time_step_min_factor", &cdtmin);
    cs_gui_time_parameters("time_step_max_factor", &cdtmax);
    cs_gui_time_parameters("max_courant_num", coumax);
    cs_gui_time_parameters("max_fourier_num", foumax);
    cs_gui_time_parameters("time_step_var", varrdt);

    *dtmin = cdtmin*(*dtref);
    *dtmax = cdtmax*(*dtref);

    /* We keep these two lines in case we read an old XML file... */
    cs_gui_time_parameters("time_step_min", dtmin);
    cs_gui_time_parameters("time_step_max", dtmax);

    value =(double) *ntmabs;
    cs_gui_time_parameters("iterations", &value);
    *ntmabs = (int) value;

    value =(double) *inpdt0;
    cs_gui_time_parameters("zero_time_step", &value);
    *inpdt0 = (int) value;

    value =(double) *iptlro;
    cs_gui_time_parameters("thermal_time_step", &value);
    *iptlro = (int) value;
  }

#if _XML_DEBUG_
  bft_printf("==>CSTIME\n");
  bft_printf("--idtvar = %i\n", *idtvar);
  if (*idtvar == -1) {
    bft_printf("--inpdt0 = %i\n", *inpdt0);
    bft_printf("--relxst = %f\n", *relxst);
  } else {
    bft_printf("--inpdt0 = %i\n", *inpdt0);
    bft_printf("--iptlro = %i\n", *iptlro);
    bft_printf("--ntmabs = %i\n", *ntmabs);
    bft_printf("--dtref = %f\n",  *dtref);
    bft_printf("--dtmin = %f\n",  *dtmin);
    bft_printf("--dtmax = %f\n",  *dtmax);
    bft_printf("--coumax = %f\n", *coumax);
    bft_printf("--foumax = %f\n", *foumax);
    bft_printf("--varrdt = %f\n", *varrdt);
  }
#endif
}

/*----------------------------------------------------------------------------
 * Check if a users thermal scalar is defined.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSSCA1 (ISCALT, ISCSTH)
 * *****************
 *
 * INTEGER          ISCALT  <--   number of the user thermal scalar if any
 * INTEGER          ISCSTH  <--   type of the user thermal scalar
 *----------------------------------------------------------------------------*/

void CS_PROCF (cssca1, CSSCA1) (int *const iscalt,
                                int *const iscsth)
{
  cs_gui_thermal_scalar_number(iscalt, iscsth);

#if _XML_DEBUG_
  {
    int i;
    cs_var_t  *vars = cs_glob_var;
    bft_printf("==>CSSCA1\n");
    bft_printf("--iscalt[0]=%i \n", *iscalt);
    for (i = 0 ; i < vars->nscaus ; i++)
      bft_printf("--iscsth[%i]=%i \n", i, iscsth[i]);
  }
#endif
}

/*----------------------------------------------------------------------------
 * Treatment of local numerical aspects:
 *     BLENCV, ISCHCV, ISSTPC, IRCFLU, CDTVAR, NITMAX, EPSILO
 *----------------------------------------------------------------------------*/

void CS_PROCF (uinum1, UINUM1) (const    int *const isca,
                                const    int *const iscapp,
                                      double *const blencv,
                                         int *const ischcv,
                                         int *const isstpc,
                                         int *const ircflu,
                                      double *const cdtvar,
                                         int *const nitmax,
                                      double *const epsilo,
                                         int *const iresol,
                                         int *const imgr,
                                         int *const nswrsm)
{
  int i, j, jj, k;
  double tmp;
  char* algo_choice = NULL;

  cs_var_t  *vars = cs_glob_var;

  k = vars->nvar - vars->nscaus - vars->nscapp;

  /* 1) variables from velocity_pressure and turbulence */
  /* 1-a) for pressure */
  j = vars->rtp[0];
  cs_gui_variable_value(vars->name[0], "solver_precision", &epsilo[j]);
  tmp = (double) nitmax[j];
  cs_gui_variable_value(vars->name[0], "max_iter_number", &tmp);
  nitmax[j] = (int) tmp;

  imgr[j] = 0;

  algo_choice = cs_gui_variable_choice(vars->name[0], "solver_choice");
  if (cs_gui_strcmp(algo_choice, "multigrid"))
  {
    iresol[j] = 0;
    imgr[j] = 1;
  }
  else if (cs_gui_strcmp(algo_choice, "conjugate_gradient"))
    iresol[j] = 0;
  else if (cs_gui_strcmp(algo_choice, "jacobi"))
    iresol[j] = 1;
  else if (cs_gui_strcmp(algo_choice, "bi_cgstab"))
    iresol[j] = 2;
  else if (cs_gui_strcmp(algo_choice, "gmres"))
    iresol[j] = 3;
  else if (cs_gui_strcmp(algo_choice, "automatic"))
    iresol[j] = -1;
  else //default value
  {
    iresol[j] = 0;
    imgr[j] = 1;
  }
  tmp = (double) nswrsm[j];
  cs_gui_variable_value(vars->name[0], "rhs_reconstruction", &tmp);
  nswrsm[j] = (int) tmp;

  /* 1-b) for the other variables */
  for (i=1; i < k; i++) {
    j = vars->rtp[i];
    cs_gui_variable_value(vars->name[i], "blending_factor", &blencv[j]);
    cs_gui_variable_value(vars->name[i], "solver_precision", &epsilo[j]);

    imgr[j] = 0;

    algo_choice = cs_gui_variable_choice(vars->name[i], "solver_choice");

    if (cs_gui_strcmp(algo_choice, "conjugate_gradient"))
      iresol[j] = 0;
    else if (cs_gui_strcmp(algo_choice, "jacobi"))
      iresol[j] = 1;
    else if (cs_gui_strcmp(algo_choice, "bi_cgstab"))
      iresol[j] = 2;
    else if (cs_gui_strcmp(algo_choice, "gmres"))
      iresol[j] = 3;
    else if (cs_gui_strcmp(algo_choice, "automatic"))
      iresol[j] = -1;
    else //default value
      iresol[j] = -1;

    tmp = (double) nitmax[j];
    cs_gui_variable_value(vars->name[i], "max_iter_number", &tmp);
    nitmax[j] = (int) tmp;
    cs_gui_variable_attribute(vars->name[i], "order_scheme", &ischcv[j]);
    cs_gui_variable_attribute(vars->name[i], "slope_test", &isstpc[j]);
    cs_gui_variable_attribute(vars->name[i], "flux_reconstruction", &ircflu[j]);
    tmp = (double) nswrsm[j];
    cs_gui_variable_value(vars->name[i], "rhs_reconstruction", &tmp);
    nswrsm[j] = (int) tmp;
  }

  /* 2) user scalars */

  if (vars->nscaus > 0 ) {
    for (i=0 ; i < vars->nscaus; i++) {
      j = isca[i]-1;
      cs_gui_scalar_value(vars->label[i], "blending_factor", &blencv[j]);
      cs_gui_scalar_value(vars->label[i], "solver_precision", &epsilo[j]);

      imgr[j] = 0;

      algo_choice = cs_gui_variable_choice(vars->label[i], "solver_choice");
      if (cs_gui_strcmp(algo_choice, "conjugate_gradient"))
        iresol[j] = 0;
      else if (cs_gui_strcmp(algo_choice, "jacobi"))
        iresol[j] = 1;
      else if (cs_gui_strcmp(algo_choice, "bi_cgstab"))
        iresol[j] = 2;
      else if (cs_gui_strcmp(algo_choice, "gmres"))
        iresol[j] = 3;
      else if (cs_gui_strcmp(algo_choice, "automatic"))
        iresol[j] = -1;
      else //default value
        iresol[j] = -1;

      cs_gui_scalar_value(vars->label[i], "time_step_factor", &cdtvar[j]);
      tmp = (double) nitmax[j];
      cs_gui_scalar_value(vars->label[i], "max_iter_number", &tmp);
      nitmax[j] = (int) tmp;
      cs_gui_scalar_attribute(vars->label[i], "order_scheme", &ischcv[j]);
      cs_gui_scalar_attribute(vars->label[i], "slope_test", &isstpc[j]);
      cs_gui_scalar_attribute(vars->label[i], "flux_reconstruction", &ircflu[j]);
      tmp = (double) nswrsm[j];
      cs_gui_scalar_value(vars->label[i], "rhs_reconstruction", &tmp);
      nswrsm[j] = (int) tmp;
    }
  }

  /* 3) model scalars */

  if (vars->nscapp > 0 ) {
    for (i=0 ; i < vars->nscapp ; i++) {
      j = iscapp[i] -1;
      jj = isca[j]-1;
      cs_gui_model_scalar_value(vars->model, vars->label[j], "blending_factor", &blencv[jj]);
      cs_gui_model_scalar_value(vars->model, vars->label[j], "solver_precision", &epsilo[jj]);
      cs_gui_model_scalar_value(vars->model, vars->label[j], "time_step_factor", &cdtvar[jj]);

      imgr[jj] = 0;

      algo_choice = cs_gui_model_variable_choice(vars->label[j], "solver_choice");
      if (cs_gui_strcmp(algo_choice, "conjugate_gradient"))
        iresol[jj] = 0;
      else if (cs_gui_strcmp(algo_choice, "jacobi"))
        iresol[jj] = 1;
      else if (cs_gui_strcmp(algo_choice, "bi_cgstab"))
        iresol[jj] = 2;
      else if (cs_gui_strcmp(algo_choice, "gmres"))
        iresol[jj] = 3;
      else if (cs_gui_strcmp(algo_choice, "automatic"))
        iresol[j] = -1;
      else //default value
        iresol[jj] = -1;

      tmp = (double) nitmax[jj];
      cs_gui_model_scalar_value(vars->model, vars->label[j], "max_iter_number", &tmp);
      nitmax[jj] = (int) tmp;
      cs_gui_model_scalar_output_status(vars->model, vars->label[j], "order_scheme", &ischcv[jj]);
      cs_gui_model_scalar_output_status(vars->model, vars->label[j], "slope_test", &isstpc[jj]);
      cs_gui_model_scalar_output_status(vars->model, vars->label[j], "flux_reconstruction", &ircflu[jj]);
      tmp = (double) nswrsm[jj];
      cs_gui_model_scalar_value(vars->model, vars->label[j], "rhs_reconstruction", &tmp);
      nswrsm[jj] = (int) tmp;
    }
  }

#if _XML_DEBUG_
  bft_printf("==>UINUM1\n");
  for (i=0; i < vars->nvar; i++) {
    bft_printf("-->variable[%i] = %s\n", i, vars->name[i]);
    bft_printf("--blencv = %f\n", blencv[vars->rtp[i]]);
    bft_printf("--epsilo = %g\n", epsilo[vars->rtp[i]]);
    bft_printf("--cdtvar = %g\n", cdtvar[vars->rtp[i]]);
    bft_printf("--nitmax = %i\n", nitmax[vars->rtp[i]]);
    bft_printf("--ischcv = %i\n", ischcv[vars->rtp[i]]);
    bft_printf("--isstpc = %i\n", isstpc[vars->rtp[i]]);
    bft_printf("--ircflu = %i\n", ircflu[vars->rtp[i]]);
    bft_printf("--nswrsm = %i\n", nswrsm[vars->rtp[i]]);
    bft_printf("--imgr = %i\n"  , imgr[vars->rtp[i]]);
    bft_printf("--iresol = %i\n", iresol[vars->rtp[i]]);
  }
  for (i=0 ; i < vars->nscaus + vars->nscapp ; i++) {
    bft_printf("-->scalar[%i]: %s\n", isca[i]-1, vars->label[i]);
    bft_printf("--blencv = %f\n", blencv[isca[i]-1]);
    bft_printf("--epsilo = %g\n", epsilo[isca[i]-1]);
    bft_printf("--cdtvar = %g\n", cdtvar[isca[i]-1]);
    bft_printf("--nitmax = %i\n", nitmax[isca[i]-1]);
    bft_printf("--ischcv = %i\n", ischcv[isca[i]-1]);
    bft_printf("--isstpc = %i\n", isstpc[isca[i]-1]);
    bft_printf("--ircflu = %i\n", ircflu[isca[i]-1]);
    bft_printf("--nswrsm = %i\n", nswrsm[isca[i]-1]);
    bft_printf("--imgr = %i\n"  , imgr[isca[i]-1]);
    bft_printf("--iresol = %i\n", iresol[isca[i]-1]);
  }
#endif
}

/*----------------------------------------------------------------------------
 * Global numerical parameters.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSNUM2 (IVISSE, RELAXP, IPUCOU, EXTRAG, IMRGRA, IMGRPR)
 * *****************
 * INTEGER          IVISSE  <--   gradient transposed
 * DOUBLE PRECISION RELAXP  <--   pressure relaxation
 * INTEGER          IPUCOU  <--   velocity pressure coupling
 * DOUBLE PRECISION EXTRAG  <--   wall pressure extrapolation
 * INTEGER          IMRGRA  <--   gradient reconstruction
 * INTEGER          IMGRPR  <--   multigrid algorithm for pressure
 * INTEGER          NTERUP  <--   piso sweep number
 *----------------------------------------------------------------------------*/

void CS_PROCF (csnum2, CSNUM2)(   int *const ivisse,
                               double *const relaxp,
                                  int *const ipucou,
                               double *const extrag,
                                  int *const imrgra,
                                  int *const nterup)
{
  cs_gui_numerical_int_parameters("gradient_transposed", ivisse);
  cs_gui_numerical_int_parameters("velocity_pressure_coupling", ipucou);
  cs_gui_numerical_int_parameters("gradient_reconstruction", imrgra);
  cs_gui_numerical_int_parameters("piso_sweep_number", nterup);
  cs_gui_numerical_double_parameters("wall_pressure_extrapolation", extrag);
  cs_gui_numerical_double_parameters("pressure_relaxation", relaxp);

#if _XML_DEBUG_
  bft_printf("==>CSNUM2\n");
  bft_printf("--ivisse = %i\n", *ivisse);
  bft_printf("--ipucou = %i\n", *ipucou);
  bft_printf("--imrgra = %i\n", *imrgra);
  bft_printf("--nterup = %i\n", *nterup);
  bft_printf("--extrag = %f\n", *extrag);
  bft_printf("--relaxp = %f\n", *relaxp);
#endif
}

/*----------------------------------------------------------------------------
 * Treatment of gravity and fluid physical properties
 * Initialize reference pressure and temperature if present
 *----------------------------------------------------------------------------*/

void CS_PROCF (csphys, CSPHYS)
  (
    const    int *const nmodpp,
             int *const irovar,
             int *const ivivar,
             int *const icorio,
          double *const gx,
          double *const gy,
          double *const gz,
          double *const omegax,
          double *const omegay,
          double *const omegaz,
          double *const ro0,
          double *const viscl0,
          double *const viscv0,
          double *const visls0,
          double *const cp0,
          double *const t0,
          double *const p0,
          double *const xmasmr,
             int *const itempk)

{
  int choice;

  cs_var_t  *vars = cs_glob_var;

  cs_gui_gravity_value("gravity_x", gx);
  cs_gui_gravity_value("gravity_y", gy);
  cs_gui_gravity_value("gravity_z", gz);

  cs_gui_coriolis_value("omega_x", omegax);
  cs_gui_coriolis_value("omega_y", omegay);
  cs_gui_coriolis_value("omega_z", omegaz);

  if (   cs_gui_is_equal_real(*omegax, 0.)
      && cs_gui_is_equal_real(*omegay, 0.)
      && cs_gui_is_equal_real(*omegaz, 0.))
    *icorio = 0;
  else
    *icorio = 1;

  cs_gui_properties_value("density", ro0);
  cs_gui_properties_value("molecular_viscosity", viscl0);
  cs_gui_properties_value("specific_heat", cp0);
  if (cs_gui_strcmp(vars->model, "compressible_model")) {
    cs_gui_properties_value("volumic_viscosity", viscv0);
    cs_gui_properties_value("thermal_conductivity", &visls0[*itempk -1]);
  }

  cs_gui_reference_initialization("pressure", p0);

  /* Variable rho and viscl */
  if (*nmodpp == 0) {
    if (cs_gui_properties_choice("density", &choice))
      *irovar = choice;

    if (cs_gui_properties_choice("molecular_viscosity", &choice))
      *ivivar = choice;
  }
  if (cs_gui_strcmp(vars->model, "compressible_model"))
    if (cs_gui_properties_choice("molecular_viscosity", &choice))
      *ivivar = choice;

  /* T0 if necessary */

  if (vars->model != NULL)
    cs_gui_reference_initialization("temperature", t0);
  if (cs_gui_strcmp(vars->model, "compressible_model"))
    cs_gui_reference_initialization("mass_molar", xmasmr);

#if _XML_DEBUG_
  bft_printf("==>CSPHYS\n");
  bft_printf("--gx = %f \n",*gx);
  bft_printf("--gy = %f \n",*gy);
  bft_printf("--gz = %f \n",*gz);
  bft_printf("--omegax = %f \n",*omegax);
  bft_printf("--omegay = %f \n",*omegay);
  bft_printf("--omegaz = %f \n",*omegaz);
  bft_printf("--rho = %g , variable %i\n", *ro0, *irovar);
  bft_printf("--mu = %g , variable %i \n", *viscl0, *ivivar);
  bft_printf("--icorio = %i \n", *icorio);
  bft_printf("--Cp = %g \n", *cp0);
  bft_printf("--T0 = %f \n", *t0);
  bft_printf("--P0 = %f \n", *p0);
  if (cs_gui_strcmp(vars->model, "compressible_model")) {
    bft_printf("--viscv0 = %g \n", *viscv0);
    bft_printf("--xmasmr = %f \n", *xmasmr);
  }
#endif
}

/*----------------------------------------------------------------------------
 * User scalar min and max values for clipping.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSSCA2 (ISCAVR, SCAMIN, SCAMAX)
 * *****************
 *
 * INTEGER          ISCAVR   -->  number of the related variance if any
 * DOUBLE PRECISION SCAMIN  <--   user scalar min array
 * DOUBLE PRECISION SCAMAX  <--   user scalar max array
 *----------------------------------------------------------------------------*/

void CS_PROCF (cssca2, CSSCA2) (const    int *const iscavr,
                                      double *const scamin,
                                      double *const scamax)
{
  /* Specific physics: the min max of the model scalar are not given */

  int i;
  cs_var_t  *vars = cs_glob_var;

  if (vars->nscaus > 0 ) {
    for (i=0 ; i < vars->nscaus; i++) {
      if (iscavr[i] <= 0 ) {
        cs_gui_scalar_value(vars->label[i], "min_value", &scamin[i]);
        cs_gui_scalar_value(vars->label[i], "max_value", &scamax[i]);
      }
    }

#if _XML_DEBUG_
    bft_printf("==>CSSCA2\n");
    for (i=0 ; i < vars->nscaus ; i++) {
      bft_printf("--scamin[%i] = %f\n", i, scamin[i]);
      bft_printf("--scamax[%i] = %f\n", i, scamax[i]);
    }
#endif
  }
}


/*----------------------------------------------------------------------------
 * Read reference dynamic and user scalar viscosity
 *----------------------------------------------------------------------------*/

void CS_PROCF (cssca3, CSSCA3) (const    int *const iscalt,
                                const    int *const iscsth,
                                const    int *const iscavr,
                                double *const visls0,
                                double *const t0,
                                double *const p0)
{
  int i;
  double result, coeff, density;

  cs_var_t  *vars = cs_glob_var;

  if (vars->nscaus > 0) {

    if (cs_gui_thermal_scalar()) {
      result = 0;
      cs_gui_properties_value("specific_heat", &result);
      if (result <= 0)
        bft_error(__FILE__, __LINE__, 0,
                  _("Specific heat value is zero or not found in the xml file.\n"));

      i = *iscalt-1;
      cs_gui_properties_value("thermal_conductivity", &visls0[i]);
      /* for the Temperature, the diffusivity factor is not divided by Cp */
      if (abs(iscsth[i]) != 1)
      {
        visls0[i] = visls0[i]/result;
      }
      else
      {
        visls0[i] = visls0[i];
      }
    }

    /* User scalar
       In the interface, the user gives the diffusion coefficient, whereas in
       the solver, one sets the diffusivity, thus one need to multiply
       this coefficient by the density to remain coherent */

    for (i=0 ; i < vars->nscaus; i++) {
      if (iscavr[i] <= 0 && i != *iscalt-1) {

        if (cs_gui_strcmp(vars->model, "solid_fuels")) {
          /* Air molar mass */
          result = 0.028966;
          cs_gui_reference_initialization("mass_molar", &result);
          if (result <= 0)
            bft_error(__FILE__, __LINE__, 0,
                      _("mass molar value is zero or not found in the xml file.\n"));
          density = *p0 * result / (8.31434 *(*t0));
        }
        else
          cs_gui_properties_value("density", &density);

        if (density <= 0)
          bft_error(__FILE__, __LINE__, 0,
                    _("Density value is zero or not found in the xml file.\n"));

        coeff = visls0[i] / density ;
        cs_gui_scalar_diffusion_value(i+1, &coeff);
        visls0[i] = coeff * density;
      }
    }

#if _XML_DEBUG_
    bft_printf("==>CSSCA3\n");
    for (i=0 ; i < vars->nscaus; i++)
      bft_printf("--visls0[%i] = %f\n", i, visls0[i]);
#endif
  }
}

/*----------------------------------------------------------------------------
 * Turbulence initialization parameters.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSTINI (UREF, ALMAX)
 * *****************
 *
 * INTEGER          UREF   <--   reference velocity
 * INTEGER          ALMAX  <--   reference length
 *----------------------------------------------------------------------------*/

void CS_PROCF (cstini, CSTINI) (double *const uref,
                                double *const almax)
{
  char* length_choice = NULL;

  cs_gui_reference_initialization("velocity", uref);

  length_choice = cs_gui_reference_length_initialization_choice();

  if (cs_gui_strcmp(length_choice, "prescribed"))
    cs_gui_reference_initialization("length", almax);

  BFT_FREE(length_choice);

#if _XML_DEBUG_
  bft_printf("==>CSTINI\n");
  bft_printf("--almax = %f\n", *almax);
  bft_printf("--uref  = %f\n", *uref);
#endif
}

/*----------------------------------------------------------------------------
 * Properties array used in the calculation
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiprop, UIPROP) (const int *const irom,
                                const int *const iviscl,
                                const int *const ivisct,
                                const int *const ivisls,
                                const int *const icour,
                                const int *const ifour,
                                const int *const ismago,
                                const int *const iale,
                                const int *const icp,
                                const int *const iscalt,
                                const int *const iscavr,
                                const int *const iprtot,
                                const int *const ipppro,
                                const int *const ipproc,
                                const int *const icmome,
                                const int *const ipptx,
                                const int *const ippty,
                                const int *const ipptz,
                                const int *const ippdt,
                                const int *const ivisma,
                                const int *const idtvar,
                                const int *const ipucou,
                                const int *const iappel)
{
  int itype = 0;
  int n;
  int i = 0;
  int nbp = 5;
  char *name = NULL;

  /* Compute the new size of vars->properties_name,
     vars->properties_ipp and vars->propce */

  if (*ismago != -1 ) nbp++;

  if (*icp>0) nbp++;

  if (cs_glob_var->nscaus > 0) {
    for (i=0; i < cs_glob_var->nscaus; i++)
      if (ivisls[i] > 0 && iscavr[i] <= 0)
        nbp++;
  }

  if (!cs_gui_strcmp(cs_glob_var->model, "compressible_model"))
    nbp++;

  if (*iale) {
    cs_gui_get_ale_viscosity_type(&itype);
    if (itype == 1) {
      nbp = nbp + 3;
    } else {
      nbp++;
    }
  }

  n = cs_glob_var->nprop;

  if (*iappel == 0) {

    cs_glob_var->nprop += nbp;

    /* Fisrt step : before the third call of VARPOS in INIUSI */

    BFT_REALLOC(cs_glob_var->properties_ipp,  cs_glob_var->nprop, int);
    BFT_REALLOC(cs_glob_var->propce,  cs_glob_var->nprop, int);
    BFT_REALLOC(cs_glob_var->properties_name, cs_glob_var->nprop, char*);

    cs_glob_var->properties_ipp[n] = ipppro[ ipproc[ *irom-1 ]-1 ];
    cs_glob_var->propce[n] = ipproc[ *irom -1] -1;
    BFT_MALLOC(cs_glob_var->properties_name[n], strlen("density")+1, char);
    strcpy(cs_glob_var->properties_name[n++], "density");

    cs_glob_var->properties_ipp[n] = ipppro[ ipproc[ *iviscl-1 ]-1 ];
    cs_glob_var->propce[n] = ipproc[ *iviscl -1] -1;
    BFT_MALLOC(cs_glob_var->properties_name[n], strlen("molecular_viscosity")+1, char);
    strcpy(cs_glob_var->properties_name[n++], "molecular_viscosity");

    cs_glob_var->properties_ipp[n] = ipppro[ ipproc[ *ivisct-1 ]-1 ];
    cs_glob_var->propce[n] = ipproc[*ivisct -1] -1;
    BFT_MALLOC(cs_glob_var->properties_name[n], strlen("turb_viscosity")+1, char);
    strcpy(cs_glob_var->properties_name[n++], "turb_viscosity");

    cs_glob_var->properties_ipp[n] = ipppro[ ipproc[ *icour-1 ]-1 ];
    cs_glob_var->propce[n] = ipproc[ *icour -1] -1;
    BFT_MALLOC(cs_glob_var->properties_name[n], strlen("courant_number")+1, char);
    strcpy(cs_glob_var->properties_name[n++], "courant_number");

    cs_glob_var->properties_ipp[n] = ipppro[ ipproc[ *ifour-1 ]-1 ];
    cs_glob_var->propce[n] = ipproc[ *ifour -1] -1;
    BFT_MALLOC(cs_glob_var->properties_name[n], strlen("fourier_number")+1, char);
    strcpy(cs_glob_var->properties_name[n++], "fourier_number");

    if (*ismago != -1 ) {
      cs_glob_var->properties_ipp[n] = ipppro[ ipproc[ *ismago-1 ]-1 ];
      cs_glob_var->propce[n] = ipproc[ *ismago -1] -1;
      BFT_MALLOC(cs_glob_var->properties_name[n], strlen("smagorinsky_constant")+1, char);
      strcpy(cs_glob_var->properties_name[n++], "smagorinsky_constant");
    }

    if (*icp > 0) {
      cs_glob_var->properties_ipp[n] = ipppro[ ipproc[ *icp-1 ]-1 ];
      cs_glob_var->propce[n] = ipproc[ *icp -1] -1;
      BFT_MALLOC(cs_glob_var->properties_name[n], strlen("specific_heat")+1, char);
      strcpy(cs_glob_var->properties_name[n++], "specific_heat");
    }

    if (!cs_gui_strcmp(cs_glob_var->model, "compressible_model")) {
      cs_glob_var->properties_ipp[n] = ipppro[ ipproc[ *iprtot-1 ]-1 ];
      cs_glob_var->propce[n] = ipproc[ *iprtot -1] -1;
      BFT_MALLOC(cs_glob_var->properties_name[n], strlen("total_pressure")+1, char);
      strcpy(cs_glob_var->properties_name[n++], "total_pressure");
    }

    if (*iale) {
      cs_glob_var->properties_ipp[n] = ipppro[ ipproc[ ivisma[0]-1 ]-1 ];
      cs_glob_var->propce[n] = ipproc[ivisma[0] -1] -1;
      BFT_MALLOC(cs_glob_var->properties_name[n], strlen("mesh_viscosity_1")+1, char);
      strcpy(cs_glob_var->properties_name[n++], "mesh_viscosity_1");

      if (itype == 1) {
        cs_glob_var->properties_ipp[n] = ipppro[ ipproc[ ivisma[1]-1 ]-1 ];
        cs_glob_var->propce[n] = ipproc[ivisma[1] -1] -1;
        BFT_MALLOC(cs_glob_var->properties_name[n], strlen("mesh_viscosity_2")+1, char);
        strcpy(cs_glob_var->properties_name[n++], "mesh_viscosity_2");

        cs_glob_var->properties_ipp[n] = ipppro[ ipproc[ ivisma[2]-1 ]-1 ];
        cs_glob_var->propce[n] = ipproc[ivisma[2] -1] -1;
        BFT_MALLOC(cs_glob_var->properties_name[n], strlen("mesh_viscosity_3")+1, char);
        strcpy(cs_glob_var->properties_name[n++], "mesh_viscosity_3");
      }

    }

    /* scalar diffusivity */

    if (cs_glob_var->nscaus > 0) {

      /* search lenght of first character of scalar property's suffixe : '_' */
      for (i=0; i < cs_glob_var->nscaus; i++) {

        if (iscavr[i] <= 0 && ivisls[i] > 0) {

          cs_glob_var->properties_ipp[n] = ipppro[ ipproc[ ivisls[i]-1 ]-1 ];
          cs_glob_var->propce[n] = ipproc[ *ivisls -1] -1;

          if (*iscalt == i+1) {
            BFT_MALLOC(cs_glob_var->properties_name[n], strlen("thermal_conductivity")+1, char);
            strcpy(cs_glob_var->properties_name[n++], "thermal_conductivity");
          } else {
            name = _scalar_diffusion_coefficient_name(i);
            BFT_MALLOC(cs_glob_var->properties_name[n], strlen(name)+1, char);
            strcpy(cs_glob_var->properties_name[n++], name);
            BFT_FREE(name);
          }
        }
      }
    }

  } else {

    /* Second step : before the fourth call of VARPOS in INIUSI */

    if (*idtvar == 1 || *idtvar == 2)
      cs_glob_var->nprop += 1;
    if (*ipucou == 1)
      cs_glob_var->nprop += 3;
    cs_glob_var->nprop += cs_glob_var->ntimaver;

    BFT_REALLOC(cs_glob_var->properties_ipp,  cs_glob_var->nprop, int);
    BFT_REALLOC(cs_glob_var->propce,  cs_glob_var->nprop, int);
    BFT_REALLOC(cs_glob_var->properties_name, cs_glob_var->nprop, char*);

    if (*idtvar == 1 || *idtvar == 2) {
      cs_glob_var->properties_ipp[n] = *ippdt;
      cs_glob_var->propce[n] = -1;
      BFT_MALLOC(cs_glob_var->properties_name[n], strlen("local_time_step")+1, char);
      strcpy(cs_glob_var->properties_name[n++], "local_time_step");
    }

    if (*ipucou == 1) {
      cs_glob_var->properties_ipp[n] = *ipptx;
      cs_glob_var->propce[n] = -1;
      BFT_MALLOC(cs_glob_var->properties_name[n], strlen("weight_matrix_X")+1, char);
      strcpy(cs_glob_var->properties_name[n++], "weight_matrix_X");

      cs_glob_var->properties_ipp[n] = *ippty;
      cs_glob_var->propce[n] = -1;
      BFT_MALLOC(cs_glob_var->properties_name[n], strlen("weight_matrix_Y")+1, char);
      strcpy(cs_glob_var->properties_name[n++], "weight_matrix_Y");

      cs_glob_var->properties_ipp[n] = *ipptz;
      cs_glob_var->propce[n] = -1;
      BFT_MALLOC(cs_glob_var->properties_name[n], strlen("weight_matrix_Z")+1, char);
      strcpy(cs_glob_var->properties_name[n++], "weight_matrix_Z");
    }

    for (i=0; i < cs_glob_var->ntimaver; i++) {
      cs_glob_var->properties_ipp[n] = ipppro[ ipproc[ icmome[i]-1 ]-1 ];
      cs_glob_var->propce[n] = ipproc[icmome[i] -1] -1;
      name = _get_time_average_label(i+1);
      BFT_MALLOC(cs_glob_var->properties_name[n], strlen(name)+1, char);
      strcpy(cs_glob_var->properties_name[n++], name);
      BFT_FREE(name);
    }
  }

  if (n != cs_glob_var->nprop)
    bft_error(__FILE__, __LINE__, 0,
              _("number of properties is not correct: %i instead of: %i\n"),
              n, cs_glob_var->nprop);

#if _XML_DEBUG_
  bft_printf("==>UIPROP %i\n",*iappel);
  bft_printf("-->nombre de proprietes = %i\n", cs_glob_var->nprop);
  for (i=0; i < cs_glob_var->nprop; i++) {
    bft_printf("-->properties_ipp[%i]: %i propce[%i]: %i "
        "properties_name[%i]: %s\n",
        i, cs_glob_var->properties_ipp[i],
        i, cs_glob_var->propce[i],
        i, cs_glob_var->properties_name[i]);
  }
#endif
}

/*----------------------------------------------------------------------------
 * Temporal averaging treatment
 *----------------------------------------------------------------------------*/

void CS_PROCF (uimoyt, UIMOYT) (const int *const ndgmox,
                                      int *const ntdmom,
                                      int *const imoold,
                                      int *const idfmom)
{
  int imom = 0;
  int isuite = 0;
  int i, j, n;
  char *name = NULL;

  cs_glob_var->ntimaver
    = cs_gui_get_tag_number("/analysis_control/time_averages/time_average", 1);

  /* for each average */
  for (i=0; i < cs_glob_var->ntimaver; i++) {

    imom = i + 1;

    _get_time_average_data(imom, "time_step_start", &ntdmom[i]);

    /* test on isuite */
    cs_gui_restart_parameters_status("restart", &isuite);

    if (isuite != 0) {
      _get_time_average_data(imom, "restart_from_time_average", &imoold[i]);
      if (imoold[i] == imom) imoold[i] = -2;
    }

    for (n=0; n < _get_time_average_n_variables(imom); n++) {

      name = _get_time_average_variable_name(imom, n + 1);

      for (j=0; j < cs_glob_var->nvar; j++) {
        if (cs_gui_strcmp(name,  cs_glob_var->name[j])) {
          idfmom[(imom-1)*(*ndgmox) + n] = cs_glob_var->rtp[j] +1;
        }
      }

      for (j=0 ; j < cs_glob_var->nprop; j++) {
        if (cs_gui_strcmp(name, cs_glob_var->properties_name[j]))
          idfmom[(imom-1)*(*ndgmox) + n] = -(cs_glob_var->propce[j]);
      }

      BFT_FREE(name);
    }
  }
#if _XML_DEBUG_
  bft_printf("==>UIMOYT\n");
  for (i=0; i < cs_glob_var->ntimaver; i++) {
    bft_printf("-->ntdmom =  %i\n", ntdmom[i]);
  }
#endif
}

/*----------------------------------------------------------------------------
 * Copy variable name from Fortran to C
 *----------------------------------------------------------------------------*/

void CS_PROCF(fcnmva, FCNMVA)
  (
    const char          *const fstr,    /* --> Fortran string */
    int                 *const len,     /* --> String Length  */
    int                 *const var_id   /* --> Variable Id (1 to n) */
    CS_ARGF_SUPP_CHAINE
  )
{
  int i, i1, i2, l;
  char *cstr = NULL;

  assert(*var_id > 0);

  /* Resize array if necessary */

  if (*var_id > cs_glob_label->_cs_gui_max_vars) {

    if (cs_glob_label->_cs_gui_max_vars == 0)
      cs_glob_label->_cs_gui_max_vars = 16;

    while (cs_glob_label->_cs_gui_max_vars <= *var_id)
      cs_glob_label->_cs_gui_max_vars *= 2;

    BFT_REALLOC(cs_glob_label->_cs_gui_var_name, cs_glob_label->_cs_gui_max_vars, char *);
    for (i = cs_glob_label->_cs_gui_last_var; i < cs_glob_label->_cs_gui_max_vars; i++)
      cs_glob_label->_cs_gui_var_name[i] = NULL;
  }

  /* Compute string length (removing start or end blanks) */

  for (i1 = 0;
      i1 < *len && (fstr[i1] == ' ' || fstr[i1] == '\t');
      i1++);

  for (i2 = *len - 1;
      i2 > i1 && (fstr[i2] == ' ' || fstr[i2] == '\t');
      i2--);

  l = i2 - i1 + 1;

  /* Should be called once per variable only */
  assert(cs_glob_label->_cs_gui_var_name[*var_id - 1] == NULL);

  if (l > 0) {

    /* Allocate and copy */
    BFT_MALLOC(cstr, l + 1, char);

    for (i = 0 ; i < l ; i++, i1++)
      cstr[i] = fstr[i1];

    cstr[l] = '\0';

    cs_glob_label->_cs_gui_var_name[*var_id - 1] = cstr;

  }

  /* Update variable counter */
  cs_glob_label->_cs_gui_last_var = *var_id;

}

/*----------------------------------------------------------------------------
 * Copy variable name from C to Fortran
 *----------------------------------------------------------------------------*/

void CS_PROCF(cfnmva, CFNMVA)
  (
    char          *const fstr,    /* --> Fortran string */
    int           *const len,     /* --> String Length  */
    int           *const var_id   /* --> Variable Id (1 to n) */
    CS_ARGF_SUPP_CHAINE
  )
{
  int i;
  int l = 0;
  char *cstr = NULL;

  /* Check that variable name was set */

  if (*var_id < 1 || *var_id > cs_glob_label->_cs_gui_last_var)
    bft_error(__FILE__, __LINE__, 0,
              _("Name of variable %i was never set.\n"), *var_id);

  /* Copy string */

  cstr = cs_glob_label->_cs_gui_var_name[*var_id - 1];

  if (cstr != NULL) {

    /* Compute string length (removing start or end blanks) */

    l = strlen(cstr);
    if (l > *len)
      l = *len;

    for (i = 0; i < l; i++)
      fstr[i] = cstr[i];

  }

  /* Pad with blanks if necessary */

  for (i = l; i < *len; i++)
    fstr[i] = ' ';
}

/*----------------------------------------------------------------------------
 * Clean memory for fortran name of variables
 *----------------------------------------------------------------------------*/

void CS_PROCF(nvamem, NVAMEM) (void)
{
  int i;
#if _XML_DEBUG_
  bft_printf("==>NVAMEM\n");
  for (i = 0; i < cs_glob_label->_cs_gui_max_vars; i++)
    if (cs_glob_label->_cs_gui_var_name[i])
      bft_printf("-->label[%i] = %s\n", i, cs_glob_label->_cs_gui_var_name[i]);
#endif

  for (i = 0; i < cs_glob_label->_cs_gui_max_vars; i++)
    BFT_FREE(cs_glob_label->_cs_gui_var_name[i]);

  BFT_FREE(cs_glob_label->_cs_gui_var_name);

  cs_glob_label->_cs_gui_max_vars = 0;
  cs_glob_label->_cs_gui_last_var = 0;
}

/*-----------------------------------------------------------------------------
 * Return the list of cells describing a given zone.
 *
 * parameters:
 *   zone_id   -->  volume zone id
 *   ncelet    -->  number of cells with halo
 *   faces     <--  number of selected cells
 *----------------------------------------------------------------------------*/

static int*
cs_gui_get_cells_list(const char *zone_id,
                      const int   ncelet,
                            int  *cells )
{
  int  c_id         = 0;
  int  *cells_list  = NULL;
  char *description = NULL;

  description = cs_gui_volumic_zone_localization(zone_id);

  /* build list of cells */
  BFT_MALLOC(cells_list, ncelet, int);

  c_id = fvm_selector_get_list(cs_glob_mesh->select_cells,
                               description,
                               cells,
                               cells_list);

  if (fvm_selector_n_missing(cs_glob_mesh->select_cells, c_id) > 0)
  {
    const char *missing
      = fvm_selector_get_missing(cs_glob_mesh->select_cells, c_id, 0);
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("The group or attribute \"%s\" in the selection\n"
                 "criteria:\n"
                 "\"%s\"\n"
                 " does not correspond to any cell.\n"),
                 missing, description);
  }
  BFT_FREE(description);
  return cells_list;
}

/*----------------------------------------------------------------------------
 * User momentum source terms.
 *
 * Fortran Interface:
 *
 * subroutine uitsnv (vel, tsexp, tsimp)
 * *****************
 *
 * double precision vel      -->  fluid velocity
 * double precision tsexp    <--  explicit source terms
 * double precision tsimp    <--  implicit source terms
 *----------------------------------------------------------------------------*/

void CS_PROCF(uitsnv, UITSNV)(const cs_real_3_t *restrict vel,
                              cs_real_3_t       *restrict tsexp,
                              cs_real_33_t      *restrict tsimp)
{
  const int n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  const cs_real_t *restrict cell_vol = cs_glob_mesh_quantities->cell_vol;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->cell_cen;


  int i, icel, iel;
  int zones = 0;
  int cells = 0;
  int *cells_list = NULL;
  double dSudu, dSudv, dSudw;
  double dSvdu, dSvdv, dSvdw;
  double dSwdu, dSwdv, dSwdw;
  char *path = NULL;
  char *status = NULL;
  char *zone_id = NULL;
  char *formula = NULL;
  char *labelU = NULL;
  char *labelV = NULL;
  char *labelW = NULL;

  mei_tree_t *ev_formula  = NULL;

  /* number of volumic zone */

  zones = cs_gui_get_tag_number("/solution_domain/volumic_conditions/zone\n", 1);

#if _XML_DEBUG_
  bft_printf("==>UITSNV\n");
#endif

  for (i=1; i < zones+1; i++) {
    /* momentum source term */
    path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 2, "solution_domain", "volumic_conditions");
    cs_xpath_add_element_num(&path, "zone", i);
    cs_xpath_add_attribute(&path, "momentum_source_term");
    status = cs_gui_get_attribute_value(path);
    BFT_FREE(path);

    if (cs_gui_strcmp(status, "on")) {
      zone_id = cs_gui_volumic_zone_id(i);
      cells_list = cs_gui_get_cells_list(zone_id, n_cells_ext, &cells);

      path = cs_xpath_init_path();
      cs_xpath_add_elements(&path, 1, "thermophysical_models");
      cs_xpath_add_elements(&path, 1, "source_terms");
      cs_xpath_add_elements(&path, 1, "momentum_formula");
      cs_xpath_add_test_attribute(&path, "zone_id",zone_id);
      cs_xpath_add_function_text(&path);
      formula = cs_gui_get_text_value(path);
      BFT_FREE(path);
      if (formula != NULL) {
        ev_formula = mei_tree_new(formula);
        mei_tree_insert(ev_formula,"x",0.0);
        mei_tree_insert(ev_formula,"y",0.0);
        mei_tree_insert(ev_formula,"z",0.0);
        labelU = cs_gui_variable_label("velocity_U");
        mei_tree_insert(ev_formula, labelU, 0.0);
        labelV = cs_gui_variable_label("velocity_V");
        mei_tree_insert(ev_formula, labelV, 0.0);
        labelW = cs_gui_variable_label("velocity_W");
        mei_tree_insert(ev_formula, labelW, 0.0);
        /* try to build the interpreter */
        if (mei_tree_builder(ev_formula))
          bft_error(__FILE__, __LINE__, 0,
                    _("Error: can not interpret expression: %s\n %i"),
                    ev_formula->string, mei_tree_builder(ev_formula));
        const char *symbols[] = {"Su","Sv","Sw",
                                 "dSudu","dSudv","dSudw",
                                 "dSvdu","dSvdv","dSvdw",
                                 "dSwdu","dSwdv","dSwdw"};
        if (mei_tree_find_symbols(ev_formula, 12, symbols))
          bft_error(__FILE__, __LINE__, 0,
                    _("Error: can not find the required symbol: %s\n%s\n%s\n%s\n"),
                    "Su, Sv, Sw,",
                    "dSudu, dSudv, dSudw,",
                    "dSvdu, dSvdv, dSvdw,",
                    "dSwdu, dSwdv or dSwdw");
        for (icel = 0; icel < cells; icel++) {
          iel = cells_list[icel]-1;
          mei_tree_insert(ev_formula, "x", cell_cen[iel][0]);
          mei_tree_insert(ev_formula, "y", cell_cen[iel][1]);
          mei_tree_insert(ev_formula, "z", cell_cen[iel][2]);
          mei_tree_insert(ev_formula, labelU, vel[iel][0]);
          mei_tree_insert(ev_formula, labelV, vel[iel][1]);
          mei_tree_insert(ev_formula, labelW, vel[iel][2]);
          mei_evaluate(ev_formula);

          dSudu = mei_tree_lookup(ev_formula,"dSudu");
          dSudv = mei_tree_lookup(ev_formula,"dSudv");
          dSudw = mei_tree_lookup(ev_formula,"dSudw");
          dSvdu = mei_tree_lookup(ev_formula,"dSvdu");
          dSvdv = mei_tree_lookup(ev_formula,"dSvdv");
          dSvdw = mei_tree_lookup(ev_formula,"dSvdw");
          dSwdu = mei_tree_lookup(ev_formula,"dSwdu");
          dSwdv = mei_tree_lookup(ev_formula,"dSwdv");
          dSwdw = mei_tree_lookup(ev_formula,"dSwdw");

          tsimp[iel][0][0] = cell_vol[iel]*dSudu;
          tsimp[iel][0][1] = cell_vol[iel]*dSudv;
          tsimp[iel][0][2] = cell_vol[iel]*dSudw;
          tsimp[iel][1][0] = cell_vol[iel]*dSvdu;
          tsimp[iel][1][1] = cell_vol[iel]*dSvdv;
          tsimp[iel][1][2] = cell_vol[iel]*dSvdw;
          tsimp[iel][2][0] = cell_vol[iel]*dSwdu;
          tsimp[iel][2][1] = cell_vol[iel]*dSwdv;
          tsimp[iel][2][2] = cell_vol[iel]*dSwdw;

          tsexp[iel][0] = mei_tree_lookup(ev_formula,"Su")
                        - ( dSudu*vel[iel][0]
                          + dSudv*vel[iel][1]
                          + dSudw*vel[iel][2]
                          );
          tsexp[iel][0] *= cell_vol[iel];
          tsexp[iel][1] = mei_tree_lookup(ev_formula,"Sv")
                        - ( dSvdu*vel[iel][0]
                          + dSvdv*vel[iel][1]
                          + dSvdw*vel[iel][2]
                          );
          tsexp[iel][1] *= cell_vol[iel];
          tsexp[iel][2] = mei_tree_lookup(ev_formula,"Sw")
                        - ( dSwdu*vel[iel][0]
                          + dSwdv*vel[iel][1]
                          + dSwdw*vel[iel][2]
                          );
          tsexp[iel][2] *= cell_vol[iel];
        }
        mei_tree_destroy(ev_formula);
      }
      BFT_FREE(cells_list);
      BFT_FREE(zone_id);
    }

    BFT_FREE(status);
  }
}

/*----------------------------------------------------------------------------
 * User scalar source terms.
 *
 * Fortran Interface:
 *
 * subroutine uitssc (iscal, pvar, tsexp, tsimp)
 * *****************
 *
 * integer          iscal    --> index of the corresponding scalar
 * double precision pvar     -->  scalar
 * double precision tsexp    <--  explicit source terms
 * double precision tsimp    <--  implicit source terms
 *----------------------------------------------------------------------------*/

void CS_PROCF(uitssc, UITSSC)(const int                  *iscal,
                              const cs_real_t   *restrict pvar,
                              cs_real_t         *restrict tsexp,
                              cs_real_t         *restrict tsimp)
{
  const int n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  const cs_real_t *restrict cell_vol = cs_glob_mesh_quantities->cell_vol;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->cell_cen;

  int i, icel, iel;
  int zones = 0;
  int cells = 0;
  int *cells_list = NULL;
  double dS;
  char *path = NULL;
  char *status = NULL;
  char *zone_id = NULL;
  char *formula = NULL;

  mei_tree_t *ev_formula  = NULL;

  cs_var_t  *vars = cs_glob_var;

  /* number of volumic zone */

  zones = cs_gui_get_tag_number("/solution_domain/volumic_conditions/zone\n", 1);

#if _XML_DEBUG_
  bft_printf("==>UITSSC\n");
#endif

  for (i=1; i < zones+1; i++) {

    /* species source term */
    path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 2, "solution_domain", "volumic_conditions");
    cs_xpath_add_element_num(&path, "zone", i);
    cs_xpath_add_attribute(&path, "scalar_source_term");
    status = cs_gui_get_attribute_value(path);
    BFT_FREE(path);

    if (cs_gui_strcmp(status, "on")) {
      zone_id = cs_gui_volumic_zone_id(i);
      cells_list = cs_gui_get_cells_list(zone_id, n_cells_ext, &cells);

      path = cs_xpath_init_path();
      cs_xpath_add_elements(&path, 3,
                            "thermophysical_models",
                            "source_terms",
                            "scalar_formula");
      cs_xpath_add_test_attribute(&path, "label", vars->label[*iscal-1]);
      cs_xpath_add_test_attribute(&path, "zone_id", zone_id);
      cs_xpath_add_function_text(&path);
      formula = cs_gui_get_text_value(path);
      BFT_FREE(path);

      if (formula != NULL) {
        ev_formula = mei_tree_new(formula);
        mei_tree_insert(ev_formula,"x",0.);
        mei_tree_insert(ev_formula,"y",0.);
        mei_tree_insert(ev_formula,"z",0.);
        mei_tree_insert(ev_formula, vars->label[*iscal-1], 0.0);
        /* try to build the interpreter */
        if (mei_tree_builder(ev_formula))
          bft_error(__FILE__, __LINE__, 0,
                    _("Error: can not interpret expression: %s\n %i"),
                    ev_formula->string, mei_tree_builder(ev_formula));

        const char *symbols[] = {"S","dS"};
        if (mei_tree_find_symbols(ev_formula, 2, symbols))
          bft_error(__FILE__, __LINE__, 0,
                    _("Error: can not find the required symbol: %s\n"), "S or dS");

        for (icel = 0; icel < cells; icel++) {
          iel = cells_list[icel]-1;
          mei_tree_insert(ev_formula, "x", cell_cen[iel][0]);
          mei_tree_insert(ev_formula, "y", cell_cen[iel][1]);
          mei_tree_insert(ev_formula, "z", cell_cen[iel][2]);
          mei_tree_insert(ev_formula, vars->label[*iscal-1], pvar[iel]);
          mei_evaluate(ev_formula);
          dS = mei_tree_lookup(ev_formula,"dS");
          tsimp[iel] = cell_vol[iel]*dS;
          tsexp[iel] = mei_tree_lookup(ev_formula,"S") - dS*pvar[iel];
          tsexp[iel] *= cell_vol[iel];
        }
        mei_tree_destroy(ev_formula);
      }
      BFT_FREE(cells_list);
      BFT_FREE(zone_id);
    }
    BFT_FREE(status);
  }
}

/*----------------------------------------------------------------------------
 * Thermal scalar source terms.
 *
 * Fortran Interface:
 *
 * subroutine uitsth (iscal, pvar, tsexp, tsimp)
 * *****************
 *
 * integer          iscal    --> index of the corresponding scalar
 * double precision pvar     -->  scalar
 * double precision tsexp    <--  explicit source terms
 * double precision tsimp    <--  implicit source terms
 *----------------------------------------------------------------------------*/

void CS_PROCF(uitsth, UITSTH)(const int                  *iscal,
                              const cs_real_t   *restrict pvar,
                              cs_real_t         *restrict tsexp,
                              cs_real_t         *restrict tsimp)
{
  const int n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  const cs_real_t *restrict cell_vol = cs_glob_mesh_quantities->cell_vol;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->cell_cen;

  int i, icel, iel;
  int zones = 0;
  int cells = 0;
  int *cells_list = NULL;
  double dS;
  char *path = NULL;
  char *status = NULL;
  char *zone_id = NULL;
  char *formula = NULL;

  mei_tree_t *ev_formula  = NULL;

  cs_var_t  *vars = cs_glob_var;

  /* number of volumic zone */

  zones = cs_gui_get_tag_number("/solution_domain/volumic_conditions/zone\n", 1);

#if _XML_DEBUG_
  bft_printf("==>UITSSC\n");
#endif

  for (i=1; i < zones+1; i++) {

    /* species source term */
    path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 2, "solution_domain", "volumic_conditions");
    cs_xpath_add_element_num(&path, "zone", i);
    cs_xpath_add_attribute(&path, "thermal_source_term");
    status = cs_gui_get_attribute_value(path);
    BFT_FREE(path);

    if (cs_gui_strcmp(status, "on")) {
      zone_id = cs_gui_volumic_zone_id(i);
      cells_list = cs_gui_get_cells_list(zone_id, n_cells_ext, &cells);

      path = cs_xpath_init_path();
      cs_xpath_add_elements(&path, 3,
                            "thermophysical_models",
                            "source_terms",
                            "thermal_formula");
      cs_xpath_add_test_attribute(&path, "label", vars->label[*iscal-1]);
      cs_xpath_add_test_attribute(&path, "zone_id", zone_id);
      cs_xpath_add_function_text(&path);
      formula = cs_gui_get_text_value(path);
      BFT_FREE(path);

      if (formula != NULL) {
        ev_formula = mei_tree_new(formula);
        mei_tree_insert(ev_formula,"x",0.);
        mei_tree_insert(ev_formula,"y",0.);
        mei_tree_insert(ev_formula,"z",0.);
        mei_tree_insert(ev_formula, vars->label[*iscal-1], 0.0);
        /* try to build the interpreter */
        if (mei_tree_builder(ev_formula))
          bft_error(__FILE__, __LINE__, 0,
                    _("Error: can not interpret expression: %s\n %i"),
                    ev_formula->string, mei_tree_builder(ev_formula));

        const char *symbols[] = {"S","dS"};
        if (mei_tree_find_symbols(ev_formula, 2, symbols))
          bft_error(__FILE__, __LINE__, 0,
                    _("Error: can not find the required symbol: %s\n"), "S or dS");

        for (icel = 0; icel < cells; icel++) {
          iel = cells_list[icel]-1;
          mei_tree_insert(ev_formula, "x", cell_cen[iel][0]);
          mei_tree_insert(ev_formula, "y", cell_cen[iel][1]);
          mei_tree_insert(ev_formula, "z", cell_cen[iel][2]);
          mei_tree_insert(ev_formula, vars->label[*iscal-1], pvar[iel]);
          mei_evaluate(ev_formula);
          dS = mei_tree_lookup(ev_formula,"dS");
          tsimp[iel] = cell_vol[iel]*dS;
          tsexp[iel] = mei_tree_lookup(ev_formula,"S") - dS*pvar[iel];
          tsexp[iel] *= cell_vol[iel];
        }
        mei_tree_destroy(ev_formula);
      }
      BFT_FREE(cells_list);
      BFT_FREE(zone_id);
    }
    BFT_FREE(status);
  }
}

/*-----------------------------------------------------------------------------
 * Initialize mei tree and check for symbols existence
 *
 * parameters:
 *   formula        -->  mei formula
 *   symbols        -->  array of symbol to check
 *   symbol_size    -->  number of symbol in symbols
 *----------------------------------------------------------------------------*/

static mei_tree_t *_init_mei_tree(const char *formula,
        const char *symbols)
{

  /* return an empty interpreter */
  mei_tree_t *tree = mei_tree_new(formula);

  /* add commun variables */
  mei_tree_insert(tree, "x",    0.0);
  mei_tree_insert(tree, "y",    0.0);
  mei_tree_insert(tree, "z",    0.0);

  /* try to build the interpreter */
  if (mei_tree_builder(tree))
    bft_error(__FILE__, __LINE__, 0,
              _("Error: can not interpret expression: %s\n"), tree->string);
  /* check for symbols */
  if (mei_tree_find_symbol(tree, symbols))
    bft_error(__FILE__, __LINE__, 0,
              _("Error: can not find the required symbol: %s\n"), symbols);

  return tree;
}

/*----------------------------------------------------------------------------
 * Variables and user scalars initialization.
 *
 * Fortran Interface:
 *
 * subroutine uiiniv (ncelet, isuite, isca, iscold, rtp)
 * *****************
 *
 * integer          ncelet   -->  number of cells with halo
 * integer          isuite   -->  restart indicator
 * integer          isca     -->  indirection array for scalar number
 * integer          iscold   -->  scalar number for restart
 * integer          iccfth   -->  type of initialisation(compressible model)
 * integer          ipr      -->  rtp index for pressure
 * integer          ipcrom   -->  propce index for density
 * integer          itempk   -->  rtp index for temperature (in K)
 * integer          ienerg   -->  rtp index for energy total
 * DOUBLE PRECISION RO0      -->  value of density if IROVAR=0
 * DOUBLE PRECISION CP0      -->  value of specific heat if ICP=0
 * DOUBLE PRECISION VISCL0   -->  value of viscosity if IVIVAR=0
 * DOUBLE PRECISION VISLS0   -->  value of reference molecular diffusivity
 * DOUBLE PRECISION UREF     -->  value of reference velocity
 * DOUBLE PRECISION ALMAX    -->  value of reference length
 * DOUBLE PRECISION XYZCEN   -->  cell's gravity center
 * double precision rtp     <--   variables and scalars array
 * double precision propce  <--   physical properties array
 *----------------------------------------------------------------------------*/

void CS_PROCF(uiiniv, UIINIV)(const int          *ncelet,
                              const int          *isuite,
                              const int          *isca,
                              const int          *iscold,
                                    int          *iccfth,
                              const int *const    ipr,
                              const int *const    ipcrom,
                              const int *const    itempk,
                              const int *const    ienerg,
                              const cs_real_t    *ro0,
                              const cs_real_t    *cp0,
                              const cs_real_t    *viscl0,
                              const cs_real_t    *uref,
                              const cs_real_t    *almax,
                              const double *const xyzcen,
                                    double        rtp[],
                                    double        propce[])
{
  /* Coal combustion: the initialization of the model scalar are not given */

  int i, j, icel, iel;
  int zones            = 0;
  int cells            = 0;
  int ccfth            = 0;
  int *cells_list      = NULL;
  char *choice         = NULL;
  char *buff           = NULL;
  char *path           = NULL;
  char *path1          = NULL;
  char *path_velocity  = NULL;
  char *path_turb      = NULL;
  char *path_sca       = NULL;
  char *path_meteo     = NULL;
  char *status         = NULL;
  char *zone_id        = NULL;
  char *formula_uvw    = NULL;
  char *formula_turb   = NULL;
  char *formula_sca    = NULL;
  char *formula_meteo  = NULL;
  char *formula        = NULL;
  char *model          = NULL;

  mei_tree_t *ev_formula_uvw   = NULL;
  mei_tree_t *ev_formula_turb  = NULL;
  mei_tree_t *ev_formula_sca   = NULL;
  mei_tree_t *ev_formula_meteo = NULL;
  mei_tree_t *ev_formula       = NULL;

  cs_var_t  *vars = cs_glob_var;

  /* number of volumic zone */

  zones
    = cs_gui_get_tag_number("/solution_domain/volumic_conditions/zone", 1);

#if _XML_DEBUG_
  bft_printf("==>UIINIV\n");
#endif

  for (i=1; i < zones+1; i++) {

    path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 2, "solution_domain", "volumic_conditions");
    cs_xpath_add_element_num(&path, "zone", i);
    cs_xpath_add_attribute(&path, "initialization");
    status = cs_gui_get_attribute_value(path);
    BFT_FREE(path);

    if (cs_gui_strcmp(status, "on")) {

      zone_id = cs_gui_volumic_zone_id(i);
      cells_list = cs_gui_get_cells_list(zone_id, *ncelet, &cells);

      if (*isuite == 0) {
        path_velocity = cs_xpath_init_path();
        cs_xpath_add_elements(&path_velocity, 4,
                              "thermophysical_models",
                              "velocity_pressure",
                              "initialization",
                              "formula");
        cs_xpath_add_test_attribute(&path_velocity, "zone_id", zone_id);
        cs_xpath_add_function_text(&path_velocity);
        formula_uvw = cs_gui_get_text_value(path_velocity);

        if (formula_uvw != NULL) {
          ev_formula_uvw = mei_tree_new(formula_uvw);
          mei_tree_insert(ev_formula_uvw,"x",0.0);
          mei_tree_insert(ev_formula_uvw,"y",0.0);
          mei_tree_insert(ev_formula_uvw,"z",0.0);
          /* try to build the interpreter */
          if (mei_tree_builder(ev_formula_uvw))
            bft_error(__FILE__, __LINE__, 0,
                      _("Error: can not interpret expression: %s\n %i"),
                      ev_formula_uvw->string, mei_tree_builder(ev_formula_uvw));

          const char *symbols_uvw[] = {"u", "v", "w"};
          if (mei_tree_find_symbols(ev_formula_uvw, 3, symbols_uvw))
            bft_error(__FILE__, __LINE__, 0,
                      _("Error: can not find the required symbol: %s\n"),
                      "u, v ou w");

          for (icel = 0; icel < cells; icel++) {
            iel = cells_list[icel] - 1;
            mei_tree_insert(ev_formula_uvw, "x", xyzcen[3 * iel + 0]);
            mei_tree_insert(ev_formula_uvw, "y", xyzcen[3 * iel + 1]);
            mei_tree_insert(ev_formula_uvw, "z", xyzcen[3 * iel + 2]);
            mei_evaluate(ev_formula_uvw);
            rtp[vars->rtp[1]*(*ncelet) + iel] = mei_tree_lookup(ev_formula_uvw, "u");
            rtp[vars->rtp[2]*(*ncelet) + iel] = mei_tree_lookup(ev_formula_uvw, "v");
            rtp[vars->rtp[3]*(*ncelet) + iel] = mei_tree_lookup(ev_formula_uvw, "w");
          }
          mei_tree_destroy(ev_formula_uvw);
        }
        else {
          for (j=1; j < 4; j++) {
            for (icel = 0; icel < cells; icel++) {
              iel = cells_list[icel]-1;
              rtp[vars->rtp[j] * (*ncelet) + iel] = 0.0;
            }
          }
        }
        BFT_FREE(formula_uvw);
        BFT_FREE(path_velocity);

        /* Turbulence variables initialization */
        choice = cs_gui_turbulence_initialization_choice(zone_id);

        if (cs_gui_strcmp(choice, "formula")) {

          path_turb = cs_xpath_init_path();
          cs_xpath_add_elements(&path_turb, 3,
                                "thermophysical_models",
                                "turbulence",
                                "initialization");
          cs_xpath_add_test_attribute(&path_turb, "zone_id", zone_id);
          cs_xpath_add_element(&path, "formula");
          cs_xpath_add_function_text(&path_turb);
          formula_turb = cs_gui_get_text_value(path_turb);
          BFT_FREE(path_turb);

          if (formula_turb != NULL) {
            ev_formula_turb = mei_tree_new(formula_turb);
            mei_tree_insert(ev_formula_turb, "rho0", *ro0);
            mei_tree_insert(ev_formula_turb, "mu0", *viscl0);
            mei_tree_insert(ev_formula_turb, "cp0", *cp0);
            mei_tree_insert(ev_formula_turb, "uref", *uref);
            mei_tree_insert(ev_formula_turb, "almax", *almax);
            mei_tree_insert(ev_formula_turb, "x", 0.0);
            mei_tree_insert(ev_formula_turb, "y", 0.0);
            mei_tree_insert(ev_formula_turb, "z", 0.0);

            /* try to build the interpreter */

            if (mei_tree_builder(ev_formula_turb))
              bft_error(__FILE__, __LINE__, 0,
                        _("Error: can not interpret expression: %s\n %i"),
                        ev_formula_turb->string, mei_tree_builder(ev_formula_turb));

            model = cs_gui_get_thermophysical_model("turbulence");
            if (model == NULL) return;

            if (cs_gui_strcmp(model, "k-epsilon") ||
                cs_gui_strcmp(model, "k-epsilon-PL")) {

              const char *symbols[] = {"k","eps"};
              if (mei_tree_find_symbols(ev_formula_turb, 2, symbols))
                bft_error(__FILE__, __LINE__, 0,
                          _("Error: can not find the required symbol: %s\n"),
                          "k or eps");

              for (icel = 0; icel < cells; icel++) {
                iel = cells_list[icel] - 1;
                mei_tree_insert(ev_formula_turb, "x", xyzcen[3 * iel + 0]);
                mei_tree_insert(ev_formula_turb, "y", xyzcen[3 * iel + 1]);
                mei_tree_insert(ev_formula_turb, "z", xyzcen[3 * iel + 2]);
                mei_evaluate(ev_formula_turb);
                rtp[vars->rtp[4] * (*ncelet) + iel] = mei_tree_lookup(ev_formula_turb, "k");
                rtp[vars->rtp[5] * (*ncelet) + iel] = mei_tree_lookup(ev_formula_turb, "eps");
              }
            }

            else if (cs_gui_strcmp(model, "Rij-epsilon") || cs_gui_strcmp(model, "Rij-SSG")) {
              const char *symbols[] = {"r11", "r22", "r133", "r12", "r13", "r23", "eps"};
              if (mei_tree_find_symbols(ev_formula_turb, 7, symbols))
                bft_error(__FILE__, __LINE__, 0,
                          _("Error: can not find the required symbol: %s\n"),
                          "r11, r22, r33, r12, r13, r23 or eps");

              for (icel = 0; icel < cells; icel++) {
                iel = cells_list[icel]-1;
                mei_tree_insert(ev_formula_turb, "x", xyzcen[3 * iel + 0]);
                mei_tree_insert(ev_formula_turb, "y", xyzcen[3 * iel + 1]);
                mei_tree_insert(ev_formula_turb, "z", xyzcen[3 * iel + 2]);
                mei_evaluate(ev_formula_turb);
                rtp[vars->rtp[4]  * (*ncelet) + iel] = mei_tree_lookup(ev_formula_turb, "r11");
                rtp[vars->rtp[5]  * (*ncelet) + iel] = mei_tree_lookup(ev_formula_turb, "r22");
                rtp[vars->rtp[6]  * (*ncelet) + iel] = mei_tree_lookup(ev_formula_turb, "r33");
                rtp[vars->rtp[7]  * (*ncelet) + iel] = mei_tree_lookup(ev_formula_turb, "r12");
                rtp[vars->rtp[8]  * (*ncelet) + iel] = mei_tree_lookup(ev_formula_turb, "r13");
                rtp[vars->rtp[9]  * (*ncelet) + iel] = mei_tree_lookup(ev_formula_turb, "r23");
                rtp[vars->rtp[10] * (*ncelet) + iel] = mei_tree_lookup(ev_formula_turb, "eps");
              }
            }

            else if (cs_gui_strcmp(model, "Rij-EBRSM")) {
              const char *symbols[] = {"r11", "r22", "r133", "r12", "r13", "r23", "eps", "alpha"};
              if (mei_tree_find_symbols(ev_formula_turb, 8, symbols))
                bft_error(__FILE__, __LINE__, 0,
                          _("Error: can not find the required symbol: %s\n"),
                          "r11, r22, r33, r12, r13, r23, eps or alpha");

              for (icel = 0; icel < cells; icel++) {
                iel = cells_list[icel]-1;
                mei_tree_insert(ev_formula_turb, "x", xyzcen[3 * iel + 0]);
                mei_tree_insert(ev_formula_turb, "y", xyzcen[3 * iel + 1]);
                mei_tree_insert(ev_formula_turb, "z", xyzcen[3 * iel + 2]);
                mei_evaluate(ev_formula_turb);
                rtp[vars->rtp[4]  * (*ncelet) + iel] = mei_tree_lookup(ev_formula_turb, "r11");
                rtp[vars->rtp[5]  * (*ncelet) + iel] = mei_tree_lookup(ev_formula_turb, "r22");
                rtp[vars->rtp[6]  * (*ncelet) + iel] = mei_tree_lookup(ev_formula_turb, "r33");
                rtp[vars->rtp[7]  * (*ncelet) + iel] = mei_tree_lookup(ev_formula_turb, "r12");
                rtp[vars->rtp[8]  * (*ncelet) + iel] = mei_tree_lookup(ev_formula_turb, "r13");
                rtp[vars->rtp[9]  * (*ncelet) + iel] = mei_tree_lookup(ev_formula_turb, "r23");
                rtp[vars->rtp[10] * (*ncelet) + iel] = mei_tree_lookup(ev_formula_turb, "eps");
                rtp[vars->rtp[11] * (*ncelet) + iel] = mei_tree_lookup(ev_formula_turb, "alpha");
              }
            }

            else if (cs_gui_strcmp(model, "v2f-BL-v2/k")) {
              const char *symbols[] = {"k", "eps", "phi", "alpha"};
              if (mei_tree_find_symbols(ev_formula_turb, 4, symbols))
                bft_error(__FILE__, __LINE__, 0, _("Error: can not find the required symbol: %s\n"),
                          "k, eps, phi of al");

              for (icel = 0; icel < cells; icel++) {
                iel = cells_list[icel]-1;
                mei_tree_insert(ev_formula_turb, "x", xyzcen[3 * iel + 0]);
                mei_tree_insert(ev_formula_turb, "y", xyzcen[3 * iel + 1]);
                mei_tree_insert(ev_formula_turb, "z", xyzcen[3 * iel + 2]);
                mei_evaluate(ev_formula_turb);
                rtp[vars->rtp[4] * (*ncelet) + iel] = mei_tree_lookup(ev_formula_turb, "k");
                rtp[vars->rtp[5] * (*ncelet) + iel] = mei_tree_lookup(ev_formula_turb, "eps");
                rtp[vars->rtp[6] * (*ncelet) + iel] = mei_tree_lookup(ev_formula_turb, "phi");
                rtp[vars->rtp[7] * (*ncelet) + iel] = mei_tree_lookup(ev_formula_turb, "alpha");
              }
            }

            else if (cs_gui_strcmp(model, "k-omega-SST")) {
              const char *symbols[] = {"k", "omega"};
              if (mei_tree_find_symbols(ev_formula_turb, 2, symbols))
                bft_error(__FILE__, __LINE__, 0,
                          _("Error: can not find the required symbol: %s\n"),
                          "k or omega");

              for (icel = 0; icel < cells; icel++) {
                iel = cells_list[icel]-1;
                mei_tree_insert(ev_formula_turb, "x", xyzcen[3 * iel + 0]);
                mei_tree_insert(ev_formula_turb, "y", xyzcen[3 * iel + 1]);
                mei_tree_insert(ev_formula_turb, "z", xyzcen[3 * iel + 2]);
                mei_evaluate(ev_formula_turb);
                rtp[vars->rtp[4] * (*ncelet) + iel] = mei_tree_lookup(ev_formula_turb, "k");
                rtp[vars->rtp[5] * (*ncelet) + iel] = mei_tree_lookup(ev_formula_turb, "omega");
              }
            }

            else if (cs_gui_strcmp(model, "Spalart-Allmaras")) {
              const char *symbols[] = {"nusa"};
              if (mei_tree_find_symbols(ev_formula_turb, 1, symbols))
                bft_error(__FILE__, __LINE__, 0,
                          _("Error: can not find the required symbol: %s\n"),
                          "nusa");

              for (icel = 0; icel < cells; icel++) {
                iel = cells_list[icel]-1;
                mei_tree_insert(ev_formula_turb, "x", xyzcen[3 * iel + 0]);
                mei_tree_insert(ev_formula_turb, "y", xyzcen[3 * iel + 1]);
                mei_tree_insert(ev_formula_turb, "z", xyzcen[3 * iel + 2]);
                mei_evaluate(ev_formula_turb);
                rtp[vars->rtp[4] * (*ncelet) + iel] = mei_tree_lookup(ev_formula_turb, "nusa");
              }
            }

            else
              bft_error(__FILE__, __LINE__, 0,
                        _("Invalid turbulence model: %s.\n"), model);
            mei_tree_destroy(ev_formula_turb);
            BFT_FREE(model);
          }

          BFT_FREE(formula_turb);
          BFT_FREE(choice);
        }
      }

      /* User Scalars initialization */
      for (j=0; j < vars->nscaus; j++) {
        path_sca = cs_xpath_init_path();
        cs_xpath_add_elements(&path_sca, 2,
                              "additional_scalars",
                              "scalar");
        cs_xpath_add_test_attribute(&path_sca, "label", vars->label[j]);
        cs_xpath_add_element(&path_sca, "formula");
        cs_xpath_add_test_attribute(&path_sca, "zone_id", zone_id);
        cs_xpath_add_function_text(&path_sca);
        formula_sca = cs_gui_get_text_value(path_sca);
        BFT_FREE(path_sca);

        if (formula_sca != NULL) {
          ev_formula_sca = mei_tree_new(formula_sca);
          mei_tree_insert(ev_formula_sca,"x",0.);
          mei_tree_insert(ev_formula_sca,"y",0.);
          mei_tree_insert(ev_formula_sca,"z",0.);
          /* try to build the interpreter */
          if (mei_tree_builder(ev_formula_sca))
            bft_error(__FILE__, __LINE__, 0,
                      _("Error: can not interpret expression: %s\n %i"),
                      ev_formula_sca->string, mei_tree_builder(ev_formula_sca));

          if (mei_tree_find_symbol(ev_formula_sca, vars->label[j]))
            bft_error(__FILE__, __LINE__, 0,
                      _("Error: can not find the required symbol: %s\n"),
                      vars->label[j]);

          if (*isuite == 0 || (*isuite !=0 && iscold[j] == 0)) {
            for (icel = 0; icel < cells; icel++) {
              iel = cells_list[icel]-1;
              mei_tree_insert(ev_formula_sca, "x", xyzcen[3 * iel + 0]);
              mei_tree_insert(ev_formula_sca, "y", xyzcen[3 * iel + 1]);
              mei_tree_insert(ev_formula_sca, "z", xyzcen[3 * iel + 2]);
              mei_evaluate(ev_formula_sca);
              rtp[(isca[j]-1)*(*ncelet) + iel] =
                mei_tree_lookup(ev_formula_sca,vars->label[j]);
            }
          }
          mei_tree_destroy(ev_formula_sca);
        } else {
          if (*isuite == 0 || (*isuite !=0 && iscold[j] == 0)) {
            for (icel = 0; icel < cells; icel++) {
              iel = cells_list[icel]-1;
              rtp[(isca[j]-1)*(*ncelet) + iel] = 0.0;
            }
          }
        }
        BFT_FREE(formula_sca);
      }
      /* Meteo Scalars initialization */
      if (cs_gui_strcmp(vars->model, "atmospheric_flows")) {
        for (j=vars->nscaus;j< vars->nscaus + vars->nscapp;j++) {
          path_meteo = cs_xpath_init_path();
          cs_xpath_add_elements(&path_meteo, 3,
                                 "thermophysical_models",
                                 "atmospheric_flows",
                                 "scalar");
          cs_xpath_add_test_attribute(&path_meteo, "label", vars->label[j]);
          cs_xpath_add_element(&path_meteo, "formula");
          cs_xpath_add_test_attribute(&path_meteo, "zone_id", zone_id);
          cs_xpath_add_function_text(&path_meteo);
          formula_meteo = cs_gui_get_text_value(path_meteo);
          BFT_FREE(path_meteo);
          if (formula_meteo != NULL) {
            ev_formula_meteo = mei_tree_new(formula_meteo);
            mei_tree_insert(ev_formula_meteo,"x",0.);
            mei_tree_insert(ev_formula_meteo,"y",0.);
            mei_tree_insert(ev_formula_meteo,"z",0.);
            /* try to build the interpreter */
            if (mei_tree_builder(ev_formula_meteo))
              bft_error(__FILE__, __LINE__, 0,
                        _("Error: can not interpret expression: %s\n %i"),
                        ev_formula_meteo->string, mei_tree_builder(ev_formula_meteo));

            if (mei_tree_find_symbol(ev_formula_meteo, vars->label[j]))
              bft_error(__FILE__, __LINE__, 0,
                        _("Error: can not find the required symbol: %s\n"),
                        vars->label[j]);

            if (*isuite == 0 || (*isuite !=0 && iscold[j] == 0)) {
              for (icel = 0; icel < cells; icel++) {
                iel = cells_list[icel]-1;
                mei_tree_insert(ev_formula_meteo, "x", xyzcen[3 * iel + 0]);
                mei_tree_insert(ev_formula_meteo, "y", xyzcen[3 * iel + 1]);
                mei_tree_insert(ev_formula_meteo, "z", xyzcen[3 * iel + 2]);
                mei_evaluate(ev_formula_meteo);
                rtp[(isca[j]-1)*(*ncelet) + iel] =
                         mei_tree_lookup(ev_formula_meteo,vars->label[j]);
              }
            }
            mei_tree_destroy(ev_formula_meteo);
          }
          else{

            if (*isuite == 0 || (*isuite !=0 && iscold[j] == 0)) {
              for (icel = 0; icel < cells; icel++) {
                iel = cells_list[icel]-1;
                rtp[(isca[j]-1)*(*ncelet) + iel] = 0.0;
              }
            }
          }
          BFT_FREE(formula_meteo);
        }
      }

      if (cs_gui_strcmp(vars->model, "compressible_model")) {
        ccfth = 10000;
        //              pressure initialisation
        path = cs_xpath_init_path();
        cs_xpath_add_elements(&path, 3, "thermophysical_models","velocity_pressure","variable");
        cs_xpath_add_test_attribute(&path, "name", "pressure");
        cs_xpath_add_element(&path, "formula");
        cs_xpath_add_test_attribute(&path, "zone_id",zone_id);
        BFT_MALLOC(path1, strlen(path) +1, char);
        strcpy(path1, path);
        cs_xpath_add_attribute(&path, "status");
        buff = cs_gui_get_attribute_value(path);

        if (cs_gui_strcmp(buff,"on")) {
          ccfth = ccfth * 2;
          cs_xpath_add_function_text(&path1);
          formula = cs_gui_get_text_value(path1);
          ev_formula = _init_mei_tree(formula, "P");
          if (*isuite == 0 || (*isuite !=0 && iscold[j] == 0)) {
            for (icel = 0; icel < cells; icel++) {
              iel = cells_list[icel]-1;
              mei_tree_insert(ev_formula, "x", xyzcen[3 * iel + 0]);
              mei_tree_insert(ev_formula, "y", xyzcen[3 * iel + 1]);
              mei_tree_insert(ev_formula, "z", xyzcen[3 * iel + 2]);
              mei_evaluate(ev_formula);
              rtp[(*ipr-1) * (*ncelet) + iel] = mei_tree_lookup(ev_formula, "P");
            }
          }
          mei_tree_destroy(ev_formula);
        }
        BFT_FREE(buff);
        BFT_FREE(formula);
        BFT_FREE(path);
        BFT_FREE(path1);

        //              density initialisation
        path = cs_xpath_init_path();
        cs_xpath_add_elements(&path, 3, "thermophysical_models","compressible_model","property");
        cs_xpath_add_test_attribute(&path, "name", "Rho");
        cs_xpath_add_element(&path, "formula");
        cs_xpath_add_test_attribute(&path, "zone_id",zone_id);
        BFT_MALLOC(path1, strlen(path) +1, char);
        strcpy(path1, path);
        cs_xpath_add_attribute(&path, "status");
        buff = cs_gui_get_attribute_value(path);
        if (cs_gui_strcmp(buff,"on")) {
          ccfth = ccfth * 3;
          cs_xpath_add_function_text(&path1);
          formula = cs_gui_get_text_value(path1);
          ev_formula = _init_mei_tree(formula,"rho");
          if (*isuite == 0 || (*isuite !=0 && iscold[j] == 0)) {
            for (icel = 0; icel < cells; icel++) {
              iel = cells_list[icel]-1;
              mei_tree_insert(ev_formula, "x", xyzcen[3 * iel + 0]);
              mei_tree_insert(ev_formula, "y", xyzcen[3 * iel + 1]);
              mei_tree_insert(ev_formula, "z", xyzcen[3 * iel + 2]);
              mei_evaluate(ev_formula);
              propce[(*ipcrom-1)*(*ncelet) + iel] = mei_tree_lookup(ev_formula, "rho");
            }
          }
          mei_tree_destroy(ev_formula);
        }
        BFT_FREE(buff);
        BFT_FREE(formula);
        BFT_FREE(path);
        BFT_FREE(path1);

        //              Temperature initialisation
        path = cs_xpath_init_path();
        cs_xpath_add_elements(&path, 3, "thermophysical_models","compressible_model","scalar");
        cs_xpath_add_test_attribute(&path, "name", "TempK");
        cs_xpath_add_element(&path, "formula");
        cs_xpath_add_test_attribute(&path, "zone_id",zone_id);
        BFT_MALLOC(path1, strlen(path) +1, char);
        strcpy(path1, path);
        cs_xpath_add_attribute(&path, "status");
        buff = cs_gui_get_attribute_value(path);
        if (cs_gui_strcmp(buff,"on")){
          ccfth = ccfth * 5;
          cs_xpath_add_function_text(&path1);
          formula = cs_gui_get_text_value(path1);
          ev_formula = _init_mei_tree(formula,"T");
          if (*isuite == 0 || (*isuite !=0 && iscold[j] == 0)) {
            for (icel = 0; icel < cells; icel++) {
              iel = cells_list[icel]-1;
              mei_tree_insert(ev_formula, "x", xyzcen[3 * iel + 0]);
              mei_tree_insert(ev_formula, "y", xyzcen[3 * iel + 1]);
              mei_tree_insert(ev_formula, "z", xyzcen[3 * iel + 2]);
              mei_evaluate(ev_formula);
              rtp[(isca[*itempk-1]-1)*(*ncelet) + iel] = mei_tree_lookup(ev_formula, "T");
            }
          }
          mei_tree_destroy(ev_formula);
        }
        BFT_FREE(buff);
        BFT_FREE(formula);
        BFT_FREE(path);
        BFT_FREE(path1);
        //              Total energy initialisation
        path = cs_xpath_init_path();
        cs_xpath_add_elements(&path, 3, "thermophysical_models","compressible_model","scalar");
        cs_xpath_add_test_attribute(&path, "name", "EnergieT");
        cs_xpath_add_element(&path, "formula");
        cs_xpath_add_test_attribute(&path, "zone_id",zone_id);
        BFT_MALLOC(path1, strlen(path) +1, char);
        strcpy(path1, path);
        cs_xpath_add_attribute(&path, "status");
        buff = cs_gui_get_attribute_value(path);
        if (cs_gui_strcmp(buff,"on")) {
          ccfth = ccfth * 7;
          cs_xpath_add_function_text(&path1);
          formula = cs_gui_get_text_value(path1);
          ev_formula = _init_mei_tree(formula,"E");
          if (*isuite == 0 || (*isuite !=0 && iscold[j] == 0)) {
            for (icel = 0; icel < cells; icel++) {
              iel = cells_list[icel]-1;
              mei_tree_insert(ev_formula, "x", xyzcen[3 * iel + 0]);
              mei_tree_insert(ev_formula, "y", xyzcen[3 * iel + 1]);
              mei_tree_insert(ev_formula, "z", xyzcen[3 * iel + 2]);
              mei_evaluate(ev_formula);
              rtp[(isca[*ienerg-1]-1)*(*ncelet) + iel] = mei_tree_lookup(ev_formula, "E");
            }
          }
          mei_tree_destroy(ev_formula);
        }
        BFT_FREE(buff);
        BFT_FREE(formula);
        BFT_FREE(path);
        BFT_FREE(path1);
        *iccfth = ccfth;
    }

    BFT_FREE(cells_list);

#if _XML_DEBUG_
      bft_printf("--zone zone_id: %s\n", zone_id);
      bft_printf("--zone's element number: %i\n", cells);

      if (*isuite == 0) {
        double initial_value;
        for (j=1; j < vars->nvar - vars->nscaus - vars->nscapp; j++) {
          cs_gui_variable_initial_value(vars->name[j], zone_id, &initial_value);
          bft_printf("--initial value for %s: %f\n",
                     vars->name[j], initial_value);
        }
      }
#endif
      BFT_FREE(cells_list);
      BFT_FREE(zone_id);
    }
    BFT_FREE(status);
  } /* zones+1 */
}

/*-----------------------------------------------------------------------------
 * Change the head losses matrix from the local frame to the global frame.
 *
 * parameters:
 *   a_ij     -->  change matrix from the local frame to the global frame
 *   in_ij    -->  head losses matrix in the local frame
 *   out_ij   -->  head losses matrix in the global frame
 *----------------------------------------------------------------------------*/

static void
_matrix_base_conversion(double  a11,   double  a12,   double  a13,
                        double  a21,   double  a22,   double  a23,
                        double  a31,   double  a32,   double  a33,
                        double  in11,  double  in12,  double  in13,
                        double  in21,  double  in22,  double  in23,
                        double  in31,  double  in32,  double  in33,
                        double *out11, double *out12, double *out13,
                        double *out21, double *out22, double *out23,
                        double *out31, double *out32, double *out33)
{
  int     i, j, k;
  double  tensorP[3][3], tensorA[3][3], tensorB[3][3], tensorC[3][3], tensorD[3][3];

  tensorA[0][0] = in11;
  tensorA[0][1] = in12;
  tensorA[0][2] = in13;
  tensorA[1][0] = in21;
  tensorA[1][1] = in22;
  tensorA[1][2] = in23;
  tensorA[2][0] = in31;
  tensorA[2][1] = in32;
  tensorA[2][2] = in33;

  tensorP[0][0] = a11;
  tensorP[0][1] = a12;
  tensorP[0][2] = a13;
  tensorP[1][0] = a21;
  tensorP[1][1] = a22;
  tensorP[1][2] = a23;
  tensorP[2][0] = a31;
  tensorP[2][1] = a32;
  tensorP[2][2] = a33;

  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      tensorB[i][j] = 0.;
      for (k = 0; k < 3; k++)
        tensorB[i][j] += tensorP[i][k] * tensorA[k][j];
    }
  }

  /* Inversion of a 3x3 matrix */

  tensorC[0][0] = a11;
  tensorC[0][1] = a21;
  tensorC[0][2] = a31;
  tensorC[1][0] = a12;
  tensorC[1][1] = a22;
  tensorC[1][2] = a32;
  tensorC[2][0] = a13;
  tensorC[2][1] = a23;
  tensorC[2][2] = a33;

  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      tensorD[i][j] = 0.;
      for (k = 0; k < 3; k++)
        tensorD[i][j] += tensorB[i][k] * tensorC[k][j];
    }
  }

  *out11 = tensorD[0][0];
  *out22 = tensorD[1][1];
  *out33 = tensorD[2][2];
  *out12 = tensorD[0][1];
  *out13 = tensorD[0][2];
  *out21 = tensorD[1][0];
  *out23 = tensorD[1][2];
  *out31 = tensorD[2][0];
  *out32 = tensorD[2][1];
}

/*-----------------------------------------------------------------------------
 * Return value of coefficient associated to the head losses definition.
 *
 * parameters:
 *   zone_id   -->  id of the volume zone
 *   c         -->  name of the coefficient
 *----------------------------------------------------------------------------*/

static double
_c_heads_losses(const char* zone_id, const char* c)
{
  char* path;
  double result = 0.0;
  double value  = 0.0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "thermophysical_models", "heads_losses", "head_loss");
  cs_xpath_add_test_attribute(&path, "zone_id", zone_id);
  cs_xpath_add_element(&path, c);
  cs_xpath_add_function_text(&path);
  if (cs_gui_get_double(path, &result))
    value = result;
  else
    value= 0.0;
  BFT_FREE(path);
  return value;
}

/*----------------------------------------------------------------------------
 * Head losses definition
 *
 * Fortran Interface:
 *
 * subroutine uikpdc
 * *****************
 *
 * integer          iappel   -->  number of calls during a time step
 * integer          ncelet   -->  number of cells with halo
 * integer          ncepdp  <--   number of cells with head losses
 * integer          icepdc  <--   ncepdp cells number with head losses
 * double precision ckupdc  <--   head losses matrix
 * double precision rtpa     -->  variables array at previous time step
 *----------------------------------------------------------------------------*/

void CS_PROCF(uikpdc, UIKPDC)(const int*   iappel,
                              const int*   ncelet,
                                    int    ncepdp[],
                                    int    icepdc[],
                                    double ckupdc[],
                              const double rtpa[] )
{
  int i, j, iel, ielpdc, ikpdc;
  int zones = 0;
  int cells = 0;
  int *cells_list = NULL;
  double vit;
  double a11, a12, a13, a21, a22, a23, a31, a32, a33;
  double c11, c12, c13, c21, c22, c23, c31, c32, c33;
  double k11, k22, k33;
  char *zone_id = NULL;
  char *status = NULL;
  char *path = NULL;

  cs_var_t  *vars = cs_glob_var;

  /* number of volumic zone */

  zones
    = cs_gui_get_tag_number("/solution_domain/volumic_conditions/zone", 1);


  if (*iappel == 1 || *iappel == 2)
  {
    ielpdc = 0;

    for (i=1; i < zones+1; i++)
    {
      path = cs_xpath_init_path();
      cs_xpath_add_elements(&path, 2, "solution_domain", "volumic_conditions");
      cs_xpath_add_element_num(&path, "zone", i);
      cs_xpath_add_attribute(&path, "head_losses");
      status = cs_gui_get_attribute_value(path);
      BFT_FREE(path);

      if (cs_gui_strcmp(status, "on"))
      {
        zone_id = cs_gui_volumic_zone_id(i);
        cells_list = cs_gui_get_cells_list(zone_id, *ncelet, &cells);

        for (j=0; j < cells; j++)
        {
          if (*iappel == 2)
            icepdc[ielpdc] = cells_list[j];
          ielpdc++;
        }
        BFT_FREE(cells_list);
        BFT_FREE(zone_id);
      }
      BFT_FREE(status);
    } /* zones+1 */
    if (*iappel == 1)
      *ncepdp = ielpdc;
  }

  if (*iappel == 3)
  {
    for (ikpdc = 0; ikpdc < 6; ikpdc++)
      for (ielpdc = 0; ielpdc < *ncepdp; ielpdc++)
        ckupdc[ikpdc * (*ncepdp) + ielpdc] = 0.0;

    ielpdc = 0;

    for (i=1; i < zones+1; i++)
    {
      path = cs_xpath_init_path();
      cs_xpath_add_elements(&path, 2, "solution_domain", "volumic_conditions");
      cs_xpath_add_element_num(&path, "zone", i);
      cs_xpath_add_attribute(&path, "head_losses");
      status = cs_gui_get_attribute_value(path);
      BFT_FREE(path);

      if (cs_gui_strcmp(status, "on"))
      {
        zone_id = cs_gui_volumic_zone_id(i);
        cells_list = cs_gui_get_cells_list(zone_id, *ncelet, &cells);

        k11 = _c_heads_losses(zone_id, "kxx");
        k22 = _c_heads_losses(zone_id, "kyy");
        k33 = _c_heads_losses(zone_id, "kzz");

        a11 = _c_heads_losses(zone_id, "a11");
        a12 = _c_heads_losses(zone_id, "a12");
        a13 = _c_heads_losses(zone_id, "a13");
        a21 = _c_heads_losses(zone_id, "a21");
        a22 = _c_heads_losses(zone_id, "a22");
        a23 = _c_heads_losses(zone_id, "a23");
        a31 = _c_heads_losses(zone_id, "a31");
        a32 = _c_heads_losses(zone_id, "a32");
        a33 = _c_heads_losses(zone_id, "a33");

        if (   cs_gui_is_equal_real(a12, 0.0)
            && cs_gui_is_equal_real(a13, 0.0)
            && cs_gui_is_equal_real(a23, 0.0))
        {
          c11 = k11;
          c22 = k22;
          c33 = k33;
          c12 = 0.0;
          c13 = 0.0;
          c23 = 0.0;
        }
        else
        {
          _matrix_base_conversion(a11, a12, a13, a21, a22, a23, a31, a32, a33,
                                  k11, 0.0, 0.0, 0.0, k22, 0.0, 0.0, 0.0, k33,
                                 &c11, &c12, &c13, &c21, &c22, &c23, &c31, &c32, &c33);
        }

        for (j = 0; j < cells; j++)
        {
          iel = cells_list[j];
          vit = rtpa[vars->rtp[1]*(*ncelet) + iel-1] * rtpa[vars->rtp[1]*(*ncelet) + iel-1] \
              + rtpa[vars->rtp[2]*(*ncelet) + iel-1] * rtpa[vars->rtp[2]*(*ncelet) + iel-1] \
              + rtpa[vars->rtp[3]*(*ncelet) + iel-1] * rtpa[vars->rtp[3]*(*ncelet) + iel-1] ;
          vit = sqrt(vit);
          ckupdc[0 * (*ncepdp) + ielpdc] = 0.5 * c11 * vit;
          ckupdc[1 * (*ncepdp) + ielpdc] = 0.5 * c22 * vit;
          ckupdc[2 * (*ncepdp) + ielpdc] = 0.5 * c33 * vit;
          ckupdc[3 * (*ncepdp) + ielpdc] = 0.5 * c12 * vit;
          ckupdc[4 * (*ncepdp) + ielpdc] = 0.5 * c23 * vit;
          ckupdc[5 * (*ncepdp) + ielpdc] = 0.5 * c13 * vit;
          ielpdc++;
        }
        BFT_FREE(cells_list);
        BFT_FREE(zone_id);
      }
      BFT_FREE(status);
    } /* zones+1 */
  }
#if _XML_DEBUG_
  bft_printf("==>uikpdc\n");
  if (*iappel == 1)
    bft_printf("--%i number of head losses cells: %i\n", *iappel, *ncepdp);
  if (*iappel == 3)
    bft_printf("--%i number of head losses cells: %i\n", *iappel, ielpdc);
#endif
}

/*----------------------------------------------------------------------------
 * User law for material properties
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIPHYV
 * *****************
 *
 * INTEGER          NCEL     -->  number of cells whithout halo
 * INTEGER          NCELET   -->  number of cells whith halo
 * INTEGER          NSCAUS   -->  number of user scalar including thermal scalar
 * INTEGER          IROM     -->  pointer for density rho
 * INTEGER          IVISCL   -->  pointer for mulecular viscosity mu
 * INTEGER          ICP      -->  pointer for specific heat Cp
 * INTEGER          IVISLS   -->  pointer for Lambda/Cp
 * INTEGER          IROVAR   -->  =1 if rho variable, =0 if rho constant
 * INTEGER          IVIVAR   -->  =1 if mu variable, =0 if mu constant
 * INTEGER          ISCA     -->  indirection array for scalar number
 * INTEGER          ISCALT   -->  pointer for the thermal scalar in ISCA
 * INTEGER          ISCAVR   -->  scalars that are variance
 * INTEGER          IPPROC   -->  indirection array for cell properties
 * INTEGER          IVISCV   -->  pointer for volumic viscosity viscv
 * INTEGER          ITEMPK   -->  pointer for temperature (in K)
 * DOUBLE PRECISION P0       -->  pressure reference value
 * DOUBLE PRECISION T0       -->  temperature reference value
 * DOUBLE PRECISION RO0      -->  density reference value
 * DOUBLE PRECISION CP0      -->  specific heat reference value
 * DOUBLE PRECISION VISCL0   -->  dynamic viscosity reference value
 * DOUBLE PRECISION VISLS0   -->  diffusion coefficient of the scalars
 * DOUBLE PRECISION VISCV0   -->  volumic viscosity
 * DOUBLE PRECISION RTP      -->  variables and scalars array
 * DOUBLE PRECISION PROPCE   <--  cell properties array
 *----------------------------------------------------------------------------*/

void CS_PROCF(uiphyv, UIPHYV)(const cs_int_t  *const ncel,
                              const cs_int_t  *const ncelet,
                              const cs_int_t  *const nscaus,
                              const cs_int_t         irom[],
                              const cs_int_t         iviscl[],
                              const cs_int_t         icp[],
                              const cs_int_t         ivisls[],
                              const cs_int_t         irovar[],
                              const cs_int_t         ivivar[],
                              const cs_int_t         isca[],
                              const cs_int_t         iscalt[],
                              const cs_int_t         iscsth[],
                              const cs_int_t         iscavr[],
                              const cs_int_t         ipproc[],
                              const cs_int_t         iviscv[],
                              const cs_int_t         itempk[],
                              const cs_real_t        p0[],
                              const cs_real_t        t0[],
                              const cs_real_t        ro0[],
                              const cs_real_t        cp0[],
                              const cs_real_t        viscl0[],
                              const cs_real_t        visls0[],
                              const cs_real_t        viscv0[],
                              const cs_real_t        rtp[],
                                    cs_real_t        propce[])
{
  cs_var_t  *vars = cs_glob_var;
  mei_tree_t *ev_rho = NULL;
  mei_tree_t *ev_mu  = NULL;
  mei_tree_t *ev_cp  = NULL;
  mei_tree_t *ev_la  = NULL;
  mei_tree_t *ev_Ds  = NULL;
  mei_tree_t *ev_viscv  = NULL;
  char *law_rho = NULL;
  char *law_mu  = NULL;
  char *law_cp  = NULL;
  char *law_la  = NULL;
  char *law_Ds  = NULL;
  char *law_viscv  = NULL;

  char *path = NULL;
  int i, j, iel;
  double time0;

  int user_law = 0;
  int ipcrom = ipproc[ *irom   -1 ] -1;
  int ipcvis = ipproc[ *iviscl -1 ] -1;
  int ipccp  = ipproc[ *icp    -1 ] -1;
  int ipcvsv = ipproc[ *iviscv -1 ] -1;
  int ipcvsl = -1;  /* Lambda/Cp from the current thermal scalar
                       if the thermal scalar is Enthalpy or Energy
                       Lambda if the thermal scalar is Temperature */

  /* law for density */

  user_law = 0;
  if (*irovar == 1)
  {
    char *prop_choice = _properties_choice("density");
    if (cs_gui_strcmp(prop_choice, "user_law"))
      user_law = 1;
    BFT_FREE(prop_choice);
  }

  if (user_law)
  {
    /* search the formula for the law */

    path = cs_xpath_short_path();
    cs_xpath_add_element(&path, "property");
    cs_xpath_add_test_attribute(&path, "name", "density");
    cs_xpath_add_element(&path, "formula");
    cs_xpath_add_function_text(&path);

    law_rho = cs_gui_get_text_value(path);
    BFT_FREE(path);

    /* return an empty interpreter */

    time0 = cs_timer_wtime();

    ev_rho = mei_tree_new(law_rho);

    mei_tree_insert(ev_rho, "rho0", *ro0);
    mei_tree_insert(ev_rho, "p0", *p0);

    for (i = 0; i < *nscaus; i++)
      mei_tree_insert(ev_rho, vars->label[i], 0.0);

    /* try to build the interpreter */

    if (mei_tree_builder(ev_rho))
      bft_error(__FILE__, __LINE__, 0,
                _("Error: can not interpret expression: %s\n"), ev_rho->string);

    if (mei_tree_find_symbol(ev_rho, "rho"))
      bft_error(__FILE__, __LINE__, 0,
                _("Error: can not find the required symbol: %s\n"), "rho");

    /* for each cell, update the value of the table of symbols for each scalar
       (including the thermal scalar), and evaluate the interpreter */

    for (iel = 0; iel < *ncel; iel++)
    {
      for (i = 0; i < *nscaus; i++)
        mei_tree_insert(ev_rho,
                        vars->label[i],
                        rtp[(isca[i] -1) * (*ncelet) + iel]);

      mei_evaluate(ev_rho);
      propce[ipcrom * (*ncelet) + iel] = mei_tree_lookup(ev_rho, "rho");
    }

    mei_tree_destroy(ev_rho);

    cs_gui_add_mei_time(cs_timer_wtime() - time0);
  }

  /* law for molecular viscosity */

  user_law = 0;
  if (*ivivar == 1)
  {
    char *prop_choice = _properties_choice("molecular_viscosity");
    if (cs_gui_strcmp(prop_choice, "user_law"))
      user_law = 1;
    BFT_FREE(prop_choice);
  }

  if (user_law)
  {
    /* search the formula for the law */

    path = cs_xpath_short_path();
    cs_xpath_add_element(&path, "property");
    cs_xpath_add_test_attribute(&path, "name", "molecular_viscosity");
    cs_xpath_add_element(&path, "formula");
    cs_xpath_add_function_text(&path);

    law_mu = cs_gui_get_text_value(path);
    BFT_FREE(path);

    /* return an empty interpreter */

    time0 = cs_timer_wtime();

    ev_mu = mei_tree_new(law_mu);

    mei_tree_insert(ev_mu, "rho0", *ro0);
    mei_tree_insert(ev_mu, "mu0", *viscl0);
    mei_tree_insert(ev_mu, "p0", *p0);
    mei_tree_insert(ev_mu, "rho", 0.0);
    if (cs_gui_strcmp(vars->model, "compressible_model"))
      mei_tree_insert(ev_mu, "t0", 0.0);

    for (i = 0; i < *nscaus; i++)
      mei_tree_insert(ev_mu, vars->label[i], 0.0);

    /* try to build the interpreter */

    if (mei_tree_builder(ev_mu))
      bft_error(__FILE__, __LINE__, 0,
                _("Error: can not interpret expression: %s\n"), ev_mu->string);

    if (mei_tree_find_symbol(ev_mu, "mu"))
      bft_error(__FILE__, __LINE__, 0,
                _("Error: can not find the required symbol: %s\n"), "mu");

    /* for each cell, update the value of the table of symbols for each scalar
       (including the thermal scalar) and for the density,
       then evaluate the interpreter */

    for (iel = 0; iel < *ncel; iel++)
    {
      for (i = 0; i < *nscaus; i++)
        mei_tree_insert(ev_mu,
                        vars->label[i],
                        rtp[(isca[i] -1) * (*ncelet) + iel]);

      mei_tree_insert(ev_mu,
                      "rho",
                      propce[ipcrom * (*ncelet) + iel]);

      if (cs_gui_strcmp(vars->model, "compressible_model"))
        mei_tree_insert(ev_mu, "T", rtp[(isca[*itempk -1] -1) * (*ncelet) + iel]);

      mei_evaluate(ev_mu);
      propce[ipcvis * (*ncelet) + iel] = mei_tree_lookup(ev_mu, "mu");
    }

    mei_tree_destroy(ev_mu);

    cs_gui_add_mei_time(cs_timer_wtime() - time0);
  }

  /* law for specific heat */

  user_law = 0;
  if (*icp > 0)
  {
    char *prop_choice = _properties_choice("specific_heat");
    if (cs_gui_strcmp(prop_choice, "user_law"))
      user_law = 1;
    BFT_FREE(prop_choice);
  }

  if (user_law)
  {
    /* search the formula for the law */

    path = cs_xpath_short_path();
    cs_xpath_add_element(&path, "property");
    cs_xpath_add_test_attribute(&path, "name", "specific_heat");
    cs_xpath_add_element(&path, "formula");
    cs_xpath_add_function_text(&path);

    law_cp = cs_gui_get_text_value(path);
    BFT_FREE(path);

    /* return an empty interpreter */

    time0 = cs_timer_wtime();

    ev_cp = mei_tree_new(law_cp);

    mei_tree_insert(ev_cp, "cp0", *cp0);
    mei_tree_insert(ev_cp, "p0", *p0);

    for (i = 0; i < *nscaus; i++)
      mei_tree_insert(ev_cp, vars->label[i], 0.0);

    /* try to build the interpreter */

    if (mei_tree_builder(ev_cp))
      bft_error(__FILE__, __LINE__, 0,
                _("Error: can not interpret expression: %s\n"), ev_cp->string);

    if (mei_tree_find_symbol(ev_cp, "cp"))
      bft_error(__FILE__, __LINE__, 0,
                _("Error: can not find the required symbol: %s\n"), "cp");

    /* for each cell, update the value of the table of symbols for each scalar
       (including the thermal scalar), and evaluate the interpreter */

    for (iel = 0; iel < *ncel; iel++) {
      for (i = 0; i < *nscaus; i++)
        mei_tree_insert(ev_cp,
                        vars->label[i],
                        rtp[(isca[i] -1) * (*ncelet) + iel]);

      mei_evaluate(ev_cp);
      propce[ipccp * (*ncelet) + iel] = mei_tree_lookup(ev_cp, "cp");
    }

    mei_tree_destroy(ev_cp);

    cs_gui_add_mei_time(cs_timer_wtime() - time0);
  }

  /* law for thermal conductivity */

  user_law = 0;
  if (*iscalt > 0)
  {
    if (ivisls[*iscalt -1] > 0)
    {
      char *prop_choice = _properties_choice("thermal_conductivity");
      if (cs_gui_strcmp(prop_choice, "user_law"))
        user_law = 1;
      BFT_FREE(prop_choice);
    }
  }

  if (user_law)
  {
    ipcvsl = ipproc[ ivisls[*iscalt -1 ] -1 ] -1;

    /* search the formula for the law */

    path = cs_xpath_short_path();
    cs_xpath_add_element(&path, "property");
    cs_xpath_add_test_attribute(&path, "name", "thermal_conductivity");
    cs_xpath_add_element(&path, "formula");
    cs_xpath_add_function_text(&path);

    law_la = cs_gui_get_text_value(path);
    BFT_FREE(path);

    /* return an empty interpreter */

    time0 = cs_timer_wtime();

    ev_la = mei_tree_new(law_la);

    /* for the Temperature, the diffusivity factor is not divided by Cp */
    if (abs(iscsth[*iscalt-1]) != 1)
    {
      mei_tree_insert(ev_la, "lambda0", visls0[*iscalt-1]*(*cp0));
    }
    else
    {
      mei_tree_insert(ev_la, "lambda0", visls0[*iscalt-1]);
    }
    mei_tree_insert(ev_la, "p0", *p0);

    for (i = 0; i < *nscaus; i++)
      mei_tree_insert(ev_la, vars->label[i], 0.0);

    /* try to build the interpreter */

    if (mei_tree_builder(ev_la))
      bft_error(__FILE__, __LINE__, 0,
                _("Error: can not interpret expression: %s\n"), ev_la->string);

    if (mei_tree_find_symbol(ev_la, "lambda"))
      bft_error(__FILE__, __LINE__, 0,
                _("Error: can not find the required symbol: %s\n"), "lambda");

    /* for each cell, update the value of the table of symbols for each scalar
       (including the thermal scalar), and evaluate the interpreter */

    if (*icp > 0)
    {
      for (iel = 0; iel < *ncel; iel++)
      {
        for (i = 0; i < *nscaus; i++)
          mei_tree_insert(ev_la,
                          vars->label[i],
                          rtp[(isca[i] -1) * (*ncelet) + iel]);

        mei_evaluate(ev_la);
        /* for the Temperature, the diffusivity factor is not divided by Cp */
        if (abs(iscsth[*iscalt - 1]) != 1)
        {
          propce[ipcvsl * (*ncelet) + iel] =
          mei_tree_lookup(ev_la, "lambda") / propce[ipccp * (*ncelet) + iel];
        }
        else
        {
          propce[ipcvsl * (*ncelet) + iel] =
          mei_tree_lookup(ev_la, "lambda");
        }
      }
    }
    else {
      for (iel = 0; iel < *ncel; iel++) {
        for (i = 0; i < *nscaus; i++)
          mei_tree_insert(ev_la,
                          vars->label[i],
                          rtp[(isca[i] -1) * (*ncelet) + iel]);

        mei_evaluate(ev_la);
        /* for the Temperature, the diffusivity factor is not divided by Cp */
        if (abs(iscsth[*iscalt - 1]) != 1) {
          propce[ipcvsl * (*ncelet) + iel] =
          mei_tree_lookup(ev_la, "lambda") / *cp0;
        }
        else {
          propce[ipcvsl * (*ncelet) + iel] =
          mei_tree_lookup(ev_la, "lambda");
        }
      }
    }
    mei_tree_destroy(ev_la);

    cs_gui_add_mei_time(cs_timer_wtime() - time0);
  }

  /* law for thermal conductivity (compressible model) */

  if (cs_gui_strcmp(vars->model, "compressible_model")) {
    if (ivisls[*itempk -1] > 0) {
      char *prop_choice = _properties_choice("thermal_conductivity");
      if (cs_gui_strcmp(prop_choice, "user_law"))
        user_law = 1;
      BFT_FREE(prop_choice);
    }

    if (user_law) {
      ipcvsl = ipproc[ ivisls[*itempk -1 ] -1 ] -1;

      /* search the formula for the law */

      path = cs_xpath_short_path();
      cs_xpath_add_element(&path, "property");
      cs_xpath_add_test_attribute(&path, "name", "thermal_conductivity");
      cs_xpath_add_element(&path, "formula");
      cs_xpath_add_function_text(&path);

      law_la = cs_gui_get_text_value(path);
      BFT_FREE(path);

      /* return an empty interpreter */

      time0 = cs_timer_wtime();

      ev_la = mei_tree_new(law_la);

      mei_tree_insert(ev_la, "lambda0", visls0[*itempk -1]);
      mei_tree_insert(ev_la, "p0", *p0);
      mei_tree_insert(ev_la, "t0", *t0);
      mei_tree_insert(ev_la, "rho0", *ro0);

      for (i = 0; i < *nscaus; i++)
        mei_tree_insert(ev_la, vars->label[i], 0.0);

      /* try to build the interpreter */

      if (mei_tree_builder(ev_la))
        bft_error(__FILE__, __LINE__, 0,
                  _("Error: can not interpret expression: %s\n"), ev_la->string);

      if (mei_tree_find_symbol(ev_la, "lambda"))
        bft_error(__FILE__, __LINE__, 0,
                  _("Error: can not find the required symbol: %s\n"), "lambda");

      /* for each cell, update the value of the table of symbols for each scalar
         (including the thermal scalar), and evaluate the interpreter */

      for (iel = 0; iel < *ncel; iel++) {
        for (i = 0; i < *nscaus; i++)
          mei_tree_insert(ev_la,
                          vars->label[i],
                          rtp[(isca[i] -1) * (*ncelet) + iel]);

        mei_tree_insert(ev_la, "T", rtp[(isca[*itempk -1] -1) * (*ncelet) + iel]);

        mei_evaluate(ev_la);
        propce[ipcvsl * (*ncelet) + iel] = mei_tree_lookup(ev_la, "lambda");
      }
      mei_tree_destroy(ev_la);

      cs_gui_add_mei_time(cs_timer_wtime() - time0);
    }
  }

  /* law for scalar diffusivity */

  for (j = 0; j < *nscaus; j++) {
    char *name = _scalar_diffusion_coefficient_name(j);

    user_law = 0;
    if (j != *iscalt -1 && iscavr[j] <= 0 && ivisls[j] > 0) {
      char *prop_choice = _properties_choice("name");
      if (cs_gui_strcmp(prop_choice, "user_law"))
        user_law = 1;
      BFT_FREE(prop_choice);
    }

    if (user_law) {
      ipcvsl = ipproc[ ivisls[j] -1 ] -1;

      /* search the formula for the law */

      path = cs_xpath_init_path();
      cs_xpath_add_element(&path, "additional_scalars");
      cs_xpath_add_element_num(&path, "scalar", j+1);
      cs_xpath_add_element(&path, "property");
      cs_xpath_add_element(&path, "formula");
      cs_xpath_add_function_text(&path);

      law_Ds = cs_gui_get_text_value(path);
      BFT_FREE(path);

      /* return an empty interpreter */

      time0 = cs_timer_wtime();

      ev_Ds = mei_tree_new(law_Ds);
      BFT_FREE(law_Ds);
      for (i = 0; i < *nscaus; i++)
        mei_tree_insert(ev_Ds, vars->label[i], 0.0);

      /* try to build the interpreter */

      if (mei_tree_builder(ev_Ds))
        bft_error(__FILE__, __LINE__, 0,
                  _("Error: can not interpret expression: %s\n"), ev_Ds->string);

      if (mei_tree_find_symbol(ev_Ds, "diffusivity"))
        bft_error(__FILE__, __LINE__, 0,
                  _("Error: can not find the required symbol: %s\n"),
                  "diffusivity");

      /* for each cell, update the value of the table of symbols for each scalar
         (including the thermal scalar), and evaluate the interpreter */

      if (*irovar == 1) {
        for (iel = 0; iel < *ncel; iel++) {
          for (i = 0; i < *nscaus; i++)
            mei_tree_insert(ev_Ds,
                            vars->label[i],
                            rtp[(isca[i] -1) * (*ncelet) + iel]);

          mei_evaluate(ev_Ds);
          propce[ipcvsl * (*ncelet) + iel]
            =    mei_tree_lookup(ev_Ds, "diffusivity")
               * propce[ipcrom * (*ncelet) + iel];
        }
      }
      else {
        for (iel = 0; iel < *ncel; iel++) {
          for (i = 0; i < *nscaus; i++)
            mei_tree_insert(ev_Ds,
                            vars->label[i],
                            rtp[(isca[i] -1) * (*ncelet) + iel]);

          mei_evaluate(ev_Ds);
          propce[ipcvsl * (*ncelet) + iel] =
            mei_tree_lookup(ev_Ds, "diffusivity") * (*ro0);
        }
      }
      mei_tree_destroy(ev_Ds);

      cs_gui_add_mei_time(cs_timer_wtime() - time0);

    }
    BFT_FREE(name);
  }

  /* law for volumic viscosity (compressible model) */
  if (cs_gui_strcmp(vars->model, "compressible_model")) {
    user_law = 0;
    if (*iviscv > 0) {
      char *prop_choice = _properties_choice("volumic_viscosity");
      if (cs_gui_strcmp(prop_choice, "user_law"))
        user_law = 1;
      BFT_FREE(prop_choice);
    }

    if (user_law) {
      /* search the formula for the law */

      path = cs_xpath_short_path();
      cs_xpath_add_element(&path, "property");
      cs_xpath_add_test_attribute(&path, "name", "volumic_viscosity");
      cs_xpath_add_element(&path, "formula");
      cs_xpath_add_function_text(&path);

      law_viscv = cs_gui_get_text_value(path);
      BFT_FREE(path);

      /* return an empty interpreter */

      ev_viscv = mei_tree_new(law_viscv);
      mei_tree_insert(ev_viscv,"viscv0", *viscv0);
      mei_tree_insert(ev_viscv,"p0", *p0);
      mei_tree_insert(ev_viscv,"T", 0.);
      mei_tree_insert(ev_viscv,"t0", *t0);

      /* try to build the interpreter */

      if (mei_tree_builder(ev_viscv))
        bft_error(__FILE__, __LINE__, 0,
            _("Error: can not interpret expression: %s\n"), ev_viscv->string);

      if (mei_tree_find_symbol(ev_viscv, "viscv"))
        bft_error(__FILE__, __LINE__, 0,
            _("Error: can not find the required symbol: %s\n"), "viscv");

      /* for each cell, update the value of the table of symbols for each scalar
         (including the thermal scalar), and evaluate the interpreter */

      for (iel = 0; iel < *ncel; iel++) {
        mei_tree_insert(ev_viscv, "T",
                        rtp[(isca[*itempk -1] -1) * (*ncelet) + iel]);

        mei_evaluate(ev_viscv);
        propce[ipcvsv * (*ncelet) + iel] = mei_tree_lookup(ev_viscv, "viscv");
      }

      mei_tree_destroy(ev_viscv);
    }
  }

#if _XML_DEBUG_
  bft_printf("==>UIPHYV\n");
  if (*irovar == 1)
    bft_printf("--law for density: %s\n", law_rho);

  if (*ivivar == 1)
    bft_printf("--law for viscosity: %s\n", law_mu);

  if (*icp > 0)
    bft_printf("--law for specific heat: %s\n", law_cp);

  if (*iscalt > 0) {
    if (ivisls[*iscalt -1] > 0)
      bft_printf("--law for thermal conductivity: %s\n", law_la);
  }

  for (j = 0; j < *nscaus; j++) {
    if (j != *iscalt -1 && iscavr[j] <= 0 && ivisls[j] > 0) {
      ipcvsl = ipproc[ ivisls[j] -1 ] -1;

      path = cs_xpath_init_path();
      cs_xpath_add_element(&path, "additional_scalars");
      cs_xpath_add_element_num(&path, "scalar", j+1);
      cs_xpath_add_element(&path, "property");
      cs_xpath_add_element(&path, "formula");
      cs_xpath_add_function_text(&path);

      law_Ds = cs_gui_get_text_value(path);
      bft_printf("--law for the coefficient of diffusity of the scalar %s: %s\n",
                 vars->label[j], law_Ds);
      BFT_FREE(path);
      BFT_FREE(law_Ds);
    }
  }
#endif

  BFT_FREE(law_rho);
  BFT_FREE(law_mu);
  BFT_FREE(law_cp);
  BFT_FREE(law_la);
}

/*----------------------------------------------------------------------------
 * 1D profile postprocessing
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIPROF
 * *****************
 *
 * INTEGER          NCELET   -->  number of cells with halo
 * INTEGER          NCEL     -->  number of cells without halo
 * INTEGER          NTMABS   -->  max iterations numbers
 * INTEGER          NTCABS   -->  current iteration number
 * DOUBLE PRECISION TTCABS   -->  current physical time
 * DOUBLE PRECISION TTMABS   -->  max physical time
 * DOUBLE PRECISION TTPABS   -->  physical time at calculation beginning
 * DOUBLE PRECISION XYZCEN   -->  cell's gravity center
 * DOUBLE PRECISION RTP      -->  variables and scalars array
 * DOUBLE PRECISION PROPCE   -->  property array
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiprof, UIPROF) (const int    *const ncelet,
                                const int    *const ncel,
                                const int    *const ntmabs,
                                const int    *const ntcabs,
                                const double *const ttcabs,
                                const double *const ttmabs,
                                const double *const ttpabs,
                                const double *const xyzcen,
                                const double *const rtp,
                                const double *const propce)
{
  FILE *file = NULL;
  char *filename = NULL;
  char *title = NULL;
  char *name = NULL;
  char *path = NULL;
  char *formula = NULL;
  char *output_type = NULL;

  int output_format = 0;
  int fic_nbr = 0;
  int i, ii, iii, j;
  int npoint, iel1, irang1, iel, irangv;
  int nvar_prop, nvar_prop4, output_frequency;
  double time_output;
  double x1 = 0., y1 = 0., z1 = 0.;
  double xx, yy, zz, xyz[3];
  double a, aa;
  double *array;
  static int ipass = 0;
  bool status;

  mei_tree_t *ev_formula  = NULL;

  cs_var_t  *vars = cs_glob_var;

  /* get the number of 1D profile file to write*/

  fic_nbr = cs_gui_get_tag_number("/analysis_control/profiles/profile", 1);

  if (!fic_nbr) return;

  for (i = 0 ; i < fic_nbr ; i++) {

    /* for each profile, check the output frequency */

    output_format = _get_profile_format(i);
    output_type = _get_profile_output_type(i);
    time_output = 0.;
    status = false;

    if (cs_gui_strcmp(output_type, "time_value")) {
      time_output = _get_profile_coordinate(i, "output_frequency");
      int ifreqs = (int)((*ttcabs - *ttpabs) / time_output);
      if ((ifreqs > ipass) || (*ttcabs >= *ttmabs && *ttmabs > 0.))
        status = true;
    } else {
      output_frequency = (int) _get_profile_coordinate(i, "output_frequency");
      if ((*ntmabs == *ntcabs) || (output_frequency > 0 && (*ntcabs % output_frequency) == 0))
        status = true;
    }
    BFT_FREE(output_type);

    if (status) {

      ipass++;
      path = cs_xpath_init_path();
      cs_xpath_add_elements(&path, 2, "analysis_control", "profiles");
      cs_xpath_add_element_num(&path, "profile", i+1);
      cs_xpath_add_element(&path, "formula");
      cs_xpath_add_function_text(&path);
      formula = cs_gui_get_text_value(path);
      ev_formula = mei_tree_new(formula);
      mei_tree_insert(ev_formula, "s", 0.0);
      BFT_FREE(formula);
      BFT_FREE(path);

      /* try to build the interpreter */

      if (mei_tree_builder(ev_formula))
        bft_error(__FILE__, __LINE__, 0,
                  _("Error: can not interpret expression: %s\n %i"),
                  ev_formula->string, mei_tree_builder(ev_formula));

      const char *coord[] = {"x","y","z"};

      if (mei_tree_find_symbols(ev_formula, 3, coord))
        bft_error(__FILE__, __LINE__, 0,
                  _("Error: can not find the required symbol: %s\n"),
                  "x, y or z");

      nvar_prop = _get_profile_names_number(i);
      nvar_prop4 = nvar_prop + 4;
      BFT_MALLOC(array, nvar_prop4, double);

      /* Only the first processor rank opens the file */

      if (cs_glob_rank_id <= 0) {

        filename = _get_profile("label", i);
        title    = _get_profile("title", i);

        if (output_frequency > 0 || time_output > 0.) {

          char buf1[5];
          const char buffer[5] = "%.4i";

          /* Extension creation : format stored in 'buffer' */

          sprintf(buf1, buffer, *ntcabs);

          BFT_REALLOC(filename, strlen(filename) + 1 + 4 + 1, char);

          strcat(filename, "_");
          strcat(filename, buf1);

        }

        BFT_REALLOC(filename, strlen(filename) + 4 + 1, char);
        if (output_format == 0)
          strcat(filename, ".dat");
        else
          strcat(filename, ".csv");
        file = fopen(filename, "w");

        if (file ==  NULL) {
          cs_base_warn(__FILE__, __LINE__);
          bft_printf( _("Unable to open the file: %s\n"), filename);
          break;
        }

        if (output_format == 0) {
          fprintf(file, "# Code_Saturne results 1D profile\n#\n");
          fprintf(file, "# Iteration output: %i\n", *ntcabs);
          fprintf(file, "# Time output:     %12.5e\n#\n", *ttcabs);
          fprintf(file, "#TITLE: %s\n", title);
          fprintf(file, "#COLUMN_TITLES: Distance | X | Y | Z");
          for (ii = 0 ; ii < nvar_prop ; ii++) {
            char *buffer = _get_profile_label_name(i, ii);
            fprintf(file, " | %s", buffer);
            BFT_FREE(buffer);
          }
          fprintf(file, "\n");
        }
        else {
          fprintf(file, "s, x, y, z");
          for (ii = 0 ; ii < nvar_prop ; ii++) {
            char *buffer = _get_profile_label_name(i, ii);
            fprintf(file, ", %s", buffer);
            BFT_FREE(buffer);
          }
          fprintf(file, "\n");
        }
        BFT_FREE(filename);
        BFT_FREE(title);
      }

      path = cs_xpath_init_path();
      cs_xpath_add_elements(&path, 2, "analysis_control", "profiles");
      cs_xpath_add_element_num(&path, "profile", i+1);
      cs_xpath_add_element(&path, "points");
      cs_xpath_add_function_text(&path);
      if (!cs_gui_get_int(path, &npoint))
        bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);
      BFT_FREE(path);

      iel1   = -999;
      irang1 = -999;

      a = 1. / (double) (npoint-1);

      for (ii = 0; ii < npoint; ii++) {

        aa = ii*a;
        mei_tree_insert(ev_formula,"s",aa);
        mei_evaluate(ev_formula);

        xyz[0] = mei_tree_lookup(ev_formula,"x");
        xyz[1] = mei_tree_lookup(ev_formula,"y");
        xyz[2] = mei_tree_lookup(ev_formula,"z");

        if (ii == 0) {
          x1 = xyz[0];
          y1 = xyz[1];
          z1 = xyz[2];
        }

        CS_PROCF(findpt, FINDPT)(ncelet,  ncel,    xyzcen,
                                &xyz[0], &xyz[1], &xyz[2],
                                &iel,    &irangv);

        if ((iel != iel1) || (irangv != irang1)) {
          iel1 = iel;
          irang1 = irangv;

          if (cs_glob_rank_id == irangv) {

            iel--;
            xx = xyzcen[3 * iel + 0];
            yy = xyzcen[3 * iel + 1];
            zz = xyzcen[3 * iel + 2];
            array[1] = xx;
            array[2] = yy;
            array[3] = zz;
            xx = xx - x1;
            yy = yy - y1;
            zz = zz - z1;
            array[0] = sqrt(xx*xx + yy*yy + zz*zz);

            for (iii=0; iii < nvar_prop; iii++) {

              name = _get_profile_name(i, iii);

              for (j=0; j < vars->nvar; j++) {
                if (cs_gui_strcmp(name,  vars->name[j]))
                  array[iii+4] = rtp[vars->rtp[j] * (*ncelet) + iel];
              }

              for (j=0; j < vars->nprop; j++) {
                if (cs_gui_strcmp(name, vars->properties_name[j]))
                  array[iii+4]
                    = propce[vars->propce[j] * (*ncelet) + iel];
              }

              BFT_FREE(name);
            }

          } else {

            for (iii=0; iii < nvar_prop4; iii++)
              array[iii] = 0.0;
          }

          /* Send to other processors if parallel */
#if defined(HAVE_MPI)
          if (cs_glob_rank_id >= 0) {
            MPI_Bcast(array,
                      nvar_prop4,
                      CS_MPI_REAL,
                      irangv,
                      cs_glob_mpi_comm);
          }
#endif

          if (cs_glob_rank_id <= 0) {
            if (output_format == 0) {
              for (iii=0; iii < nvar_prop4; iii++)
                fprintf(file, "%12.5e ", array[iii]);
              fprintf(file, "\n");
            }
            else {
              if (nvar_prop > 0) {
                for (iii=0; iii < nvar_prop4 - 1; iii++)
                  fprintf(file, "%12.5e, ", array[iii]);
                fprintf(file, "%12.5e ", array[nvar_prop4 - 1]);
              }
              fprintf(file, "\n");
            }
          }
        }
      }
      mei_tree_destroy(ev_formula);

      if (cs_glob_rank_id <= 0) fclose(file);

      BFT_FREE(array);
    }
  }
}

/*----------------------------------------------------------------------------
 * Free memory: clean global private variables and libxml2 variables
 *
 * Fortran Interface:
 *
 * SUBROUTINE MEMUI1
 * *****************
 *
 * INTEGER          NCHARB  --> number of coal
 *----------------------------------------------------------------------------*/

void CS_PROCF (memui1, MEMUI1) (const int *const ncharb)
{
  cs_gui_boundary_conditions_free_memory(ncharb);
  cs_gui_clean_memory();
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Get initial value from property markup.
 *
 * parameters:
 *   property_name      <--   name of the property
 *   value              -->   new initial value of the property
 *----------------------------------------------------------------------------*/

void
cs_gui_properties_value(const char  *property_name,
                        double      *value)
{
  char   *path = NULL;
  double  result;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "property");
  cs_xpath_add_test_attribute(&path, "name", property_name);
  cs_xpath_add_element(&path, "initial_value");
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &result))
    *value = result;

  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Initialization choice of the reference variables parameters.
 *
 * parameters:
 *   name            <--   parameter name
 *   value           -->   parameter value
 *----------------------------------------------------------------------------*/

void
cs_gui_reference_initialization(const char  *param,
                                double      *value)
{
  char *path = NULL;
  double  result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3,
                        "thermophysical_models",
                        "reference_values",
                        param);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &result))
    *value = result;
  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Set partitioning options.
 *----------------------------------------------------------------------------*/

void
cs_gui_partition(void)
{
  char  *path = NULL;
  char  *part_name = NULL;
  char  *s_perio = NULL;
  char  *s_output = NULL;
  char  *s_list = NULL;

  cs_partition_algorithm_t a = CS_PARTITION_DEFAULT;
  bool ignore_perio = false;
  int  rank_step = 1;
  int  write_level = 1;
  int  n_add_parts = 0;
  int  *add_parts = NULL;

  if (!cs_gui_file_is_loaded())
    return;

  /* Partitioning type */

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3,
                        "calculation_management", "partitioning", "type");
  cs_xpath_add_function_text(&path);

  part_name = cs_gui_get_text_value(path);

  if (part_name != NULL) {
    if (!strcmp(part_name, "default"))
      a = CS_PARTITION_DEFAULT;
    else if (!strcmp(part_name, "morton sfc"))
      a = CS_PARTITION_SFC_MORTON_BOX;
    else if (!strcmp(part_name, "morton sfc cube"))
      a = CS_PARTITION_SFC_MORTON_CUBE;
    else if (!strcmp(part_name, "hilbert sfc"))
      a = CS_PARTITION_SFC_HILBERT_BOX;
    else if (!strcmp(part_name, "hilbert sfc cube"))
      a = CS_PARTITION_SFC_HILBERT_CUBE;
    else if (!strcmp(part_name, "scotch"))
      a = CS_PARTITION_SCOTCH;
    else if (!strcmp(part_name, "metis"))
      a = CS_PARTITION_METIS;
    else if (!strcmp(part_name, "block"))
      a = CS_PARTITION_BLOCK;
    BFT_FREE(part_name);
  }

  BFT_FREE(path);

  /* Rank step */

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3,
                        "calculation_management",
                        "partitioning",
                        "rank_step");
  cs_xpath_add_function_text(&path);
  cs_gui_get_int(path, &rank_step);

  BFT_FREE(path);

  /* Ignore periodicity option */

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3,
                        "calculation_management",
                        "partitioning",
                        "ignore_periodicity");
  cs_xpath_add_attribute(&path, "status");
  s_perio = cs_gui_get_attribute_value(path);
  if (s_perio != NULL) {
    if (cs_gui_strcmp(s_perio, "on"))
      ignore_perio = true;
    BFT_FREE(s_perio);
  }

  BFT_FREE(path);

  /* Output option */

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3,
                        "calculation_management", "partitioning", "output");
  cs_xpath_add_function_text(&path);

  s_output = cs_gui_get_text_value(path);

  if (s_output != NULL) {
    if (!strcmp(s_output, "no"))
      write_level = 0;
    else if (!strcmp(s_output, "default"))
      write_level = 1;
    else if (!strcmp(s_output, "yes"))
      write_level = 2;
    BFT_FREE(s_output);
  }

  BFT_FREE(path);

  /* List of partitions to output */

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3,
                        "calculation_management",
                        "partitioning",
                        "partition_list");
  cs_xpath_add_function_text(&path);

  s_list = cs_gui_get_text_value(path);

  if (s_list != NULL) {
    char *p = strtok(s_list, " \t,;");
    while (p != NULL) {
      int np = atoi(p);
      if (np > 1) {
        BFT_REALLOC(add_parts, n_add_parts + 1, int);
        add_parts[n_add_parts] = np;
        n_add_parts += 1;
      }
      p = strtok(NULL, " \t,;");
    }
    BFT_FREE(s_list);
  }

  BFT_FREE(path);

  /* Set options */

  cs_partition_set_algorithm
    (CS_PARTITION_MAIN,
     a,
     rank_step,
     ignore_perio);

  cs_partition_set_write_level(write_level);

  if (n_add_parts > 0) {
    cs_partition_add_partitions(n_add_parts, add_parts);
    BFT_FREE(add_parts);
  }
}

/*-----------------------------------------------------------------------------
 * Define parallel IO settings.
 *----------------------------------------------------------------------------*/

void
cs_gui_parallel_io(void)
{
  int op_id;
  char  *path = NULL;

  int  rank_step = 0, block_size = -1;

  cs_file_mode_t op_mode[2] = {CS_FILE_MODE_READ, CS_FILE_MODE_WRITE};
  const char *op_name[2] = {"read_method", "write_method"};

  if (!cs_gui_file_is_loaded())
    return;

  /* Block IO read and write method */

  for (op_id = 0; op_id < 2; op_id++) {

    cs_file_access_t  m = CS_FILE_DEFAULT;
    char  *method_name = NULL;

    path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 3,
                          "calculation_management", "block_io", op_name[op_id]);
    cs_xpath_add_function_text(&path);

    method_name = cs_gui_get_text_value(path);

    if (method_name != NULL) {
      if (!strcmp(method_name, "default"))
        m = CS_FILE_DEFAULT;
      else if (!strcmp(method_name, "stdio serial"))
        m = CS_FILE_STDIO_SERIAL;
      else if (!strcmp(method_name, "stdio parallel"))
        m = CS_FILE_STDIO_PARALLEL;
      else if (!strcmp(method_name, "mpi independent"))
        m = CS_FILE_MPI_INDEPENDENT;
      else if (!strcmp(method_name, "mpi noncollective"))
        m = CS_FILE_MPI_NON_COLLECTIVE;
      else if (!strcmp(method_name, "mpi collective"))
        m = CS_FILE_MPI_COLLECTIVE;
#if defined(HAVE_MPI)
      cs_file_set_default_access(op_mode[op_id], m, MPI_INFO_NULL);
#else
      cs_file_set_default_access(op_mode[op_id], m);
#endif
      BFT_FREE(method_name);
    }

    BFT_FREE(path);

  }

#if defined(HAVE_MPI)

  /* Rank step and block buffer size */

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3,
                        "calculation_management", "block_io", "rank_step");
  cs_xpath_add_function_text(&path);
  cs_gui_get_int(path, &rank_step);

  BFT_FREE(path);

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3,
                        "calculation_management",
                        "block_io",
                        "min_block_size");
  cs_xpath_add_function_text(&path);
  cs_gui_get_int(path, &block_size);

  BFT_FREE(path);

  if (rank_step > 0 || block_size > -1) {
    int def_rank_step, def_block_size;
    cs_file_get_default_comm(&def_rank_step, &def_block_size, NULL, NULL);
    if (rank_step < 1)
      rank_step = def_rank_step;
    if (block_size < 0)
      block_size = def_block_size;
    cs_file_set_default_comm(rank_step, block_size, cs_glob_mpi_comm);
  }

#endif /* defined(HAVE_MPI) */
}

/*-----------------------------------------------------------------------------
 * Free memory: clean global private variables and libxml2 variables
 *----------------------------------------------------------------------------*/

void
cs_gui_clean_memory(void)
{
  int i;

  /* clean memory for global private structure vars */

  if (cs_glob_var->type != NULL)
    for (i = 0; i < cs_glob_var->nvar; i++)
      BFT_FREE(cs_glob_var->type[i]);
  BFT_FREE(cs_glob_var->type);

  if (cs_glob_var->head != NULL)
    for (i = 0; i < cs_glob_var->nvar; i++)
      BFT_FREE(cs_glob_var->head[i]);
  BFT_FREE(cs_glob_var->head);

  if (cs_glob_var->name != NULL)
    for (i = 0; i < cs_glob_var->nvar; i++)
      BFT_FREE(cs_glob_var->name[i]);
  BFT_FREE(cs_glob_var->name);

  if (cs_glob_var->label != NULL)
    for (i = 0; i < cs_glob_var->nscaus + cs_glob_var->nscapp; i++)
      BFT_FREE(cs_glob_var->label[i]);
  BFT_FREE(cs_glob_var->label);

  if (cs_glob_var->properties_name != NULL)
    for (i = 0; i < cs_glob_var->nprop; i++)
      BFT_FREE(cs_glob_var->properties_name[i]);
  BFT_FREE(cs_glob_var->properties_name);

  BFT_FREE(cs_glob_var->model);
  BFT_FREE(cs_glob_var->model_value);
  BFT_FREE(cs_glob_var->rtp);
  BFT_FREE(cs_glob_var->properties_ipp);
  BFT_FREE(cs_glob_var->propce);
  BFT_FREE(cs_glob_var);

  for (i = 0; i < cs_glob_label->_cs_gui_max_vars; i++)
    BFT_FREE(cs_glob_label->_cs_gui_var_name[i]);

  BFT_FREE(cs_glob_label->_cs_gui_var_name);
  BFT_FREE(cs_glob_label);

  mei_data_free();

  /* clean memory for xml document */

#if defined(HAVE_LIBXML2)
  if (xpathCtx != NULL) xmlXPathFreeContext(xpathCtx);
  if (node != NULL) xmlFreeNode(node);
#endif

  /* Shutdown libxml */

#if defined(HAVE_LIBXML2)
  xmlCleanupParser();
  xmlMemoryDump();
#endif
}

/*----------------------------------------------------------------------------
 * Logging output for MEI usage.
 *----------------------------------------------------------------------------*/

void
cs_gui_usage_log(void)
{
  double mei_wtime = cs_gui_get_mei_times();

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {
    double _wtime_loc = mei_wtime;
    MPI_Allreduce(&_wtime_loc, &mei_wtime, 1, MPI_DOUBLE, MPI_MAX,
                   cs_glob_mpi_comm);
  }

#endif

  if (mei_wtime > 0.0)
    bft_printf(_("\nTime elapsed defining values using MEI: %12.5f\n"),
        mei_wtime);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
