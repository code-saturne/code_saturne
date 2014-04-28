/*============================================================================
 * Management of the GUI parameters file: main parameters
 *============================================================================*/

/*
   This file is part of Code_Saturne, a general-purpose CFD tool.

   Copyright (C) 1998-2014 EDF S.A.

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
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_parameters.h"
#include "cs_partition.h"
#include "cs_prototypes.h"
#include "cs_timer.h"
#include "cs_physical_properties.h"

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
 * Return the label of a scalar field, given by its scalar Id
 *
 * parameters:
 *   id <-- field id
 *----------------------------------------------------------------------------*/

static inline const char *
_scalar_label(const int id)
{
  cs_field_t  *f = cs_field_by_id(id);
  return cs_field_get_label(f);
}

/*----------------------------------------------------------------------------
 * Turbulence model parameters.
 *
 * parameters:
 *   param                <--  name of the parameters
 *   keyword             -->   turbulence model parameter
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
 * Return the name of the related scalar if the scalar "name" is a variance
 *
 * parameter:
 *   name           <--  scalar name
 *----------------------------------------------------------------------------*/

static char *
cs_gui_scalar_variance(const char *const name)
{
  char *path = NULL;
  char *variance = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "additional_scalars");
  cs_xpath_add_element(&path, "variable");
  cs_xpath_add_test_attribute(&path, "name", name);
  cs_xpath_add_element(&path, "variance");
  cs_xpath_add_function_text(&path);

  variance = cs_gui_get_text_value(path);

  BFT_FREE(path);

  return variance;
}

/*----------------------------------------------------------------------------
 * Get thermal scalar model.
 *
 * return:
 *   value of itherm
 *----------------------------------------------------------------------------*/

int
gui_thermal_model(void)
{
  char *model_name = NULL;
  int   test = 0;

  model_name = cs_gui_get_thermophysical_model("thermal_scalar");

  if (cs_gui_strcmp(model_name, "off"))
    test = 0;
  else {
    if (cs_gui_strcmp(model_name, "enthalpy"))
      test = 20;
    else if (cs_gui_strcmp(model_name, "temperature_kelvin"))
      test = 11;
    else if (cs_gui_strcmp(model_name, "temperature_celsius"))
      test = 10;
    else if (cs_gui_strcmp(model_name, "potential_temperature"))
      test = 12;
    else if (cs_gui_strcmp(model_name, "liquid_potential_temperature"))
      test = 13;
    else if (cs_gui_strcmp(model_name, "total_energy"))
      test = 30;
    else
      bft_error(__FILE__, __LINE__, 0,
          _("Invalid thermal model: %s\n"), model_name);
  }

  BFT_FREE(model_name);

  return test;
}

/*-----------------------------------------------------------------------------
 * Return the name of the diffusion_coefficient property for a scalar
 *
 * parameters:
 *   scalar_index   <-- index of the scalar
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

/*----------------------------------------------------------------------------
 * Return the value of the choice attribute for material, method, ...
 *
 * parameters:
 *   name        <--  name of the property
 *----------------------------------------------------------------------------*/

static char*
_thermal_table_choice(const char *const name)
{
  char *path   = NULL;
  char *choice = NULL;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "fluid_properties");
  cs_xpath_add_element(&path, name);
  cs_xpath_add_attribute(&path, "choice");
  choice = cs_gui_get_attribute_value(path);
  BFT_FREE(path);
  return choice;
}

/*----------------------------------------------------------------------------
 * Return the value of the choice attribute from a property name.
 *
 * parameters:
 *   property_name        <--  name of the property
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
 * Return 0 if default value is needed
 *
 * parameters:
 *   name        <--  name of the property
 *----------------------------------------------------------------------------*/

static int
_thermal_table_needed(const char *const name)
{
  int choice = 0;

  char *prop_choice = _properties_choice(name);
  if (cs_gui_strcmp(prop_choice, "thermal_law"))
    choice = 1;
  BFT_FREE(prop_choice);
  return choice;
}

/*-----------------------------------------------------------------------------
 * use MEI for physical property
 *----------------------------------------------------------------------------*/
static void
_physical_property(const char *const param,
                   const char *const symbol,
                   const cs_int_t  *ncel,
                   const cs_int_t  *ncelet,
                   const cs_int_t  *itherm,
                   const cs_int_t  *iscalt,
                   const cs_int_t  *icp,
                   const cs_real_t *p0,
                   const cs_real_t *ro0,
                   const cs_real_t *cp0,
                   const cs_real_t *viscl0,
                   const cs_real_t *visls0,
                   double values[])
{
  cs_var_t  *vars = cs_glob_var;

  int user_law = 0;
  char *law = NULL;
  double time0;
  char *path = NULL;
  mei_tree_t *ev_law = NULL;
  int i, iel;

  char *prop_choice = _properties_choice(param);

  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->cell_cen;

  if (cs_gui_strcmp(prop_choice, "variable"))
    user_law = 1;

  if (user_law)
  {
    /* search the formula for the law */
    path = cs_xpath_short_path();
    cs_xpath_add_element(&path, "property");
    cs_xpath_add_test_attribute(&path, "name", param);
    cs_xpath_add_element(&path, "formula");
    cs_xpath_add_function_text(&path);

    law = cs_gui_get_text_value(path);
    BFT_FREE(path);

    if (law != NULL) {
      time0 = cs_timer_wtime();

      ev_law = mei_tree_new(law);

      mei_tree_insert(ev_law, "x", 0.0);
      mei_tree_insert(ev_law, "y", 0.0);
      mei_tree_insert(ev_law," z", 0.0);

      mei_tree_insert(ev_law, "p0", *p0);

      if (cs_gui_strcmp(param, "density"))
      {
        mei_tree_insert(ev_law, "rho0", *ro0);
      }
      else if (cs_gui_strcmp(param, "molecular_viscosity")) {
        mei_tree_insert(ev_law, "rho0", *ro0);
        mei_tree_insert(ev_law, "mu0", *viscl0);
        mei_tree_insert(ev_law, "rho", 0.0);
        if (cs_gui_strcmp(vars->model, "compressible_model"))
          mei_tree_insert(ev_law, "t0", 0.0);
      }
      else if (cs_gui_strcmp(param, "specific_heat")) {
        mei_tree_insert(ev_law, "cp0", *cp0);
      }
      else if (cs_gui_strcmp(param, "thermal_conductivity")) {
        /* for the Temperature, the diffusivity factor is not divided by Cp */
        if (*itherm != 1)
        {
          mei_tree_insert(ev_law, "lambda0", visls0[*iscalt-1]*(*cp0));
        }
        else {
          mei_tree_insert(ev_law, "lambda0", visls0[*iscalt-1]);
        }
      }

      for (int f_id2 = 0; f_id2 < cs_field_n_fields(); f_id2++) {
        const cs_field_t  *f2 = cs_field_by_id(f_id2);
        if (f2->type & CS_FIELD_USER)
          mei_tree_insert(ev_law, _scalar_label(f_id2), 0.0);
      }

      char *buff = NULL;
      int mdl = gui_thermal_model();

      if (mdl < 20) {
        BFT_MALLOC(buff, 12, char);
        strcpy(buff, "temperature");
      }
      else if (mdl < 30) {
        BFT_MALLOC(buff, 9, char);
        strcpy(buff, "enthalpy");
      }
      else {
        BFT_MALLOC(buff, 13, char);
        strcpy(buff, "total_energy");
      }
      cs_field_t *fth = cs_field_by_name_try(buff);
      int fth_id = fth->id;

      mei_tree_insert(ev_law, _scalar_label(fth_id), 0.0);

      /* try to build the interpreter */

      if (mei_tree_builder(ev_law))
        bft_error(__FILE__, __LINE__, 0,
                  _("Error: can not interpret expression: %s\n"), ev_law->string);

      if (mei_tree_find_symbol(ev_law, symbol))
        bft_error(__FILE__, __LINE__, 0,
                  _("Error: can not find the required symbol: %s\n"), symbol);

      /* for each cell, update the value of the table of symbols for each scalar
         (including the thermal scalar), and evaluate the interpreter */

      cs_field_t *c_cp = CS_F_(cp);
      cs_field_t *c_rho = CS_F_(rho);
      cs_field_t *c_t = CS_F_(t);

      for (iel = 0; iel < *ncel; iel++)
      {
        mei_tree_insert(ev_law, "x", cell_cen[iel][0]);
        mei_tree_insert(ev_law, "y", cell_cen[iel][1]);
        mei_tree_insert(ev_law, "z", cell_cen[iel][2]);
        for (int f_id = 0; f_id < cs_field_n_fields(); f_id++) {
          cs_field_t  *f = cs_field_by_id(f_id);
          if (f->type & CS_FIELD_USER)
            mei_tree_insert(ev_law,
                           _scalar_label(f_id),
                            f->val[iel]);
        }

        mei_tree_insert(ev_law,
                       _scalar_label(fth_id),
                        fth->val[iel]);

        if (cs_gui_strcmp(param, "molecular_viscosity")) {
          mei_tree_insert(ev_law, "rho", c_rho->val[iel]);
          if (cs_gui_strcmp(vars->model, "compressible_model"))
            mei_tree_insert(ev_law, "T", c_t->val[iel]);
          }

        mei_evaluate(ev_law);

        if (cs_gui_strcmp(param, "thermal_conductivity")) {
          if (*itherm == 1)
            values[iel] = mei_tree_lookup(ev_law, symbol);
          else if (*icp > 0)
            values[iel] = mei_tree_lookup(ev_law, symbol) / c_cp->val[iel];
          else
            values[iel] = mei_tree_lookup(ev_law, symbol) / *cp0;
        }
        else {
          values[iel] = mei_tree_lookup(ev_law, symbol);
        }
      }

      mei_tree_destroy(ev_law);
      BFT_FREE(buff);

      cs_gui_add_mei_time(cs_timer_wtime() - time0);
    }
  }
  else if (cs_gui_strcmp(prop_choice, "thermal_law")) {
    cs_phys_prop_type_t property;
    cs_field_t *c_prop;

    if (cs_gui_strcmp(param, "density")) {
      property = CS_PHYS_PROP_DENSITY;
      c_prop = CS_F_(rho);
    }
    else if (cs_gui_strcmp(param, "molecular_viscosity")) {
      property = CS_PHYS_PROP_DYNAMIC_VISCOSITY;
      c_prop = CS_F_(mu);
    }
    else if (cs_gui_strcmp(param, "specific_heat")) {
      property = CS_PHYS_PROP_ISOBARIC_HEAT_CAPACITY;
      c_prop = CS_F_(cp);
    }
    else if (cs_gui_strcmp(param, "thermal_conductivity")) {
      property = CS_PHYS_PROP_THERMAL_CONDUCTIVITY;

      cs_field_t *_th_f[] = {CS_F_(t), CS_F_(h), CS_F_(energy)};

      for (i = 0; i < 3; i++)
        if (_th_f[i]) {
          if ((_th_f[i])->type & CS_FIELD_VARIABLE) {
            int k = cs_field_key_id("scalar_diffusivity_id");
            int cond_diff_id = cs_field_get_key_int(_th_f[i], k);
            if (cond_diff_id > -1)
              c_prop = cs_field_by_id(cond_diff_id);
            break;
          }
        }
    }
    else {
      bft_error(__FILE__, __LINE__, 0,
                _("Error: can not use evaluate property: %s\n"), prop_choice);
    }

    cs_field_t *c_pres = CS_F_(p);

    cs_real_t *ptot;
    BFT_MALLOC(ptot, *ncelet, cs_real_t);
    for (iel = 0; iel < *ncelet; iel++)
      ptot[iel] = c_pres->val[iel] + *p0;

    cs_field_t *_th_f[] = {CS_F_(t), CS_F_(h), CS_F_(energy)};

    for (i = 0; i < 3; i++)
      if (_th_f[i]) {
        if ((_th_f[i])->type & CS_FIELD_VARIABLE) {
          cs_phys_prop_compute(property,
                              *ncel, ptot, _th_f[i]->val, c_prop->val);
          break;
        }
      }
    BFT_FREE(ptot);
  }
  BFT_FREE(prop_choice);
  BFT_FREE(law);

}

/*-----------------------------------------------------------------------------
 * use MEI for compressible physical property
 *----------------------------------------------------------------------------*/

static void
_compressible_physical_property(const char *const param,
                                     const char *const symbol,
                                     const cs_int_t   idx,
                                     const cs_int_t  *ncel,
                                     const cs_int_t  *itempk,
                                     const cs_real_t *p0,
                                     const cs_real_t *t0,
                                     const cs_real_t *ro0,
                                     const cs_real_t *visls0,
                                     const cs_real_t *viscv0)
{
  int variable = 0;
  char *law = NULL;
  double time0;
  char *path = NULL;
  mei_tree_t *ev_law = NULL;

  char *prop_choice = _properties_choice(param);
  int n_fields = cs_field_n_fields();

  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->cell_cen;

  if (cs_gui_strcmp(prop_choice, "variable"))
    variable = 1;
  BFT_FREE(prop_choice);

  if (variable) {
    /* search the formula for the law */
    path = cs_xpath_short_path();
    cs_xpath_add_element(&path, "property");
    cs_xpath_add_test_attribute(&path, "name", param);
    cs_xpath_add_element(&path, "formula");
    cs_xpath_add_function_text(&path);

    law = cs_gui_get_text_value(path);

    BFT_FREE(path);

    if (law != NULL) {
      time0 = cs_timer_wtime();

      ev_law = mei_tree_new(law);

      mei_tree_insert(ev_law, "x", 0.0);
      mei_tree_insert(ev_law, "y", 0.0);
      mei_tree_insert(ev_law," z", 0.0);

      mei_tree_insert(ev_law, "p0", *p0);
      mei_tree_insert(ev_law, "t0", *t0);

      if (cs_gui_strcmp(param, "thermal_conductivity")) {
        mei_tree_insert(ev_law, "lambda0", visls0[*itempk -1]);
        mei_tree_insert(ev_law, "rho0", *ro0);
      }
      else if (cs_gui_strcmp(param, "volume_viscosity")) {
        mei_tree_insert(ev_law,"viscv0", *viscv0);
        mei_tree_insert(ev_law,"T", 0.);
      }

      if (cs_gui_strcmp(param, "thermal_conductivity")) {
        for (int f_id2 = 0; f_id2 < n_fields; f_id2++) {
          const cs_field_t  *f2 = cs_field_by_id(f_id2);
          if (f2->type & CS_FIELD_USER)
            mei_tree_insert(ev_law, f2->name, 0.0);
        }
      }

      /* try to build the interpreter */

      if (mei_tree_builder(ev_law))
        bft_error(__FILE__, __LINE__, 0,
                  _("Error: can not interpret expression: %s\n"), ev_law->string);

      if (mei_tree_find_symbol(ev_law, symbol))
        bft_error(__FILE__, __LINE__, 0,
                  _("Error: can not find the required symbol: %s\n"), symbol);

      /* for each cell, update the value of the table of symbols for each scalar
         (including the thermal scalar), and evaluate the interpreter */

      cs_field_t *c = cs_field_by_id(idx);
      char *buff = NULL;
      int mdl = gui_thermal_model();

      if (mdl < 20) {
        BFT_MALLOC(buff, 12, char);
        strcpy(buff, "temperature");
      }
      else if (mdl < 30) {
        BFT_MALLOC(buff, 9, char);
        strcpy(buff, "enthalpy");
      }
      else {
        BFT_MALLOC(buff, 13, char);
        strcpy(buff, "total_energy");
      }

      cs_field_t *f = cs_field_by_name_try(buff);

      for (int iel = 0; iel < *ncel; iel++) {
        mei_tree_insert(ev_law, "x", cell_cen[iel][0]);
        mei_tree_insert(ev_law, "y", cell_cen[iel][1]);
        mei_tree_insert(ev_law, "z", cell_cen[iel][2]);
        if (cs_gui_strcmp(param, "thermal_conductivity")) {
            for (int f_id2 = 0; f_id2 < n_fields; f_id2++) {
              const cs_field_t  *f2 = cs_field_by_id(f_id2);
              if (f2->type & CS_FIELD_USER)
                mei_tree_insert(ev_law,
                                f2->name,
                                f2->val[iel]);
            }
        }

        mei_tree_insert(ev_law, buff, f->val[iel]);

        mei_evaluate(ev_law);
        c->val[iel] = mei_tree_lookup(ev_law, symbol);
      }
      mei_tree_destroy(ev_law);
      BFT_FREE(buff);

      cs_gui_add_mei_time(cs_timer_wtime() - time0);
    }
  }
}

/*-----------------------------------------------------------------------------
 * Return the value of choice for user scalar's property
 *
 * parameters:
 *   scalar_num     <-- number of scalar
 *   choice         <-- choice for property
 *----------------------------------------------------------------------------*/

static int
cs_gui_scalar_properties_choice(const int scalar_num, int *const choice)
{
  char *path = NULL;
  char *buff = NULL;
  int   ichoice;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "additional_scalars");
  cs_xpath_add_element_num(&path, "variable", scalar_num);
  cs_xpath_add_element(&path, "property");
  cs_xpath_add_attribute(&path, "choice");

  buff = cs_gui_get_attribute_value(path);

  if (buff == NULL) {
    ichoice = 0;

  } else {
    ichoice = 1;

    if (cs_gui_strcmp(buff, "variable"))
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
 *   num_sca  <-- number of scalar
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
  cs_xpath_add_element_num(&path, "variable", num_sca);
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
 *   keyword         -->  if 1 unsteady management else steady management
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
 *   param           <--  steady parameter
 *   keyword         -->  new value for the steady parameter
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
 *   param              <--  time parameter
 *   keyword            -->  new value of the time parameter
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
 *   param               <--  restart parameter
 *   keyword            <<--  new value of the restart parameter
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
 *   variable_name  <-- name of variable
 *   value_type     <-- name of numerical parameter parkup
 *   value          --> value of numerical parameter
 *----------------------------------------------------------------------------*/

static void
cs_gui_variable_value(const char   *const variable_name,
                      const char   *const value_type,
                            double *const value)
{
  char  *path = NULL;
  double result;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "variable");
  cs_xpath_add_test_attribute(&path, "name", variable_name);
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
 *   path          <-- path for xpath query
 *   child         <-- child markup
 *   keyword      -->  value of attribute node
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
 *   name          <--  name of the variable markup
 *   child         <--  child markup
 *   keyword      -->   value of attribute node contained in the child markup
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
 *   name          <--  name of the variable markup
 *   child         <--  child markup
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

/*-----------------------------------------------------------------------------
 * Modify double numerical parameters.
 *
 * parameters:
 *   param               <--  label of the numerical parameter
 *   keyword            <<--  value of the numerical parameter
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
 *   param               <--  label of the numerical parameter
 *   keyword            <<--  value of the numerical parameter
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
 *   param               <--  gravity parameter (GX, GY, GZ)
 *   keyword            <<--  new value of the gravity parameter
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
 *   param               <--  coriolis parameter (OMEGAX, OMEGAY, OMEGAZ)
 *   keyword            <<--  new value of the coriolis parameter
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
 * Get the value of the choice attribute from a property markup.
 * Return 1 if the xpath request has succeeded, 0 otherwise.
 *
 * parameters:
 *   property_name        <--  name of the property
 *   choice              -->   value of the attribute choice
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
    if (cs_gui_strcmp(buff, "variable") || cs_gui_strcmp(buff, "thermal_law"))
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
 *   param                <--  name of the parameters
 *   keyword             -->   turbulence model parameter
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
 *   zone_id        <--  zone number
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
 *   id           <--  time average number (imom)
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
 * Return the component of variables or properties or scalar for a given time average
 *
 * parameters:
 *   id           -->  number of 1D profile
 *   nm           -->  number of the variable name of the idst 1D profile
 *----------------------------------------------------------------------------*/

static int _get_time_average_component(const int id, const int nm)
{
  char *path = NULL;
  char *comp = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "analysis_control", "time_averages");
  cs_xpath_add_element_num(&path, "time_average", id);
  cs_xpath_add_element_num(&path, "var_prop", nm);
  cs_xpath_add_attribute(&path, "component");

  comp = cs_gui_get_attribute_value(path);
  if (comp == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Invalid xpath: %s\n component not found"), path);
  BFT_FREE(path);

  int compId = atoi(comp);

  BFT_FREE(comp);

  return compId;
}

/*----------------------------------------------------------------------------
 * Get value of a parameter for a given time average.
 *
 * parameters:
 *   id              <--  time average number (imom)
 *   param           <--  name of the parameter
 *   data           -->   value of the parameter
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
 *   id           <--  time average number (imom)
 *   nb           <--  variable or property number
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
 *   id              <--  time average number (imom)
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
 *   variable   <-- name of variable
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
 *   property_name        <--  name of the property
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
 *   kw                   <--  keyword: 'label' or 'name'
 *   scalar_num          -->   number of the searching scalar
 *----------------------------------------------------------------------------*/

static char *_scalar_name_label(const char *kw, const int scalar_num)
{
  char *path = NULL;
  char *str  = NULL;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "additional_scalars");
  cs_xpath_add_element_num(&path, "variable", scalar_num);
  cs_xpath_add_attribute(&path, kw);

  str = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return str;
}

/*-----------------------------------------------------------------------------
 * Return the name tor thermal scalar.
 *
 * parameters:
 *   kw                   <--  scalar name
 *----------------------------------------------------------------------------*/

static char *_thermal_scalar_name_label(const char *kw)
{
  char *path = NULL;
  char *str  = NULL;

  path = cs_xpath_short_path();
  cs_xpath_add_elements(&path, 3,
                        "thermophysical_models",
                        "thermal_scalar",
                        "variable");
  cs_xpath_add_attribute(&path, kw);

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
 *   ith_zone        <--  id of volumic zone
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
 *   zone_id      <--  volumic zone id
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
 *   variable_name    <--  name of variable
 *   zone_id          <--  id of volumic zone
 *   initial_value    -->  initial value
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
 *   id           <--  number of order in list of 1D profile
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
 * Return the component of variables or properties or scalar for 1D profile
 *
 * parameters:
 *   id           -->  number of 1D profile
 *   nm           -->  number of the variable name of the idst 1D profile
 *----------------------------------------------------------------------------*/

static int _get_profile_component(const int id, const int nm)
{
  char *path = NULL;
  char *comp = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "analysis_control", "profiles");
  cs_xpath_add_element_num(&path, "profile", id + 1);
  cs_xpath_add_element_num(&path, "var_prop", nm + 1);
  cs_xpath_add_attribute(&path, "component");

  comp = cs_gui_get_attribute_value(path);
  if (comp == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Invalid xpath: %s\n component not found"), path);
  BFT_FREE(path);

  int compId = atoi(comp);

  BFT_FREE(comp);

  return compId;
}

/*----------------------------------------------------------------------------
 * Get number of variables or properties or scalar for 1D profile
 *
 * parameters:
 *   id           <--  number of 1D profile
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
 *   id           <--  number of 1D profile
 *   nm           <--  number of the variable name of the idst 1D profile
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
 *   id           <--  number of 1D profile
 *   nm           <--  number of the variable name of the idst 1D profile
 *----------------------------------------------------------------------------*/

static char *_get_profile_label_name(const int id, const int nm)
{
  char *path = NULL;
  char *name = NULL;
  char *label = NULL;

  name = _get_profile_name(id, nm);
  int idim = _get_profile_component(id, nm);

  for (int f_id = 0; f_id < cs_field_n_fields(); f_id++) {
    const cs_field_t *f = cs_field_by_id(f_id);
    if (cs_gui_strcmp(name,  f->name)) {
      if (f->type & CS_FIELD_VARIABLE)
      {
        label = cs_gui_variable_label(name);
        if (f->dim > 1)
        {
          int len = strlen(label) + 4;
          char *tmp = NULL;
          char *snumpp = NULL;
          BFT_MALLOC(snumpp, 2, char);
          sprintf(snumpp, "%1.1i", idim);
          BFT_MALLOC(tmp, len, char);
          strcpy(tmp, label);
          strcat(tmp, "[");
          strcat(tmp, snumpp);
          strcat(tmp, "]");
          BFT_FREE(label);
          BFT_MALLOC(label, len, char);
          strcpy(label, tmp);
          BFT_FREE(snumpp);
          BFT_FREE(tmp);
        }
      }
      else if (f->type & CS_FIELD_PROPERTY)
      {
        label = cs_gui_properties_label(name);
        if (f->dim > 1)
        {
          int len = strlen(label) + 4;
          char *tmp = NULL;
          char *snumpp = NULL;
          BFT_MALLOC(snumpp, 2, char);
          sprintf(snumpp, "%1.1i", idim);
          BFT_MALLOC(tmp, len, char);
          strcpy(tmp, label);
          strcat(tmp, "[");
          strcat(tmp, snumpp);
          strcat(tmp, "]");
          BFT_FREE(label);
          BFT_MALLOC(label, len, char);
          strcpy(label, tmp);
          BFT_FREE(snumpp);
          BFT_FREE(tmp);
        }
      }
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
 *   id           <--  number of 1D profile
 *    x          -->   name of the coordinate (x1, y1, z1, x2, y2, z2)
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
 *   id           <--  number of average
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
 *   id           <--  number of 1D profile
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

  cs_glob_var->model            = NULL;
  cs_glob_var->model_value      = NULL;

  BFT_MALLOC(cs_glob_label, 1, cs_label_t);

  cs_glob_label->_cs_gui_max_vars = 0;
  cs_glob_label->_cs_gui_last_var = 0;
  cs_glob_label->_cs_gui_var_name = NULL;
}

/*----------------------------------------------------------------------------
 * Thermal model.
 *
 * Fortran Interface:
 *
 * subroutine csther (itherm, itpscl)
 * *****************
 *
 * integer          itherm  --> thermal model
 * integer          itpscl  --> temperature scale if itherm = 1
 *----------------------------------------------------------------------------*/


void CS_PROCF (csther, CSTHER) (int  *itherm,
                                int  *itpscl)
{
  switch(gui_thermal_model()) {
  case 10:
    *itherm = 1;
    *itpscl = 2;
    break;
  case 11:
    *itherm = 1;
    *itpscl = 1;
    break;
  case 12:
    *itherm = 1;
    *itpscl = 2;
    break;
  case 13:
    *itherm = 1;
    *itpscl = 2;
    break;
  case 20:
    *itherm = 2;
    *itpscl = 1;
    break;
  case 30:
    *itherm = 3;
    *itpscl = 1;
    break;
  default:
    *itherm = 0;
    *itpscl = 0;
    break;
  }
}

/*----------------------------------------------------------------------------
 * Turbulence model.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSTURB (ITURB, IDEUCH, IGRAKE, IGRAKI, XLOMLG)
 * *****************
 *
 * INTEGER          ITURB   -->   turbulence model
 * INTEGER          IDEUCH  -->   wall law treatment
 * INTEGER          IGRAKE  -->   k-eps gravity effects
 * INTEGER          IGRAKI  -->   Rij-eps gravity effects
 * DOUBLE PRECISION XLOMLG  -->   mixing_length_scale
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
 * INTEGER          ICP     -->   specific heat variable or constant indicator
 *----------------------------------------------------------------------------*/

void CS_PROCF (cscpva, CSCPVA) (int *const icp)
{
  int choice;

  if (cs_gui_properties_choice("specific_heat", &choice))
    *icp = choice;

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
 * INTEGER          IVISCV     -->   specific heat variable or constant indicator
 *----------------------------------------------------------------------------*/

void CS_PROCF (csvvva, CSVVVA) (int *const iviscv)
{
  int choice;

  if (cs_gui_properties_choice("volume_viscosity", &choice))
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
 * INTEGER          NSCAUS     -->   user scalars number
 *----------------------------------------------------------------------------*/

void CS_PROCF (csnsca, CSNSCA) (int *const nscaus)
{
  *nscaus = cs_gui_get_tag_number("/additional_scalars/variable", 1);

#if _XML_DEBUG_
  bft_printf("==>CSNSCA\n");
  bft_printf("--number of user scalars: %i\n", *nscaus);
#endif
}

/*----------------------------------------------------------------------------
 * User scalars labels
 *
 * Fortran Interface:
 *
 * SUBROUTINE UISCAU
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiscau, UISCAU) (void)
{
  int   i     = 0;
  char *label = NULL;

  const int n_fields = cs_field_n_fields();
  const int keysca = cs_field_key_id("scalar_id");
  const int keylbl = cs_field_key_id("label");

#if _XML_DEBUG_
  bft_printf("==>UISCAU\n");
#endif

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    i = cs_field_get_key_int(f, keysca) - 1;
    if (f->type & CS_FIELD_USER) {
      label = _scalar_name_label("label", i+1);
      cs_field_set_key_str(f, keylbl, label);
#if _XML_DEBUG_
      bft_printf("--label of scalar[%i]: %s\n", i, label);
#endif
      BFT_FREE(label);
    }
  }
}

/*----------------------------------------------------------------------------
 * User thermal scalar.
 *
 * Fortran Interface:
 *
 * SUBROUTINE UITHSC
 * *****************
 *
 * INTEGER          ISCALT     -->   thermal scalars number
 *----------------------------------------------------------------------------*/

void CS_PROCF (uithsc, UITHSC) (int *const iscalt)
{
  cs_var_t  *vars = cs_glob_var;
  char *label = NULL;

  const int n_fields = cs_field_n_fields();
  const int keysca = cs_field_key_id("scalar_id");
  const int keylbl = cs_field_key_id("label");

  label = _thermal_scalar_name_label("label");

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    int i = cs_field_get_key_int(f, keysca) - 1;
    if (i == *iscalt - 1) {
#if _XML_DEBUG_
      bft_printf("--label of thermal scalar: %s\n", label);
#endif
      cs_field_set_key_str(f, keylbl, label);
      break;
    }
  }

  BFT_FREE(label);

  BFT_MALLOC(vars->model, strlen("thermal_scalar")+1, char);
  strcpy(vars->model, "thermal_scalar");
}

/*----------------------------------------------------------------------------
 * User scalars which are variance.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSISCA (ISCAVR)
 * *****************
 *
 * INTEGER          ISCAVR  -->  user scalars variance array
 * integer          itherm  <--  type of thermal model
 *----------------------------------------------------------------------------*/

void CS_PROCF (csisca, CSISCA) (      int *const iscavr,
                                      int *const itherm,
                                const int *const iscapp)
{
  char *variance = NULL;

  const int keysca = cs_field_key_id("scalar_id");

#if _XML_DEBUG_
  bft_printf("==>CSISCA\n");
#endif

  for (int f_id = 0; f_id < cs_field_n_fields(); f_id++) {
    cs_field_t  *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_USER) {
      variance = cs_gui_scalar_variance(f->name);

      if (variance != NULL) {
        int i = cs_field_get_key_int(f, keysca) - 1;
        for (int f_id2 = 0; f_id2 < cs_field_n_fields(); f_id2++) {
          cs_field_t  *f2 = cs_field_by_id(f_id2);

          if (f2->type & CS_FIELD_USER) {
            if (cs_gui_strcmp(variance, _scalar_label(f2->id))) {
              if (f_id == f_id2)
                bft_error(__FILE__, __LINE__, 0,
                          _("Scalar: %s and its variance: %s are the same.\n"),
                          f->name, f2->name);
              int j = cs_field_get_key_int(f2, keysca);
              iscavr[i] = j;
            }
          }
        }

        if (*itherm && iscavr[i] == 0) {
          for (int f_id2 = 0; f_id2 < cs_field_n_fields(); f_id2++) {
            cs_field_t  *f2 = cs_field_by_id(f_id2);
            if (f2->type & CS_FIELD_VARIABLE && !(f2->type & CS_FIELD_USER)) {
              if (cs_gui_strcmp(variance, _scalar_label(f2->id))) {
                int j = cs_field_get_key_int(f2, keysca);
                iscavr[i] = j;
              }
            }
          }
        }
      }
      BFT_FREE(variance);

#if _XML_DEBUG_
      bft_printf("--iscavr[%i] = %i \n", f_id, iscavr[f_id]);
#endif
    }
  }

  return;
}

/*----------------------------------------------------------------------------
 * Constant or variable indicator for the user scalar laminar viscosity.
 *
 * Fortran Interface:
 *
 * subroutine csivis (iscavr, ivisls, iscalt, itherm, itempk)
 * *****************
 *
 * integer          iscavr  <-->  number of the related variance if any
 * integer          ivisls  <--   indicator for the user scalar viscosity
 * integer          iscalt  <-->  number of the user thermal scalar if any
 * integer          itherm  <-->  type of thermal model
 * integer          itempk   -->  rtp index for temperature (in K)
 *----------------------------------------------------------------------------*/

void CS_PROCF (csivis, CSIVIS) (int *const iscavr,
                                int *const ivisls,
                                int *const iscalt,
                                int *const itherm,
                                int *const itempk)
{
  int choice1, choice2;
  int test1, test2;

  cs_var_t  *vars = cs_glob_var;

  const int keysca = cs_field_key_id("scalar_id");

#if _XML_DEBUG_
  bft_printf("==>CSIVIS\n");
#endif

  if (vars->model != NULL)
    if (*itherm) {
      test1 = cs_gui_properties_choice("thermal_conductivity", &choice1);
      test2 = cs_gui_properties_choice("specific_heat", &choice2);

      if (test1 && test2) {
        if (choice1 || choice2)
          ivisls[*iscalt-1] = 1;
        else
          ivisls[*iscalt-1] = 0;
      }
#if _XML_DEBUG_
      bft_printf("--ivisls[%i] = %i\n", i, ivisls[i]);
#endif
    }

  int n_fields = cs_field_n_fields();
  for (int f_id = 0; f_id < n_fields; f_id++) {
    const cs_field_t  *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_USER) {
      int i = cs_field_get_key_int(f, keysca) - 1;
      if (iscavr[i] <= 0 ) {
        if (cs_gui_scalar_properties_choice(i+1, &choice1))
          if (*iscalt != i+1)
            ivisls[i] = choice1;
      }
#if _XML_DEBUG_
    bft_printf("--ivisls[%i] = %i\n", i, ivisls[i]);
#endif
    }
  }

  if (cs_gui_strcmp(vars->model, "compressible_model"))
  {
    ivisls[*itempk -1] = 0;

    char *prop_choice = _properties_choice("thermal_conductivity");
    if (cs_gui_strcmp(prop_choice, "variable"))
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
 * INTEGER          IDTVAR  -->   fixed or variable time step
 *----------------------------------------------------------------------------*/

void CS_PROCF (csidtv, CSIDTV) (int *const idtvar)
{
  double param;
  int steady = 0;

  cs_gui_get_steady_status(&steady);
  if (steady) {
    char *algo_choice = cs_gui_velocity_pressure_algo_choice();
    if (cs_gui_strcmp(algo_choice, "simple"))
      *idtvar = -1;
    else
      *idtvar = 2;
    BFT_FREE(algo_choice);
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
 * INTEGER          IPHYDR  -->   hydrostatic pressure
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
 * INTEGER          icfgrp  -->   hydrostatic equilibrium
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
 * Restart parameters.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSISUI (NTSUIT, ILEAUX, ICCVFG)
 * *****************
 *
 * INTEGER          NTSUIT  -->   checkpoint frequency
 * INTEGER          ILEAUX  -->   restart with auxiliary
 * INTEGER          ICCFVG  -->   restart with frozen field
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
 * INTEGER          INPDT0  -->   zero time step
 * INTEGER          IPTLTO  -->   thermal time step control
 * INTEGER          NTMABS  -->   iterations numbers
 * INTEGER          IDTVAR  -->   time steps'options
 * DOUBLE PRECISION DTREF   -->   time step
 * DOUBLE PRECISION DTMIN   -->   minimal time step
 * DOUBLE PRECISION DTMAX   -->   maximal time step
 * DOUBLE PRECISION COUMAX  -->   maximal courant number
 * DOUBLE PRECISION FOUMAX  -->   maximal fournier number
 * DOUBLE PRECISION VARRDT  -->   max time step variation between 2 iterations
 * DOUBLE PRECISION RELXST  -->   relaxation coefficient id idtvar = -1
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
 * Treatment of local numerical aspects:
 *     BLENCV, ISCHCV, ISSTPC, IRCFLU, CDTVAR, NITMAX, EPSILO
 *----------------------------------------------------------------------------*/

void CS_PROCF (uinum1, UINUM1) (double *const blencv,
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
  double tmp;
  char* algo_choice = NULL;
  int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  int var_key_id = cs_field_key_id("variable_id");
  cs_var_cal_opt_t var_cal_opt;

  /* 1) variables from velocity_pressure and turbulence */
  /* 1-a) for pressure */
  cs_field_t *c_pres = cs_field_by_name("pressure");
  cs_field_get_key_struct(c_pres, key_cal_opt_id, &var_cal_opt);
  int j = cs_field_get_key_int(c_pres, var_key_id) -1;

  cs_gui_variable_value(c_pres->name, "solver_precision", &epsilo[j]);
  tmp = (double) nitmax[j];
  cs_gui_variable_value(c_pres->name, "max_iter_number", &tmp);
  nitmax[j] = (int) tmp;

  imgr[j] = 0;

  algo_choice = cs_gui_variable_choice(c_pres->name, "solver_choice");
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
  cs_gui_variable_value(c_pres->name, "rhs_reconstruction", &tmp);
  nswrsm[j] = (int) tmp;
  BFT_FREE(algo_choice);

  // Set Field calculation options in the field structure
  var_cal_opt.epsilo = epsilo[j];
  // TODO add nitmax, imgr, iresol
  var_cal_opt.nswrsm = nswrsm[j];
  cs_field_set_key_struct(c_pres, key_cal_opt_id, &var_cal_opt);

  /* 1-b) for the other variables */
  int n_fields = cs_field_n_fields();
  for (int f_id = 0; f_id < n_fields; f_id++) {
    const cs_field_t  *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE && !cs_gui_strcmp(f->name, "pressure")) {
      j = cs_field_get_key_int(f, var_key_id) -1;
      cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);

      cs_gui_variable_value(f->name, "blending_factor", &blencv[j]);
      cs_gui_variable_value(f->name, "solver_precision", &epsilo[j]);

      imgr[j] = 0;

      algo_choice = cs_gui_variable_choice(f->name, "solver_choice");

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

      // only for nscaus and model scalar
      cs_gui_variable_value(f->name, "time_step_factor", &cdtvar[j]);

      tmp = (double) nitmax[j];
      cs_gui_variable_value(f->name, "max_iter_number", &tmp);
      nitmax[j] = (int) tmp;
      cs_gui_variable_attribute(f->name, "order_scheme", &ischcv[j]);
      cs_gui_variable_attribute(f->name, "slope_test", &isstpc[j]);
      cs_gui_variable_attribute(f->name, "flux_reconstruction", &ircflu[j]);
      tmp = (double) nswrsm[j];
      cs_gui_variable_value(f->name, "rhs_reconstruction", &tmp);
      nswrsm[j] = (int) tmp;

      // Set Field calculation options in the field structure
      var_cal_opt.blencv = blencv[j];
      var_cal_opt.epsilo = epsilo[j];
      // TODO add nitmax, imgr, iresol, cdtvar
      var_cal_opt.nswrsm = nswrsm[j];
      cs_field_set_key_struct(f, key_cal_opt_id, &var_cal_opt);
    }
  }

#if _XML_DEBUG_
  bft_printf("==>UINUM1\n");
  for (f_id = 0; f_id < n_fields; f_id++) {
    const cs_field_t  *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      j = cs_field_get_key_int(f, var_key_id) -1;
      bft_printf("-->variable[%i] = %s\n", i, f->name);
      bft_printf("--blencv = %f\n", blencv[i]);
      bft_printf("--epsilo = %g\n", epsilo[i]);
      bft_printf("--cdtvar = %g\n", cdtvar[i]);
      bft_printf("--nitmax = %i\n", nitmax[i]);
      bft_printf("--ischcv = %i\n", ischcv[i]);
      bft_printf("--isstpc = %i\n", isstpc[i]);
      bft_printf("--ircflu = %i\n", ircflu[i]);
      bft_printf("--nswrsm = %i\n", nswrsm[i]);
      bft_printf("--imgr = %i\n"  , imgr[i]);
      bft_printf("--iresol = %i\n", iresol[i]);
    }
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
 * INTEGER          IVISSE  -->   gradient transposed
 * DOUBLE PRECISION RELAXP  -->   pressure relaxation
 * INTEGER          IPUCOU  -->   velocity pressure coupling
 * DOUBLE PRECISION EXTRAG  -->   wall pressure extrapolation
 * INTEGER          IMRGRA  -->   gradient reconstruction
 * INTEGER          IMGRPR  -->   multigrid algorithm for pressure
 * INTEGER          NTERUP  -->   piso sweep number
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
             int *const itempk,
             int *const itherm,
             int *const itpscl)

{
  int choice;
  char *material = NULL;
  char *phas = NULL;

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
  if (vars->model != NULL && !cs_gui_strcmp(vars->model, "thermal_scalar"))
    cs_gui_reference_initialization("temperature", t0);
  if (cs_gui_strcmp(vars->model, "compressible_model"))
    cs_gui_reference_initialization("mass_molar", xmasmr);

  if (cs_gui_strcmp(vars->model, "thermal_scalar")) {
    material = _thermal_table_choice("material");
    if (material != NULL) {
      if (!(cs_gui_strcmp(material, "user_material"))) {
        cs_gui_reference_initialization("temperature", t0);
        phas = _thermal_table_choice("phas");

        if (!phas) {
          BFT_MALLOC(phas, 6, char);
          strcpy(phas, "undef");
        }

        cs_phys_prop_thermo_plane_type_t thermal_plane = CS_PHYS_PROP_PLANE_PH;
        if (*itherm == 1)
          thermal_plane = CS_PHYS_PROP_PLANE_PT;
        //else if (*itherm == 3)
        //  // TODO compressible
        //  thermal_plane = CS_PHYS_PROP_PLANE_PS;

        cs_thermal_table_set(material,
                             _thermal_table_choice("method"),
                             phas,
                             _thermal_table_choice("reference"),
                             thermal_plane,
                             *itpscl);
      }
      BFT_FREE(material);
    }
  }

  /* ro0, viscl0, cp0, isls0[*iscalt-1] si tables*/
  if (_thermal_table_needed("density") == 0)
    cs_gui_properties_value("density", ro0);
  else
    cs_phys_prop_compute(CS_PHYS_PROP_DENSITY,
                         1, p0, t0, ro0);

  if (_thermal_table_needed("molecular_viscosity") == 0)
    cs_gui_properties_value("molecular_viscosity", viscl0);
  else
    cs_phys_prop_compute(CS_PHYS_PROP_DYNAMIC_VISCOSITY,
                         1, p0, t0, viscl0);

  if (_thermal_table_needed("specific_heat") == 0)
    cs_gui_properties_value("specific_heat", cp0);
  else
    cs_phys_prop_compute(CS_PHYS_PROP_ISOBARIC_HEAT_CAPACITY,
                         1, p0, t0, cp0);

  if (cs_gui_strcmp(vars->model, "compressible_model")) {
    cs_gui_properties_value("volume_viscosity", viscv0);
    cs_gui_properties_value("thermal_conductivity", &visls0[*itempk -1]);
  }

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
 * subroutine cssca2 (iscalt, iscavr, scamin, scamax)
 * *****************
 *
 * integer          iscavr   <--  number of the related variance if any
 *----------------------------------------------------------------------------*/

void CS_PROCF (cssca2, CSSCA2) (const int  *iscavr)
{
#if _XML_DEBUG_
  bft_printf("==>CSSCA2\n");
#endif

  cs_var_t  *vars = cs_glob_var;

  const int kscmin = cs_field_key_id("min_scalar_clipping");
  const int kscmax = cs_field_key_id("max_scalar_clipping");

  /* Specific physics: the min max of the model scalar are not given */
  const int keysca = cs_field_key_id("scalar_id");

  for (int f_id = 0; f_id < cs_field_n_fields(); f_id++) {
    cs_field_t  *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_USER) {
      int i = cs_field_get_key_int(f, keysca) - 1;

      if (iscavr[i] <= 0 ) {
        double scal_min = cs_field_get_key_double(f, kscmin);
        double scal_max = cs_field_get_key_double(f, kscmax);
        cs_gui_variable_value(f->name, "min_value", &scal_min);
        cs_gui_variable_value(f->name, "max_value", &scal_max);
        cs_field_set_key_double(f, kscmin, scal_min);
        cs_field_set_key_double(f, kscmax, scal_max);
#if _XML_DEBUG_
        bft_printf("--min_scalar_clipping[%i] = %f\n", i, scal_min);
        bft_printf("--max_scalar_clipping[%i] = %f\n", i, scal_max);
#endif
      }
    }
  }

  if (cs_gui_strcmp(vars->model, "thermal_scalar")) {
    /* thermal model with no specific physics */
    char *buff = NULL;
    int mdl = gui_thermal_model();

    if (mdl < 20) {
      BFT_MALLOC(buff, 12, char);
      strcpy(buff, "temperature");
    }
    else if (mdl < 30) {
      BFT_MALLOC(buff, 9, char);
      strcpy(buff, "enthalpy");
    }
    else {
      BFT_MALLOC(buff, 13, char);
      strcpy(buff, "total_energy");
    }

    cs_field_t *f = cs_field_by_name(buff);
    BFT_FREE(buff);
    double scal_min = cs_field_get_key_double(f, kscmin);
    double scal_max = cs_field_get_key_double(f, kscmax);
    cs_gui_variable_value(f->name, "min_value", &scal_min);
    cs_gui_variable_value(f->name, "max_value", &scal_max);
    cs_field_set_key_double(f, kscmin, scal_min);
    cs_field_set_key_double(f, kscmax, scal_max);
#if _XML_DEBUG_
    bft_printf("--min_scalar_clipping[%i] = %f\n", i, scal_min);
    bft_printf("--max_scalar_clipping[%i] = %f\n", i, scal_max);
#endif
  }
}


/*----------------------------------------------------------------------------
 * Read reference dynamic and user scalar viscosity
 *----------------------------------------------------------------------------*/

void CS_PROCF (cssca3, CSSCA3) (const    int *const itherm,
                                const    int *const iscalt,
                                const    int *const iscavr,
                                      double *const visls0,
                                      double *const t0,
                                      double *const p0,
                                      double *const cp0)
{
  double result, coeff, density;

  cs_var_t  *vars = cs_glob_var;

  const int keysca = cs_field_key_id("scalar_id");

  if (vars->model != NULL) {

    if (gui_thermal_model()) {
      int i = *iscalt-1;

      if (_thermal_table_needed("thermal_conductivity") == 0)
        cs_gui_properties_value("thermal_conductivity", &visls0[i]);
      else
        cs_phys_prop_compute(CS_PHYS_PROP_THERMAL_CONDUCTIVITY,
                             1, p0, t0, &visls0[i]);

      /* for the Temperature, the diffusivity factor is not divided by Cp */
      if (*itherm != 1)
        visls0[i] = visls0[i] / *cp0;
    }
  }

  /* User scalar
     In the interface, the user gives the diffusion coefficient, whereas in
     the solver, one sets the diffusivity, thus one need to multiply
     this coefficient by the density to remain coherent */

  int n_fields = cs_field_n_fields();
  for (int f_id = 0; f_id < n_fields; f_id++) {
    const cs_field_t  *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_USER) {
      int i = cs_field_get_key_int(f, keysca) - 1;
      if (iscavr[i] <= 0) {

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
#if _XML_DEBUG_
      bft_printf("--visls0[%i] = %f\n", i, visls0[i]);
#endif
    }
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
 * INTEGER          UREF   -->   reference velocity
 * INTEGER          ALMAX  -->   reference length
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

void CS_PROCF (uiprop, UIPROP) (const int *const ivisls,
                                const int *const ismago,
                                const int *const iale,
                                const int *const icp,
                                const int *const iscavr)
{
  int itype = 0;
  int nbp = 5;

  const int keysca = cs_field_key_id("scalar_id");

  if (*ismago>0) nbp++;

  if (*icp>0) nbp++;

  int n_fields = cs_field_n_fields();
  for (int f_id = 0; f_id < n_fields; f_id++) {
    const cs_field_t  *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_USER) {
      int i = cs_field_get_key_int(f, keysca) - 1;
      if (ivisls[i] > 0 && iscavr[i] <= 0)
        nbp++;
    }
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

#if _XML_DEBUG_
  bft_printf("==>UIPROP %i\n",*iappel);
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
  char *name = NULL;

  int ntimaver = cs_gui_get_tag_number("/analysis_control/time_averages/time_average", 1);

  /* for each average */
  for (int i = 0; i < ntimaver; i++) {

    imom = i + 1;

    _get_time_average_data(imom, "time_step_start", &ntdmom[i]);

    /* test on isuite */
    cs_gui_restart_parameters_status("restart", &isuite);

    if (isuite != 0) {
      _get_time_average_data(imom, "restart_from_time_average", &imoold[i]);
      if (imoold[i] == imom) imoold[i] = -2;
    }

    for (int n = 0; n < _get_time_average_n_variables(imom); n++) {

      name = _get_time_average_variable_name(imom, n + 1);
      int idim = _get_time_average_component(imom, n + 1);

      cs_field_t *f = cs_field_by_name_try(name);
      idfmom[((imom-1)*(*ndgmox) + n)*2 + 0] = f->id;
      idfmom[((imom-1)*(*ndgmox) + n)*2 + 1] = idim;
      BFT_FREE(name);
    }
  }
#if _XML_DEBUG_
  bft_printf("==>UIMOYT\n");
  for (i = 0; i < ntimaver; i++) {
    bft_printf("-->ntdmom =  %i\n", ntdmom[i]);
  }
#endif
}

/*-----------------------------------------------------------------------------
 * Return the list of cells describing a given zone.
 *
 * parameters:
 *   zone_id   <--  volume zone id
 *   ncelet    <--  number of cells with halo
 *   faces     -->  number of selected cells
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
 * double precision vel      <--  fluid velocity
 * double precision tsexp    -->  explicit source terms
 * double precision tsimp    -->  implicit source terms
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
  char *label  = NULL;

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
        label = cs_gui_variable_label("velocity");
        BFT_MALLOC(labelU, strlen(label) + 6, char);
        strcpy(labelU, label);
        strcat(labelU, "[0]");
        mei_tree_insert(ev_formula, labelU, 0.0);
        BFT_MALLOC(labelV, strlen(label) + 6, char);
        strcpy(labelV, label);
        strcat(labelV, "[1]");
        mei_tree_insert(ev_formula, labelV, 0.0);
        BFT_MALLOC(labelW, strlen(label) + 6, char);
        strcpy(labelW, label);
        strcat(labelW, "[2]");
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
        BFT_FREE(label);
        BFT_FREE(labelU);
        BFT_FREE(labelV);
        BFT_FREE(labelW);
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
 * subroutine uitssc (f_id, pvar, tsexp, tsimp)
 * *****************
 *
 * integer          f_id     <--  field id
 * double precision pvar     <--  scalar
 * double precision tsexp    -->  explicit source terms
 * double precision tsimp    -->  implicit source terms
 *----------------------------------------------------------------------------*/

void CS_PROCF(uitssc, UITSSC)(const int                  *f_id,
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
      cs_xpath_add_test_attribute(&path, "label", _scalar_label(*f_id));
      cs_xpath_add_test_attribute(&path, "zone_id", zone_id);
      cs_xpath_add_function_text(&path);
      formula = cs_gui_get_text_value(path);
      BFT_FREE(path);

      if (formula != NULL) {
        ev_formula = mei_tree_new(formula);
        mei_tree_insert(ev_formula,"x",0.);
        mei_tree_insert(ev_formula,"y",0.);
        mei_tree_insert(ev_formula,"z",0.);
        mei_tree_insert(ev_formula, _scalar_label(*f_id), 0.0);
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
          mei_tree_insert(ev_formula, _scalar_label(*f_id), pvar[iel]);
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
 * subroutine uitsth (f_id, pvar, tsexp, tsimp)
 * *****************
 *
 * integer          f_id     <--  field id
 * double precision pvar     <--  scalar
 * double precision tsexp    -->  explicit source terms
 * double precision tsimp    -->  implicit source terms
 *----------------------------------------------------------------------------*/

void CS_PROCF(uitsth, UITSTH)(const int                  *f_id,
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
      cs_xpath_add_test_attribute(&path, "label", _scalar_label(*f_id));
      cs_xpath_add_test_attribute(&path, "zone_id", zone_id);
      cs_xpath_add_function_text(&path);
      formula = cs_gui_get_text_value(path);
      BFT_FREE(path);

      if (formula != NULL) {
        ev_formula = mei_tree_new(formula);
        mei_tree_insert(ev_formula,"x",0.);
        mei_tree_insert(ev_formula,"y",0.);
        mei_tree_insert(ev_formula,"z",0.);
        mei_tree_insert(ev_formula, _scalar_label(*f_id), 0.0);
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
          mei_tree_insert(ev_formula, _scalar_label(*f_id), pvar[iel]);
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
 *   formula        <--  mei formula
 *   symbols        <--  array of symbol to check
 *   symbol_size    <--  number of symbol in symbols
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
 * subroutine uiiniv
 * *****************
 *
 * integer          ncelet   <--  number of cells with halo
 * integer          isuite   <--  restart indicator
 * integer          iccfth   <--  type of initialisation(compressible model)
 * double precision ro0      <--  value of density if IROVAR=0
 * double precision cp0      <--  value of specific heat if ICP=0
 * double precision viscl0   <--  value of viscosity if IVIVAR=0
 * double precision uref     <--  value of reference velocity
 * double precision almax    <--  value of reference length
 * double precision xyzcen   <--  cell's gravity center
 *----------------------------------------------------------------------------*/

void CS_PROCF(uiiniv, UIINIV)(const int          *ncelet,
                              const int          *isuite,
                                    int          *iccfth,
                              const cs_real_t    *ro0,
                              const cs_real_t    *cp0,
                              const cs_real_t    *viscl0,
                              const cs_real_t    *uref,
                              const cs_real_t    *almax,
                              const double *const xyzcen)
{
  /* Coal combustion: the initialization of the model scalar are not given */

  int icel, iel;
  int zones            = 0;
  int cells            = 0;
  int ccfth            = 0;
  int *cells_list      = NULL;
  char *path           = NULL;
  char *path1          = NULL;
  char *status         = NULL;
  char *zone_id        = NULL;

  cs_var_t  *vars = cs_glob_var;

  /* number of volumic zone */

  zones
    = cs_gui_get_tag_number("/solution_domain/volumic_conditions/zone", 1);

#if _XML_DEBUG_
  bft_printf("==>UIINIV\n");
#endif

  for (int i = 1; i < zones + 1; i++) {

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
        char *path_velocity = cs_xpath_init_path();
        cs_xpath_add_elements(&path_velocity, 4,
            "thermophysical_models",
            "velocity_pressure",
            "initialization",
            "formula");
        cs_xpath_add_test_attribute(&path_velocity, "zone_id", zone_id);
        cs_xpath_add_function_text(&path_velocity);
        char *formula_uvw = cs_gui_get_text_value(path_velocity);

        cs_field_t *c_vel = cs_field_by_name("velocity");

        if (formula_uvw != NULL) {
          mei_tree_t *ev_formula_uvw = mei_tree_new(formula_uvw);
          mei_tree_insert(ev_formula_uvw, "x", 0.);
          mei_tree_insert(ev_formula_uvw, "y", 0.);
          mei_tree_insert(ev_formula_uvw, "z", 0.);
          /* try to build the interpreter */
          if (mei_tree_builder(ev_formula_uvw))
            bft_error(__FILE__, __LINE__, 0,
                      _("Error: can not interpret expression: %s\n %i"),
                      ev_formula_uvw->string, mei_tree_builder(ev_formula_uvw));

          const char *symbols_uvw[] = {"velocity[0]", "velocity[1]", "velocity[2]"};
          if (mei_tree_find_symbols(ev_formula_uvw, 3, symbols_uvw))
            bft_error(__FILE__, __LINE__, 0,
                      _("Error: can not find the required symbol: %s\n"),
                      "velocity[0], velocity[1] ou velocity[2]");

          for (icel = 0; icel < cells; icel++) {
            iel = cells_list[icel] - 1;
            mei_tree_insert(ev_formula_uvw, "x", xyzcen[3 * iel + 0]);
            mei_tree_insert(ev_formula_uvw, "y", xyzcen[3 * iel + 1]);
            mei_tree_insert(ev_formula_uvw, "z", xyzcen[3 * iel + 2]);
            mei_evaluate(ev_formula_uvw);
            if (c_vel->interleaved) {
              c_vel->val[3 * iel    ] = mei_tree_lookup(ev_formula_uvw, "velocity[0]");
              c_vel->val[3 * iel + 1] = mei_tree_lookup(ev_formula_uvw, "velocity[1]");
              c_vel->val[3 * iel + 2] = mei_tree_lookup(ev_formula_uvw, "velocity[2]");
            }
            else {
              c_vel->val[                iel] = mei_tree_lookup(ev_formula_uvw, "velocity[0]");
              c_vel->val[    (*ncelet) + iel] = mei_tree_lookup(ev_formula_uvw, "velocity[1]");
              c_vel->val[2 * (*ncelet) + iel] = mei_tree_lookup(ev_formula_uvw, "velocity[2]");
            }
          }
          mei_tree_destroy(ev_formula_uvw);
        }
        else {
          for (int j=0; j < 3; j++) {
            for (icel = 0; icel < cells; icel++) {
              iel = cells_list[icel]-1;
              if (c_vel->interleaved)
                c_vel->val[3 * iel + j] = 0.0;
              else
                c_vel->val[j * (*ncelet) + iel] = 0.0;
            }
          }
        }
        BFT_FREE(formula_uvw);
        BFT_FREE(path_velocity);

        /* Turbulence variables initialization */
        char *choice = cs_gui_turbulence_initialization_choice(zone_id);

        if (cs_gui_strcmp(choice, "formula")) {
          char *path_turb = cs_xpath_init_path();
          cs_xpath_add_elements(&path_turb, 3,
                                "thermophysical_models",
                                "turbulence",
                                "initialization");
          cs_xpath_add_test_attribute(&path_turb, "zone_id", zone_id);
          cs_xpath_add_element(&path_turb, "formula");
          cs_xpath_add_function_text(&path_turb);
          char *formula_turb = cs_gui_get_text_value(path_turb);
          BFT_FREE(path_turb);

          if (formula_turb != NULL) {
            mei_tree_t *ev_formula_turb = mei_tree_new(formula_turb);
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

            char *model = cs_gui_get_thermophysical_model("turbulence");
            if (model == NULL)
              break;

            if (cs_gui_strcmp(model, "k-epsilon") ||
                cs_gui_strcmp(model, "k-epsilon-PL")) {

              const char *symbols[] = {"k","epsilon"};
              if (mei_tree_find_symbols(ev_formula_turb, 2, symbols))
                bft_error(__FILE__, __LINE__, 0,
                          _("Error: can not find the required symbol: %s\n"),
                          "k or epsilon");

              cs_field_t *c_k   = cs_field_by_name("k");
              cs_field_t *c_eps = cs_field_by_name("epsilon");

              for (icel = 0; icel < cells; icel++) {
                iel = cells_list[icel] - 1;
                mei_tree_insert(ev_formula_turb, "x", xyzcen[3 * iel + 0]);
                mei_tree_insert(ev_formula_turb, "y", xyzcen[3 * iel + 1]);
                mei_tree_insert(ev_formula_turb, "z", xyzcen[3 * iel + 2]);
                mei_evaluate(ev_formula_turb);
                c_k->val[iel]   = mei_tree_lookup(ev_formula_turb, "k");
                c_eps->val[iel] = mei_tree_lookup(ev_formula_turb, "epsilon");
              }
            }

            else if (cs_gui_strcmp(model, "Rij-epsilon") || cs_gui_strcmp(model, "Rij-SSG")) {
              const char *symbols[] = {"r11", "r22", "r33", "r12", "r13", "r23", "epsilon"};
              if (mei_tree_find_symbols(ev_formula_turb, 7, symbols))
                bft_error(__FILE__, __LINE__, 0,
                          _("Error: can not find the required symbol: %s\n"),
                          "r11, r22, r33, r12, r13, r23 or epsilon");

              cs_field_t *c_r11 = cs_field_by_name("r11");
              cs_field_t *c_r22 = cs_field_by_name("r22");
              cs_field_t *c_r33 = cs_field_by_name("r33");
              cs_field_t *c_r12 = cs_field_by_name("r12");
              cs_field_t *c_r13 = cs_field_by_name("r13");
              cs_field_t *c_r23 = cs_field_by_name("r23");
              cs_field_t *c_eps = cs_field_by_name("epsilon");

              for (icel = 0; icel < cells; icel++) {
                iel = cells_list[icel]-1;
                mei_tree_insert(ev_formula_turb, "x", xyzcen[3 * iel + 0]);
                mei_tree_insert(ev_formula_turb, "y", xyzcen[3 * iel + 1]);
                mei_tree_insert(ev_formula_turb, "z", xyzcen[3 * iel + 2]);
                mei_evaluate(ev_formula_turb);
                c_r11->val[iel] = mei_tree_lookup(ev_formula_turb, "r11");
                c_r22->val[iel] = mei_tree_lookup(ev_formula_turb, "r22");
                c_r33->val[iel] = mei_tree_lookup(ev_formula_turb, "r33");
                c_r12->val[iel] = mei_tree_lookup(ev_formula_turb, "r12");
                c_r13->val[iel] = mei_tree_lookup(ev_formula_turb, "r13");
                c_r23->val[iel] = mei_tree_lookup(ev_formula_turb, "r23");
                c_eps->val[iel] = mei_tree_lookup(ev_formula_turb, "epsilon");
              }
            }

            else if (cs_gui_strcmp(model, "Rij-EBRSM")) {
              const char *symbols[] = {"r11", "r22", "r33", "r12", "r13", "r23", "epsilon", "alpha"};
              if (mei_tree_find_symbols(ev_formula_turb, 8, symbols))
                bft_error(__FILE__, __LINE__, 0,
                          _("Error: can not find the required symbol: %s\n"),
                          "r11, r22, r33, r12, r13, r23, epsilon or alpha");

              cs_field_t *c_r11 = cs_field_by_name("r11");
              cs_field_t *c_r22 = cs_field_by_name("r22");
              cs_field_t *c_r33 = cs_field_by_name("r33");
              cs_field_t *c_r12 = cs_field_by_name("r12");
              cs_field_t *c_r13 = cs_field_by_name("r13");
              cs_field_t *c_r23 = cs_field_by_name("r23");
              cs_field_t *c_eps = cs_field_by_name("epsilon");
              cs_field_t *c_alp = cs_field_by_name("alpha");

              for (icel = 0; icel < cells; icel++) {
                iel = cells_list[icel]-1;
                mei_tree_insert(ev_formula_turb, "x", xyzcen[3 * iel + 0]);
                mei_tree_insert(ev_formula_turb, "y", xyzcen[3 * iel + 1]);
                mei_tree_insert(ev_formula_turb, "z", xyzcen[3 * iel + 2]);
                mei_evaluate(ev_formula_turb);
                c_r11->val[iel] = mei_tree_lookup(ev_formula_turb, "r11");
                c_r22->val[iel] = mei_tree_lookup(ev_formula_turb, "r22");
                c_r33->val[iel] = mei_tree_lookup(ev_formula_turb, "r33");
                c_r12->val[iel] = mei_tree_lookup(ev_formula_turb, "r12");
                c_r13->val[iel] = mei_tree_lookup(ev_formula_turb, "r13");
                c_r23->val[iel] = mei_tree_lookup(ev_formula_turb, "r23");
                c_eps->val[iel] = mei_tree_lookup(ev_formula_turb, "epsilon");
                c_alp->val[iel] = mei_tree_lookup(ev_formula_turb, "alpha");
              }
            }

            else if (cs_gui_strcmp(model, "v2f-BL-v2/k")) {
              const char *symbols[] = {"k", "epsilon", "phi", "alpha"};
              if (mei_tree_find_symbols(ev_formula_turb, 4, symbols))
                bft_error(__FILE__, __LINE__, 0, _("Error: can not find the required symbol: %s\n"),
                          "k, epsilon, phi of al");

              cs_field_t *c_k   = cs_field_by_name("k");
              cs_field_t *c_eps = cs_field_by_name("epsilon");
              cs_field_t *c_phi = cs_field_by_name("phi");
              cs_field_t *c_alp = cs_field_by_name("alpha");

              for (icel = 0; icel < cells; icel++) {
                iel = cells_list[icel]-1;
                mei_tree_insert(ev_formula_turb, "x", xyzcen[3 * iel + 0]);
                mei_tree_insert(ev_formula_turb, "y", xyzcen[3 * iel + 1]);
                mei_tree_insert(ev_formula_turb, "z", xyzcen[3 * iel + 2]);
                mei_evaluate(ev_formula_turb);
                c_k->val[iel]   = mei_tree_lookup(ev_formula_turb, "k");
                c_eps->val[iel] = mei_tree_lookup(ev_formula_turb, "epsilon");
                c_phi->val[iel] = mei_tree_lookup(ev_formula_turb, "phi");
                c_alp->val[iel] = mei_tree_lookup(ev_formula_turb, "alpha");
              }
            }

            else if (cs_gui_strcmp(model, "k-omega-SST")) {
              const char *symbols[] = {"k", "omega"};
              if (mei_tree_find_symbols(ev_formula_turb, 2, symbols))
                bft_error(__FILE__, __LINE__, 0,
                          _("Error: can not find the required symbol: %s\n"),
                          "k or omega");

              cs_field_t *c_k   = cs_field_by_name("k");
              cs_field_t *c_ome = cs_field_by_name("omega");

              for (icel = 0; icel < cells; icel++) {
                iel = cells_list[icel]-1;
                mei_tree_insert(ev_formula_turb, "x", xyzcen[3 * iel + 0]);
                mei_tree_insert(ev_formula_turb, "y", xyzcen[3 * iel + 1]);
                mei_tree_insert(ev_formula_turb, "z", xyzcen[3 * iel + 2]);
                mei_evaluate(ev_formula_turb);
                c_k->val[iel]   = mei_tree_lookup(ev_formula_turb, "k");
                c_ome->val[iel] = mei_tree_lookup(ev_formula_turb, "omega");
              }
            }

            else if (cs_gui_strcmp(model, "Spalart-Allmaras")) {
              const char *symbols[] = {"nu_tilda"};
              if (mei_tree_find_symbols(ev_formula_turb, 1, symbols))
                bft_error(__FILE__, __LINE__, 0,
                          _("Error: can not find the required symbol: %s\n"),
                          "nu_tilda");

              cs_field_t *c_nu = cs_field_by_name("nu_tilda");

              for (icel = 0; icel < cells; icel++) {
                iel = cells_list[icel]-1;
                mei_tree_insert(ev_formula_turb, "x", xyzcen[3 * iel + 0]);
                mei_tree_insert(ev_formula_turb, "y", xyzcen[3 * iel + 1]);
                mei_tree_insert(ev_formula_turb, "z", xyzcen[3 * iel + 2]);
                mei_evaluate(ev_formula_turb);
                c_nu->val[iel] = mei_tree_lookup(ev_formula_turb, "nu_tilda");
              }
            }

            else
              bft_error(__FILE__, __LINE__, 0,
                        _("Invalid turbulence model: %s.\n"), model);
            mei_tree_destroy(ev_formula_turb);
            BFT_FREE(model);
          }
          BFT_FREE(formula_turb);
        }
        BFT_FREE(choice);
      }

      /* Thermal scalar initialization */
      if (gui_thermal_model()) {
        char *path_sca       = NULL;
        char *formula_sca    = NULL;
        mei_tree_t *ev_formula_sca   = NULL;

        path_sca = cs_xpath_init_path();
        cs_xpath_add_elements(&path_sca, 3,
                              "thermophysical_models",
                              "thermal_scalar",
                              "variable");
        cs_xpath_add_element(&path_sca, "formula");
        cs_xpath_add_test_attribute(&path_sca, "zone_id", zone_id);
        cs_xpath_add_function_text(&path_sca);
        formula_sca = cs_gui_get_text_value(path_sca);
        BFT_FREE(path_sca);

        char *buff = NULL;
        int mdl = gui_thermal_model();

        if (mdl < 20) {
          BFT_MALLOC(buff, 12, char);
          strcpy(buff, "temperature");
        }
        else if (mdl < 30) {
          BFT_MALLOC(buff, 9, char);
          strcpy(buff, "enthalpy");
        }
        else {
          BFT_MALLOC(buff, 13, char);
          strcpy(buff, "total_energy");
        }

        cs_field_t *c = cs_field_by_name(buff);

        if (formula_sca != NULL) {
          ev_formula_sca = mei_tree_new(formula_sca);
          mei_tree_insert(ev_formula_sca, "x", 0.);
          mei_tree_insert(ev_formula_sca, "y", 0.);
          mei_tree_insert(ev_formula_sca, "z", 0.);
          /* try to build the interpreter */
          if (mei_tree_builder(ev_formula_sca))
            bft_error(__FILE__, __LINE__, 0,
                      _("Error: can not interpret expression: %s\n %i"),
                      ev_formula_sca->string, mei_tree_builder(ev_formula_sca));

          if (mei_tree_find_symbol(ev_formula_sca, buff))
            bft_error(__FILE__, __LINE__, 0,
                      _("Error: can not find the required symbol: %s\n"),
                      buff);

          if (*isuite == 0) {
            for (icel = 0; icel < cells; icel++) {
              iel = cells_list[icel]-1;
              mei_tree_insert(ev_formula_sca, "x", xyzcen[3 * iel + 0]);
              mei_tree_insert(ev_formula_sca, "y", xyzcen[3 * iel + 1]);
              mei_tree_insert(ev_formula_sca, "z", xyzcen[3 * iel + 2]);
              mei_evaluate(ev_formula_sca);
              c->val[iel] = mei_tree_lookup(ev_formula_sca, buff);
            }
          }
          mei_tree_destroy(ev_formula_sca);
        } else {
          if (*isuite == 0) {
            for (icel = 0; icel < cells; icel++) {
              iel = cells_list[icel]-1;
              c->val[iel] = 0.0;
            }
          }
        }
        BFT_FREE(formula_sca);
        BFT_FREE(buff);
      }

      /* User Scalars initialization */
      int n_fields = cs_field_n_fields();
      for (int f_id = 0; f_id < n_fields; f_id++) {
        const cs_field_t  *f = cs_field_by_id(f_id);
        if (f->type & CS_FIELD_USER) {
          char *path_sca       = NULL;
          char *formula_sca    = NULL;
          mei_tree_t *ev_formula_sca   = NULL;

          path_sca = cs_xpath_init_path();
          cs_xpath_add_elements(&path_sca, 2,
                                "additional_scalars",
                                "variable");
          cs_xpath_add_element(&path_sca, "formula");
          cs_xpath_add_test_attribute(&path_sca, "zone_id", zone_id);
          cs_xpath_add_function_text(&path_sca);
          formula_sca = cs_gui_get_text_value(path_sca);
          BFT_FREE(path_sca);

          if (formula_sca != NULL) {
            ev_formula_sca = mei_tree_new(formula_sca);
            mei_tree_insert(ev_formula_sca, "x", 0.);
            mei_tree_insert(ev_formula_sca, "y", 0.);
            mei_tree_insert(ev_formula_sca, "z", 0.);
            /* try to build the interpreter */
            if (mei_tree_builder(ev_formula_sca))
              bft_error(__FILE__, __LINE__, 0,
                        _("Error: can not interpret expression: %s\n %i"),
                        ev_formula_sca->string, mei_tree_builder(ev_formula_sca));

            if (mei_tree_find_symbol(ev_formula_sca, f->name))
                bft_error(__FILE__, __LINE__, 0,
                        _("Error: can not find the required symbol: %s\n"),
                        f->name);

            if (*isuite == 0) {
              for (icel = 0; icel < cells; icel++) {
                iel = cells_list[icel]-1;
                mei_tree_insert(ev_formula_sca, "x", xyzcen[3 * iel + 0]);
                mei_tree_insert(ev_formula_sca, "y", xyzcen[3 * iel + 1]);
                mei_tree_insert(ev_formula_sca, "z", xyzcen[3 * iel + 2]);
                mei_evaluate(ev_formula_sca);
                f->val[iel] = mei_tree_lookup(ev_formula_sca, f->name);
              }
            }
            mei_tree_destroy(ev_formula_sca);
          } else {
            if (*isuite == 0) {
              for (icel = 0; icel < cells; icel++) {
                iel = cells_list[icel]-1;
                f->val[iel] = 0.0;
              }
            }
          }
          BFT_FREE(formula_sca);
        }
      }
      /* Meteo Scalars initialization */
      if (cs_gui_strcmp(vars->model, "atmospheric_flows")) {
        int    size;
        char *name       = NULL;
        char *path_meteo = NULL;
        char *formula_meteo  = NULL;
        mei_tree_t *ev_formula_meteo = NULL;

        path_meteo = cs_xpath_init_path();
        cs_xpath_add_elements(&path_meteo, 3,
                               "thermophysical_models",
                               "atmospheric_flows",
                               "variable");

        cs_gui_get_text_values(path_meteo, &size);
        BFT_FREE(path_meteo);

        for (int j = 0; j < size; j++)
        {
          path_meteo = cs_xpath_init_path();
          cs_xpath_add_elements(&path_meteo, 2,
                                 "thermophysical_models",
                                 "atmospheric_flows");
          cs_xpath_add_element_num(&path_meteo, "variable", j);
          cs_xpath_add_element(&path_meteo, "name");
          cs_xpath_add_function_text(&path);
          name = cs_gui_get_text_value(path_meteo);

          cs_field_t *c = cs_field_by_name_try(name);
          BFT_FREE(path_meteo);

          path_meteo = cs_xpath_init_path();
          cs_xpath_add_elements(&path_meteo, 2,
                                 "thermophysical_models",
                                 "atmospheric_flows");
          cs_xpath_add_element_num(&path_meteo, "variable", j);
          name = cs_gui_get_text_value(path_meteo);
          cs_xpath_add_element(&path_meteo, "formula");
          cs_xpath_add_test_attribute(&path_meteo, "zone_id", zone_id);
          cs_xpath_add_function_text(&path_meteo);
          formula_meteo = cs_gui_get_text_value(path_meteo);
          BFT_FREE(path_meteo);

          if (formula_meteo != NULL) {
            ev_formula_meteo = mei_tree_new(formula_meteo);
            mei_tree_insert(ev_formula_meteo, "x", 0.);
            mei_tree_insert(ev_formula_meteo, "y", 0.);
            mei_tree_insert(ev_formula_meteo, "z", 0.);
            /* try to build the interpreter */
            if (mei_tree_builder(ev_formula_meteo))
              bft_error(__FILE__, __LINE__, 0,
                        _("Error: can not interpret expression: %s\n %i"),
                        ev_formula_meteo->string, mei_tree_builder(ev_formula_meteo));

            if (mei_tree_find_symbol(ev_formula_meteo, name))
              bft_error(__FILE__, __LINE__, 0,
                        _("Error: can not find the required symbol: %s\n"),
                        name);

            if (*isuite == 0) {
              for (icel = 0; icel < cells; icel++) {
                iel = cells_list[icel]-1;
                mei_tree_insert(ev_formula_meteo, "x", xyzcen[3 * iel + 0]);
                mei_tree_insert(ev_formula_meteo, "y", xyzcen[3 * iel + 1]);
                mei_tree_insert(ev_formula_meteo, "z", xyzcen[3 * iel + 2]);
                mei_evaluate(ev_formula_meteo);
                c->val[iel] = mei_tree_lookup(ev_formula_meteo, name);
              }
            }
            mei_tree_destroy(ev_formula_meteo);
          }
          else {
            if (*isuite == 0) {
              for (icel = 0; icel < cells; icel++) {
                iel = cells_list[icel]-1;
                c->val[iel] = 0.0;
              }
            }
          }
          BFT_FREE(formula_meteo);
          BFT_FREE(name);
        }
      }

      if (cs_gui_strcmp(vars->model, "compressible_model")) {
        char *formula        = NULL;
        char *buff           = NULL;
        mei_tree_t *ev_formula       = NULL;
        char *name[] = {"pressure", "temperature", "total_energy", "density"};

        ccfth = 10000;
        for (int j = 0; j < 4; j++) {
          path = cs_xpath_short_path();
          if (j < 3)
            cs_xpath_add_element(&path, "variable");
          else
            cs_xpath_add_element(&path, "property");
          cs_xpath_add_test_attribute(&path, "name", name[j]);
          cs_xpath_add_element(&path, "formula");
          cs_xpath_add_test_attribute(&path, "zone_id", zone_id);
          BFT_MALLOC(path1, strlen(path) +1, char);
          strcpy(path1, path);
          cs_xpath_add_attribute(&path, "status");
          buff = cs_gui_get_attribute_value(path);

          if (cs_gui_strcmp(buff, "on")) {
            if (j == 0)
              ccfth = ccfth * 2;
            else if (j == 1)
              ccfth = ccfth * 5;
            else if (j == 2)
              ccfth = ccfth * 7;
            else if (j == 3)
              ccfth = ccfth * 3;

            cs_field_t *c = cs_field_by_name_try(name[j]);

            cs_xpath_add_function_text(&path1);
            formula = cs_gui_get_text_value(path1);
            ev_formula = _init_mei_tree(formula, name[j]);
            if (*isuite == 0) {
              for (icel = 0; icel < cells; icel++) {
                iel = cells_list[icel]-1;
                mei_tree_insert(ev_formula, "x", xyzcen[3 * iel + 0]);
                mei_tree_insert(ev_formula, "y", xyzcen[3 * iel + 1]);
                mei_tree_insert(ev_formula, "z", xyzcen[3 * iel + 2]);
                mei_evaluate(ev_formula);
                c->val[iel] = mei_tree_lookup(ev_formula, name[j]);
              }
            }
            mei_tree_destroy(ev_formula);
          }
          BFT_FREE(buff);
          BFT_FREE(formula);
          BFT_FREE(path);
          BFT_FREE(path1);
        }
        *iccfth = ccfth;
      }
      BFT_FREE(cells_list);

#if _XML_DEBUG_
      bft_printf("--zone zone_id: %s\n", zone_id);
      bft_printf("--zone's element number: %i\n", cells);

      if (*isuite == 0) {
        double initial_value;
        for (int f_id = 0; f_id < cs_field_n_fields(); f_id++) {
          const cs_field_t *f = cs_field_by_id(f_id);
          cs_gui_variable_initial_value(f->name, zone_id, &initial_value);
          bft_printf("--initial value for %s: %f\n", f->name, initial_value);
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
 *   a_ij     <--  change matrix from the local frame to the global frame
 *   in_ij    <--  head losses matrix in the local frame
 *   out_ij   <--  head losses matrix in the global frame
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
 *   zone_id   <--  id of the volume zone
 *   c         <--  name of the coefficient
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
 * integer          iappel   <--  number of calls during a time step
 * integer          ncelet   <--  number of cells with halo
 * integer          ncepdp  -->   number of cells with head losses
 * integer          icepdc  -->   ncepdp cells number with head losses
 * double precision ckupdc  -->   head losses matrix
 *----------------------------------------------------------------------------*/

void CS_PROCF(uikpdc, UIKPDC)(const int*   iappel,
                              const int*   ncelet,
                                    int    ncepdp[],
                                    int    icepdc[],
                                    double ckupdc[])
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

    cs_field_t *c_vel = cs_field_by_name("velocity");

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
          iel = cells_list[j] - 1;
          if (c_vel->interleaved) {
            vit = c_vel->val_pre[3 * iel    ] * c_vel->val_pre[3 * iel    ] +
                  c_vel->val_pre[3 * iel + 1] * c_vel->val_pre[3 * iel + 1] +
                  c_vel->val_pre[3 * iel + 2] * c_vel->val_pre[3 * iel + 2];
          }
          else {
            vit = c_vel->val_pre[                iel] * c_vel->val_pre[                iel] +
                  c_vel->val_pre[    (*ncelet) + iel] * c_vel->val_pre[    (*ncelet) + iel] +
                  c_vel->val_pre[2 * (*ncelet) + iel] * c_vel->val_pre[2 * (*ncelet) + iel];
          }
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
 * INTEGER          NCEL     <--  number of cells whithout halo
 * INTEGER          NCELET   <--  number of cells whith halo
 * INTEGER          NSCAUS   <--  number of user scalar including thermal scalar
 * INTEGER          IROM     <--  pointer for density rho
 * INTEGER          IVISCL   <--  pointer for mulecular viscosity mu
 * INTEGER          ICP      <--  pointer for specific heat Cp
 * INTEGER          IVISLS   <--  pointer for Lambda/Cp
 * INTEGER          IROVAR   <--  =1 if rho variable, =0 if rho constant
 * INTEGER          IVIVAR   <--  =1 if mu variable, =0 if mu constant
 * INTEGER          ISCA     <--  indirection array for scalar number
 * INTEGER          ISCALT   <--  pointer for the thermal scalar in ISCA
 * INTEGER          ISCAVR   <--  scalars that are variance
 * INTEGER          IPPROC   <--  indirection array for cell properties
 * INTEGER          IVISCV   <--  pointer for volumic viscosity viscv
 * INTEGER          ITEMPK   <--  pointer for temperature (in K)
 * DOUBLE PRECISION P0       <--  pressure reference value
 * DOUBLE PRECISION T0       <--  temperature reference value
 * DOUBLE PRECISION RO0      <--  density reference value
 * DOUBLE PRECISION CP0      <--  specific heat reference value
 * DOUBLE PRECISION VISCL0   <--  dynamic viscosity reference value
 * DOUBLE PRECISION VISLS0   <--  diffusion coefficient of the scalars
 * DOUBLE PRECISION VISCV0   <--  volumic viscosity
 *----------------------------------------------------------------------------*/

void CS_PROCF(uiphyv, UIPHYV)(const cs_int_t  *ncel,
                              const cs_int_t  *ncelet,
                              const cs_int_t  *nscaus,
                              const cs_int_t  *itherm,
                              const cs_int_t  *iviscl,
                              const cs_int_t  *icp,
                              const cs_int_t   ivisls[],
                              const cs_int_t  *irovar,
                              const cs_int_t  *ivivar,
                              const cs_int_t   isca[],
                              const cs_int_t  *iscalt,
                              const cs_int_t   iscavr[],
                              const cs_int_t   ipproc[],
                              const cs_int_t  *iviscv,
                              const cs_int_t  *itempk,
                              const cs_real_t *p0,
                              const cs_real_t *t0,
                              const cs_real_t *ro0,
                              const cs_real_t *cp0,
                              const cs_real_t *viscl0,
                              const cs_real_t *visls0,
                              const cs_real_t *viscv0)
{
  char *path = NULL;
  char *law = NULL;
  double time0;
  mei_tree_t *ev_law = NULL;
  int i, iel;

  cs_var_t  *vars = cs_glob_var;

  /* law for density */
  if (!cs_gui_strcmp(vars->model, "compressible_model") ||
       cs_glob_time_step->nt_cur == 1) {
      if (*irovar == 1) {
          cs_field_t *c_rho = CS_F_(rho);
          _physical_property("density", "density",
                             ncel, ncelet, itherm, iscalt, icp,
                             p0, ro0, cp0, viscl0, visls0,
                             c_rho->val);
      }
  }

  /* law for molecular viscosity */
  if (*ivivar == 1) {
    cs_field_t *c_mu = CS_F_(mu);
    _physical_property("molecular_viscosity", "molecular_viscosity",
                       ncel, ncelet, itherm, iscalt, icp,
                       p0, ro0, cp0, viscl0, visls0,
                       c_mu->val);
  }

  /* law for specific heat */
  if (*icp > 0) {
    cs_field_t *c_cp = CS_F_(cp);
    _physical_property("specific_heat", "specific_heat",
                       ncel, ncelet, itherm, iscalt, icp,
                       p0, ro0, cp0, viscl0, visls0,
                       c_cp->val);
  }

  /* law for thermal conductivity */
  if (*iscalt > 0)
    if (ivisls[*iscalt -1] > 0) {
      cs_field_t  *cond_dif = NULL;
      bool         th_f_use_cp = false;

      cs_field_t *_th_f[] = {CS_F_(t), CS_F_(h), CS_F_(energy)};
      bool _th_f_use_cp[] = {true, false, false};

      for (i = 0; i < 3; i++)
        if (_th_f[i]) {
          if ((_th_f[i])->type & CS_FIELD_VARIABLE) {
            th_f_use_cp = _th_f_use_cp[i];
            int k = cs_field_key_id("scalar_diffusivity_id");
            int cond_diff_id = cs_field_get_key_int(_th_f[i], k);
            if (cond_diff_id > -1)
              cond_dif = cs_field_by_id(cond_diff_id);
            break;
          }
        }
      _physical_property("thermal_conductivity", "thermal_conductivity",
                         ncel, ncelet, itherm, iscalt, icp,
                         p0, ro0, cp0, viscl0, visls0,
                         cond_dif->val);
    }

  /* law for volumic viscosity (compressible model) */
  if (cs_gui_strcmp(vars->model, "compressible_model")) {
    if (*iviscv > 0) {
      cs_field_t *c = cs_field_by_name_try("volume_viscosity");
      _compressible_physical_property("volume_viscosity",
                                      "volume_viscosity", c->id,
                                      ncel,
                                      itempk,
                                      p0, t0, ro0,
                                      visls0, viscv0);
    }
  }

  /* law for scalar diffusivity */
  int user_id = -1;
  int n_fields = cs_field_n_fields();
  for (int f_id = 0; f_id < n_fields; f_id++) {
    const cs_field_t  *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_USER) {
      user_id++;

      int user_law = 0;

      if (iscavr[user_id] <= 0 && ivisls[user_id] > 0) {
        char *prop_choice = _properties_choice("name");
        if (cs_gui_strcmp(prop_choice, "variable"))
          user_law = 1;
        BFT_FREE(prop_choice);
      }

      if (user_law) {
        int k = cs_field_key_id("scalar_diffusivity_id");
        int cond_diff_id = cs_field_get_key_int(f, k);
        cs_field_t  *cond_dif = cs_field_by_id(cond_diff_id);

        /* search the formula for the law */
        path = cs_xpath_init_path();
        cs_xpath_add_element(&path, "additional_scalars");
        cs_xpath_add_element_num(&path, "variable", user_id +1);
        cs_xpath_add_element(&path, "property");
        cs_xpath_add_element(&path, "formula");
        cs_xpath_add_function_text(&path);

        law = cs_gui_get_text_value(path);
        BFT_FREE(path);

        if (law != NULL) {
          /* return an empty interpreter */
          time0 = cs_timer_wtime();

          ev_law = mei_tree_new(law);
          BFT_FREE(law);

          for (int f_id2 = 0; f_id2 < n_fields; f_id2++) {
            const cs_field_t  *f2 = cs_field_by_id(f_id2);
            if (f2->type & CS_FIELD_USER)
              mei_tree_insert(ev_law, f2->name, 0.0);
          }

          /* try to build the interpreter */

          if (mei_tree_builder(ev_law))
            bft_error(__FILE__, __LINE__, 0,
                      _("Error: can not interpret expression: %s\n"), ev_law->string);

          if (mei_tree_find_symbol(ev_law, "diffusivity"))
            bft_error(__FILE__, __LINE__, 0,
                      _("Error: can not find the required symbol: %s\n"),
                      "diffusivity");

          /* for each cell, update the value of the table of symbols for each scalar
             (including the thermal scalar), and evaluate the interpreter */

          if (*irovar == 1) {
            cs_field_t *c_rho = CS_F_(rho);
            for (iel = 0; iel < *ncel; iel++) {
              for (int f_id2 = 0; f_id2 < n_fields; f_id2++) {
                const cs_field_t  *f2 = cs_field_by_id(f_id2);
                if (f2->type & CS_FIELD_USER)
                  mei_tree_insert(ev_law,
                                  f2->name,
                                  f2->val[iel]);
              }

              mei_evaluate(ev_law);
              f->val[iel] = mei_tree_lookup(ev_law, "diffusivity") * c_rho->val[iel];
            }
          }
          else {
            for (iel = 0; iel < *ncel; iel++) {
              for (int f_id2 = 0; f_id2 < n_fields; f_id2++) {
                const cs_field_t  *f2 = cs_field_by_id(f_id2);
                if (f2->type & CS_FIELD_USER)
                  mei_tree_insert(ev_law,
                                  f2->name,
                                  f2->val[iel]);
              }

              mei_evaluate(ev_law);
              f->val[iel] = mei_tree_lookup(ev_law, "diffusivity") * (*ro0);
            }
          }
          mei_tree_destroy(ev_law);

          cs_gui_add_mei_time(cs_timer_wtime() - time0);
        }
      }
    }
  }

#if _XML_DEBUG_
  bft_printf("==>UIPHYV\n");
  user_id = -1;
  for (int f_id = 0; f_id < n_fields; f_id++) {
    const cs_field_t  *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_USER) {
      user_id++;
      if (iscavr[user_id] <= 0 && ivisls[user_id] > 0) {
        path = cs_xpath_init_path();
        cs_xpath_add_element(&path, "additional_scalars");
        cs_xpath_add_element_num(&path, "variable", user_id +1);
        cs_xpath_add_element(&path, "property");
        cs_xpath_add_element(&path, "formula");
        cs_xpath_add_function_text(&path);

        law_Ds = cs_gui_get_text_value(path);
        bft_printf("--law for the coefficient of diffusity of the scalar %s: %s\n",
                   f->name, law_Ds);
        BFT_FREE(path);
        BFT_FREE(law_Ds);
      }
    }
#endif
}

/*----------------------------------------------------------------------------
 * 1D profile postprocessing
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIPROF
 * *****************
 *
 * INTEGER          NCELET   <--  number of cells with halo
 * INTEGER          NCEL     <--  number of cells without halo
 * INTEGER          NTMABS   <--  max iterations numbers
 * INTEGER          NTCABS   <--  current iteration number
 * DOUBLE PRECISION TTCABS   <--  current physical time
 * DOUBLE PRECISION TTMABS   <--  max physical time
 * DOUBLE PRECISION TTPABS   <--  physical time at calculation beginning
 * DOUBLE PRECISION XYZCEN   <--  cell's gravity center
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiprof, UIPROF) (const int    *const ncelet,
                                const int    *const ncel,
                                const int    *const ntmabs,
                                const int    *const ntcabs,
                                const double *const ttcabs,
                                const double *const ttmabs,
                                const double *const ttpabs,
                                const double *const xyzcen)
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
  int i, ii, iii, idim;
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
              idim = _get_profile_component(i, iii);

              cs_field_t *f = cs_field_by_name_try(name);

              if (f != NULL) {
                if (f->type & CS_FIELD_VARIABLE) {
                  if (f->interleaved)
                    array[iii+4] = f->val[3 * iel + idim];
                  else
                    array[iii+4] = f->val[iel + idim * (*ncelet)];
                }
                else
                  array[iii+4] = f->val[iel];
              }
              else {
                char *label = _get_profile_label_name(i, iii);
                const int keylbl = cs_field_key_id("label");
                for (int f_id = 0; f_id < cs_field_n_fields(); f_id++) {
                  f = cs_field_by_id(f_id);
                  char *flab = cs_field_get_key_str(f, keylbl);
                  if (cs_gui_strcmp(label, flab))
                    array[iii+4] = f->val[iel];
                }
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
 * INTEGER          NCHARB  <-- number of coal
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
  BFT_FREE(cs_glob_var->model);
  BFT_FREE(cs_glob_var->model_value);
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
