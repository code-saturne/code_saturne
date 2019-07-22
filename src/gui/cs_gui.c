/*============================================================================
 * Management of the GUI parameters file: main parameters
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
#include "cs_boundary_zone.h"
#include "cs_equation.h"
#include "cs_equation_param.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_file.h"
#include "cs_log.h"
#include "cs_gui_util.h"
#include "cs_gui_variables.h"
#include "cs_gui_boundary_conditions.h"
#include "cs_gui_specific_physics.h"
#include "cs_gui_mobile_mesh.h"
#include "cs_geom.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_location.h"
#include "cs_multigrid.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_partition.h"
#include "cs_physical_model.h"
#include "cs_prototypes.h"
#include "cs_rotation.h"
#include "cs_selector.h"
#include "cs_timer.h"
#include "cs_time_moment.h"
#include "cs_thermal_model.h"
#include "cs_physical_properties.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_sles.h"
#include "cs_sles_it.h"
#include "cs_turbulence_model.h"
#include "cs_wall_functions.h"
#include "cs_physical_constants.h"
#include "cs_stokes_model.h"
#include "cs_balance_by_zone.h"
#include "cs_fan.h"
#include "cs_volume_zone.h"
#include "cs_gwf_physical_properties.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gui.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

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

/*============================================================================
 * Private global variables
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer on the main variable structure */

cs_var_t    *cs_glob_var = NULL;

/*============================================================================
 * Static local variables
 *============================================================================*/

static char _rij_c_names[6][4] = {"r11", "r22", "r33", "r12", "r23", "r13"};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Remove a node inside a tree.
 *
 * The node is removed from the tree but not destroyed, so it may be inserted
 * in another position.
 *
 * \param[in, out]  tn    tree node to remove
 *----------------------------------------------------------------------------*/

static void
_tree_node_remove(cs_tree_node_t  *tn)
{
  if (tn->prev != NULL)
    tn->prev->next = tn->next;
  if (tn->next != NULL)
    tn->next->prev = tn->prev;

  if (tn->parent != NULL) {
    if (tn->parent->children == tn)
      tn->parent->children = tn->next;
  }

  tn->prev = NULL;
  tn->next = NULL;
}

/*-----------------------------------------------------------------------------
 * Modify double numerical parameters.
 *
 * parameters:
 *   param    <-- label of the numerical parameter
 *   value    <-- value of the numerical parameter
 *----------------------------------------------------------------------------*/

static void
_numerical_double_parameters(const char  *param,
                             double      *value)
{
  cs_tree_node_t *tn = NULL;
  char *path0 = NULL;
  BFT_MALLOC(path0, strlen("numerical_parameters/") + strlen(param) + 1, char);
  strcpy(path0, "numerical_parameters/");
  strcat(path0, param);
  tn = cs_tree_get_node(cs_glob_tree, path0);
  BFT_FREE(path0);

  /* value not changed if path not found */
  cs_gui_node_get_real(tn, value);
}

/*----------------------------------------------------------------------------
 * Return the value of the choice attribute for material, method, ...
 *
 * parameters:
 *   name        <--  name of the property
 *----------------------------------------------------------------------------*/

static const char*
_thermal_table_choice(const char *name)
{
  cs_tree_node_t *tn
    = cs_tree_get_node(cs_glob_tree, "physical_properties/fluid_properties");
  tn = cs_tree_node_get_child(tn, name);

  const char *choice = cs_tree_node_get_tag(tn, "choice");

  return choice;
}

/*----------------------------------------------------------------------------
 * Return the value of the reference table for a given method, ...
 *
 * parameters:
 *   name        <--  name of the option
 *----------------------------------------------------------------------------*/

static const char*
_thermal_table_option(const char *name)
{
  cs_tree_node_t *tn
    = cs_tree_get_node(cs_glob_tree,
                       "physical_properties/fluid_properties/method");

  const char *option = cs_tree_node_get_child_value_str(tn, name);

  return option;
}

/*----------------------------------------------------------------------------
 * Return the value of the choice attribute from a property name.
 *
 * parameters:
 *   property_name        <--  name of the property
 *----------------------------------------------------------------------------*/

static const char*
_properties_choice(const char *property_name)
{
  const char *choice = NULL;
  cs_tree_node_t *tn
    = cs_tree_get_node(cs_glob_tree,
                       "physical_properties/fluid_properties/property");
  tn = cs_tree_node_get_sibling_with_tag(tn, "name", property_name);

  choice = cs_tree_node_get_child_value_str(tn, "choice");

  return choice;
}

/*----------------------------------------------------------------------------
 * Return 0 if default value is needed
 *
 * parameters:
 *   name        <--  name of the property
 *----------------------------------------------------------------------------*/

static int
_thermal_table_needed(const char *name)
{
  int choice = 0;

  const char *prop_choice = _properties_choice(name);
  if (cs_gui_strcmp(prop_choice, "thermal_law"))
    choice = 1;

  return choice;
}

/*-----------------------------------------------------------------------------
 * use MEI for physical property
 *----------------------------------------------------------------------------*/

static void
_physical_property(cs_field_t          *c_prop,
                   const cs_zone_t     *z)
{
  int user_law = 0;
  const char *law = NULL;
  const char *prop_choice = _properties_choice(c_prop->name);

  if (cs_gui_strcmp(prop_choice, "user_law"))
    user_law = 1;

  if (user_law) {

    /* search the formula for the law */
    cs_tree_node_t *tn = cs_tree_find_node(cs_glob_tree, "property");
    while (tn != NULL) {
      const char *name = cs_tree_node_get_child_value_str(tn, "name");
      if (cs_gui_strcmp(name, c_prop->name))
        break;
      else
        tn = cs_tree_find_node_next(cs_glob_tree, tn, "property");
    }
    tn = cs_tree_get_node(tn, "formula");
    law = cs_tree_node_get_value_str(tn);

    if (law != NULL) {
      cs_field_t *fmeg[1] = {c_prop};
      cs_meg_volume_function(z, fmeg);
    }

  }
  else if (cs_gui_strcmp(prop_choice, "thermal_law")) {
    cs_phys_prop_type_t property = -1;

    if (cs_gui_strcmp(c_prop->name, "density"))
      property = CS_PHYS_PROP_DENSITY;

    else if (cs_gui_strcmp(c_prop->name, "molecular_viscosity"))
      property = CS_PHYS_PROP_DYNAMIC_VISCOSITY;

    else if (cs_gui_strcmp(c_prop->name, "specific_heat"))
      property = CS_PHYS_PROP_ISOBARIC_HEAT_CAPACITY;

    else if (cs_gui_strcmp(c_prop->name, "thermal_conductivity"))
      property = CS_PHYS_PROP_THERMAL_CONDUCTIVITY;

    else
      bft_error(__FILE__, __LINE__, 0,
                _("Error: can not evaluate property: %s using a thermal law\n"),
                c_prop->name);

    /* For incompressible flows, the thermodynamic pressure is constant over
     * time and is the reference pressure. */

    cs_lnum_t thermodynamic_pressure_stride = 0;
    cs_lnum_t thermal_f_val_stride = 1;
    cs_real_t _p0 = cs_glob_fluid_properties->p0;
    cs_real_t _t0 = cs_glob_fluid_properties->t0;

    const cs_real_t *thermodynamic_pressure = &_p0;
    const cs_real_t *_thermal_f_val = NULL;

    if (CS_F_(t) != NULL) {
      if (CS_F_(t)->type & CS_FIELD_VARIABLE)
        _thermal_f_val = CS_F_(t)->val;
    }
    else if (CS_F_(h) != NULL) {
      if (CS_F_(h)->type & CS_FIELD_VARIABLE)
        _thermal_f_val = CS_F_(h)->val;
    }
    else if (CS_F_(e_tot) != NULL) {
      if (CS_F_(h)->type & CS_FIELD_VARIABLE) {
        _thermal_f_val = CS_F_(e_tot)->val;
        thermodynamic_pressure = CS_F_(p)->val;
        thermodynamic_pressure_stride = 1;
      }
    }
    else {
      thermal_f_val_stride = 0;
      _thermal_f_val = &_t0;
    }

    const cs_lnum_t ncel = z->n_elts;

    cs_phys_prop_compute(property,
                         ncel,
                         thermodynamic_pressure_stride,
                         thermal_f_val_stride,
                         thermodynamic_pressure,
                         _thermal_f_val,
                         c_prop->val);

  }
}

/*-----------------------------------------------------------------------------
 * Return the value of choice for user scalar's property
 *
 * parameters:
 *   f_name     <-- field name
 *   choice     <-> choice for property
 *----------------------------------------------------------------------------*/

static int
_scalar_properties_choice(const char *f_name,
                          int *choice)
{
  const char *buff = NULL;
  int   ichoice;

  cs_tree_node_t *tn;

  for (tn = cs_tree_get_node(cs_glob_tree, "additional_scalars/variable");
      tn != NULL;) {
    if (!cs_gui_strcmp(f_name, cs_tree_node_get_tag(tn, "name")))
      tn = cs_tree_node_get_next_of_name(tn);
    else
      break;
  }

  tn = cs_tree_get_node(tn, "property/choice");
  buff = cs_tree_node_get_value_str(tn);

  if (buff == NULL) {
    ichoice = 0;

  } else {
    ichoice = 1;

    if (cs_gui_strcmp(buff, "variable"))
      *choice = 1;
    else if (cs_gui_strcmp(buff, "constant"))
      *choice = 0;
    else
      bft_error(__FILE__, __LINE__, 0, _("Invalid node in function %s\n"),
                                       __func__);
  }

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
_scalar_diffusion_value(int         num_sca,
                        cs_real_t  *value)
{
  cs_tree_node_t *tn
    = cs_tree_get_node(cs_glob_tree, "additional_scalars/variable");

  for (int i = 1;
       tn != NULL && i < num_sca ;
       i++) {
   tn = cs_tree_node_get_next_of_name(tn);
  }
  tn = cs_tree_get_node(tn, "property/initial_value");

  cs_gui_node_get_real(tn, value);
}

/*-----------------------------------------------------------------------------
 * Find the node associated with a given variable.
 *
 * parameters:
 *   variable_name  <-- name of variable
 *
 * return:
 *   pointer to associated tree node (NULL if not present)
 *----------------------------------------------------------------------------*/

static cs_tree_node_t *
_find_node_variable(const char  *variable_name)
{
  cs_tree_node_t *tn = cs_tree_find_node(cs_glob_tree, "variable");
  while (tn != NULL) {
    const char *name = cs_tree_node_get_child_value_str(tn, "name");
    if (cs_gui_strcmp(name, variable_name))
      break;
    else
      tn = cs_tree_find_node_next(cs_glob_tree, tn, "variable");
  }

  return tn;
}

/*-----------------------------------------------------------------------------
 * Return value of turbulent flux model
 *
 * parameters:
 *   tn_v   <-- tree node associated with variable
 *   value  --> value of turbulent flux model
 *----------------------------------------------------------------------------*/

static void
_variable_turbulent_flux_model(cs_tree_node_t  *tn_v,
                               int             *value)
{
  const char *result
    = cs_tree_node_get_child_value_str(tn_v, "turbulent_flux_model");

  if (cs_gui_strcmp(result, "SGDH"))
    *value = 0;
  else if (cs_gui_strcmp(result, "GGDH"))
    *value = 10;
  else if (cs_gui_strcmp(result, "EB-GGDH"))
    *value = 11;
  else if (cs_gui_strcmp(result, "AFM"))
    *value = 20;
  else if (cs_gui_strcmp(result, "EB-AFM"))
    *value = 21;
  else if (cs_gui_strcmp(result, "DFM"))
    *value = 30;
  else if (cs_gui_strcmp(result, "EB-DFM"))
    *value = 31;
  else
    *value = 0; /* assign default */
}

/*----------------------------------------------------------------------------
 * Get the attribute value from the scheme order
 *
 * parameters:
 *   tn_v    <-- node assocaited with variable
 *   keyword -->  value of attribute node
 *----------------------------------------------------------------------------*/

static void
_order_scheme_value(cs_tree_node_t  *tn_v,
                    int             *keyword)
{
  cs_tree_node_t *tn = cs_tree_node_get_child(tn_v, "order_scheme");
  const char *choice = cs_tree_node_get_child_value_str(tn, "choice");

  if (cs_gui_strcmp(choice, "centered"))
    *keyword = 1;
  else if (cs_gui_strcmp(choice, "solu"))
    *keyword = 0;
}

/*----------------------------------------------------------------------------
 * Get the attribute value from the cs_tree query.
 *
 * parameters:
 *   node          <--
 *   child         <-- child markup
 *   keyword      -->  value of attribute node
 *----------------------------------------------------------------------------*/

static void
_slope_test_value(cs_tree_node_t  *tn,
                  int             *keyword)
{
  int result = -999;
  cs_gui_node_get_child_status_int(tn, "slope_test", &result);

  if (result == 1)
    *keyword = 0;
  if (result == 0)
    *keyword = 1;
}

/*----------------------------------------------------------------------------
 * Return the attribute choice associated to a child markup from a variable.
 *
 * parameters:
 *   tn_v   <-- tree node associated with the choice
 *   child  <--  child markup
 *----------------------------------------------------------------------------*/

static const char *
_variable_choice(cs_tree_node_t   *tn_v,
                 const char       *child)
{
  cs_tree_node_t *tn = cs_tree_get_node(tn_v, child);

  const char *choice = cs_tree_node_get_child_value_str(tn, "choice");

  return choice;
}

/*-----------------------------------------------------------------------------
 * Modify integer numerical parameters.
 *
 * parameters:
 *   param     <--  label of the numerical parameter
 *   keyword   <--  value of the numerical parameter
 *----------------------------------------------------------------------------*/

static void
_numerical_int_parameters(const char  *param,
                          int         *keyword)
{
  int choice = *keyword;
  int result = *keyword;
  cs_tree_node_t *tn = cs_tree_get_node(cs_glob_tree, "numerical_parameters");

  if (cs_gui_strcmp(param, "gradient_reconstruction")) {

    tn = cs_tree_get_node(tn, param);
    tn = cs_tree_get_node(tn, "choice");
    cs_gui_node_get_int(tn, &choice);
    *keyword = choice;

  } else if (cs_gui_strcmp(param,"piso_sweep_number")) {

    tn = cs_tree_get_node(tn, "velocity_pressure_algo");
    tn = cs_tree_get_node(tn, param);
    cs_gui_node_get_int(tn, &result);
    *keyword = result;

  } else {

    tn = cs_tree_get_node(tn, param);
    cs_gui_node_get_status_int(tn, &result);
    *keyword = result;
  }
}

/*-----------------------------------------------------------------------------
 * Modify gravity parameters.
 *
 * parameters:
 *   param               <--  gravity parameter (GX, GY, GZ)
 *   keyword            <<--  new value of the gravity parameter
 *----------------------------------------------------------------------------*/

static void
_gravity_value(const char  *param,
               double      *value)
{
  cs_tree_node_t *tn
    = cs_tree_get_node(cs_glob_tree, "physical_properties/gravity");
  tn = cs_tree_get_node(tn, param);

  cs_gui_node_get_real(tn, value);
}

/*-----------------------------------------------------------------------------
 * Modify coriolis source terms parameters.
 *
 * parameters:
 *   param    <--  coriolis parameter (omegax, omegay, omegaz)
 *   value    -->  new value of the coriolis parameter
 *----------------------------------------------------------------------------*/

static void
_coriolis_value(const char *param,
                double     *value)
{
  cs_tree_node_t *tn
    = cs_tree_get_node(cs_glob_tree, "physical_properties/omega");
  tn = cs_tree_get_node(tn, param);

  cs_gui_node_get_real(tn, value);
}

/*----------------------------------------------------------------------------
 * Get the value of the choice attribute from a property markup.
 * Return 1 if the cs_tree request has succeeded, 0 otherwise.
 *
 * parameters:
 *   property_name <-- name of the property
 *   choice        --> value of the attribute choice
 *----------------------------------------------------------------------------*/

static int
_properties_choice_id(const char  *property_name,
                      int         *choice)
{
  const char *buff = NULL;
  int   iok = 0;

  buff = _properties_choice(property_name);
  *choice = 0; /* default */
  if (buff)
  {
    iok = 1;
    if (cs_gui_strcmp(buff, "user_law")       ||
        cs_gui_strcmp(buff, "predefined_law") ||
        cs_gui_strcmp(buff, "thermal_law") )
      *choice = 1;
    else if (cs_gui_strcmp(buff, "constant"))
      *choice = 0;
  }
  else
    iok = 0;
  return iok;
}

/*----------------------------------------------------------------------------
 * Return the length choice for initialize turbulence
 *----------------------------------------------------------------------------*/

static const char *
_reference_length_initialization_choice(void)
{
  cs_tree_node_t *tn
    = cs_tree_get_node(cs_glob_tree,
                       "thermophysical_models/turbulence/"
                       "reference_length/choice");

  const char *initialization_choice
    = cs_tree_node_get_value_str(tn);

  return initialization_choice;
}

/*----------------------------------------------------------------------------
 * Return the initialization choice of the turbulence variables.
 *
 * parameters:
 *   zone_id        <--  zone number
 *----------------------------------------------------------------------------*/

static const char *
_turbulence_initialization_choice(const char* zone_id)
{
  cs_tree_node_t *tn
    = cs_tree_get_node(cs_glob_tree,
                       "thermophysical_models/turbulence/initialization");
  tn = cs_tree_node_get_sibling_with_tag(tn, "zone_id", zone_id);
  tn = cs_tree_get_node(tn, "choice");

  const char *initialization_choice
    = cs_tree_node_get_value_str(tn);

  return initialization_choice;
}

/*==========================
 * FOR VOLUMIC ZONES
 *==========================*/

/*----------------------------------------------------------------------------
 * Check if a given zone type attribute is active
 *
 * parameters:
 *   z_id  <--  zone id
 *   attr  <--  attribute checked for
 *
 * return:
 *   true if attribute checked for is active, false otherwise
 *----------------------------------------------------------------------------*/

static bool
_zone_id_is_type(int          z_id,
                 const char  *attr)
{
  bool retval = false;

  const char *status = NULL;
  cs_tree_node_t *tn
    = cs_tree_get_node(cs_glob_tree, "solution_domain/volumic_conditions/zone");
  for (int i = 1; tn != NULL && i < z_id; i++) {
    tn = cs_tree_node_get_next_of_name(tn);
  }
  tn = cs_tree_get_node(tn, attr);
  status = cs_tree_node_get_value_str(tn);

  if (status != NULL) {
    if (cs_gui_strcmp(status, "on"))
      retval = true;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Indicate if a given zone type attribute is active
 *
 * parameters:
 *   tn       <-- tree node associated with the zone variable
 *   type_str <-- string describing the type
 *----------------------------------------------------------------------------*/

static bool
_zone_is_type(cs_tree_node_t  *tn,
              const char      *type_str)
{
  bool retval = false;

  const char *type_s = cs_tree_node_get_tag(tn, type_str);
  if (type_s != NULL) {
    if (strcmp(type_s, "on") == 0)
      retval = true;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Add a test attribute to an cs_tree query for a given zone id.
 *
 * This should dissapear in the future, when zone names replace ids
 * in the XML file.
 *
 * Note that ids in the zone model are 0-based, while those in the XML
 * are 1-based. The shift is handled by adding a default zone in
 * the base model.
 *
 * parameters:
 *   tn    <-- parent tree node
 *   z_id  <-- zone id (O-based)
 *
 * return:
 *   sibling node matching attribute
 *----------------------------------------------------------------------------*/

static cs_tree_node_t *
_add_zone_id_test_attribute(cs_tree_node_t *tn,
                            int             z_id)
{
  char z_id_str[32];
  snprintf(z_id_str, 31, "%d", z_id);

  return cs_tree_node_get_sibling_with_tag(tn, "zone_id", z_id_str);
}

/*-----------------------------------------------------------------------------
 * Initialize mei tree and check for symbols existence
 *
 * parameters:
 *   formula        <--  mei formula
 *   symbols        <--  array of symbol to check
 *   symbol_size    <--  number of symbol in symbols
 *----------------------------------------------------------------------------*/

static mei_tree_t *
_init_mei_tree(const char  *formula,
               const char  *symbols)
{
  /* return an empty interpreter */
  mei_tree_t *tree = mei_tree_new(formula);

  /* add commun variables */
  mei_tree_insert(tree, "x",    0.0);
  mei_tree_insert(tree, "y",    0.0);
  mei_tree_insert(tree, "z",    0.0);

  /* add variable from notebook */
  cs_gui_add_notebook_variables(tree);

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
 * Return the component of variables or properties or scalar for 1D profile
 *
 * parameters:
 *   tn <-- tree node associated with profile variable
 *----------------------------------------------------------------------------*/

static int
_get_profile_v_component(cs_tree_node_t  *tn)
{
  int comp_id = -1;

  const char *c_name = cs_tree_node_get_tag(tn, "component");

  if (c_name != NULL) {
    int n = sscanf(c_name, "%d", &comp_id);
    if (n != 1)
      bft_error(__FILE__, __LINE__, 0,
                _("Error converting profile component tag %s to integer."),
                c_name);
  }

  return comp_id;
}

/*----------------------------------------------------------------------------
 * Return associated field for a node referring to a field.
 *
 * A node referring to a field must have a "name" tag (i.e. child with a
 * string value).
 *
 * In the case of NEPTUNE_CFD, an additional "field_id" tag requiring
 * addition of the "_<field_id>" extension to match the field's name
 * may be present, so this is tested also.
 *
 * parameters:
 *   tn   <-- tree node associated with profile variable
 *   name <-- value of node's "name" tag (already determined)
 *
 * return:
 *   pointer to field if match, NULL otherwise
 *----------------------------------------------------------------------------*/

static const cs_field_t *
_tree_node_get_field(cs_tree_node_t  *tn)
{
  const cs_field_t *f = NULL;

  const char *name = cs_gui_node_get_tag(tn, "name");
  const char *id_name = cs_tree_node_get_tag(tn, "field_id");

  /* Handle phase id (field_id tag in xml) for NEPTUNE_CFD */

  if (id_name != NULL) {
    if (strcmp(id_name, "none") != 0) {
      char buffer[128];
      snprintf(buffer, 127, "%s_%s", name, id_name);
      buffer[127] = '\0';
      if (strlen(buffer) >= 127)
        bft_error(__FILE__, __LINE__, 0,
                  "Local buffer too small to assemble field name with:\n"
                  "name: %s\n"
                  "field_id: %s\n", name, id_name);
      f = cs_field_by_name_try(buffer);
    }
  }

  /* Handle segregated Reynolds tensor solver */

  cs_turb_rans_model_t *rans_model = cs_glob_turb_rans_model;
  if (f == NULL && rans_model != NULL) {
    if (rans_model->irijco == 0 && strcmp(name, "rij") == 0) {
      int idim = _get_profile_v_component(tn);
      f = cs_field_by_name_try(_rij_c_names[idim]);
    }
  }

  /* General case */
  if (f == NULL) {
    if (strcmp(name, "local_time_step") == 0)
      /* Fix time step output */
      f = CS_F_(dt);
    else
      f = cs_field_by_name_try(name);
  }

  if (f == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Field with name \"%s\" not found"), name);

  return f;
}

/*----------------------------------------------------------------------------
 * Return the label of variables or properties or scalar for 1D profile
 *
 * parameters:
 *   tn_vp <-- tree node associated with profile variable
 *----------------------------------------------------------------------------*/

static void
_write_profile_v_label_name(cs_tree_node_t  *tn_vp,
                            FILE *file,
                            int opf)
{
  char *label = NULL;

  const cs_field_t *f = _tree_node_get_field(tn_vp);
  int idim = _get_profile_v_component(tn_vp);

  if (f != NULL) {

    const char *f_label = cs_field_get_label(f);
    char buf[16] = "";
    size_t buf_len = 0;

    if(f->dim == 3)
      buf_len = 3;
    else if (f->dim > 3)
      buf_len = 4;

    size_t len = strlen(f_label) + buf_len;
    BFT_MALLOC(label, len+1, char);

    if (f->dim > 1) {
      int s_id = 0;
      int e_id = f->dim;
      if (idim > -1) {
        s_id = idim;
        e_id = s_id+1;
      }
      switch(f->dim) {
      case 3:
        for (int ldim = s_id; ldim < e_id; ldim ++) {
          strncpy(buf, cs_glob_field_comp_name_3[ldim], 15);
          sprintf(label, "%s%s", f_label, buf);
          if (opf == 0)
            fprintf(file, " | %s", label);
          else
            fprintf(file, ", %s", label);
        }
        break;
      case 6:
        for (int ldim = s_id ; ldim < e_id ; ldim ++) {
          strncpy(buf, cs_glob_field_comp_name_6[ldim], 15);
          sprintf(label, "%s%s", f_label, buf);
          if (opf == 0)
            fprintf(file, " | %s", label);
          else
            fprintf(file, ", %s", label);
        }
        break;
      case 9:
        for (int ldim = s_id ; ldim < e_id ; ldim ++) {
          strncpy(buf, cs_glob_field_comp_name_9[ldim], 15);
          sprintf(label, "%s%s", f_label, buf);
          if (opf == 0)
            fprintf(file, " | %s", label);
          else
            fprintf(file, ", %s", label);
        }
        break;
      default:
        snprintf(buf, 15, "[%d]", idim); buf[15] = '\0';
      }
    }
    /*Dimension 1*/
    else {
      sprintf(label, "%s%s", f_label, buf);
      if (opf == 0)
        fprintf(file, " | %s", label);
      else
        fprintf(file, ", %s", label);
    }

    BFT_FREE(label);
  }

}

/*----------------------------------------------------------------------------
 * Return a string value associated with a child node for a profile.
 *
 * If the matching child node is not present, an error is produced
 *
 * parameters:
 *   tn <-- tree node associated with profile
 *
 * return:
 *   pointer to matching child string
 *----------------------------------------------------------------------------*/

static const char *
_get_profile_child_str(cs_tree_node_t  *tn,
                       const char      *child_name)
{
  const char *name = cs_tree_node_get_child_value_str(tn, child_name);

  if (name == NULL) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("Incorrect setup tree definition for the following node:\n"));
    cs_tree_dump(CS_LOG_DEFAULT, 2, tn);
    bft_error(__FILE__, __LINE__, 0,
              _("Missing child node: %s"), child_name);
  }

  return name;
}

/*----------------------------------------------------------------------------
 * Return an array of integers associated with a child node for a profile.
 *
 * If the matching child node is not present, an error is produced
 *
 * parameters:
 *   tn <-- tree node associated with profile
 *
 * return:
 *   pointer to matching child values
 *----------------------------------------------------------------------------*/

static const int *
_get_profile_child_int(cs_tree_node_t  *tn,
                       const char      *child_name)
{
  const int *v = cs_tree_node_get_child_values_int(tn, child_name);

  if (v == NULL) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("Incorrect setup tree definition for the following node:\n"));
    cs_tree_dump(CS_LOG_DEFAULT, 2, tn);
    bft_error(__FILE__, __LINE__, 0,
              _("Missing child node: %s"), child_name);
  }

  return v;
}

/*----------------------------------------------------------------------------
 * Get output format for 1D profile
 *
 * parameters:
 *   tn <-- tree node associated with profile
 *
 * return:
 *   1 for CSV, 0 for DAT
 *----------------------------------------------------------------------------*/

static int
_get_profile_format(cs_tree_node_t  *tn)
{
  int format = 0;

  const char *format_s
    = cs_tree_node_get_tag(cs_tree_node_get_child(tn, "format"),
                           "name");

  if (format_s != NULL) {
    if (cs_gui_strcmp(format_s, "CSV"))
      format = 1;
    else if (cs_gui_strcmp(format_s, "DAT"))
      format = 0;
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Invalid profile format: %s"), format_s);
  }

  return format;
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
  double  tensorP[3][3], tensorA[3][3], tensorB[3][3];
  double  tensorC[3][3], tensorD[3][3];

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

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
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

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
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
 *   tn_hl <-- tree node associated with zone head losses
 *   c     <-- name of the coefficient
 *----------------------------------------------------------------------------*/

static double
_c_head_losses(cs_tree_node_t  *tn_hl,
               const char      *c)
{
  double value  = 0.0;

  const cs_real_t *v_r = cs_tree_node_get_child_values_real(tn_hl, c);
  if (v_r != NULL)
    value = v_r[0];

  return value;
}

/*-----------------------------------------------------------------------------
 * Get turbomachinery model
 *
 * parameters:
 *   model_type  -->  turbomachinery model type
 *   coupled     -->  use coupled variant
 *----------------------------------------------------------------------------*/

static void
_turbomachinery_model(cs_turbomachinery_model_t  *model_type,
                      bool                       *coupled)
{
  *model_type = CS_TURBOMACHINERY_NONE;
  *coupled = false;

  const char *model = NULL;

  if (!cs_gui_file_is_loaded())
    return;

  cs_tree_node_t *tn
    = cs_tree_get_node(cs_glob_tree,
                       "thermophysical_models/turbomachinery/model");
  model = cs_tree_node_get_value_str(tn);

  if (cs_gui_strcmp(model, "off"))
    *model_type = CS_TURBOMACHINERY_NONE;
  else if (cs_gui_strcmp(model, "transient"))
    *model_type = CS_TURBOMACHINERY_TRANSIENT;
  else if (cs_gui_strcmp(model, "frozen"))
    *model_type = CS_TURBOMACHINERY_FROZEN;
  else if (cs_gui_strcmp(model, "transient_coupled")) {
    *model_type = CS_TURBOMACHINERY_TRANSIENT;
    *coupled = true;
  }
  else if (cs_gui_strcmp(model, "frozen_coupled")) {
    *model_type = CS_TURBOMACHINERY_FROZEN;
    *coupled = true;
  }
}

/*----------------------------------------------------------------------------
 * Return the value of the choice attribute for rotor (turbomachinery)
 *
 * parameters:
 *   rotor_id    <--  id of the rotor
 *   name        <--  name of the property
 *----------------------------------------------------------------------------*/

static double
_rotor_option(int          rotor_id,
              const char  *name)
{
  double value = 0.;

  cs_tree_node_t *tn
    = cs_tree_get_node(cs_glob_tree,
                       "thermophysical_models/turbomachinery/rotor");
  for (int i = 1; tn != NULL && i < rotor_id + 1;  i++) {
    tn = cs_tree_node_get_next_of_name(tn);
  }
  tn = cs_tree_node_get_child(tn, "rotation");
  tn = cs_tree_node_get_child(tn, name);

  cs_gui_node_get_real(tn, &value);

  return value;
}

/*-----------------------------------------------------------------------------
 * Return the value to a face joining markup for turbomachinery
 *
 * parameters:
 *   keyword <-- label of the markup
 *   number  <-- joining number
 *----------------------------------------------------------------------------*/

static const char *
_get_rotor_face_joining(const char  *keyword,
                        int          number)
{
  cs_tree_node_t *tn
    = cs_tree_get_node(cs_glob_tree,
                       "thermophysical_models/turbomachinery/"
                       "joining/face_joining");
  for (int i = 1; tn != NULL && i < number; i++) {
    tn = cs_tree_node_get_next_of_name(tn);
  }
  tn = cs_tree_get_node(tn, keyword);

  const char *value = cs_tree_node_get_value_str(tn);

  return value;
}

/*----------------------------------------------------------------------------
 * Return volume zone id in tree (1 to n).
 *
 * Also checks that zone definitions are ordered if id_e > -1
 *
 * Note that the name tag for zones actually defines an integer (1 to n).
 * This tag should be removed in the future to only use the zone label
 * (the actual zone name).
 *
 * parameters:
 *   tn   <-- associated tree node
 *   id_e <-- expected zone id (0 to n-1)
 *
 * return
 *   zone's id in tree (1 to n)
 *----------------------------------------------------------------------------*/

static int
_v_zone_t_id(cs_tree_node_t  *tn,
             int              id_e)
{
  int z_t_id = id_e + 1;

  const char *id_s = cs_tree_node_get_tag(tn, "id");
  if (id_s != NULL) {
    z_t_id = atoi(id_s);
    if (id_e > -1 && z_t_id != id_e + 1)
      bft_printf(_("\n"
                   " Warning: noncontiguous %s zone ids in XML:\n"
                   "          zone with index %d has id %d.\n"),
                 tn->name, id_e, z_t_id);
  }

  return z_t_id;
}

/*----------------------------------------------------------------------------
 * Return volume zone tree node with a given id (name)
 *
 * parameters:
 *   tn_vc  <-- associated volume conditions tree node (parent)
 *   z_t_id <-- zone id in tree (1 to n)
 *
 * return
 *   pointer to associated zone node, or NULL if not found
 *----------------------------------------------------------------------------*/

static cs_tree_node_t *
_v_zone_node_by_id(cs_tree_node_t  *tn_vc,
                   int              z_t_id)
{
  cs_tree_node_t *retval = NULL;

  for (cs_tree_node_t *tn = cs_tree_node_get_child(tn_vc, "zone");
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn)) {
    if (z_t_id == _v_zone_t_id(tn, -1)) {
      retval = tn;
      break;
    }
  }
  return retval;
}

/*----------------------------------------------------------------------------
 * Ensure volume and boundary zones defined through the GUI are ordered
 * (modifying tree if necessary)
 *----------------------------------------------------------------------------*/

static void
_ensure_zones_order(void)
{
  cs_tree_node_t *tn_parent = NULL;

  /* Volume zones */
  /*------------- */

  tn_parent = cs_tree_get_node(cs_glob_tree,
                               "solution_domain/volumic_conditions");

  /* Check if volume zones are defined in increasing id order */

  bool need_reorder = false;

  int z_id_prev = -1;
  int id = 0;

  for (cs_tree_node_t *tn = cs_tree_node_get_child(tn_parent, "zone");
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn), id++) {
    if (_v_zone_t_id(tn, id) < z_id_prev)
      need_reorder = true;
  }

  const int n_v_zones = id;

  if (need_reorder) {

    cs_lnum_t *order = NULL, *z_ids = NULL;

    /* Build ordering array */

    BFT_MALLOC(z_ids, n_v_zones, cs_lnum_t);
    BFT_MALLOC(order, n_v_zones, cs_lnum_t);

    /* Loop on volume condition zones */

    id = 0;
    for (cs_tree_node_t *tn = cs_tree_node_get_child(tn_parent, "zone");
         tn != NULL;
         tn = cs_tree_node_get_next_of_name(tn), id++) {
      z_ids[id] = _v_zone_t_id(tn, id);
    }

    cs_order_lnum_allocated(NULL, z_ids, order, n_v_zones);

    /* Now loop on zones in id order */

    cs_tree_node_t *tn_head = NULL;
    cs_tree_node_t *tn_tail = NULL;

    for (int i = 0; i < n_v_zones; i++) {

      int z_id = z_ids[order[i]];

      cs_tree_node_t *tn = _v_zone_node_by_id(tn_parent, z_id);
      _tree_node_remove(tn);
      if (tn_head == NULL) {
        tn_head = tn;
        tn_tail = tn;
      }
      else {
        tn->prev = tn_tail;
        tn_tail->next = tn;
        tn_tail = tn;
      }

    }

    if (tn_parent->children != NULL)
      tn_parent->children->prev = tn_tail;
    tn_tail->next = tn_parent->children;
    tn_parent->children = tn_head;

    BFT_FREE(order);
    BFT_FREE(z_ids);
  }

  /* Boundary zones */
  /*--------------- */

  /* Loop on boundary condition zones */

  tn_parent = cs_tree_get_node(cs_glob_tree,
                               "boundary_conditions");

  need_reorder = false;
  int z_id_max = 0;

  id = 0;
  for (cs_tree_node_t *tn = cs_tree_node_get_child(tn_parent, "boundary");
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn), id++) {

    /* Zone id in tree; note that the name tag for boundary zones actually
       defines an integer (1 to n). This tag should be removed in the
       future to only use the zone label (the actual zone name). */

    const char *id_s = cs_tree_node_get_tag(tn, "name");
    if (id_s != NULL) {
      int z_t_id = atoi(id_s);
      if (z_t_id != id + 1)
        need_reorder = true;
      z_id_max = CS_MAX(z_id_max, z_t_id);
    }

  }

  const int n_b_zones = id;

  if (need_reorder) {

    cs_lnum_t *order = NULL, *z_ids = NULL;
    cs_tree_node_t **tn_bcs = NULL;

    /* Build ordering array */

    BFT_MALLOC(z_ids, n_b_zones, cs_lnum_t);
    BFT_MALLOC(order, n_b_zones, cs_lnum_t);
    BFT_MALLOC(tn_bcs, n_b_zones, cs_tree_node_t *);

    /* Loop on volume condition zones */

    id = 0;
    for (cs_tree_node_t *tn = cs_tree_node_get_child(tn_parent, "boundary");
         tn != NULL;
         tn = cs_tree_node_get_next_of_name(tn), id++) {
      const char *id_s = cs_tree_node_get_tag(tn, "name");
      if (id_s != NULL) {
        z_ids[id] = atoi(id_s);
      }
      else
        z_ids[id] = z_id_max + 1 + id;
      tn_bcs[id] = tn;
    }

    cs_order_lnum_allocated(NULL, z_ids, order, n_b_zones);

    BFT_FREE(z_ids);

    /* Now loop on zones in id order */

    cs_tree_node_t *tn_head = NULL;
    cs_tree_node_t *tn_tail = NULL;

    for (int i = 0; i < n_b_zones; i++) {

      cs_tree_node_t *tn = tn_bcs[order[i]];
      _tree_node_remove(tn);
      if (tn_head == NULL) {
        tn_head = tn;
        tn_tail = tn;
      }
      else {
        tn->prev = tn_tail;
        tn_tail->next = tn;
        tn_tail = tn;
      }

    }

    if (tn_parent->children != NULL)
      tn_parent->children->prev = tn_tail;
    tn_tail->next = tn_parent->children;
    tn_parent->children = tn_head;

    BFT_FREE(order);
    BFT_FREE(tn_bcs);
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Thermal model.
 *
 * Fortran Interface:
 *
 * subroutine csther ()
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (csther, CSTHER) (void)
{
  cs_thermal_model_t *thermal = cs_get_glob_thermal_model();

  switch(cs_gui_thermal_model()) {
  case 10:
    thermal->itherm = CS_THERMAL_MODEL_TEMPERATURE;
    thermal->itpscl = CS_TEMPERATURE_SCALE_CELSIUS;
    break;
  case 11:
    thermal->itherm = CS_THERMAL_MODEL_TEMPERATURE;
    thermal->itpscl = CS_TEMPERATURE_SCALE_KELVIN;
    break;
  case 12:
    thermal->itherm = CS_THERMAL_MODEL_TEMPERATURE;
    thermal->itpscl = CS_TEMPERATURE_SCALE_CELSIUS;
    break;
  case 13:
    thermal->itherm = CS_THERMAL_MODEL_TEMPERATURE;
    thermal->itpscl = CS_TEMPERATURE_SCALE_CELSIUS;
    break;
  case 20:
    thermal->itherm = CS_THERMAL_MODEL_ENTHALPY;
    thermal->itpscl = CS_TEMPERATURE_SCALE_KELVIN;
    break;
  case 30:
    thermal->itherm = CS_THERMAL_MODEL_TOTAL_ENERGY;
    thermal->itpscl = CS_TEMPERATURE_SCALE_KELVIN;
    break;
  default:
    thermal->itherm = CS_THERMAL_MODEL_NONE;
    thermal->itpscl = CS_TEMPERATURE_SCALE_NONE;
    break;
  }
}

/*----------------------------------------------------------------------------
 * Turbulence model
 *----------------------------------------------------------------------------*/

void cs_gui_turb_model(void)
{
  cs_tree_node_t *tn_t = cs_tree_get_node(cs_glob_tree,
                                          "thermophysical_models/turbulence");

  const char *model = cs_tree_node_get_tag(tn_t, "model");
  if (model == NULL)
    return;

  int iwallf = -1;
  cs_turb_model_t *turb_mdl = cs_get_glob_turb_model();
  cs_turb_rans_model_t *rans_mdl = cs_get_glob_turb_rans_model();

  if (cs_gui_strcmp(model, "off"))
    turb_mdl->iturb = 0;
  else if (cs_gui_strcmp(model, "mixing_length")) {
    turb_mdl->iturb = 10;
    cs_gui_node_get_child_real(tn_t, "mixing_length_scale", &(rans_mdl->xlomlg));
  }
  else if (cs_gui_strcmp(model, "k-epsilon")) {
    turb_mdl->iturb = 20;
    cs_gui_node_get_child_int(tn_t, "wall_function", &iwallf);
    cs_gui_node_get_child_status_int(tn_t, "gravity_terms", &(rans_mdl->igrake));
  }
  else if (cs_gui_strcmp(model, "k-epsilon-PL")) {
    turb_mdl->iturb = 21;
    cs_gui_node_get_child_int(tn_t, "wall_function", &iwallf);
    cs_gui_node_get_child_status_int(tn_t, "gravity_terms", &(rans_mdl->igrake));
  }
  else if (cs_gui_strcmp(model, "Rij-epsilon")) {
    turb_mdl->iturb = 30;
    cs_gui_node_get_child_int(tn_t, "wall_function", &iwallf);
    cs_gui_node_get_child_status_int(tn_t, "gravity_terms", &(rans_mdl->igrari));
  }
  else if (cs_gui_strcmp(model, "Rij-SSG")) {
    turb_mdl->iturb = 31;
    cs_gui_node_get_child_int(tn_t, "wall_function", &iwallf);
    cs_gui_node_get_child_status_int(tn_t, "gravity_terms", &(rans_mdl->igrari));
  }
  else if (cs_gui_strcmp(model, "Rij-EBRSM")) {
    turb_mdl->iturb = 32;
    cs_gui_node_get_child_int(tn_t, "wall_function", &iwallf);
    cs_gui_node_get_child_status_int(tn_t, "gravity_terms", &(rans_mdl->igrari));
  }
  else if (cs_gui_strcmp(model, "LES_Smagorinsky")) {
    turb_mdl->iturb = 40;
  }
  else if (cs_gui_strcmp(model, "LES_dynamique")) {
    turb_mdl->iturb = 41;
  }
  else if (cs_gui_strcmp(model, "LES_WALE")) {
    turb_mdl->iturb = 42;
  }
  else if (cs_gui_strcmp(model, "v2f-phi")) {
    turb_mdl->iturb = 50;
    cs_gui_node_get_child_int(tn_t, "wall_function", &iwallf);
    cs_gui_node_get_child_status_int(tn_t, "gravity_terms", &(rans_mdl->igrake));
  }
  else if (cs_gui_strcmp(model, "v2f-BL-v2/k")) {
    turb_mdl->iturb = 51;
    cs_gui_node_get_child_int(tn_t, "wall_function", &iwallf);
    cs_gui_node_get_child_status_int(tn_t, "gravity_terms", &(rans_mdl->igrake));
  }
  else if (cs_gui_strcmp(model, "k-omega-SST")) {
    turb_mdl->iturb = 60;
    cs_gui_node_get_child_int(tn_t, "wall_function", &iwallf);
    cs_gui_node_get_child_status_int(tn_t, "gravity_terms", &(rans_mdl->igrake));
  }
  else if (cs_gui_strcmp(model, "Spalart-Allmaras")) {
    turb_mdl->iturb = 70;
  }
  else
    bft_error(__FILE__, __LINE__, 0,
        _("Invalid turbulence model: %s.\n"), model);

  if (iwallf != -1) {
    cs_wall_functions_t *wall_fnt = cs_get_glob_wall_functions();
    wall_fnt->iwallf = (cs_wall_f_type_t)iwallf;
  }

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--model: %s\n", model);
  bft_printf("--iturb = %i\n", turb_mdl->iturb);
  bft_printf("--igrake = %i\n", rans_mdl->igrake);
  bft_printf("--igrari = %i\n", rans_mdl->igrari);
  bft_printf("--iwallf = %i\n", wall_fnt->iwallf);
  bft_printf("--xlomlg = %f\n", rans_mdl->xlomlg);
#endif
}

/*----------------------------------------------------------------------------
 * Define reference length and reference velocity for initialization of the
 * turbulence variables
 *----------------------------------------------------------------------------*/

void cs_gui_turb_ref_values(void)
{
  cs_tree_node_t *tn_t = cs_tree_get_node(cs_glob_tree,
                                          "thermophysical_models/turbulence");

  cs_turb_model_t *turb_mdl = cs_get_glob_turb_model();

  if (turb_mdl->iturb != 0) {
    const char* length_choice = NULL;
    cs_turb_ref_values_t *ref_values = cs_get_glob_turb_ref_values();

    ref_values->uref = 1.; /* default if not specified */

    cs_gui_node_get_child_real(tn_t,
                               "reference_velocity",
                               &(ref_values->uref));

    length_choice = _reference_length_initialization_choice();

    if (length_choice != NULL) {
      if (cs_gui_strcmp(length_choice, "prescribed"))
        cs_gui_node_get_child_real(tn_t,
                                   "reference_length",
                                   &(ref_values->almax));
    }
  }

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--almax = %f\n", ref_values->almax);
  bft_printf("--uref  = %f\n", ref_values->uref);
#endif
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

void CS_PROCF (cscpva, CSCPVA) (void)
{
  int choice;
  cs_fluid_properties_t *phys_pp = cs_get_glob_fluid_properties();

  if (_properties_choice_id("specific_heat", &choice))
    phys_pp->icp = (choice > 0) ? 0 : -1;

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--icp = %i\n", phys_pp->icp);
#endif
}

/*----------------------------------------------------------------------------
 * Volumic viscosity variable or constant indicator.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSCVVVA (IVISCV)
 * *****************
 *
 * INTEGER        IVISCV  --> volumic viscosity variable or constant indicator
 *----------------------------------------------------------------------------*/

void CS_PROCF (csvvva, CSVVVA) (int *iviscv)
{
  int choice;

  if (_properties_choice_id("volume_viscosity", &choice))
    *iviscv = (choice > 0) ? 0 : -1;

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--iviscv = %i\n", *iviscv);
#endif
}

/*----------------------------------------------------------------------------
 * User thermal scalar.
 *
 * Fortran Interface:
 *
 * SUBROUTINE UITHSC
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (uithsc, UITHSC) (void)
{
  cs_var_t  *vars = cs_glob_var;

  BFT_REALLOC(vars->model, strlen("thermal_scalar")+1, char);
  strcpy(vars->model, "thermal_scalar");
}

/*----------------------------------------------------------------------------
 * Constant or variable indicator for the user scalar molecular diffusivity
 *
 * Fortran Interface:
 *
 * subroutine csivis
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (csivis, CSIVIS) (void)
{
  int choice1, choice2;
  int test1, test2;

  cs_var_t  *vars = cs_glob_var;

  const int keysca = cs_field_key_id("scalar_id");
  const int kivisl = cs_field_key_id("diffusivity_id");
  const int kscavr = cs_field_key_id("first_moment_id");
  const int n_fields = cs_field_n_fields();
  const int itherm = cs_glob_thermal_model->itherm;
  const int iscalt = cs_glob_thermal_model->iscalt;

  if (vars->model != NULL && itherm != CS_THERMAL_MODEL_NONE) {
    test1 = _properties_choice_id("thermal_conductivity", &choice1);
    test2 = _properties_choice_id("specific_heat", &choice2);

    if (strcmp(vars->model, "thermal_scalar") == 0 && test1 && test2) {

      for (int f_id = 0; f_id < n_fields; f_id++) {
        cs_field_t  *f = cs_field_by_id(f_id);
        if (f->type & CS_FIELD_VARIABLE) {
          if (cs_field_get_key_int(f, keysca) == iscalt) {
            if (choice1 || choice2)
              cs_field_set_key_int(f, kivisl, 0);
            else
              cs_field_set_key_int(f, kivisl, -1);
          }
        }
      }
    }
  }

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t  *f = cs_field_by_id(f_id);

    if (   (f->type & CS_FIELD_VARIABLE)
        && (f->type & CS_FIELD_USER)) {
      int iscal = cs_field_get_key_int(f, keysca);
      if (iscal > 0) {
        if (cs_field_get_key_int(f, kscavr) < 0) {
          if (_scalar_properties_choice(f->name, &choice1))
            if (iscalt != iscal)
              cs_field_set_key_int(f, kivisl, choice1 - 1);
          // for groundwater we impose variable property
          if (cs_gui_strcmp(vars->model, "groundwater_model"))
            if (iscalt != iscal)
              cs_field_set_key_int(f, kivisl, 0);
        }
      }
    }
  }

  if (cs_gui_strcmp(vars->model, "compressible_model")) {
    int d_f_id = -1;
    const char *prop_choice = _properties_choice("thermal_conductivity");
    if (cs_gui_strcmp(prop_choice, "user_law") ||
        cs_gui_strcmp(prop_choice, "predefined_law"))
      d_f_id = 0;
    cs_field_t *c_temp = cs_field_by_name("temperature");
    cs_field_set_key_int(c_temp, kivisl, d_f_id);
  }
}

/*----------------------------------------------------------------------------
 * Time passing parameter.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSIDTV ()
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (csidtv, CSIDTV) (void)
{
  cs_time_step_options_t *time_opt = cs_get_glob_time_step_options();

  cs_tree_node_t *tn
    = cs_tree_get_node(cs_glob_tree, "analysis_control/time_parameters");
  cs_gui_node_get_child_int(tn, "time_passing", &time_opt->idtvar);

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--idtvar = %i\n", time_opt->idtvar);
#endif
}

/*----------------------------------------------------------------------------
 * Hydrostatic pressure parameter.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSIPHY ()
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (csiphy, CSIPHY) (void)
{
  cs_stokes_model_t *stokes = cs_get_glob_stokes_model();
  int result = stokes->iphydr;
  cs_tree_node_t *tn
    = cs_tree_find_node(cs_glob_tree,
                        "numerical_parameters/hydrostatic_pressure");
  cs_gui_node_get_status_int(tn, &result);
  stokes->iphydr = result;

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--iphydr = %i\n", stokes->iphydr);
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

void CS_PROCF (cscfgp, CSCFGP) (int *icfgrp)
{
  int result = *icfgrp;
  cs_tree_node_t *tn
    = cs_tree_find_node(cs_glob_tree,
                        "numerical_parameters/hydrostatic_equilibrium/");
  cs_gui_node_get_status_int(tn, &result);
  *icfgrp = result;

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
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

void CS_PROCF (csisui, CSISUI) (int *ntsuit,
                                int *ileaux,
                                int *iccvfg)
{
  cs_tree_node_t *tn = cs_tree_get_node(cs_glob_tree,
                                        "calculation_management/start_restart");

  cs_gui_node_get_child_int(tn, "restart_rescue", ntsuit);

  cs_gui_node_get_child_status_int(tn, "restart_with_auxiliary", ileaux);
  cs_gui_node_get_child_status_int(tn, "frozen_field",           iccvfg);

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
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
 * SUBROUTINE CSTIME ()
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (cstime, CSTIME) (void)
{
  /* Default values for time step factor */
  double cdtmin = 0.1, cdtmax = 1000.;

  cs_tree_node_t *tn
    = cs_tree_get_node(cs_glob_tree, "analysis_control/time_parameters");

  cs_time_step_options_t *time_opt = cs_get_glob_time_step_options();
  cs_time_step_t *time_stp = cs_get_glob_time_step();

  cs_gui_node_get_child_real(tn, "time_step_ref", &(time_stp->dt_ref));
  cs_gui_node_get_child_real(tn, "time_step_min_factor", &cdtmin);
  cs_gui_node_get_child_real(tn, "time_step_max_factor", &cdtmax);
  cs_gui_node_get_child_real(tn, "max_courant_num", &(time_opt->coumax));
  cs_gui_node_get_child_real(tn, "max_fourier_num", &(time_opt->foumax));
  cs_gui_node_get_child_real(tn, "time_step_var", &(time_opt->varrdt));
  cs_gui_node_get_child_real(tn, "relaxation_coefficient", &(time_opt->relxst));

  time_opt->dtmin = cdtmin * time_stp->dt_ref;
  time_opt->dtmax = cdtmax * time_stp->dt_ref;

  /* We keep these two lines in case we read an old XML file... */
  cs_gui_node_get_child_real(tn, "time_step_min", &(time_opt->dtmin));
  cs_gui_node_get_child_real(tn, "time_step_max", &(time_opt->dtmax));

  /* Stop criterion */

  cs_real_t  _t_max = -1;

  cs_gui_node_get_child_real(tn, "maximum_time", &_t_max);
  if (_t_max >= 0)
    time_stp->t_max = _t_max;
  else {
    cs_gui_node_get_child_real(tn, "maximum_time_add", &_t_max);
    if (_t_max >= 0)
      time_stp->t_max = time_stp->t_prev + _t_max;
  }

  if (_t_max < 0) {
    int _nt_max = -1;
    cs_gui_node_get_child_int(tn, "iterations", &_nt_max);
    if (_nt_max > -1)
      time_stp->nt_max = _nt_max;
    else {
      cs_gui_node_get_child_int(tn, "iterations_add", &_nt_max);
      if (_nt_max > -1)
        time_stp->nt_max = time_stp->nt_prev + _nt_max;
    }
  }

  cs_gui_node_get_child_status_int(tn, "thermal_time_step", &(time_opt->iptlro));

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--idtvar = %i\n", time_opt->idtvar);
  bft_printf("--iptlro = %i\n", time_opt->iptlro);
  bft_printf("--ntmabs = %i\n", time_stp->nt_max);
  bft_printf("--dtref = %f\n",  time_opt->dtref);
  bft_printf("--dtmin = %f\n",  time_opt->dtmin);
  bft_printf("--dtmax = %f\n",  time_opt->dtmax);
  bft_printf("--coumax = %f\n", time_opt->coumax);
  bft_printf("--foumax = %f\n", time_opt->foumax);
  bft_printf("--varrdt = %f\n", time_opt->varrdt);
  bft_printf("--relxst = %f\n", time_opt->relxst);
#endif
}

/*----------------------------------------------------------------------------
 * Treatment of local numerical aspects:
 *     BLENCV, ISCHCV, ISSTPC, IRCFLU, CDTVAR, NITMAX, EPSILO
 *----------------------------------------------------------------------------*/

void CS_PROCF (uinum1, UINUM1) (double  *cdtvar)
{
  int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  int var_key_id = cs_field_key_id("variable_id");
  cs_var_cal_opt_t var_cal_opt;

  /* 1) variables from velocity_pressure and turbulence */
  /* 1-a) for pressure or hydraulic head */
  cs_field_t *f = NULL;
  if (cs_glob_physical_model_flag[CS_GROUNDWATER] > -1) {
    f = cs_field_by_name("hydraulic_head");
  }
  else {
    f = cs_field_by_name("pressure");
  }
  cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);

  cs_tree_node_t *tn_v = _find_node_variable(f->name);

  cs_gui_node_get_child_real(tn_v, "solver_precision", &var_cal_opt.epsilo);
  cs_gui_node_get_child_int(tn_v, "rhs_reconstruction", &var_cal_opt.nswrsm);
  cs_gui_node_get_child_int(tn_v, "verbosity", &var_cal_opt.iwarni);
  cs_gui_node_get_child_status_int(tn_v, "flux_reconstruction",
                                   &var_cal_opt.ircflu);


  /* For CDO equation */
  cs_equation_param_t *eqp = cs_equation_param_by_name(f->name);
  if (eqp != NULL) {
    eqp->sles_param.eps = var_cal_opt.epsilo;
  }

  /* Set Field calculation options in the field structure */
  cs_field_set_key_struct(f, key_cal_opt_id, &var_cal_opt);

  /* 1-b) for the other variables */
  int n_fields = cs_field_n_fields();
  for (int f_id = 0; f_id < n_fields; f_id++) {
    f = cs_field_by_id(f_id);
    if (   f->type & CS_FIELD_VARIABLE
        && !cs_gui_strcmp(f->name, "pressure")
        && !cs_gui_strcmp(f->name, "hydraulic_head")) {

      int j = cs_field_get_key_int(f, var_key_id) -1;
      cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);

      const char *ref_name = f->name;

      if (   cs_gui_strcmp(f->name, "r11")
          || cs_gui_strcmp(f->name, "r22")
          || cs_gui_strcmp(f->name, "r33")
          || cs_gui_strcmp(f->name, "r12")
          || cs_gui_strcmp(f->name, "r23")
          || cs_gui_strcmp(f->name, "r13"))
        ref_name = "rij";

      tn_v = _find_node_variable(ref_name);

      cs_gui_node_get_child_real(tn_v, "blending_factor", &var_cal_opt.blencv);
      cs_gui_node_get_child_real(tn_v, "solver_precision", &var_cal_opt.epsilo);

      /* For CDO equation, if non-automatic value ie != -1 */
      eqp = cs_equation_param_by_name(f->name);
      if (eqp != NULL && var_cal_opt.epsilo != -1)
        eqp->sles_param.eps = var_cal_opt.epsilo;

      // only for nscaus and model scalar
      cs_gui_node_get_child_real(tn_v, "time_step_factor", &cdtvar[j]);

      _order_scheme_value(tn_v, &var_cal_opt.ischcv);
      _slope_test_value(tn_v, &var_cal_opt.isstpc);
      cs_gui_node_get_child_status_int(tn_v, "flux_reconstruction",
                                       &var_cal_opt.ircflu);

      cs_gui_node_get_child_int(tn_v, "rhs_reconstruction", &var_cal_opt.nswrsm);
      cs_gui_node_get_child_int(tn_v, "verbosity", &var_cal_opt.iwarni);

      // Set Field calculation options in the field structure
      cs_field_set_key_struct(f, key_cal_opt_id, &var_cal_opt);
    }
  }

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  for (int f_id = 0; f_id < n_fields; f_id++) {
    const cs_field_t  *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      j = cs_field_get_key_int(f, var_key_id) -1;
      bft_printf("-->variable[%i] = %s\n", j, f->name);
      bft_printf("--blencv = %f\n", var_cal_opt.blencv);
      bft_printf("--epsilo = %g\n", var_cal_opt.epsilo);
      bft_printf("--cdtvar = %g\n", cdtvar[j]);
      bft_printf("--ischcv = %i\n", var_cal_opt.ischcv);
      bft_printf("--isstpc = %i\n", var_cal_opt.isstpc);
      bft_printf("--ircflu = %i\n", var_cal_opt.ircflu);
      bft_printf("--nswrsm = %i\n", var_cal_opt.nswrsm);
    }
  }
#endif
}

/*----------------------------------------------------------------------------
 * Global numerical parameters.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSNUM2 (RELAXP, EXTRAG, IMRGRA)
 * *****************
 * DOUBLE PRECISION RELAXP  -->   pressure relaxation
 * DOUBLE PRECISION EXTRAG  -->   wall pressure extrapolation
 * INTEGER          IMRGRA  -->   gradient reconstruction
 *----------------------------------------------------------------------------*/

void CS_PROCF (csnum2, CSNUM2)(double *relaxp,
                               double *extrag,
                                  int *imrgra)
{
  cs_piso_t *piso = cs_get_glob_piso();
  cs_stokes_model_t *stokes = cs_get_glob_stokes_model();
  _numerical_int_parameters("gradient_transposed", &(stokes->ivisse));
  _numerical_int_parameters("velocity_pressure_coupling", &(stokes->ipucou));
  _numerical_int_parameters("gradient_reconstruction", imrgra);
  _numerical_int_parameters("piso_sweep_number", &(piso->nterup));
  _numerical_double_parameters("wall_pressure_extrapolation", extrag);
  _numerical_double_parameters("pressure_relaxation", relaxp);

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--ivisse = %i\n", stokes->ivisse);
  bft_printf("--ipucou = %i\n", stokes->ipucou);
  bft_printf("--imrgra = %i\n", *imrgra);
  bft_printf("--nterup = %i\n", piso->nterup);
  bft_printf("--extrag = %f\n", *extrag);
  bft_printf("--relaxp = %f\n", *relaxp);
#endif
}

/*----------------------------------------------------------------------------
 * Treatment of gravity and fluid physical properties
 * Initialize reference pressure and temperature if present
 *----------------------------------------------------------------------------*/

void CS_PROCF (csphys, CSPHYS) (double     *viscv0,
                                double     *visls0,
                                const int  *itempk)
{
  int choice;
  const char *material = NULL;

  cs_var_t  *vars = cs_glob_var;

  const int itherm = cs_glob_thermal_model->itherm;

  cs_physical_constants_t *phys_cst = cs_get_glob_physical_constants();

  _gravity_value("gravity_x", &(phys_cst->gravity[0]));
  _gravity_value("gravity_y", &(phys_cst->gravity[1]));
  _gravity_value("gravity_z", &(phys_cst->gravity[2]));

  cs_real_t w_x, w_y, w_z;
  w_x = 0.;
  w_y = 0.;
  w_z = 0.;

  _coriolis_value("omega_x", &w_x);
  _coriolis_value("omega_y", &w_y);
  _coriolis_value("omega_z", &w_z);

  if (w_x*w_x + w_y*w_y + w_z*w_z > 0.) {
    cs_rotation_define(w_x, w_y, w_z, 0, 0, 0);
    phys_cst->icorio = 1;
  }
  else
    phys_cst->icorio = 0;

  cs_fluid_properties_t *phys_pp = cs_get_glob_fluid_properties();
  cs_gui_fluid_properties_value("reference_pressure", &(phys_pp->p0));

  /* Variable rho and viscl */
  if (_properties_choice_id("density", &choice))
    phys_pp->irovar = choice;

  if (_properties_choice_id("molecular_viscosity", &choice))
    phys_pp->ivivar = choice;
  if (cs_gui_strcmp(vars->model, "compressible_model"))
    if (_properties_choice_id("molecular_viscosity", &choice))
      phys_pp->ivivar = choice;

  /* Read T0 in each case for user */
  cs_gui_fluid_properties_value("reference_temperature", &(phys_pp->t0));

  if (cs_gui_strcmp(vars->model, "compressible_model"))
    cs_gui_fluid_properties_value("reference_molar_mass", &(phys_pp->xmasmr));

  material = _thermal_table_choice("material");
  if (material != NULL) {
    if (!(cs_gui_strcmp(material, "user_material"))) {
      cs_phys_prop_thermo_plane_type_t thermal_plane = CS_PHYS_PROP_PLANE_PH;
      if (itherm <= CS_THERMAL_MODEL_TEMPERATURE)
        thermal_plane = CS_PHYS_PROP_PLANE_PT;
      //else if (itherm == CS_THERMAL_MODEL_TOTAL_ENERGY)
      //  // TODO compressible
      //  thermal_plane = CS_PHYS_PROP_PLANE_PS;

      const int itpscl = cs_glob_thermal_model->itpscl;

      cs_thermal_table_set(material,
                           _thermal_table_choice("method"),
                           _thermal_table_option("reference"),
                           thermal_plane,
                           itpscl);
    }
  }

  /* ro0, viscl0, cp0, isls0[iscalt-1] si tables */
  if (_thermal_table_needed("density") == 0)
    cs_gui_properties_value("density", &phys_pp->ro0);
  else
    cs_phys_prop_compute(CS_PHYS_PROP_DENSITY,
                         1,
                         0,
                         0,
                         &phys_pp->p0,
                         &phys_pp->t0,
                         &phys_pp->ro0);

  if (_thermal_table_needed("molecular_viscosity") == 0)
    cs_gui_properties_value("molecular_viscosity", &phys_pp->viscl0);
  else
    cs_phys_prop_compute(CS_PHYS_PROP_DYNAMIC_VISCOSITY,
                         1,
                         0,
                         0,
                         &phys_pp->p0,
                         &phys_pp->t0,
                         &phys_pp->viscl0);

  if (_thermal_table_needed("specific_heat") == 0)
    cs_gui_properties_value("specific_heat", &phys_pp->cp0);
  else
    cs_phys_prop_compute(CS_PHYS_PROP_ISOBARIC_HEAT_CAPACITY,
                         1,
                         0,
                         0,
                         &phys_pp->p0,
                         &phys_pp->t0,
                         &phys_pp->cp0);

  if (cs_gui_strcmp(vars->model, "compressible_model")) {
    cs_gui_properties_value("volume_viscosity", viscv0);
    cs_gui_properties_value("thermal_conductivity", &visls0[*itempk -1]);
  }

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--gx = %f \n", phys_cst->gravity[0]);
  bft_printf("--gy = %f \n", phys_cst->gravity[1]);
  bft_printf("--gz = %f \n", phys_cst->gravity[2]);
  bft_printf("--icorio = %i \n", cs_glob_physical_constants->icorio);
  bft_printf("--rho = %g , variable %i\n",
             cs_glob_fluid_properties->ro0,
             cs_glob_fluid_properties->irovar);
  bft_printf("--mu = %g , variable %i \n",
             cs_glob_fluid_properties->viscl0,
             cs_glob_fluid_properties->ivivar);
  bft_printf("--Cp = %g \n", cs_glob_fluid_properties->cp0);
  bft_printf("--T0 = %f \n", cs_glob_fluid_properties->t0);
  bft_printf("--P0 = %f \n", cs_glob_fluid_properties->p0);
  if (cs_gui_strcmp(vars->model, "compressible_model")) {
    bft_printf("--viscv0 = %g \n", *viscv0);
    bft_printf("--xmasmr = %f \n", cs_glob_fluid_properties->xmasmr);
  }
#endif
}

/*----------------------------------------------------------------------------
 * User scalar min and max values for clipping.
 *
 * Fortran Interface:
 *
 * subroutine cssca2 (iturt)
 * *****************
 *
 * integer          iturt    -->  turbulent flux model
 *----------------------------------------------------------------------------*/

void CS_PROCF (cssca2, CSSCA2) (int        *iturt)
{
#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
#endif

  cs_var_t  *vars = cs_glob_var;

  const int kscmin = cs_field_key_id("min_scalar_clipping");
  const int kscmax = cs_field_key_id("max_scalar_clipping");

  /* Specific physics: the min max of the model scalar are not given */
  const int keysca = cs_field_key_id("scalar_id");
  const int kscavr = cs_field_key_id("first_moment_id");

  for (int f_id = 0; f_id < cs_field_n_fields(); f_id++) {
    cs_field_t  *f = cs_field_by_id(f_id);
    if (   (f->type & CS_FIELD_VARIABLE)
        && (f->type & CS_FIELD_USER)) {
      int i = cs_field_get_key_int(f, keysca) - 1;
      if (i > -1) {
        if (cs_field_get_key_int(f, kscavr) < 0) {
          double scal_min = cs_field_get_key_double(f, kscmin);
          double scal_max = cs_field_get_key_double(f, kscmax);
          cs_tree_node_t *tn_v = _find_node_variable(f->name);
          cs_gui_node_get_child_real(tn_v, "min_value", &scal_min);
          cs_gui_node_get_child_real(tn_v, "max_value", &scal_max);
          cs_field_set_key_double(f, kscmin, scal_min);
          cs_field_set_key_double(f, kscmax, scal_max);

          if (cs_glob_turb_model->iturb/10 == 3) {
            int turb_mdl;
            _variable_turbulent_flux_model(tn_v, &turb_mdl);
            iturt[i] = turb_mdl;
          }
#if _XML_DEBUG_
          bft_printf("--min_scalar_clipping[%i] = %f\n", i, scal_min);
          bft_printf("--max_scalar_clipping[%i] = %f\n", i, scal_max);
#endif
        }
      }
    }
  }

  if (cs_gui_strcmp(vars->model, "thermal_scalar")) {

    /* thermal model with no specific physics */

    const int itherm = cs_glob_thermal_model->itherm;
    assert(itherm > CS_THERMAL_MODEL_NONE);

    const char *t_names[] = {"temperature", "enthalpy", "total_energy"};

    cs_field_t *f = cs_field_by_name(t_names[itherm-1]);

    double scal_min = cs_field_get_key_double(f, kscmin);
    double scal_max = cs_field_get_key_double(f, kscmax);

    cs_tree_node_t *tn_v = _find_node_variable(f->name);
    cs_gui_node_get_child_real(tn_v, "min_value", &scal_min);
    cs_gui_node_get_child_real(tn_v, "max_value", &scal_max);
    cs_field_set_key_double(f, kscmin, scal_min);
    cs_field_set_key_double(f, kscmax, scal_max);
    int i = cs_field_get_key_int(f, keysca) - 1;

    if (cs_glob_turb_model->iturb/10 == 3) {
      _variable_turbulent_flux_model(tn_v, &(iturt[i]));
    }
#if _XML_DEBUG_
    bft_printf("--min_scalar_clipping[%i] = %f\n", i, scal_min);
    bft_printf("--max_scalar_clipping[%i] = %f\n", i, scal_max);
#endif
  }
}

/*----------------------------------------------------------------------------
 * Read reference dynamic and user scalar viscosity
 *----------------------------------------------------------------------------*/

void CS_PROCF (cssca3, CSSCA3) (double     *visls0)
{
  double result, coeff, density;

  cs_var_t  *vars = cs_glob_var;

  const int keysca = cs_field_key_id("scalar_id");
  const int kscavr = cs_field_key_id("first_moment_id");

  const int itherm = cs_glob_thermal_model->itherm;
  const int iscalt = cs_glob_thermal_model->iscalt;

  if (vars->model != NULL) {

    if (itherm != CS_THERMAL_MODEL_NONE) {
      int i = iscalt-1;

      if (_thermal_table_needed("thermal_conductivity") == 0)
        cs_gui_properties_value("thermal_conductivity", &visls0[i]);
      else
        cs_phys_prop_compute(CS_PHYS_PROP_THERMAL_CONDUCTIVITY,
                             1,
                             0,
                             0,
                             &(cs_glob_fluid_properties->p0),
                             &(cs_glob_fluid_properties->t0),
                             &visls0[i]);

      /* for the Temperature, the diffusivity factor is not divided by Cp */
      if (itherm != CS_THERMAL_MODEL_TEMPERATURE)
        visls0[i] = visls0[i] / cs_glob_fluid_properties->cp0;
    }
  }

  /* User scalar
     In the interface, the user gives the diffusion coefficient, whereas in
     the solver, one sets the diffusivity, thus one need to multiply
     this coefficient by the density to remain coherent */

  if (!cs_gui_strcmp(vars->model, "groundwater_model")) {
    int n_fields = cs_field_n_fields();
    for (int f_id = 0; f_id < n_fields; f_id++) {
      const cs_field_t  *f = cs_field_by_id(f_id);
      if (   (f->type & CS_FIELD_VARIABLE)
          && (f->type & CS_FIELD_USER)) {
        int i = cs_field_get_key_int(f, keysca) - 1;
        if (cs_field_get_key_int(f, kscavr) < 0) {

          if (cs_gui_strcmp(vars->model, "solid_fuels")) {
            /* Air molar mass */
            result = 0.028966;
            cs_gui_fluid_properties_value("reference_molar_mass", &result);
            if (result <= 0)
              bft_error(__FILE__, __LINE__, 0,
                        _("mass molar value is zero or not found in the xml file.\n"));
            density = cs_glob_fluid_properties->p0 *
                      result / (8.31446 *(cs_glob_fluid_properties->t0));
          }
          else
            cs_gui_properties_value("density", &density);

          if (density <= 0)
            bft_error(__FILE__, __LINE__, 0,
                      _("Density value is zero or not found in the xml file.\n"));

          coeff = visls0[i] / density ;
          _scalar_diffusion_value(i+1, &coeff);
          visls0[i] = coeff * density;
        }
#if _XML_DEBUG_
        bft_printf("--visls0[%i] = %f\n", i, visls0[i]);
#endif
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Define porosity.
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIPORO
 *----------------------------------------------------------------------------*/

void CS_PROCF(uiporo, UIPORO)(void)
{
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->cell_cen;

  mei_tree_t *ev_formula  = NULL;

  int n_zones = cs_volume_zone_n_zones();

  /* Porosity fields */
  cs_field_t *fporo = CS_F_(poro);
  cs_field_t *ftporo = CS_F_(t_poro);

  cs_real_t   *porosi = NULL;
  cs_real_6_t *porosf = NULL;

  if (fporo != NULL) {
    porosi = fporo->val;
    if (ftporo != NULL) {
      porosf = (cs_real_6_t *)ftporo->val;
    }
  }

  for (cs_lnum_t iel = 0; iel < n_cells_ext; iel++) {
    porosi[iel] = 1.;
    if (ftporo != NULL) {
      porosf[iel][0] = 1.;
      porosf[iel][1] = 1.;
      porosf[iel][2] = 1.;
      porosf[iel][3] = 0.;
      porosf[iel][4] = 0.;
      porosf[iel][5] = 0.;
    }
  }

  cs_tree_node_t *tn_p
    = cs_tree_get_node(cs_glob_tree, "thermophysical_models/porosities/porosity");

  for (int z_id = 0; z_id < n_zones; z_id++) {
    const cs_zone_t *z = cs_volume_zone_by_id(z_id);

    if (z->type & CS_VOLUME_ZONE_POROSITY) {

      cs_lnum_t  n_cells = z->n_elts;
      const cs_lnum_t *cell_ids = z->elt_ids;

      cs_tree_node_t *tn_zp = _add_zone_id_test_attribute(tn_p, z->id);
      const char *mdl = cs_tree_node_get_child_value_str(tn_zp, "model");
      const char *formula = cs_tree_node_get_child_value_str(tn_zp, "formula");

      if (formula != NULL) {
        if (cs_gui_strcmp(mdl, "anisotropic")) {
          cs_field_t *fmeg[2] = {fporo, ftporo};
          cs_meg_volume_function(z, fmeg);

        } else {
          cs_field_t *fmeg[1] = {fporo};
          cs_meg_volume_function(z, fmeg);
        }

      }
    }
  }
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

void CS_PROCF(uitsnv, UITSNV)(const cs_real_3_t  *restrict vel,
                              cs_real_3_t        *restrict tsexp,
                              cs_real_33_t       *restrict tsimp)
{

  const cs_real_t *restrict cell_f_vol = cs_glob_mesh_quantities->cell_f_vol;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->cell_cen;

  double Su, Sv, Sw;
  double dSudu, dSudv, dSudw;
  double dSvdu, dSvdv, dSvdw;
  double dSwdu, dSwdv, dSwdw;

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
#endif

  cs_field_t *c_rho = CS_F_(rho);

  int n_zones = cs_volume_zone_n_zones();

  cs_tree_node_t *tn_mf
    = cs_tree_get_node(cs_glob_tree,
                       "thermophysical_models/source_terms/momentum_formula");

  for (int z_id = 0; z_id < n_zones; z_id++) {
    const cs_zone_t *z = cs_volume_zone_by_id(z_id);

    if (! (z->type & CS_VOLUME_ZONE_SOURCE_TERM))
      continue;

    if (_zone_id_is_type(z->id, "momentum_source_term")) {
      const cs_lnum_t n_cells = z->n_elts;
      const cs_lnum_t *cell_ids = z->elt_ids;

      cs_tree_node_t *tn = _add_zone_id_test_attribute(tn_mf, z->id);
      const char *formula = cs_tree_node_get_value_str(tn);

      if (formula != NULL) {

        cs_real_t *st_vals = cs_meg_source_terms(z,
                                                 "momentum",
                                                 "momentum_source_term");

        for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
          cs_lnum_t c_id = cell_ids[e_id];

          /* Read values from the newly created array */
          Su = st_vals[12*e_id];
          Sv = st_vals[12*e_id + 1];
          Sw = st_vals[12*e_id + 2];

          dSudu = st_vals[12*e_id + 3];
          dSudv = st_vals[12*e_id + 4];
          dSudw = st_vals[12*e_id + 5];

          dSvdu = st_vals[12*e_id + 6];
          dSvdv = st_vals[12*e_id + 7];
          dSvdw = st_vals[12*e_id + 8];

          dSwdu = st_vals[12*e_id + 9];
          dSwdv = st_vals[12*e_id + 10];
          dSwdw = st_vals[12*e_id + 11];

          /* Fill the explicit and implicit source terms' arrays */
          tsexp[c_id][0] = cell_f_vol[c_id]
                         * ( Su
                           - dSudu * vel[c_id][0]
                           - dSudv * vel[c_id][1]
                           - dSudw * vel[c_id][2] );

          tsexp[c_id][1] = cell_f_vol[c_id]
                         * ( Sv
                           - dSvdu * vel[c_id][0]
                           - dSvdv * vel[c_id][1]
                           - dSvdw * vel[c_id][2] );

          tsexp[c_id][2] = cell_f_vol[c_id]
                         * ( Sw
                           - dSwdu * vel[c_id][0]
                           - dSwdv * vel[c_id][1]
                           - dSwdw * vel[c_id][2] );

          tsimp[c_id][0][0] = cell_f_vol[c_id]*dSudu;
          tsimp[c_id][0][1] = cell_f_vol[c_id]*dSudv;
          tsimp[c_id][0][2] = cell_f_vol[c_id]*dSudw;
          tsimp[c_id][1][0] = cell_f_vol[c_id]*dSvdu;
          tsimp[c_id][1][1] = cell_f_vol[c_id]*dSvdv;
          tsimp[c_id][1][2] = cell_f_vol[c_id]*dSvdw;
          tsimp[c_id][2][0] = cell_f_vol[c_id]*dSwdu;
          tsimp[c_id][2][1] = cell_f_vol[c_id]*dSwdv;
          tsimp[c_id][2][2] = cell_f_vol[c_id]*dSwdw;

        }
        if (st_vals != NULL)
          BFT_FREE(st_vals);
      }
    }
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
 * integer          idarcy   <--  groundwater module activation
 * integer          f_id     <--  field id
 * double precision pvar     <--  scalar
 * double precision tsexp    -->  explicit source terms
 * double precision tsimp    -->  implicit source terms
 *----------------------------------------------------------------------------*/

void CS_PROCF(uitssc, UITSSC)(const int                  *idarcy,
                              const int                  *f_id,
                              const cs_real_t   *restrict pvar,
                              cs_real_t         *restrict tsexp,
                              cs_real_t         *restrict tsimp)
{
  const cs_real_t *restrict cell_f_vol = cs_glob_mesh_quantities->cell_f_vol;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->cell_cen;

  const char *formula = NULL;

  cs_field_t *f = cs_field_by_id(*f_id);

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
#endif

  int n_zones = cs_volume_zone_n_zones();

  for (int z_id = 0; z_id < n_zones; z_id++) {
    const cs_zone_t *z = cs_volume_zone_by_id(z_id);

    if (! (z->type & CS_VOLUME_ZONE_SOURCE_TERM))
      continue;

    /* species source term */
    if (_zone_id_is_type(z->id, "scalar_source_term")) {
      const cs_lnum_t n_cells = z->n_elts;
      const cs_lnum_t *cell_ids = z->elt_ids;

      cs_tree_node_t *tn
        = cs_tree_get_node(cs_glob_tree,
                           "thermophysical_models/source_terms/scalar_formula");
      char z_id_str[32];
      snprintf(z_id_str, 31, "%d", z->id);
      while (tn != NULL){
        const char *name = cs_gui_node_get_tag(tn, "name");
        const char *zone_id = cs_gui_node_get_tag(tn, "zone_id");
        if (cs_gui_strcmp(name, f->name) && cs_gui_strcmp(zone_id, z_id_str))
          break;
        tn = cs_tree_node_get_next_of_name(tn);
      }
      formula = cs_tree_node_get_value_str(tn);

      if (formula != NULL) {
        if (*idarcy == -1) {
          cs_real_t *st_vals = cs_meg_source_terms(z,
                                                   f->name,
                                                   "scalar_source_term");

          for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
            cs_lnum_t c_id = cell_ids[e_id];
            tsimp[c_id] = cell_f_vol[c_id] * st_vals[2 * e_id + 1];
            tsexp[c_id] = cell_f_vol[c_id] * st_vals[2 * e_id]
                        - tsimp[c_id] * pvar[c_id];
          }
          if (st_vals != NULL)
            BFT_FREE(st_vals);
        }
      }
    }
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
  const cs_real_t *restrict cell_f_vol = cs_glob_mesh_quantities->cell_f_vol;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->cell_cen;

  double dS;
  const char *formula = NULL;

  cs_field_t *f = cs_field_by_id(*f_id);

  /* number of volumic zone */

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
#endif

  int n_zones = cs_volume_zone_n_zones();

  for (int z_id = 0; z_id < n_zones; z_id++) {
    const cs_zone_t *z = cs_volume_zone_by_id(z_id);

    if (!(z->type & CS_VOLUME_ZONE_SOURCE_TERM))
      continue;

    /* species source term */
    if (_zone_id_is_type(z->id, "thermal_source_term")) {
      const cs_lnum_t n_cells = z->n_elts;
      const cs_lnum_t *cell_ids = z->elt_ids;

      cs_tree_node_t *tn
        = cs_tree_get_node(cs_glob_tree,
                           "thermophysical_models/source_terms/thermal_formula");
      char z_id_str[32];
      snprintf(z_id_str, 31, "%d", z->id);
      while (tn != NULL){

        const char *name = cs_gui_node_get_tag(tn, "name");
        const char *zone_id = cs_gui_node_get_tag(tn, "zone_id");
        if (cs_gui_strcmp(name, f->name) && cs_gui_strcmp(zone_id, z_id_str))
          break;
        tn = cs_tree_node_get_next_of_name(tn);
      }
      formula = cs_tree_node_get_value_str(tn);

      if (formula != NULL) {

        cs_real_t *st_vals = cs_meg_source_terms(z,
                                                 f->name,
                                                 "thermal_source_term");

        for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
          cs_lnum_t c_id = cell_ids[e_id];

          tsimp[c_id] = cell_f_vol[c_id] * st_vals[2 * e_id + 1];
          tsexp[c_id] = cell_f_vol[c_id] * st_vals[2 * e_id]
                      - tsimp[c_id] * pvar[c_id];
        }
        if (st_vals != NULL)
          BFT_FREE(st_vals);
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Variables and user scalars initialization.
 *
 * Fortran Interface:
 *
 * subroutine uiiniv
 * *****************
 *
 * integer          isuite   <--  restart indicator
 * integer          idarcy   <--  groundwater module activation
 * integer          iccfth   -->  type of initialization (compressible model)
 *----------------------------------------------------------------------------*/

void CS_PROCF(uiiniv, UIINIV)(const int          *isuite,
                              const int          *idarcy,
                              int                *iccfth)
{
  /* Coal combustion: the initialization of the model scalar are not given */
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->cell_cen;

  int ccfth             = 0;

  cs_var_t  *vars = cs_glob_var;

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
#endif

  const int n_zones = cs_volume_zone_n_zones();

  for (int z_id = 0; z_id < n_zones; z_id++) {
    const cs_zone_t *z = cs_volume_zone_by_id(z_id);

    if (z->type & CS_VOLUME_ZONE_INITIALIZATION) {
      const cs_lnum_t n_cells = z->n_elts;
      const cs_lnum_t *cell_ids = z->elt_ids;

      char z_id_str[32];
      snprintf(z_id_str, 31, "%d", z_id);

      if (*isuite == 0) {

        cs_tree_node_t *tn_velocity
          = cs_tree_get_node(cs_glob_tree,
                             "thermophysical_models/velocity_pressure"
                             "/initialization/formula");
        tn_velocity = _add_zone_id_test_attribute(tn_velocity, z->id);
        const char *formula_uvw = cs_tree_node_get_value_str(tn_velocity);

        cs_field_t *c_vel = cs_field_by_name("velocity");

        if (formula_uvw != NULL) {
          cs_real_t *ini_vals = cs_meg_initialization(z, "velocity");
          if (ini_vals != NULL) {
            for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
              cs_lnum_t c_id = cell_ids[e_id];
              for (int d = 0; d < 3; d++)
                c_vel->val[3 * c_id + d] = ini_vals[3 * e_id + d];
            }
            BFT_FREE(ini_vals);
          }
        }
        else {
          for (cs_lnum_t icel = 0; icel < n_cells; icel++) {
            cs_lnum_t iel = cell_ids[icel];
            for (cs_lnum_t j = 0; j < 3; j++)
              c_vel->val[3 * iel + j] = 0.0;
          }
        }

        /* pressure initialization for groundwater model */
        if (*idarcy > 0) {
          const char *formula = NULL;

          cs_tree_node_t *tn = _find_node_variable("hydraulic_head");
          tn = cs_tree_find_node(tn, "formula");
          tn = _add_zone_id_test_attribute(tn, z->id);
          formula = cs_tree_node_get_value_str(tn);

          cs_field_t *c = cs_field_by_name_try("hydraulic_head");

          if (formula != NULL) {
            cs_real_t *ini_vals = cs_meg_initialization(z, "hydraulic_head");
            if (ini_vals != NULL) {
              for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
                cs_lnum_t c_id = cell_ids[e_id];
                c->val[c_id] = ini_vals[e_id];
              }
              BFT_FREE(ini_vals);
            }
          }

        }

        /* Turbulence variables initialization */
        const char *choice = _turbulence_initialization_choice(z_id_str);

        if (cs_gui_strcmp(choice, "formula")) {

          const char *formula_turb = NULL;
          cs_tree_node_t *tn_turb
            = cs_tree_get_node(cs_glob_tree,
                               "thermophysical_models/turbulence/initialization");
          tn_turb = _add_zone_id_test_attribute(tn_turb, z->id);
          tn_turb = cs_tree_get_node(tn_turb, "formula");
          formula_turb = cs_tree_node_get_value_str(tn_turb);

          if (formula_turb != NULL) {
            const char *model = cs_gui_get_thermophysical_model("turbulence");
            if (model == NULL)
              break;
            if (cs_gui_strcmp(model, "off"))
              break;

            cs_real_t *ini_vals = cs_meg_initialization(z, "turbulence");

            if (ini_vals != NULL) {

              if (   cs_gui_strcmp(model, "k-epsilon")
                  || cs_gui_strcmp(model, "k-epsilon-PL")) {

                cs_field_t *c_k   = cs_field_by_name("k");
                cs_field_t *c_eps = cs_field_by_name("epsilon");

                for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
                  cs_lnum_t c_id = cell_ids[e_id];
                  c_k->val[c_id]   = ini_vals[2 * e_id];
                  c_eps->val[c_id] = ini_vals[2 * e_id + 1];
                }
              }
              else if (   cs_gui_strcmp(model, "Rij-epsilon")
                       || cs_gui_strcmp(model, "Rij-SSG")) {

                cs_field_t *c_rij = cs_field_by_name_try("rij");
                cs_field_t *c_eps = cs_field_by_name("epsilon");

                if (c_rij != NULL) {
                  for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
                    cs_lnum_t c_id = cell_ids[e_id];
                    for (int drij = 0; drij < 6; drij++) {
                      c_rij->val[6*c_id + drij] = ini_vals[7*e_id + drij];
                      c_eps->val[c_id] = ini_vals[7 * e_id + 6];
                    }
                  }
                }
                else {
                  cs_field_t *c_r11 = cs_field_by_name("r11");
                  cs_field_t *c_r22 = cs_field_by_name("r22");
                  cs_field_t *c_r33 = cs_field_by_name("r33");
                  cs_field_t *c_r12 = cs_field_by_name("r12");
                  cs_field_t *c_r13 = cs_field_by_name("r13");
                  cs_field_t *c_r23 = cs_field_by_name("r23");
                  for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
                    cs_lnum_t c_id = cell_ids[e_id];
                    c_r11->val[c_id] = ini_vals[7 * e_id];
                    c_r22->val[c_id] = ini_vals[7 * e_id + 1];
                    c_r33->val[c_id] = ini_vals[7 * e_id + 2];
                    c_r12->val[c_id] = ini_vals[7 * e_id + 3];
                    c_r23->val[c_id] = ini_vals[7 * e_id + 4];
                    c_r13->val[c_id] = ini_vals[7 * e_id + 5];
                    c_eps->val[c_id] = ini_vals[7 * e_id + 6];
                  }
                }
              }
              else if (cs_gui_strcmp(model, "Rij-EBRSM")) {
                cs_field_t *c_rij = cs_field_by_name_try("rij");
                cs_field_t *c_eps = cs_field_by_name("epsilon");
                cs_field_t *c_alp = cs_field_by_name("alpha");

                if (c_rij != NULL) {
                  for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
                    cs_lnum_t c_id = cell_ids[e_id];
                    for (int drij = 0; drij < 6; drij++) {
                      c_rij->val[6*c_id + drij] = ini_vals[8*e_id + drij];
                      c_eps->val[c_id] = ini_vals[8 * e_id + 6];
                      c_alp->val[c_id] = ini_vals[8 * e_id + 7];
                    }
                  }
                }
                else {
                  cs_field_t *c_r11 = cs_field_by_name("r11");
                  cs_field_t *c_r22 = cs_field_by_name("r22");
                  cs_field_t *c_r33 = cs_field_by_name("r33");
                  cs_field_t *c_r12 = cs_field_by_name("r12");
                  cs_field_t *c_r13 = cs_field_by_name("r13");
                  cs_field_t *c_r23 = cs_field_by_name("r23");
                  for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
                    cs_lnum_t c_id = cell_ids[e_id];
                    c_r11->val[c_id] = ini_vals[8 * e_id];
                    c_r22->val[c_id] = ini_vals[8 * e_id + 1];
                    c_r33->val[c_id] = ini_vals[8 * e_id + 2];
                    c_r12->val[c_id] = ini_vals[8 * e_id + 3];
                    c_r23->val[c_id] = ini_vals[8 * e_id + 4];
                    c_r13->val[c_id] = ini_vals[8 * e_id + 5];
                    c_eps->val[c_id] = ini_vals[8 * e_id + 6];
                    c_alp->val[c_id] = ini_vals[8 * e_id + 7];
                  }
                }
              }
              else if (cs_gui_strcmp(model, "v2f-BL-v2/k")) {

                cs_field_t *c_k   = cs_field_by_name("k");
                cs_field_t *c_eps = cs_field_by_name("epsilon");
                cs_field_t *c_phi = cs_field_by_name("phi");
                cs_field_t *c_alp = cs_field_by_name("alpha");

                for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
                  cs_lnum_t c_id = cell_ids[e_id];

                  c_k->val[c_id]   = ini_vals[4 * e_id];
                  c_eps->val[c_id] = ini_vals[4 * e_id + 1];
                  c_phi->val[c_id] = ini_vals[4 * e_id + 2];
                  c_alp->val[c_id] = ini_vals[4 * e_id + 3];
                }
              }
              else if (cs_gui_strcmp(model, "k-omega-SST")) {

                cs_field_t *c_k   = cs_field_by_name("k");
                cs_field_t *c_ome = cs_field_by_name("omega");

                for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
                  cs_lnum_t c_id = cell_ids[e_id];

                  c_k->val[c_id]   = ini_vals[2 * e_id];
                  c_ome->val[c_id] = ini_vals[2 * e_id + 1];
                }
              }
              else if (cs_gui_strcmp(model, "Spalart-Allmaras")) {
                cs_field_t *c_nu = cs_field_by_name("nu_tilda");

                for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
                  cs_lnum_t c_id = cell_ids[e_id];
                  c_nu->val[c_id] = ini_vals[e_id];
                }
              }

              else
                bft_error(__FILE__, __LINE__, 0,
                          _("Invalid turbulence model: %s.\n"), model);

              BFT_FREE(ini_vals);

            }
          }
        }

        /* Thermal scalar initialization */
        if (cs_gui_thermal_model()) {

          const char *formula_sca    = NULL;
          mei_tree_t *ev_formula_sca   = NULL;
          cs_tree_node_t *tn_sca
            = cs_tree_get_node(cs_glob_tree,
                               "thermophysical_models/thermal_scalar/variable/formula");
          tn_sca = _add_zone_id_test_attribute(tn_sca, z->id);
          formula_sca = cs_tree_node_get_value_str(tn_sca);

          /* For non-specific physics defined with the GUI,
             the thermal variable can only be temperature or enthalpy
             (as the thermal model is on) */

          cs_field_t *c = cs_thermal_model_field();

          assert(c != NULL);

          if (formula_sca != NULL) {
            cs_real_t *ini_vals = cs_meg_initialization(z, "thermal");
            if (ini_vals != NULL) {
              for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
                cs_lnum_t c_id = cell_ids[e_id];
                c->val[c_id]   = ini_vals[e_id];
              }
              BFT_FREE(ini_vals);
            }
          }
          else {
            for (cs_lnum_t icel = 0; icel < n_cells; icel++) {
              cs_lnum_t iel = cell_ids[icel];
              c->val[iel] = 0.0;
            }
          }
        }

        /* User Scalars initialization */
        int n_fields = cs_field_n_fields();

        for (int f_id = 0; f_id < n_fields; f_id++) {

          const cs_field_t  *f = cs_field_by_id(f_id);

          if (f->type & CS_FIELD_USER && f->location_id == CS_MESH_LOCATION_CELLS) {

            const char *formula_sca    = NULL;

            cs_tree_node_t *tn_sca = NULL;
            tn_sca = cs_tree_get_node(cs_glob_tree, "additional_scalars/variable");
            tn_sca = cs_tree_node_get_sibling_with_tag(tn_sca, "name", f->name);
            tn_sca = cs_tree_get_node(tn_sca, "formula");
            tn_sca = _add_zone_id_test_attribute(tn_sca, z->id);
            formula_sca = cs_tree_node_get_value_str(tn_sca);

            if (formula_sca != NULL) {
              cs_real_t *ini_vals = cs_meg_initialization(z, f->name);
              if (ini_vals != NULL) {
                for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
                  cs_lnum_t c_id = cell_ids[e_id];
                  f->val[c_id] = ini_vals[e_id];
                }
                BFT_FREE(ini_vals);
              }
            }
          }
        }

        /* Meteo Scalars initialization */
        if (cs_gui_strcmp(vars->model, "atmospheric_flows")) {

          cs_tree_node_t *tn_m0
            = cs_tree_get_node(cs_glob_tree,
                               "thermophysical_models/atmospheric_flows");

          const char *name       = NULL;
          const char *formula_meteo  = NULL;
          mei_tree_t *ev_formula_meteo = NULL;

          int size = cs_tree_get_sub_node_count_simple(tn_m0, "variable");

          for (int j = 0; j < size; j++) {
            cs_tree_node_t *tn_meteo = cs_tree_get_node(tn_m0, "variable");
            for (int i = 1;
                 tn_meteo != NULL && i < j + 1;
                 i++) {
             tn_meteo = cs_tree_node_get_next_of_name(tn_meteo);
            }
            cs_tree_node_t *tn_meteo2 = tn_meteo;
            tn_meteo = cs_tree_get_node(tn_meteo, "name");
            name = cs_tree_node_get_value_str(tn_meteo);

            cs_field_t *c = cs_field_by_name_try(name);

            snprintf(z_id_str, 31, "%d", z_id);
            const char *zone_id
              = cs_tree_node_get_child_value_str(tn_meteo2, "zone_id");

            if (cs_gui_strcmp(zone_id, z_id_str))
              tn_meteo2 = cs_tree_get_node(tn_meteo2, "formula");
            else
              tn_meteo2 = NULL;

            formula_meteo = cs_tree_node_get_value_str(tn_meteo2);

            if (formula_meteo != NULL) {
              cs_real_t *ini_vals = cs_meg_initialization(z, c->name);
              if (ini_vals != NULL) {
                for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
                  cs_lnum_t c_id = cell_ids[e_id];
                  c->val[c_id] = ini_vals[e_id];
                }
                BFT_FREE(ini_vals);
              }
            }
            else {
              for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
                cs_lnum_t c_id = cell_ids[e_id];
                c->val[c_id] = 0.0;
              }
            }
          }
        }

        if (cs_gui_strcmp(vars->model, "compressible_model")) {
          const char *formula        = NULL;
          const char *buff           = NULL;
          mei_tree_t *ev_formula       = NULL;
          const char *name[] = {"pressure", "temperature", "total_energy",
                                "density"};

          ccfth = 10000;
          for (int j = 0; j < 4; j++) {

            cs_tree_node_t *tn = NULL;
            if (j < 3) {
              tn = cs_tree_find_node(cs_glob_tree, "variable");
              while (tn != NULL) {
                const char *name_tn = cs_tree_node_get_child_value_str(tn, "name");
                if (cs_gui_strcmp(name_tn, name[j]))
                  break;
                else
                  tn = cs_tree_find_node_next(cs_glob_tree, tn, "variable");
              }
            }
            else {
              tn = cs_tree_find_node(cs_glob_tree, "property");
              while (tn != NULL) {
                const char *name_tn = cs_tree_node_get_child_value_str(tn, "name");
                if (cs_gui_strcmp(name_tn, name[j]))
                  break;
                else
                  tn = cs_tree_find_node_next(cs_glob_tree, tn, "property");
              }
            }
            tn = cs_tree_get_node(tn, "formula");
            tn =_add_zone_id_test_attribute(tn, z->id);
            buff = cs_tree_node_get_child_value_str(tn, "status");

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

              formula = cs_tree_node_get_value_str(tn);

              if (formula != NULL) {
                cs_real_t *ini_vals = cs_meg_initialization(z, c->name);
                if (ini_vals != NULL) {
                  for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
                    cs_lnum_t c_id = cell_ids[e_id];
                    c->val[c_id] = ini_vals[e_id];
                  }
                  BFT_FREE(ini_vals);
                }
              }
            }

          }
          *iccfth = ccfth;
        }

      } /* END OF ISUITE == 0 */
    }
  } /* zones+1 */
}

/*----------------------------------------------------------------------------
 * User law for material properties
 *
 * Fortran Interface:
 *
 * subroutine uiphyv
 * *****************
 *
 * integer          iviscv   <--  pointer for volumic viscosity viscv
 * integer          itempk   <--  pointer for temperature (in K)
 * double precision visls0   <--  diffusion coefficient of the scalars
 * double precision viscv0   <--  volumic viscosity
 *----------------------------------------------------------------------------*/

void CS_PROCF(uiphyv, UIPHYV)(const cs_int_t  *iviscv,
                              const cs_int_t  *itempk,
                              const cs_real_t *visls0,
                              const cs_real_t *viscv0)
{
  CS_UNUSED(itempk);
  CS_UNUSED(visls0);
  CS_UNUSED(viscv0);

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const char *law = NULL;
  double time0 = cs_timer_wtime();

  cs_lnum_t i;

  cs_var_t  *vars = cs_glob_var;
  const int iscalt = cs_glob_thermal_model->iscalt;

  const cs_zone_t *z_all = cs_volume_zone_by_name("all_cells");

  /* law for density (built-in for all current integrated physical models) */
  if (cs_glob_fluid_properties->irovar == 1) {
    cs_field_t *c_rho = CS_F_(rho);
    _physical_property(c_rho, z_all);
  }

  /* law for molecular viscosity */
  if (cs_glob_fluid_properties->ivivar == 1) {
    cs_field_t *c_mu = CS_F_(mu);
    _physical_property(c_mu, z_all);
  }

  /* law for specific heat */
  if (cs_glob_fluid_properties->icp > 0) {
    cs_field_t *c_cp = CS_F_(cp);
    _physical_property(c_cp, z_all);
  }

  /* law for thermal conductivity */
  if (iscalt > 0) {

    cs_field_t  *cond_dif = NULL;

    cs_field_t *_th_f[] = {CS_F_(t), CS_F_(h), CS_F_(e_tot)};

    for (i = 0; i < 3; i++)
      if (_th_f[i]) {
        if ((_th_f[i])->type & CS_FIELD_VARIABLE) {
          int k = cs_field_key_id("diffusivity_id");
          int cond_diff_id = cs_field_get_key_int(_th_f[i], k);
          if (cond_diff_id > -1) {
            cond_dif = cs_field_by_id(cond_diff_id);
            _physical_property(cond_dif, z_all);
          }
          break;
        }
      }
    }

  /* law for volumic viscosity (compressible model) */
  if (cs_gui_strcmp(vars->model, "compressible_model")) {
    if (*iviscv > 0) {
      cs_field_t *c = cs_field_by_name_try("volume_viscosity");
      _physical_property(c, z_all);
    }
  }

  /* law for scalar diffusivity */
  int user_id = -1;
  int n_fields = cs_field_n_fields();
  const int kivisl = cs_field_key_id("diffusivity_id");
  const int kscavr = cs_field_key_id("first_moment_id");

  for (int f_id = 0; f_id < n_fields; f_id++) {

    const cs_field_t  *f = cs_field_by_id(f_id);

    if (   (f->type & CS_FIELD_VARIABLE)
        && (f->type & CS_FIELD_USER)) {
      user_id++;

      int user_law = 0;

      if (   cs_field_get_key_int(f, kscavr) < 0
          && cs_field_get_key_int(f, kivisl) >= 0) {

        char *tmp = NULL;
        BFT_MALLOC(tmp, strlen(f->name) + 13, char);
        strcpy(tmp, f->name);
        strcat(tmp, "_diffusivity");

        const char *prop_choice = _properties_choice(tmp);
        if (cs_gui_strcmp(prop_choice, "user_law"))
          user_law = 1;
        BFT_FREE(tmp);
      }

      if (user_law) {

        int diff_id = cs_field_get_key_int(f, kivisl);
        cs_field_t *c_prop = NULL;

        if (diff_id > -1)
          c_prop = cs_field_by_id(diff_id);

        /* search the formula for the law */
        cs_tree_node_t *tn
          = cs_tree_get_node(cs_glob_tree, "additional_scalars/variable");
        for (int n = 1;
             tn != NULL && n < user_id+1 ;
             n++) {
         tn = cs_tree_node_get_next_of_name(tn);
        }
        tn = cs_tree_get_node(tn, "property/formula");
        law = cs_tree_node_get_value_str(tn);

        if (law != NULL) {
          _physical_property(c_prop, z_all);
          if (cs_glob_fluid_properties->irovar == 1) {
            cs_real_t *c_rho = CS_F_(rho)->val;
            for (int c_id = 0; c_id < n_cells; c_id++)
              c_prop->val[c_id] *= c_rho[c_id];

          } else {
            for (int c_id = 0; c_id < n_cells; c_id++)
              c_prop->val[c_id] *= cs_glob_fluid_properties->ro0;
          }
        }
        cs_gui_add_mei_time(cs_timer_wtime() - time0);
      }
    }
  }

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  user_id = -1;
  for (int f_id = 0; f_id < n_fields; f_id++) {

    const cs_field_t  *f = cs_field_by_id(f_id);

    if (   (f->type & CS_FIELD_VARIABLE)
        && (f->type & CS_FIELD_USER)) {
      user_id++;

      if (   cs_field_get_key_int(f, kscavr) < 0
          && cs_field_get_key_int(f, kivisl) >= 0) {

        i = 1;
        for (tn = cs_tree_get_node(cs_glob_tree, "additional_scalars/variable");
             tn != NULL && i < user_id+1 ;
             i++) {
         tn = cs_tree_node_get_next_of_name(tn);
        }
        tn = cs_tree_get_node(tn, "property/formula");
        law = cs_tree_node_get_value_str(tn);

        bft_printf("--law for the coefficient of diffusity of the scalar %s: %s\n",
                   f->name, law);
      }
    }
  }
#endif
}

/*----------------------------------------------------------------------------
 * extra operations
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIEXOP
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiexop, UIEXOP)(void)
{
  cs_gui_balance_by_zone();
  cs_gui_pressure_drop_by_zone();
}

/*----------------------------------------------------------------------------
 * groundwater model : read laws for capacity, saturation and permeability
 *
 * Fortran Interface:
 *
 * subroutine uidapp
 * *****************
 * integer         permeability    <--  permeability type
 * integer         diffusion       <--  diffusion type
 * integer         gravity         <--  check if gravity is taken into account
 * double          gravity_x       <--  gravity direction
 * double          gravity_y       <--  gravity direction
 * double          gravity_z       <--  gravity direction
 * integer         unsaturated     <--  unsaturated zone taken into account
 *----------------------------------------------------------------------------*/

void CS_PROCF (uidapp, UIDAPP) (const int       *permeability,
                                const int       *diffusion,
                                const int       *gravity,
                                const cs_real_t *gravity_x,
                                const cs_real_t *gravity_y,
                                const cs_real_t *gravity_z,
                                const int       *unsaturated)
{
  mei_tree_t *ev_formula  = NULL;

  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->cell_cen;

  const cs_real_3_t *vel = (const cs_real_3_t *)(CS_F_(vel)->val);

  cs_field_t *fsaturation   = cs_field_by_name_try("saturation");
  cs_field_t *fcapacity     = cs_field_by_name_try("capacity");
  cs_field_t *fpermeability = cs_field_by_name_try("permeability");
  cs_field_t *fhhead     = CS_F_(head);
  cs_field_t *fsoil_density = cs_field_by_name_try("soil_density");

  cs_real_t   *saturation_field = fsaturation->val;
  cs_real_t   *capacity_field   = fcapacity->val;
  cs_real_t   *h_head_field   = fhhead->val;
  cs_real_t   *soil_density   = fsoil_density->val;

  cs_real_t     *permeability_field = NULL;
  cs_real_6_t   *permeability_field_v = NULL;

  cs_gnum_t cw[3];

  if (*permeability == 0)
    permeability_field = fpermeability->val;
  else
    permeability_field_v = (cs_real_6_t *)fpermeability->val;

  cs_tree_node_t *tn_gw
    = cs_tree_get_node(cs_glob_tree, "thermophysical_models/groundwater");

  /* number of volumic zone */

  int n_zones = cs_volume_zone_n_zones();

  for (int z_id = 0; z_id < n_zones; z_id++) {

    const cs_zone_t *z = cs_volume_zone_by_id(z_id);

    if (_zone_id_is_type(z->id, "groundwater_law")) {

      cs_tree_node_t *tn_zl = cs_tree_get_node(tn_gw, "groundwater_law");
      tn_zl = _add_zone_id_test_attribute(tn_zl, z->id);

      const cs_lnum_t n_cells = z->n_elts;
      const cs_lnum_t *cell_ids = z->elt_ids;

      char z_id_str[32];
      snprintf(z_id_str, 31, "%d", z_id);

      /* get ground properties for each zone */

      /* get soil density by zone */
      cs_real_t rhosoil = 0.;

      cs_gui_node_get_child_real(tn_zl, "soil_density", &rhosoil);

      for (cs_lnum_t icel = 0; icel < n_cells; icel++) {
        cs_lnum_t iel = cell_ids[icel];
        soil_density[iel] = rhosoil;
      }

      const char *mdl = cs_tree_node_get_child_value_str(tn_zl, "model");

      /* law for permeability */
      /* TODO: rename it in GUI, it is not Van Genuchten if saturated */

      if (cs_gui_strcmp(mdl, "VanGenuchten")) {

        cs_real_t alpha_param, ks_param, l_param, n_param;
        cs_real_t thetas_param, thetar_param;
        cs_real_t ks_xx, ks_yy, ks_zz, ks_xy, ks_xz, ks_yz;

        cs_tree_node_t *tn_m
          = cs_tree_node_get_child(tn_zl, "VanGenuchten_parameters");

        /* Van Genuchten parameters */
        if (*unsaturated) {
          cs_gui_node_get_child_real(tn_m, "alpha",  &alpha_param);
          cs_gui_node_get_child_real(tn_m, "l",      &l_param);
          cs_gui_node_get_child_real(tn_m, "n",      &n_param);
          cs_gui_node_get_child_real(tn_m, "thetar", &thetar_param);
        }

        cs_gui_node_get_child_real(tn_m, "thetas", &thetas_param);

        if (*permeability == 0)
          cs_gui_node_get_child_real(tn_m, "ks", &ks_param);
        else {
          cs_gui_node_get_child_real(tn_m, "ks_xx", &ks_xx);
          cs_gui_node_get_child_real(tn_m, "ks_yy", &ks_yy);
          cs_gui_node_get_child_real(tn_m, "ks_zz", &ks_zz);
          cs_gui_node_get_child_real(tn_m, "ks_xy", &ks_xy);
          cs_gui_node_get_child_real(tn_m, "ks_yz", &ks_yz);
          cs_gui_node_get_child_real(tn_m, "ks_xz", &ks_xz);
        }

        /* unsaturated zone considered */
        if (*unsaturated) {
          for (cs_lnum_t icel = 0; icel < n_cells; icel++) {
            cs_lnum_t iel = cell_ids[icel];
            cs_real_t head = h_head_field[iel];

            if (*gravity == 1)
              head -= (  cell_cen[iel][0] * *gravity_x
                       + cell_cen[iel][1] * *gravity_y
                       + cell_cen[iel][2] * *gravity_z);

            if (head >= 0) {
              capacity_field[iel] = 0.;
              saturation_field[iel] = thetas_param;

              if (*permeability == 0)
                permeability_field[iel] = ks_param;
              else {
                permeability_field_v[iel][0] = ks_xx;
                permeability_field_v[iel][1] = ks_yy;
                permeability_field_v[iel][2] = ks_zz;
                permeability_field_v[iel][3] = ks_xy;
                permeability_field_v[iel][4] = ks_yz;
                permeability_field_v[iel][5] = ks_xz;
              }
            }
            else {
              cs_real_t m_param = 1 - 1 / n_param;
              cs_real_t tmp1 = pow(fabs(alpha_param * head), n_param);
              cs_real_t tmp2 = 1. / (1. + tmp1);
              cs_real_t se_param = pow(tmp2, m_param);
              cs_real_t perm = pow(se_param, l_param) *
                               pow((1. - pow((1. - tmp2), m_param)), 2);

              capacity_field[iel] = -m_param * n_param * tmp1 *
                                    (thetas_param - thetar_param) *
                                     se_param * tmp2 / head;
              saturation_field[iel] = thetar_param +
                                      se_param * (thetas_param - thetar_param);

              if (*permeability == 0)
                permeability_field[iel] = perm * ks_param;
              else {
                permeability_field_v[iel][0] = perm * ks_xx;
                permeability_field_v[iel][1] = perm * ks_yy;
                permeability_field_v[iel][2] = perm * ks_zz;
                permeability_field_v[iel][3] = perm * ks_xy;
                permeability_field_v[iel][4] = perm * ks_yz;
                permeability_field_v[iel][5] = perm * ks_xz;
              }
            }
          }
        }
        else { /* saturated */
          for (cs_lnum_t icel = 0; icel < n_cells; icel++) {
            cs_lnum_t iel = cell_ids[icel];
            capacity_field[iel] = 0.;
            saturation_field[iel] = thetas_param;

            if (*permeability == 0)
              permeability_field[iel] = ks_param;
            else {
              permeability_field_v[iel][0] = ks_xx;
              permeability_field_v[iel][1] = ks_yy;
              permeability_field_v[iel][2] = ks_zz;
              permeability_field_v[iel][3] = ks_xy;
              permeability_field_v[iel][4] = ks_yz;
              permeability_field_v[iel][5] = ks_xz;
            }
          }
        }

      } else {
      /* user law for permeability */
        const char *formula
          = cs_tree_node_get_child_value_str(tn_zl, "formula");

        if (formula != NULL) {
          cs_field_t *fmeg[3] = {fcapacity, fsaturation, fpermeability};
          cs_meg_volume_function(z, fmeg);
        }
      }

      const int kivisl = cs_field_key_id("diffusivity_id");
      int n_fields = cs_field_n_fields();


      /* get diffusivity and Kd for each scalar defined by the user on current zone
         (and kplus and kminus only for scalars with kinetic model) */
      for (int f_id = 0; f_id < n_fields; f_id++) {

        cs_field_t *f = cs_field_by_id(f_id);

        if (   (f->type & CS_FIELD_VARIABLE)
            && (f->type & CS_FIELD_USER)) {

          /* get kd for current scalar and current zone */
          char *kdname = NULL;
          int len = strlen(f->name) + 4;
          BFT_MALLOC(kdname, len, char);
          strcpy(kdname, f->name);
          strcat(kdname, "_kd");
          cs_field_t *fkd = cs_field_by_name_try(kdname);
          BFT_FREE(kdname);

          cs_real_t kd_val = 0., diff_val = 0.;
          cs_tree_node_t *tn_s = cs_tree_get_node(tn_zl, "scalar");
          tn_s = cs_tree_node_get_sibling_with_tag(tn_s, "name", f->name);

          cs_gui_node_get_child_real(tn_s, "kd", &kd_val);
          cs_gui_node_get_child_real(tn_s, "diffusivity", &diff_val);

          for (cs_lnum_t icel = 0; icel < n_cells; icel++) {
            cs_lnum_t iel = cell_ids[icel];
            fkd->val[iel] = kd_val;
          }

          /* get diffusivity for current scalar and current zone */
          int diff_id = cs_field_get_key_int(f, kivisl);
          cs_field_t *fdiff = cs_field_by_id(diff_id);

          for (cs_lnum_t icel = 0; icel < n_cells; icel++) {
            cs_lnum_t iel = cell_ids[icel];
            fdiff->val[iel] = saturation_field[iel]*diff_val;
          }

          /* get kplus and kminus for current scalar and current zone
             (if EK model is chosen) */
          cs_gwf_soilwater_partition_t sorption_scal;
          int key_part = cs_field_key_id("gwf_soilwater_partition");
          cs_field_t *kp, *km;
          cs_field_get_key_struct(f, key_part, &sorption_scal);

          if (sorption_scal.kinetic == 1) {

            kp = cs_field_by_id(sorption_scal.ikp);
            km = cs_field_by_id(sorption_scal.ikm);

            cs_real_t kp_val = 0., km_val = 0.;
            cs_gui_node_get_child_real(tn_s, "kplus", &kp_val);
            cs_gui_node_get_child_real(tn_s, "kminus", &km_val);

            for (cs_lnum_t icel = 0; icel < n_cells; icel++) {
              cs_lnum_t iel = cell_ids[icel];
              kp->val[iel] = kp_val;
              km->val[iel] = km_val;
            }

          }

        }
      }

      /* get dispersion coefficient */
      if (*diffusion == 1) { /* anisotropic dispersion */
        /* TODO use a dedicated tensor field by species */

        cs_field_t *fturbvisco
          = cs_field_by_name_try("anisotropic_turbulent_viscosity");
        cs_real_6_t  *visten_v = (cs_real_6_t *)fturbvisco->val;

        cs_real_t long_diffus;
        double trans_diffus;

        cs_tree_node_t *tn = cs_tree_node_get_child(tn_zl, "diffusion_coefficient");
        cs_gui_node_get_child_real(tn, "longitudinal", &long_diffus);
        cs_gui_node_get_child_real(tn, "transverse", &trans_diffus);

        for (cs_lnum_t icel = 0; icel < n_cells; icel++) {
          cs_lnum_t iel = cell_ids[icel];
          double norm = sqrt(vel[iel][0] * vel[iel][0] +
                             vel[iel][1] * vel[iel][1] +
                             vel[iel][2] * vel[iel][2]);
          double tmp = trans_diffus * norm;
          double diff = long_diffus - trans_diffus;
          double denom = norm + 1.e-15;
          visten_v[iel][0] = tmp + diff * vel[iel][0] * vel[iel][0] / denom;
          visten_v[iel][1] = tmp + diff * vel[iel][1] * vel[iel][1] / denom;
          visten_v[iel][2] = tmp + diff * vel[iel][2] * vel[iel][2] / denom;
          visten_v[iel][3] =       diff * vel[iel][1] * vel[iel][0] / denom;
          visten_v[iel][4] =       diff * vel[iel][1] * vel[iel][2] / denom;
          visten_v[iel][5] =       diff * vel[iel][2] * vel[iel][0] / denom;
        }

      }
      else { /* isotropic dispersion */
        /* - same value of isotropic dispersion for each species
           - assigned to diffusivity field of each species
           TODO: allow to specifiy one value by species in GUI */

        double diffus;

        cs_tree_node_t *tn = cs_tree_node_get_child(tn_zl,
                                                    "diffusion_coefficient");
        cs_gui_node_get_child_real(tn, "isotropic", &diffus);

        for (int f_id = 0; f_id < n_fields; f_id++) {
          cs_field_t *f = cs_field_by_id(f_id);
          if (   (f->type & CS_FIELD_VARIABLE)
              && (f->type & CS_FIELD_USER)) {
            int diff_id = cs_field_get_key_int(f, kivisl);
            cs_field_t *fdiff = cs_field_by_id(diff_id);
            cs_real_t *visten = fdiff->val;

            /* WARNING: dispersion adds up to diffusivity
               already assigned above */
            for (cs_lnum_t l_id = 0; l_id < n_cells; l_id++) {
              cs_lnum_t iel = cell_ids[l_id];
              cs_real_t norm = sqrt(vel[iel][0] * vel[iel][0] +
                                    vel[iel][1] * vel[iel][1] +
                                    vel[iel][2] * vel[iel][2]);
              visten[iel] = visten[iel] + diffus * norm;
            }
          }
        }

      }
    }
  }

  /* check values */

  {
    const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

    cw[0] = 0; cw[1] = 0; cw[2] = 0;

    for (cs_lnum_t iel = 0; iel < n_cells; iel++) {
      if (saturation_field[iel] > 1. || saturation_field[iel] < 0.)
        cw[0] += 1;

      if (capacity_field[iel] < 0.)
        cw[1] += 1;

      if (*permeability == 0) {
        if (permeability_field[iel] < 0.)
          cw[2] += 1;
      }
    }

    cs_parall_counter(cw, 3);

    if (cw[0] > 0)
      bft_printf(_("soil_tracer_law, WARNING:\n"
                   "  saturation is outside [0, 1] in %llu cells.\n"),
                 (unsigned long long)(cw[0]));

    if (cw[1] > 0)
      bft_printf(_("soil_tracer_law, WARNING:\n"
                   "  capacity is < 0 in %llu cells.\n"),
                 (unsigned long long)(cw[1]));

    if (cw[2] > 0)
      bft_printf(_("soil_tracer_law, WARNING:\n"
                   "  isotropic permeability is < 0 in %llu cells.\n"),
                 (unsigned long long)(cw[2]));
  }
}

/*----------------------------------------------------------------------------
 * Define fans with GUI
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIFANS
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uifans, UIFANS) (void)
{
  cs_gui_define_fans();
}

/*----------------------------------------------------------------------------
 * Define error estimators
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIERES
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uieres, UIERES) (int *iescal,
                                int *iespre,
                                int *iesder,
                                int *iescor,
                                int *iestot)
{
  cs_gui_error_estimator(iescal, iespre, iesder, iescor, iestot);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize GUI reader structures.
 *----------------------------------------------------------------------------*/

void
cs_gui_init(void)
{
  if (!cs_gui_file_is_loaded())
    return;

  assert(cs_glob_var == NULL);

  BFT_MALLOC(cs_glob_var, 1, cs_var_t);

  cs_glob_var->model       = NULL;
  cs_glob_var->model_value = NULL;
}

/*-----------------------------------------------------------------------------
 * Free memory: clean global private variables.
 *----------------------------------------------------------------------------*/

void
cs_gui_finalize(void)
{
  if (!cs_gui_file_is_loaded())
    return;

  cs_gui_boundary_conditions_free_memory();

  /* clean memory for global private structure vars */

  if (cs_glob_var != NULL) {
    BFT_FREE(cs_glob_var->model);
    BFT_FREE(cs_glob_var->model_value);
    BFT_FREE(cs_glob_var);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add notebook variable to a formula.
 *
 * \param[in, out]  ev_law  pointer to MEI formula structure
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_add_notebook_variables(void  *ev_law)
{
  mei_tree_t  *_ev_law = ev_law;

  const char path0[] = "physical_properties/notebook/var";

  for (cs_tree_node_t *tn = cs_tree_get_node(cs_glob_tree, path0);
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn)) {

    const char *name = cs_tree_node_get_tag(tn, "name");
    const char *c_value = cs_tree_node_get_tag(tn, "value");

    assert(name != NULL);
    assert(c_value != NULL);

    if (name != NULL && c_value != NULL) {
      cs_real_t val = atof(c_value);
      mei_tree_insert(_ev_law, name, val);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute GUI-defined head losses for a given volume zone.
 *
 * Head loss tensor coefficients for each cell are organized as follows:
 * cku11, cku22, cku33, cku12, cku13, cku23.
 *
 * \param[in]       zone       pointer to zone structure
 * \param[in]       cvara_vel  pointer to the velocity values of the previous time step
 * \param[in, out]  cku        head loss coefficients
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_head_losses(const cs_zone_t   *zone,
                   const cs_real_3_t *cvara_vel,
                   cs_real_t          cku[][6])
{
  if (!cs_gui_file_is_loaded())
    return;

  if (! (zone->type & CS_VOLUME_ZONE_HEAD_LOSS))
    return;

  double c11, c12, c13, c21, c22, c23, c31, c32, c33;

  const cs_lnum_t n_cells = zone->n_elts;
  const cs_lnum_t *cell_ids = zone->elt_ids;

  char z_id_str[32];
  snprintf(z_id_str, 31, "%d", zone->id);

  cs_tree_node_t *tn
    = cs_tree_get_node(cs_glob_tree,
                       "thermophysical_models/head_losses/head_loss");
  tn = cs_tree_node_get_sibling_with_tag(tn, "zone_id", z_id_str);

  double k11 = _c_head_losses(tn, "kxx");
  double k22 = _c_head_losses(tn, "kyy");
  double k33 = _c_head_losses(tn, "kzz");

  double a11 = _c_head_losses(tn, "a11");
  double a12 = _c_head_losses(tn, "a12");
  double a13 = _c_head_losses(tn, "a13");
  double a21 = _c_head_losses(tn, "a21");
  double a22 = _c_head_losses(tn, "a22");
  double a23 = _c_head_losses(tn, "a23");
  double a31 = _c_head_losses(tn, "a31");
  double a32 = _c_head_losses(tn, "a32");
  double a33 = _c_head_losses(tn, "a33");

  if (   cs_gui_is_equal_real(a12, 0.0)
      && cs_gui_is_equal_real(a13, 0.0)
      && cs_gui_is_equal_real(a23, 0.0)) {

    c11 = k11;
    c22 = k22;
    c33 = k33;
    c12 = 0.0;
    c13 = 0.0;
    c23 = 0.0;

  }
  else
    _matrix_base_conversion(a11, a12, a13, a21, a22, a23, a31, a32, a33,
                            k11, 0.0, 0.0, 0.0, k22, 0.0, 0.0, 0.0, k33,
                            &c11, &c12, &c13,
                            &c21, &c22, &c23,
                            &c31, &c32, &c33);

  for (cs_lnum_t j = 0; j < n_cells; j++) {
    cs_lnum_t c_id = cell_ids[j];
    cs_real_t v = cs_math_3_norm(cvara_vel[c_id]);
    cku[j][0] = 0.5 * c11 * v;
    cku[j][1] = 0.5 * c22 * v;
    cku[j][2] = 0.5 * c33 * v;
    cku[j][3] = 0.5 * c12 * v;
    cku[j][4] = 0.5 * c23 * v;
    cku[j][5] = 0.5 * c13 * v;
  }
}

/*-----------------------------------------------------------------------------
 * Selection of linear solvers.
 *----------------------------------------------------------------------------*/

void
cs_gui_linear_solvers(void)
{
  if (!cs_gui_file_is_loaded())
    return;

  bool multigrid = false;
  cs_sles_it_type_t sles_it_type = CS_SLES_N_IT_TYPES;
  cs_multigrid_type_t mg_type = CS_MULTIGRID_N_TYPES;

  const char* algo_choice = NULL;
  const char* precond_choice = NULL;

  const int n_max_iter_default = 10000;

  int n_fields = cs_field_n_fields();

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t  *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {

      const char *ref_name = f->name;

      if (   cs_gui_strcmp(f->name, "r11")
          || cs_gui_strcmp(f->name, "r22")
          || cs_gui_strcmp(f->name, "r33")
          || cs_gui_strcmp(f->name, "r12")
          || cs_gui_strcmp(f->name, "r23")
          || cs_gui_strcmp(f->name, "r13"))
        ref_name = "rij";

      cs_tree_node_t *tn_v = _find_node_variable(ref_name);

      int n_max_iter = n_max_iter_default;
      cs_gui_node_get_child_int(tn_v, "max_iter_number", &n_max_iter);

      multigrid = false;
      sles_it_type = CS_SLES_N_IT_TYPES;

      algo_choice = _variable_choice(tn_v, "solver_choice");
      precond_choice = _variable_choice(tn_v, "preconditioning_choice");

      if (cs_gui_strcmp(algo_choice, "multigrid_k_cycle")) {
        multigrid = true;
        mg_type = CS_MULTIGRID_K_CYCLE;
      }
      else if (cs_gui_strcmp(algo_choice, "multigrid")) {
        multigrid = true;
        mg_type = CS_MULTIGRID_V_CYCLE;
      }
      else if (cs_gui_strcmp(algo_choice, "conjugate_gradient"))
        sles_it_type = CS_SLES_PCG;
      else if (cs_gui_strcmp(algo_choice, "flexible_conjugate_gradient"))
        sles_it_type = CS_SLES_FCG;
      else if (cs_gui_strcmp(algo_choice, "inexact_conjugate_gradient"))
        sles_it_type = CS_SLES_IPCG;
      else if (cs_gui_strcmp(algo_choice, "jacobi"))
        sles_it_type = CS_SLES_JACOBI;
      else if (cs_gui_strcmp(algo_choice, "bi_cgstab"))
        sles_it_type = CS_SLES_BICGSTAB;
      else if (cs_gui_strcmp(algo_choice, "bi_cgstab2"))
        sles_it_type = CS_SLES_BICGSTAB2;
      else if (cs_gui_strcmp(algo_choice, "gmres"))
        sles_it_type = CS_SLES_GMRES;
      else if (cs_gui_strcmp(algo_choice, "gauss_seidel"))
        sles_it_type = CS_SLES_P_GAUSS_SEIDEL;
      else if (cs_gui_strcmp(algo_choice, "symmetric_gauss_seidel"))
        sles_it_type = CS_SLES_P_SYM_GAUSS_SEIDEL;
      else if (cs_gui_strcmp(algo_choice, "PCR3"))
        sles_it_type = CS_SLES_PCR3;

      /* If choice is "automatic" or unspecified, delay
         choice to cs_sles_default, so do nothing here */

      if (sles_it_type < CS_SLES_N_IT_TYPES) {

        int poly_degree = 0;
        bool pc_multigrid = false;

        if (cs_gui_strcmp(precond_choice, "jacobi"))
          poly_degree = 0;
        else if (cs_gui_strcmp(precond_choice, "none"))
          poly_degree = -1;
        else if (cs_gui_strcmp(precond_choice, "polynomial"))
          poly_degree = 1;
        else if (cs_gui_strcmp(precond_choice, "multigrid_k_cycle")) {
          pc_multigrid = true;
          mg_type = CS_MULTIGRID_K_CYCLE;
          poly_degree = -1;
        }
        else if (cs_gui_strcmp(precond_choice, "multigrid")) {
          pc_multigrid = true;
          mg_type = CS_MULTIGRID_V_CYCLE;
          poly_degree = -1;
        }
        else { /* "automatic" or undefined */
          if (sles_it_type == CS_SLES_PCG) {
            pc_multigrid = true;
            mg_type = CS_MULTIGRID_V_CYCLE;
            poly_degree = -1;
          }
        }

        cs_sles_it_t *c = cs_sles_it_define(f->id, NULL, sles_it_type,
                                            poly_degree, n_max_iter);

        if (pc_multigrid) {
          cs_sles_pc_t *pc = cs_multigrid_pc_create(mg_type);
          cs_sles_it_transfer_pc(c, &pc);
        }

      }

      else if (multigrid == true) {
        cs_multigrid_t *mg = cs_multigrid_define(f->id, NULL, mg_type);

        /* If we have convection, set appropriate options */
        if (f_id >= 0) {
          cs_var_cal_opt_t var_cal_opt;
          cs_field_get_key_struct(cs_field_by_id(f_id),
                                  cs_field_key_id("var_cal_opt"),
                                  &var_cal_opt);
          if (var_cal_opt.iconv > 0)
            cs_multigrid_set_solver_options
              (mg,
               CS_SLES_P_SYM_GAUSS_SEIDEL,
               CS_SLES_P_SYM_GAUSS_SEIDEL,
               CS_SLES_P_SYM_GAUSS_SEIDEL,
               100, /* n max cycles */
               3,   /* n max iter for descent (default 2) */
               2,   /* n max iter for ascent (default 10) */
               100, /* n max iter coarse solver */
               0, 0, 0,  /* precond degree */
               -1, -1, 1); /* precision multiplier */
        }
      }
    }
  }
}

/*-----------------------------------------------------------------------------
 * Define parallel IO settings.
 *----------------------------------------------------------------------------*/

void
cs_gui_parallel_io(void)
{
  int op_id;
  int  rank_step = 0, block_size = -1;

  cs_file_mode_t op_mode[2] = {CS_FILE_MODE_READ, CS_FILE_MODE_WRITE};
  const char *op_name[2] = {"read_method", "write_method"};

  if (!cs_gui_file_is_loaded())
    return;

  cs_tree_node_t *tn_bio = cs_tree_get_node(cs_glob_tree,
                                            "calculation_management/block_io");
  /* Block IO read and write method */

  for (op_id = 0; op_id < 2; op_id++) {

    cs_file_access_t  m = CS_FILE_DEFAULT;
    const char  *method_name
      = cs_tree_node_get_child_value_str(tn_bio, op_name[op_id]);

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
    }

  }

#if defined(HAVE_MPI)

  /* Rank step and block buffer size */

  cs_gui_node_get_child_int(tn_bio, "rank_step", &rank_step);
  cs_gui_node_get_child_int(tn_bio, "min_block_size", &block_size);

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
 * Set partitioning options.
 *----------------------------------------------------------------------------*/

void
cs_gui_partition(void)
{
  cs_partition_algorithm_t a = CS_PARTITION_DEFAULT;
  bool ignore_perio = false;
  int  rank_step = 1;
  int  write_level = 1;
  int  n_add_parts = 0;
  int  *add_parts = NULL;

  if (!cs_gui_file_is_loaded())
    return;

  cs_tree_node_t *tn_p
    = cs_tree_get_node(cs_glob_tree, "calculation_management/partitioning");

  /* Partitioning type */
  const char  *part_name = cs_tree_node_get_child_value_str(tn_p, "type");

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
  }

  /* Rank step */

  cs_gui_node_get_child_int(tn_p, "rank_step", &rank_step);

  /* Ignore periodicity option */

  cs_gui_node_get_child_status_bool(tn_p, "ignore_periodicity", &ignore_perio);

  /* Output option */

  const char *s_output = cs_tree_node_get_child_value_str(tn_p, "output");

  if (s_output != NULL) {
    if (!strcmp(s_output, "no"))
      write_level = 0;
    else if (!strcmp(s_output, "default"))
      write_level = 1;
    else if (!strcmp(s_output, "yes"))
      write_level = 2;
  }

  /* List of partitions to output */

  const char *s_list = cs_tree_node_get_child_value_str(tn_p, "partition_list");

  if (s_list != NULL) {
    char *buf;
    BFT_MALLOC(buf, strlen(s_list)+1, char);
    strcpy(buf, s_list);
    char *p = strtok(buf, " \t,;");
    while (p != NULL) {
      int np = atoi(p);
      if (np > 1) {
        BFT_REALLOC(add_parts, n_add_parts + 1, int);
        add_parts[n_add_parts] = np;
        n_add_parts += 1;
      }
      p = strtok(NULL, " \t,;");
    }
    BFT_FREE(buf);
  }

  /* Set options */

  cs_partition_set_algorithm(CS_PARTITION_MAIN,
                             a,
                             rank_step,
                             ignore_perio);

  cs_partition_set_write_level(write_level);

  if (n_add_parts > 0) {
    cs_partition_add_partitions(n_add_parts, add_parts);
    BFT_FREE(add_parts);
  }
}

/*----------------------------------------------------------------------------
 * Determine porosity model type
 *----------------------------------------------------------------------------*/

void
cs_gui_porous_model(void)
{
  int n_zones = cs_volume_zone_n_zones();

  cs_tree_node_t *tn_p
    = cs_tree_get_node(cs_glob_tree,
                       "thermophysical_models/porosities/porosity");

  for (int z_id = 0; z_id < n_zones; z_id++) {
    const cs_zone_t *z = cs_volume_zone_by_id(z_id);

    if (z->type & CS_VOLUME_ZONE_POROSITY) {
      cs_tree_node_t *tn = _add_zone_id_test_attribute(tn_p, z->id);
      tn = cs_tree_get_node(tn, "model");
      const char *mdl = cs_tree_node_get_value_str(tn);

      cs_glob_porous_model = CS_MAX(1, cs_glob_porous_model);
      if (mdl) {
        if (cs_gui_strcmp(mdl, "anisotropic"))
          cs_glob_porous_model = 2;
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * 1D profile postprocessing
 *----------------------------------------------------------------------------*/

void
cs_gui_profile_output(void)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  const cs_real_3_t *cell_cen
    = (const cs_real_3_t *)(cs_glob_mesh_quantities->cell_cen);

  const cs_time_step_t *ts = cs_glob_time_step;

  static cs_real_t *t_prev = NULL;

  /* Loop on 1D profile definitions */

  int profile_id = 0;

  const char path0[] = "analysis_control/profiles/profile";

  for (cs_tree_node_t *tn = cs_tree_get_node(cs_glob_tree, path0);
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn), profile_id++) {

    const char *label = cs_gui_node_get_tag(tn, "label");

    /* for each profile, check the output frequency */

    int output_format = _get_profile_format(tn);
    const char *output_type = _get_profile_child_str(tn, "output_type");
    int output_frequency = -1;
    cs_real_t time_output = -1.;
    bool active = false;

    if (ts->nt_cur == ts->nt_prev + 1) {
      BFT_REALLOC(t_prev, profile_id+1, cs_real_t);
      t_prev[profile_id] = ts->t_prev;
    }

    if (cs_gui_strcmp(output_type, "time_value")) {
      const cs_real_t *v
        = cs_tree_node_get_child_values_real(tn, "output_frequency");
      if (v != NULL)
        time_output = v[0];
      if (   (ts->t_cur >= t_prev[profile_id] + time_output)
          || (ts->t_cur >= ts->t_max && ts->t_max > 0.)) {
        t_prev[profile_id] = ts->t_cur;
        active = true;
      }
    }
    else if (cs_gui_strcmp(output_type, "frequency")) {
      const int *v
        = cs_tree_node_get_child_values_int(tn, "output_frequency");
      if (v != NULL)
        output_frequency = v[0];
      else
        output_frequency = 1;
      if (   (ts->nt_max == ts->nt_cur)
          || (output_frequency > 0 && (ts->nt_cur % output_frequency) == 0))
        active = true;
    }
    else if (cs_gui_strcmp(output_type, "end")) {
      if (ts->nt_max == ts->nt_cur)
        active = true;
    }

    if (active) {

      FILE *file = NULL;

      cs_real_t *array = NULL;

      const char *formula = _get_profile_child_str(tn, "formula");
      mei_tree_t *ev_formula  = mei_tree_new(formula);

      mei_tree_insert(ev_formula, "s", 0.0);

      /* add variable from notebook */
      cs_gui_add_notebook_variables(ev_formula);

      /* try to build the interpreter */

      if (mei_tree_builder(ev_formula))
        bft_error(__FILE__, __LINE__, 0,
                  _("Error: can not interpret expression: %s\n %i"),
                  ev_formula->string, mei_tree_builder(ev_formula));

      const char *coord[] = {"x", "y", "z"};

      if (mei_tree_find_symbols(ev_formula, 3, coord))
        bft_error(__FILE__, __LINE__, 0,
                  _("Error: can not find the required symbol: %s\n"),
                  "x, y or z");

      int nvar_prop = cs_tree_get_node_count(tn, "var_prop");

      for (cs_tree_node_t *tn_vp = cs_tree_node_get_child(tn, "var_prop");
             tn_vp != NULL;
             tn_vp = cs_tree_node_get_next_of_name(tn_vp)) {

        const cs_field_t *f = _tree_node_get_field(tn_vp);
        int idim = _get_profile_v_component(tn_vp);
        if (f->dim > 1 && idim == -1) nvar_prop = nvar_prop + f->dim - 1;
      }

      int nvar_prop4 = nvar_prop + 4;
      BFT_MALLOC(array, nvar_prop4, cs_real_t);

      /* Only the first processor rank opens the file */

      if (cs_glob_rank_id <= 0) {

        char ts_str[16] = "";
        char fmt_str[8] = "";

        if (output_frequency > 0 || time_output > 0.) {
          snprintf(ts_str, 15, "_%.4i", ts->nt_cur); ts_str[15] = '\0';
        }
        if (output_format == 0)
          strncpy(fmt_str, ".dat", 7);
        else
          strncpy(fmt_str, ".csv", 7);

        char *filename;
        BFT_MALLOC(filename,
                   strlen(label) + strlen(ts_str) + strlen(fmt_str) + 1,
                   char);

        sprintf(filename, "%s%s%s", label, ts_str, fmt_str);

        file = fopen(filename, "w");

        if (file ==  NULL) {
          cs_base_warn(__FILE__, __LINE__);
          bft_printf( _("Unable to open the file: %s\n"), filename);
          break;
        }

        if (output_format == 0) {
          fprintf(file, "# Code_Saturne plot output (1D profile)\n#\n");
          fprintf(file, "# Iteration output: %i\n", ts->nt_cur);
          fprintf(file, "# Time output:     %12.5e\n#\n", ts->t_cur);
          fprintf(file, "#COLUMN_TITLES: Distance | X | Y | Z");
        }
        else {
          fprintf(file, "s, x, y, z");
        }

        for (cs_tree_node_t *tn_vp = cs_tree_node_get_child(tn, "var_prop");
             tn_vp != NULL;
             tn_vp = cs_tree_node_get_next_of_name(tn_vp)) {
          _write_profile_v_label_name(tn_vp, file, output_format);
        }

        fprintf(file, "\n");

        BFT_FREE(filename);
      }

      cs_lnum_t npoint = 0;
      const int *v_np = _get_profile_child_int(tn, "points");
      if (v_np != NULL)
        npoint = v_np[0];

      cs_lnum_t  c_id1 = -999, c_id = -999;
      int        rank1 = -999, c_rank = -999;
      double x1 = 0., y1 = 0., z1 = 0.;

      cs_real_t a = 1. / (double) (npoint-1);

      for (int ii = 0; ii < npoint; ii++) {

        double xx, yy, zz, xyz[3];

        cs_real_t aa = ii*a;
        mei_tree_insert(ev_formula, "s", aa);
        mei_evaluate(ev_formula);

        xyz[0] = mei_tree_lookup(ev_formula, "x");
        xyz[1] = mei_tree_lookup(ev_formula, "y");
        xyz[2] = mei_tree_lookup(ev_formula, "z");

        if (ii == 0) {
          x1 = xyz[0];
          y1 = xyz[1];
          z1 = xyz[2];
        }

        cs_geom_closest_point(n_cells, cell_cen, xyz,
                              &c_id, &c_rank);

        cs_parall_bcast(c_rank, 1, CS_LNUM_TYPE, &c_id);

        if ((c_id != c_id1) || (c_rank != rank1)) {
          c_id1 = c_id;
          rank1 = c_rank;

          if (cs_glob_rank_id == c_rank) {

            xx = cell_cen[c_id][0];
            yy = cell_cen[c_id][1];
            zz = cell_cen[c_id][2];
            array[1] = xx;
            array[2] = yy;
            array[3] = zz;
            xx = xx - x1;
            yy = yy - y1;
            zz = zz - z1;
            array[0] = sqrt(xx*xx + yy*yy + zz*zz);

            int vp_id = 4;
            for (cs_tree_node_t *tn_vp = cs_tree_node_get_child(tn, "var_prop");
                 tn_vp != NULL;
                 tn_vp = cs_tree_node_get_next_of_name(tn_vp)) {

              const cs_field_t *f = _tree_node_get_field(tn_vp);

              if (f->dim > 1) {
                int idim = _get_profile_v_component(tn_vp);
                if (idim > -1) {
                  array[vp_id] = f->val[f->dim * c_id + idim];
                  vp_id++;
                }
                else {
                  for (int ldim = 0; ldim < f->dim ; ldim++) {
                    array[vp_id] = f->val[f->dim * c_id + ldim];
                    vp_id++;
                  }
                }
              }
              else {
                array[vp_id] = f->val[c_id];
                vp_id++;
              }
            }

          }
          else {
            for (int vp_id = 0; vp_id < nvar_prop4; vp_id++)
              array[vp_id] = 0.0;
          }

          /* Send to other processors if parallel */
#if defined(HAVE_MPI)
          if (cs_glob_rank_id >= 0) {
            MPI_Bcast(array,
                      nvar_prop4,
                      CS_MPI_REAL,
                      c_rank,
                      cs_glob_mpi_comm);
          }
#endif

          if (cs_glob_rank_id <= 0) {
            if (output_format == 0) {
              for (int vp_id = 0; vp_id < nvar_prop4; vp_id++)
                fprintf(file, "%12.5e ", array[vp_id]);
              fprintf(file, "\n");
            }
            else {
              if (nvar_prop > 0) {
                for (int vp_id = 0; vp_id < nvar_prop4 - 1; vp_id++)
                  fprintf(file, "%12.5e, ", array[vp_id]);
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

  if (   (ts->nt_max == ts->nt_cur)
      || (ts->t_cur >= ts->t_max && ts->t_max > 0.)) {
    BFT_FREE(t_prev);
  }
}

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
  cs_tree_node_t *tn = cs_tree_find_node(cs_glob_tree, "property");
  while (tn != NULL) {
    const char *name_tn = cs_tree_node_get_child_value_str(tn, "name");
    if (cs_gui_strcmp(name_tn, property_name))
      break;
    else
      tn = cs_tree_find_node_next(cs_glob_tree, tn, "property");
  }
  tn = cs_tree_get_node(tn, "initial_value");

  cs_gui_node_get_real(tn, value);
}

/*-----------------------------------------------------------------------------
 * Get value of reference fluid properties parameter.
 *
 * parameters:
 *   name            <--   parameter name
 *   value           -->   parameter value
 *----------------------------------------------------------------------------*/

void
cs_gui_fluid_properties_value(const char  *param,
                              double      *value)
{
  cs_tree_node_t *tn
    = cs_tree_get_node(cs_glob_tree, "physical_properties/fluid_properties");
  tn = cs_tree_get_node(tn, param);

  cs_gui_node_get_real(tn, value);
}

/*----------------------------------------------------------------------------
 * Get thermal scalar model.
 *
 * return:
 *   value of itherm*10 + (temperature variant flag)
 *----------------------------------------------------------------------------*/

int
cs_gui_thermal_model(void)
{
  int   test = 0;

  const char *model = cs_gui_get_thermophysical_model("thermal_scalar");

  if (cs_gui_strcmp(model, "off"))
    test = 0;
  else {
    if (cs_gui_strcmp(model, "enthalpy"))
      test = 20;
    else if (cs_gui_strcmp(model, "temperature_kelvin"))
      test = 11;
    else if (cs_gui_strcmp(model, "temperature_celsius"))
      test = 10;
    else if (cs_gui_strcmp(model, "potential_temperature"))
      test = 12;
    else if (cs_gui_strcmp(model, "liquid_potential_temperature"))
      test = 13;
    else if (cs_gui_strcmp(model, "total_energy"))
      test = 30;
    else
      bft_error(__FILE__, __LINE__, 0,
          _("Invalid thermal model: %s\n"), model);
  }

  return test;
}

/*----------------------------------------------------------------------------
 * Time moments definition
 *----------------------------------------------------------------------------*/

void
cs_gui_time_moments(void)
{
  if (!cs_gui_file_is_loaded())
    return;

  int imom = 1;
  int restart = cs_restart_present();

  /* Loop on time average definitions */

  const char path0[] = "/analysis_control/time_averages/time_average";

  for (cs_tree_node_t *tn = cs_tree_get_node(cs_glob_tree, path0);
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn), imom++) {

    const char *restart_name;
    cs_time_moment_restart_t  restart_mode = CS_TIME_MOMENT_RESTART_AUTO;

    const int *v_i;
    const cs_real_t *v_r;

    /* Older files used "label", now use "name", so try both */

    const char *m_name = cs_tree_node_get_tag(tn, "name");

    if (m_name == NULL) {
      m_name = cs_tree_node_get_tag(tn, "label");
      if (m_name == NULL) /* if neither found, force error case */
        m_name = cs_gui_node_get_tag(tn, "name");
    }

    v_i = cs_tree_node_get_child_values_int(tn, "time_step_start");
    int nt_start = (v_i != NULL) ? v_i[0] : 0;

    v_r = cs_tree_node_get_child_values_real(tn, "time_start");
    double t_start = (v_r != NULL) ? v_r[0] : -1;

    /* test on restart */

    if (restart != 0) {
      v_i = cs_tree_node_get_child_values_int(tn, "restart_from_time_average");
      int restart_id = (v_i != NULL) ? v_i[0] : -2;
      cs_time_moment_restart_options_by_id(restart_id,
                                           &restart_mode,
                                           &restart_name);
    }

    int n_m_fields = cs_tree_get_node_count(tn, "var_prop");

    int *m_f_id;
    BFT_MALLOC(m_f_id, n_m_fields*2, int);
    int *m_c_id = m_f_id + n_m_fields;

    int j = 0;
    for (cs_tree_node_t *tn_vp = cs_tree_node_get_child(tn, "var_prop");
         tn_vp != NULL;
         tn_vp = cs_tree_node_get_next_of_name(tn_vp), j++) {

      const char *f_name = cs_gui_node_get_tag(tn_vp, "name");
      v_i = cs_tree_node_get_child_values_int(tn_vp, "component");
      int idim = (v_i != NULL) ? v_i[0] : -1;

      cs_field_t *f = cs_field_by_name_try(f_name);

      /* If we failed to find Rij, we search for Rxx.
       * This test is needed for the case where irijco = 0 */

      if (f == NULL && cs_gui_strcmp(f_name, "rij")) {
        switch(idim) {
        case 0:
          f = CS_F_(r11);
          break;
        case 1:
          f = CS_F_(r22);
          break;
        case 2:
          f = CS_F_(r33);
          break;
        case 3:
          f = CS_F_(r12);
          break;
        case 4:
          f = CS_F_(r23);
          break;
        case 5:
          f = CS_F_(r13);
          break;
        }
        m_f_id[j] = f->id;
        m_c_id[j] = 0;
      }

      else {
        m_f_id[j] = f->id;
        m_c_id[j] = idim;
      }

    }

    cs_time_moment_define_by_field_ids(m_name,
                                       n_m_fields,
                                       m_f_id,
                                       m_c_id,
                                       CS_TIME_MOMENT_MEAN,
                                       nt_start,
                                       t_start,
                                       restart_mode,
                                       restart_name);

    m_c_id = NULL;
    BFT_FREE(m_f_id);

  }

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
#endif
}

/*-----------------------------------------------------------------------------
 * Set turbomachinery model
 *----------------------------------------------------------------------------*/

void
cs_gui_turbomachinery(void)
{
  cs_turbomachinery_model_t  model_type;
  bool coupled;

  _turbomachinery_model(&model_type, &coupled);

  cs_turbomachinery_set_model(model_type);
}

/*-----------------------------------------------------------------------------
 * Set turbomachinery options.
 *----------------------------------------------------------------------------*/

void
cs_gui_turbomachinery_rotor(void)
{
  if (!cs_gui_file_is_loaded())
    return;

  cs_turbomachinery_model_t  model_type;
  bool coupled;

  _turbomachinery_model(&model_type, &coupled);

  if (model_type != CS_TURBOMACHINERY_NONE) {

    cs_tree_node_t *tn = NULL;
    cs_tree_node_t *tn2 = NULL;
    int n_rotors =
      cs_tree_get_node_count(cs_glob_tree,
                             "/thermophysical_models/turbomachinery/rotor");

    for (int rotor_id = 0; rotor_id < n_rotors; rotor_id++) {

      double rotation_axis[3];
      double rotation_invariant[3];
      double rotation_velocity;

      const char *cell_criteria;

      rotation_axis[0] = _rotor_option(rotor_id, "axis_x");
      rotation_axis[1] = _rotor_option(rotor_id, "axis_y");
      rotation_axis[2] = _rotor_option(rotor_id, "axis_z");

      rotation_invariant[0] = _rotor_option(rotor_id, "invariant_x");
      rotation_invariant[1] = _rotor_option(rotor_id, "invariant_y");
      rotation_invariant[2] = _rotor_option(rotor_id, "invariant_z");

      tn = cs_tree_get_node(cs_glob_tree,
                              "thermophysical_models/turbomachinery/rotor");
      int i;
      for (i = 1;
           tn != NULL && i < rotor_id + 1 ;
           i++) {
       tn = cs_tree_node_get_next_of_name(tn);
      }
      tn2 = tn;
      tn = cs_tree_get_node(tn, "velocity/value");
      cs_gui_node_get_real(tn, &rotation_velocity);

      tn2 = cs_tree_get_node(tn2, "criteria");
      cell_criteria = cs_tree_node_get_value_str(tn2);

      cs_turbomachinery_add_rotor(cell_criteria,
                                  rotation_velocity,
                                  rotation_axis,
                                  rotation_invariant);

    }

    int n_join = cs_tree_get_node_count(cs_glob_tree,
                                        "/thermophysical_models"
                                        "/turbomachinery/joining/face_joining");

    for (int join_id = 0; join_id < n_join; join_id++) {

      const char *selector_s  =  _get_rotor_face_joining("selector", join_id+1);
      const char *fraction_s  =  _get_rotor_face_joining("fraction", join_id+1);
      const char *plane_s     =  _get_rotor_face_joining("plane", join_id+1);
      const char *verbosity_s = _get_rotor_face_joining("verbosity", join_id+1);
      const char *visu_s =  _get_rotor_face_joining("visualization", join_id+1);

      double fraction = (fraction_s != NULL) ? atof(fraction_s) : 0.1;
      double plane = (plane_s != NULL) ? atof(plane_s) : 25.0;
      int verbosity = (verbosity_s != NULL) ? atoi(verbosity_s) : 0;
      int visualization = (visu_s != NULL) ? atoi(visu_s) : 0;

      if (coupled == false)
        (void)cs_turbomachinery_join_add(selector_s,
                                         fraction,
                                         plane,
                                         verbosity,
                                         visualization);
      else
        (void)cs_turbomachinery_coupling_add(selector_s,
                                             fraction,
                                             verbosity);

    }

  }
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

  if (mei_wtime > 0.0) {
    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("\nTime elapsed defining values using MEI: %12.5f\n"),
                  mei_wtime);
    cs_log_printf(CS_LOG_PERFORMANCE, "\n");
    cs_log_separator(CS_LOG_PERFORMANCE);
  }
}

/*----------------------------------------------------------------------------
 * Define user variables through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_user_variables(void)
{
  int i = 0;

  const char *t_scalar_name = NULL; /* thermal scalar name if present */

  const char path_s[] = "additional_scalars/variable";
  cs_tree_node_t *tn_s = cs_tree_get_node(cs_glob_tree, path_s);

  for (cs_tree_node_t *tn = tn_s;
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn), i++) {

    if (i == 0 && cs_glob_thermal_model->itherm != CS_THERMAL_MODEL_NONE) {
      const char path_t[] = "thermophysical_models/thermal_scalar/variable";
      t_scalar_name = cs_tree_node_get_tag
                        (cs_tree_get_node(cs_glob_tree, path_t), "name");
    }

    const char *name = cs_gui_node_get_tag(tn, "name");

    const char *variance_name = cs_tree_node_get_child_value_str(tn, "variance");

    /* In case of variance, check for presence of matching field
       in thermal and user scalars */

    if (variance_name != NULL) {

      bool found = false;
      if (t_scalar_name != NULL) {
        if (strcmp(t_scalar_name, variance_name) == 0)
          found = true;
      }
      for (cs_tree_node_t *tn_c = tn_s;
           tn_c != NULL && found == false;
           tn_c = cs_tree_node_get_next_of_name(tn_c), i++) {
        const char *cmp_name = cs_tree_node_get_tag(tn_c, "name");
        if (cmp_name != NULL) {
          if (strcmp(cmp_name, variance_name) == 0)
            found = true;
        }
      }

      if (found)
        cs_parameters_add_variable_variance(name, variance_name);

    }

    /* If not a variance, we have a regular variable */

    else
      cs_parameters_add_variable(name, 1);
  }
}

/*----------------------------------------------------------------------------
 * Define user arrays through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_user_arrays(void)
{
  const char path_s[] = "additional_scalars/users/property";
  cs_tree_node_t *tn_s = cs_tree_get_node(cs_glob_tree, path_s);

  for (cs_tree_node_t *tn = tn_s;
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn)) {

    const char *name = cs_gui_node_get_tag(tn, "name");

    int array_dim = 1;
    cs_tree_node_t *dtn = cs_tree_get_node(tn, "dimension");
    cs_gui_node_get_int(dtn, &array_dim);

    const char *location_name = cs_gui_node_get_tag(tn, "support");

    if (strcmp(location_name, "cells") == 0)
      cs_parameters_add_property(name, array_dim, CS_MESH_LOCATION_CELLS);

    else if (strcmp(location_name, "internal") == 0)
      cs_parameters_add_property(name, array_dim, CS_MESH_LOCATION_INTERIOR_FACES);

    else if (strcmp(location_name, "boundary") == 0)
      cs_parameters_add_property(name, array_dim, CS_MESH_LOCATION_BOUNDARY_FACES);

    else if (strcmp(location_name, "vertices") == 0)
      cs_parameters_add_property(name, array_dim, CS_MESH_LOCATION_VERTICES);

  }
}

/*----------------------------------------------------------------------------
 * Define volume and boundary zones through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_zones(void)
{
  if (!cs_gui_file_is_loaded())
    return;

  /* Ensure zones ordering for safety (should be removed in the future)*/

  _ensure_zones_order();

  int id = 0;

  const char default_criteria[] = "all[]";

  /* Volume zones */
  /*------------- */

  cs_tree_node_t *tn_vc = cs_tree_get_node(cs_glob_tree,
                                           "solution_domain/volumic_conditions");

  const int n_v_zones = cs_tree_get_node_count(tn_vc, "zone");

  /* Build ordering array to check zones are defined in increasing id order */

  cs_lnum_t *order = NULL, *z_ids = NULL;
  BFT_MALLOC(order, n_v_zones, cs_lnum_t);
  BFT_MALLOC(z_ids, n_v_zones, cs_lnum_t);

  /* Loop on volume condition zones */

  id = 0;
  for (cs_tree_node_t *tn = cs_tree_node_get_child(tn_vc, "zone");
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn), id++) {
    z_ids[id] = _v_zone_t_id(tn, id);
  }

  assert(id == n_v_zones);

  cs_order_lnum_allocated(NULL, z_ids, order, n_v_zones);

  /* Now loop on zones in id order */

  for (int i = 0; i < n_v_zones; i++) {

    int type_flag = 0, z_id = z_ids[order[i]];

    cs_tree_node_t *tn = _v_zone_node_by_id(tn_vc, z_id);

    /* zone name */

    const char *name = cs_tree_node_get_tag(tn, "label");

    /* location criteria */

    const char *_criteria = cs_tree_node_get_value_str(tn);
    const char *criteria = (_criteria != NULL) ? _criteria : default_criteria;

    /* Check for initialization */

    if (_zone_is_type(tn, "initialization"))
      type_flag = type_flag | CS_VOLUME_ZONE_INITIALIZATION;

    /* Check for porosity */

    if (_zone_is_type(tn, "porosity"))
      type_flag = type_flag | CS_VOLUME_ZONE_POROSITY;

    /* Check for head losses */

    if (_zone_is_type(tn, "head_losses"))
      type_flag = type_flag | CS_VOLUME_ZONE_HEAD_LOSS;

    /* Check for source terms */

    if (_zone_is_type(tn, "momentum_source_term"))
      type_flag = type_flag | CS_VOLUME_ZONE_SOURCE_TERM;
    if (_zone_is_type(tn, "scalar_source_term"))
      type_flag = type_flag | CS_VOLUME_ZONE_SOURCE_TERM;
    if (_zone_is_type(tn, "thermal_source_term"))
      type_flag = type_flag | CS_VOLUME_ZONE_SOURCE_TERM;

    /* Finally, define zone */

    cs_volume_zone_define(name, criteria, type_flag);
  }

  BFT_FREE(order);
  BFT_FREE(z_ids);

  /* Boundary zones */
  /*--------------- */

  /* Loop on boundary condition zones */

  cs_tree_node_t *tn_bc = cs_tree_get_node(cs_glob_tree,
                                           "boundary_conditions");

  id = 0;
  for (cs_tree_node_t *tn = cs_tree_node_get_child(tn_bc, "boundary");
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn), id++) {

    /* Zone id in tree; note that the name tag for boundary zones actually
       defines an integer (1 to n). This tag should be removed in the
       future to only use the zone label (the actual zone name). */

    const char *id_s = cs_tree_node_get_tag(tn, "name");
    if (id_s != NULL) {
      int z_t_id = atoi(id_s);
      if (z_t_id != id + 1)
        bft_printf(_("\n"
                     " Warning: noncontiguous %s zone ids in XML:\n"
                     "          zone with index %d has id %d.\n"),
                   tn->name, id, z_t_id);
    }

    /* Zone name */

    const char *name = cs_tree_node_get_tag(tn, "label");

    /* location criteria */

    const char *_criteria = cs_tree_node_get_value_str(tn);
    const char *criteria = (_criteria != NULL) ? _criteria : default_criteria;

    int type_flag = 0;

    /* Define zone */

    cs_boundary_zone_define(name, criteria, type_flag);
  }

}

/*----------------------------------------------------------------------------
 * Define balance by zone through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_balance_by_zone(void)
{
  const char path0[] = "/analysis_control/scalar_balances/scalar_balance";

  for (cs_tree_node_t *tn = cs_tree_get_node(cs_glob_tree, path0);
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn)) {

    const char _default_criteria[] = "all[]";

    const char *criteria = cs_tree_node_get_child_value_str(tn, "criteria");
    if (criteria == NULL) criteria = _default_criteria;

    for (cs_tree_node_t *tn_v = cs_tree_node_get_child(tn, "var_prop");
         tn_v != NULL;
         tn_v = cs_tree_node_get_next_of_name(tn_v)) {

      const char *name = cs_gui_node_get_tag(tn_v, "name");
      cs_balance_by_zone(criteria, name);

    }
  }
}

/*----------------------------------------------------------------------------
 * Define pressure drop through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_pressure_drop_by_zone(void)
{
  const char path0[] = "/analysis_control/scalar_balances/pressure_drop";

  for (cs_tree_node_t *tn = cs_tree_get_node(cs_glob_tree, path0);
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn)) {

    const char _default_criteria[] = "all[]";

    const char *criteria = cs_tree_node_get_child_value_str(tn, "criteria");
    if (criteria == NULL) criteria = _default_criteria;

    cs_pressure_drop_by_zone(criteria);
  }
}

/*----------------------------------------------------------------------------
 * Define fans through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_define_fans(void)
{
  if (!cs_gui_file_is_loaded())
    return;

  const char path0[] = "thermophysical_models/fans/fan";

  for (cs_tree_node_t *tn = cs_tree_get_node(cs_glob_tree, path0);
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn)) {

    const int *v_i;
    const cs_real_t *v_r;

    const char *i_axis_s[] = {"inlet_axis_x", "inlet_axis_y", "inlet_axis_z"};
    const char *o_axis_s[] = {"outlet_axis_x", "outlet_axis_y", "outlet_axis_z"};
    const char *p_coeff_s[]
      = {"curve_coeffs_x", "curve_coeffs_y", "curve_coeffs_z"};

    v_i = cs_tree_node_get_child_values_int(tn, "mesh_dimension");
    int dim = (v_i != NULL) ? v_i[0] : 3;

    cs_real_t inlet_axis_coords[3] = {0, 0, 0};
    cs_real_t outlet_axis_coords[3] = {0.1, 0, 0};
    cs_real_t pressure_curve_coeffs[3] = {0.6, -0.1, -0.05};

    for (int i = 0; i < 3; i++) {
      v_r = cs_tree_node_get_child_values_real(tn, i_axis_s[i]);
      if (v_r != NULL) inlet_axis_coords[i] = v_r[0];
    }
    for (int i = 0; i < 3; i++) {
      v_r = cs_tree_node_get_child_values_real(tn, o_axis_s[i]);
      if (v_r != NULL) outlet_axis_coords[i] = v_r[0];
    }

    v_r = cs_tree_node_get_child_values_real(tn, "fan_radius");
    cs_real_t fan_radius = (v_r != NULL) ? v_r[0] : 0.7;

    v_r = cs_tree_node_get_child_values_real(tn, "blades_radius");
    cs_real_t blades_radius = (v_r != NULL) ? v_r[0] : 0.5;

    v_r = cs_tree_node_get_child_values_real(tn, "hub_radius");
    cs_real_t hub_radius = (v_r != NULL) ? v_r[0] : 0.1;

    v_r = cs_tree_node_get_child_values_real(tn, "axial_torque");
    cs_real_t axial_torque = (v_r != NULL) ? v_r[0] : 0.01;

    for (int i = 0; i < 3; i++) {
      v_r = cs_tree_node_get_child_values_real(tn, p_coeff_s[i]);
      if (v_r != NULL) pressure_curve_coeffs[i] = v_r[0];
    }

    cs_fan_define(dim, /* fan (mesh) dimension (2D or 3D) */
                  inlet_axis_coords,
                  outlet_axis_coords,
                  fan_radius,
                  blades_radius,
                  hub_radius,
                  pressure_curve_coeffs,
                  axial_torque);
  }
}

/*----------------------------------------------------------------------------
 * Define error estimator through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_error_estimator(int *iescal,
                       int *iespre,
                       int *iesder,
                       int *iescor,
                       int *iestot)
{
  cs_tree_node_t *tn_ee
    = cs_tree_get_node(cs_glob_tree, "analysis_control/error_estimator");

  cs_tree_node_t *tn = cs_tree_get_node(tn_ee, "Correction/model");

  const char *result = cs_tree_node_get_value_str(tn);

  if (cs_gui_strcmp(result, "1"))
    iescal[*iescor -1] = 1;
  else if (cs_gui_strcmp(result, "2"))
    iescal[*iescor -1] = 2;
  else
    iescal[*iescor -1] = 0;

  tn = cs_tree_get_node(tn_ee, "Drift/model");

  result = cs_tree_node_get_value_str(tn);

  if (cs_gui_strcmp(result, "1"))
    iescal[*iesder -1] = 1;
  else if (cs_gui_strcmp(result, "2"))
    iescal[*iesder -1] = 2;
  else
    iescal[*iesder -1] = 0;

  tn = cs_tree_get_node(tn_ee, "Prediction/model");

  result = cs_tree_node_get_value_str(tn);

  if (cs_gui_strcmp(result, "1"))
    iescal[*iespre -1] = 1;
  else if (cs_gui_strcmp(result, "2"))
    iescal[*iespre -1] = 2;
  else
    iescal[*iespre -1] = 0;

  tn = cs_tree_get_node(tn_ee, "Total/model");

  result = cs_tree_node_get_value_str(tn);

  if (cs_gui_strcmp(result, "1"))
    iescal[*iestot -1] = 1;
  else if (cs_gui_strcmp(result, "2"))
    iescal[*iestot -1] = 2;
  else
    iescal[*iestot -1] = 0;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
