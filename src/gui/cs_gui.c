/*============================================================================
 * Management of the GUI parameters file: main parameters
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

#include "cs_all_to_all.h"
#include "cs_array.h"
#include "cs_base.h"
#include "cs_boundary.h"
#include "cs_boundary_zone.h"
#include "cs_cf_model.h"
#include "cs_equation.h"
#include "cs_equation_param.h"
#include "cs_ext_neighborhood.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_file.h"
#include "cs_function.h"
#include "cs_log.h"
#include "cs_gui_util.h"
#include "cs_gui_boundary_conditions.h"
#include "cs_gui_specific_physics.h"
#include "cs_gui_mobile_mesh.h"
#include "cs_geom.h"
#include "cs_internal_coupling.h"
#include "cs_math.h"
#include "cs_meg_prototypes.h"
#include "cs_meg_xdef_wrapper.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_location.h"
#include "cs_multigrid.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_porous_model.h"
#include "cs_parameters.h"
#include "cs_param_sles.h"
#include "cs_partition.h"
#include "cs_physical_model.h"
#include "cs_rotation.h"
#include "cs_selector.h"
#include "cs_timer.h"
#include "cs_time_moment.h"
#include "cs_time_table.h"
#include "cs_thermal_model.h"
#include "cs_physical_properties.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_sles.h"
#include "cs_sles_it.h"
#include "cs_turbulence_model.h"
#include "cs_wall_functions.h"
#include "cs_physical_constants.h"
#include "cs_balance_by_zone.h"
#include "cs_fan.h"
#include "cs_velocity_pressure.h"
#include "cs_vof.h"
#include "cs_volume_zone.h"

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

/* MEG contexts */

static int                            _n_v_meg_contexts = 0;
static cs_gui_volume_meg_context_t  **_v_meg_contexts = NULL;

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Static local variables
 *============================================================================*/

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
 * Return the tree node associated to a given property and zone.
 *----------------------------------------------------------------------------*/

static cs_tree_node_t *
_get_property_node(const char *property_name,
                   const char *zone_name)
{
  cs_tree_node_t *tn = cs_tree_find_node(cs_glob_tree, "property");
  while (tn != NULL) {
    const char *name = cs_tree_node_get_child_value_str(tn, "name");
    if (cs_gui_strcmp(name, property_name))
      break;
    else
      tn = cs_tree_find_node_next(cs_glob_tree, tn, "property");
  }

  /* FIXME: Currently zone definitions (other than all_cells) are
   * handled using sub-nodes. This will be changed for v7.1
   */
  if (zone_name != NULL) {
    if (strcmp(zone_name, "all_cells")) {
      /* Get zone_id from xml */
      cs_tree_node_t *id_n =
        cs_tree_get_node(cs_glob_tree,
                         "solution_domain/volumic_conditions/zone");
      const char *id_s = NULL;
      while (id_n != NULL) {
        const char *zname = cs_tree_node_get_tag(id_n, "label");
        if (cs_gui_strcmp(zname, zone_name)) {
          id_s = cs_tree_node_get_tag(id_n, "id");
          break;
        }
        id_n = cs_tree_node_get_next_of_name(id_n);
      }
      if (id_s != NULL) {
        tn = cs_tree_node_get_child(tn, "zone");
        tn = cs_tree_node_get_sibling_with_tag(tn, "zone_id", id_s);
      }
    }
  }

  return tn;
}

/*----------------------------------------------------------------------------
 * Return the value of the choice attribute from a property name.
 *
 * parameters:
 *   property_name        <--  name of the property
 *   zone_name            <--  name of zone. If NULL we use
 *                             default definition
 *----------------------------------------------------------------------------*/

static const char*
_properties_choice(const char *property_name,
                   const char *zone_name)
{
  const char *choice = NULL;

  cs_tree_node_t *node = _get_property_node(property_name, zone_name);

  choice = cs_tree_node_get_child_value_str(node, "choice");

  return choice;
}

/*----------------------------------------------------------------------------
 * Return the formula associated to a given property and zone.
 *----------------------------------------------------------------------------*/

static const char*
_property_formula(const char *property_name,
                  const char *zone_name)
{
  const char *law = NULL;

  cs_tree_node_t *node = _get_property_node(property_name, zone_name);

  node = cs_tree_node_get_child(node, "formula");

  law = cs_tree_node_get_value_str(node);

  return law;
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

  const char *prop_choice = _properties_choice(name, NULL);
  if (cs_gui_strcmp(prop_choice, "thermal_law"))
    choice = 1;

  return choice;
}

/*-----------------------------------------------------------------------------
 * Compute a physical property based on a thermal law.
 *----------------------------------------------------------------------------*/

static void
_physical_property_thermal_law(cs_field_t           *c_prop,
                               const cs_zone_t      *z,
                               cs_phys_prop_type_t   property)
{
  /* For incompressible flows, the thermodynamic pressure is constant over
   * time and is the reference pressure. */

  cs_lnum_t thermodynamic_pressure_stride = 0;
  cs_lnum_t thermal_f_val_stride = 1;
  cs_real_t _p0 = cs_glob_fluid_properties->p0;
  cs_real_t _t0 = cs_glob_fluid_properties->t0;

  const cs_real_t *thermodynamic_pressure = &_p0;
  const cs_real_t *_thermal_f_val = NULL;

  /* For variable density flows with pressure dependent density (idilat > 1)
   * we use total pressure values insted of the reference pressue P0. */
  cs_field_t *p_tot_field = cs_field_by_name_try("total_pressure");
  if (p_tot_field != NULL && cs_glob_velocity_pressure_model->idilat > 1) {
    thermodynamic_pressure = p_tot_field->val;
    thermodynamic_pressure_stride = 1;
  }

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

  cs_phys_prop_compute(property,
                       z->n_elts,
                       thermodynamic_pressure_stride,
                       thermal_f_val_stride,
                       thermodynamic_pressure,
                       _thermal_f_val,
                       c_prop->val);
}

/*-----------------------------------------------------------------------------
 * use MEG for thermal diffusivity
 *
 * This is a special case of _physical_property_, where the thermal
 * diffusivity computation is based on that of the thermal conductivity.
 *----------------------------------------------------------------------------*/

static void
_physical_property_th_diffusivity(cs_field_t          *c_prop,
                                  const cs_zone_t     *z)
{
  const char prop_name[] = "thermal_conductivity";

  const cs_lnum_t n_cells = z->n_elts;
  const cs_lnum_t *cell_ids = z->elt_ids;

  int user_law = 0;
  const char *prop_choice = _properties_choice(prop_name, NULL);

  if (cs_gui_strcmp(prop_choice, "user_law"))
    user_law = 1;
  else if (z->id > 1 && z->type & CS_VOLUME_ZONE_PHYSICAL_PROPERTIES) {
    const char *law = _property_formula(prop_name, z->name);
    if (law != NULL)
      user_law = 1;
  }

  if (user_law) {

    const char *law = _property_formula(prop_name, z->name);

    if (law != NULL) {
      /* Pass "field overlay": shallow copy of field with modified name
       so as to be handled by the MEG volume function. */
      cs_field_t _c_prop = *c_prop;
      _c_prop.name = prop_name;

      const cs_real_3_t *restrict cell_cen =
        (const cs_real_3_t *restrict)cs_glob_mesh_quantities->cell_cen;

      cs_meg_volume_function(z->name,
                             n_cells,
                             cell_ids,
                             cell_cen,
                             _c_prop.name,
                             &(_c_prop.val));
    }

  }
  else if (cs_gui_strcmp(prop_choice, "thermal_law") &&
           cs_gui_strcmp(z->name, "all_cells")) {

    _physical_property_thermal_law(c_prop,
                                   z,
                                   CS_PHYS_PROP_THERMAL_CONDUCTIVITY);

  }

  /* If property is globaly variable but defined as a constant
     for this zone, set from constant */

  else if (z->id <= 1 && cs_gui_strcmp(prop_choice, "constant")) {

    const cs_fluid_properties_t *phys_pp = cs_glob_fluid_properties;
    const cs_real_t p_val = phys_pp->lambda0;

    for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
      cs_lnum_t c_id = cell_ids[e_id];
      c_prop->val[c_id] = p_val;
    }

  }

  /* If property is predefined by the model, do nothing and return. */

  else if (cs_gui_strcmp(prop_choice, "predefined_law"))
    return;

  /* Finalize special case for conduction to diffusion conversion */

  if (CS_F_(cp) == NULL) {
    const cs_real_t cp0 = cs_glob_fluid_properties->cp0;
    for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
      cs_lnum_t c_id = cell_ids[e_id];
        c_prop->val[c_id] /= cp0;
    }
  }
  else {
    const cs_real_t *cp = CS_F_(cp)->val;
    for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
      cs_lnum_t c_id = cell_ids[e_id];
      c_prop->val[c_id] /= cp[c_id];
    }
  }
}

/*-----------------------------------------------------------------------------
 * use MEI for physical property
 *----------------------------------------------------------------------------*/

static void
_physical_property(cs_field_t          *c_prop,
                   const cs_zone_t     *z)
{
  int user_law = 0;
  const char *prop_choice = _properties_choice(c_prop->name, NULL);

  /* Special case: "thermal_conductivity" is defined by GUI, but
     in case of Enthalpy thermal model, the matching field is
     "thermal_diffusivity". */

  if (prop_choice == NULL) {
    if (   cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_ENTHALPY
        && cs_gui_strcmp(c_prop->name, "thermal_diffusivity")) {
      _physical_property_th_diffusivity(c_prop, z);
      return;
    }
  }

  if (cs_gui_strcmp(prop_choice, "user_law"))
    user_law = 1;
  else if (z->id > 1 && z->type & CS_VOLUME_ZONE_PHYSICAL_PROPERTIES) {
    const char *law = _property_formula(c_prop->name, z->name);
    if (law != NULL)
      user_law = 1;
  }

  if (user_law) {

    const char *law = _property_formula(c_prop->name, z->name);

    if (law != NULL) {
      const cs_real_3_t *restrict cell_cen =
        (const cs_real_3_t *restrict)cs_glob_mesh_quantities->cell_cen;
      cs_meg_volume_function(z->name,
                             z->n_elts,
                             z->elt_ids,
                             cell_cen,
                             c_prop->name,
                             &(c_prop->val));
    }

  }
  else if (cs_gui_strcmp(prop_choice, "thermal_law") &&
           cs_gui_strcmp(z->name, "all_cells")) {
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

    _physical_property_thermal_law(c_prop, z, property);
  }

  /* If property is globaly variable but defined as a constant
     for this zone, set from constant */

  else if (z->id <= 1 && cs_gui_strcmp(prop_choice, "constant")) {
    cs_real_t p_val = -1;
    const cs_fluid_properties_t *phys_pp = cs_glob_fluid_properties;

    if (cs_gui_strcmp(c_prop->name, "density"))
      p_val = phys_pp->ro0;
    else if (cs_gui_strcmp(c_prop->name, "molecular_viscosity"))
      p_val = phys_pp->viscl0;
    else if (cs_gui_strcmp(c_prop->name, "specific_heat"))
      p_val = phys_pp->cp0;
    else if (cs_gui_strcmp(c_prop->name, "thermal_conductivity")) {
      p_val = phys_pp->lambda0;
    }

    const cs_lnum_t n_cells = z->n_elts;
    const cs_lnum_t *cell_ids = z->elt_ids;
    for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
      cs_lnum_t c_id = cell_ids[e_id];
      c_prop->val[c_id] = p_val;
    }
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

    if (cs_gui_strcmp(buff, "user_law"))
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
 *   f       <-- pointer to field
 *   value   <-- value of diffusion coefficient
 *----------------------------------------------------------------------------*/

static void
_scalar_diffusion_value(const cs_field_t  *f,
                        cs_real_t         *value)
{
  cs_tree_node_t *tn
    = cs_tree_get_node(cs_glob_tree, "additional_scalars/variable");

  while (tn != NULL) {
    const char *name = cs_tree_node_get_child_value_str(tn, "name");
    if (cs_gui_strcmp(name, f->name))
      break;
    else
      tn = cs_tree_find_node_next(cs_glob_tree, tn, "variable");
  }

  tn = cs_tree_get_node(tn, "property/initial_value");

  cs_gui_node_get_real(tn, value);
}

/*----------------------------------------------------------------------------
 * Check if at least one user property is based on a thermal law.
 *
 * Fortran Interface:
 *
 * subroutine uiphyv
 * *****************
 *
 * integer          iviscv   <--  pointer for volumic viscosity viscv
 *----------------------------------------------------------------------------*/

static int
_count_thermal_laws(void)
{
  int n_thl = 0;

  int n_zones = cs_volume_zone_n_zones();
  for (int z_id = 0; z_id < n_zones; z_id++) {
    const cs_zone_t *z = cs_volume_zone_by_id(z_id);
    if (! (z->type & CS_VOLUME_ZONE_PHYSICAL_PROPERTIES))
      continue;
    if (strcmp(z->name, "all_cells") != 0)
      continue;

    const char *names[] = {
      "density", "molecular_viscosity", "specific_heat",
      "thermal_conductivity"
    };

    for (int i = 0; i < 4; i++) {
      const char *prop_choice = _properties_choice(names[i], NULL);
      if (cs_gui_strcmp(prop_choice, "thermal_law"))
        n_thl += 1;
    }

    /* volumic viscosity (compressible model) */
    if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1) {
      const char *prop_choice = _properties_choice("volume_viscosity", NULL);
      if (cs_gui_strcmp(prop_choice, "thermal_law"))
        n_thl += 1;
    }
  }

  return n_thl;
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

  if (cs_gui_strcmp(choice, "upwind")) /* for old xml compatibility */
    *keyword = 1;
  if (cs_gui_strcmp(choice, "centered"))
    *keyword = 1;
  else if (cs_gui_strcmp(choice, "solu"))
    *keyword = 0;
  else if (cs_gui_strcmp(choice, "solu_upwind_gradient"))
    *keyword = 2;
  else if (cs_gui_strcmp(choice, "blending"))
    *keyword = 3;
  else if (cs_gui_strcmp(choice, "nvd_tvd"))
    *keyword = 4;
}

/*----------------------------------------------------------------------------
 * Get the attribute value from the NVD limiter
 *
 * parameters:
 *   tn_v    <-- node assocaited with variable
 *   keyword -->  value of attribute node
 *----------------------------------------------------------------------------*/

static void
_nvd_limiter_scheme_value(cs_tree_node_t  *tn_v,
                          int             *keyword)
{
  cs_tree_node_t *tn = cs_tree_node_get_child(tn_v, "nvd_limiter");
  const char *choice = cs_tree_node_get_child_value_str(tn, "choice");

  static const char *names[] = {"gamma", "smart", "cubista", "superbee",
                                "muscl", "minmod", "clam", "stoic",
                                "osher", "waseb", "hric", "cicsam", "stacs"};

  if (choice != NULL) {
    for (int i = 0; i < 13; i++) {
      if (strcmp(choice, names[i]) == 0) {
        *keyword = i;
        break;
      }
    }
  }
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
  cs_tree_node_t *tn_c = cs_tree_node_get_child(tn, "slope_test");
  if (tn_c != NULL) {
    cs_gui_node_get_status_int(tn_c, &result);

    if (result == 1)
      *keyword = 0;
    if (result == 0)
      *keyword = 1;

    else {
      const char *choice = cs_tree_node_get_tag(tn_c, "choice");
      if (choice != NULL) {
        if (strcmp(choice, "beta_limiter") == 0)
          *keyword = 2;
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Get the attribute value for a variable's cell gradient
 *
 * parameters:
 *   tn_v    <-- node assocaited with variable
 *   kw      <-- associated keyword
 *   keyword -->  value of attribute node
 *----------------------------------------------------------------------------*/

static void
_var_gradient_type(cs_tree_node_t  *tn_v,
                   const char      *kw,
                   int             *keyword)
{
  const char *choice = cs_tree_node_get_child_value_str(tn_v, kw);

  if (choice != NULL) {
    /* Default: "global" for cells, "automatic" for boundary faces, or none */

    if (strcmp(choice, "green_iter") == 0)
      *keyword = 0;
    else if (strcmp(choice, "lsq") == 0)
      *keyword = 1;
    else if (strcmp(choice, "lsq_ext") == 0)
      *keyword = 2;
    else if (strcmp(choice, "green_lsq") == 0)
      *keyword = 4;
    else if (strcmp(choice, "green_lsq_ext") == 0)
      *keyword = 5;
    else if (strcmp(choice, "green_vtx") == 0)
      *keyword = 7;
  }
}

/*----------------------------------------------------------------------------
 * Get the attribute value for a variable's cell gradient
 *
 * parameters:
 *   tn_v    <-- node assocaited with variable
 *   keyword -->  value of attribute node
 *----------------------------------------------------------------------------*/

static void
_var_gradient_limiter_type(cs_tree_node_t  *tn_v,
                           int             *keyword)
{
  const char *choice = cs_tree_node_get_child_value_str(tn_v,
                                                        "gradient_limiter_type");

  if (choice != NULL) {
    /* Default: "none" for -1 */

    if (strcmp(choice, "cell") == 0)
      *keyword = 0;
    else if (strcmp(choice, "face") == 0)
      *keyword = 1;
  }
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

  } else if (cs_gui_strcmp(param, "piso_sweep_number")) {

    tn = cs_tree_get_node(tn, "velocity_pressure_algo");
    const char *algo_choice = cs_tree_node_get_child_value_str(tn, "choice");
    if (cs_gui_strcmp(algo_choice, "piso")) {
      tn = cs_tree_get_node(tn, param);
      cs_gui_node_get_int(tn, &result);
      *keyword = result;
    }

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

  buff = _properties_choice(property_name, NULL);
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

/*----------------------------------------------------------------------------
 * Read reference dynamic and user scalar viscosity
 *----------------------------------------------------------------------------*/

static void
_read_diffusivity(void)
{
  double result, density;

  const int kscavr = cs_field_key_id("first_moment_id");
  const int kvisls0 = cs_field_key_id("diffusivity_ref");

  const int itherm = cs_glob_thermal_model->itherm;

  cs_fluid_properties_t *fprops = cs_get_glob_fluid_properties();

  if (itherm != CS_THERMAL_MODEL_NONE) {

    if (_thermal_table_needed("thermal_conductivity") == 0)
      cs_gui_properties_value("thermal_conductivity", &(fprops->lambda0));
    else
      cs_phys_prop_compute(CS_PHYS_PROP_THERMAL_CONDUCTIVITY,
                           1,
                           0,
                           0,
                           &(cs_glob_fluid_properties->p0),
                           &(cs_glob_fluid_properties->t0),
                           &(fprops->lambda0));

    double visls_0 = fprops->lambda0;

    /* for the Temperature, the diffusivity factor is not divided by Cp
     * i.e. it remains lambda */
    if (itherm != CS_THERMAL_MODEL_TEMPERATURE)
      visls_0 /= cs_glob_fluid_properties->cp0;

    cs_field_t *tf = cs_thermal_model_field();
    cs_field_set_key_double(tf, kvisls0, visls_0);

  }

  /* User scalar
     In the interface, the user gives the diffusion coefficient, whereas in
     the solver, one sets the diffusivity, thus one need to multiply
     this coefficient by the density to remain coherent */

  if (cs_glob_physical_model_flag[CS_GROUNDWATER] < 0) {
    int n_fields = cs_field_n_fields();
    for (int f_id = 0; f_id < n_fields; f_id++) {
      cs_field_t  *f = cs_field_by_id(f_id);
      if (   (f->type & CS_FIELD_VARIABLE)
          && (f->type & CS_FIELD_USER)) {
        if (cs_field_get_key_int(f, kscavr) < 0) {

          if (cs_glob_physical_model_flag[CS_COMBUSTION_COAL] > -1) {
            /* Air molar mass */
            result = 0.028966;
            cs_gui_fluid_properties_value("reference_molar_mass", &result);
            if (result <= 0)
              bft_error
                (__FILE__, __LINE__, 0,
                 _("mass molar value is zero or not found in the xml file.\n"));
            density = cs_glob_fluid_properties->p0 *
                      result / (8.31446 *(cs_glob_fluid_properties->t0));
          }
          else
            density = cs_glob_fluid_properties->ro0;

          double visls_0 = cs_field_get_key_double(f, kvisls0);
          double coeff = visls_0 / density;
          _scalar_diffusion_value(f, &coeff);
          visls_0 = coeff * density;

          cs_field_set_key_double(f, kvisls0, visls_0);
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Get cs_mesh_location type from string
 *----------------------------------------------------------------------------*/
static cs_mesh_location_type_t
_cs_mesh_location_type_from_str(const char *location_name)
{
  cs_mesh_location_type_t loc = CS_MESH_LOCATION_NONE;

  if (strcmp(location_name, "cells") == 0)
    loc = CS_MESH_LOCATION_CELLS;

  else if (strcmp(location_name, "internal") == 0)
    loc = CS_MESH_LOCATION_INTERIOR_FACES;

  else if (strcmp(location_name, "boundary") == 0)
    loc = CS_MESH_LOCATION_BOUNDARY_FACES;

  else if (strcmp(location_name, "vertices") == 0)
    loc = CS_MESH_LOCATION_VERTICES;

  return loc;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Specific heat variable or constant indicator.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CFPPVA
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (csfpva, CSFPVA) (void)
{
  int choice;
  cs_fluid_properties_t *phys_pp = cs_get_glob_fluid_properties();

  if (_properties_choice_id("specific_heat", &choice))
    phys_pp->icp = (choice > 0) ? 0 : -1;

  if (_properties_choice_id("volume_viscosity", &choice))
    phys_pp->iviscv = (choice > 0) ? 0 : -1;

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--icp = %i\n", phys_pp->icp);
  bft_printf("--iviscv = %i\n", phys_pp->iviscv);
#endif
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

  const int keysca = cs_field_key_id("scalar_id");
  const int kivisl = cs_field_key_id("diffusivity_id");
  const int kscavr = cs_field_key_id("first_moment_id");
  const int n_fields = cs_field_n_fields();

  cs_field_t *tf = cs_thermal_model_field();

  if (   cs_glob_physical_model_flag[CS_PHYSICAL_MODEL_FLAG] <= 0
      && tf != NULL) {
    test1 = _properties_choice_id("thermal_conductivity", &choice1);
    test2 = _properties_choice_id("specific_heat", &choice2);

    if (test1 && test2) {

      for (int f_id = 0; f_id < n_fields; f_id++) {
        cs_field_t  *f = cs_field_by_id(f_id);
        if (f->type & CS_FIELD_VARIABLE) {
          if (f == tf) {
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
        && (f->type & CS_FIELD_USER)
        && (f != tf)) {
      int iscal = cs_field_get_key_int(f, keysca);
      if (iscal > 0) {
        if (cs_field_get_key_int(f, kscavr) < 0) {
          if (_scalar_properties_choice(f->name, &choice1))
            cs_field_set_key_int(f, kivisl, choice1 - 1);
          // for groundwater we impose variable property
          if (cs_glob_physical_model_flag[CS_GROUNDWATER] > -1)
            cs_field_set_key_int(f, kivisl, 0);
        }
      }
    }
  }

  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1) {
    int d_f_id = -1;
    /*FIXME:CK*/
    const char *prop_choice = _properties_choice("thermal_conductivity", NULL);
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
  cs_velocity_pressure_param_t *vp_param
    = cs_get_glob_velocity_pressure_param();
  int result = vp_param->iphydr;
  cs_tree_node_t *tn
    = cs_tree_find_node(cs_glob_tree,
                        "numerical_parameters/hydrostatic_pressure");
  cs_gui_node_get_status_int(tn, &result);
  vp_param->iphydr = result;

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--iphydr = %i\n", vp_model->iphydr);
#endif
}

/*----------------------------------------------------------------------------
 * Hydrostatic equilibrium parameter.
 *
 * Fortran Interface:
 *
 * subroutine cscfgp
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (cscfgp, CSCFGP) (void)
{
  cs_cf_model_t *cf_model = cs_get_glob_cf_model();

  int result = cf_model->icfgrp;
  cs_tree_node_t *tn
    = cs_tree_find_node(cs_glob_tree,
                        "numerical_parameters/hydrostatic_equilibrium/");
  cs_gui_node_get_status_int(tn, &result);

  cf_model->icfgrp = result;

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--icfgrp = %i\n", cf_model->icfgrp);
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
  /* Default, forbidden values for time step factor */
  double cdtmin = -1., cdtmax = -1.;

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

  if (cdtmin > 0)
    time_opt->dtmin = cdtmin * time_stp->dt_ref;
  if (cdtmax > 0)
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

  cs_gui_node_get_child_status_int(tn,
                                   "thermal_time_step",
                                   &(time_opt->iptlro));

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
 * Define porosity.
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIPORO
 *----------------------------------------------------------------------------*/

void CS_PROCF(uiporo, UIPORO)(void)
{
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  int n_zones = cs_volume_zone_n_zones();

  /* Porosity fields */
  cs_field_t *fporo = CS_F_(poro);
  cs_field_t *ftporo = CS_F_(t_poro);

  if (fporo != NULL)
    cs_array_real_set_scalar(n_cells_ext, 1., fporo->val);

  if (ftporo != NULL) {
    cs_real_6_t *porosf = (cs_real_6_t *)ftporo->val;
    for (cs_lnum_t iel = 0; iel < n_cells_ext; iel++) {
      porosf[iel][0] = 1.;
      porosf[iel][1] = 1.;
      porosf[iel][2] = 1.;
      porosf[iel][3] = 0.;
      porosf[iel][4] = 0.;
      porosf[iel][5] = 0.;
    }
  }

  cs_tree_node_t *tn_p
    = cs_tree_get_node(cs_glob_tree,
                       "thermophysical_models/porosities/porosity");

  for (int z_id = 0; z_id < n_zones; z_id++) {
    const cs_zone_t *z = cs_volume_zone_by_id(z_id);

    if (z->type & CS_VOLUME_ZONE_POROSITY) {

      cs_tree_node_t *tn_zp = _add_zone_id_test_attribute(tn_p, z->id);
      const char *mdl = cs_tree_node_get_child_value_str(tn_zp, "model");
      const char *formula = cs_tree_node_get_child_value_str(tn_zp, "formula");

      if (formula != NULL) {
        const cs_real_3_t *restrict cell_cen =
          (const cs_real_3_t *restrict)cs_glob_mesh_quantities->cell_cen;

        if (cs_gui_strcmp(mdl, "anisotropic")) {
          char _name_string[512];
          snprintf(_name_string, 511, "%s+%s", fporo->name, ftporo->name);
          _name_string[511] = '\0';
          cs_real_t *fvals[2] = {fporo->val, ftporo->val};
          cs_meg_volume_function(z->name,
                                 z->n_elts,
                                 z->elt_ids,
                                 cell_cen,
                                 _name_string,
                                 fvals);

        }
        else {
          cs_meg_volume_function(z->name,
                                 z->n_elts,
                                 z->elt_ids,
                                 cell_cen,
                                 fporo->name,
                                 &(fporo->val));
        }

      }
    }
  }

  cs_porous_model_auto_face_porosity();
}

/*----------------------------------------------------------------------------
 * User law for material properties
 *
 * Fortran Interface:
 *
 * subroutine uiphyv
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF(uiphyv, UIPHYV)(void)
{
  double time0 = cs_timer_wtime();

  int n_zones_pp
    = cs_volume_zone_n_type_zones(CS_VOLUME_ZONE_PHYSICAL_PROPERTIES);
  int n_zones = cs_volume_zone_n_zones();
  /* law for density (built-in for all current integrated physical models) */
  if (cs_glob_fluid_properties->irovar == 1) {
    cs_field_t *c_rho = CS_F_(rho);
    if (n_zones_pp > 0) {
      for (int z_id = 0; z_id < n_zones; z_id++) {
        const cs_zone_t *z = cs_volume_zone_by_id(z_id);
        if (z->type & CS_VOLUME_ZONE_PHYSICAL_PROPERTIES)
          _physical_property(c_rho, z);
      }
    }
  }

  /* law for molecular viscosity */
  if (cs_glob_fluid_properties->ivivar == 1) {
    cs_field_t *c_mu = CS_F_(mu);
    if (n_zones_pp > 0) {
      for (int z_id = 0; z_id < n_zones; z_id++) {
        const cs_zone_t *z = cs_volume_zone_by_id(z_id);
        if (z->type & CS_VOLUME_ZONE_PHYSICAL_PROPERTIES)
          _physical_property(c_mu, z);
      }
    }
  }

  /* law for specific heat */
  if (cs_glob_fluid_properties->icp > 0) {
    cs_field_t *c_cp = CS_F_(cp);
    if (n_zones_pp > 0) {
      for (int z_id = 0; z_id < n_zones; z_id++) {
        const cs_zone_t *z = cs_volume_zone_by_id(z_id);
        if (z->type & CS_VOLUME_ZONE_PHYSICAL_PROPERTIES)
          _physical_property(c_cp, z);
      }
    }
  }

  /* law for thermal conductivity */
  if (cs_glob_thermal_model->itherm != CS_THERMAL_MODEL_NONE) {

    cs_field_t  *cond_dif = NULL;

    cs_field_t *_th_f[] = {CS_F_(t), CS_F_(h), CS_F_(e_tot)};

    for (int i = 0; i < 3; i++)
      if (_th_f[i]) {
        if ((_th_f[i])->type & CS_FIELD_VARIABLE) {
          int k = cs_field_key_id("diffusivity_id");
          int cond_diff_id = cs_field_get_key_int(_th_f[i], k);
          if (cond_diff_id > -1) {
            cond_dif = cs_field_by_id(cond_diff_id);
            if (n_zones_pp > 0) {
              for (int z_id = 0; z_id < n_zones; z_id++) {
                const cs_zone_t *z = cs_volume_zone_by_id(z_id);
                if (z->type & CS_VOLUME_ZONE_PHYSICAL_PROPERTIES)
                  _physical_property(cond_dif, z);
              }
            }
          }
          break;
        }
      }
  }

  /* law for volumic viscosity (compressible model) */
  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1) {
    if (cs_glob_fluid_properties->iviscv > 0) {
      cs_field_t *c = cs_field_by_name_try("volume_viscosity");
      if (n_zones_pp > 0 && c != NULL) {
        for (int z_id = 0; z_id < n_zones; z_id++) {
          const cs_zone_t *z = cs_volume_zone_by_id(z_id);
          if (z->type & CS_VOLUME_ZONE_PHYSICAL_PROPERTIES)
            _physical_property(c, z);
        }
      }
    }
  }

  /* law for scalar diffusivity */
  int n_fields = cs_field_n_fields();
  const int kivisl = cs_field_key_id("diffusivity_id");
  const int kscavr = cs_field_key_id("first_moment_id");

  for (int f_id = 0; f_id < n_fields; f_id++) {

    const cs_field_t  *f = cs_field_by_id(f_id);

    if (   (f->type & CS_FIELD_VARIABLE)
        && (f->type & CS_FIELD_USER)) {

      if (   cs_field_get_key_int(f, kscavr) < 0
          && cs_field_get_key_int(f, kivisl) >= 0) {

        /* Get diffusivity pointer.
         * diff_id is >= 0 since it is a part of the if test
         * above!
         */
        int diff_id = cs_field_get_key_int(f, kivisl);
        cs_field_t *c_prop = cs_field_by_id(diff_id);

        if (n_zones_pp > 0) {
          for (int z_id = 0; z_id < n_zones; z_id++) {
            const cs_zone_t *z = cs_volume_zone_by_id(z_id);
            if (z->type & CS_VOLUME_ZONE_PHYSICAL_PROPERTIES) {
              const char *law = _property_formula(c_prop->name, z->name);
              if (law != NULL) {
                _physical_property(c_prop, z);
                if (cs_glob_fluid_properties->irovar == 1) {
                  cs_real_t *c_rho = CS_F_(rho)->val;
                  for (cs_lnum_t e_id = 0; e_id < z->n_elts; e_id++) {
                    cs_lnum_t c_id = z->elt_ids[e_id];
                    c_prop->val[c_id] *= c_rho[c_id];
                  }
                }
                else {
                  for (cs_lnum_t e_id = 0; e_id < z->n_elts; e_id++) {
                    cs_lnum_t c_id = z->elt_ids[e_id];
                    c_prop->val[c_id] *= cs_glob_fluid_properties->ro0;
                  }
                }
                cs_gui_add_mei_time(cs_timer_wtime() - time0);
#if _XML_DEBUG
                bft_printf("==> %s\n", __func__);
                bft_printf("--law for the diffusivity coefficient "
                           "of the scalar '%s' over zone '%s':\n  %s\n",
                           f->name, z->name, law);
#endif
              }
            }
          }
        }
      }
    }
  }

}

/*----------------------------------------------------------------------------
 * extra operations
 *
 * Fortran Interface:
 *
 * subroutine uiexop
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiexop, UIEXOP)(void)
{
  cs_gui_balance_by_zone();
  cs_gui_pressure_drop_by_zone();
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read GUi-defined Checkpoint parameters.
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_checkpoint_parameters(void)
{
  int   nt_interval = 0;
  double  t_interval, wt_interval = -1;

  cs_restart_checkpoint_get_intervals(&nt_interval, &t_interval, &wt_interval);

  cs_tree_node_t *tn = cs_tree_get_node(cs_glob_tree,
                                        "calculation_management/start_restart");

  cs_gui_node_get_child_int(tn, "restart_rescue", &nt_interval);

  cs_restart_checkpoint_set_interval(nt_interval, t_interval, wt_interval);

  cs_gui_node_get_child_status_int
    (tn, "restart_with_auxiliary",
     &(cs_glob_restart_auxiliary->read_auxiliary));

  cs_time_scheme_t *t_sch = cs_get_glob_time_scheme();

  cs_gui_node_get_child_status_int(tn, "frozen_field", &(t_sch->iccvfg));

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--nt_interval = %d\n", nt_interval);
  bft_printf("--read_auxiliary = %d\n",
             cs_glob_restart_auxiliary->read_auxiliary);
  bft_printf("--iccvfg = %d\n", t_sch->iccvfg);
#endif
}

/*----------------------------------------------------------------------------
 * Space scheme options, linear solver precision and time step factor
 *----------------------------------------------------------------------------*/

void
cs_gui_equation_parameters(void)
{
  const int keysca = cs_field_key_id("scalar_id");
  const int k_ts_fact = cs_field_key_id("time_step_factor");

  int n_fields = cs_field_n_fields();
  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      cs_equation_param_t *eqp = cs_field_get_equation_param(f);
      cs_equation_param_t *eqp_eq = NULL;

      cs_tree_node_t *tn_v = _find_node_variable(f->name);

      cs_gui_node_get_child_real(tn_v, "solver_precision", &(eqp->epsilo));
      cs_gui_node_get_child_status_int(tn_v, "flux_reconstruction",
                                       &(eqp->ircflu));
      cs_gui_node_get_child_int(tn_v, "rhs_reconstruction",
                                &(eqp->nswrsm));
      cs_gui_node_get_child_int(tn_v, "verbosity", &(eqp->verbosity));

      /* For CDO equation, if non-automatic value ie != -1 */
      if (eqp_eq == NULL)
        eqp_eq = cs_equation_param_by_name(f->name);

      if (   eqp_eq != NULL && cs_gui_is_equal_real(eqp->epsilo, -1) == 0
          && eqp_eq->sles_param != NULL)
        eqp_eq->sles_param->cvg_param.rtol = eqp->epsilo;

      /* convection scheme options */
      if (eqp->iconv > 0) {
        cs_gui_node_get_child_real(tn_v, "blending_factor",
                                   &(eqp->blencv));
        _order_scheme_value(tn_v, &(eqp->ischcv));
        if (eqp->ischcv == 4) {
          int nvd_limiter = -1;
          _nvd_limiter_scheme_value(tn_v, &nvd_limiter);
          if (nvd_limiter >= 0) {
            cs_field_set_key_int(f,
                                 cs_field_key_id("limiter_choice"),
                                 nvd_limiter);
          }
        }
        _slope_test_value(tn_v, &(eqp->isstpc));
      }

      /* Gradient options */
      {
        _var_gradient_type(tn_v, "cell_gradient_type",
                           &(eqp->imrgra));
        _var_gradient_type(tn_v, "boundary_gradient_type",
                           &(eqp->b_gradient_r));
        cs_gui_node_get_child_real(tn_v, "gradient_epsilon",
                                   &(eqp->epsrgr));
        _var_gradient_limiter_type(tn_v, &(eqp->imligr));
        cs_gui_node_get_child_real(tn_v, "gradient_limiter_factor",
                                   &(eqp->climgr));
      }

      /* only for additional variables (user or model) */
      int isca = cs_field_get_key_int(f, keysca);
      if (isca > 0) {
        /* time step factor */
        double cdtvar = cs_field_get_key_double(f, k_ts_fact);
        cs_gui_node_get_child_real(tn_v, "time_step_factor", &cdtvar);
        cs_field_set_key_double(f, k_ts_fact, cdtvar);
      }
    }
  }

#if _XML_DEBUG_
  const int var_key_id = cs_field_key_id("variable_id");
  bft_printf("==> %s\n", __func__);
  for (int f_id = 0; f_id < n_fields; f_id++) {
    const cs_field_t  *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      j = cs_field_get_key_int(f, var_key_id) -1;
      bft_printf("-->variable[%i] = %s\n", j, f->name);
      bft_printf("--blencv = %f\n", eqp->blencv);
      bft_printf("--epsilo = %g\n", eqp->epsilo);
      bft_printf("--ischcv = %i\n", eqp->ischcv);
      bft_printf("--isstpc = %i\n", eqp->isstpc);
      bft_printf("--ircflu = %i\n", eqp->ircflu);
      bft_printf("--nswrsm = %i\n", eqp->nswrsm);
      if (cs_field_get_key_int(f, keysca) > 0) {
        double cdtvar = cs_field_get_key_double(f, k_ts_fact);
        cs_gui_node_get_child_real(tn_v, "time_step_factor", &cdtvar);
        bft_printf("--cdtvar = %g\n", cdtvar);
      }
    }
  }
#endif
}

/*-----------------------------------------------------------------------------
 * Free memory: clean global private variables.
 *----------------------------------------------------------------------------*/

void
cs_gui_finalize(void)
{
  cs_gui_boundary_conditions_free_memory();

  /* Clean MEG contexts */

  for (int i = 0; i < _n_v_meg_contexts; i++)
    BFT_FREE(_v_meg_contexts[i]);

  BFT_FREE(_v_meg_contexts);
}

/*----------------------------------------------------------------------------
 * Return a pointer to equation parameters based on a field or equation name.
 *
 * parameters:
 *   name <-- field or equation name
 *
 * return:
 *   pointer to matching child string
 *----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_gui_get_equation_param(const char  *name)
{
  cs_equation_param_t *eqp = NULL;

  cs_field_t *f = cs_field_by_name_try(name);
  if (f != NULL)
    eqp = cs_field_get_equation_param(f);

  else
    eqp = cs_equation_param_by_name(name);

  return eqp;
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
 * Groundwater model : read laws for capacity, saturation and permeability
 *
 * parameters:
 *   permeability  <--  permeability type
 *   diffusion     <--  diffusion type
 *   unsaturated   <--  unsaturated zone taken into account
 *----------------------------------------------------------------------------*/

void
cs_gui_groundwater_property_laws(int  permeability,
                                 int  diffusion,
                                 int  unsaturated)
{
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

  if (permeability == 0)
    permeability_field = fpermeability->val;
  else
    permeability_field_v = (cs_real_6_t *)fpermeability->val;

  cs_tree_node_t *tn_gw
    = cs_tree_get_node(cs_glob_tree, "thermophysical_models/groundwater");

  /* number of volumic zones */

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
        if (unsaturated) {
          cs_gui_node_get_child_real(tn_m, "alpha",  &alpha_param);
          cs_gui_node_get_child_real(tn_m, "l",      &l_param);
          cs_gui_node_get_child_real(tn_m, "n",      &n_param);
          cs_gui_node_get_child_real(tn_m, "thetar", &thetar_param);
        }

        cs_gui_node_get_child_real(tn_m, "thetas", &thetas_param);

        if (permeability == 0)
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
        if (unsaturated) {

          const cs_physical_constants_t *phys_cst = cs_glob_physical_constants;
          cs_real_t g_norm = cs_math_3_norm(phys_cst->gravity);
          bool gravity = (g_norm > cs_math_epzero) ? true : false;
          const cs_real_t g_n[3] = {phys_cst->gravity[0]/g_norm,
                                    phys_cst->gravity[1]/g_norm,
                                    phys_cst->gravity[2]/g_norm};

          for (cs_lnum_t icel = 0; icel < n_cells; icel++) {
            cs_lnum_t iel = cell_ids[icel];
            cs_real_t head = h_head_field[iel];

            if (gravity)
              head -= cs_math_3_dot_product(cell_cen[iel], g_n);

            if (head >= 0) {
              capacity_field[iel] = 0.;
              saturation_field[iel] = thetas_param;

              if (permeability == 0)
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

              if (permeability == 0)
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

            if (permeability == 0)
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

      }
      else {
      /* user law for permeability */
        const char *formula
          = cs_tree_node_get_child_value_str(tn_zl, "formula");

        if (formula != NULL) {
          char _name_string[512];
          snprintf(_name_string, 511, "%s+%s+%s",
                   fcapacity->name, fsaturation->name, fpermeability->name);
          _name_string[511] = '\0';
          cs_real_t *fvals[3] = {fcapacity->val,
                                 fsaturation->val,
                                 fpermeability->val};
          cs_meg_volume_function(z->name,
                                 n_cells,
                                 cell_ids,
                                 cell_cen,
                                 _name_string,
                                 fvals);
        }
      }

      const int kivisl = cs_field_key_id("diffusivity_id");
      int n_fields = cs_field_n_fields();

      /* get diffusivity and Kd for each scalar defined by the user on
         current zone (and clsat only for scalars
         with kinetic model) */
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

          /* get clsat for current scalar and current zone */
          cs_field_t *kp, *km;

          if (false) {  /* TODO: kd + precipitation */

            cs_real_t kp_val = 0., km_val = 0.;
            cs_gui_node_get_child_real(tn_s, "clsat", &kp_val);

            for (cs_lnum_t icel = 0; icel < n_cells; icel++) {
              cs_lnum_t iel = cell_ids[icel];
              kp->val[iel] = kp_val;
              km->val[iel] = km_val;
            }

          }

        }
      }

      /* get dispersion coefficient */
      if (diffusion == 1) { /* anisotropic dispersion */
        /* TODO use a dedicated tensor field by species */

        cs_field_t *fturbvisco
          = cs_field_by_name_try("anisotropic_turbulent_viscosity");

        if (fturbvisco != NULL) {
          cs_real_6_t  *visten_v = (cs_real_6_t *)fturbvisco->val;

          double long_diffus = 0, trans_diffus = 0;

          cs_tree_node_t *tn
            = cs_tree_node_get_child(tn_zl, "diffusion_coefficient");
          cs_gui_node_get_child_real(tn, "longitudinal", &long_diffus);
          cs_gui_node_get_child_real(tn, "transverse", &trans_diffus);

          const double diff = long_diffus - trans_diffus;

          for (cs_lnum_t icel = 0; icel < n_cells; icel++) {
            cs_lnum_t iel = cell_ids[icel];
            double norm = cs_math_3_square_norm(vel[iel]);
            double tmp = trans_diffus * norm;
            double denom = norm + 1.e-15;
            visten_v[iel][0] = tmp + diff * vel[iel][0] * vel[iel][0] / denom;
            visten_v[iel][1] = tmp + diff * vel[iel][1] * vel[iel][1] / denom;
            visten_v[iel][2] = tmp + diff * vel[iel][2] * vel[iel][2] / denom;
            visten_v[iel][3] =       diff * vel[iel][1] * vel[iel][0] / denom;
            visten_v[iel][4] =       diff * vel[iel][1] * vel[iel][2] / denom;
            visten_v[iel][5] =       diff * vel[iel][2] * vel[iel][0] / denom;
          }
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

      if (permeability == 0) {
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute GUI-defined head losses for a given volume zone.
 *
 * Head loss tensor coefficients for each cell are organized as follows:
 * cku11, cku22, cku33, cku12, cku13, cku23.
 *
 * \param[in]       zone       pointer to zone structure
 * \param[in]       cvara_vel  velocity values at the  previous time step
 * \param[in, out]  cku        head loss coefficients
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_head_losses(const cs_zone_t   *zone,
                   const cs_real_3_t *cvara_vel,
                   cs_real_t          cku[][6])
{
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply initial conditions based on GUI-defined settings.
 */
/*----------------------------------------------------------------------------*/

void cs_gui_initial_conditions(void)
{
  /* Coal combustion: the initialization of the model scalar are not given */

  int restart = cs_restart_present();
  int ccfth = 0;

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
#endif

  const int n_zones = cs_volume_zone_n_zones();

  const cs_real_3_t *restrict cell_cen =
    (const cs_real_3_t *restrict)cs_glob_mesh_quantities->cell_cen;

  for (int z_id = 0; z_id < n_zones; z_id++) {
    const cs_zone_t *z = cs_volume_zone_by_id(z_id);

    if (z->type & CS_VOLUME_ZONE_INITIALIZATION) {
      const cs_lnum_t n_cells = z->n_elts;
      const cs_lnum_t *cell_ids = z->elt_ids;

      char z_id_str[32];
      snprintf(z_id_str, 31, "%d", z_id);

      if (restart == 0) {

        cs_tree_node_t *tn_velocity
          = cs_tree_get_node(cs_glob_tree,
                             "thermophysical_models/velocity_pressure"
                             "/initialization/formula");
        tn_velocity = _add_zone_id_test_attribute(tn_velocity, z->id);
        const char *formula_uvw = cs_tree_node_get_value_str(tn_velocity);

        cs_field_t *c_vel = cs_field_by_name("velocity");

        if (formula_uvw != NULL) {
          cs_real_t *ini_vals = NULL;
          BFT_MALLOC(ini_vals, c_vel->dim * n_cells, cs_real_t);

          cs_meg_initialization(z->name,
                                n_cells,
                                cell_ids,
                                cell_cen,
                                "velocity",
                                ini_vals);

          cs_array_real_copy_subset(n_cells,
                                    c_vel->dim,
                                    cell_ids,
                                    CS_ARRAY_SUBSET_OUT,
                                    ini_vals,
                                    c_vel->val);

          BFT_FREE(ini_vals);
        }

        /* Void fraction initialization for VoF approach */
        if (cs_glob_vof_parameters->vof_model > 0) {

          cs_tree_node_t *tn
            = cs_tree_get_node(cs_glob_tree,
                               "thermophysical_models/"
                               "hgn_model/initialization");
          tn = cs_tree_find_node(tn, "formula");
          tn = _add_zone_id_test_attribute(tn, z->id);
          const char *formula = cs_tree_node_get_value_str(tn);

          cs_field_t *c = cs_field_by_name_try("void_fraction");

          if (c != NULL && formula != NULL) {
            cs_real_t *ini_vals = NULL;
            BFT_MALLOC(ini_vals, c->dim*n_cells, cs_real_t);

            cs_meg_initialization(z->name,
                                  n_cells,
                                  cell_ids,
                                  cell_cen,
                                  "void_fraction",
                                  ini_vals);

            cs_array_real_copy_subset(n_cells,
                                      c->dim,
                                      cell_ids,
                                      CS_ARRAY_SUBSET_OUT,
                                      ini_vals,
                                      c->val);
            BFT_FREE(ini_vals);
          }

        }

        /* Pressure initialization for groundwater model */
        if (cs_glob_physical_model_flag[CS_GROUNDWATER] > 0) {

          cs_tree_node_t *tn = _find_node_variable("hydraulic_head");
          tn = cs_tree_find_node(tn, "formula");
          tn = _add_zone_id_test_attribute(tn, z->id);

          const char *formula = cs_tree_node_get_value_str(tn);

          cs_field_t *c = cs_field_by_name_try("hydraulic_head");

          if (formula != NULL) {
            cs_real_t *ini_vals = NULL;
            BFT_MALLOC(ini_vals, c->dim*n_cells, cs_real_t);
            cs_meg_initialization(z->name,
                                  n_cells,
                                  cell_ids,
                                  cell_cen,
                                  "hydraulic_head",
                                  ini_vals);

            cs_array_real_copy_subset(n_cells,
                                      c->dim,
                                      cell_ids,
                                      CS_ARRAY_SUBSET_OUT,
                                      ini_vals,
                                      c->val);
            BFT_FREE(ini_vals);
          }

        }

        /* Turbulence variables initialization */
        const char *choice = _turbulence_initialization_choice(z_id_str);

        if (cs_gui_strcmp(choice, "formula")) {

          const char *formula_turb = NULL;
          cs_tree_node_t *tn_turb
            = cs_tree_get_node(cs_glob_tree,
                               "thermophysical_models/"
                               "turbulence/initialization");
          tn_turb = _add_zone_id_test_attribute(tn_turb, z->id);
          tn_turb = cs_tree_get_node(tn_turb, "formula");
          formula_turb = cs_tree_node_get_value_str(tn_turb);

          if (formula_turb != NULL) {
            const char *model = cs_gui_get_thermophysical_model("turbulence");
            if (model == NULL)
              break;
            if (cs_gui_strcmp(model, "off"))
              break;

            /* Set number of values to compute */
            int n_ini_vals = 0;
            if (   cs_gui_strcmp(model, "k-epsilon")
                || cs_gui_strcmp(model, "k-epsilon-PL")
                || cs_gui_strcmp(model, "k-omega-SST")) {
              n_ini_vals = 2;
            }
            else if (   cs_gui_strcmp(model, "Rij-epsilon")
                     || cs_gui_strcmp(model, "Rij-SSG")) {
              n_ini_vals = 7;
            }
            else if (cs_gui_strcmp(model, "Rij-EBRSM")) {
              n_ini_vals = 8;
            }
            else if (cs_gui_strcmp(model, "v2f-BL-v2/k")) {
              n_ini_vals = 4;
            }
            else if (cs_gui_strcmp(model, "Spalart-Allmaras")) {
              n_ini_vals = 1;
            }
            else
              bft_error(__FILE__, __LINE__, 0,
                        _("Invalid turbulence model: %s.\n"), model);


            cs_real_t *ini_vals = NULL;
            BFT_MALLOC(ini_vals, n_ini_vals*n_cells, cs_real_t);
            cs_meg_initialization(z->name,
                                  n_cells,
                                  cell_ids,
                                  cell_cen,
                                  "turbulence",
                                  ini_vals);


            if (   cs_gui_strcmp(model, "k-epsilon")
                || cs_gui_strcmp(model, "k-epsilon-PL")) {

              cs_field_t *c_k   = cs_field_by_name("k");
              cs_field_t *c_eps = cs_field_by_name("epsilon");

              for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
                cs_real_t *_vals = ini_vals + n_ini_vals * e_id;

                cs_lnum_t c_id = cell_ids[e_id];
                c_k->val[c_id]   = _vals[0];
                c_eps->val[c_id] = _vals[1];
              }
            }
            else if (   cs_gui_strcmp(model, "Rij-epsilon")
                     || cs_gui_strcmp(model, "Rij-SSG")) {

              cs_field_t *c_rij = cs_field_by_name_try("rij");
              cs_field_t *c_eps = cs_field_by_name("epsilon");

              for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
                cs_real_t *_vals = ini_vals + n_ini_vals * e_id;

                cs_lnum_t c_id = cell_ids[e_id];
                cs_real_t *_rij = c_rij->val + 6*c_id;
                for (int drij = 0; drij < 6; drij++)
                  _rij[drij] = _vals[drij];

                c_eps->val[c_id] = _vals[6];
              }
            }
            else if (cs_gui_strcmp(model, "Rij-EBRSM")) {
              cs_field_t *c_rij = cs_field_by_name_try("rij");
              cs_field_t *c_eps = cs_field_by_name("epsilon");
              cs_field_t *c_alp = cs_field_by_name("alpha");

              for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
                cs_real_t *_vals = ini_vals + n_ini_vals * e_id;

                cs_lnum_t c_id = cell_ids[e_id];

                cs_real_t *_rij = c_rij->val + 6*c_id;
                for (int drij = 0; drij < 6; drij++)
                  _rij[drij] = _vals[drij];

                c_eps->val[c_id] = _vals[6];
                c_alp->val[c_id] = _vals[7];
              }
            }
            else if (cs_gui_strcmp(model, "v2f-BL-v2/k")) {

              cs_field_t *c_k   = cs_field_by_name("k");
              cs_field_t *c_eps = cs_field_by_name("epsilon");
              cs_field_t *c_phi = cs_field_by_name("phi");
              cs_field_t *c_alp = cs_field_by_name("alpha");

              for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
                cs_real_t *_vals = ini_vals + n_ini_vals * e_id;

                cs_lnum_t c_id = cell_ids[e_id];

                c_k->val[c_id]   = _vals[0];
                c_eps->val[c_id] = _vals[1];
                c_phi->val[c_id] = _vals[2];
                c_alp->val[c_id] = _vals[3];
              }
            }
            else if (cs_gui_strcmp(model, "k-omega-SST")) {

              cs_field_t *c_k   = cs_field_by_name("k");
              cs_field_t *c_ome = cs_field_by_name("omega");

              for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
                cs_real_t *_vals = ini_vals + n_ini_vals * e_id;

                cs_lnum_t c_id = cell_ids[e_id];

                c_k->val[c_id]   = _vals[0];
                c_ome->val[c_id] = _vals[1];
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

        /* Thermal scalar initialization */
        if (cs_gui_thermal_model_code() > 0) {

          const char *formula_sca    = NULL;
          cs_tree_node_t *tn_sca
            = cs_tree_get_node(cs_glob_tree,
                               "thermophysical_models/"
                               "thermal_scalar/variable/formula");
          tn_sca = _add_zone_id_test_attribute(tn_sca, z->id);
          formula_sca = cs_tree_node_get_value_str(tn_sca);

          /* For non-specific physics defined with the GUI,
             the thermal variable can only be temperature or enthalpy
             (as the thermal model is on) */

          cs_field_t *c = cs_thermal_model_field();

          assert(c != NULL);

          if (formula_sca != NULL) {
            cs_real_t *ini_vals = NULL;
            BFT_MALLOC(ini_vals, n_cells, cs_real_t);
            cs_meg_initialization(z->name,
                                  n_cells,
                                  cell_ids,
                                  cell_cen,
                                  "thermal",
                                  ini_vals);

            cs_array_real_copy_subset(n_cells,
                                      c->dim,
                                      cell_ids,
                                      CS_ARRAY_SUBSET_OUT,
                                      ini_vals,
                                      c->val);
            BFT_FREE(ini_vals);
          }
          /* If no formula was provided, the previous field values are
             kept (allowing mode-specific automatic initialization). */
        }

        /* User Scalars initialization */
        int n_fields = cs_field_n_fields();

        for (int f_id = 0; f_id < n_fields; f_id++) {

          const cs_field_t  *f = cs_field_by_id(f_id);

          if (   f->type & CS_FIELD_USER
              && f->location_id == CS_MESH_LOCATION_CELLS) {

            const char *formula_sca    = NULL;

            cs_tree_node_t *tn_sca = NULL;
            tn_sca = cs_tree_get_node(cs_glob_tree,
                                      "additional_scalars/variable");
            tn_sca = cs_tree_node_get_sibling_with_tag(tn_sca, "name", f->name);
            tn_sca = cs_tree_get_node(tn_sca, "formula");
            tn_sca = _add_zone_id_test_attribute(tn_sca, z->id);
            formula_sca = cs_tree_node_get_value_str(tn_sca);

            if (formula_sca != NULL) {
              cs_real_t *ini_vals = NULL;
              BFT_MALLOC(ini_vals, n_cells*f->dim, cs_real_t);
              cs_meg_initialization(z->name,
                                    n_cells,
                                    cell_ids,
                                    cell_cen,
                                    f->name,
                                    ini_vals);

              cs_array_real_copy_subset(n_cells,
                                        f->dim,
                                        cell_ids,
                                        CS_ARRAY_SUBSET_OUT,
                                        ini_vals,
                                        f->val);
              BFT_FREE(ini_vals);
            }
          }
        }

        /* Meteo Scalars initialization */
        if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] > -1) {

          cs_tree_node_t *tn_m0
            = cs_tree_get_node(cs_glob_tree,
                               "thermophysical_models/atmospheric_flows");

          const char *name       = NULL;
          const char *formula_meteo  = NULL;

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
              cs_real_t *ini_vals = NULL;
              BFT_MALLOC(ini_vals, c->dim*n_cells, cs_real_t);
              cs_meg_initialization(z->name,
                                    n_cells,
                                    cell_ids,
                                    cell_cen,
                                    c->name,
                                    ini_vals);

              cs_array_real_copy_subset(n_cells,
                                        c->dim,
                                        cell_ids,
                                        CS_ARRAY_SUBSET_OUT,
                                        ini_vals,
                                        c->val);
              BFT_FREE(ini_vals);
            }

            /* else:
               - Do not overwrite possible prior automatic initialization;
               - All fields are otherwise set to 0 by default when
                 allocated, so no need to assign a value here. */
          }

        }

        /* Combustion Scalars initialization */
        if (cs_glob_physical_model_flag[CS_COMBUSTION_3PT] > -1) {

          cs_tree_node_t *tn_gas
            = cs_tree_get_node(cs_glob_tree,
                               "thermophysical_models/gas_combustion");

          const char *name       = NULL;
          const char *formula_comb  = NULL;

          int size = cs_tree_get_sub_node_count_simple(tn_gas, "variable");

          for (int j = 0; j < size; j++) {
            cs_tree_node_t *tn_combustion = cs_tree_get_node(tn_gas,
                                                             "variable");
            for (int i = 1;
                 tn_combustion != NULL && i < j + 1;
                 i++) {
             tn_combustion = cs_tree_node_get_next_of_name(tn_combustion);
            }

            cs_tree_node_t *tn_combustion2 = tn_combustion;
            tn_combustion = cs_tree_get_node(tn_combustion, "name");
            name = cs_tree_node_get_value_str(tn_combustion);

            tn_combustion2 = cs_tree_get_node(tn_combustion2, "formula");
            tn_combustion2 = _add_zone_id_test_attribute(tn_combustion2, z->id);

            cs_field_t *c_comb = cs_field_by_name_try(name);

            formula_comb = cs_tree_node_get_value_str(tn_combustion2);

            if (formula_comb != NULL) {
              cs_real_t *ini_vals = NULL;
              BFT_MALLOC(ini_vals, c_comb->dim*n_cells, cs_real_t);
              cs_meg_initialization(z->name,
                                    n_cells,
                                    cell_ids,
                                    cell_cen,
                                    c_comb->name,
                                    ini_vals);

            cs_array_real_copy_subset(n_cells,
                                      c_comb->dim,
                                      cell_ids,
                                      CS_ARRAY_SUBSET_OUT,
                                      ini_vals,
                                      c_comb->val);
            BFT_FREE(ini_vals);
            }
          }

        }

        if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1) {
          const char *formula        = NULL;
          const char *buff           = NULL;
          const char *name[] = {"pressure", "temperature", "total_energy",
                                "density"};

          ccfth = 10000;
          for (int j = 0; j < 4; j++) {

            cs_tree_node_t *tn = NULL;
            if (j < 3) {
              tn = cs_tree_find_node(cs_glob_tree, "variable");
              while (tn != NULL) {
                const char *name_tn
                  = cs_tree_node_get_child_value_str(tn, "name");
                if (cs_gui_strcmp(name_tn, name[j]))
                  break;
                else
                  tn = cs_tree_find_node_next(cs_glob_tree, tn, "variable");
              }
            }
            else {
              tn = cs_tree_find_node(cs_glob_tree, "property");
              while (tn != NULL) {
                const char *name_tn
                  = cs_tree_node_get_child_value_str(tn, "name");
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

              if (formula != NULL && c != NULL) {
                cs_real_t *ini_vals = NULL;
                BFT_MALLOC(ini_vals, c->dim*n_cells, cs_real_t);

                cs_meg_initialization(z->name,
                                      n_cells,
                                      cell_ids,
                                      cell_cen,
                                      c->name,
                                      ini_vals);

                cs_array_real_copy_subset(n_cells,
                                          c->dim,
                                          cell_ids,
                                          CS_ARRAY_SUBSET_OUT,
                                          ini_vals,
                                          c->val);
                BFT_FREE(ini_vals);
              }
            }

          }
          cs_cf_model_t *cf_model = cs_get_glob_cf_model();
          cf_model->ithvar = ccfth;
        }

      } /* End for test on restart */
    }
  } /* zones+1 */
}

/*-----------------------------------------------------------------------------
 * Selection of linear solvers.
 *----------------------------------------------------------------------------*/

void
cs_gui_linear_solvers(void)
{
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

      cs_tree_node_t *tn_v = _find_node_variable(f->name);

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
      else if (cs_gui_strcmp(algo_choice, "gcr"))
        sles_it_type = CS_SLES_GCR;
      else if (cs_gui_strcmp(algo_choice, "gauss_seidel"))
        sles_it_type = CS_SLES_P_GAUSS_SEIDEL;
      else if (cs_gui_strcmp(algo_choice, "symmetric_gauss_seidel"))
        sles_it_type = CS_SLES_P_SYM_GAUSS_SEIDEL;
      else if (cs_gui_strcmp(algo_choice, "PCR3"))
        sles_it_type = CS_SLES_PCR3;

      /* If preconditioning choice is set, we need to define a solver
         (otherwise, we would need to pass at least this partial
         information to cs_sles_default) */

      if (sles_it_type >= CS_SLES_N_IT_TYPES) {
        if (   precond_choice != NULL
            && !cs_gui_strcmp(precond_choice, "automatic")) {

          bool symmetric = false;
          cs_equation_param_t *eqp = cs_field_get_equation_param(f);
          if (eqp != NULL) {
            if (eqp->iconv == 0)
              symmetric = true;
          }

          if (symmetric)
            sles_it_type = CS_SLES_FCG;
          else {
            int coupling_id
              = cs_field_get_key_int(f, cs_field_key_id("coupling_entity"));
            if (coupling_id < 0 && cs_gui_strcmp(precond_choice, "none"))
              sles_it_type = CS_SLES_P_SYM_GAUSS_SEIDEL;
            else
              sles_it_type = CS_SLES_BICGSTAB;
          }

        }
      }

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
        else if (cs_gui_strcmp(precond_choice, "multigrid_k_cycle_hpc")) {
          pc_multigrid = true;
          mg_type = CS_MULTIGRID_K_CYCLE_HPC;
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

/*----------------------------------------------------------------------------
 * User momentum source terms.
 *
 * parameters:
 *   vel      <--  fluid velocity
 *   tsexp    -->  explicit source terms
 *   tsimp    -->  implicit source terms
 *----------------------------------------------------------------------------*/

void
cs_gui_momentum_source_terms(const cs_real_3_t  *restrict vel,
                             cs_real_3_t        *restrict tsexp,
                             cs_real_33_t       *restrict tsimp)
{
  const cs_real_t *restrict cell_f_vol = cs_glob_mesh_quantities->cell_f_vol;

  double Su, Sv, Sw;
  double dSudu, dSudv, dSudw;
  double dSvdu, dSvdv, dSvdw;
  double dSwdu, dSwdv, dSwdw;

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
#endif

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
        cs_real_t *st_vals = NULL;
        BFT_MALLOC(st_vals, 12*n_cells, cs_real_t);

        const cs_real_3_t *restrict cell_cen =
          (const cs_real_3_t *restrict)cs_glob_mesh_quantities->cell_cen;

        cs_meg_source_terms(z->name,
                            z->n_elts,
                            z->elt_ids,
                            cell_cen,
                            "momentum",
                            "momentum_source_term",
                            st_vals);

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

        BFT_FREE(st_vals);
      }
    }
  }
}

/*-----------------------------------------------------------------------------
 * Define global numerical options.
 *----------------------------------------------------------------------------*/

void
cs_gui_numerical_options(void)
{
  cs_velocity_pressure_param_t *vp_param = cs_get_glob_velocity_pressure_param();
  cs_velocity_pressure_model_t *vp_model = cs_get_glob_velocity_pressure_model();
  cs_space_disc_t *space_disc = cs_get_glob_space_disc();

  const char *choice = NULL;

  cs_tree_node_t *tn_n = cs_tree_get_node(cs_glob_tree, "numerical_parameters");

  int _imrgra = -1;

  cs_ext_neighborhood_type_t enh_type = cs_ext_neighborhood_get_type();
  bool enh_boundary = cs_ext_neighborhood_get_boundary_complete();

  choice = cs_tree_node_get_tag(cs_tree_get_node(tn_n,
                                                 "gradient_reconstruction"),
                                "choice");
  if (cs_gui_strcmp(choice, "green_iter"))
    _imrgra = 0;
  else if (cs_gui_strcmp(choice, "lsq"))
    _imrgra = 1;
  else if (cs_gui_strcmp(choice, "green_lsq"))
    _imrgra = 4;
  else if (cs_gui_strcmp(choice, "green_vtx"))
    _imrgra = 7;

  if (_imrgra != 0) {
    int _imrgra_add = 0;
    choice = cs_tree_node_get_tag(cs_tree_get_node(tn_n,
                                                   "extended_neighborhood"),
                                  "choice");
    if (cs_gui_strcmp(choice, "none")) {
      enh_type = CS_EXT_NEIGHBORHOOD_NONE;
    }
    else if (cs_gui_strcmp(choice, "boundary")) {
      enh_type = CS_EXT_NEIGHBORHOOD_NONE;
      enh_boundary = true;
    }
    else if (cs_gui_strcmp(choice, "complete")) {
      enh_type = CS_EXT_NEIGHBORHOOD_COMPLETE;
      _imrgra_add = 1;
    }
    else if (cs_gui_strcmp(choice, "optimized")) {
      enh_type = CS_EXT_NEIGHBORHOOD_OPTIMIZED;
      _imrgra_add = 2;
    }
    else if (cs_gui_strcmp(choice, "optimized_with_boundary")) {
      enh_type = CS_EXT_NEIGHBORHOOD_OPTIMIZED;
      _imrgra_add = 2;
      enh_boundary = true;
    }
    else if (cs_gui_strcmp(choice, "cell_center_opposite")) {
      enh_type = CS_EXT_NEIGHBORHOOD_CELL_CENTER_OPPOSITE;
      _imrgra_add = 2;
    }
    else if (cs_gui_strcmp(choice, "cell_center_opposite_with_boundary")) {
      enh_type = CS_EXT_NEIGHBORHOOD_CELL_CENTER_OPPOSITE;
      _imrgra_add = 2;
      enh_boundary = true;
    }
    else if (cs_gui_strcmp(choice, "non_ortho_max")) {
      enh_type = CS_EXT_NEIGHBORHOOD_NON_ORTHO_MAX;
      _imrgra_add = 2;
    }

    if (_imrgra_add > 0) {
      if (_imrgra < 0)
        _imrgra = space_disc->imrgra;
      _imrgra += _imrgra_add;
    }
  }

  cs_ext_neighborhood_set_type(enh_type);
  cs_ext_neighborhood_set_boundary_complete(enh_boundary);

  if (_imrgra > -1)
    space_disc->imrgra = _imrgra;

  int _idilat = -1;

  choice = cs_tree_node_get_tag(cs_tree_get_node(tn_n, "algo_density_variation"),
                                "choice");
  if (cs_gui_strcmp(choice, "boussi"))
    _idilat = 0;
  else if (cs_gui_strcmp(choice, "dilat_std"))
    _idilat = 1;
  else if (cs_gui_strcmp(choice, "dilat_unstd"))
    _idilat = 2;
  else if (cs_gui_strcmp(choice, "low_mach"))
    _idilat = 3;
  else if (cs_gui_strcmp(choice, "algo_fire"))
    _idilat = 4;

  if (_idilat > -1)
    vp_model->idilat = _idilat;

  _numerical_int_parameters("gradient_transposed", &(vp_model->ivisse));
  _numerical_int_parameters("velocity_pressure_coupling", &(vp_param->ipucou));
  _numerical_int_parameters("piso_sweep_number", &(vp_param->nterup));

  if (cs_glob_time_step_options->idtvar > -1) {
    double _relaxp = -1.;
    _numerical_double_parameters("pressure_relaxation", &_relaxp);
    if (_relaxp > -1.0 && CS_F_(p) != NULL) {
      cs_equation_param_t *eqp = cs_field_get_equation_param(CS_F_(p));
      eqp->relaxv = _relaxp;
    }
  }

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--ivisse = %i\n", vp_model->ivisse);
  bft_printf("--ipucou = %i\n", vp_model->ipucou);
  bft_printf("--imrgra = %i\n", *imrgra);
  bft_printf("--nterup = %i\n", vp_param->nterup);
  bft_printf("--relaxp = %f\n", _relaxp);
#endif
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
    int def_rank_step;
    cs_file_get_default_comm(&def_rank_step, NULL, NULL);
    size_t def_block_size = cs_parall_get_min_coll_buf_size();
    if (rank_step < 1)
      rank_step = def_rank_step;
    if (block_size < 0)
      block_size = def_block_size;
    cs_file_set_default_comm(rank_step, cs_glob_mpi_comm);
    cs_parall_set_min_coll_buf_size(block_size);
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

/*-----------------------------------------------------------------------------
 * Set MPI related algorithm options
 *----------------------------------------------------------------------------*/

void
cs_gui_mpi_algorithms(void)
{
  cs_all_to_all_type_t a = CS_ALL_TO_ALL_MPI_DEFAULT;

  cs_tree_node_t *tn_cm
    = cs_tree_get_node(cs_glob_tree, "calculation_management");

  /* Partitioning type */
  const char  *all_to_all_name
    = cs_tree_node_get_child_value_str(tn_cm, "all_to_all");

  if (all_to_all_name != NULL) {
    if (!strcmp(all_to_all_name, "default"))
      a = CS_ALL_TO_ALL_MPI_DEFAULT;
    else if (!strcmp(all_to_all_name, "crystal router"))
      a = CS_ALL_TO_ALL_CRYSTAL_ROUTER;
    cs_all_to_all_set_type(a);
  }
}

/*----------------------------------------------------------------------------
 * Treatment of physical constants (gravity and Coriolis).
 *----------------------------------------------------------------------------*/

void
cs_gui_physical_constants(void)
{
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

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--gx = %f \n", phys_cst->gravity[0]);
  bft_printf("--gy = %f \n", phys_cst->gravity[1]);
  bft_printf("--gz = %f \n", phys_cst->gravity[2]);
  bft_printf("--icorio = %i \n", phy_cst->icorio);
#endif
}

/*----------------------------------------------------------------------------
 * Treatment of fluid physical properties
 * Initialize reference pressure and temperature if present
 *----------------------------------------------------------------------------*/

void
cs_gui_physical_properties(void)
{
  int choice;
  const char *material = NULL;

  const int itherm = cs_glob_thermal_model->itherm;

  cs_fluid_properties_t *phys_pp = cs_get_glob_fluid_properties();
  cs_gui_fluid_properties_value("reference_pressure", &(phys_pp->p0));

  /* Variable rho and viscl */
  if (_properties_choice_id("density", &choice))
    phys_pp->irovar = choice;

  if (_properties_choice_id("molecular_viscosity", &choice))
    phys_pp->ivivar = choice;
  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1)
    if (_properties_choice_id("molecular_viscosity", &choice))
      phys_pp->ivivar = choice;

  /* Read T0 in each case for user */
  cs_gui_fluid_properties_value("reference_temperature", &(phys_pp->t0));

  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1)
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

      const cs_temperature_scale_t itpscl = cs_glob_thermal_model->itpscl;
      const char *method = _thermal_table_choice("method");

      if (   strcmp(method, "CoolProp") == 0
          && itherm > CS_THERMAL_MODEL_NONE) {
        if (cs_physical_properties_get_coolprop_backend() == NULL) {
          int n_th_laws = _count_thermal_laws();
          if (n_th_laws > 0)
            cs_physical_properties_set_coolprop_backend("TTSE&HEOS");
        }
      }

      cs_thermal_table_set(material,
                           method,
                           _thermal_table_option("reference"),
                           thermal_plane,
                           itpscl);
    }
  }

  cs_vof_parameters_t *vof_param = cs_get_glob_vof_parameters();

  if (_thermal_table_needed("density") == 0) {
    cs_gui_properties_value("density", &phys_pp->ro0);
    if (vof_param->vof_model & CS_VOF_ENABLED) {
      cs_gui_properties_value_by_fluid_id(0, "density", &vof_param->rho1);
      cs_gui_properties_value_by_fluid_id(1, "density", &vof_param->rho2);
    }
  }
  else {
    cs_phys_prop_compute(CS_PHYS_PROP_DENSITY,
                         1,
                         0,
                         0,
                         &phys_pp->p0,
                         &phys_pp->t0,
                         &phys_pp->ro0);
  }

  const char *mv_name = "molecular_viscosity";
  if (_thermal_table_needed(mv_name) == 0) {
    cs_gui_properties_value(mv_name, &phys_pp->viscl0);
    if (vof_param->vof_model & CS_VOF_ENABLED) {
      cs_gui_properties_value_by_fluid_id(0, mv_name, &vof_param->mu1);
      cs_gui_properties_value_by_fluid_id(1, mv_name, &vof_param->mu2);
    }
  }
  else {
    cs_phys_prop_compute(CS_PHYS_PROP_DYNAMIC_VISCOSITY,
                         1,
                         0,
                         0,
                         &phys_pp->p0,
                         &phys_pp->t0,
                         &phys_pp->viscl0);
  }

  if (vof_param->vof_model & CS_VOF_ENABLED) {
    const char *st_name = "surface_tension";
    cs_gui_properties_value(st_name, &vof_param->sigma_s);
  }

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

  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1) {
    cs_gui_properties_value("volume_viscosity", &phys_pp->viscv0);
    double visls_0 = -1;
    cs_gui_properties_value("thermal_conductivity", &visls_0);
    cs_field_set_key_double(cs_field_by_name("temperature"),
                            cs_field_key_id("diffusivity_ref"),
                            visls_0);
  }

  int n_zones_pp
    = cs_volume_zone_n_type_zones(CS_VOLUME_ZONE_PHYSICAL_PROPERTIES);
  if (n_zones_pp > 1) {
    phys_pp->irovar = 1;
    phys_pp->ivivar = 1;

    cs_field_t *tf = cs_thermal_model_field();
    if (tf != NULL) {
      phys_pp->icp = 1;
      int k = cs_field_key_id("diffusivity_id");
      int cond_diff_id = cs_field_get_key_int(tf, k);
      if (cond_diff_id < 0)
        cs_field_set_key_int(tf, k, 0);
    }
  }

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
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
  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1) {
    bft_printf("--viscv0 = %g \n", *viscv0);
    bft_printf("--xmasmr = %f \n", cs_glob_fluid_properties->xmasmr);
  }
#endif
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
        else if (cs_gui_strcmp(mdl, "integral"))
          cs_glob_porous_model = 3;
      }
    }
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
 * Get value of property markup for fluid of given id
 *
 * parameters:
 *   fluid_id           <--   fluid index
 *   property_name      <--   name of the property
 *   value              -->   new initial value of the property
 *----------------------------------------------------------------------------*/

void
cs_gui_properties_value_by_fluid_id(const int    fluid_id,
                                    const char  *property_name,
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

  char *label = NULL;

  const char *prefix = "value_";
  size_t len = strlen(prefix) + 1;
  BFT_MALLOC(label, len+1, char);

  sprintf(label, "%s%1i", prefix, fluid_id);

  tn = cs_tree_get_node(tn, label);
  cs_gui_node_get_real(tn, value);

  BFT_FREE(label);
}

/*----------------------------------------------------------------------------
 * Read minimum / maximum values (used in clipping) and turbulent flux model
 * for additional user or model variables.
 *
 * Also read reference dynamic and user scalar viscosity
 *----------------------------------------------------------------------------*/

void
cs_gui_scalar_model_settings(void)
{
#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
#endif

  const cs_turb_model_t  *turb_model = cs_get_glob_turb_model();
  assert(turb_model != NULL);

  const int kscmin = cs_field_key_id("min_scalar_clipping");
  const int kscmax = cs_field_key_id("max_scalar_clipping");

  /* Specific physics: the min max of the model scalar are not given */
  const int keysca = cs_field_key_id("scalar_id");
  const int kturt  = cs_field_key_id("turbulent_flux_model");

  /* Values used by optional beta-limiters: do not need to be the same
     as clipping values, but this seesm a good initial setting
     (an independent value in the GUI would be better). */
  const int k_beta_max = cs_field_key_id("max_scalar");
  const int k_beta_min = cs_field_key_id("min_scalar");

  for (int f_id = 0; f_id < cs_field_n_fields(); f_id++) {
    cs_field_t  *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) { /* variable ? */
      int i = cs_field_get_key_int(f, keysca) - 1;
      if (i > -1) { /* additional user or model variable ? */

        double scal_min = cs_field_get_key_double(f, kscmin);
        double scal_max = cs_field_get_key_double(f, kscmax);

        cs_tree_node_t *tn_v = _find_node_variable(f->name);
        if (tn_v != NULL) { /* variable is in xml ? */
          cs_gui_node_get_child_real(tn_v, "min_value", &scal_min);
          cs_gui_node_get_child_real(tn_v, "max_value", &scal_max);
          cs_field_set_key_double(f, kscmin, scal_min);
          cs_field_set_key_double(f, kscmax, scal_max);

          cs_field_set_key_double(f, k_beta_max, scal_max);
          cs_field_set_key_double(f, k_beta_min, scal_min);

#if _XML_DEBUG_
          bft_printf("--min_scalar_clipping[%i] = %f\n", i, scal_min);
          bft_printf("--max_scalar_clipping[%i] = %f\n", i, scal_max);
#endif

          if (turb_model->order == CS_TURB_SECOND_ORDER) {
            int turb_mdl;
            _variable_turbulent_flux_model(tn_v, &turb_mdl);
            cs_field_set_key_int(f, kturt, turb_mdl);
#if _XML_DEBUG_
            bft_printf("--turb_model[%i] = %d\n", i, turb_mdl);
#endif
          }

        }
      }
    }
  }

  /* Also read/set diffusivity values */

  _read_diffusivity();
}

/*----------------------------------------------------------------------------
 * Thermal model.
 *----------------------------------------------------------------------------*/

void
cs_gui_thermal_model(void)
{
  cs_thermal_model_t *thermal = cs_get_glob_thermal_model();

  switch(cs_gui_thermal_model_code()) {
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
 * Get thermal scalar model.
 *
 * return:
 *   value of itherm*10 + (temperature variant flag), or -1 if not defined
 *----------------------------------------------------------------------------*/

int
cs_gui_thermal_model_code(void)
{
  int   test = -1;

  const char *model = cs_gui_get_thermophysical_model("thermal_scalar");

  if (model == NULL)
    return test;

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

      if (f != NULL) {
        m_f_id[j] = f->id;
        m_c_id[j] = idim;
      }

      else
        bft_error(__FILE__, __LINE__, 0,
                  _("Time moment \"%s\"\n"
                    "requires undefined field \"%s\"."),
                  m_name, f_name);

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
 * Turbulence model
 *----------------------------------------------------------------------------*/

void
cs_gui_turb_model(void)
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
    turb_mdl->iturb = CS_TURB_NONE;
  else if (cs_gui_strcmp(model, "mixing_length")) {
    turb_mdl->iturb = CS_TURB_MIXING_LENGTH;
    cs_gui_node_get_child_real(tn_t, "mixing_length_scale", &(rans_mdl->xlomlg));
  }
  else if (cs_gui_strcmp(model, "k-epsilon")) {
    turb_mdl->iturb = CS_TURB_K_EPSILON;
    cs_gui_node_get_child_int(tn_t, "wall_function", &iwallf);
    cs_gui_node_get_child_status_int(tn_t, "gravity_terms", &(rans_mdl->igrake));
  }
  else if (cs_gui_strcmp(model, "k-epsilon-PL")) {
    turb_mdl->iturb = CS_TURB_K_EPSILON_LIN_PROD;
    cs_gui_node_get_child_int(tn_t, "wall_function", &iwallf);
    cs_gui_node_get_child_status_int(tn_t, "gravity_terms", &(rans_mdl->igrake));
  }
  else if (cs_gui_strcmp(model, "Rij-epsilon")) {
    turb_mdl->iturb = CS_TURB_RIJ_EPSILON_LRR;
    cs_gui_node_get_child_int(tn_t, "wall_function", &iwallf);
    cs_gui_node_get_child_status_int(tn_t, "gravity_terms", &(rans_mdl->igrari));
  }
  else if (cs_gui_strcmp(model, "Rij-SSG")) {
    turb_mdl->iturb = CS_TURB_RIJ_EPSILON_SSG;
    cs_gui_node_get_child_int(tn_t, "wall_function", &iwallf);
    cs_gui_node_get_child_status_int(tn_t, "gravity_terms", &(rans_mdl->igrari));
    cs_gui_node_get_child_status_int(tn_t, "coupled_rij", &(rans_mdl->irijco));
  }
  else if (cs_gui_strcmp(model, "Rij-EBRSM")) {
    turb_mdl->iturb = CS_TURB_RIJ_EPSILON_EBRSM;
    cs_gui_node_get_child_int(tn_t, "wall_function", &iwallf);
    cs_gui_node_get_child_status_int(tn_t, "gravity_terms", &(rans_mdl->igrari));
    cs_gui_node_get_child_status_int(tn_t, "coupled_rij", &(rans_mdl->irijco));
  }
  else if (cs_gui_strcmp(model, "LES_Smagorinsky")) {
    turb_mdl->iturb = CS_TURB_LES_SMAGO_CONST;
  }
  else if (cs_gui_strcmp(model, "LES_dynamique")) {
    turb_mdl->iturb = CS_TURB_LES_SMAGO_DYN;
  }
  else if (cs_gui_strcmp(model, "LES_WALE")) {
    turb_mdl->iturb = CS_TURB_LES_WALE;
  }
  else if (cs_gui_strcmp(model, "v2f-phi")) {
    turb_mdl->iturb = CS_TURB_V2F_PHI;
    cs_gui_node_get_child_int(tn_t, "wall_function", &iwallf);
    cs_gui_node_get_child_status_int(tn_t, "gravity_terms", &(rans_mdl->igrake));
  }
  else if (cs_gui_strcmp(model, "v2f-BL-v2/k")) {
    turb_mdl->iturb = CS_TURB_V2F_BL_V2K;
    cs_gui_node_get_child_int(tn_t, "wall_function", &iwallf);
    cs_gui_node_get_child_status_int(tn_t, "gravity_terms", &(rans_mdl->igrake));
  }
  else if (cs_gui_strcmp(model, "k-omega-SST")) {
    turb_mdl->iturb = CS_TURB_K_OMEGA;
    cs_gui_node_get_child_int(tn_t, "wall_function", &iwallf);
    cs_gui_node_get_child_status_int(tn_t, "gravity_terms", &(rans_mdl->igrake));
  }
  else if (cs_gui_strcmp(model, "Spalart-Allmaras")) {
    turb_mdl->iturb = CS_TURB_SPALART_ALLMARAS;
  }
  else
    bft_error(__FILE__, __LINE__, 0,
        _("Invalid turbulence model: %s.\n"), model);

  if (iwallf != -1) {
    cs_wall_functions_t *wall_fnt = cs_get_glob_wall_functions();
    wall_fnt->iwallf = (cs_wall_f_type_t)iwallf;
  }

  if (   turb_mdl->iturb >= CS_TURB_RIJ_EPSILON_LRR
      && turb_mdl->iturb <= CS_TURB_RIJ_EPSILON_EBRSM) {
    const char *s
      = cs_tree_node_get_child_value_str(tn_t, "turbulent_diffusion_model");
    if (s != NULL) {
      if (cs_gui_strcmp(s, "shir"))
        rans_mdl->idirsm = 0;
      else if (cs_gui_strcmp(s, "daly_harlow"))
        rans_mdl->idirsm = 1;
    }
  }

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--model: %s\n", model);
  bft_printf("--iturb = %i\n", turb_mdl->iturb);
  bft_printf("--igrake = %i\n", rans_mdl->igrake);
  bft_printf("--igrari = %i\n", rans_mdl->igrari);
  bft_printf("--iwallf = %i\n", wall_fnt->iwallf);
  bft_printf("--xlomlg = %f\n", rans_mdl->xlomlg);
  bft_printf("--idirsm = %f\n", rans_mdl->idirsm);
#endif
}

/*----------------------------------------------------------------------------
 * Define reference length and reference velocity for initialization of the
 * turbulence variables
 *----------------------------------------------------------------------------*/

void
cs_gui_turb_ref_values(void)
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
    cs_mesh_location_type_t _loc
      = _cs_mesh_location_type_from_str(location_name);

    cs_parameters_add_property(name, array_dim, _loc);
  }
}

/*----------------------------------------------------------------------------
 * Define user calculator functions through the GUI
 *----------------------------------------------------------------------------*/

void
cs_gui_calculator_functions(void)
{
  const char path_s[] = "user_functions/calculator/function";
  cs_tree_node_t *tn_s = cs_tree_get_node(cs_glob_tree, path_s);

  for (cs_tree_node_t *tn = tn_s;
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn)) {

    const char *name = cs_gui_node_get_tag(tn, "name");

    int dim = 1;
    cs_tree_node_t *dtn = cs_tree_get_node(tn, "dimension");
    cs_gui_node_get_int(dtn, &dim);

    const char *location_name = cs_gui_node_get_tag(tn, "support");
    cs_mesh_location_type_t _loc
      = _cs_mesh_location_type_from_str(location_name);

    /* Define the cs_function_t based on MEG function*/
    cs_tree_node_t *n_f = cs_tree_get_node(tn, "formula");
    if (n_f != NULL) {
      cs_meg_xdef_input_t *_input
        = cs_meg_xdef_wrapper_add_input(CS_MEG_CALCULATOR_FUNC,
                                        -1,
                                        dim,
                                        name,
                                        NULL);

      cs_function_t *f
        = cs_function_define_by_analytic_func(name,
                                              _loc,
                                              dim,
                                              true,
                                              cs_meg_xdef_wrapper,
                                              _input);
      cs_function_set_label(f, name);
      f->log      = 1;
    }
  }
}

/*----------------------------------------------------------------------------
 * Define volume and boundary zones through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_zones(void)
{
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

    /* Check for solid zones */

    if (_zone_is_type(tn, "solid"))
      type_flag = type_flag | CS_VOLUME_ZONE_SOLID;

    /* Check if zone is used to define variable physical properties */
    if (_zone_is_type(tn, "physical_properties"))
      type_flag = type_flag | CS_VOLUME_ZONE_PHYSICAL_PROPERTIES;

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
                  0,   /* mode */
                  inlet_axis_coords,
                  outlet_axis_coords,
                  fan_radius,
                  blades_radius,
                  hub_radius,
                  pressure_curve_coeffs,
                  axial_torque);
  }
}

/*----------------------------------------------é------------------------------
 * Define error estimator through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_error_estimator(void)
{
  const char *name[] = {"est_error_pre_2",
                        "est_error_der_2",
                        "est_error_cor_2",
                        "est_error_tot_2"};

  const char *label[] = {"EsPre2",
                         "EsDer2",
                         "EsCor",
                         "EsTot2"};

  const char *node_name[] = {"Correction/model",
                             "Drift/model",
                             "Prediction/model",
                             "Total/model"};

  cs_tree_node_t *tn_ee
    = cs_tree_get_node(cs_glob_tree, "analysis_control/error_estimator");

  for (int i = 0; i < 4; i++) {

    cs_tree_node_t *tn = cs_tree_get_node(tn_ee, node_name[i]);

    const char *result = cs_tree_node_get_value_str(tn);

    if (   cs_gui_strcmp(result, "1")
        || cs_gui_strcmp(result, "2")) {

      int field_type = CS_FIELD_INTENSIVE | CS_FIELD_POSTPROCESS;

      cs_field_t *f = cs_field_create(name[i],
                                      field_type,
                                      CS_MESH_LOCATION_CELLS,
                                      1,
                                      false);

      cs_field_set_key_int(f, cs_field_key_id("log"), 1);
      cs_field_set_key_int(f, cs_field_key_id("post_vis"), 1);

      cs_field_set_key_str(f, cs_field_key_id("label"), label[i]);

    }

  }
}

/*----------------------------------------------------------------------------
 * Define internal coupling through the GUI.
 *----------------------------------------------------------------------------*/

void
cs_gui_internal_coupling(void)
{
  int n_zones = cs_volume_zone_n_zones();

  int n_solid_zones = 0;
  for (int i = 0; i < n_zones; i++) {
    const cs_zone_t  *z = cs_volume_zone_by_id(i);
    if (z->type & CS_VOLUME_ZONE_SOLID)
      n_solid_zones += 1;
  }

  if (n_solid_zones < 1)
    return;
  if (cs_glob_param_cdo_mode == CS_PARAM_CDO_MODE_ONLY)
    return;  /* The current rationale is not compatible with CDO up to now */

  cs_tree_node_t *node_int_cpl
    = cs_tree_get_node(cs_glob_tree, "thermophysical_models/internal_coupling");

  if (node_int_cpl != NULL) {

    if (node_int_cpl->children == NULL)
      return;

    int j = 0;
    int *z_ids;
    BFT_MALLOC(z_ids, n_solid_zones, int);

    for (int i = 0; i < n_zones; i++) {
      const cs_zone_t  *z = cs_volume_zone_by_id(i);
      if (z->type & CS_VOLUME_ZONE_SOLID)
        z_ids[j++] = z->id;
    }

    int coupling_id = cs_internal_coupling_n_couplings();

    cs_internal_coupling_add_volume_zones(n_solid_zones, z_ids);
    BFT_FREE(z_ids);

    {
      cs_internal_coupling_t *cpl = cs_internal_coupling_by_id(coupling_id);

      char i_name[64], e_name[64];
      snprintf(i_name, 63, "auto:internal_coupling_%d_fluid", cpl->id);
      i_name[63] = '\0';
      snprintf(e_name, 63, "auto:internal_coupling_%d_solid", cpl->id);
      e_name[63] = '\0';

      cs_internal_coupling_add_boundary_groups(cpl, i_name, e_name);
    }

    if (n_solid_zones > 0) {
      cs_tree_node_t *ns
        = cs_tree_node_get_child(node_int_cpl, "coupled_scalars");
      /* Add the coupled scalars defined in the GUI */
      for (cs_tree_node_t *tn = cs_tree_node_get_child(ns, "scalar");
           tn != NULL;
           tn = cs_tree_node_get_next_of_name(tn)) {

        const char *scalar_name = cs_tree_node_get_tag(tn, "name");
        int f_id = cs_field_id_by_name(scalar_name);
        if (f_id < 0)
          bft_error(__FILE__, __LINE__, 0,
                    _("Scalar %s does not exist!.\n"), scalar_name);

        cs_internal_coupling_add_entity(f_id);
      }
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add new MEG context info.
 *
 * \param[in]  zone      pointer to associated zone
 * \param[in]  fields    array of field pointers
 * \param[in]  n_fields  number gof field pointers
 *
 * \return: pointer to MEG context info
 */
/*----------------------------------------------------------------------------*/

const cs_gui_volume_meg_context_t *
cs_gui_add_volume_meg_context(const  cs_zone_t   *zone,
                              const  cs_field_t  *fields[],
                              const  int          n_fields)
{
  BFT_REALLOC(_v_meg_contexts,
              _n_v_meg_contexts+1,
              cs_gui_volume_meg_context_t *);

  /* Allocate field pointers at end of structure (be careful of alignment) */

  size_t f_size = (n_fields+1)*sizeof(cs_field_t *);
  size_t b_size = sizeof(cs_gui_volume_meg_context_t);

  int n_contexts = 1 + f_size/b_size;
  if (f_size % b_size)
    n_contexts += 1;

  cs_gui_volume_meg_context_t  *meg_context = NULL;
  BFT_MALLOC(meg_context, n_contexts, cs_gui_volume_meg_context_t);

  meg_context->zone = zone;

  /* Now set field pointers */

  void *f_p = meg_context + 1;
  meg_context->fields = (const cs_field_t **)f_p;

  for (int i = 0; i < n_fields; i++)
    meg_context->fields[i] = fields[i];

  meg_context->fields[n_fields] = NULL;

  /* Now set in structure */

  _v_meg_contexts[_n_v_meg_contexts] = meg_context;
  _n_v_meg_contexts += 1;

  return meg_context;
}

/*----------------------------------------------------------------------------
 * Define user scalar source terms.
 *----------------------------------------------------------------------------*/

void
cs_gui_scalar_source_terms(cs_field_t        *f,
                           const cs_real_t   *restrict pvar,
                           cs_real_t         *restrict tsexp,
                           cs_real_t         *restrict tsimp)
{
  const cs_real_t *restrict cell_f_vol = cs_glob_mesh_quantities->cell_f_vol;
  const cs_real_3_t *restrict cell_cen =
    (const cs_real_3_t *restrict)cs_glob_mesh_quantities->cell_cen;

  const char *formula = NULL;

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
#endif

  const int idarcy = cs_glob_physical_model_flag[CS_GROUNDWATER];

  const int n_zones = cs_volume_zone_n_zones();

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
        cs_real_t *st_vals = NULL;
        BFT_MALLOC(st_vals, 2*n_cells, cs_real_t);

        cs_meg_source_terms(z->name,
                            n_cells,
                            cell_ids,
                            cell_cen,
                            f->name,
                            "scalar_source_term",
                            st_vals);

        cs_real_t sign = 1.0;
        cs_real_t non_linear = 1.0;
        /* for groundwater flow, the user filled in the positive radioactive
           decay rate (lambda) - this source term is always linear:
           -lambda Y^{n+1} */
        if (idarcy > -1) {
          sign = -1.0;
          non_linear = 0.;
        }

        for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
          cs_lnum_t c_id = cell_ids[e_id];
          tsimp[c_id] = cell_f_vol[c_id] * sign * st_vals[2 * e_id + 1];
          tsexp[c_id] = cell_f_vol[c_id] * st_vals[2 * e_id]
                        - non_linear * tsimp[c_id] * pvar[c_id];
        }

        BFT_FREE(st_vals);
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Define user thermal scalar source terms
 *----------------------------------------------------------------------------*/

void
cs_gui_thermal_source_terms(cs_field_t                 *f,
                            const cs_real_t   *restrict pvar,
                            cs_real_t         *restrict tsexp,
                            cs_real_t         *restrict tsimp)
{

  const cs_real_t *restrict cell_f_vol = cs_glob_mesh_quantities->cell_f_vol;

  const cs_real_3_t *restrict cell_cen =
    (const cs_real_3_t *restrict)cs_glob_mesh_quantities->cell_cen;

  const char *formula = NULL;

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
      while (tn != NULL) {

        const char *name = cs_gui_node_get_tag(tn, "name");
        const char *zone_id = cs_gui_node_get_tag(tn, "zone_id");
        if (cs_gui_strcmp(name, f->name) && cs_gui_strcmp(zone_id, z_id_str))
          break;
        tn = cs_tree_node_get_next_of_name(tn);
      }
      formula = cs_tree_node_get_value_str(tn);

      if (formula != NULL) {
        cs_real_t *st_vals = NULL;
        BFT_MALLOC(st_vals, 2*n_cells, cs_real_t);

        cs_meg_source_terms(z->name,
                            n_cells,
                            cell_ids,
                            cell_cen,
                            f->name,
                            "thermal_source_term",
                            st_vals);

        for (cs_lnum_t e_id = 0; e_id < n_cells; e_id++) {
          cs_lnum_t c_id = cell_ids[e_id];

          tsimp[c_id] = cell_f_vol[c_id] * st_vals[2 * e_id + 1];
          tsexp[c_id] = cell_f_vol[c_id] * st_vals[2 * e_id]
                      - tsimp[c_id] * pvar[c_id];
        }

        BFT_FREE(st_vals);
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read GUI defined time tables
 */
/*----------------------------------------------------------------------------*/

void
cs_gui_time_tables(void)
{

  cs_tree_node_t *tt_n = cs_tree_get_node(cs_glob_tree,
                                          "physical_properties/time_tables");

  for (cs_tree_node_t *n = cs_tree_find_node(tt_n, "table");
       n != NULL;
       n = cs_tree_node_get_next_of_name(n)) {

    const char *name  = cs_tree_node_get_tag(n, "name");
    const char *fpath = cs_tree_node_get_tag(n, "file_name");
    const char *delim = cs_tree_node_get_tag(n, "delimiter");

    /* Get number of rows to skip if needed */
    int n_rows_to_skip = 0;
    cs_tree_node_t *node = cs_tree_node_get_child(n, "skip_rows");
    cs_gui_node_get_int(node, &n_rows_to_skip);

    /* Columns */
    int  n_columns = -1; /* Default mode = all */
    int *col_ids   = NULL; /* Default mode : NULL (all columns) */
    node = cs_tree_node_get_child(n, "cols2import");

    if (node != NULL) {
      const char *mode = cs_tree_node_get_value_str(node);
      if (mode != NULL && cs_gui_strcmp(mode, "subset")) {
        n_columns = 0;

        node = cs_tree_node_get_child(n, "col_ids");
        const char *ids_r = cs_tree_node_get_value_str(node);
        char *ids;
        BFT_MALLOC(ids, strlen(ids_r) + 1, char);
        strcpy(ids, ids_r);

        /* Parse line and count number of columns */
        char *token = strtok(ids, ",");
        while (token != NULL) {
          n_columns += 1;
          BFT_REALLOC(col_ids, n_columns, int);
          sscanf(token, "%d", col_ids + (n_columns - 1));

          token = strtok(NULL, ",");
        }

        /* GUI numbering of columns' ids starts at 1 */
        for (int col_id = 0; col_id < n_columns; col_id++)
          col_ids[col_id] -= 1;

        BFT_FREE(ids);

      }
    }

    cs_time_table_t *new_table = cs_time_table_from_csv_file(name,
                                                             fpath,
                                                             delim,
                                                             n_rows_to_skip,
                                                             n_columns,
                                                             col_ids,
                                                             true);

    BFT_FREE(col_ids);

    /* Time offset */
    node = cs_tree_node_get_child(n, "time_offset_choice");
    const char *toffset_choice = cs_tree_node_get_value_str(node);
    if (cs_gui_strcmp(toffset_choice, "yes")) {
      cs_real_t _t_offset = 0.;

      node = cs_tree_node_get_child(n, "time_offset_value");
      cs_gui_node_get_real(node, &_t_offset);
      cs_time_table_set_offset(new_table, _t_offset);
    }

    /* Set headers */
    node = cs_tree_node_get_child(n, "headers_list");
    if (node != NULL) {
      const char *hl_r = cs_tree_node_get_value_str(node);

      char *hl;
      BFT_MALLOC(hl, strlen(hl_r)+1, char);
      strcpy(hl, hl_r);

      int n_headers = 0;
      char **headers = NULL;

      /* Parse line */
      char *token = strtok(hl, ",");
      while (token != NULL) {
        n_headers += 1;
        BFT_REALLOC(headers, n_headers, char *);
        BFT_MALLOC(headers[n_headers - 1], strlen(token) + 1, char);
        strcpy(headers[n_headers - 1], token);

        token = strtok(NULL, ",");
      }

      cs_time_table_set_headers(new_table, n_headers, (const char **)headers);

      /* Free */
      for (int i = 0; i < n_headers; i++)
        BFT_FREE(headers[i]);
      BFT_FREE(headers);

      BFT_FREE(hl);

    }
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
