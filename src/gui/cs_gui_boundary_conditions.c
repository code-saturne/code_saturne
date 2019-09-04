/*============================================================================
 * Management of the GUI parameters file: boundary conditions
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
#include "cs_boundary_conditions.h"
#include "cs_boundary_zone.h"
#include "cs_parameters.h"
#include "cs_gui_util.h"
#include "cs_gui.h"
#include "cs_gui_specific_physics.h"
#include "cs_gui_variables.h"
#include "cs_mesh.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_prototypes.h"
#include "cs_thermal_model.h"
#include "cs_timer.h"
#include "cs_turbulence_model.h"
#include "cs_parall.h"
#include "cs_elec_model.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gui_boundary_conditions.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* debugging switch */
#define _XML_DEBUG_ 0

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* Enum for boundary conditions */

typedef enum {
  DIRICHLET,
  FLOW1,
  HYDRAULIC_DIAMETER,
  TURBULENT_INTENSITY,
  NEUMANN,
  EXCHANGE_COEFF,
  COALFLOW,
  WALL_FUNCTION,
  DIRICHLET_FORMULA,
  DIRICHLET_IMPLICIT,
  NEUMANN_FORMULA,
  NEUMANN_IMPLICIT,
  EXCHANGE_COEFF_FORMULA
} cs_boundary_value_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer on the main boundaries structure */

cs_boundary_t *boundaries = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Return the choice for the scalar of boundary condition type
 *
 * parameters:
 *   nature      <--  nature of boundary condition (inlet, wall, symmetry ..)
 *   label       <--  label of boundary condition
 *   var_sca     <--  name of variable(velocity_pressure, turbulence ...)
 *----------------------------------------------------------------------------*/

static char*
_boundary_choice(const char  *nature,
                 const char  *label,
                 const char  *var_sca,
                 const char  *choice)
{
  char *path = NULL;
  char *c = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "boundary_conditions", nature);
  cs_xpath_add_test_attribute(&path, "label", label);
  cs_xpath_add_element(&path, var_sca);
  cs_xpath_add_attribute(&path, choice);

  c = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return c;
}

/*-----------------------------------------------------------------------------
 * Get status of data for inlet or outlet information.
 *
 * parameters:
 *   nature      <--  nature of the boundary
 *   label       <--  label of the inlet or the outlet
 *   tag         <--  name of researched data
 *   data       -->   value associated to the data
 *----------------------------------------------------------------------------*/

static void
_boundary_status(const char  *nature,
                 const char  *label,
                 const char  *tag,
                 int         *data)
{
  char  *path = NULL;
  int  result = 0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "boundary_conditions", nature);
  cs_xpath_add_test_attribute(&path, "label", label);
  cs_xpath_add_element(&path, tag);
  cs_xpath_add_attribute(&path, "status");

  if (cs_gui_get_status(path, &result))
    *data = result;
  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Get status of data for velocity_pressure node of inlet or outlet information.
 *
 * parameters:
 *   nature      <--  nature of the boundary
 *   label       <--  label of the inlet or the outlet
 *   tag         <--  name of researched data
 *   data       -->   value associated to the data
 *----------------------------------------------------------------------------*/

static void
_boundary_status_vp(const char  *nature,
                    const char  *label,
                    const char  *tag,
                    int         *data)
{
  char  *path = NULL;
  int  result = 0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "boundary_conditions", nature);
  cs_xpath_add_test_attribute(&path, "label", label);
  cs_xpath_add_elements(&path, 2, "velocity_pressure", tag);
  cs_xpath_add_attribute(&path, "status");

  if (cs_gui_get_status(path, &result))
    *data = result;
  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Check if a zone uses a mapped inlet, and return the locator from the
 * associated mapping if this is the case.
 *
 * parameters:
 *   label    <-- label of wall boundary condition
 *   n_faces  <-- number of selected boundary faces
 *   faces    <-- list of selected boundary faces (0 to n-1),
 *                or NULL if no indirection is needed
 *
 * returns:
 *   pointer to mapping locator, or NULL
 *----------------------------------------------------------------------------*/

static  ple_locator_t *
_mapped_inlet(const char                *label,
              cs_lnum_t                  n_faces,
              const cs_lnum_t           *faces)
{
  ple_locator_t *bl = NULL;

  int mapped_inlet = 0;
  _boundary_status("inlet", label, "mapped_inlet", &mapped_inlet);

  if (mapped_inlet) {
    cs_real_t coord_shift[3] = {0., 0., 0.};
    const char *tname[] = {"translation_x",
                           "translation_y",
                           "translation_z"};

    for (int i = 0; i < 3; i++) {
      char *path = cs_xpath_init_path();
      double value = 0.;
      cs_xpath_add_element(&path, "boundary_conditions");
      cs_xpath_add_element(&path, "inlet");
      cs_xpath_add_element(&path, "mapped_inlet");
      cs_xpath_add_element(&path, tname[i]);
      cs_xpath_add_function_text(&path);
      cs_gui_get_double(path, &value);
      BFT_FREE(path);
      coord_shift[i] = value;
    }

    bl = cs_boundary_conditions_map(CS_MESH_LOCATION_CELLS,
                                    cs_glob_mesh->n_cells,
                                    n_faces,
                                    NULL,
                                    faces,
                                    &coord_shift,
                                    0,
                                    0.1);
  }

  return bl;
}

/*-----------------------------------------------------------------------------
 * Value of velocity for sliding wall.
 *
 * parameters:
 *   label       <--  label of wall boundary condition
 *   izone       <--  number of zone
 *----------------------------------------------------------------------------*/

static void
_sliding_wall(const char  *label,
              const int    izone)
{
  char suf[2];
  char *path = NULL;
  double result = 0.0;

  const cs_field_t  *f = cs_field_by_name("velocity");

  for (int i = 0; i < f->dim; i++) {
    path = cs_xpath_init_path();
    cs_xpath_add_element(&path, "boundary_conditions");
    cs_xpath_add_element(&path, "wall");
    cs_xpath_add_test_attribute(&path, "label", label);
    cs_xpath_add_element(&path, "velocity_pressure");
    cs_xpath_add_test_attribute(&path, "choice", "on");
    cs_xpath_add_element(&path, "dirichlet");
    cs_xpath_add_test_attribute(&path, "name", "velocity");
    sprintf(suf, "%i", i);
    cs_xpath_add_test_attribute(&path, "component", suf);
    cs_xpath_add_function_text(&path);

    if (cs_gui_get_double(path, &result)) {
      boundaries->type_code[f->id][izone] = WALL_FUNCTION;
      boundaries->values[f->id][f->dim * izone + i].val1 = result;
    }
    BFT_FREE(path);
  }
}

/*-----------------------------------------------------------------------------
 * Value of roughness for wall
 *
 * parameters:
 *   label       <--  label of boundary condition
 *   izone       <--  number of zone
 *----------------------------------------------------------------------------*/

static void
_wall_roughness(const char  *label,
                const int    izone)
{
  char *path = NULL;
  double result = 0.0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "boundary_conditions", "wall");
  cs_xpath_add_test_attribute(&path, "label", label);
  cs_xpath_add_elements(&path, 2, "velocity_pressure", "roughness");
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &result))
    boundaries->rough[izone] = result;
  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * add notebook variable to formula
 *----------------------------------------------------------------------------*/

static void
_add_notebook_variables(mei_tree_t *ev_law)
{
  char *path = NULL;

  /* number of variable */
  int nbvar = cs_gui_get_tag_count("/physical_properties/notebook/var\n", 1);

  if (nbvar == 0)
    return;

  for (int ivar = 0; ivar < nbvar; ivar++) {
    path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 2, "physical_properties", "notebook");
    cs_xpath_add_element_num(&path, "var", ivar +1);
    cs_xpath_add_attribute(&path, "name");
    char *name = cs_gui_get_attribute_value(path);
    BFT_FREE(path);

    path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 2, "physical_properties", "notebook");
    cs_xpath_add_element_num(&path, "var", ivar +1);
    cs_xpath_add_attribute(&path, "value");
    char *value = cs_gui_get_attribute_value(path);
    double val = atof(value);
    BFT_FREE(path);

    mei_tree_insert(ev_law, name, val);
    BFT_FREE(name);
    BFT_FREE(value);
  }
}

/*-----------------------------------------------------------------------------
 * Get value of data for inlet velocity.
 *
 * parameters:
 *   label       <--  label of the inlet
 *   tag         <--  name of researched data
 *   data       -->   value associated to the data
 *----------------------------------------------------------------------------*/

static void
_inlet_data(const char  *label,
            const char  *tag,
            double      *data)
{
  char  *path = NULL;
  double result = 0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "boundary_conditions", "inlet");
  cs_xpath_add_test_attribute(&path, "label", label);
  cs_xpath_add_elements(&path, 2, "velocity_pressure", tag);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &result))
    *data = result;
  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Formula inlet velocity.
 *
 * parameters:
 *   label       <--  label of the inlet
 *----------------------------------------------------------------------------*/

static char*
_inlet_formula(const char  *label,
               const char  *choice)
{
  char *path = NULL;
  char *form = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "boundary_conditions", "inlet");
  cs_xpath_add_test_attribute(&path, "label", label);
  cs_xpath_add_elements(&path, 2, "velocity_pressure", choice);
  cs_xpath_add_function_text(&path);

  form = cs_gui_get_text_value(path);
  BFT_FREE(path);

  if (form == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: inlet formula for %s, %s is not defined."),
              label, choice);

  return form;
}

/*-----------------------------------------------------------------------------
 * Formula for head loss.
 *
 * parameters:
 *   label       <--  label of the inlet
 *----------------------------------------------------------------------------*/

static char*
_head_loss_formula(const char  *label)
{
  char *path = NULL;
  char *form = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "boundary_conditions", "free_inlet_outlet");
  cs_xpath_add_test_attribute(&path, "label", label);
  cs_xpath_add_elements(&path, 2, "headLoss", "formula");
  cs_xpath_add_function_text(&path);

  form = cs_gui_get_text_value(path);
  BFT_FREE(path);
  return form;
}

/*-----------------------------------------------------------------------------
 * Values for turbulence variable for the current inlet.
 *
 * parameters:
 *   choice      <--  type of choice to calculate turbulence
 *   izone       <--  number of zone
 *----------------------------------------------------------------------------*/

static void
_inlet_turbulence(const char  *choice,
                  int          izone)
{
  char *path1 = NULL;
  char *path2 = NULL;
  double result;

  if (cs_gui_strcmp(choice, "hydraulic_diameter"))
    boundaries->icalke[izone] = 1  ;
  else if(cs_gui_strcmp(choice, "turbulent_intensity"))
    boundaries->icalke[izone] = 2  ;
  else if(cs_gui_strcmp(choice, "formula")) {
    boundaries->icalke[izone] = 0  ;
    return;
  }
  else
    return;

  path1 = cs_xpath_init_path();
  cs_xpath_add_elements(&path1, 2, "boundary_conditions", "inlet");
  cs_xpath_add_test_attribute(&path1, "label", boundaries->label[izone]);
  cs_xpath_add_element(&path1, "turbulence");

  BFT_MALLOC(path2, strlen(path1) + 1, char);
  strcpy(path2, path1);

  cs_xpath_add_element(&path1, "hydraulic_diameter");
  cs_xpath_add_function_text(&path1);

  if (cs_gui_get_double(path1, &result))
    boundaries->dh[izone] = result;

  BFT_FREE(path1);

  if(cs_gui_strcmp(choice, "turbulent_intensity")) {
    cs_xpath_add_element(&path2, "turbulent_intensity");
    cs_xpath_add_function_text(&path2);

    if (cs_gui_get_double(path2, &result))
      boundaries->xintur[izone] = result * 0.01;
  }

  BFT_FREE(path2);
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
_boundary_scalar_init_mei_tree(const char   *formula,
                               const char   *symbols[],
                               int           symbol_size)

{
  int i = 0;
  double ttcabs = cs_glob_time_step->t_cur;
  int    ntcabs = cs_glob_time_step->nt_cur;
  double dtref  = cs_glob_time_step_options->dtref;

  /* return an empty interpreter */
  mei_tree_t *tree = mei_tree_new(formula);

  /* add commun variables */
  mei_tree_insert(tree, "x", 0.0);
  mei_tree_insert(tree, "y", 0.0);
  mei_tree_insert(tree, "z", 0.0);
  mei_tree_insert(tree, "t", ttcabs);
  mei_tree_insert(tree, "dt", dtref);
  mei_tree_insert(tree, "iter", ntcabs);

  /* add variable from notebook */
  _add_notebook_variables(tree);

  /* try to build the interpreter */
  if (mei_tree_builder(tree))
    bft_error(__FILE__, __LINE__, 0,
              _("Error: can not interpret expression: %s\n"), tree->string);

  /* check for symbols */
  for (i = 0; i < symbol_size; ++i)
    if (mei_tree_find_symbol(tree, symbols[i]))
      bft_error(__FILE__, __LINE__, 0,
          _("Error: can not find the required symbol: %s\n"), symbols[i]);

  return tree;
}

/*-----------------------------------------------------------------------------
 * get scalar's values
 *
 * parameters:
 *   nature      <--  nature of boundary condition
 *   izone       <--  number of zone
 *   f_id        <--  field id
 *----------------------------------------------------------------------------*/

static void
_boundary_scalar(const char   *nature,
                 int           izone,
                 int           f_id)
{
  char *path = NULL;
  char *path_commun = NULL;
  char *path2 = NULL;
  char *choice = NULL;
  double result;

  const cs_field_t  *f = cs_field_by_id(f_id);

  int dim = f->dim;

  for (int i = 0; i < dim; i++) {
    path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 2, "boundary_conditions", nature);
    cs_xpath_add_test_attribute(&path, "label", boundaries->label[izone]);
    cs_xpath_add_element(&path, "scalar");
    cs_xpath_add_test_attribute(&path, "name", f->name);

    if (dim > 1)
      cs_xpath_add_element_num(&path, "component", i);

    BFT_MALLOC(path_commun, strlen(path)+1, char);
    strcpy(path_commun, path);

    BFT_MALLOC(path2, strlen(path)+1, char);
    strcpy(path2, path);

    cs_xpath_add_attribute(&path_commun, "choice");
    choice = cs_gui_get_attribute_value(path_commun);

    if (choice != NULL) {
      if (cs_gui_strcmp(choice, "dirichlet")) {
        cs_xpath_add_element(&path, "dirichlet");
        cs_xpath_add_function_text(&path);
        if (cs_gui_get_double(path, &result)) {
          if (cs_gui_strcmp(nature, "wall")) {
            boundaries->type_code[f_id][izone] = WALL_FUNCTION;
          }
          else {
            boundaries->type_code[f_id][izone] = DIRICHLET;
          }
          boundaries->values[f_id][izone * dim + i].val1 = result;
        }
      } else if(cs_gui_strcmp(choice, "neumann")) {
        cs_xpath_add_element(&path, "neumann");
        cs_xpath_add_function_text(&path);
        if (cs_gui_get_double(path, &result)) {
          boundaries->type_code[f_id][izone] = NEUMANN;
          boundaries->values[f_id][izone * dim + i].val3 = result;
        }
      } else if (cs_gui_strcmp(choice, "dirichlet_formula")) {
        cs_xpath_add_element(&path, "dirichlet_formula");
        cs_xpath_add_function_text(&path);
        char *t_value = cs_gui_get_text_value(path);
        if (t_value) {
          const char *sym[] = {f->name};
          boundaries->type_code[f_id][izone] = DIRICHLET_FORMULA;
          boundaries->scalar[f_id][izone * dim + i] =
            _boundary_scalar_init_mei_tree(t_value,
                                           sym,
                                           1);
          BFT_FREE(t_value);
        }
      } else if (cs_gui_strcmp(choice, "neumann_formula")) {
        cs_xpath_add_element(&path, "neumann_formula");
        cs_xpath_add_function_text(&path);
        char *t_value = cs_gui_get_text_value(path);
        if (t_value != NULL) {
          const char *sym[] = {"flux"};
          boundaries->type_code[f_id][izone] = NEUMANN_FORMULA;
          boundaries->scalar[f_id][izone * dim + i] =
            _boundary_scalar_init_mei_tree(t_value,
                                           sym,
                                           1);
          BFT_FREE(t_value);
        }
      } else if (cs_gui_strcmp(choice, "exchange_coefficient_formula")) {
        cs_xpath_add_element (&path, "exchange_coefficient_formula");
        cs_xpath_add_function_text(&path);
        char *t_value = cs_gui_get_text_value(path);
        if (t_value != NULL) {
          const char *sym[] = {f->name,"hc"};
          boundaries->type_code[f_id][izone] = EXCHANGE_COEFF_FORMULA;
          boundaries->scalar[f_id][izone * dim + i] =
            _boundary_scalar_init_mei_tree(t_value,
                                           sym,
                                           2);
          BFT_FREE(t_value);
        }
      } else if (cs_gui_strcmp(choice, "exchange_coefficient")) {
	cs_xpath_add_element(&path, "dirichlet");
        cs_xpath_add_function_text(&path);
	if (cs_gui_get_double(path, &result)) {
	  boundaries->values[f_id][izone * dim + i].val1 = result;
	}
        cs_xpath_add_element(&path2, "exchange_coefficient");
        cs_xpath_add_function_text(&path2);
        if (cs_gui_get_double(path2, &result)) {
          boundaries->type_code[f_id][izone] = EXCHANGE_COEFF;
          boundaries->values[f_id][izone * dim + i].val2 = result;
        }
      } else if (cs_gui_strcmp(choice, "dirichlet_implicit")) {
          boundaries->type_code[f_id][izone] = DIRICHLET_IMPLICIT;
      } else if (cs_gui_strcmp(choice, "neumann_implicit")) {
          boundaries->type_code[f_id][izone] = NEUMANN_IMPLICIT;
      }
    }

    BFT_FREE(choice);
  }

  BFT_FREE(path);
  BFT_FREE(path2);
  BFT_FREE(path_commun);
}

/*-----------------------------------------------------------------------------
 * Get coal's data for inlet. Check if the current zone is an inlet only
 * for an oxydant, of for oxydant and coal.
 *
 * parameters:
 *   izone       <--  number of the current zone
 *   ncharb      <--  number of coals (1 to 3)
 *   nclpch      <--  number of class for eah coal
 *----------------------------------------------------------------------------*/

static void
_inlet_coal(const int         izone,
            const int  *const ncharb,
            const int  *const nclpch)
{
  int    icoal;
  int    iclass;
  int    _n_coals = 0;
  double value;
  char  *path0 = NULL;
  char  *path1 = NULL;
  char  *path2 = NULL;
  char  *path3 = NULL;
  char  *path4 = NULL;
  char  *path5 = NULL;
  char  **list_of_coals = NULL;
  char  **list_of_classes = NULL;

  path0 = cs_xpath_init_path();
  cs_xpath_add_elements(&path0, 2, "boundary_conditions", "inlet");
  cs_xpath_add_test_attribute(&path0, "label", boundaries->label[izone]);
  cs_xpath_add_elements(&path0, 2, "velocity_pressure", "coal");

  BFT_MALLOC(path1, strlen(path0) + 1, char);
  strcpy(path1, path0);
  cs_xpath_add_attribute(&path1, "name");
  list_of_coals = cs_gui_get_attribute_values(path1, &_n_coals);
  BFT_FREE(path1);

  /* if there is no coal, it is an inlet only for oxydant */
  if (_n_coals == 0) {
    boundaries->ientat[izone] = 1;
    boundaries->ientcp[izone] = 0;
    BFT_FREE(path0);
  }
  else {
    if (_n_coals != *ncharb)
      bft_error(__FILE__, __LINE__, 0,
          _("Invalid number of coal-> dp_FCP: %i xml: %i\n"),
          *ncharb, _n_coals);

    boundaries->ientat[izone] = 0;
    boundaries->ientcp[izone] = 1;

    /* loop on number of coals */

    for (icoal = 0; icoal < *ncharb; icoal++) {
      BFT_MALLOC(path2, strlen(path0) + 1, char);
      strcpy(path2, path0);
      /* sprintf(coalname, "%.4s%2.2i", "coal", icoal+1); */
      cs_xpath_add_test_attribute(&path2, "name", list_of_coals[icoal]);

      BFT_MALLOC(path3, strlen(path2) + 1, char);
      strcpy(path3, path2);

      BFT_MALLOC(path4, strlen(path2) + 1, char);
      strcpy(path4, path2);

      /* mass flow rate of coal */

      cs_xpath_add_element(&path3, "flow1");
      cs_xpath_add_function_text(&path3);
      if (cs_gui_get_double(path3, &value))
        boundaries->qimpcp[izone][icoal] = value;

      /* temperature of coal */

      cs_xpath_add_element(&path4, "temperature");
      cs_xpath_add_function_text(&path4);
      if (cs_gui_get_double(path4, &value))
        boundaries->timpcp[izone][icoal] = value;

      /* loop on number of class by coal for ratio (%) stored in distch */

      cs_xpath_add_element(&path2, "ratio");
      BFT_MALLOC(path1, strlen(path2) + 1, char);
      strcpy(path1, path2);
      cs_xpath_add_attribute(&path1, "name");
      int _n_classes = 0;
      list_of_classes = cs_gui_get_attribute_values(path1, &_n_classes);
      BFT_FREE(path1);

      for (iclass = 0; iclass < nclpch[icoal]; iclass++) {
        BFT_MALLOC(path5, strlen(path2) + 1, char);
        strcpy(path5, path2);

        /* sprintf(classname, "%.5s%2.2i", "class", iclass+1); */
        cs_xpath_add_test_attribute(&path5, "name", list_of_classes[iclass]);
        cs_xpath_add_function_text(&path5);

        if (cs_gui_get_double(path5, &value))
          boundaries->distch[izone][icoal][iclass] = value;

        BFT_FREE(path5);
      }

      for (iclass = 0; iclass < _n_classes; iclass++)
        BFT_FREE(list_of_classes[iclass]);
      BFT_FREE(list_of_classes);

      BFT_FREE(path2);
      BFT_FREE(path3);
      BFT_FREE(path4);
    }

    BFT_FREE(path0);
  }

  for (icoal = 0; icoal < _n_coals; icoal++)
    BFT_FREE(list_of_coals[icoal]);
  BFT_FREE(list_of_coals);
}

/*-----------------------------------------------------------------------------
 * Get gas combustion's data for inlet. Check if the current zone is an inlet only
 * for an oxydant, of for oxydant and coal.
 *
 * parameters:
 *   izone       -->  number of the current zone
 *----------------------------------------------------------------------------*/

static void
_inlet_gas(const int         izone)
{
  double value;
  char  *buff  = NULL;
  char  *path0 = NULL;
  char  *path1 = NULL;
  char  *path2 = NULL;

  path0 = cs_xpath_init_path();
  cs_xpath_add_elements(&path0, 2, "boundary_conditions", "inlet");
  cs_xpath_add_test_attribute(&path0, "label", boundaries->label[izone]);
  cs_xpath_add_element(&path0, "velocity_pressure");
  BFT_MALLOC(path1, strlen(path0) + 1, char);
  strcpy(path1, path0);
  BFT_MALLOC(path2, strlen(path0) + 1, char);
  strcpy(path2, path0);
  cs_xpath_add_element(&path0, "gas_type");
  cs_xpath_add_attribute(&path0, "choice");
  buff = cs_gui_get_attribute_value(path0);

  if (cs_gui_strcmp(buff, "oxydant")) {
    boundaries->ientox[izone] = 1;
  }
  else if (cs_gui_strcmp(buff, "fuel")) {
    boundaries->ientfu[izone] = 1;
  }
  else if (cs_gui_strcmp(buff,"unburned")) {
    boundaries->ientgf[izone] = 1;
    cs_xpath_add_element(&path1, "temperature");
    cs_xpath_add_function_text(&path1);
    if (cs_gui_get_double(path1, &value))
      boundaries->tkent[izone] = value;
    cs_xpath_add_element(&path2, "fraction");
    cs_xpath_add_function_text(&path2);
    if (cs_gui_get_double(path2, &value))
      boundaries->fment[izone] = value;
  }
  else if (cs_gui_strcmp(buff,"burned")) {
    boundaries->ientgb[izone] = 1;
    cs_xpath_add_element(&path1, "temperature");
    cs_xpath_add_function_text(&path1);
    if (cs_gui_get_double(path1, &value))
      boundaries->tkent[izone] = value;
    cs_xpath_add_element(&path2, "fraction");
    cs_xpath_add_function_text(&path2);
    if (cs_gui_get_double(path2, &value))
      boundaries->fment[izone] = value;
  }
  BFT_FREE(path0);
  BFT_FREE(path1);
  BFT_FREE(path2);
  BFT_FREE(buff);
}


/*-----------------------------------------------------------------------------
 * Get compressible data for inlet.
 *
 * parameters:
 *   izone       -->  number of the current zone
 *----------------------------------------------------------------------------*/

static void
_inlet_compressible(int   izone)
{
  double value;
  char  *buff  = NULL;
  char  *status  = NULL;
  char  *path0 = NULL;
  char  *path1 = NULL;
  char  *path_value = NULL;
  char  *path2 = NULL;
  char  *path3 = NULL;
  char  *path4 = NULL;

  path0 = cs_xpath_init_path();
  cs_xpath_add_elements(&path0, 2, "boundary_conditions", "inlet");
  cs_xpath_add_test_attribute(&path0, "label", boundaries->label[izone]);
  cs_xpath_add_element(&path0, "velocity_pressure");
  BFT_MALLOC(path1, strlen(path0) + 1, char);
  strcpy(path1, path0);
  BFT_MALLOC(path2, strlen(path0) + 1, char);
  strcpy(path2, path0);
  BFT_MALLOC(path3, strlen(path0) + 1, char);
  strcpy(path3, path0);
  BFT_MALLOC(path4, strlen(path0) + 1, char);
  strcpy(path4, path0);
  cs_xpath_add_element(&path0, "compressible_type");
  cs_xpath_add_attribute(&path0, "choice");
  buff = cs_gui_get_attribute_value(path0);

  if (cs_gui_strcmp(buff, "imposed_inlet")) {
    boundaries->itype[izone] = CS_ESICF;
    cs_xpath_add_element(&path1, "pressure");
    BFT_MALLOC(path_value, strlen(path1) + 1, char);
    strcpy(path_value, path1);
    cs_xpath_add_attribute(&path1, "status");
    status = cs_gui_get_attribute_value(path1);
    if (cs_gui_strcmp(status, "on")) {
      cs_xpath_add_function_text(&path_value);
      if (cs_gui_get_double(path_value, &value))
        boundaries->prein[izone] = value;
    }
    BFT_FREE(status);
    BFT_FREE(path_value);

    cs_xpath_add_element(&path2, "density");
    BFT_MALLOC(path_value, strlen(path2) + 1, char);
    strcpy(path_value, path2);
    cs_xpath_add_attribute(&path2, "status");
    status = cs_gui_get_attribute_value(path2);
    if (cs_gui_strcmp(status, "on")) {
      cs_xpath_add_function_text(&path_value);
      if (cs_gui_get_double(path_value, &value))
        boundaries->rhoin[izone] = value;
    }
    BFT_FREE(status);
    BFT_FREE(path_value);

    cs_xpath_add_element(&path3, "temperature");
    BFT_MALLOC(path_value, strlen(path3) + 1, char);
    strcpy(path_value, path3);
    cs_xpath_add_attribute(&path3, "status");
    status = cs_gui_get_attribute_value(path3);
    if (cs_gui_strcmp(status, "on")) {
      cs_xpath_add_function_text(&path_value);
      if (cs_gui_get_double(path_value, &value))
        boundaries->tempin[izone] = value;
    }
    BFT_FREE(status);
    BFT_FREE(path_value);

    cs_xpath_add_element(&path4, "energy");
    BFT_MALLOC(path_value, strlen(path4) + 1, char);
    strcpy(path_value, path4);
    cs_xpath_add_attribute(&path4, "status");
    status = cs_gui_get_attribute_value(path4);
    if (cs_gui_strcmp(status, "on")) {
      cs_xpath_add_function_text(&path_value);
      if (cs_gui_get_double(path_value, &value))
        boundaries->entin[izone] = value;
    }
    BFT_FREE(status);
    BFT_FREE(path_value);
  }
  else if (cs_gui_strcmp(buff, "subsonic_inlet_PH")) {
    boundaries->itype[izone] = CS_EPHCF;

    cs_xpath_add_element(&path1, "total_pressure");
    BFT_MALLOC(path_value, strlen(path1) + 1, char);
    strcpy(path_value, path1);
    cs_xpath_add_function_text(&path_value);
    if (cs_gui_get_double(path_value, &value))
      boundaries->prein[izone] = value;
    BFT_FREE(status);
    BFT_FREE(path_value);

    cs_xpath_add_element(&path2, "enthalpy");
    BFT_MALLOC(path_value, strlen(path2) + 1, char);
    strcpy(path_value, path2);
    cs_xpath_add_function_text(&path_value);
    if (cs_gui_get_double(path_value, &value))
      boundaries->entin[izone] = value;
    BFT_FREE(status);
    BFT_FREE(path_value);
  }

  BFT_FREE(path0);
  BFT_FREE(path1);
  BFT_FREE(path2);
  BFT_FREE(path3);
  BFT_FREE(path4);
  BFT_FREE(buff);
}

/*-----------------------------------------------------------------------------
 * Get compressible data for inlet.
 *
 * parameters:
 *   izone       -->  number of the current zone
 *----------------------------------------------------------------------------*/

static void
_outlet_compressible(int         izone)
{
  double value;
  char  *buff  = NULL;
  char  *path0 = NULL;
  char  *path1 = NULL;

  path0 = cs_xpath_init_path();
  cs_xpath_add_elements(&path0, 2, "boundary_conditions", "outlet");
  cs_xpath_add_test_attribute(&path0, "label", boundaries->label[izone]);
  BFT_MALLOC(path1, strlen(path0) + 1, char);
  strcpy(path1, path0);
  cs_xpath_add_element(&path0, "compressible_type");
  cs_xpath_add_attribute(&path0, "choice");
  buff = cs_gui_get_attribute_value(path0);

  if (cs_gui_strcmp(buff, "supersonic_outlet")) {
    boundaries->itype[izone] = CS_SSPCF;
  }
  else if (cs_gui_strcmp(buff, "subsonic_outlet")) {
    boundaries->itype[izone] = CS_SOPCF;
    cs_xpath_add_element(&path1, "dirichlet");
    cs_xpath_add_test_attribute(&path1, "name", "pressure");
    cs_xpath_add_function_text(&path1);
    if (cs_gui_get_double(path1, &value))
      boundaries->preout[izone] = value;
  }

  BFT_FREE(path0);
  BFT_FREE(path1);
  BFT_FREE(buff);
}

/*-----------------------------------------------------------------------------
 * Initialize mei tree and check for symbols existence
 *
 * parameters:
 *   formula        <--  mei formula
 *   symbols        <--  array of symbol to check
 *   symbol_size    <--  number of symbol in symbols
 *----------------------------------------------------------------------------*/

static mei_tree_t *_boundary_init_mei_tree(const char *formula,
                                           const char *symbols[],
                                           const int   symbol_size)
{
  int i = 0;

  /* return an empty interpreter */
  mei_tree_t *tree = mei_tree_new(formula);

  /* add commun variables */
  mei_tree_insert(tree, "dt",   0.0);
  mei_tree_insert(tree, "t",    0.0);
  mei_tree_insert(tree, "iter", 0.0);
  mei_tree_insert(tree, "x",    0.0);
  mei_tree_insert(tree, "y",    0.0);
  mei_tree_insert(tree, "z",    0.0);

  /* add variable from notebook */
  _add_notebook_variables(tree);

  /* try to build the interpreter */
  if (mei_tree_builder(tree))
    bft_error(__FILE__, __LINE__, 0,
        _("Error: can not interpret expression: %s\n"), tree->string);

  /* check for symbols */
  for (i = 0; i < symbol_size; ++i)
    if (mei_tree_find_symbol(tree, symbols[i]))
      bft_error(__FILE__, __LINE__, 0,
          _("Error: can not find the required symbol: %s\n"),
          symbols[i]);

  return tree;
}

/*-----------------------------------------------------------------------------
 * Formula for hydraulic head.
 *
 * parameters:
 *   label       <--  label of the inlet
 *----------------------------------------------------------------------------*/

static char*
_hydraulic_head_formula(const char *nature,
                        const char *label)
{
  char *path = NULL;
  char *form = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "boundary_conditions", nature);
  cs_xpath_add_test_attribute(&path, "label", label);
  cs_xpath_add_elements(&path, 1, "dirichlet_formula");
  cs_xpath_add_test_attribute(&path, "name", "hydraulicHead");
  cs_xpath_add_elements(&path, 1, "formula");
  cs_xpath_add_function_text(&path);

  form = cs_gui_get_text_value(path);
  BFT_FREE(path);
  return form;
}

/*-----------------------------------------------------------------------------
 * Get pressure value for darcy (inlet/outlet/groundwater).
 *
 * parameters:
 *   izone       -->  number of the current zone
 *----------------------------------------------------------------------------*/

static void
_boundary_darcy(const char *nature,
                const char *label,
                      int   izone)
{
  double value;
  char  *path = NULL;
  char *choice = NULL;

  choice = _boundary_choice(nature, label, "hydraulicHead", "choice");

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "boundary_conditions", nature);
  cs_xpath_add_test_attribute(&path, "label", boundaries->label[izone]);

  if (cs_gui_strcmp(choice, "dirichlet") || cs_gui_strcmp(choice, "neumann"))
  {
    cs_xpath_add_element(&path, choice);
    cs_xpath_add_test_attribute(&path, "name", "hydraulic_head");
    cs_xpath_add_function_text(&path);
    if (cs_gui_get_double(path, &value))
      boundaries->preout[izone] = value;
  }
  else if (cs_gui_strcmp(choice, "dirichlet_formula"))
  {
    const char *sym[] = {"H"};
    if (_hydraulic_head_formula(nature, label) != NULL)
      boundaries->groundwat[izone] =
          _boundary_init_mei_tree(_hydraulic_head_formula(nature, label), sym, 1);
    else
    {
      bft_printf("Warning : groundwater flow boundary conditions\n");
      bft_printf("          without formula for hydraulic head\n");
    }
  }

  BFT_FREE(path);
  BFT_FREE(choice);
}

/*-----------------------------------------------------------------------------
 * Get pressure value for imposed pressure boundary
 *
 * parameters:
 *   izone       -->  number of the current zone
 *----------------------------------------------------------------------------*/

static void
_boundary_imposed_pressure(const char *nature,
                           const char *label,
                                 int   izone)
{
  double value;
  char  *path = NULL;
  char *choice = NULL;

  choice = _boundary_choice(nature, label, "pressure", "choice");

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "boundary_conditions", nature);
  cs_xpath_add_test_attribute(&path, "label", boundaries->label[izone]);
  cs_xpath_add_element(&path, "dirichlet");
  cs_xpath_add_test_attribute(&path, "name", "pressure");
  cs_xpath_add_function_text(&path);
  if (cs_gui_get_double(path, &value))
    boundaries->preout[izone] = value;

  BFT_FREE(path);
  BFT_FREE(choice);
}

/*----------------------------------------------------------------------------
 * Boundary conditions treatment: global structure initialization
 *
 * parameters:
 *   n_b_faces            <-- number of boundary faces
 *   nozppm               <-- max number of boundary conditions zone
 *   ncharb               <-- number of simulated coals
 *   nclpch               <-- number of simulated class per coals
 *   ncharb               <-- number of simulated gaz for electrical models
 *   izfppp               <-- zone number for each boundary face
 *   idarcy               <-- darcy module activate or not
 *----------------------------------------------------------------------------*/

static void
_init_boundaries(const cs_lnum_t   n_b_faces,
                 const int        *nozppm,
                 const int        *ncharb,
                 const int        *nclpch,
                 int              *izfppp,
                 const int        *idarcy)
{
  cs_lnum_t faces = 0;
  cs_lnum_t zones = 0;
  cs_lnum_t ifac, izone, ith_zone, zone_nbr;
  cs_lnum_t ifbr, i;
  int icharb, iclass;
  char *choice = NULL;
  char *choice_v = NULL;
  char *choice_d = NULL;
  char *nature = NULL;
  char *label = NULL;

  cs_var_t  *vars = cs_glob_var;
  assert(boundaries == NULL);
  int n_fields = cs_field_n_fields();

  zones = cs_gui_boundary_zones_number();

  BFT_MALLOC(boundaries,            1,          cs_boundary_t);
  BFT_MALLOC(boundaries->label,     zones,      char*        );
  BFT_MALLOC(boundaries->nature,    zones,      char*        );
  BFT_MALLOC(boundaries->type_code, n_fields,   int*         );
  BFT_MALLOC(boundaries->values,    n_fields,   cs_val_t*    );
  BFT_MALLOC(boundaries->iqimp,     zones,      int          );
  BFT_MALLOC(boundaries->qimp,      zones,      double       );
  BFT_MALLOC(boundaries->icalke,    zones,      int          );
  BFT_MALLOC(boundaries->dh,        zones,      double       );
  BFT_MALLOC(boundaries->xintur,    zones,      double       );
  BFT_MALLOC(boundaries->rough,     zones,      double       );
  BFT_MALLOC(boundaries->norm,      zones,      double       );
  BFT_MALLOC(boundaries->dirx,      zones,      double       );
  BFT_MALLOC(boundaries->diry,      zones,      double       );
  BFT_MALLOC(boundaries->dirz,      zones,      double       );
  BFT_MALLOC(boundaries->locator,   zones,      ple_locator_t*);

  BFT_MALLOC(boundaries->velocity,  zones,      mei_tree_t*  );
  BFT_MALLOC(boundaries->direction, zones,      mei_tree_t*  );
  BFT_MALLOC(boundaries->headLoss,  zones,      mei_tree_t*  );
  BFT_MALLOC(boundaries->scalar,    n_fields,   mei_tree_t** );
  BFT_MALLOC(boundaries->preout,    zones,      double);

  if (cs_gui_strcmp(vars->model, "solid_fuels"))
  {
    BFT_MALLOC(boundaries->ientat, zones, int      );
    BFT_MALLOC(boundaries->inmoxy, zones, int      );
    BFT_MALLOC(boundaries->timpat, zones, double   );
    BFT_MALLOC(boundaries->ientcp, zones, int      );
    BFT_MALLOC(boundaries->qimpcp, zones, double*  );
    BFT_MALLOC(boundaries->timpcp, zones, double*  );
    BFT_MALLOC(boundaries->distch, zones, double** );

    for (izone=0; izone < zones; izone++)
    {
      BFT_MALLOC(boundaries->qimpcp[izone], *ncharb, double );
      BFT_MALLOC(boundaries->timpcp[izone], *ncharb, double );
      BFT_MALLOC(boundaries->distch[izone], *ncharb, double*);

      for (icharb = 0; icharb < *ncharb; icharb++)
        BFT_MALLOC(boundaries->distch[izone][icharb],
                   nclpch[icharb], double);
    }
  }
  else if (cs_gui_strcmp(vars->model, "gas_combustion")) {
    BFT_MALLOC(boundaries->ientfu,  zones, int   );
    BFT_MALLOC(boundaries->ientox,  zones, int   );
    BFT_MALLOC(boundaries->ientgb,  zones, int   );
    BFT_MALLOC(boundaries->ientgf,  zones, int   );
    BFT_MALLOC(boundaries->tkent,   zones, double);
    BFT_MALLOC(boundaries->fment,   zones, double);
  }
  else if (cs_gui_strcmp(vars->model, "compressible_model")) {
    BFT_MALLOC(boundaries->itype,   zones, int   );
    BFT_MALLOC(boundaries->prein,   zones, double);
    BFT_MALLOC(boundaries->rhoin,   zones, double);
    BFT_MALLOC(boundaries->tempin,  zones, double);
    BFT_MALLOC(boundaries->entin,   zones, double);
  }
  else if (cs_gui_strcmp(vars->model, "groundwater_model")) {
    BFT_MALLOC(boundaries->groundwat, zones,      mei_tree_t*  );
  }
  else {
    boundaries->ientat = NULL;
    boundaries->ientfu = NULL;
    boundaries->ientox = NULL;
    boundaries->ientgb = NULL;
    boundaries->ientgf = NULL;
    boundaries->timpat = NULL;
    boundaries->tkent  = NULL;
    boundaries->fment  = NULL;
    boundaries->ientcp = NULL;
    boundaries->qimpcp = NULL;
    boundaries->timpcp = NULL;
    boundaries->distch = NULL;
  }

  if (cs_gui_strcmp(vars->model, "atmospheric_flows"))
    BFT_MALLOC(boundaries->meteo, zones, cs_meteo_t);
  else
    boundaries->meteo = NULL;

  for (int f_id = 0; f_id < n_fields; f_id++) {
    const cs_field_t  *f = cs_field_by_id(f_id);

    if (f->type & CS_FIELD_VARIABLE) {
      i = f->id;
      BFT_MALLOC(boundaries->type_code[i], zones, int);
      BFT_MALLOC(boundaries->values[i], zones * f->dim, cs_val_t);
      BFT_MALLOC(boundaries->scalar[i], zones * f->dim, mei_tree_t *);
    }
  }

  /* initialize for each zone */
  for (izone = 0; izone < zones; izone++)
  {
    boundaries->iqimp[izone]     = 0;
    boundaries->qimp[izone]      = 0;
    boundaries->norm[izone]      = 0;
    boundaries->icalke[izone]    = 0;
    boundaries->dh[izone]        = 0;
    boundaries->xintur[izone]    = 0;
    boundaries->rough[izone]     = -999;
    boundaries->velocity[izone]  = NULL;
    boundaries->direction[izone] = NULL;
    boundaries->headLoss[izone]  = NULL;
    boundaries->preout[izone]    = 0;
    boundaries->locator[izone]   = NULL;

    if (cs_gui_strcmp(vars->model, "solid_fuels"))
    {
      boundaries->ientat[izone] = 0;
      boundaries->inmoxy[izone] = 1;
      boundaries->ientcp[izone] = 0;
      boundaries->timpat[izone] = 0;

      for (icharb = 0; icharb < *ncharb; icharb++) {
        boundaries->qimpcp[izone][icharb] = 0;
        boundaries->timpcp[izone][icharb] = 0;

        for (iclass = 0; iclass < nclpch[icharb]; iclass++)
          boundaries->distch[izone][icharb][iclass] = 0;
      }
    }

    if (cs_gui_strcmp(vars->model, "gas_combustion")) {
      boundaries->ientfu[izone]  = 0;
      boundaries->ientox[izone]  = 0;
      boundaries->ientgb[izone]  = 0;
      boundaries->ientgf[izone]  = 0;
      boundaries->tkent[izone]   = 0;
      boundaries->fment[izone]   = 0;
    }

    if (cs_gui_strcmp(vars->model, "compressible_model")) {
      boundaries->itype[izone]     = 0;
      boundaries->prein[izone]     = 0;
      boundaries->rhoin[izone]     = 0;
      boundaries->tempin[izone]    = 0;
      boundaries->entin[izone]     = 0;
    }

    if (cs_gui_strcmp(vars->model, "groundwater_model")) {
      boundaries->groundwat[izone] = NULL;
    }

    if (cs_gui_strcmp(vars->model, "atmospheric_flows")) {
      boundaries->meteo[izone].read_data = 0;
      boundaries->meteo[izone].automatic = 0;
    }
  }

  /* Initialization of boundary->type_code and boundary->values */
  for (int f_id = 0; f_id < n_fields; f_id++) {
    const cs_field_t  *f = cs_field_by_id(f_id);

    if (f->type & CS_FIELD_VARIABLE) {
      i = f->id;
      for (izone = 0; izone < zones; izone++)
      {
        boundaries->type_code[i][izone] = -1;
        for (int ii = 0; ii < f->dim; ii++) {
          boundaries->values[i][izone * f->dim + ii].val1 = 1.e30;
          boundaries->values[i][izone * f->dim + ii].val2 = 1.e30;
          boundaries->values[i][izone * f->dim + ii].val3 = 0.;
          boundaries->scalar[i][izone * f->dim + ii] = NULL;
        }
      }
    }
  }

  for (ifac = 0; ifac < n_b_faces; ifac++)
    izfppp[ifac] = 0;

  /* filling of the "boundaries" structure */

  for (izone = 0; izone < zones; izone++) {

    /* nature, label of the ith initialization zone */

    ith_zone = izone + 1;
    nature = cs_gui_boundary_zone_nature(ith_zone);
    label = cs_gui_boundary_zone_label(ith_zone);

    BFT_MALLOC(boundaries->label[izone], strlen(label)+1, char);
    strcpy(boundaries->label[izone], label);

    BFT_MALLOC(boundaries->nature[izone], strlen(nature)+1, char);
    strcpy(boundaries->nature[izone], nature);

    if (cs_gui_strcmp(nature, "inlet")) {

      if (*idarcy < 0) {

        char *i_formula = NULL;

        /* Inlet: VELOCITY */
        choice_v = _boundary_choice("inlet", label, "velocity_pressure", "choice");
        choice_d = _boundary_choice("inlet", label, "velocity_pressure", "direction");

        if (cs_gui_strcmp(choice_v, "norm"))
        {
          _inlet_data(label, "norm", &boundaries->norm[izone]);
        }
        else if (cs_gui_strcmp(choice_v, "flow1")) {
          _inlet_data(label, "flow1", &boundaries->qimp[izone]);
          boundaries->iqimp[izone] = 1;
        }
        else if (cs_gui_strcmp(choice_v, "flow2")) {
          _inlet_data(label, "flow2", &boundaries->qimp[izone]);
          boundaries->iqimp[izone] = 2;
        }
        else if (cs_gui_strcmp(choice_v, "norm_formula")) {
          const char *sym[] = {"u_norm"};
          i_formula = _inlet_formula(label, choice_v);
          boundaries->velocity[izone]
            = _boundary_init_mei_tree(i_formula, sym, 1);
        }
        else if (cs_gui_strcmp(choice_v, "flow1_formula")) {
          const char *sym[] = {"q_m"};
          i_formula = _inlet_formula(label, choice_v);
          boundaries->velocity[izone]
            = _boundary_init_mei_tree(i_formula, sym, 1);
          boundaries->iqimp[izone] = 1;
        }
        else if (cs_gui_strcmp(choice_v, "flow2_formula")) {
          const char *sym[] = {"q_v"};
          i_formula = _inlet_formula(label, choice_v);
          boundaries->velocity[izone]
            = _boundary_init_mei_tree(i_formula, sym, 1);
          boundaries->iqimp[izone] = 2;
        }

        if (   cs_gui_strcmp(choice_d, "coordinates")
            || cs_gui_strcmp(choice_d, "translation"))
        {
          _inlet_data(label, "direction_x", &boundaries->dirx[izone]);
          _inlet_data(label, "direction_y", &boundaries->diry[izone]);
          _inlet_data(label, "direction_z", &boundaries->dirz[izone]);
        }
        else if (cs_gui_strcmp(choice_d, "formula")) {
          const char *sym[] = {"dir_x", "dir_y", "dir_z"};
          i_formula = _inlet_formula(label, "direction_formula");
          boundaries->direction[izone]
            = _boundary_init_mei_tree(i_formula, sym, 3);
        }
        BFT_FREE(i_formula);
        BFT_FREE(choice_v);
        BFT_FREE(choice_d);
      }

      /* Inlet: data for COAL COMBUSTION */
      if (cs_gui_strcmp(vars->model, "solid_fuels"))
      {
        _inlet_data(label, "temperature", &boundaries->timpat[izone]);
        double inmoxy;
        _inlet_data(label, "oxydant", &inmoxy);
        boundaries->inmoxy[izone] = inmoxy;
        _inlet_coal(izone, ncharb, nclpch);
      }

      /* Inlet: data for GAS COMBUSTION */
      if (cs_gui_strcmp(vars->model, "gas_combustion"))
        _inlet_gas(izone);

      /* Inlet: data for COMPRESSIBLE MODEL */
      if (cs_gui_strcmp(vars->model, "compressible_model"))
        _inlet_compressible(izone);

      /* Inlet: data for ATMOSPHERIC FLOWS */
      if (cs_gui_strcmp(vars->model, "atmospheric_flows"))
      {
        _boundary_status_vp("inlet", label, "meteo_data",
                            &boundaries->meteo[izone].read_data);
        _boundary_status_vp("inlet", label, "meteo_automatic",
                            &boundaries->meteo[izone].automatic);
      }

      /* Inlet: data for darcy */
      if (cs_gui_strcmp(vars->model, "groundwater_model"))
        _boundary_darcy(nature, label, izone);

      /* Inlet: TURBULENCE */
      choice = _boundary_choice("inlet", label, "turbulence", "choice");
      _inlet_turbulence(choice, izone);
      BFT_FREE(choice);

    }
    else if (cs_gui_strcmp(nature, "wall")) {
      /* Wall: VELOCITY */
      choice = _boundary_choice("wall", label, "velocity_pressure", "choice");

      if (cs_gui_strcmp(choice, "on"))
        _sliding_wall(label, izone);

      BFT_FREE(choice);

      /* Wall: ROUGH */
      _wall_roughness(label, izone);
    }
    else if (cs_gui_strcmp(nature, "outlet")) {
      /* Outlet: data for ATMOSPHERIC FLOWS */
      if (cs_gui_strcmp(vars->model, "atmospheric_flows")) {
        _boundary_status_vp("outlet", label, "meteo_data",
                            &boundaries->meteo[izone].read_data);
        _boundary_status_vp("outlet", label, "meteo_automatic",
                            &boundaries->meteo[izone].automatic);
      }

      /* Outlet: data for COMPRESSIBLE MODEL */
      if (cs_gui_strcmp(vars->model, "compressible_model"))
        _outlet_compressible(izone);

      /* Inlet: data for darcy */
      if (cs_gui_strcmp(vars->model, "groundwater_model"))
        _boundary_darcy(nature, label, izone);
    }
    else if (cs_gui_strcmp(nature, "free_inlet_outlet")) {
      const char *sym[] = {"K"};
      if (_head_loss_formula(label) != NULL)
        boundaries->headLoss[izone] =
            _boundary_init_mei_tree(_head_loss_formula(label), sym, 1);
      else
      {
        bft_printf("Warning : free inlet outlet boundary conditions\n");
        bft_printf("          without external head loss definition\n");
      }
    }
    else if (cs_gui_strcmp(nature, "imposed_p_outlet")) {
      _boundary_imposed_pressure(nature, label, izone);
    }
    else if (cs_gui_strcmp(nature, "groundwater")) {
      _boundary_darcy(nature, label, izone);
    }

    /* for each zone */
    if (!cs_gui_strcmp(nature, "symmetry")) {

      /* Thermal scalar */

      cs_field_t *f = cs_thermal_model_field();

      if (f != NULL) {
        if (boundaries->meteo == NULL)
          _boundary_scalar(nature,
                           izone,
                           f->id);
        else if (boundaries->meteo[izone].read_data == 0)
          _boundary_scalar(nature,
                           izone,
                           f->id);
      }

      /* Meteo scalars */
      if (cs_gui_strcmp(vars->model, "atmospheric_flows")) {
        char *path_meteo = cs_xpath_init_path();
        int size = -1;
        cs_xpath_add_elements(&path_meteo, 3,
                               "thermophysical_models",
                               "atmospheric_flows",
                               "variable");
        size = cs_gui_get_nb_element(path_meteo);
        BFT_FREE(path_meteo);

        if (boundaries->meteo[izone].read_data == 0) {
          for (int j = 0; j < size; j++)
          {
            path_meteo = cs_xpath_init_path();
            cs_xpath_add_elements(&path_meteo, 2,
                                   "thermophysical_models",
                                   "atmospheric_flows");
            cs_xpath_add_element_num(&path_meteo, "variable", j +1);
            cs_xpath_add_element(&path_meteo, "name");
            cs_xpath_add_function_text(&path_meteo);
            char *name = cs_gui_get_text_value(path_meteo);

            cs_field_t *c = cs_field_by_name_try(name);
            BFT_FREE(path_meteo);

            _boundary_scalar(nature,
                             izone,
                             c->id);
          }
        }
      }

      /* Electric arc scalars */
      if (cs_gui_strcmp(vars->model, "joule_effect")) {
        char *path_elec = cs_xpath_init_path();
        int size = -1;
        cs_xpath_add_elements(&path_elec, 3,
                               "thermophysical_models",
                               "joule_effect",
                               "variable");
        size = cs_gui_get_nb_element(path_elec);
        BFT_FREE(path_elec);

        for (int j = 0; j < size; j++)
        {
          path_elec = cs_xpath_init_path();
          cs_xpath_add_elements(&path_elec, 2,
                                 "thermophysical_models",
                                 "joule_effect");
          cs_xpath_add_element_num(&path_elec, "variable", j +1);
          cs_xpath_add_attribute(&path_elec, "name");
          char *name = cs_gui_get_attribute_value(path_elec);

          cs_field_t *c = cs_field_by_name_try(name);

          BFT_FREE(path_elec);
          _boundary_scalar(nature,
                           izone,
                           c->id);
        }
      }

      /* User scalars */
      for (int f_id = 0; f_id < n_fields; f_id++) {
        const cs_field_t  *fu = cs_field_by_id(f_id);

        if (   (fu->type & CS_FIELD_VARIABLE)
            && (fu->type & CS_FIELD_USER))
          _boundary_scalar(nature,
                           izone,
                           fu->id);
      }
    }

    BFT_FREE(nature);
    BFT_FREE(label);
  }  /* for izones */

  int overlap_error[4] = {0, -1, -1, -1};

  for (izone = 0 ; izone < zones ; izone++) {

    zone_nbr = cs_gui_boundary_zone_number(izone + 1);

    if (nozppm && zone_nbr > *nozppm)
      bft_error(__FILE__, __LINE__, 0,
                _("zone's label number %i is greater than %i,"
                  " the maximum allowed \n"), zone_nbr, *nozppm);

    const cs_lnum_t *face_ids
      = cs_gui_get_boundary_faces(boundaries->label[izone],
                                  &faces);

    /* check if faces are already marked with a zone number */

    for (ifac = 0; ifac < faces; ifac++) {
      ifbr = face_ids[ifac];
      if (izfppp[ifbr] > 0) {
        if (overlap_error[0] == 0) {
          overlap_error[0] = 1;
          overlap_error[1] = izfppp[ifbr];
          overlap_error[2] = zone_nbr;
          overlap_error[3] = izone;
        }
        izfppp[ifbr] = - izfppp[ifbr];
      }
      else {
        izfppp[ifbr] = zone_nbr;
      }
    } /* for ifac */
  } /*  for izone */

  /* Check for zone overlap errors */

  cs_parall_max(1, CS_INT_TYPE, overlap_error);

  if (overlap_error[0] > 0) {

    int old_ref, err_ref, err_izone;

    err_ref = -1, err_izone = -1;

    old_ref = overlap_error[1];
    cs_parall_max(1, CS_INT_TYPE, &old_ref);
    if (old_ref == overlap_error[1])
      err_ref = overlap_error[2];
    cs_parall_max(1, CS_INT_TYPE, &err_ref);
    if (err_ref == overlap_error[2])
      err_izone = overlap_error[3];
    cs_parall_max(1, CS_INT_TYPE, &err_izone);

    if (cs_glob_rank_id < 1)
      bft_printf(_("\n"
                   "Error: boundary face zone overlap\n"
                   "======\n\n"
                   "Zone %i (\"%s\") contains at least\n"
                   "one face already marked with zone number %i.\n"),
                 err_ref, boundaries->label[err_izone], old_ref);

    cs_boundary_conditions_error(izfppp, _("zone number"));

  }

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Boundary conditions treatment
 *
 * Remember: rdoccl[k][j][i] = rcodcl[ k * dim1 *dim2 + j *dim1 + i]
 *
 * Fortran Interface:
 *
 * subroutine uiclim
 * *****************
 *
 * integer          idarcy           <-- darcy module activate or not
 * integer          nozppm           <-- max number of boundary conditions zone
 * integer          ncharm           <-- maximal number of coals
 * integer          ncharb           <-- number of simulated coals
 * integer          nclpch           <-- number of simulated class per coals
 * integer          ngasg            <-- number of simulated gas for electrical models
 * integer          iqimp            <-- 1 if flow rate is applied
 * integer          icalke           <-- 1 for automatic turbulent boundary conditions
 * integer          ientat           <-- 1 for air temperature boundary conditions (coal)
 * integer          ientcp           <-- 1 for coal temperature boundary conditions (coal)
 * INTEGER          INMOXY           <-- coal: number of oxydant for the current inlet
 * integer          ientox           <-- 1 for an air fow inlet (gas combustion)
 * integer          ientfu           <-- 1 for fuel flow inlet (gas combustion)
 * integer          ientgb           <-- 1 for burned gas inlet (gas combustion)
 * integer          ientgf           <-- 1 for unburned gas inlet (gas combustion)
 * integer          iprofm           <-- atmospheric flows: on/off for profile from data
 * integer          iautom           <-- atmospheric flows: auto inlet/outlet flag
 * integer          itypfb           <-- type of boundary for each face
 * integer          izfppp           <-- zone number for each boundary face
 * integer          icodcl           <-- boundary conditions array type
 * double precision surfbo           <-- boundary faces surface
 * double precision cgdfbo           <-- boundary faces center of gravity
 * double precision qimp             <-- inlet flow rate
 * double precision qimpat           <-- inlet air flow rate (coal)
 * double precision qimpcp           <-- inlet coal flow rate (coal)
 * double precision dh               <-- hydraulic diameter
 * double precision xintur           <-- turbulent intensity
 * double precision timpat           <-- air temperature boundary conditions (coal)
 * double precision timpcp           <-- inlet coal temperature (coal)
 * double precision tkent            <-- inlet temperature (gas combustion)
 * double precision fment            <-- Mean Mixture Fraction at Inlet (gas combustion)
 * double precision distch           <-- ratio for each coal
 * integer          nvar             <-- dimension for rcodcl
 * double precision rcodcl           <-- boundary conditions array value
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiclim, UICLIM)(const int  *idarcy,
                               const int  *nozppm,
                               const int  *ncharm,
                               const int  *ncharb,
                               const int  *nclpch,
                               int        *iqimp,
                               int        *icalke,
                               int        *ientat,
                               int        *ientcp,
                               int        *inmoxy,
                               int        *ientox,
                               int        *ientfu,
                               int        *ientgf,
                               int        *ientgb,
                               int        *iprofm,
                               int        *iautom,
                               int        *itypfb,
                               int        *izfppp,
                               int        *icodcl,
                               double     *surfbo,
                               double     *cdgfbo,
                               double     *qimp,
                               double     *qimpat,
                               double     *qimpcp,
                               double     *dh,
                               double     *xintur,
                               double     *timpat,
                               double     *timpcp,
                               double     *tkent,
                               double     *fment,
                               double     *distch,
                               int        *nvar,
                               double     *rcodcl)
{
  cs_lnum_t faces = 0;
  int zones = 0;
  int zone_nbr;
  cs_lnum_t ifbr;
  cs_lnum_t i, ivar;
  double t0;
  double norm = 0.;
  double X[3];

  cs_var_t  *vars = cs_glob_var;
  const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const cs_time_step_t *ts = cs_glob_time_step;

  zones = cs_gui_boundary_zones_number();

  /* First iteration only: memory allocation */

  if (boundaries == NULL)
    _init_boundaries(n_b_faces, nozppm, ncharb, nclpch,
                     izfppp, idarcy);

  /* At each time-step, loop on boundary faces:
     One sets itypfb, rcodcl and icodcl thanks to the arrays
     of the structures "conditions.limites" defined
     in the first part ofthe function */

  /* initialize number of variable field */
  int n_fields = cs_field_n_fields();

  for (int izone = 0; izone < zones; izone++) {

    int ith_zone = izone + 1;

    zone_nbr = cs_gui_boundary_zone_number(ith_zone);

    const cs_lnum_t *face_ids
      = cs_gui_get_boundary_faces(boundaries->label[izone],
                                  &faces);

    /* Mapped inlet? */

    if (   cs_gui_strcmp(boundaries->nature[izone], "inlet")
        && boundaries->locator[izone] == NULL)
      boundaries->locator[izone] = _mapped_inlet(boundaries->label[izone],
                                                 faces,
                                                 face_ids);

    /* for each field */
    for (int f_id = 0; f_id < n_fields; f_id++) {
      const cs_field_t  *f = cs_field_by_id(f_id);
      const int var_key_id = cs_field_key_id("variable_id");
      ivar = cs_field_get_key_int(f, var_key_id) -1;

      if (f->type & CS_FIELD_VARIABLE) {
        switch (boundaries->type_code[f->id][izone]) {
          case NEUMANN:
            for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
              ifbr = face_ids[ifac];
              for (i = 0; i < f->dim; i++) {
                icodcl[(ivar + i)*n_b_faces + ifbr] = 3;
                rcodcl[2 * n_b_faces * (*nvar) + (ivar + i) * n_b_faces + ifbr]
                  = boundaries->values[f->id][izone * f->dim + i].val3;
              }
            }
            break;

          case DIRICHLET:
            for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
              ifbr = face_ids[ifac];
              /* if wall_function <-- icodcl[ivar *n_b_faces + ifbr] = 1; */
              for (i = 0; i < f->dim; i++) {
                icodcl[(ivar + i) * n_b_faces + ifbr] = 1;
                rcodcl[0 * n_b_faces * (*nvar) + (ivar + i) * n_b_faces + ifbr]
                  = boundaries->values[f->id][izone * f->dim + i].val1;
              }
            }
            break;

          case WALL_FUNCTION:
            for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
              ifbr = face_ids[ifac];
              for (i = 0; i < f->dim; i++) {
                icodcl[(ivar + i) * n_b_faces + ifbr] = 5;
                if (boundaries->rough[izone] >= 0.0)
                  icodcl[(ivar + i) * n_b_faces + ifbr] = 6;
                rcodcl[0 * n_b_faces * (*nvar) + (ivar + i) * n_b_faces + ifbr]
                  = boundaries->values[f->id][izone * f->dim + i].val1;
              }
            }
            break;

          case EXCHANGE_COEFF:
            for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
              ifbr = face_ids[ifac];
              for (i = 0; i < f->dim; i++) {
                icodcl[(ivar + i) * n_b_faces + ifbr] = 5;
                rcodcl[0 * n_b_faces * (*nvar) + (ivar + i) * n_b_faces + ifbr]
                  = boundaries->values[f->id][izone * f->dim + i].val1;
                rcodcl[1 * n_b_faces * (*nvar) + (ivar + i) * n_b_faces + ifbr]
                  = boundaries->values[f->id][izone * f->dim + i].val2;
              }
            }
            break;

          case DIRICHLET_FORMULA:
            for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
              ifbr = face_ids[ifac];
              for (i = 0; i < f->dim; i++) {
                mei_tree_t *ev_formula = boundaries->scalar[f->id][izone * f->dim + i];
                mei_tree_insert(ev_formula, "t", ts->t_cur);
                mei_tree_insert(ev_formula, "dt", cs_glob_time_step_options->dtref);
                mei_tree_insert(ev_formula, "iter", ts->nt_cur);
                mei_tree_insert(ev_formula, "x", cdgfbo[3 * ifbr + 0]);
                mei_tree_insert(ev_formula, "y", cdgfbo[3 * ifbr + 1]);
                mei_tree_insert(ev_formula, "z", cdgfbo[3 * ifbr + 2]);

                /* add variable from notebook */
                _add_notebook_variables(ev_formula);

                icodcl[(ivar + i) *n_b_faces + ifbr] = 1;
                mei_evaluate(ev_formula);
                if (f->dim > 1)
                {
                  char *name = NULL;
                  BFT_MALLOC(name, strlen(f->name) + 4, char);
                  sprintf(name, "%s[%d]", f->name, i);
                  rcodcl[0 * n_b_faces * (*nvar) + (ivar + i) * n_b_faces + ifbr]
                    = mei_tree_lookup(ev_formula, name);
                  BFT_FREE(name);
                } else {
                  rcodcl[0 * n_b_faces * (*nvar) + (ivar + i) * n_b_faces + ifbr]
                    = mei_tree_lookup(ev_formula, f->name);
                }
              }
            }
            break;

          case NEUMANN_FORMULA:
            for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
              ifbr = face_ids[ifac];
              for (i = 0; i < f->dim; i++) {
                icodcl[(ivar + i) *n_b_faces + ifbr] = 3;
                mei_tree_t *ev_formula = boundaries->scalar[f->id][izone * f->dim + i];
                mei_tree_insert(ev_formula, "t", ts->t_cur);
                mei_tree_insert(ev_formula, "dt", cs_glob_time_step_options->dtref);
                mei_tree_insert(ev_formula, "iter", ts->nt_cur);

                /* add variable from notebook */
                _add_notebook_variables(ev_formula);

                mei_evaluate(ev_formula);
                rcodcl[2 * n_b_faces * (*nvar) + (ivar + i) * n_b_faces + ifbr]
                  = mei_tree_lookup(ev_formula, "flux");
              }
            }
            break;

          case EXCHANGE_COEFF_FORMULA:
            for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
              ifbr = face_ids[ifac];
              for (i = 0; i < f->dim; i++) {
                icodcl[(ivar + i) *n_b_faces + ifbr] = 5;
                mei_tree_t *ev_formula = boundaries->scalar[f->id][izone * f->dim + i];
                mei_tree_insert(ev_formula, "t", ts->t_cur);
                mei_tree_insert(ev_formula, "dt", cs_glob_time_step_options->dtref);
                mei_tree_insert(ev_formula, "iter", ts->nt_cur);

                /* add variable from notebook */
                _add_notebook_variables(ev_formula);

                mei_evaluate(ev_formula);
                if (f->dim > 1)
                {
                  char *name = NULL;
                  BFT_MALLOC(name, strlen(f->name) + 4, char);
                  sprintf(name, "%s[%d]", f->name, i);
                  rcodcl[0 * n_b_faces * (*nvar) + (ivar + i) * n_b_faces + ifbr]
                      = mei_tree_lookup(ev_formula, name);
                  BFT_FREE(name);
                } else {
                  rcodcl[0 * n_b_faces * (*nvar) + (ivar + i) * n_b_faces + ifbr]
                      = mei_tree_lookup(ev_formula, f->name);
                }
                rcodcl[1 * n_b_faces * (*nvar) + (ivar + i) * n_b_faces + ifbr]
                    = mei_tree_lookup(ev_formula, "hc");
              }
            }
            break;
        }
      } /* switch */
    }

    if (cs_gui_strcmp(vars->model_value, "joule")) {
      if (cs_glob_elec_option->ielcor == 1) {
        const cs_field_t  *f = CS_F_(potr);
        const int var_key_id = cs_field_key_id("variable_id");
        ivar = cs_field_get_key_int(f, var_key_id) -1;

        for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
          ifbr = face_ids[ifac];
          rcodcl[ivar * n_b_faces + ifbr] *= cs_glob_elec_option->coejou;
        }

        if (cs_glob_elec_option->ieljou == 2 || cs_glob_elec_option->ieljou == 4) {
          const cs_field_t  *fi = CS_F_(poti);
          ivar = cs_field_get_key_int(fi, var_key_id) -1;
          for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
            ifbr = face_ids[ifac];
            rcodcl[ivar * n_b_faces + ifbr] *= cs_glob_elec_option->coejou;
          }
        }
      }
    }

    if (cs_gui_strcmp(vars->model_value, "arc")) {
      const cs_field_t  *f = CS_F_(potr);
      const int var_key_id = cs_field_key_id("variable_id");
      ivar = cs_field_get_key_int(f, var_key_id) -1;

      if (   boundaries->type_code[f->id][izone] == DIRICHLET_IMPLICIT
          && cs_glob_elec_option->ielcor == 1)
        for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
          ifbr = face_ids[ifac];
          icodcl[ivar * n_b_faces + ifbr] = 5;
          rcodcl[ivar * n_b_faces + ifbr] = cs_glob_elec_option->pot_diff;
        }

      const cs_field_t  *fp = cs_field_by_name_try("vec_potential");
      ivar = cs_field_get_key_int(fp, var_key_id) -1;

      if (boundaries->type_code[fp->id][izone] == NEUMANN_IMPLICIT)
        for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
          ifbr = face_ids[ifac];
          cs_lnum_t iel = b_face_cells[ifbr];
          icodcl[ivar *n_b_faces + ifbr] = 5;
          icodcl[(ivar+1) *n_b_faces + ifbr] = 5;
          icodcl[(ivar+2) *n_b_faces + ifbr] = 5;
          rcodcl[ivar * n_b_faces + ifbr] = fp->val_pre[3*iel];
          rcodcl[(ivar+1) * n_b_faces + ifbr] = fp->val_pre[3*iel+1];
          rcodcl[(ivar+2) * n_b_faces + ifbr] = fp->val_pre[3*iel+2];
        }
    }

    /* velocity and specific physics */
    if (cs_gui_strcmp(boundaries->nature[izone], "inlet")) {
      char *choice_v = _boundary_choice("inlet",
                                        boundaries->label[izone],
                                        "velocity_pressure",
                                        "choice");
      char *choice_d = _boundary_choice("inlet",
                                        boundaries->label[izone],
                                        "velocity_pressure",
                                        "direction");

      /* Update the depending zone's arrays (iqimp, dh, xintur, icalke, qimp,...)
         because they are initialized at each time step
         in PRECLI and PPPRCL routines */

      /* data by zone */
      iqimp[zone_nbr-1]  = boundaries->iqimp[izone];
      dh[zone_nbr-1]     = boundaries->dh[izone];
      xintur[zone_nbr-1] = boundaries->xintur[izone];
      icalke[zone_nbr-1] = boundaries->icalke[izone];

      if (cs_gui_strcmp(vars->model, "solid_fuels")) {
        ientat[zone_nbr-1] = boundaries->ientat[izone];
        inmoxy[zone_nbr-1] = boundaries->inmoxy[izone];
        ientcp[zone_nbr-1] = boundaries->ientcp[izone];
        qimpat[zone_nbr-1] = boundaries->qimp[izone];
        timpat[zone_nbr-1] = boundaries->timpat[izone];

        for (int icharb = 0; icharb < *ncharb; icharb++) {
          int ich = icharb *(*nozppm)+zone_nbr-1;
          qimpcp[ich] = boundaries->qimpcp[izone][icharb];
          timpcp[ich] = boundaries->timpcp[izone][icharb];

          for (int iclass = 0; iclass < nclpch[icharb]; iclass++) {
            int icl = iclass * (*nozppm) * (*ncharm) + ich;
            distch[icl] = boundaries->distch[izone][icharb][iclass];
          }
        }
      }
      else if (cs_gui_strcmp(vars->model, "gas_combustion")) {
        ientfu[zone_nbr-1] = boundaries->ientfu[izone];
        ientox[zone_nbr-1] = boundaries->ientox[izone];
        ientgb[zone_nbr-1] = boundaries->ientgb[izone];
        ientgf[zone_nbr-1] = boundaries->ientgf[izone];
        tkent[zone_nbr-1]  = boundaries->tkent[izone];
        fment[zone_nbr-1]  = boundaries->fment[izone];

        if (   cs_gui_strcmp(choice_v, "flow1_formula")
            || cs_gui_strcmp(choice_v, "flow2_formula")) {
          mei_tree_insert(boundaries->velocity[izone], "t", ts->t_cur);
          mei_tree_insert(boundaries->velocity[izone], "dt", cs_glob_time_step_options->dtref);
          mei_tree_insert(boundaries->velocity[izone], "iter", ts->nt_cur);

          /* add variable from notebook */
          _add_notebook_variables(boundaries->velocity[izone]);

          mei_evaluate(boundaries->velocity[izone]);

          if (cs_gui_strcmp(choice_v, "flow1_formula"))
            qimp[zone_nbr-1] = mei_tree_lookup(boundaries->velocity[izone], "q_m");
          else if (cs_gui_strcmp(choice_v, "flow2_formula"))
            qimp[zone_nbr-1] = mei_tree_lookup(boundaries->velocity[izone], "q_v");
        }
        else {
          qimp[zone_nbr-1] = boundaries->qimp[izone];
        }
      }
      else if (cs_gui_strcmp(vars->model, "compressible_model")) {
        const int var_key_id = cs_field_key_id("variable_id");

        if (boundaries->itype[izone] == CS_ESICF ||
            boundaries->itype[izone] == CS_EPHCF) {
          const cs_field_t  *fp = cs_field_by_name_try("pressure");
          int ivarp = cs_field_get_key_int(fp, var_key_id) -1;
          const cs_field_t  *fe = cs_field_by_name_try("total_energy");
          int ivare = cs_field_get_key_int(fe, var_key_id) -1;

          for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
            ifbr = face_ids[ifac];
            rcodcl[ivarp * n_b_faces + ifbr] = boundaries->prein[izone];
            rcodcl[ivare * n_b_faces + ifbr] = boundaries->entin[izone];
          }
        }

        if (boundaries->itype[izone] == CS_ESICF) {
          cs_field_t *b_rho = cs_field_by_name_try("boundary_density");
          const cs_field_t  *ft = cs_field_by_name_try("temperature");
          int ivart = cs_field_get_key_int(ft, var_key_id) -1;

          for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
            ifbr = face_ids[ifac];
            rcodcl[ivart * n_b_faces + ifbr] = boundaries->tempin[izone];
            b_rho->val[ifbr] = boundaries->rhoin[izone];
          }
        }
      }
      else {
        if (   cs_gui_strcmp(choice_v, "flow1_formula")
            || cs_gui_strcmp(choice_v, "flow2_formula")) {
          mei_tree_insert(boundaries->velocity[izone], "t", ts->t_cur);
          mei_tree_insert(boundaries->velocity[izone], "dt",
                          cs_glob_time_step_options->dtref);
          mei_tree_insert(boundaries->velocity[izone], "iter", ts->nt_cur);

          /* add variable from notebook */
          _add_notebook_variables(boundaries->velocity[izone]);

          mei_evaluate(boundaries->velocity[izone]);

          if (cs_gui_strcmp(choice_v, "flow1_formula"))
            qimp[zone_nbr-1] = mei_tree_lookup(boundaries->velocity[izone], "q_m");
          else if (cs_gui_strcmp(choice_v, "flow2_formula"))
            qimp[zone_nbr-1] = mei_tree_lookup(boundaries->velocity[izone], "q_v");
        }
        else {
          qimp[zone_nbr-1] = boundaries->qimp[izone];
        }
      }

      /* data by boundary faces */

      int inlet_type = CS_INLET;

      if (cs_gui_strcmp(vars->model, "compressible_model"))
        inlet_type = boundaries->itype[izone];
      else {
        int convective_inlet = 0;
        _boundary_status("inlet", boundaries->label[izone],
                         "convective_inlet", &convective_inlet);
        if (convective_inlet)
          inlet_type = CS_CONVECTIVE_INLET;
      }

      for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
        ifbr = face_ids[ifac];

        /* zone number and nature of boundary */
        izfppp[ifbr] = zone_nbr;
        itypfb[ifbr] = inlet_type;
      }

      if (cs_gui_strcmp(vars->model, "atmospheric_flows")) {
        iprofm[zone_nbr-1] = boundaries->meteo[izone].read_data;
        if (iprofm[zone_nbr-1] == 1) {
          BFT_FREE(choice_v);
          BFT_FREE(choice_d);
        }
        if (boundaries->meteo[izone].automatic) {
          for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
            ifbr = face_ids[ifac];
            iautom[ifbr] = 1;
          }
        }
      }

      /* Dirichlet for velocity */

      const cs_field_t  *fv = cs_field_by_name_try("velocity");
      const int var_key_id = cs_field_key_id("variable_id");
      int ivarv = cs_field_get_key_int(fv, var_key_id) -1;
      if (cs_gui_strcmp(choice_d, "coordinates")) {
        if (cs_gui_strcmp(choice_v, "norm")) {
          norm = boundaries->norm[izone] /
              sqrt(  boundaries->dirx[izone] * boundaries->dirx[izone]
                   + boundaries->diry[izone] * boundaries->diry[izone]
                   + boundaries->dirz[izone] * boundaries->dirz[izone]);

          for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
            ifbr = face_ids[ifac];
            rcodcl[ ivarv      * n_b_faces + ifbr] = boundaries->dirx[izone] * norm;
            rcodcl[(ivarv + 1) * n_b_faces + ifbr] = boundaries->diry[izone] * norm;
            rcodcl[(ivarv + 2) * n_b_faces + ifbr] = boundaries->dirz[izone] * norm;
          }
        }
        else if (cs_gui_strcmp(choice_v, "flow1")
            || cs_gui_strcmp(choice_v, "flow2")
            || cs_gui_strcmp(choice_v, "flow1_formula")
            || cs_gui_strcmp(choice_v, "flow2_formula")) {
          for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
            ifbr = face_ids[ifac];
            rcodcl[ ivarv      * n_b_faces + ifbr] = boundaries->dirx[izone];
            rcodcl[(ivarv + 1) * n_b_faces + ifbr] = boundaries->diry[izone];
            rcodcl[(ivarv + 2) * n_b_faces + ifbr] = boundaries->dirz[izone];
          }
        }
        else if (cs_gui_strcmp(choice_v, "norm_formula")) {
          t0 = cs_timer_wtime();

          mei_tree_insert(boundaries->velocity[izone], "t", ts->t_cur);
          mei_tree_insert(boundaries->velocity[izone], "dt",
                          cs_glob_time_step_options->dtref);
          mei_tree_insert(boundaries->velocity[izone], "iter", ts->nt_cur);

          /* add variable from notebook */
          _add_notebook_variables(boundaries->velocity[izone]);

          for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
            ifbr = face_ids[ifac];

            mei_tree_insert(boundaries->velocity[izone], "x", cdgfbo[3 * ifbr + 0]);
            mei_tree_insert(boundaries->velocity[izone], "y", cdgfbo[3 * ifbr + 1]);
            mei_tree_insert(boundaries->velocity[izone], "z", cdgfbo[3 * ifbr + 2]);

            cs_real_t dir[3] = {boundaries->dirx[izone],
                                boundaries->diry[izone],
                                boundaries->dirz[izone]};

            cs_real_t x_norm = cs_math_3_norm(dir);
            if (x_norm <= 0.)
              bft_error(__FILE__, __LINE__, 0,
                  _("Error in the GUI boundary conditions: the normal direction is of norm 0.\n "));

            mei_evaluate(boundaries->velocity[izone]);

            norm = mei_tree_lookup(boundaries->velocity[izone], "u_norm")
              / x_norm;

            rcodcl[ ivarv      * n_b_faces + ifbr] = boundaries->dirx[izone] * norm;
            rcodcl[(ivarv + 1) * n_b_faces + ifbr] = boundaries->diry[izone] * norm;
            rcodcl[(ivarv + 2) * n_b_faces + ifbr] = boundaries->dirz[izone] * norm;
          }

          cs_gui_add_mei_time(cs_timer_wtime() - t0);

        }
        if (cs_gui_strcmp(vars->model, "compressible_model")) {
          if (boundaries->itype[izone] == CS_EPHCF) {
            for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
              ifbr = face_ids[ifac];
              rcodcl[ ivarv      * n_b_faces + ifbr] = boundaries->dirx[izone];
              rcodcl[(ivarv + 1) * n_b_faces + ifbr] = boundaries->diry[izone];
              rcodcl[(ivarv + 2) * n_b_faces + ifbr] = boundaries->dirz[izone];
            }
          }
        }
      }
      else if (   cs_gui_strcmp(choice_d, "normal")
               || cs_gui_strcmp(choice_d, "translation")) {
        if (cs_gui_strcmp(choice_v, "norm")) {
          for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
            ifbr = face_ids[ifac];

            norm =   boundaries->norm[izone]
            / sqrt(  surfbo[3 * ifbr + 0] * surfbo[3 * ifbr + 0]
                   + surfbo[3 * ifbr + 1] * surfbo[3 * ifbr + 1]
                   + surfbo[3 * ifbr + 2] * surfbo[3 * ifbr + 2]);

            for (i = 0; i < 3; i++)
              rcodcl[(ivarv + i) * n_b_faces + ifbr] = -surfbo[3 * ifbr + i] * norm;
          }
        }
        else if (   cs_gui_strcmp(choice_v, "flow1")
                 || cs_gui_strcmp(choice_v, "flow2")
                 || cs_gui_strcmp(choice_v, "flow1_formula")
                 || cs_gui_strcmp(choice_v, "flow2_formula")) {
          for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
            ifbr = face_ids[ifac];
            norm = sqrt(  surfbo[3 * ifbr + 0] * surfbo[3 * ifbr + 0]
                        + surfbo[3 * ifbr + 1] * surfbo[3 * ifbr + 1]
                        + surfbo[3 * ifbr + 2] * surfbo[3 * ifbr + 2]);

            for (i = 0; i < 3; i++)
              rcodcl[(ivarv + i) * n_b_faces + ifbr] = -surfbo[3 * ifbr + i] / norm;
          }
        }
        else if (cs_gui_strcmp(choice_v, "norm_formula")) {
          t0 = cs_timer_wtime();

          mei_tree_insert(boundaries->velocity[izone], "t", ts->t_cur);
          mei_tree_insert(boundaries->velocity[izone], "dt",
                          cs_glob_time_step_options->dtref);
          mei_tree_insert(boundaries->velocity[izone], "iter", ts->nt_cur);

          /* add variable from notebook */
          _add_notebook_variables(boundaries->velocity[izone]);

          for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
            ifbr = face_ids[ifac];

            mei_tree_insert(boundaries->velocity[izone], "x", cdgfbo[3 * ifbr + 0]);
            mei_tree_insert(boundaries->velocity[izone], "y", cdgfbo[3 * ifbr + 1]);
            mei_tree_insert(boundaries->velocity[izone], "z", cdgfbo[3 * ifbr + 2]);

            mei_evaluate(boundaries->velocity[izone]);

            norm =   mei_tree_lookup(boundaries->velocity[izone], "u_norm")
            / sqrt(  surfbo[3 * ifbr + 0] * surfbo[3 * ifbr + 0]
                   + surfbo[3 * ifbr + 1] * surfbo[3 * ifbr + 1]
                   + surfbo[3 * ifbr + 2] * surfbo[3 * ifbr + 2]);

            for (i = 0; i < 3; i++)
              rcodcl[(ivarv + i) * n_b_faces + ifbr] = -surfbo[3 * ifbr + i] * norm;
          }
          cs_gui_add_mei_time(cs_timer_wtime() - t0);
        }

        if (cs_gui_strcmp(vars->model, "compressible_model")) {
          if (boundaries->itype[izone] == CS_EPHCF) {
            for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
              ifbr = face_ids[ifac];

              for (i = 0; i < 3; i++)
                rcodcl[(ivarv + i) * n_b_faces + ifbr] = -surfbo[3 * ifbr + i];
            }
          }
        }
      }
      else if (cs_gui_strcmp(choice_d, "formula")) {
        t0 = cs_timer_wtime();

        mei_tree_insert(boundaries->direction[izone], "t", ts->t_cur);
        mei_tree_insert(boundaries->direction[izone], "dt",
                        cs_glob_time_step_options->dtref);
        mei_tree_insert(boundaries->direction[izone], "iter", ts->nt_cur);

          /* add variable from notebook */
          _add_notebook_variables(boundaries->direction[izone]);

        if (cs_gui_strcmp(choice_v, "norm")) {
          for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
            ifbr = face_ids[ifac];

            mei_tree_insert(boundaries->direction[izone], "x", cdgfbo[3 * ifbr + 0]);
            mei_tree_insert(boundaries->direction[izone], "y", cdgfbo[3 * ifbr + 1]);
            mei_tree_insert(boundaries->direction[izone], "z", cdgfbo[3 * ifbr + 2]);

            mei_evaluate(boundaries->direction[izone]);

            X[0] = mei_tree_lookup(boundaries->direction[izone], "dir_x");
            X[1] = mei_tree_lookup(boundaries->direction[izone], "dir_y");
            X[2] = mei_tree_lookup(boundaries->direction[izone], "dir_z");

            cs_real_t x_norm = cs_math_3_norm(X);
            if (x_norm <= 0.)
              bft_error(__FILE__, __LINE__, 0,
                  _("Error in the GUI boundary conditions: the normal direction is of norm 0.\n "));

            norm = boundaries->norm[izone] / x_norm;

            for (i = 0; i < 3; i++)
              rcodcl[(ivarv + i) * n_b_faces + ifbr] = X[i] * norm;
          }
        }
        else if (   cs_gui_strcmp(choice_v, "flow1")
                 || cs_gui_strcmp(choice_v, "flow2")
                 || cs_gui_strcmp(choice_v, "flow1_formula")
                 || cs_gui_strcmp(choice_v, "flow2_formula")) {
          for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
            ifbr = face_ids[ifac];

            mei_tree_insert(boundaries->direction[izone], "x", cdgfbo[3 * ifbr + 0]);
            mei_tree_insert(boundaries->direction[izone], "y", cdgfbo[3 * ifbr + 1]);
            mei_tree_insert(boundaries->direction[izone], "z", cdgfbo[3 * ifbr + 2]);

            mei_evaluate(boundaries->direction[izone]);

            X[0] = mei_tree_lookup(boundaries->direction[izone], "dir_x");
            X[1] = mei_tree_lookup(boundaries->direction[izone], "dir_y");
            X[2] = mei_tree_lookup(boundaries->direction[izone], "dir_z");

            cs_real_t x_norm = cs_math_3_norm(X);
            if (x_norm <= 0.)
              bft_error(__FILE__, __LINE__, 0,
                  _("Error in the GUI boundary conditions: the normal direction is of norm 0.\n "));

            for (i = 0; i < 3; i++)
              rcodcl[(ivarv + i) * n_b_faces + ifbr] = X[i];
          }
        }
        else if (cs_gui_strcmp(choice_v, "norm_formula")) {
          mei_tree_insert(boundaries->velocity[izone], "t", ts->t_cur);
          mei_tree_insert(boundaries->velocity[izone], "dt",
                          cs_glob_time_step_options->dtref);
          mei_tree_insert(boundaries->velocity[izone], "iter", ts->nt_cur);

          for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
            ifbr = face_ids[ifac];

            mei_tree_insert(boundaries->velocity[izone], "x", cdgfbo[3 * ifbr + 0]);
            mei_tree_insert(boundaries->velocity[izone], "y", cdgfbo[3 * ifbr + 1]);
            mei_tree_insert(boundaries->velocity[izone], "z", cdgfbo[3 * ifbr + 2]);

            mei_tree_insert(boundaries->direction[izone], "x", cdgfbo[3 * ifbr + 0]);
            mei_tree_insert(boundaries->direction[izone], "y", cdgfbo[3 * ifbr + 1]);
            mei_tree_insert(boundaries->direction[izone], "z", cdgfbo[3 * ifbr + 2]);

            mei_evaluate(boundaries->velocity[izone]);
            mei_evaluate(boundaries->direction[izone]);

            X[0] = mei_tree_lookup(boundaries->direction[izone], "dir_x");
            X[1] = mei_tree_lookup(boundaries->direction[izone], "dir_y");
            X[2] = mei_tree_lookup(boundaries->direction[izone], "dir_z");

            cs_real_t x_norm = cs_math_3_norm(X);
            if (x_norm <= 0.)
              bft_error(__FILE__, __LINE__, 0,
                  _("Error in the GUI boundary conditions: the normal direction is of norm 0.\n "));

            norm = mei_tree_lookup(boundaries->velocity[izone], "u_norm")
              / x_norm;

            for (i = 0; i < 3; i++)
              rcodcl[(ivarv + i) * n_b_faces + ifbr] = X[i] * norm;
          }
        }

        if (cs_gui_strcmp(vars->model, "compressible_model")) {
          if (boundaries->itype[izone] == CS_EPHCF) {
            for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
              ifbr = face_ids[ifac];

              mei_tree_insert(boundaries->direction[izone], "x", cdgfbo[3 * ifbr + 0]);
              mei_tree_insert(boundaries->direction[izone], "y", cdgfbo[3 * ifbr + 1]);
              mei_tree_insert(boundaries->direction[izone], "z", cdgfbo[3 * ifbr + 2]);

              mei_evaluate(boundaries->direction[izone]);

              X[0] = mei_tree_lookup(boundaries->direction[izone], "dir_x");
              X[1] = mei_tree_lookup(boundaries->direction[izone], "dir_y");
              X[2] = mei_tree_lookup(boundaries->direction[izone], "dir_z");

              for (i = 0; i < 3; i++)
                  rcodcl[(ivarv + i) * n_b_faces + ifbr] = X[i];
            }
          }
        }
        cs_gui_add_mei_time(cs_timer_wtime() - t0);

      }

      if (cs_gui_strcmp(vars->model, "groundwater_model")) {
        for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
          ifbr = face_ids[ifac];
          for (i = 0; i < 3; i++)
          {
            icodcl[(ivarv + i)* n_b_faces + ifbr] = 3;
            rcodcl[(ivarv + i)* n_b_faces + ifbr] = 0.;
          }
        }
      }
      BFT_FREE(choice_v);
      BFT_FREE(choice_d);

      if (cs_gui_strcmp(vars->model, "groundwater_model")) {
        const cs_field_t  *fp1 = cs_field_by_name_try("hydraulic_head");
        int ivar1 = cs_field_get_key_int(fp1, var_key_id) -1;
        char *_choice_d = _boundary_choice(boundaries->nature[izone],
                                           boundaries->label[izone],
                                           "hydraulicHead", "choice");

        if (cs_gui_strcmp(_choice_d, "dirichlet")) {
          for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
            ifbr = face_ids[ifac];
            icodcl[ivar1 * n_b_faces + ifbr] = 1;
            rcodcl[ivar1 * n_b_faces + ifbr] = boundaries->preout[izone];
          }
        }
        else if (cs_gui_strcmp(_choice_d, "neumann")) {
          for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
            ifbr = face_ids[ifac];
            icodcl[ivar1 * n_b_faces + ifbr] = 3;
            rcodcl[2 * n_b_faces * (*nvar) + ivar1 * n_b_faces + ifbr]
              = boundaries->preout[izone];
          }
        }
        else if (cs_gui_strcmp(_choice_d, "dirichlet_formula")) {
          for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
            ifbr = face_ids[ifac];
            icodcl[ivar1 * n_b_faces + ifbr] = 1;

            mei_tree_t *ev_formula = boundaries->groundwat[izone];
            mei_tree_insert(ev_formula, "t", ts->t_cur);
            mei_tree_insert(ev_formula, "dt", cs_glob_time_step_options->dtref);
            mei_tree_insert(ev_formula, "iter", ts->nt_cur);
            mei_tree_insert(ev_formula, "x", cdgfbo[3 * ifbr + 0]);
            mei_tree_insert(ev_formula, "y", cdgfbo[3 * ifbr + 1]);
            mei_tree_insert(ev_formula, "z", cdgfbo[3 * ifbr + 2]);

            /* add variable from notebook */
            _add_notebook_variables(ev_formula);

            mei_evaluate(ev_formula);
            rcodcl[ivar1 * n_b_faces + ifbr]
              = mei_tree_lookup(ev_formula, "H");
          }
        }
        BFT_FREE(_choice_d);
      }

      /* turbulent inlet, with formula */
      if (icalke[zone_nbr-1] == 0) {
        t0 = cs_timer_wtime();
        char *path1 = cs_xpath_init_path();
        cs_xpath_add_elements(&path1, 2, "boundary_conditions", "inlet");
        cs_xpath_add_test_attribute(&path1, "label", boundaries->label[izone]);
        cs_xpath_add_element(&path1, "turbulence");
        cs_xpath_add_element(&path1, "formula");
        cs_xpath_add_function_text(&path1);
        char *formula = cs_gui_get_text_value(path1);
        if (formula != NULL) {
          mei_tree_t *ev_formula = mei_tree_new(formula);
          mei_tree_insert(ev_formula,"x",0.0);
          mei_tree_insert(ev_formula,"y",0.0);
          mei_tree_insert(ev_formula,"z",0.0);
          mei_tree_insert(ev_formula, "t", ts->t_cur);
          mei_tree_insert(ev_formula, "dt", cs_glob_time_step_options->dtref);
          mei_tree_insert(ev_formula, "iter", ts->nt_cur);

          /* add variable from notebook */
          _add_notebook_variables(ev_formula);

          /* try to build the interpreter */

          if (mei_tree_builder(ev_formula))
            bft_error(__FILE__, __LINE__, 0,
                      _("Error: can not interpret expression: %s\n %i"),
                      ev_formula->string, mei_tree_builder(ev_formula));

          char *model = cs_gui_get_thermophysical_model("turbulence");
          if (model == NULL)
            return;

          if (   cs_gui_strcmp(model, "k-epsilon")
              || cs_gui_strcmp(model, "k-epsilon-PL")) {
            const char *symbols[] = {"k","epsilon"};
            if (mei_tree_find_symbols(ev_formula, 2, symbols))
              bft_error(__FILE__, __LINE__, 0,
                        _("Error: can not find the required symbol: %s\n"),
                        "k or epsilon");

            cs_field_t *c_k   = cs_field_by_name("k");
            cs_field_t *c_eps = cs_field_by_name("epsilon");
            int ivark = cs_field_get_key_int(c_k, var_key_id) -1;
            int ivare = cs_field_get_key_int(c_eps, var_key_id) -1;

            for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
              ifbr = face_ids[ifac];
              mei_tree_insert(ev_formula, "x", cdgfbo[3 * ifbr + 0]);
              mei_tree_insert(ev_formula, "y", cdgfbo[3 * ifbr + 1]);
              mei_tree_insert(ev_formula, "z", cdgfbo[3 * ifbr + 2]);
              mei_evaluate(ev_formula);
              rcodcl[ivark * n_b_faces + ifbr] = mei_tree_lookup(ev_formula, "k");
              rcodcl[ivare * n_b_faces + ifbr] = mei_tree_lookup(ev_formula, "epsilon");
            }
          }
          else if (  cs_gui_strcmp(model, "Rij-epsilon")
                   ||cs_gui_strcmp(model, "Rij-SSG")) {
            const char *symbols[] = {"r11", "r22", "r33", "r12", "r13", "r23", "epsilon"};
            if (mei_tree_find_symbols(ev_formula, 7, symbols))
              bft_error(__FILE__, __LINE__, 0,
                        _("Error: can not find the required symbol: %s\n"),
                        "r11, r22, r33, r12, r13, r23 or epsilon");

            cs_field_t *cfld_rij;
            if (cs_glob_turb_rans_model->irijco == 1)
              cfld_rij = cs_field_by_name("rij");
            else
              cfld_rij = cs_field_by_name("r11");

            cs_field_t *c_eps = cs_field_by_name("epsilon");
            int ivarrij = cs_field_get_key_int(cfld_rij, var_key_id) - 1;
            int ivare   = cs_field_get_key_int(c_eps, var_key_id) -1;

            for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
              ifbr = face_ids[ifac];
              mei_tree_insert(ev_formula, "x", cdgfbo[3 * ifbr + 0]);
              mei_tree_insert(ev_formula, "y", cdgfbo[3 * ifbr + 1]);
              mei_tree_insert(ev_formula, "z", cdgfbo[3 * ifbr + 2]);
              mei_evaluate(ev_formula);
              rcodcl[ ivarrij      * n_b_faces + ifbr]  = mei_tree_lookup(ev_formula, "r11");
              rcodcl[(ivarrij + 1) * n_b_faces + ifbr]  = mei_tree_lookup(ev_formula, "r22");
              rcodcl[(ivarrij + 2) * n_b_faces + ifbr]  = mei_tree_lookup(ev_formula, "r33");
              rcodcl[(ivarrij + 3) * n_b_faces + ifbr]  = mei_tree_lookup(ev_formula, "r12");
              rcodcl[(ivarrij + 4) * n_b_faces + ifbr]  = mei_tree_lookup(ev_formula, "r23");
              rcodcl[(ivarrij + 5) * n_b_faces + ifbr]  = mei_tree_lookup(ev_formula, "r13");
              rcodcl[ivare   * n_b_faces + ifbr] = mei_tree_lookup(ev_formula, "epsilon");
            }
          }
          else if (cs_gui_strcmp(model, "Rij-EBRSM")) {
            const char *symbols[] = {"r11", "r22", "r33", "r12", "r13", "r23",
                                     "epsilon", "alpha"};
            if (mei_tree_find_symbols(ev_formula, 8, symbols))
              bft_error(__FILE__, __LINE__, 0,
                        _("Error: can not find the required symbol: %s\n"),
                        "R11, R22, R33, R12, R13, R23, eps or alpha");

            cs_field_t *cfld_rij;
            if (cs_glob_turb_rans_model->irijco == 1)
              cfld_rij = cs_field_by_name("rij");
            else
              cfld_rij = cs_field_by_name("r11");
            cs_field_t *c_eps = cs_field_by_name("epsilon");
            cs_field_t *c_a   = cs_field_by_name("alpha");

            int ivarrij = cs_field_get_key_int(cfld_rij, var_key_id) - 1;
            int ivare   = cs_field_get_key_int(c_eps, var_key_id) -1;
            int ivara   = cs_field_get_key_int(c_a, var_key_id) -1;

            for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
              ifbr = face_ids[ifac];
              mei_tree_insert(ev_formula, "x", cdgfbo[3 * ifbr + 0]);
              mei_tree_insert(ev_formula, "y", cdgfbo[3 * ifbr + 1]);
              mei_tree_insert(ev_formula, "z", cdgfbo[3 * ifbr + 2]);
              mei_evaluate(ev_formula);
              rcodcl[ ivarrij      * n_b_faces + ifbr]  = mei_tree_lookup(ev_formula, "r11");
              rcodcl[(ivarrij + 1) * n_b_faces + ifbr]  = mei_tree_lookup(ev_formula, "r22");
              rcodcl[(ivarrij + 2) * n_b_faces + ifbr]  = mei_tree_lookup(ev_formula, "r33");
              rcodcl[(ivarrij + 3) * n_b_faces + ifbr]  = mei_tree_lookup(ev_formula, "r12");
              rcodcl[(ivarrij + 4) * n_b_faces + ifbr]  = mei_tree_lookup(ev_formula, "r23");
              rcodcl[(ivarrij + 5) * n_b_faces + ifbr]  = mei_tree_lookup(ev_formula, "r13");
              rcodcl[ivare   * n_b_faces + ifbr] = mei_tree_lookup(ev_formula, "epsilon");
              rcodcl[ivara   * n_b_faces + ifbr] = mei_tree_lookup(ev_formula, "alpha");
            }
          }
          else if (cs_gui_strcmp(model, "v2f-BL-v2/k")) {
            const char *symbols[] = {"k", "epsilon", "phi", "alpha"};

            if (mei_tree_find_symbols(ev_formula, 4, symbols))
              bft_error(__FILE__, __LINE__, 0,
                        _("Error: can not find the required symbol: %s\n"),
                        "k, eps, phi of alpha");

            cs_field_t *c_k   = cs_field_by_name("k");
            cs_field_t *c_eps = cs_field_by_name("epsilon");
            cs_field_t *c_phi = cs_field_by_name("phi");
            cs_field_t *c_a   = cs_field_by_name("alpha");
            int ivark = cs_field_get_key_int(c_k,   var_key_id) -1;
            int ivare = cs_field_get_key_int(c_eps, var_key_id) -1;
            int ivarp = cs_field_get_key_int(c_phi, var_key_id) -1;
            int ivara = cs_field_get_key_int(c_a,   var_key_id) -1;

            for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
              ifbr = face_ids[ifac];
              mei_tree_insert(ev_formula, "x", cdgfbo[3 * ifbr + 0]);
              mei_tree_insert(ev_formula, "y", cdgfbo[3 * ifbr + 1]);
              mei_tree_insert(ev_formula, "z", cdgfbo[3 * ifbr + 2]);
              mei_evaluate(ev_formula);
              rcodcl[ivark * n_b_faces + ifbr] = mei_tree_lookup(ev_formula, "k");
              rcodcl[ivare * n_b_faces + ifbr] = mei_tree_lookup(ev_formula, "epsilon");
              rcodcl[ivarp * n_b_faces + ifbr] = mei_tree_lookup(ev_formula, "phi");
              rcodcl[ivara * n_b_faces + ifbr] = mei_tree_lookup(ev_formula, "alpha");
            }
          }
          else if (cs_gui_strcmp(model, "k-omega-SST")) {
            const char *symbols[] = {"k", "omega"};

            if (mei_tree_find_symbols(ev_formula, 2, symbols))
              bft_error(__FILE__, __LINE__, 0,
                        _("Error: can not find the required symbol: %s\n"),
                        "k or omega");

            cs_field_t *c_k = cs_field_by_name("k");
            cs_field_t *c_o = cs_field_by_name("omega");
            int ivark = cs_field_get_key_int(c_k,   var_key_id) -1;
            int ivaro = cs_field_get_key_int(c_o,   var_key_id) -1;

            for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
              ifbr = face_ids[ifac];
              mei_tree_insert(ev_formula, "x", cdgfbo[3 * ifbr + 0]);
              mei_tree_insert(ev_formula, "y", cdgfbo[3 * ifbr + 1]);
              mei_tree_insert(ev_formula, "z", cdgfbo[3 * ifbr + 2]);
              mei_evaluate(ev_formula);
              rcodcl[ivark * n_b_faces + ifbr] = mei_tree_lookup(ev_formula, "k");
              rcodcl[ivaro * n_b_faces + ifbr] = mei_tree_lookup(ev_formula, "omega");
            }
          }
          else if (cs_gui_strcmp(model, "Spalart-Allmaras")) {
            const char *symbols[] = {"nu_tilda"};

            if (mei_tree_find_symbols(ev_formula, 1, symbols))
              bft_error(__FILE__, __LINE__, 0,
                        _("Error: can not find the required symbol: %s\n"),
                        "nu_tilda");

            cs_field_t *c_nu = cs_field_by_name("nu_tilda");
            int ivarnu = cs_field_get_key_int(c_nu, var_key_id) -1;

            for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
              ifbr = face_ids[ifac];
              mei_tree_insert(ev_formula, "x", cdgfbo[3 * ifbr + 0]);
              mei_tree_insert(ev_formula, "y", cdgfbo[3 * ifbr + 1]);
              mei_tree_insert(ev_formula, "z", cdgfbo[3 * ifbr + 2]);
              mei_evaluate(ev_formula);
              rcodcl[ivarnu * n_b_faces + ifbr] = mei_tree_lookup(ev_formula, "nu_tilda");
            }
          }
          else
            bft_error(__FILE__, __LINE__, 0,
                      _("Invalid turbulence model: %s.\n"), model);
          mei_tree_destroy(ev_formula);
          BFT_FREE(formula);
          BFT_FREE(model);
        }
        BFT_FREE(path1);
        cs_gui_add_mei_time(cs_timer_wtime() - t0);
      }
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "wall")) {
      int iwall = CS_SMOOTHWALL;

      if (boundaries->rough[izone] >= 0.0) {
        const cs_field_t  *fv = cs_field_by_name_try("velocity");
        const int var_key_id = cs_field_key_id("variable_id");
        int ivarv = cs_field_get_key_int(fv, var_key_id) -1;

        iwall = CS_ROUGHWALL;

        /* Roughness value is stored in Velocity_U (z0) */
        /* Remember: rcodcl(ifac, ivar, 1) -> rcodcl[k][j][i]
           = rcodcl[ k*dim1*dim2 + j*dim1 + i] */
        for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
          ifbr = face_ids[ifac];
          cs_lnum_t idx = 2 * n_b_faces * (*nvar) + ivarv * n_b_faces + ifbr;
          rcodcl[idx] = boundaries->rough[izone];

          /* Roughness value is also stored in Velocity_V for eventual scalar
           * (even if there is no scalar). In this case rugd = rugt. */
          cs_lnum_t idx2 = 2 * n_b_faces * (*nvar)
                         + (ivarv + 1) * n_b_faces + ifbr;
          rcodcl[idx2] = boundaries->rough[izone];
        }
      }

      for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
        ifbr = face_ids[ifac];
        izfppp[ifbr] = zone_nbr;
        itypfb[ifbr] = iwall;
      }

    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "outlet")) {
      for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
        ifbr = face_ids[ifac];
        izfppp[ifbr] = zone_nbr;
        if (cs_gui_strcmp(vars->model, "compressible_model"))
          itypfb[ifbr] = boundaries->itype[izone];
        else
          itypfb[ifbr] = CS_OUTLET;
      }

      if (cs_gui_strcmp(vars->model, "atmospheric_flows")) {
        iprofm[zone_nbr-1] = boundaries->meteo[izone].read_data;
        if (boundaries->meteo[izone].automatic) {
          for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
            ifbr = face_ids[ifac];
            itypfb[ifbr] = 0;
          }
        }
      }
      else if (cs_gui_strcmp(vars->model, "compressible_model")) {
        if (boundaries->itype[izone] == CS_SOPCF) {
          const cs_field_t  *fp1 = cs_field_by_name_try("pressure");
          const int var_key_id = cs_field_key_id("variable_id");
          int ivar1 = cs_field_get_key_int(fp1, var_key_id) -1;
          for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
            ifbr = face_ids[ifac];
            rcodcl[ivar1 * n_b_faces + ifbr] = boundaries->preout[izone];
          }
        }
      }
      else if (cs_gui_strcmp(vars->model, "groundwater_model")) {
        const cs_field_t  *fp1 = cs_field_by_name_try("hydraulic_head");
        const int var_key_id = cs_field_key_id("variable_id");
        int ivar1 = cs_field_get_key_int(fp1, var_key_id) -1;
        char *choice_d = _boundary_choice(boundaries->nature[izone],
                                          boundaries->label[izone],
                                          "hydraulicHead", "choice");

        if (cs_gui_strcmp(choice_d, "dirichlet")) {
          for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
            ifbr = face_ids[ifac];
            icodcl[ivar1 * n_b_faces + ifbr] = 1;
            rcodcl[ivar1 * n_b_faces + ifbr] = boundaries->preout[izone];
          }
        }
        else if (cs_gui_strcmp(choice_d, "neumann")) {
          for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
            ifbr = face_ids[ifac];
            icodcl[ivar1 * n_b_faces + ifbr] = 3;
            rcodcl[2 * n_b_faces * (*nvar) + ivar1 * n_b_faces + ifbr]
              = boundaries->preout[izone];
          }
        }
        else if (cs_gui_strcmp(choice_d, "dirichlet_formula")) {
          for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
            ifbr = face_ids[ifac];
            icodcl[ivar1 * n_b_faces + ifbr] = 1;

            mei_tree_t *ev_formula = boundaries->groundwat[izone];
            mei_tree_insert(ev_formula, "t", ts->t_cur);
            mei_tree_insert(ev_formula, "dt", cs_glob_time_step_options->dtref);
            mei_tree_insert(ev_formula, "iter", ts->nt_cur);
            mei_tree_insert(ev_formula, "x", cdgfbo[3 * ifbr + 0]);
            mei_tree_insert(ev_formula, "y", cdgfbo[3 * ifbr + 1]);
            mei_tree_insert(ev_formula, "z", cdgfbo[3 * ifbr + 2]);

            /* add variable from notebook */
            _add_notebook_variables(ev_formula);

            mei_evaluate(ev_formula);
            rcodcl[ivar1 * n_b_faces + ifbr]
              = mei_tree_lookup(ev_formula, "H");

          }
        }
        BFT_FREE(choice_d);
      }
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "imposed_p_outlet")) {
      for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
        ifbr = face_ids[ifac];
        izfppp[ifbr] = zone_nbr;
        itypfb[ifbr] = CS_OUTLET;
      }

      /* imposed outlet pressure */
      const cs_field_t  *fp1 = cs_field_by_name_try("pressure");
      const int var_key_id = cs_field_key_id("variable_id");
      int ivar1 = cs_field_get_key_int(fp1, var_key_id) -1;
      char *_choice_d = _boundary_choice(boundaries->nature[izone],
                                         boundaries->label[izone],
                                         "pressure", "choice");

      if (cs_gui_strcmp(_choice_d, "dirichlet")) {
        for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
          ifbr = face_ids[ifac];
          icodcl[ivar1 * n_b_faces + ifbr] = 1;
          rcodcl[ivar1 * n_b_faces + ifbr] = boundaries->preout[izone];
        }
      }

    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "symmetry")) {
      for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
        ifbr = face_ids[ifac];
        izfppp[ifbr] = zone_nbr;
        itypfb[ifbr] = CS_SYMMETRY;
      }
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "free_inlet_outlet")) {
      for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
        ifbr = face_ids[ifac];
        izfppp[ifbr] = zone_nbr;
        itypfb[ifbr] = CS_FREE_INLET;
      }

      if (boundaries->headLoss[izone] != NULL) {
        t0 = cs_timer_wtime();

        const cs_field_t  *fp = cs_field_by_name_try("pressure");
        const int var_key_id = cs_field_key_id("variable_id");
        int ivarp = cs_field_get_key_int(fp, var_key_id) -1;

        mei_tree_insert(boundaries->headLoss[izone], "t", ts->t_cur);
        mei_tree_insert(boundaries->headLoss[izone], "dt", cs_glob_time_step_options->dtref);
        mei_tree_insert(boundaries->headLoss[izone], "iter", ts->nt_cur);

        /* add variable from notebook */
        _add_notebook_variables(boundaries->headLoss[izone]);

        for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
          ifbr = face_ids[ifac];

          mei_tree_insert(boundaries->headLoss[izone], "x", cdgfbo[3 * ifbr + 0]);
          mei_tree_insert(boundaries->headLoss[izone], "y", cdgfbo[3 * ifbr + 1]);
          mei_tree_insert(boundaries->headLoss[izone], "z", cdgfbo[3 * ifbr + 2]);

          mei_evaluate(boundaries->headLoss[izone]);
          rcodcl[1 * n_b_faces * (*nvar) + ivarp * n_b_faces + ifbr] =
              mei_tree_lookup(boundaries->headLoss[izone], "K");
        }
        cs_gui_add_mei_time(cs_timer_wtime() - t0);
      }
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "free_surface")) {
      for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
        ifbr = face_ids[ifac];
        izfppp[ifbr] = zone_nbr;
        itypfb[ifbr] = CS_FREE_SURFACE;
      }
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "groundwater")) {
      const cs_field_t  *fp1 = cs_field_by_name_try("hydraulic_head");
      const int var_key_id = cs_field_key_id("variable_id");
      int ivar1 = cs_field_get_key_int(fp1, var_key_id) -1;
      char *choice_d = _boundary_choice(boundaries->nature[izone],
                                        boundaries->label[izone],
                                        "hydraulicHead",
                                        "choice");

      for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
        ifbr = face_ids[ifac];
        izfppp[ifbr] = zone_nbr;
        itypfb[ifbr] = CS_INDEF;
      }

      // TODO set velocity to 0
      const cs_field_t  *fp2 = cs_field_by_name_try("velocity");
      int ivar2 = cs_field_get_key_int(fp2, var_key_id) -1;

      for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
        ifbr = face_ids[ifac];
        for (i = 0; i < 3; i++) {
          icodcl[(ivar2 + i) * n_b_faces + ifbr] = 3;
          rcodcl[(ivar2 + i) * n_b_faces + ifbr] = 0.;
        }
      }

      if (cs_gui_strcmp(choice_d, "dirichlet"))
      {
        for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
          ifbr = face_ids[ifac];
          icodcl[ivar1 * n_b_faces + ifbr] = 1;
          rcodcl[ivar1 * n_b_faces + ifbr] = boundaries->preout[izone];
        }
      }
      else if (cs_gui_strcmp(choice_d, "neumann"))
      {
        for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
          ifbr = face_ids[ifac];
          icodcl[ivar1 * n_b_faces + ifbr] = 3;
          rcodcl[2 * n_b_faces * (*nvar) + ivar1 * n_b_faces + ifbr]
            = boundaries->preout[izone];
        }
      }
      else if (cs_gui_strcmp(choice_d, "dirichlet_formula"))
      {
        for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
          ifbr = face_ids[ifac];
          icodcl[ivar1 * n_b_faces + ifbr] = 1;

          mei_tree_t *ev_formula = boundaries->groundwat[izone];
          mei_tree_insert(ev_formula, "t", ts->t_cur);
          mei_tree_insert(ev_formula, "dt", cs_glob_time_step_options->dtref);
          mei_tree_insert(ev_formula, "iter", ts->nt_cur);
          mei_tree_insert(ev_formula, "x", cdgfbo[3 * ifbr + 0]);
          mei_tree_insert(ev_formula, "y", cdgfbo[3 * ifbr + 1]);
          mei_tree_insert(ev_formula, "z", cdgfbo[3 * ifbr + 2]);

          /* add variable from notebook */
          _add_notebook_variables(ev_formula);

          mei_evaluate(ev_formula);
          rcodcl[ivar1 * n_b_faces + ifbr]
            = mei_tree_lookup(ev_formula, "H");

        }
      }
      BFT_FREE(choice_d);
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "undefined")) {
      for (cs_lnum_t ifac = 0; ifac < faces; ifac++) {
        ifbr = face_ids[ifac];
        izfppp[ifbr] = zone_nbr;
        itypfb[ifbr] = CS_INDEF;
      }

    }
    else {
      bft_error(__FILE__, __LINE__, 0,
          _("boundary nature %s is unknown \n"),
          boundaries->nature[izone]);
    }

    /* treatment of mapped inlet for each field */

    if (boundaries->locator[izone] != NULL && ts->nt_cur > 1) {
      icalke[zone_nbr-1] = 0;

      for (int f_id = 0; f_id < n_fields; f_id++) {
        const cs_field_t  *f = cs_field_by_id(f_id);
        const int var_key_id = cs_field_key_id("variable_id");
        ivar = cs_field_get_key_int(f, var_key_id) -1;

        if (f->type & CS_FIELD_VARIABLE) {
          int interpolate = 0;
          int normalize = 0;
          if (f == CS_F_(u))
            normalize = 1;
          else {
            const int keysca = cs_field_key_id("scalar_id");
            if (cs_field_get_key_int(f, keysca) > 0)
              normalize = 1;
          }
          if (f != CS_F_(p))
            cs_boundary_conditions_mapped_set(f, boundaries->locator[izone],
                                              CS_MESH_LOCATION_CELLS,
                                              normalize, interpolate,
                                              faces, face_ids,
                                              NULL, *nvar, rcodcl);
        }

      }
    }

  } /*  for (izone=0 ; izone < zones ; izone++) */

#if _XML_DEBUG_
  bft_printf("==>UICLIM\n");
  bft_printf("--boundary zones number: %i\n", zones);

  for (izone = 0 ; izone < zones ; izone++) {

    const cs_lnum_t *face_ids
      = cs_gui_get_boundary_faces(boundaries->label[izone], &faces);

    zone_nbr = cs_gui_boundary_zone_number(izone+1);

    bft_printf("\n---zone %i label: %s\n", zone_nbr, boundaries->label[izone]);
    bft_printf("---zone %i nature: %s\n", zone_nbr, boundaries->nature[izone]);
    bft_printf("---zone %i number of faces: %i\n", zone_nbr, faces);
    bft_printf("----localization: %s\n",
               cs_gui_boundary_zone_localization(boundaries->label[izone]));

    if (cs_gui_strcmp(boundaries->nature[izone], "inlet")) {
      char * choice_v = _boundary_choice("inlet",
                                         boundaries->label[izone],
                                         "velocity_pressure",
                                         "choice");
      char * choice_d = _boundary_choice("inlet",
                                         boundaries->label[izone],
                                         "velocity_pressure",
                                         "direction");

      if (cs_gui_strcmp(choice_v, "norm"))
        bft_printf("-----velocity: %s => %12.5e \n", choice_v, boundaries->norm[izone]);
      if (cs_gui_strcmp(choice_v, "flow1") || cs_gui_strcmp(choice_v, "flow2"))
        bft_printf("-----velocity: %s => %12.5e \n", choice_v, boundaries->qimp[izone]);
      if (   cs_gui_strcmp(choice_v, "norm_formula")
          || cs_gui_strcmp(choice_v, "flow1_formula")
          || cs_gui_strcmp(choice_v, "flow2_formula"))
        bft_printf("-----velocity: %s => %s \n",
            choice_v, boundaries->velocity[izone]->string);
      if (cs_gui_strcmp(choice_d, "coordinates") || cs_gui_strcmp(choice_d, "translation"))
        bft_printf("-----direction: %s => %12.5e %12.5e %12.5e \n",
            choice_v, boundaries->dirx[izone],
            boundaries->diry[izone], boundaries->dirz[izone]);
      else if (cs_gui_strcmp(choice_d, "formula"))
        bft_printf("-----direction: %s => %s \n", choice_d,
            boundaries->direction[izone]->string);

      BFT_FREE(choice_v);
      BFT_FREE(choice_d);

      if (cs_gui_strcmp(vars->model, "solid_fuels")) {
        bft_printf("-----iqimp=%i, qimpat=%12.5e \n",
            iqimp[zone_nbr-1], qimpat[zone_nbr-1]);
        bft_printf("-----ientat=%i, ientcp=%i, timpat=%12.5e \n",
            ientat[zone_nbr-1], ientcp[zone_nbr-1], timpat[zone_nbr-1]);

        for (int icharb = 0; icharb < *ncharb; icharb++) {
          bft_printf("-----coal=%i, qimpcp=%12.5e, timpcp=%12.5e \n",
              icharb+1, qimpcp[icharb *(*nozppm)+zone_nbr-1],
              timpcp[icharb *(*nozppm)+zone_nbr-1]);

          for (int iclass = 0; iclass < nclpch[icharb]; k++)
            bft_printf("-----coal=%i, class=%i, distch=%f \n",
                icharb+1, iclass=1,
                distch[iclass * (*nozppm) * (*ncharm) +icharb * (*nozppm) +zone_nbr-1]);
        }
      }
      else if (cs_gui_strcmp(vars->model, "gas_combustion")) {
        bft_printf("-----iqimp=%i \n",
            iqimp[zone_nbr-1]);
        bft_printf("-----ientox=%i, ientfu=%i, ientgf=%i, ientgb=%i \n",
            ientox[zone_nbr-1], ientfu[zone_nbr-1],
            ientgf[zone_nbr-1], ientgb[zone_nbr-1]);
      }
      else if (cs_gui_strcmp(vars->model, "compressible_model")) {
        if (boundaries->itype[izone] == CS_ESICF) {
          bft_printf("-----imposed_inlet\n");
          bft_printf("-----premin=%i \n",boundaries->prein[zone_nbr-1]);
          bft_printf("-----rhoin=%i \n",boundaries->rhoin[zone_nbr-1]);
          bft_printf("-----tempin=%i \n",boundaries->tempin[zone_nbr-1]);
          bft_printf("-----entin=%i \n",boundaries->entin[zone_nbr-1]);
        }
        if (boundaries->itype[izone] == CS_EPHCF) {
          bft_printf("-----subsonic_inlet_PH\n");
          bft_printf("-----prein=%i \n",boundaries->prein[zone_nbr-1]);
          bft_printf("-----entin=%i \n",boundaries->entin[zone_nbr-1]);
        }
      }
      else {
        bft_printf("-----iqimp=%i, qimp=%12.5e \n",
            iqimp[zone_nbr-1], qimp[zone_nbr-1]);
      }
      bft_printf("-----icalke=%i, dh=%12.5e, xintur=%12.5e \n",
          icalke[zone_nbr-1], dh[zone_nbr-1], xintur[zone_nbr-1]);

      if (cs_gui_strcmp(vars->model, "atmospheric_flows"))
        bft_printf("-----iprofm=%i, automatic=%i \n",
            iprofm[zone_nbr-1], boundaries->meteo[izone].automatic);
    }

    if (faces > 0) {
      ifbr = face_ids[0];

      for (int f_id = 0; f_id < n_fields; f_id++) {
        const cs_field_t  *f = cs_field_by_id(f_id);
        const int var_key_id = cs_field_key_id("variable_id");
        ivar = cs_field_get_key_int(f, var_key_id) -1;
        if (f->type & CS_FIELD_VARIABLE) {
          bft_printf("------%s: icodcl=%i, "
                     "rcodcl(1)=%12.5e, rcodcl(2)=%12.5e, rcodcl(3)=%12.5e\n",
                     f->name,
                     icodcl[ivar *n_b_faces +ifbr ],
                     rcodcl[0 * n_b_faces * (*nvar) +ivar * n_b_faces +ifbr],
                     rcodcl[1 * n_b_faces * (*nvar) +ivar * n_b_faces +ifbr],
                     rcodcl[2 * n_b_faces * (*nvar) +ivar * n_b_faces +ifbr]);
        }
      }
    }
  }
#endif
}

/*----------------------------------------------------------------------------
 * Boundary conditions input verification
 *
 * Fortran Interface:
 *
 * SUBROUTINE UICLVE
 * *****************
 *
 * integer          nozppm  <-- max number of boundary conditions zone
 * integer          iale    <-- ale module activated
 * integer          itypfb  <-- type of boundary for each face
 * integer          izfppp  <-- zone number
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiclve, UICLVE)(const int  *nozppm,
                               const int  *iale,
                               int        *itypfb,
                               int        *izfppp)
{
  cs_lnum_t ifbr, ifac;
  int izone, zones, zone_nbr;
  int inature = 0;
  int inature2 = 0;
  cs_lnum_t faces = 0;

  zones   = cs_gui_boundary_zones_number();

  for (izone=0 ; izone < zones ; izone++) {
    if (cs_gui_strcmp(boundaries->nature[izone], "inlet")) {
      inature = CS_INLET;
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "wall")) {
      inature = CS_ROUGHWALL;
      if (boundaries->rough[izone] < 0.0)
        inature = CS_SMOOTHWALL;
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "outlet")
        || cs_gui_strcmp(boundaries->nature[izone], "imposed_p_outlet")) {
      inature = CS_OUTLET;
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "symmetry")) {
      inature = CS_SYMMETRY;
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "free_inlet_outlet")) {
      inature = CS_FREE_INLET;
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "free_surface") && *iale) {
      inature = CS_FREE_SURFACE;
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "undefined")) {
      inature = CS_INDEF;
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "groundwater")) {
      inature = CS_INDEF;
    }
    else
      bft_error(__FILE__, __LINE__, 0,
          _("boundary nature %s is unknown \n"),
          boundaries->nature[izone]);

    zone_nbr = cs_gui_boundary_zone_number(izone + 1);

    if (nozppm && zone_nbr > *nozppm)
      bft_error(__FILE__, __LINE__, 0,
                _("zone's label number %i is greater than %i,"
                  " the maximum allowed \n"), zone_nbr, *nozppm);

    /* Boundary types compatibilty:
       the nature of the boundary can be changed from smooth wall to
       rough wall or vice-versa between the GUI and the user code.
       It can also be changed from inlet to convective inlet and
       vice-versa */

    int enature = inature;
    if (inature == CS_ROUGHWALL)
      enature = CS_SMOOTHWALL;
    else if (inature == CS_CONVECTIVE_INLET)
      enature = CS_INLET;

    int atmo_auto = 0;
    int compr_auto = 0;

    if (cs_gui_strcmp(cs_glob_var->model, "atmospheric_flows")) {
      if (boundaries->meteo[izone].automatic) {
        if (inature == CS_INLET || inature == CS_OUTLET)
          atmo_auto = inature;
      }
    }
    else if (cs_gui_strcmp(cs_glob_var->model, "compressible_model")) {
      if (inature == CS_INLET || inature == CS_OUTLET)
        compr_auto = inature;
    }

    const cs_lnum_t *face_ids
      = cs_gui_get_boundary_faces(boundaries->label[izone], &faces);

    for (ifac = 0; ifac < faces; ifac++) {
      ifbr = face_ids[ifac];

      if (izfppp[ifbr] != zone_nbr)
        bft_error
          (__FILE__, __LINE__, 0,
           _("@                                                            \n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@                                                            \n"
             "@ @@ WARNING: BOUNDARY CONDITIONS ERROR                      \n"
             "@    *******                                                 \n"
             "@                                                            \n"
             "@    The zone %s does not have the same id number            \n"
             "@    in the GUI and in the user subroutine.                  \n"
             "@                                                            \n"
             "@    GUI zone number:             %i                         \n"
             "@    USER SUBROUTINE zone number: %i                         \n"
             "@                                                            \n"
             "@    The id number given in the GUI cannot be modified       \n"
             "@    in the user subroutine (fortran array IZFPPP).          \n"
             "@                                                            \n"
             "@    The calculation will stop.                              \n"
             "@                                                            \n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@                                                            \n"),
           boundaries->label[izone], zone_nbr, izfppp[ifbr]);

      inature2 = itypfb[ifbr];

      int enature2 = inature;
      if (inature2 == CS_ROUGHWALL)
        enature2 = CS_SMOOTHWALL;
      else if (inature2 == CS_CONVECTIVE_INLET)
        inature2 = CS_INLET;

      if (atmo_auto  && inature2 == 0)
        inature2 = inature;

      else if (compr_auto) {
        if (   (compr_auto == CS_INLET  && (   inature2 == CS_ESICF
                                            || inature2 == CS_EPHCF))
            || (compr_auto = CS_OUTLET  &&  (   inature2 == CS_SSPCF
                                             || inature2 == CS_SOPCF)))
          inature2 = inature;
      }

      if (enature2 != enature)
        bft_error
          (__FILE__, __LINE__, 0,
           _("@                                                            \n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@                                                            \n"
             "@ @@ WARNING: BOUNDARY CONDITIONS ERROR                      \n"
             "@    *******                                                 \n"
             "@                                                            \n"
             "@    The zone %s does not have the same nature               \n"
             "@    in the GUI and in the user subroutine.                  \n"
             "@                                                            \n"
             "@    GUI zone nature:             %s                         \n"
             "@    USER SUBROUTINE ITYPFB:      %i                         \n"
             "@                                                            \n"
             "@    The nature given in the GUI cannot be modified          \n"
             "@    in the user subroutine (fortran array ITYPFB).          \n"
             "@                                                            \n"
             "@    The calculation will stop.                              \n"
             "@                                                            \n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@                                                            \n"),
           boundaries->label[izone], boundaries->nature[izone], inature2);
    }
  } /*  for izone */
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Return number of boundary regions definition
 *----------------------------------------------------------------------------*/

int
cs_gui_boundary_zones_number(void)
{
  int zones = 0;
  char *path = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "boundary_conditions");
  cs_xpath_add_element(&path, "boundary");

  zones = cs_gui_get_nb_element(path);

  BFT_FREE(path);

  return zones;
}

/*-----------------------------------------------------------------------------
 * Return the nature of boundary condition for the given zone
 *----------------------------------------------------------------------------*/

char *
cs_gui_boundary_zone_nature(const int ith_zone)
{
  char *path = NULL;
  char *nature = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "boundary_conditions");
  cs_xpath_add_element_num(&path, "boundary", ith_zone);
  cs_xpath_add_attribute(&path, "nature");

  nature = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return nature;
}

/*-----------------------------------------------------------------------------
 * Return the label of boundary condition for the given zone
 *----------------------------------------------------------------------------*/

char *cs_gui_boundary_zone_label(const int ith_zone)
{
  char *path = NULL;
  char *label = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "boundary_conditions");
  cs_xpath_add_element_num(&path, "boundary", ith_zone);
  cs_xpath_add_attribute(&path, "label");

  label = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return label;
}


/*-----------------------------------------------------------------------------
 * Return the zone number of boundary condition for the given zone
 *----------------------------------------------------------------------------*/

int
cs_gui_boundary_zone_number(const int ith_zone)
{
  char *path = NULL;
  char *czone = NULL;
  int zone;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "boundary_conditions");
  cs_xpath_add_element_num(&path, "boundary", ith_zone);
  cs_xpath_add_attribute(&path, "name");

  czone = cs_gui_get_attribute_value(path);
  zone = atoi(czone);

  BFT_FREE(path);
  BFT_FREE(czone);

  return zone;
}

/*-----------------------------------------------------------------------------
 * Return the description of a boundary zone
 *
 * parameters:
 *   label                   <--  label of boundary zone
 *----------------------------------------------------------------------------*/

char *
cs_gui_boundary_zone_localization(const char  *label)
{
  char *path = NULL;
  char *localization = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "boundary_conditions",  "boundary");
  cs_xpath_add_test_attribute(&path, "label", label);
  cs_xpath_add_function_text(&path);

  localization = cs_gui_get_text_value(path);

  BFT_FREE(path);

  return localization;
}

/*-----------------------------------------------------------------------------
 * Helper to get the face list for the izone
 *
 * parameters:
 *   label     <--  boundary label
 *   n_faces   -->  number of faces
 *
 * returns:
 *   pointer to face list
 *----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_gui_get_boundary_faces(const char   *label,
                          cs_lnum_t    *n_faces)
{
  const cs_lnum_t *face_ids = NULL;

  const cs_boundary_zone_t *z = cs_boundary_zone_by_name(label);

  *n_faces = z->n_faces;
  face_ids = z->face_ids;

  return face_ids;
}

/*----------------------------------------------------------------------------
 * Free boundary conditions structures
 *
 * parameters:
 *   ncharb  <-- number of coals
 *----------------------------------------------------------------------------*/

void
cs_gui_boundary_conditions_free_memory(const int  *ncharb)
{
  int izone;
  int zones;
  int icharb;

  cs_var_t  *vars = cs_glob_var;

  /* clean memory for global private structure boundaries */

  if (boundaries != NULL) {

    zones = cs_gui_boundary_zones_number();
    for (izone=0 ; izone < zones ; izone++) {
      BFT_FREE(boundaries->label[izone]);
      BFT_FREE(boundaries->nature[izone]);
      mei_tree_destroy(boundaries->velocity[izone]);
      mei_tree_destroy(boundaries->direction[izone]);
      mei_tree_destroy(boundaries->headLoss[izone]);
      for (int f_id = 0; f_id < cs_field_n_fields(); f_id++) {
        const cs_field_t  *f = cs_field_by_id(f_id);

        if (f->type & CS_FIELD_VARIABLE)
          for (int i = 0; i < f->dim; i++)
            mei_tree_destroy(boundaries->scalar[f->id][izone * f->dim + i]);
      }
    }

    for (int f_id = 0; f_id < cs_field_n_fields(); f_id++) {
      const cs_field_t  *f = cs_field_by_id(f_id);
      if (f->type & CS_FIELD_VARIABLE) {
        BFT_FREE(boundaries->type_code[f->id]);
        BFT_FREE(boundaries->values[f->id]);
        BFT_FREE(boundaries->scalar[f->id]);
      }
    }

    if (cs_gui_strcmp(vars->model, "solid_fuels")) {
      for (izone = 0; izone < zones; izone++) {
        BFT_FREE(boundaries->qimpcp[izone]);
        BFT_FREE(boundaries->timpcp[izone]);
        for (icharb = 0; icharb < *ncharb; icharb++)
          BFT_FREE(boundaries->distch[izone][icharb]);
        BFT_FREE(boundaries->distch[izone]);
      }
      BFT_FREE(boundaries->ientat);
      BFT_FREE(boundaries->ientcp);
      BFT_FREE(boundaries->inmoxy);
      BFT_FREE(boundaries->timpat);
      BFT_FREE(boundaries->qimpcp);
      BFT_FREE(boundaries->timpcp);
      BFT_FREE(boundaries->distch);
    }
    if (cs_gui_strcmp(vars->model, "gas_combustion")) {
      BFT_FREE(boundaries->ientfu);
      BFT_FREE(boundaries->ientox);
      BFT_FREE(boundaries->ientgb);
      BFT_FREE(boundaries->ientgf);
      BFT_FREE(boundaries->tkent);
      BFT_FREE(boundaries->fment);
    }
    if (cs_gui_strcmp(vars->model, "compressible_model")) {
      BFT_FREE(boundaries->itype);
      BFT_FREE(boundaries->prein);
      BFT_FREE(boundaries->rhoin);
      BFT_FREE(boundaries->tempin);
      BFT_FREE(boundaries->entin);
    }
    if (cs_gui_strcmp(vars->model, "groundwater_model")) {
      for (izone=0 ; izone < zones ; izone++)
        if (boundaries->groundwat[izone] != NULL)
          mei_tree_destroy(boundaries->groundwat[izone]);
      BFT_FREE(boundaries->groundwat);
    }
    if (cs_gui_strcmp(vars->model, "atmospheric_flows"))
      BFT_FREE(boundaries->meteo);

    for (izone=0 ; izone < zones ; izone++) {
      if (boundaries->locator[izone] != NULL)
        boundaries->locator[izone]
          = ple_locator_destroy(boundaries->locator[izone]);
    }

    BFT_FREE(boundaries->label);
    BFT_FREE(boundaries->nature);
    BFT_FREE(boundaries->iqimp);
    BFT_FREE(boundaries->icalke);
    BFT_FREE(boundaries->qimp);
    BFT_FREE(boundaries->dh);
    BFT_FREE(boundaries->xintur);
    BFT_FREE(boundaries->type_code);
    BFT_FREE(boundaries->values);
    BFT_FREE(boundaries->rough);
    BFT_FREE(boundaries->norm);
    BFT_FREE(boundaries->dirx);
    BFT_FREE(boundaries->diry);
    BFT_FREE(boundaries->dirz);
    BFT_FREE(boundaries->velocity);
    BFT_FREE(boundaries->direction);
    BFT_FREE(boundaries->headLoss);
    BFT_FREE(boundaries->scalar);
    BFT_FREE(boundaries->preout);
    BFT_FREE(boundaries->locator);
    BFT_FREE(boundaries);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
