/*============================================================================
 * Management of the GUI parameters file: boundary conditions
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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_selector.h"

#include "mei_evaluate.h"

#include "cs_base.h"
#include "cs_gui_util.h"
#include "cs_gui_specific_physics.h"
#include "cs_gui_variables.h"
#include "cs_mesh.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_prototypes.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gui_boundary_conditions.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

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

/*----------------------------------------------------------------------------
 * Return the label of a scalar field, given by its scalar Id
 *
 * parameters:
 *   id <-- scalar field id
 *----------------------------------------------------------------------------*/

static inline const char *
_scalar_label(const int id)
{
  cs_field_t  *f = cs_field_by_id(cs_glob_var->scal_f_id[id]);
  return cs_field_get_label(f);
}

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
 * Value of velocity for sliding wall.
 *
 * parameters:
 *   label       <--  label of wall boundary condition
 *   izone       <--  number of zone
 *   ivar        <--  number of variable
 *----------------------------------------------------------------------------*/

static void
_sliding_wall(const char  *label,
              const int    izone,
              const int    ivar)
{
  char *path = NULL;
  double result = 0.0;

  cs_var_t  *vars = cs_glob_var;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "boundary_conditions");
  cs_xpath_add_element(&path, "wall");
  cs_xpath_add_test_attribute(&path, "label", label);
  cs_xpath_add_element(&path, "velocity_pressure");
  cs_xpath_add_test_attribute(&path, "choice", "on");
  cs_xpath_add_element(&path, "dirichlet");
  cs_xpath_add_test_attribute(&path, "name", vars->name[ivar]);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &result)) {
    boundaries->type_code[vars->rtp[ivar]][izone] = DIRICHLET;
    boundaries->values[vars->rtp[ivar]][izone].val1 = result;
  }
  BFT_FREE(path);
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
  cs_xpath_add_elements(&path, 2, "velocity_pressure", tag);
  cs_xpath_add_attribute(&path, "status");

  if (cs_gui_get_status(path, &result))
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
_boundary_scalar_init_mei_tree(const char  *formula,
                               const char  *symbols[],
                               int          symbol_size)
{
  int i = 0;

  /* return an empty interpreter */
  mei_tree_t *tree = mei_tree_new(formula);

  /* add commun variables */
  mei_tree_insert(tree, "x", 0.0);
  mei_tree_insert(tree, "y", 0.0);
  mei_tree_insert(tree, "z", 0.0);

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
 *   nsca        <--  number of user scalar
 *----------------------------------------------------------------------------*/

static void
_boundary_scalar(const char  *nature,
                 int          izone,
                 int          nsca)
{
  int numvar;
  char *path = NULL;
  char *path_commun = NULL;
  char *path2 = NULL;
  char *choice = NULL;
  double result;

  cs_var_t  *vars = cs_glob_var;

  numvar  = vars->nvar - vars->nscaus - vars->nscapp;
  numvar = numvar + nsca;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "boundary_conditions", nature);
  cs_xpath_add_test_attribute(&path, "label", boundaries->label[izone]);
  cs_xpath_add_element(&path, "scalar");
  cs_xpath_add_test_attribute(&path, "label", _scalar_label(nsca));

  BFT_MALLOC(path_commun, strlen(path)+1, char);
  strcpy(path_commun, path);

  BFT_MALLOC(path2, strlen(path)+1, char);
  strcpy(path2, path);

  cs_xpath_add_attribute(&path_commun, "choice");
  choice = cs_gui_get_attribute_value(path_commun);

  if (choice != NULL) {

    if (cs_gui_strcmp(choice, "dirichlet") || cs_gui_strcmp(choice, "exchange_coefficient") || cs_gui_strcmp(choice, "wall_function")) {
      cs_xpath_add_element(&path, "dirichlet");
      cs_xpath_add_function_text(&path);
      if (cs_gui_get_double(path, &result)) {
        if (cs_gui_strcmp(choice, "wall_function")) {
          boundaries->type_code[vars->rtp[numvar]][izone] = WALL_FUNCTION;
        } else {
          boundaries->type_code[vars->rtp[numvar]][izone] = DIRICHLET;
        }
        boundaries->values[vars->rtp[numvar]][izone].val1 = result;
      }

    } else if(cs_gui_strcmp(choice, "neumann")) {
      cs_xpath_add_element(&path, "neumann");
      cs_xpath_add_function_text(&path);
      if (cs_gui_get_double(path, &result)) {
        boundaries->type_code[vars->rtp[numvar]][izone] = NEUMANN;
        boundaries->values[vars->rtp[numvar]][izone].val3 = result;
      }
    } else if (cs_gui_strcmp(choice, "dirichlet_formula")) {
      cs_xpath_add_element(&path, "dirichlet_formula");
      cs_xpath_add_function_text(&path);
      if (cs_gui_get_text_value(path)){
        const char *sym[] = {_scalar_label(nsca)};
        boundaries->type_code[vars->rtp[numvar]][izone] = DIRICHLET_FORMULA;
        boundaries->scalar[vars->rtp[numvar]][izone] = _boundary_scalar_init_mei_tree(cs_gui_get_text_value(path), sym, 1);
      }
    } else if (cs_gui_strcmp(choice, "neumann_formula")) {
      cs_xpath_add_element(&path, "neumann_formula");
      cs_xpath_add_function_text(&path);
      if (cs_gui_get_text_value(path)){
        const char *sym[] = {"flux"};
        boundaries->type_code[vars->rtp[numvar]][izone] = NEUMANN_FORMULA;
        boundaries->scalar[vars->rtp[numvar]][izone] = _boundary_scalar_init_mei_tree(cs_gui_get_text_value(path), sym, 1);
      }
    } else if (cs_gui_strcmp(choice, "exchange_coefficient_formula")){
      cs_xpath_add_element (&path, "exchange_coefficient_formula");
      cs_xpath_add_function_text(&path);
      if (cs_gui_get_text_value(path)){
        const char *sym[] = {_scalar_label(nsca),"hc"};
        boundaries->type_code[vars->rtp[numvar]][izone] = EXCHANGE_COEFF_FORMULA;
        boundaries->scalar[vars->rtp[numvar]][izone] = _boundary_scalar_init_mei_tree(cs_gui_get_text_value(path), sym, 2);
      }
    }

    if (cs_gui_strcmp(choice, "exchange_coefficient")) {
      cs_xpath_add_element(&path2, "exchange_coefficient");
      cs_xpath_add_function_text(&path2);
      if (cs_gui_get_double(path2, &result)) {
        boundaries->type_code[vars->rtp[numvar]][izone] = EXCHANGE_COEFF;
        boundaries->values[vars->rtp[numvar]][izone].val2 = result;
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
  int    size = 0;
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
  list_of_coals = cs_gui_get_attribute_values(path1, &size);
  BFT_FREE(path1);

  /* if there is no coal, it is an inlet only for oxydant */
  if (size == 0) {
    boundaries->ientat[izone] = 1;
    boundaries->ientcp[izone] = 0;
    BFT_FREE(path0);
  }
  else {
    if (size != *ncharb)
      bft_error(__FILE__, __LINE__, 0,
                _("Invalid number of coal-> dp_FCP: %i xml: %i\n"),
                *ncharb, size);

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
      list_of_classes = cs_gui_get_attribute_values(path1, &size);
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

      for (iclass = 0; iclass < nclpch[icoal]; iclass++)
        BFT_FREE(list_of_classes[iclass]);
      BFT_FREE(list_of_classes);

      BFT_FREE(path2);
      BFT_FREE(path3);
      BFT_FREE(path4);
    }

    for (icoal = 0; icoal < *ncharb; icoal++)
      BFT_FREE(list_of_coals[icoal]);
    BFT_FREE(list_of_coals);

    BFT_FREE(path0);
  }
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
 * Get electrical's data for boundary.
 *
 * parameters:
 *   izone       -->  number of the current zone
 *----------------------------------------------------------------------------*/

static void
_boundary_elec(const char  *nature,
               const int    izone,
               int          nsca)
{
  int numvar;
  char *path = NULL;
  char *path_commun = NULL;
  char *path2 = NULL;
  char *choice = NULL;
  double result;

  cs_var_t  *vars = cs_glob_var;

  numvar  = vars->nvar - vars->nscaus - vars->nscapp;
  numvar = numvar + nsca;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "boundary_conditions", nature);
  cs_xpath_add_test_attribute(&path, "label", boundaries->label[izone]);
  cs_xpath_add_element(&path, "scalar");
  cs_xpath_add_test_attribute(&path, "label", _scalar_label(nsca));

  BFT_MALLOC(path_commun, strlen(path)+1, char);
  strcpy(path_commun, path);

  BFT_MALLOC(path2, strlen(path)+1, char);
  strcpy(path2, path);

  cs_xpath_add_attribute(&path_commun, "choice");
  choice = cs_gui_get_attribute_value(path_commun);

  if (choice != NULL) {

    if (   cs_gui_strcmp(choice, "dirichlet")
        || cs_gui_strcmp(choice, "exchange_coefficient")
        || cs_gui_strcmp(choice, "wall_function")) {
      cs_xpath_add_element(&path, "dirichlet");
      cs_xpath_add_function_text(&path);
      if (cs_gui_get_double(path, &result)) {
        if (cs_gui_strcmp(choice, "wall_function")) {
          boundaries->type_code[vars->rtp[numvar]][izone] = WALL_FUNCTION;
        } else {
          boundaries->type_code[vars->rtp[numvar]][izone] = DIRICHLET;
        }
        boundaries->values[vars->rtp[numvar]][izone].val1 = result;
      }

    } else if(cs_gui_strcmp(choice, "neumann")) {
      cs_xpath_add_element(&path, "neumann");
      cs_xpath_add_function_text(&path);
      if (cs_gui_get_double(path, &result)) {
        boundaries->type_code[vars->rtp[numvar]][izone] = NEUMANN;
        boundaries->values[vars->rtp[numvar]][izone].val3 = result;
      }
    } else if (cs_gui_strcmp(choice, "dirichlet_formula")) {
      cs_xpath_add_element(&path, "dirichlet_formula");
      cs_xpath_add_function_text(&path);
      if (cs_gui_get_text_value(path)){
        const char *sym[] = {_scalar_label(nsca)};
        boundaries->type_code[vars->rtp[numvar]][izone] = DIRICHLET_FORMULA;
        boundaries->scalar[vars->rtp[numvar]][izone] = _boundary_scalar_init_mei_tree(cs_gui_get_text_value(path), sym, 1);
      }
    } else if (cs_gui_strcmp(choice, "neumann_formula")) {
      cs_xpath_add_element(&path, "neumann_formula");
      cs_xpath_add_function_text(&path);
      if (cs_gui_get_text_value(path)){
        const char *sym[] = {"flux"};
        boundaries->type_code[vars->rtp[numvar]][izone] = NEUMANN_FORMULA;
        boundaries->scalar[vars->rtp[numvar]][izone] = _boundary_scalar_init_mei_tree(cs_gui_get_text_value(path), sym, 1);
      }
    } else if (cs_gui_strcmp(choice, "exchange_coefficient_formula")){
      cs_xpath_add_element (&path, "exchange_coefficient_formula");
      cs_xpath_add_function_text(&path);
      if (cs_gui_get_text_value(path)){
        const char *sym[] = {_scalar_label(nsca),"hc"};
        boundaries->type_code[vars->rtp[numvar]][izone] = EXCHANGE_COEFF_FORMULA;
        boundaries->scalar[vars->rtp[numvar]][izone] = _boundary_scalar_init_mei_tree(cs_gui_get_text_value(path), sym, 2);
      }
    } else if (cs_gui_strcmp(choice, "dirichlet_implicit")){
        boundaries->type_code[vars->rtp[numvar]][izone] = DIRICHLET_IMPLICIT;
    } else if (cs_gui_strcmp(choice, "neumann_implicit")){
        boundaries->type_code[vars->rtp[numvar]][izone] = NEUMANN_IMPLICIT;
    }

    if (cs_gui_strcmp(choice, "exchange_coefficient")) {
      cs_xpath_add_element(&path2, "exchange_coefficient");
      cs_xpath_add_function_text(&path2);
      if (cs_gui_get_double(path2, &result)) {
        boundaries->type_code[vars->rtp[numvar]][izone] = EXCHANGE_COEFF;
        boundaries->values[vars->rtp[numvar]][izone].val2 = result;
      }
    }

    BFT_FREE(choice);
  }

  BFT_FREE(path);
  BFT_FREE(path2);
  BFT_FREE(path_commun);
}


/*-----------------------------------------------------------------------------
 * Get compressible data for inlet.
 *
 * parameters:
 *   izone       -->  number of the current zone
 *----------------------------------------------------------------------------*/

static void
_inlet_compressible(int         izone,
                    const  int *iesicf,
                    const  int *iephcf)
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
    boundaries->itype[izone] = *iesicf;
    cs_xpath_add_element(&path1, "pressure");
    BFT_MALLOC(path_value, strlen(path1) + 1, char);
    strcpy(path_value, path1);
    cs_xpath_add_attribute(&path1, "status");
    status = cs_gui_get_attribute_value(path1);
    if (cs_gui_strcmp(status, "on")){
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
    if (cs_gui_strcmp(status, "on")){
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
    if (cs_gui_strcmp(status, "on")){
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
    if (cs_gui_strcmp(status, "on")){
      cs_xpath_add_function_text(&path_value);
      if (cs_gui_get_double(path_value, &value))
        boundaries->entin[izone] = value;
    }
    BFT_FREE(status);
    BFT_FREE(path_value);
  }
  else if (cs_gui_strcmp(buff, "subsonic_inlet_PH")){
    boundaries->itype[izone] = *iephcf;

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
_outlet_compressible(int         izone,
                     const int  *isspcf,
                     const int  *isopcf)
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
    boundaries->itype[izone] = *isspcf;
  }
  else if (cs_gui_strcmp(buff, "subsonic_outlet")) {
    boundaries->itype[izone] = *isopcf;
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

/*----------------------------------------------------------------------------
 * Boundary conditions treatment: global structure initialization
 *
 * parameters:
 *   nfabor               <-- number of boundary faces
 *   nozppm               <-- max number of boundary conditions zone
 *   iscalt               <-- index of the thermal scalar
 *   ncharb               <-- number of simulated coals
 *   nclpch               <-- number of simulated class per coals
 *   ncharb               <-- number of simulated gaz for electrical models
 *   izfppp               <-- zone number for each boundary face
 *   iesicf               <-- type of boundary: imposed inlet (compressible)
 *   isspcf               <-- type of boundary: supersonic outlet (compressible)
 *   iephcf               <-- type of boundary: subsonic inlet imposed total pressure and total enthalpy (compressible)
 *   isopcf               <-- type of boundary: subsonic outlet (compressible)
 *----------------------------------------------------------------------------*/

static void
_init_boundaries(const int  *nfabor,
                 const int  *nozppm,
                 const int  *iscalt,
                 const int  *ncharb,
                 const int  *nclpch,
                       int  *izfppp,
                 const int  *iesicf,
                 const int  *isspcf,
                 const int  *iephcf,
                 const int  *isopcf)
{
    int faces = 0;
    int zones = 0;
    int ifac, izone, ith_zone, zone_nbr;
    int ifbr, i;
    int ivar, isca, icharb, iclass;
    char *choice = NULL;
    char *choice_v = NULL;
    char *choice_d = NULL;
    char *nature = NULL;
    char *label = NULL;
    int  *faces_list = NULL;

    cs_var_t  *vars = cs_glob_var;
    assert(boundaries == NULL);

    zones = cs_gui_boundary_zones_number();

    BFT_MALLOC(boundaries,            1,          cs_boundary_t);
    BFT_MALLOC(boundaries->label,     zones,      char*        );
    BFT_MALLOC(boundaries->nature,    zones,      char*        );
    BFT_MALLOC(boundaries->type_code, vars->nvar, int*         );
    BFT_MALLOC(boundaries->values,    vars->nvar, cs_val_t*    );
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

    BFT_MALLOC(boundaries->velocity,  zones,      mei_tree_t*  );
    BFT_MALLOC(boundaries->direction, zones,      mei_tree_t*  );
    BFT_MALLOC(boundaries->scalar,    vars->nvar, mei_tree_t** );

    if (cs_gui_strcmp(vars->model, "solid_fuels"))
    {
        BFT_MALLOC(boundaries->ientat, zones, int      );
        BFT_MALLOC(boundaries->inmoxy, zones, double   );
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
                           nclpch[icharb],
                           double);
        }
    }
    else if (cs_gui_strcmp(vars->model, "gas_combustion"))
    {
        BFT_MALLOC(boundaries->ientfu,  zones, int   );
        BFT_MALLOC(boundaries->ientox,  zones, int   );
        BFT_MALLOC(boundaries->ientgb,  zones, int   );
        BFT_MALLOC(boundaries->ientgf,  zones, int   );
        BFT_MALLOC(boundaries->tkent,   zones, double);
        BFT_MALLOC(boundaries->fment,   zones, double);
    }
    else if (cs_gui_strcmp(vars->model, "compressible_model"))
    {
        BFT_MALLOC(boundaries->itype,   zones, int   );
        BFT_MALLOC(boundaries->prein,   zones, double);
        BFT_MALLOC(boundaries->rhoin,   zones, double);
        BFT_MALLOC(boundaries->tempin,  zones, double);
        BFT_MALLOC(boundaries->entin,   zones, double);
        BFT_MALLOC(boundaries->preout,  zones, double);
    }
    else
    {
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
    {
       BFT_MALLOC(boundaries->meteo, zones, cs_meteo_t);
    }
    else
    {
        boundaries->meteo = NULL;
    }

    for (ivar = 0; ivar < vars->nvar; ivar++)
    {
        i = vars->rtp[ivar];
        BFT_MALLOC(boundaries->type_code[i], zones, int);
        BFT_MALLOC(boundaries->values[i], zones, cs_val_t);
        BFT_MALLOC(boundaries->scalar[i], zones, mei_tree_t *);
    }

    for (izone = 0; izone < zones; izone++)
    {
        boundaries->iqimp[izone]  = 0;
        boundaries->qimp[izone]   = 0;
        boundaries->norm[izone]   = 0;
        boundaries->icalke[izone] = 0;
        boundaries->dh[izone]     = 0;
        boundaries->xintur[izone] = 0;
        boundaries->rough[izone]  = -999;
        boundaries->velocity[izone] = NULL;
        boundaries->direction[izone] = NULL;

        if (cs_gui_strcmp(vars->model, "solid_fuels"))
        {
            boundaries->ientat[izone] = 0;
            boundaries->inmoxy[izone] = 1;
            boundaries->ientcp[izone] = 0;
            boundaries->timpat[izone] = 0;

            for (icharb = 0; icharb < *ncharb; icharb++)
            {
                boundaries->qimpcp[izone][icharb] = 0;
                boundaries->timpcp[izone][icharb] = 0;

                for (iclass = 0; iclass < nclpch[icharb]; iclass++)
                    boundaries->distch[izone][icharb][iclass] = 0;
            }
        }
        else if (cs_gui_strcmp(vars->model, "gas_combustion"))
        {
            boundaries->ientfu[izone]  = 0;
            boundaries->ientox[izone]  = 0;
            boundaries->ientgb[izone]  = 0;
            boundaries->ientgf[izone]  = 0;
            boundaries->tkent[izone]   = 0;
            boundaries->fment[izone]   = 0;
        }
        else if (cs_gui_strcmp(vars->model, "compressible_model"))
        {
            boundaries->itype[izone]     = 0;
            boundaries->prein[izone]     = 0;
            boundaries->rhoin[izone]     = 0;
            boundaries->tempin[izone]    = 0;
            boundaries->preout[izone]    = 0;
            boundaries->entin[izone]     = 0;
        }
        else if (cs_gui_strcmp(vars->model, "atmospheric_flows"))
        {
           boundaries->meteo[izone].read_data = 0;
           boundaries->meteo[izone].automatic = 0;
        }
    }

    /* Initialization of boundary->type_code and boundary->values */

    for (ivar = 0; ivar < vars->nvar; ivar++)
    {
        i = vars->rtp[ivar];
        for (izone = 0; izone < zones; izone++)
        {
            boundaries->type_code[i][izone] = -1;
            boundaries->values[i][izone].val1 = 1.e30;
            boundaries->values[i][izone].val2 = 1.e30;
            boundaries->values[i][izone].val3 = 0.;
            boundaries->scalar[i][izone] = NULL;
        }
    }

    for (ifac = 0; ifac < *nfabor; ifac++)
        izfppp[ifac] = 0;

    /* filling of the "boundaries" structure */

    for (izone = 0; izone < zones; izone++)
    {
        /* nature, label of the ith initialization zone */

        ith_zone = izone + 1;
        nature = cs_gui_boundary_zone_nature(ith_zone);
        label = cs_gui_boundary_zone_label(ith_zone);

        BFT_MALLOC(boundaries->label[izone], strlen(label)+1, char);
        strcpy(boundaries->label[izone], label);

        BFT_MALLOC(boundaries->nature[izone], strlen(nature)+1, char);
        strcpy(boundaries->nature[izone], nature);

        if (cs_gui_strcmp(nature, "inlet"))
        {
            /* Inlet: VELOCITY */
            choice_v = _boundary_choice("inlet", label, "velocity_pressure", "choice");
            choice_d = _boundary_choice("inlet", label, "velocity_pressure", "direction");

            if (cs_gui_strcmp(choice_v, "norm"))
            {
                _inlet_data(label, "norm", &boundaries->norm[izone]);
            }
            else if (cs_gui_strcmp(choice_v, "flow1"))
            {
                _inlet_data(label, "flow1", &boundaries->qimp[izone]);
                boundaries->iqimp[izone] = 1;
            }
            else if (cs_gui_strcmp(choice_v, "flow2"))
            {
                _inlet_data(label, "flow2", &boundaries->qimp[izone]);
                boundaries->iqimp[izone] = 2;
            }
            else if (cs_gui_strcmp(choice_v, "norm_formula"))
            {
                const char *sym[] = {"u_norm"};
                boundaries->velocity[izone] = _boundary_init_mei_tree(_inlet_formula(label, choice_v), sym, 1);
            }
            else if (cs_gui_strcmp(choice_v, "flow1_formula"))
            {
                const char *sym[] = {"q_m"};
                boundaries->velocity[izone] = _boundary_init_mei_tree(_inlet_formula(label, choice_v), sym, 1);
                boundaries->iqimp[izone] = 1;
            }
            else if (cs_gui_strcmp(choice_v, "flow2_formula"))
            {
                const char *sym[] = {"q_v"};
                boundaries->velocity[izone] = _boundary_init_mei_tree(_inlet_formula(label, choice_v), sym, 1);
                boundaries->iqimp[izone] = 2;
            }
            if (cs_gui_strcmp(choice_d, "coordinates"))
            {
                _inlet_data(label, "direction_x", &boundaries->dirx[izone]);
                _inlet_data(label, "direction_y", &boundaries->diry[izone]);
                _inlet_data(label, "direction_z", &boundaries->dirz[izone]);
            }
            else if (cs_gui_strcmp(choice_d, "formula"))
            {
                const char *sym[] = {"dir_x", "dir_y", "dir_z"};
                boundaries->direction[izone] =
                    _boundary_init_mei_tree(_inlet_formula(label, "direction_formula"), sym, 3);
            }
            BFT_FREE(choice_v);
            BFT_FREE(choice_d);

            /* Inlet: data for COAL */
            if (cs_gui_strcmp(vars->model, "solid_fuels"))
            {
                _inlet_data(label, "temperature", &boundaries->timpat[izone]);
                _inlet_data(label, "oxydant", &boundaries->inmoxy[izone]);
                _inlet_coal(izone, ncharb, nclpch);
            }

            if (cs_gui_strcmp(vars->model, "gas_combustion"))
                _inlet_gas(izone);

            /* Inlet: data for COMPRESSIBLE MODEL */
            if (cs_gui_strcmp(vars->model, "compressible_model"))
            {
                _inlet_compressible(izone, iesicf, iephcf);
            }

            /* Inlet: data for ATMOSPHERIC FLOWS */
            if (cs_gui_strcmp(vars->model, "atmospheric_flows"))
            {
                _boundary_status("inlet", label, "meteo_data", &boundaries->meteo[izone].read_data);
                _boundary_status("inlet", label, "meteo_automatic", &boundaries->meteo[izone].automatic);
            }

            /* Inlet: TURBULENCE */
            choice = _boundary_choice("inlet", label, "turbulence", "choice");
            _inlet_turbulence(choice, izone);
            BFT_FREE(choice);

            /* Inlet: Thermal scalar */
            if (*iscalt > 0)
                _boundary_scalar("inlet", izone, *iscalt-1);

            /* Inlet: USER SCALARS */
            for (isca = 0; isca < vars->nscaus; isca++)
                _boundary_scalar("inlet", izone, isca);

            /* Inlet: METEO SCALARS */
            if (cs_gui_strcmp(vars->model, "atmospheric_flows"))
                for (isca = 0; isca < vars->nscapp; isca++)
                    _boundary_scalar("inlet", izone, isca+vars->nscaus);

        }
        else if (cs_gui_strcmp(nature, "wall"))
        {
            /* Wall: VELOCITY */
            choice = _boundary_choice("wall", label, "velocity_pressure", "choice");

            if (cs_gui_strcmp(choice, "on"))
            {
                for (ivar = 1; ivar < 4; ivar++)
                    _sliding_wall(label, izone, ivar);
            }
            BFT_FREE(choice);

            /* Wall: ROUGH */
            _wall_roughness(label, izone);

            /* Wall: Thermal scalar */
            if (*iscalt > 0)
                _boundary_scalar("wall", izone, *iscalt-1);

            /* Wall: USER SCALARS */
            for (isca = 0; isca < vars->nscaus; isca++)
                _boundary_scalar("wall", izone, isca);

            /* Wall: METEO SCALARS */
            if (cs_gui_strcmp(vars->model, "atmospheric_flows"))
                for (isca = 0; isca < vars->nscapp; isca++)
                    _boundary_scalar("wall", izone, isca+vars->nscaus);

        }
        else if (cs_gui_strcmp(nature, "outlet"))
        {
            /* Outlet: data for ATMOSPHERIC FLOWS */
            if (cs_gui_strcmp(vars->model, "atmospheric_flows"))
            {
                _boundary_status("outlet", label, "meteo_data", &boundaries->meteo[izone].read_data);
                _boundary_status("outlet", label, "meteo_automatic", &boundaries->meteo[izone].automatic);
            }

            /* Outlet: Thermal scalar */
            if (*iscalt > 0)
                _boundary_scalar("outlet", izone, *iscalt-1);

            /* Outlet: USER SCALARS */
            for (isca = 0; isca < vars->nscaus; isca++)
                _boundary_scalar("outlet", izone, isca);

            /* Outlet: METEO SCALARS */
            if (cs_gui_strcmp(vars->model, "atmospheric_flows"))
                for (isca = 0; isca < vars->nscapp; isca++)
                    _boundary_scalar("outlet", izone, isca+vars->nscaus);

            /* Inlet: data for COMPRESSIBLE MODEL */
            if (cs_gui_strcmp(vars->model, "compressible_model"))
                _outlet_compressible(izone, isspcf, isopcf);
        }  /* if (cs_gui_strcmp(nature, "outlet")) */

        if (cs_gui_strcmp(vars->model, "joule_effect"))
          for (isca = 0; isca < vars->nscapp; isca++)
            _boundary_elec(nature, izone, isca+vars->nscaus);
        BFT_FREE(nature);
        BFT_FREE(label);

    }  /* for izones */

    for (izone = 0 ; izone < zones ; izone++)
    {
        zone_nbr = cs_gui_boundary_zone_number(izone + 1);
        faces_list = cs_gui_get_faces_list(izone,
                                           boundaries->label[izone],
                                           *nfabor,
                                           *nozppm,
                                           &faces);

        /* check if faces are already marked with a zone number */

        for (ifac = 0; ifac < faces; ifac++)
        {
            ifbr = faces_list[ifac] -1;

            if (izfppp[ifbr] != 0)
            {
                bft_error(__FILE__, __LINE__, 0,
                _("@                                                            \n"
                  "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                  "@                                                            \n"
                  "@ @@ WARNING: BOUNDARY CONDITIONS ERROR                      \n"
                  "@    *******                                                 \n"
                  "@                                                            \n"
                  "@    In the zone %s has a face already marked                \n"
                  "@    with a zone number.                                     \n"
                  "@                                                            \n"
                  "@    new zone number:             %i                         \n"
                  "@    previous zone number:        %i                         \n"
                  "@                                                            \n"
                  "@    It seems that zones definitions are overlapping.        \n"
                  "@                                                            \n"
                  "@    The calculation will stop.                              \n"
                  "@                                                            \n"
                  "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                  "@                                                            \n"),
                  boundaries->label[izone], zone_nbr, izfppp[ifbr]);

            }
            else
            {
                izfppp[ifbr] = zone_nbr;
            }
        } /* for ifac */
        BFT_FREE(faces_list);
    } /*  for izone */
}

/*============================================================================
 * C API public functions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Return number of boundary regions definition
 *----------------------------------------------------------------------------*/

int cs_gui_boundary_zones_number(void)
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

char *cs_gui_boundary_zone_nature(const int ith_zone)
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

int cs_gui_boundary_zone_number(const int ith_zone)
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

char *cs_gui_boundary_zone_localization(const char  *label)
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
 *   izone   <--  zone index
 *   label   <--  boundary label
 *   nfabor  <--  number of boundary faces
 *   nozppm  <--  max number of boundary zone for preefined physics
 *   faces   -->  number of faces
 *----------------------------------------------------------------------------*/

int*
cs_gui_get_faces_list(int          izone,
                      const char  *label,
                      int          nfabor,
                      int          nozppm,
                      int         *faces )
{
  int  c_id        = 0;
  int  *faces_list = NULL;
  char *description = NULL;

  int  boundary_zones = cs_gui_boundary_zone_number(izone + 1);

  if (nozppm && boundary_zones >  nozppm)
    bft_error(__FILE__, __LINE__, 0,
              _("zone's label number %i is greater than %i,"
                " the maximum allowed \n"), boundary_zones , nozppm);

  description = cs_gui_boundary_zone_localization(label);

  /* list of faces building */
  BFT_MALLOC(faces_list, nfabor, int);

  c_id = fvm_selector_get_list(cs_glob_mesh->select_b_faces,
                               description,
                               faces,
                               faces_list);

  if (fvm_selector_n_missing(cs_glob_mesh->select_b_faces, c_id) > 0) {
    const char *missing
      = fvm_selector_get_missing(cs_glob_mesh->select_b_faces, c_id, 0);
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("The group or attribute \"%s\" in the selection\n"
                 "criteria:\n"  "\"%s\"\n"
                 " does not correspond to any boundary face.\n"),
               missing, description);
  }
  BFT_FREE(description);
  return faces_list;
}

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
 * integer          ntcabs  <-- current iteration number
 * integer          nfabor  <-- number of boundary faces
 * integer          nozppm  <-- max number of boundary conditions zone
 * integer          ncharm  <-- maximal number of coals
 * integer          ncharb  <-- number of simulated coals
 * integer          nclpch  <-- number of simulated class per coals
 * integer          ngasg   <-- number of simulated gas for electrical models
 * integer          iindef  <-- type of boundary: not defined
 * integer          ientre  <-- type of boundary: inlet
 * integer          iesicf  <-- type of boundary: imposed inlet (compressible)
 * integer          isspcf  <-- type of boundary: supersonic outlet (compressible)
 * integer          iephcf  <-- type of boundary: subsonic inlet imposed total pressure and total enthalpy (compressible)
 * integer          isopcf  <-- type of boundary: subsonic outlet (compressible)
 * integer          iparoi  <-- type of boundary: smooth wall
 * integer          iparug  <-- type of boundary: rough wall
 * integer          isymet  <-- type of boundary: symetry
 * integer          isolib  <-- type of boundary: outlet
 * integer          isca    <-- indirection array for scalar number
 * integer          ipr     <-- rtp index for pressure
 * integer          iscalt  <-- index of the thermal scalar
 * integer          itempk  <-- rtp index for temperature (in K)
 * integer          ienerg  <-- rtp index for energy total
 * integer          iqimp   <-- 1 if flow rate is applied
 * integer          icalke  <-- 1 for automatic turbulent boundary conditions
 * integer          ientat  <-- 1 for air temperature boundary conditions (coal)
 * integer          ientcp  <-- 1 for coal temperature boundary conditions (coal)
 * INTEGER          INMOXY  <-- coal: number of oxydant for the current inlet
 * integer          ientox  <-- 1 for an air fow inlet (gas combustion)
 * integer          ientfu  <-- 1 for fuel flow inlet (gas combustion)
 * integer          ientgb  <-- 1 for burned gas inlet (gas combustion)
 * integer          ientgf  <-- 1 for unburned gas inlet (gas combustion)
 * integer          iprofm  <-- atmospheric flows: on/off for profile from data
 * double precision coejou  <-- electric arcs
 * double precision dpot    <-- electric arcs : potential difference
 * double precision rtpa    <-- rtpa for implicit flux
 * integer          itypfb  <-- type of boundary for each face
 * integer          izfppp  <-- zone number for each boundary face
 * integer          icodcl  <-- boundary conditions array type
 * double precision dtref   <-- time step
 * double precision ttcabs  <-- current time
 * double precision surfbo  <-- boundary faces surface
 * double precision cgdfbo  <-- boundary faces center of gravity
 * double precision qimp    <-- inlet flow rate
 * double precision qimpat  <-- inlet air flow rate (coal)
 * double precision qimpcp  <-- inlet coal flow rate (coal)
 * double precision dh      <-- hydraulic diameter
 * double precision xintur  <-- turbulent intensity
 * double precision timpat  <-- air temperature boundary conditions (coal)
 * double precision timpcp  <-- inlet coal temperature (coal)
 * double precision tkent   <-- inlet temperature (gas combustion)
 * double precision fment   <-- Mean Mixture Fraction at Inlet (gas combustion)
 * double precision distch  <-- ratio for each coal
 * double precision rcodcl  <-- boundary conditions array value
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiclim, UICLIM)(const int  *ntcabs,
                               const int  *nfabor,
                               const int  *nozppm,
                               const int  *ncharm,
                               const int  *ncharb,
                               const int  *nclpch,
                               const int  *iindef,
                               const int  *ientre,
                               const int  *iesicf,
                               const int  *isspcf,
                               const int  *iephcf,
                               const int  *isopcf,
                               const int  *iparoi,
                               const int  *iparug,
                               const int  *isymet,
                               const int  *isolib,
                               const int  *isca,
                               const int  *ipr,
                               const int  *iscalt,
                               const int  *itempk,
                               const int  *ienerg,
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
                               double     *coejou,
                               double     *dpot,
                               cs_real_t   rtpa[],
                               int        *ielcor,
                               int        *ipotr,
                               int        *ipoti,
                               int         ipotva[],
                               int        *ncelet,
                               int        *itypfb,
                               int        *izfppp,
                               int        *icodcl,
                               double     *dtref,
                               double     *ttcabs,
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
                               double     *rcodcl)
{
  int faces = 0;
  int zones = 0;
  int izone, ith_zone, zone_nbr;
  int ifac, ifbr;
  int i, ivar, icharb, iclass, iwall;
  double t0;
  double norm = 0.;
  double X[3];
  int *faces_list = NULL;
  char *choice_v = NULL;
  char *choice_d = NULL;
  char *path1 = NULL;
  char *formula = NULL;
  char *model = NULL;
  mei_tree_t *ev_formula  = NULL;

  cs_var_t  *vars = cs_glob_var;
  int *ifabor = cs_glob_mesh->b_face_cells;

  zones = cs_gui_boundary_zones_number();

  /* First iteration only: memory allocation */

  if (boundaries == NULL)
    _init_boundaries(nfabor, nozppm, iscalt, ncharb, nclpch,
                     izfppp, iesicf, isspcf, iephcf, isopcf);

  /* At each time-step, loop on boundary faces:
     One sets itypfb, rcodcl and icodcl thanks to the arrays
     of the structures "conditions.limites" defined
     in the first part ofthe function */

  for (izone=0; izone < zones; izone++) {
    ith_zone = izone + 1;
    zone_nbr = cs_gui_boundary_zone_number(ith_zone);

    faces_list = cs_gui_get_faces_list(izone,
                                       boundaries->label[izone],
                                       *nfabor, *nozppm, &faces);

    if (cs_gui_strcmp(boundaries->nature[izone], "inlet")) {
      choice_v = _boundary_choice("inlet", boundaries->label[izone],
                                  "velocity_pressure", "choice");
      choice_d = _boundary_choice("inlet", boundaries->label[izone],
                                  "velocity_pressure", "direction");

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
        inmoxy[zone_nbr-1] = (int) boundaries->inmoxy[izone];
        ientcp[zone_nbr-1] = boundaries->ientcp[izone];
        qimpat[zone_nbr-1] = boundaries->qimp[izone];
        timpat[zone_nbr-1] = boundaries->timpat[izone];

        for (icharb = 0; icharb < *ncharb; icharb++) {
          qimpcp[icharb *(*nozppm)+zone_nbr-1]
            = boundaries->qimpcp[izone][icharb];
          timpcp[icharb *(*nozppm)+zone_nbr-1]
            = boundaries->timpcp[izone][icharb];

          for (iclass = 0; iclass < nclpch[icharb]; iclass++)
            distch[iclass * (*nozppm) * (*ncharm) +icharb * (*nozppm) +zone_nbr-1]
              = boundaries->distch[izone][icharb][iclass];
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
          mei_tree_insert(boundaries->velocity[izone], "t", *ttcabs);
          mei_tree_insert(boundaries->velocity[izone], "dt", *dtref);
          mei_tree_insert(boundaries->velocity[izone], "iter", *ntcabs);
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
        if (boundaries->itype[izone] == *iesicf) {
          cs_field_t *b_rho = CS_F_(rho_b);
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac] -1;
            rcodcl[(*ipr-1) * (*nfabor) + ifbr] = boundaries->prein[izone];
            rcodcl[(isca[*itempk-1]-1) * (*nfabor) + ifbr] = boundaries->tempin[izone];
            rcodcl[(isca[*ienerg-1]-1) * (*nfabor) + ifbr] = boundaries->entin[izone];
            b_rho->val[ifbr] = boundaries->rhoin[izone];
          }
        }
        if (boundaries->itype[izone] == *iephcf) {
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac] -1;
            rcodcl[(*ipr-1) * (*nfabor) + ifbr] = boundaries->prein[izone];
            rcodcl[(isca[*ienerg-1]-1) * (*nfabor) + ifbr] = boundaries->entin[izone];
          }
        }
      }
      else {
        if (   cs_gui_strcmp(choice_v, "flow1_formula")
            || cs_gui_strcmp(choice_v, "flow2_formula")) {
          mei_tree_insert(boundaries->velocity[izone], "t", *ttcabs);
          mei_tree_insert(boundaries->velocity[izone], "dt", *dtref);
          mei_tree_insert(boundaries->velocity[izone], "iter", *ntcabs);
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

      for (ifac = 0; ifac < faces; ifac++) {
        ifbr = faces_list[ifac]- 1;

        /* zone number and nature of boundary */
        izfppp[ifbr] = zone_nbr;
        if (cs_gui_strcmp(vars->model, "compressible_model"))
          itypfb[ifbr] = boundaries->itype[izone];
        else
          itypfb[ifbr] = *ientre;
      }

      /* dirichlet for turbulent variables and scalars */

      for (i = 0; i < vars->nvar; i++) {
        ivar = vars->rtp[i];
        int scal_id;
        int numvar  = vars->nvar - vars->nscaus - vars->nscapp;
        if (i < vars->nvar - vars->nscapp)
          scal_id = ivar - numvar;
        else
          scal_id = ivar - numvar + vars->nscaus;

        switch (boundaries->type_code[ivar][izone]) {

        case DIRICHLET:
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac]-1;
            rcodcl[ivar * (*nfabor) + ifbr]
              = boundaries->values[ivar][izone].val1;
          }
          break;

        case DIRICHLET_FORMULA:
          mei_tree_insert(boundaries->velocity[izone], "t", *ttcabs);
          mei_tree_insert(boundaries->velocity[izone], "dt", *dtref);
          mei_tree_insert(boundaries->velocity[izone], "iter", *ntcabs);
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac]-1;
            mei_evaluate(boundaries->scalar[ivar][izone]);
            rcodcl[ivar * (*nfabor) + ifbr]
              = mei_tree_lookup(boundaries->scalar[ivar][izone],
                                _scalar_label(scal_id));
          }
          break;
        }
      }

      if (cs_gui_strcmp(vars->model_value, "joule")) {
        if (*ielcor == 1) {
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac]- 1;
            rcodcl[(isca[*ipotr -1] -1) * (*nfabor) + ifbr] = boundaries->values[isca[*ipotr -1] -1][izone].val1 * (*coejou);
            if (*ipoti > 0)
              rcodcl[(isca[*ipoti -1] -1) * (*nfabor) + ifbr] = boundaries->values[isca[*ipoti -1] -1][izone].val1 * (*coejou);
          }
        }
      }

      if (cs_gui_strcmp(vars->model_value, "arc")) {
        if (boundaries->type_code[(isca[*ipotr -1] -1)][izone] == DIRICHLET_IMPLICIT && *ielcor == 1)
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac]- 1;
            icodcl[(isca[*ipotr-1]-1) *(*nfabor) + ifbr] = 5;
            rcodcl[(isca[*ipotr-1]-1) * (*nfabor) + ifbr] = *dpot;
          }
        if (boundaries->type_code[(isca[ipotva[0] -1] -1)][izone] == NEUMANN_IMPLICIT)
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac]- 1;
            int iel = ifabor[ifbr];
            icodcl[(isca[ipotva[0] -1] -1) *(*nfabor) + ifbr] = 5;
            icodcl[(isca[ipotva[1] -1] -1) *(*nfabor) + ifbr] = 5;
            icodcl[(isca[ipotva[2] -1] -1) *(*nfabor) + ifbr] = 5;
            rcodcl[(isca[ipotva[0] -1] -1) * (*nfabor) + ifbr] = rtpa[(isca[ipotva[0] -1] -1) * (*ncelet) + iel-1];
            rcodcl[(isca[ipotva[1] -1] -1) * (*nfabor) + ifbr] = rtpa[(isca[ipotva[1] -1] -1) * (*ncelet) + iel-1];
            rcodcl[(isca[ipotva[2] -1] -1) * (*nfabor) + ifbr] = rtpa[(isca[ipotva[2] -1] -1) * (*ncelet) + iel-1];
          }
      }

      if (cs_gui_strcmp(vars->model, "atmospheric_flows")) {
        iprofm[zone_nbr-1] = boundaries->meteo[izone].read_data;
        if (iprofm[zone_nbr-1] == 1) {
          BFT_FREE(choice_v);
          BFT_FREE(choice_d);
        }
        if (boundaries->meteo[izone].automatic) {
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac]- 1;
            itypfb[ifbr] = 0;
          }
        }
      }

      /* dirichlet for velocity */

      if (cs_gui_strcmp(choice_d, "coordinates")) {
        if (cs_gui_strcmp(choice_v, "norm")) {
          norm = boundaries->norm[izone] /
            sqrt(  boundaries->dirx[izone] * boundaries->dirx[izone]
                 + boundaries->diry[izone] * boundaries->diry[izone]
                 + boundaries->dirz[izone] * boundaries->dirz[izone]);

          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac] -1;
            rcodcl[vars->rtp[1] * (*nfabor) + ifbr] = boundaries->dirx[izone] * norm;
            rcodcl[vars->rtp[2] * (*nfabor) + ifbr] = boundaries->diry[izone] * norm;
            rcodcl[vars->rtp[3] * (*nfabor) + ifbr] = boundaries->dirz[izone] * norm;
          }
        }
        else if (   cs_gui_strcmp(choice_v, "flow1")
                 || cs_gui_strcmp(choice_v, "flow2")
                 || cs_gui_strcmp(choice_v, "flow1_formula")
                 || cs_gui_strcmp(choice_v, "flow2_formula")) {
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac] -1;
            rcodcl[vars->rtp[1] * (*nfabor) + ifbr] = boundaries->dirx[izone];
            rcodcl[vars->rtp[2] * (*nfabor) + ifbr] = boundaries->diry[izone];
            rcodcl[vars->rtp[3] * (*nfabor) + ifbr] = boundaries->dirz[izone];
          }
        }
        else if (cs_gui_strcmp(choice_v, "norm_formula")) {
          t0 = cs_timer_wtime();

          mei_tree_insert(boundaries->velocity[izone], "t", *ttcabs);
          mei_tree_insert(boundaries->velocity[izone], "dt", *dtref);
          mei_tree_insert(boundaries->velocity[izone], "iter", *ntcabs);

          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac] -1;

            mei_tree_insert(boundaries->velocity[izone], "x", cdgfbo[3 * ifbr + 0]);
            mei_tree_insert(boundaries->velocity[izone], "y", cdgfbo[3 * ifbr + 1]);
            mei_tree_insert(boundaries->velocity[izone], "z", cdgfbo[3 * ifbr + 2]);

            mei_evaluate(boundaries->velocity[izone]);

            norm =   mei_tree_lookup(boundaries->velocity[izone], "u_norm")
                   / sqrt(  boundaries->dirx[izone] * boundaries->dirx[izone]
                          + boundaries->diry[izone] * boundaries->diry[izone]
                          + boundaries->dirz[izone] * boundaries->dirz[izone]);

            rcodcl[vars->rtp[1] * (*nfabor) + ifbr] = boundaries->dirx[izone] * norm;
            rcodcl[vars->rtp[2] * (*nfabor) + ifbr] = boundaries->diry[izone] * norm;
            rcodcl[vars->rtp[3] * (*nfabor) + ifbr] = boundaries->dirz[izone] * norm;
          }

          cs_gui_add_mei_time(cs_timer_wtime() - t0);

        }
        if (cs_gui_strcmp(vars->model, "compressible_model")) {
          if (boundaries->itype[izone] == *iephcf) {
            for (ifac = 0; ifac < faces; ifac++) {
              ifbr = faces_list[ifac] -1;
              rcodcl[vars->rtp[1] * (*nfabor) + ifbr] = boundaries->dirx[izone];
              rcodcl[vars->rtp[2] * (*nfabor) + ifbr] = boundaries->diry[izone];
              rcodcl[vars->rtp[3] * (*nfabor) + ifbr] = boundaries->dirz[izone];
            }
          }
        }
      }
      else if (cs_gui_strcmp(choice_d, "normal")) {
        if (cs_gui_strcmp(choice_v, "norm")) {
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac] -1;

            norm =   boundaries->norm[izone]
                   / sqrt(  surfbo[3 * ifbr + 0] * surfbo[3 * ifbr + 0]
                          + surfbo[3 * ifbr + 1] * surfbo[3 * ifbr + 1]
                          + surfbo[3 * ifbr + 2] * surfbo[3 * ifbr + 2]);

            for (i = 1; i < 4; i++)
              rcodcl[vars->rtp[i] * (*nfabor) + ifbr] = -surfbo[3 * ifbr + vars->rtp[i]-1] * norm;
          }
        }
        else if (   cs_gui_strcmp(choice_v, "flow1")
                 || cs_gui_strcmp(choice_v, "flow2")
                 || cs_gui_strcmp(choice_v, "flow1_formula")
                 || cs_gui_strcmp(choice_v, "flow2_formula")) {
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac] -1;
            norm = sqrt(  surfbo[3 * ifbr + 0] * surfbo[3 * ifbr + 0]
                        + surfbo[3 * ifbr + 1] * surfbo[3 * ifbr + 1]
                        + surfbo[3 * ifbr + 2] * surfbo[3 * ifbr + 2]);
            rcodcl[vars->rtp[1] * (*nfabor) + ifbr] = -surfbo[3 * ifbr + vars->rtp[1] -1]/norm;
            rcodcl[vars->rtp[2] * (*nfabor) + ifbr] = -surfbo[3 * ifbr + vars->rtp[2] -1]/norm;
            rcodcl[vars->rtp[3] * (*nfabor) + ifbr] = -surfbo[3 * ifbr + vars->rtp[3] -1]/norm;
          }
        }
        else if (cs_gui_strcmp(choice_v, "norm_formula")) {
          t0 = cs_timer_wtime();

          mei_tree_insert(boundaries->velocity[izone], "t", *ttcabs);
          mei_tree_insert(boundaries->velocity[izone], "dt", *dtref);
          mei_tree_insert(boundaries->velocity[izone], "iter", *ntcabs);

          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac] -1;

            mei_tree_insert(boundaries->velocity[izone], "x", cdgfbo[3 * ifbr + 0]);
            mei_tree_insert(boundaries->velocity[izone], "y", cdgfbo[3 * ifbr + 1]);
            mei_tree_insert(boundaries->velocity[izone], "z", cdgfbo[3 * ifbr + 2]);

            mei_evaluate(boundaries->velocity[izone]);

            norm =   mei_tree_lookup(boundaries->velocity[izone], "u_norm")
                   / sqrt(  surfbo[3 * ifbr + 0] * surfbo[3 * ifbr + 0]
                          + surfbo[3 * ifbr + 1] * surfbo[3 * ifbr + 1]
                          + surfbo[3 * ifbr + 2] * surfbo[3 * ifbr + 2]);

            for (i = 1; i < 4; i++)
              rcodcl[vars->rtp[i] * (*nfabor) + ifbr] = -surfbo[3 * ifbr + vars->rtp[i]-1] * norm;
          }

          cs_gui_add_mei_time(cs_timer_wtime() - t0);

        }
        if (cs_gui_strcmp(vars->model, "compressible_model")) {
          if (boundaries->itype[izone] == *iephcf) {
            for (ifac = 0; ifac < faces; ifac++) {
              ifbr = faces_list[ifac] -1;
              for (i = 1; i < 4; i++)
                rcodcl[vars->rtp[i] * (*nfabor) + ifbr] = -surfbo[3 * ifbr + vars->rtp[i]-1];
            }
          }
        }
      }
      else if (cs_gui_strcmp(choice_d, "formula")) {
        t0 = cs_timer_wtime();

        mei_tree_insert(boundaries->direction[izone], "t", *ttcabs);
        mei_tree_insert(boundaries->direction[izone], "dt", *dtref);
        mei_tree_insert(boundaries->direction[izone], "iter", *ntcabs);

        if (cs_gui_strcmp(choice_v, "norm")) {
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac] -1;

            mei_tree_insert(boundaries->direction[izone], "x", cdgfbo[3 * ifbr + 0]);
            mei_tree_insert(boundaries->direction[izone], "y", cdgfbo[3 * ifbr + 1]);
            mei_tree_insert(boundaries->direction[izone], "z", cdgfbo[3 * ifbr + 2]);

            mei_evaluate(boundaries->direction[izone]);

            X[0] = mei_tree_lookup(boundaries->direction[izone], "dir_x");
            X[1] = mei_tree_lookup(boundaries->direction[izone], "dir_y");
            X[2] = mei_tree_lookup(boundaries->direction[izone], "dir_z");

            norm =   boundaries->norm[izone]
                   / sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2]);

            for (i = 1; i < 4; i++)
              rcodcl[vars->rtp[i] * (*nfabor) + ifbr] = X[i-1] * norm;
          }
        }
        else if (   cs_gui_strcmp(choice_v, "flow1")
                 || cs_gui_strcmp(choice_v, "flow2")
                 || cs_gui_strcmp(choice_v, "flow1_formula")
                 || cs_gui_strcmp(choice_v, "flow2_formula")) {
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac] -1;

            mei_tree_insert(boundaries->direction[izone], "x", cdgfbo[3 * ifbr + 0]);
            mei_tree_insert(boundaries->direction[izone], "y", cdgfbo[3 * ifbr + 1]);
            mei_tree_insert(boundaries->direction[izone], "z", cdgfbo[3 * ifbr + 2]);

            mei_evaluate(boundaries->direction[izone]);

            X[0] = mei_tree_lookup(boundaries->direction[izone], "dir_x");
            X[1] = mei_tree_lookup(boundaries->direction[izone], "dir_y");
            X[2] = mei_tree_lookup(boundaries->direction[izone], "dir_z");

            rcodcl[vars->rtp[1] * (*nfabor) + ifbr] = X[0];
            rcodcl[vars->rtp[2] * (*nfabor) + ifbr] = X[1];
            rcodcl[vars->rtp[3] * (*nfabor) + ifbr] = X[2];
          }
        }
        else if (cs_gui_strcmp(choice_v, "norm_formula")) {
          mei_tree_insert(boundaries->velocity[izone], "t", *ttcabs);
          mei_tree_insert(boundaries->velocity[izone], "dt", *dtref);
          mei_tree_insert(boundaries->velocity[izone], "iter", *ntcabs);

          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac] -1;

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

            norm =   mei_tree_lookup(boundaries->velocity[izone], "u_norm")
                   / sqrt( X[0] * X[0] + X[1] * X[1] + X[2] * X[2]);

            for (i = 1; i < 4; i++)
              rcodcl[vars->rtp[i] * (*nfabor) + ifbr] = X[i-1] * norm;
          }
        }

        if (cs_gui_strcmp(vars->model, "compressible_model")) {
          if (boundaries->itype[izone] == *iephcf) {
            for (ifac = 0; ifac < faces; ifac++) {
              ifbr = faces_list[ifac] -1;

              mei_tree_insert(boundaries->direction[izone], "x", cdgfbo[3 * ifbr + 0]);
              mei_tree_insert(boundaries->direction[izone], "y", cdgfbo[3 * ifbr + 1]);
              mei_tree_insert(boundaries->direction[izone], "z", cdgfbo[3 * ifbr + 2]);

              mei_evaluate(boundaries->direction[izone]);

              X[0] = mei_tree_lookup(boundaries->direction[izone], "dir_x");
              X[1] = mei_tree_lookup(boundaries->direction[izone], "dir_y");
              X[2] = mei_tree_lookup(boundaries->direction[izone], "dir_z");

              rcodcl[vars->rtp[1] * (*nfabor) + ifbr] = X[0];
              rcodcl[vars->rtp[2] * (*nfabor) + ifbr] = X[1];
              rcodcl[vars->rtp[3] * (*nfabor) + ifbr] = X[2];
            }
          }
        }

        cs_gui_add_mei_time(cs_timer_wtime() - t0);

      }
      BFT_FREE(choice_v);
      BFT_FREE(choice_d);

      /* turbulent inlet, with formula */
      if (icalke[zone_nbr-1] == 0) {
        t0 = cs_timer_wtime();
        path1 = cs_xpath_init_path();
        cs_xpath_add_elements(&path1, 2, "boundary_conditions", "inlet");
        cs_xpath_add_test_attribute(&path1, "label", boundaries->label[izone]);
        cs_xpath_add_element(&path1, "turbulence");
        cs_xpath_add_element(&path1, "formula");
        cs_xpath_add_function_text(&path1);
        formula = cs_gui_get_text_value(path1);
        if (formula != NULL) {
          ev_formula = mei_tree_new(formula);
          mei_tree_insert(ev_formula,"x",0.0);
          mei_tree_insert(ev_formula,"y",0.0);
          mei_tree_insert(ev_formula,"z",0.0);
          mei_tree_insert(ev_formula, "t", *ttcabs);
          mei_tree_insert(ev_formula, "dt", *dtref);
          mei_tree_insert(ev_formula, "iter", *ntcabs);
          /* try to build the interpreter */

          if (mei_tree_builder(ev_formula))
            bft_error(__FILE__, __LINE__, 0,
                      _("Error: can not interpret expression: %s\n %i"),
                      ev_formula->string, mei_tree_builder(ev_formula));

          model = cs_gui_get_thermophysical_model("turbulence");
          if (model == NULL)
            return;

          if (   cs_gui_strcmp(model, "k-epsilon")
              || cs_gui_strcmp(model, "k-epsilon-PL")) {
            const char *symbols[] = {"k","eps"};
            if (mei_tree_find_symbols(ev_formula, 2, symbols))
              bft_error(__FILE__, __LINE__, 0,
                        _("Error: can not find the required symbol: %s\n"),
                        "k or eps");

            for (ifac = 0; ifac < faces; ifac++) {
              ifbr = faces_list[ifac] - 1;
              mei_tree_insert(ev_formula, "x", cdgfbo[3 * ifbr + 0]);
              mei_tree_insert(ev_formula, "y", cdgfbo[3 * ifbr + 1]);
              mei_tree_insert(ev_formula, "z", cdgfbo[3 * ifbr + 2]);
              mei_evaluate(ev_formula);
              rcodcl[vars->rtp[4] * (*nfabor) + ifbr] = mei_tree_lookup(ev_formula, "k");
              rcodcl[vars->rtp[5] * (*nfabor) + ifbr] = mei_tree_lookup(ev_formula, "eps");
            }
          }
          else if (  cs_gui_strcmp(model, "Rij-epsilon")
                   ||cs_gui_strcmp(model, "Rij-SSG")) {
            const char *symbols[] = {"R11", "R22", "R133", "R12", "R13", "R23", "eps"};
            if (mei_tree_find_symbols(ev_formula, 7, symbols))
              bft_error(__FILE__, __LINE__, 0, _("Error: can not find the required symbol: %s\n"), "R11, R22, R33, R12, R13, R23 or eps");
            for (ifac = 0; ifac < faces; ifac++) {
              ifbr = faces_list[ifac]- 1;
              mei_tree_insert(ev_formula, "x", cdgfbo[3 * ifbr + 0]);
              mei_tree_insert(ev_formula, "y", cdgfbo[3 * ifbr + 1]);
              mei_tree_insert(ev_formula, "z", cdgfbo[3 * ifbr + 2]);
              mei_evaluate(ev_formula);
              rcodcl[vars->rtp[4] * (*nfabor) + ifbr]  = mei_tree_lookup(ev_formula, "R11");
              rcodcl[vars->rtp[5] * (*nfabor) + ifbr]  = mei_tree_lookup(ev_formula, "R22");
              rcodcl[vars->rtp[6] * (*nfabor) + ifbr]  = mei_tree_lookup(ev_formula, "R33");
              rcodcl[vars->rtp[7] * (*nfabor) + ifbr]  = mei_tree_lookup(ev_formula, "R12");
              rcodcl[vars->rtp[8] * (*nfabor) + ifbr]  = mei_tree_lookup(ev_formula, "R13");
              rcodcl[vars->rtp[9] * (*nfabor) + ifbr]  = mei_tree_lookup(ev_formula, "R23");
              rcodcl[vars->rtp[10] * (*nfabor) + ifbr] = mei_tree_lookup(ev_formula, "eps");
            }
          }
          else if (cs_gui_strcmp(model, "Rij-EBRSM")) {
            const char *symbols[] = {"R11","R22","R133","R12","R13","R23","eps", "alpha"};
            if (mei_tree_find_symbols(ev_formula, 8, symbols))
              bft_error(__FILE__, __LINE__, 0, _("Error: can not find the required symbol: %s\n"), "R11, R22, R33, R12, R13, R23, eps or alpha");
            for (ifac = 0; ifac < faces; ifac++) {
              ifbr = faces_list[ifac]- 1;
              mei_tree_insert(ev_formula, "x", cdgfbo[3 * ifbr + 0]);
              mei_tree_insert(ev_formula, "y", cdgfbo[3 * ifbr + 1]);
              mei_tree_insert(ev_formula, "z", cdgfbo[3 * ifbr + 2]);
              mei_evaluate(ev_formula);
              rcodcl[vars->rtp[4] * (*nfabor) + ifbr]  = mei_tree_lookup(ev_formula,"R11");
              rcodcl[vars->rtp[5] * (*nfabor) + ifbr]  = mei_tree_lookup(ev_formula,"R22");
              rcodcl[vars->rtp[6] * (*nfabor) + ifbr]  = mei_tree_lookup(ev_formula,"R33");
              rcodcl[vars->rtp[7] * (*nfabor) + ifbr]  = mei_tree_lookup(ev_formula,"R12");
              rcodcl[vars->rtp[8] * (*nfabor) + ifbr]  = mei_tree_lookup(ev_formula,"R13");
              rcodcl[vars->rtp[9] * (*nfabor) + ifbr]  = mei_tree_lookup(ev_formula,"R23");
              rcodcl[vars->rtp[10] * (*nfabor) + ifbr] = mei_tree_lookup(ev_formula,"eps");
              rcodcl[vars->rtp[11] * (*nfabor) + ifbr] = mei_tree_lookup(ev_formula,"alpha");
            }
          }
          else if (cs_gui_strcmp(model, "v2f-BL-v2/k")) {
            const char *symbols[] = {"k", "eps", "phi", "alpha"};

            if (mei_tree_find_symbols(ev_formula, 4, symbols))
              bft_error(__FILE__, __LINE__, 0,
                        _("Error: can not find the required symbol: %s\n"),
                        "k, eps, phi of alpha");

            for (ifac = 0; ifac < faces; ifac++) {
              ifbr = faces_list[ifac]- 1;
              mei_tree_insert(ev_formula, "x", cdgfbo[3 * ifbr + 0]);
              mei_tree_insert(ev_formula, "y", cdgfbo[3 * ifbr + 1]);
              mei_tree_insert(ev_formula, "z", cdgfbo[3 * ifbr + 2]);
              mei_evaluate(ev_formula);
              rcodcl[vars->rtp[4] * (*nfabor) + ifbr] = mei_tree_lookup(ev_formula, "k");
              rcodcl[vars->rtp[5] * (*nfabor) + ifbr] = mei_tree_lookup(ev_formula, "eps");
              rcodcl[vars->rtp[6] * (*nfabor) + ifbr] = mei_tree_lookup(ev_formula, "phi");
              rcodcl[vars->rtp[7] * (*nfabor) + ifbr] = mei_tree_lookup(ev_formula, "alpha");
            }
          }
          else if (cs_gui_strcmp(model, "k-omega-SST")) {
            const char *symbols[] = {"k", "omega"};

            if (mei_tree_find_symbols(ev_formula, 2, symbols))
              bft_error(__FILE__, __LINE__, 0,
                        _("Error: can not find the required symbol: %s\n"),
                        "k or omega");

            for (ifac = 0; ifac < faces; ifac++) {
              ifbr = faces_list[ifac]- 1;
              mei_tree_insert(ev_formula, "x", cdgfbo[3 * ifbr + 0]);
              mei_tree_insert(ev_formula, "y", cdgfbo[3 * ifbr + 1]);
              mei_tree_insert(ev_formula, "z", cdgfbo[3 * ifbr + 2]);
              mei_evaluate(ev_formula);
              rcodcl[vars->rtp[4] * (*nfabor) + ifbr] = mei_tree_lookup(ev_formula, "k");
              rcodcl[vars->rtp[5] * (*nfabor) + ifbr] = mei_tree_lookup(ev_formula, "omega");
            }
          }
          else if (cs_gui_strcmp(model, "Spalart-Allmaras")) {
            const char *symbols[] = {"nusa"};

            if (mei_tree_find_symbols(ev_formula, 1, symbols))
              bft_error(__FILE__, __LINE__, 0,
                        _("Error: can not find the required symbol: %s\n"),
                        "nusa");

            for (ifac = 0; ifac < faces; ifac++) {
              ifbr = faces_list[ifac]- 1;
              mei_tree_insert(ev_formula, "x", cdgfbo[3 * ifbr + 0]);
              mei_tree_insert(ev_formula, "y", cdgfbo[3 * ifbr + 1]);
              mei_tree_insert(ev_formula, "z", cdgfbo[3 * ifbr + 2]);
              mei_evaluate(ev_formula);
              rcodcl[vars->rtp[4] * (*nfabor) + ifbr] = mei_tree_lookup(ev_formula, "nusa");
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
      if (boundaries->rough[izone] >= 0.0) {
        /* roughness value is only stored in Velocity_U (z0) */
        /* Remember: rcodcl(ifac, ivar, 1) -> rcodcl[k][j][i]
                   = rcodcl[ k*dim1*dim2 + j*dim1 + i] */
        iwall = *iparug;
        for (ifac = 0; ifac < faces; ifac++) {
          ifbr = faces_list[ifac] -1;
          rcodcl[2 * (*nfabor * (vars->nvar)) + vars->rtp[1] * (*nfabor) + ifbr]
            = boundaries->rough[izone];

          /* if atmospheric flows and "dry" or "humid atmosphere" mode,
             roughness value also stored in Velocity_V (z0t)*/
          if (cs_gui_strcmp(vars->model, "atmospheric_flows")) {
            if (   cs_gui_strcmp(vars->model_value, "dry")
                || cs_gui_strcmp(vars->model_value, "humid")) {
              rcodcl[2 * (*nfabor * (vars->nvar)) + vars->rtp[2] * (*nfabor) + ifbr]
                = boundaries->rough[izone];
            }
          }
        }
      }
      else {
        iwall = *iparoi;
      }

      for (ifac = 0; ifac < faces; ifac++) {
        ifbr = faces_list[ifac]-1;
        izfppp[ifbr] = zone_nbr;
        itypfb[ifbr] = iwall;
      }

      for (i = 0; i < vars->nvar; i++) {
        ivar = vars->rtp[i];

        switch (boundaries->type_code[ivar][izone]) {

        case NEUMANN:
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac]-1;
            icodcl[ivar *(*nfabor) + ifbr] = 3;
            rcodcl[2 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
              = boundaries->values[ivar][izone].val3;
          }
          break;

        case DIRICHLET:
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac]-1;
            icodcl[ivar *(*nfabor) + ifbr] = 5;
            /* if wall_function <-- icodcl[ivar *(*nfabor) + ifbr] = 1; */
            rcodcl[0 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
              = boundaries->values[ivar][izone].val1;
          }
          break;

        case WALL_FUNCTION:
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac]-1;
            icodcl[ivar *(*nfabor) + ifbr] = 5;
            rcodcl[0 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
              = boundaries->values[ivar][izone].val1;
          }
          break;

        case EXCHANGE_COEFF:
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac]-1;
            icodcl[ivar *(*nfabor) + ifbr] = 5;
            rcodcl[0 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
              = boundaries->values[ivar][izone].val1;
            rcodcl[1 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
              = boundaries->values[ivar][izone].val2;
          }
          break;

        case DIRICHLET_FORMULA:
          mei_tree_insert(boundaries->velocity[izone], "t", *ttcabs);
          mei_tree_insert(boundaries->velocity[izone], "dt", *dtref);
          mei_tree_insert(boundaries->velocity[izone], "iter", *ntcabs);
          for (ifac = 0; ifac < faces; ifac++) {
            int scal_id;
            int numvar  = vars->nvar - vars->nscaus - vars->nscapp;
            if (i < vars->nvar - vars->nscapp)
              scal_id = ivar - numvar;
            else
              scal_id = ivar - numvar + vars->nscaus;
            ifbr = faces_list[ifac]-1;
            mei_evaluate(boundaries->scalar[ivar][izone]);
            icodcl[ivar *(*nfabor) + ifbr] = 5;
            rcodcl[0 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
              = mei_tree_lookup(boundaries->scalar[ivar][izone],
                                _scalar_label(scal_id));
          }
          break;

        case NEUMANN_FORMULA:
          mei_tree_insert(boundaries->velocity[izone], "t", *ttcabs);
          mei_tree_insert(boundaries->velocity[izone], "dt", *dtref);
          mei_tree_insert(boundaries->velocity[izone], "iter", *ntcabs);
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac]-1;
            mei_evaluate(boundaries->scalar[ivar][izone]);
            icodcl[ivar *(*nfabor) + ifbr] = 3;
            rcodcl[2 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
              = mei_tree_lookup(boundaries->scalar[ivar][izone], "flux");
          }
          break;

        case EXCHANGE_COEFF_FORMULA:
          mei_tree_insert(boundaries->velocity[izone], "t", *ttcabs);
          mei_tree_insert(boundaries->velocity[izone], "dt", *dtref);
          mei_tree_insert(boundaries->velocity[izone], "iter", *ntcabs);
          for (ifac = 0; ifac < faces; ifac++) {
            int scal_id;
            int numvar  = vars->nvar - vars->nscaus - vars->nscapp;
            if (i < vars->nvar - vars->nscapp)
              scal_id = ivar - numvar;
            else
              scal_id = ivar - numvar + vars->nscaus;
            ifbr = faces_list[ifac]-1;
            mei_evaluate(boundaries->scalar[ivar][izone]);
            icodcl[ivar *(*nfabor) + ifbr] = 5;
            rcodcl[0 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
              = mei_tree_lookup(boundaries->scalar[ivar][izone],
                                _scalar_label(scal_id));
            rcodcl[1 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
              = mei_tree_lookup(boundaries->scalar[ivar][izone], "hc");
          }
          break;
        } /* switch */
      } /* for (i = 0; i < vars->nvar; i++) */

      if (cs_gui_strcmp(vars->model_value, "joule")) {
        if (*ielcor == 1) {
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac]- 1;
            rcodcl[(isca[*ipotr -1] -1) * (*nfabor) + ifbr] = boundaries->values[isca[*ipotr -1] -1][izone].val1 * (*coejou);
            if (*ipoti > 0)
              rcodcl[(isca[*ipoti -1] -1) * (*nfabor) + ifbr] = boundaries->values[isca[*ipoti -1] -1][izone].val1 * (*coejou);
          }
        }
      }

      if (cs_gui_strcmp(vars->model_value, "arc")) {
        if (    boundaries->type_code[(isca[*ipotr-1]-1)][izone] == DIRICHLET_IMPLICIT
            && *ielcor == 1)
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac]- 1;
            icodcl[(isca[*ipotr-1]-1) *(*nfabor) + ifbr] = 5;
            rcodcl[(isca[*ipotr-1]-1) * (*nfabor) + ifbr] = *dpot;
          }
        if (boundaries->type_code[(isca[ipotva[0] -1] -1)][izone] == NEUMANN_IMPLICIT)
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac]- 1;
            int iel = ifabor[ifbr];
            icodcl[(isca[ipotva[0] -1] -1) *(*nfabor) + ifbr] = 5;
            icodcl[(isca[ipotva[1] -1] -1) *(*nfabor) + ifbr] = 5;
            icodcl[(isca[ipotva[2] -1] -1) *(*nfabor) + ifbr] = 5;
            rcodcl[(isca[ipotva[0] -1] -1) * (*nfabor) + ifbr] = rtpa[(isca[ipotva[0] -1] -1) * (*ncelet) + iel-1];
            rcodcl[(isca[ipotva[1] -1] -1) * (*nfabor) + ifbr] = rtpa[(isca[ipotva[1] -1] -1) * (*ncelet) + iel-1];
            rcodcl[(isca[ipotva[2] -1] -1) * (*nfabor) + ifbr] = rtpa[(isca[ipotva[2] -1] -1) * (*ncelet) + iel-1];
          }
      }
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "outlet")) {

      for (ifac = 0; ifac < faces; ifac++) {
        ifbr = faces_list[ifac]-1;
        izfppp[ifbr] = zone_nbr;
        if (cs_gui_strcmp(vars->model, "compressible_model"))
          itypfb[ifbr] = boundaries->itype[izone];
        else
          itypfb[ifbr] = *isolib;
      }

      for (i = 0; i < vars->nvar; i++) {
        ivar = vars->rtp[i];

        switch (boundaries->type_code[ivar][izone]) {

        case DIRICHLET:
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac]-1;
            icodcl[ivar *(*nfabor) + ifbr] = 1;
            rcodcl[0 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
              = boundaries->values[ivar][izone].val1;
          }
          break;

        case DIRICHLET_FORMULA:
          mei_tree_insert(boundaries->velocity[izone], "t", *ttcabs);
          mei_tree_insert(boundaries->velocity[izone], "dt", *dtref);
          mei_tree_insert(boundaries->velocity[izone], "iter", *ntcabs);
          for (ifac = 0; ifac < faces; ifac++) {
            int scal_id;
            int numvar  = vars->nvar - vars->nscaus - vars->nscapp;
            if (i < vars->nvar - vars->nscapp)
              scal_id = ivar - numvar;
            else
              scal_id = ivar - numvar + vars->nscaus;
            ifbr = faces_list[ifac]-1;
            mei_evaluate(boundaries->scalar[ivar][izone]);
            icodcl[ivar *(*nfabor) + ifbr] = 1;
            rcodcl[0 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
              = mei_tree_lookup(boundaries->scalar[ivar][izone],
                                _scalar_label(scal_id));
          }
          break;
        }
      }

      if (cs_gui_strcmp(vars->model, "atmospheric_flows")) {
        iprofm[zone_nbr-1] = boundaries->meteo[izone].read_data;
        if (boundaries->meteo[izone].automatic) {
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac]- 1;
            itypfb[ifbr] = 0;
          }
        }
      }
      else if (cs_gui_strcmp(vars->model, "compressible_model")) {
        if (boundaries->itype[izone] == *isopcf) {
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac] -1;
            rcodcl[(*ipr-1) * (*nfabor) + ifbr] = boundaries->preout[izone];
          }
        }
      }

      if (cs_gui_strcmp(vars->model_value, "joule")) {
        if (*ielcor == 1) {
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac]- 1;
            rcodcl[(isca[*ipotr -1] -1) * (*nfabor) + ifbr] = boundaries->values[isca[*ipotr -1] -1][izone].val1 * (*coejou);
            if (*ipoti > 0)
              rcodcl[(isca[*ipoti -1] -1) * (*nfabor) + ifbr] = boundaries->values[isca[*ipoti -1] -1][izone].val1 * (*coejou);
          }
        }
      }

      if (cs_gui_strcmp(vars->model_value, "arc")) {
        if (boundaries->type_code[(isca[*ipotr-1]-1)][izone] == DIRICHLET_IMPLICIT && *ielcor == 1)
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac]- 1;
            icodcl[(isca[*ipotr-1]-1) *(*nfabor) + ifbr] = 5;
            rcodcl[(isca[*ipotr-1]-1) * (*nfabor) + ifbr] = *dpot;
          }
        if (boundaries->type_code[(isca[ipotva[0] -1] -1)][izone] == NEUMANN_IMPLICIT)
          for (ifac = 0; ifac < faces; ifac++) {
            ifbr = faces_list[ifac]- 1;
            int iel = ifabor[ifbr];
            icodcl[(isca[ipotva[0] -1] -1) *(*nfabor) + ifbr] = 5;
            icodcl[(isca[ipotva[1] -1] -1) *(*nfabor) + ifbr] = 5;
            icodcl[(isca[ipotva[2] -1] -1) *(*nfabor) + ifbr] = 5;
            rcodcl[(isca[ipotva[0] -1] -1) * (*nfabor) + ifbr] = rtpa[(isca[ipotva[0] -1] -1) * (*ncelet) + iel-1];
            rcodcl[(isca[ipotva[1] -1] -1) * (*nfabor) + ifbr] = rtpa[(isca[ipotva[1] -1] -1) * (*ncelet) + iel-1];
            rcodcl[(isca[ipotva[2] -1] -1) * (*nfabor) + ifbr] = rtpa[(isca[ipotva[2] -1] -1) * (*ncelet) + iel-1];
          }
      }
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "symmetry")) {
      for (ifac = 0; ifac < faces; ifac++) {
        ifbr = faces_list[ifac]-1;
        izfppp[ifbr] = zone_nbr;
        itypfb[ifbr] = *isymet;
      }

    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "undefined")) {
      for (ifac = 0; ifac < faces; ifac++) {
        ifbr = faces_list[ifac]-1;
        izfppp[ifbr] = zone_nbr;
        itypfb[ifbr] = *iindef;
      }

    }
    else {
      bft_error(__FILE__, __LINE__, 0,
                _("boundary nature %s is unknown \n"),
                boundaries->nature[izone]);
    }
    BFT_FREE(faces_list);

  } /*  for (izone=0 ; izone < zones ; izone++) */

#if _XML_DEBUG_
  bft_printf("==>UICLIM\n");
  bft_printf("--boundary zones number: %i\n", zones);

  for (izone = 0 ; izone < zones ; izone++) {

    faces_list = cs_gui_get_faces_list(izone,
                                       boundaries->label[izone],
                                       *nfabor, *nozppm, &faces);

    zone_nbr = cs_gui_boundary_zone_number(izone+1);

    bft_printf("\n---zone %i label: %s\n", zone_nbr, boundaries->label[izone]);
    bft_printf("---zone %i nature: %s\n", zone_nbr, boundaries->nature[izone]);
    bft_printf("---zone %i number of faces: %i\n", zone_nbr, faces);
    bft_printf("----localization: %s\n", cs_gui_boundary_zone_localization(boundaries->label[izone]));

    if (cs_gui_strcmp(boundaries->nature[izone], "inlet")) {
      choice_v = _boundary_choice("inlet", boundaries->label[izone], "velocity_pressure", "choice");
      choice_d = _boundary_choice("inlet", boundaries->label[izone], "velocity_pressure", "direction");

      if (cs_gui_strcmp(choice_v, "norm"))
        bft_printf("-----velocity: %s => %12.5e \n", choice_v, boundaries->norm[izone]);
      if (cs_gui_strcmp(choice_v, "flow1") || cs_gui_strcmp(choice_v, "flow2"))
        bft_printf("-----velocity: %s => %12.5e \n", choice_v, boundaries->qimp[izone]);
      if (   cs_gui_strcmp(choice_v, "norm_formula")
          || cs_gui_strcmp(choice_v, "flow1_formula")
          || cs_gui_strcmp(choice_v, "flow2_formula"))
        bft_printf("-----velocity: %s => %s \n",
                   choice_v, boundaries->velocity[izone]->string);
      if (cs_gui_strcmp(choice_d, "coordinates"))
        bft_printf("-----direction: %s => %12.5e %12.5e %12.5e \n",
                   choice_v, boundaries->dirx[izone],
                   boundaries->diry[izone], boundaries->dirz[izone]);
      if (cs_gui_strcmp(choice_d, "formula"))
        bft_printf("-----direction: %s => %s \n", choice_d,
                   boundaries->direction[izone]->string);
      BFT_FREE(choice_v);
      BFT_FREE(choice_d);

      if (cs_gui_strcmp(vars->model, "solid_fuels")) {
        bft_printf("-----iqimp=%i, qimpat=%12.5e \n",
                   iqimp[zone_nbr-1], qimpat[zone_nbr-1]);
        bft_printf("-----ientat=%i, ientcp=%i, timpat=%12.5e \n",
                   ientat[zone_nbr-1], ientcp[zone_nbr-1], timpat[zone_nbr-1]);

        for (icharb = 0; icharb < *ncharb; icharb++) {
          bft_printf("-----coal=%i, qimpcp=%12.5e, timpcp=%12.5e \n",
                     icharb+1, qimpcp[icharb *(*nozppm)+zone_nbr-1],
                     timpcp[icharb *(*nozppm)+zone_nbr-1]);

          for (iclass = 0; iclass < nclpch[icharb]; k++)
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
        if (boundaries->itype[izone] == *iesicf) {
          bft_printf("-----imposed_inlet\n");
          bft_printf("-----premin=%i \n",boundaries->prein[zone_nbr-1]);
          bft_printf("-----rhoin=%i \n",boundaries->rhoin[zone_nbr-1]);
          bft_printf("-----tempin=%i \n",boundaries->tempin[zone_nbr-1]);
          bft_printf("-----entin=%i \n",boundaries->entin[zone_nbr-1]);
        }
        if (boundaries->itype[izone] == *iephcf) {
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
      ifbr = faces_list[0]-1;

      for (i = 0; i < vars->nvar ; i++) {
        ivar = vars->rtp[i];
        bft_printf("------%s: icodcl=%i, "
                   "rcodcl(1)=%12.5e, rcodcl(2)=%12.5e, rcodcl(3)=%12.5e\n",
                   vars->name[ivar],
                   icodcl[ivar *(*nfabor) +ifbr ],
                   rcodcl[0 * (*nfabor * (vars->nvar)) +ivar * (*nfabor) +ifbr],
                   rcodcl[1 * (*nfabor * (vars->nvar)) +ivar * (*nfabor) +ifbr],
                   rcodcl[2 * (*nfabor * (vars->nvar)) +ivar * (*nfabor) +ifbr]);
      }
    }
    BFT_FREE(faces_list);
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
 * INTEGER          NFABOR  <-- number of boundary faces
 * INTEGER          NOZPPM  <-- max number of boundary conditions zone
 * INTEGER          IINDEF  <-- type of boundary: not defined
 * INTEGER          IENTRE  <-- type of boundary: inlet
 * INTEGER          IESICF  --> type of boundary: imposed inlet (compressible)
 * INTEGER          ISSPCF  --> type of boundary: supersonic outlet (compressible)
 * INTEGER          IEPHCF  --> type of boundary: subsonic inlet at given total pressure
 *                                                and enthalpy (compressible)
 * INTEGER          ISOPCF  --> type of boundary: subsonic outlet (compressible)
 * INTEGER          IPAROI  <-- type of boundary: smooth wall
 * INTEGER          IPARUG  <-- type of boundary: rough wall
 * INTEGER          ISYMET  <-- type of boundary: symetry
 * INTEGER          ISOLIB  <-- type of boundary: outlet
 * INTEGER          ITYPFB  <-- type of boundary for each face
 * INTEGER          IZFPPP  <-- zone number
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiclve, UICLVE)(const int  *nfabor,
                               const int  *nozppm,
                               const int  *iindef,
                               const int  *ientre,
                               const int  *iesicf,
                               const int  *iephcf,
                               const int  *isspcf,
                               const int  *isopcf,
                               const int  *iparoi,
                               const int  *iparug,
                               const int  *isymet,
                               const int  *isolib,
                                     int  *itypfb,
                                     int  *izfppp)
{
  int ifbr, ifac;
  int izone, zones, zone_nbr;
  int inature = 0;
  int inature2 = 0;
  int *faces_list = NULL;
  int faces = 0;

  zones   = cs_gui_boundary_zones_number();

  for (izone=0 ; izone < zones ; izone++) {
    if (cs_gui_strcmp(boundaries->nature[izone], "inlet")) {
      inature = *ientre;
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "wall")) {
      inature = *iparug;
      if (boundaries->rough[izone] < 0.0)
        inature = *iparoi;
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "outlet")) {
      inature = *isolib;
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "symmetry")) {
      inature = *isymet;
    }
    else if (cs_gui_strcmp(boundaries->nature[izone], "undefined")) {
      inature = *iindef;
    }
    else
      bft_error(__FILE__, __LINE__, 0,
                _("boundary nature %s is unknown \n"),
                boundaries->nature[izone]);

    zone_nbr =  cs_gui_boundary_zone_number(izone + 1);
    faces_list = cs_gui_get_faces_list(izone,
                                       boundaries->label[izone],
                                       *nfabor, *nozppm, &faces);
    for (ifac = 0; ifac < faces; ifac++) {
      ifbr = faces_list[ifac]-1;

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

      /* The nature of the boundary can be changed from smooth wall to
         rough wall or vice-versa between the GUI and the FORTRAN */
      if (inature2 == *iparug ) inature2 = *iparoi;
      if (inature == *iparug ) inature = *iparoi;

      if (cs_gui_strcmp(cs_glob_var->model, "atmospheric_flows")) {
        if (boundaries->meteo[izone].automatic) {
          if ((inature == *ientre || inature == *isolib) && inature2 == 0)
            inature2 = inature;
        }
      }
      else if (cs_gui_strcmp(cs_glob_var->model, "compressible_model")) {
        if (inature == *ientre  && inature2 == *iesicf)
          inature2 = inature;
        if (inature == *ientre  && inature2 == *iephcf)
          inature2 = inature;
        if (inature == *isolib  && inature2 == *isspcf)
          inature2 = inature;
        if (inature == *isolib  && inature2 == *isopcf)
          inature2 = inature;
      }

      if (inature2 != inature)
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
    BFT_FREE(faces_list);
  } /*  for izone */
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Free boundary conditions structures
 *
 * parameters:
 *   ncharb  <-- number of coals
 *----------------------------------------------------------------------------*/

void
cs_gui_boundary_conditions_free_memory(const int  *ncharb)
{
  int i;
  int ivar;
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
      for (i = 0 ; i < vars->nvar ; i++) {
        ivar = vars->rtp[i];
        mei_tree_destroy(boundaries->scalar[ivar][izone]);
      }
    }

    for (i=0; i < vars->nvar; i++) {
      ivar = vars->rtp[i];
      BFT_FREE(boundaries->type_code[ivar]);
      BFT_FREE(boundaries->values[ivar]);
      BFT_FREE(boundaries->scalar[ivar]);
    }

    if (cs_gui_strcmp(vars->model, "solid_fuels")) {
      for (izone=0 ; izone < zones ; izone++) {
        BFT_FREE(boundaries->qimpcp[izone]);
        BFT_FREE(boundaries->timpcp[izone]);
        for (icharb=0; icharb < *ncharb; icharb++)
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
      BFT_FREE(boundaries->preout);
    }
    if (cs_gui_strcmp(vars->model, "atmospheric_flows"))
        BFT_FREE(boundaries->meteo);

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
    BFT_FREE(boundaries->scalar);
    BFT_FREE(boundaries);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
