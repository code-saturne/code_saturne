/*============================================================================
 * Notebook management.
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_printf.h"

#include "cs_gui_util.h"
#include "cs_parameters.h"


/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_notebook.h"

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/


/*============================================================================
 * Static global variables
 *============================================================================*/

static cs_tree_node_t *_notebook = NULL;

/*----------------------------------------------------------------------------*/
/*!
 *  \brief Initialize the notebook object (based on cs_tree_node_t)
 *
 *  The name used to identify the object is "cs_notebook".
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_notebook_load_from_file(void)
{
  _notebook = cs_tree_get_node(cs_glob_tree, "physical_properties/notebook/var");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a parameter value (real).
 *
 * The name used is the same as the one in the GUI
 *
 * \param[in] name  name of the parameter
 *
 * \return value of the given parameter
 *
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_notebook_parameter_value_by_name(const char *name)
{

  cs_tree_node_t *n = cs_tree_node_get_sibling_with_tag(_notebook, "name", name);

  if (n == NULL) {
    bft_error(__FILE__, __LINE__, 0,
              "Variable %s was not defined in the notebook\n",
              name);
  }

  const char *c_value = cs_tree_node_get_tag(n, "value");

  cs_real_t val = atof(c_value);
  return val;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Indicate whether the parameter is used for a study with openturns.
 *
 * Returns an int flag to indicate whether this paramter is used for an
 * OpenTurns study.
 *  -1 : The parameter is not used with OpenTurns
 *   0 : The parameter is used as an input from OpenTurns
 *   1 : The parameter is used as an output to OpenTurns
 *
 * \param[in] name  name of the parameter
 *
 * \return  an int flag value
 */
/*----------------------------------------------------------------------------*/

int
cs_notebook_parameter_get_openturns_status(char *name)
{
  cs_tree_node_t *n = cs_tree_node_get_sibling_with_tag(_notebook,
                                                        "name",
                                                        name);

  if (n == NULL) {
    bft_error(__FILE__, __LINE__, 0,
              "Variable %s was not defined in the notebook\n",
              name);
  }

  const char *c_value = cs_tree_node_get_tag(n, "oturns");

  int status = -1;
  if (c_value != NULL) {
    if (strcmp(c_value, "Yes: Input") == 0)
      status = 0;
    else if (strcmp(c_value, "Yes: Output") == 0)
      status = 1;
  }

  return status;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Returns the description of the parameter (GUI defined)
 *
 * \param[in] name  name of the parameter
 *
 * \return  a const char pointer containing the description.
 *
 */
/*----------------------------------------------------------------------------*/

const char *
cs_notebook_parameter_get_description(char *name)
{
  cs_tree_node_t *n = cs_tree_node_get_sibling_with_tag(_notebook,
                                                        "name",
                                                        name);

  if (n == NULL) {
    bft_error(__FILE__, __LINE__, 0,
              "Variable %s was not defined in the notebook\n",
              name);
  }

  const char *d = cs_tree_node_get_tag(n, "description");
  return d;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
