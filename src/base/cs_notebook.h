#ifndef __CS_NOTEBOOK_H__
#define __CS_NOTEBOOK_H__

/*============================================================================
 * Notebook management.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*
 * \brief Initialize the notebook object (based on cs_tree_node_t).
 */
/*----------------------------------------------------------------------------*/

void
cs_notebook_load_from_file(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Check if a parameter value is present.
 *
 * \param[in]   name      name of the parameter
 * \param[out]  editable  1 if the value is editable, 0 otherwise (optional)
 *
 * \return 0 if not present, 1 if present
 */
/*----------------------------------------------------------------------------*/

int
cs_notebook_parameter_is_present(const char  *name,
                                 int         *editable);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return a parameter value (real).
 *
 * The name used is the same as the one in the GUI.
 *
 * \param[in] name  name of the parameter
 *
 * \return value of the given parameter
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_notebook_parameter_value_by_name(const char *name);

/*----------------------------------------------------------------------------*/
/*
 * \brief Set a parameter value (real) for an editable parameter.
 *
 * The name used is the same as the one in the GUI.
 *
 * \param[in] name  name of the parameter
 * \param[in] val   value of the parameter
 */
/*----------------------------------------------------------------------------*/

void
cs_notebook_parameter_set_value(const char *name,
                                cs_real_t   val);

/*----------------------------------------------------------------------------*/
/*
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
cs_notebook_parameter_get_openturns_status(char *name);

/*----------------------------------------------------------------------------*/
/*
 * \brief Returns the description of the parameter (GUI defined).
 *
 * \param[in] name  name of the parameter
 *
 * \return  a const char pointer containing the description.
 */
/*----------------------------------------------------------------------------*/

const char *
cs_notebook_parameter_get_description(char *name);

/*----------------------------------------------------------------------------*/
/*
 * \brief Get id associated with a notebook parameter.
 *
 * \param[in]   name      name of the parameter
 *
 * \return -1 if not present, id if present
 */
/*----------------------------------------------------------------------------*/

int
cs_notebook_parameter_get_id(const char  *name);

/*----------------------------------------------------------------------------*/
/*
 * \brief Get a group of notebook variable values
 *
 * \param[in]   n       number of notebook variables to query
 * \param[in]   ids     ids of notebook variables to query
 *                      (value set to 0 where id < 0)
 * \param[out]  values  values of notebook variables to query
 */
/*----------------------------------------------------------------------------*/

void
cs_notebook_get_values(int        n,
                       const int  ids[],
                       double     values[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Set a group of notebook variable values
 *
 * \param[in]  n       number of notebook variables to set
 * \param[in]  ids     ids of notebook variables to set
 *                     (ignored where id < 0)
 * \param[in]  values  values of notebook variables to set
 */
/*----------------------------------------------------------------------------*/

void
cs_notebook_set_values(int           n,
                       const int     ids[],
                       const double  values[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Destroy the notebook structure.
 *
 * Destroys the structures related to the notebook.
 */
/*----------------------------------------------------------------------------*/

void
cs_notebook_destroy_all(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Output the notebook info to the setup log.
 */
/*----------------------------------------------------------------------------*/

void
cs_notebook_log_setup(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Number of notebook variables
 *
 * \returns number of notebook variables (int)
 */
/*----------------------------------------------------------------------------*/

int
cs_notebook_nb_var(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Indicate if the notebook parameter is editable
 *
 * Returns a boolean to indicate wheter this parameter is editable
 *
 * \param[in]   id   Id of the notebook parameter
 *
 * \returns true is variable can be edited, false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_notebook_var_is_editable(int  id);

/*----------------------------------------------------------------------------*/
/*
 * \brief Indicate if the notebook parameter is read at restart
 *
 * Returns a boolean to indicate wheter this parameter is read at restart
 *
 * \param[in]   id   Id of the notebook parameter
 *
 * \returns true if variable should be read from checkpoint file, false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_notebook_var_is_read_from_checkpoint(int  id);

/*----------------------------------------------------------------------------*/
/*
 * \brief Change the editable property of the notebook parameter
 *
 * \param[in]   id   Id of the notebook parameter
 * \param[in]   val  flag (bool) indicating if the value is set to editable
 */
/*----------------------------------------------------------------------------*/

void
cs_notebook_var_change_editable(int   id,
                                bool  val);

/*----------------------------------------------------------------------------*/
/*
 * \brief Get name of a notebook parameter based on its id
 *
 * \param[in]   id    Id of the notebook parameter
 * \param[out]  name  Name of the notebook parameter
 *
 * \returns name of variable (char *)
 */
/*----------------------------------------------------------------------------*/

const char *
cs_notebook_name_by_id(int  id);

/*----------------------------------------------------------------------------*/
/*
 * \brief Write uncertain values to output file.
 *
 * If input and output uncertain variables are provided, output values
 * are written to an output file: cs_uncertain_output.dat
 * Results are ordered in the definition order in the notebook.
 */
/*----------------------------------------------------------------------------*/

void
cs_notebook_uncertain_output(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_NOTEBOOK_H__ */
