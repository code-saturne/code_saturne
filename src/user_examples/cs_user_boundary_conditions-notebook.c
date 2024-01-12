/*============================================================================
 * User definition of boundary conditions.
 *============================================================================*/

/* VERS */

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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_dof_func_t function to compute an inlet value based
 *        on the value of a notebook variable named "t_inlet" .
 *
 * For the calling function, elt_ids is optional. If not NULL, array(s) should
 * be accessed with an indirection. The same indirection can be applied to fill
 * retval if dense_output is set to false.
 * In the current case, retval is allocated to mesh->n_b_faces.
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         NULL or pointer to a structure cast on-the-fly
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

/*! [notebook_dof_func] */
static void
_notebook_t_inlet(cs_lnum_t         n_elts,
                  const cs_lnum_t  *elt_ids,
                  bool              dense_output,
                  void             *input,
                  cs_real_t        *retval)
{
  /* Get a user parameter defined in the GUI notebook */

  cs_real_t t_bnd = cs_notebook_parameter_value_by_name("t_inlet");

  /* Apply values at selected locations */

  for (cs_lnum_t i = 0; i < n_elts; i++) {
    const cs_lnum_t  face_id = (elt_ids == NULL) ? i : elt_ids[i];
    const cs_lnum_t  j = dense_output ? i : face_id;
    retval[j] = t_bnd;
  }
}
/*! [notebook_dof_func] */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set boundary conditions to be applied.
 *
 * This function is called just before \ref cs_user_finalize_setup, and
 * boundary conditions can be defined in either of those functions,
 * depending on whichever is considered more readable or practical for a
 * given use.
 *
 * \param[in, out]  domain  pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_boundary_conditions_setup(cs_domain_t  *domain)
{
  {
    /*! [notebook_dof_inlet] */
    cs_equation_param_t  *eqp = cs_equation_param_by_name("scalar_1");

    cs_equation_add_bc_by_dof_func(eqp,
                                   CS_PARAM_BC_DIRICHLET,
                                   "inlet",                // zone name
                                   cs_flag_boundary_face,  // location flag
                                   _notebook_t_inlet,      // callback function
                                   NULL);                  // input structure
    /*! [notebook_dof_inlet] */
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
