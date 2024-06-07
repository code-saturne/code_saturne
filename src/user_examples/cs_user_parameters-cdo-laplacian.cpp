/*============================================================================
 * User functions for input of calculation parameters.
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
#include <math.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_parameters-cdo-laplacian.c
 *
 * \brief User functions for input of calculation parameters.
 *
 * See \ref parameters for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select physical model options, including user fields.
 *
 * This function is called at the earliest stages of the data setup,
 * so field ids are not available yet.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_model(void)
{
  /*! [param_cdo_laplacian_init] */
  {
    /* Activate CDO/HHO mode */

    cs_param_cdo_mode_set(CS_PARAM_CDO_MODE_ONLY);

    /* Add a user-defined equation */

    cs_equation_add_user("Laplacian", /* equation name */
                         "potential", /* associated variable field name */
                         1,           /* dimension of the unknown */
                      CS_BC_SYMMETRY); /* default boundary condition */
  }
  /*! [param_cdo_laplacian_init] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify output user parameters.
 *
 * For CDO schemes, this function concludes the setup of properties,
 * equations, source terms...
 *
 * \param[in, out] domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_finalize_setup(cs_domain_t   *domain)
{
  CS_NO_WARN_IF_UNUSED(domain);

  /*! [param_cdo_laplacian_finalize] */
  {
    cs_equation_param_t  *eqp = cs_equation_param_by_name("Laplacian");

    /* The property named "unity" is defined by default. Associate this
       property to the diffusion term and add it to the equation settings. */

    cs_equation_add_diffusion(eqp, cs_property_by_name("unity"));

    /* Boundary conditions (One assumes that two boundary zones named "X0" and
       "X1" exist. This is the case for instance if a Cartesian mesh is
       generated from the GUI. The two boundary zones are added in the
       GUI. Label is set to "X0" (resp. "X1") and selection criterion "X0"
       (resp. "X1").
    */

    cs_real_t  T0 = 0, T1 = 1;
    cs_equation_add_bc_by_value(eqp, CS_BC_DIRICHLET, "X0", &T0);
    cs_equation_add_bc_by_value(eqp, CS_BC_DIRICHLET, "X1", &T1);
  }
  /*! [param_cdo_laplacian_finalize] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
