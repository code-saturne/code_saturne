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
#include <math.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*
 * \file cs_user_boundary_conditions.c
 *
 * \brief User functions for boundary condition definitions.
 */
/*----------------------------------------------------------------------------*/

/*=============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
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
  CS_UNUSED(domain);

  /* Air inlet with pulverized coal, e.g. secondary or tertiary air. */

  /*! [example_1] */
  {
    const cs_zone_t *z = cs_boundary_zone_by_name("inlet");

    cs_coal_bc_inlet_t *ci = cs_coal_boundary_conditions_get_inlet(z);

    /* Inlet velocity direction based on constant vector;
       value based on oxidizer mass flow rate (in kg/s);*/

    cs_real_t u[] = {0, 0, 0.5};
    cs_boundary_conditions_open_set_velocity_by_value(z, u);

    cs_boundary_conditions_open_set_mass_flow_rate_by_value(z, 1.46e-03);

    ci->ientat = 1;

    /* Oxidizer's number (1 to 3) */
    ci->inmoxy = 1;

    /* Oxidizer Temperature in K */
    ci->t_air = 400.0 + cs_physical_constants_celsius_to_kelvin;

    /* Turbulence variables BC's calculated based on hydraulic diameter */

    cs_boundary_conditions_inlet_set_turbulence_hyd_diam(z, 0.032);

    /* Automatic treatment of non-user scalars */

    /* Treatment of user-defined scalars */

    const char *names[] = {"scalar1", "scalar2"};

    for (int i = 0; i < 2; i++) {

      cs_equation_param_t  *eqp = cs_equation_param_by_name(names[i]);
      if (eqp == NULL)
        continue;

      cs_real_t val[1] = {1.0};
      cs_equation_add_bc_by_value(eqp,
                                  CS_BC_DIRICHLET,
                                  z->name,
                                  val);

    }
  }
  /*! [example_1] */

  /* Primary air inlet with pulverized coal. */

  /*! [example_2] */
  {
    const cs_zone_t *z = cs_boundary_zone_by_name("primary_inlet");

    cs_coal_bc_inlet_t *ci = cs_coal_boundary_conditions_get_inlet(z);

    /* Inlet velocity direction based on constant vector;
       value based on oxidizer mass flow rate (in kg/s);*/

    cs_real_t u[] = {0, 0, 0.5};
    cs_boundary_conditions_open_set_velocity_by_value(z, u);

    cs_coal_boundary_conditions_inlet_set_air_mass_flow_rate_by_value
      (z, 1.46e-03);

    /* Turbulence variables BC's calculated based on hydraulic diameter */

    cs_boundary_conditions_inlet_set_turbulence_hyd_diam(z, 0.1);

    /* Air inlet with pulverized coal */
    ci->ientcp = 1;

    /* Oxidizer's number (1 to 3) */
    ci->inmoxy = 1;

    /* Oxidizer Temperature in K */
    ci->t_air = 800.0 + cs_physical_constants_celsius_to_kelvin;

    /* code_saturne deals with ncha different coals (component of blend)
       every coal is described by n_classes_per_coal[icha] class of particles
       (each of them described by an inlet diameter) */

    cs_coal_model_t *coal = cs_glob_coal_model;

    /* Treatment for the first coal */

    int icha = 0;

    /* Coal mass flow rate in kg/s */
    ci->qimpcp[icha] = 1.46e-4;

    /* Percentage mass fraction of each granulometric class */
    for (int iclapc = 0; iclapc < coal->n_classes_per_coal[icha]; iclapc++) {
      ci->distch[icha][iclapc] = 100.0/(double)coal->n_classes_per_coal[icha];
    }

    /* Inlet temperature for coal & primary air */
    ci->timpcp[icha] = 800.0 + cs_physical_constants_celsius_to_kelvin;

    /* Automatic treatment of non-user scalars */
  }
  /*! [example_2] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
