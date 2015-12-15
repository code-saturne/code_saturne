/*============================================================================
 * Set main parameters for the current simulation when the CDO kernel is used
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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

#include <errno.h>
#include <locale.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_mesh_location.h"
#include "cs_cdo_toolbox.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_user_cdo.c

  \brief  Set main parameters for the current simulation when the CDO kernel
          is used

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*! \endcond (end ignore by Doxygen) */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

static const double  one6 = 1/6.;
static const double  L = 200;

/* Solution of the TRACY 1D verification testcase
        F.T. Tracy, "1D, 2D, 3D analytical solutions of unsaturated flow in
        groundwater", Journal of Hydrology, 170, pp. 199--214 (1995)
*/
static void
get_sol(cs_real_t           time,
        const cs_real_3_t   xyz,
        cs_get_t           *get)
{
  /* Physical parameters */
  const double  ks = 1.15741e-4;
  const double  theta_r = 0.15, theta_s = 0.45, dtheta = theta_s - theta_r;
  const double  hr = -100;
  const double  td = -5*L*L*dtheta/(6*hr*ks);

  /* Space-dependent part */
  const double  xll = (xyz[0] - L)/L, beta = xll*xll;
  /* Time-dependent part */
  const double  alpha = 6 - 5*time/td;

  (*get).val = hr*(1 - beta/alpha);
}

/* Same as get_sol but optimize for time=0 */
static void
get_ic(cs_real_t           time,
       const cs_real_3_t   xyz,
       cs_get_t           *get)
{
  const double  x = xyz[0], xll = (x - L)/L;
  const double  hr = -100;

  (*get).val = 1-one6*xll*xll;
  (*get).val *= hr;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate or not the CDO module
 */
/*----------------------------------------------------------------------------*/

bool
cs_user_cdo_activated(void)
{
  return  true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Specify additional mesh locations
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_add_mesh_locations(void)
{
  /* ===========================
     Define mesh locations
     ===========================

     By default several mesh locations are predefined
     >> "cells"
     >> "interior_faces"
     >> "boundary_faces"
     >> "vertices"

 */

  cs_mesh_location_add("left", CS_MESH_LOCATION_BOUNDARY_FACES, "x < 1e-3");

  char cmd[20];
  const double  tol = 1e-5;

  sprintf(cmd, "x > %10.7e", L-tol);
  cs_mesh_location_add("right", CS_MESH_LOCATION_BOUNDARY_FACES, cmd);

  return;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Specify for the computational domain:
 *         -- which type of boundaries closed the computational domain
 *         -- the settings for the time step
 *         -- activate predefined equations or modules
 *         -- add user-defined properties and/or advection fields
 *         -- add user-defined equations
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_init_domain(cs_domain_t   *domain)
{

  /* =========================================
     Define boundary of the domain
     =========================================

     Choose a boundary by default
     boundary keyword is one of the following keyword
     >> wall or symmetry
  */

  cs_domain_set_default_boundary(domain, "wall");

  /* Add a boundary
     >>   cs_domain_add_boundary(domain,
                                 mesh location name,
                                 boundary keyword)

     mesh location name is either a predefined mesh location or one defined
     by user

     boundary keyword is one of the following keyword
     >> wall, inlet, outlet, symmetry
  */

  cs_domain_add_boundary(domain, "left", "inlet");
  cs_domain_add_boundary(domain, "right", "outlet");

  /* =========================================
     Time step management
     =========================================

     If there is an inconsistency between the max. number of iteration in
     time and the final physical time, the first condition encountered stops
     the calculation.

     Type of definition is among the following choices:
     >> "value", "time_func", "user"
     By value, the time step is constant

  */

  cs_domain_set_time_step(domain,
                          864000.,  /* Final time of the simulation */
                          200,      /* Max. number of time iteration */
                          "value",  /* How time step is define */
                          "4320");  /* Value of the time step (2160) */

  /* Rk: Final time is 10 days = 864000 and dt = 0.05 day i.e 20 iters
     for one day */

  /* ================================
     User-defined material properties
     ================================

     By default, one material property is defined:
     >> "unity" (isotropic and value equal 1.0)

     Users can also define additional material properties
     cs_domain_add_property(domain,
                            "name_of_property",
                            "type_keyword");

      type_keyword has predefined values among:
        >> "isotropic", "orthotropic" or "anisotropic"
  */

  /* =========================================
     Activate predefined equations
     =========================================

     For the wall distance:
       cs_domain_activate_wall_distance(domain);

  */

  /* =========================================
     Activate groundwater module
     =========================================

     For the groundwater module:
       cs_domain_activate_groundwater(domain, "model_keyword");

     where "model_keyword" is chosen among the following list:
     >> "saturated", "tracy", "genutchen"

     This call implies:
     - add an equation named "Richards" along with an associated field named
     "hydraulic_head". Default boundary condition is set to "zero_flux".
     - definition of an advection field named "darcian_flux"
     - definition of a property called "permeability".
     - definition a property called "soil_capacity" used to take into account
     the unsteady phenomena if the model is not "satured".

  */

  cs_domain_activate_groundwater(domain, "tracy");

  /* Set parameters associated to the groundwater module */
  cs_groundwater_t  *gw = cs_domain_get_groundwater(domain);

  cs_groundwater_set_param(gw, "saturated_permeability", "1.15741e-4");
  cs_groundwater_set_param(gw, "tracy_hr", "-100");
  cs_groundwater_set_param(gw, "max_moisture", "0.45");
  cs_groundwater_set_param(gw, "residual_moisture", "0.15");
  cs_groundwater_set_param(gw, "post_freq", "10");
  cs_groundwater_set_param(gw, "output_moisture", "true");

  /* Then add a tracer equation convected by the darcean velocity
     >> cs_groundwater_add_tracer_equation(domain,
                                           "eq_name",
                                           "variable_name",
                                           dispersivity vector,
                                           bulk density,
                                           distribution coefficient,
                                           first order rate reaction);
   */





}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  - Specify the elements such as properties, advection fields,
 *           user-defined equations and modules which have been previously
 *           added.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
*/
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_set_domain(cs_domain_t   *domain)
{
  /* ==========
     Properties
     ==========

     Retrieve the property to set
     cs_property_t  *pty = cs_domain_get_property(domain, "pty_name");

     Several ways exist to define a property
      >> cs_property_def_by_value(pty, value);
         -- pty is the structure related to the property to set
         -- value is "1.0" for instance for an isotropic property
            or "0.5 0.1 1." for instance for an orthotropic property

      >> cs_property_def_by_analytic(pty, func);
         -- pty is the structure related to the property to set
         -- func is a function with a predefined prototype

      >> cs_property_def_by_law(pty, func);
         -- pty is the structure related to the property to set
         -- func is a function with a predefined prototype

  */

  /* ===============
     Advection field
     ===============

     Retrieve the advection field to set
     cs_adv_field_t  *adv = cs_domain_get_advection_field(domain, "adv_name");

  */

  cs_adv_field_t  *adv = cs_domain_get_advection_field(domain, "darcian_flux");

  cs_advection_field_set_option(adv, "post_freq", "0");
  cs_advection_field_set_option(adv, "cell_field", "true");

  /* =========
     Equations
     =========

     Retrieve the equation to set
     cs_equation_t  *eq = cs_domain_get_equation(domain, "eq_name");

     Define the boundary conditions
     >> cs_equation_add_bc(eq,
                           "mesh_location_name",
                           "bc_type_keyword",
                           "definition_type_keyword",
                           pointer to the definition);

     -- eq is the structure related to the equation to set
     -- Keyword related to the boundary condition type is a choice among:
        >> "dirichlet", "neumann" or "robin"
     -- Keyword related to the type of definition is a choice among:
        >> "value", "analytic"

  */

  cs_equation_t  *eq = cs_domain_get_equation(domain, "Richards");

  /* Define the boundary conditions
     >> cs_equation_add_bc(eq,
                           "mesh_location_name",
                           "bc_type_keyword",
                           "definition_type_keyword",
                           pointer to the definition);

     -- eq is the structure related to the equation to set
     -- Keyword related to the boundary condition type is a choice among:
        >> "dirichlet", "neumann" or "robin"
     -- Keyword related to the type of definition is a choice among:
        >> "value", "analytic"

  */

  cs_equation_add_bc(eq,           // equation
                     "left",       // name of the mesh location
                     "dirichlet",  // BC type
                     "analytic",   // type of definition
                     get_sol);     // pointer to the analytic function

  cs_equation_add_bc(eq,           // equation
                     "right",      // name of the mesh location
                     "dirichlet",  // BC type
                     "value",      // type of definition
                     "-100");      // pointer to the analytic function

  /* Define the initial condition (By default: zero is set) */
  cs_equation_set_ic(eq,         // equation
                     "analytic", // type of definition
                     get_ic);    // pointer to the analytic function

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
