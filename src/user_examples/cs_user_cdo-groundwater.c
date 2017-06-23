/*============================================================================
 * Set main parameters for the current simulation when the CDO kernel is used
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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
#include "cs_property.h"
#include "cs_advection_field.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_cdo-groundwater.c
 *
 * \brief Main user function for setting of a calculation with CDO for the
 *        groundwater flow module
 */
/*----------------------------------------------------------------------------*/

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
        cs_lnum_t           n_pts,
        const cs_real_t    *xyz,
        cs_real_t          *retval)
{
  /* Physical parameters */
  const double  ks = 1.15741e-4;
  const double  theta_r = 0.15, theta_s = 0.45, dtheta = theta_s - theta_r;
  const double  hr = -100;
  const double  td = -5*L*L*dtheta/(6*hr*ks);

  /* Time-dependent part */
  const double  alpha = 6 - 5*time/td;

  for (cs_lnum_t p = 0; p < n_pts; p++) {

    /* Space-dependent part */
    const double  xll = (xyz[3*p] - L)/L, beta = xll*xll;

    retval[p] = hr*(1 - beta/alpha);

  }
}

/* Same as get_sol but optimize for time=0 */
static void
get_ic(cs_real_t           time,
       cs_lnum_t           n_pts,
       const cs_real_t    *xyz,
       cs_real_t          *retval)
{
  CS_UNUSED(time);

  const double  hr = -100;

  for (cs_lnum_t p = 0; p < n_pts; p++) {

    const double  x = xyz[3*p], xll = (x - L)/L;

    retval[p] = 1-one6*xll*xll;
    retval[p] *= hr;

  }
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

int
cs_user_cdo_activated(void)
{
  /* CS_CDO_OFF     = -1 --> CDO schemes are not used (no activation)
     CS_CDO_WITH_FV =  0 --> CDO schemes are used as well as finite volume
     CS_CDO_ONLY    =  1 --> CDO schemes are exclusively used */

  return  CS_CDO_ONLY;
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
  /* ======================
     Boundary of the domain
     ====================== */

  /* Choose a boundary by default.
     Valid choice is CS_PARAM_BOUNDARY_WALL or CS_PARAM_BOUNDARY_SYMMETRY */
  cs_domain_set_default_boundary(domain, CS_PARAM_BOUNDARY_SYMMETRY);

  /* Add a new boundary
     >> cs_domain_add_boundary(domain,
                               type_of_boundary,
                               mesh_location_name);

     * mesh_location_name is either a predefined mesh location or one defined
     by the user
     * type_of_boundary is one of the following keyword
        CS_PARAM_BOUNDARY_WALL,
        CS_PARAM_BOUNDARY_INLET,
        CS_PARAM_BOUNDARY_OUTLET,
        CS_PARAM_BOUNDARY_SYMMETRY
  */

  cs_domain_add_boundary(domain, CS_PARAM_BOUNDARY_INLET, "left");
  cs_domain_add_boundary(domain, CS_PARAM_BOUNDARY_OUTLET, "right");

  /* =========================
     Generic output management
     ========================= */

  cs_domain_set_output_param(domain,
                             10,     // output log frequency
                             2);     // verbosity (-1: no, 0, ...)

  /* ====================
     Time step management
     ====================

     If there is an inconsistency between the max. number of iteration in
     time and the final physical time, the first condition encountered stops
     the calculation.
  */

  cs_domain_set_time_param(domain,
                           200,       // nt_max or -1 (automatic)
                           86400.);   // t_max or < 0. (automatic)

  /* Define the value of the time step
     >> cs_domain_def_time_step_by_value(domain, dt_val);
     >> cs_domain_def_time_step_by_func(domain, dt_func);

     The second way to define the time step enable complex definitions.
     dt_func must have the following prototype:

     double dt_func(int  nt_cur, double  time_cur)
  */

  cs_domain_def_time_step_by_value(domain, 4320);

  /* Rk: Final time is 10 days = 864000 and dt = 0.05 day i.e 20 iters
     for one day */

  /* ================================
     Activate groundwater flow module
     ================================

     For the groundwater flow module:
     >> cs_domain_activate_gwf(domain,
                               permeability_type,
                               Richards_time,
                               n_soils,
                               n_tracers);

     * permeability_type is one of the following keywords:
       "isotropic", "orthotropic" or "anisotropic"
     * Richards_time is one of the following keywords:
       "steady" or "unsteady"
     * n_soils should be at least equal to 1.

     * Consequences of the activation of the groundwater flow module are:
     - add a new equation named "Richards" along with an associated field named
       "hydraulic_head". Default boundary condition is set to "zero_flux".
     - define a new advection field named "darcian_flux"
     - define a new property called "permeability".
     - define a new property called "soil_capacity" if "unsteady" is chosen
  */

  cs_gwf_t  *gwf = cs_domain_activate_gwf(domain,
                                          "isotropic", // type of permeability
                                          "unsteady",  // steady or unsteady
                                          1,           // number of soils
                                          0);          // number of tracers

  /* =========
     Add soils
     ========= */

  cs_gwf_add_iso_soil_by_value(gwf,
                               CS_GWF_HYDRAULIC_TRACY, // Hydraulic model to use
                               "cells",                // mesh location
                               1.15741e-4,             // saturated permeability
                               0.45,                   // saturated moisture
                               1.0);                   // bulk density (useless)

  /* Set additional parameters defining this soil
     >> cs_gwf_set_soil_param(gw, mesh_location_name, key, keyval);

     If mesh_location_name is set to NULL, all soils are set.

     Available keys are:
     CS_SOILKEY_SAT_MOISTURE,  // Set the saturated moisture content
     CS_SOILKEY_RES_MOISTURE,  // Set the residual moisture content

     Keys specific to the Tracy model
     CS_SOILKEY_TRACY_SAT_H,   // Head related to the saturated moisture content
     CS_SOILKEY_TRACY_RES_H,   // Head related to the residual moisture content
  */

  cs_gwf_set_soil_param(gwf, "cells", CS_SOILKEY_TRACY_RES_H, -100);
  cs_gwf_set_soil_param(gwf, NULL, CS_SOILKEY_RES_MOISTURE, 0.15);

  /* ====================
     Add tracer equations
     ====================

     Add a tracer equation which is unsteady and convected by the darcean flux
     >> cs_domain_add_gwf_tracer_eq(domain, eqname, varname);

     This implies the creation of a new equation called eqname and a new
     field called varname.
  */

  /* Set parameters related to each tracer equation in each soil
     >> cs_domain_set_gwf_tracer_eq(domain,
                                    eqname,
                                    mesh_location_name,
                                    water_diff,
                                    alpha_l,
                                    alpha_t,
                                    kd,
                                    lambda);

     According to the setting, additional properties can be created which are
     associated to the diffusion and/or reaction terms.

     If mesh_location_name is set to NULL, all soils are set.
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
  /* Retrieve the equation to set
     cs_equation_t  *eq = cs_domain_get_equation(domain, "eq_name");
  */

  cs_equation_t  *eq = NULL;

  /* =================
     Richards equation
     ================= */

  eq = cs_domain_get_equation(domain, "Richards");

  /* Define the boundary conditions
     >> cs_equation_add_bc_by_analytic(eq,
                                       bc_type,
                                       "mesh_location_name",
                                       analytic_function);

     -> eq is the structure related to the equation to set
     -> type of boundary condition:
        CS_PARAM_BC_DIRICHLET, CS_PARAM_BC_HMG_DIRICHLET,
        CS_PARAM_BC_NEUMANN, CS_PARAM_BC_HMG_NEUMANN, CS_PARAM_BC_ROBIN

     >> cs_equation_add_bc_by_value(eq,
                                    bc_type,
                                    "mesh_location_name",
                                    get);

     -> get : accessor to the value
  */

  cs_equation_add_bc_by_analytic(eq,
                                 CS_PARAM_BC_DIRICHLET,
                                 "left",    // name of the mesh location
                                 get_sol);  // analytic function

  /* Value to set */
  cs_get_t  get_bc = {.val = -100};
  cs_equation_add_bc_by_value(eq,
                              CS_PARAM_BC_DIRICHLET,
                              "right",  // name of the related mesh location
                              get_bc);     // value to set

  /* Define the initial condition by an analytical function
     (By default: zero is set) */
  cs_equation_set_ic_by_analytic(eq,       // equation
                                 NULL,     // NULL --> all cells
                                 get_ic);  // pointer to the analytic function
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
