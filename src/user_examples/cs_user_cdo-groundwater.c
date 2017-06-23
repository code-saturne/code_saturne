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
#include "cs_boundary_zone.h"
#include "cs_volume_zone.h"
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
 * \brief  Specify additional zones and/or mesh locations
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

  cs_boundary_zone_define("left", "x < 1e-3", 0);

  char cmd[20];
  sprintf(cmd, "x > %10.7e", L - 1e-5);
  cs_boundary_zone_define("right", cmd, 0);

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
cs_user_cdo_init_setup(cs_domain_t   *domain)
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
     >> cs_gwf_activate(permeability_type,
                        richards_flag);

     * permeability_type is one of the following keywords:
       CS_PROPERTY_ISO, CS_POPERTY_ORTHO or CS_PROPERTY_ANISO

     * richards_flag are
     CS_GWF_GRAVITATION, CS_GWF_RICHARDS_UNSTEADY, CS_GWF_SOIL_PROPERTY_UNSTEADY
     CS_GWF_SOIL_ALL_SATURATED

     * Consequences of the activation of the groundwater flow module are:
     - add a new equation named "Richards" along with an associated field named
       "hydraulic_head". The default boundary condition set is a homogeneous
       Neumann.
     - define a new advection field named "darcian_flux"
     - define a new property called "permeability".
     - define a new property called "soil_capacity" if "unsteady" is chosen
  */

  cs_gwf_activate(CS_PROPERTY_ISO,
                  CS_GWF_RICHARDS_UNSTEADY);

  /* =========
     Add soils
     ========= */

  cs_gwf_soil_add("cells", CS_GWF_SOIL_USER);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Specify for each soil and tracer how is defined each term of the
 *         the tracer equation. Soils and tracer equations have to be added
 *         previously
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
*/
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_setup_gwf(cs_domain_t   *domain)
{
  CS_UNUSED(domain);

  /* =========
     Set soils
     ========= */

  cs_gwf_soil_t  *soil = cs_gwf_soil_by_name("cells");
  cs_gwf_set_iso_saturated_soil(soil,
                                1.15741e-4,  // saturated permeability
                                0.45,        // saturated moisture
                                1.0);        // bulk density (useless)

  /* ===========
     Set tracers
     =========== */

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
cs_user_cdo_finalize_setup(cs_domain_t   *domain)
{
  CS_UNUSED(domain);

  /* =================
     Richards equation
     ================= */

  cs_equation_t  *eq = cs_equation_by_name("Richards");

  /* Define the boundary conditions
     -> type of boundary condition:
        CS_PARAM_BC_DIRICHLET, CS_PARAM_BC_HMG_DIRICHLET,
        CS_PARAM_BC_NEUMANN, CS_PARAM_BC_HMG_NEUMANN, CS_PARAM_BC_ROBIN

     >> cs_equation_add_bc_by_value(eq,
                                    bc_type,
                                    "mesh_location_name",
                                    val);  // pointer to cs_real_t  */

  cs_equation_add_bc_by_analytic(eq,
                                 CS_PARAM_BC_DIRICHLET,
                                 "left",            // name of the boundary zone
                                 (void *)get_sol);  // analytic function

  /* Value to set */
  cs_real_t  bc_val = -100;
  cs_equation_add_bc_by_value(eq,
                              CS_PARAM_BC_DIRICHLET,
                              "right",  // name of the related mesh location
                              &bc_val); // value to set

  /* Define the initial condition by an analytical function
     (By default: zero is set) */
  cs_equation_add_ic_by_analytic(eq,       // equation
                                 NULL,     // NULL --> all cells
                                 (void *)get_ic);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the bulk density related to a soil structure
 *
 * \param[in]  soil      pointer to a cs_gwf_soil_t structure
 * \param[out] density   return value for the density
 */
/*----------------------------------------------------------------------------*/

void
cs_user_gwf_get_soil_density(const cs_gwf_soil_t   *soil,
                             cs_real_t             *density)
{
  CS_UNUSED(soil);
  *density = 1.0;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
