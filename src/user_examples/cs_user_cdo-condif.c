/*============================================================================
 * Set main parameters for the current simulation when the CDO kernel is used
 *============================================================================*/

/* VERS */

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
 * \file cs_user_cdo-condif.c
 *
 * \brief Main user function for setting of a calculation with CDO.
 *
 */
/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*! \endcond (end ignore by Doxygen) */

static const double  one_third = 1./3.;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/* ---------------------------------------------------------------------------
 * TEST 1 -- Boundary conditions
 * -------------------------------------------------------------------------- */

static void
_define_adv_field(cs_real_t           time,
                  cs_lnum_t           n_pts,
                  const cs_real_t    *xyz,
                  cs_real_t          *res)
{
  CS_UNUSED(time);

  for (cs_lnum_t p = 0; p < n_pts; p++) {

    const cs_real_t  *pxyz = xyz + 3*p;
    cs_real_t  *pres = res + 3*p;

    pres[0] = pxyz[1] - 0.5;
    pres[1] = 0.5 - pxyz[0];
    pres[2] = pxyz[2];

  }
}

/* ---------------------------------------------------------------------------
 * TEST 1 -- Boundary conditions
 * -------------------------------------------------------------------------- */

static void
_define_bcs(cs_real_t           time,
            cs_lnum_t           n_pts,
            const cs_real_t    *xyz,
            cs_real_t          *res)
{
  CS_UNUSED(time);

  const double  pi = 4.0*atan(1.0);
  for (cs_lnum_t p = 0; p < n_pts; p++) {

    const cs_real_t  *pxyz = xyz + 3*p;

    res[p] = 1 +
      sin(pi*pxyz[0]) * sin(pi*(pxyz[1]+0.5)) * sin(pi*(pxyz[2]+one_third));

  }
}

/* ---------------------------------------------------------------------------
 * TEST 1 -- Source term
 * -------------------------------------------------------------------------- */

static void
_define_source(cs_real_t           time,
               cs_lnum_t           n_pts,
               const cs_real_t    *xyz,
               cs_real_t          *res)
{
  CS_UNUSED(time);
  const double  pi = 4.0*atan(1.0), pi2 = pi*pi;

  for (cs_lnum_t p = 0; p < n_pts; p++) {

    const double  x = xyz[3*p], y = xyz[3*p+1], z = xyz[3*p+2];
    const double  cpx = cos(pi*x), spx = sin(pi*x);
    const double  cpy = cos(pi*(y+0.5)), spy = sin(pi*(y+0.5));
    const double  cpz = cos(pi*(z+one_third)), spz = sin(pi*(z+one_third));

    /* first derivatives */
    cs_real_t gx = pi*cpx*spy*spz;
    cs_real_t gy = pi*spx*cpy*spz;
    cs_real_t gz = pi*spx*spy*cpz;

    /* second derivatives */
    cs_real_t gxx, gyy, gzz, gxy, gxz, gyz;
    gxx = gyy = gzz = -pi2*spx*spy*spz;
    gxy = pi2*cpx*cpy*spz, gxz = pi2*cpx*spy*cpz, gyz = pi2*spx*cpy*cpz;

    /* Material property */
    cs_real_33_t  cond;
    cond[0][0] = 1.0, cond[0][1] = 0.5, cond[0][2] = 0.0;
    cond[1][0] = 0.5, cond[1][1] = 1.0, cond[1][2] = 0.5;
    cond[2][0] = 0.0, cond[2][1] = 0.5, cond[2][2] = 1.0;

    /* Contribution of the diffusive part */
    res[p] = cond[0][0]*gxx + cond[1][1]*gyy + cond[2][2]*gzz +
      2*( cond[0][1]*gxy + cond[0][2]*gxz + cond[1][2]*gyz);
    res[p] *= -1;

    /* Contribution of the advection term */
    res[p] += (y - 0.5)*gx + (0.5 - x)*gy + z*gz + 1 + spx*spy*spz;

  } // Loop on evaluation points
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

  cs_mesh_location_add("in", CS_MESH_LOCATION_BOUNDARY_FACES, "x < 1e-5");
  cs_mesh_location_add("out", CS_MESH_LOCATION_BOUNDARY_FACES, "x > 0.9999");
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
  cs_domain_set_default_boundary(domain, CS_PARAM_BOUNDARY_WALL);

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

  cs_domain_add_boundary(domain, CS_PARAM_BOUNDARY_INLET, "in");
  cs_domain_add_boundary(domain, CS_PARAM_BOUNDARY_OUTLET, "out");

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
                           100,     // nt_max or -1 (automatic)
                           -1.);    // t_max or < 0. (automatic)

  /* Define the value of the time step
     >> cs_domain_def_time_step_by_value(domain, dt_val);
     >> cs_domain_def_time_step_by_func(domain, dt_func);

     The second way to define the time step enable complex definitions.
     dt_func must have the following prototype:

     double dt_func(int  nt_cur, double  time_cur)
  */

  cs_domain_def_time_step_by_value(domain, 1.0);

  /* =============================
     Activate predefined equations
     =============================

     For the wall distance:
       cs_domain_activate_wall_distance(domain);

  */

  cs_domain_activate_wall_distance(domain);

  /* =========================================
     Define additional user equations to solve
     =========================================

     cs_domain_add_user_equation(...)

     >> available type of equation: "scalar", "vector" or "tensor"
     >> default_bc: "zero_value" or "zero_flux"

     By default, initial values are set to zero (or the value given by the
     restart file in case of restart).
  */

  cs_domain_add_user_equation(domain,
                              "AdvDiff",     // equation name
                              "Potential",   // associated variable field name
                              "scalar",      // type of equation
                              "zero_value"); // default boundary condition

  /* ================================
     User-defined material properties
     ================================

     By default, one material property is defined:
     >> "unity" (isotropic and value equal 1.0)

     Users can also define additional material properties
     cs_domain_add_property(domain,
                            "name_of_property",
                            "type_keyword",
                            n_subdomains);

      type_keyword has predefined values among:
        >> "isotropic", "orthotropic" or "anisotropic"

      n_subdomains corresponds to the number of subdomains used to define
      this property (if 1 is given, use "cells" as mesh location, otherwise
      give the name of a mesh locations based on a selection of cells which
      has been previously defined).
  */

  cs_domain_add_property(domain,
                         "conductivity",  // property name
                         "anisotropic",   // type of material property
                         1);              // definition in n_subdomains

  cs_domain_add_property(domain,
                         "rho.cp",       // property name
                         "isotropic",    // type of material property
                         1);             // definition in n_subdomains

  /* =============================
     User-defined advection fields
     =============================

     Users can also define advection fields
     cs_domain_add_advection_field(domain,
                                   "name_of_advection_field");

  */

  cs_domain_add_advection_field(domain,
                                "adv_field");
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
  /* =======================
     User-defined properties
     =======================

     Retrieve the property to set
     cs_property_t  *pty = cs_domain_get_property(domain, "pty_name");

    Several ways exist to define a property
      >> cs_property_def_by_value(pty, "location_name", value);
      >> cs_property_def_by_analytic(pty, "location_name", func);
      >> cs_property_def_by_law(pty, "location_name", func);

      -- pty is the structure related to the property to set
      -- "location_name" is the name of the mesh location (based on cells)
      where the given definition has to be applied

      For a definition by value:
      -- value is "1.0" for instance for an isotropic property
         or "0.5 0.1 1." for instance for an orthotropic property
      For a definition by analytical function or law
      -- func is a function with a predefined prototype

  */

  cs_property_t  *conductivity = cs_domain_get_property(domain, "conductivity");

  cs_property_def_by_value(conductivity,     // property structure
                           "cells",          // name of the mesh location
                           "1.0  0.5  0.0\n" // values of the property
                           "0.5  1.0  0.5\n"
                           "0.0  0.5  1.0\n");

  cs_property_t  *rhocp = cs_domain_get_property(domain, "rho.cp");

  cs_property_def_by_value(rhocp,    // property structure
                           "cells",  // name of the mesh location
                           "1.0");   // value of the property

  /* =============================
     User-defined advection fields
     =============================

     Retrieve the advection field to set
     cs_adv_field_t  *adv = cs_domain_get_advection_field(domain, "adv_name");

     Several ways exist to define an advection field
      >> cs_advection_field_def_by_value(adv, values);
         -- adv is the structure related to the advection field to set
         -- values is "0.5 0.1 1." for instance

      >> cs_property_def_by_analytic(adv, func);
         -- adv is the structure related to the advection field to set
         -- func is a function with a predefined prototype

  */

  cs_adv_field_t  *adv = cs_domain_get_advection_field(domain, "adv_field");

  cs_advection_field_def_by_analytic(adv, _define_adv_field);

  /* Enable also the defintion of the advection field at mesh vertices */
  cs_advection_field_set_option(adv, CS_ADVKEY_DEFINE_AT, "vertices");

  /* Activate the post-processing of the advection field */
  cs_advection_field_set_option(adv, CS_ADVKEY_POST, "field");

  /* Activate the post-processing of the related Courant number */
  cs_advection_field_set_option(adv, CS_ADVKEY_POST, "courant");

  /* ======================
     User-defined equations
     ====================== */

  /* Retrieve the equation to set
     >> cs_equation_t  *eq = cs_domain_get_equation(domain, "eq_name");
  */

  cs_equation_t  *eq = cs_domain_get_equation(domain, "AdvDiff");

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
                                 "boundary_faces",  // name of the mesh location
                                 _define_bcs);      // pointer to the function

  /* Link properties to different terms of this equation
     >> cs_equation_link(eq,
                         "term_keyword",
                         structure_to_link);

     - eq is the structure related to the equation to set
     - Keyword related to the term to set is a choice among:
       >> "diffusion", "time" or "advection"
     - If keyword is "time" or "diffusion", the structure to link is a
       property.
       If keyword is "advection", the structure to link is an advection field
  */

  /* Activate unsteady effect */
  cs_equation_link(eq, "time", rhocp);
  /* Activate diffusion effect */
  cs_equation_link(eq, "diffusion", conductivity);
  /* Activate advection effect */
  cs_equation_link(eq, "advection", adv);

  /* Add a source term:
     >> cs_equation_add_source_term_by_val(eq,
                                           label,
                                           ml_name,
                                           val);
     or
     >> cs_equation_add_source_term_by_analytic(eq,
                                                label,
                                                ml_name,
                                                analytic_func);

     where label of the source term is optional (i.e. NULL is possible)
     This label is mandatory if additional settings are requested only for
     this specific source term.

     where ml_name is the name of the mesh location where this source term is
     defined (a selection of cells)

     where val is the value of the source term by m^3
     or where analytic_func is the name of the analytical function
   */

  cs_source_term_t  *st
    = cs_equation_add_source_term_by_analytic(eq,
                                              "SourceTerm",    // label
                                              "cells",         // ml_name
                                              _define_source); // function

  /* Optional: specify the quadrature used for computing a source term

     >> key takes value among
     CS_QUADRATURE_BARY     used the barycenter approximation
     CS_QUADRATURE_HIGHER   used 4 Gauss points for approximating the integral
     CS_QUADRATURE_HIGHEST  used 5 Gauss points for approximating the integral
  */

  cs_source_term_set_quadrature(st, CS_QUADRATURE_BARY);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
