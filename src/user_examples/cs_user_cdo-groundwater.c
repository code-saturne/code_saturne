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
#include "cs_sdm.h"
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

/* Permeability in each subdomain */
static const double k1 = 1e5, k2 = 1;

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the explicit definition of the problem for the Richards eq.
 *         pt_ids is optional. If not NULL, it enables to access in coords
 *         at the right location and the same thing to fill retval if compact
 *         is set to false
 *         Rely on a generic function pointer for an analytic function
 *
 * \param[in]      time      when ?
 * \param[in]      n_elts    number of elements to consider
 * \param[in]      pt_ids    list of elements ids (to access coords and fill)
 * \param[in]      coords    where ?
 * \param[in]      compact   true:no indirection, false:indirection for filling
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 * \param[in, out] retval    result of the function
 */
/*----------------------------------------------------------------------------*/

static inline void
get_tracer_sol(cs_real_t          time,
               cs_lnum_t          n_points,
               const cs_lnum_t   *pt_ids,
               const cs_real_t   *xyz,
               bool               compact,
               void              *input,
               cs_real_t         *retval)
{
  CS_UNUSED(input);

  /* Physical parameters */
  const double  magnitude = 2*k1/(k1 + k2);
  const double  x_front = magnitude * time;

  if (pt_ids != NULL && !compact) {

    for (cs_lnum_t  i = 0; i < n_points; i++) {

      const cs_lnum_t  id = pt_ids[i];
      const double  x = xyz[3*id];
      if (x <= x_front)
        retval[id] = 1;
      else
        retval[id] = 0;
    }

  }
  else if (pt_ids != NULL && compact) {

    for (cs_lnum_t  i = 0; i < n_points; i++) {
      const double  x = xyz[3*pt_ids[i]];
      if (x <= x_front)
        retval[i] = 1;
      else
        retval[i] = 0;
    }

  }
  else {

    for (cs_lnum_t  i = 0; i < n_points; i++) {
      const double  x = xyz[3*i];
      if (x <= x_front)
        retval[i] = 1;
      else
        retval[i] = 0;
    }

  }

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate or not the CDO module
 */
/*----------------------------------------------------------------------------*/

int
cs_user_cdo_activated(void)
{
  /* CS_PARAM_CDO_MODE_OFF     = -1 --> CDO schemes are not used (no activation)
     CS_PARAM_CDO_MODE_WITH_FV =  0 --> Both CDO schemes and FV activated
     CS_PARAM_CDO_MODE_ONLY    =  1 --> CDO schemes are exclusively used */

  return  CS_PARAM_CDO_MODE_ONLY;
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
                             100,     // output log frequency
                             3);     // verbosity (-1: no, 0, ...)

  /* ====================
     Time step management
     ====================

     If there is an inconsistency between the max. number of iteration in
     time and the final physical time, the first condition encountered stops
     the calculation.
  */

  cs_domain_set_time_param(domain,
                           500,    // nt_max or -1 (automatic)
                           -1.);   // t_max or < 0. (automatic)

  /* Define the value of the time step
     >> cs_domain_def_time_step_by_value(domain, dt_val);
     >> cs_domain_def_time_step_by_func(domain, dt_func);

     The second way to define the time step enable complex definitions.
     dt_func must have the following prototype:

     double dt_func(int  nt_cur, double  time_cur)
  */

  cs_domain_def_time_step_by_value(domain, 1e-3);

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
     or 0 if there is no flag to set

     * Consequences of the activation of the groundwater flow module are:
     - add a new equation named "Richards" along with an associated field named
       "hydraulic_head". The default boundary condition set is a homogeneous
       Neumann.
     - define a new advection field named "darcian_flux"
     - define a new property called "permeability".
     - define a new property called "soil_capacity" if "unsteady" is chosen
  */

  cs_gwf_activate(CS_PROPERTY_ISO, 0); // no flag to set

  /* =========
     Add soils (must be done before adding tracers)
     ========= */

  cs_gwf_soil_add("soil1", CS_GWF_SOIL_SATURATED);
  cs_gwf_soil_add("soil2", CS_GWF_SOIL_SATURATED);

  /* ====================
     Add tracer equations
     ====================

     Add a tracer equation which is unsteady and convected by the darcean flux
     This implies the creation of a new equation called eqname along with a new
     field called varname.

     For standard tracer:
       cs_gwf_add_tracer(eqname, varname);
     For user-defined tracer
       cs_gwf_add_tracer_user(eqname, varname, setup_func);
  */

  cs_gwf_add_tracer("Tracer_01","C1");

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

  cs_gwf_soil_t  *s1 = cs_gwf_soil_by_name("soil1");
  cs_gwf_set_iso_saturated_soil(s1,
                                k1,     // saturated permeability
                                1.0,    // saturated moisture
                                1.0);   // bulk density (useless)

  cs_gwf_soil_t  *s2 = cs_gwf_soil_by_name("soil2");
  cs_gwf_set_iso_saturated_soil(s2,
                                k2,     // saturated permeability
                                1.0,    // saturated moisture
                                1.0);   // bulk density (useless)

  /* ===========
     Set tracers (soil have to be defined first)
     ===========

     Set parameters related to each tracer equation in each soil */

  cs_gwf_tracer_t *tr = cs_gwf_tracer_by_name("Tracer_01");
  cs_gwf_set_standard_tracer(tr,
                             NULL,        // zone name or NULL for all
                             0.,          // water molecular diffusivity
                             0., 0.,      // alpha (longi. and transvesal)
                             0.,          // distribution coef.
                             0.);         // 1st order decay coef.
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

  /*=============
    Set equations
    ============= */

  cs_equation_t  *eq = NULL;

  /* Richards equation */
  eq = cs_equation_by_name("Richards");

  /* Define the boundary conditions  */
  cs_real_t  val0 = 0.0, val1 = 1.0;

  cs_equation_add_bc_by_value(eq,
                              CS_PARAM_BC_DIRICHLET,
                              "left", // boundary zone name
                              &val1);  // value to set

  cs_equation_add_bc_by_value(eq,
                              CS_PARAM_BC_DIRICHLET,
                              "right", // boundary zone name
                              &val0);  // value to set

  /* Tracer equation */
  eq = cs_equation_by_name("Tracer_01");

  cs_equation_add_bc_by_value(eq,
                              CS_PARAM_BC_DIRICHLET,
                              "left",  // boundary zone name
                              &val1);  // value to set

  /* Define the initial condition (By default: zero is set) */
  cs_equation_add_ic_by_analytic(eq, "cells", get_tracer_sol, NULL);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
