/*============================================================================
 * User functions for input of calculation parameters.
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

#include <assert.h>
#include <math.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_advection_field.h"
#include "cs_base.h"
#include "cs_domain.h"
#include "cs_field.h"
#include "cs_gwf.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_param.h"
#include "cs_property.h"
#include "cs_prototypes.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_parameters-cdo-gwf.c
 *
 * \brief User functions for setting a calcultion using the groundwater flow
 *        module with CDO schemes
 *
 * See \subpage parameters for examples.
 */
/*----------------------------------------------------------------------------*/

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
  /* Activate CDO/HHO module so that main additional structure are built */
  cs_domain_t  *domain = cs_glob_domain;

  cs_domain_set_cdo_mode(domain, CS_DOMAIN_CDO_MODE_ONLY);

  /* ======================
     Boundary of the domain
     ====================== */

  /* Choose a boundary by default */
  cs_domain_set_default_boundary(domain, CS_DOMAIN_BOUNDARY_SYMMETRY);

  /* Add new boundaries */
  cs_domain_add_boundary(domain, CS_DOMAIN_BOUNDARY_INLET, "left");
  cs_domain_add_boundary(domain, CS_DOMAIN_BOUNDARY_OUTLET, "right");

  /* =========================
     Generic output management
     ========================= */

  cs_domain_set_output_param(domain,
                             100,     // output log frequency
                             3);      // verbosity (-1: no, 0, ...)

  /* ====================
     Time step management
     ==================== */

  /*
    If there is an inconsistency between the max. number of iteration in
    time and the final physical time, the first condition encountered stops
    the calculation.
  */

  cs_domain_set_time_param(domain,
                           500,     // nt_max or -1 (automatic)
                           -1.);    // t_max or < 0. (automatic)

  /* Define the value of the time step
     >> cs_domain_def_time_step_by_value(domain, dt_val);
     >> cs_domain_def_time_step_by_func(domain, dt_func, context);
     context may be NULL.
     This second way to define the time step enable complex definitions.
  */

  cs_domain_def_time_step_by_value(domain, 1.0);

  /* ================================
     Activate groundwater flow module
     ================================

     For the groundwater flow module:
     >> cs_gwf_activate(permeability_type,
                        richards_flag);

     * permeability_type is one of the following keywords:
       CS_PROPERTY_ISO, CS_PROPERTY_ORTHO or CS_PROPERTY_ANISO

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

  /*! [param_cdo_activate_gwf] */
  {

  cs_gwf_activate(CS_PROPERTY_ISO, 0); // no flag to set

  }
  /*! [param_cdo_activate_gwf] */

  /* =========
     Add soils (must be done before adding tracers)
     ========= */

  /*! [param_cdo_gwf_add_soil] */
  {
    cs_gwf_soil_add("soil1", CS_GWF_SOIL_SATURATED);
    cs_gwf_soil_add("soil2", CS_GWF_SOIL_SATURATED);
  }
  /*! [param_cdo_gwf_add_soil] */

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

  /*! [param_cdo_gwf_add_tracer] */
  {
    cs_gwf_add_tracer("Tracer_01","C1");
  }
  /*! [param_cdo_gwf_add_tracer] */

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

  /* Richards equation */
  cs_equation_param_t  *r_eqp = cs_equation_param_by_name("Richards");

  /* Define the boundary conditions  */
  cs_real_t  val0 = 0.0, val1 = 1.0;

  cs_equation_add_bc_by_value(r_eqp,
                              CS_PARAM_BC_DIRICHLET,
                              "left", // boundary zone name
                              &val1);  // value to set

  cs_equation_add_bc_by_value(r_eqp,
                              CS_PARAM_BC_DIRICHLET,
                              "right", // boundary zone name
                              &val0);  // value to set

  /* Tracer equation */
  cs_equation_param_t  *t_eqp = cs_equation_param_by_name("Tracer_01");

  cs_equation_add_bc_by_value(t_eqp,
                              CS_PARAM_BC_DIRICHLET,
                              "left",  // boundary zone name
                              &val1);  // value to set

  /* Define the initial condition (By default: zero is set) */
  cs_equation_add_ic_by_analytic(t_eqp, "cells", get_tracer_sol, NULL);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
