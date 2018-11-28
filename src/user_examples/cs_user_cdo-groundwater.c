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

#include "cs_advection_field.h"
#include "cs_boundary_zone.h"
#include "cs_gwf.h"
#include "cs_mesh_location.h"
#include "cs_property.h"
#include "cs_sdm.h"
#include "cs_volume_zone.h"

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

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

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

  cs_boundary_t  *boundaries = domain->boundaries;

  /* Choose a boundary by default.
     Valid choice is CS_BOUNDARY_WALL or CS_BOUNDARY_SYMMETRY */
  cs_boundary_set_default(boundaries, CS_BOUNDARY_SYMMETRY);

  /* Add a new boundary
     >> cs_domain_add_boundary(domain,
                               type_of_boundary,
                               mesh_location_name);

     * mesh_location_name is either a predefined mesh location or one defined
     by the user
     * type_of_boundary is one of the following keyword
        CS_BOUNDARY_WALL,
        CS_BOUNDARY_INLET,
        CS_BOUNDARY_OUTLET,
        CS_BOUNDARY_SYMMETRY
  */

  cs_boundary_add(boundaries, CS_BOUNDARY_INLET, "left");
  cs_boundary_add(boundaries, CS_BOUNDARY_OUTLET, "right");

  /* =========================
     Generic output management
     ========================= */

  cs_domain_set_output_param(domain,
                             -1,     /* restart frequency */
                             100,    /* output log frequency */
                             3);     /* verbosity (-1: no, 0, ...) */

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

END_C_DECLS
