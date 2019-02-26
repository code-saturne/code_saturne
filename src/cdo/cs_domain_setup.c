/* ===========================================================================
 * Routines to handle the setup of a computational domain
 * High level interface
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

#include "cs_ale.h"
#include "cs_boundary_zone.h"
#include "cs_evaluate.h"
#include "cs_equation.h"
#include "cs_equation_common.h"
#include "cs_equation_param.h"
#include "cs_gwf.h"
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_log_iteration.h"
#include "cs_mesh_deform.h"
#include "cs_mesh_location.h"
#include "cs_navsto_system.h"
#include "cs_parall.h"
#include "cs_prototypes.h"
#include "cs_source_term.h"
#include "cs_time_step.h"
#include "cs_walldistance.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_domain_setup.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
 *  \file cs_domain_setup.c
 *
 * \brief  Routines to handle the setup of a computational domain
 *         High level interface for handling the computation.
 */

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local static variables
 *============================================================================*/

static const char _err_empty_domain[] =
  " Stop setting an empty cs_domain_t structure.\n"
  " Please check your settings.\n";

static const char _err_empty_cdo_context[] =
  " Stop setting an empty cs_domain_cdo_context_t structure.\n"
  " Please check your settings.\n";

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the scheme flags for the current computational domain
 *         Requirement: domain->cdo_context is alloctated
 *
 * \param[in, out]  domain       pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_set_scheme_flags(cs_domain_t    *domain)
{
  if (domain == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_domain);
  if (domain->cdo_context == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_cdo_context);

  cs_domain_cdo_context_t  *cc = domain->cdo_context;

  /* Define a scheme flag for the current domain */
  const int  n_equations = cs_equation_get_n_equations();
  for (int eq_id = 0; eq_id < n_equations; eq_id++) {

    cs_equation_t  *eq = cs_equation_by_id(eq_id);
    cs_param_space_scheme_t  scheme = cs_equation_get_space_scheme(eq);
    int  vardim = cs_equation_get_var_dim(eq);

    switch (scheme) {

    case CS_SPACE_SCHEME_CDOVB:
      cc->vb_scheme_flag |= CS_FLAG_SCHEME_POLY0;
      if (vardim == 1)
        cc->vb_scheme_flag |= CS_FLAG_SCHEME_SCALAR;
      else if (vardim == 3)
        cc->vb_scheme_flag |= CS_FLAG_SCHEME_VECTOR;
      else
        bft_error(__FILE__, __LINE__, 0, "Invalid case");
      break;

    case CS_SPACE_SCHEME_CDOVCB:
      cc->vcb_scheme_flag |= CS_FLAG_SCHEME_POLY0;
      if (vardim == 1)
        cc->vcb_scheme_flag |= CS_FLAG_SCHEME_SCALAR;
      else if (vardim == 3)
        cc->vcb_scheme_flag |= CS_FLAG_SCHEME_VECTOR;
      else
        bft_error(__FILE__, __LINE__, 0, "Invalid case");
      break;

    case CS_SPACE_SCHEME_CDOFB:
      cc->fb_scheme_flag |= CS_FLAG_SCHEME_POLY0;
      if (vardim == 1)
        cc->fb_scheme_flag |= CS_FLAG_SCHEME_SCALAR;
      else if (vardim == 3)
        cc->fb_scheme_flag |= CS_FLAG_SCHEME_VECTOR;
      else
        bft_error(__FILE__, __LINE__, 0, "Invalid case");
      break;

    case CS_SPACE_SCHEME_HHO_P0:
      assert(cs_equation_get_space_poly_degree(eq) == 0);
      cc->hho_scheme_flag |= CS_FLAG_SCHEME_POLY0;
      if (vardim == 1)
        cc->hho_scheme_flag |= CS_FLAG_SCHEME_SCALAR;
      else if (vardim == 3)
        cc->hho_scheme_flag |= CS_FLAG_SCHEME_VECTOR;
      else
        bft_error(__FILE__, __LINE__, 0, "Invalid case");
      break;

    case CS_SPACE_SCHEME_HHO_P1:
      cc->hho_scheme_flag |= CS_FLAG_SCHEME_POLY1;
      assert(cs_equation_get_space_poly_degree(eq) == 1);
      if (vardim == 1)
        cc->hho_scheme_flag |= CS_FLAG_SCHEME_SCALAR;
      else if (vardim == 3)
        cc->hho_scheme_flag |= CS_FLAG_SCHEME_VECTOR;
      else
        bft_error(__FILE__, __LINE__, 0, "Invalid case");
      break;

    case CS_SPACE_SCHEME_HHO_P2:
      cc->hho_scheme_flag |= CS_FLAG_SCHEME_POLY2;
      assert(cs_equation_get_space_poly_degree(eq) == 2);
      if (vardim == 1)
        cc->hho_scheme_flag |= CS_FLAG_SCHEME_SCALAR;
      else if (vardim == 3)
        cc->hho_scheme_flag |= CS_FLAG_SCHEME_VECTOR;
      else
        bft_error(__FILE__, __LINE__, 0, "Invalid case");
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Undefined type of scheme to solve for eq. %s."
                  " Please check your settings."), cs_equation_get_name(eq));
    }

  } /* Loop on equations */

  /* Navier-Stokes sytem */
  if (cs_navsto_system_is_activated()) {

    cs_navsto_param_t  *nsp = cs_navsto_system_get_param();

    switch (nsp->space_scheme) {

    case CS_SPACE_SCHEME_CDOVB:
      cc->vb_scheme_flag |= CS_FLAG_SCHEME_NAVSTO;
      break;

    case CS_SPACE_SCHEME_CDOVCB:
      cc->vcb_scheme_flag |= CS_FLAG_SCHEME_NAVSTO;
      break;

    case CS_SPACE_SCHEME_CDOFB:
      cc->fb_scheme_flag |= CS_FLAG_SCHEME_NAVSTO;
      break;

    case CS_SPACE_SCHEME_HHO_P0:
    case CS_SPACE_SCHEME_HHO_P1:
    case CS_SPACE_SCHEME_HHO_P2:
      cc->hho_scheme_flag |= CS_FLAG_SCHEME_NAVSTO;
      break;

    default:
      break;

    }

  } /* NavSto is activated */

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set to true the automatic update of all advection fields
 *
 * \param[in, out]  domain    pointer to a \ref cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_update_advfield(cs_domain_t       *domain)
{
  if (domain == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_domain);
  if (domain->cdo_context == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_cdo_context);

  domain->cdo_context->force_advfield_update = true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set auxiliary parameters related to the way output is done
 *
 * \param[in, out]  domain       pointer to a cs_domain_t structure
 * \param[in]       nt_interval  frequency for the restart process
 * \param[in]       nt_list      output frequency into the log
 * \param[in]       verbosity    level of information displayed
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_set_output_param(cs_domain_t       *domain,
                           int                nt_interval,
                           int                nt_list,
                           int                verbosity)
{
  if (domain == NULL) bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  domain->restart_nt = nt_interval;
  domain->output_nt = nt_list;
  if (domain->output_nt == 0)
    domain->output_nt = -1;

  domain->verbosity = verbosity;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set parameters for unsteady computations: the max number of time
 *         steps or the final physical time of the simulation
 *
 * \param[in, out]  domain    pointer to a \ref cs_domain_t structure
 * \param[in]       nt_max    max. number of time step iterations
 * \param[in]       t_max     final physical time of the simulation
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_set_time_param(cs_domain_t       *domain,
                         int                nt_max,
                         double             t_max)
{
  if (domain == NULL)  bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  domain->time_step->nt_max = nt_max;
  domain->time_step->t_max = t_max;

  /* Add a property related to the time step (in case of use when building a
     linear system) */
  cs_property_add("time_step", CS_PROPERTY_ISO);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the value of the time step thanks to a predefined function
 *
 * \param[in, out] domain      pointer to a cs_domain_t structure
 * \param[in]      func        pointer to a cs_time_func_t function
 * \param[in]      func_input  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_def_time_step_by_function(cs_domain_t        *domain,
                                    cs_time_func_t     *func,
                                    void               *func_input)
{
  if (domain == NULL) bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  domain->time_step->is_variable = 1; /* not constant time step */
  domain->time_options.idtvar = 1;    /* uniform in space but can change
                                         from one time step to the other */

  cs_xdef_time_func_input_t  def = {.input = func_input,
                                    .func = func};

  domain->time_step_def = cs_xdef_timestep_create(CS_XDEF_BY_TIME_FUNCTION,
                                                  0,     /* state flag */
                                                  0,     /* meta flag */
                                                  &def);

  /* Set the property related to the time step if used for building a system */
  cs_property_def_by_time_func(cs_property_by_name("time_step"),
                               NULL, /* all cells are selected */
                               func, func_input);

  /* Default initialization.
     To be changed at first call to cs_domain_time_step_increment() */

  domain->time_step->dt[0] = domain->time_step->t_max;
  domain->time_step->dt_ref =  domain->time_step->t_max;
  domain->time_options.dtmin =  domain->time_step->t_max;
  domain->time_options.dtmax = 0.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the value of the time step.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 * \param[in]        dt        value of the constant time step
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_def_time_step_by_value(cs_domain_t   *domain,
                                 double         dt)
{
  if (domain == NULL) bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  domain->time_step->is_variable = 0; /* constant time step */
  domain->time_options.idtvar = 0;    /* constant time step by default */

  domain->time_step_def = cs_xdef_timestep_create(CS_XDEF_BY_VALUE,
                                                  0, /* state flag */
                                                  0, /* meta flag */
                                                  &dt);

  domain->time_step->dt[0] = dt;
  domain->time_step->dt_ref = dt;
  domain->time_options.dtmin = dt;
  domain->time_options.dtmax = dt;

  /* Set the property related to the time step if used for building a system */
  cs_property_def_iso_by_value(cs_property_by_name("time_step"), NULL, dt);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  First setup stage of the cs_domain_t structure
 *         Define extra domain boundaries
 *         Setup predefined equations
 *         Create fields
 *         Define cs_sles_t structures for variable fields
 *
 * \param[in, out]  domain    pointer to a cs_domain_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_initialize_setup(cs_domain_t    *domain)
{
  /* Setup predefined equations which are activated. At this stage,
   * no equation is added. Space discretization scheme and the related
   * numerical parameters are set.
   */

  /* Wall distance */
  if (cs_walldistance_is_activated())
    cs_walldistance_setup();

  /* Mesh deformation */
  if (cs_mesh_deform_is_activated())
    cs_mesh_deform_setup(domain);

  /* Groundwater flow module */
  if (cs_gwf_is_activated())
    cs_gwf_init_setup();

  /* Navier-Stokes system */
  if (cs_navsto_system_is_activated())
    cs_navsto_system_init_setup();

  /* ALE mesh velocity */
  if (cs_ale_is_activated())
    cs_ale_init_setup(domain);

  /* Add variables related to user-defined and predefined equations */
  cs_equation_create_fields();
  cs_advection_field_create_fields();

  /* Set the scheme flag for the computational domain */
  _set_scheme_flags(domain);

  /* Last step: Proceed to the settings of a cs_equation_t structure. Setup the
     structure related to cs_sles_*
  */

  /* Navier-Stokes system */
  if (cs_navsto_system_is_activated())
    cs_navsto_system_set_sles();

  /* Set the remaining equations */
  cs_equation_set_sles();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Last setup stage of the cs_domain_t structure
 *
 * \param[in, out]  domain            pointer to a cs_domain_t struct.
 * \param[in, out]  mesh              pointer to a cs_mesh_t struct.
 * \param[in]       mesh_quantities   pointer to a cs_mesh_quantities_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_finalize_setup(cs_domain_t                 *domain,
                         cs_mesh_t                   *mesh,
                         const cs_mesh_quantities_t  *mesh_quantities)
{
  if (domain == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_domain);
  if (domain->cdo_context == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_cdo_context);

  /* Manage checkpoint/restart settings
   * Use the same default values for t_interval and wt_interval as the FV */
  double  t_interval = -1.0, wt_interval = -1.0;
  cs_restart_checkpoint_set_defaults(domain->restart_nt,
                                     t_interval,
                                     wt_interval);

  /* Build additional connectivity structures
     Update mesh structure with range set structures */
  cs_domain_cdo_context_t  *cc = domain->cdo_context;
  domain->connect = cs_cdo_connect_init(mesh,
                                        cc->vb_scheme_flag,
                                        cc->vcb_scheme_flag,
                                        cc->fb_scheme_flag,
                                        cc->hho_scheme_flag);

  /* Build additional mesh quantities in a seperate structure */
  domain->cdo_quantities =  cs_cdo_quantities_build(mesh,
                                                    mesh_quantities,
                                                    domain->connect);

  /* Shared main generic structure
     Avoid the declaration of global variables by sharing pointers */
  cs_source_term_set_shared_pointers(domain->cdo_quantities,
                                     domain->connect);
  cs_evaluate_set_shared_pointers(domain->cdo_quantities,
                                  domain->connect);
  cs_property_set_shared_pointers(domain->cdo_quantities,
                                  domain->connect);
  cs_advection_field_set_shared_pointers(domain->cdo_quantities,
                                         domain->connect);

  /* Groundwater flow module */
  if (cs_gwf_is_activated()) {

    /* Setup for the soil structures and the tracer equations */
    cs_user_gwf_setup(domain);

    /* Add if needed new terms (as diffusion or reaction) to tracer equations
       according to the settings */
    cs_gwf_add_tracer_terms();

  }

  /* Allocate all fields created during the setup stage */
  cs_field_allocate_or_map_all();

  /* Allocate common structures for solving equations */
  cs_equation_allocate_structures(domain->connect,
                                  domain->cdo_quantities,
                                  domain->time_step,
                                  cc->vb_scheme_flag,
                                  cc->vcb_scheme_flag,
                                  cc->fb_scheme_flag,
                                  cc->hho_scheme_flag);

  /* Set the range set structure for synchronization in parallel computing */
  cs_equation_set_range_set(domain->connect);

  cs_equation_set_shared_structures(domain->connect,
                                    domain->cdo_quantities,
                                    domain->time_step,
                                    cc->vb_scheme_flag,
                                    cc->vcb_scheme_flag,
                                    cc->fb_scheme_flag,
                                    cc->hho_scheme_flag);


  /* Set the definition of user-defined properties and/or advection
   * fields (no more fields are created at this stage)
   * Last setting stage for equations: Associate properties to activate or not
   * terms such as:
   * --> Unsteady terms (link with a property)
   * --> Diffusion terms (link with a property)
   * --> Advection term  (link with an advection field)
   * --> Reaction term (link with a property)
   * --> Source term
   */

  cs_user_finalize_setup(cs_glob_domain);

  /* Assign to a cs_equation_t structure a list of function to manage this
   * structure during the computation.
   * The set of functions chosen for each equation depends on the parameters
   * specifying the cs_equation_t structure
   */

  domain->only_steady = cs_equation_set_functions();
  if (domain->only_steady)
    domain->is_last_iter = true;

  /* Last stage for the settings for each predefined set of equations:
     - wall distance computation
     - groundwater flow module
     - Navier-Stokes system
   */

  if (cs_walldistance_is_activated())
    cs_walldistance_finalize_setup(domain->connect, domain->cdo_quantities);

  if (cs_gwf_is_activated())
    cs_gwf_finalize_setup(domain->connect, domain->cdo_quantities);

  if (cs_navsto_system_is_activated())
    cs_navsto_system_finalize_setup(domain->mesh,
                                    domain->connect,
                                    domain->cdo_quantities,
                                    domain->time_step);

  if (cs_ale_is_activated())
    cs_ale_finalize_setup(domain);

  /* Last stage to define properties (when complex definition is requested) */
  cs_property_finalize_setup();

  /* Last stage to define properties (when complex definition is requested) */
  cs_advection_field_finalize_setup();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize systems of equations and their related field values
 *         according to the user settings
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_initialize_systems(cs_domain_t   *domain)
{
  /* Initialize system before resolution for all equations
     - create system builder
     - initialize field according to initial conditions
     - initialize source term
     - set the initial condition to all variable fields */
  cs_equation_initialize(domain->mesh,
                         domain->connect,
                         domain->cdo_quantities,
                         domain->time_step);

  /* Set the initial condition for all advection fields */
  cs_advection_field_update(domain->time_step->t_cur,
                            false); /* operate current to previous ? */

  /* Set the initial state for the groundawater flow module */
  if (cs_navsto_system_is_activated())
    cs_navsto_system_initialize(domain->mesh,
                                domain->connect,
                                domain->cdo_quantities,
                                domain->time_step);

  /* Set the initial state for the groundawater flow module */
  if (cs_gwf_is_activated())
    cs_gwf_update(domain->mesh,
                  domain->connect,
                  domain->cdo_quantities,
                  domain->time_step,
                  false); /* operate current to previous ? */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
