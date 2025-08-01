/* ===========================================================================
 * Functions to handle the setup of a computational domain
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_mem.h"

#include "base/cs_ale.h"
#include "base/cs_boundary_zone.h"
#include "base/cs_log.h"
#include "base/cs_log_iteration.h"
#include "base/cs_parall.h"
#include "base/cs_pressure_correction.h"
#include "base/cs_prototypes.h"
#include "base/cs_time_step.h"
#include "cdo/cs_cdo_blas.h"
#include "cdo/cs_cdo_system.h"
#include "cdo/cs_equation.h"
#include "cdo/cs_equation_param.h"
#include "cdo/cs_evaluate.h"
#include "cdo/cs_hodge.h"
#include "cdo/cs_maxwell.h"
#include "cdo/cs_mesh_deform.h"
#include "cdo/cs_navsto_system.h"
#include "cdo/cs_solidification.h"
#include "cdo/cs_source_term.h"
#include "cdo/cs_thermal_system.h"
#include "cdo/cs_walldistance.h"
#include "gui/cs_gui.h"
#include "gwf/cs_gwf.h"
#include "mesh/cs_mesh_location.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cdo/cs_domain_setup.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
 *  \file cs_domain_setup.cpp
 *
 * \brief  Functions to handle the setup of a computational domain
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
  if (domain == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_domain);
  if (domain->cdo_context == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_cdo_context);

  cs_flag_t  quant_flag = 0;
  cs_domain_cdo_context_t  *cc = domain->cdo_context;

  /* Define a scheme flag for the current domain */

  const int  n_equations = cs_equation_get_n_equations();
  for (int eq_id = 0; eq_id < n_equations; eq_id++) {

    cs_equation_t  *eq = cs_equation_by_id(eq_id);
    cs_param_space_scheme_t  scheme = cs_equation_get_space_scheme(eq);
    int  vardim = cs_equation_get_var_dim(eq);

    switch (scheme) {

    case CS_SPACE_SCHEME_CDOVB:
      quant_flag |= CS_CDO_QUANTITIES_VB_SCHEME;
      cc->vb_scheme_flag |= CS_FLAG_SCHEME_POLY0;
      if (vardim == 1)
        cc->vb_scheme_flag |= CS_FLAG_SCHEME_SCALAR;
      else if (vardim == 3)
        cc->vb_scheme_flag |= CS_FLAG_SCHEME_VECTOR;
      else
        bft_error(__FILE__, __LINE__, 0, "Invalid case");
      break;

    case CS_SPACE_SCHEME_CDOVCB:
      quant_flag |= CS_CDO_QUANTITIES_VCB_SCHEME;
      cc->vcb_scheme_flag |= CS_FLAG_SCHEME_POLY0;
      if (vardim == 1)
        cc->vcb_scheme_flag |= CS_FLAG_SCHEME_SCALAR;
      else if (vardim == 3)
        cc->vcb_scheme_flag |= CS_FLAG_SCHEME_VECTOR;
      else
        bft_error(__FILE__, __LINE__, 0, "Invalid case");
      break;

    case CS_SPACE_SCHEME_CDOEB:
      assert(vardim == 3);
      quant_flag |= CS_CDO_QUANTITIES_EB_SCHEME;

      /* vardim should equal to 3 but each edge is associated to a scalar-valued
         quantity */

      cc->eb_scheme_flag |= CS_FLAG_SCHEME_POLY0 | CS_FLAG_SCHEME_SCALAR;
      break;

    case CS_SPACE_SCHEME_CDOFB:
      quant_flag |= CS_CDO_QUANTITIES_FB_SCHEME;
      cc->fb_scheme_flag |= CS_FLAG_SCHEME_POLY0;

      /* Always build quantities related to scalar-valued equations (in
         particular the scalar-valued interface can be useful */

      cc->fb_scheme_flag |= CS_FLAG_SCHEME_SCALAR;
      if (vardim == 3)
        cc->fb_scheme_flag |= CS_FLAG_SCHEME_VECTOR;
      else if (vardim > 3)
        bft_error(__FILE__, __LINE__, 0, "Invalid case");
      break;

    case CS_SPACE_SCHEME_CDOCB:
      quant_flag |= CS_CDO_QUANTITIES_CB_SCHEME;
      cc->cb_scheme_flag |= CS_FLAG_SCHEME_POLY0;

      /* Always build quantities related to scalar-valued equations (in
         particular the scalar-valued interface can be useful */

      cc->cb_scheme_flag |= CS_FLAG_SCHEME_SCALAR;
      if (vardim > 1)
        bft_error(__FILE__, __LINE__, 0, "Invalid case");
      break;

    case CS_SPACE_SCHEME_HHO_P0:
      assert(cs_equation_get_space_poly_degree(eq) == 0);
      quant_flag |= CS_CDO_QUANTITIES_HHO_SCHEME;
      cc->hho_scheme_flag |= CS_FLAG_SCHEME_POLY0;
      if (vardim == 1)
        cc->hho_scheme_flag |= CS_FLAG_SCHEME_SCALAR;
      else if (vardim == 3)
        cc->hho_scheme_flag |= CS_FLAG_SCHEME_VECTOR;
      else
        bft_error(__FILE__, __LINE__, 0, "Invalid case");
      break;

    case CS_SPACE_SCHEME_HHO_P1:
      quant_flag |= CS_CDO_QUANTITIES_HHO_SCHEME;
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
      quant_flag |= CS_CDO_QUANTITIES_HHO_SCHEME;
      cc->hho_scheme_flag |= CS_FLAG_SCHEME_POLY2;
      assert(cs_equation_get_space_poly_degree(eq) == 2);
      if (vardim == 1)
        cc->hho_scheme_flag |= CS_FLAG_SCHEME_SCALAR;
      else if (vardim == 3)
        cc->hho_scheme_flag |= CS_FLAG_SCHEME_VECTOR;
      else
        bft_error(__FILE__, __LINE__, 0, "Invalid case");
      break;

    case CS_SPACE_SCHEME_MACFB:
      quant_flag |= CS_CDO_QUANTITIES_MAC_SCHEME;
      cc->mac_scheme_flag |= CS_FLAG_SCHEME_POLY0;

      /* Always build quantities related to scalar-valued equations (in
         particular the scalar-valued interface can be useful */

      cc->mac_scheme_flag |= CS_FLAG_SCHEME_SCALAR;
      if (vardim == 3) {
        cc->mac_scheme_flag |= CS_FLAG_SCHEME_VECTOR;
      }
      else if (vardim > 3) {
        bft_error(__FILE__, __LINE__, 0, "Invalid case");
      }
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Undefined type of scheme to solve for eq. %s."
                  " Please check your settings."), cs_equation_get_name(eq));
    }

  } /* Loop on equations */

  /* Navier-Stokes system */

  if (cs_navsto_system_is_activated()) {

    cs_navsto_param_t  *nsp = cs_navsto_system_get_param();

    switch (nsp->space_scheme) {

    case CS_SPACE_SCHEME_CDOVB:
      quant_flag |= CS_CDO_QUANTITIES_VB_SCHEME;
      cc->vb_scheme_flag |= CS_FLAG_SCHEME_NAVSTO;
      break;

    case CS_SPACE_SCHEME_CDOVCB:
      quant_flag |= CS_CDO_QUANTITIES_VCB_SCHEME;
      cc->vcb_scheme_flag |= CS_FLAG_SCHEME_NAVSTO;
      break;

    case CS_SPACE_SCHEME_CDOEB:
      quant_flag |= CS_CDO_QUANTITIES_EB_SCHEME;
      cc->eb_scheme_flag |= CS_FLAG_SCHEME_NAVSTO;
      break;

    case CS_SPACE_SCHEME_CDOFB:
      quant_flag |= CS_CDO_QUANTITIES_FB_SCHEME;
      cc->fb_scheme_flag |= CS_FLAG_SCHEME_NAVSTO;
      break;

    case CS_SPACE_SCHEME_HHO_P0:
    case CS_SPACE_SCHEME_HHO_P1:
    case CS_SPACE_SCHEME_HHO_P2:
      quant_flag |= CS_CDO_QUANTITIES_HHO_SCHEME;
      cc->hho_scheme_flag |= CS_FLAG_SCHEME_NAVSTO;
      break;

    case CS_SPACE_SCHEME_MACFB:
      quant_flag |= CS_CDO_QUANTITIES_MAC_SCHEME;
      cc->mac_scheme_flag |= CS_FLAG_SCHEME_NAVSTO;
      break;

    default:
      break;

    }

  } /* NavSto is activated */

  /* Update the flag storing which geometrical quantities have to be
     computed */

  cs_cdo_quantities_set(quant_flag);
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

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
  if (domain == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  domain->time_step->nt_max = nt_max;
  domain->time_step->t_max = t_max;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set time step parameters for unsteady computations when this is not
 *         already done. This situation should occur when the GUI is used to
 *         set a constant time step.
 *
 * \param[in, out]  domain    pointer to a \ref cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_automatic_time_step_settings(cs_domain_t       *domain)
{
  if (domain == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  cs_time_step_t  *ts = domain->time_step;

  if (ts->t_max < 0 && ts->nt_max < 1)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Please check your settings.\n"
              " Unsteady computation but no definition available.\n",
              __func__);
  if (ts->dt_ref < 0)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Please check your settings.\n"
              " Unsteady computation but no dt_ref available.\n",
              __func__);

  cs_domain_set_time_param(domain, ts->nt_max, ts->t_max);
  cs_domain_def_time_step_by_value(domain, ts->dt_ref);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the value of the time step thanks to a predefined function
 *
 * \param[in, out] domain      pointer to a cs_domain_t structure
 * \param[in]      func        pointer to a cs_time_func_t function
 * \param[in]      func_input  pointer to a structure cast on-the-fly
 *
 * \return a pointer to the created definition (\ref cs_xdef_t structure)
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_domain_def_time_step_by_function(cs_domain_t    *domain,
                                    cs_time_func_t *func,
                                    void           *func_input)
{
  if (domain == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  domain->time_step->is_variable = 1; /* not constant time step */

  /* Uniform in space but can change from one time step to the other */

  domain->time_options->idtvar = CS_TIME_STEP_ADAPTIVE;

  /* Set the property related to the time step if used for building a system */

  cs_property_t  *dt_pty = cs_property_by_name("time_step");
  assert(dt_pty != nullptr);

  cs_property_set_reference_value(dt_pty, domain->time_step->t_max);

  cs_xdef_t *def = cs_property_def_by_time_func(dt_pty,
                                                nullptr, /* all cells */
                                                func,
                                                func_input);

  /* Default initialization.
     To be changed at first call to cs_domain_time_step_increment() */

  domain->time_step->dt[0] = domain->time_step->t_max;
  domain->time_step->dt[1] = domain->time_step->t_max;
  domain->time_step->dt[2] = domain->time_step->t_max;
  domain->time_step->dt_ref = -1; // Updated at initialization
  domain->time_options->dtmin = domain->time_step->t_max;
  domain->time_options->dtmax = 0.;

  return def;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the value of the time step.
 *
 * \param[in, out] domain  pointer to a cs_domain_t structure
 * \param[in]      dt      value of the constant time step
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_def_time_step_by_value(cs_domain_t *domain,
                                 double       dt)
{
  if (domain == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  domain->time_step->is_variable = 0; /* constant time step */

  /* Constant time step by default */

  domain->time_options->idtvar = CS_TIME_STEP_CONSTANT;

  domain->time_step->dt[0] = dt;    /* time step n */
  domain->time_step->dt[1] = dt;    /* time step n-1 */
  domain->time_step->dt[2] = dt;    /* time step n-2 */
  domain->time_step->dt_ref = dt;
  domain->time_step->dt_next = dt;
  domain->time_options->dtmin = dt;
  domain->time_options->dtmax = dt;

  /* Set the property related to the time step if used for building a system */

  cs_property_t  *dt_pty = cs_property_by_name("time_step");
  assert(dt_pty != nullptr);

  cs_property_def_constant_value(dt_pty, dt); // ref. value set inside
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief After having read the mesh and the first setup stage build the
 *        connectivities and mesh quantities related to CDO/HHO schemes
 *
 * \param[in, out] domain  pointer to a cs_domain_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_setup_init_cdo_structures(cs_domain_t *domain)
{
  if (domain == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_domain);
  if (domain->cdo_context == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_cdo_context);

  /* Manage checkpoint/restart settings
   * Use the same default values for t_interval and wt_interval as the FV */

  int nt_interval = -1;
  double  t_interval = -1.0, wt_interval = -1.0;

  cs_restart_checkpoint_get_intervals(&nt_interval, &t_interval, &wt_interval);
  if (cs_glob_time_step->t_max > -1.) {
    /* Computation stop depends on time, we set default interval based on that
     * instead of number of iterations.
     */
    if (nt_interval == CS_RESTART_INTERVAL_DEFAULT) {
      /* default behavior => nt_interval is set to default by the GUI! */
      t_interval = (cs_glob_time_step->t_max - cs_glob_time_step->t_prev) * 0.25;
      nt_interval = -1;
    }
  }
  cs_restart_checkpoint_set_interval(nt_interval, t_interval, wt_interval);

  cs_restart_checkpoint_set_interval(domain->restart_nt,
                                     t_interval,
                                     wt_interval);

  /* Build additional connectivity structures
     Update mesh structure with range set structures */

  cs_domain_cdo_context_t  *cc = domain->cdo_context;
  domain->connect = cs_cdo_connect_build(domain->mesh,
                                         cc->eb_scheme_flag,
                                         cc->fb_scheme_flag,
                                         cc->cb_scheme_flag,
                                         cc->vb_scheme_flag,
                                         cc->vcb_scheme_flag,
                                         cc->hho_scheme_flag,
                                         cc->mac_scheme_flag);

  /* Build additional mesh quantities in a separate structure */

  cs_flag_t  cdo_quantities_flag = 0;
  if (cc->eb_scheme_flag)
    cdo_quantities_flag |= CS_CDO_QUANTITIES_EB_SCHEME;
  if (cc->fb_scheme_flag)
    cdo_quantities_flag |= CS_CDO_QUANTITIES_FB_SCHEME;
  if (cc->cb_scheme_flag)
    cdo_quantities_flag |= CS_CDO_QUANTITIES_CB_SCHEME;
  if (cc->hho_scheme_flag)
    cdo_quantities_flag |= CS_CDO_QUANTITIES_HHO_SCHEME;
  if (cc->vb_scheme_flag)
    cdo_quantities_flag |= CS_CDO_QUANTITIES_VB_SCHEME;
  if (cc->vcb_scheme_flag)
    cdo_quantities_flag |= CS_CDO_QUANTITIES_VCB_SCHEME;
  if (cc->mac_scheme_flag)
    cdo_quantities_flag |= CS_CDO_QUANTITIES_MAC_SCHEME;

  cs_cdo_quantities_set(cdo_quantities_flag);

  domain->cdo_quantities = cs_cdo_quantities_build(domain->mesh,
                                                   domain->mesh_quantities,
                                                   domain->connect);

  /* Main generic structures are shared with low-level files.
     Avoid the declaration of global variables by sharing pointers */

  cs_advection_field_init_sharing(domain->cdo_quantities, domain->connect);

  cs_cdo_blas_init_sharing(domain->cdo_quantities, domain->connect);

  cs_cdo_system_init_sharing(domain->mesh, domain->connect);

  cs_evaluate_init_sharing(domain->cdo_quantities,
                           domain->connect,
                           domain->mesh);

  cs_equation_init_sharing(domain->mesh,
                           domain->connect,
                           domain->cdo_quantities,
                           domain->time_step,
                           cc->cb_scheme_flag,
                           cc->eb_scheme_flag,
                           cc->fb_scheme_flag,
                           cc->vb_scheme_flag,
                           cc->vcb_scheme_flag,
                           cc->hho_scheme_flag,
                           cc->mac_scheme_flag);

  cs_equation_system_init_sharing(domain->mesh,
                                  domain->connect,
                                  domain->cdo_quantities,
                                  domain->time_step);

  cs_property_init_sharing(domain->mesh,
                           domain->cdo_quantities,
                           domain->connect);

  cs_source_term_init_sharing(domain->cdo_quantities, domain->connect);

  cs_gwf_init_sharing(domain->cdo_quantities,
                      domain->connect);

  /* Allocate common local structures for building/solving equations */

  cs_cdo_toolbox_init(domain->connect,
                      cc->eb_scheme_flag,
                      cc->fb_scheme_flag,
                      cc->cb_scheme_flag,
                      cc->vb_scheme_flag,
                      cc->vcb_scheme_flag,
                      cc->hho_scheme_flag,
                      cc->mac_scheme_flag);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief First setup stage of the cs_domain_t structure
 *        Define extra domain boundaries
 *        Setup predefined equations
 *        Create fields (already done in the FV part)
 *        Define cs_sles_t structures for variable fields
 *
 * \param[in, out] domain  pointer to a cs_domain_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_setup_initialize(cs_domain_t *domain)
{
  /* Setup predefined equations which are activated. At this stage, no equation
   * is added. Space discretization scheme and the related numerical parameters
   * are set.
   */

  /* User equations:
   *  Add a variable field related to each user-defined equation.
   */

  cs_equation_user_create_fields();

  /* For predefined equations (equations associated to a module), adding a
   * variable field is done inside each setup function of the module */

  /* Hybrid velcoity/pressure coupling between FV and CDO */

  if (cs_pressure_correction_cdo_is_activated())
    cs_pressure_correction_cdo_init_setup();

  /* Wall distance */

  if (cs_walldistance_is_activated())
    cs_walldistance_setup();

  /* Mesh deformation */

  if (cs_mesh_deform_is_activated())
    cs_mesh_deform_setup(domain);

  /* Thermal module */

  if (cs_thermal_system_is_activated())
    cs_thermal_system_init_setup();

  /* Groundwater flow module */

  if (cs_gwf_is_activated())
    cs_gwf_init_setup();

  /* ALE mesh velocity */

  if (cs_ale_is_activated())
    cs_ale_init_setup(domain);

  /* Maxwell module */

  if (cs_maxwell_is_activated())
    cs_maxwell_init_setup();

  /* Navier-Stokes system */

  if (cs_navsto_system_is_activated()) {

    /* To make more easy the settings for the end-user, one may have to ensure
     * that the Navier-Stokes system has the sufficient knowledge of what is
     * requested */

    if (cs_thermal_system_needs_navsto())
      cs_navsto_system_update_model(true); /* true = with thermal */

    cs_navsto_system_init_setup();

  }
  else {

    /* Switch off the turbulence modelling if in CDO mode only */

    if (cs_param_cdo_has_cdo_only()) {
      cs_turb_model_t  *turb = cs_get_glob_turb_model();

      turb->model = CS_TURB_NONE;          /* laminar flow */
      turb->itytur = 0;                    /* deprecated */
      turb->hybrid_turb = CS_HYBRID_NONE;
      turb->type        = CS_TURB_NONE;
    }
  }

  if (cs_solidification_is_activated())
    cs_solidification_init_setup();

  /* Add fields associated to advection fields */

  cs_advection_field_create_fields();

  /* Set the scheme flag for the computational domain */

  _set_scheme_flags(domain);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Last GUI and then user setup stage of the cs_domain_t structure
 *
 * \param[in, out] domain  pointer to a cs_domain_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_setup_finalize(cs_domain_t  *domain)
{
  if (domain == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_domain);
  if (domain->cdo_context == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_cdo_context);

  /* Allocate all fields created during the setup stage */

  cs_field_allocate_or_map_all();

  /* GUI initial conditions - Only Conjugate Heat Transfer for the moment */

  if (cs_thermal_system_is_activated())
    cs_gui_initial_conditions();

  // User-defined initialization

  cs_user_initialization(domain);

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

  cs_gui_thermal_source_terms_setup();
  cs_user_boundary_conditions_setup(domain);
  cs_user_finalize_setup(domain);

  /* Setup boundary conditions for CHT coupling if needed */

  if (cs_thermal_system_is_activated())
    cs_thermal_system_cht_boundary_conditions_setup();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Last user setup stage of the cs_domain_t structure
 *
 * \param[in, out] domain  pointer to a cs_domain_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_setup_finalize_module(cs_domain_t *domain)
{
  if (domain == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_domain);
  if (domain->cdo_context == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_cdo_context);

  /* Last stage for the settings for each predefined set of equations:
     - thermal module
     - groundwater flow module
     - Maxwell equations
     - Navier-Stokes system
     - ALE equation
   */

  if (cs_thermal_system_is_activated())
    cs_thermal_system_finalize_setup(domain->connect,
                                     domain->cdo_quantities,
                                     domain->time_step);

  if (cs_gwf_is_activated())
    cs_gwf_finalize_setup(domain->connect, domain->cdo_quantities);

  if (cs_maxwell_is_activated())
    cs_maxwell_finalize_setup(domain->connect, domain->cdo_quantities);

  if (cs_navsto_system_is_activated())
    cs_navsto_system_finalize_setup(domain->mesh,
                                    domain->connect,
                                    domain->cdo_quantities,
                                    domain->time_step);

  if (cs_ale_is_activated())
    cs_ale_finalize_setup(domain);

  if (cs_solidification_is_activated())
    cs_solidification_finalize_setup(domain->connect,
                                     domain->cdo_quantities);

  if (cs_pressure_correction_cdo_is_activated())
    cs_pressure_correction_cdo_finalize_setup(domain);

  /* Last stage to define properties (when complex definition is requested) */

  cs_property_finalize_setup();

  /* Last stage to define properties (when complex definition is requested) */

  cs_advection_field_finalize_setup();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize systems of equations and their related field/properties
 *        values according to the user settings
 *
 * \param[in, out] domain  pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_setup_init_state(cs_domain_t *domain)
{
  assert(domain != nullptr);

  /* Set the initial condition to all variable fields */

  cs_equation_init_field_values(domain->mesh, domain->time_step);

  /* Set the initial condition for all advection fields */

  cs_advection_field_update(domain->time_step->t_cur,
                            false); /* operate current to previous ? */

  /* Set the initial state for the thermal module */

  if (cs_thermal_system_is_activated())
    cs_thermal_system_update(domain->mesh,
                             domain->connect,
                             domain->cdo_quantities,
                             domain->time_step,
                             false); /* operate current to previous ? */

  /* Set the initial state for the Navier-Stokes system */

  if (cs_navsto_system_is_activated())
    cs_navsto_system_init_values(domain->mesh,
                                 domain->connect,
                                 domain->cdo_quantities,
                                 domain->time_step);

  /* Set the initial state for the Maxwell module */

  if (cs_maxwell_is_activated())
    cs_maxwell_update(domain->mesh,
                      domain->connect,
                      domain->cdo_quantities,
                      domain->time_step,
                      false); /* operate current to previous ? */

  /* Additional initializations (Navier-Stokes system has already been
     initialized) */

  if (cs_solidification_is_activated())
    cs_solidification_init_values(domain->mesh,
                                  domain->connect,
                                  domain->cdo_quantities,
                                  domain->time_step);

  /* Set the initial state for the groundwater flow module */

  if (cs_gwf_is_activated())
    cs_gwf_init_values(domain->mesh,
                       domain->connect,
                       domain->cdo_quantities,
                       domain->time_step);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of the main domain settings
 *
 * \param[in]   domain    pointer to the cs_domain_t structure to summarize
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_setup_log(const cs_domain_t   *domain)
{
  cs_log_printf(CS_LOG_SETUP, "\n Working directory names\n");
  cs_log_printf(CS_LOG_SETUP, "%s\n", cs_sep_h1);

  cs_log_printf(CS_LOG_SETUP, " %10s: %s\n", "Study", domain->study_name);
  cs_log_printf(CS_LOG_SETUP, " %10s: %s\n", "Case", domain->case_name);
  cs_log_printf(CS_LOG_SETUP, " %10s: %s\n", "Run_id", domain->run_id);

  cs_log_printf(CS_LOG_SETUP, "\nSummary of the CDO domain settings\n");
  cs_log_printf(CS_LOG_SETUP, "%s\n", cs_sep_h1);

  cs_param_cdo_setup_log();

  cs_cdo_assembly_setup_log();

  /* CDO main structure count */

  cs_log_printf(CS_LOG_SETUP, "\n## CDO main structures\n");

  int  n_equations, n_predef_equations, n_user_equations;
  cs_equation_get_count(&n_equations, &n_predef_equations, &n_user_equations);

  int  n_systems = cs_equation_system_get_n_systems();

  cs_log_printf(CS_LOG_SETUP, " * Number of systems of equations       %3d\n",
                n_systems);
  cs_log_printf(CS_LOG_SETUP, " * Number of equations                  %3d\n",
                n_equations);
  cs_log_printf(CS_LOG_SETUP, "  - with number of predefined equations %3d\n",
                n_predef_equations);
  cs_log_printf(CS_LOG_SETUP, "  - with number of user equations       %3d\n",
                n_user_equations);
  cs_log_printf(CS_LOG_SETUP, " * Number of properties                 %3d\n",
                cs_property_get_n_properties());
  cs_log_printf(CS_LOG_SETUP, " * Number of advection fields           %3d\n",
                cs_advection_field_get_n_fields());

  if (cs_param_cdo_has_fv_only())
    return;

  cs_domain_cdo_context_t  *cc = domain->cdo_context;

  cs_cdo_connect_log_summary(domain->connect,
                             cc->eb_scheme_flag,
                             cc->vb_scheme_flag,
                             cc->vcb_scheme_flag);

  cs_cdo_quantities_log_summary(domain->cdo_quantities);

  /* Time step summary */

  cs_log_printf(CS_LOG_SETUP, "\n## Time step information\n");
  if (domain->only_steady)
    cs_log_printf(CS_LOG_SETUP, " * Steady-state computation\n");

  else { /* Time information */

    cs_log_printf(CS_LOG_SETUP, " * Unsteady computation\n");

    if (domain->time_step->t_max > 0.)
      cs_log_printf(CS_LOG_SETUP, "%-30s %5.3e\n",
                    " * Final simulation time:", domain->time_step->t_max);
    if (domain->time_step->nt_max > 0)
      cs_log_printf(CS_LOG_SETUP, "%-30s %9d\n",
                    " * Final time step:", domain->time_step->nt_max);

    if (domain->time_options->idtvar == 0)
      cs_log_printf(CS_LOG_SETUP, " * Time step *constant*\n\n");
    else if (domain->time_options->idtvar == 1)
      cs_log_printf(CS_LOG_SETUP, " * Time step *variable in time*\n\n");
    else {
      if (cs_glob_param_cdo_mode != CS_PARAM_CDO_MODE_WITH_FV)
        bft_error(__FILE__, __LINE__, 0,
                  _(" Invalid idtvar value for the CDO module.\n"));
    }

  }
}


/*----------------------------------------------------------------------------*/

END_C_DECLS
