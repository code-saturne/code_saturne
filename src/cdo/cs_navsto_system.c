/*============================================================================
 * Functions to handle cs_navsto_system_t structure
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_cdofb_ac.h"
#include "cs_cdofb_monolithic.h"
#include "cs_cdofb_monolithic_sles.h"
#include "cs_cdofb_navsto.h"
#include "cs_cdofb_predco.h"
#include "cs_hho_stokes.h"
#include "cs_evaluate.h"
#include "cs_flag.h"
#include "cs_log.h"
#include "cs_navsto_coupling.h"
#include "cs_post.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_navsto_system.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 *  \file cs_navsto_system.c
 *
 *  \brief  Functions to handle the cs_navsto_system_t structure which is
 *          the high-level structure to manage the Navier-Stokes system of
 *          equations.
 */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define CS_NAVSTO_SYSTEM_DBG  0

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private variables
 *============================================================================*/

static const char _err_empty_ns[] =
  " Stop execution. The structure related to the Navier-Stokes system is"
  " empty.\n Please check your settings.\n";

static const char _err_invalid_coupling[] =
  " %s: Invalid case for the coupling algorithm.\n";

static int  cs_navsto_system_solid_enforcement_id = -1;
static cs_navsto_system_t  *cs_navsto_system = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the \ref cs_equation_param_t structure related to the
 *         momentum equation according to the type of coupling
 *
 * \param[in]  nsp       pointer to a \ref cs_navsto_param_t structure
 *
 * \return a pointer to the corresponding \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

static inline bool
_handle_non_linearities(cs_navsto_param_t    *nsp)
{
  if (nsp == NULL)
    return false;

  switch (nsp->model) {

  case CS_NAVSTO_MODEL_OSEEN:
  case CS_NAVSTO_MODEL_STOKES:
    return false;

  case CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES:
    if (nsp->adv_strategy == CS_PARAM_ADVECTION_IMPLICIT_FULL)
      return true;
    else {
      assert(nsp->adv_strategy == CS_PARAM_ADVECTION_IMPLICIT_LINEARIZED ||
             nsp->adv_strategy == CS_PARAM_ADVECTION_EXPLICIT);
      return false;
    }
    break;

  default:
    return true;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate an empty Navier-Stokes system
 *
 * \return a pointer to a new allocated groundwater flow structure
 */
/*----------------------------------------------------------------------------*/

static cs_navsto_system_t *
_allocate_navsto_system(void)
{
  cs_navsto_system_t  *navsto = NULL;

  BFT_MALLOC(navsto, 1, cs_navsto_system_t);

  navsto->param = NULL;

  /* Array of boundary type */

  navsto->bf_type = NULL;

  /* Main set of variables */

  navsto->velocity = NULL;
  navsto->pressure = NULL;

  /* Advection field */

  navsto->adv_field = NULL;
  navsto->mass_flux_array = NULL;
  navsto->mass_flux_array_pre = NULL;

  /* Related modules */

  navsto->turbulence = NULL;

  /* Post-processing fields */

  navsto->plot_writer = NULL;
  navsto->velocity_divergence = NULL;
  navsto->pressure_gradient = NULL;
  navsto->mass_density = NULL;
  navsto->mass_flux_balance = NULL;
  navsto->kinetic_energy = NULL;
  navsto->velocity_gradient = NULL;
  navsto->vorticity = NULL;
  navsto->helicity = NULL;
  navsto->enstrophy = NULL;

  /* Stream function is associated to the variable field of an equation
     So the treatment is different */

  navsto->stream_function_eq = NULL;

  /* Additional data fitting the choice of the coupling model */

  navsto->coupling_context = NULL;
  navsto->scheme_context = NULL;

  /* Function pointers */

  navsto->init_scheme_context = NULL;
  navsto->free_scheme_context = NULL;
  navsto->init_velocity = NULL;
  navsto->init_pressure = NULL;
  navsto->compute_steady = NULL;
  navsto->compute= NULL;

  return navsto;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if the resolution of the Navier-Stokes system has been
 *        activated
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_navsto_system_is_activated(void)
{
  if (cs_navsto_system == NULL)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the flag associated to the modelling options
 *
 * \param[in]   with_thermal     true or false
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_update_model(bool   with_thermal)
{
  if (cs_navsto_system == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s: The main structure for the Navier-Stokes equations is not"
              " allocated.\n"
              " Please check your settings", __func__);

  cs_navsto_param_t  *nsp = cs_navsto_system->param;

  if (with_thermal) { /* Thermal system is switch on and relies on the mass flux
                         for the advection */

    if ((nsp->model_flag & (CS_NAVSTO_MODEL_PASSIVE_THERMAL_TRACER |
                            CS_NAVSTO_MODEL_BOUSSINESQ)) == 0) {

      /* Thermal system is linked to the Navier-Stokes one but nothing has been
       * set. Add the "minimal" flag. */

      nsp->model_flag |= CS_NAVSTO_MODEL_PASSIVE_THERMAL_TRACER;

    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize the Navier-Stokes (NS) system
 *
 * \param[in] boundaries     pointer to the domain boundaries
 * \param[in] model          type of model related to the NS system
 * \param[in] model_flag     additional high-level model options
 * \param[in] algo_coupling  algorithm used for solving the NS system
 * \param[in] post_flag      predefined post-processings options
 *
 * \return a pointer to a new allocated cs_navsto_system_t structure
 */
/*----------------------------------------------------------------------------*/

cs_navsto_system_t *
cs_navsto_system_activate(const cs_boundary_t           *boundaries,
                          cs_navsto_param_model_t        model,
                          cs_navsto_param_model_flag_t   model_flag,
                          cs_navsto_param_coupling_t     algo_coupling,
                          cs_navsto_param_post_flag_t    post_flag)
{
  if (model == CS_NAVSTO_N_MODELS)
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid model for Navier-Stokes.\n",
              __func__);

  /* Allocate an empty structure */

  cs_navsto_system_t  *navsto = _allocate_navsto_system();

  /* Initialize the set of parameters */

  navsto->param = cs_navsto_param_create(boundaries,
                                         model,
                                         model_flag,
                                         algo_coupling,
                                         post_flag);

  /* Set the default boundary condition for the equations of the Navier-Stokes
     system according to the default domain boundary */

  cs_param_bc_type_t  default_bc = CS_PARAM_N_BC_TYPES;
  switch (boundaries->default_type) {

  case CS_BOUNDARY_WALL:
    default_bc = CS_PARAM_BC_HMG_DIRICHLET;
    break;
  case CS_BOUNDARY_SYMMETRY:
    default_bc = CS_PARAM_BC_SLIDING;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid boundary default type\n",
              __func__);
    break;

  } /* End of switch */

  /* Advection field related to the resolved velocity */

  cs_advection_field_status_t  adv_status =
    CS_ADVECTION_FIELD_NAVSTO | CS_ADVECTION_FIELD_TYPE_SCALAR_FLUX;

  if (cs_navsto_param_is_steady(navsto->param))
    adv_status |= CS_ADVECTION_FIELD_STEADY;

  navsto->adv_field = cs_advection_field_add("mass_flux", adv_status);

  /* Additional initialization fitting the choice of model */

  switch (navsto->param->coupling) {

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
    navsto->coupling_context =
      cs_navsto_ac_create_context(default_bc, navsto->param);
    break;
  case CS_NAVSTO_COUPLING_MONOLITHIC:
    navsto->coupling_context =
      cs_navsto_monolithic_create_context(default_bc, navsto->param);
    break;
  case CS_NAVSTO_COUPLING_PROJECTION:
    navsto->coupling_context =
      cs_navsto_projection_create_context(default_bc, navsto->param);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _err_invalid_coupling, __func__);
    break;

  }

  if (post_flag & CS_NAVSTO_POST_STREAM_FUNCTION) {

    navsto->stream_function_eq = cs_equation_add(CS_NAVSTO_STREAM_EQNAME,
                                                 "stream_function",
                                                 CS_EQUATION_TYPE_NAVSTO,
                                                 1,
                                                 CS_PARAM_BC_HMG_NEUMANN);

    cs_equation_param_t  *eqp =
      cs_equation_get_param(navsto->stream_function_eq);
    assert(eqp != NULL);

    /* Default settings for this equation */

    cs_equation_param_set(eqp, CS_EQKEY_SPACE_SCHEME, "cdo_vb");
    cs_equation_param_set(eqp, CS_EQKEY_HODGE_DIFF_COEF, "dga");
    cs_equation_param_set(eqp, CS_EQKEY_PRECOND, "amg");
    cs_equation_param_set(eqp, CS_EQKEY_AMG_TYPE, "k_cycle");
    cs_equation_param_set(eqp, CS_EQKEY_ITSOL, "cg");

    /* This is for post-processing purpose, so, there is no need to have
     * a restrictive convergence tolerance on the resolution of the linear
     * system */

    cs_equation_param_set(eqp, CS_EQKEY_ITSOL_EPS, "1e-6");

  }

  /* Create the main structure to handle the turbulence modelling */

  navsto->turbulence = cs_turbulence_create(navsto->param->turbulence);

  /* Set the static variable */

  cs_navsto_system = navsto;

  return navsto;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the main structure related to the Navier-Stokes system
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_destroy(void)
{
  cs_navsto_system_t  *navsto = cs_navsto_system;

  if (navsto == NULL)
    return;

  BFT_FREE(navsto->bf_type);

  /*
   * Properties, advection fields, equations and fields are all destroyed
   * respectively inside cs_property_destroy_all(),
   * cs_advection_field_destroy_all(), cs_equation_destroy_all() and
   * cs_field_destroy_all()
   */

  BFT_FREE(navsto->mass_flux_array);
  BFT_FREE(navsto->mass_flux_array_pre);

  /* Free the plot writer */

  if (navsto->plot_writer != NULL)
    cs_time_plot_finalize(&navsto->plot_writer);

  cs_navsto_param_t  *nsp = navsto->param;

  /* Free the context according to the model choice */

  switch (nsp->coupling) {

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
    navsto->coupling_context =
      cs_navsto_ac_free_context(navsto->coupling_context);
    break;

  case CS_NAVSTO_COUPLING_MONOLITHIC:
    navsto->coupling_context =
      cs_navsto_monolithic_free_context(navsto->coupling_context);
    if (nsp->space_scheme == CS_SPACE_SCHEME_CDOFB) {
      cs_cdofb_monolithic_finalize_common(nsp);
      cs_cdofb_monolithic_sles_finalize();
    }
    break;

  case CS_NAVSTO_COUPLING_PROJECTION:
    navsto->coupling_context =
      cs_navsto_projection_free_context(navsto->coupling_context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _err_invalid_coupling, __func__);
    break;
  }

  /* Destroy the context related to the discretization scheme */

  navsto->free_scheme_context(navsto->scheme_context);

  /* Free the context and the structure related to the turbulence modelling */

  cs_turbulence_free(&(navsto->turbulence));

  /* Set of numerical parameters */

  navsto->param = cs_navsto_param_free(nsp);

  BFT_FREE(navsto);
  cs_navsto_system = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the structure storing the parameters for the Navier--Stokes
 *         system
 *
 * \return NULL or the pointer to a \ref cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_navsto_param_t *
cs_navsto_system_get_param(void)
{
  cs_navsto_system_t  *navsto = cs_navsto_system;

  if (navsto == NULL)
    return NULL;

  return navsto->param;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a pointer to the equation related to the momentum equation
 *
 * \return NULL or the pointer
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_navsto_system_get_momentum_eq(void)
{
  cs_navsto_system_t  *navsto = cs_navsto_system;

  if (navsto == NULL)
    return NULL;

  cs_navsto_param_t  *nsp = navsto->param;
  cs_equation_t  *eq = NULL;

  switch (nsp->coupling) {

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
    eq = cs_navsto_ac_get_momentum_eq(navsto->coupling_context);
    break;
  case CS_NAVSTO_COUPLING_MONOLITHIC:
    eq = cs_navsto_monolithic_get_momentum_eq(navsto->coupling_context);
    break;
  case CS_NAVSTO_COUPLING_PROJECTION:
    eq = cs_navsto_projection_get_momentum_eq(navsto->coupling_context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _err_invalid_coupling, __func__);
    break;

  }

  return eq;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the advection field structure (the mass flux) related to
 *         the Navier-Stokes system.
 *
 * \return a pointer to the advection field structure
 */
/*----------------------------------------------------------------------------*/

cs_adv_field_t *
cs_navsto_get_adv_field(void)
{
  cs_navsto_system_t  *ns = cs_navsto_system;

  if (ns == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_ns));

  return ns->adv_field;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the mass flux array related to the Navier-Stokes system.
 *
 * \param[in]  previous    if true return the previous state otherwise the
 *                         current state.
 *
 * \return a pointer to an array of cs_real_t
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_navsto_get_mass_flux(bool   previous)
{
  cs_navsto_system_t  *ns = cs_navsto_system;

  if (ns == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_ns));

  cs_real_t  *mass_flux = ns->mass_flux_array;
  if (previous)
    mass_flux = ns->mass_flux_array_pre;

  return mass_flux;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes system
 *         At this stage, numerical settings should be completely determined
 *         but connectivity and geometrical information is not yet available.
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_init_setup(void)
{
  cs_navsto_system_t  *ns = cs_navsto_system;

  if (ns == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_ns));

  cs_navsto_param_t  *nsp = ns->param;

  /* Plotter output:
   * By default: there is at least the norm of the velocity field divergence */

  const char  **labels;
  int  n_plotter_outputs = 1;

  /* Set field metadata */

  const int  log_key = cs_field_key_id("log");
  const int  post_key = cs_field_key_id("post_vis");
  const bool  has_previous = cs_navsto_param_is_steady(nsp) ? false : true;
  int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE | CS_FIELD_CDO;

  /* Set the location id to define a mesh location support */

  int  location_id = -1;
  switch (nsp->space_scheme) {

  case CS_SPACE_SCHEME_CDOFB:
  case CS_SPACE_SCHEME_HHO_P0:
  case CS_SPACE_SCHEME_HHO_P1:
  case CS_SPACE_SCHEME_HHO_P2:
    location_id = cs_mesh_location_get_id_by_name("cells");
    break; /* Face-based scheme family */

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid space discretization scheme.", __func__);
  }

  /* Create if needed velocity and pressure fields */

  const int  field_post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;

  /* Handle the velocity field */

  ns->velocity = cs_field_find_or_create("velocity",
                                         field_mask,
                                         location_id,
                                         3, /* dimension */
                                         has_previous);

  /* Set the default value for keys related to log and post-processing */

  cs_field_set_key_int(ns->velocity, log_key, 1);
  cs_field_set_key_int(ns->velocity, post_key, field_post_flag);

  /* Handle the pressure field */

  ns->pressure = cs_field_find_or_create("pressure",
                                         field_mask,
                                         location_id,
                                         1, /* dimension */
                                         has_previous);

  /* Set the default value for keys related to log and post-processing */

  cs_field_set_key_int(ns->pressure, log_key, 1);
  cs_field_set_key_int(ns->pressure, post_key, field_post_flag);

  /* Handle the divergence of the velocity field.
   * Up to now, always defined the divergence of the velocity field. This
   * should be changed in the future */

  int  p_mask = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY | CS_FIELD_CDO;

  nsp->post_flag |= CS_NAVSTO_POST_VELOCITY_DIVERGENCE;

  /* Always post-process the velocity divergence */

  ns->velocity_divergence = cs_field_find_or_create("velocity_divergence",
                                                    p_mask,
                                                    location_id,
                                                    1, /* dimension */
                                                    has_previous);

  /* Set default value for keys related to log and post-processing */

  cs_field_set_key_int(ns->velocity_divergence, log_key, 1);
  cs_field_set_key_int(ns->velocity_divergence, post_key, field_post_flag);

  if (nsp->post_flag & CS_NAVSTO_POST_KINETIC_ENERGY) {

    n_plotter_outputs += 1;   /* Integral of the kinetic energy is monitored */

    ns->kinetic_energy = cs_field_find_or_create("kinetic_energy",
                                                 p_mask,
                                                 location_id,
                                                 1, /* dimension */
                                                 has_previous);

    /* Set the default value for keys related to log and post-processing */

    cs_field_set_key_int(ns->kinetic_energy, log_key, 1);
    cs_field_set_key_int(ns->kinetic_energy, post_key, field_post_flag);

  }

  if (nsp->post_flag & CS_NAVSTO_POST_MASS_DENSITY) {

    n_plotter_outputs += 1;   /* Integral of the mass density is monitored */

    ns->mass_density = cs_field_find_or_create("mass_density",
                                               p_mask,
                                               location_id,
                                               1, /* dimension */
                                               has_previous);

    /* Set the default value for keys related to log and post-processing */

    cs_field_set_key_int(ns->mass_density, log_key, 1);
    cs_field_set_key_int(ns->mass_density, post_key, field_post_flag);

  }

  if (nsp->post_flag & CS_NAVSTO_POST_CELL_MASS_FLUX_BALANCE) {

    ns->mass_flux_balance = cs_field_find_or_create("mass_flux_balance",
                                                    p_mask,
                                                    location_id,
                                                    1, /* dimension */
                                                    has_previous);

    /* Set the default value for keys related to log and post-processing */

    cs_field_set_key_int(ns->mass_flux_balance, log_key, 1);
    cs_field_set_key_int(ns->mass_flux_balance, post_key, field_post_flag);

  }

  if (nsp->post_flag & CS_NAVSTO_POST_PRESSURE_GRADIENT) {

    ns->pressure_gradient = cs_field_find_or_create("pressure_gradient",
                                                    p_mask,
                                                    location_id,
                                                    3, /* dimension */
                                                    has_previous);

    /* Set the default value for keys related to post-processing */

    cs_field_set_key_int(ns->pressure_gradient, post_key, field_post_flag);

  }

  if (nsp->post_flag & CS_NAVSTO_POST_STREAM_FUNCTION)
    nsp->post_flag |= CS_NAVSTO_POST_VORTICITY; /* automatic */

  if (nsp->post_flag & CS_NAVSTO_POST_HELICITY) {

    n_plotter_outputs += 1;   /* Integral of the helicity is monitored */

    nsp->post_flag |= CS_NAVSTO_POST_VORTICITY; /* automatic */
    ns->helicity = cs_field_find_or_create("helicity",
                                           p_mask,
                                           location_id,
                                           1, /* dimension */
                                           has_previous);

    /* Set default value for keys related to log and post-processing */

    cs_field_set_key_int(ns->helicity, log_key, 1);
    cs_field_set_key_int(ns->helicity, post_key, field_post_flag);

  }

  if (nsp->post_flag & CS_NAVSTO_POST_ENSTROPHY) {

    n_plotter_outputs += 1;   /* Integral of the enstrophy is monitored */

    nsp->post_flag |= CS_NAVSTO_POST_VORTICITY; /* automatic */
    ns->enstrophy = cs_field_find_or_create("enstrophy",
                                            p_mask,
                                            location_id,
                                            1, /* dimension */
                                            has_previous);

    /* Set default value for keys related to log and post-processing */

    cs_field_set_key_int(ns->enstrophy, log_key, 1);
    cs_field_set_key_int(ns->enstrophy, post_key, field_post_flag);

  }

  if (nsp->post_flag & CS_NAVSTO_POST_VORTICITY) {

    ns->vorticity = cs_field_find_or_create("vorticity",
                                            p_mask,
                                            location_id,
                                            3, /* dimension */
                                            has_previous);

    /* Set default value for keys related to log and post-processing */

    cs_field_set_key_int(ns->vorticity, log_key, 1);
    cs_field_set_key_int(ns->vorticity, post_key, field_post_flag);

  }

  if (nsp->post_flag & CS_NAVSTO_POST_VELOCITY_GRADIENT) {

    ns->velocity_gradient = cs_field_find_or_create("velocity_gradient",
                                                    p_mask,
                                                    location_id,
                                                    9, /* dimension */
                                                    has_previous);

    /* Set default value for keys related to log and post-processing */

    cs_field_set_key_int(ns->velocity_gradient, log_key, 1);
    cs_field_set_key_int(ns->velocity_gradient, post_key, field_post_flag);

  }

  /* Time plot monitor of global quantities predefined by the automatic
     post-processing */

  if (cs_glob_rank_id < 1) {

    assert(n_plotter_outputs > 0);
    BFT_MALLOC(labels, n_plotter_outputs, const char *);

    int  n_cols = 0;

    if (nsp->post_flag & CS_NAVSTO_POST_VELOCITY_DIVERGENCE)
      labels[n_cols++] = "vel_div_norm";
    if (nsp->post_flag & CS_NAVSTO_POST_MASS_DENSITY)
      labels[n_cols++] = "mass";
    if (nsp->post_flag & CS_NAVSTO_POST_KINETIC_ENERGY)
      labels[n_cols++] = "kinetic_energy";
    if (nsp->post_flag & CS_NAVSTO_POST_ENSTROPHY)
      labels[n_cols++] = "enstrophy";
    if (nsp->post_flag & CS_NAVSTO_POST_HELICITY)
      labels[n_cols++] = "helicity";

    ns->plot_writer = cs_time_plot_init_probe("navsto_monitor",
                                              "",
                                              CS_TIME_PLOT_DAT,
                                              false, /* use iteration */
                                              300,   /* flush time */
                                              -1,
                                              n_cols,
                                              NULL,
                                              NULL,
                                              labels);

    BFT_FREE(labels);

  } /* monitoring */

  /* Setup data according to the type of coupling (add the advection field in
     case of Navier-Stokes equations) */

  switch (nsp->coupling) {

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
    cs_navsto_ac_init_setup(nsp, ns->adv_field, ns->coupling_context);
    break;
  case CS_NAVSTO_COUPLING_MONOLITHIC:
    cs_navsto_monolithic_init_setup(nsp, ns->adv_field, ns->coupling_context);
    break;
  case CS_NAVSTO_COUPLING_PROJECTION:
    cs_navsto_projection_init_setup(nsp,
                                    ns->adv_field,
                                    location_id,
                                    has_previous,
                                    ns->coupling_context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _err_invalid_coupling, __func__);
    break;

  }

  /* Initialize the turbulence modelling */

  cs_equation_t *mom_eq = cs_navsto_system_get_momentum_eq();
  cs_turbulence_init_setup(ns->turbulence, mom_eq);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the settings for SLES related to the Navier-Stokes system
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_set_sles(void)
{
  cs_navsto_system_t  *ns = cs_navsto_system;

  if (ns == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_ns));

  cs_navsto_param_t *nsp = ns->param;

  switch (nsp->space_scheme) {

  case CS_SPACE_SCHEME_CDOFB:
  case CS_SPACE_SCHEME_HHO_P0:
    switch (nsp->coupling) {

    case CS_NAVSTO_COUPLING_MONOLITHIC:
      cs_cdofb_monolithic_set_sles(nsp, ns->scheme_context);
      break;

    case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
      cs_cdofb_ac_set_sles(nsp, ns->coupling_context);
      break;

    case CS_NAVSTO_COUPLING_PROJECTION:
      cs_cdofb_predco_set_sles(nsp, ns->coupling_context);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, _err_invalid_coupling, __func__);
      break;

    } /* Switch algo. for coupling velocity/pressure */
    break; /* Face-based scheme family */

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid space discretization scheme.", __func__);

  } /* Switch space scheme */

  if (nsp->post_flag & CS_NAVSTO_POST_STREAM_FUNCTION) {

    cs_equation_param_t  *eqp = cs_equation_get_param(ns->stream_function_eq);
    assert(eqp != NULL);

    /* Equation related to Navier-Stokes do not follow the classical setup
       stage */
    cs_equation_param_set_sles(eqp);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Last step of the setup of the Navier-Stokes system
 *
 * \param[in]  mesh       pointer to a cs_mesh_t structure
 * \param[in]  connect    pointer to a cs_cdo_connect_t structure
 * \param[in]  quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_step  pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_finalize_setup(const cs_mesh_t            *mesh,
                                const cs_cdo_connect_t     *connect,
                                const cs_cdo_quantities_t  *quant,
                                const cs_time_step_t       *time_step)
{
  cs_navsto_system_t  *ns = cs_navsto_system;

  if (ns == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_ns));

  cs_navsto_param_t  *nsp = ns->param;
  assert(connect != NULL && quant != NULL && nsp != NULL);

  /* Setup checkings */

  if ((nsp->model_flag & (CS_NAVSTO_MODEL_BOUSSINESQ |
                          CS_NAVSTO_MODEL_PASSIVE_THERMAL_TRACER)) > 0) {

    if (cs_thermal_system_is_activated() == false)
      bft_error(__FILE__, __LINE__, 0,
                " %s: The Navier-Stokes module is activated with options"
                " that imply a thermal module.\n"
                " Please check that cs_thermal_system_activate() has been"
                " called.\n", __func__);

    /* Check that an advection term has been set */

    cs_equation_param_t  *eqp = cs_equation_param_by_name(CS_THERMAL_EQNAME);
    if (cs_equation_param_has_convection(eqp) == false)
      bft_error(__FILE__, __LINE__, 0,
                " %s: No advection field is associated with the thermal"
                " equation\n whereas the Navier-Stokes is associated with"
                " the thermal equation.\n"
                " Please check your settings.", __func__);

  }

  /* Remaining boundary conditions:
   * 1. Walls
   * 2. Symmetries
   * 3. Outlets
   */

  cs_navsto_set_fixed_walls(nsp);
  cs_navsto_set_symmetries(nsp);
  cs_navsto_set_outlets(nsp);

  /* Define the advection field which relies on the mass flux. */

  switch (nsp->space_scheme) {

  case CS_SPACE_SCHEME_CDOFB:
  case CS_SPACE_SCHEME_HHO_P0:
    {
      BFT_MALLOC(ns->mass_flux_array, quant->n_faces, cs_real_t);
      memset(ns->mass_flux_array, 0, sizeof(cs_real_t)*quant->n_faces);

      BFT_MALLOC(ns->mass_flux_array_pre, quant->n_faces, cs_real_t);
      memset(ns->mass_flux_array_pre, 0, sizeof(cs_real_t)*quant->n_faces);

      cs_flag_t loc_flag =
        CS_FLAG_FULL_LOC | CS_FLAG_SCALAR | cs_flag_primal_face;

      cs_advection_field_def_by_array(ns->adv_field,
                                      loc_flag,
                                      ns->mass_flux_array,
                                      false, /* advection field is not owner */
                                      NULL, NULL); /* no index/ids */
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid space discretization scheme.", __func__);

  } /* End of switch on space_scheme */

  /* Last setup stage according to the type of coupling (not related to
     space discretization scheme */

  switch (nsp->coupling) {

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
    cs_navsto_ac_last_setup(nsp, ns->coupling_context);
    break;
  case CS_NAVSTO_COUPLING_MONOLITHIC:
    cs_navsto_monolithic_last_setup(nsp, ns->coupling_context);
    break;
  case CS_NAVSTO_COUPLING_PROJECTION:
    cs_navsto_projection_last_setup(quant, nsp, ns->coupling_context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _err_invalid_coupling, __func__);
    break;

  }

  /* Settings with respect to the discretization scheme */

  switch (nsp->space_scheme) {

  case CS_SPACE_SCHEME_CDOFB:
  case CS_SPACE_SCHEME_HHO_P0:

    /* Setup functions and data according to the type of coupling */

    switch (nsp->coupling) {

    case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
      /* ============================================= */

      ns->init_scheme_context = cs_cdofb_ac_init_scheme_context;
      ns->free_scheme_context = cs_cdofb_ac_free_scheme_context;
      ns->init_velocity = NULL;
      ns->init_pressure = cs_cdofb_navsto_init_pressure;
      ns->compute_steady = NULL;

      if (_handle_non_linearities(nsp)) {

        /* Same function. The difference will be taken care of by mean of a
         * build sub-function */

        if (nsp->time_scheme == CS_TIME_SCHEME_EULER_IMPLICIT)
          ns->compute = cs_cdofb_ac_compute_implicit_nl;
        else
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Invalid time scheme for the AC and Navier-Stokes",
                    __func__);
      }
      else {

        switch (nsp->time_scheme) {

        case CS_TIME_SCHEME_EULER_IMPLICIT:
          ns->compute = cs_cdofb_ac_compute_implicit;
          break;

        case CS_TIME_SCHEME_THETA:
        case CS_TIME_SCHEME_CRANKNICO:
        case CS_TIME_SCHEME_BDF2:
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Time scheme not implemented for the AC coupling",
                    __func__);
          break;

        case CS_TIME_SCHEME_STEADY:
        default:
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Invalid time scheme for the AC coupling", __func__);
          break;

        } /* Switch */

      } /* If nonlinear */

      cs_cdofb_ac_init_common(quant, connect, time_step);
      break;

    case CS_NAVSTO_COUPLING_MONOLITHIC:
      /* ============================= */

      ns->init_scheme_context = cs_cdofb_monolithic_init_scheme_context;
      ns->free_scheme_context = cs_cdofb_monolithic_free_scheme_context;
      ns->init_velocity = NULL;
      ns->init_pressure = cs_cdofb_navsto_init_pressure;
      if (_handle_non_linearities(nsp))
        ns->compute_steady = cs_cdofb_monolithic_steady_nl;
      else
        ns->compute_steady = cs_cdofb_monolithic_steady;

      switch (nsp->time_scheme) {

      case CS_TIME_SCHEME_STEADY:
        if (_handle_non_linearities(nsp))
          ns->compute = cs_cdofb_monolithic_steady_nl;
        else
          ns->compute = cs_cdofb_monolithic_steady;
        break; /* Nothing to set */

      case CS_TIME_SCHEME_EULER_IMPLICIT:
      case CS_TIME_SCHEME_THETA:
      case CS_TIME_SCHEME_CRANKNICO:
        if (_handle_non_linearities(nsp))
          ns->compute = cs_cdofb_monolithic_nl;
        else
          ns->compute = cs_cdofb_monolithic;
        break;

      case CS_TIME_SCHEME_BDF2:
      default:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid time scheme for the monolithic coupling",
                  __func__);
        break;

      } /* Switch */

      cs_cdofb_monolithic_init_sharing(nsp, mesh, quant, connect, time_step);
      break;

    case CS_NAVSTO_COUPLING_PROJECTION:
      /* ============================= */

      ns->init_scheme_context = cs_cdofb_predco_init_scheme_context;
      ns->free_scheme_context = cs_cdofb_predco_free_scheme_context;
      ns->init_velocity = NULL;
      ns->init_pressure = cs_cdofb_navsto_init_pressure;
      ns->compute_steady = NULL;

      switch (nsp->time_scheme) {

      case CS_TIME_SCHEME_STEADY:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: The projection coupling algorithm "
                  "can be used only in unsteady problems", __func__);
        break;

      case CS_TIME_SCHEME_EULER_IMPLICIT:
        ns->compute = cs_cdofb_predco_compute_implicit;
        break;

      case CS_TIME_SCHEME_THETA:
      case CS_TIME_SCHEME_CRANKNICO:
      case CS_TIME_SCHEME_BDF2:
      default:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid time scheme for the projection coupling"
                  " algorithm", __func__);
        break;

      } /* Switch */

      cs_cdofb_predco_init_common(quant, connect, time_step);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, _err_invalid_coupling, __func__);
      break;

    }
    break; /* Lowest-order face-based schemes */

  case CS_SPACE_SCHEME_HHO_P1:
  case CS_SPACE_SCHEME_HHO_P2:
    /* TODO: set function pointers */
    break; /* HHO schemes */

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid space discretization scheme.", __func__);
  }

  if (fabs(nsp->reference_pressure) > 0 && nsp->n_pressure_ic_defs == 0) {

    /* Initialize the initial pressure to the reference pressure */

    cs_navsto_add_pressure_ic_by_value(nsp, NULL, &(nsp->reference_pressure));

  }

  /* Add default post-processing related to the Navier-Stokes system */

  cs_post_add_time_mesh_dep_output(cs_navsto_system_extra_post, ns);

  if (nsp->post_flag & CS_NAVSTO_POST_STREAM_FUNCTION) {

    cs_equation_param_t  *eqp = cs_equation_get_param(ns->stream_function_eq);
    assert(eqp != NULL);
    cs_field_t  *w = cs_field_by_name("vorticity");

    /* Add a laplacian term: -div.grad */

    cs_equation_add_diffusion(eqp, cs_property_by_name("unity"));

    /* Add source term as the vorticity w.r.t. the z-axis */

    cs_equation_add_source_term_by_dof_func(eqp,
                                            NULL,
                                            cs_flag_primal_cell,
                                            cs_cdofb_navsto_stream_source_term,
                                            (void *)w->val);

  } /* Post-processing of the stream function is requested */

  /* Finalize the setting of the turbulence modelling */

  cs_turbulence_finalize_setup(mesh, connect, quant, time_step,
                               ns->turbulence);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the scheme context structure used to build the algebraic
 *         system. This is done after the setup step.
 *
 * \param[in]  mesh        pointer to a cs_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_define_context(const cs_mesh_t             *mesh)
{
  cs_navsto_system_t  *ns = cs_navsto_system;

  if (ns == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_ns));

  const cs_navsto_param_t *nsp = ns->param;
  assert(nsp != NULL);
  if (nsp->space_scheme != CS_SPACE_SCHEME_CDOFB)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid space discretization scheme.", __func__);

  /* Allocate then define an array of boundary types for each boundary face */

  BFT_MALLOC(ns->bf_type, mesh->n_b_faces, cs_boundary_type_t);
  cs_boundary_build_type_array(nsp->boundaries, mesh->n_b_faces, ns->bf_type);

  /* Allocate and initialize the scheme context structure */

  ns->scheme_context = ns->init_scheme_context(nsp,
                                               ns->adv_field,
                                               ns->mass_flux_array,
                                               ns->mass_flux_array_pre,
                                               ns->bf_type,
                                               ns->coupling_context);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set an initial value for the velocity and pressure fields as well
 *         as mass fluxes and tubulent quantities if needed
 *
 * \param[in]  mesh        pointer to a cs_mesh_t structure
 * \param[in]  connect     pointer to a cs_cdo_connect_t structure
 * \param[in]  quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_step   pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_init_values(const cs_mesh_t             *mesh,
                             const cs_cdo_connect_t      *connect,
                             const cs_cdo_quantities_t   *quant,
                             const cs_time_step_t        *time_step)
{
  cs_navsto_system_t  *ns = cs_navsto_system;

  if (ns == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_ns));

  const cs_navsto_param_t *nsp = ns->param;
  assert(nsp != NULL);
  if (nsp->space_scheme != CS_SPACE_SCHEME_CDOFB)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid space discretization scheme.", __func__);

  if (time_step->nt_cur < 1) {

    /* Initial conditions for the velocity */

    if (ns->init_velocity != NULL)
      ns->init_velocity(nsp, quant, time_step, ns->scheme_context);

    /* Initial conditions for the pressure */

    if (ns->init_pressure != NULL)
      ns->init_pressure(nsp, quant, time_step, ns->pressure);

  }

  if (nsp->space_scheme == CS_SPACE_SCHEME_CDOFB) {

    if (nsp->coupling == CS_NAVSTO_COUPLING_PROJECTION) {

      /* The call to the initialization of the cell pressure should be done
         before */

      cs_real_t  *pr_f = cs_cdofb_predco_get_face_pressure(ns->scheme_context);

      cs_cdofb_navsto_init_face_pressure(nsp, connect, time_step, pr_f);

    }

    /* Initialize the mass flux values */

    const cs_equation_t  *mom_eq = cs_navsto_system_get_momentum_eq();
    const cs_real_t  *face_vel = cs_equation_get_face_values(mom_eq, false);
    cs_cdofb_navsto_mass_flux(nsp, quant, face_vel, ns->mass_flux_array);

  } /* Face-based schemes */

  /* Initialize the values of the tubulent quantities */

  cs_turbulence_init_values(mesh, connect, quant, time_step, ns->turbulence);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a solid zone related to the Navier-Stokes equations
 *
 * \param[in] n_solid_cells    number of solid cells
 * \param[in] solid_cell_ids   list of cell ids (local numbering)
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_set_solid_cells(cs_lnum_t          n_solid_cells,
                                 cs_lnum_t          solid_cell_ids[])
{
  cs_navsto_system_t  *ns = cs_navsto_system;

  if (ns == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_ns));

  /* Do not exit the function if n_elts = 0 since there is a parallel
     synchronization to count the global number of cells to enforce in the next
     call */

  /* The momentum equation has to enforce a zero-velocity */

  cs_equation_t  *mom_eq = cs_navsto_system_get_momentum_eq();
  cs_equation_param_t  *mom_eqp = cs_equation_get_param(mom_eq);
  cs_real_t  zero_velocity[3] = {0, 0, 0};

  if (cs_navsto_system_solid_enforcement_id < 0)
    cs_navsto_system_solid_enforcement_id = mom_eqp->n_enforcements;
  cs_equation_add_or_replace_cell_enforcement(mom_eqp,
                      cs_navsto_system_solid_enforcement_id,
                                              n_solid_cells,
                                              solid_cell_ids,
                                              zero_velocity,
                                              NULL);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update variables and related quantities when a new state of the
 *         Navier-Stokes system has been computed
 *
 * \param[in] mesh       pointer to a cs_mesh_t structure
 * \param[in] connect    pointer to a cs_cdo_connect_t structure
 * \param[in] quant      pointer to a cs_cdo_quantities_t structure
 * \param[in] time_step  structure managing the time stepping
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_update(const cs_mesh_t             *mesh,
                        const cs_cdo_connect_t      *connect,
                        const cs_cdo_quantities_t   *quant,
                        const cs_time_step_t        *time_step)
{
  CS_UNUSED(mesh);

  cs_navsto_system_t  *ns = cs_navsto_system;

  if (ns == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_ns));

  /* Synchronize the cell pressure before doing something */

  cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, ns->pressure->val);

  /* Update quantities relative to turbulence (e.g. turbulence
   * viscosity...) */

  cs_turbulence_t  *tbs = ns->turbulence;
  if (tbs->update != NULL)
    tbs->update(mesh,
                connect,
                quant,
                time_step,
                tbs);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build, solve and update the Navier-Stokes system in case of a
 *         steady-state approach
 *
 * \param[in] mesh       pointer to a cs_mesh_t structure
 * \param[in] connect    pointer to a cs_cdo_connect_t structure
 * \param[in] quant      pointer to a cs_cdo_quantities_t structure
 * \param[in] time_step  structure managing the time stepping
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_compute_steady_state(const cs_mesh_t             *mesh,
                                      const cs_cdo_connect_t      *connect,
                                      const cs_cdo_quantities_t   *quant,
                                      const cs_time_step_t        *time_step)
{
  cs_navsto_system_t  *ns = cs_navsto_system;

  if (ns == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_ns));

  cs_navsto_param_t  *nsp = ns->param;
  if (!cs_navsto_param_is_steady(nsp))
    return;

  cs_turbulence_t  *tbs = ns->turbulence;

  /* First resolve the thermal system if needed */

  cs_equation_t  *th_eq = cs_thermal_system_get_equation();

  if (nsp->model_flag & CS_NAVSTO_MODEL_PASSIVE_THERMAL_TRACER) {

    assert(th_eq != NULL);

    /* Update variable, properties according to the computed variables */

    cs_navsto_system_update(mesh, connect, quant, time_step);

    /* Build and solve the Navier-Stokes system */

    ns->compute_steady(mesh, nsp, ns->scheme_context);

    /* Build and solve the turbulence variable system */

    if (tbs->compute_steady != NULL)
      tbs->compute_steady(mesh,
                          connect,
                          quant,
                          time_step,
                          tbs);

    /* Solve the thermal equation */

    if (cs_equation_param_has_time(cs_equation_get_param(th_eq)) == false)
      cs_thermal_system_compute_steady_state(mesh, connect, quant, time_step);

  }
  else if (cs_flag_test(nsp->model_flag, CS_NAVSTO_MODEL_BOUSSINESQ) &&
           cs_flag_test(nsp->model_flag,
                        CS_NAVSTO_MODEL_WITH_SOLIDIFICATION) == false) {

    /* Remark: The "solidification" case is handled in a dedicated module */

    assert(th_eq != NULL);
    if (cs_equation_param_has_time(cs_equation_get_param(th_eq)))
      bft_error(__FILE__, __LINE__, 0,
                " %s: Steady-state computation for Navier-Stokes with a"
                " Boussinesq approximation\n"
                " whereas an unsteady thermal equation is set.\n"
                " Please check your settings.", __func__);

    cs_real_t  *th_var = cs_equation_get_cell_values(th_eq, false);

    cs_real_t  *th_var_iter_prev = NULL;

    BFT_MALLOC(th_var_iter_prev, quant->n_cells, cs_real_t);
    memcpy(th_var_iter_prev, th_var, quant->n_cells*sizeof(cs_real_t));

    cs_real_t  inv_tref = cs_thermal_system_get_reference_temperature();

    if (fabs(inv_tref) > cs_math_zero_threshold)
      inv_tref = 1./inv_tref;
    else
      inv_tref = 1.;

    cs_real_t  delta_th_tolerance = 1 + nsp->delta_thermal_tolerance;
    int  iter = 0;

    do {

      /* Update variable, properties according to the computed variables */

      cs_navsto_system_update(mesh, connect, quant, time_step);

      /* Build and solve the thermal system */

      cs_thermal_system_compute_steady_state(mesh, connect, quant, time_step);

      /* Build and solve the Navier-Stokes system */

      ns->compute_steady(mesh, nsp, ns->scheme_context);

      /* Build and solve the turbulence variable system */

      if (tbs->compute_steady != NULL)
        tbs->compute_steady(mesh,
                            connect,
                            quant,
                            time_step,
                            tbs);

      /* Check convergence */

      cs_real_t  delta = -1;
      for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {
        cs_real_t  _delta = fabs(th_var[c_id] - th_var_iter_prev[c_id]);
        th_var_iter_prev[c_id] = th_var[c_id];
        if (_delta > delta)
          delta = _delta;
      }

      /* Update the stopping criteria */

      delta_th_tolerance = delta * inv_tref;
      iter++;

      if (nsp->verbosity > 0)
        cs_log_printf(CS_LOG_DEFAULT,
                      "### Boussinesq.Iteration: %2d | delta_th_var= %5.3e\n",
                      iter, delta_th_tolerance);

    }
    while ( (delta_th_tolerance > nsp->delta_thermal_tolerance) &&
            (iter < nsp->n_max_outer_iter) );

    cs_log_printf(CS_LOG_DEFAULT, " Steady algorithm exits.\n"
                  " Number of iterations: %d\n"
                  " Tolerance on the thermal variable: %5.3e\n",
                  iter, delta_th_tolerance);

    BFT_FREE(th_var_iter_prev);

  }
  else {

    /* Build and solve the Navier-Stokes system */

    ns->compute_steady(mesh, nsp, ns->scheme_context);

    /* Build and solve the turbulence variable system
     * TODO, we should loop over NS+turbulence and add a
     * stopping criteria on the NS matrix for instance
     * */

    if (tbs->compute_steady != NULL)
      tbs->compute_steady(mesh,
                          connect,
                          quant,
                          time_step,
                          tbs);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build, solve and update the Navier-Stokes system
 *
 * \param[in] mesh       pointer to a cs_mesh_t structure
 * \param[in] connect    pointer to a cs_cdo_connect_t structure
 * \param[in] quant       pointer to a cs_cdo_quantities_t structure
 * \param[in] time_step  structure managing the time stepping
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_compute(const cs_mesh_t             *mesh,
                         const cs_cdo_connect_t      *connect,
                         const cs_cdo_quantities_t   *quant,
                         const cs_time_step_t        *time_step)
{
  cs_navsto_system_t  *ns = cs_navsto_system;

  if (ns == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_ns));

  const cs_navsto_param_t  *nsp = ns->param;

  cs_turbulence_t  *tbs = ns->turbulence;

  /* Update variable, properties according to the new computed variables */

  cs_navsto_system_update(mesh, connect, quant, time_step);

  if (nsp->model_flag & CS_NAVSTO_MODEL_PASSIVE_THERMAL_TRACER) {

    /* First resolve the thermal system if needed */

    cs_equation_t  *th_eq = cs_thermal_system_get_equation();

    assert(th_eq != NULL);

    /* Build and solve the Navier-Stokes system */

    if (cs_navsto_param_is_steady(nsp) == false) {
      ns->compute(mesh, nsp, ns->scheme_context);

      /* Build and solve the turbulence variable system */

      if (tbs->compute != NULL)
        tbs->compute(mesh,
                     connect,
                     quant,
                     time_step,
                     tbs);

    }

    /* Solve the thermal equation */

    if (cs_equation_param_has_time(cs_equation_get_param(th_eq)))
      cs_thermal_system_compute(true, mesh, connect, quant, time_step);

  }
  else if (cs_flag_test(nsp->model_flag, CS_NAVSTO_MODEL_BOUSSINESQ) &&
           cs_flag_test(nsp->model_flag,
                        CS_NAVSTO_MODEL_WITH_SOLIDIFICATION) == false) {

    /* Remark: The "solidification" case is handled in a dedicated module */

    if (cs_navsto_param_is_steady(nsp))
      return;

    /* First resolve the thermal system if needed */

    cs_equation_t  *th_eq = cs_thermal_system_get_equation();

    assert(th_eq != NULL);

    if (cs_equation_param_has_time(cs_equation_get_param(th_eq))) {

      /* Build and solve the thermal system */

      cs_thermal_system_compute(true, mesh, connect, quant, time_step);

      /* Build and solve the Navier-Stokes system */

      ns->compute(mesh, nsp, ns->scheme_context);

      /* Build and solve the turbulence variable system */

      if (tbs->compute != NULL)
        tbs->compute(mesh,
                     connect,
                     quant,
                     time_step,
                     tbs);


    }
    else { /* Thermal system is declared as steady. So, this is equivalent to a
              standard Navier-Stokes iteration */

      /* Build and solve the Navier-Stokes system */

      ns->compute(mesh, nsp, ns->scheme_context);

    }

  }
  else {

    if (cs_navsto_param_is_steady(nsp))
      return;

    /* Build and solve the Navier-Stokes system */

    ns->compute(mesh, nsp, ns->scheme_context);

    /* Build and solve the turbulence variable system */

    if (tbs->compute != NULL)
      tbs->compute(mesh,
                   connect,
                   quant,
                   time_step,
                   tbs);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations for the Navier-Stokes system
 *
 * \param[in]  mesh        pointer to a cs_mesh_t structure
 * \param[in]  connect     pointer to a cs_cdo_connect_t structure
 * \param[in]  quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]  time_step   pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_extra_op(const cs_mesh_t             *mesh,
                          const cs_cdo_connect_t      *connect,
                          const cs_cdo_quantities_t   *quant,
                          const cs_time_step_t        *time_step)
{
  cs_navsto_system_t  *navsto = cs_navsto_system;

  if (navsto == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_ns));

  const cs_navsto_param_t  *nsp = navsto->param;

  switch (nsp->space_scheme) {

  case CS_SPACE_SCHEME_CDOFB:
    {
      /* Get the current values not the previous one */

      bool  need_prev = false;

      /* Mass flux is a scalar associated to each face (first interior faces
         then border faces */

      const cs_real_t  *mass_flux = cs_navsto_get_mass_flux(need_prev);

      const cs_adv_field_t  *adv = cs_navsto_get_adv_field();
      const cs_equation_t  *eq = cs_navsto_system_get_momentum_eq();
      const cs_real_t  *u_face = cs_equation_get_face_values(eq, need_prev);
      const cs_real_t  *u_cell = navsto->velocity->val;
      cs_real_t  *p_cell = navsto->pressure->val;

      /* Synchronize the cell pressure before doing something */

      if (nsp->post_flag & CS_NAVSTO_POST_PRESSURE_GRADIENT)
        cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, p_cell);

      cs_cdofb_navsto_extra_op(nsp, mesh, quant, connect, time_step,
                               navsto->plot_writer,
                               adv, mass_flux, p_cell, u_cell, u_face);
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid space discretization scheme.", __func__);
    break;

  } /* End of switch */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for the Navier-Stokes system.
 *         The prototype of this function is fixed since it is a function
 *         pointer defined in cs_post.h (\ref cs_post_time_mesh_dep_output_t)
 *
 * \param[in, out] input        pointer to a optional structure (here a
 *                              cs_navsto_system_t structure)
 * \param[in]      mesh_id      id of the output mesh for the current call
 * \param[in]      cat_id       category id of the output mesh for this call
 * \param[in]      ent_flag     indicate global presence of cells (ent_flag[0]),
 *                              interior faces (ent_flag[1]), boundary faces
 *                              (ent_flag[2]), particles (ent_flag[3]) or probes
 *                              (ent_flag[4])
 * \param[in]      n_cells      local number of cells of post_mesh
 * \param[in]      n_i_faces    local number of interior faces of post_mesh
 * \param[in]      n_b_faces    local number of boundary faces of post_mesh
 * \param[in]      cell_ids     list of cells (0 to n-1)
 * \param[in]      i_face_ids   list of interior faces (0 to n-1)
 * \param[in]      b_face_ids   list of boundary faces (0 to n-1)
 * \param[in]      time_step    pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_extra_post(void                      *input,
                            int                        mesh_id,
                            int                        cat_id,
                            int                        ent_flag[5],
                            cs_lnum_t                  n_cells,
                            cs_lnum_t                  n_i_faces,
                            cs_lnum_t                  n_b_faces,
                            const cs_lnum_t            cell_ids[],
                            const cs_lnum_t            i_face_ids[],
                            const cs_lnum_t            b_face_ids[],
                            const cs_time_step_t      *time_step)
{
  CS_UNUSED(n_cells);
  CS_UNUSED(n_i_faces);
  CS_UNUSED(n_b_faces);
  CS_UNUSED(cell_ids);
  CS_UNUSED(i_face_ids);
  CS_UNUSED(b_face_ids);

  cs_navsto_system_t  *ns = (cs_navsto_system_t *)input;
  if (ns == NULL)
    return;

  const cs_navsto_param_t  *nsp = ns->param;

  /* Post-processing of the boundary mass flux */
  if (cat_id == CS_POST_MESH_BOUNDARY && ent_flag[2] > 0) {

    switch (nsp->space_scheme) {

    case CS_SPACE_SCHEME_CDOFB:
    case CS_SPACE_SCHEME_HHO_P0:
      {
        /* Get the current values not the previous one */

        bool  need_prev = false;

        /* In case of postprocessing of the border faces, one has to check if
           there is a mesh modification. In particular, a removal of 2D
           extruded border faces*/

        bool  use_parent =
          (cs_glob_mesh->n_g_b_faces_all > cs_glob_mesh->n_g_b_faces) ?
          false : true;

        /* Mass flux is a scalar associated to each face (first interior faces
           then border faces */
        const cs_real_t  *mass_flux = cs_navsto_get_mass_flux(need_prev);
        const cs_real_t  *mass_bflux = mass_flux + cs_glob_mesh->n_i_faces;

        cs_post_write_var(mesh_id,
                          CS_POST_WRITER_DEFAULT,
                          "boundary_mass_flux",
                          1,
                          false,                  /* interlace */
                          use_parent,             /* true = original mesh */
                          CS_POST_TYPE_cs_real_t,
                          NULL,                   /* values on cells */
                          NULL,                   /* values at internal faces */
                          mass_bflux,             /* values at border faces */
                          time_step);             /* time step structure */

      }
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                "%s: Invalid space scheme\n", __func__);
      break;

    }

  } /* Boundary mesh type */

  if (cat_id == CS_POST_MESH_VOLUME && ent_flag[0] > 0) {
    /* Velocity-Pressure coupling algorithm */
    switch (nsp->coupling) {

    case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
    case CS_NAVSTO_COUPLING_MONOLITHIC:
      /* Nothing to do up to now */
      break;

    case CS_NAVSTO_COUPLING_PROJECTION:
      {
        cs_navsto_projection_t  *cc
          = (cs_navsto_projection_t *)ns->coupling_context;

        const cs_field_t  *velp = cc->predicted_velocity;

        /* Post-process the predicted velocity */
        cs_post_write_var(mesh_id,
                          CS_POST_WRITER_DEFAULT,
                          velp->name,
                          3,
                          true,           // interlace
                          true,           // true = original mesh
                          CS_POST_TYPE_cs_real_t,
                          velp->val,      // values on cells
                          NULL,           // values at internal faces
                          NULL,           // values at border faces
                          time_step);     // time step management struct.

        /* Post-process the source term of the correction equation on the
           pressure increment (-div(velp_f) */
        cs_post_write_var(mesh_id,
                          CS_POST_WRITER_DEFAULT,
                          "-DivVelPred",
                          1,
                          true,       // interlace
                          true,       // true = original mesh
                          CS_POST_TYPE_cs_real_t,
                          cc->div_st, // values on cells
                          NULL,       // values at internal faces
                          NULL,       // values at border faces
                          time_step); // time step management struct.
      }
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, _err_invalid_coupling, __func__);
      break;
    }

  } /* cat_id == CS_POST_MESH_VOLUME */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of the main cs_navsto_system_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_log_setup(void)
{
  cs_navsto_system_t  *ns = cs_navsto_system;

  if (ns == NULL)
    return;

  cs_log_printf(CS_LOG_SETUP, "\n");
  cs_log_printf(CS_LOG_SETUP, "%s", cs_sep_h1);
  cs_log_printf(CS_LOG_SETUP, "\tSummary of the Navier-Stokes system\n");
  cs_log_printf(CS_LOG_SETUP, "%s", cs_sep_h1);

  /* Main set of numerical parameters */
  cs_navsto_param_log(ns->param);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
