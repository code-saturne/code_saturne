/*============================================================================
 * Routines to handle cs_navsto_system_t structure
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
#include "cs_cdofb_navsto.h"
#include "cs_cdofb_uzawa.h"
#include "cs_hho_stokes.h"
#include "cs_equation.h"
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

/*!
 *  \file cs_navsto_system.c
 *
 *  \brief  Routines to handle the cs_navsto_system_t structure
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

static cs_navsto_system_t  *cs_navsto_system = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

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

  /* Velocity in the case of Navier-Stokes or Stokes,
     Wind or advection field in the case of the Oseen problem */
  navsto->adv_field = NULL;

  /* Main set of variables */
  navsto->velocity = NULL;
  navsto->velocity_divergence = NULL;
  navsto->pressure = NULL;
  navsto->temperature = NULL;

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
 * \brief  Allocate and initialize the Navier-Stokes (NS) system
 *
 * \param[in] boundaries     pointer to the domain boundaries
 * \param[in] model          type of model related to the NS system
 * \param[in] time_state     state of the time for the NS equations
 * \param[in] algo_coupling  algorithm used for solving the NS system
 *
 * \return a pointer to a new allocated cs_navsto_system_t structure
 */
/*----------------------------------------------------------------------------*/

cs_navsto_system_t *
cs_navsto_system_activate(const cs_boundary_t           *boundaries,
                          cs_navsto_param_model_t        model,
                          cs_navsto_param_time_state_t   time_state,
                          cs_navsto_param_coupling_t     algo_coupling)
{
  /* Sanity checks */
  if (model == CS_NAVSTO_N_MODELS)
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid model for Navier-Stokes.\n",
              __func__);

  /* Allocate an empty structure */
  cs_navsto_system_t  *navsto = _allocate_navsto_system();

  /* Initialize the set of parameters */
  navsto->param
    = cs_navsto_param_create(boundaries, model, time_state, algo_coupling);

  /* Advection field related to the resolved velocity */
  navsto->adv_field = cs_advection_field_add("velocity_field",
                                             CS_ADVECTION_FIELD_NAVSTO);

  /* Normal flux at the boundary faces */
  cs_advection_field_set_option(navsto->adv_field,
                                CS_ADVKEY_DEFINE_AT_BOUNDARY_FACES);

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

  } /* Switch */

  /* Additional initialization fitting the choice of model */
  switch (navsto->param->coupling) {

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
    navsto->coupling_context = cs_navsto_ac_create_context(navsto->param,
                                                           default_bc);
    break;
  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY_VPP:
    navsto->coupling_context = cs_navsto_ac_vpp_create_context(navsto->param,
                                                               default_bc);
    break;
  case CS_NAVSTO_COUPLING_MONOLITHIC:
    navsto->coupling_context
      = cs_navsto_monolithic_create_context(navsto->param, default_bc);
    break;
  case CS_NAVSTO_COUPLING_PROJECTION:
    navsto->coupling_context =
      cs_navsto_projection_create_context(navsto->param, default_bc);
    break;
  case CS_NAVSTO_COUPLING_UZAWA:
    navsto->coupling_context = cs_navsto_uzawa_create_context(navsto->param,
                                                              default_bc);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _err_invalid_coupling, __func__);
    break;

  }

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
    Properties, advection fields, equations and fields are all destroyed
    respectively inside cs_property_destroy_all(),
    cs_advection_field_destroy_all(), cs_equation_destroy_all() and
    cs_field_destroy_all()
  */

  cs_navsto_param_t  *nsp = navsto->param;

  /* Free the context according to the model choice */
  switch (nsp->coupling) {

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
    navsto->coupling_context =
      cs_navsto_ac_free_context(nsp, navsto->coupling_context);
    break;
  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY_VPP:
    navsto->coupling_context =
      cs_navsto_ac_vpp_free_context(nsp, navsto->coupling_context);
    break;
  case CS_NAVSTO_COUPLING_MONOLITHIC:
    navsto->coupling_context =
      cs_navsto_monolithic_free_context(nsp, navsto->coupling_context);
    break;
  case CS_NAVSTO_COUPLING_PROJECTION:
    navsto->coupling_context =
      cs_navsto_projection_free_context(nsp, navsto->coupling_context);
    break;
  case CS_NAVSTO_COUPLING_UZAWA:
    navsto->coupling_context =
      cs_navsto_uzawa_free_context(nsp, navsto->coupling_context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _err_invalid_coupling, __func__);
    break;
  }

  /* Destroy the context related to the discretization scheme */
  navsto->free_scheme_context(navsto->scheme_context);

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

  /* Set field metadata */
  const bool  has_previous = cs_navsto_param_is_steady(nsp) ? false : true;
  int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE;

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
  const int  post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;

  /* Handle the velocity field */
  ns->velocity = cs_field_find_or_create("velocity",
                                         field_mask,
                                         location_id,
                                         3, /* dimension */
                                         has_previous);

  /* Set default value for keys related to log and post-processing */
  cs_field_set_key_int(ns->velocity, cs_field_key_id("log"), 1);
  cs_field_set_key_int(ns->velocity, cs_field_key_id("post_vis"), post_flag);

  /* Handle the pressure field */
  ns->pressure = cs_field_find_or_create("pressure",
                                         field_mask,
                                         location_id,
                                         1, /* dimension */
                                         has_previous);

  /* Set default value for keys related to log and post-processing */
  cs_field_set_key_int(ns->pressure, cs_field_key_id("log"), 1);
  cs_field_set_key_int(ns->pressure, cs_field_key_id("post_vis"), post_flag);

  /* Handle the divergence of the velocity field */
  ns->velocity_divergence = cs_field_find_or_create("velocity_divergence",
                                                    CS_FIELD_INTENSIVE,
                                                    location_id,
                                                    1, /* dimension */
                                                    has_previous);

  /* Set default value for keys related to log and post-processing */
  cs_field_set_key_int(ns->velocity_divergence, cs_field_key_id("log"), 1);
  cs_field_set_key_int(ns->velocity_divergence, cs_field_key_id("post_vis"),
                       post_flag);

  /* TODO: temperature for the energy equation */

  /* Setup data according to the type of coupling */
  switch (nsp->coupling) {

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
    cs_navsto_ac_init_setup(nsp, ns->coupling_context);
    break;
  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY_VPP:
    cs_navsto_ac_vpp_init_setup(nsp, ns->coupling_context);
    break;
  case CS_NAVSTO_COUPLING_MONOLITHIC:
    cs_navsto_monolithic_init_setup(nsp, ns->coupling_context);
    break;
  case CS_NAVSTO_COUPLING_PROJECTION:
    cs_navsto_projection_init_setup(nsp, ns->coupling_context);
    break;
  case CS_NAVSTO_COUPLING_UZAWA:
    cs_navsto_uzawa_init_setup(nsp, ns->coupling_context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _err_invalid_coupling, __func__);
    break;

  }

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
  void  *nscc = ns->coupling_context;

  const cs_navsto_param_t *nsp = ns->param;

  switch (nsp->space_scheme) {

  case CS_SPACE_SCHEME_CDOFB:
  case CS_SPACE_SCHEME_HHO_P0:
    switch (nsp->coupling) {

    case CS_NAVSTO_COUPLING_MONOLITHIC:
      cs_cdofb_monolithic_set_sles(nsp, nscc);
      break;

    case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
      cs_cdofb_ac_set_sles(nsp, nscc);
      break;

    case CS_NAVSTO_COUPLING_UZAWA:
      cs_cdofb_uzawa_set_sles(nsp, nscc);
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

  assert(connect != NULL && quant != NULL);
  if (ns == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_ns));

  cs_navsto_param_t  *nsp = ns->param;

  /* Avoid an error if no definition is given for the mandatory physical
     properties */
  cs_real_t  one = 1.0;
  if (nsp->density->n_definitions == 0) /* Not set by the user */
    cs_property_def_iso_by_value(nsp->density,
                                 NULL, /* all cells */
                                 one);

  if (nsp->lami_viscosity->n_definitions == 0) /* Not set by the user */
    cs_property_def_iso_by_value(nsp->lami_viscosity,
                                 NULL, /* all cells */
                                 one);

  /* Remaining boundary conditions:
   * 1. Walls
   * 2. Symmetries
   * 3. Outlets
   */
  cs_navsto_set_fixed_walls(nsp);
  cs_navsto_set_symmetries(nsp);
  cs_navsto_set_outlets(nsp);

  /* Last setup stage according to the type of coupling (not related to
     space discretization scheme */
  switch (nsp->coupling) {

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
    cs_navsto_ac_last_setup(connect, quant, nsp, ns->coupling_context);
    break;
  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY_VPP:
    cs_navsto_ac_vpp_last_setup(connect, quant, nsp, ns->coupling_context);
    break;
  case CS_NAVSTO_COUPLING_MONOLITHIC:
    cs_navsto_monolithic_last_setup(connect, quant, nsp, ns->coupling_context);
    break;
  case CS_NAVSTO_COUPLING_PROJECTION:
    cs_navsto_projection_last_setup(connect, quant, nsp, ns->coupling_context);
    break;
  case CS_NAVSTO_COUPLING_UZAWA:
    cs_navsto_uzawa_last_setup(connect, quant, nsp, ns->coupling_context);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _err_invalid_coupling, __func__);
    break;

  }

  /* Set functions according to the discretization scheme */
  switch (nsp->space_scheme) {

  case CS_SPACE_SCHEME_CDOFB:
  case CS_SPACE_SCHEME_HHO_P0:

    /* Setup data according to the type of coupling */
    switch (nsp->coupling) {

    case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
      ns->init_scheme_context = cs_cdofb_ac_init_scheme_context;
      ns->free_scheme_context = cs_cdofb_ac_free_scheme_context;
      ns->init_velocity = NULL;
      ns->init_pressure = cs_cdofb_navsto_init_pressure;
      ns->compute_steady = NULL;

      switch (nsp->time_scheme) {

      case CS_TIME_SCHEME_STEADY:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: The Artificial Compressibility "
                  "can be used only in unsteady problems", __func__);
        break;
      case CS_TIME_SCHEME_EULER_IMPLICIT:
        ns->compute = cs_cdofb_ac_compute_implicit;
        break;
      case CS_TIME_SCHEME_THETA:
      case CS_TIME_SCHEME_CRANKNICO:
        ns->compute = cs_cdofb_ac_compute_theta;
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid time scheme for the "
                  " Artificial Compressibility coupling", __func__);
        break;

      } /* Switch */

      cs_cdofb_ac_init_common(quant, connect, time_step);
      break;

    case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY_VPP:
      /* ns->init = cs_cdofb_navsto_init_ac_vpp_context; */
      /* ns->compute = cs_cdofb_navsto_ac_vpp_compute; */
      break;

    case CS_NAVSTO_COUPLING_MONOLITHIC:
      ns->init_scheme_context = cs_cdofb_monolithic_init_scheme_context;
      ns->free_scheme_context = cs_cdofb_monolithic_free_scheme_context;
      ns->init_velocity = NULL;
      ns->init_pressure = cs_cdofb_navsto_init_pressure;
      ns->compute_steady = cs_cdofb_monolithic_compute_steady;

      switch (nsp->time_scheme) {

      case CS_TIME_SCHEME_STEADY:
        ns->compute = NULL;
        break; /* Nothing to set */

      case CS_TIME_SCHEME_EULER_IMPLICIT:
        ns->compute = cs_cdofb_monolithic_compute_implicit;
        break;
      case CS_TIME_SCHEME_THETA:
      case CS_TIME_SCHEME_CRANKNICO:
        ns->compute = cs_cdofb_monolithic_compute_theta;
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid time scheme for the monolithic coupling",
                  __func__);
        break;

      } /* Switch */

      cs_cdofb_monolithic_init_common(mesh, quant, connect, time_step);
      break;

    case CS_NAVSTO_COUPLING_PROJECTION:
      /* ns->init = cs_cdofb_navsto_init_proj_context; */
      /* ns->compute = cs_cdofb_navsto_proj_compute; */
      break;

    case CS_NAVSTO_COUPLING_UZAWA:
      ns->init_scheme_context = cs_cdofb_uzawa_init_scheme_context;
      ns->free_scheme_context = cs_cdofb_uzawa_free_scheme_context;
      ns->init_velocity = NULL;
      ns->init_pressure = cs_cdofb_navsto_init_pressure;
      ns->compute_steady = cs_cdofb_uzawa_compute_steady;

      switch (nsp->time_scheme) {

      case CS_TIME_SCHEME_EULER_IMPLICIT:
        ns->compute = cs_cdofb_uzawa_compute_implicit;
        break;
      case CS_TIME_SCHEME_THETA:
      case CS_TIME_SCHEME_CRANKNICO:
        ns->compute = cs_cdofb_uzawa_compute_theta;
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid time scheme for the Uzawa coupling", __func__);
        break;

      } /* Switch */

      cs_cdofb_uzawa_init_common(quant, connect, time_step);
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

  /* Add default post-processing related to the Navier-Stokes system */
  cs_post_add_time_mesh_dep_output(cs_navsto_system_extra_post, ns);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the context structure used to build the algebraic system
 *         This is done after the setup step.
 *
 * \param[in]  mesh      pointer to a cs_mesh_t structure
 * \param[in]  connect   pointer to a cs_cdo_connect_t structure
 * \param[in]  quant     pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts        pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_initialize(const cs_mesh_t             *mesh,
                            const cs_cdo_connect_t      *connect,
                            const cs_cdo_quantities_t   *quant,
                            const cs_time_step_t        *ts)
{
  CS_UNUSED(connect);

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
                                               ns->bf_type,
                                               ns->coupling_context);

  /* Initial conditions for the velocity */
  if (ns->init_velocity != NULL)
    ns->init_velocity(nsp, quant, ts, ns->scheme_context);

  /* Initial conditions for the pressure */
  if (ns->init_pressure != NULL)
    ns->init_pressure(nsp, quant, ts, ns->pressure);

  /* Define the advection field. Since one links the advection field to the
     face velocity this is only available for Fb schemes and should be done
     after initializing the context structure */
  cs_real_t *face_vel = NULL;
  cs_field_t  *bd_nflux = NULL;

  switch (nsp->coupling) {

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
  case CS_NAVSTO_COUPLING_MONOLITHIC:
  case CS_NAVSTO_COUPLING_UZAWA:
    {
      cs_equation_t  *mom_eq = cs_equation_by_name("momentum");
      face_vel = cs_equation_get_face_values(mom_eq);
    }
    break;

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY_VPP:
  case CS_NAVSTO_COUPLING_PROJECTION:
  default:
    bft_error(__FILE__, __LINE__, 0, _err_invalid_coupling, __func__);
    break;

  }

  const cs_flag_t loc_flag
    = CS_FLAG_FULL_LOC | cs_flag_primal_face | CS_FLAG_VECTOR;

  cs_advection_field_def_by_array(ns->adv_field, loc_flag, face_vel,
                                  false, /* the advection field is not owner */
                                  NULL);

  /* Assign the velocity boundary flux to the boundary flux for the advection
     field*/
  if (bd_nflux != NULL)
    ns->adv_field->bdy_field_id = bd_nflux->id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build, solve and update the Navier-Stokes system in case of a
 *         steady-state approach
 *
 * \param[in] mesh       pointer to a cs_mesh_t structure
 * \param[in] time_step  structure managing the time stepping
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_compute_steady_state(const cs_mesh_t        *mesh,
                                      const cs_time_step_t   *time_step)
{
  cs_navsto_system_t  *ns = cs_navsto_system;

  if (ns == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_ns));

  cs_navsto_param_t  *nsp = ns->param;

  /* Build and solve the Navier-Stokes system */
  if (nsp->time_state == CS_NAVSTO_TIME_STATE_FULL_STEADY)
    ns->compute_steady(mesh, ns->param, ns->scheme_context);

  /* Retrieve the boundary velocity flux (mass flux) and perform the update */
  cs_field_t  *nflx
    = cs_advection_field_get_field(ns->adv_field,
                                   CS_MESH_LOCATION_BOUNDARY_FACES);

  assert(nflx != NULL);
  cs_advection_field_across_boundary(ns->adv_field,
                                     time_step->t_cur,
                                     nflx->val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build, solve and update the Navier-Stokes system
 *
 * \param[in] mesh       pointer to a cs_mesh_t structure
 * \param[in] time_step  structure managing the time stepping
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_compute(const cs_mesh_t         *mesh,
                         const cs_time_step_t    *time_step)
{
  cs_navsto_system_t  *ns = cs_navsto_system;

  if (ns == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_ns));

  const cs_navsto_param_t  *nsp = ns->param;
  if (nsp->time_state == CS_NAVSTO_TIME_STATE_FULL_STEADY)
    return;

  /* Build and solve the Navier-Stokes system */
  ns->compute(mesh, nsp, ns->scheme_context);
  /* Retrieve the boundary velocity flux (mass flux) and perform the update */
  cs_field_t  *nflx
    = cs_advection_field_get_field(ns->adv_field,
                                   CS_MESH_LOCATION_BOUNDARY_FACES);

  assert(nflx != NULL);
  cs_advection_field_across_boundary(ns->adv_field, time_step->t_cur,
                                     nflx->val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations for the Navier-Stokes system
 *
 * \param[in]  connect   pointer to a cs_cdo_connect_t structure
 * \param[in]  cdoq      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_extra_op(const cs_cdo_connect_t      *connect,
                          const cs_cdo_quantities_t   *cdoq)
{
  cs_navsto_system_t  *navsto = cs_navsto_system;

  if (navsto == NULL) bft_error(__FILE__, __LINE__, 0, _(_err_empty_ns));

  const cs_navsto_param_t  *nsp = navsto->param;

  switch (nsp->space_scheme) {

  case CS_SPACE_SCHEME_CDOFB:
    cs_cdofb_navsto_extra_op(nsp, cdoq, connect, navsto->adv_field);
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
 *                              cs_gwf_t structure)
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
  CS_UNUSED(mesh_id);
  CS_UNUSED(cat_id);
  CS_UNUSED(ent_flag);
  CS_UNUSED(n_cells);
  CS_UNUSED(n_i_faces);
  CS_UNUSED(n_b_faces);
  CS_UNUSED(cell_ids);
  CS_UNUSED(i_face_ids);
  CS_UNUSED(b_face_ids);
  CS_UNUSED(time_step);

  cs_navsto_system_t  *ns = (cs_navsto_system_t *)input;

  /* TODO */
  CS_UNUSED(ns);
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
  cs_log_printf(CS_LOG_SETUP, "%s", lsepline);
  cs_log_printf(CS_LOG_SETUP, "\tSummary of the Navier-Stokes system\n");
  cs_log_printf(CS_LOG_SETUP, "%s", lsepline);

  /* Main set of numerical parameters */
  cs_navsto_param_log(ns->param);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
