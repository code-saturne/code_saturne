/* ===========================================================================
 * Routines to handle the setup of a computational domain
 * High level interface
 *============================================================================*/

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

#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

#include "cs_boundary_zone.h"
#include "cs_domain_post.h"
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the list of boundary faces attached to a wall boundary
 *         condition
 *         Function pointer to mesh location elements selection definition.
 *
 * If non-empty and not containing all elements, a list of elements
 * of the parent mesh belonging to the location should be allocated
 * (using BFT_MALLOC) and defined by this function when called.
 * This list's lifecycle is then managed by the mesh location object.
 *
 * \param[in]   input        pointer to a structure cast on-the-fly
 * \param[in]   m            pointer to associated mesh structure.
 * \param[in]   location_id  id of associated location.
 * \param[out]  n_elts       number of selected elements
 * \param[out]  elt_list     list of selected elements.
 */
/*----------------------------------------------------------------------------*/

static void
_wall_boundary_selection(void              *input,
                         const cs_mesh_t   *m,
                         int                location_id,
                         cs_lnum_t         *n_elts,
                         cs_lnum_t        **elt_ids)
{
  CS_UNUSED(location_id);

  cs_domain_t  *domain = (cs_domain_t *)input;

  /* Handle case where global domain has a temporary existence
     and has been destroyed */

  if (domain == NULL) {
    *n_elts = 0;
    *elt_ids = NULL;
    return;
  }

  cs_domain_boundary_t  *db = domain->boundary;
  cs_lnum_t  n_wall_elts = 0;
  cs_lnum_t *wall_elts = NULL;
  bool  *is_wall = NULL;

  BFT_MALLOC(is_wall, m->n_b_faces, bool);

  if (db->default_type == CS_DOMAIN_BOUNDARY_WALL) {

    for (cs_lnum_t i = 0; i < m->n_b_faces; i++)
      is_wall[i] = true;

    for (int i = 0; i < db->n_zones; i++) {
      if (db->zone_type[i] != CS_DOMAIN_BOUNDARY_WALL) {

        int z_id = db->zone_ids[i];
        const cs_boundary_zone_t  *z = cs_boundary_zone_by_id(z_id);
        const cs_lnum_t  _n_faces = z->n_faces;
        const cs_lnum_t  *_face_ids = z->face_ids;

        for (cs_lnum_t j = 0; j < _n_faces; j++)
          is_wall[_face_ids[j]] = false;

      }
    }

  }
  else { /* Wall is not the default boundary */

    for (cs_lnum_t i = 0; i < m->n_b_faces; i++)
      is_wall[i] = false;

    for (int i = 0; i < db->n_zones; i++) {
      if (db->zone_type[i] == CS_DOMAIN_BOUNDARY_WALL) {

        int z_id = db->zone_ids[i];
        const cs_boundary_zone_t  *z = cs_boundary_zone_by_id(z_id);
        const cs_lnum_t  _n_faces = z->n_faces;
        const cs_lnum_t  *_face_ids = z->face_ids;

        for (cs_lnum_t j = 0; j < _n_faces; j++)
          is_wall[_face_ids[j]] = true;

      }
    }

  } /* Which default ? */

  /* Count  */
  for (cs_lnum_t i = 0; i < m->n_b_faces; i++)
    if (is_wall[i])
      n_wall_elts++;

  if (n_wall_elts < m->n_b_faces) {

    /* Fill list  */
    BFT_MALLOC(wall_elts, n_wall_elts, cs_lnum_t);

    cs_lnum_t shift = 0;
    for (cs_lnum_t i = 0; i < m->n_b_faces; i++)
      if (is_wall[i])
        wall_elts[shift++] = i;

    assert(shift == n_wall_elts);

  }

  BFT_FREE(is_wall);

  /* Return pointers */
  *n_elts = n_wall_elts;
  *elt_ids = wall_elts;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set to true the automatic update of all advection fields
 *
 * \param[in, out]  domain    pointer to a cs_domain_t structure
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
 * \brief  Add a new zone gathering all CS_DOMAIN_BOUNDARY_WALL zone type
 *
 * \param[in]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_update_wall_zones(cs_domain_t   *domain)
{
  /* Add a new boundary zone (and also a new mesh location) related to all
     wall boundary faces */
  const char  zone_name[] = "domain_walls";

  int flag = CS_BOUNDARY_ZONE_WALL | CS_BOUNDARY_ZONE_PRIVATE;

  int  z_id = cs_boundary_zone_define_by_func(zone_name,
                                              _wall_boundary_selection,
                                              domain,
                                              flag);

  /* Allow overlay with other boundary zones used to set BCs on transport
     equations for instance (not really needed since zone is private) */
  cs_boundary_zone_set_overlay(z_id, true);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup predefined equations which are activated.
 *         At this stage, no equation is added and the space discretization
 *         scheme and the related numerical parameters are set.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_setup_predefined_equations(cs_domain_t   *domain)
{
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
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the scheme flags for the current computational domain
 *         Requirement: domain->cdo_context is alloctated
 *
 * \param[in, out]  domain       pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_set_scheme_flags(cs_domain_t    *domain)
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
                _(" Undefined type of schme to solve for eq. %s."
                  " Please check your settings."), cs_equation_get_name(eq));
    }

  } // Loop on equations

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a cs_domain_t structure
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

  domain->mesh = mesh;
  domain->mesh_quantities = mesh_quantities;

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
                                     domain->connect,
                                     domain->time_step);

  cs_evaluate_set_shared_pointers(domain->cdo_quantities,
                                  domain->connect,
                                  domain->time_step);

  cs_property_set_shared_pointers(domain->cdo_quantities,
                                  domain->connect,
                                  domain->time_step);

  cs_advection_field_set_shared_pointers(domain->cdo_quantities,
                                         domain->connect,
                                         domain->time_step);

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
  cs_equation_common_allocate(domain->connect,
                              domain->cdo_quantities,
                              domain->time_step,
                              domain->cdo_context);

  /* Set the definition of user-defined properties and/or advection
     fields (no more fields are created at this stage) */
  cs_user_cdo_finalize_setup(cs_glob_domain);

  if (cs_walldistance_is_activated())
    cs_walldistance_finalize_setup(domain->connect, domain->cdo_quantities);

  if (cs_gwf_is_activated())
    cs_gwf_finalize_setup(domain->connect, domain->cdo_quantities);

  /* Navier-Stokes system */
  if (cs_navsto_system_is_activated())
    cs_navsto_system_finalize_setup(domain->connect, domain->cdo_quantities);

  /* Last stage to define properties (when complex definition is requested) */
  cs_property_finalize_setup();

  /* Proceed to the last settings of a cs_equation_t structure
     - Assign to a cs_equation_t structure a list of function to manage this
       structure during the computation.
     - The set of functions chosen for each equation depends on the parameters
       specifying the cs_equation_t structure
     - Setup the structure related to cs_sles_*
  */

  domain->only_steady = cs_equation_finalize_setup(domain->connect,
                                                   domain->profiling);

  if (domain->only_steady)
    domain->is_last_iter = true;
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
  cs_advection_field_update(false); // operate current to previous ?

  /* Set the initial state for the groundawater flow module */
  if (cs_gwf_is_activated())
    cs_gwf_update(domain->mesh,
                  domain->connect,
                  domain->cdo_quantities,
                  domain->time_step,
                  false); // operate current to previous ?
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
