/*============================================================================
 * Manage a computational domain within the CDO framework
 *  - Physical boundary conditions attached to a domain
 *  - Properties and advection fields attached to this domain
 *  - Equations to solve on this domain
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

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <locale.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_boundary_zone.h"
#include "cs_domain_post.h"
#include "cs_evaluate.h"
#include "cs_equation_common.h"
#include "cs_gwf.h"
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_log_iteration.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"
#include "cs_prototypes.h"
#include "cs_restart.h"
#include "cs_restart_default.h"
#include "cs_source_term.h"
#include "cs_time_step.h"
#include "cs_walldistance.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_domain.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Static global variables
 *============================================================================*/

cs_domain_t  *cs_glob_domain = NULL; /* Pointer to the main computational
                                        domain used in CDO/HHO schemes */

/*============================================================================
 * Local structure definitions
 *============================================================================*/


/*============================================================================
 * Local variables
 *============================================================================*/

static const char _err_empty_domain[] =
  " Stop setting an empty cs_domain_t structure.\n"
  " Please check your settings.\n";

static double  cs_domain_kahan_time_compensation = 0.0;

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
 * \param [in]   m            pointer to associated mesh structure.
 * \param [in]   location_id  id of associated location.
 * \param [out]  n_elts       number of selected elements
 * \param [out]  elt_list     list of selected elements.
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

  cs_domain_boundary_t  *db = (cs_domain_boundary_t *)input;

  cs_lnum_t  n_wall_elts = 0;
  cs_lnum_t *wall_elts = NULL;
  bool  *is_wall = NULL;

  BFT_MALLOC(is_wall, m->n_b_faces, bool);

  if (db->default_type == CS_PARAM_BOUNDARY_WALL) {

    for (cs_lnum_t i = 0; i < m->n_b_faces; i++)
      is_wall[i] = true;

    for (int i = 0; i < db->n_zones; i++) {
      if (db->type_by_zone[i] != CS_PARAM_BOUNDARY_WALL) {

        int z_id = db->zone_ids[i];
        const cs_boundary_zone_t  *z = cs_boundary_zone_by_id(z_id);
        const cs_lnum_t  *_n_elts = cs_mesh_location_get_n_elts(z->location_id);
        const cs_lnum_t  *_elt_ids =
          cs_mesh_location_get_elt_list(z->location_id);

        for (cs_lnum_t j = 0; j < _n_elts[0]; j++)
          is_wall[_elt_ids[j]] = false;

      }
    }

  }
  else { /* Wall is not the default boundary */

    for (cs_lnum_t i = 0; i < m->n_b_faces; i++)
      is_wall[i] = false;

    for (int i = 0; i < db->n_zones; i++) {
      if (db->type_by_zone[i] == CS_PARAM_BOUNDARY_WALL) {

        int z_id = db->zone_ids[i];
        const cs_boundary_zone_t  *z = cs_boundary_zone_by_id(z_id);
        const cs_lnum_t  *_n_elts = cs_mesh_location_get_n_elts(z->location_id);
        const cs_lnum_t  *_elt_ids =
          cs_mesh_location_get_elt_list(z->location_id);

        for (cs_lnum_t j = 0; j < _n_elts[0]; j++)
          is_wall[_elt_ids[j]] = true;

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute equations which user-defined and steady-state
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_steady_user_equations(cs_domain_t   *domain)
{
  int  n_equations = cs_equation_get_n_equations();

  for (int eq_id = 0; eq_id < n_equations; eq_id++) {

    cs_equation_t  *eq = cs_equation_by_id(eq_id);

    if (cs_equation_is_steady(eq)) {

      cs_equation_type_t  type = cs_equation_get_type(eq);

      if (type == CS_EQUATION_TYPE_USER) {

        /* Define the algebraic system */
        cs_equation_build_system(domain->mesh,
                                 domain->time_step,
                                 domain->dt_cur,
                                 eq);

        /* Solve the algebraic system */
        cs_equation_solve(eq);

      } /* User-defined equation */

    } /* Steady-state equation */

  } // Loop on equations

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute user-defined equation which are time-dependent
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 * \param[in]       nt_cur     current number of iteration done
 */
/*----------------------------------------------------------------------------*/

static void
_compute_unsteady_user_equations(cs_domain_t   *domain,
                                 int            nt_cur)
{
  const int  n_equations = cs_equation_get_n_equations();

  if (nt_cur > 0) {

    for (int eq_id = 0; eq_id < n_equations; eq_id++) {

      cs_equation_t  *eq = cs_equation_by_id(eq_id);

      if (!cs_equation_is_steady(eq)) {

        cs_equation_type_t  type = cs_equation_get_type(eq);

        if (type == CS_EQUATION_TYPE_USER) {

          /* Define the algebraic system */
          if (cs_equation_needs_build(eq))
            cs_equation_build_system(domain->mesh,
                                     domain->time_step,
                                     domain->dt_cur,
                                     eq);

          /* Solve domain */
          cs_equation_solve(eq);

        } /* User-defined equation */

      } /* Unsteady equations */

    } /* Loop on equations */

  } /* nt_cur > 0 */

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize by default a cs_domain_t structure
 *
 * \return a pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

cs_domain_t *
cs_domain_create(void)
{
  cs_real_t  default_time_step = -1e13;
  cs_domain_t  *domain = NULL;

  BFT_MALLOC(domain, 1, cs_domain_t);

  domain->mesh = NULL;
  domain->mesh_quantities = NULL;
  domain->connect = NULL;
  domain->cdo_quantities = NULL;

  /* Default initialization of the time step */
  domain->is_last_iter = false;
  domain->dt_cur = default_time_step;
  domain->time_step_def = NULL;

  domain->time_step = cs_get_glob_time_step();
  domain->time_step->is_variable = 0;
  domain->time_step->is_local = 0; // In CDO, this is always equal to 0
  domain->time_step->nt_prev = 0;  // Changed in case of restart
  domain->time_step->nt_cur =  0;  // Do not modify this value
  domain->time_step->nt_max = 0;
  domain->time_step->nt_ini = 2;   // Not useful in CDO module
  domain->time_step->t_prev = 0.;
  domain->time_step->t_cur = 0.;   // Assume initial time is 0.0
  domain->time_step->t_max = 0.;

  domain->time_options.inpdt0 = 0; // standard calculation
  domain->time_options.iptlro = 0;
  domain->time_options.idtvar = 0; // constant time step by default
  domain->time_options.dtref = default_time_step;
  domain->time_options.coumax = 1.;
  domain->time_options.cflmmx = 0.99;
  domain->time_options.foumax = 10.;
  domain->time_options.varrdt = 0.1;
  domain->time_options.dtmin = default_time_step;
  domain->time_options.dtmax = default_time_step;
  domain->time_options.relxst = 0.7; // Not useful in CDO schemes

  domain->scheme_flag = 0;
  domain->only_steady = true;
  domain->force_advfield_update = false;

  /* Other options */
  domain->output_nt = -1;
  domain->verbosity = 1;
  domain->profiling = false;

  /* Add predefined properties */
  cs_property_t  *pty = cs_property_add("unity", CS_PROPERTY_ISO);

  cs_property_def_iso_by_value(pty, "cells", 1.0);

  /* Allocate domain builder */
  BFT_MALLOC(domain->boundary_def, 1, cs_domain_boundary_t);
  domain->boundary_def->default_type = CS_PARAM_BOUNDARY_WALL; // Set by default
  domain->boundary_def->n_zones = 0;
  domain->boundary_def->zone_ids = NULL;
  domain->boundary_def->type_by_zone = NULL;

  /* Monitoring */
  CS_TIMER_COUNTER_INIT(domain->tcp); // domain post
  CS_TIMER_COUNTER_INIT(domain->tcs); // domain setup

  return domain;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_domain_t structure
 *
 * \param[in, out]   domain    pointer to the cs_domain_t structure to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_domain_t *
cs_domain_free(cs_domain_t   *domain)
{
  if (domain == NULL)
    return domain;

  cs_domain_post_finalize();

  /* cs_mesh_t and cs_mesh_quantities_t structure are not freed since they
     are only shared */
  domain->mesh = NULL;
  domain->mesh_quantities = NULL;

  BFT_FREE(domain->boundary_def->zone_ids);
  BFT_FREE(domain->boundary_def->type_by_zone);
  BFT_FREE(domain->boundary_def);

  domain->time_step_def = cs_xdef_free(domain->time_step_def);
  domain->time_step = NULL;

  /* Print monitoring information */
  cs_equation_log_monitoring();

  /* Free memory related to equations */
  cs_equation_destroy_all();

  /* Free memory related to advection fields */
  cs_advection_field_destroy_all();

  /* Free memory related to properties */
  cs_property_destroy_all();

  /* Free memory related to the groundwater flow module */
  cs_gwf_destroy_all();

  /* Free common structures relatated to equations */
  cs_equation_free_common_structures(domain->scheme_flag);

  /* Free CDO structures related to geometric quantities and connectivity */
  domain->cdo_quantities = cs_cdo_quantities_free(domain->cdo_quantities);
  domain->connect = cs_cdo_connect_free(domain->connect);

  BFT_FREE(domain);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the default boundary related to this domain
 *
 * \param[in, out]   domain       pointer to a cs_domain_t structure
 * \param[in]        type         type of boundary to set
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_set_default_boundary(cs_domain_t                *domain,
                               cs_param_boundary_type_t    type)
{
  if (domain == NULL) bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  if (type == CS_PARAM_BOUNDARY_WALL || type == CS_PARAM_BOUNDARY_SYMMETRY)
    domain->boundary_def->default_type = type;
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of boundary by default.\n"
                " Valid choice is CS_PARAM_BOUNDARY_WALL or"
                " CS_PARAM_BOUNDARY_SYMMETRY."));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a boundary type defined on a mesh location
 *
 * \param[in, out]  domain       pointer to a cs_domain_t structure
 * \param[in]       type         type of boundary to set
 * \param[in]       zone_name    name of the zone related to this boundary
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_add_boundary(cs_domain_t                *domain,
                       cs_param_boundary_type_t    type,
                       const char                 *zone_name)
{
  if (domain == NULL) bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  const cs_boundary_zone_t  *zone = cs_boundary_zone_by_name(zone_name);

  if (zone == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid zone name %s.\n"
                " This zone is not already defined.\n"), zone_name);

  int  new_id = domain->boundary_def->n_zones;

  domain->boundary_def->n_zones++;

  BFT_REALLOC(domain->boundary_def->zone_ids,
              domain->boundary_def->n_zones, int);
  domain->boundary_def->zone_ids[new_id] = zone->id;

  BFT_REALLOC(domain->boundary_def->type_by_zone,
              domain->boundary_def->n_zones, cs_param_boundary_type_t);
  domain->boundary_def->type_by_zone[new_id] = type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set parameters for unsteady computations: the max number of time
 *         steps or the final physical time of the simulation
 *
 * \param[in, out]  domain    pointer to a cs_domain_t structure
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
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set auxiliary parameters related to the way output is done
 *
 * \param[in, out]  domain      pointer to a cs_domain_t structure
 * \param[in]       nt_list     output frequency into the listing
 * \param[in]       verbosity   level of information displayed
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_set_output_param(cs_domain_t       *domain,
                           int                nt_list,
                           int                verbosity)
{
  if (domain == NULL) bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  domain->output_nt = nt_list;
  if (domain->output_nt == 0)
    domain->output_nt = -1;

  domain->verbosity = verbosity;
}

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
  if (domain == NULL) bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  domain->force_advfield_update = true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set auxiliary parameters related to a cs_domain_t structure
 *
 * \param[in, out]  domain    pointer to a cs_domain_t structure
 * \param[in]       key       key related to the parameter to set
 * \param[in]       keyval    value related to the parameter to set
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_set_advanced_param(cs_domain_t       *domain,
                             cs_domain_key_t    key,
                             const char        *keyval)
{
  CS_UNUSED(keyval);

  if (domain == NULL) bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  switch(key) {

  case CS_DOMAIN_PROFILING:
    domain->profiling = true;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key for setting a cs_domain_t structure."));

  } /* Switch on keys */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the value of the time step thanks to a predefined function
 *
 * \param[in, out] domain      pointer to a cs_domain_t structure
 * \param[in]      func        pointer to a cs_timestep_func_t function
 * \param[in]      func_input  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_def_time_step_by_function(cs_domain_t          *domain,
                                    cs_timestep_func_t   *func,
                                    void                 *func_input)
{
  if (domain == NULL) bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  domain->time_step->is_variable = 1; // not constant time step
  domain->time_options.idtvar = 1;    /* uniform in space but can change
                                         from one time step to the other */

  cs_xdef_timestep_input_t  def = {.input = func_input,
                                   .func = func};

  domain->time_step_def = cs_xdef_timestep_create(CS_XDEF_BY_TIME_FUNCTION,
                                                  0, // state flag
                                                  0, // meta flag
                                                  &def);

  /* Default initialization.
     To be changed at first call to cs_domain_time_step_increment() */

  domain->dt_cur = domain->time_step->t_max;
  domain->time_options.dtref =  domain->time_step->t_max;
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

  domain->time_step->is_variable = 0; // constant time step
  domain->time_options.idtvar = 0;    // constant time step by default

  domain->time_step_def = cs_xdef_timestep_create(CS_XDEF_BY_VALUE,
                                                  0, // state flag
                                                  0, // meta flag
                                                  &dt);

  domain->dt_cur = dt;
  domain->time_options.dtref = domain->dt_cur;
  domain->time_options.dtmin = domain->dt_cur;
  domain->time_options.dtmax = domain->dt_cur;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add new mesh locations related to domain boundaries from existing
 *         mesh locations
 *
 * \param[in]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_update_mesh_locations(cs_domain_t   *domain)
{
  /* Add a new boundary zone (and also a new mesh location) related to all
     wall boundary faces */
  const char *zone_name =
    cs_param_get_boundary_domain_name(CS_PARAM_BOUNDARY_WALL);

  int  z_id = cs_boundary_zone_define_by_func(zone_name,
                                              _wall_boundary_selection,
                                              domain->boundary_def,
                                              CS_BOUNDARY_ZONE_CDO_DOMAIN);

  /* Allow overlay with other boundary zones used to set BCs on transport
     equations for instance */
  cs_boundary_zone_set_overlay(z_id, true);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup predefined equations which are activated
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_setup_predefined_equations(cs_domain_t   *domain)
{
  CS_UNUSED(domain);

  /* Wall distance */
  if (cs_walldistance_is_activated())
    cs_walldistance_setup();

  /* Groundwater flow module */
  if (cs_gwf_is_activated())
    cs_gwf_init_setup();

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the scheme flag for the current computational domain
 *
 * \param[in, out]  domain            pointer to a cs_domain_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_set_scheme_flag(cs_domain_t    *domain)
{
  if (domain == NULL) bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  /* Define a scheme flag for the current domain */
  const int  n_equations = cs_equation_get_n_equations();
  for (int eq_id = 0; eq_id < n_equations; eq_id++) {

    cs_equation_t  *eq = cs_equation_by_id(eq_id);
    cs_space_scheme_t  scheme = cs_equation_get_space_scheme(eq);
    int  vardim = cs_equation_get_var_dim(eq);

    if (vardim == 1)
      domain->scheme_flag |= CS_SCHEME_FLAG_SCALAR;
    else if (vardim == 3)
      domain->scheme_flag |= CS_SCHEME_FLAG_VECTOR;

    switch (scheme) {

    case CS_SPACE_SCHEME_CDOVB:
      domain->scheme_flag |= CS_SCHEME_FLAG_CDOVB | CS_SCHEME_FLAG_POLY0;
      break;
    case CS_SPACE_SCHEME_CDOVCB:
      domain->scheme_flag |= CS_SCHEME_FLAG_CDOVCB | CS_SCHEME_FLAG_POLY0;
      break;
    case CS_SPACE_SCHEME_CDOFB:
      domain->scheme_flag |= CS_SCHEME_FLAG_CDOFB | CS_SCHEME_FLAG_POLY0;
      break;
    case CS_SPACE_SCHEME_HHO_P0:
      domain->scheme_flag |= CS_SCHEME_FLAG_HHO | CS_SCHEME_FLAG_POLY0;
      assert(cs_equation_get_space_poly_degree(eq) == 0);
      break;
    case CS_SPACE_SCHEME_HHO_P1:
      domain->scheme_flag |= CS_SCHEME_FLAG_HHO | CS_SCHEME_FLAG_POLY1;
      assert(cs_equation_get_space_poly_degree(eq) == 1);
      break;
    case CS_SPACE_SCHEME_HHO_P2:
      domain->scheme_flag |= CS_SCHEME_FLAG_HHO | CS_SCHEME_FLAG_POLY2;
      assert(cs_equation_get_space_poly_degree(eq) == 2);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Undefined type of schme to solve for eq. %s."
                  " Please check your settings."), cs_equation_get_name(eq));
    }

  } // Loop on equations

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
  if (domain == NULL) bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  domain->mesh = mesh;
  domain->mesh_quantities = mesh_quantities;

  /* Build additional connectivity structures
     Update mesh structure with range set structures */
  domain->connect = cs_cdo_connect_init(mesh, domain->scheme_flag);

  /* Default = CS_CDO_CC_SATURNE but can be modify by the user */
  cs_cdo_cell_center_algo_t  cc_algo =
    cs_user_cdo_geometric_settings();

  /* Build additional mesh quantities in a seperate structure */
  domain->cdo_quantities =  cs_cdo_quantities_build(cc_algo,
                                                    mesh,
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
    cs_user_cdo_setup_gwf(domain);

    /* Add if needed new terms (as diffusion or reaction) to tracer equations
       according to the settings */
    cs_gwf_add_tracer_terms();

  }

  /* Allocate all fields created during the setup stage */
  cs_field_allocate_or_map_all();

  /* Initialization default post-processing for the computational domain */
  cs_domain_post_init(domain->dt_cur, domain->cdo_quantities);

  /* Allocate common structures for solving equations */
  cs_equation_allocate_common_structures(domain->connect,
                                         domain->cdo_quantities,
                                         domain->time_step,
                                         domain->scheme_flag);

  /* Set the definition of user-defined properties and/or advection
     fields (no more fields are created at this stage) */
  cs_user_cdo_finalize_setup(cs_glob_domain);

  if (cs_walldistance_is_activated())
    cs_walldistance_finalize_setup(domain->connect, domain->cdo_quantities);

  if (cs_gwf_is_activated())
    cs_gwf_finalize_setup(domain->connect, domain->cdo_quantities);

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
/*!
 * \brief  Check if one needs to continue iterations in time
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_domain_needs_iteration(cs_domain_t  *domain)
{
  bool  one_more_iter = true;

  cs_time_step_t  *ts = domain->time_step;

  if (ts->nt_max > 0) // nt_max has been set
    if (ts->nt_cur > ts->nt_max)
      one_more_iter = false;

  if (ts->t_max > 0) // t_max has been set
    if (ts->t_cur > ts->t_max)
      one_more_iter = false;

  if (domain->only_steady && ts->nt_cur > 0)
    one_more_iter = false;

  if (!domain->only_steady && ts->nt_max <= 0 && ts->t_max <= 0)
    one_more_iter = false;

  return one_more_iter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if an ouput is requested according to the domain setting
 *
 * \param[in]   domain    pointer to a cs_domain_t structure
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_domain_needs_log(const cs_domain_t      *domain)
{
  const cs_time_step_t  *ts = domain->time_step;

  if (domain->verbosity < 0)
    return false;

  if (domain->only_steady)
    return true;

  if (domain->output_nt > 0)
    if (ts->nt_cur % domain->output_nt == 0)
      return true;

  if (domain->is_last_iter)
    return true;

  return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the current time step for this new time iteration
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_define_current_time_step(cs_domain_t   *domain)
{
  if (domain == NULL) bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  if (domain->only_steady)
    return;

  const cs_time_step_t  *ts = domain->time_step;
  const double  t_cur = ts->t_cur;
  const int  nt_cur = ts->nt_cur;

  cs_xdef_t  *ts_def = domain->time_step_def;

  if (ts_def->type != CS_XDEF_BY_VALUE) { /* dt_cur may change */

    if (ts_def->type == CS_XDEF_BY_TIME_FUNCTION) {

      /* Compute current time step */
      cs_xdef_timestep_input_t  *param =
        (cs_xdef_timestep_input_t *)ts_def->input;
      domain->dt_cur = param->func(nt_cur, t_cur, param->input);

      /* Update time_options */
      double  dtmin = CS_MIN(domain->time_options.dtmin, domain->dt_cur);
      double  dtmax = CS_MAX(domain->time_options.dtmax, domain->dt_cur);

      domain->time_options.dtmin = dtmin;
      domain->time_options.dtmax = dtmax;
      // TODO: Check how the following value is set in FORTRAN
      // domain->time_options.dtref = 0.5*(dtmin + dtmax);
      if (domain->time_options.dtref < 0) // Should be the initial val.
        domain->time_options.dtref = domain->dt_cur;

    }
    else
      bft_error(__FILE__, __LINE__, 0,
                " Invalid way of defining the current time step.\n"
                " Please modify your settings.");

    cs_domain_post_update(domain->dt_cur);
  }

  /* Check if this is the last iteration */
  if (ts->t_max > 0) // t_max has been set
    if (t_cur + domain->dt_cur > ts->t_max)
      domain->is_last_iter = true;
  if (ts->nt_max > 0) // nt_max has been set
    if (nt_cur + 1 > ts->nt_max)
      domain->is_last_iter = true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update time step after one temporal iteration
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_increment_time(cs_domain_t  *domain)
{
  cs_time_step_t  *ts = domain->time_step;

  /* Increment time iteration */
  ts->nt_cur++;
  ts->t_prev = ts->t_cur;

  /* Use Kahan's trick to limit the truncation error */
  double  z = domain->dt_cur - cs_domain_kahan_time_compensation;
  double  t = ts->t_cur + z;

  cs_domain_kahan_time_compensation = (t - ts->t_cur) - z;
  ts->t_cur = t;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve all the equations of a computational domain for one time step
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_solve(cs_domain_t  *domain)
{
  int  nt_cur = domain->time_step->nt_cur;
  bool  do_output = cs_domain_needs_log(domain);

  /* Setup step for all equations */
  if (nt_cur == 0) {

    /* Output information */
    if (domain->only_steady) {
      cs_log_printf(CS_LOG_DEFAULT, "\n%s", lsepline);
      cs_log_printf(CS_LOG_DEFAULT, "#      Solve steady-state problem(s)\n");
      cs_log_printf(CS_LOG_DEFAULT, "%s", lsepline);
    }
    else if (do_output) {
      cs_log_printf(CS_LOG_DEFAULT, "\n%s", lsepline);
      cs_log_printf(CS_LOG_DEFAULT,
                    "-ite- %5d; time= %5.3e s; dt= %5.3e >> Solve domain\n",
                    nt_cur, domain->time_step->t_cur, domain->dt_cur);
      cs_log_printf(CS_LOG_DEFAULT, "%s", lsepline);
    }

    /* Predefined equation for the computation of the wall distance */
    if (cs_walldistance_is_activated())
      cs_walldistance_compute(domain->mesh,
                              domain->time_step,
                              domain->dt_cur,
                              domain->connect,
                              domain->cdo_quantities);

    /* If unsteady, only initialization is done, otherwise one makes the whole
       computation */
    if (cs_gwf_is_activated())
      cs_gwf_compute(domain->mesh,
                     domain->time_step,
                     domain->dt_cur,
                     domain->connect,
                     domain->cdo_quantities);

    /* User-defined equations */
    _compute_steady_user_equations(domain);

    /* Only initialization is done */
    _compute_unsteady_user_equations(domain, nt_cur);

  }
  else { /* nt_cur > 0: solve unsteady problems */

    /* Output information */
    if (do_output) {
      cs_log_printf(CS_LOG_DEFAULT, "\n%s", lsepline);
      cs_log_printf(CS_LOG_DEFAULT,
                    "-ite- %5d; time = %5.3e s >> Solve domain\n",
                    nt_cur, domain->time_step->t_cur);
      cs_log_printf(CS_LOG_DEFAULT, "%s", lsepline);
    }

    if (cs_gwf_is_activated())
      cs_gwf_compute(domain->mesh,
                     domain->time_step,
                     domain->dt_cur,
                     domain->connect,
                     domain->cdo_quantities);

    /* User-defined equations */
    _compute_unsteady_user_equations(domain, nt_cur);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Process the computational domain after the resolution
 *
 * \param[in]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_process_after_solve(cs_domain_t  *domain)
{
  cs_timer_t  t0 = cs_timer_time();

  /* Pre-stage for post-processing for the current time step */
  cs_domain_post_activate(domain->time_step);

  /* Extra-operations */
  /* ================ */

  /* Predefined extra-operations related to advection fields */
  if (domain->force_advfield_update)
    cs_advection_field_update(true);

  /* User-defined extra operations */
  cs_user_cdo_extra_op(domain);

  /* Log output */
  if (cs_domain_needs_log(domain))
    cs_log_iteration();

  /* Post-processing */
  /* =============== */

  cs_domain_post(domain->time_step);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(domain->tcp), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read a restart file for the CDO module
 *
 * \param[in]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_read_restart(const cs_domain_t  *domain)
{
  CS_UNUSED(domain);

  cs_restart_t  *restart = cs_restart_create("main", // restart file name
                                             NULL,   // directory name
                                             CS_RESTART_MODE_READ);

  const char err_i_val[] = N_("Restart mismatch for: %s\n"
                              "read: %d\n"
                              "expected: %d.");
  int i_val;

  /* Read a new section: version */
  int  version = 400000;
  cs_restart_read_section(restart,
                          "code_saturne:checkpoint:main:version", // secname
                          CS_MESH_LOCATION_NONE,                  // ml_id
                          1,                                      // nb. values
                          CS_TYPE_cs_int_t,                       // val. type
                          &i_val);                                // value(s)

  if (i_val != version)
    bft_error(__FILE__, __LINE__, 0, _(err_i_val),
              "code_saturne:checkpoint:main:version", version, i_val);

  /* Read a new section: field information */
  cs_map_name_to_id_t  *old_field_map = NULL;

  cs_restart_read_field_info(restart, &old_field_map);

  /* Read a new section */
  int  n_equations = cs_equation_get_n_equations();
  cs_restart_read_section(restart,
                          "cdo:n_equations",
                          CS_MESH_LOCATION_NONE,
                          1,
                          CS_TYPE_cs_int_t,
                          &i_val);

  if (i_val != n_equations)
    bft_error(__FILE__, __LINE__, 0, _(err_i_val),
              "cdo:n_equations", n_equations, i_val);

  /* Read a new section */
  int  n_properties = cs_property_get_n_properties();
  cs_restart_read_section(restart,
                          "cdo:n_properties",
                          CS_MESH_LOCATION_NONE,
                          1,
                          CS_TYPE_cs_int_t,
                          &i_val);

  if (i_val != n_properties)
    bft_error(__FILE__, __LINE__, 0, _(err_i_val),
              "cdo:n_properties", n_properties, i_val);

  /* Read a new section */
  int  n_adv_fields = cs_advection_field_get_n_fields();
  cs_restart_read_section(restart,
                          "cdo:n_adv_fields",
                          CS_MESH_LOCATION_NONE,
                          1,
                          CS_TYPE_cs_int_t,
                          &i_val);

  if (i_val != n_adv_fields)
    bft_error(__FILE__, __LINE__, 0, _(err_i_val),
              "cdo:n_adv_fields", n_adv_fields, i_val);

  /* Read a new section: activation or not of the groundwater flow module */
  int  igwf = 0; // not activated by default
  if (cs_gwf_is_activated()) igwf = 1;
  cs_restart_read_section(restart,
                          "groundwater_flow_module",
                          CS_MESH_LOCATION_NONE,
                          1,
                          CS_TYPE_cs_int_t,
                          &i_val);

  if (i_val != igwf)
    bft_error(__FILE__, __LINE__, 0, _(err_i_val),
              "groundwater_flow_module", igwf, i_val);

  /* Read a new section: computation or not of the wall distance */
  int  iwall = 0;
  if (cs_walldistance_is_activated()) iwall = 1;
  cs_restart_read_section(restart,
                          "wall_distance",
                          CS_MESH_LOCATION_NONE,
                          1,
                          CS_TYPE_cs_int_t,
                          &i_val);

  if (i_val != iwall)
    bft_error(__FILE__, __LINE__, 0, _(err_i_val),
              "wall_distance", iwall, i_val);

  /* Read a new section: number of computed time steps */
  int  nt_cur = 0;
  cs_restart_read_section(restart,
                          "cur_time_step",
                          CS_MESH_LOCATION_NONE,
                          1,
                          CS_TYPE_cs_int_t,
                          &nt_cur);

  /* Read a new section: number of computed time steps */
  cs_real_t  t_cur = 0;
  cs_restart_read_section(restart,
                          "cur_time",
                          CS_MESH_LOCATION_NONE,
                          1,
                          CS_TYPE_cs_real_t,
                          &t_cur);

  cs_time_step_redefine_cur(nt_cur, t_cur);

  /* Main variables */
  int  t_id_flag = 0; // Only current values
  cs_restart_read_variables(restart, old_field_map, t_id_flag, NULL);

  cs_map_name_to_id_destroy(&old_field_map);

  // TODO: read field values for previous time step if needed

  int n_fields = cs_field_n_fields();
  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    cs_field_current_to_previous(f);
  }

  /* Finalize restart process */
  cs_restart_destroy(&restart);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Write a restart file for the CDO module
 *
 * \param[in]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_write_restart(const cs_domain_t  *domain)
{
  cs_restart_t  *restart = cs_restart_create("main", // restart file name
                                             NULL,   // directory name
                                             CS_RESTART_MODE_WRITE);

  /* Write a new section: version */
  int  version = 400000;
  cs_restart_write_section(restart,
                           "code_saturne:checkpoint:main:version", // secname
                           CS_MESH_LOCATION_NONE,                  // ml_id
                           1,                                      // nb. values
                           CS_TYPE_cs_int_t,                       // val. type
                           &version);                              // value(s)

  /* Write a new section: field information */
  cs_restart_write_field_info(restart);

  /* Write a new section */
  int  n_equations = cs_equation_get_n_equations();
  cs_restart_write_section(restart,
                           "cdo:n_equations",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_cs_int_t,
                           &n_equations);

  /* Write a new section */
  int  n_properties = cs_property_get_n_properties();
  cs_restart_write_section(restart,
                           "cdo:n_properties",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_cs_int_t,
                           &n_properties);

  /* Write a new section */
  int  n_adv_fields = cs_advection_field_get_n_fields();
  cs_restart_write_section(restart,
                           "cdo:n_adv_fields",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_cs_int_t,
                           &n_adv_fields);

  /* Write a new section: activation or not of the groundwater flow module */
  int  igwf = 0; // not activated by default
  if (cs_gwf_is_activated()) igwf = 1;
  cs_restart_write_section(restart,
                           "groundwater_flow_module",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_cs_int_t,
                           &igwf);

  /* Write a new section: computation or not of the wall distance */
  int  iwall = 0;
  if (cs_walldistance_is_activated()) iwall = 1;
  cs_restart_write_section(restart,
                           "wall_distance",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_cs_int_t,
                           &iwall);

  /* Write a new section: number of computed time steps */
  int  ntcabs = domain->time_step->nt_cur;
  cs_restart_write_section(restart,
                           "cur_time_step",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_cs_int_t,
                           &ntcabs);

  /* Read a new section: number of computed time steps */
  cs_real_t  ttcabs = domain->time_step->t_cur;
  cs_restart_write_section(restart,
                           "cur_time",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_cs_real_t,
                           &ttcabs);

  /* Main variables */
  int  t_id_flag = 0; // Only current values
  cs_restart_write_variables(restart, t_id_flag, NULL);

  // TODO: write field values for previous time step if needed

  /* Finalize restart process */
  cs_restart_destroy(&restart);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of a cs_domain_t structure
 *
 * \param[in]   domain    pointer to the cs_domain_t structure to summarize
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_summary(const cs_domain_t   *domain)
{
  if (domain == NULL)
    return;

  /* Output information */
  cs_log_printf(CS_LOG_SETUP, "\n%s", lsepline);
  cs_log_printf(CS_LOG_SETUP, "\tSummary of domain settings\n");
  cs_log_printf(CS_LOG_SETUP, "%s", lsepline);

  /* Boundary */
  cs_domain_boundary_t  *bdy = domain->boundary_def;

  cs_log_printf(CS_LOG_SETUP, "\n  Domain boundary by default: ");
  switch (bdy->default_type) {
  case CS_PARAM_BOUNDARY_WALL:
    cs_log_printf(CS_LOG_SETUP, " wall\n");
    break;
  case CS_PARAM_BOUNDARY_SYMMETRY:
    cs_log_printf(CS_LOG_SETUP, " symmetry\n");
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid boundary by default.\n"
                " Please modify your settings."));
  }

  /* Number of border faces for each type of boundary */
  for (int i = 0; i < bdy->n_zones; i++) {

    const cs_boundary_zone_t *z = cs_boundary_zone_by_id(bdy->zone_ids[i]);

    cs_gnum_t  n_g_elts = (cs_gnum_t)z->n_faces;
    if (cs_glob_n_ranks > 1)
      cs_parall_counter(&n_g_elts, 1);

    cs_log_printf(CS_LOG_SETUP, " %s: %u,", z->name, (unsigned int)n_g_elts);
    switch (bdy->type_by_zone[i]) {
    case CS_PARAM_BOUNDARY_INLET:
      cs_log_printf(CS_LOG_SETUP, " inlet\n");
      break;
    case CS_PARAM_BOUNDARY_OUTLET:
      cs_log_printf(CS_LOG_SETUP, " outlet\n");
      break;
    case CS_PARAM_BOUNDARY_WALL:
      cs_log_printf(CS_LOG_SETUP, " wall\n");
      break;
    case CS_PARAM_BOUNDARY_SYMMETRY:
      cs_log_printf(CS_LOG_SETUP, " symmetry\n");
      break;
    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid boundary by default.\n"
                  " Please modify your settings."));
    }

  }

  /* Time step summary */
  cs_log_printf(CS_LOG_SETUP, "\n  Time step information\n");
  if (domain->only_steady)
    cs_log_printf(CS_LOG_SETUP, "  >> Steady-state computation");

  else { /* Time information */

    cs_log_printf(CS_LOG_SETUP, "  >> Time step status:");
    if (domain->time_options.idtvar == 0)
      cs_log_printf(CS_LOG_SETUP, "  constant\n");
    else if (domain->time_options.idtvar == 1)
      cs_log_printf(CS_LOG_SETUP, "  variable in time\n");
    else
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid idtvar value for the CDO module.\n"));

    cs_xdef_log(domain->time_step_def);

    if (domain->time_step->t_max > 0.)
      cs_log_printf(CS_LOG_SETUP, "%-30s %5.3e\n",
                    "  >> Final simulation time:", domain->time_step->t_max);
    if (domain->time_step->nt_max > 0)
      cs_log_printf(CS_LOG_SETUP, "%-30s %9d\n",
                    "  >> Final time step:", domain->time_step->nt_max);

  }
  cs_log_printf(CS_LOG_SETUP, "\n");

  /* Summary for each equation */
  cs_equation_log_setup();

  if (domain->verbosity > 0) {

    /* Properties */
    cs_property_log_setup();

    /* Advection fields */
    cs_advection_field_log_setup();

    /* Summary of the groundwater module */
    cs_gwf_log_setup();

  } /* Domain->verbosity > 0 */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
