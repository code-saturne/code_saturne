/*============================================================================
 * Manage a computational domain within the CDO framework
 *  - Physical boundary conditions attached to a domain
 *  - Properties and advection fields attached to this domain
 *  - Equations to solve on this domain
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

struct _cs_domain_boundary_t {

  cs_param_boundary_type_t    default_boundary; // boundary set by default

  cs_lnum_t                   n_b_faces;        // number of boundary faces
  cs_param_boundary_type_t   *b_face_types;     // type of each boundary face

  /*
     A mesh location attached to a domain boundary can be compound of several
     existing mesh locations, hereafter called sub mesh locations
  */

  /* Mesh location id related to a domain boundary (-1 if not useful)
     Mesh locations are generated automatically with the following reserved
     names: "domain_walls",
            "domain_inlets",
            "domain_outlets",
            "domain_symmetries"
   */
  int        autogen_ml_ids[CS_PARAM_N_BOUNDARY_TYPES];
  int        n_sub_ml_ids[CS_PARAM_N_BOUNDARY_TYPES];
  int       *sub_ml_id_lst[CS_PARAM_N_BOUNDARY_TYPES];

  /* Number of border faces related to each type of boundary */
  cs_lnum_t  n_type_elts[CS_PARAM_N_BOUNDARY_TYPES];

};

/*============================================================================
 * Local variables
 *============================================================================*/

static const char _err_empty_domain[] =
  " Stop setting an empty cs_domain_t structure.\n"
  " Please check your settings.\n";

static const char
_domain_boundary_ml_name[CS_PARAM_N_BOUNDARY_TYPES][CS_BASE_STRING_LEN] =
  { N_("domain_walls"),
    N_("domain_inlets"),
    N_("domain_outlets"),
    N_("domain_symmetries") };

static double  cs_domain_kahan_time_compensation = 0.0;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a new cs_domain_boundary_t structure
 *
 * \return a pointer to a new allocated cs_domain_boundary_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_domain_boundary_t *
_create_domain_boundaries(void)
{
  cs_domain_boundary_t  *dby = NULL;

  BFT_MALLOC(dby, 1, cs_domain_boundary_t);

  dby->default_boundary = CS_PARAM_N_BOUNDARY_TYPES;
  dby->n_b_faces = -1;
  dby->b_face_types = NULL;

  for (int i = 0; i < CS_PARAM_N_BOUNDARY_TYPES; i++) {
    dby->autogen_ml_ids[i] = -1;
    dby->n_sub_ml_ids[i] = 0;
    dby->sub_ml_id_lst[i] = NULL;
    dby->n_type_elts[i] = 0;
  }

  return dby;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a cs_domain_boundary_t structure
 *
 * \param[in]   n_b_faces   number of boundary faces
 *
 * \return a pointer to a new allocated cs_domain_boundary_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_build_domain_boundaries(cs_lnum_t               n_b_faces,
                         cs_domain_boundary_t   *dby)
{
  /* Sanity check */
  assert(dby !=NULL);

  dby->n_b_faces = n_b_faces;

  /* Define the b_face_types array */
  BFT_MALLOC(dby->b_face_types, n_b_faces, cs_param_boundary_type_t);
  for (int i = 0; i < n_b_faces; i++)
    dby->b_face_types[i] = dby->default_boundary;

  for (unsigned int type = 0; type < CS_PARAM_N_BOUNDARY_TYPES; type++) {
    if (dby->autogen_ml_ids[type] > -1 && type != dby->default_boundary) {

      const int  ml_id = dby->autogen_ml_ids[type];
      const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(ml_id);
      const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(ml_id);

      if (elt_ids == NULL)
        for (cs_lnum_t i = 0; i < n_elts[0]; i++)
          dby->b_face_types[i] = type;
      else
        for (cs_lnum_t i = 0; i < n_elts[0]; i++)
          dby->b_face_types[elt_ids[i]] = type;

    }
  }

  /* Count how many border faces are associated to each type of boundary */
  cs_gnum_t  error_count = 0;

  for (cs_lnum_t i = 0; i < n_b_faces; i++) {
    if (dby->b_face_types[i] == CS_PARAM_N_BOUNDARY_TYPES)
      error_count++;
    else
      dby->n_type_elts[dby->b_face_types[i]] += 1;
  }

  if (cs_glob_n_ranks > 1)
    cs_parall_counter(&error_count, 1);

  if (error_count > 0)
    bft_error(__FILE__, __LINE__, 0,
              _(" Problem detected during the setup.\n"
                " %lu boundary faces have no boundary type."), error_count);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_domain_boundary_t structure
 *
 * \param[in,out]  dby       pointer to the cs_domain_t structure to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

static cs_domain_boundary_t *
_free_domain_boundaries(cs_domain_boundary_t   *dby)
{
  if (dby == NULL)
    return dby;

  BFT_FREE(dby->b_face_types);

  for (int type = 0; type < CS_PARAM_N_BOUNDARY_TYPES; type++)
    BFT_FREE(dby->sub_ml_id_lst[type]);

  BFT_FREE(dby);

  return NULL;
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
  for (int eq_id = 0; eq_id < domain->n_equations; eq_id++) {

    cs_equation_t  *eq = domain->equations[eq_id];

    if (cs_equation_is_steady(eq)) {

      cs_equation_type_t  type = cs_equation_get_type(eq);

      if (type == CS_EQUATION_TYPE_USER) {

        /* Initialize system before resolution for all equations
           - create system builder
           - initialize field according to initial conditions
           - initialize source term */
        cs_equation_init_system(domain->mesh, eq);

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
  if (nt_cur == 0) { /* Initialization */

    for (int eq_id = 0; eq_id < domain->n_equations; eq_id++) {

      cs_equation_t  *eq = domain->equations[eq_id];

      if (!cs_equation_is_steady(eq)) { // Unsteady eq.

        cs_equation_type_t  type = cs_equation_get_type(eq);

        if (type == CS_EQUATION_TYPE_USER) {

          /* Initialize system before resolution for all equations
             - create system builder
             - initialize field according to initial conditions
             - initialize source term */
          cs_equation_init_system(domain->mesh, eq);

        } /* User-defined equation */

      } /* Unsteady equations */

    } /* Loop on equations */

  } /* nt_cur == 0 */
  else {

    for (int eq_id = 0; eq_id < domain->n_equations; eq_id++) {

      cs_equation_t  *eq = domain->equations[eq_id];

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
  domain->time_step_def_type = CS_PARAM_DEF_BY_VALUE;
  domain->time_step_def.get.val = default_time_step;

  BFT_MALLOC(domain->time_step, 1, cs_time_step_t);
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

  /* Equations */
  domain->n_equations = 0;
  domain->n_predef_equations = 0;
  domain->n_user_equations = 0;
  domain->equations = NULL;

  domain->scheme_flag = 0;
  domain->only_steady = true;
  domain->force_advfield_update = false;

  /* Other options */
  domain->output_nt = -1;
  domain->verbosity = 1;
  domain->profiling = false;

  /* Predefined equations or modules */
  domain->richards_eq_id = -1;
  domain->wall_distance_eq_id = -1;
  domain->gwf = NULL;

  /* Specify the "physical" domain boundaries. Define a mesh location for
     each boundary type */
  domain->boundaries = _create_domain_boundaries();

  /* Initialize properties */
  domain->n_properties = 0;
  domain->properties = NULL;

  /* Add predefined properties */
  cs_property_t  *pty = cs_domain_add_property(domain, "unity", "isotropic", 1);
  cs_property_iso_def_by_value(pty, "cells", 1.0);

  /* Advection fields */
  domain->n_adv_fields = 0;
  domain->adv_fields = NULL;

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

  domain->boundaries = _free_domain_boundaries(domain->boundaries);

  BFT_FREE(domain->time_step);

  if (domain->gwf != NULL)
    domain->gwf = cs_gwf_finalize(domain->gwf);

  /* Free properties */
  for (int i = 0; i < domain->n_properties; i++)
    domain->properties[i] = cs_property_free(domain->properties[i]);
  BFT_FREE(domain->properties);

  /* Free advection fields */
  if (domain->n_adv_fields > 0) {
    for (int i = 0; i < domain->n_adv_fields; i++)
      domain->adv_fields[i] = cs_advection_field_free(domain->adv_fields[i]);
    BFT_FREE(domain->adv_fields);
  }

  /* Print monitoring information */
  cs_log_printf(CS_LOG_PERFORMANCE,
                "%-36s %9s %9s %9s %9s %9s %9s\n",
                " ", "SysBuild", "Diffusion", "Advection", "Reaction",
                "Source", "Extra");
  for (int i = 0; i < domain->n_equations; i++)
    cs_equation_print_monitoring(domain->equations[i]);

  /* Free memory related to equations */
  for (int i = 0; i < domain->n_equations; i++)
    domain->equations[i] = cs_equation_free(domain->equations[i]);

  BFT_FREE(domain->equations);

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
 * \brief  Add a boundary type defined on a mesh location
 *
 * \param[in, out]   domain       pointer to a cs_domain_t structure
 * \param[in]        type         type of boundary to set
 * \param[in]        ml_name      mesh location name
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_set_default_boundary(cs_domain_t                *domain,
                               cs_param_boundary_type_t    type)
{
  if (domain == NULL) bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  if (type == CS_PARAM_BOUNDARY_WALL || type == CS_PARAM_BOUNDARY_SYMMETRY)
    domain->boundaries->default_boundary = type;
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
 * \param[in, out]   domain       pointer to a cs_domain_t structure
 * \param[in]        type         type of boundary to set
 * \param[in]        ml_name      mesh location name
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_add_boundary(cs_domain_t                *domain,
                       cs_param_boundary_type_t    type,
                       const char                 *ml_name)
{
  if (domain == NULL) bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  int  ml_id = cs_mesh_location_get_id_by_name(ml_name);

  if (ml_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid mesh location name %s.\n"
                " This mesh location is not already defined.\n"), ml_name);

  /* Sanity check */
  assert(cs_mesh_location_get_type(ml_id) == CS_MESH_LOCATION_BOUNDARY_FACES);

  /* Add this mesh location id to a list of mesh location ids of
     the same type */
  cs_domain_boundary_t  *dby = domain->boundaries;

  /* Number of mesh locations already defined for this type of boundary */
  int  n_ml_ids = dby->n_sub_ml_ids[type];

  /* Add a new mesh location for this type of boundary */
  BFT_REALLOC(dby->sub_ml_id_lst[type], n_ml_ids + 1, int);
  dby->sub_ml_id_lst[type][n_ml_ids] = ml_id;
  dby->n_sub_ml_ids[type] += 1;
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
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 * \param[in]        func      pointer to a cs_timestep_func_t function
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_def_time_step_by_function(cs_domain_t          *domain,
                                    cs_timestep_func_t   *func)
{
  if (domain == NULL) bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  domain->time_step->is_variable = 1; // not constant time step
  domain->time_options.idtvar = 1;    /* uniform in space but can change
                                         from one time step to the other */

  domain->time_step_def_type = CS_PARAM_DEF_BY_TIME_FUNCTION;
  domain->time_step_def.time_func = func;

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

  domain->time_step_def_type = CS_PARAM_DEF_BY_VALUE;
  domain->time_step->is_variable = 0; // constant time step
  domain->time_options.idtvar = 0;    // constant time step by default

  domain->time_step_def.get.val = dt;
  domain->dt_cur = domain->time_step_def.get.val;
  domain->time_options.dtref = domain->dt_cur;
  domain->time_options.dtmin = domain->dt_cur;
  domain->time_options.dtmax = domain->dt_cur;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new property to the current computational domain
 *
 * \param[in, out]  domain        pointer to a cs_domain_t structure
 * \param[in]       pty_name      name of the property to add
 * \param[in]       type_name     key name related to the type of property
 * \param[in]       n_subdomains  specify a definition in n_subdomains
 *
 * \return a pointer to the new cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_domain_add_property(cs_domain_t     *domain,
                       const char      *pty_name,
                       const char      *type_name,
                       int              n_subdomains)
{
  if (domain == NULL)
    return NULL;

  cs_property_t  *pty = cs_domain_get_property(domain, pty_name);

  if (pty != NULL) {
    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_DEFAULT,
                  _(" An existing property has already the name %s.\n"
                    " Stop adding this property.\n"), pty_name);
    return NULL;
  }

  int  pty_id = domain->n_properties;

  domain->n_properties += 1;
  BFT_REALLOC(domain->properties, domain->n_properties, cs_property_t *);

  pty = cs_property_create(pty_name, type_name, n_subdomains);
  domain->properties[pty_id] = pty;

  return pty;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new advection field to the current computational domain
 *
 * \param[in, out]   domain       pointer to a cs_domain_t structure
 * \param[in]        adv_name     name of the advection field to add
 *
 * \return a pointer to the new cs_adv_field_t structure
 */
/*----------------------------------------------------------------------------*/

cs_adv_field_t *
cs_domain_add_advection_field(cs_domain_t     *domain,
                              const char      *adv_name)
{
  if (domain == NULL)
    return NULL;

  cs_adv_field_t  *adv = cs_domain_get_advection_field(domain, adv_name);

  if (adv != NULL) {
    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_DEFAULT,
                  _(" An existing advection field has already the name %s.\n"
                    " Stop adding this advection field.\n"), adv_name);
    return adv;
  }

  int  adv_id = domain->n_adv_fields;

  domain->n_adv_fields += 1;
  BFT_REALLOC(domain->adv_fields, domain->n_adv_fields, cs_adv_field_t *);

  adv = cs_advection_field_create(adv_name);
  domain->adv_fields[adv_id] = adv;

  return adv;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate the computation of the wall distance
 *         Define a new equation called "WallDistance" and its related field
 *         named "WallDistance"
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 *
 * \return a pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_domain_activate_wall_distance(cs_domain_t   *domain)
{
  if (domain == NULL) bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  domain->wall_distance_eq_id = domain->n_equations;
  domain->n_predef_equations += 1;
  domain->n_equations += 1;
  BFT_REALLOC(domain->equations, domain->n_equations, cs_equation_t *);

  cs_equation_t  *eq =
    cs_equation_create("WallDistance",              // equation name
                       "WallDistance",              // variable name
                       CS_EQUATION_TYPE_PREDEFINED, // type of equation
                       CS_PARAM_VAR_SCAL,           // type of variable
                       CS_PARAM_BC_HMG_NEUMANN);    // default BC

  domain->equations[domain->wall_distance_eq_id] = eq;

  return eq;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate the computation of the Richards' equation
 *         This activation yields severals consequences:
 *         * Add a new equation named "Richards" along with an associated field
 *           named "hydraulic_head". Default boundary condition is set to
 *          "zero_flux".
 *         * Define a new advection field named "darcian_flux"
 *         * Define a new property called "permeability".
 *         * Define a new property called "soil_capacity" if "unsteady" is
 *           chosen
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 * \param[in]       kw_type    "isotropic", "orthotropic or "anisotropic"
 * \param[in]       kw_time    Richards equation is "steady" or "unsteady"
 * \param[in]       n_soils    number of soils to consider
 * \param[in]       n_tracers  number of tracer equations
 *
 * \return a pointer to a cs_grounwater_t structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_t *
cs_domain_activate_gwf(cs_domain_t   *domain,
                       const char    *kw_type,
                       const char    *kw_time,
                       int            n_soils,
                       int            n_tracers)
{
  if (domain == NULL) bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  int  richards_eq_id = domain->n_equations;

  /* Allocate a new structure for managing groundwater module */
  domain->gwf = cs_gwf_create();

  /* Add a property related to the diffusion term of the Richards eq. */
  cs_property_t  *permeability = cs_domain_add_property(domain,
                                                        "permeability",
                                                        kw_type,
                                                        n_soils);

  /* Add a property related to the unsteady term of the Richards eq. */
  cs_property_t  *soil_capacity = NULL;
  if (strcmp(kw_time, "unsteady") == 0)
    soil_capacity = cs_domain_add_property(domain,
                                           "soil_capacity",
                                           "isotropic",
                                           n_soils);

  /* Add an advection field related to the darcian flux stemming from the
     Richards equation */
  cs_adv_field_t  *adv_field = cs_domain_add_advection_field(domain,
                                                             "darcian_flux");

  cs_advection_field_set_option(adv_field, CS_ADVKEY_DEFINE_AT, "cells");

  /* Create a new equation for solving the Richards equation */
  cs_equation_t  *richards_eq = cs_gwf_initialize(richards_eq_id,
                                                  n_soils,
                                                  n_tracers,
                                                  permeability,
                                                  soil_capacity,
                                                  adv_field,
                                                  domain->gwf);

  if (richards_eq == NULL) {
    bft_error(__FILE__, __LINE__, 0,
              " The module dedicated to groundwater flows is activated but"
              " the Richards' equation is not set.\n"
              " Please check your settings.");
    return NULL; // Avoid a warning
  }

  /* Update cs_domain_t structure */
  domain->richards_eq_id = richards_eq_id;
  domain->n_predef_equations += 1;
  domain->n_equations += 1;
  BFT_REALLOC(domain->equations, domain->n_equations, cs_equation_t *);
  domain->equations[richards_eq_id] = richards_eq;

  return domain->gwf;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new equation related to the groundwater flow module
 *         This equation is a particular type of unsteady advection-diffusion
 *         reaction eq.
 *         Tracer is advected thanks to the darcian velocity and
 *         diffusion/reaction parameters result from a physical modelling.
 *
 * \param[in, out]  domain         pointer to a cs_domain_t structure
 * \param[in]       eqname         name of the equation
 * \param[in]       varname        name of the related variable
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_add_gwf_tracer_eq(cs_domain_t   *domain,
                            const char    *eq_name,
                            const char    *var_name)
{
  /* Sanity checks */
  if (domain == NULL) bft_error(__FILE__, __LINE__, 0, _err_empty_domain);
  if (domain->gwf == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" Groundwater module is requested but is not activated.\n"
                " Please first activate this module."));

  /* Add a new equation */
  BFT_REALLOC(domain->equations, domain->n_equations + 1, cs_equation_t *);

  cs_equation_t  *tracer_eq = cs_gwf_add_tracer(domain->gwf,
                                                domain->n_equations,
                                                eq_name,
                                                var_name);

  if (tracer_eq == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " Problem during the definition of the tracer equation %s for"
              " the groundwater module.", eq_name);

  /* Add a new property related to the time-depedent term */
  char  *pty_name = NULL;
  int  len = strlen(eq_name) + strlen("_time") + 1;
  BFT_MALLOC(pty_name, len, char);
  sprintf(pty_name, "%s_time", eq_name);

  const int  n_soils = cs_gwf_get_n_soils(domain->gwf);
  cs_property_t  *time_pty = cs_domain_add_property(domain,
                                                    pty_name,
                                                    "isotropic",
                                                    n_soils);

  cs_equation_link(tracer_eq, "time", time_pty);

  domain->equations[domain->n_equations] = tracer_eq;
  domain->n_predef_equations += 1;
  domain->n_equations += 1;

  BFT_FREE(pty_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the parameters related to a tracer equation used in the
 *         groundwater flow module
 *
 * \param[in, out]  domain         pointer to a cs_domain_t structure
 * \param[in]       eqname         name of the equation
 * \param[in]       ml_name        name of the related mesh location
 * \param[in]       wmd            value of the water molecular diffusivity
 * \param[in]       alpha_l        value of the longitudinal dispersivity
 * \param[in]       alpha_t        value of the transversal dispersivity
 * \param[in]       distrib_coef   value of the distribution coefficient
 * \param[in]       reaction_rate  value of the first order rate of reaction
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_set_gwf_tracer_eq(cs_domain_t   *domain,
                            const char    *eq_name,
                            const char    *ml_name,
                            double         wmd,
                            double         alpha_l,
                            double         alpha_t,
                            double         distrib_coef,
                            double         reaction_rate)
{
  /* Sanity checks */
  if (domain == NULL) bft_error(__FILE__, __LINE__, 0, _err_empty_domain);
  if (domain->gwf == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" Groundwater module is requested but is not activated.\n"
                " Please first activate this module."));

  int  eq_id = -1;
  bool found = false;
  for (int i = 0; i < domain->n_equations && found == false; i++) {
    cs_equation_t  *_eq = domain->equations[i];
    if (strcmp(eq_name, cs_equation_get_name(_eq)) == 0)
      eq_id = i, found = true;
  }

  if (found == false)
    bft_error(__FILE__, __LINE__, 0,
              _(" Stop setting a tracer equation %s.\n"
                " No equation with this name has been found.\n"
                " Please check your settings."), eq_name);

  cs_gwf_set_tracer_param(domain->gwf, eq_id, ml_name,
                          wmd,
                          alpha_l, alpha_t,
                          distrib_coef,
                          reaction_rate);
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
  cs_domain_boundary_t  *dby = domain->boundaries;
  cs_param_boundary_type_t  default_type = dby->default_boundary;

  /* Store all mesh locations not associated to the default type */
  int  n_sub_ids = 0;
  int  *sub_ids = NULL;

  for (unsigned int type = 0; type < CS_PARAM_N_BOUNDARY_TYPES; type++) {

    int  _n_sub_ids = dby->n_sub_ml_ids[type];
    int  *_sub_ids = dby->sub_ml_id_lst[type];

    if (_n_sub_ids > 0 && type != default_type) {

      // There is at least one mesh location
      assert(_sub_ids != NULL);

      dby->autogen_ml_ids[type] =
        cs_mesh_location_add_by_union(_domain_boundary_ml_name[type],
                                      CS_MESH_LOCATION_BOUNDARY_FACES,
                                      _n_sub_ids,
                                      _sub_ids,
                                      false);  // complement ?

      /* Increment the list of mesh locations */
      BFT_REALLOC(sub_ids, n_sub_ids + _n_sub_ids, int);
      for (int i = 0; i < _n_sub_ids; i++)
        sub_ids[n_sub_ids + i] = _sub_ids[i];
      n_sub_ids += _n_sub_ids;

    }

  } /* Loop on boundary types */

  /* Treatment of the default boundary type */
  dby->autogen_ml_ids[default_type] =
    cs_mesh_location_add_by_union(_domain_boundary_ml_name[default_type],
                                  CS_MESH_LOCATION_BOUNDARY_FACES,
                                  n_sub_ids,
                                  sub_ids,
                                  true);

  BFT_FREE(sub_ids);
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
  assert(domain != NULL);

  cs_domain_boundary_t  *bdy = domain->boundaries;
  cs_equation_t  *eq = NULL;

  /* Wall distance */
  /* ------------- */

  if (domain->wall_distance_eq_id > -1) {

    eq = domain->equations[domain->wall_distance_eq_id];

    cs_walldistance_setup(eq,
                          cs_domain_get_property(domain, "unity"),
                          bdy->autogen_ml_ids[CS_PARAM_BOUNDARY_WALL]);

  } // Wall distance is activated

  /* Groundwater flow module */
  /* ----------------------- */

  if (domain->richards_eq_id > -1) {

    int  len = 0, max_len = 0;
    char  *pty_name = NULL;

    /* Automatic settings for the Richards equation */
    eq = domain->equations[domain->richards_eq_id];
    cs_gwf_richards_setup(domain->gwf, eq);

    /* Automatic settings for the tracer equations */
    for (int eq_id = 0; eq_id < domain->n_equations; eq_id++) {

      if (eq_id != domain->richards_eq_id) {

        eq = domain->equations[eq_id];

        if (cs_equation_get_type(eq) == CS_EQUATION_TYPE_GROUNDWATER) {

          const int  n_soils = cs_gwf_get_n_soils(domain->gwf);

          if (cs_gwf_tracer_needs_diffusion(domain->gwf, eq_id)) {

            /* Add a new property related to the diffusion property */
            const char *eq_name = cs_equation_get_name(eq);

            len = strlen(eq_name) + strlen("_diffusivity") + 1;
            if (len > max_len)
              max_len = len, BFT_REALLOC(pty_name, len, char);
            sprintf(pty_name, "%s_diffusivity", eq_name);

            cs_property_t *diff_pty = cs_domain_add_property(domain,
                                                             pty_name,
                                                             "anisotropic",
                                                             n_soils);

            cs_equation_link(eq, "diffusion", diff_pty);

          } /* Add a diffusion property for this equation ? */

          if (cs_gwf_tracer_needs_reaction(domain->gwf, eq_id)) {

            /* Add a new property related to the reaction property */
            const char *eq_name = cs_equation_get_name(eq);

            len = strlen(eq_name) + strlen("_reaction") + 1;
            if (len > max_len)
              max_len = len, BFT_REALLOC(pty_name, len, char);
            sprintf(pty_name, "%s_reaction", eq_name);

            cs_property_t *r_pty = cs_domain_add_property(domain,
                                                          pty_name,
                                                          "isotropic",
                                                          n_soils);

            cs_equation_add_linear_reaction(eq, r_pty, "decay");

          } /* Add a reaction property for this equation ? */

          /* Automatic settings for this tracer equation */
          cs_gwf_tracer_setup(eq_id, eq, domain->gwf);

        } /* Tracer equation related to the groundwater flow module */

      } /* Not the Richards equation */

    } /* Loop on equations */

    BFT_FREE(pty_name);

  } /* Groundwater flow module is activated */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new user equation to a domain
 *
 * \param[in, out] domain         pointer to a cs_domain_t structure
 * \param[in]      eqname         name of the equation
 * \param[in]      varname        name of the related variable
 * \param[in]      key_type       type of equation: "scalar", "vector", "tensor"
 * \param[in]      key_bc         type of boundary condition set by default
 *                                "zero_value" or "zero_flux"
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_add_user_equation(cs_domain_t         *domain,
                            const char          *eqname,
                            const char          *varname,
                            const char          *key_type,
                            const char          *key_bc)
{
  /* Sanity checks */
  if (domain == NULL) bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  cs_param_var_type_t  var_type = CS_PARAM_N_VAR_TYPES;
  cs_param_bc_type_t  default_bc = CS_PARAM_N_BC_TYPES;

  BFT_REALLOC(domain->equations, domain->n_equations + 1, cs_equation_t *);

  /* Define the type of equation */
  if (strcmp(key_type, "scalar") == 0)
    var_type = CS_PARAM_VAR_SCAL;
  else if (strcmp(key_type, "vector") == 0)
    var_type = CS_PARAM_VAR_VECT;
  else if (strcmp(key_type, "tensor") == 0)
    var_type = CS_PARAM_VAR_TENS;
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of equation: %s\n"
                " Choices are scalar, vector or tensor."), key_type);

  /* Define a boundary condition by default */
  if (strcmp(key_bc, "zero_value") == 0)
    default_bc = CS_PARAM_BC_HMG_DIRICHLET;
  else if (strcmp(key_bc, "zero_flux") == 0)
    default_bc = CS_PARAM_BC_HMG_NEUMANN;
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of boundary condition by default: %s\n"
                " Choices are zero_value or zero_flux."), key_bc);

  domain->equations[domain->n_equations] =
    cs_equation_create(eqname,                // equation name
                       varname,               // variable name
                       CS_EQUATION_TYPE_USER, // type of equation
                       var_type,              // type of equation
                       default_bc);           // default BC

  domain->n_user_equations += 1;
  domain->n_equations += 1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Find the related property definition from its name
 *
 * \param[in]  domain      pointer to a domain structure
 * \param[in]  ref_name    name of the property to find
 *
 * \return NULL if not found otherwise the associated pointer
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_domain_get_property(const cs_domain_t    *domain,
                       const char           *ref_name)
{
  for (int i = 0; i < domain->n_properties; i++) {

    cs_property_t  *pty = domain->properties[i];
    if (cs_property_check_name(pty, ref_name))
      return pty;

  }

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Find the related advection field definition from its name
 *
 * \param[in]  domain      pointer to a domain structure
 * \param[in]  ref_name    name of the adv_field to find
 *
 * \return NULL if not found otherwise the associated pointer
 */
/*----------------------------------------------------------------------------*/

cs_adv_field_t *
cs_domain_get_advection_field(const cs_domain_t    *domain,
                              const char           *ref_name)
{
  for (int i = 0; i < domain->n_adv_fields; i++) {

    cs_adv_field_t  *adv = domain->adv_fields[i];
    if (cs_advection_field_check_name(adv, ref_name))
      return adv;

  }

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Find the cs_equation_t structure whith name eqname
 *         Return NULL if not find
 *
 * \param[in]  domain    pointer to a cs_domain_t structure
 * \param[in]  eqname    name of the equation to find
 *
 * \return a pointer to a cs_equation_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_domain_get_equation(const cs_domain_t  *domain,
                       const char         *eqname)
{
  cs_equation_t  *eq = NULL;

  for (int i = 0; i < domain->n_equations; i++) {

    cs_equation_t  *_eq = domain->equations[i];
    if (strcmp(eqname, cs_equation_get_name(_eq)) == 0) {
      eq = _eq;
      break;
    }

  }

  return eq;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the pointer to a cs_gwf_t structure related to this
 *         domain
 *
 * \param[in]   domain         pointer to a cs_domain_t structure
 *
 * \return a pointer to a cs_gwf_t structure
 */
/*----------------------------------------------------------------------------*/

cs_gwf_t *
cs_domain_get_gwf_struct(const cs_domain_t    *domain)
{
  if (domain == NULL)
    return NULL;

  return domain->gwf;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create cs_field_t structures attached to equation unknowns and
 *         advection fields
 *
 * \param[in, out]  domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_create_fields(cs_domain_t  *domain)
{
  /* Loop on all equations */
  for (int eq_id = 0; eq_id < domain->n_equations; eq_id++)
    cs_equation_create_field(domain->equations[eq_id]);

  /* Loop on all advection fields */
  for (int adv_id = 0; adv_id < domain->n_adv_fields; adv_id++)
    cs_advection_field_create_field(domain->adv_fields[adv_id]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the scheme flag for the current computational domain
 *
 * \param[in, out]  domain            pointer to a cs_domain_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_set_scheme_flag(cs_domain_t                 *domain)
{
  if (domain == NULL) bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  /* Define a scheme flag for the current domain */
  for (int eq_id = 0; eq_id < domain->n_equations; eq_id++) {

    cs_equation_t  *eq = domain->equations[eq_id];
    cs_space_scheme_t  scheme = cs_equation_get_space_scheme(eq);
    cs_param_var_type_t  vartype = cs_equation_get_var_type(eq);

    if (vartype == CS_PARAM_VAR_SCAL)
      domain->scheme_flag |= CS_SCHEME_FLAG_SCALAR;
    else if (vartype == CS_PARAM_VAR_VECT)
      domain->scheme_flag |= CS_SCHEME_FLAG_VECTOR;

    if (scheme == CS_SPACE_SCHEME_CDOVB)
      domain->scheme_flag |= CS_SCHEME_FLAG_CDOVB | CS_SCHEME_FLAG_POLY0;
    else if (scheme == CS_SPACE_SCHEME_CDOVCB)
      domain->scheme_flag |= CS_SCHEME_FLAG_CDOVCB | CS_SCHEME_FLAG_POLY0;
    else if (scheme == CS_SPACE_SCHEME_CDOFB)
      domain->scheme_flag |= CS_SCHEME_FLAG_CDOFB | CS_SCHEME_FLAG_POLY0;
    else if (scheme == CS_SPACE_SCHEME_HHO)
      domain->scheme_flag |= CS_SCHEME_FLAG_HHO;
    else
      bft_error(__FILE__, __LINE__, 0,
                _(" Undefined type of equation to solve for eq. %s."
                  " Please check your settings."), cs_equation_get_name(eq));

    switch (cs_equation_get_space_poly_degree(eq)) {
    case 0:
      domain->scheme_flag |= CS_SCHEME_FLAG_POLY0;
      break;
    case 1:
      domain->scheme_flag |= CS_SCHEME_FLAG_POLY1;
      break;
    case 2:
      domain->scheme_flag |= CS_SCHEME_FLAG_POLY2;
      break;
    default:
      bft_error(__FILE__, __LINE__, 0, " Undefined space polynomial order.");
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
cs_domain_initialize(cs_domain_t                 *domain,
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

  /* Allocate all fields created during the setup stage */
  cs_field_allocate_or_map_all();

  /* Specify the "physical" domain boundaries. Define a mesh location for
     each boundary type */
  _build_domain_boundaries(mesh->n_b_faces, domain->boundaries);

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

  /* Last stage to define properties (when complex definition is requested) */
  for (int i = 0; i < domain->n_properties; i++)
    cs_property_last_definition_stage(domain->properties[i]);

  /* Initialization default post-processing for the computational domain */
  cs_domain_post_init(domain->dt_cur,
                      domain->cdo_quantities,
                      domain->n_adv_fields, domain->adv_fields,
                      domain->n_properties, domain->properties,
                      domain->n_equations, domain->equations);

  /* Allocate common structures for solving equations */
  cs_equation_allocate_common_structures(domain->connect,
                                         domain->cdo_quantities,
                                         domain->time_step,
                                         domain->scheme_flag);

  if (domain->richards_eq_id > -1)
    cs_gwf_final_initialization(domain->connect,
                                domain->n_equations,
                                domain->equations,
                                domain->gwf);

  /* Proceed to the last settings of a cs_equation_t structure
     - Assign to a cs_equation_t structure a list of function to manage this
       structure during the computation.
     - The set of functions chosen for each equation depends on the parameters
       specifying the cs_equation_t structure
     - Setup the structure related to cs_sles_*
  */

  for (int eq_id = 0; eq_id < domain->n_equations; eq_id++) {

    cs_equation_t  *eq = domain->equations[eq_id];

    if (domain->profiling)
      cs_equation_set_timer_stats(eq);

    cs_equation_last_setup(domain->connect, eq);

    if (!cs_equation_is_steady(eq))
      domain->only_steady = false;

  } // Loop on domain equations

  if (domain->only_steady)
    domain->is_last_iter = true;
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
cs_domain_needs_iterate(cs_domain_t  *domain)
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

  cs_param_def_type_t  def_type = domain->time_step_def_type;

  if (def_type != CS_PARAM_DEF_BY_VALUE) {

    if (def_type == CS_PARAM_DEF_BY_TIME_FUNCTION) {

      domain->dt_cur = domain->time_step_def.time_func(nt_cur, t_cur);

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
    if (domain->wall_distance_eq_id > -1) {

      cs_equation_t  *walld_eq = domain->equations[domain->wall_distance_eq_id];

      /* Compute the wall distance */
      cs_walldistance_compute(domain->mesh,
                              domain->time_step,
                              domain->dt_cur,
                              domain->connect,
                              domain->cdo_quantities,
                              walld_eq);

    } // wall distance

    /* If unsteady, only initialization is done, otherwise one makes the whole
       computation */
    if (domain->richards_eq_id > -1)
      cs_gwf_compute(domain->mesh,
                     domain->time_step,
                     domain->dt_cur,
                     domain->connect,
                     domain->cdo_quantities,
                     domain->equations,
                     domain->gwf);

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

    if (domain->richards_eq_id > -1)
      cs_gwf_compute(domain->mesh,
                     domain->time_step,
                     domain->dt_cur,
                     domain->connect,
                     domain->cdo_quantities,
                     domain->equations,
                     domain->gwf);

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
  if (domain->force_advfield_update) {
    for (int adv_id = 0; adv_id < domain->n_adv_fields; adv_id++)
      cs_advection_field_update(domain->adv_fields[adv_id]);
  }

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
 * \brief  Write a restart file for the CDO module
 *
 * \param[in]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_write_restart(const cs_domain_t  *domain)
{
  cs_restart_t  *restart = cs_restart_create("cdo", // restart file name
                                             NULL,  // directory name
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
  cs_restart_write_section(restart,
                           "cdo:n_equations",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_cs_int_t,
                           &(domain->n_equations));

  /* Write a new section */
  cs_restart_write_section(restart,
                           "cdo:n_properties",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_cs_int_t,
                           &(domain->n_properties));

  /* Write a new section */
  cs_restart_write_section(restart,
                           "cdo:n_adv_fields",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_cs_int_t,
                           &(domain->n_adv_fields));

  /* Write a new section: activation or not of the groundwater flow module */
  int  igwf = (domain->gwf != NULL) ? 1 : 0;
  cs_restart_write_section(restart,
                           "groundwater_flow_module",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_cs_int_t,
                           &igwf);

  /* Write a new section: computation or not of the wall distance */
  int  iwall = (domain->wall_distance_eq_id != -1) ? 1 : 0;
  cs_restart_write_section(restart,
                           "wall_distance",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_cs_int_t,
                           &iwall);

  /* Write a new section: number of computed time steps */
  int  ntcabs = domain->time_step->nt_cur;
  cs_restart_write_section(restart,
                           "nbre_pas_de_temps",
                           CS_MESH_LOCATION_NONE,
                           1,
                           CS_TYPE_cs_int_t,
                           &ntcabs);

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

  cs_log_printf(CS_LOG_SETUP, " -msg- n_cdo_equations          %d\n",
                domain->n_equations);
  cs_log_printf(CS_LOG_SETUP, " -msg- n_predefined_equations   %d\n",
                domain->n_predef_equations);
  cs_log_printf(CS_LOG_SETUP, " -msg- n_user_equations         %d\n",
                domain->n_user_equations);
  cs_log_printf(CS_LOG_SETUP, " -msg- n_properties             %d\n",
                domain->n_properties);

  /* Boundary */
  cs_domain_boundary_t  *bdy = domain->boundaries;

  cs_log_printf(CS_LOG_SETUP, "\n  Domain boundary by default: ");
  switch (bdy->default_boundary) {
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
  cs_gnum_t  counter[4] = { bdy->n_type_elts[CS_PARAM_BOUNDARY_WALL],
                            bdy->n_type_elts[CS_PARAM_BOUNDARY_INLET],
                            bdy->n_type_elts[CS_PARAM_BOUNDARY_OUTLET],
                            bdy->n_type_elts[CS_PARAM_BOUNDARY_SYMMETRY] };
  if (cs_glob_n_ranks > 1)
    cs_parall_counter(counter, 4);

  cs_log_printf(CS_LOG_SETUP,
                "  >> Number of faces with a wall boundary:      %lu\n",
                counter[0]);
  cs_log_printf(CS_LOG_SETUP,
                "  >> Number of faces with an inlet boundary:    %lu\n",
                counter[1]);
  cs_log_printf(CS_LOG_SETUP,
                "  >> Number of faces with an outlet boundary:   %lu\n",
                counter[2]);
  cs_log_printf(CS_LOG_SETUP,
                "  >> Number of faces with a symmetry boundary:  %lu\n",
                counter[3]);

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
    cs_log_printf(CS_LOG_SETUP, "  >> Type of definition: %s",
                  cs_param_get_def_type_name(domain->time_step_def_type));
    if (domain->time_step_def_type == CS_PARAM_DEF_BY_VALUE)
      cs_log_printf(CS_LOG_SETUP, " => %5.3e\n", domain->dt_cur);
    else
      cs_log_printf(CS_LOG_SETUP, "\n");

    if (domain->time_step->t_max > 0.)
      cs_log_printf(CS_LOG_SETUP, "%-30s %5.3e\n",
                    "  >> Final simulation time:", domain->time_step->t_max);
    if (domain->time_step->nt_max > 0)
      cs_log_printf(CS_LOG_SETUP, "%-30s %9d\n",
                    "  >> Final time step:", domain->time_step->nt_max);

  }
  cs_log_printf(CS_LOG_SETUP, "\n");

  if (domain->verbosity > 0) {

    /* Properties */
    cs_log_printf(CS_LOG_SETUP, "\n%s", lsepline);
    cs_log_printf(CS_LOG_SETUP, "\tSummary of the definition of properties\n");
    cs_log_printf(CS_LOG_SETUP, "%s", lsepline);

    for (int i = 0; i < domain->n_properties; i++)
      cs_property_summary(domain->properties[i]);

    /* Advection fields */
    if (domain->n_adv_fields > 0) {

      cs_log_printf(CS_LOG_SETUP, "\n%s", lsepline);
      cs_log_printf(CS_LOG_SETUP, "\tSummary of the advection field\n");
      cs_log_printf(CS_LOG_SETUP, "%s", lsepline);

      for (int i = 0; i < domain->n_adv_fields; i++)
        cs_advection_field_summary(domain->adv_fields[i]);

    }

    /* Summary of the groundwater module */
    cs_gwf_summary(domain->gwf);

    /* Summary for each equation */
    for (int  eq_id = 0; eq_id < domain->n_equations; eq_id++)
      cs_equation_summary(domain->equations[eq_id]);

  } /* Domain->verbosity > 0 */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
