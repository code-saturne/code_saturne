#ifndef __CS_DOMAIN_H__
#define __CS_DOMAIN_H__

/*============================================================================
 * Manage a computational domain
 *  - Physical boundary conditions attached to a domain
 *  - Equations
 *============================================================================*/

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_advection_field.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"
#include "cs_equation.h"
#include "cs_gwf.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_param.h"
#include "cs_property.h"
#include "cs_time_step.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Store information about the kind of boundary attached to the computational
   domain: inlet, outlet, wall, symmetry... */

typedef struct _cs_domain_boundary_t cs_domain_boundary_t;

typedef struct {

  /* Code_Saturne mesh and mesh quantities structures already computed */
  const  cs_mesh_t              *mesh;
  const  cs_mesh_quantities_t   *mesh_quantities;

  /* CDO structures:
     - cs_cdo_connect_t contains additional information about connectivity
     - cs_cdo_quantities_t contains additional information on mesh quantities
  */
  cs_cdo_connect_t              *connect;
  cs_cdo_quantities_t           *cdo_quantities;

  /* Physical boundary conditions on the computational domain */
  cs_domain_boundary_t          *boundaries;

  /* Time step management */
  bool                     is_last_iter;       // true or false
  double                   dt_cur;             // current time step
  cs_param_def_type_t      time_step_def_type; // Way of defining the time step
  cs_def_t                 time_step_def;      // Definition of the time_step
  cs_time_step_t          *time_step;          // time step descriptor
  cs_time_step_options_t   time_options;       // time step options

  /* Properties attached to the computational domain.
     "unity" is created by default
  */
  int               n_properties;
  cs_property_t   **properties;

  /* Advection fields attached to the computational domain */
  int               n_adv_fields;
  cs_adv_field_t  **adv_fields;

  /* Number of equations defined on this domain splitted into
     predefined equations and user equations.
     Predefined equations are stored first.
  */
  int              n_equations;
  int              n_predef_equations;
  int              n_user_equations;
  cs_equation_t  **equations;

  bool             only_steady;
  bool             force_advfield_update;

  /* Flag to know if scalar or vector equations are requested and which kind
     of numerical schemes is requested to solve these equations */
  cs_flag_t        scheme_flag;

  /* Pre-defined equations to solve
     If xxxxx_eq_id = -1, then this equation is not activated */

  /* GROUNDWATER FLOW MODULE */
  int        richards_eq_id;   // Main equation of the groundwater flow module
  cs_gwf_t  *gwf;              // NULL is not activated

  /* WALL DISTANCE */
  int   wall_distance_eq_id;  // Wall distance computation

  /* TODO: NAVIER-STOKES */

  /* Output options */
  int        output_nt;  /* Log information every nt iteration(s) */
  int        verbosity;  /* Level of details given in log */
  bool       profiling;  /* Activate a set of timer statistics (details differ
                            according to the verbosity level) */

  /* Monitoring */
  cs_timer_counter_t    tcp; /* Cumulated elapsed time for extra-operations
                                and post-processing */
  cs_timer_counter_t    tcs; /* Cumulated elapsed time for setup operations */

} cs_domain_t;

/* List of available keys for setting a property */
typedef enum {

  CS_DOMAIN_PROFILING,        // activate the profiling
  CS_DOMAIN_N_KEYS

} cs_domain_key_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

extern cs_domain_t *cs_glob_domain; /* Pointer to main computational domain
                                       used in CDO/HHO schmes */

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
cs_domain_create(void);

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
cs_domain_free(cs_domain_t   *domain);

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
                               cs_param_boundary_type_t    type);

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
                       const char                 *ml_name);

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
                         double             t_max);

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
                           int                verbosity);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set to true the automatic update of all advection fields
 *
 * \param[in, out]  domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_update_advfield(cs_domain_t       *domain);

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
                             const char        *keyval);

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
                                    cs_timestep_func_t   *func);

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
                                 double         dt);

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
                       int              n_subdomains);

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
                              const char      *adv_name);

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
cs_domain_activate_wall_distance(cs_domain_t   *domain);

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
                       int            n_tracers);

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
                            const char    *var_name);

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
                            double         reaction_rate);

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
                            const char          *key_bc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add new mesh locations related to domain boundaries from existing
 *         mesh locations
 *
 * \param[in]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_update_mesh_locations(cs_domain_t   *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup predefined equations which are activated
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_setup_predefined_equations(cs_domain_t   *domain);

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
                       const char           *ref_name);

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
                              const char           *ref_name);

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
                       const char         *eqname);

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
cs_domain_get_gwf_struct(const cs_domain_t    *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a cs_field_t structure for each equation defined in the
 *         domain
 *
 * \param[in, out]  domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_create_fields(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a cs_domain_t structure
 *
 * \param[in, out]  domain            pointer to a cs_domain_t struct.
 * \param[in, out]  mesh              pointer to a cs_mesh_t struct.
 * \param[in]       mesh_quantities   pointer to a cs_mesh_quantities_t struct.
 *
 * \return a pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_initialize(cs_domain_t                 *domain,
                     cs_mesh_t                   *mesh,
                     const cs_mesh_quantities_t  *mesh_quantities);

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
cs_domain_needs_iterate(cs_domain_t  *domain);

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
cs_domain_needs_log(const cs_domain_t      *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the current time step for this new time iteration
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_define_current_time_step(cs_domain_t   *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update time step after one temporal iteration
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_increment_time(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve all the equations of a computational domain for one time step
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_solve(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Process the computational domain after the resolution
 *
 * \param[in]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_process_after_solve(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Write a restart file for the CDO module
 *
 * \param[in]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_write_restart(const cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of a cs_domain_t structure
 *
 * \param[in]   domain    pointer to the cs_domain_t structure to summarize
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_summary(const cs_domain_t   *domain);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_DOMAIN_H__ */
