#ifndef __CS_DOMAIN_H__
#define __CS_DOMAIN_H__

/*============================================================================
 * Manage a computational domain
 *  - Physical boundary conditions attached to a domain
 *  - Equations
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

#include "cs_time_step.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"
#include "cs_param.h"
#include "cs_equation.h"
#include "cs_property.h"
#include "cs_advection_field.h"
#include "cs_groundwater.h"

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

  /* Flag to know if scalar or vector equations are requested and which kind
     of numerical schemes is requested to solve these equations */
  cs_flag_t        scheme_flag;


  /* Pre-defined equations to solve
     If xxxxx_eq_id = -1, then this equation is not activated */

  /* GROUNDWATER FLOW MODULE */

  int   richards_eq_id;   // Main equation of the groundwater module
  cs_groundwater_t  *gw;  // NULL is not activated

  /* WALL DISTANCE */
  int   wall_distance_eq_id;  // Wall distance computation

  /* TODO: NAVIER-STOKES */

  /* Output options */
  int        output_nt;  /* Log information every nt iteration(s) */
  cs_real_t  output_dt;  /* Log information every dt elapsed simulated time */
  int        verbosity;  /* Level of details given in log */
  bool       profiling;  /* Activate a set of timer statistics (details differ
                            according to the verbosity level) */

} cs_domain_t;

/* List of available keys for setting a property */
typedef enum {

  CS_DOMAIN_DEFAULT_BOUNDARY, // set the default boundary
  CS_DOMAIN_OUTPUT_NT,        // set the output frequency (number of iterations)
  CS_DOMAIN_OUTPUT_DT,        // set the output frequency (simulation time)
  CS_DOMAIN_PROFILING,        // activate the profiling
  CS_DOMAIN_NTMAX,            // set the max number of iterations to perform
  CS_DOMAIN_TMAX,             // set the maximum simulated time
  CS_DOMAIN_VERBOSITY,        // set the verbosity
  CS_DOMAIN_N_KEYS

} cs_domain_key_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize of cs_domain_t structure
 *
 * \param[in]   mesh              pointer to a cs_mesh_t struct.
 * \param[in]   mesh_quantities   pointer to a cs_mesh_quantities_t struct.
 *
 * \return a pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

cs_domain_t *
cs_domain_init(const cs_mesh_t             *mesh,
               const cs_mesh_quantities_t  *mesh_quantities);

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
 * \brief  Set auxiliary parameters related to a cs_domain_t structure
 *
 * \param[in, out]  domain    pointer to a cs_domain_t structure
 * \param[in]       key       key related to the parameter to set
 * \param[in]       keyval    value related to the parameter to set
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_set_param(cs_domain_t       *domain,
                    cs_domain_key_t    key,
                    const char        *keyval);

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
 * \brief  Summary of a cs_domain_t structure
 *
 * \param[in]   domain    pointer to the cs_domain_t structure to summarize
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_summary(const cs_domain_t   *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Proceed to the last settings of a cs_domain_t structure
 *
 * \param[in, out]  domain    pointer to the cs_domain_t structure to set
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_last_setup(cs_domain_t    *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a boundary type defined on a mesh location
 *
 * \param[in, out]   domain       pointer to a cs_domain_t structure
 * \param[in]        ml_name      mesh location name
 * \param[in]        bdy_name     key name of boundary to set
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_add_boundary(cs_domain_t               *domain,
                       const char                *ml_name,
                       const char                *bdy_name);

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
 * \brief  Add a new property to the current computational domain
 *
 * \param[in, out]  domain        pointer to a cs_domain_t structure
 * \param[in]       pty_name      name of the property to add
 * \param[in]       type_name     key name related to the type of property
 * \param[in]       n_subdomains  specify a definition in n_subdomains
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_add_property(cs_domain_t     *domain,
                       const char      *pty_name,
                       const char      *type_name,
                       int              n_subdomains);

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
 * \brief  Add a new advection field to the current computational domain
 *
 * \param[in, out]   domain       pointer to a cs_domain_t structure
 * \param[in]        adv_name     name of the advection field to add
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_add_advection_field(cs_domain_t     *domain,
                              const char      *adv_name);

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
 * \brief  Activate the computation of the wall distance
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_activate_wall_distance(cs_domain_t   *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate the computation of the Richards' equation
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 * \param[in]       kw_type    "isotropic", "orthotropic or "anisotropic"
 * \param[in]       kw_time    Richards equation is "steady" or "unsteady"
 * \param[in]       n_soils    number of soils to consider
 * \param[in]       n_tracers  number of tracer equations
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_activate_groundwater(cs_domain_t   *domain,
                               const char    *kw_type,
                               const char    *kw_time,
                               int            n_soils,
                               int            n_tracers);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the pointer to a cs_groundwater_t structure related to this
 *         domain
 *
 * \param[in]   domain         pointer to a cs_domain_t structure
 *
 * \return a pointer to a cs_groundwater_t structure
 */
/*----------------------------------------------------------------------------*/

cs_groundwater_t *
cs_domain_get_groundwater(const cs_domain_t    *domain);

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
cs_domain_add_groundwater_tracer(cs_domain_t   *domain,
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
 * \param[in]       bulk_density   value of the bulk density
 * \param[in]       distrib_coef   value of the distribution coefficient
 * \param[in]       reaction_rate  value of the first order rate of reaction
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_set_groundwater_tracer(cs_domain_t   *domain,
                                 const char    *eq_name,
                                 const char    *ml_name,
                                 double         wmd,
                                 double         alpha_l,
                                 double         alpha_t,
                                 double         bulk_density,
                                 double         distrib_coef,
                                 double         reaction_rate);

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
 * \brief  Solve all the equations of a computational domain for one time step
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_solve(cs_domain_t  *domain);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Process the computed solution
 *
 * \param[in]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_postprocess(cs_domain_t  *domain);

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

END_C_DECLS

#endif /* __CS_DOMAIN_H__ */
