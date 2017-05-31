#ifndef __CS_EQUATION_H__
#define __CS_EQUATION_H__

/*============================================================================
 * Routines to handle cs_equation_t structure and its related structures
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

#include "cs_cdo_quantities.h"
#include "cs_equation_param.h"
#include "cs_field.h"
#include "cs_param.h"
#include "cs_mesh.h"
#include "cs_source_term.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _cs_equation_t cs_equation_t;

/* List of available keys for setting an equation */
typedef enum {

  CS_EQKEY_ADV_FORMULATION,
  CS_EQKEY_ADV_SCHEME,
  CS_EQKEY_ADV_FLUX_QUADRA,
  CS_EQKEY_BC_ENFORCEMENT,
  CS_EQKEY_BC_QUADRATURE,
  CS_EQKEY_EXTRA_OP,
  CS_EQKEY_HODGE_DIFF_ALGO,
  CS_EQKEY_HODGE_DIFF_COEF,
  CS_EQKEY_HODGE_TIME_ALGO,
  CS_EQKEY_HODGE_TIME_COEF,
  CS_EQKEY_HODGE_REAC_ALGO,
  CS_EQKEY_HODGE_REAC_COEF,
  CS_EQKEY_ITSOL,
  CS_EQKEY_ITSOL_EPS,
  CS_EQKEY_ITSOL_MAX_ITER,
  CS_EQKEY_ITSOL_RESNORM,
  CS_EQKEY_PRECOND,
  CS_EQKEY_SLES_VERBOSITY,
  CS_EQKEY_SOLVER_FAMILY,
  CS_EQKEY_SPACE_SCHEME,
  CS_EQKEY_TIME_SCHEME,
  CS_EQKEY_TIME_THETA,
  CS_EQKEY_VERBOSITY,
  CS_EQKEY_N_KEYS

} cs_equation_key_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a pointer to a buffer of size at least the 2*n_cells for
 *         managing temporary usage of memory when dealing with equations
 *         Call specific structure allocation related to a numerical scheme
 *         according the scheme flag
 *         The size of the temporary buffer can be bigger according to the
 *         numerical settings
 *         Set also shared pointers from the main domain members
 *
 * \param[in]  connect       pointer to a cs_cdo_connect_t structure
 * \param[in]  quant         pointer to additional mesh quantities struct.
 * \param[in]  time_step     pointer to a time step structure
 * \param[in]  scheme_flag   flag to identify which kind of numerical scheme is
 *                           requested to solve the computational domain
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_allocate_common_structures(const cs_cdo_connect_t     *connect,
                                       const cs_cdo_quantities_t  *quant,
                                       const cs_time_step_t       *time_step,
                                       cs_flag_t                   scheme_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a pointer to a buffer of size at least the 2*n_cells for
 *         managing temporary usage of memory when dealing with equations
 *         Call specific structure allocation related to a numerical scheme
 *         according the scheme flag
 *         The size of the temporary buffer can be bigger according to the
 *         numerical settings
 *
 * \param[in]  scheme_flag   flag to identify which kind of numerical scheme is
 *                           requested to solve the computational domain
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_free_common_structures(cs_flag_t   scheme_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a pointer to a buffer of size at least the 2*n_cells
 *         The size of the temporary buffer can be bigger according to the
 *         numerical settings
 *
 * \return  a pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_equation_get_tmpbuf(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the allocation size of the temporary buffer
 *
 * \return  the size of the temporary buffer
 */
/*----------------------------------------------------------------------------*/

size_t
cs_equation_get_tmpbuf_size(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to store parameters related
 *         to an equation
 *
 * \param[in] eqname           name of the equation
 * \param[in] varname          name of the variable associated to this equation
 * \param[in] eqtype           type of equation (user, predefined...)
 * \param[in] vartype          type of variable (scalar, vector, tensor...)
 * \param[in] default_bc       type of boundary condition set by default
 *
 * \return  a pointer to the new allocated cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_equation_create(const char            *eqname,
                   const char            *varname,
                   cs_equation_type_t     eqtype,
                   cs_param_var_type_t    vartype,
                   cs_param_bc_type_t     default_bc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_equation_t structure
 *
 * \param[in, out] eq    pointer to a cs_equation_t structure
 *
 * \return  a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_equation_free(cs_equation_t  *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print a synthesis of the monitoring information in the performance
 *         file
 *
 * \param[in] eq    pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_print_monitoring(const cs_equation_t  *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of a cs_equation_t structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_summary(const cs_equation_t  *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create timer statistics structures to enable a "home-made" profiling
 *
 * \param[in, out]  eq       pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_timer_stats(cs_equation_t  *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign a set of pointer functions for managing the cs_equation_t
 *         structure during the computation
 *
 * \param[in]       connect  pointer to a cs_cdo_connect_t structure
 * \param[in, out]  eq       pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_last_setup(const cs_cdo_connect_t   *connect,
                       cs_equation_t            *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a parameter in a cs_equation_t structure attached to keyname
 *
 * \param[in, out]  eq       pointer to a cs_equation_t structure
 * \param[in]       key      key related to the member of eq to set
 * \param[in]       keyval   accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_param(cs_equation_t       *eq,
                      cs_equation_key_t    key,
                      const char          *keyval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Associate a material property or an advection field with an equation
 *         for a given term (diffusion, time, convection)
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       keyword   "time", "diffusion", "advection"
 * \param[in]       pointer   pointer to a given structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_link(cs_equation_t       *eq,
                 const char          *keyword,
                 void                *pointer);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the unknown related to this equation
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here a constant value is set to all the entities belonging to the
 *         given mesh location
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       ml_name   name of the associated mesh location (if NULL or
 *                            "" all cells are considered)
 * \param[in]       val       pointer to the value
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_ic_by_value(cs_equation_t    *eq,
                            const char       *ml_name,
                            cs_get_t          get);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the unknown related to this equation
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the value related to all the entities belonging to the
 *         given mesh location is such that the integral over these cells
 *         returns the requested quantity
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       ml_name   name of the associated mesh location (if NULL or
 *                            "" all cells are considered)
 * \param[in]       quantity  quantity to distribute over the mesh location
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_ic_by_qov(cs_equation_t    *eq,
                          const char       *ml_name,
                          double            quantity);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the unknown related to this equation
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set according to an analytical function
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       ml_name   name of the associated mesh location (if NULL or
 *                            "" all cells are considered)
 * \param[in]       analytic  pointer to an analytic function
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_ic_by_analytic(cs_equation_t        *eq,
                               const char           *ml_name,
                               cs_analytic_func_t   *analytic);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to set a boundary condition
 *         related to the givan equation structure
 *         ml_name corresponds to the name of a pre-existing cs_mesh_location_t
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       bc_type   type of boundary condition to add
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       get       pointer to a cs_get_t structure storing the value
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_bc_by_value(cs_equation_t              *eq,
                            const cs_param_bc_type_t    bc_type,
                            const char                 *ml_name,
                            const cs_get_t              get);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to set a boundary condition
 *         related to the givan equation structure
 *         ml_name corresponds to the name of a pre-existing cs_mesh_location_t
 *
 * \param[in, out] eq        pointer to a cs_equation_t structure
 * \param[in]      bc_type   type of boundary condition to add
 * \param[in]      ml_name   name of the related mesh location
 * \param[in]      analytic  pointer to an analytic function defining the value
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_bc_by_analytic(cs_equation_t              *eq,
                               const cs_param_bc_type_t    bc_type,
                               const char                 *ml_name,
                               cs_analytic_func_t         *analytic);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to store parameters related
 *         to a reaction term
 *
 * \param[in, out] eq         pointer to a cs_equation_t structure
 * \param[in]      property   pointer to a cs_property_t struct.
 * \param[in]      r_name     name of the reaction term (optional, i.e. NULL)
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_linear_reaction(cs_equation_t   *eq,
                                cs_property_t   *property,
                                const char      *r_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize by value a new structure to store parameters
 *         related to a source term defined by a user
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       st_name   name of the source term or NULL
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       val       pointer to the value
 *
 * \return a pointer to the new cs_source_term_t structure
 */
/*----------------------------------------------------------------------------*/

cs_source_term_t *
cs_equation_add_source_term_by_val(cs_equation_t   *eq,
                                   const char      *st_name,
                                   const char      *ml_name,
                                   const void      *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize by an analytical function a new structure
 *         related to a source term defined by a user
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       st_name   name of the source term or NULL
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       ana       pointer to an analytical function
 *
 * \return a pointer to the new cs_source_term_t structure
 */
/*----------------------------------------------------------------------------*/

cs_source_term_t *
cs_equation_add_source_term_by_analytic(cs_equation_t        *eq,
                                        const char           *st_name,
                                        const char           *ml_name,
                                        cs_analytic_func_t   *ana);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a field structure related to this cs_equation_t structure
 *         to an equation
 *
 * \param[in, out]  eq       pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_create_field(cs_equation_t      *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the values of a field according to the initial condition
 *         related to its equation
 *
 * \param[in]       mesh       pointer to the mesh structure
 * \param[in, out]  eq         pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_init_system(const cs_mesh_t            *mesh,
                        cs_equation_t              *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if one has to build the linear system
 *
 * \param[in]  eq        pointer to a cs_equation_t structure
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_needs_build(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system for this equation
 *
 * \param[in]       m           pointer to a cs_mesh_t structure
 * \param[in]       time_step   pointer to a time step structure
 * \param[in]       dt_cur      value of the current time step
 * \param[in, out]  eq          pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_build_system(const cs_mesh_t            *m,
                         const cs_time_step_t       *time_step,
                         double                      dt_cur,
                         cs_equation_t              *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the linear system for this equation
 *
 * \param[in, out]  eq          pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_solve(cs_equation_t   *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return true is the given equation is steady otherwise false
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_is_steady(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the values at each face of the mesh for the field unknowns
 *         related to this equation.
 *
 * \param[in]   eq        pointer to a cs_equation_t structure
 *
 * \return a pointer to the face values
 */
/*----------------------------------------------------------------------------*/

const cs_real_t *
cs_equation_get_face_values(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the values at each cell centers for the field unknowns
 *         related to this equation.
 *
 * \param[in]   eq        pointer to a cs_equation_t structure
 *
 * \return a pointer to the cell values
 */
/*----------------------------------------------------------------------------*/

const cs_real_t *
cs_equation_get_cell_values(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the name related to the given cs_equation_t structure
 *         to an equation
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a name or NULL if not found
 */
/*----------------------------------------------------------------------------*/

const char *
cs_equation_get_name(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the field structure associated to a cs_equation_t structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a cs_field_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_field_t *
cs_equation_get_field(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the flag associated to an equation
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a flag (cs_flag_t type)
 */
/*----------------------------------------------------------------------------*/

cs_flag_t
cs_equation_get_flag(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the cs_equation_param_t structure associated to a
 *         cs_equation_t structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a cs_equation_param_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

const cs_equation_param_t *
cs_equation_get_param(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a pointer to the cs_property_t structure associated to the
 *         diffusion term for this equation.
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a pointer to a cs_property_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_equation_get_diffusion_property(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a pointer to the cs_property_t structure associated to the
 *         unsteady term for this equation.
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a pointer to a cs_property_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_equation_get_time_property(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a pointer to the cs_property_t structure associated to the
 *         reaction term called r_name and related to this equation
 *
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a pointer to a cs_property_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_equation_get_reaction_property(const cs_equation_t    *eq,
                                  const char             *r_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the type of numerical scheme used for the discretization in
 *         space
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return  a cs_space_scheme_t variable
 */
/*----------------------------------------------------------------------------*/

cs_space_scheme_t
cs_equation_get_space_scheme(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the max. degree used in the polynomial basis for the space
 *         discretization
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return  the polynomial order
 */
/*----------------------------------------------------------------------------*/

int
cs_equation_get_space_poly_degree(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the type of variable solved by this equation
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return  the type of variable (sclar, vector...) associated to this equation
 */
/*----------------------------------------------------------------------------*/

cs_param_var_type_t
cs_equation_get_var_type(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the type of equation for the given equation structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return  the type of the given equation
 */
/*----------------------------------------------------------------------------*/

cs_equation_type_t
cs_equation_get_type(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the diffusive and convective flux accross a plane defined
 *         by a mesh location structure attached to the name ml_name.
 *
 * \param[in]      eq          pointer to a cs_equation_t structure
 * \param[in]      ml_name     name of the related mesh location
 * \param[in]      direction   vector indicating in which direction flux is > 0
 * \param[in, out] diff_flux   value of the diffusive part of the flux
 * \param[in, out] conv_flux   value of the convective part of the flux
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_flux_across_plane(const cs_equation_t   *eq,
                                      const char            *ml_name,
                                      const cs_real_3_t      direction,
                                      cs_real_t             *diff_flux,
                                      cs_real_t             *conv_flux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the diffusive flux across all cell faces.
 *         Primal or dual faces are considered according to the space scheme.
 *
 * \param[in]      eq          pointer to a cs_equation_t structure
 * \param[in]      location    indicate where the flux has to be computed
 * \param[in, out] diff_flux   value of the diffusive flux
  */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_diff_flux_cellwise(const cs_equation_t   *eq,
                                       cs_flag_t              location,
                                       cs_real_t             *diff_flux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the diffusive flux across all cell faces.
 *         Primal or dual faces are considered according to the space scheme.
 *
 * \param[in]      eq          pointer to a cs_equation_t structure
 * \param[in, out] diff_flux   value of the diffusive flux
  */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_diff_flux(const cs_equation_t   *eq,
                              cs_real_t             *diff_flux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]  eq      pointer to a cs_equation_t structure
 * \param[in]  ts      pointer to a cs_time_step_t struct.
 * \param[in]  dt      value of the current time step
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_extra_post(const cs_equation_t     *eq,
                       const cs_time_step_t    *ts,
                       double                   dt);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EQUATION_H__ */
