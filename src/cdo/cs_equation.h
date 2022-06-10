#ifndef __CS_EQUATION_H__
#define __CS_EQUATION_H__

/*============================================================================
 * Functions to handle the cs_equation_t structure and its related structures
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"
#include "cs_cdo_toolbox.h"
#include "cs_equation_param.h"
#include "cs_equation_builder.h"
#include "cs_field.h"
#include "cs_mesh.h"
#include "cs_restart.h"
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

/*!
 * \struct cs_equation_core_t
 * \brief Main structures composiing an equation structure
 *
 * Intermediate structure useful to manipulate an array of (sub-)structures.
 * Especially, the scheme context relies on the space discretization and it is
 * not easy to manipulate void ** object. This is a way around to this
 * operation.
 *
 * These three structures allow one to use nearly all operations related to
 * an equation without having to build an equation structure. This is useful
 * when handling extra-diagonal block in systems of equations.
 *
 * \var param
 *      Set of parameters to specifiy the settings
 *
 * \var builder
 *      Part of the quantities useful to build/manipulate an equation. All
 *      quantities that are shared among all discretizations are in this
 *      structure.
 *
 * \var scheme_context
 *      Part of the quantities useful to build/manipulate an equation. All
 *      quantities that are specific to the discrization are in this structure.
 */

typedef struct {

  cs_equation_param_t     *param;
  cs_equation_builder_t   *builder;
  void                    *scheme_context;

} cs_equation_core_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the number of equations
 *
 * \return the current number of cs_equation_t structure allocated
 */
/*----------------------------------------------------------------------------*/

int
cs_equation_get_n_equations(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Find the cs_equation_t structure with name eqname
 *         Return NULL if not find
 *
 * \param[in]  eqname    name of the equation to find
 *
 * \return a pointer to a cs_equation_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_equation_by_name(const char    *eqname);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the pointer to a cs_equation_t structure thanks to the field
 *         name of the variable field associated to a cs_equation_t structure
 *
 * \param[in]  field_name       name of the field
 *
 * \return a pointer to a cs_equation_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_equation_by_field_name(const char    *field_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the asociated field to a \ref cs_equation_t structure
 *         has name equal to fld_name
 *
 * \param[in]  eq          pointer to a \ref cs_equation_t structure to test
 * \param[in]  fld_name    name of the field
 *
 * \return true if the \ref cs_equation_t structure has an associated field
 *         named fld_name, otherwise false
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_has_field_name(const cs_equation_t  *eq,
                           const char           *fld_name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the cs_equation_param_t structure associated to a
 *         cs_equation_t structure based on the equation name
 *
 * If no equation matches the given name but a field does, equation
 * parameter structure associated to the field will be returned instead.
 * This allows using this function with non-CDO (legacy) fields.
 *
 * \param[in]  eqname       name of the equation
 *
 * \return a cs_equation_param_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_equation_param_by_name(const char    *eqname);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the cs_equation_param_t structure related to a
 *         cs_equation_t structure thanks to the field name of the variable
 *         field associated to a cs_equation_t structure
 *
 * \param[in]  field_name       name of the field
 *
 * \return a cs_equation_param_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_equation_param_by_field_name(const char    *field_name);

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

cs_equation_param_t *
cs_equation_get_param(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Find the \ref cs_equation_t structure with id eq_id
 *         Return NULL if not find
 *
 * \param[in]  eq_id    id of the equation to find
 *
 * \return a pointer to a \ref cs_equation_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_equation_by_id(int   eq_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the name related to the given cs_equation_t structure
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
 * \brief  Return the id number related to the given cs_equation_t structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return an id (0 ... n-1) or -1 if not found
 */
/*----------------------------------------------------------------------------*/

int
cs_equation_get_id(const cs_equation_t    *eq);

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
 * \brief  Return the id related to the variable field structure associated to
 *         the cs_equation_t structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return an integer (-1 if the field is not defined)
 */
/*----------------------------------------------------------------------------*/

int
cs_equation_get_field_id(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the range set structure associated to a cs_equation_t
 *         structure. One assumes that there is only one block (it could be a
 *         split block) otherwise this means that one handles systems of
 *         equations.
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a cs_range_set_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

const cs_range_set_t *
cs_equation_get_range_set(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the global number of degrees of freedom associated to this
 *         cs_equation_t structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 * \param[in]  cdoq     pointer to a cs_cdo_quantities_t structure
 *
 * \return a global number of degrees of freedom (DoFs)
 */
/*----------------------------------------------------------------------------*/

cs_gnum_t
cs_equation_get_global_n_dofs(const cs_equation_t         *eq,
                              const cs_cdo_quantities_t   *cdoq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the field structure for the (normal) boundary flux associated
 *         to a cs_equation_t structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a cs_field_t structure or NULL
 */
/*----------------------------------------------------------------------------*/

cs_field_t *
cs_equation_get_boundary_flux(const cs_equation_t    *eq);

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
 * \brief  Redefine the flag associated to an equation
 *
 * \param[in, out]  eq       pointer to a cs_equation_t structure
 * \param[in]       flag     new flag to set
*/
/*----------------------------------------------------------------------------*/

void
cs_equation_set_flag(cs_equation_t    *eq,
                     cs_flag_t         flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a hook function to enable an advanced control during the
 *         cellwise system building.
 *         Only for an advanced usage. The context may be set to NULL if there
 *         is no need to get additional information.
 *
 * \param[in, out] eq        pointer to the cs_equation_t stucture to update
 * \param[in]      context   pointer to a structure for additional information
 * \param[in]      func      pointer to the build function
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_build_hook(cs_equation_t               *eq,
                           void                        *context,
                           cs_equation_build_hook_t    *func);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the cs_equation_builder_t structure associated to a
 *         cs_equation_t structure. Only for an advanced usage.
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a cs_equation_builder_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_equation_builder_t *
cs_equation_get_builder(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a pointer to a structure useful to handle low-level
 *         operations for the given equation
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a pointer to a structure to cast on-the-fly or NULL if not found
 */
/*----------------------------------------------------------------------------*/

void *
cs_equation_get_scheme_context(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a pointer to a structure useful to handle low-level
 *         operations for the given equation
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return  structure storing the main structure associated to an equation
 */
/*----------------------------------------------------------------------------*/

cs_equation_core_t
cs_equation_get_core(const cs_equation_t    *eq);

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
 * \brief  Return true is the given equation is steady otherwise false
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_uses_new_mechanism(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new equation structure and set a first set of parameters
 *
 * \param[in] eqname        name of the equation
 * \param[in] varname       name of the variable associated to this equation
 * \param[in] eqtype        type of equation (user, predefined...)
 * \param[in] dim           dimension of the unknow attached to this equation
 * \param[in] default_bc    type of boundary condition set by default
 *
 * \return  a pointer to the new allocated cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_equation_add(const char            *eqname,
                const char            *varname,
                cs_equation_type_t     eqtype,
                int                    dim,
                cs_param_bc_type_t     default_bc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new user equation structure and set a first set of parameters
 *
 * \param[in] eqname        name of the equation
 * \param[in] varname       name of the variable associated to this equation
 * \param[in] dim           dimension of the unknow attached to this equation
 * \param[in] default_bc    type of boundary condition set by default
 *
 * \return  a pointer to the new allocated cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_equation_add_user(const char            *eqname,
                     const char            *varname,
                     int                    dim,
                     cs_param_bc_type_t     default_bc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new user transport equation and set a first set of parameters
 *         If time_pty is NULL, then no unsteady term is added.
 *         If adv is NULL, then no advection term is added.
 *         If diff_pty is NULL, then no diffusion term is added.
 *
 * \param[in] eqname       name of the equation
 * \param[in] varname      name of the variable associated to this equation
 * \param[in] dim          dimension of the unknow attached to this equation
 * \param[in] default_bc   type of boundary condition set by default
 * \param[in] time_pty     property related to the unsteady term
 * \param[in] adv          advection field
 * \param[in] diff_pty     property related to the diffusion term
 *
 * \return  a pointer to the new allocated cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_equation_add_user_tracer(const char            *eqname,
                            const char            *varname,
                            int                    dim,
                            cs_param_bc_type_t     default_bc,
                            cs_property_t         *time_pty,
                            cs_adv_field_t        *adv,
                            cs_property_t         *diff_pty);
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy all cs_equation_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_destroy_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if a steady-state computation is requested according to the
 *         setting
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_needs_steady_state_solve(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print a synthesis of the monitoring information in the performance
 *         file
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_log_monitoring(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the count of equations of each macro type
 *
 * \param[out]  n_equations          total number of equations
 * \param[out]  n_predef_equations   number of predefined equations
 * \param[out]  n_user_equations     number of user equations
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_get_count(int      *n_equations,
                      int      *n_predef_equations,
                      int      *n_user_equations);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summarize all cs_equation_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_log_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a parameter attached to a keyname for the default settigns
 *
 * \param[in] key      key related to the member of eq to set
 * \param[in] keyval   accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_default_param(cs_equation_key_t      key,
                              const char            *keyval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup the linear algebra requirements
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_sles(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers to the main structures. Associate these
 *         structures among the activated class of discretization schemes
 *
 * \param[in]  connect          pointer to a cs_cdo_connect_t structure
 * \param[in]  quant            pointer to additional mesh quantities struct.
 * \param[in]  time_step        pointer to a time step structure
 * \param[in]  eb_scheme_flag   metadata for Eb schemes
 * \param[in]  fb_scheme_flag   metadata for Fb schemes
 * \param[in]  vb_scheme_flag   metadata for Vb schemes
 * \param[in]  vcb_scheme_flag  metadata for V+C schemes
 * \param[in]  hho_scheme_flag  metadata for HHO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_init_sharing(const cs_cdo_connect_t      *connect,
                         const cs_cdo_quantities_t   *quant,
                         const cs_time_step_t        *time_step,
                         cs_flag_t                    eb_scheme_flag,
                         cs_flag_t                    fb_scheme_flag,
                         cs_flag_t                    vb_scheme_flag,
                         cs_flag_t                    vcb_scheme_flag,
                         cs_flag_t                    hho_scheme_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free shared local structures among the discretization schemes
 *
 * \param[in]  vb_scheme_flag   metadata for Vb schemes
 * \param[in]  vcb_scheme_flag  metadata for V+C schemes
 * \param[in]  eb_scheme_flag   metadata for Eb schemes
 * \param[in]  fb_scheme_flag   metadata for Fb schemes
 * \param[in]  hho_scheme_flag  metadata for HHO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_finalize_sharing(cs_flag_t    vb_scheme_flag,
                             cs_flag_t    vcb_scheme_flag,
                             cs_flag_t    eb_scheme_flag,
                             cs_flag_t    fb_scheme_flag,
                             cs_flag_t    hho_scheme_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign a set of pointer functions for managing the cs_equation_t
 *         structure during the computation
 *         After this call, parameters related to an equation are set once for
 *         all
 *
 * \return true if all equations are steady-state otherwise false
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_set_functions(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a field structure related to the predefined equation given as
 *         parameter. This includes an equation associated to all modules and
 *         also wall distance or mesh deformation for instance
 *
 *         When an automatic behavior is asked then one checks the flag
 *         CS_EQUATION_UNSTEADY to decide. One can force the behavior when
 *         handling predefined equations since more complex situations can
 *         occur such as a steady computation with non-linearities (in which
 *         case one wants a field with a previous state)
 *
 * \param[in]       n_previous     number of previous states to keep
 *                                 -1 means automatic
 * \param[in, out]  eq             pointer to an equation structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_predefined_create_field(int               n_previous,
                                    cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a field structure related to all user-defined equations
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_user_create_fields(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and define the builder structure
 *
 * \param[in]  mesh      pointer to a cs_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_define_builders(const cs_mesh_t     *mesh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and define the context structure associated to each equation
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_define_context_structures(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a pointer to a core structure. If the input core structure is
 *         not allocated, then one allocates the structure.
 *
 * \param[in]       eq       pointer to a cs_equation_t structure
 * \param[in, out]  p_core   double pointer to a core structure to build
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_define_core(const cs_equation_t    *eq,
                        cs_equation_core_t    **p_core);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the initialize condition to all variable fields associated to
 *         each cs_equation_t structure.
 *
 * \param[in]  mesh      pointer to a cs_mesh_t structure
 * \param[in]  ts        pointer to a cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_init_field_values(const cs_mesh_t             *mesh,
                              const cs_time_step_t        *ts);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and then solve the linear system for this equation when the
 *         goal is to find the steady state
 *
 * \param[in]       mesh        pointer to a cs_mesh_t structure
 * \param[in, out]  eq          pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_solve_steady_state(const cs_mesh_t            *mesh,
                               cs_equation_t              *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and then solve the linear system for an equation with an
 *         unsteady term
 *
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in, out] eq         pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_solve(bool                        cur2prev,
                  const cs_mesh_t            *mesh,
                  cs_equation_t              *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and then solve the linear system for a steady-state equation.
 *         This is wrapper for the FORTRAN interface (limitation of the
 *         parameters to simple types).
 *
 * \param[in] eqname     name of the equation to solve
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_solve_steady_state_wrapper(const char                 *eqname);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and then solve the linear system for an equation with an
 *         unsteady term. This is wrapper for the FORTRAN interface (limitation
 *         of the parameters to simple types)
 *
 * \param[in] cur2prev   true="current to previous" operation is performed
 * \param[in] eqname     name of the equation to solve
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_solve_wrapper(bool                        cur2prev,
                          const char                 *eqname);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system for this equation (deprecated)
 *
 * \param[in]       mesh        pointer to a cs_mesh_t structure
 * \param[in, out]  eq          pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_build_system(const cs_mesh_t            *mesh,
                         cs_equation_t              *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the linear system for this equation (deprecated)
 *
 * \param[in, out]  eq          pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_solve_deprecated(cs_equation_t   *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a pointer to the cs_property_t structure associated to the
 *         diffusion term for this equation (NULL if not activated).
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_equation_get_diffusion_property(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a pointer to the cs_property_t structure associated to the
 *         unsteady term for this equation (NULL if not activated).
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_equation_get_time_property(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a pointer to the cs_property_t structure associated to the
 *         reaction term with id equal to reaction_id and related to this
 *         equation
 *
 *
 * \param[in]  eq            pointer to a cs_equation_t structure
 * \param[in]  reaction_id   id related to this reaction term
 *
 * \return a pointer to a cs_property_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_equation_get_reaction_property(const cs_equation_t    *eq,
                                  const int               reaction_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the type of numerical scheme used for the discretization in
 *         time
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return  a cs_param_time_scheme_t variable
 */
/*----------------------------------------------------------------------------*/

cs_param_time_scheme_t
cs_equation_get_time_scheme(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the value of the theta parameter in theta time scheme
 *         discretization
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return  the value of the theta coefficient. -1 if the time scheme is not
 *          a theta time scheme
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_equation_get_theta_time_val(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the type of numerical scheme used for the discretization in
 *         space
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return  a cs_param_space_scheme_t variable
 */
/*----------------------------------------------------------------------------*/

cs_param_space_scheme_t
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
 * \brief  Return the dimension of the variable solved by this equation
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return  an integer corresponding to the dimension of the variable
 */
/*----------------------------------------------------------------------------*/

int
cs_equation_get_var_dim(const cs_equation_t    *eq);

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
 * \brief  Apply the current to previous to all fields (and potentially arrays)
 *         related to an equation. This function fas to be called when a solve
 *         step is called with the parameter: cur2prev = false
 *
 * \param[in]   eq       pointer to a \ref cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_current_to_previous(const cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  For a given equation, retrieve the related cellwise builder
 *         structures: cs_cell_builder_t and cs_cell_system_t structures
 *
 * \param[in]   eq       pointer to a \ref cs_equation_t structure
 * \param[out]  cb       pointer to a pointer on a cs_cell_sys_t structure
 * \param[out]  csys     pointer to a pointer on a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_get_cellwise_builders(const cs_equation_t    *eq,
                                  cs_cell_sys_t         **csys,
                                  cs_cell_builder_t     **cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  For a given equation, retrieve an array of values related to each
 *         cell of the mesh for the unknowns
 *
 * \param[in]   eq        pointer to a \ref cs_equation_t structure
 * \param[in]   previous  retrieve the previous state (true/false)
 *
 * \return a pointer to an array of cell values
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_equation_get_cell_values(const cs_equation_t    *eq,
                            bool                    previous);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  For a given equation, retrieve an array of values related to each
 *         face of the mesh for the unknowns
 *
 * \param[in]   eq        pointer to a \ref cs_equation_t structure
 * \param[in]   previous  retrieve the previous state (true/false)
 *
 * \return a pointer to an array of face values
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_equation_get_face_values(const cs_equation_t    *eq,
                            bool                    previous);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  For a given equation, retrieve an array of values related to each
 *         edge of the mesh for the unknowns
 *
 * \param[in]   eq        pointer to a \ref cs_equation_t structure
 * \param[in]   previous  retrieve the previous state (true/false)
 *
 * \return a pointer to an array of edge values
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_equation_get_edge_values(const cs_equation_t    *eq,
                            bool                    previous);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  For a given equation, retrieve an array of values related to each
 *         vertex of the mesh for the unknowns
 *
 * \param[in]   eq        pointer to a \ref cs_equation_t structure
 * \param[in]   previous  retrieve the previous state (true/false)
 *
 * \return a pointer to an array of vertex values
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_equation_get_vertex_values(const cs_equation_t    *eq,
                              bool                    previous);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over the domain of the variable field
 *         associated to the given equation.
 *
 * \param[in]      connect    pointer to a \ref cs_cdo_connect_t structure
 * \param[in]      cdoq       pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]      eq         pointer to a \ref cs_equation_t structure
 * \param[in, out] integral   result of the computation
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_integrate_variable(const cs_cdo_connect_t     *connect,
                               const cs_cdo_quantities_t  *cdoq,
                               const cs_equation_t        *eq,
                               cs_real_t                  *result);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the diffusive flux across all boundary faces
 *         According to the space discretization scheme, the size of the
 *         resulting array differs.
 *         For Vb and VCb schemes, this array relies on the bf2v adjacency.
 *
 * \param[in]      t_eval     time at which one performs the property evaluation
 * \param[in]      eq         pointer to a cs_equation_t structure
 * \param[in, out] diff_flux  value of the diffusive part of the flux
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_boundary_diff_flux(cs_real_t              t_eval,
                                       const cs_equation_t   *eq,
                                       cs_real_t             *diff_flux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the diffusive and convective flux across a plane defined
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
 * \brief Cellwise computation of the diffusive flux across the requested
 *        location. If the location is not the "natural" one (which depends on
 *        the space discretization scheme) then the diffusive flux is only an
 *        approximation.
 *
 * \param[in]      eq          pointer to a cs_equation_t structure
 * \param[in]      location    indicate where the flux has to be computed
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in, out] diff_flux   value of the diffusive flux (must be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_diffusive_flux(const cs_equation_t   *eq,
                                   cs_flag_t              location,
                                   cs_real_t              t_eval,
                                   cs_real_t             *diff_flux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the discrete gradient at vertices
 *
 * \param[in]      eq          pointer to a cs_equation_t structure
 * \param[in, out] v_gradient  gradient at vertices
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_vtx_field_gradient(const cs_equation_t   *eq,
                                       cs_real_t             *v_gradient);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute and post-process Peclet number if requested
 *
 * \param[in]      eq       pointer to a cs_equation_t structure
 * \param[in]      ts       pointer to a cs_time_step_t struct.
 * \param[in, out] peclet   pointer to an array storing the resulting Peclet
 *                          number in each cell
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_peclet(const cs_equation_t        *eq,
                           const cs_time_step_t       *ts,
                           cs_real_t                   peclet[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Write into the restart file additionnal arrays (not defined as
 *         fields) but useful for the checkpoint/restart process
 *
 * \param[in, out]  restart    pointer to a \ref cs_restart_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_read_extra_restart(cs_restart_t   *restart);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Write into the restart file additionnal arrays (not defined as
 *         fields) but useful for the checkpoint/restart process
 *
 * \param[in, out]  restart    pointer to a \ref cs_restart_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_write_extra_restart(cs_restart_t   *restart);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to all equations
 *
 * \param[in]  mesh      pointer to a cs_mesh_t structure
 * \param[in]  connect   pointer to a cs_cdo_connect_t structure
 * \param[in]  cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in]  ts        pointer to a cs_time_step_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_post_balance(const cs_mesh_t            *mesh,
                         const cs_cdo_connect_t     *connect,
                         const cs_cdo_quantities_t  *cdoq,
                         const cs_time_step_t       *ts);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the cellwise stiffness matrix associated to the property
 *         given as a parameter and apply it to the pot array to define
 *         the resulting array associated to entities defined at loc_res
 *
 * \param[in]      eq        pointer to a \ref cs_equation_t structure
 * \param[in]      property  pointer to the property to consider
 * \param[in]      pot       array to multiply with the stiffness matrix
 * \param[in]      loc_res   location of entities in the resulting array
 * \param[in, out] res       resulting array
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_apply_stiffness(cs_equation_t          *eq,
                            const cs_property_t    *property,
                            const cs_real_t        *pot,
                            cs_flag_t               loc_res,
                            cs_real_t              *res);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to equations according to the
 *         type of numerical scheme (for the space discretization)
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_extra_post(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EQUATION_H__ */
