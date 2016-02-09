#ifndef __CS_EQUATION_H__
#define __CS_EQUATION_H__

/*============================================================================
 * Routines to handle cs_equation_t structure and its related structures
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

#include "cs_mesh.h"
#include "cs_time_step.h"

#include "cs_param.h"
#include "cs_field.h"
#include "cs_cdo_quantities.h"
#include "cs_equation_priv.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _cs_equation_t cs_equation_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

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
 * \brief  Summary of a cs_equation_t structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_summary(const cs_equation_t  *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign a set of pointer functions for managing the cs_equation_t
 *         structure during the computation
 *
 * \param[in, out]  eq       pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_last_setup(cs_equation_t  *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a parameter in a cs_equation_t structure attached to keyname
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       keyname   name of key related to the member of eq to set
 * \param[in]       val       accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_option(cs_equation_t       *eq,
                       const char          *keyname,
                       const void          *val);

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
 * \brief  Define the initial condition of the unknown related to this equation
 *         def_key is among "value", "analytic"
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       def_key   way of defining the value of the bc
 * \param[in]       val       pointer to the value
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_ic(cs_equation_t    *eq,
                   const char       *def_key,
                   void             *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to store parameters related
 *         to an equation
 *         bc_key among "dirichlet", "neumann" or "robin"
 *         def_key among "value", "analytic", "user"
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       bc_key    type of boundary condition to add
 * \param[in]       def_key   way of defining the value of the bc
 * \param[in]       val       pointer to the value
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_bc(cs_equation_t    *eq,
                   const char       *ml_name,
                   const char       *bc_key,
                   const char       *def_key,
                   const void       *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to store parameters related
 *         to a reaction term
 *
 * \param[in, out] eq         pointer to a cs_equation_t structure
 * \param[in]      r_name     name of the reaction term or NULL
 * \param[in]      type_name  type of reaction term to add
 * \param[in]      property   pointer to a cs_property_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_reaction(cs_equation_t   *eq,
                         const char      *r_name,
                         const char      *type_name,
                         cs_property_t   *property);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set advanced parameters related to a reaction term
 *         keyname among "lumping", "hodge_algo", "hodge_coef"...
 *         If r_name is NULL, all reaction terms of the given equation are set
 *         according to the couple (keyname, keyval)
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       r_name    name of the reaction term
 * \param[in]       keyname   name of the key
 * \param[in]       keyval    pointer to the value to set to the key
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_reaction_option(cs_equation_t    *eq,
                                const char       *r_name,
                                const char       *keyname,
                                const char       *keyval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to store parameters related
 *         to a source term
 *         def_key among "value", "analytic", "user"...
 *
 * \param[in, out]  eq             pointer to a cs_equation_t structure
 * \param[in]       ml_id          id related to a mesh location
 * \param[in]       array_support  indicate where the values are defined
 * \param[in]       array_values   pointer to the array values
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_gravity_source_term(cs_equation_t   *eq,
                                    int              ml_id,
                                    cs_flag_t        array_support,
                                    cs_real_t       *array_values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize by value a new structure to store parameters
 *         related to a source term defined by a user
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       st_name   name of the source term or NULL
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       val       pointer to the value
 */
/*----------------------------------------------------------------------------*/

void
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
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_source_term_by_analytic(cs_equation_t        *eq,
                                        const char           *st_name,
                                        const char           *ml_name,
                                        cs_analytic_func_t   *ana);


/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set members defined by default in a source term structure
 *         keyname among "quadrature", "post"...
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       st_name   name of the source term
 * \param[in]       keyname   name of the key
 * \param[in]       keyval    pointer to the value to set to the key
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_source_term_option(cs_equation_t    *eq,
                                   const char       *st_name,
                                   const char       *keyname,
                                   const char       *keyval);

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
 * \param[in]       connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]       cdoq       pointer to a cs_cdo_quantities_t struct.
 * \param[in]       time_step  pointer to a time step structure
 * \param[in, out]  eq         pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_init_system(const cs_mesh_t            *mesh,
                        const cs_cdo_connect_t     *connect,
                        const cs_cdo_quantities_t  *cdoq,
                        const cs_time_step_t       *time_step,
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
 * \param[in]       do_logcvg   output information on convergence or not
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_solve(cs_equation_t   *eq,
                  bool             do_logcvg);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]  time_step  pointer to a time step structure
 * \param[in]  eq         pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_extra_op(const cs_time_step_t       *time_step,
                     const cs_equation_t        *eq);

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
 * \brief  Compute the values of the associated field at each face of the mesh
 *         If the pointer storing the values is NULL, it is alloacted inside the
 *         function
 *
 * \param[in]       eq       pointer to a cs_equation_t structure
 *
 * \return a pointer to the values
 */
/*----------------------------------------------------------------------------*/

const cs_real_t *
cs_equation_get_face_values(const cs_equation_t    *eq);

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

END_C_DECLS

#endif /* __CS_EQUATION_H__ */
