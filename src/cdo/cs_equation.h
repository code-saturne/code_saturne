#ifndef __CS_EQUATION_H__
#define __CS_EQUATION_H__

/*============================================================================
 * Routines to handle cs_equation_t structure and its related structures
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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
 * \param[in] status           user or predefined equations
 * \param[in] type             type of equation (scalar, vector, tensor...)
 * \param[in] is_steady        add an unsteady term or not
 * \param[in] do_convection    add a convection term
 * \param[in] do_diffusion     add a diffusion term
 * \param[in] default_bc       type of boundary condition set by default
 *
 * \return  a pointer to the new allocated cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_equation_create(const char            *name,
                   const char            *varname,
                   cs_equation_status_t   status,
                   cs_equation_type_t     type,
                   bool                   is_steady,
                   bool                   do_convection,
                   bool                   do_diffusion,
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
 * \brief  Resume a cs_equation_t structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_resume(const cs_equation_t  *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a builder structure
 *
 * \param[in]       mesh     pointer to a cs_mesh_t structure
 * \param[in, out]  eq       pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_create_builder(const cs_mesh_t   *mesh,
                           cs_equation_t     *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign a set of pointer functions for managing the cs_equation_t
 *         structure during the computation
 *
 * \param[in, out]  eq       pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_last_init(cs_equation_t  *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to store parameters related
 *         to an equation
 *
 * \param[in, out]  eq       pointer to a cs_equation_t structure
 * \param[in]       key      key related to the member of eq to set
 * \param[in]       val      accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set(cs_equation_t       *eq,
                const char          *key,
                const void          *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to store parameters related
 *         to an equation
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       pty_key   "time", "diffusion"...
 * \param[in]       pty_name  name of the material property to associate
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_pty(cs_equation_t       *eq,
                    const char          *pty_key,
                    const char          *pty_name);

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
 *         to a source term
 *         st_key among "implicit", "explicit", "imex"...
 *         def_key among "value", "analytic", "user"...
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       st_name   name of the source term or NULL
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       st_key    type of boundary condition to add
 * \param[in]       def_key   way of defining the value of the bc
 * \param[in]       val       pointer to the value
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_source_term(cs_equation_t    *eq,
                            const char       *st_name,
                            const char       *ml_name,
                            const char       *st_key,
                            const char       *def_key,
                            const void       *val);

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
cs_equation_source_term_set(cs_equation_t    *eq,
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
cs_equation_create_field(cs_equation_t    *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if one has to solve this equation now
 *
 * \param[in]  eq          pointer to a cs_equation_t structure
 * \param[in]  time_iter   id of the time iteration
 * \param[in]  tcur        current physical time
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_needs_solve(const cs_equation_t   *eq,
                        int                    time_iter,
                        double                 tcur);

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
 * \param[in]       m        pointer to a cs_mesh_t structure
 * \param[in]       connect  pointer to a cs_cdo_connect_t structure
 * \param[in]       cdoq     pointer to a cs_cdo_quantities_t structure
 * \param[in]       tcur     current physical time of the simulation
 * \param[in, out]  eq       pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_build_system(const cs_mesh_t            *m,
                         const cs_cdo_connect_t     *connect,
                         const cs_cdo_quantities_t  *cdoq,
                         double                      tcur,
                         cs_equation_t              *eq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the linear system for this equation
 *
 * \param[in]       connect  pointer to a cs_cdo_connect_t structure
 * \param[in]       cdoq     pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  eq       pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_solve(const cs_cdo_connect_t     *connect,
                  const cs_cdo_quantities_t  *cdoq,
                  cs_equation_t              *eq);

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

END_C_DECLS

#endif /* __CS_EQUATION_H__ */
