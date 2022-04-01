/*============================================================================
 * User definition of boundary conditions.
 *============================================================================*/

/* VERS */

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

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_boundary_conditions.c
 *
 * \brief User functions for input of calculation parameters.
 *
 * See \ref parameters for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the velocity inlet profile as an analytic function.
 *         elt_ids is optional. If not NULL, it enables to access to the coords
 *         array with an indirection. The same indirection can be applied to
 *         the val array if dense_output is set to false.
 *
 * \param[in]      time          when ?
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids (in coords and retval)
 * \param[in]      coords        where ? Coordinates array
 * \param[in]      dense_output  perform an indirection in res or not
 * \param[in]      input         NULL or pointer to a structure cast on-the-fly
 * \param[in, out] val           resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

/*! [inlet_vel_analytic_func] */
static void
_vel_profile(cs_real_t           time,
             cs_lnum_t           n_elts,
             const cs_lnum_t    *elt_ids,
             const cs_real_t    *coords,
             bool                dense_output,
             void               *input,
             cs_real_t          *val)
{
  CS_UNUSED(time);
  CS_UNUSED(input);

  const cs_real_3_t *elt_coords = (const cs_real_3_t *)coords;

  cs_real_3_t  *v = (cs_real_3_t *)val;

  for (cs_lnum_t i = 0; i < n_elts; i++) {

    const cs_lnum_t  elt_id = (elt_ids == NULL) ? i : elt_ids[i];
    const cs_lnum_t  j = dense_output ? i : elt_id;

    cs_real_t  y = elt_coords[elt_id][1] - 0.5;
    cs_real_t  z = elt_coords[elt_id][2] - 0.5;
    cs_real_t r = sqrt(y*y + z*z);

    v[j][0] = 1.5 * (1 - r);
    v[j][1] = 0;
    v[j][2] = 0;

  }
}
/*! [inlet_vel_analytic_func] */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the scalar inlet profile as an analytic function.
 *         elt_ids is optional. If not NULL, it enables to access to the coords
 *         array with an indirection. The same indirection can be applied to
 *         the val array if dense_output is set to false.
 *
 * \param[in]      time          when ?
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids (in coords and retval)
 * \param[in]      coords        where ? Coordinates array
 * \param[in]      dense_output  perform an indirection in res or not
 * \param[in]      input         NULL or pointer to a structure cast on-the-fly
 * \param[in, out] val           resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

/*! [inlet_scalar_analytic_func] */
static void
_scalar_inlet_profile(cs_real_t           time,
                      cs_lnum_t           n_elts,
                      const cs_lnum_t    *elt_ids,
                      const cs_real_t    *coords,
                      bool                dense_output,
                      void               *input,
                      cs_real_t          *val)
{
  CS_UNUSED(time);

  const cs_real_3_t *elt_coords = (const cs_real_3_t *)coords;

  const cs_field_t *f = input;  /* field pointer passed as input
                                   upon assignment */

  if (strcmp(f->name, "scalar1") == 0) {

    for (cs_lnum_t i = 0; i < n_elts; i++) {

      const cs_lnum_t  elt_id = (elt_ids == NULL) ? i : elt_ids[i];
      const cs_lnum_t  j = dense_output ? i : elt_id;

      cs_real_t  y = elt_coords[elt_id][1] - 0.5;
      cs_real_t  z = elt_coords[elt_id][2] - 0.5;
      cs_real_t r = sqrt(y*y + z*z);

      val[j] = 4.0 * (1 - r);

    }
  }
}
/*! [inlet_scalar_analytic_func] */

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_dof_func_t function to compute a profile at boundary faces
 *        using a MEG generated function for exchange coefficients.
 *
 * For the calling function, elt_ids is optional. If not NULL, array(s) should
 * be accessed with an indirection. The same indirection can be applied to fill
 * retval if dense_output is set to false.
 * In the current case, retval is allocated to mesh->n_b_faces * stride.
 *
 * The retval values can be decomposed into the alpha, u0, and g values in:
 *
 * K du/dn + alpha*(u - u0) = g
 *
 * For multidimentsional cases, we assume scalar alpha, variable dimension u0,
 * and dimension^2 g (g = 0 here but storage required), so we have a stride
 * of 1 + dim + dim*dim, with values alpha, u0, and g in order.
 *
 * This example assumes a scalar variable, with a resulting stride of 3.
 * Here we set a constant exchange coefficient: alpha = -10.5,
 *             a variable external value:       u0    = 25. + 0.1*x
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         NULL or pointer to a structure cast on-the-fly
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

/*! [wall_exchange_dof_func] */
static void
_scalar_exchange_profile(cs_lnum_t         n_elts,
                         const cs_lnum_t  *elt_ids,
                         bool              dense_output,
                         void             *input,
                         cs_real_t        *retval)
{
  CS_UNUSED(input);

  const cs_lnum_t dim = 1;
  const cs_lnum_t stride = 1 + dim + dim*dim;

  /* Face center coordinates */

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_real_3_t *f_cog = (const cs_real_3_t *)mq->b_face_cog;

  /* Exchange coefficient first, Dirichlet values second. */

  for (cs_lnum_t i = 0; i < n_elts; i++) {
    const cs_lnum_t  face_id = (elt_ids == NULL) ? i : elt_ids[i];
    const cs_lnum_t  j = dense_output ? i : face_id;
    retval[j*stride]   = - 10.5;
    retval[j*stride+1] = 25. + 0.1*f_cog[face_id][0];
    retval[j*stride*3+2] = 0;
  }
}
/*! [wall_exchange_dof_func] */

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_dof_func_t function to compute a wall flux normalized by the
 *        associated zone's surface.
 *
 * For the calling function, elt_ids is optional. If not NULL, array(s) should
 * be accessed with an indirection. The same indirection can be applied to fill
 * retval if dense_output is set to false.
 * In the current case, retval is allocated to mesh->n_b_faces.
 *
 * In this example, the input variable is assumed to be a pointer
 * to the associated zone.
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         pointer to associated zone
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

/*! [wall_top_flux_dof_func] */
static void
_w_flux_top(cs_lnum_t         n_elts,
            const cs_lnum_t  *elt_ids,
            bool              dense_output,
            void             *input,
            cs_real_t        *retval)
{
  CS_NO_WARN_IF_UNUSED(input);

  const cs_zone_t *z = cs_boundary_zone_by_name("wall_top");

  /* Get the fluid measure (i.e. surface) of the zone */
  cs_real_t flux = 1. / z->f_measure;

  /* Exchange coefficient first, Dirichlet values second. */

  for (cs_lnum_t i = 0; i < n_elts; i++) {
    const cs_lnum_t  face_id = (elt_ids == NULL) ? i : elt_ids[i];
    const cs_lnum_t  j = dense_output ? i : face_id;
    retval[j] = flux;
  }
}
/*! [wall_top_flux_dof_func] */

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_dof_func_t function to compute a wall flux normalized by the
 *        associated zone's surface.
 *
 * For the calling function, elt_ids is optional. If not NULL, array(s) should
 * be accessed with an indirection. The same indirection can be applied to fill
 * retval if dense_output is set to false.
 * In the current case, retval is allocated to mesh->n_b_faces.
 *
 * In this example, the input variable is assumed to be a pointer
 * to the associated zone.
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         pointer to associated zone
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

/*! [wall_side_flux_dof_func] */
static void
_w_flux_side(cs_lnum_t         n_elts,
             const cs_lnum_t  *elt_ids,
             bool              dense_output,
             void             *input,
             cs_real_t        *retval)
{
  const cs_zone_t *z = (const cs_zone_t *)input;

  /* Get the fluid measure (i.e. surface) of the zone */
  cs_real_t flux = 1. / z->f_measure;

  /* Exchange coefficient first, Dirichlet values second. */

  for (cs_lnum_t i = 0; i < n_elts; i++) {
    const cs_lnum_t  face_id = (elt_ids == NULL) ? i : elt_ids[i];
    const cs_lnum_t  j = dense_output ? i : face_id;
    retval[j] = flux;
  }
}
/*! [wall_side_flux_dof_func] */

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set boundary conditions to be applied.
 *
 * This function is called just before \ref cs_user_finalize_setup, and
 * boundary conditions can be defined in either of those functions,
 * depending on whichever is considered more readable or practical for a
 * given use.
 *
 * \param[in, out]  domain  pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_boundary_conditions_setup(cs_domain_t  *domain)
{
  CS_UNUSED(domain);

  /* Example: define constant inlet velocity with vector (1,0,0) */
  {
    /*! [inlet_vel_value] */
    cs_equation_param_t  *eqp = cs_equation_param_by_name("velocity");

    cs_real_t inlet_velocity[] = {1.5, 0, 0};

    cs_equation_add_bc_by_value(eqp,
                                CS_PARAM_BC_DIRICHLET,
                                "inlet",           // zone name
                                inlet_velocity);

    /*! [inlet_vel_value] */
  }

  /* Example: define inlet velocity using user-defined callback function */

  {
    /*! [inlet_vel_analytic] */
    cs_equation_param_t  *eqp = cs_equation_param_by_name("velocity");

    cs_equation_add_bc_by_analytic(eqp,
                                   CS_PARAM_BC_DIRICHLET,
                                   "inlet",           // zone name
                                   _vel_profile,      // callback function
                                   NULL);             // input structure
    /*! [inlet_vel_analytic] */
  }

  /* Example: define inlet scalar value using user-defined callback function */

  {
    /*! [inlet_scalar_analytic] */
    cs_field_t  *f = cs_field_by_name("scalar1");
    cs_equation_param_t  *eqp = cs_equation_param_by_name("scalar1");

    cs_equation_add_bc_by_analytic(eqp,
                                   CS_PARAM_BC_DIRICHLET,
                                   "inlet",                 // zone name
                                   _scalar_inlet_profile,   // callback function
                                   f);                      // input structure
    /*! [inlet_scalar_analytic] */
  }

  /* Example: define constant wall exchange coefficient + external value */

  {
    /*! [wall_scalar_exchange_const] */

    /* Robin BC: K du/dn + alpha*(u - u0) = g
       with alpha = -10.5 and u0 = 0.1 */

    cs_real_t robin_values[3] = {-10.5, 0.1, 0};

    cs_equation_param_t  *eqp = cs_equation_param_by_name("scalar1");

    cs_equation_add_bc_by_value(eqp,
                                CS_PARAM_BC_ROBIN,
                                "exchanger_wall",         // zone name
                                robin_values);
    /*! [wall_scalar_exchange_const] */
  }

  /* Example: define wall exchange coefficient + external value
     with user-defined callback function. */
  {
    /*! [wall_scalar_exchange_dof] */
    cs_field_t  *f = cs_field_by_name("scalar1");
    cs_equation_param_t  *eqp = cs_equation_param_by_name("scalar1");

    cs_equation_add_bc_by_dof_func(eqp,
                                   CS_PARAM_BC_ROBIN,
                                   "exchanger_wall",         // zone name
                                   cs_flag_boundary_face,    // location flag
                                   _scalar_exchange_profile, // callback function
                                   f);                       // input structure
    /*! [wall_scalar_exchange_dof] */
  }

  /* Example: define wall flux based on total flux on zone. */

  {
    /*! [wall_top_flux_dof] */
    cs_equation_param_t  *eqp = cs_equation_param_by_name("scalar_1");

    cs_equation_add_bc_by_dof_func(eqp,
                                   CS_PARAM_BC_NEUMANN,
                                   "wall_top",             // zone name
                                   cs_flag_boundary_face,  // location flag
                                   _w_flux_top,            // callback function
                                   NULL);                  // input structure
    /*! [wall_top_flux_dof] */
  }

  /* Example: define wall flux based on total flux on zone (variant) */

  {
    /*! [wall_side_flux_dof] */
    const cs_zone_t *z = cs_boundary_zone_by_name("wall_side");
    cs_equation_param_t  *eqp = cs_equation_param_by_name("scalar_1");

    cs_equation_add_bc_by_dof_func(eqp,
                                   CS_PARAM_BC_NEUMANN,
                                   z->name,                // zone name
                                   cs_flag_boundary_face,  // location flag
                                   _w_flux_side,           // callback function
                                   (void *)z);             // input structure
    /*! [wall_side_flux_dof] */
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
