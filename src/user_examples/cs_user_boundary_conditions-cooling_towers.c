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

/*=============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of boundary conditions
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in, out]  bc_type  boundary face types
 *
 * The icodcl and rcodcl arrays are pre-initialized based on default
 * and GUI-defined definitions, and may be modified here.
 *
 * For a given variable field f, and a given face "face_id", these arrays
 * may be used as follows:
 *
 * - Boundary condition type code given at:
 *   f->bc_coeffs->icodcl[face_id]
 *
 * - Dirichlet value defined at:
 *   f->bc_coeffs->rcodcl1[face_id]
 *
 * - Interior exchange coefficient (infinite if no exchange) at:
 *   f->bc_coeffs->rcodcl2[face_id]
 *
 * - Flux density defined at:
 *   f->bc_coeffs->rcodcl3[face_id]
 *
 * For vector or tensor fields, these arrays are not interleaved,
 * so for a given face "face_id" and field component "comp_id", acess
 * is as follows (where n_b_faces is domain->mesh->n_b_faces):
 *
 *   f->bc_coeffs->rcodcl1[n_b_faces*comp_id + face_id]
 *   f->bc_coeffs->rcodcl2[n_b_faces*comp_id + face_id]
 *   f->bc_coeffs->rcodcl3[n_b_faces*comp_id + face_id]
 *
 * Only the icodcl code values from the first component are used in the case
 * of vector or tensor fields, so the icodcl values can be defined as for
 * a scalar.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_boundary_conditions(cs_domain_t  *domain,
                            int           bc_type[])
{
  const cs_real_3_t *b_face_cog
    = (const cs_real_3_t *)domain->mesh_quantities->b_face_cog;

  const cs_real_t ro0 = cs_glob_fluid_properties->ro0;
  const cs_real_t *gravity = cs_glob_physical_constants->gravity;
  const cs_real_t *xyzp0 = (const cs_real_t *)cs_glob_fluid_properties->xyzp0;

  const cs_zone_t  *zn = NULL;

  /* Assign a free outlet for faces of group "OUTLET" */

  /*![example_1]*/
  /* Direct pointers (without multiple indirections) */

  int       *p_icodcl  = CS_F_(p)->bc_coeffs->icodcl;
  cs_real_t *p_rcodcl1 = CS_F_(p)->bc_coeffs->rcodcl1;

  /* Set BC's over zone. */

  zn = cs_boundary_zone_by_name("OUTLET");

  for (cs_lnum_t ilelt = 0; ilelt < zn->n_elts; ilelt++) {

    /* outlet: zero flux for velocity and temperature, prescribed pressure
     *         note that pressure will be set to P0 on the free outlet face
     *         (CS_OUTLET) closest to xyz0 */
    const cs_lnum_t face_id = zn->elt_ids[ilelt];

    bc_type[face_id] = CS_OUTLET;

    /* Precribe a pressure profile for all faces
     * Warning: the pressure has to be specified in term of TOTAL pressure
     * i.e. including ro0, gravity ... */

    p_icodcl[face_id] = 1;
    p_rcodcl1[face_id]
      = ro0 * cs_math_3_distance_dot_product(xyzp0,
                                             b_face_cog[face_id],
                                             gravity);
  }
  /*![example_1]*/
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
