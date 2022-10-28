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
 * \file cs_user_boundary_conditions-advanced..c
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
  /*! [loc_var_dec] */
  const cs_lnum_t n_b_faces = domain->mesh->n_b_faces;

  const int n_fields = cs_field_n_fields();

  /* Example of specific boundary conditions fully defined by the user,
   * on the basis of wall conditions.
   *
   * We prescribe for zone 'wall_s' a wall, with in addition:
   *   - a Dirichlet condition on velocity (sliding wall with no-slip condition)
   *   - a Dirichlet condition on the first scalar. */

  cs_zone_t  *zn = NULL;

  cs_field_t *scal = cs_field_by_name("scalar1");
  /*! [loc_var_dec] */

  /*! [example_1] */
  zn = cs_boundary_zone_by_name("wall_s");

  cs_real_t *vel_rcodcl1 = CS_F_(vel)->bc_coeffs->rcodcl1;
  cs_real_t *vel_rcodcl2 = CS_F_(vel)->bc_coeffs->rcodcl2;
  cs_real_t *vel_rcodcl3 = CS_F_(vel)->bc_coeffs->rcodcl3;

  for (cs_lnum_t ilelt = 0; ilelt < zn->n_elts; ilelt++) {

    const cs_lnum_t face_id = zn->elt_ids[ilelt];

    bc_type[face_id] = CS_SMOOTHWALL;

    scal->bc_coeffs->icodcl[face_id] = 1;
    CS_F_(vel)->bc_coeffs->icodcl[face_id] = 1;

    /* Dirichlet value */

    scal->bc_coeffs->rcodcl1[face_id] = 10.;

    vel_rcodcl1[n_b_faces*0 + face_id] = 1.;
    vel_rcodcl1[n_b_faces*1 + face_id] = 0.;
    vel_rcodcl1[n_b_faces*2 + face_id] = 0.;

    /* No exchange coefficient */

    scal->bc_coeffs->rcodcl2[face_id] = cs_math_infinite_r;

    vel_rcodcl2[n_b_faces*0 + face_id] = cs_math_infinite_r;
    vel_rcodcl2[n_b_faces*1 + face_id] = cs_math_infinite_r;
    vel_rcodcl2[n_b_faces*2 + face_id] = cs_math_infinite_r;

    /* Flux density at 0 */

    scal->bc_coeffs->rcodcl3[face_id] = 0;

    vel_rcodcl3[n_b_faces*0 + face_id] = 0;
    vel_rcodcl3[n_b_faces*1 + face_id] = 0;
    vel_rcodcl3[n_b_faces*2 + face_id] = 0;

  }
  /*! [example_1] */

  /* Example of specific boundary conditions fully defined by the user,
   * with no definition of a specific type.
   * We prescribe at zone 'surf_h' a homogeneous Neumann condition for
   * all variables. */

  /*! [example_2] */
  zn = cs_boundary_zone_by_name("surf_h");

  for (cs_lnum_t ilelt = 0; ilelt < zn->n_elts; ilelt++) {

    const cs_lnum_t face_id = zn->elt_ids[ilelt];

    /* CAUTION: the value of bc_type must be assigned to CS_INDEF */
    bc_type[face_id] = CS_INDEF;

    for (int f_id = 0; f_id < n_fields; f_id++) {

      cs_field_t *f = cs_field_by_id(f_id);
      if (!(f->type & CS_FIELD_VARIABLE))
        continue;

      f->bc_coeffs->icodcl[face_id] = 3;

      for (cs_lnum_t ii = 0; ii < f->dim; ii++) {
        f->bc_coeffs->rcodcl1[n_b_faces*ii + face_id] = 0.;
        f->bc_coeffs->rcodcl2[n_b_faces*ii + face_id] = cs_math_infinite_r;
        f->bc_coeffs->rcodcl3[n_b_faces*ii + face_id] = 0.;
      }
    }
  }
  /*! [example_2] */

  /* Example of wall boundary condition with automatic continuous switch
   * between rough and smooth.
   * Here the boundary_roughness is the length scale so that
   * "y+ = log (y/boundary_roughness)" in rough regime.
   * So different from the Sand grain roughness */

  /*! [example_3] */
  zn = cs_boundary_zone_by_name("r_wall");

  if (cs_field_by_name_try("boundary_roughness") != NULL) {

    cs_real_t *bpro_roughness = cs_field_by_name_try("boundary_roughness")->val;

    for (cs_lnum_t ilelt = 0; ilelt < zn->n_elts; ilelt++) {

      const cs_lnum_t face_id = zn->elt_ids[ilelt];

      bc_type[face_id] = CS_SMOOTHWALL;

      /* Boundary roughtness (in meter) */
      bpro_roughness[face_id] = 0.05;
    }
  }
  /*! [example_3] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
