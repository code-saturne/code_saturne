/*============================================================================
 * User definition of boundary conditions.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include "cs_headers.h"

/*----------------------------------------------------------------------------
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*
 * \file cs_user_boundary_conditions-advanced.cpp
 *
 * \brief User functions for boundary condition definitions.
 */
/*----------------------------------------------------------------------------*/

/*=============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * User definition of boundary conditions.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in, out]  bc_type  boundary face types
 */
/*----------------------------------------------------------------------------*/

void
cs_user_boundary_conditions([[maybe_unused]] cs_domain_t  *domain,
                            [[maybe_unused]] int           bc_type[])
{
  /*! [loc_var_dec] */
  const int n_fields = cs_field_n_fields();

  /* Example of specific boundary conditions fully defined by the user,
   * on the basis of wall conditions.
   *
   * We prescribe for zone 'wall_s' a wall, with in addition:
   *   - a Dirichlet condition on velocity (sliding wall with no-slip condition)
   *   - a Dirichlet condition on the first scalar. */

  cs_field_t *scal = cs_field("scalar1");
  /*! [loc_var_dec] */

  /*! [example_1] */
  {
    const cs_zone_t *zn = cs_boundary_zone_by_name("wall_s");

    const cs_lnum_t n_elts = zn->n_elts;
    const cs_lnum_t *elt_ids = zn->elt_ids;

    int *vel_icodcl = CS_F_(vel)->bc_coeffs->icodcl;
    auto vel_val_ext = CS_F_(vel)->bc_coeffs->get_val_ext_v();
    auto vel_h_ext = CS_F_(vel)->bc_coeffs->get_h_ext_v();
    auto vel_q_ext = CS_F_(vel)->bc_coeffs->get_q_ext_v();

    int *scal_icodcl = scal->bc_coeffs->icodcl;
    auto scal_val_ext = scal->bc_coeffs->get_val_ext();
    auto scal_h_ext = scal->bc_coeffs->get_h_ext();
    auto scal_q_ext = scal->bc_coeffs->get_q_ext();

    for (cs_lnum_t elt_idx = 0; elt_idx < n_elts; elt_idx++) {

      const cs_lnum_t face_id = elt_ids[elt_idx];

      bc_type[face_id] = CS_SMOOTHWALL;

      scal_icodcl[face_id] = CS_BC_DIRICHLET;
      vel_icodcl[face_id] = CS_BC_DIRICHLET;

      /* Dirichlet value */

      scal_val_ext[face_id] = 10.;

      vel_val_ext(face_id, 0) = 1.;
      vel_val_ext(face_id, 1) = 0.;
      vel_val_ext(face_id, 2) = 0.;

      /* No exchange coefficient */

      scal_h_ext[face_id] = cs_math_infinite_r;

      vel_h_ext(face_id, 0) = cs_math_infinite_r;
      vel_h_ext(face_id, 1) = cs_math_infinite_r;
      vel_h_ext(face_id, 2) = cs_math_infinite_r;

      /* Flux density at 0 */

      scal_q_ext[face_id] = 0;

      vel_q_ext(face_id, 0) = 0;
      vel_q_ext(face_id, 1) = 0;
      vel_q_ext(face_id, 1) = 0;
    }
  }
  /*! [example_1] */

  /* Example of specific boundary conditions fully defined by the user,
   * with no definition of a specific type.
   * We prescribe at zone 'surf_h' a homogeneous Neumann condition for
   * all variables. */

  /*! [example_2] */

  {
    const cs_zone_t *zn = cs_boundary_zone_by_name("surf_h");

    const cs_lnum_t n_elts = zn->n_elts;
    const cs_lnum_t *elt_ids = zn->elt_ids;

    /* CAUTION: the value of bc_type must be assigned to CS_UNDEF */
    for (cs_lnum_t elt_idx = 0; elt_idx < n_elts; elt_idx++) {
      const cs_lnum_t face_id = elt_ids[elt_idx];
      bc_type[face_id] = CS_UNDEF;
    }

    for (int f_id = 0; f_id < n_fields; f_id++) {

      cs_field_t *f = cs_field(f_id);
      if (!(f->type & CS_FIELD_VARIABLE))
        continue;

      int *icodcl = f->bc_coeffs->icodcl;
      auto val_ext = f->bc_coeffs->get_val_ext_2d();
      auto h_ext = f->bc_coeffs->get_h_ext_2d();
      auto q_ext = f->bc_coeffs->get_q_ext_2d();

      for (cs_lnum_t elt_idx = 0; elt_idx < n_elts; elt_idx++) {
        const cs_lnum_t face_id = elt_ids[elt_idx];

        icodcl[face_id] = CS_BC_NEUMANN;
        for (cs_lnum_t ii = 0; ii < f->dim; ii++) {
          val_ext(face_id, ii) = 0.;
          h_ext(face_id, ii) = cs_math_infinite_r;
          q_ext(face_id, ii) = 0.;
        }
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
  {
    const cs_zone_t *zn = cs_boundary_zone_by_name("r_wall");

    const cs_lnum_t n_elts = zn->n_elts;
    const cs_lnum_t *elt_ids = zn->elt_ids;

    if (cs_field_try("boundary_roughness") != nullptr) {

      auto bpro_roughness = cs_field("boundary_roughness")->get_val_s();

      for (cs_lnum_t elt_idx = 0; elt_idx < n_elts; elt_idx++) {

        const cs_lnum_t face_id = elt_ids[elt_idx];

        bc_type[face_id] = CS_SMOOTHWALL;

        /* Boundary roughtness (in meter) */
        bpro_roughness[face_id] = 0.05;
      }
    }
  }
  /*! [example_3] */
}

/*----------------------------------------------------------------------------*/
