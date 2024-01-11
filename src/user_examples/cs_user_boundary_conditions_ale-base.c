/*============================================================================
 * User definition of boundary conditions.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
 * \brief User definition of boundary conditions for ALE
 *
 * \param[in, out]  domain       pointer to a cs_domain_t structure
 * \param[in, out]  bc_type      boundary face types
 * \param[in, out]  ale_bc_type  boundary face types for mesh velocity
 * \param[in]       impale       indicator for fixed node displacement
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
 * so for a given face "face_id" and field component "comp_id", access
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
cs_user_boundary_conditions_ale(cs_domain_t  *domain,
                                int           bc_type[],
                                int           ale_bc_type[],
                                int           impale[])
{
  /* Initialization
   * -------------- */

  /*![loc_var]*/

  const cs_lnum_t n_b_faces = domain->mesh->n_b_faces;
  const cs_lnum_t *b_face_vtx_idx = domain->mesh->b_face_vtx_idx;
  const cs_lnum_t *b_face_vtx_lst = domain->mesh->b_face_vtx_lst;
  const cs_lnum_t *b_face_cells = domain->mesh->b_face_cells;

  const int nt_cur = domain->time_step->nt_cur;

  const cs_real_t *dt = CS_F_(dt)->val;

  /* nodes displacement */
  cs_real_3_t *disale
    = (cs_real_3_t*)cs_field_by_name("mesh_displacement")->val;

  const cs_zone_t  *zn = NULL;

  /*![loc_var]*/

  /* Assign boundary conditions to boundary faces here

   *     One may use selection criteria to filter boundary case subsets
   *       Loop on faces from a subset
   *         Set the boundary condition for each face */

  /*![example_1]*/

  /* Example: For boundary faces of zone 'fv' assign a fixed velocity */
  zn = cs_boundary_zone_by_name("fv");

  cs_field_t *mesh_u = CS_F_(mesh_u);

  /* Calculation of displacement at current time step */
  const cs_real_t deltaa = sin(3.141596*(nt_cur-1)/50);
  const cs_real_t delta  = sin(3.141596*nt_cur/50.);

  for (cs_lnum_t ilelt = 0; ilelt < zn->n_elts; ilelt++) {

    const cs_lnum_t face_id = zn->elt_ids[ilelt];
    const cs_lnum_t c_id = b_face_cells[face_id];

    ale_bc_type[face_id] = CS_BOUNDARY_ALE_IMPOSED_VEL;
    mesh_u->bc_coeffs->rcodcl1[n_b_faces*0 + face_id] = 0;
    mesh_u->bc_coeffs->rcodcl1[n_b_faces*1 + face_id] = 0;
    mesh_u->bc_coeffs->rcodcl1[n_b_faces*2 + face_id] = (delta-deltaa)/dt[c_id];

  }
  /*![example_1]*/

  /* Example: for boundary faces zone "fd" assign a fixed displacement on nodes */

  /*![example_2]*/
  zn = cs_boundary_zone_by_name("fd");

  for (cs_lnum_t ilelt = 0; ilelt < zn->n_elts; ilelt++) {

    const cs_lnum_t face_id = zn->elt_ids[ilelt];

    for (cs_lnum_t ii = b_face_vtx_idx[face_id];
         ii < b_face_vtx_idx[face_id+1];
         ii++) {
      const cs_lnum_t vtx_id = b_face_vtx_lst[ii];
      if (impale[vtx_id] == 0) {
        disale[vtx_id][0] = 0;
        disale[vtx_id][1] = 0;
        disale[vtx_id][2] = delta;
        impale[vtx_id] = 1;
      }
    }
  }

  /*![example_2]*/

  /* Example: For boundary faces of zone "sb" assign a sliding boundary */

  /*![example_3]*/
  zn = cs_boundary_zone_by_name("sb");

  for (cs_lnum_t ilelt = 0; ilelt < zn->n_elts; ilelt++) {

    const cs_lnum_t face_id = zn->elt_ids[ilelt];

    ale_bc_type[face_id] = CS_BOUNDARY_ALE_SLIDING;

  }
  /*![example_3]*/

  /* Example: prescribe a fixed boundary for zone "fixed" */
  /*![example_4]*/

  zn = cs_boundary_zone_by_name("fixed");

  for (cs_lnum_t ilelt = 0; ilelt < zn->n_elts; ilelt++) {

    const cs_lnum_t face_id = zn->elt_ids[ilelt];

    ale_bc_type[face_id] = CS_BOUNDARY_ALE_FIXED;

  }
  /*![example_4]*/
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
