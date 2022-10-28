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
 * \file cs_user_boundary_conditions-atmospheric.c
 *
 * \brief Atmospheric example of cs_user_boundary_conditions function.
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
  /* For boundary faces of zone "inlet_3",
   * assign an inlet boundary condition.
   * Here, all other variables prescribed from the meteo profile
   * are assumed to be managed by the GUI, except for
   * dynamic variables which are prescribed with a rough log law */

  /*![example_3]*/
  {
    const cs_lnum_t n_b_faces = domain->mesh->n_b_faces;
    const cs_real_3_t *restrict b_face_cog
      = (const cs_real_3_t *restrict)domain->mesh_quantities->b_face_cog;

    const cs_real_t d2o3 = 2./3;

    /* Parameters for the analytical rough wall law (neutral) */
    const cs_real_t rugd = 0.10;
    const cs_real_t zref = 10.0;
    const cs_real_t xuref = 10.0;

    const cs_zone_t *zn = cs_boundary_zone_by_name("inlet_3");

    for (cs_lnum_t e_idx = 0; e_idx < zn->n_elts; e_idx++) {

      const cs_lnum_t face_id = zn->elt_ids[e_idx];
      bc_type[face_id] = CS_INLET;

      /* Dynamic variables are prescribed with a rough log law;
         note: using functions from the `cs_turbulence_bc` series
         is preferrable when the appropriate function is available. */
      const cs_real_t zent = b_face_cog[face_id][2];

      const cs_real_t ustar = cs_turb_xkappa*xuref/log((zref+rugd)/rugd);
      const cs_real_t xuent = ustar/cs_turb_xkappa*log((zent+rugd)/rugd);
      const cs_real_t xvent = 0.0;
      const cs_real_t xkent = cs_math_pow2(ustar)/sqrt(cs_turb_cmu);
      const cs_real_t xeent = cs_math_pow3(ustar)/cs_turb_xkappa/(zent+rugd);

      if (cs_glob_turb_model->itytur == 2) {
        CS_F_(k)->bc_coeffs->rcodcl1[face_id] = xkent;
        CS_F_(eps)->bc_coeffs->rcodcl1[face_id] = xeent;
      }

      else if (cs_glob_turb_model->itytur == 3) {
        for (int ii = 0; ii< 3; ii++)
          CS_F_(rij)->bc_coeffs->rcodcl1[n_b_faces*ii + face_id] = d2o3*xkent;
        for (int ii = 3; ii< 6; ii++)
          CS_F_(rij)->bc_coeffs->rcodcl1[n_b_faces*ii + face_id] = 0;
        CS_F_(eps)->bc_coeffs->rcodcl1[face_id] = xeent;
      }

      else if (cs_glob_turb_model->iturb == CS_TURB_V2F_PHI) {
        CS_F_(k)->bc_coeffs->rcodcl1[face_id] = xkent;
        CS_F_(eps)->bc_coeffs->rcodcl1[face_id] = xeent;
        CS_F_(phi)->bc_coeffs->rcodcl1[face_id] = d2o3;
        CS_F_(f_bar)->bc_coeffs->rcodcl1[face_id] = 0.0;
      }

      else if (cs_glob_turb_model->iturb == CS_TURB_K_OMEGA) {
        CS_F_(k)->bc_coeffs->rcodcl1[face_id] = xkent;
        CS_F_(omg)->bc_coeffs->rcodcl1[face_id] = xeent/cs_turb_cmu/xkent;
      }

      else if (cs_glob_turb_model->iturb ==  CS_TURB_SPALART_ALLMARAS) {
        CS_F_(nusa)->bc_coeffs->rcodcl1[face_id]
          = cs_turb_cmu*cs_math_pow2(xkent)/xeent;
      }
    }
  }
  /*! [example_3] */

  /* Rough wall at boundary faces of zone "b_5". */

  /*! [example_4] */
  {
    /* Parameters for the analytical rough wall law (neutral) */
    const cs_real_t rugd = 0.10;

    cs_real_t *bpro_roughness = NULL;
    cs_real_t *bpro_roughness_t = NULL;

    if (cs_field_by_name_try("boundary_roughness") != NULL)
      bpro_roughness = cs_field_by_name_try("boundary_roughness")->val;

    if (cs_field_by_name_try("boundary_thermal_roughness") != NULL)
      bpro_roughness = cs_field_by_name_try("boundary_thermal_roughness")->val;

    const cs_zone_t *zn = cs_boundary_zone_by_name("b_5");

    for (cs_lnum_t e_idx = 0; e_idx < zn->n_elts; e_idx++) {

      const cs_lnum_t face_id = zn->elt_ids[e_idx];
      bc_type[face_id] = CS_ROUGHWALL;

      if (bpro_roughness != NULL)
        bpro_roughness[face_id] = rugd;

      if (bpro_roughness_t != NULL)
        bpro_roughness_t[face_id] = 0.01;
    }
  }
  /*! [example_4] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
