/*============================================================================
 * User definition of boundary conditions.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
/*
 * \file cs_user_boundary_conditions.c
 *
 * \brief User functions for boundary condition definitions.
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
cs_user_boundary_conditions(cs_domain_t  *domain,
                            int           bc_type[])
{
  /*! [init] */
  const cs_lnum_t n_b_faces = domain->mesh->n_b_faces;
  const cs_lnum_t *b_face_cells = domain->mesh->b_face_cells;

  const int n_fields = cs_field_n_fields();

  cs_real_t viscl0 = cs_glob_fluid_properties->viscl0;

  const cs_real_t *bpro_rho = CS_F_(rho_b)->val;

  const int keysca = cs_field_key_id("scalar_id");
  const cs_zone_t  *zn = NULL;

  cs_real_t *vel_rcodcl1 = CS_F_(vel)->bc_coeffs->rcodcl1;

  /*! [init] */

  /* Assign  boundary conditions of boundary faces.

   * For each subset:
   *   - use selection criteria to filter boundary faces of a given subset
   *   - loop on faces from a subset
   *   - set the boundary condition for each face */

  /*! [example_1] */

  zn = cs_boundary_zone_by_name("inlet_1");

  for (cs_lnum_t e_idx = 0; e_idx < zn->n_elts; e_idx++) {

    const cs_lnum_t face_id = zn->elt_ids[e_idx];

    bc_type[face_id] = CS_INLET;

    vel_rcodcl1[n_b_faces*0 + face_id] = 0;     /* vel_x */
    vel_rcodcl1[n_b_faces*1 + face_id] = 1.1;   /* vel_y */
    vel_rcodcl1[n_b_faces*2 + face_id] = 0;     /* vel_z */

    cs_real_t uref2 = 0;
    for (cs_lnum_t ii = 0; ii < 3; ii++)
      uref2 += cs_math_pow2(vel_rcodcl1[n_b_faces*ii + face_id]);
    uref2 = cs_math_fmax(uref2, 1e-12);

    /* Turbulence example computed using equations valid for a pipe.
     *
     * We will be careful to specify a hydraulic diameter adapted
     *   to the current inlet.

     * We will also be careful if necessary to use a more precise
     *   formula for the dynamic viscosity use in the calculation of
     *   the Reynolds number (especially if it is variable, it may be
     *   useful to take the law from 'cs_user_physical_properties'
     *   Here, we use by default the 'viscl0" value.
     *   Regarding the density, we have access to its value at boundary
     *   faces (b_rho) so this value is the one used here (specifically,
     *   it is consistent with the processing in 'cs_user_physical_properties',
     *   in case of variable density) */

    /* Hydraulic diameter */
    cs_real_t xdh = cs_notebook_parameter_value_by_name("hydraulic_diam");

    /* Calculation of turbulent inlet conditions using
       the turbulence intensity and standard laws for a circular pipe
       (their initialization is not needed here but is good practice) */

    cs_real_t b_rho = bpro_rho[face_id];

    cs_turbulence_bc_inlet_hyd_diam(face_id,
                                    uref2,
                                    xdh,
                                    b_rho,
                                    viscl0);
  }

  /* Set user scalar values to 1 */

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (cs_field_get_key_int(f, keysca) > 0) {
      for (cs_lnum_t e_idx = 0; e_idx < zn->n_elts; e_idx++) {
        const cs_lnum_t face_id = zn->elt_ids[e_idx];
        f->bc_coeffs->rcodcl1[face_id] = 1.0;
      }
    }
  }
  /*! [example_1] */

  /* Assign an inlet to boundary faces of zone "inlet_2" */

  /*! [example_2] */
  zn = cs_boundary_zone_by_name("inlet_2");

  /* Hydraulic diameter: to be defined in the GUI notebook. */
  cs_real_t xdh = cs_notebook_parameter_value_by_name("hydraulic_diam");

  /* Turbulence intensity: to be defined in the GUI notebook. */
  cs_real_t xitur = cs_notebook_parameter_value_by_name("turb_intensity");

  for (cs_lnum_t e_idx = 0; e_idx < zn->n_elts; e_idx++) {

    const cs_lnum_t face_id = zn->elt_ids[e_idx];

    bc_type[face_id] = CS_INLET;

    vel_rcodcl1[n_b_faces*0 + face_id] = 0;     /* vel_x */
    vel_rcodcl1[n_b_faces*1 + face_id] = 1.1;   /* vel_y */
    vel_rcodcl1[n_b_faces*2 + face_id] = 0;     /* vel_z */

    cs_real_t uref2 = 0;
    for (cs_lnum_t ii = 0; ii < 3; ii++)
      uref2 += cs_math_pow2(vel_rcodcl1[n_b_faces*ii + face_id]);
    uref2 = cs_math_fmax(uref2, 1e-12);

    /* Calculation of turbulent inlet conditions using
     *   the turbulence intensity and standard laws for a circular pipe
     *   (their initialization is not needed here but is good practice) */

    cs_turbulence_bc_inlet_turb_intensity(face_id, uref2, xitur, xdh);
  }

  /* Set user scalar values to 1 */

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (cs_field_get_key_int(f, keysca) > 0) {
      for (cs_lnum_t e_idx = 0; e_idx < zn->n_elts; e_idx++) {
        const cs_lnum_t face_id = zn->elt_ids[e_idx];
        f->bc_coeffs->rcodcl1[face_id] = 1.0;
      }
    }
  }
  /*![example_2]*/

  /* Assign an outlet to boundary faces of zone outlet */

  /*![example_3]*/
  zn = cs_boundary_zone_by_name("outlet");

  for (cs_lnum_t e_idx = 0; e_idx < zn->n_elts; e_idx++) {

    const cs_lnum_t face_id = zn->elt_ids[e_idx];
    bc_type[face_id] = CS_OUTLET;

  }
  /*![example_3]*/

  /* Assign a wall to boundary faces of zone '5' */

  /*![example_4]*/
  zn = cs_boundary_zone_by_name("5");

  /* Temperature */
  cs_field_t *th_f = cs_thermal_model_field();

  for (cs_lnum_t e_idx = 0; e_idx < zn->n_elts; e_idx++) {

    const cs_lnum_t face_id = zn->elt_ids[e_idx];
    bc_type[face_id] =  CS_SMOOTHWALL;

    if (th_f == NULL)
      continue;

    /* If temperature prescribed to 20 with wall law */
    th_f->bc_coeffs->icodcl[face_id] = 5;
    th_f->bc_coeffs->rcodcl1[face_id] = 20.;

    /* If temperature prescribed to 50 with no wall law (simple Dirichlet)
     *  with exchange coefficient 8 */
    th_f->bc_coeffs->icodcl[face_id] = 1;
    th_f->bc_coeffs->rcodcl1[face_id] = 50.;
    th_f->bc_coeffs->rcodcl2[face_id] = 8.;

    /* If flux prescribed to 4. */
    th_f->bc_coeffs->icodcl[face_id] = 3;
    th_f->bc_coeffs->rcodcl3[face_id] = 4.;
  }
  /*![example_4]*/

  /* Assign a rough wall to boundary faces of zone "7" */

  /*![example_5]*/
  cs_real_t *bpro_roughness = NULL, *bpro_roughness_t = NULL;

  if (cs_field_by_name_try("boundary_roughness") != NULL)
    bpro_roughness = cs_field_by_name_try("boundary_roughness")->val;

  if (cs_field_by_name_try("boundary_thermal_roughness") != NULL)
    bpro_roughness_t = cs_field_by_name_try("boundary_thermal_roughness")->val;

  zn = cs_boundary_zone_by_name("7");

  for (cs_lnum_t e_idx = 0; e_idx < zn->n_elts; e_idx++) {

    const cs_lnum_t face_id = zn->elt_ids[e_idx];

    /* Wall: zero flow (zero flux for pressure)
     *       rough friction for velocities (+ turbulent variables)
     *       zero flux for scalars*/

    bc_type[face_id] = CS_ROUGHWALL;

    /* Roughness for velocity: 1cm */
    if (bpro_roughness != NULL)
      bpro_roughness[face_id] = 0.01;

    /* Roughness for temperature (if required): 1cm */
    if (bpro_roughness_t != NULL)
      bpro_roughness_t[face_id] = 0.01;

    /* If sliding wall with velocity */
    CS_F_(vel)->bc_coeffs->rcodcl1[face_id] = 1.;
  }
  /*![example_5]*/

  /* Assign a symmetry to boundary faces of zone "s" */

  /*![example_6]*/
  zn = cs_boundary_zone_by_name("s");

  for (cs_lnum_t e_idx = 0; e_idx < zn->n_elts; e_idx++) {

    const cs_lnum_t face_id = zn->elt_ids[e_idx];
    bc_type[face_id] = CS_SYMMETRY;

  }
  /*![example_6]*/
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
