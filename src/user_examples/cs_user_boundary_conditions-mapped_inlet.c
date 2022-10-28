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
 * \file cs_user_boundary_conditions-mapped_inlet.c
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
  /*![init]*/

  const cs_lnum_t n_cells = domain->mesh->n_cells;
  const cs_lnum_t n_b_faces = domain->mesh->n_b_faces;
  const cs_lnum_t n_cells_ext = domain->mesh->n_cells_with_ghosts;

  const cs_lnum_t *b_face_cells = domain->mesh->b_face_cells;

  const int n_fields = cs_field_n_fields();

  const int nt_cur = domain->time_step->nt_cur;
  const int nt_max = domain->time_step->nt_max;
  const int nt_prev = domain->time_step->nt_prev;

  const cs_real_t viscl0 = cs_glob_fluid_properties->viscl0;

  const cs_real_t *bpro_rho = CS_F_(rho_b)->val;

  const cs_zone_t  *zn = NULL;

  const int keysca = cs_field_key_id("scalar_id");

  cs_real_t *vel_rcodcl1 = CS_F_(vel)->bc_coeffs->rcodcl1;

  /*![init]*/

  /* Assign a pseudo-periodic channel type inlet to a set of boundary faces.

   * For each subset:
   *   - use selection criteria to filter boundary faces of a given subset
   *   - loop on faces from a subset
   *   - set the boundary condition for each face
   *
   * A feedback loop is used so as to progressively reach a state similar
   * to that of a periodic channel at the inlet. */

  /*![example_1_base]*/

  zn = cs_boundary_zone_by_name("inlet");

  const cs_real_t fmprsc = 1.; /* mean prescribed velocity */
  const cs_real_t xdh = 1.0;   /* Hydraulic diameter */

  for (cs_lnum_t e_idx = 0; e_idx < zn->n_elts; e_idx++) {

    const cs_lnum_t face_id = zn->elt_ids[e_idx];
    const cs_lnum_t c_id = b_face_cells[face_id];

    bc_type[face_id] = CS_INLET;

    vel_rcodcl1[n_b_faces*0 + face_id] = -fmprsc;

    cs_real_t uref2 = 0;
    for (cs_lnum_t ii = 0; ii< CS_F_(vel)->dim; ii++)
      uref2 += cs_math_pow2(vel_rcodcl1[n_b_faces*ii + face_id]);
    uref2 = cs_math_fmax(uref2, 1e-12);

    /* Turbulence example computed using equations valid for a pipe.

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

    /* Calculation of turbulent inlet conditions using
       the turbulence intensity and standard laws for a circular pipe
       (their initialization is not needed here but is good practice) */

    cs_real_t b_rho = bpro_rho[c_id];

    cs_turbulence_bc_inlet_hyd_diam(face_id,
                                    uref2,
                                    xdh,
                                    b_rho,
                                    viscl0);

    for (int f_id = 0; f_id < n_fields; f_id++) {
      cs_field_t *f = cs_field_by_id(f_id);

      /* Here we only handle user scalar */
     if (cs_field_get_key_int(f, keysca) > 1)
       f->bc_coeffs->rcodcl1[face_id] = 1;
    }
  }
  /*![example_1_base]*/

  /* Create locator at initialization */

  /*![example_1_map_init]*/

  ple_locator_t *inlet_l = NULL;
  zn = cs_boundary_zone_by_name("inlet");

  if (nt_cur == nt_prev) {

    cs_real_t coord_shift[1][3] = {{5.95, 0, 0}};

    cs_lnum_t *cells_ids;
    BFT_MALLOC(cells_ids, n_cells_ext, cs_lnum_t);
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      cells_ids[c_id] = c_id;

    inlet_l = cs_boundary_conditions_map(CS_MESH_LOCATION_CELLS,
                                         n_cells,
                                         zn->n_elts,
                                         cells_ids,
                                         zn->elt_ids,
                                         coord_shift,
                                         0,
                                         0.10);
    BFT_FREE(cells_ids);
  }

  /*![example_1_map_init]*/

  /* Subsequent time steps
   *----------------------*/

  /*![example_1_map_apply] */
  if (nt_cur == 1) {

    const int interpolate = 0;

    for (int f_id = 0; f_id < n_fields; f_id++) {

      cs_field_t *f = cs_field_by_id(f_id);

      if (!(f->type & CS_FIELD_VARIABLE))
        continue;

      cs_real_t normalize = 0;
      const int sc_id = cs_field_get_key_int(f, keysca) - 1;
      if ((f == CS_F_(vel)) || (sc_id > -1))
        normalize = 1;

      cs_boundary_conditions_mapped_set(f,
                                        inlet_l,
                                        CS_MESH_LOCATION_CELLS,
                                        normalize,
                                        interpolate,
                                        zn->n_elts,
                                        zn->elt_ids,
                                        NULL);
    }
  }
  /*![example_1_map_apply] */

  /* Destroy locator at end */

  /*![example_1_map_free]*/
  if (nt_cur == nt_max)
    inlet_l = ple_locator_destroy(inlet_l);
  /*![example_1_map_free]*/
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
