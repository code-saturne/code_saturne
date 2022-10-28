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
 * \file cs_user_boundary_conditions-richards.c
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
  const cs_lnum_t n_b_faces = domain->mesh->n_b_faces;
  const cs_lnum_t *b_face_cells = domain->mesh->b_face_cells;

  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)domain->mesh_quantities->b_face_cog;
  const cs_real_t *restrict b_face_surf
    = ( const cs_real_t *restrict)domain->mesh_quantities->b_face_surf;

  const cs_real_t t_cur = domain->time_step->t_cur;

  const int n_fields = cs_field_n_fields();
  const int keysca = cs_field_key_id("scalar_id");

  const cs_real_t *dt = CS_F_(dt)->val;

  const cs_zone_t *zn = NULL;

  /* Boundary mass flux */
  int iflmab = cs_field_get_key_int(CS_F_(p),
                                    cs_field_key_id("boundary_mass_flux_id"));
  const cs_real_t *b_mass_flux = cs_field_by_id(iflmab)->val;

  /* For groundwater flow module, the undefined type face iindef is always used.
   * BC for the flow part are on the pressure head H=h+z */

  /* Example 1: Dirichlet on hydraulic head and solutes */

  /*! [richards_bc_ex1] */
  zn = cs_boundary_zone_by_name("FACE1");

  for (cs_lnum_t e_idx = 0; e_idx < zn->n_elts; e_idx++) {

    const cs_lnum_t face_id = zn->elt_ids[e_idx];

    /* Undefined type face */
    bc_type[face_id] = CS_INDEF;

    /* Velocity is not solved but deduced from darcy law, so no BC is required.
     * However, no flux BC are imposed for safety. */

    CS_F_(vel)->bc_coeffs->icodcl[face_id] = 3;

    for (int ii = 0; ii< CS_F_(vel)->dim; ii++)
      CS_F_(vel)->bc_coeffs->rcodcl3[n_b_faces*ii + face_id] = 0;

    /* Dirichlet BC on hydraulic head (H = h + z) to impose a constant value */
    CS_F_(p)->bc_coeffs->icodcl[face_id] = 1;
    CS_F_(p)->bc_coeffs->rcodcl1[face_id] = 10.;

    /* Dirichlet BC on centration C (g/m^3) */
    for (int f_id = 0; f_id < n_fields; f_id++) {
      cs_field_t *f = cs_field_by_id(f_id);

      /* Here we only handle user scalars */
      if (cs_field_get_key_int(f, keysca) > 0) {
        f->bc_coeffs->icodcl[face_id] = 1;
        f->bc_coeffs->rcodcl1[face_id] = 1.0;
      }
    }
  }
  /*! [richards_bc_ex1] */

  /*! [richards_bc_ex2] */
  zn = cs_boundary_zone_by_name("FACE2");

  for (cs_lnum_t e_idx = 0; e_idx < zn->n_elts; e_idx++) {

    const cs_lnum_t face_id = zn->elt_ids[e_idx];

    /* Undefined type face */
    bc_type[face_id] = CS_INDEF;

    /* Dirichlet BC on hydraulic head (H = h + z) to impose a constant
       gradient over z axis here \f$ \grad H \cdot \vect{z} = -0.1 \f$*/
    CS_F_(p)->bc_coeffs->icodcl[face_id] = 1;
    CS_F_(p)->bc_coeffs->rcodcl1[face_id] = -1.*b_face_cog[face_id][2];
  }
  /*! [richards_bc_ex2] */

  /* Neumann on hydraulic head and solutes */

  /*! [richards_bc_ex3] */
  zn = cs_boundary_zone_by_name("FACE3");

  for (cs_lnum_t e_idx = 0; e_idx < zn->n_elts; e_idx++) {

    const cs_lnum_t face_id = zn->elt_ids[e_idx];

    /* Undefined type face */
    bc_type[face_id] = CS_INDEF;

    /* Neumann BC on hydraulic head (H = h + z) to impose a gradient along
       the surface normal
       Here \f$ \grad H \cdot \vect{n} = 0 \f$ */
    CS_F_(p)->bc_coeffs->icodcl[face_id] = 3;
    CS_F_(p)->bc_coeffs->rcodcl1[face_id] = 0.;

    /* Dirichlet BC on centration C (g/m^3) */
    for (int f_id = 0; f_id < n_fields; f_id++) {
      cs_field_t *f = cs_field_by_id(f_id);

      /* Here we only handle user scalars */
      if (cs_field_get_key_int(f, keysca) > 0) {

        /* Neumann BC on centration C for boundary surface with outward
         *  or null normal flow $V_out$:
         * It allows to impose a gradient along the surface normal for a
         * diffusive flux as the
         * convective flux is defined by the $C*V_out$.
         * Here \f$ \grad C \cdot \vect{n} = 0 \f$ */
        f->bc_coeffs->icodcl[face_id] = 3;
        f->bc_coeffs->rcodcl3[face_id] = 0.;

      }
    }
  }
  /*! [richards_bc_ex3] */

  /* Mixed (or Robin) BC to impose a total flux (diffusive and convective flux)*/

  /*! [richards_bc_ex4] */
  zn = cs_boundary_zone_by_name("FACE4");

  for (cs_lnum_t e_idx = 0; e_idx < zn->n_elts; e_idx++) {

    const cs_lnum_t face_id = zn->elt_ids[e_idx];
    const cs_lnum_t c_id = b_face_cells[face_id];

    /* Undefined type face */
    bc_type[face_id] = CS_INDEF;

    /* Dirichlet BC on hydraulic head (H = h + z) to impose a constant value */
    CS_F_(p)->bc_coeffs->icodcl[face_id] = 1;
    CS_F_(p)->bc_coeffs->rcodcl1[face_id] = 10.;

    /* For impose a radioactive activity R_act (Bq) at a surface S_act (m^2)
       during a period of time T_act, a mixed (or Robin) BC is used.
       We need two quantities:
         the velocity at an entrance is Vel = - mass_flux / surf
         or Vel = mass_flux / surf at an exit
         the reference concentration Cref = R_act / (S_act * dt * Vel) */

    const cs_real_t R_act = 0.5;
    const cs_real_t S_act = 50.0;
    const cs_real_t T_act = 10.0;

    /* The total flux is imposed from 0 to T_act */
    const cs_real_t Vel = -b_mass_flux[face_id]/b_face_surf[face_id];
    cs_real_t  Cref = R_act / (S_act * dt[c_id] * Vel);

    if (t_cur > T_act)
      Cref = 0.;

    for (int f_id = 0; f_id < n_fields; f_id++) {
      cs_field_t *f = cs_field_by_id(f_id);

      /* Here we only handle user scalar */
      if (cs_field_get_key_int(f, keysca) > 0) {
        f->bc_coeffs->icodcl[face_id] = 1;
        f->bc_coeffs->rcodcl1[face_id] = Cref;
        f->bc_coeffs->rcodcl2[face_id] = Vel;
      }
    }
  }
  /*! [richards_bc_ex4] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
