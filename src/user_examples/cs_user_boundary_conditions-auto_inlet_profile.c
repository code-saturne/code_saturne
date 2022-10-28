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
  /*! [init] */
  const cs_lnum_t n_cells_ext = domain->mesh->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces = domain->mesh->n_b_faces;

  const int n_fields = cs_field_n_fields();
  const int keysca = cs_field_key_id("scalar_id");

  const int nt_cur = domain->time_step->nt_cur;

  const cs_lnum_t *b_face_cells = domain->mesh->b_face_cells;

  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)domain->mesh_quantities->b_face_normal;
  const cs_real_t *restrict b_face_surf
    = ( const cs_real_t *restrict)domain->mesh_quantities->b_face_surf;

  const cs_real_t viscl0 = cs_glob_fluid_properties->viscl0;

  const cs_real_t *bpro_rho = CS_F_(rho_b)->val;
  const cs_real_3_t *cvar_vel = (const cs_real_3_t *)CS_F_(vel)->val;

  const cs_zone_t  *zn = NULL;
  /*! [init] */

  /* Assign a pseudo-periodic channel type inlet to a set of boundary faces.

   * For each subset:
   *   - use selection criteria to filter boundary faces of a given subset
   *   - loop on faces from a subset
   *   - set the boundary condition for each face
   *
   * A feedback loop is used so as to progressively reach a state similar
   * to that of a periodic channel at the inlet. */

  /*! [example_1] */
  zn = cs_boundary_zone_by_name("inlet");

  const cs_real_t fmprsc = 1.; // mean prescribed velocity

  if (nt_cur == 1) {

    /* For the Rij-EBRSM model (and possibly V2f), we need a non-flat profile,
     * so as to ensure turbulent production, and avoid laminarization;
     * here, we simply divide the initial velocity by 10 for inlet
     * faces adjacent to the wall.

     * The loop below assumes wall conditions have been defined first
     * (in the GUI, or in this file, before the current test). */

    int *mrkcel = NULL;
    if ((cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM) ||
        (cs_glob_turb_model->itytur == 5)){

      BFT_MALLOC(mrkcel, n_cells_ext, int);
      for (cs_lnum_t i = 0; i < n_cells_ext; i++)
        mrkcel[i] = 0;

      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
        if (bc_type[f_id] == CS_SMOOTHWALL) {
          const cs_lnum_t c_id = b_face_cells[f_id];
          mrkcel[c_id] = 1;
        }
      }
    }

    for (cs_lnum_t e_idx = 0; e_idx < zn->n_elts; e_idx++) {

      const cs_lnum_t face_id = zn->elt_ids[e_idx];
      const cs_lnum_t c_id = b_face_cells[face_id];

      bc_type[face_id] = CS_INLET;

      for (int ii = 0; ii< CS_F_(vel)->dim; ii++)
        CS_F_(vel)->bc_coeffs->rcodcl1[n_b_faces*ii + face_id]
          = -fmprsc*b_face_normal[face_id][ii]/b_face_surf[face_id];

      if (mrkcel[c_id] == 1)
        CS_F_(vel)->bc_coeffs->rcodcl1[n_b_faces*0 + face_id] = fmprsc/10;

      cs_real_t uref2 = 0;
      for (int ii = 0; ii< CS_F_(vel)->dim; ii++)
        uref2 += pow(CS_F_(vel)->bc_coeffs->rcodcl1[n_b_faces*ii + face_id], 2);
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

      /* Hydraulic diameter */
      cs_real_t xdh = 1.0;

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
        cs_field_t *fld = cs_field_by_id(f_id);

        /* Here we only handle user scalar */
        int sc_id = cs_field_get_key_int(fld, keysca) - 1;
        if (sc_id < 0)
          continue;

        fld->bc_coeffs->rcodcl1[face_id] = 1;
      }
    }
    BFT_FREE(mrkcel);
  }
  else {

    /* Subsequent time steps
     *----------------------*/

    cs_real_2_t acc = {0, 0};

    /* Estimate multiplier */

    for (cs_lnum_t e_idx = 0; e_idx < zn->n_elts; e_idx++) {

      const cs_lnum_t face_id = zn->elt_ids[e_idx];
      const cs_lnum_t c_id = b_face_cells[face_id];

      const cs_real_t vnrm = cs_math_3_norm(cvar_vel[c_id]);

      acc[0] += vnrm*b_face_surf[face_id];
      acc[1] += b_face_surf[face_id];
    }

    cs_parall_sum(2, CS_REAL_TYPE, acc);

    cs_real_t fmul = 0; /* zero velocity in bulk domain */
    if (acc[0] > cs_math_epzero)
      fmul = fmprsc/(acc[0]/acc[1]); /* estimate flow multiplier */

    /* Apply BC */
    for (cs_lnum_t e_idx = 0; e_idx < zn->n_elts; e_idx++) {

      const cs_lnum_t face_id = zn->elt_ids[e_idx];
      const cs_lnum_t c_id = b_face_cells[face_id];

      const cs_real_t vnrm = cs_math_3_norm(cvar_vel[c_id]);

      bc_type[face_id] = CS_INLET;

      for (int ii = 0; ii< CS_F_(vel)->dim; ii++)
        CS_F_(vel)->bc_coeffs->rcodcl1[n_b_faces*ii + face_id]
          = -fmul*vnrm*b_face_normal[face_id][ii]/b_face_surf[face_id];

      if (cs_glob_turb_model->itytur == 2) {
        CS_F_(k)->bc_coeffs->rcodcl1[face_id] = CS_F_(k)->val[c_id];
        CS_F_(eps)->bc_coeffs->rcodcl1[face_id] = CS_F_(eps)->val[c_id];
      }

      else if (cs_glob_turb_model->itytur == 3) {
        for (cs_lnum_t ii = 0; ii< CS_F_(rij)->dim; ii++)
          CS_F_(rij)->bc_coeffs->rcodcl1[n_b_faces*ii + face_id]
            = CS_F_(rij)->val[c_id];

        CS_F_(eps)->bc_coeffs->rcodcl1[face_id] = CS_F_(eps)->val[c_id];

        if (cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM)
          CS_F_(alp_bl)->bc_coeffs->rcodcl1[face_id] = CS_F_(alp_bl)->val[c_id];
      }

      else if (cs_glob_turb_model->itytur == 5) {

        CS_F_(k)->bc_coeffs->rcodcl1[face_id] = CS_F_(k)->val[c_id];
        CS_F_(eps)->bc_coeffs->rcodcl1[face_id] = CS_F_(eps)->val[c_id];
        CS_F_(phi)->bc_coeffs->rcodcl1[face_id] = CS_F_(phi)->val[c_id];

        if (cs_glob_turb_model->iturb == CS_TURB_V2F_PHI)
          CS_F_(f_bar)->bc_coeffs->rcodcl1[face_id] = CS_F_(f_bar)->val[c_id];
        else if (cs_glob_turb_model->iturb == CS_TURB_V2F_BL_V2K)
          CS_F_(alp_bl)->bc_coeffs->rcodcl1[face_id] = CS_F_(alp_bl)->val[c_id];

      }

      else if (cs_glob_turb_model->iturb == CS_TURB_K_OMEGA) {

        CS_F_(k)->bc_coeffs->rcodcl1[face_id] = CS_F_(k)->val[c_id];
        CS_F_(omg)->bc_coeffs->rcodcl1[face_id] = CS_F_(omg)->val[c_id];
      }

      else if (cs_glob_turb_model->iturb ==  CS_TURB_SPALART_ALLMARAS) {

        CS_F_(nusa)->bc_coeffs->rcodcl1[face_id] = CS_F_(nusa)->val[c_id];

      }

      for (int f_id = 0; f_id < n_fields; f_id++) {
        cs_field_t *fld = cs_field_by_id(f_id);

        /* Here we only handle user scalar */
        int sc_id = cs_field_get_key_int(fld, keysca) - 1;
        if (sc_id < 0)
          continue;

        fld->bc_coeffs->rcodcl1[face_id] = fld->val[c_id];
      }

    }
  }
  /*!example_1*/
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
