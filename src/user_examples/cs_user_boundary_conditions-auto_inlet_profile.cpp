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
 * \file cs_user_boundary_conditions-auto_inlet_profile.cpp
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
  /*! [init] */
  const cs_lnum_t n_cells_ext = domain->mesh->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces = domain->mesh->n_b_faces;

  const int n_fields = cs_field_n_fields();
  const int keysca = cs_field_key_id("scalar_id");

  const int nt_cur = domain->time_step->nt_cur;

  const cs_lnum_t *b_face_cells = domain->mesh->b_face_cells;

  const cs_nreal_3_t *restrict b_face_u_normal
    = domain->mesh_quantities->b_face_u_normal;
  const cs_real_t *restrict b_face_surf
    = domain->mesh_quantities->b_face_surf;

  const cs_real_t viscl0 = cs_glob_fluid_properties->viscl0;

  const auto bpro_rho = CS_F_(rho_b)->get_val_s();
  const auto cvar_vel = CS_F_(vel)->get_val_v();

  const cs_zone_t  *zn = nullptr;
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

  constexpr cs_real_t fmprsc = 1.; // mean prescribed velocity

  auto vel_val_ext = CS_F_(vel)->bc_coeffs->get_val_ext_v();

  if (nt_cur == 1) {

    /* For the Rij-EBRSM model (and possibly V2f), we need a non-flat profile,
     * so as to ensure turbulent production, and avoid laminarization;
     * here, we simply divide the initial velocity by 10 for inlet
     * faces adjacent to the wall.

     * The loop below assumes wall conditions have been defined first
     * (in the GUI, or in this file, before the current test). */

    cs_array<int> mrkcel;
    if ((cs_glob_turb_model->model == CS_TURB_RIJ_EPSILON_EBRSM) ||
        (cs_glob_turb_model->itytur == 5)){

      mrkcel.reshape(n_cells_ext);
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

      for (int ii = 0; ii < 3; ii++)
        vel_val_ext(face_id, ii) = -fmprsc*b_face_u_normal[face_id][ii];

      if (mrkcel[c_id] == 1)
        vel_val_ext(face_id, 0) = fmprsc/10;

      cs_real_t uref2 = 0;
      for (int ii = 0; ii< CS_F_(vel)->dim; ii++)
        uref2 += cs_math_pow2(vel_val_ext(face_id, ii));
      uref2 = cs::max(uref2, 1e-12);

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

      cs_real_t b_rho = bpro_rho[face_id];

      cs_turbulence_bc_inlet_hyd_diam(face_id,
                                      uref2,
                                      xdh,
                                      b_rho,
                                      viscl0);
    }

    /* Initialize user scalars inlet */

    for (int f_id = 0; f_id < n_fields; f_id++) {
      cs_field_t *fld = cs_field(f_id);

      if (! (fld->type & (CS_FIELD_VARIABLE | CS_FIELD_USER)))
        continue;
      int sc_id = fld->get_key_int(keysca) - 1;
      if (sc_id < 0)
        continue;

      auto val_ext = fld->bc_coeffs->get_val_ext();

      for (cs_lnum_t e_idx = 0; e_idx < zn->n_elts; e_idx++) {
        const cs_lnum_t face_id = zn->elt_ids[e_idx];
        val_ext[face_id] = 1;
      }

    }
  }

  else {

    /* Subsequent time steps
     *----------------------*/

    cs_real_2_t acc = {0, 0};

    /* Accessors for BC values */

    cs_span<cs_real_t> k_val_ext;
    cs_span<cs_real_t> eps_val_ext;
    cs::mdspan<cs_real_t, 2, cs::layout::left> rij_val_ext;
    cs_span<cs_real_t> alp_bl_val_ext;
    cs_span<cs_real_t> phi_val_ext;
    cs_span<cs_real_t> f_bar_val_ext;
    cs_span<cs_real_t> omg_val_ext;
    cs_span<cs_real_t> nusa_val_ext;

    cs_span<cs_real_t> cvar_k;
    cs_span<cs_real_t> cvar_eps;
    cs_span_2d<cs_real_t> cvar_rij;
    cs_span<cs_real_t> cvar_alp_bl;
    cs_span<cs_real_t> cvar_phi;
    cs_span<cs_real_t> cvar_f_bar;
    cs_span<cs_real_t> cvar_omg;
    cs_span<cs_real_t> cvar_nusa;

    if (cs_glob_turb_model->itytur == 2) {
      k_val_ext = CS_F_(k)->bc_coeffs->get_val_ext();
      eps_val_ext = CS_F_(eps)->bc_coeffs->get_val_ext();
      cvar_k = CS_F_(k)->get_val_s();
      cvar_eps = CS_F_(eps)->get_val_s();
    }
    else if (cs_glob_turb_model->itytur == 3) {
      rij_val_ext = CS_F_(rij)->bc_coeffs->get_val_ext_t();
      eps_val_ext = CS_F_(eps)->bc_coeffs->get_val_ext();
      cvar_rij = CS_F_(rij)->get_val_t();
      cvar_eps = CS_F_(eps)->get_val_s();
      if (cs_glob_turb_model->model == CS_TURB_RIJ_EPSILON_EBRSM) {
        alp_bl_val_ext = CS_F_(alp_bl)->bc_coeffs->get_val_ext();
        cvar_alp_bl= CS_F_(alp_bl)->get_val_s();
      }
    }
    if (cs_glob_turb_model->itytur == 5) {
      k_val_ext = CS_F_(k)->bc_coeffs->get_val_ext();
      eps_val_ext = CS_F_(eps)->bc_coeffs->get_val_ext();
      phi_val_ext = CS_F_(phi)->bc_coeffs->get_val_ext();
      cvar_k = CS_F_(k)->get_val_s();
      cvar_eps = CS_F_(eps)->get_val_s();
      cvar_phi = CS_F_(phi)->get_val_s();
      if (cs_glob_turb_model->model == CS_TURB_V2F_PHI) {
        f_bar_val_ext = CS_F_(f_bar)->bc_coeffs->get_val_ext();
        cvar_f_bar = CS_F_(f_bar)->get_val_s();
      }
      else if (cs_glob_turb_model->model == CS_TURB_V2F_BL_V2K) {
        alp_bl_val_ext = CS_F_(alp_bl)->bc_coeffs->get_val_ext();
        cvar_alp_bl= CS_F_(alp_bl)->get_val_s();
      }
    }
    else if (cs_glob_turb_model->model == CS_TURB_K_OMEGA) {
      k_val_ext = CS_F_(k)->bc_coeffs->get_val_ext();
      omg_val_ext = CS_F_(omg)->bc_coeffs->get_val_ext();
      cvar_k = CS_F_(k)->get_val_s();
      cvar_omg = CS_F_(omg)->get_val_s();
    }
    else if (cs_glob_turb_model->model ==  CS_TURB_SPALART_ALLMARAS) {
      nusa_val_ext = CS_F_(nusa)->bc_coeffs->get_val_ext();
      cvar_nusa = CS_F_(nusa)->get_val_s();
    }

    /* Estimate multiplier */

    for (cs_lnum_t e_idx = 0; e_idx < zn->n_elts; e_idx++) {

      const cs_lnum_t face_id = zn->elt_ids[e_idx];
      const cs_lnum_t c_id = b_face_cells[face_id];

      const cs_real_t vnrm = cs_math_3_norm(cvar_vel.sub_array(c_id));

      acc[0] += vnrm*b_face_surf[face_id];
      acc[1] += b_face_surf[face_id];
    }

    cs_parall_sum(2, CS_REAL_TYPE, acc);

    cs_real_t fmul = 0; /* zero velocity in bulk domain */
    if (acc[0] > cs_math_epzero)
      fmul = fmprsc/(acc[0]/acc[1]); /* estimate flow multiplier */

    /* Apply BCs */
    for (cs_lnum_t e_idx = 0; e_idx < zn->n_elts; e_idx++) {
      const cs_lnum_t face_id = zn->elt_ids[e_idx];
      const cs_lnum_t c_id = b_face_cells[face_id];

      const cs_real_t vnrm = cs_math_3_norm(cvar_vel.sub_array(c_id));

      assert(bc_type[face_id] == CS_INLET);

      for (cs_lnum_t ii = 0; ii < 3; ii++)
        vel_val_ext(face_id, ii) = -fmul*vnrm*b_face_u_normal[face_id][ii];

      if (cs_glob_turb_model->itytur == 2) {
        k_val_ext[face_id] = cvar_k[c_id];
        eps_val_ext[face_id] = cvar_eps[c_id];
      }

      else if (cs_glob_turb_model->itytur == 3) {
        for (cs_lnum_t ii = 0; ii< 6; ii++)
          rij_val_ext(face_id, ii) = cvar_rij(c_id, ii);
        eps_val_ext[face_id] = cvar_eps[c_id];
        if (cs_glob_turb_model->model == CS_TURB_RIJ_EPSILON_EBRSM)
          alp_bl_val_ext[face_id] = cvar_alp_bl[c_id];
      }

      else if (cs_glob_turb_model->itytur == 5) {
        k_val_ext[face_id] = cvar_k[c_id];
        eps_val_ext[face_id] = cvar_eps[c_id];
        phi_val_ext[face_id] = cvar_phi[c_id];
        if (cs_glob_turb_model->model == CS_TURB_V2F_PHI)
          f_bar_val_ext[face_id] = cvar_f_bar[c_id];
        else if (cs_glob_turb_model->model == CS_TURB_V2F_BL_V2K)
          alp_bl_val_ext[face_id] = cvar_alp_bl[c_id];
      }

      else if (cs_glob_turb_model->model == CS_TURB_K_OMEGA) {
        k_val_ext[face_id] = cvar_k[c_id];
        omg_val_ext[face_id] = cvar_omg[c_id];
      }
      else if (cs_glob_turb_model->model ==  CS_TURB_SPALART_ALLMARAS) {
        nusa_val_ext[face_id] = cvar_nusa[c_id];
      }
    }

    for (int f_id = 0; f_id < n_fields; f_id++) {
      cs_field_t *fld = cs_field(f_id);

      if (! (fld->type & (CS_FIELD_VARIABLE | CS_FIELD_USER)))
        continue;
      int sc_id = fld->get_key_int(keysca) - 1;
      if (sc_id < 0)
        continue;

      auto val_ext = fld->bc_coeffs->get_val_ext();
      auto cvar = fld->get_val_s();

      for (cs_lnum_t e_idx = 0; e_idx < zn->n_elts; e_idx++) {
        const cs_lnum_t face_id = zn->elt_ids[e_idx];
        const cs_lnum_t c_id = b_face_cells[face_id];
        val_ext[face_id] = cvar[c_id];
      }
    }
  }
  /*! [example_1] */
}

/*----------------------------------------------------------------------------*/
