/*============================================================================
 * This function is called at the end of each time step, and has a very
 *  general purpose
 *  (i.e. anything that does not have another dedicated user function)
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_prototypes.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_selector.h"

#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_extra_operations-scalar_balance.c
 *
 * \brief This function is called at the end of each time step, and has a very
 * general purpose (i.e. anything that does not have another dedicated
 * user function).
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Example for scalar balance.
 *----------------------------------------------------------------------------*/

void
cs_user_extra_operations(cs_domain_t     *domain)
{

  /* Local variables */
  cs_lnum_t n_faces;
  cs_lnum_t *face_list;

  int cell_id, cell_id1, cell_id2, face_id;
  int nt_cur = domain->time_step->nt_cur;

  const cs_mesh_t *m = domain->mesh;
  const cs_mesh_quantities_t *mq = domain->mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)m->i_face_cells;
  const cs_lnum_t *b_face_cells = (const cs_lnum_t *)m->b_face_cells;

  const cs_real_t *cell_vol = mq->cell_vol;
  const cs_real_3_t *diipb = (const cs_real_3_t *)mq->diipb;
  const cs_real_t *b_face_surf = (const cs_real_t *)mq->b_face_surf;

  /* Get physical fields */
  const cs_real_t *dt = CS_F_(dt)->val;
  const cs_real_t *rho = CS_F_(rho)->val;
  const cs_field_t *h = cs_field_by_name_try("enthalpy");

  /*-------------------------------------------------------------------------
   * This example computes energy balance relative to enthalpy
   * We assume that we want to compute balances (convective and diffusive)
   * at the boundaries of the calculation domain represented below
   * (with boundaries marked by colors).
   *
   * The scalar considered if the enthalpy. We will also use the
   * specific heat (to obtain balances in Joules)
   *
   *
   * Domain and associated boundary colors:
   * - 2, 3, 4, 7 : adiabatic walls
   * - 6          : wall with fixed enthalpy
   * - 1          : inlet
   * - 5          : outlet
   * - 8, 9       : symmetry
   *-------------------------------------------------------------------------*/

  /* 1. Initialization
     =================

    --> Local variables
        ---------------

    vol_balance   : volume contribution of unsteady terms
    div_balance   : volume contribution due to to term in div(rho u)
    a_wall_balance: contribution from adiabatic walls
    h_wall_balance: contribution from walls with fixed temperature
    sym_balance   : contribution from symmetry boundaries
    in_balance    : contribution from inlets
    out_balance   : contribution from outlets
    mass_i_balance: contribution from mass injections
    mass_o_balance: constribution from mass suctions
    tot_balance   : total balance */

  double vol_balance = 0.;
  double div_balance = 0.;
  double a_wall_balance = 0.;
  double h_wall_balance = 0.;
  double sym_balance = 0.;
  double in_balance = 0.;
  double out_balance = 0.;
  double mass_i_balance = 0.;
  double mass_o_balance = 0.;
  double tot_balance = 0.;

  /* If the scalar enthalpy is not computed, return */
  if (h == NULL)
    return;

  /* Boundary condition coefficient for h */
  const cs_real_t *a_H = h->bc_coeffs->a;
  const cs_real_t *b_H = h->bc_coeffs->b;
  const cs_real_t *af_H = h->bc_coeffs->af;
  const cs_real_t *bf_H = h->bc_coeffs->bf;

  /* Convective mass fluxes for inner and boundary faces */
  int iflmas = cs_field_get_key_int(h, cs_field_key_id("inner_mass_flux_id"));
  const cs_real_t *i_mass_flux = cs_field_by_id(iflmas)->val;

  int iflmab = cs_field_get_key_int(h, cs_field_key_id("boundary_mass_flux_id"));
  const cs_real_t *b_mass_flux = cs_field_by_id(iflmab)->val;

  /* Allocate temporary array */
  cs_real_t *h_reconstructed;
  BFT_MALLOC(h_reconstructed, n_b_faces, cs_real_t);

  /* Reconstructed value */
  if (false) {
    cs_real_3_t *grad;
    BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);

    cs_field_gradient_scalar(h,
                             true, /* use_previous_t */
                             1, /* inc */
                             true, /* _recompute_cocg */
                             grad);

    for (face_id = 0; face_id < n_b_faces; face_id++) {
      cell_id = b_face_cells[face_id]; // associated boundary cell
      h_reconstructed[face_id] = h->val[cell_id]
                               + grad[cell_id][0]*diipb[face_id][0]
                               + grad[cell_id][1]*diipb[face_id][1]
                               + grad[cell_id][2]*diipb[face_id][2];
    }

    BFT_FREE(grad);

  /* Non-reconstructed value */
  } else {
    for (face_id = 0; face_id < n_b_faces; face_id++) {
      cell_id = b_face_cells[face_id]; // associated boundary cell
      h_reconstructed[face_id] = h->val[cell_id];
    }
  }

  /* 2. Compute the balance at time step n
    ======================================

    --> Balance on interior volumes
        --------------------------- */

  for (cell_id = 0; cell_id < n_cells; cell_id++) {
    vol_balance += cell_vol[cell_id] * rho[cell_id]
                 * (h->val_pre[cell_id] - h->val[cell_id]);
  }

  /*
    --> Balance on all faces (interior and boundary), for div(rho u)
        ------------------------------------------------------------
   */

  for (face_id = 0; face_id < n_i_faces; face_id++) {

    cell_id1 = i_face_cells[face_id][0]; // associated boundary cell
    cell_id2 = i_face_cells[face_id][1]; // associated boundary cell

    /* Contribution to flux from the two cells of the current face
      (The cell is count only once in parallel by checking that
       the cell_id is not in the halo) */

    if (cell_id1 < n_cells)
      div_balance += i_mass_flux[face_id] * dt[cell_id1] * h->val[cell_id1];

    if (cell_id2 < n_cells)
      div_balance -= i_mass_flux[face_id] * dt[cell_id2] * h->val[cell_id2];
  }

  for (face_id = 0; face_id < n_b_faces; face_id++) {

    cell_id = b_face_cells[face_id]; // associated boundary cell

    /* Contribution to flux from the current face */
    div_balance += b_mass_flux[face_id] * dt[cell_id] * h->val[cell_id];

  }

  // TODO mass source terms and mass accumulation term
  // In case of a mass source term, add contribution from Gamma*Tn+1

  /*
    --> Balance on boundary faces
        -------------------------

    We handle different types of boundary faces separately to better
    analyze the information, but this is not mandatory. */

  /*
    Compute the contribution from walls with colors 2, 3, 4 and 7
    (adiabatic here, so flux should be 0)
  */
  BFT_MALLOC(face_list, n_b_faces, cs_lnum_t);

  cs_selector_get_b_face_list("2 or 3 or 4 or 7", &n_faces, face_list);

  for (cs_lnum_t i = 0; i < n_faces; i++) {

    face_id = face_list[i];
    cell_id = b_face_cells[face_id]; // associated boundary cell

    /* Contribution to flux from the current face
      (diffusion and convection flux, negative if incoming) */

    a_wall_balance += - b_face_surf[face_id] * dt[cell_id]
                        * (af_H[face_id] + bf_H[face_id] * h_reconstructed[face_id])
                      - b_mass_flux[face_id] * dt[cell_id]
                        * (a_H[face_id] + b_H[face_id] * h_reconstructed[face_id]);

  }

  /*
    Contribution from walls with color 6
    (here at fixed enthalpy; the convective flux should be 0)
  */
  cs_selector_get_b_face_list("6", &n_faces, face_list);

  for (cs_lnum_t i = 0; i < n_faces; i++) {

    face_id = face_list[i];
    cell_id = b_face_cells[face_id]; // associated boundary cell

    /* Contribution to flux from the current face
      (diffusion and convection flux, negative if incoming) */

    h_wall_balance += - b_face_surf[face_id] * dt[cell_id]
                        * (af_H[face_id] + bf_H[face_id] * h_reconstructed[face_id])
                      - b_mass_flux[face_id] * dt[cell_id]
                        * (a_H[face_id] + b_H[face_id] * h_reconstructed[face_id]);

  }

  /*
    Contribution from symmetries (should be 0).
  */
  cs_selector_get_b_face_list("8 or 9", &n_faces, face_list);

  for (cs_lnum_t i = 0; i < n_faces; i++) {

    face_id = face_list[i];
    cell_id = b_face_cells[face_id]; // associated boundary cell

    /* Contribution to flux from the current face
      (diffusion and convection flux, negative if incoming) */

    sym_balance += - b_face_surf[face_id] * dt[cell_id]
                     * (af_H[face_id] + bf_H[face_id] * h_reconstructed[face_id])
                   - b_mass_flux[face_id] * dt[cell_id]
                     * (a_H[face_id] + b_H[face_id] * h_reconstructed[face_id]);

  }

  /*
    Contribution from inlet (color 1, diffusion and convection flux)
  */
  cs_selector_get_b_face_list("1", &n_faces, face_list);

  for (cs_lnum_t i = 0; i < n_faces; i++) {

    face_id = face_list[i];
    cell_id = b_face_cells[face_id]; // associated boundary cell

    /* Contribution to flux from the current face
      (diffusion and convection flux, negative if incoming) */

    in_balance += - b_face_surf[face_id] * dt[cell_id]
                    * (af_H[face_id] + bf_H[face_id] * h_reconstructed[face_id])
                  - b_mass_flux[face_id] * dt[cell_id]
                    * (a_H[face_id] + b_H[face_id] * h_reconstructed[face_id]);

  }

  /*
    Contribution from outlet (color 5, diffusion and convection flux)
  */
  cs_selector_get_b_face_list("5", &n_faces, face_list);

  for (cs_lnum_t i = 0; i < n_faces; i++) {

    face_id = face_list[i];
    cell_id = b_face_cells[face_id]; // associated boundary cell

    /* Contribution to flux from the current face
      (diffusion and convection flux, negative if incoming) */

    out_balance += - b_face_surf[face_id] * dt[cell_id]
                     * (af_H[face_id] + bf_H[face_id] * h_reconstructed[face_id])
                   - b_mass_flux[face_id] * dt[cell_id]
                     * (a_H[face_id] + b_H[face_id] * h_reconstructed[face_id]);

  }

  /* Free memory */
  BFT_FREE(face_list);
  BFT_FREE(h_reconstructed);

  /* Sum of values on all ranks (parallel calculations) */

  cs_parall_sum(1, CS_DOUBLE, &vol_balance);
  cs_parall_sum(1, CS_DOUBLE, &div_balance);
  cs_parall_sum(1, CS_DOUBLE, &a_wall_balance);
  cs_parall_sum(1, CS_DOUBLE, &h_wall_balance);
  cs_parall_sum(1, CS_DOUBLE, &sym_balance);
  cs_parall_sum(1, CS_DOUBLE, &in_balance);
  cs_parall_sum(1, CS_DOUBLE, &out_balance);
  cs_parall_sum(1, CS_DOUBLE, &mass_i_balance);
  cs_parall_sum(1, CS_DOUBLE, &mass_o_balance);

  /* --> Total balance
         ------------- */

  /* We add the different contributions calculated above */

  tot_balance = vol_balance + div_balance + a_wall_balance + h_wall_balance
              + sym_balance + in_balance + out_balance + mass_i_balance
              + mass_o_balance;

  /* 3. Write the balance at time step n
    ==================================== */

  bft_printf("\n   ** Enthalpy balance **\n"
             "      ----------------\n"
             "-----------"
             "----------------------------------------------------------\n"
             "bt   Iter"
             "   Volume     Divergence  Adia Wall   Fixed_H Wall  Symmetry"
             "      Inlet       Outlet  Inj. Mass.  Suc. Mass.  Total\n"
             "bt %6i %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e "
             "%12.4e %12.4e %12.4e\n"
             "-----------"
             "----------------------------------------------------------\n",
    nt_cur, vol_balance, div_balance, a_wall_balance, h_wall_balance,
    sym_balance, in_balance, out_balance,
    mass_i_balance, mass_o_balance, tot_balance);

}

END_C_DECLS
