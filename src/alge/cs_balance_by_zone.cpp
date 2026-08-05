/*============================================================================
 * Scalar balance on zones.
 *============================================================================*/

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_array.h"
#include "base/cs_base.h"
#include "base/cs_boundary_conditions.h"
#include "base/cs_boundary_conditions_set_coeffs.h"
#include "alge/cs_convection_diffusion.h"
#include "alge/cs_convection_diffusion_priv.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_pointer.h"
#include "base/cs_field_operator.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_halo.h"
#include "base/cs_halo_perio.h"
#include "base/cs_log.h"
#include "base/cs_parall.h"
#include "base/cs_parameters.h"
#include "base/cs_post.h"
#include "base/cs_time_step.h"
#include "base/cs_turbomachinery.h"
#include "base/cs_selector.h"
#include "alge/cs_face_viscosity.h"
#include "base/cs_physical_constants.h"
#include "base/cs_thermal_model.h"
#include "base/cs_volume_mass_injection.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "alge/cs_balance_by_zone.h"

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file  cs_balance_by_zone.cpp

*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add convective flux to flux at boundary face.
 *
 * The convective flux is a pure upwind flux.
 *
 * \param[in]     bc_type      type of boundary face
 * \param[in]     pi           value at cell i
 * \param[in]     coefap       explicit boundary coefficient for convection operator
 * \param[in]     coefbp       implicit boundary coefficient for convection operator
 * \param[in]     b_massflux   mass flux at boundary face
 * \param[in,out] flux         flux at boundary face
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline void
_b_upwind_flux(const int        bc_type,
               const cs_real_t  pi,
               const cs_real_t  coefap,
               const cs_real_t  coefbp,
               const cs_real_t  b_massflux,
               cs_real_t       *flux)
{
  cs_real_t flui, fluj;

  /* Remove decentering for coupled faces */
  if (bc_type == CS_COUPLED_FD) {
    flui = 0.0;
    fluj = b_massflux;
  }
  else {
    flui = (signbit(b_massflux)) ? 0. : b_massflux;
    fluj = (signbit(b_massflux)) ? b_massflux : 0.;
  }

  cs_real_t pfac = coefap + coefbp*pi;
  *flux += (flui*pi + fluj*pfac);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add convective fluxes to fluxes at face ij.
 *
 * \param[in]     pi           value at cell i
 * \param[in]     pj           value at cell j
 * \param[in]     i_massflux   mass flux at face ij
 * \param[in,out] fluxij       fluxes at face ij
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE inline void
_i_conv_flux(const cs_real_t  pi,
             const cs_real_t  pj,
             const cs_real_t  i_massflux,
             cs_real_2_t      fluxij)
{
  cs_real_t flui = (signbit(i_massflux)) ? 0. : i_massflux;
  cs_real_t fluj = (signbit(i_massflux)) ? i_massflux : 0.;

  fluxij[0] += (flui*pi + fluj*pj);
  fluxij[1] += (flui*pi + fluj*pj);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the different terms of the balance of a given scalar,
 *        on a volume zone defined by selected cell ids/
 *
 * This function computes the balance relative to a given scalar
 * on a selected zone of the mesh.
 * We assume that we want to compute balances (convective and diffusive)
 * at the boundaries of the calculation domain represented below
 * (with different boundary types).
 *
 * In the case of the temperature, the energy balance in Joules will be
 * computed by multiplying by the specific heat.
 *
 * \param[in]     scalar_name         scalar name
 * \param[in]     n_cells_sel         number of selected cells
 * \param[in]     cell_sel_ids        ids of selected cells
 * \param[out]    balance             array of computed balance terms
 *                                    (see \ref cs_balance_term_t)
 */
/*----------------------------------------------------------------------------*/

void
cs_balance_by_zone_compute(const char      *scalar_name,
                           cs_lnum_t        n_cells_sel,
                           const cs_lnum_t  cell_sel_ids[],
                           cs_real_t        balance[CS_BALANCE_N_TERMS])
{
  const cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_halo_t  *halo = m->halo;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;

  const int *bc_type = cs_glob_bc_type;

  /* Parallel or device dispatch */
  cs_dispatch_context ctx;
  ctx.set_use_gpu(false);  /* balance_by_zone case not ported to GPU */

  /* initialize output */

  for (int i = 0; i < CS_BALANCE_N_TERMS; i++)
    balance[i] = 0;

  /* Get physical fields */
  const cs_real_t *dt = CS_F_(dt)->val;
  const cs_real_t *rho = CS_F_(rho)->val;
  cs_field_t *f = cs_field_by_name_try(scalar_name);
  const int field_id = cs_field_id_by_name(scalar_name);

  /* If the requested scalar field is not computed, return */
  if (field_id == -1) {
    bft_printf("Scalar field does not exist. Balance will not be computed.\n");
    return;
  }

  /* Get the calculation option from the field */
  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);

  /* Specific heat (CP) */
  const cs_real_t *xcpp = nullptr;
  cs_real_t *_xcpp = nullptr;
  cs_real_t cp0 = 1;
  if (f->get_key_int("is_temperature")) {
    const int icp = cs_field_id_by_name("specific_heat");
    if (icp != -1)
      xcpp = CS_F_(cp)->val;
    else
      cp0 = cs_glob_fluid_properties->cp0;
  }
  if (xcpp == nullptr) {
    CS_MALLOC_HD(_xcpp, n_cells_ext, cs_real_t, cs_alloc_mode);
    cs_arrays_set_value<cs_real_t, 1>(n_cells_ext, cp0, _xcpp);
    xcpp = _xcpp;
  }

  /* Zone cells selection variables*/
  cs_lnum_t n_i_faces_sel = 0;
  cs_lnum_t *i_face_sel_ids = nullptr;
  cs_lnum_t n_bb_faces_sel = 0;
  cs_lnum_t *bb_face_sel_ids = nullptr;
  cs_lnum_t n_bi_faces_sel = 0;
  cs_lnum_t *bi_face_sel_ids = nullptr;
  cs_lnum_2_t *bi_face_cells = nullptr;
  cs_lnum_t *cells_tag_ids = nullptr;

  /* Initialize balance contributions
    ---------------------------------

    vol_balance   : volume contribution of unsteady terms
    div_balance   : volume contribution due to to term in div(rho u)
    mass_i_balance: contribution from mass injections
    mass_o_balance: contribution from mass suctions
    bi_i_balance  : contribution from inlet boundary faces of the selected zone
                    which are internal in the total mesh
    bi_o_balance  : contribution from outlet boundary faces of the selected zone
                    which are internal in the total mesh
    in_balance    : contribution from inlets
    out_balance   : contribution from outlets
    sym_balance   : contribution from symmetry boundaries
    s_wall_balance: contribution from smooth walls
    r_wall_balance: contribution from rough walls
    cpl_balance   : contribution from coupled faces
    i_cpl_balance : contribution from internal coupled faces
    ndef_balance  : contribution from undefined faces
    tot_balance   : total balance */

  double vol_balance = 0.;
  double tot_vol_balance2 = 0.;
  double div_balance = 0.;
  double mass_i_balance = 0.;
  double mass_o_balance = 0.;
  double bi_i_balance = 0.;
  double bi_o_balance = 0.;
  double in_balance = 0.;
  double out_balance = 0.;
  double sym_balance = 0.;
  double s_wall_balance = 0.;
  double r_wall_balance = 0.;
  double cpl_balance = 0.;
  double i_cpl_balance = 0.;
  double ndef_balance = 0.;

  /* Convective mass fluxes for inner and boundary faces */
  int iflmas = f->get_key_int("inner_mass_flux_id");
  const cs_real_t *i_mass_flux = cs_field_by_id(iflmas)->val;

  int iflmab = f->get_key_int("boundary_mass_flux_id");
  const cs_real_t *b_mass_flux = cs_field_by_id(iflmab)->val;

  if (field_id != -1) {
    f = cs_field_by_id(field_id);

    /* Update bc_coeffs val_f and flux_diff if nullptr.
       (It may occur when scalars are not computed,
       for example when nt_max = 0) */

    if (f->bc_coeffs->val_f == nullptr || f->bc_coeffs->flux_diff == nullptr)
      cs_boundary_conditions_update_bc_coeff_face_values
        (ctx, f, eqp,
         true, true,
         0, nullptr, // hyd_p_flag, f_ext
         nullptr, // visel
         nullptr, // weighb
         f->val);
  }

  /* Face viscosity */
  int imvisf = eqp->imvisf;
  cs_real_t *i_visc;
  cs_real_t *b_visc;
  CS_MALLOC_HD(i_visc, n_i_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(b_visc, n_b_faces, cs_real_t, cs_alloc_mode);

  cs_real_t *c_visc = nullptr;
  CS_MALLOC_HD(c_visc, n_cells_ext, cs_real_t, cs_alloc_mode);
  const int kivisl = f->get_key_int("diffusivity_id");
  if (kivisl != -1) {
    cs_array_copy(n_cells_ext, cs_field_by_id(kivisl)->val, c_visc);
  }
  else {
    const double visls0 = f->get_key_double("diffusivity_ref");
    cs_arrays_set_value<cs_real_t, 1>(n_cells_ext, visls0, c_visc);
  }

  /* Turbulent part */
  cs_real_t *c_visct = cs_field_by_name("turbulent_viscosity")->val;

  if (eqp->idifft == 1) {
    const int ksigmas = cs_field_key_id("turbulent_schmidt");
    cs_real_t turb_schmidt = f->get_key_double(ksigmas);
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      c_visc[c_id] += xcpp[c_id] * c_visct[c_id]/turb_schmidt;
  }

  cs_face_viscosity(m, fvq, imvisf, c_visc, i_visc, b_visc);

  /* Get user-selected zone
     ====================== */

  /* Initialize arrays */

  /* Internal faces of the selected zone */
  CS_MALLOC_HD(i_face_sel_ids, n_i_faces, cs_lnum_t, cs_alloc_mode);
  /* Boundary faces of the selected zone,
     which are internal faces of the global mesh.
     Faces -> cells connectivity */
  CS_MALLOC_HD(bi_face_sel_ids, n_i_faces, cs_lnum_t, cs_alloc_mode);
  CS_MALLOC_HD(bi_face_cells, n_i_faces, cs_lnum_2_t, cs_alloc_mode);
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    i_face_sel_ids[f_id] = -1;
    bi_face_sel_ids[f_id] = -1;
    bi_face_cells[f_id][0] = -999;
    bi_face_cells[f_id][1] = -999;
  }

  /* Boundary faces of the selected zone,
     which are also boundary faces of the global mesh */
  CS_MALLOC(bb_face_sel_ids, n_b_faces, cs_lnum_t);
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    bb_face_sel_ids[f_id] = -1;
  }

  /* Synchronization for parallelism */
  CS_MALLOC_HD(cells_tag_ids, n_cells_ext, cs_lnum_t, cs_alloc_mode);
  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
    cells_tag_ids[c_id] = 0;
  }
  for (cs_lnum_t c_id = 0; c_id < n_cells_sel; c_id++) {
    cs_lnum_t c_id_sel = cell_sel_ids[c_id];
    cells_tag_ids[c_id_sel] = 1;
  }
  if (halo != nullptr) {
    cs_halo_sync_num(halo, CS_HALO_STANDARD, cells_tag_ids);
  }

  /* Classify mesh faces with respect to the selected zone */

  /* Check boundary faces:
     if they are in the selected zone, they are boundary as well */
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t c_id = b_face_cells[f_id];

    if (cells_tag_ids[c_id] == 1) {
      n_bb_faces_sel++;
      bb_face_sel_ids[n_bb_faces_sel-1] = f_id;
    }
  }

  /* Check internal faces:
     if they are in the selected zone, they can be either
     internal or boundary faces */
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    cs_lnum_t c_id1 = i_face_cells[f_id][0];
    cs_lnum_t c_id2 = i_face_cells[f_id][1];

    bool indic1 = false;
    bool indic2 = false;

    if (cells_tag_ids[c_id1] == 1)
      indic1 = true;
    if (cells_tag_ids[c_id2] == 1)
      indic2 = true;

    if (indic1 && indic2) {
      n_i_faces_sel++;
      i_face_sel_ids[n_i_faces_sel-1] = f_id;
    }
    else if (indic1 || indic2) {
      bi_face_sel_ids[n_bi_faces_sel] = f_id;
      n_bi_faces_sel++;
      /* Build the faces -> cells connectivity as done in
         i_face_cells */
      if (indic1)
        bi_face_cells[f_id][0] = c_id1;
      else
        bi_face_cells[f_id][1] = c_id2;
    }

  }

  /* Compute the balance at time step n
     ==================================

    --> Balance on interior volumes and
        total quantity on interior volumes
        ---------------------------------- */

  for (cs_lnum_t c_id = 0; c_id < n_cells_sel; c_id++) {

    cs_lnum_t c_id_sel = cell_sel_ids[c_id];

    vol_balance += cell_vol[c_id_sel] * rho[c_id_sel]
                 * xcpp[c_id_sel]
                 * (f->val_pre[c_id_sel] - f->val[c_id_sel]);

    cs_real_t rho_y_dt =  rho[c_id_sel] * xcpp[c_id_sel]
                        * f->val_pre[c_id_sel] * dt[c_id_sel];
    tot_vol_balance2 += cell_vol[c_id_sel] * rho_y_dt * rho_y_dt;
  }

  /* Balance on all faces (interior and boundary), for div(rho u)
     ------------------------------------------------------------ */

  /* Interior faces */
  for (cs_lnum_t f_id = 0; f_id < n_i_faces_sel; f_id++) {

    cs_lnum_t f_id_sel = i_face_sel_ids[f_id];
    /* Associated internal cells */
    cs_lnum_t c_id1 = i_face_cells[f_id_sel][0];
    cs_lnum_t c_id2 = i_face_cells[f_id_sel][1];

    /* Contribution to flux from the two cells of the current face
      (The cell is counted only once in parallel by checking that
       the c_id is not in the halo) */

    if (c_id1 < n_cells)
      div_balance += i_mass_flux[f_id_sel] * dt[c_id1] * f->val[c_id1]
                   * xcpp[c_id1];

    if (c_id2 < n_cells)
      div_balance -= i_mass_flux[f_id_sel] * dt[c_id2] * f->val[c_id2]
                   * xcpp[c_id2];

  }

  /* Boundary faces which are internal in the total mesh */
  for (cs_lnum_t f_id = 0; f_id < n_bi_faces_sel; f_id++) {

    cs_lnum_t f_id_sel = bi_face_sel_ids[f_id];
    /* Associated boundary-internal cells */
    cs_lnum_t c_id1 = bi_face_cells[f_id_sel][0];
    cs_lnum_t c_id2 = bi_face_cells[f_id_sel][1];

    /* Contribution to flux from the only cell of the current face
       lying inside the selected zone
      (The cell is counted only once in parallel by checking that
       the c_id is not in the halo) */

    if (c_id1 >= 0) {

      if (c_id1 < n_cells)
        div_balance += i_mass_flux[f_id_sel] * dt[c_id1] * f->val[c_id1]
                     * xcpp[c_id1];
    }

    else {

      if (c_id2 < n_cells)
        div_balance -= i_mass_flux[f_id_sel] * dt[c_id2] * f->val[c_id2]
                     * xcpp[c_id2];
    }

  }

  /* Boundary faces which are also boundary in the total mesh */
  for (cs_lnum_t f_id = 0; f_id < n_bb_faces_sel; f_id++) {

    cs_lnum_t f_id_sel = bb_face_sel_ids[f_id];
    /* Associated boundary cell */
    cs_lnum_t c_id = b_face_cells[f_id_sel];

    /* Contribution to flux from the current face */
      div_balance += b_mass_flux[f_id_sel] * dt[c_id] * f->val[c_id]
                   * xcpp[c_id];

  }

  /* Mass source terms and mass accumulation term.
     In case of a mass source term, add contribution from Gamma*Tn+1 */

  cs_lnum_t ncesmp = 0;
  const cs_lnum_t *icetsm = nullptr;
  int *itpsmp = nullptr;
  cs_real_t *smcelp, *gamma = nullptr;

  cs_volume_mass_injection_get_arrays(f, &ncesmp, &icetsm, &itpsmp,
                                      &smcelp, &gamma);

  if (ncesmp > 0) {

    for (cs_lnum_t c_idx = 0; c_idx < ncesmp; c_idx++) {
      cs_lnum_t c_id_sel = icetsm[c_idx];

      if (cells_tag_ids[c_id_sel]) {

        cs_real_t vg = gamma[c_idx];
        cs_real_t v;
        if (itpsmp[c_idx] == 0 || vg < 0)
          v = f->val[c_id_sel];
        else
          v = smcelp[c_idx];

        cs_real_t c_st = cell_vol[c_id_sel] * dt[c_id_sel]* vg * v;

        c_st *= xcpp[c_id_sel];

        if (vg < 0)
          mass_o_balance += c_st;
        else
          mass_i_balance += c_st;
      }

    }
  }

  cs_real_t *c_weight = nullptr;
  if (f->type & CS_FIELD_VARIABLE && eqp->iwgrec == 1) {
    if (eqp->idiff > 0) {
      int key_id = cs_field_key_id("gradient_weighting_id");
      int diff_id = f->get_key_int(key_id);
      if (diff_id > -1) {
        cs_field_t *weight_f = cs_field_by_id(diff_id);
        c_weight = weight_f->val;
      }
    }
  }

  cs_array_2d<cs_real_t> i_flux(n_bi_faces_sel, 2, cs_alloc_mode);
  cs_array<cs_real_t> b_flux(n_bb_faces_sel, cs_alloc_mode);
  cs_array<int> icvfli(n_bb_faces_sel, cs_alloc_mode);

  i_flux.zero();
  b_flux.zero();
  cs_arrays_set_value<int, 1>(n_bb_faces_sel, 0, icvfli.data());

  cs_convection_diffusion_scalar_at_faces(f,
                                          *eqp,
                                          0, // icvflb
                                          n_bi_faces_sel,
                                          n_bb_faces_sel,
                                          bi_face_sel_ids,
                                          bb_face_sel_ids,
                                          f->val,
                                          icvfli.data(),
                                          f->bc_coeffs,
                                          i_mass_flux,
                                          b_mass_flux,
                                          i_visc,
                                          b_visc,
                                          c_weight,
                                          xcpp,
                                          (cs_real_2_t *)i_flux.data(),
                                          b_flux.data());

  /* Balance on boundary faces
     -------------------------

     We handle different types of boundary faces separately to better
     analyze the information, but this is not mandatory. */

  for (cs_lnum_t f_id = 0; f_id < n_bb_faces_sel; f_id++) {

    cs_lnum_t f_id_sel = bb_face_sel_ids[f_id];
    /* Associated boundary cell */
    cs_lnum_t c_id = b_face_cells[f_id_sel];

    cs_real_t term_balance = b_flux(f_id);

    if (bc_type[f_id_sel] == CS_INLET ||
        bc_type[f_id_sel] == CS_CONVECTIVE_INLET ||
        bc_type[f_id_sel] == CS_FREE_INLET ||
        bc_type[f_id_sel] == CS_ESICF ||
        bc_type[f_id_sel] == CS_EPHCF)
      in_balance -= term_balance*dt[c_id];
    else if (bc_type[f_id_sel] == CS_OUTLET ||
             bc_type[f_id_sel] == CS_SSPCF ||
             bc_type[f_id_sel] == CS_SOPCF)
      out_balance -= term_balance*dt[c_id];
    else if (bc_type[f_id_sel] == CS_SYMMETRY)
      sym_balance -= term_balance*dt[c_id];
    else if (bc_type[f_id_sel] == CS_SMOOTHWALL)
      s_wall_balance -= term_balance*dt[c_id];
    else if (bc_type[f_id_sel] == CS_ROUGHWALL)
      r_wall_balance -= term_balance*dt[c_id];
    else if (   bc_type[f_id_sel] == CS_COUPLED
             || bc_type[f_id_sel] == CS_COUPLED_FD)
      cpl_balance -= term_balance*dt[c_id];
    else
      ndef_balance -= term_balance*dt[c_id];

  }

  /* Note: Balance on coupled faces (internal coupling) is treated by
   * standard BCs */

  /* Balance on boundary faces of the selected zone
     that are interior to the total mesh
     ---------------------------------------------- */

  for (cs_lnum_t f_id = 0; f_id < n_bi_faces_sel; f_id++) {

    cs_lnum_t f_id_sel = bi_face_sel_ids[f_id];
    /* Associated boundary-internal cells */
    cs_lnum_t c_id1 = i_face_cells[f_id_sel][0];
    cs_lnum_t c_id2 = i_face_cells[f_id_sel][1];

    /* (The cell is counted only once in parallel by checking that
       the c_id is not in the halo) */
    /* Face normal well oriented (check bi_face_cells array) */
    if (bi_face_cells[f_id_sel][0] >= 0) {
      if (c_id1 < n_cells) {
        if (i_mass_flux[f_id_sel] > 0)
          bi_o_balance -= i_flux(f_id, 0)*dt[c_id1];
        else
          bi_i_balance -= i_flux(f_id, 0)*dt[c_id1];
      }
    }
    /* Face normal direction reversed */
    else {
      if (c_id2 < n_cells) {
        if (i_mass_flux[f_id_sel] > 0)
          bi_i_balance += i_flux(f_id, 1)*dt[c_id2];
        else
          bi_o_balance += i_flux(f_id, 1)*dt[c_id2];
      }
    }

  }

  /* Free memory */

  CS_FREE(_xcpp);
  CS_FREE(c_visc);
  CS_FREE(i_visc);
  CS_FREE(b_visc);

  CS_FREE(cells_tag_ids);
  CS_FREE(bi_face_cells);
  CS_FREE(i_face_sel_ids);
  CS_FREE(bb_face_sel_ids);
  CS_FREE(bi_face_sel_ids);

  /* Sum of values on all ranks (parallel calculations) */

  balance[CS_BALANCE_TOTAL_NORMALIZED] = tot_vol_balance2; /* temporary */

  balance[CS_BALANCE_VOLUME] = vol_balance;
  balance[CS_BALANCE_DIV] = div_balance;
  balance[CS_BALANCE_UNSTEADY] = vol_balance + div_balance;
  balance[CS_BALANCE_MASS] = mass_i_balance + mass_o_balance;
  balance[CS_BALANCE_MASS_IN] = mass_i_balance;
  balance[CS_BALANCE_MASS_OUT] = mass_o_balance;
  balance[CS_BALANCE_INTERIOR_IN] = bi_i_balance;
  balance[CS_BALANCE_INTERIOR_OUT] = bi_o_balance;
  balance[CS_BALANCE_BOUNDARY_IN] = in_balance;
  balance[CS_BALANCE_BOUNDARY_OUT] = out_balance;
  balance[CS_BALANCE_BOUNDARY_SYM] = sym_balance;
  balance[CS_BALANCE_BOUNDARY_WALL] = s_wall_balance + r_wall_balance;
  balance[CS_BALANCE_BOUNDARY_WALL_S] = s_wall_balance;
  balance[CS_BALANCE_BOUNDARY_WALL_R] = r_wall_balance;
  balance[CS_BALANCE_BOUNDARY_COUPLED] = cpl_balance + i_cpl_balance;
  balance[CS_BALANCE_BOUNDARY_COUPLED_E] = cpl_balance;
  balance[CS_BALANCE_BOUNDARY_COUPLED_I] = i_cpl_balance;
  balance[CS_BALANCE_BOUNDARY_OTHER] = ndef_balance;

  cs_parall_sum(CS_BALANCE_N_TERMS, CS_REAL_TYPE, balance);

  /* Total balance: add the different contributions calculated above */

  balance[CS_BALANCE_TOTAL]
    =   balance[CS_BALANCE_UNSTEADY] + balance[CS_BALANCE_MASS]
      + balance[CS_BALANCE_INTERIOR_IN] + balance[CS_BALANCE_INTERIOR_OUT]
      + balance[CS_BALANCE_BOUNDARY_IN] + balance[CS_BALANCE_BOUNDARY_OUT]
      + balance[CS_BALANCE_BOUNDARY_SYM] + balance[CS_BALANCE_BOUNDARY_WALL]
      + balance[CS_BALANCE_BOUNDARY_COUPLED]
      + balance[CS_BALANCE_BOUNDARY_OTHER];

  tot_vol_balance2 = balance[CS_BALANCE_TOTAL_NORMALIZED]; /* from temporary above */
  balance[CS_BALANCE_TOTAL_NORMALIZED] = balance[CS_BALANCE_TOTAL];

  if (tot_vol_balance2 > 0.)
    balance[CS_BALANCE_TOTAL_NORMALIZED] /= sqrt(tot_vol_balance2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute and log the different terms of the balance of a given scalar,
 *        on a volumic zone defined by selection criteria.
 *        The different contributions to the balance are printed in the
 *        run_solver.log.
 *
 * This function computes the balance relative to a given scalar
 * on a selected zone of the mesh.
 * We assume that we want to compute balances (convective and diffusive)
 * at the boundaries of the calculation domain represented below
 * (with different boundary types).
 *
 * The scalar and the zone are selected at the top of the routine
 * by the user.
 * In the case of the temperature, the energy balance in Joules will be
 * computed by multiplying by the specific heat.
 *
 * \param[in]     selection_crit      zone selection criterion
 * \param[in]     scalar_name         scalar name
 */
/*----------------------------------------------------------------------------*/

void
cs_balance_by_zone(const char  *selection_crit,
                   const char  *scalar_name)
{
  cs_real_t balance[CS_BALANCE_N_TERMS];

  const cs_mesh_t *m = cs_glob_mesh;
  const int nt_cur = cs_glob_time_step->nt_cur;

  /* Select cells */

  cs_lnum_t n_cells_sel = 0;
  cs_lnum_t *cells_sel_ids = nullptr;

  CS_MALLOC(cells_sel_ids, m->n_cells, cs_lnum_t);
  cs_selector_get_cell_list(selection_crit, &n_cells_sel, cells_sel_ids);

  /* Compute balance */

  cs_balance_by_zone_compute(scalar_name,
                             n_cells_sel,
                             cells_sel_ids,
                             balance);

  CS_FREE(cells_sel_ids);

  /* Log results at time step n */

  bft_printf
    (_("   ** Scalar balance by zone at iteration %d\n"
       "   -------------------------\n"
       "------------------------------------------------------------\n"
       "   Scalar: %s\n"
       "   Zone selection criteria: \"%s\"\n"
       "------------------------------------------------------------\n"
       "   Unst. term   Inj. Mass.   Suc. Mass.\n"
       "  %12.4e %12.4e %12.4e\n"
       "------------------------------------------------------------\n"
       "   IB inlet     IB outlet\n"
       "  %12.4e %12.4e\n"
       "------------------------------------------------------------\n"
       "   Inlet        Outlet\n"
       "  %12.4e %12.4e\n"
       "------------------------------------------------------------\n"
       "   Sym.         Smooth W.    Rough W.\n"
       "  %12.4e %12.4e %12.4e\n"
       "------------------------------------------------------------\n"
       "   Coupled      Int. Coupling    Undef. BC\n"
       "  %12.4e %12.4e     %12.4e\n"
       "------------------------------------------------------------\n"
       "   Total        Instant. norm. total\n"
       "  %12.4e %12.4e\n"
       "------------------------------------------------------------\n\n"),
     nt_cur, scalar_name, selection_crit,
     balance[CS_BALANCE_UNSTEADY],
     balance[CS_BALANCE_MASS_IN], balance[CS_BALANCE_MASS_OUT],
     balance[CS_BALANCE_INTERIOR_IN], balance[CS_BALANCE_INTERIOR_OUT],
     balance[CS_BALANCE_BOUNDARY_IN], balance[CS_BALANCE_BOUNDARY_OUT],
     balance[CS_BALANCE_BOUNDARY_SYM],
     balance[CS_BALANCE_BOUNDARY_WALL_S], balance[CS_BALANCE_BOUNDARY_WALL_R],
     balance[CS_BALANCE_BOUNDARY_COUPLED_E],
     balance[CS_BALANCE_BOUNDARY_COUPLED_I],
     balance[CS_BALANCE_BOUNDARY_OTHER],
     balance[CS_BALANCE_TOTAL], balance[CS_BALANCE_TOTAL_NORMALIZED]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computes one term of the head loss balance (pressure drop) on a
 *        on a volume zone defined by selected cell ids/
 *
 * \param[in]     n_cells_sel         number of selected cells
 * \param[in]     cell_sel_ids        ids of selected cells
 * \param[out]    balance             array of computed balance terms
 *                                    (see \ref cs_balance_p_term_t)
 */
/*----------------------------------------------------------------------------*/

void
cs_pressure_drop_by_zone_compute(cs_lnum_t        n_cells_sel,
                                 const cs_lnum_t  cell_sel_ids[],
                                 cs_real_t        balance[CS_BALANCE_P_N_TERMS])
{
  const cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_halo_t  *halo = m->halo;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog = fvq->b_face_cog;

  const int *bc_type = cs_glob_bc_type;

  /* initialize output */

  for (int i = 0; i < CS_BALANCE_P_N_TERMS; i++)
    balance[i] = 0;

  /* Get physical fields */
  const cs_real_t *rho = CS_F_(rho)->val;
  const cs_field_t *f_pres = CS_F_(p);
  const cs_real_t *pressure = f_pres->val;
  const cs_field_t *f_vel = CS_F_(vel);
  const cs_real_3_t *velocity =  (const cs_real_3_t *)f_vel->val;
  const cs_real_t gravity[3] = {cs_glob_physical_constants->gravity[0],
                                cs_glob_physical_constants->gravity[1],
                                cs_glob_physical_constants->gravity[2]};

  /* Zone cells selection variables*/
  cs_lnum_t n_i_faces_sel = 0;
  cs_lnum_t *i_face_sel_ids = nullptr;
  cs_lnum_t n_bb_faces_sel = 0;
  cs_lnum_t *bb_face_sel_ids = nullptr;
  cs_lnum_t n_bi_faces_sel = 0;
  cs_lnum_t *bi_face_sel_ids = nullptr;
  cs_lnum_2_t *bi_face_cells = nullptr;
  cs_lnum_t *cells_tag_ids = nullptr;

  /* Initialization of balance contributions
     ---------------------------------------

    in_pressure   : contribution from inlets
    out_pressure  : contribution from outlets
    in_u2         : contribution from inlets
    out_u2        : contribution from outlets
    in_rhogx      : contribution from inlets
    out_rhogx     : contribution from outlets
    in_debit      : debit from inlets
    out_debit     : debit from outlets
    in_m_debit    : mass flow from inlets
    out_m_debit   : mass flow from outlets

  */

  double in_pressure= 0.;
  double out_pressure= 0.;
  double in_u2 = 0.;
  double out_u2 = 0.;
  double in_rhogx = 0.;
  double out_rhogx = 0.;
  double in_debit = 0.;
  double out_debit = 0.;
  double in_m_debit = 0.;
  double out_m_debit = 0.;

  /* Boundary condition coefficient for p */
  const cs_real_t *a_p = f_pres->bc_coeffs->a;
  const cs_real_t *b_p = f_pres->bc_coeffs->b;

  /* Boundary condition coefficient for u */
  const cs_real_3_t *a_u = (const cs_real_3_t *)f_vel->bc_coeffs->a;
  const cs_real_33_t *b_u = (const cs_real_33_t *)f_vel->bc_coeffs->b;

  /* Convective mass fluxes for inner and boundary faces */
  int iflmas = f_pres->get_key_int("inner_mass_flux_id");
  const cs_real_t *i_mass_flux = cs_field_by_id(iflmas)->val;

  int iflmab = f_pres->get_key_int("boundary_mass_flux_id");
  const cs_real_t *b_mass_flux = cs_field_by_id(iflmab)->val;

  /* Get user-selected zone
     ====================== */

  /* Initialize arrays */

  /* Internal faces of the selected zone */
  CS_MALLOC(i_face_sel_ids, n_i_faces, cs_lnum_t);
  /* Boundary faces of the selected zone,
     which are internal faces of the global mesh.
     Faces -> cells connectivity */
  CS_MALLOC(bi_face_sel_ids, n_i_faces, cs_lnum_t);
  CS_MALLOC(bi_face_cells, n_i_faces, cs_lnum_2_t);
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    i_face_sel_ids[f_id] = -1;
    bi_face_sel_ids[f_id] = -1;
    bi_face_cells[f_id][0] = -999;
    bi_face_cells[f_id][1] = -999;
  }

  /* Boundary faces of the selected zone,
     which are also boundary faces of the global mesh */
  CS_MALLOC(bb_face_sel_ids, n_b_faces, cs_lnum_t);
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    bb_face_sel_ids[f_id] = -1;
  }

  /* Synchronization for parallelism */
  CS_MALLOC(cells_tag_ids, n_cells_ext, cs_lnum_t);
  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
    cells_tag_ids[c_id] = 0;
  }

  for (cs_lnum_t c_id = 0; c_id < n_cells_sel; c_id++) {
    cs_lnum_t c_id_sel = cell_sel_ids[c_id];
    cells_tag_ids[c_id_sel] = 1;
  }

  if (halo != nullptr) {
    cs_halo_sync_num(halo, CS_HALO_STANDARD, cells_tag_ids);
  }

  /* Classify mesh faces with respect to the selected zone */

  /* Check boundary faces:
     if they are in the selected zone, they are boundary as well */
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    cs_lnum_t c_id = b_face_cells[f_id];

    if (cells_tag_ids[c_id] == 1) {
      n_bb_faces_sel++;
      bb_face_sel_ids[n_bb_faces_sel-1] = f_id;
    }
  }

  /* Check internal faces:
     if they are in the selected zone, they can be either
     internal or boundary faces */
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    cs_lnum_t c_id1 = i_face_cells[f_id][0];
    cs_lnum_t c_id2 = i_face_cells[f_id][1];

    bool indic1 = false;
    bool indic2 = false;

    if (cells_tag_ids[c_id1] == 1)
      indic1 = true;
    if (cells_tag_ids[c_id2] == 1)
      indic2 = true;

    if (indic1 && indic2) {
      n_i_faces_sel++;
      i_face_sel_ids[n_i_faces_sel-1] = f_id;
    }
    else if (indic1 || indic2) {
      n_bi_faces_sel++;
      bi_face_sel_ids[n_bi_faces_sel-1] = f_id;
      /* Build the faces -> cells connectivity as done in
         i_face_cells */
      if (indic1)
        bi_face_cells[f_id][0] = c_id1;
      else
        bi_face_cells[f_id][1] = c_id2;
    }

  }

  /* Balance computation
     =================== */

  /* Compute the balance at time step n */

  /* Balance on boundary faces
     -------------------------

     We handle different types of boundary faces separately to better
     analyze the information, but this is not mandatory. */

  for (cs_lnum_t f_id = 0; f_id < n_bb_faces_sel; f_id++) {

    cs_lnum_t f_id_sel = bb_face_sel_ids[f_id];
    /* Associated boundary cell */
    cs_lnum_t c_id = b_face_cells[f_id_sel];

    /* Pressure term FIXME rho0*gravity*(X-X0) should be added */
    cs_real_t p_rho = pressure[c_id] / rho[c_id];
    cs_real_t a_p_rho = a_p[f_id_sel] / rho[c_id];
    cs_real_t b_p_rho = b_p[f_id_sel];

    cs_real_t term_balance = 0.;

    _b_upwind_flux(bc_type[f_id_sel],
                   p_rho,
                   a_p_rho,
                   b_p_rho,
                   b_mass_flux[f_id_sel],
                   &term_balance);

    if (b_mass_flux[f_id_sel] > 0) {
      out_debit += b_mass_flux[f_id_sel]/rho[c_id];
      out_m_debit += b_mass_flux[f_id_sel];
      out_pressure += term_balance;
    }
    else {
      in_debit += b_mass_flux[f_id_sel]/rho[c_id];
      in_m_debit += b_mass_flux[f_id_sel];
      in_pressure += term_balance;
    }

    /* Kinematic term */
    cs_real_t u2 = 0.5 * cs_math_3_square_norm(velocity[c_id]);
    cs_real_t a_u2 = 0.5 * cs_math_3_square_norm(a_u[f_id_sel]);
    /* Approximation of u^2 BC */
    cs_real_t b_u2 = 1./6.*(  b_u[f_id_sel][0][0] * b_u[f_id_sel][0][0]
                            + b_u[f_id_sel][1][1] * b_u[f_id_sel][1][1]
                            + b_u[f_id_sel][2][2] * b_u[f_id_sel][2][2]);

    term_balance = 0.;

    _b_upwind_flux(bc_type[f_id_sel],
                   u2,
                   a_u2,
                   b_u2,
                   b_mass_flux[f_id_sel],
                   &term_balance);

    if (b_mass_flux[f_id_sel] > 0) {
      out_u2 += term_balance;
    }
    else {
      in_u2 += term_balance;
    }

    /* Gravity term */
    cs_real_t gx = - cs_math_3_dot_product(gravity, b_face_cog[f_id_sel]);
    /* Trivial BCs */
    cs_real_t a_gx = gx;
    cs_real_t b_gx = 0.;

    term_balance = 0.;

    _b_upwind_flux(bc_type[f_id_sel],
                   gx,
                   a_gx,
                   b_gx,
                   b_mass_flux[f_id_sel],
                   &term_balance);

    if (b_mass_flux[f_id_sel] > 0) {
      out_rhogx += term_balance;
    } else {
      in_rhogx += term_balance;
    }

  }

  /* Balance on boundary faces of the selected zone
     that are internal of the total mesh
     ---------------------------------------------- */

  for (cs_lnum_t f_id = 0; f_id < n_bi_faces_sel; f_id++) {

    cs_lnum_t f_id_sel = bi_face_sel_ids[f_id];
    /* Associated boundary-internal cells */
    cs_lnum_t c_id1 = i_face_cells[f_id_sel][0];
    cs_lnum_t c_id2 = i_face_cells[f_id_sel][1];

    cs_real_2_t bi_bterms = {0.,0.};

    /* Pressure term */
    cs_real_t p_rho_id1 = pressure[c_id1] / rho[c_id1];
    cs_real_t p_rho_id2 = pressure[c_id2] / rho[c_id2];

    _i_conv_flux(p_rho_id1,
                 p_rho_id2,
                 i_mass_flux[f_id_sel],
                 bi_bterms);

    /* (The cell is counted only once in parallel by checking that
       the c_id is not in the halo) */
    /* Face normal well oriented (check bi_face_cells array) */
    if (bi_face_cells[f_id_sel][0] >= 0) {
      if (c_id1 < n_cells) {
        if (i_mass_flux[f_id_sel] > 0) {
          out_pressure += bi_bterms[0];
          out_debit += i_mass_flux[f_id_sel] / rho[c_id1];
          out_m_debit += i_mass_flux[f_id_sel];
        }
        else {
          in_pressure += bi_bterms[0];
          in_debit += i_mass_flux[f_id_sel] / rho[c_id1];
          in_m_debit += i_mass_flux[f_id_sel];
        }
      }
    }
    /* Face normal direction reversed */
    else {
      if (c_id2 < n_cells) {
        if (i_mass_flux[f_id_sel] > 0) {
          in_pressure -= bi_bterms[1];
          in_debit -= i_mass_flux[f_id_sel] / rho[c_id2];
          in_m_debit -= i_mass_flux[f_id_sel];
        }
        else {
          out_pressure -= bi_bterms[1];
          out_debit -= i_mass_flux[f_id_sel] / rho[c_id2];
          out_m_debit -= i_mass_flux[f_id_sel];
        }
      }
    }

    /* Kinematic term */
    bi_bterms[0] = 0.;
    bi_bterms[1] = 0.;

    cs_real_t u2_id1 = 0.5 * cs_math_3_square_norm(velocity[c_id1]);
    cs_real_t u2_id2 = 0.5 * cs_math_3_square_norm(velocity[c_id2]);

    _i_conv_flux(u2_id1,
                 u2_id2,
                 i_mass_flux[f_id_sel],
                 bi_bterms);

    /* (The cell is counted only once in parallel by checking that
       the c_id is not in the halo) */
    /* Face normal well oriented (check bi_face_cells array) */
    if (bi_face_cells[f_id_sel][0] >= 0) {
      if (c_id1 < n_cells) {
        if (i_mass_flux[f_id_sel] > 0) {
          out_u2 += bi_bterms[0];
        }
        else {
          in_u2 += bi_bterms[0];
        }
      }
    }
    /* Face normal direction reversed */
    else {
      if (c_id2 < n_cells) {
        if (i_mass_flux[f_id_sel] > 0) {
          in_u2 -= bi_bterms[1];
        }
        else {
          out_u2 -= bi_bterms[1];
        }
      }
    }

    /* Gravity term */
    bi_bterms[0] = 0.;
    bi_bterms[1] = 0.;

    cs_real_t gx_id1 = - cs_math_3_dot_product(gravity, i_face_cog[f_id_sel]);
    cs_real_t gx_id2 = - cs_math_3_dot_product(gravity, i_face_cog[f_id_sel]);

    _i_conv_flux(gx_id1,
                 gx_id2,
                 i_mass_flux[f_id_sel],
                 bi_bterms);

    /* (The cell is counted only once in parallel by checking that
       the c_id is not in the halo) */
    /* Face normal well oriented (check bi_face_cells array) */
    if (bi_face_cells[f_id_sel][0] >= 0) {
      if (c_id1 < n_cells) {
        if (i_mass_flux[f_id_sel] > 0) {
          out_rhogx += bi_bterms[0];
        }
        else {
          in_rhogx += bi_bterms[0];
        }
      }
    }
    /* Face normal direction reversed */
    else {
      if (c_id2 < n_cells) {
        if (i_mass_flux[f_id_sel] > 0) {
          in_rhogx -= bi_bterms[1];
        }
        else {
          out_rhogx -= bi_bterms[1];
        }
      }
    }

  }

  /* Free memory */

  CS_FREE(cells_tag_ids);
  CS_FREE(bi_face_cells);
  CS_FREE(i_face_sel_ids);
  CS_FREE(bb_face_sel_ids);
  CS_FREE(bi_face_sel_ids);

  /* Sum of values on all ranks (parallel calculations) */

  balance[CS_BALANCE_P_IN] = in_pressure;
  balance[CS_BALANCE_P_OUT] = out_pressure;
  balance[CS_BALANCE_P_U2_IN] = in_u2;
  balance[CS_BALANCE_P_U2_OUT] = out_u2;
  balance[CS_BALANCE_P_RHOGX_IN] = in_rhogx;
  balance[CS_BALANCE_P_RHOGX_OUT] = out_rhogx;
  balance[CS_BALANCE_P_U_IN] = in_debit;
  balance[CS_BALANCE_P_U_OUT] = out_debit;
  balance[CS_BALANCE_P_RHOU_IN] = in_m_debit;
  balance[CS_BALANCE_P_RHOU_OUT] = out_m_debit;

  cs_parall_sum(CS_BALANCE_P_N_TERMS, CS_REAL_TYPE, balance);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computes one term of the head loss balance (pressure drop) on a
 * volumic zone defined by the criterion also given as argument.
 * The different contributions are printed in the run_solver.log.
 *
 * \param[in]     selection_crit      zone selection criterion
 */
/*----------------------------------------------------------------------------*/

void
cs_pressure_drop_by_zone(const char * selection_crit)
{
  cs_real_t balance[CS_BALANCE_P_N_TERMS];

  const cs_mesh_t *m = cs_glob_mesh;
  const int nt_cur = cs_glob_time_step->nt_cur;

  /* Select cells */

  cs_lnum_t n_cells_sel = 0;
  cs_lnum_t *cells_sel_ids = nullptr;

  CS_MALLOC(cells_sel_ids, m->n_cells, cs_lnum_t);
  cs_selector_get_cell_list(selection_crit, &n_cells_sel, cells_sel_ids);

  /* Compute pressure drop terms */

  cs_pressure_drop_by_zone_compute(n_cells_sel,
                                   cells_sel_ids,
                                   balance);

  CS_FREE(cells_sel_ids);

  /* Log results at time step n */

  bft_printf(_("   ** Pressure drop by zone (at iteration %d)\n"
               "   ------------------------\n"
               "------------------------------------------------------------\n"
               "   Zone selection criteria: \"%s\"\n"
               "------------------------------------------------------------\n"
               "  |                 |\n"
               "  | p u . dS        | p u . dS\n"
               "  |   -    -        |   -    -\n"
               "  |                 |\n"
               "  | inlet           | outlet\n"
               "  %12.4e      %12.4e\n"
               "------------------------------------------------------------\n"
               "  |                 |\n"
               "  | u^2/2 rho u . dS| u^2/2 rho u . dS\n"
               "  | -         -    -| -         -    -\n"
               "  |                 |\n"
               "  | inlet           | outlet\n"
               "  %12.4e      %12.4e\n"
               "------------------------------------------------------------\n"
               "  |                 |\n"
               "  |-rho(g . x)u . dS|-rho(g . x)u . dS\n"
               "  |     -   - -    -|     -   - -    -\n"
               "  |                 |\n"
               "  | inlet           | outlet\n"
               "  %12.4e      %12.4e\n"
               "------------------------------------------------------------\n"
               "  |                 |\n"
               "  | u . dS          | u . dS\n"
               "  | -    -          | -    -\n"
               "  |                 |\n"
               "  | inlet           | outlet\n"
               "  %12.4e      %12.4e\n"
               "------------------------------------------------------------\n"
               "  |                 |\n"
               "  | rho u . dS      | rho u . dS\n"
               "  |     -    -      |     -    -\n"
               "  |                 |\n"
               "  | inlet           | outlet\n"
               "  %12.4e      %12.4e\n"
               "------------------------------------------------------------\n\n"),
             nt_cur, selection_crit,
             balance[CS_BALANCE_P_IN], balance[CS_BALANCE_P_OUT],
             balance[CS_BALANCE_P_U2_IN], balance[CS_BALANCE_P_U2_OUT],
             balance[CS_BALANCE_P_RHOGX_IN], balance[CS_BALANCE_P_RHOGX_OUT],
             balance[CS_BALANCE_P_U_IN], balance[CS_BALANCE_P_U_OUT],
             balance[CS_BALANCE_P_RHOU_IN], balance[CS_BALANCE_P_RHOU_OUT]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the surface balance of a given scalar.
 *
 * For interior faces, the flux is counted negatively relative to the given
 * normal (as neighboring interior faces may have differently-aligned normals).
 *
 * For boundary faces, the flux is counted negatively in the outwards-facing
 * direction.
 *
 * \param[in]     selection_crit      zone selection criterion
 * \param[in]     scalar_name         scalar name
 * \param[in]     normal              outwards normal direction
 */
/*----------------------------------------------------------------------------*/

void
cs_surface_balance(const char       *selection_crit,
                   const char       *scalar_name,
                   const cs_real_t   normal[3])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;

  const int nt_cur = cs_glob_time_step->nt_cur;

  /* Faces selection */

  cs_lnum_t n_b_faces_sel = 0;
  cs_lnum_t *b_face_sel_ids = nullptr;
  cs_lnum_t n_i_faces_sel = 0;
  cs_lnum_t *i_face_sel_ids = nullptr;

  CS_MALLOC(i_face_sel_ids, m->n_i_faces, cs_lnum_t);
  CS_MALLOC(b_face_sel_ids, m->n_b_faces, cs_lnum_t);

  cs_selector_get_i_face_list(selection_crit, &n_i_faces_sel, i_face_sel_ids);
  cs_selector_get_b_face_list(selection_crit, &n_b_faces_sel, b_face_sel_ids);

  /* Balance on selected faces */

  cs_real_t  balance[CS_BALANCE_N_TERMS];

  cs_flux_through_surface(scalar_name,
                          normal,
                          n_b_faces_sel,
                          n_i_faces_sel,
                          b_face_sel_ids,
                          i_face_sel_ids,
                          balance,
                          nullptr,   /* flux_b_faces */
                          nullptr);  /* flux_i_faces */

  /* Recount selected interior faces (parallel test) */

  cs_gnum_t n_sel[2] = {(cs_gnum_t)n_b_faces_sel, 0};

  for (cs_lnum_t i = 0; i < n_i_faces_sel; i++) {
    cs_lnum_t f_id = i_face_sel_ids[i];
    if (i_face_cells[f_id][0] < n_cells)
      n_sel[1] += 1;
  }

  cs_parall_sum(2, CS_GNUM_TYPE, n_sel);

  /* Free memory */

  CS_FREE(i_face_sel_ids);
  CS_FREE(b_face_sel_ids);

  /* Compute some sums */

  cs_real_t flux_b_faces
    = balance[CS_BALANCE_BOUNDARY_IN] + balance[CS_BALANCE_BOUNDARY_OUT]
    + balance[CS_BALANCE_BOUNDARY_SYM] + balance[CS_BALANCE_BOUNDARY_WALL]
    + balance[CS_BALANCE_BOUNDARY_COUPLED_E]
    + balance[CS_BALANCE_BOUNDARY_OTHER];

  cs_real_t flux_i_faces
    = balance[CS_BALANCE_INTERIOR_IN] + balance[CS_BALANCE_INTERIOR_OUT];

  /* Log balance */

  bft_printf
    (_("\n   ** Surface balance at iteration %d\n"
       "     ------------------\n"
       "------------------------------------------------------------\n"
       "   Scalar: %s\n"
       "   Zone selection criteria: \"%s\"\n"
       "   Outgoing normal: [%.2e, %.2e, %.2e] \n"
       "------------------------------------------------------------\n"
       "   Interior faces selected: %llu of %llu \n"
       "   Boundary faces selected: %llu of %llu \n"
       "------------------------------------------------------------\n"
       "    Boundary faces:        %12.4e \n"
       "    Int. Coupling faces:   %12.4e \n"
       "    Interior faces:        \n"
       "      In:                  %12.4e \n"
       "      Out:                 %12.4e \n"
       "      Balance:             %12.4e \n"
       "------------------------------------------------------------\n"),
     nt_cur, scalar_name, selection_crit,
     normal[0], normal[1], normal[2],
     (unsigned long long)n_sel[1], (unsigned long long)(m->n_g_i_faces),
     (unsigned long long)n_sel[0], (unsigned long long)(m->n_g_b_faces),
     flux_b_faces, balance[CS_BALANCE_BOUNDARY_COUPLED_E],
     balance[CS_BALANCE_INTERIOR_IN], balance[CS_BALANCE_INTERIOR_OUT],
     flux_i_faces);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the face by face surface flux of a given scalar, through a
 *        surface area defined by the given face ids.
 *
 * For interior faces, the flux is counted negatively relative to the given
 * normal (as neighboring interior faces may have differently-aligned normals).
 *
 * For boundary faces, the flux is counted negatively in the outwards-facing
 * direction.
 *
 * \param[in]   scalar_name       scalar name
 * \param[in]   normal            outwards normal direction
 * \param[in]   n_b_faces_sel     number of selected boundary faces
 * \param[in]   n_i_faces_sel     number of selected internal faces
 * \param[in]   b_face_sel_ids    ids of selected boundary faces
 * \param[in]   i_face_sel_ids    ids of selected internal faces
 * \param[out]  balance           optional array of computed balance terms
 *                                (see \ref cs_balance_term_t), of
 *                                size CS_BALANCE_N_TERMS, or nullptr
 * \param[out]  flux_b_faces      optional surface flux through selected
 *                                boundary faces, or nullptr
 * \param[out]  flux_i_faces      optional surface flux through selected
 *                                interior faces, or nullptr
 */
/*----------------------------------------------------------------------------*/

void
cs_flux_through_surface(const char         *scalar_name,
                        const cs_real_t     normal[3],
                        cs_lnum_t           n_b_faces_sel,
                        cs_lnum_t           n_i_faces_sel,
                        const cs_lnum_t     b_face_sel_ids[],
                        const cs_lnum_t     i_face_sel_ids[],
                        cs_real_t          *balance,
                        cs_real_t          *flux_b_faces,
                        cs_real_t          *flux_i_faces)
{
  const cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_nreal_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;

  /* Parallel or device dispatch */
  cs_dispatch_context ctx;
  ctx.set_use_gpu(false);  /* not yet ported to GPU */

  const int *bc_type = cs_glob_bc_type;

  cs_field_t *f = cs_field_by_name_try(scalar_name);
  if (f == nullptr)
    return;

  const int field_id = cs_field_id_by_name(scalar_name);

  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);

  /* initialize output */

  cs_real_t  _balance[CS_BALANCE_N_TERMS];
  for (int i = 0; i < CS_BALANCE_N_TERMS; i++)
    _balance[i] = 0;

 /* Physical properties
    ------------------- */

  /* Specific heat (CP) */
  const cs_real_t *xcpp = nullptr;
  cs_real_t *_xcpp = nullptr;
  cs_real_t cp0 = 1;
  if (f->get_key_int("is_temperature")) {
    const int icp = cs_field_id_by_name("specific_heat");
    if (icp != -1)
      xcpp = CS_F_(cp)->val;
    else
      cp0 = cs_glob_fluid_properties->cp0;
  }
  if (xcpp == nullptr) {
    CS_MALLOC_HD(_xcpp, n_cells_ext, cs_real_t, cs_alloc_mode);
    cs_arrays_set_value<cs_real_t, 1>(n_cells_ext, cp0, _xcpp);
    xcpp = _xcpp;
  }

  /* Convective mass fluxes for inner and boundary faces */
  int iflmas = f->get_key_int("inner_mass_flux_id");
  const cs_real_t *i_mass_flux = cs_field_by_id(iflmas)->val;

  int iflmab = f->get_key_int("boundary_mass_flux_id");
  const cs_real_t *b_mass_flux = cs_field_by_id(iflmab)->val;

  if (field_id != -1) {
    f = cs_field_by_id(field_id);

    /* Update bc_coeffs val_f and flux_diff if nullptr.
       (It may occur when scalars are not computed,
       for example when nt_max = 0) */

    if (f->bc_coeffs->val_f == nullptr || f->bc_coeffs->flux_diff == nullptr)
      cs_boundary_conditions_update_bc_coeff_face_values
        (ctx, f, eqp,
         true, true,
         0, nullptr, // hyd_p_flag, f_ext
         nullptr, // visel
         nullptr, // weighb
         f->val);
  }

  /* Face viscosity */
  int imvisf = eqp->imvisf;
  cs_real_t *i_visc;
  cs_real_t *b_visc;
  CS_MALLOC_HD(i_visc, n_i_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(b_visc, n_b_faces, cs_real_t, cs_alloc_mode);

  cs_real_t *c_visc = nullptr;
  CS_MALLOC_HD(c_visc, n_cells_ext, cs_real_t, cs_alloc_mode);
  const int kivisl = f->get_key_int("diffusivity_id");
  if (kivisl != -1) {
    cs_array_copy(n_cells_ext, cs_field_by_id(kivisl)->val, c_visc);
  }
  else {
    const double visls0 = f->get_key_double("diffusivity_ref");
    cs_arrays_set_value<cs_real_t, 1>(n_cells_ext, visls0, c_visc);
  }

  /* Turbulent part */
  cs_real_t *c_visct = cs_field_by_name("turbulent_viscosity")->val;

  if (eqp->idifft == 1) {
    const int ksigmas = cs_field_key_id("turbulent_schmidt");
    const cs_real_t turb_schmidt = f->get_key_double(ksigmas);
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      c_visc[c_id] += xcpp[c_id] * c_visct[c_id]/turb_schmidt;
  }

  cs_face_viscosity(m, fvq, imvisf, c_visc, i_visc, b_visc);

  /* Choose gradient type */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  const int imrgra = eqp->imrgra;
  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  /* Limiters */

  if (field_id != -1) {
    f = cs_field_by_id(field_id);
  }

  /* Faces selection
     --------------- */

  cs_lnum_2_t *bi_face_cells = nullptr;

  if (n_i_faces_sel > 0) {

    CS_MALLOC(bi_face_cells, n_i_faces, cs_lnum_2_t);
    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
      bi_face_cells[f_id][0] = -999;
      bi_face_cells[f_id][1] = -999;
    }

    for (cs_lnum_t f_id = 0; f_id < n_i_faces_sel; f_id++) {
      cs_lnum_t f_id_sel = i_face_sel_ids[f_id];
      cs_lnum_t c_id1 = i_face_cells[f_id_sel][0];
      cs_lnum_t c_id2 = i_face_cells[f_id_sel][1];

      cs_real_t dot_pro = cs_math_3_dot_product(normal,
                                                i_face_u_normal[f_id_sel]);
      if (fabs(dot_pro) < 1.0e-14)//FIXME
        dot_pro = 0;
      if (dot_pro > 0.)
        bi_face_cells[f_id_sel][0] = c_id1;
      else if (dot_pro < 0.)
        bi_face_cells[f_id_sel][1] = c_id2;
    }

    if (flux_i_faces != nullptr) {
      for (cs_lnum_t f_id = 0; f_id < n_i_faces_sel; f_id++)
        flux_i_faces[f_id] = 0.;
    }

    if (flux_b_faces != nullptr) {
      for (cs_lnum_t f_id = 0; f_id < n_b_faces_sel; f_id++)
        flux_b_faces[f_id] = 0.;
    }

  }

  cs_real_t *c_weight = nullptr;
  if (f->type & CS_FIELD_VARIABLE && eqp->iwgrec == 1) {
    if (eqp->idiff > 0) {
      int key_id = cs_field_key_id("gradient_weighting_id");
      int diff_id = f->get_key_int(key_id);
      if (diff_id > -1) {
        cs_field_t *weight_f = cs_field_by_id(diff_id);
        c_weight = weight_f->val;
      }
    }
  }

  cs_array_2d<cs_real_t> i_flux(n_i_faces_sel, 2, cs_alloc_mode);
  cs_array<cs_real_t> b_flux(n_b_faces_sel, cs_alloc_mode);
  cs_array<int> icvfli(n_b_faces_sel, cs_alloc_mode);

  i_flux.zero();
  b_flux.zero();
  cs_arrays_set_value<int, 1>(n_b_faces_sel, 0, icvfli.data());

  cs_convection_diffusion_scalar_at_faces(f,
                                          *eqp,
                                          0, // icvflb
                                          n_i_faces_sel,
                                          n_b_faces_sel,
                                          i_face_sel_ids,
                                          b_face_sel_ids,
                                          f->val,
                                          icvfli.data(),
                                          f->bc_coeffs,
                                          i_mass_flux,
                                          b_mass_flux,
                                          i_visc,
                                          b_visc,
                                          c_weight,
                                          xcpp,
                                          (cs_real_2_t *)i_flux.data(),
                                          b_flux.data());

  /* Boundary faces contribution
     --------------------------- */

  for (cs_lnum_t f_id = 0; f_id < n_b_faces_sel; f_id++) {

    cs_lnum_t f_id_sel = (b_face_sel_ids != nullptr) ?
      b_face_sel_ids[f_id] : f_id;

    cs_real_t term_balance = b_flux(f_id);

    if (flux_b_faces != nullptr)
      flux_b_faces[f_id] = term_balance;

    if (bc_type[f_id_sel] == CS_INLET ||
        bc_type[f_id_sel] == CS_CONVECTIVE_INLET ||
        bc_type[f_id_sel] == CS_FREE_INLET ||
        bc_type[f_id_sel] == CS_ESICF ||
        bc_type[f_id_sel] == CS_EPHCF)
      _balance[CS_BALANCE_BOUNDARY_IN] -= term_balance;
    else if (bc_type[f_id_sel] == CS_OUTLET ||
             bc_type[f_id_sel] == CS_SSPCF ||
             bc_type[f_id_sel] == CS_SOPCF)
      _balance[CS_BALANCE_BOUNDARY_OUT] -= term_balance;
    else if (bc_type[f_id_sel] == CS_SYMMETRY)
      _balance[CS_BALANCE_BOUNDARY_SYM] -= term_balance;
    else if (bc_type[f_id_sel] == CS_SMOOTHWALL)
      _balance[CS_BALANCE_BOUNDARY_WALL_S] -= term_balance;
    else if (bc_type[f_id_sel] == CS_ROUGHWALL)
      _balance[CS_BALANCE_BOUNDARY_WALL_R] -= term_balance;
    else if (   bc_type[f_id_sel] == CS_COUPLED
             || bc_type[f_id_sel] == CS_COUPLED_FD)
      _balance[CS_BALANCE_BOUNDARY_COUPLED_E] -= term_balance;
    else
      _balance[CS_BALANCE_BOUNDARY_OTHER] -= term_balance;

  }

  /* Balance on selected interior faces */

  for (cs_lnum_t f_id = 0; f_id < n_i_faces_sel; f_id++) {

    cs_lnum_t f_id_sel = i_face_sel_ids[f_id];
    /* Associated boundary-internal cells */
    cs_lnum_t c_id1 = i_face_cells[f_id_sel][0];
    cs_lnum_t c_id2 = i_face_cells[f_id_sel][1];

    /* (The cell is counted only once in parallel by checking that
       the c_id is not in the halo) */
    /* Face normal well oriented (check i_face_cells array) */
    if (bi_face_cells[f_id_sel][0] >= 0) {
      if (c_id1 < n_cells) {
        if (flux_i_faces != nullptr)
          flux_i_faces[f_id] -= i_flux(f_id, 0);
        if (i_mass_flux[f_id_sel] > 0)
          _balance[CS_BALANCE_INTERIOR_IN] -= i_flux(f_id, 0);
        else
          _balance[CS_BALANCE_INTERIOR_OUT] -= i_flux(f_id, 0);
      }
    }
    /* Face normal direction reversed */
    else if (bi_face_cells[f_id_sel][1] >= 0) {
      if (c_id2 < n_cells) {
        if (flux_i_faces != nullptr)
          flux_i_faces[f_id] += i_flux(f_id, 1);
        if (i_mass_flux[f_id_sel] > 0)
          _balance[CS_BALANCE_INTERIOR_IN] += i_flux(f_id, 1);
        else
          _balance[CS_BALANCE_INTERIOR_OUT] += i_flux(f_id, 1);
      }
    }
  }

  if (balance != nullptr) {
    _balance[CS_BALANCE_BOUNDARY_WALL] =   _balance[CS_BALANCE_BOUNDARY_WALL_S]
                                         + _balance[CS_BALANCE_BOUNDARY_WALL_R];
    _balance[CS_BALANCE_BOUNDARY_COUPLED]
      =   _balance[CS_BALANCE_BOUNDARY_COUPLED_E]
        + _balance[CS_BALANCE_BOUNDARY_COUPLED_I];

    for (int i = 0; i < CS_BALANCE_N_TERMS; i++)
      balance[i] = _balance[i];
    cs_parall_sum(CS_BALANCE_N_TERMS, CS_REAL_TYPE, balance);
  }

  /* Free memory */

  CS_FREE(bi_face_cells);

  CS_FREE(_xcpp);
  CS_FREE(c_visc);
  CS_FREE(i_visc);
  CS_FREE(b_visc);
}

/*----------------------------------------------------------------------------*/
