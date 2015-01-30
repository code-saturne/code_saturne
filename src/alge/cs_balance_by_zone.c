/*============================================================================
 * Scalar balance on zones.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_writer.h"

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

#include "cs_face_viscosity.h"
#include "cs_gradient_perio.h"
#include "cs_physical_constants.h"
#include "cs_thermal_model.h"
#include "cs_convection_diffusion.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_balance_by_zone.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file  cs_balance_by_zone.c

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computes the different terms of the balance of a scalar which name is
 * given as argument, on a volumic zone defined by the criterium also given as
 * argument. The different contributions to the balance are printed in the
 * listing.
 *
 * \param[in]     bc_type             boundary condition type
 * \param[in]     selection_crit      zone selection criterium
 * \param[in]     scalar_name         scalar name
 */
/*----------------------------------------------------------------------------*/

void
cs_balance_by_zone(const int  bc_type[],
                   const char *selection_crit,
                   const char *scalar_name)
{
  int nt_cur = cs_glob_time_step->nt_cur;
  int idtvar = cs_glob_time_step_options->idtvar;

  const cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_halo_t  *halo = m->halo;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict i_face_surf = fvq->i_face_surf;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)fvq->i_face_normal;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;
  const cs_real_3_t *restrict dijpf
    = (const cs_real_3_t *restrict)fvq->dijpf;
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;

  /* all boundary convective fluxes are upwind */
  int icvflb = 0;
  int *icvfli = NULL;

  /* coupling is disabled */
  int ifaccp = 0;

  /* Get physical fields */
  const cs_real_t *dt = CS_F_(dt)->val;
  const cs_real_t *rho = CS_F_(rho)->val;
  const cs_field_t *f = cs_field_by_name_try(scalar_name);
  const int field_id = cs_field_id_by_name(scalar_name);

  /* If the requested scalar field is not computed, return */
  if (field_id == -1) {
    bft_printf("Scalar field does not exist. Balance will not be computed.\n");
    return;
  }

  /* Temperature indicator.
     Will multiply by CP in order to have energy. */
  bool itemperature = false;
  const int scal_id = cs_field_get_key_int(f, cs_field_key_id("scalar_id"));
  if (scal_id == cs_glob_thermal_model->iscalt) {
    if (cs_glob_thermal_model->itherm == 1)
      itemperature = true;
  }

  /* Specific heat (CP) */
  cs_real_t *cpro_cp = NULL;
  const int icp = cs_field_id_by_name("specific_heat");
  if (itemperature) {
    if (icp != -1)
      cpro_cp = CS_F_(cp)->val;
    else {
      const double cp0 = cs_glob_fluid_properties->cp0;
      BFT_MALLOC(cpro_cp, n_cells, cs_real_t);
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cpro_cp[c_id] = cp0;
      }
    }
  }
  else {
    BFT_MALLOC(cpro_cp, n_cells, cs_real_t);
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      cpro_cp[c_id] = 1.;
    }
  }

  /* Zone cells selection variables*/
  cs_lnum_t n_cells_sel = 0;
  cs_lnum_t *cells_sel_list = NULL;
  cs_lnum_t n_i_faces_sel = 0;
  cs_lnum_t *i_face_sel_list = NULL;
  cs_lnum_t n_bb_faces_sel = 0;
  cs_lnum_t *bb_face_sel_list = NULL;
  cs_lnum_t n_bi_faces_sel = 0;
  cs_lnum_t *bi_face_sel_list = NULL;
  cs_lnum_2_t *bi_face_cells = NULL;
  cs_lnum_t *cells_tag_list = NULL;

  /*--------------------------------------------------------------------------
   * This example computes the balance relative to a given scalar
   * on a selected zone of the mesh.
   * We assume that we want to compute balances (convective and diffusive)
   * at the boundaries of the calculation domain represented below
   * (with different boundary types).
   *
   * The scalar and the zone are selected at the top of the routine
   * by the user.
   * In the case of the temperature, the energy balance in Joules will be
   * computed by multiplying by the specific heat.
   *--------------------------------------------------------------------------*/

  /* 1. Initialization
     =================

    --> List of balance contributions
        -----------------------------

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
    ndef_balance  : contribution from undefined faces
    tot_balance   : total balance */

  double vol_balance = 0.;
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
  double ndef_balance = 0.;
  double tot_balance = 0.;

  /* Boundary condition coefficient for h */
  const cs_real_t *a_F = f->bc_coeffs->a;
  const cs_real_t *b_F = f->bc_coeffs->b;
  const cs_real_t *af_F = f->bc_coeffs->af;
  const cs_real_t *bf_F = f->bc_coeffs->bf;

  /* Convective mass fluxes for inner and boundary faces */
  int iflmas = cs_field_get_key_int(f, cs_field_key_id("inner_mass_flux_id"));
  const cs_real_t *i_mass_flux = cs_field_by_id(iflmas)->val;

  int iflmab = cs_field_get_key_int(f, cs_field_key_id("boundary_mass_flux_id"));
  const cs_real_t *b_mass_flux = cs_field_by_id(iflmab)->val;

  /* Allocate temporary array */
  cs_real_t *f_reconstructed;
  BFT_MALLOC(f_reconstructed, n_b_faces, cs_real_t);

  /* Reconstructed value */
  cs_real_3_t *grad;
  BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);

  int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_var_cal_opt_t var_cal_opt;

  /* Get the calculation option from the field */
  cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);

  cs_halo_type_t halo_type;
  cs_gradient_type_t gradient_type;

  cs_gradient_type_by_imrgra(var_cal_opt.imrgra,
                             &gradient_type,
                             &halo_type);

  cs_field_gradient_scalar(f,
                           true, /* use_previous_t */
                           gradient_type,
                           halo_type,
                           1, /* inc */
                           true, /* _recompute_cocg */
                           grad);

  if (false) {
    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      /* Associated boundary cell */
      cs_lnum_t c_id = b_face_cells[f_id];
      f_reconstructed[f_id] = f->val[c_id]
                               + grad[c_id][0]*diipb[f_id][0]
                               + grad[c_id][1]*diipb[f_id][1]
                               + grad[c_id][2]*diipb[f_id][2];
    }

  /* Non-reconstructed value */
  } else {
    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      /* Associated boundary cell */
      cs_lnum_t c_id = b_face_cells[f_id];
      f_reconstructed[f_id] = f->val[c_id];
    }
  }

  int inc = 1;

  /* Gradient for slope test */
  cs_real_3_t *grdpa = NULL;
  if (var_cal_opt.blencv > 0 && var_cal_opt.isstpc == 0) {
    BFT_MALLOC(grdpa, n_cells_ext, cs_real_3_t);

    for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
      grdpa[c_id][0] = 0.;
      grdpa[c_id][1] = 0.;
      grdpa[c_id][2] = 0.;
    }

    if (var_cal_opt.iconv > 0)
      cs_slope_test_gradient(field_id,
                             inc,
                             halo_type,
                             grad,
                             grdpa,
                             f->val,
                             a_F,
                             b_F,
                             i_mass_flux);
  }

  /* Face viscosity */
  int imvisf = cs_glob_space_disc->imvisf;
  cs_real_t *i_visc;
  cs_real_t *b_visc;
  BFT_MALLOC(i_visc, n_i_faces, cs_real_t);
  BFT_MALLOC(b_visc, n_b_faces, cs_real_t);

  cs_real_t *c_visc = NULL;
  const int kivisl =
    cs_field_get_key_int(f, cs_field_key_id("scalar_diffusivity_id"));
  if (kivisl != -1)
    c_visc = cs_field_by_id(kivisl)->val;
  else {
    const double visls0 =
      cs_field_get_key_double(f, cs_field_key_id("scalar_diffusivity_ref"));
    BFT_MALLOC(c_visc, n_cells_ext, cs_real_t);
    for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
      c_visc[c_id] = visls0;
    }
  }

  cs_face_viscosity(m, fvq, imvisf, c_visc, i_visc, b_visc);

  /* =========================================================================
     ---> Get user-selected zone
     =========================================================================*/

  /* Initialise arrays */

  /* Internal faces of the selected zone */
  BFT_MALLOC(i_face_sel_list, n_i_faces, cs_lnum_t);
  /* Boundary faces of the selected zone,
     which are internal faces of the global mesh.
     Faces -> cells connectivity */
  BFT_MALLOC(bi_face_sel_list, n_i_faces, cs_lnum_t);
  BFT_MALLOC(bi_face_cells, n_i_faces, cs_lnum_2_t);
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    i_face_sel_list[f_id] = -1;
    bi_face_sel_list[f_id] = -1;
    bi_face_cells[f_id][0] = -999;
    bi_face_cells[f_id][1] = -999;
  }

  /* Boundary faces of the selected zone,
     which are also boundary faces of the global mesh */
  BFT_MALLOC(bb_face_sel_list, n_b_faces, cs_lnum_t);
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    bb_face_sel_list[f_id] = -1;
  }

  /* Select cells */
  BFT_MALLOC(cells_sel_list, n_cells, cs_lnum_t);
  cs_selector_get_cell_list(selection_crit, &n_cells_sel, cells_sel_list);

  /* Synchronization for parallelism */
  BFT_MALLOC(cells_tag_list, n_cells_ext, cs_lnum_t);
  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++){
    cells_tag_list[c_id] = 0;
  }
  for (cs_lnum_t c_id = 0; c_id < n_cells_sel; c_id++){
    cs_lnum_t c_id_sel = cells_sel_list[c_id];
    cells_tag_list[c_id_sel] = 1;
  }
  if (halo != NULL) {
    cs_halo_sync_num(halo, halo_type, cells_tag_list);
  }

  /* Classify mesh faces with respect to the selected zone */

  /* Check boundary faces:
     if they are in the selected zone, they are boundary as well */
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    cs_lnum_t c_id = b_face_cells[f_id];

    if (cells_tag_list[c_id] == 1) {
      n_bb_faces_sel++;
      bb_face_sel_list[n_bb_faces_sel-1] = f_id;
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

    if (cells_tag_list[c_id1] == 1)
      indic1 = true;
    if (cells_tag_list[c_id2] == 1)
      indic2 = true;

    if (indic1 && indic2) {
      n_i_faces_sel++;
      i_face_sel_list[n_i_faces_sel-1] = f_id;
    }
    else if (indic1 || indic2) {
      n_bi_faces_sel++;
      bi_face_sel_list[n_bi_faces_sel-1] = f_id;
      /* Build the faces -> cells connectivity as done in
         i_face_cells */
      if (indic1)
        bi_face_cells[f_id][0] = c_id1;
      else
        bi_face_cells[f_id][1] = c_id2;
    }

  }

  /* =========================================================================
     ---> Balance computation
     =========================================================================*/

  /* 2. Compute the balance at time step n
    ======================================

    --> Balance on interior volumes
        --------------------------- */

  for (cs_lnum_t c_id = 0; c_id < n_cells_sel; c_id++) {

    cs_lnum_t c_id_sel = cells_sel_list[c_id];

    vol_balance += cell_vol[c_id_sel] * rho[c_id_sel]
                 * cpro_cp[c_id_sel]
                 * (f->val_pre[c_id_sel] - f->val[c_id_sel]);
  }

  /*
    --> Balance on all faces (interior and boundary), for div(rho u)
        ------------------------------------------------------------
   */

  /* Interior faces */
  for (cs_lnum_t f_id = 0; f_id < n_i_faces_sel; f_id++) {

    cs_lnum_t f_id_sel = i_face_sel_list[f_id];
    /* Associated internal cells */
    cs_lnum_t c_id1 = i_face_cells[f_id_sel][0];
    cs_lnum_t c_id2 = i_face_cells[f_id_sel][1];

    /* Contribution to flux from the two cells of the current face
      (The cell is counted only once in parallel by checking that
       the c_id is not in the halo) */

    if (c_id1 < n_cells)
      div_balance += i_mass_flux[f_id_sel] * dt[c_id1] * f->val[c_id1]
                   * cpro_cp[c_id1];

    if (c_id2 < n_cells)
      div_balance -= i_mass_flux[f_id_sel] * dt[c_id2] * f->val[c_id2]
                   * cpro_cp[c_id2];

  }

  /* Boundary faces which are internal in the total mesh */
  for (cs_lnum_t f_id = 0; f_id < n_bi_faces_sel; f_id++) {

    cs_lnum_t f_id_sel = bi_face_sel_list[f_id];
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
                     * cpro_cp[c_id1];
    }

    else {

      if (c_id2 < n_cells)
        div_balance -= i_mass_flux[f_id_sel] * dt[c_id2] * f->val[c_id2]
                     * cpro_cp[c_id2];
    }

  }

  /* Boundary faces which are also boundary in the total mesh */
  for (cs_lnum_t f_id = 0; f_id < n_bb_faces_sel; f_id++) {

    cs_lnum_t f_id_sel = bb_face_sel_list[f_id];
    /* Associated boundary cell */
    cs_lnum_t c_id = b_face_cells[f_id_sel];

    /* Contribution to flux from the current face */
      div_balance += b_mass_flux[f_id_sel] * dt[c_id] * f->val[c_id]
                   * cpro_cp[c_id];

  }

  // TODO mass source terms and mass accumulation term
  // In case of a mass source term, add contribution from Gamma*Tn+1

  int iconvp = var_cal_opt.iconv;
  int idiffp = var_cal_opt.idiff;
  int ircflp = var_cal_opt.ircflu;
  double relaxp = var_cal_opt.relaxv;
  /*
    --> Balance on boundary faces
        -------------------------

    We handle different types of boundary faces separately to better
    analyze the information, but this is not mandatory. */

  cs_real_t pir, pip, pjp, pipr, pjpr;
  cs_real_t pifri, pjfri, pifrj, pjfrj;

  if (icvflb == 0) {
    /* ====================
       ---> Upwind
       ====================*/

    /* Steady */
    if (idtvar < 0) {

      for (cs_lnum_t f_id = 0; f_id < n_bb_faces_sel; f_id++) {

        cs_lnum_t f_id_sel = bb_face_sel_list[f_id];
        /* Associated boundary cell */
        cs_lnum_t c_id = b_face_cells[f_id_sel];

        cs_b_cd_steady(ircflp,
                       relaxp,
                       diipb[f_id_sel],
                       grad[c_id],
                       f->val[c_id],
                       f->val_pre[c_id],
                       &pir,
                       &pipr);

        cs_real_t term_balance = 0.;

        cs_b_upwind_flux_cons(iconvp,
                              inc,
                              ifaccp,
                              bc_type[f_id_sel],
                              pir,
                              pipr,
                              a_F[f_id_sel],
                              b_F[f_id_sel],
                              b_mass_flux[f_id_sel],
                              cpro_cp[c_id],
                              &term_balance);

        cs_b_diff_flux(idiffp,
                       inc,
                       pipr,
                       af_F[f_id_sel],
                       bf_F[f_id_sel],
                       b_visc[f_id_sel],
                       &term_balance);

        if (bc_type[f_id_sel] == CS_INLET ||
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
        else if (bc_type[f_id_sel] == CS_COUPLED)
          cpl_balance -= term_balance*dt[c_id];
        else
          ndef_balance -= term_balance*dt[c_id];
      }
    /* Unsteady */
    } else {

      for (cs_lnum_t f_id = 0; f_id < n_bb_faces_sel; f_id++) {

        cs_lnum_t f_id_sel = bb_face_sel_list[f_id];
        /* Associated boundary cell */
        cs_lnum_t c_id = b_face_cells[f_id_sel];

        cs_b_cd_unsteady(ircflp,
                         diipb[f_id_sel],
                         grad[c_id],
                         f->val[c_id],
                         &pir,
                         &pipr);

        cs_real_t term_balance = 0.;

        cs_b_upwind_flux_cons(iconvp,
                              inc,
                              ifaccp,
                              bc_type[f_id_sel],
                              pir,
                              pipr,
                              a_F[f_id_sel],
                              b_F[f_id_sel],
                              b_mass_flux[f_id_sel],
                              cpro_cp[c_id],
                              &term_balance);

        cs_b_diff_flux(idiffp,
                       inc,
                       pipr,
                       af_F[f_id_sel],
                       bf_F[f_id_sel],
                       b_visc[f_id_sel],
                       &term_balance);

        if (bc_type[f_id_sel] == CS_INLET ||
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
        else if (bc_type[f_id_sel] == CS_COUPLED)
          cpl_balance -= term_balance*dt[c_id];
        else
          ndef_balance -= term_balance*dt[c_id];
      }
    }

  /* Boundary convective flux is imposed at some faces
     (tagged in icvfli array) */
  } else { /* icvflb =1 */
    /* ====================
       ---> Imposed
       ====================*/

    const cs_real_t *ac_F = f->bc_coeffs->ac;
    const cs_real_t *bc_F = f->bc_coeffs->bc;

    /* Steady */
    if (idtvar < 0) {

      for (cs_lnum_t f_id = 0; f_id < n_bb_faces_sel; f_id++) {

        cs_lnum_t f_id_sel = bb_face_sel_list[f_id];
        /* Associated boundary cell */
        cs_lnum_t c_id = b_face_cells[f_id_sel];

        cs_b_cd_steady(ircflp,
                       relaxp,
                       diipb[f_id_sel],
                       grad[c_id],
                       f->val[c_id],
                       f->val_pre[c_id],
                       &pir,
                       &pipr);

        cs_real_t term_balance = 0.;

        cs_b_imposed_conv_flux_cons(iconvp,
                                    inc,
                                    ifaccp,
                                    bc_type[f_id_sel],
                                    icvfli[f_id_sel],
                                    pir,
                                    pipr,
                                    a_F[f_id_sel],
                                    b_F[f_id_sel],
                                    ac_F[f_id_sel],
                                    bc_F[f_id_sel],
                                    b_mass_flux[f_id_sel],
                                    cpro_cp[c_id],
                                    &term_balance);

        cs_b_diff_flux(idiffp,
                       inc,
                       pipr,
                       af_F[f_id_sel],
                       bf_F[f_id_sel],
                       b_visc[f_id_sel],
                       &term_balance);

        if (bc_type[f_id_sel] == CS_INLET ||
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
        else if (bc_type[f_id_sel] == CS_COUPLED)
          cpl_balance -= term_balance*dt[c_id];
        else
          ndef_balance -= term_balance*dt[c_id];
      }
    /* Unsteady */
    } else {

      for (cs_lnum_t f_id = 0; f_id < n_bb_faces_sel; f_id++) {

        cs_lnum_t f_id_sel = bb_face_sel_list[f_id];
        /* Associated boundary cell */
        cs_lnum_t c_id = b_face_cells[f_id_sel];

        cs_b_cd_unsteady(ircflp,
                         diipb[f_id_sel],
                         grad[c_id],
                         f->val[c_id],
                         &pir,
                         &pipr);

        cs_real_t term_balance = 0.;

        cs_b_imposed_conv_flux_cons(iconvp,
                                    inc,
                                    ifaccp,
                                    bc_type[f_id_sel],
                                    icvfli[f_id_sel],
                                    pir,
                                    pipr,
                                    a_F[f_id_sel],
                                    b_F[f_id_sel],
                                    ac_F[f_id_sel],
                                    bc_F[f_id_sel],
                                    b_mass_flux[f_id_sel],
                                    cpro_cp[c_id],
                                    &term_balance);

        cs_b_diff_flux(idiffp,
                       inc,
                       pipr,
                       af_F[f_id_sel],
                       bf_F[f_id_sel],
                       b_visc[f_id_sel],
                       &term_balance);

        if (bc_type[f_id_sel] == CS_INLET ||
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
        else if (bc_type[f_id_sel] == CS_COUPLED)
          cpl_balance -= term_balance*dt[c_id];
        else
          ndef_balance -= term_balance*dt[c_id];
      }
    }
  }

  /*
    --> Balance on boundary faces of the selected zone
        that are internal of the total mesh
        ------------------------------------------------------------
   */

  int isstpp = var_cal_opt.isstpc;
  int ischcp = var_cal_opt.ischcv;
  double blencp = var_cal_opt.blencv;
  int iupwin = (blencp > 0.) ? 0 : 1;

  /* Upwind indicator (useless here) */
  bool indic = false;

  if (iupwin == 1) {
    /* ====================
       ---> Upwind
       ====================*/

    /* Steady */
    if (idtvar < 0) {

      for (cs_lnum_t f_id = 0; f_id < n_bi_faces_sel; f_id++) {

        cs_lnum_t f_id_sel = bi_face_sel_list[f_id];
        /* Associated boundary-internal cells */
        cs_lnum_t c_id1 = i_face_cells[f_id_sel][0];
        cs_lnum_t c_id2 = i_face_cells[f_id_sel][1];

        cs_real_2_t bi_bterms = {0.,0.};

        cs_i_cd_steady_upwind(ircflp,
                              relaxp,
                              weight[f_id],
                              cell_cen[c_id1],
                              cell_cen[c_id2],
                              i_face_cog[f_id_sel],
                              dijpf[f_id_sel],
                              grad[c_id1],
                              grad[c_id2],
                              f->val[c_id1],
                              f->val[c_id2],
                              f->val_pre[c_id1],
                              f->val_pre[c_id2],
                              &pifri,
                              &pifrj,
                              &pjfri,
                              &pjfrj,
                              &pip,
                              &pjp,
                              &pipr,
                              &pjpr);

        cs_i_conv_flux_cons(iconvp,
                            pifri,
                            pifrj,
                            pjfri,
                            pjfrj,
                            i_mass_flux[f_id_sel],
                            cpro_cp[c_id1],
                            cpro_cp[c_id2],
                            bi_bterms);

        cs_i_diff_flux(idiffp,
                       pip,
                       pjp,
                       pipr,
                       pjpr,
                       i_visc[f_id_sel],
                       bi_bterms);

        /* (The cell is counted only once in parallel by checking that
           the c_id is not in the halo) */
        /* Face normal well oriented (check bi_face_cells array) */
        if (bi_face_cells[f_id_sel][0] >= 0) {
          if (c_id1 < n_cells) {
            if (i_mass_flux[f_id_sel] > 0)
              bi_o_balance -= bi_bterms[0]*dt[c_id1];
            else
              bi_i_balance -= bi_bterms[0]*dt[c_id1];
          }
        }
        /* Face normal direction reversed */
        else {
          if (c_id2 < n_cells) {
            if (i_mass_flux[f_id_sel] > 0)
              bi_i_balance += bi_bterms[1]*dt[c_id2];
            else
              bi_o_balance += bi_bterms[1]*dt[c_id2];
          }
        }
      }
    /* Unsteady */
    } else {

      for (cs_lnum_t f_id = 0; f_id < n_bi_faces_sel; f_id++) {

        cs_lnum_t f_id_sel = bi_face_sel_list[f_id];
        /* Associated boundary-internal cells */
        cs_lnum_t c_id1 = i_face_cells[f_id_sel][0];
        cs_lnum_t c_id2 = i_face_cells[f_id_sel][1];

        cs_real_2_t bi_bterms = {0.,0.};

        cs_i_cd_unsteady_upwind(ircflp,
                                weight[f_id],
                                cell_cen[c_id1],
                                cell_cen[c_id2],
                                i_face_cog[f_id_sel],
                                dijpf[f_id_sel],
                                grad[c_id1],
                                grad[c_id2],
                                f->val[c_id1],
                                f->val[c_id2],
                                &pifri,
                                &pifrj,
                                &pjfri,
                                &pjfrj,
                                &pip,
                                &pjp,
                                &pipr,
                                &pjpr);

        cs_i_conv_flux_cons(iconvp,
                            pifri,
                            pifrj,
                            pjfri,
                            pjfrj,
                            i_mass_flux[f_id_sel],
                            cpro_cp[c_id1],
                            cpro_cp[c_id2],
                            bi_bterms);

        cs_i_diff_flux(idiffp,
                       pip,
                       pjp,
                       pipr,
                       pjpr,
                       i_visc[f_id_sel],
                       bi_bterms);

        /* (The cell is counted only once in parallel by checking that
           the c_id is not in the halo) */
        /* Face normal well oriented (check bi_face_cells array) */
        if (bi_face_cells[f_id_sel][0] >= 0) {
          if (c_id1 < n_cells) {
            if (i_mass_flux[f_id_sel] > 0)
              bi_o_balance -= bi_bterms[0]*dt[c_id1];
            else
              bi_i_balance -= bi_bterms[0]*dt[c_id1];
          }
        }
        /* Face normal direction reversed */
        else {
          if (c_id2 < n_cells) {
            if (i_mass_flux[f_id_sel] > 0)
              bi_i_balance += bi_bterms[1]*dt[c_id2];
            else
              bi_o_balance += bi_bterms[1]*dt[c_id2];
          }
        }
      }
    }
  } else if (isstpp == 1) {
    /* ====================
       ---> No slope
       ====================*/

    /* Steady */
    if (idtvar < 0) {

      for (cs_lnum_t f_id = 0; f_id < n_bi_faces_sel; f_id++) {

        cs_lnum_t f_id_sel = bi_face_sel_list[f_id];
        /* Associated boundary-internal cells */
        cs_lnum_t c_id1 = i_face_cells[f_id_sel][0];
        cs_lnum_t c_id2 = i_face_cells[f_id_sel][1];

        cs_real_2_t bi_bterms = {0.,0.};

        cs_i_cd_steady(ircflp,
                       ischcp,
                       relaxp,
                       blencp,
                       weight[f_id],
                       cell_cen[c_id1],
                       cell_cen[c_id2],
                       i_face_cog[f_id_sel],
                       dijpf[f_id_sel],
                       grad[c_id1],
                       grad[c_id2],
                       f->val[c_id1],
                       f->val[c_id2],
                       f->val_pre[c_id1],
                       f->val_pre[c_id2],
                       &pifri,
                       &pifrj,
                       &pjfri,
                       &pjfrj,
                       &pip,
                       &pjp,
                       &pipr,
                       &pjpr);

        cs_i_conv_flux_cons(iconvp,
                            pifri,
                            pifrj,
                            pjfri,
                            pjfrj,
                            i_mass_flux[f_id_sel],
                            cpro_cp[c_id1],
                            cpro_cp[c_id2],
                            bi_bterms);

        cs_i_diff_flux(idiffp,
                       pip,
                       pjp,
                       pipr,
                       pjpr,
                       i_visc[f_id_sel],
                       bi_bterms);

        /* (The cell is counted only once in parallel by checking that
           the c_id is not in the halo) */
        /* Face normal well oriented (check bi_face_cells array) */
        if (bi_face_cells[f_id_sel][0] >= 0) {
          if (c_id1 < n_cells) {
            if (i_mass_flux[f_id_sel] > 0)
              bi_o_balance -= bi_bterms[0]*dt[c_id1];
            else
              bi_i_balance -= bi_bterms[0]*dt[c_id1];
          }
        }
        /* Face normal direction reversed */
        else {
          if (c_id2 < n_cells) {
            if (i_mass_flux[f_id_sel] > 0)
              bi_i_balance += bi_bterms[1]*dt[c_id2];
            else
              bi_o_balance += bi_bterms[1]*dt[c_id2];
          }
        }
      }
    /* Unsteady */
    } else {

      for (cs_lnum_t f_id = 0; f_id < n_bi_faces_sel; f_id++) {

        cs_lnum_t f_id_sel = bi_face_sel_list[f_id];
        /* Associated boundary-internal cells */
        cs_lnum_t c_id1 = i_face_cells[f_id_sel][0];
        cs_lnum_t c_id2 = i_face_cells[f_id_sel][1];

        cs_real_2_t bi_bterms = {0.,0.};

        cs_i_cd_unsteady(ircflp,
                         ischcp,
                         blencp,
                         weight[f_id],
                         cell_cen[c_id1],
                         cell_cen[c_id2],
                         i_face_cog[f_id_sel],
                         dijpf[f_id_sel],
                         grad[c_id1],
                         grad[c_id2],
                         f->val[c_id1],
                         f->val[c_id2],
                         &pifri,
                         &pifrj,
                         &pjfri,
                         &pjfrj,
                         &pip,
                         &pjp,
                         &pipr,
                         &pjpr);

        cs_i_conv_flux_cons(iconvp,
                            pifri,
                            pifrj,
                            pjfri,
                            pjfrj,
                            i_mass_flux[f_id_sel],
                            cpro_cp[c_id1],
                            cpro_cp[c_id2],
                            bi_bterms);

        cs_i_diff_flux(idiffp,
                       pip,
                       pjp,
                       pipr,
                       pjpr,
                       i_visc[f_id_sel],
                       bi_bterms);

        /* (The cell is counted only once in parallel by checking that
           the c_id is not in the halo) */
        /* Face normal well oriented (check bi_face_cells array) */
        if (bi_face_cells[f_id_sel][0] >= 0) {
          if (c_id1 < n_cells) {
            if (i_mass_flux[f_id_sel] > 0)
              bi_o_balance -= bi_bterms[0]*dt[c_id1];
            else
              bi_i_balance -= bi_bterms[0]*dt[c_id1];
          }
        }
        /* Face normal direction reversed */
        else {
          if (c_id2 < n_cells) {
            if (i_mass_flux[f_id_sel] > 0)
              bi_i_balance += bi_bterms[1]*dt[c_id2];
            else
              bi_o_balance += bi_bterms[1]*dt[c_id2];
          }
        }
      }
    }
  } else {
    /* ====================
       ---> Slope test
       ====================*/

    /* Steady */
    if (idtvar < 0) {

      for (cs_lnum_t f_id = 0; f_id < n_bi_faces_sel; f_id++) {

        cs_lnum_t f_id_sel = bi_face_sel_list[f_id];
        /* Associated boundary-internal cells */
        cs_lnum_t c_id1 = i_face_cells[f_id_sel][0];
        cs_lnum_t c_id2 = i_face_cells[f_id_sel][1];

        cs_real_2_t bi_bterms = {0.,0.};

        cs_i_cd_steady_slope_test(&indic,
                                  ircflp,
                                  ischcp,
                                  relaxp,
                                  blencp,
                                  weight[f_id],
                                  i_dist[f_id_sel],
                                  i_face_surf[f_id_sel],
                                  cell_cen[c_id1],
                                  cell_cen[c_id2],
                                  i_face_normal[f_id_sel],
                                  i_face_cog[f_id_sel],
                                  dijpf[f_id_sel],
                                  i_mass_flux[f_id_sel],
                                  grad[c_id1],
                                  grad[c_id2],
                                  grdpa[c_id1],
                                  grdpa[c_id2],
                                  f->val[c_id1],
                                  f->val[c_id2],
                                  f->val_pre[c_id1],
                                  f->val_pre[c_id2],
                                  &pifri,
                                  &pifrj,
                                  &pjfri,
                                  &pjfrj,
                                  &pip,
                                  &pjp,
                                  &pipr,
                                  &pjpr);

        cs_i_conv_flux_cons(iconvp,
                            pifri,
                            pifrj,
                            pjfri,
                            pjfrj,
                            i_mass_flux[f_id_sel],
                            cpro_cp[c_id1],
                            cpro_cp[c_id2],
                            bi_bterms);

        cs_i_diff_flux(idiffp,
                       pip,
                       pjp,
                       pipr,
                       pjpr,
                       i_visc[f_id_sel],
                       bi_bterms);

        /* (The cell is counted only once in parallel by checking that
           the c_id is not in the halo) */
        /* Face normal well oriented (check bi_face_cells array) */
        if (bi_face_cells[f_id_sel][0] >= 0) {
          if (c_id1 < n_cells) {
            if (i_mass_flux[f_id_sel] > 0)
              bi_o_balance -= bi_bterms[0]*dt[c_id1];
            else
              bi_i_balance -= bi_bterms[0]*dt[c_id1];
          }
        }
        /* Face normal direction reversed */
        else {
          if (c_id2 < n_cells) {
            if (i_mass_flux[f_id_sel] > 0)
              bi_i_balance += bi_bterms[1]*dt[c_id2];
            else
              bi_o_balance += bi_bterms[1]*dt[c_id2];
          }
        }
      }
    /* Unsteady */
    } else {

      for (cs_lnum_t f_id = 0; f_id < n_bi_faces_sel; f_id++) {

        cs_lnum_t f_id_sel = bi_face_sel_list[f_id];
        /* Associated boundary-internal cells */
        cs_lnum_t c_id1 = i_face_cells[f_id_sel][0];
        cs_lnum_t c_id2 = i_face_cells[f_id_sel][1];

        cs_real_2_t bi_bterms = {0.,0.};

        cs_i_cd_unsteady_slope_test(&indic,
                                    ircflp,
                                    ischcp,
                                    blencp,
                                    weight[f_id],
                                    i_dist[f_id_sel],
                                    i_face_surf[f_id_sel],
                                    cell_cen[c_id1],
                                    cell_cen[c_id2],
                                    i_face_normal[f_id_sel],
                                    i_face_cog[f_id_sel],
                                    dijpf[f_id_sel],
                                    i_mass_flux[f_id_sel],
                                    grad[c_id1],
                                    grad[c_id2],
                                    grdpa[c_id1],
                                    grdpa[c_id2],
                                    f->val[c_id1],
                                    f->val[c_id2],
                                    &pifri,
                                    &pifrj,
                                    &pjfri,
                                    &pjfrj,
                                    &pip,
                                    &pjp,
                                    &pipr,
                                    &pjpr);

        cs_i_conv_flux_cons(iconvp,
                            pifri,
                            pifrj,
                            pjfri,
                            pjfrj,
                            i_mass_flux[f_id_sel],
                            cpro_cp[c_id1],
                            cpro_cp[c_id2],
                            bi_bterms);

        cs_i_diff_flux(idiffp,
                       pip,
                       pjp,
                       pipr,
                       pjpr,
                       i_visc[f_id_sel],
                       bi_bterms);

        /* (The cell is counted only once in parallel by checking that
           the c_id is not in the halo) */
        /* Face normal well oriented (check bi_face_cells array) */
        if (bi_face_cells[f_id_sel][0] >= 0) {
          if (c_id1 < n_cells) {
            if (i_mass_flux[f_id_sel] > 0)
              bi_o_balance -= bi_bterms[0]*dt[c_id1];
            else
              bi_i_balance -= bi_bterms[0]*dt[c_id1];
          }
        }
        /* Face normal direction reversed */
        else {
          if (c_id2 < n_cells) {
            if (i_mass_flux[f_id_sel] > 0)
              bi_i_balance += bi_bterms[1]*dt[c_id2];
            else
              bi_o_balance += bi_bterms[1]*dt[c_id2];
          }
        }
      }
    }
  }

  /* Free memory */

  BFT_FREE(grad);
  if (grdpa != NULL)
    BFT_FREE(grdpa);
  BFT_FREE(f_reconstructed);

  if (!itemperature || icp == -1)
    BFT_FREE(cpro_cp);
  if (kivisl == -1)
    BFT_FREE(c_visc);
  BFT_FREE(i_visc);
  BFT_FREE(b_visc);

  BFT_FREE(cells_sel_list);
  BFT_FREE(cells_tag_list);
  BFT_FREE(bi_face_cells);
  BFT_FREE(i_face_sel_list);
  BFT_FREE(bb_face_sel_list);
  BFT_FREE(bi_face_sel_list);

  /* Sum of values on all ranks (parallel calculations) */

  cs_parall_sum(1, CS_DOUBLE, &vol_balance);
  cs_parall_sum(1, CS_DOUBLE, &div_balance);
  cs_parall_sum(1, CS_DOUBLE, &mass_i_balance);
  cs_parall_sum(1, CS_DOUBLE, &mass_o_balance);
  cs_parall_sum(1, CS_DOUBLE, &bi_i_balance);
  cs_parall_sum(1, CS_DOUBLE, &bi_o_balance);
  cs_parall_sum(1, CS_DOUBLE, &in_balance);
  cs_parall_sum(1, CS_DOUBLE, &out_balance);
  cs_parall_sum(1, CS_DOUBLE, &sym_balance);
  cs_parall_sum(1, CS_DOUBLE, &s_wall_balance);
  cs_parall_sum(1, CS_DOUBLE, &r_wall_balance);
  cs_parall_sum(1, CS_DOUBLE, &cpl_balance);
  cs_parall_sum(1, CS_DOUBLE, &ndef_balance);

  /* --> Total balance
         ------------- */

  /* We add the different contributions calculated above */

  tot_balance = vol_balance + div_balance
              + bi_i_balance + bi_o_balance
              + in_balance + out_balance + sym_balance
              + s_wall_balance + r_wall_balance
              + mass_i_balance + mass_o_balance
              + cpl_balance + ndef_balance;

  /* 3. Write the balance at time step n
    ==================================== */

  bft_printf(_("   ** SCALAR BALANCE BY ZONE at iteration %6i\n"
               "   ---------------------------------------------\n"
               "------------------------------------------------------------\n"
               "   SCALAR: %s\n"
               "   ZONE SELECTION CRITERIA: \"%s\"\n"
               "------------------------------------------------------------\n"
               "   Volume       Divergence   Inj. Mass.   Suc. Mass.\n"
               "  %12.4e %12.4e %12.4e %12.4e\n"
               "------------------------------------------------------------\n"
               "   IB inlet     IB outlet    Inlet        Outlet\n"
               "  %12.4e %12.4e %12.4e %12.4e\n"
               "------------------------------------------------------------\n"
               "   Sym.         Smooth W.    Rough W.     Coupled      Undef. BC\n"
               "  %12.4e %12.4e %12.4e %12.4e %12.4e\n"
               "------------------------------------------------------------\n"
               "   Total\n"
               "  %12.4e\n"
               "------------------------------------------------------------\n\n"),
             nt_cur, scalar_name, selection_crit,
             vol_balance, div_balance, mass_i_balance, mass_o_balance,
             bi_i_balance, bi_o_balance, in_balance, out_balance, sym_balance,
             s_wall_balance, r_wall_balance, cpl_balance, ndef_balance,
             tot_balance);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
