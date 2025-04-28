/*============================================================================
 * Turbulent viscosity for LES models.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_printf.h"

#include "base/cs_array.h"
#include "base/cs_base.h"
#include "base/cs_dispatch.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_pointer.h"
#include "base/cs_field_operator.h"
#include "alge/cs_gradient.h"
#include "base/cs_halo.h"
#include "base/cs_halo_perio.h"
#include "turb/cs_les_filter.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_location.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_parall.h"
#include "base/cs_physical_constants.h"
#include "pprt/cs_physical_model.h"
#include "base/cs_prototypes.h"
#include "base/cs_reducers.h"
#include "turb/cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "turb/cs_les_mu_t.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Tensor to vector (t2v) and vector to tensor (v2t) mask arrays */

#undef _IV2T
#undef _JV2T
#define _IV2T {0, 1, 2, 0, 1, 0};
#define _JV2T {0, 1, 2, 1, 2, 2};

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of  Mij:Mij and Mij:Lij for dynamic Smagorinsky model
 *
 * Please refer to the
 * <a href="../../theory.pdf#dynsmago"><b>dynamic Smagorinsky model</b></a>
 * section of the theory guide for more informations.
 *
 * \param[out]    s_n           strain rate (sqrt(2SijSij))
 * \param[out]    sf_n          filtered strain rate
 * \param[out]    f_vel         filtered velocity
 * \param[out]    mijmij        Mij:Mij
 * \param[out]    mijlij        Mij:Lij
 *
 !*/
/*----------------------------------------------------------------------------*/

void
cs_les_mu_t_smago_dyn_prepare(cs_real_t  s_n[],
                              cs_real_t  sf_n[],
                              cs_real_t  f_vel[][3],
                              cs_real_t  mijmij[],
                              cs_real_t  mijlij[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_real_t *cell_vol = cs_glob_mesh_quantities->cell_vol;

  const int irovar = cs_glob_fluid_properties->irovar;

  constexpr cs_real_t c_1ov3 = 1./3.;

  const cs_real_t xlesfl = cs_turb_xlesfl;
  const cs_real_t ales = cs_turb_ales;
  const cs_real_t bles = cs_turb_bles;
  const cs_real_t xlesfd = cs_turb_xlesfd;

  const int iv2t[6] = _IV2T;
  const int jv2t[6] = _JV2T;

  cs_dispatch_context ctx;

  /*  Allocate some work arrays */

  cs_real_t *w0, *xro, *xrof;
  cs_real_6_t *mij;
  cs_real_33_t *gradv;

  CS_MALLOC_HD(gradv, n_cells_ext, cs_real_33_t, cs_alloc_mode);
  CS_MALLOC_HD(w0, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(mij, n_cells_ext, cs_real_6_t, cs_alloc_mode);
  CS_MALLOC_HD(xro, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(xrof, n_cells_ext, cs_real_t, cs_alloc_mode);

  /* Take into account variable density case: Favre filtering
   * Constant density case: Reynolds filtering */

  const cs_real_t *crom  = CS_F_(rho)->val;
  const cs_real_3_t *vel = (const cs_real_3_t *)CS_F_(vel)->val;

  if (irovar == 1)
    cs_array_copy<cs_real_t>(n_cells, crom, xro);
  else
    cs_arrays_set_value<cs_real_t, 1>(n_cells, 1.0, xro);

  /* In case of constant density, xrof is 1.0 */

  cs_les_filter_scalar(xro, xrof);

  /* Computation of velocity gradient and of
   * S11^2+S22^2+S33^2+2*(S12^2+S13^2+S23^2)
   *======================================== */

  /* Allocate temporary arrays for gradients computation */

  cs_real_6_t *w61, *w62;

  CS_MALLOC_HD(w61, n_cells_ext, cs_real_6_t, cs_alloc_mode);
  CS_MALLOC_HD(w62, n_cells_ext, cs_real_6_t, cs_alloc_mode);

  cs_field_gradient_vector(CS_F_(vel),
                           false, // no use_previous_t
                           1,     // inc
                           gradv);

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    /* gradv[c_id][uvw][xyz] */

    cs_real_t divu = cs_math_33_trace(gradv[c_id]);

    /* In the case of constant density, s11+s22+s33 is zero
     * TODO: this is not exactly true, we should always remove the trace  */

    mij[c_id][0] = gradv[c_id][0][0]-irovar*c_1ov3*divu;
    mij[c_id][1] = gradv[c_id][1][1]-irovar*c_1ov3*divu;
    mij[c_id][2] = gradv[c_id][2][2]-irovar*c_1ov3*divu;
    mij[c_id][3] = 0.5*(gradv[c_id][0][1]+gradv[c_id][1][0]);
    mij[c_id][4] = 0.5*(gradv[c_id][1][2]+gradv[c_id][2][1]);
    mij[c_id][5] = 0.5*(gradv[c_id][0][2]+gradv[c_id][2][0]);

    s_n[c_id] =
      sqrt(2. * cs_math_sym_33_sym_33_product_trace(mij[c_id], mij[c_id]));

    for (cs_lnum_t ij = 0; ij < 6; ij++)
      w62[c_id][ij] = xro[c_id] * mij[c_id][ij];

  });

  ctx.wait();

  CS_FREE(gradv);

  /* w62 temporarily contains rho*S */

  cs_les_filter_strided<6>(w62, w61);

  /* w61 <rho*S>/<rho>, sf_n is ||<rho*S>/<rho>|| */

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (cs_lnum_t i = 0; i < 6; i++)
      w61[c_id][i] /= xrof[c_id];

    sf_n[c_id] =
      sqrt(2. * cs_math_sym_33_sym_33_product_trace(w61[c_id], w61[c_id]));

    /* Here mij contains Sij^d
     *   S_n  contains ||S||
     *       sqrt(2)*sqrt(S11^2+S22^2+S33^2+2(S12^2+S13^2+S23^2))
     *   Sf_n contains ||SF||
     *       sqrt(2)*sqrt(S11F^2+S22F^2+S33F^2+2(S12F^2+S13F^2+S23F^2)) */

    /* Computation of Mij
     * ================== */

    const cs_real_t delta = xlesfl * pow(ales*cell_vol[c_id], bles);
    w0[c_id] = delta;
    for (cs_lnum_t ij = 0; ij < 6; ij++)
      mij[c_id][ij] *= -2.0 * xro[c_id] * cs_math_pow2(delta) * s_n[c_id];
  });

  ctx.wait();

  /* w62 now contains <-2*rho*delta^2*||S||*S> */
  cs_les_filter_strided<6>(mij, w62);

  /* Now compute final mij value: M_ij = alpha_ij - beta_ij */

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    const cs_real_t delta = w0[c_id];
    const cs_real_t deltaf = xlesfd * delta;
    for (cs_lnum_t ij = 0; ij < 6; ij++)
      mij[c_id][ij] = -2.0 * xro[c_id] * cs_math_pow2(deltaf)
                           * sf_n[c_id] * w61[c_id][ij]
                      - w62[c_id][ij];
  });

  /* Allocate work arrays */
  cs_real_6_t *lij;
  cs_real_6_t *rho_ui_uj;
  cs_real_6_t *w_t;
  cs_real_3_t *w_v;

  CS_MALLOC_HD(rho_ui_uj, n_cells_ext, cs_real_6_t, cs_alloc_mode);
  CS_MALLOC_HD(lij, n_cells_ext, cs_real_6_t, cs_alloc_mode);
  CS_MALLOC_HD(w_t, n_cells_ext, cs_real_6_t, cs_alloc_mode);
  CS_MALLOC_HD(w_v, n_cells_ext, cs_real_3_t, cs_alloc_mode);

  /* Filtering the velocity and its square */

  /* Second order moment  <rho u_i u_j> */
  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (cs_lnum_t ij = 0; ij < 6; ij++)
      w_t[c_id][ij] = xro[c_id]*vel[c_id][iv2t[ij]]*vel[c_id][jv2t[ij]];
  });

  ctx.wait();

  cs_les_filter_strided<6>(w_t, rho_ui_uj);

  /* <rho u_i>/rho */
  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (cs_lnum_t i = 0; i < 3; i++)
      w_v[c_id][i] = xro[c_id]*vel[c_id][i];
  });

  ctx.wait();

  cs_les_filter_strided<3>(w_v, f_vel);

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    for (cs_lnum_t i = 0; i < 3; i++)
      f_vel[c_id][i] /= xrof[c_id];

    /* Computation of Lij */
    for (cs_lnum_t ij = 0; ij < 6; ij++)
      lij[c_id][ij] = rho_ui_uj[c_id][ij]
        - xrof[c_id] * f_vel[c_id][iv2t[ij]] * f_vel[c_id][jv2t[ij]];

    /* Computation of Mij :: Lij */
    mijlij[c_id] = cs_math_sym_33_sym_33_product_trace(mij[c_id], lij[c_id]);

    /* Computation of Mij :: Mij */
    mijmij[c_id] = cs_math_sym_33_sym_33_product_trace(mij[c_id], mij[c_id]);

  });

  ctx.wait();

  /* Free memory */
  CS_FREE(w61);
  CS_FREE(w62);
  CS_FREE(gradv);
  CS_FREE(w0);
  CS_FREE(mij);
  CS_FREE(xro);
  CS_FREE(xrof);
  CS_FREE(rho_ui_uj);
  CS_FREE(lij);
  CS_FREE(w_t);
  CS_FREE(w_v);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of the turbulent viscosity for
 *        a dynamic Smagorinsky LES model
 *
 * \f[ smago = \dfrac{L_{ij}M_{ij}}{M_{ij}M_{ij}} \f]
 *
 * \f[ \mu_T = \rho smago L^2  \sqrt{2 S_{ij}S_{ij}} \f]
 * \f[ S_{ij} = \dfrac{\der{u_i}{x_j} + \der{u_j}{x_i}}{2}\f]
 *
 * Please refer to the
 * <a href="../../theory.pdf#dynsmago"><b>dynamic Smagorinsky model</b></a>
 * section of the theory guide for more informations.
 *
 !*/
/*----------------------------------------------------------------------------*/

void
cs_les_mu_t_smago_dyn(void)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_real_t *cell_vol = cs_glob_mesh_quantities->cell_vol;

  const int irovar = cs_glob_fluid_properties->irovar;

  cs_dispatch_context ctx;

  /* Initialization
   * ============== */

  /*  Map field arrays */

  cs_real_t *visct = CS_F_(mu_t)->val;
  const cs_real_t *crom = CS_F_(rho)->val;
  const cs_real_3_t *vel = (const cs_real_3_t *)CS_F_(vel)->val;

  cs_real_t *cpro_smago
    = cs_field_by_name("smagorinsky_constant^2")->val;

  /* For the calculation of the viscosity of the sub-mesh */

  const cs_real_t ales = cs_turb_ales;
  const cs_real_t bles = cs_turb_bles;
  const cs_real_t xsmgmx = cs_turb_csmago_max;
  const cs_real_t xsmgmn = cs_turb_csmago_min;
  const cs_real_t xlesfl = cs_turb_xlesfl;
  const cs_real_t xlesfd = cs_turb_xlesfd;
  const cs_real_t c0 = cs_turb_crij_c0;

  /*  Allocate some work arrays */

  cs_real_t *w0, *w1, *xro, *xrof;
  CS_MALLOC_HD(w0, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(w1, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(xro, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(xrof, n_cells_ext, cs_real_t, cs_alloc_mode);

  /* Take into account variable density case: Favre filtering
   * Constant density case: Reynolds filtering */

  if (irovar == 1)
    cs_array_copy<cs_real_t>(n_cells, crom, xro);
  else
    cs_arrays_set_value<cs_real_t, 1>(n_cells, 1.0, xro);

  /* In case of constant density, xrof always 1.0 */

  cs_les_filter_scalar(xro, xrof);

  /* Allocate work arrays */
  cs_real_t *s_n, *sf_n;
  cs_real_t *w2, *w3, *w4;
  cs_real_3_t *f_vel;

  CS_MALLOC_HD(s_n, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(sf_n, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(w2, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(w3, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(w4, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(f_vel, n_cells_ext, cs_real_3_t, cs_alloc_mode);

  /* Compute:
   *   s_n (aka sqrt(2SijSij))
   *   sf_n (aka ||<rho*S>/<rho>||),
   *   filtered vel
   *   Mij:Mij
   *   Lij:Mij
   */

  cs_les_mu_t_smago_dyn_prepare(s_n, sf_n, f_vel, w2, w1);

  /* By default we compute a local average of the numerator and of the
     denominator, then only compute  the quotient. The user can overwrite
     this in cs_user_physical_properties_turb_viscosity. */

  cs_les_filter_scalar(w1, w3);
  cs_les_filter_scalar(w2, w4);

  cs_gnum_t iclipc = 0;

  ctx.parallel_for_reduce_sum
    (n_cells, iclipc, [=] CS_F_HOST_DEVICE
     (cs_lnum_t c_id,
      CS_DISPATCH_REDUCER_TYPE(cs_gnum_t) &sum) {

    if (cs_math_fabs(w4[c_id]) <= cs_math_epzero)
      cpro_smago[c_id] = xsmgmx;
    else
      cpro_smago[c_id] = w3[c_id]/w4[c_id];

    /* Clipping of cpro_smago */

    if (cpro_smago[c_id] >= xsmgmx) {
      cpro_smago[c_id] = xsmgmx;
      sum += 1;
    }
    else if (cpro_smago[c_id] <= xsmgmn) {
      cpro_smago[c_id] = xsmgmn;
      sum += 1;
    }

    /* Computation of (dynamic) viscosity
     * ================================== */

    const cs_real_t coef = cpro_smago[c_id];
    const cs_real_t delta = xlesfl * pow(ales*cell_vol[c_id], bles);

    visct[c_id] = crom[c_id]*coef*cs_math_pow2(delta)*s_n[c_id];
  });

  /* Compute k_SGS and its dissipation if need (e.g. Lagrangian module) */
  cs_field_t *f_k = cs_field_by_name_try("k_sgs");
  cs_field_t *f_eps = cs_field_by_name_try("epsilon_sgs");

  if (f_k != nullptr && f_eps != nullptr) {
    cs_real_t *cvar_k = f_k->val;
    cs_real_t *cvar_eps = f_eps->val;
    const cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();
    cs_real_t viscl0 = phys_pro->viscl0; /* reference molecular viscosity */

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      const cs_real_t nu_t =   cs::max(visct[c_id], cs_math_epzero*viscl0)
                             / crom[c_id];

      const cs_real_t s = s_n[c_id];

      cvar_eps[c_id] = nu_t * cs_math_pow2(s);
      cvar_k[c_id] = nu_t * s * (1. + 1.5 * c0)/sqrt(c0);
    });

    ctx.wait();

    /* Parallel sync */
    cs_halo_sync_var(m->halo, CS_HALO_STANDARD, f_eps->val);
    cs_halo_sync_var(m->halo, CS_HALO_STANDARD, f_k->val);

    /* Update previous values */
    cs_field_current_to_previous(f_eps);
    cs_field_current_to_previous(f_k);
  }

  const cs_equation_param_t *eqp
    = cs_field_get_equation_param_const(CS_F_(vel));

  if (eqp->verbosity >= 1) {

    struct cs_data_1double_2float rd;
    struct cs_reduce_sum1double_min1float_max1float reducer;

    ctx.parallel_for_reduce
      (n_cells, rd, reducer,
       [=] CS_F_HOST_DEVICE (cs_lnum_t c_id, cs_data_1double_2float &res) {

      res.d[0] = cpro_smago[c_id]*cell_vol[c_id];
      res.r[0] = cpro_smago[c_id];
      res.r[1] = cpro_smago[c_id];
    });

    ctx.wait();

    cs_parall_min(1, CS_FLOAT, &rd.r[0]);
    cs_parall_max(1, CS_FLOAT, &rd.r[1]);
    cs_parall_sum(1, CS_DOUBLE, &rd.d[0]);
    cs_parall_counter(&iclipc, 1);

    cs_real_t smagmy = rd.d[0] / cs_glob_mesh_quantities->tot_vol;
    cs_real_t smagmi = rd.r[0];
    cs_real_t smagma = rd.r[1];

    cs_log_printf(CS_LOG_DEFAULT,
                  _("N. clippings of the Smagorinsky constant %lu\n"
                    " --- Information on the squared Smagorinsky constant\n"
                    " --------------------------------\n"
                    " Mean value  Min value  Max value\n"
                    " --------------------------------\n"
                    " %12.4e, %12.4e, %12.4e\n"
                    " --------------------------------\n"),
                  (unsigned long)iclipc, smagmy, smagmi, smagma);
  }

  /* Clipping of the turbulent viscosity
   * =================================== */

  /* Clip turbulent viscosity so that it is always positive. */

  cs_gnum_t n_clip_visct = 0;

  ctx.parallel_for_reduce_sum
    (n_cells, n_clip_visct, [=] CS_F_HOST_DEVICE
     (cs_lnum_t c_id,
      CS_DISPATCH_REDUCER_TYPE(cs_gnum_t) &sum) {

    if (visct[c_id] < 0.0) {
      visct[c_id] = 0.0;
      sum += 1;
    }
  });

  ctx.wait();

  if (eqp->verbosity >= 1) {
    cs_parall_counter(&n_clip_visct, 1);
    cs_log_printf(CS_LOG_DEFAULT,
                  _("N. clippings for turbulent viscosity (mu_t>0): %lu\n"),
                  (unsigned long)n_clip_visct);
  }

  /* Scalar turbulent model
   *======================= */

  /* In case of gas combustion, the SGS scalar flux constant and the turbulent
   * diffusivity are only evaluated with the mixture fraction, then applied
   * automatically to the other scalar equations */

  const int key_turb_diff = cs_field_key_id("turbulent_diffusivity_id");
  const int key_sgs_sca_coef = cs_field_key_id("sgs_scalar_flux_coef_id");

  const int n_fields = cs_field_n_fields();
  const int keysca = cs_field_key_id("scalar_id");
  const int kscavr = cs_field_key_id("first_moment_id");

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *fld = cs_field_by_id(f_id);

    int sc_id = cs_field_get_key_int(fld, keysca)-1;
    int iscavr = cs_field_get_key_int(fld, kscavr);
    /* For variance of a scalar, the turbulent diffusivity is not computed */
    if ((sc_id < 0) || (iscavr > -1))
      continue;

    /* For any scalar other than the mixture fraction in diffusion flames,
     * Dt is not computed either. TODO Soot may be an exception */
    if (cs_glob_physical_model_flag[CS_COMBUSTION_SLFM] > -1) {
      if (fld != CS_F_(fm))
        continue;
    }

    int t_dif_id = cs_field_get_key_int(fld, key_turb_diff);
    int sca_dync_id = cs_field_get_key_int(fld, key_sgs_sca_coef);

    if ((t_dif_id >= 0) && (sca_dync_id >= 0)) {
      cs_real_t *cpro_turb_diff = cs_field_by_id(t_dif_id)->val;
      cs_real_t *cpro_sca_dync = cs_field_by_id(sca_dync_id)->val;

      /* Compute the Mi for scalar
       * ========================= */

      cs_real_3_t *grads, *gradsf;
      CS_MALLOC_HD(grads, n_cells_ext, cs_real_3_t, cs_alloc_mode);
      CS_MALLOC_HD(gradsf, n_cells_ext, cs_real_3_t, cs_alloc_mode);

      const cs_field_bc_coeffs_t *bc_coeffs = fld->bc_coeffs;

      cs_real_t *cvar_sca = fld->val;
      cs_field_gradient_scalar(fld,
                               false,     /* use previous t   */
                               1,         /* not on increment */
                               grads);

      /* compute grad (<rho.Y>/<rho>) */
      cs_real_3_t *scami, *scamif;
      CS_MALLOC_HD(scami, n_cells_ext, cs_real_3_t, cs_alloc_mode);
      CS_MALLOC_HD(scamif, n_cells_ext, cs_real_3_t, cs_alloc_mode);

      cs_real_3_t *w_v, *f_sca_vel;
      CS_MALLOC_HD(w_v, n_cells_ext, cs_real_3_t, cs_alloc_mode);
      CS_MALLOC_HD(f_sca_vel, n_cells_ext, cs_real_3_t, cs_alloc_mode);

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        w0[c_id] = cvar_sca[c_id]*xro[c_id];
      });

      ctx.wait();

      cs_les_filter_scalar(w0, w4);

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        w4[c_id] /= xrof[c_id];
      });

      ctx.wait();

      const cs_equation_param_t *eqp_fld
        = cs_field_get_equation_param_const(fld);

      cs_halo_type_t halo_type = CS_HALO_STANDARD;
      cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

      cs_gradient_type_by_imrgra(eqp_fld->imrgra,
                                 &gradient_type,
                                 &halo_type);

      cs_gradient_scalar("Work array",
                         gradient_type,
                         halo_type,
                         1, /* inc */
                         eqp_fld->nswrgr,
                         0,
                         1, /* w_stride */
                         eqp_fld->verbosity,
                         static_cast<cs_gradient_limit_t>(eqp_fld->imligr),
                         eqp_fld->epsrgr,
                         eqp_fld->climgr,
                         nullptr,
                         bc_coeffs,
                         w4,
                         nullptr,
                         nullptr, /* internal coupling */
                         gradsf);

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        const cs_real_t delta = xlesfl * pow(ales*cell_vol[c_id], bles);
        for (cs_lnum_t i = 0; i < 3; i++)
          scami[c_id][i]
            = -xro[c_id]*cs_math_pow2(delta)*s_n[c_id]*grads[c_id][i];
      });

      ctx.wait();

      cs_les_filter_strided<3>(scami, scamif);

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        const cs_real_t deltaf  = xlesfd * pow(ales*cell_vol[c_id], bles);
        for (cs_lnum_t i = 0; i < 3; i++) {
          scami[c_id][i] = - cs_math_pow2(deltaf)*xrof[c_id]
                              * sf_n[c_id]*gradsf[c_id][i]
                           - scamif[c_id][i];

          /* Compute the Li for scalar
           * ========================= */

          /* rho*Y*vel */
          w_v[c_id][i] = xro[c_id]*vel[c_id][i]*cvar_sca[c_id];
        }

      });

      ctx.wait();

      cs_les_filter_strided<3>(w_v, f_sca_vel);

      CS_FREE(w_v);

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        /* filter(rho Y vel) - rho filter(rho vel)/rho filter(Y)/rho */
        cs_real_t variance[3];
        for (cs_lnum_t i = 0; i < 3; i++)
          variance[i] = f_sca_vel[c_id][i] - xrof[c_id]*f_vel[c_id][i]*w4[c_id];

        w1[c_id] = cs_math_3_dot_product(variance, scami[c_id]);
        w2[c_id] = cs_math_3_square_norm(scami[c_id]);
      });

      ctx.wait();

      CS_FREE(f_sca_vel);

      cs_les_filter_scalar(w1, w3);
      cs_les_filter_scalar(w2, w4);

      /*
       * Compute the SGS flux coefficient and SGS diffusivity
       * Cs >= 0, Dt >=0
       */

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        if (cs_math_fabs(w4[c_id]) <= cs_math_epzero)
          cpro_sca_dync[c_id] = 0.0;
        else
          cpro_sca_dync[c_id] = cs_math_fmax(w3[c_id]/w4[c_id], 0.0);

        const cs_real_t delta = xlesfl * pow(ales*cell_vol[c_id], bles);

        cpro_turb_diff[c_id] =   crom[c_id] * cpro_sca_dync[c_id]
                               * cs_math_pow2(delta) * s_n[c_id];
      });

      ctx.wait();

      CS_FREE(scami);
      CS_FREE(scamif);
      CS_FREE(grads);
      CS_FREE(gradsf);
    }
  }

  /* Free memory */
  CS_FREE(s_n);
  CS_FREE(sf_n);

  CS_FREE(f_vel);
  CS_FREE(w4);
  CS_FREE(w3);
  CS_FREE(w2);
  CS_FREE(w1);
  CS_FREE(w0);

  CS_FREE(xro);
  CS_FREE(xrof);
}

/*----------------------------------------------------------------------------*/
/*! \brief Computation of turbulent viscosity for a Smagorinsky LES model
 *
 * \f[ \mu_T = \rho (C_{S} l)^2  \sqrt{2 S_{ij}S_{ij}} \f]
 * \f[ S_{ij} = \dfrac{\der{u_i}{x_j} + \der{u_j}{x_i}}{2}\f]
 *
 !*/
/*----------------------------------------------------------------------------*/

void
cs_les_mu_t_smago_const(void)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_real_t *cell_vol = cs_glob_mesh_quantities->cell_vol;

  const cs_real_t xlesfl = cs_turb_xlesfl;
  const cs_real_t ales = cs_turb_ales;
  const cs_real_t bles = cs_turb_bles;
  const cs_real_t csmago = cs_turb_csmago;
  const cs_real_t c0 = cs_turb_crij_c0;

  cs_dispatch_context ctx;

  /* Initialization
   * ============== */

  cs_real_33_t *gradv;
  CS_MALLOC_HD(gradv, n_cells_ext, cs_real_33_t, cs_alloc_mode);

  cs_real_t *visct = CS_F_(mu_t)->val;
  const cs_real_t *crom = CS_F_(rho)->val;

  /* We need the velocity gradient */

  cs_field_gradient_vector(CS_F_(vel),
                           false,  /* no use_previous_t */
                           1,      /* inc */
                           gradv);

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    /* S2 = 2 S:S = 2 Sij Sij
     * Note: should be the deviatoric part */
    const cs_real_t s2 = 2. * cs_math_33_main_invariant_2(gradv[c_id]);

    /* Computation of (dynamic) viscosity */

    const cs_real_t delta = xlesfl * pow(ales*cell_vol[c_id], bles);

    visct[c_id] = crom[c_id] * cs_math_pow2(csmago * delta) * sqrt(s2);
  });

  /* Compute k_SGS and its dissipation if need (e.g. Lagrangian module) */
  cs_field_t *f_k = cs_field_by_name_try("k_sgs");
  cs_field_t *f_eps = cs_field_by_name_try("epsilon_sgs");

  if (f_k != nullptr && f_eps != nullptr) {
    cs_real_t *cvar_k = f_k->val;
    cs_real_t *cvar_eps = f_eps->val;
    const cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();
    cs_real_t viscl0 = phys_pro->viscl0; /* reference molecular viscosity */

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      const cs_real_t nu_t =   cs::max(visct[c_id], cs_math_epzero*viscl0)
                             / crom[c_id];

      const cs_real_t s = sqrt(2. * cs_math_33_main_invariant_2(gradv[c_id]));

      cvar_eps[c_id] = nu_t * cs_math_pow2(s);
      cvar_k[c_id] = nu_t * s * (1. + 1.5 * c0)/sqrt(c0);
    });

    ctx.wait();

    /* Parallel sync */
    cs_halo_sync_var(m->halo, CS_HALO_STANDARD, f_eps->val);
    cs_halo_sync_var(m->halo, CS_HALO_STANDARD, f_k->val);

    /* Update previous values */
    cs_field_current_to_previous(f_eps);
    cs_field_current_to_previous(f_k);
  }

  ctx.wait();

  /* Free memory */
  CS_FREE(gradv);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the turbulent viscosity for the WALE LES model.
 *
 * The turbulent viscosity is:
 * \f$ \mu_T = \rho (C_{wale} L)^2 * \dfrac{(\tens{S}:\tens{Sd})^{3/2}}
 *                                         {(\tens{S} :\tens{S})^(5/2)
 *                                         +(\tens{Sd}:\tens{Sd})^(5/4)} \f$
 * with \f$ \tens{S}  = \frac{1}{2}(\gradt \vect{u} + \transpose{\gradt \vect{u}})\f$
 * and  \f$ \tens{Sd} = \deviator{(\symmetric{(\tens{S}^2)})}\f$
 *
 !*/
/*----------------------------------------------------------------------------*/

void
cs_les_mu_t_wale(void)

{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_real_t *cell_vol = cs_glob_mesh_quantities->cell_vol;

  constexpr cs_real_t c_1ov3 = 1./3.;

  const cs_real_t xlesfl = cs_turb_xlesfl;
  const cs_real_t ales = cs_turb_ales;
  const cs_real_t bles = cs_turb_bles;
  const cs_real_t c0 = cs_turb_crij_c0;
  const cs_real_t cwale = cs_turb_cwale;

  cs_dispatch_context ctx;

  /* Initialization
   * ============== */

  cs_real_33_t *gradv;
  CS_MALLOC_HD(gradv, n_cells_ext, cs_real_33_t, cs_alloc_mode);

  cs_real_t *visct = CS_F_(mu_t)->val;
  const cs_real_t *crom = CS_F_(rho)->val;

  /* Computation of the velocity gradient
   * ==================================== */

  cs_field_gradient_vector(CS_F_(vel),
                           false,  /* no use_previous_t */
                           1,      /* inc */
                           gradv);

  /* Kronecker delta Dij */

  const cs_real_t tdeltij[3][3] = {{1., 0., 0.},
                                   {0., 1., 0.},
                                   {0., 0., 1.}};

  /* Compute k_SGS and its dissipation if need (e.g. Lagrangian module) */

  cs_field_t *f_k = cs_field_by_name_try("k_sgs");
  cs_field_t *f_eps = cs_field_by_name_try("epsilon_sgs");
  cs_real_t *s_eq = nullptr;

  /* Store inverse of the time scale if needed */
  if (f_k != nullptr && f_eps != nullptr)
    CS_MALLOC_HD(s_eq, n_cells, cs_real_t, cs_alloc_mode);

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    /* Dudx is interleaved, but not gradv...
     * gradv[c_id][xyz][uvw]
     * ====================================== */

    cs_real_t dudx[3][3], g2[3][3];

    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        dudx[i][j] = gradv[c_id][i][j];
    }

    cs_real_t s = 0., trace_g2 = 0.;

    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        /* s = 1/4 * (dUi/dXj + dUj/dXi) * (dUi/dXj + dUj/dXi) */
        s += 0.25 * cs_math_pow2(dudx[i][j] + dudx[j][i]);

        /* g2 is the square tensor of the velocity gradient */
        g2[i][j] = 0.0;
        for (cs_lnum_t k = 0; k < 3; k++)
          g2[i][j] += dudx[i][k]*dudx[k][j];
      }

      trace_g2 += g2[i][i];
    }

    cs_real_t sd = 0.;
    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        /* Traceless symmetric part of the square of the velocity gradient tensor
         *   Sijd =   0.5*(dUi/dXk dUk/dXj + dUj/dXk dUk/dXi)
         *          - 1/3 Dij dUk/dXl dUl/dXk */
        cs_real_t sijd = 0.5*(g2[i][j]+g2[j][i])-c_1ov3*tdeltij[i][j]*trace_g2;
        sd += cs_math_pow2(sijd);
      }
    }

    /* Turbulent viscosity
     * =================== */

    /* Turbulent inverse time scale
     *   = (Sijd Sijd)^3/2 / [ (Sij Sij)^5/2 + (Sijd Sijd)^5/4 ] */

    cs_real_t sinv = pow(s, 2.5) + pow(sd, 1.25);
    cs_real_t con = 0;
    if (sinv > 0)
      con = sqrt(2.) * pow(sd, 1.5)/sinv;

    if (s_eq != nullptr)
      s_eq[c_id] = con;

    cs_real_t delta = xlesfl * pow(ales*cell_vol[c_id], bles);

    visct[c_id] = crom[c_id] * cs_math_pow2(cwale * delta) * con;

  });

  ctx.wait();

  /* Compute k_SGS and its dissipation if need (e.g. Lagrangian module) */

  if (f_k != nullptr && f_eps != nullptr) {
    cs_real_t *cvar_k = f_k->val;
    cs_real_t *cvar_eps = f_eps->val;
    const cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();
    cs_real_t viscl0 = phys_pro->viscl0; /* reference molecular viscosity */

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      const cs_real_t nu_t =   cs::max(visct[c_id], cs_math_epzero*viscl0)
                             / crom[c_id];

      cs_real_t s = s_eq[c_id];

      cvar_eps[c_id] = nu_t * cs_math_pow2(s);
      cvar_k[c_id] = nu_t * s * (1. + 1.5 * c0) / sqrt(c0);
    });

    ctx.wait();

    /* Parallel sync */
    cs_halo_sync_var(m->halo, CS_HALO_STANDARD, f_eps->val);
    cs_halo_sync_var(m->halo, CS_HALO_STANDARD, f_k->val);

    /* Update previous values */
    cs_field_current_to_previous(f_eps);
    cs_field_current_to_previous(f_k);

    CS_FREE(s_eq);
  }

  /* Free memory */
  CS_FREE(gradv);
}


/*----------------------------------------------------------------------------*/
/*! \brief Computation of turbulent viscosity for
 *        a k-SGS LES model
 !*/
/*----------------------------------------------------------------------------*/

void
cs_les_mu_t_ksgs(void)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_real_t *cell_vol = cs_glob_mesh_quantities->cell_vol;

  const cs_real_t xlesfl = cs_turb_xlesfl;
  const cs_real_t ales = cs_turb_ales;
  const cs_real_t bles = cs_turb_bles;
  const cs_real_t csmago = cs_turb_csmago;

  cs_dispatch_context ctx;

  /* Initialization
   * ============== */

  cs_real_t *visct = CS_F_(mu_t)->val;
  const cs_real_t *crom = CS_F_(rho)->val;
  cs_real_t *cvar_k = CS_F_(k)->val;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    /* Computation of (dynamic) viscosity */
    const cs_real_t delta = xlesfl * pow(ales*cell_vol[c_id], bles);

    visct[c_id] = crom[c_id] * csmago * delta * sqrt(cvar_k[c_id]);
  });

  ctx.wait();
}

/*----------------------------------------------------------------------------*/
/*! \brief Computation of turbulent viscosity for
 *        a tau-SGS LES model
 !*/
/*----------------------------------------------------------------------------*/

void
cs_les_mu_t_tausgs(void)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_real_t *cell_vol = cs_glob_mesh_quantities->cell_vol;

  const cs_real_t xlesfl = cs_turb_xlesfl;
  const cs_real_t ales = cs_turb_ales;
  const cs_real_t bles = cs_turb_bles;
  const cs_real_t csmago = cs_turb_csmago;

  cs_dispatch_context ctx;

  /* Initialization
   * ============== */

  cs_real_t *visct =  CS_F_(mu_t)->val;
  const cs_real_t *crom  = CS_F_(rho)->val;
  cs_real_6_t *cvar_rij = (cs_real_6_t *)(CS_F_(rij)->val);

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    /* Computation of (dynamic) viscosity */

    const cs_real_t delta = xlesfl * pow(ales*cell_vol[c_id], bles);

    cs_real_t tke = 0.5 * cs_math_6_trace(cvar_rij[c_id]);
    visct[c_id] = crom[c_id] * csmago * delta * sqrt(tke);
  });

  ctx.wait();
}

/*----------------------------------------------------------------------------*/

#endif /* cplusplus */
