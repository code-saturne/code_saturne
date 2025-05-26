/*============================================================================
 * Convection-diffusion operators.
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

/*-----------------------------------------------------------------------------
 * Standard C and C++ library headers
 *----------------------------------------------------------------------------*/

#include <cmath>
#include <chrono>

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_array.h"
#include "alge/cs_blas.h"
#include "alge/cs_bad_cells_regularisation.h"
#include "base/cs_boundary_conditions.h"
#include "base/cs_dispatch.h"
#include "cdo/cs_equation_param.h"
#include "base/cs_halo.h"
#include "base/cs_log.h"
#include "base/cs_internal_coupling.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_operator.h"
#include "base/cs_field_pointer.h"
#include "alge/cs_gradient.h"
#include "alge/cs_gradient_boundary.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_parall.h"
#include "base/cs_parameters.h"
#include "base/cs_porous_model.h"
#include "base/cs_profiling.h"
#include "base/cs_prototypes.h"
#include "base/cs_timer.h"
#include "base/cs_velocity_pressure.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "alge/cs_convection_diffusion.h"
#include "alge/cs_convection_diffusion_priv.h"

/*----------------------------------------------------------------------------*/

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file  cs_convection_diffusion.cpp
 *
 * \brief Convection-diffusion operators.
 *
 * Please refer to the
 * <a href="../../theory.pdf#conv-diff"><b>convection-diffusion</b></a> section
 *  of the theory guide for more informations.
 */
/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return the equivalent heat transfer coefficient. If both terms are
 * below a given tolerance, 0. is returned.
 *
 * parameters:
 *   h1     <-- first exchange coefficient
 *   h2     <-- second exchange coefficient
 *
 * return:
 *   value of equivalent exchange coefficient
 *----------------------------------------------------------------------------*/

static inline cs_real_t
_calc_heq(cs_real_t h1,
          cs_real_t h2)
{
  const cs_real_t h_eps = 1.e-12;

  cs_real_t heq = 0.;
  if (h1 + h2 > h_eps)
    heq = h1 * h2 / (h1 + h2);

  return heq;
}

/*----------------------------------------------------------------------------
 * Return the denominator to build the beta blending coefficient of the
 * beta limiter (ensuring preservation of a given min/max pair of values).
 *
 * parameters:
 *   f           <-- pointer to field
 *   ctx         <-- Reference to dispatch context
 *   eqp         <-- associated equation parameters
 *   inc         <-- 0 if an increment, 1 otherwise
 *   denom_inf   --> computed denominator for the lower bound
 *   denom_sup   --> computed denominator for the upper bound
 *
 * return:
 *   pointer to local values array, or nullptr;
 *----------------------------------------------------------------------------*/

static void
_beta_limiter_denom(cs_field_t                 *f,
                    cs_dispatch_context        &ctx,
                    const cs_equation_param_t  *eqp,
                    const int                   inc,
                    cs_real_t         *restrict denom_inf,
                    cs_real_t         *restrict denom_sup)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_nreal_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_rreal_3_t *restrict diipf = fvq->diipf;
  const cs_rreal_3_t *restrict djjpf = fvq->djjpf;

  const cs_real_t *restrict pvar  = f->val;
  const cs_real_t *restrict pvara = f->val_pre;

  const int ischcp = eqp->ischcv;
  const int ircflp = eqp->ircflu;
  const cs_real_t thetap = eqp->theta;
  const cs_real_t blencp = eqp->blencv;

  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  const cs_real_t *restrict i_massflux =
    cs_field_by_id(cs_field_get_key_int(f, kimasf))->val;
  const cs_real_t *restrict b_massflux =
    cs_field_by_id(cs_field_get_key_int(f, kbmasf))->val;

  const cs_real_t *hybrid_blend;
  if (CS_F_(hybrid_blend) != nullptr)
    hybrid_blend = CS_F_(hybrid_blend)->val;
  else
    hybrid_blend = nullptr;

  /* select halo type according to field gradient method */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  cs_gradient_type_by_imrgra(eqp->imrgra,
                             &gradient_type,
                             &halo_type);

  /* for NVD / TVD schemes (including VoF schemes) */

  const int key_lim_choice = cs_field_key_id("limiter_choice");
  cs_nvd_type_t limiter_choice = CS_NVD_N_TYPES;

  cs_real_t *local_min = nullptr;
  cs_real_t *local_max = nullptr;
  cs_real_t *courant = nullptr;

  if (ischcp == 4) {
    /* get limiter choice */
    limiter_choice = (cs_nvd_type_t)(cs_field_get_key_int(f, key_lim_choice));

    /* local extrema computation */
    CS_MALLOC_HD(local_max, n_cells_ext, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(local_min, n_cells_ext, cs_real_t, cs_alloc_mode);
    cs_field_local_extrema_scalar(f->id,
                                  halo_type,
                                  local_max,
                                  local_min);

    /* cell Courant number computation */
    if (limiter_choice >= CS_NVD_VOF_HRIC) {
      CS_MALLOC_HD(courant, n_cells_ext, cs_real_t, cs_alloc_mode);
      cs_cell_courant_number(f, ctx, courant);
    }
  }

  cs_real_t *df_limiter = nullptr;
  int df_limiter_id =
    cs_field_get_key_int(f, cs_field_key_id("diffusion_limiter_id"));
  if (df_limiter_id > -1)
    df_limiter = cs_field_by_id(df_limiter_id)->val;

  /* step 1: gradient computation */
  cs_real_3_t *grdpa;  /* for the implicit part */
  cs_real_3_t *grdpaa; /* for the explicit part */

  CS_MALLOC_HD(grdpa, n_cells_ext, cs_real_3_t, cs_alloc_mode);
  CS_MALLOC_HD(grdpaa, n_cells_ext, cs_real_3_t, cs_alloc_mode);

  /* TODO check if we can remove this initialization, as it seems
     redundant with assignment in following calls */
  ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    grdpa[c_id][0] = 0.;
    grdpa[c_id][1] = 0.;
    grdpa[c_id][2] = 0.;

    grdpaa[c_id][0] = 0.;
    grdpaa[c_id][1] = 0.;
    grdpaa[c_id][2] = 0.;
  });
  ctx.wait();

  /* legacy SOLU Scheme or centered scheme */
  if (ischcp != 2) {

    cs_field_gradient_scalar(f,
                             false, /* use_previous_t */
                             inc,
                             grdpa);

    cs_field_gradient_scalar(f,
                             true, /* use_previous_t */
                             inc,
                             grdpaa);

  }
  /* pure SOLU scheme (upwind gradient reconstruction) */
  else if (ischcp == 2) {

    cs_upwind_gradient(f->id,
                       ctx,
                       inc,
                       halo_type,
                       f->bc_coeffs,
                       i_massflux,
                       b_massflux,
                       pvar,
                       grdpa);

    cs_upwind_gradient(f->id,
                       ctx,
                       inc,
                       halo_type,
                       f->bc_coeffs,
                       i_massflux,
                       b_massflux,
                       pvara,
                       grdpaa);

  }

  /* Step 2: Building of denominator */
  ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    denom_inf[c_id] = 0.;
    denom_sup[c_id] = 0.;
  });
  ctx.wait();

  /* ---> Contribution from interior faces */

  ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

    cs_lnum_t ii = i_face_cells[face_id][0];
    cs_lnum_t jj = i_face_cells[face_id][1];

    cs_real_t pi = pvar[ii];
    cs_real_t pj = pvar[jj];
    cs_real_t pia = pvara[ii];
    cs_real_t pja = pvara[jj];
    cs_real_t pif, pjf, pip, pjp;
    cs_real_t pifa, pjfa, pipa, pjpa;

    cs_real_t _i_massflux = i_massflux[face_id];

    cs_real_t hybrid_coef_ii = 0., hybrid_coef_jj = 0.;
    if (ischcp == 3) {
      hybrid_coef_ii = hybrid_blend[ii];
      hybrid_coef_jj = hybrid_blend[jj];
    }

    cs_real_t bldfrp = (cs_real_t) ircflp;
    /* Local limitation of the reconstruction */
    if (df_limiter != nullptr && ircflp > 0)
      bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]), 0.);

    if (ischcp == 4) {
      /* NVD/TVD family of high accuracy schemes */

      cs_lnum_t ic, id;

      /* Determine central and downwind sides w.r.t. current face */
      cs_central_downwind_cells(ii,
                                jj,
                                _i_massflux,
                                &ic,  /* central cell id */
                                &id); /* downwind cell id */

      cs_real_t courant_c = -1.;
      if (courant != nullptr)
        courant_c = courant[ic];

      cs_i_cd_unsteady_nvd(limiter_choice,
                           blencp,
                           cell_cen[ic],
                           cell_cen[id],
                           i_face_u_normal[face_id],
                           i_face_cog[face_id],
                           grdpa[ic],
                           pvar[ic],
                           pvar[id],
                           local_max[ic],
                           local_min[ic],
                           courant_c,
                           &pif,
                           &pjf);

      cs_i_cd_unsteady_nvd(limiter_choice,
                           blencp,
                           cell_cen[ic],
                           cell_cen[id],
                           i_face_u_normal[face_id],
                           i_face_cog[face_id],
                           grdpaa[ic],
                           pvara[ic],
                           pvara[id],
                           local_max[ic],
                           local_min[ic],
                           courant_c,
                           &pifa,
                           &pjfa);
    }
    else {
      /* Value at time n */
      cs_i_cd_unsteady(bldfrp,
                       ischcp,
                       blencp,
                       weight[face_id],
                       cell_cen[ii],
                       cell_cen[jj],
                       i_face_cog[face_id],
                       hybrid_coef_ii,
                       hybrid_coef_jj,
                       diipf[face_id],
                       djjpf[face_id],
                       grdpa[ii], /* Std gradient when needed */
                       grdpa[jj], /* Std gradient when needed */
                       grdpa[ii], /* Upwind gradient when needed */
                       grdpa[jj], /* Upwind gradient when needed */
                       pi,
                       pj,
                       &pif,
                       &pjf,
                       &pip,
                       &pjp);

      /* Value at time n-1 */
      cs_i_cd_unsteady(bldfrp,
                       ischcp,
                       blencp,
                       weight[face_id],
                       cell_cen[ii],
                       cell_cen[jj],
                       i_face_cog[face_id],
                       hybrid_coef_ii, /* FIXME use previous values
                                          of blending function */
                       hybrid_coef_jj, /* FIXME use previous values
                                          of blending function */
                       diipf[face_id],
                       djjpf[face_id],
                       grdpaa[ii], /* Std gradient when needed */
                       grdpaa[jj], /* Std gradient when needed */
                       grdpaa[ii], /* Upwind gradient when needed */
                       grdpaa[jj], /* Upwind gradient when needed */
                       pia,
                       pja,
                       &pifa,
                       &pjfa,
                       &pipa,
                       &pjpa);
    }

    cs_real_t flui = 0.5*(_i_massflux + cs::abs(_i_massflux));
    cs_real_t fluj = 0.5*(_i_massflux - cs::abs(_i_massflux));

    cs_real_t flux =   thetap  * (  (pif  - pi )*flui
                                  + (pjf  - pj )*fluj)
                     + (1.-thetap) * (  (pifa - pia)*flui
                                      + (pjfa - pja)*fluj);

    /* blending to prevent lower bound violation
       We need to take the positive part*/
    cs_real_t partii = 0.5*(flux + cs::abs(flux));
    cs_real_t partjj = 0.5*(flux - cs::abs(flux));

    cs_dispatch_sum(&denom_inf[ii],  partii, i_sum_type);
    cs_dispatch_sum(&denom_inf[jj], -partjj, i_sum_type);

    /* blending to prevent upper bound violation
       We need to take the negative part
       Note: the swap between ii and jj is due to the fact
       that an upwind value on (Y-bound) is equivalent to a
       downwind value on (bound-Y) */
    cs_dispatch_sum(&denom_sup[ii], -partjj, i_sum_type);
    cs_dispatch_sum(&denom_sup[jj],  partii, i_sum_type);
  });

  ctx.wait();

  //Free Gradient arrays
  CS_FREE_HD(local_min);
  CS_FREE_HD(local_max);
  CS_FREE_HD(courant);
  CS_FREE_HD(grdpa);
  CS_FREE_HD(grdpaa);
}

/*----------------------------------------------------------------------------
 * Return the diagonal part of the numerator to build the Min/max limiter.
 *
 * parameters:
 *   f           <-- pointer to field
 *   ctx         <-- Reference to dispatch context
 *   eqp         <-- associated equation parameters
 *   rovsdt      <-- rho * volume / dt
 *   num_inf     --> computed numerator for the lower bound
 *   num_sup     --> computed numerator for the upper bound
 *
 *----------------------------------------------------------------------------*/

static void
_beta_limiter_num(cs_field_t                 *f,
                  cs_dispatch_context        &ctx,
                  const cs_equation_param_t  *eqp,
                  const int                   inc,
                  const cs_real_t             rovsdt[],
                  cs_real_t                  *num_inf,
                  cs_real_t                  *num_sup)
{
  const cs_mesh_t  *m = cs_glob_mesh;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  /* Get option from the field */

  const cs_real_t *restrict pvara = f->val_pre;

  const cs_real_t *restrict coefap = f->bc_coeffs->a;
  const cs_real_t *restrict coefbp = f->bc_coeffs->b;

  int key_scamax_id = cs_field_key_id("max_scalar");
  int key_scamin_id = cs_field_key_id("min_scalar");

  cs_real_t scalar_max = cs_field_get_key_double(f, key_scamax_id);
  cs_real_t scalar_min = cs_field_get_key_double(f, key_scamin_id);

  const cs_real_t thetex =  1. - eqp->theta;

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  const cs_real_t *restrict i_massflux =
    cs_field_by_id( cs_field_get_key_int(f, kimasf) )->val;
  const cs_real_t *restrict b_massflux =
    cs_field_by_id( cs_field_get_key_int(f, kbmasf) )->val;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    num_inf[c_id] = rovsdt[c_id] * cs::max(pvara[c_id] - scalar_min, 0.);
    num_sup[c_id] = rovsdt[c_id] * cs::max(scalar_max - pvara[c_id], 0.);
  });
  //ctx.wait();

  /* ---> Contribution from interior faces */

  ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

    cs_lnum_t ii = i_face_cells[face_id][0];
    cs_lnum_t jj = i_face_cells[face_id][1];

    cs_real_t _i_massflux = i_massflux[face_id];
    cs_real_t flui = 0.5*(_i_massflux + cs::abs(_i_massflux));
    cs_real_t fluj = 0.5*(_i_massflux - cs::abs(_i_massflux));

    cs_real_t pi = cs::max(pvara[ii] - scalar_min, 0.);
    cs_real_t pj = cs::max(pvara[jj] - scalar_min, 0.);

    cs_real_t pii = cs::max(scalar_max - pvara[ii], 0.);
    cs_real_t pjj = cs::max(scalar_max - pvara[jj], 0.);

    cs_real_t flux_inf = thetex * (pi  * flui + pj  * fluj);
    cs_real_t flux_sup = thetex * (pii * flui + pjj * fluj);

    if (ii < n_cells) {
      cs_dispatch_sum(&num_inf[ii], -flux_inf, i_sum_type);
      cs_dispatch_sum(&num_sup[ii], -flux_sup, i_sum_type);
    }
    if (jj < n_cells) {
      cs_dispatch_sum(&num_inf[jj], flux_inf, i_sum_type);
      cs_dispatch_sum(&num_sup[jj], flux_sup, i_sum_type);
    }
  });

  //ctx.wait();

  /* ---> Contribution from boundary faces */

  ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

    cs_lnum_t ii = b_face_cells[face_id];

    cs_real_t _b_massflux = b_massflux[face_id];
    cs_real_t flui = 0.5*(_b_massflux + cs::abs(_b_massflux));
    cs_real_t fluf = 0.5*(_b_massflux - cs::abs(_b_massflux));

    cs_real_t face_value = inc*coefap[face_id] + coefbp[face_id]*pvara[ii];

    cs_real_t flux_inf
      = thetex * (  cs::max(pvara[ii]  - scalar_min, 0.) * flui
                  + cs::max(face_value - scalar_min, 0.) * fluf);

    cs_real_t flux_sup
      = thetex *(  cs::max(scalar_max - pvara[ii],  0.) * flui
                 + cs::max(scalar_max - face_value, 0.) * fluf);

    cs_dispatch_sum(&num_inf[ii], -flux_inf, b_sum_type);
    cs_dispatch_sum(&num_sup[ii], -flux_sup, b_sum_type);
  });

  ctx.wait();
}

/*----------------------------------------------------------------------------
 * Synchronize strided gradient ghost cell values.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   on_device,     <-- is data on device (GPU) ?
 *   halo_type      <-- halo type (extended or not)
 *   grad           --> gradient of a variable
 *----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
static void
_sync_strided_gradient_halo(const cs_mesh_t         *m,
                            cs_halo_type_t           halo_type,
                            [[maybe_unused]] bool    on_device,
                            cs_real_t (*restrict grad)[stride][3])
{
#if defined(HAVE_ACCEL)
  if (on_device)
    cs_halo_sync_d(m->halo, halo_type, CS_REAL_TYPE, stride*3,
                   (cs_real_t *)grad);
  else
#endif
    cs_halo_sync(m->halo, halo_type, CS_REAL_TYPE, stride*3,
                 (cs_real_t *)grad);

  if (m->have_rotation_perio) {
#if defined(HAVE_ACCEL)
    if (on_device)
      cs_sync_d2h((void  *)grad);
#endif
    if (stride == 1)
      cs_halo_perio_sync_var_vect(m->halo, halo_type, (cs_real_t *)grad, 3);
    else if (stride == 3)
      cs_halo_perio_sync_var_tens(m->halo, halo_type, (cs_real_t *)grad);
    else if (stride == 6)
      cs_halo_perio_sync_var_sym_tens_grad(m->halo,
                                           halo_type,
                                           (cs_real_t *)grad);
#if defined(HAVE_ACCEL)
    if (on_device)
      cs_sync_h2d((void  *)grad);
#endif
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the upwind gradient used in the slope tests, on host
 *
 * This function assumes the input gradient and pvar values have already
 * been synchronized.
 *
 * \param[in]     f_id         field id
 * \param[in]     ctx          Reference to dispatch context
 * \param[in]     inc          Not an increment flag
 * \param[in]     grad         standard gradient
 * \param[out]    grdpa        upwind gradient
 * \param[in]     pvar         values
 * \param[in]     bc_coeffs    boundary condition structure for the variable
 * \param[in]     i_massflux   mass flux at interior faces
 */
/*----------------------------------------------------------------------------*/

static void
_slope_test_gradient_h
  (cs_host_context            &ctx,
   int                         inc,
   const cs_real_3_t          *grad,
   cs_real_3_t                *restrict grdpa,
   const cs_real_t            *pvar,
   const cs_field_bc_coeffs_t *bc_coeffs,
   const cs_real_t            *i_massflux)
{
  cs_real_t *coefap = bc_coeffs->a;
  cs_real_t *coefbp = bc_coeffs->b;

  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_real_t *restrict i_face_surf = fvq->i_face_surf;
  const cs_real_t *restrict b_face_surf = fvq->b_face_surf;
  const cs_nreal_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;
  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;

  /* Cast to the parent class to obtain info on the parent class.
     On host, we know this is not necessary and we use a simple
     sum, but with this precaution, we can enable this function on
     just by changing the "ctx" argument type. */

  cs_dispatch_context &p_ctx = static_cast<cs_dispatch_context&>(ctx);
  cs_dispatch_sum_type_t i_sum_type = p_ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = p_ctx.get_parallel_for_b_faces_sum_type(m);

  ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

    cs_lnum_t ii = i_face_cells[face_id][0];
    cs_lnum_t jj = i_face_cells[face_id][1];

    cs_real_t difv[3] = {i_face_cog[face_id][0] - cell_cen[ii][0],
                         i_face_cog[face_id][1] - cell_cen[ii][1],
                         i_face_cog[face_id][2] - cell_cen[ii][2]};

    cs_real_t djfv[3] = {i_face_cog[face_id][0] - cell_cen[jj][0],
                         i_face_cog[face_id][1] - cell_cen[jj][1],
                         i_face_cog[face_id][2] - cell_cen[jj][2]};

    cs_real_t pif = pvar[ii] + cs_math_3_dot_product(difv, grad[ii]);
    cs_real_t pjf = pvar[jj] + cs_math_3_dot_product(djfv, grad[jj]);

    cs_real_t pfac = (i_massflux[face_id] > 0.) ? pif : pjf;
    pfac *= i_face_surf[face_id];

    cs_real_t vfac_i[3], vfac_j[3];
    for (cs_lnum_t k = 0; k < 3; k++) {
      vfac_i[k] = pfac*i_face_u_normal[face_id][k];
      vfac_j[k] = - vfac_i[k];
    }

    cs_dispatch_sum<3>(grdpa[ii], vfac_i, i_sum_type);
    cs_dispatch_sum<3>(grdpa[jj], vfac_j, i_sum_type);

  });

  ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

    cs_lnum_t ii = b_face_cells[face_id];

    cs_real_t pfac = inc*coefap[face_id] + coefbp[face_id]
      * (pvar[ii] + cs_math_3_dot_product(grad[ii], diipb[face_id]));

    cs_real_t vfac[3] = {pfac*b_face_surf[face_id]*b_face_u_normal[face_id][0],
                         pfac*b_face_surf[face_id]*b_face_u_normal[face_id][1],
                         pfac*b_face_surf[face_id]*b_face_u_normal[face_id][2]};

    cs_dispatch_sum<3>(grdpa[ii], vfac, b_sum_type);

  });

  const int *restrict c_disable_flag = (fvq->has_disable_flag) ?
    fvq->c_disable_flag : nullptr;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
    cs_real_t unsvol = cs_mq_cell_vol_inv(cell_id, c_disable_flag, cell_vol);

    grdpa[cell_id][0] *= unsvol;
    grdpa[cell_id][1] *= unsvol;
    grdpa[cell_id][2] *= unsvol;
  });
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the upwind gradient used in the slope tests, on device
 *
 * This function assumes the input gradient and pvar values have already
 * been synchronized.
 *
 * \param[in]     ctx          Reference to dispatch context
 * \param[in]     inc          Not an increment flag
 * \param[in]     grad         standard gradient
 * \param[out]    grdpa        upwind gradient
 * \param[in]     pvar         values
 * \param[in]     bc_coeffs    boundary condition structure for the variable
 * \param[in]     i_massflux   mass flux at interior faces
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

static void
_slope_test_gradient_d
  (cs_device_context          &ctx,
   int                         inc,
   const cs_real_3_t          *grad,
   cs_real_3_t                *restrict grdpa,
   const cs_real_t            *pvar,
   const cs_field_bc_coeffs_t *bc_coeffs,
   const cs_real_t            *i_massflux)
{
  cs_real_t *coefap = bc_coeffs->a;
  cs_real_t *coefbp = bc_coeffs->b;

  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;

  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_real_t *restrict i_face_surf = fvq->i_face_surf;
  const cs_real_t *restrict b_face_surf = fvq->b_face_surf;
  const cs_nreal_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;
  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  cs_mesh_adjacencies_update_cell_i_faces();
  const cs_lnum_t *c2c = ma->cell_cells;
  const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
  const short int *c2f_sgn = ma->cell_i_faces_sgn;
  const cs_lnum_t *cell_i_faces = ma->cell_i_faces;

  const int *restrict c_disable_flag = (fvq->has_disable_flag) ?
    fvq->c_disable_flag : nullptr;

  /* Cast to the parent class to obtain info on the parent class.
     On device, we know this is not necessary and we use an atomic
     sum, but with this precaution, we can enable this function on
     just by changing the "ctx" argument type. */
  cs_dispatch_context &p_ctx = static_cast<cs_dispatch_context&>(ctx);
  cs_dispatch_sum_type_t b_sum_type = p_ctx.get_parallel_for_b_faces_sum_type(m);

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {

    cs_real_t grdpa_c[3] = {0, 0, 0};

    /* Loop on interior faces */
    const cs_lnum_t s_id_i = c2c_idx[cell_id];
    const cs_lnum_t e_id_i = c2c_idx[cell_id + 1];

    for (cs_lnum_t cidx = s_id_i; cidx < e_id_i; cidx++) {
      const cs_lnum_t face_id = cell_i_faces[cidx];

      /* Which cell is upwind ? */
      cs_lnum_t u_cell_id = cell_id;
      short int f_sgn = c2f_sgn[cidx];
      if (f_sgn*i_massflux[face_id] <= 0.)
        u_cell_id = c2c[cidx];

      cs_real_t dufv[3];
      for (cs_lnum_t isou = 0; isou < 3; isou++)
        dufv[isou] = i_face_cog[face_id][isou] - cell_cen[u_cell_id][isou];

      cs_real_t pfac =   pvar[u_cell_id]
                       + cs_math_3_dot_product(dufv, grad[u_cell_id]);

      pfac *= i_face_surf[face_id] * f_sgn;

      for (cs_lnum_t isou = 0; isou < 3; isou++)
        grdpa_c[isou] += pfac*i_face_u_normal[face_id][isou];

      /* Scale now to avoid second loop */

      cs_real_t unsvol = cs_mq_cell_vol_inv(cell_id, c_disable_flag, cell_vol);
      for (cs_lnum_t isou = 0; isou < 3; isou++) {
        grdpa[cell_id][isou] = grdpa_c[isou]*unsvol;
      }

    }
  });

  /* Contribution from boundary faces

     TODO! test folding this into previous kernel when BC coefficients
     are handled using the new scheme. For now, we keep this loop
     separate, as handling of BC coefficients is more costly so this
     could add a lot of imbalance between cells in the main kernel */

  ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

    cs_lnum_t ii = b_face_cells[face_id];

    cs_real_t unsvol = cs_mq_cell_vol_inv(ii, c_disable_flag, cell_vol);
    const cs_real_t _b_face_surf_o_v = b_face_surf[face_id] * unsvol;

    cs_real_t pfac = inc*coefap[face_id] + coefbp[face_id]
      * (pvar[ii] + cs_math_3_dot_product(grad[ii], diipb[face_id]));

    cs_real_t vfac[3] = {pfac*_b_face_surf_o_v*b_face_u_normal[face_id][0],
                         pfac*_b_face_surf_o_v*b_face_u_normal[face_id][1],
                         pfac*_b_face_surf_o_v*b_face_u_normal[face_id][2]};

    cs_dispatch_sum<3>(grdpa[ii], vfac, b_sum_type);

  });
}

#endif /* defined(HAVE_ACCEL) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the upwind gradient used in the slope tests, on host
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * This function assumes the input gradient and pvar values have already
 * been synchronized.
 *
 * \param[in]     ctx          Reference to dispatch context
 * \param[in]     halo_type    halo type
 * \param[in]     grad         standard gradient
 * \param[out]    grdpa        upwind gradient
 * \param[in]     pvar         values
 * \param[in]     val_f        face values for gradient
 * \param[in]     i_massflux   mass flux at interior faces
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
static void
_slope_test_gradient_strided_h
  (cs_host_context             &ctx,
   const cs_real_t              grad[][stride][3],
   cs_real_t                  (*restrict grdpa)[stride][3],
   const cs_real_t              pvar[][stride],
   const cs_real_t              val_f[][stride],
   const cs_real_t             *i_massflux)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_nreal_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;
  const cs_real_t *restrict i_f_face_surf = fvq->i_face_surf;
  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
  const cs_real_t *restrict b_f_face_surf = fvq->b_face_surf;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;

  /* Cast to the parent class to obtain info on the parent class.
     On host, we know this is not necessary and we use a simple
     sum, but with this precaution, we can enable this function on
     just by changing the "ctx" argument type. */

  cs_dispatch_context &p_ctx = static_cast<cs_dispatch_context&>(ctx);
  cs_dispatch_sum_type_t i_sum_type = p_ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = p_ctx.get_parallel_for_b_faces_sum_type(m);

  ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
    for (cs_lnum_t isou = 0; isou < stride; isou++) {
      for (cs_lnum_t jsou = 0; jsou < 3; jsou++)
        grdpa[cell_id][isou][jsou] = 0.;
    }
  });

  ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

    cs_real_t difv[3], djfv[3];

    cs_lnum_t ii = i_face_cells[face_id][0];
    cs_lnum_t jj = i_face_cells[face_id][1];

    for (cs_lnum_t jsou = 0; jsou < 3; jsou++) {
      difv[jsou] = i_face_cog[face_id][jsou] - cell_cen[ii][jsou];
      djfv[jsou] = i_face_cog[face_id][jsou] - cell_cen[jj][jsou];
    }

    /* x-y-z component, p = u, v, w */

    for (cs_lnum_t isou = 0; isou < stride; isou++) {
      cs_real_t pif = pvar[ii][isou];
      cs_real_t pjf = pvar[jj][isou];
      for (cs_lnum_t jsou = 0; jsou < 3; jsou++) {
        pif = pif + grad[ii][isou][jsou]*difv[jsou];
        pjf = pjf + grad[jj][isou][jsou]*djfv[jsou];
      }

      cs_real_t pfac = pjf;
      if (i_massflux[face_id] > 0.) pfac = pif;

      /* U gradient */

      pfac *= i_f_face_surf[face_id];
      cs_real_t vfac_i[3], vfac_j[3];

      for (cs_lnum_t jsou = 0; jsou < 3; jsou++) {
        vfac_i[jsou] = pfac*i_face_u_normal[face_id][jsou];
        vfac_j[jsou] = - vfac_i[jsou];
      }

      cs_dispatch_sum<3>(grdpa[ii][isou], vfac_i, i_sum_type);
      cs_dispatch_sum<3>(grdpa[jj][isou], vfac_j, i_sum_type);
    }

  });

  ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

    cs_real_t vfac[3];
    cs_lnum_t ii = b_face_cells[face_id];

    const cs_real_t &_b_f_face_surf = b_f_face_surf[face_id];

    for (cs_lnum_t isou = 0; isou < stride; isou++) {
      for (cs_lnum_t jsou =  0; jsou < 3; jsou++)
        vfac[jsou] =   val_f[face_id][isou] * _b_f_face_surf
                     * b_face_u_normal[face_id][jsou];

      cs_dispatch_sum<3>(grdpa[ii][isou], vfac, b_sum_type);
    }

  });

  const int *restrict c_disable_flag = (fvq->has_disable_flag) ?
    fvq->c_disable_flag : nullptr;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
    cs_real_t unsvol = cs_mq_cell_vol_inv(cell_id, c_disable_flag, cell_vol);
    for (cs_lnum_t isou = 0; isou < stride; isou++) {
      for (cs_lnum_t jsou = 0; jsou < 3; jsou++)
        grdpa[cell_id][isou][jsou] = grdpa[cell_id][isou][jsou]*unsvol;
    }
  });
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the upwind gradient used in the slope tests, on device
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * This function assumes the input gradient and pvar values have already
 * been synchronized.
 *
 * \param[in]     ctx          Reference to dispatch context
 * \param[in]     halo_type    halo type
 * \param[in]     grad         standard gradient
 * \param[out]    grdpa        upwind gradient
 * \param[in]     pvar         values
 * \param[in]     val_f        face values for gradient
 * \param[in]     i_massflux   mass flux at interior faces
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

template <cs_lnum_t stride>
static void
_slope_test_gradient_strided_d
  (cs_device_context           &ctx,
   const cs_real_t              grad[][stride][3],
   cs_real_t                  (*restrict grdpa)[stride][3],
   const cs_real_t              pvar[][stride],
   const cs_real_t              val_f[][stride],
   const cs_real_t             *i_massflux)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;

  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_nreal_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;
  const cs_real_t *restrict i_f_face_surf = fvq->i_face_surf;
  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
  const cs_real_t *restrict b_f_face_surf = fvq->b_face_surf;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const int *restrict c_disable_flag = (fvq->has_disable_flag) ?
    fvq->c_disable_flag : nullptr;

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  cs_mesh_adjacencies_update_cell_i_faces();
  const cs_lnum_t *c2c = ma->cell_cells;
  const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
  const short int *c2f_sgn = ma->cell_i_faces_sgn;
  const cs_lnum_t *cell_i_faces = ma->cell_i_faces;

  /* Cast to the parent class to obtain info on the parent class.
     On device, we know this is not necessary and we use an atomic
     sum, but with this precaution, we can enable this function on
     just by changing the "ctx" argument type. */
  cs_dispatch_context &p_ctx = static_cast<cs_dispatch_context&>(ctx);
  cs_dispatch_sum_type_t b_sum_type = p_ctx.get_parallel_for_b_faces_sum_type(m);

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {

    cs_real_t grdpa_c[stride][3];

    for (cs_lnum_t isou = 0; isou < stride; isou++) {
      for (cs_lnum_t jsou = 0; jsou < 3; jsou++)
        grdpa_c[isou][jsou] = 0.;
    }

    /* Loop on interior faces */
    const cs_lnum_t s_id_i = c2c_idx[cell_id];
    const cs_lnum_t e_id_i = c2c_idx[cell_id + 1];

    for (cs_lnum_t cidx = s_id_i; cidx < e_id_i; cidx++) {
      const cs_lnum_t face_id = cell_i_faces[cidx];

      /* Which cell is upwind ? */
      cs_lnum_t u_cell_id = cell_id;
      short int f_sgn = c2f_sgn[cidx];
      if (f_sgn*i_massflux[face_id] <= 0.)
        u_cell_id = c2c[cidx];

      cs_real_t dufv[3];
      for (cs_lnum_t jsou = 0; jsou < 3; jsou++)
        dufv[jsou] = i_face_cog[face_id][jsou] - cell_cen[u_cell_id][jsou];

      /* For each component */

      for (cs_lnum_t isou = 0; isou < stride; isou++) {
        cs_real_t pfac = pvar[u_cell_id][isou];
        for (cs_lnum_t jsou = 0; jsou < 3; jsou++) {
          pfac += grad[u_cell_id][isou][jsou]*dufv[jsou];
        }

        /* U gradient */

        pfac *= i_f_face_surf[face_id] * f_sgn;

        for (cs_lnum_t jsou = 0; jsou < 3; jsou++)
          grdpa_c[isou][jsou] += pfac*i_face_u_normal[face_id][jsou];
      }
    }

    /* Scale now to avoid second loop */

    cs_real_t unsvol = cs_mq_cell_vol_inv(cell_id, c_disable_flag, cell_vol);
    for (cs_lnum_t isou = 0; isou < stride; isou++) {
      for (cs_lnum_t jsou = 0; jsou < 3; jsou++)
        grdpa[cell_id][isou][jsou] = grdpa_c[isou][jsou]*unsvol;
    }

  });

  /* Contribution from boundary faces

     TODO! test folding this into previous kernel when BC coefficients
     are handled using the new scheme. For now, we keep this loop
     separate, as handling of BC coefficients is more costly so this
     could add a lot of imbalance between cells in the main kernel */

  ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

    cs_lnum_t ii = b_face_cells[face_id];
    cs_real_t vfac[3];

    /* x-y-z components, p = u, v, w */

    cs_real_t unsvol = cs_mq_cell_vol_inv(ii, c_disable_flag, cell_vol);
    const cs_real_t _b_f_face_surf_o_v = b_f_face_surf[face_id] * unsvol;

    for (cs_lnum_t isou = 0; isou < stride; isou++) {
      for (cs_lnum_t jsou = 0; jsou < 3; jsou++)
        vfac[jsou] =   val_f[face_id][isou] * _b_f_face_surf_o_v
                     * b_face_u_normal[face_id][jsou];

      cs_dispatch_sum<3>(grdpa[ii][isou], vfac, b_sum_type);
    }

  });
}

#endif /* defined(HAVE_ACCEL) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the upwind gradient used in the slope tests.
 *
 * template parameters:
 *   stride        1 for scalars, 3 for vectors, 6 for symmetric tensors
 *
 * This function assumes the input gradient and pvar values have already
 * been synchronized.
 *
 * \param[in]     ctx          Reference to dispatch context
 * \param[in]     grad         standard gradient
 * \param[out]    grdpa        upwind gradient
 * \param[in]     pvar         values
 * \param[in]     val_f        face values for gradient
 * \param[in]     i_massflux   mass flux at interior faces
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
static void
_slope_test_gradient_strided
  (cs_dispatch_context         &ctx,
   const cs_real_t              grad[][stride][3],
   cs_real_t                  (*restrict grdpa)[stride][3],
   const cs_real_t              pvar[][stride],
   const cs_real_t              val_f[][stride],
   const cs_real_t             *i_massflux)
{
  bool use_gpu = ctx.use_gpu();
  const cs_mesh_t  *m = cs_glob_mesh;

  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

#if defined(HAVE_ACCEL)

  if (use_gpu) {
    cs_device_context &d_ctx = static_cast<cs_device_context&>(ctx);

    _slope_test_gradient_strided_d<stride>
      (d_ctx,
       grad,
       grdpa,
       pvar,
       val_f,
       i_massflux);
  }

#endif

  if (use_gpu == false) {
    cs_host_context &h_ctx = static_cast<cs_host_context&>(ctx);

    _slope_test_gradient_strided_h<stride>
      (h_ctx,
       grad,
       grdpa,
       pvar,
       val_f,
       i_massflux);
  }

  ctx.wait();

  /* Handle parallelism and periodicity */

  if (m->halo != nullptr)
    _sync_strided_gradient_halo<stride>(m,
                                        CS_HALO_STANDARD,
                                        use_gpu,
                                        grdpa);

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();

    std::chrono::microseconds elapsed;
    printf("%d: %s<%d>", cs_glob_rank_id, __func__, stride);

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_stop - t_start);
    printf(", total = %ld\n", elapsed.count());
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the convection/diffusion terms of a
 * standard transport equation of a scalar field \f$ \varia \f$.
 *
 * More precisely, the right hand side \f$ Rhs \f$ is updated as
 * follows:
 * \f[
 * Rhs = Rhs - \sum_{\fij \in \Facei{\celli}}      \left(
 *        \dot{m}_\ij \left( \varia_\fij - \varia_\celli \right)
 *      - \mu_\fij \gradv_\fij \varia \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ Rhs \f$ has already been initialized before calling bilsc2!
 * - mind the sign minus
 *
 * Please refer to the
 * <a href="../../theory.pdf#bilsc2"><b>bilsc2</b></a> section of the
 * theory guide for more informations.
 *
 * \param[in]     f             pointer to field
 * \param[in]     eqp           equation parameters
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in]     bc_coeffs     boundary condition structure for the variable
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     xcpp          array of specific heat (\f$ C_p \f$), or nullptr
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

static void
_convection_diffusion_scalar_steady(const cs_field_t           *f,
                                    const cs_equation_param_t   eqp,
                                    bool                        icvflb,
                                    int                         inc,
                                    cs_real_t         *restrict pvar,
                                    const cs_real_t   *restrict pvara,
                                    const int                   icvfli[],
                                    const cs_field_bc_coeffs_t *bc_coeffs,
                                    const cs_real_t             i_massflux[],
                                    const cs_real_t             b_massflux[],
                                    const cs_real_t             i_visc[],
                                    const cs_real_t             b_visc[],
                                    const cs_real_t             xcpp[],
                                    cs_real_t         *restrict rhs)
{

  const cs_real_t *coefap = bc_coeffs->a;
  const cs_real_t *coefbp = bc_coeffs->b;
  const cs_real_t *cofafp = bc_coeffs->af;
  const cs_real_t *cofbfp = bc_coeffs->bf;

  const int iconvp = eqp.iconv;
  const int idiffp = eqp.idiff;
  const int ircflp = eqp.ircflu;
  const int ircflb = (ircflp > 0) ? eqp.b_diff_flux_rc : 0;
  const int ischcp = eqp.ischcv;
  const int isstpp = eqp.isstpc;
  const int iwarnp = eqp.verbosity;
  const int icoupl = eqp.icoupl;
  const double blencp = eqp.blencv;
  const double blend_st = eqp.blend_st;
  const double relaxp = eqp.relaxv;
  const double thetap = eqp.theta;

  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict b_face_surf = fvq->b_face_surf;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_nreal_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_rreal_3_t *restrict diipf = fvq->diipf;
  const cs_rreal_3_t *restrict djjpf = fvq->djjpf;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;

  const int *bc_type = cs_glob_bc_type;
  const int *restrict c_disable_flag = (fvq->has_disable_flag) ?
    fvq->c_disable_flag : nullptr;

  /* Parallel or device dispatch */
  cs_dispatch_context ctx;
  ctx.set_use_gpu(false);  /* steady case not ported to GPU */

  /* Local variables */

  const int f_id = (f != nullptr) ? f->id : -1;
  char var_name[64];

  int iupwin = 0;
  int w_stride = 1;

  /* Default value for cp */
  cs_real_t cpi = 1.0, cpj = 1.0;

  cs_real_t *coface = nullptr, *cofbce = nullptr;

  cs_real_3_t *grad;
  cs_real_3_t *gradup = nullptr;
  cs_real_3_t *gradst = nullptr;

  cs_real_t *df_limiter = nullptr;

  cs_real_t *gweight = nullptr;

  cs_real_t  *v_slope_test = cs_get_v_slope_test(f_id,  eqp);

  /* Internal coupling variables */
  cs_real_t *pvar_local = nullptr;
  cs_real_t *pvar_distant = nullptr;
  cs_real_t *df_limiter_local = nullptr;
  int coupling_id = -1;
  cs_lnum_t n_local = 0, n_distant = 0;
  const cs_lnum_t *faces_local = nullptr, *faces_distant = nullptr;
  cs_internal_coupling_t *cpl = nullptr;

  /* Initialization */

  /* Allocate work arrays */

  CS_MALLOC(grad, n_cells_ext, cs_real_3_t);

  /* Choose gradient type */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  cs_gradient_type_by_imrgra(eqp.d_gradient_r,
                             &gradient_type,
                             &halo_type);

  /* Handle cases where only the previous values (already synchronized)
     or current values are provided */

  if (pvar != nullptr)
    cs_sync_scalar_halo(m, halo_type, pvar);
  if (pvara == nullptr)
    pvara = (const cs_real_t *)pvar;

  const cs_real_t  *restrict _pvar = (pvar != nullptr) ? pvar : pvara;

  /* Limiters */

  if (f_id != -1) {
    f = cs_field_by_id(f_id);

    int df_limiter_id =
      cs_field_get_key_int(f, cs_field_key_id("diffusion_limiter_id"));
    if (df_limiter_id > -1)
      df_limiter = cs_field_by_id(df_limiter_id)->val;

    snprintf(var_name, 63, "%s", f->name);
  }
  else if (isstpp > 1) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of isstpp for a work array"));
  }
  else {
    strncpy(var_name, "[scalar convection-diffusion]", 63);
  }
  var_name[63] = '\0';

  if (iwarnp >= 2) {
    if (ischcp == 1) {
      bft_printf(
        _(" %s: Convection in centered blending with %f percent of upwind\n"),
        var_name, (1.-blencp)*100.);
    }
    else {
      bft_printf(
        _(" %s: Convection in 2nd order blending with %f percent of upwind\n"),
        var_name, (1.-blencp)*100.);
    }
  }

  iupwin = (blencp > 0.) ? 0 : 1;

  if (icoupl > 0) {
    assert(f_id != -1);
    const int coupling_key_id = cs_field_key_id("coupling_entity");
    coupling_id = cs_field_get_key_int(f, coupling_key_id);
    cpl = cs_internal_coupling_by_id(coupling_id);
    cs_internal_coupling_coupled_faces(cpl,
                                       &n_local,
                                       &faces_local,
                                       &n_distant,
                                       &faces_distant);
  }

  /* Compute the balance with reconstruction */

  /* Compute the gradient of the variable */

  /* The gradient (grad) is used in the flux reconstruction and the slope test.
     Thus we must compute it:
         - when we have diffusion and we reconstruct the fluxes,
         - when the convection scheme is the legacy SOLU,
         - when we have convection, we are not in pure upwind
           and we reconstruct the fluxes,
         - when we have convection, we are not in pure upwind
           and we have not shunted the slope test,
         - when we use NVD / TVD schemes.
  */

  bool is_ischcp
    = (   (xcpp != nullptr && ischcp == 0)
       || (xcpp == nullptr && (ischcp == 0 || ischcp == 3 || ischcp == 4)));

  if (   (idiffp != 0 && ircflp == 1)
      || (   iconvp != 0 && iupwin == 0
          && (ircflp == 1 || isstpp == 0 || is_ischcp))) {

    if (f_id != -1) {
      /* Get the calculation option from the field */
      if (f->type & CS_FIELD_VARIABLE && eqp.iwgrec == 1) {
        if (eqp.idiff > 0) {
          int key_id = cs_field_key_id("gradient_weighting_id");
          int diff_id = cs_field_get_key_int(f, key_id);
          if (diff_id > -1) {
            cs_field_t *weight_f = cs_field_by_id(diff_id);
            gweight = weight_f->val;
            w_stride = weight_f->dim;
            cs_field_synchronize(weight_f, halo_type);
          }
        }
      }
    }

    cs_gradient_scalar_synced_input(var_name,
                                    gradient_type,
                                    halo_type,
                                    inc,
                                    eqp.nswrgr,
                                    0, /* hyd_p_flag */
                                    w_stride,
                                    eqp.verbosity,
                                    (cs_gradient_limit_t)eqp.imligr,
                                    eqp.epsrgr,
                                    eqp.d_climgr,
                                    nullptr, /* f_ext exterior force */
                                    bc_coeffs,
                                    _pvar,
                                    gweight, /* Weighted gradient */
                                    cpl,
                                    grad);

  }
  else {

#   pragma omp parallel for
    for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
      grad[cell_id][0] = 0.;
      grad[cell_id][1] = 0.;
      grad[cell_id][2] = 0.;
    }
  }

  /* Compute gradients used in convection schemes */

  if (iconvp > 0 && iupwin == 0) {

    /* Compute cell gradient used in slope test */
    if (isstpp == 0) {

      CS_MALLOC(gradst, n_cells_ext, cs_real_3_t);

#     pragma omp parallel for
      for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
        gradst[cell_id][0] = 0.;
        gradst[cell_id][1] = 0.;
        gradst[cell_id][2] = 0.;
      }

      cs_slope_test_gradient(f_id,
                             ctx,
                             inc,
                             (const cs_real_3_t *)grad,
                             gradst,
                             _pvar,
                             bc_coeffs,
                             i_massflux);

    }

    /* Pure SOLU scheme */
    if (ischcp == 2 || (xcpp != nullptr && ischcp == 4)) {

      CS_MALLOC(gradup, n_cells_ext, cs_real_3_t);

#     pragma omp parallel for
      for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
        gradup[cell_id][0] = 0.;
        gradup[cell_id][1] = 0.;
        gradup[cell_id][2] = 0.;
      }

      cs_upwind_gradient(f_id,
                         ctx,
                         inc,
                         halo_type,
                         bc_coeffs,
                         i_massflux,
                         b_massflux,
                         _pvar,
                         gradup);

    }

  }

  /* ======================================================================
    ---> Contribution from interior faces
    ======================================================================*/

  cs_gnum_t n_upwind = 0;

  if (n_cells_ext>n_cells) {
#   pragma omp parallel for if(n_cells_ext - n_cells > CS_THR_MIN)
    for (cs_lnum_t cell_id = n_cells; cell_id < n_cells_ext; cell_id++) {
      rhs[cell_id] = 0.;
    }
  }

  /* --> Pure upwind flux
    =====================*/

  if (iupwin == 1) {

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for reduction(+:n_upwind)
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          if (xcpp != nullptr) {
            cpi = xcpp[ii];
            cpj = xcpp[jj];
          }

          /* in parallel, face will be counted by one and only one rank */
          if (ii < n_cells) {
            n_upwind++;
          }

          cs_real_2_t fluxij = {0.,0.};

          cs_real_t pifri, pjfri, pifrj, pjfrj;
          cs_real_t pip, pjp, pipr, pjpr;

          cs_real_t bldfrp = (cs_real_t) ircflp;
          /* Local limitation of the reconstruction */
          if (df_limiter != nullptr && ircflp > 0)
            bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]), 0.);

          cs_i_cd_steady_upwind(bldfrp,
                                relaxp,
                                diipf[face_id],
                                djjpf[face_id],
                                grad[ii],
                                grad[jj],
                                _pvar[ii],
                                _pvar[jj],
                                pvara[ii],
                                pvara[jj],
                                &pifri,
                                &pifrj,
                                &pjfri,
                                &pjfrj,
                                &pip,
                                &pjp,
                                &pipr,
                                &pjpr);

          cs_i_conv_flux(iconvp,
                         1.,
                         1,
                         _pvar[ii],
                         _pvar[jj],
                         pifri,
                         pifrj,
                         pjfri,
                         pjfrj,
                         i_massflux[face_id],
                         cpi,
                         cpj,
                         fluxij);

          cs_i_diff_flux(idiffp,
                         1.,
                         pip,
                         pjp,
                         pipr,
                         pjpr,
                         i_visc[face_id],
                         fluxij);

          rhs[ii] -= fluxij[0];
          rhs[jj] += fluxij[1];

        }
      }
    }

    /* --> Flux with no slope test or Min/Max Beta limiter
       ====================================================*/

  }
  else if (isstpp == 1 || isstpp == 2) {

    if (ischcp < 0 || ischcp > 4) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcv"));
    }

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          if (xcpp != nullptr) {
            cpi = xcpp[ii];
            cpj = xcpp[jj];
          }

          cs_real_2_t fluxij = {0.,0.};

          cs_real_t pifri, pjfri, pifrj, pjfrj;
          cs_real_t pip, pjp, pipr, pjpr;

          cs_real_t bldfrp = (cs_real_t) ircflp;
          /* Local limitation of the reconstruction */
          if (df_limiter != nullptr && ircflp > 0)
            bldfrp = cs::max(cs::min(df_limiter[ii],
                                     df_limiter[jj]), 0.);

          cs_i_cd_steady(bldfrp,
                         ischcp,
                         relaxp,
                         blencp,
                         weight[face_id],
                         cell_cen[ii],
                         cell_cen[jj],
                         i_face_cog[face_id],
                         diipf[face_id],
                         djjpf[face_id],
                         grad[ii],
                         grad[jj],
                         gradup[ii],
                         gradup[jj],
                         _pvar[ii],
                         _pvar[jj],
                         pvara[ii],
                         pvara[jj],
                         &pifri,
                         &pifrj,
                         &pjfri,
                         &pjfrj,
                         &pip,
                         &pjp,
                         &pipr,
                         &pjpr);

          cs_i_conv_flux(iconvp,
                         1.,
                         1,
                         _pvar[ii],
                         _pvar[jj],
                         pifri,
                         pifrj,
                         pjfri,
                         pjfrj,
                         i_massflux[face_id],
                         cpi,
                         cpj,
                         fluxij);

          cs_i_diff_flux(idiffp,
                         1.,
                         pip,
                         pjp,
                         pipr,
                         pjpr,
                         i_visc[face_id],
                         fluxij);

          rhs[ii] -= fluxij[0];
          rhs[jj] += fluxij[1];

        }
      }
    }

    /* --> Flux with slope test
       ============================================*/

  }
  else { /* isstpp = 0 */

    if (ischcp < 0 || ischcp > 2) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcv"));
    }

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for reduction(+:n_upwind)
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          if (xcpp != nullptr) {
            cpi = xcpp[ii];
            cpj = xcpp[jj];
          }

          cs_real_2_t fluxij = {0., 0.};

          bool upwind_switch = false;
          cs_real_t pifri, pjfri, pifrj, pjfrj;
          cs_real_t pip, pjp, pipr, pjpr;

          cs_real_t bldfrp = (cs_real_t) ircflp;
          /* Local limitation of the reconstruction */
          if (df_limiter != nullptr && ircflp > 0)
            bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                             0.);

          cs_i_cd_steady_slope_test(&upwind_switch,
                                    iconvp,
                                    bldfrp,
                                    ischcp,
                                    relaxp,
                                    blencp,
                                    blend_st,
                                    weight[face_id],
                                    i_dist[face_id],
                                    cell_cen[ii],
                                    cell_cen[jj],
                                    i_face_u_normal[face_id],
                                    i_face_cog[face_id],
                                    diipf[face_id],
                                    djjpf[face_id],
                                    i_massflux[face_id],
                                    grad[ii],
                                    grad[jj],
                                    gradup[ii],
                                    gradup[jj],
                                    gradst[ii],
                                    gradst[jj],
                                    _pvar[ii],
                                    _pvar[jj],
                                    pvara[ii],
                                    pvara[jj],
                                    &pifri,
                                    &pifrj,
                                    &pjfri,
                                    &pjfrj,
                                    &pip,
                                    &pjp,
                                    &pipr,
                                    &pjpr);

          cs_i_conv_flux(iconvp,
                         1.,
                         1,
                         _pvar[ii],
                         _pvar[jj],
                         pifri,
                         pifrj,
                         pjfri,
                         pjfrj,
                         i_massflux[face_id],
                         cpi,
                         cpj,
                         fluxij);

          cs_i_diff_flux(idiffp,
                         1.,
                         pip,
                         pjp,
                         pipr,
                         pjpr,
                         i_visc[face_id],
                         fluxij);

          if (upwind_switch) {

            /* in parallel, face will be counted by one and only one rank */
            if (ii < n_cells)
              n_upwind++;
            if (v_slope_test != nullptr) {
              cs_real_t q_d_vol_ii
                =   std::abs(i_massflux[face_id])
                  * cs_mq_cell_vol_inv(ii, c_disable_flag, cell_vol);
              cs_real_t q_d_vol_jj
                =   std::abs(i_massflux[face_id])
                  * cs_mq_cell_vol_inv(jj, c_disable_flag, cell_vol);

              v_slope_test[ii] += q_d_vol_ii;
              v_slope_test[jj] += q_d_vol_jj;
            }

          }

          rhs[ii] -= fluxij[0];
          rhs[jj] += fluxij[1];

        }
      }
    }

  } /* End test on pure upwind */

  if (iwarnp >= 2 && iconvp == 1) {

    /* Sum number of clippings */
    cs_parall_counter(&n_upwind, 1);

    bft_printf(_(" %s: %llu Faces with upwind on %llu interior faces\n"),
               var_name, (unsigned long long)n_upwind,
               (unsigned long long)m->n_g_i_c_faces);
  }

  /* ======================================================================
    ---> Contribution from boundary faces
    ======================================================================*/

  /* Boundary convective flux are all computed with an upwind scheme */
  if (icvflb == false || xcpp != nullptr) {

#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
    for (int t_id = 0; t_id < n_b_threads; t_id++) {
      for (cs_lnum_t face_id = b_group_index[t_id*2];
           face_id < b_group_index[t_id*2 + 1];
           face_id++) {

        cs_lnum_t ii = b_face_cells[face_id];

        if (xcpp != nullptr)
          cpi = xcpp[ii];

        cs_real_t fluxi = 0.;
        cs_real_t pir, pipr;

        cs_real_t bldfrp = (cs_real_t) ircflb;
        /* Local limitation of the reconstruction */
        if (df_limiter != nullptr && ircflb > 0)
          bldfrp = cs::max(df_limiter[ii], 0.);

        cs_b_cd_steady(bldfrp,
                       relaxp,
                       diipb[face_id],
                       grad[ii],
                       _pvar[ii],
                       pvara[ii],
                       &pir,
                       &pipr);

        cs_b_upwind_flux(iconvp,
                         1.,
                         1,
                         inc,
                         bc_type[face_id],
                         _pvar[ii],
                         pir,
                         pipr,
                         coefap[face_id],
                         coefbp[face_id],
                         b_massflux[face_id],
                         cpi,
                         &fluxi);

        cs_b_diff_flux(idiffp,
                       1., /* thetap */
                       inc,
                       pipr,
                       cofafp[face_id],
                       cofbfp[face_id],
                       b_visc[face_id],
                       &fluxi);

        rhs[ii] -= fluxi;

      }
    }

    /* The scalar is internal_coupled and an implicit contribution
     * is required */
    if (icoupl > 0) {
      /* Prepare data for sending */
      CS_MALLOC(pvar_distant, n_distant, cs_real_t);

      for (cs_lnum_t ii = 0; ii < n_distant; ii++) {
        cs_lnum_t face_id = faces_distant[ii];
        cs_lnum_t jj = b_face_cells[face_id];
        cs_real_t pip, pipr;

        cs_real_t bldfrp = (cs_real_t) ircflb;
        /* Local limitation of the reconstruction */
        if (df_limiter != nullptr && ircflb > 0)
          bldfrp = cs::max(df_limiter[jj], 0.);

        cs_b_cd_steady(bldfrp,
                       relaxp,
                       diipb[face_id],
                       grad[jj],
                       _pvar[jj],
                       pvara[jj],
                       &pip,
                       &pipr);
        pvar_distant[ii] = pipr;
      }

      /* Receive data */
      CS_MALLOC(pvar_local, n_local, cs_real_t);
      cs_internal_coupling_exchange_var(cpl,
                                        1, /* Dimension */
                                        pvar_distant,
                                        pvar_local);

      /* Exchange diffusion limiter */
      if (df_limiter != nullptr) {
        CS_MALLOC(df_limiter_local, n_local, cs_real_t);
        cs_internal_coupling_exchange_var(cpl,
                                          1, /* Dimension */
                                          df_limiter,
                                          df_limiter_local);
      }

      /* Flux contribution */
      assert(f != nullptr);
      cs_real_t *hintp = f->bc_coeffs->hint;
      cs_real_t *hextp = f->bc_coeffs->rcodcl2;
      for (cs_lnum_t ii = 0; ii < n_local; ii++) {
        cs_lnum_t face_id = faces_local[ii];
        cs_lnum_t jj = b_face_cells[face_id];
        cs_real_t surf = b_face_surf[face_id];
        cs_real_t pip, pipr, pjpr;
        cs_real_t fluxi = 0.;

        cs_real_t bldfrp = (cs_real_t) ircflb;
        /* Local limitation of the reconstruction */
        if (df_limiter != nullptr && ircflb > 0)
          bldfrp = cs::max(cs::min(df_limiter_local[ii],
                                   df_limiter[jj]),
                           0.);

        cs_b_cd_steady(bldfrp,
                       relaxp,
                       diipb[face_id],
                       grad[jj],
                       _pvar[jj],
                       pvara[jj],
                       &pip,
                       &pipr);

        pjpr = pvar_local[ii];

        cs_real_t hint = hintp[face_id];
        cs_real_t hext = hextp[face_id];
        cs_real_t heq = _calc_heq(hint, hext)*surf;

        cs_b_diff_flux_coupling(idiffp,
                                pipr,
                                pjpr,
                                heq,
                                &fluxi);

        rhs[jj] -= thetap * fluxi;
      }

      CS_FREE(pvar_local);
      /* Sending structures are no longer needed */
      CS_FREE(pvar_distant);
      if (df_limiter != nullptr)
        CS_FREE(df_limiter_local);
    }

    /* Boundary convective flux is imposed at some faces
       (tagged in icvfli array) */
  }
  else if (icvflb && xcpp == nullptr) {

    /* Retrieve the value of the convective flux to be imposed */
    if (f_id != -1) {
      coface = f->bc_coeffs->ac;
      cofbce = f->bc_coeffs->bc;
    }
    else {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of icvflb and f_id"));
    }

#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
    for (int t_id = 0; t_id < n_b_threads; t_id++) {
      for (cs_lnum_t face_id = b_group_index[t_id*2];
           face_id < b_group_index[t_id*2 + 1];
           face_id++) {

        cs_lnum_t ii = b_face_cells[face_id];

        cs_real_t fluxi = 0.;
        cs_real_t pir, pipr;

        cs_real_t bldfrp = (cs_real_t) ircflb;
        /* Local limitation of the reconstruction */
        if (df_limiter != nullptr && ircflb > 0)
          bldfrp = cs::max(df_limiter[ii], 0.);

        cs_b_cd_steady(bldfrp,
                       relaxp,
                       diipb[face_id],
                       grad[ii],
                       _pvar[ii],
                       pvara[ii],
                       &pir,
                       &pipr);

        cs_b_imposed_conv_flux(iconvp,
                               1.,
                               1,
                               inc,
                               bc_type[face_id],
                               icvfli[face_id],
                               _pvar[ii],
                               pir,
                               pipr,
                               coefap[face_id],
                               coefbp[face_id],
                               coface[face_id],
                               cofbce[face_id],
                               b_massflux[face_id],
                               1., /* xcpp */
                               &fluxi);

        cs_b_diff_flux(idiffp,
                       1., /* thetap */
                       inc,
                       pipr,
                       cofafp[face_id],
                       cofbfp[face_id],
                       b_visc[face_id],
                       &fluxi);

        rhs[ii] -= fluxi;

      }
    }
  }

  /* Free memory */
  CS_FREE(grad);
  CS_FREE(gradup);
  CS_FREE(gradst);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update face flux with convection contribution of a standard transport
 * equation of a scalar field \f$ \varia \f$.
 *
 * <a name="cs_face_convection_scalar"></a>
 *
 * \f[
 * C_\ij = \dot{m}_\ij \left( \varia_\fij - \varia_\celli \right)
 * \f]
 *
 * \param[in]     f_id          pointer to field id, or nullptr
 * \param[in]     eqp           equation parameters
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imasac        take mass accumulation into account?
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in]     bc_coeffs     boundary condition structure for the variable
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in,out] i_conv_flux   scalar convection flux at interior faces
 * \param[in,out] b_conv_flux   scalar convection flux at boundary faces
 */
/*----------------------------------------------------------------------------*/

static void
_face_convection_scalar_steady(const cs_field_t           *f,
                               const cs_equation_param_t   eqp,
                               int                         icvflb,
                               int                         inc,
                               cs_real_t         *restrict pvar,
                               const cs_real_t   *restrict pvara,
                               const int                   icvfli[],
                               const cs_field_bc_coeffs_t *bc_coeffs,
                               const cs_real_t             i_massflux[],
                               const cs_real_t             b_massflux[],
                               cs_real_t                   i_conv_flux[][2],
                               cs_real_t                   b_conv_flux[])
{
  cs_real_t *coefap = bc_coeffs->a;
  cs_real_t *coefbp = bc_coeffs->b;

  const int iconvp = eqp.iconv;
  const int nswrgp = eqp.nswrgr;
  const int ircflp = eqp.ircflu;
  const int ircflb = (ircflp > 0) ? eqp.b_diff_flux_rc : 0;
  const int ischcp = eqp.ischcv;
  const int isstpp = eqp.isstpc;
  const int iwarnp = eqp.verbosity;
  const int icoupl = eqp.icoupl;
  const double blencp = eqp.blencv;
  const double blend_st = eqp.blend_st;
  const double relaxp = eqp.relaxv;

  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_nreal_3_t *restrict i_face_u_normal  = fvq->i_face_u_normal;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_rreal_3_t *restrict diipf = fvq->diipf;
  const cs_rreal_3_t *restrict djjpf = fvq->djjpf;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;

  const int *bc_type = cs_glob_bc_type;
  const int *restrict c_disable_flag = (fvq->has_disable_flag) ?
    fvq->c_disable_flag : nullptr;

  /* Local variables */

  char var_name[64];
  const int f_id = (f != nullptr) ? f->id : -1;

  /* Parallel or device dispatch */
  cs_dispatch_context ctx;
  ctx.set_use_gpu(false);  /* steady case not ported to GPU */

  int w_stride = 1;

  cs_real_t *coface = nullptr, *cofbce = nullptr;

  cs_real_3_t *grad;
  cs_real_3_t *gradup = nullptr;
  cs_real_3_t *gradst = nullptr;

  cs_real_t *local_min = nullptr;
  cs_real_t *local_max = nullptr;

  cs_real_t *df_limiter = nullptr;
  cs_real_t *gweight = nullptr;

  cs_real_t  *v_slope_test = cs_get_v_slope_test(f_id,  eqp);

  /* Internal coupling variables */
  int coupling_id;
  cs_internal_coupling_t *cpl = nullptr;

  /* Allocate work arrays */

  CS_MALLOC(grad, n_cells_ext, cs_real_3_t);

  /* Choose gradient type */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  cs_gradient_type_by_imrgra(eqp.d_gradient_r,
                             &gradient_type,
                             &halo_type);

  /* Handle cases where only the previous values (already synchronized)
     or current values are provided */

  if (pvar != nullptr)
    cs_sync_scalar_halo(m, halo_type, pvar);
  if (pvara == nullptr)
    pvara = (const cs_real_t *)pvar;

  const cs_real_t  *restrict _pvar = (pvar != nullptr) ? pvar : pvara;

  /* Slope limiters */

  if (f_id != -1) {
    f = cs_field_by_id(f_id);

    /* NVD/TVD limiters */
    if (ischcp == 4) {
      CS_MALLOC(local_max, n_cells_ext, cs_real_t);
      CS_MALLOC(local_min, n_cells_ext, cs_real_t);
      cs_field_local_extrema_scalar(f_id,
                                    halo_type,
                                    local_max,
                                    local_min);
    }

    int df_limiter_id =
      cs_field_get_key_int(f, cs_field_key_id("diffusion_limiter_id"));
    if (df_limiter_id > -1)
      df_limiter = cs_field_by_id(df_limiter_id)->val;

    snprintf(var_name, 63, "%s", f->name);
  }
  else if (isstpp > 1) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of isstpp for a work array"));
  }
  else {
    strncpy(var_name, "[scalar face flux from convection]", 63);
  }
  var_name[63] = '\0';

  if (iwarnp >= 2) {
    if (ischcp == 1) {
      bft_printf(
        _(" %s: Convection in centered blending with %f percent of upwind\n"),
        var_name, (1.-blencp)*100.);
    }
    else {
      bft_printf(
        _(" %s: Convection in 2nd order blending with %f percent of upwind\n"),
        var_name, (1.-blencp)*100.);
    }
  }

  bool pure_upwind = (blencp > 0.) ? false : true;

  if (icoupl > 0) {
    assert(f_id != -1);
    const int coupling_key_id = cs_field_key_id("coupling_entity");
    coupling_id = cs_field_get_key_int(f, coupling_key_id);
    cpl = cs_internal_coupling_by_id(coupling_id);
  }

  /* Compute the balance with reconstruction */

  /* Compute the gradient of the variable */

  /* The gradient (grad) is used in the flux reconstruction and the slope test.
     Thus we must compute it:
         - when the convection scheme is the legacy SOLU,
         - when we have convection, we are not in pure upwind
           and we reconstruct the fluxes,
         - when we have convection, we are not in pure upwind
           and we have not shunted the slope test,
         - when we use NVD / TVD schemes.
  */

  if (   iconvp != 0 && pure_upwind == false
      && (ischcp == 0 || ircflp == 1 || isstpp == 0 || ischcp == 4)) {

    if (f_id != -1) {
      /* Get the calculation option from the field */
      if (f->type & CS_FIELD_VARIABLE && eqp.iwgrec == 1) {
        if (eqp.idiff > 0) {
          int key_id = cs_field_key_id("gradient_weighting_id");
          int diff_id = cs_field_get_key_int(f, key_id);
          if (diff_id > -1) {
            cs_field_t *weight_f = cs_field_by_id(diff_id);
            gweight = weight_f->val;
            w_stride = weight_f->dim;
            cs_field_synchronize(weight_f, halo_type);
          }
        }
      }
    }

    cs_gradient_scalar_synced_input(var_name,
                                    gradient_type,
                                    halo_type,
                                    inc,
                                    nswrgp,
                                    0, /* hyd_p_flag */
                                    w_stride,
                                    eqp.verbosity,
                                    (cs_gradient_limit_t)(eqp.imligr),
                                    eqp.epsrgr,
                                    eqp.climgr,
                                    nullptr, /* f_ext exterior force */
                                    bc_coeffs,
                                    _pvar,
                                    gweight, /* Weighted gradient */
                                    cpl,
                                    grad);

  }
  else {

#   pragma omp parallel for
    for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
      grad[cell_id][0] = 0.;
      grad[cell_id][1] = 0.;
      grad[cell_id][2] = 0.;
    }
  }

  /* Compute gradients used in convection schemes */

  if (iconvp > 0 && pure_upwind == false) {

    /* Compute cell gradient used in slope test */
    if (isstpp == 0) {

      CS_MALLOC(gradst, n_cells_ext, cs_real_3_t);

#     pragma omp parallel for
      for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
        gradst[cell_id][0] = 0.;
        gradst[cell_id][1] = 0.;
        gradst[cell_id][2] = 0.;
      }

      cs_slope_test_gradient(f_id,
                             ctx,
                             inc,
                             (const cs_real_3_t *)grad,
                             gradst,
                             _pvar,
                             bc_coeffs,
                             i_massflux);

    }

    /* Pure SOLU scheme */
    if (ischcp == 2) {

      CS_MALLOC(gradup, n_cells_ext, cs_real_3_t);

#     pragma omp parallel for
      for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++) {
        gradup[cell_id][0] = 0.;
        gradup[cell_id][1] = 0.;
        gradup[cell_id][2] = 0.;
      }

      cs_upwind_gradient(f_id,
                         ctx,
                         inc,
                         halo_type,
                         bc_coeffs,
                         i_massflux,
                         b_massflux,
                         _pvar,
                         gradup);

    }

  }

  /* ======================================================================
    ---> Contribution from interior faces
    ======================================================================*/

  cs_gnum_t n_upwind = 0;

  /* --> Pure upwind flux
    =====================*/

  if (pure_upwind) {

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for reduction(+:n_upwind)
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          /* in parallel, face will be counted by one and only one rank */
          if (ii < n_cells) {
            n_upwind++;
          }

          cs_real_t pifri, pjfri, pifrj, pjfrj;
          cs_real_t pip, pjp, pipr, pjpr;

          cs_real_t bldfrp = (cs_real_t) ircflp;
          /* Local limitation of the reconstruction */
          if (df_limiter != nullptr && ircflp > 0)
            bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]), 0.);

          cs_i_cd_steady_upwind(bldfrp,
                                relaxp,
                                diipf[face_id],
                                djjpf[face_id],
                                grad[ii],
                                grad[jj],
                                _pvar[ii],
                                _pvar[jj],
                                pvara[ii],
                                pvara[jj],
                                &pifri,
                                &pifrj,
                                &pjfri,
                                &pjfrj,
                                &pip,
                                &pjp,
                                &pipr,
                                &pjpr);

          cs_i_conv_flux(iconvp,
                         1.,
                         1,
                         _pvar[ii],
                         _pvar[jj],
                         pifri,
                         pifrj,
                         pjfri,
                         pjfrj,
                         i_massflux[face_id],
                         1., /* xcpp */
                         1., /* xcpp */
                         i_conv_flux[face_id]);
        }
      }
    }
  }

  /* --> Flux with no slope test or Min/Max Beta limiter
    ====================================================*/

  else if (isstpp == 1 || isstpp == 2) {

    if (ischcp < 0 || ischcp > 4) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcv"));
    }

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

            cs_real_t pifri, pjfri, pifrj, pjfrj;
            cs_real_t pip, pjp, pipr, pjpr;

            cs_real_t bldfrp = (cs_real_t) ircflp;
            /* Local limitation of the reconstruction */
            if (df_limiter != nullptr && ircflp > 0)
              bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                               0.);

            cs_i_cd_steady(bldfrp,
                           ischcp,
                           relaxp,
                           blencp,
                           weight[face_id],
                           cell_cen[ii],
                           cell_cen[jj],
                           i_face_cog[face_id],
                           diipf[face_id],
                           djjpf[face_id],
                           grad[ii],
                           grad[jj],
                           gradup[ii],
                           gradup[jj],
                           _pvar[ii],
                           _pvar[jj],
                           pvara[ii],
                           pvara[jj],
                           &pifri,
                           &pifrj,
                           &pjfri,
                           &pjfrj,
                           &pip,
                           &pjp,
                           &pipr,
                           &pjpr);

            cs_i_conv_flux(iconvp,
                           1.,
                           1,
                           _pvar[ii],
                           _pvar[jj],
                           pifri,
                           pifrj,
                           pjfri,
                           pjfrj,
                           i_massflux[face_id],
                           1., /* xcpp */
                           1., /* xcpp */
                           i_conv_flux[face_id]);

        }
      }
    }

    /* --> Flux with slope test or NVD/TVD limiter
    ============================================*/

  }
  else { /* isstpp = 0 */

    if (ischcp < 0 || ischcp > 2) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcv"));
    }

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#       pragma omp parallel for reduction(+:n_upwind)
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          bool upwind_switch = false;
          cs_real_t pifri, pjfri, pifrj, pjfrj;
          cs_real_t pip, pjp, pipr, pjpr;

          cs_real_t bldfrp = (cs_real_t) ircflp;
          /* Local limitation of the reconstruction */
          if (df_limiter != nullptr && ircflp > 0)
            bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                             0.);

          cs_i_cd_steady_slope_test(&upwind_switch,
                                    iconvp,
                                    bldfrp,
                                    ischcp,
                                    relaxp,
                                    blencp,
                                    blend_st,
                                    weight[face_id],
                                    i_dist[face_id],
                                    cell_cen[ii],
                                    cell_cen[jj],
                                    i_face_u_normal[face_id],
                                    i_face_cog[face_id],
                                    diipf[face_id],
                                    djjpf[face_id],
                                    i_massflux[face_id],
                                    grad[ii],
                                    grad[jj],
                                    gradup[ii],
                                    gradup[jj],
                                    gradst[ii],
                                    gradst[jj],
                                    _pvar[ii],
                                    _pvar[jj],
                                    pvara[ii],
                                    pvara[jj],
                                    &pifri,
                                    &pifrj,
                                    &pjfri,
                                    &pjfrj,
                                    &pip,
                                    &pjp,
                                    &pipr,
                                    &pjpr);

          cs_i_conv_flux(iconvp,
                         1.,
                         1,
                         _pvar[ii],
                         _pvar[jj],
                         pifri,
                         pifrj,
                         pjfri,
                         pjfrj,
                         i_massflux[face_id],
                         1., /* xcpp */
                         1., /* xcpp */
                         i_conv_flux[face_id]);

          if (upwind_switch) {

            /* in parallel, face will be counted by one and only one rank */
            if (ii < n_cells)
              n_upwind++;
            if (v_slope_test != nullptr) {
              cs_real_t q_d_vol_ii
                =   std::abs(i_massflux[face_id])
                  * cs_mq_cell_vol_inv(ii, c_disable_flag, cell_vol);
              cs_real_t q_d_vol_jj
                =   std::abs(i_massflux[face_id])
                  * cs_mq_cell_vol_inv(jj, c_disable_flag, cell_vol);

              v_slope_test[ii] += q_d_vol_ii;
              v_slope_test[jj] += q_d_vol_jj;
            }

          }

        }
      }
    }

  } /* End pure_upwind */

  if (iwarnp >= 2 && iconvp == 1) {

    /* Sum number of clippings */
    cs_parall_counter(&n_upwind, 1);

    bft_printf(_(" %s: %llu Faces with upwind on %llu interior faces\n"),
               var_name, (unsigned long long)n_upwind,
               (unsigned long long)m->n_g_i_c_faces);
  }

  /* ======================================================================
    ---> Contribution from boundary faces
    ======================================================================*/

  /* Boundary convective flux are all computed with an upwind scheme */
  if (icvflb == 0) {

#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
    for (int t_id = 0; t_id < n_b_threads; t_id++) {
      for (cs_lnum_t face_id = b_group_index[t_id*2];
           face_id < b_group_index[t_id*2 + 1];
           face_id++) {

        cs_lnum_t ii = b_face_cells[face_id];

        cs_real_t pir, pipr;

        cs_real_t bldfrp = (cs_real_t) ircflb;
        /* Local limitation of the reconstruction */
        if (df_limiter != nullptr && ircflb > 0)
          bldfrp = cs::max(df_limiter[ii], 0.);

        cs_b_cd_steady(bldfrp,
                       relaxp,
                       diipb[face_id],
                       grad[ii],
                       _pvar[ii],
                       pvara[ii],
                       &pir,
                       &pipr);

        cs_b_upwind_flux(iconvp,
                         1.,
                         1,
                         inc,
                         bc_type[face_id],
                         _pvar[ii],
                         pir,
                         pipr,
                         coefap[face_id],
                         coefbp[face_id],
                         b_massflux[face_id],
                         1., /* xcpp */
                         &(b_conv_flux[face_id]));

      }
    }

  /* Boundary convective flux is imposed at some faces
     (tagged in icvfli array) */
  }
  else if (icvflb == 1) {

    /* Retrieve the value of the convective flux to be imposed */
    if (f_id != -1) {
      coface = f->bc_coeffs->ac;
      cofbce = f->bc_coeffs->bc;
    }
    else {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of icvflb and f_id"));
    }

#     pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
    for (int t_id = 0; t_id < n_b_threads; t_id++) {
      for (cs_lnum_t face_id = b_group_index[t_id*2];
           face_id < b_group_index[t_id*2 + 1];
           face_id++) {

        cs_lnum_t ii = b_face_cells[face_id];

        cs_real_t pir, pipr;

        cs_real_t bldfrp = (cs_real_t) ircflb;
        /* Local limitation of the reconstruction */
        if (df_limiter != nullptr && ircflb > 0)
          bldfrp = cs::max(df_limiter[ii], 0.);

        cs_b_cd_steady(bldfrp,
                       relaxp,
                       diipb[face_id],
                       grad[ii],
                       _pvar[ii],
                       pvara[ii],
                       &pir,
                       &pipr);

        cs_b_imposed_conv_flux(iconvp,
                               1.,
                               1,
                               inc,
                               bc_type[face_id],
                               icvfli[face_id],
                               _pvar[ii],
                               pir,
                               pipr,
                               coefap[face_id],
                               coefbp[face_id],
                               coface[face_id],
                               cofbce[face_id],
                               b_massflux[face_id],
                               1., /* xcpp */
                               &(b_conv_flux[face_id]));

      }
    }
  }

  /* Free memory */
  CS_FREE(grad);
  CS_FREE(gradup);
  CS_FREE(gradst);
  CS_FREE(local_max);
  CS_FREE(local_min);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the convection/diffusion terms of a
 * standard transport equation of a scalar field \f$ \varia \f$.
 *
 * More precisely, the right hand side \f$ Rhs \f$ is updated as
 * follows:
 * \f[
 * Rhs = Rhs - \sum_{\fij \in \Facei{\celli}}      \left(
 *        \dot{m}_\ij \left( \varia_\fij - \varia_\celli \right)
 *      - \mu_\fij \gradv_\fij \varia \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ Rhs \f$ has already been initialized before calling bilsc2!
 * - mind the sign minus
 *
 * Please refer to the
 * <a href="../../theory.pdf#bilsc2"><b>bilsc2</b></a> section of the
 * theory guide for more informations.
 *
 * \param[in]     f             pointer to field
 * \param[in]     eqp           equation parameters
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imasac        take mass accumulation into account?
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in]     bc_coeffs     boundary condition structure for the variable
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     xcpp          array of specific heat (\f$ C_p \f$), or nullptr
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 * \param[in,out] i_flux        interior flux (or nullptr)
 * \param[in,out] b_flux        boundary flux (or nullptr)
 */
/*----------------------------------------------------------------------------*/

template <bool is_thermal, bool store_flux>
static void
_convection_diffusion_scalar_unsteady
  (const cs_field_t           *f,
   const cs_equation_param_t  &eqp,
   bool                        icvflb,
   int                         inc,
   int                         imasac,
   cs_real_t         *restrict pvar,
   const cs_real_t   *restrict pvara,
   const int                   icvfli[],
   const cs_field_bc_coeffs_t *bc_coeffs,
   const cs_real_t             i_massflux[],
   const cs_real_t             b_massflux[],
   const cs_real_t             i_visc[],
   const cs_real_t             b_visc[],
   const cs_real_t             xcpp[],
   cs_real_t         *restrict rhs,
   cs_real_2_t       *restrict i_flux,
   cs_real_t         *restrict b_flux)
{
  const cs_real_t *coefap = bc_coeffs->a;
  const cs_real_t *coefbp = bc_coeffs->b;
  const cs_real_t *cofafp = bc_coeffs->af;
  const cs_real_t *cofbfp = bc_coeffs->bf;

  const int iconvp = eqp.iconv;
  const int idiffp = eqp.idiff;
  const int ircflp = eqp.ircflu;
  const int ircflb = (ircflp > 0) ? eqp.b_diff_flux_rc : 0;
  const int ischcp = eqp.ischcv;
  const int isstpp = eqp.isstpc;
  const int iwarnp = eqp.verbosity;
  const int icoupl = eqp.icoupl;
  cs_nvd_type_t limiter_choice = CS_NVD_N_TYPES;
  const double blencp = eqp.blencv;
  const double blend_st = eqp.blend_st;
  const double thetap = eqp.theta;

  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict b_face_surf = fvq->b_face_surf;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_nreal_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_rreal_3_t *restrict diipf = fvq->diipf;
  const cs_rreal_3_t *restrict djjpf = fvq->djjpf;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;

  const int *bc_type = cs_glob_bc_type;
  const int *restrict c_disable_flag = (fvq->has_disable_flag) ?
    fvq->c_disable_flag : nullptr;

  /* Parallel or device dispatch */
  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  /* Local variables */

  char var_name[64];
  const int f_id = (f != nullptr) ? f->id : -1;

  int w_stride = 1;

  bool pure_upwind = (blencp > 0.) ? false : true;
  cs_real_t *coface = nullptr, *cofbce = nullptr;

  cs_real_3_t *grad;
  cs_real_3_t *gradup = nullptr;
  cs_real_3_t *gradst = nullptr;

  cs_real_t *local_min = nullptr;
  cs_real_t *local_max = nullptr;
  cs_real_t *courant = nullptr;

  cs_real_t *cv_limiter = nullptr;
  cs_real_t *df_limiter = nullptr;

  cs_real_t *gweight = nullptr;

  const int key_lim_choice = cs_field_key_id("limiter_choice");

  cs_real_t  *v_slope_test = cs_get_v_slope_test(f_id,  eqp);

  /* Internal coupling variables */
  cs_real_t *pvar_local = nullptr;
  cs_real_t *pvar_distant = nullptr;
  cs_real_t *df_limiter_local = nullptr;
  int coupling_id = -1;
  cs_lnum_t n_local = 0, n_distant = 0;
  const cs_lnum_t *faces_local = nullptr, *faces_distant = nullptr;
  cs_internal_coupling_t *cpl = nullptr;

  /* Initialization */

  /* Allocate work arrays */

  CS_MALLOC_HD(grad, n_cells_ext, cs_real_3_t, cs_alloc_mode);

  /* Choose gradient type */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  cs_gradient_type_by_imrgra(eqp.d_gradient_r,
                             &gradient_type,
                             &halo_type);

  /* Handle cases where only the previous values (already synchronized)
     or current values are provided */

  if (pvar != nullptr)
    cs_sync_scalar_halo(m, halo_type, pvar);
  if (pvara == nullptr)
    pvara = (const cs_real_t *)pvar;

  const cs_real_t  *restrict _pvar = (pvar != nullptr) ? pvar : pvara;

  /* Limiters */

  if (f != nullptr) {

    /* NVD/TVD limiters */
    if (ischcp == 4) {
      limiter_choice = (cs_nvd_type_t)(cs_field_get_key_int(f, key_lim_choice));
      CS_MALLOC_HD(local_max, n_cells_ext, cs_real_t, cs_alloc_mode);
      CS_MALLOC_HD(local_min, n_cells_ext, cs_real_t, cs_alloc_mode);
      cs_field_local_extrema_scalar(f_id,
                                    halo_type,
                                    local_max,
                                    local_min);
      if (is_thermal) {
        assert(limiter_choice < CS_NVD_VOF_HRIC); /* VOF scheme forbidden here */
      }
      else if (limiter_choice >= CS_NVD_VOF_HRIC) {
        CS_MALLOC_HD(courant, n_cells_ext, cs_real_t, cs_alloc_mode);
        cs_cell_courant_number(f, ctx, courant);
      }
    }

    int cv_limiter_id =
      cs_field_get_key_int(f, cs_field_key_id("convection_limiter_id"));
    if (cv_limiter_id > -1)
      cv_limiter = cs_field_by_id(cv_limiter_id)->val;

    int df_limiter_id =
      cs_field_get_key_int(f, cs_field_key_id("diffusion_limiter_id"));
    if (df_limiter_id > -1)
      df_limiter = cs_field_by_id(df_limiter_id)->val;

    snprintf(var_name, 63, "%s", f->name);
  }
  else if (isstpp > 1) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of isstpp for a work array"));
  }
  else {
    strncpy(var_name, "[scalar convection-diffusion]", 63);
  }
  var_name[63] = '\0';

  if (iwarnp >= 2) {
    if (ischcp == 1) {
      bft_printf(
        _(" %s: Convection in centered blending with %f percent of upwind\n"),
        var_name, (1.-blencp)*100.);
    }
    else {
      bft_printf(
        _(" %s: Convection in 2nd order blending with %f percent of upwind\n"),
        var_name, (1.-blencp)*100.);
    }
  }

  if (icoupl > 0) {
    assert(f_id != -1);
    const int coupling_key_id = cs_field_key_id("coupling_entity");
    coupling_id = cs_field_get_key_int(f, coupling_key_id);
    cpl = cs_internal_coupling_by_id(coupling_id);
    cs_internal_coupling_coupled_faces(cpl,
                                       &n_local,
                                       &faces_local,
                                       &n_distant,
                                       &faces_distant);
  }

  /* Compute the balance with reconstruction */

  /* Compute the gradient of the variable */

  /* The gradient (grad) is used in the flux reconstruction and the slope test.
     Thus we must compute it:
         - when we have diffusion and we reconstruct the fluxes,
         - when the convection scheme is the legacy SOLU,
         - when we have convection, we are not in pure upwind
           and we reconstruct the fluxes,
         - when we have convection, we are not in pure upwind
           and we have not shunted the slope test,
         - when we use NVD / TVD schemes.
  */

  bool is_ischcp
    = (   (is_thermal && ischcp == 0)
       || (is_thermal == false && (ischcp == 0 || ischcp == 3 || ischcp == 4)));

  if (   (idiffp != 0 && ircflp == 1)
      || (   iconvp != 0 && pure_upwind == false
          && (ircflp == 1 || isstpp == 0 || is_ischcp))) {

    if (f != nullptr) {
      /* Get the calculation option from the field */
      if (f->type & CS_FIELD_VARIABLE && eqp.iwgrec == 1) {
        if (eqp.idiff > 0) {
          int key_id = cs_field_key_id("gradient_weighting_id");
          int diff_id = cs_field_get_key_int(f, key_id);
          if (diff_id > -1) {
            cs_field_t *weight_f = cs_field_by_id(diff_id);
            gweight = weight_f->val;
            w_stride = weight_f->dim;
            cs_field_synchronize(weight_f, halo_type);
          }
        }
      }
    }

    cs_gradient_scalar_synced_input(var_name,
                                    gradient_type,
                                    halo_type,
                                    inc,
                                    eqp.nswrgr,
                                    0, /* hyd_p_flag */
                                    w_stride,
                                    eqp.verbosity,
                                    (cs_gradient_limit_t)(eqp.imligr),
                                    eqp.epsrgr,
                                    eqp.d_climgr,
                                    nullptr, /* f_ext exterior force */
                                    bc_coeffs,
                                    _pvar,
                                    gweight, /* Weighted gradient */
                                    cpl,
                                    grad);

  }
  else {
    ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      grad[cell_id][0] = 0.;
      grad[cell_id][1] = 0.;
      grad[cell_id][2] = 0.;
    });
    ctx.wait();
  }

  /* Compute gradients used in convection schemes */

  if (iconvp > 0 && pure_upwind == false) {

    /* Compute cell gradient used in slope test */
    if (isstpp == 0) {

      CS_MALLOC_HD(gradst, n_cells_ext, cs_real_3_t, cs_alloc_mode);

      ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        gradst[cell_id][0] = 0.;
        gradst[cell_id][1] = 0.;
        gradst[cell_id][2] = 0.;
      });

      ctx.wait();

      cs_slope_test_gradient(f_id,
                             ctx,
                             inc,
                             (const cs_real_3_t *)grad,
                             gradst,
                             _pvar,
                             bc_coeffs,
                             i_massflux);

    }

    /* Pure SOLU scheme */
    if (ischcp == 2 || (is_thermal && ischcp == 4)) {

      CS_MALLOC_HD(gradup, n_cells_ext, cs_real_3_t, cs_alloc_mode);

      ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        gradup[cell_id][0] = 0.;
        gradup[cell_id][1] = 0.;
        gradup[cell_id][2] = 0.;
      });

      ctx.wait();

      cs_upwind_gradient(f_id,
                         ctx,
                         inc,
                         halo_type,
                         bc_coeffs,
                         i_massflux,
                         b_massflux,
                         _pvar,
                         gradup);

    }

  }

  /* ======================================================================
    ---> Contribution from interior faces
    ======================================================================*/

  short *i_upwind = nullptr;
  if (eqp.verbosity >= 2 && iconvp == 1 && pure_upwind == false) {
    CS_MALLOC_HD(i_upwind, n_i_faces, short, cs_alloc_mode);
    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      i_upwind[face_id] = 0;
    });
    ctx.wait();
  }

  /* --> Pure upwind flux
     =====================*/

  if (pure_upwind) {

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t cpi = 1.0, cpj = 1.0;
      if (is_thermal) {
        cpi = xcpp[ii];
        cpj = xcpp[jj];
      }

      cs_real_t fluxi = 0., fluxj = 0.;
      cs_real_t pip, pjp;

      if (ircflp == 1) {
        cs_real_t bldfrp = 1.;
        /* Local limitation of the reconstruction */
        if (df_limiter != nullptr && ircflp > 0)
          bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                           0.);

        cs_real_t recoi, recoj;
        cs_i_compute_quantities(bldfrp,
                                diipf[face_id], djjpf[face_id],
                                grad[ii], grad[jj],
                                _pvar[ii], _pvar[jj],
                                &recoi, &recoj,
                                &pip, &pjp);
      }
      else {
        pip = _pvar[ii];
        pjp = _pvar[jj];
      }

      /* No relaxation in following expressions (pipr = pip, pjpr = pjp) */

      /* Convective flux (as we are upwind, pif == _pvar[ii],
                                            pjf == _pvar[jj] */

      if (iconvp == 1) {
        cs_real_t _i_massflux = i_massflux[face_id];
        cs_real_t flui = 0.5*(_i_massflux + cs::abs(_i_massflux));
        cs_real_t fluj = 0.5*(_i_massflux - cs::abs(_i_massflux));

        fluxi += cpi*(  thetap * (flui*_pvar[ii] + fluj*_pvar[jj])
                      - imasac * _i_massflux*_pvar[ii]);

        fluxj += cpj*(  thetap * (flui*_pvar[ii] + fluj*_pvar[jj])
                      - imasac *  _i_massflux*_pvar[jj]);

        /* Fluxes without mass accumulation */
        if (store_flux) {
          i_flux[face_id][0] +=  cpi*thetap*(flui*_pvar[ii] + fluj*_pvar[jj]);
          i_flux[face_id][1] +=  cpj*thetap*(flui*_pvar[ii] + fluj*_pvar[jj]);
        }
      }

      /* Diffusive flux */

      cs_real_t diff_contrib = idiffp*thetap*i_visc[face_id]*(pip - pjp);
      fluxi += diff_contrib;
      fluxj += diff_contrib;
      /* Fluxes if needed */
      if (store_flux) {
        i_flux[face_id][0] += diff_contrib;
        i_flux[face_id][1] += diff_contrib;
      }

      if (ii < n_cells)
        cs_dispatch_sum(&rhs[ii], -fluxi, i_sum_type);
      if (jj < n_cells)
        cs_dispatch_sum(&rhs[jj],  fluxj, i_sum_type);

    });
  } /* End pure upwind */

  /* Flux with no slope test or Min/Max Beta limiter
     =============================================== */

  else if (isstpp == 1 || isstpp == 2) {

    if (ischcp < 0 || ischcp > 4) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcv"));
    }

    const cs_real_t *hybrid_blend = nullptr;
    if (CS_F_(hybrid_blend) != nullptr)
      hybrid_blend = CS_F_(hybrid_blend)->val;

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t cpi = 1.0, cpj = 1.0;
      if (is_thermal) {
        cpi = xcpp[ii];
        cpj = xcpp[jj];
      }

      cs_real_t beta = blencp;

      cs_real_t pif, pjf;
      cs_real_t pip, pjp;

#if defined(__INTEL_LLVM_COMPILER)
      // Silence unitialized variables warning due do compiler ignoring
      // initializations in inlined functions.
      pif = 0., pjf = 0.;
#endif

      /* Beta blending coefficient ensuring positivity of the scalar */
      if (isstpp == 2) {
        beta = cs::max(cs::min(cv_limiter[ii], cv_limiter[jj]),
                       0.);
      }

      cs_real_t fluxi = 0., fluxj = 0.;

      if (ircflp == 1) {
        cs_real_t bldfrp = 1.;
        if (df_limiter != nullptr)  /* Local limiter of the reconstruction */
          bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                           0.);

        cs_real_t recoi, recoj;
        cs_i_compute_quantities(bldfrp,
                                diipf[face_id], djjpf[face_id],
                                grad[ii], grad[jj],
                                _pvar[ii], _pvar[jj],
                                &recoi, &recoj,
                                &pip, &pjp);
      }
      else {
        pip = _pvar[ii];
        pjp = _pvar[jj];
      }

      const cs_real_t *cell_ceni = cell_cen[ii];
      const cs_real_t *cell_cenj = cell_cen[jj];
      const cs_real_t w_f = weight[face_id];

      if (ischcp != 4) {

        if (ischcp == 0) {

          /* Legacy SOLU
             -----------*/

          cs_solu_f_val(cell_ceni,
                        i_face_cog[face_id],
                        grad[ii],
                        _pvar[ii],
                        &pif);
          cs_solu_f_val(cell_cenj,
                        i_face_cog[face_id],
                        grad[jj],
                        _pvar[jj],
                        &pjf);
        }
        else if (ischcp == 1) {

          /* Centered
             --------*/

          pif = w_f*pip + (1.-w_f)*pjp;
          pjf = pif;

        }
        else if (ischcp == 2) {

          /* SOLU
             ---- */

          cs_solu_f_val(cell_ceni,
                        i_face_cog[face_id],
                        gradup[ii],
                        _pvar[ii],
                        &pif);
          cs_solu_f_val(cell_cenj,
                        i_face_cog[face_id],
                        gradup[jj],
                        _pvar[jj],
                        &pjf);

        }
        else {

          assert(ischcp == 3);

          /* Centered
             -------- */

          pif = w_f*pip + (1.-w_f)*pjp;
          pjf = pif;

          /* Legacy SOLU
             -----------*/
          cs_real_t pif_up, pjf_up;

          cs_solu_f_val(cell_ceni,
                        i_face_cog[face_id],
                        grad[ii],
                        _pvar[ii],
                        &pif_up);

          cs_solu_f_val(cell_cenj,
                        i_face_cog[face_id],
                        grad[jj],
                        _pvar[jj],
                        &pjf_up);

          cs_real_t hybrid_blend_interp
            = cs::min(hybrid_blend[ii], hybrid_blend[jj]);

          pif = hybrid_blend_interp*pif + (1. - hybrid_blend_interp)*pif_up;
          pjf = hybrid_blend_interp*pjf + (1. - hybrid_blend_interp)*pjf_up;

        }

        /* Blending
           -------- */

        pif = beta * pif + (1. - beta) * _pvar[ii];
        pjf = beta * pjf + (1. - beta) * _pvar[jj];

      }
      else if (ischcp == 4) {

        /* NVD/TVD family of high accuracy schemes */

        cs_lnum_t ic, id;

        /* Determine central and downwind sides w.r.t. current face */
        if (i_massflux[face_id] >= 0.) {
          ic = ii;
          id = jj;
        } else {
          ic = jj;
          id = ii;
        }

        cs_real_t courant_c = -1.;
        if (courant != nullptr && is_thermal == false)
          courant_c = courant[ic];

        cs_i_cd_unsteady_nvd(limiter_choice,
                             beta,
                             cell_cen[ic],
                             cell_cen[id],
                             i_face_u_normal[face_id],
                             i_face_cog[face_id],
                             grad[ic],
                             _pvar[ic],
                             _pvar[id],
                             local_max[ic],
                             local_min[ic],
                             courant_c,
                             &pif,
                             &pjf);

      }

      // Convective flux

      if (iconvp == 1) {
        cs_real_t _i_massflux = i_massflux[face_id];
        cs_real_t flui = 0.5*(_i_massflux + cs::abs(_i_massflux));
        cs_real_t fluj = 0.5*(_i_massflux - cs::abs(_i_massflux));

        fluxi += cpi*(  thetap*(flui*pif + fluj*pjf)
                      - imasac*_i_massflux*_pvar[ii]);

        fluxj += cpj*(  thetap*(flui*pif + fluj*pjf)
                      - imasac*_i_massflux*_pvar[jj]);

        /* Fluxes without mass accumulation if needed */
        if (store_flux) {
          i_flux[face_id][0] +=  cpi*thetap*(flui*pif + fluj*pjf);
          i_flux[face_id][1] +=  cpj*thetap*(flui*pif + fluj*pjf);
        }
      }

      // Diffusive flux (no relaxation)

      cs_real_t diff_contrib = idiffp*thetap*i_visc[face_id]*(pip - pjp);
      fluxi += diff_contrib;
      fluxj += diff_contrib;

      /* Fluxes if needed */
      if (store_flux) {
        i_flux[face_id][0] += diff_contrib;
        i_flux[face_id][1] += diff_contrib;
      }

      if (ii < n_cells)
        cs_dispatch_sum(&rhs[ii], -fluxi, i_sum_type);
      if (jj < n_cells)
        cs_dispatch_sum(&rhs[jj],  fluxj, i_sum_type);

    });

  /* --> Flux with slope test
     ============================================*/

  }
  else { /* isstpp = 0 */

    if (ischcp < 0 || ischcp > 2) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcv"));
    }

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t cpi = 1., cpj = 1.;
      if (is_thermal) {
        cpi = xcpp[ii];
        cpj = xcpp[jj];
      }

      bool upwind_switch = false;

      cs_real_t fluxi = 0., fluxj = 0.;

      cs_real_t pif, pjf;
      cs_real_t pip, pjp;

      if (ircflp == 1) {
        cs_real_t bldfrp = 1.;
        if (df_limiter != nullptr)  /* Local limiter */
          bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                           0.);

        cs_real_t recoi, recoj;
        cs_i_compute_quantities(bldfrp,
                                diipf[face_id], djjpf[face_id],
                                grad[ii], grad[jj],
                                _pvar[ii], _pvar[jj],
                                &recoi, &recoj,
                                &pip, &pjp);
      }
      else {
        pip = _pvar[ii];
        pjp = _pvar[jj];
      }

      const cs_real_t *cell_ceni = cell_cen[ii];
      const cs_real_t *cell_cenj = cell_cen[jj];
      const cs_real_t w_f = weight[face_id];

      /* Slope test is needed with convection */
      if (iconvp > 0) {
        cs_real_t testij, tesqck;

        cs_slope_test(_pvar[ii],
                      _pvar[jj],
                      i_dist[face_id],
                      i_face_u_normal[face_id],
                      grad[ii],
                      grad[jj],
                      gradst[ii],
                      gradst[jj],
                      i_massflux[face_id],
                      &testij,
                      &tesqck);

        if (ischcp == 0) {

          /* Original SOLU
             --------------*/

          cs_solu_f_val(cell_ceni,
                        i_face_cog[face_id],
                        grad[ii],
                        _pvar[ii],
                        &pif);
          cs_solu_f_val(cell_cenj,
                        i_face_cog[face_id],
                        grad[jj],
                        _pvar[jj],
                        &pjf);
        }
        else if (ischcp == 1) {

          /* Centered
             --------*/

          pif = w_f*pip + (1.-w_f)*pjp;
          pjf = pif;

        }
        else {

          /* SOLU
             -----*/

          cs_solu_f_val(cell_ceni,
                        i_face_cog[face_id],
                        gradup[ii],
                        _pvar[ii],
                        &pif);
          cs_solu_f_val(cell_cenj,
                        i_face_cog[face_id],
                        gradup[jj],
                        _pvar[jj],
                        &pjf);

        }

        /* Slope test: percentage of upwind
           -------------------------------- */

        if (tesqck <= 0. || testij <= 0.) {

          cs_blend_f_val(blend_st, _pvar[ii], &pif);
          cs_blend_f_val(blend_st, _pvar[jj], &pjf);

          upwind_switch = true;

        }

        /* Blending
           -------- */

        cs_blend_f_val(blencp, _pvar[ii], &pif);
        cs_blend_f_val(blencp, _pvar[jj], &pjf);

      }
      else { /* If iconv=0 p*fr* are useless */

        pif = _pvar[ii];
        pjf = _pvar[jj];

      } /* End for slope test */

      // Convective flux

      if (iconvp == 1) {
        cs_real_t _i_massflux = i_massflux[face_id];
        cs_real_t flui = 0.5*(_i_massflux + cs::abs(_i_massflux));
        cs_real_t fluj = 0.5*(_i_massflux - cs::abs(_i_massflux));

        fluxi += cpi*(  thetap*(flui*pif + fluj*pjf)
                      - imasac*_i_massflux*_pvar[ii]);

        fluxj += cpj*(  thetap*(flui*pif + fluj*pjf)
                      - imasac*_i_massflux*_pvar[jj]);

        /* Fluxes without mass accumulation if needed */
        if (store_flux) {
          i_flux[face_id][0] +=  cpi*thetap*(flui*pif + fluj*pjf);
          i_flux[face_id][1] +=  cpj*thetap*(flui*pif + fluj*pjf);
        }

      }

      // Diffusive flux (no relaxation)

      cs_real_t diff_contrib = idiffp*thetap*i_visc[face_id]*(pip - pjp);
      fluxi += diff_contrib;
      fluxj += diff_contrib;

      /* Fluxes if needed */
      if (store_flux) {
        i_flux[face_id][0] += diff_contrib;
        i_flux[face_id][1] += diff_contrib;
      }

      if (upwind_switch) {
        /* in parallel, face will be counted by one and only one rank */
        if (i_upwind != nullptr && ii < n_cells)
          i_upwind[face_id] = 1;

        if (v_slope_test != nullptr) {
          cs_real_t q_d_vol_ii
            =   cs::abs(i_massflux[face_id])
              * cs_mq_cell_vol_inv(ii, c_disable_flag, cell_vol);
          cs_real_t q_d_vol_jj
            =   cs::abs(i_massflux[face_id])
              * cs_mq_cell_vol_inv(jj, c_disable_flag, cell_vol);

          cs_dispatch_sum(&v_slope_test[ii], q_d_vol_ii, i_sum_type);
          cs_dispatch_sum(&v_slope_test[jj], q_d_vol_jj, i_sum_type);
        }

      }

      if (ii < n_cells)
        cs_dispatch_sum(&rhs[ii], -fluxi, i_sum_type);
      if (jj < n_cells)
        cs_dispatch_sum(&rhs[jj],  fluxj, i_sum_type);

    });
  } /* pure upwind, without slope test, with slope test */

  ctx.wait();

  if (eqp.verbosity >= 2) {
    cs_gnum_t n_upwind = 0;

    if (i_upwind != nullptr) {
      ctx.parallel_for_reduce_sum
        (n_i_faces, n_upwind, [=] CS_F_HOST_DEVICE
         (cs_lnum_t i,
          CS_DISPATCH_REDUCER_TYPE(cs_gnum_t) &sum) {
        sum += (cs_gnum_t)i_upwind[i];
      });

      ctx.wait();
      cs_parall_counter(&n_upwind, 1);
    }
    else if (pure_upwind)
      n_upwind = m->n_g_i_faces;

    bft_printf(_(" %s: %llu Faces with upwind on %llu interior faces\n"),
               var_name, (unsigned long long)n_upwind,
               (unsigned long long)m->n_g_i_c_faces);

    CS_FREE_HD(i_upwind);
  }

  /* ======================================================================
     ---> Contribution from boundary faces
     ======================================================================*/

  /* Boundary convective flux are all computed with an upwind scheme */
  if (icvflb == false || is_thermal) {

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = b_face_cells[face_id];

      cs_real_t cpi = 1.;
      if (is_thermal)
        cpi = xcpp[ii];

      cs_real_t fluxi = 0.;
      cs_real_t pip;

      cs_real_t bldfrp = (cs_real_t) ircflb;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflb > 0)
        bldfrp = cs::max(df_limiter[ii], 0.);

      cs_b_cd_unsteady(bldfrp,
                       diipb[face_id],
                       grad[ii],
                       _pvar[ii],
                       &pip);

      cs_b_upwind_flux(iconvp,
                       thetap,
                       imasac,
                       inc,
                       bc_type[face_id],
                       _pvar[ii],
                       _pvar[ii], /* no relaxation */
                       pip,
                       coefap[face_id],
                       coefbp[face_id],
                       b_massflux[face_id],
                       cpi,
                       &fluxi);

      cs_b_diff_flux(idiffp,
                     thetap,
                     inc,
                     pip,
                     cofafp[face_id],
                     cofbfp[face_id],
                     b_visc[face_id],
                     &fluxi);

      cs_dispatch_sum(&rhs[ii], -fluxi, b_sum_type);

      /* Fluxes without mass accumulation if needed */
      if (store_flux) {
        fluxi = 0.;

        cs_b_upwind_flux(iconvp,
                         thetap,
                         0, //imasac,
                         inc,
                         bc_type[face_id],
                         _pvar[ii],
                         _pvar[ii], /* no relaxation */
                         pip,
                         coefap[face_id],
                         coefbp[face_id],
                         b_massflux[face_id],
                         cpi,
                         &fluxi);

        cs_b_diff_flux(idiffp,
                       thetap,
                       inc,
                       pip,
                       cofafp[face_id],
                       cofbfp[face_id],
                       b_visc[face_id],
                       &fluxi);

        b_flux[face_id] += fluxi;
      }

    });

    /* The scalar is internally coupled and an implicit contribution
     * is required */
    if (icoupl > 0) {
      /* Prepare data for sending */
      CS_MALLOC_HD(pvar_distant, n_distant, cs_real_t, cs_alloc_mode);

      for (cs_lnum_t ii = 0; ii < n_distant; ii++) {
        cs_lnum_t face_id = faces_distant[ii];
        cs_lnum_t jj = b_face_cells[face_id];
        cs_real_t pip;

        cs_real_t bldfrp = (cs_real_t) ircflb;
        /* Local limitation of the reconstruction */
        if (df_limiter != nullptr && ircflb > 0)
          bldfrp = cs::max(df_limiter[jj], 0.);

        cs_b_cd_unsteady(bldfrp,
                         diipb[face_id],
                         grad[jj],
                         _pvar[jj],
                         &pip);
        pvar_distant[ii] = pip;
      }

      /* Receive data */
      CS_MALLOC_HD(pvar_local, n_local, cs_real_t, cs_alloc_mode);
      cs_internal_coupling_exchange_var(cpl,
                                        1, /* Dimension */
                                        pvar_distant,
                                        pvar_local);

      /* Exchange diffusion limiter */
      if (df_limiter != nullptr) {
        CS_MALLOC_HD(df_limiter_local, n_local, cs_real_t, cs_alloc_mode);
        cs_internal_coupling_exchange_var(cpl,
                                          1, /* Dimension */
                                          df_limiter,
                                          df_limiter_local);
      }

      ctx.wait(); // Next loop will contribute to rhs.

      /* Flux contribution */
      assert(f != nullptr);
      cs_real_t *hintp = f->bc_coeffs->hint;
      cs_real_t *hextp = f->bc_coeffs->rcodcl2;
      for (cs_lnum_t ii = 0; ii < n_local; ii++) {
        cs_lnum_t face_id = faces_local[ii];
        cs_lnum_t jj = b_face_cells[face_id];
        cs_real_t surf = b_face_surf[face_id];
        cs_real_t pip, pjp;
        cs_real_t fluxi = 0.;

        cs_real_t bldfrp = (cs_real_t) ircflb;
        /* Local limitation of the reconstruction */
        if (df_limiter != nullptr && ircflb > 0)
          bldfrp = cs::max(cs::min(df_limiter_local[ii],
                                   df_limiter[jj]),
                           0.);

        cs_b_cd_unsteady(bldfrp,
                         diipb[face_id],
                         grad[jj],
                         _pvar[jj],
                         &pip);

        pjp = pvar_local[ii];

        cs_real_t hint = hintp[face_id];
        cs_real_t hext = hextp[face_id];
        cs_real_t heq = _calc_heq(hint, hext)*surf;

        cs_b_diff_flux_coupling(idiffp,
                                pip,
                                pjp,
                                heq,
                                &fluxi);

        if (store_flux)
          b_flux[face_id] += thetap * fluxi;

        rhs[jj] -= thetap * fluxi;
      }

      CS_FREE_HD(pvar_local);
      /* Sending structures are no longer needed */
      CS_FREE_HD(pvar_distant);
      if (df_limiter != nullptr)
        CS_FREE_HD(df_limiter_local);
    }
  }

  /* Boundary convective flux is imposed at some faces
     (tagged in icvfli array) */

  else if (icvflb && is_thermal == false) {

    /* Retrieve the value of the convective flux to be imposed */
    if (f_id != -1) {
      coface = f->bc_coeffs->ac;
      cofbce = f->bc_coeffs->bc;
    }
    else {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of icvflb and f_id"));
    }

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = b_face_cells[face_id];

      cs_real_t fluxi = 0.;
      cs_real_t pip;

      cs_real_t bldfrp = (cs_real_t) ircflb;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflb > 0)
        bldfrp = cs::max(df_limiter[ii], 0.);

      cs_b_cd_unsteady(bldfrp,
                       diipb[face_id],
                       grad[ii],
                       _pvar[ii],
                       &pip);

      cs_b_imposed_conv_flux(iconvp,
                             thetap,
                             imasac,
                             inc,
                             bc_type[face_id],
                             icvfli[face_id],
                             _pvar[ii],
                             _pvar[ii], /* no relaxation */
                             pip,
                             coefap[face_id],
                             coefbp[face_id],
                             coface[face_id],
                             cofbce[face_id],
                             b_massflux[face_id],
                             1., /* xcpp */
                             &fluxi);

      cs_b_diff_flux(idiffp,
                     thetap,
                     inc,
                     pip,
                     cofafp[face_id],
                     cofbfp[face_id],
                     b_visc[face_id],
                     &fluxi);

      cs_dispatch_sum(&rhs[ii], -fluxi, b_sum_type);

      /* Fluxes without mass accumulation if needed */
      if (store_flux) {
        fluxi = 0.;

        cs_b_imposed_conv_flux(iconvp,
                               thetap,
                               0, //imasac,
                               inc,
                               bc_type[face_id],
                               icvfli[face_id],
                               _pvar[ii],
                               _pvar[ii], /* no relaxation */
                               pip,
                               coefap[face_id],
                               coefbp[face_id],
                               coface[face_id],
                               cofbce[face_id],
                               b_massflux[face_id],
                               1., /* xcpp */
                               &fluxi);

        cs_b_diff_flux(idiffp,
                       thetap,
                       inc,
                       pip,
                       cofafp[face_id],
                       cofbfp[face_id],
                       b_visc[face_id],
                       &fluxi);

        b_flux[face_id] += fluxi;
      }

    });

  }

  ctx.wait();

  /* Free memory */
  CS_FREE_HD(grad);
  CS_FREE_HD(gradup);
  CS_FREE_HD(gradst);
  CS_FREE_HD(local_max);
  CS_FREE_HD(local_min);
  CS_FREE_HD(courant);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update face flux with convection contribution of a standard transport
 * equation of a scalar field \f$ \varia \f$.
 *
 * <a name="cs_face_convection_scalar"></a>
 *
 * \f[
 * C_\ij = \dot{m}_\ij \left( \varia_\fij - \varia_\celli \right)
 * \f]
 *
 * \param[in]     f_id          pointer to field, or nullptr
 * \param[in]     eqp           equation parameters
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imasac        take mass accumulation into account?
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in]     bc_coeffs     boundary condition structure for the variable
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in,out] i_conv_flux   scalar convection flux at interior faces
 * \param[in,out] b_conv_flux   scalar convection flux at boundary faces
 */
/*----------------------------------------------------------------------------*/

static void
_face_convection_scalar_unsteady(const cs_field_t           *f,
                                 const cs_equation_param_t   eqp,
                                 int                         icvflb,
                                 int                         inc,
                                 int                         imasac,
                                 cs_real_t         *restrict pvar,
                                 const cs_real_t   *restrict pvara,
                                 const int                   icvfli[],
                                 const cs_field_bc_coeffs_t *bc_coeffs,
                                 const cs_real_t             i_massflux[],
                                 const cs_real_t             b_massflux[],
                                 cs_real_t                   i_conv_flux[][2],
                                 cs_real_t                   b_conv_flux[])
{

  cs_real_t *coefap = bc_coeffs->a;
  cs_real_t *coefbp = bc_coeffs->b;

  const int iconvp = eqp.iconv;
  const int nswrgp = eqp.nswrgr;
  const int ircflp = eqp.ircflu;
  const int ircflb = (ircflp > 0) ? eqp.b_diff_flux_rc : 0;
  const int ischcp = eqp.ischcv;
  const int isstpp = eqp.isstpc;
  const int iwarnp = eqp.verbosity;
  const int icoupl = eqp.icoupl;
  cs_nvd_type_t limiter_choice = CS_NVD_N_TYPES;
  const double blencp = eqp.blencv;
  const double blend_st = eqp.blend_st;
  const double thetap = eqp.theta;

  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_faces = m->n_i_faces;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_nreal_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_rreal_3_t *restrict diipf = fvq->diipf;
  const cs_rreal_3_t *restrict djjpf = fvq->djjpf;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;

  const int *bc_type = cs_glob_bc_type;

  const int *restrict c_disable_flag = (fvq->has_disable_flag) ?
    fvq->c_disable_flag : nullptr;

  /* Parallel or device dispatch */
  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);

  /* Local variables */

  const int f_id = (f != nullptr) ? f->id : -1;
  char var_name[64];

  int w_stride = 1;

  cs_real_t *coface = nullptr, *cofbce = nullptr;

  cs_real_3_t *grad;
  cs_real_3_t *gradup = nullptr;
  cs_real_3_t *gradst = nullptr;

  cs_real_t *local_min = nullptr;
  cs_real_t *local_max = nullptr;
  cs_real_t *courant = nullptr;

  cs_real_t *cv_limiter = nullptr;
  cs_real_t *df_limiter = nullptr;

  cs_real_t *gweight = nullptr;

  const int key_lim_choice = cs_field_key_id("limiter_choice");

  cs_real_t  *v_slope_test = cs_get_v_slope_test(f_id,  eqp);

  /* Internal coupling variables */
  int coupling_id;
  cs_internal_coupling_t *cpl = nullptr;

  /* Allocate work arrays */

  CS_MALLOC_HD(grad, n_cells_ext, cs_real_3_t, cs_alloc_mode);

  /* Choose gradient type */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  cs_gradient_type_by_imrgra(eqp.d_gradient_r,
                             &gradient_type,
                             &halo_type);

  /* Handle cases where only the previous values (already synchronized)
     or current values are provided */

  if (pvar != nullptr)
    cs_sync_scalar_halo(m, halo_type, pvar);
  if (pvara == nullptr)
    pvara = (const cs_real_t *)pvar;

  const cs_real_t  *restrict _pvar = (pvar != nullptr) ? pvar : pvara;

  /* Slope limiters */

  if (f != nullptr) {

    /* NVD/TVD limiters */
    if (ischcp == 4) {
      limiter_choice = (cs_nvd_type_t)cs_field_get_key_int(f, key_lim_choice);
      CS_MALLOC_HD(local_max, n_cells_ext, cs_real_t, cs_alloc_mode);
      CS_MALLOC_HD(local_min, n_cells_ext, cs_real_t, cs_alloc_mode);
      cs_field_local_extrema_scalar(f_id,
                                    halo_type,
                                    local_max,
                                    local_min);
      if (limiter_choice >= CS_NVD_VOF_HRIC) {
        CS_MALLOC_HD(courant, n_cells_ext, cs_real_t, cs_alloc_mode);
        cs_cell_courant_number(f, ctx, courant);
      }
    }

    int cv_limiter_id =
      cs_field_get_key_int(f, cs_field_key_id("convection_limiter_id"));
    if (cv_limiter_id > -1)
      cv_limiter = cs_field_by_id(cv_limiter_id)->val;

    int df_limiter_id =
      cs_field_get_key_int(f, cs_field_key_id("diffusion_limiter_id"));
    if (df_limiter_id > -1)
      df_limiter = cs_field_by_id(df_limiter_id)->val;

    snprintf(var_name, 63, "%s", f->name);
  }
  else if (isstpp > 1) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of isstpp for a work array"));
  }
  else {
    strncpy(var_name, "[scalar face flux from convection]", 63);
  }
  var_name[63] = '\0';

  if (iwarnp >= 2) {
    if (ischcp == 1) {
      bft_printf(
        _(" %s: Convection in centered blending with %f percent of upwind\n"),
        var_name, (1.-blencp)*100.);
    }
    else {
      bft_printf(
        _(" %s: Convection in 2nd order blending with %f percent of upwind\n"),
        var_name, (1.-blencp)*100.);
    }
  }

  bool pure_upwind = (blencp > 0.) ? false : true;

  if (icoupl > 0) {
    assert(f_id != -1);
    const int coupling_key_id = cs_field_key_id("coupling_entity");
    coupling_id = cs_field_get_key_int(f, coupling_key_id);
    cpl = cs_internal_coupling_by_id(coupling_id);
  }

  /* Compute the balance with reconstruction */

  /* Compute the gradient of the variable */

  /* The gradient (grad) is used in the flux reconstruction and the slope test.
     Thus we must compute it:
         - when the convection scheme is the legacy SOLU,
         - when we have convection, we are not in pure upwind
           and we reconstruct the fluxes,
         - when we have convection, we are not in pure upwind
           and we have not shunted the slope test,
         - when we use NVD / TVD schemes.
  */

  if (   iconvp != 0 && pure_upwind == false
      && (ischcp == 0 || ircflp == 1 || isstpp == 0 || ischcp == 4)) {

    if (f != nullptr) {
      /* Get the calculation option from the field */
      if (f->type & CS_FIELD_VARIABLE && eqp.iwgrec == 1) {
        if (eqp.idiff > 0) {
          int key_id = cs_field_key_id("gradient_weighting_id");
          int diff_id = cs_field_get_key_int(f, key_id);
          if (diff_id > -1) {
            cs_field_t *weight_f = cs_field_by_id(diff_id);
            gweight = weight_f->val;
            w_stride = weight_f->dim;
            cs_field_synchronize(weight_f, halo_type);
          }
        }
      }
    }

    cs_gradient_scalar_synced_input(var_name,
                                    gradient_type,
                                    halo_type,
                                    inc,
                                    nswrgp,
                                    0, /* hyd_p_flag */
                                    w_stride,
                                    eqp.verbosity,
                                    (cs_gradient_limit_t)eqp.imligr,
                                    eqp.epsrgr,
                                    eqp.d_climgr,
                                    nullptr, /* f_ext exterior force */
                                    bc_coeffs,
                                    _pvar,
                                    gweight, /* Weighted gradient */
                                    cpl,
                                    grad);

  }
  else {

    ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      grad[cell_id][0] = 0.;
      grad[cell_id][1] = 0.;
      grad[cell_id][2] = 0.;
    });
    ctx.wait();
  }

  /* Compute gradients used in convection schemes */

  if (iconvp > 0 && pure_upwind == false) {

    /* Compute cell gradient used in slope test */
    if (isstpp == 0) {

      CS_MALLOC_HD(gradst, n_cells_ext, cs_real_3_t, cs_alloc_mode);

      ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        gradst[cell_id][0] = 0.;
        gradst[cell_id][1] = 0.;
        gradst[cell_id][2] = 0.;
      });
      ctx.wait();

      cs_slope_test_gradient(f_id,
                             ctx,
                             inc,
                             (const cs_real_3_t *)grad,
                             gradst,
                             _pvar,
                             bc_coeffs,
                             i_massflux);

    }

    /* Pure SOLU scheme */
    if (ischcp == 2) {

      CS_MALLOC_HD(gradup, n_cells_ext, cs_real_3_t, cs_alloc_mode);

      ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        gradup[cell_id][0] = 0.;
        gradup[cell_id][1] = 0.;
        gradup[cell_id][2] = 0.;
      });
      ctx.wait();

      cs_upwind_gradient(f_id,
                         ctx,
                         inc,
                         halo_type,
                         bc_coeffs,
                         i_massflux,
                         b_massflux,
                         _pvar,
                         gradup);

    }

  }

  /* ======================================================================
    ---> Contribution from interior faces
    ======================================================================*/

  short *i_upwind = nullptr;
  if (eqp.verbosity >= 2 && iconvp == 1 && pure_upwind == false) {
    CS_MALLOC_HD(i_upwind, n_i_faces, short, cs_alloc_mode);
    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      i_upwind[face_id] = 0;
    });
    ctx.wait();
  }

  /* --> Pure upwind flux
    =====================*/

  if (pure_upwind) {

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t pip, pjp;

      if (ircflp == 1) {
        cs_real_t bldfrp = 1.;
        /* Local limitation of the reconstruction */
        if (df_limiter != nullptr)
          bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                           0.);

        cs_real_t recoi, recoj;
        cs_i_compute_quantities(bldfrp,
                                diipf[face_id], djjpf[face_id],
                                grad[ii], grad[jj],
                                _pvar[ii], _pvar[jj],
                                &recoi, &recoj,
                                &pip, &pjp);
      }
      else {
        pip = _pvar[ii];
        pjp = _pvar[jj];
      }

      /* No relaxation in following expressions (pipr = pip, pjpr = pjp) */

      /* Convective flux (as we are upwind, pif == _pvar[ii],
                                            pjf == _pvar[jj] */

      if (iconvp == 1) {
        cs_real_t _i_massflux = i_massflux[face_id];
        cs_real_t flui = 0.5*(_i_massflux + cs::abs(_i_massflux));
        cs_real_t fluj = 0.5*(_i_massflux - cs::abs(_i_massflux));

        i_conv_flux[face_id][0] +=   thetap*(flui*_pvar[ii] + fluj*_pvar[jj])
                                   - imasac * _i_massflux*_pvar[ii];

        i_conv_flux[face_id][1] +=   thetap * (flui*_pvar[ii] + fluj*_pvar[jj])
                                   - imasac * _i_massflux*_pvar[jj];

      }

    });
  } /* End pure upwind */

  /* Flux with no slope test or Min/Max Beta limiter
     =============================================== */

  else if (isstpp == 1 || isstpp == 2) {

    if (ischcp < 0 || ischcp > 4) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcv"));
    }

    const cs_real_t *hybrid_blend = nullptr;
    if (CS_F_(hybrid_blend) != nullptr)
      hybrid_blend = CS_F_(hybrid_blend)->val;

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t beta = blencp;

      cs_real_t pif = 0., pjf = 0.;
      cs_real_t pip = 0., pjp = 0.;

      /* Beta blending coefficient ensuring positivity of the scalar */
      if (isstpp == 2) {
        beta = cs::max(cs::min(cv_limiter[ii], cv_limiter[jj]),
                       0.);
      }

      if (ircflp == 1) {
        cs_real_t bldfrp = 1.;
        if (df_limiter != nullptr)  /* Local limiter of the reconstruction */
          bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                           0.);

        cs_real_t recoi, recoj;
        cs_i_compute_quantities(bldfrp,
                                diipf[face_id], djjpf[face_id],
                                grad[ii], grad[jj],
                                _pvar[ii], _pvar[jj],
                                &recoi, &recoj,
                                &pip, &pjp);
      }
      else {
        pip = _pvar[ii];
        pjp = _pvar[jj];
      }

      const cs_real_t *cell_ceni = cell_cen[ii];
      const cs_real_t *cell_cenj = cell_cen[jj];
      const cs_real_t w_f = weight[face_id];

      if (ischcp != 4) {

        if (ischcp == 0) {

          /* Legacy SOLU
             -----------*/

          cs_solu_f_val(cell_ceni,
                        i_face_cog[face_id],
                        grad[ii],
                        _pvar[ii],
                        &pif);
          cs_solu_f_val(cell_cenj,
                        i_face_cog[face_id],
                        grad[jj],
                        _pvar[jj],
                        &pjf);
        }
        else if (ischcp == 1) {

          /* Centered
             --------*/

          pif = w_f*pip + (1.-w_f)*pjp;
          pjf = pif;

        }
        else if (ischcp == 2) {

          /* SOLU
             ---- */

          cs_solu_f_val(cell_ceni,
                        i_face_cog[face_id],
                        gradup[ii],
                        _pvar[ii],
                        &pif);
          cs_solu_f_val(cell_cenj,
                        i_face_cog[face_id],
                        gradup[jj],
                        _pvar[jj],
                        &pjf);

        }
        else if (ischcp == 3) {

          /* Centered
             -------- */

          pif = w_f*pip + (1.-w_f)*pjp;
          pjf = pif;

          /* Legacy SOLU
             -----------*/
          cs_real_t pif_up, pjf_up;

          cs_solu_f_val(cell_ceni,
                        i_face_cog[face_id],
                        grad[ii],
                        _pvar[ii],
                        &pif_up);

          cs_solu_f_val(cell_cenj,
                        i_face_cog[face_id],
                        grad[jj],
                        _pvar[jj],
                        &pjf_up);

          cs_real_t hybrid_blend_interp
            = cs::min(hybrid_blend[ii], hybrid_blend[jj]);

          pif = hybrid_blend_interp*pif + (1. - hybrid_blend_interp)*pif_up;
          pjf = hybrid_blend_interp*pjf + (1. - hybrid_blend_interp)*pjf_up;

        }

        /* Blending
           -------- */

        pif = beta * pif + (1. - beta) * _pvar[ii];
        pjf = beta * pjf + (1. - beta) * _pvar[jj];

      }
      else if (ischcp == 4) {

        /* NVD/TVD family of high accuracy schemes */

        cs_lnum_t ic, id;

        /* Determine central and downwind sides w.r.t. current face */
        if (i_massflux[face_id] >= 0.) {
          ic = ii;
          id = jj;
        }
        else {
          ic = jj;
          id = ii;
        }

        cs_real_t courant_c = -1.;
        if (courant != nullptr)
          courant_c = courant[ic];

        cs_i_cd_unsteady_nvd(limiter_choice,
                             beta,
                             cell_cen[ic],
                             cell_cen[id],
                             i_face_u_normal[face_id],
                             i_face_cog[face_id],
                             grad[ic],
                             _pvar[ic],
                             _pvar[id],
                             local_max[ic],
                             local_min[ic],
                             courant_c,
                             &pif,
                             &pjf);

      }

      // Convective flux

      if (iconvp == 1) {
        cs_real_t _i_massflux = i_massflux[face_id];
        cs_real_t flui = 0.5*(_i_massflux + cs::abs(_i_massflux));
        cs_real_t fluj = 0.5*(_i_massflux - cs::abs(_i_massflux));

        i_conv_flux[face_id][0] +=   thetap*(flui*pif + fluj*pjf)
                                   - imasac*_i_massflux*_pvar[ii];

        i_conv_flux[face_id][1] +=   thetap*(flui*pif + fluj*pjf)
                                   - imasac*_i_massflux*_pvar[jj];

      }

    });

  }
  /* Flux with slope test or NVD/TVD limiter
     ======================================= */

  else { /* isstpp = 0 */

    if (ischcp < 0 || ischcp > 2) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcv"));
    }

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      bool upwind_switch = false;

      cs_real_t pif, pjf;
      cs_real_t pip, pjp;

      if (ircflp == 1) {
        cs_real_t bldfrp = 1.;
        if (df_limiter != nullptr)  /* Local limiter */
          bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                           0.);

        cs_real_t recoi, recoj;
        cs_i_compute_quantities(bldfrp,
                                diipf[face_id], djjpf[face_id],
                                grad[ii], grad[jj],
                                _pvar[ii], _pvar[jj],
                                &recoi, &recoj,
                                &pip, &pjp);
      }
      else {
        pip = _pvar[ii];
        pjp = _pvar[jj];
      }

      const cs_real_t *cell_ceni = cell_cen[ii];
      const cs_real_t *cell_cenj = cell_cen[jj];
      const cs_real_t w_f = weight[face_id];

      /* Slope test is needed with convection */
      if (iconvp > 0) {
        cs_real_t testij, tesqck;

        cs_slope_test(_pvar[ii],
                      _pvar[jj],
                      i_dist[face_id],
                      i_face_u_normal[face_id],
                      grad[ii],
                      grad[jj],
                      gradst[ii],
                      gradst[jj],
                      i_massflux[face_id],
                      &testij,
                      &tesqck);

        if (ischcp == 0) {

          /* Original SOLU
             --------------*/

          cs_solu_f_val(cell_ceni,
                        i_face_cog[face_id],
                        grad[ii],
                        _pvar[ii],
                        &pif);
          cs_solu_f_val(cell_cenj,
                        i_face_cog[face_id],
                        grad[jj],
                        _pvar[jj],
                        &pjf);
        }
        else if (ischcp == 1) {

          /* Centered
             --------*/

          pif = w_f*pip + (1.-w_f)*pjp;
          pjf = pif;

        }
        else {

          /* SOLU
             -----*/

          cs_solu_f_val(cell_ceni,
                        i_face_cog[face_id],
                        gradup[ii],
                        _pvar[ii],
                        &pif);
          cs_solu_f_val(cell_cenj,
                        i_face_cog[face_id],
                        gradup[jj],
                        _pvar[jj],
                        &pjf);

        }

        /* Slope test: percentage of upwind
           -------------------------------- */

        if (tesqck <= 0. || testij <= 0.) {

          cs_blend_f_val(blend_st, _pvar[ii], &pif);
          cs_blend_f_val(blend_st, _pvar[jj], &pjf);

          upwind_switch = true;

        }

        /* Blending
           -------- */

        cs_blend_f_val(blencp, _pvar[ii], &pif);
        cs_blend_f_val(blencp, _pvar[jj], &pjf);

      }
      else { /* If iconv=0 p*fr* are useless */

        pif = _pvar[ii];
        pjf = _pvar[jj];

      } /* End for slope test */

      if (iconvp == 1) {
        cs_real_t _i_massflux = i_massflux[face_id];
        cs_real_t flui = 0.5*(_i_massflux + cs::abs(_i_massflux));
        cs_real_t fluj = 0.5*(_i_massflux - cs::abs(_i_massflux));

        i_conv_flux[face_id][0] +=   thetap*(flui*pif + fluj*pjf)
                                   - imasac*_i_massflux*_pvar[ii];

        i_conv_flux[face_id][1] +=   thetap*(flui*pif + fluj*pjf)
                                   - imasac*_i_massflux*_pvar[jj];
      }

      if (upwind_switch) {
        /* in parallel, face will be counted by one and only one rank */
        if (i_upwind != nullptr && ii < n_cells)
          i_upwind[face_id] = 1;

        if (v_slope_test != nullptr) {
          cs_real_t q_d_vol_ii
            =   cs::abs(i_massflux[face_id])
              * cs_mq_cell_vol_inv(ii, c_disable_flag, cell_vol);
          cs_real_t q_d_vol_jj
            =   cs::abs(i_massflux[face_id])
              * cs_mq_cell_vol_inv(jj, c_disable_flag, cell_vol);

          cs_dispatch_sum(&v_slope_test[ii], q_d_vol_ii, i_sum_type);
          cs_dispatch_sum(&v_slope_test[jj], q_d_vol_jj, i_sum_type);
        }
      }
    });

  } /* pure upwind, without slope test, with slope test */

  if (eqp.verbosity >= 2 && iconvp == 1) {
    cs_gnum_t n_upwind = 0;

    if (i_upwind != nullptr) {
      ctx.parallel_for_reduce_sum
        (n_i_faces, n_upwind, [=] CS_F_HOST_DEVICE
         (cs_lnum_t i,
          CS_DISPATCH_REDUCER_TYPE(cs_gnum_t) &sum) {
        sum += (cs_gnum_t)i_upwind[i];
      });

      ctx.wait();
      cs_parall_counter(&n_upwind, 1);
    }
    else if (pure_upwind)
      n_upwind = m->n_g_i_faces;

    bft_printf(_(" %s: %llu Faces with upwind on %llu interior faces\n"),
               var_name, (unsigned long long)n_upwind,
               (unsigned long long)m->n_g_i_c_faces);

    CS_FREE_HD(i_upwind);
  }

  /* ======================================================================
    ---> Contribution from boundary faces
    ======================================================================*/

  /* Boundary convective flux are all computed with an upwind scheme */
  if (icvflb == false) {

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = b_face_cells[face_id];
      cs_real_t pip;

      cs_real_t bldfrp = (cs_real_t) ircflb;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflb > 0)
        bldfrp = cs::max(df_limiter[ii], 0.);

      cs_b_cd_unsteady(bldfrp,
                       diipb[face_id],
                       grad[ii],
                       _pvar[ii],
                       &pip);

      cs_b_upwind_flux(iconvp,
                       thetap,
                       imasac,
                       inc,
                       bc_type[face_id],
                       _pvar[ii],
                       _pvar[ii], /* no relaxation */
                       pip,
                       coefap[face_id],
                       coefbp[face_id],
                       b_massflux[face_id],
                       1., /* xcpp */
                       &(b_conv_flux[face_id]));

    });
  }

  /* Boundary convective flux is imposed at some faces
     (tagged in icvfli array) */

  else if (icvflb) {

    /* Retrieve the value of the convective flux to be imposed */
    if (f_id != -1) {
      coface = f->bc_coeffs->ac;
      cofbce = f->bc_coeffs->bc;
    }
    else {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of icvflb and f_id"));
    }

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = b_face_cells[face_id];
      cs_real_t pip;

      cs_real_t bldfrp = (cs_real_t) ircflb;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflb > 0)
        bldfrp = cs::max(df_limiter[ii], 0.);

      cs_b_cd_unsteady(bldfrp,
                       diipb[face_id],
                       grad[ii],
                       _pvar[ii],
                       &pip);

      cs_b_imposed_conv_flux(iconvp,
                             thetap,
                             imasac,
                             inc,
                             bc_type[face_id],
                             icvfli[face_id],
                             _pvar[ii],
                             _pvar[ii], /* no relaxation */
                             pip,
                             coefap[face_id],
                             coefbp[face_id],
                             coface[face_id],
                             cofbce[face_id],
                             b_massflux[face_id],
                             1., /* xcpp */
                             &(b_conv_flux[face_id]));

    });
  }

  ctx.wait();

  /* Free memory */
  CS_FREE_HD(grad);
  CS_FREE_HD(gradup);
  CS_FREE_HD(gradst);
  CS_FREE_HD(local_max);
  CS_FREE_HD(local_min);
  CS_FREE_HD(courant);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the convection/diffusion terms of a transport
 *  equation of a vector field \f$ \vect{\varia} \f$.
 *
 * More precisely, the right hand side \f$ \vect{Rhs} \f$ is updated as
 * follows:
 * \f[
 *  \vect{Rhs} = \vect{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
 *         \dot{m}_\ij \left( \vect{\varia}_\fij - \vect{\varia}_\celli \right)
 *       - \mu_\fij \gradt_\fij \vect{\varia} \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * \deprecated
 *
 * Remark:
 * if ivisep = 1, then we also take \f$ \mu \transpose{\gradt\vect{\varia}}
 * + \lambda \trace{\gradt\vect{\varia}} \f$, where \f$ \lambda \f$ is
 * the secondary viscosity, i.e. usually \f$ -\frac{2}{3} \mu \f$.
 *
 * Warning:
 * - \f$ \vect{Rhs} \f$ has already been initialized before calling bilsc!
 * - mind the sign minus
 *
 * \param[in]     f             pointer to field, or nullptr
 * \param[in]     name          pointer to associated field or array name
 * \param[in]     eqp           equation parameters
 * \param[in]     cpl           structure associated with internal coupling
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     pvar          solved velocity (current time step)
 * \param[in]     pvara         solved velocity (previous time step)
 * \param[in]     icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in]     bc_coeffs_v   boundary conditions structure for the variable
 * \param[in]     bc_coeffs_solve_v   sweep loop boundary conditions structure
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     grad          associated gradient
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

static void
_convection_diffusion_vector_steady(cs_field_t                 *f,
                                    const char                 *var_name,
                                    const cs_equation_param_t  &eqp,
                                    int                         icvflb,
                                    int                         inc,
                                    cs_real_3_t                *pvar,
                                    const cs_real_3_t          *pvara,
                                    const int                   icvfli[],
                                    const cs_field_bc_coeffs_t *bc_coeffs_v,
                                    const cs_bc_coeffs_solve_t *bc_coeffs_solve_v,
                                    const cs_real_t             i_massflux[],
                                    const cs_real_t             b_massflux[],
                                    const cs_real_t             i_visc[],
                                    const cs_real_t             b_visc[],
                                    cs_real_33_t               *grad,
                                    cs_real_3_t       *restrict rhs)
{
  const cs_real_3_t  *coefav = (const cs_real_3_t  *)bc_coeffs_v->a;
  const cs_real_33_t *coefbv = (const cs_real_33_t *)bc_coeffs_v->b;
  const cs_real_3_t  *cofafv = (const cs_real_3_t  *)bc_coeffs_v->af;
  const cs_real_33_t *cofbfv = (const cs_real_33_t *)bc_coeffs_v->bf;

  const cs_real_3_t *val_f_g = (bc_coeffs_solve_v == nullptr) ?
                               (const cs_real_3_t *)bc_coeffs_v->val_f :
                               (const cs_real_3_t *)bc_coeffs_solve_v->val_f;

  const int iconvp = eqp.iconv;
  const int idiffp = eqp.idiff;
  const int ircflp = eqp.ircflu;
  const int ircflb = (ircflp > 0) ? eqp.b_diff_flux_rc : 0;
  const int ischcp = eqp.ischcv;
  const int isstpp = eqp.isstpc;
  const int iwarnp = eqp.verbosity;
  const double blencp = eqp.blencv;
  const double blend_st = eqp.blend_st;
  const double relaxp = eqp.relaxv;

  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_nreal_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;
  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_rreal_3_t *restrict diipf = fvq->diipf;
  const cs_rreal_3_t *restrict djjpf = fvq->djjpf;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;

  const int *bc_type = cs_glob_bc_type;
  cs_real_2_t *i_f_face_factor = nullptr;
  cs_real_t *b_f_face_factor = nullptr;

  /* Flux limiter */
  cs_real_t *df_limiter = nullptr;
  if (f != nullptr) {
    int df_limiter_id
      = cs_field_get_key_int(f, cs_field_key_id("diffusion_limiter_id"));
    if (df_limiter_id > -1)
      df_limiter = cs_field_by_id(df_limiter_id)->val;
  }

  /* Discontinuous porous treatment */
  if (cs_glob_porous_model == 3 && f == CS_F_(vel)) {
    i_f_face_factor = fvq->i_f_face_factor;
    b_f_face_factor = fvq->b_f_face_factor;
  }

  const cs_real_3_t *coface = nullptr;
  const cs_real_33_t *cofbce = nullptr;

  cs_dispatch_context ctx;
  ctx.set_use_gpu(false);  /* Steady case deprecated, so not ported to GPU */

  /*==========================================================================*/

  /* 1. Initialization */

  /* Allocate work arrays */

  cs_real_33_t *grdpa = nullptr;

  const cs_real_3_t  *restrict _pvar
    = (pvar != nullptr) ? (const cs_real_3_t  *)pvar : pvara;

  /* Slope limiters */

  if (iwarnp >= 2 && iconvp == 1) {
    if (ischcp == 1) {
      bft_printf
        (
         _(" %s: Convection in centered blending with %f percent of upwind\n"),
         var_name, (1.-blencp)*100.);
    }
    else {
      bft_printf
        (
         _(" %s: Convection in 2nd order blending with %f percent of upwind\n"),
         var_name, (1.-blencp)*100.);
    }
  }

  int pure_upwind = (blencp > 0.) ? 0 : 1;

  /* ======================================================================
     ---> Compute uncentered gradient grdpa for the slope test
     ======================================================================*/

  if (iconvp > 0 && pure_upwind == false && isstpp == 0) {
    CS_MALLOC_HD(grdpa, n_cells_ext, cs_real_33_t, cs_alloc_mode);

    _slope_test_gradient_strided<3>(ctx,
                                    (const cs_real_33_t *)grad,
                                    grdpa,
                                    _pvar,
                                    val_f_g,
                                    i_massflux);
    ctx.wait();
  }

  /* ======================================================================
     ---> Contribution from interior faces
     ======================================================================*/

  short *i_upwind = nullptr;
  if (iwarnp >= 2 && iconvp == 1) {
    cs_lnum_t n_i_faces = m->n_i_faces;
    CS_MALLOC_HD(i_upwind, n_i_faces, short, cs_alloc_mode);
    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      i_upwind[face_id] = 0;
    });
    ctx.wait();
  }

  if (n_cells_ext > n_cells) {
#   pragma omp parallel for if(n_cells_ext -n_cells > CS_THR_MIN)
    for (cs_lnum_t cell_id = n_cells; cell_id < n_cells_ext; cell_id++) {
      for (cs_lnum_t isou = 0; isou < 3; isou++)
        rhs[cell_id][isou] = 0.;
    }
  }

  /* Pure upwind flux
     ================ */

  if (pure_upwind == 1) {

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          /* in parallel, face will be counted by one and only one rank */
          if (i_upwind != nullptr && ii < n_cells) {
            i_upwind[face_id] = 1;
          }

          cs_real_t fluxi[3], fluxj[3] ;
          for (cs_lnum_t isou =  0; isou < 3; isou++) {
            fluxi[isou] = 0;
            fluxj[isou] = 0;
          }

          cs_real_3_t pip, pjp, pipr, pjpr;
          cs_real_3_t pifri, pifrj, pjfri, pjfrj;
          cs_real_3_t _pi, _pj, _pia, _pja;

          for (cs_lnum_t i = 0; i < 3; i++) {
            _pi[i]  = _pvar[ii][i];
            _pj[i]  = _pvar[jj][i];
            _pia[i] = pvara[ii][i];
            _pja[i] = pvara[jj][i];
          }

          /* Scaling due to mass balance in porous modelling */
          if (i_f_face_factor != nullptr) {
            const cs_nreal_t *n = i_face_u_normal[face_id];
            cs_math_3_normal_scaling(n, i_f_face_factor[face_id][0], _pi);
            cs_math_3_normal_scaling(n, i_f_face_factor[face_id][0], _pia);
            cs_math_3_normal_scaling(n, i_f_face_factor[face_id][1], _pj);
            cs_math_3_normal_scaling(n, i_f_face_factor[face_id][1], _pja);
          }

          cs_real_t bldfrp = (cs_real_t) ircflp;
          /* Local limitation of the reconstruction */
          if (df_limiter != nullptr && ircflp > 0)
            bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                             0.);

          cs_i_cd_steady_upwind_strided<3>(bldfrp,
                                           relaxp,
                                           diipf[face_id],
                                           djjpf[face_id],
                                           grad[ii],
                                           grad[jj],
                                           _pi,
                                           _pj,
                                           _pia,
                                           _pja,
                                           pifri,
                                           pifrj,
                                           pjfri,
                                           pjfrj,
                                           pip,
                                           pjp,
                                           pipr,
                                           pjpr);

          cs_i_conv_flux_strided<3>(iconvp,
                                    1.,
                                    1,
                                    _pvar[ii],
                                    _pvar[jj],
                                    pifri,
                                    pifrj,
                                    pjfri,
                                    pjfrj,
                                    i_massflux[face_id],
                                    fluxi,
                                    fluxj);

          cs_i_diff_flux_strided<3>(idiffp,
                                    1.,
                                    pip,
                                    pjp,
                                    pipr,
                                    pjpr,
                                    i_visc[face_id],
                                    fluxi,
                                    fluxj);

          for (cs_lnum_t isou = 0; isou < 3; isou++) {
            rhs[ii][isou] -= fluxi[isou];
            rhs[jj][isou] += fluxj[isou];
          }

        }
      }

    }

  }

  /* --> Flux with no slope test
     ============================*/

  else if (isstpp == 1) {

    if (ischcp < 0 || ischcp == 2 || ischcp > 3) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcv"));
    }

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          cs_real_t fluxi[3], fluxj[3] ;
          for (cs_lnum_t isou =  0; isou < 3; isou++) {
            fluxi[isou] = 0;
            fluxj[isou] = 0;
          }
          cs_real_3_t pip, pjp, pipr, pjpr;
          cs_real_3_t pifri, pifrj, pjfri, pjfrj;
          cs_real_3_t _pi, _pj, _pia, _pja;

          for (cs_lnum_t i = 0; i < 3; i++) {
            _pi[i]  = _pvar[ii][i];
            _pj[i]  = _pvar[jj][i];
            _pia[i] = pvara[ii][i];
            _pja[i] = pvara[jj][i];
          }

          /* Scaling due to mass balance in porous modelling */
          if (i_f_face_factor != nullptr) {
            const cs_nreal_t *n = i_face_u_normal[face_id];
            cs_math_3_normal_scaling(n, i_f_face_factor[face_id][0], _pi);
            cs_math_3_normal_scaling(n, i_f_face_factor[face_id][0], _pia);
            cs_math_3_normal_scaling(n, i_f_face_factor[face_id][1], _pj);
            cs_math_3_normal_scaling(n, i_f_face_factor[face_id][1], _pja);
          }

          cs_real_t bldfrp = (cs_real_t) ircflp;
          /* Local limitation of the reconstruction */
          if (df_limiter != nullptr && ircflp > 0)
            bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                             0.);

          cs_i_cd_steady_strided<3>(bldfrp,
                                    ischcp,
                                    relaxp,
                                    blencp,
                                    weight[face_id],
                                    cell_cen[ii],
                                    cell_cen[jj],
                                    i_face_cog[face_id],
                                    diipf[face_id],
                                    djjpf[face_id],
                                    grad[ii],
                                    grad[jj],
                                    _pi,
                                    _pj,
                                    _pia,
                                    _pja,
                                    pifri,
                                    pifrj,
                                    pjfri,
                                    pjfrj,
                                    pip,
                                    pjp,
                                    pipr,
                                    pjpr);

          cs_i_conv_flux_strided<3>(iconvp,
                                    1.,
                                    1,
                                    _pvar[ii],
                                    _pvar[jj],
                                    pifri,
                                    pifrj,
                                    pjfri,
                                    pjfrj,
                                    i_massflux[face_id],
                                    fluxi,
                                    fluxj);

          cs_i_diff_flux_strided<3>(idiffp,
                                    1.,
                                    pip,
                                    pjp,
                                    pipr,
                                    pjpr,
                                    i_visc[face_id],
                                    fluxi,
                                    fluxj);

          for (cs_lnum_t isou = 0; isou < 3; isou++) {
            rhs[ii][isou] -= fluxi[isou];
            rhs[jj][isou] += fluxj[isou];
          }

        }
      }
    }

  }

  /* Flux with slope test
     ==================== */

  else {

    if (ischcp < 0 || ischcp > 1) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcv"));
    }

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          cs_real_t fluxi[3], fluxj[3] ;
          for (cs_lnum_t isou =  0; isou < 3; isou++) {
            fluxi[isou] = 0;
            fluxj[isou] = 0;
          }
          cs_real_3_t pip, pjp, pipr, pjpr;
          cs_real_3_t pifri, pifrj, pjfri, pjfrj;
          bool upwind_switch = false;
          cs_real_3_t _pi, _pj, _pia, _pja;

          for (cs_lnum_t i = 0; i < 3; i++) {
            _pi[i]  = _pvar[ii][i];
            _pj[i]  = _pvar[jj][i];
            _pia[i] = pvara[ii][i];
            _pja[i] = pvara[jj][i];
          }

          /* Scaling due to mass balance in porous modelling */
          if (i_f_face_factor != nullptr) {
            const cs_nreal_t *n = i_face_u_normal[face_id];
            cs_math_3_normal_scaling(n, i_f_face_factor[face_id][0], _pi);
            cs_math_3_normal_scaling(n, i_f_face_factor[face_id][0], _pia);
            cs_math_3_normal_scaling(n, i_f_face_factor[face_id][1], _pj);
            cs_math_3_normal_scaling(n, i_f_face_factor[face_id][1], _pja);
          }

          cs_real_t bldfrp = (cs_real_t) ircflp;
          /* Local limitation of the reconstruction */
          if (df_limiter != nullptr && ircflp > 0)
            bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                             0.);

          cs_i_cd_steady_slope_test_strided<3>(&upwind_switch,
                                               iconvp,
                                               bldfrp,
                                               ischcp,
                                               relaxp,
                                               blencp,
                                               blend_st,
                                               weight[face_id],
                                               i_dist[face_id],
                                               cell_cen[ii],
                                               cell_cen[jj],
                                               i_face_u_normal[face_id],
                                               i_face_cog[face_id],
                                               diipf[face_id],
                                               djjpf[face_id],
                                               i_massflux[face_id],
                                               grad[ii],
                                               grad[jj],
                                               grdpa[ii],
                                               grdpa[jj],
                                               _pi,
                                               _pj,
                                               _pia,
                                               _pja,
                                               pifri,
                                               pifrj,
                                               pjfri,
                                               pjfrj,
                                               pip,
                                               pjp,
                                               pipr,
                                               pjpr);

          cs_i_conv_flux_strided<3>(iconvp,
                                    1.,
                                    1,
                                    _pvar[ii],
                                    _pvar[jj],
                                    pifri,
                                    pifrj,
                                    pjfri,
                                    pjfrj,
                                    i_massflux[face_id],
                                    fluxi,
                                    fluxj);

          cs_i_diff_flux_strided<3>(idiffp,
                                    1.,
                                    pip,
                                    pjp,
                                    pipr,
                                    pjpr,
                                    i_visc[face_id],
                                    fluxi,
                                    fluxj);

          for (cs_lnum_t isou = 0; isou < 3; isou++) {
            rhs[ii][isou] -= fluxi[isou];
            rhs[jj][isou] += fluxj[isou];
          }

        }
      }
    }

  } /* iupwin */

  if (iwarnp >= 2 && i_upwind != nullptr) {
    cs_gnum_t n_upwind = 0;
    const cs_lnum_t n_i_faces = m->n_i_faces;

    if (i_upwind != nullptr) {
      ctx.parallel_for_reduce_sum
        (n_i_faces, n_upwind, [=] CS_F_HOST_DEVICE
         (cs_lnum_t i,
          CS_DISPATCH_REDUCER_TYPE(cs_gnum_t) &sum) {
        sum += (cs_gnum_t)i_upwind[i];
      });

      ctx.wait();
      cs_parall_counter(&n_upwind, 1);
    }

    bft_printf(_(" %s: %llu Faces with upwind on %llu interior faces\n"),
               var_name, (unsigned long long)n_upwind,
               (unsigned long long)m->n_g_i_c_faces);

    CS_FREE_HD(i_upwind);
  }

  /* ======================================================================
     ---> Contribution from boundary faces
     ======================================================================*/

  /* Boundary convective flux are all computed with an upwind scheme */
  if (icvflb == 0) {

#   pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
    for (int t_id = 0; t_id < n_b_threads; t_id++) {
      for (cs_lnum_t face_id = b_group_index[t_id*2];
           face_id < b_group_index[t_id*2 + 1];
           face_id++) {

        cs_lnum_t ii = b_face_cells[face_id];

        cs_real_t fluxi[3];
        for (cs_lnum_t isou =  0; isou < 3; isou++) {
          fluxi[isou] = 0;
        }
        cs_real_3_t pir, pipr, val_f, val_f_d;
        cs_real_3_t _pi, _pia;
        for (int i = 0; i < 3; i++) {
          _pi[i]  = _pvar[ii][i];
          _pia[i] = pvara[ii][i];
        }

        /* Scaling due to mass balance in porous modelling */
        if (b_f_face_factor != nullptr) {
          const cs_nreal_t *n = i_face_u_normal[face_id];
          cs_math_3_normal_scaling(n, b_f_face_factor[face_id], _pi);
          cs_math_3_normal_scaling(n, b_f_face_factor[face_id], _pia);
        }

        cs_real_t bldfrp = (cs_real_t) ircflb;
        /* Local limitation of the reconstruction */
        if (df_limiter != nullptr && ircflb > 0)
          bldfrp = cs::max(df_limiter[ii], 0.);

        cs_b_cd_steady_strided<3>(bldfrp,
                                  relaxp,
                                  diipb[face_id],
                                  (const cs_real_3_t *)grad[ii],
                                  _pi,
                                  _pia,
                                  pir,
                                  pipr);

        /* Compute face value for gradient and diffusion for the
           steady case (relaxation value in iprime) */
        for (cs_lnum_t i = 0; i < 3; i++) {
          val_f[i] = inc * coefav[face_id][i];
          val_f_d[i] = inc * cofafv[face_id][i];

          for (cs_lnum_t j = 0; j < 3; j++) {
            val_f[i] += coefbv[face_id][j][i] * pipr[j];
            val_f_d[i] += cofbfv[face_id][j][i] * pipr[j];
          }
        }

        cs_b_upwind_flux_strided<3>(iconvp,
                                    1., /* thetap */
                                    1, /* imasac */
                                    bc_type[face_id],
                                    _pi,
                                    pir,
                                    b_massflux[face_id],
                                    val_f,
                                    fluxi);

        cs_b_diff_flux_strided<3>(idiffp,
                                  1., /* thetap */
                                  b_visc[face_id],
                                  val_f_d,
                                  fluxi);

        for (cs_lnum_t isou = 0; isou < 3; isou++) {
          rhs[ii][isou] -= fluxi[isou];
        } /* isou */

      }
    }

  }

  /* Boundary convective flux imposed at some faces (tags in icvfli array) */

  else if (icvflb == 1) {

    /* Retrieve the value of the convective flux to be imposed */
    if (f != nullptr) {
      coface = (const cs_real_3_t *)(f->bc_coeffs->ac);
      cofbce = (const cs_real_33_t *)(f->bc_coeffs->bc);
    }
    else {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of icvflb and f_id"));
    }

#   pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
    for (int t_id = 0; t_id < n_b_threads; t_id++) {
      for (cs_lnum_t face_id = b_group_index[t_id*2];
           face_id < b_group_index[t_id*2 + 1];
           face_id++) {

        cs_lnum_t ii = b_face_cells[face_id];

        cs_real_t fluxi[3], pir[3], pipr[3], val_f[3], val_f_d[3];

        for (cs_lnum_t isou =  0; isou < 3; isou++) {
          fluxi[isou] = 0;
        }
        cs_real_3_t _pi, _pia;

        for (cs_lnum_t i = 0; i < 3; i++) {
          _pi[i]  = _pvar[ii][i];
          _pia[i] = pvara[ii][i];
        }

        /* Scaling due to mass balance in porous modelling */
        if (b_f_face_factor != nullptr) {
          const cs_nreal_t *n = b_face_u_normal[face_id];
          cs_math_3_normal_scaling(n, b_f_face_factor[face_id], _pi);
          cs_math_3_normal_scaling(n, b_f_face_factor[face_id], _pia);
        }

        cs_real_t bldfrp = (cs_real_t) ircflb;
        /* Local limitation of the reconstruction */
        if (df_limiter != nullptr && ircflb > 0)
          bldfrp = cs::max(df_limiter[ii], 0.);

        cs_b_cd_steady_strided<3>(bldfrp,
                                  relaxp,
                                  diipb[face_id],
                                  (const cs_real_3_t *)grad[ii],
                                  _pi,
                                  _pia,
                                  pir,
                                  pipr);

        /* Compute face value for gradient and diffusion for the
           steady case (relaxation value in iprime) */
        for (cs_lnum_t i = 0; i < 3; i++) {
          val_f[i] = inc * coefav[face_id][i];
          val_f_d[i] = inc * cofafv[face_id][i];

          for (cs_lnum_t j = 0; j < 3; j++) {
            val_f[i] += coefbv[face_id][j][i] * pipr[j];
            val_f_d[i] += cofbfv[face_id][j][i] * pipr[j];
          }
        }

        cs_b_imposed_conv_flux_strided<3>(iconvp,
                                          1., /* thetap */
                                          1., /* imasac */
                                          inc,
                                          bc_type[face_id],
                                          icvfli[face_id],
                                          _pvar[ii],
                                          pir,
                                          pipr,
                                          coface[face_id],
                                          cofbce[face_id],
                                          b_massflux[face_id],
                                          val_f,
                                          fluxi);

        cs_b_diff_flux_strided<3>(idiffp,
                                  1., /* thetap */
                                  b_visc[face_id],
                                  val_f_d,
                                  fluxi);

        for (cs_lnum_t isou = 0; isou < 3; isou++) {
          rhs[ii][isou] -= fluxi[isou];
        }

      }
    }
  }

  /* Free memory */
  CS_FREE_HD(grdpa);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the convection/diffusion terms of a transport
 *  equation of a tensor field \f$ \tens{\varia} \f$.
 *
 * More precisely, the right hand side \f$ \tens{Rhs} \f$ is updated as
 * follows:
 * \f[
 *  \tens{Rhs} = \tens{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
 *         \dot{m}_\ij \left( \tens{\varia}_\fij - \tens{\varia}_\celli \right)
 *       - \mu_\fij \gradt_\fij \tens{\varia} \cdot \tens{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ \tens{Rhs} \f$ has already been initialized before calling bilsc!
 * - mind the sign minus
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          index of the current variable
 * \param[in]     eqp           equation parameters
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     pvar          solved velocity (current time step)
 * \param[in]     pvara         solved velocity (previous time step)
 * \param[in]     bc_coeffs_ts  boundary condition structure for the variable
 * \param[in]     bc_coeffs_solve_ts   sweep loop boundary conditions structure
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     grad          associated gradient
 * \param[in,out] rhs           right hand side \f$ \tens{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

static void
_convection_diffusion_tensor_steady
(
 cs_field_t                  *f,
 const char                  *var_name,
 const cs_equation_param_t   &eqp,
 int                          icvflb,
 int                          inc,
 cs_real_6_t                 *pvar,
 const cs_real_6_t           *pvara,
 const cs_field_bc_coeffs_t  *bc_coeffs_ts,
 const cs_bc_coeffs_solve_t  *bc_coeffs_solve_ts,
 const cs_real_t              i_massflux[],
 const cs_real_t              b_massflux[],
 const cs_real_t              i_visc[],
 const cs_real_t              b_visc[],
 cs_real_63_t                *grad,
 cs_real_6_t                 *rhs
)
{
  const cs_real_6_t  *coefa = (const cs_real_6_t  *)bc_coeffs_ts->a;
  const cs_real_66_t *coefb = (const cs_real_66_t *)bc_coeffs_ts->b;
  const cs_real_6_t  *cofaf = (const cs_real_6_t  *)bc_coeffs_ts->af;
  const cs_real_66_t *cofbf = (const cs_real_66_t *)bc_coeffs_ts->bf;

  const cs_real_6_t *val_f_g = (bc_coeffs_solve_ts == nullptr) ?
                               (const cs_real_6_t *)bc_coeffs_ts->val_f :
                               (const cs_real_6_t *)bc_coeffs_solve_ts->val_f;

  const int iconvp = eqp.iconv;
  const int idiffp = eqp.idiff;
  const int ircflp = eqp.ircflu;
  const int ircflb = (ircflp > 0) ? eqp.b_diff_flux_rc : 0;
  const int ischcp = eqp.ischcv;
  const int isstpp = eqp.isstpc;
  const int iwarnp = eqp.verbosity;
  const double blencp = eqp.blencv;
  const double blend_st = eqp.blend_st;
  const double relaxp = eqp.relaxv;

  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_nreal_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_rreal_3_t *restrict diipf = fvq->diipf;
  const cs_rreal_3_t *restrict djjpf = fvq->djjpf;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;

  const int *bc_type = cs_glob_bc_type;

  /* Flux limiter */
  cs_real_t *df_limiter = nullptr;
  if (f != nullptr) {
    int df_limiter_id
      = cs_field_get_key_int(f, cs_field_key_id("diffusion_limiter_id"));
    if (df_limiter_id > -1)
      df_limiter = cs_field_by_id(df_limiter_id)->val;
  }

  cs_dispatch_context ctx;
  ctx.set_use_gpu(false);  /* Steady case deprecated, so not ported to GPU */

  /*==========================================================================*/

  /* 1. Initialization */

  /* Allocate work arrays */

  cs_real_63_t *grdpa = nullptr;

  /* Handle cases where only the previous values (already synchronized)
     or current values are provided */

  const cs_real_6_t  *restrict _pvar
    = (pvar != nullptr) ? (const cs_real_6_t  *)pvar : pvara;

  /* Logging info */

  if (iwarnp >= 2 && iconvp == 1) {
    if (ischcp == 1) {
      bft_printf
        (
         _(" %s: Convection in centered blending with %f percent of upwind\n"),
         var_name, (1.-blencp)*100.);
    }
    else {
      bft_printf
        (
         _(" %s: Convection in 2nd order blending with %f percent of upwind\n"),
         var_name, (1.-blencp)*100.);
    }
  }

  int iupwin = (blencp > 0.) ? 0 : 1;

  /* ======================================================================
     ---> Compute uncentered gradient grdpa for the slope test
     ======================================================================*/

  if (iconvp > 0 && iupwin == 0 && isstpp == 0) {
    CS_MALLOC(grdpa, n_cells_ext, cs_real_63_t);

    _slope_test_gradient_strided<6>(ctx,
                                    (const cs_real_63_t *)grad,
                                    grdpa,
                                    _pvar,
                                    val_f_g,
                                    i_massflux);
    ctx.wait();
  }

  /* ======================================================================
     ---> Contribution from interior faces
     ======================================================================*/

  cs_gnum_t n_upwind = 0;

  if (n_cells_ext > n_cells) {
#   pragma omp parallel for
    for (cs_lnum_t cell_id = n_cells; cell_id < n_cells_ext; cell_id++) {
      for (cs_lnum_t isou = 0; isou < 6; isou++)
        rhs[cell_id][isou] = 0.;
    }
  }

  /* --> Pure upwind flux
     =====================*/

  if (iupwin == 1) {

    n_upwind = m->n_g_i_faces;

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          cs_real_t fluxi[6], fluxj[6] ;
          for (cs_lnum_t isou =  0; isou < 6; isou++) {
            fluxi[isou] = 0;
            fluxj[isou] = 0;
          }
          cs_real_t pip[6], pjp[6], pipr[6], pjpr[6];
          cs_real_t pifri[6], pifrj[6], pjfri[6], pjfrj[6];

          cs_real_t bldfrp = (cs_real_t) ircflp;
          /* Local limitation of the reconstruction */
          if (df_limiter != nullptr && ircflp > 0)
            bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                             0.);

          cs_i_cd_steady_upwind_strided<6>(bldfrp,
                                           relaxp,
                                           diipf[face_id],
                                           djjpf[face_id],
                                           (const cs_real_3_t *)grad[ii],
                                           (const cs_real_3_t *)grad[jj],
                                           _pvar[ii],
                                           _pvar[jj],
                                           pvara[ii],
                                           pvara[jj],
                                           pifri,
                                           pifrj,
                                           pjfri,
                                           pjfrj,
                                           pip,
                                           pjp,
                                           pipr,
                                           pjpr);


          cs_i_conv_flux_strided<6>(iconvp,
                                    1.,
                                    1,
                                    _pvar[ii],
                                    _pvar[jj],
                                    pifri,
                                    pifrj,
                                    pjfri,
                                    pjfrj,
                                    i_massflux[face_id],
                                    fluxi,
                                    fluxj);

          cs_i_diff_flux_strided<6>(idiffp,
                                    1.,
                                    pip,
                                    pjp,
                                    pipr,
                                    pjpr,
                                    i_visc[face_id],
                                    fluxi,
                                    fluxj);

          for (cs_lnum_t isou = 0; isou < 6; isou++) {
            rhs[ii][isou] -= fluxi[isou];
            rhs[jj][isou] += fluxj[isou];
          }

        }
      }
    }

  }

  /* --> Flux with no slope test
     ============================*/

  else if (isstpp == 1) {

    if (ischcp < 0 || ischcp > 1) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcp"));
    }

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          double fluxi[6], fluxj[6] ;
          for (cs_lnum_t isou =  0; isou < 6; isou++) {
            fluxi[isou] = 0;
            fluxj[isou] = 0;
          }
          cs_real_t pip[6], pjp[6], pipr[6], pjpr[6];
          cs_real_t pifri[6], pifrj[6], pjfri[6], pjfrj[6];

          cs_real_t bldfrp = (cs_real_t) ircflp;
          /* Local limitation of the reconstruction */
          if (df_limiter != nullptr && ircflp > 0)
            bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                             0.);

          cs_i_cd_steady_strided<6>(bldfrp,
                                    ischcp,
                                    relaxp,
                                    blencp,
                                    weight[face_id],
                                    cell_cen[ii],
                                    cell_cen[jj],
                                    i_face_cog[face_id],
                                    diipf[face_id],
                                    djjpf[face_id],
                                    grad[ii],
                                    grad[jj],
                                    _pvar[ii],
                                    _pvar[jj],
                                    pvara[ii],
                                    pvara[jj],
                                    pifri,
                                    pifrj,
                                    pjfri,
                                    pjfrj,
                                    pip,
                                    pjp,
                                    pipr,
                                    pjpr);

          cs_i_conv_flux_strided<6>(iconvp,
                                    1.,
                                    1,
                                    _pvar[ii],
                                    _pvar[jj],
                                    pifri,
                                    pifrj,
                                    pjfri,
                                    pjfrj,
                                    i_massflux[face_id],
                                    fluxi,
                                    fluxj);

          cs_i_diff_flux_strided<6>(idiffp,
                                    1.,
                                    pip,
                                    pjp,
                                    pipr,
                                    pjpr,
                                    i_visc[face_id],
                                    fluxi,
                                    fluxj);

          for (cs_lnum_t isou = 0; isou < 6; isou++) {
            rhs[ii][isou] -= fluxi[isou];
            rhs[jj][isou] += fluxj[isou];
          }

        }
      }

    }

  }

  /* --> Flux with slope test
     =========================*/

  else {

    if (ischcp < 0 || ischcp > 1) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcp"));
    }

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for reduction(+:n_upwind)
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          double fluxi[6], fluxj[6] ;
          for (cs_lnum_t isou =  0; isou < 6; isou++) {
            fluxi[isou] = 0;
            fluxj[isou] = 0;
          }
          cs_real_t pip[6], pjp[6], pipr[6], pjpr[6];
          cs_real_t pifri[6], pifrj[6], pjfri[6], pjfrj[6];
          bool upwind_switch = false;

          cs_real_t bldfrp = (cs_real_t) ircflp;
          /* Local limitation of the reconstruction */
          if (df_limiter != nullptr && ircflp > 0)
            bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                             0.);

          cs_i_cd_steady_slope_test_strided<6>(&upwind_switch,
                                               iconvp,
                                               bldfrp,
                                               ischcp,
                                               relaxp,
                                               blencp,
                                               blend_st,
                                               weight[face_id],
                                               i_dist[face_id],
                                               cell_cen[ii],
                                               cell_cen[jj],
                                               i_face_u_normal[face_id],
                                               i_face_cog[face_id],
                                               diipf[face_id],
                                               djjpf[face_id],
                                               i_massflux[face_id],
                                               grad[ii],
                                               grad[jj],
                                               grdpa[ii],
                                               grdpa[jj],
                                               _pvar[ii],
                                               _pvar[jj],
                                               pvara[ii],
                                               pvara[jj],
                                               pifri,
                                               pifrj,
                                               pjfri,
                                               pjfrj,
                                               pip,
                                               pjp,
                                               pipr,
                                               pjpr);

          cs_i_conv_flux_strided<6>(iconvp,
                                    1.,
                                    1,
                                    _pvar[ii],
                                    _pvar[jj],
                                    pifri,
                                    pifrj,
                                    pjfri,
                                    pjfrj,
                                    i_massflux[face_id],
                                    fluxi,
                                    fluxj);

          cs_i_diff_flux_strided<6>(idiffp,
                                    1.,
                                    pip,
                                    pjp,
                                    pipr,
                                    pjpr,
                                    i_visc[face_id],
                                    fluxi,
                                    fluxj);

          for (cs_lnum_t isou = 0; isou < 6; isou++) {
            rhs[ii][isou] -= fluxi[isou];
            rhs[jj][isou] += fluxj[isou];
          }

        }
      }
    }

    if (iwarnp >= 2 && iconvp == 1)
      cs_parall_counter(&n_upwind, 1);

  } /* iupwin */

  if (iwarnp >= 2 && iconvp == 1)
    bft_printf(_(" %s: %llu Faces with upwind on %llu interior faces\n"),
               var_name, (unsigned long long)n_upwind,
               (unsigned long long)m->n_g_i_c_faces);

  /* ======================================================================
     Contribution from boundary faces
     ======================================================================*/

  /* Boundary convective flux are all computed with an upwind scheme */
  if (icvflb == 0) {

#   pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
    for (int t_id = 0; t_id < n_b_threads; t_id++) {
      for (cs_lnum_t face_id = b_group_index[t_id*2];
           face_id < b_group_index[t_id*2 + 1];
           face_id++) {

        cs_lnum_t ii = b_face_cells[face_id];

        double fluxi[6] ;
        for (cs_lnum_t isou =  0; isou < 6; isou++) {
          fluxi[isou] = 0;
        }
        cs_real_t pir[6], pipr[6], val_f[6], val_f_d[6];

        cs_real_t bldfrp = (cs_real_t) ircflb;
        /* Local limitation of the reconstruction */
        if (df_limiter != nullptr && ircflb > 0)
          bldfrp = cs::max(df_limiter[ii], 0.);

        cs_b_cd_steady_strided<6>(bldfrp,
                                  relaxp,
                                  diipb[face_id],
                                  (const cs_real_3_t *)grad[ii],
                                  _pvar[ii],
                                  pvara[ii],
                                  pir,
                                  pipr);

        /* Compute face value for gradient and diffusion for the
           steady case (relaxation value in iprime) */
        for (cs_lnum_t i = 0; i < 6; i++) {
          val_f[i] = inc * coefa[face_id][i];
          val_f_d[i] = inc * cofaf[face_id][i];

          for (cs_lnum_t j = 0; j < 6; j++) {
            val_f[i] += coefb[face_id][j][i] * pipr[j];
            val_f_d[i] += cofbf[face_id][j][i] * pipr[j];
          }
        }

        cs_b_upwind_flux_strided<6>(iconvp,
                                    1., /* thetap */
                                    1, /* imasac */
                                    bc_type[face_id],
                                    _pvar[ii],
                                    pir,
                                    b_massflux[face_id],
                                    val_f,
                                    fluxi);

        cs_b_diff_flux_strided<6>(idiffp,
                                  1., /* thetap */
                                  b_visc[face_id],
                                  val_f_d,
                                  fluxi);

        for (cs_lnum_t isou = 0; isou < 6; isou++) {
          rhs[ii][isou] -= fluxi[isou];
        }
      }
    }

  }

  /* Free memory */
  CS_FREE(grdpa);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the convection/diffusion terms of a transport
 *  equation of a vector field \f$ \vect{\varia} \f$.
 *
 * More precisely, the right hand side \f$ \vect{Rhs} \f$ is updated as
 * follows:
 * \f[
 *  \vect{Rhs} = \vect{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
 *         \dot{m}_\ij \left( \vect{\varia}_\fij - \vect{\varia}_\celli \right)
 *       - \mu_\fij \gradt_\fij \vect{\varia} \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Remark:
 * if ivisep = 1, then we also take \f$ \mu \transpose{\gradt\vect{\varia}}
 * + \lambda \trace{\gradt\vect{\varia}} \f$, where \f$ \lambda \f$ is
 * the secondary viscosity, i.e. usually \f$ -\frac{2}{3} \mu \f$.
 *
 * Warning:
 * - \f$ \vect{Rhs} \f$ has already been initialized before calling bilsc!
 * - mind the sign minus
 *
 * \param[in]     context       reference to dispatch context
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f             pointer to current field, or nullptr
 * \param[in]     name          pointer to associated field or array name
 * \param[in]     eqp           equation parameters
 * \param[in]     cpl           structure associated with internal coupling
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imasac        take mass accumulation into account?
 * \param[in]     pvar          solved velocity (current time step)
 * \param[in]     pvara         solved velocity (previous time step)
 * \param[in]     icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in]     bc_coeffs   boundary conditions structure for the variable
 * \param[in]     bc_coeffs_solve   sweep loop boundary conditions structure
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     i_pvar        velocity at interior faces
 * \param[in]     b_pvar        velocity at boundary faces
 * \param[in]     grad          associated gradient
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride, bool porous_vel>
static void
_convection_diffusion_unsteady_strided
  (cs_dispatch_context         &ctx,
   const cs_field_t            *f,
   const char                  *var_name,
   const cs_equation_param_t   &eqp,
   int                          icvflb,
   int                          inc,
   int                          imasac,
   cs_real_t                  (*pvar)[stride],
   const cs_real_t            (*pvara)[stride],
   const int                    icvfli[],
   const cs_field_bc_coeffs_t  *bc_coeffs,
   const cs_bc_coeffs_solve_t  *bc_coeffs_solve,
   const cs_real_t              i_massflux[],
   const cs_real_t              b_massflux[],
   const cs_real_t              i_visc[],
   const cs_real_t              b_visc[],
   cs_real_t         (*restrict i_pvar)[stride],
   cs_real_t         (*restrict b_pvar)[stride],
   cs_real_t         (*restrict grad)[stride][3],
   cs_real_t         (*restrict rhs)[stride])
{
  using grad_t = cs_real_t[stride][3];
  using var_t = cs_real_t[stride];
  using b_t = cs_real_t[stride][stride];

  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  const var_t *val_f = (bc_coeffs_solve == nullptr) ?
                       (const var_t *)bc_coeffs->val_f :
                       (const var_t *)bc_coeffs_solve->val_f;

  const var_t *val_f_d_lim = (bc_coeffs_solve == nullptr) ?
                             (const var_t *)bc_coeffs->val_f_d_lim :
                             (const var_t *)bc_coeffs_solve->val_f_d_lim;

  const int iconvp = eqp.iconv;
  const int idiffp = eqp.idiff;
  const int ircflp = eqp.ircflu;
  const int ircflb = (ircflp > 0) ? eqp.b_diff_flux_rc : 0;
  const int ischcp = eqp.ischcv;
  const int isstpp = eqp.isstpc;
  const double blencp = eqp.blencv;
  const double blend_st = eqp.blend_st;
  const double thetap = eqp.theta;

  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_halo_t  *halo = m->halo;
  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_nreal_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;
  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_rreal_3_t *restrict diipf = fvq->diipf;
  const cs_rreal_3_t *restrict djjpf = fvq->djjpf;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;

  const int *restrict c_disable_flag = (fvq->has_disable_flag) ?
    fvq->c_disable_flag : nullptr;

  const int *bc_type = cs_glob_bc_type;
  cs_real_2_t *i_f_face_factor = nullptr;
  cs_real_t *b_f_face_factor = nullptr;

  /* Local variables */

  cs_real_t *df_limiter = nullptr;

  if (f != nullptr) {
    int df_limiter_id
      = cs_field_get_key_int(f, cs_field_key_id("diffusion_limiter_id"));
    if (df_limiter_id > -1)
      df_limiter = cs_field_by_id(df_limiter_id)->val;
  }

  /* Discontinuous porous treatment */
  if (porous_vel) {
    i_f_face_factor = fvq->i_f_face_factor;
    b_f_face_factor = fvq->b_f_face_factor;
  }

  const int f_id = (f != nullptr) ? f->id : -1;

  const var_t *coface = nullptr;
  const b_t *cofbce = nullptr;

  cs_real_t  *v_slope_test = cs_get_v_slope_test(f_id,  eqp);

  /* Parallel or device dispatch */
  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  /*==========================================================================*/

  /* Initialization */

  cs_alloc_mode_t amode = ctx.alloc_mode(true);

  /* Handle cases where only the previous values (already synchronized)
     or current values are provided */

  if (pvar != nullptr && halo != nullptr) {
    cs_halo_sync_r(m->halo, CS_HALO_STANDARD, ctx.use_gpu(), pvar);
  }
  if (pvara == nullptr)
    pvara = (const var_t *)pvar;

  const var_t *_pvar = (pvar != nullptr) ? (const var_t *)pvar : pvara;

  /* Slope limiters */

  if (eqp.verbosity >= 2 && eqp.iconv == 1) {
    if (ischcp == 1) {
      bft_printf
        (
         _(" %s: Convection in centered blending with %f percent of upwind\n"),
         var_name, (1.-blencp)*100.);
    }
    else {
      bft_printf
        (
         _(" %s: Convection in 2nd order blending with %f percent of upwind\n"),
         var_name, (1.-blencp)*100.);
    }
  }

  bool pure_upwind = (blencp > 0.) ? false : true;

  /* Compute the balance with reconstruction */

  /* ======================================================================
     Compute uncentered gradient grdpa for the slope test
     ======================================================================*/

  grad_t *grdpa = nullptr;

  ctx.wait();

  if (iconvp > 0 && pure_upwind == false && isstpp == 0) {
    CS_MALLOC_HD(grdpa, n_cells_ext, grad_t, amode);

    _slope_test_gradient_strided<stride>(ctx,
                                         (const grad_t *)grad,
                                         grdpa,
                                         _pvar,
                                         val_f,
                                         i_massflux);
  }

  /* ======================================================================
     Contribution from interior faces
     ======================================================================*/

  short *i_upwind = nullptr;
  if (eqp.verbosity >= 2 && iconvp == 1 && pure_upwind == false) {
    CS_MALLOC_HD(i_upwind, n_i_faces, short, cs_alloc_mode);
    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      i_upwind[face_id] = 0;
    });
    ctx.wait();
  }

  /* Pure upwind flux
     ================ */

  if (pure_upwind) {

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t fluxi[stride], fluxj[stride];
      for (cs_lnum_t i = 0; i < stride; i++) {
        fluxi[i] = 0;
        fluxj[i] = 0;
      }
      cs_real_t pip[stride], pjp[stride];
      cs_real_t _pi[stride], _pj[stride];

      const var_t &pvar_i = _pvar[ii];
      const var_t &pvar_j = _pvar[jj];

      for (cs_lnum_t i = 0; i < stride; i++) {
        _pi[i]  = _pvar[ii][i];
        _pj[i]  = _pvar[jj][i];
      }

      /* Scaling due to mass balance in porous modelling */
      if (porous_vel) {
        const cs_nreal_t *n = i_face_u_normal[face_id];
        cs_math_3_normal_scaling(n, i_f_face_factor[face_id][0], _pi);
        cs_math_3_normal_scaling(n, i_f_face_factor[face_id][1], _pj);
      }

      if (ircflp == 1) {
        cs_real_t bldfrp = 1.;
        if (df_limiter != nullptr)  /* Local limiter */
          bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                           0.);

        cs_real_t recoi[stride], recoj[stride];
        cs_i_compute_quantities_strided<stride>(bldfrp,
                                                diipf[face_id], djjpf[face_id],
                                                grad[ii], grad[jj],
                                                _pi, _pj,
                                                recoi, recoj,
                                                pip, pjp);
      }
      else {
        for (cs_lnum_t isou = 0; isou < stride; isou++) {
          pip[isou] = _pi[isou];
          pjp[isou] = _pj[isou];
        }
      }

      /* No relaxation in following expressions */

      // Convective flux (pif == _pi, pjf == _pj as we are upwind)

      if (iconvp) {
        cs_real_t _i_massflux = i_massflux[face_id];
        cs_real_t flui = 0.5*(_i_massflux + cs::abs(_i_massflux));
        cs_real_t fluj = 0.5*(_i_massflux - cs::abs(_i_massflux));

        for (cs_lnum_t isou = 0; isou < stride; isou++) {
          fluxi[isou] +=   thetap*(flui*_pi[isou] + fluj*_pj[isou])
                         - imasac*_i_massflux*pvar_i[isou];
          fluxj[isou] +=   thetap*(flui*_pi[isou] + fluj*_pj[isou])
                         - imasac*_i_massflux*pvar_j[isou];
        }
      }

      // Diffusive flux

      cs_real_t t_i_visc = idiffp * i_visc[face_id] * thetap;
      for (cs_lnum_t isou = 0; isou < stride; isou++) {
        fluxi[isou] += t_i_visc*(pip[isou] -pjp[isou]);
        fluxj[isou] += t_i_visc*(pip[isou] -pjp[isou]);
      }

      /* Saving values at faces, if needed */
      if (i_pvar != nullptr) {
        if (i_massflux[face_id] >= 0.) {
          for (cs_lnum_t i = 0; i < stride; i++)
            i_pvar[face_id][i] += thetap * _pi[i];
        }
        else {
          for (cs_lnum_t i = 0; i < stride; i++)
            i_pvar[face_id][i] += thetap * _pj[i];
        }
      }

      for (cs_lnum_t isou = 0; isou < stride; isou++)
        fluxi[isou] *= -1;  /* For substraction in dispatch sum */

      if (ii < n_cells)
        cs_dispatch_sum<stride>(rhs[ii], fluxi, i_sum_type);
      if (jj < n_cells)
        cs_dispatch_sum<stride>(rhs[jj], fluxj, i_sum_type);
    });
  }

  /* Flux with no slope test
     ======================= */

  else if (isstpp == 1) {

    const cs_real_t *hybrid_blend = nullptr;
    if (CS_F_(hybrid_blend) != nullptr)
      hybrid_blend = CS_F_(hybrid_blend)->val;

    switch(ischcp) {
    case 0:
      [[fallthrough]];
    case 1:
      [[fallthrough]];
    case 3:
      ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

        cs_lnum_t ii = i_face_cells[face_id][0];
        cs_lnum_t jj = i_face_cells[face_id][1];

        cs_real_t fluxi[stride], fluxj[stride] ;
        for (cs_lnum_t isou =  0; isou < stride; isou++) {
          fluxi[isou] = 0;
          fluxj[isou] = 0;
        }

        cs_real_t pip[stride], pjp[stride];
        cs_real_t pif[stride], pjf[stride];
        cs_real_t _pi[stride], _pj[stride];

        const var_t &pvar_i = _pvar[ii];
        const var_t &pvar_j = _pvar[jj];

        for (cs_lnum_t i = 0; i < stride; i++) {
          _pi[i]  = _pvar[ii][i];
          _pj[i]  = _pvar[jj][i];
        }

        /* Scaling due to mass balance in porous modelling */
        if (porous_vel) {
          const cs_nreal_t *n = i_face_u_normal[face_id];
          cs_math_3_normal_scaling(n, i_f_face_factor[face_id][0], _pi);
          cs_math_3_normal_scaling(n, i_f_face_factor[face_id][1], _pj);
        }

        if (ircflp == 1) {
          cs_real_t bldfrp = 1.;
          if (df_limiter != nullptr)  /* Local limiter */
            bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                             0.);

          cs_real_t recoi[stride], recoj[stride];
          cs_i_compute_quantities_strided<stride>(bldfrp,
                                                  diipf[face_id], djjpf[face_id],
                                                  grad[ii], grad[jj],
                                                  _pi, _pj,
                                                  recoi, recoj,
                                                  pip, pjp);
        }
        else {
          for (cs_lnum_t isou = 0; isou < stride; isou++) {
            pip[isou] = _pi[isou];
            pjp[isou] = _pj[isou];
          }
        }

        const cs_real_t *cell_ceni = cell_cen[ii];
        const cs_real_t *cell_cenj = cell_cen[jj];
        const cs_real_t w_f = weight[face_id];

        if (ischcp == 1) {

          /* Centered
             --------*/

          for (cs_lnum_t isou = 0; isou < stride; isou++) {
            pif[isou] = w_f*pip[isou] + (1.-w_f)*pjp[isou];
            pjf[isou] = pif[isou];
          }

        }
        else if (ischcp == 3) {

          cs_real_t hybrid_blend_interp
            = cs::min(hybrid_blend[ii], hybrid_blend[jj]);

          /* Centered
             -------- */

          for (cs_lnum_t isou = 0; isou < stride; isou++) {
            pif[isou] = w_f*pip[isou] + (1.-w_f)*pjp[isou];
            pjf[isou] = pif[isou];
          }

          /* SOLU
             -----*/
          cs_real_t pif_up[stride], pjf_up[stride];

          cs_solu_f_val_strided<stride>(cell_ceni,
                                        i_face_cog[face_id],
                                        grad[ii],
                                        _pi,
                                        pif_up);
          cs_solu_f_val_strided<stride>(cell_cenj,
                                        i_face_cog[face_id],
                                        grad[jj],
                                        _pj,
                                        pjf_up);

          for (cs_lnum_t isou = 0; isou < stride; isou++) {
            pif[isou] =         hybrid_blend_interp *pif[isou]
                        + (1. - hybrid_blend_interp)*pif_up[isou];
            pjf[isou] =         hybrid_blend_interp *pjf[isou]
                        + (1. - hybrid_blend_interp)*pjf_up[isou];
          }
        }
        else {

          /* Second order
             ------------*/

          cs_solu_f_val_strided<stride>(cell_ceni,
                                        i_face_cog[face_id],
                                        grad[ii],
                                        _pi,
                                        pif);
          cs_solu_f_val_strided<stride>(cell_cenj,
                                        i_face_cog[face_id],
                                        grad[jj],
                                        _pj,
                                        pjf);

        }

        /* Blending
           --------*/

        for (cs_lnum_t isou = 0; isou < stride; isou++) {
          pif[isou] = blencp*(pif[isou])+(1.-blencp)*_pi[isou];
          pjf[isou] = blencp*(pjf[isou])+(1.-blencp)*_pj[isou];
        }

        // Convective flux

        if (iconvp) {
          cs_real_t _i_massflux = i_massflux[face_id];
          cs_real_t flui = 0.5*(_i_massflux + cs::abs(_i_massflux));
          cs_real_t fluj = 0.5*(_i_massflux - cs::abs(_i_massflux));

          for (cs_lnum_t isou = 0; isou < stride; isou++) {
            fluxi[isou] +=   thetap*(flui*pif[isou] + fluj*pjf[isou])
                           - imasac*_i_massflux*pvar_i[isou];
            fluxj[isou] +=   thetap*(flui*pif[isou] + fluj*pjf[isou])
                           - imasac*_i_massflux*pvar_j[isou];
          }
        }

        // Diffusive flux

        cs_real_t t_i_visc = idiffp * i_visc[face_id] * thetap;
        for (cs_lnum_t isou = 0; isou < stride; isou++) {
          fluxi[isou] += t_i_visc*(pip[isou] -pjp[isou]);
          fluxj[isou] += t_i_visc*(pip[isou] -pjp[isou]);
        }

        /* Save values at internal faces, if needed */
        if (i_pvar != nullptr) {
          if (i_massflux[face_id] >= 0.) {
            for (cs_lnum_t i = 0; i < stride; i++)
              i_pvar[face_id][i] += thetap * pif[i];
          }
          else {
            for (cs_lnum_t i = 0; i < stride; i++)
              i_pvar[face_id][i] += thetap * pjf[i];
          }
        }

        for (cs_lnum_t isou = 0; isou < stride; isou++)
          fluxi[isou] *= -1;  /* For substraction in dispatch sum */

        if (ii < n_cells)
          cs_dispatch_sum<stride>(rhs[ii], fluxi, i_sum_type);
        if (jj < n_cells)
          cs_dispatch_sum<stride>(rhs[jj], fluxj, i_sum_type);
      });
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _("Invalid value of ischcv with isstpp == 1:\n"
                  "0, 1 or 3 expected, %d here."), ischcp);
    }
  }

  /* Flux with slope test
     ==================== */

  else {

    switch(ischcp) {
    case 0:
      [[fallthrough]];
    case 1:
      ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

        cs_lnum_t ii = i_face_cells[face_id][0];
        cs_lnum_t jj = i_face_cells[face_id][1];

        cs_real_t fluxi[stride], fluxj[stride] ;
        for (cs_lnum_t isou =  0; isou < stride; isou++) {
          fluxi[isou] = 0;
          fluxj[isou] = 0;
        }
        cs_real_t pip[stride], pjp[stride];
        cs_real_t pif[stride], pjf[stride];
        bool upwind_switch = false;
        cs_real_t _pi[stride], _pj[stride];

        const var_t &pvar_i = _pvar[ii];
        const var_t &pvar_j = _pvar[jj];

        for (cs_lnum_t i = 0; i < stride; i++) {
          _pi[i]  = _pvar[ii][i];
          _pj[i]  = _pvar[jj][i];
        }

        /* Scaling due to mass balance in porous modelling */
        if (porous_vel) {
          const cs_nreal_t *n = i_face_u_normal[face_id];
          cs_math_3_normal_scaling(n, i_f_face_factor[face_id][0], _pi);
          cs_math_3_normal_scaling(n, i_f_face_factor[face_id][1], _pj);
        }

        if (ircflp == 1) {
          cs_real_t bldfrp = 1.;
          if (df_limiter != nullptr)  /* Local limiter */
            bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                             0.);

          cs_real_t recoi[stride], recoj[stride];
          cs_i_compute_quantities_strided<stride>(bldfrp,
                                                  diipf[face_id], djjpf[face_id],
                                                  grad[ii], grad[jj],
                                                  _pi, _pj,
                                                  recoi, recoj,
                                                  pip, pjp);
        }
        else {
          for (cs_lnum_t isou = 0; isou < stride; isou++) {
            pip[isou] = _pi[isou];
            pjp[isou] = _pj[isou];
          }
        }

        const cs_real_t *cell_ceni = cell_cen[ii];
        const cs_real_t *cell_cenj = cell_cen[jj];
        const cs_real_t w_f = weight[face_id];

        /* Slope test is needed with convection */
        if (iconvp > 0) {
          cs_real_t testij, tesqck;

          const grad_t &gradi = grad[ii];
          const grad_t &gradj = grad[jj];

          cs_slope_test_strided<stride>(_pi, _pj,
                                        i_dist[face_id],
                                        i_face_u_normal[face_id],
                                        gradi, gradj,
                                        grdpa[ii], grdpa[jj],
                                        i_massflux[face_id],
                                        &testij, &tesqck);

          for (cs_lnum_t isou = 0; isou < stride; isou++) {

            if (ischcp == 1) {

              /* Centered
                 --------*/

              cs_centered_f_val(w_f,
                                pip[isou], pjp[isou],
                                &pif[isou]);
              cs_centered_f_val(w_f,
                                pip[isou], pjp[isou],
                                &pjf[isou]);

            }
            else {

              /* Second order
                 ------------*/

              cs_solu_f_val(cell_ceni,
                            i_face_cog[face_id],
                            gradi[isou],
                            _pi[isou],
                            &pif[isou]);
              cs_solu_f_val(cell_cenj,
                            i_face_cog[face_id],
                            gradj[isou],
                            _pj[isou],
                            &pjf[isou]);
            }

          }

          /* Slope test activated: percentage of upwind */
          if (tesqck <= 0. || testij <= 0.) {

            /* Upwind
               ------ */

            cs_blend_f_val_strided<stride>(blend_st, _pi, pif);
            cs_blend_f_val_strided<stride>(blend_st, _pj, pjf);

            upwind_switch = true;
          }

          /* Blending
             -------- */

          cs_blend_f_val_strided<stride>(blencp, _pi, pif);
          cs_blend_f_val_strided<stride>(blencp, _pj, pjf);

        }
        else { /* If iconv=0 p*fr* are useless */

          for (cs_lnum_t isou = 0; isou < stride; isou++) {
            pif[isou] = _pi[isou];
            pjf[isou] = _pj[isou];
          }

        } /* End for slope test */

        // Convective flux

        if (iconvp) {
          cs_real_t _i_massflux = i_massflux[face_id];
          cs_real_t flui = 0.5*(_i_massflux + cs::abs(_i_massflux));
          cs_real_t fluj = 0.5*(_i_massflux - cs::abs(_i_massflux));

          for (cs_lnum_t isou = 0; isou < stride; isou++) {
            fluxi[isou] +=   thetap*(flui*pif[isou] + fluj*pjf[isou])
                           - imasac*_i_massflux*pvar_i[isou];
            fluxj[isou] +=   thetap*(flui*pif[isou] + fluj*pjf[isou])
                           - imasac*_i_massflux*pvar_j[isou];
          }
        }

        // Diffusive flux

        cs_real_t t_i_visc = idiffp * i_visc[face_id] * thetap;
        for (cs_lnum_t isou = 0; isou < stride; isou++) {
          fluxi[isou] += t_i_visc*(pip[isou] -pjp[isou]);
          fluxj[isou] += t_i_visc*(pip[isou] -pjp[isou]);
        }

        if (upwind_switch) {
          /* in parallel, face will be counted by one and only one rank */
          if (i_upwind != nullptr && ii < n_cells) {
            i_upwind[face_id] = 1;
          }

          if (v_slope_test != nullptr) {
            cs_real_t q_d_vol_ii
              =   std::abs(i_massflux[face_id])
                * cs_mq_cell_vol_inv(ii, c_disable_flag, cell_vol);
            cs_real_t q_d_vol_jj
              =   std::abs(i_massflux[face_id])
                * cs_mq_cell_vol_inv(jj, c_disable_flag, cell_vol);

            cs_dispatch_sum(&v_slope_test[ii], q_d_vol_ii, i_sum_type);
            cs_dispatch_sum(&v_slope_test[jj], q_d_vol_jj, i_sum_type);
          }
        }

        /* Save values at internal faces, if needed */
        if (i_pvar != nullptr) {
          if (i_massflux[face_id] >= 0.) {
            for (cs_lnum_t i = 0; i < stride; i++)
              i_pvar[face_id][i] += thetap * pif[i];
          }
          else {
            for (cs_lnum_t i = 0; i < stride; i++)
              i_pvar[face_id][i] += thetap * pjf[i];
          }
        }

        for (cs_lnum_t isou = 0; isou < stride; isou++)
          fluxi[isou] *= -1;  /* For substraction in dispatch sum */

        if (ii < n_cells)
          cs_dispatch_sum<stride>(rhs[ii], fluxi, i_sum_type);
        if (jj < n_cells)
          cs_dispatch_sum<stride>(rhs[jj], fluxj, i_sum_type);
      });
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcv"));

    }

  } /* pure upwind, without slope test, with slope test */

  ctx.wait();

  if (eqp.verbosity >= 2) {
    cs_gnum_t n_upwind = 0;

    if (i_upwind != nullptr) {
      ctx.parallel_for_reduce_sum
        (n_i_faces, n_upwind, [=] CS_F_HOST_DEVICE
         (cs_lnum_t i,
          CS_DISPATCH_REDUCER_TYPE(cs_gnum_t) &sum) {
        sum += (cs_gnum_t)i_upwind[i];
      });

      ctx.wait();
      cs_parall_counter(&n_upwind, 1);
    }
    else if (pure_upwind)
      n_upwind = m->n_g_i_faces;

    bft_printf(_(" %s: %llu Faces with upwind on %llu interior faces\n"),
               var_name, (unsigned long long)n_upwind,
               (unsigned long long)m->n_g_i_c_faces);

    CS_FREE_HD(i_upwind);
  }

  /* ======================================================================
     Contribution from boundary faces
     ======================================================================*/

  /* Boundary convective fluxes are all computed with an upwind scheme */
  if (icvflb == 0) {

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = b_face_cells[face_id];

      cs_real_t fluxi[stride], b_val_g[stride], b_val_d[stride];
      for (cs_lnum_t isou = 0; isou < stride; isou++) {
        fluxi[isou] = 0;
        b_val_g[isou] = val_f[face_id][isou];
        b_val_d[isou] = val_f_d_lim[face_id][isou];
      }
      cs_real_t _pi[stride];

      for (cs_lnum_t i = 0; i < stride; i++) {
        _pi[i]  = _pvar[ii][i];
      }

      /* Scaling due to mass balance in porous modeling */
      if (porous_vel) {
        const cs_nreal_t *n = b_face_u_normal[face_id];
        cs_math_3_normal_scaling(n, b_f_face_factor[face_id], _pi);
      }

      cs_b_upwind_flux_strided<stride>(iconvp,
                                       thetap,
                                       imasac,
                                       bc_type[face_id],
                                       _pi,
                                       _pi, /* no relaxation */
                                       b_massflux[face_id],
                                       b_val_g,
                                       fluxi);

      /* Save values on boundary faces */
      if (b_pvar != nullptr) {
        if (b_massflux[face_id] >= 0.) {
          for (cs_lnum_t i = 0; i < stride; i++)
            b_pvar[face_id][i] += thetap * _pi[i];
        }
        else {
          for (cs_lnum_t i = 0; i < stride; i++) {
            b_pvar[face_id][i] += thetap * b_val_g[i];
          }
        }
      }

      cs_b_diff_flux_strided<stride>(idiffp,
                                     thetap,
                                     b_visc[face_id],
                                     b_val_d,
                                     fluxi);

      for (cs_lnum_t isou = 0; isou < stride; isou++)
        fluxi[isou] *= -1;  /* For substraction in dispatch sum */
      cs_dispatch_sum<stride>(rhs[ii], fluxi, b_sum_type);

    });
  }

  /* Boundary convective flux imposed at some faces (tags in icvfli array) */
  else if (icvflb == 1) {

    /* Retrieve the value of the convective flux to be imposed */
    if (f_id != -1) {
      coface = (const var_t *)(f->bc_coeffs->ac);
      cofbce = (const b_t *)(f->bc_coeffs->bc);
    }
    else {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of icvflb and f_id"));
    }

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = b_face_cells[face_id];

      cs_real_t fluxi[stride], b_val_g[stride], b_val_d[stride];
      for (cs_lnum_t isou = 0; isou < stride; isou++) {
        fluxi[isou] = 0;
        b_val_g[isou] = val_f[face_id][isou];
        b_val_d[isou] = val_f_d_lim[face_id][isou];
      }

      cs_real_t pip[stride], _pi[stride];
      for (cs_lnum_t i = 0; i < stride; i++) {
        _pi[i]  = _pvar[ii][i];
      }

      /* Scaling due to mass balance in porous modelling */
      if (porous_vel) {
        const cs_nreal_t *n = b_face_u_normal[face_id];
        cs_math_3_normal_scaling(n, b_f_face_factor[face_id], _pi);
      }

      { /* This part is to compute val_iprime to compute
           face_value (local b_val_g) with coeffs AC and BC */

        cs_real_t bldfrp = (cs_real_t) ircflb;
        /* Local limitation of the reconstruction */
        if (df_limiter != nullptr && ircflb > 0)
          bldfrp = cs::max(df_limiter[ii], 0.);

        cs_b_cd_unsteady_strided<stride>(bldfrp,
                                         diipb[face_id],
                                         grad[ii],
                                         _pi,
                                         pip);
      }

      cs_b_imposed_conv_flux_strided<stride>(iconvp,
                                             thetap,
                                             imasac,
                                             inc,
                                             bc_type[face_id],
                                             icvfli[face_id],
                                             _pi,
                                             _pi, /* no relaxation */
                                             pip,
                                             coface[face_id],
                                             cofbce[face_id],
                                             b_massflux[face_id],
                                             b_val_g,
                                             fluxi);

      cs_b_diff_flux_strided<stride>(idiffp,
                                     thetap,
                                     b_visc[face_id],
                                     b_val_d,
                                     fluxi);

      /* Saving velocity on boundary faces if needed */
      if (b_pvar != nullptr) {
        if (b_massflux[face_id] >= 0.) {
          for (cs_lnum_t i = 0; i < stride; i++)
            b_pvar[face_id][i] += thetap * _pi[i];
        }
        else {
          for (cs_lnum_t i = 0; i < stride; i++)
            b_pvar[face_id][i] += thetap * b_val_g[i];
        }
      }

      for (cs_lnum_t isou = 0; isou < stride; isou++)
        fluxi[isou] *= -1;  /* For substraction in dispatch sum */
      cs_dispatch_sum<stride>(rhs[ii], fluxi, b_sum_type);

    });

  }

  ctx.wait();

  /* Free memory */
  CS_FREE_HD(grdpa);

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();

    std::chrono::microseconds elapsed;
    printf("%d: %s<%d>", cs_glob_rank_id, __func__, stride);

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_stop - t_start);
    printf(", total = %ld\n", elapsed.count());
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return pointer to slope test indicator field values if active.
 *
 * parameters:
 *   f_id  <-- field id (or -1)
 *   eqp   <-- equation parameters
 *
 * return:
 *   pointer to local values array, or nullptr;
 *----------------------------------------------------------------------------*/

cs_real_t *
cs_get_v_slope_test(int                        f_id,
                    const cs_equation_param_t  eqp)
{
  const int iconvp = eqp.iconv;
  const int isstpp = eqp.isstpc;
  const double blencp = eqp.blencv;

  cs_real_t  *v_slope_test = nullptr;

  if (f_id > -1 && iconvp > 0 && blencp > 0. && isstpp == 0) {

    static int _k_slope_test_f_id = -1;

    cs_field_t *f = cs_field_by_id(f_id);

    int f_track_slope_test_id = -1;

    if (_k_slope_test_f_id < 0)
      _k_slope_test_f_id = cs_field_key_id_try("slope_test_upwind_id");
    if (_k_slope_test_f_id > -1 && isstpp == 0)
      f_track_slope_test_id = cs_field_get_key_int(f, _k_slope_test_f_id);

    if (f_track_slope_test_id > -1)
      v_slope_test = (cs_field_by_id(f_track_slope_test_id))->val;

    if (v_slope_test != nullptr) {
      const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
#     pragma omp parallel for
      for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++)
        v_slope_test[cell_id] = 0.;
    }

  }

  return v_slope_test;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the beta blending coefficient of the
 * beta limiter (ensuring preservation of a given min/max pair of values).
 *
 * \param[in]     f_id         field id
 * \param[in]     inc          "not an increment" flag
 * \param[in]     rovsdt       rho * volume / dt
 */
/*----------------------------------------------------------------------------*/

void
cs_beta_limiter_building(int                   f_id,
                         int                   inc,
                         const cs_real_t       rovsdt[])
{
  /* Parallel or device dispatch */
  cs_dispatch_context ctx;

  /* Get options from the field */

  cs_field_t *f = cs_field_by_id(f_id);
  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);

  if (eqp->isstpc != 2)
    return;

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_halo_t *halo = m->halo;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  cs_real_t* cv_limiter = cs_field_by_id(
      cs_field_get_key_int(f, cs_field_key_id("convection_limiter_id")))->val;

  cs_real_t *denom_inf, *num_inf, *denom_sup, *num_sup;
  CS_MALLOC_HD(denom_inf, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(denom_sup, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(num_inf,   n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(num_sup,   n_cells_ext, cs_real_t, cs_alloc_mode);

  /* computation of the denominator for the inferior and upper bound */

  _beta_limiter_denom(f,
                      ctx,
                      eqp,
                      inc,
                      denom_inf,
                      denom_sup);

  /* computation of the numerator for the inferior and upper bound */

  _beta_limiter_num(f,
                    ctx,
                    eqp,
                    inc,
                    rovsdt,
                    num_inf,
                    num_sup);

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    /* Treatment of the lower bound */
    cs_real_t beta_inf;
    if (denom_inf[c_id] <= num_inf[c_id]) {
      beta_inf = 1.;
    }
    else if (denom_inf[c_id] <= cs::abs(num_inf[c_id])) {
      beta_inf = -1.;
    }
    else {
      beta_inf = num_inf[c_id]/denom_inf[c_id]; //FIXME division by 0
      beta_inf = cs::min(beta_inf, 1.);
    }

    /* Treatment of the upper bound */
    cs_real_t beta_sup;
    if (denom_sup[c_id] <= num_sup[c_id]) {
      beta_sup = 1.;
    }
    else if (denom_sup[c_id] <= cs::abs(num_sup[c_id])) {
      beta_sup = -1.;
    }
    else {
      beta_sup = num_sup[c_id]/denom_sup[c_id]; //FIXME division by 0
      beta_sup = cs::min(beta_sup, 1.);
    }

    cv_limiter[c_id] = cs::min(beta_inf, beta_sup);
  });

  ctx.wait();

  /* Synchronize variable */
  if (halo != nullptr)
    cs_halo_sync(halo, CS_HALO_STANDARD, ctx.use_gpu(), cv_limiter);

  CS_FREE_HD(denom_inf);
  CS_FREE_HD(denom_sup);
  CS_FREE_HD(num_inf);
  CS_FREE_HD(num_sup);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the convection/diffusion terms of a
 * standard transport equation of a scalar field \f$ \varia \f$.
 *
 * More precisely, the right hand side \f$ Rhs \f$ is updated as
 * follows:
 * \f[
 * Rhs = Rhs - \sum_{\fij \in \Facei{\celli}}      \left(
 *        \dot{m}_\ij \left( \varia_\fij - \varia_\celli \right)
 *      - \mu_\fij \gradv_\fij \varia \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ Rhs \f$ has already been initialized before calling bilsc2!
 * - mind the sign minus
 *
 * Please refer to the
 * <a href="../../theory.pdf#bilsc2"><b>bilsc2</b></a> section of the
 * theory guide for more informations.
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          field id (or -1)
 * \param[in]     eqp           equation parameters
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imasac        take mass accumulation into account?
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in]     bc_coeffs     boundary condition structure for the variable
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 * \param[in,out] i_flux        interior flux (or nullptr)
 * \param[in,out] b_flux        boundary flux (or nullptr)
 */
/*----------------------------------------------------------------------------*/

void
cs_convection_diffusion_scalar(int                         idtvar,
                               int                         f_id,
                               const cs_equation_param_t   eqp,
                               int                         icvflb,
                               int                         inc,
                               int                         imasac,
                               cs_real_t         *restrict pvar,
                               const cs_real_t   *restrict pvara,
                               const int                   icvfli[],
                               const cs_field_bc_coeffs_t *bc_coeffs,
                               const cs_real_t             i_massflux[],
                               const cs_real_t             b_massflux[],
                               const cs_real_t             i_visc[],
                               const cs_real_t             b_visc[],
                               cs_real_t                  *rhs,
                               cs_real_2_t                 i_flux[],
                               cs_real_t                   b_flux[])
{
  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  cs_field_t *f = nullptr;
  if (f_id != -1)
    f = cs_field_by_id(f_id);

  if (idtvar >= 0) {
    if (i_flux == nullptr)
      _convection_diffusion_scalar_unsteady<false, false>
        (f, eqp, icvflb, inc, imasac,
         pvar, pvara,
         icvfli,
         bc_coeffs,
         i_massflux, b_massflux,
         i_visc, b_visc,
         nullptr, rhs, i_flux, b_flux);
    else
      _convection_diffusion_scalar_unsteady<false, true>
        (f, eqp, icvflb, inc, imasac,
         pvar, pvara,
         icvfli,
         bc_coeffs,
         i_massflux, b_massflux,
         i_visc, b_visc,
         nullptr, rhs, i_flux, b_flux);
  }
  else {

  /* The following function, for the steady (idtvar < 0 case) is deprecated
     ---------------------------------------------------------------------- */

    _convection_diffusion_scalar_steady
      (f, eqp, icvflb, inc,
       pvar, pvara,
       icvfli,
       bc_coeffs,
       i_massflux, b_massflux,
       i_visc, b_visc,
       nullptr, rhs);
  }

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();

    std::chrono::microseconds elapsed;
    printf("%d: %s", cs_glob_rank_id, __func__);

    elapsed = std::chrono::duration_cast
      <std::chrono::microseconds>(t_stop - t_start);
    printf(", total = %ld\n", elapsed.count());
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update face flux with convection contribution of a standard transport
 * equation of a scalar field \f$ \varia \f$.
 *
 * <a name="cs_face_convection_scalar"></a>
 *
 * \f[
 * C_\ij = \dot{m}_\ij \left( \varia_\fij - \varia_\celli \right)
 * \f]
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          pointer to field
 * \param[in]     eqp           equation parameters
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imasac        take mass accumulation into account?
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in]     bc_coeffs     boundary condition structure for the variable
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in,out] i_conv_flux   scalar convection flux at interior faces
 * \param[in,out] b_conv_flux   scalar convection flux at boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_face_convection_scalar(int                         idtvar,
                          int                         f_id,
                          const cs_equation_param_t   eqp,
                          int                         icvflb,
                          int                         inc,
                          int                         imasac,
                          cs_real_t         *restrict pvar,
                          const cs_real_t   *restrict pvara,
                          const int                   icvfli[],
                          const cs_field_bc_coeffs_t *bc_coeffs,
                          const cs_real_t             i_massflux[],
                          const cs_real_t             b_massflux[],
                          cs_real_t                   i_conv_flux[][2],
                          cs_real_t                   b_conv_flux[])
{
  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  cs_field_t *f = nullptr;
  if (f_id != -1)
    f = cs_field_by_id(f_id);

  if (idtvar >= 0) {
    _face_convection_scalar_unsteady
      (f, eqp, icvflb, inc, imasac,
       pvar, pvara,
       icvfli,
       bc_coeffs,
       i_massflux, b_massflux,
       i_conv_flux, b_conv_flux);
  }
  else {

    /* The following function, for the steady (idtvar < 0 case) is deprecated
       ---------------------------------------------------------------------- */

    _face_convection_scalar_steady
      (f, eqp, icvflb, inc,
       pvar, pvara,
       icvfli,
       bc_coeffs,
       i_massflux, b_massflux,
       i_conv_flux, b_conv_flux);
  }

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();

    std::chrono::microseconds elapsed;
    printf("%d: %s", cs_glob_rank_id, __func__);

    elapsed = std::chrono::duration_cast
      <std::chrono::microseconds>(t_stop - t_start);
    printf(", total = %ld\n", elapsed.count());
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the convection/diffusion terms of a transport
 *  equation of a vector field \f$ \vect{\varia} \f$.
 *
 * More precisely, the right hand side \f$ \vect{Rhs} \f$ is updated as
 * follows:
 * \f[
 *  \vect{Rhs} = \vect{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
 *         \dot{m}_\ij \left( \vect{\varia}_\fij - \vect{\varia}_\celli \right)
 *       - \mu_\fij \gradt_\fij \vect{\varia} \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Remark:
 * if ivisep = 1, then we also take \f$ \mu \transpose{\gradt\vect{\varia}}
 * + \lambda \trace{\gradt\vect{\varia}} \f$, where \f$ \lambda \f$ is
 * the secondary viscosity, i.e. usually \f$ -\frac{2}{3} \mu \f$.
 *
 * Warning:
 * - \f$ \vect{Rhs} \f$ has already been initialized before calling bilsc!
 * - mind the sign minus
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          index of the current variable
 * \param[in]     eqp           equation parameters
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     ivisep        indicator to take \f$ \divv
 *                               \left(\mu \gradt \transpose{\vect{a}} \right)
 *                               -2/3 \grad\left( \mu \dive \vect{a} \right)\f$
 *                               - 1 take into account,
 *                               - 0 otherwise
 * \param[in]     imasac        take mass accumulation into account?
 * \param[in]     pvar          solved velocity (current time step)
 * \param[in]     pvara         solved velocity (previous time step)
 * \param[in]     icvfli        boundary face indicator array of convection flux
 *                               - 0 upwind scheme
 *                               - 1 imposed flux
 * \param[in]     bc_coeffs_v   boundary conditions structure for the variable
 * \param[in]     bc_coeffs_solve_v   sweep loop boundary conditions structure
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     i_secvis      secondary viscosity at interior faces
 * \param[in]     b_secvis      secondary viscosity at boundary faces
 * \param[in]     i_pvar        velocity at interior faces
 * \param[in]     b_pvar        velocity at boundary faces
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_convection_diffusion_vector(int                         idtvar,
                               int                         f_id,
                               const cs_equation_param_t   eqp,
                               int                         icvflb,
                               int                         inc,
                               int                         ivisep,
                               int                         imasac,
                               cs_real_3_t       *restrict pvar,
                               const cs_real_3_t *restrict pvara,
                               const int                   icvfli[],
                               const cs_field_bc_coeffs_t *bc_coeffs_v,
                               const cs_bc_coeffs_solve_t *bc_coeffs_solve_v,
                               const cs_real_t             i_massflux[],
                               const cs_real_t             b_massflux[],
                               const cs_real_t             i_visc[],
                               const cs_real_t             b_visc[],
                               const cs_real_t             i_secvis[],
                               const cs_real_t             b_secvis[],
                               cs_real_3_t       *restrict i_pvar,
                               cs_real_3_t       *restrict b_pvar,
                               cs_real_3_t       *restrict rhs)
{
  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  /* Initialization
     -------------- */

  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  static const char _novar_name[] = "[convection-diffusion, vector]";

  cs_field_t *f = nullptr;
  const char *var_name = (const char *)_novar_name;
  if (f_id != -1) {
    f = cs_field_by_id(f_id);
    var_name = (const char *)(f->name);
  }

  /* Parallel or device dispatch */

  cs_dispatch_context ctx;
  if (idtvar < 0)
    ctx.set_use_gpu(false);  /* steady case not ported to GPU */

  cs_alloc_mode_t amode = ctx.alloc_mode(true);

  /* Halo (and gradient) type */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;
  cs_gradient_type_by_imrgra(eqp.d_gradient_r,
                             &gradient_type,
                             &halo_type);

  /* Handle cases where only the previous values (already synchronized)
     or current values are provided */

  if (pvar != nullptr && m->halo != nullptr) {
    cs_halo_sync_r(m->halo, halo_type, ctx.use_gpu(), pvar);
  }
  if (pvara == nullptr)
    pvara = (const cs_real_3_t *)pvar;

  /* Marker for inlet/outlet cells for transpose gradient
     ----------------------------------------------------

     We start the computation now, so that on a GPU, this can be done
     in a separate stream, concurrently with the heavier computation
     below. This is done so as to compensate for the kernel launch
     latency of this very light operation.
  */

  short *bndcel = nullptr;
  cs_dispatch_context ctx_c;
  ctx_c.set_use_gpu(ctx.use_gpu()); /* Follows behavior of main context */
#if defined(HAVE_CUDA)
  if (ctx_c.use_gpu())
    ctx_c.set_cuda_stream(cs_cuda_get_stream(1));
#endif

  if (ivisep == 1 && eqp.idiff == 1) {
    const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
    const int *bc_type = cs_glob_bc_type;

    /* We do not know what condition to put in the inlets and the outlets, so we
       assume that there is an equilibrium. Moreover, cells containing a
       coupled face (velocity-wise) are ignored also. */

    /* Allocate a temporary array */
    CS_MALLOC_HD(bndcel, n_cells_ext, short, amode);

    ctx_c.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      bndcel[cell_id] = 1;
    });

    ctx_c.parallel_for(m->n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      int ityp = bc_type[face_id];
      if (   ityp == CS_OUTLET
          || ityp == CS_INLET
          || ityp == CS_FREE_INLET
          || ityp == CS_CONVECTIVE_INLET
          || ityp == CS_COUPLED_FD)
        bndcel[b_face_cells[face_id]] = 0;
    });
  }

  /* Compute the gradient of the vector (velocity)
     ----------------------------------

     The gradient (grad) is used in the flux reconstruction and the slope test.
     Thus we must compute it:
         - when we have diffusion and we reconstruct the fluxes,
         - when the convection scheme is SOLU,
         - when we have convection, we are not in pure upwind
           and we reconstruct the fluxes,
         - when we have convection, we are not in pure upwind
           and we have not shunted the slope test. */

  using grad_t = cs_real_t[3][3];
  grad_t *grad;
  CS_MALLOC_HD(grad, n_cells_ext, grad_t, amode);

  if (  (eqp.idiff != 0 && eqp.ircflu == 1)
      || ivisep == 1
      || (    eqp.iconv != 0 && eqp.blencv > 0.
          && (eqp.ischcv == 0 || eqp.ircflu == 1 || eqp.isstpc == 0))) {

    cs_real_t *gweight = nullptr;
    if (f_id != -1) {
      /* Get the calculation option from the field */
      if (f->type & CS_FIELD_VARIABLE && eqp.iwgrec == 1) {
        if (eqp.idiff > 0) {
          int key_id = cs_field_key_id("gradient_weighting_id");
          int diff_id = cs_field_get_key_int(f, key_id);
          if (diff_id > -1) {
            cs_field_t *weight_f = cs_field_by_id(diff_id);
            gweight = weight_f->val;
            cs_field_synchronize(weight_f, halo_type);
          }
        }
      }
    }

    const cs_real_3_t  *restrict _pvar
      = (pvar != nullptr) ? (const cs_real_3_t  *)pvar : pvara;

    const cs_real_3_t *val_f = (bc_coeffs_solve_v == nullptr) ?
                               (const cs_real_3_t *)bc_coeffs_v->val_f :
                               (const cs_real_3_t *)bc_coeffs_solve_v->val_f;

    cs_gradient_vector_synced_input(var_name,
                                    gradient_type,
                                    halo_type,
                                    inc,
                                    eqp.nswrgr,
                                    eqp.verbosity,
                                    (cs_gradient_limit_t)(eqp.imligr),
                                    eqp.epsrgr,
                                    eqp.d_climgr,
                                    bc_coeffs_v,
                                    (const cs_real_3_t *)_pvar,
                                    val_f,
                                    gweight, /* weighted gradient */
                                    (cs_real_33_t *)grad);

  }
  else {
    ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      for (cs_lnum_t isou = 0; isou < 3; isou++) {
        for (cs_lnum_t jsou = 0; jsou < 3; jsou++)
          grad[cell_id][isou][jsou] = 0.;
      }
    });
    ctx.wait();
  }

  /* Main convection/diffusion block
     ------------------------------- */

  if (idtvar >= 0) {
    bool porous_vel = false;
    if (cs_glob_porous_model == 3 && f == CS_F_(vel))
      porous_vel = true;
    if (porous_vel == false)
      _convection_diffusion_unsteady_strided<3, false>
        (ctx, f, var_name, eqp, icvflb, inc, imasac,
         pvar, pvara,
         icvfli,
         bc_coeffs_v,
         bc_coeffs_solve_v,
         i_massflux, b_massflux,
         i_visc, b_visc,
         i_pvar, b_pvar, grad, rhs);
    else
      _convection_diffusion_unsteady_strided<3, true>
        (ctx, f, var_name, eqp, icvflb, inc, imasac,
         pvar, pvara,
         icvfli,
         bc_coeffs_v,
         bc_coeffs_solve_v,
         i_massflux, b_massflux,
         i_visc, b_visc,
         i_pvar, b_pvar, grad, rhs);
  }

  else {
    _convection_diffusion_vector_steady
      (f, var_name, eqp, icvflb, inc,
       pvar, pvara,
       icvfli,
       bc_coeffs_v,
       bc_coeffs_solve_v,
       i_massflux, b_massflux,
       i_visc, b_visc,
       grad, rhs);
  }

  /* Computation of the transpose grad(vel) term and grad(-2/3 div(vel))
     ------------------------------------------------------------------- */

  if (ivisep == 1 && eqp.idiff == 1) {

    cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

    const cs_lnum_t n_cells = m->n_cells;

    const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
    const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

    const cs_real_t *restrict weight = fvq->weight;
    const cs_real_t *restrict i_dist = fvq->i_dist;
    const cs_nreal_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;
    const cs_real_t *restrict i_f_face_surf = fvq->i_face_surf;
    const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
    const cs_real_t *restrict b_f_face_surf = fvq->b_face_surf;

    const double thetap = eqp.theta;

    cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
    cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

    ctx_c.wait();  /* We now need bndcel, computed by ctx_c */

    /* Interior faces */

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t w_f = weight[face_id];
      cs_real_t secvis = i_secvis[face_id]; /* - 2/3 * mu */
      cs_real_t visco = i_visc[face_id]; /* mu S_ij / d_ij */

      cs_real_t grdtrv =      w_f  * cs_math_33_trace(grad[ii])
                        + (1.-w_f) * cs_math_33_trace(grad[jj]);

      cs_real_t ipjp[3];
      const cs_nreal_t *n = i_face_u_normal[face_id];

      /* I'J' */
      for (cs_lnum_t i = 0; i < 3; i++)
        ipjp[i] = n[i] * i_dist[face_id];

      cs_real_t _i_f_face_surf = i_f_face_surf[face_id];
      /* We need to compute trans_grad(u).I'J' which is equal to I'J'.grad(u) */

      cs_real_t  flux_p[3], flux_m[3];
      cs_real_t  theta_bc_mult_i = thetap * bndcel[ii];
      cs_real_t  theta_bc_mult_j = thetap * bndcel[jj];

      for (cs_lnum_t isou = 0; isou < 3; isou++) {

        cs_real_t tgrdfl =   ipjp[0] * (      w_f*grad[ii][0][isou]
                                        + (1.-w_f)*grad[jj][0][isou])
                           + ipjp[1] * (      w_f*grad[ii][1][isou]
                                        + (1.-w_f)*grad[jj][1][isou])
                           + ipjp[2] * (      w_f*grad[ii][2][isou]
                                        + (1.-w_f)*grad[jj][2][isou]);

        cs_real_t flux =   visco*tgrdfl
                         + secvis*grdtrv*i_face_u_normal[face_id][isou]
                                        *_i_f_face_surf;

        flux_p[isou] =  theta_bc_mult_i * flux;
        flux_m[isou] = -theta_bc_mult_j * flux;

      }

      if (ii < n_cells)
        cs_dispatch_sum<3>(rhs[ii], flux_p, i_sum_type);
      if (jj < n_cells)
        cs_dispatch_sum<3>(rhs[jj], flux_m, i_sum_type);

    });

    /* Boundary faces
       the whole flux term of the stress tensor is already taken into account
       (so, no corresponding term in forbr)
       TODO in theory we should take the normal component into account (the
       tangential one is modeled by the wall law) */

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = b_face_cells[face_id];

      /* Trace of velocity gradient */
      cs_real_t grdtrv = (grad[ii][0][0]+grad[ii][1][1]+grad[ii][2][2]);
      cs_real_t secvis = b_secvis[face_id]; /* - 2/3 * mu */
      cs_real_t mult = thetap * bndcel[ii]
        * secvis * grdtrv * b_f_face_surf[face_id];

      cs_real_t flux[3];
      for (cs_lnum_t isou = 0; isou < 3; isou++) {
        flux[isou] = mult * b_face_u_normal[face_id][isou];
      }
      cs_dispatch_sum<3>(rhs[ii], flux, b_sum_type);

    });

  }

  ctx.wait();

  /* Free memory */
  CS_FREE_HD(bndcel);
  CS_FREE_HD(grad);

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();

    std::chrono::microseconds elapsed;
    printf("%d: %s", cs_glob_rank_id, __func__);

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_stop - t_start);
    printf(", total = %ld\n", elapsed.count());
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the convection/diffusion terms of a transport
 *  equation of a tensor field \f$ \tens{\varia} \f$.
 *
 * More precisely, the right hand side \f$ \tens{Rhs} \f$ is updated as
 * follows:
 * \f[
 *  \tens{Rhs} = \tens{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
 *         \dot{m}_\ij \left( \tens{\varia}_\fij - \tens{\varia}_\celli \right)
 *       - \mu_\fij \gradt_\fij \tens{\varia} \cdot \tens{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ \tens{Rhs} \f$ has already been initialized before calling bilsc!
 * - mind the sign minus
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          index of the current variable
 * \param[in]     eqp           equation parameters
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imasac        take mass accumulation into account?
 * \param[in]     pvar          solved velocity (current time step)
 * \param[in]     pvara         solved velocity (previous time step)
 * \param[in]     bc_coeffs_ts  boundary condition structure for the variable
 * \param[in]     bc_coeffs_solve_ts   sweep loop boundary conditions structure
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in,out] rhs           right hand side \f$ \tens{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_convection_diffusion_tensor(int                          idtvar,
                               int                          f_id,
                               const cs_equation_param_t    eqp,
                               int                          icvflb,
                               int                          inc,
                               int                          imasac,
                               cs_real_6_t        *restrict pvar,
                               const cs_real_6_t  *restrict pvara,
                               const cs_field_bc_coeffs_t  *bc_coeffs_ts,
                               const cs_bc_coeffs_solve_t  *bc_coeffs_solve_ts,
                               const cs_real_t              i_massflux[],
                               const cs_real_t              b_massflux[],
                               const cs_real_t              i_visc[],
                               const cs_real_t              b_visc[],
                               cs_real_6_t        *restrict rhs)
{
  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  /* Initialization
     -------------- */

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  static const char _novar_name[] = "[convection-diffusion, tensor]";

  cs_field_t *f = nullptr;
  const char *var_name = (const char *)_novar_name;
  if (f_id != -1) {
    f = cs_field_by_id(f_id);
    var_name = (const char *)(f->name);
  }

  /* Parallel or device dispatch */

  cs_dispatch_context ctx;
  if (idtvar < 0)
    ctx.set_use_gpu(false);  /* steady case not ported to GPU */

  cs_alloc_mode_t amode = ctx.alloc_mode(true);

  /* Halo (and gradient) type */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;
  cs_gradient_type_by_imrgra(eqp.d_gradient_r,
                             &gradient_type,
                             &halo_type);

  /* Handle cases where only the previous values (already synchronized)
     or current values are provided */

  if (pvar != nullptr && m->halo != nullptr) {
    cs_halo_sync_r(m->halo, halo_type, ctx.use_gpu(), pvar);
  }
  if (pvara == nullptr)
    pvara = (const cs_real_6_t *)pvar;

  /* Compute the gradient of the tensor

     The gradient (grad) is used in the flux reconstruction and the slope test.
     Thus we must compute it:
         - when we have diffusion and we reconstruct the fluxes,
         - when the convection scheme is SOLU,
         - when we have convection, we are not in pure upwind
           and we reconstruct the fluxes,
         - when we have convection, we are not in pure upwind
           and we have not shunted the slope test. */

  using grad_t = cs_real_t[6][3];
  grad_t *grad;
  CS_MALLOC_HD(grad, n_cells_ext, grad_t, amode);

  if (  (eqp.idiff != 0 && eqp.ircflu == 1)
      || (   eqp.iconv != 0 && eqp.blencv > 0.
          && (eqp.ischcv == 0 || eqp.ircflu == 1 || eqp.isstpc == 0))) {

    const cs_gradient_limit_t imligp = (cs_gradient_limit_t)(eqp.imligr);
    const cs_real_6_t  *restrict _pvar
      = (pvar != nullptr) ? (const cs_real_6_t  *)pvar : pvara;

    const cs_real_6_t *val_f = (bc_coeffs_solve_ts == nullptr) ?
                               (const cs_real_6_t *)bc_coeffs_ts->val_f :
                               (const cs_real_6_t *)bc_coeffs_solve_ts->val_f;

    cs_gradient_tensor_synced_input(var_name,
                                    gradient_type,
                                    halo_type,
                                    inc,
                                    eqp.nswrgr,
                                    eqp.verbosity,
                                    imligp,
                                    eqp.epsrgr,
                                    eqp.d_climgr,
                                    bc_coeffs_ts,
                                    _pvar,
                                    val_f,
                                    (cs_real_63_t *)grad);

  }
  else {
    ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      for (cs_lnum_t isou = 0; isou < 6; isou++) {
        for (cs_lnum_t jsou = 0; jsou < 3; jsou++)
          grad[cell_id][isou][jsou] = 0.;
      }
    });
    ctx.wait();
  }

  /* Main convection/diffusion block
     ------------------------------- */

  if (idtvar >= 0) {
    assert(icvflb == 0);
    _convection_diffusion_unsteady_strided<6, false>
      (ctx, f, var_name, eqp, 0, inc, imasac,
       pvar, pvara,
       nullptr, // icvfli,
       bc_coeffs_ts,
       bc_coeffs_solve_ts,
       i_massflux, b_massflux,
       i_visc, b_visc,
       nullptr, nullptr, grad, rhs);
  }
  else {
    _convection_diffusion_tensor_steady
      (f, var_name, eqp, icvflb, inc,
       pvar, pvara,
       bc_coeffs_ts,
       bc_coeffs_solve_ts,
       i_massflux, b_massflux,
       i_visc, b_visc,
       grad, rhs);
  }

  /* Free memory */
  CS_FREE_HD(grad);

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();

    std::chrono::microseconds elapsed;
    printf("%d: %s", cs_glob_rank_id, __func__);

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_stop - t_start);
    printf(", total = %ld\n", elapsed.count());
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the convection/diffusion terms of a transport
 * equation of a scalar field \f$ \varia \f$ such as the temperature.
 *
 * More precisely, the right hand side \f$ Rhs \f$ is updated as
 * follows:
 * \f[
 * Rhs = Rhs + \sum_{\fij \in \Facei{\celli}}      \left(
 *        C_p\dot{m}_\ij \varia_\fij
 *      - \lambda_\fij \gradv_\fij \varia \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * \f$ Rhs \f$ has already been initialized before calling bilsct!
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          index of the current variable
 * \param[in]     eqp           equation parameters)
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imasac        take mass accumulation into account?
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     bc_coeffs     boundary condition structure for the variable
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     xcpp          array of specific heat (\f$ C_p \f$)
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_convection_diffusion_thermal(int                         idtvar,
                                int                         f_id,
                                const cs_equation_param_t   eqp,
                                int                         inc,
                                int                         imasac,
                                cs_real_t        *restrict  pvar,
                                const cs_real_t  *restrict  pvara,
                                const cs_field_bc_coeffs_t *bc_coeffs,
                                const cs_real_t             i_massflux[],
                                const cs_real_t             b_massflux[],
                                const cs_real_t             i_visc[],
                                const cs_real_t             b_visc[],
                                const cs_real_t             xcpp[],
                                cs_real_t        *restrict  rhs)
{
  cs_field_t *f = nullptr;
  if (f_id != -1)
    f = cs_field_by_id(f_id);

  if (idtvar >= 0) {
    _convection_diffusion_scalar_unsteady<true, false>
      (f, eqp,
       false, /* icvflb */
       inc, imasac,
       pvar, pvara,
       nullptr, /* icvfli */
       bc_coeffs,
       i_massflux, b_massflux,
       i_visc, b_visc,
       xcpp, rhs, nullptr, nullptr);
    return;
  }

  /* The following function, for the steady (idtvar < 0 case) is deprecated
     ---------------------------------------------------------------------- */

  _convection_diffusion_scalar_steady
    (f, eqp,
     false, /* icvflb */
     inc,
     pvar, pvara,
     nullptr, /* icvfli */
     bc_coeffs,
     i_massflux, b_massflux,
     i_visc, b_visc,
     xcpp, rhs);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the diffusion terms with a symmetric tensor
 * diffusivity for a transport equation of a scalar field \f$ \varia \f$.
 *
 * More precisely, the right hand side \f$ Rhs \f$ is updated as
 * follows:
 * \f[
 * Rhs = Rhs - \sum_{\fij \in \Facei{\celli}}      \left(
 *      - \tens{\mu}_\fij \gradv_\fij \varia \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ Rhs \f$ has already been initialized before
 *   calling cs_anisotropic_diffusion_scalar!
 * - mind the sign minus
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          index of the current variable
 * \param[in]     eqp           equation parameters
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     bc_coeffs     boundary condition structure for the variable
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]     weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_anisotropic_diffusion_scalar(int                         idtvar,
                                int                         f_id,
                                const cs_equation_param_t   eqp,
                                int                         inc,
                                cs_real_t        *restrict  pvar,
                                const cs_real_t  *restrict  pvara,
                                const cs_field_bc_coeffs_t *bc_coeffs,
                                const cs_real_t             i_visc[],
                                const cs_real_t             b_visc[],
                                cs_real_6_t      *restrict  viscel,
                                const cs_real_2_t           weighf[],
                                const cs_real_t             weighb[],
                                cs_real_t        *restrict  rhs)
{
  const cs_real_t *cofafp = bc_coeffs->af;
  const cs_real_t *cofbfp = bc_coeffs->bf;

  const int nswrgp = eqp.nswrgr;
  const cs_gradient_limit_t imligp = (cs_gradient_limit_t)(eqp.imligr);
  const int ircflp = eqp.ircflu;
  const int ircflb = (ircflp > 0) ? eqp.b_diff_flux_rc : 0;
  const int iwarnp = eqp.verbosity;
  const int icoupl = eqp.icoupl;
  const double epsrgp = eqp.epsrgr;
  const double relaxp = eqp.relaxv;
  const double thetap = eqp.theta;

  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_halo_t  *halo = m->halo;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_nreal_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;
  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
  const cs_real_t *restrict i_face_surf = fvq->i_face_surf;
  const cs_real_t *restrict b_face_surf = fvq->b_face_surf;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog = fvq->b_face_cog;

  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  /* Local variables */

  cs_real_t *df_limiter = nullptr;

  char var_name[64];

  int w_stride = 1;

  cs_real_6_t *viscce = nullptr;
  cs_real_6_t *w2 = nullptr;
  cs_real_3_t *grad;

  cs_field_t *f = nullptr;

  cs_real_t *gweight = nullptr;

  /* Internal coupling variables */
  cs_real_t *pvar_local = nullptr;
  cs_real_3_t *grad_local = nullptr;
  cs_real_t *df_limiter_local = nullptr;
  cs_real_6_t *viscce_local = nullptr;
  cs_real_t *weighb_local = nullptr;
  const cs_lnum_t *faces_local = nullptr;
  cs_lnum_t n_local = 0;
  int coupling_id = -1;
  cs_internal_coupling_t *cpl = nullptr;

  /* 1. Initialization */

  /* Allocate work arrays */
  CS_MALLOC_HD(grad, n_cells_ext, cs_real_3_t, cs_alloc_mode);

  /* Choose gradient type */
  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  cs_gradient_type_by_imrgra(eqp.d_gradient_r,
                             &gradient_type,
                             &halo_type);

  /* Handle cases where only the previous values (already synchronized)
     or current values are provided */

  if (pvar != nullptr)
    cs_sync_scalar_halo(m, halo_type, pvar);
  if (pvara == nullptr)
    pvara = (const cs_real_t *)pvar;

  const cs_real_t  *restrict _pvar = (pvar != nullptr) ? pvar : pvara;

  /* Logging */

  if (f_id != -1) {
    f = cs_field_by_id(f_id);

    int df_limiter_id =
      cs_field_get_key_int(f, cs_field_key_id("diffusion_limiter_id"));
    if (df_limiter_id > -1)
      df_limiter = cs_field_by_id(df_limiter_id)->val;

    snprintf(var_name, 63, "%s", f->name);
  }
  else
    strncpy(var_name, "[anisotropic diffusion, scalar]", 63);
  var_name[63] = '\0';

  /* Porosity fields */
  cs_field_t *fporo = cs_field_by_name_try("porosity");
  cs_field_t *ftporo = cs_field_by_name_try("tensorial_porosity");

  cs_real_t *porosi = nullptr;
  cs_real_6_t *porosf = nullptr;

  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2) {
    porosi = fporo->val;
    if (ftporo != nullptr) {
      porosf = (cs_real_6_t *)ftporo->val;
    }
  }

  /* Without porosity */
  if (porosi == nullptr) {
    viscce = viscel;

    /* With porosity */
  }
  else if (porosi != nullptr && porosf == nullptr) {
    CS_MALLOC_HD(w2, n_cells_ext, cs_real_6_t, cs_alloc_mode);
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      for (cs_lnum_t isou = 0; isou < 6; isou++) {
        w2[cell_id][isou] = porosi[cell_id]*viscel[cell_id][isou];
      }
    });
    ctx.wait();
    viscce = w2;

    /* With tensorial porosity */
  }
  else if (porosi != nullptr && porosf != nullptr) {
    CS_MALLOC_HD(w2, n_cells_ext, cs_real_6_t, cs_alloc_mode);
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      cs_math_sym_33_product(porosf[cell_id],
                             viscel[cell_id],
                             w2[cell_id]);
    });
    ctx.wait();
    viscce = w2;
  }

  /* ---> Periodicity and parallelism treatment of symmetric tensors */
  cs_halo_sync_r(halo, halo_type, false, viscce);

  if (icoupl > 0) {
    assert(f_id != -1);
    const int coupling_key_id = cs_field_key_id("coupling_entity");
    coupling_id = cs_field_get_key_int(f, coupling_key_id);
    cpl = cs_internal_coupling_by_id(coupling_id);
    cs_internal_coupling_coupled_faces(cpl,
                                       &n_local,
                                       &faces_local,
                                       nullptr,
                                       nullptr);
  }

  /* 2. Compute the diffusive part with reconstruction technics */

  /* ======================================================================
     Compute the gradient of the current variable if needed
     ======================================================================*/

  if (ircflp == 1) {

    if (f_id != -1) {
      /* Get the calculation option from the field */
      if (f->type & CS_FIELD_VARIABLE && eqp.iwgrec == 1) {
        if (eqp.idifft > 0) {
          int key_id = cs_field_key_id("gradient_weighting_id");
          int diff_id = cs_field_get_key_int(f, key_id);
          if (diff_id > -1) {
            cs_field_t *weight_f = cs_field_by_id(diff_id);
            gweight = weight_f->val;
            w_stride = weight_f->dim;
            cs_field_synchronize(weight_f, halo_type);
          }
        }
      }
    }

    cs_gradient_scalar_synced_input(var_name,
                                    gradient_type,
                                    halo_type,
                                    inc,
                                    nswrgp,
                                    0, /* hyd_p_flag */
                                    w_stride,
                                    iwarnp,
                                    imligp,
                                    epsrgp,
                                    eqp.d_climgr,
                                    nullptr, /* f_ext exterior force */
                                    bc_coeffs,
                                    _pvar,
                                    gweight, /* Weighted gradient */
                                    cpl, /* internal coupling */
                                    grad);

  }
  else {
    ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      grad[cell_id][0] = 0.;
      grad[cell_id][1] = 0.;
      grad[cell_id][2] = 0.;
    });
  }

  /* ======================================================================
     ---> Contribution from interior faces
     ======================================================================*/

  /* Steady */
  if (idtvar < 0) {

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t pi = _pvar[ii];
      cs_real_t pj = _pvar[jj];
      cs_real_t pia = pvara[ii];
      cs_real_t pja = pvara[jj];

      cs_real_t bldfrp = (cs_real_t) ircflp;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflp > 0)
        bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                         0.);

      /* Recompute II" and JJ"
         ----------------------*/

      cs_real_t visci[3][3], viscj[3][3];

      visci[0][0] = viscce[ii][0];
      visci[1][1] = viscce[ii][1];
      visci[2][2] = viscce[ii][2];
      visci[1][0] = viscce[ii][3];
      visci[0][1] = viscce[ii][3];
      visci[2][1] = viscce[ii][4];
      visci[1][2] = viscce[ii][4];
      visci[2][0] = viscce[ii][5];
      visci[0][2] = viscce[ii][5];

      /* IF.Ki.S / ||Ki.S||^2 */
      cs_real_t fikdvi = weighf[face_id][0];

      cs_real_t diippf[3], djjppf[3];

      /* II" = IF + FI" */
      for (cs_lnum_t i = 0; i < 3; i++) {
        diippf[i] = i_face_cog[face_id][i]-cell_cen[ii][i]
                    - fikdvi*(  visci[0][i]*i_face_u_normal[face_id][0]
                              + visci[1][i]*i_face_u_normal[face_id][1]
                              + visci[2][i]*i_face_u_normal[face_id][2])
                            *i_face_surf[face_id];
      }

      viscj[0][0] = viscce[jj][0];
      viscj[1][1] = viscce[jj][1];
      viscj[2][2] = viscce[jj][2];
      viscj[1][0] = viscce[jj][3];
      viscj[0][1] = viscce[jj][3];
      viscj[2][1] = viscce[jj][4];
      viscj[1][2] = viscce[jj][4];
      viscj[2][0] = viscce[jj][5];
      viscj[0][2] = viscce[jj][5];

      /* FJ.Kj.S / ||Kj.S||^2 */
      cs_real_t fjkdvi = weighf[face_id][1];

      /* JJ" = JF + FJ" */
      for (cs_lnum_t i = 0; i < 3; i++) {
        djjppf[i] =   i_face_cog[face_id][i]-cell_cen[jj][i]
                    + fjkdvi*(  viscj[0][i]*i_face_u_normal[face_id][0]
                              + viscj[1][i]*i_face_u_normal[face_id][1]
                              + viscj[2][i]*i_face_u_normal[face_id][2])
                            *i_face_surf[face_id];
      }

      /* p in I" and J" */
      cs_real_t pipp = pi + bldfrp*(  grad[ii][0]*diippf[0]
                                    + grad[ii][1]*diippf[1]
                                    + grad[ii][2]*diippf[2]);
      cs_real_t pjpp = pj + bldfrp*(  grad[jj][0]*djjppf[0]
                                    + grad[jj][1]*djjppf[1]
                                    + grad[jj][2]*djjppf[2]);

      cs_real_t pir = pi/relaxp - (1.-relaxp)/relaxp * pia;
      cs_real_t pjr = pj/relaxp - (1.-relaxp)/relaxp * pja;

      /* pr in I" and J" */
      cs_real_t pippr = pir + bldfrp*(  grad[ii][0]*diippf[0]
                                      + grad[ii][1]*diippf[1]
                                      + grad[ii][2]*diippf[2]);
      cs_real_t pjppr = pjr + bldfrp*(  grad[jj][0]*djjppf[0]
                                      + grad[jj][1]*djjppf[1]
                                      + grad[jj][2]*djjppf[2]);


      cs_real_t fluxi = i_visc[face_id]*(pippr - pjpp);
      cs_real_t fluxj = i_visc[face_id]*(pipp - pjppr);

      if (ii < n_cells)
        cs_dispatch_sum(&rhs[ii], -fluxi, i_sum_type);
      if (jj < n_cells)
        cs_dispatch_sum(&rhs[jj], fluxj, i_sum_type);

    });

    /* Unsteady */
  }
  else {
    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t pi = _pvar[ii];
      cs_real_t pj = _pvar[jj];

      cs_real_t bldfrp = (cs_real_t) ircflp;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflp > 0)
        bldfrp = cs::max(df_limiter[jj], 0.);

      /* Recompute II" and JJ"
         ----------------------*/

      cs_real_t visci[3][3], viscj[3][3];

      visci[0][0] = viscce[ii][0];
      visci[1][1] = viscce[ii][1];
      visci[2][2] = viscce[ii][2];
      visci[1][0] = viscce[ii][3];
      visci[0][1] = viscce[ii][3];
      visci[2][1] = viscce[ii][4];
      visci[1][2] = viscce[ii][4];
      visci[2][0] = viscce[ii][5];
      visci[0][2] = viscce[ii][5];

      /* IF.Ki.S / ||Ki.S||^2 */
      cs_real_t fikdvi = weighf[face_id][0];

      cs_real_t diippf[3], djjppf[3];

      /* II" = IF + FI" */
      for (cs_lnum_t i = 0; i < 3; i++) {
        diippf[i] =   i_face_cog[face_id][i]-cell_cen[ii][i]
                    - fikdvi*(  visci[0][i]*i_face_u_normal[face_id][0]
                              + visci[1][i]*i_face_u_normal[face_id][1]
                              + visci[2][i]*i_face_u_normal[face_id][2])
                            *i_face_surf[face_id];
      }

      viscj[0][0] = viscce[jj][0];
      viscj[1][1] = viscce[jj][1];
      viscj[2][2] = viscce[jj][2];
      viscj[1][0] = viscce[jj][3];
      viscj[0][1] = viscce[jj][3];
      viscj[2][1] = viscce[jj][4];
      viscj[1][2] = viscce[jj][4];
      viscj[2][0] = viscce[jj][5];
      viscj[0][2] = viscce[jj][5];

      /* FJ.Kj.S / ||Kj.S||^2 */
      cs_real_t fjkdvi = weighf[face_id][1];

      /* JJ" = JF + FJ" */
      for (cs_lnum_t i = 0; i < 3; i++) {
        djjppf[i] =   i_face_cog[face_id][i]-cell_cen[jj][i]
                    + fjkdvi*( viscj[0][i]*i_face_u_normal[face_id][0]
                              + viscj[1][i]*i_face_u_normal[face_id][1]
                              + viscj[2][i]*i_face_u_normal[face_id][2])
                            *i_face_surf[face_id];
      }

      /* p in I" and J" */
      cs_real_t pipp = pi + bldfrp*(  grad[ii][0]*diippf[0]
                                    + grad[ii][1]*diippf[1]
                                    + grad[ii][2]*diippf[2]);
      cs_real_t pjpp = pj + bldfrp*(  grad[jj][0]*djjppf[0]
                                    + grad[jj][1]*djjppf[1]
                                    + grad[jj][2]*djjppf[2]);

      cs_real_t flux = i_visc[face_id]*(pipp -pjpp);

      if (ii < n_cells)
        cs_dispatch_sum(&rhs[ii], -thetap*flux, i_sum_type);
      if (jj < n_cells)
      cs_dispatch_sum(&rhs[jj], thetap*flux, i_sum_type);

    });
  }

  /* ======================================================================
     ---> Contribution from boundary faces
     ======================================================================*/

  /* Steady */
  if (idtvar < 0) {

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {

      cs_lnum_t ii = b_face_cells[face_id];

      cs_real_t pi = _pvar[ii];
      cs_real_t pia = pvara[ii];

      cs_real_t pir = pi/relaxp - (1.-relaxp)/relaxp*pia;

      cs_real_t bldfrp = (cs_real_t) ircflb;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflb > 0)
        bldfrp = cs::max(df_limiter[ii], 0.);

      /* Recompute II"
         --------------*/

      cs_real_t visci[3][3];

      visci[0][0] = viscce[ii][0];
      visci[1][1] = viscce[ii][1];
      visci[2][2] = viscce[ii][2];
      visci[1][0] = viscce[ii][3];
      visci[0][1] = viscce[ii][3];
      visci[2][1] = viscce[ii][4];
      visci[1][2] = viscce[ii][4];
      visci[2][0] = viscce[ii][5];
      visci[0][2] = viscce[ii][5];

      /* IF.Ki.S / ||Ki.S||^2 */
      cs_real_t fikdvi = weighb[face_id];

      cs_real_t diippf[3];

      /* II" = IF + FI" */
      for (cs_lnum_t i = 0; i < 3; i++) {
        diippf[i] = b_face_cog[face_id][i] - cell_cen[ii][i]
                    - fikdvi*(  visci[0][i]*b_face_u_normal[face_id][0]
                              + visci[1][i]*b_face_u_normal[face_id][1]
                              + visci[2][i]*b_face_u_normal[face_id][2])
                            *b_face_surf[face_id];
      }

      cs_real_t pippr = pir + bldfrp*(  grad[ii][0]*diippf[0]
                                      + grad[ii][1]*diippf[1]
                                      + grad[ii][2]*diippf[2]);

      cs_real_t pfacd = inc*cofafp[face_id] + cofbfp[face_id]*pippr;

      cs_real_t flux = b_visc[face_id]*pfacd;

      cs_dispatch_sum(&rhs[ii], -flux, b_sum_type);
    });

    /* Unsteady */
  }
  else {

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {

      cs_lnum_t ii = b_face_cells[face_id];

      cs_real_t pi = _pvar[ii];

      cs_real_t bldfrp = (cs_real_t) ircflb;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflb > 0)
        bldfrp = cs::max(df_limiter[ii], 0.);

      /* Recompute II"
         --------------*/

      cs_real_t visci[3][3];

      visci[0][0] = viscce[ii][0];
      visci[1][1] = viscce[ii][1];
      visci[2][2] = viscce[ii][2];
      visci[1][0] = viscce[ii][3];
      visci[0][1] = viscce[ii][3];
      visci[2][1] = viscce[ii][4];
      visci[1][2] = viscce[ii][4];
      visci[2][0] = viscce[ii][5];
      visci[0][2] = viscce[ii][5];

      /* IF.Ki.S / ||Ki.S||^2 */
      cs_real_t fikdvi = weighb[face_id];

      cs_real_t diippf[3];

      /* II" = IF + FI" */
      for (cs_lnum_t i = 0; i < 3; i++) {
        diippf[i] = b_face_cog[face_id][i] - cell_cen[ii][i]
                    - fikdvi*(  visci[0][i]*b_face_u_normal[face_id][0]
                              + visci[1][i]*b_face_u_normal[face_id][1]
                              + visci[2][i]*b_face_u_normal[face_id][2])
                            *b_face_surf[face_id];
      }

      cs_real_t pipp = pi + bldfrp*(  grad[ii][0]*diippf[0]
                                    + grad[ii][1]*diippf[1]
                                    + grad[ii][2]*diippf[2]);

      cs_real_t pfacd = inc*cofafp[face_id] + cofbfp[face_id]*pipp;

      cs_real_t flux = b_visc[face_id]*pfacd;
      cs_dispatch_sum(&rhs[ii], -thetap*flux, b_sum_type);

    });

    ctx.wait();

    /* The scalar is internal_coupled and an implicit contribution
     * is required */
    if (icoupl > 0) {

      /* Exchange pvar */
      CS_MALLOC(pvar_local, n_local, cs_real_t);
      cs_internal_coupling_exchange_by_cell_id(cpl,
                                               1, /* Dimension */
                                               _pvar,
                                               pvar_local);

      /* Exchange grad */
      CS_MALLOC(grad_local, n_local, cs_real_3_t);
      cs_internal_coupling_exchange_by_cell_id(cpl,
                                               3, /* Dimension */
                                               (const cs_real_t*)grad,
                                               (cs_real_t *)grad_local);

      /* Exchange viscce */
      CS_MALLOC(viscce_local, n_local, cs_real_6_t);
      cs_internal_coupling_exchange_by_cell_id(cpl,
                                               6, /* Dimension */
                                               (const cs_real_t*)viscce,
                                               (cs_real_t *)viscce_local);

      /* Exchange weighb */
      CS_MALLOC(weighb_local, n_local, cs_real_t);
      cs_internal_coupling_exchange_by_face_id(cpl,
                                               1, /* Dimension */
                                               (const cs_real_t*)weighb,
                                               (cs_real_t *)weighb_local);

      /* Exchange diffusion limiter */
      if (df_limiter != nullptr) {
        CS_MALLOC(df_limiter_local, n_local, cs_real_t);
        cs_internal_coupling_exchange_var(cpl,
                                          1, /* Dimension */
                                          df_limiter,
                                          df_limiter_local);
      }

      /* Flux contribution */
      for (cs_lnum_t jj = 0; jj < n_local; jj++) {
        cs_lnum_t face_id = faces_local[jj];
        cs_lnum_t ii = b_face_cells[face_id];

        cs_real_t pi = _pvar[ii];
        cs_real_t pj = pvar_local[jj];

        cs_real_t bldfrp = (cs_real_t) ircflb;
        /* Local limitation of the reconstruction */
        if (df_limiter != nullptr && ircflb > 0)
          bldfrp = cs::max(cs::min(df_limiter_local[ii],
                                   df_limiter[jj]),
                           0.);

        /* Recompute II" and JJ" */
        cs_real_t visci[3][3], viscj[3][3];
        cs_real_t diippf[3], djjppf[3];

        visci[0][0] = viscce[ii][0];
        visci[1][1] = viscce[ii][1];
        visci[2][2] = viscce[ii][2];
        visci[1][0] = viscce[ii][3];
        visci[0][1] = viscce[ii][3];
        visci[2][1] = viscce[ii][4];
        visci[1][2] = viscce[ii][4];
        visci[2][0] = viscce[ii][5];
        visci[0][2] = viscce[ii][5];

        /* IF.Ki.S / ||Ki.S||^2 */
        cs_real_t fikdvi = weighb[face_id];

        /* II" = IF + FI" */
        for (cs_lnum_t i = 0; i < 3; i++) {
          diippf[i] = b_face_cog[face_id][i]-cell_cen[ii][i]
                    - fikdvi*( visci[0][i]*b_face_u_normal[face_id][0]
                             + visci[1][i]*b_face_u_normal[face_id][1]
                             + visci[2][i]*b_face_u_normal[face_id][2])
                            *b_face_surf[face_id];
        }

        viscj[0][0] = viscce_local[jj][0];
        viscj[1][1] = viscce_local[jj][1];
        viscj[2][2] = viscce_local[jj][2];
        viscj[1][0] = viscce_local[jj][3];
        viscj[0][1] = viscce_local[jj][3];
        viscj[2][1] = viscce_local[jj][4];
        viscj[1][2] = viscce_local[jj][4];
        viscj[2][0] = viscce_local[jj][5];
        viscj[0][2] = viscce_local[jj][5];

        /* FJ.Kj.S / ||Kj.S||^2
         * weighb_local defined with vector JF and surface -S */
        cs_real_t fjkdvi = weighb_local[jj];

        /* JJ" = JF + FJ"
         *   */
        for (cs_lnum_t i = 0; i < 3; i++) {
          djjppf[i] =    b_face_cog[face_id][i]
                      - cell_cen[ii][i]-cpl->ci_cj_vect[jj][i]
                      + fjkdvi*(  viscj[0][i]*b_face_u_normal[face_id][0]
                                + viscj[1][i]*b_face_u_normal[face_id][1]
                                + viscj[2][i]*b_face_u_normal[face_id][2])
                              *b_face_surf[face_id];
        }

        /* p in I" and J" */
        cs_real_t pipp = pi + bldfrp*(  grad[ii][0]*diippf[0]
                                      + grad[ii][1]*diippf[1]
                                      + grad[ii][2]*diippf[2]);
        cs_real_t pjpp = pj + bldfrp*(  grad_local[jj][0]*djjppf[0]
                                      + grad_local[jj][1]*djjppf[1]
                                      + grad_local[jj][2]*djjppf[2]);

        /* Reproduce multiplication by i_visc[face_id] */
        cs_real_t flux = (pipp - pjpp) / (weighb[face_id] + weighb_local[jj]);

        rhs[ii] -= thetap*flux;

      }

      /* Remote data no longer needed */
      CS_FREE(pvar_local);
      if (df_limiter != nullptr)
        CS_FREE(df_limiter_local);
      CS_FREE(grad_local);
      CS_FREE(viscce_local);
      CS_FREE(weighb_local);
    }

  }

  ctx.wait();

  /* Free memory */
  CS_FREE(grad);
  CS_FREE(w2);
}

/*-----------------------------------------------------------------------------*/
/*!
 * \brief Add explicit part of the terms of diffusion by a left-multiplying
 * symmetric tensorial diffusivity for a transport equation of a vector field
 * \f$ \vect{\varia} \f$.
 *
 * More precisely, the right hand side \f$ \vect{Rhs} \f$ is updated as
 * follows:
 * \f[
 * \vect{Rhs} = \vect{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
 *      - \gradt_\fij \vect{\varia} \tens{\mu}_\fij  \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Remark:
 * if ivisep = 1, then we also take \f$ \mu \transpose{\gradt\vect{\varia}}
 * + \lambda \trace{\gradt\vect{\varia}} \f$, where \f$ \lambda \f$ is
 * the secondary viscosity, i.e. usually \f$ -\frac{2}{3} \mu \f$.
 *
 * Warning:
 * - \f$ \vect{Rhs} \f$ has already been initialized before calling the present
 *   function
 * - mind the sign minus
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          index of the current variable
 * \param[in]     eqp           equation parameters
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     ivisep        indicator to take \f$ \divv
 *                               \left(\mu \gradt \transpose{\vect{a}} \right)
 *                               -2/3 \grad\left( \mu \dive \vect{a} \right)\f$
 *                               - 1 take into account,
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     bc_coeffs_v   boundary condition structure for the variable
 * \param[in]     bc_coeffs_solve_v   sweep loop boundary conditions structure
 * \param[in]     i_visc        \f$ \tens{\mu}_\fij \dfrac{S_\fij}{\ipf\jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     i_secvis      secondary viscosity at interior faces
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_anisotropic_left_diffusion_vector
  (int                         idtvar,
   int                         f_id,
   const cs_equation_param_t   eqp,
   int                         inc,
   int                         ivisep,
   cs_real_3_t       *restrict pvar,
   const cs_real_3_t *restrict pvara,
   const cs_field_bc_coeffs_t *bc_coeffs_v,
   const cs_bc_coeffs_solve_t *bc_coeffs_solve_v,
   const cs_real_33_t          i_visc[],
   const cs_real_t             b_visc[],
   const cs_real_t             i_secvis[],
   cs_real_3_t       *restrict rhs)
{
  const int idiffp = eqp.idiff;
  const int ircflp = eqp.ircflu;
  const int ircflb = (ircflp > 0) ? eqp.b_diff_flux_rc : 0;
  const double relaxp = eqp.relaxv;
  const double thetap = eqp.theta;

  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_halo_t  *halo = m->halo;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict i_face_surf = fvq->i_face_surf;
  const cs_nreal_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_rreal_3_t *restrict diipf = fvq->diipf;
  const cs_rreal_3_t *restrict djjpf = fvq->djjpf;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;

  const int *bc_type = cs_glob_bc_type;

  const cs_real_3_t *val_f = (bc_coeffs_solve_v == nullptr) ?
                             (const cs_real_3_t *)bc_coeffs_v->val_f :
                             (const cs_real_3_t *)bc_coeffs_solve_v->val_f;

  const cs_real_3_t *val_f_d_lim = (bc_coeffs_solve_v == nullptr) ?
                                   (const cs_real_3_t *)bc_coeffs_v->val_f_d_lim :
                                   (const cs_real_3_t *)bc_coeffs_solve_v->val_f_d_lim;

  /* Local variables */

  cs_real_t *df_limiter = nullptr;

  char var_name[64];

  cs_real_33_t *gradv;

  cs_field_t *f = nullptr;

  /* Dispatch */
  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  /* 1. Initialization */

  /* Allocate work arrays */
  CS_MALLOC_HD(gradv, n_cells_ext, cs_real_33_t, cs_alloc_mode);

  /* Choose gradient type */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  cs_gradient_type_by_imrgra(eqp.d_gradient_r,
                             &gradient_type,
                             &halo_type);

  /* Handle cases where only the previous values (already synchronized)
     or current values are provided */

  if (pvar != nullptr && halo != nullptr) {
    cs_halo_sync_r(halo, halo_type, false, pvar);
  }
  if (pvara == nullptr)
    pvara = (const cs_real_3_t *)pvar;

  const cs_real_3_t  *restrict _pvar
    = (pvar != nullptr) ? (const cs_real_3_t  *)pvar : pvara;

  /* logging info */

  if (f_id != -1) {
    f = cs_field_by_id(f_id);

    int df_limiter_id =
      cs_field_get_key_int(f, cs_field_key_id("diffusion_limiter_id"));
    if (df_limiter_id > -1)
      df_limiter = cs_field_by_id(df_limiter_id)->val;

    snprintf(var_name, 63, "%s", f->name);
  }
  else
    strncpy(var_name, "[anisotropic left diffusion, vector]", 63);
  var_name[63] = '\0';

  /* 2. Compute the diffusive part with reconstruction technics */

  /* Compute the gradient of the current variable if needed */

  if (ircflp == 1 || ivisep == 1) {

    cs_gradient_vector_synced_input(var_name,
                                    gradient_type,
                                    halo_type,
                                    inc,
                                    eqp.nswrgr,
                                    eqp.verbosity,
                                    (cs_gradient_limit_t)(eqp.imligr),
                                    eqp.epsrgr,
                                    eqp.d_climgr,
                                    bc_coeffs_v,
                                    _pvar,
                                    val_f,
                                    nullptr, /* weighted gradient */
                                    gradv);

  }
  else {
    ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      for (cs_lnum_t isou = 0; isou < 3; isou++) {
        for (cs_lnum_t jsou = 0; jsou < 3; jsou++)
          gradv[cell_id][isou][jsou] = 0.;
      }
    });
  }

  /* ======================================================================
     ---> Contribution from interior faces
     ======================================================================*/

  /* Steady */
  if (idtvar < 0) {

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          cs_real_t pip[3], pjp[3], pipr[3], pjpr[3];

          cs_real_t bldfrp = (cs_real_t) ircflp;
          /* Local limitation of the reconstruction */
          if (df_limiter != nullptr && ircflp > 0)
            bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                             0.);

          /*-----------------
            X-Y-Z components, p = u, v, w */
          for (cs_lnum_t isou = 0; isou < 3; isou++) {

            cs_real_t dpvf[3];

            for (cs_lnum_t jsou = 0; jsou < 3; jsou++) {
              dpvf[jsou] = 0.5*(gradv[ii][isou][jsou] + gradv[jj][isou][jsou]);
            }

            cs_real_t pi  = _pvar[ii][isou];
            cs_real_t pj  = _pvar[jj][isou];

            cs_real_t pia = pvara[ii][isou];
            cs_real_t pja = pvara[jj][isou];

            /* reconstruction only if IRCFLP = 1 */
            pip[isou] = pi + bldfrp * (cs_math_3_dot_product(dpvf,
                                                             diipf[face_id]));
            pjp[isou] = pj + bldfrp * (cs_math_3_dot_product(dpvf,
                                                             djjpf[face_id]));

            pipr[isou] = pi /relaxp - (1.-relaxp)/relaxp * pia
                         + bldfrp * (cs_math_3_dot_product(dpvf,
                                                           diipf[face_id]));
            pjpr[isou] = pj /relaxp - (1.-relaxp)/relaxp * pja
                         + bldfrp * (cs_math_3_dot_product(dpvf,
                                                           djjpf[face_id]));

          }

          for (cs_lnum_t isou = 0; isou < 3; isou++) {

            cs_real_t fluxi =   i_visc[face_id][0][isou]*(pipr[0] - pjp[0])
                              + i_visc[face_id][1][isou]*(pipr[1] - pjp[1])
                              + i_visc[face_id][2][isou]*(pipr[2] - pjp[2]);
            cs_real_t fluxj =   i_visc[face_id][0][isou]*(pip[0] - pjpr[0])
                              + i_visc[face_id][1][isou]*(pip[1] - pjpr[1])
                              + i_visc[face_id][2][isou]*(pip[2] - pjpr[2]);

            if (ii < n_cells)
              rhs[ii][isou] = rhs[ii][isou] - fluxi;
            if (jj < n_cells)
              rhs[jj][isou] = rhs[jj][isou] + fluxj;

          }

        }
      }
    }

    /* Unsteady */
  }
  else {

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t pip[3], pjp[3];

      cs_real_t bldfrp = (cs_real_t) ircflp;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflp > 0)
        bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                         0.);

      /*-----------------
        X-Y-Z components, p = u, v, w */
      for (cs_lnum_t isou = 0; isou < 3; isou++) {
        cs_real_t dpvf[3];

        for (cs_lnum_t jsou = 0; jsou < 3; jsou++) {
          dpvf[jsou] = 0.5*(gradv[ii][isou][jsou] + gradv[jj][isou][jsou]);
        }

        cs_real_t pi = _pvar[ii][isou];
        cs_real_t pj = _pvar[jj][isou];

        pip[isou] = pi + bldfrp * (cs_math_3_dot_product(dpvf,
                                                         diipf[face_id]));
        pjp[isou] = pj + bldfrp * (cs_math_3_dot_product(dpvf,
                                                         djjpf[face_id]));
      }

      cs_real_t fluxi[3], fluxj[3];
      for (cs_lnum_t isou = 0; isou < 3; isou++) {
        cs_real_t flux =   i_visc[face_id][0][isou]*(pip[0] - pjp[0])
                         + i_visc[face_id][1][isou]*(pip[1] - pjp[1])
                         + i_visc[face_id][2][isou]*(pip[2] - pjp[2]);

        fluxj[isou] = thetap * flux;
        fluxi[isou] = - fluxj[isou];
      }

      if (ii < n_cells)
        cs_dispatch_sum<3>(rhs[ii], fluxi, i_sum_type);
      if (jj < n_cells)
        cs_dispatch_sum<3>(rhs[jj], fluxj, i_sum_type);
    });
  }

  /* ======================================================================
     ---> Contribution from boundary faces
     ======================================================================*/

  /* Steady */
  if (idtvar < 0) {

    const cs_real_3_t  *cofafv = (const cs_real_3_t  *)bc_coeffs_v->af;
    const cs_real_33_t *cofbfv = (const cs_real_33_t *)bc_coeffs_v->bf;

#   pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
    for (int t_id = 0; t_id < n_b_threads; t_id++) {
      for (cs_lnum_t face_id = b_group_index[t_id*2];
           face_id < b_group_index[t_id*2 + 1];
           face_id++) {

        cs_lnum_t cell_id = b_face_cells[face_id];

        cs_real_t bldfrp = (cs_real_t) ircflb;
        /* Local limitation of the reconstruction */
        if (df_limiter != nullptr && ircflb > 0)
          bldfrp = cs::max(df_limiter[cell_id], 0.);

        const cs_rreal_t *diipbv = diipb[face_id];
        cs_real_t pipr[3];

        for (cs_lnum_t k = 0; k < 3; k++) {
          cs_real_t pir  =   _pvar[cell_id][k]/relaxp
                           - (1.-relaxp)/relaxp*pvara[cell_id][k];

          pipr[k] = pir + bldfrp*(  gradv[cell_id][k][0]*diipbv[0]
                                  + gradv[cell_id][k][1]*diipbv[1]
                                  + gradv[cell_id][k][2]*diipbv[2]);

        }

        /* X-Y-Z components, p = u, v, w */
        for (cs_lnum_t i = 0; i < 3; i++) {

          cs_real_t pfacd = inc*cofafv[face_id][i];

          /*coefu and cofuf and b_visc are matrices */
          for (cs_lnum_t j = 0; j < 3; j++) {
            pfacd += cofbfv[face_id][i][j]*pipr[j];
          }

          rhs[cell_id][i] -= b_visc[face_id] * pfacd;

        } /* i */

      }
    }

    /* Unsteady */
  }
  else {

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
      cs_lnum_t cell_id = b_face_cells[face_id];

      cs_real_t fluxi[3];
      /* X-Y-Z components, p = u, v, w */
      for (cs_lnum_t i = 0; i < 3; i++) {
        cs_real_t flux = thetap * b_visc[face_id] * val_f_d_lim[face_id][i];
        fluxi[i] = -flux;
      }

      cs_dispatch_sum<3>(rhs[cell_id], fluxi, b_sum_type);
    });
  } /* idtvar */

  /* 3. Computation of the transpose grad(vel) term and grad(-2/3 div(vel)) */
  short *bndcel = nullptr;

  if (ivisep == 1 && idiffp == 1) {

    /* We do not know what condition to put in the inlets and the outlets, so we
       assume that there is an equilibrium. Moreover, cells containing a coupled
       are removed. */

    /* Allocate a temporary array */
    CS_MALLOC_HD(bndcel, n_cells_ext, short, cs_alloc_mode);

    ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      bndcel[cell_id] = 1;
    });

    ctx.parallel_for(m->n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      int ityp = bc_type[face_id];
      if (   ityp == CS_OUTLET
          || ityp == CS_INLET
          || ityp == CS_FREE_INLET
          || ityp == CS_CONVECTIVE_INLET
          || ityp == CS_COUPLED_FD)
        bndcel[b_face_cells[face_id]] = 0;
    });

    /* ---> Interior faces */

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t pnd = weight[face_id];
      cs_real_t grdtrv =       pnd  * cs_math_33_trace(gradv[ii])
                         + (1.-pnd) * cs_math_33_trace(gradv[jj]);

      cs_real_t flux0 = i_secvis[face_id] * grdtrv * i_face_surf[face_id];

      cs_real_t ipjp[3];
      const cs_nreal_t *n = i_face_u_normal[face_id];

      /* I'J' */
      for (cs_lnum_t i = 0; i < 3; i++)
        ipjp[i] = n[i] * i_dist[face_id];

      cs_real_t fluxi[3], fluxj[3];
      for (cs_lnum_t i = 0; i < 3; i++) {
        cs_real_t flux = flux0*i_face_u_normal[face_id][i];

        /* We need to compute (K grad(u)^T) .I'J'
           which is equal to I'J' . (grad(u) . K^T)
           But: (I'J' . (grad(u) . K^T))_i = I'J'_k grad(u)_kj K_ij
        */

        for (cs_lnum_t j = 0; j < 3; j++) {
          for (cs_lnum_t k = 0; k < 3; k++) {
            flux +=   ipjp[k]
                    * (pnd*gradv[ii][k][j]+(1.-pnd)*gradv[jj][k][j])
                    * i_visc[face_id][i][j];
          }
        }

        fluxi[i] =  flux * bndcel[ii];
        fluxj[i] = -flux * bndcel[jj];
      }

      if (ii < n_cells)
        cs_dispatch_sum<3>(rhs[ii], fluxi, i_sum_type);
      if (jj < n_cells)
        cs_dispatch_sum<3>(rhs[jj], fluxj, i_sum_type);
    });

    /* ---> Boundary FACES
       the whole flux term of the stress tensor is already taken into account
       (so, no corresponding term in forbr)
       TODO in theory we should take the normal component into account (the
       tangential one is modeled by the wall law) */
  }

  ctx.wait();

  /* Free memory */
  CS_FREE(bndcel);
  CS_FREE(gradv);
}

/*-----------------------------------------------------------------------------*/
/*!
 * \brief Add explicit part of the terms of diffusion by a right-multiplying
 * symmetric tensorial diffusivity for a transport equation of a vector field
 * \f$ \vect{\varia} \f$.
 *
 * More precisely, the right hand side \f$ \vect{Rhs} \f$ is updated as
 * follows:
 * \f[
 * \vect{Rhs} = \vect{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
 *      - \gradt_\fij \vect{\varia} \tens{\mu}_\fij  \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ \vect{Rhs} \f$ has already been initialized before calling the present
 *   function
 * - mind the sign minus
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          index of the current variable
 * \param[in]     eqp           equation parameters
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     bc_coeffs_v   boundary condition structure for the variable
 * \param[in]     bc_coeffs_solve_v   sweep loop boundary conditions structure
 * \param[in]     i_visc        \f$ \tens{\mu}_\fij \dfrac{S_\fij}{\ipf\jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]     weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_anisotropic_right_diffusion_vector
  (int                          idtvar,
   int                          f_id,
   const cs_equation_param_t    eqp,
   int                          inc,
   cs_real_3_t        *restrict pvar,
   const cs_real_3_t  *restrict pvara,
   const cs_field_bc_coeffs_t  *bc_coeffs_v,
   const cs_bc_coeffs_solve_t  *bc_coeffs_solve_v,
   const cs_real_t              i_visc[],
   const cs_real_t              b_visc[],
   cs_real_6_t        *restrict viscel,
   const cs_real_2_t            weighf[],
   const cs_real_t              weighb[],
   cs_real_3_t        *restrict rhs)
{
  const int ircflp = eqp.ircflu;
  const int ircflb = (ircflp > 0) ? eqp.b_diff_flux_rc : 0;
  const int icoupl = eqp.icoupl;
  const double relaxp = eqp.relaxv;
  const double thetap = eqp.theta;

  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_halo_t  *halo = m->halo;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const int n_b_threads = m->b_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_nreal_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;
  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
  const cs_real_t *restrict i_face_surf = fvq->i_face_surf;
  const cs_real_t *restrict b_face_surf = fvq->b_face_surf;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog = fvq->b_face_cog;

  const cs_real_3_t *val_f = (bc_coeffs_solve_v == nullptr) ?
                             (const cs_real_3_t *)bc_coeffs_v->val_f :
                             (const cs_real_3_t *)bc_coeffs_solve_v->val_f;

  const cs_real_3_t  *cofafv = (const cs_real_3_t  *)bc_coeffs_v->af;
  const cs_real_33_t *cofbfv = (const cs_real_33_t *)bc_coeffs_v->bf;

  /* Dispatch */
  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  /* Internal coupling variables */
  cs_real_3_t *pvar_local = nullptr;
  cs_real_33_t *grad_local = nullptr;
  cs_real_t *df_limiter_local = nullptr;
  cs_real_6_t *viscce_local = nullptr;
  cs_real_t *weighb_local = nullptr;
  const cs_lnum_t *faces_local = nullptr, *faces_distant = nullptr;
  cs_lnum_t n_local = 0, n_distant = 0;
  int coupling_id = -1;
  cs_internal_coupling_t *cpl = nullptr;

  /* Local variables */

  cs_real_t *df_limiter = nullptr;

  char var_name[64];

  cs_real_6_t *viscce;
  cs_real_33_t *grad;

  cs_field_t *f = nullptr;

  /* 1. Initialization */

  viscce = nullptr;

  /* Allocate work arrays */
  CS_MALLOC_HD(grad, n_cells_ext, cs_real_33_t, cs_alloc_mode);

  /* Choose gradient type */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  cs_gradient_type_by_imrgra(eqp.d_gradient_r,
                             &gradient_type,
                             &halo_type);

  /* Handle cases where only the previous values (already synchronized)
     or current values are provided */

  if (pvar != nullptr && halo != nullptr) {
    cs_halo_sync_r(halo, halo_type, false, pvar);
  }
  if (pvara == nullptr)
    pvara = (const cs_real_3_t *)pvar;

  const cs_real_3_t  *restrict _pvar
    = (pvar != nullptr) ? (const cs_real_3_t  *)pvar : pvara;

  /* logging info */

  if (f_id != -1) {
    f = cs_field_by_id(f_id);

    int df_limiter_id =
      cs_field_get_key_int(f, cs_field_key_id("diffusion_limiter_id"));
    if (df_limiter_id > -1)
      df_limiter = cs_field_by_id(df_limiter_id)->val;

    snprintf(var_name, 63, "%s", f->name);
  }
  else
    strncpy(var_name, "[anisotropic right diffusion, vector]", 63);
  var_name[63] = '\0';

  viscce = viscel;

  if (icoupl > 0) {
    assert(f_id != -1);
    const int coupling_key_id = cs_field_key_id("coupling_entity");
    coupling_id = cs_field_get_key_int(f, coupling_key_id);
    cpl = cs_internal_coupling_by_id(coupling_id);
    cs_internal_coupling_coupled_faces(cpl,
                                       &n_local,
                                       &faces_local,
                                       &n_distant,
                                       &faces_distant);
  }

  /* 2. Compute the diffusive part with reconstruction technics */

  /* Compute the gradient of the current variable if needed */

  if (ircflp == 1) {

    cs_gradient_vector_synced_input(var_name,
                                    gradient_type,
                                    halo_type,
                                    inc,
                                    eqp.nswrgr,
                                    eqp.verbosity,
                                    (cs_gradient_limit_t)(eqp.imligr),
                                    eqp.epsrgr,
                                    eqp.d_climgr,
                                    bc_coeffs_v,
                                    _pvar,
                                    val_f,
                                    nullptr, /* weighted gradient */
                                    grad);

  }
  else {
    ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      for (cs_lnum_t isou = 0; isou < 3; isou++) {
        for (cs_lnum_t jsou = 0; jsou < 3; jsou++)
          grad[cell_id][isou][jsou] = 0.;
      }
    });
  }

  /* ======================================================================
     ---> Contribution from interior faces
     ======================================================================*/

  /* Steady */
  if (idtvar < 0) {

    for (int g_id = 0; g_id < n_i_groups; g_id++) {
#     pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t face_id = i_group_index[(t_id*n_i_groups + g_id)*2];
             face_id < i_group_index[(t_id*n_i_groups + g_id)*2 + 1];
             face_id++) {

          cs_lnum_t ii = i_face_cells[face_id][0];
          cs_lnum_t jj = i_face_cells[face_id][1];

          cs_real_t visci[3][3], viscj[3][3];
          cs_real_t diippf[3], djjppf[3], pipp[3], pjpp[3];
          cs_real_t  pir[3], pjr[3], pippr[3], pjppr[3];
          cs_real_t  pi[3], pj[3], pia[3], pja[3];

          for (cs_lnum_t isou = 0; isou < 3; isou++) {
            pi[isou] = _pvar[ii][isou];
            pj[isou] = _pvar[jj][isou];
            pia[isou] = pvara[ii][isou];
            pja[isou] = pvara[jj][isou];
          }

          cs_real_t bldfrp = (cs_real_t) ircflp;
          /* Local limitation of the reconstruction */
          if (df_limiter != nullptr && ircflp > 0)
            bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                             0.);

          /* Recompute II' and JJ' at this level */

          visci[0][0] = viscce[ii][0];
          visci[1][1] = viscce[ii][1];
          visci[2][2] = viscce[ii][2];
          visci[1][0] = viscce[ii][3];
          visci[0][1] = viscce[ii][3];
          visci[2][1] = viscce[ii][4];
          visci[1][2] = viscce[ii][4];
          visci[2][0] = viscce[ii][5];
          visci[0][2] = viscce[ii][5];

          /* IF.Ki.S / ||Ki.S||^2 */
          cs_real_t fikdvi_s = weighf[face_id][0] * i_face_surf[face_id];

          /* II" = IF + FI" */
          for (cs_lnum_t i = 0; i < 3; i++) {
            diippf[i] =   i_face_cog[face_id][i]-cell_cen[ii][i]
                        - fikdvi_s*(  visci[0][i]*i_face_u_normal[face_id][0]
                                    + visci[1][i]*i_face_u_normal[face_id][1]
                                    + visci[2][i]*i_face_u_normal[face_id][2]);
          }

          viscj[0][0] = viscce[jj][0];
          viscj[1][1] = viscce[jj][1];
          viscj[2][2] = viscce[jj][2];
          viscj[1][0] = viscce[jj][3];
          viscj[0][1] = viscce[jj][3];
          viscj[2][1] = viscce[jj][4];
          viscj[1][2] = viscce[jj][4];
          viscj[2][0] = viscce[jj][5];
          viscj[0][2] = viscce[jj][5];

          /* FJ.Kj.S / ||Kj.S||^2 */
          cs_real_t fjkdvi_s = weighf[face_id][1] * i_face_surf[face_id];

          /* JJ" = JF + FJ" */
          for (cs_lnum_t i = 0; i < 3; i++) {
            djjppf[i] = i_face_cog[face_id][i]-cell_cen[jj][i]
                      + fjkdvi_s*(  viscj[0][i]*i_face_u_normal[face_id][0]
                                  + viscj[1][i]*i_face_u_normal[face_id][1]
                                  + viscj[2][i]*i_face_u_normal[face_id][2]);
          }

          for (cs_lnum_t isou = 0; isou < 3; isou++) {
            /* p in I" and J" */
            pipp[isou] = pi[isou] + bldfrp*( grad[ii][isou][0]*diippf[0]
                                           + grad[ii][isou][1]*diippf[1]
                                           + grad[ii][isou][2]*diippf[2]);
            pjpp[isou] = pj[isou] + bldfrp*( grad[jj][isou][0]*djjppf[0]
                                           + grad[jj][isou][1]*djjppf[1]
                                           + grad[jj][isou][2]*djjppf[2]);

            pir[isou] = pi[isou]/relaxp - (1.-relaxp)/relaxp * pia[isou];
            pjr[isou] = pj[isou]/relaxp - (1.-relaxp)/relaxp * pja[isou];


            /* pr in I" and J" */
            pippr[isou] = pir[isou] + bldfrp*( grad[ii][isou][0]*diippf[0]
                                             + grad[ii][isou][1]*diippf[1]
                                             + grad[ii][isou][2]*diippf[2]);
            pjppr[isou] = pjr[isou] + bldfrp*( grad[jj][isou][0]*djjppf[0]
                                             + grad[jj][isou][1]*djjppf[1]
                                             + grad[jj][isou][2]*djjppf[2]);

            cs_real_t fluxi = i_visc[face_id]*(pippr[isou] - pjpp[isou]);
            cs_real_t fluxj = i_visc[face_id]*(pipp[isou] - pjppr[isou]);

            if (ii < n_cells)
              rhs[ii][isou] -= fluxi;
            if (jj < n_cells)
              rhs[jj][isou] += fluxj;
          }
        }
      }
    }

    /* Unsteady */
  }
  else {

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t visci[3][3], viscj[3][3];
      cs_real_t diippf[3], djjppf[3], pipp[3], pjpp[3];
      cs_real_t pi[3], pj[3];

      for (cs_lnum_t isou = 0; isou < 3; isou++) {
        pi[isou] = _pvar[ii][isou];
        pj[isou] = _pvar[jj][isou];
      }

      cs_real_t bldfrp = (cs_real_t) ircflp;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflp > 0)
        bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                         0.);

      /* Recompute II' and JJ' at this level */

      visci[0][0] = viscce[ii][0];
      visci[1][1] = viscce[ii][1];
      visci[2][2] = viscce[ii][2];
      visci[1][0] = viscce[ii][3];
      visci[0][1] = viscce[ii][3];
      visci[2][1] = viscce[ii][4];
      visci[1][2] = viscce[ii][4];
      visci[2][0] = viscce[ii][5];
      visci[0][2] = viscce[ii][5];

      /* IF.Ki.S / ||Ki.S||^2 */
      cs_real_t fikdvi_s = weighf[face_id][0] * i_face_surf[face_id];

      /* II" = IF + FI" */
      for (cs_lnum_t i = 0; i < 3; i++) {
        diippf[i] =   i_face_cog[face_id][i]-cell_cen[ii][i]
                    - fikdvi_s*(  visci[0][i]*i_face_u_normal[face_id][0]
                                + visci[1][i]*i_face_u_normal[face_id][1]
                                + visci[2][i]*i_face_u_normal[face_id][2]);
      }

      viscj[0][0] = viscce[jj][0];
      viscj[1][1] = viscce[jj][1];
      viscj[2][2] = viscce[jj][2];
      viscj[1][0] = viscce[jj][3];
      viscj[0][1] = viscce[jj][3];
      viscj[2][1] = viscce[jj][4];
      viscj[1][2] = viscce[jj][4];
      viscj[2][0] = viscce[jj][5];
      viscj[0][2] = viscce[jj][5];

      /* FJ.Kj.S / ||Kj.S||^2 */
      cs_real_t fjkdvi_s = weighf[face_id][1] * i_face_surf[face_id];

      /* JJ" = JF + FJ" */
      for (cs_lnum_t i = 0; i < 3; i++) {
        djjppf[i] =   i_face_cog[face_id][i]-cell_cen[jj][i]
                    + fjkdvi_s*(  viscj[0][i]*i_face_u_normal[face_id][0]
                                + viscj[1][i]*i_face_u_normal[face_id][1]
                                + viscj[2][i]*i_face_u_normal[face_id][2]);
      }

      cs_real_t fluxi[3], fluxj[3];

      for (cs_lnum_t isou = 0; isou < 3; isou++) {
        /* p in I" and J" */
        pipp[isou] =  pi[isou]
                    + bldfrp * cs_math_3_dot_product(grad[ii][isou], diippf);

        pjpp[isou] =  pj[isou]
                    + bldfrp * cs_math_3_dot_product(grad[jj][isou], djjppf);

        cs_real_t flux = i_visc[face_id]*(pipp[isou] - pjpp[isou]);

        fluxj[isou] =   thetap * flux;
        fluxi[isou] = - fluxj[isou];
      }

      if (ii < n_cells)
        cs_dispatch_sum<3>(rhs[ii], fluxi, i_sum_type);
      if (jj < n_cells)
        cs_dispatch_sum<3>(rhs[jj], fluxj, i_sum_type);
    });
  }

  /* ======================================================================
     ---> Contribution from boundary faces
     ======================================================================*/

  /* Steady */
  if (idtvar < 0) {

#   pragma omp parallel for if(m->n_b_faces > CS_THR_MIN)
    for (int t_id = 0; t_id < n_b_threads; t_id++) {
      for (cs_lnum_t face_id = b_group_index[t_id*2];
           face_id < b_group_index[t_id*2 + 1];
           face_id++) {

        cs_lnum_t ii = b_face_cells[face_id];

        cs_real_t pi[3], pia[3], pir[3], pippr[3];
        cs_real_t visci[3][3];
        cs_real_t diippf[3];

        for (cs_lnum_t isou = 0; isou < 3; isou++) {
          pi[isou] = pvar[ii][isou];
          pia[isou] = pvara[ii][isou];
          pir[isou] = pi[isou]/relaxp - (1.-relaxp)/relaxp*pia[isou];
        }

        cs_real_t bldfrp = (cs_real_t) ircflb;
        /* Local limitation of the reconstruction */
        if (df_limiter != nullptr && ircflb > 0)
          bldfrp = cs::max(df_limiter[ii], 0.);

        /* Recompute II"
           --------------*/

        visci[0][0] = viscce[ii][0];
        visci[1][1] = viscce[ii][1];
        visci[2][2] = viscce[ii][2];
        visci[1][0] = viscce[ii][3];
        visci[0][1] = viscce[ii][3];
        visci[2][1] = viscce[ii][4];
        visci[1][2] = viscce[ii][4];
        visci[2][0] = viscce[ii][5];
        visci[0][2] = viscce[ii][5];

        /* IF.Ki.S / ||Ki.S||^2 */
        cs_real_t fikdvi_s = weighb[face_id] * b_face_surf[face_id];

        /* II" = IF + FI" */
        for (cs_lnum_t i = 0; i < 3; i++) {
          diippf[i] = b_face_cog[face_id][i] - cell_cen[ii][i]
                      - fikdvi_s*(  visci[0][i]*b_face_u_normal[face_id][0]
                                  + visci[1][i]*b_face_u_normal[face_id][1]
                                  + visci[2][i]*b_face_u_normal[face_id][2]);
        }
        for (cs_lnum_t isou = 0; isou < 3; isou++) {
          pippr[isou] = pir[isou] + bldfrp*(  grad[ii][isou][0]*diippf[0]
                                            + grad[ii][isou][1]*diippf[1]
                                            + grad[ii][isou][2]*diippf[2]);
        }
        for (cs_lnum_t isou = 0; isou < 3; isou++) {
          cs_real_t pfacd = inc*cofafv[face_id][isou];
          for (cs_lnum_t jsou = 0; jsou < 3; jsou++) {
            pfacd += cofbfv[face_id][isou][jsou]*pippr[jsou];
          }

          cs_real_t flux = b_visc[face_id]*pfacd;
          rhs[ii][isou] -= flux;

        } /* isou */

      }
    }

    /* Unsteady */
  }
  else {

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = b_face_cells[face_id];

      cs_real_t visci[3][3];
      cs_real_t diippf[3], pi[3], pipp[3];

      for (cs_lnum_t isou = 0; isou < 3; isou++) {
        pi[isou] = pvar[ii][isou];
      }

      cs_real_t bldfrp = (cs_real_t) ircflb;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflb > 0)
        bldfrp = cs::max(df_limiter[ii], 0.);

      /* Recompute II"
         --------------*/

      visci[0][0] = viscce[ii][0];
      visci[1][1] = viscce[ii][1];
      visci[2][2] = viscce[ii][2];
      visci[1][0] = viscce[ii][3];
      visci[0][1] = viscce[ii][3];
      visci[2][1] = viscce[ii][4];
      visci[1][2] = viscce[ii][4];
      visci[2][0] = viscce[ii][5];
      visci[0][2] = viscce[ii][5];

      /* IF.Ki.S / ||Ki.S||^2 */
      cs_real_t fikdvi_s = weighb[face_id] * b_face_surf[face_id];

      /* II" = IF + FI" */
      for (cs_lnum_t i = 0; i < 3; i++) {
        diippf[i] = b_face_cog[face_id][i] - cell_cen[ii][i]
                    - fikdvi_s*(  visci[0][i]*b_face_u_normal[face_id][0]
                                + visci[1][i]*b_face_u_normal[face_id][1]
                                + visci[2][i]*b_face_u_normal[face_id][2]);
      }
      for (cs_lnum_t isou = 0; isou < 3; isou++) {
        pipp[isou] =  pi[isou]
                    + bldfrp * cs_math_3_dot_product(grad[ii][isou], diippf);
      }

      cs_real_t fluxi[3];

      for (cs_lnum_t isou = 0; isou < 3; isou++) {
        cs_real_t pfacd = inc*cofafv[face_id][isou];

        for (cs_lnum_t jsou = 0; jsou < 3; jsou++) {
          pfacd += cofbfv[face_id][isou][jsou]*pipp[jsou];
        }

        cs_real_t flux = b_visc[face_id]*pfacd;
        fluxi[isou] = -thetap * flux;
      } /* isou */

      cs_dispatch_sum<3>(rhs[ii], fluxi, b_sum_type);
    });

    ctx.wait();

    /* The vector is internal_coupled and an implicit contribution
     * is required */
    if (icoupl > 0) {

      /* Exchange pvar */
      CS_MALLOC(pvar_local, n_local, cs_real_3_t);
      cs_internal_coupling_exchange_by_cell_id(cpl,
                                               3, /* Dimension */
                                               (const cs_real_t*)_pvar,
                                               (cs_real_t *)pvar_local);

      /* Exchange grad */
      CS_MALLOC(grad_local, n_local, cs_real_33_t);
      cs_internal_coupling_exchange_by_cell_id(cpl,
                                               9, /* Dimension */
                                               (const cs_real_t*)grad,
                                               (cs_real_t *)grad_local);

      /* Exchange viscce */
      CS_MALLOC(viscce_local, n_local, cs_real_6_t);
      cs_internal_coupling_exchange_by_cell_id(cpl,
                                               6, /* Dimension */
                                               (const cs_real_t*)viscce,
                                               (cs_real_t *)viscce_local);

      /* Exchange weighb */
      CS_MALLOC(weighb_local, n_local, cs_real_t);
      cs_internal_coupling_exchange_by_face_id(cpl,
                                               1, /* Dimension */
                                               (const cs_real_t*)weighb,
                                               (cs_real_t *)weighb_local);

      /* Exchange diffusion limiter */
      if (df_limiter != nullptr) {
        CS_MALLOC(df_limiter_local, n_local, cs_real_t);
        cs_internal_coupling_exchange_var(cpl,
                                          1, /* Dimension */
                                          df_limiter,
                                          df_limiter_local);
      }

      /* Flux contribution */
      for (cs_lnum_t jj = 0; jj < n_local; jj++) {
        cs_lnum_t face_id = faces_local[jj];
        cs_lnum_t ii = b_face_cells[face_id];

        cs_real_t pipp[3], pjpp[3];
        cs_real_t pi[3], pj[3];

        for (cs_lnum_t isou = 0; isou < 3; isou++) {
          pi[isou] = _pvar[ii][isou];
          pj[isou] = pvar_local[jj][isou];
        }

        cs_real_t bldfrp = (cs_real_t) ircflb;
        /* Local limitation of the reconstruction */
        if (df_limiter != nullptr && ircflb > 0)
          bldfrp = cs::max(cs::min(df_limiter_local[ii],
                                   df_limiter[jj]),
                           0.);

        /* Recompute II" and JJ" */
        cs_real_t visci[3][3], viscj[3][3];
        cs_real_t diippf[3], djjppf[3];

        visci[0][0] = viscce[ii][0];
        visci[1][1] = viscce[ii][1];
        visci[2][2] = viscce[ii][2];
        visci[1][0] = viscce[ii][3];
        visci[0][1] = viscce[ii][3];
        visci[2][1] = viscce[ii][4];
        visci[1][2] = viscce[ii][4];
        visci[2][0] = viscce[ii][5];
        visci[0][2] = viscce[ii][5];

        /* IF.Ki.S / ||Ki.S||^2 */
        cs_real_t fikdvi_s = weighb[face_id] * b_face_surf[face_id];

        /* II" = IF + FI" */
        for (cs_lnum_t i = 0; i < 3; i++) {
          diippf[i] = b_face_cog[face_id][i]-cell_cen[ii][i]
                    - fikdvi_s*(  visci[0][i]*b_face_u_normal[face_id][0]
                                + visci[1][i]*b_face_u_normal[face_id][1]
                                + visci[2][i]*b_face_u_normal[face_id][2]);
        }

        viscj[0][0] = viscce_local[jj][0];
        viscj[1][1] = viscce_local[jj][1];
        viscj[2][2] = viscce_local[jj][2];
        viscj[1][0] = viscce_local[jj][3];
        viscj[0][1] = viscce_local[jj][3];
        viscj[2][1] = viscce_local[jj][4];
        viscj[1][2] = viscce_local[jj][4];
        viscj[2][0] = viscce_local[jj][5];
        viscj[0][2] = viscce_local[jj][5];

        /* FJ.Kj.S / ||Kj.S||^2
         * weighb_local defined with vector JF and surface -S */
        cs_real_t fjkdvi_s = weighb_local[jj] * b_face_surf[face_id];

        /* JJ" = JF + FJ"
         *   */
        for (cs_lnum_t i = 0; i < 3; i++) {
          djjppf[i] =   b_face_cog[face_id][i]
                      - cell_cen[ii][i] - cpl->ci_cj_vect[jj][i]
                       + fjkdvi_s*(  viscj[0][i]*b_face_u_normal[face_id][0]
                                   + viscj[1][i]*b_face_u_normal[face_id][1]
                                   + viscj[2][i]*b_face_u_normal[face_id][2]);
        }

        for (cs_lnum_t isou = 0; isou < 3; isou++) {
          /* p in I" and J" */
          pipp[isou] = pi[isou] + bldfrp*(  grad[ii][isou][0]*diippf[0]
                                          + grad[ii][isou][1]*diippf[1]
                                          + grad[ii][isou][2]*diippf[2]);
          pjpp[isou] = pj[isou] + bldfrp*(  grad_local[jj][isou][0]*djjppf[0]
                                          + grad_local[jj][isou][1]*djjppf[1]
                                          + grad_local[jj][isou][2]*djjppf[2]);

          /* Reproduce multiplication by i_visc[face_id] */
          cs_real_t flux =   (pipp[isou] - pjpp[isou])
                           / (weighb[face_id] + weighb_local[jj]);

          rhs[ii][isou] -= thetap*flux;
        }

      }

      /* Remote data no longer needed */
      CS_FREE(pvar_local);
      if (df_limiter != nullptr)
        CS_FREE(df_limiter_local);
      CS_FREE(grad_local);
      CS_FREE(viscce_local);
      CS_FREE(weighb_local);

    }

  } /* idtvar */

  /* Free memory */
  CS_FREE(grad);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the diffusion terms with a symmetric tensor
 * diffusivity for a transport equation of a scalar field \f$ \varia \f$.
 *
 * More precisely, the right hand side \f$ Rhs \f$ is updated as
 * follows:
 * \f[
 * Rhs = Rhs - \sum_{\fij \in \Facei{\celli}}      \left(
 *      - \tens{\mu}_\fij \gradv_\fij \varia \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * Warning:
 * - \f$ Rhs \f$ has already been initialized before
 *   calling cs_anisotropic_diffusion_scalar!
 * - mind the sign minus
 *
 * \param[in]     idtvar        indicator of the temporal scheme
 * \param[in]     f_id          index of the current variable
 * \param[in]     eqp           equation parameters
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     pvara         solved variable (previous time step)
 * \param[in]     bc_coeffs_ts  boundary condition structure for the variable
 * \param[in]     bc_coeffs_solve_ts   sweep loop boundary conditions structure
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]     weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_anisotropic_diffusion_tensor(int                          idtvar,
                                int                          f_id,
                                const cs_equation_param_t    eqp,
                                int                          inc,
                                cs_real_6_t        *restrict pvar,
                                const cs_real_6_t  *restrict pvara,
                                const cs_field_bc_coeffs_t  *bc_coeffs_ts,
                                const cs_bc_coeffs_solve_t  *bc_coeffs_solve_ts,
                                const cs_real_t              i_visc[],
                                const cs_real_t              b_visc[],
                                cs_real_6_t        *restrict viscel,
                                const cs_real_2_t            weighf[],
                                const cs_real_t              weighb[],
                                cs_real_6_t        *restrict rhs)
{
  const cs_real_6_t *val_f = (bc_coeffs_solve_ts == nullptr) ?
                             (const cs_real_6_t *)bc_coeffs_ts->val_f :
                             (const cs_real_6_t *)bc_coeffs_solve_ts->val_f;

  const cs_real_6_t  *cofaf = (const cs_real_6_t  *)bc_coeffs_ts->af;
  const cs_real_66_t *cofbf = (const cs_real_66_t *)bc_coeffs_ts->bf;

  const int ircflp = eqp.ircflu;
  const int ircflb = (ircflp > 0) ? eqp.b_diff_flux_rc : 0;
  const double relaxp = eqp.relaxv;
  const double thetap = eqp.theta;

  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_halo_t  *halo = m->halo;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_real_t *restrict i_face_surf = fvq->i_face_surf;
  const cs_real_t *restrict b_face_surf = fvq->b_face_surf;
  const cs_nreal_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;
  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog = fvq->b_face_cog;

  /* Dispatch */
  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  /* Local variables */

  cs_real_t *df_limiter = nullptr;

  char var_name[64];

  cs_real_6_t *viscce = nullptr;
  cs_real_6_t *w2 = nullptr;

  cs_field_t *f = nullptr;

  /* Allocate work arrays */
  cs_real_63_t *grad;
  CS_MALLOC_HD(grad, n_cells_ext, cs_real_63_t, cs_alloc_mode);

  /* Choose gradient type */
  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  cs_gradient_type_by_imrgra(eqp.d_gradient_r,
                             &gradient_type,
                             &halo_type);

  /* Handle cases where only the previous values (already synchronized)
     or current values are provided */

  if (pvar != nullptr && m->halo != nullptr) {
    cs_halo_sync(m->halo, halo_type, false, pvar);
  }
  if (pvara == nullptr)
    pvara = (const cs_real_6_t *)pvar;

  const cs_real_6_t  *restrict _pvar
    = (pvar != nullptr) ? (const cs_real_6_t  *)pvar : pvara;

  /* Logging info */

  if (f_id != -1) {
    f = cs_field_by_id(f_id);

    int df_limiter_id
      = cs_field_get_key_int(f, cs_field_key_id("diffusion_limiter_id"));
    if (df_limiter_id > -1)
      df_limiter = cs_field_by_id(df_limiter_id)->val;

    snprintf(var_name, 63, "%s", f->name);
  }
  else
    strncpy(var_name, "[anisotropic diffusion, tensor]", 63);
  var_name[63] = '\0';

  /* Porosity fields */
  cs_field_t *fporo = cs_field_by_name_try("porosity");
  cs_field_t *ftporo = cs_field_by_name_try("tensorial_porosity");

  cs_real_t *porosi = nullptr;
  cs_real_6_t *porosf = nullptr;

  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2) {
    porosi = fporo->val;
    if (ftporo != nullptr) {
      porosf = (cs_real_6_t *)ftporo->val;
    }
  }

  /* Without porosity */
  if (porosi == nullptr) {
    viscce = viscel;

    /* With porosity */
  }
  else if (porosi != nullptr && porosf == nullptr) {
    CS_MALLOC_HD(w2, n_cells_ext, cs_real_6_t, cs_alloc_mode);
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      for (cs_lnum_t isou = 0; isou < 6; isou++) {
        w2[cell_id][isou] = porosi[cell_id]*viscel[cell_id][isou];
      }
    });
    ctx.wait();
    viscce = w2;
    /* With tensorial porosity */
  }
  else if (porosi != nullptr && porosf != nullptr) {
    CS_MALLOC_HD(w2, n_cells_ext, cs_real_6_t, cs_alloc_mode);
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      cs_math_sym_33_product(porosf[cell_id],
                             viscel[cell_id],
                             w2[cell_id]);
    });
    ctx.wait();
    viscce = w2;
  }

  /* ---> Periodicity and parallelism treatment of symmetric tensors */
  cs_halo_sync_r(halo, halo_type, false, viscce);

  /* 2. Compute the diffusive part with reconstruction technics */

  /* ======================================================================
     ---> Compute the gradient of the current variable if needed
     ======================================================================*/

  if (ircflp == 1) {

    cs_gradient_tensor_synced_input(var_name,
                                    gradient_type,
                                    halo_type,
                                    inc,
                                    eqp.nswrgr,
                                    eqp.verbosity,
                                    (cs_gradient_limit_t)(eqp.imligr),
                                    eqp.epsrgr,
                                    eqp.d_climgr,
                                    bc_coeffs_ts,
                                    _pvar,
                                    val_f,
                                    grad);

  }
  else {
    ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
      for (cs_lnum_t isou = 0; isou < 6; isou++) {
        for (cs_lnum_t jsou = 0; jsou < 3; jsou++)
          grad[cell_id][isou][jsou] = 0.;
      }
    });
  }

  /* ======================================================================
     ---> Contribution from interior faces
     ======================================================================*/

  /* Steady */
  if (idtvar < 0) {

    cs_host_context &h_ctx = static_cast<cs_host_context&>(ctx);
    i_sum_type = CS_DISPATCH_SUM_SIMPLE;

    h_ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t visci[3][3], viscj[3][3];
      cs_real_t diippf[3], djjppf[3], pipp[6], pjpp[6];
      cs_real_t  pir[6], pjr[6], pippr[6], pjppr[6];
      cs_real_t  pi[6], pj[6], pia[6], pja[6];
      cs_real_t fluxi[6], fluxj[6];

      for (cs_lnum_t isou = 0; isou < 6; isou++) {
        pi[isou] = _pvar[ii][isou];
        pj[isou] = _pvar[jj][isou];
        pia[isou] = pvara[ii][isou];
        pja[isou] = pvara[jj][isou];
      }

      cs_real_t bldfrp = (cs_real_t) ircflp;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflp > 0)
        bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                         0.);

      /* Recompute II" and JJ"
         ----------------------*/

      visci[0][0] = viscce[ii][0];
      visci[1][1] = viscce[ii][1];
      visci[2][2] = viscce[ii][2];
      visci[1][0] = viscce[ii][3];
      visci[0][1] = viscce[ii][3];
      visci[2][1] = viscce[ii][4];
      visci[1][2] = viscce[ii][4];
      visci[2][0] = viscce[ii][5];
      visci[0][2] = viscce[ii][5];

      /* IF.Ki.S / ||Ki.S||^2 */
      cs_real_t fikdvi_s = weighf[face_id][0] * i_face_surf[face_id];

      /* II" = IF + FI" */
      for (cs_lnum_t i = 0; i < 3; i++) {
        diippf[i] =   i_face_cog[face_id][i]-cell_cen[ii][i]
                    - fikdvi_s*(  visci[0][i]*i_face_u_normal[face_id][0]
                                + visci[1][i]*i_face_u_normal[face_id][1]
                                + visci[2][i]*i_face_u_normal[face_id][2]);
      }

      viscj[0][0] = viscce[jj][0];
      viscj[1][1] = viscce[jj][1];
      viscj[2][2] = viscce[jj][2];
      viscj[1][0] = viscce[jj][3];
      viscj[0][1] = viscce[jj][3];
      viscj[2][1] = viscce[jj][4];
      viscj[1][2] = viscce[jj][4];
      viscj[2][0] = viscce[jj][5];
      viscj[0][2] = viscce[jj][5];

      /* FJ.Kj.S / ||Kj.S||^2 */
      cs_real_t fjkdvi_s = weighf[face_id][1] * i_face_surf[face_id];

      /* JJ" = JF + FJ" */
      for (cs_lnum_t i = 0; i < 3; i++) {
        djjppf[i] =   i_face_cog[face_id][i]-cell_cen[jj][i]
                    + fjkdvi_s*(  viscj[0][i]*i_face_u_normal[face_id][0]
                                + viscj[1][i]*i_face_u_normal[face_id][1]
                                + viscj[2][i]*i_face_u_normal[face_id][2]);
      }

      for (cs_lnum_t isou = 0; isou < 6; isou++) {
        /* p in I" and J" */
        pipp[isou] = pi[isou] + bldfrp*(  grad[ii][isou][0]*diippf[0]
                                        + grad[ii][isou][1]*diippf[1]
                                        + grad[ii][isou][2]*diippf[2]);
        pjpp[isou] = pj[isou] + bldfrp*(  grad[jj][isou][0]*djjppf[0]
                                        + grad[jj][isou][1]*djjppf[1]
                                        + grad[jj][isou][2]*djjppf[2]);

        pir[isou] = pi[isou]/relaxp - (1.-relaxp)/relaxp * pia[isou];
        pjr[isou] = pj[isou]/relaxp - (1.-relaxp)/relaxp * pja[isou];

        /* pr in I" and J" */
        pippr[isou] = pir[isou] + bldfrp*(  grad[ii][isou][0]*diippf[0]
                                          + grad[ii][isou][1]*diippf[1]
                                          + grad[ii][isou][2]*diippf[2]);
        pjppr[isou] = pjr[isou] + bldfrp*(  grad[jj][isou][0]*djjppf[0]
                                          + grad[jj][isou][1]*djjppf[1]
                                          + grad[jj][isou][2]*djjppf[2]);

        fluxi[isou] = - i_visc[face_id]*(pippr[isou] - pjpp[isou]);
        fluxj[isou] =   i_visc[face_id]*(pipp[isou] - pjppr[isou]);

      } // loop on isou

      if (ii < n_cells)
        cs_dispatch_sum<6>(rhs[ii], fluxi, i_sum_type);
      if (jj < n_cells)
        cs_dispatch_sum<6>(rhs[jj], fluxj, i_sum_type);

    });

    /* Unsteady */
  }
  else {

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t visci[3][3], viscj[3][3];
      cs_real_t diippf[3], djjppf[3], pipp[6], pjpp[6];
      cs_real_t pi[6], pj[6];

      for (cs_lnum_t isou = 0; isou < 6; isou++) {
        pi[isou] = _pvar[ii][isou];
        pj[isou] = _pvar[jj][isou];
      }

      cs_real_t bldfrp = (cs_real_t) ircflp;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflp > 0)
        bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                         0.);

      /* Recompute II" and JJ"
         ----------------------*/

      visci[0][0] = viscce[ii][0];
      visci[1][1] = viscce[ii][1];
      visci[2][2] = viscce[ii][2];
      visci[1][0] = viscce[ii][3];
      visci[0][1] = viscce[ii][3];
      visci[2][1] = viscce[ii][4];
      visci[1][2] = viscce[ii][4];
      visci[2][0] = viscce[ii][5];
      visci[0][2] = viscce[ii][5];

      /* IF.Ki.S / ||Ki.S||^2 */
      cs_real_t fikdvi_s = weighf[face_id][0] * i_face_surf[face_id];

      /* II" = IF + FI" */
      for (cs_lnum_t i = 0; i < 3; i++) {
        diippf[i] =   i_face_cog[face_id][i]-cell_cen[ii][i]
                    - fikdvi_s*(  visci[0][i]*i_face_u_normal[face_id][0]
                                + visci[1][i]*i_face_u_normal[face_id][1]
                                + visci[2][i]*i_face_u_normal[face_id][2]);
      }

      viscj[0][0] = viscce[jj][0];
      viscj[1][1] = viscce[jj][1];
      viscj[2][2] = viscce[jj][2];
      viscj[1][0] = viscce[jj][3];
      viscj[0][1] = viscce[jj][3];
      viscj[2][1] = viscce[jj][4];
      viscj[1][2] = viscce[jj][4];
      viscj[2][0] = viscce[jj][5];
      viscj[0][2] = viscce[jj][5];

      /* FJ.Kj.S / ||Kj.S||^2 */
      cs_real_t fjkdvi_s = weighf[face_id][1] * i_face_surf[face_id];

      /* JJ" = JF + FJ" */
      for (cs_lnum_t i = 0; i < 3; i++) {
        djjppf[i] = i_face_cog[face_id][i]-cell_cen[jj][i]
                    + fjkdvi_s*(  viscj[0][i]*i_face_u_normal[face_id][0]
                                + viscj[1][i]*i_face_u_normal[face_id][1]
                                + viscj[2][i]*i_face_u_normal[face_id][2]);
      }

      cs_real_t fluxi[6], fluxj[6];

      for (cs_lnum_t isou = 0; isou < 6; isou++) {
        /* p in I" and J" */
        pipp[isou] =  pi[isou]
                    + bldfrp * cs_math_3_dot_product(grad[ii][isou], diippf);

        pjpp[isou] =  pj[isou]
                    + bldfrp * cs_math_3_dot_product(grad[jj][isou], djjppf);

        cs_real_t flux = i_visc[face_id] * (pipp[isou] - pjpp[isou]);

        fluxj[isou] = thetap * flux;
        fluxi[isou] = - fluxj[isou];
      }

      if (ii < n_cells)
        cs_dispatch_sum<6>(rhs[ii], fluxi, i_sum_type);
      if (jj < n_cells)
        cs_dispatch_sum<6>(rhs[jj], fluxj, i_sum_type);
    });
  }

  /* ======================================================================
     ---> Contribution from boundary faces
     ======================================================================*/

  /* Steady */
  if (idtvar < 0) {

    cs_host_context &h_ctx = static_cast<cs_host_context&>(ctx);
    b_sum_type = CS_DISPATCH_SUM_SIMPLE;

    h_ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = b_face_cells[face_id];

      cs_real_t pi[6], pia[6], pir[6], pippr[6];
      cs_real_t visci[3][3];
      cs_real_t diippf[3];

      for (cs_lnum_t isou = 0; isou < 6; isou++) {
        pi[isou] = _pvar[ii][isou];
        pia[isou] = pvara[ii][isou];
        pir[isou] = pi[isou]/relaxp - (1.-relaxp)/relaxp*pia[isou];
      }

      cs_real_t bldfrp = (cs_real_t) ircflb;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflb > 0)
        bldfrp = cs::max(df_limiter[ii], 0.);

      /* Recompute II"
         --------------*/

      visci[0][0] = viscce[ii][0];
      visci[1][1] = viscce[ii][1];
      visci[2][2] = viscce[ii][2];
      visci[1][0] = viscce[ii][3];
      visci[0][1] = viscce[ii][3];
      visci[2][1] = viscce[ii][4];
      visci[1][2] = viscce[ii][4];
      visci[2][0] = viscce[ii][5];
      visci[0][2] = viscce[ii][5];

      /* IF.Ki.S / ||Ki.S||^2 */
      cs_real_t fikdvi_s = weighb[face_id] * b_face_surf[face_id];

      /* II" = IF + FI" */
      for (cs_lnum_t i = 0; i < 3; i++) {
        diippf[i] =   b_face_cog[face_id][i] - cell_cen[ii][i]
                    - fikdvi_s*(  visci[0][i]*b_face_u_normal[face_id][0]
                                + visci[1][i]*b_face_u_normal[face_id][1]
                                + visci[2][i]*b_face_u_normal[face_id][2]);
      }
      for (cs_lnum_t isou = 0; isou < 6; isou++) {
        pippr[isou] = pir[isou] + bldfrp*(  grad[ii][isou][0]*diippf[0]
                                          + grad[ii][isou][1]*diippf[1]
                                          + grad[ii][isou][2]*diippf[2]);
      }

      cs_real_t flux[6];

      for (cs_lnum_t isou = 0; isou < 6; isou++) {
        cs_real_t pfacd = inc*cofaf[face_id][isou];
        for (cs_lnum_t jsou = 0; jsou < 6; jsou++) {
          pfacd += cofbf[face_id][isou][jsou]*pippr[jsou];
        }

        flux[isou] = - b_visc[face_id]*pfacd;
      }

      cs_dispatch_sum<6>(rhs[ii], flux, b_sum_type);
    });

  }

  /* Unsteady */

  else {

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = b_face_cells[face_id];

      cs_real_t visci[3][3];
      cs_real_t diippf[3], pi[6], pipp[6];

      for (cs_lnum_t isou = 0; isou < 6; isou++) {
        pi[isou] = _pvar[ii][isou];
      }

      cs_real_t bldfrp = (cs_real_t) ircflb;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflb > 0)
        bldfrp = cs::max(df_limiter[ii], 0.);

      /* Recompute II"
         --------------*/

      visci[0][0] = viscce[ii][0];
      visci[1][1] = viscce[ii][1];
      visci[2][2] = viscce[ii][2];
      visci[1][0] = viscce[ii][3];
      visci[0][1] = viscce[ii][3];
      visci[2][1] = viscce[ii][4];
      visci[1][2] = viscce[ii][4];
      visci[2][0] = viscce[ii][5];
      visci[0][2] = viscce[ii][5];

      /* IF.Ki.S / ||Ki.S||^2 */
      cs_real_t fikdvi_s = weighb[face_id] * b_face_surf[face_id];

      /* II" = IF + FI" */
      for (cs_lnum_t i = 0; i < 3; i++) {
        diippf[i] = b_face_cog[face_id][i] - cell_cen[ii][i]
                    - fikdvi_s*(  visci[0][i]*b_face_u_normal[face_id][0]
                                + visci[1][i]*b_face_u_normal[face_id][1]
                                + visci[2][i]*b_face_u_normal[face_id][2]);
      }
      for (cs_lnum_t isou = 0; isou < 6; isou++) {
        pipp[isou] =  pi[isou]
                    + bldfrp * cs_math_3_dot_product(grad[ii][isou], diippf);
      }

      cs_real_t fluxi[6];

      for (cs_lnum_t isou = 0; isou < 6; isou++) {
        cs_real_t pfacd = inc*cofaf[face_id][isou];
        for (cs_lnum_t jsou = 0; jsou < 6; jsou++) {
          pfacd += cofbf[face_id][isou][jsou]*pipp[jsou];
        }

        cs_real_t flux = b_visc[face_id]*pfacd;
        fluxi[isou] = -thetap * flux;
      }

      cs_dispatch_sum<6>(rhs[ii], fluxi, b_sum_type);
    });
  }

  ctx.wait();

  /* Free memory */
  CS_FREE(grad);
  CS_FREE(w2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the face mass flux with the face pressure (or pressure
 * increment, or pressure double increment) gradient.
 *
 * \f[
 * \dot{m}_\ij = \dot{m}_\ij
 *             - \Delta t \grad_\fij \delta p \cdot \vect{S}_\ij
 * \f]
 *
 * Please refer to the
 * <a href="../../theory.pdf#cs_face_diffusion_potential">
     <b>cs_face_diffusion_potential/cs_diffusion_potential</b></a>
 * section of the theory guide for more information.
 *
 * \param[in]     f_id          field id (or -1)
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
 * \param[in]     init          indicator
 *                               - 1 initialize the mass flux to 0
 *                               - 0 otherwise
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imrgra        indicator
 *                               - 0 iterative gradient
 *                               - 1 least squares gradient
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     imligp        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thank to neighbooring gradients
 *                               - = 1 thank to the mean gradient
 * \param[in]     iphydp        hydrostatic pressure indicator
 * \param[in]     iwgrp         indicator
 *                               - 1 weight gradient by vicosity*porosity
 *                               - weighting determined by field options
 * \param[in]     verbosity     verbosity
 * \param[in]     epsrgp        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgp        clipping coeffecient for the computation of
 *                               the gradient
 * \param[in]     frcxt         body force creating the hydrostatic pressure
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     bc_coeffs     boundary condition structure for the variable
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     visel         viscosity by cell
 * \param[in,out] i_massflux    mass flux at interior faces
 * \param[in,out] b_massflux    mass flux at boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_face_diffusion_potential(const int                   f_id,
                            const cs_mesh_t            *m,
                            cs_mesh_quantities_t       *fvq,
                            int                         init,
                            int                         inc,
                            int                         imrgra,
                            int                         nswrgp,
                            int                         imligp,
                            int                         iphydp,
                            int                         iwgrp,
                            int                         verbosity,
                            double                      epsrgp,
                            double                      climgp,
                            cs_real_3_t       *restrict frcxt,
                            cs_real_t         *restrict pvar,
                            const cs_field_bc_coeffs_t *bc_coeffs,
                            const cs_real_t             i_visc[],
                            const cs_real_t             b_visc[],
                            cs_real_t         *restrict visel,
                            cs_real_t         *restrict i_massflux,
                            cs_real_t         *restrict b_massflux)
{
  CS_PROFILE_FUNC_RANGE();

  cs_real_t *cofafp = bc_coeffs->af;
  cs_real_t *cofbfp = bc_coeffs->bf;

  const cs_halo_t  *halo = m->halo;

  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict i_f_face_surf = fvq->i_face_surf;
  const cs_rreal_3_t *restrict diipf = fvq->diipf;
  const cs_rreal_3_t *restrict djjpf = fvq->djjpf;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;

  /* Parallel or device dispatch */
  cs_dispatch_context ctx_i, ctx_b;
#if defined(HAVE_CUDA)
  if (ctx_b.use_gpu())
    ctx_b.set_cuda_stream(cs_cuda_get_stream(1));
#endif

  const bool on_device = ctx_i.use_gpu();

  /* Local variables */

  int w_stride = 1;
  cs_real_t *gweight = nullptr;

  /*===========================================================================*/

  /* i_visc and carry similar information */

  /*===========================================================================
    1. Initialization
    ==========================================================================*/

  if (init >= 1) {
    ctx_i.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      i_massflux[face_id] = 0.;
    });

    ctx_b.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      b_massflux[face_id] = 0.;
    });
  }
  else if (init != 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of init"));
  }

  /* Use iterative gradient */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  if (imrgra < 0)
    imrgra = 0;

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  cs_field_t *f = nullptr;
  char var_name[64];

  if (f_id > -1) {
    f = cs_field_by_id(f_id);
    snprintf(var_name, 63, "%s", f->name);
  }
  else
    strncpy(var_name, "[face mass flux update]", 63);
  var_name[63] = '\0';

  /* Handle parallelism and periodicity */

  if (halo != nullptr)
    cs_halo_sync(halo, halo_type, on_device, pvar);

  /*==========================================================================
    2. Update mass flux without reconstruction
    ==========================================================================*/

  if (nswrgp <= 1) {

    /* Mass flow through interior faces */

    ctx_i.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      i_massflux[face_id] += i_visc[face_id]*(pvar[ii] - pvar[jj]);
    });

    /* Mass flow through boundary faces */

    ctx_b.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      cs_lnum_t ii = b_face_cells[face_id];
      cs_real_t pfac = inc*cofafp[face_id] + cofbfp[face_id]*pvar[ii];

      b_massflux[face_id] += b_visc[face_id]*pfac;
    });

    ctx_i.wait();
    ctx_b.wait();

  }

  /*==========================================================================
    3. Update mass flux with reconstruction technique if the mesh is non
       orthogonal
    ==========================================================================*/

  if (nswrgp > 1) {

    /* Allocate a work array for the gradient calculation */
    cs_real_3_t *grad;
    CS_MALLOC_HD(grad, n_cells_ext, cs_real_3_t, cs_alloc_mode);

    /* Compute gradient */
    if (iwgrp > 0) {
      gweight = visel;
      if (halo != nullptr)
        cs_halo_sync(halo, halo_type, on_device, gweight);
    }

    else if (f_id > -1) {
      /* Get the calculation option from the field */
      const cs_equation_param_t *eqp
        = cs_field_get_equation_param_const(f);
      if (f->type & CS_FIELD_VARIABLE && eqp->iwgrec == 1) {
        if (eqp->idiff > 0) {
          int key_id = cs_field_key_id("gradient_weighting_id");
          int diff_id = cs_field_get_key_int(f, key_id);
          if (diff_id > -1) {
            cs_field_t *weight_f = cs_field_by_id(diff_id);
            gweight = weight_f->val;
            w_stride = weight_f->dim;
            cs_field_synchronize(weight_f, halo_type);
          }
        }
      }
    }

    cs_gradient_scalar_synced_input(var_name,
                                    gradient_type,
                                    halo_type,
                                    inc,
                                    nswrgp,
                                    iphydp,
                                    w_stride,
                                    verbosity,
                                    (cs_gradient_limit_t)imligp,
                                    epsrgp,
                                    climgp,
                                    frcxt,
                                    bc_coeffs,
                                    (const cs_real_t *)pvar,
                                    gweight, /* Weighted gradient */
                                    nullptr, /* internal coupling */
                                    grad);

    /* Handle parallelism and periodicity */

    if (halo != nullptr)
      cs_halo_sync(halo, halo_type, on_device, visel);

    /* Mass flow through interior faces */

    ctx_i.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t dpxf = 0.5*(  visel[ii]*grad[ii][0]
                            + visel[jj]*grad[jj][0]);
      cs_real_t dpyf = 0.5*(  visel[ii]*grad[ii][1]
                            + visel[jj]*grad[jj][1]);
      cs_real_t dpzf = 0.5*(  visel[ii]*grad[ii][2]
                            + visel[jj]*grad[jj][2]);

      /*---> Dij = IJ - (IJ.N) N = II' - JJ' */
      cs_real_t dijx = diipf[face_id][0] - djjpf[face_id][0];
      cs_real_t dijy = diipf[face_id][1] - djjpf[face_id][1];
      cs_real_t dijz = diipf[face_id][2] - djjpf[face_id][2];

      cs_real_t f_mass_flux =   i_visc[face_id]*(pvar[ii] - pvar[jj])
                              + (  dpxf *dijx + dpyf*dijy + dpzf*dijz)
                                 * i_f_face_surf[face_id]/i_dist[face_id];

      i_massflux[face_id] += f_mass_flux;
    });

    /* Mass flow through boundary faces */

    ctx_b.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      cs_lnum_t ii = b_face_cells[face_id];

      cs_real_t pip = pvar[ii] + cs_math_3_dot_product(grad[ii], diipb[face_id]);
      cs_real_t pfac = inc*cofafp[face_id] + cofbfp[face_id]*pip;

      b_massflux[face_id] += b_visc[face_id]*pfac;
    });

    ctx_i.wait();
    ctx_b.wait();

    /* Free memory */
    CS_FREE_HD(grad);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the pressure gradient term to the mass flux
 * in case of anisotropic diffusion of the pressure field \f$ P \f$.
 *
 * More precisely, the mass flux side \f$ \dot{m}_\fij \f$ is updated as
 * follows:
 * \f[
 * \dot{m}_\fij = \dot{m}_\fij -
 *              \left( \tens{\mu}_\fij \gradv_\fij P \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * \param[in]     f_id          field id (or -1)
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
 * \param[in]     init           indicator
 *                               - 1 initialize the mass flux to 0
 *                               - 0 otherwise
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imrgra        indicator
 *                               - 0 iterative gradient
 *                               - 1 least squares gradient
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     imligp        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thank to neighbooring gradients
 *                               - = 1 thank to the mean gradient
 * \param[in]     ircflp        indicator
 *                               - 1 flux reconstruction,
 *                               - 0 otherwise
 * \param[in]     ircflb        indicator
 *                               - 1 boundary flux reconstruction,
 *                               - 0 otherwise
 * \param[in]     iphydp        indicator
 *                               - 1 hydrostatic pressure taken into account
 *                               - 0 otherwise
 * \param[in]     iwgrp         indicator
 *                               - 1 weight gradient by vicosity*porosity
 *                               - weighting determined by field options
 * \param[in]     iwarnp        verbosity
 * \param[in]     epsrgp        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgp        clipping coeffecient for the computation of
 *                               the gradient
 * \param[in]     frcxt         body force creating the hydrostatic pressure
 * \param[in]     pvar          solved variable (pressure)
 * \param[in]     bc_coeffs     boundary condition structure for the variable
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]     weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in,out] i_massflux    mass flux at interior faces
 * \param[in,out] b_massflux    mass flux at boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_face_anisotropic_diffusion_potential(const int                   f_id,
                                        const cs_mesh_t            *m,
                                        cs_mesh_quantities_t       *fvq,
                                        int                         init,
                                        int                         inc,
                                        int                         imrgra,
                                        int                         nswrgp,
                                        int                         imligp,
                                        int                         ircflp,
                                        int                         ircflb,
                                        int                         iphydp,
                                        int                         iwgrp,
                                        int                         iwarnp,
                                        double                      epsrgp,
                                        double                      climgp,
                                        cs_real_3_t       *restrict frcxt,
                                        cs_real_t         *restrict pvar,
                                        const cs_field_bc_coeffs_t *bc_coeffs,
                                        const cs_real_t             i_visc[],
                                        const cs_real_t             b_visc[],
                                        cs_real_6_t       *restrict viscel,
                                        const cs_real_2_t           weighf[],
                                        const cs_real_t             weighb[],
                                        cs_real_t         *restrict i_massflux,
                                        cs_real_t         *restrict b_massflux)
{
  cs_real_t *cofafp = bc_coeffs->af;
  cs_real_t *cofbfp = bc_coeffs->bf;

  const cs_halo_t  *halo = m->halo;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_real_t *restrict i_face_surf = fvq->i_face_surf;
  const cs_real_t *restrict b_face_surf = fvq->b_face_surf;
  const cs_nreal_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;
  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog = fvq->b_face_cog;

  /* Local variables */

  cs_real_t *df_limiter = nullptr;

  char var_name[64];
  int w_stride = 6;

  cs_field_t *f = nullptr;

  cs_real_t *gweight = nullptr;

  /* Parallel or device dispatch */
  cs_dispatch_context ctx, ctx_b;
#if defined(HAVE_CUDA)
  if (ctx_b.use_gpu())
    ctx_b.set_cuda_stream(cs_cuda_get_stream(1));
#endif
  bool on_device = ctx.use_gpu();

  /*==========================================================================
    1. Initialization
    ==========================================================================*/

  if (init >= 1) {
    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      i_massflux[face_id] = 0.;
    });
    ctx_b.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      b_massflux[face_id] = 0.;
    });
  }
  else if (init != 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of init"));
  }

  /* Use iterative gradient */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  if (imrgra < 0)
    imrgra = 0;

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  if (f_id > -1) {
    f = cs_field_by_id(f_id);

    int df_limiter_id =
      cs_field_get_key_int(f, cs_field_key_id("diffusion_limiter_id"));
    if (df_limiter_id > -1)
      df_limiter = cs_field_by_id(df_limiter_id)->val;

    snprintf(var_name, 63, "%s", f->name);
  }
  else
    strncpy(var_name, "[face mass flux update]", 63);
  var_name[63] = '\0';

  /* Porosity fields */
  cs_field_t *fporo = cs_field_by_name_try("porosity");
  cs_field_t *ftporo = cs_field_by_name_try("tensorial_porosity");

  cs_real_t *porosi = nullptr;
  cs_real_6_t *porosf = nullptr;

  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2) {
    porosi = fporo->val;
    if (ftporo != nullptr) {
      porosf = (cs_real_6_t *)ftporo->val;
    }
  }

  /* Handle parallelism and periodicity */

  cs_halo_sync(halo, halo_type, on_device, pvar);

  /*==========================================================================
    2. Update mass flux without reconstruction technique
    ==========================================================================*/

  if (nswrgp <= 1) {

    /* Contribution from interior faces */

    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      i_massflux[face_id] += i_visc[face_id]*(pvar[ii] - pvar[jj]);
    });

    /* Contribution from boundary faces */

    ctx_b.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
      cs_lnum_t ii = b_face_cells[face_id];
      double pfac = inc*cofafp[face_id] + cofbfp[face_id]*pvar[ii];

      b_massflux[face_id] += b_visc[face_id]*pfac;
    });

    ctx.wait();
    ctx_b.wait();

  }

  /*==========================================================================
    3. Update mass flux with reconstruction technique
    ==========================================================================*/

  if (nswrgp > 1) {

    cs_real_6_t *viscce = nullptr;
    cs_real_6_t *w2 = nullptr;

    /* Without porosity */
    if (porosi == nullptr) {
      viscce = viscel;

      /* With porosity */
    }
    else if (porosi != nullptr && porosf == nullptr) {
      CS_MALLOC_HD(w2, n_cells_ext, cs_real_6_t, cs_alloc_mode);
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t  cell_id) {
        for (cs_lnum_t isou = 0; isou < 6; isou++) {
          w2[cell_id][isou] = porosi[cell_id]*viscel[cell_id][isou];
        }
      });
      viscce = w2;

      /* With tensorial porosity */
    }
    else if (porosi != nullptr && porosf != nullptr) {
      CS_MALLOC_HD(w2, n_cells_ext, cs_real_6_t, cs_alloc_mode);
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t  cell_id) {
        cs_math_sym_33_product(porosf[cell_id],
                               viscel[cell_id],
                               w2[cell_id]);
      });
      viscce = w2;
    }

    ctx.wait();

    /* Periodicity and parallelism treatment of symmetric tensors */
    cs_halo_sync_r(halo, CS_HALO_STANDARD, on_device, viscce);

    /* Allocate a work array for the gradient calculation */
    cs_real_3_t *grad;
    CS_MALLOC_HD(grad, n_cells_ext, cs_real_3_t, cs_alloc_mode);

    /* Compute gradient */
    if (iwgrp > 0) {
      gweight = (cs_real_t *)viscce;
    }

    else if (f_id > -1) {
      /* Get the calculation option from the field */
      const cs_equation_param_t *eqp
        = cs_field_get_equation_param_const(f);
      if (f->type & CS_FIELD_VARIABLE && eqp->iwgrec == 1) {
        if (eqp->idifft > 0) {
          int key_id = cs_field_key_id("gradient_weighting_id");
          int diff_id = cs_field_get_key_int(f, key_id);
          if (diff_id > -1) {
            cs_field_t *weight_f = cs_field_by_id(diff_id);
            gweight = weight_f->val;
            w_stride = weight_f->dim;
            cs_field_synchronize(weight_f, halo_type, on_device);
          }
        }
      }
    }

    /* Compute gradient */
    cs_gradient_scalar_synced_input(var_name,
                                    gradient_type,
                                    halo_type,
                                    inc,
                                    nswrgp,
                                    iphydp,
                                    w_stride,
                                    iwarnp,
                                    (cs_gradient_limit_t)imligp,
                                    epsrgp,
                                    climgp,
                                    frcxt,
                                    bc_coeffs,
                                    pvar,
                                    gweight, /* Weighted gradient */
                                    nullptr, /* internal coupling */
                                    grad);

    /* Mass flow through interior faces */

    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t pi = pvar[ii];
      cs_real_t pj = pvar[jj];

      /* Recompute II" and JJ"
         -------------------- */

      cs_real_t visci[3][3], viscj[3][3];
      cs_real_t diippf[3], djjppf[3];

      visci[0][0] = viscce[ii][0];
      visci[1][1] = viscce[ii][1];
      visci[2][2] = viscce[ii][2];
      visci[1][0] = viscce[ii][3];
      visci[0][1] = viscce[ii][3];
      visci[2][1] = viscce[ii][4];
      visci[1][2] = viscce[ii][4];
      visci[2][0] = viscce[ii][5];
      visci[0][2] = viscce[ii][5];

      cs_real_t bldfrp = (cs_real_t) ircflp;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflp > 0)
        bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]), 0.);

      /* IF.Ki.S / ||Ki.S||^2 */
      cs_real_t fikdvi_s = weighf[face_id][0] * i_face_surf[face_id];

      /* II" = IF + FI" */
      for (cs_lnum_t i = 0; i < 3; i++) {
        diippf[i] =   i_face_cog[face_id][i]-cell_cen[ii][i]
                    - fikdvi_s * (  visci[0][i]*i_face_u_normal[face_id][0]
                                  + visci[1][i]*i_face_u_normal[face_id][1]
                                  + visci[2][i]*i_face_u_normal[face_id][2]);
      }

      viscj[0][0] = viscce[jj][0];
      viscj[1][1] = viscce[jj][1];
      viscj[2][2] = viscce[jj][2];
      viscj[1][0] = viscce[jj][3];
      viscj[0][1] = viscce[jj][3];
      viscj[2][1] = viscce[jj][4];
      viscj[1][2] = viscce[jj][4];
      viscj[2][0] = viscce[jj][5];
      viscj[0][2] = viscce[jj][5];

      /* FJ.Kj.S / ||Kj.S||^2 */
      cs_real_t fjkdvi_s = weighf[face_id][1] * i_face_surf[face_id];

      /* JJ" = JF + FJ" */
      for (cs_lnum_t i = 0; i < 3; i++) {
        djjppf[i] =   i_face_cog[face_id][i]-cell_cen[jj][i]
                    + fjkdvi_s*(  viscj[0][i]*i_face_u_normal[face_id][0]
                                + viscj[1][i]*i_face_u_normal[face_id][1]
                                + viscj[2][i]*i_face_u_normal[face_id][2]);
      }

      /* p in I" and J" */
      cs_real_t pipp = pi + bldfrp * cs_math_3_dot_product(grad[ii], diippf);
      cs_real_t pjpp = pj + bldfrp * cs_math_3_dot_product(grad[jj], djjppf);

      i_massflux[face_id] += i_visc[face_id]*(pipp - pjpp);
    });

    /* Contribution from boundary faces */

    ctx_b.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
      cs_lnum_t ii = b_face_cells[face_id];

      cs_real_t pi = pvar[ii];

      cs_real_t bldfrp = (cs_real_t) ircflb;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflb > 0)
        bldfrp = cs::max(df_limiter[ii], 0.);

      /* Recompute II"
         ------------- */

      cs_real_t visci[3][3];
      cs_real_t diippf[3];

      visci[0][0] = viscce[ii][0];
      visci[1][1] = viscce[ii][1];
      visci[2][2] = viscce[ii][2];
      visci[1][0] = viscce[ii][3];
      visci[0][1] = viscce[ii][3];
      visci[2][1] = viscce[ii][4];
      visci[1][2] = viscce[ii][4];
      visci[2][0] = viscce[ii][5];
      visci[0][2] = viscce[ii][5];

      /* IF.Ki.S / ||Ki.S||^2 */
      cs_real_t fikdvi_s = weighb[face_id] * b_face_surf[face_id];

      /* II" = IF + FI" */
      for (cs_lnum_t i = 0; i < 3; i++) {
        diippf[i] =   b_face_cog[face_id][i] - cell_cen[ii][i]
                    - fikdvi_s*(  visci[0][i]*b_face_u_normal[face_id][0]
                                + visci[1][i]*b_face_u_normal[face_id][1]
                                + visci[2][i]*b_face_u_normal[face_id][2]);
      }

      cs_real_t pipp = pi + bldfrp * cs_math_3_dot_product(grad[ii], diippf);

      cs_real_t pfac = inc*cofafp[face_id] + cofbfp[face_id]*pipp;

      b_massflux[face_id] += b_visc[face_id]*pfac;
    });

    ctx.wait();
    ctx_b.wait();

    /* Free memory */
    CS_FREE(grad);
    CS_FREE(w2);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the cell mass flux divergence with the face pressure (or
 * pressure increment, or pressure double increment) gradient.
 *
 * \f[
 * \dot{m}_\ij = \dot{m}_\ij
 *             - \sum_j \Delta t \grad_\fij p \cdot \vect{S}_\ij
 * \f]
 *
 * \param[in]     f_id          field id (or -1)
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
 * \param[in]     init          indicator
 *                               - 1 initialize the mass flux to 0
 *                               - 0 otherwise
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imrgra        indicator
 *                               - 0 iterative gradient
 *                               - 1 least squares gradient
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     imligp        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thank to neighbooring gradients
 *                               - = 1 thank to the mean gradient
 * \param[in]     iphydp        hydrostatic pressure indicator
 * \param[in]     iwarnp        verbosity
 * \param[in]     iwgrp         indicator
 *                               - 1 weight gradient by vicosity*porosity
 *                               - weighting determined by field options
 * \param[in]     epsrgp        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgp        clipping coeffecient for the computation of
 *                               the gradient
 * \param[in]     frcxt         body force creating the hydrostatic pressure
 * \param[in]     pvar          solved variable (current time step)
 * \param[in]     bc_coeffs     boundary condition structure for the variable
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     visel         viscosity by cell
 * \param[in,out] diverg        mass flux divergence
 */
/*----------------------------------------------------------------------------*/

void
cs_diffusion_potential(const int                   f_id,
                       const cs_mesh_t            *m,
                       cs_mesh_quantities_t       *fvq,
                       int                         init,
                       int                         inc,
                       int                         imrgra,
                       int                         nswrgp,
                       int                         imligp,
                       int                         iphydp,
                       int                         iwgrp,
                       int                         iwarnp,
                       double                      epsrgp,
                       double                      climgp,
                       cs_real_3_t       *restrict frcxt,
                       cs_real_t         *restrict pvar,
                       const cs_field_bc_coeffs_t *bc_coeffs,
                       const cs_real_t             i_visc[],
                       const cs_real_t             b_visc[],
                       cs_real_t                   visel[],
                       cs_real_t         *restrict diverg)
{
  CS_PROFILE_FUNC_RANGE();

  cs_real_t *cofafp = bc_coeffs->af;
  cs_real_t *cofbfp = bc_coeffs->bf;

  const cs_halo_t  *halo = m->halo;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict i_f_face_surf = fvq->i_face_surf;
  const cs_rreal_3_t *restrict diipf = fvq->diipf;
  const cs_rreal_3_t *restrict djjpf = fvq->djjpf;
  const cs_rreal_3_t *restrict diipb = fvq->diipb;

  /* Parallel or device dispatch */
  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);
  bool on_device = ctx.use_gpu();

  /* Local variables */

  char var_name[64];
  int mass_flux_rec_type = cs_glob_velocity_pressure_param->irecmf;
  int w_stride = 1;

  cs_field_t *f = nullptr;
  cs_real_3_t *grad = nullptr;
  cs_real_t *gweight = nullptr;

  /*==========================================================================
    1. Initialization
    ==========================================================================*/

  if (init >= 1) {
    ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
      diverg[ii] = 0.;
    });
  }
  else if (init < 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of init"));
  }

  /* Use iterative gradient */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  if (imrgra < 0)
    imrgra = 0;

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);
  if (f_id != -1) {
    f = cs_field_by_id(f_id);
    snprintf(var_name, 63, "%s", f->name);
  }
  else
    strncpy(var_name, "[cell mass flux divergence update]", 63);
  var_name[63] = '\0';

  /* Handle parallelism and periodicity */

  if (halo != nullptr)
    cs_halo_sync(halo, halo_type, ctx.use_gpu(), pvar);

  /*==========================================================================
    2. Update mass flux without reconstruction techniques
    ==========================================================================*/

  if (nswrgp <= 1) {

    /* Mass flow through interior faces */

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t i_massflux = i_visc[face_id]*(pvar[ii] - pvar[jj]);

      if (ii < n_cells)
        cs_dispatch_sum(&diverg[ii],  i_massflux, i_sum_type);
      if (jj < n_cells)
        cs_dispatch_sum(&diverg[jj], -i_massflux, i_sum_type);

    });

    /* Mass flow through boundary faces */

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {

      cs_lnum_t ii = b_face_cells[face_id];
      cs_real_t pfac = inc*cofafp[face_id] +cofbfp[face_id]*pvar[ii];

      cs_real_t b_massflux = b_visc[face_id]*pfac;
      cs_dispatch_sum(&diverg[ii], b_massflux, b_sum_type);

    });

  }

  /*==========================================================================
    3. Update mass flux with reconstruction if the mesh is non-orthogonal
    ==========================================================================*/

  if (nswrgp > 1) {

    /* Allocate a work array for the gradient calculation */
    CS_MALLOC_HD(grad, n_cells_ext, cs_real_3_t, cs_alloc_mode);

    /* Compute gradient */
    if (iwgrp > 0) {
      gweight = visel;
      if (halo != nullptr)
        cs_halo_sync(halo, halo_type, ctx.use_gpu(), gweight);
    }

    else if (f_id > -1) {
      /* Get the calculation option from the field */
      const cs_equation_param_t *eqp
        = cs_field_get_equation_param_const(f);
      if (f->type & CS_FIELD_VARIABLE && eqp->iwgrec == 1) {
        if (eqp->idiff > 0) {
          int key_id = cs_field_key_id("gradient_weighting_id");
          int diff_id = cs_field_get_key_int(f, key_id);
          if (diff_id > -1) {
            cs_field_t *weight_f = cs_field_by_id(diff_id);
            gweight = weight_f->val;
            w_stride = weight_f->dim;
            cs_field_synchronize(weight_f, halo_type);
          }
        }
      }
    }

    cs_real_t *_pvar = pvar;

    if (cs_glob_mesh_quantities_flag & CS_BAD_CELLS_REGULARISATION) {
      CS_MALLOC(_pvar, n_cells_ext, cs_real_t);
      cs_array_real_copy(n_cells_ext, pvar, _pvar);

      cs_bad_cells_regularisation_scalar(_pvar);
    }

    cs_gradient_scalar_synced_input(var_name,
                                    gradient_type,
                                    halo_type,
                                    inc,
                                    nswrgp,
                                    iphydp,
                                    w_stride,
                                    iwarnp,
                                    (cs_gradient_limit_t)imligp,
                                    epsrgp,
                                    climgp,
                                    frcxt,
                                    bc_coeffs,
                                    _pvar,
                                    gweight, /* Weighted gradient */
                                    nullptr, /* internal coupling */
                                    grad);

    if (_pvar != pvar)
      CS_FREE(_pvar);

    /* Handle parallelism and periodicity */

    cs_halo_sync(halo, halo_type, on_device, visel);

    /* Mass flow through interior faces */

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t i_massflux = i_visc[face_id]*(pvar[ii] - pvar[jj]);

      if (mass_flux_rec_type == 0) {
        /*---> Dij = IJ - (IJ.N) N = II' - JJ' */
        cs_real_t dijx = diipf[face_id][0] - djjpf[face_id][0];
        cs_real_t dijy = diipf[face_id][1] - djjpf[face_id][1];
        cs_real_t dijz = diipf[face_id][2] - djjpf[face_id][2];

        cs_real_t dpxf = 0.5*(  visel[ii]*grad[ii][0]
                              + visel[jj]*grad[jj][0]);
        cs_real_t dpyf = 0.5*(  visel[ii]*grad[ii][1]
                              + visel[jj]*grad[jj][1]);
        cs_real_t dpzf = 0.5*(  visel[ii]*grad[ii][2]
                              + visel[jj]*grad[jj][2]);

        i_massflux +=  (dpxf*dijx + dpyf*dijy + dpzf*dijz)
                      * i_f_face_surf[face_id]/i_dist[face_id];
      }
      else {
        i_massflux +=   i_visc[face_id]
                      * (  cs_math_3_dot_product(grad[ii], diipf[face_id])
                         - cs_math_3_dot_product(grad[jj], djjpf[face_id]));
      }

      if (ii < n_cells)
        cs_dispatch_sum(&diverg[ii],  i_massflux, i_sum_type);
      if (jj < n_cells)
        cs_dispatch_sum(&diverg[jj], -i_massflux, i_sum_type);

    });

    /* Mass flow through boundary faces */

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {

      cs_lnum_t ii = b_face_cells[face_id];

      cs_real_t pip = pvar[ii] + cs_math_3_dot_product(grad[ii],
                                                       diipb[face_id]);
      cs_real_t pfac = inc*cofafp[face_id] +cofbfp[face_id]*pip;
      cs_real_t b_massflux = b_visc[face_id]*pfac;

      cs_dispatch_sum(&diverg[ii], b_massflux, b_sum_type);

    });

  }
  ctx.wait();

  /* Free memory */
  CS_FREE_HD(grad);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the explicit part of the divergence of the mass flux due to the
 * pressure gradient (analog to cs_anisotropic_diffusion_scalar).
 *
 * More precisely, the divergence of the mass flux side
 * \f$ \sum_{\fij \in \Facei{\celli}} \dot{m}_\fij \f$ is updated as follows:
 * \f[
 * \sum_{\fij \in \Facei{\celli}} \dot{m}_\fij
 *  = \sum_{\fij \in \Facei{\celli}} \dot{m}_\fij
 *  - \sum_{\fij \in \Facei{\celli}}
 *    \left( \tens{\mu}_\fij \gradv_\fij P \cdot \vect{S}_\ij  \right)
 * \f]
 *
 * \param[in]     f_id          field id (or -1)
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
 * \param[in]     init           indicator
 *                               - 1 initialize the mass flux to 0
 *                               - 0 otherwise
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imrgra        indicator
 *                               - 0 iterative gradient
 *                               - 1 least squares gradient
 * \param[in]     nswrgp        number of reconstruction sweeps for the
 *                               gradients
 * \param[in]     imligp        clipping gradient method
 *                               - < 0 no clipping
 *                               - = 0 thank to neighbooring gradients
 *                               - = 1 thank to the mean gradient
 * \param[in]     ircflp        indicator
 *                               - 1 flux reconstruction,
 *                               - 0 otherwise
 * \param[in]     ircflb        indicator
 *                               - 1 boundary flux reconstruction,
 *                               - 0 otherwise
 * \param[in]     iphydp        indicator
 *                               - 1 hydrostatic pressure taken into account
 *                               - 0 otherwise
 * \param[in]     iwgrp         indicator
 *                               - 1 weight gradient by vicosity*porosity
 *                               - weighting determined by field options
 * \param[in]     iwarnp        verbosity
 * \param[in]     epsrgp        relative precision for the gradient
 *                               reconstruction
 * \param[in]     climgp        clipping coeffecient for the computation of
 *                               the gradient
 * \param[in]     frcxt         body force creating the hydrostatic pressure
 * \param[in]     pvar          solved variable (pressure)
 * \param[in]     bc_coeffs     boundary condition structure for the variable
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     viscel        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
 * \param[in]     weighf        internal face weight between cells i j in case
 *                               of tensor diffusion
 * \param[in]     weighb        boundary face weight for cells i in case
 *                               of tensor diffusion
 * \param[in,out] diverg        divergence of the mass flux
 */
/*----------------------------------------------------------------------------*/

void
cs_anisotropic_diffusion_potential(const int                   f_id,
                                   const cs_mesh_t            *m,
                                   cs_mesh_quantities_t       *fvq,
                                   int                         init,
                                   int                         inc,
                                   int                         imrgra,
                                   int                         nswrgp,
                                   int                         imligp,
                                   int                         ircflp,
                                   int                         ircflb,
                                   int                         iphydp,
                                   int                         iwgrp,
                                   int                         iwarnp,
                                   double                      epsrgp,
                                   double                      climgp,
                                   cs_real_3_t       *restrict frcxt,
                                   cs_real_t         *restrict pvar,
                                   const cs_field_bc_coeffs_t *bc_coeffs,
                                   const cs_real_t             i_visc[],
                                   const cs_real_t             b_visc[],
                                   cs_real_6_t       *restrict viscel,
                                   const cs_real_2_t           weighf[],
                                   const cs_real_t             weighb[],
                                   cs_real_t         *restrict diverg)
{
  cs_real_t *cofafp = bc_coeffs->af;
  cs_real_t *cofbfp = bc_coeffs->bf;

  const cs_halo_t  *halo = m->halo;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_3_t *restrict cell_cen = fvq->cell_cen;
  const cs_real_t *restrict i_face_surf = fvq->i_face_surf;
  const cs_real_t *restrict b_face_surf = fvq->b_face_surf;
  const cs_nreal_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;
  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
  const cs_real_3_t *restrict i_face_cog = fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog = fvq->b_face_cog;

  /* Local variables */

  cs_real_t *df_limiter = nullptr;

  char var_name[64];
  int w_stride = 6;

  cs_real_6_t *viscce = nullptr;
  cs_real_6_t *w2 = nullptr;
  cs_real_3_t *grad = nullptr;
  cs_field_t *f = nullptr;

  cs_real_t *gweight = nullptr;

  /* Parallel or device dispatch */
  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);
  bool on_device = ctx.use_gpu();

  /*==========================================================================
    1. Initialization
    ==========================================================================*/

  if (init >= 1) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
      diverg[ii] = 0.;
    });
  }
  else if (init != 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("invalid value of init"));
  }

  /* Use iterative gradient */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  if (imrgra < 0)
    imrgra = 0;

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  if (f_id != -1) {
    f = cs_field_by_id(f_id);

    int df_limiter_id
      = cs_field_get_key_int(f, cs_field_key_id("diffusion_limiter_id"));
    if (df_limiter_id > -1)
      df_limiter = cs_field_by_id(df_limiter_id)->val;

    snprintf(var_name, 63, "%s", f->name);
  }
  else
    strncpy(var_name, "[cell mass flux divergence update]", 63);
  var_name[63] = '\0';

  /* Porosity fields */
  cs_field_t *fporo = cs_field_by_name_try("porosity");
  cs_field_t *ftporo = cs_field_by_name_try("tensorial_porosity");

  cs_real_t *porosi = nullptr;
  cs_real_6_t *porosf = nullptr;

  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2) {
    porosi = fporo->val;
    if (ftporo != nullptr) {
      porosf = (cs_real_6_t *)ftporo->val;
    }
  }

  /* Handle parallelism and periodicity */

  cs_halo_sync(m->halo, halo_type, on_device, pvar);

  /*==========================================================================
    2. Update mass flux without reconstruction technics
    ==========================================================================*/

  if (nswrgp <= 1) {

    /* Mass flow through interior faces */

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      double flux = i_visc[face_id]*(pvar[ii] - pvar[jj]);

      if (ii < n_cells)
        cs_dispatch_sum(&diverg[ii], flux, i_sum_type);
      if (jj < n_cells)
        cs_dispatch_sum(&diverg[jj], -flux, i_sum_type);

    });

    /* Mass flow though boundary faces */

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = b_face_cells[face_id];
      double pfac = inc*cofafp[face_id] + cofbfp[face_id]*pvar[ii];

      double flux = b_visc[face_id]*pfac;
      cs_dispatch_sum(&diverg[ii], flux, b_sum_type);

    });

  }

  /*==========================================================================
    3. Update mass flux with reconstruction techniques
    ==========================================================================*/

  if (nswrgp > 1) {

    viscce = nullptr;
    w2 = nullptr;

    /* Without porosity */
    if (porosi == nullptr) {
      viscce = viscel;

      /* With porosity */
    }
    else if (porosi != nullptr && porosf == nullptr) {
      CS_MALLOC_HD(w2, n_cells_ext, cs_real_6_t, cs_alloc_mode);
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t  cell_id) {
        for (cs_lnum_t isou = 0; isou < 6; isou++) {
          w2[cell_id][isou] = porosi[cell_id]*viscel[cell_id][isou];
        }
      });
      viscce = w2;

      /* With tensorial porosity */
    }
    else if (porosi != nullptr && porosf != nullptr) {
      CS_MALLOC_HD(w2, n_cells_ext, cs_real_6_t, cs_alloc_mode);
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t  cell_id) {
        cs_math_sym_33_product(porosf[cell_id],
                               viscel[cell_id],
                               w2[cell_id]);
      });
      viscce = w2;
    }

    ctx.wait();

    /* Periodicity and parallelism treatment of symmetric tensors */
    cs_halo_sync_r(halo, CS_HALO_STANDARD, on_device, viscce);

    /* Allocate a work array for the gradient calculation */
    CS_MALLOC_HD(grad, n_cells_ext, cs_real_3_t, cs_alloc_mode);

    /* Compute gradient */
    if (iwgrp > 0) {
      gweight = (cs_real_t *)viscce;
    }

    else if (f_id > -1) {
      /* Get the calculation option from the field */
      const cs_equation_param_t *eqp
        = cs_field_get_equation_param_const(f);
      if (f->type & CS_FIELD_VARIABLE && eqp->iwgrec == 1) {
        if (eqp->idifft > 0) {
          int key_id = cs_field_key_id("gradient_weighting_id");
          int diff_id = cs_field_get_key_int(f, key_id);
          if (diff_id > -1) {
            cs_field_t *weight_f = cs_field_by_id(diff_id);
            gweight = weight_f->val;
            w_stride = weight_f->dim;
            cs_field_synchronize(weight_f, halo_type, on_device);
          }
        }
      }
    }

    /* Compute gradient */
    cs_gradient_scalar_synced_input(var_name,
                                    gradient_type,
                                    halo_type,
                                    inc,
                                    nswrgp,
                                    iphydp,
                                    w_stride,
                                    iwarnp,
                                    (cs_gradient_limit_t)imligp,
                                    epsrgp,
                                    climgp,
                                    frcxt,
                                    bc_coeffs,
                                    pvar,
                                    gweight, /* Weighted gradient */
                                    nullptr, /* internal coupling */
                                    grad);

    /* Mass flow through interior faces */

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t pi = pvar[ii];
      cs_real_t pj = pvar[jj];

      cs_real_t bldfrp = (cs_real_t) ircflp;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflp > 0)
        bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                         0.);

      /* Recompute II" and JJ"
         ----------------------*/

      cs_real_t visci[3][3], viscj[3][3];

      visci[0][0] = viscce[ii][0];
      visci[1][1] = viscce[ii][1];
      visci[2][2] = viscce[ii][2];
      visci[1][0] = viscce[ii][3];
      visci[0][1] = viscce[ii][3];
      visci[2][1] = viscce[ii][4];
      visci[1][2] = viscce[ii][4];
      visci[2][0] = viscce[ii][5];
      visci[0][2] = viscce[ii][5];

      /* IF.Ki.S / ||Ki.S||^2 */
      cs_real_t fikdvi = weighf[face_id][0] * i_face_surf[face_id];

      /* II" = IF + FI" */

      cs_real_t diippf[3], djjppf[3];

      for (cs_lnum_t i = 0; i < 3; i++) {
        diippf[i] =   i_face_cog[face_id][i]-cell_cen[ii][i]
                    - fikdvi*(  visci[0][i]*i_face_u_normal[face_id][0]
                              + visci[1][i]*i_face_u_normal[face_id][1]
                              + visci[2][i]*i_face_u_normal[face_id][2]);
      }

      viscj[0][0] = viscce[jj][0];
      viscj[1][1] = viscce[jj][1];
      viscj[2][2] = viscce[jj][2];
      viscj[1][0] = viscce[jj][3];
      viscj[0][1] = viscce[jj][3];
      viscj[2][1] = viscce[jj][4];
      viscj[1][2] = viscce[jj][4];
      viscj[2][0] = viscce[jj][5];
      viscj[0][2] = viscce[jj][5];

      /* FJ.Kj.S / ||Kj.S||^2 */
      cs_real_t fjkdvi = weighf[face_id][1] * i_face_surf[face_id];

      /* JJ" = JF + FJ" */
      for (cs_lnum_t i = 0; i < 3; i++) {
        djjppf[i] =   i_face_cog[face_id][i]-cell_cen[jj][i]
                    + fjkdvi*(  viscj[0][i]*i_face_u_normal[face_id][0]
                              + viscj[1][i]*i_face_u_normal[face_id][1]
                              + viscj[2][i]*i_face_u_normal[face_id][2]);
      }

      /* p in I" and J" */
      cs_real_t pipp = pi + bldfrp*(  grad[ii][0]*diippf[0]
                                    + grad[ii][1]*diippf[1]
                                    + grad[ii][2]*diippf[2]);
      cs_real_t pjpp = pj + bldfrp*(  grad[jj][0]*djjppf[0]
                                    + grad[jj][1]*djjppf[1]
                                    + grad[jj][2]*djjppf[2]);

      cs_real_t flux = i_visc[face_id]*(pipp - pjpp);

      if (ii < n_cells)
        cs_dispatch_sum(&diverg[ii], flux, i_sum_type);
      if (jj < n_cells)
        cs_dispatch_sum(&diverg[jj], -flux, i_sum_type);

    });

    /* Mass flow though boundary faces */

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = b_face_cells[face_id];

      cs_real_t pi = pvar[ii];

      cs_real_t bldfrp = (cs_real_t)ircflb;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflb > 0)
        bldfrp = cs::max(df_limiter[ii], 0.);

      /* Recompute II"
         --------------*/

      cs_real_t visci[3][3];

      visci[0][0] = viscce[ii][0];
      visci[1][1] = viscce[ii][1];
      visci[2][2] = viscce[ii][2];
      visci[1][0] = viscce[ii][3];
      visci[0][1] = viscce[ii][3];
      visci[2][1] = viscce[ii][4];
      visci[1][2] = viscce[ii][4];
      visci[2][0] = viscce[ii][5];
      visci[0][2] = viscce[ii][5];

      /* IF.Ki.S / ||Ki.S||^2 */
      cs_real_t fikdvi = weighb[face_id] * b_face_surf[face_id];

      cs_real_t diippf[3];

      /* II" = IF + FI" */
      for (cs_lnum_t i = 0; i < 3; i++) {
        diippf[i] =   b_face_cog[face_id][i] - cell_cen[ii][i]
                    - fikdvi*(  visci[0][i]*b_face_u_normal[face_id][0]
                              + visci[1][i]*b_face_u_normal[face_id][1]
                              + visci[2][i]*b_face_u_normal[face_id][2]);
      }

      cs_real_t pipp = pi + bldfrp*(  grad[ii][0]*diippf[0]
                                    + grad[ii][1]*diippf[1]
                                    + grad[ii][2]*diippf[2]);

      cs_real_t pfac = inc*cofafp[face_id] + cofbfp[face_id]*pipp;

      cs_real_t flux = b_visc[face_id]*pfac;
      cs_dispatch_sum(&diverg[ii], flux, b_sum_type);

    });

    /* Free memory */
    CS_FREE(grad);
    CS_FREE(w2);

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

/*----------------------------------------------------------------------------
 * Compute the local cell Courant number as the maximum of all cell face based
 * Courant number at each cell.
 *
 * parameters:
 *   f_id        <-- pointer to field
 *   ctx         <-- reference to dispatch context
 *   courant     --> cell Courant number
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_courant_number(const cs_field_t    *f,
                       cs_dispatch_context &ctx,
                       cs_real_t           *courant)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;

  cs_mesh_adjacencies_update_cell_i_faces();

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
  const cs_lnum_t *cell_i_faces = ma->cell_i_faces;
  const cs_lnum_t *cell_b_faces_idx = ma->cell_b_faces_idx;
  const cs_lnum_t *cell_b_faces = ma->cell_b_faces;

  const cs_real_t *restrict vol = fvq->cell_vol;

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  const cs_real_t *restrict i_massflux
    = cs_field_by_id( cs_field_get_key_int(f, kimasf) )->val;
  const cs_real_t *restrict b_massflux
    = cs_field_by_id( cs_field_get_key_int(f, kbmasf) )->val;

  const cs_real_t *restrict dt = CS_F_(dt)->val;

  const int *restrict c_disable_flag = (fvq->has_disable_flag) ?
    fvq->c_disable_flag : nullptr;

  /* ---> Contribution from interior and boundary faces */

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    /* Initialization */
    courant[c_id] = 0.;

    cs_real_t dvol = cs_mq_cell_vol_inv(c_id, c_disable_flag, vol);

    /* Loop on interior faces */
    const cs_lnum_t s_id_i = c2c_idx[c_id];
    const cs_lnum_t e_id_i = c2c_idx[c_id + 1];

    for (cs_lnum_t cidx = s_id_i; cidx < e_id_i; cidx++) {
      const cs_lnum_t face_id = cell_i_faces[cidx];
      cs_real_t cnt = cs::abs(i_massflux[face_id])*dt[c_id]*dvol;
      courant[c_id] = cs::max(courant[c_id], cnt); //FIXME may contain rho
    }

    /* Loop on boundary faces */
    const cs_lnum_t s_id_b = cell_b_faces_idx[c_id];
    const cs_lnum_t e_id_b = cell_b_faces_idx[c_id + 1];

    for (cs_lnum_t cidx = s_id_b; cidx < e_id_b; cidx++) {
      const cs_lnum_t face_id = cell_b_faces[cidx];
      cs_real_t cnt = cs::abs(b_massflux[face_id])*dt[c_id]*dvol;
      courant[c_id] = cs::max(courant[c_id], cnt);
    }
  });

  ctx.wait();

  if (m->halo != nullptr)
    cs_halo_sync(m->halo, CS_HALO_STANDARD, ctx.use_gpu(), courant);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the upwind gradient used in the slope tests.
 *
 * This function assumes the input gradient and pvar values have already
 * been synchronized.
 *
 * \param[in]     f_id         field id
 * \param[in]     ctx          Reference to dispatch context
 * \param[in]     inc          Not an increment flag
 * \param[in]     grad         standard gradient
 * \param[out]    grdpa        upwind gradient
 * \param[in]     pvar         values
 * \param[in]     bc_coeffs    boundary condition structure for the variable
 * \param[in]     i_massflux   mass flux at interior faces
 */
/*----------------------------------------------------------------------------*/

void
cs_slope_test_gradient(int                         f_id,
                       cs_dispatch_context        &ctx,
                       int                         inc,
                       const cs_real_3_t          *grad,
                       cs_real_3_t                *grdpa,
                       const cs_real_t            *pvar,
                       const cs_field_bc_coeffs_t *bc_coeffs,
                       const cs_real_t            *i_massflux)
{
  CS_UNUSED(f_id);

  bool use_gpu = ctx.use_gpu();
  const cs_mesh_t  *m = cs_glob_mesh;

  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

#if defined(HAVE_ACCEL)

  if (use_gpu) {
    cs_device_context &d_ctx = static_cast<cs_device_context&>(ctx);

    _slope_test_gradient_d
      (d_ctx,
       inc,
       grad,
       grdpa,
       pvar,
       bc_coeffs,
       i_massflux);
  }

#endif

  if (use_gpu == false) {
    cs_host_context &h_ctx = static_cast<cs_host_context&>(ctx);

    _slope_test_gradient_h
      (h_ctx,
       inc,
       grad,
       grdpa,
       pvar,
       bc_coeffs,
       i_massflux);
  }

  ctx.wait();

  /* Synchronization for parallelism or periodicity */

  if (m->halo != nullptr) {
    cs_halo_sync_r(m->halo, CS_HALO_STANDARD, use_gpu, grdpa);
  }

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();

    std::chrono::microseconds elapsed;
    printf("%d: %s", cs_glob_rank_id, __func__);

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_stop - t_start);
    printf(", total = %ld\n", elapsed.count());
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the upwind gradient used in the pure SOLU schemes
 *        (observed in the litterature).
 *
 * \param[in]     f_id         field index
 * \param[in]     ctx          Reference to dispatch context
 * \param[in]     inc          Not an increment flag
 * \param[in]     halo_type    halo type
 * \param[in]     bc_coeffs    boundary condition structure for the variable
 * \param[in]     i_massflux   mass flux at interior faces
 * \param[in]     b_massflux   mass flux at boundary faces
 * \param[in]     pvar         values
 * \param[out]    grdpa        upwind gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_upwind_gradient(const int                     f_id,
                   cs_dispatch_context          &ctx,
                   const int                     inc,
                   const cs_halo_type_t          halo_type,
                   const cs_field_bc_coeffs_t   *bc_coeffs,
                   const cs_real_t               i_massflux[],
                   const cs_real_t               b_massflux[],
                   const cs_real_t     *restrict pvar,
                   cs_real_3_t         *restrict grdpa)
{
  CS_UNUSED(f_id);

  cs_real_t *coefap = bc_coeffs->a;
  cs_real_t *coefbp = bc_coeffs->b;

  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_halo_t  *halo = m->halo;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;
  const cs_real_t *restrict cell_vol = fvq->cell_vol;
  const cs_real_t *restrict i_face_surf = fvq->i_face_surf;
  const cs_real_t *restrict b_face_surf = fvq->b_face_surf;
  const cs_nreal_3_t *restrict i_face_u_normal = fvq->i_face_u_normal;
  const cs_nreal_3_t *restrict b_face_u_normal = fvq->b_face_u_normal;
  const int *restrict c_disable_flag = (fvq->has_disable_flag) ?
    fvq->c_disable_flag : nullptr;

  /* Parallel or device dispatch */
  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

    cs_lnum_t ii = i_face_cells[face_id][0];
    cs_lnum_t jj = i_face_cells[face_id][1];

    cs_real_t pif = pvar[ii];
    cs_real_t pjf = pvar[jj];

    cs_real_t pfac = (i_massflux[face_id] > 0.) ? pif : pjf;
    pfac *= i_face_surf[face_id];

    cs_real_t vfac_i[3], vfac_j[3];
    for (cs_lnum_t k = 0; k < 3; k++) {
      vfac_i[k] = pfac*i_face_u_normal[face_id][k];
      vfac_j[k] = - vfac_i[k];
    }

    cs_dispatch_sum<3>(grdpa[ii], vfac_i, i_sum_type);
    cs_dispatch_sum<3>(grdpa[jj], vfac_j, i_sum_type);

  });

  ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

    cs_lnum_t ii = b_face_cells[face_id];

    cs_real_t pfac = (b_massflux[face_id] < 0) ?
                      inc*coefap[face_id] + coefbp[face_id] * pvar[ii] :
                      pvar[ii];

    cs_real_t vfac[3];
    for (cs_lnum_t k = 0; k < 3; k++)
      vfac[k] = pfac * b_face_surf[face_id] * b_face_u_normal[face_id][k];

    cs_dispatch_sum<3>(grdpa[ii], vfac, b_sum_type);
  });

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {

    cs_real_t unsvol = cs_mq_cell_vol_inv(cell_id, c_disable_flag, cell_vol);

    grdpa[cell_id][0] *= unsvol;
    grdpa[cell_id][1] *= unsvol;
    grdpa[cell_id][2] *= unsvol;

  });

  ctx.wait();

  /* Synchronization for parallelism or periodicity */

  if (halo != nullptr)
    cs_halo_sync_r(halo, halo_type, ctx.use_gpu(), grdpa);
}

/*----------------------------------------------------------------------------*/
