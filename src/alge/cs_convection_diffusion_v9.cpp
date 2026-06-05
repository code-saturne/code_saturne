/*============================================================================
 * Convection-diffusion operators (v9.0 and older schemes).
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

/*-----------------------------------------------------------------------------
 * Standard C and C++ library headers
 *----------------------------------------------------------------------------*/

#include <chrono>

#include <assert.h>
#include <stdio.h>
#include <string.h>

#if defined(HAVE_MPI)
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_algorithm.h"
#include "base/cs_boundary_conditions.h"
#include "base/cs_boundary_conditions_set_coeffs.h"
#include "base/cs_dispatch.h"
#include "cdo/cs_equation_param.h"
#include "base/cs_halo.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "base/cs_field.h"
#include "base/cs_field_operator.h"
#include "base/cs_field_pointer.h"
#include "base/cs_parall.h"
#include "base/cs_porous_model.h"
#include "base/cs_profiling.h"
#include "base/cs_timer.h"

#include "alge/cs_gradient.h"
#include "mesh/cs_mesh_quantities.h"

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
  \brief Convection-diffusion operators, v9 scheme

  Please refer to the
  <a href="../../theory.pdf#conv-diff"><b>convection-diffusion</b></a> section
  of the theory guide for more information.

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
 * theory guide for more information.
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
 * \param[in]     pvar          solved variable
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
 * \param[in]     c_weight      diffusion gradient weighting
 * \param[in]     xcpp          array of specific heat (\f$ C_p \f$), or nullptr
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 * \param[in,out] i_flux        interior flux (or nullptr)
 * \param[in,out] b_flux        boundary flux (or nullptr)
 */
/*----------------------------------------------------------------------------*/

template <bool is_thermal, bool store_flux>
static bool
_convection_diffusion_scalar_unsteady_v91
  (const cs_field_t           *f,
   const cs_equation_param_t  &eqp,
   bool                        icvflb,
   int                         inc,
   int                         imasac,
   const cs_real_t            *restrict pvar,
   const int                   icvfli[],
   cs_field_bc_coeffs_t       *bc_coeffs,
   const cs_real_t             i_massflux[],
   const cs_real_t             b_massflux[],
   const cs_real_t             i_visc[],
   const cs_real_t             b_visc[],
   const cs_real_t            *c_weight,
   const cs_real_t             xcpp[],
   cs_real_t         *restrict rhs,
   cs_real_2_t       *restrict i_flux,
   cs_real_t         *restrict b_flux)
{
  const int iconvp = eqp.iconv;
  const int idiffp = eqp.idiff;
  const int ircflp = eqp.ircflu;
  const int ircflb = (ircflp > 0  && icvflb && is_thermal == false) ?
    eqp.b_diff_flux_rc : 0;
  const int ischcp = eqp.ischcv;
  const int isstpp = eqp.isstpc;
  const int iwarnp = eqp.verbosity;
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

  assert(eqp.ischcv < 3);

  /* Parallel or device dispatch */
  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  /* Local variables */

  char var_name[64];
  const int f_id = (f != nullptr) ? f->id : -1;

  bool pure_upwind = (blencp > 0.) ? false : true;
  cs_real_t *coface = nullptr, *cofbce = nullptr;

  cs_rreal_3_t *grad = nullptr;
  cs_rreal_3_t *gradup = nullptr;
  cs_rreal_3_t *gradst = nullptr;
  cs_real_2_t *bounds = nullptr;

  cs_real_t *local_min = nullptr;
  cs_real_t *local_max = nullptr;
  cs_real_t *courant = nullptr;

  cs_real_t *df_limiter = nullptr;

  cs_real_t  *v_slope_test = cs_get_v_slope_test(f_id,  eqp);

  int w_stride = 1;

  /* Initialization */

  /* Allocate work arrays */

  CS_MALLOC_HD(grad, n_cells_ext, cs_rreal_3_t, cs_alloc_mode);

  const cs_real_t rc_clip_factor
    = (eqp.ircflu != 0) ? eqp.rc_clip_factor : -1;
  if (rc_clip_factor >= 0)
    CS_MALLOC_HD(bounds, n_cells_ext, cs_real_2_t, cs_alloc_mode);

  /* Choose gradient type */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  cs_gradient_type_by_imrgra(eqp.d_gradient_r,
                             &gradient_type,
                             &halo_type);

  /* Limiters */

  if (f != nullptr) {
    int df_limiter_id = eqp.diffusion_limiter_id;
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

  /* Update BC coeffs */

  // Update BCs

  bool is_ischcp
    = (   (is_thermal && ischcp == 0)
       || (is_thermal == false && ischcp == 0));

  bool gradient_reco
    =    (   (idiffp != 0 && ircflp == 1)
      || (   iconvp != 0 && pure_upwind == false
          && (ircflp == 1 || isstpp == 0 || is_ischcp)));

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

  if (gradient_reco) {

    cs_gradient_scalar_synced_input(var_name,
                                    gradient_type,
                                    halo_type,
                                    inc,
                                    eqp.nswrgr,
                                    0, /* hyd_p_flag */
                                    w_stride,
                                    eqp.verbosity,
                                    (cs_gradient_limit_t)(eqp.d_imligr),
                                    eqp.epsrgr,
                                    eqp.d_climgr,
                                    nullptr, /* f_ext exterior force */
                                    bc_coeffs,
                                    pvar,
                                    c_weight, /* Weighted gradient */
                                    grad,
                                    bounds);

    /* Adjust bounds if reconstruction clip factor is >= 0 and != 1. */
    if (rc_clip_factor >= 0) {
      cs_convection_diffusion_adjust_and_check_bounds_scalar
        (ctx,
         var_name,
         eqp,
         true, /* face_gradient */
         ircflb,
         m,
         fvq,
         pvar,
         grad,
         df_limiter,
         bounds);
    }

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

      CS_MALLOC_HD(gradst, n_cells_ext, cs_rreal_3_t, cs_alloc_mode);

      /* Precomputed boundary face value */
      const cs_real_t *val_f = bc_coeffs->val_f;

      cs_slope_test_gradient(f_id,
                             ctx,
                             (const cs_rreal_3_t *)grad,
                             gradst,
                             pvar,
                             val_f,
                             i_massflux);

    }

    /* Pure SOLU scheme */
    if (ischcp == 2) {

      CS_MALLOC_HD(gradup, n_cells_ext, cs_rreal_3_t, cs_alloc_mode);

      ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t cell_id) {
        gradup[cell_id][0] = 0.;
        gradup[cell_id][1] = 0.;
        gradup[cell_id][2] = 0.;
      });

      cs_upwind_gradient(ctx,
                         inc,
                         halo_type,
                         bc_coeffs,
                         i_massflux,
                         b_massflux,
                         pvar,
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
        cs_rreal_t bldfrp = 1.;
        if (df_limiter != nullptr)  /* Local limitation of the reconstruction */
          bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                           0.);

        cs_rreal_t recoi, recoj;
        cs_i_compute_quantities(bldfrp,
                                diipf[face_id], djjpf[face_id],
                                grad[ii], grad[jj],
                                pvar[ii], pvar[jj],
                                &recoi, &recoj,
                                &pip, &pjp);

        if (bounds != nullptr) {
          cs_clip_quantity(bounds[ii], pip);
          cs_clip_quantity(bounds[jj], pjp);
        }
      }
      else {
        pip = pvar[ii];
        pjp = pvar[jj];
      }

      /* No relaxation in following expressions (pipr = pip, pjpr = pjp) */

      /* Convective flux (as we are upwind, pif == pvar[ii],
                                            pjf == pvar[jj] */

      if (iconvp == 1) {
        cs_real_t _i_massflux = i_massflux[face_id];
        cs_real_t flui = 0.5*(_i_massflux + cs::abs(_i_massflux));
        cs_real_t fluj = 0.5*(_i_massflux - cs::abs(_i_massflux));

        fluxi += cpi*(  thetap * (flui*pvar[ii] + fluj*pvar[jj])
                      - imasac * _i_massflux*pvar[ii]);

        fluxj += cpj*(  thetap * (flui*pvar[ii] + fluj*pvar[jj])
                      - imasac *  _i_massflux*pvar[jj]);

        /* Fluxes without mass accumulation */
        if (store_flux) {
          i_flux[face_id][0] += cpi*thetap*(flui*pvar[ii] + fluj*pvar[jj]);
          i_flux[face_id][1] += cpj*thetap*(flui*pvar[ii] + fluj*pvar[jj]);
        }
      }

      /* Diffusive flux */

      if (idiffp > 0) {
        cs_real_t diff_contrib = thetap*i_visc[face_id]*(pip - pjp);
        fluxi += diff_contrib;
        fluxj += diff_contrib;

        /* Fluxes if needed */
        if (store_flux) {
          i_flux[face_id][0] += diff_contrib;
          i_flux[face_id][1] += diff_contrib;
        }
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

    if (ischcp < 0 || ischcp > 3) {
      bft_error(__FILE__, __LINE__, 0,
                _("invalid value of ischcv"));
    }

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t cpi = 1.0, cpj = 1.0;
      if (is_thermal) {
        cpi = xcpp[ii];
        cpj = xcpp[jj];
      }

      cs_real_t pif, pjf;
      cs_real_t pip, pjp;

#if defined(__INTEL_LLVM_COMPILER)
      // Silence unitialized variables warning due do compiler ignoring
      // initializations in inlined functions.
      pif = 0., pjf = 0.;
#endif

      cs_real_t fluxi = 0., fluxj = 0.;

      if (ircflp == 1) {
        cs_rreal_t bldfrp = 1.;
        if (df_limiter != nullptr)  /* Local limiter of the reconstruction */
          bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                           0.);

        cs_rreal_t recoi, recoj;
        cs_i_compute_quantities(bldfrp,
                                diipf[face_id], djjpf[face_id],
                                grad[ii], grad[jj],
                                pvar[ii], pvar[jj],
                                &recoi, &recoj,
                                &pip, &pjp);

        if (bounds != nullptr) {
          cs_clip_quantity(bounds[ii], pip);
          cs_clip_quantity(bounds[jj], pjp);
        }
      }
      else {
        pip = pvar[ii];
        pjp = pvar[jj];
      }

      const cs_real_t *cell_ceni = cell_cen[ii];
      const cs_real_t *cell_cenj = cell_cen[jj];
      const cs_real_t w_f = weight[face_id];

      if (ischcp == 0) {

        /* Legacy SOLU
           -----------*/

        cs_solu_f_val(cell_ceni,
                      i_face_cog[face_id],
                      grad[ii],
                      pvar[ii],
                      &pif);
        cs_solu_f_val(cell_cenj,
                      i_face_cog[face_id],
                      grad[jj],
                      pvar[jj],
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
                      pvar[ii],
                      &pif);
        cs_solu_f_val(cell_cenj,
                      i_face_cog[face_id],
                      gradup[jj],
                      pvar[jj],
                      &pjf);

      }

      // Convective flux

      if (iconvp == 1) {
        cs_real_t _i_massflux = i_massflux[face_id];
        cs_real_t flui = 0.5*(_i_massflux + cs::abs(_i_massflux));
        cs_real_t fluj = 0.5*(_i_massflux - cs::abs(_i_massflux));

        fluxi += cpi*(  thetap*(flui*pif + fluj*pjf)
                      - imasac*_i_massflux*pvar[ii]);

        fluxj += cpj*(  thetap*(flui*pif + fluj*pjf)
                      - imasac*_i_massflux*pvar[jj]);

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
        cs_rreal_t bldfrp = 1.;
        if (df_limiter != nullptr)  /* Local limiter */
          bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                           0.);

        cs_rreal_t recoi, recoj;
        cs_i_compute_quantities(bldfrp,
                                diipf[face_id], djjpf[face_id],
                                grad[ii], grad[jj],
                                pvar[ii], pvar[jj],
                                &recoi, &recoj,
                                &pip, &pjp);

        if (bounds != nullptr) {
          cs_clip_quantity(bounds[ii], pip);
          cs_clip_quantity(bounds[jj], pjp);
        }
      }
      else {
        pip = pvar[ii];
        pjp = pvar[jj];
      }

      const cs_real_t *cell_ceni = cell_cen[ii];
      const cs_real_t *cell_cenj = cell_cen[jj];
      const cs_real_t w_f = weight[face_id];

      /* Slope test is needed with convection */
      if (iconvp > 0) {
        cs_rreal_t testij, tesqck;

        cs_slope_test(pvar[ii],
                      pvar[jj],
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
                        pvar[ii],
                        &pif);
          cs_solu_f_val(cell_cenj,
                        i_face_cog[face_id],
                        grad[jj],
                        pvar[jj],
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
                        pvar[ii],
                        &pif);
          cs_solu_f_val(cell_cenj,
                        i_face_cog[face_id],
                        gradup[jj],
                        pvar[jj],
                        &pjf);

        }

        /* Slope test: percentage of upwind
           -------------------------------- */

        if (tesqck <= 0. || testij <= 0.) {

          cs_blend_f_val(blend_st, pvar[ii], &pif);
          cs_blend_f_val(blend_st, pvar[jj], &pjf);

          upwind_switch = true;

        }

        /* Blending
           -------- */

        cs_blend_f_val(blencp, pvar[ii], &pif);
        cs_blend_f_val(blencp, pvar[jj], &pjf);

      }
      else { /* If iconv=0 p*fr* are useless */

        pif = pvar[ii];
        pjf = pvar[jj];

      } /* End for slope test */

      // Convective flux

      if (iconvp == 1) {
        cs_real_t _i_massflux = i_massflux[face_id];
        cs_real_t flui = 0.5*(_i_massflux + cs::abs(_i_massflux));
        cs_real_t fluj = 0.5*(_i_massflux - cs::abs(_i_massflux));

        fluxi += cpi*(  thetap*(flui*pif + fluj*pjf)
                      - imasac*_i_massflux*pvar[ii]);

        fluxj += cpj*(  thetap*(flui*pif + fluj*pjf)
                      - imasac*_i_massflux*pvar[jj]);

        /* Fluxes without mass accumulation if needed */
        if (store_flux) {
          i_flux[face_id][0] +=  cpi*thetap*(flui*pif + fluj*pjf);
          i_flux[face_id][1] +=  cpj*thetap*(flui*pif + fluj*pjf);
        }

      }

      // Diffusive flux (no relaxation)

      if (idiffp > 0) {
        cs_real_t diff_contrib = thetap*i_visc[face_id]*(pip - pjp);
        fluxi += diff_contrib;
        fluxj += diff_contrib;

        /* Fluxes if needed */
        if (store_flux) {
          i_flux[face_id][0] += diff_contrib;
          i_flux[face_id][1] += diff_contrib;
        }
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
      n_upwind = cs::algorithm::count_reduce_sum(ctx, n_i_faces, i_upwind);
      cs_parall_counter(&n_upwind, 1);
    }
    else if (pure_upwind)
      n_upwind = m->n_g_i_faces;

    bft_printf(_(" %s: %llu Faces with upwind on %llu interior faces\n"),
               var_name, (unsigned long long)n_upwind,
               (unsigned long long)m->n_g_i_c_faces);

    CS_FREE(i_upwind);
  }

  /* ======================================================================
     ---> Contribution from boundary faces
     ======================================================================*/

   const cs_real_t *val_f_g = bc_coeffs->val_f;
   const cs_real_t *flux_d = bc_coeffs->flux_diff;

   /* Boundary convective flux are all computed with an upwind scheme */
  if (icvflb == false || is_thermal) {

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {

      cs_lnum_t ii = b_face_cells[face_id];

      cs_real_t cpi = 1.;
      if (is_thermal)
        cpi = xcpp[ii];

      cs_real_t fluxi = 0.;

      cs_b_upwind_flux(iconvp,
                       thetap,
                       imasac,
                       bc_type[face_id],
                       pvar[ii],
                       pvar[ii], /* no relaxation */
                       val_f_g[face_id],
                       b_massflux[face_id],
                       cpi,
                       &fluxi);

      if (idiffp > 0) {
        cs_b_diff_flux(1.0, // idiffp
                       thetap,
                       flux_d[face_id],
                       b_visc[face_id],
                       &fluxi);
      }

      cs_dispatch_sum(&rhs[ii], -fluxi, b_sum_type);

      /* Fluxes without mass accumulation if needed */
      if (store_flux) {
        fluxi = 0.;

        cs_b_upwind_flux(iconvp,
                         thetap,
                         0, //imasac,
                         bc_type[face_id],
                         pvar[ii],
                         pvar[ii], /* no relaxation */
                         val_f_g[face_id],
                         b_massflux[face_id],
                         cpi,
                         &fluxi);

        if (idiffp > 0) {
          cs_b_diff_flux(1.0, // idiffp,
                         thetap,
                         flux_d[face_id],
                         b_visc[face_id],
                         &fluxi);
        }

        b_flux[face_id] += fluxi;
      }

    });
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

      cs_rreal_t bldfrp = (cs_real_t) ircflb;
      /* Local limitation of the reconstruction */
      if (df_limiter != nullptr && ircflb > 0)
        bldfrp = cs::max(df_limiter[ii], 0.);

      cs_b_cd_unsteady(bldfrp,
                       diipb[face_id],
                       grad[ii],
                       pvar[ii],
                       &pip);

      if (bounds != nullptr)
        cs_clip_quantity(bounds[ii], pip);

      cs_b_imposed_conv_flux(iconvp,
                             thetap,
                             imasac,
                             inc,
                             bc_type[face_id],
                             icvfli[face_id],
                             pvar[ii],
                             pvar[ii], /* no relaxation */
                             pip,
                             coface[face_id],
                             cofbce[face_id],
                             b_massflux[face_id],
                             1., /* xcpp */
                             val_f_g[face_id],
                             &fluxi);

      cs_b_diff_flux(idiffp,
                     thetap,
                     flux_d[face_id],
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
                               pvar[ii],
                               pvar[ii], /* no relaxation */
                               pip,
                               coface[face_id],
                               cofbce[face_id],
                               b_massflux[face_id],
                               1., /* xcpp */
                               val_f_g[face_id],
                               &fluxi);

        cs_b_diff_flux(idiffp,
                       thetap,
                       flux_d[face_id],
                       b_visc[face_id],
                       &fluxi);

        b_flux[face_id] += fluxi;
      }

    });

  }

  ctx.wait();

  /* Free memory */
  CS_FREE(bounds);
  CS_FREE(grad);
  CS_FREE(gradup);
  CS_FREE(gradst);
  CS_FREE(local_max);
  CS_FREE(local_min);
  CS_FREE(courant);

  return true;
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
 * \param[in]     bc_coeffs     boundary conditions structure for the variable
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     i_pvar        velocity at interior faces
 * \param[in]     b_pvar        velocity at boundary faces
 * \param[in]     grad          associated gradient
 * \param[in]     bounds        optional bounds (square distance)
 *                              of values in adjacent cells and faces, or null
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride, bool porous_vel, typename T>
static void
_convection_diffusion_unsteady_strided_v91
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
   const cs_real_t              i_massflux[],
   const cs_real_t              b_massflux[],
   const cs_real_t              i_visc[],
   const cs_real_t              b_visc[],
   cs_real_t         (*restrict i_pvar)[stride],
   cs_real_t         (*restrict b_pvar)[stride],
   T                 (*restrict grad)[stride][3],
   cs_real_t          *restrict bounds,
   cs_real_t         (*restrict rhs)[stride])
{
  using rgrad_t = T[stride][3];
  using var_t = cs_real_t[stride];
  using b_t = cs_real_t[stride][stride];

  /* Precomputed boundary face value */
  const var_t *val_f = (var_t *)bc_coeffs->val_f;
  const var_t *flux = (var_t *)bc_coeffs->flux_diff;

  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

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

  int df_limiter_id = eqp.diffusion_limiter_id;
  if (df_limiter_id > -1)
    df_limiter = cs_field_by_id(df_limiter_id)->val;

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

  cs_alloc_mode_t amode = ctx.alloc_mode();

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

  /* Adjust bounds if reconstruction clip factor is >= 0 and != 1. */
  if (eqp.rc_clip_factor >= 0) {
    cs_convection_diffusion_adjust_and_check_bounds_strided
      (ctx,
       var_name,
       eqp,
       true, // face gradient reconstruction
       ircflb,
       m,
       fvq,
       grad,
       df_limiter,
       bounds);
  }

  /* Compute the balance with reconstruction */

  /* ======================================================================
     Compute uncentered gradient grdpa for the slope test
     ======================================================================*/

  rgrad_t *grdpa = nullptr;

  if (iconvp > 0 && pure_upwind == false && isstpp == 0) {
    CS_MALLOC_HD(grdpa, n_cells_ext, rgrad_t, amode);

    cs_slope_test_gradient_strided<stride>(ctx,
                                           (const rgrad_t *)grad,
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

  ctx.wait();

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
        T bldfrp = T{1};
        if (df_limiter != nullptr)  /* Local limiter */
          bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                           0.);

        T recoi[stride], recoj[stride];
        cs_i_compute_quantities_strided<stride>(bldfrp,
                                                diipf[face_id], djjpf[face_id],
                                                grad[ii], grad[jj],
                                                _pi, _pj,
                                                recoi, recoj,
                                                pip, pjp);

        if (bounds != nullptr) {
          cs_clip_quantity_strided<stride>(bounds[ii], _pi, pip);
          cs_clip_quantity_strided<stride>(bounds[jj], _pj, pjp);
        }
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
          T bldfrp = T{1};
          if (df_limiter != nullptr)  /* Local limiter */
            bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                             0.);

          T recoi[stride], recoj[stride];
          cs_i_compute_quantities_strided<stride>(bldfrp,
                                                  diipf[face_id], djjpf[face_id],
                                                  grad[ii], grad[jj],
                                                  _pi, _pj,
                                                  recoi, recoj,
                                                  pip, pjp);

          if (bounds != nullptr) {
            cs_clip_quantity_strided<stride>(bounds[ii], _pi, pip);
            cs_clip_quantity_strided<stride>(bounds[jj], _pj, pjp);
          }
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
          T bldfrp = T{1};
          if (df_limiter != nullptr)  /* Local limiter */
            bldfrp = cs::max(cs::min(df_limiter[ii], df_limiter[jj]),
                             0.);

          T recoi[stride], recoj[stride];
          cs_i_compute_quantities_strided<stride>(bldfrp,
                                                  diipf[face_id], djjpf[face_id],
                                                  grad[ii], grad[jj],
                                                  _pi, _pj,
                                                  recoi, recoj,
                                                  pip, pjp);

          if (bounds != nullptr) {
            cs_clip_quantity_strided<stride>(bounds[ii], _pi, pip);
            cs_clip_quantity_strided<stride>(bounds[jj], _pj, pjp);
          }
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
          T testij, tesqck;

          const rgrad_t &gradi = grad[ii];
          const rgrad_t &gradj = grad[jj];

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
          if (cs::min(tesqck, testij) <= T{0}) {
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
      n_upwind = cs::algorithm::count_reduce_sum(ctx, n_i_faces, i_upwind);
      cs_parall_counter(&n_upwind, 1);
    }
    else if (pure_upwind)
      n_upwind = m->n_g_i_faces;

    bft_printf(_(" %s: %llu Faces with upwind on %llu interior faces\n"),
               var_name, (unsigned long long)n_upwind,
               (unsigned long long)m->n_g_i_c_faces);

    CS_FREE(i_upwind);
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
        b_val_d[isou] = flux[face_id][isou];
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
        b_val_d[isou] = flux[face_id][isou];
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

        T bldfrp = (T)ircflb;
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
  CS_FREE(grdpa);

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

/*============================================================================
 * Public function definitions
 *============================================================================*/

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
 * theory guide for more information.
 *
 * \param[in]     f             pointer to field, or null
 * \param[in]     eqp           equation parameters
 * \param[in]     icvflb        global indicator of boundary convection flux
 *                               - 0 upwind scheme at all boundary faces
 *                               - 1 imposed flux at some boundary faces
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imasac        take mass accumulation into account?
 * \param[in]     pvar          solved variable
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
 * \param[in]     c_weight      diffusion gradient weighting
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 * \param[in,out] i_flux        interior flux (or nullptr)
 * \param[in,out] b_flux        boundary flux (or nullptr)
 *
 * \return true if v9 algorithm called, false if nothing done
 */
/*----------------------------------------------------------------------------*/

bool
cs_convection_diffusion_scalar_v9(const cs_field_t           *f,
                                  const cs_equation_param_t   eqp,
                                  int                         icvflb,
                                  int                         inc,
                                  int                         imasac,
                                  const cs_real_t            *restrict pvar,
                                  const int                   icvfli[],
                                  cs_field_bc_coeffs_t       *bc_coeffs,
                                  const cs_real_t             i_massflux[],
                                  const cs_real_t             b_massflux[],
                                  const cs_real_t             i_visc[],
                                  const cs_real_t             b_visc[],
                                  const cs_real_t            *c_weight,
                                  cs_real_t                  *rhs,
                                  cs_real_2_t                 i_flux[],
                                  cs_real_t                   b_flux[])
{
  /* TVD/NVD schemes rarely used except for specific VoF case, so we
     do not maintain compatibility for those.
     Hybrid blend not commonly used either. */

  if (eqp.ischcv == 3 || eqp.ischcv == 4)
    return false;

  CS_PROFILE_FUNC_RANGE();

  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  bool retval = 0;

  if (i_flux == nullptr)
    retval = _convection_diffusion_scalar_unsteady_v91<false, false>
               (f, eqp, icvflb, inc, imasac,
                pvar,
                icvfli,
                bc_coeffs,
                i_massflux, b_massflux,
                i_visc, b_visc, c_weight,
                nullptr, rhs, i_flux, b_flux);
  else
    retval = _convection_diffusion_scalar_unsteady_v91<false, true>
               (f, eqp, icvflb, inc, imasac,
                pvar,
                icvfli,
                bc_coeffs,
                i_massflux, b_massflux,
                i_visc, b_visc, c_weight,
                nullptr, rhs, i_flux, b_flux);

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();

    std::chrono::microseconds elapsed;
    printf("%d: %s", cs_glob_rank_id, __func__);

    elapsed = std::chrono::duration_cast
      <std::chrono::microseconds>(t_stop - t_start);
    printf(", total = %ld\n", elapsed.count());
  }

  return retval;
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
 * \warning \f$ Rhs \f$ must have been initialized before calling this function!
 * \warning The ghost cell values of pvar must already have been synchronized.
 *
 * \param[in]     f             pointer to field, or null
 * \param[in]     eqp           equation parameters)
 * \param[in]     inc           indicator
 *                               - 0 when solving an increment
 *                               - 1 otherwise
 * \param[in]     imasac        take mass accumulation into account?
 * \param[in]     pvar          solved variable
 * \param[in]     bc_coeffs     boundary condition structure for the variable
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in]     w_stride      stride for weighting coefficient
 * \param[in]     c_weight      diffusion gradient weighting
 * \param[in]     xcpp          array of specific heat (\f$ C_p \f$)
 * \param[in,out] rhs           right hand side \f$ \vect{Rhs} \f$
 *
 * \return true if v9 algorithm called, false if nothing done
 */
/*----------------------------------------------------------------------------*/

bool
cs_convection_diffusion_thermal_v9(const cs_field_t           *f,
                                   const cs_equation_param_t   eqp,
                                   int                         inc,
                                   int                         imasac,
                                   const cs_real_t           * pvar,
                                   cs_field_bc_coeffs_t       *bc_coeffs,
                                   const cs_real_t             i_massflux[],
                                   const cs_real_t             b_massflux[],
                                   const cs_real_t             i_visc[],
                                   const cs_real_t             b_visc[],
                                   const cs_real_t            *c_weight,
                                   const cs_real_t             xcpp[],
                                   cs_real_t        *restrict  rhs)
{
  /* TVD/NVD schemes rarely used except for specific VoF case, so we
     do not maintain compatibility for those.
     Hybrid blend not commonly used either. */

  if (eqp.ischcv == 3 || eqp.ischcv == 4)
    return false;

  CS_PROFILE_FUNC_RANGE();

  bool retval = _convection_diffusion_scalar_unsteady_v91<true, false>
                  (f, eqp,
                   false, /* icvflb */
                   inc, imasac,
                   pvar,
                   nullptr, /* icvfli */
                   bc_coeffs,
                   i_massflux, b_massflux,
                   i_visc, b_visc, c_weight,
                   xcpp, rhs, nullptr, nullptr);

  return retval;
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
 * \param[in]     bc_coeffs     boundary conditions structure for the variable
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
 *
 * \return true if v9 algorithm called, false if nothing done
 */
/*----------------------------------------------------------------------------*/

bool
cs_convection_diffusion_vector_v9(int                         idtvar,
                                  int                         f_id,
                                  const cs_equation_param_t   eqp,
                                  int                         icvflb,
                                  int                         inc,
                                  int                         ivisep,
                                  int                         imasac,
                                  cs_real_3_t       *restrict pvar,
                                  const cs_real_3_t *restrict pvara,
                                  const int                   icvfli[],
                                  cs_field_bc_coeffs_t       *bc_coeffs,
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
  CS_PROFILE_FUNC_RANGE();

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

  cs_alloc_mode_t amode = ctx.alloc_mode();

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

  using grad_t = cs_rreal_t[3][3];
  grad_t *grad;
  CS_MALLOC_HD(grad, n_cells_ext, grad_t, amode);

  cs_real_t *bounds = nullptr;

  const cs_real_3_t  *restrict _pvar
    = (pvar != nullptr) ? (const cs_real_3_t  *)pvar : pvara;

  /* Update BC coeffs */
  cs_boundary_conditions_update_bc_coeff_face_values_strided<3>
    (ctx, f, bc_coeffs, inc, &eqp, _pvar);

  if (  (eqp.idiff != 0 && eqp.ircflu == 1)
      || ivisep == 1
      || (    eqp.iconv != 0 && eqp.blencv > 0.
          && (eqp.ischcv == 0 || eqp.ircflu == 1 || eqp.isstpc == 0))) {

    const cs_real_t rc_clip_factor
      = (eqp.ircflu != 0) ? eqp.rc_clip_factor : -1;
    if (rc_clip_factor >= 0)
      CS_MALLOC_HD(bounds, n_cells_ext, cs_real_t, cs_alloc_mode);

    cs_real_t *gweight = nullptr;
    if (f_id != -1) {
      /* Get the calculation option from the field */
      if (f->type & CS_FIELD_VARIABLE && eqp.iwgrec == 1) {
        if (eqp.idiff > 0) {
          int key_id = cs_field_key_id("gradient_weighting_id");
          int diff_id = f->get_key_int(key_id);
          if (diff_id > -1) {
            cs_field_t *weight_f = cs_field_by_id(diff_id);
            gweight = weight_f->val;
            cs_field_synchronize(weight_f, halo_type);
          }
        }
      }
    }

    cs_gradient_vector_synced_input(var_name,
                                    gradient_type,
                                    halo_type,
                                    inc,
                                    eqp.nswrgr,
                                    eqp.verbosity,
                                    (cs_gradient_limit_t)(eqp.d_imligr),
                                    eqp.epsrgr,
                                    eqp.d_climgr,
                                    bc_coeffs,
                                    (const cs_real_3_t *)_pvar,
                                    gweight, /* weighted gradient */
                                    grad,
                                    bounds);
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
      _convection_diffusion_unsteady_strided_v91<3, false>
        (ctx, f, var_name, eqp, icvflb, inc, imasac,
         pvar, pvara,
         icvfli,
         bc_coeffs,
         i_massflux, b_massflux,
         i_visc, b_visc,
         i_pvar, b_pvar, grad, bounds, rhs);
    else
      _convection_diffusion_unsteady_strided_v91<3, true>
        (ctx, f, var_name, eqp, icvflb, inc, imasac,
         pvar, pvara,
         icvfli,
         bc_coeffs,
         i_massflux, b_massflux,
         i_visc, b_visc,
         i_pvar, b_pvar, grad, bounds, rhs);
  }

  else {
    cs_convection_diffusion_steady_strided<3>
      (f, var_name, eqp, icvflb, inc,
       pvar, pvara,
       icvfli,
       bc_coeffs,
       i_massflux, b_massflux,
       i_visc, b_visc,
       (cs_real_33_t *)grad, rhs);
  }

  /* Computation of the transpose grad(vel) term and grad(-2/3 div(vel))
     ------------------------------------------------------------------- */

  if (ivisep == 1 && eqp.idiff == 1)
    cs_convection_diffusion_secvis(ctx,
                                   m,
                                   cs_glob_mesh_quantities,
                                   eqp.theta,
                                   i_visc,
                                   i_secvis,
                                   b_secvis,
                                   grad,
                                   rhs);

  ctx.wait();

  /* Free memory */
  CS_FREE(bounds);
  CS_FREE(grad);

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();

    std::chrono::microseconds elapsed;
    printf("%d: %s", cs_glob_rank_id, __func__);

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_stop - t_start);
    printf(", total = %ld\n", elapsed.count());
  }

  return true;
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
 * \param[in]     bc_coeffs     boundary condition structure for the variable
 * \param[in]     i_massflux    mass flux at interior faces
 * \param[in]     b_massflux    mass flux at boundary faces
 * \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
 *                               at interior faces for the r.h.s.
 * \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
 *                               at border faces for the r.h.s.
 * \param[in,out] rhs           right hand side \f$ \tens{Rhs} \f$
 *
 * \return true if v9 algorithm called, false if nothing done
 */
/*----------------------------------------------------------------------------*/

bool
cs_convection_diffusion_tensor_v9(int                          idtvar,
                                  int                          f_id,
                                  const cs_equation_param_t    eqp,
                                  int                          icvflb,
                                  int                          inc,
                                  int                          imasac,
                                  cs_real_6_t        *restrict pvar,
                                  const cs_real_6_t  *restrict pvara,
                                  cs_field_bc_coeffs_t        *bc_coeffs,
                                  const cs_real_t              i_massflux[],
                                  const cs_real_t              b_massflux[],
                                  const cs_real_t              i_visc[],
                                  const cs_real_t              b_visc[],
                                  cs_real_6_t        *restrict rhs)
{
  CS_PROFILE_FUNC_RANGE();

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

  cs_alloc_mode_t amode = ctx.alloc_mode();

  /* Halo (and gradient) type */

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;
  cs_gradient_type_by_imrgra(eqp.d_gradient_r,
                             &gradient_type,
                             &halo_type);

  /* Handle cases where only the previous values (already synchronized)
     or current values are provided */

  if (pvar != nullptr)
    cs_halo_sync_r(m->halo, halo_type, ctx.use_gpu(), pvar);
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

  cs_real_t *bounds = nullptr;

  if (  (eqp.idiff != 0 && eqp.ircflu == 1)
      || (   eqp.iconv != 0 && eqp.blencv > 0.
          && (eqp.ischcv == 0 || eqp.ircflu == 1 || eqp.isstpc == 0))) {

    if (eqp.rc_clip_factor >= 0)
      CS_MALLOC_HD(bounds, n_cells_ext, cs_real_t, cs_alloc_mode);

    const cs_real_6_t  *restrict _pvar
      = (pvar != nullptr) ? (const cs_real_6_t *)pvar : pvara;

    cs_gradient_tensor_synced_input(var_name,
                                    gradient_type,
                                    halo_type,
                                    inc,
                                    eqp.nswrgr,
                                    eqp.verbosity,
                                    (cs_gradient_limit_t)(eqp.d_imligr),
                                    eqp.epsrgr,
                                    eqp.d_climgr,
                                    bc_coeffs,
                                    _pvar,
                                    (cs_real_63_t *)grad,
                                    bounds);

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
    _convection_diffusion_unsteady_strided_v91<6, false>
      (ctx, f, var_name, eqp, 0, inc, imasac,
       pvar, pvara,
       nullptr, // icvfli,
       bc_coeffs,
       i_massflux, b_massflux,
       i_visc, b_visc,
       nullptr, nullptr, grad, bounds, rhs);
  }
  else {
    cs_convection_diffusion_steady_strided<6>
      (f, var_name, eqp, icvflb, inc,
       pvar, pvara, nullptr,
       bc_coeffs,
       i_massflux, b_massflux,
       i_visc, b_visc,
       grad, rhs);
  }

  /* Free memory */
  CS_FREE(bounds);
  CS_FREE(grad);

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();

    std::chrono::microseconds elapsed;
    printf("%d: %s", cs_glob_rank_id, __func__);

    elapsed = std::chrono::duration_cast
                <std::chrono::microseconds>(t_stop - t_start);
    printf(", total = %ld\n", elapsed.count());
  }

  return true;
}

/*----------------------------------------------------------------------------*/