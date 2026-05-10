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

#include "base/cs_algorithm.h"
#include "base/cs_array.h"
#include "base/cs_boundary_conditions.h"
#include "base/cs_boundary_conditions_set_coeffs.h"
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
#include "base/cs_parall.h"
#include "base/cs_parameters.h"
#include "base/cs_porous_model.h"
#include "base/cs_profiling.h"
#include "base/cs_reducers.h"
#include "base/cs_timer.h"
#include "base/cs_velocity_pressure.h"

#include "alge/cs_blas.h"
#include "alge/cs_bad_cells_regularisation.h"
#include "alge/cs_gradient.h"
#include "alge/cs_gradient_boundary.h"
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
  \brief Convection-diffusion operators.

  Please refer to the
  <a href="../../theory.pdf#conv-diff"><b>convection-diffusion</b></a> section
  of the theory guide for more information.

  \section rc_bounds Bounds for diffusion scheme value reconstruction.

  For non-orthogonal meshes, the diffusion scheme involves reconstruction
  of cell variables at points orthogonal to the associated face centers.
  Deactivating this reconstruction improves robustness, at the cost of the
  loss of the scheme's second order spatial accuracy and even consistency.

  Since the reconstruction (based on 1st order Taylor expansion using
  the cell gradients or face gradients (mean or adjacent cell gradients)
  can lead to reconstructed values outside the range of values in neighboring
  cells, diffusion from a cell _I_ with a given field value to a cell _J_
  with a higher value is possible if the reconstructed value at the _I'_
  location associated to face _IJ_ is higher than the value in J, leading
  to loss of monotonicity.

  To avoid this, using gradient limiters can reduce the reconstructed
  values range where computed gradients are detected as being too strong.
  It is now also possible to use *bounded reconstruction*, where the gradient
  value is not modified, but values outside the range of adjacent cell
  (and boundary face) values are clipped to values within this range.

  For scalars, bounds are naturally based on the minimum and maximum values
  at adjacent cells and boundary faces. the actual upper bound  applied in
  cell _I_ is equal to
  \f$ \varia_clip = \varia_\celli + rc{\_}clip{\_}factor . (\varia_{max} - \varia_\celli) \f$
  so using a reconstruction clipping factor of 0 is equivalent to disabling
  reconstruction, and using a clipping factor of 1 limits the reconstructed
  value to the exact bounds of adjacent cells (using face adjacency only
  or additional cell neighbors based on the cell neighborhood used for
  the matching reconstruction gradient. The lower bound is defined in a
  similar manner.

  For vectors and tensors, bounds are not so straightforward, as
  component-based bounds are not invariant by rotation.
  So we choose a simpled dispersion bound for vectors, for which the
  Euclidean norm of the difference between the reconstructed vector and
  the non-reconstructed vector in cell I does not exceed that of if the
  norm of this difference with all adjacent cell and boundary face values.
  If this bound is exceeded, the reconstruction gradient is scaled so
  that the reconstructed value lies withing this bound
  (here again a reconstruction clipping factor multiplier may be used,
  so a factor of 0 disables reconstruction entirely, while a factor of 1
  limits the reconstruction so as to exactly fit the bound).

  This is illustrated in the following figures:
  \image html vector_bounds_a.svg "Vector values in cell _I_ and adjacent cells." width=40%

  Placing vectors at the same origin, we can compute their dispersion:
  \image html vector_bounds_b.svg "Vector values in cell _I_ and adjacent cells."  width=40%

  Finally, the value reconstructed at _I'_ should fir inside both circles:
  \image html vector_bounds_c.svg "Vector value bounds at _I'_." width=40%

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

  bool pure_upwind = (blencp > 0.) ? false : true;
  cs_real_t *coface = nullptr, *cofbce = nullptr;

  cs_rreal_3_t *grad = nullptr;
  cs_rreal_3_t *gradup = nullptr;
  cs_rreal_3_t *gradst = nullptr;
  cs_real_2_t *bounds = nullptr;

  cs_real_t *local_min = nullptr;
  cs_real_t *local_max = nullptr;
  cs_real_t *courant = nullptr;

  cs_real_t *cv_limiter = nullptr;
  cs_real_t *df_limiter = nullptr;

  const int key_lim_choice = cs_field_key_id("limiter_choice");

  cs_real_t  *v_slope_test = cs_get_v_slope_test(f_id,  eqp);

  int w_stride = 1;

  /* Initialization */

  /* Allocate work arrays */

  CS_MALLOC_HD(grad, n_cells_ext, cs_real_3_t, cs_alloc_mode);

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

    /* NVD/TVD limiters */
    if (ischcp == 4) {
      limiter_choice = (cs_nvd_type_t)(cs_field_get_key_int(f, key_lim_choice));
      CS_MALLOC_HD(local_max, n_cells_ext, cs_real_t, cs_alloc_mode);
      CS_MALLOC_HD(local_min, n_cells_ext, cs_real_t, cs_alloc_mode);
      // FIXME replace this with bounds such as those computed with gradients.
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

    int cv_limiter_id = eqp.convection_limiter_id;
    if (cv_limiter_id > -1)
      cv_limiter = cs_field_by_id(cv_limiter_id)->val;

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
       || (is_thermal == false && (ischcp == 0 || ischcp == 3 || ischcp == 4)));

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
         ircflb,
         true, /* face_gradient */
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
    if (ischcp == 2 || (is_thermal && ischcp == 4)) {

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

      if (ischcp != 4) {

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
                        pvar[ii],
                        &pif_up);

          cs_solu_f_val(cell_cenj,
                        i_face_cog[face_id],
                        grad[jj],
                        pvar[jj],
                        &pjf_up);

          cs_real_t hybrid_blend_interp
            = cs::min(hybrid_blend[ii], hybrid_blend[jj]);

          pif = hybrid_blend_interp*pif + (1. - hybrid_blend_interp)*pif_up;
          pjf = hybrid_blend_interp*pjf + (1. - hybrid_blend_interp)*pjf_up;

        }

        /* Blending
           -------- */

        pif = beta * pif + (1. - beta) * pvar[ii];
        pjf = beta * pjf + (1. - beta) * pvar[jj];

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
        if (courant != nullptr && is_thermal == false)
          courant_c = courant[ic];

        cs_i_cd_unsteady_nvd(limiter_choice,
                             beta,
                             cell_cen[ic],
                             cell_cen[id],
                             i_face_u_normal[face_id],
                             i_face_cog[face_id],
                             grad[ic],
                             pvar[ic],
                             pvar[id],
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
