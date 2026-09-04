/*============================================================================
 * Turbulence transport equation.
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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>

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
#include "base/cs_dispatch.h"
#include "alge/cs_divergence.h"
#include "base/cs_equation_iterative_solve.h"
#include "cdo/cs_equation_param.h"
#include "alge/cs_face_viscosity.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_operator.h"
#include "base/cs_field_pointer.h"
#include "base/cs_log_iteration.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_location.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_parall.h"
#include "base/cs_physical_constants.h"
#include "base/cs_prototypes.h"
#include "base/cs_thermal_model.h"
#include "base/cs_solid_zone.h"
#include "base/cs_time_step.h"
#include "turb/cs_turbulence_bc.h"
#include "turb/cs_turbulence_model.h"
#include "base/cs_velocity_pressure.h"

#include "turb/cs_turbulence_rij.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "turb/cs_turbulence_rit.h"

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * source_time_stepping -- exact exponential integration of the frozen-tau
 * Rotta-Monin/buoyancy source subsystem (thesis Chapter 6, Section
 * "Exact solution of the frozen-time-scale source subsystem").
 *
 * Restricted to pure Rotta closure (crij2 == 0): the "rapid" pressure-
 * strain contribution is not part of this source step and is left to
 * the standard treatment (identically zero for crij2 == 0 anyway).
 *
 * See the engine definition further below for the actual resonance
 * handling (built from a single J(dt;rate) primitive, robust through
 * every resonance without a separate branch for each one).
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Overflow-safe exp(): clamps its argument to [-700, 700] before
 *        calling the standard exp().
 *
 * All the closed-form kernels below assume tau = k/eps > 0 and CR > 1,
 * as required for a physically well-posed Rotta-Monin source step. A
 * transient, purely numerical excursion of R or epsilon outside their
 * realizable range elsewhere in the solve (e.g. a momentarily negative
 * k or eps before the next clipping pass) can make tau non-positive or
 * near-zero, which turns one of the many rate*dt exponents below into
 * a large positive value and makes the plain exp() overflow to inf,
 * raising SIGFPE (observed in practice after ~900 time steps of an
 * otherwise-converged run). Clamping the exponent is a deliberately
 * blunt safeguard -- it does not attempt to fix the underlying
 * transient realizability excursion, only to prevent it from crashing
 * the whole run; exp(700) ~ 1e304 is already far outside any physical
 * source-step magnitude, so the clamp has no effect in the normal,
 * well-posed regime.
 */
/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the turbulent flux source terms
 *
 * \param[in]     name       name of current field
 * \param[in]     f_ut       scalar turbulent flux field
 * \param[in]     f_tv       variance of the thermal scalar field, or nullptr
 * \param[in]     n_cells    number of cells
 * \param[in]     xcpp       \f$ C_p \f$
 * \param[in]     viscl      Molecular viscosity
 * \param[in]     viscls     variable diffusivity field
 * \param[in]     xuta       calculated variables at cell centers
 *                           (at current and previous time steps)
 * \param[in]     gradv      mean velocity gradient
 * \param[in]     gradt      mean scalar gradient
 * \param[in]     grad_al    alpha scalar gradient
 * \param[in]     grav       gravity
 * \param[out]    fimp       implicit part of source term
 * \param[out]    rhs_ut     right-hand side part source term
 */
/*----------------------------------------------------------------------------*/

static void
_turb_flux_st(const char          *name,
              const cs_field_t    *f_ut,
              const cs_field_t    *f_tv,
              const cs_lnum_t      n_cells,
              const cs_real_t      xcpp[],
              const cs_real_t      viscl[],
              const cs_real_t      viscls[],
              const cs_real_33_t   gradv[],
              const cs_real_3_t    gradt[],
              const cs_real_3_t    grad_al[],
              cs_real_33_t         fimp[],
              cs_real_3_t          rhs_ut[])
{
  const cs_real_t *cell_f_vol = cs_glob_mesh_quantities->cell_vol;

  cs_field_t *f = cs_field(name);

  const cs_real_t *crom = CS_F_(rho)->val;
  const cs_real_t *cvar_ep = CS_F_(eps)->val;
  const cs_real_6_t *cvar_rij = (const cs_real_6_t *)CS_F_(rij)->val;

  const cs_real_3_t *xuta = (const cs_real_3_t *)f_ut->val_pre;

  cs_real_t *cpro_beta = nullptr;
  cs_field_t *f_beta = cs_field_try("thermal_expansion");
  if (f_beta != nullptr)
    cpro_beta = f_beta->val;

  /* source_time_stepping: R and epsilon AT THE START OF THE TIME STEP
   * (val_pre), consistent with q_theta_pre/theta2_pre (xuta/cvara_tt
   * above) -- NOT cvar_rij/cvar_ep (->val) just above, which already
   * hold this time step's updated values (R and epsilon are solved
   * before q_theta in the standard segregated order). The exact
   * source-step integration requires all four quantities to start
   * from the SAME state N. */
  const int st_scheme = cs_glob_turb_rans_model->source_time_stepping;
  const cs_real_t ce2   = cs_turb_ce2;
  const cs_real_t ce3   = cs_turb_ce3;
  const int rij_scheme = cs_glob_turb_rans_model->rij_discretization_scheme;
  const bool source_time_stepping_active =
    (st_scheme != CS_TURB_RIJ_SOURCE_TS_IMEX);
  const cs_real_6_t *cvara_rij = source_time_stepping_active ?
    (const cs_real_6_t *)CS_F_(rij)->val_pre : nullptr;
  const cs_real_t *cvara_ep = source_time_stepping_active ?
    CS_F_(eps)->val_pre : nullptr;
  const cs_real_t *dt = source_time_stepping_active ?
    CS_F_(dt)->val : nullptr;

  if (source_time_stepping_active) {
    assert(f_tv != nullptr);
  }

  const int krvarfl = cs_field_key_id("variance_dissipation");
  const cs_real_t ctheta = source_time_stepping_active ?
    1. / f_tv->get_key_double(krvarfl) : 1.;

  /* source_time_stepping: C_theta is derived exclusively from rvarfl (the
   * SAME relation used for theta2's own dissipation rate just above,
   * epsilon_theta = epsilon/(k*rvarfl)*theta2 => tau_theta = rvarfl*tau
   * => C_theta = tau/tau_theta = 1/rvarfl, exact for the homogeneous
   * case alpha_theta == 1, i.e. non-EBRSM).
   *
   * IMPORTANT: this deliberately does NOT involve c1trit (the model
   * constant otherwise used for q_theta's own relaxation rate,
   * c1trit/xttdrbt below). The exact source-step trajectory computed
   * by cs_turbulence_rit_source_step_frozen_tau internally uses the
   * document's own consistent relaxation rate (CR+Ctheta)/(2*tau),
   * derived from crij1 and this same Ctheta -- it does NOT blend with
   * or partially reuse c1trit. c1trit is therefore left UNUSED for
   * the relaxation term when source_time_stepping != 0 (the whole phiit_relax
   * term is replaced, not adjusted): whatever value c1trit is set to
   * in the case setup has no effect on the source step in that
   * configuration. If c1trit differs substantially from
   * (crij1+1/rvarfl)/2, this reflects a genuine modeling choice
   * difference between code_saturne's general DFM closure and the
   * pure Rotta-Monin theory source_time_stepping implements -- not a bug to be
   * silently reconciled. */

  const cs_real_t *cvar_tt = nullptr, *cvara_tt = nullptr, *cvar_al = nullptr;

  const cs_turb_rans_model_t *rans_mdl = cs_glob_turb_rans_model;
  const cs_turb_model_type_t model
    = (cs_turb_model_type_t)cs_glob_turb_model->model;

  /* Get the turbulent flux model */
  int turb_flux_model = f->get_key_int("turbulent_flux_model");

  if (f_tv != nullptr) {
    cvar_tt = f_tv->val;
    cvara_tt = f_tv->val_pre;
  }

  /* Save production terms if required */

  cs_real_3_t *prod_ut = nullptr;
  cs_field_t *f_ut_prod = cs_field_by_double_composite_name_try
                            ("algo:", f->name, "_turbulent_flux_production");

  if (f_ut_prod != nullptr)
    prod_ut = (cs_real_3_t *)f_ut_prod->val;

  cs_real_3_t *phi_ut = nullptr;
  cs_field_t *f_phi_ut = cs_field_by_double_composite_name_try
                           ("algo:", f->name, "_turbulent_flux_scrambling");
  if (f_phi_ut != nullptr)
    phi_ut = (cs_real_3_t *)f_phi_ut->val;

  cs_real_3_t *prod_by_vel_grad_ut = nullptr;
  cs_field_t *f_ut_prod_by_vel
    = cs_field_by_double_composite_name_try
        ("algo:", f->name, "_turbulent_flux_production_by_velocity_gradient");
  if (f_ut_prod_by_vel != nullptr)
    prod_by_vel_grad_ut = (cs_real_3_t *)f_ut_prod_by_vel->val;

  cs_real_3_t *prod_by_scal_grad_ut = nullptr;
  cs_field_t *f_ut_prod_by_scal
    = cs_field_by_double_composite_name_try
        ("algo:", f->name, "_turbulent_flux_production_by_scalar_gradient");
  if (f_ut_prod_by_scal != nullptr)
    prod_by_scal_grad_ut = (cs_real_3_t *)f_ut_prod_by_scal->val;

  cs_real_3_t *buo_ut = nullptr;
  cs_field_t *f_buo_ut = cs_field_by_double_composite_name_try
                           ("algo:", f->name, "_turbulent_flux_buoyancy");
  if (f_buo_ut != nullptr)
    buo_ut = (cs_real_3_t *)f_buo_ut->val;

  cs_real_3_t *dissip_ut = nullptr;
  cs_field_t *f_dissip_ut = cs_field_by_double_composite_name_try
                              ("algo:", f_ut->name, "_dissipation");
  if (f_dissip_ut != nullptr)
    dissip_ut = (cs_real_3_t *)f_dissip_ut->val;

  if (turb_flux_model == 31)
    cvar_al = cs_field_by_composite_name_try(f->name, "alpha")->val;

  const cs_real_t rhebdfm = 0.5;
  const cs_real_t grav[3] = {cs_glob_physical_constants->gravity[0],
                             cs_glob_physical_constants->gravity[1],
                             cs_glob_physical_constants->gravity[2]};

  cs_field_t * f_beta2 = cs_field_try("algo:rij_beta2");
  cs_real_t * v_beta2 = nullptr;
  if (f_beta2 != nullptr)
    v_beta2 = f_beta2->val;

  const cs_real_t c1trit = cs_turb_c1trit;
  const cs_real_t crij1  = cs_turb_crij1;
  const cs_real_t c2trit = cs_turb_c2trit;
  const cs_real_t c3trit = cs_turb_c3trit;
  const cs_real_t c4trit = cs_turb_c4trit;

  const int has_buoyant_term = rans_mdl->has_buoyant_term;

  cs_real_t _visls_0 = -1;
  if (viscls == nullptr)
    _visls_0 = f->get_key_double("diffusivity_ref");

  cs_dispatch_context ctx;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    cs_real_t xrij[3][3];
    xrij[0][0] = cvar_rij[c_id][0];
    xrij[0][1] = cvar_rij[c_id][3];
    xrij[0][2] = cvar_rij[c_id][5];
    xrij[1][0] = cvar_rij[c_id][3];
    xrij[1][1] = cvar_rij[c_id][1];
    xrij[1][2] = cvar_rij[c_id][4];
    xrij[2][0] = cvar_rij[c_id][5];
    xrij[2][1] = cvar_rij[c_id][4];
    xrij[2][2] = cvar_rij[c_id][2];

    /* source_time_stepping: exact frozen-tau source-step integration,
     * computed once per cell (all four quantities together, from the
     * SAME state-N inputs cvara_rij/xuta/cvara_tt/cvara_ep), and
     * reused below for phiit_relax (relaxation+buoyancy combined --
     * buoyancy_i is left at 0 for source_time_stepping != 0, its
     * contribution is already folded in here). */
    cs_real_3_t qtheta1 = {0., 0., 0.};
    if (source_time_stepping_active && cvara_tt != nullptr) {
      const cs_real_t beta_c = (cpro_beta != nullptr) ? cpro_beta[c_id] : 0.;
      cs_real_6_t r1;
      cs_real_t theta2_1, eps1;
      if (st_scheme == CS_TURB_RIJ_SOURCE_TS_VAR_TAU)
        cs_turbulence_rit_source_step_variable_tau(
          crij1, ctheta, ce2, ce3, beta_c,
          grav, dt[c_id],
          cvara_rij[c_id], xuta[c_id], cvara_tt[c_id], cvara_ep[c_id],
          r1, qtheta1, &theta2_1, &eps1);
      else
        cs_turbulence_rit_source_step_frozen_tau(
          crij1, ctheta, ce2, beta_c, grav, dt[c_id],
          cvara_rij[c_id], xuta[c_id], cvara_tt[c_id], cvara_ep[c_id],
          r1, qtheta1, &theta2_1, &eps1);
    }

    cs_real_t prdtl = viscl[c_id]*xcpp[c_id];
    if (viscls != nullptr)
      prdtl /= viscls[c_id];
    else
      prdtl /= _visls_0;

    const cs_real_t tke = 0.5 * cs_math_6_trace(cvar_rij[c_id]);

    /* Compute Durbin time scheme */
    const cs_real_t xttke = tke / cvar_ep[c_id];

    cs_real_t alpha = 1., xttdrbt = xttke, xttdrbw = xttke;
    cs_real_t xxc1 = 0, xxc2 = 0, xxc3 = 0;
    cs_real_t xnal[3] = {0, 0, 0};

    /* EB-DFM */
    if (turb_flux_model == 31) {
      alpha = cvar_al[c_id];
      /* FIXME Warning / rhebdfm**0.5 compared to F Dehoux
       * And so multiplied by (R/Prandt)^0.5 */
      xttdrbt = xttke * sqrt((1.0-alpha)*prdtl/rhebdfm + alpha);
      xttdrbw = xttdrbt * sqrt(rhebdfm/prdtl);

      /* Compute the unit normal vector */
      cs_real_t xnoral = cs_math_3_norm(grad_al[c_id]);
      const cs_real_t eps = cs_math_epzero/cbrt(cell_f_vol[c_id]);
      if (xnoral > eps) {
        for (cs_lnum_t i = 0; i < 3; i++)
          xnal[i] = grad_al[c_id][i] / xnoral;
      }

      /* Production and buoyancy for TKE */
      cs_real_t pk = 0;
      for (cs_lnum_t i = 0; i < 3; i++) {
        for (cs_lnum_t j = 0; j < 3; j++)
          pk -= xrij[i][j]*gradv[c_id][i][j];
      }

      cs_real_t gk = 0;
      /* FIXME make buoyant term coherent elsewhere */
      if (cpro_beta != nullptr && has_buoyant_term == 1)
        gk = - cpro_beta[c_id] * cs_math_3_dot_product(xuta[c_id], grav);

      xxc1 = 1.+2.*(1.-cvar_al[c_id])*(pk+gk)/cvar_ep[c_id];
      xxc2 = 0.5*(1.+1./prdtl)*(1.-0.3*(1.-cvar_al[c_id])
                                      *(pk+gk)/cvar_ep[c_id]);
      xxc3 = xxc2;
    }

    cs_real_t phiith[3], phiitw[3];
    cs_real_t phiit[3];


    /* Pressure/thermal fluctuation correlation term
     * --------------------------------------------- */

    /* Dynamic model must impose dynamic part for the thermal model
     * See BFH */
    if (model == CS_TURB_RIJ_EPSILON_BFH) {

      cs_real_t beta2 = c2trit;
      /* Production of TKE */
      cs_real_t pk = 0;
      for (cs_lnum_t i = 0; i < 3; i++) {
        for (cs_lnum_t j = 0; j < 3; j++)
          pk -= xrij[i][j]*gradv[c_id][i][j];
      }

      if (v_beta2 != nullptr)
        beta2 = v_beta2[c_id];

      for (cs_lnum_t i = 0; i < 3; i++) {
        /* Pope 1994:
         * - 0.5 (C_theta / Tt - 2 alpha1 / Td) * T'u'i
         *   with Tt the thermal time scale and Td the dynamique time scale
         *
         *   -2 alpha1 correspond to C1 (Rotta constant)
         *
         *   Rapid term writes:
         *   beta2 (gradu + gradu^T - Pk / k Id)_ij T'u'j
         *
         *   Note that the last term is added to term proportional to u'T'
         *   (in factor)
         * */
        cs_real_t factor = 0.5 * (c1trit / xttdrbt + crij1 /xttke)
                         - beta2 *pk/tke;

        phiit[i] = - factor * xuta[c_id][i];
        for (cs_lnum_t j = 0; j < 3; j++)
          phiit[i] += beta2
            * (gradv[c_id][i][j]+gradv[c_id][j][i])* xuta[c_id][j];

         if ((cvar_tt != nullptr) && (cpro_beta != nullptr)
             && has_buoyant_term == 1)
           phiit[i] += c3trit*(cpro_beta[c_id] * grav[i] * cvar_tt[c_id]);

         if (f_phi_ut != nullptr) /* Save it if needed */
           phi_ut[c_id][i] = phiit[i];

         cs_real_t imp_term = cell_f_vol[c_id] * crom[c_id] * factor;

         fimp[c_id][i][i] += cs::max(imp_term, 0);

      }

    }
    /* Phi_T for other models */
    else {
      for (cs_lnum_t i = 0; i < 3; i++) {
        /* source_time_stepping: the Monin relaxation piece
         * (-c1trit/xttdrbt * xuta[i], matching the document's
         * -(CR+Ctheta)/(2tau) * theta'u' term) is isolated here so it
         * can be neutralized under source_time_stepping != 0 (handled instead
         * by the exact/quadrature source-step integration in this file)
         * without touching the rapid-redistribution pieces below
         * (c2trit, c4trit), which are a separate closure mechanism,
         * unrelated to the frozen/variable-tau source step, and are
         * therefore left untouched regardless of source_time_stepping. */
        const cs_real_t phiit_relax =
          (st_scheme == CS_TURB_RIJ_SOURCE_TS_IMEX) ?
          - c1trit / xttdrbt * xuta[c_id][i] :
          (qtheta1[i] - xuta[c_id][i]) / dt[c_id];

        phiith[i] = phiit_relax
                    + c2trit * cs_math_3_dot_product( gradv[c_id][i], xuta[c_id])
                    + c4trit * (-xrij[0][i] * gradt[c_id][0]
                                -xrij[1][i] * gradt[c_id][1]
                                -xrij[2][i] * gradt[c_id][2]);

         if ((cvar_tt != nullptr) && (cpro_beta != nullptr)
             && has_buoyant_term == 1)
           phiith[i] += c3trit*(cpro_beta[c_id] * grav[i] * cvar_tt[c_id]);

         phiitw[i] =   -1. / xttdrbw *xxc1   /* FIXME full implicit */
                     * (  xuta[c_id][0]*xnal[0]*xnal[i]
                        + xuta[c_id][1]*xnal[1]*xnal[i]
                        + xuta[c_id][2]*xnal[2]*xnal[i]);

         phiit[i] = alpha * phiith[i] + (1.-alpha) * phiitw[i];
         if (f_phi_ut != nullptr) /* Save it if needed */
           phi_ut[c_id][i] = phiit[i];

         /* source_time_stepping: drop the c1trit/xttdrbt implicit
          * stabilization when the corresponding explicit relaxation
          * term (phiit_relax above) has itself been dropped. */
         const cs_real_t c1_impl_term =
           (st_scheme == 0) ? c1trit/xttdrbt : 0.;

         cs_real_t imp_term
           =   cell_f_vol[c_id] * crom[c_id]
             * (      alpha  * (c1_impl_term - c2trit*gradv[c_id][i][i])
                 // TODO All the following matrix can be implicit
                + (1.-alpha) * (xxc1*xnal[i]*xnal[i]/xttdrbw));

         fimp[c_id][i][i] += cs::max(imp_term, 0);

      }
    }

    for (cs_lnum_t i = 0; i < 3; i++) {
      /* Production terms
       *----------------- */

      /* Production term due to the mean velocity */
      const cs_real_t prod_by_vel_grad_i =
        - cs_math_3_dot_product(gradv[c_id][i], xuta[c_id]);
      if (prod_by_vel_grad_ut != nullptr) /* Save it if needed */
        prod_by_vel_grad_ut[c_id][i] = prod_by_vel_grad_i;

      /* Production term due to the mean temperature */
      const cs_real_t prod_by_scal_grad_i =  - (   xrij[i][0]*gradt[c_id][0]
                                                 + xrij[i][1]*gradt[c_id][1]
                                                 + xrij[i][2]*gradt[c_id][2]);
      if (prod_by_scal_grad_ut != nullptr) /* Save it if needed */
        prod_by_scal_grad_ut[c_id][i] = prod_by_scal_grad_i;

      /* Production term due to the gravity */
      /* source_time_stepping: this is the primary buoyancy mechanism for
       * q_theta (matches the document's -beta_theta*theta2*g term
       * exactly); handled by the source-step integration instead when
       * source_time_stepping != 0. */
      cs_real_t buoyancy_i = 0.;
      if ((cvar_tt != nullptr) && (cpro_beta != nullptr)
          && has_buoyant_term == 1
          && st_scheme == 0)
        buoyancy_i = -grav[i] * cpro_beta[c_id] * cvara_tt[c_id];

      if (buo_ut != nullptr) /* Save it if needed */
        buo_ut[c_id][i] = buoyancy_i;

      /* Dissipation (Wall term only because "h" term is zero */
      const cs_real_t dissip_i =  (1.-alpha)/xttdrbw
                                * (  xxc2 * xuta[c_id][i]
                                   + xxc3 * (  xuta[c_id][0]*xnal[0]*xnal[i]
                                             + xuta[c_id][1]*xnal[1]*xnal[i]
                                             + xuta[c_id][2]*xnal[2]*xnal[i]));
      if (dissip_ut != nullptr)/* Save it if needed */
        dissip_ut[c_id][i] = dissip_i;

      /* Save production terms for post-processing */
      if (prod_ut != nullptr)
        prod_ut[c_id][i] = prod_by_vel_grad_i + prod_by_scal_grad_i
                         + buoyancy_i - dissip_i;

      /* GODUNOV scheme: for CS_RIJ_SCHEME_GODUNOV, the mechanical
       * production terms (by mean velocity gradient and by mean
       * temperature gradient) are already captured by the exact
       * Riemann interface state, added explicitly as cross terms in
       * _solve_rit's divqtheta assembly (built from i_velocity,
       * i_reynolds_stress, i_temperature, i_turbulent_heat_flux).
       * Adding the standard gradient-based values here as well would
       * double-count them, exactly as pij for R. Diagnostics
       * (prod_by_vel_grad_ut, prod_by_scal_grad_ut, prod_ut just
       * above) are left showing the standard gradient-based value for
       * reference; only the contribution actually added to rhs_ut is
       * zeroed here. */
      const cs_real_t mech_prod_vel =
        (rij_scheme == CS_RIJ_SCHEME_GODUNOV) ? 0. : prod_by_vel_grad_i;
      const cs_real_t mech_prod_scal =
        (rij_scheme == CS_RIJ_SCHEME_GODUNOV) ? 0. : prod_by_scal_grad_i;

      rhs_ut[c_id][i] += (  mech_prod_vel + mech_prod_scal
                          + buoyancy_i + phiit[i] - dissip_i)
                        * cell_f_vol[c_id]*crom[c_id];

      /* TODO we can implicit more terms */
      cs_real_t imp_term =   cell_f_vol[c_id] * crom[c_id]
                 * (1.-alpha)/xttdrbw * (xxc2+xxc3*xnal[i]*xnal[i]);

      fimp[c_id][i][i] += cs::max(imp_term, 0);

      if ((cvar_tt != nullptr) && (cpro_beta != nullptr)
          && has_buoyant_term == 1) {

        /* Stable if negative w'T' */
        cs_real_t mez[3];
        cs_math_3_normalize(grav, mez);
        cs_real_t wptp = -cs_math_3_dot_product(mez, xuta[c_id]);
        cs_real_t w2 = cs_math_3_sym_33_3_dot_product(mez,
                                                      cvar_rij[c_id],
                                                      mez);

        if (wptp < - cs_math_epzero * sqrt(cvara_tt[c_id] * w2)) {

          /* Note Cauchy Schwarz implies that
           * T'2/|w'T'| > |w'T'| / w'2
           * */
          imp_term =   cell_f_vol[c_id] * crom[c_id]
            * grav[i] * cpro_beta[c_id] * cvara_tt[c_id] / wptp;

          fimp[c_id][i][i] += cs::max(imp_term, 0);
        }
      }
    }
  });

  ctx.wait();

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief compute the thermal fluxes and Diffusivity.
 *
 * \param[in]     f                  Current field
 * \param[in]     f_tv               variance of the thermal scalar field,
 *                                   or nullptr
 * \param[in]     n_cells            number of cells
 * \param[in]     n_b_faces          number od boundary faces
 * \param[in]     n_cells_ext        number of cells + gost
 * \param[in]     turb_flux_model    turb_flux_model
 * \param[in]     xcpp               \f$ C_p \f$
 * \param[in]     gradv              mean velocity gradient
 * \param[in]     gradt              mean scalar gradient
 * \param[in]     grad_al            alpha scalar gradient
 * \param[out]    xut                calculated variables at cell centers
 *                                    (at current and previous time steps)
 * \param[out]    thflxf             thermal flux on interior faces
 * \param[out]    thflxb             thermal flux on boundary faces
 * \param[out]    vistet             Diffusivity tensor
 */
/*----------------------------------------------------------------------------*/

static void
_thermal_flux_and_diff(cs_field_t         *f,
                       const cs_field_t   *f_tv,
                       cs_lnum_t           n_cells,
                       cs_lnum_t           n_cells_ext,
                       cs_lnum_t           n_b_faces,
                       int                 turb_flux_model,
                       const cs_real_t     xcpp[],
                       const cs_real_33_t  gradv[],
                       const cs_real_3_t   gradt[],
                       const cs_real_3_t   grad_al[],
                       cs_real_3_t         xut[],
                       cs_real_t           thflxf[],
                       cs_real_t           thflxb[],
                       cs_real_6_t         vistet[])
{
  const cs_real_t *cell_f_vol = cs_glob_mesh_quantities->cell_vol;

  const cs_real_t *crom = CS_F_(rho)->val;
  const cs_real_t *viscl  = CS_F_(mu)->val;
  const cs_real_t *brom = CS_F_(rho_b)->val;

  const cs_real_t *cvara_ep = CS_F_(eps)->val_pre;
  const cs_real_6_t *cvara_rij = (const cs_real_6_t *)CS_F_(rij)->val_pre;

  const cs_field_t *f_beta = cs_field_try("thermal_expansion");
  const cs_turb_rans_model_t *rans_mdl = cs_glob_turb_rans_model;
  const cs_real_t *cpro_beta = nullptr, *cvara_tt = nullptr;
  if (f_beta != nullptr)
    cpro_beta = f_beta->val;

  if (f_tv != nullptr)
    cvara_tt = f_tv->val_pre;

  cs_real_t _visls_0 = -1;
  const cs_real_t *viscls = nullptr;
  {
    int ifcvsl = f->get_key_int("diffusivity_id");
    if (ifcvsl > -1) {
      viscls = cs_field(ifcvsl)->val;
    }
    else {
      _visls_0 = f->get_key_double("diffusivity_ref");
    }
  }

  cs_real_t *cvar_al = nullptr;
  if (   (turb_flux_model == 11)
      || (turb_flux_model == 21)
      || (turb_flux_model == 31))
    cvar_al = cs_field_by_composite_name(f->name, "alpha")->val;

  const cs_real_t *grav = cs_glob_physical_constants->gravity;

  const cs_real_t ctheta_c = f->get_key_double("turbulent_flux_ctheta");

  const int has_buoyant_term = rans_mdl->has_buoyant_term;

  const cs_real_t etaafm = cs_turb_etaafm;
  const cs_real_t xiafm = cs_turb_xiafm;

  cs_array_2d<cs_real_t> w1(n_cells_ext, 3, cs_alloc_mode);

  /* loop on cells */

  cs_dispatch_context ctx;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    cs_real_t xnal[3] = {0, 0, 0}, temp[3] = {0, 0, 0};
    cs_real_t ctheta = ctheta_c;

    /* Rij (copy to local 3x3 tensor to allow loops */

    cs_real_t xrij[3][3];
    xrij[0][0] = cvara_rij[c_id][0];
    xrij[0][1] = cvara_rij[c_id][3];
    xrij[0][2] = cvara_rij[c_id][5];
    xrij[1][0] = cvara_rij[c_id][3];
    xrij[1][1] = cvara_rij[c_id][1];
    xrij[1][2] = cvara_rij[c_id][4];
    xrij[2][0] = cvara_rij[c_id][5];
    xrij[2][1] = cvara_rij[c_id][4];
    xrij[2][2] = cvara_rij[c_id][2];

    /* Epsilon */
    const cs_real_t xe = cvara_ep[c_id];

    /* Kinetic turbulent energy */
    const cs_real_t xk = 0.5 * cs_math_6_trace(cvara_rij[c_id]);

    /* Turbulent time-scale (constant in AFM) */
    const cs_real_t xtt = xk/xe;
    cs_real_t alpha_theta = 0, xpk = 0., xgk = 0;
    cs_real_t eta_ebafm = 0, xi_ebafm = 0, gamma_eb = 0;

    if ((turb_flux_model == 11) || (turb_flux_model == 21)) {

      alpha_theta = cvar_al[c_id];

      /* Production and buoyancy */
      xpk = 0;
      for (cs_lnum_t ii = 0; ii < 3; ii++) {
        for (cs_lnum_t jj = 0; jj < 3; jj++)
          xpk -= xrij[jj][ii] * gradv[c_id][jj][ii];
      }
      if (cpro_beta != nullptr)
        xgk = -cpro_beta[c_id] * cs_math_3_dot_product(xut[c_id], grav);

      /* Thermo-mecanical scales ratio R */
      cs_real_t prdtl = viscl[c_id] * xcpp[c_id];
      if (viscls != nullptr)
        prdtl /= viscls[c_id];
      else
        prdtl /= _visls_0;

      cs_real_t xr_h = 0.5;
      cs_real_t xr = (1.0-alpha_theta)*prdtl + alpha_theta*xr_h;

      /* Unit normal vector */
      cs_real_t xnoral = cs_math_3_norm(grad_al[c_id]);

      const cs_real_t eps = cs_math_epzero/cbrt(cell_f_vol[c_id]);
      if (xnoral > eps) {
        for (cs_lnum_t i = 0; i < 3; i++)
          xnal[i] = grad_al[c_id][i] / xnoral;
      }

      /* Constants for EB-GGDH and EB-AFM */

      cs_real_t xxc1 = 1. + 2.*(1.-alpha_theta)*(xpk+xgk)/cvara_ep[c_id];
      cs_real_t xxc2 = 0.5*(1.+1./prdtl)*(1. -0.3*(1.-alpha_theta)
                                                 *(xpk+xgk)/cvara_ep[c_id]);

      ctheta =   (0.97*sqrt(xr))/(alpha_theta*(4.15*sqrt(0.5))
               + (1.-alpha_theta)*(sqrt(prdtl))*xxc2);
      gamma_eb = (1.-alpha_theta)*(xxc1 + xxc2);

      /* Constants for EB-AFM */
      if (turb_flux_model == 21) {
        eta_ebafm = 1.0 - alpha_theta*0.6;
        xi_ebafm  = 1.0 - alpha_theta*0.3;
      }

    }

    /* Compute thermal flux u'T' */

    for (cs_lnum_t ii = 0; ii < 3; ii++) {
      temp[ii] = 0;

      /* AFM model
         "-C_theta*k/eps*( xi* uT'.Grad u + eta*beta*g_i*T'^2)" */
      if (turb_flux_model == 20) {
        if ((cvara_tt != nullptr) && (cpro_beta != nullptr)
            && has_buoyant_term == 1)
          temp[ii] -=   ctheta * xtt * etaafm
                      * cpro_beta[c_id] * grav[ii] * cvara_tt[c_id];

        for (cs_lnum_t jj = 0; jj < 3; jj++) {
          /*  Partial implicitation of "-C_theta*k/eps*(xi* uT'.Grad u)"
           *  Only the i != j  components are added. */
          if (ii != jj)
            temp[ii] -=  ctheta*xtt*xiafm
                        *xut[c_id][jj]*gradv[c_id][ii][jj];
          else
            temp[ii] -= cs::min(  ctheta*xtt*xiafm*xut[c_id][jj]
                                * gradv[c_id][ii][jj],
                                0.);
        }
      }

      /* EB-AFM model
       *  "-C_theta*k/eps*(  xi*uT'.Gradu+eta*beta*g_i*T'^2
       *                   + eps/k gamma uT' ni nj )"
       */
      if (turb_flux_model == 21) {
        if ((cvara_tt != nullptr) && (cpro_beta != nullptr)
            && has_buoyant_term == 1)
          temp[ii] -=   ctheta * xtt * eta_ebafm
                      * cpro_beta[c_id] * grav[ii] * cvara_tt[c_id];
        for (cs_lnum_t jj = 0; jj < 3; jj++) {
          /* Partial implicitation of
           * "-C_theta*k/eps*( xi* uT'.Grad u + eps/k gamma uT' ni nj)"
           * Only the i.ne.j  components are added. */
          cs_real_t tmp1 = xtt * xi_ebafm * gradv[c_id][ii][jj] * xut[c_id][jj];
          if (ii != jj)
            temp[ii] -=   ctheta * tmp1
                        + ctheta*gamma_eb*xnal[ii]*xnal[jj]*xut[c_id][jj];
          else
            temp[ii] -=   ctheta
                        * cs::min(tmp1 +  gamma_eb*xnal[ii]
                                         *xnal[jj]*xut[c_id][jj],
                                  0);
        }
      }

      /* EB-GGDH model
       *  "-C_theta*k/eps*( eps/k gamma uT' ni nj)" */
      if (turb_flux_model == 11) {
        for (cs_lnum_t jj = 0; jj < 3; jj++) {
          /* Partial implicitation of "-C_theta*k/eps*( eps/k gamma uT' ni nj)"
           * Only the i.ne.j  components are added. */
          if (ii != jj)
            temp[ii] -= ctheta*gamma_eb*xnal[ii]*xnal[jj]*xut[c_id][jj];
        }
      }
    }

    cs_real_t coeff_imp;

    for (cs_lnum_t ii = 0; ii < 3; ii++) {
      /* Add the term in "grad T" which is implicited by the GGDH part in
         cs_solve_equation_scalar.
       *  "-C_theta*k/eps* R.grad T"
       * The resulting XUT array is only use for post processing purpose in
       * (EB)GGDH & (EB)AFM */
      xut[c_id][ii] =   temp[ii]
                      - ctheta*xtt*(  xrij[0][ii]*gradt[c_id][0]
                                    + xrij[1][ii]*gradt[c_id][1]
                                    + xrij[2][ii]*gradt[c_id][2]);
      /* Partial implicitation of "-C_theta*k/eps*( xi* uT'.Grad u )" for
       * EB-GGDH & (EB)-AFM
       * if positive
       * X_i = C*Y_ij*X_j -> X_i = Coeff_imp * Y_ij * X_j for i.ne.j
       * with Coeff_imp = C/(1+C*Y_ii) */

      /* AFM */
      if (turb_flux_model == 20) {
        coeff_imp
          = 1.+cs::max(ctheta*xtt*xiafm*gradv[c_id][ii][ii], 0);

        xut[c_id][ii] = xut[c_id][ii]/coeff_imp;
        temp[ii] = temp[ii]/coeff_imp;
        /* Calculation of the diffusion tensor for the implicited part
         * of the model computed in cs_convection_diffusion_solve.c */
        vistet[c_id][ii] = crom[c_id]*ctheta*xtt*xrij[ii][ii]/coeff_imp;

      }

      /* EB-AFM */
      else if (turb_flux_model == 21) {
        coeff_imp = 1. + cs::max(  ctheta*xtt*xi_ebafm*gradv[c_id][ii][ii]
                                 + ctheta*gamma_eb*xnal[ii]*xnal[ii],
                                 0.);

        xut[c_id][ii] = xut[c_id][ii]/coeff_imp;
        temp[ii] = temp[ii] / coeff_imp;
        /* Calculation of the diffusion tensor for the implicited part
         * of the model computed in cs_convection_diffusion_solve.c */
        vistet[c_id][ii] = crom[c_id]*ctheta*xtt*xrij[ii][ii]/coeff_imp;

      }

      /*! EB-GGDH */
      else if (turb_flux_model == 11) {
        coeff_imp = 1. + ctheta*gamma_eb*xnal[ii]*xnal[ii];

        xut[c_id][ii] = xut[c_id][ii]/coeff_imp;
        temp[ii] = temp[ii]/coeff_imp;
        /* Calculation of the diffusion tensor for the implicited part
         * of the model computed in cs_convection_diffusion_solve.c */
        vistet[c_id][ii] = crom[c_id]*ctheta*xtt*xrij[ii][ii]/coeff_imp;
      }

      /* In the next step, we compute the divergence of:
       *  "-Cp*C_theta*k/eps*( xi* uT'.Grad u + eta*beta*g_i*T'^2)"
       *  The part "-C_theta*k/eps* R.Grad T" is computed by the GGDH part */
      w1(c_id, ii) = xcpp[c_id]*temp[ii];
    }

    /*  Extra diag part of the diffusion tensor
        for cs_convection_diffusion_solve.c */
    if (   (turb_flux_model == 11)
        || (turb_flux_model == 20)
        || (turb_flux_model == 21)) {
      vistet[c_id][3] = crom[c_id]*ctheta*xtt*xrij[1][0];
      vistet[c_id][4] = crom[c_id]*ctheta*xtt*xrij[2][1];
      vistet[c_id][5] = crom[c_id]*ctheta*xtt*xrij[2][0];
    }

  }); /* End loop on cells */
  ctx.wait();

  cs_solid_zone_set_zero_on_cells(3, (cs_real_t *)xut);

  /* FIXME the line below would reproduce the previous behavior, which
     cannot be correct (see issue #387). Either we should consider
     ctheta here purely local, or we must use an associated field to save it. */

  /* f->set_key_double(kctheta, ctheta); */

  cs_field_bc_coeffs_t bc_coeffs_v_loc(3);
  CS_MALLOC_HD(bc_coeffs_v_loc.a, 3*n_b_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(bc_coeffs_v_loc.b, 9*n_b_faces, cs_real_t, cs_alloc_mode);

  cs_real_3_t  *coefat = (cs_real_3_t  *)bc_coeffs_v_loc.a;
  cs_real_33_t *coefbt = (cs_real_33_t *)bc_coeffs_v_loc.b;

  const cs_real_t kr_33[3][3] = {{1., 0., 0.},
                                 {0., 1., 0.},
                                 {0., 0., 1.}};

  ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t f_id) {
    for (cs_lnum_t ii = 0; ii < 3; ii++) {
      coefat[f_id][ii] = 0.;
    }
    for (cs_lnum_t ii = 0; ii < 3; ii++) {
      for (cs_lnum_t jj = 0; jj < 3; jj++)
        coefbt[f_id][ii][jj] = kr_33[ii][jj];
    }
  });
  ctx.wait();

  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);;

  cs_mass_flux(cs_glob_mesh,
               cs_glob_mesh_quantities,
               -1,
               1,
               1,
               1,
               1,
               eqp->imrgra,
               eqp->nswrgr,
               (cs_gradient_limit_t)(eqp->imligr),
               eqp->verbosity,
               eqp->epsrgr,
               eqp->climgr,
               crom,
               brom,
               w1.data<cs_real_3_t>(),
               &bc_coeffs_v_loc,
               thflxf,
               thflxb);

  cs_field_bc_coeffs_clear(&bc_coeffs_v_loc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function perform the solving of the transport equation
 * of the turbulent heat fluxes.
 *
 * \param[in]     f             pointer to scalar field
 * \param[in,out] f_ut          pointer to turbulent flux field
 * \param[in]     xcpp          \f$ C_p \f$
 * \param[in]     gradv         mean velocity gradient
 * \param[in]     gradt         mean scalar gradient
 * \param[in]     grad_al       alpha scalar gradient
 */
/*----------------------------------------------------------------------------*/

static void
_solve_rit(const cs_field_t     *f,
           cs_field_t           *f_ut,
           const cs_real_t       xcpp[],
           const cs_real_33_t    gradv[],
           const cs_real_3_t     gradt[],
           const cs_real_3_t     grad_al[])
{
  if (cs_glob_turb_model->order == CS_TURB_FIRST_ORDER)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: use an Rij model with thermal model."),
              __func__);

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_real_t *cell_f_vol = mq->cell_vol;

  const cs_real_t *dt = CS_F_(dt)->val;
  const cs_real_t *crom = CS_F_(rho)->val;
  const cs_real_t *viscl  = CS_F_(mu)->val;
  const cs_real_t *visct = CS_F_(mu_t)->val;
  const cs_real_6_t *visten
    = (const cs_real_6_t *)cs_field
                             ("anisotropic_turbulent_viscosity")->val;

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  const int iflmas = CS_F_(vel)->get_key_int(kimasf);
  const int iflmab = CS_F_(vel)->get_key_int(kbmasf);

  const cs_real_t *imasfl = cs_field(iflmas)->val;
  const cs_real_t *bmasfl = cs_field(iflmab)->val;

  const cs_real_3_t *xuta = (cs_real_3_t *)f_ut->val_pre;
  cs_real_3_t *xut = (cs_real_3_t *)f_ut->val;

  /* vcopt */
  const cs_equation_param_t *eqp
    = cs_field_get_equation_param_const(f);

  /* vcopt_ut */
  const cs_equation_param_t *eqp_ut
    = cs_field_get_equation_param_const(f_ut);

  if (eqp->verbosity >= 1)
    bft_printf(" Solving variable %s\n", f_ut->name);

  int st_prv_id = f_ut->get_key_int("source_term_prev_id");
  cs_real_3_t *c_st_prv = nullptr;
  if (st_prv_id > -1)
    c_st_prv = (cs_real_3_t *)cs_field(st_prv_id)->val;

  const int rij_scheme
    = cs_glob_turb_rans_model->rij_discretization_scheme;

  cs_real_t _visls_0 = -1;
  const cs_real_t *viscls = nullptr;
  {
    int ifcvsl = f->get_key_int("diffusivity_id");
    if (ifcvsl > -1) {
      viscls = cs_field(ifcvsl)->val;
    }
    else {
      _visls_0 = f->get_key_double("diffusivity_ref");
    }
  }

  cs_dispatch_context ctx;

  cs_array_3d<cs_real_t> fimp(n_cells_ext, 3, 3, cs_alloc_mode);
  cs_array_2d<cs_real_t> rhs_ut(n_cells_ext, 3, cs_alloc_mode);

  fimp.zero(ctx);
  rhs_ut.zero(ctx);
  ctx.wait();

  /* Find the corresponding variance of the scalar */

  const cs_real_t *grav = cs_glob_physical_constants->gravity;

  const cs_field_t *f_tv = nullptr;

  if (cs_math_3_norm(grav) > cs_math_epzero)
    f_tv = cs_field_get_variance(f);
  else
    grav = nullptr;

  /* User source terms
     ----------------- */

  cs_user_source_terms(cs_glob_domain,
                       f_ut->id,
                       rhs_ut,
                       fimp);

  const cs_real_t thetv = eqp->theta;

  if (st_prv_id > -1) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t i = 0; i < 3; i++) {
        for (cs_lnum_t j = 0; j < 3; j++) {
          rhs_ut(c_id, i) = fimp(c_id, i, j)*xuta[c_id][j];
          fimp(c_id, i, j) = -thetv*fimp(c_id, i, j);
        }
      }
    });
  }

  /* If we do not extrapolate the source terms */
  else {
    const cs_real_t zero_threshold = cs_math_zero_threshold;
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t i = 0; i < 3; i++) {
        for (cs_lnum_t j = 0; j < 3; j++) {
          /* User source term */
          rhs_ut(c_id, i) += fimp(c_id, i, j)*xuta[c_id][j];
        }
        /* Diagonal */
        fimp(c_id, i, i) = cs::max(-fimp(c_id, i, i),
                                   zero_threshold);
      }
    });
  }

  /* Mass source terms FIXME
   * ----------------------- */

  /* Unsteady term
   * ------------- */

  if (eqp->istat == 1) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t i = 0; i < 3; i++)
        fimp(c_id, i, i) += (crom[c_id] / dt[c_id]) * cell_f_vol[c_id];
    });
  }

  ctx.wait();

  /* Right Hand Side of the turbulent fluxes:
   *     rho*(Pit + Git + Phi*_it - eps_it)
   * -------------------------------------- */

  _turb_flux_st(f->name, f_ut, f_tv, n_cells,
                xcpp, viscl, viscls, gradv,
                gradt, grad_al,
                fimp.data<cs_real_33_t>(),
                rhs_ut.data<cs_real_3_t>());

  /* Tensor diffusion
   * ---------------- */

  cs_array<cs_real_t> w1(n_cells_ext, cs_alloc_mode);
  cs_array<cs_real_t> viscf(n_i_faces, cs_alloc_mode);
  cs_array<cs_real_t> viscb(n_b_faces, cs_alloc_mode);
  cs_array<cs_real_t> weighb(n_b_faces, cs_alloc_mode);
  cs_array_2d<cs_real_t> weighf(n_i_faces, 2, cs_alloc_mode);
  cs_array_2d<cs_real_t> viscce(n_cells_ext, 6, cs_alloc_mode);

  cs_real_t mdifft = (cs_real_t)(eqp_ut->idifft);

  const cs_real_t ctheta = f->get_key_double("turbulent_flux_ctheta");

  /* Symmetric tensor diffusivity (GGDH) */
  if (eqp_ut->idiff > 0) {
    if (eqp_ut->idften & CS_ANISOTROPIC_RIGHT_DIFFUSION) {
      const cs_real_t a = mdifft * ctheta / cs_turb_csrij;

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        cs_real_t prdtl = viscl[c_id]*xcpp[c_id];
        if (viscls != nullptr)
          prdtl /= viscls[c_id];
        else
          prdtl /= _visls_0;

        for (cs_lnum_t i = 0; i < 3; i++)
          viscce(c_id, i) =   0.5*(viscl[c_id]*(1.+1./prdtl))
                            + a*visten[c_id][i];
        for (cs_lnum_t i = 3; i < 6; i++)
          viscce(c_id, i) = a*visten[c_id][i];
      });
      ctx.wait();

      cs_face_anisotropic_viscosity_scalar(m,
                                           mq,
                                           viscce.data<cs_real_6_t>(),
                                           eqp->verbosity,
                                           weighf.data<cs_real_2_t>(),
                                           weighb,
                                           viscf,
                                           viscb);
    }

    /* Scalar diffusivity */
    else {

      cs_real_t cmu = cs_turb_cmu;
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        w1[c_id] = viscl[c_id] + mdifft*(ctheta*visct[c_id]/cmu);
      });
      ctx.wait();

      cs_face_viscosity(m,
                        mq,
                        eqp->imvisf,
                        w1,
                        viscf,
                        viscb);
    }
  }
  /* No diffusion */
  else {
    viscf.zero(ctx);
    viscb.zero(ctx);
    ctx.wait();
  }

  /* Add Rusanov fluxes */
  if (rij_scheme == CS_RIJ_SCHEME_RUSANOV) {
    cs_real_t *ipro_rusanov = cs_field("i_rusanov_diff")->val;
    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      viscf[face_id] = cs::max(viscf[face_id], 0.5 * ipro_rusanov[face_id]);
    });

    const cs_nreal_3_t *restrict b_face_u_normal = mq->b_face_u_normal;
    cs_real_t *bpro_rusanov = cs_field("b_rusanov_diff")->val;

    //cs_real_3_t *coefap = (cs_real_3_t *)f_ut->bc_coeffs->a;
    cs_real_33_t *cofbfp = (cs_real_33_t *)f_ut->bc_coeffs->bf;
    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      const cs_nreal_t *n = b_face_u_normal[face_id];

      for (cs_lnum_t i = 0; i < 3; i++) {
        for (cs_lnum_t j = 0; j < 3; j++) {
          cofbfp[face_id][i][j] +=  bpro_rusanov[face_id] * n[i]*n[j];
          //TODO ?cofafp[face_id][i] -= bf[i][j] * coefap[face_id][j];
        }
      }
    });

    ctx.wait();
  }

  /* Vectorial solving of the turbulent thermal fluxes
   * ------------------------------------------------- */

  if (st_prv_id > -1) {
    const cs_time_scheme_t *time_scheme = cs_glob_time_scheme;
    const cs_real_t thets = time_scheme->thetst;
    const cs_real_t thetp1 = 1.0+thets;
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t i = 0; i < 3; i++)
        rhs_ut(c_id, i) += thetp1*c_st_prv[c_id][i];
    });
    ctx.wait();
  }

  cs_equation_param_t eqp_loc = *eqp;
  eqp_loc.iwgrec = 0;     /* Warning, may be overwritten if a field */
  eqp_loc.theta = thetv;
  eqp_loc.blend_st = 0;   /* Warning, may be overwritten if a field */

  /* GODUNOV scheme: explicit convection + cross-production
   * of q_theta itself, from the continuous equation
   *   d_t q_theta + u.grad(q_theta) + grad(u).q_theta + R.grad(T) = ...
   * Godunov deferred-correction form, built from the exact Riemann
   * interface state shared with the {u,R} and T/variance equations
   * (i_velocity, i_reynolds_stress, i_temperature,
   * i_turbulent_heat_flux). iconv is disabled below so that the
   * standard convective operator is not assembled a second time on
   * top of this explicit contribution -- same rationale as
   * GODUNOV scheme notes for R in cs_turbulence_rij.cpp. */
  if (rij_scheme == CS_RIJ_SCHEME_GODUNOV) {

    const cs_lnum_2_t *restrict i_face_cells
      = (const cs_lnum_2_t *)m->i_face_cells;
    const cs_lnum_t *restrict b_face_cells
      = (const cs_lnum_t *)m->b_face_cells;
    const cs_real_3_t *restrict i_face_normal_g
      = (const cs_real_3_t *)mq->i_face_normal;
    const cs_real_3_t *restrict b_face_normal_g
      = (const cs_real_3_t *)mq->b_face_normal;

    const cs_real_3_t *c_vel_g = (const cs_real_3_t *)CS_F_(vel)->val;
    const cs_real_t   *c_temp_g = f->val;

    std::string i_name = std::string("i_") + f->name;
    std::string b_name = std::string("b_") + f->name;
    std::string i_tf_name =
      std::string("i_") + f->name + "_turbulent_flux";
    std::string b_tf_name =
      std::string("b_") + f->name + "_turbulent_flux";

    const cs_real_t *i_temp_g =
      cs_field(i_name.c_str())->val;
    const cs_real_t *b_temp_g =
      cs_field(b_name.c_str())->val;
    const cs_real_3_t *i_qtheta_g =
      (const cs_real_3_t *)cs_field(i_tf_name.c_str())->val;
    const cs_real_3_t *b_qtheta_g =
      (const cs_real_3_t *)cs_field(b_tf_name.c_str())->val;
    const cs_real_3_t *i_vel_g =
      (const cs_real_3_t *) cs_field("i_velocity")->val;
    const cs_real_3_t *b_vel_g =
      (const cs_real_3_t *) cs_field("b_velocity")->val;
    const cs_real_6_t *i_rij_g =
      (const cs_real_6_t *) cs_field("i_reynolds_stress")->val;
    const cs_real_6_t *b_rij_g =
      (const cs_real_6_t *) cs_field("b_reynolds_stress")->val;

    cs_real_3_t *divqtheta;
    CS_MALLOC_HD(divqtheta, n_cells_ext, cs_real_3_t, cs_alloc_mode);
    cs_arrays_set_value<cs_real_t, 1>(3*n_cells_ext, 0., (cs_real_t *)divqtheta);

    cs_dispatch_sum_type_t i_sum_type_g =
      ctx.get_parallel_for_i_faces_sum_type(m);
    cs_dispatch_sum_type_t b_sum_type_g =
      ctx.get_parallel_for_b_faces_sum_type(m);

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {

      const cs_lnum_t c_id_l = i_face_cells[face_id][0];
      const cs_lnum_t c_id_r = i_face_cells[face_id][1];

      const cs_real_t qn_f =
        cs_math_3_dot_product(i_qtheta_g[face_id], i_face_normal_g[face_id]);

      cs_real_3_t rin_f;
      cs_math_sym_33_3_product(i_rij_g[face_id], i_face_normal_g[face_id], rin_f);

      const cs_real_t dtheta_l = i_temp_g[face_id] - c_temp_g[c_id_l];
      const cs_real_t dtheta_r = i_temp_g[face_id] - c_temp_g[c_id_r];

      cs_real_3_t flux_q_l, flux_q_r;
      for (cs_lnum_t i = 0; i < 3; i++) {

        flux_q_l[i] =
            (i_qtheta_g[face_id][i] - xut[c_id_l][i]) * imasfl[face_id]
          + rin_f[i] * dtheta_l
          + qn_f * (i_vel_g[face_id][i] - c_vel_g[c_id_l][i]);

        flux_q_r[i] =
          -(  (i_qtheta_g[face_id][i] - xut[c_id_r][i]) * imasfl[face_id]
            + rin_f[i] * dtheta_r
            + qn_f * (i_vel_g[face_id][i] - c_vel_g[c_id_r][i]));
      }

      if (c_id_l < n_cells)
        cs_dispatch_sum<3>(divqtheta[c_id_l], flux_q_l, i_sum_type_g);
      if (c_id_r < n_cells)
        cs_dispatch_sum<3>(divqtheta[c_id_r], flux_q_r, i_sum_type_g);
    });

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {

      const cs_lnum_t c_id_l = b_face_cells[face_id];

      const cs_real_t qn_f =
        cs_math_3_dot_product(b_qtheta_g[face_id], b_face_normal_g[face_id]);

      cs_real_3_t rin_f;
      cs_math_sym_33_3_product(b_rij_g[face_id], b_face_normal_g[face_id], rin_f);

      const cs_real_t dtheta = b_temp_g[face_id] - c_temp_g[c_id_l];

      cs_real_3_t flux_q;
      for (cs_lnum_t i = 0; i < 3; i++) {
        flux_q[i] =
            (b_qtheta_g[face_id][i] - xut[c_id_l][i]) * bmasfl[face_id]
          + rin_f[i] * dtheta
          + qn_f * (b_vel_g[face_id][i] - c_vel_g[c_id_l][i]);
      }

      if (c_id_l < n_cells)
        cs_dispatch_sum<3>(divqtheta[c_id_l], flux_q, b_sum_type_g);
    });

    ctx.wait();

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t i = 0; i < 3; i++)
        rhs_ut(c_id, i) -= divqtheta[c_id][i];
    });
    ctx.wait();

    CS_FREE(divqtheta);

    /* GODUNOV scheme notes: the standard convective operator MUST be
     * disabled, or convection would be assembled twice for q_theta as
     * well -- same double-counting risk as for R and T. */
    eqp_loc.iconv = 0;
  }

  cs_equation_iterative_solve_vector(cs_glob_time_step_options->idtvar,
                                     1, // init
                                     f_ut->id,
                                     nullptr,
                                     0,
                                     0,
                                     &eqp_loc,
                                     xuta,
                                     xuta,
                                     f_ut->bc_coeffs,
                                     imasfl,
                                     bmasfl,
                                     viscf,
                                     viscb,
                                     viscf,
                                     viscb,
                                     nullptr,
                                     nullptr,
                                     viscce.data<cs_real_6_t>(),
                                     weighf.data<cs_real_2_t>(),
                                     weighb,
                                     0,
                                     nullptr,
                                     fimp.data<cs_real_33_t>(),
                                     rhs_ut.data<cs_real_3_t>(),
                                     xut,
                                     nullptr);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the divergence of turbulent flux to a scalar transport equation.
 *
 * \param[in]   field_id  transported field id
 * \param[in]   xcpp      Cp
 * \param[out]  vistet    diffusivity tensor
 * \param[out]  rhs       right hand side to update
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rit_div(const int        field_id,
                      const cs_real_t  xcpp[],
                      cs_real_t        vistet[][6],
                      cs_real_t        rhs[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  /* TODO: declare field as const when ctheta issue (#387) is solved */
  cs_field_t *f = cs_field(field_id);

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const int turb_flux_model = f->get_key_int("turbulent_flux_model");
  const int turb_flux_model_type = turb_flux_model / 10;

  cs_dispatch_context ctx;

  /* Value of the corresponding turbulent flux */
  cs_field_t *f_ut = cs_field_by_composite_name(f->name, "turbulent_flux");
  cs_real_3_t *xut = (cs_real_3_t *)f_ut->val;

  cs_field_t *f_vel = CS_F_(vel);

  /* Compute velocity gradient */

  cs_array_3d<cs_real_t> gradv;
  {
    cs_field_t *f_vg = cs_field_try("algo:velocity_gradient");

    if (f_vel->grad != nullptr)
      gradv = cs_array_3d<cs_real_t>(f_vel->grad,
                                     n_cells_ext, 3, 3);
    else if (f_vg != nullptr)
      gradv = cs_array_3d<cs_real_t>(f_vg->val,
                                     n_cells_ext, 3, 3);
    else {
      gradv = cs_array_3d<cs_real_t>(n_cells_ext, 3, 3, cs_alloc_mode);
    }
  }

  cs_field_gradient_vector(f_vel, false, 1, gradv.data<cs_real_33_t>());

  /* Compute scalar gradient */

  cs_array_2d<cs_real_t> gradt;
  {
    cs_field_t *f_tg = cs_field_by_double_composite_name_try
                         ("algo:", f->name, "_gradient");

    if (f_tg != nullptr)
      gradt = cs_array_2d<cs_real_t>(f_tg->val, n_cells_ext, 3);
    else {
      gradt = cs_array_2d<cs_real_t>(n_cells_ext, 3, cs_alloc_mode);
    }
  }

  cs_field_gradient_scalar(f,
                           true,     /* use previous t   */
                           1,        /* not on increment */
                           gradt.data<cs_real_3_t>());


  /* EB- AFM or EB-DFM: compute the gradient of alpha of the scalar */

  cs_array_2d<cs_real_t> grad_al;

  if (   (turb_flux_model == 11)
      || (turb_flux_model == 21)
      || (turb_flux_model == 31)) {

    grad_al.set_alloc_mode(cs_alloc_mode);
    grad_al.reshape(n_cells_ext, 3);

    cs_field_gradient_scalar(cs_field_by_composite_name(f->name, "alpha"),
                             false,       /* use previous t */
                             1,           /* not on increment */
                             grad_al.data<cs_real_3_t>());
  }

  /* Find the corresponding variance of the scalar */

  const cs_field_t *f_tv = nullptr;

  const int irovar = cs_glob_fluid_properties->irovar;
  const int idilat = cs_glob_velocity_pressure_model->idilat;
  const cs_real_t *grav = cs_glob_physical_constants->gravity;
  const cs_turb_rans_model_t *rans_mdl = cs_glob_turb_rans_model;

  const cs_real_t mod_grav = cs_math_3_norm(grav);
  if (   (mod_grav > cs_math_epzero)
      && ((irovar > 0) || (idilat == 0))
      && ((turb_flux_model_type == 2) || (turb_flux_model_type == 3))
      && rans_mdl->has_buoyant_term == 1) {

    f_tv = cs_field_get_variance(f);

    if (f_tv == nullptr)
      bft_error(__FILE__, __LINE__, 0,
                _("%s: the variance field required for\n"
                  "the turbulent transport of \"%s\" is not available."),
                __func__, f->name);

  }
  else
    grav = nullptr;

  /* Agebraic models AFM
   * ------------------- */

  cs_array<cs_real_t> thflxf(n_i_faces, cs_alloc_mode);
  cs_array<cs_real_t> thflxb(n_b_faces, cs_alloc_mode);

  if (turb_flux_model_type != 3) {

    thflxf.zero(ctx);
    thflxb.zero(ctx);
    ctx.wait();

    _thermal_flux_and_diff(f,
                           f_tv,
                           n_cells,
                           n_cells_ext,
                           n_b_faces,
                           turb_flux_model,
                           xcpp,
                           gradv.data<cs_real_33_t>(),
                           gradt.data<cs_real_3_t>(),
                           grad_al.data<cs_real_3_t>(),
                           xut,
                           thflxf,
                           thflxb,
                           vistet);

  }
  else {

    /* Transport equation on turbulent thermal fluxes (DFM)
     * ---------------------------------------------------- */

    _solve_rit(f, f_ut, xcpp, gradv.data<cs_real_33_t>(),
               gradt.data<cs_real_3_t>(), grad_al.data<cs_real_3_t>());

    /*  Clipping of the turbulence flux vector */
    if ((f_tv != nullptr) && (cs_glob_time_step->nt_cur > 1)) {
      const int clprit = f_ut->get_key_int("is_clipped");
      if (clprit > 0)
        cs_clip_turbulent_fluxes(f_ut->id,
                                 f_tv->id);
    }

    const cs_real_t *crom = CS_F_(rho)->val;
    const cs_real_t *brom = CS_F_(rho_b)->val;

    cs_array_2d<cs_real_t> w1(n_cells_ext, 3, cs_alloc_mode);

    ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t ii = 0; ii < 3; ii ++)
        w1(c_id, ii) = xcpp[c_id] * xut[c_id][ii];
    });
    ctx.wait();

    /* Boundary Conditions on T'u' for the divergence term of
     * the thermal transport equation */

    cs_field_bc_coeffs_t bc_coeffs(1);

    bc_coeffs.a = f_ut->bc_coeffs->ad;
    bc_coeffs.b = f_ut->bc_coeffs->bd;

    const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);

    cs_mass_flux(m,
                 mq,
                 -1, /*f_id */
                 1,
                 1,
                 1,
                 1,
                 eqp->imrgra,
                 eqp->nswrgr,
                 static_cast<cs_gradient_limit_t>(eqp->imligr),
                 eqp->verbosity,
                 eqp->epsrgr,
                 eqp->climgr,
                 crom,
                 brom,
                 w1.data<cs_real_3_t>(),
                 &bc_coeffs,
                 thflxf,
                 thflxb);

    bc_coeffs.a = nullptr;
    bc_coeffs.b = nullptr;
    cs_field_bc_coeffs_clear(&bc_coeffs);

  }

  /* Add the divergence of the thermal flux to the thermal transport equation
     ------------------------------------------------------------------------ */

  if (   turb_flux_model == 11
      || turb_flux_model_type == 2
      || turb_flux_model_type == 3) {

    cs_field_t *f_dut = cs_field_by_double_composite_name_try
                          ("algo:", f_ut->name, "_divergence");

    cs_array<cs_real_t> divut;
    if (f_dut != nullptr)
      divut = cs_array<cs_real_t>(f_dut->val, n_cells_ext);
    else {
      divut = cs_array<cs_real_t>(n_cells_ext, cs_alloc_mode);
    }

    cs_divergence(m, 1, thflxf, thflxb, divut);

    ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      rhs[c_id] -= divut[c_id];
    });

    /* For post-processing intensive quantities */
    if (f_dut != nullptr) {
      int has_disable_flag = mq->has_disable_flag;
      int *c_disable_flag = mq->c_disable_flag;
      const cs_real_t *cell_f_vol = mq->cell_vol;

      ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        cs_real_t dvol = 0;
        const int ind = has_disable_flag * c_id;
        const int c_act = (1 - (has_disable_flag * c_disable_flag[ind]));
        if (c_act == 1)
          dvol = 1.0/cell_f_vol[c_id];
        divut[c_id] *= dvol;
      });
    }
  }
  ctx.wait();
}

