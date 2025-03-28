/*============================================================================
 * Modify the k-epsilon turbulence model for the atmospheric module.
 *============================================================================*/

/* This file is part of code_saturne, a general-purpose CFD tool.

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
   Street, Fifth Floor, Boston, MA 02110-1301, USA. */

/*----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_mem.h"

#include "atmo/cs_atmo.h"
#include "atmo/cs_air_props.h"
#include "base/cs_base.h"
#include "base/cs_field.h"
#include "base/cs_field_operator.h"
#include "base/cs_field_pointer.h"
#include "base/cs_math.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_physical_constants.h"
#include "pprt/cs_physical_model.h"
#include "base/cs_thermal_model.h"
#include "base/cs_time_step.h"
#include "turb/cs_turbulence_model.h"
#include "base/cs_velocity_pressure.h"

#include "atmo/cs_atmo_profile_std.h"
#include "atmo/cs_intprf.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "atmo/cs_atprke.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_atprke.cpp

  Modify the \f$k-\varepsilon\f$ turbulence model formulation
  for the atmospheric module.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Compute etheta and eq variable knowing the saturation.
 *
 * \param[in]       pphy        pressure [Pa]
 * \param[in]       thetal      Liquid potential temperature
 * \param[in]       qw          total water amount
 * \param[in]       qldia       mean liquid water content
 * \param[in]       xnebdia     nebulosity after the diagnostic
 * \param[in]       xnn         second order moment "n" <s'ql'>/(2*sigma_s**2)
 * \param[out]      etheta      sensible heat part of buoyancy flux
 * \param[out]      eq          latent heat part of buoyancy flux
 */
/*----------------------------------------------------------------------------*/

static void
_etheq(cs_real_t   pphy,
       cs_real_t   thetal,
       cs_real_t   qw,
       cs_real_t   qldia,
       cs_real_t   xnebdia,
       cs_real_t   xnn,
       cs_real_t  *etheta,
       cs_real_t  *eq)
{
  const cs_fluid_properties_t *phys_pro = cs_glob_fluid_properties;

  cs_real_t rair = phys_pro->r_pg_cnst;
  cs_real_t p0 = phys_pro->p0;
  cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

  /* Most probable and simplest case
     =============================== */

  cs_real_t rvsra = phys_pro->rvsra;
  if (qldia <= 0. || cs_glob_atmo_option->subgrid_model == 0) {
    *etheta = 1.;
    *eq     = (rvsra - 1.)*thetal;
    return;
  }

  /* Initialization
     ============== */

  cs_real_t cp = phys_pro->cp0;
  cs_real_t clatev = phys_pro->clatev;

  cs_real_t lscp = clatev/cp;

  *etheta = 1.;
  *eq     = (rvsra - 1.)*thetal;

  /* Sub-grid model
     ============== */

  /* standard k-eps */

  /* Compute "liquid" temperature */
  cs_real_t tl = thetal*pow(p0/pphy, -rair/cp);

  /* Constants */
  cs_real_t beta1 = cs_math_pow2(clatev)/(cp*(rair*rvsra)*cs_math_pow2(tl));

  /* Compute saturation ql */
  cs_real_t qsl = cs_air_yw_sat(tl-tkelvi, pphy);

  cs_real_t a = 1./(1. + beta1*qsl);

  cs_real_t alpha1 = qsl*beta1*pow(pphy/p0, rair/cp)/lscp;

  /* Potential temperature (dry) */
  cs_real_t theta = thetal + (clatev/cp)*pow(p0/pphy, rair/cp)*qldia;

  /* Compute "thermodynamic" temperature */
  cs_real_t t = theta*pow(p0/pphy, -rair/cp);

  /* Compute q */
  cs_real_t q = qw - qldia;

  /* Compute saturation q */
  cs_real_t qs = cs_air_yw_sat(t-tkelvi, pphy);

  /* Compute de */
  cs_real_t de = (clatev/cp)*pow(p0/pphy, rair/cp) - rvsra*theta;

  /* Constants for moddis = 3 */
  cs_real_t beta2  = cs_math_pow2(clatev)/(cp*(rair*rvsra)*cs_math_pow2(t));
  cs_real_t aa     = 1./(1. + beta2*qs);
  cs_real_t alpha2 = qs*(beta2*cp/clatev)*pow(pphy/p0, rair/cp);

  /* New computation of d for moddis = 3 */
  cs_real_t dd =   (clatev/cp)*pow(p0/pphy, rair/cp)
                 * (1. + (rvsra - 1. )*q - qldia) - rvsra*theta;

  /* Nebulosity */
  cs_real_t rn = xnebdia;

  if (cs_glob_atmo_option->subgrid_model == 1) {
    /* Bechtold et al. 1995 */
    *etheta = 1. - a*alpha1*de*xnn;
    *eq     = (rvsra - 1.)*theta + a*de*xnn;
  }
  else if (cs_glob_atmo_option->subgrid_model == 2) {
    /* Bouzereau et al. 2004 */
    *etheta = 1. + (rvsra - 1.)*q - qldia - a*alpha1*dd*xnn;
    *eq     = (rvsra - 1.)*theta + a*dd*xnn;
  }
  else if (cs_glob_atmo_option->subgrid_model == 3) {
    /* Cuijpers et Duynkerke 1993, etc.
     * Linear interpolation between saturated and non-saturated cases
     * (coefficient of partial nebulosity r) */
    *etheta = 1. + (rvsra - 1.)*q - rn*(qldia + aa*alpha2*dd);
    *eq     = (rvsra - 1.)*theta + rn*aa*dd;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Modify the \f$k-\varepsilon\f$ turbulence model
 *        formulation for the atmospheric module with a dry atmosphere.
 *
 * Adjunction of a production term for buyancy in the \f$k-\varepsilon\f$
 * model in the context of the atmospheric module
 * g = g*grad(theta)/prdtur/theta
 *
 * \param[in]       cromo       density
 * \param[in]       cpro_pcvto  turbulent viscosity
 * \param[in, out]  gk          explicit part of the buoyancy term (for k)
 */
/*----------------------------------------------------------------------------*/

static void
_dry_atmosphere(const cs_real_t  cromo[],
                const cs_real_t  cpro_pcvto[],
                cs_real_t        gk[])
{
  const cs_lnum_t n_cells     = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  const cs_real_t *grav = cs_glob_physical_constants->gravity;

  cs_field_t *f_thm = cs_thermal_model_field();

  /* Initialization
     ============== */

  cs_real_t prdtur = 1;

  if (f_thm != nullptr) {
    prdtur = cs_field_get_key_double(f_thm,
                                     cs_field_key_id("turbulent_schmidt"));
  }

  /* Allocate work arrays */
  cs_real_3_t *grad;
  CS_MALLOC(grad, n_cells_ext, cs_real_3_t);

  /* Compute potential temperature derivatives
     ======================================== */

  /* Computation of the gradient of the potential temperature
   *
   * compute the turbulent production/destruction terms:
   * dry atmo: (1/turb_schmidt*theta)*(dtheta/dz)*gz */

  cs_field_gradient_scalar(f_thm,
                           false, // use_previous_t
                           1,    // inc
                           grad);

  /* TKE Production by gravity term G */

  cs_real_t rho, visct, gravke;
  cs_field_t *f_tke_buoy = cs_field_by_name_try("tke_buoyancy");
  cs_field_t *f_beta = cs_field_by_name("thermal_expansion");

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    rho    = cromo[c_id];
    visct  = cpro_pcvto[c_id];
    cs_real_t beta = f_beta->val[c_id];

    gravke = beta * cs_math_3_dot_product(grad[c_id], grav) / prdtur;

    /* Explicit Buoyant terms */
    gk[c_id] = visct*gravke;

    /* Save for post processing */
    if (f_tke_buoy != nullptr)
      f_tke_buoy->val[c_id] = visct*gravke/rho;
  }

  CS_FREE(grad);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Modify the \f$k-\varepsilon\f$ turbulence model
 *        formulation for the atmospheric module with a humid atmosphere.
 *
 * Adjunction of a production term for buyancy in the \f$k-\varepsilon\f$
 * model in the context of the atmospheric module
 * g = g*grad(theta)/prdtur/theta
 *
 * \param[in]       cpro_pcvto  turbulent viscosity
 * \param[in, out]  gk          explicit part of the buoyancy term (for k)
 */
/*----------------------------------------------------------------------------*/

static void
_humid_atmosphere(const cs_real_t  cpro_pcvto[],
                  cs_real_t        gk[])
{
  cs_mesh_quantities_t *fvq   = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells     = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  const cs_real_3_t *cell_cen = (const cs_real_3_t *)fvq->cell_cen;

  const cs_real_t *grav = cs_glob_physical_constants->gravity;
  const cs_fluid_properties_t *fluid_props = cs_glob_fluid_properties;

  cs_field_t *f_thm = cs_thermal_model_field();

  /* Initialization
     ============== */

  cs_real_t prdtur = 1;

  if (f_thm != nullptr) {
    prdtur = cs_field_get_key_double(f_thm,
                                     cs_field_key_id("turbulent_schmidt"));
  }

  cs_field_t *f_beta = cs_field_by_name("thermal_expansion");

  /* Allocate work arrays */
  cs_real_3_t *grad;
  CS_MALLOC(grad, n_cells_ext, cs_real_3_t);

  /* Compute potential temperature derivatives
     ======================================== */

  /* compute the production term in case of humid atmosphere
   * ie. when iatmos eq 2 */
  cs_real_t *etheta, *eq, *gravke_theta, *gravke_qw;
  CS_MALLOC(etheta, n_cells_ext, cs_real_t);
  CS_MALLOC(eq, n_cells_ext, cs_real_t);
  CS_MALLOC(gravke_theta, n_cells_ext, cs_real_t);
  CS_MALLOC(gravke_qw, n_cells_ext, cs_real_t);

  /* Computation of the gradient of the potentalp_bl temperature */

  cs_real_t *cvar_tpp  = f_thm->val;
  cs_real_t *cvar_qw   = cs_field_by_name("ym_water")->val;
  cs_real_t *cpro_pcliq = cs_field_by_name("liquid_water")->val;

  /* compute the coefficients etheta,eq */

  cs_real_t pphy, dum;

  cs_real_t *diag_neb = cs_field_by_name("nebulosity_diag")->val;
  cs_real_t *frac_neb = cs_field_by_name("nebulosity_frac")->val;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    /* calculate the physical pressure 'pphy' */

    if (cs_glob_atmo_option->meteo_profile == 0) {
      cs_atmo_profile_std(0., /* z_ref */
                          fluid_props->p0,
                          fluid_props->t0,
                          cell_cen[c_id][2], &pphy, &dum, &dum);
    }
    else if (cs_glob_atmo_option->meteo_profile == 1) {
      int nbmett = cs_glob_atmo_option->met_1d_nlevels_t;
      int nbmetm = cs_glob_atmo_option->met_1d_ntimes;
      pphy = cs_intprf(nbmett,
                       nbmetm,
                       cs_glob_atmo_option->z_temp_met,
                       cs_glob_atmo_option->time_met,
                       cs_glob_atmo_option->hyd_p_met,
                       cell_cen[c_id][2],
                       cs_glob_time_step->t_cur);
    }
    else {
      pphy = cs_field_by_name("meteo_pressure")->val[c_id];
    }

    _etheq(pphy,
           cvar_tpp[c_id],
           cvar_qw[c_id],
           cpro_pcliq[c_id],
           diag_neb[c_id],
           frac_neb[c_id],
           &(etheta[c_id]),
           &(eq[c_id]));
  }

  /* compute the turbulent production/destruction terms:
   * humid atmo: (beta/turb_schmidt*theta_v)*(dtheta_l/dz)*gz
   * with beta = 1/theta_v */

  cs_field_gradient_scalar(f_thm,
                           false, // use_previous_t
                           1,    // inc
                           grad);

  /* Production and gravity terms
   * RHS k  = P + G et RHS ep = ce1 P + ce3 * G */

  /* Now store the production term due to theta_liq in gravke_theta */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    gravke_theta[c_id] =   etheta[c_id] * f_beta->val[c_id]
                         * cs_math_3_dot_product(grad[c_id], grav)
                         / prdtur;
  }

  /* Gradient of humidity and it's associated production term
     -------------------------------------------------------- */

  /* Compute the turbulent production/destruction terms:
   * humid atmo: (1/turb_schmidt*theta_v)*(dtheta_l/dz)*gz */

  cs_field_gradient_scalar(cs_field_by_name("ym_water"),
                           false, // use_previous_t
                           1,    // inc
                           grad);

  /* Production and gravity terms
   * RHS k  = P + G et RHS ep = ce1 P + ce3 * G */

  /* Store the production term due to qw in gravke_qw */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    gravke_qw[c_id] =   eq[c_id] * f_beta->val[c_id]
                      * cs_math_3_dot_product(grad[c_id], grav)
                      / prdtur;
  }

  /* Finalization */

  cs_real_t visct, gravke;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    visct = cpro_pcvto[c_id];

    gravke  = gravke_theta[c_id] + gravke_qw[c_id];

    /* Explicit Buoyant terms */
    gk[c_id] = visct*gravke;

  }

  CS_FREE(etheta);
  CS_FREE(eq);
  CS_FREE(gravke_theta);
  CS_FREE(gravke_qw);

  CS_FREE(grad);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Modify the \f$k-\varepsilon\f$ turbulence model
 *        formulation for the atmospheric module
 *
 * Adjunction of a production term for buoyancy in the \f$k-\varepsilon\f$
 * model in the context of the atmospheric module
 * g = g*grad(theta)/prdtur/theta
 *
 * \param[in, out]  gk        Explicit part of the buoyancy term (for k)
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_buoyancy_ke_prod(cs_real_t  gk[])
{
  /* Time extrapolation? */
  int key_t_ext_id = cs_field_key_id("time_extrapolated");

  /* Pointer to density and turbulent viscosity */
  cs_real_t *cromo      =  (cs_real_t *)CS_F_(rho)->val;
  cs_real_t *cpro_pcvto =  (cs_real_t *)CS_F_(mu_t)->val;

  const cs_time_scheme_t *time_scheme = cs_glob_time_scheme;

  if (time_scheme->isto2t > 0) {
    if (cs_field_get_key_int(CS_F_(rho), key_t_ext_id) > 0) {
      cromo = CS_F_(rho)->val_pre;
    }
    if (cs_field_get_key_int(CS_F_(mu), key_t_ext_id) > 0) {
      cpro_pcvto = CS_F_(mu_t)->val_pre;
    }
  }

  /* Compute potential temperature derivatives
     ======================================== */

  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_DRY)
    _dry_atmosphere(cromo,
                    cpro_pcvto,
                    gk);

  else if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_HUMID)
    _humid_atmosphere(cpro_pcvto,
                      gk);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
