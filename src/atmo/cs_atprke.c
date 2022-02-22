/*============================================================================
 * Modify the k-epsilon turbulence model for the atmospheric module.
 *============================================================================*/

/* This file is part of code_saturne, a general-purpose CFD tool.

   Copyright (C) 1998-2022 EDF S.A.

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

#include "cs_defs.h"

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

#include "bft_mem.h"

#include "cs_atmo.h"
#include "cs_air_props.h"
#include "cs_base.h"
#include "cs_field.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_thermal_model.h"
#include "cs_time_step.h"
#include "cs_turbulence_model.h"
#include "cs_velocity_pressure.h"

#include "cs_atmo_profile_std.h"
#include "cs_intprf.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_atprke.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_atprke.c

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

  /* Computesaturation q */
  cs_real_t qs = cs_air_yw_sat(t-tkelvi, pphy);

  /* Calcul de de */
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
 * \param[in, out]  tinstk      implicit part of the buoyancy term (for k)
 * \param[in, out]  smbrk       explicit part of the buoyancy term (for k)
 * \param[in, out]  smbre       explicit part of the buoyancy term (for eps)
 */
/*----------------------------------------------------------------------------*/

static void
_dry_atmosphere(const cs_real_t  cromo[],
                const cs_real_t  cpro_pcvto[],
                cs_real_t        tinstk[],
                cs_real_t        smbrk[],
                cs_real_t        smbre[])
{
  cs_mesh_quantities_t *fvq   = cs_glob_mesh_quantities;
  const cs_real_t *cell_f_vol = fvq->cell_vol;
  const cs_lnum_t n_cells     = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  const cs_real_t *grav = cs_glob_physical_constants->gravity;

  cs_field_t *f_thm = cs_thermal_model_field();

  /* Initialization
     ============== */

  cs_real_t prdtur = 1;

  if (f_thm != NULL) {
    prdtur = cs_field_get_key_double(f_thm,
                                     cs_field_key_id("turbulent_schmidt"));
  }

  /* Allocate work arrays */
  cs_real_3_t *grad;
  BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);

  cs_real_t *cvara_k = NULL, *cvara_ep = NULL;

  if (CS_F_(k) != NULL) {
    cvara_k  =  (cs_real_t *)CS_F_(k)->val_pre;
    cvara_ep =  (cs_real_t *)CS_F_(eps)->val_pre;
  }

  /* Compute potential temperature derivatives
     ======================================== */

  cs_real_t *cvara_tpp = f_thm->val_pre;

  /* Computation of the gradient of the potential temperature
   *
   * compute the turbulent production/destruction terms:
   * dry atmo: (1/turb_schmidt*theta)*(dtheta/dz)*gz */

  cs_field_gradient_scalar(f_thm,
                           true, // use_previous_t
                           1,    // inc
                           true, // iccocg
                           grad);

  /* Production and gravity terms
     TINSTK = P + G et TINSTE = P + (1 - CE3)*G */

  cs_real_t rho, visct, xeps, xk, ttke, gravke;
  cs_field_t *f_tke_buoy = cs_field_by_name_try("tke_buoyancy");
  cs_real_t *cpro_beta = NULL;
  cs_velocity_pressure_model_t *vp_model = cs_get_glob_velocity_pressure_model();
  if (vp_model->idilat == 0)
    cpro_beta = cs_field_by_name("thermal_expansion")->val;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    rho    = cromo[c_id];
    visct  = cpro_pcvto[c_id];
    xeps   = cvara_ep[c_id];
    xk     = cvara_k[c_id];
    ttke   = xk / xeps;
    cs_real_t beta = (cpro_beta == NULL) ? 1./cvara_tpp[c_id] : cpro_beta[c_id];

    gravke = beta *
      cs_math_3_dot_product(grad[c_id], grav) / prdtur;

    /* Implicit part (no implicit part for epsilon because the source
     * term is positive) */

    tinstk[c_id] += fmax(-rho*cell_f_vol[c_id]*cs_turb_cmu*ttke*gravke, 0.);

    /* Explicit part */
    smbre[c_id] += visct*fmax(gravke, 0.);
    smbrk[c_id] += visct*gravke;

    /* Save for post processing */
    if (f_tke_buoy != NULL)
      f_tke_buoy->val[c_id] = visct*gravke/rho;
  }

  BFT_FREE(grad);
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
 * \param[in]       cromo       density
 * \param[in]       cpro_pcvto  turbulent viscosity
 * \param[in, out]  tinstk      implicit part of the buoyancy term (for k)
 * \param[in, out]  smbrk       explicit part of the buoyancy term (for k)
 * \param[in, out]  smbre       explicit part of the buoyancy term (for eps)
 */
/*----------------------------------------------------------------------------*/

static void
_humid_atmosphere(const cs_real_t  cromo[],
                  const cs_real_t  cpro_pcvto[],
                  cs_real_t        tinstk[],
                  cs_real_t        smbrk[],
                  cs_real_t        smbre[])
{
  cs_mesh_quantities_t *fvq   = cs_glob_mesh_quantities;
  const cs_real_t *cell_f_vol = fvq->cell_vol;
  const cs_lnum_t n_cells     = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  const cs_real_3_t *cell_cen = (const cs_real_3_t *)fvq->cell_cen;

  const cs_real_t *grav = cs_glob_physical_constants->gravity;

  cs_field_t *f_thm = cs_thermal_model_field();

  /* Initialization
     ============== */

  cs_real_t prdtur = 1;

  if (f_thm != NULL) {
    prdtur = cs_field_get_key_double(f_thm,
                                     cs_field_key_id("turbulent_schmidt"));
  }

  /* Allocate work arrays */
  cs_real_3_t *grad;
  BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);

  cs_real_t *cvara_k = NULL, *cvara_ep = NULL;

  if (CS_F_(k) != NULL) {
    cvara_k  =  (cs_real_t *)CS_F_(k)->val_pre;
    cvara_ep =  (cs_real_t *)CS_F_(eps)->val_pre;
  }

  /* Compute potential temperature derivatives
     ======================================== */

  const cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();
  cs_real_t rvsra = phys_pro->rvsra;

  /* compute the production term in case of humid atmosphere
   * ie. when iatmos eq 2 */
  cs_real_t *etheta, *eq, *gravke_theta, *gravke_qw;
  BFT_MALLOC(etheta, n_cells_ext, cs_real_t);
  BFT_MALLOC(eq, n_cells_ext, cs_real_t);
  BFT_MALLOC(gravke_theta, n_cells_ext, cs_real_t);
  BFT_MALLOC(gravke_qw, n_cells_ext, cs_real_t);

  /* Computation of the gradient of the potentalp_bl temperature */

  cs_real_t *cvara_tpp  = f_thm->val_pre;
  cs_real_t *cvara_qw   = cs_field_by_name("ym_water")->val_pre;
  cs_real_t *cpro_pcliq = cs_field_by_name("liquid_water")->val;

  /* compute the coefficients etheta,eq */

  cs_real_t pphy, dum;

  cs_real_t *diag_neb = cs_field_by_name("nebulosity_diag")->val;
  cs_real_t *frac_neb = cs_field_by_name("nebulosity_frac")->val;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    /* calculate the physical pressure 'pphy' */

    if (cs_glob_atmo_option->meteo_profile == 0) {
      cs_atmo_profile_std(cell_cen[c_id][2], &pphy, &dum, &dum);
    }
    else if (cs_glob_atmo_option->meteo_profile == 1) {
      int nbmett = cs_glob_atmo_option->nbmett;
      int nbmetm = cs_glob_atmo_option->nbmetm;
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
           cvara_tpp[c_id],
           cvara_qw[c_id],
           cpro_pcliq[c_id],
           diag_neb[c_id],
           frac_neb[c_id],
           &(etheta[c_id]),
           &(eq[c_id]));
  }

  /* compute the turbulent production/destruction terms:
   * humid atmo: (1/turb_schmidt*theta_v)*(dtheta_l/dz)*gz */

  cs_field_gradient_scalar(f_thm,
                           true, // use_previous_t
                           1,    // inc
                           true, // iccocg
                           grad);

  /* Production and gravity terms
   * TINSTK = P + G et TINSTE = P + (1 - CE3)*G */

  /* Now store the production term due to theta_liq in gravke_theta */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_real_t qw         = cvara_qw[c_id];    /* total water content */
    cs_real_t qldia      = cpro_pcliq[c_id];  /* liquid water content */
    cs_real_t theta_virt = cvara_tpp[c_id]*(1. + (rvsra - 1.)*qw - rvsra*qldia);

    gravke_theta[c_id] =   etheta[c_id]
                         * (cs_math_3_dot_product(grad[c_id], grav)
                         / (theta_virt*prdtur));
  }

  /* Gradient of humidity and it's associated production term
     -------------------------------------------------------- */

  /* Compute the turbulent production/destruction terms:
   * humid atmo: (1/turb_schmidt*theta_v)*(dtheta_l/dz)*gz */

  cs_field_gradient_scalar(cs_field_by_name("ym_water"),
                           true, // use_previous_t
                           1,    // inc
                           true, // iccocg
                           grad);

  /* Production and gravity terms
   * TINSTK = P + G et TINSTE = P + (1 - CE3)*G */

  /* Store the production term due to qw in gravke_qw */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_real_t qw         = cvara_qw[c_id];    /* total water content*/
    cs_real_t qldia      = cpro_pcliq[c_id];  /* liquid water content */
    cs_real_t theta_virt = cvara_tpp[c_id]*(1. + (rvsra - 1.)*qw - rvsra*qldia);

    gravke_qw[c_id] =   eq[c_id]
                      * (cs_math_3_dot_product(grad[c_id], grav)
                      / (theta_virt*prdtur));
  }

  /* Finalization */

  cs_real_t rho, visct, xeps, xk, ttke, gravke;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    rho   = cromo[c_id];
    visct = cpro_pcvto[c_id];
    xeps  = cvara_ep[c_id];
    xk    = cvara_k[c_id];
    ttke  = xk / xeps;

    gravke  = gravke_theta[c_id] + gravke_qw[c_id];

    /* Implicit part (no implicit part for epsilon because the source
     * term is positive) */

    tinstk[c_id] += fmax(-rho*cell_f_vol[c_id]*cs_turb_cmu*ttke*gravke, 0.);

    /* Explicit part */
    smbre[c_id] += visct*fmax(gravke, 0.);
    smbrk[c_id] += visct*gravke;
  }

  BFT_FREE(etheta);
  BFT_FREE(eq);
  BFT_FREE(gravke_theta);
  BFT_FREE(gravke_qw);

  BFT_FREE(grad);
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
 * Adjunction of a production term for buyancy in the \f$k-\varepsilon\f$
 * model in the context of the atmospheric module
 * g = g*grad(theta)/prdtur/theta
 *
 * \param[in, out]  tinstk    Implicit part of the buoyancy term (for k)
 * \param[in, out]  smbrk     Explicit part of the buoyancy term (for k)
 * \param[in, out]  smbre     Explicit part of the buoyancy term (for eps)
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_buoyancy_ke_prod(cs_real_t  tinstk[],
                         cs_real_t  smbrk[],
                         cs_real_t  smbre[])
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
                    tinstk,
                    smbrk,
                    smbre);

  else if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_HUMID)
    _humid_atmosphere(cromo,
                      cpro_pcvto,
                      tinstk,
                      smbrk,
                      smbre);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
