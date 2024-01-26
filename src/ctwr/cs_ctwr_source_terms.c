/*============================================================================
 * Cooling towers related functions
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_locator.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_nodal_extract.h"

#include "cs_air_props.h"
#include "cs_atmo.h"
#include "cs_base.h"
#include "cs_boundary_conditions.h"
#include "cs_boundary_zone.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_parameters.h"
#include "cs_parall.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_post.h"
#include "cs_prototypes.h"
#include "cs_restart.h"
#include "cs_selector.h"
#include "cs_thermal_model.h"
#include "cs_velocity_pressure.h"
#include "cs_volume_zone.h"

#include "cs_ctwr.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ctwr_source_terms.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------
 * Compute the Lewis factor used for the evaluation of the heat transfer
 * phase change source terms
 *
 * parameters:
 *   evap_model  <-- Evaporation model: CS_CTWR_POPPE or CS_CTWR_MERKEL
 *   molmassrat  <-- Dry air to water vapor molecular mass ratio
 *   x           <-- Humidity
 *   x_s_tl      <-- Saturation humidity at the temperature of the liquid
 *
 * returns:
 *   le_f        --> Lewis factor
 *----------------------------------------------------------------------------*/

static cs_real_t
_lewis_factor(const int        evap_model,
              const cs_real_t  molmassrat,
              const cs_real_t  x,
              const cs_real_t  x_s_tl)
{
  /* Merkel Model
     Hypothesis of unity Lewis factor */
  cs_real_t le_f = 1.;

  if (evap_model == CS_CTWR_POPPE) {
    /* Poppe evaporation model
       Compute Lewis factor using Bosnjakovic hypothesis
       NB: clippings ensuring xi > 1 and le_f > 0 */
    cs_real_t xi = (molmassrat + x_s_tl)/(molmassrat + CS_MIN(x, x_s_tl));
    if ((xi - 1.) < 1.e-15)
      le_f = pow(0.866,(2./3.));
    else
      le_f = pow(0.866,(2./3.))*(xi-1.)/log(xi);
  }

  return le_f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Phase change source terms - Exchange terms between the injected
 *        liquid and the water vapor phase in the bulk, humid air
 *
 * \param[in]     f_id          field id
 * \param[in,out] exp_st        Explicit source term
 * \param[in,out] imp_st        Implicit source term
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_source_term(int              f_id,
                    cs_real_t        exp_st[],
                    cs_real_t        imp_st[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_2_t *i_face_cells
    = (const cs_lnum_2_t *)(m->i_face_cells);
  const cs_lnum_t n_i_faces = m->n_i_faces;

  const cs_real_t *cell_f_vol = cs_glob_mesh_quantities->cell_f_vol;

  cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();
  cs_air_fluid_props_t *air_prop = cs_glob_air_props;

  /* Water / air molar mass ratio */
  const cs_real_t molmassrat = air_prop->molmass_rat;

  cs_real_t  *rho_h = (cs_real_t *)CS_F_(rho)->val; /* Humid air
                                                       (bulk) density */
  cs_real_3_t *vel_h = (cs_real_3_t *)CS_F_(vel)->val; /* Humid air (bulk)
                                                          velocity*/

  cs_real_t *ym_w = (cs_real_t *)CS_F_(ym_w)->val; /* Water mass fraction
                                                     in humid air */

  cs_real_t *t_h = cs_field_by_name("temperature")->val; /* Humid air
                                                            temperature */
  cs_real_t *x = cs_field_by_name("humidity")->val; /* Humidity in air (bulk) */
  cs_real_t *x_s = cs_field_by_name("x_s")->val;
  cs_real_t *vel_l = cs_field_by_name("vertvel_l")->val;  /* Liquid velocity
                                                             in packing */

  cs_real_t *t_l_p = cs_field_by_name("temp_l_packing")->val;
  cs_real_t *h_l_p = cs_field_by_name("h_l_packing")->val;
  cs_real_t *y_l_p = CS_F_(y_l_pack)->val_pre;

  /* Variable and properties for rain drops */
  cs_field_t *cfld_yp = cs_field_by_name("y_l_r");     /* Rain mass fraction */
  cs_field_t *cfld_yh_rain = cs_field_by_name("yh_l_r"); /* Yp times Tp */
  cs_field_t *cfld_drift_vel = cs_field_by_name("drift_vel_y_l_r"); /* Rain drift
                                                                     velocity */

  /* Rain inner mass flux */
  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  cs_real_t *imasfl_r = cs_field_by_id(cs_field_get_key_int(cfld_yp, kimasf))->val;

  cs_real_t vertical[3], horizontal[3];

  const cs_ctwr_option_t *ct_opt = cs_glob_ctwr_option;

  int evap_model = ct_opt->evap_model;


  /* Cooling tower zones */
  cs_ctwr_zone_t **_ct_zone = cs_get_glob_ctwr_zone();
  const int *_n_ct_zones = cs_get_glob_ctwr_n_zones();

  /* Need to cook up the cell value of the liquid mass flux
     In the old code, it seems to be taken as the value of the
     face mass flux upstream of the cell */

  cs_real_t v_air = 0.;

  cs_real_t mass_flux_h = 0.; /* Highly suspicious for rain
                                 zones - not recomputed */

  /* Identify the source term formulation for the required field */

  const cs_field_t *f = cs_field_by_id(f_id);

  cs_real_t *f_var = f->val;  /* Field variable */

  /* Compute the source terms */

  /* Fields for source terms post-processing */
  cs_real_t *evap_rate_pack = NULL;
  cs_real_t *evap_rate_rain = NULL;
  evap_rate_pack = cs_field_by_name("evaporation_rate_packing")->val;
  evap_rate_rain = cs_field_by_name("evaporation_rate_rain")->val;

  cs_real_t *thermal_power_pack = NULL;
  cs_real_t *thermal_power_rain = NULL;
  thermal_power_pack = cs_field_by_name("thermal_power_packing")->val;
  thermal_power_rain = cs_field_by_name("thermal_power_rain")->val;

  /* Table to track cells belonging to packing zones */
  cs_lnum_t  *packing_cell;
  BFT_MALLOC(packing_cell, m->n_cells_with_ghosts, int);
#   pragma omp parallel for if (m->n_cells> CS_THR_MIN)
  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells_with_ghosts; cell_id++)
    packing_cell[cell_id] = -1;

  /* Air / fluid properties */
  cs_real_t cp_v = air_prop->cp_v;
  cs_real_t cp_l = air_prop->cp_l;
  cs_real_t hv0 = air_prop->hv0;
  cs_real_t rho_l = air_prop->rho_l;
  cs_real_t visc = fp->viscl0;
  cs_real_t  p0 = fp->p0;
  cs_real_t lambda_h = air_prop->lambda_h;
  cs_real_t droplet_diam  = air_prop->droplet_diam;
  cs_real_t sigma = air_prop->sigma;

  if (evap_model != CS_CTWR_NONE) {

    /* =========================
       Phase change source terms
       ========================= */

    cs_math_3_normalize(cs_glob_physical_constants->gravity, vertical);

    vertical[0] *= -1;
    vertical[1] *= -1;
    vertical[2] *= -1;

    horizontal[0] = vertical[0] -1.;
    horizontal[1] = vertical[1] -1.;
    horizontal[2] = vertical[2] -1.;

    /* Phase change in the packing zones
       Between the liquid film and the humid air
       ========================================= */

    for (int ict = 0; ict < *_n_ct_zones; ict++) {

      cs_ctwr_zone_t *ct = _ct_zone[ict];

      /* We skip this if we are in injection zone */
      if (ct->type == CS_CTWR_INJECTION)
        continue;

      /* Packing zone characteristics */
      cs_real_t a_0 = ct->xap;
      cs_real_t xnp = ct->xnp;
      int zone_type = ct->type;

      const cs_lnum_t *ze_cell_ids = cs_volume_zone_by_name(ct->name)->elt_ids;

      for (cs_lnum_t j = 0; j < ct->n_cells; j++) {

        cs_lnum_t cell_id = ze_cell_ids[j];
        /* Identify packing cells ids */
        if (ct->type != CS_CTWR_INJECTION)
          packing_cell[cell_id] = ict;

        /* For correlations, T_h cannot be greater than T_l */
        cs_real_t temp_h = CS_MIN(t_h[cell_id], t_l_p[cell_id]);

        /* Saturation humidity at humid air temperature */
        cs_real_t x_s_th = cs_air_x_sat(temp_h, p0);

        /* Saturation humidity at injected liquid temperature */
        cs_real_t x_s_tl = cs_air_x_sat(t_l_p[cell_id], p0);

        /*--------------------------------------------*
         * Counter or cross flow packing zone         *
         *--------------------------------------------*/

        if (zone_type == CS_CTWR_COUNTER_CURRENT) {
          /* Counter flow packing */
          v_air = CS_ABS(cs_math_3_dot_product(vel_h[cell_id], vertical));
        }
        else if (zone_type == CS_CTWR_CROSS_CURRENT) {
          /* Cross flow packing */
          v_air = CS_ABS(cs_math_3_dot_product(vel_h[cell_id], horizontal));
        }

        /* Dry air flux */
        mass_flux_h = rho_h[cell_id] * v_air * (1. - ym_w[cell_id]);

        /* Liquid mass flux */
        cs_real_t mass_flux_l = rho_h[cell_id] * y_l_p[cell_id] * vel_l[cell_id];
        cs_real_t mass_flux_l_oy = rho_h[cell_id] * vel_l[cell_id];

        /* Evaporation coefficient 'Beta_x' (kg/m2/s) times exchange surface
         * per unit of volume 'ai' (m2/m3)*/
        cs_real_t beta_x_ai = 0.;
        cs_real_t beta_x_ai_oy = 0.;

        /* There is evaporation only if we have an injected liquid flow */
        if (mass_flux_l > 0.){
          beta_x_ai = a_0 * mass_flux_l * pow((mass_flux_h/mass_flux_l), xnp);
          beta_x_ai_oy = a_0 * mass_flux_l_oy * pow((mass_flux_h/mass_flux_l), xnp);
        }

        /* Source terms for the different equations */

        /* Humid air mass source term */
        cs_real_t mass_source = 0.0;
        cs_real_t mass_source_oy = 0.0;

        if (x[cell_id] <= x_s_th) {
          mass_source = beta_x_ai * (x_s_tl - x[cell_id]);
          mass_source_oy = beta_x_ai_oy * (x_s_tl - x[cell_id]);
        }
        else {
          mass_source = beta_x_ai * (x_s_tl - x_s_th);
          mass_source_oy = beta_x_ai_oy * (x_s_tl - x_s_th);
        }
        mass_source = CS_MAX(mass_source, 0.);
        mass_source_oy = CS_MAX(mass_source_oy, 0.);

        cs_real_t vol_mass_source = mass_source * cell_f_vol[cell_id];
        cs_real_t vol_mass_source_oy = mass_source_oy * cell_f_vol[cell_id];
        cs_real_t vol_beta_x_ai = beta_x_ai * cell_f_vol[cell_id];

        /* Humid air thermal source term */
        cs_real_t cp_h = cs_air_cp_humidair(x[cell_id], x_s[cell_id]);

        /* Global mass source term for continuity (pressure) equation
         * Note that rain is already considered in the bulk, so inner
         * mass transfer between liquid and vapor disappears */
        if (f_id == (CS_F_(p)->id)) {
          /* Warning: not multiplied by Cell volume! no addition neither */
          exp_st[cell_id] = mass_source;

          /* Saving evaporation rate for post-processing */
          evap_rate_pack[cell_id] = mass_source;
        }

        /* Water mass fraction equation except rain */
        else if (f_id == (CS_F_(ym_w)->id)) {
          //TODO add mass_from_rain
          exp_st[cell_id] += vol_mass_source * (1. - f_var[cell_id]);
          imp_st[cell_id] += vol_mass_source;
        }

        /* Injected liquid mass equation (solve in drift model form) */
        else if (f_id == (CS_F_(y_l_pack)->id)) {
          exp_st[cell_id] -= vol_mass_source_oy * y_l_p[cell_id];
          imp_st[cell_id] += vol_mass_source_oy;
        }

        /* Humid air temperature equation */
        else if (f_id == (CS_F_(t)->id)) {
          // FIXME source term for theta_l instead...
          /* Because the writing is in a non-conservative form */
          cs_real_t l_imp_st = vol_mass_source * cp_h;
          cs_real_t l_exp_st = 0.;
          cs_real_t le_f = _lewis_factor(evap_model, molmassrat,
                                         x[cell_id], x_s_tl);
          if (x[cell_id] <= x_s_th) {
            /* Implicit term */
            l_imp_st += vol_beta_x_ai * (le_f * cp_h
                                          + (x_s_tl - x[cell_id]) * cp_v
                                            / (1. + x[cell_id]));
            l_exp_st = l_imp_st * (t_l_p[cell_id] - f_var[cell_id]);
          }
          else {
            cs_real_t coeft = le_f * cp_h;
            /* Implicit term */
            l_imp_st += vol_beta_x_ai * (coeft + (x_s_tl - x_s_th) * cp_l
                                                 / (1. + x[cell_id]));
            /* Explicit term */
            //FIXME : why does a hv0 term appears while we work on temperature ?
            l_exp_st = vol_beta_x_ai * (coeft * t_l_p[cell_id]
                                        + (x_s_tl - x_s_th)
                                          * (cp_v * t_l_p[cell_id] + hv0)
                                          / (1. + x[cell_id]))
                       - l_imp_st * f_var[cell_id];
          }
          imp_st[cell_id] += CS_MAX(l_imp_st, 0.);
          exp_st[cell_id] += l_exp_st;
        }

        /* Injected liquid enthalpy equation (solve in drift model form)
         * NB: it is in fact "y_l_p x h_l" */
        else if (f_id == (CS_F_(yh_l_pack)->id)) {
          /* Liquid temperature in Kelvin */
          cs_real_t t_l_k = t_l_p[cell_id]
                            + cs_physical_constants_celsius_to_kelvin;
          cs_real_t l_exp_st = 0.;

          /* Note: the solved variable is yl_p.hl_p so the source term associated
           * to evaporation is: Gamma/y_lp * (yl_p.h_lp) */
          /* Implicit term */
          cs_real_t l_imp_st = vol_mass_source_oy;

          l_exp_st -= l_imp_st * f_var[cell_id];

          cs_real_t le_f = _lewis_factor(evap_model,
                                         molmassrat,
                                         x[cell_id],
                                         x_s_tl);
          /* Under saturated */
          if (x[cell_id] <= x_s_th) {
            /* Explicit term */
            l_exp_st -= vol_beta_x_ai * ((x_s_tl - x[cell_id])
                                         * (cp_v * t_l_k + hv0)
                                         + le_f * cp_h
                                           * (t_l_p[cell_id] - t_h[cell_id]));
          }
          /* Over saturated */
          else {
            cs_real_t coefh = le_f * cp_h;
            /* Explicit term */
            l_exp_st += vol_beta_x_ai * (coefh * (t_h[cell_id] - t_l_p[cell_id])
                                         + (x_s_tl - x_s_th) / (1. + x[cell_id])
                                           * (cp_l * t_h[cell_id]
                                              - (cp_v * t_l_p[cell_id] + hv0)));
          }
          /* Because we deal with an increment */
          exp_st[cell_id] += l_exp_st;
          imp_st[cell_id] += CS_MAX(l_imp_st, 0.);

          /* Saving thermal power for post-processing */
          thermal_power_pack[cell_id] = -(l_exp_st + l_imp_st * f_var[cell_id])
                                         / cell_f_vol[cell_id];
        }
      } /* end loop over the cells of a packing zone */

    } /* end loop over all the packing zones */

    /* Phase change in the whole domain
       Between the rain drops and the humid air
       ======================================== */

    if (cfld_yp != NULL) {
      cs_real_t *y_rain = (cs_real_t *)cfld_yp->val;
      cs_real_t *t_l_r = (cs_real_t *)cs_field_by_name("temp_l_r")->val;

      for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {

        if (y_rain[cell_id] > 0.) {

          /* For correlations, T_h cannot be greater than T_p */
          cs_real_t temp_h = CS_MIN(t_h[cell_id], t_l_r[cell_id]);

          /* Saturation humidity at the temperature of the humid air */
          cs_real_t x_s_th = cs_air_x_sat(temp_h, p0);

          /* Saturation humidity at the temperature of the rain drop  */
          cs_real_t x_s_tl = cs_air_x_sat(t_l_r[cell_id], p0);

          cs_real_3_t *drift_vel_rain
            = (cs_real_3_t *restrict)(cfld_drift_vel->val);
          cs_real_t drift_vel_mag = cs_math_3_norm(drift_vel_rain[cell_id]);

          /* Lewis factor computation */
          cs_real_t le_f = _lewis_factor(evap_model,
                                         molmassrat,
                                         x[cell_id],
                                         x_s_tl);

          cs_real_t cp_h = cs_air_cp_humidair(x[cell_id], x_s[cell_id]);

          /* Rain droplets Reynolds number */
          cs_real_t rey = rho_h[cell_id] * drift_vel_mag * droplet_diam / visc;

          /* Prandtl number */
          cs_real_t pr = cp_h * visc / lambda_h;

          /* Nusselt number correlations */
          /* Ranz-Marshall or Hughmark when rey <= 776.06 && pr <= 250. */
          cs_real_t nusselt = 2. + 0.6 * sqrt(rey) * pow(pr,(1./3.));
          /* Hughmark when rey > 776.06 && pr <= 250. */
          if (rey > 776.06 && pr <= 250.) {
            nusselt = 2. + 0.27 * pow(rey, 0.62) * pow(pr,(1./3.));
          }

          /* Convective exchange coefficient 'a_c' */
          cs_real_t a_c = (nusselt * lambda_h) / droplet_diam;

          /* beta_x coefficient */
          cs_real_t beta_x = a_c / (le_f * cp_h);

          /* Exchange surface area per unit volume based on the total droplets
           * surface in the cell
           * NOTE: Use rho_h to compute the number of particles per unit volume
           * since conservation equation for Y_p based on rho_h
           *   --> Should really be rho_mixture!?
           * Use the symmetric relationship:
           *   a_i = 6*alpha_p*(1.-alpha_p)/droplet_diam
           * where alpha_p is the droplets volume fraction
           * - this kills transfer when there is only one phase (pure humid air
           *   or pure rain) */
          cs_real_t vol_frac_rain = y_rain[cell_id] * rho_h[cell_id] / rho_l;
          cs_real_t vol_frac_rain_oy = rho_h[cell_id] / rho_l;

          if (vol_frac_rain >= 1.0)
            vol_frac_rain = 1.0;
          cs_real_t a_i =  6.0 * vol_frac_rain * (1.0 - vol_frac_rain)
                           / droplet_diam;
          cs_real_t a_i_oy = 6.0 * vol_frac_rain_oy * (1.0 - vol_frac_rain)
                           / droplet_diam;

          /* Evaporation coefficient 'Beta_x' times exchange surface 'ai' */
          cs_real_t beta_x_ai = beta_x * a_i;
          cs_real_t beta_x_ai_oy = beta_x * a_i_oy;

          /* Source terms for the different equations */

          /* Humid air mass source term */
          cs_real_t mass_source = 0.0;
          cs_real_t mass_source_oy = 0.0;

          if (x[cell_id] <= x_s_th) {
            mass_source = beta_x_ai * (x_s_tl - x[cell_id]);
            mass_source_oy = beta_x_ai_oy * (x_s_tl - x[cell_id]);
          }
          else {
            mass_source = beta_x_ai * (x_s_tl - x_s_th);
            mass_source_oy = beta_x_ai_oy * (x_s_tl - x_s_th);
          }
          mass_source = CS_MAX(mass_source, 0.);
          mass_source_oy = CS_MAX(mass_source_oy, 0.);

          cs_real_t vol_mass_source = mass_source * cell_f_vol[cell_id];
          cs_real_t vol_mass_source_oy = mass_source_oy * cell_f_vol[cell_id];

          cs_real_t vol_beta_x_ai = beta_x_ai * cell_f_vol[cell_id];

          /* Global mass source term for continuity (pressure) equation
           * Note that if rain were already considered in the bulk, then inner
           * mass transfer between liquid and vapor would disappear */

          //FIXME: Consider putting rain in the bulk
          //       --> implication on 'c' variables different from 'h'
          if (f_id == (CS_F_(p)->id)) {
            /* Warning: not multiplied by Cell volume! no addition neither */
            // FIXME: Addition needed to avoid deleting packing mass
            //        source term ?
            exp_st[cell_id] += mass_source;

            /* Saving evaporation rate for post-processing */
            evap_rate_rain[cell_id] = mass_source;
          }

          /* Water (vapor + condensate) in gas mass fraction equation
             except rain */
          else if (f_id == (CS_F_(ym_w)->id)) {
            exp_st[cell_id] += vol_mass_source * (1. - f_var[cell_id]);
            imp_st[cell_id] += vol_mass_source;
          }

          /* Rain drop mass equation (solve in drift model form) */
          else if (f_id == cfld_yp->id) {
            exp_st[cell_id] -= vol_mass_source_oy * y_rain[cell_id];
            imp_st[cell_id] += vol_mass_source_oy;
          }

          /* Humid air temperature equation */
          else if (f_id == (CS_F_(t)->id)) {
            /* Because the writing is in a non-conservative form */
            cs_real_t l_imp_st = vol_mass_source * cp_h;
            cs_real_t l_exp_st = 0.;
            if (x[cell_id] <= x_s_th) {
              /* Implicit term */
              l_imp_st += vol_beta_x_ai * (le_f * cp_h
                                           + (x_s_tl - x[cell_id]) * cp_v
                                             / (1. + x[cell_id]));
              /* Explicit term */
              l_exp_st = l_imp_st * (t_l_r[cell_id] - f_var[cell_id]);
            }
            else {
              cs_real_t coeft = le_f * cp_h;
              /* Implicit term */
              l_imp_st += vol_beta_x_ai * (coeft + (x_s_tl - x_s_th)
                                                   * cp_l / (1. + x[cell_id]));
              /* Explicit term */
              l_exp_st = vol_beta_x_ai * (coeft * t_l_r[cell_id]
                         + (x_s_tl - x_s_th) * (cp_v * t_l_r[cell_id] + hv0)
                           / (1. + x[cell_id])) - l_imp_st * f_var[cell_id];
            }
            imp_st[cell_id] += CS_MAX(l_imp_st, 0.);
            exp_st[cell_id] += l_exp_st;
          }

          /* Rain enthalpy equation (solve in drift model form)
           * NB: The actual variable being solved is y_rain x h_rain */
          else if (f_id == cfld_yh_rain->id) {
          /* Liquid temperature in Kelvin */
          cs_real_t t_l_k = t_l_r[cell_id]
                            + cs_physical_constants_celsius_to_kelvin;

          cs_real_t l_exp_st = 0.;

          /* Note: the solved variable is yl_r.hl_r so the source term associated
           * to evaporation is: Gamma/y_lr * (yl_r.h_lr) */
          /* Implicit term */
          cs_real_t l_imp_st = vol_mass_source_oy;

          l_exp_st -= l_imp_st * f_var[cell_id];

          /* Under saturated */
          if (x[cell_id] <= x_s_th) {
            /* Explicit term */
            l_exp_st -= vol_beta_x_ai * ((x_s_tl - x[cell_id])
                                         * (cp_v * t_l_k + hv0)
                                         + le_f * cp_h
                                           * (t_l_r[cell_id] - t_h[cell_id]));
          }
          /* Over saturated */
          else {
            cs_real_t coefh = le_f * cp_h;
            /* Explicit term */
            l_exp_st += vol_beta_x_ai * (coefh * (t_h[cell_id] - t_l_r[cell_id])
                                         + (x_s_tl - x_s_th) / (1. + x[cell_id])
                                           * (cp_l * t_h[cell_id]
                                              - (cp_v * t_l_r[cell_id] + hv0)));
          }
          /* Because we deal with an increment */
          exp_st[cell_id] += l_exp_st;
          imp_st[cell_id] += CS_MAX(l_imp_st, 0.);

            /* Saving thermal power for post-processing */
            if (t_l_r[cell_id] > 0.){
              thermal_power_rain[cell_id] = -(l_exp_st + l_imp_st * f_var[cell_id])
                                             / cell_f_vol[cell_id];
            }
          }

        } /* End if (y_rain > 0) */

      } /* End loop over all the cells of the domain */
    }
  } /* End evaporation model active */

  if (ct_opt->has_rain) {

    /* Generate rain from packing zones which are leaking
       ================================================== */

    cs_real_t *liq_mass_frac = CS_F_(y_l_pack)->val; /* Liquid mass fraction
                                                        in packing */
    /* Inner mass flux of liquids (in the packing) */
    cs_real_t *liq_mass_flow
      = cs_field_by_name("inner_mass_flux_y_l_packing")->val;

    for (int ict = 0; ict < *_n_ct_zones; ict++) {

      cs_ctwr_zone_t *ct = _ct_zone[ict];

      if (ct->xleak_fac > 0.0 && ct->type != CS_CTWR_INJECTION) {

        /* Rain generation source terms
           ============================ */

        for (cs_lnum_t i = 0; i < ct->n_outlet_faces; i++) {

          /* Leak face_id */
          cs_lnum_t face_id = ct->outlet_faces_ids[i];
          cs_lnum_t cell_id_leak, cell_id_rain;

          /* Convention: outlet is positive mass flux
           * Then upwind cell for liquid is i_face_cells[][0] */
          int sign = 1;
          if (liq_mass_flow[face_id] < 0) {
            sign = -1;
            cell_id_leak = i_face_cells[face_id][1];
            cell_id_rain = i_face_cells[face_id][0];
          }
          else {
            cell_id_leak = i_face_cells[face_id][0];
            cell_id_rain = i_face_cells[face_id][1];
          }

          /* Note: vol_mass_source must not be multiplied by
           * cell_f_vol[cell_id_rain]
           * because mass source computed from liq_mass_flow is
           * already in kg/s associated to the facing rain cell */
          cs_real_t vol_mass_source = ct->xleak_fac
            * liq_mass_frac[cell_id_leak] * sign * liq_mass_flow[face_id];

          /* Global bulk mass - continuity */
          //FIXME - Ignore for now because drops are not in the bulk
          if (f_id == (CS_F_(p)->id)) {
            /* Warning: not multiplied by Cell volume! */
            //exp_st[cell_id_rain] = mass_source;
          }

          /* Injected liquid mass equation for rain zones
             (solve in drift model form) */
          if (f_id == cfld_yp->id) {
            /* Because we deal with an increment */
            exp_st[cell_id_rain] += vol_mass_source;
                                    //* (1. - f_var[cell_id_rain]);
            //imp_st[cell_id_rain] += vol_mass_source;
          }
          /* Rain enthalpy */
          else if (f_id == cfld_yh_rain->id) {
            // FIXME: There should be a y_p factor in there so that
            // mass and enthalpy are compatible
            /* The transported variable is y_rain * T_rain */
            /* Since it is treated as a scalar, no multiplication by cp_l is
             * required */
            /* For temperature equation of the rain */
            exp_st[cell_id_rain] += vol_mass_source * h_l_p[cell_id_leak];
                                    //- f_var[cell_id_rain];
            //imp_st[cell_id_rain] += vol_mass_source;
          }

        } /* End of loop through outlet cells of the packing zone */

      } /* End of leaking zone test */

      /* Testing if we are in an rain injection zone */
      else if (ct->xleak_fac > 0.0 && ct->type == CS_CTWR_INJECTION) {
        const cs_lnum_t *ze_cell_ids = cs_volume_zone_by_name(ct->name)->elt_ids;
        cs_real_t inj_vol = ct->vol_f;

        for (cs_lnum_t j = 0; j < ct->n_cells; j++) {
          cs_lnum_t cell_id = ze_cell_ids[j];

          cs_real_t vol_mass_source = cell_f_vol[cell_id] / inj_vol
                                      * ct->q_l_bc * ct->xleak_fac;
          cs_real_t t_inj = ct->t_l_bc;

          /* Injected liquid mass equation for rain zones
             (solve in drift model form) */
          if (f_id == cfld_yp->id) {
            /* Because we deal with an increment */
            exp_st[cell_id] += vol_mass_source * (1. - f_var[cell_id]);
            imp_st[cell_id] += vol_mass_source;
          }
          /* Rain enthalpy */
          else if (f_id == cfld_yh_rain->id) {
            // FIXME: There should be a y_p factor in there so that
            // mass and enthalpy are compatible
            /* The transported variable is y_rain * h_rain */
            cs_real_t h_inj = cs_liq_t_to_h(t_inj);
            exp_st[cell_id] += vol_mass_source * (h_inj - f_var[cell_id]);
            imp_st[cell_id] += vol_mass_source;
          }
        }
      }
    } /* End of loop through the packing zones */

    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
      cs_lnum_t cell_id_0 = i_face_cells[face_id][0];
      cs_lnum_t cell_id_1 = i_face_cells[face_id][1];

      /* one of neigh. cells is in packing */
      if (packing_cell[cell_id_0] != -1 || packing_cell[cell_id_1] != -1) {

        /* Rain sink term in packing zones */
        //TODO : Add rain leak portion inside packing
        if (f_id == cfld_yp->id) {
          if (packing_cell[cell_id_0] != -1) {
            imp_st[cell_id_0] += CS_MAX(imasfl_r[face_id], 0.);
            exp_st[cell_id_0] -= CS_MAX(imasfl_r[face_id],0) * f_var[cell_id_0];
          }
          if (packing_cell[cell_id_1] != -1) {
            imp_st[cell_id_1] += CS_MAX(-imasfl_r[face_id], 0.);
            exp_st[cell_id_1] -= CS_MAX(-imasfl_r[face_id],0) * f_var[cell_id_1];
          }
        }

        if (f_id == cfld_yh_rain->id) {
          if (packing_cell[cell_id_0] != -1) {
            imp_st[cell_id_0] += CS_MAX(imasfl_r[face_id], 0.);
            exp_st[cell_id_0] -= CS_MAX(imasfl_r[face_id],0) * f_var[cell_id_0];
          }
          if (packing_cell[cell_id_1] != -1) {
            imp_st[cell_id_1] += CS_MAX(-imasfl_r[face_id], 0.);
            exp_st[cell_id_1] -= CS_MAX(-imasfl_r[face_id],0) * f_var[cell_id_1];
          }
        }

        /* Liquid source term in packing zones from rain */
        if (f_id == CS_F_(y_l_pack)->id) {
          if (packing_cell[cell_id_0] != -1)
          {
            exp_st[cell_id_0] += CS_MAX(imasfl_r[face_id],0) * cfld_yp->val[cell_id_0];
          }
          if (packing_cell[cell_id_1] != -1)
          {
            exp_st[cell_id_1] += CS_MAX(-imasfl_r[face_id],0) * cfld_yp->val[cell_id_1];
          }
        }

        if (f_id == CS_F_(yh_l_pack)->id) {
          if (packing_cell[cell_id_0] != -1) {
           exp_st[cell_id_0] += CS_MAX(imasfl_r[face_id],0)
                                 * cfld_yh_rain->val[cell_id_0];
          }

          if (packing_cell[cell_id_1] != -1) {
           exp_st[cell_id_1] += CS_MAX(imasfl_r[face_id],0)
                                * cfld_yh_rain->val[cell_id_1];
          }
        }
      }
    }

  } /* End of test on whether to generate rain */

  /* Source terms for rain drops velocity
   * ==================================== */
  if (ct_opt->solve_rain_velocity) {
    int class_id = 1;
    char vg_lim_name[80];
    sprintf(vg_lim_name, "vg_lim_p_%02d", class_id);

    /* Drops terminal velocity fields */
    cs_real_3_t *vg_lim_p = (cs_real_3_t *)cs_field_by_name(vg_lim_name)->val;
    cs_field_t *cfld_taup = cs_field_by_name_try("drift_tau_y_l_r");
    cs_real_t *cpro_taup = NULL;
    if (cfld_taup != NULL)
      cpro_taup = cfld_taup->val;

    /* Continuous phase drift velocity */
    cs_real_3_t *vd_c = (cs_real_3_t *)cs_field_by_name("vd_c")->val;

    /* Rain drift velocity variables */
    char f_name[80];
    sprintf(f_name, "v_p_%02d", class_id);
    cs_field_t *f_vp = cs_field_by_name(f_name);

    if (f_id == f_vp->id) {
      cs_real_33_t *_imp_st = (cs_real_33_t *)imp_st;
      cs_real_3_t *_exp_st = (cs_real_3_t *)exp_st;
      cs_real_3_t *vel = (cs_real_3_t *)CS_F_(vel)->val;
      cs_real_3_t *vp = (cs_real_3_t *)f_vp->val;

      for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {
        for (cs_lnum_t i = 0; i < 3; i++) {
          /* Explicit source term */
          _exp_st[cell_id][i] += rho_h[cell_id] * cell_f_vol[cell_id]
                              * 1. / cpro_taup[cell_id]
                              * (vel[cell_id][i] + vd_c[cell_id][i]
                                  + vg_lim_p[cell_id][i] - vp[cell_id][i]);
          /* Implicit source term: only diagonal terms */
          _imp_st[cell_id][i][i] += rho_h[cell_id] * cell_f_vol[cell_id]
                                       / cpro_taup[cell_id];
        }
      }
    }
    /* Interfacial pressure drop due to air / rain friction */
    if (f_id == (CS_F_(vel)->id)) {
      /* Casting implicit source term on a 3x3 symmetric matrix */
      cs_real_33_t *_imp_st = (cs_real_33_t *)imp_st;
      cs_real_3_t *_exp_st = (cs_real_3_t *)exp_st;

      /* Rain mass fraction field */
      cs_real_t *y_rain = (cs_real_t *)cfld_yp->val;
      /* Rain drift and velocity fields */
      sprintf(f_name, "vd_p_%02d", class_id);
      cs_real_3_t *cfld_drift
        = (cs_real_3_t *)cs_field_by_name("drift_vel_y_l_r")->val;
      cs_real_3_t *vp = (cs_real_3_t *)f_vp->val;

      /* Gravity norm */
      cs_real_t g = cs_math_3_norm(cs_glob_physical_constants->gravity);
      for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {
        if (y_rain[cell_id] > 0.){
          /* Droplet drift and absolute velocity */
          cs_real_t drift = cs_math_3_norm(cfld_drift[cell_id]);
          cs_real_t v_drop = cs_math_3_norm(vp[cell_id]);

          /* Droplet Reynolds and Eotvos number */
          cs_real_t re_d = rho_h[cell_id] * drift * droplet_diam / visc;
          cs_real_t e_o = g * droplet_diam * (rho_l - rho_h[cell_id]) / sigma;
          /* Sphere drag coefficient */
          cs_real_t cd = 0.;
          if (re_d > 0.){
            cd = (24. / re_d) * (1. + 0.15 * pow(re_d, 0.685));
          }
          /* Droplet terminal velocity */
          cs_real_t v_term = pow((4. * rho_l * droplet_diam * g)
              / (3. * cd * rho_h[cell_id]), 0.5);
          /* Droplet deformation / elongation */
          cs_real_t e_tau = 1. / (1. + 0.148 * pow(e_o, 0.85));
          //FIXME : check positivity of E
          cs_real_t E =   1. - cs_math_pow2(CS_MIN(v_drop / v_term, 1.))
                        * (1. - e_tau);

          /* Total drag coefficient for deformed drop */
          cs_real_t cd_tot = cd * (1. - 0.17185 * (1. - E)
              + 6.692 * cs_math_pow2(1. - E)
              - 6.605 * cs_math_pow3(1. - E));

          /* Air / droplets interfacial area density calculation */
          cs_real_t vol_frac_rain = y_rain[cell_id] * rho_h[cell_id] / rho_l;
          if (vol_frac_rain >= 1.0)
            vol_frac_rain = 1.0;
          cs_real_t a_i =  6.0 * vol_frac_rain * (1.0 - vol_frac_rain)
            / droplet_diam;

          /* Droplet relaxation time */
          cs_real_t tau_d = rho_l * cs_math_pow2(droplet_diam) / (18. * visc);
          /* Final head loss coefficient */
          cs_real_t k_drop = rho_l * (cd_tot * re_d / 24.) * droplet_diam * a_i
            / (6. * tau_d);
          for (int k = 0; k < 3; k++){
            _imp_st[cell_id][k][k] += -cell_f_vol[cell_id] * k_drop;
            _exp_st[cell_id][k] += cell_f_vol[cell_id] * k_drop * vp[cell_id][k];
          }
        }
      }
    }

  } /* End of solve_rain variable check */

  BFT_FREE(packing_cell);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Phase change mass source term from the evaporating liquid to the
 *        bulk, humid air.
 *
 * Careful, this is different from an injection source term, which would
 * normally be handled with a 'cs_equation_add_volume_mass_injection_' function.
 *
 * \param[out]  mass_source     Mass source term
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_bulk_mass_source_term(cs_real_t         mass_source[])
{
  cs_lnum_t n_cells_with_ghosts = cs_glob_mesh->n_cells_with_ghosts;
  /* Compute the mass exchange term */
  cs_real_t *imp_st;

  BFT_MALLOC(imp_st, n_cells_with_ghosts, cs_real_t);

  for (cs_lnum_t cell_id = 0; cell_id < n_cells_with_ghosts; cell_id++) {
    mass_source[cell_id] = 0.0;
    imp_st[cell_id] = 0.0;
  }

  /* Bulk mass source term is stored for pressure */

  cs_ctwr_source_term(CS_F_(p)->id,
                      mass_source,
                      imp_st);

  BFT_FREE(imp_st);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
