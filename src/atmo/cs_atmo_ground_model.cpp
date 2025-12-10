/*============================================================================
 * Atmospheric ground module - Build constants and variables to describe ground model
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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "atmo/cs_air_props.h"
#include "atmo/cs_atmo.h"
#include "atmo/cs_intprf.h"
#include "atmo/cs_atmo_profile_std.h"
#include "atmo/cs_atmo_1d_rad.h"

#include "base/cs_boundary_zone.h"
#include "base/cs_physical_constants.h"
#include "base/cs_field.h"
#include "base/cs_field_pointer.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "base/cs_parall.h"
#include "base/cs_volume_zone.h"
#include "bft/bft_error.h"
#include "bft/bft_printf.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "pprt/cs_physical_model.h"
#include "rayt/cs_rad_transfer.h"

/*----------------------------------------------------------------------------
*  Header for the current file
*----------------------------------------------------------------------------*/

#include "atmo/cs_atmo_ground_model.h"

/*----------------------------------------------------------------------------*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* global plant options structure */
static cs_plant_option_t  _plant_option = {
  .h_canopy = 0.6,
  .h_canopy_exch = 0.5,
  .d_leaf = 0.01,
  .leaf_area_index = 3,
  .k_ext_coef = 0.7,
  .albedo_plant = 0.22,
  .emi_plant = 0.94,
  .c_co2_air = 420.e-6,
  .f_co2_0_ground = 3.0e-6,
  .temp_ref_co2 = 298,
  .ea_ground = 53000.,
  .a_pho = 2.0,
  .g0_pho = 3.0e-4,
  .sf_psiv = 2.,
  .psif_psiv = -2.8,
  .oi = 210.0e-3,
  .koref = 256.0e-3,
  .kcref = 302.0e-6,
  .eao = 36000.0,
  .eac = 59430.0,
  .gamma_0 = 28.0e-6,
  .gamma_1 = 0.0509,
  .gamma_2 = 0.001,
  .tem_ref_pho = 293.2,
  .vlmaxref = 100.0e-6,
  .knitro = 1,
  .jlmaxref = 200.0e-6,
  .ea_vlmax = 58520.0,
  .ed_vlmax = 2.2e5,
  .s_vlmax = 700.0,
  .ea_jlmax = 37000.0,
  .ed_jlmax = 2.2e5,
  .s_jlmax = 700.0,
  .alpha_pho = 0.2,
  .theta_pho = 0.9,
  .v_gco2_incr = 1000.0,
  .v_gco2_decr = -1000.0,
  .xiv = 1.0e7,
  .r_root = 0.1e-3,
  .z_root = 1,
  .ld_root = 5.e3,
  .water_potential_model = 1,
  .psi_e = -1.30e-3,
  .psi_pow = 3.5,
  .kref = 250.0e-6,
  .canopy_mix_l = 0.08,
  .turb_prandtl = 0.9,
  .cdrag_leaf = 0.01,
  .zos = 0.001,
  .zref_plant = 0.8,
  .k_eddy_hc = 1.0,
  .tort_factor = 1.5,
  .ground_exch_l = 5.0e-3,
  .dry_ground_porosity = 0.4,
  .dv = 2.3e-5
};

/*============================================================================
 * Static global variables
 *============================================================================*/

cs_plant_option_t *cs_glob_plant_option = &_plant_option;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize ground model variables
 */
/*----------------------------------------------------------------------------*/

static void
_ground_model_variable_init(void)
{
  /* Local variables */
  cs_real_t rscp;
  cs_real_t esaini, qsaini, huini, psini;

  cs_lnum_t n_elts;
  int n_ground_cat;
  const cs_lnum_t *elt_ids;
  cs_atmo_get_ground_zone(&n_elts, &n_ground_cat, &elt_ids);

  /* Constants */
  const cs_real_t cpvcpa_value = cs_glob_air_props->cp_v / cs_glob_air_props->cp_a;
  const cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;
  const cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();
  cs_atmo_option_t *at_opt = cs_glob_atmo_option;

  /* Get pointers to field values */
  cs_field_t *ground_temperature = cs_field_by_name("ground_temperature");
  cs_field_t *ground_pot_temperature = cs_field_by_name("ground_pot_temperature");
  cs_field_t *ground_total_water = cs_field_by_name("ground_total_water");
  cs_field_t *ground_w1 = cs_field_by_name("ground_w1");
  cs_field_t *ground_w2 = cs_field_by_name("ground_w2");

  cs_real_t *bvar_ground_temp = ground_temperature->val;
  cs_real_t *bvar_potentiel_temp = ground_pot_temperature->val;
  cs_real_t *bvar_total_water = ground_total_water->val;
  cs_real_t *bvar_w1 = ground_w1->val;
  cs_real_t *bvar_w2 = ground_w2->val;

  /* Initialization of t and qv at z0 */
  if (at_opt->ground_humidity > 1.0) {
    /* If qvsini > 1, qvsini represents relative humidity in % */
    esaini = 610.78 * exp(17.2694 * at_opt->ground_surf_temp
                          / (at_opt->ground_surf_temp + tkelvi - 35.86));
    qsaini = esaini / (phys_pro->rvsra * phys_pro->p0
                       + esaini * (1.0 - phys_pro->rvsra));
    at_opt->ground_humidity = at_opt->ground_humidity * qsaini / 100.0;
  }

  /* Loop over ground models */
  for (cs_lnum_t ground_id = 0; ground_id < n_elts; ground_id++) {
    psini = phys_pro->p0;
    if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_HUMID) {
      rscp = (phys_pro->r_pg_cnst / phys_pro->cp0)
             * (1.0 + (phys_pro->rvsra - cpvcpa_value) * at_opt->ground_humidity);
    }
    else {
      rscp = (phys_pro->r_pg_cnst / phys_pro->cp0);
    }

    bvar_ground_temp[ground_id] = at_opt->ground_surf_temp;
    bvar_potentiel_temp[ground_id] = (at_opt->ground_surf_temp + tkelvi)
                                   * pow((cs_glob_atmo_constants->ps / psini),
                                         rscp);

    /*TODO A bit too restrictive; it's possible to have water in the ground
      without a humid atmosphere*/
    if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_HUMID) {
      bvar_total_water[ground_id] = at_opt->ground_humidity;
      bvar_w1[ground_id] = at_opt->ground_w1_ini;
      bvar_w2[ground_id] = at_opt->ground_w2_ini;

      /* If not initialized, compute an approximation of the
         initial water content of the top reservoir */
      if (bvar_w1[ground_id] < 1.0e-20) {
        esaini = 610.78 * exp(17.2694 * at_opt->ground_surf_temp
                              / (at_opt->ground_surf_temp + tkelvi - 35.86));
        qsaini = esaini / (phys_pro->rvsra * psini
                           + esaini * (1.0 - phys_pro->rvsra));
        huini = at_opt->ground_humidity / qsaini;
        huini = cs::min(huini, 1.0);
        bvar_w1[ground_id] = acos(1.0 - 2.0 * huini) / acos(-1.0);
      }
      /* If the deep reservoir is uninitialized,
         set the ratio w1/w2 to 1 (layer equilibrium) */
      if (bvar_w2[ground_id] < 1.0e-20) {
        bvar_w2[ground_id] = bvar_w1[ground_id];
      }
    }
    else {
      bvar_total_water[ground_id] = 0.0;
      bvar_w1[ground_id] = 0.0;
      bvar_w2[ground_id] = 0.0;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute mean values
 */
/*----------------------------------------------------------------------------*/

static void
cs_ground_compute_mean(void)
{
  /* Local variables */
  cs_real_t rugdij, rugtij, albeij, emisij;
  cs_real_t vegeij, c1wij, c2wij, csolij;
  cs_real_t r1ij, r2ij;
  cs_real_t surf_zone;
  cs_real_t groundmax[10], groundmean[10], groundmin[10];
  const char *groundname[10] = {
    "z0 dynamique", "z0 thermique", "albedo      ", "emissivite  ",
    "csol (e-6)  ", "vegetation  ", "c1w         ", "c2w         ",
    "r1          ", "r2          "
  };

  cs_lnum_t n_elts;
  int n_ground_cat;
  const cs_lnum_t *elt_ids;
  cs_lnum_t face_id;
  cs_atmo_get_ground_zone(&n_elts, &n_ground_cat, &elt_ids);

  /* Pointers to field values */
  cs_field_t *boundary_roughness = cs_field_by_name("boundary_roughness");
  cs_field_t *boundary_thermal_roughness
    = cs_field_by_name("boundary_thermal_roughness");
  cs_field_t *boundary_albedo = cs_field_by_name("boundary_albedo");
  cs_field_t *emissivity = cs_field_by_name("emissivity");
  cs_field_t *boundary_vegetation = cs_field_by_name("boundary_vegetation");
  cs_field_t *ground_water_capacity = cs_field_by_name("ground_water_capacity");
  cs_field_t *ground_water_ratio = cs_field_by_name("ground_water_ratio");
  cs_field_t *ground_thermal_capacity = cs_field_by_name("ground_thermal_capacity");
  cs_field_t *ground_r1 = cs_field_by_name("ground_r1");
  cs_field_t *ground_r2 = cs_field_by_name("ground_r2");
  cs_field_t *ground_temperature_deep = cs_field_by_name("ground_temperature_deep");
  cs_field_t *atmo_ground_percentages = cs_field_by_name("atmo_ground_percentages");

  cs_real_t *bvar_dyn_rough = boundary_roughness->val;
  cs_real_t *bvar_therm_rough = boundary_thermal_roughness->val;
  cs_real_t *bvar_albedo = boundary_albedo->val;
  cs_real_t *bvar_emissi = emissivity->val;
  cs_real_t *bvar_vegeta = boundary_vegetation->val;
  cs_real_t *bvar_water_capacity = ground_water_capacity->val;
  cs_real_t *bvar_water_ratio = ground_water_ratio->val;
  cs_real_t *bvar_thermal_capacity = ground_thermal_capacity->val;
  cs_real_t *bvar_r1 = ground_r1->val;
  cs_real_t *bvar_r2 = ground_r2->val;
  cs_real_t *bvar_temperature_deep = ground_temperature_deep->val;
  cs_real_t *bvar_ground_percentages = atmo_ground_percentages->val;

  /* Initializations */
  const cs_real_t codinv = -999.0;

  for (cs_lnum_t ground_id = 0; ground_id < n_elts ; ground_id++) {
    face_id =elt_ids[ground_id];
    bvar_dyn_rough[face_id] = codinv;
    bvar_therm_rough[face_id] = codinv;
    bvar_albedo[face_id] = codinv;
    bvar_emissi[face_id] = codinv;

    bvar_vegeta[ground_id] = codinv;
    bvar_water_capacity[ground_id] = codinv;
    bvar_water_ratio[ground_id] = codinv;
    bvar_thermal_capacity[ground_id] = codinv;
    bvar_r1[ground_id] = codinv;
    bvar_r2[ground_id] = codinv;
  }

  const cs_atmo_option_t *aopt = cs_glob_atmo_option;
  cs_lnum_t as_dim = atmo_ground_percentages->dim;

  /* Calculate coefficients for each cell */
  for (cs_lnum_t ground_id = 0; ground_id < n_elts ; ground_id++) {
    rugdij = 0.0;
    rugtij = 0.0;
    albeij = 0.0;
    emisij = 0.0;
    csolij = 0.0;
    vegeij = 0.0;
    c1wij = 0.0;
    c2wij = 0.0;
    r1ij = 0.0;
    r2ij = 0.0;

    for (cs_lnum_t n = 1; n < as_dim; n++) {
      cs_real_t fac = bvar_ground_percentages[n+ground_id*as_dim] / 100.0;
      /* shifting of 1 for the percentage because starting by default value */
      rugdij += aopt->ground_cat_roughness[n] * fac;
      rugtij += aopt->ground_cat_thermal_roughness[n] * fac;
      albeij += aopt->ground_cat_albedo[n] * fac;
      emisij += aopt->ground_cat_emissi[n] * fac;
      csolij += aopt->ground_cat_thermal_inertia[n] * fac;
      vegeij += aopt->ground_cat_vegeta[n] * fac;
      c1wij += aopt->ground_cat_w1[n] * fac;
      c2wij += aopt->ground_cat_w2[n] * fac;
      r1ij += aopt->ground_cat_r1[n] * fac;
      r2ij += aopt->ground_cat_r2[n] * fac;
    }

    face_id = elt_ids[ground_id];
    bvar_dyn_rough[face_id] = rugdij;
    bvar_therm_rough[face_id] = rugtij;
    bvar_albedo[face_id] = albeij;
    bvar_emissi[face_id] = emisij;

    bvar_thermal_capacity[ground_id] = csolij;
    bvar_vegeta[ground_id] = vegeij;
    bvar_water_capacity[ground_id] = c1wij;
    bvar_water_ratio[ground_id] = c2wij;
    bvar_r1[ground_id] = r1ij;
    bvar_r2[ground_id] = r2ij;

    /* Initialize deep temperatures to tprini */
    bvar_temperature_deep[ground_id] = aopt->ground_temperature;
  }

  /* Error checking */
  cs_gnum_t error_count = 0;
  for (cs_lnum_t ground_id = 0; ground_id < n_elts; ground_id++) {
    face_id = elt_ids[ground_id];
    if (bvar_dyn_rough[face_id] <= codinv ) error_count++;
    if (bvar_therm_rough[face_id] <= codinv) error_count++;
    if (bvar_albedo[face_id] <= codinv) error_count++;
    if (bvar_emissi[face_id] <= codinv) error_count++;

    if (bvar_thermal_capacity[ground_id] <= codinv) error_count++;
    if (bvar_vegeta[ground_id] <= codinv) error_count++;
    if (bvar_water_capacity[ground_id] <= codinv) error_count++;
    if (bvar_water_ratio[ground_id] <= codinv) error_count++;
    if (bvar_r1[ground_id] <= codinv) error_count++;
    if (bvar_r2[ground_id] <= codinv) error_count++;
  }

  cs_parall_counter(&error_count, 1);

  /* Print error message if necessary */
  if (error_count != 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("%s: incorrect initialization of coefficients of the\n"
                "ground-atmosphere interface\n"
                "For %d quantities, %llu local values are not initialized."),
              __func__, as_dim,
              (unsigned long long)error_count);
  }
  else {
    /* Initialize control variables */
    for (cs_lnum_t n = 0; n < 10; n++) {
      groundmin[n] = 999999.0;
      groundmean[n] = 0.0;
      groundmax[n] = -999999.0;
    }
    surf_zone = 0.0;

    cs_real_t *b_face_surf = cs_glob_mesh_quantities->b_face_surf;

    for (cs_lnum_t ground_id = 0; ground_id < n_elts ; ground_id++) {
      face_id = elt_ids[ground_id];
      if (bvar_dyn_rough[face_id] > groundmax[0])
        groundmax[0] = bvar_dyn_rough[face_id];
      if (bvar_therm_rough[face_id] > groundmax[1])
        groundmax[1] = bvar_therm_rough[face_id];
      if (bvar_albedo[face_id] > groundmax[2])
        groundmax[2] = bvar_albedo[face_id];
      if (bvar_emissi[face_id] > groundmax[3])
        groundmax[3] = bvar_emissi[face_id];

      if (bvar_thermal_capacity[ground_id] > groundmax[4])
        groundmax[4] = bvar_thermal_capacity[ground_id];
      if (bvar_vegeta[ground_id] > groundmax[5])
        groundmax[5] = bvar_vegeta[ground_id];
      if (bvar_water_capacity[ground_id] > groundmax[6])
        groundmax[6] = bvar_water_capacity[ground_id];
      if (bvar_water_ratio[ground_id] > groundmax[7])
        groundmax[7] = bvar_water_ratio[ground_id];
      if (bvar_r1[ground_id] > groundmax[8])
        groundmax[8] = bvar_r1[ground_id];
      if (bvar_r2[ground_id] > groundmax[9])
        groundmax[9] = bvar_r2[ground_id];

      if (bvar_dyn_rough[face_id] < groundmin[0])
        groundmin[0] = bvar_dyn_rough[face_id];
      if (bvar_therm_rough[face_id] < groundmin[1])
        groundmin[1] = bvar_therm_rough[face_id];
      if (bvar_albedo[face_id] < groundmin[2])
        groundmin[2] = bvar_albedo[face_id];
      if (bvar_emissi[face_id] < groundmin[3])
        groundmin[3] = bvar_emissi[face_id];

      if (bvar_thermal_capacity[ground_id] < groundmin[4])
        groundmin[4] = bvar_thermal_capacity[ground_id];
      if (bvar_vegeta[ground_id] < groundmin[5])
        groundmin[5] = bvar_vegeta[ground_id];
      if (bvar_water_capacity[ground_id] < groundmin[6])
        groundmin[6] = bvar_water_capacity[ground_id];
      if (bvar_water_ratio[ground_id] < groundmin[7])
        groundmin[7] = bvar_water_ratio[ground_id];
      if (bvar_r1[ground_id] < groundmin[8])
        groundmin[8] = bvar_r1[ground_id];
      if (bvar_r2[ground_id] < groundmin[9])
        groundmin[9] = bvar_r2[ground_id];

      groundmean[0] += b_face_surf[face_id] * bvar_dyn_rough[face_id];
      groundmean[1] += b_face_surf[face_id] * bvar_therm_rough[face_id];
      groundmean[2] += b_face_surf[face_id] * bvar_albedo[face_id];
      groundmean[3] += b_face_surf[face_id] * bvar_emissi[face_id];

      groundmean[4] += b_face_surf[face_id] * bvar_thermal_capacity[ground_id];
      groundmean[5] += b_face_surf[face_id] * bvar_vegeta[ground_id];
      groundmean[6] += b_face_surf[face_id] * bvar_water_capacity[ground_id];
      groundmean[7] += b_face_surf[face_id] * bvar_water_ratio[ground_id];
      groundmean[8] += b_face_surf[face_id] * bvar_r1[ground_id];
      groundmean[9] += b_face_surf[face_id] * bvar_r2[ground_id];

      /* Surface of the zone, could use it directly in C */
      surf_zone += b_face_surf[face_id];
    }

    cs_parall_sum(10, CS_DOUBLE, groundmean);
    cs_parall_min(10, CS_DOUBLE, groundmin);
    cs_parall_max(10, CS_DOUBLE, groundmax);
    cs_parall_sum(1, CS_DOUBLE, &surf_zone);

    for (cs_lnum_t n = 0; n < 10; n++) {
      groundmean[n] = groundmean[n] / surf_zone;
    }

    bft_printf(" ** ========================================= **\n"
               " ** Soil/atmosphere interface                 **\n"
               " ** Array of constants                        **\n"
               " ** ========================================= **\n");
    bft_printf(" *            * minimum*    mean* maximum*\n");
    for (cs_lnum_t n = 0; n < 10; n++) {
      if (n == 4) {
        bft_printf(" *%-12s*%8.4f*%8.4f*%8.4f*\n",
                   groundname[n], groundmin[n] * 1e6, groundmean[n] * 1e6,
                   groundmax[n] * 1e6);
      }
      else {
        bft_printf(" *%-12s*%8.4f*%8.4f*%8.4f*\n",
                   groundname[n], groundmin[n], groundmean[n], groundmax[n]);
      }
    }
    bft_printf(" ** ========================================= **\n");
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Soil Plant Atmo Continuum model - compute plant and ground convective
 * exchange resistances.
 *
 * \param[in] face_id and ground_id.
 */
/*----------------------------------------------------------------------------*/

static void
_compute_convective_exch_resistances(cs_lnum_t face_id,
                                     cs_lnum_t ground_id)
{
  cs_plant_option_t *plant_opt = cs_glob_plant_option;

  cs_field_t *u_k_atm = cs_field_by_name_try("boundary_k");
  cs_field_t *u_star_atm = cs_field_by_name_try("boundary_ustar");

  cs_field_t *plant_air_resistance = cs_field_by_name_try("plant_air_resistance");
  cs_field_t *ground_air_resistance = cs_field_by_name_try("ground_air_resistance");

  cs_real_t eta_ext_coef = plant_opt->h_canopy
    * pow(plant_opt->cdrag_leaf
    * plant_opt->leaf_area_index
    / (plant_opt->h_canopy
    * 2*pow(plant_opt->canopy_mix_l,2)),1./3.);

  cs_real_t u_star_ground = u_star_atm->val[face_id] * exp(-eta_ext_coef);;

  cs_real_t u_star_plant = pow(pow(u_star_atm->val[face_id],2)
    - pow(u_star_ground,2),0.5);

  cs_real_t canopy_exch_height = plant_opt->h_canopy
    * (1+log(1-exp(-eta_ext_coef)) / 2 / eta_ext_coef);

  cs_real_t g_plant_forced = u_star_plant * plant_opt->canopy_mix_l / canopy_exch_height / plant_opt->turb_prandtl;
  cs_real_t g_plant_nat = 0.005; /* TODO use expression in term of T_plant - T_air */

  cs_real_t g_ground_forced = u_star_ground * plant_opt->canopy_mix_l / canopy_exch_height / plant_opt->turb_prandtl;
  cs_real_t g_ground_nat = 0.005; /* TODO use expression in term of T_plant - T_air */

  /* Results */
  plant_air_resistance->val[ground_id] = 1. /(g_plant_nat + g_plant_forced) ;
  ground_air_resistance->val[ground_id] = 1. /(g_ground_nat + g_ground_forced) ;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Soil Plant Atmo Continuum model - compute water stress.
 *
 * \param[in] face_id and ground_id.
 */
/*----------------------------------------------------------------------------*/

static void
_compute_water_stress(cs_lnum_t ground_id)
{
  cs_plant_option_t *plant_opt = cs_glob_plant_option;
  cs_air_fluid_props_t *ct_prop = cs_glob_air_props;

  cs_field_t *ground_water_capacity = cs_field_by_name("ground_water_capacity");
  cs_field_t *ground_w1 = cs_field_by_name("ground_w1");
  cs_field_t *et_plant = cs_field_by_name_try("et_plant");
  cs_field_t *ground_water_potential = cs_field_by_name_try("ground_water_potential");
  cs_field_t *root_water_potential = cs_field_by_name_try("root_water_potential");
  cs_field_t *plant_water_potential = cs_field_by_name_try("plant_water_potential");
  cs_field_t *root_p_water_volumetric_capacity = cs_field_by_name_try("root_p_water_volumetric_capacity");
  cs_field_t *ground_p_water_volumetric_capacity = cs_field_by_name_try("ground_p_water_volumetric_capacity");
  cs_field_t *water_stress_factor = cs_field_by_name_try("water_stress_factor");

  /* Resistive approach */

  cs_real_t krground = 7.85e-10/(2*cs_math_pi
    * plant_opt->r_root);
  cs_real_t xis = 1/(plant_opt->z_root*2*cs_math_pi
    * plant_opt->ld_root
    * plant_opt->r_root*krground);

  cs_real_t psis = plant_opt->psi_e
    * pow(1./ground_w1->val[ground_id] , plant_opt->psi_pow);

  cs_real_t psir = psis - xis*et_plant->val[ground_id]/ct_prop->rho_l;

  cs_real_t psiv = psir - plant_opt->xiv
    * et_plant->val[ground_id]/ct_prop->rho_l;

  cs_real_t fpsiv = (1+exp(plant_opt->sf_psiv
    * plant_opt->psif_psiv))
    / (1+exp(plant_opt->sf_psiv
    * (plant_opt->psif_psiv - psiv)));

  /* Results water stress factor */
  water_stress_factor->val[ground_id] = fpsiv;

  /* Post processing */
  ground_water_potential->val[ground_id] = psis;
  root_water_potential->val[ground_id] = psir;
  plant_water_potential->val[ground_id] = psiv;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Soil Plant Atmo Continuum model - compute stomatal conductance.
 *
 * \param[in] face_id and ground_id.
 */
/*----------------------------------------------------------------------------*/

static void
_compute_stomatal_conductance(cs_lnum_t ground_id,
                              cs_real_t photo_active_rad,
                              cs_real_t air_leaf_temp,
                              cs_real_t ps,
                              cs_real_t dtref)
{
  cs_plant_option_t *plant_opt = cs_glob_plant_option;

  cs_field_t *ci_pho = cs_field_by_name_try("ci_pho");
  cs_field_t *assimilation_rate = cs_field_by_name_try("assimilation_rate");
  cs_field_t *gco2 = cs_field_by_name_try("gco2");
  cs_field_t *plant_air_resistance = cs_field_by_name_try("plant_air_resistance");
  cs_field_t *ground_air_resistance = cs_field_by_name_try("ground_air_resistance");
  cs_field_t *water_stress_factor = cs_field_by_name_try("water_stress_factor");
  cs_field_t *ground_temperature = cs_field_by_name("ground_temperature");
  cs_field_t *leaf_temp = cs_field_by_name_try("leaf_temp");

  cs_real_t f_co2_ground = plant_opt->f_co2_0_ground
    * exp(plant_opt->ea_ground
    * (ground_temperature->val[ground_id]
    + cs_physical_constants_celsius_to_kelvin
    - plant_opt->temp_ref_co2)
    / (cs_physical_constants_r*(ground_temperature->val[ground_id]
    + cs_physical_constants_celsius_to_kelvin)
    * plant_opt->temp_ref_co2));

  cs_real_t eta_ext_coef = plant_opt->h_canopy
    * pow(plant_opt->cdrag_leaf
    * plant_opt->leaf_area_index
    / (plant_opt->h_canopy
    * 2. * pow(plant_opt->canopy_mix_l,2.)),1./3.);

  cs_real_t rv = -1./plant_opt->k_eddy_hc
    * exp(eta_ext_coef)*plant_opt->h_canopy/eta_ext_coef
    * (exp(-eta_ext_coef*plant_opt->h_canopy
    / plant_opt->h_canopy)
    - exp(-eta_ext_coef*plant_opt->h_canopy_exch
    / plant_opt->h_canopy));

  cs_real_t ra = -1./plant_opt->k_eddy_hc
    * exp(eta_ext_coef)*plant_opt->h_canopy/eta_ext_coef
    * (exp(-eta_ext_coef*plant_opt->zref_plant
    / plant_opt->h_canopy)
    - exp(-eta_ext_coef*plant_opt->h_canopy
    / plant_opt->h_canopy));

  cs_real_t ci_pho_no_a = (plant_opt->c_co2_air
    + f_co2_ground*(rv+ra)
    / ps*cs_physical_constants_r
    * (leaf_temp->val[ground_id]
    + cs_physical_constants_celsius_to_kelvin));

  cs_real_t gamma_pho = plant_opt->gamma_0
    * (1. + plant_opt->gamma_1*(leaf_temp->val[ground_id]
    - plant_opt->tem_ref_pho
    + cs_physical_constants_celsius_to_kelvin)
    + plant_opt->gamma_2*pow(leaf_temp->val[ground_id]
    + cs_physical_constants_celsius_to_kelvin
    - plant_opt->tem_ref_pho,2));

  cs_real_t vlmax = plant_opt->vlmaxref
    * exp(plant_opt->ea_vlmax/cs_physical_constants_r
    / plant_opt->tem_ref_pho
    * (1. - plant_opt->tem_ref_pho/(leaf_temp->val[ground_id]
    + cs_physical_constants_celsius_to_kelvin )))
    / (1. + exp((plant_opt->s_vlmax*(leaf_temp->val[ground_id]
    + cs_physical_constants_celsius_to_kelvin )
    - plant_opt->ed_vlmax)/cs_physical_constants_r
    / (leaf_temp->val[ground_id]
    + cs_physical_constants_celsius_to_kelvin )));

  /* TODO Check the value for knitro Leuning 1995 */

  cs_real_t vcmax = plant_opt->leaf_area_index * vlmax
    * (1. - exp(-plant_opt->knitro)) / plant_opt->knitro;

  cs_real_t jlmax = plant_opt->jlmaxref
    * exp(plant_opt->ea_jlmax/cs_physical_constants_r
    / plant_opt->tem_ref_pho
    * (1. - plant_opt->tem_ref_pho/(leaf_temp->val[ground_id]
    + cs_physical_constants_celsius_to_kelvin )))
    / (1. + exp((plant_opt->s_jlmax*(leaf_temp->val[ground_id]
    + cs_physical_constants_celsius_to_kelvin )
    - plant_opt->ed_jlmax)/cs_physical_constants_r
    / (leaf_temp->val[ground_id]
    + cs_physical_constants_celsius_to_kelvin )));

  cs_real_t kc =  plant_opt->kcref
    * exp( plant_opt->eac/cs_physical_constants_r
    / plant_opt->tem_ref_pho
    * (1. -  plant_opt->tem_ref_pho
    / (leaf_temp->val[ground_id]
    + cs_physical_constants_celsius_to_kelvin )));

  cs_real_t ko =  plant_opt->koref
    * exp( plant_opt->eao/cs_physical_constants_r
    / plant_opt->tem_ref_pho
    * (1. -  plant_opt->tem_ref_pho
    / (leaf_temp->val[ground_id]
    + cs_physical_constants_celsius_to_kelvin )));

  cs_real_t discr_j = (pow(plant_opt->alpha_pho*photo_active_rad
    + jlmax,2.) - 4.*photo_active_rad*plant_opt->theta_pho
    * plant_opt->alpha_pho*jlmax);

  cs_real_t j = (plant_opt->alpha_pho*photo_active_rad + jlmax
    - pow(discr_j,0.5))/2./plant_opt->theta_pho;

  cs_real_t gco2_previous = gco2->val[ground_id];

  /* Newton algorithm resolution for photosynthesis */

  cs_real_t ci_pho_0 = 0.0003;
  cs_real_t ci_pho_1 = 0.0002;
  cs_real_t f_ci_pho_0 = 0.;
  cs_real_t f_ci_pho_1 = 0.;
  cs_real_t assimilation_rate_0 = 0.;
  cs_real_t assimilation_rate_1 = 0.;
  cs_real_t gco2_0 = 0.;
  cs_real_t gco2_1 = 0.;
  cs_real_t vc_0 = 0.;
  cs_real_t vc_1 = 0.;
  cs_real_t vj_0 = 0.;
  cs_real_t vj_1 = 0.;
  cs_real_t delta_ci = 0.;
  int i;
  for (i=0; i<30; i++){

    vc_0 = vcmax*(ci_pho_0
      - gamma_pho)/(ci_pho_0 + kc
      * (1 + plant_opt->oi/ko));
    vc_1 = vcmax*(ci_pho_1
      - gamma_pho)/(ci_pho_1 + kc
      * (1 + plant_opt->oi/ko));

    vj_0 = j/4.*(ci_pho_0 - gamma_pho)
      / (ci_pho_0+ 2.*gamma_pho);
    vj_1 = j/4.*(ci_pho_1 - gamma_pho)
      / (ci_pho_1 + 2.*gamma_pho);

    assimilation_rate_0 =cs::max(0., fmin(vc_0 , vj_0));
    assimilation_rate_1 =cs::max(0., fmin(vc_1 , vj_1));

    gco2_0 =  plant_opt->g0_pho
      + plant_opt->a_pho*assimilation_rate_0
      / (ci_pho_0 - gamma_pho)* water_stress_factor->val[ground_id];
    gco2_0 = cs::max( plant_opt->g0_pho
      , gco2_0);
    gco2_1 =  plant_opt->g0_pho
      + plant_opt->a_pho*assimilation_rate_1
      / (ci_pho_1 - gamma_pho)* water_stress_factor->val[ground_id];
    gco2_1 = cs::max( plant_opt->g0_pho
      , gco2_1);

    f_ci_pho_0 = ci_pho_0 - plant_opt->c_co2_air
      - f_co2_ground*(ra + rv) / ps*cs_physical_constants_r
      * (air_leaf_temp + cs_physical_constants_celsius_to_kelvin)
      + assimilation_rate_0*((rv+ra+1.6*plant_air_resistance->val[ground_id])
      / ps*cs_physical_constants_r
      * (air_leaf_temp+ cs_physical_constants_celsius_to_kelvin)
      + 1./gco2_0);

    f_ci_pho_1 = ci_pho_1 - plant_opt->c_co2_air
      - f_co2_ground*(ra + rv) / ps*cs_physical_constants_r
      * (air_leaf_temp + cs_physical_constants_celsius_to_kelvin)
      + assimilation_rate_1*((rv+ra+1.6*plant_air_resistance->val[ground_id])
      / ps*cs_physical_constants_r
      * (air_leaf_temp + cs_physical_constants_celsius_to_kelvin)
      + 1./gco2_1);

    if (( f_ci_pho_0 == f_ci_pho_1) && (f_ci_pho_1 != 0.)){
      ci_pho->val[ground_id] = (ci_pho_0 + ci_pho_1)/2;
    }
    else if (( f_ci_pho_0 == f_ci_pho_1) && (f_ci_pho_1 == 0.)){
      ci_pho->val[ground_id] = ci_pho_1;
    }
    else if (f_ci_pho_0 == 0. ){
      ci_pho->val[ground_id] = ci_pho_0;
    }
    else if (f_ci_pho_1 == 0. ){
      ci_pho->val[ground_id] = ci_pho_1;
    }
    else {
      ci_pho->val[ground_id] = ci_pho_1 - (ci_pho_1-ci_pho_0)
        / (f_ci_pho_1-f_ci_pho_0)*f_ci_pho_1;
    }
    if(( ci_pho->val[ground_id] < 0.) || (ci_pho->val[ground_id] > 0.0006)){
      bft_printf("PROBLEM ci equal to %f\n", ci_pho->val[ground_id]);
      ci_pho->val[ground_id] = (ci_pho_0 + ci_pho_1)/2;
    }
    ci_pho_0 = ci_pho_1;
    ci_pho_1 = ci_pho->val[ground_id];
    delta_ci = 1000*(ci_pho->val[ground_id] - ci_pho_1);
  }

  if (ci_pho->val[ground_id] == 0){
    bft_printf("PROBLEM ci can't be equal to 0 \n");
  }
  else if (delta_ci/ci_pho->val[ground_id]>0.2) {
    bft_printf("PROBLEM NEW ALOGRITHM HAS NOT CONVERGED \n");
  }

  cs_real_t vc = vcmax*(ci_pho->val[ground_id]
    - gamma_pho)/(ci_pho->val[ground_id] + kc
    * (1. + plant_opt->oi/ko));

  cs_real_t vj = j/4.*(ci_pho->val[ground_id] - gamma_pho)
    / (ci_pho->val[ground_id] + 2.*gamma_pho);

  /* Assumption photorespiration is more important
   * than the other processes */
  assimilation_rate->val[ground_id] = cs::max(0., cs::min(vc , vj));
  assimilation_rate->val[ground_id] = cs::max(0.,vj);

  gco2->val[ground_id] =  plant_opt->g0_pho
    + plant_opt->a_pho*assimilation_rate->val[ground_id]
    / (ci_pho->val[ground_id] - gamma_pho)* water_stress_factor->val[ground_id];

  gco2->val[ground_id] = cs::max(plant_opt->g0_pho , gco2->val[ground_id]);

  /* To account for plant time response delay to light fluctuations*/
  cs_real_t gco2_max = gco2_previous + plant_opt->v_gco2_incr*dtref;
  cs_real_t gco2_min = gco2_previous + plant_opt->v_gco2_decr*dtref;
  gco2->val[ground_id] = cs::min(gco2->val[ground_id],gco2_max);
  gco2->val[ground_id] = cs::max(gco2->val[ground_id],gco2_min);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Soil Plant Atmo Continuum model - compute energy exchanges and plant temperature.
 *
 * \param[in] face_id and ground_id.
 */
/*----------------------------------------------------------------------------*/

static void
_compute_le_h_and_leaf_temp(cs_lnum_t ground_id,
                            cs_real_t rsh2o,
                            cs_real_t air_leaf_temp,
                            cs_real_t pvsat_dewpoint,
                            cs_real_t wv_molar_mass,
                            cs_real_t rho_moist_air,
                            cs_real_t cp_moist_air,
                            cs_real_t rn_thermal_atm_plant,
                            cs_real_t rn_thermal_ground_plant)
{
  cs_plant_option_t *plant_opt = cs_glob_plant_option;
  cs_air_fluid_props_t *ct_prop = cs_glob_air_props;

  const cs_real_t stephn = cs_physical_constants_stephan;

  cs_field_t *cover_geometry_ratio
    = cs_field_by_name_try("cover_geometry_ratio");
  cs_field_t *gco2 = cs_field_by_name_try("gco2");
  cs_field_t *plant_air_resistance
    = cs_field_by_name_try("plant_air_resistance");
  cs_field_t *leaf_temp = cs_field_by_name_try("leaf_temp");
  cs_field_t *et_plant = cs_field_by_name_try("et_plant");
  cs_field_t *le_plant = cs_field_by_name_try("latent_heat_plant_to_air");
  cs_field_t *hs_plant = cs_field_by_name_try("sensible_heat_plant_to_air");
  cs_field_t *plant_water_potential
    = cs_field_by_name_try("plant_water_potential");

  cs_field_t *ray_net_plant = cs_field_by_name_try("ray_net_plant");
  cs_field_t *ray_net_ir_plant = cs_field_by_name_try("ray_net_ir_plant");
  cs_field_t *ray_net_solar_plant
    = cs_field_by_name_try("ray_net_solar_plant");
  cs_field_t *ray_ir_plant_to_ground
    = cs_field_by_name_try("ray_ir_plant_to_ground");
  cs_field_t *ray_ir_atm_to_plant
    = cs_field_by_name_try("ray_ir_atm_to_plant");

  cs_real_t leaf_temp_0 = leaf_temp->val[ground_id]+2.;
  cs_real_t leaf_temp_1 = leaf_temp->val[ground_id]-2.;

  cs_real_t pvsat_plant_0 = 0.;
  cs_real_t pvsat_plant_1 = 0.;

  cs_real_t hi = exp(wv_molar_mass * plant_water_potential->val[ground_id]
    / (ct_prop->rho_l*cs_physical_constants_r
    * (air_leaf_temp + cs_physical_constants_celsius_to_kelvin)));

  /* Init of the Newton algo */
  cs_real_t rn_thermal_plant_0 = 0.;
  cs_real_t rn_plant_0 = 0.;
  cs_real_t rn_thermal_plant_1 = 0.;
  cs_real_t rn_plant_1 = 0.;
  cs_real_t et_plant_0 = 0.;
  cs_real_t et_plant_1 = 0.;
  cs_real_t le_plant_0 = 0.;
  cs_real_t le_plant_1 = 0.;
  cs_real_t hs_plant_0 = 0.;
  cs_real_t hs_plant_1 = 0.;
  cs_real_t f_leaf_temp_0 = 0.;
  cs_real_t f_leaf_temp_1 = 0.;
  cs_real_t delta_leaf_temp = 0.;

  int nb;

  for (nb=0; nb<30; nb++){

    pvsat_plant_0 = cs_air_pwv_sat(leaf_temp_0);
    pvsat_plant_1 = cs_air_pwv_sat(leaf_temp_1);

    rn_thermal_plant_0 = 2.*plant_opt->emi_plant*stephn
      * cs_math_pow4(leaf_temp_0
      + cs_physical_constants_celsius_to_kelvin)
      * cover_geometry_ratio->val[ground_id];
    rn_plant_0 = ray_net_solar_plant->val[ground_id] + rn_thermal_atm_plant
      - rn_thermal_plant_0 + rn_thermal_ground_plant;

    rn_thermal_plant_1 = 2.*plant_opt->emi_plant*stephn
      * cs_math_pow4(leaf_temp_1
      + cs_physical_constants_celsius_to_kelvin)
      * cover_geometry_ratio->val[ground_id];
    rn_plant_1 = ray_net_solar_plant->val[ground_id] + rn_thermal_atm_plant
      - rn_thermal_plant_1 + rn_thermal_ground_plant;

    et_plant_0 = cs::max(0.,wv_molar_mass/cs_physical_constants_r
      * (1./(rsh2o + plant_air_resistance->val[ground_id]))
      * (hi*pvsat_plant_0/(leaf_temp_0
      + cs_physical_constants_celsius_to_kelvin)
      - pvsat_dewpoint/(air_leaf_temp
      + cs_physical_constants_celsius_to_kelvin)));
    et_plant_1 = cs::max(0.,wv_molar_mass/cs_physical_constants_r
      * (1./(rsh2o + plant_air_resistance->val[ground_id]))
      * (hi*pvsat_plant_1/(leaf_temp_1
      + cs_physical_constants_celsius_to_kelvin)
      - pvsat_dewpoint/(air_leaf_temp
      + cs_physical_constants_celsius_to_kelvin)));
    le_plant_0 = ct_prop->hv0*et_plant_0;
    le_plant_1 = ct_prop->hv0*et_plant_1;

    hs_plant_0 = rho_moist_air*cp_moist_air/plant_air_resistance->val[ground_id]
      *(leaf_temp_0 - air_leaf_temp);
    hs_plant_1 = rho_moist_air*cp_moist_air/plant_air_resistance->val[ground_id]
      *(leaf_temp_1 - air_leaf_temp);

    f_leaf_temp_0 = rn_plant_0 - le_plant_0 - hs_plant_0;
    f_leaf_temp_1 = rn_plant_1 - le_plant_1 - hs_plant_1;

    if (( f_leaf_temp_0 == f_leaf_temp_1) && (f_leaf_temp_1 != 0.)){
      leaf_temp->val[ground_id] = (leaf_temp_0 + leaf_temp_1)/2.;
    }
    else if (( f_leaf_temp_0 == f_leaf_temp_1) && (f_leaf_temp_1 == 0.)){
      leaf_temp->val[ground_id] = leaf_temp_1;
    }
    else if (f_leaf_temp_0 == 0. ){
      leaf_temp->val[ground_id] = leaf_temp_0;
    }
    else if (f_leaf_temp_1 == 0. ){
      leaf_temp->val[ground_id] = leaf_temp_1;
    }
    else {
      leaf_temp->val[ground_id] = leaf_temp_1 - (leaf_temp_1-leaf_temp_0)
        / (f_leaf_temp_1-f_leaf_temp_0)*f_leaf_temp_1;
    }
    if(( leaf_temp->val[ground_id] < -20.) || (leaf_temp->val[ground_id] > 80.)){
      bft_printf("PROBLEM leaf temp equal to %f\n", leaf_temp->val[ground_id]);
      leaf_temp->val[ground_id] = (leaf_temp_0 + leaf_temp_1)/2.;
    }
    leaf_temp_0 = leaf_temp_1;
    leaf_temp_1 = leaf_temp->val[ground_id];
    delta_leaf_temp = 1000.*(leaf_temp->val[ground_id] - leaf_temp_1);
  }

  if (delta_leaf_temp/(0.1 + CS_ABS(leaf_temp->val[ground_id]))>0.2) {
    bft_printf("PROBLEM NEW ALOGRITHM HAS NOT CONVERGED \n");
  }
  cs_real_t pvsat_plant = cs_air_pwv_sat(leaf_temp->val[ground_id]);

  et_plant->val[ground_id] = cs::max(0.,wv_molar_mass/cs_physical_constants_r
    * (1./(rsh2o + plant_air_resistance->val[ground_id]))
    * (hi*pvsat_plant/(leaf_temp->val[ground_id]
    + cs_physical_constants_celsius_to_kelvin)
    - pvsat_dewpoint/(air_leaf_temp
    + cs_physical_constants_celsius_to_kelvin)));

  cs_real_t rn_thermal_plant = 2.*plant_opt->emi_plant*stephn
    * cs_math_pow4(leaf_temp->val[ground_id]
    + cs_physical_constants_celsius_to_kelvin)
    * cover_geometry_ratio->val[ground_id];
  cs_real_t rn_plant = ray_net_solar_plant->val[ground_id] + rn_thermal_atm_plant
    - rn_thermal_plant + rn_thermal_ground_plant;

  le_plant->val[ground_id] = ct_prop->hv0*et_plant->val[ground_id];

  hs_plant->val[ground_id] = rho_moist_air*cp_moist_air
    / plant_air_resistance->val[ground_id]
    *(leaf_temp->val[ground_id] - air_leaf_temp);

  ray_net_plant->val[ground_id] = rn_plant;
  ray_net_ir_plant->val[ground_id] = rn_thermal_atm_plant
    - rn_thermal_plant + rn_thermal_ground_plant;
  ray_ir_plant_to_ground->val[ground_id] = plant_opt->emi_plant*stephn
    * cs_math_pow4(leaf_temp->val[ground_id]
    + cs_physical_constants_celsius_to_kelvin)
    * cover_geometry_ratio->val[ground_id];
  ray_ir_atm_to_plant->val[ground_id] = rn_thermal_atm_plant;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Soil Plant Atmo Continuum model - compute source term plant and ground
 * to air.
 *
 * \param[in] face_id and ground_id.
 */
/*----------------------------------------------------------------------------*/

static void
_compute_plant_to_air_source_terms(cs_lnum_t ground_id,
                                   cs_real_t dtref,
                                   cs_real_t dz_max,
                                   cs_real_t rho_moist_air)
{

  cs_air_fluid_props_t *ct_prop = cs_glob_air_props;

  cs_field_t *leaf_temp = cs_field_by_name_try("leaf_temp");
  cs_field_t *et_plant = cs_field_by_name_try("et_plant");
  cs_field_t *sensible_heat_plant_to_air
    = cs_field_by_name_try("sensible_heat_plant_to_air");
  cs_field_t *st_rho_p = cs_field_by_name_try("source_term_vapor_mass_plant");
  cs_field_t *st_q_p_exp
    = cs_field_by_name_try("source_term_exp_specific_humidity_plant");
  cs_field_t *st_q_p_imp
    = cs_field_by_name_try("source_term_imp_specific_humidity_plant");
  cs_field_t *st_conv_energy_p
    = cs_field_by_name_try("source_term_convective_energy_plant");
  cs_field_t *st_mass_energy_p
    = cs_field_by_name_try("source_term_mass_energy_plant");

  cs_real_t drho_p_to_atm  = et_plant->val[ground_id] * dtref / dz_max;

  st_rho_p->val[ground_id] = drho_p_to_atm / dtref;

  cs_real_t dq_p_to_atm = drho_p_to_atm / rho_moist_air ;

  st_q_p_exp->val[ground_id] = dq_p_to_atm / (1. + dq_p_to_atm) / dtref;
  st_q_p_imp->val[ground_id] = - dq_p_to_atm / (1. + dq_p_to_atm) / dtref;

  st_conv_energy_p->val[ground_id] = sensible_heat_plant_to_air->val[ground_id] / dz_max ;
  st_mass_energy_p->val[ground_id] = drho_p_to_atm
    * ct_prop->cp_v * ( leaf_temp->val[ground_id]
    + cs_physical_constants_celsius_to_kelvin ) / dtref;

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Deardorff force restore model
 */
/*----------------------------------------------------------------------------*/

void
cs_ground_model(void)
{
  cs_atmo_option_t *at_opt = cs_glob_atmo_option;
  cs_atmo_1d_rad_t *at_1d_rad = cs_glob_atmo_1d_rad;

  int z_id = at_opt->ground_zone_id;
  if (z_id > -1) {
    int micro_scale_option = at_opt->ground_meb_model;

    cs_rad_transfer_params_t *rt_params = cs_glob_rad_transfer_params;
    cs_plant_option_t *plant_opt = cs_glob_plant_option;

    const cs_real_t stephn = cs_physical_constants_stephan;
    const cs_zone_t *z = cs_boundary_zone_by_id(z_id);
    const cs_lnum_t *elt_ids = z->elt_ids;
    const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;
    cs_mesh_quantities_t *fvq   = cs_glob_mesh_quantities;
    const cs_real_3_t *cell_cen = (const cs_real_3_t *)fvq->cell_cen;
    cs_real_t *dt = cs_field_by_name("dt")->val;
    /* Post treatment fields */
    cs_field_t *ground_sensible_heat
      = cs_field_by_name_try("ground_sensible_heat");
    cs_field_t *ground_latent_heat
      = cs_field_by_name_try("ground_latent_heat");
    cs_field_t *ground_thermal_rad_upward
      = cs_field_by_name_try("ground_thermal_rad_upward");
    cs_field_t *ground_thermal_rad_downward
      = cs_field_by_name_try("ground_thermal_rad_downward");
    cs_field_t *ground_visible_rad_absorbed
      = cs_field_by_name_try("ground_visible_rad_absorbed");

    /* Soil related fields */
    cs_field_t *ground_temperature = cs_field_by_name("ground_temperature");
    cs_field_t *ground_pot_temperature = cs_field_by_name("ground_pot_temperature");
    cs_field_t *ground_total_water = cs_field_by_name("ground_total_water");
    cs_field_t *ground_w1 = cs_field_by_name("ground_w1");
    cs_field_t *ground_w2 = cs_field_by_name("ground_w2");
    cs_real_t *f_fos = cs_field_by_name("ground_solar_incident_flux")->val;
    cs_real_t *f_foir = cs_field_by_name("ground_infrared_incident_flux")->val;
    cs_field_t *ground_temperature_deep
      = cs_field_by_name("ground_temperature_deep");
    cs_field_t *ground_r1 = cs_field_by_name("ground_r1");
    cs_field_t *ground_r2 = cs_field_by_name("ground_r2");
    cs_field_t *ground_water_capacity = cs_field_by_name("ground_water_capacity");
    cs_field_t *ground_water_ratio = cs_field_by_name("ground_water_ratio");
    cs_field_t *ground_thermal_capacity
      = cs_field_by_name("ground_thermal_capacity");
    cs_field_t *ground_percentages = cs_field_by_name("atmo_ground_percentages");
    cs_field_t *boundary_vegetation = cs_field_by_name("boundary_vegetation");
    /* Fields related to all faces */
    cs_field_t *boundary_albedo = cs_field_by_name_try("boundary_albedo");
    cs_field_t *emissivity = cs_field_by_name_try("emissivity");
    /* Cell fields  */
    cs_field_t *density_moist_air = cs_field_by_name_try("density");
    cs_field_t *atm_total_water = cs_field_by_name_try("ym_water");
    cs_field_t *atm_temp = CS_F_(t);
    cs_field_t *meteo_pressure = cs_field_by_name_try("meteo_pressure");
    /* Radiative tables */
    cs_real_t *sold = (cs_real_t *)at_1d_rad->sold;
    cs_real_t *ird = (cs_real_t *) at_1d_rad->ird;
    /* Pointer to the spectral flux density field */
    cs_field_t *f_qinspe = nullptr;
    if (rt_params->atmo_model != CS_RAD_ATMO_3D_NONE)
      f_qinspe = cs_field_by_name_try("spectral_rad_incident_flux");

    /* Exchange coefficients*/
    cs_real_t *h_t = CS_F_(t)->bc_coeffs->bf;
    cs_real_t *h_q = atm_total_water->bc_coeffs->bf;

    const cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();
    /* Deardorff parameterisation */
    const cs_real_t tau_1 = 86400.;

    cs_air_fluid_props_t *ct_prop = cs_glob_air_props;

    /* Update previous values */
    cs_field_current_to_previous(ground_pot_temperature);
    cs_field_current_to_previous(ground_total_water);

    /* In case of multi energy balance (MEB) models including PV */
    cs_field_t *cover_geometry_ratio
      = cs_field_by_name_try("cover_geometry_ratio");
    cs_field_t *cover_reflectivity
      = cs_field_by_name_try("cover_reflectivity");
    cs_field_t *cover_temperature_radiative
      = cs_field_by_name_try("cover_temperature_radiative");

    /* Plant related fields */
    cs_field_t *ray_net_plant = cs_field_by_name_try("ray_net_plant");
    cs_field_t *ray_net_ir_plant = cs_field_by_name_try("ray_net_ir_plant");
    cs_field_t *ray_net_solar_plant = cs_field_by_name_try("ray_net_solar_plant");
    cs_field_t *ray_ir_plant_to_ground = cs_field_by_name_try("ray_ir_plant_to_ground");
    cs_field_t *ray_ir_atm_to_plant = cs_field_by_name_try("ray_ir_atm_to_plant");
    cs_field_t *et_plant = cs_field_by_name_try("et_plant");
    cs_field_t *et_ground = cs_field_by_name_try("et_ground");
    cs_field_t *assimilation_rate = cs_field_by_name_try("assimilation_rate");
    cs_field_t *gco2 = cs_field_by_name_try("gco2");
    cs_field_t *ground_w1_transpirated_by_the_plant = cs_field_by_name_try("ground_w1_transpirated_by_the_plant");
    cs_field_t *ground_w1_evaporated_by_the_ground = cs_field_by_name_try("ground_w1_evaporated_by_the_ground");
    cs_field_t *q_emitted_by_the_plant = cs_field_by_name_try("q_emitted_by_the_plant");
    cs_field_t *leaf_temp = cs_field_by_name_try("leaf_temp");
    cs_field_t *sensible_heat_plant_to_air = cs_field_by_name_try("sensible_heat_plant_to_air");
    cs_field_t *root_p_water_volumetric_capacity = cs_field_by_name_try("root_p_water_volumetric_capacity");
    cs_field_t *ground_p_water_volumetric_capacity = cs_field_by_name_try("ground_p_water_volumetric_capacity");
    cs_field_t *water_stress_factor = cs_field_by_name_try("water_stress_factor");
    cs_field_t *sensible_heat_plant = cs_field_by_name_try("sensible_heat_plant");
    cs_field_t *latent_heat_plant = cs_field_by_name_try("latent_heat_plant");
    cs_field_t *plant_air_resistance = cs_field_by_name_try("plant_air_resistance");
    cs_field_t *ground_air_resistance = cs_field_by_name_try("ground_air_resistance");

    cs_field_t *st_rho_s = cs_field_by_name_try("source_term_vapor_mass_ground");
    cs_field_t *st_conv_energy_s = cs_field_by_name_try("source_term_convective_energy_ground");
    cs_field_t *st_mass_energy_s = cs_field_by_name_try("source_term_mass_energy_ground");

    /*
    plant_opt->h_canopy = cs_notebook_parameter_value_by_name("h_canopy");
    plant_opt->h_canopy_exch = cs_notebook_parameter_value_by_name("h_canopy_exch");
    plant_opt->leaf_area_index = cs_notebook_parameter_value_by_name("leaf_area_index");
    plant_opt->albedo_plant = cs_notebook_parameter_value_by_name("albedo_plant");

    plant_opt->a_pho = cs_notebook_parameter_value_by_name("a_pho");
    plant_opt->sf_psiv = cs_notebook_parameter_value_by_name("sf_psiv");
    plant_opt->psif_psiv = cs_notebook_parameter_value_by_name("psif_psiv");
    plant_opt->canopy_mix_l = cs_notebook_parameter_value_by_name("canopy_mix_l");
    */

    for (cs_lnum_t ground_id = 0; ground_id < z->n_elts; ground_id++) {
      cs_real_t ray2 = 0;
      cs_real_t chas2 = 0;
      cs_real_t chal2 = 0;
      cs_real_t rapp2 = 0;
      cs_real_t secmem = 0.;
      cs_real_t premem = 0.;
      cs_real_t pphy = 0.;
      cs_real_t dum = 0.;
      cs_real_t qvs_new = 0.;
      cs_real_t w1_min = 0.;
      cs_real_t w1_max = 1.;
      cs_real_t w2_min = 0.;
      cs_real_t w2_max = 1.;
      cs_real_t precip = 0.;
      cs_real_t transpire = 0.;
      cs_real_t tseuil = 16. + cs_physical_constants_celsius_to_kelvin;
      cs_real_t ts_new = 0.;
      cs_real_t ts_c_new = 0.; /* in Celsius */
      cs_real_t cpvcpa = ct_prop->cp_v / ct_prop->cp_a;
      cs_real_t rvsra = phys_pro->rvsra;
      cs_real_t clatev = phys_pro->clatev;
      cs_real_t cp0 = phys_pro->cp0;
      cs_real_t rair = phys_pro->r_pg_cnst;
      cs_real_t wv_molar_mass = phys_pro->xmasmr*ct_prop->molmass_rat;
      cs_real_t ps = cs_glob_atmo_constants->ps;
      cs_real_t esat = 0.;

      cs_lnum_t face_id = elt_ids[ground_id];
      cs_lnum_t cell_id = b_face_cells[face_id];
      cs_real_t foir = 0.;
      cs_real_t fos = 0.;

      cs_real_t dtref  = dt[cell_id];
      cs_real_t dz_max = cell_cen[cell_id][2]*2.;

      cs_real_t h_t_ground = h_t[face_id];
      cs_real_t h_q_ground = h_q[face_id];
      cs_real_t albedo_s = cs_glob_atmo_option->ground_cat_albedo[3];
      cs_real_t emi_s = cs_glob_atmo_option->ground_cat_emissi[3];
      cs_real_t rho_moist_air = density_moist_air->val[cell_id];
      cs_real_t cp_moist_air = cp0*(1+0.856*atm_total_water->val[cell_id]);
      cs_real_t coef_flux_to_hum = dtref / rho_moist_air
        / dz_max;
      cs_real_t gcr = 0.;
      cs_real_t refl = 0.;
      if (cover_geometry_ratio != nullptr)
        gcr = cover_geometry_ratio->val[ground_id];
      if (cover_reflectivity != nullptr)
        refl = cover_reflectivity->val[ground_id];

      /* Infrared and Solar radiative fluxes
       * Warning: should be adapted for many verticales */
      if (at_1d_rad->radiative_model_1d == 1
         && rt_params->atmo_model == CS_RAD_ATMO_3D_NONE) {
        foir = ird[0];
        fos = sold[0];
      }
      /* In case of 3D, take the incident flux */
      else if (rt_params->atmo_model != CS_RAD_ATMO_3D_NONE) {

        /* Number of bands (stride) */
        cs_lnum_t stride = rt_params->nwsgg;

        /* Compute fos */
        int gg_id = rt_params->atmo_dr_id;
        if (gg_id > -1)
          fos += f_qinspe->val[gg_id + face_id * stride];

        gg_id = rt_params->atmo_dr_o3_id;
        if (gg_id > -1)
          fos += f_qinspe->val[gg_id + face_id * stride];

        gg_id = rt_params->atmo_df_id;
        if (gg_id > -1)
          fos += f_qinspe->val[gg_id + face_id * stride];

        gg_id = rt_params->atmo_df_o3_id;
        if (gg_id > -1)
          fos += f_qinspe->val[gg_id + face_id * stride];

        /* Compute foir */
        gg_id = rt_params->atmo_ir_id;
        if (gg_id > -1)
          foir += f_qinspe->val[gg_id + face_id * stride];
      }

      /* f_fos and f_foir store the radiations coming ontop of
       * the boundary cell */
      f_fos[ground_id] = fos;
      f_foir[ground_id] = foir;

      f_fos[ground_id] = 400.;
      f_foir[ground_id] = 300.;

      if (at_opt->meteo_profile == 0) {
        cs_atmo_profile_std(0., /* z_ref */
                            phys_pro->p0,
                            phys_pro->t0,
                            cell_cen[cell_id][2], &pphy, &dum, &dum);
      }
      else if (at_opt->meteo_profile == 1) {
        int nbmett = at_opt->met_1d_nlevels_t;
        int nbmetm = at_opt->met_1d_ntimes;
        pphy = cs_intprf(nbmett,
                         nbmetm,
                         at_opt->z_temp_met,
                         at_opt->time_met,
                         at_opt->hyd_p_met,
                         cell_cen[cell_id][2],
                         cs_glob_time_step->t_cur);
      }
      else {
        pphy = meteo_pressure->val[cell_id];
      }

      /* ====================================
       * Specific case for water faces (sea/lake)
       * ==================================== */

      /* Water is second component of ground_percentages */
      cs_lnum_t cat_id = 1 + ground_percentages->dim * ground_id;
      if (ground_percentages->val[cat_id] > 50.) {
        /* NB: ground_temperature in °C */
        esat = cs_air_pwv_sat(ground_temperature->val[ground_id]);
        qvs_new = esat / ( rvsra * pphy
          + esat * (1.0 - rvsra));
        /* Constant temperature at the moment
         * TODO integrate the lake model GLM to match with LNHE models */
        ts_c_new = ground_temperature->val[ground_id];
        ts_new = ts_c_new + cs_physical_constants_celsius_to_kelvin;
      }
      else {

      cs_real_t rscp1 = (rair / cp0) * (1. + (rvsra - cpvcpa)
          * ground_total_water->val[ground_id] );

      /* ====================================
       * Specific case for vegetation
       * ==================================== */

      if (micro_scale_option == CS_ATMO_SOIL_VEGETATION) {

        /* Check if latent heat of vaporization of water (J/kg) is defined */
        if (ct_prop->hv0 < 100000.){
          ct_prop->hv0 = 2260000.;
        }
        /* Check if density of liquid water (kg/m3) is defined */
        if (ct_prop->rho_l < 900.){
          ct_prop->rho_l = 1000.;
        }

        /* ============================
         * Parameters for the plant layer
         * ============================ */

        cover_geometry_ratio->val[ground_id] = (1-exp(
          - plant_opt->k_ext_coef
          * plant_opt->leaf_area_index));

        cover_reflectivity->val[ground_id] = 1 - plant_opt->emi_plant;

        /* ============================
         * Radiation on the canopy
         * ============================ */

        ray_net_solar_plant->val[ground_id] = (1.-plant_opt->albedo_plant)
          * fos*cover_geometry_ratio->val[ground_id];

        /* Note that the plant acts as a black body for PAR */
        cs_real_t photo_active_rad = ray_net_solar_plant->val[ground_id]
          / (1 - plant_opt->albedo_plant) *1.e-6/0.48;

        cs_real_t rn_thermal_atm_plant = cs::max(plant_opt->emi_plant*foir
          * cover_geometry_ratio->val[ground_id],50.);

        cs_real_t rn_thermal_ground_plant = emissivity->val[face_id] * stephn
          * cs_math_pow4(ground_temperature->val[ground_id]
          + cs_physical_constants_celsius_to_kelvin)
          * cover_geometry_ratio->val[ground_id];

        /* ============================
         * Air close to the canopy
         * ============================ */

        /* Assumption no liquid water */
        cs_real_t rho_vapor = rho_moist_air*atm_total_water->val[cell_id];

        cs_real_t air_leaf_temp = atm_temp->val[cell_id]
          - cs_physical_constants_celsius_to_kelvin;
        if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == 2){
          air_leaf_temp = atm_temp->val[cell_id] * pow(pphy / ps , rscp1)
            - cs_physical_constants_celsius_to_kelvin;
        }
        else if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == 1){
          air_leaf_temp = atm_temp->val[cell_id] * pow(pphy / ps , rscp1
            / (1. + (rvsra - cpvcpa)* ground_total_water->val[ground_id] ))
            - cs_physical_constants_celsius_to_kelvin;
        }

        /* ============================
         * Exchange resistances
         * ============================ */

        _compute_convective_exch_resistances(face_id,
                                             ground_id);

        h_t_ground = rho_moist_air * cp0 / ground_air_resistance->val[ground_id];
        h_q_ground = rho_moist_air / ground_air_resistance->val[ground_id];

        /* ============================
         * Saturation vapor pressure
         * ============================ */

        cs_real_t pvsat_airleaf = cs_air_pwv_sat(air_leaf_temp);

        cs_real_t pv_airleaf = rho_vapor*rair*(air_leaf_temp
          + cs_physical_constants_celsius_to_kelvin)/ct_prop->molmass_rat;

        cs_real_t h_rel = pv_airleaf/pvsat_airleaf;
        cs_real_t t_dewpoint = 237.7*(17.27*(air_leaf_temp)
          / (237.7 + air_leaf_temp) + log(h_rel))
          / (17.27 - 17.27*air_leaf_temp
          / (237.7 + air_leaf_temp) + log(h_rel));

        cs_real_t pvsat_dewpoint = cs_air_pwv_sat(t_dewpoint);

        cs_real_t pvsat_ground = cs_air_pwv_sat(ground_temperature->val[ground_id]);

        /* ============================
         * Water potential
         * ============================ */

        _compute_water_stress(ground_id);

        /* ============================
         * Stomatal conductance
         * ============================ */

        _compute_stomatal_conductance(ground_id,
                                      photo_active_rad,
                                      air_leaf_temp,
                                      ps,
                                      dtref);

        cs_real_t rsh2o = ps/1.6/gco2->val[ground_id]/cs_physical_constants_r
          / (air_leaf_temp + cs_physical_constants_celsius_to_kelvin);

        /* ============================
         * Energy exchange and leaf temperature
         * ============================ */

        _compute_le_h_and_leaf_temp(ground_id,
                                    rsh2o,
                                    air_leaf_temp,
                                    pvsat_dewpoint,
                                    wv_molar_mass,
                                    rho_moist_air,
                                    cp_moist_air,
                                    rn_thermal_atm_plant,
                                    rn_thermal_ground_plant);

        cover_temperature_radiative->val[ground_id] = leaf_temp->val[ground_id];

        /* ============================
         * Compute water exchanges
         * ============================ */

        transpire = et_plant->val[ground_id];

        et_ground->val[ground_id] = cs::max(0. , wv_molar_mass
          / cs_physical_constants_r*(1/ground_air_resistance->val[ground_id])
          *(pvsat_ground/(ground_temperature->val[ground_id]
          + cs_physical_constants_celsius_to_kelvin)
          - pvsat_dewpoint/(air_leaf_temp
          + cs_physical_constants_celsius_to_kelvin)));

        ground_w1_transpirated_by_the_plant->val[ground_id] = dtref * transpire;
        ground_w1_evaporated_by_the_ground->val[ground_id] = dtref
          * et_ground->val[ground_id];   /* To compare with the ground model*/

        q_emitted_by_the_plant->val[ground_id] = dtref * transpire
          / plant_opt->h_canopy / rho_moist_air;

        /* ============================
         * Compute source terms plant to air
         * ============================ */

        _compute_plant_to_air_source_terms(ground_id,
                                           dtref,
                                           dz_max,
                                           rho_moist_air);

      }

      /* ===============================
       * Compute coefficients for heat and latent heat fluxes
       * =============================== */

      /* ratio specific heat of humide air/ specidif heat of dry air
       * Cph/Cpd */
      cs_real_t cph_dcpd = (1. + (cpvcpa - 1.)
        * ground_total_water->val[ground_id] );
      /* Conversion theta -> T */
      cs_real_t cht = h_t[face_id] * pow(ps / pphy, rscp1) * cph_dcpd;

      cs_real_t chq = h_q[face_id]
        * (clatev - 2370.* (ground_temperature->val[ground_id]) );

      /* ===============================
       * Compute reservoirs water (shallow and deep)
       * =============================== */

      cs_real_t evapor = - h_q[face_id] * (atm_total_water->val[cell_id]
        - ground_total_water->val[ground_id]);

      cs_real_t w1_num = ground_w1->val[ground_id]
        + dtref * (precip - evapor - transpire)
        / ground_water_capacity->val[ground_id]
        + ground_w2->val[ground_id] * dtref
          / (tau_1 + ground_water_ratio->val[ground_id] * dtref);
      cs_real_t w1_den = 1.
        + 1. / (tau_1/dtref + ground_water_ratio->val[ground_id]);
      cs_real_t w1_new = w1_num / w1_den;
      w1_new = cs::max(w1_new, w1_min);
      w1_new = cs::min(w1_new, w1_max);
      cs_real_t w2_num = ground_w2->val[ground_id] * tau_1
        + w1_new * dtref * ground_water_ratio->val[ground_id];
      cs_real_t w2_den = tau_1 + dtref * ground_water_ratio->val[ground_id];
      cs_real_t w2_new = w2_num / w2_den;
      w2_new = cs::max(w2_new, w2_min);
      w2_new = cs::min(w2_new, w2_max);

      ground_w1->val[ground_id] = w1_new;
      ground_w2->val[ground_id] = w2_new;

      cs_real_t hu = 0.5 * (1. - cos(cs_math_pi * w1_new));

      /* ============================
       * Compute saturated pressure and DL1
       * ============================ */

      esat = cs_air_pwv_sat(ground_temperature->val[ground_id]);
      cs_real_t rapsat = rvsra * pphy + esat * (1. - rvsra);
      cs_real_t qsat = esat / rapsat;
      cs_real_t cstder = 17.438 * 239.78; //#TODO transfer in cs_air_props.c
      cs_real_t dqsat = pphy * rvsra * cstder * esat
        / cs_math_pow2(rapsat * (ground_temperature->val[ground_id] + 239.78));

      /* Compute equivalent emissivity and reflexivity du to pannels/plants */
      cs_real_t c_refl = 1. - refl;
      cs_real_t emi = emissivity->val[face_id];
      cs_real_t emi_eq = gcr * emi / (emi * refl + c_refl);
      cs_real_t c_gcr = 1. - gcr;
      /* ============================
       * Compute the first member of ground_temperature equation
       * ============================ */

      cs_lnum_t iseuil = 0;

      /* Soil temp is in Celius.*/
      if (ground_temperature->val[ground_id]
          + cs_physical_constants_celsius_to_kelvin < tseuil)
        iseuil = 1;

      cs_lnum_t ichal = 1;

      cs_real_t ray1 = 4. * c_gcr * emissivity->val[face_id] * stephn
        * cs_math_pow3(ground_temperature->val[ground_id]
        + cs_physical_constants_celsius_to_kelvin);
      cs_real_t chas1 = cht;
      cs_real_t chal1 = chq * hu * dqsat;
      cs_real_t rapp1 = 2.* cs_math_pi / tau_1;

      premem = ground_thermal_capacity->val[ground_id] * ( ray1
        + chas1 + chal1 * ichal + ground_r2->val[ground_id] * iseuil ) + rapp1;

      /* ============================
       * Compute the second member of ground_temperature equation
       * ============================ */

      ray2 = c_gcr * fos*(1. - boundary_albedo->val[face_id])
        + emi * c_gcr * foir
        + 3. * (emi_eq * c_refl  + c_gcr) * emi * stephn
        * cs_math_pow4(ground_temperature->val[ground_id]
        + cs_physical_constants_celsius_to_kelvin);

      if ((micro_scale_option == CS_ATMO_SOIL_PHOTOVOLTAICS) ||
          (micro_scale_option == CS_ATMO_SOIL_VEGETATION)) {
        ray2 += cs_math_pow4(cover_temperature_radiative->val[ground_id]
            + cs_physical_constants_celsius_to_kelvin) * stephn
          * emi_eq * c_refl * refl;
      }
      chas2 = cht * atm_temp->val[cell_id]
        * pow(pphy / ps , rscp1);
      chal2 = chq * (atm_total_water->val[cell_id]
        * (1. - boundary_vegetation->val[ground_id] * ( 1. - hu))
        - hu * (qsat - (ground_temperature->val[ground_id]
              + cs_physical_constants_celsius_to_kelvin)
              * dqsat));
      rapp2 = 2.*cs_math_pi * (ground_temperature_deep->val[ground_id]
        + cs_physical_constants_celsius_to_kelvin) / tau_1;

      secmem = ground_thermal_capacity->val[ground_id] * ( ray2
        + chas2 + chal2 * ichal + ground_r1->val[ground_id]
        + tseuil * ground_r2->val[ground_id] * iseuil ) + rapp2;

      /* ============================
       * Compute new ground variables
       * ============================ */

      ts_new = (ground_temperature->val[ground_id]
        + cs_physical_constants_celsius_to_kelvin + dtref * secmem )
        / ( 1. + dtref * premem );
      ts_c_new = (ground_temperature->val[ground_id]
        + dtref * (secmem - cs_physical_constants_celsius_to_kelvin * premem))
        / ( 1. + dtref * premem );

      qvs_new = hu * ( qsat + dqsat * ( ts_c_new
        - ground_temperature->val[ground_id] ))
        + boundary_vegetation->val[ground_id] * atm_total_water->val[cell_id]
        * ( 1. - hu );

      /* TODO filling ground fields should be done one for loop below
       * by computing cht and chq for the ""lake model""
       * At the moment the allocation is performed ONLY for ""ground model""  */

      if (ground_latent_heat != nullptr)
        ground_latent_heat->val[ground_id] = chq * (qvs_new
            - atm_total_water->val[cell_id]);
      if (ground_sensible_heat != nullptr)
        ground_sensible_heat->val[ground_id] = cht * ((ts_new
              * (pow(ps/pphy,(rair/cp0)
                  * (1. + (rvsra - cpvcpa) * qvs_new ))))
            - atm_temp->val[cell_id]);
      if (ground_thermal_rad_upward != nullptr)
        ground_thermal_rad_upward->val[ground_id] = stephn * emi
          * cs_math_pow4(ts_new);
      if (ground_thermal_rad_downward != nullptr)
        ground_thermal_rad_downward->val[ground_id] = emi * foir;
      if (ground_visible_rad_absorbed != nullptr)
        ground_visible_rad_absorbed->val[ground_id] =
          (1. - boundary_albedo->val[face_id]) * fos;

      }

      /* ============================
       * Update new ground variables
       * ============================ */

      ground_temperature->val[ground_id] = ts_c_new;

      ground_pot_temperature->val[ground_id] = ts_new
        * (pow(ps/pphy,(rair/cp0)
        * (1. + (rvsra - cpvcpa) * qvs_new )));
      ground_total_water->val[ground_id] = qvs_new;
    }
  }
  else {
    bft_error(__FILE__, __LINE__, 0,
              _("at_opt->ground_zone_id is missing."));
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build constants and variables to describe ground model
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_ground_initialize(void)
{
  /* Local variables */
  cs_lnum_t n_elts;
  int n_ground_cat;
  const cs_lnum_t *elt_ids;
  cs_atmo_get_ground_zone(&n_elts, &n_ground_cat, &elt_ids);

  /* Get the number of elements in the ground zone */
  const cs_atmo_option_t *at_opt = cs_glob_atmo_option;

  cs_parall_sum(1, CS_INT_TYPE,&n_elts);

  /* There are some ground faces on some ranks */
  /* Note: we can use ground categories without ground model */
  /* (which solve temperature and humidity) */
  if (n_elts > 0) {
    /* Second pass, print and check ground categories parameters */
    cs_atmo_ground_cat(2);
    cs_ground_compute_mean();

    /* Initialization of ground variables */
    /* Only if ground is activated */
    if (at_opt->ground_model >= 0) {
      _ground_model_variable_init();
    }
  } /* End of second call */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Soil - atmosphere parameters computed from a "Land use" file
 *
 * \param[in] call_stage  first pass to set default values,
 *                        second pass to perform some checks and log
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_ground_cat(int call_stage)
{
  /* Local variables */
  int water, forest, diverse, mineral, diffu, mixed, dense, building;
  const int inityp = -9;
  const double codinv = -999.0;
  char nomcat[8][10];

  cs_atmo_option_t *at_opt = cs_glob_atmo_option;
  cs_lnum_t n_elts;
  int n_ground_cat;
  const cs_lnum_t *elt_ids;

  cs_atmo_get_ground_zone(&n_elts, &n_ground_cat, &elt_ids);

  /* Initialization */

  water = inityp;
  forest = inityp;
  diverse = inityp;
  mineral = inityp;
  diffu = inityp;
  mixed = inityp;
  dense = inityp;
  building = inityp;

  /* Case of a file based on IGN data in 7 categories */
  if (n_ground_cat == 7) {
    water = 1;
    forest = 2;
    diverse = 3;
    mineral = 4;
    diffu = 5;
    mixed = 6;
    dense = 7;
  }
  /* Case of a file based on IGN data in 5 categories */
  else if (n_ground_cat == 5) {
    water = 1;
    forest = 2;
    diverse = 3;
    mineral = 4;
    building = 5;
  }

  /* Soil categories names */
  if (water != inityp) strncpy(nomcat[water], "water", 6);
  if (forest != inityp) strncpy(nomcat[forest], "forest", 7);
  if (diverse != inityp) strncpy(nomcat[diverse], "diverse", 8);
  if (mineral != inityp) strncpy(nomcat[mineral], "mineral", 8);
  if (diffu != inityp) strncpy(nomcat[diffu], "bg diffu", 9);
  if (mixed != inityp) strncpy(nomcat[mixed], "bg mixed", 9);
  if (dense != inityp) strncpy(nomcat[dense], "bg dense", 9);
  if (building != inityp) strncpy(nomcat[building], "building", 9);

  /* First pass, default values according to the choice of number of grounds */
  /*----------------------------------------------------------------------*/

  if (call_stage == 1) {
    cs_atmo_ground_init_arrays(&n_ground_cat,
                             &at_opt->ground_cat_thermal_inertia,
                             &at_opt->ground_cat_roughness,
                             &at_opt->ground_cat_thermal_roughness,
                             &at_opt->ground_cat_albedo,
                             &at_opt->ground_cat_emissi,
                             &at_opt->ground_cat_vegeta,
                             &at_opt->ground_cat_w1,
                             &at_opt->ground_cat_w2,
                             &at_opt->ground_cat_r1,
                             &at_opt->ground_cat_r2);

    for (cs_lnum_t n = 1; n <= n_ground_cat; n++) {
      at_opt->ground_cat_thermal_inertia[n] = codinv;
      at_opt->ground_cat_roughness[n] = codinv;
      at_opt->ground_cat_thermal_roughness[n] = codinv;
      at_opt->ground_cat_albedo[n] = codinv;
      at_opt->ground_cat_emissi[n] = codinv;
      at_opt->ground_cat_vegeta[n] = codinv;
      at_opt->ground_cat_w1[n] = codinv;
      at_opt->ground_cat_w2[n] = codinv;
      at_opt->ground_cat_r1[n] = codinv;
      at_opt->ground_cat_r2[n] = codinv;
    }

    /* Standard values of parameters */
    if (water != inityp) at_opt->ground_cat_roughness[water] = 0.0005;
    if (forest != inityp) at_opt->ground_cat_roughness[forest] = 0.800;
    if (diverse != inityp) at_opt->ground_cat_roughness[diverse] = 0.100;
    if (mineral != inityp) at_opt->ground_cat_roughness[mineral] = 0.0012;
    if (diffu != inityp) at_opt->ground_cat_roughness[diffu] = 0.250;
    if (mixed != inityp) at_opt->ground_cat_roughness[mixed] = 0.600;
    if (dense != inityp) at_opt->ground_cat_roughness[dense] = 1.000;
    if (building != inityp) at_opt->ground_cat_roughness[building] = 0.600;

    if (water != inityp) {
      at_opt->ground_cat_thermal_roughness[water]
        = at_opt->ground_cat_roughness[water];
    }
    if (forest != inityp) {
      at_opt->ground_cat_thermal_roughness[forest]
        = at_opt->ground_cat_roughness[forest] * exp(-2.0);
    }
    if (diverse != inityp) {
      at_opt->ground_cat_thermal_roughness[diverse]
        = at_opt->ground_cat_roughness[diverse] * exp(-2.0);
    }
    if (mineral != inityp) {
      at_opt->ground_cat_thermal_roughness[mineral]
        = at_opt->ground_cat_roughness[mineral] * exp(-2.0);
    }
    if (diffu != inityp) {
      at_opt->ground_cat_thermal_roughness[diffu]
        = at_opt->ground_cat_roughness[diffu] * exp(-2.0);
    }
    if (mixed != inityp) {
      at_opt->ground_cat_thermal_roughness[mixed]
        = at_opt->ground_cat_roughness[mixed] * exp(-2.0);
    }
    if (dense != inityp) {
      at_opt->ground_cat_thermal_roughness[dense]
        = at_opt->ground_cat_roughness[dense] * exp(-2.0);
    }
    if (building != inityp) {
      at_opt->ground_cat_thermal_roughness[building]
        = at_opt->ground_cat_roughness[building] * exp(-2.0);
    }

    if (water != inityp) at_opt->ground_cat_albedo[water] = 0.08;
    if (forest != inityp) at_opt->ground_cat_albedo[forest] = 0.16;
    if (diverse != inityp) at_opt->ground_cat_albedo[diverse] = 0.20;
    if (mineral != inityp) at_opt->ground_cat_albedo[mineral] = 0.25;
    if (diffu != inityp) at_opt->ground_cat_albedo[diffu] = 0.18;
    if (mixed != inityp) at_opt->ground_cat_albedo[mixed] = 0.18;
    if (dense != inityp) at_opt->ground_cat_albedo[dense] = 0.18;
    if (building != inityp) at_opt->ground_cat_albedo[building] = 0.18;

    if (water != inityp) at_opt->ground_cat_emissi[water] = 0.980;
    if (forest != inityp) at_opt->ground_cat_emissi[forest] = 0.950;
    if (diverse != inityp) at_opt->ground_cat_emissi[diverse] = 0.940;
    if (mineral != inityp) at_opt->ground_cat_emissi[mineral] = 0.965;
    if (diffu != inityp) at_opt->ground_cat_emissi[diffu] = 0.880;
    if (mixed != inityp) at_opt->ground_cat_emissi[mixed] = 0.880;
    if (dense != inityp) at_opt->ground_cat_emissi[dense] = 0.880;
    if (building != inityp) at_opt->ground_cat_emissi[building] = 0.880;

    if (water != inityp) at_opt->ground_cat_vegeta[water] = 0.00;
    if (forest != inityp) at_opt->ground_cat_vegeta[forest] = 1.00;
    if (diverse != inityp) at_opt->ground_cat_vegeta[diverse] = 1.00;
    if (mineral != inityp) at_opt->ground_cat_vegeta[mineral] = 0.00;
    if (diffu != inityp) at_opt->ground_cat_vegeta[diffu] = 0.50;
    if (mixed != inityp) at_opt->ground_cat_vegeta[mixed] = 0.25;
    if (dense != inityp) at_opt->ground_cat_vegeta[dense] = 0.00;
    if (building != inityp) at_opt->ground_cat_vegeta[building] = 0.25;

    if (water != inityp) at_opt->ground_cat_thermal_inertia[water] = 7.6e-06;
    if (forest != inityp) at_opt->ground_cat_thermal_inertia[forest] = 11.0e-06;
    if (diverse != inityp) at_opt->ground_cat_thermal_inertia[diverse] = 11.0e-06;
    if (mineral != inityp) at_opt->ground_cat_thermal_inertia[mineral] = 5.0e-06;
    if (dense != inityp) at_opt->ground_cat_thermal_inertia[dense] = 3.9e-06;
    if (diffu != inityp) {
      at_opt->ground_cat_thermal_inertia[diffu]
        =      at_opt->ground_cat_thermal_inertia[forest]
             * at_opt->ground_cat_vegeta[diffu]
           +   at_opt->ground_cat_thermal_inertia[dense]
             * (1.0 - at_opt->ground_cat_vegeta[diffu]);
    }
    if (mixed != inityp) {
      at_opt->ground_cat_thermal_inertia[mixed]
        =     at_opt->ground_cat_thermal_inertia[forest]
            * at_opt->ground_cat_vegeta[mixed]
          +   at_opt->ground_cat_thermal_inertia[dense]
            * (1.0 - at_opt->ground_cat_vegeta[mixed]);
    }
    if (building != inityp) {
      at_opt->ground_cat_thermal_inertia[building]
        =     at_opt->ground_cat_thermal_inertia[forest]
            * at_opt->ground_cat_vegeta[building]
          +   3.9e-06 * (1.0 - at_opt->ground_cat_vegeta[building]);
    }

    if (water != inityp) at_opt->ground_cat_w1[water] = 100.0;
    if (forest != inityp) {
      at_opt->ground_cat_w1[forest] = 18.0 * at_opt->ground_cat_vegeta[forest] + 2.0;
    }
    if (diverse != inityp) {
      at_opt->ground_cat_w1[diverse]
        = 18.0 * at_opt->ground_cat_vegeta[diverse] + 2.0;
    }
    if (mineral != inityp) {
      at_opt->ground_cat_w1[mineral]
        = 18.0 * at_opt->ground_cat_vegeta[mineral] + 2.0;
    }
    if (diffu != inityp) {
      at_opt->ground_cat_w1[diffu] = 18.0 * at_opt->ground_cat_vegeta[diffu] + 2.0;
    }
    if (mixed != inityp) {
      at_opt->ground_cat_w1[mixed] = 18.0 * at_opt->ground_cat_vegeta[mixed] + 2.0;
    }
    if (dense != inityp) {
      at_opt->ground_cat_w1[dense] = 18.0 * at_opt->ground_cat_vegeta[dense] + 2.0;
    }
    if (building != inityp) {
      at_opt->ground_cat_w1[building]
        = 18.0 * at_opt->ground_cat_vegeta[building] + 2.0;
    }

    if (water != inityp) at_opt->ground_cat_w2[water] = 1.00;
    if (forest != inityp) at_opt->ground_cat_w2[forest] = 0.20;
    if (diverse != inityp) at_opt->ground_cat_w2[diverse] = 0.20;
    if (mineral != inityp) at_opt->ground_cat_w2[mineral] = 0.20;
    if (diffu != inityp) at_opt->ground_cat_w2[diffu] = 0.20;
    if (mixed != inityp) at_opt->ground_cat_w2[mixed] = 0.20;
    if (dense != inityp) at_opt->ground_cat_w2[dense] = 0.20;
    if (building != inityp) at_opt->ground_cat_w2[building] = 0.20;

    if (water != inityp) at_opt->ground_cat_r1[water] = 0.0;
    if (forest != inityp) at_opt->ground_cat_r1[forest] = 0.0;
    if (diverse != inityp) at_opt->ground_cat_r1[diverse] = 0.0;
    if (mineral != inityp) at_opt->ground_cat_r1[mineral] = 0.0;
    if (dense != inityp) at_opt->ground_cat_r1[dense] = 30.0;
    if (diffu != inityp) at_opt->ground_cat_r1[diffu] = 10.0;
    if (mixed != inityp) at_opt->ground_cat_r1[mixed] = 15.0;
    if (building != inityp) at_opt->ground_cat_r1[building] = 15.0;

    if (water != inityp) at_opt->ground_cat_r2[water] = 0.0;
    if (forest != inityp) at_opt->ground_cat_r2[forest] = 0.0;
    if (diverse != inityp) at_opt->ground_cat_r2[diverse] = 0.0;
    if (mineral != inityp) at_opt->ground_cat_r2[mineral] = 0.0;
    if (dense != inityp) at_opt->ground_cat_r2[dense] = 2.0;
    if (diffu != inityp) at_opt->ground_cat_r2[diffu] = 2.0 / 3.0;
    if (mixed != inityp) at_opt->ground_cat_r2[mixed] = 1.0;
    if (building != inityp) at_opt->ground_cat_r2[building] = 1.0;
  }

  /* Second pass: log and checks
     --------------------------- */

  if (call_stage == 2) {
    /* Log */
    bft_printf(" Soil-atmosphere interface model\n");
    bft_printf(" Values of tabulated coefficients\n");
    bft_printf(" --------------------------------\n");
    bft_printf(" Name      z0 dyn   z0 th    albedo   emissi   "
               "Cs(e-6)  vegeta   c1w      c2w      r1       r2\n");

    for (cs_lnum_t n = 1; n <= n_ground_cat; n++) {
      bft_printf
        ("%s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
         nomcat[n], at_opt->ground_cat_roughness[n],
         at_opt->ground_cat_thermal_roughness[n],
        at_opt->ground_cat_albedo[n],
         at_opt->ground_cat_emissi[n],
         1.e+06 * at_opt->ground_cat_thermal_inertia[n],
        at_opt->ground_cat_vegeta[n],
         at_opt->ground_cat_w1[n], at_opt->ground_cat_w2[n],
         at_opt->ground_cat_r1[n], at_opt->ground_cat_r2[n]);
    }
    bft_printf(" --------------------------------\n\n");

    /* Check */
    int n_errors = n_ground_cat;

    for (cs_lnum_t n = 1; n <= n_ground_cat; n++) {
      if (   at_opt->ground_cat_roughness[n] > codinv
          && at_opt->ground_cat_thermal_roughness[n] > codinv
          && at_opt->ground_cat_albedo[n] > codinv
          && at_opt->ground_cat_emissi[n] > codinv
          && at_opt->ground_cat_w1[n] > codinv
          && at_opt->ground_cat_w2[n] > codinv
          && at_opt->ground_cat_thermal_inertia[n] > codinv
          && at_opt->ground_cat_r1[n] > codinv
          && at_opt->ground_cat_r2[n] > codinv) {
        n_errors--;
      }
    }

    /* Error message if needed */
    if (n_errors != 0) {
      bft_error(__FILE__, __LINE__, 0,
                _("%s: number of grounds with incorrect ground coefficients = %d."),
                __func__, n_errors);
    }
  }
}

/*----------------------------------------------------------------------------*/
