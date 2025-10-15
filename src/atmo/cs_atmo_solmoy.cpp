/*============================================================================
 * Atmospheric soil module - Initialize ground level parameters from land use
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

#include "base/cs_physical_constants.h"
#include "atmo/cs_atmo.h"
#include "base/cs_field.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "atmo/cs_air_props.h"
#include "base/cs_math.h"
#include "bft/bft_error.h"
#include "bft/bft_printf.h"
#include "base/cs_parall.h"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "atmo/cs_atmo_solmoy.h"

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize ground level parameters from land use
 */
/*----------------------------------------------------------------------------*/

void
cs_solmoy(void)
{
  /* Local variables */
  cs_real_t rugdij, rugtij, albeij, emisij;
  cs_real_t vegeij, c1wij, c2wij, csolij;
  cs_real_t r1ij, r2ij;
  cs_real_t surf_zone;
  cs_real_t soilmax[10], soilmean[10], soilmin[10];
  const char *soilname[10] = {
    "z0 dynamique", "z0 thermique", "albedo      ", "emissivite  ",
    "csol (e-6)  ", "vegetation  ", "c1w         ", "c2w         ",
    "r1          ", "r2          "
  };

  cs_lnum_t n_elts;
  int n_soil_cat;
  const cs_lnum_t *elt_ids;
  cs_lnum_t face_id;
  cs_f_atmo_get_soil_zone(&n_elts, &n_soil_cat, &elt_ids);

  /* Pointers to field values */
  cs_field_t *boundary_roughness = cs_field_by_name("boundary_roughness");
  cs_field_t *boundary_thermal_roughness
    = cs_field_by_name("boundary_thermal_roughness");
  cs_field_t *boundary_albedo = cs_field_by_name("boundary_albedo");
  cs_field_t *emissivity = cs_field_by_name("emissivity");
  cs_field_t *boundary_vegetation = cs_field_by_name("boundary_vegetation");
  cs_field_t *soil_water_capacity = cs_field_by_name("soil_water_capacity");
  cs_field_t *soil_water_ratio = cs_field_by_name("soil_water_ratio");
  cs_field_t *soil_thermal_capacity = cs_field_by_name("soil_thermal_capacity");
  cs_field_t *soil_r1 = cs_field_by_name("soil_r1");
  cs_field_t *soil_r2 = cs_field_by_name("soil_r2");
  cs_field_t *soil_temperature_deep = cs_field_by_name("soil_temperature_deep");
  cs_field_t *atmo_soil_percentages = cs_field_by_name("atmo_soil_percentages");

  cs_real_t *bvar_dyn_rough = boundary_roughness->val;
  cs_real_t *bvar_therm_rough = boundary_thermal_roughness->val;
  cs_real_t *bvar_albedo = boundary_albedo->val;
  cs_real_t *bvar_emissi = emissivity->val;
  cs_real_t *bvar_vegeta = boundary_vegetation->val;
  cs_real_t *bvar_water_capacity = soil_water_capacity->val;
  cs_real_t *bvar_water_ratio = soil_water_ratio->val;
  cs_real_t *bvar_thermal_capacity = soil_thermal_capacity->val;
  cs_real_t *bvar_r1 = soil_r1->val;
  cs_real_t *bvar_r2 = soil_r2->val;
  cs_real_t *bvar_temperature_deep = soil_temperature_deep->val;
  cs_real_t *bvar_soil_percentages = atmo_soil_percentages->val;

  /* Initializations */
  const cs_real_t codinv = -999.0;

  for (cs_lnum_t soil_id = 0; soil_id < n_elts ; soil_id++) {
    face_id =elt_ids[soil_id];
    bvar_dyn_rough[face_id] = codinv;
    bvar_therm_rough[face_id] = codinv;
    bvar_albedo[face_id] = codinv;
    bvar_emissi[face_id] = codinv;

    bvar_vegeta[soil_id] = codinv;
    bvar_water_capacity[soil_id] = codinv;
    bvar_water_ratio[soil_id] = codinv;
    bvar_thermal_capacity[soil_id] = codinv;
    bvar_r1[soil_id] = codinv;
    bvar_r2[soil_id] = codinv;
  }

  const cs_atmo_option_t *aopt = cs_glob_atmo_option;
  cs_lnum_t as_dim = atmo_soil_percentages->dim;

  /* Calculate coefficients for each cell */
  for (cs_lnum_t soil_id = 0; soil_id < n_elts ; soil_id++) {
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
      cs_real_t fac = bvar_soil_percentages[n+soil_id*as_dim] / 100.0;
      /* shifting of 1 for the percentage because starting by default value */
      rugdij += aopt->soil_cat_roughness[n] * fac;
      rugtij += aopt->soil_cat_thermal_roughness[n] * fac;
      albeij += aopt->soil_cat_albedo[n] * fac;
      emisij += aopt->soil_cat_emissi[n] * fac;
      csolij += aopt->soil_cat_thermal_inertia[n] * fac;
      vegeij += aopt->soil_cat_vegeta[n] * fac;
      c1wij += aopt->soil_cat_w1[n] * fac;
      c2wij += aopt->soil_cat_w2[n] * fac;
      r1ij += aopt->soil_cat_r1[n] * fac;
      r2ij += aopt->soil_cat_r2[n] * fac;
    }

    face_id = elt_ids[soil_id];
    bvar_dyn_rough[face_id] = rugdij;
    bvar_therm_rough[face_id] = rugtij;
    bvar_albedo[face_id] = albeij;
    bvar_emissi[face_id] = emisij;

    bvar_thermal_capacity[soil_id] = csolij;
    bvar_vegeta[soil_id] = vegeij;
    bvar_water_capacity[soil_id] = c1wij;
    bvar_water_ratio[soil_id] = c2wij;
    bvar_r1[soil_id] = r1ij;
    bvar_r2[soil_id] = r2ij;

    /* Initialize deep temperatures to tprini */
    bvar_temperature_deep[soil_id] = aopt->soil_temperature;
  }

  /* Error checking */
  cs_gnum_t error_count = 0;
  for (cs_lnum_t soil_id = 0; soil_id < n_elts; soil_id++) {
    face_id = elt_ids[soil_id];
    if (bvar_dyn_rough[face_id] <= codinv ) error_count++;
    if (bvar_therm_rough[face_id] <= codinv) error_count++;
    if (bvar_albedo[face_id] <= codinv) error_count++;
    if (bvar_emissi[face_id] <= codinv) error_count++;

    if (bvar_thermal_capacity[soil_id] <= codinv) error_count++;
    if (bvar_vegeta[soil_id] <= codinv) error_count++;
    if (bvar_water_capacity[soil_id] <= codinv) error_count++;
    if (bvar_water_ratio[soil_id] <= codinv) error_count++;
    if (bvar_r1[soil_id] <= codinv) error_count++;
    if (bvar_r2[soil_id] <= codinv) error_count++;
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
      soilmin[n] = 999999.0;
      soilmean[n] = 0.0;
      soilmax[n] = -999999.0;
    }
    surf_zone = 0.0;

    cs_real_t *b_face_surf = cs_glob_mesh_quantities->b_face_surf;

    for (cs_lnum_t soil_id = 0; soil_id < n_elts ; soil_id++) {
      face_id = elt_ids[soil_id];
      if (bvar_dyn_rough[face_id] > soilmax[0])
        soilmax[0] = bvar_dyn_rough[face_id];
      if (bvar_therm_rough[face_id] > soilmax[1])
        soilmax[1] = bvar_therm_rough[face_id];
      if (bvar_albedo[face_id] > soilmax[2])
        soilmax[2] = bvar_albedo[face_id];
      if (bvar_emissi[face_id] > soilmax[3])
        soilmax[3] = bvar_emissi[face_id];

      if (bvar_thermal_capacity[soil_id] > soilmax[4])
        soilmax[4] = bvar_thermal_capacity[soil_id];
      if (bvar_vegeta[soil_id] > soilmax[5])
        soilmax[5] = bvar_vegeta[soil_id];
      if (bvar_water_capacity[soil_id] > soilmax[6])
        soilmax[6] = bvar_water_capacity[soil_id];
      if (bvar_water_ratio[soil_id] > soilmax[7])
        soilmax[7] = bvar_water_ratio[soil_id];
      if (bvar_r1[soil_id] > soilmax[8])
        soilmax[8] = bvar_r1[soil_id];
      if (bvar_r2[soil_id] > soilmax[9])
        soilmax[9] = bvar_r2[soil_id];

      if (bvar_dyn_rough[face_id] < soilmin[0])
        soilmin[0] = bvar_dyn_rough[face_id];
      if (bvar_therm_rough[face_id] < soilmin[1])
        soilmin[1] = bvar_therm_rough[face_id];
      if (bvar_albedo[face_id] < soilmin[2])
        soilmin[2] = bvar_albedo[face_id];
      if (bvar_emissi[face_id] < soilmin[3])
        soilmin[3] = bvar_emissi[face_id];

      if (bvar_thermal_capacity[soil_id] < soilmin[4])
        soilmin[4] = bvar_thermal_capacity[soil_id];
      if (bvar_vegeta[soil_id] < soilmin[5])
        soilmin[5] = bvar_vegeta[soil_id];
      if (bvar_water_capacity[soil_id] < soilmin[6])
        soilmin[6] = bvar_water_capacity[soil_id];
      if (bvar_water_ratio[soil_id] < soilmin[7])
        soilmin[7] = bvar_water_ratio[soil_id];
      if (bvar_r1[soil_id] < soilmin[8])
        soilmin[8] = bvar_r1[soil_id];
      if (bvar_r2[soil_id] < soilmin[9])
        soilmin[9] = bvar_r2[soil_id];

      soilmean[0] += b_face_surf[face_id] * bvar_dyn_rough[face_id];
      soilmean[1] += b_face_surf[face_id] * bvar_therm_rough[face_id];
      soilmean[2] += b_face_surf[face_id] * bvar_albedo[face_id];
      soilmean[3] += b_face_surf[face_id] * bvar_emissi[face_id];

      soilmean[4] += b_face_surf[face_id] * bvar_thermal_capacity[soil_id];
      soilmean[5] += b_face_surf[face_id] * bvar_vegeta[soil_id];
      soilmean[6] += b_face_surf[face_id] * bvar_water_capacity[soil_id];
      soilmean[7] += b_face_surf[face_id] * bvar_water_ratio[soil_id];
      soilmean[8] += b_face_surf[face_id] * bvar_r1[soil_id];
      soilmean[9] += b_face_surf[face_id] * bvar_r2[soil_id];

      /* Surface of the zone, could use it directly in C */
      surf_zone += b_face_surf[face_id];
    }

    cs_parall_sum(10, CS_DOUBLE, soilmean);
    cs_parall_min(10, CS_DOUBLE, soilmin);
    cs_parall_max(10, CS_DOUBLE, soilmax);
    cs_parall_sum(1, CS_DOUBLE, &surf_zone);

    for (cs_lnum_t n = 0; n < 10; n++) {
      soilmean[n] = soilmean[n] / surf_zone;
    }

    bft_printf(" ** ========================================= **\n"
               " ** Soil/atmosphere interface                 **\n"
               " ** Array of constants                        **\n"
               " ** ========================================= **\n");
    bft_printf(" *            * minimum*    mean* maximum*\n");
    for (cs_lnum_t n = 0; n < 10; n++) {
      if (n == 4) {
        bft_printf(" *%-12s*%8.4f*%8.4f*%8.4f*\n",
                   soilname[n], soilmin[n] * 1e6, soilmean[n] * 1e6,
                   soilmax[n] * 1e6);
      }
      else {
        bft_printf(" *%-12s*%8.4f*%8.4f*%8.4f*\n",
                   soilname[n], soilmin[n], soilmean[n], soilmax[n]);
      }
    }
    bft_printf(" ** ========================================= **\n");
  }
}

/*----------------------------------------------------------------------------*/
