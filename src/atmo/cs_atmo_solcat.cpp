/*============================================================================
 * Atmospheric soil module
 * Soil - atmosphere parameters computed from a "Land use" file
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
#include "mesh/cs_mesh.h"
#include "atmo/cs_air_props.h"
#include "base/cs_math.h"
#include "alge/cs_divergence.h"
#include "bft/bft_mem.h"
#include "bft/bft_error.h"
#include "bft/bft_printf.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "atmo/cs_atmo_solcat.h"

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Soil - atmosphere parameters computed from a "Land use" file
 *
 * \param[in] call_stage  first pass to set default values,
 *                        second pass to perform some checks and log
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_soil_cat(int call_stage)
{
  /* Local variables */
  int water, forest, diverse, mineral, diffu, mixed, dense, building;
  const int inityp = -9;
  const double codinv = -999.0;
  char nomcat[8][10];

  cs_atmo_option_t *at_opt = cs_glob_atmo_option;
  cs_lnum_t n_elts;
  int n_soil_cat;
  const cs_lnum_t *elt_ids;

  cs_f_atmo_get_soil_zone(&n_elts, &n_soil_cat, &elt_ids);

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
  if (n_soil_cat == 7) {
    water = 1;
    forest = 2;
    diverse = 3;
    mineral = 4;
    diffu = 5;
    mixed = 6;
    dense = 7;
  }
  /* Case of a file based on IGN data in 5 categories */
  else if (n_soil_cat == 5) {
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

  /* First pass, default values according to the choice of number of soils */
  /*----------------------------------------------------------------------*/

  if (call_stage == 1) {
    cs_atmo_soil_init_arrays(&n_soil_cat,
                             &at_opt->soil_cat_thermal_inertia,
                             &at_opt->soil_cat_roughness,
                             &at_opt->soil_cat_thermal_roughness,
                             &at_opt->soil_cat_albedo,
                             &at_opt->soil_cat_emissi,
                             &at_opt->soil_cat_vegeta,
                             &at_opt->soil_cat_w1,
                             &at_opt->soil_cat_w2,
                             &at_opt->soil_cat_r1,
                             &at_opt->soil_cat_r2);

    for (cs_lnum_t n = 1; n <= n_soil_cat; n++) {
      at_opt->soil_cat_thermal_inertia[n] = codinv;
      at_opt->soil_cat_roughness[n] = codinv;
      at_opt->soil_cat_thermal_roughness[n] = codinv;
      at_opt->soil_cat_albedo[n] = codinv;
      at_opt->soil_cat_emissi[n] = codinv;
      at_opt->soil_cat_vegeta[n] = codinv;
      at_opt->soil_cat_w1[n] = codinv;
      at_opt->soil_cat_w2[n] = codinv;
      at_opt->soil_cat_r1[n] = codinv;
      at_opt->soil_cat_r2[n] = codinv;
    }

    /* Standard values of parameters */
    if (water != inityp) at_opt->soil_cat_roughness[water] = 0.0005;
    if (forest != inityp) at_opt->soil_cat_roughness[forest] = 0.800;
    if (diverse != inityp) at_opt->soil_cat_roughness[diverse] = 0.100;
    if (mineral != inityp) at_opt->soil_cat_roughness[mineral] = 0.0012;
    if (diffu != inityp) at_opt->soil_cat_roughness[diffu] = 0.250;
    if (mixed != inityp) at_opt->soil_cat_roughness[mixed] = 0.600;
    if (dense != inityp) at_opt->soil_cat_roughness[dense] = 1.000;
    if (building != inityp) at_opt->soil_cat_roughness[building] = 0.600;

    if (water != inityp) {
      at_opt->soil_cat_thermal_roughness[water]
        = at_opt->soil_cat_roughness[water];
    }
    if (forest != inityp) {
      at_opt->soil_cat_thermal_roughness[forest]
        = at_opt->soil_cat_roughness[forest] * exp(-2.0);
    }
    if (diverse != inityp) {
      at_opt->soil_cat_thermal_roughness[diverse]
        = at_opt->soil_cat_roughness[diverse] * exp(-2.0);
    }
    if (mineral != inityp) {
      at_opt->soil_cat_thermal_roughness[mineral]
        = at_opt->soil_cat_roughness[mineral] * exp(-2.0);
    }
    if (diffu != inityp) {
      at_opt->soil_cat_thermal_roughness[diffu]
        = at_opt->soil_cat_roughness[diffu] * exp(-2.0);
    }
    if (mixed != inityp) {
      at_opt->soil_cat_thermal_roughness[mixed]
        = at_opt->soil_cat_roughness[mixed] * exp(-2.0);
    }
    if (dense != inityp) {
      at_opt->soil_cat_thermal_roughness[dense]
        = at_opt->soil_cat_roughness[dense] * exp(-2.0);
    }
    if (building != inityp) {
      at_opt->soil_cat_thermal_roughness[building]
        = at_opt->soil_cat_roughness[building] * exp(-2.0);
    }

    if (water != inityp) at_opt->soil_cat_albedo[water] = 0.08;
    if (forest != inityp) at_opt->soil_cat_albedo[forest] = 0.16;
    if (diverse != inityp) at_opt->soil_cat_albedo[diverse] = 0.20;
    if (mineral != inityp) at_opt->soil_cat_albedo[mineral] = 0.25;
    if (diffu != inityp) at_opt->soil_cat_albedo[diffu] = 0.18;
    if (mixed != inityp) at_opt->soil_cat_albedo[mixed] = 0.18;
    if (dense != inityp) at_opt->soil_cat_albedo[dense] = 0.18;
    if (building != inityp) at_opt->soil_cat_albedo[building] = 0.18;

    if (water != inityp) at_opt->soil_cat_emissi[water] = 0.980;
    if (forest != inityp) at_opt->soil_cat_emissi[forest] = 0.950;
    if (diverse != inityp) at_opt->soil_cat_emissi[diverse] = 0.940;
    if (mineral != inityp) at_opt->soil_cat_emissi[mineral] = 0.965;
    if (diffu != inityp) at_opt->soil_cat_emissi[diffu] = 0.880;
    if (mixed != inityp) at_opt->soil_cat_emissi[mixed] = 0.880;
    if (dense != inityp) at_opt->soil_cat_emissi[dense] = 0.880;
    if (building != inityp) at_opt->soil_cat_emissi[building] = 0.880;

    if (water != inityp) at_opt->soil_cat_vegeta[water] = 0.00;
    if (forest != inityp) at_opt->soil_cat_vegeta[forest] = 1.00;
    if (diverse != inityp) at_opt->soil_cat_vegeta[diverse] = 1.00;
    if (mineral != inityp) at_opt->soil_cat_vegeta[mineral] = 0.00;
    if (diffu != inityp) at_opt->soil_cat_vegeta[diffu] = 0.50;
    if (mixed != inityp) at_opt->soil_cat_vegeta[mixed] = 0.25;
    if (dense != inityp) at_opt->soil_cat_vegeta[dense] = 0.00;
    if (building != inityp) at_opt->soil_cat_vegeta[building] = 0.25;

    if (water != inityp) at_opt->soil_cat_thermal_inertia[water] = 7.6e-06;
    if (forest != inityp) at_opt->soil_cat_thermal_inertia[forest] = 11.0e-06;
    if (diverse != inityp) at_opt->soil_cat_thermal_inertia[diverse] = 11.0e-06;
    if (mineral != inityp) at_opt->soil_cat_thermal_inertia[mineral] = 5.0e-06;
    if (dense != inityp) at_opt->soil_cat_thermal_inertia[dense] = 3.9e-06;
    if (diffu != inityp) {
      at_opt->soil_cat_thermal_inertia[diffu]
        =      at_opt->soil_cat_thermal_inertia[forest]
             * at_opt->soil_cat_vegeta[diffu]
           +   at_opt->soil_cat_thermal_inertia[dense]
             * (1.0 - at_opt->soil_cat_vegeta[diffu]);
    }
    if (mixed != inityp) {
      at_opt->soil_cat_thermal_inertia[mixed]
        =     at_opt->soil_cat_thermal_inertia[forest]
            * at_opt->soil_cat_vegeta[mixed]
          +   at_opt->soil_cat_thermal_inertia[dense]
            * (1.0 - at_opt->soil_cat_vegeta[mixed]);
    }
    if (building != inityp) {
      at_opt->soil_cat_thermal_inertia[building]
        =     at_opt->soil_cat_thermal_inertia[forest]
            * at_opt->soil_cat_vegeta[building]
          +   3.9e-06 * (1.0 - at_opt->soil_cat_vegeta[building]);
    }

    if (water != inityp) at_opt->soil_cat_w1[water] = 100.0;
    if (forest != inityp) {
      at_opt->soil_cat_w1[forest] = 18.0 * at_opt->soil_cat_vegeta[forest] + 2.0;
    }
    if (diverse != inityp) {
      at_opt->soil_cat_w1[diverse]
        = 18.0 * at_opt->soil_cat_vegeta[diverse] + 2.0;
    }
    if (mineral != inityp) {
      at_opt->soil_cat_w1[mineral]
        = 18.0 * at_opt->soil_cat_vegeta[mineral] + 2.0;
    }
    if (diffu != inityp) {
      at_opt->soil_cat_w1[diffu] = 18.0 * at_opt->soil_cat_vegeta[diffu] + 2.0;
    }
    if (mixed != inityp) {
      at_opt->soil_cat_w1[mixed] = 18.0 * at_opt->soil_cat_vegeta[mixed] + 2.0;
    }
    if (dense != inityp) {
      at_opt->soil_cat_w1[dense] = 18.0 * at_opt->soil_cat_vegeta[dense] + 2.0;
    }
    if (building != inityp) {
      at_opt->soil_cat_w1[building]
        = 18.0 * at_opt->soil_cat_vegeta[building] + 2.0;
    }

    if (water != inityp) at_opt->soil_cat_w2[water] = 1.00;
    if (forest != inityp) at_opt->soil_cat_w2[forest] = 0.20;
    if (diverse != inityp) at_opt->soil_cat_w2[diverse] = 0.20;
    if (mineral != inityp) at_opt->soil_cat_w2[mineral] = 0.20;
    if (diffu != inityp) at_opt->soil_cat_w2[diffu] = 0.20;
    if (mixed != inityp) at_opt->soil_cat_w2[mixed] = 0.20;
    if (dense != inityp) at_opt->soil_cat_w2[dense] = 0.20;
    if (building != inityp) at_opt->soil_cat_w2[building] = 0.20;

    if (water != inityp) at_opt->soil_cat_r1[water] = 0.0;
    if (forest != inityp) at_opt->soil_cat_r1[forest] = 0.0;
    if (diverse != inityp) at_opt->soil_cat_r1[diverse] = 0.0;
    if (mineral != inityp) at_opt->soil_cat_r1[mineral] = 0.0;
    if (dense != inityp) at_opt->soil_cat_r1[dense] = 30.0;
    if (diffu != inityp) at_opt->soil_cat_r1[diffu] = 10.0;
    if (mixed != inityp) at_opt->soil_cat_r1[mixed] = 15.0;
    if (building != inityp) at_opt->soil_cat_r1[building] = 15.0;

    if (water != inityp) at_opt->soil_cat_r2[water] = 0.0;
    if (forest != inityp) at_opt->soil_cat_r2[forest] = 0.0;
    if (diverse != inityp) at_opt->soil_cat_r2[diverse] = 0.0;
    if (mineral != inityp) at_opt->soil_cat_r2[mineral] = 0.0;
    if (dense != inityp) at_opt->soil_cat_r2[dense] = 2.0;
    if (diffu != inityp) at_opt->soil_cat_r2[diffu] = 2.0 / 3.0;
    if (mixed != inityp) at_opt->soil_cat_r2[mixed] = 1.0;
    if (building != inityp) at_opt->soil_cat_r2[building] = 1.0;
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

    for (cs_lnum_t n = 1; n <= n_soil_cat; n++) {
      bft_printf
        ("%s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
         nomcat[n], at_opt->soil_cat_roughness[n],
         at_opt->soil_cat_thermal_roughness[n],
        at_opt->soil_cat_albedo[n],
         at_opt->soil_cat_emissi[n],
         1.e+06 * at_opt->soil_cat_thermal_inertia[n],
        at_opt->soil_cat_vegeta[n],
         at_opt->soil_cat_w1[n], at_opt->soil_cat_w2[n],
         at_opt->soil_cat_r1[n], at_opt->soil_cat_r2[n]);
    }
    bft_printf(" --------------------------------\n\n");

    /* Check */
    int n_errors = n_soil_cat;

    for (cs_lnum_t n = 1; n <= n_soil_cat; n++) {
      if (   at_opt->soil_cat_roughness[n] > codinv
          && at_opt->soil_cat_thermal_roughness[n] > codinv
          && at_opt->soil_cat_albedo[n] > codinv
          && at_opt->soil_cat_emissi[n] > codinv
          && at_opt->soil_cat_w1[n] > codinv
          && at_opt->soil_cat_w2[n] > codinv
          && at_opt->soil_cat_thermal_inertia[n] > codinv
          && at_opt->soil_cat_r1[n] > codinv
          && at_opt->soil_cat_r2[n] > codinv) {
        n_errors--;
      }
    }

    /* Error message if needed */
    if (n_errors != 0) {
      bft_error(__FILE__, __LINE__, 0,
                _("%s: number of soils with incorrect soil coefficients = %d."),
                __func__, n_errors);
    }
  }
}

/*----------------------------------------------------------------------------*/
