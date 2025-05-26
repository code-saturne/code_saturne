/*============================================================================
 * Atmospheric soil module - soil variables initialization
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
#include "pprt/cs_physical_model.h"

/*----------------------------------------------------------------------------
*  Header for the current file
*----------------------------------------------------------------------------*/

#include "atmo/cs_atmo_soliva.h"

/*----------------------------------------------------------------------------*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize soil model variables
 */
/*----------------------------------------------------------------------------*/

void cs_soliva(void)
{
    /* Local variables */
    int iphysi;
    cs_real_t rscp;
    cs_real_t esaini, qsaini, huini, psini;

    cs_lnum_t n_elts;
    int n_soil_cat;
    const cs_lnum_t *elt_ids;
    cs_f_atmo_get_soil_zone(&n_elts, &n_soil_cat, &elt_ids);

    /* Constants */
    const cs_real_t cpvcpa_value = cs_glob_air_props->cp_v / cs_glob_air_props->cp_a;
    const cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;
    const cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();
    cs_atmo_option_t *at_opt = cs_glob_atmo_option;

    /* Get pointers to field values */
    cs_field_t *soil_temperature = cs_field_by_name("soil_temperature");
    cs_field_t *soil_pot_temperature = cs_field_by_name("soil_pot_temperature");
    cs_field_t *soil_total_water = cs_field_by_name("soil_total_water");
    cs_field_t *soil_w1 = cs_field_by_name("soil_w1");
    cs_field_t *soil_w2 = cs_field_by_name("soil_w2");

    cs_real_t *bvar_soil_temp = soil_temperature->val;
    cs_real_t *bvar_potentiel_temp = soil_pot_temperature->val;
    cs_real_t *bvar_total_water = soil_total_water->val;
    cs_real_t *bvar_w1 = soil_w1->val;
    cs_real_t *bvar_w2 = soil_w2->val;

    /* Initialization of t and qv at z0 */
    if (at_opt->soil_humidity > 1.0) {
        /* If qvsini > 1, qvsini represents relative humidity in % */
        esaini = 610.78 * exp(17.2694 * at_opt->soil_surf_temp
                              / (at_opt->soil_surf_temp + tkelvi - 35.86));
        qsaini = esaini / (phys_pro->rvsra * phys_pro->p0
                            + esaini * (1.0 - phys_pro->rvsra));
        at_opt->soil_humidity = at_opt->soil_humidity * qsaini / 100.0;
    }

    /* Loop over soil models */
    for (cs_lnum_t soil_id = 0; soil_id < n_elts; soil_id++) {
        psini = phys_pro->p0;
        if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_HUMID) {
            rscp = (phys_pro->r_pg_cnst / phys_pro->cp0)
                   * (1.0 + (phys_pro->rvsra - cpvcpa_value) * at_opt->soil_humidity);
        } else {
            rscp = (phys_pro->r_pg_cnst / phys_pro->cp0);
        }

        bvar_soil_temp[soil_id] = at_opt->soil_surf_temp;
        bvar_potentiel_temp[soil_id] = (at_opt->soil_surf_temp + tkelvi)
                            * pow((cs_glob_atmo_constants->ps / psini), rscp);
        
        /*TODO A bit too restrictive; it's possible to have water in the soil without a humid atmosphere*/
        if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_HUMID) {
            bvar_total_water[soil_id] = at_opt->soil_humidity;
            bvar_w1[soil_id] = at_opt->soil_w1_ini;
            bvar_w2[soil_id] = at_opt->soil_w2_ini;

            /* If not initialized, compute an approximation of the initial water content of the top reservoir */
            if (bvar_w1[soil_id] < 1.0e-20) {
                esaini = 610.78 * exp(17.2694 * at_opt->soil_surf_temp
                                      / (at_opt->soil_surf_temp + tkelvi - 35.86));
                qsaini = esaini / (phys_pro->rvsra * psini
                                    + esaini * (1.0 - phys_pro->rvsra));
                huini = at_opt->soil_humidity / qsaini;
                huini = cs::min(huini, 1.0);
                bvar_w1[soil_id] = acos(1.0 - 2.0 * huini) / acos(-1.0);
            }
            /* If the deep reservoir is uninitialized, set the ratio w1/w2 to 1 (layer equilibrium) */
            if (bvar_w2[soil_id] < 1.0e-20) {
                bvar_w2[soil_id] = bvar_w1[soil_id];
            }
        } else {
            bvar_total_water[soil_id] = 0.0;
            bvar_w1[soil_id] = 0.0;
            bvar_w2[soil_id] = 0.0;
        }
    }
}
