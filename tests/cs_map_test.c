/*============================================================================
 * Unit test for cs_map.c;
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#include <stdlib.h>
#include <string.h>

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_map.h"
#include "cs_timer.h"

/*---------------------------------------------------------------------------*/

const char *k_base[] =
  {"boundary_mass_flux_id",
   "boundary_value_id",
   "coupled",
   "diffusivity_tensor",
   "drift_scalar_model",
   "first_moment_id",
   "gradient_weighting_id",
   "inner_mass_flux_id",
   "label",
   "log",
   "max_scalar_clipping",
   "min_scalar_clipping",
   "moment_id",
   "post_id",
   "post_probes",
   "post_vis",
   "property_id",
   "scalar_class",
   "diffusivity_id",
   "diffusivity_ref",
   "scalar_id",
   "slope_test_upwind_id",
   "source_term_prev_id",
   "turbulent_flux_id",
   "turbulent_flux_model",
   "variable_id",
   "solving_info",
   "var_cal_opt"};

const char *f_base[] =
  {"velocity",
   "pressure",
   "k",
   "epsilon",
   "enthalpy",
   "n_p_01",
   "x_p_coal_01",
   "x_p_char_01",
   "x_p_h_01",
   "x_c_h",
   "fr_mv1_01",
   "fr_mv2_01",
   "fr_het_o2",
   "f1f2_variance",
   "x_c_co2",
   "x_c_hcn",
   "x_c_nh3",
   "x_c_no",
   "x_c_h_ox",
   "density",
   "boundary_density",
   "molecular_viscosity",
   "turbulent_viscosity",
   "courant_number",
   "fourier_number",
   "total_pressure",
   "t_gas",
   "rho_gas",
   "ym_chx1m",
   "ym_chx2m",
   "ym_co",
   "ym_h2s",
   "ym_h2",
   "ym_hcn",
   "ym_nh3",
   "ym_o2",
   "ym_co2",
   "ym_h2o",
   "ym_so2",
   "ym_n2",
   "xm",
   "exp1",
   "exp2",
   "exp3",
   "exp4",
   "exp5",
   "f_hcn_dev",
   "f_hcn_het",
   "f_nh3_dev",
   "f_nh3_het",
   "f_no_hcn",
   "f_no_nh3",
   "f_no_het",
   "f_no_the",
   "c_no_hcn",
   "c_no_nh3",
   "f_hcn_rb",
   "c_no_rb",
   "exp_rb",
   "t_p_01",
   "x_p_01",
   "rho_p_01",
   "diam_p_01",
   "dissapear_rate_p_01",
   "m_transfer_v1_p_01",
   "m_transfer_v2_p_01",
   "het_ts_o2_p_01",
   "imp_m_transfer_to_g_p_01",
   "x_c",
   "b_x_c",
   "x_h_c_exp_st",
   "x_h_c_imp_st",
   "x_carbone",
   "x_oxygen",
   "x_hydrogen",
   "luminance",
   "radiative_flux",
   "rad_st",
   "rad_st_implicit",
   "rad_absorption",
   "rad_emission",
   "rad_absorption_coeff",
   "rad_st_02",
   "rad_st_implicit_02",
   "rad_absorption_02",
   "rad_emission_02",
   "rad_absorption_coeff_02",
   "wall_temperature",
   "rad_incident_flux",
   "wall_thermal_conductivity",
   "wall_thickness",
   "emissivity",
   "rad_net_flux",
   "rad_convective_flux",
   "rad_exchange_coefficient",
   "inner_mass_flux",
   "boundary_mass_flux",
   "boundary_ym_fuel",
   "boundary_ym_oxydizer",
   "boundary_ym_product"};

/*---------------------------------------------------------------------------*/

int
main (int argc, char *argv[])
{
  CS_UNUSED(argc);
  CS_UNUSED(argv);

  bft_mem_init(getenv("CS_MEM_LOG"));

  /* Cont elements */

  int n_keys = 0, n_fields = 0;
  size_t tot;

  for (n_keys = 0, tot= 0; tot < sizeof(k_base); n_keys++)
    tot += strlen(k_base[n_keys]);
  for (n_fields = 0, tot = 0; tot < sizeof(f_base); n_fields++)
    tot += strlen(f_base[n_fields]);

  bft_printf("n f %d %d\n", n_keys, (int)sizeof(k_base));
  bft_printf("n f %d %d\n", n_fields, (int)sizeof(f_base));

  cs_map_name_to_id_t *k_map = cs_map_name_to_id_create();
  cs_map_name_to_id_t *f_map = cs_map_name_to_id_create();

  for (int i = 0; i < n_keys; i++)
    cs_map_name_to_id(k_map, k_base[i]);

  for (int i = 0; i < n_fields; i++)
    cs_map_name_to_id(f_map, f_base[i]);

  /* Timings */

  cs_timer_counter_t k_time, f_time;

  CS_TIMER_COUNTER_INIT(k_time);
  CS_TIMER_COUNTER_INIT(f_time);

  int n_passes = 500;
  int l = 0;

  for (int k = 0; k < n_passes; k++) {

    /* Time keys */

    cs_timer_t t0 = cs_timer_time();

    for (int i = 0; i < n_keys; i++)
      l += cs_map_name_to_id(k_map, k_base[i]);

    cs_timer_t t1 = cs_timer_time();

    cs_timer_counter_add_diff(&k_time, &t0, &t1);

    /* Time fields */

    t0 = cs_timer_time();

    for (int i = 0; i < n_fields; i++)
      l += cs_map_name_to_id(f_map, f_base[i]);

    t1 = cs_timer_time();

    cs_timer_counter_add_diff(&f_time, &t0, &t1);

  }

  cs_real_t k_query = (k_time.wall_nsec / n_passes / n_keys) * 1.e-9;
  cs_real_t f_query = (f_time.wall_nsec / n_passes / n_fields) * 1.e-9;

  bft_printf("mean key   id query time: %g\n", k_query);
  bft_printf("mean field id query time: %g\n", f_query);

  cs_map_name_to_id_destroy(&k_map);
  cs_map_name_to_id_destroy(&f_map);

  bft_mem_end();

  exit (EXIT_SUCCESS);
}
