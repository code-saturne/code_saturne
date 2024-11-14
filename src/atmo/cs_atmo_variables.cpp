/*============================================================================
 * Functions relative to atmospheric model fields.
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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"

#include "cs_atmo.h"
#include "cs_atmo_aerosol.h"
#include "cs_atmo_profile_std.h"
#include "cs_air_props.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_intprf.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_physical_properties.h"
#include "cs_post.h"
#include "cs_prototypes.h"
#include "cs_thermal_model.h"
#include "cs_turbulence_model.h"
#include "cs_velocity_pressure.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_atmo_variables.h"

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_atmo_variables.cpp
        Add atmospheric model fields.
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Type and macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Internal function -
 *        subgrid condensation scheme assuming
 *        a gaussian distribution for the
 *        fluctuations of both qw and thetal.
 */
/*----------------------------------------------------------------------------*/

static void
_gaussian(const cs_mesh_t             *m,
          const cs_mesh_quantities_t  *mq,
          const cs_atmo_option_t      *at_opt,
          const cs_fluid_properties_t *fluid_props,
          cs_real_t                   *crom,
          cs_real_t                   *cpro_tempc,
          cs_real_t                   *cpro_liqwt,
          const cs_real_t             *cpro_met_p,
          const cs_real_t             *cvar_totwt)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_real_3_t *cell_cen = (const cs_real_3_t *)mq->cell_cen;

  const cs_real_6_t *cvar_rij = nullptr;
  const cs_real_t *cvar_k = nullptr, *cvar_ep = nullptr;
  const cs_real_t *cvar_omg = nullptr, *cvar_nusa = nullptr;

  const cs_field_t *th_f = cs_thermal_model_field();
  const cs_turb_model_t *turb_model = cs_glob_turb_model;

  if (turb_model->itytur == 2 || turb_model->model == CS_TURB_V2F_PHI) {
    cvar_k = CS_F_(k)->val;
    cvar_ep = CS_F_(eps)->val;
  }
  else if (turb_model->order == CS_TURB_SECOND_ORDER) {
    cvar_ep = CS_F_(eps)->val;
    cvar_rij = (const cs_real_6_t *)(CS_F_(rij)->val);
  }
  else if (turb_model->model == CS_TURB_K_OMEGA) {
    cvar_k = CS_F_(k)->val;
    cvar_omg = CS_F_(omg)->val;
  }
  else if (turb_model->model == CS_TURB_SPALART_ALLMARAS)
    cvar_nusa = CS_F_(nusa)->val;

  cs_real_t *nn = cs_field_by_name("nebulosity_frac")->val;
  cs_real_t *nebdia = cs_field_by_name("nebulosity_diag")->val;

  /* Gradients are used for estimating standard
     deviations of the subgrid fluctuations */

  cs_real_t a_coeff = 0.0;
  const cs_real_t cp0 = fluid_props->cp0;
  const cs_real_t rvsra = fluid_props->rvsra;
  const cs_real_t clatev = fluid_props->clatev;
  const cs_real_t rair = fluid_props->r_pg_cnst;
  const cs_real_t ps = cs_glob_atmo_constants->ps;
  const cs_real_t a_const = 2.0*cs_turb_cmu/2.3;
  const cs_real_t rscp = fluid_props->r_pg_cnst/fluid_props->cp0;
  const cs_real_t rvap = fluid_props->r_pg_cnst*fluid_props->rvsra;
  const cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

  const cs_real_t *cvar_vart = th_f->val;

  cs_real_3_t *dqsd = nullptr, *dtlsd = nullptr;
  CS_MALLOC_HD(dqsd, n_cells_ext, cs_real_3_t, cs_alloc_mode);
  CS_MALLOC_HD(dtlsd, n_cells_ext, cs_real_3_t, cs_alloc_mode);

  cs_field_gradient_scalar(th_f, true, 1, dtlsd);
  cs_field_gradient_scalar(cs_field_by_name("ym_water"), true, 1, dqsd);

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    // coeff = 2 cmu/c2 * k^3 / eps^2
    if (turb_model->itytur == 2 || turb_model->model == CS_TURB_V2F_PHI)
      a_coeff = a_const*cs_math_pow3(cvar_k[c_id])
              / cs_math_pow2(cvar_ep[c_id]);
    else if (turb_model->order == CS_TURB_SECOND_ORDER) {
      const cs_real_t ek = 0.5*cs_math_6_trace(cvar_rij[c_id]);
      a_coeff = a_const*cs_math_pow3(ek)/cs_math_pow2(cvar_ep[c_id]);
    }
    else if (turb_model->model == CS_TURB_K_OMEGA) {
      const cs_real_t ep = cvar_omg[c_id]*cvar_k[c_id]*cs_turb_cmu;
      a_coeff = a_const*cs_math_pow3(cvar_k[c_id])/cs_math_pow2(ep);
    }
    else if (cvar_nusa != nullptr)
      /* using cvar_nusa[c_id] = cmu*xkent^2/xeent
       * FIXME: There is no good way to calculate tke and eps from nusa.
       * For the moment we use tke^4/eps^2 instead of tk^3/eps^2
       * Need to return WARNING that in case of Spalart-Allmaras we use bad assumpltion
       * or RETURN error for this case. */
      a_coeff = a_const*cs_math_pow2(cvar_nusa[c_id])/cs_math_pow2(cs_turb_cmu);

    cs_real_t pp = 0.0, dum = 0.0;
    const cs_real_t zent = cell_cen[c_id][2];
    if (at_opt->meteo_profile == 0)
      cs_atmo_profile_std(zent, &pp, &dum, &dum);
    else if (at_opt->meteo_profile == 1)
      pp = cs_intprf(at_opt->met_1d_nlevels_t,
                     at_opt->met_1d_ntimes,
                     at_opt->z_temp_met,
                     at_opt->time_met,
                     at_opt->hyd_p_met,
                     zent,
                     cs_glob_time_step->t_cur);
    else
      pp = cpro_met_p[c_id];

     const cs_real_t xvart = cvar_vart[c_id]; // thermal scalar: liquid potential temperature
     const cs_real_t tliq = xvart*pow(pp/ps, rscp); // liquid temperature
     const cs_real_t qsl = cs_air_yw_sat(tliq-tkelvi, pp); // saturated vapor content
     const cs_real_t alpha = (clatev*qsl/(rvap*pow(tliq,2)))*pow(pp/ps, rscp);
     const cs_real_t var_q_tl = a_coeff *
                                (  pow(dqsd[c_id][0] - alpha * dtlsd[c_id][0], 2)
                                 + pow(dqsd[c_id][1] - alpha * dtlsd[c_id][1], 2)
                                 + pow(dqsd[c_id][2] - alpha * dtlsd[c_id][2], 2)  );

     const cs_real_t sig_flu = fmax(sqrt(var_q_tl), 1e-30);
     const cs_real_t qwt = cvar_totwt[c_id]; // total water content
     cs_real_t deltaq = qwt - qsl;
     const cs_real_t q1 = deltaq/sig_flu;

     nebdia[c_id] = 0.5*(1.0 + erf(q1/sqrt(2.0)));

     // FIXME MF : put in input of the global function...
     cs_real_t yw_liq = (sig_flu /(1.0 + qsl*pow(clatev, 2)/(rvap*cp0*pow(tliq, 2))))
                      * (nebdia[c_id]*q1 + exp(-pow(q1, 2)/2.0)/sqrt(2.0*cs_math_pi));
     yw_liq   = fmax(yw_liq, 0.0);
     nn[c_id] = nebdia[c_id] - (  nebdia[c_id]*q1
                                + exp(-pow(q1, 2)/2.0)/sqrt(2.0*cs_math_pi))
                                * exp(-pow(q1, 2)/2.0)/sqrt(2.0*cs_math_pi);

     // go back to all or nothing
     if (qwt < yw_liq) {
       nn[c_id] = 0.0;
       // deltaq set to 0 if unsaturted air parcel
       if (deltaq < 0.0) {
         deltaq = 0.0;
         nebdia[c_id] = 0.0;
       }
       else {
         nebdia[c_id] = 1.0;
       }

       /* TODO input ?
        * 0 if unsaturated air parcel */
       yw_liq = deltaq / (1.0 + qsl*pow(clatev, 2)/(rvap*cp0*pow(tliq, 2)));
     }

     //Celcius temperature of the air parcel
     cpro_tempc[c_id] = tliq + (clatev/cp0)*yw_liq - tkelvi;
     // liquid water content
     cpro_liqwt[c_id] = yw_liq;
     //density
     const cs_real_t lrhum = rair*(1.0 - yw_liq + (rvsra - 1.0)*(qwt - yw_liq));
     crom[c_id] = pp/(lrhum*(tliq + (clatev/cp0)*yw_liq));

  } // end loop on cells

  CS_FREE_HD(dqsd);
  CS_FREE_HD(dtlsd);

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * Add atmospheric variables fields
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_add_variable_fields(void)
{
  /* Key ids for clipping */
  const int kscmin = cs_field_key_id("min_scalar_clipping");
  const int kscmax = cs_field_key_id("max_scalar_clipping");

  cs_atmo_option_t *at_opt = cs_glob_atmo_option;
  bool rain = at_opt->rain;

  if (rain == true) {
    cs_field_t *f;

    /*Create the mass fraction of rain */
    int f_id = cs_variable_field_create("ym_l_r",
                                        "Mass frac rain",
                                        CS_MESH_LOCATION_CELLS,
                                        1);
    f = cs_field_by_id(f_id);

    /* Clipping of rain mass fraction 0 <= ym_l_r <=1 */
    cs_field_set_key_double(f, kscmin, 0.e0);
    cs_field_set_key_double(f, kscmax, 1.e0);

    /* Add the rain mass fraction to the index of fields */
    cs_add_model_field_indexes(f->id);

    /* Equation parameters */
    cs_equation_param_t *eqp = cs_field_get_equation_param(f);
    /* Set beta limiter to maintain y_p in the limits */
    eqp->isstpc = 2;
    /* Full upwind scheme */
    eqp->blencv = 0.0;


    /* Create the concentration of rain droplets */
    f_id = cs_variable_field_create("n_r",
                                    "Number of rain drops",
                                    CS_MESH_LOCATION_CELLS,
                                    1);
    f = cs_field_by_id(f_id);
    eqp = cs_field_get_equation_param(f);

    /* Clipping of rain mass fraction 0 <= ym_l_r <=1 */
    cs_field_set_key_double(f, kscmin, 0.e0);
    cs_field_set_key_double(f, kscmax, 1.e10);

    /* Add the rain mass fraction to the index of fields */
    cs_add_model_field_indexes(f->id);

    /* Set beta limiter to maintain y_p in the limits */
    eqp->isstpc = 2;
    /* Full upwind scheme */
    eqp->blencv = 0.0;
  }

  cs_thermal_model_t *thermal = cs_get_glob_thermal_model();

  /* Add variables
     -------------*/

  // Dry atmosphere
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_DRY) {
    // Potential temperature, in Kelvin
    thermal->temperature_scale = CS_TEMPERATURE_SCALE_KELVIN;
    thermal->thermal_variable = CS_THERMAL_MODEL_TEMPERATURE;

    int f_id = cs_variable_field_create("temperature",
                                        "PotTemp",
                                        CS_MESH_LOCATION_CELLS,
                                        1);
    cs_field_t *f = cs_field_by_id(f_id);
    cs_add_model_thermal_field_indexes(f->id);
    cs_field_set_key_double(f, kscmin, 0.0);
  }

  // Humid atmosphere
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_HUMID) {
    // Potential temperature, in Kelvin
    thermal->temperature_scale = CS_TEMPERATURE_SCALE_KELVIN;
    thermal->thermal_variable = CS_THERMAL_MODEL_TEMPERATURE;

    int f_id = cs_variable_field_create("temperature",
                                        "LqPotTmp",
                                        CS_MESH_LOCATION_CELLS,
                                        1);
    cs_field_t *f = cs_field_by_id(f_id);
    cs_add_model_thermal_field_indexes(f->id);
    cs_field_set_key_double(f, kscmin, 200.0);

    // total water content
    f_id = cs_variable_field_create("ym_water",
                                    "Ym water",
                                    CS_MESH_LOCATION_CELLS,
                                    1);
    f = cs_field_by_id(f_id);
    cs_add_model_field_indexes(f->id);
    cs_field_set_key_double(f, kscmin, 0.0);

    // total number droplets
    f_id = cs_variable_field_create("number_of_droplets",
                                    "TotDrop",
                                    CS_MESH_LOCATION_CELLS,
                                    1);
    f = cs_field_by_id(f_id);
    cs_add_model_field_indexes(f->id);
    cs_field_set_key_double(f, kscmin, 0.0);

  }

  /* Chemistry variables
     ------------------- */
  cs_atmo_chemistry_t *at_chem = cs_glob_atmo_chemistry;

  /* Atmospheric gaseous chemistry
   * Do not change this order */
  if (at_chem->aerosol_model == CS_ATMO_AEROSOL_SSH)
    at_chem->model = 4;
  if (at_chem->frozen_gas_chem && at_chem->aerosol_model == CS_ATMO_AEROSOL_OFF)
    at_chem->model = 0;

  if (at_chem->model > 0) {
    // Set the name of the chemical profiles file
    cs_atmo_set_chem_conc_file_name("chemistry");

    /* Initialization of the chemical scheme */
    cs_atmo_init_chemistry();
  }

  // Atmospheric aerosol chemistryc
  if (at_chem->aerosol_model != CS_ATMO_AEROSOL_OFF) {
    // Set the name of the chemical profiles file
    cs_atmo_set_aero_conc_file_name("aerosols");

    // Verification
    if (at_chem->model != 4)
      bft_error(__FILE__, __LINE__, 0,
                "    WARNING:   STOP WHILE READING INPUT DATA\n"
                "    =========\n"
                "  When aerosol chemistry model is used\n"
                "   a full gaseous scheme (CB05) is automatically used\n"
                " The user cannot specify any other scheme (ichemistry)\n"
                "  Computation CAN NOT run.\n\n"
                "  Check the input data given through the User Interfac\n"
                "  or cs_user_parameters.c\n");


    /* Load shared library
     * Initialise external aerosol code
     * Create variables */
    cs_atmo_aerosol_initialize();

  }

  // Set clippings for gas aerosol species
  if (at_chem->model > 0) {
    for (int ii = 0; ii  < at_chem->n_species; ii++) {
      cs_field_t *f = cs_field_by_id(at_chem->species_to_field_id[ii]);
      cs_field_set_key_double(f, kscmin, 0.0);
    }
  }

  if (at_chem->aerosol_model != CS_ATMO_AEROSOL_OFF) {
    const int n_end
      = at_chem->n_species + at_chem->n_size*(at_chem->n_layer + 1);
    for (int ii = at_chem->n_species; ii < n_end; ii++) {
      cs_field_t *f = cs_field_by_id(at_chem->species_to_field_id[ii]);
      cs_field_set_key_double(f, kscmin, 0.0);
    }
    // Allow large aerosol numbers
    const int n_start = at_chem->n_species + at_chem->n_size*at_chem->n_layer;
    for (int ii = n_start; ii < n_end; ii++) {
      cs_field_t *f = cs_field_by_id(at_chem->species_to_field_id[ii]);
      cs_field_set_key_double(f, kscmax, 1.0e40);
    }
  }

  /* Map to fields and GUI
     --------------------- */
  cs_field_pointer_map_atmospheric(at_chem->n_species,
                                   at_chem->species_to_field_id);

  /* General field and physical properties
     ------------------------------------- */
  cs_get_glob_fluid_properties()->icp = -1;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add if needed the variables fields
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_add_property_fields(void)
{

  cs_atmo_option_t *at_opt = cs_glob_atmo_option;
  bool rain = at_opt->rain;

  if (rain == true) {
   {
    cs_field_t *f;
    int field_type = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY;
    bool has_previous = false;
    const int klbl   = cs_field_key_id("label");
    const int keyvis = cs_field_key_id("post_vis");
    const int keylog = cs_field_key_id("log");
    const int post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;

    /* Continuous phase density (humid air) */
    f = cs_field_create("rho_humid_air",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Density humid air");
    }

  }
  const int klbl   = cs_field_key_id("label");
  const int keyvis = cs_field_key_id("post_vis");
  const int keylog = cs_field_key_id("log");

  int field_type = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY;
  const int post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;

  cs_field_t *f = nullptr;

  /* Momentum source terms */
  if (cs_glob_atmo_option->open_bcs_treatment > 0) {
    f = cs_field_create("momentum_source_terms",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        3,
                        false);
    cs_field_set_key_int(f, keyvis, 1);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "MomentumSourceTerms");
  }

  /* Boundary roughness */
  f = cs_field_create("boundary_roughness",
                      field_type,
                      CS_MESH_LOCATION_BOUNDARY_FACES,
                      1,
                      false);
  cs_field_set_key_int(f, keyvis, 0);
  cs_field_set_key_int(f, keylog, 1);
  cs_field_set_key_str(f, klbl, "Boundary Roughness");

  f = cs_field_create("boundary_thermal_roughness",
                      field_type,
                      CS_MESH_LOCATION_BOUNDARY_FACES,
                      1,
                      false);
  cs_field_set_key_int(f, keyvis, 0);
  cs_field_set_key_int(f, keylog, 1);
  cs_field_set_key_str(f, klbl, "Boundary Thermal Roughness");

  /* Temperature for DRY or HUMID */
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] > CS_ATMO_CONSTANT_DENSITY) {

    f = cs_field_by_name_try("real_temperature");
    cs_physical_property_define_from_field("real_temperature",
                                           field_type,
                                           CS_MESH_LOCATION_CELLS,
                                           1,
                                           false);
    cs_field_set_key_int(f, keyvis, 1);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "RealTemp");

    f = cs_field_create("non_neutral_scalar_correction",
                        field_type,
                        CS_MESH_LOCATION_BOUNDARY_FACES,
                        1,
                        false);
    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 0);
    cs_field_set_key_str(f, klbl, "Non Neutral Scalar Correction");

    f = cs_field_by_name_try("thermal_expansion");
    cs_physical_property_define_from_field("thermal_expansion",
                                           field_type,
                                           CS_MESH_LOCATION_CELLS,
                                           1,
                                           false);
    cs_field_set_key_int(f, keyvis, 1);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Beta");
  }

  /* Liquid water content HUMID */
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_HUMID) {

    f = cs_field_by_name_try("liquid_water");
    cs_physical_property_define_from_field("liquid_water",
                                           field_type,
                                           CS_MESH_LOCATION_CELLS,
                                           1,
                                           false);
    cs_field_set_key_int(f, keyvis, 1);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "LiqWater");

    if (at_opt->sedimentation_model > 0 || at_opt->deposition_model > 0) {
      cs_field_find_or_create("boundary_ustar",
                              field_type,
                              CS_MESH_LOCATION_BOUNDARY_FACES,
                              1,
                              false);
    }

    // Radiative cooling
    if (at_opt->sedimentation_model == 1 || at_opt->radiative_model_1d > 0) {
      cs_physical_property_define_from_field("radiative_cooling",
                                             field_type,
                                             CS_MESH_LOCATION_CELLS,
                                             1,
                                             false);
      f = cs_field_by_name_try("radiative_cooling");
      cs_field_set_key_int(f, keyvis, 1);
      cs_field_set_key_int(f, keylog, 1);
      cs_field_set_key_str(f, klbl, "Radiative cooling");
    }

    // fractional nebulosity
    cs_physical_property_define_from_field("nebulosity_frac",
                                           field_type,
                                           CS_MESH_LOCATION_CELLS,
                                           1,
                                           false);
    f = cs_field_by_name_try("nebulosity_frac");
    cs_field_set_key_int(f, keyvis, 1);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Nebulo frac");

    // Diagnosed nebulosity
    cs_physical_property_define_from_field("nebulosity_diag",
                                           field_type,
                                           CS_MESH_LOCATION_CELLS,
                                           1,
                                           false);
    f = cs_field_by_name_try("nebulosity_diag");
    cs_field_set_key_int(f, keyvis, 1);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Nebulo diag");

    cs_physical_property_define_from_field("droplet_eq_radius",
                                           field_type,
                                           CS_MESH_LOCATION_CELLS,
                                           1,
                                           false);
    f = cs_field_by_name_try("droplet_eq_radius");
    cs_field_set_key_int(f, keyvis, 1);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Drop eq radius3");
  }

  int z_id = cs_glob_atmo_option->soil_zone_id;

  if (z_id > -1) {
    const cs_zone_t *z = cs_boundary_zone_by_id(z_id);
    int soil_num = 5;
    switch (cs_glob_atmo_option->soil_cat) {
    case CS_ATMO_SOIL_5_CAT:
      soil_num = 5;
      break;
    case CS_ATMO_SOIL_7_CAT:
      soil_num = 7;
      break;
    case CS_ATMO_SOIL_23_CAT:
      soil_num = 23;
      break;
    }

    /* Note the size is soil_num+1, first value is undef.
     * */
    f = cs_field_create("atmo_soil_percentages",
                        field_type,
                        z->location_id,
                        soil_num + 1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil percentages");

    /* Boundary variable fields for the soil */
    /*---------------------------------------*/

    /* Soil surface temperature in Celcius */
    f = cs_field_create("soil_temperature",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil T");

    /* Soil surface potential temperature (K) */
    f = cs_field_create("soil_pot_temperature",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        true); /* has_previous */
    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil theta");

    /* Soil total water content */
    f = cs_field_create("soil_total_water",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        true); /* has_previous */
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil qw");

    /* ratio of the shallow reservoir water content to its maximum capacity */
    f = cs_field_create("soil_w1",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil w1");

    /* ratio of the deep reservoir water content to its maximum capacity */
    f = cs_field_create("soil_w2",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil w2");

    /* Incident solar radiative flux */
    f = cs_field_create("soil_solar_incident_flux",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil solar incid flux");

    /* Incident solar radiative flux */
    f = cs_field_create("soil_infrared_incident_flux",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil IR incid flux");

    /* Boundary parameters fields characterizing soil */
    /*------------------------------------------------*/

    f = cs_field_create("boundary_albedo",
                        field_type,
                        /* Note: as for boundary_roughness,
                         * location can be reduced in the future */
                        CS_MESH_LOCATION_BOUNDARY_FACES,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Albedo");

    f = cs_field_by_name_try("emissivity");
    if (f == nullptr)
      f = cs_field_create("emissivity",
                          field_type,
                          /* Note: as for boundary_roughness,
                           * location can be reduced in the future */
                          CS_MESH_LOCATION_BOUNDARY_FACES,
                          1, /* dim */
                          false); /* has_previous */

    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Emissivity");

    f = cs_field_create("boundary_vegetation",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Vegetation");

    /* maximum water capacity of shallow reservoir */
    f = cs_field_create("soil_water_capacity",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil c1w");

    /* ratio of the maximum water capacity of the shallow reservoir to the deep
     * reservoir [0,1] */
    f = cs_field_create("soil_water_ratio",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil c2w");

    /* Thermal inertia of the soil */
    f = cs_field_create("soil_thermal_capacity",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil Cp");

    f = cs_field_create("soil_r1",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil r1");

    f = cs_field_create("soil_r2",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil r2");

    /* Deep soil temperature (in Celsius)
     * FIXME potential value? */
    f = cs_field_create("soil_temperature_deep",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil deep T");

    /* Fields usefull for heat budget plot on the soil boundary */
    f = cs_field_create("soil_sensible_heat",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil sensible heat");

    f = cs_field_create("soil_latent_heat",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil latent heat");

    f = cs_field_create("soil_thermal_rad_upward",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil thermal radiation upward");

    f = cs_field_create("soil_thermal_rad_downward",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil thermal radiation downward");

    f = cs_field_create("soil_visible_rad_absorbed",
                        field_type,
                        z->location_id,
                        1, /* dim */
                        false); /* has_previous */
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Soil visible radiation absorbed");

    if (cs_glob_atmo_option->soil_meb_model > CS_ATMO_SOIL_GENUINE) {
      f = cs_field_create("cover_geometry_ratio",
                          field_type,
                          z->location_id,
                          1,
                          false);
      cs_field_set_key_int(f, keyvis, post_flag);
      cs_field_set_key_int(f, keylog, 1);
      cs_field_set_key_str(f, klbl, "Cover geometry ratio");

      f = cs_field_create("cover_reflectivity",
                          field_type,
                          z->location_id,
                          1,
                          false);
      cs_field_set_key_int(f, keyvis, post_flag);
      cs_field_set_key_int(f, keylog, 1);
      cs_field_set_key_str(f, klbl, "Cover reflectivity");

      f = cs_field_create("cover_temperature_radiative",
                          field_type,
                          z->location_id,
                          1,
                          false);
      cs_field_set_key_int(f, keyvis, post_flag);
      cs_field_set_key_int(f, keylog, 1);
      cs_field_set_key_str(f, klbl, "Cover temperature radiation");
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the thermo physical properties fields for the humid air and
 *        the liquid \n
 *        Remarques :
 *        This function  is called at the beginning of each time step
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_physical_properties_update(void)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_fluid_properties_t *fluid_props = cs_glob_fluid_properties;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_real_3_t *cell_cen = (const cs_real_3_t *)mq->cell_cen;

  const cs_atmo_option_t *at_opt = cs_glob_atmo_option;

  cs_real_t *cpro_beta = nullptr;
  cs_real_t *cpro_met_p = nullptr;
  cs_real_t *cpro_met_rho = nullptr;

  if (cs_field_by_name_try("thermal_expansion") != nullptr)
    cpro_beta = cs_field_by_name_try("thermal_expansion")->val;

  if (at_opt->meteo_profile > 1) {
    cpro_met_p = cs_field_by_name_try("meteo_pressure")->val;
    cpro_met_rho = cs_field_by_name_try("meteo_density")->val;
  }

  /* This routine computes the density and the thermodynamic temperature.
     The computations may require the pressure profile which is here taken from
     the meteo file. If no meteo file is used, the user can
     give the laws for RHO and T in cs_user_physical_properties */

  if (cs_thermal_model_field() == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "@                                                            \n"
              "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
              "@                                                            \n"
              "@ @@ WARNING : STOP WHEN CALCULATING PHYSICAL QUANTITIES     \n"
              "@    =========                                               \n"
              "@    The thermal field is not defined check its definition in"
              "@    the GUI, cs_user_parameters or cs_user_physical_properties"
              "@                                                            \n"
              "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
              "@                                                            \n");

  cs_real_t *cvar_totwt = nullptr, *cpro_liqwt = nullptr;
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_HUMID) {
    cvar_totwt = cs_field_by_name("ym_water")->val;
    cpro_liqwt = cs_field_by_name("liquid_water")->val;
  }
  cs_real_t *crom = CS_F_(rho)->val;
  const cs_real_t *cvar_vart = cs_thermal_model_field()->val;
  cs_real_t *cpro_tempc = cs_field_by_name("real_temperature")->val;

  /* From potential temperature, compute:
   * - Temperature in Celsius
   * - Density
   * -------------------------------------
   * Computes the perfect gas constants according to the physics */

  const cs_real_t rscp = fluid_props->r_pg_cnst/fluid_props->cp0;
  const cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;
  // Adiabatic (constant) potential temperature
  const cs_real_t theta0
    = fluid_props->t0 * pow(fluid_props->p0/cs_glob_atmo_constants->ps, rscp);

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    cs_real_t pp = 0.0, dum = 0.0, qwt = 0.0;
    const cs_real_t zent = cell_cen[c_id][2];

    if (at_opt->meteo_profile == 0)
      cs_atmo_profile_std(zent, &pp, &dum, &dum);
    else if (at_opt->meteo_profile == 1)
      pp = cs_intprf(at_opt->met_1d_nlevels_t,
                     at_opt->met_1d_ntimes,
                     at_opt->z_temp_met,
                     at_opt->time_met,
                     at_opt->hyd_p_met,
                     zent,
                     cs_glob_time_step->t_cur);
    else
      pp = cpro_met_p[c_id];

    // Potential temperature
    // or liquid potential temperature for humid atmosphere
    const cs_real_t xvart = cvar_vart[c_id];

    /* (liquid) temperature
     * law: T = theta * (p/ps) ** (Rair/Cp0) */
    const cs_real_t tliq = xvart*pow(pp/cs_glob_atmo_constants->ps, rscp);

    if (cvar_totwt != nullptr)
      qwt = cvar_totwt[c_id];

    /*  Density in cell centers:
     * ------------------------
     * law: rho = P / ( R_mixture * T_mixture(K) )
     * Boussinesq / anelastic approximation */
    if (cs_glob_velocity_pressure_model->idilat == 0) {
      // Compute T in Celisus
      cpro_tempc[c_id] = tliq - tkelvi;

      //  Boussinesq with respect to the adiabatic density
      //  (so called anelastic approximation)
      if (cpro_met_rho != nullptr)
        crom[c_id] = cpro_met_rho[c_id];
      else
        crom[c_id] = fluid_props->ro0;
      /* "delta rho = - beta0 rho_a delta theta" gives
         "beta0 = 1 / theta0" */
      cpro_beta[c_id] = 1.0 / theta0;
    }
    else {
      cs_real_t beta = 0.0, yw_liq = 0.0;
      cs_rho_humidair(qwt,
                      xvart,
                      pp,
                      &yw_liq,
                      &cpro_tempc[c_id],
                      &crom[c_id],
                      &beta);
      // Humid atmosphere
      if (cpro_liqwt != nullptr)
        cpro_liqwt[c_id] = yw_liq;
      /* Thermal expansion for turbulent production
       * "delta rho = - beta rho delta theta" gives
       *"beta = 1 / theta_v", theta_v the virtual temperature */
      if (cpro_beta != nullptr)
        cpro_beta[c_id] = beta;
    }

  } // end loop

  if (   at_opt->distribution_model == 2
      && cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_HUMID  )
    _gaussian(m,
              mq,
              at_opt,
              fluid_props,
              crom,
              cpro_tempc,
              cpro_liqwt,
              cpro_met_p,
              cvar_totwt);

  /* Update the thermo physical properties
     fields for the humid air and he liquid.
     -----------------------------------------*/

  if (!at_opt->rain)
    return;

  cs_real_t *ym_w = (cs_real_t *)CS_F_(ym_w)->val;     // Water mass fraction
  cs_real_t *yr = cs_field_by_name_try("ym_l_r")->val;   // Rain mass fraction
  cs_real_t *rho_h = cs_field_by_name("rho_humid_air")->val; // Humid air density
  cs_real_t *theta_liq = cs_field_by_name("temperature")->val; // Liq. pot. temp.

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_rho_humidair(ym_w[c_id],
                    theta_liq[c_id],
                    cpro_met_p[c_id],
                    &(cpro_liqwt[c_id]),
                    &(cpro_tempc[c_id]),
                    &(rho_h[c_id]),
                    &(cpro_beta[c_id]));
    /* Homogeneous mixture density */
    crom[c_id] = 1.0 / ((1.0 - yr[c_id])/rho_h[c_id] + yr[c_id]/1000);
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
