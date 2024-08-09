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
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_physical_properties.h"
#include "cs_post.h"
#include "cs_prototypes.h"
#include "cs_thermal_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_atmo_variables.h"

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_atmo_variables.c
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

  cs_field_t *f = NULL;

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
    if (f == NULL)
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

END_C_DECLS
