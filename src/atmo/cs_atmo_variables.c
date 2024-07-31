/*============================================================================
 * Functions relative to fields atmospheric
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
        Add fields atmospheric.

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
 * Add atmospheric property fields
 */
/*----------------------------------------------------------------------------*/

/*void
cs_atmo_add_property_fields(void)
{

} */

/*----------------------------------------------------------------------------*/

END_C_DECLS
