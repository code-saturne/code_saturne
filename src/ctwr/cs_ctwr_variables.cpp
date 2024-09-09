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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

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
#include "cs_thermal_model.h"
#include "cs_velocity_pressure.h"
#include "cs_volume_zone.h"

#include "cs_ctwr.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ctwr_variables.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS


/*----------------------------------------------------------------------------
 * Add variables fields
 *----------------------------------------------------------------------------*/

void
cs_ctwr_add_variable_fields(void)
{
  /* Key id of the scalar class */
  const int keyccl = cs_field_key_id("scalar_class");

  /* Key id for drift scalar */
  const int keydri = cs_field_key_id("drift_scalar_model");

  /* Key ids for clipping */
  const int kscmin = cs_field_key_id("min_scalar_clipping");
  const int kscmax = cs_field_key_id("max_scalar_clipping");

  /* Key id for the diffusivity */
  const int kivisl = cs_field_key_id("diffusivity_id");

  /* Fluid properties and physical variables */
  cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();

  /* Set fluid properties parameters */

  /* Variable density */
  fp->irovar = 1;
  /* Activate compressibility */
  cs_velocity_pressure_model_t *vp_model =
    cs_get_glob_velocity_pressure_model();
  vp_model->idilat = 2;
  /* Constant molecular viscosity */
  fp->ivivar = 0;

  /* 1. Definition of fields
   * --------------------------------------------------------------------------
   *  Bulk definition - For cooling towers, the bulk is the humid air.
   *  By definition, humid air is composed of two species: dry air and water
   *  vapor (whether in gas or condensate form)
   *  -------------------------------------------------------------------------
   */

  cs_field_t *f;

  {
    /* CONTINUOUS PHASE - HUMID AIR VARIABLES
     * ====================================== */

    int class_id = -1;
    /* NB : 'c' stands for continuous and 'p' for particles */

    /* Total water mass fraction in bulk humid air */
    {
      f = cs_field_by_name_try("ym_water");
      /* If not using the atmospheric module, we create the field */
      if (f == NULL) {
        int f_id = cs_variable_field_create("ym_water",
            "Mass frac water air",
            CS_MESH_LOCATION_CELLS,
            1);
        f = cs_field_by_id(f_id);
        cs_add_model_field_indexes(f->id);
      }
      /* Clipping : 0 < ym < 1 */
      cs_field_set_key_double(f, kscmin, 0.e0);
      cs_field_set_key_double(f, kscmax,  1.e0);

      /* Set the class index for the field */
      cs_field_set_key_int(f, keyccl, class_id);

      /* Set constant diffusivity for the dry air mass fraction.
       * The diffusivity used in the transport equation will be the cell value
       * of the" diffusivity_ref" for field f */
      int ifcvsl = -1;
      cs_field_set_key_int(f, kivisl, ifcvsl);

      /* Activate the drift for all scalars with key "drift" > 0 */
      int drift = CS_DRIFT_SCALAR_ON + CS_DRIFT_SCALAR_ADD_DRIFT_FLUX;

      /* Activated drift. As it is the continuous phase class (class_id = -1),
       * the convective flux is deduced for classes > 0
       * and bulk class (class_id = 0) */
      cs_field_set_key_int(f, keydri, drift);


      /* Equation parameters */
      cs_equation_param_t *eqp = cs_field_get_equation_param(f);
      /* Full upwind scheme */
      eqp->blencv = 0.0;
    }

    {
      /* Bulk humid air temperature */
      f = cs_field_by_name_try("temperature");
      /* If no temperature field exist, we create the field */
      if (f == NULL) {
        /* Thermal model - Set parameters of calculations (module optcal) */
        cs_thermal_model_t *thermal = cs_get_glob_thermal_model();

        /* Solve for temperature of bulk humid air */
        thermal->thermal_variable = CS_THERMAL_MODEL_TEMPERATURE;

        /* Temperature treated in Celsius */
        thermal->temperature_scale = CS_TEMPERATURE_SCALE_CELSIUS;
        /* The thermal transported scalar is the temperature of the bulk.
         * If the atmospheric module is switched off (i.e., iatmos!= 2)
         * , we create the field. */
        int f_id = cs_variable_field_create("temperature",
                                            "Temperature humid air",
                                            CS_MESH_LOCATION_CELLS,
                                            1);

        f = cs_field_by_id(f_id);
        cs_add_model_thermal_field_indexes(f->id);
      }

      /* Variable cp (0 = variable, -1 = constant) since it changed with humidity
       * Needs to be specified here because the automated creation and
       * initialization of the cell array for cp in 'iniva0' depends on its value
       * (unlike the cell arrays for the density and viscosity which are
       * initialized irrespective of the values of irovar and ivivar) */
      fp->icp = 0;

      /* Set variable diffusivity for the humid air temperature.
       * The diffusivity used in the transport equation will be the cell value
       * of the viscls array for f_id.
       * This value is updated at the top of each time step in 'ctphyv' along
       * with the other variable properties */
      int ifcvsl = 0;
      cs_field_set_key_int(f, kivisl, ifcvsl);

      /* Associate temperature to continuous phase -> class_id = -1 */
      cs_field_set_key_int(f, keyccl, class_id);

      /* Activate the drift versus the mixture velocity */
      int drift = CS_DRIFT_SCALAR_ON + CS_DRIFT_SCALAR_ADD_DRIFT_FLUX;

      /* Activated drift. As it is the continuous phase class (class_id = -1),
       * the convective flux is deduced for classes > 0
       * and bulk class (class_id = 0) */
      cs_field_set_key_int(f, keydri, drift);
    }
  }
  {
    /* RAIN VARIABLES
     * ============== */

    /* Associate liquid water rain with class 1 */
    int class_id = 1;

    int f_id = cs_variable_field_create("ym_l_r",
                                        "Mass frac rain",
                                        CS_MESH_LOCATION_CELLS,
                                        1);
    f = cs_field_by_id(f_id);

    /* Clipping of rain mass fraction 0 <= ym_l_r <=1 */
    cs_field_set_key_double(f, kscmin, 0.e0);
    cs_field_set_key_double(f, kscmax, 1.e0);

    /* Set the class index for the field */
    cs_field_set_key_int(f, keyccl, class_id);

    /* Scalar with drift: create additional mass flux.
     * This flux will then be reused for all scalars associated to this class
     * (here : rain liquid water variables)
     * We set the bit corresponding to drift flux computation to 1.
     * TODO : make it optional ?*/
    int drift = CS_DRIFT_SCALAR_ON + CS_DRIFT_SCALAR_ADD_DRIFT_FLUX;

    cs_field_set_key_int(f, keydri, drift);

    /* Set constant diffusivity for injected liquid mass fraction */
    int ifcvsl = -1;
    cs_field_set_key_int(f, kivisl, ifcvsl);
    cs_add_model_field_indexes(f->id);

    /* Equation parameters */
    cs_equation_param_t *eqp = cs_field_get_equation_param(f);
    /* Set beta limiter to maintain y_p in the limits */
    eqp->isstpc = 2;
    /* Full upwind scheme */
    eqp->blencv = 0.0;

    /* Transport and solve for the enthalpy of the liquid rain-with the same
     * drift as the mass y_l_r. The transported variable is actually Ylr.hlr.
     * NB : Enthalpy of the liquid must be transported after the bulk
     * enthalpy. */

    f_id = cs_variable_field_create("ymh_l_r",
                                    "ym.hl rain",
                                    CS_MESH_LOCATION_CELLS,
                                    1);
    f = cs_field_by_id(f_id);
    cs_field_set_key_int(f, keyccl, class_id);

    /* Scalar with drift, but do not create an additional mass flux for the
     * enthalpy (use ^= to reset the bit for drift flux calculation).
     * It reuses the mass flux already identified with the mass fraction. */
    drift = CS_DRIFT_SCALAR_ON;

    cs_field_set_key_int(f, keydri, drift);

    /* Set variable diffusivity for the injected liquid enthalpy transport.
     * The diffusivity used in the transport equation will be the cell value
     * of the viscls array for field f */

    ifcvsl = 0;
    cs_field_set_key_int(f, kivisl, ifcvsl);
    cs_add_model_field_indexes(f->id);

    /* Equation parameters */
    eqp = cs_field_get_equation_param(f);
    /* Full upwind scheme */
    eqp->blencv = 0.0;

    /* Variable fields creation for rain drops velocities if we want to solve
     * rain fall velocity */
    cs_ctwr_option_t *ct_opt = cs_get_glob_ctwr_option();

    if (ct_opt->solve_rain_velocity) {
      char f_name[80];
      char f_label[80];

      /* Rain drops velocities --> treated as particles */
      sprintf(f_name, "v_p_%02d", class_id);
      sprintf(f_label, "Velocity rain");
      f_id = cs_variable_field_create(f_name, f_label,
                                      CS_MESH_LOCATION_CELLS, 3);
      f = cs_field_by_id(f_id);
      cs_field_set_key_int(f, keyccl, class_id);
      cs_add_model_field_indexes(f_id);

      /* Scalar with drift, but do not create an additional mass flux */
      drift = CS_DRIFT_SCALAR_ON + CS_DRIFT_SCALAR_NO_MASS_AGGREGATION;
      cs_field_set_key_int(f, keydri, drift);

      /* Equation parameters */
      eqp = cs_field_get_equation_param(f);
      /*Full upwind scheme */
      eqp->blencv = 0.0;

      //TODO : Check equation parameters to set for v_p
      //       + what equation is exactly solved for v_p ??
      //       --> it is the equation on y_p*v_p
    }
  }

  {
    /* PACKING ZONE VARIABLES
     * ====================== */

    /* Associate injected liquid water in packing with class 2 */
    int class_id = 2;

    /* Mass of injected liquid */
    int f_id = cs_variable_field_create("y_l_packing",
                                        "Yl packing",
                                        CS_MESH_LOCATION_CELLS,
                                        1);
    f = cs_field_by_id(f_id);

    /* Clipping of packing liquid mass 0 < y_l_packing */
    cs_field_set_key_double(f, kscmin, 0.e0);

    /* Set the class index for the field */
    cs_field_set_key_int(f, keyccl, class_id);

    /* Scalar with drift: create additional mass flux.
     * This flux will then be reused for all scalars associated to this class
     * (here : injected liquid water variables in packing)
     * We set the bit corresponding to drift flux computation to 1.
     * We impose the mass flux in the packing because liquid film velocity
     * is nearly constant.
     * TODO : make it optional ?*/
    int drift = CS_DRIFT_SCALAR_ON + CS_DRIFT_SCALAR_ADD_DRIFT_FLUX
                + CS_DRIFT_SCALAR_IMPOSED_MASS_FLUX;

    cs_field_set_key_int(f, keydri, drift);

    /* Set constant diffusivity for injected liquid mass fraction */
    int ifcvsl = -1;
    cs_field_set_key_int(f, kivisl, ifcvsl);
    cs_add_model_field_indexes(f->id);

    /* Equation parameters */
    cs_equation_param_t *eqp = cs_field_get_equation_param(f);
    /* Upwind schemes for scalars in packing zone */
    eqp->blencv = 0.;
    eqp->idiff  = 0;
    eqp->idifft = 0;

    /* Do not show y_l_packing in post-processing */
    cs_field_set_key_int(f, cs_field_key_id("post_vis"), 0);

    /* Transport and solve for the enthalpy of the liquid in packing - with the
     * same drift as the mass y_l_packing in the packing zones.
     * The solved variable is actually Ylp.hlp.
     * NB : Enthalpy of the liquid must be transported after the bulk
     * enthalpy. */

    f_id = cs_variable_field_create("yh_l_packing",
                                    "Yl.hl_packing",
                                    CS_MESH_LOCATION_CELLS,
                                    1);
    /* TODO (from ctvarp.f90) : x_p_h_l or y_p_h_2 */

    f = cs_field_by_id(f_id);
    cs_field_set_key_int(f, keyccl, class_id);

    /* Scalar with drift, but do not create an additional mass flux for the
     * enthalpy (use ^= to reset the bit for drift flux calculation).
     * It reuses the mass flux already identified with the mass variable. */
    drift = CS_DRIFT_SCALAR_ON + CS_DRIFT_SCALAR_IMPOSED_MASS_FLUX;

    cs_field_set_key_int(f, keydri, drift);

    /* Set variable diffusivity for the injected liquid enthalpy transport.
     * The diffusivity used in the transport equation will be the cell value
     * of the viscls array for field f */

    ifcvsl = 0;
    cs_field_set_key_int(f, kivisl, ifcvsl);
    cs_add_model_field_indexes(f->id);

    /* Equation parameters */
    eqp = cs_field_get_equation_param(f);
    /* Upwind schemes for scalars in packing zone */
    eqp->blencv = 0.;
    eqp->idiff  = 0;
    eqp->idifft = 0;
  }

}

/*----------------------------------------------------------------------------
 * Add property fields
 *----------------------------------------------------------------------------*/

void
cs_ctwr_add_property_fields(void)
{
  cs_field_t *f;
  int class_id = 1;
  int field_type = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY;
  bool has_previous = false;
  const int klbl   = cs_field_key_id("label");
  const int keyvis = cs_field_key_id("post_vis");
  const int keylog = cs_field_key_id("log");
  const int post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;
  cs_ctwr_option_t *ct_opt = cs_get_glob_ctwr_option();

  if (cs_field_by_name_try("thermal_expansion") == NULL) {
    /* Humid air thermal expansion coefficient */
    f = cs_field_create("thermal_expansion",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Beta");
  }
  /* CONTINUOUS FIELD (HUMID AIR) PROPERTIES */
  /* NB: 'c' stands for continuous and 'p' for particles */
 {
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
  {
    /* Humidity field */
    f = cs_field_create("humidity",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Humidity");
  }

  {
    /* Saturated humidity field */
    f = cs_field_create("x_s",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Humidity sat");
  }

  {
    /* Relative humidity field */
    f = cs_field_create("x_rel",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Humidity rel");
  }


  {
    /* Humid air enthalpy field */
    f = cs_field_create("enthalpy",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Enthalpy humid air");
  }

  {
    /* Mass fraction of the continuous phase (X1) */
    f = cs_field_create("x_c",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Gas mass fraction");
  }

  {
    /* Mass fraction of the continuous phase (X1) BOUNDARY VALUE */
    f = cs_field_create("b_x_c",
                        field_type,
                        CS_MESH_LOCATION_BOUNDARY_FACES,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Boundary gas mass fraction");
  }

  {
    /* Continuous phase volume fraction */
    f = cs_field_create("vol_f_c",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Vol. frac. air");
  }

  /* Properties to create for rain velocity equation solving */
  if (ct_opt->solve_rain_velocity) {
  }


  /* RAIN FIELD PROPERTIES */
  {
    /* Rain temperature */
    f = cs_field_create("temp_l_r",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Temperature rain");
  }

  {
    /* Rain enthalpy */
    f = cs_field_create("h_l_r",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Enthalpy rain");
  }

  {
    /* Rain volume fraction */
    f = cs_field_create("vol_f_r",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Vol. frac. rain");
  }

  /* RAIN VELOCITY TRANSPORT EQUATION FOR HOMOGENEOUS MODEL */
  /* Properties to create for rain velocity equation solving */
  if (ct_opt->solve_rain_velocity) {

    /* Continuous phase drift velocity */
    f = cs_field_create("vd_c",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        3,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Drift velocity gas phase");

    /* Continuous phase velocity */
    f = cs_field_create("v_c",
        field_type,
        CS_MESH_LOCATION_CELLS,
        3,
        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Velocity continuous phase");

    char f_name[80];
    char f_label[80];

    /* Particle limit velocity */
    sprintf(f_name, "vg_lim_p_%02d", class_id);
    sprintf(f_label, "Terminal velocity rain");
    f = cs_field_create(f_name,
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        3,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, f_label);

    /* Drift velocity for rain drops */
    sprintf(f_name, "vd_p_%02d", class_id);
    sprintf(f_label, "Drift velocity rain");
    f = cs_field_create(f_name,
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        3,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, f_label);
  }

  /* LIQUID IN PACKING FIELD PROPERTIES */
  {
    /* Liquid temperature in packing */
    f = cs_field_create("temp_l_packing",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Temperature liq packing");
  }
  {
    /* Liquid enthalpy in packing */
    f = cs_field_create("h_l_packing",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Enthalpy liq packing");
  }
  {
    /* True liquid mass fraction in packing */
    f = cs_field_create("ym_l_packing",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Mass frac liq packing");
  }


  {
    /* Liquid vertical velocity in packing */
    f = cs_field_create("vertvel_l",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Velocity liq packing");
  }

  {
    /* Liquid mass flux in packing */
    f = cs_field_create("mass_flux_l",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Mass flux liq packing");
  }


  /* Mass and energy exchange terms in packing and rain for evaporation rate
   * and thermal power post-processing, they are updated
   * in cs_ctwr_source_term */
  {
    /* Evaporation rate in packing */
    f = cs_field_create("evaporation_rate_packing",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Evaporation packing");
  }

  {
    /* Evaporation rate in rain */
    f = cs_field_create("evaporation_rate_rain",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Evaporation rain");
  }

  {
    /* Thermal power in packing */
    f = cs_field_create("thermal_power_packing",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Thermal power packing");
  }

  {
    /* Thermal power in rain */
    f = cs_field_create("thermal_power_rain",
                        field_type,
                        CS_MESH_LOCATION_CELLS,
                        1,
                        has_previous);
    cs_field_set_key_int(f, keyvis, post_flag);
    cs_field_set_key_int(f, keylog, 1);
    cs_field_set_key_str(f, klbl, "Thermal power rain");
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Map fields used by the cooling tower module to pointers.
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_field_pointer_map(void)
{
  /* No need to redefine the temperature and enthalpy for humid air as they
     have already been defined in 'cs_field_pointer_map',
     which comes after 'ctvarp' */
  cs_field_pointer_map(CS_ENUMF_(humid), cs_field_by_name_try("humidity"));
  cs_field_pointer_map(CS_ENUMF_(ym_w),
                       cs_field_by_name_try("ym_water"));
  cs_field_pointer_map(CS_ENUMF_(t_l_pack),
                       cs_field_by_name_try("temp_l_packing"));
  cs_field_pointer_map(CS_ENUMF_(yh_l_pack),
                       cs_field_by_name_try("yh_l_packing"));
  cs_field_pointer_map(CS_ENUMF_(y_l_pack),
                       cs_field_by_name_try("y_l_packing"));
  cs_field_pointer_map(CS_ENUMF_(thermal_diff_h),
                       cs_field_by_name_try("thermal_conductivity"));
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
