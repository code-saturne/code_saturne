/*============================================================================
 * Log field and other array statistics at relevant time steps.
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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_ale.h"
#include "cs_atmo.h"
#include "cs_base.h"
#include "cs_boundary.h"
#include "cs_boundary_zone.h"
#include "cs_cf_model.h"
#include "cs_ctwr.h"
#include "cs_combustion_model.h"
#include "cs_domain.h"
#include "cs_fan.h"
#include "cs_field.h"
#include "cs_function.h"
#include "cs_log.h"
#include "cs_mesh_quantities.h"
#include "cs_mobile_structures.h"
#include "cs_notebook.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_restart.h"
#include "cs_sles.h"
#include "cs_sles_default.h"
#include "cs_thermal_model.h"
#include "cs_time_moment.h"
#include "cs_turbomachinery.h"
#include "cs_rad_transfer_options.h"
#include "cs_rotation.h"
#include "cs_turbulence_model.h"
#include "cs_lagr_log.h"
#include "cs_velocity_pressure.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_log_setup.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_log_setup.c

  \brief Setup info at the end of the setup stage.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type and macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Log various global model options.
 *----------------------------------------------------------------------------*/

static void
_log_error_estimators(void)
{
  int ee_count = 0;

  const char *name[] = {"est_error_pre_2",
                        "est_error_der_2",
                        "est_error_cor_2",
                        "est_error_tot_2"};

  const char *desc[] = {"prediction",
                        "drift",
                        "correction",
                        "total"};

  for (int i = 0; i < 4; i++) {

    const cs_field_t *f = cs_field_by_name_try(name[i]);
    if (f != NULL) {
      if (ee_count == 0)
        cs_log_printf(CS_LOG_SETUP,
                      _("\n"
                        "Error estimators for Navier-Stokes\n"
                        "----------------------------------\n\n"));

      cs_log_printf(CS_LOG_SETUP,
                    _("  %s: %s\n"), name[i], desc[i]);

      ee_count += 1;
    }
  }
}

/*----------------------------------------------------------------------------
 * Log various global model options.
 *----------------------------------------------------------------------------*/

static void
_log_global_model_options(void)
{
  /* Mesh quantity options */

  cs_mesh_quantities_log_setup();

  /* Notebook parameters */

  cs_notebook_log_setup();

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "Physical model options\n"
                  "----------------------\n"));

  /* Physical properties */

  cs_physical_constants_log_setup();
  cs_fluid_properties_log_setup();

  /* Thermal model */

  cs_thermal_model_log_setup();

  /* Turbulence */

  cs_turb_model_log_setup();
  cs_turb_constants_log_setup();

  /* Time discretization */

  cs_time_step_log_setup();
  cs_time_scheme_log_setup();

  /* Velocity-pressure coupling */

  cs_velocity_pressure_model_log_setup();
  cs_velocity_pressure_param_log_setup();

  _log_error_estimators();

  /* Compressible model */

  cs_cf_model_log_setup();

  /* Atmospheric */

  cs_atmo_log_setup();

  /* Atmospheric chemistry */

  cs_atmo_chemistry_log_setup();

  /* Atmospheric aerosols */

  cs_atmo_aerosol_log_setup();

  /* Combustion */

  cs_combustion_log_setup();

  /* TODO: iroext, etc... */

  /* Face viscosity */

  cs_space_disc_log_setup();

  /* ALE */

  if (cs_glob_ale > CS_ALE_NONE)
    cs_mobile_structures_log_setup();

  /* Rotation info */

  if (cs_turbomachinery_get_model() == CS_TURBOMACHINERY_NONE) {
    const cs_rotation_t  *r = cs_glob_rotation;

    cs_log_printf(CS_LOG_SETUP, _("\nSubdomain rotation\n"
                                  "------------------\n\n"));

    cs_log_printf(CS_LOG_SETUP,
                  _("  Global domain rotation:\n"
                    "    axis:             [%g, %g, %g]\n"
                    "    invariant point:  [%g, %g, %g]\n"
                    "    angular velocity:  %g radians/s\n"),
                  r->axis[0], r->axis[1], r->axis[2],
                  r->invariant[0], r->invariant[1], r->invariant[2],
                  r->omega);

  }

  /* Zone information */

  cs_volume_zone_log_setup();
  cs_boundary_zone_log_setup();

  /* BC information */

  cs_boundary_log_setup(cs_glob_domain->boundaries);
  cs_boundary_log_setup(cs_glob_domain->ale_boundaries);
}


/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log setup options and define the default setup for SLES.
 */
/*----------------------------------------------------------------------------*/

void
cs_log_setup(void)
{
  cs_field_log_defs();
  cs_field_log_key_defs();
  cs_field_log_all_key_vals(false);

  cs_time_moment_log_setup();

  cs_function_log_defs();
  cs_function_log_all_settings();

  cs_sles_default_setup();

  cs_restart_log_setup();
  cs_log_printf(CS_LOG_SETUP,
                _("  read auxiliary:       %d\n"),
                cs_glob_restart_auxiliary->read_auxiliary);
  cs_log_printf(CS_LOG_SETUP,
                _("  write auxiliary:      %d\n"),
                cs_glob_restart_auxiliary->write_auxiliary);

  _log_global_model_options();

  cs_rad_transfer_log_setup();

  cs_lagr_log_setup();

  cs_fan_log_setup();

  cs_ctwr_log_setup();

  cs_log_printf_flush(CS_LOG_SETUP);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
