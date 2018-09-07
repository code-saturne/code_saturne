/*============================================================================
 * Log field and other array statistics at relevant time steps.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include "cs_base.h"
#include "cs_boundary_zone.h"
#include "cs_ctwr.h"
#include "cs_fan.h"
#include "cs_field.h"
#include "cs_log.h"
#include "cs_mesh_quantities.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_sles.h"
#include "cs_sles_default.h"
#include "cs_stokes_model.h"
#include "cs_thermal_model.h"
#include "cs_time_moment.h"
#include "cs_turbomachinery.h"
#include "cs_rad_transfer_options.h"
#include "cs_rotation.h"
#include "cs_turbulence_model.h"
#include "cs_lagr_log.h"
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
_log_global_model_options(void)
{
  /* Mesh quantity options */

  cs_mesh_quantities_log_setup();

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "Physical model options\n"
                  "----------------------\n"));

  /* Physical properties */

  cs_physical_constants_log_setup();

  cs_fluid_properties_log_setup();

  /* TODO : Add diftl0 printing */

  /* Thermal model */

  cs_thermal_model_log_setup();

  /* Turbulence */
  cs_turb_model_log_setup();

  cs_turb_constants_log_setup();

  cs_time_step_log_setup();

  /* Stokes model*/
  cs_stokes_model_log_setup();

  /* TODO : Partie iroext etc... */

  /* Face viscosity */
  cs_space_disc_log_setup();

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
 * \brief Log setup options.
 */
/*----------------------------------------------------------------------------*/

void
cs_log_setup(void)
{
  cs_field_log_defs();
  cs_field_log_key_defs();
  cs_field_log_all_key_vals(false);

  cs_time_moment_log_setup();

  cs_sles_default_setup();

  _log_global_model_options();

  cs_rad_transfer_log_setup();

  cs_lagr_log_setup();

  cs_fan_log_setup();

  cs_ctwr_log_setup();

  cs_log_printf_flush(CS_LOG_SETUP);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
