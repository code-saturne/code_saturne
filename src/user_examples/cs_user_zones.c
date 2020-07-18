/*============================================================================
 * Boundary and volume zone definitions.
 *============================================================================*/

/* VERS */

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define volume and surface zones.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_zones(void)
{
  /* Example:
     a volume zone for head losses based on geometrical criteria
     described by a rule. */

  /*! [user_zones_head_loss_1] */
  {
    cs_volume_zone_define("head_loss_1",
                          "box[4.0, 6.0, -1e6, 2.0, 8.0, 1e6]",
                          CS_VOLUME_ZONE_HEAD_LOSS);
  }
  /*! [user_zones_head_loss_1] */

  /* Example:
     define 3 spherical source term zones based on geometrical criteria
     described by a rule. */

  /*! [user_zones_volume_1] */
  {
    char name[128], criteria[128];

    for (int i = 0; i < 3; i++) {

      double s_coords[] = {0, 0, 2.0*i};

      snprintf(name, 127, "source_%d", i);
      snprintf(criteria, 127, "sphere[%f, %f, %f, 0.5]",
               s_coords[0], s_coords[1], s_coords[2]);

      cs_volume_zone_define(name, criteria, CS_VOLUME_ZONE_SOURCE_TERM);
    }
  }
  /*! [user_zones_volume_1] */

   /* Example:
     Define zones from a text file defining turbines modelled
     with the actuator disk approach. */

  {
    char name[128], criteria[128];

    FILE* file = fopen("turbines", "rt");
    int n_turbines = 0;

    /* Some parameters */
    cs_real_t turb_lenght = 1.;
    cs_real_t radius = 126./2.;

    cs_real_t wind_dir = 0.;

    if (fscanf(file, "%d\n", &n_turbines) != 1)
      bft_error(__FILE__,__LINE__, 0, _("Could not read the number of turbines."));

    for (int i = 0; i < n_turbines; i++) {
      int turbine_id = 0;
      if (fscanf(file, "%d\n", &turbine_id) != 1)
        bft_error(__FILE__,__LINE__, 0, _("Could not read turbine %d."), i);

      float turb_cen[3];
      if (fscanf(file, "%f\n", &(turb_cen[0])) != 1)
        bft_error(__FILE__,__LINE__, 0, _("Could not read turbine x."));
      if (fscanf(file, "%f\n", &(turb_cen[1])) != 1)
        bft_error(__FILE__,__LINE__, 0, _("Could not read turbine y."));
      if (fscanf(file, "%f\n", &(turb_cen[2])) != 1)
        bft_error(__FILE__,__LINE__, 0, _("Could not read turbine z."));

      double s_coords[] = {
	      turb_cen[0] - 0.5 * turb_lenght * cos(wind_dir),
	      turb_cen[1] - 0.5 * turb_lenght * sin(wind_dir),
	      turb_cen[2]};
      double e_coords[] = {
	      turb_cen[0] + 0.5 * turb_lenght * cos(wind_dir),
	      turb_cen[1] + 0.5 * turb_lenght * sin(wind_dir),
	      turb_cen[2]};

      snprintf(name, 127, "turbine_%d", i);
      snprintf(criteria, 127, "cylinder[%f, %f, %f, %f, %f, %f, %f]",
               s_coords[0], s_coords[1], s_coords[2],/* Start */
               e_coords[0], e_coords[1], e_coords[2],/* End */
               radius /* Radius */
	       );

      cs_volume_zone_define(name, criteria, CS_VOLUME_ZONE_SOURCE_TERM);
    }
    if (fclose(file) != 0)
      bft_error(__FILE__,__LINE__, 0, _("Could not close the file."));
  }

  /* Example:
     define simple boundary zones, allowing all faces not in the
     "INLET" or "OUTLET" groups to be considered as walls */

  /*! [user_zones_boundary_1] */
  {
    int z_id = cs_boundary_zone_define("wall", "all[]", 0);
    cs_boundary_zone_set_overlay(z_id, true);

    cs_boundary_zone_define("inlet", "INLET", 0);
    cs_boundary_zone_define("outlet", "OUTLET", 0);
  }
  /*! [user_zones_boundary_1] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
