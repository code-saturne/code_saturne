/*============================================================================
 * Boundary and volume zone definitions.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_mesh_location.h"
#include "cs_base.h"
#include "cs_boundary_zone.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_log.h"
#include "cs_parameters.h"
#include "cs_time_step.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

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
               s_coords[0], s_coords[0], s_coords[0]);

      cs_volume_zone_define(name, criteria, CS_VOLUME_ZONE_SOURCE_TERM);
    }
  }
  /*! [user_zones_volume_1] */

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
