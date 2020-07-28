/*============================================================================
 * Lagrangian volume injection definitions.
 *============================================================================*/

/* Code_Saturne version 5.1 */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_volume_zone.h"
#include "cs_math.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_prototypes.h"
#include "cs_random.h"

#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_log.h"

#include "cs_lagr.h"
#include "cs_lagr_new.h"
#include "cs_lagr_tracking.h"
#include "cs_lagr_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define particle volume conditions.
 *
 * This is used for the definition of volume injections,
 * based on predefined volume zones (\ref cs_volume_zone_t).
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_volume_conditions(void)
{
  cs_lagr_zone_data_t *lagr_vol_conds = cs_lagr_get_volume_conditions();

  /* The volume zone containing all cells always has id 0;
     a given zone may otherwise be selected using cs_volume_zone_by_name() */

  const cs_volume_zone_t *z = cs_volume_zone_by_name("particle_injection");

  /* Inject particle set */

  int set_id = 0;

  cs_lagr_injection_set_t *zis
    = cs_lagr_get_injection_set(lagr_vol_conds, z->id, set_id);

  zis->n_inject = 1000;
  zis->injection_frequency = 1;

  zis->velocity_profile = -1; /* fluid velocity */

  zis->stat_weight = 1.0;

  zis->diameter = 5e-6;
  zis->diameter_variance = 1e-6;

  zis->density = 2475.;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
