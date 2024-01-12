/*============================================================================
 * Lagrangian volume injection definitions.
 *============================================================================*/

/* VERS */

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

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
 * based on predefined volume zones (\ref cs_zone_t).
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_volume_conditions(void)
{
  /*! [lagr_vol_cond_variables] */
  cs_lagr_zone_data_t *lagr_vol_conds = cs_lagr_get_volume_conditions();
  /*! [lagr_vol_cond_variables] */

  /* Example for a uniform injection at computation initialization */

  /*! [lagr_vol_define_injection_1] */

  {
    /* The volume zone named "particle_injection" is created in the GUI */

    const cs_zone_t *z = cs_volume_zone_by_name("particle_injection");

    /* Inject 1 particle set every time step */

    int set_id = 0;

    cs_lagr_injection_set_t *zis
      = cs_lagr_get_injection_set(lagr_vol_conds, z->id, set_id);

    zis->n_inject = 1000;
    zis->injection_frequency = 1; /* if <= 0, injection at
                                     initialization only */
    zis->velocity_profile = -1; /* fluid velocity */

    zis->stat_weight = 1.0;

    zis->diameter = 5e-6;
    zis->diameter_variance = 1e-6;

    zis->density = 2475.;

    zis->fouling_index = 100.0;

    /* Activation of agglomeration */
    if (cs_glob_lagr_model->agglomeration == 1 ||
        cs_glob_lagr_model->fragmentation == 1 ) {
      zis->aggregat_class_id = 1;
      zis->aggregat_fractal_dim = 3.;
      zis->diameter = cs_glob_lagr_agglomeration_model->base_diameter
                      * pow((cs_real_t)zis->aggregat_class_id,
                            1./zis->aggregat_fractal_dim);
    }

  }
  /*! [lagr_vol_define_injection_1] */


  /* Example for a uniform injection at computation initialization */

  /*! [lagr_vol_define_injection_2] */

  {
    /* The volume zone containing all cells always has id 0;
       a given zone may otherwise be selected using cs_volume_zone_by_name() */

    const cs_zone_t *z = cs_volume_zone_by_id(0);

    /* Inject 2 particle sets of different diameters */

    cs_gnum_t n_inject[] = {500, 500};

    cs_real_t diam[] = {1e-3, 1e-2};
    cs_real_t diam_dev[] = {1e-6, 1e-5};
    cs_real_t density[] = {2500., 1500.};

    for (int set_id = 0; set_id < 2; set_id++) {

      cs_lagr_injection_set_t *zis
        = cs_lagr_get_injection_set(lagr_vol_conds, z->id, set_id);

      zis->n_inject = n_inject[set_id];
      zis->injection_frequency = 0; /* if <= 0, injection at
                                       initialization only */

      if (cs_glob_lagr_model->n_stat_classes > 0)
        zis->cluster = set_id + 1;

      zis->velocity_profile = -1; /* fluid velocity */

      zis->stat_weight = 1.0;
      zis->flow_rate   = 0;

      zis->diameter = diam[set_id];
      zis->diameter_variance = diam_dev[set_id];

      zis->density = density[set_id];
      zis->fouling_index = 100.0;

      zis->temperature_profile = 0; /* fluid temperature */

      zis->cp = 1400.0;
      zis->emissivity = 0.7;

    }
  }

  /*! [lagr_vol_define_injection_2] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
