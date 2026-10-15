/*============================================================================
 * User definition of boundary conditions.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include "cs_headers.h"

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * User definition of boundary conditions.
 *
 * Configure the integral wall function on the four cavity walls using the
 * characteristic lengths and activation thresholds of the Sebilleau case.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in, out]  bc_type  boundary face types
 */
/*----------------------------------------------------------------------------*/

void
cs_user_boundary_conditions([[maybe_unused]] cs_domain_t  *domain,
                            [[maybe_unused]] int           bc_type[])
{
  const char *zone_names[] = {
    "X0",
    "X1",
    "Z0",
    "Z1"
  };
  const cs_real_t characteristic_lengths[] = {1., 1., 0.5, 0.5};

  cs_wall_functions_iwf_options_t options = {
    cs_notebook_parameter_value_by_name("IWF_Ri_min"),
    cs_notebook_parameter_value_by_name("IWF_y_cross_min"),
    1.e-4, /* tolerance */
    1000,  /* max_iterations */
    4      /* first_time_step */
  };

  cs_wall_functions_iwf_set_options(&options);
  cs_wall_functions_iwf_reset_faces();

  for (int zone_id = 0; zone_id < 4; zone_id++) {
    const cs_zone_t *z = cs_boundary_zone_by_name(zone_names[zone_id]);
    cs_wall_functions_iwf_set_faces(z->n_elts,
                                    z->elt_ids,
                                    characteristic_lengths[zone_id]);
  }
}

/*----------------------------------------------------------------------------*/
