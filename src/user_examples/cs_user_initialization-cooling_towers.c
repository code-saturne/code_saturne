/*============================================================================
 * User initialization prior to solving time steps.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_initialization-cooling_towers.c
 *
 * \brief Initialization prior to solving time steps.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_initialization.c
 *
 * \brief Initialize variables.
 *
 * This function is called at beginning of the computation
 * (restart or not) before the time step loop.
 *
 * This is intended to initialize or modify (when restarted)
 * variable and time step values.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_initialization(cs_domain_t     *domain)
{
  /*![init]*/

  const cs_lnum_t n_cells = domain->mesh->n_cells;

  /* Apply only at the true computation start, not on restarts */
  if (domain->time_step->nt_prev > 0)
    return;

  /* Map field arrays */
  cs_real_3_t *vel = (cs_real_3_t *)CS_F_(vel)->val;
  cs_real_t *cvar_temp = cs_thermal_model_field()->val;
  cs_real_t *cpro_humid = cs_field_by_name("humidity")->val;

  /* Initialize temperature of humid air at 11 deg Celsius
   * and of humidity at 0.0063 */

   cs_array_set_value_real(n_cells, 1, 11.0, cvar_temp);
   cs_array_set_value_real(n_cells, 1, 0.0063, cpro_humid);

   /* Initialize temperature of humid air at 20 deg Celsius
    * and of humidity at 0.012
    * and of velocity at 0.5 m/s
    * for cells of zone "6" */

   const cs_zone_t *zn = cs_volume_zone_by_name("6");

   for (cs_lnum_t ii = 0; ii < zn->n_elts; ii++) {

      const cs_lnum_t c_id = zn->elt_ids[ii];

      vel[c_id][0] = 0.5;
      cvar_temp[c_id] = 20.;
      cpro_humid[c_id] = 0.012;
   }

   /*![init]*/
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
