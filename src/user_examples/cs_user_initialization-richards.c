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
 * \file cs_user_initialization-richards.c
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
  /* Apply only at the true computation start, not on restarts */
  if (domain->time_step->nt_prev > 0)
    return;

  const cs_lnum_t n_cells = domain->mesh->n_cells;
  const cs_real_3_t *cell_cen
    = (const cs_real_3_t *)domain->mesh_quantities->cell_cen;

  /*![richards_init_cell]*/
  /* Initialization of the hydraulic pressure H with a contant gradient of 1
   * along z axis and -0.01 along x axis */

  cs_real_t *cvar_pr = CS_F_(p)->val;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cvar_pr[c_id] = 1.0*cell_cen[c_id][2] - 1.e-2*cell_cen[c_id][0];
  }
  /*![richards_init_cell]*/

  /*![richards_init_grp]*/
  cs_real_t *cvar_scal_1 = cs_field_by_name("scalar1")->val;

  /*  Initialization per group of volume */
  const cs_zone_t *zn = cs_volume_zone_by_name("SOURCE");

  for (cs_lnum_t ii = 0; ii < zn->n_elts; ii++) {
    const cs_lnum_t c_id = zn->elt_ids[ii];
    cvar_scal_1[c_id] = 1.0;
  }
  /*![richards_init_grp]*/

  /*![richards_init_sorb]*/
  /* Index field of sorbed concentration */
  const int keysrb = cs_field_key_id("gwf_sorbed_concentration_id");
  const int isorb = cs_field_get_key_int(cs_field_by_name("scalar1"), keysrb);
  cs_real_t *cpro_sorb = cs_field_by_id(isorb)->val;

  /* no initial contamination of sorbed phase */
  cs_array_set_value_real(n_cells, 1, 0., cpro_sorb);
  /*![richards_init_sorb]*/

  /*![richards_init_precip]*/
  /* Index field of precipitated concentration */
  const int keypre = cs_field_key_id("gwf_precip_concentration_id");
  const int igwfpr = cs_field_get_key_int(cs_field_by_name("scalar1"), keypre);
  cs_real_t *cpro_precip = cs_field_by_id(igwfpr)->val;

  /* no initial precipitation phase */
  cs_array_set_value_real(n_cells, 1, 0., cpro_precip);
  /*![richards_init_precip]*/
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
