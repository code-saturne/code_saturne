/*============================================================================
 * User source terms for a scalar in a channel example.
 *============================================================================*/

/* code_saturne version 6.2-alpha */

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
#include <math.h>

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_source_terms-scalar_in_a_channel.c
 *
 * \brief User source terms for a scalar in a channel example.
 *
 * See the reference \ref cs_user_source_terms.c for documentation.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function called at each time step to define source terms.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in]       f_id     field id of the variable
 * \param[out]      st_exp   explicit source term
 * \param[out]      st_imp   implicit part of the source term
 */
/*----------------------------------------------------------------------------*/

void
cs_user_source_terms(cs_domain_t  *domain,
                     int           f_id,
                     cs_real_t    *st_exp,
                     cs_real_t    *st_imp)
{
  CS_NO_WARN_IF_UNUSED(domain);

  /* field structure */
  const cs_field_t  *f = cs_field_by_id(f_id);

  /* local number of mesh cells */
  const cs_lnum_t  n_cells = cs_glob_mesh->n_cells;

  /* mesh quantities */
  const cs_real_t  *cell_f_vol = cs_glob_mesh_quantities->cell_vol;

  if (f != cs_thermal_model_field())
    return;

  /* velocity */
  const cs_real_3_t  *cvar_vel = (const cs_real_3_t *)(CS_F_(vel)->val);

  /* bulk mean velocity (x component) */
  cs_real_t  ubulk = 0;
  for (cs_lnum_t i = 0; i < n_cells; i++)
    ubulk += cvar_vel[i][0] * cell_f_vol[i];

  cs_parall_sum(1, CS_DOUBLE, &ubulk);  /* sum across processes if needed */

  ubulk /= cs_glob_mesh_quantities->tot_vol;

  /* Flux x Total surface / (rho Cp) */
  cs_real_t tot_flux = 1.;

  for (cs_lnum_t i = 0; i < n_cells; i++) {
    st_imp[i] = 0.;
    st_exp[i] = cell_f_vol[i] * cvar_vel[i][0] * tot_flux / ubulk;
  }
}
