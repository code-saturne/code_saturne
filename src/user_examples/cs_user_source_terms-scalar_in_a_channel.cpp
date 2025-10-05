/*============================================================================
 * Additional user-defined source terms for variable equations.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Standard library headers
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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * Additional user-defined source terms for variable equations
 * (momentum, scalars, turbulence...).
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in]       f_id     field id of the variable
 * \param[out]      st_exp   explicit source term
 * \param[out]      st_imp   implicit part of the source term
 */
/*----------------------------------------------------------------------------*/

void
cs_user_source_terms
(
  [[maybe_unused]] cs_domain_t  *domain,
  [[maybe_unused]] int           f_id,
  [[maybe_unused]] cs_real_t    *st_exp,
  [[maybe_unused]] cs_real_t    *st_imp
)
{
  /*! [st_meta] */
  /* field structure */
  const cs_field_t  *f = cs_field_by_id(f_id);

  /* local number of mesh cells */
  const cs_lnum_t  n_cells = cs_glob_mesh->n_cells;

  /* mesh quantities */
  const cs_real_t  *cell_f_vol = cs_glob_mesh_quantities->cell_vol;
  /*! [st_meta] */

  /*! [thermal_scalar_only] */
  if (f != cs_thermal_model_field())
    return;
  /*! [thermal_scalar_only] */

  /*! [map_fields] */
  /* velocity */
  const cs_real_3_t  *cvar_vel = (const cs_real_3_t *)(CS_F_(vel)->val);
  /*! [map_fields] */

  /*! [bulk_mean_velocity] */
  /* bulk mean velocity (x component) */
  cs_real_t  ubulk = 0;
  for (cs_lnum_t i = 0; i < n_cells; i++)
    ubulk += cvar_vel[i][0] * cell_f_vol[i];

  cs_parall_sum(1, CS_REAL_TYPE, &ubulk);  /* sum across processes if needed */

  ubulk /= cs_glob_mesh_quantities->tot_vol;
  /*! [bulk_mean_velocity] */

  /*! [scalar_st] */
  /* Flux x Total surface / (rho Cp) */
  cs_real_t tot_flux = 1.;

  for (cs_lnum_t i = 0; i < n_cells; i++) {
    st_imp[i] = 0.;
    st_exp[i] = cell_f_vol[i] * cvar_vel[i][0] * tot_flux / ubulk;
  }
  /*! [scalar_st] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
