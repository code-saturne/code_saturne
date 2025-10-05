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
  /* mesh quantities */
  const cs_lnum_t  n_cells = cs_glob_mesh->n_cells;
  const cs_real_t  *cell_f_vol = cs_glob_mesh_quantities->cell_vol;
  /*! [st_meta] */

  /* Define a cpro_rom pointer to the density */
  /*! [dens_array_3] */
  const cs_real_t  *cpro_rom = CS_F_(rho)->val;
  /*! [dens_array_3] */

  /*! [current_turb_3] */
  /* Define pointer f to the current turbulent variable field */
  const cs_field_t  *f = cs_field_by_id(f_id);

  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);

  /*! [current_turb_3] */

  /* Example of arbitrary source term for turbulence models
   * (Source term on the TKE 'k' here)
   *
   *  Source term for cvar_var:
   *  rho cell_f_vol d(cvar_var)/dt       = ...
   *                 ... - rho*cell_f_vol*ff - rho*cell_f_vol*cvar_var/tau
   *
   *  With ff=3.0 and tau = 4.0
   */

  /*! [rem_code_3] */
  if (f == CS_F_(k)) {

    if (eqp->verbosity >= 1)
      bft_printf(" User source terms for turbulence variable %s\n",
                 cs_field_get_label(f));

    const cs_real_t ff  = 3.0;
    const cs_real_t tau = 4.0;

    for (cs_lnum_t i = 0; i < n_cells; i++) {

      /* explicit source terms */
      st_exp[i] = -cpro_rom[i] * cell_f_vol[i] * ff;

      /* implicit source terms */
      st_imp[i] = -cpro_rom[i] * cell_f_vol[i] / tau;
    }
  }
  /*! [rem_code_3] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
