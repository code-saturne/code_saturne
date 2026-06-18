/*============================================================================
 * User initialization prior to solving time steps.
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

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * Initialize variables.
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
cs_user_initialization([[maybe_unused]] cs_domain_t  *domain)
{
  /*! [loc_var_dec] */
  const cs_lnum_t n_cells = domain->mesh->n_cells;
  const cs_real_3_t *cell_cen
    = (const cs_real_3_t *)domain->mesh_quantities->cell_cen;

  const cs_real_t hhot = 0.4075;
  /*! [loc_var_dec] */

  /*! [init2] */
  /* Apply only at the true computation start, not on restarts */
  if (domain->time_step->nt_prev > 0)
    return;

  /* Mass fraction = 1 for first gas, 0 for others */

  int n_gasses = cs_glob_elec_properties->n_gas;

  if (n_gasses > 1) {

    for (int sp_id = 0; sp_id < n_gasses-1; sp_id++) {
      auto cvar_ycoel = CS_FI_(ycoel, sp_id)->get_val_s();
      if (sp_id == 0)
        cvar_ycoel.set_to_val(1.0, n_cells);
      else
        cvar_ycoel.set_to_val(0.0, n_cells);
    }
  }

  /* Enthalpy = H(T0) or 0

     For electric arc,
     for the whole compution domain enthalpy is set to H(T) using
     the model conversion function.

     For Joule jeating by direct conduction,
     enthalpy is set to zero, and the user will enter his H(T) function
     tabulation.
  */

  /* Enthaly reinitialization */

  auto cpro_t = CS_F_(t)->get_val_s();
  auto cvar_h = CS_F_(h)->get_val_s();

  if (cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] >= 1) {
    cs_elec_convert_t_to_h_cells(cpro_t, cvar_h);
  }
  else {
    cs_user_physical_properties_t_to_h(domain,
                                       cs_volume_zone_by_id(0),
                                       false,
                                       cpro_t,
                                       cvar_h);
  }

  /* Enthalpy modification in given area */

  for (cs_lnum_t i = 0; i < n_cells; i++) {
    cs_real_t rayo = sqrt(  cell_cen[i][0]*cell_cen[i][0]
                          + cell_cen[i][1]*cell_cen[i][1]);
    if (rayo <= 0.8e-3)
      cvar_h[i] = hhot;
  }

  /* Set electric potential to 0 (real component and imaginary component -- for
     Joule heating by direct conduction) */

  auto cvar_potr = CS_F_(potr)->get_val_s();
  cvar_potr.set_to_val(0.0, n_cells);

  if (CS_F_(poti) != nullptr) {
    auto cvar_poti = CS_F_(poti)->get_val_s();
    cvar_poti.set_to_val(0.0, n_cells);
  }

  /* Vector potential (electric arcs, 3d) */

  if (CS_F_(potva) != nullptr) {
    auto cvar_potva = CS_F_(potva)->get_val_v();
    cvar_potva.set_to_val(0.0, 3*n_cells);
  }
  /*! [init2] */
}

/*----------------------------------------------------------------------------*/
