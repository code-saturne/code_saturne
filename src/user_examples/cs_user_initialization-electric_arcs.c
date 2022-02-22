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
 * \file cs_user_initialization-electric-arcs.c
 *
 * \brief Electric arcs fields initialization example.
 *
 * See \ref cs_user_initialization for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
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
  /*! [loc_var_dec] */
  const cs_lnum_t n_cells = domain->mesh-> n_cells;
  const cs_real_3_t *cell_cen
    = (const cs_real_3_t *)domain->mesh_quantities->cell_cen;

  const cs_real_t hhot = 0.4075;
  /*! [loc_var_dec] */

  /*! [init2] */
  /* Apply only at the true computation start, not on restarts */
  if (domain->time_step->nt_prev > 0)
    return;

  /* Mass fraction = 1 for first gas, 0 for others */

  int n_gasses = cs_glob_elec_properties->ngaz;

  if (n_gasses > 1) {
    cs_array_set_value_real(n_cells, 1, 1.0, CS_FI_(ycoel, 0)->val);

    for (int sp_id = 1; sp_id < n_gasses-1; sp_id++) {
      cs_array_set_value_real(n_cells, 1, 0.0, CS_FI_(ycoel, sp_id)->val);
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

  cs_real_t  *cpro_t = CS_F_(t)->val;
  cs_real_t  *cvar_h = CS_F_(h)->val;

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

  /* Set electric potential to 0 */

  /* real component */
  cs_array_set_value_real(n_cells, 1, 0.0, CS_F_(potr)->val);

  /* imaginary component (for Joule heating by direct conduction) */
  if (CS_F_(poti) != NULL)
    cs_array_set_value_real(n_cells, 1, 0.0, CS_F_(poti)->val);

  /* Vector potential (electric arcs, 3d) */
  if (CS_F_(potva) != NULL)
    cs_array_set_value_real(n_cells, 3, 0.0, CS_F_(potva)->val);
  /*! [init2] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
