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
 * \file cs_user_initialization-gas_3ptchem.c
 *
 * \brief Initialization prior to solving time steps.
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
  /*![init]*/
  /* If this is restarted computation, do not reinitialize values */
  if (domain->time_step->nt_prev > 0)
    return;

  const cs_lnum_t n_cells = domain->mesh->n_cells;
  const cs_combustion_model_t *cm = cs_glob_combustion_model;

  cs_real_t *cvar_fm = CS_F_(fm)->val;
  cs_real_t *cvar_scalt = cs_thermal_model_field()->val;

  const cs_real_t fs = cm->gas.fs[0];

  /* Mean Mixture Fraction */
  cs_array_set_value_real(n_cells, 1, fs, cvar_fm);

  /* Enthalpy */

  if (cs_glob_physical_model_flag[CS_COMBUSTION_3PT] == 1) {
    const cs_real_t hinoxy = cm->hinoxy;
    const cs_real_t hinfue = cm->fuel.hinfue;

    cs_real_t h_ini = hinfue*fs + hinoxy*(1.0-fs);
    cs_array_set_value_real(n_cells, 1, h_ini, cvar_scalt);
  }
  /*![init]*/
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
