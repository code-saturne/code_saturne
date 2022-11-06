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
 * \file cs_user_initialization-base.c
 *
 * \brief Initialization prior to solving time steps.
 *        Basic examples
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
  /*! [init] */

  const cs_lnum_t n_cells = domain->mesh->n_cells;

  /* If this is restarted computation, do not reinitialize values */
  if (domain->time_step->nt_prev > 0)
    return;

  /* Initialize "scalar1" field to 25 only if it exists  */
  cs_field_t *fld = cs_field_by_name_try("scalar1");

  if (fld != NULL) {
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
      fld->val[cell_id] = 25;
  }

  /* In the case of the EBU pre-mixed flame module the user can initialise (by example 25):
   * the mixing rate
   * the fresh gas mass fraction
   * the mixture enthalpy */

  cs_real_t *cvar_fm = CS_F_(fm)->val;
  cs_real_t *cvar_ygfm = CS_F_(ygfm)->val;
  cs_real_t *cvar_scalt = cs_thermal_model_field()->val;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cvar_fm[c_id] = 25;
    cvar_ygfm[c_id] = 25;
    cvar_scalt[c_id] = 25;
  }
  /*! [init] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
