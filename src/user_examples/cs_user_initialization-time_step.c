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
 * \file cs_user_initialization-time_step.c
 *
 * \brief Time step example
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
  /*
    Time step modification

    We do a computation restart with an adaptive (variable in time and constant
    in space) or local (variable in time and space) time step.
    We want to modify the time step read from the restart file
    (in order to overcome a too slow evolution for instance).
  */

  /*! [init] */
  if (! cs_restart_present())      /* Useful only for a restart */
    return;

  const cs_lnum_t n_cells = domain->mesh->n_cells;

  cs_field_t *f = cs_field_by_name("dt");

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    f->val[cell_id] *= 10.;
  }
  /*! [init] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
