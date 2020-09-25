/*============================================================================
 * User initialization prior to solving time steps.
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
  const cs_mesh_t *m = domain->mesh;

  cs_field_t *f = CS_F_(k);

  cs_real_t c0 = 2.1; /* To be coherent with the Lagrangian module */
  cs_real_t cl = 1. / (0.5 + 0.75 * c0);

  cs_real_t tlag_x = cs_notebook_parameter_value_by_name("tlag_x");
  cs_real_t tlag_y = cs_notebook_parameter_value_by_name("tlag_y");
  cs_real_t tlag_z = cs_notebook_parameter_value_by_name("tlag_z");
  cs_real_t bx     = cs_notebook_parameter_value_by_name("bx");
  
  cs_real_t tlag = (tlag_x + tlag_y + tlag_z)/3.0;

  /* We want bx = sqrt(c0 eps) so eps = bx^2/c0 */
  cs_real_t eps = bx * bx / c0;

  /* We want tlag = cl k / eps  so k = tlag / cl * eps */
  cs_real_t tke = tlag * eps / cl;

  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {
    CS_F_(k)->val[cell_id] = tke;
    CS_F_(eps)->val[cell_id] = eps;
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
