/*============================================================================
 * General-purpose user-defined functions called before time stepping, at
 * the end of each time step, and after time-stepping.
 *
 * These can be used for operations which do not fit naturally in any other
 * dedicated user function.
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_extra_operations.c
 *
 * \brief This function is called at the end of each time step, and has a very
 * general purpose (i.e. anything that does not have another dedicated
 * user function)
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is called at the end of each time step.
 *
 * It has a very general purpose, although it is recommended to handle
 * mainly postprocessing or data-extraction type operations.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations(cs_domain_t     *domain)
{
  CS_UNUSED(domain);

  /* ! [vorticity_d] */
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  cs_real_33_t *gradv;
  /* ! [vorticity_d] */

  /* ! [vorticity_a] */
  BFT_MALLOC(gradv, n_cells_ext, cs_real_33_t);
  /* ! [vorticity_a] */

  /* ! [vorticity_g] */
  bool use_previous_t = false;
  int inc = 1;
  cs_field_gradient_vector(CS_F_(vel),
                           use_previous_t,
                           inc,
                           gradv);
  /* ! [vorticity_g] */

  /* ! [vorticity_f] */
  cs_field_t *vort = cs_field_by_name_try("vorticity");
  /* ! [vorticity_f] */

  /* ! [vorticity_cv] */
  if (vort != NULL) {
    for (cs_lnum_t i = 0; i < n_cells; i++) {
      vort->val[i] = gradv[i][1][0] - gradv[i][0][1];
    }
  }
  /* ! [vorticity_cv] */

  /* ! [vorticity_da] */
  BFT_FREE(gradv);
  /* ! [vorticity_da] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
