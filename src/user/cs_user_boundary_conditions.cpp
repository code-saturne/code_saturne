/*============================================================================
 * User definition of boundary conditions.
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
#include <stdio.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * Setup boundary conditions to be applied.
 *
 * This function is called just before \ref cs_user_finalize_setup, and
 * boundary conditions can be defined in either of those functions,
 * depending on whichever is considered more readable or practical for a
 * given use.
 *
 * \param[in, out]  domain  pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_boundary_conditions_setup
void
cs_user_boundary_conditions_setup([[maybe_unused]] cs_domain_t  *domain)
{
}

/*----------------------------------------------------------------------------*/
/*
 * User definition of boundary conditions.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in, out]  bc_type  boundary face types
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_boundary_conditions
void
cs_user_boundary_conditions([[maybe_unused]] cs_domain_t  *domain,
                            [[maybe_unused]] int           bc_type[])
{
}

/*----------------------------------------------------------------------------*/
/*
 * User definition of boundary conditions for ALE.
 *
 * \param[in, out]  domain       pointer to a cs_domain_t structure
 * \param[in, out]  bc_type      boundary face types
 * \param[in, out]  ale_bc_type  boundary face types for mesh velocity
 *                               (see \ref cs_boundary_ale_subtype_bits_t)
 * \param[in]       impale       indicator for prescribed node displacement
 *                               (0 or 1)
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_boundary_conditions_ale
void
cs_user_boundary_conditions_ale([[maybe_unused]] cs_domain_t  *domain,
                                [[maybe_unused]] int           bc_type[],
                                [[maybe_unused]] int           ale_bc_type[],
                                [[maybe_unused]] int           impale[])
{
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
