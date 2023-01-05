/*============================================================================
 * User definitions for fluid-structure interaction using ALE.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
 * \file cs_user_fluid_structure_interaction.c
 *
 * \brief User-defined functions dedicated to Fluid-Structure interaction
 *        modeling.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define external structure ids for faces associated with external
 *        (code_aster) structures.
 *
 * Structure ids associated to a given face have the following values:
 * - -i where coupled to  i-th (1-to n) external (code_aster) structure.
 * - 0 where not coupled with an internal or external structure.
 * - i  where coupled to  i-th (1-to n) internal (mass-spring) structure.
 *
 * \param[in, out]  domain         pointer to a cs_domain_t structure
 * \param[in, out]  structure_id   structure id associated to each face
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_fsi_external_structure_id
void
cs_user_fsi_external_structure_id(cs_domain_t  *domain,
                                  int           structure_id[])
{
  CS_UNUSED(domain);
  CS_UNUSED(structure_id);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
