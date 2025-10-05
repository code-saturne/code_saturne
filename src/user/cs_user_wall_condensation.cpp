/*============================================================================
 * User source terms associated at the boundary faces and the neighboring
 * cells with surface condensation.
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
 * Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * Source terms associated at the boundary faces and the neighboring
 * cells with surface condensation.
 *
 * This function fills the condensation source terms for each variable at
 * the cell center associated to the boundary faces identifed in the mesh.
 * The fluid exchange coefficient is computed with a empirical law to be
 * imposed at the boundary face where the condensation phenomenon occurs.
 *
 * \param[in]  iappel     indicates which at which stage the routine is
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_wall_condensation
void
cs_user_wall_condensation([[maybe_unused]] int  iappel)
{
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
