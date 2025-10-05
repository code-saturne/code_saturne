/*============================================================================
 * User function. Locally modify a given porosity to take into
 *  account erosion effect (for instance).
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

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * User function. Locally modify a given porosity to take into
 * account erosion effect (for instance).
 *
 *  This function is called for each time step.
 *  Porosity will be modified if
 *  cs_ibm->porosity_user_source_term_modification = true
 *
 * \param[in]   mesh               pointer to associated mesh structure
 * \param[in]   mesh_quantities    pointer to associated mesh quantities
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_ibm_modify
void
cs_user_ibm_modify
(
  [[maybe_unused]] const cs_mesh_t *mesh,
  [[maybe_unused]] const cs_mesh_quantities_t *mesh_quantities
)
{
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
