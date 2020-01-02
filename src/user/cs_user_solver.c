/*============================================================================
 * User solver
 *============================================================================*/

/* VERS */

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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_solver.c
 *
 * \brief User solver setting and implementation
 *
 * See \subpage user_solver for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set user solver
 *
 * \return  1 if user solver is called, 0 otherwise
 */
/*----------------------------------------------------------------------------*/

int
cs_user_solver_set(void)
{
  return 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Main call to user solver
 *
 * \param[in] mesh pointer to a cs_mesh_t structure
 * \param[in,out] mesh_quantities pointer to a cs_mesh_quantities_t structure
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_user_solver(const cs_mesh_t             *mesh,
               const cs_mesh_quantities_t  *mesh_quantities)
{
  CS_UNUSED(mesh);
  CS_UNUSED(mesh_quantities);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
