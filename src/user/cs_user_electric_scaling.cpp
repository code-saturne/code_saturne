/*============================================================================
 * Define scaling parameter for electric model
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

#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_electric_scaling.cpp
 *
 * \brief Define scaling parameter for electric model.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Rescale all electro-magnetic physical fields
 *        (electric potential, current density and Joule effect).
 *
 * \param[in] mesh pointer to a cs_mesh_t structure
 * \param[in,out] mesh_quantities pointer to a cs_mesh_quantities_t structure
 * \param[in] dt pointer to a \ref cs_real_t
 *
 * These options allow defining the time step synchronization policy,
 * as well as a time step multiplier.
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_scaling_elec
void
cs_user_scaling_elec
(
  [[maybe_unused]] const cs_mesh_t             *mesh,
  [[maybe_unused]] const cs_mesh_quantities_t  *mesh_quantities,
  [[maybe_unused]] cs_real_t                   *dt
)
{
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
