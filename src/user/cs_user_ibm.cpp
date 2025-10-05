/*============================================================================
 * User function. Define immersed boundaries in time and space.
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
 * User function in which the user defines the objects to model.
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_ibm_define_objects
void
cs_user_ibm_define_objects(void)
{
}

/*----------------------------------------------------------------------------*/
/*
 * User function to set global parameters for the immersed boundaries module.
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_ibm_parameters
void
cs_user_ibm_parameters(void)
{
}

/*----------------------------------------------------------------------------*/
/*
 * User function where to apply predefined transformations to MED/STL
 * based objects.
 *
 * \param[in]  t  time value for the current time step
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_ibm_object_transformations
void
cs_user_ibm_object_transformations([[maybe_unused]] const cs_real_t time)
{
}

/*----------------------------------------------------------------------------*/
/*
 * User function which allows the definition of a 'porous' object.
 *
 * \param[in]  c_id         local cell number
 * \param[in]  xyz          x, y, z coordinates of the current position
 * \param[in]  t            time value for the current time step
 * \param[in]  num_object   num of fsi object (if fsi activated)
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_ibm_solid_por
void
cs_user_ibm_solid_por([[maybe_unused]] cs_lnum_t        c_id,
                      [[maybe_unused]] const cs_real_t  xyz[3],
                      [[maybe_unused]] cs_real_t        t,
                      [[maybe_unused]] int              num_object)
{
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
