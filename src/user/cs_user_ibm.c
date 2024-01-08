/*============================================================================
 * User function. Define immersed boundaries in time and space.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"
#include "cs_stl.h"
#include "cs_ibm.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_ibm.c
 *
 * \brief User function. Define immersed boundaries in time and space.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function in which the user defines the objects to model.
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_ibm_define_objects
void
cs_user_ibm_define_objects(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function to set global parameters for the immersed boundaries
 *         module.
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_ibm_parameters
void
cs_user_ibm_parameters(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function where to apply predefined transformations to med/stl
 *         based objects.
 *
 * \param[in]  t            time value for the current time step
 *
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_ibm_object_transformations
void
cs_user_ibm_object_transformations(const cs_real_t time)
{
  CS_UNUSED(time);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function which allows the definition of a 'porous' object.
 *
 * \param[in]  c_id         local cell number
 * \param[in]  xyz          x, y, z coordinates of the current position
 * \param[in]  t            time value for the current time step
 * \param[in]  num_object   num of fsi object (if fsi activated)
 *
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_ibm_solid_por
void
cs_user_ibm_solid_por(const cs_lnum_t    c_id,
                      const cs_real_3_t  xyz,
                      const cs_real_t    t,
                      const int          num_object)
{
  CS_UNUSED(c_id);
  CS_UNUSED(xyz);
  CS_UNUSED(t);
  CS_UNUSED(num_object);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
