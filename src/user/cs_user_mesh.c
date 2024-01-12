/*============================================================================
 * Definition of the calculation mesh.
 *
 * Mesh-related user functions (called in this order):
 *   1) Manage the exchange of data between code_saturne and the pre-processor
 *   2) Define (conforming or non-conforming) mesh joinings.
 *   3) Define (conforming or non-conforming) periodicity.
 *   4) Define thin walls.
 *   5) Modify the geometry and mesh.
 *   6) Smoothe the mesh.
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
#include <stdarg.h>
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
 * \file cs_user_mesh.c
 *
 * \brief Definition and modification of the calculation mesh.
 *
 * See \ref cs_user_mesh for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Force preprocessing behavior in case of restart.
 *
 * By default, in case of restart, if a "restart/mesh_input.csm" file
 * is present, it will be read and proprocessing will be skipped.
 *
 * This behavior may be changed in the GUI (in the "Mesh" section, unchecking
 * "Use unmodified checkpoint mesh in case of restart"), or by calling
 * \ref cs_preprocessor_data_set_restart_mode in this function.
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_mesh_restart_mode
void
cs_user_mesh_restart_mode(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define mesh files to read and optional associated transformations.
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_mesh_input
void
cs_user_mesh_input(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define mesh joinings.
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_join
void
cs_user_join(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define periodic faces.
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_periodicity
void
cs_user_periodicity(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set options for cutting of warped faces.
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_mesh_warping
void
cs_user_mesh_warping(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Insert boundaries into a mesh.
 *
 * \param[in,out] mesh  pointer to a cs_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_mesh_boundary
void
cs_user_mesh_boundary(cs_mesh_t  *mesh)
{
  CS_UNUSED(mesh);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Modify geometry and mesh.
 *
 * \param[in,out] mesh  pointer to a cs_mesh_t structure
*/
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_mesh_modify
void
cs_user_mesh_modify(cs_mesh_t  *mesh)
{
  CS_UNUSED(mesh);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Mesh smoothing.
 *
 * \param[in,out] mesh  pointer to a cs_mesh_t structure
*/
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_mesh_smoothe
void
cs_user_mesh_smoothe(cs_mesh_t  *mesh)
{
  CS_UNUSED(mesh);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Enable or disable mesh saving.
 *
 * By default, mesh is saved when modified.
 *
 * \param[in,out] mesh  pointer to a cs_mesh_t structure
*/
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_mesh_save
void
cs_user_mesh_save(cs_mesh_t  *mesh)
{
  CS_UNUSED(mesh);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Tag bad cells within the mesh based on user-defined geometric criteria.
 *
 * \param[in,out] mesh  pointer to a cs_mesh_t structure
 * \param[in,out] mesh_quantities pointer to a cs_mesh_quantities_t structure
*/
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_mesh_bad_cells_tag
void
cs_user_mesh_bad_cells_tag(cs_mesh_t             *mesh,
                           cs_mesh_quantities_t  *mesh_quantities)
{
  CS_UNUSED(mesh);
  CS_UNUSED(mesh_quantities);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply partial modifications to the mesh after the preprocessing
 *        stage, but before initial postprocessing mesh building.
 *
 * \param[in,out] mesh  pointer to a cs_mesh_t structure
 * \param[in,out] mesh_quantities pointer to a cs_mesh_quantities_t structure
*/
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_mesh_modify_partial
void
cs_user_mesh_modify_partial(cs_mesh_t             *mesh,
                            cs_mesh_quantities_t  *mesh_quantities)
{
  CS_UNUSED(mesh);
  CS_UNUSED(mesh_quantities);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a cartesian mesh.
*/
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_mesh_cartesian_define
void
cs_user_mesh_cartesian_define(void)
{

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
