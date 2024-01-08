/*============================================================================
 * Definition of the calculation mesh.
 *
 * Mesh-related user functions (called in this order):
 *   Manage the exchange of data between code_saturne and the pre-processor
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
 * \brief Mesh input definition and mesh saving examples.
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
 * This behavior may be changed in the GUI:
 * - In the "Mesh" section, unchecking
 *   "Use unmodified checkpoint mesh in case of restart" sets the mode
 *   to CS_PREPROCESSOR_DATA_RESTART_AND_MODIFY.
 * - In "Time Settings - Start/Restart", selecting "different mesh" sets the
 *   mode to CS_PREPROCESSOR_DATA_RESTART_NONE.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_mesh_restart_mode(void)
{
  /*! [mesh_restart_1] */
  cs_preprocessor_data_set_restart_mode(CS_PREPROCESSOR_DATA_RESTART_NONE);
  /*! [mesh_restart_1] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define mesh files to read and optional associated transformations.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_mesh_input(void)
{
  /*! [mesh_input_1] */
  {
    /* Determine list of files to add */
    /*--------------------------------*/

    /* Read input mesh with no modification */

    cs_preprocessor_data_add_file("mesh_input/mesh_01", 0, NULL, NULL);
  }
  /*! [mesh_input_1] */

  /*! [mesh_input_2] */
  {
    /* Add same mesh with transformations:
     * Here a translation in direction x,
     * and a rotation of theta around axe z.
     * */
    const char *renames[] = {"Inlet", "Injection_2",
                             "Group_to_remove", NULL};
    const double  theta = 0.1; /* radians */
    const double transf_matrix[3][4] = {{ cos(theta), sin(theta), 0., 5.},
                                        {-sin(theta), cos(theta), 0., 0.},
                                        {         0.,         0., 1., 0.}};

    cs_preprocessor_data_add_file("mesh_input/mesh_02",
                                  2, renames,
                                  transf_matrix);
  }
  /*! [mesh_input_2] */

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

void
cs_user_mesh_save(cs_mesh_t  *mesh)
{
  /*! [mesh_save] */
  {
    /* Disable saving of modified mesh by setting flag to 0;
       Force saving of unmodified mesh by setting it to 2. */

    mesh->save_if_modified = 0;

  }
  /*! [mesh_save] */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
