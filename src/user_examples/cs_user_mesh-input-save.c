/*============================================================================
 * Definition of the calculation mesh.
 *
 * Mesh-related user functions (called in this order):
 *   Manage the exchange of data between Code_Saturne and the pre-processor
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_mesh.c
 *
 * \brief Mesh input definition and mesh saving examples.
 *
 * See \subpage cs_user_mesh for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

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
    /* Add same mesh with transformations */
    const char *renames[] = {"Inlet", "Injection_2",
                             "Group_to_remove", NULL};
    const double transf_matrix[3][4] = {{1., 0., 0., 5.},
                                        {0., 1., 0., 0.},
                                        {0., 0., 1., 0.}};

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
    /* Mark mesh as not modified (0) to disable saving;
       Mark it as modified (> 0) to force saving */

    mesh->modified = 0;

  }
  /*! [mesh_save] */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
