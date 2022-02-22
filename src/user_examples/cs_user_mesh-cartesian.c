/*============================================================================
 * Definition of the calculation mesh.
 *
 * Mesh-related user functions (called in this order):
 *   Manage the exchange of data between code_saturne and the pre-processor
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
 * \brief Define a cartesian mesh to use during computation.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_mesh_cartesian_define(void)
{

  /*! [mesh_cartesian_1] */
  {
    /* Define a mesh with constant step size in all directions.
     * Here we create a mesh with the following parameters :
     * NX = 10, Xmin = 0.,  Xmax = 1.
     * NY = 15, Ymin = -1., Ymax = 1.
     * NZ = 6,  Zmin = 0.5, Zmax = 2.5
     */
    /* Array of size 3 containing number of cells per direction */
    int nxyz[3] = {10, 15, 6};

    /* Array of size 6, containing min values and then max values for each
     * direction coordinate.
     */
    cs_real_t xyz[6] = {0., -1., 0.5, 1., 1., 2.5};

    cs_mesh_cartesian_define_simple(nxyz, xyz);
  }
  /*! [mesh_cartesian_1] */

  /*! [mesh_cartesian_2] */
  {
    /* Create a mesh using simple parameters */
    cs_mesh_cartesian_create();

    /* Constant step for X direction */
    cs_mesh_cartesian_define_dir_params(0, /* Direction index */
                                        CS_MESH_CARTESIAN_CONSTANT_LAW,
                                        10, /* Number of cells */
                                        0., /* Min value of direction coord */
                                        5., /* Max value of direction coord */
                                        -1.); /* Progression, not used for
                                                 constant step */

    /* Geometric step for Y direction,
     * arguments are the same of the Parabolic law:
     * CS_MESH_CARTESIAN_PARABOLIC
     */
    cs_mesh_cartesian_define_dir_params(1, /* Direction index */
                                        CS_MESH_CARTESIAN_GEOMETRIC_LAW,
                                        15, /* Number of cells */
                                        0., /* Min value of direction coord */
                                        2., /* Max value of direction coord */
                                        1.2); /* Progression */

    /* User parameters for Z direction.
     * Parameters to provide are the number of cells, and an array
     * of the vertices coordinates of size ncells+1.
     */
    int nz = 7;
    cs_real_t zvtx[8] = {0., 0.1, 0.3, 0.42, 0.8, 0.9, 1.2, 1.5};
    cs_mesh_cartesian_define_dir_user(2, nz, zvtx);

  }
  /*! [mesh_cartesian_2] */

  /*! [mesh_cartesian_3] */
  {
    /* Create a mesh based on a CSV file containing vertices coordinates
     * for each direction.
     *
     *  CSV file needs to contain :
     *  (1) First line which is empty or contains a header
     *  (2) Second line containing number of vertices per direction:
     *      format is 'nx;ny;nz'
     *  (3) Third line is empty or contains a header
     *  (4) Fourth line and onwards contains vertices coordinates for each
     *      direction. Format is "X1[i];X2[i];X3[i]" for index i.
     *      If current vertex index is beyond max value for a given
     *      direction, an empty value is expected.
     *      For example, if for index 'j' the first direction is empty,
     *      format is : ';X2[j];X3[j]'
     *
     */

    cs_mesh_cartesian_define_from_csv("cartesian_vertex_coordinates.csv");
  }
  /*! [mesh_cartesian_3] */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
