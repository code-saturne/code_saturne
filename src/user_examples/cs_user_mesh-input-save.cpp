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

/*----------------------------------------------------------------------------
 * Combine transformation matrixes (c = a.b)
 *
 * parameters:
 *   a <-- first transformation matrix
 *   b <-- second transformation matrix
 *   c --> combined transformation matrix
 *---------------------------------------------------------------------------*/

static void
_combine_tr_matrixes(const double  a[3][4],
                     const double  b[3][4],
                     double  c[3][4])
{
  c[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0];
  c[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1] + a[0][2]*b[2][1];
  c[0][2] = a[0][0]*b[0][2] + a[0][1]*b[1][2] + a[0][2]*b[2][2];
  c[0][3] = a[0][0]*b[0][3] + a[0][1]*b[1][3] + a[0][2]*b[2][3] + a[0][3];

  c[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0] + a[1][2]*b[2][0];
  c[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1] + a[1][2]*b[2][1];
  c[1][2] = a[1][0]*b[0][2] + a[1][1]*b[1][2] + a[1][2]*b[2][2];
  c[1][3] = a[1][0]*b[0][3] + a[1][1]*b[1][3] + a[1][2]*b[2][3] + a[1][3];
  
  c[2][0] = a[2][0]*b[0][0] + a[2][1]*b[1][0] + a[2][2]*b[2][0];
  c[2][1] = a[2][0]*b[0][1] + a[2][1]*b[1][1] + a[2][2]*b[2][1];
  c[2][2] = a[2][0]*b[0][2] + a[2][1]*b[1][2] + a[2][2]*b[2][2];
  c[2][3] = a[2][0]*b[0][3] + a[2][1]*b[1][3] + a[2][2]*b[2][3] + a[2][3];
}

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

    cs_preprocessor_data_add_file("mesh_input/mesh_01", 0, nullptr, nullptr);
  }
  /*! [mesh_input_1] */

  /*! [mesh_input_2] */
  {
    /* Add same mesh with transformations */
    const char *renames[] = {"Inlet", "Injection_2",
                             "Group_to_remove", nullptr};

    /* Transformation matrix in homogeneous coordinates */
    /* Last raw is not represented */
    double transf_matrix[3][4];

    /* Here a rotation of -90° in direction z */
    const double theta = -0.5 * cs_math_pi; /* radians */
    const double transf_matrix_rot[3][4] = {{cos(theta), -sin(theta), 0., 0.},
                                            {sin(theta),  cos(theta), 0., 0.},
                                            {        0.,          0., 1., 0.}};

    /* Here a translation of -4m in direction x and 6m in y */
    const double transf_matrix_trans[3][4] = {{1., 0., 0., -4.},
                                              {0., 1., 0.,  6.},
                                              {0., 0., 1.,  0.}};

    /* Combine transformation matrixes */
    _combine_tr_matrixes(transf_matrix_trans, transf_matrix_rot, transf_matrix);

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
