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
 * \file cs_user_mesh-joining.c
 *
 * \brief Mesh joining example.
 *
 * See \ref cs_user_mesh for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define mesh joinings.
 *
 * This is done by calling the cs_join_add() function for each joining
 * operation to add.
 *
 * The arguments to \ref cs_join_add are:
 * \param [in] sel_criteria boundary face selection criteria string
 * \param [in] fraction value of the fraction parameter;
 *                    the initial tolerance radius associated to each vertex
 *                    is equal to the lenght of the shortest incident edge,
 *                    multiplied by this fraction.
 * \param [in] plane value of the plane parameter;
 *                    when subdividing faces, 2 faces are considered
 *                    coplanar and may be joined if angle between their
 *                    normals (in degrees) does not exceed this parameter.
 * \param [in] verbosity level of verbosity required
 *
 * The function returns a number (1 to n) associated with the
 * new joining. This number may optionnally be used to assign advanced
 * parameters to the joining.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_join(void)
{

  /*! [mesh_add_simple_joining] */

  int    join_num;

  /* Add a joining operation */
  /* ----------------------- */

  int    verbosity = 1;     /* per-task dump if > 1, debug level if >= 3 */
  int    visualization = 1; /* debug level if >= 3 */
  float  fraction = 0.10, plane = 25.;

  join_num = cs_join_add("98 or 99",
                         fraction,
                         plane,
                         verbosity,
                         visualization);

  /*! [mesh_add_simple_joining] */

  /*--------------------------------------------------------------------------*/

  /* Example with advanced parameters;
     Advanced parameters may be modified to solve errors during the
     joining step or to get a better mesh quality. */

  /*! [mesh_add_advanced_joining] */
  cs_join_set_advanced_param(join_num,
                             1.00,      /* merge tolerance factor */
                             0.10,      /* pre-merge factor */
                             1,         /* tolerance computation mode */
                             1,         /* intersection computation mode */
                             500,       /* max. nb. of equivalence breaks */
                             100,       /* maximum number of sub-faces */
                             30,        /* tree max level */
                             25,        /* tree max boxes per node */
                             5.0,       /* tree max ratio */
                             2.0);      /* distribution tree max ratio */
  /*! [mesh_add_advanced_joining] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
