/*============================================================================
 * Definition of the calculation mesh.
 *
 * Mesh-related user functions (called in this order):
 *   Manage the exchange of data between code_saturne and the pre-processor
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
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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
 * Define a cartesian mesh.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_mesh_cartesian_define(void)
{
  /*
   * If a cartesian block
   * is created with a non-null name, then all groups will be prefixed with the
   * given name. For example, "X0" face group will be called "<name>_X0" instead.
   * In this example we create a three block mesh. */

  /*! [mesh_multiblock_cartesian_1] */
  {
    /* We define three blocks : two cubic rooms linked by a square duct.*/

    /* Rooms edge length */
    const cs_real_t room_edge_size = 1.0;

    /* duct dimensions */
    const cs_real_t duct_len = 1.;
    const cs_real_t duct_edge_size = 0.1;

    cs_real_t duct_xminmax[2];
    duct_xminmax[0] = - duct_len * 0.5;
    duct_xminmax[1] =   duct_len * 0.5;

    int nxyz_rooms[3] = {40, 40, 40};
    int nxyz_duct[3] = {40, 4, 4};

    /* Define Room1 coordinates*/
    cs_real_t xyz_room1[6];
    xyz_room1[0] = duct_xminmax[0] - room_edge_size;
    xyz_room1[1] = - room_edge_size * 0.5;
    xyz_room1[2] = 0.;
    xyz_room1[3] = duct_xminmax[0];
    xyz_room1[4] = + room_edge_size * 0.5;
    xyz_room1[5] = room_edge_size;

    /* Define duct coordinates*/
    cs_real_t xyz_duct[6];
    xyz_duct[0] = duct_xminmax[0];
    xyz_duct[1] = - duct_edge_size*0.5;
    xyz_duct[2] = room_edge_size - duct_edge_size;
    xyz_duct[3] = duct_xminmax[1];
    xyz_duct[4] = + duct_edge_size*0.5;
    xyz_duct[5] = room_edge_size;

    /* Define Room2 coordinates*/
    cs_real_t xyz_room2[6];
    xyz_room2[0] = duct_xminmax[1];
    xyz_room2[1] = - room_edge_size * 0.5;
    xyz_room2[2] = 0.;
    xyz_room2[3] = duct_xminmax[1] + room_edge_size;
    xyz_room2[4] = + room_edge_size * 0.5;
    xyz_room2[5] = room_edge_size;

    /* Create the 3 cartesian blocks */
    cs_mesh_cartesian_define_simple("room1", nxyz_rooms, xyz_room1);
    cs_mesh_cartesian_define_simple("duct", nxyz_duct, xyz_duct);
    cs_mesh_cartesian_define_simple("room2", nxyz_rooms, xyz_room2);
  }
  /*! [mesh_multiblock_cartesian_1] */

}

/*----------------------------------------------------------------------------*/
/*
 * Define mesh joinings.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_join(void)
{
  /*! [mesh_multiblock_cartesian_join] */

  /* Add a joining operations for the 3 cartesian blocks defined in
   * cs_user_mesh_cartesian.
   */

  /* Add a joining operation */
  /* ----------------------- */

  int    verbosity = 1;     /* per-task dump if > 1, debug level if >= 3 */
  int    visualization = 1; /* debug level if >= 3 */
  float  fraction = 0.10, plane = 25.;

  [[maybe_unused]]
  int join_num = cs_join_add("room1_X1 or duct_X0 or duct_X1 or room2_X0",
                             fraction,
                             plane,
                             verbosity,
                             visualization);
  /*! [mesh_multiblock_cartesian_join] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
