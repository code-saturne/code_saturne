/*============================================================================
 * Unit test for some geometrical algorithms.
 *============================================================================*/

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_math.h"

/*---------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Function to check if a point is "inside" or "outside" a plane, according
 * to its normal. Basically, performs the dot product between the point and
 * the normal and check its sign.
 *
 * parameters:
 *   plane  <-- plane definition (x, y, z, nx, ny, nz)
 *   x      <-- the point considered
 *
 * returns:
 *   true if the point is "inside", false otherwise
  ----------------------------------------------------------------------------*/

static bool
_is_point_inside_plane(cs_real_t plane[6],
                       cs_real_t x[3])
{
  bool val = false;

  cs_real_t *p = &plane[0];
  cs_real_t *n = &plane[3];

  cs_real_t psca = cs_math_3_distance_dot_product(p, x, n);

  if (psca <= 0.0)
    val = true;

  return val;
}

/*----------------------------------------------------------------------------
 * Function that performs the intersection between a plane and a polygon.
 * It returns the resulting polygon at the "inner" side of the plane
 * according to its normal.
 *
 * parameters:
 *   nb_vertex     <--> number of vertices of the polygon
 *   vertex_coords <--> coordinates of the vertices (size : 3*nb_vertex)
 *   plane         <--  plane definition (point + normal)
 *
 ----------------------------------------------------------------------------*/

static void
_polygon_plane_intersection(int           *nb_vertex,
                            cs_real_3_t  **vertex_coord,
                            cs_real_t      plane[6])
{
  /* Initial number of vertices in the polygon */
  int n_vtx = *nb_vertex;
  cs_real_3_t *vtx = *vertex_coord;

  cs_real_3_t *new_vtx = NULL;
  BFT_MALLOC(new_vtx, n_vtx + 1, cs_real_3_t);
  int j = 0;

  /* Now we check which edge is intersected by the plane */
  for (int i = 0; i < n_vtx; i++) {
    cs_lnum_t v1 = i;
    cs_lnum_t v2 = (i+1) % n_vtx;

    cs_real_t xn = cs_math_3_distance_dot_product(vtx[v1], plane, plane+3);
    cs_real_t xd = cs_math_3_distance_dot_product(vtx[v1], vtx[v2], plane+3);

    // if the edge is //, we keep the vertex
    if (cs_math_fabs(xd) < 1.0e-12) {
      if (_is_point_inside_plane(plane, vtx[v2])) {
        assert(j <= n_vtx);
        for (cs_lnum_t dir = 0; dir < 3; dir ++)
          new_vtx[j][dir] = vtx[v2][dir];
        j++;
      }
    }

    // if the edge is not // to the plane
    else {
      cs_real_t t = xn/xd;

      // If intersection, new vertex
      if (t > 0 && t < 1.0) {
        assert(j <= n_vtx);
        for (cs_lnum_t dir = 0; dir < 3; dir ++)
          new_vtx[j][dir] = vtx[v1][dir] + t * (vtx[v2][dir] - vtx[v1][dir]);
        j++;
      }

      // We check if the second point is "inside" inside the plane
      // if yes, add the vertex to the new vertex list
      if (_is_point_inside_plane(plane, vtx[v2])) {
        assert(j <= n_vtx);
        for (cs_lnum_t dir = 0; dir < 3; dir ++)
          new_vtx[j][dir] = vtx[v2][dir];
        j++;
      }
    }
  }

  BFT_REALLOC(vtx, j, cs_real_3_t);

  for (cs_lnum_t i = 0; i < j; i++) {
    printf("%d: ", i);
    for (cs_lnum_t dir = 0; dir < 3; dir ++) {
      vtx[i][dir] = new_vtx[i][dir];
      printf(" %g", vtx[i][dir]);
    }
    printf("\n");
  }

  BFT_FREE(new_vtx);

  *nb_vertex = j;
  *vertex_coord = vtx;
}

/*---------------------------------------------------------------------------*/

int
main (int argc, char *argv[])
{
  CS_UNUSED(argc);
  CS_UNUSED(argv);

  int nb_vertex = 4;
  cs_real_3_t *vertex_coord;
  BFT_MALLOC(vertex_coord, 4, cs_real_3_t);

  vertex_coord[0][0] = 9;
  vertex_coord[0][1] = 1;
  vertex_coord[0][2] = 0;

  vertex_coord[1][0] = 10;
  vertex_coord[1][1] = 1;
  vertex_coord[1][2] = 0;

  vertex_coord[2][0] = 10;
  vertex_coord[2][1] = 2;
  vertex_coord[2][2] = 0;

  vertex_coord[3][0] = 9;
  vertex_coord[3][1] = 2;
  vertex_coord[3][2] = 1e-13;

  cs_real_t plane[6] = {10, 2, 0,
                        0, 0, 1};

  _polygon_plane_intersection(&nb_vertex,
                              &vertex_coord,
                              plane);

  exit(EXIT_SUCCESS);
}
