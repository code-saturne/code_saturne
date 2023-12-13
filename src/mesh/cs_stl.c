/*============================================================================
 * STL reader and writer.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_writer.h"
#include "cs_post.h"

#include "fvm_nodal_append.h"
#include "fvm_neighborhood.h"

#include "cs_math.h"
#include "cs_mesh_headers.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_rotation.h"

/*----------------------------------------------------------------------------
 * Header of the current file
 *----------------------------------------------------------------------------*/

#include "cs_stl.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Static global variables
 *============================================================================*/

static cs_stl_mesh_info_t _stl_meshes = {NULL, 0, 0};

cs_stl_mesh_info_t *cs_glob_stl_meshes = &_stl_meshes;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Function to cast 4 uint8_t value into one single uint32_t value
 *
 * parameters:
 *   val <-- array of 4 uint8_t
 *
 * returns :
 *   the uint32_t corresponding value
  ----------------------------------------------------------------------------*/

static uint32_t
_cast32(uint8_t *val)
{
  uint32_t retval;

  /* We have here 4 8bits unsigned int
   * and we cast then in a single 32bit
   * unsigned int by progessive shifts */

  retval =   (uint32_t)val[0]
         + ( (uint32_t)val[1]<<8  )
         + ( (uint32_t)val[2]<<16 )
         + ( (uint32_t)val[3]<<24 );

  return retval;
}

/*----------------------------------------------------------------------------
 * Function to cut one uint32_t value into an array of 4 uint8_t value
 *
 * parameters:
 *   out --> cut array of 4 uint8_t
 *   val <-- input uint32_t value
  ----------------------------------------------------------------------------*/

static void
_cut32(uint8_t *out, int val)
{
  out[0] = (val       ) & 0xff;
  out[1] = (val >> 8  ) & 0xff;
  out[2] = (val >> 16 ) & 0xff;
  out[3] = (val >> 24 ) & 0xff;
}

/*----------------------------------------------------------------------------
 * Function to compute the normal of a triangle defined by 3 points
 *
 * parameters:
 *   coords <-- triangle coordinates
 *   normal --> triangle normal
  ----------------------------------------------------------------------------*/

static void
_compute_normal(cs_real_t  coords[3][3],
                cs_real_t  normal[3])
{
  /* Compute and Write the face normal coordinates */
  const cs_real_t *a = coords[0], *b = coords[1], *c = coords[2];

  // ab vect ac
  normal[0] = (b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1]);
  normal[1] = (b[2]-a[2])*(c[0]-a[0]) - (b[0]-a[0])*(c[2]-a[2]);
  normal[2] = (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);

  cs_real_t norm = cs_math_3_norm(normal);

  for (int dir = 0; dir < 3; dir ++)
    normal[dir] /= norm;
}

/*----------------------------------------------------------------------------
 * Function to check if an AABB (Axis Aligned Bounding Box) is overlaped
 * by a triangle. It is based on the SAT (Separating Axis Theorem).
 *
 * Warning : the user must first be sure that the bounding boxes of
 * the AABB and the triangle intersects, otherwise, this functions might
 * give wrong results.
 *
 * parameters:
 *   box_extents <-- coordinates of the bounding box (xmin,...,xmax,...)
 *   tria_coords <-- coordinates of the triangle vertices
 *
 * returns:
 *   a boolean
  ----------------------------------------------------------------------------*/

static bool
_triangle_box_intersect(const cs_real_t  box_extents[6],
                        const cs_real_t  tria_coords[3][3])
{

  // Origin axis
  cs_real_t e[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

  // Width of the box, center of the box
  cs_real_t h[3], c[3];
  h[0] = 0.5*(box_extents[3]-box_extents[0]);
  h[1] = 0.5*(box_extents[4]-box_extents[1]);
  h[2] = 0.5*(box_extents[5]-box_extents[2]);

  c[0] = box_extents[0] + h[0];
  c[1] = box_extents[1] + h[1];
  c[2] = box_extents[2] + h[2];

  // Move triangle so that the box is virtually centered on 0
  cs_real_t v[3][3];
  for (int vi = 0; vi < 3; vi++) {
    for (int dir = 0; dir < 3; dir++)
      v[vi][dir] = tria_coords[vi][dir] - c[dir];
  }

  // Compute triangle edges
  cs_real_t f[3][3];
  for (int dir = 0; dir < 3; dir++) {
    f[0][dir] = v[1][dir] - v[0][dir];
    f[1][dir] = v[2][dir] - v[1][dir];
    f[2][dir] = v[0][dir] - v[2][dir];
  }

  // Construction of the 9 testing axis resulting from
  // cross product combinations of the edges of the
  // triangle and the edges of the AABB on which vertices
  // will be projected.

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {

      cs_real_t axis[3];
      cs_math_3_cross_product(f[i], e[j], axis);

      // Projection of the triangle's vertices on the
      // seperated axis
      cs_real_t p[3];
      for (int dir = 0; dir < 3; dir ++)
        p[dir] = cs_math_3_dot_product(axis, v[dir]);

      // Projection of the box on the seperated axis
      // computation of the "radius" r
      cs_real_t r = h[0] * CS_ABS(axis[0])
                  + h[1] * CS_ABS(axis[1])
                  + h[2] * CS_ABS(axis[2]);

      cs_real_t pmin = fmin(p[0], fmin(p[1], p[2]));
      cs_real_t pmax = fmax(p[0], fmax(p[1], p[2]));

      // Finally, we test if most extreme
      // of the triangle points intersects r
      if (pmin > r || pmax < -r)
        return false;
    }
  }

  // Construction of the last testing axis resulting from
  // the triangle normal, repeat the process
  cs_real_t n[3];
  cs_math_3_cross_product(f[0], f[1], n);

  cs_real_t p[3];
  for (int dir = 0; dir < 3; dir ++)
    p[dir] = cs_math_3_dot_product(n, v[dir]);

  cs_real_t r = h[0] * CS_ABS(n[0])
              + h[1] * CS_ABS(n[1])
              + h[2] * CS_ABS(n[2]);

  cs_real_t pmin = fmin(p[0], fmin(p[1],p[2]));
  cs_real_t pmax = fmax(p[0], fmax(p[1],p[2]));

  if ( pmin > r || pmax < -r )
    return false;


  return true;
}

/*----------------------------------------------------------------------------
 * Function to check if a point is "inside" or "outside" a plane, according
 * to its normal. Basically, performs the dot product between the point and
 * the normal and check its sign.
 *
 * parameters:
 *   plane  <-- plane definition (x,y,z, Nx,Ny,Nz)
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

  cs_real_t psca = cs_math_3_distance_dot_product(p,x,n);

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
    for (cs_lnum_t dir = 0; dir < 3; dir ++)
      vtx[i][dir] = new_vtx[i][dir];
  }

  BFT_FREE(new_vtx);

  *nb_vertex = j;
  *vertex_coord = vtx;
}

/*----------------------------------------------------------------------------
 * Function to compute the approximate surface of a triangle overlapping
 * an AABB (Axis Aligned Bounding Box). It is based on a Monte-Carlo like
 * method that randomly draw points in the triangle.
 *
 * parameters:
 *   box_extents <-- coordinates of the bounding box (xmin,...,xmax,...)
 *   tria_coords <-- coordinates of the triangle vertices
 *
 * returns:
 *   the surface of the triangle inside the bounding box
  ----------------------------------------------------------------------------*/

static cs_real_t
_triangle_box_surface_intersection(const cs_real_t  box_extents[6],
                                   const cs_real_t  tria_coords[3][3])
{
  cs_real_t surf = 0.;

  const cs_real_t *a = tria_coords[0];
  const cs_real_t *b = tria_coords[1];
  const cs_real_t *c = tria_coords[2];

  /* Computation of the total triangle surface */
  cs_real_t ab[3], ac[3], vec[3];

  for (int i = 0; i < 3; i++) {
    ab[i] = b[i]-a[i];
    ac[i] = c[i]-a[i];
  }

  cs_math_3_cross_product(ab,ac,vec);
  surf = 0.5 * cs_math_3_norm(vec);

  /* Particular case if the triangle is totally included in the box */
  cs_real_t bbox[6];
  bbox[0] =  HUGE_VAL; bbox[1] =  HUGE_VAL; bbox[2] =  HUGE_VAL;
  bbox[3] = -HUGE_VAL; bbox[4] = -HUGE_VAL; bbox[5] = -HUGE_VAL;

  for (int vtx_id = 0; vtx_id < 3; vtx_id ++) {
    bbox[0] = CS_MIN(bbox[0], tria_coords[vtx_id][0]);
    bbox[3] = CS_MAX(bbox[3], tria_coords[vtx_id][0]);
    bbox[1] = CS_MIN(bbox[1], tria_coords[vtx_id][1]);
    bbox[4] = CS_MAX(bbox[4], tria_coords[vtx_id][1]);
    bbox[2] = CS_MIN(bbox[2], tria_coords[vtx_id][2]);
    bbox[5] = CS_MAX(bbox[5], tria_coords[vtx_id][2]);
  }

  if (   bbox[0] >= box_extents[0]
      && bbox[1] >= box_extents[1]
      && bbox[2] >= box_extents[2]
      && bbox[3] <= box_extents[3]
      && bbox[4] <= box_extents[4]
      && bbox[5] <= box_extents[5])
    return surf;

  /* Approximation of the surface enclosed in the AABB */

  int max_draw = 5000;
  int nb_point_in = 0;
  cs_real_t px[3];

  for (int nb_draw = 0; nb_draw < max_draw; nb_draw ++) {
    /* Random draw in [[0:1]] */
    cs_real_t ra = (cs_real_t)rand() / (cs_real_t)RAND_MAX;
    cs_real_t rb = (cs_real_t)rand() / (cs_real_t)RAND_MAX;

    if (ra+rb >= 1.0){
      ra = 1.0-ra;
      rb = 1.0-rb;
    }

    /* Contruction of the random point */
    for (int i = 0; i < 3; i++)
      px[i] = a[i] + ra*(b[i]-a[i]) + rb*(c[i]-a[i]);

    if (   px[0] >= box_extents[0]
        && px[1] >= box_extents[1]
        && px[2] >= box_extents[2]
        && px[0] <= box_extents[3]
        && px[1] <= box_extents[4]
        && px[2] <= box_extents[5])
      nb_point_in ++;
  }

  surf *= (cs_real_t)nb_point_in / (cs_real_t)max_draw;

  return surf;
}

/*----------------------------------------------------------------------------
 * Function to compute the exact surface of a triangle overlapping
 * an AABB (Axis Aligned Bounding Box). It is based on recursive cut
 * of the triangle by the planes defined by the AABB faces.
 *
 * parameters:
 *   box_extents <-- coordinates of the bounding box (xmin,...,xmax,...)
 *   tria_coords <-- coordinates of the triangle vertices
 *
 * returns:
 *   the surface of the triangle inside the bounding box
  ----------------------------------------------------------------------------*/

static cs_real_t
_exact_triangle_box_surface_intersection(const cs_real_t  box_extents[6],
                                         const cs_real_t  tria_coords[3][3])
{
  cs_real_t surf = 0.;

  const cs_real_t *a = tria_coords[0];
  const cs_real_t *b = tria_coords[1];
  const cs_real_t *c = tria_coords[2];

  /* Computation of the total triangle surface */
  cs_real_t ab[3], ac[3], vec[3];

  for (int i = 0; i < 3; i++) {
    ab[i] = b[i]-a[i];
    ac[i] = c[i]-a[i];
  }

  cs_math_3_cross_product(ab,ac,vec);
  surf = 0.5 * cs_math_3_norm(vec);

  /* Particular case if the triangle is totally included in the box */
  cs_real_t bbox[6];
  bbox[0] =  HUGE_VAL; bbox[1] =  HUGE_VAL; bbox[2] =  HUGE_VAL;
  bbox[3] = -HUGE_VAL; bbox[4] = -HUGE_VAL; bbox[5] = -HUGE_VAL;

  for (int vtx_id = 0; vtx_id < 3; vtx_id ++) {
    bbox[0] = CS_MIN(bbox[0], tria_coords[vtx_id][0]);
    bbox[3] = CS_MAX(bbox[3], tria_coords[vtx_id][0]);
    bbox[1] = CS_MIN(bbox[1], tria_coords[vtx_id][1]);
    bbox[4] = CS_MAX(bbox[4], tria_coords[vtx_id][1]);
    bbox[2] = CS_MIN(bbox[2], tria_coords[vtx_id][2]);
    bbox[5] = CS_MAX(bbox[5], tria_coords[vtx_id][2]);
  }

  if (   bbox[0] >= box_extents[0]
      && bbox[1] >= box_extents[1]
      && bbox[2] >= box_extents[2]
      && bbox[3] <= box_extents[3]
      && bbox[4] <= box_extents[4]
      && bbox[5] <= box_extents[5])
    return surf;

  /* Successively cut the triangle by the planes
   * defined byt the box faces */
  int nv = 3;
  cs_real_3_t *coords = NULL ;
  BFT_MALLOC(coords, nv, cs_real_3_t);

  /* Polygon init */
  for (int i = 0; i < nv; i ++){
    for (int dir = 0; dir < 3; dir ++)
      coords[i][dir] = tria_coords[i][dir];
  }

  /* Plane definition */
  cs_real_t plane[6][6];
  cs_real_t e[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      plane[i][j]   = box_extents[j];
      plane[i][j+3] = -e[i][j];
    }
  }

  for (int i = 3; i < 6; i++){
    for (int j = 0; j < 3; j++){
      plane[i][j]   = box_extents[j+3];
      plane[i][j+3] = e[i-3][j];
    }
  }

  /* Recursively cut the tiangle by the planes */
  for (int f = 0; f < 6; f++)
    _polygon_plane_intersection(&nv, &coords, plane[f]);

  /* If no intersection, surface is 0 */
  if (nv ==0) {
    surf = 0;
    return surf;
  }

  /* Compute approximate face center coordinates for the polygon */
  cs_real_t a_center[3] = {0, 0, 0};
  cs_real_t f_norm[3] = {0, 0, 0};

  for (int i = 0; i < nv; i++) {
    for (int dir = 0; dir < 3; dir++)
      a_center[dir] += coords[i][dir];
  }

  for (cs_lnum_t i = 0; i < 3; i++)
    a_center[i] /= nv;

  /* Surface computation */
  cs_real_t vc0[3], vc1[3], vn[3];

  for (cs_lnum_t i = 0; i < nv; i++) {

    const int v0 = i;
    const int v1 = (i+1) % nv;

    for (cs_lnum_t dir = 0; dir < 3; dir++) {
      vc0[dir] = coords[v0][dir] - a_center[dir];
      vc1[dir] = coords[v1][dir] - a_center[dir];
    }

    cs_math_3_cross_product(vc0, vc1, vn);

    for (cs_lnum_t dir = 0; dir < 3; dir++)
      f_norm[dir] += 0.5*vn[dir];
  }

  surf = cs_math_3_norm(f_norm);

  BFT_FREE(coords);

  return surf;
}

/*----------------------------------------------------------------------------
 * Function to compute the volume of a tetrahedron
 *
 * parameters:
 *   x1      <-- coordinates of the tetrahedra (x1,y1,z1)
 *   x2      <-- coordinates of the tetrahedra (x2,y2,z2)
 *   x3      <-- coordinates of the tetrahedra (x3,y3,z3)
 *   x4      <-- coordinates of the tetrahedra (x4,y4,z4)
 *
 * returns:
 *   the volume of the tetrahedron
  ----------------------------------------------------------------------------*/

static cs_real_t
_tetra_vol(cs_real_t x1[3],
           cs_real_t x2[3],
           cs_real_t x3[3],
           cs_real_t x4[3])
{
  cs_real_t tetra_vol;

  tetra_vol =  cs_math_fabs((  (x1[0]-x4[0])*(x2[1]-x4[1])
                             - (x1[1]-x4[1])*(x2[0]-x4[0])) * (x3[2]-x4[2])
                         + (   (x1[1]-x4[1])*(x2[2]-x4[2])
                             - (x1[2]-x4[2])*(x2[1]-x4[1])) * (x3[0]-x4[0])
                         + (   (x1[2]-x4[2])*(x2[0]-x4[0])
                             - (x1[0]-x4[0])*(x2[2]-x4[2])) * (x3[1]-x4[1]));

  tetra_vol *= cs_math_1ov6;

  return tetra_vol;
}

/*----------------------------------------------------------------------------
 * Function to compute the volume of a pyramid, whose base is defined by
 * x1, x2, x3 and x4.
 *
 * parameters:
 *   x1      <-- coordinates of the pyramid (x1,y1,z1)
 *   x2      <-- coordinates of the pyramid (x2,y2,z2)
 *   x3      <-- coordinates of the pyramid (x3,y3,z3)
 *   x4      <-- coordinates of the pyramid (x4,y4,z4)
 *   x5      <-- coordinates of the pyramid (x5,y5,z5)
 *
 * returns:
 *   the volume of the pyramid
  ----------------------------------------------------------------------------*/

static cs_real_t
_pyramid_vol(cs_real_t x1[3],
             cs_real_t x2[3],
             cs_real_t x3[3],
             cs_real_t x4[3],
             cs_real_t x5[3])
{
  cs_real_t xc[3];

  for (int i = 0; i < 3; i++)
    xc[i] = 0.25 * (x1[i]+x2[i]+x3[i]+x4[i]);

  cs_real_t vol12 = _tetra_vol(x1, x2, xc, x5);
  cs_real_t vol23 = _tetra_vol(x2, x3, xc, x5);
  cs_real_t vol34 = _tetra_vol(x3, x4, xc, x5);
  cs_real_t vol41 = _tetra_vol(x4, x1, xc, x5);

  cs_real_t pyram_vol = vol12 + vol23 + vol34 + vol41;

  return pyram_vol;
}

/*----------------------------------------------------------------------------
 * Function to compute the volume of a prism, whose base is defined by
 * x1, x2, x3 and x4.
 *
 * parameters:
 *   x1      <-- coordinates of the prism (x1,y1,z1)
 *   x2      <-- coordinates of the prism (x2,y2,z2)
 *   x3      <-- coordinates of the prism (x3,y3,z3)
 *   x4      <-- coordinates of the prism (x4,y4,z4)
 *   x5      <-- coordinates of the prism (x5,y5,z5)
 *   x6      <-- coordinates of the prism (x5,y5,z5)
 *
 * returns:
 *   the volume of the prism
  ----------------------------------------------------------------------------*/

static cs_real_t
_prism_vol(cs_real_t x1[3],
           cs_real_t x2[3],
           cs_real_t x3[3],
           cs_real_t x4[3],
           cs_real_t x5[3],
           cs_real_t x6[3])
{
  cs_real_t xc[3];

  for (int i = 0; i < 3; i++)
    xc[i] = (x1[i]+x2[i]+x3[i]+x4[i]+x5[i]+x6[i]) * cs_math_1ov6;

  cs_real_t vol136c  = _tetra_vol(x1, x3, x6, xc);
  cs_real_t vol245c  = _tetra_vol(x2, x4, x5, xc);
  cs_real_t vol1256c = _pyramid_vol(x1, x2, x5, x6, xc);
  cs_real_t vol1243c = _pyramid_vol(x1, x2, x4, x3, xc);
  cs_real_t vol3456c = _pyramid_vol(x3, x4, x5, x6, xc);

  cs_real_t prisme_vol = vol136c + vol245c + vol1256c +  vol1243c + vol3456c;

  return prisme_vol;
}

/*----------------------------------------------------------------------------
 * Function to compute the "outside" length of a segment cut by a plan,
 * according to its normal
 *
 * parameters:
 *   x1      <-- coordinates of the tetrahedra (x1,y1,z1)
 *   x2      <-- coordinates of the tetrahedra (x2,y2,z2)
 *   plane   <-- coordinates of the triangle vertices and normals
 *               (x1,y1,z1,nx,ny,nz)
 *
 * returns:
 *   the volume of the tetrahedra under the plane (according to its normal)
  ----------------------------------------------------------------------------*/

static cs_real_t
_plane_segment_intersection(cs_real_t x1[3],
                            cs_real_t x2[3],
                            cs_real_t plane[6])
{
  cs_real_t length = cs_math_3_distance(x1, x2);
  cs_real_t *a = &plane[0];
  cs_real_t *n = &plane[3];

  int v1_is_inside = 0; int v2_is_inside = 0;

  if (_is_point_inside_plane(plane, x1))
    v1_is_inside = 1;

  if (_is_point_inside_plane(plane, x2))
    v2_is_inside = 1;


  if (v1_is_inside == 0 && v2_is_inside == 0)
    length = 0.0;

  if (v1_is_inside == 1 && v2_is_inside == 0) {

    cs_real_t den = cs_math_3_distance_dot_product(x2, x1, n);
    cs_real_t num = cs_math_3_distance_dot_product(x2, a, n);

    cs_real_t lambda = 1.0;
    if (cs_math_fabs(den) > 1.e-20)
      lambda = num / den;

    length *= lambda;
  }

  if (v1_is_inside == 0 && v2_is_inside == 1) {

    cs_real_t den = cs_math_3_distance_dot_product(x1, x2, n);
    cs_real_t num = cs_math_3_distance_dot_product(x1, a, n);

    cs_real_t lambda = 1.0;
    if (cs_math_fabs(den) > 1.e-20)
      lambda = num / den;

    length *= lambda;
  }

  return length;
}

/*----------------------------------------------------------------------------
 * Function to compute the volume of a tetrahedra cut by a plane surface
 *
 * parameters:
 *   x1      <-- coordinates of the tetrahedra (x1,y1,z1)
 *   x2      <-- coordinates of the tetrahedra (x2,y2,z2)
 *   x3      <-- coordinates of the tetrahedra (x3,y3,z3)
 *   x4      <-- coordinates of the tetrahedra (x4,y4,z4)
 *   plane   <-- coordinates of the triangle vertices and normals
 *               (x1,y1,z1,nx,ny,nz)
 *
 * returns:
 *   the volume of the tetrahedra under the plane (according to its normal)
  ----------------------------------------------------------------------------*/

static cs_real_t
_tetrahedron_plane_volume_intersection(cs_real_t x1[3],
                                       cs_real_t x2[3],
                                       cs_real_t x3[3],
                                       cs_real_t x4[3],
                                       cs_real_t plane[6])
{
  cs_real_t vol = 0.0;

  int v1_is_inside = 0; int v2_is_inside = 0;
  int v3_is_inside = 0; int v4_is_inside = 0;

  /* Test if the tetra vertices are "inside" or "outside" the plane,
   * according to its normal */

  if (_is_point_inside_plane(plane, x1))
    v1_is_inside = 1;

  if (_is_point_inside_plane(plane, x2))
    v2_is_inside = 1;

  if (_is_point_inside_plane(plane, x3))
    v3_is_inside = 1;

  if (_is_point_inside_plane(plane, x4))
    v4_is_inside = 1;

  int nb_point_inside_plane = v1_is_inside + v2_is_inside + v3_is_inside + v4_is_inside;

  switch(nb_point_inside_plane) {

    case 0: /* entire tetra "outside" the plane */

      vol = _tetra_vol(x1, x2, x3, x4);

      break;

    case 1: /* 3 intersection points with one point inside */
    {
      vol = _tetra_vol(x1, x2, x3, x4);

      /* Looking for the point located inside */
      cs_real_t *a = NULL, *b = NULL, *c = NULL, *d = NULL;

      if (v1_is_inside == 1) {
        a = &x1[0];
        b = &x2[0];
        c = &x3[0];
        d = &x4[0];
      }
      else if (v2_is_inside == 1) {
        a = &x2[0];
        b = &x1[0];
        c = &x3[0];
        d = &x4[0];
      }
      else if (v3_is_inside == 1) {
        a = &x3[0];
        b = &x1[0];
        c = &x2[0];
        d = &x4[0];
      }
      else if (v4_is_inside == 1) {
        a = &x4[0];
        b = &x1[0];
        c = &x2[0];
        d = &x3[0];
      }

      assert(a != NULL && b != NULL && c != NULL && d != NULL);

      cs_real_t ab[3], ac[3], ad[3];
      for (int i = 0; i < 3; i++) {
        ab[i] = a[i] - b[i];
        ac[i] = a[i] - c[i];
        ad[i] = a[i] - d[i];
      }

      cs_real_t lab = cs_math_3_norm(ab);
      cs_real_t lac = cs_math_3_norm(ac);
      cs_real_t lad = cs_math_3_norm(ad);

      cs_real_t lab_inside =  lab - _plane_segment_intersection(a, b, plane);
      cs_real_t lac_inside =  lac - _plane_segment_intersection(a, c, plane);
      cs_real_t lad_inside =  lad - _plane_segment_intersection(a, d, plane);

      cs_real_t lab_frac = lab_inside / lab;
      cs_real_t lac_frac = lac_inside / lac;
      cs_real_t lad_frac = lad_inside / lad;

      vol *= (1.0 - lab_frac * lac_frac * lad_frac);

      break;
    }

    case 3: /* 3 intersection points with 3 points inside*/
    {
      vol = _tetra_vol(x1, x2, x3, x4);

      /* Lookink for the point located outside */
      cs_real_t *a = NULL, *b = NULL, *c = NULL, *d = NULL;

      if (v1_is_inside == 0){
        a = &x1[0];
        b = &x2[0];
        c = &x3[0];
        d = &x4[0];
      }
      else if (v2_is_inside == 0) {
        a = &x2[0];
        b = &x1[0];
        c = &x3[0];
        d = &x4[0];
      }
      else if (v3_is_inside == 0) {
        a = &x3[0];
        b = &x1[0];
        c = &x2[0];
        d = &x4[0];
      }
      else if (v4_is_inside == 0) {
        a = &x4[0];
        b = &x1[0];
        c = &x2[0];
        d = &x3[0];
      }

      assert(a != NULL && b != NULL && c != NULL && d != NULL);

      cs_real_t ab[3], ac[3], ad[3];
      for (int i = 0; i < 3; i++) {
        ab[i] = a[i] - b[i];
        ac[i] = a[i] - c[i];
        ad[i] = a[i] - d[i];
      }

      cs_real_t lab = cs_math_3_norm(ab);
      cs_real_t lac = cs_math_3_norm(ac);
      cs_real_t lad = cs_math_3_norm(ad);

      cs_real_t lab_inside =  _plane_segment_intersection(a, b, plane);
      cs_real_t lac_inside =  _plane_segment_intersection(a, c, plane);
      cs_real_t lad_inside =  _plane_segment_intersection(a, d, plane);

      cs_real_t lab_frac = lab_inside / lab;
      cs_real_t lac_frac = lac_inside / lac;
      cs_real_t lad_frac = lad_inside / lad;

      vol *= lab_frac * lac_frac * lad_frac;

      break;
    }

    case 2: /* 4 intersection points */
    {
      /* Looking for the 2 inside and 2 outside points */
      cs_real_t *a = NULL, *b = NULL, *c = NULL, *d = NULL;

      if (v1_is_inside && v2_is_inside) {
        a = &x1[0];
        b = &x2[0];
        c = &x3[0];
        d = &x4[0];
      }
      else if (v1_is_inside && v3_is_inside) {
        a = &x1[0];
        b = &x3[0];
        c = &x2[0];
        d = &x4[0];
      }
      else if (v1_is_inside && v4_is_inside) {
        a = &x1[0];
        b = &x4[0];
        c = &x2[0];
        d = &x3[0];
      }
      else if (v2_is_inside && v3_is_inside) {
        a = &x2[0];
        b = &x3[0];
        c = &x1[0];
        d = &x4[0];
      }
      else if (v2_is_inside && v4_is_inside) {
        a = &x2[0];
        b = &x4[0];
        c = &x1[0];
        d = &x3[0];
      }
      else if (v3_is_inside && v4_is_inside) {
        a = &x3[0];
        b = &x4[0];
        c = &x1[0];
        d = &x2[0];
      }

      assert(a != NULL && b != NULL && c != NULL && d != NULL);

      cs_real_t ac[3], ad[3], bc[3], bd[3];
      for (int i = 0; i < 3; i++) {
        ac[i] = c[i] - a[i];
        ad[i] = d[i] - a[i];
        bc[i] = c[i] - b[i];
        bd[i] = d[i] - b[i];
      }

      cs_real_t lac = cs_math_3_norm(ac);
      cs_real_t lad = cs_math_3_norm(ad);
      cs_real_t lbc = cs_math_3_norm(bc);
      cs_real_t lbd = cs_math_3_norm(bd);

      cs_real_t lac_inside = _plane_segment_intersection(a, c, plane);
      cs_real_t lad_inside = _plane_segment_intersection(a, d, plane);
      cs_real_t lbc_inside = _plane_segment_intersection(b, c, plane);
      cs_real_t lbd_inside = _plane_segment_intersection(b, d, plane);

      cs_real_t lac_frac = lac_inside / lac;
      cs_real_t lad_frac = lad_inside / lad;
      cs_real_t lbc_frac = lbc_inside / lbc;
      cs_real_t lbd_frac = lbd_inside / lbd;

      cs_real_t ac_inter[3], ad_inter[3], bc_inter[3], bd_inter[3];
      /* Coordinates of the intersection point between the tetra and the plane */
      for (int i = 0; i < 3; i++) {
        ac_inter[i] = a[i] + lac_frac * ac[i];
        ad_inter[i] = a[i] + lad_frac * ad[i];
        bc_inter[i] = b[i] + lbc_frac * bc[i];
        bd_inter[i] = b[i] + lbd_frac * bd[i];
      }

      vol = _prism_vol(a, b, ac_inter, bc_inter, bd_inter, ad_inter);

      break;
    }

    case 4:
      vol = 0.0;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                "Error in function _tetrahedra_plane_volume_intersection\n");
      break;
  }

  return vol;

}
/*----------------------------------------------------------------------------
 * Function to update the copy of the mesh at initialization
 * (necessary is the STL mesh needs to move with time)
 *
 * parameters:
 *   stl_mesh      <-- a STL mesh structure
  ----------------------------------------------------------------------------*/

static void
_update_init_mesh(cs_stl_mesh_t *stl_mesh)
{
  cs_lnum_t n_tria = stl_mesh->n_faces;
  cs_lnum_t n_points = n_tria*3;

  for (cs_lnum_t i = 0; i < n_points; i++) {
    for (int j = 0; j < 3; j++) {
      stl_mesh->coords_ini[i][j] = stl_mesh->coords[i][j];
    }
  }

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a new STL mesh structure.
 *
 * \param[in] name name of the STL mesh
 *
 * \return pointer to the new STL mesh structure
 */
/*----------------------------------------------------------------------------*/

cs_stl_mesh_t *
cs_stl_mesh_add(const char  *name)
{
  /* First check if it already exists */
  cs_stl_mesh_t  *stl_mesh = cs_stl_mesh_get_by_name(name);

  if (stl_mesh != NULL) {
    // The STL mesh already exists
    bft_error(__FILE__, __LINE__, 0,
              _("Error creating stl mesh: mesh %s already exists."), name);

  }
  else {
    // If it does not exists create it
    _stl_meshes.n_meshes++;
    BFT_REALLOC(_stl_meshes.mesh_list, _stl_meshes.n_meshes, cs_stl_mesh_t *);

    BFT_MALLOC(stl_mesh, 1, cs_stl_mesh_t);

    if (name != NULL) {
      strncpy(stl_mesh->name, name, 19);
      stl_mesh->name[19] = '\0';
    }
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Error creating stl mesh: no name given."));

    memset(stl_mesh->header, 0, 80);
    stl_mesh->n_faces = 0;
    stl_mesh->coords = NULL;
    stl_mesh->n_seeds = 0;
    stl_mesh->seed_coords = NULL;
    stl_mesh->is_porous = false;
    stl_mesh->ext_mesh = NULL;

    _stl_meshes.mesh_list[_stl_meshes.n_meshes - 1] = stl_mesh;
  }

  return stl_mesh;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a STL mesh based on its name if present.
 *
 * \param[in] name name of the STL mesh
 *
 * If no STL mesh of the given name is defined, NULL is returned.
 *
 * \return pointer to the STL mesh structure, or NULL
 */
/*----------------------------------------------------------------------------*/

cs_stl_mesh_t *
cs_stl_mesh_get_by_name(const char  *name)
{
  cs_stl_mesh_t *ptr = NULL;

  for (int s_id = 0; s_id < _stl_meshes.n_meshes; s_id ++) {
    cs_stl_mesh_t *stl_mesh = _stl_meshes.mesh_list[s_id];
    int test = strcmp(stl_mesh->name, name);
    if (test == 0)
      ptr = stl_mesh;
  }

  return ptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all allocated STL mesh structures
 */
/*----------------------------------------------------------------------------*/

void
cs_stl_mesh_destroy_all(void)
{
  for (int s_id = 0; s_id < _stl_meshes.n_meshes; s_id ++) {
    cs_stl_mesh_t *ptr = _stl_meshes.mesh_list[s_id];
    BFT_FREE(ptr->coords);
    BFT_FREE(ptr->coords_ini);
    BFT_FREE(ptr->seed_coords);
  }

  BFT_FREE(_stl_meshes.mesh_list);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read a binary STL file and store its content in a STL mesh structure
 *
 * Each STL file composed of the following header:
 *
 * uint8[80] – Header
 * uint32 – Number of triangles
 *
 * followed by 50 byte blocks for each triangle:
 *
 * - real32[3] – Normal vector
 * - real32[3] – Vertex 1 coordinates
 * - real32[3] – Vertex 2 coordinates
 * - real32[3] – Vertex 3 coordiantes
 * - uint16    – Attribute (any other information, usually void)
 *
 * \param[in] path     path to the STL file
 * \param[in] stl_mesh pointer to the associated STL mesh structure
 */
/*----------------------------------------------------------------------------*/

void
cs_stl_file_read(cs_stl_mesh_t  *stl_mesh,
                 const char     *path)
{
  uint8_t buf[128];
  FILE *fp;
  float *loc_coords = NULL;

  cs_lnum_t n_tria = 0;
  cs_lnum_t n_tria_new = 0;

  if (cs_glob_rank_id < 1) {

    /* Open STL file */
    fp = fopen(path, "rb");
    if (fp == NULL) {
      bft_error(__FILE__, __LINE__, 0,
                _("Error opening file \"%s\":\n\n"
                  "  %s"), path, strerror(errno));
    }

    /* Check if the file is ASCII or Binary */
    char temp[6] ;
    fread(temp, 5, 1, fp);
    temp[5] = '\0';
    int test = strcmp(temp, "solid");

    if (test == 0)
      bft_error(__FILE__, __LINE__, 0,
                _("Error opening file \"%s\":\n"
                  "You are probably tyring to open an ASCII STL file\n"
                  "Please convert your file to binary :-)"), path);


    /* Read and copy header */
    rewind(fp);
    fread(buf, 84, 1, fp);
    memcpy(stl_mesh->header, buf, 80);

    /* Read number of triangles */
    uint32_t ntri;
    ntri = _cast32(buf+80);
    n_tria = (int)ntri;

    stl_mesh->n_faces = n_tria;

    /* Allocation*/
    BFT_MALLOC(stl_mesh->coords     , 3*n_tria, cs_real_3_t);
    BFT_MALLOC(loc_coords , 9, float);

    /* Loop on triangle faces
       ---------------------- */
    for (uint32_t i = 0; i < ntri; i++) {

      // Each triangle data are contained in 50bytes blocks
      fread(buf, 50, 1, fp);
      uint8_t *start = buf + 12;

      // Read the 3 vertices for the current triangle
      for (uint32_t vi = 0; vi < 3; vi ++) {

        // Read the 3 coordinates for each vertex
        for (cs_lnum_t dir = 0; dir < 3; dir ++) {

          uint32_t v_temp = _cast32(start + 3*4*vi + 4*dir);
          float *dir_coo = (float *)(&v_temp);
          loc_coords[3*vi + dir] = *dir_coo;
        }
      }

      // Check if the current triangle has a null area

      cs_real_t a[3], b[3], c[3], n[3];

      for (int dir = 0; dir < 3; dir ++) {
        a[dir] = (cs_real_t)loc_coords[3*0 + dir];
        b[dir] = (cs_real_t)loc_coords[3*1 + dir];
        c[dir] = (cs_real_t)loc_coords[3*2 + dir];
      }

      n[0] = (b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1]);
      n[1] = (b[2]-a[2])*(c[0]-a[0]) - (b[0]-a[0])*(c[2]-a[2]);
      n[2] = (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);

      cs_real_t nn = cs_math_3_norm(n);

      if (nn > 1.e-20) {
        for (cs_lnum_t dir = 0; dir < 3; dir ++) {
          for (cs_lnum_t vi = 0; vi < 3; vi ++) {
            stl_mesh->coords[3*n_tria_new + vi][dir]
              = (cs_real_t)loc_coords[3*vi + dir];
          }
        }
        n_tria_new ++;
      }

    }

    /* Some log information */
    bft_printf(_("\n"
                 " ** Reading of STL mesh \"%s\"\n"
                 "    Number of triangles: %d\n"
                 "    Number of removed triangles: %d\n\n"),
               stl_mesh->name, (int)n_tria_new, (int)n_tria-n_tria_new);

    n_tria = n_tria_new;
    stl_mesh->n_faces = n_tria;

    /* Re-allocation*/
    BFT_REALLOC(stl_mesh->coords    , 3*n_tria, cs_real_3_t);

    /* Copy coordinates to a work aray that
     * will contain all the init coordinates */
    BFT_MALLOC(stl_mesh->coords_ini , 3*n_tria, cs_real_3_t);
    for (int i = 0; i < 3*n_tria; i++)
      for (int j = 0; j < 3; j++)
        stl_mesh->coords_ini[i][j] = stl_mesh->coords[i][j];

    BFT_FREE(loc_coords);
    fclose(fp);

  }

  /* Broadcast to other ranks */

  cs_parall_bcast(0, /* root_rank */
                  1,
                  CS_LNUM_TYPE,
                  &(stl_mesh->n_faces));

  if (cs_glob_rank_id > 0) {
    BFT_MALLOC(stl_mesh->coords    , stl_mesh->n_faces*3, cs_real_3_t);
    BFT_MALLOC(stl_mesh->coords_ini, stl_mesh->n_faces*3, cs_real_3_t);
  }

  cs_parall_bcast(0, /* root_rank */
                  stl_mesh->n_faces*9,
                  CS_REAL_TYPE,
                  stl_mesh->coords);

  cs_parall_bcast(0, /* root_rank */
                  stl_mesh->n_faces*9,
                  CS_REAL_TYPE,
                  stl_mesh->coords_ini);


  /* Merge identical vertices
     ------------------------ */

#if 0
  cs_lnum_t n_coo_ini = n_tria*3;

  cs_coord_t *_face_vtx_coord = NULL;
  BFT_MALLOC(_face_vtx_coord, n_coo_ini*3, cs_coord_t);
  for (cs_lnum_t i = 0; i < n_coo_ini*3; i++)
    _face_vtx_coord[i] = stl_mesh->coords[i];

  fvm_io_num_t *v_io_num = fvm_io_num_create_from_sfc(_face_vtx_coord,
                                                      3,
                                                      n_coo_ini,
                                                      FVM_IO_NUM_SFC_MORTON_BOX);

  BFT_FREE(_face_vtx_coord);

  cs_gnum_t *vtx_gnum = fvm_io_num_transfer_global_num(v_io_num);

  cs_lnum_t *order = cs_order_gnum(NULL, vtx_gnum, n_coo_ini);

  v_io_num = fvm_io_num_destroy(v_io_num);

  cs_coord_t  *vertex_coord = NULL;
  cs_lnum_t  *vertex_num = NULL;

  if (cs_glob_rank_id < 1) {

    BFT_MALLOC(vertex_coord, n_coo_ini*3, cs_coord_t);
    BFT_MALLOC(vertex_num, n_tria*3, cs_lnum_t);

    for (cs_lnum_t i = 0; i < n_coo_ini; i++)
      vertex_num[i] = -1;

    cs_lnum_t n_coo_new = 0;
    for (cs_lnum_t i = 0; i < n_coo_ini; i++) {
      cs_lnum_t j = order[i];
      vertex_num[j] = n_coo_new + 1;
      bool identical = false;
      if (! identical) {
        for (cs_lnum_t k = 0; k < 3; k++)
          vertex_coord[n_coo_new*3 + k] = stl_mesh->coords[j*3 + k];
        n_coo_new++;
      }
    }

    BFT_REALLOC(vertex_coord, n_coo_new*3, cs_coord_t);

  }
#endif

  cs_coord_3_t *vertex_coord = NULL;
  cs_lnum_t  *vertex_num   = NULL;
  cs_gnum_t  *vertex_gnum  = NULL;
  cs_gnum_t  *faces_gnum   = NULL;

  if (cs_glob_rank_id < 1) {
    BFT_MALLOC(vertex_coord, n_tria*3, cs_coord_3_t);
    BFT_MALLOC(vertex_num  , n_tria*3, cs_lnum_t);
    BFT_MALLOC(vertex_gnum , n_tria*3, cs_gnum_t);
    BFT_MALLOC(faces_gnum  , n_tria,   cs_gnum_t);

    for (cs_lnum_t j = 0; j < n_tria*3; j++) {
      vertex_num[j] = j+1;
      vertex_gnum[j] = j+1;
      for (cs_lnum_t k = 0; k < 3; k++)
        vertex_coord[j][k] = stl_mesh->coords[j][k];
    }
    for (cs_lnum_t j = 0; j < n_tria; j++)
      faces_gnum[j] = j+1;

  }

  fvm_nodal_t *ext_mesh = fvm_nodal_create(stl_mesh->name, 3);

  /* Generate FVM structure */

  fvm_nodal_append_by_transfer(ext_mesh,
                               n_tria,
                               FVM_FACE_TRIA,
                               NULL,
                               NULL,
                               NULL,
                               vertex_num,
                               NULL);

  if (cs_glob_rank_id < 1) {
    fvm_nodal_set_shared_vertices(ext_mesh, (const cs_coord_t *)stl_mesh->coords);
  } else {
    fvm_nodal_set_shared_vertices(ext_mesh, NULL);
  }

  fvm_nodal_init_io_num(ext_mesh, faces_gnum , 2);
  fvm_nodal_init_io_num(ext_mesh, vertex_gnum, 0);

  stl_mesh->ext_mesh = ext_mesh;

  if (cs_glob_rank_id < 1) {
    BFT_FREE(vertex_gnum);
    BFT_FREE(faces_gnum);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply a transformation matrix to a STL mesh structure.
 *
 * \param[in]  stl_mesh        pointer to the associated STL mesh structure
 * \param[in]  matrix          transformation matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_stl_mesh_transform(cs_stl_mesh_t  *stl_mesh,
                      double          matrix[3][4])
{
  cs_lnum_t n_tria = stl_mesh->n_faces;
  cs_lnum_t n_points = n_tria*3;

  for (cs_lnum_t i = 0; i < n_points; i++) {

    cs_real_t *c = stl_mesh->coords[i];

    double  c_a[4] = {c[0], c[1], c[2], 1.}; /* homogeneous coords */
    double  c_b[3] = {0, 0, 0};

    for (cs_lnum_t j = 0; j < 3; j++)
      for (int k = 0; k < 4; k++)
        c_b[j] += matrix[j][k]*c_a[k];

    for (cs_lnum_t j = 0; j < 3; j++)
      c[j] = c_b[j];

  }

  _update_init_mesh(stl_mesh);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply a transformation matrix to a STL mesh structure, but use
 * \brief the initial coordinates as inputs
 *
 * \param[in]  stl_mesh        pointer to the associated STL mesh structure
 * \param[in]  matrix          transformation matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_stl_mesh_transform_from_init(cs_stl_mesh_t  *stl_mesh,
                                double          matrix[3][4])
{
  cs_lnum_t n_tria = stl_mesh->n_faces;
  cs_lnum_t n_points = n_tria*3;

  for (cs_lnum_t i = 0; i < n_points; i++) {

    cs_real_t *c = stl_mesh->coords[i];
    cs_real_t *ci = stl_mesh->coords_ini[i];

    double  c_a[4] = {ci[0], ci[1], ci[2], 1.}; /* homogeneous coords */
    double  c_b[3] = {0, 0, 0};

    for (cs_lnum_t j = 0; j < 3; j++)
      for (int k = 0; k < 4; k++)
        c_b[j] += matrix[j][k]*c_a[k];

    for (cs_lnum_t j = 0; j < 3; j++)
      c[j] = c_b[j];

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply a translation to a STL mesh structure.
 *
 * \param[in]  stl_mesh        pointer to the associated STL mesh structure
 * \param[in]  vector          translation vector
 */
/*----------------------------------------------------------------------------*/

void
cs_stl_mesh_translate(cs_stl_mesh_t  *stl_mesh,
                      cs_real_t      vector[3])
{
  cs_lnum_t n_tria = stl_mesh->n_faces;
  cs_lnum_t n_points = n_tria*3;

  for (cs_lnum_t i = 0; i < n_points; i++) {
    for (cs_lnum_t j = 0; j < 3; j++) {
      stl_mesh->coords[i][j] += vector[j];
    }
  }

  _update_init_mesh(stl_mesh);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply a rotation to a STL mesh structure.
 *
 * \param[in]  stl_mesh        pointer to the associated STL mesh structure
 * \param[in]  theta           rotation angle
 * \param[in]  axis            rotation axis
 * \param[in]  center          rotation center
 */
/*----------------------------------------------------------------------------*/

void
cs_stl_mesh_rotate(cs_stl_mesh_t  *stl_mesh,
                   double          theta,
                   double          axis[3],
                   double          center[3])
{
  double matrix[3][4];

  cs_rotation_matrix(theta, axis, center, matrix);

  cs_stl_mesh_transform(stl_mesh, matrix);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply a scaling to a STL mesh structure.
 *
 * \param[in]  stl_mesh        pointer to the associated STL mesh structure
 * \param[in]  vector          translation vector
 */
/*----------------------------------------------------------------------------*/

void
cs_stl_mesh_scale(cs_stl_mesh_t  *stl_mesh,
                  double          scale)
{
  cs_lnum_t n_tria = stl_mesh->n_faces;
  cs_lnum_t n_points = n_tria*3;

  for (cs_lnum_t i = 0; i < n_points; i++) {
    for (cs_lnum_t j = 0; j < 3; j++) {
      stl_mesh->coords[i][j] *= scale;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the points outside he STL geometry. Those points will be used as
 * \brief seeds to propagate porosity values outside the STL geometry.
 *
 * \param[in]  stl_mesh  pointer to the associated STL mesh structure
 * \param[in]  n_points  number of points
 * \param[in]  coords    coordinates (x1,y1,z1,...,xn,yn,zn) (size : 3*n_point)
 */
/*----------------------------------------------------------------------------*/

void
cs_stl_set_porosity_seed(cs_stl_mesh_t  *stl_mesh,
                         int            n_points,
                         cs_real_t      *coords)
{
  stl_mesh->n_seeds = n_points;
  BFT_REALLOC(stl_mesh->seed_coords, n_points*3, cs_real_t);

  for (int i = 0; i < 3*n_points; i++)
    stl_mesh->seed_coords[i] = coords[i];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return writer_id used for stl meshes. 0 means unused.
 */
/*----------------------------------------------------------------------------*/

int
cs_stl_post_get_writer_id(void)
{

  return _stl_meshes.writer_id;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a new writer that will contains the STL mesh added by the user
 * \brief The writer_id is stored in the cs_glob_stl_meshes structure.
 *
 * \param[in]  time_dep > 1 if the writer is transient, else writer is fixed
 */
/*----------------------------------------------------------------------------*/

void
cs_stl_post_init_writer(const char             *case_name,
                        const char             *dir_name,
                        const char             *fmt_name,
                        const char             *fmt_opts,
                        fvm_writer_time_dep_t   time_dep,
                        bool                    output_at_start,
                        bool                    output_at_end,
                        int                     frequency_n,
                        double                  frequency_t)
{
  /* We check if a writer_id has already been defined.*/
  bool writer_exists = false;
  if ( _stl_meshes.writer_id != 0)
    writer_exists = true;

  /* If a writer do no already exist, create it */
  if (!writer_exists) {
    int writer_id = cs_post_get_free_writer_id();
    _stl_meshes.writer_id = writer_id;

    /* Add writer  */
    cs_post_define_writer(writer_id,       /* writer_id */
                          case_name,       /* writer name */
                          dir_name,        /* directory name */
                          fmt_name,        /* format_name */
                          fmt_opts,        /* format_options */
                          time_dep,
                          output_at_start, /* output_at_start */
                          output_at_end,   /* output_at_end */
                          frequency_n,     /* frequency_n */
                          frequency_t);    /* frequency_t */
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Associate a STL mesh to the default writer
 *
 * \param[in]  stl_mesh  pointer to the associated STL mesh structure
 */
/*----------------------------------------------------------------------------*/

void
cs_stl_post_add_mesh(cs_stl_mesh_t  *stl_mesh)
{
  if (_stl_meshes.writer_id == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("No writer was defined for STL mesh output\n"
                "cs_stl_post_init_writer should be called first.\n"));

  int writer_ids[] = {_stl_meshes.writer_id};
  int stl_mesh_id = cs_post_get_free_mesh_id();
  cs_post_define_existing_mesh(stl_mesh_id,
                               stl_mesh->ext_mesh,
                               0,
                               true,
                               false,
                               1,
                               writer_ids);

  cs_post_write_meshes(NULL);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write a binary STL file according to a given STL mesh structure
 *
 * \param[in]  stl_mesh  pointer to the associated STL mesh structure
 * \param[in]  path      path to the STL file
 */
/*----------------------------------------------------------------------------*/

void
cs_stl_file_write(cs_stl_mesh_t  *stl_mesh,
                  const char     *path)
{
  /* Output only on rank 0 for now */
  if (cs_glob_rank_id > 0)
    return;

  uint8_t buf[128];
  FILE *fp;

  /* Open STL file */
  fp = fopen(path,"wb");
  if (fp == NULL) {
      bft_error(__FILE__, __LINE__, 0,
                _("Error opening file \"%s\":\n\n"
                  "  %s"), path, strerror(errno));
  }

  /* Write header */
  char header[80] = "Exported from code_saturne";
  memcpy(buf, header, 80);

  /* Cut number of triangles in 4 8bytes unsigned int */
  uint32_t ntri = (uint32_t)stl_mesh->n_faces;
  _cut32(buf+80,ntri);

  fwrite(buf, 84, 1, fp);

  /* Loop on each triangle */
  for (int i = 0; i < stl_mesh->n_faces; i++) {

    uint8_t *start = buf + 12;

    /* Compute and Write the face normal coordinates */
    float normals[3];
    float a[3], b[3], c[3];

    for (int dir = 0; dir < 3; dir ++) {
      a[dir] = (float)stl_mesh->coords[3*i + 0][dir];
      b[dir] = (float)stl_mesh->coords[3*i + 1][dir];
      c[dir] = (float)stl_mesh->coords[3*i + 2][dir];
    }

    // ab vect ac
    normals[0] = (b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1]);
    normals[1] = (b[2]-a[2])*(c[0]-a[0]) - (b[0]-a[0])*(c[2]-a[2]);
    normals[2] = (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);

    float norm = sqrt(  normals[0]*normals[0]
                      + normals[1]*normals[1]
                      + normals[2]*normals[2]);

    for (int dir = 0; dir < 3; dir ++)
      normals[dir] /= norm;

    uint32_t n_temp[3];
    for (int dir = 0; dir < 3; dir ++) {
      memcpy(&n_temp[dir],
             &normals[dir],
             sizeof n_temp[dir]);

      _cut32(buf + 4*dir, n_temp[dir]);
    }

    /* Write the 3 vertices for the current triangle */
    for (int vi = 0; vi < 3; vi ++) {

      uint32_t v_temp[3];
      for (int dir = 0; dir < 3; dir ++) {
        float loc_coord = (float)stl_mesh->coords[3*i + vi][dir];
        memcpy(&v_temp[dir],
               &loc_coord,
               sizeof v_temp[dir]);
        _cut32(start + 3*4*vi + 4*dir, v_temp[dir]);
      }
    }

    fwrite(buf, 50, 1, fp);
  }

  fclose(fp);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute intersection between a STL mesh and the main mesh.
 *
 * \param[in]     stl_mesh         pointer to the associated STL mesh structure
 * \param[in]     n_input          number of cells on which intersection is done
 * \param[in]     input_idx        index of input cells (size: input_idx)
 * \param[out]    n_selected_cells number of output intersecting cells
 * \param[out]    selected_cells   index of output cells (size: output_idx)
 * \param[out]    tria_in_cell_idx start index of triangle intersecting each cell
 *                                  (size: n_output)
 * \param[out]    tria_in_cell_lst list of triangles in intersecting cells
 * \param[in,out] max_size         maximum size of tria_in_cell_lst array
 */
/*----------------------------------------------------------------------------*/

void
cs_stl_intersection(cs_stl_mesh_t *stl_mesh,
                    cs_lnum_t     n_input,
                    cs_lnum_t     *input_idx,
                    cs_lnum_t     *n_selected_cells,
                    cs_lnum_t     *selected_cells,
                    cs_lnum_t     *tria_in_cell_idx,
                    cs_lnum_t     **tria_in_cell_lst,
                    cs_lnum_t     *max_size)
{
  cs_mesh_t *m = cs_glob_mesh;

  cs_lnum_t  n_cells = n_input;// Local number of cells of the main mesh
  cs_lnum_t  n_tria_stl = stl_mesh->n_faces; // Local number of triangles of the STL
  cs_lnum_t  n_boxes = n_cells + n_tria_stl;
  int dim = 3;

  cs_lnum_t *_tria_in_cell_lst = NULL;
  if (tria_in_cell_lst != NULL)
    _tria_in_cell_lst = *tria_in_cell_lst;

  cs_gnum_t *box_gnum = NULL;
  cs_coord_t *extents = NULL;

  BFT_MALLOC(box_gnum, n_boxes, cs_gnum_t);
  BFT_MALLOC(extents , 2*dim*n_boxes, cs_coord_t);

  /* Global numbering construction */
  for (cs_lnum_t i = 0; i < n_cells; i++)
    box_gnum[i] = i + 1;

  for (cs_lnum_t i = 0; i < n_tria_stl; i++)
    box_gnum[i + n_cells] = n_cells + i + 1;


  /* Compute extents */

  // For the main mesh
  cs_real_6_t *bbox = NULL;
  bbox = cs_mesh_quantities_cell_extents(m, 0.0);

  for (cs_lnum_t i = 0; i < n_cells; i++)
    for (int id = 0; id < 6; id ++)
      extents[6*i + id] = bbox[input_idx[i]][id];

  BFT_FREE(bbox);

  // For the STL mesh
  BFT_MALLOC(bbox, n_tria_stl, cs_real_6_t);

  for (cs_lnum_t i = 0; i < n_tria_stl; i++) {
    bbox[i][0] = HUGE_VAL;
    bbox[i][1] = HUGE_VAL;
    bbox[i][2] = HUGE_VAL;
    bbox[i][3] = -HUGE_VAL;
    bbox[i][4] = -HUGE_VAL;
    bbox[i][5] = -HUGE_VAL;
  }

  for (cs_lnum_t i = 0; i < n_tria_stl; i++) {
    for (int vtx_id = 0; vtx_id < 3; vtx_id ++) {
      bbox[i][0] = CS_MIN(bbox[i][0], stl_mesh->coords[3*i + vtx_id][0]);
      bbox[i][3] = CS_MAX(bbox[i][3], stl_mesh->coords[3*i + vtx_id][0]);
      bbox[i][1] = CS_MIN(bbox[i][1], stl_mesh->coords[3*i + vtx_id][1]);
      bbox[i][4] = CS_MAX(bbox[i][4], stl_mesh->coords[3*i + vtx_id][1]);
      bbox[i][2] = CS_MIN(bbox[i][2], stl_mesh->coords[3*i + vtx_id][2]);
      bbox[i][5] = CS_MAX(bbox[i][5], stl_mesh->coords[3*i + vtx_id][2]);
    }
  }

  for (cs_lnum_t i = 0; i < n_tria_stl; i++) {
    for (cs_lnum_t id = 0; id < 6; id ++)
      extents[6*(i+n_cells) + id] = bbox[i][id];
  }

  BFT_FREE(bbox);

  fvm_neighborhood_t  *cell_neighborhood = NULL;

#if defined HAVE_MPI
  cell_neighborhood = fvm_neighborhood_create(MPI_COMM_NULL);
#else
  cell_neighborhood = fvm_neighborhood_create();
#endif

  // Compute intersections by boxes (transfer ownership)
  fvm_neighborhood_by_boxes( cell_neighborhood,
                             dim,
                             n_boxes,
                             box_gnum,
                             extents,
                             NULL,
                             NULL );

  cs_lnum_t  n_elts = 0;
  cs_gnum_t *elt_num = NULL;
  cs_lnum_t *neighbor_index = NULL;
  cs_gnum_t *neighbor_num = NULL;

  // Get the data back
  fvm_neighborhood_transfer_data( cell_neighborhood,
                                  &n_elts,
                                  &elt_num,
                                  &neighbor_index,
                                  &neighbor_num );

  cs_gnum_t _n_cells = n_cells;
  cs_lnum_t _n_selected_cells = 0;
  cs_lnum_t idx_num = 0;

  /* Init of list if needed */
  if (tria_in_cell_idx != NULL)
    tria_in_cell_idx[_n_selected_cells] = 0;

  /* Loop on the elements that have neighbors */
  for (cs_lnum_t i = 0; i < n_elts; i++) {

    // Check if the element is a box from the main mesh
    if (elt_num[i] <= _n_cells) {

      cs_lnum_t cell_id = elt_num[i] - 1;
      cs_lnum_t nb_tri_in_cell = 0;

      // Loop on the neighbors
      for (cs_lnum_t j = neighbor_index[i]; j < neighbor_index[i+1]; j++) {

        // If the current neighbor is a box from the STL,
        // -> the BB from the cell intersects the BB from
        // the STL
        if (neighbor_num[j] > _n_cells) {

          // Local number of the triangle
          cs_lnum_t ntria_id = neighbor_num[j] - n_cells - 1;

          cs_real_t tri_coords[3][3];
          for (cs_lnum_t k = 0; k < 3; k++) {
            for (cs_lnum_t l = 0; l < 3; l++)
              tri_coords[k][l] = stl_mesh->coords[3*ntria_id + k][l];
          }

          // Additional tests to know weather the triangle
          // overlap the cell bounding box
          bool triangle_box_intersect
            = _triangle_box_intersect(&extents[6*cell_id],
                                      tri_coords);


          if (triangle_box_intersect) {
            nb_tri_in_cell++;
            if (tria_in_cell_lst != NULL) {
              if (idx_num >= *max_size) {
                *max_size *= 2;
                BFT_REALLOC(_tria_in_cell_lst, *max_size, cs_lnum_t);
              }
              _tria_in_cell_lst[idx_num] = ntria_id;
              idx_num ++;
            }

          }

        }

      } // End loop on neighbors

      // If at least one time a triangle neigbor was found in the
      // previous loop, tag the cell.
      if (nb_tri_in_cell > 0) {
        selected_cells[_n_selected_cells] = input_idx[cell_id];
        if (tria_in_cell_idx != NULL)
          tria_in_cell_idx[_n_selected_cells + 1] = idx_num;
        _n_selected_cells ++;
      }

    }
  }

  *n_selected_cells  = _n_selected_cells;
  if (tria_in_cell_lst != NULL)
    *tria_in_cell_lst = _tria_in_cell_lst;

  fvm_neighborhood_destroy(&cell_neighborhood);

  BFT_FREE(elt_num);
  BFT_FREE(neighbor_index);
  BFT_FREE(neighbor_num);
  BFT_FREE(box_gnum);
  BFT_FREE(extents);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Refine the mesh following a given STL mesh
 *
 * \param[in]  stl_mesh        pointer to the associated STL mesh structure
 * \param[in]  n_ref           level of refinement
 * \param[in]  n_add_layer     additional layers between two refinement stage
 */
/*----------------------------------------------------------------------------*/

void
cs_stl_refine(cs_stl_mesh_t *stl_mesh,
              int           n_ref,
              int           n_add_layer)
{
  cs_mesh_t *m = cs_glob_mesh;

  /* Initialisation of the selection criteria (first level) */
  cs_lnum_t n_input_cells = m->n_cells;
  cs_lnum_t *input_cells = NULL;
  BFT_MALLOC(input_cells, m->n_cells, cs_lnum_t);
  for (cs_lnum_t i = 0; i < m->n_cells; i++)
    input_cells[i] = i;

   /* Begining of the refinement loop */
  for (int n_level = 0; n_level < n_ref+1; n_level ++) {

    if (cs_glob_rank_id < 1)
      bft_printf("Refinement %d\n",n_level);

    cs_lnum_t n_selected_cells = 0;
    cs_lnum_t *selected_cells = NULL;
    BFT_MALLOC(selected_cells, m->n_cells, cs_lnum_t);

    /* Select previous refined cells */
    if (n_level > 0) {
       BFT_REALLOC(input_cells, m->n_cells, cs_lnum_t);

      char criteria[100];
      sprintf(criteria,"STL_refined_region_%d",n_level-1);
      cs_selector_get_cell_list(criteria,
                                &n_input_cells,
                                input_cells  );
    }

    /* Compute mesh/STL intersection  */
    cs_stl_intersection( stl_mesh,
                         n_input_cells,
                         input_cells,
                         &n_selected_cells,
                         selected_cells,
                         NULL,
                         NULL,
                         NULL);

    /* If no intersected cells, do not perform refinement */
    cs_lnum_t n_intersected_cells = n_selected_cells;
    cs_parall_sum(1, CS_LNUM_TYPE, &n_intersected_cells);
    if (n_intersected_cells == 0)
      bft_error(__FILE__, __LINE__, 0,
                "Error in function cs_stl_refine: no intersection\n"
                "detected between the given STL file and the main mesh \n");

    /* Add incremented group */
    char group_name[100];
    sprintf(group_name,"STL_refined_region_%d",n_level);
    cs_mesh_group_cells_add(m,
                            group_name,
                            n_selected_cells,
                            selected_cells);

    /* Perform refinement */
    if (n_level < n_ref) {

      int *cell_tag;
      BFT_MALLOC(cell_tag, m->n_cells_with_ghosts, int);

      for (cs_lnum_t c_id = 0; c_id < m->n_cells; c_id ++) {
        cell_tag[c_id] = 0;
        for (cs_lnum_t i = 0; i < n_selected_cells; i++) {
          if (c_id == selected_cells[i]) cell_tag[c_id] = 1;
        }
      }

      if (m->halo!=NULL) {
        cs_halo_set_use_barrier(1);
        cs_halo_sync_num(m->halo, CS_HALO_STANDARD, cell_tag);
      }

      /* Here we add additionnal layers of refined cells
       * around the original STL selected cells */
      int nn = 1;

      for (int k = 0; k < n_add_layer ; k++){
        for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++){
          cs_lnum_t c1 = m->i_face_cells[face_id][0];
          cs_lnum_t c2 = m->i_face_cells[face_id][1];
          if (cell_tag[c1] == 0 && cell_tag[c2] == nn)
            cell_tag[c1] = nn+1;
          if (cell_tag[c2] == 0 && cell_tag[c1] == nn)
            cell_tag[c2] = nn+1;
        }

        if (m->halo!=NULL)
          cs_halo_sync_num(m->halo, CS_HALO_STANDARD, cell_tag);
        nn ++;
      }

      n_selected_cells = 0;
      for (cs_lnum_t c_id = 0; c_id < m->n_cells; c_id ++) {
        if (cell_tag[c_id] > 0) {
          selected_cells[n_selected_cells] = c_id;
          n_selected_cells ++;
        }
      }

      BFT_FREE(cell_tag);

      cs_mesh_refine_simple_selected(m,
                                     false,
                                     n_selected_cells,
                                     selected_cells);
    }

    BFT_FREE(selected_cells);

    /* Every 2 refinemet stages and at the end, re-partition
     * the mesh to limit imbalance in the processors load */

    if ((n_level%2 == 0 || n_level == n_ref -1)
        && cs_glob_rank_id > -1) {

      cs_mesh_builder_destroy(&cs_glob_mesh_builder);
      cs_glob_mesh_builder = cs_mesh_builder_create();
      cs_mesh_to_builder(m, cs_glob_mesh_builder, true, NULL);

      cs_partition(m, cs_glob_mesh_builder, CS_PARTITION_MAIN);

      cs_mesh_from_builder(m, cs_glob_mesh_builder);
      cs_mesh_init_halo(m, cs_glob_mesh_builder, m->halo_type, -1, true);
    }

    cs_mesh_update_auxiliary(m);
  }

  BFT_FREE(input_cells);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute porosity field according to a given STL mesh
 *
 * \param[in]  stl_mesh      pointer to the associated STL mesh structure
 * \param[out] porosity      interpolated porosity field
 * \param[out] indic         indicator of the STL location
 */
/*----------------------------------------------------------------------------*/

void
cs_stl_compute_porosity(cs_stl_mesh_t *stl_mesh,
                        cs_real_t     *porosity,
                        int           *indic)
{
  cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  cs_stl_mesh_t *stl = stl_mesh;

  const cs_real_t *volume   = (const cs_real_t *)mq->cell_vol;
  const cs_real_3_t *xyzcen = (const cs_real_3_t *)mq->cell_cen;
  const cs_real_3_t *cdgfac = (const cs_real_3_t *)mq->i_face_cog;
  const cs_real_3_t *cdgfbo = (const cs_real_3_t *)mq->b_face_cog;
  const cs_lnum_2_t *ifacel = (const cs_lnum_2_t *)m->i_face_cells;
  const cs_lnum_t *ifabor   = (const cs_lnum_t *)m->b_face_cells;

  /* If no reference point outside the STL file is given
   * by the user, the porosity computation will fail */
  if (stl_mesh->seed_coords == NULL || stl_mesh->n_seeds == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Error in STL porosity computation: no seed points"
                "for STL mesh %s is given."), stl_mesh->name);

  /* Initialisation */
  cs_lnum_t n_input_cells      = m->n_cells;
  cs_lnum_t n_selected_cells   = 0;
  cs_lnum_t *input_cells       = NULL;
  cs_lnum_t *selected_cells    = NULL;
  cs_lnum_t *tria_in_cell_lst  = NULL;
  cs_lnum_t *tria_in_cell_idx  = NULL;
  cs_lnum_t *cell_selected_idx = NULL;
  cs_real_t *mean_plane_def    = NULL;
  cs_real_t *stl_normals       = NULL;

  cs_lnum_t max_size = stl->n_faces;

  BFT_MALLOC(input_cells      , m->n_cells      , cs_lnum_t);
  BFT_MALLOC(selected_cells   , m->n_cells      , cs_lnum_t);
  BFT_MALLOC(tria_in_cell_idx , m->n_cells      , cs_lnum_t);
  BFT_MALLOC(tria_in_cell_lst , max_size        , cs_lnum_t);
  BFT_MALLOC(cell_selected_idx, m->n_cells      , cs_lnum_t);
  BFT_MALLOC(stl_normals      , stl->n_faces*3  , cs_real_t);
  BFT_MALLOC(mean_plane_def   , m->n_cells*6    , cs_real_t);

  /* Input cells are all[] cells and initialisation */
  for (cs_lnum_t i = 0; i < m->n_cells; i++) {
    input_cells[i] = i;
    cell_selected_idx[i] = -1;
    porosity[i] = 0.0;
  }

  /* Compute normals */
  for (cs_lnum_t i = 0; i < stl->n_faces; i++)
    _compute_normal(stl_mesh->coords + (3*i), stl_normals + (3*i));

  /* Get the intersection
   * ==================== */
  cs_stl_intersection(stl_mesh,
                      n_input_cells,
                      input_cells,
                      &n_selected_cells,
                      selected_cells,
                      tria_in_cell_idx,
                      &tria_in_cell_lst,
                      &max_size);

  /* Compute the bounding boxes of the main mesh */
  cs_real_6_t *bbox = NULL;
  bbox = cs_mesh_quantities_cell_extents(m, 0.0);

  /* If a cell is overlaped by more than 1 triangle,
   * replace those triangles by a mean plane
   * ============================================== */
  for (cs_lnum_t i = 0; i < n_selected_cells; i ++) {

    cs_lnum_t start_id = tria_in_cell_idx[i];
    cs_lnum_t end_id   = tria_in_cell_idx[i+1];
    cs_lnum_t n_tria_in_cell = end_id - start_id;

    cs_lnum_t cell_id   = selected_cells[i];

    if (n_tria_in_cell > 1) {

      cs_real_t surftot  = 0;

      // Mean plane defined by normal and point
      cs_real_t n_plan[3], p_plan[3];
      for (int k = 0; k < 3; k++) {
        n_plan[k] = 0.;
        p_plan[k] = 0.;
      }

      /* Loop on triangles intersecting the current cell */
      for (cs_lnum_t tri = start_id; tri < end_id; tri ++) {

        cs_lnum_t f_tria_id = tria_in_cell_lst[tri];

        cs_real_3_t *coords  = stl_mesh->coords + (3*f_tria_id);
        cs_real_t *n         = stl_normals + (3*f_tria_id);
        cs_real_t *aabb      = &bbox[cell_id][0];
        const cs_real_t *cog = xyzcen[cell_id];

        /* Compute the surface of the triangle in the box */
        cs_real_t surf = _exact_triangle_box_surface_intersection(aabb, coords);
        surftot += surf;

        /* Center of gravity of triangle */
        cs_real_t cog_tri[3];
        for (int k = 0; k < 3; k++)
          cog_tri[k] = (coords[0][k] + coords[1][k] + coords[2][k]) / 3.0;

        double psca = cs_math_3_distance_dot_product(cog, cog_tri, n);
        psca *= surf;

        /* Weigh triangle normal and point by the
         * surface included in the box  */
        for (int k = 0; k < 3; k++) {
          n_plan[k] += surf * n[k];
          p_plan[k] += psca * n[k];

        }

      } // End loop on tri intersecting the cell

      /* Normalisation */
      cs_real_t nn = cs_math_3_norm(n_plan);

      for (int k = 0; k < 3; k++) {
        n_plan[k] /= cs_math_fmax(nn, 1.e-20);
        p_plan[k] /= cs_math_fmax(surftot, 1.e-20);
      }

      /* Mean plane definition storage */
      for (int k = 0; k < 3; k++) {
        mean_plane_def[6*i + k] = p_plan[k] + xyzcen[cell_id][k];
        mean_plane_def[6*i + 3 + k] = n_plan[k];
      }


    /* else if only one triangle in cell */
    } else {
      cs_lnum_t f_tria_id = tria_in_cell_lst[start_id];

      for (int k = 0; k < 3; k++) {
        mean_plane_def[6*i + k] = stl_mesh->coords[3*f_tria_id][k];
        mean_plane_def[6*i + 3 + k] = stl_normals[3*f_tria_id + k];
      }

    }// End test on n_tria_in_cell

  // Building indirection
  cell_selected_idx[cell_id] = i;
  } // End loop on selected_cells

  /* Porosity computation on intersecting cells
   * ========================================== */

  cs_real_t vol;

  /* Loop on internal faces */
  for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id ++) {

    /* Loop on adjacent cells */
    for (int k = 0; k < 2; k++) {
      cs_lnum_t cell_id = ifacel[face_id][k];
      cs_lnum_t idx;

      /* Test if the cell is ghost cell */
      if (cell_id < m->n_cells) {
        idx = cell_selected_idx[cell_id];

        /* Test if the cell was a previously selected
         * intersecting cell*/
        if (idx > -1) {
          cs_real_t xf[3], xc[3];

          for (int i = 0; i < 3; i++) {
            xf[i] = cdgfac[face_id][i];
            xc[i] = xyzcen[cell_id][i];
          }

          cs_real_t *tria_mean_plane = &mean_plane_def[6*idx];

          /* For each face, tetrahedra based one edge, the cog of
           * the face and the cog of the cell, are built*/
          cs_lnum_t vtx_start = m->i_face_vtx_idx[face_id];
          cs_lnum_t vtx_end   = m->i_face_vtx_idx[face_id + 1];
          cs_lnum_t nb_vtx    = vtx_end - vtx_start;

          for (cs_lnum_t vtx = vtx_start; vtx < vtx_end; vtx++) {

            // When vtx = vtx_end-1, vtx+1 has to point towards vtx_start
            cs_lnum_t vtx1 = vtx;
            cs_lnum_t vtx2 = vtx_start + (vtx+1-vtx_start) % nb_vtx;

            cs_lnum_t vtx1_id = m->i_face_vtx_lst[vtx1];
            cs_lnum_t vtx2_id = m->i_face_vtx_lst[vtx2];

            cs_real_t *x1 = &(m->vtx_coord[3*vtx1_id]);
            cs_real_t *x2 = &(m->vtx_coord[3*vtx2_id]);

            vol = _tetrahedron_plane_volume_intersection(x1, x2, xf, xc, tria_mean_plane);

            porosity[cell_id] += vol;
          }

        } // End test if it is an intersecting cell

      } // End test if cell in not a ghost

    } // End loop on adjacent cells

  } // End loop on internal faces

  /* Loop on boundary faces */
  for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id ++) {

    cs_lnum_t cell_id = ifabor[face_id];
    cs_lnum_t idx;

    /* Test if the cell is ghost cell */
    if (cell_id < m->n_cells) {
      idx = cell_selected_idx[cell_id];

      /* Test if the cell was a previously selected
       * intersecting cell*/
      if (idx > -1) {
        cs_real_t xf[3], xc[3];

        for (int i = 0; i < 3; i++) {
          xf[i] = cdgfbo[face_id][i];
          xc[i] = xyzcen[cell_id][i];
        }

        cs_real_t *tria_mean_plane = &mean_plane_def[6*idx];

        /* For each face, tetrahedra based one edge, the cog of
         * the face and the cog of the cell, are built*/
        cs_lnum_t vtx_start = m->b_face_vtx_idx[face_id];
        cs_lnum_t vtx_end   = m->b_face_vtx_idx[face_id + 1];
        cs_lnum_t nb_vtx    = vtx_end - vtx_start;

        for (cs_lnum_t vtx = vtx_start; vtx < vtx_end; vtx++) {

          // When vtx = vtx_end-1, vtx+1 has to point towards vtx_start
          cs_lnum_t vtx1 = vtx;
          cs_lnum_t vtx2 = vtx_start + (vtx+1-vtx_start) % nb_vtx;

          cs_lnum_t vtx1_id = m->b_face_vtx_lst[vtx1];
          cs_lnum_t vtx2_id = m->b_face_vtx_lst[vtx2];

          cs_real_t *x1 = &(m->vtx_coord[3*vtx1_id]);
          cs_real_t *x2 = &(m->vtx_coord[3*vtx2_id]);

          vol = _tetrahedron_plane_volume_intersection(x1, x2, xf, xc, tria_mean_plane);

          porosity[cell_id] += vol;

        }

      } // End test if it is an intersecting cell

    } // End test if cell in not a ghost

  } // End loop on boundary faces

  /* Normalisation and clipping of the porosity field */
  for (cs_lnum_t i = 0; i < n_selected_cells; i++) {
    cs_lnum_t cell_id = selected_cells[i];
    porosity[cell_id] /= volume[cell_id];
    porosity[cell_id] = CS_MAX(porosity[cell_id],0.0);
    porosity[cell_id] = CS_MIN(porosity[cell_id],1.0);

    if (porosity[cell_id] < 1.e-5)
      porosity[cell_id] = 0.0;
  }

  if (m->halo!=NULL)
    cs_halo_sync_var(m->halo, CS_HALO_STANDARD, porosity);

  /* Porosity : filling the entire domain
   * ==================================== */

  int *cell_tag = NULL;
  BFT_MALLOC(cell_tag, m->n_cells_with_ghosts, int);

  /*-------------------------
   * cell_tag is :
   *  - 1 inside solid
   *  -  0 on the solid border
   *  - -1 outside the solid
   *  ------------------------*/

  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells_with_ghosts; cell_id ++) {
    cell_tag[cell_id] = 1;
    if (cell_id < m->n_cells)
      /* If the cell in a border of the STL mesh */
      if (cell_selected_idx[cell_id] > -1)
        cell_tag[cell_id] = 0;
  }


  /* Tag outside cells according to user presribed points */
  for (int p = 0; p < stl_mesh->n_seeds; p++) {

    cs_real_t *xyz_ref = &(stl_mesh->seed_coords[3*p]);
    cs_lnum_t ref_id;
    int ref_rank;

    cs_geom_closest_point(m->n_cells,
                          (const cs_real_3_t *)(mq->cell_cen),
                          xyz_ref,
                          &ref_id,
                          &ref_rank);

    if (ref_rank == cs_glob_rank_id)
      cell_tag[ref_id] = -1;
  }

  /* Sychronize */
  if (m->halo!=NULL)
    cs_halo_sync_num(m->halo, CS_HALO_STANDARD, cell_tag);

  /* Loop on internal faces "Algo glouton" */

  bool new_cells_found = true;

  while (new_cells_found) {

    cs_gnum_t cpt = 0;

    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id ++) {

      cs_lnum_t ii = ifacel[face_id][0];
      cs_lnum_t jj = ifacel[face_id][1];

      int ti = cell_tag[ii];
      int tj = cell_tag[jj];

      if (ti == 1 && tj == -1) {
        cell_tag[ii] = -1;
        cpt ++;
      }

      if (ti == -1 && tj == 1) {
        cell_tag[jj] = -1;
        cpt ++;
      }
    }

    if (cs_glob_rank_id >= 0)
      cs_parall_counter(&cpt, 1);

    if (cpt == 0)
      new_cells_found = false;

    if (m->halo!=NULL)
      cs_halo_sync_num(m->halo, CS_HALO_STANDARD, cell_tag);
  }

  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id ++) {
    if (indic != NULL)
      indic[cell_id] = cell_tag[cell_id];
    if (cell_tag[cell_id] == -1)
      porosity[cell_id] = 1.0;
  }

  BFT_FREE(cell_tag);
  BFT_FREE(cell_selected_idx);
  BFT_FREE(mean_plane_def);
  BFT_FREE(stl_normals);
  BFT_FREE(bbox);
  BFT_FREE(input_cells);
  BFT_FREE(selected_cells);
  BFT_FREE(tria_in_cell_lst);
  BFT_FREE(tria_in_cell_idx);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
