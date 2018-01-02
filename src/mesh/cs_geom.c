/*============================================================================
 * Geometry utility functions.
 *===========================================================================*/

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
 *---------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *---------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_math.h"
#include "cs_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *---------------------------------------------------------------------------*/

#include "cs_geom.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_geom.c
        Geometry utility functions.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/* In case <math.h> does not define HUGE_VAL, use a "safe" value */
#if !defined(HUGE_VAL)
#define HUGE_VAL 1.0e+30
#endif

/*============================================================================
 * Static global variables
 *===========================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *===========================================================================*/

/*============================================================================
 * Private function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Compute whether the projection point of prev_location is inside or outside
 * the a given edge of a sub triangle
 *
 * parameters:
 *   sx0    <-- Previous location of the particle
 *   sx1    <-- Next location of the particle
 *   vtx_0  <-- First vertex of the edge (sorted by index)
 *   vtx_1  <-- Second vertex of the edge (sorted by index)
 *
 * returns:
 *   1 if point is inside, -1 if outside
 *----------------------------------------------------------------------------*/

static int
_test_edge(const cs_real_t  sx0[3],
           const cs_real_t  sx1[3],
           const cs_real_t  vtx_0[3],
           const cs_real_t  vtx_1[3])
{
  /* vO vector where the choice for v between v0 and v1 has no importance
   * so we take the smaller vertex id */
  cs_real_3_t vO = {sx0[0] - vtx_0[0],
                    sx0[1] - vtx_0[1],
                    sx0[2] - vtx_0[2]};

  cs_real_3_t edge = {vtx_1[0] - vtx_0[0],
                      vtx_1[1] - vtx_0[1],
                      vtx_1[2] - vtx_0[2]};

  cs_real_3_t disp = {sx1[0] - sx0[0],
                      sx1[1] - sx0[1],
                      sx1[2] - sx0[2]};
  /* p = edge ^ vO */
  const cs_real_3_t p = {edge[1]*vO[2] - edge[2]*vO[1],
                         edge[2]*vO[0] - edge[0]*vO[2],
                         edge[0]*vO[1] - edge[1]*vO[0]};

  return (cs_math_3_dot_product(disp, p) > 0 ? 1 : -1);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief find the closest point of a set to a given point in space.
 *
 * If the orient parameter is set to -1 or 1, intersection is only
 * considered when (sx1-sx0).normal.orient > 0.
 * If set to 0, intersection is considered in both cases.
 *
 * \param[in]   n_points      number of points
 * \param[in]   point_coords  point coordinates
 * \param[in]   query_coords  coordinates searched for
 * \param[out]  point_id      id of closest point if on the same rank,
 *                            -1 otherwise
 * \param[out]  rank_id       id of rank containing closest point
 */
/*----------------------------------------------------------------------------*/

void
cs_geom_closest_point(cs_lnum_t         n_points,
                      const cs_real_t   point_coords[][3],
                      const cs_real_t   query_coords[3],
                      cs_lnum_t        *point_id,
                      int              *rank_id)
{
  cs_lnum_t id_min = -1;
  cs_real_t d2_min = HUGE_VAL;

  for (cs_lnum_t i = 0; i < n_points; i++) {
    cs_real_t d2 = cs_math_3_distance_squared(point_coords[i], query_coords);
    if (d2 < d2_min) {
      d2_min = d2;
      id_min = i;
    }
  }

  *rank_id = cs_glob_rank_id;

  cs_parall_min_id_rank_r(&id_min, rank_id, d2_min);

  if (*rank_id != cs_glob_rank_id)
    *point_id = -1;
  else
    *point_id = id_min;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Test if a line segment intersects a face.
 *
 * If the orient parameter is set to -1 or 1, intersection is only
 * considered when (sx1-sx0).normal.orient > 0.
 * If set to 0, intersection is considered in both cases.
 *
 * \param[in]   orient         if -1 or 1, multiplies face_normal to check
 *                             for segment
 * \param[in]   n_vertices     number of face vertices
 * \param[in]   vertex_ids     ids of face vertices
 * \param[in]   vertex_coords  vertex coordinates
 * \param[in]   face_center    cordinates of face center
 * \param[in]   face_normal    face normal vector
 * \param[in]   sx0            segment start coordinates
 * \param[in]   sx1            segment end coordinates
 * \param[out]  n_crossings    number sub_face crossings
 *                             [0: in; 1: out]
 *
 * \return
 *   2 if the segment does not go through the face's plane, or minimum
 *   relative distance (in terms of barycentric coordinates)
 *   of intersection point to face.
 */
/*----------------------------------------------------------------------------*/

double
cs_geom_segment_intersect_face(int              orient,
                               cs_lnum_t        n_vertices,
                               const cs_lnum_t  vertex_ids[],
                               const cs_real_t  vertex_coords[][3],
                               const cs_real_t  face_center[3],
                               const cs_real_t  face_normal[3],
                               const cs_real_t  sx0[3],
                               const cs_real_t  sx1[3],
                               int              n_crossings[2])
{
  const double epsilon = 1.e-15;

  /* Initialization of retval to unity*/
  double retval = 2.;

  assert(sizeof(cs_real_t) == 8);

  /* Initialization */

  const cs_real_t disp[3] = {sx1[0] - sx0[0],
                             sx1[1] - sx0[1],
                             sx1[2] - sx0[2]};
  const cs_real_t vgo[3] = {sx0[0] - face_center[0],
                            sx0[1] - face_center[1],
                            sx0[2] - face_center[2]};

  int n_intersects = 0;

  /* Principle:
   *  - loop on sub-triangles of the face
   *    and test for each triangle if the intersection is inside the triangle
   *  - use of a geometric algorithm:x
   *    the intersection is in the triangle if it is on the proper side of each
   *    three edges defining the triangle (e0, e1 and e_out)
   *    This is measured calculating the sign u, v and w
   *    (and keeping them in memory to calculate each value only once)
   *
   *        e0
   *          ---------
   *          |\  xI /|        I = intersection (occurring at t such that
   *          | \   / |                           I = O + t * OD          )
   *          |  \ /  |
   *    e_out |   x G |        G = Center of gravity of the face
   *          |  / \  |
   *          | /   \ |
   *          |/     \|
   *          ---------
   *        e1
   */

  /* Initialization of triangle points and edges (vectors) */
  cs_real_3_t  e0, e1;
  int pi, pip1, p0;

  /* 1st vertex: vector e0, p0 = e0 ^ vgo  */
  cs_lnum_t vtx_id_0 = vertex_ids[0];
  const cs_real_t *vtx_0 = vertex_coords[vtx_id_0];

  p0 = _test_edge(sx0, sx1, face_center, vtx_0);
  pi = p0;
  pip1 = p0;

  /* Loop on vertices of the face */
  for (cs_lnum_t i = 0; i < n_vertices; i++) {

    vtx_id_0 = vertex_ids[i];
    cs_lnum_t vtx_id_1 = vertex_ids[(i+1)%n_vertices];

    vtx_0 = vertex_coords[vtx_id_0];
    const cs_real_t *vtx_1 = vertex_coords[vtx_id_1];
    for (int j = 0; j < 3; j++) {
      e0[j] = vtx_0[j] - face_center[j];
      e1[j] = vtx_1[j] - face_center[j];
    }

    /* P = e1^e0 */

    const cs_real_3_t pvec = {e1[1]*e0[2] - e1[2]*e0[1],
                              e1[2]*e0[0] - e1[0]*e0[2],
                              e1[0]*e0[1] - e1[1]*e0[0]};

    double det = cs_math_3_dot_product(disp, pvec);
    int sign_det = det > 0 ? 1 : -1;

    /* 2nd edge: vector ei+1, pi+1 = ei+1 ^ vgo */

    pi = - pip1;
    if (i == n_vertices - 1)
      pip1 = p0;
    else
      pip1 = _test_edge(sx0, sx1, face_center, vtx_1);

    const int u_sign = pip1 * sign_det;

    /* 1st edge: vector ei, pi = ei ^ vgo */
    const int v_sign = pi * sign_det;

    /* 3rd edge: vector e_out */

    /* Check the orientation of the edge */
    int reorient_edge = (vtx_id_0 < vtx_id_1 ? 1 : -1);

    /* Sort the vertices of the edges so that it is easier to find it after */
    cs_lnum_2_t edge_id = {vtx_id_0, vtx_id_1};
    if (reorient_edge == -1) {
      edge_id[0] = vtx_id_1;
      edge_id[1] = vtx_id_0;
    }

    vtx_0 = vertex_coords[edge_id[0]];
    vtx_1 = vertex_coords[edge_id[1]];

    int w_sign = _test_edge(sx0, sx1, vtx_0, vtx_1)
                 * reorient_edge * sign_det;

    /* The projection of point O along displacement is outside of the triangle
     * then no intersection */
    if (w_sign > 0 || u_sign  < 0 || v_sign < 0)
      continue;

    /* We have an intersection if
     * u_sign >= 0, v_sign <= 0 and w_sign <= 0
     * If det is nearly 0, consider that the particle does not cross the face */

    double go_p = - cs_math_3_dot_product(vgo, pvec);
    int sign_go_p = (go_p > 0 ? 1 : -1);

    /* We check the direction of displacement and the triangle normal
     *  to see if the directed segment enters or leaves the half-space */

    bool dir_move;
    if (orient) {
      int sign_face_orient
        = (cs_math_3_dot_product(pvec, face_normal) > 0 ? 1 : -1 );
      dir_move = (sign_face_orient * orient * sign_det > 0);
    }
    else
      dir_move = true;

    /* Same sign (meaning there is a possible intersection with t > 0). */
    if (sign_det == sign_go_p) {

      /* The particle enters or leaves the half-space */
      if (dir_move) {
        if (fabs(go_p) < fabs(det)) {
          /* There is a real intersection (outward)
             with 0 < t < 1 (n_intersect ++) */
          double t = 0.99;

          const double det_cen = cs_math_3_norm(e0)*cs_math_3_norm(e1);
          if (fabs(det/det_cen) > epsilon) {
            t = go_p / det;
          }

          n_crossings[1] += 1;
          n_intersects += 1;
          if (t < retval)
            retval = t;
        } else {
          n_crossings[1] += 1;
        }
      } else {
        n_crossings[0] += 1;
        if (fabs(go_p) < fabs(det))
          n_intersects -= 1;
        /* There is a real intersection (inward) with 0 < t < 1 (n_intersect -) */
      }
    } else {
      /* Opposite sign (meaning there is a possible intersection with t < 0). */
      if (dir_move)
        n_crossings[1] += 1;
      else
        n_crossings[0] += 1;
    }

    /* In case intersections were removed due to non-convex cases
     *  (i.e.  n_intersects < 1, but retval <1 ),
     *  the retval value is forced to 1
     *  (no intersection since the particle entered and left from this face). */

    if (n_intersects < 1 && retval < 1.) {
      retval = 2.;
    }

  }

  return retval;
}

/*---------------------------------------------------------------------------*/

END_C_DECLS
