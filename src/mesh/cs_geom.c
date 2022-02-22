/*============================================================================
 * Geometric utility functions.
 *===========================================================================*/

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
 * Compute wether the (signed) volume of the parallelepiped defined by the
 * three vectors (a,b,c) is positive or not. Its is compulted with the triple
 * vector product (a x b) . c.
 *
 * Here:
 *  a = [vtx0, vtx1]
 *  b = [vtx0, sx0]
 *  c = [sx0, sx1]
 *
 * parameters:
 *   sx0    <-- Previous location of the particle
 *   sx1    <-- Next location of the particle
 *   vtx_0  <-- First vertex of the edge (sorted by index)
 *   vtx_1  <-- Second vertex of the edge (sorted by index)
 *
 * returns:
 *   1 if positive, -1 otherwise
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
    cs_real_t d2 = cs_math_3_square_distance(point_coords[i], query_coords);
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
 *                               |
 *      x------------------------|--------x D: end coordinates
 *   O: start coordinates        |
 *                               x G: Face (Center of Gravity)
 *             x current         |
 *               cell center     |
 *                               |
 *                           Face number

 * If the orient parameter is set to -1 or 1, intersection is only
 * considered when (sx1-sx0).normal.orient > 0.
 * If set to 0, intersection is considered in both cases.
 *
 * \param[in]      orient         if -1 or 1, multiplies face_normal to check
 *                                for segment
 * \param[in]      n_vertices     number of face vertices
 * \param[in]      vertex_ids     ids of face vertices
 * \param[in]      vtx_coord      vertex coordinates
 * \param[in]      face_cog       coordinates of face center
 * \param[in]      sx0            segment start coordinates
 * \param[in]      sx1            segment end coordinates
 * \param[out]     n_inout        number of sub_face crossings with segment's
 *                                half-line. [0: in; 1: out]
 * \param[in, out] face_norm      optional local face unit normal of the
 *                                crossed sub-triangle, or NULL if unused
 *
 * \return
 *   2 if the segment does not go through the face's plane, or minimum
 *   relative distance (in terms of segment length)
 *   of intersection point to face.
 */
/*----------------------------------------------------------------------------*/

double
cs_geom_segment_intersect_face(int               orient,
                               cs_lnum_t         n_vertices,
                               const cs_lnum_t   vertex_ids[],
                               const cs_real_t   vtx_coord[][3],
                               const cs_real_t   face_cog[3],
                               const cs_real_t   sx0[3],
                               const cs_real_t   sx1[3],
                               int               n_inout[2],
                               cs_real_t        *face_norm)
{
  const double epsilon = 1.e-15;

  /* Initialization of retval to two */
  double retval = 2.;

  assert(sizeof(cs_real_t) == 8);

  /* Initialization */

  const cs_real_t disp[3] = {sx1[0] - sx0[0],
                             sx1[1] - sx0[1],
                             sx1[2] - sx0[2]};
  const cs_real_t vgo[3] = {sx0[0] - face_cog[0],
                            sx0[1] - face_cog[1],
                            sx0[2] - face_cog[2]};

  int n_intersects = 0;
  bool in_cog = true;

  /* Principle:
   *  - loop on sub-triangles of the face
   *    and test for each triangle if the intersection is inside the triangle
   *  - use of a geometric algorithm:
   *    the intersection is in the triangle if it is on the proper side of all
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
   *
   *   Note that t belongs to [0,1] if there OD crosses the sub-triangle,
   *   t is negative if the line (OD) enters the volume,
   *   t is positive if the line (OD) leaves the volume (and though there is
   *    a risk of getting out)
   *   t is given by the formula:
   *    t = (OG . e1^e0) / (OD . e1^e0)
   */

  /* Initialization of triangle points and edges (vectors) */
  cs_real_3_t e0, e1;
  int pip1, p0;

  /* 1st vertex: test edge e0 with respect to (OD)
   * in an orthogonal plane of (OD)
   * vector e0 = [G, v0], p0 = e0 ^ OG
   * return "p0 . OD > 0 ?" */
  cs_lnum_t vtx_id_0 = vertex_ids[0];
  const cs_real_t *vtx_0 = vtx_coord[vtx_id_0];

  p0 = _test_edge(sx0, sx1, face_cog, vtx_0);
  pip1 = p0;

  /* Sign of "[OD].p" for the first sub-triangle.
   * We need to store this to check if all e_out are around
   * (OD), and for that the multiplier "[OD].p" must be the same
   * Otherwise we can have troubles when [OD] is orthogonal to the face normal
   * (it can switch sign from one triangle to an other...) */
  int sign_od_p0 = 0;

  /* Loop on vertices of the face */
  for (cs_lnum_t i = 0; i < n_vertices; i++) {

    vtx_id_0 = vertex_ids[i];
    cs_lnum_t vtx_id_1 = vertex_ids[(i+1)%n_vertices];

    vtx_0 = vtx_coord[vtx_id_0];
    const cs_real_t *vtx_1 = vtx_coord[vtx_id_1];
    for (int j = 0; j < 3; j++) {
      e0[j] = vtx_0[j] - face_cog[j];
      e1[j] = vtx_1[j] - face_cog[j];
    }

    /* P = e1^e0: same value for the two neighbouring cells
     * NB: in the other direction to the face normal */

    const cs_real_3_t pvec = {e1[1]*e0[2] - e1[2]*e0[1],
                              e1[2]*e0[0] - e1[0]*e0[2],
                              e1[0]*e0[1] - e1[1]*e0[0]};

    double od_p = cs_math_3_dot_product(disp, pvec);

    /* This sign is absolute (ie same result is obtained if the face is seen
     * from the other neighbouring cell).
     */
    int sign_od_p = (od_p > 0 ? 1 : -1);

    if (sign_od_p0 == 0)
      sign_od_p0 = sign_od_p;

    /* 2nd edge: vector ei+1, pi+1 = ei+1 ^ OG */

    int pi = pip1;
    if (i == n_vertices - 1)
      pip1 = p0;
    else
      pip1 = _test_edge(sx0, sx1, face_cog, vtx_1);

    int u_sign = pip1 * sign_od_p;

    /* 1st edge: vector ei, pi = ei ^ OG */
    int v_sign = - pi * sign_od_p;

    /* 3rd edge: vector e_out */

    /* Check the orientation of the edge */
    int reorient_edge = (vtx_id_0 < vtx_id_1 ? 1 : -1);

    /* Sort the vertices of the edges so that it gives the same
     * answer for the same edge for another face */
    cs_lnum_2_t edge_id = {vtx_id_0, vtx_id_1};
    if (reorient_edge == -1) {
      edge_id[0] = vtx_id_1;
      edge_id[1] = vtx_id_0;
    }

    vtx_0 = vtx_coord[edge_id[0]];
    vtx_1 = vtx_coord[edge_id[1]];

    int w_sign = _test_edge(sx0, sx1, vtx_0, vtx_1)
                 * reorient_edge * sign_od_p;

    /* Special treatment if the intersection point is
     * in G exactly (so w_sign is < 0 for all e_out)
     * Not: all e_out should be talen with the same noram p.
     * We take the first sign_od_p
     */
    if (w_sign * sign_od_p *sign_od_p0 > 0)
      in_cog = false;

    /* Last vertex for the face,
     * in_cog is true for all e_out edges which means that
     * w_sign <= 0 for all
     * AND no triangle was selected (because the line is passing exactly
     * in the cog G) */
    if ((i == (n_vertices - 1)) && in_cog) {
      u_sign = 1;
      v_sign = 1;
    }

    /* The projection of point O along displacement (OD) is outside of the
     * triangle then no intersection */
    if (w_sign > 0 || u_sign  < 0 || v_sign < 0)
      continue;

    /* Line (OD) intersects the triangle because
     * u_sign >= 0, v_sign >= 0 and w_sign <= 0
     */
    /* No need to deal with the special case where the intersection is in G */
    in_cog = false;

    double og_p = - cs_math_3_dot_product(vgo, pvec);

    /* This sign is absolute (ie same result is obtained if the face is seen from
     * the other adjacent cell.
     */
    int sign_og_p = (og_p > 0 ? 1 : -1);

    /* Same sign (meaning there is a possible intersection with t > 0). */
    if (sign_od_p == sign_og_p) {
      /* The line (OD) enters (n_inout[0]++)
       * (it means OD points toward cell i)
       * or leaves (n_inout[1]++) the cell
       * (it means OD points toward cell j)
       * */

      /* We don't care if it is entering or outgoing
       * but we store this information correctly
       * Convention
       * t > 0: entering i
       * t < 0: outgoing i */
      if (orient == 0) {
        /* entering cell i */
        if (sign_od_p == 1)
          n_inout[0]++;
        /* going out cell i */
        if (sign_od_p == -1)
          n_inout[1]++;

        if (fabs(og_p) < fabs(od_p)) {
          /* There is a real intersection (inward or outward) with 0 <= t < 1 */
          double t = 0.;
          n_intersects++;

          const double det = cs_math_3_norm(e0)*cs_math_3_norm(pvec);
          if (fabs(od_p) > epsilon * fabs(det)) {
            t = og_p / od_p;
          }

          if (t < retval) {
            retval = t;
            /* Store the normal if needed */
            if (face_norm != NULL)
              cs_math_3_normalise(pvec, face_norm);
          }
        }
      }
      else if (orient != sign_od_p) {
        n_inout[1]++;
        if (fabs(og_p) < fabs(od_p)) {
          /* There is a real intersection (outward) with 0 <= t < 1 */
          double t = 0.;
          n_intersects++;

          const double det = cs_math_3_norm(e0)*cs_math_3_norm(pvec);
          if (fabs(od_p) > epsilon * fabs(det)) {
            t = og_p / od_p;
          }

          if (t < retval) {
            retval = t;
            /* Store the normal if needed */
            if (face_norm != NULL)
              cs_math_3_normalise(pvec, face_norm);
          }
        }
      } else {
        n_inout[0]++;
        /* Incoming intersection on segment [OD] */
        if (fabs(og_p) < fabs(od_p))
          n_intersects--;
      }

    } else {
      /* Opposite sign (meaning there is a possible intersection of the line
       * with t<0).  */

      /* We don't care if it is entering or outgoing
       * but we store this information correctly
       * Convention
       * t > 0: entering i
       * t < 0: outgoing i */
      if (orient == 0) {
        cs_real_t t = -1;
        const double det = cs_math_3_norm(e0)*cs_math_3_norm(pvec);
        if (fabs(od_p) > epsilon * fabs(det))
          t = og_p / od_p;

        if (t < retval)
          retval = t;

        if (sign_od_p == -1)
          n_inout[1]++;
        else
          n_inout[0]++;

      }
      else if (orient != sign_od_p)
        n_inout[1]++;
      else
        n_inout[0]++;
    }

  }
  /* In case intersections were removed due to non-convex cases
   *  (i.e.  n_intersects < 1 , but retval < 1),
   *  the retval value is forced to 2
   *  (no intersection since the particle entered and left from this face). */
  if ((n_intersects < 1) && retval < 1. && retval >= 0) {
    retval = 2.;
  }

  return retval;
}

/*---------------------------------------------------------------------------*/

END_C_DECLS
