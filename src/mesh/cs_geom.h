#ifndef __CS_GEOM_H__
#define __CS_GEOM_H__

/*============================================================================
 * Geometric utility functions.
 *===========================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
 *  Local headers
 *---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *===========================================================================*/

/*============================================================================
 * Type definition
 *===========================================================================*/

/*=============================================================================
 * Global variables
 *===========================================================================*/

/*=============================================================================
 * Public function prototypes
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
                      int              *rank_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Test if a line segment intersects a face.
 *
 *                               |
 *      x------------------------|--------x D: end coordinates
 *   O: start coordiantes        |
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
 * \param[out]     n_crossings    number sub_face crossings
 *                                 [0: in; 1: out]
 * \param[in, out] face_norm      local face unite noraml of the crossed sub
 *                                 triangle (if entering with something
 *                                 different from NULL)
 *
 * \return
 *   2 if the segment does not go through the face's plane, or minimum
 *   relative distance (in terms of segment length)
 *   of intersection point to face.
 */
/*----------------------------------------------------------------------------*/

double
cs_geom_segment_intersect_face(int              orient,
                               cs_lnum_t        n_vertices,
                               const cs_lnum_t  vertex_ids[],
                               const cs_real_t  vtx_coord[][3],
                               const cs_real_t  face_cog[3],
                               const cs_real_t  sx0[3],
                               const cs_real_t  sx1[3],
                               int              n_crossings[2],
                               cs_real_t        *face_norm);

/*---------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GEOM_H__ */
