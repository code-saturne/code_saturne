#ifndef __MY_PLE_POINT_LOCATION_H__
#define __MY_PLE_POINT_LOCATION_H__

/*============================================================================
 * Locate local points in a nodal representation associated with a mesh
 *============================================================================*/

/*
  This file is part of the "Parallel Location and Exchange" library,
  intended to provide mesh or particle-based code coupling services.

  Copyright (C) 2005-2018  EDF S.A.

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------*/

#include "ple_config.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "ple_defs.h"
#include "my_ple_mesh.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Query number of extents and compute extents of a mesh representation.
 *
 * The minimum required functionality for this function is to compute
 * whole mesh extents, but it could also return extents of individual
 * elements, or intermediate extents of mesh subdivisions or coarsened
 * elements. As such, it takes an argument indicating the maximum
 * local number of extents it should compute (based on the size of
 * the extents array argument), but returns the number of extents
 * really computed, which may be lower (usually 1 for mesh extents,
 * possibly even 0 if the local mesh is empty). If n_max_extents = 1,
 * the whole mesh extents should be computed.
 *
 * If n_max_extents is set to a negative value (-1), no extents are computed,
 * but the function returns the maximum number of extents it may compute.
 * This query mode allows for the caller to allocate the correct amount
 * of memory for a subsequent call.
 *
 * parameters:
 *   mesh          <-- pointer to mesh representation structure
 *   n_max_extents <-- maximum number of sub-extents (such as element extents)
 *                     to compute, or -1 to query
 *   tolerance     <-- addition to local extents of each element:
 *                     extent = base_extent * (1 + tolerance)
 *   extents       <-> extents associated with the mesh or elements (or even
 *                     aggregated elements in case of coarser representation):
 *                     x_min_0, y_min_0, ..., x_max_i, y_max_i, ...
 *                     (size: 2*dim*n_max_extents), ignored in query mode
 * returns:
 *   the number of extents computed
 *----------------------------------------------------------------------------*/

ple_lnum_t
my_ple_mesh_extents(const void  *mesh,
                    ple_lnum_t   n_max_extents,
                    double       tolerance,
                    double       extents[]);

/*----------------------------------------------------------------------------
 * Find elements in a given nodal mesh containing points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this mesh, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   mesh              <-- pointer to mesh representation structure
 *   tolerance         <-- associated tolerance
 *   locate_on_parents <-- location relative to parent element numbers if
 *                         true, id of element + 1 in concatenated sections
 *                         of same element dimension if false
 *   n_points          <-- number of points to locate
 *   point_coords      <-- point coordinates
 *   point_tag         <-- optional point tag (size: n_points, ignored here)
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated, 0 - 1 if inside,
 *                         and > 1 if outside a volume element, or absolute
 *                         distance to a surface element (size: n_points)
 *----------------------------------------------------------------------------*/

void
my_ple_point_location_contain(const void          *mesh,
                              double               tolerance,
                              _Bool                locate_on_parents,
                              ple_lnum_t           n_points,
                              const ple_coord_t    point_coords[],
                              const int            point_tag[],
                              ple_lnum_t           location[],
                              float                distance[]);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __MY_PLE_POINT_LOCATION_H__ */
