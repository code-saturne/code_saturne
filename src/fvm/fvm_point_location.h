#ifndef __FVM_POINT_LOCATION_H__
#define __FVM_POINT_LOCATION_H__

/*============================================================================
 * Locate local points in a nodal representation associated with a mesh
 *============================================================================*/

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"
#include "fvm_nodal.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

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
 * Find elements in a given nodal mesh containing points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this mesh, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh representation structure
 *   tolerance_base       <-- associated base tolerance (used for bounding
 *                            box check only, not for location test)
 *   tolerance_multiplier <-- associated fraction of element bounding boxes
 *                            added to tolerance
 *   locate_on_parents    <-- location relative to parent element numbers if 1,
 *                            id of element + 1 in concatenated sections of
 *                            same element dimension if 0
 *   n_points             <-- number of points to locate
 *   point_tag            <-- optional point tag
 *   point_coords         <-- point coordinates
 *   location             <-> number of element containing or closest to each
 *                            point (size: n_points)
 *   distance             <-> distance from point to element indicated by
 *                            location[]: < 0 if unlocated, 0 - 1 if inside,
 *                            and > 1 if outside a volume element, or absolute
 *                            distance to a surface element (size: n_points)
 *----------------------------------------------------------------------------*/

void
fvm_point_location_nodal(const fvm_nodal_t  *this_nodal,
                         float               tolerance_base,
                         float               tolerance_fraction,
                         int                 locate_on_parents,
                         cs_lnum_t           n_points,
                         const int          *point_tag,
                         const cs_coord_t    point_coords[],
                         cs_lnum_t           location[],
                         float               distance[]);

/*----------------------------------------------------------------------------
 * For each point previously located in a element, find among vertices of this
 * element the closest vertex relative to this point.
 *
 * As input, located_ent_num is an array with a numbering not using a parent
 * numbering. As output, located_ent_num may use a parent numbering
 * according to the value of locate_on_parents.
 *
 * The located_vtx_num output is also determined relative to the
 * locate_on_parents option.
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh representation structure
 *   locate_on_parents    <-- location relative to parent element numbers if 1,
 *                            id of element + 1 in concatenated sections of
 *                            same element dimension if 0
 *   n_points             <-- number of points to locate
 *   point_coords         <-- point coordinates
 *   located_ent_num      <-> input: list of elements (cells or faces according
 *                            to max entity dim) where points have been
 *                            initially located or not (size: n_points)
 *                            output: possibly modified by parent numbering
 *   located_vtx_num      <-> output: list of vertices closest to each point
 *----------------------------------------------------------------------------*/

void
fvm_point_location_closest_vertex(const fvm_nodal_t  *this_nodal,
                                  int                 locate_on_parents,
                                  cs_lnum_t           n_points,
                                  const cs_coord_t    point_coords[],
                                  cs_lnum_t           located_ent_num[],
                                  cs_lnum_t           located_vtx_num[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_POINT_LOCATION_H__ */
