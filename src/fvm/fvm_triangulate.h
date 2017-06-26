#ifndef __FVM_TRIANGULATE_H__
#define __FVM_TRIANGULATE_H__

/*============================================================================
 * Triangulation of a polygon
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*
 * Pointer to structure maintaining the state of the current triangulation;
 * the structure itself is private.
 */

typedef struct _fvm_triangulate_state_t fvm_triangulate_state_t;

/*
 * Describe how the resulting triangle connectivity is defined.
 */

typedef enum {

  FVM_TRIANGULATE_MESH_DEF,      /* Definition by mesh vertex numbers */
  FVM_TRIANGULATE_ELT_DEF        /* Definition by local (element) vertex
                                    position (1 to n) */

} fvm_triangulate_def_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a structure necessary to the polygon triangulation algorithm.
 *
 * parameters:
 *   n_vertices_max    <-- maximum expected number of vertices per polygon.
 *
 * returns:
 *   pointer to polygon triangulation state structure.
 *----------------------------------------------------------------------------*/

fvm_triangulate_state_t *
fvm_triangulate_state_create(const int  n_vertices_max);

/*----------------------------------------------------------------------------
 * Destroy a structure necessary to the polygon triangulation algorithm.
 *
 * parameters:
 *   this_state  <-> pointer to structure that should be destroyed.
 *
 * returns:
 *   NULL pointer.
 *----------------------------------------------------------------------------*/

fvm_triangulate_state_t *
fvm_triangulate_state_destroy(fvm_triangulate_state_t  *this_state);

/*----------------------------------------------------------------------------
 * Triangulate a polygonal face.
 *
 * For a polygon with n vertices, we should obtain a triangluation with
 * (n-2) triangles and (2n-3) edges. If the polygon_vertices argument
 * is NULL, 1, 2, ...,n local numbering is implied.
 *
 * parameters:
 *   dim               <-- spatial dimension (2 or 3).
 *   base              <-- base numbering (usually 0 or 1)
 *   n_vertices        <-- number of vertices defining the polygon.
 *   coords            <-- coordinates of the triangulation's vertices.
 *   parent_vertex_num <-- optional indirection to vertex coordinates
 *                         (base to n-base).
 *   polygon_vertices  <-- polygon connectivity; size: n_vertices or empty.
 *   mode              <-- triangles connectivity by vertex number or
 *                         polygon vertex index.
 *   triangle_vertices --> triangles connectivity;
 *                         size: (n_vertices - 2) * 3.
 *   state             <-> associated triangulation state structure.
 *
 * returns:
 *   number of resulting triangles.
 *----------------------------------------------------------------------------*/

int
fvm_triangulate_polygon(int                             dim,
                        int                             base,
                        int                             n_vertices,
                        const cs_coord_t                coords[],
                        const cs_lnum_t                 parent_vertex_num[],
                        const cs_lnum_t                 polygon_vertices[],
                        fvm_triangulate_def_t           mode,
                        cs_lnum_t                       triangle_vertices[],
                        fvm_triangulate_state_t  *const state);

/*----------------------------------------------------------------------------
 * Triangulate a quadrangle.
 *
 * A convex quadrangle is divided into two triangles along its shortest
 * diagonal. A non-convex quadrangle may only be divided along the diagonal
 * which lies inside the quadrangle.
 *
 * If the quadrangle_vertices argument is NULL, 1, 2, ...,n local numbering
 * is implied.
 *
 * parameters:
 *   dim                  <-- spatial dimension (2 or 3).
 *   base                 <-- base numbering (usually 0 or 1)
 *   coords               <-- coordinates of the triangulation's vertices.
 *   parent_vertex_num    <-- optional indirection to vertex coordinates
 *                            (base to n-base numbering).
 *   quadrangle_vertices  <-- polygon connectivity; size: n_vertices or empty.
 *   triangle_vertices    --> triangles connectivity; size: 2 * 3.
 *
 * returns:
 *   number of resulting triangles.
 *----------------------------------------------------------------------------*/

int
fvm_triangulate_quadrangle(int               dim,
                           int               base,
                           const cs_coord_t  coords[],
                           const cs_lnum_t   parent_vertex_num[],
                           const cs_lnum_t   quadrangle_vertices[],
                           cs_lnum_t         triangle_vertices[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_TRIANGULATE_H__ */
