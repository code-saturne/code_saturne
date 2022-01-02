#ifndef __MY_PLE_MESH_PRIV_H__
#define __MY_PLE_MESH_PRIV_H__

/*============================================================================
 * Main structure for a nodal representation associated with a mesh
 *============================================================================*/

/*
  This file is part of the "Parallel Location and Exchange" library,
  intended to provide mesh or particle-based code coupling services.

  Copyright (C) 2005-2022  EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "ple_defs.h"

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

/*----------------------------------------------------------------------------
 * Structure defining a mesh section
 *----------------------------------------------------------------------------*/

typedef struct _my_ple_mesh_section_t {

  /* Basic information */
  /*-------------------*/

  int         entity_dim;          /* Entity dimension */

  ple_lnum_t  n_elements;          /* Number of elements */

  my_ple_element_t  type;         /* Element types */

  /* Connectivity */
  /*--------------*/

  size_t      connectivity_size;   /* Size of vertex_num array;
                                      for strided elements:
                                       (n_elements * stride)
                                      for polygons:
                                       (vertex_index[n_elements])
                                      for polyhedra:
                                       (vertex_index[n_faces]) */

  int         stride;              /* Element size for regular elements
                                      (0 for polygons and polyhedra) */

  ple_lnum_t  n_faces;             /* Number of faces defining polyhedra */

  /* Pointers to connectivity arrays, which may be shared */

  const ple_lnum_t  *face_index;   /* polyhedron -> faces index (O to n-1);
                                      size: n_elements + 1 */
  const ple_lnum_t  *face_num;     /* polyhedron -> face numbers (1 to n, signed,
                                      > 0 for outwards pointing face normal
                                      < 0 for inwards pointing face normal);
                                      size: face_index[n_elements] */

  const ple_lnum_t  *vertex_index; /* polygon face -> vertices index (O to n-1);
                                      size: n_faces + 1 */

  const ple_lnum_t  *vertex_num;   /* vertex numbers (1 to n);
                                      size: connectivity_size */

  /* Pointers to local connectivity arrays, if owner */

  ple_lnum_t  *_face_index;        /* face_index if owner, NULL if shared */
  ple_lnum_t  *_face_num;          /* face_num if owner, NULL if shared */
  ple_lnum_t  *_vertex_index;      /* vertex_index if owner, NULL if shared */
  ple_lnum_t  *_vertex_num;        /* vertex numbers if owner, NULL if shared */

  /* Pointers to element numbering */

  ple_lnum_t  *element_num;        /* element numbers, or NULL */

} my_ple_mesh_section_t;

/*----------------------------------------------------------------------------
 * Structure defining a mesh in nodal definition
 *----------------------------------------------------------------------------*/

struct _my_ple_mesh_t {

  /* Global indicators */
  /*-------------------*/

  int    dim;                  /* Spatial dimension */
  int    n_sections;           /* Number of sections */

  /* Local dimensions */
  /*------------------*/

  /* Total number of cells, faces, edges, and vertices */
  ple_lnum_t  n_cells;
  ple_lnum_t  n_faces;
  ple_lnum_t  n_edges;
  ple_lnum_t  n_vertices;

  /* Vertex definitions; */
  /*---------------------*/

  const ple_coord_t  *vertex_coords;    /* pointer to  vertex coordinates
                                           (always interlaced:
                                           x1, y1, z1, x2, y2, z2, ...) */
  ple_coord_t        *_vertex_coords;   /* pointer to vertex coordinates if
                                           owner (for use with own algorithms) */

  /* Mesh connectivity */
  /*-------------------*/

  my_ple_mesh_section_t  **sections;  /* Array of section descriptions */

};

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __MY_PLE_MESH_PRIV_H__ */
