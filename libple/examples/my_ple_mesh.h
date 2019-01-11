#ifndef __MY_PLE_MESH_H__
#define __MY_PLE_MESH_H__

/*============================================================================
 * Main structure for a nodal representation associated with a mesh
 *============================================================================*/

/*

  This file is part of the "Parallel Location and Exchange" library,
  intended to provide mesh or particle-based code coupling services.

  Copyright (C) 2005-2019  EDF S.A.

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
 * Element types
 *----------------------------------------------------------------------------*/

typedef enum {

  MY_PLE_EDGE,             /* Edge */
  MY_PLE_FACE_TRIA,        /* Triangle */
  MY_PLE_FACE_QUAD,        /* Quadrangle */
  MY_PLE_FACE_POLY,        /* Simple Polygon */
  MY_PLE_CELL_TETRA,       /* Tetrahedron */
  MY_PLE_CELL_PYRAM,       /* Pyramid */
  MY_PLE_CELL_PRISM,       /* Prism (pentahedron) */
  MY_PLE_CELL_HEXA,        /* Hexahedron (brick) */
  MY_PLE_CELL_POLY,        /* Simple Polyhedron (convex or quasi-convex) */
  MY_PLE_N_ELEMENT_TYPES   /* Number of element types */

} my_ple_element_t;

/*----------------------------------------------------------------------------
 * Structure defining a mesh in nodal definition
 *----------------------------------------------------------------------------*/

typedef struct _my_ple_mesh_t my_ple_mesh_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/* Names of (single) element types */

extern const char  *my_ple_element_type_name[];

/* Number of vertices associated with each "nodal" element type */

extern const int  my_ple_mesh_n_vertices_element[];

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation of a mesh representation structure.
 *
 * parameters:
 *   dim  <-- spatial dimension
 *
 * returns:
 *  pointer to created mesh representation structure
 *----------------------------------------------------------------------------*/

my_ple_mesh_t *
my_ple_mesh_create(int  dim);

/*----------------------------------------------------------------------------
 * Destruction of a mesh representation structure.
 *
 * parameters:
 *   mesh  <-> pointer to pointer to structure that should be destroyed
 *----------------------------------------------------------------------------*/

void
my_ple_mesh_destroy(my_ple_mesh_t  **mesh);

/*----------------------------------------------------------------------------
 * Assign shared vertex coordinates to an extracted mesh.
 *
 * parameters:
 *   mesh          <-> mesh structure
 *   vertex_coords <-- vertex coordinates (interlaced)
 *----------------------------------------------------------------------------*/

void
my_ple_mesh_set_shared_vertices(my_ple_mesh_t     *mesh,
                                 const ple_coord_t   vertex_coords[]);

/*----------------------------------------------------------------------------
 * Assign private vertex coordinates to a mesh.
 *
 * Ownership of the given coordinates array is transferred to
 * the mesh representation structure.
 *
 * parameters:
 *   mesh          <-> mesh structure
 *   vertex_coords <-- vertex coordinates (interlaced)
 *----------------------------------------------------------------------------*/

void
my_ple_mesh_transfer_vertices(my_ple_mesh_t  *mesh,
                               ple_coord_t      vertex_coords[]);

/*----------------------------------------------------------------------------
 * Append a new section to an existing my_ple_mesh mesh, and transfer
 * ownership of the given connectivity to that section.
 *
 * parameters:
 *   mesh         <-> mesh structure
 *   n_elements   <-- number of elements to add
 *   type         <-- type of elements to add
 *   face_index   <-- polyhedron -> faces index (O to n-1)
 *                    size: n_elements + 1
 *   face_num     <-- polyhedron -> face numbers (1 to n, signed,
 *                    > 0 for outwards pointing face normal
 *                    < 0 for inwards pointing face normal);
 *                    size: face_index[n_elements]
 *   vertex_index <-- polygon face -> vertices index (O to n-1)
 *                    size: face_index[n_elements]
 *   vertex_num   <-- element -> vertex connectivity
 *   element_num  <-- element numbers
 *----------------------------------------------------------------------------*/

void
my_ple_mesh_append_by_transfer(my_ple_mesh_t     *mesh,
                               ple_lnum_t         n_elements,
                               my_ple_element_t   type,
                               ple_lnum_t         face_index[],
                               ple_lnum_t         face_num[],
                               ple_lnum_t         vertex_index[],
                               ple_lnum_t         vertex_num[],
                               ple_lnum_t         elt_num[]);

/*----------------------------------------------------------------------------
 * Append a new section to an existing my_ple_mesh mesh, sharing the given
 * given connectivity with the caller.
 *
 * The caller should not destroy or modify the arrays passed to this
 * function until the mesh is destroyed.
 *
 * parameters:
 *   mesh         <-> mesh structure
 *   n_elements   <-- number of elements to add
 *   type         <-- type of elements to add
 *   face_index   <-- polyhedron -> faces index (O to n-1)
 *                    size: n_elements + 1
 *   face_num     <-- polyhedron -> face numbers (1 to n, signed,
 *                    > 0 for outwards pointing face normal
 *                    < 0 for inwards pointing face normal);
 *                    size: face_index[n_elements]
 *   vertex_index <-- polygon face -> vertices index (O to n-1)
 *                    size: face_index[n_elements]
 *   vertex_num   <-- element -> vertex connectivity
 *   element_num  <-- element numbers
 *----------------------------------------------------------------------------*/

void
my_ple_mesh_append_shared(my_ple_mesh_t     *mesh,
                          ple_lnum_t         n_elements,
                          my_ple_element_t   type,
                          ple_lnum_t         face_index[],
                          ple_lnum_t         face_num[],
                          ple_lnum_t         vertex_index[],
                          ple_lnum_t         vertex_num[],
                          ple_lnum_t         element_num[]);

/*----------------------------------------------------------------------------
 * Return maximum dimension of entities in a mesh.
 *
 * parameters:
 *   mesh <-- pointer to mesh structure
 *
 * returns:
 *  maximum dimension of entities in mesh (0 to 3)
 *----------------------------------------------------------------------------*/

int
my_ple_mesh_get_max_entity_dim(const my_ple_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Dump printout of a nodal representation structure.
 *
 * parameters:
 *   mesh <-- pointer to structure that should be dumped
 *----------------------------------------------------------------------------*/

void
my_ple_mesh_dump(const my_ple_mesh_t  *mesh);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __MY_PLE_MESH_H__ */
