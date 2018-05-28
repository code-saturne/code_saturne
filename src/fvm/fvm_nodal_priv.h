#ifndef __FVM_NODAL_PRIV_H__
#define __FVM_NODAL_PRIV_H__

/*============================================================================
 * Main structure for a nodal representation associated with a mesh
 *============================================================================*/

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_mesh.h"
#include "fvm_defs.h"
#include "fvm_group.h"
#include "fvm_nodal.h"
#include "fvm_tesselation.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure defining a mesh section
 *----------------------------------------------------------------------------*/

typedef struct _fvm_nodal_section_t {

  /* Basic information */
  /*-------------------*/

  int         entity_dim;          /* Entity dimension */

  cs_lnum_t   n_elements;          /* Number of elements */

  fvm_element_t  type;             /* Element types */

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

  cs_lnum_t   n_faces;             /* Number of faces defining polyhedra */

  /* Pointers to connectivity arrays, which may be shared */

  const cs_lnum_t   *face_index;   /* polyhedron -> faces index (O to n-1);
                                      size: n_elements + 1 */
  const cs_lnum_t   *face_num;     /* polyhedron -> face numbers (1 to n, signed,
                                      > 0 for outwards pointing face normal
                                      < 0 for inwards pointing face normal);
                                      size: face_index[n_elements] */

  const cs_lnum_t   *vertex_index; /* polygon face -> vertices index (O to n-1);
                                      size: n_faces + 1 */

  const cs_lnum_t   *vertex_num;   /* vertex numbers (1 to n);
                                      size: connectivity_size */

  /* Pointers to local connectivity arrays, if owner */

  cs_lnum_t   *_face_index;        /* face_index if owner, NULL if shared */
  cs_lnum_t   *_face_num;          /* face_num if owner, NULL if shared */
  cs_lnum_t   *_vertex_index;      /* vertex_index if owner, NULL if shared */
  cs_lnum_t   *_vertex_num;        /* vertex numbers if owner, NULL if shared */

  /* Pointers to group class ids, if present */

  int         *gc_id;              /* Group class id, NULL if implicit 0 */

  /* Optional tagging */

  int         *tag;                /* Element tag */

  /* Auxiliary structure used to define subdivision of elements into
     simpler element types (usually polygons to triangles and
     polyhedra to tetrahedra and pyramids) */

  fvm_tesselation_t  *tesselation;

  /* Numbering */
  /*-----------*/

  const cs_lnum_t   *parent_element_num; /* Local numbers (1 to n) of local
                                            elements in the parent mesh,
                                            associated with the section's
                                            elements.

                                            This array is necessary to redis-
                                            tribute output fields when the
                                            section has been either associated
                                            with an unsorted mixed mesh,
                                            renumbered, or is associated with a
                                            subset of a more complete mesh,
                                            such as a clip plane. When used for
                                            a subset, it also defines the lists
                                            of elements of the parent mesh
                                            belonging to that subset.

                                            This array is present only when non
                                            "trivial" (i.e. not 1, 2, ..., n). */

  cs_lnum_t     *_parent_element_num;    /* pointer to parent_element_num if
                                            owner, NULL otherwise */

  fvm_io_num_t  *global_element_num;     /* Global element numbers */

} fvm_nodal_section_t;

/*----------------------------------------------------------------------------
 * Structure defining a mesh in nodal definition
 *----------------------------------------------------------------------------*/

struct _fvm_nodal_t {

  /* Global indicators */
  /*-------------------*/

  char  *name;                 /* Mesh name */

  int    dim;                  /* Spatial dimension */
  int    num_dom;              /* Local domain number */
  int    n_doms;               /* Global number of domains */
  int    n_sections;           /* Number of sections */

  /* Local dimensions */
  /*------------------*/

  /* Total number of cells, faces, edges, and vertices */
  cs_lnum_t   n_cells;
  cs_lnum_t   n_faces;
  cs_lnum_t   n_edges;
  cs_lnum_t   n_vertices;

  /* Vertex definitions; */
  /*---------------------*/

  const cs_coord_t  *vertex_coords;     /* pointer to  vertex coordinates
                                           (always interlaced:
                                           x1, y1, z1, x2, y2, z2, ...) */
  cs_coord_t        *_vertex_coords;    /* pointer to vertex coordinates if
                                           owner (for use with own algorithms) */

  const cs_lnum_t   *parent_vertex_num; /* Local numbers (1 to n) of local
                                           vertices in the parent mesh.

                                           This array is necessary to redis-
                                           tribute output fields when a nodal
                                           mesh has been renumbered or is
                                           associated with a subset of a more
                                           complete mesh, such as a clip plane
                                           (in which case it also defines the
                                           lists of vertices of the parent
                                           mesh in that subset).

                                           This array is present only when non
                                           "trivial" (i.e. not 1, 2, ..., n). */

  cs_lnum_t     *_parent_vertex_num;    /* pointer to parent_vertex_num if
                                           owner, NULL otherwise */

  fvm_io_num_t  *global_vertex_num;     /* Global vertex numbering */

  /* Mesh connectivity */
  /*-------------------*/

  fvm_nodal_section_t  **sections;  /* Array of section descriptions */

  /* Metadata */
  /*----------*/

  /* Group class descriptions if present */

  fvm_group_class_set_t   *gc_set;  /* Pointer to group class set, or NULL */

  /* Optional global vertex labels.
     As these are expected to be used only for small sets (i.e. probes)
     where the point set is built from a global definition and data movement
     would adds complexity and overhead, the labels refer to a global view
     on rank 0; for the same reason, only shared labels are needed */

  char   **global_vertex_labels;  /* Pointer to vertex labels, or NULL */

  /* Pointer to parent mesh, if defined */

  const cs_mesh_t  *parent;

};

/*=============================================================================
 * Semi-private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation of a nodal mesh section representation structure.
 *
 * parameters:
 *   type <-- type of element defined by this section
 *
 * returns:
 *  pointer to created nodal mesh section representation structure
 *----------------------------------------------------------------------------*/

fvm_nodal_section_t *
fvm_nodal_section_create(const fvm_element_t  type);

/*----------------------------------------------------------------------------
 * Destruction of a nodal mesh section representation structure.
 *
 * parameters:
 *   this_section <-> pointer to structure that should be destroyed
 *
 * returns:
 *  NULL pointer
 *----------------------------------------------------------------------------*/

fvm_nodal_section_t *
fvm_nodal_section_destroy(fvm_nodal_section_t  * this_section);

/*----------------------------------------------------------------------------
 * Copy selected shared connectivity information to private connectivity
 * for a nodal mesh section .
 *
 * parameters:
 *   this_section      <-> pointer to section structure
 *   copy_face_index   <-- copy face index (polyhedra only) ?
 *   copy_face_num     <-- copy face numbers (polyhedra only) ?
 *   copy_vertex_index <-- copy vertex index (polyhedra/polygons only) ?
 *   copy_vertex_num   <-- copy vertex numbers ?
 *----------------------------------------------------------------------------*/

void
fvm_nodal_section_copy_on_write(fvm_nodal_section_t  *this_section,
                                bool                  copy_face_index,
                                bool                  copy_face_num,
                                bool                  copy_vertex_index,
                                bool                  copy_vertex_num);

/*----------------------------------------------------------------------------
 * Return global number of elements associated with section.
 *
 * parameters:
 *   this_section      <-- pointer to section structure
 *
 * returns:
 *   global number of elements associated with section
 *----------------------------------------------------------------------------*/

cs_gnum_t
fvm_nodal_section_n_g_elements(const fvm_nodal_section_t  *this_section);

/*----------------------------------------------------------------------------
 * Return global number of vertices associated with nodal mesh.
 *
 * parameters:
 *   this_nodal           <-- pointer to nodal mesh structure
 *
 * returns:
 *   global number of vertices associated with nodal mesh
 *----------------------------------------------------------------------------*/

cs_gnum_t
fvm_nodal_n_g_vertices(const fvm_nodal_t  *this_nodal);

/*----------------------------------------------------------------------------
 * Define cell->face connectivity for strided cell types.
 *
 * parameters:
 *   element_type     <-- type of strided element
 *   n_faces          --> number of element faces
 *   n_face_vertices  --> number of vertices of each face
 *   face_vertices    --> face -> vertex base connectivity (0 to n-1)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_cell_face_connect(fvm_element_t   element_type,
                            int            *n_faces,
                            int             n_face_vertices[6],
                            int             face_vertices[6][4]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_NODAL_PRIV_H__ */
