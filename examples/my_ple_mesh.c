/*============================================================================
 * Main structure for a nodal representation associated with a mesh
 *============================================================================*/

/*
  This file is part of the "Parallel Location and Exchange" library,
  intended to provide mesh or particle-based code coupling services.

  Copyright (C) 2005-2024  EDF S.A.

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "ple_config_defs.h"
#include "ple_defs.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "my_ple_mesh.h"
#include "my_ple_mesh_priv.h"

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

#define _ABS(a)     ((a) <  0  ? -(a) : (a))  /* Absolute value of a */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Number of vertices associated with each "nodal" element type */

const int  my_ple_mesh_n_vertices_element[] = {2,   /* Edge */
                                               3,   /* Triangle */
                                               4,   /* Quadrangle */
                                               0,   /* Simple polygon */
                                               4,   /* Tetrahedron */
                                               5,   /* Pyramid */
                                               6,   /* Prism */
                                               8,   /* Hexahedron */
                                               0};  /* Simple polyhedron */

/* Names of (single) "nodal" element types */

const char  *my_ple_element_type_name[] = {N_("edge"),
                                           N_("triangle"),
                                           N_("quadrangle"),
                                           N_("simple polygon"),
                                           N_("tetrahedron"),
                                           N_("pyramid"),
                                           N_("prism"),
                                           N_("hexahedron"),
                                           N_("simple polyhedron")};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation of a mesh section representation structure.
 *
 * parameters:
 *   type <-- type of element defined by this section
 *
 * returns:
 *   pointer to created mesh section representation structure
 *----------------------------------------------------------------------------*/

static my_ple_mesh_section_t *
_nodal_section_create(const my_ple_element_t  type)
{
  my_ple_mesh_section_t  *this_section;

  PLE_MALLOC(this_section, 1, my_ple_mesh_section_t);

  /* Global information */

  if (type == MY_PLE_EDGE)
    this_section->entity_dim = 1;
  else if (type >= MY_PLE_FACE_TRIA && type <= MY_PLE_FACE_POLY)
    this_section->entity_dim = 2;
  else
    this_section->entity_dim = 3;

  this_section->n_elements = 0;
  this_section->type = type;

  /* Connectivity */

  this_section->connectivity_size = 0;

  if (type != MY_PLE_FACE_POLY && type != MY_PLE_CELL_POLY)
    this_section->stride = my_ple_mesh_n_vertices_element[type];
  else
    this_section->stride = 0;

  this_section->n_faces = 0;
  this_section->face_index = NULL;
  this_section->face_num   = NULL;
  this_section->vertex_index = NULL;
  this_section->vertex_num = NULL;

  this_section->_face_index = NULL;
  this_section->_face_num   = NULL;
  this_section->_vertex_index = NULL;
  this_section->_vertex_num = NULL;

  this_section->element_num = NULL;

  return (this_section);
}

/*----------------------------------------------------------------------------
 * Destruction of a mesh section representation structure.
 *
 * parameters:
 *   this_section <-> pointer to pointer to structure that should be destroyed
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

static void
_nodal_section_destroy(my_ple_mesh_section_t  **this_section)
{
  my_ple_mesh_section_t  *s = *this_section;

  /* Connectivity */

  if (s->_face_index != NULL)
    PLE_FREE(s->_face_index);
  if (s->_face_num != NULL)
    PLE_FREE(s->_face_num);

  if (s->_vertex_index != NULL)
    PLE_FREE(s->_vertex_index);
  if (s->_vertex_num != NULL)
    PLE_FREE(s->_vertex_num);

  if (s->element_num != NULL)
    PLE_FREE(s->element_num);

  /* Main structure destroyed */

  PLE_FREE(*this_section);
}

/*----------------------------------------------------------------------------
 * Create new section, transferring ownership of the given connectivity
 * and optional element number arrays to that section.
 *
 * parameters:
 *   n_elements         <-- number of elements in section
 *   type               <-- type of elements to add
 *   face_index         <-- polyhedron -> faces index (O to n-1)
 *                          size: n_elements + 1
 *   face_num           <-- polyhedron -> face numbers (1 to n, signed,
 *                          > 0 for outwards pointing face normal
 *                          < 0 for inwards pointing face normal);
 *                          size: face_index[n_elements]
 *   vertex_index       <-- polygon face -> vertices index (O to n-1)
 *                          size: face_index[n_elements]
 *   vertex_num         <-- element -> vertex connectivity
 *   element_num        <-- element number, or NULL
 *----------------------------------------------------------------------------*/

static my_ple_mesh_section_t *
_transfer_to_section(ple_lnum_t        n_elements,
                     my_ple_element_t  type,
                     ple_lnum_t        face_index[],
                     ple_lnum_t        face_num[],
                     ple_lnum_t        vertex_index[],
                     ple_lnum_t        vertex_num[],
                     ple_lnum_t        element_num[])
{
  my_ple_mesh_section_t  *this_section = NULL;

  this_section = _nodal_section_create(type);

  this_section->n_elements = n_elements;

  /* Connectivity */

  if (type == MY_PLE_CELL_POLY) {
    this_section->_face_index = face_index;
    this_section->_face_num = face_num;
  }

  if (type == MY_PLE_FACE_POLY || type == MY_PLE_CELL_POLY)
    this_section->_vertex_index = vertex_index;

  this_section->_vertex_num = vertex_num;

  /* Shared arrays */

  this_section->face_index = this_section->_face_index;
  this_section->face_num = this_section->_face_num;
  this_section->vertex_index = this_section->_vertex_index;
  this_section->vertex_num = this_section->_vertex_num;

  /* Connectivity size */

  if (this_section->stride != 0)
    this_section->connectivity_size
      = this_section->n_elements * this_section->stride;

  else if (this_section->type == MY_PLE_FACE_POLY)
    this_section->connectivity_size
      = this_section->vertex_index[this_section->n_elements];

  else if (this_section->type == MY_PLE_CELL_POLY) {
    ple_lnum_t i, _face_num;
    for (i = 0;
         i < this_section->face_index[this_section->n_elements];
         i++) {
      _face_num = _ABS(this_section->face_num[i]);
      if (_face_num > this_section->n_faces)
        this_section->n_faces = _face_num;
    }
    this_section->connectivity_size
      = this_section->vertex_index[this_section->n_faces];
  }

  this_section->element_num = element_num;

  return this_section;
}

/*----------------------------------------------------------------------------
 * Create new section, mapping the given connectivity to that section.
 *
 * The optional element element number is always transferred to the section.
 *
 * parameters:
 *   n_elements   <-- number of elements in section
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
 *   element_num  <-- element number, or NULL
 *----------------------------------------------------------------------------*/

static my_ple_mesh_section_t *
_map_to_section(ple_lnum_t         n_elements,
                my_ple_element_t   type,
                const ple_lnum_t   face_index[],
                const ple_lnum_t   face_num[],
                const ple_lnum_t   vertex_index[],
                ple_lnum_t         vertex_num[],
                ple_lnum_t         element_num[])
{
  my_ple_mesh_section_t  *this_section = NULL;

  this_section = _nodal_section_create(type);

  this_section->n_elements = n_elements;

  /* Connectivity */

  if (type == MY_PLE_CELL_POLY) {
    this_section->face_index = face_index;
    this_section->face_num = face_num;
  }

  if (type == MY_PLE_FACE_POLY || type == MY_PLE_CELL_POLY)
    this_section->vertex_index = vertex_index;

  this_section->vertex_num = vertex_num;

  /* Connectivity size */

  if (this_section->stride != 0)
    this_section->connectivity_size
      = this_section->n_elements * this_section->stride;

  else if (this_section->type == MY_PLE_FACE_POLY)
    this_section->connectivity_size
      = this_section->vertex_index[this_section->n_elements];

  else if (this_section->type == MY_PLE_CELL_POLY) {
    ple_lnum_t i, _face_num;
    for (i = 0;
         i < this_section->face_index[this_section->n_elements];
         i++) {
      _face_num = _ABS(this_section->face_num[i]);
      if (_face_num > this_section->n_faces)
        this_section->n_faces = _face_num;
    }
    this_section->connectivity_size
      = this_section->vertex_index[this_section->n_faces];
  }

  this_section->element_num = element_num;

  return this_section;
}

/*----------------------------------------------------------------------------
 * Dump printout of a nodal representation structure section.
 *
 * parameters:
 *   this_section <-- pointer to structure that should be dumped
 *----------------------------------------------------------------------------*/

static void
_my_ple_mesh_section_dump(const my_ple_mesh_section_t  *this_section)
{
  ple_lnum_t  n_elements, i, j;
  const ple_lnum_t  *idx, *num;

  /* Global indicators */
  /*--------------------*/

  ple_printf("\n"
             "Entity dimension:     %d\n"
             "Number of elements:   %ld\n"
             "Element type:         %s\n",
             this_section->entity_dim, (long)this_section->n_elements,
             my_ple_element_type_name[this_section->type]);

  ple_printf("\n"
             "Connectivity_size:     %lu\n"
             "Stride:                %d\n"
             "Number of faces:       %d\n",
             (unsigned long)(this_section->connectivity_size),
             this_section->stride,
             (long)(this_section->n_faces));

  ple_printf("\n"
             "Pointers to shareable arrays:\n"
             "  face_index:           %p\n"
             "  face_num:             %p\n"
             "  vertex_index:         %p\n"
             "  vertex_num:           %p\n",
             this_section->face_index, this_section->face_num,
             this_section->vertex_index, this_section->vertex_num);

  ple_printf("\n"
             "Pointers to local arrays:\n"
             "  _face_index:          %p\n"
             "  _face_num:            %p\n"
             "  _vertex_index:        %p\n"
             "  _vertex_num:          %p\n",
             this_section->_face_index, this_section->_face_num,
             this_section->_vertex_index, this_section->_vertex_num);

  if (this_section->face_index != NULL) {
    ple_printf("\nPolyhedra -> faces (polygons) connectivity:\n\n");
    n_elements = this_section->n_elements;
    idx = this_section->face_index;
    num = this_section->face_num;
    for (i = 0; i < n_elements; i++) {
      ple_printf("%10d (idx = %10d) %10d\n",
                 i+1, idx[i], num[idx[i]]);
      for (j = idx[i] + 1; j < idx[i + 1]; j++)
        ple_printf("                              %10d\n", num[j]);
    }
    ple_printf("      end  (idx = %10d)\n", idx[n_elements]);
  }

  if (this_section->vertex_index != NULL) {
    ple_lnum_t  n_faces = (this_section->n_faces > 0) ?
                          this_section->n_faces : this_section->n_elements;
    ple_printf("\nPolygons -> vertices connectivity:\n\n");
    idx = this_section->vertex_index;
    num = this_section->vertex_num;
    for (i = 0; i < n_faces; i++) {
      ple_printf("%10d (idx = %10d) %10d\n",
                i + 1, idx[i], num[idx[i]]);
      for (j = idx[i] + 1; j < idx[i + 1]; j++)
        ple_printf("                              %10d\n", num[j]);
    }
    ple_printf("      end  (idx = %10d)\n", idx[n_faces]);
  }

  else {
    ple_printf("\nElements -> vertices connectivity:\n\n");
    n_elements = this_section->n_elements;
    num = this_section->vertex_num;
    switch(this_section->stride) {
    case 2:
      for (i = 0; i < n_elements; i++)
        ple_printf("%10d : %10d %10d\n",
                   i+1, num[i*2], num[i*2+1]);
      break;
    case 3:
      for (i = 0; i < n_elements; i++)
        ple_printf("%10d : %10d %10d %10d\n",
                   i+1, num[i*3], num[i*3+1], num[i*3+2]);
      break;
    case 4:
      for (i = 0; i < n_elements; i++)
        ple_printf("%10d : %10d %10d %10d %10d\n",
                   i+1, num[i*4], num[i*4+1], num[i*4+2],
                   num[i*4+3]);
      break;
    case 5:
      for (i = 0; i < n_elements; i++)
        ple_printf("%10d : %10d %10d %10d %10d %10d\n",
                   i+1, num[i*5], num[i*5+1], num[i*5+2],
                   num[i*5+3], num[i*5+4]);
      break;
    case 6:
      for (i = 0; i < n_elements; i++)
        ple_printf("%10d : %10d %10d %10d %10d %10d %10d\n",
                   i+1, num[i*6], num[i*6+1], num[i*6+2],
                   num[i*6+3], num[i*6+4], num[i*6+5]);
      break;
    case 8:
      for (i = 0; i < n_elements; i++)
        ple_printf("%10d : %10d %10d %10d %10d %10d %10d %10d %10d\n",
                   i+1, num[i*8], num[i*8+1], num[i*8+2], num[i*8+3],
                   num[i*8+4], num[i*8+5], num[i*8+6], num[i*8+7]);
      break;
    default:
      for (i = 0; i < n_elements; i++) {
        ple_printf("%10d :", i+1);
        for (j = 0; j < this_section->stride; j++)
          ple_printf(" %10d", num[i*this_section->stride + j]);
        ple_printf("\n");
      }
    }
  }

  /* Numbers of associated elements */

  ple_printf("\nLocal element numbers:\n");
  if (this_section->element_num == NULL)
    ple_printf("\n  Nil\n\n");
  else {
    for (i = 0; i < this_section->n_elements; i++)
      ple_printf("  %10d %10d\n", i+1, this_section->element_num[i]);
  }
}

/*============================================================================
 * Public function definitions
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
my_ple_mesh_create(int  dim)
{
  my_ple_mesh_t  *mesh;

  PLE_MALLOC(mesh, 1, my_ple_mesh_t);

  /* Global indicators */

  mesh->dim     = dim;
  mesh->n_sections = 0;

  /* Local dimensions */

  mesh->n_cells = 0;
  mesh->n_faces = 0;
  mesh->n_edges = 0;
  mesh->n_vertices = 0;

  /* Local structures */

  mesh->vertex_coords = NULL;
  mesh->_vertex_coords = NULL;

  mesh->sections = NULL;

  return (mesh);
}

/*----------------------------------------------------------------------------
 * Destruction of a mesh representation structure.
 *
 * parameters:
 *   mesh  <-> pointer to pointer to structure that should be destroyed
 *----------------------------------------------------------------------------*/

void
my_ple_mesh_destroy(my_ple_mesh_t  **mesh)
{
  int  i;
  my_ple_mesh_t  *m = *mesh;

  /* Local structures */

  if (m->_vertex_coords != NULL)
    PLE_FREE(m->_vertex_coords);

  for (i = 0; i < m->n_sections; i++)
    _nodal_section_destroy(&(m->sections[i]));

  if (m->sections != NULL)
    PLE_FREE(m->sections);

  /* Main structure destroyed */

  PLE_FREE(*mesh);
}

/*----------------------------------------------------------------------------
 * Assign shared vertex coordinates to an extracted mesh.
 *
 * parameters:
 *   mesh              <-> mesh structure
 *   vertex_coords     <-- vertex coordinates (interlaced)
 *----------------------------------------------------------------------------*/

void
my_ple_mesh_set_shared_vertices(my_ple_mesh_t      *mesh,
                                const ple_coord_t   vertex_coords[])
{
  assert(mesh != NULL);

  /* Map vertex coordinates to array passed as argument
     (mesh->_vertex_coords remains NULL, so only
     the const pointer may be used for a shared array) */

  mesh->vertex_coords = vertex_coords;
}

/*----------------------------------------------------------------------------
 * Assign private vertex coordinates to a mesh.
 *
 * Ownership of the given coordinates array is transferred to
 * the mesh representation structure.
 *
 * parameters:
 *   mesh              <-> mesh structure
 *   vertex_coords     <-- vertex coordinates (interlaced)
 *----------------------------------------------------------------------------*/

void
my_ple_mesh_transfer_vertices(my_ple_mesh_t  *mesh,
                              ple_coord_t     vertex_coords[])
{
  assert(mesh != NULL);

  mesh->_vertex_coords = vertex_coords;
  mesh->vertex_coords = vertex_coords;
}

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
                               ple_lnum_t         elt_num[])
{
  my_ple_mesh_section_t  *new_section = NULL;
  int  n_sections = 0;

  assert(mesh != NULL);

  n_sections = mesh->n_sections;

  /* Create new section */

  PLE_REALLOC(mesh->sections, n_sections + 1, my_ple_mesh_section_t *);

  new_section = _transfer_to_section(n_elements,
                                     type,
                                     face_index,
                                     face_num,
                                     vertex_index,
                                     vertex_num,
                                     elt_num);

  mesh->sections[n_sections] = new_section;
  mesh->n_sections += 1;

  /* Update main structure information */

  switch(new_section->entity_dim) {
  case 3:
    mesh->n_cells += n_elements;
    break;
  case 2:
    mesh->n_faces += n_elements;
    break;
  case 1:
    mesh->n_edges += n_elements;
    break;
  default:
    assert(0);
  }

}

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
my_ple_mesh_append_shared(my_ple_mesh_t      *mesh,
                           ple_lnum_t         n_elements,
                           my_ple_element_t   type,
                           ple_lnum_t         face_index[],
                           ple_lnum_t         face_num[],
                           ple_lnum_t         vertex_index[],
                           ple_lnum_t         vertex_num[],
                           ple_lnum_t         element_num[])
{
  my_ple_mesh_section_t  *new_section = NULL;
  int  n_sections = 0;

  assert(mesh != NULL);

  n_sections = mesh->n_sections;

  /* Create new section */

  PLE_REALLOC(mesh->sections, n_sections + 1, my_ple_mesh_section_t *);

  new_section = _map_to_section(n_elements,
                                type,
                                face_index,
                                face_num,
                                vertex_index,
                                vertex_num,
                                element_num);

  mesh->sections[n_sections] = new_section;
  mesh->n_sections += 1;

  /* Update main structure information */

  switch(new_section->entity_dim) {
  case 3:
    mesh->n_cells += n_elements;
    break;
  case 2:
    mesh->n_faces += n_elements;
    break;
  case 1:
    mesh->n_edges += n_elements;
    break;
  default:
    assert(0);
  }

}

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
my_ple_mesh_get_max_entity_dim(const my_ple_mesh_t  *mesh)
{
  int  section_id;
  int  max_entity_dim = 0;

  assert(mesh != NULL);

  for (section_id = 0; section_id < mesh->n_sections; section_id++) {
    const my_ple_mesh_section_t  *section = mesh->sections[section_id];
    if (section->entity_dim > max_entity_dim)
      max_entity_dim = section->entity_dim;
  }

  return max_entity_dim;
}

/*----------------------------------------------------------------------------
 * Dump printout of a nodal representation structure.
 *
 * parameters:
 *   mesh <-- pointer to structure that should be dumped
 *----------------------------------------------------------------------------*/

void
my_ple_mesh_dump(const my_ple_mesh_t  *mesh)
{
  ple_lnum_t  i;
  ple_lnum_t  num_vertex = 1;
  const ple_coord_t  *coord = mesh->vertex_coords;

  /* Global indicators */
  /*--------------------*/

  ple_printf("\n"
             "Mesh dimension:               %d\n"
             "Number of sections:           %d\n",
             mesh->dim, mesh->n_sections);

  ple_printf("\n"
             "Number of cells:               %d\n"
             "Number of faces:               %d\n"
             "Number of edges:               %d\n"
             "Number of vertices:            %d\n",
            mesh->n_cells,
            mesh->n_faces,
            mesh->n_edges,
            mesh->n_vertices);

  if (mesh->n_vertices > 0) {

    ple_printf("\n"
               "Pointers to shareable arrays:\n"
               "  vertex_coords:        %p\n",
               mesh->vertex_coords);

    ple_printf("\n"
               "Pointers to local arrays:\n"
               "  _vertex_coords:       %p\n",
               mesh->_vertex_coords);

    /* Output coordinates */

    ple_printf("\nVertex coordinates:\n\n");
    switch(mesh->dim) {
    case 1:
      for (i = 0; i < mesh->n_vertices; i++)
        ple_printf("%10d : %12.5f\n",
                   num_vertex++, (double)(coord[i]));
      break;
    case 2:
      for (i = 0; i < mesh->n_vertices; i++)
        ple_printf("%10d : %12.5f %12.5f\n",
                   num_vertex++, (double)(coord[i*2]),
                   (double)(coord[i*2+1]));
      break;
    case 3:
      for (i = 0; i < mesh->n_vertices; i++)
        ple_printf("%10d : %12.5f %12.5f %12.5f\n",
                   num_vertex++, (double)(coord[i*3]),
                   (double)(coord[i*3+1]), (double)(coord[i*3+2]));
      break;
    default:
      ple_printf("coordinates not output\n"
                 "dimension = %d unsupported\n", mesh->dim);
    }
  }

  /* Dump element sections */
  /*-----------------------*/

  for (i = 0; i < mesh->n_sections; i++)
    _my_ple_mesh_section_dump(mesh->sections[i]);
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
