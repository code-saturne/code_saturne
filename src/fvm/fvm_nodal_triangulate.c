/*============================================================================
 * Triangulation of nodal mesh sections
 *============================================================================*/

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

#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_defs.h"
#include "fvm_io_num.h"
#include "fvm_nodal.h"
#include "fvm_nodal_priv.h"
#include "fvm_triangulate.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_nodal_triangulate.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Static global variables
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a new mesh section from the triangulation of a given section.
 *
 * parameters:
 *   dim               <-- spatial dimension
 *   vertex_coords     <-- associated vertex coordinates array
 *   parent_vertex_num <-- optional indirection to vertex coordinates
 *   base_section      <-- pointer to structure that should be triangulated
 *   base_element_num  <-- starting number for element numbers if no
 *                         parent element number is available
 *   error_count       --> number of triangulation errors counter (optional)
 *
 * returns:
 *  pointer to created nodal mesh section representation structure
 *----------------------------------------------------------------------------*/

static fvm_nodal_section_t *
_triangulate_section(int                         dim,
                     const cs_coord_t            vertex_coords[],
                     const cs_lnum_t             parent_vertex_num[],
                     const fvm_nodal_section_t  *base_section,
                     cs_lnum_t                   base_element_num,
                     cs_lnum_t                  *error_count)
{
  cs_lnum_t n_vertices, n_triangles, n_elements;
  cs_lnum_t n_vertices_max, n_triangles_tot;
  cs_lnum_t i, j, k, triangle_id, vertex_id;

  fvm_triangulate_state_t *state = NULL;
  cs_lnum_t *n_sub_elements = NULL;
  fvm_nodal_section_t *ret_section = NULL;

  /* Initialization */

  if (error_count != NULL)
    *error_count = 0;

  n_elements = base_section->n_elements;

  if (base_section->global_element_num != NULL)
    BFT_MALLOC(n_sub_elements, n_elements, cs_lnum_t);

  /* Count expected total and local numbers of triangles */

  n_vertices_max = 0;

  n_triangles_tot = 0;

  if (base_section->vertex_index != NULL) {
    for (i = 0 ; i < n_elements ; i++) {
      n_vertices =   base_section->vertex_index[i+1]
                   - base_section->vertex_index[i];
      n_triangles_tot += n_vertices - 2;
      if (n_vertices > n_vertices_max)
        n_vertices_max = n_vertices;
    }
  }
  else if (base_section->stride == 4) {
    n_triangles_tot = base_section->n_elements * 2;
    n_vertices_max = 4;
  }
  else if (base_section->stride == 3) {
    n_triangles_tot = base_section->n_elements;
    n_vertices_max = 3;
  }

  /* Allocate memory and state variables */

  if (n_vertices_max > 4 && base_section->vertex_index != NULL)
    state = fvm_triangulate_state_create(n_vertices_max);

  /* Create new section */
  /*--------------------*/

  ret_section = fvm_nodal_section_create(FVM_FACE_TRIA);

  ret_section->n_elements = n_triangles_tot;
  ret_section->stride = 3;
  ret_section->connectivity_size =   ret_section->stride
                                   * ret_section->n_elements;
  BFT_MALLOC(ret_section->_vertex_num,
             ret_section->connectivity_size,
             cs_lnum_t);
  ret_section->vertex_num = ret_section->_vertex_num;
  BFT_MALLOC(ret_section->_parent_element_num,
             ret_section->n_elements,
             cs_lnum_t);
  ret_section->parent_element_num = ret_section->_parent_element_num;

  triangle_id = 0;

  /* Main loop on section face elements */
  /*------------------------------------*/

  for (i = 0 ; i < n_elements ; i++) {

    if (base_section->vertex_index != NULL) {
      n_vertices =   base_section->vertex_index[i+1]
                   - base_section->vertex_index[i];
      vertex_id = base_section->vertex_index[i];
    }
    else {
      n_vertices = base_section->stride;
      vertex_id = base_section->stride * i;
    }

    n_triangles = 0;

    /* If face must be subdivided */

    if (n_vertices >= 4) {

      if (n_vertices == 4)
        n_triangles = fvm_triangulate_quadrangle(dim,
                                                 1,
                                                 vertex_coords,
                                                 parent_vertex_num,
                                                 (  base_section->vertex_num
                                                  + vertex_id),
                                                 (  ret_section->_vertex_num
                                                  + triangle_id*3));

      else {
        n_triangles = fvm_triangulate_polygon(dim,
                                              1,
                                              n_vertices,
                                              vertex_coords,
                                              parent_vertex_num,
                                              (  base_section->vertex_num
                                               + vertex_id),
                                              FVM_TRIANGULATE_MESH_DEF,
                                              (  ret_section->_vertex_num
                                               + triangle_id*3),
                                              state);

        if (n_triangles != (n_vertices - 2) && error_count != NULL)
          *error_count += 1;
      }

      if (base_section->parent_element_num != NULL) {
        for (j = 0; j < n_triangles; j++)
          ret_section->_parent_element_num[triangle_id + j]
            = base_section->parent_element_num[i];
      }
      else {
        for (j = 0; j < n_triangles; j++)
          ret_section->_parent_element_num[triangle_id + j]
            = base_element_num + i;
      }

      triangle_id += n_triangles;

    }

    /* Otherwise, face must simply be copied */

    else if (n_vertices == 3) {

      n_triangles = 1;

      for (k = 0; k < 3; k++)
        ret_section->_vertex_num[triangle_id*3 + k]
          = base_section->vertex_num[i*3 + k];

      if (base_section->parent_element_num != NULL)
        ret_section->_parent_element_num[triangle_id]
          = base_section->parent_element_num[i];
      else
        ret_section->_parent_element_num[triangle_id]
          = base_element_num + i;

      triangle_id += 1;
    }

    if (n_sub_elements != NULL)
      n_sub_elements[i] = n_triangles;

  }

  /* Free memory and state variables */

  if (n_vertices_max > 4 && base_section->vertex_index != NULL)
    state = fvm_triangulate_state_destroy(state);

  /* Now update global numbering if present */

  if (base_section->global_element_num != NULL) {
    ret_section->global_element_num
      = fvm_io_num_create_from_sub(base_section->global_element_num,
                                   n_sub_elements);
    if (n_sub_elements != NULL)
      BFT_FREE(n_sub_elements);
  }

  /* Return new (triangulated) section */

  return ret_section;
}


/*----------------------------------------------------------------------------
 * Create new mesh sections from the triangulation of a given section.
 *
 * parameters:
 *   dim               <-- spatial dimension
 *   vertex_coords     <-- associated vertex coordinates array
 *   parent_vertex_num <-- optional indirection to vertex coordinates
 *   base_section      <-- pointer to structure that should be triangulated
 *   new_sections      --> array of triangle and quadrangle sections
 *   base_element_num  <-- starting number for element numbers if no
 *                         parent element number is available
 *   error_count       --> number of triangulation errors counter (optional)
 *----------------------------------------------------------------------------*/

static void
_triangulate_section_polygons(int                         dim,
                              const cs_coord_t            vertex_coords[],
                              const cs_lnum_t             parent_vertex_num[],
                              const fvm_nodal_section_t  *base_section,
                              fvm_nodal_section_t        *new_sections[2],
                              cs_lnum_t                   base_element_num,
                              cs_lnum_t                  *error_count)
{
  int type_id;
  cs_lnum_t n_vertices, n_triangles, n_elements;
  cs_lnum_t n_vertices_max, n_triangles_max;
  cs_lnum_t i, j, k, vertex_id;

  cs_lnum_t *triangle_vertices = NULL;
  fvm_element_t element_type[2] = {FVM_FACE_TRIA,
                                   FVM_FACE_QUAD};
  cs_lnum_t n_elements_tot[2] = {0, 0}; /* New triangles/quadrangles */
  cs_lnum_t element_id[2] = {0, 0};
  cs_lnum_t *n_sub_elements[2] = {NULL, NULL};
  fvm_triangulate_state_t *state = NULL;

  /* Initialization */

  if (error_count != NULL)
    *error_count = 0;

  new_sections[0] = NULL, new_sections[1] = NULL;

  n_elements = base_section->n_elements;

  /* Count expected total and local numbers of triangles */

  n_vertices_max = 0;
  n_triangles_max = 0;

  if (base_section->vertex_index != NULL) {
    for (i = 0 ; i < n_elements ; i++) {
      n_vertices =   base_section->vertex_index[i+1]
                   - base_section->vertex_index[i];
      if (n_vertices == 4)
        n_elements_tot[1] += 1;
      else
        n_elements_tot[0] += n_vertices - 2;
      if (n_vertices > n_vertices_max)
        n_vertices_max = n_vertices;
    }
  }
  else
    return;

  n_triangles_max = n_vertices_max - 2;

  /* Allocate memory and state variables */

  if (n_vertices_max > 4 && base_section->vertex_index != NULL) {
    BFT_MALLOC(triangle_vertices, n_triangles_max*3, cs_lnum_t);
    state = fvm_triangulate_state_create(n_vertices_max);
  }

  /* Create new sections */
  /*---------------------*/

  for (type_id = 0; type_id < 2; type_id++) {

    if (n_elements_tot[type_id] > 0) {
      fvm_nodal_section_t  *_section;
      _section = fvm_nodal_section_create(element_type[type_id]);
      _section->n_elements = n_elements_tot[type_id];
      _section->stride = fvm_nodal_n_vertices_element[element_type[type_id]];
      _section->connectivity_size =   _section->stride
                                    * _section->n_elements;
      BFT_MALLOC(_section->_vertex_num,
                 _section->connectivity_size,
                 cs_lnum_t);
      _section->vertex_num = _section->_vertex_num;
      BFT_MALLOC(_section->_parent_element_num,
                 _section->n_elements,
                 cs_lnum_t);
      _section->parent_element_num = _section->_parent_element_num;
      new_sections[type_id] = _section;
      if (base_section->global_element_num != NULL)
        BFT_MALLOC(n_sub_elements[type_id], n_elements, cs_lnum_t);
    }

  }

  /* Main loop on section face elements */
  /*------------------------------------*/

  for (i = 0 ; i < n_elements ; i++) {

    fvm_nodal_section_t  *_section;

    if (base_section->vertex_index != NULL) {
      n_vertices =   base_section->vertex_index[i+1]
                   - base_section->vertex_index[i];
      vertex_id = base_section->vertex_index[i];
    }
    else {
      n_vertices = base_section->stride;
      vertex_id = base_section->stride * i;
    }

    if (n_vertices == 4) {
      type_id = 1;
    }
    else {
      type_id = 0;
    }
    _section = new_sections[type_id];

    /* If face must be subdivided */

    if (n_vertices > 4) {

      n_triangles = fvm_triangulate_polygon(dim,
                                            1,
                                            n_vertices,
                                            vertex_coords,
                                            parent_vertex_num,
                                            (  base_section->vertex_num
                                             + vertex_id),
                                            FVM_TRIANGULATE_MESH_DEF,
                                            (  _section->_vertex_num
                                             + element_id[0]*3),
                                            state);

      if (n_triangles != (n_vertices - 2) && error_count != NULL)
        *error_count += 1;

      if (base_section->parent_element_num != NULL) {
        for (j = 0; j < n_triangles; j++)
          _section->_parent_element_num[element_id[0] + j]
            = base_section->parent_element_num[i];
      }
      else {
        for (j = 0; j < n_triangles; j++)
          _section->_parent_element_num[element_id[0] + j]
            = base_element_num + i;
      }

      element_id[0] += n_triangles;

      /* Update sub-element counts */

      if (n_sub_elements[0] != NULL)
        (n_sub_elements[0])[i] = n_triangles;

      if (n_sub_elements[1] != NULL)
        (n_sub_elements[1])[i] = 0;

    }

    /* Otherwise, face must simply be copied */

    else {

      for (k = 0; k < _section->stride; k++)
        _section->_vertex_num[element_id[type_id]*_section->stride + k]
          = base_section->vertex_num[i*_section->stride + k];

      if (base_section->parent_element_num != NULL)
        _section->_parent_element_num[element_id[type_id]]
          = base_section->parent_element_num[i];
      else {
        _section->_parent_element_num[element_id[type_id]]
          = base_element_num + i;
      }

      element_id[type_id] += 1;

      /* Update sub-element counts */

      if (n_sub_elements[type_id] != NULL)
        (n_sub_elements[type_id])[i] = 1;

      if (n_sub_elements[(type_id+1)%2] != NULL)
        (n_sub_elements[(type_id+1)%2])[i] = 0;

    }
  }

  /* Free memory and state variables */

  if (n_vertices_max > 4 && base_section->vertex_index != NULL) {
    BFT_FREE(triangle_vertices);
    state = fvm_triangulate_state_destroy(state);
  }

  /* Now update global numbering if present */

  for (type_id = 0; type_id < 2; type_id++) {
    if (n_sub_elements[type_id] != NULL) {
      (new_sections[type_id])->global_element_num
        = fvm_io_num_create_from_sub(base_section->global_element_num,
                                     n_sub_elements[type_id]);
    }
    BFT_FREE(n_sub_elements[type_id]);
  }

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Triangulate all sections of a nodal mesh.
 *
 * parameters:
 *   this_nodal        <-> pointer to structure that should be triangulated
 *   error_count       --> number of triangulation errors counter (optional)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_triangulate(fvm_nodal_t  *this_nodal,
                      cs_lnum_t    *error_count)
{
  int i;
  cs_lnum_t j;

  cs_lnum_t base_element_num = 1;
  cs_lnum_t section_error_count = 0;
  cs_lnum_t n_faces = 0;

  assert(this_nodal != NULL);

  /* Now triangulate and update new section list */

  for (i = 0; i < this_nodal->n_sections; i++) {

    fvm_nodal_section_t *t_section;
    fvm_nodal_section_t *_section = this_nodal->sections[i];

    if (_section->entity_dim == 2 && _section->type != FVM_FACE_TRIA) {

      t_section = _triangulate_section(this_nodal->dim,
                                       this_nodal->vertex_coords,
                                       this_nodal->parent_vertex_num,
                                       _section,
                                       base_element_num,
                                       &section_error_count);

      if (error_count != NULL)
        *error_count += section_error_count;

      base_element_num += _section->n_elements;

      fvm_nodal_section_destroy(_section);
      this_nodal->sections[i] = t_section;

      n_faces += t_section->n_elements;

    }
    else {

      if (_section->entity_dim == 2)
        n_faces +=_section->n_elements;

      if (_section->parent_element_num == NULL) {
        BFT_MALLOC(_section->_parent_element_num,
                   _section->n_elements,
                   cs_lnum_t);
        for (j = 0; j < _section->n_elements; j++)
          _section->_parent_element_num[j] = base_element_num + j;
        _section->parent_element_num = _section->_parent_element_num;
      }

      base_element_num += _section->n_elements;

    }

  }

  this_nodal->n_faces = n_faces;
}

/*----------------------------------------------------------------------------
 * Triangulate polygonal sections of a nodal mesh.
 *
 * parameters:
 *   this_nodal        <-> pointer to structure that should be triangulated
 *   error_count       --> number of triangulation errors counter (optional)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_triangulate_polygons(fvm_nodal_t  *this_nodal,
                               cs_lnum_t    *error_count)
{
  int i, j;
  cs_lnum_t k;
  fvm_nodal_section_t *new_sections[2];

  int n_sections = 0;
  cs_lnum_t base_element_num = 1;
  cs_lnum_t section_error_count = 0;
  cs_lnum_t n_faces = 0;
  fvm_nodal_section_t **sections = NULL;

  assert(this_nodal != NULL);

  /* Best estimate for new section list size: polygonal sections
     may lead to 2 sections (triangles + quads) */

  for (i = 0; i < this_nodal->n_sections; i++) {
    fvm_nodal_section_t *_section = this_nodal->sections[i];
    if (_section->type == FVM_FACE_POLY)
      n_sections += 2;
    else
      n_sections += 1;
  }

  BFT_MALLOC(sections, n_sections, fvm_nodal_section_t *);

  /* Now triangulate and update new section list */

  n_sections = 0;
  base_element_num = 1;

  for (i = 0; i < this_nodal->n_sections; i++) {

    fvm_nodal_section_t *_section = this_nodal->sections[i];

    if (_section->type == FVM_FACE_POLY) {

      _triangulate_section_polygons(this_nodal->dim,
                                    this_nodal->vertex_coords,
                                    this_nodal->parent_vertex_num,
                                    _section,
                                    new_sections,
                                    base_element_num,
                                    &section_error_count);

      if (error_count != NULL)
        *error_count += section_error_count;

      base_element_num += _section->n_elements;

      fvm_nodal_section_destroy(_section);

      for (j = 0; j < 2; j++) {
        if (new_sections[j] != NULL) {
          sections[n_sections] = new_sections[j];
          n_faces += (new_sections[j])->n_elements;
          n_sections += 1;
        }
      }

    }
    else {

      if (_section->entity_dim == 2)
        n_faces += _section->n_elements;

      if (_section->parent_element_num == NULL) {
        BFT_MALLOC(_section->_parent_element_num,
                   _section->n_elements,
                   cs_lnum_t);
        for (k = 0; k < _section->n_elements; k++)
          _section->_parent_element_num[k] = base_element_num + k;
        _section->parent_element_num = _section->_parent_element_num;
      }

      base_element_num += _section->n_elements;

      sections[n_sections] = _section;
      n_sections += 1;

    }

  }

  /* Now replace section list */

  BFT_FREE(this_nodal->sections);
  BFT_REALLOC(sections, n_sections, fvm_nodal_section_t *);

  this_nodal->n_sections = n_sections;
  this_nodal->sections = sections;

  this_nodal->n_faces = n_faces;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
