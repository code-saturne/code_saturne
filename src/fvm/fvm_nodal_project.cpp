/*============================================================================
 * Triangulation of nodal mesh sections
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

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

#include "bft/bft_printf.h"
#include "base/cs_mem.h"

#include "fvm/fvm_defs.h"
#include "fvm/fvm_io_num.h"
#include "fvm/fvm_nodal.h"
#include "fvm/fvm_nodal_priv.h"
#include "fvm/fvm_nodal_project.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm/fvm_nodal_triangulate.h"

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
 * Create a new mesh section from projection of a given extruded section
 * to its base plane.
 *
 * The algorithm used here always leads to one edge per face.
 * The global numbering of the new section is not generated here, as it may
 * simply be transferred from the initial section.
 *
 * parameters:
 *   dim              <-- spatial dimension
 *   vertex_coords    <-- associated vertex coordinates array
 *   parent_vertex_id <-- optional indirection to vertex coordinates
 *   base_section     <-- pointer to structure that should be triangulated
 *   error_count      --> number of triangulation errors counter (optional)
 *
 * returns:
 *  pointer to created nodal mesh section representation structure
 *----------------------------------------------------------------------------*/

static fvm_nodal_section_t *
_faces_to_edges(int                         dim,
                int                         chosen_axis,
                const cs_coord_t            vertex_coords[],
                const cs_lnum_t             parent_vertex_id[],
                const fvm_nodal_section_t  *base_section,
                bool                       *selected_vertices,
                cs_lnum_t                  *error_count)
{
  cs_lnum_t n_vertices, n_elements;
  cs_lnum_t i, j, k, vertex_id;

  fvm_nodal_section_t *ret_section = nullptr;

  /* Initialization */

  if (error_count != nullptr)
    *error_count = 0;

  n_elements = base_section->n_elements;

  /* Create new section */
  /*--------------------*/

  ret_section = fvm_nodal_section_create(FVM_EDGE);

  ret_section->n_elements = base_section->n_elements;
  ret_section->stride = 2;
  ret_section->connectivity_size =   ret_section->stride
                                   * ret_section->n_elements;

  CS_MALLOC(ret_section->_vertex_num,
            ret_section->connectivity_size,
            cs_lnum_t);
  ret_section->vertex_num = ret_section->_vertex_num;

  if (base_section->parent_element_id != nullptr) {

    CS_MALLOC(ret_section->_parent_element_id,
              ret_section->n_elements,
              cs_lnum_t);
    ret_section->parent_element_id = ret_section->_parent_element_id;

  }

  /* Main loop on section face elements */
  /*------------------------------------*/

  for (i = 0 ; i < n_elements ; i++) {

    cs_lnum_t tmp_id[2], tmp_selected_edge[2];
    cs_coord_t tmp_edge_center;

    cs_lnum_t selected_edge[2] = {0, 0};
    cs_coord_t selected_edge_center = 0;

    /* Compute number of vertices on each face */

    if (base_section->vertex_index != nullptr) {
      n_vertices =   base_section->vertex_index[i+1]
                   - base_section->vertex_index[i];
      vertex_id = base_section->vertex_index[i];
    }
    else {
      n_vertices = base_section->stride;
      vertex_id = base_section->stride * i;
    }

    /* Initialize the edge selection */

    selected_edge[0] = base_section->vertex_num[vertex_id + n_vertices - 1] - 1;
    selected_edge[1] = base_section->vertex_num[vertex_id] - 1;

    if (parent_vertex_id == nullptr)
      selected_edge_center = 0.5 *
        (  vertex_coords[selected_edge[0]*dim + chosen_axis]
         + vertex_coords[selected_edge[1]*dim + chosen_axis]);

    else {

      tmp_id[0] = parent_vertex_id[selected_edge[0]];
      tmp_id[1] = parent_vertex_id[selected_edge[1]];

      selected_edge_center = 0.5 *
        ( vertex_coords[tmp_id[0] * dim + chosen_axis]
        + vertex_coords[tmp_id[1] * dim + chosen_axis] );

    }

    /* Loop on face's vertices */

    for (k = 1; k < n_vertices ; k++) {

      tmp_selected_edge[0] = base_section->vertex_num[vertex_id + k - 1] - 1;
      tmp_selected_edge[1] = base_section->vertex_num[vertex_id + k] - 1;

      if (parent_vertex_id == nullptr)
        tmp_edge_center = 0.5 *
          (  vertex_coords[tmp_selected_edge[0]*dim + chosen_axis]
           + vertex_coords[tmp_selected_edge[1]*dim + chosen_axis]);

      else {

        tmp_id[0] = parent_vertex_id[tmp_selected_edge[0]];
        tmp_id[1] = parent_vertex_id[tmp_selected_edge[1]];

        tmp_edge_center
          = 0.5 * (  vertex_coords[tmp_id[0]*dim + chosen_axis]
                   + vertex_coords[tmp_id[1]*dim + chosen_axis]);

      }

      if (tmp_edge_center < selected_edge_center) {

        selected_edge_center = tmp_edge_center;

        selected_edge[0] = tmp_selected_edge[0];
        selected_edge[1] = tmp_selected_edge[1];

      }

    } /* End of loop on element's vertices */

    selected_vertices[selected_edge[0]] = true;
    selected_vertices[selected_edge[1]] = true;

    /* Define vertex_num */

    for (j = 0; j < 2; j++)
      ret_section->_vertex_num[i*2 + j] = selected_edge[j] + 1;

    if (base_section->parent_element_id != nullptr)
      ret_section->_parent_element_id[i]
        = base_section->parent_element_id[i];

  } /* End of loop on elements of the section */

  /* Return new section */

  return ret_section;
}

/*----------------------------------------------------------------------------
 * Delete unused vertices after extracting edges from faces.
 *
 * parameters:
 *   this_nodal        <-> nodal mesh structure
 *   selected_vertices <-- number of triangulation errors counter (optional)
 *
 * returns:
 *  pointer to created nodal mesh section representation structure
 *----------------------------------------------------------------------------*/

static void
_compact_mesh(fvm_nodal_t   *this_nodal,
              bool          *selected_vertices)
{
  int i, j, section_id;

  int idx = 0;

  cs_lnum_t    *old_to_new = nullptr;
  cs_lnum_t    *new_to_old = nullptr;
  cs_lnum_t    *new_parent_vtx_id = nullptr;
  cs_coord_t  *new_coords = nullptr;
  fvm_io_num_t *new_vtx_io_num = nullptr;

  cs_lnum_t n_vertices_new = 0;
  cs_lnum_t n_vertices_old = this_nodal->n_vertices;

  const int dim = this_nodal->dim;
  const cs_gnum_t *old_global_vtx_num = nullptr;

  /* Count the new number of vertices */

  for (i = 0; i < n_vertices_old; i++) {
    if (selected_vertices[i] == true)
      idx++;
  }

  n_vertices_new = idx;


  CS_MALLOC(new_to_old, n_vertices_new, cs_lnum_t);
  CS_MALLOC(old_to_new, n_vertices_old, cs_lnum_t);

  for (i = 0, idx = 0; i < n_vertices_old; i++) {

    old_to_new[i] = -1;

    if (selected_vertices[i] == true) {
      new_to_old[idx] = i + 1;
      old_to_new[i] = idx + 1;
      idx++;
    }

  }

  if (n_vertices_new != n_vertices_old) {

    /* Adapt nodal structure if required */
    /* --------------------------------- */

    if (this_nodal->_vertex_coords != nullptr) {

      /* generate a new coordinate array */

      CS_MALLOC(new_coords, dim * n_vertices_new, cs_coord_t);

      if (this_nodal->_parent_vertex_id != nullptr) {
        CS_FREE(this_nodal->_parent_vertex_id);
        this_nodal->parent_vertex_id = nullptr;
      }

      for (i = 0, idx = 0; i < n_vertices_old; i++) {

        if (selected_vertices[i] == true) {
          for (j = 0; j < dim; j++)
            new_coords[idx*dim + j] = this_nodal->vertex_coords[i*dim + j];
          idx++;
        }

      }

    } /* End if _coords != nullptr */

    else { /* this_nodal->_coords == nullptr */

      if (this_nodal->parent_vertex_id != nullptr) {

        /* generate a new parent_vertex_id */

        CS_MALLOC(new_parent_vtx_id, n_vertices_new, cs_lnum_t);

        for (i = 0, idx = 0; i < n_vertices_old; i++) {

          if (selected_vertices[i] == true)
            new_parent_vtx_id[idx++] = this_nodal->parent_vertex_id[i];

        }

        if (this_nodal->_parent_vertex_id != nullptr)
          CS_FREE(this_nodal->_parent_vertex_id);

        this_nodal->_parent_vertex_id = new_parent_vtx_id;
        this_nodal->parent_vertex_id = this_nodal->_parent_vertex_id;


      } /* End if parent_vertex_id != nullptr */

    } /* End if this_nodal->_coords == nullptr */

    /* Adapt _vertex_num in each section */
    /* --------------------------------- */

    for (section_id = 0; section_id < this_nodal->n_sections; section_id++) {

      fvm_nodal_section_t *section = this_nodal->sections[section_id];
      const cs_lnum_t n_connect = section->stride * section->n_elements;

      if (section->type == FVM_EDGE) {

        if (section->_vertex_num == nullptr)
          CS_MALLOC(section->_vertex_num, n_connect, cs_lnum_t);

        for (j = 0; j < n_connect; j++)
          section->_vertex_num[j] = old_to_new[section->vertex_num[j] - 1];

        section->vertex_num = section->_vertex_num;
      }

    }

  } /* End if n_vertices_old != n_vertices_new */

  /* Adapt global_vertex_num if required */
  /* ----------------------------------- */

  if (this_nodal->global_vertex_num != nullptr) {

    old_global_vtx_num = fvm_io_num_get_global_num(this_nodal->global_vertex_num);

    new_vtx_io_num = fvm_io_num_create(new_to_old,
                                       old_global_vtx_num,
                                       n_vertices_new,
                                       0);

    fvm_io_num_destroy(this_nodal->global_vertex_num);

  }

  this_nodal->global_vertex_num = new_vtx_io_num;
  this_nodal->n_vertices = n_vertices_new;

  CS_FREE(old_to_new);
  CS_FREE(new_to_old);

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Project an extruded mesh to its base plane.
 *
 * This is currently only possible for 2D element topologies.
 *
 * parameters:
 *   this_nodal        <-> pointer to structure that should be cut in edges
 *   chosen_axis       <-- indicate which axis is selected to extract edges
 *   error_count       --> number of triangulation errors counter (optional)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_project(fvm_nodal_t  *this_nodal,
                  int           chosen_axis,
                  cs_lnum_t    *error_count)
{
  int i;

  cs_lnum_t n_vertices = 0;
  cs_lnum_t n_edges = 0;
  cs_lnum_t section_error_count = 0;

  bool *selected_vertices = nullptr;

  assert(this_nodal != nullptr);

  n_vertices = this_nodal->n_vertices;
  CS_MALLOC(selected_vertices, n_vertices, bool);

  for (i = 0; i < n_vertices; i++)
    selected_vertices[i] = false;

  /* Now extract edges and update new section list */

  for (i = 0; i < this_nodal->n_sections; i++) {

    fvm_nodal_section_t *e_section = nullptr;
    fvm_nodal_section_t *_section = this_nodal->sections[i];

    if (_section->entity_dim == 2) {

      e_section = _faces_to_edges(this_nodal->dim,
                                  chosen_axis,
                                  this_nodal->vertex_coords,
                                  this_nodal->parent_vertex_id,
                                  _section,
                                  selected_vertices,
                                  &section_error_count);

      if (error_count != nullptr)
        *error_count += section_error_count;

      /* Transfer global numbering if present */

      if (_section->global_element_num != nullptr) {
        e_section->global_element_num = _section->global_element_num;
        _section->global_element_num = nullptr;
      }

      fvm_nodal_section_destroy(_section);
      this_nodal->sections[i] = e_section;

      n_edges += e_section->n_elements;

    }
    else if (_section->entity_dim == 3) {

      char _name[] = "-";
      const char *name = this_nodal->name;
      if (name == nullptr)
        name = _name;
      bft_error(__FILE__, __LINE__, 0,
                _("%s: not implemented for element sections with 3D topology\n"
                  " here for mesh %s and element type %s."),
                __func__, name, fvm_element_type_name[_section->type]);

    }

  }

  /* Adapt parent_vertex_id, coords and global_vertex_num */

  _compact_mesh(this_nodal,
                selected_vertices);

  this_nodal->n_faces = 0;
  this_nodal->n_edges = n_edges;

  CS_FREE(selected_vertices);
}

/*----------------------------------------------------------------------------
 * Reduce the spatial dimension of a mesh, discarding the last coordinate
 * component.
 *
 * The mesh's spatial dimension is reduced by 1.
 *
 * parameters:
 *   this_nodal <-> pointer to structure that projected
 *   matrix     <-- projection matrix
 *                  3D -> 2D: (a11, a12, a13, a21, a22, a23)
 *                  2D -> 1D: (a11, a12)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_project_coords(fvm_nodal_t  *this_nodal,
                         double       matrix[])
{
  int i;

  int new_dim = 0, old_dim, entity_dim = 0;
  cs_lnum_t n_vertices = 0;
  cs_coord_t *new_coords = nullptr;

  assert(this_nodal != nullptr);

  n_vertices = this_nodal->n_vertices;
  new_dim = this_nodal->dim - 1;
  old_dim = this_nodal->dim;

  entity_dim = fvm_nodal_get_max_entity_dim(this_nodal);

  if (entity_dim > new_dim)
    bft_error(__FILE__, __LINE__, 0,
              _("Projecting coordinates is not allowed for a mesh\n"
                "containing entities of dimension %d, as its\n"
                "spatial dimension would be reduced to %d"),
              entity_dim, new_dim);

  CS_MALLOC(new_coords, n_vertices*new_dim, cs_coord_t);

  if (old_dim == 3) {

    if (this_nodal->parent_vertex_id  == nullptr) {
      for (i = 0; i < n_vertices; i++) {
        const double *c = this_nodal->vertex_coords + i*3;
        new_coords[i*2]     = matrix[0]*c[0] + matrix[1]*c[1] + matrix[2]*c[2];
        new_coords[i*2 + 1] = matrix[3]*c[0] + matrix[4]*c[1] + matrix[5]*c[2];
      }
    }
    else {
      for (i = 0; i < n_vertices; i++) {
        const double *c =   this_nodal->vertex_coords
                          + this_nodal->parent_vertex_id[i]*3;
        new_coords[i*2]     = matrix[0]*c[0] + matrix[1]*c[1] + matrix[2]*c[2];
        new_coords[i*2 + 1] = matrix[3]*c[0] + matrix[4]*c[1] + matrix[5]*c[2];
      }
    }

  }

  else if (old_dim == 2) {

    if (this_nodal->parent_vertex_id  == nullptr) {
      for (i = 0; i < n_vertices; i++) {
        const double *c = this_nodal->vertex_coords + i*2;
        new_coords[i] = matrix[0]*c[0] + matrix[1]*c[1];
      }
    }
    else {
      for (i = 0; i < n_vertices; i++) {
        const double *c =   this_nodal->vertex_coords
                          + this_nodal->parent_vertex_id[i]*2;
        new_coords[i] = matrix[0]*c[0] + matrix[1]*c[1];
      }
    }

  }

  else
    bft_error(__FILE__, __LINE__, 0,
              _("Projecting coordinates is only allowed for a mesh\n"
                "of initial spatial dimension %d"),
              old_dim);

  /* Update structure */

  this_nodal->dim = new_dim;

  if (this_nodal->_vertex_coords != nullptr)
    CS_FREE(this_nodal->_vertex_coords);

  this_nodal->parent_vertex_id = nullptr;
  if (this_nodal->_parent_vertex_id != nullptr)
    CS_FREE(this_nodal->_parent_vertex_id);

  this_nodal->vertex_coords = new_coords;
  this_nodal->_vertex_coords = new_coords;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
