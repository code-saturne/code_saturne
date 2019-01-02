/*============================================================================
 * Extrusion of a nodal representation associated with a mesh
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

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_nodal_extrude.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Static global variables
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Extrude strided section.
 *
 * Entity parent numbering is removed.
 *
 * parameters:
 *   this_section      <-> pointer to structure that should be extruded
 *----------------------------------------------------------------------------*/

static void
_extrude_strided_section(fvm_nodal_section_t  * this_section,
                         const cs_lnum_t        n_layers)
{
  int stride;
  size_t  connectivity_size;
  int k;
  cs_lnum_t   i, j;
  cs_lnum_t   base_vertex_id;
  cs_lnum_t   element_shift, layer_shift, bottom_shift, top_shift;
  cs_lnum_t *vertex_num;

  cs_lnum_t n_elements = this_section->n_elements;
  const cs_lnum_t   n_planes = n_layers + 1;

  /* Build new connectivity */

  stride = this_section->stride * 2;
  connectivity_size = this_section->n_elements * stride * n_layers;

  BFT_MALLOC(vertex_num, connectivity_size, cs_lnum_t);
  this_section->connectivity_size = 0;

  for (i = 0; i < this_section->n_elements; i++) {
    element_shift = n_layers * stride * i;
    for (j = 0; j < n_layers; j++) {
      layer_shift = j * stride;
      bottom_shift = element_shift + layer_shift;
      top_shift = element_shift + layer_shift + this_section->stride;
      for (k = 0; k < this_section->stride; k++) {
        base_vertex_id
          = this_section->vertex_num[this_section->stride*i + k] - 1;
        vertex_num[bottom_shift + k]
          = n_planes*base_vertex_id + j + 1;
        vertex_num[top_shift + k]
          = n_planes*base_vertex_id + j + 2;
      }
    }
  }

  this_section->connectivity_size = connectivity_size;
  /* Replace old connectivity */

  if (this_section->_vertex_num != NULL)
    BFT_FREE(this_section->_vertex_num);

  this_section->_vertex_num = vertex_num;
  this_section->vertex_num = this_section->_vertex_num;

  this_section->connectivity_size = connectivity_size;

  /* Remove old attributes */

  BFT_FREE(this_section->gc_id);
  BFT_FREE(this_section->tag);

  /* Remove old parent numbering */

  this_section->parent_element_num = NULL;
  if (this_section->_parent_element_num != NULL)
    BFT_FREE(this_section->_parent_element_num);

  /* Update global_numbering */

  if (this_section->global_element_num != NULL) {

    /* Create new global numbering */

    cs_gnum_t   *global_element_num = NULL;

    const cs_gnum_t *old_global_element_num
      = fvm_io_num_get_global_num(this_section->global_element_num);

    BFT_MALLOC(global_element_num, n_elements*n_layers, cs_gnum_t);

    for (i = 0; i < n_elements; i++) {
      cs_gnum_t   base_num = (  (old_global_element_num[i]-1)
                              * (cs_gnum_t)n_layers) + 1;
      for (j = 0; j < n_layers; j++)
        global_element_num[i*n_layers + j] = base_num + (cs_gnum_t)j;
    }

    /* Remplace old global number with new */

    fvm_io_num_destroy(this_section->global_element_num);
    this_section->global_element_num = fvm_io_num_create(NULL,
                                                         global_element_num,
                                                         n_elements * n_layers,
                                                         0);

  }

  /* Update section info */

  this_section->n_elements *= n_layers;

  switch(this_section->type) {
  case FVM_EDGE:
    this_section->type = FVM_FACE_QUAD;
    break;
  case FVM_FACE_TRIA:
    this_section->type = FVM_CELL_PRISM;
    break;
  case FVM_FACE_QUAD:
    this_section->type = FVM_CELL_HEXA;
    break;
  default:
    assert(0);
  }

  this_section->entity_dim += 1;
  this_section->stride *=2;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Extrude nodal mesh.
 *
 * Vertex and element parent numbering is removed if present.
 *
 * Note: layout of new elements in memory is such that the definitions
 *       of all elements extruded from a same ancestor are contiguous.
 *       that is, {e_1, e_2, ..., e_n} leads to
 *       {e_1_layer_1, ..., e_1_layer_m, e_2_layer_1, ... e_n_layer_m}
 *
 * parameters:
 *   this_nodal       <-> pointer to structure that should be extruded
 *   n_layers          <-> number of extruded layers
 *   extrusion_vectors <-> length and direction of extrusion for each vertex;
 *                         size: mesh_spatial_dim . n_vertices
 *   distribution      <-> optional distribution of resulting vertices
 *                         along each extrusion vector (size: n_layers + 1)
 *                         with values ranging from 0 to 1, or NULL for
 *                         a regular distribution.
 *----------------------------------------------------------------------------*/

void
fvm_nodal_extrude(fvm_nodal_t        *this_nodal,
                  const cs_lnum_t     n_layers,
                  const cs_coord_t    extrusion_vectors[],
                  const cs_coord_t    distribution[])
{
  int dim, k;
  cs_lnum_t   i, j;
  cs_lnum_t   n_vertices;
  cs_lnum_t   vertex_shift;
  cs_coord_t  *_distrib = NULL;
  cs_coord_t  *new_coords = NULL;

  const cs_lnum_t   n_planes = n_layers + 1;
  const cs_coord_t  *distrib = distribution;
  const cs_coord_t  *old_coords = NULL;

  assert(this_nodal != NULL);
  assert(extrusion_vectors != NULL || this_nodal->n_vertices == 0);

  dim = this_nodal->dim;

  /* Check that no section is of too high dimension */

  for (i = 0; i < this_nodal->n_sections; i++) {
    const fvm_nodal_section_t *_section = this_nodal->sections[i];
    if (_section->entity_dim >= dim)
      bft_error(__FILE__, __LINE__, 0,
                _("Dimension of mesh \"%s\" section %d equals %d\n"
                  "with mesh spatial dimension %d prior to extrusion\n"
                  "when it should be smaller."),
                this_nodal->name, i+1, _section->entity_dim, dim);
  }

  /* Set distribution if necessary */

  if (distribution == NULL) {
    BFT_MALLOC(_distrib, n_planes, cs_coord_t);
    for (i = 0; i < n_planes; i++)
      _distrib[i] = ((double)i) / ((double)n_layers);
    distrib = _distrib;
  }

  /* Compute new coordinates */

  n_vertices = this_nodal->n_vertices;
  old_coords = this_nodal->vertex_coords;

  BFT_MALLOC(new_coords, n_planes*n_vertices*dim, cs_coord_t);

  if (this_nodal->_parent_vertex_num != NULL) {

    for (i = 0; i < n_vertices; i++) {
      const double *_old_coords
        = old_coords + ((this_nodal->parent_vertex_num[i]-1) * dim);
      vertex_shift = n_planes * dim * i;
      for (j = 0; j < n_planes; j++) {
        for (k = 0; k < dim; k++) {
          new_coords[vertex_shift + (j*dim) + k]
            =   _old_coords[k]
              + (extrusion_vectors[i*dim + k] * distrib[j]);
        }
      }
    }

  }
  else {

    for (i = 0; i < n_vertices; i++) {
      vertex_shift = n_planes * dim * i;
      for (j = 0; j < n_planes; j++) {
        for (k = 0; k < dim; k++) {
          new_coords[vertex_shift + (j*dim) + k]
            =   old_coords[i*dim + k]
              + (extrusion_vectors[i*dim + k] * distrib[j]);
        }
      }
    }

  }

  /* Replace old coords with new */

  if (this_nodal->_vertex_coords != NULL)
    BFT_FREE(this_nodal->_vertex_coords);

  this_nodal->_vertex_coords = new_coords;
  this_nodal->vertex_coords = this_nodal->_vertex_coords;

  this_nodal->parent_vertex_num = NULL;
  if (this_nodal->_parent_vertex_num != NULL)
    BFT_FREE(this_nodal->_parent_vertex_num);

  /* Update global numbering */

  if (this_nodal->global_vertex_num != NULL) {

    /* Create new global numbering */

    cs_gnum_t   *global_vertex_num = NULL;

    const cs_gnum_t *old_global_vertex_num
      = fvm_io_num_get_global_num(this_nodal->global_vertex_num);

    BFT_MALLOC(global_vertex_num, n_planes*n_vertices, cs_gnum_t);

    for (i = 0; i < n_vertices; i++) {
      cs_gnum_t   base_num = (  (old_global_vertex_num[i]-1)
                              * (cs_gnum_t)n_planes) + 1;
      for (j = 0; j < n_planes; j++)
        global_vertex_num[i*n_planes + j] = base_num + (cs_gnum_t)j;
    }

    /* Remplace old global number with new */

    fvm_io_num_destroy(this_nodal->global_vertex_num);
    this_nodal->global_vertex_num = fvm_io_num_create(NULL,
                                                      global_vertex_num,
                                                      n_vertices * n_planes,
                                                      0);

  }

  /* We may now update the number of vertices */

  this_nodal->n_vertices = n_vertices * n_planes;

  /* Extrude element definitions */

  this_nodal->n_cells = 0;
  this_nodal->n_faces = 0;
  this_nodal->n_edges = 0;

  for (i = 0; i < this_nodal->n_sections; i++) {

    fvm_nodal_section_t *_section = this_nodal->sections[i];

    if (_section->stride > 0)
      _extrude_strided_section(_section,
                               n_layers);

    else
      bft_error(__FILE__, __LINE__, 0,
                _("Extrusion of non strided sections not implemented yet."));

    switch(_section->entity_dim) {
    case 3:
      this_nodal->n_cells += _section->n_elements;
      break;
    case 2:
      this_nodal->n_faces += _section->n_elements;
      break;
    default:
      assert(0);
    }

  }

  /* If the mesh contains only vertices and no elements, a section
     of edges linking "extruded" vertices should be created */

  if (this_nodal->n_vertices != 0 && this_nodal->n_sections == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Extrusion of vertices only to edges not implemented yet."));

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
