/*============================================================================
 * Deal with the exchange of data included in opaque structure fvm_nodal_t
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

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_nodal.h"
#include "fvm_nodal_priv.h"
#include "fvm_nodal_extract.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/* Geometric operation macros*/

#define _DOT_PRODUCT(v0, v1) \
  (v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2])

#define _MODULE(v) \
  sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

#define _CROSS_PRODUCT(cp, v0, v1)  \
  (cp[0] = v0[1] * v1[2] - v1[1] * v0[2], \
   cp[1] = v1[0] * v0[2] - v0[0] * v1[2], \
   cp[2] = v0[0] * v1[1] - v1[0] * v0[1])

/* In case <math.h> does not define HUGE_VAL, use a "safe" value */

#if !defined(HUGE_VAL)
#define HUGE_VAL 1.0e+30
#endif

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute values associated with a 3d face.
 *
 * The algorithm will hande non-convex polygons.
 * To attempt to reduce the influence of face warping, the reference point
 * is taken not at a vertex of the polygon, but at the center of gravity
 * of the element's vertices.
 *
 * parameters:
 *   n_vertices         <-- number of face vertices
 *   face_vertex_num    <-- face vertex numbers
 *   parent_vertex_num  <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords      <-- pointer to vertex coordinates
 *   center             --> face center
 *   normal             --> face normal
 *   surface            --> face surface
 *----------------------------------------------------------------------------*/

static inline void
_face_quantities_3d(cs_lnum_t           n_vertices,
                    const cs_lnum_t     face_vertex_num[],
                    const cs_lnum_t    *parent_vertex_num,
                    const cs_coord_t    vertex_coords[],
                    double              center[3],
                    double              normal[3],
                    double             *surface)
{
  cs_lnum_t   i, j;
  double v0[3], v1[3], tc[3], tn[3], ts;
  const cs_coord_t *_coords_0, *_coords_1;

  double one_third = 1./3.;
  double vertices_center[3] = {0., 0., 0.};

  /* Initialization */

  for (j = 0; j < 3; j++) {
    center[j] = 0.;
    normal[j] = 0.;
    *surface = 0;
  }

  /* Counting loop on face vertices */

  if (parent_vertex_num != NULL) {

    for (i = 0; i < n_vertices; i++) {
      _coords_0 =   vertex_coords
                  + (parent_vertex_num[face_vertex_num[i]-1] - 1)*3;
      for (j = 0; j < 3; j++)
        vertices_center[j] += _coords_0[j];
    }

  }
  else {

    for (i = 0; i < n_vertices; i++) {
      _coords_0 = vertex_coords + (face_vertex_num[i]-1)*3;
      for (j = 0; j < 3; j++)
        vertices_center[j] += _coords_0[j];
    }

  }

  for (j = 0; j < 3; j++)
    vertices_center[j] /= n_vertices;

  /* Main loop on face vertices */

  for (i = 0; i < n_vertices; i++) {

    cs_lnum_t v_id_0 = face_vertex_num[i] - 1;
    cs_lnum_t v_id_1 = face_vertex_num[(i+1)%n_vertices] - 1;

    if (parent_vertex_num != NULL) {
      _coords_0 = vertex_coords + (parent_vertex_num[v_id_0] - 1)*3;
      _coords_1 = vertex_coords + (parent_vertex_num[v_id_1] - 1)*3;
    }
    else {
      _coords_0 = vertex_coords + v_id_0*3;
      _coords_1 = vertex_coords + v_id_1*3;
    }

    /* Triangle normal and center */

    for (j = 0; j < 3; j++) {
      v0[j] = _coords_0[j] - vertices_center[j];
      v1[j] = _coords_1[j] - vertices_center[j];
      tc[j] = (_coords_0[j] + _coords_1[j] + vertices_center[j])*one_third;
    }

    _CROSS_PRODUCT(tn, v0, v1);

    for (j = 0; j < 3; j++)
      normal[j] += tn[j]*0.5;

    ts = _MODULE(tn);

    if (_DOT_PRODUCT(tn, normal) < 0.)
      ts *= -1.;

    /* Contribution to surface and center */

    *surface += ts;

    for (j = 0; j < 3; j++)
      center[j] += ts*tc[j];
  }

  for (j = 0; j < 3; j++)
    center[j] /= *surface;
}

/*----------------------------------------------------------------------------
 * Compute values associated with a 2d face.
 *
 * The algorithm will hande non-convex polygons.
 *
 * parameters:
 *   n_vertices         <-- number of face vertices
 *   face_vertex_num    <-- face vertex numbers
 *   parent_vertex_num  <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords      <-- pointer to vertex coordinates
 *   center             --> face center
 *   surface            --> face surface
 *----------------------------------------------------------------------------*/

static inline void
_face_quantities_2d(cs_lnum_t          n_vertices,
                    const cs_lnum_t    face_vertex_num[],
                    const cs_lnum_t   *parent_vertex_num,
                    const cs_coord_t   vertex_coords[],
                    double             center[2],
                    double            *surface)
{
  cs_lnum_t   i, j, v_id_1, v_id_2;
  double v0[2], v1[2], tc[2], ts;
  const cs_coord_t *_coords_0, *_coords_1, *_coords_2;

  double one_third = 1./3.;

  /* Initialization */

  for (j = 0; j < 2; j++) {
    center[j] = 0.;
    *surface = 0;
  }

  v_id_1 = face_vertex_num[0] - 1;
  v_id_2 = face_vertex_num[1] - 1;

  if (parent_vertex_num != NULL) {
    _coords_1 = vertex_coords + (parent_vertex_num[v_id_1] - 1)*2;
    _coords_2 = vertex_coords + (parent_vertex_num[v_id_2] - 1)*2;
  }
  else {
    _coords_1 = vertex_coords + v_id_1*2;
    _coords_2 = vertex_coords + v_id_2*2;
  }

  /* Loop on face vertices */

  for (i = 0; i < n_vertices; i++) {

    _coords_0 = _coords_1;
    _coords_1 = _coords_2;

    v_id_2 = face_vertex_num[(i+2)%n_vertices] - 1;

    if (parent_vertex_num != NULL)
      _coords_2 = vertex_coords + (parent_vertex_num[v_id_2] - 1)*2;
    else
      _coords_2 = vertex_coords + v_id_2*2;

    /* Triangle normal and center */

    for (j = 0; j < 2; j++) {
      v0[j] = _coords_2[j] - _coords_1[j];
      v1[j] = _coords_0[j] - _coords_1[j];
      tc[j] = (_coords_0[j] + _coords_1[j] + _coords_2[j])*one_third;
    }

    ts = v0[0]*v1[1] - v0[1]*v1[0];

    /* Contribution to surface and center */

    *surface += ts;

    for (j = 0; j < 2; j++)
      center[j] += ts*tc[j];
  }

  for (j = 0; j < 2; j++)
    center[j] /= *surface;

}

/*----------------------------------------------------------------------------
 * Copy element centers from a polyhedral mesh section to an array.
 *
 * parameters:
 *   this_section       <-- pointer to nodal mesh section
 *   n_elements         <-- number of elements of similar dimension in mesh
 *   element_count      <-> number of elements whose centers have been copied
 *   parent_vertex_num  <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords      <-- pointer to vertex coordinates
 *   interlace          <-- indicates if destination array is interlaced
 *   cell_centers       --> cell centers coordinates (pre-allocated)
 *----------------------------------------------------------------------------*/

static void
_cell_poly_section_centers(const fvm_nodal_section_t  *this_section,
                           cs_lnum_t                   n_elements,
                           cs_lnum_t                  *element_count,
                           const cs_lnum_t            *parent_vertex_num,
                           const cs_coord_t            vertex_coords[],
                           cs_interlace_t              interlace,
                           cs_coord_t                 *cell_centers)
{
  cs_lnum_t i, j;
  int k;

  cs_lnum_t _element_count = *element_count;

  assert(this_section != NULL);

  assert(this_section->type == FVM_CELL_POLY);

  /* Loop on elements */

  for (i = 0; i < this_section->n_elements; i++) {

    double denom = 0.;
    double c_center[3] = {0., 0., 0.};

    /* Loop on faces */

    for (j = this_section->face_index[i];
         j < this_section->face_index[i + 1];
         j++) {

      double f_center[3], f_normal[3], f_surface;

      cs_lnum_t face_id = CS_ABS(this_section->face_num[j]) - 1;
      cs_lnum_t v_index_start = this_section->vertex_index[face_id];
      cs_lnum_t v_index_end = this_section->vertex_index[face_id + 1];
      cs_lnum_t n_vertices = v_index_end - v_index_start;

      _face_quantities_3d(n_vertices,
                          this_section->vertex_num + v_index_start,
                          parent_vertex_num,
                          vertex_coords,
                          f_center,
                          f_normal,
                          &f_surface);

      denom += f_surface;
      for (k = 0; k < 3; k++)
        c_center[k] += f_center[k]*f_surface;

    } /* End of loop on faces */

    if (interlace == CS_INTERLACE) {
      for (k = 0; k < 3; k++)
        cell_centers[_element_count*3 + k] = c_center[k]/denom;
    }
    else {
      for (k = 0; k < 3; k++)
        cell_centers[n_elements*k + _element_count] = c_center[k]/denom;
    }

    _element_count++;

  } /* End of loop on cells */

  *element_count = _element_count;
}

/*----------------------------------------------------------------------------
 * Copy element centers from a strided cell mesh section to an array.
 *
 * parameters:
 *   this_section       <-- pointer to nodal mesh section
 *   n_elements         <-- number of elements of similar dimension in mesh
 *   element_count      <-> number of elements whose centers have been copied
 *   parent_vertex_num  <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords      <-- pointer to vertex coordinates
 *   interlace          <-- indicates if destination array is interlaced
 *   cell_centers       --> cell centers coordinates (pre-allocated)
 *----------------------------------------------------------------------------*/

static void
_cell_strided_section_centers(const fvm_nodal_section_t  *this_section,
                              cs_lnum_t                   n_elements,
                              cs_lnum_t                  *element_count,
                              const cs_lnum_t            *parent_vertex_num,
                              const cs_coord_t            vertex_coords[],
                              cs_interlace_t              interlace,
                              cs_coord_t                 *cell_centers)
{
  cs_lnum_t i, j;
  int k;

  int n_faces, stride;
  int n_face_vertices[6], face_vertices[6][4];

  cs_lnum_t _element_count = *element_count;

  assert(this_section != NULL);

  stride = this_section->stride;

  fvm_nodal_cell_face_connect(this_section->type,
                              &n_faces,
                              n_face_vertices,
                              face_vertices);

  /* Loop on elements */

  for (i = 0; i < this_section->n_elements; i++) {

    double denom = 0.;
    double c_center[3] = {0., 0., 0.};

    /* Loop on faces */

    for (j = 0; j < n_faces; j++) {

      cs_lnum_t _vertex_num[4];
      double f_center[3], f_normal[3], f_surface;

      for (k = 0; k < n_face_vertices[j]; k++)
        _vertex_num[k]
          = this_section->vertex_num[i*stride + face_vertices[j][k]];

      _face_quantities_3d(n_face_vertices[j],
                          _vertex_num,
                          parent_vertex_num,
                          vertex_coords,
                          f_center,
                          f_normal,
                          &f_surface);

      denom += f_surface;
      for (k = 0; k < 3; k++)
        c_center[k] += f_center[k]*f_surface;

    } /* End of loop on faces */

    if (interlace == CS_INTERLACE) {
      for (k = 0; k < 3; k++)
        cell_centers[_element_count*3 + k] = c_center[k]/denom;
    }
    else {
      for (k = 0; k < 3; k++)
        cell_centers[n_elements*k + _element_count] = c_center[k]/denom;
    }

    _element_count++;

  } /* End of loop on cells */

  *element_count = _element_count;
}

/*----------------------------------------------------------------------------
 * Copy element centers from a 3d face mesh section to an array.
 *
 * parameters:
 *   this_section       <-- pointer to nodal mesh section
 *   n_elements         <-- number of elements of similar dimension in mesh
 *   element_count      <-> number of elements whose centers have been copied
 *   parent_vertex_num  <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords      <-- pointer to vertex coordinates
 *   interlace          <-- indicates if destination array is interlaced
 *   cell_centers       --> cell centers coordinates (pre-allocated)
 *----------------------------------------------------------------------------*/

static void
_face_section_centers_3d(const fvm_nodal_section_t  *this_section,
                         cs_lnum_t                   n_elements,
                         cs_lnum_t                  *element_count,
                         const cs_lnum_t            *parent_vertex_num,
                         const cs_coord_t            vertex_coords[],
                         cs_interlace_t              interlace,
                         cs_coord_t                 *cell_centers)
{
  cs_lnum_t i, v_index_start, v_index_end, n_vertices;
  int k;

  cs_lnum_t _element_count = *element_count;

  /* Loop on elements */

  for (i = 0; i < this_section->n_elements; i++) {

    double f_center[3], f_normal[3], f_surface;

    if (this_section->vertex_index != NULL) {
      v_index_start = this_section->vertex_index[i];
      v_index_end = this_section->vertex_index[i + 1];
    }
    else {
      v_index_start = this_section->stride*i;
      v_index_end = this_section->stride*(i+1);
    }

    n_vertices = v_index_end - v_index_start;

    _face_quantities_3d(n_vertices,
                        this_section->vertex_num + v_index_start,
                        parent_vertex_num,
                        vertex_coords,
                        f_center,
                        f_normal,
                        &f_surface);

    if (interlace == CS_INTERLACE) {
      for (k = 0; k < 3; k++)
        cell_centers[_element_count*3 + k] = f_center[k];
    }
    else {
      for (k = 0; k < 3; k++)
        cell_centers[n_elements*k + _element_count] = f_center[k];
    }

    _element_count++;

  } /* End of loop on cells */

  *element_count = _element_count;
}

/*----------------------------------------------------------------------------
 * Copy element centers from a 2d face mesh section to an array.
 *
 * parameters:
 *   this_section       <-- pointer to nodal mesh section
 *   n_elements         <-- number of elements of similar dimension in mesh
 *   element_count      <-> number of elements whose centers have been copied
 *   parent_vertex_num  <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords      <-- pointer to vertex coordinates
 *   interlace          <-- indicates if destination array is interlaced
 *   cell_centers       --> cell centers coordinates (pre-allocated)
 *----------------------------------------------------------------------------*/

static void
_face_section_centers_2d(const fvm_nodal_section_t  *this_section,
                         cs_lnum_t                   n_elements,
                         cs_lnum_t                  *element_count,
                         const cs_lnum_t            *parent_vertex_num,
                         const cs_coord_t            vertex_coords[],
                         cs_interlace_t              interlace,
                         cs_coord_t                 *cell_centers)
{
  cs_lnum_t i, v_index_start, v_index_end, n_vertices;
  int k;

  cs_lnum_t _element_count = *element_count;

  /* Loop on elements */

  for (i = 0; i < this_section->n_elements; i++) {

    double f_center[2], f_surface;

    if (this_section->vertex_index != NULL) {
      v_index_start = this_section->vertex_index[i];
      v_index_end = this_section->vertex_index[i + 1];
    }
    else {
      v_index_start = this_section->stride*i;
      v_index_end = this_section->stride*(i+1);
    }

    n_vertices = v_index_end - v_index_start;

    _face_quantities_2d(n_vertices,
                        this_section->vertex_num + v_index_start,
                        parent_vertex_num,
                        vertex_coords,
                        f_center,
                        &f_surface);

    if (interlace == CS_INTERLACE) {
      for (k = 0; k < 2; k++)
        cell_centers[_element_count*2 + k] = f_center[k];
    }
    else {
      for (k = 0; k < 2; k++)
        cell_centers[n_elements*k + _element_count] = f_center[k];
    }

    _element_count++;

  } /* End of loop on elements */

  *element_count = _element_count;
}

/*----------------------------------------------------------------------------
 * Update element extents with a given vertex
 *
 * parameters:
 *   dim               <-- spatial (coordinates) dimension
 *   vertex_id         <-- vertex index (0 to n-1)
 *   parent_vertex_num <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords     <-- pointer to vertex coordinates
 *   elt_extents       <-> extents associated with element:
 *                         x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *   elt_initialized   <-> are extents already initialized for this vertex
 *                         (for all element vertices except the first) ?
 *----------------------------------------------------------------------------*/

inline static void
_update_elt_extents(int                 dim,
                    cs_lnum_t           vertex_id,
                    const cs_lnum_t    *parent_vertex_num,
                    const cs_coord_t    vertex_coords[],
                    double              elt_extents[],
                    bool               *elt_initialized)
{
  cs_lnum_t   i, coord_idx;

  if (parent_vertex_num == NULL)
    coord_idx = vertex_id;
  else
    coord_idx = parent_vertex_num[vertex_id] - 1;

  if (*elt_initialized == false) {
    for (i = 0; i < dim; i++) {
      elt_extents[i]       = vertex_coords[(coord_idx * dim) + i];
      elt_extents[i + dim] = vertex_coords[(coord_idx * dim) + i];
    }
    *elt_initialized = true;
  }
  else {
    for (i = 0; i < dim; i++) {
      if (elt_extents[i]       > vertex_coords[(coord_idx * dim) + i])
        elt_extents[i]       = vertex_coords[(coord_idx * dim) + i];
      if (elt_extents[i + dim] < vertex_coords[(coord_idx * dim) + i])
        elt_extents[i + dim] = vertex_coords[(coord_idx * dim) + i];
    }

  }

}

/*----------------------------------------------------------------------------
 * Adjust element extents with search tolerance and update global extents
 *
 * parameters:
 *   dim         <-- spatial (coordinates) dimension
 *   elt_dim     <-- element dimension
 *   tolerance   <-- addition to local extents of each element:
 *                   extent = base_extent * (1 + tolerance)
 *   elt_extents <-> extents associated with element:
 *                   x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *----------------------------------------------------------------------------*/

inline static void
_elt_extents_finalize(int               dim,
                      int               elt_dim,
                      double            tolerance,
                      double  *restrict elt_extents)
{
  int i;
  double delta[3];

  for (i = 0; i < dim; i++)
    delta[i] = (elt_extents[i+dim] - elt_extents[i]) * tolerance;

  if (elt_dim < dim) {
    double delta_max = delta[0];  /* for 1d or 2d elements, ensure */
    for (i = 0; i < dim; i++) {   /* search extent "thickness" */
      if (delta[i] > delta_max)
        delta_max = delta[i];
    }
    for (i = 0; i < dim; i++)
      delta[i] = delta_max;
  }

  for (i = 0; i < dim; i++) {
    elt_extents[i]     = elt_extents[i]     - delta[i];
    elt_extents[i+dim] = elt_extents[i+dim] + delta[i];
  }

}

/*----------------------------------------------------------------------------
 * Adjust extents with sub-extents
 *
 * parameters:
 *   dim         <-- spatial (coordinates) dimension
 *   sub_extents <-> extents associated with element or section:
 *                   x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *   extents     <-> optional section or mesh extents, to be updated:
 *                   x_min, y_min, ..., x_max, y_max, ... (size: 2*dim);
 *                   NULL if unused
 *----------------------------------------------------------------------------*/

inline static void
_update_extents(int               dim,
                double  *restrict sub_extents,
                double  *restrict extents)
{
  int i;

  for (i = 0; i < dim; i++) {
    if (sub_extents[i] < extents[i])
      extents[i] = sub_extents[i];
    if (sub_extents[i+dim] > extents[i+dim])
      extents[i+dim] = sub_extents[i+dim];
  }
}

/*----------------------------------------------------------------------------
 * Compute extents of a nodal mesh representation section
 *
 * parameters:
 *   this_section      <-- pointer to section structure
 *   dim               <-- spatial (coordinates) dimension
 *   parent_vertex_num <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords     <-- pointer to vertex coordinates
 *   tolerance         <-- addition to local extents of each element:
 *                         extent = base_extent * (1 + tolerance)
 *   extents           <-> extents associated with section:
 *                         x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *----------------------------------------------------------------------------*/

static void
_nodal_section_extents(const fvm_nodal_section_t  *this_section,
                       int                         dim,
                       const cs_lnum_t            *parent_vertex_num,
                       const cs_coord_t            vertex_coords[],
                       double                      tolerance,
                       double                      extents[])
{
  cs_lnum_t   i, j, k, face_id, vertex_id;
  double elt_extents[6];

  /* initialize extents in case section is empty */
  for (i = 0; i < dim; i++) {
    extents[i]       =  HUGE_VAL;
    extents[i + dim] = -HUGE_VAL;
  }

  /* Extents for polyhedra */

  if (this_section->face_index != NULL) {

    for (i = 0; i < this_section->n_elements; i++) {

      bool elt_initialized = false;

      for (j = this_section->face_index[i];
           j < this_section->face_index[i + 1];
           j++) {
        face_id = CS_ABS(this_section->face_num[j]) - 1;
        for (k = this_section->vertex_index[face_id];
             k < this_section->vertex_index[face_id + 1];
             k++) {
          vertex_id = this_section->vertex_num[k] - 1;

          _update_elt_extents(dim,
                              vertex_id,
                              parent_vertex_num,
                              vertex_coords,
                              elt_extents,
                              &elt_initialized);

        }
      }

      _elt_extents_finalize(dim, 3, tolerance, elt_extents);
      _update_extents(dim, elt_extents, extents);

    }

  }

  /* Extents for polygons */

  else if (this_section->vertex_index != NULL) {

    cs_lnum_t   n_faces = (this_section->n_faces > 0) ?
                          this_section->n_faces : this_section->n_elements;

    for (i = 0; i < n_faces; i++) {

      bool elt_initialized = false;

      for (j = this_section->vertex_index[i];
           j < this_section->vertex_index[i + 1];
           j++) {
        vertex_id = this_section->vertex_num[j] - 1;

        _update_elt_extents(dim,
                            vertex_id,
                            parent_vertex_num,
                            vertex_coords,
                            elt_extents,
                            &elt_initialized);

      }

      _elt_extents_finalize(dim, 2, tolerance, elt_extents);
      _update_extents(dim, elt_extents, extents);

    }

  }

  /* Extents for regular elements */

  else {

    for (i = 0; i < this_section->n_elements; i++) {

      bool elt_initialized = false;

      for (j = 0; j < this_section->stride; j++) {

        vertex_id = this_section->vertex_num[i*this_section->stride + j] - 1;

        _update_elt_extents(dim,
                            vertex_id,
                            parent_vertex_num,
                            vertex_coords,
                            elt_extents,
                            &elt_initialized);

      }

      _elt_extents_finalize(dim,
                            this_section->entity_dim,
                            tolerance,
                            elt_extents);
      _update_extents(dim, elt_extents, extents);

    }
  }
}

/*============================================================================
 * Semi-private function definitions (prototypes in fvm_nodal_priv.h)
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Copy global vertex numbers to an array.
 *
 * parameters:
 *   this_nodal <-- pointer to nodal mesh structure
 *   g_vtx_num  --> global vertex numbers (pre-allocated)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_get_global_vertex_num(const fvm_nodal_t  *this_nodal,
                                cs_gnum_t          *g_vtx_num)
{
  size_t size;
  cs_lnum_t   vertex_id;
  fvm_io_num_t *global_io_num = this_nodal->global_vertex_num;

  assert(g_vtx_num != NULL || this_nodal->n_vertices == 0);

  if (g_vtx_num == NULL)
    return;

  if (global_io_num == NULL) {

    for (vertex_id = 0; vertex_id < this_nodal->n_vertices; vertex_id++)
      g_vtx_num[vertex_id] = vertex_id + 1;

  }
  else {

    size = sizeof(cs_gnum_t) * fvm_io_num_get_local_count(global_io_num);
    memcpy(g_vtx_num, fvm_io_num_get_global_num(global_io_num), size);

  }
}

/*----------------------------------------------------------------------------
 * Copy global element numbers of a given element type to an array.
 *
 * Note that if the mesh contains multiple sections of the same element type,
 * the global element numbers are continued from one section to the next,
 * so to the user, all is as if the sections were concatenated.
 *
 * parameters:
 *   this_nodal   <-- pointer to nodal mesh structure
 *   element_type <-- type of elements to deal with
 *   g_elt_num    <-> pointer to global_element_num array (pre-allocated)
 *
 * returns:
 *----------------------------------------------------------------------------*/

void
fvm_nodal_get_global_element_num(const fvm_nodal_t  *this_nodal,
                                 fvm_element_t       element_type,
                                 cs_gnum_t          *g_elt_num)
{
  int element_id, section_id;

  cs_lnum_t n_elements = 0, element_count = 0;
  cs_gnum_t n_g_elements = 0, global_count = 0;
  const cs_gnum_t *g_num = NULL;

  /* Define global element numbers */

  for (section_id = 0; section_id < this_nodal->n_sections; section_id++) {

    const fvm_nodal_section_t *section = this_nodal->sections[section_id];

    if (section->type == element_type) {

      const fvm_io_num_t *io_num = section->global_element_num;

      if (io_num == NULL) {

        for (element_id = 0; element_id < section->n_elements; element_id++)
          g_elt_num[element_count + element_id] = element_id + global_count + 1;

        element_count += section->n_elements;
        global_count += section->n_elements;

      }

      else { /* io_num != NULL */

        n_elements = fvm_io_num_get_local_count(io_num);
        n_g_elements = fvm_io_num_get_global_count(io_num);
        g_num = fvm_io_num_get_global_num(io_num);

        if (global_count == 0) {
          memcpy((char *)g_elt_num,
                 g_num,
                 sizeof(cs_gnum_t) * n_elements);

        }
        else {

          for (element_id = 0; element_id < n_elements; element_id++)
            g_elt_num[element_count + element_id]
              = g_num[element_id] + global_count;

        }

        element_count += n_elements;
        global_count += n_g_elements;

      }

    } /* section->type == element_type */

  } /* End of loop on sections */

  return;
}

/*----------------------------------------------------------------------------
 * Copy vertex coordinates to an array.
 *
 * parameters:
 *   this_nodal     <-- pointer to nodal mesh structure
 *   interlace      <-- indicates if destination array is interlaced
 *   vertex_coords  --> vertices coordinates (pre-allocated)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_get_vertex_coords(const fvm_nodal_t  *this_nodal,
                            cs_interlace_t      interlace,
                            cs_coord_t         *vertex_coords)
{
  int i;
  cs_lnum_t vertex_id;

  const int dim = this_nodal->dim;
  const cs_lnum_t n_vertices = this_nodal->n_vertices;
  const cs_coord_t *coords = this_nodal->vertex_coords;
  const cs_lnum_t   *parent_num = this_nodal->parent_vertex_num;

  if (this_nodal->parent_vertex_num == NULL) {

    if (interlace == CS_INTERLACE)
      memcpy(vertex_coords, coords, sizeof(cs_coord_t) * n_vertices * dim);

    else {
      for (i = 0; i < dim; i++) {
        for (vertex_id = 0; vertex_id < n_vertices; vertex_id++)
          vertex_coords[n_vertices*i + vertex_id]
            = coords[vertex_id*dim + i];
      }

    }

  }
  else { /* parent_vertex_num != NULL */

    if (interlace == CS_INTERLACE) {

      for (i = 0; i < dim; i++) {
        for (vertex_id = 0; vertex_id < n_vertices; vertex_id++)
          vertex_coords[vertex_id * dim + i]
            = coords[(parent_num[vertex_id]-1) * dim + i];
      }

    }
    else {

      for (i = 0; i < dim; i++) {
        for (vertex_id = 0; vertex_id < n_vertices; vertex_id++)
          vertex_coords[n_vertices*i + vertex_id]
            = coords[(parent_num[vertex_id]-1) * dim + i];
      }

    }

  }

}

/*----------------------------------------------------------------------------
 * Copy element centers to an array.
 *
 * Note that if the mesh contains multiple cell element sections of, the
 * cell_centers array spans all sections, so to the user, all is as if the
 * sections were concatenated.
 *
 * parameters:
 *   this_nodal     <-- pointer to nodal mesh structure
 *   interlace      <-- indicates if destination array is interlaced
 *   entity_dim     <-- dimension of entities we want to count (0 to 3)
 *   cell_centers   --> cell centers coordinates (pre-allocated)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_get_element_centers(const fvm_nodal_t  *this_nodal,
                              cs_interlace_t      interlace,
                              int                 entity_dim,
                              cs_coord_t         *cell_centers)
{
  int i, j;
  int section_id;
  cs_lnum_t element_id;

  cs_lnum_t element_count = 0;

  const int dim = this_nodal->dim;
  const cs_coord_t *coords = this_nodal->vertex_coords;
  const cs_lnum_t   *parent_num = this_nodal->parent_vertex_num;
  const cs_lnum_t n_elements = fvm_nodal_get_n_entities(this_nodal,
                                                         entity_dim);

  for (section_id = 0; section_id < this_nodal->n_sections; section_id++) {

    const fvm_nodal_section_t *section = this_nodal->sections[section_id];

    /* Ignore sections of different dimension than that queried */

    if (section->entity_dim != entity_dim)
      continue;

    /* 3D cells */

    if (section->entity_dim == 3) {

      if (section->type == FVM_CELL_POLY)

        _cell_poly_section_centers(section,
                                   n_elements,
                                   &element_count,
                                   this_nodal->parent_vertex_num,
                                   this_nodal->vertex_coords,
                                   interlace,
                                   cell_centers);

      else {

        assert(section->stride != 0);

        _cell_strided_section_centers(section,
                                      n_elements,
                                      &element_count,
                                      this_nodal->parent_vertex_num,
                                      this_nodal->vertex_coords,
                                      interlace,
                                      cell_centers);

      }

    }

    /* 3D or 2D faces */

    else if (section->entity_dim == 2) {

      if (dim == 3)
        _face_section_centers_3d(section,
                                 n_elements,
                                 &element_count,
                                 this_nodal->parent_vertex_num,
                                 this_nodal->vertex_coords,
                                 interlace,
                                 cell_centers);

      else /* if (dim == 2) */
        _face_section_centers_2d(section,
                                 n_elements,
                                 &element_count,
                                 this_nodal->parent_vertex_num,
                                 this_nodal->vertex_coords,
                                 interlace,
                                 cell_centers);

    }

    /* 3D, 2D, or 1D edges */

    else if (section->entity_dim == 1) {

      const int stride = section->stride;

      assert(section->stride != 0);

      for (element_id = 0;
           element_id < section->n_elements;
           element_id++, element_count++) {

        double  cell_center[3] = {0., 0., 0.};
        double  denom = 0.;

        if (this_nodal->parent_vertex_num == NULL) {
          for (j = 0; j < stride; j++) {
            const cs_lnum_t vertex_id
              = section->vertex_num[(element_id * stride) + j] - 1;
            for (i = 0; i < dim; i++)
              cell_center[i] += coords[vertex_id*dim + i];
            denom += 1.;
          }
        }
        else { /* if (this_nodal->parent_vertex_num != NULL) */
          for (j = 0; j < stride; j++) {
            const cs_lnum_t vertex_id
              = section->vertex_num[(element_id * stride) + j] - 1;
            for (i = 0; i < dim; i++)
              cell_center[i] += coords[(parent_num[vertex_id]-1)*dim + i];
            denom += 1.;
          }
        }

        if (interlace == CS_INTERLACE) {
          for (i = 0; i < dim; i++)
            cell_centers[element_count*dim + i] = cell_center[i] / denom;
        }
        else {
          for (i = 0; i < dim; i++)
            cell_centers[n_elements*i + element_count]
              = cell_center[i] / denom;
        }

      }

    }

  } /* End of loop on section_id */

}

/*----------------------------------------------------------------------------
 * Copy element -> vertex connectivity of a given element type to an array.
 *
 * Note that if the mesh contains multiple sections of the same element type,
 * the connectivity spans all sections, so to the user, all is as if the
 * sections were concatenated.
 *
 * parameters:
 *   this_nodal    <-- pointer to nodal mesh structure
 *   element_type  <-- type of elements of the section to deal with
 *   connectivity  <-> pointer to connectvity (pre-allocated)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_get_strided_connect(const fvm_nodal_t  *this_nodal,
                              fvm_element_t       element_type,
                              cs_lnum_t          *connectivity)
{
  int i, section_id;
  cs_lnum_t element_id;
  cs_lnum_t element_count = 0;

  /* Verify section with "element_type" type exists and is strided */

  if (element_type == FVM_FACE_POLY || element_type == FVM_CELL_POLY)
    bft_error(__FILE__, __LINE__, 0,
              _("Elements of type : \"%s\" are not strided elements.\n"
                "Incorrect use with fvm_nodal_get_strided_connect()\n"
                "Associated nodal mesh : \"%s\"\n"),
              fvm_elements_type_name[element_type], this_nodal->name);

  /* Retrieve connectivity */

  for (section_id = 0; section_id < this_nodal->n_sections; section_id++) {

    const fvm_nodal_section_t *section = this_nodal->sections[section_id];

    if (section->type == element_type) {

      const int stride = section->stride;
      const cs_lnum_t *num = section->vertex_num;

      for (element_id = 0; element_id < section->n_elements; element_id++) {
        for (i = 0; i < stride; i++) {
          connectivity[element_id * stride + i + element_count]
            = num[element_id*stride + i];
        }
      }

      element_count += section->n_elements * stride;

    } /* Section->type == element_type */

  } /* End of loop on sections */

}

/*----------------------------------------------------------------------------
 * Build inverse vertex -> element connectivity.
 *
 * The user is responsible for freeing the returned arrays.
 * The size of element_index[] should be n_vertices + 1, where n_vertices
 * is the value returned by fvm_nodal_get_n_entities(this_nodal, entity_dim).
 * The size of element_id[] should be element_index[n_vertices].
 *
 * Note that if the mesh contains multiple cell element sections of, the
 * cell_centers array spans all sections, so to the user, all is as if the
 * sections were concatenated.
 *
 * Note also that both the vertex -> element index and connectivity arrays
 * returned use 0 to n numbering.
 *
 * parameters:
 *   this_nodal    <-- pointer to nodal mesh structure
 *   entity_dim    <-- dimension of entities we want to count (1 to 3)
 *   element_index --> vertex -> element connectivity index (O to n-1)
 *   element_id    --> vertex -> element connectivity (0 to n-1)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_get_vertex_elements(const fvm_nodal_t   *this_nodal,
                              int                  entity_dim,
                              cs_lnum_t          **element_index,
                              cs_lnum_t          **element_id)
{
  int section_id;
  cs_lnum_t i, j, k, l, element_concat_id, vertex_id;

  cs_lnum_t *element_count = NULL;
  cs_lnum_t *_element_index = NULL;
  cs_lnum_t *_element_id = NULL;

  const cs_lnum_t n_vertices = this_nodal->n_vertices;

  /* Initialize return values */

  *element_index = NULL;
  *element_id = NULL;

  /* Build count and index */
  /*-----------------------*/

  BFT_MALLOC(element_count, n_vertices, cs_lnum_t);
  for (i = 0; i < n_vertices; i++)
    element_count[i] = 0;

  /* Loop on sections */

  for (section_id = 0; section_id < this_nodal->n_sections; section_id++) {

    const fvm_nodal_section_t *section = this_nodal->sections[section_id];

    /* Ignore sections of different dimension than that queried */

    if (section->entity_dim != entity_dim)
      continue;

    if (section->type != FVM_FACE_POLY && section->type != FVM_CELL_POLY) {

      const int stride = section->stride;

      for (i = 0; i < section->n_elements; i++) {
        for (j = 0; j < stride; j++)
          element_count[section->vertex_num[i*stride + j] - 1] += 1;
      }
    }

    else if (section->type == FVM_FACE_POLY) {

      for (i = 0; i < section->n_elements; i++) {
        cs_lnum_t start_id = section->vertex_index[i];
        cs_lnum_t end_id = section->vertex_index[i+1];
        for (j = start_id; j < end_id; j++)
          element_count[section->vertex_num[j] - 1] += 1;
      }

    }

    else if (section->type == FVM_CELL_POLY) {

      for (i = 0; i < section->n_elements; i++) {
        cs_lnum_t start_face_id = section->face_index[i];
        cs_lnum_t end_face_id = section->face_index[i+1];
        for (j = start_face_id; j < end_face_id; j++) {
          cs_lnum_t face_id = CS_ABS(section->face_num[j]) - 1;
          cs_lnum_t start_id = section->vertex_index[face_id];
          cs_lnum_t end_id = section->vertex_index[face_id+1];
          for (k = start_id; k < end_id; k++)
            element_count[section->vertex_num[k] - 1] += 1;
        }
      }

    }

  } /* End of loop on sections */

  BFT_MALLOC(_element_index, n_vertices + 1, cs_lnum_t);
  _element_index[0] = 0;

  for (i = 0; i < n_vertices; i++) {
    _element_index[i+1] = _element_index[i] + element_count[i];
    element_count[i] = 0;
  }

  /* Build inverse connectivity */
  /*----------------------------*/

  BFT_MALLOC(_element_id, _element_index[n_vertices], cs_lnum_t);

  element_concat_id = 0;

  /* Loop on sections */

  for (section_id = 0; section_id < this_nodal->n_sections; section_id++) {

    const fvm_nodal_section_t *section = this_nodal->sections[section_id];

    /* Ignore sections of different dimension than that queried */

    if (section->entity_dim != entity_dim)
      continue;

    if (section->type != FVM_FACE_POLY && section->type != FVM_CELL_POLY) {

      const int stride = section->stride;

      for (i = 0; i < section->n_elements; i++, element_concat_id++) {
        for (j = 0; j < stride; j++) {
          vertex_id = section->vertex_num[i*stride + j] - 1;
          k = _element_index[vertex_id] + element_count[vertex_id];
          _element_id[k] = element_concat_id;
          element_count[vertex_id] += 1;
        }
      }
    }

    else if (section->type == FVM_FACE_POLY) {

      for (i = 0; i < section->n_elements; i++, element_concat_id++) {
        cs_lnum_t start_id = section->vertex_index[i];
        cs_lnum_t end_id = section->vertex_index[i+1];
        for (j = start_id; j < end_id; j++) {
          vertex_id = section->vertex_num[j] - 1;
          k = _element_index[vertex_id] + element_count[vertex_id];
          _element_id[k] = element_concat_id;
          element_count[vertex_id] += 1;
        }
      }

    }

    else if (section->type == FVM_CELL_POLY) {

      for (i = 0; i < section->n_elements; i++, element_concat_id++) {
        cs_lnum_t start_face_id = section->face_index[i];
        cs_lnum_t end_face_id = section->face_index[i+1];
        for (j = start_face_id; j < end_face_id; j++) {
          cs_lnum_t face_id = CS_ABS(section->face_num[j]) - 1;
          cs_lnum_t start_id = section->vertex_index[face_id];
          cs_lnum_t end_id = section->vertex_index[face_id+1];
          for (k = start_id; k < end_id; k++) {
            vertex_id = section->vertex_num[k] - 1;
            l = _element_index[vertex_id] + element_count[vertex_id];
            _element_id[l] = element_concat_id;
            element_count[vertex_id] += 1;
          }
        }
      }

    }

  } /* End of loop on sections */

  assert(element_concat_id == fvm_nodal_get_n_entities(this_nodal, entity_dim));

  /* Set return values */

  *element_index = _element_index;
  *element_id = _element_id;
}

/*----------------------------------------------------------------------------
 * Compute extents of a nodal mesh representation
 *
 * parameters:
 *   this_nodal   <-- pointer to mesh representation structure
 *   tolerance    <-- addition to local extents of each element:
 *                    extent = base_extent * (1 + tolerance)
 *   extents      <-> extents associated with mesh:
 *                    x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *----------------------------------------------------------------------------*/

void
fvm_nodal_extents(const fvm_nodal_t  *this_nodal,
                  double              tolerance,
                  double              extents[])
{
  int i, j;
  int dim;
  double section_extents[6];

  if (this_nodal == NULL)
    return;

  dim = this_nodal->dim;

  /* initialize extents in case mesh is empty or dim < 3 */
  for (i = 0; i < dim; i++) {
    extents[i]       =  HUGE_VAL;
    extents[i + dim] = -HUGE_VAL;
  }

  /* Compute extents */

  for (i = 0; i < this_nodal->n_sections; i++) {

    _nodal_section_extents(this_nodal->sections[i],
                           this_nodal->dim,
                           this_nodal->parent_vertex_num,
                           this_nodal->vertex_coords,
                           tolerance,
                           section_extents);

    for (j = 0; j < this_nodal->dim; j++) {
      if (section_extents[j] < extents[j])
        extents[j] = section_extents[j];
      if (section_extents[j+dim] > extents[j+dim])
        extents[j+dim] = section_extents[j+dim];
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
