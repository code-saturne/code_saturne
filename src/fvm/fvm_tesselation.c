/*============================================================================
 * Structure describing a mesh section's subdivision into simpler elements
 *
 * This is mostly useful to replace polygons or polyhedra by simpler
 * elements such as triangles, tetrahedra, and prisms upon data export.
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

/*
 * Notes:
 *
 * For a polygon with n vertices (and thus n edges), a compact way of
 * storing a tesselation is described here:
 *   The vertices of each of the polygon's triangles are defined in
 *   terms of the polgon's local vertex numbers. As the local number
 *   of resulting triangles should remain small, these numbers may
 *   be represented by bit fields. for example, with 32 bit integers,
 *   we may use bits 1-10, 11-20, and 21-30 for each of a triangle's
 *   vertices (allowing for polygons with up to 1024 vertices, which
 *   vastly exceeds the number of vertices we expect for polygons
 *   (usually in the range of 5-10).
 *   This requires n-2 values per polygon, while the expanded
 *   connectivity requires (n - 2) * 3 values.
 *
 * As polyhedra are defined by their polygonal faces, this encoding
 * is also used.
 *
 * For a quadrangle, a better way is to use a single diagonal_flag, equal to
 *   O if the quadrangle is split along its diagonal (1, 3), or
 *   1 if the quadrangle is split along its diagonal (2, 4).
 *
 * In the fvm_tesselation_t structure, the encoding[] array corresponds
 * to the opposite_vertex_num[] array (of the same size as vertex_num[])
 * for polygons and polyhedra, and to a diagonal_flag[] array of size
 * n_elements for quadrangles.
 *
 * If FVM model were ever extended to elements with more than 2 vertices
 * per edge, we would need to replace n_vertices by n_edges in some of
 * this file's function's logic, but the biggest difficulty would
 * be in adding extra vertices along new edges without duplicating them:
 * this would only require some extra processing for a given section,
 * (including some communication between processors), but would
 * involve difficulties mainly when using multiple tesselation structures
 * for multiple mesh sections, as extra vertices along section boundaries
 * would need to be shared, involving extra book-keeping (and not sharing
 * those vertices would make the corresponding elements non-conforming,
 * so it would not be an option). As in most known data models, general
 * simple polygons or polyhedra only have two vertices per edge,
 * restricting tesselation to linear element types would not be a problem
 * as regards its main intended use, which is to allow replacement of
 * polygons and/or polyhedra upon export to some formats or for use with
 * some tools which do not support these types of elements.
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

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_defs.h"
#include "fvm_io_num.h"
#include "fvm_triangulate.h"

#include "cs_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_tesselation.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/* Geometric primitives */

#undef _CROSS_PRODUCT_3D
#undef _DOT_PRODUCT_3D
#undef _MODULE_3D

#define _CROSS_PRODUCT_3D(cross_v1_v2, v1, v2) ( \
 cross_v1_v2[0] = v1[1]*v2[2] - v1[2]*v2[1],   \
 cross_v1_v2[1] = v1[2]*v2[0] - v1[0]*v2[2],   \
 cross_v1_v2[2] = v1[0]*v2[1] - v1[1]*v2[0]  )

#define _DOT_PRODUCT_3D(v1, v2) ( \
 v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])

#define _MODULE_3D(v) \
 sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

/* Tesselation encoding */

#undef _ENCODING_BITS
#undef _DECODE_TRIANGLE_VERTICES

#define _ENCODING_BITS (sizeof(fvm_tesselation_encoding_t)*8/3)
#define _DECODE_TRIANGLE_VERTICES(vertex_list, encoding, decoding_mask) (\
  vertex_list[0] = encoding & decoding_mask[0], \
  vertex_list[1] = (encoding & decoding_mask[1]) >> _ENCODING_BITS, \
  vertex_list[2] = (encoding & decoding_mask[2]) >> (_ENCODING_BITS*2) )

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef unsigned fvm_tesselation_encoding_t;

/*----------------------------------------------------------------------------
 * Structure defining a tesselation of a mesh section.
 *----------------------------------------------------------------------------*/

struct _fvm_tesselation_t {

  /* Parent section information */
  /*----------------------------*/

  fvm_element_t  type;             /* Element types */

  cs_lnum_t   n_elements;          /* Number of elements */

  int         dim;                 /* Spatial dimension */

  int         entity_dim;          /* Entity dimension */

  int         stride;              /* Element size for regular elements
                                      (0 for polygons and polyhedra) */

  cs_lnum_t   n_faces;             /* Number of faces defining polyhedra */

  /* Pointers to shared vertex coordinates */

  const cs_coord_t   *vertex_coords;    /* pointer to  vertex coordinates
                                           (always interlaced:
                                           x1, y1, z1, x2, y2, z2, ...) */

  const cs_lnum_t    *parent_vertex_num; /* Local numbers (1 to n) of local
                                            vertices in the parent mesh.
                                            This array is present only when non
                                            "trivial" (i.e. not 1, 2, ..., n). */

  /* Pointers to shared connectivity arrays */

  const cs_lnum_t   *face_index;   /* polyhedron -> faces index (O to n-1);
                                      dimension [n_elements + 1] */
  const cs_lnum_t   *face_num;     /* polyhedron -> face numbers (1 to n, signed,
                                      > 0 for outwards pointing face normal
                                      < 0 for inwards pointing face normal);
                                      dimension: [face_index[n_elements]] */

  const cs_lnum_t   *vertex_index; /* polygon face -> vertices index (O to n-1);
                                      dimension: [n_cell_faces + 1] */

  const cs_lnum_t   *vertex_num;   /* vertex numbers (1 to n);
                                      dimension: connectivity_size */

  /* Pointer to shared global element numbers */

  const fvm_io_num_t  *global_element_num;

  /* Basic information */
  /*-------------------*/

  int               n_sub_types;         /* Number of sub-element types
                                            (currently max 2) */
  fvm_element_t     sub_type[2];         /* Sub-element types */

  cs_lnum_t         n_sub_max[2];        /* Maximum number of sub-elements of
                                            each type per element */
  cs_lnum_t         n_sub_max_glob[2];   /* Global maximum number of
                                            sub-elements of each type
                                            per element */
  cs_lnum_t         n_sub[2];            /* Number of sub-elements of
                                            each type */
  cs_gnum_t         n_sub_glob[2];       /* Global number of sub-elements
                                            of each type */

  const fvm_tesselation_encoding_t  *encoding;    /* Compact tesselation
                                                     encoding array
                                                     (see notes above) */

  fvm_tesselation_encoding_t        *_encoding;   /* encoding if owner,
                                                     NULL if shared */

  const cs_lnum_t   *sub_elt_index[2];   /* index of sub-elements of each
                                            given type associated with
                                            each element (0 to n-1); */
  cs_lnum_t         *_sub_elt_index[2];  /* sub_elt_index if owner,
                                            NULL if shared */

};

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute coordinates of added vertices for a tesselation of a polyhedron.
 *
 * parameters:
 *   ts             <-- tesselation structure
 *   vertex_coords  --> coordinates of added vertices
 *   n_vertices_tot --> total number of vertex references (or NULL)
 *   element_id     <-- element index
 *----------------------------------------------------------------------------*/

static inline void
_added_vertex_coords(const fvm_tesselation_t  *ts,
                     cs_coord_t                vertex_coords[3],
                     int                      *n_vertices_tot,
                     cs_lnum_t                 element_id)
{
  cs_lnum_t n_vertices;
  cs_lnum_t i, j, k, face_id, vertex_id;

  cs_lnum_t _n_vertices_tot = 0;
  cs_coord_t cell_center[3] = {0., 0., 0.};
  cs_coord_t cell_center_divisor = 0.;

  /* Loop on polyhedron's faces */
  /*----------------------------*/

  for (i = ts->face_index[element_id];
       i < ts->face_index[element_id + 1];
       i++) {

    cs_coord_t sign, v1[3], v2[3];

    cs_coord_t vertices_center[3] = {0., 0., 0.};
    cs_coord_t face_center[3] = {0., 0., 0.};
    cs_coord_t face_normal[3] = {0., 0., 0.};
    cs_coord_t triangle_center[3] = {0., 0., 0.};
    cs_coord_t triangle_normal[3] = {0., 0., 0.};
    cs_coord_t face_surface = 0.;
    cs_coord_t triangle_surface = 0.;
    const cs_coord_t *current_coords = NULL;

    face_id = CS_ABS(ts->face_num[i]) - 1;

    n_vertices = ts->vertex_index[face_id+1] - ts->vertex_index[face_id];

    _n_vertices_tot += n_vertices;

    /* First loop on face vertices to estimate face center location */

    for (j = 0; j < n_vertices; j++) {

      vertex_id = ts->vertex_num[ts->vertex_index[face_id] + j] - 1;

      if (ts->parent_vertex_num != NULL)
        current_coords
          = ts->vertex_coords + ((ts->parent_vertex_num[vertex_id] - 1) * 3);
      else
        current_coords = ts->vertex_coords + (vertex_id * 3);

      for (k = 0; k < 3; k++)
        vertices_center[k] += current_coords[k];

    }

    for (k = 0; k < 3; k++)
      vertices_center[k] /= (cs_coord_t)n_vertices;

    /* At this stage, current_coords points
       to the last face vertice's coords */

    for (k = 0; k < 3; k++) {
      v1[k] = current_coords[k] - vertices_center[k];
      triangle_center[k] = current_coords[k] + vertices_center[k];
    }

    /* Second loop on face vertices to estimate face normal */

    for (j = 0; j < n_vertices; j++) {

      vertex_id = ts->vertex_num[ts->vertex_index[face_id] + j] - 1;

      if (ts->parent_vertex_num != NULL)
        current_coords
          = ts->vertex_coords + ((ts->parent_vertex_num[vertex_id] - 1) * 3);
      else
        current_coords = ts->vertex_coords + (vertex_id * 3);

      for (k = 0; k < 3; k++) {
        v2[k] = current_coords[k] - vertices_center[k];
        triangle_center[k] =   (triangle_center[k] + current_coords[k])
                             * (1./3.);
      }

      _CROSS_PRODUCT_3D(triangle_normal, v1, v2);

      for (k = 0; k < 3; k++)
        face_normal[k] += triangle_normal[k];

      triangle_surface = _MODULE_3D(triangle_normal) * 0.5;

      if (_DOT_PRODUCT_3D(triangle_normal, face_normal) > 0)
        sign = 1.;
      else
        sign = -1.;

      /* Add triangle contributions to face center and normal */

      face_surface += triangle_surface * sign;
      for (k = 0; k < 3; k++)
        face_center[k] += triangle_center[k] * triangle_surface * sign;

      /* Prepare for next face vertex */

      for (k = 0; k < 3; k++) {
        v1[k] = v2[k];
        triangle_center[k] = current_coords[k] + vertices_center[k];
      }

    } /* End of second loop on face vertices */

    /* If first normal was not oriented like the face normal, correct sign */

    if (face_surface < 0) {
      face_surface = - face_surface;
      for (k = 0; k < 3; k++)
        face_center[k] = -face_center[k];
    }

    /* We should now divide the face_center[] values by face_surface to
       finish, but the contribution of face center to element center is:
       face_center[] * face_surface;
       which the face_center[] array contains, so we use those directly */

    cell_center_divisor += face_surface;
    for (k = 0; k < 3; k++)
      cell_center[k] += face_center[k];

  } /* End of loop on polyhedron's faces */

  for (k = 0; k < 3; k++)
    vertex_coords[k] = cell_center[k] / cell_center_divisor;

  if (n_vertices_tot != NULL)
    *n_vertices_tot = _n_vertices_tot;
}

/*----------------------------------------------------------------------------
 * Solve Ax = B for a 4x4 system.
 *
 * parameters:
 *   a                  <-- matrix
 *   b                  <-- right hand side
 *   x                  --> solution of Ax = b
 *
 * returns: 0 if solution was found, 1 in case of error (zero pivot)
 *----------------------------------------------------------------------------*/

static int
_solve_ax_b_4(double            a[4][4],
              double  *restrict b,
              double  *restrict x)
{
  int i, j, k, k_pivot;

  double abs_pivot, abs_a_ki, factor;
  double _a[4][4], _b[4];
  double _a_swap[4], _b_swap;

  const double _epsilon = 1.0e-15;

  /* Copy array and RHS (avoid modifying a and b) */

  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      _a[i][j] = a[i][j];
    }
    _b[i] = b[i];
  }

  /* Forward elimination */

  for (i = 0; i < 4-1; i++) {

    /* Seek best pivot */

    k_pivot = i;
    abs_pivot = CS_ABS(_a[i][i]);

    for (k = i+1; k < 4; k++) {
      abs_a_ki = CS_ABS(_a[k][i]);
      if (abs_a_ki > abs_pivot) {
        abs_pivot = abs_a_ki;
        k_pivot = k;
      }
    }

    /* Swap if necessary */

    if (k_pivot != i) {

      for (j = 0; j < 4; j++) {
        _a_swap[j] = _a[i][j];
        _a[i][j] = _a[k_pivot][j];
        _a[k_pivot][j] = _a_swap[j];
      }

      _b_swap = _b[i];
      _b[i] = _b[k_pivot];
      _b[k_pivot] = _b_swap;

    }

    if (abs_pivot < _epsilon)
      return 1;

    /* Eliminate values */

    for (k = i+1; k < 4; k++) {

      factor = _a[k][i] / _a[i][i];

      _a[k][i] = 0.0;
      for (j = i+1; j < 4; j++) {
        _a[k][j] -= _a[i][j]*factor;
      }
      _b[k] -= _b[i]*factor;

    }

  }

  /* Solve system */

  x[3] =  _b[3]                                                   / _a[3][3];
  x[2] = (_b[2] - _a[2][3]*x[3])                                  / _a[2][2];
  x[1] = (_b[1] - _a[1][3]*x[3]  - _a[1][2]*x[2])                 / _a[1][1];
  x[0] = (_b[0] - _a[0][3]*x[3]  - _a[0][2]*x[2] - _a[0][1]*x[1]) / _a[0][0];

  return 0;
}

/*----------------------------------------------------------------------------
 * Compute coordinates of added vertices for a tesselation of a polyhedron.
 *
 * The heap and vertex_list arrays should be large enough to contain
 * all vertex references of all the polyhedron's faces.
 *
 * parameters:
 *   ts               <-- tesselation structure
 *   heap             --- working array
 *   vertex_list      --> list of polyhedron's vertices
 *   vertex_list_size --> list size
 *   element_id       <-- element index
 *----------------------------------------------------------------------------*/

static inline void
_polyhedron_vertices(const fvm_tesselation_t   *ts,
                     cs_lnum_t                  heap[],
                     cs_lnum_t                  vertex_list[],
                     int                       *vertex_list_size,
                     cs_lnum_t                  element_id)
{
  cs_lnum_t n_vertices;
  cs_lnum_t i, j, face_id, vertex_id;

  int heap_size = 0;

  /* Loop on polyhedron's faces */
  /*----------------------------*/

  for (i = ts->face_index[element_id];
       i < ts->face_index[element_id + 1];
       i++) {

    /* Loop on face vertices */
    /*-----------------------*/

    face_id = CS_ABS(ts->face_num[i]) - 1;

    n_vertices = ts->vertex_index[face_id+1] - ts->vertex_index[face_id];

    for (j = 0; j < n_vertices; j++) {

      /* Add to heap */

      int hole = heap_size;
      int parent = (hole - 1) / 2;

      vertex_id = ts->vertex_num[ts->vertex_index[face_id] + j] - 1;

      while (hole > 0 && heap[parent] > vertex_id) {
        heap[hole] = heap[parent];
        hole = parent;
        parent = (parent - 1) / 2;
      }

      heap[hole] = vertex_id;
      heap_size++;

    }

  }

  /* Now remove entries from heap */

  i = 0;

  if (heap_size > 0)
    vertex_list[i++] = heap[0];

  while (heap_size > 0) {

    if (heap[0] != vertex_list[i-1])
      vertex_list[i++] = heap[0];

    heap_size--;

    /* Insert the last item in the heap whose root is empty */

    {
      int start = 0;
      int child = 1; /* left child of 0 */

      vertex_id = heap[heap_size];

      while (child < heap_size) {

        if (child < (heap_size - 1) && heap[child] > heap[child+1])
          child++;  /* child is the smaller child of start */

        if (vertex_id < heap[child])
          break;    /* item stays in start */

        else {
          heap[start] = heap[child];
          start = child;
          child = 2*start + 1;
        }

      }

      heap[start] = vertex_id;

    }

  }

  *vertex_list_size = i;
}

/*----------------------------------------------------------------------------
 * Compute field values at added vertices for a tesselation of polyhedra.
 *
 * One additional vertex is added near the center of each polyhedra.
 * For element types other than polyhedra, there is no need for added
 * vertices, so this function returns immediately.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *   vertex_coords    <-- coordinates of added vertices
 *   src_dim          <-- dimension of source data
 *   src_dim_shift    <-- source data dimension shift (start index)
 *   dest_dim         <-- destination data dimension (1 if non interlaced)
 *   start_id         <-- added vertices start index
 *   end_id           <-- added vertices past the end index
 *   src_interlace    <-- indicates if source data is interlaced
 *   src_datatype     <-- source data type (float, double, or int)
 *   dest_datatype    <-- destination data type (float, double, or int)
 *   n_parent_lists   <-- number of parent lists (if parent_num != NULL)
 *   parent_num_shift <-- parent number to value array index shifts;
 *                        size: n_parent_lists
 *   parent_num       <-- if n_parent_lists > 0, parent entity numbers
 *   src_data         <-- array of source arrays (at least one, with one per
 *                        source dimension if non interlaced, times one per
 *                        parent list if multiple parent lists, with
 *                        x_parent_1, y_parent_1, ..., x_parent_2, ...) order
 *   dest_data        --> destination buffer
 *----------------------------------------------------------------------------*/

static void
_vertex_field_of_real_values(const fvm_tesselation_t  *this_tesselation,
                             int                       src_dim,
                             int                       src_dim_shift,
                             int                       dest_dim,
                             cs_lnum_t                 start_id,
                             cs_lnum_t                 end_id,
                             cs_interlace_t            src_interlace,
                             cs_datatype_t             src_datatype,
                             cs_datatype_t             dest_datatype,
                             int                       n_parent_lists,
                             const cs_lnum_t           parent_num_shift[],
                             const cs_lnum_t           parent_num[],
                             const void         *const src_data[],
                             void               *const dest_data)
{
  int n_vertices_tot, pl, vertex_list_size;
  cs_lnum_t i, j, parent_id, src_id, vertex_id;

  int _src_dim_shift = src_dim_shift;
  int max_list_size = 128;

  cs_lnum_t _heap[128];
  cs_lnum_t _vertex_list[128];
  cs_lnum_t *heap = _heap;
  cs_lnum_t *vertex_list = _vertex_list;

  const fvm_tesselation_t *ts = this_tesselation;
  const cs_coord_t *current_coords = NULL;

  if (ts->type != FVM_CELL_POLY)
    return;

  /* Main loop on polyhedra */
  /*------------------------*/

  for (i = start_id; i < end_id; i++) {

    for (j = 0; j < dest_dim; j++) { /* Loop on destination dimension */

      cs_coord_t vertex_coords[3];
      double coeff[4];
      double interpolated_value, v_x, v_y, v_z;
      double v_f = 0.;

      double a[4][4] = {{0., 0., 0., 0.},
                        {0., 0., 0., 0.},
                        {0., 0., 0., 0.},
                        {0., 0., 0., 0.}};
      double b[4] = {0., 0., 0., 0.};

      _added_vertex_coords(ts,
                           vertex_coords,
                           &n_vertices_tot,
                           i);

      /* Allocate or resize working arrays if too small */

      if (n_vertices_tot > max_list_size) {
        max_list_size *= 2;
        if (heap == _heap) {
          heap = NULL;
          vertex_list = NULL;
        }
        BFT_REALLOC(heap, max_list_size, cs_lnum_t);
        BFT_REALLOC(vertex_list, max_list_size, cs_lnum_t);
      }

      /* Obtain list of polyhedron's vertices */

      _polyhedron_vertices(ts,
                           heap,
                           vertex_list,
                           &vertex_list_size,
                           i);

      /* Build matrix for least squares */

      for (j = 0; j < vertex_list_size; j++) {

        vertex_id = vertex_list[j];

        if (ts->parent_vertex_num != NULL)
          current_coords
            = ts->vertex_coords + ((ts->parent_vertex_num[vertex_id] - 1) * 3);
        else
          current_coords = ts->vertex_coords + (vertex_id * 3);

        v_x = current_coords[0];
        v_y = current_coords[1];
        v_z = current_coords[2];

        /* Position in source data for current vertex value */

        if (n_parent_lists == 0) {
          pl = 0;
          parent_id = vertex_id;
        }
        else if (parent_num != NULL) {
          for (parent_id = parent_num[vertex_id] - 1, pl = n_parent_lists - 1;
               parent_id < parent_num_shift[pl];
               pl--);
          assert(pl > -1);
          parent_id -= parent_num_shift[pl];
        }
        else {
          for (parent_id = vertex_id, pl = n_parent_lists - 1 ;
               parent_id < parent_num_shift[pl] ;
               pl--);
          assert(pl > -1);
          parent_id -= parent_num_shift[pl];
        }

        if (src_interlace == CS_INTERLACE) {
          src_id = parent_id * src_dim + _src_dim_shift;
        }
        else {
          pl = src_dim*pl + _src_dim_shift;
          src_id = parent_id;
        }

        /* Now convert based on source datatype */

        switch(src_datatype) {
        case CS_FLOAT:
          v_f = ((const float *const *const)src_data)[pl][src_id];
          break;
        case CS_DOUBLE:
          v_f = ((const double *const *const)src_data)[pl][src_id];
          break;
        default:
          assert(0);
        }

        a[0][0] += v_x * v_x;
        a[0][1] += v_x * v_y;
        a[0][2] += v_x * v_z;
        a[0][3] += v_x;

        a[1][1] += v_y * v_y;
        a[1][2] += v_y * v_z;
        a[1][3] += v_y;

        a[2][2] += v_z * v_z;
        a[2][3] += v_z;

        a[3][3] += 1.;

        b[0] += v_x * v_f;
        b[1] += v_y * v_f;
        b[2] += v_z * v_f;
        b[3] += v_f;

      } /* End of loop on element vertices */

      /* Matrix is symmetric */

      a[1][0] = a[0][1];
      a[2][0] = a[0][2];
      a[3][0] = a[0][3];

      a[2][1] = a[1][2];
      a[3][1] = a[1][3];

      a[3][2] = a[2][3];

      /* Find hyperplane coefficients and compute added value */

      if (_solve_ax_b_4(a, b, coeff) == 0) {
        interpolated_value = (  coeff[0] * vertex_coords[0]
                              + coeff[1] * vertex_coords[1]
                              + coeff[2] * vertex_coords[2]
                              + coeff[3]);
      }
      else {
        interpolated_value = v_f; /* last encountered value */
      }

      /* Now convert based on destination datatype */

      switch(dest_datatype) {
      case CS_FLOAT:
        ((float *const)dest_data)[i] = interpolated_value;
        break;
      case CS_DOUBLE:
        ((double *const)dest_data)[i] = interpolated_value;
        break;
      default:
        assert(0);
      }

    } /* End of loop on destination dimension */

  } /* End of main loop on polyhedra */

  if (heap != _heap) {
    BFT_FREE(heap);
    BFT_FREE(vertex_list);
  }
}

/*----------------------------------------------------------------------------
 * Initialize decoding mask for polygon triangulations.
 *
 * parameters:
 *   decoding_mask         <-- decoding mask for triangle
 *----------------------------------------------------------------------------*/

static void
_init_decoding_mask(fvm_tesselation_encoding_t decoding_mask[3])
{
  unsigned i;

  for (i = 0; i < _ENCODING_BITS; i++)
    decoding_mask[0] = (decoding_mask[0] << 1) + 1;
  decoding_mask[1] = decoding_mask[0] << _ENCODING_BITS;
  decoding_mask[2] = decoding_mask[1] << _ENCODING_BITS;
}

/*----------------------------------------------------------------------------
 * Tesselate polygons (2D elements or 3D element faces) of a mesh section
 * referred to by an fvm_tesselation_t structure. Quadrangles are not
 * tesselated.
 *
 * parameters:
 *   this_tesselation   <-> partially initialized tesselation structure
 *   dim                <-- spatial dimension
 *   vertex_coords      <-- associated vertex coordinates array
 *   parent_vertex_num  <-- optional indirection to vertex coordinates
 *   error_count        --> number of triangulation errors counter (optional)
 *----------------------------------------------------------------------------*/

static void
_tesselate_polygons(fvm_tesselation_t  *this_tesselation,
                    int                 dim,
                    const cs_coord_t    vertex_coords[],
                    const cs_lnum_t     parent_vertex_num[],
                    cs_lnum_t          *error_count)
{
  int type_id;
  cs_lnum_t n_vertices, n_triangles, n_quads, n_elements;
  cs_lnum_t n_vertices_max, n_triangles_max;
  cs_lnum_t i, j, k, vertex_id, encoding_id;
  fvm_tesselation_encoding_t encoding_sub[3];

  cs_gnum_t n_g_elements_tot[2] = {0, 0}; /* Global new elements count */
  cs_lnum_t n_elements_tot[2] = {0, 0}; /* New triangles/quadrangles */
  cs_lnum_t n_g_elements_max[2] = {0, 0}; /* Global max triangles/quadrangles */
  cs_lnum_t n_elements_max[2] = {0, 0}; /* Max triangles/quadrangles */
  cs_lnum_t *triangle_vertices = NULL;
  fvm_triangulate_state_t *state = NULL;

  fvm_tesselation_t *ts = this_tesselation;

  /* Initialization */

  if (error_count != NULL)
    *error_count = 0;

  if (ts->type != FVM_CELL_POLY)
    n_elements = ts->n_elements;
  else
    n_elements = ts->n_faces;

  /* Count expected local numbers of triangles and quadrangles */

  n_vertices_max = 0;
  n_triangles_max = 0;

  if (ts->vertex_index != NULL) {
    for (i = 0 ; i < n_elements ; i++) {
      n_vertices = ts->vertex_index[i+1] - ts->vertex_index[i];
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

  if (n_vertices_max > 2<<_ENCODING_BITS)
    bft_error(__FILE__, __LINE__, 0,
              _("Encoding of tesselation impossible:\n"
                "maximum number of vertices per polygon: %u\n"
                "maximum integer encoded on %d/3 = %d bits: %ld\n"),
              (unsigned)n_triangles_max,
              (int)(sizeof(fvm_tesselation_encoding_t)*8),
              (int)_ENCODING_BITS, (long)(2<<_ENCODING_BITS));

  /* Destroy previous tesselation description if present */

  ts->encoding = NULL;
  if (ts->_encoding != NULL)
    BFT_FREE(ts->_encoding);

  /* Allocate memory and state variables */
  /*-------------------------------------*/

  if (n_vertices_max > 4) {
    BFT_MALLOC(ts->_encoding,
               ts->vertex_index[n_elements] - n_elements*2,
               fvm_tesselation_encoding_t);
    ts->encoding = ts->_encoding;
    BFT_MALLOC(triangle_vertices, (n_vertices_max - 2) * 3, cs_lnum_t);
    state = fvm_triangulate_state_create(n_vertices_max);
  }

  n_elements_tot[0] = 0; n_elements_tot[1] = 0; /* reset */

  /* Main loop on section face elements */
  /*------------------------------------*/

  for (i = 0 ; i < n_elements ; i++) {

    n_triangles = 0;
    n_quads = 0;
    n_vertices = ts->vertex_index[i+1] - ts->vertex_index[i];
    vertex_id = ts->vertex_index[i];

    /* We calculate the encoding index base from the polygon's
       connectivity index base, knowing that for a polygon
       with n vertices, we have at most n-2 triangles
       (exactly n-2 when no error occurs) */

    encoding_id = ts->vertex_index[i] - (i*2);

    if (n_vertices == 4)
      type_id = 1;
    else
      type_id = 0;

    /* If face must be subdivided */

    if (n_vertices > 4) {

      n_triangles = fvm_triangulate_polygon(dim,
                                            1,
                                            n_vertices,
                                            vertex_coords,
                                            parent_vertex_num,
                                            (  ts->vertex_num
                                             + vertex_id),
                                            FVM_TRIANGULATE_ELT_DEF,
                                            triangle_vertices,
                                            state);

      if (n_triangles != (n_vertices - 2) && error_count != NULL)
        *error_count += 1;

      /* Encode local triangle connectivity */

      for (j = 0; j < n_triangles; j++) {

        for (k = 0; k < 3; k++)
          encoding_sub[k]
            = (   ((fvm_tesselation_encoding_t)(triangle_vertices[j*3 + k] - 1))
               << (_ENCODING_BITS * k));

        ts->_encoding[encoding_id + j]
          = encoding_sub[0] | encoding_sub[1] | encoding_sub[2];

      }

      /* In case of incomplete tesselation due to errors,
         blank unused encoding values */

      for (j = n_triangles; j < (n_vertices - 2); j++)
        ts->_encoding[encoding_id + j] = 0;

      n_elements_tot[0] += n_triangles;

    }

    /* Otherwise, tesselation trivial or not necessary for this face */

    else {

      if (ts->_encoding != NULL) {
        for (j = 0; j < (n_vertices - 2); j++)
          ts->_encoding[encoding_id + j] = 0;
      }

      if (n_vertices == 4) {
        n_elements_tot[1] += 1;
        n_quads = 1;
      }

      else if (n_vertices == 3) {
        n_elements_tot[0] += 1;
        n_triangles = 1;
      }

    }

    if (n_triangles > n_elements_max[0])
      n_elements_max[0] = n_triangles;

    if (n_quads > n_elements_max[1])
      n_elements_max[1] = n_triangles;

  } /* End of loop on elements */

  /* Free memory and state variables */

  if (n_vertices_max > 4) {
    BFT_FREE(triangle_vertices);
    state = fvm_triangulate_state_destroy(state);
  }

  /* Update tesselation structure info */

  for (type_id = 0; type_id < 2; type_id++) {
    n_g_elements_tot[type_id] = n_elements_tot[type_id];
    n_g_elements_max[type_id] = n_elements_max[type_id];
  }

  cs_parall_counter(n_g_elements_tot, 2);
  cs_parall_counter_max(n_g_elements_max, 2);

  for (type_id = 0; type_id < 2; type_id++) {
    if (n_g_elements_tot[type_id] > 0) {
      if (type_id == 0)
        ts->sub_type[ts->n_sub_types] = FVM_FACE_TRIA;
      else
        ts->sub_type[ts->n_sub_types] = FVM_FACE_QUAD;
      ts->n_sub_max[ts->n_sub_types] = n_elements_max[type_id];
      ts->n_sub_max_glob[ts->n_sub_types] = n_g_elements_max[type_id];
      ts->n_sub[ts->n_sub_types] = n_elements_tot[type_id];
      ts->n_sub_glob[ts->n_sub_types] = n_elements_tot[type_id];
      ts->n_sub_types += 1;
    }
  }

}

/*----------------------------------------------------------------------------
 * Build indexes and global counts for polygons of a mesh section referred
 * to by an fvm_tesselation_t structure.
 *
 * parameters:
 *   this_tesselation   <-> partially initialized tesselation structure
 *   global_count       <-- indicates if global counts should be built
 *----------------------------------------------------------------------------*/

static void
_count_and_index_sub_polygons(fvm_tesselation_t  *this_tesselation,
                              _Bool               global_count)
{
  int sub_type_id, type_id;
  cs_lnum_t n_vertices, n_triangles, n_elements;
  cs_lnum_t i, j, encoding_id;

  cs_lnum_t *n_sub_elements[2] = {NULL, NULL};

  fvm_tesselation_t *ts = this_tesselation;

  /* Initial checks */

  if (ts->vertex_index == NULL)
    return;

  /* Initialization */

  n_elements = ts->n_elements;

  /* Count expected local numbers of triangles and quadrangles */

  for (sub_type_id = 0; sub_type_id < ts->n_sub_types; sub_type_id++) {

    type_id = -1;
    switch(ts->sub_type[sub_type_id]) {
    case FVM_FACE_TRIA:
      type_id = 0;
      break;
    case FVM_FACE_QUAD:
      type_id = 1;
      break;
    default:
      assert(0);
    }

    BFT_MALLOC(ts->_sub_elt_index[sub_type_id], n_elements + 1, cs_lnum_t);

    for (i = 0 ; i < n_elements + 1 ; i++)
      ts->_sub_elt_index[sub_type_id][i] = 0;

    /* The index array will first be used to hold n_sub_elements[],
       of size n_elements. To simplify future conversion to index,
       we shift shuch that n_sub_elements[i] corresponds to index[i+1]. */

    n_sub_elements[type_id] = ts->_sub_elt_index[sub_type_id] + 1;

  }

  /* Loop on elements to populate n_sub_elements[].
     Note that each n_sub_elements[] array has been initialized with zeroes,
     as it is mapped to a ts->sub_elt_index[] thus initialized. */

  for (i = 0 ; i < n_elements ; i++) {

    n_vertices = ts->vertex_index[i+1] - ts->vertex_index[i];
    n_triangles = 0;

    if (n_vertices == 3) {
      n_sub_elements[0][i] = 1;
      n_triangles += 1;
    }

    else { /* if (n_vertices > 3) */

      encoding_id = ts->vertex_index[i] - (i*2);

      for (j = 0; j < (n_vertices - 2); j++) {
        if (ts->encoding != NULL) {
          if (ts->encoding[encoding_id + j] != 0)
            n_triangles += 1;
        }
      }
      if (n_triangles > 0)
        n_sub_elements[0][i] = n_triangles;
      else if (n_vertices == 4)
        n_sub_elements[1][i] = 1;
      assert(n_triangles > 0 || n_vertices == 4);
    }

  }

  /* Now compute max global number if necessary */

  for (sub_type_id = 0; sub_type_id < ts->n_sub_types; sub_type_id++)
    ts->n_sub_glob[sub_type_id] = ts->n_sub[sub_type_id];

  if (global_count == true && ts->global_element_num != NULL) {

    for (sub_type_id = 0; sub_type_id < ts->n_sub_types; sub_type_id++) {
      if (ts->_sub_elt_index[sub_type_id] != NULL) {
        cs_lnum_t * _n_sub_elements = ts->_sub_elt_index[sub_type_id] + 1;
        if (n_elements == 0) _n_sub_elements = NULL;
        ts->n_sub_glob[sub_type_id]
          = fvm_io_num_global_sub_size(ts->global_element_num,
                                       _n_sub_elements);
      }
      else
        ts->n_sub_glob[sub_type_id] = 0;
    }

  }

  /* Finally, build index from n_sub_elements_glob */

  for (sub_type_id = 0; sub_type_id < ts->n_sub_types; sub_type_id++) {

    cs_lnum_t *sub_elt_index = ts->_sub_elt_index[sub_type_id];
    for (i = 0 ; i < n_elements ; i++)
      sub_elt_index[i+1] = sub_elt_index[i] + sub_elt_index[i+1];
    ts->sub_elt_index[sub_type_id] = ts->_sub_elt_index[sub_type_id];

  }

}

/*----------------------------------------------------------------------------
 * Count resulting number of tetrahedra and pyramids when tesselating
 * polyhedra of a mesh section referred to by an fvm_tesselation_t
 * structure. This requires that the mesh section's polygonal faces have
 * already been tesselated (i.e. _tesselate_polygons has been called).
 *
 * parameters:
 *   this_tesselation   <-> partially initialized tesselation structure
 *   error_count        --> number of tesselation errors counter (optional)
 *   global_count       <-- indicates if global counts should be built
 *----------------------------------------------------------------------------*/

static void
_count_and_index_sub_polyhedra(fvm_tesselation_t  *this_tesselation,
                               cs_lnum_t          *error_count,
                               _Bool               global_count)
{
  int sub_type_id, type_id;
  cs_lnum_t n_vertices, n_triangles, n_elements;
  cs_lnum_t n_tetras, n_pyrams;
  cs_lnum_t i, j, k, face_id, encoding_id;

  cs_gnum_t n_g_elements_tot[2] = {0, 0}; /* Global new elements count */
  cs_lnum_t n_elements_tot[2] = {0, 0}; /* New tetrahedra/pyramids */
  cs_lnum_t *n_sub_elements[2] = {NULL, NULL};
  cs_lnum_t n_g_elements_max[2] = {0, 0}; /* Global max tetrahedra/pyramids */
  cs_lnum_t n_elements_max[2] = {0, 0}; /* Max tetrahedra/pyramids */

  fvm_tesselation_t *ts = this_tesselation;

  /* Initial checks */

  assert(ts->vertex_index != NULL);

  /* Initialization */

  if (error_count != NULL)
    *error_count = 0;

  n_elements = ts->n_elements;

  /* Count expected total and local numbers of tetrahedra and pyramids;
     at this stage, the tesselation has been initialized as a
     polygon tesselation; adding a "center" vertex, triangle faces
     lead to tetrahedra, and quadrangles to pyramids */

  for (sub_type_id = 0; sub_type_id < ts->n_sub_types; sub_type_id++) {

    type_id = -1;
    switch(ts->sub_type[sub_type_id]) {
    case FVM_FACE_TRIA:
      type_id = 0;
      break;
    case FVM_FACE_QUAD:
      type_id = 1;
      break;
    default:
      assert(0);
    }

    BFT_MALLOC(ts->_sub_elt_index[sub_type_id], n_elements + 1, cs_lnum_t);

    for (i = 0 ; i < n_elements + 1 ; i++)
      ts->_sub_elt_index[sub_type_id][i] = 0;

    /* The index array will first be used to hold n_sub_elements[],
       of size n_elements. To simplify future conversion to index,
       we shift shuch that n_sub_elements[i] corresponds to index[i+1]. */

    n_sub_elements[type_id] = ts->_sub_elt_index[sub_type_id] + 1;

  }

  /* Counting loop on polyhedra elements */
  /*-------------------------------------*/

  for (i = 0 ; i < n_elements ; i++) {

    n_tetras = 0;
    n_pyrams = 0;

    for (j = ts->face_index[i];     /* Loop on element faces */
         j < ts->face_index[i+1];
         j++) {

      face_id = CS_ABS(ts->face_num[j]) - 1;

      n_vertices = ts->vertex_index[face_id+1] - ts->vertex_index[face_id];

      if (n_vertices == 3)
        n_tetras += 1;

      else { /* if (n_vertices > 3) */

        n_triangles = 0;

        encoding_id = ts->vertex_index[face_id] - (face_id*2);

        for (k = encoding_id; k < (encoding_id + n_vertices - 2); k++) {
          if (ts->encoding != NULL) {
            if (ts->encoding[k] != 0)
              n_triangles += 1;
          }
        }

        if (error_count != NULL && n_triangles < n_vertices - 2)
          *error_count += 1;

        if (n_triangles > 0)
          n_tetras += n_triangles;
        else if (n_vertices == 4)
          n_pyrams += 1;

      }

    } /* End of loop on element faces */

    n_elements_tot[0] += n_tetras;
    n_elements_tot[1] += n_pyrams;

    if (n_tetras > n_elements_max[0])
      n_elements_max[0] = n_tetras;

    if (n_pyrams > n_elements_max[1])
      n_elements_max[1] = n_pyrams;

    if (n_sub_elements[0] != NULL)
      n_sub_elements[0][i] = n_tetras;

    if (n_sub_elements[1] != NULL)
      n_sub_elements[1][i] = n_pyrams;

  }  /* End of loop on elements */

  /* Update tesselation structure info */

  for (type_id = 0; type_id < 2; type_id++) {
    n_g_elements_tot[type_id] = n_elements_tot[type_id];
    n_g_elements_max[type_id] = n_elements_max[type_id];
  }

  cs_parall_counter(n_g_elements_tot, 2);
  cs_parall_counter_max(n_g_elements_max, 2);

  ts->n_sub_types = 0;

  for (type_id = 0; type_id < 2; type_id++) {
    if (n_g_elements_tot[type_id] > 0) {
      if (type_id == 0)
        ts->sub_type[ts->n_sub_types] = FVM_CELL_TETRA;
      else
        ts->sub_type[ts->n_sub_types] = FVM_CELL_PYRAM;
      ts->n_sub_max[ts->n_sub_types] = n_elements_max[type_id];
      ts->n_sub_max_glob[ts->n_sub_types] = n_g_elements_max[type_id];
      ts->n_sub[ts->n_sub_types] = n_elements_tot[type_id];
      ts->n_sub_glob[ts->n_sub_types] = n_elements_tot[type_id];
      ts->n_sub_types += 1;
    }
  }

  /* Now compute max global number if necessary */

  for (sub_type_id = 0; sub_type_id < ts->n_sub_types; sub_type_id++)
    ts->n_sub_glob[sub_type_id] = ts->n_sub[sub_type_id];

  if (global_count == true && ts->global_element_num != NULL) {

    for (sub_type_id = 0; sub_type_id < ts->n_sub_types; sub_type_id++) {
      if (ts->_sub_elt_index[sub_type_id] != NULL) {
        cs_lnum_t * _n_sub_elements = ts->_sub_elt_index[sub_type_id] + 1;
        if (n_elements == 0) _n_sub_elements = NULL;
        ts->n_sub_glob[sub_type_id]
          = fvm_io_num_global_sub_size(ts->global_element_num,
                                       _n_sub_elements);
      }
      else
        ts->n_sub_glob[sub_type_id] = 0;
    }

  }

  /* Finally, build index from n_sub_elements_glob */

  for (sub_type_id = 0; sub_type_id < ts->n_sub_types; sub_type_id++) {

    cs_lnum_t *sub_elt_index = ts->_sub_elt_index[sub_type_id];
    for (i = 0 ; i < n_elements ; i++)
      sub_elt_index[i+1] = sub_elt_index[i] + sub_elt_index[i+1];
    ts->sub_elt_index[sub_type_id] = ts->_sub_elt_index[sub_type_id];

  }

}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Decode polygons tesselation to a connectivity buffer.
 *
 * parameters:
 *   this_tesselation   <-- tesselation structure
 *   connect_type       <-- destination element type
 *   global_vertex_num  <-- global vertex numbering
 *   vertex_num         --> sub-element (global) vertex connectivity
 *----------------------------------------------------------------------------*/

static void
_decode_polygons_tesselation_g(const fvm_tesselation_t  *this_tesselation,
                               fvm_element_t             connect_type,
                               const fvm_io_num_t       *global_vertex_num,
                               cs_gnum_t                 vertex_num[])
{
  cs_lnum_t n_vertices;
  cs_lnum_t i, k, l, vertex_id, encoding_id;
  cs_lnum_t tv[3];
  fvm_tesselation_encoding_t decoding_mask[3] = {0, 0, 0};

  cs_lnum_t n_sub_tot = 0;
  const fvm_tesselation_t *ts = this_tesselation;

  const cs_gnum_t *_global_vertex_num
                      = fvm_io_num_get_global_num(global_vertex_num);

  _init_decoding_mask(decoding_mask);

  /* Main loop on polygons */
  /*-----------------------*/

  for (i = 0; i < this_tesselation->n_elements; i++) {

    n_vertices = ts->vertex_index[i+1] - ts->vertex_index[i];
    vertex_id = ts->vertex_index[i];
    encoding_id = ts->vertex_index[i] - (i*2);

    /* Sub-elements (triangles) connectivity */
    /*---------------------------------------*/

    if (   connect_type == FVM_FACE_TRIA
        && n_vertices > 3 && ts->encoding != NULL) {

      /* Loop on triangles */

      for (k = 0; k < (n_vertices - 2); k++) {

        if (ts->encoding[encoding_id + k] != 0) {

          /* Fill connectivity array */

          /* Fill connectivity array */

          _DECODE_TRIANGLE_VERTICES(tv,
                                    ts->encoding[encoding_id + k],
                                    decoding_mask);

          if (_global_vertex_num != NULL) {
            for (l = 0; l < 3; l++)
              vertex_num[n_sub_tot*3 + l]
                = _global_vertex_num[ts->vertex_num[vertex_id + tv[l]] - 1];
          }
          else {
            for (l = 0; l < 3; l++)
              vertex_num[n_sub_tot*3 + l]
                = ts->vertex_num[vertex_id + tv[l]];
          }

          n_sub_tot += 1;

        }

      } /* End of loop on polygon vertices */

    }

    /* Non-tesselated elements connectivity */
    /*--------------------------------------*/

    else if (   (connect_type == FVM_FACE_TRIA && n_vertices == 3)
             || (connect_type == FVM_FACE_QUAD && n_vertices == 4)) {

      /* Check that polygon is not subdivided */

      k = n_vertices;
      if (ts->encoding != NULL) {
        for (k = 0; k < (n_vertices - 2); k++)
          if (ts->encoding[encoding_id + k] > 0)
            break;
      }

      if (k == n_vertices - 2 || ts->encoding == NULL) {

        /* Fill connectivity array */

        if (_global_vertex_num != NULL) {
          for (k = 0; k < n_vertices; k++)
            vertex_num[n_sub_tot*n_vertices +k]
              = _global_vertex_num[ts->vertex_num[vertex_id + k] - 1];
        }
        else {
          for (k = 0; k < n_vertices; k++)
            vertex_num[n_sub_tot*n_vertices +k]
              = ts->vertex_num[vertex_id + k];
        }

        n_sub_tot += 1;

      }

    } /* End of case where polygon is not tesselated */

  } /* End of loop on polygons */
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Decode polygons tesselation to a connectivity buffer.
 *
 * To avoid requiring huge buffers and computing unneeded element
 * connectivities this function may decode a partial connectivity range,
 * starting at polygon index start_id and ending either when the indicated
 * buffer size or the last polygon is attained.
 * It returns the effective polygon index end.
 *
 * parameters:
 *   this_tesselation   <-- tesselation structure
 *   connect_type       <-- destination element type
 *   start_id           <-- start index of polygons subset in parent section
 *   buffer_limit       <-- maximum number of sub-elements of destination
 *                          element type allowable for vertex_num[] buffer
 *   vertex_num         --> sub-element (global) vertex connectivity
 *
 * returns:
 *   polygon index end corresponding to decoded range
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_decode_polygons_tesselation_l(const fvm_tesselation_t  *this_tesselation,
                               fvm_element_t             connect_type,
                               cs_lnum_t                 start_id,
                               cs_lnum_t                 buffer_limit,
                               cs_lnum_t                 vertex_num[])
{
  cs_lnum_t n_vertices;
  cs_lnum_t i, j, k, vertex_id, encoding_id;
  cs_lnum_t tv[3];
  fvm_tesselation_encoding_t decoding_mask[3] = {0, 0, 0};

  cs_lnum_t n_sub_tot = 0;
  const fvm_tesselation_t *ts = this_tesselation;
  cs_lnum_t retval = start_id;

  _init_decoding_mask(decoding_mask);

  /* Main loop on polygons */
  /*-----------------------*/

  for (i = start_id ; i < this_tesselation->n_elements; i++) {

    n_vertices = ts->vertex_index[i+1] - ts->vertex_index[i];
    vertex_id = ts->vertex_index[i];
    encoding_id = ts->vertex_index[i] - (i*2);

    /* Sub-elements (triangles) connectivity */
    /*---------------------------------------*/

    if (   connect_type == FVM_FACE_TRIA
        && n_vertices > 3 && ts->encoding != NULL) {

      /* Loop on polygon vertices */

      for (j = 0; j < (n_vertices - 2); j++) {

        if (ts->encoding[encoding_id + j] != 0) {

          /* Fill connectivity array */
          /* Return previous element's end index if buffer size reached */

          if (n_sub_tot >= buffer_limit)
            return i;

          /* Fill connectivity array */

          _DECODE_TRIANGLE_VERTICES(tv,
                                    ts->encoding[encoding_id + j],
                                    decoding_mask);

          for (k = 0; k < 3; k++)
            vertex_num[n_sub_tot*3 + k] = ts->vertex_num[vertex_id + tv[k]];

          n_sub_tot += 1;

        }

      } /* End of loop on polygon vertices */

    }

    /* Non-tesselated elements connectivity */
    /*--------------------------------------*/

    else if (   (connect_type == FVM_FACE_TRIA && n_vertices == 3)
             || (connect_type == FVM_FACE_QUAD && n_vertices == 4)) {

      /* Check that polygon is not subdivided */

      j = n_vertices;
      if (ts->encoding != NULL) {
        for (j = 0; j < (n_vertices - 2); j++)
          if (ts->encoding[encoding_id + j] != 0)
            break;
      }

      if (j == n_vertices - 2 || ts->encoding == NULL) {

        /* Return previous element's end index if buffer size reached */

        if (n_sub_tot >= buffer_limit)
          return i;

        /* Fill connectivity array */

        for (j = 0; j < n_vertices; j++)
          vertex_num[n_sub_tot*n_vertices +j] = ts->vertex_num[vertex_id + j];

        n_sub_tot += 1;

      }

    } /* End of case where polygon is not tesselated */

    retval = i+1; /* If end of buffer has not already been reached */

  } /* End of loop on polygons */

  return retval;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Decode polyhedra tesselation to a connectivity buffer.
 *
 * parameters:
 *   this_tesselation      <-- tesselation structure
 *   connect_type          <-- destination element type
 *   extra_vertex_base     <-- starting number for added vertices
 *   global_vertex_num     <-- global vertex numbering
 *   vertex_num            --> sub-element (global) vertex connectivity
 *----------------------------------------------------------------------------*/

static void
_decode_polyhedra_tesselation_g(const fvm_tesselation_t  *this_tesselation,
                                fvm_element_t             connect_type,
                                cs_gnum_t                 extra_vertex_base_num,
                                const fvm_io_num_t       *global_vertex_num,
                                cs_gnum_t                 vertex_num[])
{
  int orient;
  cs_lnum_t n_vertices;
  cs_lnum_t i, k, l, m, base_dest_id, face_id, vertex_id, encoding_id;
  cs_lnum_t tv[3];
  fvm_tesselation_encoding_t decoding_mask[3] = {0, 0, 0};

  cs_lnum_t n_sub_tot = 0;
  const fvm_tesselation_t *ts = this_tesselation;

  const cs_gnum_t *_global_vertex_num
                      = fvm_io_num_get_global_num(global_vertex_num);
  const cs_gnum_t *global_element_num
    = fvm_io_num_get_global_num(ts->global_element_num);

  _init_decoding_mask(decoding_mask);

  /* Main loop on polyhedra */
  /*------------------------*/

  for (i = 0; i < this_tesselation->n_elements; i++) {

    for (k = ts->face_index[i];     /* Loop on element faces */
         k < ts->face_index[i+1];
         k++) {

      face_id = CS_ABS(ts->face_num[k]) - 1;

      n_vertices = ts->vertex_index[face_id+1] - ts->vertex_index[face_id];
      vertex_id = ts->vertex_index[face_id];
      encoding_id = ts->vertex_index[face_id] - (face_id*2);

      /* Invert base orientation as it points outwards in polyhedron
         definition, but inwards in tetrahedron and pyramid reference
         connectivity */

      if (ts->face_num[k] > 0)
        orient = -1;
      else
        orient = 1;

      /* Sub-elements (tetrahedra) connectivity */
      /*----------------------------------------*/

      if (   connect_type == FVM_CELL_TETRA
          && n_vertices > 3 && ts->encoding != NULL) {

        /* Loop on face vertices */

        for (l = 0; l < (n_vertices - 2); l++) {

          if (ts->encoding[encoding_id + l] != 0) {

            if (orient == 1)
              base_dest_id = n_sub_tot*4;
            else
              base_dest_id = n_sub_tot*4 + 2;

            /* Fill connectivity array */

            _DECODE_TRIANGLE_VERTICES(tv,
                                      ts->encoding[encoding_id + l],
                                      decoding_mask);

            if (_global_vertex_num != NULL) {
              for (m = 0; m < 3; m++)
                vertex_num[base_dest_id + orient*m]
                  = _global_vertex_num[ts->vertex_num[vertex_id + tv[m]] - 1];
            }
            else {
              for (m = 0; m < 3; m++)
                vertex_num[base_dest_id + orient*m]
                  = ts->vertex_num[vertex_id + tv[m]];
            }

            /* Add extra "top" vertex */

            if (global_element_num != NULL)
              vertex_num[n_sub_tot*4 + 3] =   extra_vertex_base_num
                                            + global_element_num[i] - 1;
            else
              vertex_num[n_sub_tot*4 + 3] = extra_vertex_base_num + i;

            n_sub_tot += 1;

          }

        }

      }

      /* Non-tesselated faces connectivity */
      /*-----------------------------------*/

      else if (   (connect_type == FVM_CELL_TETRA && n_vertices == 3)
               || (connect_type == FVM_CELL_PYRAM && n_vertices == 4)) {

        /* Check that face is not subdivided */

        l = n_vertices;
        if (ts->encoding != NULL) {
          for (l = 0; l < (n_vertices - 2); l++)
            if (ts->encoding[encoding_id + l] != 0)
              break;
        }

        if (l == n_vertices - 2 || ts->encoding == NULL) {

          cs_lnum_t stride = n_vertices + 1;

          if (orient == 1)
            base_dest_id = n_sub_tot*stride;
          else
            base_dest_id = n_sub_tot*stride + n_vertices-1;

          /* Fill connectivity array */

          if (_global_vertex_num != NULL) {
            for (l = 0; l < n_vertices; l++)
              vertex_num[base_dest_id + l*orient]
                = _global_vertex_num[ts->vertex_num[vertex_id + l] - 1];
          }
          else {
            for (l = 0; l < n_vertices; l++)
              vertex_num[base_dest_id + l*orient]
                = ts->vertex_num[vertex_id + l];
          }

          /* Add extra "top" vertex */

          if (global_element_num != NULL)
            vertex_num[n_sub_tot*stride + n_vertices]
              =   extra_vertex_base_num + global_element_num[i] - 1;
          else
            vertex_num[n_sub_tot*stride + n_vertices]
              = extra_vertex_base_num + i;

          n_sub_tot += 1;

        }

      } /* End of case where face is not tesselated */

    } /* End of loop on element faces */

  } /* End of main loop on polyhedra */
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Decode polyhedra tesselation to a connectivity buffer.
 *
 * To avoid requiring huge buffers, this function may decode a partial
 * connectivity range, starting at polyhedron index start_id and ending either
 * when the indicated buffer size or the last polyhedra is attained.
 * It returns the effective polyhedron index end.
 *
 * parameters:
 *   this_tesselation      <-- tesselation structure
 *   connect_type          <-- destination element type
 *   start_id              <-- start index of polyhedra subset in parent section
 *   buffer_limit          <-- maximum number of sub-elements of destination
 *                             element type allowable for vertex_num[] buffer
 *   extra_vertex_base     <-- starting number for added vertices
 *   vertex_num            --> sub-element (global) vertex connectivity
 *
 * returns:
 *   polyhedron index end corresponding to decoded range
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_decode_polyhedra_tesselation_l(const fvm_tesselation_t  *this_tesselation,
                                fvm_element_t             connect_type,
                                cs_lnum_t                 start_id,
                                cs_lnum_t                 buffer_limit,
                                cs_lnum_t                 extra_vertex_base_num,
                                cs_lnum_t                 vertex_num[])
{
  int orient;
  cs_lnum_t n_vertices;
  cs_lnum_t i, j, k, l, m, base_dest_id, face_id, vertex_id, encoding_id;
  cs_lnum_t tv[3];
  fvm_tesselation_encoding_t decoding_mask[3] = {0, 0, 0};

  cs_lnum_t n_sub_tot = 0;
  const fvm_tesselation_t *ts = this_tesselation;
  cs_lnum_t retval = start_id;

  _init_decoding_mask(decoding_mask);

  /* Main loop on polyhedra */
  /*------------------------*/

  for (i = 0, j = start_id ;
       j < this_tesselation->n_elements;
       i++, j++) {

    for (k = ts->face_index[j];     /* Loop on element faces */
         k < ts->face_index[j+1];
         k++) {

      face_id = CS_ABS(ts->face_num[k]) - 1;

      n_vertices = ts->vertex_index[face_id+1] - ts->vertex_index[face_id];
      vertex_id = ts->vertex_index[face_id];
      encoding_id = ts->vertex_index[face_id] - (face_id*2);

      /* Invert base orientation as it points outwards in polyhedron
         definition, but inwards in tetrahedron and pyramid reference
         connectivity */

      if (ts->face_num[k] > 0)
        orient = -1;
      else
        orient = 1;

      /* Sub-elements (tetrahedra) connectivity */
      /*----------------------------------------*/

      if (   connect_type == FVM_CELL_TETRA
          && n_vertices > 3 && ts->encoding != NULL) {

        /* Loop on face vertices */

        for (l = 0; l < (n_vertices - 2); l++) {

          if (ts->encoding[encoding_id + l] != 0) {

            /* Return previous element's end index if buffer size reached */

            if (n_sub_tot >= buffer_limit)
              return j;

            if (orient == 1)
              base_dest_id = n_sub_tot*4;
            else
              base_dest_id = n_sub_tot*4 + 2;

            /* Fill connectivity array */

            _DECODE_TRIANGLE_VERTICES(tv,
                                      ts->encoding[encoding_id + l],
                                      decoding_mask);

            for (m = 0; m < 3; m++)
              vertex_num[base_dest_id + orient*m]
                = ts->vertex_num[vertex_id + tv[m]];

            /* Add extra "top" vertex */

            vertex_num[n_sub_tot*4 + 3] = extra_vertex_base_num + j;

            n_sub_tot += 1;

          }

        }

      }

      /* Non-tesselated faces connectivity */
      /*-----------------------------------*/

      else if (   (connect_type == FVM_CELL_TETRA && n_vertices == 3)
               || (connect_type == FVM_CELL_PYRAM && n_vertices == 4)) {

        /* Check that face is not subdivided */

        l = n_vertices;
        if (ts->encoding != NULL) {
          for (l = 0; l < (n_vertices - 2); l++)
            if (ts->encoding[encoding_id + l] != 0)
              break;
        }

        if (l == n_vertices - 2 || ts->encoding == NULL) {

          cs_lnum_t stride = n_vertices + 1;

          /* Return previous element's end index if buffer size reached */

          if (n_sub_tot >= buffer_limit)
            return j;

          if (orient == 1)
            base_dest_id = n_sub_tot*stride;
          else
            base_dest_id = n_sub_tot*stride + n_vertices-1;

          /* Fill connectivity array */

          for (l = 0; l < n_vertices; l++)
            vertex_num[base_dest_id + l*orient]
              = ts->vertex_num[vertex_id + l];

          /* Add extra "top" vertex */

          vertex_num[n_sub_tot*stride + n_vertices] = extra_vertex_base_num + j;

          n_sub_tot += 1;

        }

      } /* End of case where face is not tesselated */

    } /* End of loop on element faces */

    retval = j+1; /* If end of buffer has not already been reached */

  } /* End of main loop on polyhedra */

  return retval;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation of a mesh section's subdivision into simpler elements.
 *
 * The structure contains pointers to the mesh section's connectivity,
 * (passed as arguments), which is not copied. This structure should thus
 * always be destroyed before the mesh section to which it relates.
 *
 * Unused connectivity array arguments should be set to NULL (such as
 * face_index[] and face_num[] for 2D or regular (strided) elements,
 * and vertex_index[] for strided elements.
 *
 * At this stage, the structure does not yet contain tesselation information.
 *
 * parameters:
 *   element_type       <-- type of elements considered
 *   n_elements         <-- number of elements
 *   face_index         <-- polyhedron -> faces index (O to n-1)
 *                          dimension [n_elements + 1]
 *   face_num           <-- element -> face numbers (1 to n, signed,
 *                           > 0 for outwards pointing face normal
 *                           < 0 for inwards pointing face normal);
 *                          dimension: [face_index[n_elements]], or NULL
 *   vertex_index       <-- element face -> vertices index (O to n-1);
 *                          dimension: [n_cell_faces + 1], [n_elements + 1],
 *                          or NULL depending on face_index and vertex_index
 *   vertex_num         <-- element -> vertex connectivity (1 to n)
 *   global_element_num <-- global element numbers (NULL in serial mode)
 *
 * returns:
 *  pointer to created mesh section tesselation structure
 *----------------------------------------------------------------------------*/

fvm_tesselation_t *
fvm_tesselation_create(fvm_element_t        element_type,
                       cs_lnum_t            n_elements,
                       const cs_lnum_t      face_index[],
                       const cs_lnum_t      face_num[],
                       const cs_lnum_t      vertex_index[],
                       const cs_lnum_t      vertex_num[],
                       const fvm_io_num_t  *global_element_num)

{
  int  i;
  int  entity_dim = 0, stride = 0;
  fvm_tesselation_t  *this_tesselation;

  /* First, check if element type is handled and initialize
     additionnal info */

  switch (element_type) {
  case FVM_FACE_QUAD:
    entity_dim = 2;
    stride = 4;
    break;
  case FVM_FACE_POLY:
    entity_dim = 2;
    stride = 0;
    break;
  case FVM_CELL_POLY:
    entity_dim = 3;
    stride = 0;
    break;
  default:
    return NULL;
  }

  /* Now, create structure */

  BFT_MALLOC(this_tesselation, 1, fvm_tesselation_t);

  /* Assign parent mesh section info */

  this_tesselation->type = element_type;
  this_tesselation->n_elements = n_elements;
  this_tesselation->dim = 0;
  this_tesselation->entity_dim = entity_dim;

  this_tesselation->stride = stride;
  this_tesselation->n_faces = 0;

  this_tesselation->vertex_coords = NULL;
  this_tesselation->parent_vertex_num = NULL;

  this_tesselation->face_index = face_index;
  this_tesselation->face_num = face_num;
  this_tesselation->vertex_index = vertex_index;
  this_tesselation->vertex_num = vertex_num;

  this_tesselation->global_element_num = global_element_num;

  /* Check argument consistency */

  if (face_index != NULL || face_num != NULL) {
    if (element_type != FVM_CELL_POLY)
      bft_error(__FILE__, __LINE__, 0,
                _("Incoherent connectivity for tesselation:\n"
                  "Connectivity face_index or face_num non NULL,\n"
                  "but element type != FVM_CELL_POLY"));
  }

  else if (vertex_index != NULL) {
    if (element_type != FVM_FACE_POLY)
      bft_error(__FILE__, __LINE__, 0,
                _("Incoherent connectivity for tesselation:\n"
                  "Connectivy vertex_index non NULL,\n"
                  "but element type != FVM_FACE_POLY"));
  }

  /* Compute number of polyhedron faces */

  if (n_elements > 0 && face_index != NULL) {
    cs_lnum_t   j, k, face_id;
    cs_lnum_t   max_face_id = 0;
    for (j = 0; j < n_elements; j++) {
      for (k = face_index[j]; k < face_index[j+1]; k++) {
        face_id = CS_ABS(face_num[k]) - 1;
        if (face_id > max_face_id)
          max_face_id = face_id;
      }
    }
    this_tesselation->n_faces = max_face_id + 1;
  }

  /* Set tesselation structures to zero */

  this_tesselation->n_sub_types = 0;

  for (i = 0; i < FVM_TESSELATION_N_SUB_TYPES_MAX; i++) {
    this_tesselation->n_sub_max[i] = 0;
    this_tesselation->n_sub_max_glob[i] = 0;
    this_tesselation->n_sub[i] = 0;
    this_tesselation->n_sub_glob[i] =0;
    this_tesselation->sub_type[i] = FVM_N_ELEMENT_TYPES;
  }

  this_tesselation->encoding = NULL;
  this_tesselation->_encoding = NULL;

  for (i = 0; i < FVM_TESSELATION_N_SUB_TYPES_MAX; i++) {
    this_tesselation->sub_elt_index[i] = NULL;
    this_tesselation->_sub_elt_index[i] = NULL;
  }

  return (this_tesselation);
}

/*----------------------------------------------------------------------------
 * Destruction of a nodal mesh polygon splitting representation structure.
 *
 * parameters:
 *   this_tesselation <-> pointer to structure that should be destroyed
 *
 * returns:
 *  NULL pointer
 *----------------------------------------------------------------------------*/

fvm_tesselation_t *
fvm_tesselation_destroy(fvm_tesselation_t  * this_tesselation)
{
  int i;

  if (this_tesselation->_encoding != NULL)
    BFT_FREE(this_tesselation->_encoding);

  for (i = 0; i < this_tesselation->n_sub_types; i++) {
    if (this_tesselation->_sub_elt_index[i] != NULL)
      BFT_FREE(this_tesselation->_sub_elt_index[i]);
  }
  BFT_FREE(this_tesselation);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Tesselate a mesh section referred to by an fvm_tesselation_t structure.
 *
 * parameters:
 *   this_tesselation   <-> partially initialized tesselation structure
 *   dim                <-- spatial dimension
 *   vertex_coords      <-- associated vertex coordinates array
 *   parent_vertex_num  <-- optional indirection to vertex coordinates
 *   error_count        --> number of elements with a tesselation error
 *                          counter (optional)
 *----------------------------------------------------------------------------*/

void
fvm_tesselation_init(fvm_tesselation_t  *this_tesselation,
                     int                 dim,
                     const cs_coord_t    vertex_coords[],
                     const cs_lnum_t     parent_vertex_num[],
                     cs_lnum_t          *error_count)
{
  assert(this_tesselation != NULL);

  this_tesselation->dim = dim;

  this_tesselation->vertex_coords = vertex_coords;
  this_tesselation->parent_vertex_num = parent_vertex_num;

  switch(this_tesselation->type) {

  case FVM_CELL_POLY:
    _tesselate_polygons(this_tesselation,
                        dim,
                        vertex_coords,
                        parent_vertex_num,
                        error_count);
    _count_and_index_sub_polyhedra(this_tesselation,
                                   error_count,
                                   true);
    break;

  case FVM_FACE_POLY:
    _tesselate_polygons(this_tesselation,
                        dim,
                        vertex_coords,
                        parent_vertex_num,
                        error_count);
    _count_and_index_sub_polygons(this_tesselation,
                                  true);
    break;

  default:
      bft_error(__FILE__, __LINE__, 0,
                _("Tesselation of element type %s not implemented."),
                fvm_elements_type_name[this_tesselation->type]);

  }
}

/*----------------------------------------------------------------------------
 * Reduction of a nodal mesh polygon splitting representation structure;
 * only the associations (numberings) necessary to redistribution of fields
 * for output are conserved, the full connectivity being no longer useful
 * once it has been output.
 *
 * parameters:
 *   this_tesselation <-> pointer to structure that should be reduced
 *----------------------------------------------------------------------------*/

void
fvm_tesselation_reduce(fvm_tesselation_t  * this_tesselation)
{
  this_tesselation->stride = 0;
  this_tesselation->n_faces = 0;

  if (this_tesselation->face_index == NULL) {
    this_tesselation->face_num = NULL;
    this_tesselation->vertex_index = NULL;
    this_tesselation->vertex_num = NULL;
  }

  this_tesselation->encoding = NULL;
  if (this_tesselation->_encoding != NULL)
    BFT_FREE(this_tesselation->_encoding);
}

/*----------------------------------------------------------------------------
 * Return number of parent elements of a tesselation.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *
 * returns:
 *   number of parent elements
 *----------------------------------------------------------------------------*/

cs_lnum_t
fvm_tesselation_n_elements(const fvm_tesselation_t  *this_tesselation)
{
  cs_lnum_t retval = 0;

  if (this_tesselation != NULL)
    retval = this_tesselation->n_elements;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return global number of added vertices associated with a tesselation.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *
 * returns:
 *   global number of added vertices associated with the tesselation
 *----------------------------------------------------------------------------*/

cs_gnum_t
fvm_tesselation_n_g_vertices_add(const fvm_tesselation_t  *this_tesselation)
{
  cs_gnum_t retval = 0;

  assert(this_tesselation != NULL);

  if (this_tesselation->type == FVM_CELL_POLY) {

    if (this_tesselation->global_element_num != NULL)
      retval = fvm_io_num_get_global_count(this_tesselation->global_element_num);
    else
      retval = this_tesselation->n_elements;

  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Return (local) number of added vertices associated with a tesselation.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *
 * returns:
 *   global number of added vertices associated with the tesselation
 *----------------------------------------------------------------------------*/

cs_lnum_t
fvm_tesselation_n_vertices_add(const fvm_tesselation_t  *this_tesselation)
{
  cs_gnum_t retval = 0;

  assert(this_tesselation != NULL);

  if (this_tesselation->type == FVM_CELL_POLY)
    retval = this_tesselation->n_elements;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return number of resulting sub-types of a tesselation.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *
 * returns:
 *   number of resulting sub-types of the tesselation
 *----------------------------------------------------------------------------*/

int
fvm_tesselation_n_sub_types(const fvm_tesselation_t  *this_tesselation)
{
  int retval = 0;

  if (this_tesselation != NULL)
    retval = this_tesselation->n_sub_types;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return given sub-types of a tesselation.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *   sub_type_id      <-- index of sub-type in tesselation (0 to n-1)
 *
 * returns:
 *   sub-types of the tesselation with the given index
 *----------------------------------------------------------------------------*/

fvm_element_t
fvm_tesselation_sub_type(const fvm_tesselation_t  *this_tesselation,
                         int                       sub_type_id)
{
  fvm_element_t retval = FVM_N_ELEMENT_TYPES;

  if (this_tesselation == NULL)
    retval = FVM_N_ELEMENT_TYPES;
  else {
    assert(sub_type_id < this_tesselation->n_sub_types);
    retval = this_tesselation->sub_type[sub_type_id];
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Return number of elements of a given sub-type of a tesselation.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *   sub_type_id      <-- index of sub-type in tesselation (0 to n-1)
 *
 * returns:
 *   sub-types of the tesselation with the given index
 *----------------------------------------------------------------------------*/

cs_lnum_t
fvm_tesselation_n_sub_elements(const fvm_tesselation_t  *this_tesselation,
                               fvm_element_t             sub_type)
{
  int id;

  cs_lnum_t retval = 0;

  if (this_tesselation != NULL) {
    for (id = 0; id < this_tesselation->n_sub_types; id++) {
      if (this_tesselation->sub_type[id] == sub_type) {
        retval = this_tesselation->n_sub[id];
        break;
      }
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Obtain the global and maximum number of elements of a given sub-type
 * of a tesselation.
 *
 * parameters:
 *   this_tesselation    <-- tesselation structure
 *   sub_type_id         <-- index of sub-type in tesselation (0 to n-1)
 *   n_sub_elements_glob --> global number of sub-elements of the given type
 *   n_sub_elements_max  --> maximum number of sub-elements per element
 *                           of the given type (for all ranks)
 *----------------------------------------------------------------------------*/

void
fvm_tesselation_get_global_size(const fvm_tesselation_t  *this_tesselation,
                                fvm_element_t             sub_type,
                                cs_gnum_t                *n_sub_elements_glob,
                                cs_lnum_t                *n_sub_elements_max)
{
  int id;

  if (n_sub_elements_max != NULL)
    *n_sub_elements_max = 0;

  if (n_sub_elements_glob != NULL)
    *n_sub_elements_glob = 0;

  if (this_tesselation != NULL) {
    for (id = 0; id < this_tesselation->n_sub_types; id++) {
      if (this_tesselation->sub_type[id] == sub_type) {
        if (n_sub_elements_max != NULL)
          *n_sub_elements_max = this_tesselation->n_sub_max_glob[id];
        if (n_sub_elements_glob != NULL)
          *n_sub_elements_glob = this_tesselation->n_sub_glob[id];
        break;
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Return global numbering of added vertices associated with a tesselation.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *
 * returns:
 *   pointer to global numbering of added vertices for this tesselation,
 *   or NULL if no added vertices are present.
 *----------------------------------------------------------------------------*/

const fvm_io_num_t *
fvm_tesselation_global_vertex_num(const fvm_tesselation_t  *this_tesselation)
{
  const fvm_io_num_t *retval = NULL;

  assert(this_tesselation != NULL);

  if (this_tesselation->type == FVM_CELL_POLY)
    retval = this_tesselation->global_element_num;

  return retval;
}

/*----------------------------------------------------------------------------
 * Compute coordinates of added vertices for a tesselation of polyhedra.
 *
 * One additional vertex is added near the center of each polyhedra.
 * For element types other than polyhedra, there is no need for added
 * vertices, so this function returns immediately.
 *
 * parameters:
 *   this_tesselation   <-- tesselation structure
 *   vertex_coords      --> coordinates of added vertices
 *----------------------------------------------------------------------------*/

void
fvm_tesselation_vertex_coords(const fvm_tesselation_t  *this_tesselation,
                              cs_coord_t                vertex_coords[])
{
  cs_lnum_t i;

  if (this_tesselation->type != FVM_CELL_POLY)
    return;

  for (i = 0; i < this_tesselation->n_elements; i++) {

    _added_vertex_coords(this_tesselation, vertex_coords + i*3, NULL, i);

  }

}

/*----------------------------------------------------------------------------
 * Return index of sub-elements associated with each element of a given
 * sub-type of a tesselation.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *   sub_type_id      <-- index of sub-type in tesselation (0 to n-1)
 *
 * returns:
 *   index of sub-elements associated with each element (0 to n-1 numbering)
 *----------------------------------------------------------------------------*/

const cs_lnum_t *
fvm_tesselation_sub_elt_index(const fvm_tesselation_t  *this_tesselation,
                              fvm_element_t             sub_type)
{
  int id;
  const cs_lnum_t *retval = NULL;

  if (this_tesselation != NULL) {
    for (id = 0; id < this_tesselation->n_sub_types; id++) {
      if (this_tesselation->sub_type[id] == sub_type) {
        retval = this_tesselation->sub_elt_index[id];
        break;
      }
    }
  }

  return retval;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Decode tesselation to a parent element->sub elements index and
 * connectivity buffer.
 *
 * parameters:
 *   this_tesselation   <-- tesselation structure
 *   connect_type       <-- destination element type
 *   extra_vertex_base  <-- starting number for added vertices
 *   global_vertex_num  <-- global vertex numbering
 *   vertex_num         --> sub-element (global) vertex connectivity
 *----------------------------------------------------------------------------*/

void
fvm_tesselation_decode_g(const fvm_tesselation_t  *this_tesselation,
                         fvm_element_t             connect_type,
                         const fvm_io_num_t       *global_vertex_num,
                         cs_gnum_t                 extra_vertex_base,
                         cs_gnum_t                 vertex_num[])
{
  switch(this_tesselation->type) {

  case FVM_CELL_POLY:
    _decode_polyhedra_tesselation_g(this_tesselation,
                                    connect_type,
                                    extra_vertex_base,
                                    global_vertex_num,
                                    vertex_num);
    break;

  case FVM_FACE_POLY:
    _decode_polygons_tesselation_g(this_tesselation,
                                   connect_type,
                                   global_vertex_num,
                                   vertex_num);
    break;

  default:
    assert(0);

  }
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Decode tesselation to a connectivity buffer.
 *
 * To avoid requiring huge buffers and computing unneeded element
 * connectivities, this function may decode a partial connectivity range,
 * starting at polygon index start_id and ending either when the indicated
 * buffer size or the last polygon is attained.
 * It returns the effective polygon index end.
 *
 * parameters:
 *   this_tesselation   <-- tesselation structure
 *   connect_type       <-- destination element type
 *   start_id           <-- start index of polygons subset in parent section
 *   buffer_limit       <-- maximum number of sub-elements of destination
 *                          element type allowable for sub_element_idx[]
 *                          and vertex_num[] buffers
 *   extra_vertex_base  <-- starting number for added vertices
 *   vertex_num         --> sub-element (global) vertex connectivity
 *
 * returns:
 *   polygon index corresponding to end of decoded range
 *----------------------------------------------------------------------------*/

cs_lnum_t
fvm_tesselation_decode(const fvm_tesselation_t  *this_tesselation,
                       fvm_element_t             connect_type,
                       cs_lnum_t                 start_id,
                       cs_lnum_t                 buffer_limit,
                       cs_lnum_t                 extra_vertex_base,
                       cs_lnum_t                 vertex_num[])
{
  cs_lnum_t retval = 0;

  switch(this_tesselation->type) {

  case FVM_CELL_POLY:
    retval = _decode_polyhedra_tesselation_l(this_tesselation,
                                             connect_type,
                                             start_id,
                                             buffer_limit,
                                             extra_vertex_base,
                                             vertex_num);
    break;

  case FVM_FACE_POLY:
    retval = _decode_polygons_tesselation_l(this_tesselation,
                                            connect_type,
                                            start_id,
                                            buffer_limit,
                                            vertex_num);
    break;

  default:
    assert(0);

  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Distribute "per element" data from the base elements to their tesselation.
 *
 * The same data array is used for input and output, so as to avoid requiring
 * excess allocation in typical use cases (extracting data from a parent mesh
 * to a buffer and distributing it as per its tesselation).
 * The data array should be at least of size:
 * [sub_elt_index[end_id] - sub_elt_index[start_id] * size
 *
 * parameters:
 *   this_tesselation   <-- tesselation structure
 *   connect_type       <-- destination element type
 *   start_id           <-- start index of elements subset in parent section
 *   end_id             <-- end index of elements subset in parent section
 *   size               <-- data size for each element (sizeof(type)*stride)
 *   data               <-> undistributed data in, distributed data out
 *----------------------------------------------------------------------------*/

#if defined(__INTEL_COMPILER) && defined(__KNC__)
#pragma optimization_level 1 /* Crash with O2 on KNC with icc 14.0.0 20130728 */
#endif

void
fvm_tesselation_distribute(const fvm_tesselation_t  *this_tesselation,
                           fvm_element_t             connect_type,
                           cs_lnum_t                 start_id,
                           cs_lnum_t                 end_id,
                           size_t                    size,
                           void                     *data)
{
  int id;
  cs_lnum_t   i, j, k, n_sub;
  size_t  l;
  char  *src, *dest;

  const cs_lnum_t *sub_elt_index = NULL;

  /* Find index, or return */

  if (this_tesselation == NULL)
    return;

  for (id = 0; id < this_tesselation->n_sub_types; id++) {
    if (this_tesselation->sub_type[id] == connect_type) {
      sub_elt_index = this_tesselation->sub_elt_index[id];
      break;
    }
  }
  if (id == this_tesselation->n_sub_types)
    return;

  /* Distribute data (starting from the end so as not to overwrite
     data at the beginning of the array) */

  for (i = end_id, j = end_id - start_id - 1; i > start_id; i--, j--) {

    src = ((char *)data) + j*size;
    dest = ((char *)data) + (sub_elt_index[i-1] - sub_elt_index[start_id])*size;
    n_sub = sub_elt_index[i] - sub_elt_index[i-1];

    for (k = 0; k < n_sub; k++) {
      for (l = 0; l < size; l++)
        dest[k*size + l] = src[l];
    }
  }

}

/*----------------------------------------------------------------------------
 * Compute field values at added vertices for a tesselation of polyhedra.
 *
 * One additional vertex is added near the center of each polyhedra.
 * For element types other than polyhedra, there is no need for added
 * vertices, so this function returns immediately.
 *
 * parameters:
 *   this_tesselation <-- tesselation structure
 *   vertex_coords    <-- coordinates of added vertices
 *   src_dim          <-- dimension of source data
 *   src_dim_shift    <-- source data dimension shift (start index)
 *   dest_dim         <-- destination data dimension (1 if non interlaced)
 *   start_id         <-- added vertices start index
 *   end_id           <-- added vertices past the end index
 *   src_interlace    <-- indicates if source data is interlaced
 *   src_datatype     <-- source data type (float, double, or int)
 *   dest_datatype    <-- destination data type (float, double, or int)
 *   n_parent_lists   <-- number of parent lists (if parent_num != NULL)
 *   parent_num_shift <-- parent number to value array index shifts;
 *                        size: n_parent_lists
 *   parent_num       <-- if n_parent_lists > 0, parent entity numbers
 *   src_data         <-- array of source arrays (at least one, with one per
 *                        source dimension if non interlaced, times one per
 *                        parent list if multiple parent lists, with
 *                        x_parent_1, y_parent_1, ..., x_parent_2, ...) order
 *   dest_data        --> destination buffer
 *----------------------------------------------------------------------------*/

void
fvm_tesselation_vertex_values(const fvm_tesselation_t  *this_tesselation,
                              int                       src_dim,
                              int                       src_dim_shift,
                              int                       dest_dim,
                              cs_lnum_t                 start_id,
                              cs_lnum_t                 end_id,
                              cs_interlace_t            src_interlace,
                              cs_datatype_t             src_datatype,
                              cs_datatype_t             dest_datatype,
                              int                       n_parent_lists,
                              const cs_lnum_t           parent_num_shift[],
                              const cs_lnum_t           parent_num[],
                              const void         *const src_data[],
                              void               *const dest_data)
{
  /* If source or destination datatype is not floating-point,
     set all return values to zero */

  if (   (src_datatype != CS_DOUBLE && src_datatype != CS_FLOAT)
      || (dest_datatype != CS_DOUBLE && dest_datatype != CS_FLOAT)) {

    unsigned char *_dest_data = dest_data;

    size_t data_shift =     start_id
                          * (dest_dim * cs_datatype_size[dest_datatype]);
    size_t data_size_c =    (end_id - start_id)
                          * (dest_dim * cs_datatype_size[dest_datatype]);

    memset(_dest_data + data_shift, 0, data_size_c);

  }

  /* Otherwise, interpolate values */

  else {

    _vertex_field_of_real_values(this_tesselation,
                                 src_dim,
                                 src_dim_shift,
                                 dest_dim,
                                 start_id,
                                 end_id,
                                 src_interlace,
                                 src_datatype,
                                 dest_datatype,
                                 n_parent_lists,
                                 parent_num_shift,
                                 parent_num,
                                 src_data,
                                 dest_data);
  }

}

/*----------------------------------------------------------------------------
 * Dump printout of a mesh section tesselation structure.
 *
 * parameters:
 *   this_tesselation  <-- pointer to structure that should be dumped
 *----------------------------------------------------------------------------*/

void
fvm_tesselation_dump(const fvm_tesselation_t  *this_tesselation)
{
  int  i;
  cs_lnum_t   n_elements, j, k;
  const cs_lnum_t   *idx;

  if (this_tesselation == NULL)
    return;

  /* Global indicators */
  /*--------------------*/

  bft_printf("\n"
             "Tesselation:\n\n"
             "Element type:         %s\n"
             "Number of elements:   %ld\n"
             "Spatial dimension:    %d\n"
             "Entity dimension:     %d\n",
             fvm_elements_type_name[this_tesselation->type],
             (long)this_tesselation->n_elements,
             this_tesselation->dim, this_tesselation->entity_dim);

  bft_printf("\n"
             "Stride:                %d\n"
             "Number of faces:       %ld\n",
             this_tesselation->stride,
             (long)(this_tesselation->n_faces));

  bft_printf("\n"
             "Pointers to shared arrays:\n"
             "  vertex_coords         %p\n"
             "  parent_vertex_num     %p\n"
             "  face_index:           %p\n"
             "  face_num:             %p\n"
             "  vertex_index:         %p\n"
             "  vertex_num:           %p\n",
             (const void *)this_tesselation->vertex_coords,
             (const void *)this_tesselation->parent_vertex_num,
             (const void *)this_tesselation->face_index,
             (const void *)this_tesselation->face_num,
             (const void *)this_tesselation->vertex_index,
             (const void *) this_tesselation->vertex_num);

  bft_printf("\n"
             "Pointers to shared global numbering:\n"
             "  global_element_num    %p\n",
             (const void *)this_tesselation->global_element_num);


  /* Basic information */
  /*-------------------*/

  bft_printf("\n"
             "Number of sub-entity types:     %d\n\n",
             this_tesselation->n_sub_types);

  for (i = 0; i < this_tesselation->n_sub_types; i++) {
    bft_printf("Maximum local number of resulting %s per element: %ld\n",
               fvm_elements_type_name[this_tesselation->sub_type[i]],
               (long)(this_tesselation->n_sub_max[i]));
  }

  for (i = 0; i < this_tesselation->n_sub_types; i++) {
    bft_printf("Maximum global number of resulting %s per element: %ld\n",
               fvm_elements_type_name[this_tesselation->sub_type[i]],
               (long)(this_tesselation->n_sub_max_glob[i]));
  }

  bft_printf("\n");

  for (i = 0; i < this_tesselation->n_sub_types; i++) {
    bft_printf("Local number of resulting %s: %ld\n",
               fvm_elements_type_name[this_tesselation->sub_type[i]],
               (long)(this_tesselation->n_sub[i]));
  }

  for (i = 0; i < this_tesselation->n_sub_types; i++) {
    bft_printf("Global number of resulting %s: %llu\n",
               fvm_elements_type_name[this_tesselation->sub_type[i]],
               (unsigned long long)(this_tesselation->n_sub_glob[i]));
  }

  bft_printf("\n"
             "Pointers to shareable arrays:\n"
             "  encoding:  %p\n",
             (const void *)this_tesselation->encoding);

  for (i = 0; i < this_tesselation->n_sub_types; i++) {
    if (this_tesselation->sub_elt_index[i] != NULL)
      bft_printf("  sub_elt_index[%s]: %p\n",
                 fvm_elements_type_name[this_tesselation->sub_type[i]],
                 (const void *)this_tesselation->sub_elt_index[i]);
  }

  bft_printf("\n"
             "Pointers to local arrays:\n"
             "  _encoding: %p\n",
             (const void *)this_tesselation->_encoding);

  for (i = 0; i < this_tesselation->n_sub_types; i++) {
    if (this_tesselation->sub_elt_index[i] != NULL)
      bft_printf("  _sub_elt_index[%s]: %p\n",
                 fvm_elements_type_name[this_tesselation->sub_type[i]],
                 (const void *)this_tesselation->_sub_elt_index[i]);
  }

  if (this_tesselation->encoding != NULL) {

    fvm_tesselation_encoding_t decoding_mask[3] = {0, 0, 0};
    cs_lnum_t tv[3];

    _init_decoding_mask(decoding_mask);

    if (this_tesselation->type != FVM_FACE_QUAD) {
      bft_printf("\nEncoding (local vertex numbers):\n\n");
      if (this_tesselation->n_faces > 0)
        n_elements = this_tesselation->n_faces;
      else
        n_elements = this_tesselation->n_elements;
      idx = this_tesselation->vertex_index;
      for (j = 0; j < n_elements; j++) {
        _DECODE_TRIANGLE_VERTICES(tv,
                                  this_tesselation->encoding[idx[j] - 2*j],
                                  decoding_mask);
        bft_printf("%10d (idx = %10d) %10d %10d %10d\n",
                   j+1, idx[j], (int)tv[0], (int)tv[1], (int)tv[2]);
        for (k = idx[j] -2*j + 1; k < idx[j+1] - 2*j; k++) {
          _DECODE_TRIANGLE_VERTICES(tv,
                                    this_tesselation->encoding[k],
                                    decoding_mask);
          bft_printf("                              %10d %10d %10d\n",
                     (int)tv[0], (int)tv[1], (int)tv[2]);
        }
      }
      bft_printf("      end  (idx = %10d)\n", idx[n_elements]);
    }
    else { /* if (this_tesselation->type != FVM_FACE_QUAD) */
      bft_printf("\nEncoding (diagonal flag):\n\n");
      n_elements = this_tesselation->n_elements;
      for (j = 0; j < n_elements; j++)
        bft_printf("%10d: %10d\n", j+1, (int)this_tesselation->encoding[j]);
    }

  }

  for (i = 0; i < this_tesselation->n_sub_types; i++) {
    if (this_tesselation->sub_elt_index[i] != NULL) {
      bft_printf("\nSub-element index [%s]:\n\n",
                 fvm_elements_type_name[this_tesselation->sub_type[i]]);
      n_elements = this_tesselation->n_elements;
      idx = this_tesselation->sub_elt_index[i];
      for (j = 0; j < n_elements; j++)
        bft_printf("%10d: idx = %10d\n", j+1, idx[j]);
      bft_printf("      end: idx = %10d\n", idx[n_elements]);
    }
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
