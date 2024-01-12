/*============================================================================
 * Locate local points in a nodal representation associated with a mesh
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

/*----------------------------------------------------------------------------*/

#include "ple_config.h"

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
#include "my_ple_mesh.h"
#include "my_ple_mesh_priv.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "my_ple_point_location.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/* In case <math.h> does not define HUGE_VAL, use a "safe" value */
#if !defined(HUGE_VAL)
#define HUGE_VAL 1.0e+30
#endif

/* Absolute, minimum, and maximum values */

#define _ABS(a)     ((a) <  0  ? -(a) : (a))  /* Absolute value of a */
#define _MIN(a,b)   ((a) > (b) ?  (b) : (a))  /* Minimum of a et b */
#define _MAX(a,b)   ((a) < (b) ?  (b) : (a))  /* Maximum of a et b */

/* Geometric operation macros*/

enum {X, Y, Z};

#define _CROSS_PRODUCT_3D(cross_v1_v2, v1, v2) ( \
 cross_v1_v2[0] = v1[1]*v2[2] - v1[2]*v2[1],   \
 cross_v1_v2[1] = v1[2]*v2[0] - v1[0]*v2[2],   \
 cross_v1_v2[2] = v1[0]*v2[1] - v1[1]*v2[0]  )

#define _DOT_PRODUCT_3D(v1, v2) ( \
 v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])

#define _DOT_PRODUCT_2D(vect1, vect2) \
  (vect1[X] * vect2[X] + vect1[Y] * vect2[Y])

#define _MODULE_3D(v) \
 sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * State of current triangulation (working structure)
 *----------------------------------------------------------------------------*/

typedef struct {

  int          *triangle_vertices; /* current triangle vertices list */
  ple_coord_t  *coords;            /* vertex coordinates */
  int          *list_previous;     /* indices of previous vertices in polygon;
                                      size:n_vertices; */
  int          *list_next;         /* indices of next vertices in polygon;
                                      size:n_vertices; */
  int          *edge_vertices;     /* edges connectivity; size: n_edges * 2 */
  int          *edge_neighbors;    /* triangles sharing a given edge
                                      (2 values per edge);
                                      - (-1,-1): edge does not exist
                                      - (x ,-1): boundary edge, triangle x
                                      - (x ,-1): internal edge, triangles x
                                                 and y */
  _Bool        *edge_is_delaunay;  /* Delaunay edge indicator */
  _Bool        *concave;           /* Concave vertex indicator */

  int           n_vertices_max;    /* Maximum number vertices; a larger
                                      size requires resizing work arrays */

} _triangulate_state_t;

/*----------------------------------------------------------------------------
 * Structure defining a local octree (3d)
 *----------------------------------------------------------------------------*/

typedef struct {

  ple_lnum_t  octant_id[8];   /* Ids of sub-octants in octree array */
  ple_lnum_t  idx[9];         /* Start index of point list for each octant */
  ple_lnum_t  n_points;       /* Number of points in octree */

} _octant_t;

typedef struct {

  size_t       n_points;      /* Number of points in octree */
  size_t       n_nodes;       /* Current number of nodes in octree */
  size_t       n_nodes_max;   /* Maximum number of nodes in octree */

  double       extents[6];    /* Associated extents */

  ple_lnum_t  *point_ids;     /* Id's of points sorted by octree
                                 (size: n_points + 1) */
  _octant_t   *nodes;         /* Array of octree nodes
                                 (size: n_nodes_max) */

} _octree_t;

/*----------------------------------------------------------------------------
 * Structure defining a local quadtree (2d)
 *----------------------------------------------------------------------------*/

typedef struct {

  ple_lnum_t  quadrant_id[4]; /* Id of sub-quadrants in quadtree array */
  ple_lnum_t  idx[5];         /* Start index of point list for each quadrant */
  ple_lnum_t  n_points;       /* Number of points in quadtree */

} _quadrant_t;

typedef struct {

  size_t        n_points;     /* Number of points in quadtree */
  size_t        n_nodes;      /* Current number of nodes in quadtree */
  size_t        n_nodes_max;  /* Maximum number of nodes in quadtree */

  double        extents[4];   /* Associated extents */

  ple_lnum_t   *point_ids;    /* Id's of points sorted by quadtree
                                 (size: n_points + 1) */
  _quadrant_t  *nodes;        /* Array of quadtree nodes
                                 (size: n_nodes_max) */

} _quadtree_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

static double      _epsilon_denom = 1.e-14; /* Minimum denominator */

static int         _octree_max_levels = 18; /* Maximum number of levels */
static ple_lnum_t  _octree_threshold = 4;   /* Number of points in octree node
                                               under which the node is final */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Test if two extents intersect
 *
 * parameters:
 *   dim             <-- spatial (coordinates) dimension
 *   extents_1       <-- extents: x_min, y_min, ..., x_max, y_max, ...
 *                       size: dim*2
 *   extents_2       <-- extents: x_min, y_min, ..., x_max, y_max, ...
 *                       size: dim*2
 *
 * returns:
 *   true if extents intersect, false otherwise
 *----------------------------------------------------------------------------*/

inline static _Bool
_intersect_extents(int           dim,
                   const double  extents_1[],
                   const double  extents_2[])
{
  int i;
  _Bool retval = true;

  for (i = 0; i < dim; i++) {
    if (   (extents_1[i] > extents_2[i + dim])
        || (extents_2[i] > extents_1[i + dim])) {
      retval = false;
      break;
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Test if a point is within given extents
 *
 * parameters:
 *   dim             <-- spatial (coordinates) dimension
 *   coords          <-- coordinates: x, y, ...
 *                       size: dim
 *   extents         <-- extents: x_min, y_min, ..., x_max, y_max, ...
 *                       size: dim*2
 *
 * returns:
 *   true if point lies within extents, false otherwise
 *----------------------------------------------------------------------------*/

inline static _Bool
_within_extents(int                dim,
                const ple_coord_t  coords[],
                const double       extents[])
{
  int i;
  _Bool retval = true;

  for (i = 0; i < dim; i++) {
    if (   (coords[i] < extents[i])
        || (coords[i] > extents[i + dim])) {
      retval = false;
      break;
    }
  }

  return retval;
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
 * Update element extents with a given vertex
 *
 * parameters:
 *   dim               <-- spatial (coordinates) dimension
 *   vertex_id         <-- vertex index (0 to n-1)
 *   vertex_coords     <-- pointer to vertex coordinates
 *   elt_extents       <-> extents associated with element:
 *                         x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *   elt_initialized   <-> are extents already initialized for this vertex
 *                         (for all element vertices except the first) ?
 *----------------------------------------------------------------------------*/

inline static void
_update_elt_extents(int                 dim,
                    ple_lnum_t          vertex_id,
                    const ple_coord_t   vertex_coords[],
                    double              elt_extents[],
                    _Bool              *elt_initialized)
{
  ple_lnum_t  i, coord_idx;

  coord_idx = vertex_id;

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
 * Compute extents of a mesh representation section
 *
 * parameters:
 *   this_section      <-- pointer to section structure
 *   dim               <-- spatial (coordinates) dimension
 *   vertex_coords     <-- pointer to vertex coordinates
 *   tolerance         <-- addition to local extents of each element:
 *                         extent = base_extent * (1 + tolerance)
 *   extents           <-> extents associated with section:
 *                         x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *----------------------------------------------------------------------------*/

static void
_mesh_section_extents(const my_ple_mesh_section_t  *this_section,
                      int                           dim,
                      const ple_coord_t             vertex_coords[],
                      double                        tolerance,
                      double                        extents[])
{
  ple_lnum_t  i, j, k, face_id, vertex_id;
  double elt_extents[6];

  /* initialize extents in case section is empty */
  for (i = 0; i < dim; i++) {
    extents[i]       =  HUGE_VAL;
    extents[i + dim] = -HUGE_VAL;
  }

  /* Extents for polyhedra */

  if (this_section->face_index != NULL) {

    for (i = 0; i < this_section->n_elements; i++) {

      _Bool elt_initialized = false;

      for (j = this_section->face_index[i];
           j < this_section->face_index[i + 1];
           j++) {
        face_id = _ABS(this_section->face_num[j]) - 1;
        for (k = this_section->vertex_index[face_id];
             k < this_section->vertex_index[face_id + 1];
             k++) {
          vertex_id = this_section->vertex_num[k] - 1;

          _update_elt_extents(dim,
                              vertex_id,
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

    ple_lnum_t  n_faces = (this_section->n_faces > 0) ?
                          this_section->n_faces : this_section->n_elements;

    for (i = 0; i < n_faces; i++) {

      _Bool elt_initialized = false;

      for (j = this_section->vertex_index[i];
           j < this_section->vertex_index[i + 1];
           j++) {
        vertex_id = this_section->vertex_num[j] - 1;

        _update_elt_extents(dim,
                            vertex_id,
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

      _Bool elt_initialized = false;

      for (j = 0; j < this_section->stride; j++) {

        vertex_id = this_section->vertex_num[i*this_section->stride + j] - 1;

        _update_elt_extents(dim,
                            vertex_id,
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

/*----------------------------------------------------------------------------
 * Compute extents of a point set
 *
 * parameters:
 *   dim          <-- space dimension of points to locate
 *   n_points     <-- number of points to locate
 *   point_index  <-- optional indirection array to point_coords
 *                    (1 to n_points numbering)
 *   point_coords <-- coordinates of points to locate
 *                    (dimension: dim * n_points)
 *   extents      --> extents associated with mesh:
 *                    x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *----------------------------------------------------------------------------*/

static void
_point_extents(const int            dim,
               const ple_lnum_t     n_points,
               const ple_lnum_t     point_index[],
               const ple_coord_t    point_coords[],
               double               extents[])
{
  int i;
  ple_lnum_t j, coord_idx;

  /* initialize extents in case mesh is empty or dim < 3 */
  for (i = 0; i < dim; i++) {
    extents[i]       =  HUGE_VAL;
    extents[i + dim] = -HUGE_VAL;
  }

  /* Compute extents */

  if (point_index != NULL) {

    for (j = 0; j < n_points; j++) {
      coord_idx = point_index[j] - 1;
      for (i = 0; i < dim; i++) {
        if (extents[i]       > point_coords[(coord_idx * dim) + i])
          extents[i]       = point_coords[(coord_idx * dim) + i];
        if (extents[i + dim] < point_coords[(coord_idx * dim) + i])
          extents[i + dim] = point_coords[(coord_idx * dim) + i];
      }
    }
  }

  else {

    for (coord_idx = 0; coord_idx < n_points; coord_idx++) {
      for (i = 0; i < dim; i++) {
        if (extents[i]       > point_coords[(coord_idx * dim) + i])
          extents[i]       = point_coords[(coord_idx * dim) + i];
        if (extents[i + dim] < point_coords[(coord_idx * dim) + i])
          extents[i + dim] = point_coords[(coord_idx * dim) + i];
      }
    }
  }

}

/*----------------------------------------------------------------------------
 * Project a face in 3d on a plane parallel to this face, which is then
 * treated as the (Oxy) plane.
 *
 * parameters:
 *   n_vertices  <-- number of vertices defining the polygon.
 *   coords      <-> coordinates of the polygon's vertices (3d in, 2d out)
 *----------------------------------------------------------------------------*/

static void
_polygon_plane_3d(const int    n_vertices,
                  ple_coord_t  coords[])
{

  int i, j;

  ple_coord_t  face_center[3], face_normal[3];
  ple_coord_t  v1[3], v2[3], cross_v1_v2[3];
  ple_coord_t  dot_v1_v2, module_v2;
  ple_coord_t  tmp_coord;

  double   cost;
  double   sint;

#define _N_VERTICES_AUTO_MAX   20 /* Size of local temporary coordinates
                                     buffer; above this size, allocation
                                     is necessary */

  ple_coord_t  _tmp_coords[_N_VERTICES_AUTO_MAX * 3];
  ple_coord_t  *_tmp_coords_p = NULL;
  ple_coord_t  *tmp_coords = _tmp_coords;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Estimate position of polygon center */

  for (i = 0; i < 3; i++)
    face_center[i] = 0.;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < n_vertices; j++)
      face_center[i] += coords[j*3 + i];
    face_center[i] /= n_vertices;
  }

  /* Estimate face normal */

  for (i = 0; i < 3; i++)
    face_normal[i] = 0.;

  for (j = 0; j < n_vertices; j++) {

    for (i = 0; i < 3; i++) {
      v1[i] = coords[j*3 + i] - face_center[i];
      if (j < n_vertices - 1)
        v2[i] = coords[(j+1)*3 + i] - face_center[i];
      else
        v2[i] = coords[          i] - face_center[i];
    }

    face_normal[0] += v1[1]*v2[2] - v1[2]*v2[1];
    face_normal[1] += v1[2]*v2[0] - v1[0]*v2[2];
    face_normal[2] += v1[0]*v2[1] - v1[1]*v2[0];

  }

  /* Project coordinates in a plane parallel to face */
  /*-------------------------------------------------*/

  /* We place the coordinate system origin at the estimated face center */

  for (j = 0; j < n_vertices; j++)
    for (i = 0; i < 3; i++)
      coords[j*3 + i] -= face_center[i];

  if (_ABS(face_normal[0]) > 1.e-12 || _ABS(face_normal[1]) > 1.e-12) {

    /* First rotation of axis (Oz) and angle (Ox, normal proj. on Oxy) */

    if (n_vertices > _N_VERTICES_AUTO_MAX) {
      PLE_MALLOC(_tmp_coords_p, n_vertices*3, ple_coord_t);
      tmp_coords = _tmp_coords_p;
    }

    v1[0] = 1.;
    v1[1] = 0.;
    v1[2] = 0.;

    v2[0] = face_normal[0];
    v2[1] = face_normal[1];
    v2[2] = 0.;

    _CROSS_PRODUCT_3D(cross_v1_v2, v1, v2);

    dot_v1_v2 = _DOT_PRODUCT_3D(v1, v2);
    module_v2 = _MODULE_3D(v2);
    cost = dot_v1_v2 / module_v2;

    if (cross_v1_v2[2] > 0.)
      sint =   _MODULE_3D(cross_v1_v2) / module_v2;
    else
      sint = - _MODULE_3D(cross_v1_v2) / module_v2;

    for (j = 0; j < n_vertices; j++) {

      tmp_coords[j*3    ] =  cost*coords[j*3] + sint*coords[j*3 + 1];
      tmp_coords[j*3 + 1] = -sint*coords[j*3] + cost*coords[j*3 + 1];
      tmp_coords[j*3 + 2] =  coords[j*3 +2];

    }

    /* Second rotation with axis (Oy) and angle (Oz', normal proj. on  Ox'z) */

    v1[0] =  0.;
    v1[1] =  0.;
    v1[2] =  1.;

    v2[0] = sqrt(face_normal[0]*face_normal[0] + face_normal[1]*face_normal[1]);
    v2[1] = 0.;
    v2[2] = face_normal[2];

    _CROSS_PRODUCT_3D(cross_v1_v2, v1, v2);

    cost = _DOT_PRODUCT_3D(v1, v2) / _MODULE_3D(v2);

    if (cross_v1_v2[2] > 0.)
      sint =  _MODULE_3D(cross_v1_v2) / _MODULE_3D(v2);
    else
      sint = -_MODULE_3D(cross_v1_v2) / _MODULE_3D(v2);

    for (j = 0; j < n_vertices; j++) {

      coords[j*3    ] = cost*tmp_coords[j*3] + sint*tmp_coords[j*3 + 2];
      coords[j*3 + 1] = tmp_coords[j*3 + 1];
      coords[j*3 + 2] = 0.;

    }

    if (_tmp_coords_p != NULL) {
      PLE_FREE(_tmp_coords_p);
      tmp_coords = NULL;
    }

  }
  else {

    /* We only need to set the vertices z coordinate to 0, possibly
       swapping the coordinates in the (Oxy) projection plane.  */

    if (face_normal[2] > 0.)
      for (j = 0; j < n_vertices; j++)
        coords[j*3 + 2] = 0.;

    else
      for (j = 0; j < n_vertices; j++) {
        tmp_coord = coords[j*3];
        coords[j*3    ] = coords[j*3 + 1];
        coords[j*3 + 1] = tmp_coord;
        coords[j*3 + 2] = 0.;
      }

  }

  /* Convert coords to 2d */

  for (j = 0; j < n_vertices; j++) {
    coords[j*2    ] = coords[j*3    ];
    coords[j*2 + 1] = coords[j*3 + 1];
  }

#undef _N_VERTICES_AUTO_MAX
}

/*----------------------------------------------------------------------------
 * Check if a (2D projection of a) polygon vertex is convex.
 *
 * parameters:
 *   previous    <-- index of previous vertex in polygon.
 *   current     <-- index of current vertex in polygon.
 *   next        <-- index of following vertex in polygon.
 *   coords      <-- coordinates of the polygon's vertices (2d).
 *----------------------------------------------------------------------------*/

static _Bool
_polygon_vertex_is_convex(const int    previous,
                          const int    current,
                          const int    next,
                          const ple_coord_t  coords[])
{
  /* sin(theta) = (v1 x v2) / (|v1| x |v2|), so the sign of the sine
     is that of the cross product's z component (as v1 and v2 lie in
     the same plane) */

  if (  (  (coords[current*2]     - coords[previous*2])
         * (coords[next*2 + 1]    - coords[current*2 + 1]))
      - (  (coords[next*2]        - coords[current*2])
         * (coords[current*2 + 1] - coords[previous*2 + 1])) > 0.0)
    return true;

  else
    return false;
}


/*----------------------------------------------------------------------------
 * Check if a (2D projection of a) polygon vertex is an ear.
 *
 * parameters:
 *   n_vertices    <-- number of vertices defining the polygon.
 *   current       <-- index of current vertex in polygon.
 *   list_previous <-- indices of previous vertices in polygon.
 *   list_next     <-- index of previous vertex in polygon.
 *   concave       <-- flag concave vertices in polygon.
 *   coords        <-- coordinates of the polygon's vertices (2d).
 *   epsilon       <-- associated relative tolerance.
 *----------------------------------------------------------------------------*/

static _Bool
_polygon_vertex_is_ear(const int          n_vertices,
                       const int          current,
                       const int          list_previous[],
                       const int          list_next[],
                       const _Bool        concave[],
                       const ple_coord_t  coords[],
                       const double       epsilon)

{
  int i, previous, next;
  double surf_2, x_iso, y_iso;
  double vect1[2], vect2[2], vect3[2];

  /* If no vertex is concave, we have an ear */

  for (i = 0; i < n_vertices && concave[i] == false; i++);

  if (i == n_vertices)
    return true;

  /* If current vertex is convex */

  else {

    if (concave[current] == false) {

      /* We check if the triangle formed by
         list_previous[current], current, and list_next[current]
         contains a concave vertex */

      previous = list_previous[current];
      next     = list_next[current];

      vect2[0] = coords[current*2    ] - coords[previous*2    ];
      vect2[1] = coords[current*2 + 1] - coords[previous*2 + 1];
      vect3[0] = coords[next*2    ]    - coords[previous*2    ];
      vect3[1] = coords[next*2 + 1]    - coords[previous*2 + 1];

      surf_2 = vect2[0]*vect3[1] - vect3[0]*vect2[1];

      for (i = list_next[next]; i != previous; i = list_next[i]) {

        if (concave[i] == true) {

          vect1[0] = coords[i*2    ] - coords[previous*2    ];
          vect1[1] = coords[i*2 + 1] - coords[previous*2 + 1];

          x_iso = (vect1[0]*vect3[1] - vect1[1]*vect3[0]) / surf_2;
          y_iso = (vect2[0]*vect1[1] - vect2[1]*vect1[0]) / surf_2;

          if (   (1.0 - x_iso - y_iso > - epsilon)
              && (      x_iso         > - epsilon)
              && (              y_iso > - epsilon))
            return false;

        }

      }

      return true;

    }
    else
      return false;

  }

}

/*----------------------------------------------------------------------------
 * Check if an edge (between two triangles) is locally Delaunay.
 *
 * We compute the power of a point compared to the circle circumscribed
 * to the triangle formed by the three others (of which two define
 * the diagonal considered).
 *
 * parameters:
 *   edge_vertex_0     <-- index of first edge vertex in coords.
 *   edge_vertex_1     <-- index of second edge vertex in coords.
 *   flip_vertex_0     <-- index of first flip edge vertex in coords.
 *   flip_vertex_1     <-- index of second flip edge vertex in coords.
 *   coords            <-- coordinates of the triangulation's vertices (2d).
 *----------------------------------------------------------------------------*/

static _Bool
_edge_is_locally_delaunay(const int          edge_vertex_0,
                          const int          edge_vertex_1,
                          const int          flip_vertex_0,
                          const int          flip_vertex_1,
                          const ple_coord_t  coords[])
{
  double   a, b, delta;
  double   lambda[4];
  double   x_center, y_center, radius;
  double   point_power;

  lambda[0] = 2*(coords[edge_vertex_1*2    ] - coords[edge_vertex_0*2    ]);
  lambda[1] = 2*(coords[edge_vertex_1*2 + 1] - coords[edge_vertex_0*2 + 1]);

  lambda[2] = 2*(coords[flip_vertex_0*2    ] - coords[edge_vertex_0*2    ]);
  lambda[3] = 2*(coords[flip_vertex_0*2 + 1] - coords[edge_vertex_0*2 + 1]);

  delta = lambda[1]*lambda[2] - lambda[0]*lambda[3];

  /* If the triangle is flat, we automatically switch diagonals to
     avoid a division by zero. */

  if (_ABS(delta) < 1.e-12)
    return false;

  a =   coords[edge_vertex_1*2    ] * coords[edge_vertex_1*2    ]
      - coords[edge_vertex_0*2    ] * coords[edge_vertex_0*2    ]
      + coords[edge_vertex_1*2 + 1] * coords[edge_vertex_1*2 + 1]
      - coords[edge_vertex_0*2 + 1] * coords[edge_vertex_0*2 + 1];

  b =   coords[flip_vertex_0*2    ] * coords[flip_vertex_0*2    ]
      - coords[edge_vertex_0*2    ] * coords[edge_vertex_0*2    ]
      + coords[flip_vertex_0*2 + 1] * coords[flip_vertex_0*2 + 1]
      - coords[edge_vertex_0*2 + 1] * coords[edge_vertex_0*2 + 1];

  /* Center and radius of the circle passing through the vertices of
     the diagonal and through a third vertex. */

  x_center = (lambda[1]*b - lambda[3]*a)/delta;
  y_center = (lambda[2]*a - lambda[0]*b)/delta;

  radius = sqrt( (  (x_center - coords[edge_vertex_0*2    ])
                  * (x_center - coords[edge_vertex_0*2    ]))
               + (  (y_center - coords[edge_vertex_0*2 + 1])
                  * (y_center - coords[edge_vertex_0*2 + 1])));

  /* Compute the power of the quadrilateral's fourth vertex compared
     to the circle. */

  point_power =   (  (coords[flip_vertex_1*2    ] - x_center)
                   * (coords[flip_vertex_1*2    ] - x_center))
                + (  (coords[flip_vertex_1*2 + 1] - y_center)
                   * (coords[flip_vertex_1*2 + 1] - y_center))
                - radius*radius;

  /* Keep a slight margin in case there is no perceptible gain
     in switching diagonals */

  if (point_power > -1.e-12)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------
 * Sort vertex indexes defining a triangle.
 *
 * parameters:
 *   triangle_vertices <-> triangles connectivity.
 *----------------------------------------------------------------------------*/

static void
_triangle_by_sorted_vertices(int triangle_vertices[])
{
  int  i_min, i_max, i_mid;

  /* First step */

  if (triangle_vertices[0] < triangle_vertices[1]) {
    i_min = triangle_vertices[0];
    i_mid = triangle_vertices[1];
    i_max = i_mid;
  }
  else {
    i_min = triangle_vertices[1];
    i_mid = triangle_vertices[0];
    i_max = i_mid;
  }

  /* Second step */

  if (triangle_vertices[2] < i_min) {
    i_mid = i_min;
    i_min = triangle_vertices[2];
  }
  else if (triangle_vertices[2] > i_max) {
    i_max = triangle_vertices[2];
  }
  else {
    i_mid = triangle_vertices[2];
  }

  /* Reordering */

  triangle_vertices[0] = i_min;
  triangle_vertices[1] = i_mid;
  triangle_vertices[2] = i_max;
}

/*----------------------------------------------------------------------------
 * Convert a given triangulation to a Delaunay triangulation using
 * edge flips.
 *
 * parameters:
 *   n_vertices        <-- number of vertices defining the triangulation.
 *   triangle_vertices <-> triangles connectivity.
 *   edge_vertices     --> edges connectivity.
 *   edge_neigbors     --> triangles sharing a given edge (2 values per edge);
 *                          - (-1,-1): edge does not exist
 *                          - (x ,-1): boundary edge, triangle x
 *                          - (x ,-1): internal edge, triangles x and y
 *   edge_is_delaunay  --> delaunay edge indicator.
 *   coords            <-- coordinates of the triangulation's vertices (2d).
 *----------------------------------------------------------------------------*/

static void
_polygon_delaunay_flip(const int          n_vertices,
                       int                triangle_vertices[],
                       int                edge_vertices[],
                       int                edge_neighbors[],
                       _Bool              edge_is_delaunay[],
                       const ple_coord_t  coords[])
{
  int    triangle_id, edge_id, vertex_id;
  int    triangle_id_0, triangle_id_1;

  int    current_edge_id, flip_edge_id;

  int    vertex_flip[2];

  int    i, i_0, i_1, i_min, i_max, j;

  _Bool   face_is_delaunay, edge_locally_delaunay;
  _Bool   restart, convex_quad;

  const int  n_edges = n_vertices*(n_vertices - 1)/2;
  const int  n_triangles = n_vertices - 2;
  int        i_previous[2] = {-1, -1}, i_next[2] = {-1, -1};

  /*
    There are (n_vertices*(n_vertices - 1)/2) possible edge combinations.
    An edge's number can be given by:
    -> SUM(k = 0, k < i_min) {n_vertices - k} + (i_max - i_min)
    i.e. n_vertices*i_min - i_min*(i_min+1)/2 + (i_max - i_min)
    An edge's index equals it's number - 1.

    Be careful with the order of arguments! arg[0]=i_min, arg[1]=i_max
  */

#undef _EDGE_INDEX
#define _EDGE_INDEX(i_min, i_max) \
(n_vertices*i_min - i_min*(i_min+1)/2 + i_max-i_min - 1)

  /* Initialization */
  /*----------------*/

  for (i_0 = 0; i_0 < n_vertices; i_0++) {
    for (i_1 = i_0 + 1; i_1 < n_vertices; i_1++) {

      edge_id = _EDGE_INDEX(i_0, i_1);

      /* Define edges */

      edge_vertices[2*edge_id    ] = i_0;
      edge_vertices[2*edge_id + 1] = i_1;

      /*
        Liste of triangles sharing an edge:
        - (-1,-1): edge does not exist
        - (x ,-1): boundary edge, triangle x
        - (x ,-1): internal edge, triangles x and y
      */

      edge_neighbors[2*edge_id    ] = -1;
      edge_neighbors[2*edge_id + 1] = -1;

      /* Initialize an array indicating if an edge is locally Delaunay */

      edge_is_delaunay[edge_id] = true;

    }
  }

  /* First traversal of triangles to determine initial neighborhood,
     as well as the list of edges which are not locally Delaunay */

  for (j = 0; j < n_triangles; j++) {
    for (i = 0; i < 3; i++) {

      i_0 = triangle_vertices[(j*3) +   i      ];
      i_1 = triangle_vertices[(j*3) + ((i+1)%3)];

      i_min = _MIN(i_0, i_1);
      i_max = _MAX(i_0, i_1);

      edge_id = _EDGE_INDEX(i_min, i_max);

      /* Update edge neighbors */

      if (edge_neighbors[2*edge_id] == -1)
        edge_neighbors[2*edge_id] = j;
      else
        edge_neighbors[2*edge_id + 1] = j;

      /* If edge is not on the boundary: */

      if (   !(i_max == i_min + 1)
          && !(i_min == 0 && i_max == n_vertices - 1))
        edge_is_delaunay[edge_id] = false;

    }
  }

  /* Main flip algorithm */
  /*---------------------*/

  edge_id = 0;

  restart = false;
  face_is_delaunay = false;

  while(face_is_delaunay == false) {

    if (edge_is_delaunay[edge_id] == false) {

      edge_is_delaunay[edge_id] = true;

      i_0 = edge_vertices[2*edge_id];
      i_1 = edge_vertices[2*edge_id + 1];

      for (j = 0; j < 2; j++) { /* Loop on triangles on each side of edge */

        triangle_id = edge_neighbors[2*edge_id + j];

        for (i = 0; i < 3; i++) { /* Loop on triangle's vertices */

          vertex_id = triangle_vertices[3*triangle_id + i];

          /* Seek opposite vertices */

          if (vertex_id != i_0 && vertex_id != i_1)
            vertex_flip[j] = vertex_id;

          /* Seek preceding and following vertices */

          if (   vertex_id == i_0
              && triangle_vertices[3*triangle_id + ((i + 1)%3)] == i_1) {
            i_previous[0] = triangle_vertices[3*triangle_id + ((i + 2)%3)];
            i_next[1] = i_previous[0];
          }
          else if (   vertex_id == i_1
                   && triangle_vertices[3*triangle_id + ((i + 1)%3)] == i_0) {
            i_next[0] = triangle_vertices[3*triangle_id + ((i + 2)%3)];
            i_previous[1] = i_next[0];
          }

        } /* End of loop on triangle's vertices */

      } /* End of loop on triangles on each side of edge */

      /* Test quadrilateral's convexity */

      if (   _polygon_vertex_is_convex(i_previous[0],
                                       i_0,
                                       i_next[0],
                                       coords) == true
          && _polygon_vertex_is_convex(i_previous[1],
                                       i_1,
                                       i_next[1],
                                       coords) == true)
        convex_quad = true;

      else
        convex_quad = false;

      /* Check if edge is locally Delaunay */

      edge_locally_delaunay =
        _edge_is_locally_delaunay(i_0,
                                  i_1,
                                  vertex_flip[0],
                                  vertex_flip[1],
                                  coords);

      /* If edge is not locally Delaunay */
      /*---------------------------------*/

      if (   edge_locally_delaunay == false
          && convex_quad == true) {

        i_min = _MIN(vertex_flip[0], vertex_flip[1]);
        i_max = _MAX(vertex_flip[0], vertex_flip[1]);

        flip_edge_id = _EDGE_INDEX(i_min, i_max);

        for (j = 0; j < 2; j++) { /* Loop on triangles on each side of edge */

          triangle_id_0 = edge_neighbors[2*edge_id + j];
          triangle_id_1 = edge_neighbors[2*edge_id + ((j + 1)%2)];

          /* Redefine adjacent triangles */

          for (i = 0; i < 3; i++) {

            vertex_id = triangle_vertices[3*triangle_id_0 + i];

            if (vertex_id == edge_vertices[2*edge_id + ((j + 1)%2)])
              triangle_vertices[3*triangle_id_0 + i] = vertex_flip[(j + 1)%2];

          }

          /* Redefine triangle so that vertices appear in increasing order */

          _triangle_by_sorted_vertices(triangle_vertices + 3*triangle_id_0);

          /* Update neighborhood and non locally Delaunay edges */

          for (i = 0; i < 3; i++) { /* Loop on triangle's vertices */

            i_0 = triangle_vertices[3*triangle_id_0 + i];
            i_1 = triangle_vertices[3*triangle_id_0 + ((i + 1)%3)];

            i_min = _MIN(i_0, i_1);
            i_max = _MAX(i_0, i_1);

            current_edge_id = _EDGE_INDEX(i_min, i_max);

            if (current_edge_id < edge_id)
              restart = true;

            /* If current edge is not the new (flip) egde ... */

            if (current_edge_id != flip_edge_id) {

              if (edge_neighbors[2*current_edge_id] == triangle_id_1)
                edge_neighbors[2*current_edge_id] = triangle_id_0;
              else if (edge_neighbors[2*current_edge_id + 1] == triangle_id_1)
                edge_neighbors[2*current_edge_id + 1] = triangle_id_0;

              /* ... and is not either a boundary edge */

              if (edge_neighbors[2*current_edge_id + 1] != -1)
                edge_is_delaunay[current_edge_id] = false;

            }

          } /* End of loop on triangle's vertices */

        } /* End of loop on triangles on each side of edge */

        triangle_id_0 = edge_neighbors[2*edge_id];
        triangle_id_1 = edge_neighbors[2*edge_id + 1];

        edge_neighbors[2*flip_edge_id] = triangle_id_0;
        edge_neighbors[2*flip_edge_id + 1] = triangle_id_1;

        edge_neighbors[2*edge_id] = -1;
        edge_neighbors[2*edge_id + 1] = -1;

      } /* End if edge was not locally Delaunay */

    } /* End if edge was initially supposed non locally Delaunay */

    if (edge_id == n_edges - 1 && restart == true) {
      restart = false;
      edge_id = 0;
    }
    else if (edge_id == n_edges - 1 && restart == false)
      face_is_delaunay = true;
    else
      edge_id++;


  } /* End of flip algorithm */

}

/*----------------------------------------------------------------------------
 * Create a structure necessary to the polygon triangulation algorithm.
 *
 * parameters:
 *   n_vertices_max    <-- maximum expected number of vertices per polygon.
 *
 * returns:
 *   pointer to polygon triangulation state structure.
 *----------------------------------------------------------------------------*/

static _triangulate_state_t *
_triangulate_state_create(const int  n_vertices_max)
{
  _triangulate_state_t  *this_state = NULL;

  int n_edges_max = (2*n_vertices_max) - 3;
  int n_edges_tot_max = n_edges_max * (n_edges_max - 1) / 2;

  PLE_MALLOC(this_state, 1, _triangulate_state_t);

  if (n_vertices_max > 3) {
    PLE_MALLOC(this_state->triangle_vertices, (n_vertices_max - 2) * 3, int);
    PLE_MALLOC(this_state->coords, n_vertices_max*3, ple_coord_t);
    PLE_MALLOC(this_state->list_previous, n_vertices_max, int);
    PLE_MALLOC(this_state->list_next, n_vertices_max, int);
    PLE_MALLOC(this_state->edge_vertices, n_edges_tot_max*2, int);
    PLE_MALLOC(this_state->edge_neighbors, n_edges_tot_max*2, int);
    PLE_MALLOC(this_state->edge_is_delaunay, n_edges_tot_max, _Bool);
    PLE_MALLOC(this_state->concave, n_vertices_max, _Bool);
  }
  else {
    this_state->triangle_vertices = NULL;
    this_state->coords = NULL;
    this_state->list_previous = NULL;
    this_state->list_next = NULL;
    this_state->edge_vertices = NULL;
    this_state->edge_neighbors = NULL;
    this_state->edge_is_delaunay = NULL;
    this_state->concave = NULL;
  }

  this_state->n_vertices_max = n_vertices_max;

  return this_state;
}

/*----------------------------------------------------------------------------
 * Destroy a structure necessary to the polygon triangulation algorithm.
 *
 * parameters:
 *   this_state  <-> pointer to pointer to structure that should be destroyed.
 *----------------------------------------------------------------------------*/

static void
_triangulate_state_destroy(_triangulate_state_t  **this_state)
{
  if (this_state != NULL) {
    _triangulate_state_t  *s = *this_state;
    if (s->triangle_vertices != NULL) {
      PLE_FREE(s->triangle_vertices);
      PLE_FREE(s->coords);
      PLE_FREE(s->list_previous);
      PLE_FREE(s->list_next);
      PLE_FREE(s->edge_vertices);
      PLE_FREE(s->edge_neighbors);
      PLE_FREE(s->edge_is_delaunay);
      PLE_FREE(s->concave);
    }
    PLE_FREE(*this_state);
  }
}

/*----------------------------------------------------------------------------
 * Triangulate a polygonal face.
 *
 * For a polygon with n vertices, we should obtain a triangluation with
 * (n-2) triangles and (2n-3) edges. If the polygon_vertices argument
 * is NULL, 1, 2, ...,n local numbering is implied.
 *
 * parameters:
 *   dim               <-- spatial dimension (2 or 3).
 *   n_vertices        <-- number of vertices defining the polygon.
 *   coords            <-- coordinates of the triangulation's vertices.
 *   polygon_vertices  <-- polygon connectivity; size: n_vertices or empty.
 *   mode              <-- triangles connectivity by vertex number or
 *                         polygon vertex index (1 to n).
 *   triangle_vertices --> triangles connectivity;
 *                         size: (n_vertices - 2) * 3.
 *   state             <-> associated triangulation state structure.
 *
 * returns:
 *   number of resulting triangles.
 *----------------------------------------------------------------------------*/

static int
_triangulate_polygon(int                    dim,
                     int                    n_vertices,
                     const ple_coord_t      coords[],
                     const ple_lnum_t       polygon_vertices[],
                     ple_lnum_t             triangle_vertices[],
                     _triangulate_state_t  *state)
{
  int i, j;

  int n_triangles = 0;
  int n_tries = 0;
  double epsilon[] = {1.0e-1, 1.0e-2, 0.0, -1.0e-2, -1.0e-1};

  int  *const list_previous = state->list_previous;
  int  *const list_next = state->list_next;
  _Bool  *const concave = state->concave;

  /* Initialize state structure */

  if (n_vertices > state->n_vertices_max) {

    int n_vertices_max = n_vertices*2;
    int n_edges_max = (2*n_vertices_max) - 3;
    int n_edges_tot_max = n_edges_max * (n_edges_max - 1) / 2;

    state->n_vertices_max = n_vertices_max;
    PLE_REALLOC(state->triangle_vertices, (n_vertices_max - 2) * 3, int);
    PLE_REALLOC(state->coords, n_vertices_max*3, ple_coord_t);
    PLE_REALLOC(state->list_previous, n_vertices_max, int);
    PLE_REALLOC(state->list_next, n_vertices_max, int);
    PLE_REALLOC(state->edge_vertices, n_edges_tot_max*2, int);
    PLE_REALLOC(state->edge_neighbors, n_edges_tot_max*2, int);
    PLE_REALLOC(state->edge_is_delaunay, n_edges_tot_max, _Bool);
    PLE_REALLOC(state->concave, n_vertices_max, _Bool);
  }

  if (polygon_vertices != NULL) {
    for (i = 0; i < n_vertices; i++) {
      for (j = 0; j < dim; j++)
        state->coords[i*dim + j] = coords[(polygon_vertices[i]-1)*dim + j];
    }
  }
  else {
    for (i = 0; i < (n_vertices * dim); i++)
      state->coords[i] = coords[i];
  }

  /* Determine the work plane (3d coords are overwritten) */

  if (dim == 3)
    _polygon_plane_3d(n_vertices, state->coords);

  /* Initialization */

  while ((n_triangles != (n_vertices - 2)) &&
         n_tries < 5) {

    n_triangles = 0; /* Reset if previous try was a failure */

    for (i = 0; i < n_vertices; i++) {
      list_previous[i] = i - 1;
      list_next[i] = i + 1;
    }
    list_previous[0] = n_vertices - 1;
    list_next[n_vertices - 1] = 0;

    for (i = 0; i < n_vertices; i++) {
      if (_polygon_vertex_is_convex(list_previous[i],
                                    i,
                                    list_next[i],
                                    state->coords) == true)
        concave[i] = false;
      else
        concave[i] = true;
    }

    i = 2;

    while (i != 0 && i != n_vertices) {

      if (_polygon_vertex_is_ear(n_vertices,
                                 list_previous[i],
                                 list_previous,
                                 list_next,
                                 state->concave,
                                 state->coords,
                                 epsilon[n_tries]) == true) {

        /* Add a triangle with vertices list_previous[list_previous[i]],
           list_previous[i], i */

        state->triangle_vertices[n_triangles*3    ] = list_previous
                                                        [list_previous[i]];
        state->triangle_vertices[n_triangles*3 + 1] = list_previous[i];
        state->triangle_vertices[n_triangles*3 + 2] = i;

        n_triangles += 1;

        /* Cut the ear corresponding to list_previous[i] */

        list_previous[i] = list_previous[list_previous[i]];
        list_next[list_previous[i]] = i;

        if (   (concave[i] == true)
            && (_polygon_vertex_is_convex(list_previous[i],
                                          i,
                                          list_next[i],
                                          state->coords) == true))
          concave[i] = false;

        if (   (concave[list_previous[i]] == true)
            && (_polygon_vertex_is_convex(list_previous[list_previous[i]],
                                          list_previous[i],
                                          i,
                                          state->coords) == true))
          concave[list_previous[i]] = false;

        if (list_previous[i] == 0)
          i = list_next[i];

      }
      else /* ! _polygon_vertex_is_ear(...) */

        i = list_next[i];

    }

    n_tries++;

  }


  /* Now that we have an initial triangulation, apply flip algorithm
     to obtain a Delaunay triangulation */

  if (n_triangles == n_vertices - 2)
    _polygon_delaunay_flip(n_vertices,
                           state->triangle_vertices,
                           state->edge_vertices,
                           state->edge_neighbors,
                           state->edge_is_delaunay,
                           state->coords);

  /* Update triangle_vertices argument */

  if (polygon_vertices != NULL) {
    for (i = 0; i < n_triangles * 3; i++)
      triangle_vertices[i] = polygon_vertices[state->triangle_vertices[i]];
  }
  else {
    for (i = 0; i < n_triangles * 3; i++)
      triangle_vertices[i] = state->triangle_vertices[i] + 1;
  }

  return n_triangles;
}

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
 *   coords               <-- coordinates of the triangulation's vertices.
 *   quadrangle_vertices  <-- polygon connectivity; size: n_vertices or empty.
 *   triangle_vertices    --> triangles connectivity; size: 2 * 3.
 *
 * returns:
 *   number of resulting triangles.
 *----------------------------------------------------------------------------*/

static int
_triangulate_quadrangle(int                dim,
                        const ple_coord_t  coords[],
                        const ple_lnum_t   quadrangle_vertices[],
                        ple_lnum_t         triangle_vertices[])
{
  int i, j;
  double d2_02, d2_13;
  int o_count = 0, o_id = 0;
  ple_lnum_t  vertex_id[4] = {0, 1, 2, 3};
  double v1[3] = {0.0, 0.0, 0.0}, v2[3] = {0.0, 0.0, 0.0};
  double n0[3] = {0.0, 0.0, 0.0}, ni[3] = {0.0, 0.0, 0.0};

  if (quadrangle_vertices != NULL) {
    for (i = 0; i < 4 ; i++)
      vertex_id[i] = quadrangle_vertices[i] - 1;
  }

  /* Check for an obtuse angle */

  for (i = 0; i < dim; i++) {
    v1[i] = coords[vertex_id[1]*dim + i] - coords[vertex_id[0]*dim + i];
    v2[i] = coords[vertex_id[3]*dim + i] - coords[vertex_id[0]*dim + i];
  }

  _CROSS_PRODUCT_3D(n0, v1, v2);

  for (j = 1; j < 4; j++) {
    for (i = 0; i < dim; i++) {
      v1[i] = coords[vertex_id[(j+1)%4]*dim + i] - coords[vertex_id[j]*dim + i];
      v2[i] = coords[vertex_id[ j-1   ]*dim + i] - coords[vertex_id[0]*dim + i];
    }
    _CROSS_PRODUCT_3D(ni, v1, v2);
    if (_DOT_PRODUCT_3D(n0, ni) < 0) {
      o_count++;
      o_id = j;
    }
  }

  /* With an obtuse angle, only one diagonal lies inside the quadrangle;
     we define it as "shorter" */

  if (o_count > 0) {

    if (o_count > 1)
      o_id = 0;

    if (o_id%2 == 0) {
      d2_02 = 0.;
      d2_13 = 1.;
    }
    else {
      d2_02 = 1.;
      d2_13 = 0.;
    }

  }

  /* With no obtuse angle, we choose the true shortest diagonal */

  else {

    for (i = 0; i < dim; i++) {
      v1[i] = coords[vertex_id[2]*dim + i] - coords[vertex_id[0]*dim + i];
      v2[i] = coords[vertex_id[3]*dim + i] - coords[vertex_id[1]*dim + i];
    }

    /* Now compute diagonal lengths (squared) */

    d2_02 = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
    d2_13 = v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2];

  }

  /* Now define triangulation */

  if (quadrangle_vertices != NULL) {
    if (d2_02 < d2_13) {
      triangle_vertices[0] = quadrangle_vertices[0]; /* 1st triangle */
      triangle_vertices[1] = quadrangle_vertices[1];
      triangle_vertices[2] = quadrangle_vertices[2];
      triangle_vertices[3] = quadrangle_vertices[2]; /* 2nd triangle */
      triangle_vertices[4] = quadrangle_vertices[3];
      triangle_vertices[5] = quadrangle_vertices[0];
    }
    else {
      triangle_vertices[0] = quadrangle_vertices[0]; /* 1st triangle */
      triangle_vertices[1] = quadrangle_vertices[1];
      triangle_vertices[2] = quadrangle_vertices[3];
      triangle_vertices[3] = quadrangle_vertices[2]; /* 2nd triangle */
      triangle_vertices[4] = quadrangle_vertices[3];
      triangle_vertices[5] = quadrangle_vertices[1];
    }
  }
  else { /* if (quadrangle_vertices == NULL) */
    if (d2_02 < d2_13) {
      triangle_vertices[0] = 1; /* 1st triangle */
      triangle_vertices[1] = 2;
      triangle_vertices[2] = 3;
      triangle_vertices[3] = 3; /* 2nd triangle */
      triangle_vertices[4] = 4;
      triangle_vertices[5] = 1;
    }
    else {
      triangle_vertices[0] = 1; /* 1st triangle */
      triangle_vertices[1] = 2;
      triangle_vertices[2] = 4;
      triangle_vertices[3] = 3; /* 2nd triangle */
      triangle_vertices[4] = 4;
      triangle_vertices[5] = 2;
    }
  }

  /* Return number of triangles (for consistency with polygon triangulation) */

  return 2;
}

/*----------------------------------------------------------------------------
 * Updates the location[] and distance[] arrays associated with a set
 * of 1d points for points that are in a given element extent, based only
 * on this extent only (temporary, unoptimzed location).
 *
 * parameters:
 *   elt_num         <-- number of element corresponding to extents
 *   extents         <-> extents associated with element:
 *                       x_min, x_max (size: 2)
 *   n_points        <-- number of points to locate
 *   point_coords    <-- point coordinates
 *   location        <-> number of element containing or closest to each
 *                       point (size: n_points)
 *   distance        <-> distance from point to element indicated by
 *                       location[]: < 0 if unlocated, 0 - 1 if inside,
 *                       > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_by_extents_1d(ple_lnum_t         elt_num,
                      const double       extents[],
                      ple_lnum_t         n_points,
                      const ple_coord_t  point_coords[],
                      ple_lnum_t         location[],
                      float              distance[])
{
  ple_lnum_t  i;

  /* For now, we base a minimal location test on the element extents */
  /* The behavior is quadradic, nothing is optimized yet */

  for (i = 0; i < n_points; i++) {

    double elt_coord_max = -1;
    double elt_coord = -1;

    double cur_coord = point_coords[i];

    elt_coord =   (cur_coord - 0.5*(extents[1] + extents[0]))
                / (            0.5*(extents[1] - extents[0]));

    elt_coord = _ABS(elt_coord);

    if (elt_coord > elt_coord_max)
      elt_coord_max = elt_coord;

    if (  (distance[i] < 0 && elt_coord_max < 1)
        || elt_coord_max < distance[i]) {

      location[i] = elt_num;
      distance[i] = elt_coord_max;

    }

  }

}

/*----------------------------------------------------------------------------
 * Updates the location[] and distance[] arrays associated with a set
 * of points for points that are in a given element extent, based on
 * a query of points with this extent.
 *
 * parameters:
 *   elt_num            <-- number of element corresponding to extents
 *   dim                <-- spatial dimension
 *   extents            <-> extents associated with element:
 *                          x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *   point_coords       <-- point coordinates (size > n_points_in_extent*dim)
 *   n_points_in_extent <-- number of points in extents
 *   points_in_extent   <-- ids of points in extents
 *   location           <-> number of element containing or closest to each
 *                          point (size: n_points)
 *   distance           <-> distance from point to element indicated by
 *                          location[]: < 0 if unlocated, 0 - 1 if inside,
 *                          > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_in_extents(const ple_lnum_t    elt_num,
                   const int           dim,
                   const double        extents[],
                   const ple_coord_t   point_coords[],
                   ple_lnum_t          n_points_in_extents,
                   const ple_lnum_t    points_in_extents[],
                   ple_lnum_t          location[],
                   float               distance[])
{
  ple_lnum_t  i, j, k;

  /* For now, we base a minimal location test on the element extents */
  /* The behavior is quadradic, nothing is optimized yet */

  for (i = 0; i < n_points_in_extents; i++) {

    double elt_coord_max = -1;
    double elt_coord = -1;

    j = points_in_extents[i];

    for (k = 0; k < dim; k++) {

      double cur_coord = point_coords[j*dim + k];

      elt_coord =   (cur_coord - 0.5*(extents[k+dim] + extents[k]))
                  / (            0.5*(extents[k+dim] - extents[k]));

      elt_coord = _ABS(elt_coord);

      if (elt_coord > elt_coord_max)
        elt_coord_max = elt_coord;

    }

    if (  (distance[j] < 0 && elt_coord_max < 1)
        || elt_coord_max < distance[j]) {

      location[j] = elt_num;
      distance[j] = elt_coord_max;

    }

  }

}

/*----------------------------------------------------------------------------
 * Build a local octree's leaves.
 *
 * parameters:
 *   level              <-- current level in octree
 *   extents            <-- extents associated with node:
 *                          x_min, y_min, z_min, x_max, y_max, z_max (size: 6)
 *   point_coords       <-- point coordinates
 *   point_ids_tmp      <-- temporary point indexes
 *   octree             <-> current octree structure
 *   point_range        <-> start and past-the end index in point_idx
 *                          for current node (size: 2)
 *----------------------------------------------------------------------------*/

static void
_build_octree_leaves(int                 level,
                     const double        extents[],
                     const ple_coord_t   point_coords[],
                     ple_lnum_t         *point_ids_tmp,
                     _octree_t          *octree,
                     ple_lnum_t          point_range[2])
{
  ple_lnum_t i, j, k, _n_nodes, _n_points, tmp_size;

  ple_lnum_t count[8], idx[9], octant_id[8];
  double mid[3], sub_extents[6];
  _octant_t  *_node;

  int octant_mask[3] = {4, 2, 1}; /* pow(2, 2), pow(2, 1), pow(2,0) */

  _n_nodes = octree->n_nodes;
  tmp_size = octree->n_nodes;

  /* Resize octree if necesary */

  if (octree->n_nodes >= octree->n_nodes_max) {
    if (octree->n_nodes == 0) {
      octree->n_nodes = 1;
      octree->n_nodes_max = 8;
    }
    octree->n_nodes_max *= 2;
    PLE_REALLOC(octree->nodes, octree->n_nodes_max, _octant_t);
  }

  /* Number of points */

  _n_points = point_range[1] - point_range[0];

  /* Extents center */

  for (j = 0; j < 3; j++)
    mid[j]= (extents[j] + extents[j + 3]) * 0.5;

  for (j = 0; j < 8; j++) {
    count[j] = 0;
    octant_id[j] = -1;
  }

  /* Count points in each octant */

  for (i = point_range[0]; i < point_range[1]; i++) {

    for (j = 0, k = 0; j < 3; j++) {
      if (point_coords[octree->point_ids[i]*3 + j] > mid[j])
        k += octant_mask[j];
    }

    count[k] += 1;
  }

  /* Build index */

  idx[0] = 0;
  for (j = 0; j < 8; j++)
    idx[j+1] = idx[j] + count[j];

  for (j = 0; j < 8; j++)
    count[j] = 0;

  for (i = point_range[0], j = 0; i < point_range[1]; i++) {

    for (j = 0, k = 0; j < 3; j++) {
      if (point_coords[octree->point_ids[i]*3 + j] > mid[j])
        k += octant_mask[j];
    }

    point_ids_tmp[idx[k] + count[k]] = octree->point_ids[i];
    count[k] += 1;
  }

  for (i = point_range[0], j = 0; i < point_range[1]; i++, j++)
    octree->point_ids[i] = point_ids_tmp[j];

  for (i = 0; i < 9; i++)
    idx[i] = point_range[0] + idx[i];

  /* Build leaves recursively */

  if (level < _octree_max_levels) {

    for (i = 0; i < 8; i++) {

      if (count[i] > _octree_threshold) {

        tmp_size++;

        octant_id[i] = tmp_size;

        if (i < 4) {
          sub_extents[0] = extents[0];
          sub_extents[3] = mid[0];
        }
        else {
          sub_extents[0] = mid[0];
          sub_extents[3] = extents[3];
        }
        /* 1.0e-12 term in assert() used to allow for
           truncation error in for xmin = xmax case */
        assert(sub_extents[0] < sub_extents[3] + 1.0e-12);

        if (i%4 < 2) {
          sub_extents[1] = extents[1];
          sub_extents[4] = mid[1];
        }
        else {
          sub_extents[1] = mid[1];
          sub_extents[4] = extents[4];
        }
        assert(sub_extents[1] < sub_extents[4] + 1.0e-12);

        if (i%2 < 1) {
          sub_extents[2] = extents[2];
          sub_extents[5] = mid[2];
        }
        else {
          sub_extents[2] = mid[2];
          sub_extents[5] = extents[5];
        }
        assert(sub_extents[2] < sub_extents[5] + 1.0e-12);

        octree->n_nodes = tmp_size;

        _build_octree_leaves(level + 1,
                             sub_extents,
                             point_coords,
                             point_ids_tmp,
                             octree,
                             idx + i);

        tmp_size = octree->n_nodes;
      }

    }

  }

  /* Finalize node */

  _node = octree->nodes + _n_nodes;

  for (i = 0; i < 9; i++)
    _node->idx[i] = idx[i];

  for (i = 0; i < 8; i++)
    _node->octant_id[i] = octant_id[i];

  _node->n_points = _n_points;
}

/*----------------------------------------------------------------------------
 * Build an octree structure to locate 3d points in mesh.
 *
 * parameters:
 *   n_points        <-- number of points to locate
 *   point_coords    <-- point coordinates
 *
 * returns:
 *   pointer to local octree structure
 *----------------------------------------------------------------------------*/

static _octree_t
_build_octree(ple_lnum_t         n_points,
              const ple_coord_t  point_coords[])
{
  size_t i;
  ple_lnum_t point_range[2];
  _octree_t _octree;

  int *point_ids_tmp = NULL;

  /* Initialization */

  point_range[0] = 0;
  point_range[1] = n_points;

  _octree.n_points = n_points;
  _octree.n_nodes = 0;
  _octree.n_nodes_max = 0;
  _octree.nodes = NULL;
  _octree.point_ids = NULL;

  if (n_points > 0) {

    _point_extents(3,
                   n_points,
                   NULL,
                   point_coords,
                   _octree.extents);

    PLE_MALLOC(_octree.point_ids, _octree.n_points, ple_lnum_t);

    for (i = 0; i < _octree.n_points; i++)
      _octree.point_ids[i] = i;

    PLE_MALLOC(point_ids_tmp, n_points, int);

    _build_octree_leaves(0,
                         _octree.extents,
                         point_coords,
                         point_ids_tmp,
                         &_octree,
                         point_range);

    PLE_FREE(point_ids_tmp);

  }

  return _octree;
}

/*----------------------------------------------------------------------------
 * Free an octree structure.
 *
 * parameters:
 *   octree <-> octree structure whose elements are to be freed
 *
 * returns:
 *   pointer to local octree structure
 *----------------------------------------------------------------------------*/

static void
_free_octree(_octree_t *octree)
{

  octree->n_points = 0;
  octree->n_nodes = 0;
  octree->n_nodes_max = 0;

  PLE_FREE(octree->nodes);
  PLE_FREE(octree->point_ids);
}

/*----------------------------------------------------------------------------
 * locate points in box defined by extents in an octant.
 *
 * parameters:
 *   search_extents  <-- extents associated with element:
 *   point_coords    <-- point coordinates
 *                       x_min, y_min, z_min, x_max, y_max, z_max (size: 6)
 *   octree          <-- point octree
 *   node_extents    <-- extents of octant
 *   node_id         <-- id of node in octree
 *   loc_point_ids   --> ids of located points (max size: octree->n_points)
 *   n_loc_points    --> number of points located
 *----------------------------------------------------------------------------*/

static void
_query_octree_node(const double        extents[],
                   const ple_coord_t   point_coords[],
                   const _octree_t    *octree,
                   const double        node_extents[],
                   int                 node_id,
                   ple_lnum_t         *loc_point_ids,
                   ple_lnum_t         *n_loc_points)
{
  _octant_t *node;
  int i, j, k;
  int dim = 3;
  double sub_extents[6], mid[3];

  node = octree->nodes + node_id;

  if (_intersect_extents(dim, node_extents, extents)) {

    for (j = 0; j < dim; j++)
      mid[j]= (node_extents[j] + node_extents[j + dim]) * 0.5;

    /* Loop on node leaves */

    for (i = 0; i < 8; i++) {

      /* Compute octant extents */

      if (i < 4) {
        sub_extents[0] = node_extents[0];
        sub_extents[3] = mid[0];
      }
      else {
        sub_extents[0] = mid[0];
        sub_extents[3] = node_extents[3];
      }

      if (i%4 < 2) {
        sub_extents[1] = node_extents[1];
        sub_extents[4] = mid[1];
      }
      else {
        sub_extents[1] = mid[1];
        sub_extents[4] = node_extents[4];
      }
      if (i%2 < 1) {
        sub_extents[2] = node_extents[2];
        sub_extents[5] = mid[2];
      }
      else {
        sub_extents[2] = mid[2];
        sub_extents[5] = node_extents[5];
      }

      /* Search recursively if octant is not final */

      if (node->octant_id[i] > -1)
        _query_octree_node(extents,
                           point_coords,
                           octree,
                           sub_extents,
                           node->octant_id[i],
                           loc_point_ids,
                           n_loc_points);

      else {

        /* Update list of points located */

        if (_intersect_extents(dim, sub_extents, extents)) {

          for (k = node->idx[i]; k < node->idx[i+1]; k++) {

            ple_lnum_t point_id = octree->point_ids[k];

            if (_within_extents(dim, point_coords + point_id*dim, extents)) {
              loc_point_ids[*n_loc_points] = point_id;
              (*n_loc_points)++;
            }

          }
        }

      }

    } /* End of loop on node leaves */

  }

}

/*----------------------------------------------------------------------------
 * locate points in box defined by extents using an octree.
 *
 * parameters:
 *   search_extents  <-- extents associated with element:
 *                       x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *   point_coords    <-- point coordinates
 *   octree          <-- point octree
 *   n_loc_points    --> number of points located
 *   loc_point_ids   --> ids of located points (max size: octree->n_points)
 *----------------------------------------------------------------------------*/

static void
_query_octree(const double        extents[],
              const ple_coord_t   point_coords[],
              const _octree_t    *octree,
              ple_lnum_t         *n_loc_points,
              ple_lnum_t          loc_point_ids[])
{
  *n_loc_points = 0;

  if (octree->n_points > 0)
    _query_octree_node(extents,
                       point_coords,
                       octree,
                       octree->extents,
                       0,
                       loc_point_ids,
                       n_loc_points);
}

/*----------------------------------------------------------------------------
 * Build a local quadtree's leaves.
 *
 * parameters:
 *   level              <-- current level in octree
 *   extents            <-- extents associated with node:
 *                          x_min, y_min, x_max, y_max, (size: 4)
 *   point_coords       <-- point coordinates
 *   point_ids_tmp      <-- temporary point indexes
 *   pos_tmp            <-- temporary point position in quadtree
 *   octree             <-> current quadtree structure
 *   point_range        <-> start and past-the end index in point_idx
 *                          for current node (size: 2)
 *----------------------------------------------------------------------------*/

static void
_build_quadtree_leaves(int                 level,
                       const double        extents[],
                       const ple_coord_t   point_coords[],
                       ple_lnum_t         *point_ids_tmp,
                       _quadtree_t        *quadtree,
                       ple_lnum_t          point_range[2])
{
  ple_lnum_t i, j, k, _n_nodes, _n_points, tmp_size;

  ple_lnum_t count[4], idx[5], quadrant_id[4];
  double mid[2], sub_extents[4];
  _quadrant_t  *_node;

  int quadrant_mask[2] = {2, 1}; /* pow(2, 1), pow(2,0) */

  _n_nodes = quadtree->n_nodes;
  tmp_size = quadtree->n_nodes;

  /* Resize quadtree if necesary */

  if (quadtree->n_nodes >= quadtree->n_nodes_max) {
    if (quadtree->n_nodes == 0) {
      quadtree->n_nodes = 1;
      quadtree->n_nodes_max = 4;
    }
    quadtree->n_nodes_max *= 2;
    PLE_REALLOC(quadtree->nodes, quadtree->n_nodes_max, _quadrant_t);
  }

  /* Number of points */

  _n_points = point_range[1] - point_range[0];

  /* Extents center */

  for (j = 0; j < 2; j++)
    mid[j]= (extents[j] + extents[j + 2]) * 0.5;

  for (j = 0; j < 4; j++) {
    count [j] = 0;
    quadrant_id[j] = -1;
  }

  /* Count points in each quadrant */

  for (i = point_range[0]; i < point_range[1]; i++) {

    for (j = 0, k = 0; j < 2; j++) {
      if (point_coords[quadtree->point_ids[i]*2 + j] > mid[j])
        k += quadrant_mask[j];
    }

    count[k] += 1;
  }

  /* Build index */

  idx[0] = 0;
  for (j = 0; j < 4; j++)
    idx[j+1] = idx[j] + count[j];

  for (j = 0; j < 4; j++)
    count[j] = 0;

  for (i = point_range[0], j = 0; i < point_range[1]; i++) {

    for (j = 0, k = 0; j < 2; j++) {
      if (point_coords[quadtree->point_ids[i]*2 + j] > mid[j])
        k += quadrant_mask[j];
    }

    point_ids_tmp[idx[k] + count[k]] = quadtree->point_ids[i];
    count[k] += 1;
  }

  for (i = point_range[0], j = 0; i < point_range[1]; i++, j++)
    quadtree->point_ids[i] = point_ids_tmp[j];

  for (i = 0; i < 5; i++)
    idx[i] = point_range[0] + idx[i];

  /* Build leaves recursively */

  if (level < _octree_max_levels) {

    for (i = 0; i < 4; i++) {

      if (count[i] > _octree_threshold) {

        tmp_size++;

        quadrant_id[i] = tmp_size;

        if (i < 2) {
          sub_extents[0] = extents[0];
          sub_extents[2] = mid[0];
        }
        else {
          sub_extents[0] = mid[0];
          sub_extents[2] = extents[2];
        }
        assert(sub_extents[0] < sub_extents[2] + 1.0e-12);

        if (i%2 < 1) {
          sub_extents[1] = extents[1];
          sub_extents[3] = mid[1];
        }
        else {
          sub_extents[1] = mid[1];
          sub_extents[3] = extents[3];
        }
        assert(sub_extents[1] < sub_extents[3] + 1.0e-12);

        quadtree->n_nodes = tmp_size;

        _build_quadtree_leaves(0,
                               sub_extents,
                               point_coords,
                               point_ids_tmp,
                               quadtree,
                               idx + i);

        tmp_size = quadtree->n_nodes;
      }

    }

  }

  /* Finalize node */

  _node = quadtree->nodes + _n_nodes;

  for (i = 0; i < 5; i++)
    _node->idx[i] = idx[i];

  for (i = 0; i < 4; i++)
    _node->quadrant_id[i] = quadrant_id[i];

  _node->n_points = _n_points;
}

/*----------------------------------------------------------------------------
 * Build a quadtree structure to locate 2d points in mesh.
 *
 * parameters:
 *   n_points        <-- number of points to locate
 *   point_coords    <-- point coordinates
 *
 * returns:
 *   pointer to local quadtree structure
 *----------------------------------------------------------------------------*/

static _quadtree_t
_build_quadtree(ple_lnum_t         n_points,
                const ple_coord_t  point_coords[])
{
  size_t i;
  ple_lnum_t point_range[2];
  _quadtree_t _quadtree;

  int *point_ids_tmp = NULL;

  /* Initialization */

  point_range[0] = 0;
  point_range[1] = n_points;

  _quadtree.n_points = n_points;
  _quadtree.n_nodes = 0;
  _quadtree.n_nodes_max = 0;
  _quadtree.nodes = NULL;
  _quadtree.point_ids = NULL;

  if (n_points > 0) {

    _point_extents(2,
                   n_points,
                   NULL,
                   point_coords,
                   _quadtree.extents);

    PLE_MALLOC(_quadtree.point_ids, _quadtree.n_points, ple_lnum_t);

    for (i = 0; i < _quadtree.n_points; i++)
      _quadtree.point_ids[i] = i;

    PLE_MALLOC(point_ids_tmp, n_points, int);

    _build_quadtree_leaves(0,
                           _quadtree.extents,
                           point_coords,
                           point_ids_tmp,
                           &_quadtree,
                           point_range);

    PLE_FREE(point_ids_tmp);

  }

  return _quadtree;
}

/*----------------------------------------------------------------------------
 * Free a quadtree structure.
 *
 * parameters:
 *   quadtree <-> quadtree structure whose elements are to be freed
 *
 * returns:
 *   pointer to local quadtree structure
 *----------------------------------------------------------------------------*/

static void
_free_quadtree(_quadtree_t *quadtree)
{

  quadtree->n_points = 0;
  quadtree->n_nodes = 0;
  quadtree->n_nodes_max = 0;

  PLE_FREE(quadtree->nodes);
  PLE_FREE(quadtree->point_ids);
}

/*----------------------------------------------------------------------------
 * locate points in box defined by extents in a quadrant.
 *
 * parameters:
 *   search_extents  <-- extents associated with element:
 *   point_coords    <-- point coordinates
 *                       x_min, y_min, x_max, y_max (size: 4)
 *   quadtree        <-- point quadtree
 *   node_extents    <-- extents of quadrant
 *   node_id         <-- if of node in octree
 *   loc_point_ids   --> ids of located points (max size: quadtree->n_points)
 *   n_loc_points    --> number of points located
 *----------------------------------------------------------------------------*/

static void
_query_quadtree_node(const double        extents[],
                     const ple_coord_t   point_coords[],
                     const _quadtree_t  *quadtree,
                     const double        node_extents[],
                     int                 node_id,
                     ple_lnum_t         *loc_point_ids,
                     ple_lnum_t         *n_loc_points)
{
  _quadrant_t *node;
  int i, j, k;
  double sub_extents[4], mid[2];

  node = quadtree->nodes + node_id;

  if (_intersect_extents(2, node_extents, extents)) {

    for (j = 0; j < 2; j++)
      mid[j]= (node_extents[j] + node_extents[j + 2]) * 0.5;

    /* Loop on node leaves */

    for (i = 0; i < 4; i++) {

      /* Compute quadrant extents */

      if (i < 2) {
        sub_extents[0] = node_extents[0];
        sub_extents[2] = mid[0];
      }
      else {
        sub_extents[0] = mid[0];
        sub_extents[2] = node_extents[2];
      }

      if (i%2 < 1) {
        sub_extents[1] = node_extents[1];
        sub_extents[3] = mid[1];
      }
      else {
        sub_extents[1] = mid[1];
        sub_extents[3] = node_extents[3];
      }

      /* Search recursively if quadrant is not final */

      if (node->quadrant_id[i] > -1)
        _query_quadtree_node(extents,
                             point_coords,
                             quadtree,
                             sub_extents,
                             node->quadrant_id[i],
                             loc_point_ids,
                             n_loc_points);

      else {

        /* Update list of points located */

        if (_intersect_extents(2, sub_extents, extents)) {

          for (k = node->idx[i]; k < node->idx[i+1]; k++) {

            ple_lnum_t point_id = quadtree->point_ids[k];

            if (_within_extents(2, point_coords + point_id*2, extents)) {
              loc_point_ids[*n_loc_points] = point_id;
              (*n_loc_points)++;
            }

          }
        }

      }

    } /* End of loop on node leaves */

  }

}

/*----------------------------------------------------------------------------
 * locate points in box defined by extents using an octree.
 *
 * parameters:
 *   search_extents  <-- extents associated with element:
 *   point_coords    <-- point coordinates
 *                       x_min, y_min, x_max, y_max (size: 4)
 *   quadtree        <-- point quadtree
 *   n_loc_points    --> number of points located
 *   loc_point_ids   --> ids of located points (max size: octree->n_points)
 *----------------------------------------------------------------------------*/

static void
_query_quadtree(const double        extents[],
                const ple_coord_t   point_coords[],
                const _quadtree_t  *quadtree,
                ple_lnum_t         *n_loc_points,
                ple_lnum_t          loc_point_ids[])
{
  *n_loc_points = 0;

  if (quadtree->n_points > 0)
    _query_quadtree_node(extents,
                         point_coords,
                         quadtree,
                         quadtree->extents,
                         0,
                         loc_point_ids,
                         n_loc_points);
}

/*----------------------------------------------------------------------------
 * Locate points on a 3d edge, and update the location[] and distance[]
 * arrays associated with the point set.
 *
 * parameters:
 *   elt_num             <-- number of element corresponding to extents
 *   element_vertex_num  <-- element vertex numbers
 *   vertex_coords       <-- pointer to vertex coordinates
 *   point_coords        <-- point coordinates
 *   n_points_in_extents <-- number of points in extents
 *   points_in_extents   <-- ids of points in extents
 *   tolerance           <-- associated tolerance (considered infinite if < 0)
 *   location            <-> number of element containing or closest to each
 *                           point (size: n_points)
 *   distance            <-> distance from point to element indicated by
 *                           location[]: < 0 if unlocated, absolute distance
 *                           to element if located (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_on_edge_3d(ple_lnum_t           elt_num,
                   const ple_lnum_t     element_vertex_num[],
                   const ple_coord_t    vertex_coords[],
                   const ple_coord_t    point_coords[],
                   ple_lnum_t           n_points_in_extents,
                   const ple_lnum_t     points_in_extents[],
                   double               tolerance,
                   ple_lnum_t           location[],
                   float                distance[])
{
  ple_lnum_t  i, j, k, coord_idx_0, coord_idx_1;

  double u[3], v[3];
  double uv, len2, isop_0;
  double dist2, epsilon2, vertex_dist2;

  /* vertex index of the edge studied */

  coord_idx_0 = element_vertex_num[0] - 1;
  coord_idx_1 = element_vertex_num[1] - 1;

  /* Calculate edge vector and length */

  for (j = 0; j < 3; j++)
    u[j] =   vertex_coords[(coord_idx_1*3) + j]
           - vertex_coords[(coord_idx_0*3) + j];

  len2 = _DOT_PRODUCT_3D(u, u);

  if (len2 < _epsilon_denom)
    return;

  else if (tolerance < 0.0)
    epsilon2 = HUGE_VAL;

  else
    epsilon2 = len2*tolerance*tolerance;

  /* Loop on points resulting from extent query */

  for (k = 0; k < n_points_in_extents; k++) {

    i =  points_in_extents[k];

    vertex_dist2 = distance[i]*distance[i];

    /* Calculate linear coordinates of projection of point on edge axis */

    for (j = 0; j < 3; j++)
      v[j] = point_coords[i*3 + j] - vertex_coords[(coord_idx_0*3) + j];

    uv = _DOT_PRODUCT_3D(u, v);

    isop_0 = uv / len2;

    /* Set v to be the vector from the point to the closest point on
       the segment (if isop_0 < 0, v is already that vector) */

    if (isop_0 >= 1) {
      for (j = 0; j < 3; j++)
        v[j] = point_coords[i*3 + j] - vertex_coords[(coord_idx_1*3) + j];
    }
    else if (isop_0 > 0) {
      for (j = 0; j < 3; j++)
        v[j] -= isop_0*u[j];
    }

    /* Distance between point to locate and its projection */

    dist2 = _DOT_PRODUCT_2D(v, v);

    if (dist2 < epsilon2 && (dist2 < vertex_dist2 || distance[i] < 0.0)) {
      location[i] = elt_num;
      distance[i] = sqrt(dist2);
    }

  } /* End of loop on points resulting from extent query */

}

/*----------------------------------------------------------------------------
 * Locate points on a 2d edge, and update the location[] and distance[]
 * arrays associated with the point set.
 *
 * parameters:
 *   elt_num             <-- number of element corresponding to extents
 *   element_vertex_num  <-- element vertex numbers
 *   vertex_coords       <-- pointer to vertex coordinates
 *   point_coords        <-- point coordinates
 *   n_points_in_extent  <-- number of points in extents
 *   points_in_extent    <-- ids of points in extents
 *   n_points_in_extent  <-- number of points in extents
 *   points_in_extent    <-- ids of points in extents
 *   tolerance           <-- associated tolerance (considered infinite if < 0)
 *   location            <-> number of element containing or closest to each
 *                           point (size: n_points)
 *   distance            <-> distance from point to element indicated by
 *                           location[]: < 0 if unlocated, absolute distance
 *                           to element if located (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_on_edge_2d(ple_lnum_t           elt_num,
                   const ple_lnum_t     element_vertex_num[],
                   const ple_coord_t    vertex_coords[],
                   const ple_coord_t    point_coords[],
                   ple_lnum_t           n_points_in_extents,
                   const ple_lnum_t     points_in_extents[],
                   double               tolerance,
                   ple_lnum_t           location[],
                   float                distance[])
{
  ple_lnum_t  i, j, k, coord_idx_0, coord_idx_1;

  double u[2], v[2];
  double uv, len2, isop_0;
  double dist2, epsilon2, vertex_dist2;

  /* vertex index of the edge studied */

  coord_idx_0 = element_vertex_num[0] - 1;
  coord_idx_1 = element_vertex_num[1] - 1;

  /* Calculate edge vector and length */

  for (j = 0; j < 2; j++)
    u[j] =   vertex_coords[(coord_idx_1*2) + j]
           - vertex_coords[(coord_idx_0*2) + j];

  len2 = _DOT_PRODUCT_2D(u, u);

  if (len2 < _epsilon_denom)
    return;

  else if (tolerance < 0.0)
    epsilon2 = HUGE_VAL;

  else
    epsilon2 = len2*tolerance*tolerance;

  /* Loop on points resulting from extent query */

  for (k = 0; k < n_points_in_extents; k++) {

    i =  points_in_extents[k];

    vertex_dist2 = distance[i]*distance[i];

    /* Calculate linear coordinates of projection of point on edge axis */

    for (j = 0; j < 2; j++)
      v[j] = point_coords[i*2 + j] - vertex_coords[(coord_idx_0*2) + j];

    uv = u[0]*v[0] + u[1]*v[1];

    isop_0 = uv / len2;

    /* Set v to be the vector from the point to the closest point on
       the segment (if isop_0 < 0, v is already that vector) */

    if (isop_0 >= 1) {
      for (j = 0; j < 2; j++)
        v[j] = point_coords[i*2 + j] - vertex_coords[(coord_idx_1*2) + j];
    }
    else if (isop_0 > 0) {
      for (j = 0; j < 2; j++)
        v[j] -= isop_0*u[j];
    }

    /* Distance between point to locate and its projection */

    dist2 = _DOT_PRODUCT_2D(v, v);

    if (dist2 < epsilon2 && (dist2 < vertex_dist2 || distance[i] < 0.0)) {
      location[i] = elt_num;
      distance[i] = sqrt(dist2);
    }

  } /* End of loop on points resulting from extent query */

}

/*----------------------------------------------------------------------------
 * Locate points in a given set of 3d triangles, and updates the location[]
 * and distance[] arrays associated with the point set.
 *
 * This function is called for sets of triangles belonging to the subdivision
 * of a given 3d face. Barycentric coordinates are used to locate the
 * projection of points.
 *
 * parameters:
 *   elt_num             <-- number of element corresponding to extents
 *   n_triangles         <-- number of triangles
 *   triangle_vertices   <-- triangles connectivity; size: 2 * 3
 *   vertex_coords       <-- pointer to vertex coordinates
 *   point_coords        <-- point coordinates
 *   n_points_in_extents <-- number of points in extents
 *   points_in_extents   <-- ids of points in extents
 *   tolerance           <-- associated tolerance (considered infinite if < 0)
 *   location            <-> number of element containing or closest to each
 *                           point (size: n_points)
 *   distance            <-> distance from point to element indicated by
 *                           location[]: < 0 if unlocated, absolute distance
 *                           to element if located (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_on_triangles_3d(ple_lnum_t           elt_num,
                        int                  n_triangles,
                        const ple_lnum_t     triangle_vertices[],
                        const ple_coord_t    vertex_coords[],
                        const ple_coord_t    point_coords[],
                        ple_lnum_t           n_points_in_extents,
                        const ple_lnum_t     points_in_extents[],
                        const double         tolerance,
                        ple_lnum_t           location[],
                        float                distance[])
{
  ple_lnum_t  i, j, k, tria_id, coord_idx_0, coord_idx_1, coord_idx_2;

  double t[3], u[3], v[3], w[3], vect_tmp[3];
  double uu, vv, uv, ut, vt, ww, det, tmp_max;
  double epsilon2, dist2, vertex_dist2, isop_0, isop_1;

  double tolerance2 = tolerance*tolerance;

  /* Loop on element's sub-triangles */

  for (tria_id = 0; tria_id < n_triangles; tria_id++) {

    /* vertex index of the triangle studied */

    coord_idx_0 = triangle_vertices[tria_id*3]     - 1;
    coord_idx_1 = triangle_vertices[tria_id*3 + 1] - 1;
    coord_idx_2 = triangle_vertices[tria_id*3 + 2] - 1;

    /* Calculate triangle-constant values for barycentric coordinates */

    for (j = 0; j < 3; j++) {
      u[j] = - vertex_coords[(coord_idx_0*3) + j]
             + vertex_coords[(coord_idx_1*3) + j];
      v[j] = - vertex_coords[(coord_idx_0*3) + j]
             + vertex_coords[(coord_idx_2*3) + j];
      w[j] =   vertex_coords[(coord_idx_1*3) + j]
             - vertex_coords[(coord_idx_2*3) + j];
    }

    uu = _DOT_PRODUCT_3D(u, u);
    vv = _DOT_PRODUCT_3D(v, v);
    ww = _DOT_PRODUCT_3D(w, w);
    uv = _DOT_PRODUCT_3D(u, v);

    det = (uu*vv - uv*uv);

    if (det < _epsilon_denom)
      continue;

    /* epsilon2 is based on maximum edge length (squared) */

    tmp_max = _MAX(vv, ww);

    if (tolerance < 0.)
      epsilon2 = HUGE_VAL;
    else
      epsilon2 = _MAX(uu, tmp_max) * tolerance2;

    /* Loop on points resulting from extent query */

    for (k = 0; k < n_points_in_extents; k++) {

      i =  points_in_extents[k];

      vertex_dist2 = distance[i]*distance[i];

      /* Calculation of the barycenter coordinates for the projected node */

      for (j = 0; j < 3; j++)
        t[j] = - vertex_coords[(coord_idx_0*3) + j]
               + point_coords[i*3 + j];

      ut = _DOT_PRODUCT_3D(u, t);
      vt = _DOT_PRODUCT_3D(v, t);

      isop_0 = (ut*vv - vt*uv) / det;
      isop_1 = (uu*vt - uv*ut) / det;

      _CROSS_PRODUCT_3D(vect_tmp, u, v);

      /* if the projected point is not on triangle, we project it
         on the nearest edge or node */

      if (isop_0 < 0.)
        isop_0 = 0.;

      if (isop_1 < 0.)
        isop_1 = 0.;

      if ((1.0 - isop_0 - isop_1) < 0.) {
        isop_0 = isop_0 / (isop_0 + isop_1);
        isop_1 = isop_1 / (isop_0 + isop_1);
      }

      /* re-use vect_tmp */

      for (j = 0; j < 3; j++)
        vect_tmp[j] =   vertex_coords[coord_idx_0*3 + j]
                      + u[j]*isop_0
                      + v[j]*isop_1
                      - point_coords[i*3 + j];

      /* distance between point to locate and its projection */

      dist2 = _DOT_PRODUCT_3D(vect_tmp, vect_tmp);

      if (dist2 < epsilon2 && (dist2 < vertex_dist2 || distance[i] < 0.0)) {
        location[i] = elt_num;
        distance[i] = sqrt(dist2);
      }

    } /* End of loop on points resulting from extent query */

  } /* End of loop on element's sub-triangles */

}

/*----------------------------------------------------------------------------
 * Locate points in a given set of 2d triangles, and updates the location[]
 * and distance[] arrays associated with the point set.
 *
 * This function is called for sets of triangles belonging to the subdivision
 * of a given 2d face. Barycentric coordinates are used to locate the
 * projection of points.
 *
 * parameters:
 *   elt_num             <-- number of element corresponding to extents
 *   n_triangles         <-- number of triangles
 *   triangle_vertices   <-- triangles connectivity; size: 2 * 3
 *   vertex_coords       <-- pointer to vertex coordinates
 *   point_coords        <-- point coordinates
 *   n_points_in_extent  <-- number of points in extents
 *   points_in_extent    <-- ids of points in extents
 *   n_points_in_extent  <-- number of points in extents
 *   points_in_extent    <-- ids of points in extents
 *   tolerance           <-- associated tolerance
 *   location            <-> number of element containing or closest to each
 *                           point (size: n_points)
 *   distance            <-> distance from point to element indicated by
 *                           location[]: < 0 if unlocated, 0 - 1 if inside,
 *                           > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_on_triangles_2d(ple_lnum_t           elt_num,
                        int                  n_triangles,
                        const ple_lnum_t     triangle_vertices[],
                        const ple_coord_t    vertex_coords[],
                        const ple_coord_t    point_coords[],
                        ple_lnum_t           n_points_in_extents,
                        const ple_lnum_t     points_in_extents[],
                        double               tolerance,
                        ple_lnum_t           location[],
                        float                distance[])
{
  ple_lnum_t  i, j, k, tria_id, coord_idx_0, coord_idx_1, coord_idx_2;

  double t[2], u[2], v[2], shapef[3];
  double uu, vv, uv, ut, vt, det;
  double dist, max_dist, isop_0, isop_1;

  /* Loop on element's sub-triangles */

  for (tria_id = 0; tria_id < n_triangles; tria_id++) {

    /* vertex index of the triangle studied */

    coord_idx_0 = triangle_vertices[tria_id*3]     - 1;
    coord_idx_1 = triangle_vertices[tria_id*3 + 1] - 1;
    coord_idx_2 = triangle_vertices[tria_id*3 + 2] - 1;

    /* Calculate triangle-constant values for barycentric coordinates */

    for (j = 0; j < 2; j++) {
      u[j] = - vertex_coords[(coord_idx_0*2) + j]
             + vertex_coords[(coord_idx_1*2) + j];
      v[j] = - vertex_coords[(coord_idx_0*2) + j]
             + vertex_coords[(coord_idx_2*2) + j];
    }

    uu = _DOT_PRODUCT_2D(u, u);
    vv = _DOT_PRODUCT_2D(v, v);
    uv = _DOT_PRODUCT_2D(u, v);

    det = (uu*vv - uv*uv);

    if (det < _epsilon_denom)
      continue;

    /* Loop on points resulting from extent query */

    for (k = 0; k < n_points_in_extents; k++) {

      i =  points_in_extents[k];

      /* Calculation of the barycenter coordinates for the projected node */

      for (j = 0; j < 2; j++)
        t[j] = - vertex_coords[(coord_idx_0*2) + j]
               + point_coords[i*2 + j];

      ut = u[0]*t[0] + u[1]*t[1];
      vt = v[0]*t[0] + v[1]*t[1];

      isop_0 = (ut*vv - vt*uv) / det;
      isop_1 = (uu*vt - uv*ut) / det;

      shapef[0] = 1. - isop_0 - isop_1;
      shapef[1] =      isop_0;
      shapef[2] =               isop_1;

      max_dist = -1.0;

      for (j = 0; j < 3; j++){

        dist = 2.*_ABS(shapef[j] - 0.5);

        if (max_dist < dist)
          max_dist = dist;
      }

      if (   (max_dist > -0.5 && max_dist < (1. + 2.*tolerance))
          && (max_dist < distance[i] || distance[i] < 0)) {
        location[i] = elt_num;
        distance[i] = max_dist;
      }

    } /* End of loop on points resulting from extent query */

  } /* End of loop on element's sub-triangles */

}

/*----------------------------------------------------------------------------
 * Locate points in a tetrahedron whose coordinates are pre-computed,
 * updating the location[] and distance[] arrays associated with a set
 * of points.
 *
 * parameters:
 *   elt_num             <-- element number
 *   tetra_coords[]      <-- tetrahedra vertex coordinates
 *   point_coords        <-- point coordinates
 *   n_points_in_extents <-- number of points in element extents
 *   points_in_extents   <-- ids of points in extents
 *   tolerance           <-- associated tolerance
 *   location            <-> number of element containing or closest to each
 *                           point (size: n_points)
 *   distance            <-> distance from point to element indicated by
 *                           location[]: < 0 if unlocated, 0 - 1 if inside,
 *                           > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_in_tetra(ple_lnum_t         elt_num,
                 ple_coord_t        tetra_coords[4][3],
                 const ple_coord_t  point_coords[],
                 ple_lnum_t         n_points_in_extents,
                 const ple_lnum_t   points_in_extents[],
                 double             tolerance,
                 ple_lnum_t         location[],
                 float              distance[])
{
  double vol6;
  double dist, max_dist;
  int i, j, k;

  double isop_0, isop_1, isop_2;
  double t00, t10, t20, t01, t02, t03, t11, t12, t13, t21, t22, t23;
  double v01[3], v02[3], v03[3], shapef[4];

  for (i = 0; i < 3; i++) {
    v01[i] = tetra_coords[1][i] - tetra_coords[0][i];
    v02[i] = tetra_coords[2][i] - tetra_coords[0][i];
    v03[i] = tetra_coords[3][i] - tetra_coords[0][i];
  }

  vol6 = fabs(  v01[0] * (v02[1]*v03[2] - v02[2]*v03[1])
              - v02[0] * (v01[1]*v03[2] - v01[2]*v03[1])
              + v03[0] * (v01[1]*v02[2] - v01[2]*v02[1]));

  if (vol6 < _epsilon_denom)
    return;

  for (k = 0; k < n_points_in_extents; k++) {

    i = points_in_extents[k];

    t00  =   point_coords[i*3]     - tetra_coords[0][0];
    t10  =   point_coords[i*3 + 1] - tetra_coords[0][1];
    t20  =   point_coords[i*3 + 2] - tetra_coords[0][2];

    t01  = - tetra_coords[0][0] + tetra_coords[1][0];
    t02  = - tetra_coords[0][0] + tetra_coords[2][0];
    t03  = - tetra_coords[0][0] + tetra_coords[3][0];

    t11  = - tetra_coords[0][1] + tetra_coords[1][1];
    t12  = - tetra_coords[0][1] + tetra_coords[2][1];
    t13  = - tetra_coords[0][1] + tetra_coords[3][1];

    t21  = - tetra_coords[0][2] + tetra_coords[1][2];
    t22  = - tetra_coords[0][2] + tetra_coords[2][2];
    t23  = - tetra_coords[0][2] + tetra_coords[3][2];

    isop_0 = (  t00 * (t12*t23 - t13*t22)
              - t10 * (t02*t23 - t22*t03)
              + t20 * (t02*t13 - t12*t03)) / vol6;
    isop_1 = (- t00 * (t11*t23 - t13*t21)
              + t10 * (t01*t23 - t21*t03)
              - t20 * (t01*t13 - t03*t11)) / vol6;
    isop_2 = (  t00 * (t11*t22 - t21*t12)
              - t10 * (t01*t22 - t21*t02)
              + t20 * (t01*t12 - t11*t02)) / vol6;

    shapef[0] = 1. - isop_0 - isop_1 - isop_2;
    shapef[1] =      isop_0;
    shapef[2] =               isop_1;
    shapef[3] =                        isop_2;

    max_dist = -1.0;

    for (j = 0; j < 4; j++){

      dist = 2.*_ABS(shapef[j] - 0.5);

      if (max_dist < dist)
        max_dist = dist;
    }

    if (   (max_dist > -0.5 && max_dist < (1. + 2.*tolerance))
        && (max_dist < distance[i] || distance[i] < 0)) {
      location[i] = elt_num;
      distance[i] = max_dist;
    }

  }

}

/*---------------------------------------------------------------------------
 * Solve the equation "matrix.x = b" with Cramer's rule.
 *
 * parameters:
 *   m[3][3] <-- equation matrix
 *   b[3]    <-- b equation right hand side
 *   x[3]    <-> equation solution (unchanged if matrix is singular)
 *
 * returns:
 *   1 if matrix is singular, 0 otherwise
 *----------------------------------------------------------------------------*/

static int
_inverse_3x3(double  m[3][3],
             double  b[3],
             double  x[3])
{
  double det, det_inv, x0, x1, x2;

  det =   m[0][0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2])
        - m[1][0]*(m[0][1]*m[2][2] - m[2][1]*m[0][2])
        + m[2][0]*(m[0][1]*m[1][2] - m[1][1]*m[0][2]);

  if (_ABS(det) < _epsilon_denom)
    return 1;
  else
    det_inv = 1./det;

  /* Use local variables to ensure no aliasing */

  x0 = (  b[0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2])
        - b[1]*(m[0][1]*m[2][2] - m[2][1]*m[0][2])
        + b[2]*(m[0][1]*m[1][2] - m[1][1]*m[0][2])) * det_inv;

  x1 = (  m[0][0]*(b[1]*m[2][2] - b[2]*m[1][2])
        - m[1][0]*(b[0]*m[2][2] - b[2]*m[0][2])
        + m[2][0]*(b[0]*m[1][2] - b[1]*m[0][2])) * det_inv;

  x2 = (  m[0][0]*(m[1][1]*b[2] - m[2][1]*b[1])
        - m[1][0]*(m[0][1]*b[2] - m[2][1]*b[0])
        + m[2][0]*(m[0][1]*b[1] - m[1][1]*b[0])) * det_inv;

  /* Copy local variables to output */

  x[0] = x0; x[1] = x1; x[2] = x2;

  return 0;
}

/*----------------------------------------------------------------------------
 * Compute 3d shape functions and their derivatives given element
 * parametric coordinates.
 *
 * This function is adapted from the CGNS interpolation tool.
 *
 * parameters:
 *   elt_type    <-- type of element
 *   uvw[]       <-- parametric coordinates
 *   shapef[]    --> barycenter's coordinates
 *   deriv [][]  --> derivative of shape function
*----------------------------------------------------------------------------*/

static void
_compute_shapef_3d(my_ple_element_t  elt_type,
                   const double      uvw[3],
                   double            shapef[8],
                   double            deriv[8][3])

{
  switch (elt_type) {

  case MY_PLE_CELL_HEXA:

    shapef[0] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[1] = uvw[0] * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[2] = uvw[0] * uvw[1] * (1.0 - uvw[2]);
    shapef[3] = (1.0 - uvw[0]) * uvw[1] * (1.0 - uvw[2]);
    shapef[4] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * uvw[2];
    shapef[5] = uvw[0] * (1.0 - uvw[1]) * uvw[2];
    shapef[6] = uvw[0] * uvw[1] * uvw[2];
    shapef[7] = (1.0 - uvw[0]) * uvw[1] * uvw[2];

    if (deriv != NULL) {
      deriv[0][0] = -(1.0 - uvw[1]) * (1.0 - uvw[2]);
      deriv[0][1] = -(1.0 - uvw[0]) * (1.0 - uvw[2]);
      deriv[0][2] = -(1.0 - uvw[0]) * (1.0 - uvw[1]);
      deriv[1][0] =  (1.0 - uvw[1]) * (1.0 - uvw[2]);
      deriv[1][1] = -uvw[0] * (1.0 - uvw[2]);
      deriv[1][2] = -uvw[0] * (1.0 - uvw[1]);
      deriv[2][0] =  uvw[1] * (1.0 - uvw[2]);
      deriv[2][1] =  uvw[0] * (1.0 - uvw[2]);
      deriv[2][2] = -uvw[0] * uvw[1];
      deriv[3][0] = -uvw[1] * (1.0 - uvw[2]);
      deriv[3][1] =  (1.0 - uvw[0]) * (1.0 - uvw[2]);
      deriv[3][2] = -(1.0 - uvw[0]) * uvw[1];
      deriv[4][0] = -(1.0 - uvw[1]) * uvw[2];
      deriv[4][1] = -(1.0 - uvw[0]) * uvw[2];
      deriv[4][2] =  (1.0 - uvw[0]) * (1.0 - uvw[1]);
      deriv[5][0] =  (1.0 - uvw[1]) * uvw[2];
      deriv[5][1] = -uvw[0] * uvw[2];
      deriv[5][2] =  uvw[0] * (1.0 - uvw[1]);
      deriv[6][0] =  uvw[1] * uvw[2];
      deriv[6][1] =  uvw[0] * uvw[2];
      deriv[6][2] =  uvw[0] * uvw[1];
      deriv[7][0] = -uvw[1] * uvw[2];
      deriv[7][1] =  (1.0 - uvw[0]) * uvw[2];
      deriv[7][2] =  (1.0 - uvw[0]) * uvw[1];
    }

    break;

  case MY_PLE_CELL_PRISM:

    shapef[0] = (1.0 - uvw[0] - uvw[1]) * (1.0 - uvw[2]);
    shapef[1] = uvw[0] * (1.0 - uvw[2]);
    shapef[2] = uvw[1] * (1.0 - uvw[2]);
    shapef[3] = (1.0 - uvw[0] - uvw[1]) * uvw[2];
    shapef[4] = uvw[0] * uvw[2];
    shapef[5] = uvw[1] * uvw[2];

    if (deriv != NULL) {
      deriv[0][0] = -(1.0 - uvw[2]);
      deriv[0][1] = -(1.0 - uvw[2]);
      deriv[0][2] = -(1.0 - uvw[0] - uvw[1]);
      deriv[1][0] =  (1.0 - uvw[2]);
      deriv[1][1] =  0.0;
      deriv[1][2] = -uvw[0];
      deriv[2][0] =  0.0;
      deriv[2][1] =  (1.0 - uvw[2]);
      deriv[2][2] = -uvw[1];
      deriv[3][0] = -uvw[2];
      deriv[3][1] = -uvw[2];
      deriv[3][2] =  (1.0 - uvw[0] - uvw[1]);
      deriv[4][0] =  uvw[2];
      deriv[4][1] =  0.0;
      deriv[4][2] =  uvw[0];
      deriv[5][0] =  0.0;
      deriv[5][1] =  uvw[2];
      deriv[5][2] =  uvw[1];
    }

    break;

  case MY_PLE_CELL_PYRAM:

    shapef[0] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[1] = uvw[0] * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[2] = uvw[0] * uvw[1] * (1.0 - uvw[2]);
    shapef[3] = (1.0 - uvw[0]) * uvw[1] * (1.0 - uvw[2]);
    shapef[4] = uvw[2];

    if (deriv != NULL) {
      deriv[0][0] = -(1.0 - uvw[1]) * (1.0 - uvw[2]);
      deriv[0][1] = -(1.0 - uvw[0]) * (1.0 - uvw[2]);
      deriv[0][2] = -(1.0 - uvw[0]) * (1.0 - uvw[1]);
      deriv[1][0] =  (1.0 - uvw[1]) * (1.0 - uvw[2]);
      deriv[1][1] = -uvw[0] * (1.0 - uvw[2]);
      deriv[1][2] = -uvw[0] * (1.0 - uvw[1]);
      deriv[2][0] =  uvw[1] * (1.0 - uvw[2]);
      deriv[2][1] =  uvw[0] * (1.0 - uvw[2]);
      deriv[2][2] = -uvw[0] * uvw[1];
      deriv[3][0] = -uvw[1] * (1.0 - uvw[2]);
      deriv[3][1] =  (1.0 - uvw[0]) * (1.0 - uvw[2]);
      deriv[3][2] = -(1.0 - uvw[0]) * uvw[1];
      deriv[4][0] =  0.0;
      deriv[4][1] =  0.0;
      deriv[4][2] =  1.0;
    }

    break;

  default:
    ple_error(__FILE__, __LINE__, 0,
              _("_compute_shapef: unhandled element type %s\n"),
              my_ple_element_type_name[elt_type]);

  }

}

/*----------------------------------------------------------------------------
 * Compute hexahedron, pyramid, or prism parametric coordinates for a
 * given point.
 *
 * This function is adapted from the CGNS interpolation tool.
 *
 * parameters:
 *   elt_type            <-- type of element
 *   point_coords        <-- point coordinates
 *   vertex_coords[]     <-- pointer to element vertex coordinates
 *   tolerance           <-- location tolerance factor
 *   uvw[]               --> parametric coordinates of point in element
*----------------------------------------------------------------------------*/
static int
_compute_uvw(my_ple_element_t    elt_type,
             const ple_coord_t   point_coords[],
             double              vertex_coords[8][3],
             double              tolerance,
             double              uvw[3])

{
  int i, j, n_elt_vertices, iter;
  int max_iter = 20;
  double dist;
  double a[3][3], b[3], x[3], shapef[8], dw[8][3];

  n_elt_vertices = my_ple_mesh_n_vertices_element[elt_type];

  assert(   elt_type == MY_PLE_CELL_HEXA
         || elt_type == MY_PLE_CELL_PRISM
         || elt_type == MY_PLE_CELL_PYRAM);

  /* Use Newton-method to determine parametric coordinates and shape function */

  for (i = 0; i < 3; i++)
    uvw[i] = 0.5;

  for (iter = 0; iter < max_iter; iter++) {

    _compute_shapef_3d(elt_type, uvw, shapef, dw);

    b[0] = - point_coords[0];
    b[1] = - point_coords[1];
    b[2] = - point_coords[2];

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++)
        a[i][j] = 0.0;
    }

    for (i = 0; i < n_elt_vertices; i++) {

      b[0] += (shapef[i] * vertex_coords[i][0]);
      b[1] += (shapef[i] * vertex_coords[i][1]);
      b[2] += (shapef[i] * vertex_coords[i][2]);

      for (j = 0; j < 3; j++) {
        a[0][j] -= (dw[i][j] * vertex_coords[i][0]);
        a[1][j] -= (dw[i][j] * vertex_coords[i][1]);
        a[2][j] -= (dw[i][j] * vertex_coords[i][2]);
      }

    }

    if (_inverse_3x3(a, b, x))
      return 0;

    dist = 0.0;

    for (i = 0; i < 3; i++) {
      dist += x[i] * x[i];
      uvw[i] += x[i];
    }

    if (dist <= (tolerance * tolerance))
      return 1;

  }

  return 0;
}

/*----------------------------------------------------------------------------
 * Locate points in a given 3d cell (other than tetrahedra or polyhedra,
 * handlesd elsewhere), updating the location[] and distance[] arrays
 * associated with a set of points.
 *
 * parameters:
 *   elt_num             <-- element number
 *   elt_type            <-- type of element
 *   element_vertex_num  <-- element vertex numbers
 *   vertex_coords[]     <-- pointer to vertex coordinates
 *   point_coords        <-- point coordinates
 *   n_points_in_extents <-- number of points in element extents
 *   points_in_extents   <-- ids of points in extents
 *   tolerance           <-- associated tolerance
 *   location            <-> number of element containing or closest to each
 *                           point (size: n_points)
 *   distance            <-> distance from point to element indicated by
 *                           location[]: < 0 if unlocated, 0 - 1 if inside,
 *                           > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_in_cell_3d(ple_lnum_t          elt_num,
                   my_ple_element_t    elt_type,
                   const ple_lnum_t    element_vertex_num[],
                   const ple_coord_t   vertex_coords[],
                   const ple_coord_t   point_coords[],
                   ple_lnum_t          n_points_in_extents,
                   const ple_lnum_t    points_in_extents[],
                   double              tolerance,
                   ple_lnum_t          location[],
                   float               distance[])
{
  int i, j, k, n_vertices;
  ple_lnum_t coord_idx, vertex_id;

  double uvw[3], dist, shapef[8],max_dist;
  double  _vertex_coords[8][3];

  n_vertices = my_ple_mesh_n_vertices_element[elt_type];

  /* Initialize local element coordinates copy */

  for (vertex_id = 0; vertex_id < n_vertices; vertex_id++) {

    coord_idx = element_vertex_num[vertex_id] -1;

    for (j = 0; j < 3; j++)
      _vertex_coords[vertex_id][j] = vertex_coords[(coord_idx * 3) + j];

  }

  /* Shape functions may be computed directly with tetrahedra */

  if (elt_type == MY_PLE_CELL_TETRA)

    _locate_in_tetra(elt_num,
                     _vertex_coords,
                     point_coords,
                     n_points_in_extents,
                     points_in_extents,
                     tolerance,
                     location,
                     distance);

  /* For cell shapes other than tetrahedra, find shape functions iteratively */

  else {

    for (k = 0; k < n_points_in_extents; k++) {

      i = points_in_extents[k];

      if (_compute_uvw(elt_type,
                       point_coords + 3*i,
                       _vertex_coords,
                       tolerance,
                       uvw)) {

        max_dist = -1.0;

        /* For hexahedra, no need to compute shape functions, as
           the 3 parametric coordinates are simpler to use */

        if (elt_type == MY_PLE_CELL_HEXA) {

          for (j = 0; j < 3; j++){

            dist = 2.*_ABS(uvw[j] - 0.5);

            if (max_dist < dist)
              max_dist = dist;
          }

        }

        /* For pyramids ands prisms, we need to compute shape functions */

        else {

          _compute_shapef_3d(elt_type, uvw, shapef, NULL);

          for (j = 0; j < n_vertices; j++){

            dist = 2.*_ABS(shapef[j] - 0.5);

          if (max_dist < dist)
            max_dist = dist;
          }

        }

        /* For all element types, update location and distance arrays */

        if (   (max_dist > -0.5 && max_dist < (1. + 2.*tolerance))
            && (max_dist < distance[i] || distance[i] < 0)) {
          location[i] = elt_num;
          distance[i] = max_dist;
        }

      }

    } /* End of loop on points in extents */

  }

}

/*----------------------------------------------------------------------------
 * Find elements in a given polyhedral section containing points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this section, or closer to one
 * than to previously encountered elements.
 *
 * Space and element dimensions are both equal to 3 here.
 *
 * parameters:
 *   this_section      <-- pointer to mesh section representation structure
 *   vertex_coords     <-- pointer to vertex coordinates
 *   tolerance         <-- associated tolerance
 *   base_element_num  <-- < 0 for location relative to parent element numbers,
 *                         number of elements in preceding sections of same
 *                         element dimension + 1 otherwise
 *   point_coords      <-- point coordinates
 *   octree            <-- point octree
 *   points_in_extents <-> array for query of ids of points in extents
 *                         (size: n_points, less usually needed)
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated, 0 - 1 if inside,
 *                         > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_polyhedra_section_locate(const my_ple_mesh_section_t  *this_section,
                          const ple_coord_t             vertex_coords[],
                          double                        tolerance,
                          ple_lnum_t                    base_element_num,
                          const ple_coord_t             point_coords[],
                          _octree_t                    *octree,
                          ple_lnum_t                    points_in_extents[],
                          ple_lnum_t                    location[],
                          float                         distance[])
{
  ple_lnum_t  i, j, k, n_vertices, face_id, vertex_id, elt_num;
  ple_coord_t  center[3];
  double elt_extents[6];

  double _tolerance = tolerance * 2; /* double tolerance, as polyhedra is
                                        split into tetrahedra, whose extents
                                        are 1/2 the polyhedron extents */

  ple_lnum_t n_vertices_max = 0;
  ple_lnum_t n_points_in_extents = 0;
  ple_lnum_t *triangle_vertices = NULL;
  _triangulate_state_t *state = NULL;

  assert(this_section->face_index != NULL);

  /* Return immediately if nothing to do for this rank */

  if (this_section->n_elements == 0)
    return;

  /* Counting loop on faces */

  for (i = 0; i < this_section->n_faces; i++) {
    n_vertices =   this_section->vertex_index[i + 1]
                 - this_section->vertex_index[i];
    if (n_vertices > n_vertices_max)
      n_vertices_max = n_vertices;
  }

  if (n_vertices_max < 3)
    return;

  PLE_MALLOC(triangle_vertices, (n_vertices_max-2)*3, int);
  state = _triangulate_state_create(n_vertices_max);

  /* Loop on elements */

  for (i = 0; i < this_section->n_elements; i++) {

    _Bool elt_initialized = false;

    /* Compute extents */

    for (j = this_section->face_index[i];
         j < this_section->face_index[i + 1];
         j++) {
      face_id = _ABS(this_section->face_num[j]) - 1;
      for (k = this_section->vertex_index[face_id];
           k < this_section->vertex_index[face_id + 1];
           k++) {
        vertex_id = this_section->vertex_num[k] - 1;

        _update_elt_extents(3,
                            vertex_id,
                            vertex_coords,
                            elt_extents,
                            &elt_initialized);

      }
    }

    _elt_extents_finalize(3, 3, tolerance, elt_extents);

    if (base_element_num < 0) {
      if (this_section->element_num != NULL)
        elt_num = this_section->element_num[i];
      else
        elt_num = i + 1;
    }
    else
      elt_num = base_element_num + i;

    _query_octree(elt_extents,
                  point_coords,
                  octree,
                  &n_points_in_extents,
                  points_in_extents);

    if (n_points_in_extents < 1)
      continue;

    /* Compute psuedo-element center */

    for (j = 0; j < 3; j++)
      center[j] = (elt_extents[j] + elt_extents[j + 3]) * 0.5;

    /* Loop on element faces */

    for (j = this_section->face_index[i];
         j < this_section->face_index[i + 1];
         j++) {

      ple_lnum_t n_triangles;

      const ple_lnum_t *_vertex_num;

      face_id = _ABS(this_section->face_num[j]) - 1;

      n_vertices = (  this_section->vertex_index[face_id + 1]
                    - this_section->vertex_index[face_id]);

      _vertex_num = (  this_section->vertex_num
                     + this_section->vertex_index[face_id]);

      if (n_vertices == 4)
        n_triangles = _triangulate_quadrangle(3,
                                              vertex_coords,
                                              _vertex_num,
                                              triangle_vertices);

      else if (n_vertices > 4)
        n_triangles = _triangulate_polygon(3,
                                           n_vertices,
                                           vertex_coords,
                                           _vertex_num,
                                           triangle_vertices,
                                           state);

      else { /* n_vertices == 3 */

        n_triangles = 1;
        for (k = 0; k < 3; k++)
          triangle_vertices[k] = _vertex_num[k];

      }

      /* Loop on face triangles so as to loop on tetrahedra
         built by joining face triangles and psuedo-center */

      for (k = 0; k < n_triangles; k++) {

        ple_lnum_t l, coord_id[3];
        ple_coord_t tetra_coords[4][3];

        coord_id[0] = triangle_vertices[k*3    ] - 1;
        coord_id[1] = triangle_vertices[k*3 + 2] - 1;
        coord_id[2] = triangle_vertices[k*3 + 1] - 1;

        for (l = 0; l < 3; l++) {
          tetra_coords[0][l] = vertex_coords[3*coord_id[0] + l];
          tetra_coords[1][l] = vertex_coords[3*coord_id[1] + l];
          tetra_coords[2][l] = vertex_coords[3*coord_id[2] + l];
          tetra_coords[3][l] = center[l];
        }

        _locate_in_tetra(elt_num,
                         tetra_coords,
                         point_coords,
                         n_points_in_extents,
                         points_in_extents,
                         _tolerance,
                         location,
                         distance);

      } /* End of loop on face triangles */

    } /* End of loop on element faces */

    _locate_in_extents(elt_num,
                       3,
                       elt_extents,
                       point_coords,
                       n_points_in_extents,
                       points_in_extents,
                       location,
                       distance);

  } /* End of loop on elements */

  PLE_FREE(triangle_vertices);
  _triangulate_state_destroy(&state);
}

/*----------------------------------------------------------------------------
 * Find elements in a given polygonal section containing 3d points: updates
 * the location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this section, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   this_section      <-- pointer to mesh section representation structure
 *   vertex_coords     <-- pointer to vertex coordinates
 *   tolerance         <-- associated tolerance
 *   base_element_num  <-- < 0 for location relative to parent element numbers,
 *                         number of elements in preceding sections of same
 *                         element dimension + 1 otherwise
 *   n_points          <-- number of points to locate
 *   point_coords      <-- point coordinates
 *   octree            <-- point octree
 *   points_in_extents <-- array for query of ids of points in extents
 *                         (size: n_points, less usually needed)
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated, absolute distance
 *                         to element if located (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_polygons_section_locate_3d(const my_ple_mesh_section_t   *this_section,
                            const ple_coord_t              vertex_coords[],
                            const double                   tolerance,
                            ple_lnum_t                     base_element_num,
                            const ple_coord_t              point_coords[],
                            _octree_t                     *octree,
                            ple_lnum_t                     points_in_extents[],
                            ple_lnum_t                     location[],
                            float                          distance[])
{
  ple_lnum_t  i, j, n_vertices, vertex_id, elt_num;
  int n_triangles;
  double elt_extents[6];

  int n_vertices_max = 0;
  ple_lnum_t n_points_in_extents = 0;

  ple_lnum_t *triangle_vertices = NULL;
  _triangulate_state_t *state = NULL;

  /* Return immediately if nothing to do for this rank */

  if (this_section->n_elements == 0)
    return;

  /* Counting loop on elements */

  for (i = 0; i < this_section->n_elements; i++) {

    n_vertices = (  this_section->vertex_index[i + 1]
                  - this_section->vertex_index[i]);

    if (n_vertices > n_vertices_max)
      n_vertices_max = n_vertices;

  }

  if (n_vertices_max < 3)
    return;

  PLE_MALLOC(triangle_vertices, (n_vertices_max-2)*3, int);
  state = _triangulate_state_create(n_vertices_max);

  /* Main loop on elements */

  for (i = 0; i < this_section->n_elements; i++) {

    _Bool elt_initialized = false;

    for (j = this_section->vertex_index[i];
         j < this_section->vertex_index[i + 1];
         j++) {
      vertex_id = this_section->vertex_num[j] - 1;

      _update_elt_extents(3,
                          vertex_id,
                          vertex_coords,
                          elt_extents,
                          &elt_initialized);

    }

    _elt_extents_finalize(3, 2, tolerance, elt_extents);

    if (base_element_num < 0) {
      if (this_section->element_num != NULL)
        elt_num = this_section->element_num[i];
      else
        elt_num = i + 1;
    }
    else
      elt_num = base_element_num + i;

    _query_octree(elt_extents,
                  point_coords,
                  octree,
                  &n_points_in_extents,
                  points_in_extents);

    /* Triangulate polygon */

    n_vertices = (  this_section->vertex_index[i + 1]
                  - this_section->vertex_index[i]);
    vertex_id = this_section->vertex_index[i];

    n_triangles = _triangulate_polygon(3,
                                       n_vertices,
                                       vertex_coords,
                                       (  this_section->vertex_num
                                        + vertex_id),
                                       triangle_vertices,
                                       state);

    /* Locate on triangulated polygon */

    _locate_on_triangles_3d(elt_num,
                            n_triangles,
                            triangle_vertices,
                            vertex_coords,
                            point_coords,
                            n_points_in_extents,
                            points_in_extents,
                            tolerance,
                            location,
                            distance);

  } /* End of loop on elements */

  PLE_FREE(triangle_vertices);
  _triangulate_state_destroy(&state);
}

/*----------------------------------------------------------------------------
 * Find elements in a given section containing 3d points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this section, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   this_section      <-- pointer to mesh section representation structure
 *   vertex_coords     <-- pointer to vertex coordinates
 *   tolerance         <-- associated tolerance
 *   base_element_num  <-- < 0 for location relative to parent element numbers,
 *                         number of elements in preceding sections of same
 *                         element dimension + 1 otherwise
 *   point_coords      <-- point coordinates
 *   octree            <-- point octree
 *   points_in_extents <-- array for query of ids of points in extents
 *                         (size: octree->n_points, less usually needed)
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated; 0 - 1 if inside,
 *                         and > 1 if outside a volume element, or absolute
 *                         distance to a surface element (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_mesh_section_locate_3d(const my_ple_mesh_section_t  *this_section,
                        const ple_coord_t             vertex_coords[],
                        double                        tolerance,
                        ple_lnum_t                    base_element_num,
                        const ple_coord_t             point_coords[],
                        _octree_t                    *octree,
                        ple_lnum_t                    points_in_extents[],
                        ple_lnum_t                    location[],
                        float                         distance[])
{
  ple_lnum_t  i, j, vertex_id, elt_num, triangle_vertices[6];
  int n_triangles;
  double elt_extents[6];

  ple_lnum_t n_points_in_extents = 0;

  /* If section contains polyhedra */

  if (this_section->type == MY_PLE_CELL_POLY)

    _polyhedra_section_locate(this_section,
                              vertex_coords,
                              tolerance,
                              base_element_num,
                              point_coords,
                              octree,
                              points_in_extents,
                              location,
                              distance);

  /* If section contains polygons */

  else if (this_section->type == MY_PLE_FACE_POLY)

    _polygons_section_locate_3d(this_section,
                                vertex_coords,
                                tolerance,
                                base_element_num,
                                point_coords,
                                octree,
                                points_in_extents,
                                location,
                                distance);

  /* If section contains regular elements */

  else {

    for (i = 0; i < this_section->n_elements; i++) {

      _Bool elt_initialized = false;

      if (base_element_num < 0) {
        if (this_section->element_num != NULL)
          elt_num = this_section->element_num[i];
        else
          elt_num = i + 1;
      }
      else
        elt_num = base_element_num + i;

      for (j = 0; j < this_section->stride; j++) {

        vertex_id = this_section->vertex_num[i*this_section->stride + j] - 1;

        _update_elt_extents(3,
                            vertex_id,
                            vertex_coords,
                            elt_extents,
                            &elt_initialized);

      }

      _elt_extents_finalize(3,
                            this_section->entity_dim,
                            tolerance,
                            elt_extents);

      _query_octree(elt_extents,
                    point_coords,
                    octree,
                    &n_points_in_extents,
                    points_in_extents);

      if (this_section->entity_dim == 3)

        _locate_in_cell_3d(elt_num,
                           this_section->type,
                           this_section->vertex_num + i*this_section->stride,
                           vertex_coords,
                           point_coords,
                           n_points_in_extents,
                           points_in_extents,
                           tolerance,
                           location,
                           distance);

      else if (this_section->entity_dim == 2) {

        if (this_section->type == MY_PLE_FACE_QUAD)

          n_triangles = _triangulate_quadrangle(3,
                                                vertex_coords,
                                                (  this_section->vertex_num
                                                 + i*this_section->stride),
                                                triangle_vertices);

        else {

          assert(this_section->type == MY_PLE_FACE_TRIA);

          n_triangles = 1;
          for (j = 0; j < 3; j++)
            triangle_vertices[j]
              = this_section->vertex_num[i*this_section->stride + j];


        }

        _locate_on_triangles_3d(elt_num,
                                n_triangles,
                                triangle_vertices,
                                vertex_coords,
                                point_coords,
                                n_points_in_extents,
                                points_in_extents,
                                tolerance,
                                location,
                                distance);
      }

      else if (this_section->entity_dim == 1) {

        assert(this_section->type == MY_PLE_EDGE);

        _locate_on_edge_3d(elt_num,
                           this_section->vertex_num + i*this_section->stride,
                           vertex_coords,
                           point_coords,
                           n_points_in_extents,
                           points_in_extents,
                           tolerance,
                           location,
                           distance);

      }
    }

  }
}

/*----------------------------------------------------------------------------
 * Find elements in a given section containing 2d points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this section, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   this_section      <-- pointer to mesh section representation structure
 *   vertex_coords     <-- pointer to vertex coordinates
 *   tolerance         <-- associated tolerance
 *   base_element_num  <-- < 0 for location relative to parent element numbers,
 *                         number of elements in preceding sections of same
 *                         element dimension + 1 otherwise
 *   point_coords      <-- point coordinates
 *   quadtree          <-- point quadtree
 *   points_in_extents <-- array for query of ids of points in extents
 *                         (size: quadtree->n_points, less usually needed)
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated, 0 - 1 if inside,
 *                         > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_mesh_section_locate_2d(const my_ple_mesh_section_t  *this_section,
                        const ple_coord_t             vertex_coords[],
                        double                        tolerance,
                        ple_lnum_t                    base_element_num,
                        const ple_coord_t             point_coords[],
                        _quadtree_t                  *quadtree,
                        ple_lnum_t                    points_in_extents[],
                        ple_lnum_t                    location[],
                        float                         distance[])
{
  ple_lnum_t  i, j, vertex_id, elt_num, _triangle_vertices[6];
  int n_vertices;
  double elt_extents[4];

  int n_triangles = 0;
  int n_vertices_max = 0;
  ple_lnum_t n_points_in_extents = 0;
  ple_lnum_t *triangle_vertices = _triangle_vertices;
  _triangulate_state_t *state = NULL;

  /* Return immediately if nothing to do for this rank */

  if (this_section->n_elements == 0)
    return;

  /* Count maximum number of vertices */

  if (this_section->type == MY_PLE_FACE_POLY) {

    for (i = 0; i < this_section->n_elements; i++) {

      n_vertices = (  this_section->vertex_index[i + 1]
                    - this_section->vertex_index[i]);

      if (n_vertices > n_vertices_max)
        n_vertices_max = n_vertices;

    }

    if (n_vertices_max < 3)
      return;

    PLE_MALLOC(triangle_vertices, (n_vertices_max-2)*3, int);
    state = _triangulate_state_create(n_vertices_max);

  }

  else if (this_section->type == MY_PLE_FACE_QUAD)
    n_vertices_max = 4;

  else if (this_section->type == MY_PLE_FACE_TRIA)
    n_vertices_max = 3;

  /* Main loop on elements */

  for (i = 0; i < this_section->n_elements; i++) {

    _Bool elt_initialized = false;

    if (this_section->type == MY_PLE_FACE_POLY) {

      for (j = this_section->vertex_index[i];
           j < this_section->vertex_index[i + 1];
           j++) {
        vertex_id = this_section->vertex_num[j] - 1;

        _update_elt_extents(2,
                            vertex_id,
                            vertex_coords,
                            elt_extents,
                            &elt_initialized);

      }

    }
    else {

      for (j = 0; j < this_section->stride; j++) {

        vertex_id = this_section->vertex_num[i*this_section->stride + j] - 1;

        _update_elt_extents(2,
                            vertex_id,
                            vertex_coords,
                            elt_extents,
                            &elt_initialized);

      }

    }

    _elt_extents_finalize(2,
                          this_section->entity_dim,
                          tolerance,
                          elt_extents);

    if (base_element_num < 0) {
      if (this_section->element_num != NULL)
        elt_num = this_section->element_num[i];
      else
        elt_num = i + 1;
    }
    else
      elt_num = base_element_num + i;

    _query_quadtree(elt_extents,
                    point_coords,
                    quadtree,
                    &n_points_in_extents,
                    points_in_extents);

    /* Divide all face types into triangles */

    if (this_section->type == MY_PLE_FACE_POLY) {

      /* Triangulate polygon */

      n_vertices = (  this_section->vertex_index[i + 1]
                    - this_section->vertex_index[i]);
      vertex_id = this_section->vertex_index[i];

      n_triangles = _triangulate_polygon(2,
                                         n_vertices,
                                         vertex_coords,
                                         (  this_section->vertex_num
                                          + vertex_id),
                                         triangle_vertices,
                                         state);

    }
    else if (this_section->type == MY_PLE_FACE_QUAD) {

      /* Triangulate quadrangle */

      n_triangles = _triangulate_quadrangle(2,
                                            vertex_coords,
                                            (  this_section->vertex_num
                                             + i*this_section->stride),
                                            triangle_vertices);

    }

    else if (this_section->type == MY_PLE_FACE_TRIA) {

      /* Already a triangle */

      n_triangles = 1;

      for (j = 0; j < 3; j++)
        triangle_vertices[j]
          = this_section->vertex_num[i*this_section->stride + j];

    }

    /* Locate on triangulated face */

    if (this_section->entity_dim == 2)

      _locate_on_triangles_2d(elt_num,
                              n_triangles,
                              triangle_vertices,
                              vertex_coords,
                              point_coords,
                              n_points_in_extents,
                              points_in_extents,
                              tolerance,
                              location,
                              distance);

    else if (this_section->entity_dim == 1) {

      assert(this_section->type == MY_PLE_EDGE);

      _locate_on_edge_2d(elt_num,
                         this_section->vertex_num + i*this_section->stride,
                         vertex_coords,
                         point_coords,
                         n_points_in_extents,
                         points_in_extents,
                         tolerance,
                         location,
                         distance);

    }

  } /* End of loop on elements */

  /* Free axiliary arrays and structures */

  if (triangle_vertices != _triangle_vertices)
    PLE_FREE(triangle_vertices);

  if (state != NULL)
    _triangulate_state_destroy(&state);
}

/*----------------------------------------------------------------------------
 * Find elements in a given section containing 1d points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this section, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   this_section      <-- pointer to mesh section representation structure
 *   vertex_coords     <-- pointer to vertex coordinates
 *   tolerance         <-- associated tolerance
 *   base_element_num  <-- < 0 for location relative to parent element numbers,
 *                         number of elements in preceding sections of same
 *                         element dimension + 1 otherwise
 *   n_points          <-- number of points to locate
 *   point_coords      <-- point coordinates
 *   points_in_extents <-- array for query of ids of points in extents
 *                         (size: n_points, less usually needed)
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated, 0 - 1 if inside,
 *                         > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_mesh_section_locate_1d(const my_ple_mesh_section_t  *this_section,
                        const ple_coord_t             vertex_coords[],
                        double                        tolerance,
                        ple_lnum_t                    base_element_num,
                        ple_lnum_t                    n_points,
                        const ple_coord_t             point_coords[],
                        ple_lnum_t                    location[],
                        float                         distance[])
{
  ple_lnum_t  i, j, vertex_id, elt_num;
  ple_coord_t edge_coords[2];
  double delta, elt_extents[2];

  for (i = 0; i < this_section->n_elements; i++) {

    if (base_element_num < 0) {
      if (this_section->element_num != NULL)
        elt_num = this_section->element_num[i];
      else
        elt_num = i + 1;
    }
    else
      elt_num = base_element_num + i;

    for (j = 0; j < 2; j++) {

      vertex_id = this_section->vertex_num[i*this_section->stride] - 1;

      edge_coords[j] = vertex_coords[vertex_id];

    }

    if (edge_coords[0] < edge_coords[1]) {
      elt_extents[0] = edge_coords[0];
      elt_extents[1] = edge_coords[1];
    }
    else {
      elt_extents[0] = edge_coords[1];
      elt_extents[1] = edge_coords[0];
    }

    delta = (elt_extents[1] - elt_extents[0]) * tolerance;

    elt_extents[0] -= delta;
    elt_extents[1] += delta;

    _locate_by_extents_1d(elt_num,
                          elt_extents,
                          n_points,
                          point_coords,
                          location,
                          distance);

  }

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Query number of extents and compute extents of a mesh representation.
 *
 * The minimum required functionality for this function is to compute
 * whole mesh extents, but it could also return extents of individual
 * elements, or intermediate extents of mesh subdivisions or coarsened
 * elements. As such, it takes an argument indicating the maximum
 * local number of extents it should compute (based on the size of
 * the extents array argument), but returns the number of extents
 * really computed, which may be lower (usually 1 for mesh extents,
 * possibly even 0 if the local mesh is empty). If n_max_extents = 1,
 * the whole mesh extents should be computed.
 *
 * If n_max_extents is set to a negative value (-1), no extents are computed,
 * but the function returns the maximum number of extents it may compute.
 * This query mode allows for the caller to allocate the correct amount
 * of memory for a subsequent call.
 *
 * parameters:
 *   mesh          <-- pointer to mesh representation structure
 *   n_max_extents <-- maximum number of sub-extents (such as element extents)
 *                     to compute, or -1 to query
 *   tolerance     <-- addition to local extents of each element:
 *                     extent = base_extent * (1 + tolerance)
 *   extents       <-> extents associated with the mesh or elements (or even
 *                     aggregated elements in case of coarser representation):
 *                     x_min_0, y_min_0, ..., x_max_i, y_max_i, ...
 *                     (size: 2*dim*n_max_extents), ignored in query mode
 * returns:
 *   the number of extents computed
 *----------------------------------------------------------------------------*/

ple_lnum_t
my_ple_mesh_extents(const void  *mesh,
                    ple_lnum_t   n_max_extents,
                    double       tolerance,
                    double       extents[])
{
  int i, j;
  int dim;
  double section_extents[6];

  const my_ple_mesh_t  *m = mesh;

  if (mesh == NULL)
    return 0;

  /* In query mode, return 1 (currently, global extents only) */

  if (n_max_extents < 0)
    return 1;

  dim = m->dim;

  /* initialize extents in case mesh is empty or dim < 3 */
  for (i = 0; i < dim; i++) {
    extents[i]       =  HUGE_VAL;
    extents[i + dim] = -HUGE_VAL;
  }

  /* Compute extents */

  for (i = 0; i < m->n_sections; i++) {

    _mesh_section_extents(m->sections[i],
                          m->dim,
                          m->vertex_coords,
                          tolerance,
                          section_extents);

    for (j = 0; j < m->dim; j++) {
      if (section_extents[j] < extents[j])
        extents[j] = section_extents[j];
      if (section_extents[j+dim] > extents[j+dim])
        extents[j+dim] = section_extents[j+dim];
    }

  }

  return 1; /* 1 global mesh extent */
}

/*----------------------------------------------------------------------------
 * Find elements in a given nodal mesh containing points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this mesh, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   mesh              <-- pointer to mesh representation structure
 *   tolerance         <-- associated tolerance
 *   locate_on_parents <-- location relative to parent element numbers if
 *                         true, id of element + 1 in concatenated sections
 *                         of same element dimension if false
 *   n_points          <-- number of points to locate
 *   point_coords      <-- point coordinates
 *   point_tag         <-- optional point tag (size: n_points, ignored here)
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated; 0 - 1 if inside,
 *                         and > 1 if outside a volume element, or absolute
 *                         distance to a surface element (size: n_points)
 *----------------------------------------------------------------------------*/

void
my_ple_point_location_contain(const void          *mesh,
                              double               tolerance,
                              _Bool                locate_on_parents,
                              ple_lnum_t           n_points,
                              const ple_coord_t    point_coords[],
                              const int            point_tag[],
                              ple_lnum_t           location[],
                              float                distance[])
{
  int i;
  int max_entity_dim;
  ple_lnum_t   base_element_num;
  ple_lnum_t  *points_in_extents = NULL;
  const my_ple_mesh_t  *m = mesh;

  if (mesh == NULL)
    return;

  if (locate_on_parents == true)
    base_element_num = -1;
  else
    base_element_num = 1;

  max_entity_dim = my_ple_mesh_get_max_entity_dim(m);

  /* Build point query list
     (max size: n_points, usually much less) */

  PLE_MALLOC(points_in_extents, n_points, ple_lnum_t);

  /* Use octree for 3d point location */

  if (m->dim == 3) {

    _octree_t  octree = _build_octree(n_points, point_coords);

    /* Locate for all sections */

    for (i = 0; i < m->n_sections; i++) {

      const my_ple_mesh_section_t  *this_section = m->sections[i];

      if (this_section->entity_dim == max_entity_dim) {

        _mesh_section_locate_3d(this_section,
                                m->vertex_coords,
                                tolerance,
                                base_element_num,
                                point_coords,
                                &octree,
                                points_in_extents,
                                location,
                                distance);

        if (base_element_num > -1)
          base_element_num += this_section->n_elements;

      }

    }

    _free_octree(&octree);
  }

  /* Use quadtree for 2d point location */

  else if (m->dim == 2) {

    _quadtree_t  quadtree = _build_quadtree(n_points, point_coords);

    /* Locate for all sections */

    for (i = 0; i < m->n_sections; i++) {

      const my_ple_mesh_section_t  *this_section = m->sections[i];

      if (this_section->entity_dim == max_entity_dim) {

        _mesh_section_locate_2d(this_section,
                                m->vertex_coords,
                                tolerance,
                                base_element_num,
                                point_coords,
                                &quadtree,
                                points_in_extents,
                                location,
                                distance);

        if (base_element_num > -1)
          base_element_num += this_section->n_elements;

      }

    }

    _free_quadtree(&quadtree);
  }

  /* Use brute force for 1d point location */

  else if (m->dim == 1) {

    for (i = 0; i < m->n_sections; i++) {

      const my_ple_mesh_section_t  *this_section = m->sections[i];

      if (this_section->entity_dim == 1) {

        _mesh_section_locate_1d(this_section,
                                m->vertex_coords,
                                tolerance,
                                base_element_num,
                                n_points,
                                point_coords,
                                location,
                                distance);

        if (base_element_num > -1)
          base_element_num += this_section->n_elements;

      }
    }
  }

  PLE_FREE(points_in_extents);
}

/*----------------------------------------------------------------------------*/

#undef _DOT_PRODUCT_3D
#undef _DOT_PRODUCT_2D
#undef _MODULE_3D
#undef _CROSS_PRODUCT_3D
#undef _CROSS_PRODUCT_2D

#ifdef __cplusplus
}
#endif /* __cplusplus */
