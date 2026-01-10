/*============================================================================
 * Mesh refinement.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include <float.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

#include <new>  // For placement new

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "fvm/fvm_io_num.h"
#include "fvm/fvm_triangulate.h"

#include "base/cs_assert.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "base/cs_order.h"
#include "base/cs_parall.h"

#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_adjacencies.h"
#include "mesh/cs_mesh_quantities.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "mesh/cs_mesh_refine.h"

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_mesh_refine.cpp
        Mesh refinement.

  Cells are refined based on predefined templates for given (autodetected)
  cell types. Refinement may optionally be propagated to neighbor cells.

  The basic face refinement templates are given here:

  \image html r_templates_2d_base.png "Basic face refinement templates"

  Polygonal faces do not use a specific template, but their edges
  may be subdivided.

  Tetrahedra are divided into 8 tetrahedra, as follows:
  \image html r_template_tetra.png "Basic tetrahedron refinement template and subcells"

  Pyramids are divided into 6 pyramids and 4 tetrahedra, as follows
  (with the pyramids a darker shade and tetrahedra lighter shade
  on the right view):
  \image html r_template_pyram.png "Basic pyramid refinement template and subcells"

  Prisms are divided into 8 prisms, as follows:
  \image html r_template_prism.png "Basic prism refinement template and subcells"

  Hexahedra are divided into 8 hexahedra, as follows:
  \image html r_template_hexa.png "Basic hexahedron refinement template and subcells"

  Polyhedra are subdivided by adding a vertex to the cell's center,
  and subdividing the cells into n-faced cones joing the cell faces
  and cell center (tetrahedra for triangle faces, pyramids for
  quadrangle faces, generic cone-shaped polyhedra otherwise).
  The following view provides and example for one hexahedral cell
  with one face subdivided (as may occur for hexahedral cells adjacent
  to refined hexahadra when a conforming refinement is required):
  \image html r_template_polyhedron_e1.png "Polyhedron refinement and subcells example"
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* Refinement template used (should be extensible and allow checking) */

typedef enum {

  CS_REFINE_NONE,        /*!< no refinement */

  CS_REFINE_DEFAULT,     /*!< default refinement */

  CS_REFINE_TRIA,        /*!< simple triangle subdivision into 4 triangles
                           using mid-edge vertices, no added face vertex */
  CS_REFINE_TRIA_Q,      /*!< simple triangle subdivision into 3 quadrangles
                           using mid-edge vertices and a face center */
  CS_REFINE_QUAD,        /*!< simple quadrangle subdivision into 4 quadrangles
                           using mid-edge vertices and a face center */

  CS_REFINE_POLYGON,     /*!< refine polygon by tesselation, no vertices added */
  CS_REFINE_POLYGON_T,   /*!< refine polygon by triangulation with added
                           center vertex */
  CS_REFINE_POLYGON_Q,   /*!< refine polygon by quadrangulation with added
                           center vertex */

  CS_REFINE_TETRA,       /*!< simple tetrahedron subdivision scheme */
  CS_REFINE_TETRA_H,     /*!< simple tetrahedron to hexahedron scheme */
  CS_REFINE_PYRAM,       /*!< simple pyramid subdivision scheme */
  CS_REFINE_PRISM,       /*!< simple prism subdivision scheme */
  CS_REFINE_HEXA,        /*!< simple hexahedron subdivision scheme */

  CS_REFINE_POLYHEDRON,     /*!< refine polyhedron by addition of center vertex
                              and tetrahedra/pyramids joining faces to center,
                              with an attempt at */
  CS_REFINE_POLYHEDRON_P,  /*!< refine polyhedron by addition of center vertex
                             and tetrahedra/pyramids joining faces to center */

  CS_REFINE_N_TYPES      /*!< number of available refinement templates */

} cs_mesh_refine_type_t;

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

typedef struct {

  cs_lnum_t v0;    // id of first vertex
  cs_lnum_t v1;    // id of second vertex
  int16_t   f0;    // local id of first face
  int16_t   f1;    // local id of second face
  int16_t   c0;    // local id of first sub-cell
  int16_t   c1;    // local id of second sub-cell

} _mrh_edge_t;

/*----------------------------------------------------------------------------
 * Compare mesh refinement helper edges (qsort function).
 *
 * parameters:
 *   a <-> pointer to first edge
 *   b <-> pointer to second edge
 *
 * returns:
 *   result of lexicographical comparison of edge vertices
 *----------------------------------------------------------------------------*/

static int
_rh_edge_compare(const void  *a,
                 const void  *b)
{
  int retval = 1;

  const _mrh_edge_t *e_a = (const _mrh_edge_t *)a;
  const _mrh_edge_t *e_b = (const _mrh_edge_t *)b;

  if (e_a->v0 < e_b->v0)
    retval = -1;
  else if (e_a->v0 == e_b->v0) {
    if (e_a->v1 < e_b->v1)
      retval = -1;
    else if (e_a->v1 == e_b->v1)
      retval = 0;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Compare mesh refinement helper edges-sub-cells connectivity
 * (qsort function).
 *
 * parameters:
 *   a <-> pointer to first edge
 *   b <-> pointer to second edge
 *
 * returns:
 *   result of lexicographical comparison of edge vertices
 *----------------------------------------------------------------------------*/

static int
_rh_edge_compare_s_cells(const void  *a,
                         const void  *b)
{
  int retval = 1;

  const _mrh_edge_t *e_a = (const _mrh_edge_t *)a;
  const _mrh_edge_t *e_b = (const _mrh_edge_t *)b;

  int16_t c_a[2], c_b[2];
  if (e_a->c0 < e_a->c1) {
    c_a[0] = e_a->c0;
    c_a[1] = e_a->c1;
  }
  else {
    c_a[0] = e_a->c1;
    c_a[1] = e_a->c0;
  }
  if (e_b->c0 < e_b->c1) {
    c_b[0] = e_b->c0;
    c_b[1] = e_b->c1;
  }
  else {
    c_b[0] = e_b->c1;
    c_b[1] = e_b->c0;
  }

  if (c_a[0] < c_b[0])
    retval = -1;
  else if (c_a[0] == c_b[0]) {
    if (c_a[1] < c_b[1])
      retval = -1;
    else if (c_a[1] == c_b[1]) {
      retval = 0;
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * State of current face merging (working structure)
 *----------------------------------------------------------------------------*/

class cs_cell_refine_helper {

public:

  int16_t       n_sub_cells;      /* current triangle vertices list */

  int16_t       n_faces;          /* number of exterior faces */
  int           n_edges;          /* number of exterior edges */
  int           n_vtx;            /* number of vertices */
  int           c_i_f2v_size;     /* size of interior faces to vertex
                                     connectivity (includes extra -1 per face) */

  int           n_faces_max;      /* current max face allocation */
  int           n_edges_max;      /* current max edges allocation */
  int           n_vtx_max;        /* max number of vertices */

  cs_lnum_t     c_vtx_id;         /* center vertex id, or -1 */

  char          f_r_range[2];     /* min and max face refinement levels */

  int16_t      *sc_equiv;         /* equivalence for sub-cells */
  int16_t      *e_f2c;            /* exterior face to sub-cell id */
  _mrh_edge_t  *edges;            /* edges -> vertices and exterior faces */

  cs_lnum_t    *vtx_id;           /* vertices list */
  int          *vtx_valence;      /* vertices valence */

  /*--------------------------------------------------------------------------
   * Constructor
   *--------------------------------------------------------------------------*/

  cs_cell_refine_helper()
  {
    n_sub_cells = 0;
    n_faces = 0;
    n_edges = 0;
    n_vtx = 0;
    c_i_f2v_size = 0;
    n_faces_max = 32;
    n_edges_max = 64;
    n_vtx_max = 32;
    f_r_range[0] = 0, f_r_range[1] = 0;

    c_vtx_id = -1;

    CS_MALLOC(sc_equiv, n_faces_max, int16_t);
    CS_MALLOC(e_f2c, n_faces_max, int16_t);

    CS_MALLOC(edges, n_edges_max, _mrh_edge_t);

    CS_MALLOC(vtx_id, n_vtx_max, cs_lnum_t);
    CS_MALLOC(vtx_valence, n_vtx_max, int);
  }

  /*--------------------------------------------------------------------------
   * Destructor
   *--------------------------------------------------------------------------*/

  ~cs_cell_refine_helper()
  {
    CS_FREE(vtx_valence);
    CS_FREE(vtx_id);

    CS_FREE(edges);

    CS_FREE(e_f2c);
    CS_FREE(sc_equiv);
  }

  /*--------------------------------------------------------------------------
   * Reset for new cell
   *--------------------------------------------------------------------------*/

  void
  reset_for_new_cell()
  {
    n_sub_cells = 0;
    n_faces = 0;
    n_edges = 0;
    n_vtx = 0;
    c_i_f2v_size = 0;
    f_r_range[0] = 0, f_r_range[1] = 0;
  }

  /*--------------------------------------------------------------------------
   * Set central vertex id
   *--------------------------------------------------------------------------*/

  inline void
  set_center_vtx(cs_lnum_t  id)
  {
    c_vtx_id = id;
  }

  /*--------------------------------------------------------------------------
   * Add a face
   *--------------------------------------------------------------------------*/

  void
  add_face(int        orient,       // 1 for outwards normal, -1 for inwards
           char       r_level,      // face refinement level
           cs_lnum_t  n_f_vtx,      // number of face vertices
           const cs_lnum_t  f_vtx_ids[])  // face vertex ids
  {
    // Reallocate if needed (with margin)

    if (n_faces + n_f_vtx > n_faces_max) {
      n_faces_max *= 2;
      CS_REALLOC(e_f2c, n_faces_max, int16_t);
      CS_MALLOC(sc_equiv, n_faces_max, int16_t);
    }
    if (n_edges + n_f_vtx > n_edges_max) {
      n_edges_max *= 2;
      CS_REALLOC(edges, n_edges_max, _mrh_edge_t);
    }
    if (n_vtx + n_f_vtx > n_vtx_max) {
      n_vtx_max *= 2;
      CS_REALLOC(vtx_id, n_vtx_max, cs_lnum_t);
      CS_REALLOC(vtx_valence, n_vtx_max, int);
    }

    e_f2c[n_faces] = n_faces;

    if (n_faces == 0) {
      f_r_range[0] = r_level;
      f_r_range[1] = r_level;
    }
    else {
      if (r_level < f_r_range[0])
        f_r_range[0] = r_level;
      if (r_level > f_r_range[1])
        f_r_range[1] = r_level;
    }

    cs_lnum_t s_id = n_edges;
    cs_lnum_t e_id = s_id;
    if (orient == 1) {
      for (cs_lnum_t i = 0; i < n_f_vtx; i++) {
        edges[e_id].v0 = f_vtx_ids[i];
        edges[e_id].v1 = f_vtx_ids[(i+1)%n_f_vtx];
        edges[e_id].f0 = n_faces;
        edges[e_id].f1 = -1;
        edges[e_id].c0 = n_faces;
        edges[e_id].c1 = -1;
        e_id++;
      }
    }
    else {
      for (cs_lnum_t i = n_f_vtx-1; i > -1; i--) {
        edges[e_id].v0 = f_vtx_ids[(i+1)%n_f_vtx];
        edges[e_id].v1 = f_vtx_ids[i];
        edges[e_id].f0 = n_faces;
        edges[e_id].f1 = -1;
        edges[e_id].c0 = n_faces;
        edges[e_id].c1 = -1;
        e_id++;
      }
    }
    n_edges = e_id;

    for (cs_lnum_t i = s_id; i < e_id; i++) {
      if (edges[i].v1 < edges[i].v0) {
        _mrh_edge_t es = edges[i];
        edges[i].v0 = edges[i].v1;
        edges[i].v1 = es.v0;
        edges[i].f0 = edges[i].f1;
        edges[i].f1 = es.f0;
        edges[i].c0 = edges[i].c1;
        edges[i].c1 = es.c0;
      }
    }

    n_faces += 1;
    n_sub_cells += 1;
  }

  /*--------------------------------------------------------------------------
   * Add interior face for cell refinement with center vertex scheme
   *
   * \param[in]   s_id             start edge id
   * \param[in]   e_id             past the end edge id
   * \param[out]  c_n_s_cells      number of sub-cells
   * \param[out]  c_n_i_faces      number of added interior faces
   * \param[out]  c_i_faces_size   size of added interior faces connectivity
   * \param[out]  i_f2c_pre        temp. cell to interior faces connect
   * \param[out]  i_f2v_pre        temp. interior faces to vertices connect,
   *                               with each face "closed" with a -1 value.
   *--------------------------------------------------------------------------*/

  void
  build_interior_face_with_center_vertex(cs_lnum_t   s_id,
                                         cs_lnum_t   e_id,
                                         cs_lnum_t   &c_n_i_faces,
                                         cs_lnum_t   &c_i_faces_size,
                                         int16_t      i_f2c_pre[][2],
                                         cs_lnum_t    i_f2v_pre[])
  {
    cs_lnum_t *f2v_pre = i_f2v_pre + c_i_f2v_size;

    int e0 = 0;
    cs_lnum_t i_face_size = e_id - s_id + 2;

    cs_lnum_t v0 = cs::min(edges[s_id].v0, edges[s_id].v1);
    f2v_pre[0] = v0;

    /* In case of single edge, simply start face with smallest vertex id */

    if (i_face_size == 3) {
      _mrh_edge_t *e = edges + s_id;
      if (e->v0 == v0) {
        f2v_pre[1] = e->v1;
      }
      else {
        f2v_pre[1] = e->v0;
      }
      f2v_pre[2] = c_vtx_id;
    }

    /* In case of multiple edges sharing faces, link edges from s_id to e_id;
       first (reference) edge is that with lowest vertex id */

    else {
      cs_lnum_t j = 1;
      int n_e = e_id - s_id;
      constexpr int n_f_vtx_max = 16;
      cs_lnum_t edge_queue[n_f_vtx_max];

      /* Search for vertex with smallest id */
      for (cs_lnum_t i = s_id+1; i < e_id; i++) {
        cs_lnum_t v1 = cs::min(edges[i].v0, edges[i].v1);
        if (v1 < v0) {
          v0 = v1;
          e0 = i - s_id;
        }
      }

      cs_lnum_t v_s = v0;

      int edge_queue_size = n_e;
      for (int i = 0; i < n_e; i++)
        edge_queue[i] = s_id + (e0 + i)%n_e;

      /* Now build list of vertices, searching for next vertex in queue */
      while (edge_queue_size > 0) {
        cs_lnum_t v1 = -1;
        for (int i = 0; i < edge_queue_size; i++) {
          _mrh_edge_t *e = edges + edge_queue[i];
          if (e->v0 == v0)
            v1 = e->v1;
          else if (e->v1 == v0)
            v1 = e->v0;
          if (v1 != -1) {
            f2v_pre[j] = v1;
            j++;
            edge_queue[i] = edge_queue[edge_queue_size - 1];
            edge_queue_size -= 1;
            break;
          }
        }
        v0 = v1;
        if (v1 == -1)
          break;
      }

      /* Now continue contour with cell center */

      f2v_pre[j] = c_vtx_id;

      /* If some edges are still not connected, they should
         loop back to the starting point */

      v0 = v_s;
      cs_lnum_t k = i_face_size -1;
      while (edge_queue_size > 0) {
        cs_lnum_t v1 = -1;
        for (int i = 0; i < edge_queue_size; i++) {
          _mrh_edge_t *e = edges + edge_queue[i];
          if (e->v0 == v0)
            v1 = e->v1;
          else if (e->v1 == v0)
            v1 = e->v0;
          if (v1 != -1) {
            f2v_pre[k] = v1;
            k--;
            edge_queue[i] = edge_queue[edge_queue_size - 1];
            edge_queue_size -= 1;
            break;
          }
        }
        v0 = v1;
        if (v1 == -1)
          bft_error(__FILE__, __LINE__, 0,
                    _("Error connecting edges with adjacent faces."));
      }
      assert(k == j);
    }

    /* Mark end of face */

    f2v_pre[i_face_size] = -1;

    /* Update indexes and metadata */

    _mrh_edge_t *e = edges + s_id + (cs_lnum_t)e0;
    cs_lnum_t k = (e->v1 > e->v0) ? 1 : 0;

    i_f2c_pre[c_n_i_faces][(0+k)%2] = e->c0;
    i_f2c_pre[c_n_i_faces][(1+k)%2] = e->c1;

    n_sub_cells = cs::max(n_sub_cells,
                          (int16_t)(cs::max(e->c0, e->c1) + 1));

    c_i_f2v_size += i_face_size+1; // size in temporary connectivity
    c_i_faces_size += i_face_size; // size in final connectivity
    c_n_i_faces += 1;
  }

  /*--------------------------------------------------------------------------
   * Merge equivalent sub-cells.
   *
   * \param[in]  vtx_r_gen  refinedment level associated to each vertex.
   *--------------------------------------------------------------------------*/

  void
  merge_equivalent_sub_cells()
  {
    /* First pass: ensure equivalences point to anchor for each sub-cell */

    int16_t count = 0;

    for (int16_t i = 0; i < n_sub_cells; i++) {
      int16_t j = sc_equiv[i];
      if (j == i)
        count++;
      else {
        int16_t k = sc_equiv[j];
        while (k < j) {
          j = k;
          sc_equiv[i] = j;
          k = sc_equiv[j];
        }
      }
    }

    /* Return if no cells need to be merged */

    if (count == n_sub_cells)
      return;

    /* Second pass: compact numbering for anchor (sub-cell with no lower
       id equivalent) of each sub-cell. Use negative numbers to mark
       sub-cell as anchor. */

    count = 1;
    for (int16_t i = 1; i < n_sub_cells; i++) {
      if (sc_equiv[i] == i) {
        sc_equiv[i] = - count - 1;
        count++;
      }
    }

    /* Third pass: update equivalences (point to root) for each sub-cell */

    for (int16_t i = 0; i < n_sub_cells; i++) {
      int16_t j = sc_equiv[i];
      if (j < 0)
        sc_equiv[i] = -j -1;
      else {
        j = sc_equiv[i]; // points to anchor
        assert(j > -1);
        sc_equiv[i] = sc_equiv[j];
      }
    }

    /* Now update edges, removing those pointing to equivalent cells,
       since no interior face will be based on them. */

    cs_lnum_t j = 0;
    for (cs_lnum_t i = 0; i < n_edges; i++) {
      _mrh_edge_t *e = edges + i;
      e->c0 = sc_equiv[e->c0];
      e->c1 = sc_equiv[e->c1];
      if (e->c0 != e->c1) {
        edges[j] = edges[i];
        j++;
      }
    }

    n_edges = j;
    n_sub_cells = count;

    /* Finally, update exterior faces to cells adjacency */

    for (int16_t i = 0; i < n_faces; i++)
      e_f2c[i] = sc_equiv[e_f2c[i]];
  }

  /*--------------------------------------------------------------------------
   * Remove edges connecting only vertices of higher refinement level
   * than local cells.
   *
   * \param[in]  vtx_r_gen  refinedment level associated to each vertex.
   *--------------------------------------------------------------------------*/

  void
  filter_high_refinement_level_edges(const char  vtx_r_gen[])
  {
    char r_threshold = f_r_range[0];

    for (cs_lnum_t i = 0; i < n_edges; i++) {
      _mrh_edge_t *e = edges + i;
      if (   vtx_r_gen[e->v0] > r_threshold
          && vtx_r_gen[e->v1] > r_threshold) {
        if (e->c0 > e->c1)
          sc_equiv[e->c0] = sc_equiv[e->c1];
        else
          sc_equiv[e->c1] = sc_equiv[e->c0];
      }
    }
  }

  /*--------------------------------------------------------------------------
   * Remove edges connected to corners present at lower refinement level
   * of local cells.
   *
   * \param[in]  vtx_r_gen  refinedment level associated to each vertex.
   *--------------------------------------------------------------------------*/

  void
  filter_corner_edges(const char  vtx_r_gen[])
  {
    char r_threshold = f_r_range[0];

    for (cs_lnum_t i = 0; i < n_edges; i++) {
      _mrh_edge_t *e = edges + i;
      if (   vtx_r_gen[e->v0] < r_threshold
          || vtx_r_gen[e->v1] < r_threshold) {
        if (e->c0 > e->c1)
          sc_equiv[e->c0] = sc_equiv[e->c1];
        else
          sc_equiv[e->c1] = sc_equiv[e->c0];
      }
    }
  }

  /*--------------------------------------------------------------------------
   * Refine cell using scheme with central vertex
   *
   * \param[in]   c_r_flag         cell refinement type flag
   * \param[in]   vtx_r_gen        refinement level associated to each vertex.
   * \param[out]  c_n_s_cells      number of sub-cells
   * \param[out]  c_n_i_faces      number of added interior faces
   * \param[out]  c_i_faces_size   size of added interior faces connectivity
   * \param[out]  e_f2c_pre        temp. exterior faces to cell connect
   *                               (size: n_faces member value)
   * \param[out]  i_f2c_pre        temp. interior faces to cells connect
   * \param[out]  i_f2v_pre        temp. interior faces to vertices connect,
   *                               with each face "closed" with a -1 value.
   *--------------------------------------------------------------------------*/

  void
  refine_cell_with_center_vertex(cs_mesh_refine_type_t   c_r_flag,
                                 const char              vtx_r_gen[],
                                 cs_lnum_t              &c_n_s_cells,
                                 cs_lnum_t              &c_n_i_faces,
                                 cs_lnum_t              &c_i_faces_size,
                                 int16_t                 e_f2c_pre[],
                                 int16_t                 i_f2c_pre[][2],
                                 cs_lnum_t               i_f2v_pre[])
  {
    cs_lnum_t c_n_s_cells_max = c_n_s_cells;
    cs_lnum_t c_n_i_faces_max = c_n_i_faces;
    cs_lnum_t c_i_faces_size_max = c_i_faces_size;

    c_n_s_cells = 0;
    c_n_i_faces = 0;
    c_i_faces_size = 0;

    // Reorder and merge edges
    // -----------------------

    qsort(edges, n_edges, sizeof(_mrh_edge_t), &_rh_edge_compare);

    cs_lnum_t count = (n_edges > 0) ? 1 : 0;

    for (cs_lnum_t i = 1; i < n_edges; i++) {
      _mrh_edge_t *e_p = edges + (count - 1);
      if (edges[i].v0 == e_p->v0 && edges[i].v1 == e_p->v1) {
        if (edges[i].f0 != -1) {
          assert(edges[i].f1 == -1);
          cs_assert(e_p->f0 == -1 && e_p->f1 != -1);
          e_p->f0 = edges[i].f0;
          e_p->c0 = edges[i].c0;
        }
        else {
          assert(edges[i].f1 != -1);
          cs_assert(e_p->f0 != -1 && e_p->f1 == -1);
          e_p->f1 = edges[i].f1;
          e_p->c1 = edges[i].c1;
        }
      }
      else {
        edges[count] = edges[i];
        count++;
      }
    }
    n_edges = count;

    // Filter edges which for which no faces are generated
    // ---------------------------------------------------
    // Associated cell indexes need to be updated, as this
    // amounts to merging 2 adjacent pyramidal cells.

    for (int16_t i = 0; i < n_sub_cells; i++)
      sc_equiv[i] = i;

    filter_high_refinement_level_edges(vtx_r_gen);

    if (c_r_flag == CS_REFINE_POLYHEDRON) // not CS_REFINE_POLYHEDRON_P
      filter_corner_edges(vtx_r_gen);

    merge_equivalent_sub_cells();

    // Now reorder edges based on sub-cells adjacency
    // -----------------------------------------------

    qsort(edges, n_edges, sizeof(_mrh_edge_t), &_rh_edge_compare_s_cells);

    // Now buid faces

    c_i_faces_size = 0;

    cs_lnum_t s_id = 0;
    for (cs_lnum_t i = 1; i < n_edges; i++) {
      _mrh_edge_t *e_c = edges + i;
      _mrh_edge_t *e_p = edges + (i - 1);
      if (_rh_edge_compare_s_cells(e_c, e_p) != 0) {
        build_interior_face_with_center_vertex(s_id, i,
                                               c_n_i_faces,
                                               c_i_faces_size,
                                               i_f2c_pre,
                                               i_f2v_pre);
        s_id = i;
      }
    }
    if (s_id < n_edges-1)
      build_interior_face_with_center_vertex(s_id, n_edges,
                                             c_n_i_faces,
                                             c_i_faces_size,
                                             i_f2c_pre,
                                             i_f2v_pre);

    c_n_s_cells = n_sub_cells;

    for (int16_t i = 0; i < n_faces; i++)
      e_f2c_pre[i] = e_f2c[i];

    assert(c_n_s_cells_max >= c_n_s_cells);
    assert(c_n_i_faces_max >= c_n_i_faces);
    assert(c_i_faces_size_max >= c_i_faces_size);
  }
};

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer to subdivide a given cell, updating faces.
 *
 * The cells->faces adjacency is updated as faces are locally renumbered
 * in canonical order for known cell types.
 *
 * \param[in]  m                  pointer to mesh structure
 * \param[in]  cell_id            (old) cell id (0 to n-1)
 * \param[in]  n_b_f_ini          old number of boundary faces
 * \param[in]  c_o2n_idx          old to new cells index
 * \param[in]  i_face_o2n_idx     old to new interior faces index
 * \param[in]  b_face_o2n_idx     old to new boundary faces index
 * \param[in]  c2f                cells->faces adjacency (boundary faces first)
 * \param[in]  c2f2v_start        start id for adjacent face vertices
 *                                definitions (initialized to 0)
 * \param[in]  c_v_idx            for each cell, start index of added vertices
 * \param[in]  c_f_n_idx          cells to new faces index (count in),
 *                                or nullptr
 */
/*----------------------------------------------------------------------------*/

typedef void
(subdivide_cell_t)(const cs_mesh_t              *m,
                   cs_lnum_t                     cell_id,
                   cs_lnum_t                     n_b_f_ini,
                   const cs_lnum_t               c_o2n_idx[],
                   const cs_lnum_t               i_face_o2n_idx[],
                   const cs_lnum_t               b_face_o2n_idx[],
                   cs_adjacency_t               *c2f,
                   const cs_lnum_t               c2f2v_start[],
                   const cs_lnum_t               c_v_idx[],
                   const cs_lnum_t               c_f_n_idx[]);

/*============================================================================
 * Static global variables.
 *============================================================================*/

constexpr int _n_f_vtx_max = 96;

static cs_mesh_refine_type_t _refine_tria_type = CS_REFINE_TRIA;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Transform a counts array to an index.
 *
 * The counts array values are assumed shifted by one, that is
 * count[i] is stored in elt_idx[i+1]
 *
 * The value of elt_idx[0] is unchanged, so can be set to 0 or another
 * (shift) value by the caller.
 *
 * \param[in]       n_elts      number of elements
 * \param[in, out]  elt_idx     counts in, index out
 */
/*----------------------------------------------------------------------------*/

static void
_counts_to_index(cs_lnum_t   n_elts,
                 cs_lnum_t   elt_idx[])
{
  for (cs_lnum_t i = 0; i < n_elts; i++)
    elt_idx[i+1] += elt_idx[i];
}

/*----------------------------------------------------------------------------
 * Print information on a mesh structure.
 *
 * parameters:
 *   m     <--  pointer to mesh structure.
 *   name  <--  associated name.
 *----------------------------------------------------------------------------*/

static void
_print_mesh_counts(const cs_mesh_t  *m,
                   const char       *name)
{
  cs_log_printf(CS_LOG_DEFAULT, "\n");
  cs_log_printf(CS_LOG_DEFAULT,
                _(" %s\n"
                  "     Number of cells:          %llu\n"
                  "     Number of interior faces: %llu\n"
                  "     Number of boundary faces: %llu\n"
                  "     Number of vertices:       %llu\n"),
             name,
             (unsigned long long)(m->n_g_cells),
             (unsigned long long)(m->n_g_i_faces),
             (unsigned long long)(m->n_g_b_faces - m->n_g_free_faces),
             (unsigned long long)(m->n_g_vertices));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Determine face refinement level based on associated vertex
 *        levels.
 *
 * A refinement level is different from a refinement generation:
 * - The generation corresponds to the cell refinement level at which a face
 *   is added.
 * - The level corresponds to the number of subdivisions of the root
 *   ancestor face.
 *
 * Faces based on previous refinement should usually be triangles (with
 * 1 edge/2 consecutive vertices of matching refinement level) or quadrangles
 * (with 2 edges/3 consecutive vertices of matching refinement level).
 *
 * When a face is polygonal, we must try to identify corners, as some edges
 * may have been split due to adjacent face/cell refinement. When an edge
 * is split in this way, the added vertex is of a higher level than both
 * adjacent vertices. We use this property to try to reconstruct the original
 * topology.
 *
 * \param[in, out]  n_v        number of face vertices
 * \param[in]       face_vtx   face vertex idx
 * \param[in]       vtx_r_gen  refinedment level associated to each vertex.
 */
/*----------------------------------------------------------------------------*/

static char
_compute_face_r_level(cs_lnum_t        n_v,
                      const cs_lnum_t  face_vtx[],
                      const char       vtx_r_gen[])
{
  char retval = 0;

  char r_gen[256];
  char *_r_gen = r_gen;

  if (n_v > 256) {
    /* If we ever need to do this, mesh quality is probably quite bad... */
    CS_MALLOC(_r_gen, 256, char);
  }

  for (cs_lnum_t i = 0; i < n_v; i++)
    _r_gen[i] = vtx_r_gen[face_vtx[i]];

  cs_lnum_t _n_v = n_v, v_count = n_v;
  do {
    _n_v = v_count;
    v_count = 0;
    bool prev_up = (_r_gen[0] > _r_gen[n_v - 1]) ? true : false;

    for (cs_lnum_t i = 0; i < _n_v; i++) {
      cs_lnum_t i_n = (i+1)%_n_v;
      if (prev_up && _r_gen[i] > _r_gen[i_n]) { /* Remove vertex */
        prev_up = false;
      }
      else {
        prev_up = (_r_gen[i_n] > _r_gen[i]) ? true : false;
        _r_gen[v_count++] = _r_gen[i]; /* Keep vertex */
      }
    }
  } while (v_count < _n_v);

  for (cs_lnum_t i = 0; i < _n_v; i++)
    retval = cs::max(retval, _r_gen[i]);

  if (_r_gen != r_gen)
    CS_FREE(_r_gen);

  return retval;
}

/*----------------------------------------------------------------------------
 * Check if a polygon is convex or at least star-shaped.
 *
 * parameters:
 *   n_face_vertices <-- number of polygon vertices
 *   f2v_lst         <-- polygon connectivity
 *   vtx_coords      <-- vertex coordinates
 *
 * returns:
 *   true if quadrangle is convex, false otherwise
 *----------------------------------------------------------------------------*/

static bool
_polygon_is_nearly_convex(cs_lnum_t         n_face_vertices,
                          const cs_lnum_t   f2v_lst[],
                          const cs_coord_t  vtx_coords[][3])
{
  bool is_convex = true;

  /* Compute approximate face center coordinates for the polygon */

  cs_real_t a_center[3] = {0, 0, 0};
  cs_real_t f_norm[3] = {0, 0, 0};

  for (cs_lnum_t j = 0; j < n_face_vertices; j++) {
    const cs_lnum_t v0 = f2v_lst[j];
    for (cs_lnum_t i = 0; i < 3; i++)
      a_center[i] += vtx_coords[v0][i];
  }

  for (cs_lnum_t i = 0; i < 3; i++)
    a_center[i] /= n_face_vertices;

  /* first loop on edges, with implied subdivision into triangles
     defined by edge and cell center */

  cs_real_t vc0[3], vc1[3], vn[3];

  for (cs_lnum_t tri_id = 0; tri_id < n_face_vertices; tri_id++) {

    const cs_lnum_t v0 = f2v_lst[tri_id];
    const cs_lnum_t v1 = f2v_lst[(tri_id+1)%n_face_vertices];

    for (cs_lnum_t i = 0; i < 3; i++) {
      vc0[i] = vtx_coords[v0][i] - a_center[i];
      vc1[i] = vtx_coords[v1][i] - a_center[i];
    }

    cs_math_3_cross_product(vc0, vc1, vn);

    for (cs_lnum_t i = 0; i < 3; i++)
      f_norm[i] += vn[i];

  }

  /* second loop on edges, to check consistency of orientation */

  for (cs_lnum_t tri_id = 0; tri_id < n_face_vertices; tri_id++) {

    const cs_lnum_t v0 = f2v_lst[tri_id];
    const cs_lnum_t v1 = f2v_lst[(tri_id+1)%n_face_vertices];

    for (cs_lnum_t i = 0; i < 3; i++) {
      vc0[i] = vtx_coords[v0][i] - a_center[i];
      vc1[i] = vtx_coords[v1][i] - a_center[i];
    }

    cs_math_3_cross_product(vc0, vc1, vn);

    if (cs_math_3_dot_product(vn, f_norm) < 0) {
      is_convex = false;
      break;
    }

  }

  return is_convex;
}

/*----------------------------------------------------------------------------
 * Tetrahedron reconstruction.
 *
 * parameters:
 *   vtx_tria    <-- triangular faces connectivity
 *   f_id_tria   <-> matching (triangle) face ids, in initial order
 *                   set to -1 when matched
 *   f_sgn_tria  <-- matching triangle face orientation relative to cell
 *   f_ids       --> matching face ids, in template order
 *   f_sgn       --> matching face orientation relative to cell
 *   c2f2v_start --> matching initial vertex ids
 *
 * return:
 *   true if reconstruction successful, false otherwise
 *----------------------------------------------------------------------------*/

static bool
_match_tetra(cs_lnum_t  vtx_tria[4][3],
             cs_lnum_t  f_id_tria[4],
             int        f_sgn_tria[4],
             cs_lnum_t  f_ids[4],
             int        f_sgn[4],
             cs_lnum_t  c2f2v_start[4])
{
  bool matched = true;

  /*
    Base : vertices of the first face; we take the vertices in opposite
    order ("bottom" face numbering of the tetrahedron with outward
    pointing normal).
                                                  Vertices
    *        x 3                         Face 0:  0 2 1
    *       /|\                          Face 1:  0 1 3
    *      / | \                         Face 2:  1 2 3
    *     /  |  \                        Face 3:  2 0 3
    *  0 x- -|- -x 2
    *     \  |  /
    *      \ | /
    *       \|/
    *        x 1
  */

  cs_lnum_t  cell_vtx[4];

  cell_vtx[0] = vtx_tria[0][0];
  cell_vtx[1] = vtx_tria[0][2];
  cell_vtx[2] = vtx_tria[0][1];
  cell_vtx[3] = -1;

  f_ids[0] = f_id_tria[0];
  f_sgn[0] = f_sgn_tria[0];
  c2f2v_start[0] = 0;
  f_id_tria[0] = -1;

  /*
    We start with 3 of 4 vertices found; all other triangles should share
    vertex 4, and one edge of the base triangle.
  */

  for (cs_lnum_t i = 1; i < 4; i++) {

    cs_lnum_t v_id_0 = cell_vtx[i-1];
    cs_lnum_t v_id_1 = cell_vtx[i%3];

    f_ids[i] = -1;

    for (cs_lnum_t f_id = 1; f_id < 4 && f_ids[i] == -1; f_id++) {
      if (f_id_tria[f_id] < 0)
        continue;
      cs_lnum_t v_id_3 = -1;
      for (cs_lnum_t j = 0; j < 3; j++) {
        if (vtx_tria[f_id][j] == v_id_0) {
          if (vtx_tria[f_id][((j+1) % 3)] == v_id_1) {
            v_id_3 = vtx_tria[f_id][((j+2) % 3)];
            if (i == 1)
              cell_vtx[3] = v_id_3;
            else if (cell_vtx[3] != v_id_3)
              matched = false; /* error case */
            f_ids[i] = f_id_tria[f_id];
            f_sgn[i] = f_sgn_tria[f_id];
            f_id_tria[f_id] = -1;
            c2f2v_start[i] = j;
            break;
          }
        }
      }
    }

    if (f_ids[i] < 0)
      matched = false;
  }

  return matched;
}

/*----------------------------------------------------------------------------
 * Pyramid reconstruction.
 *
 * parameters:
 *   vtx_tria    <-- triangular faces connectivity
 *   vtx_quad    <-- quadrangle faces connectivity
 *   f_id_tria   <-> matching (triangle) face ids, in initial order
 *                   set to -1 when matched
 *   f_id_quad   <-> matching (quadrangle) face ids, in initial order
 *                   set to -1 when matched
 *   f_sgn_tria  <-- matching triangle face orientation relative to cell
 *   f_sgn_quad  <-- matching quadrangle face orientation relative to cell
 *   f_ids       --> matching face ids, in template order
 *   f_sgn       --> matching face orientation relative to cell
 *   c2f2v_start --> matching initial vertex ids
 *
 * return:
 *   true if reconstruction successful, false otherwise
 *----------------------------------------------------------------------------*/

static bool
_match_pyram(cs_lnum_t  vtx_tria[4][3],
             cs_lnum_t  vtx_quad[1][4],
             cs_lnum_t  f_id_tria[4],
             cs_lnum_t  f_id_quad[1],
             int        f_sgn_tria[4],
             int        f_sgn_quad[4],
             cs_lnum_t  f_ids[5],
             int        f_sgn[5],
             cs_lnum_t  c2f2v_start[5])
{
  bool matched = true;

  /*
    Base : vertices of the quadrangle; we take the vertices in opposite
    order ("bottom" face numbering of the pyramid with outward
    pointing normal).
                                                  Vertices
    *         4 x                        Face 0:  0 3 2 1
    *          /|\                       Face 1:  0 1 4
    *         //| \                      Face 2:  1 2 4
    *        // |  \                     Face 3:  2 3 4
    *     3 x/--|---x 2                  Face 4:  3 0 4
    *      //   |  /
    *     //    | /
    *  0 x-------x 1
  */

  cs_lnum_t  cell_vtx[5];

  cell_vtx[0] = vtx_quad[0][0];
  cell_vtx[1] = vtx_quad[0][3];
  cell_vtx[2] = vtx_quad[0][2];
  cell_vtx[3] = vtx_quad[0][1];
  cell_vtx[4] = -1;

  f_ids[0] = f_id_quad[0];
  f_sgn[0] = f_sgn_quad[0];
  c2f2v_start[0] = 0;
  f_id_quad[0] = -1;

  /*
    We have found 4 out of 5 vertices; all 4 triangles should share
    vertex 4, and one edge of the base quadrangle.
  */

  for (cs_lnum_t i = 1; i < 5; i++) {

    cs_lnum_t v_id_0 = cell_vtx[i-1];
    cs_lnum_t v_id_1 = cell_vtx[i%4];

    f_ids[i] = -1;

    for (cs_lnum_t f_id = 1; f_id < 5 && f_ids[i] == -1; f_id++) {
      cs_lnum_t t_f_id = f_id-1;
      if (f_id_tria[t_f_id] < 0)
        continue;
      cs_lnum_t v_id_4 = -1;
      for (cs_lnum_t j = 0; j < 3; j++) {
        if (vtx_tria[t_f_id][j] == v_id_0) {
          if (vtx_tria[t_f_id][((j+1) % 3)] == v_id_1) {
            v_id_4 = vtx_tria[t_f_id][((j+2) % 3)];
            if (i == 1)
              cell_vtx[4] = v_id_4;
            else if (cell_vtx[4] != v_id_4)
              matched = false; /* error case */
            f_ids[i] = f_id_tria[t_f_id];
            f_sgn[i] = f_sgn_tria[t_f_id];
            c2f2v_start[i] = j;
            f_id_tria[t_f_id] = -1;
            break;
          }
        }
      }
    }

    if (f_ids[i] < 0)
      matched = false;
  }

  return matched;
}

/*----------------------------------------------------------------------------
 * Prism (pentahedron) reconstruction
 *
 * parameters:
 *   vtx_tria    <-- triangular faces connectivity
 *   vtx_quad    <-- quadrangle faces connectivity
 *   f_id_tria   <-> matching (triangle) face ids, in initial order
 *                   set to -1 when matched
 *   f_id_quad   <-> matching (quadrangle) face ids, in initial order
 *                   set to -1 when matched
 *   f_sgn_tria  <-- matching triangle face orientation relative to cell
 *   f_sgn_quad  <-- matching quadrangle face orientation relative to cell
 *   f_ids       --> matching face ids, in template order
 *   f_sgn       --> matching face orientation relative to cell
 *   c2f2v_start --> matching initial vertex ids
 *
 * return:
 *   true if reconstruction successful, false otherwise
 *----------------------------------------------------------------------------*/

static bool
_match_prism(cs_lnum_t  vtx_tria[2][3],
             cs_lnum_t  vtx_quad[3][4],
             cs_lnum_t  f_id_tria[2],
             cs_lnum_t  f_id_quad[3],
             int        f_sgn_tria[4],
             int        f_sgn_quad[4],
             cs_lnum_t  f_ids[5],
             int        f_sgn[5],
             cs_lnum_t  c2f2v_start[5])
{
  bool matched = true;

  /*
    Base : vertices of the first triangle; we take the vertices in opposite
    order ("bottom" face numbering of the prism with outward
    pointing normal).

                                                  Vertices
    *  3 x-------x 5                     Face 0:  0 1 4 3
    *    |\     /|                       Face 1:  1 2 5 4
    *    | \   / |                       Face 2:  2 0 3 5
    *  0 x- \-/ -x 2                     Face 3:  0 2 1
    *     \ 4x  /                        Face 4:  3 4 5
    *      \ | /
    *       \|/
    *        x 1
  */

  cs_lnum_t  cell_vtx[6];

  cell_vtx[0] = vtx_tria[0][0]; /* base: Face 3 */
  cell_vtx[1] = vtx_tria[0][2];
  cell_vtx[2] = vtx_tria[0][1];
  cell_vtx[3] = -1;
  cell_vtx[4] = -1;
  cell_vtx[5] = -1;

  f_ids[3] = f_id_tria[0];
  f_sgn[3] = f_sgn_tria[0];
  c2f2v_start[3] = 0;
  f_id_tria[0] = -1;

  /*
    We have found 3 out of 6 vertices; all 3 quadrangles should share
    one edge of the base quadrangle.
  */

  for (cs_lnum_t i = 0; i < 3; i++) {

    cs_lnum_t v_id_0 = cell_vtx[i];
    cs_lnum_t v_id_1 = cell_vtx[(i+1)%3];

    f_ids[i] = -1;

    for (cs_lnum_t f_id = 0; f_id < 3 && f_ids[i] == -1; f_id++) {
      if (f_id_quad[f_id] < 0)
        continue;
      cs_lnum_t v_id_2 = -1;
      cs_lnum_t v_id_3 = -1;
      for (cs_lnum_t j = 0; j < 4; j++) {
        if (vtx_quad[f_id][j] == v_id_0) {
          if (vtx_quad[f_id][((j+1) % 4)] == v_id_1) {
            v_id_2 = vtx_quad[f_id][((j+2) % 4)];
            v_id_3 = vtx_quad[f_id][((j+3) % 4)];
            if (i == 0) {
              cell_vtx[3] = v_id_3;
              cell_vtx[4] = v_id_2;
            }
            else if (i == 1) {
              if (cell_vtx[4] != v_id_3)
                matched = false; /* error case */
              cell_vtx[5] = v_id_2;
            }
            else if (i == 2) {
              if (cell_vtx[3] != v_id_2 || cell_vtx[5] != v_id_3)
                matched = false; /* error case */
            }
            f_ids[i] = f_id_quad[f_id];
            f_sgn[i] = f_sgn_quad[f_id];
            c2f2v_start[i] = j;
            f_id_quad[f_id] = -1;
            break;
          }
        }
      }
    }

    if (f_ids[i] < 0)
      matched = false;
  }

  /* Now match "top" triangle */

  f_ids[4] = -1;

  for (cs_lnum_t j = 0; j < 3; j++) {

    if (   vtx_tria[1][j]       == cell_vtx[3]
        && vtx_tria[1][(j+1)%3] == cell_vtx[4]
        && vtx_tria[1][(j+2)%3] == cell_vtx[5]) {
      f_ids[4] = f_id_tria[1];
      f_sgn[4] = f_sgn_tria[1];
      c2f2v_start[4] = j;
      f_id_tria[1] = -1;
    }

  }

  if (f_ids[4] < 0)
    matched = false;

  return matched;
}

/*----------------------------------------------------------------------------
 * Hexahedron reconstruction
 *
 * parameters:
 *   vtx_quad    <-- quadrangle faces connectivity
 *   f_id_quad   <-> matching (quadrangle) face ids, in initial order
 *                   set to -1 when matched
 *   f_sgn_quad  <-- matching quadrangle face orientation relative to cell
 *   f_ids       --> matching face ids, in template order
 *   f_sgn       --> matching face orientation relative to cell
 *   c2f2v_start --> matching initial vertex ids
 *
 * return:
 *   true if reconstruction successful, false otherwise
 *----------------------------------------------------------------------------*/

static bool
_match_hexa(cs_lnum_t  vtx_quad[6][4],
            cs_lnum_t  f_id_quad[6],
            int        f_sgn_quad[4],
            cs_lnum_t  f_ids[6],
            int        f_sgn[6],
            cs_lnum_t  c2f2v_start[6])
{
  bool matched = true;

  /*
    Base : vertices of the first face; we take the vertices in opposite
    order ("bottom" face numbering of the hexahedron with outward
    pointing normal).

                                                  Vertices
    *     7 x-------x 6                  Face 0:  0 3 2 1
    *      /|      /|                    Face 1:  0 1 5 4
    *     / |     / |                    Face 2:  1 2 6 5
    *  4 x-------x5 |                    Face 3:  2 3 7 6
    *    | 3x----|--x 2                  Face 4:  3 0 4 7
    *    | /     | /                     Face 5:  4 5 6 7
    *    |/      |/
    *  0 x-------x 1
  */

  cs_lnum_t  cell_vtx[8];

  cell_vtx[0] = vtx_quad[0][0];
  cell_vtx[1] = vtx_quad[0][3];
  cell_vtx[2] = vtx_quad[0][2];
  cell_vtx[3] = vtx_quad[0][1];

  f_ids[0] = f_id_quad[0];
  f_sgn[0] = f_sgn_quad[0];
  c2f2v_start[0] = 0;
  f_id_quad[0] = -1;

  /*
    We have found 4 out of 8 vertices; 4 other quadrangles should share
    one edge of the base quadrangle.
  */

  for (cs_lnum_t i = 1; i < 5; i++) {

    cs_lnum_t v_id_0 = cell_vtx[i-1];
    cs_lnum_t v_id_1 = cell_vtx[i%4];

    f_ids[i] = -1;

    for (cs_lnum_t f_id = 1; f_id < 6 && f_ids[i] == -1; f_id++) {
      if (f_id_quad[f_id] < 0)
        continue;
      cs_lnum_t v_id_2 = -1;
      cs_lnum_t v_id_3 = -1;
      for (cs_lnum_t j = 0; j < 4; j++) {
        if (vtx_quad[f_id][j] == v_id_0) {
          if (vtx_quad[f_id][((j+1) % 4)] == v_id_1) {
            v_id_2 = vtx_quad[f_id][((j+2) % 4)];
            v_id_3 = vtx_quad[f_id][((j+3) % 4)];
            if (i == 1) {
              cell_vtx[4] = v_id_3;
              cell_vtx[5] = v_id_2;
            }
            else if (i < 4) {  /* i = 2, 3 */
              if (cell_vtx[i+3] != v_id_3)
                matched = false; /* error case */
              cell_vtx[i+4] = v_id_2;
            }
            else { /* i = 4 */
              if (cell_vtx[4] != v_id_2 || cell_vtx[7] != v_id_3)
                matched = false; /* error case */
            }
            f_ids[i] = f_id_quad[f_id];
            f_sgn[i] = f_sgn_quad[f_id];
            c2f2v_start[i] = j;
            f_id_quad[f_id] = -1;
            break;
          }
        }
      }
    }

    if (f_ids[i] < 0)
      matched = false;

  }

  /* Now match "top" quadrangle */

  f_ids[5] = -1;

  for (cs_lnum_t i = 1; i < 6; i++) {
    if (f_id_quad[i] > -1) {

      for (cs_lnum_t j = 0; j < 4; j++) {

        if (   vtx_quad[i][j]       == cell_vtx[4]
            && vtx_quad[i][(j+1)%4] == cell_vtx[5]
            && vtx_quad[i][(j+2)%4] == cell_vtx[6]
            && vtx_quad[i][(j+3)%4] == cell_vtx[7]) {
          f_ids[5] = f_id_quad[i];
          f_sgn[5] = f_sgn_quad[i];
          f_id_quad[i] = -1;
          c2f2v_start[5] = j;
        }

      }

    }
  }

  if (f_ids[5] < 0)
    matched = false;

  return matched;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Determination of a given cell's refinement type.
 *
 * The cells->faces adjacency is updated as faces are locally renumbered
 * in canonical order for known cell types.
 *
 * \param[in]       m            pointer to mesh structure
 * \param[in]       cell_id      cell id (0 to n-1)
 * \param[in]       conforming   force conforming option ?
 * \param[in, out]  c2f          cells->faces adjacency (boundary faces first)
 * \param[in, out]  c2f2v_start  start id for adjacent face vertices
 *                               definitions (initialized to 0)
 * \param[in, out]  c_r_flag     cell refinement type
 * \param[in]       f_r_flag     face refinement type
 */
/*----------------------------------------------------------------------------*/

inline static void
_cell_r_type(const cs_mesh_t              *m,
             cs_lnum_t                     cell_id,
             bool                          conforming,
             cs_adjacency_t               *c2f,
             cs_lnum_t                     c2f2v_start[],
             cs_mesh_refine_type_t         c_r_flag[],
             const cs_mesh_refine_type_t   f_r_flag[])
{
  const cs_lnum_t s_id = c2f->idx[cell_id];
  const cs_lnum_t e_id = c2f->idx[cell_id+1];

  const cs_lnum_t n_cell_faces = e_id - s_id;

  /* Cells with more than 6 faces handled as general polyhedra.
     TODO check face refinement levels and re-associated faces
     of highest levels if possible to see if a pattern of cells with some
     faces already subdivided can be identified. */

  if (n_cell_faces > 6 && c_r_flag[cell_id] == CS_REFINE_DEFAULT) {
    c_r_flag[cell_id] = CS_REFINE_POLYHEDRON;
    for (cs_lnum_t i = s_id; i < e_id; i++) {
      if (f_r_flag[c2f->ids[i]] == CS_REFINE_POLYGON)
        c_r_flag[cell_id] = CS_REFINE_POLYHEDRON_P;
    }
  }
  else if (conforming) {
    for (cs_lnum_t i = s_id; i < e_id; i++) {
      if (f_r_flag[c2f->ids[i]] > CS_REFINE_NONE) {
        c_r_flag[cell_id] = CS_REFINE_DEFAULT;
        break;
      }
    }
  }

  if (   c_r_flag[cell_id] == CS_REFINE_NONE
      || c_r_flag[cell_id] >= CS_REFINE_POLYHEDRON)
    return;

  /*
    First guess at element type, based on matching face refinement templates;
    assume we have a "classical" element; for example, with 6 faces, we
    probably have a hexahedron, though we could also have a tetrahedron
    with one of its edges and thus two of its faces split. We count vertices
    per face to check.
  */

  cs_mesh_refine_type_t  _f_r_flag[6]; /* local flags */

  if (n_cell_faces <= 6) {
    for (cs_lnum_t i = s_id; i < e_id; i++) {
      cs_lnum_t l_id = i - s_id;
      cs_lnum_t f_id = c2f->ids[i];
      _f_r_flag[l_id] = f_r_flag[f_id];
      if (_f_r_flag[l_id] == CS_REFINE_POLYGON_T) {
        c_r_flag[cell_id] = CS_REFINE_POLYHEDRON_P;
        return;
      }
    }
  }

  switch (n_cell_faces) {

  case 4: /* tetrahedron ? */
    {
      for (cs_lnum_t i = 1; i < 4; i++) {
        if (_f_r_flag[i] != _f_r_flag[0]) {
          c_r_flag[cell_id] = CS_REFINE_POLYHEDRON;
          return;
        }
      }
      if (_f_r_flag[0] == CS_REFINE_TRIA) {
        c_r_flag[cell_id] = CS_REFINE_TETRA;
      }
      else if (_f_r_flag[0] == CS_REFINE_TRIA_Q) {
        c_r_flag[cell_id] = CS_REFINE_TETRA_H;
      }
      else if (_f_r_flag[0] == CS_REFINE_POLYGON_Q) {
        c_r_flag[cell_id] = CS_REFINE_POLYHEDRON;
      }
    }
    break;

  case 5: /* prism or pyramid ? */
    {
      int n_tria = 0, n_quad = 0;
      for (cs_lnum_t i = 0; i < 5; i++) {
        if (_f_r_flag[i] == CS_REFINE_TRIA)
          n_tria += 1;
        else if (_f_r_flag[i] == CS_REFINE_QUAD)
          n_quad += 1;
        else {
          c_r_flag[cell_id] = CS_REFINE_POLYHEDRON;
          return;
        }
        if (n_tria == 2 && n_quad == 3) {
          c_r_flag[cell_id] = CS_REFINE_PRISM;
        }
        else if (n_tria == 4 && n_quad == 1) {
          c_r_flag[cell_id] = CS_REFINE_PYRAM;
        }
      }
    }
    break;

  case 6: /* hexahedron ? */
    {
      c_r_flag[cell_id] = CS_REFINE_HEXA;
      for (cs_lnum_t i = 0; i < 6; i++) {
        if (_f_r_flag[i] != CS_REFINE_QUAD) {
          c_r_flag[cell_id] = CS_REFINE_POLYHEDRON;
          return;
        }
      }
    }
    break;

  default:
    c_r_flag[cell_id] = CS_REFINE_POLYHEDRON;
    return;
  }

  /* Now build base nodal connectivity */

  cs_lnum_t  n_trias = 0, n_quads = 0;
  cs_lnum_t  f_id_tria[4], f_id_quad[6];
  int        f_sgn_tria[4], f_sgn_quad[6];
  cs_lnum_t  vtx_tria[4][3];
  cs_lnum_t  vtx_quad[6][4];

  for (cs_lnum_t i = 0; i < n_cell_faces; i++) {

    cs_lnum_t f_id = c2f->ids[s_id + i];
    cs_lnum_t f_vtx_s, f_vtx_e;
    const cs_lnum_t *f_vtx;

    if (f_id < m->n_b_faces) {
      f_vtx_s = m->b_face_vtx_idx[f_id];
      f_vtx_e = m->b_face_vtx_idx[f_id+1];
      f_vtx = m->b_face_vtx_lst;
    }
    else {
      f_vtx_s = m->i_face_vtx_idx[f_id - m->n_b_faces];
      f_vtx_e = m->i_face_vtx_idx[f_id - m->n_b_faces + 1];
      f_vtx = m->i_face_vtx_lst;
    }

    cs_lnum_t n_fv = f_vtx_e - f_vtx_s;
    short int sgn = c2f->sgn[s_id + i];

    if (n_fv == 3) {
      assert(n_trias < 4);

      if (sgn > 0) {
        for (cs_lnum_t j = 0; j < 3; j++)
          vtx_tria[n_trias][j] = f_vtx[f_vtx_s + j];
      }
      else {
        for (cs_lnum_t j = 0; j < 3; j++)
          vtx_tria[n_trias][j] = f_vtx[f_vtx_e - 1 - j];
      }
      f_id_tria[n_trias] = f_id;
      f_sgn_tria[n_trias] = sgn;

      n_trias += 1;
    }
    else if (n_fv == 4) {
      assert(n_quads < 6);

      if (sgn > 0) {
        for (cs_lnum_t j = 0; j < 4; j++)
          vtx_quad[n_quads][j] = f_vtx[f_vtx_s + j];
      }
      else {
        for (cs_lnum_t j = 0; j < 4; j++)
          vtx_quad[n_quads][j] = f_vtx[f_vtx_e - 1 - j];
      }
      f_id_quad[n_quads] = f_id;
      f_sgn_quad[n_quads] = sgn;

      n_quads += 1;
    }

  }

  /* Now build nodal representation */

  cs_lnum_t  _f_ids[6], _c2f2v_start[6];
  int        _f_sgn[6];
  bool matched = false;

  switch(n_quads) {

  case 0:
    if (n_trias == 4) {
      matched = _match_tetra(vtx_tria, f_id_tria, f_sgn_tria,
                             _f_ids, _f_sgn, _c2f2v_start);
    }
    break;

  case 1:
    if (n_trias == 4) {
      matched = _match_pyram(vtx_tria, vtx_quad,
                             f_id_tria, f_id_quad,
                             f_sgn_tria, f_sgn_quad,
                             _f_ids, _f_sgn, _c2f2v_start);
    }
    break;

  case 3:
    if (n_trias == 2) {
      matched = _match_prism(vtx_tria, vtx_quad,
                             f_id_tria, f_id_quad,
                             f_sgn_tria, f_sgn_quad,
                             _f_ids, _f_sgn, _c2f2v_start);
    }
    break;

  case 6:
    if (n_trias == 0) {
      matched = _match_hexa(vtx_quad, f_id_quad, f_sgn_quad,
                            _f_ids, _f_sgn, _c2f2v_start);
    }
    break;

  default:
    break;
  }

  if (matched) {

    assert(n_cell_faces < 7);

    for (cs_lnum_t i = 0; i < n_cell_faces; i++) {
      c2f->ids[s_id + i] = _f_ids[i];
      c2f->sgn[s_id + i] = _f_sgn[i];
      c2f2v_start[s_id + i] = _c2f2v_start[i];
    }

  }
  else
    c_r_flag[cell_id] = CS_REFINE_POLYHEDRON;

  assert(c_r_flag[cell_id] != CS_REFINE_DEFAULT);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Determination of cell refinement flags.
 *
 * This assumes a first pass to determine adjacent face refinement flags
 * and an initial cell refinement flag is done.
 *
 * The cells->faces adjacency is updated as faces are locally renumbered
 * in canonical order for known cell types.
 *
 * \param[in]       m            pointer to mesh structure
 * \param[in]       conforming   force conforming option ?
 * \param[in, out]  c2f          cells->faces adjacency
 * \param[out]      c2f2v_start  start id for adjacent face vertices
 *                               definitions (initialized to 0)
 * \param[in, out]  c_r_flag     cell refinement type
 * \param[in]       f_r_flag     face refinement type
 */
/*----------------------------------------------------------------------------*/

static void
_cell_r_types(const cs_mesh_t              *m,
              bool                          conforming,
              cs_adjacency_t               *c2f,
              cs_lnum_t                     c2f2v_start[],
              cs_mesh_refine_type_t         c_r_flag[],
              const cs_mesh_refine_type_t   f_r_flag[])
{
# pragma omp for schedule(dynamic, CS_CL_SIZE)
  for (cs_lnum_t c_id = 0; c_id < m->n_cells; c_id++) {

    cs_lnum_t e_id = c2f->idx[c_id+1];
    for (cs_lnum_t j = c2f->idx[c_id]; j < e_id; j++)
      c2f2v_start[j] = 0;

    if (conforming || c_r_flag[c_id] != CS_REFINE_NONE)
      _cell_r_type(m, c_id, conforming, c2f, c2f2v_start, c_r_flag, f_r_flag);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Count number of added faces and associated connectivity
 *        size for new cells.
 *
 * The cells->faces adjacency is updated as faces are locally renumbered
 * in canonical order for known cell types.
 *
 * For polyhedral refinement schemes, the computation of these sizes is more
 * complex, so an overestimate is used here, and a separate pass will be
 * done do subdivide those cells. To mark such cells, the maximum estimation
 * of the number of sub-cells is multiplied by -1.
 *
 * \param[in]       m               pointer to mesh structure
 * \param[in, out]  c2f             cells->faces adjacency (boundary faces
 *                                  first, interior face ids shifted)
 * \param[in]       n_b_f_ini       old number of boundary faces
 * \param[in]       c_r_flag        cell refinement type
 * \param[in]       b_face_o2n_idx  boundary subface index
 * \param[in]       i_face_o2n_idx  interior subface index
 * \param[out]      c_poly_s_idx    index in polyhedral division
 *                                  (0 for simple cells, 1 for poly)
 * \param[out]      c_poly_f_idx    index in polyhedral division
 *                                  (0 for simple cells, n_faces for poly)
 * \param[out]      c_n_s_cells     number of sub-cells
 * \param[out]      c_n_i_faces     number of added interior faces
 * \param[out]      c_i_faces_size  size of added interior faces connectivity
 */
/*----------------------------------------------------------------------------*/

static void
_new_cells_i_faces_count(const cs_mesh_t            *m,
                         cs_adjacency_t             *c2f,
                         cs_lnum_t                   n_b_f_ini,
                         const cs_mesh_refine_type_t c_r_flag[],
                         const cs_lnum_t             b_face_o2n_idx[],
                         const cs_lnum_t             i_face_o2n_idx[],
                         cs_lnum_t         *restrict c_poly_s_idx,
                         cs_lnum_t         *restrict c_n_s_cells,
                         cs_lnum_t         *restrict c_n_i_faces,
                         cs_lnum_t         *restrict c_i_faces_size)
{
  const cs_lnum_t n_cells = m->n_cells;

  /* Counts for standard subdivition types */

  cs_lnum_t n_sub_cells[CS_REFINE_N_TYPES];    /* sub-elements */
  cs_lnum_t n_type_faces[CS_REFINE_N_TYPES];   /* additional faces */
  cs_lnum_t n_type_size[CS_REFINE_N_TYPES];    /* added connectivity size */

  for (int i = 0; i < CS_REFINE_N_TYPES; i++) {
    n_sub_cells[i] = 0;
    n_type_faces[i] = 0;
    n_type_size[i] = 0;
  }

  n_sub_cells[CS_REFINE_NONE] = 1;

  n_sub_cells[CS_REFINE_TETRA] = 8;     /* tetrahedra */
  n_type_faces[CS_REFINE_TETRA] = 8;
  n_type_size[CS_REFINE_TETRA] = 8*3;

  n_sub_cells[CS_REFINE_TETRA_H] = 4;   /* tetrahedra to hexahedra */
  n_type_faces[CS_REFINE_TETRA_H] = 6;
  n_type_size[CS_REFINE_TETRA_H] = 6*4;

  n_sub_cells[CS_REFINE_PYRAM] = 6 + 4; /* pyramids + tetrahedra */
  n_type_faces[CS_REFINE_PYRAM] = 1 + 12;
  n_type_size[CS_REFINE_PYRAM] = 4 + 12*3;

  n_sub_cells[CS_REFINE_PRISM] = 8;     /* prisms */
  n_type_faces[CS_REFINE_PRISM] = 6 + 4;
  n_type_size[CS_REFINE_PRISM] = 6*4 + 4*3;

  n_sub_cells[CS_REFINE_HEXA] = 8;      /* hexahedra */
  n_type_faces[CS_REFINE_HEXA] = 12;
  n_type_size[CS_REFINE_HEXA] = 12*4;

  /* Now loop on faces to determine counts */

# pragma omp for schedule(dynamic, CS_CL_SIZE)
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    const cs_mesh_refine_type_t rf = c_r_flag[cell_id];

    c_poly_s_idx[cell_id] = 0;

    if (rf < CS_REFINE_POLYHEDRON) {
      c_n_s_cells[cell_id] = n_sub_cells[rf];
      c_n_i_faces[cell_id] = n_type_faces[rf];
      c_i_faces_size[cell_id] = n_type_size[rf];
    }
    else {

      cs_lnum_t n_s_cells = 0;
      cs_lnum_t n_i_faces = 0;

      const cs_lnum_t s_id = c2f->idx[cell_id];
      const cs_lnum_t e_id = c2f->idx[cell_id+1];

      const cs_lnum_t n_cell_faces = e_id - s_id;

      for (cs_lnum_t i = 0; i < n_cell_faces; i++) {
        cs_lnum_t f_id = c2f->ids[s_id + i];
        const cs_lnum_t *o2n_idx, *vtx_idx;
        if (f_id < n_b_f_ini) {
          o2n_idx = b_face_o2n_idx;
          vtx_idx = m->b_face_vtx_idx;
        }
        else {
          f_id -= n_b_f_ini;
          o2n_idx = i_face_o2n_idx;
          vtx_idx = m->i_face_vtx_idx;
        }
        cs_lnum_t s_f_idx = o2n_idx[f_id];
        cs_lnum_t e_f_idx = o2n_idx[f_id+1];

        n_s_cells += e_f_idx - s_f_idx;
        // sub-faces are contiguous, so sum of face sizes does not need loop
        n_i_faces += vtx_idx[e_f_idx] - vtx_idx[s_f_idx];
      }

      assert(n_i_faces%2 == 0);

      c_poly_s_idx[cell_id] = 1;
      c_n_s_cells[cell_id] = n_s_cells;
      c_n_i_faces[cell_id] = n_i_faces/2;         /* each face shared */
      c_i_faces_size[cell_id] = n_i_faces/2 * 3;  /* triangle faces */

    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate polyhedral cell refinement temporary connnectivity arrays.
 *
 * The cells->faces adjacency is updated as faces are locally renumbered
 * in canonical order for known cell types.
 *
 * For polyhedral refinement schemes, the computation of these sizes is more
 * complex, so an overestimate is used here, and a separate pass will be
 * done do subdivide those cells. To mark such cells, the maximum estimation
 * of the number of sub-cells is multiplied by -1.
 *
 * \param[in]   n_cells_ini      number of initial cells
 * \param[in]   c_poly_s_idx     index in polyhedral division
 * \param[in]   c_n_e_faces      max number of exterior faces per cell
 *                               (also max number of sub-cells)
 * \param[in]   c_n_i_faces      max number of added interior faces per cell
 * \param[in]   c_i_faces_size   max size of added interior faces connectivity
 * \param[out]  c_o_id           old cell id
 * \param[out]  e_f2c_pre_idx    temp. exterior faces to cell index
 * \param[out]  e_f2c_pre        temp. exterior faces to cell connect
 * \param[out]  i_f2c_pre_idx    temp. interior faces to cells index
 * \param[out]  i_f2c_pre        temp. interior faces to cells connect
 * \param[out]  i_f2v_pre_idx    temp interior faces to vertices index
 * \param[out]  i_f2v_pre        temp interior faces to vertices connect.
 */
/*----------------------------------------------------------------------------*/

static void
_new_cells_poly_pre_alloc(cs_lnum_t          n_cells_ini,
                          const cs_lnum_t   *restrict c_poly_s_idx,
                          const cs_lnum_t   *restrict c_n_e_faces,
                          const cs_lnum_t   *restrict c_n_i_faces,
                          const cs_lnum_t   *restrict c_i_faces_size,
                          cs_lnum_t         *&c_o_id,
                          cs_lnum_t         *&e_f2c_pre_idx,
                          int16_t           *&e_f2c_pre,
                          cs_lnum_t         *&i_f2c_pre_idx,
                          int16_t           *&i_f2c_pre,
                          cs_lnum_t         *&i_f2v_pre_idx,
                          cs_lnum_t         *&i_f2v_pre
)
{
  cs_lnum_t n_poly_r_cells = c_poly_s_idx[n_cells_ini];

  CS_MALLOC(c_o_id, n_poly_r_cells, cs_lnum_t);
  CS_MALLOC(e_f2c_pre_idx, n_poly_r_cells+1, cs_lnum_t);
  CS_MALLOC(i_f2c_pre_idx, n_poly_r_cells+1, cs_lnum_t);
  CS_MALLOC(i_f2v_pre_idx, n_poly_r_cells+1, cs_lnum_t);

  e_f2c_pre_idx[0] = 0;
  i_f2c_pre_idx[0] = 0;
  i_f2v_pre_idx[0] = 0;

  int n_threads = cs_parall_n_threads(n_cells_ini, CS_THR_MIN);

# pragma omp parallel for  num_threads(n_threads)
  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ini; cell_id++) {
    cs_lnum_t s_idx = c_poly_s_idx[cell_id];
    cs_lnum_t e_idx = c_poly_s_idx[cell_id+1];

    /* For temporary cell interior faces to vertices connectivity,
       each face's connectivity will be closed by a -1 value,
       so one extra value per face is needed (hence + c_n_i_faces[cell_id]). */

    if (e_idx > s_idx) {
      c_o_id[s_idx] = cell_id;
      e_f2c_pre_idx[s_idx+1] = c_n_e_faces[cell_id];
      i_f2c_pre_idx[s_idx+1] = c_n_i_faces[cell_id] + 1;
      i_f2v_pre_idx[s_idx+1] = c_i_faces_size[cell_id] + c_n_i_faces[cell_id];
    }
  }

  _counts_to_index(n_poly_r_cells, e_f2c_pre_idx);
  _counts_to_index(n_poly_r_cells, i_f2c_pre_idx);
  _counts_to_index(n_poly_r_cells, i_f2v_pre_idx);

  CS_MALLOC(e_f2c_pre, e_f2c_pre_idx[n_poly_r_cells], int16_t);
  CS_MALLOC(i_f2c_pre, i_f2c_pre_idx[n_poly_r_cells]*2, int16_t);
  CS_MALLOC(i_f2v_pre, i_f2v_pre_idx[n_poly_r_cells], cs_lnum_t);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief return edge id matching 2 vertices
 *
 * \param[in]   v0   first vertex id
 * \param[in]   v1   second vertex id
 * \param[in]   v2v  vertex adjacency
 *
 * \return edge if matching 2 vertices, or -1
 */
/*----------------------------------------------------------------------------*/

static inline cs_lnum_t
_v2v_edge_id(cs_lnum_t              v0,
             cs_lnum_t              v1,
             const cs_adjacency_t  *v2v)
{
  cs_lnum_t edge_id = -1;

  if (v1 < v0) {
    cs_lnum_t vt = v0;
    v0 = v1;
    v1 = vt;
  }

  cs_lnum_t s_id = v2v->idx[v0];
  cs_lnum_t e_id = v2v->idx[v0+1];

  for (cs_lnum_t j = s_id; j < e_id; j++) {
    if (v2v->ids[j] == v1) {
      edge_id = j;
      break;
    }
  }

  return edge_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Sync edges flag for parallelism and determine associated
 *        added vertices global numbers.
 *
 * \param[in]       m            pointer to mesh structure
 * \param[in]       v2v          vertex adjacency
 * \param[in, out]  e_v_flag     for each edge, flag (count) for added vertices
 * \param[out]      g_edges_num  global edges number, or nullptr
 *
 * \return: global number of edges
 */
/*----------------------------------------------------------------------------*/

static cs_gnum_t
_sync_edges_flag(const cs_mesh_t        *m,
                 const cs_adjacency_t   *v2v,
                 cs_lnum_t               e_v_flag[],
                 cs_gnum_t              *g_edges_num)
{
  const cs_lnum_t n_vertices = v2v->n_elts;
  const cs_lnum_t n_edges = v2v->idx[v2v->n_elts];

  cs_gnum_t n_g_edges = n_edges;

  /* Build global edge numbering and edges interface */

  cs_gnum_t *g_e_vtx;
  CS_MALLOC(g_e_vtx, n_edges*2, cs_gnum_t);

  cs_lnum_t edge_id = 0;

  for (cs_lnum_t i = 0; i < n_vertices; i++) {
    cs_gnum_t g_v0 = m->global_vtx_num[i];
    cs_lnum_t e_id = v2v->idx[i+1];
    for (cs_lnum_t j = v2v->idx[i]; j < e_id; j++) {
      cs_gnum_t g_v1 = m->global_vtx_num[v2v->ids[j]];
      if (g_v0 < g_v1) {
        g_e_vtx[edge_id*2]   = g_v0;
        g_e_vtx[edge_id*2+1] = g_v1;
      }
      else {
        g_e_vtx[edge_id*2]   = g_v1;
        g_e_vtx[edge_id*2+1] = g_v0;
      }
      edge_id++;
    }
  }

  fvm_io_num_t *edge_io_num
    = fvm_io_num_create_from_adj_s(nullptr, g_e_vtx, n_edges, 2);

  CS_FREE(g_e_vtx);

  if (cs_glob_n_ranks > 1 || g_edges_num != nullptr) {
    n_g_edges = fvm_io_num_get_global_count(edge_io_num);
    const cs_gnum_t *_g_num =  fvm_io_num_get_global_num(edge_io_num);
    for (cs_lnum_t i = 0; i < n_edges; i++)
      g_edges_num[i] = _g_num[i];
    /* Rebuild as shared to free a bit of memory */
    edge_io_num = fvm_io_num_destroy(edge_io_num);
    edge_io_num = fvm_io_num_create_shared(g_edges_num, n_g_edges, n_edges);
  }

  cs_interface_set_t *e_if
    = cs_interface_set_create(n_edges,
                              nullptr,
                              fvm_io_num_get_global_num(edge_io_num),
                              nullptr, 0, nullptr, nullptr, nullptr);

  /* Synchronize added vertex counts */

  cs_interface_set_max(e_if, n_edges, 1, true, CS_LNUM_TYPE, e_v_flag);

  cs_interface_set_destroy(&e_if);

  edge_io_num = fvm_io_num_destroy(edge_io_num);

  return n_g_edges;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Synchronize interior faces flag for parallelism.
 *
 * \param[in]       m         pointer to mesh structure
 * \param[in, out]  f_v_flag  for each face, flag (count) for added vertices
 * \param[in, out]  f_r_flag  face refinement type flag
 *
 * \return: global number of added vertices
 */
/*----------------------------------------------------------------------------*/

static void
_sync_i_faces_flag(const cs_mesh_t        *m,
                   cs_lnum_t               f_v_flag[],
                   cs_mesh_refine_type_t   f_r_flag[])
{
  /* Build global interior faces interface */

  cs_interface_set_t *f_if
    = cs_interface_set_create(m->n_i_faces,
                              nullptr,
                              m->global_i_face_num,
                              nullptr, 0, nullptr, nullptr, nullptr);

  /* Synchronize added vertex counts */

  cs_interface_set_max(f_if, m->n_i_faces, 1, true, CS_LNUM_TYPE, f_v_flag);

  /* Synchronize refinement templates;
     NOTE: for anisotropic refinement, using the "max" enum
     value would not be sufficient, unless a spcific enum numbering
     is used */

  cs_interface_set_max(f_if, m->n_i_faces, 1, true, CS_INT_TYPE, f_r_flag);

  cs_interface_set_destroy(&f_if);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Flag faces and associated edges for subdivision.
 *
 * Triangle and quadrangle edges are subdivided, but not those of general
 * polygons (unless required by other faces and cells), as the chosen
 * subdivision scheme is different.
 *
 * Edges already refined at the current level are not refined.
 * If necessary, the list of vertices is rotated so that it
 * starts with a vertex of the lowest refinement level.
 *
 * TODO: for polygonal faces, check for angle between edges, and do not
 *       subdivide edges which are nearly aligned, as the vertex between
 *       the two could be considered to already be the result of a
 *       subdivision; this should allow improving mesh quality as we refine.
 *
 * \param[in]       f_id          id of this face
 * \param[in]       n_fv          number of vertices for this face
 * \param[in]       check_convex  check if faces are convex ?
 * \param[in]       f_r_level     initial face refinement level
 * \param[in, out]  f2v_lst       face vertices list
 * \param[in]       v2v           vertex adjacency
 * \param[in]       vtx_coords    vertex coordinates
 * \param[in]       vtx_r_gen     vertices refinement generation
 * \param[in, out]  e_v_idx       for each edge, start index of added vertices
 * \param[in, out]  f_v_idx       for each face, start index of added vertices
 * \param[in, out]  f_r_flag      face refinement type flag
 */
/*----------------------------------------------------------------------------*/

static void
_flag_faces_and_edges(cs_lnum_t               f_id,
                      cs_lnum_t               n_fv,
                      bool                    check_convex,
                      char                    f_r_level,
                      cs_lnum_t               f2v_lst[],
                      const cs_adjacency_t   *v2v,
                      const cs_coord_t        vtx_coords[][3],
                      const char              vtx_r_gen[],
                      cs_lnum_t               e_v_idx[],
                      cs_lnum_t               f_v_idx[],
                      cs_mesh_refine_type_t   f_r_flag[])
{
  if (f_r_flag[f_id] == CS_REFINE_NONE)
    return;

  cs_mesh_refine_type_t _f_flag = f_r_flag[f_id];

  char _vtx_r_gen_f[64];
  char *vtx_r_gen_f = _vtx_r_gen_f;
  if (n_fv > 64) {
    CS_MALLOC(vtx_r_gen_f, n_fv, char);
  }

  int last_corner_idx = 0;
  cs_lnum_t v0 = f2v_lst[0];
  char r_lv0 = vtx_r_gen[v0];

  /* Copy to local array for cheaper access, as multiple
     passes may be needed */

  for (cs_lnum_t i = 0; i < n_fv; i++) {
    cs_lnum_t j = f2v_lst[i];
    char r_lv = vtx_r_gen[j];
    vtx_r_gen_f[i] = r_lv;
    if (r_lv < r_lv0) {
      last_corner_idx = i;
      r_lv0 = r_lv;
    }
  }

  /* Rotate list if needed */

  if (last_corner_idx > 0) {
    cs_lnum_t _f2v_lst_o[64];
    cs_lnum_t *f2v_lst_o = _f2v_lst_o;
    if (n_fv > 64)
      CS_MALLOC(f2v_lst_o, n_fv, cs_lnum_t);
    memcpy(f2v_lst_o, f2v_lst, n_fv*sizeof(cs_lnum_t));
    for (cs_lnum_t i = 0; i < n_fv; i++) {
      cs_lnum_t j = f2v_lst_o[(i+last_corner_idx)%n_fv];
      f2v_lst[i] = j;
    }
    if (f2v_lst_o != _f2v_lst_o)
      CS_FREE(f2v_lst_o);
  }

  /* Flag edges */

  cs_lnum_t n_corner = 0, n_mid = 0, n_fv_c = n_fv;

  v0 = f2v_lst[0];
  r_lv0 = vtx_r_gen[v0];

  last_corner_idx = 0;  // Last corner

  for (cs_lnum_t i = 0; i < n_fv; i++) {
    int idx1 = (i+1)%n_fv, idx2 = (i+2)%n_fv;
    cs_lnum_t v1 = f2v_lst[idx1];
    char r_lv1 = vtx_r_gen_f[idx1], r_lv2 = vtx_r_gen_f[idx2];

    if (r_lv1 <= f_r_level) {  // Corner
      n_corner += 1;
      if (last_corner_idx == i) {
        cs_lnum_t vr, vc;
        if (v0 < v1) {
          vr = v0; vc = v1;
        }
        else {
          vr = v1; vc = v0;
        }
        cs_lnum_t edge_id = _v2v_edge_id(vr, vc, v2v);
        assert(edge_id > -1);
        e_v_idx[edge_id + 1] = 1;
        n_fv_c += 1;
      }
      last_corner_idx = idx1;
    }

    else {
      if (r_lv1 == f_r_level + 1) // Pre-existing mid-point
        n_mid += 1;
    }

    // Prepare for next vertex.
    v0 = v1;
    r_lv0 = r_lv1;
  }

  if (vtx_r_gen_f != _vtx_r_gen_f)
    CS_FREE(vtx_r_gen_f);

  /* Determine specific template */

  n_mid += (n_fv_c - n_fv);  // Add new mid-points to old ones

  if (_f_flag == CS_REFINE_DEFAULT) {
    _f_flag = CS_REFINE_POLYGON_Q;
    if (n_corner == 3 && n_fv_c == 6) {
      _f_flag = _refine_tria_type;
    }
    else if (n_corner == 4 && n_fv_c == 8) {
      _f_flag = CS_REFINE_QUAD;
    }
    if (_f_flag == CS_REFINE_POLYGON_Q) {
      if (n_mid != n_corner) {
        _f_flag = CS_REFINE_POLYGON;
      }
    }
  }

  if (check_convex && n_fv >= 4) {
    if (_polygon_is_nearly_convex(n_fv, f2v_lst, vtx_coords) == false) {
      _f_flag = CS_REFINE_POLYGON;
    }
  }

  /* Flag element centers if required */

  if (_f_flag != CS_REFINE_TRIA && _f_flag != CS_REFINE_POLYGON)
    f_v_idx[f_id+1] = 1;

  f_r_flag[f_id] = _f_flag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define vertices to add at refined edges and faces
 *
 * \param[in]       m             mesh
 * \param[in]       check_convex  check if faces are convex ?
 * \param[in]       v2v           vertex->vertex adjacency
 * \param[in]       vtx_r_gen     vertex refinement generation
 * \param[in]       c_r_level     cell refinement level
 * \param[in]       c_r_flag      cell refinement flag
 * \param[in, out]  f_r_level     face refinement level
 * \param[in, out]  f_r_flag      face refinement flag
 * \param[in, out]  e_v_idx       for each edge, start index of added vertices
 * \param[in, out]  f_v_idx       for each face, start index of added vertices
 * \param[in, out]  g_edges_num   for each edge, start global number of
 *                                added vertices (1 to n, to shift) or nullptr
 * \param[out]      n_add_vtx     number of added vertices for edges and faces
 *
 * \return global number of edges
 */
/*----------------------------------------------------------------------------*/

static cs_gnum_t
_new_edge_and_face_vertex_ids(cs_mesh_t                    *m,
                              bool                          check_convex,
                              const cs_adjacency_t         *v2v,
                              const char                    vtx_r_gen[],
                              const char                    c_r_level[],
                              const cs_mesh_refine_type_t   c_r_flag[],
                              char                          f_r_level[],
                              cs_mesh_refine_type_t         f_r_flag[],
                              cs_lnum_t                     e_v_idx[],
                              cs_lnum_t                     f_v_idx[],
                              cs_gnum_t                    *g_edges_num,
                              cs_lnum_t                     n_add_vtx[2])
{
  cs_gnum_t n_g_edges = 0;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells = m->b_face_cells;

  const cs_real_3_t *vtx_coords = (const cs_real_3_t *)m->vtx_coord;

  /* Loop on boundary faces */

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    bool subdivide = false;

    cs_lnum_t c_id = b_face_cells[f_id];
    if (c_id < n_cells) {
      if (c_r_flag[c_id] == CS_REFINE_DEFAULT) {
        subdivide = true;
        f_r_flag[f_id] = CS_REFINE_DEFAULT;
      }
      f_r_level[f_id] = c_r_level[c_id];
    }

    if (subdivide) {
      cs_lnum_t n_fv = m->b_face_vtx_idx[f_id+1] - m->b_face_vtx_idx[f_id];
      _flag_faces_and_edges(f_id,
                            n_fv,
                            check_convex,
                            f_r_level[f_id],
                            m->b_face_vtx_lst + m->b_face_vtx_idx[f_id],
                            v2v,
                            vtx_coords,
                            vtx_r_gen,
                            e_v_idx,
                            f_v_idx,
                            f_r_flag);
    }

  }

  /* Loop on interior faces */

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    bool subdivide = false;

    char f_c_r_level_ini[2], f_c_r_level_next[2];
    for (cs_lnum_t i = 0; i < 2; i++) {
      cs_lnum_t c_id = i_face_cells[f_id][i];
      f_c_r_level_ini[i] = c_r_level[c_id];
      f_c_r_level_next[i] = f_c_r_level_ini[i];
      if (c_r_flag[c_id] > CS_REFINE_NONE)
        f_c_r_level_next[i] += 1;
    }

    /* If a face is shared with a cell with higher refinement level,
       which is not further refined here, no need to subdivide it a
       second time. */

    char f_r_level_ini = cs::max(f_c_r_level_ini[0], f_c_r_level_ini[1]);
    char f_r_level_next = cs::max(f_c_r_level_next[0], f_c_r_level_next[1]);

    f_r_level[n_b_faces + f_id] = f_r_level_ini;
    if (f_r_level_next > f_r_level_ini) {
      subdivide = true;
      f_r_flag[n_b_faces + f_id] = CS_REFINE_DEFAULT;
    }

    if (subdivide) {
      cs_lnum_t n_fv = m->i_face_vtx_idx[f_id+1] - m->i_face_vtx_idx[f_id];
      _flag_faces_and_edges(f_id,
                            n_fv,
                            check_convex,
                            f_r_level_ini,
                            m->i_face_vtx_lst + m->i_face_vtx_idx[f_id],
                            v2v,
                            vtx_coords,
                            vtx_r_gen,
                            e_v_idx,
                            f_v_idx + m->n_b_faces,
                            f_r_flag + m->n_b_faces);
    }

  }

  cs_lnum_t n_edges = v2v->idx[v2v->n_elts];

  /* Parallel synchronization */

  if (cs_glob_n_ranks > 1) {
    n_g_edges = _sync_edges_flag(m, v2v, e_v_idx+1, g_edges_num);
    _sync_i_faces_flag(m,
                       f_v_idx + m->n_b_faces + 1,
                       f_r_flag + m->n_b_faces);
  }
  else
    n_g_edges = n_edges;

  cs_lnum_t n_base_vertices = m->n_vertices;

  /* Transform counts to index */

  for (cs_lnum_t i = 0; i < n_edges; i++) {
    e_v_idx[i+1] += e_v_idx[i];
    e_v_idx[i] += n_base_vertices;
  }
  e_v_idx[n_edges] += n_base_vertices;

  n_add_vtx[0] = e_v_idx[n_edges] - e_v_idx[0];
  n_base_vertices += n_add_vtx[0];

  /* Additional vertices at face centers */

  const cs_lnum_t n_faces = m->n_b_faces + m->n_i_faces;

  for (cs_lnum_t i = 0; i < n_faces; i++) {
    f_v_idx[i+1] += f_v_idx[i];
    f_v_idx[i] += n_base_vertices;
  }
  f_v_idx[n_faces] += n_base_vertices;

  n_add_vtx[1] = f_v_idx[n_faces] - f_v_idx[0];

  return n_g_edges;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define vertices to add inside refined cells.
 *
 * The number of added vertices for edges and faces must have been determined
 * first.
 *
 * \param[in]       m          mesh
 * \param[in]       c_r_flag   cell refinement flag
 * \param[in, out]  c_v_idx    for each cell, start index of added vertices
 * \param[in, out]  n_add_vtx  number of added vertices for edges, faces,
 *                             and cells
 */
/*----------------------------------------------------------------------------*/

static void
_new_cell_vertex_ids(cs_mesh_t                    *m,
                     const cs_mesh_refine_type_t   c_r_flag[],
                     cs_lnum_t                     c_v_idx[],
                     cs_lnum_t                     n_add_vtx[3])
{
  const cs_lnum_t n_cells = m->n_cells;

  /* Additional vertices at cell centers */

  c_v_idx[0] = 0;

  for (cs_lnum_t c_id = 0; c_id < m->n_cells; c_id++) {

    if (   c_r_flag[c_id] == CS_REFINE_HEXA
        || c_r_flag[c_id] == CS_REFINE_POLYHEDRON
        || c_r_flag[c_id] == CS_REFINE_POLYHEDRON_P
        || c_r_flag[c_id] == CS_REFINE_TETRA_H)
      c_v_idx[c_id+1] = 1;
    else
      c_v_idx[c_id+1] = 0;
  }

  /* Transfer to index */

  cs_lnum_t n_base_vertices = m->n_vertices + n_add_vtx[0] + n_add_vtx[1];

  for (cs_lnum_t i = 0; i < n_cells; i++) {
    c_v_idx[i+1] += c_v_idx[i];
    c_v_idx[i] += n_base_vertices;
  }
  c_v_idx[n_cells] += n_base_vertices;

  n_add_vtx[2] = c_v_idx[n_cells] - c_v_idx[0];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build new vertices on selected edges.
 *
 * These vertices are appended at the end of the initial vertex definitions.
 * The coordinates and numbering arrays should be resized before calling
 * this function (to allow for vertices inserted on edges, faces, and
 * cells with a single resize).
 *
 * \param[in]  m          mesh
 * \param[in]  v2v        vertex->vertex adjacency
 * \param[in]  n_e_vtx    local number of vertices added on edges
 * \param[in]  e_v_idx    for each edge, start index of added vertices
 * \param[in]  g_e_v_num  for each edge, start global number of
 *                        added vertices (1 to n, to shift) or nullptr
 *
 * \return  local number of vertices added on mid edges
 */
/*----------------------------------------------------------------------------*/

static void
_build_edge_vertices(cs_mesh_t             *m,
                     const cs_adjacency_t  *v2v,
                     const cs_lnum_t        n_e_vtx,
                     const cs_lnum_t        e_v_idx[],
                     const cs_gnum_t       *g_e_v_num)
{
  const cs_lnum_t n_edges = v2v->idx[v2v->n_elts];

  if (m->global_vtx_num != nullptr && g_e_v_num != nullptr) {

    cs_gnum_t gnum_shift = m->n_g_vertices - 1;

    for (cs_lnum_t id0 = 0; id0 < v2v->n_elts; id0++) {
      for (cs_lnum_t j = v2v->idx[id0]; j < v2v->idx[id0+1]; j++) {
        cs_lnum_t s_id = e_v_idx[j];
        cs_lnum_t e_id = e_v_idx[j+1];
        cs_lnum_t n_sub = e_id - s_id;
        cs_lnum_t id1 = v2v->ids[j];
        for (cs_lnum_t k = 0; k < n_sub; k++) {
          cs_lnum_t id2 = k + s_id;
          cs_real_t r = (1. + k) / (n_sub+1);
          for (cs_lnum_t l = 0; l < 3; l++)
            m->vtx_coord[(id2+k)*3 + l] = r * (  m->vtx_coord[id0*3 + l]
                                               + m->vtx_coord[id1*3 + l]);
          m->global_vtx_num[id2+k] = gnum_shift + g_e_v_num[j] + (cs_gnum_t)k;
          m->vtx_r_gen[id2+k] = cs::max(m->vtx_r_gen[id0],
                                        m->vtx_r_gen[id1]) + 1;
        }
      }
    }

  }
  else {

    for (cs_lnum_t id0 = 0; id0 < v2v->n_elts; id0++) {
      for (cs_lnum_t j = v2v->idx[id0]; j < v2v->idx[id0+1]; j++) {
        cs_lnum_t s_id = e_v_idx[j];
        cs_lnum_t e_id = e_v_idx[j+1];
        cs_lnum_t n_sub = e_id - s_id;
        cs_lnum_t id1 = v2v->ids[j];
        for (cs_lnum_t k = 0; k < n_sub; k++) {
          cs_lnum_t id2 = k + s_id;
          cs_real_t r = (1. + k) / (n_sub+1);
          for (cs_lnum_t l = 0; l < 3; l++)
            m->vtx_coord[(id2+k)*3 + l] = r * (  m->vtx_coord[id0*3 + l]
                                               + m->vtx_coord[id1*3 + l]);
          m->vtx_r_gen[id2+k] = cs::max(m->vtx_r_gen[id0],
                                        m->vtx_r_gen[id1]) + 1;
        }
      }
    }

    if (m->global_vtx_num != nullptr) {
      for (cs_lnum_t i = 0; i < n_e_vtx; i++)
        m->global_vtx_num[m->n_vertices + i] = m->n_vertices + i + 1;
    }

  }

  m->n_vertices = cs::max(m->n_vertices, e_v_idx[n_edges]);
}

/*----------------------------------------------------------------------------
 * Compute added vertex centers for quadrangle faces.
 *
 * parameters:
 *   vtx_ids     <-- face vertex ids
 *   vtx_coords  <-- vertex coordinates
 *   center      --> face center
 *----------------------------------------------------------------------------*/

static inline void
_quad_face_center(const cs_lnum_t    vtx_ids[],
                  const cs_real_t    vtx_coords[],
                  cs_real_t          center[3])
{
  /* mid-segment coordinates */

  cs_real_t ec[4][3];

  for (cs_lnum_t i = 0; i < 4; i++) {
    cs_lnum_t v0 = vtx_ids[i], v1 = vtx_ids[(i+1)%4];
    for (cs_lnum_t j = 0; j < 3; j++) {
      ec[i][j] = 0.5*(vtx_coords[v0*3 + j] + vtx_coords[v1*3 + j]);
    }
  }

  /* centers of segments joining opposite mid-edge coordinates */

  cs_real_t mc[2][3];

  for (cs_lnum_t i = 0; i < 2; i++) {
    cs_lnum_t v0 = i, v1 = i+2;
    for (cs_lnum_t j = 0; j < 3; j++) {
      mc[i][j] = 0.5*(ec[v0][j] + ec[v1][j]);
    }
  }

  /* center approximation (unweighted) */

  for (cs_lnum_t j = 0; j < 3; j++) {
    center[j] = 0.5*(mc[0][j] + mc[1][j]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Modify some face centers for the initial mesh.
 *
 * This focuses on better balance of sub-element shapes.
 *
 * \param[in]       n_faces    number of faces
 * \param[in]       f_r_flag   face refinement flag
 * \param[in]       f2v_idx    face-vertices index
 * \param[in]       f2v_lst    face-vertices connectivity
 * \param[in]       vtx_coord  vertex coordinates
 * \param[in, out]  f_center   face centers
 */
/*----------------------------------------------------------------------------*/

static void
_adjust_face_centers(cs_lnum_t                     n_faces,
                     const cs_mesh_refine_type_t   f_r_flag[],
                     const cs_lnum_t               f2v_idx[],
                     const cs_lnum_t               f2v_lst[],
                     const cs_real_t               vtx_coord[],
                     cs_real_t                     f_center[][3])
{
  /* Loop on faces */

# pragma omp for schedule(dynamic, CS_CL_SIZE)
  for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++) {

    if (f_r_flag[f_id] == CS_REFINE_QUAD) {
      cs_lnum_t n_vtx = f2v_idx[f_id+1] - f2v_idx[f_id];
      if (n_vtx == 4)
        _quad_face_center(f2v_lst + f2v_idx[f_id],
                          vtx_coord,
                          f_center[f_id]);
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute approximate tetrahedron center as the mean of its face
 *         centers weighted by the associated surfaces.
 *
 *           n-1
 *           Sum   G(Fi)
 *           i=0
 *  G(C) = -------------
 *               n
 *
 * \param[in]   m            pointer to mesh structure
 * \param[in]   n_b_f_ini    old number of boundary faces
 * \param[in]   f_ids        ids of tetrahedral faces for this tetrahedron
 * \param[out]  cell_cen     cell center
 */
/*----------------------------------------------------------------------------*/

static void
_cell_tetra_isocog(const cs_mesh_t  *m,
                   cs_lnum_t         n_b_f_ini,
                   const cs_lnum_t   f_ids[4],
                   cs_real_t         cell_cen[3])
{
  assert(cell_cen != nullptr);

  const cs_real_3_t *vtx_coord = (const cs_real_3_t *)m->vtx_coord;

  /* Initialization */

  cs_real_t one_third = 1./3.;

  for (cs_lnum_t i = 0; i < 3; i++)
    cell_cen[i] = 0.;

  /* Loop on faces */

  const cs_lnum_t *f_vtx;

  for (cs_lnum_t f_idx = 0; f_idx < 4; f_idx++) {

    cs_lnum_t f_id = f_ids[f_idx];

    if (f_id < n_b_f_ini) {
      cs_lnum_t s_id = m->b_face_vtx_idx[f_id];
      f_vtx = m->b_face_vtx_lst + s_id;
    }
    else {
      cs_lnum_t s_id = m->i_face_vtx_idx[f_id - n_b_f_ini];
      f_vtx = m->i_face_vtx_lst + s_id;
    }

    /* Computation of the area of the face */

    const cs_lnum_t v0 = f_vtx[0];
    const cs_lnum_t v1 = f_vtx[1];
    const cs_lnum_t v2 = f_vtx[2];

    for (cs_lnum_t i = 0; i < 3; i++)
      cell_cen[i] += one_third * (  vtx_coord[v0][i]
                                  + vtx_coord[v1][i]
                                  + vtx_coord[v2][i]);

  } /* End of loop on faces */

  for (cs_lnum_t i = 0; i < 3; i++)
    cell_cen[i] *= 0.25;
}

/*----------------------------------------------------------------------------
 * Compute cell and face centers adapted for refinement.
 *
 * This uses a slightly different algorithm for some element centers
 * than the algorithm used for mesh quantities, focusing on better balance
 * of sub-element shapes.
 *
 * The corresponding areas are allocated by this function, and it is the
 * caller's responsability to free them when they are no longer needed.
 *
 * \param[in]   m           mesh
 * \param[in]   f_r_flag    face refinement flag
 * \param[out]  cell_cen    pointer to the cell centers array
 * \param[out]  i_face_cen  pointer to the interior face centers array
 * \param[out]  b_face_cen  pointer to the boundary face centers array
 *----------------------------------------------------------------------------*/

static void
_element_centers(const cs_mesh_t              *m,
                 const cs_mesh_refine_type_t   f_r_flag[],
                 cs_real_3_t                  *cell_cen[],
                 cs_real_3_t                  *i_face_cen[],
                 cs_real_3_t                  *b_face_cen[])
{
  cs_lnum_t  n_cells_with_ghosts = m->n_cells_with_ghosts;

  cs_real_3_t  *_cell_cen = nullptr;
  CS_MALLOC(_cell_cen, n_cells_with_ghosts, cs_real_3_t);

  /* Modify some cell centers and compute surfaces;
     Note that this could be delayed to cell center creation,
     avoiding computations for non-refined elements.*/

  cs_real_3_t  *_i_face_cen = nullptr, *i_face_normal = nullptr;
  cs_real_3_t  *_b_face_cen = nullptr, *b_face_normal = nullptr;
  CS_MALLOC_HD(_i_face_cen, m->n_i_faces, cs_real_3_t, cs_alloc_mode);
  CS_MALLOC_HD(i_face_normal, m->n_i_faces, cs_real_3_t, cs_alloc_mode);
  CS_MALLOC_HD(_b_face_cen, m->n_b_faces, cs_real_3_t, cs_alloc_mode);
  CS_MALLOC_HD(b_face_normal, m->n_b_faces, cs_real_3_t, cs_alloc_mode);

  cs_mesh_quantities_compute_face_cog_sn
    (m->n_i_faces,
     reinterpret_cast<const cs_real_3_t *>(m->vtx_coord),
     m->i_face_vtx_idx,
     m->i_face_vtx_lst,
     _i_face_cen,
     i_face_normal);

  cs_mesh_quantities_compute_face_cog_sn
    (m->n_b_faces,
     reinterpret_cast<const cs_real_3_t *>(m->vtx_coord),
     m->b_face_vtx_idx,
     m->b_face_vtx_lst,
     _b_face_cen,
     b_face_normal);

  _adjust_face_centers(m->n_i_faces,
                       f_r_flag + m->n_b_faces,
                       m->i_face_vtx_idx,
                       m->i_face_vtx_lst,
                       m->vtx_coord,
                       _i_face_cen);

  _adjust_face_centers(m->n_b_faces,
                       f_r_flag,
                       m->b_face_vtx_idx,
                       m->b_face_vtx_lst,
                       m->vtx_coord,
                       _b_face_cen);

  cs_mesh_quantities_cell_faces_cog(m,
                                    i_face_normal,
                                    _i_face_cen,
                                    b_face_normal,
                                    _b_face_cen,
                                    _cell_cen);

  CS_FREE(b_face_normal);
  CS_FREE(i_face_normal);

  *cell_cen = (cs_real_3_t *)_cell_cen;
  *i_face_cen = (cs_real_3_t *)_i_face_cen;
  *b_face_cen = (cs_real_3_t *)_b_face_cen;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build new vertices on faces.
 *
 * These vertices are appended at the end of the initial vertex definitions.
 * The coordinates and numbering arrays should be resized before calling
 * this function (to allow for vertices inserted on edges, faces, and possibly
 * cells with a single resize).
 *
 * \param[in, out]  m          mesh
 * \param[in]       n_faces    number of faces
 * \param[in]       f_vtx_idx  face -> vertices index
 * \param[in]       f_vtx_lst  face -> vertices list
 * \param[in]       f_r_flag   face refinement flag
 * \param[in]       f_v_idx    for each face, start index of added vertices
 * \param[in]       f_center   face centers
 *
 * \return  local number of vertices added on faces
 */
/*----------------------------------------------------------------------------*/

static void
_build_face_vertices(cs_mesh_t                    *m,
                     cs_lnum_t                     n_faces,
                     const cs_lnum_t               f_vtx_idx[],
                     const cs_lnum_t               f_vtx_lst[],
                     const cs_mesh_refine_type_t   f_r_flag[],
                     const cs_lnum_t               f_v_idx[],
                     const cs_real_t               f_center[][3])
{
  /* Loop on faces */

# pragma omp for schedule(dynamic, CS_CL_SIZE)
  for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++) {

    cs_lnum_t n_add = f_v_idx[f_id+1] - f_v_idx[f_id];

    if (n_add == 1) {
      cs_lnum_t s_id = f_vtx_idx[f_id];
      cs_lnum_t e_id = f_vtx_idx[f_id+1];
      char r_level = _compute_face_r_level(e_id - s_id,
                                           f_vtx_lst + s_id,
                                           m->vtx_r_gen);

      if (f_r_flag[f_id] == CS_REFINE_QUAD && (e_id - s_id == 4)) {
        _quad_face_center(f_vtx_lst + s_id,
                          m->vtx_coord,
                          m->vtx_coord + (f_v_idx[f_id]*3));
      }
      else {
        for (cs_lnum_t i = 0; i < 3; i++)
          m->vtx_coord[(f_v_idx[f_id]*3) + i] = f_center[f_id][i];
      }
      m->vtx_r_gen[f_v_idx[f_id]] = r_level + 1;

    }
    else if (n_add > 1) {
      CS_UNUSED(f_r_flag);
      assert(0); /* Handle other cases when other templates are added */
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build new vertices inside cells.
 *
 * These vertices are appended at the end of the initial vertex definitions.
 * The coordinates and numbering arrays should be resized before calling
 * this function (to allow for vertices inserted on edges, faces, and possibly
 * cells with a single resize).
 *
 * \param[in, out]  m          mesh
 * \param[in]       c2f        cells->faces adjacency (boundary faces first)
 * \param[in]       n_b_f_ini  old number of boundary faces
 * \param[in]       n_cells    number of cells
 * \param[in]       c_r_flag   cell refinement flag
 * \param[in]       c_r_gen    cell refinement generation
 * \param[in]       c_v_idx    for each cell, start index of added vertices
 * \param[in]       c_center   cell centers
 *
 * \return  local number of vertices added on mid edges
 */
/*----------------------------------------------------------------------------*/

static void
_build_cell_vertices(cs_mesh_t                    *m,
                     const cs_adjacency_t         *c2f,
                     cs_lnum_t                     n_b_f_ini,
                     cs_lnum_t                     n_cells,
                     const cs_mesh_refine_type_t   c_r_flag[],
                     const char                    c_r_gen[],
                     const cs_lnum_t               c_v_idx[],
                     const cs_real_t               c_center[][3])
{
  /* Loop on cells */

# pragma omp for schedule(dynamic, CS_CL_SIZE)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    if (c_v_idx[c_id+1] - c_v_idx[c_id] > 0) {
      switch(c_r_flag[c_id]) {
      case CS_REFINE_HEXA:
      case CS_REFINE_PRISM:
      case CS_REFINE_POLYHEDRON:
      case CS_REFINE_POLYHEDRON_P:
        {
          const cs_lnum_t v_id = c_v_idx[c_id];
          for (cs_lnum_t i = 0; i < 3; i++)
            m->vtx_coord[v_id*3 + i] = c_center[c_id][i];
          m->vtx_r_gen[v_id] = c_r_gen[c_id] + 1;
        }
        break;
      case CS_REFINE_TETRA_H:
        {
          const cs_lnum_t s_id_c = c2f->idx[c_id];
          const cs_lnum_t v_id = c_v_idx[c_id];

          _cell_tetra_isocog(m, n_b_f_ini, c2f->ids+s_id_c,
                             m->vtx_coord + v_id*3);
          m->vtx_r_gen[v_id] = c_r_gen[c_id] + 1;
        }
        break;
      default:
        assert(0);
        break;
      }
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build global numbers for new vertices on edges, faces, or cells.
 *
 * These vertices are appended at the end of the initial vertex definitions.
 * The numbering arrays should be resized before calling this function
 * (to allow for vertices inserted on edges, faces, and possibly
 * cells with a single resize).
 *
 * Each call of this function updates the global number of vertices
 * member of the mesh structure.
 *
 * \param[in, out]  m          mesh
 * \param[in]       n_elts     number of parent elements
 * \param[in]       n_g_elts   global number of parent elements
 * \param[in]       elt_v_idx  for each element, start index of added vertices
 * \param[in]       g_elt_num  global number of each element
 */
/*----------------------------------------------------------------------------*/

static void
_build_vertices_gnum(cs_mesh_t       *m,
                     cs_lnum_t        n_elts,
                     cs_gnum_t        n_g_elts,
                     const cs_lnum_t  elt_v_idx[],
                     const cs_gnum_t  g_elt_num[])
{
  cs_gnum_t n_g_add_vtx = 0;

  /* Loop on elements */

  if (cs_glob_n_ranks == 1 && g_elt_num == nullptr) {

    if (m->global_vtx_num != nullptr) {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        for (cs_lnum_t j = elt_v_idx[i]; j < elt_v_idx[i+1]; j++)
          m->global_vtx_num[j] = j+1;
      }
    }

    n_g_add_vtx = elt_v_idx[n_elts] - m->n_g_vertices;

  }
  else {

    /* Build associated global numbering */

    fvm_io_num_t *elt_io_num
      = fvm_io_num_create_shared(g_elt_num, n_g_elts, n_elts);

    cs_lnum_t *n_sub;
    CS_MALLOC(n_sub, n_elts, cs_lnum_t);
    cs_lnum_t *restrict _n_sub = n_sub;
    for (cs_lnum_t i = 0; i < n_elts; i++)
      _n_sub[i] = elt_v_idx[i+1] - elt_v_idx[i];
    _n_sub = nullptr;

    fvm_io_num_t *vtx_io_num
      = fvm_io_num_create_from_sub(elt_io_num, n_sub);

    elt_io_num = fvm_io_num_destroy(elt_io_num);

    CS_FREE(n_sub);

    const cs_gnum_t *add_vtx_gnum = fvm_io_num_get_global_num(vtx_io_num);
    n_g_add_vtx = fvm_io_num_get_global_count(vtx_io_num);

    assert(   elt_v_idx[n_elts] - elt_v_idx[0]
           == fvm_io_num_get_local_count(vtx_io_num));

    if (m->global_vtx_num != nullptr) {
      cs_lnum_t k = 0;
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        for (cs_lnum_t j = elt_v_idx[i]; j < elt_v_idx[i+1]; j++, k++)
          m->global_vtx_num[j] = add_vtx_gnum[k] + m->n_g_vertices;
      }
    }

    vtx_io_num = fvm_io_num_destroy(vtx_io_num);

  }

  m->n_g_vertices += n_g_add_vtx;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Count new faces and face vertices for non-subdivided face
 *
 * \param[in] n_fv     number of vertices for this face
 * \param[in] f2v_lst  face vertices list
 * \param[in] v2v      vertex adjacency
 * \param[in] e_v_idx  for each edge, start index of added vertices
 *
 * \return new number of vertices for a given face
 */
/*----------------------------------------------------------------------------*/

static inline cs_lnum_t
_count_face_edge_vertices_new(cs_lnum_t              n_fv,
                              const cs_lnum_t        f2v_lst[],
                              const cs_adjacency_t  *v2v,
                              const cs_lnum_t        e_v_idx[])
{
  cs_lnum_t n_new = n_fv;

  for (cs_lnum_t i = 0; i < n_fv; i++) {
    cs_lnum_t v0 = f2v_lst[i];
    cs_lnum_t v1 = f2v_lst[(i+1)%n_fv];

    cs_lnum_t edge_id = _v2v_edge_id(v0, v1, v2v);
    n_new += (e_v_idx[edge_id+1] - e_v_idx[edge_id]);
  }

  return n_new;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Count new face vertices and edge mid-points.
 *
 * Whenever a edge is refined, the added vertex has a level higher than its
 * adjacent vertices, so we can identify "non-corner" vertices in this way.
 *
 * \param[in]   n_fv       number of vertices for this face
 * \param[in]   f2v_lst    face vertices list
 * \param[in]   v2v        vertex adjacency
 * \param[in]   vtx_r_gen  vertices refinement generation
 * \param[in]   e_v_idx    for each edge, start index of added vertices
 *
 * \param[out]  n_new      number of vertices in subdivided perimeter
 * \param[out]  n_corner   number of corner points in subdivided perimeter
 * \param[out]  n_mid      number of mid-points in subdivided perimeter
 */
/*----------------------------------------------------------------------------*/

static inline void
_count_face_edge_vertices_sub(int                    n_fv,
                              char                   f_r_level,
                              const cs_lnum_t        f2v_lst[],
                              const cs_adjacency_t  *v2v,
                              const char             vtx_r_gen[],
                              const cs_lnum_t        e_v_idx[],
                              cs_lnum_t             &n_new,
                              cs_lnum_t             &n_corner,
                              cs_lnum_t             &n_mid)
{
  n_new = n_fv, n_corner = 0, n_mid = 0;

  cs_lnum_t v0 = f2v_lst[0];
  char r_lv0 = vtx_r_gen[v0];

  /* Corners should have a refinement generation lower or equal
     to the initial face refinement level, so can be identified
     using this property. In case edges are subdivided multiple
     times, a mid-point should be the vertex with highest
     refinement generation between 2 corners.

     Note that no face vertex should have a lower refinement
     generation than the first one, as we have rotated its vertices
     to ensure this in a prior step (in _flag_faces_and_edges).
  */

  cs_lnum_t last_corner_idx = 0;

  for (cs_lnum_t i = 0; i < n_fv; i++) {
    cs_lnum_t v1 = f2v_lst[(i+1)%n_fv];
    char r_lv1 = vtx_r_gen[v1];

    cs_lnum_t edge_id = _v2v_edge_id(v0, v1, v2v);
    assert(edge_id > -1);

    cs_lnum_t n_ev = (e_v_idx[edge_id+1] - e_v_idx[edge_id]);
    n_new += n_ev;

    if (r_lv0 <= f_r_level) {
      n_corner += 1;
      if (r_lv1 <= f_r_level && n_ev > 0)
        n_mid += 1;
      last_corner_idx = i;
    }
    else {
      /* Actual midpoint might not directly follow corner,
         but there can only be one midpoint between 2 corners. */
      if (i == last_corner_idx+1)
        n_mid += 1;
    }

    r_lv0 = r_lv1; // for next loop iteration
    v0 = v1;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief compute new connectivity size for a given face.
 *
 * \param[in]       n_fv           number of vertices for this face
 * \param[in, out]  f_r_flag       face refinement type flag
 * \param[in]       f_r_level_o    old face refinement level
 * \param[in]       f2v_lst_o      old face vertices list
 * \param[in]       v2v            vertex adjacency
 * \param[in]       vtx_r_gen      vertices refinement generation
 * \param[in]       e_v_idx        for each edge, start index of added vertices
 * \param[out]      n_sub          number of sub-faces
 * \param[out]      connect_size   associated connectivity size
 */
/*----------------------------------------------------------------------------*/

static void
_subdivided_face_sizes(const cs_lnum_t          n_fv,
                       cs_mesh_refine_type_t   &f_r_flag,
                       const char               f_r_level_o,
                       const cs_lnum_t          f2v_lst_o[],
                       const cs_adjacency_t    *v2v,
                       const char               vtx_r_gen[],
                       const cs_lnum_t          e_v_idx[],
                       cs_lnum_t               *n_sub,
                       cs_lnum_t               *connect_size)
{
  if (f_r_flag < CS_REFINE_POLYGON) {
    cs_lnum_t n_fv_c = _count_face_edge_vertices_new(n_fv,
                                                     f2v_lst_o,
                                                     v2v,
                                                     e_v_idx);
    switch(f_r_flag) {
    case CS_REFINE_NONE:
      *n_sub = 1;
      *connect_size = n_fv_c;
      break;
    case CS_REFINE_TRIA:
      *n_sub = 4;
      *connect_size = n_fv_c + 6;
      if (n_fv_c != 6)
        f_r_flag = CS_REFINE_POLYGON_Q;
      break;
    case CS_REFINE_TRIA_Q:
      *n_sub = 3;
      *connect_size = n_fv_c + 6;
      if (n_fv_c != 6) {
        f_r_flag = CS_REFINE_POLYGON_Q;
        /* TODO: the setting above modifies the behavior for trianges, as this
         * is the equivalent of replacing CS_REFINE_TRIA with CS_REFINE_TRIA_Q.
         * For this reason, handling subdivided edges in the CS_REFINE_TRIA case
         * in _subdivide_face and in the CS_REFINE_TETRA, CS_REFINE_PYRAM, and
         * CS_REFINE_PRISME matching and subdivision functions should be done. */
      }
      break;
    case CS_REFINE_QUAD:
      *n_sub = 4;
      *connect_size = n_fv_c + 8;
      if (n_fv_c != 8)
        f_r_flag = CS_REFINE_POLYGON_Q;
      break;
    default:
      break;
    }
  }

  else {
    cs_lnum_t n_fv_tot, n_corner, n_mid;
    _count_face_edge_vertices_sub(n_fv,
                                  f_r_level_o,
                                  f2v_lst_o,
                                  v2v,
                                  vtx_r_gen,
                                  e_v_idx,
                                  n_fv_tot,
                                  n_corner,
                                  n_mid);
    switch(f_r_flag) {
    case CS_REFINE_POLYGON_T:
      *n_sub = n_fv_tot;
      *connect_size = n_fv_tot*3;
      break;
    case CS_REFINE_POLYGON_Q:
      *n_sub = n_corner;
      *connect_size = n_fv_tot + n_mid*2;
      break;
    case CS_REFINE_POLYGON:
      *n_sub = n_fv_tot-2;
      *connect_size = (n_fv_tot-2)*3;
      break;
    default:
      *n_sub = 1;
      *connect_size = n_fv_tot;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief return vertex id ad mid-edge for simple edge subdivision patterns
 *
 * \param[in]   v0       first vertex id
 * \param[in]   v1       second vertex id
 * \param[in]   v2v      vertex adjacency
 * \param[in]   e_v_idx  for each edge, start index of added vertices

 * \return mid-edge vertex id
 */
/*----------------------------------------------------------------------------*/

static inline cs_lnum_t
_v2v_mid_vtx_id(cs_lnum_t              v0,
                cs_lnum_t              v1,
                const cs_adjacency_t  *v2v,
                const cs_lnum_t        e_v_idx[])
{
  cs_lnum_t edge_id = _v2v_edge_id(v0, v1, v2v);

  assert(edge_id > -1);

  return e_v_idx[edge_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief prepare new connectivity for a given face.
 *
 * If the f2v_count_n and f2v_lst_n arguments are nullptr, only sub-face
 * and connectivity size counts are returned.
 *
 * \param[in]      f_id           id for this face
 * \param[in]      n_fv           number of vertices for this face
 * \param[in]      f_r_flag       face refinement type flag
 * \param[in]      f_r_level_o    old face refinement level
 * \param[in]      f2v_lst_o      old face vertices list
 * \param[in]      v2v            vertex adjacency
 * \param[in]      e_v_idx        for each edge, start index of added vertices
 * \param[in]      f_v_idx        for each face, start index of added vertices
 * \param[in]      vtx_coords     mesh vertex coordinates
 * \param[in]      vtx_r_gen      vertices refinement generation
 * \param[in,out]  t_state        triangulation state if needed
 * \param[in,out]  f2v_idx_n      new faces to vertices index, initialized for
 *                                first vertex
 * \param[out]     f2v_lst_n      new faces connectivity
 */
/*----------------------------------------------------------------------------*/

static void
_subdivide_face(cs_lnum_t                 f_id,
                cs_lnum_t                 n_fv,
                cs_mesh_refine_type_t     f_r_flag,
                char                      f_r_level_o,
                const cs_lnum_t           f2v_lst_o[],
                const cs_adjacency_t     *v2v,
                const cs_lnum_t           e_v_idx[],
                const cs_lnum_t           f_v_idx[],
                const cs_real_t           vtx_coords[],
                const char                vtx_r_gen[],
                fvm_triangulate_state_t  *t_state,
                cs_lnum_t                 f2v_idx_n[],
                cs_lnum_t                 f2v_lst_n[])
{
  cs_lnum_t n_fv_max = 32;
  cs_lnum_t _f2v_lst_c[32];
  cs_lnum_t *f2v_lst_c = _f2v_lst_c;

  cs_lnum_t n_fv_c = 0;

  /* Build list of vertices on base face contour */

  for (cs_lnum_t i = 0; i < n_fv; i++) {

    if (n_fv_c+2 >= n_fv_max) {
      CS_MALLOC(f2v_lst_c, n_fv*2, cs_lnum_t);
      memcpy(f2v_lst_c, _f2v_lst_c, n_fv_c*sizeof(cs_lnum_t));
    }

    cs_lnum_t v0 = f2v_lst_o[i];
    cs_lnum_t v1 = f2v_lst_o[(i+1)%n_fv];

    f2v_lst_c[n_fv_c++] = v0;

    cs_lnum_t vr, vc;
    if (v0 < v1) {
      vr = v0; vc = v1;
    }
    else {
      vr = v1; vc = v0;
    }
    cs_lnum_t edge_id = _v2v_edge_id(vr, vc, v2v);
    assert(edge_id > -1);

    cs_lnum_t v_idx = e_v_idx[edge_id];
    if (e_v_idx[edge_id+1] > v_idx) {
      f2v_lst_c[n_fv_c++] = v_idx;
    }

  }

  switch(f_r_flag) {

  case CS_REFINE_NONE:
    {
      f2v_idx_n[1] = f2v_idx_n[0] + n_fv_c;
      memcpy(f2v_lst_n, f2v_lst_c, n_fv_c*sizeof(cs_lnum_t));
    }
    break;

  case CS_REFINE_TRIA:
    {
      assert(n_fv_c == 6);

      f2v_idx_n[1] = f2v_idx_n[0] + 3;
      f2v_idx_n[2] = f2v_idx_n[0] + 6;
      f2v_idx_n[3] = f2v_idx_n[0] + 9;

      f2v_lst_n[0] = f2v_lst_c[0];
      f2v_lst_n[1] = f2v_lst_c[1];
      f2v_lst_n[2] = f2v_lst_c[5];

      f2v_lst_n[3] = f2v_lst_c[2];
      f2v_lst_n[4] = f2v_lst_c[3];
      f2v_lst_n[5] = f2v_lst_c[1];

      f2v_lst_n[6] = f2v_lst_c[4];
      f2v_lst_n[7] = f2v_lst_c[5];
      f2v_lst_n[8] = f2v_lst_c[3];

      f2v_lst_n[9]  = f2v_lst_c[1];
      f2v_lst_n[10] = f2v_lst_c[3];
      f2v_lst_n[11] = f2v_lst_c[5];
    }
    break;

  case CS_REFINE_TRIA_Q:
    {
      assert(n_fv_c == 6);

      f2v_idx_n[1] = f2v_idx_n[0] + 4;
      f2v_idx_n[2] = f2v_idx_n[0] + 8;

      f2v_lst_n[0] = f2v_lst_c[0];
      f2v_lst_n[1] = f2v_lst_c[1];
      f2v_lst_n[2] = f_v_idx[f_id];
      f2v_lst_n[3] = f2v_lst_c[5];

      f2v_lst_n[4] = f2v_lst_c[2];
      f2v_lst_n[5] = f2v_lst_c[3];
      f2v_lst_n[6] = f_v_idx[f_id];
      f2v_lst_n[7] = f2v_lst_c[1];

      f2v_lst_n[8]  = f2v_lst_c[4];
      f2v_lst_n[9]  = f2v_lst_c[5];
      f2v_lst_n[10] = f_v_idx[f_id];
      f2v_lst_n[11] = f2v_lst_c[3];
    }
    break;

  case CS_REFINE_QUAD:
    {
      assert(n_fv_c == 8);

      f2v_idx_n[1] = f2v_idx_n[0] + 4;
      f2v_idx_n[2] = f2v_idx_n[0] + 8;
      f2v_idx_n[3] = f2v_idx_n[0] + 12;

      f2v_lst_n[0] = f2v_lst_c[0];
      f2v_lst_n[1] = f2v_lst_c[1];
      f2v_lst_n[2] = f_v_idx[f_id];
      f2v_lst_n[3] = f2v_lst_c[7];

      f2v_lst_n[4] = f2v_lst_c[2];
      f2v_lst_n[5] = f2v_lst_c[3];
      f2v_lst_n[6] = f_v_idx[f_id];
      f2v_lst_n[7] = f2v_lst_c[1];

      f2v_lst_n[8]  = f2v_lst_c[4];
      f2v_lst_n[9]  = f2v_lst_c[5];
      f2v_lst_n[10] = f_v_idx[f_id];
      f2v_lst_n[11] = f2v_lst_c[3];

      f2v_lst_n[12] = f2v_lst_c[6];
      f2v_lst_n[13] = f2v_lst_c[7];
      f2v_lst_n[14] = f_v_idx[f_id];
      f2v_lst_n[15] = f2v_lst_c[5];
    }
    break;

  default:
    {
      switch(f_r_flag) {
      case CS_REFINE_POLYGON:
        {
          int n_tria = fvm_triangulate_polygon(3, /* dim */
                                               0, /* base */
                                               n_fv_c,
                                               vtx_coords,
                                               nullptr,
                                               f2v_lst_c,
                                               FVM_TRIANGULATE_MESH_DEF,
                                               f2v_lst_n,
                                               t_state);

          if (n_tria != n_fv_c-2)
            bft_error(__FILE__, __LINE__, 0,
                      _("Error triangulating polygonal face.\n"
                        "This may be due to excessive warping."));

          for (cs_lnum_t i = 1; i < n_tria; i++)
            f2v_idx_n[i] = f2v_idx_n[0] + 3*i;
        }
        break;

      case CS_REFINE_POLYGON_T:
        {
          cs_lnum_t vc = f_v_idx[f_id];

          for (cs_lnum_t i = 0; i < n_fv_c; i++) {
            f2v_idx_n[i+1]     = f2v_idx_n[0] + i*3;

            f2v_lst_n[i*3]     = f2v_lst_c[i];
            f2v_lst_n[i*3 + 1] = f2v_lst_c[(i+1)%n_fv_c];
            f2v_lst_n[i*3 + 2] = vc;
          }
        }
        break;

      case CS_REFINE_POLYGON_Q:
        {
          /* Corners should have a refinement generation lower or equal
             to the initial face refinement level, so can be identified
             using this property. In case edges are subdivided multiple
             times, a mid-point should be the vertex with lowest
             refinement generation between 2 corners.

             Note that no face vertex should have a lower refinement
             generation than the first one, as we have rotated its vertices
             to ensure this in a prior step (in _flag_faces_and_edges).
          */

          cs_lnum_t mid_idx_first = n_fv_c-1;
          char r_lv_min = 127;  // Maximum allowable value (for unlocated).

          /* First corner: find mid-point of previous segment */

          {
            for (cs_lnum_t i = n_fv_c-1; i > 0; i--) {
              cs_lnum_t v_id = f2v_lst_c[i];
              char r_lv = vtx_r_gen[v_id];
              if (r_lv <= f_r_level_o)
                break;
              else if (r_lv < r_lv_min) {
                r_lv_min = r_lv;
                mid_idx_first = i;
              }
            }
          }

          cs_lnum_t idx_c0 = 0;
          cs_lnum_t idx_mid0 = mid_idx_first;
          cs_lnum_t idx_mid1 = 1;

          cs_lnum_t n_corner = 0;
          cs_lnum_t n = 0;
          r_lv_min = 127; // Reset for next search.

          /* Now loop on face vertices */

          for (cs_lnum_t i = 1; i < mid_idx_first; i++) {

            cs_lnum_t v_id = f2v_lst_c[i];
            char r_lv = vtx_r_gen[v_id];

            // Reached corner: mid-point between this corner and
            // The previous one should be known
            if (r_lv <= f_r_level_o) {
              cs_lnum_t idx_c1 = i;
              cs_lnum_t idx_c2 = (idx_c0 > 0) ? idx_c0 : n_fv_c;

              for (cs_lnum_t j = idx_c0; j <= idx_mid1; j++) {
                f2v_lst_n[n++] = f2v_lst_c[j];
              }
              f2v_lst_n[n++] = f_v_idx[f_id];
              for (cs_lnum_t j = idx_mid0; j < idx_c2; j++) {
                f2v_lst_n[n++] = f2v_lst_c[j];
              }
              f2v_idx_n[n_corner + 1] = f2v_idx_n[0] + n;
              n_corner += 1;

              // Prepare for next corner
              r_lv_min = 127;
              idx_mid0 = idx_mid1;
              idx_c0 = idx_c1;
            }

            // Non-corner point; might be a mid-point.
            else if (r_lv < r_lv_min) {
              r_lv_min = r_lv;
              idx_mid1 = i;
            }

          }

          // Last corner
          {
            idx_mid1  = mid_idx_first;

            for (cs_lnum_t j = idx_c0; j <= idx_mid1; j++) {
              f2v_lst_n[n++] = f2v_lst_c[j];
            }
            f2v_lst_n[n++] = f_v_idx[f_id];
            for (cs_lnum_t j = idx_mid0; j < idx_c0; j++) {
              f2v_lst_n[n++] = f2v_lst_c[j];
            }
            f2v_idx_n[n_corner + 1] = f2v_idx_n[0] + n;
            n_corner += 1;
          }

        }
        break;

      default:
        break;

      }
    }
  }

  if (f2v_lst_c != _f2v_lst_c)
    CS_FREE(f2v_lst_c);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute connectivity sizes for a given faces set.
 *
 * The face refinement flag may also be adjusted for faces with simple
 * refinement templates containing some edges which have further been marked
 * for subdivision due to a neighboring faces subdivision.
 *
 * For example, for a simple quad subdivision with 8 vertices on the contour:
 *   3---6---2
 *   |   |   |
 *   7---8---5
 *   |   |   |
 *   0---4---1
 *
 * if edge 0-4 was already subdivided into edges 0-4 and 4-1 du to a higher
 * refinement level in the adjecent face, and one of the adjacent faces is
 * further refined (so the edge is further subdivided), the face contour
 * will now contain at least one additional vertex.
 *
 * \param[in]       v2v          vertex->vertex adjacency
 * \param[in]       vtx_r_gen    vertices refinement generation
 * \param[in]       f_r_level_o  old faces refinement level
 * \param[in]       e_v_idx      for each edge, start index of added vertices
 * \param[in]       f_v_idx      for each face, start index of added vertices
 * \param[in]       n_faces      number of faces
 * \param[in, out]  f_r_flag     face refinement type flag
 * \param[in]       f2v_idx      face->vertices index
 * \param[in]       f2v_lst      face->vertices connectivity
 * \param[out]      f_o2n_idx    old to new faces index
 * \param[out]      f_o2n_connect_size   old to new faces connectivity size
 */
/*----------------------------------------------------------------------------*/

static void
_subdivided_faces_sizes(const cs_adjacency_t       *v2v,
                        const char                  vtx_r_gen[],
                        const char                  f_r_level_o[],
                        const cs_lnum_t             e_v_idx[],
                        cs_lnum_t                   n_faces,
                        cs_mesh_refine_type_t       f_r_flag[],
                        const cs_lnum_t             f2v_idx[],
                        const cs_lnum_t             f2v_lst[],
                        cs_lnum_t         *restrict f_o2n_idx,
                        cs_lnum_t         *restrict f_o2n_connect_size)
{
  f_o2n_idx[0] = 0;

# pragma omp parallel for  if(n_faces > CS_THR_MIN)
  for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++) {

    cs_lnum_t s_id = f2v_idx[f_id];
    cs_lnum_t n_fv = f2v_idx[f_id+1] - s_id;

    _subdivided_face_sizes(n_fv,
                           f_r_flag[f_id],
                           f_r_level_o[f_id],
                           f2v_lst + s_id,
                           v2v,
                           vtx_r_gen,
                           e_v_idx,
                           f_o2n_idx + f_id + 1,
                           f_o2n_connect_size + f_id +1);

  }

  _counts_to_index(n_faces, f_o2n_idx);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build subdivided triangle face lookup.
 *
 *       2
 *      / \
 *     5-- 4
 *    / \ / \
 *   0---3---1
 *
 * Numberings are based on the sub-face numbering relation based on
 * CS_REFINE_TRIA in _subdivide_face
 *
 * \param[in]   s_id          start id of first face in subset
 *                            (4 subfaces are adjacent)
 * \param[in]   start_vertex  vertex position (0-3) of first vertex
 *                            matching reference (for permutation)
 * \param[in]   sgn           parent face orientation (-1 if inverted)
 * \param[in]   f_vtx_idx     face->vertices connectivity index
 * \param[out]  tria_vtx      lookup face->vertices connectivity
 */
/*----------------------------------------------------------------------------*/

static void
_subdivided_tria(cs_lnum_t        s_id,
                 cs_lnum_t        start_vertex,
                 int              sgn,
                 const cs_lnum_t  f_vtx_idx[],
                 const cs_lnum_t  f_vtx[],
                 cs_lnum_t        tria_vtx[6])
{
  cs_lnum_t s_id_f;
  cs_lnum_t l = (3 - start_vertex)%3;

  if (sgn > 0) {

    s_id_f = f_vtx_idx[s_id];
    tria_vtx[l%3]     = f_vtx[s_id_f];
    tria_vtx[l%3 + 3] = f_vtx[s_id_f + 1];

    s_id_f = f_vtx_idx[s_id + 1];
    tria_vtx[(l+1)%3]     = f_vtx[s_id_f];
    tria_vtx[(l+1)%3 + 3] = f_vtx[s_id_f + 1];

    s_id_f = f_vtx_idx[s_id + 2];
    tria_vtx[(l+2)%3]     = f_vtx[s_id_f];
    tria_vtx[(l+2)%3 + 3] = f_vtx[s_id_f + 1];

  }
  else {

    s_id_f = f_vtx_idx[s_id + 2];
    tria_vtx[l%3]     = f_vtx[s_id_f];
    tria_vtx[l%3 + 3] = f_vtx[s_id_f + 2];

    s_id_f = f_vtx_idx[s_id + 1];
    tria_vtx[(l+1)%3]     = f_vtx[s_id_f];
    tria_vtx[(l+1)%3 + 3] = f_vtx[s_id_f + 2];

    s_id_f = f_vtx_idx[s_id];
    tria_vtx[(l+2)%3]     = f_vtx[s_id_f];
    tria_vtx[(l+2)%3 + 3] = f_vtx[s_id_f + 2];

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build subdivided triangle face lookup.
 *
 *       2
 *      / \
 *     5-6-4
 *    /  |  \
 *   0---3---1
 *
 * Numberings are based on the sub-face numbering relation based on
 * CS_REFINE_TRIA_H in _subdivide_face
 *
 * \param[in]   s_id          start id of first face in subset
 *                            (3 subfaces are adjacent)
 * \param[in]   start_vertex  vertex position (0-3) of first vertex
 *                            matching reference (for permutation)
 * \param[in]   sgn           parent face orientation (-1 if inverted)
 * \param[in]   f_vtx_idx     face->vertices connectivity index
 * \param[out]  quad_vtx      lookup face->vertices connectivity
 */
/*----------------------------------------------------------------------------*/

static void
_subdivided_tria_h(cs_lnum_t        s_id,
                   cs_lnum_t        start_vertex,
                   int              sgn,
                   const cs_lnum_t  f_vtx_idx[],
                   const cs_lnum_t  f_vtx[],
                   cs_lnum_t        quad_vtx[7])
{
  cs_lnum_t s_id_f;
  cs_lnum_t l = (3 - start_vertex)%3;

  if (sgn > 0) {

    s_id_f = f_vtx_idx[s_id];
    quad_vtx[6] = f_vtx[s_id_f + 2];
    quad_vtx[l%3]     = f_vtx[s_id_f];
    quad_vtx[l%3 + 3] = f_vtx[s_id_f + 1];

    s_id_f = f_vtx_idx[s_id + 1];
    quad_vtx[(l+1)%3]     = f_vtx[s_id_f];
    quad_vtx[(l+1)%3 + 3] = f_vtx[s_id_f + 1];

    s_id_f = f_vtx_idx[s_id + 2];
    quad_vtx[(l+2)%3]     = f_vtx[s_id_f];
    quad_vtx[(l+2)%3 + 3] = f_vtx[s_id_f + 1];

  }
  else {

    s_id_f = f_vtx_idx[s_id + 2];
    quad_vtx[6] = f_vtx[s_id_f + 2];
    quad_vtx[l%3]     = f_vtx[s_id_f];
    quad_vtx[l%3 + 3] = f_vtx[s_id_f + 3];

    s_id_f = f_vtx_idx[s_id + 1];
    quad_vtx[(l+1)%3]     = f_vtx[s_id_f];
    quad_vtx[(l+1)%3 + 3] = f_vtx[s_id_f + 3];

    s_id_f = f_vtx_idx[s_id];
    quad_vtx[(l+2)%3]     = f_vtx[s_id_f];
    quad_vtx[(l+2)%3 + 3] = f_vtx[s_id_f + 3];

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build subdivided quad face lookup.
 *
 *   3---6---2
 *   |   |   |
 *   7---8---5
 *   |   |   |
 *   0---4---1
 *
 * Numberings are based on the sub-face numbering relation based on
 * CS_REFINE_QUAD in _subdivide_face
 *
 * \param[in]   s_id          start id of first face in subset
 *                            (4 subfaces are adjacent)
 * \param[in]   start_vertex  vertex position (0-3) of first vertex
 *                            matching reference (for permutation)
 * \param[in]   sgn           parent face orientation (-1 if inverted)
 * \param[in]   f_vtx_idx     face->vertices connectivity index
 * \param[out]  quad_vtx      lookup face->vertices connectivity
 */
/*----------------------------------------------------------------------------*/

static void
_subdivided_quad(cs_lnum_t        s_id,
                 cs_lnum_t        start_vertex,
                 int              sgn,
                 const cs_lnum_t  f_vtx_idx[],
                 const cs_lnum_t  f_vtx[],
                 cs_lnum_t        quad_vtx[9])
{
  cs_lnum_t s_id_f;
  cs_lnum_t l = (4 - start_vertex)%4;

  if (sgn > 0) {

    s_id_f = f_vtx_idx[s_id];
    quad_vtx[8] = f_vtx[s_id_f + 2];
    quad_vtx[l%4] = f_vtx[s_id_f];
    quad_vtx[l%4 + 4] = f_vtx[s_id_f + 1];

    s_id_f = f_vtx_idx[s_id + 1];
    quad_vtx[(l+1)%4]     = f_vtx[s_id_f];
    quad_vtx[(l+1)%4 + 4] = f_vtx[s_id_f + 1];

    s_id_f = f_vtx_idx[s_id + 2];
    quad_vtx[(l+2)%4]     = f_vtx[s_id_f];
    quad_vtx[(l+2)%4 + 4] = f_vtx[s_id_f + 1];

    s_id_f = f_vtx_idx[s_id + 3];
    quad_vtx[(l+3)%4]     = f_vtx[s_id_f];
    quad_vtx[(l+3)%4 + 4] = f_vtx[s_id_f + 1];

  }
  else {

    s_id_f = f_vtx_idx[s_id + 3];
    quad_vtx[8] = f_vtx[s_id_f + 2];
    quad_vtx[l%4] = f_vtx[s_id_f];
    quad_vtx[l%4 + 4] = f_vtx[s_id_f + 3];

    s_id_f = f_vtx_idx[s_id + 2];
    quad_vtx[(l+1)%4]     = f_vtx[s_id_f];
    quad_vtx[(l+1)%4 + 4] = f_vtx[s_id_f + 3];

    s_id_f = f_vtx_idx[s_id + 1];
    quad_vtx[(l+2)%4]     = f_vtx[s_id_f];
    quad_vtx[(l+2)%4 + 4] = f_vtx[s_id_f + 3];

    s_id_f = f_vtx_idx[s_id];
    quad_vtx[(l+3)%4]     = f_vtx[s_id_f];
    quad_vtx[(l+3)%4 + 4] = f_vtx[s_id_f + 3];

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define interior triangle face and associate with cells.
 *
 * \param[in, out]  m           pointer to mesh structure
 * \param[in]       c_id_0      cell id with face normal out
 * \param[in]       c_id_1      cell id with face normal in
 * \param[in]       f_id        assigned face id
 * \param[in, out]  vertex_ids  associated vertex ids (may be permuted)
 * \param[in]       c_f_ranges  cell added faces range
 */
/*----------------------------------------------------------------------------*/

inline static void
_add_interior_face_tria(const cs_mesh_t   *m,
                        cs_lnum_t          c_id_0,
                        cs_lnum_t          c_id_1,
                        cs_lnum_t          f_id,
                        cs_lnum_t          vertex_ids[],
                        const cs_lnum_t    c_f_range[2])
{
  cs_lnum_t *_vtx_lst = m->i_face_vtx_lst + m->i_face_vtx_idx[f_id];

  for (cs_lnum_t i = 0; i < 3; i++)
    _vtx_lst[i] = vertex_ids[i];

  if (f_id < c_f_range[1]-1)
    m->i_face_vtx_idx[f_id+1] = m->i_face_vtx_idx[f_id] + 3;
  else
    assert(m->i_face_vtx_idx[f_id+1] == m->i_face_vtx_idx[f_id] + 3);

  m->i_face_cells[f_id][0] = c_id_0;
  m->i_face_cells[f_id][1] = c_id_1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define interior quadrangle face and associate with cells.
 *
 * \param[in, out]  m           pointer to mesh structure
 * \param[in]       c_id_0      cell id with face normal out
 * \param[in]       c_id_1      cell id with face normal in
 * \param[in]       f_id        assigned face id
 * \param[in, out]  vertex_ids  associated vertex ids (may be permuted)
 * \param[in]       c_f_ranges  cell added faces range
 */
/*----------------------------------------------------------------------------*/

inline static void
_add_interior_face_quad(const cs_mesh_t   *m,
                        cs_lnum_t          c_id_0,
                        cs_lnum_t          c_id_1,
                        cs_lnum_t          f_id,
                        cs_lnum_t          vertex_ids[],
                        const cs_lnum_t    c_f_range[2])
{
  cs_lnum_t *_vtx_lst = m->i_face_vtx_lst + m->i_face_vtx_idx[f_id];

  for (cs_lnum_t i = 0; i < 4; i++)
    _vtx_lst[i] = vertex_ids[i];

  if (f_id < c_f_range[1]-1)
    m->i_face_vtx_idx[f_id+1] = m->i_face_vtx_idx[f_id] + 4;
  else
    assert(m->i_face_vtx_idx[f_id+1] == m->i_face_vtx_idx[f_id] + 4);

  m->i_face_cells[f_id][0] = c_id_0;
  m->i_face_cells[f_id][1] = c_id_1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Match a cell's triangle faces subdivision and partially update
 *        face->cell connectivity.
 *
 * The connectivity of exisiting (subdivided) boundary of interior
 * faces (on the exterior of the cell, not added interior faces)
 * is updated, assuming a mapping defined by the c_id_shift argument.
 *
 * The cells->faces adjacency is updated as faces are locally renumbered
 * in canonical order for known cell types.
 *
 * \param[in, out]  m             pointer to mesh structure
 * \param[in]       c_f_id_s      start id of cell's local faces to handle
 * \param[in]       c_f_id_e      past-end id of cell's local faces to handle
 * \param[in]       n_b_f_ini     old number of boundary faces
 * \param[in]       c_o2n_idx     old to new cells index
 * \param[in]       i_face_o2n_idx  old to new interior faces index
 * \param[in]       b_face_o2n_idx  old to new boundary faces index
 * \param[in]       c2f           cells->faces adjacency (boundary faces first)
 * \param[in]       c2f2v_start   start id for adjacent face vertices
 *                                definitions (initialized to 0)
 * \param[in]       c_v_idx       for each cell, start index of added vertices
 * \param[in]       c_f_n_idx     cells to new faces index
 * \param[in]       c_id_shift    subface to subcell mapping
 * \param[out]      tria_vtx      subdivided triangle connectivity
 *                                (see \ref _subdivided_tria)
 */
/*----------------------------------------------------------------------------*/

static void
_subdivide_cell_tria_faces(const cs_mesh_t    *m,
                           cs_lnum_t           c_f_id_s,
                           cs_lnum_t           c_f_id_e,
                           cs_lnum_t           cell_id,
                           cs_lnum_t           n_b_f_ini,
                           const cs_lnum_t     c_o2n_idx[],
                           const cs_lnum_t     i_face_o2n_idx[],
                           const cs_lnum_t     b_face_o2n_idx[],
                           cs_adjacency_t     *c2f,
                           const cs_lnum_t     c2f2v_start[],
                           const cs_lnum_t     c_id_shift[][4],
                           cs_lnum_t           tria_vtx[][6])
{
  const cs_lnum_t s_id_c = c2f->idx[cell_id];

  /* Loop on faces; face->sub-face numbering relations based
     on CS_REFINE_TRIA in _subdivide_face */

  for (cs_lnum_t i = c_f_id_s; i < c_f_id_e; i++) {

    cs_lnum_t s_id, e_id;
    cs_lnum_t f_id_o = c2f->ids[s_id_c + i];
    cs_lnum_t sgn = c2f->sgn[s_id_c + i];

    cs_lnum_t *_tria_vtx = tria_vtx[i-c_f_id_s];

    const cs_lnum_t *f_vtx_idx, *f_vtx;

    cs_lnum_t s_f_id_shift = (3 - c2f2v_start[s_id_c + i])%3;

    if (f_id_o < n_b_f_ini) {
      s_id = b_face_o2n_idx[f_id_o];
      e_id = b_face_o2n_idx[f_id_o + 1];
      assert(e_id - s_id == 4);
      f_vtx_idx = m->b_face_vtx_idx;
      f_vtx = m->b_face_vtx_lst;
      assert(sgn > 0);
      for (cs_lnum_t k = 0; k < 3; k++) {
        cs_lnum_t f_id = s_id + k;
        assert(m->b_face_cells[f_id] == c_o2n_idx[cell_id]);
        m->b_face_cells[f_id] += c_id_shift[i][(s_f_id_shift+k)%3];
      }
      m->b_face_cells[s_id+3] += c_id_shift[i][3];
    }
    else {
      s_id = i_face_o2n_idx[f_id_o - n_b_f_ini];
      e_id = i_face_o2n_idx[f_id_o - n_b_f_ini + 1];
      assert(e_id - s_id == 4);
      f_vtx_idx = m->i_face_vtx_idx;
      f_vtx = m->i_face_vtx_lst;
      if (sgn > 0) {
        for (cs_lnum_t k = 0; k < 3; k++) {
          cs_lnum_t f_id = s_id + k;
          assert(m->i_face_cells[f_id][0] == c_o2n_idx[cell_id]);
          m->i_face_cells[f_id][0] += c_id_shift[i][(s_f_id_shift+k)%3];
        }
        m->i_face_cells[s_id+3][0] += c_id_shift[i][3];
      }
      else {
        for (cs_lnum_t k = 0; k < 3; k++) {
          cs_lnum_t f_id = s_id + 2 - k;
          assert(m->i_face_cells[f_id][1] == c_o2n_idx[cell_id]);
          m->i_face_cells[f_id][1] += c_id_shift[i][(s_f_id_shift+k)%3];
        }
        m->i_face_cells[s_id+3][1] += c_id_shift[i][3];
      }
    }

    _subdivided_tria(s_id,
                     c2f2v_start[s_id_c + i],
                     sgn,
                     f_vtx_idx,
                     f_vtx,
                     _tria_vtx);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Match a cell's triangle faces subdivision and partially update
 *        face->cell connectivity.
 *
 * The connectivity of exisiting (subdivided) boundary of interior
 * faces (on the exterior of the cell, not added interior faces)
 * is updated, assuming a mapping defined by the c_id_shift argument.
 *
 * The cells->faces adjacency is updated as faces are locally renumbered
 * in canonical order for known cell types.
 *
 * \param[in, out]  m             pointer to mesh structure
 * \param[in]       c_f_id_s      start id of cell's local faces to handle
 * \param[in]       c_f_id_e      past-end id of cell's local faces to handle
 * \param[in]       n_b_f_ini     old number of boundary faces
 * \param[in]       c_o2n_idx     old to new cells index
 * \param[in]       i_face_o2n_idx  old to new interior faces index
 * \param[in]       b_face_o2n_idx  old to new boundary faces index
 * \param[in]       c2f           cells->faces adjacency (boundary faces first)
 * \param[in]       c2f2v_start   start id for adjacent face vertices
 *                                definitions (initialized to 0)
 * \param[in]       c_v_idx       for each cell, start index of added vertices
 * \param[in]       c_f_n_idx     cells to new faces index
 * \param[in]       c_id_shift    subface to subcell mapping
 * \param[out]      tria_vtx      subdivided triangle connectivity
 *                                (see \ref _subdivided_tria)
 */
/*----------------------------------------------------------------------------*/

static void
_subdivide_cell_tria_q_faces(const cs_mesh_t  *m,
                             cs_lnum_t           c_f_id_s,
                             cs_lnum_t           c_f_id_e,
                             cs_lnum_t           cell_id,
                             cs_lnum_t           n_b_f_ini,
                             const cs_lnum_t     c_o2n_idx[],
                             const cs_lnum_t     i_face_o2n_idx[],
                             const cs_lnum_t     b_face_o2n_idx[],
                             cs_adjacency_t     *c2f,
                             const cs_lnum_t     c2f2v_start[],
                             const cs_lnum_t     c_id_shift[][3],
                             cs_lnum_t           tria_vtx[][7])
{
  const cs_lnum_t s_id_c = c2f->idx[cell_id];

  /* Loop on faces; face->sub-face numbering relations based
     on CS_REFINE_TRIA in _subdivide_face */

  for (cs_lnum_t i = c_f_id_s; i < c_f_id_e; i++) {

    cs_lnum_t s_id, e_id;
    cs_lnum_t f_id_o = c2f->ids[s_id_c + i];
    cs_lnum_t sgn = c2f->sgn[s_id_c + i];

    cs_lnum_t *_tria_vtx = tria_vtx[i-c_f_id_s];

    const cs_lnum_t *f_vtx_idx, *f_vtx;

    cs_lnum_t s_f_id_shift = (3 - c2f2v_start[s_id_c + i])%3;

    if (f_id_o < n_b_f_ini) {
      s_id = b_face_o2n_idx[f_id_o];
      e_id = b_face_o2n_idx[f_id_o + 1];
      assert(e_id - s_id == 3);
      f_vtx_idx = m->b_face_vtx_idx;
      f_vtx = m->b_face_vtx_lst;
      assert(sgn > 0);
      for (cs_lnum_t k = 0; k < 3; k++) {
        cs_lnum_t f_id = s_id + k;
        assert(m->b_face_cells[f_id] == c_o2n_idx[cell_id]);
        m->b_face_cells[f_id] += c_id_shift[i][(s_f_id_shift+k)%3];
      }

    }
    else {
      s_id = i_face_o2n_idx[f_id_o - n_b_f_ini];
      e_id = i_face_o2n_idx[f_id_o - n_b_f_ini + 1];
      assert(e_id - s_id == 3);
      f_vtx_idx = m->i_face_vtx_idx;
      f_vtx = m->i_face_vtx_lst;
      if (sgn > 0) {
        for (cs_lnum_t k = 0; k < 3; k++) {
          cs_lnum_t f_id = s_id + k;
          assert(m->i_face_cells[f_id][0] == c_o2n_idx[cell_id]);
          m->i_face_cells[f_id][0] += c_id_shift[i][(s_f_id_shift+k)%3];
        }

      }
      else {
        for (cs_lnum_t k = 0; k < 3; k++) {
          cs_lnum_t f_id = s_id + 2 - k;
          assert(m->i_face_cells[f_id][1] == c_o2n_idx[cell_id]);
          m->i_face_cells[f_id][1] += c_id_shift[i][(s_f_id_shift+k)%3];
        }

      }
    }

    _subdivided_tria_h(s_id,
                     c2f2v_start[s_id_c + i],
                     sgn,
                     f_vtx_idx,
                     f_vtx,
                     _tria_vtx);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Match a cell's quadrangle faces subdivision and partially update
 *        face->cell connectivity.
 *
 * The connectivity of exisiting (subdivided) boundary of interior
 * faces (on the exterior of the cell, not added interior faces)
 * is updated, assuming a mapping defined by the c_id_shift argument.
 *
 * The cells->faces adjacency is updated as faces are locally renumbered
 * in canonical order for known cell types.
 *
 * \param[in, out]  m             pointer to mesh structure
 * \param[in]       c_f_id_s      start id of cell's local faces to handle
 * \param[in]       c_f_id_e      past-end id of cell's local faces to handle
 * \param[in]       n_b_f_ini     old number of boundary faces
 * \param[in]       c_o2n_idx     old to new cells index
 * \param[in]       i_face_o2n_idx  old to new interior faces index
 * \param[in]       b_face_o2n_idx  old to new boundary faces index
 * \param[in]       c2f           cells->faces adjacency (boundary faces first)
 * \param[in]       c2f2v_start   start id for adjacent face vertices
 *                                definitions (initialized to 0)
 * \param[in]       c_v_idx       for each cell, start index of added vertices
 * \param[in]       c_f_n_idx     cells to new faces index
 * \param[in]       c_id_shift    subface to subcell mapping
 * \param[out]      quad_vtx      subdivided quadrangle connectivity
 *                                (see \ref _subdivided_quad)
 */
/*----------------------------------------------------------------------------*/

static void
_subdivide_cell_quad_faces(const cs_mesh_t    *m,
                           cs_lnum_t           c_f_id_s,
                           cs_lnum_t           c_f_id_e,
                           cs_lnum_t           cell_id,
                           cs_lnum_t           n_b_f_ini,
                           const cs_lnum_t     c_o2n_idx[],
                           const cs_lnum_t     i_face_o2n_idx[],
                           const cs_lnum_t     b_face_o2n_idx[],
                           cs_adjacency_t     *c2f,
                           const cs_lnum_t     c2f2v_start[],
                           const cs_lnum_t     c_id_shift[][4],
                           cs_lnum_t           quad_vtx[][9])
{
  const cs_lnum_t s_id_c = c2f->idx[cell_id];

  /* Loop on faces; face->sub-face numbering relations based
     on CS_REFINE_QUAD in _subdivide_face */

  for (cs_lnum_t i = c_f_id_s; i < c_f_id_e; i++) {

    cs_lnum_t s_id, e_id;
    cs_lnum_t f_id_o = c2f->ids[s_id_c + i];
    cs_lnum_t sgn = c2f->sgn[s_id_c + i];

    cs_lnum_t *_quad_vtx = quad_vtx[i-c_f_id_s];

    const cs_lnum_t *f_vtx_idx, *f_vtx;

    cs_lnum_t s_f_id_shift = (4 - c2f2v_start[s_id_c + i])%4;

    if (f_id_o < n_b_f_ini) {
      s_id = b_face_o2n_idx[f_id_o];
      e_id = b_face_o2n_idx[f_id_o + 1];
      assert(e_id - s_id == 4);
      f_vtx_idx = m->b_face_vtx_idx;
      f_vtx = m->b_face_vtx_lst;
      assert(sgn > 0);
      for (cs_lnum_t k = 0; k < 4; k++) {
        cs_lnum_t f_id = s_id + k;
        assert(m->b_face_cells[f_id] == c_o2n_idx[cell_id]);
        m->b_face_cells[f_id] += c_id_shift[i][(s_f_id_shift+k)%4];
      }
    }
    else {
      s_id = i_face_o2n_idx[f_id_o - n_b_f_ini];
      e_id = i_face_o2n_idx[f_id_o - n_b_f_ini + 1];
      assert(e_id - s_id == 4);
      f_vtx_idx = m->i_face_vtx_idx;
      f_vtx = m->i_face_vtx_lst;
      if (sgn > 0) {
        for (cs_lnum_t k = 0; k < 4; k++) {
          cs_lnum_t f_id = s_id + k;
          assert(m->i_face_cells[f_id][0] == c_o2n_idx[cell_id]);
          m->i_face_cells[f_id][0] += c_id_shift[i][(s_f_id_shift+k)%4];
        }
      }
      else {
        for (cs_lnum_t k = 0; k < 4; k++) {
          cs_lnum_t f_id = s_id + 3 - k;
          assert(m->i_face_cells[f_id][1] == c_o2n_idx[cell_id]);
          m->i_face_cells[f_id][1] += c_id_shift[i][(s_f_id_shift+k)%4];
        }
      }
    }

    _subdivided_quad(s_id,
                     c2f2v_start[s_id_c + i],
                     sgn,
                     f_vtx_idx,
                     f_vtx,
                     _quad_vtx);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Subdivide a given tetrahedron, updating faces.
 *
 * The cells->faces adjacency is updated as faces are locally renumbered
 * in canonical order for known cell types.
 *
 * \param[in, out]  m             pointer to mesh structure
 * \param[in]       cell_id       (old) cell id (0 to n-1)
 * \param[in]       n_b_f_ini     old number of boundary faces
 * \param[in]       c_o2n_idx     old to new cells index
 * \param[in]       i_face_o2n_idx  old to new interior faces index
 * \param[in]       b_face_o2n_idx  old to new boundary faces index
 * \param[in]       c2f           cells->faces adjacency (boundary faces first)
 * \param[in]       c2f2v_start   start id for adjacent face vertices
 *                                definitions (initialized to 0)
 * \param[in]       c_v_idx       for each cell, start index of added vertices
 * \param[in]       c_f_n_idx     cells to new faces index
 */
/*----------------------------------------------------------------------------*/

static void
_subdivide_cell_tetra(const cs_mesh_t              *m,
                      cs_lnum_t                     cell_id,
                      cs_lnum_t                     n_b_f_ini,
                      const cs_lnum_t               c_o2n_idx[],
                      const cs_lnum_t               i_face_o2n_idx[],
                      const cs_lnum_t               b_face_o2n_idx[],
                      cs_adjacency_t               *c2f,
                      const cs_lnum_t               c2f2v_start[],
                      const cs_lnum_t               c_v_idx[],
                      const cs_lnum_t               c_f_n_idx[])
{
  CS_UNUSED(c_v_idx);

  assert(c2f->idx[cell_id+1] - c2f->idx[cell_id] == 4);

  static const cs_lnum_t c_id_shift[4][4] = {{0, 2, 1, 4},
                                             {0, 1, 3, 5},
                                             {1, 2, 3, 6},
                                             {2, 0, 3, 7}};

  /* Case of tetrahedron split into 8 smaller prisms */

  cs_lnum_t new_cell_id = c_o2n_idx[cell_id];

  cs_lnum_t tria_vtx[4][6];

  /* Loop on faces; face->sub-face numbering relations based
     on CS_REFINE_TRIA in _subdivide_face */

  _subdivide_cell_tria_faces(m,
                             0, 4, /* tria range */
                             cell_id,
                             n_b_f_ini,
                             c_o2n_idx,
                             i_face_o2n_idx,
                             b_face_o2n_idx,
                             c2f,
                             c2f2v_start,
                             c_id_shift,
                             tria_vtx);

  cs_lnum_t v_ids[3];

  v_ids[0] = tria_vtx[0][5];
  v_ids[1] = tria_vtx[0][3];
  v_ids[2] = tria_vtx[1][5];
  _add_interior_face_tria(m,
                          new_cell_id,
                          new_cell_id+5,
                          c_f_n_idx[cell_id],
                          v_ids,
                          c_f_n_idx + cell_id);

  v_ids[0] = tria_vtx[1][4];
  v_ids[1] = tria_vtx[2][3];
  v_ids[2] = tria_vtx[1][3];
  _add_interior_face_tria(m,
                          new_cell_id+1,
                          new_cell_id+4,
                          c_f_n_idx[cell_id] + 1,
                          v_ids,
                          c_f_n_idx + cell_id);

  v_ids[0] = tria_vtx[2][4];
  v_ids[1] = tria_vtx[3][3];
  v_ids[2] = tria_vtx[2][3];
  _add_interior_face_tria(m,
                          new_cell_id+2,
                          new_cell_id+6,
                          c_f_n_idx[cell_id] + 2,
                          v_ids,
                          c_f_n_idx + cell_id);

  v_ids[0] = tria_vtx[3][5];
  v_ids[1] = tria_vtx[1][4];
  v_ids[2] = tria_vtx[3][4];
  _add_interior_face_tria(m,
                          new_cell_id+3,
                          new_cell_id+7,
                          c_f_n_idx[cell_id] + 3,
                          v_ids,
                          c_f_n_idx + cell_id);

  v_ids[0] = tria_vtx[0][3];
  v_ids[1] = tria_vtx[0][5];
  v_ids[2] = tria_vtx[1][4];
  _add_interior_face_tria(m,
                          new_cell_id+4,
                          new_cell_id+5,
                          c_f_n_idx[cell_id] + 4,
                          v_ids,
                          c_f_n_idx + cell_id);

  v_ids[0] = tria_vtx[0][4];
  v_ids[1] = tria_vtx[0][3];
  v_ids[2] = tria_vtx[2][5];
  _add_interior_face_tria(m,
                          new_cell_id+4,
                          new_cell_id+6,
                          c_f_n_idx[cell_id] + 5,
                          v_ids,
                          c_f_n_idx + cell_id);

  v_ids[0] = tria_vtx[1][5];
  v_ids[1] = tria_vtx[1][4];
  v_ids[2] = tria_vtx[3][3];
  _add_interior_face_tria(m,
                          new_cell_id+5,
                          new_cell_id+7,
                          c_f_n_idx[cell_id] + 6,
                          v_ids,
                          c_f_n_idx + cell_id);

  v_ids[0] = tria_vtx[2][5];
  v_ids[1] = tria_vtx[2][4];
  v_ids[2] = tria_vtx[3][3];
  _add_interior_face_tria(m,
                          new_cell_id+6,
                          new_cell_id+7,
                          c_f_n_idx[cell_id] + 7,
                          v_ids,
                          c_f_n_idx + cell_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Subdivide a given tetrahedron to hexahedra, updating faces.
 *
 * The cells->faces adjacency is updated as faces are locally renumbered
 * in canonical order for known cell types.
 *
 * \param[in, out]  m             pointer to mesh structure
 * \param[in]       cell_id       (old) cell id (0 to n-1)
 * \param[in]       n_b_f_ini     old number of boundary faces
 * \param[in]       c_o2n_idx     old to new cells index
 * \param[in]       i_face_o2n_idx  old to new interior faces index
 * \param[in]       b_face_o2n_idx  old to new boundary faces index
 * \param[in]       c2f           cells->faces adjacency (boundary faces first)
 * \param[in]       c2f2v_start   start id for adjacent face vertices
 *                                definitions (initialized to 0)
 * \param[in]       c_v_idx       for each cell, start index of added vertices
 * \param[in]       c_f_n_idx     cells to new faces index
 */
/*----------------------------------------------------------------------------*/

static void
_subdivide_cell_tetra_h(const cs_mesh_t              *m,
                        cs_lnum_t                     cell_id,
                        cs_lnum_t                     n_b_f_ini,
                        const cs_lnum_t               c_o2n_idx[],
                        const cs_lnum_t               i_face_o2n_idx[],
                        const cs_lnum_t               b_face_o2n_idx[],
                        cs_adjacency_t               *c2f,
                        const cs_lnum_t               c2f2v_start[],
                        const cs_lnum_t               c_v_idx[],
                        const cs_lnum_t               c_f_n_idx[])
{
  assert(c2f->idx[cell_id+1] - c2f->idx[cell_id] == 4);

  static const cs_lnum_t c_id_shift[4][3] = {{0, 2, 1},
                                             {0, 1, 3},
                                             {1, 2, 3},
                                             {2, 0, 3}};

  /* Case of tetrahedron split into 4 hexahedra */

  const cs_lnum_t c_vtx_id = c_v_idx[cell_id];

  cs_lnum_t new_cell_id = c_o2n_idx[cell_id];

  cs_lnum_t tria_vtx[4][7];

  /* Loop on faces; face->sub-face numbering relations based
     on CS_REFINE_TRIA in _subdivide_face */

  _subdivide_cell_tria_q_faces(m,
                               0, 4, /* tria range */
                               cell_id,
                               n_b_f_ini,
                               c_o2n_idx,
                               i_face_o2n_idx,
                               b_face_o2n_idx,
                               c2f,
                               c2f2v_start,
                               c_id_shift,
                               tria_vtx);

  cs_lnum_t v_ids[4];

  v_ids[0] = tria_vtx[0][5];
  v_ids[1] = tria_vtx[0][6];
  v_ids[2] = c_vtx_id;
  v_ids[3] = tria_vtx[1][6];
  _add_interior_face_quad(m,
                          new_cell_id,
                          new_cell_id+1,
                          c_f_n_idx[cell_id],
                          v_ids,
                          c_f_n_idx + cell_id);

  v_ids[0] = tria_vtx[2][3];
  v_ids[1] = tria_vtx[0][6];
  v_ids[2] = c_vtx_id;
  v_ids[3] = tria_vtx[2][6];
  _add_interior_face_quad(m,
                          new_cell_id+1,
                          new_cell_id+2,
                          c_f_n_idx[cell_id] + 1,
                          v_ids,
                          c_f_n_idx + cell_id);

  v_ids[0] = tria_vtx[3][3];
  v_ids[1] = tria_vtx[3][6];
  v_ids[2] = c_vtx_id;
  v_ids[3] = tria_vtx[0][6];
  _add_interior_face_quad(m,
                          new_cell_id,
                          new_cell_id+2,
                          c_f_n_idx[cell_id] + 2,
                          v_ids,
                          c_f_n_idx + cell_id);

  v_ids[0] = tria_vtx[1][5];
  v_ids[1] = tria_vtx[1][6];
  v_ids[2] = c_vtx_id;
  v_ids[3] = tria_vtx[3][6];
  _add_interior_face_quad(m,
                          new_cell_id,
                          new_cell_id+3,
                          c_f_n_idx[cell_id] + 3,
                          v_ids,
                          c_f_n_idx + cell_id);

  v_ids[0] = tria_vtx[2][5];
  v_ids[1] = tria_vtx[2][6];
  v_ids[2] = c_vtx_id;
  v_ids[3] = tria_vtx[1][6];
  _add_interior_face_quad(m,
                          new_cell_id+1,
                          new_cell_id+3,
                          c_f_n_idx[cell_id] + 4,
                          v_ids,
                          c_f_n_idx + cell_id);

  v_ids[0] = tria_vtx[3][5];
  v_ids[1] = tria_vtx[3][6];
  v_ids[2] = c_vtx_id;
  v_ids[3] = tria_vtx[2][6];
  _add_interior_face_quad(m,
                          new_cell_id+2,
                          new_cell_id+3,
                          c_f_n_idx[cell_id] + 5,
                          v_ids,
                          c_f_n_idx + cell_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Subdivide a given pyramid, updating faces.
 *
 * The cells->faces adjacency is updated as faces are locally renumbered
 * in canonical order for known cell types.
 *
 * \param[in, out]  m             pointer to mesh structure
 * \param[in]       cell_id       (old) cell id (0 to n-1)
 * \param[in]       n_b_f_ini     old number of boundary faces
 * \param[in]       c_o2n_idx     old to new cells index
 * \param[in]       i_face_o2n_idx  old to new interior faces index
 * \param[in]       b_face_o2n_idx  old to new boundary faces index
 * \param[in]       c2f           cells->faces adjacency (boundary faces first)
 * \param[in]       c2f2v_start   start id for adjacent face vertices
 *                                definitions (initialized to 0)
 * \param[in]       c_v_idx       for each cell, start index of added vertices
 * \param[in]       c_f_n_idx     cells to new faces index
 */
/*----------------------------------------------------------------------------*/

static void
_subdivide_cell_pyram(const cs_mesh_t              *m,
                      cs_lnum_t                     cell_id,
                      cs_lnum_t                     n_b_f_ini,
                      const cs_lnum_t               c_o2n_idx[],
                      const cs_lnum_t               i_face_o2n_idx[],
                      const cs_lnum_t               b_face_o2n_idx[],
                      cs_adjacency_t               *c2f,
                      const cs_lnum_t               c2f2v_start[],
                      const cs_lnum_t               c_v_idx[],
                      const cs_lnum_t               c_f_n_idx[])
{
  CS_UNUSED(c_v_idx);

  assert(c2f->idx[cell_id+1] - c2f->idx[cell_id] == 5);

  static const cs_lnum_t c_id_shift[5][4] = {{0, 3, 2, 1},
                                             {0, 1, 4, 6},
                                             {1, 2, 4, 7},
                                             {2, 3, 4, 8},
                                             {3, 0, 4, 9}};

  /* Case of pyramid split into 6 pyramids and 4 tetrahedra */

  cs_lnum_t new_cell_id = c_o2n_idx[cell_id];

  cs_lnum_t quad_vtx[1][9];
  cs_lnum_t tria_vtx[4][6];

  /* Loop on faces; face->sub-face numbering relations based
     on CS_REFINE_TRIA and CS_REFINE_QUAD in _subdivide_face */

  _subdivide_cell_quad_faces(m,
                             0, 1, /* quad range */
                             cell_id,
                             n_b_f_ini,
                             c_o2n_idx,
                             i_face_o2n_idx,
                             b_face_o2n_idx,
                             c2f,
                             c2f2v_start,
                             c_id_shift,
                             quad_vtx);

  _subdivide_cell_tria_faces(m,
                             1, 5, /* tria range */
                             cell_id,
                             n_b_f_ini,
                             c_o2n_idx,
                             i_face_o2n_idx,
                             b_face_o2n_idx,
                             c2f,
                             c2f2v_start,
                             c_id_shift,
                             tria_vtx);

  cs_lnum_t v_ids[4];

  v_ids[3] = tria_vtx[0][4];
  v_ids[2] = tria_vtx[1][4];
  v_ids[1] = tria_vtx[2][4];
  v_ids[0] = tria_vtx[3][4];
  _add_interior_face_quad(m,
                          new_cell_id+4,
                          new_cell_id+5,
                          c_f_n_idx[cell_id],
                          v_ids,
                          c_f_n_idx + cell_id);

  for (cs_lnum_t i = 0; i < 4; i++) {

    v_ids[0] = tria_vtx[i][3];
    v_ids[1] = quad_vtx[0][8];
    v_ids[2] = tria_vtx[i][5];
    _add_interior_face_tria(m,
                            new_cell_id+i,
                            new_cell_id+6+i,
                            c_f_n_idx[cell_id] + 1 + 3*i,
                            v_ids,
                            c_f_n_idx + cell_id);

    v_ids[0] = tria_vtx[i][3];
    v_ids[1] = tria_vtx[i][4];
    v_ids[2] = quad_vtx[0][8];
    _add_interior_face_tria(m,
                            new_cell_id+(1+i)%4,
                            new_cell_id+6+i,
                            c_f_n_idx[cell_id] + 1 + 3*i + 1,
                            v_ids,
                            c_f_n_idx + cell_id);

    v_ids[0] = tria_vtx[i][5];
    v_ids[1] = quad_vtx[0][8];
    v_ids[2] = tria_vtx[i][4];
    _add_interior_face_tria(m,
                            new_cell_id+5,
                            new_cell_id+6+i,
                            c_f_n_idx[cell_id] + 1 + 3*i + 2,
                            v_ids,
                            c_f_n_idx + cell_id);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Subdivide a given prism, updating faces.
 *
 * The cells->faces adjacency is updated as faces are locally renumbered
 * in canonical order for known cell types.
 *
 * \param[in, out]  m             pointer to mesh structure
 * \param[in]       cell_id       (old) cell id (0 to n-1)
 * \param[in]       n_b_f_ini     old number of boundary faces
 * \param[in]       c_o2n_idx     old to new cells index
 * \param[in]       i_face_o2n_idx  old to new interior faces index
 * \param[in]       b_face_o2n_idx  old to new boundary faces index
 * \param[in]       c2f           cells->faces adjacency (boundary faces first)
 * \param[in]       c2f2v_start   start id for adjacent face vertices
 *                                definitions (initialized to 0)
 * \param[in]       c_v_idx       for each cell, start index of added vertices
 * \param[in]       c_f_n_idx     cells to new faces index
 */
/*----------------------------------------------------------------------------*/

static void
_subdivide_cell_prism(const cs_mesh_t              *m,
                      cs_lnum_t                     cell_id,
                      cs_lnum_t                     n_b_f_ini,
                      const cs_lnum_t               c_o2n_idx[],
                      const cs_lnum_t               i_face_o2n_idx[],
                      const cs_lnum_t               b_face_o2n_idx[],
                      cs_adjacency_t               *c2f,
                      const cs_lnum_t               c2f2v_start[],
                      const cs_lnum_t               c_v_idx[],
                      const cs_lnum_t               c_f_n_idx[])
{
  CS_UNUSED(c_v_idx);

  assert(c2f->idx[cell_id+1] - c2f->idx[cell_id] == 5);

  static const cs_lnum_t c_id_shift[5][4] = {{0, 1, 5, 4},
                                             {1, 2, 6, 5},
                                             {2, 0, 4, 6},
                                             {0, 2, 1, 3},
                                             {4, 5, 6, 7}};

  /* Case of prism split into 8 smaller prisms */

  cs_lnum_t new_cell_id = c_o2n_idx[cell_id];

  cs_lnum_t quad_vtx[3][9];
  cs_lnum_t tria_vtx[2][6];

  /* Loop on faces; face->sub-face numbering relations based
     on CS_REFINE_TRIA and CS_REFINE_QUAD in _subdivide_face */

  _subdivide_cell_quad_faces(m,
                             0, 3, /* quad range */
                             cell_id,
                             n_b_f_ini,
                             c_o2n_idx,
                             i_face_o2n_idx,
                             b_face_o2n_idx,
                             c2f,
                             c2f2v_start,
                             c_id_shift,
                             quad_vtx);

  _subdivide_cell_tria_faces(m,
                             3, 5, /* tria range */
                             cell_id,
                             n_b_f_ini,
                             c_o2n_idx,
                             i_face_o2n_idx,
                             b_face_o2n_idx,
                             c2f,
                             c2f2v_start,
                             c_id_shift,
                             tria_vtx);

  cs_lnum_t v_ids[4];

  for (cs_lnum_t i = 0; i < 3; i++) {
    v_ids[0] = quad_vtx[(i+2)%3][4];
    v_ids[1] = quad_vtx[(i+2)%3][8];
    v_ids[2] = quad_vtx[i][8];
    v_ids[3] = quad_vtx[i][4];
    _add_interior_face_quad(m,
                            new_cell_id+i,
                            new_cell_id+3,
                            c_f_n_idx[cell_id] + i,
                            v_ids,
                            c_f_n_idx + cell_id);
  }

  for (cs_lnum_t i = 0; i < 3; i++) {
    v_ids[0] = quad_vtx[i][6];
    v_ids[1] = quad_vtx[i][8];
    v_ids[2] = quad_vtx[(i+2)%3][8];
    v_ids[3] = quad_vtx[(i+2)%3][6];
    _add_interior_face_quad(m,
                            new_cell_id+4+i,
                            new_cell_id+7,
                            c_f_n_idx[cell_id] + 3 + i,
                            v_ids,
                            c_f_n_idx + cell_id);
  }

  for (cs_lnum_t i = 0; i < 3; i++) {
    v_ids[0] = quad_vtx[i][7];
    v_ids[1] = quad_vtx[i][8];
    v_ids[2] = quad_vtx[(i+2)%3][8];
    _add_interior_face_tria(m,
                            new_cell_id+i,
                            new_cell_id+4+i,
                            c_f_n_idx[cell_id] + 6 + i,
                            v_ids,
                            c_f_n_idx + cell_id);
  }

  v_ids[0] = quad_vtx[0][8];
  v_ids[1] = quad_vtx[1][8];
  v_ids[2] = quad_vtx[2][8];
  _add_interior_face_tria(m,
                          new_cell_id+3,
                          new_cell_id+7,
                          c_f_n_idx[cell_id] + 9,
                          v_ids,
                          c_f_n_idx + cell_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Subdivide a given hexahedron, updating faces.
 *
 * The cells->faces adjacency is updated as faces are locally renumbered
 * in canonical order for known cell types.
 *
 * \param[in, out]  m             pointer to mesh structure
 * \param[in]       cell_id       (old) cell id (0 to n-1)
 * \param[in]       n_b_f_ini     old number of boundary faces
 * \param[in]       c_o2n_idx     old to new cells index
 * \param[in]       i_face_o2n_idx  old to new interior faces index
 * \param[in]       b_face_o2n_idx  old to new boundary faces index
 * \param[in]       c2f           cells->faces adjacency (boundary faces first)
 * \param[in]       c2f2v_start   start id for adjacent face vertices
 *                                definitions (initialized to 0)
 * \param[in]       c_v_idx       for each cell, start index of added vertices
 * \param[in]       c_f_n_idx     cells to new faces index
 */
/*----------------------------------------------------------------------------*/

static void
_subdivide_cell_hexa(const cs_mesh_t              *m,
                     cs_lnum_t                     cell_id,
                     cs_lnum_t                     n_b_f_ini,
                     const cs_lnum_t               c_o2n_idx[],
                     const cs_lnum_t               i_face_o2n_idx[],
                     const cs_lnum_t               b_face_o2n_idx[],
                     cs_adjacency_t               *c2f,
                     const cs_lnum_t               c2f2v_start[],
                     const cs_lnum_t               c_v_idx[],
                     const cs_lnum_t               c_f_n_idx[])
{
  static const cs_lnum_t c_id_shift[6][4] = {{0, 3, 2, 1},
                                             {0, 1, 5, 4},
                                             {1, 2, 6, 5},
                                             {2, 3, 7, 6},
                                             {3, 0, 4, 7},
                                             {4, 5, 6, 7}};

  assert(c2f->idx[cell_id+1] - c2f->idx[cell_id] == 6);

  /* Case of hexahedra split into 8 smaller hexahedra */

  const cs_lnum_t c_vtx_id = c_v_idx[cell_id];

  cs_lnum_t new_cell_id = c_o2n_idx[cell_id];

  cs_lnum_t quad_vtx[6][9];

  /* Loop on faces; face->sub-face numbering relation based
     on CS_REFINE_QUAD in _subdivide_face */

  _subdivide_cell_quad_faces(m,
                             0, 6, /* quad range */
                             cell_id,
                             n_b_f_ini,
                             c_o2n_idx,
                             i_face_o2n_idx,
                             b_face_o2n_idx,
                             c2f,
                             c2f2v_start,
                             c_id_shift,
                             quad_vtx);

  cs_lnum_t v_ids[4];

  v_ids[0] = quad_vtx[0][7];
  v_ids[1] = quad_vtx[0][8];
  v_ids[2] = c_vtx_id;
  v_ids[3] = quad_vtx[1][8];
  _add_interior_face_quad(m,
                          new_cell_id,
                          new_cell_id+1,
                          c_f_n_idx[cell_id] + 0,
                          v_ids,
                          c_f_n_idx + cell_id);

  v_ids[0] = quad_vtx[5][4];
  v_ids[1] = quad_vtx[1][8];
  v_ids[2] = c_vtx_id;
  v_ids[3] = quad_vtx[5][8];
  _add_interior_face_quad(m,
                          new_cell_id+4,
                          new_cell_id+5,
                          c_f_n_idx[cell_id] + 1,
                          v_ids,
                          c_f_n_idx + cell_id);

  v_ids[0] = quad_vtx[5][6];
  v_ids[1] = quad_vtx[5][8];
  v_ids[2] = c_vtx_id;
  v_ids[3] = quad_vtx[3][8];
  _add_interior_face_quad(m,
                          new_cell_id+7,
                          new_cell_id+6,
                          c_f_n_idx[cell_id] + 2,
                          v_ids,
                          c_f_n_idx + cell_id);

  v_ids[0] = quad_vtx[3][4];
  v_ids[1] = quad_vtx[3][8];
  v_ids[2] = c_vtx_id;
  v_ids[3] = quad_vtx[0][8];
  _add_interior_face_quad(m,
                          new_cell_id+3,
                          new_cell_id+2,
                          c_f_n_idx[cell_id] + 3,
                          v_ids,
                          c_f_n_idx + cell_id);

  v_ids[0] = quad_vtx[4][4];
  v_ids[1] = quad_vtx[4][8];
  v_ids[2] = c_vtx_id;
  v_ids[3] = quad_vtx[0][8];
  _add_interior_face_quad(m,
                          new_cell_id,
                          new_cell_id+3,
                          c_f_n_idx[cell_id] + 4,
                          v_ids,
                          c_f_n_idx + cell_id);

  v_ids[0] = quad_vtx[2][4];
  v_ids[1] = quad_vtx[0][8];
  v_ids[2] = c_vtx_id;
  v_ids[3] = quad_vtx[2][8];
  _add_interior_face_quad(m,
                          new_cell_id+1,
                          new_cell_id+2,
                          c_f_n_idx[cell_id] + 5,
                          v_ids,
                          c_f_n_idx + cell_id);

  v_ids[0] = quad_vtx[5][5];
  v_ids[1] = quad_vtx[2][8];
  v_ids[2] = c_vtx_id;
  v_ids[3] = quad_vtx[5][8];
  _add_interior_face_quad(m,
                          new_cell_id+5,
                          new_cell_id+6,
                          c_f_n_idx[cell_id] + 6,
                          v_ids,
                          c_f_n_idx + cell_id);

  v_ids[0] = quad_vtx[5][7];
  v_ids[1] = quad_vtx[5][8];
  v_ids[2] = c_vtx_id;
  v_ids[3] = quad_vtx[4][8];
  _add_interior_face_quad(m,
                          new_cell_id+4,
                          new_cell_id+7,
                          c_f_n_idx[cell_id] + 7,
                          v_ids,
                          c_f_n_idx + cell_id);

  v_ids[0] = quad_vtx[4][5];
  v_ids[1] = quad_vtx[1][8];
  v_ids[2] = c_vtx_id;
  v_ids[3] = quad_vtx[4][8];
  _add_interior_face_quad(m,
                          new_cell_id+0,
                          new_cell_id+4,
                          c_f_n_idx[cell_id] + 8,
                          v_ids,
                          c_f_n_idx + cell_id);


  v_ids[0] = quad_vtx[1][5];
  v_ids[1] = quad_vtx[2][8];
  v_ids[2] = c_vtx_id;
  v_ids[3] = quad_vtx[1][8];
  _add_interior_face_quad(m,
                          new_cell_id+1,
                          new_cell_id+5,
                          c_f_n_idx[cell_id] + 9,
                          v_ids,
                          c_f_n_idx + cell_id);

  v_ids[0] = quad_vtx[2][5];
  v_ids[1] = quad_vtx[3][8];
  v_ids[2] = c_vtx_id;
  v_ids[3] = quad_vtx[2][8];
  _add_interior_face_quad(m,
                          new_cell_id+2,
                          new_cell_id+6,
                          c_f_n_idx[cell_id] + 10,
                          v_ids,
                          c_f_n_idx + cell_id);

  v_ids[0] = quad_vtx[3][5];
  v_ids[1] = quad_vtx[4][8];
  v_ids[2] = c_vtx_id;
  v_ids[3] = quad_vtx[3][8];
  _add_interior_face_quad(m,
                          new_cell_id+3,
                          new_cell_id+7,
                          c_f_n_idx[cell_id] + 11,
                          v_ids,
                          c_f_n_idx + cell_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Pre-refine polyhedral cell.
 *
 * Here, the actual refinement is determined, but stored to temporary
 * structures, as the total sizes and associated indexes for the refined
 * mesh are not known yet.
 *
 * \param[in]       m               pointer to mesh structure
 * \param[in]       c2f             cells->faces adjacency (boundary faces
 *                                  first, interior face ids shifted)
 * \param[in, out]  crh             cell refinement helper
 * \param[in]       c_id            cell id
 * \param[in]       n_b_f_ini       old number of boundary faces
 * \param[in]       c_r_flag        cell refinement type flag
 * \param[in]       b_face_o2n_idx  boundary subface index
 * \param[in]       i_face_o2n_idx  interior subface index
 * \param[in]       c_v_idx         for each cell, start index of added vertices
 * \param[in]       f_r_level_o  old faces refinement level
 * \param[out]      c_n_s_cells     number of sub-cells
 * \param[out]      c_n_i_faces     number of added interior faces
 * \param[out]      c_i_faces_size  size of added interior faces connectivity
 * \param[out]      e_f2c_pre       temp. exterior faces to cell connect
 * \param[out]      i_f2c_pre       temp. cell to interior faces connect
 * \param[out]      i_f2v_pre       temp interior faces to vertices connect.
 */
/*----------------------------------------------------------------------------*/

static void
_pre_refine_poly_cell(const cs_mesh_t            *m,
                      const cs_adjacency_t       *c2f,
                      cs_cell_refine_helper      *crh,
                      cs_lnum_t                   c_id,
                      cs_lnum_t                   n_b_f_ini,
                      cs_mesh_refine_type_t       c_r_flag,
                      const cs_lnum_t             b_face_o2n_idx[],
                      const cs_lnum_t             i_face_o2n_idx[],
                      const cs_lnum_t             c_v_idx[],
                      const char                  f_r_level_o[],
                      cs_lnum_t                  &c_n_s_cells,
                      cs_lnum_t                  &c_n_i_faces,
                      cs_lnum_t                  &c_i_faces_size,
                      int16_t                     e_f2c_pre[],
                      int16_t                     i_f2c_pre[][2],
                      cs_lnum_t                   i_f2v_pre[])
{
  crh->reset_for_new_cell();
  crh->set_center_vtx(c_v_idx[c_id]);

  const cs_lnum_t s_id_c = c2f->idx[c_id];
  const cs_lnum_t e_id_c = c2f->idx[c_id+1];

  const cs_lnum_t n_cell_faces = e_id_c - s_id_c;

  cs_lnum_t e_f2c_count = 0;  // counter for number of exterior faces

  for (cs_lnum_t i = 0; i < n_cell_faces; i++) {
    cs_lnum_t f_id_o = c2f->ids[s_id_c + i];
    int sgn = c2f->sgn[s_id_c + i];
    char f_r_level = f_r_level_o[f_id_o];

    const cs_lnum_t *o2n_idx, *vtx_idx, *vtx_lst;
    if (f_id_o < n_b_f_ini) {
      o2n_idx = b_face_o2n_idx;
      vtx_idx = m->b_face_vtx_idx;
      vtx_lst = m->b_face_vtx_lst;
    }
    else {
      f_id_o -= n_b_f_ini;
      o2n_idx = i_face_o2n_idx;
      vtx_lst = m->i_face_vtx_lst;
      vtx_idx = m->i_face_vtx_idx;
    }
    cs_lnum_t s_f_idx = o2n_idx[f_id_o];
    cs_lnum_t e_f_idx = o2n_idx[f_id_o+1];

    if (e_f_idx - s_f_idx > 1)
      f_r_level += 1;

    for (cs_lnum_t j = s_f_idx; j < e_f_idx; j++) {

      cs_lnum_t f_id_n = j;
      cs_lnum_t v_s_id = vtx_idx[f_id_n];
      cs_lnum_t v_e_id = vtx_idx[f_id_n+1];

      crh->add_face(sgn,
                    f_r_level,
                    v_e_id - v_s_id,
                    vtx_lst + v_s_id);

      e_f2c_count++;

    }

  }

  if (   c_r_flag == CS_REFINE_POLYHEDRON
      || c_r_flag == CS_REFINE_POLYHEDRON_P)
    crh->refine_cell_with_center_vertex(c_r_flag,
                                        m->vtx_r_gen,
                                        c_n_s_cells,
                                        c_n_i_faces,
                                        c_i_faces_size,
                                        e_f2c_pre,
                                        i_f2c_pre,
                                        i_f2v_pre);
  else
    cs_assert(0);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Refine polyhedral cell.
 *
 * Here, the refinement stored in temporary structures, is applied.
 *
 * \param[in, out]  m             pointer to mesh structure
 * \param[in]       cell_id       (old) cell id (0 to n-1)
 * \param[in]       n_b_f_ini     old number of boundary faces
 * \param[in]       c_o2n_idx     old to new cells index
 * \param[in]       i_face_o2n_idx  old to new interior faces index
 * \param[in]       b_face_o2n_idx  old to new boundary faces index
 * \param[in]       c2f           cells->faces adjacency (boundary faces first)
 * \param[in]       c_v_idx       for each cell, start index of added vertices
 * \param[in]       c_f_n_idx     cells to new faces index
 * \param[in]       e_f2c_pre     temp. exterior faces to cell connect
 * \param[in]       i_f2c_pre     temp. interior faces to cell connect
 * \param[in]       i_f2v_pre     temp interior faces to vertices connect.
 */
/*----------------------------------------------------------------------------*/

static void
_refine_poly_cell(const cs_mesh_t            *m,
                  cs_lnum_t                   cell_id,
                  cs_lnum_t                   n_b_f_ini,
                  const cs_lnum_t             c_o2n_idx[],
                  const cs_lnum_t             i_face_o2n_idx[],
                  const cs_lnum_t             b_face_o2n_idx[],
                  cs_adjacency_t             *c2f,
                  const cs_lnum_t             c_f_n_idx[],
                  cs_lnum_t                   c_n_e_faces,
                  const int16_t               e_f2c_pre[],
                  const int16_t               i_f2c_pre[][2],
                  const cs_lnum_t             i_f2v_pre[])
{
  const cs_lnum_t s_id_c = c2f->idx[cell_id];
  const cs_lnum_t e_id_c = c2f->idx[cell_id+1];

  const cs_lnum_t n_cell_faces = e_id_c - s_id_c;

  const cs_lnum_t c_id_n = c_o2n_idx[cell_id];

  /* Update exterior faces to cell connectivity
     ------------------------------------------ */

  cs_lnum_t e_f2c_count = 0;  // counter for number of exterior faces

  for (cs_lnum_t i = 0; i < n_cell_faces; i++) {
    cs_lnum_t f_id_o = c2f->ids[s_id_c + i];
    int sgn = c2f->sgn[s_id_c + i];

    if (f_id_o < n_b_f_ini) {
      cs_lnum_t s_f_idx = b_face_o2n_idx[f_id_o];
      cs_lnum_t e_f_idx = b_face_o2n_idx[f_id_o+1];
      for (cs_lnum_t f_id = s_f_idx; f_id < e_f_idx; f_id++) {
        m->b_face_cells[f_id] = c_id_n + (cs_lnum_t)(e_f2c_pre[e_f2c_count]);
        e_f2c_count++;
      }
    }
    else {
      f_id_o -= n_b_f_ini;
      cs_lnum_t s_f_idx = i_face_o2n_idx[f_id_o];
      cs_lnum_t e_f_idx = i_face_o2n_idx[f_id_o+1];
      cs_lnum_t j = (sgn > 0) ? 0 : 1;
      for (cs_lnum_t f_id = s_f_idx; f_id < e_f_idx; f_id++) {
        m->i_face_cells[f_id][j] = c_id_n + (cs_lnum_t)(e_f2c_pre[e_f2c_count]);
        e_f2c_count++;
      }
    }

  }

  assert(e_f2c_count == c_n_e_faces);

  /* Now add new interior faces to cell connectivity
     ----------------------------------------------- */

  cs_lnum_t *i_face_vtx_lst = m->i_face_vtx_lst;
  cs_lnum_t j = 0;

  cs_lnum_t cf_s_id = c_f_n_idx[cell_id];
  cs_lnum_t cf_e_id = c_f_n_idx[cell_id+1];

  for (cs_lnum_t f_id = cf_s_id, i_f_idx = 0;
       f_id < cf_e_id;
       f_id++, i_f_idx++) {

    cs_lnum_t s_id = m->i_face_vtx_idx[f_id];

    for (cs_lnum_t i = s_id; true; i++) {
      cs_lnum_t k = i_f2v_pre[j++];
      if (k > -1)
        i_face_vtx_lst[i] = k;
      else {
#if defined(DEBUG) && !defined(NDEBUG)
        if (m->i_face_vtx_idx[f_id+1] == -1)
          m->i_face_vtx_idx[f_id+1] = i; // end of face
        else {
          assert(m->i_face_vtx_idx[f_id+1] == i);
        }
#else
        m->i_face_vtx_idx[f_id+1] = i; // end of face
#endif
        break;
      }
    }

    for (cs_lnum_t k = 0; k < 2; k++)
      m->i_face_cells[f_id][k] = c_id_n + i_f2c_pre[i_f_idx][k];

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Subdivide cells, updating faces.
 *
 * The cells->faces adjacency is updated as faces are locally renumbered
 * in canonical order for known cell types.
 *
 * \param[in]  m                  pointer to mesh structure
 * \param[in]  n_c_ini            old number of cells
 * \param[in]  n_b_f_ini          old number of boundary faces
 * \param[in]  c_o2n_idx          old to new cells index
 * \param[in]  i_face_o2n_idx     old to new interior faces index
 * \param[in]  b_face_o2n_idx     old to new boundary faces index
 * \param[in]  c2f                cells->faces adjacency (boundary faces first)
 * \param[in]  c2f2v_start        start id for adjacent face vertices
 *                                definitions (initialized to 0)
 * \param[in]  c_v_idx            for each cell, start index of added vertices
 * \param[in]  c_f_n_idx          cells to new faces index (count in),
 *                                or nullptr
 * \param[in]  c_r_flag           cell refinement type
 * \param[in]  c_r_level          cell refinement level
 * \param[in]  e_f2c_pre_idx      temp. exterior faces to cell index
 * \param[in]  e_f2c_pre          temp. exterior faces to cell connect
 * \param[in]  i_f2c_pre_idx      temp. interior faces to cells index
 * \param[in]  i_f2c_pre          temp. interior faces to cells connect
 * \param[in]  i_f2v_pre_idx      temp interior faces to vertices index
 * \param[in]  i_f2v_pre          temp interior faces to vertices connect.
 */
/*----------------------------------------------------------------------------*/

inline static void
_subdivide_cells(const cs_mesh_t              *m,
                 cs_lnum_t                     n_c_ini,
                 cs_lnum_t                     n_b_f_ini,
                 const cs_lnum_t               c_o2n_idx[],
                 const cs_lnum_t               i_face_o2n_idx[],
                 const cs_lnum_t               b_face_o2n_idx[],
                 cs_adjacency_t               *c2f,
                 cs_lnum_t                     c2f2v_start[],
                 const cs_lnum_t               c_v_idx[],
                 const cs_lnum_t               c_f_n_idx[],
                 const cs_mesh_refine_type_t   c_r_flag[],
                 const char                    c_r_level[],
                 cs_lnum_t                    *c_poly_s_idx,
                 const cs_lnum_t              *e_f2c_pre_idx,
                 const int16_t                *e_f2c_pre,
                 const cs_lnum_t              *i_f2c_pre_idx,
                 const int16_t                *i_f2c_pre,
                 const cs_lnum_t              *i_f2v_pre_idx,
                 const cs_lnum_t              *i_f2v_pre)
{
  subdivide_cell_t  *c_r_func[CS_REFINE_N_TYPES];

  for (int i = 0; i < CS_REFINE_N_TYPES; i++)
    c_r_func[i] = nullptr;

  c_r_func[CS_REFINE_TETRA] = _subdivide_cell_tetra;
  c_r_func[CS_REFINE_TETRA_H] = _subdivide_cell_tetra_h;
  c_r_func[CS_REFINE_PYRAM] = _subdivide_cell_pyram;
  c_r_func[CS_REFINE_PRISM] = _subdivide_cell_prism;
  c_r_func[CS_REFINE_HEXA] = _subdivide_cell_hexa;

  using c2f_pre_t = int16_t[2];
  const c2f_pre_t *_i_f2c_pre = reinterpret_cast<const c2f_pre_t *>(i_f2c_pre);

# pragma omp for schedule(dynamic, CS_CL_SIZE)
  for (cs_lnum_t c_id = 0; c_id < n_c_ini; c_id++) {

    subdivide_cell_t  *_c_r_func = c_r_func[c_r_flag[c_id]];
    if (_c_r_func != nullptr) {
      _c_r_func(m,
                c_id,
                n_b_f_ini,
                c_o2n_idx,
                i_face_o2n_idx,
                b_face_o2n_idx,
                c2f,
                c2f2v_start,
                c_v_idx,
                c_f_n_idx);

      cs_lnum_t e_id = c_f_n_idx[c_id+1];
      for (cs_lnum_t f_id = c_f_n_idx[c_id]; f_id < e_id; f_id++)
        m->i_face_r_gen[f_id] = c_r_level[c_id] + 1;
    }
    else if (c_r_flag[c_id] != CS_REFINE_NONE) {
      cs_lnum_t c_idx = c_poly_s_idx[c_id];

      if (c_poly_s_idx[c_id+1] != c_idx+1) { // true if helper is used
        bft_error(__FILE__, __LINE__, 0,
                  _("Unhandled cell refinement type: %d"), c_r_flag[c_id]);
        continue;
      }

      cs_lnum_t c_n_e_faces = e_f2c_pre_idx[c_idx+1] - e_f2c_pre_idx[c_idx];
      assert(   i_f2c_pre_idx[c_idx+1] - i_f2c_pre_idx[c_idx]
             >= c_f_n_idx[c_id+1] - c_f_n_idx[c_id]);

      _refine_poly_cell(m,
                        c_id,
                        n_b_f_ini,
                        c_o2n_idx,
                        i_face_o2n_idx,
                        b_face_o2n_idx,
                        c2f,
                        c_f_n_idx,
                        c_n_e_faces,
                        e_f2c_pre + e_f2c_pre_idx[c_idx],
                        _i_f2c_pre + i_f2c_pre_idx[c_idx],
                        i_f2v_pre + i_f2v_pre_idx[c_idx]);

      cs_lnum_t e_id = c_f_n_idx[c_id+1];
      for (cs_lnum_t f_id = c_f_n_idx[c_id]; f_id < e_id; f_id++)
        m->i_face_r_gen[f_id] = c_r_level[c_id] + 1;
    }

  }

  for (cs_lnum_t i = 0; i < m->n_i_faces+1; i++) {
    if (m->i_face_vtx_idx[i] < -1) printf("a\n");
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update connectivity for a given faces set
 *
 * This function also transforms some counts to indexes
 * (grouping some conversions so as to reduce the number of required loops).
 *
 * \param[in]       v2v                vertex->vertex adjacency
 * \param[in]       vtx_coords         mesh vertex coordinates
 * \param[in]       vtx_r_gen          vertices refinement generation
 * \param[in]       e_v_idx            for each edge, start index of
 *                                     added vertices
 * \param[in]       f_v_idx            for each face, start index of
 *                                     added vertices
 * \param[in]       n_faces            number of faces
 * \param[in]       f_r_flag           face refinement type flag
 * \param[in]       f_r_level_o        old faces refinement level
 * \param[in, out]  f_o2n_idx          old to new faces index
 * \param[in, out]  f_o2n_connect_idx  old to new faces connectivity index
 *                                     (count in)
 * \param[in, out]  f2v_idx            face->vertices index
 * \param[in, out]  f2v_lst            face->vertices connectivity
 *
 * \return  new number of faces
 */
/*----------------------------------------------------------------------------*/

static cs_lnum_t
_update_face_connectivity(const cs_adjacency_t         *v2v,
                          const cs_real_t               vtx_coords[],
                          const char                    vtx_r_gen[],
                          const cs_lnum_t               e_v_idx[],
                          const cs_lnum_t               f_v_idx[],
                          cs_lnum_t                     n_faces,
                          const cs_mesh_refine_type_t   f_r_flag[],
                          const char                    f_r_level_o[],
                          const cs_lnum_t               f_o2n_idx[],
                          cs_lnum_t                     f_o2n_connect_idx[],
                          cs_lnum_t                    *f2v_idx[],
                          cs_lnum_t                    *f2v_lst[])
{
  cs_lnum_t *f2v_lst_o = *f2v_lst, *f2v_idx_o = *f2v_idx;
  cs_lnum_t n_faces_new = 0;

  /* Convert counts to indexes
     ------------------------- */

  n_faces_new = f_o2n_idx[n_faces];

  cs_lnum_t idx_size_new = n_faces_new + 1;

  cs_lnum_t *_f2v_idx;
  CS_MALLOC(_f2v_idx, idx_size_new, cs_lnum_t);

  /* Partial index creation (allows threading following part) */

  _f2v_idx[0] = 0;
  f_o2n_connect_idx[0] = 0;
  for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++) {
    cs_lnum_t s_id = f_o2n_idx[f_id];
    cs_lnum_t e_id = f_o2n_idx[f_id+1];
    cs_lnum_t connect_size = f_o2n_connect_idx[f_id+1];
    f_o2n_connect_idx[f_id+1] += f_o2n_connect_idx[f_id];
    _f2v_idx[e_id] = _f2v_idx[s_id] + connect_size;
  }
  cs_lnum_t connect_size_new = f_o2n_connect_idx[n_faces];

  cs_lnum_t *_f2v_lst;
  CS_MALLOC(_f2v_lst, connect_size_new, cs_lnum_t);

  /* Build new faces */

# pragma omp parallel  if(n_faces > CS_THR_MIN)
  {
    fvm_triangulate_state_t *t_state = fvm_triangulate_state_create(8);

#   pragma omp for
    for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++) {

      cs_lnum_t s_id = f2v_idx_o[f_id];
      cs_lnum_t n_fv = f2v_idx_o[f_id+1] - s_id;

      _subdivide_face(f_id,
                      n_fv,
                      f_r_flag[f_id],
                      f_r_level_o[f_id],
                      f2v_lst_o + s_id,
                      v2v,
                      e_v_idx,
                      f_v_idx,
                      vtx_coords,
                      vtx_r_gen,
                      t_state,
                      _f2v_idx + f_o2n_idx[f_id],
                      _f2v_lst + f_o2n_connect_idx[f_id]);

    }

    t_state = fvm_triangulate_state_destroy(t_state);
  }

  CS_FREE(f2v_idx_o);
  CS_FREE(f2v_lst_o);

  *f2v_idx = _f2v_idx;
  *f2v_lst = _f2v_lst;

  return n_faces_new;
}

/*----------------------------------------------------------------------------
 * Update a global numbering array in case of entity renumbering
 *
 * parameters:
 *   n_old      <-- old number of elements
 *   n_g_old    <-- old global number of elements
 *   o2n_idx    <-- old to new index
 *   global_num <-> global numbering (allocated if initially nullptr)
 *
 * returns:
 *   new global number of elements
 *----------------------------------------------------------------------------*/

static cs_gnum_t
_o2n_idx_update_global_num(cs_lnum_t          n_old,
                           cs_gnum_t          n_g_old,
                           const cs_lnum_t    o2n_idx[],
                           cs_gnum_t        **global_num)
{
  cs_gnum_t n_g_new = o2n_idx[n_old];

  if (cs_glob_n_ranks == 1 && *global_num == nullptr)
    return n_g_new;

  fvm_io_num_t *o_io_num
    = fvm_io_num_create_shared(*global_num, n_g_old, n_old);

  cs_lnum_t *n_sub;
  CS_MALLOC(n_sub, n_old, cs_lnum_t);
  cs_lnum_t *restrict _n_sub = n_sub;
  for (cs_lnum_t i = 0; i < n_old; i++)
    _n_sub[i] = o2n_idx[i+1] - o2n_idx[i];
  _n_sub = nullptr;

  fvm_io_num_t *n_io_num
    = fvm_io_num_create_from_sub(o_io_num, n_sub);

  o_io_num = fvm_io_num_destroy(o_io_num);

  CS_FREE(n_sub);
  CS_FREE(*global_num);

  *global_num = fvm_io_num_transfer_global_num(n_io_num);

  n_g_new = fvm_io_num_get_global_count(n_io_num);

  n_io_num = fvm_io_num_destroy(n_io_num);

  return n_g_new;
}

/*----------------------------------------------------------------------------
 * Complete a global numbering array in case of entity renumbering
 *
 * parameters:
 *   n_old            <-- old number of referencing elements
 *   n_g_old          <-- old global number of elements
 *   global_num_shift <-- global numbering shift for added elements
 *   o2n_idx          <-- old to new index
 *   old_global_num   <-- old global numbering
 *   new_global_num   --> new global numbering
 *
 * returns:
 *   new global number of elements
 *----------------------------------------------------------------------------*/

static cs_gnum_t
_o2n_idx_complete_global_num(cs_lnum_t          n_old,
                             cs_gnum_t          n_g_old,
                             cs_gnum_t          global_num_shift,
                             const cs_lnum_t    o2n_idx[],
                             const cs_gnum_t   *old_global_num,
                             cs_gnum_t         *new_global_num)
{
  cs_gnum_t n_g_new = global_num_shift;
  cs_gnum_t *_old_global_num = nullptr;
  const cs_gnum_t *_new_global_num = nullptr;

  const cs_lnum_t n_new = o2n_idx[n_old] - o2n_idx[0];

  if (old_global_num == nullptr) {
    CS_MALLOC(_old_global_num, n_old, cs_gnum_t);
    for (cs_lnum_t i = 0; i < n_old; i++)
      _old_global_num[i] = i+1;
    old_global_num = _old_global_num;
  }

  fvm_io_num_t *o_io_num
    = fvm_io_num_create_shared(old_global_num, n_g_old, n_old);

  cs_lnum_t *n_sub;
  CS_MALLOC(n_sub, n_old, cs_lnum_t);
  cs_lnum_t *restrict _n_sub = n_sub;
  for (cs_lnum_t i = 0; i < n_old; i++)
    _n_sub[i] = o2n_idx[i+1] - o2n_idx[i];
  _n_sub = nullptr;

  fvm_io_num_t *n_io_num
    = fvm_io_num_create_from_sub(o_io_num, n_sub);

  o_io_num = fvm_io_num_destroy(o_io_num);

  CS_FREE(n_sub);
  CS_FREE(_old_global_num);

  n_g_new += fvm_io_num_get_global_count(n_io_num);

  _new_global_num = fvm_io_num_get_global_num(n_io_num);
  for (cs_lnum_t i = 0; i < n_new; i++)
    new_global_num[i] = _new_global_num[i] + global_num_shift;

  n_io_num = fvm_io_num_destroy(n_io_num);

  return n_g_new;
}

/*----------------------------------------------------------------------------
 * Update arrays related to cells, as well as cell counts.
 *
 * Also updates face references to the first cell of a given subset
 * (which is sufficient for non-refined cells, and allows using local
 * shifts for refined cells).
 *
 * Must be called after face array updates.
 *
 * parameters:
 *   m         <-> pointer to global mesh structure
 *   n_c_ini   <-- old number of cells
 *   o2n_idx   <-- old to new index
 *----------------------------------------------------------------------------*/

static void
_o2n_idx_update_cell_arrays(cs_mesh_t        *m,
                            cs_lnum_t         n_c_ini,
                            const cs_lnum_t  *o2n_idx)
{
  const cs_lnum_t n_new = o2n_idx[n_c_ini];

  /* Allocate new arrays */

  int *cell_family;
  CS_MALLOC(cell_family, n_new, int);

  for (cs_lnum_t o_id = 0; o_id < n_c_ini; o_id++) {
    for (cs_lnum_t n_id = o2n_idx[o_id]; n_id < o2n_idx[o_id+1]; n_id++) {
      /* update family */
      cell_family[n_id] = m->cell_family[o_id];
    }
  }

  CS_FREE(m->cell_family);
  m->cell_family = cell_family;

  /* Update global numbering */

  m->n_g_cells
    = _o2n_idx_update_global_num(n_c_ini, m->n_g_cells,
                                 o2n_idx, &(m->global_cell_num));
  m->n_cells = n_new;
  m->n_cells_with_ghosts = n_new;

  /* Update face references to first cell of subset */

  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

# pragma omp for schedule(dynamic, CS_CL_SIZE)
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t i0 = m->i_face_cells[f_id][0];
    cs_lnum_t i1 = m->i_face_cells[f_id][1];
    if (i0 >= n_c_ini)
      m->i_face_cells[f_id][0] = -1;
    else if (i0 > -1)
      m->i_face_cells[f_id][0] = o2n_idx[i0];
    if (i1 >= n_c_ini)
      m->i_face_cells[f_id][1] = -1;
    else if (i1 > -1)
      m->i_face_cells[f_id][1] = o2n_idx[i1];
  }

# pragma omp for schedule(dynamic, CS_CL_SIZE)
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t i = m->b_face_cells[f_id];
    if (i > -1)
      m->b_face_cells[f_id] = o2n_idx[i];
  }
}

/*----------------------------------------------------------------------------
 * Update arrays related to interior faces, as well as face counts.
 *
 * parameters:
 *   m       <-> pointer to global mesh structure
 *   o2n_idx <-- old to new interior faces index
 *----------------------------------------------------------------------------*/

static void
_o2n_idx_update_i_face_arrays(cs_mesh_t        *m,
                              const cs_lnum_t   o2n_idx[])
{
  const cs_lnum_t n_old = m->n_i_faces;
  const cs_lnum_t n_new = o2n_idx[n_old];

  /* Allocate new arrays */

  cs_lnum_2_t *i_face_cells;
  int *i_face_family;
  char *i_face_r_gen;

  CS_MALLOC(i_face_cells, n_new, cs_lnum_2_t);
  CS_MALLOC(i_face_family, n_new, int);
  CS_MALLOC(i_face_r_gen, n_new, char);

  for (cs_lnum_t o_id = 0; o_id < n_old; o_id++) {
    for (cs_lnum_t n_id = o2n_idx[o_id]; n_id < o2n_idx[o_id+1]; n_id++) {
      /* update faces -> cells connectivity */
      i_face_cells[n_id][0] = m->i_face_cells[o_id][0];
      i_face_cells[n_id][1] = m->i_face_cells[o_id][1];
      /* update family */
      i_face_family[n_id] = m->i_face_family[o_id];
      /* update generation */
      i_face_r_gen[n_id] = m->i_face_r_gen[o_id];
    }
  }

  CS_FREE(m->i_face_r_gen);
  CS_FREE(m->i_face_family);
  CS_FREE(m->i_face_cells);
  m->i_face_r_gen = i_face_r_gen;
  m->i_face_cells = i_face_cells;
  m->i_face_family = i_face_family;

  /* Update global numbering */

  m->n_g_i_faces
    = _o2n_idx_update_global_num(n_old, m->n_g_i_faces,
                                 o2n_idx, &(m->global_i_face_num));

  m->n_i_faces = n_new;
  m->i_face_vtx_connect_size = m->i_face_vtx_idx[n_new];
}

/*----------------------------------------------------------------------------
 * Update arrays related to boundary faces, as well as face counts.
 *
 * parameters:
 *   m       <-> pointer to global mesh structure
 *   o2n_idx <-- old to new index
 *----------------------------------------------------------------------------*/

static void
_o2n_idx_update_b_face_arrays(cs_mesh_t        *m,
                              const cs_lnum_t  *o2n_idx)
{
  const cs_lnum_t n_old = m->n_b_faces;
  const cs_lnum_t n_new = o2n_idx[n_old];

  /* Allocate new arrays */

  cs_lnum_t *b_face_cells;
  int *b_face_family;
  char *b_face_r_c_idx;

  CS_MALLOC(b_face_cells, n_new, cs_lnum_t);
  CS_MALLOC(b_face_family, n_new, int);
  CS_MALLOC(b_face_r_c_idx, n_new, char);

  for (cs_lnum_t o_id = 0; o_id < n_old; o_id++) {
    for (cs_lnum_t n_id = o2n_idx[o_id]; n_id < o2n_idx[o_id+1]; n_id++) {
      /* update faces -> cells connectivity */
      b_face_cells[n_id] = m->b_face_cells[o_id];
      /* update family */
      b_face_family[n_id] = m->b_face_family[o_id];
      /* update local index of boundary faces */
      b_face_r_c_idx[n_id] = m->b_face_r_c_idx[o_id];
    }
  }

  CS_FREE(m->b_face_family);
  CS_FREE(m->b_face_cells);
  CS_FREE(m->b_face_r_c_idx);
  m->b_face_cells = b_face_cells;
  m->b_face_family = b_face_family;
  m->b_face_r_c_idx = b_face_r_c_idx;

  /* Update global numbering */

  m->n_g_b_faces
    = _o2n_idx_update_global_num(n_old, m->n_g_b_faces,
                                 o2n_idx, &(m->global_b_face_num));
  m->n_b_faces = n_new;
  m->b_face_vtx_connect_size = m->b_face_vtx_idx[n_new];
}

/*----------------------------------------------------------------------------
 * Update arrays related to interior faces, as well as face counts.
 *
 * parameters:
 *   m              <-> pointer to global mesh structure
 *   n_c_ini        <-- old number of cells
 *   c_i_face_idx   <-- cells to new faces index, or nullptr
 *   c_i_faces_size <-- connectivity size for faces of a given cell
 *----------------------------------------------------------------------------*/

static void
_o2n_idx_update_cell_i_face_arrays(cs_mesh_t        *m,
                                   cs_lnum_t         n_c_ini,
                                   const cs_lnum_t   c_i_face_idx[],
                                   const cs_lnum_t   c_i_faces_size[])
{
  const cs_lnum_t n_old = m->n_i_faces;
  const cs_lnum_t n_new = c_i_face_idx[n_c_ini];

  const int default_family_id = 1;

  /* Resize arrays */

  CS_REALLOC(m->i_face_vtx_idx, n_new+1, cs_lnum_t);

  CS_REALLOC(m->i_face_cells, n_new, cs_lnum_2_t);
  CS_REALLOC(m->i_face_family, n_new, int);
  CS_REALLOC(m->i_face_r_gen, n_new, char);

  cs_lnum_2_t *i_face_cells = m->i_face_cells;
  int *i_face_family = m->i_face_family;
  char *i_face_r_gen = m->i_face_r_gen;
  cs_lnum_t *i_face_vtx_idx = m->i_face_vtx_idx;

  for (cs_lnum_t n_id = n_old; n_id < n_new; n_id++) {
    /* initialize faces -> cells connectivity */
    i_face_cells[n_id][0] = -1;
    i_face_cells[n_id][1] = -1;
    /* initialize family */
    i_face_family[n_id] = default_family_id;
    /* initialize generation */
    i_face_r_gen[n_id] = 0;
  }

  for (cs_lnum_t c_id = 0; c_id < n_c_ini; c_id++) {
    cs_lnum_t s_id = c_i_face_idx[c_id];
    cs_lnum_t e_id = c_i_face_idx[c_id+1];
    cs_lnum_t connect_size = c_i_faces_size[c_id];
#if defined(DEBUG) && !defined(NDEBUG)
    for (cs_lnum_t i = s_id+1; i < e_id; i++)
      i_face_vtx_idx[i] = -1;
#endif
    i_face_vtx_idx[e_id] = i_face_vtx_idx[s_id] + connect_size;
  }

  m->i_face_vtx_connect_size = m->i_face_vtx_idx[n_new];
  CS_REALLOC(m->i_face_vtx_lst, m->i_face_vtx_connect_size, cs_lnum_t);

  /* Update global numbering */

  if (cs_glob_n_ranks > 1 || m->global_i_face_num != nullptr) {
    CS_REALLOC(m->global_i_face_num, n_new, cs_gnum_t);
    m->n_g_i_faces
      = _o2n_idx_complete_global_num(n_c_ini,
                                     m->n_g_cells,
                                     m->n_g_i_faces,
                                     c_i_face_idx,
                                     m->global_cell_num,
                                     m->global_i_face_num + n_old);
  }
  else
    m->n_g_i_faces = n_new;

  m->n_i_faces = n_new;
  m->i_face_vtx_connect_size = m->i_face_vtx_idx[n_new];
}

/*----------------------------------------------------------------------------
 * Free cell refinement helpers.
 *
 * parameters:
 *   n        <--  number of mesh refinement helpers
 *   helpers  <->  reference to helpers array
 *----------------------------------------------------------------------------*/

static void
_free_cell_refine_helpers(int                      n,
                          cs_cell_refine_helper  *&helpers)
{
  for (int i = 0; i < n; i++) {
    // Explicit destructor call to free associated arrays.
    cs_cell_refine_helper *crh = helpers + i;
    crh->~cs_cell_refine_helper();
  }

  CS_FREE(helpers);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Refine flagged mesh cells.
 *
 * \param[in, out]  m           mesh
 * \param[in]       conforming  if true, propagate refinement to ensure
 *                              subdivision is conforming
 * \param[in]       cell_flag   subdivision type for each cell
 *                              (0: none; 1: isotropic)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_refine_simple(cs_mesh_t  *m,
                      bool        conforming,
                      const int   cell_flag[])
{
  /* Timers:
     0: total
     1: update mesh structure (ghosts, ...)
     2: build edge and face indexing
     3: compute sub-face indexing
     4: build cell->faces connectivity and identify refined cell types
     5: build new vertices
     6: subdivide faces connectivity
     7: compute added interior faces indexing
     8: subdivide cells
  */

  cs_timer_counter_t  timers[9];
  for (int i = 0; i < 9; i++)
    CS_TIMER_COUNTER_INIT(timers[i]);

  cs_lnum_t n_v_ini = m->n_vertices;
  cs_lnum_t n_f_ini = m->n_b_faces + m->n_i_faces;
  cs_lnum_t n_c_ini = m->n_cells;

  cs_lnum_t n_b_f_ini = m->n_b_faces;

  bool check_convex = true;

  cs_timer_t t0 = cs_timer_time();

  /* Build ghosts in case they are not present */

  int mv_save = m->verbosity;
  m->verbosity = -1;

  if ((m->n_domains > 1 || m->n_init_perio > 0) && m->halo == nullptr) {
    cs_halo_type_t halo_type = m->halo_type;
    cs_mesh_builder_t *mb = (m == cs_glob_mesh) ? cs_glob_mesh_builder : nullptr;
    cs_mesh_init_halo(m, mb, halo_type, -1, true);
    cs_mesh_update_auxiliary(m);
  }

  m->verbosity = mv_save;

  if (m->verbosity > 0) {
    cs_log_printf(CS_LOG_DEFAULT, "\n");
    cs_log_separator(CS_LOG_DEFAULT);
    _print_mesh_counts(m, _("Mesh before refinement"));
  }

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(timers[1]), &t0, &t1);

  /* Compute some mesh quantities */

  /* Build vertex and face generation if not available yet */

  if (m->vtx_r_gen == nullptr) {
    CS_MALLOC(m->vtx_r_gen, m->n_vertices, char);
    for (cs_lnum_t i = 0; i < m->n_vertices; i++)
      m->vtx_r_gen[i] = 0;
  }

  if (m->i_face_r_gen == nullptr) {
    CS_MALLOC(m->i_face_r_gen, m->n_i_faces, char);
    for (cs_lnum_t i = 0; i < m->n_i_faces; i++)
      m->i_face_r_gen[i] = 0;
  }

  if (m->b_face_r_c_idx == nullptr) {
    CS_MALLOC(m->b_face_r_c_idx, m->n_b_faces, char);

    /* Build cell->faces connectivity to group associated faces */
    cs_adjacency_t *c2f = cs_mesh_adjacency_c2f_boundary(m);

    for (cs_lnum_t c_id = 0; c_id < m->n_cells; c_id++) {
      cs_lnum_t s_id = c2f->idx[c_id];
      cs_lnum_t e_id = c2f->idx[c_id+1];
      cs_lnum_t n_loc_b_faces = e_id - s_id;
      for (cs_lnum_t f_b_loc = 0; f_b_loc < n_loc_b_faces; f_b_loc++){
        cs_lnum_t f_id = c2f->ids[s_id + f_b_loc];
        m->b_face_r_c_idx[f_id] = f_b_loc;
      }
    }

    cs_adjacency_destroy(&c2f);
  }

  /* Determine current cell refinement level */

  char *c_r_level = nullptr;
  {
    CS_MALLOC(c_r_level, m->n_cells_with_ghosts, char);
    for (cs_lnum_t i = 0; i < m->n_cells_with_ghosts; i++)
      c_r_level[i] = 0;

    const cs_lnum_2_t *restrict i_face_cells = m->i_face_cells;

    for (cs_lnum_t f_id = 0; f_id < m->n_i_faces; f_id++) {
      for (cs_lnum_t i = 0; i < 2; i++) {
        cs_lnum_t c_id = i_face_cells[f_id][i];
        if (m->i_face_r_gen[f_id] > c_r_level[c_id])
          c_r_level[c_id] = m->i_face_r_gen[f_id];
      }
    }

    if (m->halo != nullptr)
      cs_halo_sync_untyped(m->halo,
                           CS_HALO_STANDARD,
                           1,
                           c_r_level);
  }

  /* Number of added vertices for edges, faces, and cells */

  cs_lnum_t n_add_vtx[3] = {0, 0, 0};

  /* Refinement type flags for tracking and propagating refinement
     for cells and faces */

  cs_mesh_refine_type_t *c_r_flag, *f_r_flag;

  CS_MALLOC(c_r_flag, m->n_cells_with_ghosts, cs_mesh_refine_type_t);

  for (cs_lnum_t i = 0; i < m->n_cells; i++)
    c_r_flag[i] = static_cast<cs_mesh_refine_type_t>(cs::max(0, cell_flag[i]));
  for (cs_lnum_t i = m->n_cells; i < m->n_cells_with_ghosts; i++)
    c_r_flag[i] = CS_REFINE_NONE;

  CS_MALLOC(f_r_flag, n_f_ini, cs_mesh_refine_type_t);
  for (cs_lnum_t i = 0; i < n_f_ini; i++)
    f_r_flag[i] = CS_REFINE_NONE;

  cs_lnum_t n_refined_cells = 0;
  for (cs_lnum_t i = 0; i < m->n_cells; i++) {
    if (c_r_flag[i] > CS_REFINE_NONE)
      n_refined_cells += 1;
  }

  cs_lnum_t *refined_cell_id = nullptr;
  if (n_refined_cells < m->n_cells) {
    CS_MALLOC(refined_cell_id, m->n_cells, cs_lnum_t);
    n_refined_cells = 0;
    for (cs_lnum_t i = 0; i < m->n_cells; i++) {
      if (c_r_flag[i] > CS_REFINE_NONE)
        refined_cell_id[n_refined_cells++] = i;
      else
        refined_cell_id[n_refined_cells++] = -1;
    }
  }

  t1 = cs_timer_time();

  /* Build additional edge and face indexing and flagging arrays
     ----------------------------------------------------------- */

  /* Build vertices to vertices (edges) graph */

  cs_adjacency_t *v2v = cs_mesh_adjacency_v2v(m);

  cs_lnum_t n_edges = v2v->idx[n_v_ini];

  /* For each edge, e_v_idx will contain the starting index of
     vertices inserted on edges requiring subdivision */

  cs_lnum_t *e_v_idx;
  CS_MALLOC(e_v_idx, n_edges+1, cs_lnum_t);
  for (cs_lnum_t i = 0; i < n_edges+1; i++)
    e_v_idx[i] = 0;

  /* For each face, f_v_idx will contain the starting index of
     vertices inserted on faces requiring subdivision */

  cs_lnum_t *f_v_idx;
  CS_MALLOC(f_v_idx, n_f_ini+1, cs_lnum_t);
  for (cs_lnum_t i = 0; i < n_f_ini+1; i++)
    f_v_idx[i] = 0;

  /* Mark required edges and faces,
     determining face and edge subdivision types */

  char *f_r_level_ini;
  CS_MALLOC(f_r_level_ini, n_f_ini, char);

  cs_gnum_t *g_edges_num = nullptr;
  if (cs_glob_n_ranks > 1)
    CS_MALLOC(g_edges_num, n_edges, cs_gnum_t);

  cs_gnum_t n_g_edges
    = _new_edge_and_face_vertex_ids(m, check_convex, v2v,
                                    m->vtx_r_gen, c_r_level, c_r_flag,
                                    f_r_level_ini, f_r_flag,
                                    e_v_idx, f_v_idx, g_edges_num,
                                    n_add_vtx);

  cs_timer_t t2 = cs_timer_time();
  cs_timer_counter_add_diff(&(timers[2]), &t1, &t2);
  t1 = t2;

  /* Compute number of sub-faces and associated connectivity.
     Size will be transformed to index later so named as index,
     and values for f_id placed in position f_id+1 */

  cs_lnum_t *b_face_o2n_idx, *b_face_o2n_connect_idx;

  CS_MALLOC(b_face_o2n_idx, m->n_b_faces + 1, cs_lnum_t);
  CS_MALLOC(b_face_o2n_connect_idx, m->n_b_faces + 1, cs_lnum_t);

  _subdivided_faces_sizes(v2v,
                          m->vtx_r_gen,
                          f_r_level_ini,
                          e_v_idx,
                          m->n_b_faces,
                          f_r_flag,
                          m->b_face_vtx_idx,
                          m->b_face_vtx_lst,
                          b_face_o2n_idx,
                          b_face_o2n_connect_idx);

  cs_lnum_t *i_face_o2n_idx, *i_face_o2n_connect_idx;

  CS_MALLOC(i_face_o2n_idx, m->n_i_faces + 1, cs_lnum_t);
  CS_MALLOC(i_face_o2n_connect_idx, m->n_i_faces + 1, cs_lnum_t);

  _subdivided_faces_sizes(v2v,
                          m->vtx_r_gen,
                          f_r_level_ini + m->n_b_faces,
                          e_v_idx,
                          m->n_i_faces,
                          f_r_flag + m->n_b_faces,
                          m->i_face_vtx_idx,
                          m->i_face_vtx_lst,
                          i_face_o2n_idx,
                          i_face_o2n_connect_idx);

  t2 = cs_timer_time();
  cs_timer_counter_add_diff(&(timers[3]), &t1, &t2);
  t1 = t2;

  /* Build cell->faces connectivity and identify refined cell types
     ------------------------------------------------------------- */

  cs_adjacency_t *c2f = cs_mesh_adjacency_c2f(m, 0);

  /* Compute element centers (adjusted for subdivision);
     Note that for optimization, this could be done locally
     when adding face and cell centers */

  cs_real_3_t *cell_cen_o, *i_face_cen_o, *b_face_cen_o;

  _element_centers(m, f_r_flag, &cell_cen_o, &i_face_cen_o, &b_face_cen_o);

  /* Now determine cell subdivision types
     ------------------------------------ */

  cs_lnum_t *c2f2v_start;
  CS_MALLOC(c2f2v_start, c2f->idx[c2f->n_elts], cs_lnum_t);

  _cell_r_types(m, conforming, c2f, c2f2v_start, c_r_flag, f_r_flag);

  /* For each cell, c_v_idx will contain the starting index of
     vertices inserted inside cells requiring subdivision */

  cs_lnum_t *c_v_idx;
  CS_MALLOC(c_v_idx, n_c_ini+1, cs_lnum_t);
  for (cs_lnum_t i = 0; i < n_c_ini+1; i++)
    c_v_idx[i] = 0;

  /* Mark required cell interior vertex additions */

  _new_cell_vertex_ids(m, c_r_flag, c_v_idx, n_add_vtx);

  t2 = cs_timer_time();
  cs_timer_counter_add_diff(&(timers[4]), &t1, &t2);
  t1 = t2;

  /* Free mesh data that will be rebuilt, as it would become
     inconsistent once the mesh is modified. */

  cs_mesh_free_rebuildable(m, true);

  /* Now build new vertices
     ---------------------- */

  cs_lnum_t n_vtx_new =   m->n_vertices
                        + n_add_vtx[0] + n_add_vtx[1] + n_add_vtx[2];

  cs_lnum_t n_g_vtx_new =  m->n_g_vertices;

  CS_REALLOC(m->vtx_coord, n_vtx_new*3, cs_real_t);
  CS_REALLOC(m->vtx_r_gen, n_vtx_new, char);
  if (m->global_vtx_num != nullptr)
    CS_REALLOC(m->global_vtx_num, n_vtx_new, cs_gnum_t);

  _build_edge_vertices(m, v2v, n_add_vtx[0], e_v_idx, g_edges_num);
  _build_vertices_gnum(m, n_edges, n_g_edges, e_v_idx, g_edges_num);

  CS_FREE(g_edges_num);

  _build_face_vertices(m,
                       m->n_b_faces,
                       m->b_face_vtx_idx,
                       m->b_face_vtx_lst,
                       f_r_flag,
                       f_v_idx,
                       (const cs_real_3_t *)b_face_cen_o);
  _build_face_vertices(m,
                       m->n_i_faces,
                       m->i_face_vtx_idx,
                       m->i_face_vtx_lst,
                       f_r_flag + m->n_b_faces,
                       f_v_idx + m->n_b_faces,
                       (const cs_real_3_t *)i_face_cen_o);

  CS_FREE(b_face_cen_o);
  CS_FREE(i_face_cen_o);

  _build_vertices_gnum(m, m->n_b_faces, m->n_g_b_faces,
                       f_v_idx, m->global_b_face_num);
  _build_vertices_gnum(m, m->n_i_faces, m->n_g_i_faces,
                       f_v_idx + m->n_b_faces, m->global_i_face_num);

  _build_cell_vertices(m,
                       c2f,
                       n_b_f_ini,
                       m->n_cells,
                       c_r_flag,
                       c_r_level,
                       c_v_idx,
                       (const cs_real_3_t *)cell_cen_o);

  CS_FREE(cell_cen_o);

  _build_vertices_gnum(m, m->n_cells, m->n_g_cells,
                       c_v_idx, m->global_cell_num);

  /* Update counts */

  m->n_vertices = n_vtx_new;
  m->n_g_vertices = n_g_vtx_new;

  t2 = cs_timer_time();
  cs_timer_counter_add_diff(&(timers[5]), &t1, &t2);
  t1 = t2;

  /* Update faces
     ------------ */

  _update_face_connectivity(v2v,
                            m->vtx_coord,
                            m->vtx_r_gen,
                            e_v_idx,
                            f_v_idx,
                            m->n_b_faces,
                            f_r_flag,
                            f_r_level_ini,
                            b_face_o2n_idx,
                            b_face_o2n_connect_idx,
                            &(m->b_face_vtx_idx),
                            &(m->b_face_vtx_lst));

  CS_FREE(b_face_o2n_connect_idx);

  _update_face_connectivity(v2v,
                            m->vtx_coord,
                            m->vtx_r_gen,
                            e_v_idx,
                            f_v_idx + m->n_b_faces,
                            m->n_i_faces,
                            f_r_flag + m->n_b_faces,
                            f_r_level_ini + m->n_b_faces,
                            i_face_o2n_idx,
                            i_face_o2n_connect_idx,
                            &(m->i_face_vtx_idx),
                            &(m->i_face_vtx_lst));

  CS_FREE(i_face_o2n_connect_idx);

  /* Update arrays and counts based on faces (families) and number of faces */

  _o2n_idx_update_b_face_arrays(m, b_face_o2n_idx);
  _o2n_idx_update_i_face_arrays(m, i_face_o2n_idx);

  t2 = cs_timer_time();
  cs_timer_counter_add_diff(&(timers[6]), &t1, &t2);
  t1 = t2;

  /* Now subdivide cells */

  /* Count number of sub-cells, added interior faces and their connectivity
     size due to refinement (note face_o2_counts arrays are shifted by 1
     because they were built to be transformed as indexes, with
     initial values shifted by 1). */

  cs_lnum_t *c_poly_s_idx;
  cs_lnum_t *c_o2n_idx, *c_i_face_idx, *c_i_faces_size;

  CS_MALLOC(c_poly_s_idx, n_c_ini+1, cs_lnum_t);
  CS_MALLOC(c_o2n_idx, n_c_ini + 1, cs_lnum_t);
  CS_MALLOC(c_i_face_idx, n_c_ini + 1, cs_lnum_t);
  CS_MALLOC(c_i_faces_size, n_c_ini, cs_lnum_t);

  c_poly_s_idx[0] = 0;
  c_o2n_idx[0] = 0;
  c_i_face_idx[0] = m->n_i_faces;

  _new_cells_i_faces_count(m,
                           c2f,
                           n_b_f_ini,
                           c_r_flag,
                           b_face_o2n_idx,
                           i_face_o2n_idx,
                           c_poly_s_idx + 1,
                           c_o2n_idx + 1,
                           c_i_face_idx + 1,
                           c_i_faces_size);

  _counts_to_index(n_c_ini, c_poly_s_idx);

  /* If more complex subdivision is needed, allocate required
     structures and precompute division. */

  cs_lnum_t *c_o_id = nullptr;
  cs_lnum_t *e_f2c_pre_idx = nullptr, *i_f2c_pre_idx = nullptr;
  int16_t *e_f2c_pre = nullptr, *i_f2c_pre = nullptr;
  cs_lnum_t *i_f2v_pre_idx = nullptr, *i_f2v_pre = nullptr;

  cs_lnum_t n_poly_r_cells = c_poly_s_idx[n_c_ini];
  if (n_poly_r_cells > 0) {
    _new_cells_poly_pre_alloc(n_c_ini,
                              c_poly_s_idx,
                              c_o2n_idx + 1,    // also number of exterior faces
                              c_i_face_idx + 1, // still a count at this stage
                              c_i_faces_size,
                              c_o_id,
                              e_f2c_pre_idx,
                              e_f2c_pre,
                              i_f2c_pre_idx,
                              i_f2c_pre,
                              i_f2v_pre_idx,
                              i_f2v_pre);

    int n_threads_r = cs_parall_n_threads(n_poly_r_cells, 1);

    cs_cell_refine_helper *r_helpers;
    CS_MALLOC(r_helpers, n_threads_r, cs_cell_refine_helper);
    for (int i = 0; i < n_threads_r; i++)
      new(r_helpers+i) cs_cell_refine_helper();

    using c2f_pre_t = int16_t[2];

#   pragma omp parallel for schedule(dynamic) num_threads(n_threads_r)
    for (cs_lnum_t c_idx = 0; c_idx < n_poly_r_cells; c_idx++) {
      cs_cell_refine_helper *crh = r_helpers + cs_get_thread_id();

      cs_lnum_t c_id = c_o_id[c_idx];

      _pre_refine_poly_cell(m,
                            c2f,
                            crh,
                            c_id,
                            n_b_f_ini,
                            c_r_flag[c_id],
                            b_face_o2n_idx,
                            i_face_o2n_idx,
                            c_v_idx,
                            f_r_level_ini,
                            c_o2n_idx[c_id+1],     // c_n_s_cells
                            c_i_face_idx[c_id+1],  // c_n_i_faces
                            c_i_faces_size[c_id],
                            e_f2c_pre + e_f2c_pre_idx[c_idx],
                            ((c2f_pre_t *)i_f2c_pre) + i_f2c_pre_idx[c_idx],
                            i_f2v_pre + i_f2v_pre_idx[c_idx]);
    }

    _free_cell_refine_helpers(n_threads_r, r_helpers);
  }

  _counts_to_index(n_c_ini, c_o2n_idx);
  _counts_to_index(n_c_ini, c_i_face_idx);

  _o2n_idx_update_cell_arrays(m, n_c_ini, c_o2n_idx);

  t2 = cs_timer_time();
  cs_timer_counter_add_diff(&(timers[7]), &t1, &t2);
  t1 = t2;

  _o2n_idx_update_cell_i_face_arrays(m, n_c_ini, c_i_face_idx, c_i_faces_size);
  CS_FREE(c_i_faces_size);

  _subdivide_cells(m,
                   n_c_ini,
                   n_b_f_ini,
                   c_o2n_idx,
                   i_face_o2n_idx,
                   b_face_o2n_idx,
                   c2f,
                   c2f2v_start,
                   c_v_idx,
                   c_i_face_idx,
                   c_r_flag,
                   c_r_level,
                   c_poly_s_idx,
                   e_f2c_pre_idx, e_f2c_pre,
                   i_f2c_pre_idx, i_f2c_pre,
                   i_f2v_pre_idx, i_f2v_pre);

  t2 = cs_timer_time();
  cs_timer_counter_add_diff(&(timers[8]), &t1, &t2);
  t1 = t2;

  /* Cleanup*/

  CS_FREE(i_f2v_pre);
  CS_FREE(i_f2v_pre_idx);
  CS_FREE(i_f2c_pre);
  CS_FREE(i_f2c_pre_idx);
  CS_FREE(e_f2c_pre);
  CS_FREE(e_f2c_pre_idx);

  CS_FREE(c_o_id);
  CS_FREE(c_poly_s_idx);

  CS_FREE(f_r_level_ini);

  CS_FREE(c2f2v_start);

  CS_FREE(c_o2n_idx);
  CS_FREE(c_i_face_idx);

  CS_FREE(refined_cell_id);
  cs_adjacency_destroy(&c2f);

  CS_FREE(c_v_idx);
  CS_FREE(f_v_idx);
  CS_FREE(e_v_idx);
  cs_adjacency_destroy(&v2v);

  CS_FREE(i_face_o2n_idx);
  CS_FREE(b_face_o2n_idx);

  CS_FREE(c_r_level);
  CS_FREE(c_r_flag);
  CS_FREE(f_r_flag);

  m->have_r_gen = true;
  m->modified |= (CS_MESH_MODIFIED | CS_MESH_MODIFIED_BALANCE);

  /* Rebuild ghosts */

  mv_save = m->verbosity;
  m->verbosity = -1;

  if (   m->n_domains > 1 || m->n_init_perio > 0
      || m->halo_type == CS_HALO_EXTENDED) {
    cs_halo_type_t halo_type = m->halo_type;
    cs_mesh_builder_t *mb = (m == cs_glob_mesh) ? cs_glob_mesh_builder : nullptr;
    cs_mesh_init_halo(m, mb, halo_type, -1, true);
  }

  cs_mesh_update_auxiliary(m);

  m->verbosity = mv_save;

  t2 = cs_timer_time();
  cs_timer_counter_add_diff(&(timers[1]), &t1, &t2);
  t1 = t2;

  cs_timer_counter_add_diff(&(timers[0]), &t0, &t2);

  if (m->verbosity > 0) {

    _print_mesh_counts(m, _("Mesh after refinement"));
    cs_log_printf(CS_LOG_DEFAULT, "\n");
    cs_log_separator(CS_LOG_DEFAULT);

    cs_log_printf
      (CS_LOG_PERFORMANCE,
       _("\nMesh refinement:\n\n"
         "  Pre and post update of mesh structure:        %.3g\n"
         "  Edge and face indexing:                       %.3g\n"
         "  Build sub-face indexing:                      %.3g\n"
         "  Cell-faces indexing and type identification:  %.3g\n"
         "  Build new vertices:                           %.3g\n"
         "  Subdivide faces:                              %.3g\n"
         "  Build added interior faces indexing:          %.3g\n"
         "  Subdivide cells                               %.3g\n\n"
         "  Total:                                        %.3g\n"),
       (double)(timers[1].nsec*1.e-9),
       (double)(timers[2].nsec*1.e-9),
       (double)(timers[3].nsec*1.e-9),
       (double)(timers[4].nsec*1.e-9),
       (double)(timers[5].nsec*1.e-9),
       (double)(timers[6].nsec*1.e-9),
       (double)(timers[7].nsec*1.e-9),
       (double)(timers[8].nsec*1.e-9),
       (double)(timers[0].nsec*1.e-9));
    cs_log_printf(CS_LOG_PERFORMANCE, "\n");
    cs_log_separator(CS_LOG_PERFORMANCE);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Refine selected mesh cells.
 *
 * \param[in, out]  m           mesh
 * \param[in]       conforming  if true, propagate refinement to ensure
 *                              subdivision is conforming
 * \param[in]       n_cells     number of selected cells
 * \param[in]       cells       list of selected cells (0 to n-1)
 *                              or nullptr if no indirection is needed
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_refine_simple_selected(cs_mesh_t        *m,
                               bool              conforming,
                               cs_lnum_t         n_cells,
                               const cs_lnum_t   cells[])
{
  cs_lnum_t n_c_ini = m->n_cells;

  int *cell_flag;
  CS_MALLOC(cell_flag, n_c_ini, int);
  for (cs_lnum_t i = 0; i < n_c_ini; i++)
    cell_flag[i] = 0;

  if (cells != nullptr) {
    for (cs_lnum_t i = 0; i < n_cells; i++)
      cell_flag[cells[i]] = 1;
  }
  else {
    for (cs_lnum_t i = 0; i < n_cells; i++)
      cell_flag[i] = 1;
  }

  cs_mesh_refine_simple(m, conforming, cell_flag);

  CS_FREE(cell_flag);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set refinement options, using key/value pairs.
 *
 * Accepted keys and values:
 *
 * - "triangle_subdivision"
 *   - "triangle" (default)
 *   - "quadrangle"
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_refine_set_option(const char  *key,
                          const char  *value)
{
  bool wrong_value = false;

  if (key == nullptr || value == nullptr)
    return;

  if (strcmp(key, "triangle_subdivision") == 0) {
    if (strcmp(value, "triangle") == 0)
      _refine_tria_type = CS_REFINE_TRIA;
    else if (strcmp(value, "quadrangle") == 0)
      _refine_tria_type = CS_REFINE_TRIA_Q;
    else
      wrong_value = true;
  }

  else
    bft_error(__FILE__, __LINE__, 0,
              _("%s: unknown key \"%s\"."), __func__, key);

  if (wrong_value)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: unknown value \"%s\" for key \"%s\"."),
              __func__, value, key);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
