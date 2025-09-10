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
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

#include <cmath>
#include <chrono>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_printf.h"
#include "base/cs_math.h"

#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_adjacencies.h"
#include "mesh/cs_mesh_boundary.h"
#include "mesh/cs_mesh_group.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "mesh/cs_mesh_cut.h"

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

constexpr cs_real_t _plane_tol = 1e-12;
constexpr int       _default_family_id = 1;

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*------------------------------------------------------------------------------
 * This structure holds all the necessary buffers for the algorithm.
 * The buffers are allocated once, with enough memory to successfully handle
 * all the cells to be cut.
 * Unless specified, all the buffers are reset and filled on a cell-by-cell
 * basis.
 *
 * Members:
 *   sd       : signed distances of cell vertices to the cut plane.
 *   occurs   : number of occurrances of vertices in a polyline.
 *   e_stride : maximum number of edges per face.
 *   f_size   : maximum size of the face connectivity array induced by the
 *              cut.
 *   edges    : edge-vertex connectivity.
 *   e2f      : edge-face connectivity.
 *   xyz      : coordinates of the intersection points.
 *   indices  : indirection table, useful for eliminating duplicate vertices.
 *   compact  : maps unique vertices to their indices within the mesh.
 *   n_polys  : number of closing polygons created by the cutting algorithm.
 *              This variable gets incremented after every cell cut.
 *   polys    : the indices of the closing polygons are appended to this array.
 *----------------------------------------------------------------------------*/

struct _cut_data {
  cs_real_t    *sd;
  uint8_t      *occurs;
  cs_lnum_t     e_stride;
  cs_lnum_t     f_size;
  cs_lnum_2_t  *edges;
  cs_lnum_t    *e2f;
  cs_lnum_t    *f2e;
  cs_real_t    *xyz;
  cs_lnum_t    *indices;
  cs_lnum_t    *compact;
  cs_lnum_t     n_polys;
  cs_lnum_t    *polys;
};

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*------------------------------------------------------------------------------
 * Reverse the order of the elements in the range [arr, arr+n]
 *----------------------------------------------------------------------------*/

static void
_reverse_array(cs_lnum_t *arr, cs_lnum_t n)
{
  cs_lnum_t start = 0, end = n-1, temp;
  while (start < end) {
    temp = arr[start];
    arr[start] = arr[end];
    arr[end] = temp;
    start++;
    end--;
  }
}

/*------------------------------------------------------------------------------
 * A |      | D
 *   |      |
 * B |______| C
 *
 * If a polyline is open, find the start and end vertices (e.g. A and D), and
 * return true.
 * Else, return false.
 *----------------------------------------------------------------------------*/

static bool
_get_open_polyline(cs_lnum_t    fe[],
                   cs_lnum_t    e_stride,
                   uint8_t     *occurs,
                   cs_lnum_t   *start,
                   cs_lnum_t   *end,
                   cs_lnum_2_t  edges[])
{
  for (cs_lnum_t i = 0; i < e_stride; i++) {
    if (fe[i] == -1) continue;
    cs_lnum_t e_id = fe[i];
    occurs[edges[e_id][0]]++;
    occurs[edges[e_id][1]]++;
  }

  *start = -1;
  *end = -1;

  for (cs_lnum_t i = 0; i < e_stride; i++) {
    if (fe[i] == -1) continue;
    cs_lnum_t e_id = fe[i];
    cs_lnum_t p = edges[e_id][0], q = edges[e_id][1];

    if (occurs[p] == 1) {
      if (*start == -1) *start = p; else if (*end == -1) *end = p;
    }
    if (occurs[q] == 1) {
      if (*start == -1) *start = q; else if (*end == -1) *end = q;
    }
  }

  return *start != -1;
}

static inline void
_get_b_face_vertices(cs_mesh_t *mesh,
                     cs_lnum_t  f_id,
                     cs_lnum_t *fv[],
                     cs_lnum_t *nv)
{
  *nv = mesh->b_face_vtx_idx[f_id+1] - mesh->b_face_vtx_idx[f_id];
  *fv = mesh->b_face_vtx_lst + mesh->b_face_vtx_idx[f_id];
}

/*------------------------------------------------------------------------------
 * Remove edge e_id from the face-edge connectivity of face f_id.
 *----------------------------------------------------------------------------*/

static void
_remove_edge_from_face(cs_lnum_t f2e[],
                       cs_lnum_t e_stride,
                       cs_lnum_t f_id,
                       cs_lnum_t e_id,
                       cs_lnum_t e2f[])
{
  cs_lnum_t *fe = f2e + e_stride*f_id;
  cs_lnum_t i;
  for (i = 0; i < e_stride; i++) {
    if (fe[i] == e_id) {
      fe[i] = -1;
      e2f[e_id] = -1;
      break;
    }
  }
  assert(i != e_stride);
}

/*------------------------------------------------------------------------------
 * Add edge e_id to the face-edge connectivity of face f_id.
 *----------------------------------------------------------------------------*/

static void
_add_edge_to_face(cs_lnum_t f2e[],
                  cs_lnum_t e_stride,
                  cs_lnum_t f_id,
                  cs_lnum_t e_id,
                  cs_lnum_t e2f[])
{
  cs_lnum_t *fe = f2e + e_stride*f_id;
  cs_lnum_t i;
  for (i = 0; i < e_stride; i++) {
    if (fe[i] == -1) {
      fe[i] = e_id;
      e2f[e_id] = f_id;
      break;
    }
  }
  assert(i != e_stride);
}

/*------------------------------------------------------------------------------
 * Lexicographical comparison of two vertex coordinates.
 *----------------------------------------------------------------------------*/

static int
_cmp_crd(const cs_real_t x[3], const cs_real_t y[3])
{
  for (int i = 0; i < 3; i++) {
    if (x[i] < y[i]) return -1;
    if (x[i] > y[i]) return 1;
  }
  return 0;
}

/*-----------------------------------------------------------------------------
 * Deduce the face-vertex connectivity from the face-edge connectivity.
 *----------------------------------------------------------------------------*/

static void
_fill_connectivity(cs_lnum_t   b_face_vtx_lst[],
                   cs_lnum_t   e_stride,
                   cs_lnum_t   f_id,
                   cs_lnum_t   fe[],
                   cs_lnum_2_t edges[],
                   cs_lnum_t   b_face_vtx_idx[],
                   cs_lnum_t   ne)
{
  cs_lnum_t *fv = b_face_vtx_lst + e_stride*f_id;
  for (cs_lnum_t j = 0; j < e_stride; j++) fv[j] = -1;

  cs_lnum_t count;
  for (count = 0; count < 2; count++) fv[count] = edges[fe[0]][count];

  while (count != ne) {
    for (cs_lnum_t j = 1; j < e_stride; j++) {
      cs_lnum_t e_id = fe[j];
      if (fe[j] == -1) continue;

      for (int k = 0; k < 2; k++) {
        if (edges[e_id][k] == fv[count-1]) {
          fv[count++] = edges[e_id][(k+1)%2];
          fe[j] = -1;
          break;
        }
      }
      if (fe[j] == -1) break;
    }
  }

  b_face_vtx_idx[f_id+1] = ne;
}

/*------------------------------------------------------------------------------
 * Compute the normal of a face given its vertex connectivity.
 * For now, approximates it with the normal of the triangle formed by the
 * first three vertices.
 * TODO: more robust computation.
 *----------------------------------------------------------------------------*/

static void
_get_normal(const cs_lnum_t  *fv,
            const cs_real_t  *crd,
            cs_real_t         normal[3])
{
  /* Get the normal to the first triangle. */
  cs_lnum_t a = fv[0], b = fv[1], c = fv[2];
  cs_real_t u[3], v[3];
  for (int i = 0; i < 3; i++) {
    u[i] = crd[3*b+i] - crd[3*a+i];
    v[i] = crd[3*c+i] - crd[3*a+i];
  }
  cs_math_3_cross_product(u, v, normal);
}

/*-----------------------------------------------------------------------------
 * After getting the face-vertex connectivity of a new child face,
 * recover its correct orientation using the orientation of its parent face.
 *----------------------------------------------------------------------------*/

static void
_fill_connectivity_and_reorient(cs_mesh_t   *mesh,
                                cs_lnum_t    b_face_vtx_lst[],
                                cs_lnum_t    e_stride,
                                cs_lnum_t    f_id,
                                cs_lnum_t    fe[],
                                cs_lnum_2_t  edges[],
                                cs_lnum_t    b_face_vtx_idx[],
                                cs_lnum_t    ne,
                                cs_lnum_t    o_fid)
{
  _fill_connectivity(b_face_vtx_lst,
                     e_stride,
                     f_id,
                     fe,
                     edges,
                     b_face_vtx_idx,
                     ne);

  /* Get the original normal. */
  cs_real_t o_normal[3] = {0, 0, 0};
  cs_lnum_t *fv = mesh->b_face_vtx_lst + mesh->b_face_vtx_idx[o_fid];
  _get_normal(fv, mesh->vtx_coord, o_normal);

  /* Compute current normal. */
  cs_real_t c_normal[3] = {0.0, 0.0, 0.0};
  fv = b_face_vtx_lst + f_id*e_stride;
  _get_normal(fv, mesh->vtx_coord, c_normal);

  /* Swap the orientation of the current face if necessary. */
  cs_real_t dp = cs_math_3_dot_product(o_normal, c_normal);
  if (dp >= 0.0) return;
  _reverse_array(fv+1, ne-1);
}

/*-----------------------------------------------------------------------------
 * Allocate the buffers necessary for the algorithm with enough memory to be
 * reused across all the cells.
 *----------------------------------------------------------------------------*/

static _cut_data
_allocate_scratch_data(cs_mesh_t *mesh,
                       cs_lnum_t  max_nf,
                       cs_lnum_t  max_nv,
                       cs_lnum_t  n_new_vertices,
                       cs_lnum_t  n_new_cells)
{
  _cut_data cd;
  CS_MALLOC(cd.sd, mesh->n_vertices + n_new_vertices, cs_real_t);
  CS_MALLOC(cd.occurs, mesh->n_vertices + n_new_vertices, uint8_t);
  cd.e_stride = CS_MAX(max_nv+1, max_nf);
  cd.f_size = (2*max_nf+2)*cd.e_stride;
  CS_MALLOC(cd.edges, cd.f_size, cs_lnum_2_t);
  CS_MALLOC(cd.e2f, cd.f_size, cs_lnum_t);
  CS_MALLOC(cd.f2e, cd.f_size, cs_lnum_t);
  CS_MALLOC(cd.xyz, 3*n_new_vertices, cs_real_t);
  CS_MALLOC(cd.indices, n_new_vertices, cs_lnum_t);
  CS_MALLOC(cd.compact, n_new_vertices, cs_lnum_t);
  cd.n_polys = 0;
  CS_MALLOC(cd.polys, 2*n_new_cells, cs_lnum_t);
  return cd;
}

/*-----------------------------------------------------------------------------
 * Cut the cell c_id with the plane (p_normal, p_origin).
 * A cut cell is replaced by the cell lying in the negative half-space
 * spanned by the plane. The cell lying in the positive half-space is appended
 * to the mesh connectivity arrays.
 * Set the new face-vertex connectivity data inside b_face_vtx_idx and
 * b_face_vtx_ids.
 *----------------------------------------------------------------------------*/

static void
_cut_cell(cs_mesh_t            *mesh,
          _cut_data            *cd,
          cs_lnum_t             c_id,
          const cs_real_t       p_normal[],
          const cs_real_t       p_origin[],
          const cs_adjacency_t *c2v,
          const cs_adjacency_t *c2f_b,
          cs_lnum_t            *b_face_vtx_idx,
          cs_lnum_t            *b_face_vtx_lst)
{
  cs_real_t *sd = cd->sd;
  uint8_t *occurs = cd->occurs;
  cs_lnum_t e_stride = cd->e_stride;
  cs_lnum_t f_size = cd->f_size;
  cs_lnum_2_t *edges = cd->edges;
  cs_lnum_t *e2f = cd->e2f;
  cs_lnum_t *f2e = cd->f2e;
  cs_real_t *xyz = cd->xyz;
  cs_lnum_t *indices = cd->indices;
  cs_lnum_t *compact = cd->compact;
  cs_lnum_t *polys = cd->polys;

  /* Reset. */
  for (cs_lnum_t i = 0; i < f_size; i++) {
    edges[i][0] = edges[i][1] = -1;
    e2f[i] = -1;
    f2e[i] = -1;
  }

  /* Compute the signed distances of the cell vertices to the cut plane. */
  int positive = 0, negative = 0;

  for (cs_lnum_t i = c2v->idx[c_id]; i < c2v->idx[c_id+1]; i++) {
    cs_lnum_t p = c2v->ids[i];
    const cs_real_t *crd = mesh->vtx_coord + 3*p;
    sd[p] = 0.0;
    for (int j = 0; j < 3; j++) sd[p] += (crd[j] - p_origin[j]) * p_normal[j];
    if (sd[p] >= _plane_tol) positive++;
    else if (sd[p] <= -_plane_tol) negative++;
    else sd[p] = 0.0;
  }

  /* If the plane does not cut the cell, do nothing. */
  if (positive == 0 || negative == 0) return;

  /* Make the edge-face connectivity. */
  cs_lnum_t n_edges = 0;
  cs_lnum_t nf = c2f_b->idx[c_id+1] - c2f_b->idx[c_id];
  const cs_lnum_t *cf = c2f_b->ids + c2f_b->idx[c_id];

  for (cs_lnum_t i = 0; i < nf; i++) {
    cs_lnum_t f_id = cf[i];
    cs_lnum_t *fv, nv;
    _get_b_face_vertices(mesh, f_id, &fv, &nv);

    cs_lnum_t *fe = f2e + e_stride*i;

    for (cs_lnum_t j = 0; j < nv; j++) {
      cs_lnum_t p = fv[j];
      cs_lnum_t q = fv[(j+1)%nv];
      edges[n_edges][0] = p;
      edges[n_edges][1] = q;
      e2f[n_edges] = i;
      fe[j] = n_edges;
      n_edges++;
    }
  }

  /* Process the edges. */
  cs_lnum_t n_cut_e = 0;

  for (cs_lnum_t e_id = 0; e_id < n_edges; e_id++) {
    cs_lnum_t p = edges[e_id][0];
    cs_lnum_t q = edges[e_id][1];

    cs_real_t dp = sd[p];
    cs_real_t dq = sd[q];

    cs_lnum_t f = e2f[e_id];

    if (dp <= 0 && dq <= 0) {
      /* Edge on the negative side. Do nothing. */
      continue;
    }

    if (dp >= 0 && dq >= 0) {
      /* Edge is on the positive side. */
      _remove_edge_from_face(f2e, e_stride, f, e_id, e2f);
      _add_edge_to_face(f2e, e_stride, f + nf, e_id, e2f);
      continue;
    }

    /* Split the edge by the plane. */
    cs_lnum_t a = p, b = q;
    if (a > b) { cs_lnum_t t = a; a = b; b = t; }
    cs_real_t da = sd[a];
    cs_real_t db = sd[b];
    cs_real_t t = da / (da - db);
    cs_real_t *crd = xyz + 3*n_cut_e;
    const cs_real_t *a_crd = mesh->vtx_coord + 3*a;
    const cs_real_t *b_crd = mesh->vtx_coord + 3*b;
    for (int k = 0; k < 3; k++) crd[k] = (1.0-t)*a_crd[k] + t*b_crd[k];

    cs_lnum_t new_e_id = n_edges + n_cut_e;

    if (dp > 0) {
      edges[e_id][0] = q;
      edges[e_id][1] = mesh->n_vertices + n_cut_e;
      edges[new_e_id][0] = mesh->n_vertices + n_cut_e;
      edges[new_e_id][1] = p;
    }
    else {
      edges[e_id][0] = p;
      edges[e_id][1] = mesh->n_vertices + n_cut_e;
      edges[new_e_id][0] = mesh->n_vertices + n_cut_e;
      edges[new_e_id][1] = q;
    }

    _add_edge_to_face(f2e, e_stride, f + nf, new_e_id, e2f);

    n_cut_e++;
  }

  n_edges += n_cut_e;

  /* Get rid of duplicate vertices */
  for (cs_lnum_t i = 0; i < n_cut_e; i++) indices[i] = i;

  cs_lnum_t n_dups = 0;
  for (cs_lnum_t i = 0; i < n_cut_e; i++) {
    cs_real_t *i_crd = xyz + 3*i;
    for (cs_lnum_t j = i+1; j < n_cut_e; j++) {
      cs_real_t *j_crd = xyz + 3*j;
      if (_cmp_crd(i_crd, j_crd) == 0) {
        indices[j] = indices[i];
        n_dups++;
      }
    }
  }

  cs_lnum_t unique_vtx = 0;
  for (int i = 0; i < n_cut_e; i++) compact[i] = -1;

  for (cs_lnum_t i = 0; i < n_cut_e; i++) {
    if (indices[i] == i) {
      compact[i] = unique_vtx++;

      /* Insert the coordinates into the mesh. */
      cs_real_t *crd = mesh->vtx_coord + 3*(mesh->n_vertices+compact[i]);
      const cs_real_t *_crd = xyz + 3*i;
      for (int k = 0; k < 3; k++) crd[k] = _crd[k];
    }
  }

  assert(unique_vtx == n_cut_e-n_dups);

  for (cs_lnum_t i = 0; i < n_edges; i++) {
    for (int j = 0; j < 2; j++) {
      cs_lnum_t p = edges[i][j];
      if (p >= mesh->n_vertices) {
        cs_lnum_t delta = p - mesh->n_vertices;
        cs_lnum_t first = indices[delta];
        cs_lnum_t index = compact[first];
        assert(index != -1);
        edges[i][j] = mesh->n_vertices + index;
      }
    }
  }

  /* Process the faces. */

  cs_lnum_t n_close_e = 0;
  cs_lnum_t seeds[2] = {-1, -1};

  for (cs_lnum_t i = 0; i < nf; i++) {
    cs_lnum_t *fe = f2e + i*e_stride;

    cs_lnum_t old_f = i;
    cs_lnum_t new_f = i + nf;
    cs_lnum_t closing_f = 2*nf;

    for (cs_lnum_t j = 0; j < e_stride; j++) {
      if (fe[j] == -1) continue;
      cs_lnum_t e_id = fe[j];
      occurs[edges[e_id][0]] = 0;
      occurs[edges[e_id][1]] = 0;
    }

    cs_lnum_t start, end;
    if (_get_open_polyline(fe, e_stride, occurs, &start, &end, edges)) {
      /* Polyline is open, close it. */
      cs_lnum_t e_id = n_edges + n_close_e;
      edges[e_id][0] = start;
      edges[e_id][1] = end;

      /* Add the edge to the old and new faces and to the two closing faces. */
      _add_edge_to_face(f2e, e_stride, old_f, e_id, e2f);
      _add_edge_to_face(f2e, e_stride, new_f, e_id, e2f);
      _add_edge_to_face(f2e, e_stride, closing_f, e_id, e2f);
      _add_edge_to_face(f2e, e_stride, closing_f+1, e_id, e2f);

      n_close_e++;

      if (seeds[0] == -1) {
        assert(seeds[1] == -1);
        /* Link the closing polygons with the first faces that share one of its
         * edges. Useful to set their correct orientations later.
         */
        seeds[0] = old_f;
        seeds[1] = new_f;
      }
    }
  }

  assert(seeds[0] != -1);
  assert(seeds[1] != -1);

  /* Fill in the new face connectivities */

  cs_lnum_t f_incr = 0;
  bool seeds_updated = false;

  for (cs_lnum_t i = 0; i < nf; i++) {
    /* Old face. */
    cs_lnum_t *fe = f2e + i*e_stride;

    cs_lnum_t ne = 0;
    for (cs_lnum_t j = 0; j < e_stride; j++) {
      if (fe[j] != -1) ne++;
    }

    cs_lnum_t f_id = cf[i];

    if (ne == 0) {
      /* Old face belongs to the positive side.
       * Copy its connectivity but do not increment n_b_faces.
       */
      fe = f2e + (i+nf)*e_stride;
      ne = 0;
      for (cs_lnum_t j = 0; j < e_stride; j++) {
        if (fe[j] != -1) ne++;
      }
      assert(ne != 0);

      _fill_connectivity_and_reorient(mesh,
                                      b_face_vtx_lst,
                                      e_stride,
                                      f_id,
                                      fe,
                                      edges,
                                      b_face_vtx_idx,
                                      ne,
                                      f_id);

      mesh->b_face_cells[f_id] = mesh->n_cells;

      continue;
    }

    assert(ne <= e_stride);

    _fill_connectivity_and_reorient(mesh,
                                    b_face_vtx_lst,
                                    e_stride,
                                    f_id,
                                    fe,
                                    edges,
                                    b_face_vtx_idx,
                                    ne,
                                    f_id);

    /* New face. */
    fe = f2e + (i+nf)*e_stride;

    ne = 0;
    for (cs_lnum_t j = 0; j < e_stride; j++) {
      if (fe[j] != -1) ne++;
    }

    if (ne == 0) continue;

    cs_lnum_t f_id_new = mesh->n_b_faces + f_incr++;

    if (!seeds_updated) {
      if (i == seeds[0]) {
        seeds_updated = true;
        seeds[0] = f_id;
        seeds[1] = f_id_new;
      }
    }

    _fill_connectivity_and_reorient(mesh,
                                    b_face_vtx_lst,
                                    e_stride,
                                    f_id_new,
                                    fe,
                                    edges,
                                    b_face_vtx_idx,
                                    ne,
                                    f_id);

    assert(mesh->b_face_family[f_id_new] == _default_family_id);
    assert(mesh->b_face_cells[f_id_new] == -1);

    mesh->b_face_cells[f_id_new] = mesh->n_cells;
    mesh->b_face_family[f_id_new] = mesh->b_face_family[f_id];
  }

  assert(seeds_updated);

  for (cs_lnum_t i = 0; i < 2; i++) {
    cs_lnum_t *fe = f2e + (i+2*nf)*e_stride;

    cs_lnum_t ne = 0;
    for (cs_lnum_t j = 0; j < e_stride; j++) {
      if (fe[j] != -1) ne++;
    }

    assert(ne != 0);
    assert(ne <= e_stride);

    cs_lnum_t f_id = mesh->n_b_faces + f_incr++;

    _fill_connectivity(b_face_vtx_lst,
                       e_stride,
                       f_id,
                       fe,
                       edges,
                       b_face_vtx_idx,
                       ne);

    if (i == 0) {
      mesh->b_face_cells[f_id] = c_id;
    } else {
      mesh->b_face_cells[f_id] = mesh->n_cells;
    }

    mesh->b_face_family[f_id] = _default_family_id;

    polys[cd->n_polys++] = f_id;

    /* Set its correct orientation. */
    bool reverse = false;
    cs_lnum_t *fv = b_face_vtx_lst + f_id * e_stride;

    const cs_lnum_t *_fv = b_face_vtx_lst + seeds[i] * e_stride;
    cs_lnum_t _ne = b_face_vtx_idx[seeds[i]+1];

    bool found = false;
    for (cs_lnum_t j = 0; j < _ne; j++) {
      cs_lnum_t _p = _fv[j];
      cs_lnum_t _q = _fv[(j+1)%_ne];
      for (cs_lnum_t k = 0; k < ne && !found; k++) {
        cs_lnum_t p = fv[k];
        cs_lnum_t q = fv[(k+1)%ne];

        if (p == _p && q == _q) {
          found = true;
          reverse = true;
          break;
        } else if (p == _q && q == _p) {
          found = true;
          reverse = false;
          break;
        }
      }
    }
    assert(found);

    if (reverse) {
      _reverse_array(fv+1, ne-1);
    }
  }

  /* Update cell family. */
  assert(mesh->cell_family[mesh->n_cells] == _default_family_id);
  mesh->cell_family[mesh->n_cells] = mesh->cell_family[c_id];

  /* Increment mesh counts. */
  mesh->n_cells++;
  mesh->n_b_cells++;
  mesh->n_b_faces += f_incr;
  mesh->n_vertices += unique_vtx;
}

/*-----------------------------------------------------------------------------
 * Adds the cells lying in the positive/negative half-spaces to a group.
 * Useful for further mesh processing.
 *----------------------------------------------------------------------------*/

static void
_mark_cut_cells(cs_mesh_t *mesh, cs_lnum_t n_new_cells, const cs_lnum_t cells[])
{
  cs_mesh_group_cells_set(mesh, "auto:negative_cells", n_new_cells, cells);

  cs_lnum_t *sel_cells = nullptr;
  CS_MALLOC(sel_cells, n_new_cells, cs_lnum_t);
  cs_lnum_t n_c_ini = mesh->n_cells - n_new_cells;
  for (cs_lnum_t i = 0; i < n_new_cells; i++) {
    sel_cells[i] = n_c_ini + i;
  }
  cs_mesh_group_cells_set(mesh, "auto:positive_cells", n_new_cells, sel_cells);
  CS_FREE(sel_cells);
}

/*-----------------------------------------------------------------------------
 * Update a global numbering array.
 *
 * Since the algorithm only creates new unique entities, a simple scan is
 * enough.
 *
 * \param[in, out]  global_num  the global numbering array to update.
 * \param[in]       n_local     the local number of entities post-cut.
 * \param[in]       delta       the number of new entities created by the cut.
 * \param[in, out]  n_global    the global number of entities post-cut.
 *----------------------------------------------------------------------------*/

static void
_update_global_num(cs_gnum_t *global_num[],
                   cs_lnum_t  n_local,
                   cs_lnum_t  delta,
                   cs_gnum_t *n_global)
{
  cs_gnum_t local_new = (cs_gnum_t)delta;
  cs_gnum_t scan_new = local_new;
  cs_gnum_t total_new = local_new;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {
    MPI_Scan(&local_new, &scan_new, 1, CS_MPI_GNUM, MPI_SUM, cs_glob_mpi_comm);
    MPI_Allreduce(&local_new, &total_new, 1, CS_MPI_GNUM, MPI_SUM,
        cs_glob_mpi_comm);
  }

#endif // defined(HAVE_MPI)

  if (*global_num) {
    CS_REALLOC(*global_num, n_local, cs_gnum_t);
    cs_gnum_t *_global_num = *global_num;
    cs_gnum_t start = scan_new - local_new;

    for (cs_lnum_t i = 0; i < delta; i++)
      _global_num[n_local-delta+i] = *n_global + start + i + 1;
  }

  *n_global = *n_global + total_new;
}

static void
_update_global_vertices(cs_mesh_t *mesh, cs_lnum_t n_new_vertices)
{
  _update_global_num(&mesh->global_vtx_num,
                     mesh->n_vertices,
                     n_new_vertices,
                     &mesh->n_g_vertices);
}

static void
_update_global_b_faces(cs_mesh_t *mesh, cs_lnum_t n_new_b_faces)
{
  _update_global_num(&mesh->global_b_face_num,
                     mesh->n_b_faces,
                     n_new_b_faces,
                     &mesh->n_g_b_faces);
}

static void
_update_global_cells(cs_mesh_t *mesh, cs_lnum_t n_new_cells)
{
  _update_global_num(&mesh->global_cell_num,
                     mesh->n_cells,
                     n_new_cells,
                     &mesh->n_g_cells);
}

/*-----------------------------------------------------------------------------
 * Update parallel mesh structures and counts.
 *
 * Since the algorithm only creates new unique entities, a simple scan is
 * enough.
 *
 * \param[in, out]   mesh            pointer to mesh structure.
 * \param[in]        n_new_vertices  new number of local vertices
 * \param[in]        n_new_b_faces   new number of local boundary faces
 * \param[in]        n_new_cells     new number of local cells
 *----------------------------------------------------------------------------*/

static void
_update_parallelism(cs_mesh_t *mesh,
                    cs_lnum_t  n_new_vertices,
                    cs_lnum_t  n_new_b_faces,
                    cs_lnum_t  n_new_cells)
{
  mesh->n_cells_with_ghosts = mesh->n_cells;

  _update_global_vertices(mesh, n_new_vertices);
  _update_global_b_faces(mesh, n_new_b_faces);
  _update_global_cells(mesh, n_new_cells);

  if (cs_glob_n_ranks > 1) {
    cs_halo_type_t halo_type = mesh->halo_type;
    cs_mesh_builder_t *mb = (mesh == cs_glob_mesh) ?
                            cs_glob_mesh_builder :
                            nullptr;
    cs_mesh_init_halo(mesh, mb, halo_type, -1, true);
  }
}

/*-----------------------------------------------------------------------------
 * Prepare the cutting of cell faces.
 *
 * - Determine which cells are to be cut,
 * - Transform interior faces to boundary faces for those cells.
 *
 * \param[in, out]   mesh        pointer to mesh structure.
 * \param[in]        p_normals   new number of local vertices
 * \param[in, out]   n_c_cut     new number of local boundary faces
 * \param[out]       c_cut       new number of local cells
 *----------------------------------------------------------------------------*/

static void
_prepare_cut_cell_faces(cs_mesh_t        *mesh,
                        const cs_real_t   p_normals[][3],
                        cs_lnum_t        &n_c_cut,
                        cs_lnum_t         c_cut[])
{
  cs_real_t *p_norms = nullptr;
  CS_MALLOC(p_norms, mesh->n_cells_with_ghosts, cs_real_t);

  n_c_cut = 0;
  for (cs_lnum_t c_id = 0; c_id < mesh->n_cells; c_id++) {
    p_norms[c_id] = cs_math_3_square_norm(p_normals[c_id]);
    if (p_norms[c_id] > 0.0)
      c_cut[n_c_cut++] = c_id;
  }

  if (mesh->halo) {
    cs_halo_sync_untyped(mesh->halo,
                         mesh->halo_type,
                         sizeof(cs_real_t),
                         p_norms);
  }

  /* Transform the cut cells internal faces into boundary faces */
  cs_lnum_t *sel_faces = nullptr;
  CS_MALLOC(sel_faces, mesh->n_i_faces, cs_lnum_t);
  cs_lnum_t n_sel = 0;
  for (cs_lnum_t f_id = 0; f_id < mesh->n_i_faces; f_id++) {
    for (int j = 0; j < 2; j++) {
      if (p_norms[mesh->i_face_cells[f_id][j]] > 0.0) {
        sel_faces[n_sel++] = f_id;
        break;
      }
    }
  }

  cs_mesh_group_i_faces_set(mesh,
                            "auto:transformed_internal_faces",
                            n_sel,
                            sel_faces);

  cs_mesh_boundary_insert_with_shared_vertices(mesh,
                                               n_sel,
                                               sel_faces);

  CS_FREE(sel_faces);

  cs_mesh_free_rebuildable(mesh, true);

  // TODO: move this inside cs_mesh_free_rebuildable
  for (cs_lnum_t f_id = 0; f_id < mesh->n_i_faces; f_id++) {
    for (int j = 0; j < 2; j++) {
      if (mesh->i_face_cells[f_id][j] >= mesh->n_cells) {
        mesh->i_face_cells[f_id][j] = -1;
        break;
      }
    }
  }

  CS_FREE(p_norms);
}

/*------------------------------------------------------------------------------
 * Update the mesh face connectivity post cut.
 *----------------------------------------------------------------------------*/

static void
_update_face_connectivity(cs_mesh_t  *mesh,
                          cs_lnum_t   b_face_vtx_idx[],
                          cs_lnum_t   b_face_vtx_lst_size,
                          cs_lnum_t  *b_face_vtx_lst[])
{
  CS_FREE(mesh->b_face_vtx_idx);

  /* Counts to indices */
  for (cs_lnum_t i = 0; i < mesh->n_b_faces; i++) {
    b_face_vtx_idx[i+1] += b_face_vtx_idx[i];
  }
  mesh->b_face_vtx_idx = b_face_vtx_idx;
  mesh->b_face_vtx_connect_size = b_face_vtx_idx[mesh->n_b_faces];

  CS_FREE(mesh->b_face_vtx_lst);

  cs_lnum_t *_b_face_vtx_lst = nullptr;
  CS_MALLOC(_b_face_vtx_lst, mesh->b_face_vtx_connect_size, cs_lnum_t);
  cs_lnum_t j = 0;

  cs_lnum_t *lst = *b_face_vtx_lst;

  for (cs_lnum_t i = 0; i < b_face_vtx_lst_size; i++) {
    if (lst[i] != -1)
      _b_face_vtx_lst[j++] = lst[i];
  }

  assert(j == mesh->b_face_vtx_connect_size);

  CS_FREE(lst);

  mesh->b_face_vtx_lst = _b_face_vtx_lst;
}

/*-----------------------------------------------------------------------------
 * Allocate the new face-vertex connectivity arrays.
 *
 * \param[in, out]   mesh            pointer to mesh structure.
 * \param[out]       b_face_vtx_idx  boundary face->vertices index
 * \param[out]       b_face_vtx_lst  boundary face->vertices list
 * \param[in]        n_new_faces  number of new faces
 * \param[in]        e_stride     maximum number of edges per face.
 *
 * \return size of allocated b_face_vtx_lst connectivity array
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_init_b_face_connectivity_arrays(cs_mesh_t  *mesh,
                                 cs_lnum_t  *b_face_vtx_idx[],
                                 cs_lnum_t  *b_face_vtx_lst[],
                                 cs_lnum_t   n_new_faces,
                                 cs_lnum_t   e_stride)
{
  CS_MALLOC(*b_face_vtx_idx, mesh->n_b_faces+n_new_faces+1, cs_lnum_t);
  cs_lnum_t *idx = *b_face_vtx_idx;

  for (cs_lnum_t i = 0; i < mesh->n_b_faces+n_new_faces+1; i++) {
    idx[i] = 0;
  }
  for (cs_lnum_t i = 0; i < mesh->n_b_faces; i++) {
    idx[i+1] = mesh->b_face_vtx_idx[i+1] - mesh->b_face_vtx_idx[i];
  }

  cs_lnum_t b_face_vtx_lst_size = (mesh->n_b_faces+n_new_faces)*e_stride;
  CS_MALLOC(*b_face_vtx_lst, b_face_vtx_lst_size, cs_lnum_t);
  cs_lnum_t *ids = *b_face_vtx_lst;

  for (cs_lnum_t i = 0; i < b_face_vtx_lst_size; i++) ids[i] = -1;

  for (cs_lnum_t i = 0; i < mesh->n_b_faces; i++) {
    cs_lnum_t *_fv = ids + i*e_stride;

    cs_lnum_t start = mesh->b_face_vtx_idx[i];
    cs_lnum_t end = mesh->b_face_vtx_idx[i+1];
    cs_lnum_t nv = end - start;
    const cs_lnum_t *fv = mesh->b_face_vtx_lst + start;

    for (cs_lnum_t j = 0; j < e_stride && j < nv; j++) {
      _fv[j] = fv[j];
    }
  }

  return b_face_vtx_lst_size;
}

/*------------------------------------------------------------------------------
 * Resize the mesh connectivity arrays before mesh modification.
 *----------------------------------------------------------------------------*/

static void
_resize_mesh(cs_mesh_t  *mesh,
             cs_lnum_t   n_new_vertices,
             cs_lnum_t   n_new_faces,
             cs_lnum_t   n_new_cells)
{
  CS_REALLOC(mesh->vtx_coord, 3*(mesh->n_vertices+n_new_vertices), cs_real_t);

  CS_REALLOC(mesh->b_face_cells, mesh->n_b_faces+n_new_faces, cs_lnum_t);

  for (cs_lnum_t i = 0; i < n_new_faces; i++)
    mesh->b_face_cells[mesh->n_b_faces+i] = -1;

  CS_REALLOC(mesh->cell_family, mesh->n_cells+n_new_cells, cs_lnum_t);

  for (cs_lnum_t i = mesh->n_cells; i < mesh->n_cells+n_new_cells; i++)
    mesh->cell_family[i] = _default_family_id;

  CS_REALLOC(mesh->b_face_family, mesh->n_b_faces+n_new_faces, cs_lnum_t);

  for (cs_lnum_t i = mesh->n_b_faces; i < mesh->n_b_faces+n_new_faces; i++)
    mesh->b_face_family[i] = _default_family_id;
}

/*------------------------------------------------------------------------------
 * Free scratch data used for cell cut algorithm.
 *
 * \param[in, out]  cd  cut helper structure whose members are to be freed.
 *----------------------------------------------------------------------------*/

static void
_free_cut_data(_cut_data  *cd)
{
  CS_FREE(cd->sd);
  CS_FREE(cd->occurs);
  CS_FREE(cd->xyz);
  CS_FREE(cd->indices);
  CS_FREE(cd->compact);
  CS_FREE(cd->polys);
  CS_FREE(cd->f2e);
  CS_FREE(cd->e2f);
  CS_FREE(cd->edges);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Cut cells with planes.
 *
 * Each cell can be cut by a single plane, defined by its normal and an origin
 * (i.e. any point in the plane). Cells whose assigned normals are null
 * vectors are not cut.
 *
 * The polygons created by the cut are added to a new group,
 * "auto:closing_polygons".
 *
 * This function should be followed by applying a joining on the group,
 * "auto:transformed_internal_faces".
 *
 * \param[in, out]  mesh      mesh to cut
 * \param[in]       p_normals array of plane_normals of size mesh->n_cells
 * \param[in]       p_origins array of plane origins of size mesh->n_cells
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cut(cs_mesh_t       *mesh,
            const cs_real_t  p_normals[][3],
            const cs_real_t  p_origins[][3])
{
  std::chrono::high_resolution_clock::time_point t0
    = std::chrono::high_resolution_clock::now();

  bft_printf("\nStart of cell-plane cut\n\n");

  cs_lnum_t *sel_cells = nullptr;
  CS_MALLOC(sel_cells, mesh->n_cells, cs_lnum_t);
  cs_lnum_t n_new_cells;

  _prepare_cut_cell_faces(mesh, p_normals, n_new_cells, sel_cells);

  /* Create useful cell connectivity arrays. */
  const cs_lnum_t n_c_ini = mesh->n_cells;
  const cs_lnum_t n_v_ini = mesh->n_vertices;
  const cs_lnum_t n_b_ini = mesh->n_b_faces;
  const cs_adjacency_t *c2v = cs_mesh_adjacencies_cell_vertices();
  cs_adjacency_t *c2f_b = cs_mesh_adjacency_c2f_boundary(mesh);

  /* Size-up the problem. */
  cs_lnum_t max_nv = 0, max_nf = 0;
  for (cs_lnum_t f_id = 0; f_id < mesh->n_b_faces; f_id++) {
    cs_lnum_t c_id = mesh->b_face_cells[f_id];
    cs_lnum_t nf = c2f_b->idx[c_id+1] - c2f_b->idx[c_id];
    if (nf > max_nf) max_nf = nf;

    cs_lnum_t nv = mesh->b_face_vtx_idx[f_id+1] - mesh->b_face_vtx_idx[f_id];
    if (nv > max_nv) max_nv = nv;
  }

  cs_lnum_t n_new_faces = n_new_cells*(max_nf + 2);
  cs_lnum_t n_new_vertices = (n_new_faces - 2*n_new_cells)*2;

  _cut_data cd = _allocate_scratch_data(mesh,
                                        max_nf,
                                        max_nv,
                                        n_new_vertices,
                                        n_new_cells);

  /* Init b_face connectivity arrays. */
  cs_lnum_t *b_face_vtx_idx = nullptr;
  cs_lnum_t *b_face_vtx_lst = nullptr;
  cs_lnum_t b_face_vtx_lst_size;
  b_face_vtx_lst_size = _init_b_face_connectivity_arrays(mesh,
                                                         &b_face_vtx_idx,
                                                         &b_face_vtx_lst,
                                                         n_new_faces,
                                                         cd.e_stride);

  /* Resize mesh. */
  _resize_mesh(mesh, n_new_vertices, n_new_faces, n_new_cells);

  std::chrono::high_resolution_clock::time_point t1
    = std::chrono::high_resolution_clock::now();
  std::chrono::microseconds te_prepare
    = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);

  /* Cut. */
  for (cs_lnum_t i = 0; i < n_new_cells; i++) {
    cs_lnum_t c_id = sel_cells[i];
    const cs_real_t *p_normal = p_normals[c_id];
    const cs_real_t *p_origin = p_origins[c_id];
    _cut_cell(mesh,
              &cd,
              c_id,
              p_normal,
              p_origin,
              c2v,
              c2f_b,
              b_face_vtx_idx,
              b_face_vtx_lst);
  }

  std::chrono::high_resolution_clock::time_point t2
    = std::chrono::high_resolution_clock::now();
  std::chrono::microseconds te_cut
    = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);

  n_new_cells = mesh->n_cells - n_c_ini;

  /* Add the closing polygons to a group. */
  cs_mesh_group_b_faces_set(mesh,
                            "auto:closing_polygons",
                            cd.n_polys,
                            cd.polys);

  /* Free cut data. */
  _free_cut_data(&cd);
  cs_adjacency_destroy(&c2f_b);
  cs_mesh_adjacencies_finalize();

  /* Set up the end b_face-vertex connectivity arrays. */
  _update_face_connectivity(mesh,
                            b_face_vtx_idx,
                            b_face_vtx_lst_size,
                            &b_face_vtx_lst);

  /* Optional: mark the cut cells for further post-processing. */
  _mark_cut_cells(mesh, n_new_cells, sel_cells);
  CS_FREE(sel_cells);

  /* Update parallel data. */
  n_new_vertices = mesh->n_vertices - n_v_ini;
  n_new_faces = mesh->n_b_faces - n_b_ini;
  _update_parallelism(mesh, n_new_vertices, n_new_faces, n_new_cells);
  cs_mesh_update_auxiliary(mesh);

  /* Tag mesh for repartitionning. */
  mesh->modified |= CS_MESH_MODIFIED;
  mesh->modified |= CS_MESH_MODIFIED_BALANCE;

  bft_printf("\nEnd of cell-plane cut.\n");

  t1 = std::chrono::high_resolution_clock::now();
  std::chrono::microseconds te_update
    = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t2);

  if (mesh->verbosity > 0) {

    cs_mesh_print_element_counts(mesh,
                                 _("Mesh after cells cut"));

    cs_log_printf(CS_LOG_DEFAULT, "\n");
    cs_log_separator(CS_LOG_DEFAULT);

    cs_log_printf
      (CS_LOG_PERFORMANCE,
       _("\nMesh cells cut:\n\n"
         "  Preparation:                                  %.3g\n"
         "  Cells cut:                                    %.3g\n"
         "  Mesh update:                                  %.3g\n"),
       (double)(te_prepare.count()*1.e-6),
       (double)(te_cut.count()*1.e-6),
       (double)(te_update.count()*1.e-6));
    cs_log_printf(CS_LOG_PERFORMANCE, "\n");
    cs_log_separator(CS_LOG_PERFORMANCE);

  }
}

/*----------------------------------------------------------------------------*/
