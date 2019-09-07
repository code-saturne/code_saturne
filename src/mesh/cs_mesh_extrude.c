/*============================================================================
 * Mesh extrusion.
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

#include <float.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_io_num.h"

#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_order.h"
#include "cs_parall.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_extrude.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_mesh_extrude.c
        Mesh extrusion.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#undef _CS_MODULE

#define _CS_MODULE(vect) \
  sqrt(vect[0] * vect[0] + vect[1] * vect[1] + vect[2] * vect[2])

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Determine angle between two face edges adjacent to a given vertex.
 *
 * \param[in]       v_ids          vertex ids (e0v0, e0v1, e1v1)
 * \param[in]       f_n            face normal
 * \param[in]       vertex_coords  vertex coordinates
 *
 * \return angle
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_angle(const cs_lnum_t  v_ids[3],
       const cs_real_t  f_n[3],
       const cs_real_t  vtx_coords[][3])
{
  const cs_real_t *c0 = vtx_coords[v_ids[0]];
  const cs_real_t *c1 = vtx_coords[v_ids[1]];
  const cs_real_t *c2 = vtx_coords[v_ids[2]];

  const cs_real_t u[3] = {c1[0]-c0[0], c1[1]-c0[1], c1[2]-c0[2]};
  const cs_real_t v[3] = {c2[0]-c0[0], c2[1]-c0[1], c2[2]-c0[2]};

  const cs_real_t d    = cs_math_3_dot_product(u, v);
  const cs_real_t lsq0 = cs_math_3_square_norm(u);
  const cs_real_t lsq1 = cs_math_3_square_norm(v);

  double r = d/sqrt(lsq0*lsq1);
  if (r > 1)
    r = 1;
  else if (r < -1)
    r = -1;

  double theta = acos(r);

  cs_real_t uv[3];

  cs_math_3_cross_product(u, v, uv);

  /* Check the sign */
  if (cs_math_3_dot_product(uv, f_n) < 0)
    theta = 2.*cs_math_pi - theta;

  return theta;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Generate a vertices selection from adjacent boundary faces.
 *
 * The calling code is responsible for freeing the allocated vertices
 * array when not needed anymore.
 *
 * \param[in, out]  m           mesh
 * \param[in]       n_b_faces   number of selected boundary faces
 * \param[in]       b_faces     list of selected boundary faces (0 to n-1),
 *                              or NULL if no indirection is needed
 * \param[out]      n_vertices  number of selected vertices
 * \param[out]      vertices    list of vertices (0 to n-1)
 */
/*----------------------------------------------------------------------------*/

static void
_select_vertices_from_adj_b_faces(const cs_mesh_t   *m,
                                  cs_lnum_t          n_b_faces,
                                  const cs_lnum_t    b_faces[],
                                  cs_lnum_t         *n_vertices,
                                  cs_lnum_t        **vertices)
{
  cs_lnum_t _n_vertices = 0;
  cs_lnum_t *_vertices = NULL;

  char *v_flag = NULL;

  /* Build vertices indirection */

  BFT_MALLOC(v_flag, m->n_vertices, char);

  for (cs_lnum_t i = 0; i < m->n_vertices; i++)
    v_flag[i] = 0;

  /* Mark based on selected faces */

  if (b_faces != NULL) {
    for (cs_lnum_t i = 0; i < n_b_faces; i++) {
      cs_lnum_t f_id = b_faces[i];
      cs_lnum_t s_id = m->b_face_vtx_idx[f_id];
      cs_lnum_t e_id = m->b_face_vtx_idx[f_id+1];
      for (cs_lnum_t j = s_id; j < e_id; j++)
        v_flag[m->b_face_vtx_lst[j]] = 1;
    }
  }
  else {
    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      cs_lnum_t s_id = m->b_face_vtx_idx[f_id];
      cs_lnum_t e_id = m->b_face_vtx_idx[f_id+1];
      for (cs_lnum_t j = s_id; j < e_id; j++)
        v_flag[m->b_face_vtx_lst[j]] = 1;
    }
  }

  /* Synchronize across ranks if necessary */

  if (m->vtx_interfaces != NULL)
    cs_interface_set_max(m->vtx_interfaces,
                         m->n_vertices,
                         1,
                         true,
                         CS_CHAR,
                         v_flag);

  /* Transform to list */

  for (cs_lnum_t i = 0; i < m->n_vertices; i++) {
    if (v_flag[i] != 0)
      _n_vertices++;
  }

  BFT_MALLOC(_vertices, _n_vertices, cs_lnum_t);

  _n_vertices = 0;
  for (cs_lnum_t i = 0; i < m->n_vertices; i++) {
    if (v_flag[i] != 0) {
      _vertices[_n_vertices] = i;
      _n_vertices++;
    }
  }

  /* Free temporary arrays and return result array */

  BFT_FREE(v_flag);

  *n_vertices = _n_vertices;
  *vertices = _vertices;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compact indexed array, removing unused entries
 *
 * This function is used when an initial array is oversized, and has been
 * built with a first pass.
 *
 * Arrays are not freed or reallocated here, so the caller is responsible
 * for freeing the count array and possibly reallocating data to the final size
 *
 * \param[in]       elt_size    element size
 * \param[in]       n_elts      number of elements
 * \param[in, out]  idx         index of data to compact (size: n_elts + 1)
 * \param[in, out]  count       count of used data entries (size: n_elts)
 * \param[in,out]   data        indexed data
 *
 * \return  new size of data
 */
/*----------------------------------------------------------------------------*/

static cs_lnum_t
_compact_indexed_data(size_t      elt_size,
                      cs_lnum_t   n_elts,
                      cs_lnum_t  *idx,
                      cs_lnum_t  *count,
                      void       *data)
{
  cs_lnum_t data_size = 0;

  unsigned char *_data = data;

  for (cs_lnum_t i = 0; i < n_elts; i++) {
    cs_lnum_t s_id0 = idx[i]*elt_size;
    cs_lnum_t s_id1 = (idx[i] + count[i])*elt_size;
    cs_lnum_t d_p = data_size*elt_size;
    idx[i] = data_size;
    data_size += count[i];
    for (cs_lnum_t j = s_id0; j < s_id1; j++)
      _data[d_p++] = _data[j];
  }
  idx[n_elts] = data_size;

  return data_size;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build edges adjacent to selected boundary faces
 *
 * Marked vertices (belonging to selected faces) are used here rather
 * than faces, so as to also build edges adjacent to selected faces,
 * even on other ranks (the marking of vertices having been done
 * and synchronized previously).
 *
 * The calling code is responsible for freeing the allocated arrays
 * when not needed anymore.
 *
 * \param[in, out]  m             mesh
 * \param[in]       n_faces       number of selected boundary faces
 * \param[in]       n_vertices    number of selected vertices
 * \param[in]       faces         list of selected boundary faces (0 to n-1)
 * \param[in]       vertices      ids of selected vertices (0 to n-1)
 * \param[in]       n_layers      number of layers per selected vertex
 * \param[in]       vertices      ids of selected vertices (0 to n-1)
 * \param[out]      n_i_edges     number of interior edges
 * \param[out]      n_b_edges     number of boundary edges
 * \param[out]      i_e2sf        interior edges to selected faces
 *                                connectivity
 * \param[out]      i_e2v         interior edge to vertices connectivity
 * \param[out]      i_e_io_num    interior edge global number (or NULL)
 * \param[out]      b_e2sf        boundary edges to selected faces
 *                                connectivity
 * \param[out]      b_e2v         boundary edge to vertices connectivity
 * \param[out]      b_e_gc        boundary edge group class
 * \param[out]      b_e_io_num    boundary edge global number (or NULL)
 *
 * \return  pointer modified number of layers par selected vertex array if
*           locking of some vertices is required, NULL otherwise
 */
/*----------------------------------------------------------------------------*/

static cs_lnum_t *
_build_face_edges(cs_mesh_t         *m,
                  cs_lnum_t          n_faces,
                  cs_lnum_t          n_vertices,
                  const cs_lnum_t    faces[],
                  const cs_lnum_t    vertices[],
                  const cs_lnum_t    n_layers[],
                  cs_lnum_t         *n_i_edges,
                  cs_lnum_t         *n_b_edges,
                  cs_lnum_2_t       *i_e2sf[],
                  cs_lnum_2_t       *i_e2v[],
                  fvm_io_num_t     **i_e_io_num,
                  cs_lnum_t         *b_e2sf[],
                  cs_lnum_2_t       *b_e2v[],
                  int               *b_e_gc[],
                  fvm_io_num_t     **b_e_io_num)
{
  const int default_family_id = 1;

  /* Build vertices flag */

  char *v_flag;
  BFT_MALLOC(v_flag, m->n_vertices, char);
  for (cs_lnum_t i = 0; i < m->n_vertices; i++)
    v_flag[i] = 0;
  for (cs_lnum_t i = 0; i < n_vertices; i++)
    v_flag[vertices[i]] = 1;

  /* Build vertices indirection */

  cs_lnum_t *v2v_idx;
  BFT_MALLOC(v2v_idx, m->n_vertices+1, cs_lnum_t);

  for (cs_lnum_t i = 0; i < m->n_vertices+1; i++)
    v2v_idx[i] = 0;

  /* Counting loop for b_f2e_idx */

  for (cs_lnum_t f_id = 0; f_id < m->n_b_faces; f_id++) {
    cs_lnum_t s_id = m->b_face_vtx_idx[f_id];
    cs_lnum_t n_f_vtx = m->b_face_vtx_idx[f_id+1] - s_id;
    for (cs_lnum_t j = 0; j < n_f_vtx; j++) {
      cs_lnum_t vid0 = m->b_face_vtx_lst[s_id + j];
      cs_lnum_t vid1 = m->b_face_vtx_lst[s_id + ((j+1)%n_f_vtx)];
      if (v_flag[vid0] && v_flag[vid1]) {
        if (vid0 < vid1)
          v2v_idx[vid0+1] += 1;
        else
          v2v_idx[vid1+1] += 1;
      }
    }
  }

  /* Convert b_f2e_idx to true index */

  for (cs_lnum_t i = 0; i < m->n_vertices; i++)
    v2v_idx[i+1] += v2v_idx[i];

  /* Now build index and edges */

  cs_lnum_t *v2v, *v2v_count;

  BFT_MALLOC(v2v, v2v_idx[m->n_vertices], cs_lnum_t);
  BFT_MALLOC(v2v_count, m->n_vertices, cs_lnum_t);

  for (cs_lnum_t i = 0; i < m->n_vertices; i++)
    v2v_count[i] = 0;

  for (cs_lnum_t f_id = 0; f_id < m->n_b_faces; f_id++) {
    cs_lnum_t s_id = m->b_face_vtx_idx[f_id];
    cs_lnum_t n_f_vtx = m->b_face_vtx_idx[f_id+1] - s_id;
    for (cs_lnum_t j = 0; j < n_f_vtx; j++) {
      cs_lnum_t vid0 = m->b_face_vtx_lst[s_id + j];
      cs_lnum_t vid1 = m->b_face_vtx_lst[s_id + ((j+1)%n_f_vtx)];
      if (v_flag[vid0] && v_flag[vid1]) {
        if (vid0 > vid1) { /* swap vertices, smallest id first */
          cs_lnum_t vid_tmp = vid1;
          vid1 = vid0;
          vid0 = vid_tmp;
        }
        /* search for couple in index */
        cs_lnum_t k = v2v_idx[vid0], l = 0;
        while (l < v2v_count[vid0]) {
          if (v2v[k] == vid1)
            break;
          else
            k++, l++;
        }
        /* add couple if not previously found */
        if (l >= v2v_count[vid0]) {
          v2v[k] = vid1;
          v2v_count[vid0] += 1;
        }
      }
    }
  }

  BFT_FREE(v_flag);

  /* Make indexed structure compact */

  cs_lnum_t _n_edges = _compact_indexed_data(sizeof(cs_lnum_t),
                                             m->n_vertices,
                                             v2v_idx,
                                             v2v_count,
                                             v2v);

  BFT_REALLOC(v2v, _n_edges, cs_lnum_t);

  BFT_FREE(v2v_count);

  /* Generate boundary face -> selected face indirection */

  cs_lnum_t *f_s_id;
  BFT_MALLOC(f_s_id, m->n_b_faces, cs_lnum_t);

  for (cs_lnum_t f_id = 0; f_id < m->n_b_faces; f_id++)
    f_s_id[f_id] = -1;

  for (cs_lnum_t i = 0; i < n_faces; i++)
    f_s_id[faces[i]] = i;

  /* Now build edges;
     ---------------

     we will define edges by their adjacent faces, where e2sf[i][0] is
     the boundary face for which the directed edge is seen in the "positive"
     direction when looping over the face's edges in the trigonometric
     direction, and e2sf[i][1] is the face for which the edge is seen
     in the negative direction. If no matching face is in the selection,
     e2sf[i][j] = -1.

     Also, we want to base group classes (families) on the faces
     just outside the selection (to prolong those naturally),
     so the group class is added only for edges visited from
     faces which are not selected.
  */

  cs_lnum_2_t *e2v, *e2sf;
  int *e_gc, *e_nf;

  BFT_MALLOC(e2v, _n_edges, cs_lnum_2_t);
  BFT_MALLOC(e2sf, _n_edges, cs_lnum_2_t);
  BFT_MALLOC(e_gc, _n_edges, int);
  BFT_MALLOC(e_nf, _n_edges, int);

  for (cs_lnum_t i = 0; i < _n_edges; i++) {
    e2v[i][0] = -1;
    e2v[i][1] = -1;
    e2sf[i][0] = -1;
    e2sf[i][1] = -1;
    e_gc[i] = default_family_id;
    e_nf[i] = 0;
  }

  for (cs_lnum_t f_id = 0; f_id < m->n_b_faces; f_id++) {
    cs_lnum_t s_id = m->b_face_vtx_idx[f_id];
    cs_lnum_t n_f_vtx = m->b_face_vtx_idx[f_id+1] - s_id;
    for (cs_lnum_t j = 0; j < n_f_vtx; j++) {
      cs_lnum_t vid0 = m->b_face_vtx_lst[s_id + j];
      cs_lnum_t vid1 = m->b_face_vtx_lst[s_id + ((j+1)%n_f_vtx)];
      int e_orient;
      if (vid0 > vid1) { /* swap vertices, smallest id first */
        cs_lnum_t vid_tmp = vid1;
        vid1 = vid0;
        vid0 = vid_tmp;
        e_orient = 1;
      }
      else
        e_orient = 0;
      /* search for couple in index */
      for (cs_lnum_t k = v2v_idx[vid0]; k < v2v_idx[vid0+1]; k++) {
        if (v2v[k] == vid1) {
          e2v[k][0] = vid0;
          e2v[k][1] = vid1;
          e2sf[k][e_orient] = f_s_id[f_id];
          if (f_s_id[f_id] > -1)
            e_nf[k] += 1;
          else
            e_gc[k] = m->b_face_family[f_id];
          break;
        }
      }
    }
  }

  /* Free work arrays */

  BFT_FREE(v2v);
  BFT_FREE(v2v_idx);

  BFT_FREE(f_s_id);

  /* Build initial global numbering if applicable */
  /* -------------------------------------------- */

  cs_gnum_t *e_gnum = NULL;

  if (m->global_vtx_num != NULL || cs_glob_n_ranks > 1) {

    const cs_gnum_t  *global_vtx_num = m->global_vtx_num;
    cs_gnum_t *g_e2v;

    BFT_MALLOC(g_e2v, _n_edges*2, cs_gnum_t);

    /* Orient and order edges based on global (not local) numbers */

    for (cs_lnum_t i = 0; i < _n_edges; i++) {
      cs_gnum_t gv0 = global_vtx_num[e2v[i][0]];
      cs_gnum_t gv1 = global_vtx_num[e2v[i][1]];
      if (gv1 < gv0) {
        cs_lnum_t vid_tmp = e2v[i][0];
        e2v[i][0] = e2v[i][1];
        e2v[i][1] = vid_tmp;
        cs_lnum_t fid_tmp = e2sf[i][0];
        e2sf[i][0] = e2sf[i][1];
        e2sf[i][1] = fid_tmp;
        cs_gnum_t gvid_tmp = gv0;
        gv0 = gv1;
        gv1 = gvid_tmp;
      }
      g_e2v[i*2]     = gv0;
      g_e2v[i*2 + 1] = gv1;
    }

    cs_lnum_t *order = cs_order_gnum_s(NULL, g_e2v, 2, _n_edges);

    cs_order_reorder_data(_n_edges, 2*sizeof(cs_gnum_t), order, g_e2v);
    cs_order_reorder_data(_n_edges, 2*sizeof(cs_lnum_t), order, e2v);
    cs_order_reorder_data(_n_edges, 2*sizeof(cs_lnum_t), order, e2sf);
    cs_order_reorder_data(_n_edges, sizeof(int), order, e_gc);
    cs_order_reorder_data(_n_edges, sizeof(int), order, e_nf);

    BFT_FREE(order);

    fvm_io_num_t  *e_io_num
      = fvm_io_num_create_from_adj_s(NULL, g_e2v, _n_edges, 2);

    e_gnum = fvm_io_num_transfer_global_num(e_io_num);

    e_io_num = fvm_io_num_destroy(e_io_num);

    BFT_FREE(g_e2v);

    /* No need to handle periodicity here, as periodic edges could
       have independent group classes. */

    cs_interface_set_t *e_if
      = cs_interface_set_create(_n_edges, NULL, e_gnum,
                                NULL, 0, NULL, NULL, NULL);

    const cs_datatype_t int_type
      = (sizeof(int) == 8) ? CS_INT64 : CS_INT32;

    cs_interface_set_sum(e_if, _n_edges, 1, true, int_type, e_nf);

    cs_interface_set_destroy(&e_if);

  }

  /* Separate interior from boundary edges */
  /*--------------------------------------*/

  /* Note that edges with no local adjacent face (which can only
     be edges both on the selection boundary and on a parallel or
     periodic boundary) are discarded here, as they are not needed
     anymore. */

  /* Counting loop */

  cs_lnum_t *lock_mult = NULL;

  cs_lnum_t _n_i_edges = 0, _n_b_edges = 0;

  for (cs_lnum_t i = 0; i < _n_edges; i++) {
    if (e_nf[i] == 2)
      _n_i_edges += 1;
    else if (e_nf[i] == 1) {
      /* Ignore boundary edge on rank not owning the adjacent face */
      if (e2sf[i][0] + e2sf[i][1] < -1)
        e_nf[i] = 0;
      else
        _n_b_edges += 1;
    }
    else if (e_nf[i] > 2) {
      /* Too many incident faces: no extrusion along this edge */
      if (lock_mult == NULL) {
        BFT_MALLOC(lock_mult,  m->n_vertices, cs_lnum_t);
        for (cs_lnum_t j = 0; j < m->n_vertices; j++)
          lock_mult[j] = 1;
      }
      lock_mult[e2v[i][0]] = 0;
      lock_mult[e2v[i][1]] = 0;
    }
  }

  cs_lnum_t     *_b_e2sf;
  cs_lnum_2_t   *_b_e2v;
  int           *_b_e_gc;
  BFT_MALLOC(_b_e2sf, _n_b_edges, cs_lnum_t);
  BFT_MALLOC(_b_e2v, _n_b_edges, cs_lnum_2_t);
  BFT_MALLOC(_b_e_gc, _n_b_edges, int);

  /* Distribution loop */

  _n_i_edges = 0; _n_b_edges = 0;

  for (cs_lnum_t i = 0; i < _n_edges; i++) {
    if (e_nf[i] == 2) {
      e2sf[_n_i_edges][0] = e2sf[i][0];
      e2sf[_n_i_edges][1] = e2sf[i][1];
      e2v[_n_i_edges][0] = e2v[i][0];
      e2v[_n_i_edges][1] = e2v[i][1];
      _n_i_edges += 1;
    }
    else if (e_nf[i] == 1) {
      /* Reorient boundary edge if needed */
      if (e2sf[i][0] < 0) {
        assert(e2sf[i][1] >= 0);
        _b_e2sf[_n_b_edges] = e2sf[i][1];
        _b_e2v[_n_b_edges][0] = e2v[i][1];
        _b_e2v[_n_b_edges][1] = e2v[i][0];
      }
      else {
        _b_e2sf[_n_b_edges] = e2sf[i][0];
        _b_e2v[_n_b_edges][0] = e2v[i][0];
        _b_e2v[_n_b_edges][1] = e2v[i][1];
      }
      _b_e_gc[_n_b_edges] = e_gc[i];
      _n_b_edges += 1;
    }
  }

  BFT_REALLOC(e2sf, _n_i_edges, cs_lnum_2_t);
  BFT_REALLOC(e2v, _n_i_edges, cs_lnum_2_t);
  BFT_FREE(e_gc);

  /* Compact global numberings if required */

  if (m->global_vtx_num != NULL || cs_glob_n_ranks > 1) {

    cs_lnum_t ki = 0, kb = 0;
    cs_lnum_t *i_edges, *b_edges;
    BFT_MALLOC(i_edges, _n_i_edges, cs_lnum_t);
    BFT_MALLOC(b_edges, _n_b_edges, cs_lnum_t);

    for (cs_lnum_t i = 0; i < _n_edges; i++) {
      if (e_nf[i] == 2)
        i_edges[ki++] = i;
      else if (e_nf[i] == 1)
        b_edges[kb++] = i;
    }

    *i_e_io_num  = fvm_io_num_create_from_select(i_edges,
                                                 e_gnum,
                                                 _n_i_edges,
                                                 0);

    *b_e_io_num  = fvm_io_num_create_from_select(b_edges,
                                                 e_gnum,
                                                 _n_b_edges,
                                                 0);

    BFT_FREE(i_edges);
    BFT_FREE(b_edges);

  }
  else {

    *i_e_io_num = NULL;
    *b_e_io_num = NULL;

  }

  BFT_FREE(e_gnum);
  BFT_FREE(e_nf);

  /* Other return values */

  *n_i_edges = _n_i_edges;
  *i_e2sf = e2sf;
  *i_e2v = e2v;

  *n_b_edges = _n_b_edges;
  *b_e2sf = _b_e2sf;
  *b_e2v = _b_e2v;
  *b_e_gc = _b_e_gc;

  /* Restrict _n_layers to selected vertices */

  cs_lnum_t  *_n_layers = NULL;

  if (lock_mult != NULL) {
    BFT_MALLOC(_n_layers, n_vertices, cs_lnum_t);
    for (cs_lnum_t j = 0; j < n_vertices; j++) {
      _n_layers[j] = n_layers[j] * lock_mult[vertices[j]];
    }
    BFT_FREE(lock_mult);
  }

  return _n_layers;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create extruded vertices.
 *
 * Extrusion is defined on selected boundary faces, and the number of layers
 * for each associated vertex may be (slightly) variable, to account for
 * cluttered areas where extrusion may be constrained, or more complex
 * extrusions.
 *
 * \param[in, out]  m             mesh
 * \param[in]       n_vertices    number of selected vertices
 * \param[in]       vertices      ids of selected vertices (0 to n-1)
 * \param[in]       n_layers      number of layers for each vertex
 * \param[in]       n_layers_ini  predicted number of layers for each vertex
 * \param[in]       coord_shift   extrusion vector for each vertex
 * \param[in]       distribution  optional distribution of resulting vertices
 *                                along each extrusion vector
 *                                (size: n_cumulative_layers) with values
 *                                in range ]0, 1].
 *
 * \return  array defining the starting position of added layer vertices for
 *          selected vertex.
 */
/*----------------------------------------------------------------------------*/

static cs_lnum_t *
_add_extruded_vertices(cs_mesh_t          *m,
                       cs_lnum_t           n_vertices,
                       const cs_lnum_t     vertices[],
                       const cs_lnum_t     n_layers[],
                       const cs_lnum_t     n_layers_ini[],
                       const cs_coord_3_t  coord_shift[],
                       const float         distribution[])
{
  const cs_lnum_t n_vertices_ini = m->n_vertices;

  /* Count added vertices */

  cs_lnum_t *v_shift;

  BFT_MALLOC(v_shift, n_vertices+1, cs_lnum_t);

  v_shift[0] = 0;
  for (cs_lnum_t i = 0; i < n_vertices; i++)
    v_shift[i+1] = v_shift[i] + n_layers[i];

  cs_lnum_t n_vertices_add = v_shift[n_vertices];

  /* Update coordinates */

  BFT_REALLOC(m->vtx_coord, (n_vertices_ini + n_vertices_add)*3, cs_real_t);

  if (distribution != NULL) {
    for (cs_lnum_t i = 0; i < n_vertices; i++) {
      cs_lnum_t v_id = vertices[i];
      const cs_real_t *s_coo = m->vtx_coord + 3*v_id;
      const cs_lnum_t s_id = v_shift[i];
      cs_real_t *d_coo = m->vtx_coord + (n_vertices_ini + s_id)*3;
      for (cs_lnum_t j = 0; j < n_layers[i]; j++) {
        for (cs_lnum_t k = 0; k < 3; k++)
          d_coo[j*3 + k] = s_coo[k] + distribution[s_id + j]*coord_shift[i][k];
      }
    }
  }

  else {
    for (cs_lnum_t i = 0; i < n_vertices; i++) {
      cs_lnum_t v_id = vertices[i];
      const cs_real_t *s_coo = m->vtx_coord + 3*v_id;
      const cs_lnum_t s_id = v_shift[i];
      cs_real_t *d_coo = m->vtx_coord + (n_vertices_ini + s_id)*3;
      cs_real_t d_f = 1./n_layers_ini[i];
      for (cs_lnum_t j = 0; j < n_layers[i]; j++) {
        for (cs_lnum_t k = 0; k < 3; k++)
          d_coo[j*3 + k] = s_coo[k] + d_f*(j+1)*coord_shift[i][k];
      }
    }
  }

  /* Update global numbering */

  if (m->global_vtx_num != NULL || cs_glob_n_ranks > 1) {

    fvm_io_num_t *v_io_num
      = fvm_io_num_create_from_select(vertices,
                                      m->global_vtx_num,
                                      n_vertices,
                                      0);

    fvm_io_num_t *v_add_io_num
      = fvm_io_num_create_from_sub(v_io_num, n_layers);

    v_io_num = fvm_io_num_destroy(v_io_num);

    cs_gnum_t  v_add_gcount = fvm_io_num_get_global_count(v_add_io_num);
    const cs_gnum_t *v_add_gnum = fvm_io_num_get_global_num(v_add_io_num);

    assert(n_vertices_add == fvm_io_num_get_local_count(v_add_io_num));

    BFT_REALLOC(m->global_vtx_num, n_vertices_ini + n_vertices_add, cs_gnum_t);

    for (cs_lnum_t i = 0; i < n_vertices_add; i++)
      m->global_vtx_num[n_vertices_ini + i] = v_add_gnum[i] + m->n_g_vertices;

    v_add_io_num = fvm_io_num_destroy(v_add_io_num);

    m->n_g_vertices += v_add_gcount;

  }
  else
    m->n_g_vertices += n_vertices_add;

  m->n_vertices += n_vertices_add;

  /* Destroy vertex interfaces for future rebuild
     FIXME: handle periodicity */

  if (m->vtx_interfaces)
    cs_interface_set_destroy(&(m->vtx_interfaces));

  return v_shift;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create extruded cells.
 *
 * Cells inherit the group class of cells adjacent to selected faces.
 *
 * \param[in, out]  m             mesh
 * \param[in]       n_faces       number of selected boundary faces
 * \param[in]       faces         list of selected boundary faces (0 to n-1)
 * \param[in]       n_v_sub       number of layers for each vertex in mesh
 *                                (O for vertices not on selected faces)
 *
 * \return  array defining the starting position of added layer cells for
 *          selected face.
 */
/*----------------------------------------------------------------------------*/

static cs_lnum_t *
_add_extruded_cells(cs_mesh_t          *m,
                    cs_lnum_t           n_faces,
                    const cs_lnum_t     faces[],
                    const cs_lnum_t     n_v_sub[])
{
  const cs_lnum_t n_cells_ini = m->n_cells;

  /* Determine number of extruded cells per face */

  cs_lnum_t *c_shift;
  int *c_family;

  BFT_MALLOC(c_shift, n_faces+1, cs_lnum_t);
  BFT_MALLOC(c_family, n_faces, int);

  c_shift[0] = 0;

  for (cs_lnum_t i = 0; i < n_faces; i++) {
    cs_lnum_t f_id = faces[i];
    cs_lnum_t s_id = m->b_face_vtx_idx[f_id];
    cs_lnum_t e_id = m->b_face_vtx_idx[f_id+1];
    cs_lnum_t ns = 0;
    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t vid = m->b_face_vtx_lst[j];
      if (n_v_sub[vid] > ns)
        ns = n_v_sub[vid];
    }
    c_shift[i+1] = ns;
    c_family[i] = m->cell_family[m->b_face_cells[f_id]];
  }

  /* Transform count to index */

  for (cs_lnum_t i = 0; i < n_faces; i++)
    c_shift[i+1] += c_shift[i];

  const cs_lnum_t n_cells_add = c_shift[n_faces];

  /* Update global numbering */

  if (m->global_b_face_num != NULL || cs_glob_n_ranks > 1) {

    if (m->global_cell_num == NULL) {
      BFT_MALLOC(m->global_cell_num, n_cells_ini + n_cells_add, cs_gnum_t);
      for (cs_lnum_t i = 0; i < n_cells_ini; i++)
        m->global_cell_num[i] = i+1;
    }
    else
      BFT_REALLOC(m->global_cell_num, n_cells_ini + n_cells_add, cs_gnum_t);

    fvm_io_num_t *c_io_num
      = fvm_io_num_create_from_select(faces,
                                      m->global_b_face_num,
                                      n_faces,
                                      0);

    cs_lnum_t *n_f_sub;
    BFT_MALLOC(n_f_sub, n_faces, cs_lnum_t);
    cs_lnum_t *restrict _n_f_sub = n_f_sub;
    for (cs_lnum_t i = 0; i < n_faces; i++)
      _n_f_sub[i] = c_shift[i+1] - c_shift[i];
    _n_f_sub = NULL;

    fvm_io_num_t *c_add_io_num
      = fvm_io_num_create_from_sub(c_io_num, n_f_sub);

    c_io_num = fvm_io_num_destroy(c_io_num);

    BFT_FREE(n_f_sub);

    const cs_gnum_t n_g_cells_ini = m->n_g_cells;

    m->n_g_cells += fvm_io_num_get_global_count(c_add_io_num);

    const cs_gnum_t *c_add_gnum = fvm_io_num_get_global_num(c_add_io_num);

    assert(n_cells_add == fvm_io_num_get_local_count(c_add_io_num));

    for (cs_lnum_t i = 0; i < n_faces; i++) {
      for (cs_lnum_t j = c_shift[i]; j < c_shift[i+1]; j++)
        m->global_cell_num[n_cells_ini + j] = n_g_cells_ini + c_add_gnum[j];
    }

    c_add_io_num = fvm_io_num_destroy(c_add_io_num);

  }
  else
    m->n_g_cells += n_cells_add;

  BFT_REALLOC(m->cell_family, n_cells_ini + n_cells_add, int);

  for (cs_lnum_t i = 0; i < n_faces; i++) {
    for (cs_lnum_t j = c_shift[i]; j < c_shift[i+1]; j++)
      m->cell_family[n_cells_ini + j] = c_family[i];
  }

  /* Reset ghost cell connectivity */

  for (cs_lnum_t i = 0; i < m->n_i_faces; i++) {
    for (cs_lnum_t j = 0; j < 2; j++) {
      if (m->i_face_cells[i][j] >= m->n_cells)
        m->i_face_cells[i][j] = -1;
    }
  }

  m->n_cells += n_cells_add;
  m->n_cells_with_ghosts += n_cells_add;

  BFT_FREE(c_family);

  /* Return shift array */

  return c_shift;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create faces defining extrusion layers.
 *
 * The added boundary faces replace the initial boundary faces, while
 * the added interior faces are inserted.
 *
 * This function may only be called after \ref _add_extruded_cells, as
 * it bases the global numbering scheme for added faces on that for
 * added cells.
 *
 * \param[in, out]  m              mesh
 * \param[in]       interior_gc    if true, maintain group classes of
 *                                 interior faces previously on boundary
 * \param[in]       n_faces        number of selected boundary faces
 * \param[in]       n_vtx_ini      number of initial vertices
 * \param[in]       n_cells_ini    number of initial cells
 * \param[in]       n_g_cells_ini  global number of initial cells
 * \param[in]       faces          list of selected boundary faces (0 to n-1),
 *                                 or NULL if no indirection is needed
 * \param[in]       n_c_shift      shift for each added cell
 * \param[in]       n_v_shift      shift for each added vertex
 * \param[in]       v_s_id         id of a given mesh vertex in the selected
 *                                 vertices array
 */
/*----------------------------------------------------------------------------*/

static void
_add_layer_faces(cs_mesh_t        *m,
                 bool              interior_gc,
                 cs_lnum_t         n_faces,
                 cs_lnum_t         n_vtx_ini,
                 cs_lnum_t         n_cells_ini,
                 cs_gnum_t         n_g_cells_ini,
                 const cs_lnum_t   faces[],
                 const cs_lnum_t   n_c_shift[],
                 const cs_lnum_t   n_v_shift[],
                 const cs_lnum_t   v_s_id[])
{
  const int default_family_id = 1;

  /* Compute added size of interior faces -> vertices connectivity */

  cs_lnum_t n_add_faces = n_c_shift[n_faces];
  cs_lnum_t connect_add_size = 0;

  if (n_add_faces > 0) {

    cs_gnum_t *_g_cell_num = m->global_cell_num;

    /* Counting stage */

    for (cs_lnum_t i = 0; i < n_faces; i++) {
      cs_lnum_t f_id = faces[i];
      cs_lnum_t n_f_vtx = m->b_face_vtx_idx[f_id+1] - m->b_face_vtx_idx[f_id];
      cs_lnum_t n_f_sub = n_c_shift[i+1] - n_c_shift[i];
      connect_add_size += n_f_vtx*n_f_sub;
    }

    /* Preparation stage */

    BFT_REALLOC(m->i_face_cells, m->n_i_faces + n_add_faces, cs_lnum_2_t);

    BFT_REALLOC(m->i_face_vtx_idx, m->n_i_faces + n_add_faces + 1, cs_lnum_t);
    BFT_REALLOC(m->i_face_vtx_lst,
                m->i_face_vtx_idx[m->n_i_faces] + connect_add_size,
                cs_lnum_t);

    BFT_REALLOC(m->i_face_family, m->n_i_faces + n_add_faces, cs_lnum_t);

    bool have_g_i_face_num = (m->global_i_face_num != NULL) ? true : false;
    bool have_g_cell_num = (m->global_cell_num != NULL) ? true : false;

    /* Ensure we have all the necessary global numberings, at least as
       temporary arrays, if either global cell or interior face
       numbers are required */

    if (have_g_i_face_num || have_g_cell_num) {

      BFT_REALLOC(m->global_i_face_num, m->n_i_faces + n_add_faces, cs_gnum_t);

      if (! have_g_i_face_num) {
        for (cs_lnum_t i = 0; i < m->n_i_faces; i++)
          m->global_i_face_num[i] = i+1;
      }

      if (_g_cell_num == NULL) {
        BFT_MALLOC(_g_cell_num, m->n_cells, cs_gnum_t);
        for (cs_lnum_t i = 0; i < m->n_cells; i++)
          _g_cell_num[i] = i+1;
      }

    }

    /* Update face -> vertex connectivity */

    cs_lnum_t *vtx_idx_p = m->i_face_vtx_idx + m->n_i_faces;
    cs_lnum_t *vtx_lst_p = m->i_face_vtx_lst + m->i_face_vtx_idx[m->n_i_faces];

    for (cs_lnum_t i = 0; i < n_faces; i++) {

      cs_lnum_t f_id = faces[i];
      cs_lnum_t s_id = m->b_face_vtx_idx[f_id];
      cs_lnum_t e_id = m->b_face_vtx_idx[f_id+1];

      cs_lnum_t n_f_vtx = e_id - s_id;
      cs_lnum_t n_f_sub = n_c_shift[i+1] - n_c_shift[i];

      if (n_f_sub > 0) {

        /* First face copies previous boundary face */

        vtx_idx_p[1] = vtx_idx_p[0] + n_f_vtx;
        vtx_idx_p++;

        for (cs_lnum_t k = s_id; k < e_id; k++) {
          *vtx_lst_p = m->b_face_vtx_lst[k];
          vtx_lst_p++;
        }

        /* Other faces shifted */

        for (cs_lnum_t j = 0; j < n_f_sub-1; j++) {

          vtx_idx_p[1] = vtx_idx_p[0] + n_f_vtx;
          vtx_idx_p++;

          for (cs_lnum_t k = s_id; k < e_id; k++) {
            cs_lnum_t v_id = m->b_face_vtx_lst[k];
            cs_lnum_t l = v_s_id[v_id];
            cs_lnum_t n_v_sub = n_v_shift[l+1] - n_v_shift[l];
            if (j >= n_f_sub - n_v_sub)
              *vtx_lst_p = n_vtx_ini + n_v_shift[l] + j + n_v_sub - n_f_sub;
            else
              *vtx_lst_p = m->b_face_vtx_lst[k];
            vtx_lst_p++;
          }

        }

        /* Boundary face updated to new boundary */

        for (cs_lnum_t k = s_id; k < e_id; k++) {
          cs_lnum_t v_id = m->b_face_vtx_lst[k];
          cs_lnum_t l = v_s_id[v_id];
          cs_lnum_t n_v_sub = n_v_shift[l+1] - n_v_shift[l];
          if (n_v_sub > 0)
            m->b_face_vtx_lst[k] = n_vtx_ini + n_v_shift[l] + n_v_sub - 1;
        }

      }

    }

    m->i_face_vtx_connect_size += connect_add_size;

    /* Update face -> cells connectivity */

    cs_lnum_t n_i_faces_ini = m->n_i_faces;

    for (cs_lnum_t i = 0; i < n_faces; i++) {

      cs_lnum_t f_id = faces[i];
      cs_lnum_t s_id = n_c_shift[i];
      cs_lnum_t e_id = n_c_shift[i+1];
      cs_lnum_t n_f_sub = e_id - s_id;

      if (n_f_sub > 0) {

        /* First face links previous boundary cell to new cell */

        m->i_face_cells[n_i_faces_ini + s_id][0] = m->b_face_cells[f_id];
        m->i_face_cells[n_i_faces_ini + s_id][1] = n_cells_ini + s_id;
        if (!interior_gc)
          m->i_face_family[n_i_faces_ini + s_id] = default_family_id;
        else
          m->i_face_family[n_i_faces_ini + s_id]
            = m->b_face_family[f_id];

        /* Other faces shifted */

        for (cs_lnum_t j = 1; j < n_f_sub; j++) {
          cs_lnum_t l = n_i_faces_ini + s_id + j;
          m->i_face_cells[l][0] = n_cells_ini + s_id + j - 1;
          m->i_face_cells[l][1] = n_cells_ini + s_id + j;
          m->i_face_family[l] = default_family_id;
        }

        /* Boundary face updated to new boundary cell */

        m->b_face_cells[f_id] = n_cells_ini + e_id - 1;

      }

      if (_g_cell_num != NULL) {

        for (cs_lnum_t j = 0; j < n_f_sub; j++) {
          cs_lnum_t l = n_i_faces_ini + s_id + j;
          cs_gnum_t g_diff = _g_cell_num[n_cells_ini + s_id] - n_g_cells_ini;
          m->global_i_face_num[l] = m->n_g_i_faces + g_diff + j;
        }

      }

    }

    m->n_i_faces += n_add_faces;
    m->i_face_vtx_connect_size += connect_add_size;

    m->n_g_i_faces += (m->n_g_cells - n_g_cells_ini);

    if (_g_cell_num != m->global_cell_num)
      BFT_FREE(_g_cell_num);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create faces extruded from edges.
 *
 * This function may only be called after \ref _add_extruded_cells, as
 * it bases the global numbering scheme for added faces on that for
 * added cells.
 *
 * \param[in, out]  m              mesh
 * \param[in]       n_vtx_ini      number of initial vertices
 * \param[in]       n_cells_ini    number of initial cells
 * \param[in]       faces          list of selected boundary faces (0 to n-1),
 *                                 or NULL if no indirection is needed
 * \param[in]       e2f_stride     2 for interior edges, 1 for boundary edges
 * \param[out]      n_edges        number of edges
 * \param[out]      e2sf           edges to faces connectivity
 *                                 (size: n_edges*e2f_stride)
 * \param[out]      e2v            edge to vertices connectivity
 * \param[out]      e_gc           edge group class (or NULL)
 * \param[out]      e_io_num       edge global numbering (or NULL)
 * \param[in]       n_c_shift      shift for each added cell
 * \param[in]       n_v_shift      shift for each added vertex
 * \param[in]       v_s_id         id of a given mesh vertex in the selected
 *                                 vertices array
 */
/*----------------------------------------------------------------------------*/

static void
_add_side_faces(cs_mesh_t           *m,
                cs_lnum_t            n_vtx_ini,
                cs_lnum_t            n_cells_ini,
                cs_lnum_t            e2f_stride,
                cs_lnum_t            n_edges,
                const cs_lnum_t     *e2sf,
                const cs_lnum_2_t    e2v[],
                const int            e_gc[],
                const fvm_io_num_t  *e_io_num,
                const cs_lnum_t      n_c_shift[],
                const cs_lnum_t      n_v_shift[],
                const cs_lnum_t      v_s_id[])
{
  const int default_family_id = 1;

  /* Determine number of generated faces per edge */

  cs_lnum_t *f_shift;

  BFT_MALLOC(f_shift, n_edges+1, cs_lnum_t);

  f_shift[0] = 0;

  for (cs_lnum_t i = 0; i < n_edges; i++) {
    cs_lnum_t ns = 0;
    for (cs_lnum_t j = 0; j < 2; j++) {
      cs_lnum_t vsid = v_s_id[e2v[i][j]];
      cs_lnum_t n_v_sub = n_v_shift[vsid+1] - n_v_shift[vsid];
      if (n_v_sub > ns)
        ns = n_v_sub;
    }
    f_shift[i+1] = ns;
  }

  /* Transform count to index */

  for (cs_lnum_t i = 0; i < n_edges; i++)
    f_shift[i+1] += f_shift[i];

  const cs_lnum_t n_faces_add = f_shift[n_edges];

  /* Build associated global numbering */

  cs_gnum_t  f_ini_gcount = 0, f_add_gcount = 0;

  if (e2f_stride == 2)
    f_ini_gcount = m->n_g_i_faces;
  else if (e2f_stride == 1)
    f_ini_gcount = m->n_g_b_faces;

  cs_gnum_t *f_add_gnum = NULL;

  if (e_io_num != NULL || cs_glob_n_ranks > 1) {

    cs_lnum_t *n_f_sub;
    BFT_MALLOC(n_f_sub, n_edges, cs_lnum_t);
    cs_lnum_t *restrict _n_f_sub = n_f_sub;
    for (cs_lnum_t i = 0; i < n_edges; i++)
      _n_f_sub[i] = f_shift[i+1] - f_shift[i];
    _n_f_sub = NULL;

    fvm_io_num_t *f_add_io_num = fvm_io_num_create_from_sub(e_io_num, n_f_sub);

    f_add_gcount = fvm_io_num_get_global_count(f_add_io_num);
    f_add_gnum = fvm_io_num_transfer_global_num(f_add_io_num);

    f_add_io_num = fvm_io_num_destroy(f_add_io_num);

    BFT_FREE(n_f_sub);

  }

  /* Now generate faces */

  if (n_faces_add > 0) {

    /* Compute added face->vertices connectivity size.

       Generated faces are usualy quadrangles, but may be triangles
       if one of an edge's vertices has less extrusion layers than the
       others. The number of such faces is twice the number of extruded
       faces, minus the sum of the number of extruded vertices (the number
       of extruded faces being identical to the number of extruded
       vertices for at east one of the vertices. */

    cs_lnum_t n_tria = 0;
    cs_lnum_t n_quad = 0;

    for (cs_lnum_t i = 0; i < n_edges; i++) {
      cs_lnum_t n_f_sub = f_shift[i+1] - f_shift[i];
      if (n_f_sub < 1)
        continue;
      cs_lnum_t vid0 = e2v[i][0], vid1 = e2v[i][1];
      cs_lnum_t vsid0 = v_s_id[vid0], vsid1 = v_s_id[vid1];
      cs_lnum_t n_v_sub0 = n_v_shift[vsid0+1] - n_v_shift[vsid0];
      cs_lnum_t n_v_sub1 = n_v_shift[vsid1+1] - n_v_shift[vsid1];
      cs_lnum_t n_f_tria = CS_ABS(n_v_sub1 - n_v_sub0);
      n_tria += n_f_tria;
      n_quad += n_f_sub - n_f_tria;
    }

    /* Reallocate accordingly */

    cs_lnum_t f2v_size_add = n_quad*4 + n_tria*3;
    cs_lnum_t *a_face_cell = NULL;
    cs_lnum_t *p_face_vtx_idx = NULL;
    cs_lnum_t *p_face_vtx_lst = NULL;
    cs_lnum_t *a_face_gc = NULL;
    cs_gnum_t *a_face_gnum = NULL;

    if (e2f_stride == 2) {
      cs_lnum_t f2v_size_ini = m->i_face_vtx_idx[m->n_i_faces];
      f_ini_gcount = m->n_g_i_faces;
      BFT_REALLOC(m->i_face_cells, m->n_i_faces + n_faces_add, cs_lnum_2_t);
      BFT_REALLOC(m->i_face_vtx_idx, m->n_i_faces + n_faces_add + 1, cs_lnum_t);
      BFT_REALLOC(m->i_face_vtx_lst, f2v_size_ini + f2v_size_add, cs_lnum_t);
      BFT_REALLOC(m->i_face_family, m->n_i_faces + n_faces_add, cs_lnum_t);
      a_face_cell = (cs_lnum_t *)(m->i_face_cells + m->n_i_faces);
      p_face_vtx_idx = m->i_face_vtx_idx + m->n_i_faces;
      p_face_vtx_lst = m->i_face_vtx_lst + f2v_size_ini;
      a_face_gc = m->i_face_family + m->n_i_faces;
      if (e_io_num != NULL) {
        BFT_REALLOC(m->global_i_face_num, m->n_i_faces + n_faces_add, cs_gnum_t);
        a_face_gnum = m->global_i_face_num +  m->n_i_faces;
      }
    }
    else if (e2f_stride == 1) {
      cs_lnum_t f2v_size_ini = m->b_face_vtx_idx[m->n_b_faces];
      f_ini_gcount = m->n_g_b_faces;
      BFT_REALLOC(m->b_face_cells, m->n_b_faces + n_faces_add, cs_lnum_t);
      BFT_REALLOC(m->b_face_vtx_idx, m->n_b_faces + n_faces_add + 1, cs_lnum_t);
      BFT_REALLOC(m->b_face_vtx_lst, f2v_size_ini + f2v_size_add, cs_lnum_t);
      BFT_REALLOC(m->b_face_family, m->n_b_faces + n_faces_add, cs_lnum_t);
      a_face_cell = m->b_face_cells + m->n_b_faces;
      p_face_vtx_idx = m->b_face_vtx_idx + m->n_b_faces;
      p_face_vtx_lst = m->b_face_vtx_lst + f2v_size_ini;
      a_face_gc = m->b_face_family + m->n_b_faces;
      if (e_io_num != NULL) {
        BFT_REALLOC(m->global_b_face_num, m->n_b_faces + n_faces_add, cs_gnum_t);
        a_face_gnum = m->global_b_face_num +  m->n_b_faces;
      }
    }

    /* Now generate new faces */

    for (cs_lnum_t i = 0; i < n_edges; i++) {

      cs_lnum_t n_f_sub = f_shift[i+1] - f_shift[i];

      if (n_f_sub < 1)
        continue;

      cs_lnum_t vid0 = e2v[i][0], vid1 = e2v[i][1];
      cs_lnum_t vsid0 = v_s_id[vid0], vsid1 = v_s_id[vid1];

      cs_lnum_t n_v_sub0 = n_v_shift[vsid0+1] - n_v_shift[vsid0];
      cs_lnum_t n_v_sub1 = n_v_shift[vsid1+1] - n_v_shift[vsid1];

      cs_lnum_t n_v_diff = CS_ABS(n_v_sub1 - n_v_sub0);

      /* Face -> vertices connectivity;
         when edge vertices do not have the same extrusion count,
         arrange for triangular faces first (closest to the interior
         of the mesh), then regular quadrangle faces last (near the
         boundary). */

      cs_lnum_t l0 = -1, l1 = -1;

      /* Triangles */

      for (cs_lnum_t j = 0; j < n_v_diff; j++) {

        p_face_vtx_idx[1] = p_face_vtx_idx[0] + 3;

        if (n_v_sub0 < n_v_sub1) {
          p_face_vtx_lst[0] = vid0;
          if (l1 < 0)
            p_face_vtx_lst[1] = vid1;
          else
            p_face_vtx_lst[1] = n_vtx_ini + n_v_shift[vsid1] + l1;
          p_face_vtx_lst[2] = n_vtx_ini + n_v_shift[vsid1] + l1 + 1;
          l1++;
        }
        else { /* n_v_sub0 > n_v_sub1 */
          if (l0 < 0)
            p_face_vtx_lst[0] = vid0;
          else
            p_face_vtx_lst[0] = n_vtx_ini + n_v_shift[vsid0] + l0;
          p_face_vtx_lst[1] = vid1;
          p_face_vtx_lst[2] = n_vtx_ini + n_v_shift[vsid0] + l0 + 1;
          l0++;
        }

        p_face_vtx_idx += 1;
        p_face_vtx_lst += 3;

      }

      /* Quadrangles */

      for (cs_lnum_t j = 0; j < n_f_sub; j++) {

        if (j < n_v_diff)
          continue;

        p_face_vtx_idx[1] = p_face_vtx_idx[0] + 4;

        if (l0 < 0)
          p_face_vtx_lst[0] = vid0;
        else
          p_face_vtx_lst[0] = n_vtx_ini + n_v_shift[vsid0] + l0;
        if (l1 < 0)
          p_face_vtx_lst[1] = vid1;
        else
          p_face_vtx_lst[1] = n_vtx_ini + n_v_shift[vsid1] + l1;
        p_face_vtx_lst[2] = n_vtx_ini + n_v_shift[vsid1] + l1 + 1;
        p_face_vtx_lst[3] = n_vtx_ini + n_v_shift[vsid0] + l0 + 1;
        l0++;
        l1++;

        p_face_vtx_idx += 1;
        p_face_vtx_lst += 4;

      }

      /* Face -> cells connectivity;
         when faces have less extrusion layers than adjacent cells
         (due to some other faces of the cell having more layers),
         we arrange for the face to be adjacent to the outermost
         layers (so as to have regular prism cells near to the boundary). */

      for (cs_lnum_t j = 0; j < e2f_stride; j++) {

        cs_lnum_t f_id = e2sf[i*e2f_stride + j];

        if (f_id > -1) {
          cs_lnum_t n_c_sub = n_c_shift[f_id+1] - n_c_shift[f_id];
          for (cs_lnum_t k = 0, l = n_c_sub - n_f_sub;
               k < n_f_sub; k++, l++) {
            cs_lnum_t nf = f_shift[i] + k;
            a_face_cell[e2f_stride*nf + j] = n_cells_ini + n_c_shift[f_id] + l;
          }
        }
        else {
          for (cs_lnum_t k = 0; k < n_f_sub; k++) {
            cs_lnum_t nf = f_shift[i] + k;
            a_face_cell[e2f_stride*nf + j] = -1;
          }
        }

      }

      /* Face family */

      if (e2f_stride != 1) {
        for (cs_lnum_t k = 0; k < n_f_sub; k++)
          a_face_gc[f_shift[i] + k] = default_family_id;
      }
      else {
        for (cs_lnum_t k = 0; k < n_f_sub; k++)
          a_face_gc[f_shift[i] + k] = e_gc[i];
      }

      /* Global number */

      if (a_face_gnum != NULL) {
        for (cs_lnum_t k = 0; k < n_f_sub; k++)
          a_face_gnum[f_shift[i] + k]
            = f_ini_gcount + f_add_gnum[f_shift[i] + k];
      }
    }
  }

  /* Update metadata */

  if (e2f_stride == 2) {
    m->n_i_faces += n_faces_add;
    m->i_face_vtx_connect_size = m->i_face_vtx_idx[m->n_i_faces];
    m->n_g_i_faces += f_add_gcount;
  }
  else if (e2f_stride == 1) {
    m->n_b_faces += n_faces_add;
    m->b_face_vtx_connect_size = m->b_face_vtx_idx[m->n_b_faces];
    m->n_g_b_faces += f_add_gcount;
  }

  /* Free temporary arrays */

  BFT_FREE(f_add_gnum);
  BFT_FREE(f_shift);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define extrusion vectors by face info.
 *
 * \param[in, out]  e    extrusion vectors structure
 * \param[in]       efi  extrusion face info structure
 */
/*----------------------------------------------------------------------------*/

static void
_cs_mesh_extrude_vectors_by_face_info(cs_mesh_extrude_vectors_t          *e,
                                      const cs_mesh_extrude_face_info_t  *efi)
{
  cs_mesh_t *m = cs_glob_mesh;

  cs_mesh_quantities_t *mq = cs_mesh_quantities_create();

  cs_mesh_quantities_compute_preprocess(m, mq);

  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_vertices = m->n_vertices;

  cs_lnum_t *_n_layers = NULL;
  float *_expansion = NULL;
  cs_coord_3_t *_coord_shift = NULL;
  cs_real_2_t *_thickness_se = NULL;

  /* Determine vertices to extrude */

  cs_real_2_t *w = NULL;
  int *c = NULL;

  BFT_MALLOC(_n_layers, n_vertices, cs_lnum_t);
  BFT_MALLOC(_expansion, n_vertices, float);
  BFT_MALLOC(_thickness_se, n_vertices, cs_real_2_t);
  BFT_MALLOC(_coord_shift, n_vertices, cs_coord_3_t);

  BFT_MALLOC(w, n_vertices, cs_real_2_t);
  BFT_MALLOC(c, n_vertices, int);

  for (cs_lnum_t i = 0; i < n_vertices; i++) {
    _n_layers[i] = 0;
    _expansion[i] = 0;
    _thickness_se[i][0] = 0;
    _thickness_se[i][1] = 0;
    _coord_shift[i][0] = 0;
    _coord_shift[i][1] = 0;
    _coord_shift[i][2] = 0;
    w[i][0] = 0;
    w[i][1] = 0;
    c[i] = 0;
  }

  /* Global check for zone thickness specification; compute default
     thickness based on smoothed cell size it thickness not defined. */

  int z_thickness_spec = 1;

  for (cs_lnum_t i = 0; i < n_b_faces; i++) {
    if (efi->n_layers[i] != 0 && efi->distance[i] < 0)
      z_thickness_spec = 0;
  }

  cs_parall_min(1, CS_INT_TYPE, &z_thickness_spec);

  cs_real_t *_distance;
  BFT_MALLOC(_distance, n_b_faces, cs_real_t);

  if (z_thickness_spec < 1)
    cs_mesh_quantities_b_thickness_f(m,
                                     mq,
                                     3, /* n_passes */
                                     _distance);

  for (cs_lnum_t i = 0; i < n_b_faces; i++) {
    if (efi->n_layers[i] == 0)
      _distance[i] = 0;
    else if (efi->distance[i] < 0)
      _distance[i] *= -efi->distance[i];
    else
      _distance[i] = efi->distance[i];
  }

  /* Build selected faces ids array */

  {
    cs_lnum_t n_faces = 0;

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      if (efi->n_layers[f_id] > -1)
        n_faces++;
    }

    BFT_REALLOC(e->face_ids, n_faces, cs_lnum_t);

    e->n_faces = n_faces;

    n_faces = 0;

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      if (efi->n_layers[f_id] > -1)
        e->face_ids[n_faces++] = f_id;
    }

    assert(n_faces == e->n_faces);
  }

  /* Now determine other parameters */

  for (cs_lnum_t j = 0; j < e->n_faces; j++) {

    cs_lnum_t f_id = e->face_ids[j];

    const cs_lnum_t n_layers = efi->n_layers[f_id];
    const cs_real_t distance = _distance[f_id];
    const cs_real_t expansion_factor = efi->expansion_factor[f_id];
    const cs_real_t thickness_s = n_layers > 2 ? efi->thickness_s[f_id] : 0;
    const cs_real_t thickness_e = n_layers > 1 ? efi->thickness_e[f_id] : 0;

    cs_lnum_t s_id = m->b_face_vtx_idx[f_id];
    cs_lnum_t e_id = m->b_face_vtx_idx[f_id+1];
    const cs_real_t *f_n = mq->b_face_normal + f_id*3;
    const cs_real_t f_s = cs_math_3_norm(f_n);

    for (cs_lnum_t k = s_id; k < e_id; k++) {
      cs_lnum_t v_id = m->b_face_vtx_lst[k];
      cs_lnum_t k_0 = (k < e_id-1) ? k+1 : s_id;
      cs_lnum_t k_1 = (k > s_id) ? k-1 : e_id-1;
      cs_lnum_t v_ids[3] = {v_id,
                            m->b_face_vtx_lst[k_0],
                            m->b_face_vtx_lst[k_1]};
      cs_real_t a = _angle(v_ids, f_n, (const cs_real_3_t *)(m->vtx_coord));
      _n_layers[v_id] += n_layers;
      _expansion[v_id] += expansion_factor * a;
      _thickness_se[v_id][0] += thickness_s * a;
      _thickness_se[v_id][1] += thickness_e * a;
      for (cs_lnum_t l = 0; l < 3; l++)
        _coord_shift[v_id][l] += a * f_n[l]/f_s;
      w[v_id][0] += a;
      w[v_id][1] += distance * a;
      c[v_id] += 1;
    }
  }

  BFT_FREE(_distance);

  /* Handle parallelism */

  if (m->vtx_interfaces != NULL) {
    cs_interface_set_sum(m->vtx_interfaces,
                         m->n_vertices,
                         1,
                         true,
                         CS_LNUM_TYPE,
                         _n_layers);
    cs_interface_set_sum(m->vtx_interfaces,
                         m->n_vertices,
                         1,
                         true,
                         CS_FLOAT,
                         _expansion);
    cs_interface_set_sum(m->vtx_interfaces,
                         m->n_vertices,
                         2,
                         true,
                         CS_REAL_TYPE,
                         _thickness_se);
    cs_interface_set_sum(m->vtx_interfaces,
                         m->n_vertices,
                         3,
                         true,
                         CS_COORD_TYPE,
                         _coord_shift);
    cs_interface_set_sum(m->vtx_interfaces,
                         m->n_vertices,
                         2,
                         true,
                         CS_REAL_TYPE,
                         w);
    cs_interface_set_sum(m->vtx_interfaces,
                         m->n_vertices,
                         1,
                         true,
                         CS_INT_TYPE,
                         c);
  }

  for (cs_lnum_t i = 0; i < m->n_vertices; i++) {
    if (c[i] > 0) {
      _n_layers[i] /= c[i];
      _expansion[i] /= w[i][0];
      _thickness_se[i][0] /= w[i][0];
      _thickness_se[i][1] /= w[i][0];
      cs_real_t nn = cs_math_3_square_norm(_coord_shift[i]);
      if (nn > FLT_MIN) {
        for (cs_lnum_t l = 0; l < 3; l++)
          _coord_shift[i][l] *= (w[i][1] / nn);
      }
      else {
        _n_layers[i] = 0;
      }
    }
    else
      _n_layers[i] = 0;
  }

  /* Check for opposing normal directions (may occur at edges of thin
     boundaries; block extrusion there) */

  for (cs_lnum_t j = 0; j < e->n_faces; j++) {

    cs_lnum_t f_id = e->face_ids[j];

    cs_lnum_t s_id = m->b_face_vtx_idx[f_id];
    cs_lnum_t e_id = m->b_face_vtx_idx[f_id+1];

    const cs_real_t *f_n = mq->b_face_normal + f_id*3;

    for (cs_lnum_t k = s_id; k < e_id; k++) {
      cs_lnum_t v_id = m->b_face_vtx_lst[k];
      if (cs_math_3_dot_product(_coord_shift[v_id], f_n) < 0) {
        _n_layers[v_id] = 0;
        _coord_shift[v_id][0] = 0;
        _coord_shift[v_id][1] = 0;
        _coord_shift[v_id][2] = 0;
      }
    }
  }

  if (m->vtx_interfaces != NULL)
    cs_interface_set_min(m->vtx_interfaces,
                         m->n_vertices,
                         1,
                         true,
                         CS_LNUM_TYPE,
                         _n_layers);

  /* Limit jumps in number of layers */

  cs_gnum_t layer_limiter = 0;

  do {

    layer_limiter = 0;

    for (cs_lnum_t j = 0; j < e->n_faces; j++) {

      cs_lnum_t f_id = e->face_ids[j];

      cs_lnum_t s_id = m->b_face_vtx_idx[f_id];
      cs_lnum_t e_id = m->b_face_vtx_idx[f_id+1];

      cs_lnum_t n_l_min = _n_layers[m->b_face_vtx_lst[s_id]];
      int up = 0;
      int down = 0;
      for (cs_lnum_t k = s_id+1; k < e_id; k++) {
        cs_lnum_t l = (k+1) < e_id ? k+1 : s_id;
        cs_lnum_t v_id0 = m->b_face_vtx_lst[k];
        cs_lnum_t v_id1 = m->b_face_vtx_lst[l];
        if (_n_layers[v_id1] > _n_layers[v_id0])
          up += 1;
        else if (_n_layers[v_id1] < _n_layers[v_id0])
          down += 1;
        if (_n_layers[v_id0] < n_l_min)
          n_l_min = _n_layers[v_id0];
      }
      if (up > 1 || down > 1) {
        layer_limiter += 1;
        for (cs_lnum_t k = s_id; k < e_id; k++) {
          cs_lnum_t v_id = m->b_face_vtx_lst[k];
          _n_layers[v_id] = n_l_min;
        }
      }

    }

    cs_parall_counter(&layer_limiter, 1);

    if (layer_limiter > 0 && m->vtx_interfaces != NULL)
      cs_interface_set_min(m->vtx_interfaces,
                           m->n_vertices,
                           1,
                           true,
                           CS_LNUM_TYPE,
                           _n_layers);

  } while (layer_limiter > 0);

  /* Now determine other parameters */

  /* Free temporaries */

  BFT_FREE(c);
  BFT_FREE(w);

  mq = cs_mesh_quantities_destroy(mq);

  /* Build vertex selection list */

  _select_vertices_from_adj_b_faces(m,
                                    e->n_faces,
                                    e->face_ids,
                                    &(e->n_vertices),
                                    &(e->vertex_ids));

  BFT_REALLOC(e->n_layers, e->n_vertices, cs_lnum_t);
  BFT_REALLOC(e->coord_shift, e->n_vertices, cs_coord_3_t);

  for (cs_lnum_t i = 0; i < e->n_vertices; i++) {
    cs_lnum_t v_id = e->vertex_ids[i];
    e->n_layers[i] = _n_layers[v_id];
    for (cs_lnum_t l = 0; l < 3; l++)
      e->coord_shift[i][l] = _coord_shift[v_id][l];
  }

  BFT_FREE(_n_layers);
  BFT_FREE(_coord_shift);

  BFT_REALLOC(e->distribution_idx, e->n_vertices+1, cs_lnum_t);

  e->distribution_idx[0] = 0;
  for (cs_lnum_t i = 0; i < e->n_vertices; i++)
    e->distribution_idx[i+1] = e->distribution_idx[i] + e->n_layers[i];

  BFT_REALLOC(e->distribution, e->distribution_idx[e->n_vertices], float);

  /* Compute distribution for each extruded vertex */

  for (cs_lnum_t i = 0; i < e->n_vertices; i++) {

    int n_l = e->n_layers[i];

    if (n_l > 0) {

      cs_lnum_t v_id = e->vertex_ids[i];
      float *_d = e->distribution + e->distribution_idx[i];

      cs_lnum_t s_id = 0;
      cs_lnum_t e_id = n_l;
      double d_expansion = 1;

      /* Handle optional start and end thickness:
         - transform to relative values, ensure they are within bounds
           of extrusion distance (reducing them if necessary)
         - adjust start and end layers for expansion factor-based layers */

      for (int j = 0; j < 2; j++) {
        if (_thickness_se[v_id][j] < 0)
          _thickness_se[v_id][j] = 0;
      }

      if (_thickness_se[v_id][0] + _thickness_se[v_id][1] > 0) {
        cs_real_t vn = cs_math_3_norm(e->coord_shift[i]);
        if (vn > 0) {
          _thickness_se[v_id][0] /= vn;
          _thickness_se[v_id][1] /= vn;
        }
        cs_real_t f = 0;
        if (_thickness_se[v_id][1] > 0) {
          f += _thickness_se[v_id][1];
          n_l -= 1;
          e_id -= 1;
        }
        if (n_l > 0 && _thickness_se[v_id][1] > 0) {
          f += _thickness_se[v_id][0];
          n_l -= 1;
          s_id += 1;
        }
        if (f > 1) {
          if (n_l > 0)
            f *= 1.2;
          _thickness_se[v_id][0] /= f;
          _thickness_se[v_id][1] /= f;
        }

        d_expansion -= (_thickness_se[v_id][0] + _thickness_se[v_id][1]);
      }

      /* Estimation to reference distance */

      _d[s_id] = 1;
      for (cs_lnum_t l_id = s_id + 1; l_id < e_id; l_id++)
        _d[l_id] = _d[l_id-1]*_expansion[v_id];
      double d_tot = 0;
      for (cs_lnum_t l_id = s_id; l_id < e_id; l_id++)
        d_tot += _d[l_id];

      /* Now set distances */

      assert(s_id == 0);
      if (s_id > 0)
        _d[0] = _thickness_se[v_id][0];

      _d[s_id] = d_expansion/d_tot;
      for (cs_lnum_t l_id = s_id+1; l_id < e_id-1; l_id++)
        _d[l_id] = _d[l_id-1] + _d[l_id]/d_tot;

      if (e_id < n_l)
        _d[e_id-1] = _thickness_se[v_id][1];

      _d[n_l-1] = 1.0;
    }

    else { /* n_l = 0 */
      for (cs_lnum_t l = 0; l < 3; l++)
        e->coord_shift[i][l] = 0.;
    }

  }

  BFT_FREE(_expansion);
  BFT_FREE(_thickness_se);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Extrude mesh boundary faces in the normal direction.
 *
 * Extrusion is defined on selected boundary faces, and the number of layers
 * for each associated vertex may be (slightly) variable, to account for
 * cluttered areas where extrusion may be constrained, or more complex
 * extrusions.
 *
 * \param[in, out]  m             mesh
 * \param[in]       e             extrusion vector definitions
 * \param[in]       interior_gc   if true, maintain group classes of
 *                                interior faces previously on boundary
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_extrude(cs_mesh_t                        *m,
                const cs_mesh_extrude_vectors_t  *e,
                bool                              interior_gc)
{
  cs_lnum_t *l_faces = NULL, *l_vertices = NULL;

  const cs_lnum_t *_faces = e->face_ids;
  const cs_lnum_t *_vertices = e->vertex_ids;

  cs_lnum_t  n_cells_ini = m->n_cells;
  cs_lnum_t  n_vtx_ini = m->n_vertices;
  cs_gnum_t  n_g_cells_ini = m->n_g_cells;

  /* Check we have something to do */

  cs_gnum_t n_g_sel_faces = e->n_faces;
  cs_parall_counter(&n_g_sel_faces, 1);

  bft_printf(_("\n"
               " Extrusion: %llu boundary faces selected.\n"),
             (unsigned long long)n_g_sel_faces);

  if (n_g_sel_faces < 1)
    return;

  cs_mesh_free_rebuildable(m, false);

  /* Local names for parameters */

  const cs_lnum_t      n_faces = e->n_faces;
  const cs_lnum_t      n_vertices = e->n_vertices;
  const cs_lnum_t     *faces = e->face_ids;
  const cs_lnum_t     *vertices = e->vertex_ids;
  const cs_lnum_t     *n_layers_ini = e->n_layers;
  const cs_coord_3_t  *coord_shift = (const cs_coord_3_t *)e->coord_shift;
  const float         *distribution = e->distribution;

  /* Ensure we have explicit selections */

  if (faces == NULL) {
    BFT_MALLOC(l_faces, n_faces, cs_lnum_t);
    for (cs_lnum_t i = 0; i < n_faces; i++)
      l_faces[i] = i;
    _faces = l_faces;
  }

  if (vertices == NULL) {
    BFT_MALLOC(l_vertices, n_vertices, cs_lnum_t);
    for (cs_lnum_t i = 0; i < n_vertices; i++)
      l_vertices[i] = i;
    _vertices = l_vertices;
  }

  /* Build local edges */

  cs_lnum_t n_i_edges = 0, n_b_edges = 0;
  cs_lnum_2_t  *i_e2sf = NULL, *i_e2v = NULL, *b_e2v = NULL;
  cs_lnum_t  *b_e2sf = NULL;
  int  *b_e_gc = NULL;
  fvm_io_num_t  *i_e_io_num = NULL, *b_e_io_num = NULL;

  cs_lnum_t  *_n_layers
    = _build_face_edges(m,
                        n_faces,
                        n_vertices,
                        _faces,
                        _vertices,
                        n_layers_ini,
                        &n_i_edges,
                        &n_b_edges,
                        &i_e2sf,
                        &i_e2v,
                        &i_e_io_num,
                        &b_e2sf,
                        &b_e2v,
                        &b_e_gc,
                        &b_e_io_num);

  const cs_lnum_t *n_layers = (const cs_lnum_t *)_n_layers;
  if (_n_layers == NULL)
    n_layers = n_layers_ini;

  /* Now generate added vertices */

  cs_lnum_t *n_v_shift = _add_extruded_vertices(m,
                                                n_vertices,
                                                _vertices,
                                                n_layers_ini,
                                                n_layers,
                                                coord_shift,
                                                distribution);

  /* Now determine number of extruded faces per edge and cells per face */

  cs_lnum_t *n_v_sub;
  BFT_MALLOC(n_v_sub, m->n_vertices, cs_lnum_t);

  for (cs_lnum_t v_id = 0; v_id < m->n_vertices; v_id++)
    n_v_sub[v_id] = 0;

  for (cs_lnum_t i = 0; i < n_vertices; i++)
    n_v_sub[_vertices[i]] = n_layers[i];

  cs_lnum_t *n_c_shift = _add_extruded_cells(m,
                                             n_faces,
                                             _faces,
                                             n_v_sub);

  BFT_FREE(n_v_sub);

  bft_printf(_("            %llu cells added.\n"),
             (unsigned long long)(m->n_g_cells - n_g_cells_ini));

  /* Add faces whose normals are parallel to the extrusion direction */

  cs_lnum_t *v_s_id;
  BFT_MALLOC(v_s_id, m->n_vertices, cs_lnum_t);

  for (cs_lnum_t v_id = 0; v_id < m->n_vertices; v_id++)
    v_s_id[v_id] = -1;

  for (cs_lnum_t i = 0; i < n_vertices; i++)
    v_s_id[_vertices[i]] = i;

  _add_layer_faces(m,
                   interior_gc,
                   n_faces,
                   n_vtx_ini,
                   n_cells_ini,
                   n_g_cells_ini,
                   _faces,
                   n_c_shift,
                   n_v_shift,
                   v_s_id);

  /* Extrude edges for interior and boundary faces */

  _add_side_faces(m,
                  n_vtx_ini,
                  n_cells_ini,
                  2, /* e2f_stride */
                  n_i_edges,
                  (cs_lnum_t *)i_e2sf,
                  (const cs_lnum_2_t *)i_e2v,
                  NULL, /* e_gc */
                  i_e_io_num,
                  n_c_shift,
                  n_v_shift,
                  v_s_id);

  _add_side_faces(m,
                  n_vtx_ini,
                  n_cells_ini,
                  1, /* e2f_stride */
                  n_b_edges,
                  b_e2sf,
                  (const cs_lnum_2_t *)b_e2v,
                  b_e_gc,
                  b_e_io_num,
                  n_c_shift,
                  n_v_shift,
                  v_s_id);

  /* Free numberings */

  if (i_e_io_num != NULL)
    i_e_io_num = fvm_io_num_destroy(i_e_io_num);

  if (b_e_io_num != NULL)
    b_e_io_num = fvm_io_num_destroy(b_e_io_num);

  /* Free local arrays */

  n_layers = NULL;
  BFT_FREE(_n_layers);

  BFT_FREE(v_s_id);

  BFT_FREE(n_c_shift);
  BFT_FREE(n_v_shift);

  BFT_FREE(i_e2sf);
  BFT_FREE(b_e2sf);
  BFT_FREE(i_e2v);
  BFT_FREE(b_e2v);
  BFT_FREE(b_e_gc);

  BFT_FREE(l_faces);
  BFT_FREE(l_vertices);

  m->modified = CS_MAX(m->modified, 1);

  /* Rebuild ghosts */

  if (m->halo != NULL || m->halo_type == CS_HALO_EXTENDED) {
    cs_halo_type_t halo_type = m->halo_type;
    cs_halo_destroy(&(m->halo));
    assert(m == cs_glob_mesh);
    cs_mesh_builder_t *mb = (m == cs_glob_mesh) ? cs_glob_mesh_builder : NULL;
    cs_mesh_init_halo(m, mb, halo_type);
  }

  cs_mesh_update_auxiliary(cs_glob_mesh);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Extrude mesh boundary faces in the normal direction by a constant
 *        thickness.
 *
 * \param[in, out]  m                 mesh
 * \param[in]       interior_gc       if true, maintain group classes of
 *                                    interior faces previously on boundary
 * \param[in]       n_layers          number of layers
 * \param[in]       thickness         extrusion thickness
 * \param[in]       expansion_factor  geometric expansion factor for
 *                                    extrusion refinement
 * \param[in]       n_faces           number of selected boundary faces
 * \param[in]       faces             list of selected boundary faces (0 to n-1),
 *                                    or NULL if no indirection is needed
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_extrude_constant(cs_mesh_t        *m,
                         bool              interior_gc,
                         cs_lnum_t         n_layers,
                         double            thickness,
                         double            expansion_factor,
                         cs_lnum_t         n_faces,
                         const cs_lnum_t   faces[])
{
  cs_mesh_extrude_face_info_t *efi = cs_mesh_extrude_face_info_create(m);

  cs_mesh_extrude_set_info_by_zone(efi,
                                   n_layers,
                                   thickness,
                                   expansion_factor,
                                   n_faces,
                                   faces);

  cs_mesh_extrude_vectors_t *e= cs_mesh_extrude_vectors_create(efi);

  cs_mesh_extrude_face_info_destroy(&efi);

  /* Build selection list */

  /* Now call lower-level function */

  cs_mesh_extrude(m, e, interior_gc);

  /* Free work arrays */

  cs_mesh_extrude_vectors_destroy(&e);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a mesh extrusion face information structure.
 *
 * \param[in]  m  mesh
 *
 * \return pointer to new mesh extrusion face information structure.
 */
/*----------------------------------------------------------------------------*/

cs_mesh_extrude_face_info_t *
cs_mesh_extrude_face_info_create(const cs_mesh_t  *m)
{
  cs_mesh_extrude_face_info_t *efi;

  const cs_lnum_t n_faces = m->n_b_faces;

  BFT_MALLOC(efi, 1, cs_mesh_extrude_face_info_t);

  BFT_MALLOC(efi->n_layers, n_faces, cs_lnum_t);
  BFT_MALLOC(efi->distance, n_faces, cs_real_t);
  BFT_MALLOC(efi->expansion_factor, n_faces, float_t);
  BFT_MALLOC(efi->thickness_s, n_faces, cs_real_t);
  BFT_MALLOC(efi->thickness_e, n_faces, cs_real_t);

  for (cs_lnum_t i = 0; i < n_faces; i++) {
    efi->n_layers[i] = -1;
    efi->distance[i] = -1;
    efi->expansion_factor[i] = 0.8;
    efi->thickness_s[i] = 0;
    efi->thickness_e[i] = 0;
  }

  return efi;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a mesh extrusion face information structure.
 *
 * \param[in, out]  efi  pointer to pointer to mesh extrusion face information.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_extrude_face_info_destroy(cs_mesh_extrude_face_info_t  **efi)
{
  if (efi != NULL) {
    cs_mesh_extrude_face_info_t *_efi = *efi;
    if (_efi != NULL) {
      BFT_FREE(_efi->n_layers);
      BFT_FREE(_efi->distance);
      BFT_FREE(_efi->expansion_factor);
      BFT_FREE(_efi->thickness_s);
      BFT_FREE(_efi->thickness_e);
      BFT_FREE(*efi);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set face extrusion information by zone.
 *
 * \param[in, out]  efi               mesh extrusion face information
 * \param[in]       n_layers          number of layers for selected faces
 * \param[in]       distance          extrusion distance for selected faces
 *                                    (if < 0, absolute value used as
 *                                    multiplier for boundary cell thickness)
 * \param[in]       expansion_factor  expansion factor for selected faces
 * \param[in]       n_faces           number of selected faces
 * \param[in]       face_ids          ids of selected faces, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_extrude_set_info_by_zone(cs_mesh_extrude_face_info_t  *efi,
                                 int                           n_layers,
                                 double                        distance,
                                 float                         expansion_factor,
                                 const cs_lnum_t               n_faces,
                                 const cs_lnum_t               face_ids[])
{
  if (efi == NULL)
    return;

  if (face_ids != NULL) {
    for (cs_lnum_t i = 0; i < n_faces; i++) {
      cs_lnum_t f_id = face_ids[i];
      efi->n_layers[f_id] = n_layers;
      efi->distance[f_id] = distance;
      efi->expansion_factor[f_id] = expansion_factor;
      efi->thickness_s[f_id] = 0;
      efi->thickness_e[f_id] = 0;
    }
  }

  else {
    for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++) {
      efi->n_layers[f_id] = n_layers;
      efi->distance[f_id] = distance;
      efi->expansion_factor[f_id] = expansion_factor;
      efi->thickness_s[f_id] = 0;
      efi->thickness_e[f_id] = 0;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and build a mesh extrusion vectors definition.
 *
 * Extrusion vectors will be computed based on the provided extrusion
 * face information structure. If no such structure is provided, an empty
 * structure is returned.
 *
 * \param[in]  efi  mesh extrusion face information, or NULL
 *
 * \return pointer to created mesh extrusion vectors definition.
 */
/*----------------------------------------------------------------------------*/

cs_mesh_extrude_vectors_t *
cs_mesh_extrude_vectors_create(const cs_mesh_extrude_face_info_t  *efi)
{
  cs_mesh_extrude_vectors_t *e;

  BFT_MALLOC(e, 1, cs_mesh_extrude_vectors_t);

  e->n_faces = 0;
  e->n_vertices = 0;
  e->face_ids = NULL;
  e->vertex_ids = NULL;
  e->n_layers = NULL;
  e->coord_shift = NULL;
  e->distribution_idx = NULL;
  e->distribution = NULL;

  if (efi != NULL)
    _cs_mesh_extrude_vectors_by_face_info(e, efi);

  return e;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a mesh extrusion vectors definition.
 *
 *
 * \param[in, out]  e  pointer to pointer to mesh extrusion vectors definition.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_extrude_vectors_destroy(cs_mesh_extrude_vectors_t **e)
{
  if (e != NULL) {
    cs_mesh_extrude_vectors_t *_e = *e;
    if (_e != NULL) {
      BFT_FREE(_e->face_ids);
      BFT_FREE(_e->vertex_ids);
      BFT_FREE(_e->n_layers);
      BFT_FREE(_e->coord_shift);
      BFT_FREE(_e->distribution_idx);
      BFT_FREE(_e->distribution);
      BFT_FREE(*e);
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
