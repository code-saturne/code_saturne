/*============================================================================
 * Manage connectivity (Topological features of the mesh)
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_order.h"
#include "cs_sort.h"
#include "cs_cdo.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_connect.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local macro and structure definitions
 *============================================================================*/

#define INNOV_CONNECT_DEBUG 0

/* Temporary structure to build mesh connectivities */
typedef struct {

  cs_lnum_t  n_vertices;
  cs_lnum_t  n_edges;
  cs_lnum_t  n_faces;  /* border + interior faces */
  cs_lnum_t  n_cells;

  int  *e2v_lst;  /* Edge ref. definition (2*n_edges) */

  /* Vertex -> edge connect based on vertex -> vertices connect. */

  int  *v2v_idx;
  int  *v2v_lst;
  int  *v2v_edge_lst;

  /* Cell -> faces connect. (int + border faces) */

  int  *c2f_idx;
  int  *c2f_lst;

} cs_cdo_connect_builder_t;

/*============================================================================
 * Private function prototypes
 *============================================================================*/


static int
_get_edge_num(int  v1_num,
              int  v2_num,
              const cs_cdo_connect_builder_t  *cbuilder)
{
  int  i;
  int  edge_num = 0;

  assert(cbuilder != NULL);
  assert(cbuilder->v2v_idx != NULL);
  assert(v1_num > 0);
  assert(v2_num > 0);

  if (cbuilder->v2v_idx[v1_num] - cbuilder->v2v_idx[v1_num-1] == 0)
    bft_error(__FILE__, __LINE__, 0,
              _(" The given vertex number: %d is not defined"
                " in the edge structure (edges->vtx_idx).\n"), v1_num);

  for (i = cbuilder->v2v_idx[v1_num-1]; i < cbuilder->v2v_idx[v1_num]; i++) {
    if (cbuilder->v2v_lst[i] == v2_num) {
      edge_num = cbuilder->v2v_edge_lst[i];
      break;
    }
  }

  if (edge_num == 0)
    bft_error(__FILE__, __LINE__, 0,
              _(" The given couple of vertex numbers :\n"
                "   vertex 1 : %d\n"
                "   vertex 2 : %d\n"
                " is not defined in the edge structure.\n"), v1_num, v2_num);

  return edge_num;
}

/*----------------------------------------------------------------------------*/

static void
_create_c2f_connect(const cs_mesh_t     *mesh,
                    cs_cdo_connect_builder_t  *cbuilder)
{
  int  i, cid1, cid2, shift;

  int  *cell_shift = NULL, *c2f_idx = NULL, *c2f_lst = NULL;

  const int  n_i_faces = mesh->n_i_faces;
  const int  n_b_faces = mesh->n_b_faces;
  const int  n_cells = mesh->n_cells;

  assert(c2f_idx == NULL && c2f_lst == NULL);

  BFT_MALLOC(c2f_idx, n_cells + 1, int);
  BFT_MALLOC(cell_shift, n_cells, int);

  c2f_idx[0] = 0;
  for (i = 0; i < n_cells; i++) {
    c2f_idx[i+1] = 0;
    cell_shift[i] = 0;
  }

  /* Count */
  for (i = 0; i < n_i_faces; i++) {

    cid1 = mesh->i_face_cells[i][0];
    cid2 = mesh->i_face_cells[i][1];

    if (cid1 < n_cells)
      c2f_idx[cid1+1] += 1;

    if (cid2 < n_cells)
      c2f_idx[cid2+1] += 1;

  } /* End of loop on internal faces */

  for (i = 0; i < n_b_faces; i++) {

    cid1 = mesh->b_face_cells[i];
    assert(cid1 < n_cells);
    c2f_idx[cid1+1] += 1;

  } /* End of loop on border faces */

  /* Build index */
  for (i = 0; i < n_cells; i++)
    c2f_idx[i+1] += c2f_idx[i];

  BFT_MALLOC(c2f_lst, c2f_idx[n_cells], int);

  for (i = 0; i < n_i_faces; i++) {

    cid1 = mesh->i_face_cells[i][0];
    cid2 = mesh->i_face_cells[i][1];

    if (cid1 < n_cells) {
      shift = c2f_idx[cid1] + cell_shift[cid1];
      c2f_lst[shift] = i+1;
      cell_shift[cid1] += 1;
    }

    if (cid2 < n_cells) {
      shift = c2f_idx[cid2] + cell_shift[cid2];
      c2f_lst[shift] = i+1;
      cell_shift[cid2] += 1;
    }

  } /* End of loop on internal faces */

  for (i = 0; i < n_b_faces; i++) {

    cid1 = mesh->b_face_cells[i];
    shift = c2f_idx[cid1] + cell_shift[cid1];
    c2f_lst[shift] = n_i_faces+i+1;
    cell_shift[cid1] += 1;

  } /* End of loop on border faces */

  /* Free memory */
  BFT_FREE(cell_shift);

  cbuilder->c2f_idx = c2f_idx;
  cbuilder->c2f_lst = c2f_lst;

}

/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_build_c2f_connect(const cs_mesh_t   *mesh)
{
  int  i, shift;

  int  idx_size = 0;
  int  *cell_shift = NULL;
  cs_sla_matrix_t  *m = NULL;

  const int  n_cells = mesh->n_cells;
  const int  n_i_faces = mesh->n_i_faces;
  const int  n_b_faces = mesh->n_b_faces;
  const int  n_faces = n_i_faces + n_b_faces;

  m = cs_sla_matrix_create(n_cells, n_faces, 1, CS_SLA_MAT_DEC, false);

  BFT_MALLOC(cell_shift, n_cells, int);
  for (i = 0; i < n_cells; i++)
    cell_shift[i] = 0;

  for (i = 0; i < n_b_faces; i++) {
    m->idx[mesh->b_face_cells[i]+1] += 1;
    idx_size += 1;
  }

  for (i = 0; i < n_i_faces; i++) {

    int  c1_id = mesh->i_face_cells[i][0];
    int  c2_id = mesh->i_face_cells[i][1];

    if (c1_id < n_cells)
      m->idx[c1_id+1] += 1, idx_size += 1;
    if (c2_id < n_cells)
      m->idx[c2_id+1] += 1, idx_size += 1;
  }

  for (i = 0; i < n_cells; i++)
    m->idx[i+1] += m->idx[i];

  assert(m->idx[n_cells] == idx_size);

  BFT_MALLOC(m->col, idx_size, int);
  BFT_MALLOC(m->sgn, idx_size, short int);

  for (i = 0; i < n_i_faces; i++) {

    int  c1_id = mesh->i_face_cells[i][0];
    int  c2_id = mesh->i_face_cells[i][1];

    if (c1_id < n_cells) { /* Don't want ghost cells */

      shift = m->idx[c1_id] + cell_shift[c1_id];
      m->col[shift] = i + 1;
      m->sgn[shift] = 1;
      cell_shift[c1_id] += 1;

    }

    if (c2_id < n_cells) { /* Don't want ghost cells */

      shift = m->idx[c2_id] + cell_shift[c2_id];
      m->col[shift] = i + 1;
      m->sgn[shift] = -1;
      cell_shift[c2_id] += 1;

    }

  } /* End of loop on internal faces */

  for (i = 0; i < n_b_faces; i++) {

    int  cid = mesh->b_face_cells[i];

    shift = m->idx[cid] + cell_shift[cid];
    m->col[shift] = n_i_faces + i + 1;
    m->sgn[shift] = 1;
    cell_shift[cid] += 1;

  } /* End of loop on border faces */

  /* Free memory */
  BFT_FREE(cell_shift);

  return m;
}

/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_build_f2e_connect(const cs_mesh_t     *mesh,
                   cs_cdo_connect_builder_t  *cbuilder)
{
  int  i, j, s, e, shift, v1_num, v2_num, e_num;

  int  *face_shift = NULL;
  cs_sla_matrix_t  *m = NULL;

  const int  n_i_faces = mesh->n_i_faces;
  const int  n_b_faces = mesh->n_b_faces;
  const int  n_faces = n_i_faces + n_b_faces;
  const int  n_edges = cbuilder->n_edges;

  m = cs_sla_matrix_create(n_faces, n_edges, 1, CS_SLA_MAT_DEC, false);

  /* Build index */
  BFT_MALLOC(face_shift, n_faces, int);
  for (i = 0; i < n_faces; i++)
    face_shift[i] = 0;

  for (i = 0; i < n_b_faces; i++) {
    s = mesh->b_face_vtx_idx[i], e = mesh->b_face_vtx_idx[i+1];
    m->idx[n_i_faces+i+1] += e-s;
  }

  for (i = 0; i < n_i_faces; i++) {
    s = mesh->i_face_vtx_idx[i], e = mesh->i_face_vtx_idx[i+1];
    m->idx[i+1] += e-s;
  }

  for (i = 0; i < n_faces; i++)
    m->idx[i+1] += m->idx[i];

  assert(m->idx[n_faces] ==
         mesh->i_face_vtx_idx[n_i_faces] + mesh->b_face_vtx_idx[n_b_faces]);

  /* Build matrix */
  BFT_MALLOC(m->col, m->idx[n_faces], int);
  BFT_MALLOC(m->sgn, m->idx[n_faces], short int);

  for (i = 0; i < n_b_faces; i++) {

    int  fid = n_i_faces + i;

    s = mesh->b_face_vtx_idx[i], e = mesh->b_face_vtx_idx[i+1];

    shift = m->idx[fid] + face_shift[fid];
    face_shift[fid] += 1;

    v1_num = mesh->b_face_vtx_lst[e-1] + 1;
    v2_num = mesh->b_face_vtx_lst[s] + 1;
    e_num = _get_edge_num(v1_num, v2_num, cbuilder);
    if (e_num < 0)
      m->col[shift] = -e_num, m->sgn[shift] = -1;
    else
      m->col[shift] = e_num, m->sgn[shift] = 1;

    for (j = s; j < e-1; j++) {

      shift++;
      face_shift[fid] += 1;
      v1_num = mesh->b_face_vtx_lst[j] + 1;
      v2_num = mesh->b_face_vtx_lst[j+1] + 1;
      e_num = _get_edge_num(v1_num, v2_num, cbuilder);
      if (e_num < 0)
        m->col[shift] = -e_num, m->sgn[shift] = -1;
      else
        m->col[shift] = e_num, m->sgn[shift] = 1;

    }

  } /* End of loop on border faces */

  for (i = 0; i < n_i_faces; i++) {

    s = mesh->i_face_vtx_idx[i], e = mesh->i_face_vtx_idx[i+1];

    shift = m->idx[i] + face_shift[i];
    face_shift[i] += 1;

    v1_num = mesh->i_face_vtx_lst[e-1] + 1;
    v2_num = mesh->i_face_vtx_lst[s] + 1;
    e_num = _get_edge_num(v1_num, v2_num, cbuilder);
    if (e_num < 0)
      m->col[shift] = -e_num, m->sgn[shift] = -1;
    else
      m->col[shift] = e_num, m->sgn[shift] = 1;

    for (j = s; j < e-1; j++) {

      shift++;
      face_shift[i] += 1;
      v1_num = mesh->i_face_vtx_lst[j] + 1;
      v2_num = mesh->i_face_vtx_lst[j+1] + 1;
      e_num = _get_edge_num(v1_num, v2_num, cbuilder);
      if (e_num < 0)
        m->col[shift] = -e_num, m->sgn[shift] = -1;
      else
        m->col[shift] = e_num, m->sgn[shift] = 1;

    }

  } /* End of loop on internal faces */

  /* Free memory */
  BFT_FREE(face_shift);

  return m;
}

/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_build_e2v_connect(cs_cdo_connect_builder_t  *cbuilder)
{
  int  i;

  cs_sla_matrix_t  *m = NULL;

  const int  n_vertices = cbuilder->n_vertices;
  const int  n_edges = cbuilder->n_edges;

  m = cs_sla_matrix_create(n_edges, n_vertices, 1, CS_SLA_MAT_DEC, false);

  /* Build index */
  m->idx[0] = 0;
  for (i = 0; i < n_edges; i++)
    m->idx[i+1] = m->idx[i] + 2;

  assert(m->idx[n_edges] == 2*n_edges);

  /* Build matrix */
  BFT_MALLOC(m->col, m->idx[n_edges], int);
  BFT_MALLOC(m->sgn, m->idx[n_edges], short int);

  for (i = 0; i < n_edges; i++) {

    m->col[2*i] = cbuilder->e2v_lst[2*i];
    m->sgn[2*i] = -1;
    m->col[2*i+1] = cbuilder->e2v_lst[2*i+1];
    m->sgn[2*i+1] = 1;

  }

  return m;
}

/*----------------------------------------------------------------------------*/

static void
_build_edge_connect(const cs_mesh_t      *mesh,
                    cs_cdo_connect_builder_t   *cbuilder)
{
  int  i, j, k, v1, v2, o1, o2, s, e, nfv, shift, s1, s2;

  int  n_edges = 0, n_init_edges = 0;
  int  n_max_face_vertices = 0;
  int  *f_vertices = NULL, *vtx_shift = NULL;
  int  *v2v_idx = NULL, *v2v_lst = NULL, *v2e_lst = NULL;
  int  *e2v_ref_lst = NULL;
  cs_gnum_t  *e2v_lst = NULL; /* Only because we have to use cs_order */
  cs_lnum_t  *order = NULL;

  const int n_vertices = mesh->n_vertices;
  const int n_i_faces = mesh->n_i_faces;
  const int n_b_faces = mesh->n_b_faces;


  for (i = 0; i < n_b_faces; i++)
    n_max_face_vertices = CS_MAX(n_max_face_vertices,
                                 mesh->b_face_vtx_idx[i+1]
                                 - mesh->b_face_vtx_idx[i]);

  for (i = 0; i < n_i_faces; i++)
    n_max_face_vertices = CS_MAX(n_max_face_vertices,
                                 mesh->i_face_vtx_idx[i+1]
                                 - mesh->i_face_vtx_idx[i]);

  BFT_MALLOC(f_vertices, n_max_face_vertices + 1, int);

  n_init_edges = mesh->b_face_vtx_idx[n_b_faces];
  n_init_edges += mesh->i_face_vtx_idx[n_i_faces];

  /* Build e2v_lst */
  BFT_MALLOC(e2v_lst, 2*n_init_edges, cs_gnum_t);

  shift = 0;
  for (i = 0; i < n_b_faces; i++) {

    s = mesh->b_face_vtx_idx[i], e = mesh->b_face_vtx_idx[i+1];
    nfv = e - s;

    for (j = s, k = 0; j < e; j++, k++)
      f_vertices[k] = mesh->b_face_vtx_lst[j] + 1;
    f_vertices[nfv] = mesh->b_face_vtx_lst[s] + 1;

    for (k = 0; k < nfv; k++) {

      v1 = f_vertices[k], v2 = f_vertices[k+1];
      if (v1 < v2)
        e2v_lst[2*shift] = v1, e2v_lst[2*shift+1] = v2;
      else
        e2v_lst[2*shift] = v2, e2v_lst[2*shift+1] = v1;
      shift++;

    }

  } /* End of loop on border faces */

  for (i = 0; i < n_i_faces; i++) {

    s = mesh->i_face_vtx_idx[i], e = mesh->i_face_vtx_idx[i+1];
    nfv = e - s;

    for (j = s, k = 0; j < e; j++, k++)
      f_vertices[k] = mesh->i_face_vtx_lst[j] + 1;
    f_vertices[nfv] = mesh->i_face_vtx_lst[s] + 1;

    for (k = 0; k < nfv; k++) {

      v1 = f_vertices[k], v2 = f_vertices[k+1];
      if (v1 < v2)
        e2v_lst[2*shift] = v1, e2v_lst[2*shift+1] = v2;
      else
        e2v_lst[2*shift] = v2, e2v_lst[2*shift+1] = v1;
      shift++;

    }

  } /* End of loop on interior faces */

  assert(shift == n_init_edges);

  BFT_MALLOC(order, n_init_edges, cs_lnum_t);
  cs_order_gnum_allocated_s(NULL, e2v_lst, 2, order, n_init_edges);

  BFT_MALLOC(v2v_idx, n_vertices + 1, int);
  for (i = 0; i < n_vertices + 1; i++)
    v2v_idx[i] = 0;

  if (n_init_edges > 0) {

    BFT_MALLOC(e2v_ref_lst, 2*n_init_edges, int);

    o1 = order[0];
    v1 = e2v_lst[2*o1];
    v2 = e2v_lst[2*o1+1];

    e2v_ref_lst[0] = v1;
    e2v_ref_lst[1] = v2;
    v2v_idx[v1] += 1;
    v2v_idx[v2] += 1;
    shift = 1;

    for (i = 1; i < n_init_edges; i++) {

      o1 = order[i-1];
      o2 = order[i];

      if (   e2v_lst[2*o1]   != e2v_lst[2*o2]
          || e2v_lst[2*o1+1] != e2v_lst[2*o2+1]) {

        v2v_idx[e2v_lst[2*o2]] += 1;
        v2v_idx[e2v_lst[2*o2+1]] += 1;
        e2v_ref_lst[2*shift] = e2v_lst[2*o2];
        e2v_ref_lst[2*shift+1] = e2v_lst[2*o2+1];
        shift++;

      }

    } /* End of loop on edges */

  } /* n_init_edges > 0 */

  n_edges = shift;

  for (i = 0; i < n_vertices; i++)
    v2v_idx[i+1] += v2v_idx[i];

  /* Free memory */
  BFT_FREE(e2v_lst);
  BFT_FREE(order);
  BFT_FREE(f_vertices);

  if (n_edges > 0) {

    BFT_MALLOC(v2v_lst, v2v_idx[n_vertices], int);
    BFT_MALLOC(v2e_lst, v2v_idx[n_vertices], int);
    BFT_MALLOC(vtx_shift, n_vertices, int);

    for (i = 0; i < n_vertices; i++)
      vtx_shift[i] = 0;

    for (i = 0; i < n_edges; i++) {

      v1 = e2v_ref_lst[2*i] - 1;
      v2 = e2v_ref_lst[2*i+1] - 1;
      s1 = v2v_idx[v1] + vtx_shift[v1];
      s2 = v2v_idx[v2] + vtx_shift[v2];
      vtx_shift[v1] += 1;
      vtx_shift[v2] += 1;
      v2v_lst[s1] = v2 + 1;
      v2v_lst[s2] = v1 + 1;

      if (v1 < v2)
        v2e_lst[s1] = i+1, v2e_lst[s2] = -(i+1);
      else
        v2e_lst[s1] = -(i+1), v2e_lst[s2] = i+1;

    } /* End of loop on edges */

    BFT_FREE(vtx_shift);

  } /* n_edges > 0 */

  /* Return pointers */
  cbuilder->n_edges = n_edges;
  cbuilder->e2v_lst = e2v_ref_lst;
  cbuilder->v2v_idx = v2v_idx;
  cbuilder->v2v_lst = v2v_lst;
  cbuilder->v2v_edge_lst = v2e_lst;

}

/*----------------------------------------------------------------------------*/

static cs_cdo_connect_builder_t *
_define_connect_builder(const cs_mesh_t  *mesh)
{
  cs_cdo_connect_builder_t  *connect_builder = NULL;

  BFT_MALLOC(connect_builder, 1, cs_cdo_connect_builder_t);

  connect_builder->n_vertices = mesh->n_vertices;
  connect_builder->n_cells = mesh->n_cells;
  connect_builder->n_faces = mesh->n_b_faces + mesh->n_i_faces;

  _create_c2f_connect(mesh, connect_builder);
  _build_edge_connect(mesh, connect_builder);

  return connect_builder;
}

/*----------------------------------------------------------------------------*/

static void
_free_connect_builder(cs_cdo_connect_builder_t  **pconnect_builder)
{
  cs_cdo_connect_builder_t  *connect_builder = *pconnect_builder;

  if (connect_builder == NULL)
    return;

  BFT_FREE(connect_builder->e2v_lst);
  BFT_FREE(connect_builder->v2v_idx);
  BFT_FREE(connect_builder->v2v_lst);
  BFT_FREE(connect_builder->v2v_edge_lst);

  BFT_FREE(connect_builder->c2f_idx);
  BFT_FREE(connect_builder->c2f_lst);

  BFT_FREE(connect_builder);

  *pconnect_builder = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build a cs_cdo_connect_t structure
 *
 * \param[in]  connect     pointer to the cs_cdo_connect_t struct.
 */
/*----------------------------------------------------------------------------*/

static void
_build_cdo_connect(cs_cdo_connect_t  *connect)
{

  cs_connect_index_t  *c2f = NULL, *f2e = NULL, *e2v = NULL;

  /* Mapping */
  c2f = cs_index_map(connect->c2f->n_rows,
                     connect->c2f->idx,
                     connect->c2f->col);

  f2e = cs_index_map(connect->f2e->n_rows,
                     connect->f2e->idx,
                     connect->f2e->col);


  e2v = cs_index_map(connect->e2v->n_rows,
                     connect->e2v->idx,
                     connect->e2v->col);

  /* Build new connectivity */
  connect->c2e = cs_index_convol(connect->e2v->n_rows, c2f, f2e);
  connect->c2v = cs_index_convol(connect->v2e->n_rows, connect->c2e, e2v);

  /* Sort list for each entry */
  cs_index_sort(connect->c2v);
  cs_index_sort(connect->c2e);

  /* Free mapped structures */
  cs_index_free(&c2f);
  cs_index_free(&f2e);
  cs_index_free(&e2v);

#if INNOV_CONNECT_DEBUG
  /* Dump for debugging pruposes */
  cs_index_dump("Connect-c2e.log", NULL, connect->c2e);
  cs_index_dump("Connect-c2v.log", NULL, connect->c2v);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute max number of entities by cell
 *
 * \param[in]  connect     pointer to the cs_cdo_connect_t struct.
 */
/*----------------------------------------------------------------------------*/

static void
_compute_max_ent(cs_cdo_connect_t  *connect)
{
  int  i, n_ent;

  /* Max number of faces for a cell */
  connect->n_max_fbyc = 0;
  if (connect->c2f != NULL) {
    for (i = 0; i < connect->c2f->n_rows; i++) {
      n_ent = connect->c2f->idx[i+1] - connect->c2f->idx[i];
      if (n_ent > connect->n_max_fbyc)
        connect->n_max_fbyc = n_ent;
    }
  }

  /* Max number of edges for a cell */
  connect->n_max_ebyc = 0;
  if (connect->c2e != NULL) {
    for (i = 0; i < connect->c2e->n; i++) {
      n_ent = connect->c2e->idx[i+1] - connect->c2e->idx[i];
      if (n_ent > connect->n_max_ebyc)
        connect->n_max_ebyc = n_ent;
    }
  }

  /* Max number of vertices for a cell */
  connect->n_max_vbyc = 0;
  if (connect->c2v != NULL) {
    for (i = 0; i < connect->c2v->n; i++) {
      n_ent = connect->c2v->idx[i+1] - connect->c2v->idx[i];
      if (n_ent > connect->n_max_vbyc)
        connect->n_max_vbyc = n_ent;
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocated and initialize a cs_connect_info_t structure
 *
 * \param[in]     n_elts    Size of the maximal set of entities related to
 *                          this structure
 *
 * \return  a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

static cs_connect_info_t *
_connect_info_create(cs_lnum_t     n_elts)
{
  cs_lnum_t  i;

  cs_connect_info_t  *info = NULL;

  if (n_elts < 1)
    return NULL;

  BFT_MALLOC(info, 1, cs_connect_info_t);

  BFT_MALLOC(info->flag, n_elts, short int);
  for (i = 0; i < n_elts; i++)
    info->flag[i] = 0;

  info->n = n_elts;
  info->n_in = 0;
  info->n_bd = 0;
  info->n_ii = 0;
  info->n_ib = 0;
  info->n_bb = 0;
  info->n_bi = 0;

  return info;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocated and initialize a cs_cdo_connect_info_t structure
 *
 * \param[in]     n_elts    Size of the maximal set of entities related to
 *                          this structure
 *
 * \return  a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

static cs_connect_info_t *
_connect_info_free(cs_connect_info_t    *info)
{
  if (info == NULL)
    return info;

  BFT_FREE(info->flag);
  BFT_FREE(info);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a status Int/Border 1st and 2nd level for sets of vertices,
 *         edges and faces
 *
 * \param[inout]     connect    pointer to a cs_cdo_connect_t struct.
 */
/*----------------------------------------------------------------------------*/

static void
_define_connect_info(cs_cdo_connect_t     *connect)
{
  int  i, j, nn, vid, eid, fid, cid, count;
  short int  flag1, flag2;

  cs_connect_index_t  *v2v = NULL, *v2e = NULL, *e2v = NULL;
  cs_connect_info_t  *vi = NULL, *ei = NULL, *fi = NULL, *ci = NULL;

  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_lnum_t  n_vertices = connect->v2e->n_rows;
  const cs_lnum_t  n_edges = connect->e2f->n_rows;
  const cs_lnum_t  n_faces = connect->f2e->n_rows;
  const cs_lnum_t  n_cells = connect->c2f->n_rows;

  /* Allocate info structures */
  vi = _connect_info_create(n_vertices);
  ei = _connect_info_create(n_edges);
  fi = _connect_info_create(n_faces);
  ci = _connect_info_create(n_cells);

  /* By default all entities are set "interior" */
  for (i = 0; i < n_vertices; i++)
    vi->flag[i] = CS_CDO_CONNECT_IN;

  for (i = 0; i < n_edges; i++)
    ei->flag[i] = CS_CDO_CONNECT_IN;

  for (i = 0; i < n_faces; i++)
    fi->flag[i] = CS_CDO_CONNECT_IN;

  for (i = 0; i < n_cells; i++)
    ci->flag[i] = CS_CDO_CONNECT_IN;

  /* Loop on border faces => flag all border entities */
  for (fid = m->n_i_faces; fid < n_faces; fid++) {

    fi->flag[fid] = CS_CDO_CONNECT_BD;
    assert(connect->f2c->idx[fid+1]-connect->f2c->idx[fid]==1);
    cid = connect->f2c->col[connect->f2c->idx[fid]]-1;
    ci->flag[cid] = CS_CDO_CONNECT_BD;

    for (i = connect->f2e->idx[fid]; i < connect->f2e->idx[fid+1]; i++) {

      eid = connect->f2e->col[i]-1;
      ei->flag[eid] = CS_CDO_CONNECT_BD;
      for (j = connect->e2v->idx[eid]; j < connect->e2v->idx[eid+1]; j++) {

        vid = connect->e2v->col[j]-1;
        vi->flag[vid] = CS_CDO_CONNECT_BD;

      } /* Loop on border vertices */

    } /* Loop on border edges */

  } /* Loop on border faces */

  /* Count number of border vertices */
  for (i = 0; i < n_vertices; i++)
    if (vi->flag[i] & CS_CDO_CONNECT_BD)
      vi->n_bd++;
  vi->n_in = vi->n - vi->n_bd;

  /* Count number of border edges */
  for (i = 0; i < n_edges; i++)
    if (ei->flag[i] & CS_CDO_CONNECT_BD)
      ei->n_bd++;
  ei->n_in = ei->n - ei->n_bd;

  /* Count number of border faces */
  for (i = 0; i < n_faces; i++)
    if (fi->flag[i] & CS_CDO_CONNECT_BD)
      fi->n_bd++;
  fi->n_in = fi->n - fi->n_bd;
  assert(m->n_i_faces == fi->n_in);

  /* Count number of border cells */
  for (i = 0; i < n_cells; i++)
    if (ci->flag[i] & CS_CDO_CONNECT_BD)
      ci->n_bd++;
  ci->n_in = ci->n - ci->n_bd;

  /* Build v -> v connectivity */
  v2e = cs_index_map(n_vertices, connect->v2e->idx, connect->v2e->col);
  e2v = cs_index_map(n_edges, connect->e2v->idx, connect->e2v->col);
  v2v = cs_index_convol(n_vertices, v2e, e2v);

  /* Compute second level interior/border for vertices */
  for (i = 0; i < n_vertices; i++) {

    nn = v2v->idx[i+1] - v2v->idx[i], count = 0;
    for (j = v2v->idx[i]; j < v2v->idx[i+1]; j++)
      if (vi->flag[v2v->lst[j]-1] & CS_CDO_CONNECT_BD)
        count++;

    if (vi->flag[i] & CS_CDO_CONNECT_BD) { /* Border vertices */
      if (count == nn) vi->flag[i] |= CS_CDO_CONNECT_BB, vi->n_bb++;
      else             vi->flag[i] |= CS_CDO_CONNECT_BI;
    }
    else if (vi->flag[i] & CS_CDO_CONNECT_IN) { /* Interior vertices */
      if (count == 0) vi->flag[i] |= CS_CDO_CONNECT_II, vi->n_ii++;
      else            vi->flag[i] |= CS_CDO_CONNECT_IB;
    }
    else
      bft_error(__FILE__, __LINE__, 0,
                _(" Vertex %d is neither interior nor border.\n"
                  " Stop execution\n"), i+1);

  } /* End of loop on vertices */

  vi->n_bi = vi->n_bd - vi->n_bb;
  vi->n_ib = vi->n_in - vi->n_ii;

  cs_index_free(&v2v);
  cs_index_free(&v2e);
  cs_index_free(&e2v);

  /* Set of edges */
  for (i = 0; i < n_edges; i++) {

    j = connect->e2v->idx[i];
    flag1 = vi->flag[connect->e2v->col[j  ]-1];
    flag2 = vi->flag[connect->e2v->col[j+1]-1];

    if (ei->flag[i] == CS_CDO_CONNECT_IN) {
      if ( (flag1 & CS_CDO_CONNECT_II) && (flag2 & CS_CDO_CONNECT_II) )
        ei->flag[i] |= CS_CDO_CONNECT_II, ei->n_ii++;
      else
        ei->flag[i] |= CS_CDO_CONNECT_IB;
    }
    else if (ei->flag[i] == CS_CDO_CONNECT_BD) {
      if ( (flag1 & CS_CDO_CONNECT_BB) && (flag2 & CS_CDO_CONNECT_BB) )
        ei->flag[i] |= CS_CDO_CONNECT_BB, ei->n_bb++;
      else
        ei->flag[i] |= CS_CDO_CONNECT_BI;
    }
    else
      bft_error(__FILE__, __LINE__, 0,
                _(" Edge %d is neither interior nor border.\n"
                  " Stop execution\n"), i+1);
  } /* Loop on edges */

  ei->n_ib = ei->n_in - ei->n_ii;
  ei->n_bi = ei->n_bd - ei->n_bb;

  /* Set of faces */
  for (fid = 0; fid < n_faces; fid++) {

    nn = connect->f2e->idx[fid+1] - connect->f2e->idx[fid];
    count = 0;
    for (j = connect->f2e->idx[fid]; j < connect->f2e->idx[fid+1]; j++) {
      if (ei->flag[connect->f2e->col[j]-1] & CS_CDO_CONNECT_BD)
        count++;
    }

    if (fi->flag[fid] & CS_CDO_CONNECT_IN) {
      if (count == 0) fi->flag[fid] |= CS_CDO_CONNECT_II;
      else            fi->flag[fid] |= CS_CDO_CONNECT_IB;
    }
    else
      /* Border faces are built from only border edges. Second level tag
         is therfore not useful */
      assert(fi->flag[fid] & CS_CDO_CONNECT_BD && count == nn);

} /* Loop on faces */

  fi->n_ib = fi->n_in - fi->n_ii;

  /* Return pointers */
  connect->v_info = vi;
  connect->e_info = ei;
  connect->f_info = fi;
  connect->c_info = ci;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  String related to flag in cs_cdo_connect_info_t
 *
 * \param[in]  flag     retrieve name for this flag
 */
/*----------------------------------------------------------------------------*/

const char *
cs_cdo_connect_flagname(short int  flag)
{
  short int  _flag = 0;

  /* Second level is prior */
  if (flag & CS_CDO_CONNECT_II) _flag = CS_CDO_CONNECT_II;
  if (flag & CS_CDO_CONNECT_IB) _flag = CS_CDO_CONNECT_IB;
  if (flag & CS_CDO_CONNECT_BB) _flag = CS_CDO_CONNECT_BB;
  if (flag & CS_CDO_CONNECT_BI) _flag = CS_CDO_CONNECT_BI;

  if (_flag == 0) { /* Second level */
    if (flag & CS_CDO_CONNECT_IN) _flag = CS_CDO_CONNECT_IN;
    if (flag & CS_CDO_CONNECT_BD) _flag = CS_CDO_CONNECT_BD;
  }

  switch (_flag) {

  case CS_CDO_CONNECT_BD:
    return " Bd ";
    break;
  case CS_CDO_CONNECT_IN:
    return " In ";
    break;
  case CS_CDO_CONNECT_II:
    return "InIn";
    break;
  case CS_CDO_CONNECT_IB:
    return "InBd";
    break;
  case CS_CDO_CONNECT_BI:
    return "BdIn";
    break;
  case CS_CDO_CONNECT_BB:
    return "BdBd";
    break;
  default:
    return "Full";

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a cs_cdo_connect_t structure
 *
 * \param[in]  m            pointer to a cs_mesh_t structure
 *
 * \return  a cs_cdo_connect_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_connect_t *
cs_cdo_connect_build(const cs_mesh_t      *m)
{
  cs_cdo_connect_t  *connect = NULL;
  cs_cdo_connect_builder_t  *connect_builder = NULL;

  /* Build the connectivity structure */
  BFT_MALLOC(connect, 1, cs_cdo_connect_t);

  /* Build DEC matrices related to connectivity */
  connect->c2f = _build_c2f_connect(m);
  connect->f2c = cs_sla_matrix_transpose(connect->c2f);

  connect_builder = _define_connect_builder(m);

  connect->f2e = _build_f2e_connect(m, connect_builder);
  connect->e2f = cs_sla_matrix_transpose(connect->f2e);
  connect->e2v = _build_e2v_connect(connect_builder);
  connect->v2e = cs_sla_matrix_transpose(connect->e2v);

  _free_connect_builder(&connect_builder);

  /* Build additional connectivity c2e, c2v */
  _build_cdo_connect(connect);

  /* Build status flag: interior/border and related connection to
     interior/border entities */
  _define_connect_info(connect);

  /* Max number of entities (vertices, edges and faces) by cell */
  _compute_max_ent(connect);

  connect->max_set_size = CS_MAX(connect->v2e->n_rows, connect->e2v->n_rows);
  connect->max_set_size = CS_MAX(connect->f2e->n_rows, connect->max_set_size);
  connect->max_set_size = CS_MAX(connect->c2f->n_rows, connect->max_set_size);

  return connect;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a cs_cdo_connect_t structure
 *
 * \param[in]  connect     pointer to the cs_cdo_connect_t struct. to destroy
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_cdo_connect_t *
cs_cdo_connect_free(cs_cdo_connect_t   *connect)
{
  if (connect == NULL)
    return connect;

  connect->v2e = cs_sla_matrix_free(connect->v2e);
  connect->e2f = cs_sla_matrix_free(connect->e2f);
  connect->e2v = cs_sla_matrix_free(connect->e2v);
  connect->f2e = cs_sla_matrix_free(connect->f2e);
  connect->f2c = cs_sla_matrix_free(connect->f2c);
  connect->c2f = cs_sla_matrix_free(connect->c2f);

  /* Specific CDO connectivity */
  cs_index_free(&(connect->c2e));
  cs_index_free(&(connect->c2v));

  connect->v_info = _connect_info_free(connect->v_info);
  connect->e_info = _connect_info_free(connect->e_info);
  connect->f_info = _connect_info_free(connect->f_info);
  connect->c_info = _connect_info_free(connect->c_info);

  BFT_FREE(connect);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Resume connectivity information
 *
 * \param[in]  connect     pointer to cs_cdo_connect_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_connect_resume(const cs_cdo_connect_t  *connect)
{
  cs_connect_info_t  *i = NULL;

  /* Output */
  bft_printf("\n Connectivity information:\n");
  bft_printf("  --dim-- max. number of faces by cell:    %4d\n",
             connect->n_max_fbyc);
  bft_printf("  --dim-- max. number of edges by cell:    %4d\n",
             connect->n_max_ebyc);
  bft_printf("  --dim-- max. number of vertices by cell: %4d\n",
             connect->n_max_vbyc);

  if (connect->v_info != NULL) {
    i = connect->v_info;
    bft_printf("\n");
    bft_printf("                     |   full  |  intern |  border |  in/in  |"
               "  in/bd  |  bd/bd  |  bd/in  |\n");
    bft_printf("  --dim-- n_vertices |"
               " %7d | %7d | %7d | %7d | %7d | %7d | %7d |\n",
               i->n, i->n_in, i->n_bd, i->n_ii, i->n_ib, i->n_bb, i->n_bi);
  }
  if (connect->e_info != NULL) {
    i = connect->e_info;
    bft_printf("  --dim-- n_edges    |"
               " %7d | %7d | %7d | %7d | %7d | %7d | %7d |\n",
               i->n, i->n_in, i->n_bd, i->n_ii, i->n_ib, i->n_bb, i->n_bi);
  }
  if (connect->f_info != NULL) {
    i = connect->f_info;
    bft_printf("  --dim-- n_faces    |"
               " %7d | %7d | %7d | %7d | %7d | %7d | %7d |\n",
               i->n, i->n_in, i->n_bd, i->n_ii, i->n_ib, i->n_bb, i->n_bi);
  }
  if (connect->c_info != NULL) {
    i = connect->c_info;
    bft_printf("  --dim-- n_cells    |"
               " %7d | %7d | %7d |\n", i->n, i->n_in, i->n_bd);
  }

  bft_printf("\n");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Dump a cs_cdo_connect_t structure
 *
 * \param[in]  connect     pointer to cs_cdo_connect_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_connect_dump(const cs_cdo_connect_t  *connect)
{
  FILE  *fdump = NULL;

  fdump = fopen("Innov_connect_dump.dat", "w");

  if (connect == NULL) {
    fprintf(fdump, "Empty structure.\n");
    fclose(fdump);
    return;
  }

  fprintf(fdump, "\n Connect structure: %p\n", (const void *)connect);

  /* Dump CONNECT matrices */
  cs_sla_matrix_dump("Connect c2f mat", fdump, connect->c2f);
  cs_sla_matrix_dump("Connect f2c mat", fdump, connect->f2c);
  cs_sla_matrix_dump("Connect f2e mat", fdump, connect->f2e);
  cs_sla_matrix_dump("Connect e2f mat", fdump, connect->e2f);
  cs_sla_matrix_dump("Connect e2v mat", fdump, connect->e2v);
  cs_sla_matrix_dump("Connect v2e mat", fdump, connect->v2e);

  /* Dump specific CDO connectivity */
  cs_index_dump("Connect c2e", fdump, connect->c2e);
  cs_index_dump("Connect c2v", fdump, connect->c2v);

  fclose(fdump);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Create an index structure of size n
 *
 * \param[in]  n     number of entries of the indexed list
 *
 * \return  a pointer to a cs_connect_index_t
 */
/*----------------------------------------------------------------------------*/

cs_connect_index_t *
cs_index_create(int  n)
{
  int  i;

  cs_connect_index_t  *x = NULL;

  BFT_MALLOC(x, 1, cs_connect_index_t);

  x->n = n;
  x->owner = true;
  x->lst = NULL;

  BFT_MALLOC(x->idx, n+1, int);
  for (i = 0; i < x->n + 1; i++)
    x->idx[i] = 0;

  return x;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Map arrays into an index structure of size n (owner = false)
 *
 * \param[in]  n     number of entries of the indexed list
 * \param[in]  idx   array of size n+1
 * \param[in]  lst   array of size idx[n]
 *
 * \return  a pointer to a cs_connect_index_t
 */
/*----------------------------------------------------------------------------*/

cs_connect_index_t *
cs_index_map(int    n,
             int   *idx,
             int   *lst)
{
  cs_connect_index_t  *x = NULL;

  BFT_MALLOC(x, 1, cs_connect_index_t);

  x->n = n;
  x->owner = false;
  x->idx = idx;
  x->lst = lst;

  return x;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Destroy a cs_connect_index_t structure
 *
 * \param[in]  pidx     pointer of pointer to a cs_connect_index_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_index_free(cs_connect_index_t   **pidx)
{
  cs_connect_index_t  *x = *pidx;

  if (x == NULL)
    return;

  if (x->owner) {
    BFT_FREE(x->idx);
    BFT_FREE(x->lst);
  }

  BFT_FREE(x);
  *pidx = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   From 2 indexes : A -> B and B -> C create a new index A -> C
 *
 * \param[in]  nc      number of elements in C set
 * \param[in]  xab     pointer to the index A -> B
 * \param[in]  xbc     pointer to the index B -> C
 *
 *\return  a pointer to the cs_connect_index_t structure A -> C
 */
/*----------------------------------------------------------------------------*/

cs_connect_index_t *
cs_index_convol(int                        nc,
                const cs_connect_index_t  *xab,
                const cs_connect_index_t  *xbc)
{
  int  i, j, k, aid, bid, cid, shift;

  int  *ctag = NULL;
  cs_connect_index_t  *xac = NULL;

  xac = cs_index_create(xab->n);

  BFT_MALLOC(ctag, nc, int);
  for (i = 0; i < nc; i++)
    ctag[i] = -1;

  /* Build index */
  for (aid = 0; aid < xab->n; aid++) {

    for (j = xab->idx[aid]; j < xab->idx[aid+1]; j++) {

      bid = xab->lst[j] - 1;
      for (k = xbc->idx[bid]; k < xbc->idx[bid+1]; k++) {

        cid = xbc->lst[k] - 1;
        if (ctag[cid] != aid) { /* Not tagged yet */
          ctag[cid] = aid;
          xac->idx[aid+1] += 1;
        }

      } /* End of loop on C elements */

    } /* End of loop on B elements */

  } /* End of loop on A elements */

  for (i = 0; i < xac->n; i++)
    xac->idx[i+1] = xac->idx[i+1] +  xac->idx[i];

  BFT_MALLOC(xac->lst, xac->idx[xac->n], int);

  /* Reset ctag */
  for (i = 0; i < nc; i++)
    ctag[i] = -1;

  /* Fill lst */
  shift = 0;
  for (aid = 0; aid < xab->n; aid++) {

    for (j = xab->idx[aid]; j < xab->idx[aid+1]; j++) {

      bid = xab->lst[j] - 1;
      for (k = xbc->idx[bid]; k < xbc->idx[bid+1]; k++) {

        cid = xbc->lst[k] - 1;
        if (ctag[cid] != aid) { /* Not tagged yet */
          ctag[cid] = aid;
          xac->lst[shift++] = cid + 1;
        }

      } /* End of loop on C elements */

    } /* End of loop on B elements */

  } /* End of loop on A elements */

  assert(shift == xac->idx[xac->n]);

  BFT_FREE(ctag);

  return xac;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   From a cs_connect_index_t A -> B create a new index B -> A
 *
 * \param[in]  nb     size of the "b" set
 * \param[in]  ab     pointer to the index A -> B
 *
 * \return  a new pointer to the cs_connect_index_t structure B -> A
 */
/*----------------------------------------------------------------------------*/

cs_connect_index_t *
cs_index_transpose(int                        nb,
                   const cs_connect_index_t  *ab)
{
  int  i, j, b_id, shift;
  int  *count = NULL;

  cs_connect_index_t  *ba = cs_index_create(nb);

  if (nb == 0)
    return ba;

  /* Build idx */
  for (i = 0; i < ab->n; i++)
    for (j = ab->idx[i]; j < ab->idx[i+1]; j++)
      ba->idx[ab->lst[j]] += 1;

  for (i = 0; i < ba->n; i++)
    ba->idx[i+1] += ba->idx[i];

  /* Build lst */
  BFT_MALLOC(ba->lst, ba->idx[ba->n], int);

  /* Allocate temporary buffer */
  BFT_MALLOC(count, nb, int);
  for (i = 0; i < nb; i++) count[i] = 0;

  for (i = 0; i < ab->n; i++) {
    for (j = ab->idx[i]; j < ab->idx[i+1]; j++) {
      b_id = ab->lst[j] - 1;
      shift = count[b_id] + ba->idx[b_id];
      ba->lst[shift] = i+1;
      count[b_id] += 1;
    }
  }

  /* Free temporary buffer */
  BFT_FREE(count);

  return ba;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Sort each sub-list related to an entry in a cs_connect_index_t
 *          structure
 *
 * \param[in]  x     pointer to a cs_connect_index_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_index_sort(cs_connect_index_t   *x)
{
  int  i;

  if (x == NULL)
    return;

  for (i = 0; i < x->n; i++)
    cs_sort_shell(x->idx[i], x->idx[i+1], x->lst);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Dump a cs_connect_index_t structure to a file or into the
 *          standard output
 *
 * \param[in]  name  name of the dump file. Can be set to NULL
 * \param[in]  _f    pointer to a FILE structure. Can be set to NULL.
 * \param[in]  x     pointer to a cs_connect_index_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_index_dump(const char           *name,
              FILE                 *_f,
              cs_connect_index_t   *x)
{
  int  i, j;

  FILE  *f = _f;
  _Bool  close_file = false;

  if (f == NULL) {
    if (name == NULL)
      f = stdout;
    else {
      f = fopen(name,"w");
      close_file = true;
    }
  }

  fprintf(f, "\n Dump cs_connect_index_t struct: %p (%s)\n",
          (const void *)x, name);

  if (x == NULL) {
    if (close_file) fclose(f);
    return;
  }

  fprintf(f, "  owner:             %6d\n", x->owner);
  fprintf(f, "  n_elts:            %6d\n", x->n);
  fprintf(f, "  lst_size:          %6d\n", x->idx[x->n]);

  for (i = 0; i < x->n; i++) {

    fprintf(f, "\n[%4d] ", i+1);
    for (j = x->idx[i]; j < x->idx[i+1]; j++)
      fprintf(f, "%5d |", x->lst[j]);

  } /* End of loop on elements */

  if (close_file)
    fclose(f);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
