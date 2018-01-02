/*============================================================================
 * Manage connectivity (Topological features of the mesh)
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

#include "bft_mem.h"

#include "cs_cdo.h"
#include "cs_log.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_sort.h"


/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_connect.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local macro and structure definitions
 *============================================================================*/

#define CS_CDO_CONNECT_DBG 0

/* Temporary structure to build edge/vertices connectivities */
typedef struct {

  cs_lnum_t  n_vertices;
  cs_lnum_t  n_edges;

  int  *e2v_lst;  /* Edge ref. definition (2*n_edges) */
  int  *v2v_idx;
  int  *v2v_lst;
  int  *v2v_edge_lst;

} _edge_builder_t;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a entry in the face --> edges connectivity
 *
 * \param[in]      shift     position where to add the new entry
 * \param[in]      v1_num    number of the first vertex
 * \param[in]      v2_num    number of the second vertex
 * \param[in]      builder   pointer to a _edge_builder_t structure
 * \param[in, out] f2e       face --> edges connectivity
 */
/*----------------------------------------------------------------------------*/

static void
_add_f2e_entry(cs_lnum_t                   shift,
               cs_lnum_t                   v1_num,
               cs_lnum_t                   v2_num,
               const _edge_builder_t      *builder,
               cs_sla_matrix_t            *f2e)
{
  cs_lnum_t  i;

  /* Sanity check */
  assert(v1_num > 0);
  assert(v2_num > 0);
  assert(builder != NULL);
  assert(builder->v2v_idx != NULL);
  assert(builder->v2v_idx[v1_num] > builder->v2v_idx[v1_num-1]);

  /* Get edge number */
  cs_lnum_t  edge_sgn_num = 0;
  cs_lnum_t  *v2v_idx = builder->v2v_idx, *v2v_lst = builder->v2v_lst;

  for (i = v2v_idx[v1_num-1]; i < v2v_idx[v1_num]; i++) {
    if (v2v_lst[i] == v2_num) {
      edge_sgn_num = builder->v2v_edge_lst[i];
      break;
    }
  }

  if (edge_sgn_num == 0)
    bft_error(__FILE__, __LINE__, 0,
              _(" The given couple of vertices (number): [%d, %d]\n"
                " is not defined in the edge structure.\n"), v1_num, v2_num);

  f2e->col_id[shift] = CS_ABS(edge_sgn_num) - 1;
  if (edge_sgn_num < 0)
    f2e->sgn[shift] = -1;
  else
    f2e->sgn[shift] = 1;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the face -> edges connectivity which is stored in a
 *         cs_sla_matrix_t structure
 *
 * \param[in]  m         pointer to a cs_mesh_t structure
 * \param[in]  builder   pointer to the _edge_builder_t structure
 *
 * \return a pointer to a new allocated cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_build_f2e_connect(const cs_mesh_t         *m,
                   const _edge_builder_t   *builder)
{
  int  i, j, s, e, shift, v1_num, v2_num;

  cs_sla_matrix_t  *f2e = NULL;

  const int  n_i_faces = m->n_i_faces;
  const int  n_b_faces = m->n_b_faces;
  const int  n_faces = n_i_faces + n_b_faces;
  const int  n_edges = builder->n_edges;

  f2e = cs_sla_matrix_create(n_faces, n_edges, 1, CS_SLA_MAT_DEC, false);

  /* Build index */
  for (i = 0; i < n_i_faces; i++)
    f2e->idx[i+1] += m->i_face_vtx_idx[i+1] - m->i_face_vtx_idx[i];
  for (i = 0, j=n_i_faces+1; i < n_b_faces; i++, j++)
    f2e->idx[j] += m->b_face_vtx_idx[i+1] - m->b_face_vtx_idx[i];
  for (i = 0; i < n_faces; i++)
    f2e->idx[i+1] += f2e->idx[i];

  assert(f2e->idx[n_faces]
         == m->i_face_vtx_idx[n_i_faces] + m->b_face_vtx_idx[n_b_faces]);

  /* Build matrix */
  BFT_MALLOC(f2e->col_id, f2e->idx[n_faces], cs_lnum_t);
  BFT_MALLOC(f2e->sgn, f2e->idx[n_faces], short int);

  /* Border faces */
  for (i = 0; i < n_b_faces; i++) {

    s = m->b_face_vtx_idx[i], e = m->b_face_vtx_idx[i+1];

    cs_lnum_t  f_id = n_i_faces + i;

    shift = f2e->idx[f_id];
    v1_num = m->b_face_vtx_lst[e-1] + 1;
    v2_num = m->b_face_vtx_lst[s] + 1;
    _add_f2e_entry(shift, v1_num, v2_num, builder, f2e);

    for (j = s; j < e-1; j++) {

      shift++;
      v1_num = m->b_face_vtx_lst[j] + 1;
      v2_num = m->b_face_vtx_lst[j+1] + 1;
      _add_f2e_entry(shift, v1_num, v2_num, builder, f2e);

    }

  } /* End of loop on border faces */

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    s = m->i_face_vtx_idx[f_id], e = m->i_face_vtx_idx[f_id+1];
    shift = f2e->idx[f_id];
    v1_num = m->i_face_vtx_lst[e-1] + 1;
    v2_num = m->i_face_vtx_lst[s] + 1;
    _add_f2e_entry(shift, v1_num, v2_num, builder, f2e);

    for (j = s; j < e-1; j++) {

      shift++;
      v1_num = m->i_face_vtx_lst[j] + 1;
      v2_num = m->i_face_vtx_lst[j+1] + 1;
      _add_f2e_entry(shift, v1_num, v2_num, builder, f2e);

    }

  } /* End of loop on internal faces */

  return f2e;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the edge -> vertices connectivity which is stored in a
 *         cs_sla_matrix_t structure
 *
 * \param[in]  builder   pointer to the _edge_builder_t structure
 *
 * \return a pointer to a new allocated cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_build_e2v_connect(const _edge_builder_t  *builder)
{
  int  i;

  cs_sla_matrix_t  *e2v = NULL;

  const int  n_vertices = builder->n_vertices;
  const int  n_edges = builder->n_edges;

  e2v = cs_sla_matrix_create(n_edges, n_vertices, 1, CS_SLA_MAT_DEC, false);

  /* Build index */
  e2v->idx[0] = 0;
  for (i = 0; i < n_edges; i++)
    e2v->idx[i+1] = e2v->idx[i] + 2;

  assert(e2v->idx[n_edges] == 2*n_edges);

  /* Build matrix */
  BFT_MALLOC(e2v->col_id, e2v->idx[n_edges], cs_lnum_t);
  BFT_MALLOC(e2v->sgn, e2v->idx[n_edges], short int);

  for (i = 0; i < n_edges; i++) {

    e2v->col_id[2*i] = builder->e2v_lst[2*i] - 1;
    e2v->sgn[2*i] = -1;
    e2v->col_id[2*i+1] = builder->e2v_lst[2*i+1] - 1;
    e2v->sgn[2*i+1] = 1;

  }

  return e2v;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and define a _edge_builder_t structure
 *
 * \param[in]  m   pointer to the cs_mesh_t structure
 *
 * \return a pointer to a new allocated _edge_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static _edge_builder_t *
_create_edge_builder(const cs_mesh_t  *m)
{
  int  i, j, k, v1, v2, o1, o2, s, e, nfv, shift, s1, s2;

  int  n_edges = 0, n_init_edges = 0;
  int  n_max_face_vertices = 0;
  cs_lnum_t  *f_vertices = NULL, *vtx_shift = NULL;
  cs_lnum_t  *v2v_idx = NULL, *v2v_lst = NULL, *v2e_lst = NULL;
  int  *e2v_ref_lst = NULL;
  cs_gnum_t  *e2v_lst = NULL; /* Only because we have to use cs_order */
  cs_lnum_t  *order = NULL;

  _edge_builder_t  *builder = NULL;

  const int n_vertices = m->n_vertices;
  const int n_i_faces = m->n_i_faces;
  const int n_b_faces = m->n_b_faces;

  /* Compute max. number of vertices by face */
  for (i = 0; i < n_b_faces; i++)
    n_max_face_vertices = CS_MAX(n_max_face_vertices,
                                 m->b_face_vtx_idx[i+1] - m->b_face_vtx_idx[i]);
  for (i = 0; i < n_i_faces; i++)
    n_max_face_vertices = CS_MAX(n_max_face_vertices,
                                 m->i_face_vtx_idx[i+1] - m->i_face_vtx_idx[i]);

  BFT_MALLOC(f_vertices, n_max_face_vertices + 1, cs_lnum_t);

  n_init_edges = m->b_face_vtx_idx[n_b_faces];
  n_init_edges += m->i_face_vtx_idx[n_i_faces];

  /* Build e2v_lst */
  BFT_MALLOC(e2v_lst, 2*n_init_edges, cs_gnum_t);

  shift = 0;
  for (i = 0; i < n_b_faces; i++) {

    s = m->b_face_vtx_idx[i], e = m->b_face_vtx_idx[i+1];
    nfv = e - s;

    for (j = s, k = 0; j < e; j++, k++)
      f_vertices[k] = m->b_face_vtx_lst[j] + 1;
    f_vertices[nfv] = m->b_face_vtx_lst[s] + 1;

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

    s = m->i_face_vtx_idx[i], e = m->i_face_vtx_idx[i+1];
    nfv = e - s;

    for (j = s, k = 0; j < e; j++, k++)
      f_vertices[k] = m->i_face_vtx_lst[j] + 1;
    f_vertices[nfv] = m->i_face_vtx_lst[s] + 1;

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

    BFT_MALLOC(v2v_lst, v2v_idx[n_vertices], cs_lnum_t);
    BFT_MALLOC(v2e_lst, v2v_idx[n_vertices], cs_lnum_t);
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
  BFT_MALLOC(builder, 1, _edge_builder_t);

  builder->n_vertices = n_vertices;
  builder->n_edges = n_edges;
  builder->e2v_lst = e2v_ref_lst;
  builder->v2v_idx = v2v_idx;
  builder->v2v_lst = v2v_lst;
  builder->v2v_edge_lst = v2e_lst;

  return builder;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a _edge_builder structure
 *
 * \param[in]  p_builder   pointer to the _edge_builder structure pointer
 */
/*----------------------------------------------------------------------------*/

static void
_free_edge_builder(_edge_builder_t  **p_builder)
{
  _edge_builder_t  *_builder = *p_builder;

  if (_builder == NULL)
    return;

  BFT_FREE(_builder->e2v_lst);
  BFT_FREE(_builder->v2v_idx);
  BFT_FREE(_builder->v2v_lst);
  BFT_FREE(_builder->v2v_edge_lst);

  BFT_FREE(_builder);

  *p_builder = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the cell -> faces connectivity which is stored in a
 *         cs_sla_matrix_t structure
 *
 * \param[in]  mesh      pointer to a cs_mesh_t structure
 *
 * \return a pointer to a new allocated cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_build_c2f_connect(const cs_mesh_t   *mesh)
{
  int  idx_size = 0;
  int  *cell_shift = NULL;
  cs_sla_matrix_t  *c2f = NULL;

  const int  n_cells = mesh->n_cells;
  const int  n_i_faces = mesh->n_i_faces;
  const int  n_b_faces = mesh->n_b_faces;
  const int  n_faces = n_i_faces + n_b_faces;

  c2f = cs_sla_matrix_create(n_cells, n_faces, 1, CS_SLA_MAT_DEC, false);

  BFT_MALLOC(cell_shift, n_cells, int);
# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_cells; i++)
    cell_shift[i] = 0;

  for (cs_lnum_t i = 0; i < n_b_faces; i++) {
    c2f->idx[mesh->b_face_cells[i]+1] += 1;
    idx_size += 1;
  }

  for (cs_lnum_t i = 0; i < n_i_faces; i++) {

    int  c1_id = mesh->i_face_cells[i][0];
    int  c2_id = mesh->i_face_cells[i][1];

    if (c1_id < n_cells) // cell owned by the local rank
      c2f->idx[c1_id+1] += 1, idx_size += 1;
    if (c2_id < n_cells) // cell owned by the local rank
      c2f->idx[c2_id+1] += 1, idx_size += 1;
  }

  for (cs_lnum_t i = 0; i < n_cells; i++)
    c2f->idx[i+1] += c2f->idx[i];

  assert(c2f->idx[n_cells] == idx_size);

  BFT_MALLOC(c2f->col_id, idx_size, cs_lnum_t);
  BFT_MALLOC(c2f->sgn, idx_size, short int);

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    const cs_lnum_t  c1_id = mesh->i_face_cells[f_id][0];
    const cs_lnum_t  c2_id = mesh->i_face_cells[f_id][1];

    if (c1_id < n_cells) { /* Don't want ghost cells */

      const cs_lnum_t  shift = c2f->idx[c1_id] + cell_shift[c1_id];
      c2f->col_id[shift] = f_id;
      c2f->sgn[shift] = 1;
      cell_shift[c1_id] += 1;

    }

    if (c2_id < n_cells) { /* Don't want ghost cells */

      const cs_lnum_t  shift = c2f->idx[c2_id] + cell_shift[c2_id];
      c2f->col_id[shift] = f_id;
      c2f->sgn[shift] = -1;
      cell_shift[c2_id] += 1;

    }

  } /* End of loop on internal faces */

  for (cs_lnum_t  f_id = 0; f_id < n_b_faces; f_id++) {

    const cs_lnum_t  c_id = mesh->b_face_cells[f_id];
    const cs_lnum_t  shift = c2f->idx[c_id] + cell_shift[c_id];

    c2f->col_id[shift] = n_i_faces + f_id;
    c2f->sgn[shift] = 1;
    cell_shift[c_id] += 1;

  } /* End of loop on border faces */

  /* Free memory */
  BFT_FREE(cell_shift);

  return c2f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build additional connectivities for accessing geometrical quantities
 *        c2e: cell --> edges connectivity
 *        c2v: cell --> vertices connectivity
 *
 * \param[in, out]  connect     pointer to the cs_cdo_connect_t struct.
 */
/*----------------------------------------------------------------------------*/

static void
_build_additional_connect(cs_cdo_connect_t  *connect)
{
  cs_connect_index_t  *c2f = cs_index_map(connect->c2f->n_rows,
                                          connect->c2f->idx,
                                          connect->c2f->col_id);
  cs_connect_index_t  *f2e = cs_index_map(connect->f2e->n_rows,
                                          connect->f2e->idx,
                                          connect->f2e->col_id);
  cs_connect_index_t  *e2v = cs_index_map(connect->e2v->n_rows,
                                          connect->e2v->idx,
                                          connect->e2v->col_id);

  /* Build new connectivity */
  connect->c2e = cs_index_compose(connect->e2v->n_rows, c2f, f2e);
  connect->c2v = cs_index_compose(connect->v2e->n_rows, connect->c2e, e2v);

  /* Sort list for each entry */
  cs_index_sort(connect->c2v);
  cs_index_sort(connect->c2e);

  /* Free mapped structures */
  cs_index_free(&c2f);
  cs_index_free(&f2e);
  cs_index_free(&e2v);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute max number of entities by cell and the max range between
 *        the min. id and the max.id for edges and vertices
 *
 * \param[in]       m         pointer to a cs_mesh_t structure
 * \param[in, out]  connect   pointer to the cs_cdo_connect_t struct.
 */
/*----------------------------------------------------------------------------*/

static void
_compute_max_ent(const cs_mesh_t      *m,
                 cs_cdo_connect_t     *connect)
{
  assert(connect != NULL && m != NULL);
  assert(connect->c2v != NULL && connect->c2e != NULL && connect->f2c != NULL);
  assert(connect->f2e != NULL);

  const cs_lnum_t  n_vertices = connect->n_vertices;

  short int  *v_count = NULL;
  BFT_MALLOC(v_count, n_vertices, short int);
# pragma omp parallel for if (n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_vertices; i++) v_count[i] = 0;

  const cs_lnum_t  n_cells = connect->n_cells;
  const cs_lnum_t  n_edges = connect->n_edges;
  const cs_connect_index_t  *c2v = connect->c2v;
  const cs_connect_index_t  *c2e = connect->c2e;
  const cs_sla_matrix_t  *e2v = connect->e2v;

  assert(c2e->n == n_cells);
  assert(c2v->n == n_cells);

  int  n_max_vc = 0, n_max_ec = 0, n_max_fc = 0;
  int  n_max_v2ec = 0, n_max_v2fc = 0, n_max_vf = 0;

  connect->e_max_cell_range = 0;
  connect->v_max_cell_range = 0;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    /* Vertices */
    const cs_lnum_t  *c2v_idx = c2v->idx + c_id;
    const cs_lnum_t  *c2v_ids = c2v->ids + c2v_idx[0];
    const int n_vc = c2v_idx[1] - c2v_idx[0];

    cs_lnum_t  min_id = n_vertices, max_id = 0;
    for (short int v = 0; v < n_vc; v++) {
      if (c2v_ids[v] < min_id) min_id = c2v_ids[v];
      if (c2v_ids[v] > max_id) max_id = c2v_ids[v];
    }
    cs_lnum_t  _range = max_id - min_id;

    if (n_vc > n_max_vc) n_max_vc = n_vc;
    if (connect->v_max_cell_range < _range)
      connect->v_max_cell_range = _range;

    /* Edges */
    const cs_lnum_t  *c2e_idx = c2e->idx + c_id;
    const cs_lnum_t  *c2e_ids = c2e->ids + c2e_idx[0];
    const int n_ec = c2e_idx[1] - c2e_idx[0];

    min_id = n_edges, max_id = 0;
    for (short int e = 0; e < n_ec; e++) {
      if (c2e_ids[e] < min_id) min_id = c2e_ids[e];
      if (c2e_ids[e] > max_id) max_id = c2e_ids[e];
    }
    _range = max_id - min_id;

    if (n_ec > n_max_ec) n_max_ec = n_ec;
    if (connect->e_max_cell_range < _range)
      connect->e_max_cell_range = _range;

    for (short int e = 0; e < n_ec; e++) {

      const cs_lnum_t  *e2v_ids = e2v->col_id + 2*c2e_ids[e];

      v_count[e2v_ids[0]] += 1;
      v_count[e2v_ids[1]] += 1;

    }

    /* Update n_max_v2ec and reset v_count */
    for (short int v = 0; v < n_vc; v++) {

      const cs_lnum_t  v_id = c2v_ids[v];
      if (v_count[v_id] > n_max_v2ec) n_max_v2ec = v_count[v_id];
      v_count[v_id] = 0; // reset

    }

    const cs_lnum_t  *c2f_idx = connect->c2f->idx + c_id;
    const cs_lnum_t  *c2f_ids = connect->c2f->col_id + c2f_idx[0];
    const int  n_fc = c2f_idx[1] - c2f_idx[0];

    if (n_fc > n_max_fc) n_max_fc = n_fc;

    for (short int f = 0; f < n_fc; f++) {

      if (c2f_ids[f] < m->n_i_faces) { // Interior face

        const cs_lnum_t  *f2v_idx = m->i_face_vtx_idx + f;
        const cs_lnum_t  *f2v_ids = m->i_face_vtx_lst + f2v_idx[0];
        const int  n_vf = f2v_idx[1] - f2v_idx[0];

        if (n_vf > n_max_vf) n_max_vf = n_vf;
        for (short int v = 0; v < n_vf; v++)
          v_count[f2v_ids[v]] += 1;

      }
      else { // Border face

        const cs_lnum_t  *f2v_idx =
          m->b_face_vtx_idx + c2f_ids[f] - m->n_i_faces;
        const cs_lnum_t  *f2v_ids = m->b_face_vtx_lst + f2v_idx[0];
        const int  n_vf = f2v_idx[1] - f2v_idx[0];

        if (n_vf > n_max_vf) n_max_vf = n_vf;
        for (short int v = 0; v < n_vf; v++)
          v_count[f2v_ids[v]] += 1;

      }

    } // Loop on cell faces

    /* Update n_max_v2fc and reset  v_count */
    for (short int v = 0; v < n_vc; v++) {

      const cs_lnum_t  v_id = c2v_ids[v];
      if (v_count[v_id] > n_max_v2ec) n_max_v2fc = v_count[v_id];
      v_count[v_id] = 0; // reset

    }

  } // Loop on cells

  BFT_FREE(v_count);

  /* Store computed values */
  connect->n_max_vbyc = n_max_vc;   // Max number of vertices for a cell
  connect->n_max_ebyc = n_max_ec;   // Max number of edges for a cell
  connect->n_max_fbyc = n_max_fc;   // Max number of faces for a cell
  connect->n_max_v2ec = n_max_v2ec;
  connect->n_max_v2fc = n_max_v2fc;
  connect->n_max_vbyf = n_max_vf;   // Max number of vertices in a face
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Associate to each cell a type of element (fvm_element_t)
 *
 * \param[in]  c_id      cell id
 * \param[in]  connect   pointer to a cs_cdo_connect_t struct.
 *
 * \return  type of element for this cell
 */
/*----------------------------------------------------------------------------*/

static fvm_element_t
_get_cell_type(cs_lnum_t                 c_id,
               const cs_cdo_connect_t   *connect)
{
  fvm_element_t  ret_type = FVM_CELL_POLY; // Default value

  int  n_vc = connect->c2v->idx[c_id+1] - connect->c2v->idx[c_id];
  int  n_ec = connect->c2e->idx[c_id+1] - connect->c2e->idx[c_id];
  int  n_fc = connect->c2f->idx[c_id+1] - connect->c2f->idx[c_id];

  /* Tetrahedron */
  if (n_vc == 4 && n_ec == 6 && n_fc == 4)
    ret_type = FVM_CELL_TETRA;

  /* Pyramid */
  else if (n_vc == 5 && n_ec == 8 && n_fc == 5)
    ret_type = FVM_CELL_PYRAM;

  /* Prism ? */
  else if (n_vc == 6 && n_ec == 9 && n_fc == 5) { // Potentially a prism

    int  count[2] = {0, 0};

    /* Loop on cell faces */
    for (cs_lnum_t i = connect->c2f->idx[c_id]; i < connect->c2f->idx[c_id+1];
         i++) {

      cs_lnum_t  f_id = connect->c2f->col_id[i];

      if (connect->f2e->idx[f_id+1] - connect->f2e->idx[f_id] == 4) // Quad
        count[1] += 1;
      if (connect->f2e->idx[f_id+1] - connect->f2e->idx[f_id] == 3) // Tria
        count[0] += 1;

      if (count[0] == 2 && count[1] == 3)
        ret_type = FVM_CELL_PRISM;
    }

  }

  /* Hexahedron ? */
  else if (n_vc == 8 && n_ec == 12 && n_fc == 6) { // Potentially a hexahedron

    _Bool  is_hexa = true;

    /* Loop on cell faces */
    for (cs_lnum_t i = connect->c2f->idx[c_id]; i < connect->c2f->idx[c_id+1];
         i++) {

      cs_lnum_t  f_id = connect->c2f->col_id[i];

      if (connect->f2e->idx[f_id+1] - connect->f2e->idx[f_id] != 4) {
        is_hexa = false;
        break;
      }

    }

    if (is_hexa)
      ret_type = FVM_CELL_HEXA;

  }

  return ret_type;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and define a new cs_cdo_connect_t structure
 *        Range sets related to vertices and faces are computed inside and
 *        set as members of the cs_mesh_t structure
 *
 * \param[in, out]  mesh          pointer to a cs_mesh_t structure
 * \param[in]       scheme_flag   flag storing requested space schemes
 *
 * \return  a pointer to a cs_cdo_connect_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_connect_t *
cs_cdo_connect_init(cs_mesh_t      *mesh,
                    cs_flag_t       scheme_flag)
{
  cs_timer_t t0 = cs_timer_time();

  cs_cdo_connect_t  *connect = NULL;

  /* Build the connectivity structure */
  BFT_MALLOC(connect, 1, cs_cdo_connect_t);

  /* Build the cell --> faces connectivity */
  connect->c2f = _build_c2f_connect(mesh);

  /* Build the face --> cells connectivity */
  connect->f2c = cs_sla_matrix_transpose(connect->c2f);

  /* Build the face --> edges connectivity */
  _edge_builder_t  *builder = _create_edge_builder(mesh);
  connect->f2e = _build_f2e_connect(mesh, builder);

  /* Build the edge --> faces connectivity */
  connect->e2f = cs_sla_matrix_transpose(connect->f2e);

  /* Build the edge --> vertices connectivity */
  connect->e2v = _build_e2v_connect(builder);

  /* Build the vertex --> edges connectivity */
  connect->v2e = cs_sla_matrix_transpose(connect->e2v);

  _free_edge_builder(&builder);

  /* Build additional connectivity c2e, c2v (used to access dual faces and
     dual volumes) */
  _build_additional_connect(connect);

  const cs_lnum_t  n_vertices = connect->v2e->n_rows;
  const cs_lnum_t  n_edges = connect->e2f->n_rows;
  const cs_lnum_t  n_faces = connect->f2e->n_rows;
  const cs_lnum_t  n_cells = connect->c2f->n_rows;

  connect->n_vertices = n_vertices;
  connect->n_edges = n_edges;
  connect->n_faces[0] = n_faces;
  connect->n_faces[1] = mesh->n_b_faces;
  connect->n_faces[2] = mesh->n_i_faces;
  connect->n_cells = n_cells;

  /* Build the cell flag and associate a cell type to each cell */
  BFT_MALLOC(connect->cell_flag, n_cells, cs_flag_t);
  BFT_MALLOC(connect->cell_type, n_cells, fvm_element_t);

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    connect->cell_flag[c_id] = 0;
    connect->cell_type[c_id] = _get_cell_type(c_id, connect);
  }

  /* Loop on border faces and flag boundary cells */
  cs_lnum_t  *c_ids = connect->f2c->col_id + connect->f2c->idx[mesh->n_i_faces];
  for (cs_lnum_t f_id = 0; f_id < mesh->n_b_faces; f_id++)
    connect->cell_flag[c_ids[f_id]] = CS_FLAG_BOUNDARY;

  /* Max number of entities (vertices, edges and faces) by cell */
  _compute_max_ent(mesh, connect);

  connect->v_rs = NULL;
  connect->f_rs = NULL;

  if ((scheme_flag & CS_SCHEME_FLAG_CDOVB) ||
      (scheme_flag & CS_SCHEME_FLAG_CDOVCB)) {

    /* Vertex range set */
    cs_range_set_t  *v_rs = cs_range_set_create(mesh->vtx_interfaces,
                                                NULL,
                                                n_vertices,
                                                false,   // TODO: Ask Yvan
                                                0);      // g_id_base
    mesh->vtx_range_set = v_rs;
    connect->v_rs = v_rs;

  }

  if ((scheme_flag & CS_SCHEME_FLAG_CDOFB) ||
      (scheme_flag & CS_SCHEME_FLAG_HHO)) {

    /* Face range set */
    connect->f_rs = NULL; // TODO

  }

  /* Monitoring */
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_t  time_count = cs_timer_diff(&t0, &t1);
  cs_log_printf(CS_LOG_PERFORMANCE, " %-35s %9.3f s\n",
                "<CDO/Connectivity> Runtime", time_count.wall_nsec*1e-9);

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

  BFT_FREE(connect->cell_type);
  BFT_FREE(connect->cell_flag);

  /* Structures for parallelism */
  cs_range_set_destroy(&(connect->f_rs));

  BFT_FREE(connect);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of connectivity information
 *
 * \param[in]  connect     pointer to cs_cdo_connect_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_connect_summary(const cs_cdo_connect_t  *connect)
{
  cs_lnum_t  n_max_entbyc[5] = {connect->n_max_fbyc,
                                connect->n_max_ebyc,
                                connect->n_max_vbyc,
                                connect->v_max_cell_range,
                                connect->e_max_cell_range};

  if (cs_glob_n_ranks > 1)
    cs_parall_max(5, CS_LNUM_TYPE, n_max_entbyc);

  /* Output */
  cs_log_printf(CS_LOG_DEFAULT, "\n Connectivity information:\n");
  cs_log_printf(CS_LOG_DEFAULT,
                " --dim-- max. number of faces by cell:    %4d\n",
                n_max_entbyc[0]);
  cs_log_printf(CS_LOG_DEFAULT,
                " --dim-- max. number of edges by cell:    %4d\n",
                n_max_entbyc[1]);
  cs_log_printf(CS_LOG_DEFAULT,
                " --dim-- max. number of vertices by cell: %4d\n",
                n_max_entbyc[2]);
  cs_log_printf(CS_LOG_DEFAULT,
                " --dim-- max. vertex range for a cell:      %d\n",
                n_max_entbyc[3]);
  cs_log_printf(CS_LOG_DEFAULT,
                " --dim-- max. edge range for a cell:        %d\n\n",
                n_max_entbyc[4]);

  /* Information about the element types */
  cs_gnum_t  n_type_cells[FVM_N_ELEMENT_TYPES];
  for (int i = 0; i < FVM_N_ELEMENT_TYPES; i++)
    n_type_cells[i] = 0;

  for (cs_lnum_t i = 0; i < connect->n_cells; i++)
    n_type_cells[connect->cell_type[i]] += 1;

  if (cs_glob_n_ranks > 1)
    cs_parall_sum(FVM_N_ELEMENT_TYPES, CS_GNUM_TYPE, n_type_cells);

  cs_log_printf(CS_LOG_DEFAULT,
                " --dim-- number of tetrahedra: %8lu\n",
                n_type_cells[FVM_CELL_TETRA]);
  cs_log_printf(CS_LOG_DEFAULT,
                " --dim-- number of pyramids:   %8lu\n",
                n_type_cells[FVM_CELL_PYRAM]);
  cs_log_printf(CS_LOG_DEFAULT,
                " --dim-- number of prisms:     %8lu\n",
                n_type_cells[FVM_CELL_PRISM]);
  cs_log_printf(CS_LOG_DEFAULT,
                " --dim-- number of hexahedra:  %8lu\n",
                n_type_cells[FVM_CELL_HEXA]);
  cs_log_printf(CS_LOG_DEFAULT,
                " --dim-- number of polyhedra:  %8lu\n\n",
                n_type_cells[FVM_CELL_POLY]);

#if CS_CDO_CONNECT_DBG > 0 && defined(DEBUG) && !defined(NDEBUG)
  cs_cdo_connect_dump(connect);
#endif
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
  int  lname = strlen("DumpConnect.dat") + 1;

  /* Define the name of the dump file */
  char *fname = NULL;
  if (cs_glob_n_ranks > 1) {
    lname += 6;
    BFT_MALLOC(fname, lname, char);
    sprintf(fname, "DumpConnect.%05d.dat", cs_glob_rank_id);
  }
  else {
    BFT_MALLOC(fname, lname, char);
    sprintf(fname, "DumpConnect.dat");
  }
  FILE  *fdump = fopen(fname, "w");

  if (connect == NULL) {
    fprintf(fdump, "Empty structure.\n");
    fclose(fdump);
    return;
  }

  fprintf(fdump, "\n Connect structure: %p\n", (const void *)connect);

  cs_sla_matrix_dump("Connect c2f mat", fdump, connect->c2f);
  cs_sla_matrix_dump("Connect f2e mat", fdump, connect->f2e);
  cs_sla_matrix_dump("Connect e2v mat", fdump, connect->e2v);

  cs_sla_matrix_dump("Connect f2c mat", fdump, connect->f2c);

  cs_sla_matrix_dump("Connect e2f mat", fdump, connect->e2f);
  cs_sla_matrix_dump("Connect v2e mat", fdump, connect->v2e);

  /* Dump specific CDO connectivity */
  cs_index_dump("Connect c2e", fdump, connect->c2e);
  cs_index_dump("Connect c2v", fdump, connect->c2v);

  fclose(fdump);
  BFT_FREE(fname);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
