/*============================================================================
 * Manage connectivity (Topological features of the mesh)
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include "fvm_io_num.h"

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

  cs_lnum_t  *e2v_lst;  /* Edge ref. definition (2*n_edges) */
  cs_lnum_t  *v2v_idx;
  cs_lnum_t  *v2v_lst;
  cs_lnum_t  *v2v_edge_lst;

} _edge_builder_t;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the cell -> faces connectivity which is stored in a
 *         cs_adjacency_t structure
 *         Interior and border faces are both considered. Border face ids are
 *         shift by n_i_faces.
 *
 * \param[in]  mesh      pointer to a cs_mesh_t structure
 *
 * \return a pointer to a new allocated cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_adjacency_t *
_build_c2f_connect(const cs_mesh_t   *mesh)
{
  cs_lnum_t  *cell_shift = NULL;
  cs_adjacency_t  *c2f = NULL;

  const cs_lnum_t  n_cells = mesh->n_cells;
  const cs_lnum_t  n_i_faces = mesh->n_i_faces;
  const cs_lnum_t  n_b_faces = mesh->n_b_faces;

  c2f = cs_adjacency_create(CS_ADJACENCY_SIGNED, // flag
                            -1,                  // no stride
                            n_cells);

  BFT_MALLOC(cell_shift, n_cells, cs_lnum_t);
# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_cells; i++) cell_shift[i] = 0;

  /* Build index */
  for (cs_lnum_t i = 0; i < n_b_faces; i++)
    c2f->idx[mesh->b_face_cells[i]+1] += 1;

  for (cs_lnum_t i = 0; i < n_i_faces; i++) {

    cs_lnum_t  c1_id = mesh->i_face_cells[i][0];
    cs_lnum_t  c2_id = mesh->i_face_cells[i][1];

    if (c1_id < n_cells) // cell owned by the local rank
      c2f->idx[c1_id+1] += 1;
    if (c2_id < n_cells) // cell owned by the local rank
      c2f->idx[c2_id+1] += 1;
  }

  for (cs_lnum_t i = 0; i < n_cells; i++)
    c2f->idx[i+1] += c2f->idx[i];

  const cs_lnum_t  idx_size = c2f->idx[n_cells];

  /* Fill arrays */
  BFT_MALLOC(c2f->ids, idx_size, cs_lnum_t);
  BFT_MALLOC(c2f->sgn, idx_size, short int);

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    const cs_lnum_t  c1_id = mesh->i_face_cells[f_id][0];
    const cs_lnum_t  c2_id = mesh->i_face_cells[f_id][1];

    if (c1_id < n_cells) { /* Do not consider ghost cells */

      const cs_lnum_t  shift = c2f->idx[c1_id] + cell_shift[c1_id];
      c2f->ids[shift] = f_id;
      c2f->sgn[shift] = 1;
      cell_shift[c1_id] += 1;

    }

    if (c2_id < n_cells) { /* Do not consider ghost cells */

      const cs_lnum_t  shift = c2f->idx[c2_id] + cell_shift[c2_id];
      c2f->ids[shift] = f_id;
      c2f->sgn[shift] = -1;
      cell_shift[c2_id] += 1;

    }

  } /* End of loop on internal faces */

  for (cs_lnum_t  f_id = 0; f_id < n_b_faces; f_id++) {

    const cs_lnum_t  c_id = mesh->b_face_cells[f_id];
    const cs_lnum_t  shift = c2f->idx[c_id] + cell_shift[c_id];

    c2f->ids[shift] = n_i_faces + f_id;
    c2f->sgn[shift] = 1;
    cell_shift[c_id] += 1;

  } /* End of loop on border faces */

  /* Free memory */
  BFT_FREE(cell_shift);

  return c2f;
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
  const cs_lnum_t n_vertices = m->n_vertices;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  /* Compute max. number of vertices by face */
  int  n_max_face_vertices = 0;
  for (cs_lnum_t i = 0; i < n_b_faces; i++)
    n_max_face_vertices = CS_MAX(n_max_face_vertices,
                                 m->b_face_vtx_idx[i+1] - m->b_face_vtx_idx[i]);
  for (cs_lnum_t i = 0; i < n_i_faces; i++)
    n_max_face_vertices = CS_MAX(n_max_face_vertices,
                                 m->i_face_vtx_idx[i+1] - m->i_face_vtx_idx[i]);

  cs_lnum_t  *f_vertices = NULL;
  BFT_MALLOC(f_vertices, n_max_face_vertices + 1, cs_lnum_t);

  const cs_lnum_t  n_init_edges =
    m->b_face_vtx_idx[n_b_faces] + m->i_face_vtx_idx[n_i_faces];

  if (n_init_edges == 0)
    bft_error(__FILE__, __LINE__, 0,
              " Empty mesh connectivity. Stop build edge connectivity.\n");

  /* Build e2v_lst */
  cs_gnum_t  *e2v_lst = NULL; /* Only because we have to use cs_order */
  BFT_MALLOC(e2v_lst, 2*n_init_edges, cs_gnum_t);

  cs_lnum_t  shift = 0;
  for (cs_lnum_t i = 0; i < n_b_faces; i++) {

    const cs_lnum_t  s = m->b_face_vtx_idx[i], e = m->b_face_vtx_idx[i+1];
    const int  nfv = e - s;

    for (int k = 0; k < nfv; k++)
      f_vertices[k] = m->b_face_vtx_lst[s+k] + 1;
    f_vertices[nfv] = m->b_face_vtx_lst[s] + 1;

    for (int k = 0; k < nfv; k++) {

      const cs_lnum_t  v1 = f_vertices[k], v2 = f_vertices[k+1];
      if (v1 < v2)
        e2v_lst[2*shift] = v1, e2v_lst[2*shift+1] = v2;
      else
        e2v_lst[2*shift] = v2, e2v_lst[2*shift+1] = v1;
      shift++;

    }

  } /* End of loop on border faces */

  for (cs_lnum_t i = 0; i < n_i_faces; i++) {

    const cs_lnum_t  s = m->i_face_vtx_idx[i], e = m->i_face_vtx_idx[i+1];
    const int  nfv = e - s;

    for (int k = 0; k < nfv; k++)
      f_vertices[k] = m->i_face_vtx_lst[s+k] + 1;
    f_vertices[nfv] = m->i_face_vtx_lst[s] + 1;

    for (int k = 0; k < nfv; k++) {

      const cs_lnum_t v1 = f_vertices[k], v2 = f_vertices[k+1];
      if (v1 < v2)
        e2v_lst[2*shift] = v1, e2v_lst[2*shift+1] = v2;
      else
        e2v_lst[2*shift] = v2, e2v_lst[2*shift+1] = v1;
      shift++;

    }

  } /* End of loop on interior faces */

  assert(shift == n_init_edges);

  /* Remove duplicated edges */
  cs_lnum_t  *order = NULL;
  BFT_MALLOC(order, n_init_edges, cs_lnum_t);
  cs_order_gnum_allocated_s(NULL, e2v_lst, 2, order, n_init_edges);

  cs_lnum_t  *v2v_idx = NULL;
  BFT_MALLOC(v2v_idx, n_vertices + 1, cs_lnum_t);
# pragma omp parallel for if (n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_vertices + 1; i++) v2v_idx[i] = 0;

  cs_lnum_t  *e2v_ref_lst = NULL;
  BFT_MALLOC(e2v_ref_lst, 2*n_init_edges, cs_lnum_t);

  {
    cs_lnum_t  o1 = order[0], o2 = order[0];
    cs_lnum_t  v1 = e2v_lst[2*o1],  v2 = e2v_lst[2*o1+1];

    /* Add the first edge */
    e2v_ref_lst[0] = v1;
    e2v_ref_lst[1] = v2;
    v2v_idx[v1] += 1;
    v2v_idx[v2] += 1;
    shift = 1;

    /* Add new edges which are different from the one already added */
    for (cs_lnum_t i = 1; i < n_init_edges; i++) {

      o1 = order[i-1], o2 = order[i];
      if (e2v_lst[2*o1]   != e2v_lst[2*o2] ||
          e2v_lst[2*o1+1] != e2v_lst[2*o2+1]) { /* New edge */

        v2v_idx[e2v_lst[2*o2]] += 1;
        v2v_idx[e2v_lst[2*o2+1]] += 1;
        e2v_ref_lst[2*shift] = e2v_lst[2*o2];
        e2v_ref_lst[2*shift+1] = e2v_lst[2*o2+1];
        shift++;

      }

    } /* End of loop on edges */

  }

  /* Free memory */
  BFT_FREE(e2v_lst);
  BFT_FREE(order);
  BFT_FREE(f_vertices);

  const cs_lnum_t  n_edges = shift;
  for (cs_lnum_t i = 0; i < n_vertices; i++)
    v2v_idx[i+1] += v2v_idx[i];

  assert(n_edges > 0);

  /* Allocate and set members of the edge builder structure */
  _edge_builder_t  *builder = NULL;
  BFT_MALLOC(builder, 1, _edge_builder_t);

  builder->n_vertices = n_vertices;
  builder->n_edges = n_edges;
  builder->v2v_idx = v2v_idx;
  builder->e2v_lst = e2v_ref_lst;

  cs_lnum_t  *v2v_lst = NULL, *v2e_lst = NULL;
  BFT_MALLOC(v2v_lst, v2v_idx[n_vertices], cs_lnum_t);
  BFT_MALLOC(v2e_lst, v2v_idx[n_vertices], cs_lnum_t);

  cs_lnum_t  *vtx_shift = NULL;
  BFT_MALLOC(vtx_shift, n_vertices, cs_lnum_t);
#   pragma omp parallel for if (n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_vertices; i++) vtx_shift[i] = 0;

  for (cs_lnum_t e_id = 0; e_id < n_edges; e_id++) {

    const cs_lnum_t  v1 = e2v_ref_lst[2*e_id] - 1;
    const cs_lnum_t  v2 = e2v_ref_lst[2*e_id+1] - 1;
    const cs_lnum_t  s1 = v2v_idx[v1] + vtx_shift[v1];
    const cs_lnum_t  s2 = v2v_idx[v2] + vtx_shift[v2];

    vtx_shift[v1] += 1;
    vtx_shift[v2] += 1;
    v2v_lst[s1] = v2 + 1;
    v2v_lst[s2] = v1 + 1;

    if (v1 < v2)
      v2e_lst[s1] = e_id+1, v2e_lst[s2] = -(e_id+1);
    else
      v2e_lst[s1] = -(e_id+1), v2e_lst[s2] = e_id+1;

  } /* End of loop on edges */

  /* Assign last members */
  builder->v2v_edge_lst = v2e_lst;
  builder->v2v_lst = v2v_lst;

  BFT_FREE(vtx_shift);

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
_add_f2e_entry(cs_lnum_t                  shift,
               cs_lnum_t                  v1_num,
               cs_lnum_t                  v2_num,
               const _edge_builder_t     *builder,
               cs_adjacency_t            *f2e)
{
  /* Sanity check */
  assert(v1_num > 0);
  assert(v2_num > 0);
  assert(builder->v2v_idx[v1_num] > builder->v2v_idx[v1_num-1]);

  /* Get edge number */
  cs_lnum_t  edge_sgn_num = 0;
  cs_lnum_t  *v2v_idx = builder->v2v_idx, *v2v_lst = builder->v2v_lst;

  for (cs_lnum_t i = v2v_idx[v1_num-1]; i < v2v_idx[v1_num]; i++) {
    if (v2v_lst[i] == v2_num) {
      edge_sgn_num = builder->v2v_edge_lst[i];
      break;
    }
  }

  if (edge_sgn_num == 0)
    bft_error(__FILE__, __LINE__, 0,
              _(" The given couple of vertices (number): [%d, %d]\n"
                " is not defined in the edge structure.\n"), v1_num, v2_num);

  f2e->ids[shift] = CS_ABS(edge_sgn_num) - 1;
  if (edge_sgn_num < 0)
    f2e->sgn[shift] = -1;
  else
    f2e->sgn[shift] = 1;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the face -> edges connectivity which is stored in a
 *         cs_adjacency_t structure
 *
 * \param[in]  m         pointer to a cs_mesh_t structure
 * \param[in]  builder   pointer to the _edge_builder_t structure
 *
 * \return a pointer to a new allocated cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_adjacency_t *
_build_f2e_connect(const cs_mesh_t         *m,
                   const _edge_builder_t   *builder)
{
  /* Sanity checks */
  assert(builder != NULL);

  const cs_lnum_t  n_i_faces = m->n_i_faces;
  const cs_lnum_t  n_b_faces = m->n_b_faces;
  const cs_lnum_t  n_faces = n_i_faces + n_b_faces;

  cs_adjacency_t  *f2e = NULL;
  f2e = cs_adjacency_create(CS_ADJACENCY_SIGNED, -1, n_faces);

  /* Build index */
  for (cs_lnum_t i = 0; i < n_i_faces; i++)
    f2e->idx[i+1] += m->i_face_vtx_idx[i+1] - m->i_face_vtx_idx[i];
  for (cs_lnum_t i = 0; i < n_b_faces; i++)
    f2e->idx[n_i_faces+i+1] += m->b_face_vtx_idx[i+1] - m->b_face_vtx_idx[i];
  for (cs_lnum_t i = 0; i < n_faces; i++)
    f2e->idx[i+1] += f2e->idx[i];

  assert(f2e->idx[n_faces] ==
         m->i_face_vtx_idx[n_i_faces] + m->b_face_vtx_idx[n_b_faces]);

  /* Build matrix */
  BFT_MALLOC(f2e->ids, f2e->idx[n_faces], cs_lnum_t);
  BFT_MALLOC(f2e->sgn, f2e->idx[n_faces], short int);

  /* Border faces */
  for (cs_lnum_t i = 0; i < n_b_faces; i++) {

    const cs_lnum_t  s = m->b_face_vtx_idx[i], e = m->b_face_vtx_idx[i+1];
    const cs_lnum_t  f_id = n_i_faces + i;

    cs_lnum_t  v1_num = m->b_face_vtx_lst[e-1] + 1;
    cs_lnum_t  v2_num = m->b_face_vtx_lst[s] + 1;
    cs_lnum_t  shift = f2e->idx[f_id];

    _add_f2e_entry(shift, v1_num, v2_num, builder, f2e);

    for (cs_lnum_t j = s; j < e-1; j++) {

      shift++;
      v1_num = m->b_face_vtx_lst[j] + 1;
      v2_num = m->b_face_vtx_lst[j+1] + 1;
      _add_f2e_entry(shift, v1_num, v2_num, builder, f2e);

    }

  } /* End of loop on border faces */

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    const cs_lnum_t  s = m->i_face_vtx_idx[f_id], e = m->i_face_vtx_idx[f_id+1];

    cs_lnum_t  shift = f2e->idx[f_id];
    cs_lnum_t  v1_num = m->i_face_vtx_lst[e-1] + 1;
    cs_lnum_t  v2_num = m->i_face_vtx_lst[s] + 1;

    _add_f2e_entry(shift, v1_num, v2_num, builder, f2e);

    for (cs_lnum_t j = s; j < e-1; j++) {

      shift++;
      v1_num = m->i_face_vtx_lst[j] + 1;
      v2_num = m->i_face_vtx_lst[j+1] + 1;
      _add_f2e_entry(shift, v1_num, v2_num, builder, f2e);

    }

  } /* End of loop on interior faces */

  return f2e;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the edge -> vertices connectivity which is stored in a
 *         cs_adjacency_t structure
 *
 * \param[in]  builder   pointer to the _edge_builder_t structure
 *
 * \return a pointer to a new allocated cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_adjacency_t *
_build_e2v_connect(const _edge_builder_t  *builder)
{
  const cs_lnum_t  n_edges = builder->n_edges;

  /* Arrays ids and sgn are allocated during the creation of the structure */
  cs_adjacency_t  *e2v = cs_adjacency_create(CS_ADJACENCY_SIGNED,
                                             2, // stride
                                             n_edges);

  /* Fill arrays */
# pragma omp parallel for if (n_edges > CS_THR_MIN)
  for (cs_lnum_t shift = 0; shift < 2*n_edges; shift += 2) {

    e2v->ids[shift] = builder->e2v_lst[shift] - 1;
    e2v->sgn[shift] = -1;
    e2v->ids[shift+1] = builder->e2v_lst[shift+1] - 1;
    e2v->sgn[shift+1] = 1;

  }

  return e2v;
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
  const cs_adjacency_t  *c2v = connect->c2v;
  const cs_adjacency_t  *c2e = connect->c2e;
  const cs_adjacency_t  *e2v = connect->e2v;

  assert(c2e->n_elts == n_cells);
  assert(c2v->n_elts == n_cells);

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

      const cs_lnum_t  *e2v_ids = e2v->ids + 2*c2e_ids[e];

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
    const cs_lnum_t  *c2f_ids = connect->c2f->ids + c2f_idx[0];
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

      cs_lnum_t  f_id = connect->c2f->ids[i];

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

      cs_lnum_t  f_id = connect->c2f->ids[i];

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and define a cs_range_set_t structure and a
 *        cs_interface_set_t structure for schemes with DoFs at faces.
 *
 * \param[in]       mesh          pointer to a cs_mesh_t structure
 * \param[in]       n_faces       number of faces (interior + border)
 * \param[in]       n_face_dofs   number of DoFs per face
 * \param[in, out]  p_ifs         pointer of  pointer to a cs_interface_set_t
 * \param[in, out]  p_ifs         pointer of  pointer to a cs_interface_set_t
 */
/*----------------------------------------------------------------------------*/

static void
_assign_ifs_rs(const cs_mesh_t       *mesh,
               cs_lnum_t              n_faces,
               int                    n_face_dofs,
               cs_interface_set_t   **p_ifs,
               cs_range_set_t       **p_rs)
{
  fvm_io_num_t  *i_face_io_num = NULL, *b_face_io_num = NULL;

  const cs_lnum_t  n_elts = n_faces * n_face_dofs;

  cs_gnum_t *face_gnum = NULL;
  BFT_MALLOC(face_gnum, n_elts, cs_gnum_t);

  if (cs_glob_n_ranks > 1) {

    const cs_gnum_t  *i_face_gnum;
    const cs_gnum_t  *b_face_gnum;

    if (n_face_dofs == 1) {

      i_face_gnum = mesh->global_i_face_num;
      b_face_gnum = mesh->global_b_face_num;

    }
    else {

      cs_lnum_t  *order = NULL;
      cs_gnum_t  *gnum_ordered = NULL;

      BFT_MALLOC(order, CS_MAX(mesh->n_i_faces,
                               mesh->n_b_faces), cs_lnum_t);
      BFT_MALLOC(gnum_ordered, CS_MAX(mesh->n_i_faces,
                                      mesh->n_b_faces), cs_gnum_t);

      /* Handle interior faces */
      cs_order_gnum_allocated(NULL,
                              mesh->global_i_face_num,
                              order,
                              mesh->n_i_faces);

      for (cs_lnum_t i = 0; i < mesh->n_i_faces; i++)
        gnum_ordered[i] = mesh->global_i_face_num[order[i]];

      i_face_io_num = fvm_io_num_create_from_adj_s(NULL,
                                                   gnum_ordered,
                                                   mesh->n_i_faces,
                                                   n_face_dofs);
      i_face_gnum = fvm_io_num_get_global_num(i_face_io_num);

      /* Handle border faces */
      cs_order_gnum_allocated(NULL,
                              mesh->global_b_face_num,
                              order,
                              mesh->n_b_faces);

      for (cs_lnum_t i = 0; i < mesh->n_b_faces; i++)
        gnum_ordered[i] = mesh->global_b_face_num[order[i]];

      b_face_io_num = fvm_io_num_create_from_adj_s(NULL,
                                                   gnum_ordered,
                                                   mesh->n_b_faces,
                                                   n_face_dofs);
      b_face_gnum = fvm_io_num_get_global_num(b_face_io_num);

      BFT_FREE(order);
      BFT_FREE(gnum_ordered);
    }

    cs_gnum_t  i_shift = (cs_gnum_t)mesh->n_i_faces * (cs_gnum_t)n_face_dofs;
    cs_gnum_t  global_i_shift = mesh->n_g_i_faces * (cs_gnum_t)n_face_dofs;

    assert(i_face_gnum != NULL);
#   pragma omp parallel for if (i_shift > CS_THR_MIN)
    for (cs_gnum_t i = 0; i < i_shift; i++)
      face_gnum[i] = i_face_gnum[i];

#   pragma omp parallel for if (mesh->n_b_faces*n_face_dofs > CS_THR_MIN)
    for (cs_gnum_t i = 0; i < (cs_gnum_t)(mesh->n_b_faces*n_face_dofs); i++)
      face_gnum[i + i_shift] = b_face_gnum[i] + global_i_shift;

  }
  else {

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_gnum_t i = 0; i < (cs_gnum_t)(n_elts); i++)
      face_gnum[i] = i + 1;

  } /* Sequential or parallel run */

  /* Do not consider periodicity up to now. Should split the face interface
     into interior and border faces to do this, since only boundary faces
     can be associated to a periodicity */
  cs_interface_set_t  *ifs = NULL;
  ifs = cs_interface_set_create(n_elts,
                                NULL,
                                face_gnum,
                                mesh->periodicity, 0, NULL, NULL, NULL);

  cs_range_set_t  *rs = cs_range_set_create(ifs,
                                            NULL,
                                            n_elts,
                                            false, //TODO: Ask Yvan
                                            0);    // g_id_base

  /* Free memory */
  BFT_FREE(face_gnum);
  if (n_face_dofs > 1 && cs_glob_n_ranks > 1) {
    i_face_io_num = fvm_io_num_destroy(i_face_io_num);
    b_face_io_num = fvm_io_num_destroy(b_face_io_num);
  }

  /* Return pointers */
  *p_ifs = ifs;
  *p_rs = rs;
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
  const cs_lnum_t  n_cells = mesh->n_cells;
  connect->c2f = _build_c2f_connect(mesh);

  /* Build the face --> cells connectivity */
  const cs_lnum_t  n_faces = mesh->n_i_faces + mesh->n_b_faces;
  connect->f2c = cs_adjacency_transpose(n_faces, connect->c2f);

  /* Build the face --> edges connectivity */
  _edge_builder_t  *builder = _create_edge_builder(mesh);
  connect->f2e = _build_f2e_connect(mesh, builder);

  /* Build the edge --> vertices connectivity */
  const cs_lnum_t  n_edges = builder->n_edges;
  const cs_lnum_t  n_vertices = mesh->n_vertices;
  connect->e2v = _build_e2v_connect(builder);

  _free_edge_builder(&builder);

  connect->n_vertices = n_vertices;
  connect->n_edges = n_edges;
  connect->n_faces[0] = n_faces;
  connect->n_faces[1] = mesh->n_b_faces;
  connect->n_faces[2] = mesh->n_i_faces;
  connect->n_cells = n_cells;

  /* Build additional connectivity cell --> edges and cell --> vertices
     Useful for accessing dual faces and dual volumes for instance */

  connect->c2e = cs_adjacency_compose(n_edges, connect->c2f, connect->f2e);
  cs_adjacency_sort(connect->c2e);

  connect->c2v = cs_adjacency_compose(n_vertices, connect->c2e, connect->e2v);
  cs_adjacency_sort(connect->c2v);

  /* Build the cell flag and associate a cell type to each cell */
  BFT_MALLOC(connect->cell_flag, n_cells, cs_flag_t);
  BFT_MALLOC(connect->cell_type, n_cells, fvm_element_t);

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    connect->cell_flag[c_id] = 0;
    connect->cell_type[c_id] = _get_cell_type(c_id, connect);
  }

  /* Loop on border faces and flag boundary cells */
  cs_lnum_t  *c_ids = connect->f2c->ids + connect->f2c->idx[mesh->n_i_faces];
  for (cs_lnum_t f_id = 0; f_id < mesh->n_b_faces; f_id++)
    connect->cell_flag[c_ids[f_id]] = CS_FLAG_BOUNDARY;

  /* Max number of entities (vertices, edges and faces) by cell */
  _compute_max_ent(mesh, connect);

  /* Members to handle assembly process and parallel sync. */
  connect->v_rs = NULL;

  connect->f_rs = NULL;
  connect->f_ifs = NULL;

  connect->hho1_rs = NULL;
  connect->hho2_rs = NULL;
  connect->hho1_ifs = NULL;
  connect->hho2_ifs = NULL;

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

  /* CDO face-based schemes or HHO schemes with k=0 */
  if (scheme_flag & CS_SCHEME_FLAG_CDOFB ||
      cs_test_flag(scheme_flag, CS_SCHEME_FLAG_HHO | CS_SCHEME_FLAG_POLY0))
    _assign_ifs_rs(mesh, n_faces, 1, &(connect->f_ifs), &(connect->f_rs));

  /* HHO schemes with k=1 */
  if (cs_test_flag(scheme_flag, CS_SCHEME_FLAG_HHO | CS_SCHEME_FLAG_POLY1))
    _assign_ifs_rs(mesh, n_faces, 3, &(connect->hho1_ifs), &(connect->hho1_rs));

  /* HHO schemes with k=2 */
  if (cs_test_flag(scheme_flag, CS_SCHEME_FLAG_HHO | CS_SCHEME_FLAG_POLY2))
    _assign_ifs_rs(mesh, n_faces, 6, &(connect->hho2_ifs), &(connect->hho2_rs));

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

  cs_adjacency_free(&(connect->e2v));
  cs_adjacency_free(&(connect->f2e));
  cs_adjacency_free(&(connect->f2c));

  cs_adjacency_free(&(connect->c2f));
  cs_adjacency_free(&(connect->c2e));
  cs_adjacency_free(&(connect->c2v));

  BFT_FREE(connect->cell_type);
  BFT_FREE(connect->cell_flag);

  /* Structures for parallelism */
  cs_range_set_destroy(&(connect->f_rs));
  cs_range_set_destroy(&(connect->hho1_rs));
  cs_range_set_destroy(&(connect->hho2_rs));

  cs_interface_set_destroy(&(connect->f_ifs));
  cs_interface_set_destroy(&(connect->hho1_ifs));
  cs_interface_set_destroy(&(connect->hho2_ifs));

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

  cs_adjacency_dump("Cell   --> Faces",    fdump, connect->c2f);
  cs_adjacency_dump("Face   --> Edges",    fdump, connect->f2e);
  cs_adjacency_dump("Edge   --> Vertices", fdump, connect->e2v);
  cs_adjacency_dump("Face   --> Cells",    fdump, connect->f2c);
  cs_adjacency_dump("Cell   --> Edges",    fdump, connect->c2e);
  cs_adjacency_dump("Cell   --> Vertices", fdump, connect->c2v);

  fclose(fdump);
  BFT_FREE(fname);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
