/*============================================================================
 * Manage connectivity (Topological features of the mesh)
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

#include "cs_array.h"
#include "cs_flag.h"
#include "cs_log.h"
#include "cs_mesh_adjacencies.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_param_types.h"
#include "cs_param_cdo.h"
#include "cs_sort.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_connect.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_cdo_connect.c

  \brief Build additional connectivities (or adjacencies) useful for building
         CDO or HHO schemes.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local macro and structure definitions
 *============================================================================*/

#define CS_CDO_CONNECT_DBG 0

static long long cs_cdo_connect_time = 0;

/* Auxiliary buffers related to local (cell-wise) problems. These buffers can
 * also be used for computing quantities related to a cs_cell_mesh_t (there are
 * as many buffers as threads since a call to these buffers can be done inside
 * an OpenMP directive) */

static int       cs_cdo_connect_cw_buffer_size = 0;
static double  **cs_cdo_connect_cw_buffer = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a entry in the face --> edges connectivity
 *
 * \param[in]      shift     position where to add the new entry
 * \param[in]      v1_id     id of the first vertex
 * \param[in]      v2_id     id of the second vertex
 * \param[in]      v2v       pointer to a cs_adjacency_t structure
 * \param[in, out] f2e       face --> edges connectivity
 */
/*----------------------------------------------------------------------------*/

static inline void
_add_f2e_entry(cs_lnum_t                  shift,
               cs_lnum_t                  v1_id,
               cs_lnum_t                  v2_id,
               const cs_adjacency_t      *v2v,
               cs_adjacency_t            *f2e)
{
  /* Convention:  sgn = -1 => v2 < v1 otherwise sgn = 1
     Edge id corresponds to the position in v2v->idx */

  cs_lnum_t  vidx, vref;
  if (v1_id < v2_id)
    f2e->sgn[shift] = 1, vidx = v1_id, vref = v2_id;
  else
    f2e->sgn[shift] = -1, vidx = v2_id, vref = v1_id;

#if defined(DEBUG) && !defined(NDEBUG)
  f2e->ids[shift] = -1;
#endif

  for (cs_lnum_t i = v2v->idx[vidx]; i < v2v->idx[vidx+1]; i++) {
    if (v2v->ids[i] == vref) {
      f2e->ids[shift] = i;
      break;
    }
  }

#if defined(DEBUG) && !defined(NDEBUG)
  if (f2e->ids[shift] == -1)
    bft_error(__FILE__, __LINE__, 0, " %s: edge not found (v1: %ld, v2: %ld)\n",
              __func__, (long)v1_id, (long)v2_id);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the face -> edges connectivity which is stored in a
 *         cs_adjacency_t structure
 *
 * \param[in]  m      pointer to a cs_mesh_t structure
 * \param[in]  v2v    pointer to the cs_adjacency_t structure
 *
 * \return a pointer to a new allocated cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_adjacency_t *
_build_f2e_connect(const cs_mesh_t         *m,
                   const cs_adjacency_t    *v2v)
{
  assert(v2v != NULL);

  const cs_lnum_t  n_i_faces = m->n_i_faces;
  const cs_lnum_t  n_b_faces = m->n_b_faces;
  const cs_lnum_t  n_faces = n_i_faces + n_b_faces;

  cs_adjacency_t  *f2e = cs_adjacency_create(CS_ADJACENCY_SIGNED, -1, n_faces);

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

  /* Interior faces */

# pragma omp parallel for if (n_i_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_i_faces; i++) {

    const cs_lnum_t  s = m->i_face_vtx_idx[i];
    const int  n_vf = m->i_face_vtx_idx[i+1] - s;
    const cs_lnum_t  *f2v_lst = m->i_face_vtx_lst + s;

    cs_lnum_t  shift = f2e->idx[i];
    for (int j = 0; j < n_vf - 1; j++) {
      _add_f2e_entry(shift, f2v_lst[j], f2v_lst[j+1], v2v, f2e);
      shift++;
    }
    _add_f2e_entry(shift, f2v_lst[n_vf-1], f2v_lst[0], v2v, f2e);

  }

  /* Boundary faces */

# pragma omp parallel for if (n_b_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_b_faces; i++) {

    const cs_lnum_t  s = m->b_face_vtx_idx[i];
    const int  n_vf = m->b_face_vtx_idx[i+1] - s;
    const cs_lnum_t  *f2v_lst = m->b_face_vtx_lst + s;

    cs_lnum_t  shift = f2e->idx[i + n_i_faces];
    for (int j = 0; j < n_vf - 1; j++) {
      _add_f2e_entry(shift, f2v_lst[j], f2v_lst[j+1], v2v, f2e);
      shift++;
    }
    _add_f2e_entry(shift, f2v_lst[n_vf-1], f2v_lst[0], v2v, f2e);

  } /* End of loop on border faces */

  return f2e;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the edge -> vertices connectivity which is stored in a
 *         cs_adjacency_t structure
 *
 * \param[in]  v2v   pointer to a cs_adjacency_t structure
 *
 * \return a pointer to a new allocated cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_adjacency_t *
_build_e2v_connect(const cs_adjacency_t  *v2v)
{
  const cs_lnum_t  n_edges = v2v->idx[v2v->n_elts];

  /* Arrays ids and sgn are allocated during the creation of the structure */

  cs_adjacency_t  *e2v = cs_adjacency_create(CS_ADJACENCY_SIGNED,
                                             2, /* stride */
                                             n_edges);

  /* Fill arrays */

# pragma omp parallel for if (n_edges > CS_THR_MIN)
  for (cs_lnum_t v1_id = 0; v1_id < v2v->n_elts; v1_id++) {
    for (cs_lnum_t j = v2v->idx[v1_id]; j < v2v->idx[v1_id+1]; j++) {

      e2v->ids[2*j] = v1_id;         /* v1_id */
      e2v->ids[2*j+1] = v2v->ids[j]; /* v2_id */
      e2v->sgn[2*j] = -1;            /* orientation v1 -> v2 */
      e2v->sgn[2*j+1] = 1;

      /* Assumption made when building a discrete Hodge operator */

      assert(e2v->ids[2*j+1] > e2v->ids[2*j]);

    }
  } /* Loop on vertices */

  return e2v;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a connectivity vertex -> vertices through cells
 *
 * \param[in]  connect       pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_adjacency_t *
_build_v2v_through_cell(const cs_cdo_connect_t     *connect)
{
  /* Build a (sorted) v2v connectivity index */

  const cs_lnum_t  n_vertices = connect->n_vertices;
  const cs_adjacency_t  *c2v = connect->c2v;

  cs_adjacency_t  *v2c = cs_adjacency_transpose(n_vertices, c2v);
  cs_adjacency_t  *v2v = cs_adjacency_compose(n_vertices, v2c, c2v);

  cs_adjacency_sort(v2v);

  /* Update index (v2v has a diagonal entry. We remove it since we have in
     mind an matrix structure stored using the MSR format (with diagonal terms
     counted outside the index) */

  cs_adjacency_remove_self_entries(v2v);

  /* Free temporary buffers */

  cs_adjacency_destroy(&v2c);

  return v2v;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a connectivity face -> faces through cells
 *
 * \param[in]  connect       pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_adjacency_t *
_build_f2f_through_cell(const cs_cdo_connect_t     *connect)
{
  cs_adjacency_t  *f2f = NULL;

  const cs_lnum_t  n_faces = connect->n_faces[CS_ALL_FACES];

  /* Build a face -> face connectivity */

  f2f = cs_adjacency_compose(n_faces, connect->f2c, connect->c2f);
  cs_adjacency_sort(f2f);

  /* Update index (f2f has a diagonal entry. We remove it since we have in
     mind an index structure for a matrix stored using the MSR format */

  cs_adjacency_remove_self_entries(f2f);

  return f2f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a connectivity edge -> edges through cells
 *
 * \param[in]  connect       pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_adjacency_t *
_build_e2e_through_cell(const cs_cdo_connect_t     *connect)
{
  const cs_lnum_t  n_edges = connect->n_edges;

  /* Build a edge -> edge connectivity */

  cs_adjacency_t  *e2c = cs_adjacency_transpose(n_edges, connect->c2e);
  assert(e2c != NULL);
  cs_adjacency_t  *e2e = cs_adjacency_compose(n_edges, e2c, connect->c2e);

  cs_adjacency_sort(e2e);

  /* Update index (e2e has a diagonal entry. We remove it since we have in
     mind an index structure for a matrix stored using the MSR format */

  cs_adjacency_remove_self_entries(e2e);

  /* Free temporary buffers */

  cs_adjacency_destroy(&e2c);

  return e2e;
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
      v_count[v_id] = 0; /* reset */

    }

    const cs_lnum_t  *c2f_idx = connect->c2f->idx + c_id;
    const cs_lnum_t  *c2f_ids = connect->c2f->ids + c2f_idx[0];
    const int  n_fc = c2f_idx[1] - c2f_idx[0];

    if (n_fc > n_max_fc) n_max_fc = n_fc;

    for (short int f = 0; f < n_fc; f++) {

      const cs_lnum_t  f_id = c2f_ids[f];
      if (f_id < m->n_i_faces) { /* Interior face */

        const cs_lnum_t  *f2v_idx = m->i_face_vtx_idx + f_id;
        const cs_lnum_t  *f2v_ids = m->i_face_vtx_lst + f2v_idx[0];
        const int  n_vf = f2v_idx[1] - f2v_idx[0];

        if (n_vf > n_max_vf) n_max_vf = n_vf;
        for (short int v = 0; v < n_vf; v++)
          v_count[f2v_ids[v]] += 1;

      }
      else { /* Border face */

        const cs_lnum_t  *f2v_idx =
          m->b_face_vtx_idx + f_id - m->n_i_faces;
        const cs_lnum_t  *f2v_ids = m->b_face_vtx_lst + f2v_idx[0];
        const int  n_vf = f2v_idx[1] - f2v_idx[0];

        if (n_vf > n_max_vf) n_max_vf = n_vf;
        for (short int v = 0; v < n_vf; v++)
          v_count[f2v_ids[v]] += 1;

      }

    } /* Loop on cell faces */

    /* Update n_max_v2fc and reset v_count */

    for (short int v = 0; v < n_vc; v++) {

      const cs_lnum_t  v_id = c2v_ids[v];
      if (v_count[v_id] > n_max_v2fc) n_max_v2fc = v_count[v_id];
      v_count[v_id] = 0; /* reset */

    }

  } /* Loop on cells */

  BFT_FREE(v_count);

  /* Store computed values */

  connect->n_max_vbyc = n_max_vc;   /* Max number of vertices for a cell */
  connect->n_max_ebyc = n_max_ec;   /* Max number of edges for a cell */
  connect->n_max_fbyc = n_max_fc;   /* Max number of faces for a cell */
  connect->n_max_v2ec = n_max_v2ec;
  connect->n_max_v2fc = n_max_v2fc;
  connect->n_max_vbyf = n_max_vf;   /* Max number of vertices in a face */
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
  fvm_element_t  ret_type = FVM_CELL_POLY; /* Default value */

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

  else if (n_vc == 6 && n_ec == 9 && n_fc == 5) { /* Could be a prism */

    int  count[2] = {0, 0};

    /* Loop on cell faces */

    for (cs_lnum_t i = connect->c2f->idx[c_id]; i < connect->c2f->idx[c_id+1];
         i++) {

      cs_lnum_t  f_id = connect->c2f->ids[i];

      if (connect->f2e->idx[f_id+1] - connect->f2e->idx[f_id] == 4) /* Quad */
        count[1] += 1;
      if (connect->f2e->idx[f_id+1] - connect->f2e->idx[f_id] == 3) /* Tria */
        count[0] += 1;

      if (count[0] == 2 && count[1] == 3)
        ret_type = FVM_CELL_PRISM;
    }

  }

  /* Hexahedron ? */

  else if (n_vc == 8 && n_ec == 12 && n_fc == 6) { /* Could be a hexahedron */

    bool  is_hexa = true;

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
 * \brief  Associate a flag to each cell.
 *
 * \param[in, out]  connect           pointer to a cs_cdo_connect_t structure
 * \param[in]       eb_scheme_flag    metadata for Edge-based schemes
 * \param[in]       vb_scheme_flag    metadata for Vertex-based schemes
 * \param[in]       vcb_scheme_flag   metadata for Vertex+Cell-based schemes
 */
/*----------------------------------------------------------------------------*/

static void
_build_cell_flag(cs_cdo_connect_t   *connect,
                 cs_flag_t           eb_scheme_flag,
                 cs_flag_t           vb_scheme_flag,
                 cs_flag_t           vcb_scheme_flag)
{
  const cs_lnum_t  n_i_faces = connect->n_faces[CS_INT_FACES];
  const cs_lnum_t  n_b_faces = connect->n_faces[CS_BND_FACES];
  const cs_lnum_t  n_cells = connect->n_cells;

  BFT_MALLOC(connect->cell_flag, n_cells, cs_flag_t);
  cs_array_flag_fill_zero(n_cells, connect->cell_flag);

  /* Loop on volume zones to find solid zones and then tag cells */

  for (int z_id = 0; z_id < cs_volume_zone_n_zones(); z_id++) {

    const cs_zone_t  *z = cs_volume_zone_by_id(z_id);
    if (z->type & CS_VOLUME_ZONE_SOLID) {
      for (cs_lnum_t i = 0; i < z->n_elts; i++)
        connect->cell_flag[z->elt_ids[i]] |= CS_FLAG_SOLID_CELL;
    }

  } /* Loop on volume zones */

  /* Loop on border faces and flag boundary cells */

  const cs_lnum_t  *c_ids = connect->f2c->ids + connect->f2c->idx[n_i_faces];
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++)
    connect->cell_flag[c_ids[f_id]] |= CS_FLAG_BOUNDARY_CELL_BY_FACE;

  if (vb_scheme_flag > 0 || vcb_scheme_flag > 0) {

    const cs_lnum_t  n_vertices = connect->n_vertices;

    cs_flag_t  *is_border_vtx = NULL;
    BFT_MALLOC(is_border_vtx, n_vertices, cs_flag_t);
    cs_array_flag_fill_zero(n_vertices, is_border_vtx);

    const cs_adjacency_t  *bf2v = connect->bf2v;
    for (cs_lnum_t bf_id = 0; bf_id < n_b_faces; bf_id++) {
      for (cs_lnum_t j = bf2v->idx[bf_id]; j < bf2v->idx[bf_id+1]; j++) {
        assert(bf2v->ids[j] < n_vertices);
        is_border_vtx[bf2v->ids[j]] = 1;
      }
    } /* Loop on border faces */

    /* Synchronization needed in parallel or periodic computations */

    if (connect->vtx_ifs != NULL)
      cs_interface_set_max(connect->vtx_ifs,
                           n_vertices,
                           1,             /* stride */
                           false,         /* interlace (not useful here) */
                           CS_FLAG_TYPE,  /* unsigned short int */
                           is_border_vtx);

    const cs_adjacency_t  *c2v = connect->c2v;
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {
        if (is_border_vtx[c2v->ids[j]] > 0)
          connect->cell_flag[c_id] |= CS_FLAG_BOUNDARY_CELL_BY_VERTEX;
      }
    } /* Loop on cells */

    BFT_FREE(is_border_vtx);

  } /* At least one equation with DoFs at vertices */

  if (eb_scheme_flag > 0) {

    const cs_lnum_t  n_edges = connect->n_edges;

    cs_flag_t  *is_border_edge = NULL;
    BFT_MALLOC(is_border_edge, n_edges, cs_flag_t);
    cs_array_flag_fill_zero(n_edges, is_border_edge);

    const cs_adjacency_t  *f2e = connect->f2e;
    for (cs_lnum_t bf_id = n_i_faces; bf_id < n_i_faces + n_b_faces; bf_id++) {
      for (cs_lnum_t j = f2e->idx[bf_id]; j < f2e->idx[bf_id+1]; j++) {
        is_border_edge[f2e->ids[j]] = 1;
      }
    } /* Loop on border faces */

    /* Synchronization needed in parallel or periodic computations */

    if (connect->edge_ifs != NULL)
      cs_interface_set_max(connect->edge_ifs,
                           n_edges,
                           1,             /* stride */
                           false,         /* interlace (not useful here) */
                           CS_FLAG_TYPE,  /* unsigned short int */
                           is_border_edge);

    const cs_adjacency_t  *c2e = connect->c2e;
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      for (cs_lnum_t j = c2e->idx[c_id]; j < c2e->idx[c_id+1]; j++) {
        if (is_border_edge[c2e->ids[j]] > 0)
          connect->cell_flag[c_id] |= CS_FLAG_BOUNDARY_CELL_BY_EDGE;
      }
    } /* Loop on cells */

    BFT_FREE(is_border_edge);

  } /* edge interfaces */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and define a \ref cs_range_set_t structure and a
 *        \ref cs_interface_set_t structure for schemes with DoFs at edges
 *        otherwise set the number of global edges
 *
 * \param[in]       mesh            pointer to a cs_mesh_t structure
 * \param[in]       connect         pointer to a cs_cdo_connect_t structure
 * \param[in]       eb_scheme_flag  metadata for edge-based schemes
 * \param[in, out]  p_ifs           pointer of pointer to a cs_interface_set_t
 * \param[in, out]  p_rs            pointer of pointer to a cs_range_set_t
 */
/*----------------------------------------------------------------------------*/

static void
_assign_edge_ifs_rs(const cs_mesh_t       *mesh,
                    cs_cdo_connect_t      *connect,
                    cs_flag_t              eb_scheme_flag,
                    cs_interface_set_t   **p_ifs,
                    cs_range_set_t       **p_rs)
{
  const cs_lnum_t  n_edges = connect->n_edges;
  cs_gnum_t  n_g_edges = n_edges;
  cs_gnum_t  *edge_gnum = NULL;

  BFT_MALLOC(edge_gnum, n_edges, cs_gnum_t);

  if (cs_glob_n_ranks > 1) {

    const cs_adjacency_t  *e2v = connect->e2v;

    /* Build global edge numbering and edges interface */

    cs_gnum_t *e2v_gnum = NULL;
    BFT_MALLOC(e2v_gnum, n_edges*2, cs_gnum_t);

#   pragma omp parallel for if (n_edges > CS_THR_MIN)
    for (cs_lnum_t e_id = 0; e_id < n_edges; e_id++) {

      cs_gnum_t  *v_gids = e2v_gnum + 2*e_id;
      const cs_lnum_t  *_v_ids = e2v->ids + 2*e_id;
      const cs_gnum_t  v0_gid = mesh->global_vtx_num[_v_ids[0]];
      const cs_gnum_t  v1_gid = mesh->global_vtx_num[_v_ids[1]];

      if (v0_gid < v1_gid)
        v_gids[0] = v0_gid, v_gids[1] = v1_gid;
      else
        v_gids[0] = v1_gid, v_gids[1] = v0_gid;

    } /* Loop on edges */

    cs_lnum_t  *order = NULL;
    BFT_MALLOC(order, n_edges, cs_lnum_t);
    cs_order_gnum_allocated_s(NULL, e2v_gnum, 2, order, n_edges);

    cs_gnum_t  *order_couples = NULL;
    BFT_MALLOC(order_couples, 2*n_edges, cs_gnum_t);
#   pragma omp parallel for if (n_edges > CS_THR_MIN)
    for (cs_lnum_t e = 0; e < n_edges; e++) {
      const cs_lnum_t  o_id = 2*order[e];
      order_couples[2*e] = e2v_gnum[o_id];
      order_couples[2*e+1] = e2v_gnum[o_id+1];
    }

    fvm_io_num_t *edge_io_num
      = fvm_io_num_create_from_adj_s(NULL, order_couples, n_edges, 2);

    BFT_FREE(order);
    BFT_FREE(order_couples);
    BFT_FREE(e2v_gnum);

    const cs_gnum_t *_g_num =  fvm_io_num_get_global_num(edge_io_num);
    memcpy(edge_gnum, _g_num, n_edges*sizeof(cs_gnum_t));

    connect->n_g_edges = fvm_io_num_get_global_count(edge_io_num);
    edge_io_num = fvm_io_num_destroy(edge_io_num);

  }
  else {

#   pragma omp parallel for if (n_edges > CS_THR_MIN)
    for (cs_gnum_t i = 0; i < n_g_edges; i++) edge_gnum[i] = i + 1;

  } /* Sequential or parallel run */

  cs_interface_set_t  *ifs = NULL;
  cs_range_set_t  *rs = NULL;

  if (eb_scheme_flag > 0) {

  /* Do not consider periodicity up to now. Should split the face interface
     into interior and border faces to do this, since only boundary faces
     can be associated to a periodicity */

    ifs = cs_interface_set_create(n_edges,
                                  NULL,
                                  edge_gnum,
                                  mesh->periodicity,
                                  0,
                                  NULL, NULL, NULL);

    rs = cs_range_set_create(ifs,   /* interface set */
                             NULL,  /* halo */
                             n_edges,
                             false, /* TODO: add option for balance */
                             2,     /* tr_ignore */
                             0);    /* g_id_base */

  }

  /* Free memory */

  BFT_FREE(edge_gnum);

  /* Return pointers */

  *p_ifs = ifs;
  *p_rs = rs;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and define a new cs_interface_set_t structure on faces
 *
 * \param[in]  mesh          pointer to a cs_mesh_t structure
 *
 * \return a pointer to a new allocated cs_interface_set_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_interface_set_t *
_define_face_interface(const cs_mesh_t  *mesh)
{
  cs_interface_set_t  *ifs = NULL;

  const cs_gnum_t *face_gnum = mesh->global_i_face_num;
  cs_gnum_t *_face_gnum = NULL;

  if (face_gnum == NULL) {
    cs_lnum_t _n_faces = mesh->n_i_faces;
    BFT_MALLOC(_face_gnum, _n_faces, cs_gnum_t);

#   pragma omp parallel for if (_n_faces > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < _n_faces; i++)
      _face_gnum[i] = i+1;

    face_gnum = _face_gnum;
  }

  int n_perio = mesh->n_init_perio;
  int *perio_num = NULL;
  cs_lnum_t *n_perio_face_couples = NULL;
  cs_gnum_t **perio_face_couples = NULL;

  if (n_perio > 0) {
    BFT_MALLOC(perio_num, n_perio, int);
    for (int i = 0; i < n_perio; i++)
      perio_num[i] = i+1;
  }

  cs_mesh_get_perio_faces(mesh,
                          &n_perio_face_couples,
                          &perio_face_couples);

  ifs = cs_interface_set_create(mesh->n_i_faces,
                                NULL,
                                face_gnum,
                                mesh->periodicity,
                                n_perio,
                                perio_num,
                                n_perio_face_couples,
                                (const cs_gnum_t *const *)perio_face_couples);

  if (n_perio > 0) {
    for (int i = 0; i < n_perio; i++)
      BFT_FREE(perio_face_couples[i]);
    BFT_FREE(perio_face_couples);
    BFT_FREE(n_perio_face_couples);
    BFT_FREE(perio_num);
  }

  BFT_FREE(_face_gnum);

  return ifs;
}

/*! \endcond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and define a new cs_cdo_connect_t structure
 *        Range sets and interface sets are allocated and defined according to
 *        the value of the different scheme flags.
 *        cs_range_set_t structure related to vertices is shared the cs_mesh_t
 *        structure (the global one)
 *
 * \param[in, out] mesh              pointer to a cs_mesh_t structure
 * \param[in]      eb_scheme_flag    metadata for Edge-based schemes
 * \param[in]      fb_scheme_flag    metadata for Face-based schemes
 * \param[in]      cb_scheme_flag    metadata for Cell-based schemes
 * \param[in]      vb_scheme_flag    metadata for Vertex-based schemes
 * \param[in]      vcb_scheme_flag   metadata for Vertex+Cell-based schemes
 * \param[in]      hho_scheme_flag   metadata for HHO schemes
 *
 * \return a pointer to a cs_cdo_connect_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_connect_t *
cs_cdo_connect_build(cs_mesh_t      *mesh,
                     cs_flag_t       eb_scheme_flag,
                     cs_flag_t       fb_scheme_flag,
                     cs_flag_t       cb_scheme_flag,
                     cs_flag_t       vb_scheme_flag,
                     cs_flag_t       vcb_scheme_flag,
                     cs_flag_t       hho_scheme_flag)
{
  cs_timer_t t0 = cs_timer_time();

  cs_cdo_connect_t  *connect = NULL;

  const cs_lnum_t  n_vertices = mesh->n_vertices;
  const cs_lnum_t  n_faces = mesh->n_i_faces + mesh->n_b_faces;
  const cs_lnum_t  n_cells = mesh->n_cells;

  /* Build the connectivity structure */

  BFT_MALLOC(connect, 1, cs_cdo_connect_t);

  /* Map the boundary face --> vertices connectivity */

  connect->bf2v = cs_adjacency_create_from_i_arrays(mesh->n_b_faces,
                                                    mesh->b_face_vtx_idx,
                                                    mesh->b_face_vtx_lst,
                                                    NULL);

  /* Map the interior face --> vertices connectivity */

  connect->if2v = cs_adjacency_create_from_i_arrays(mesh->n_i_faces,
                                                    mesh->i_face_vtx_idx,
                                                    mesh->i_face_vtx_lst,
                                                    NULL);

  /* Build the cell --> faces connectivity */

  connect->c2f = cs_mesh_adjacency_c2f(mesh, 1);

  /* Build the face --> cells connectivity */

  connect->f2c = cs_adjacency_transpose(n_faces, connect->c2f);

  /* Build the face --> edges connectivity */

  cs_adjacency_t  *v2v_e = cs_mesh_adjacency_v2v(mesh);
  const cs_lnum_t  n_edges = v2v_e->idx[n_vertices];
  connect->f2e = _build_f2e_connect(mesh, v2v_e);

  /* Build the edge --> vertices connectivity */

  connect->e2v = _build_e2v_connect(v2v_e);

  cs_adjacency_destroy(&v2v_e);

  connect->n_vertices = n_vertices;
  connect->n_edges = n_edges;
  connect->n_g_edges = n_edges;
  connect->n_faces[CS_ALL_FACES] = n_faces;
  connect->n_faces[CS_BND_FACES] = mesh->n_b_faces;
  connect->n_faces[CS_INT_FACES] = mesh->n_i_faces;
  connect->n_cells = n_cells;
  connect->n_cells_with_ghosts = mesh->n_cells_with_ghosts;

  /* Build additional connectivity cell --> edges and cell --> vertices
     Useful for accessing dual faces and dual volumes for instance */

  connect->c2e = cs_adjacency_compose(n_edges, connect->c2f, connect->f2e);
  cs_adjacency_sort(connect->c2e);

  connect->c2v = cs_adjacency_compose(n_vertices, connect->c2e, connect->e2v);
  cs_adjacency_sort(connect->c2v);

  /* Max number of entities (vertices, edges and faces) by cell */

  _compute_max_ent(mesh, connect);

  /* Build the adjacency needed to define the matrix structure of the
     linear systems either for Vertex-based or Face-based schemes */

  if (vb_scheme_flag > 0 || vcb_scheme_flag > 0)
    connect->v2v = _build_v2v_through_cell(connect);
  else
    connect->v2v = NULL;

  if (fb_scheme_flag > 0 || cb_scheme_flag > 0 || hho_scheme_flag > 0)
    connect->f2f = _build_f2f_through_cell(connect);
  else
    connect->f2f = NULL;

  if (eb_scheme_flag > 0)
    connect->e2e = _build_e2e_through_cell(connect);
  else
    connect->e2e = NULL;

  /* Members to handle assembly process and parallel sync. */

  connect->vtx_rset = NULL;
  connect->vtx_ifs = mesh->vtx_interfaces;  /* Already defined */

  connect->edge_rset = NULL;
  connect->edge_ifs = NULL;

  connect->face_rset = NULL;
  connect->face_ifs = NULL;

  /* At least one CDO vertex- or vertex+cell-based schemes has been requested */

  if ( vb_scheme_flag > 0 || vcb_scheme_flag > 0) {

    connect->vtx_rset = cs_range_set_create(connect->vtx_ifs,
                                            NULL,
                                            n_vertices,
                                            false,  /* option for balance */
                                            1,      /* tr_ignore */
                                            0);     /* g_id_base */

    if (mesh->vtx_range_set != NULL)
      cs_range_set_destroy(&(mesh->vtx_range_set));

    mesh->vtx_range_set = connect->vtx_rset;

  }

  /* CDO face-based schemes or HHO schemes with k=0 */

  if (fb_scheme_flag > 0 || cb_scheme_flag > 0 || hho_scheme_flag > 0) {

    connect->face_ifs = _define_face_interface(mesh);

    connect->face_rset = cs_range_set_create(connect->face_ifs,
                                             NULL,   /* halo */
                                             n_faces,
                                             false,  /* option for balance */
                                             1,      /* tr_ignore */
                                             0);     /* g_id_base */

  }

  /* CDO edge-based schemes for vector-valued variables if needed (DoF is
     scalar-valued) */

  _assign_edge_ifs_rs(mesh, connect, eb_scheme_flag,
                      &(connect->edge_ifs),
                      &(connect->edge_rset));

  /* Build the cell type for each cell */

  BFT_MALLOC(connect->cell_type, n_cells, fvm_element_t);
# pragma omp parallel if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    connect->cell_type[c_id] = _get_cell_type(c_id, connect);

  /* Build the cell flag for each cell */

  _build_cell_flag(connect, eb_scheme_flag, vb_scheme_flag, vcb_scheme_flag);

  /* Monitoring */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_t  time_count = cs_timer_diff(&t0, &t1);
  cs_cdo_connect_time += time_count.nsec;

  return connect;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a cs_cdo_connect_t structure
 *
 * \param[in]      mesh        pointer to a mesh structure
 * \param[in, out] connect     pointer to additional CDO connectivities
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_cdo_connect_t *
cs_cdo_connect_free(const cs_mesh_t    *mesh,
                    cs_cdo_connect_t   *connect)
{
  if (connect == NULL)
    return connect;

  cs_adjacency_destroy(&(connect->bf2v));
  cs_adjacency_destroy(&(connect->if2v));

  cs_adjacency_destroy(&(connect->e2v));
  cs_adjacency_destroy(&(connect->f2e));
  cs_adjacency_destroy(&(connect->f2c));

  cs_adjacency_destroy(&(connect->c2f));
  cs_adjacency_destroy(&(connect->c2e));
  cs_adjacency_destroy(&(connect->c2v));

  cs_adjacency_destroy(&(connect->v2v));
  cs_adjacency_destroy(&(connect->f2f));
  cs_adjacency_destroy(&(connect->e2e));

  BFT_FREE(connect->cell_type);
  BFT_FREE(connect->cell_flag);

  /* Structures for parallelism */

  if (mesh->vtx_range_set != connect->vtx_rset)
    cs_range_set_destroy(&(connect->vtx_rset));
  else
    connect->vtx_rset = NULL;

  cs_range_set_destroy(&(connect->edge_rset));
  cs_range_set_destroy(&(connect->face_rset));

  /* Done in cs_mesh.c since it is shared with the mesh structure :
     cs_interface_set_destroy(&(connect->vtx_ifs)); */

  cs_interface_set_destroy(&(connect->edge_ifs));
  cs_interface_set_destroy(&(connect->face_ifs));

  BFT_FREE(connect);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Summary of the connectivity information
 *
 * \param[in] connect           pointer to cs_cdo_connect_t structure
 * \param[in] eb_scheme_flag    metadata for Edge-based schemes
 * \param[in] vb_scheme_flag    metadata for Vertex-based schemes
 * \param[in] vcb_scheme_flag   metadata for Vertex+Cell-based schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_connect_log_summary(const cs_cdo_connect_t  *connect,
                           cs_flag_t                eb_scheme_flag,
                           cs_flag_t                vb_scheme_flag,
                           cs_flag_t                vcb_scheme_flag)
{
  /* Information about the element types */

  cs_gnum_t  n_type_cells[FVM_N_ELEMENT_TYPES];
  for (int i = 0; i < FVM_N_ELEMENT_TYPES; i++)
    n_type_cells[i] = 0;

  for (cs_lnum_t i = 0; i < connect->n_cells; i++)
    n_type_cells[connect->cell_type[i]] += 1;

  cs_parall_sum(FVM_N_ELEMENT_TYPES, CS_GNUM_TYPE, n_type_cells);

  cs_log_printf(CS_LOG_DEFAULT, "\n Connectivity information:\n");
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


  cs_lnum_t  n_max_entbyc[5] = {connect->n_max_fbyc,
                                connect->n_max_ebyc,
                                connect->n_max_vbyc,
                                connect->v_max_cell_range,
                                connect->e_max_cell_range};

  cs_parall_max(5, CS_LNUM_TYPE, n_max_entbyc);

  /* Output */

  cs_log_printf(CS_LOG_DEFAULT,
                " --dim-- max. number of faces by cell:    %4d\n",
                (int)n_max_entbyc[0]);
  cs_log_printf(CS_LOG_DEFAULT,
                " --dim-- max. number of edges by cell:    %4d\n",
                (int)n_max_entbyc[1]);
  cs_log_printf(CS_LOG_DEFAULT,
                " --dim-- max. number of vertices by cell: %4d\n",
                (int)n_max_entbyc[2]);
  cs_log_printf(CS_LOG_DEFAULT,
                "\n --dim-- max. vertex range for a cell:    %ld\n",
                (long)n_max_entbyc[3]);
  cs_log_printf(CS_LOG_DEFAULT,
                " --dim-- max. edge range for a cell:      %ld\n",
                (long)n_max_entbyc[4]);

  /* Information about special case where vertices are lying on the boundary
     but not a face (for instance a tetrahedron) */

  if (vb_scheme_flag > 0 || vcb_scheme_flag > 0) {

    cs_lnum_t  v_count = 0;
    for (cs_lnum_t c = 0; c < connect->n_cells; c++) {
      if ((connect->cell_flag[c] & CS_FLAG_BOUNDARY_CELL_BY_VERTEX) &&
          !(connect->cell_flag[c] & CS_FLAG_BOUNDARY_CELL_BY_FACE))
        v_count++;
    }

    cs_gnum_t  v_counter = v_count;
    cs_parall_counter(&v_counter, 1);

    cs_log_printf(CS_LOG_DEFAULT,
                  " --dim-- number of boundary cells through a vertex only"
                  " %lu\n\n", v_counter);

  } /* At least one equation with a scheme with DoFs at vertices */

  if (eb_scheme_flag > 0) {

    /* Information about special case where edges are lying on the boundary
       but not a face (for instance a tetrahedron) */

    cs_lnum_t  e_count = 0;
    for (cs_lnum_t c = 0; c < connect->n_cells; c++) {
      if ((connect->cell_flag[c] & CS_FLAG_BOUNDARY_CELL_BY_EDGE) &&
          !(connect->cell_flag[c] & CS_FLAG_BOUNDARY_CELL_BY_FACE))
        e_count++;
    }

    cs_gnum_t  e_counter = e_count;
    cs_parall_counter(&e_counter, 1);

    cs_log_printf(CS_LOG_DEFAULT,
                  " --dim-- number of boundary cells through an edge only"
                  " %lu\n\n", e_counter);

  } /* At least one equation with a scheme with DoFs at edges */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and initialize the cell-wise buffer(s)
 *
 * \param[in] connect     pointer to additional connectivities
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_connect_allocate_cw_buffer(const cs_cdo_connect_t      *connect)
{
  int  n_vc = connect->n_max_vbyc;
  int  max_ent = 3*CS_MAX(n_vc, CS_MAX(connect->n_max_ebyc,
                                       connect->n_max_fbyc));

  cs_cdo_connect_cw_buffer_size = CS_MAX(n_vc*(n_vc+1)/2, max_ent);

  BFT_MALLOC(cs_cdo_connect_cw_buffer, cs_glob_n_threads, double *);

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    BFT_MALLOC(cs_cdo_connect_cw_buffer[t_id],
               cs_cdo_connect_cw_buffer_size, double);
    memset(cs_cdo_connect_cw_buffer[t_id], 0,
           cs_cdo_connect_cw_buffer_size*sizeof(double));
  }
#else
  assert(cs_glob_n_threads == 1);

  BFT_MALLOC(cs_cdo_connect_cw_buffer[0],
             cs_cdo_connect_cw_buffer_size, double);
  memset(cs_cdo_connect_cw_buffer[0], 0,
         cs_cdo_connect_cw_buffer_size*sizeof(double));
#endif /* openMP */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the cell-wise buffer(s)
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_connect_free_cw_buffer(void)
{
  cs_cdo_connect_cw_buffer_size = 0;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    BFT_FREE(cs_cdo_connect_cw_buffer[t_id]);
  }
#else
  assert(cs_glob_n_threads == 1);
  BFT_FREE(cs_cdo_connect_cw_buffer[0]);
#endif /* openMP */

  BFT_FREE(cs_cdo_connect_cw_buffer);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the cell-wise buffer associated to the current thread. Use
 *        \ref cs_get_thread_id when the thread id is not directly available.
 *
 * \param[in] thr_id   thread_id (0 if no openMP)
 *
 * \return a pointer to an allocated auxiliary cell-wise buffer
 */
/*----------------------------------------------------------------------------*/

double *
cs_cdo_connect_get_cw_buffer(int   thr_id)
{
  assert(thr_id > -1 && thr_id < cs_glob_n_threads);
  return cs_cdo_connect_cw_buffer[thr_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the time elapsed to build the cs_cdo_connect_t structure
 *
 * \return the value of the time elapsed in ns
 */
/*----------------------------------------------------------------------------*/

long long
cs_cdo_connect_get_time_perfo(void)
{
  return cs_cdo_connect_time;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and define a \ref cs_range_set_t structure and a
 *        \ref cs_interface_set_t structure for schemes with DoFs at vertices
 *
 * \param[in]       mesh          pointer to a cs_mesh_t structure
 * \param[in]       n_vtx_dofs    number of DoFs per vertex
 * \param[in]       interlaced    false means a block viewpoint
 * \param[in, out]  p_ifs         pointer of  pointer to a cs_interface_set_t
 * \param[in, out]  p_rs          pointer of  pointer to a cs_range_set_t
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_connect_assign_vtx_ifs_rs(const cs_mesh_t       *mesh,
                                 int                    n_vtx_dofs,
                                 bool                   interlaced,
                                 cs_interface_set_t   **p_ifs,
                                 cs_range_set_t       **p_rs)
{
  assert(mesh != NULL);

  const cs_lnum_t  n_vertices = mesh->n_vertices;
  const cs_lnum_t  n_elts = n_vertices * n_vtx_dofs;

  /* Structure to define */

  cs_range_set_t  *rs =  NULL;
  cs_interface_set_t  *ifs = NULL;

  switch (n_vtx_dofs) {

  case 1: /* Scalar-valued */
    ifs = *p_ifs;
    assert(ifs != NULL || cs_glob_n_ranks < 2);  /* Should be already set */
    break;

  default:
    if (interlaced)
      ifs = cs_interface_set_dup(mesh->vtx_interfaces, n_vtx_dofs);
    else
      ifs = cs_interface_set_dup_blocks(mesh->vtx_interfaces,
                                        n_vertices,
                                        n_vtx_dofs);
    break;

  } /* End of switch on dimension */

  rs = cs_range_set_create(ifs,
                           NULL,
                           n_elts,
                           false,  /* TODO: add option for balance */
                           1,      /* tr_ignore */
                           0);     /* g_id_base */

  /* Return pointers */

  *p_ifs = ifs;
  *p_rs = rs;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the discrete curl operator across each primal faces.
 *        From an edge-based array (seen as circulations) compute a face-based
 *        array (seen as fluxes)
 *
 * \param[in]      connect        pointer to a cs_cdo_connect_t struct.
 * \param[in]      edge_values    array of values at edges
 * \param[in, out] p_curl_values  array storing the curl across faces. This
 *                                array is allocated if necessary.
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_connect_discrete_curl(const cs_cdo_connect_t    *connect,
                             const cs_real_t           *edge_values,
                             cs_real_t                **p_curl_values)
{
  if (connect == NULL || edge_values == NULL)
    return;
  const cs_lnum_t  n_faces = connect->n_faces[CS_ALL_FACES];

  cs_real_t  *curl_values = *p_curl_values;
  if (curl_values == NULL)
    BFT_MALLOC(curl_values, n_faces, cs_real_t);

  const cs_adjacency_t  *f2e = connect->f2e;
  assert(f2e != NULL && f2e->sgn != NULL); /* Sanity checks */

# pragma omp parallel for if (n_faces > CS_THR_MIN)
  for (cs_lnum_t f = 0; f < n_faces; f++) {

    const cs_lnum_t  start = f2e->idx[f], end = f2e->idx[f+1];
    const cs_lnum_t  *ids = f2e->ids + start;
    const short int  *sgn = f2e->sgn + start;

    curl_values[f] = 0;
    for (int e = 0; e < end-start; e++)
      curl_values[f] += sgn[e]*edge_values[ids[e]];

  } /* Loop on faces */

  /* Return pointer */

  *p_curl_values = curl_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Dump a cs_cdo_connect_t structure for debugging purpose
 *
 * \param[in] connect     pointer to cs_cdo_connect_t structure
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
  cs_adjacency_dump("Bd Face--> Vertices", fdump, connect->bf2v);
  cs_adjacency_dump("In Face--> Vertices", fdump, connect->if2v);
  cs_adjacency_dump("Edge   --> Vertices", fdump, connect->e2v);
  cs_adjacency_dump("Face   --> Cells",    fdump, connect->f2c);
  cs_adjacency_dump("Cell   --> Edges",    fdump, connect->c2e);
  cs_adjacency_dump("Cell   --> Vertices", fdump, connect->c2v);

  fclose(fdump);
  BFT_FREE(fname);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
