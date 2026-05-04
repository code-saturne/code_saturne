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

#include "fvm/fvm_io_num.h"

#include "bft/bft_printf.h"
#include "base/cs_math.h"

#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_adaptive_refinement.h"
#include "mesh/cs_mesh_adjacencies.h"
#include "cs_mesh_to_builder.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "mesh/cs_mesh_algorithm.h"

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Update arrays related to interior faces.
 *
 * parameters:
 *   m     <-> pointer to global mesh structure
 *   n_new <-- number of new faces
 *   f_n2o <-- new to old interior faces array
 *----------------------------------------------------------------------------*/

static void
_update_i_face_arrays(cs_mesh_t        *m,
                      cs_lnum_t         n_new,
                      const cs_lnum_t   f_n2o[])
{
  /* Allocate new arrays */

  cs_lnum_2_t *i_face_cells;
  int *i_face_family;
  char *i_face_r_gen;

  CS_MALLOC(i_face_cells, n_new, cs_lnum_2_t);
  CS_MALLOC(i_face_family, n_new, int);
  CS_MALLOC(i_face_r_gen, n_new, char);

# pragma omp parallel for if (n_new > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_new; i++) {
    cs_lnum_t j = f_n2o[i];
    /* update faces -> cells connectivity */
    i_face_cells[i][0] = m->i_face_cells[j][0];
    i_face_cells[i][1] = m->i_face_cells[j][1];
    /* update family */
    i_face_family[i] = m->i_face_family[j];
    /* update generation */
    i_face_r_gen[i] = m->i_face_r_gen[j];
  }

  CS_FREE(m->i_face_r_gen);
  CS_FREE(m->i_face_family);
  CS_FREE(m->i_face_cells);
  m->i_face_r_gen = i_face_r_gen;
  m->i_face_cells = i_face_cells;
  m->i_face_family = i_face_family;

  /* Update global numbering */

  m->n_g_i_faces
    = cs_mesh_algorithm_n2o_update_global_num(n_new, f_n2o,
                                              &(m->global_i_face_num));

  m->n_i_faces = n_new;

  /* Update connectivity */

  cs_lnum_t *i_face_vtx_idx, *i_face_vtx;
  CS_MALLOC(i_face_vtx_idx, n_new+1, cs_lnum_t);
  CS_MALLOC(i_face_vtx, m->i_face_vtx_connect_size, cs_lnum_t);

  i_face_vtx_idx[0] = 0;
  for (cs_lnum_t i = 0; i < n_new; i++) {
    cs_lnum_t j = f_n2o[i];
    cs_lnum_t dst = i_face_vtx_idx[i];
    cs_lnum_t src = m->i_face_vtx_idx[j];
    cs_lnum_t n_f_vtx = m->i_face_vtx_idx[j+1] - src;
    for (cs_lnum_t k = 0; k < n_f_vtx; k++)
      i_face_vtx[dst+k] = m->i_face_vtx_lst[src+k];
    i_face_vtx_idx[i+1] = i_face_vtx_idx[i] + n_f_vtx;
  }

  CS_FREE(m->i_face_vtx_idx);
  CS_FREE(m->i_face_vtx_lst);

  m->i_face_vtx_idx = i_face_vtx_idx;
  m->i_face_vtx_lst = i_face_vtx;

  i_face_vtx_idx = nullptr;
  i_face_vtx = nullptr;

  m->i_face_vtx_connect_size = m->i_face_vtx_idx[n_new];
}

/*----------------------------------------------------------------------------
 * Build new to old array from old to new array
 *
 * The caller is responsible for freeing the returned array.
 *
 * parameters:
 *   n_old      <-- old number of elements
 *   n_new      <-- new number of elements
 *   o2n        <-- old to new array
 *
 * returns:
 *   new to old numbering
 *----------------------------------------------------------------------------*/

static cs_lnum_t *
_build_n2o(cs_lnum_t          n_old,
           cs_lnum_t          n_new,
           const cs_lnum_t    o2n[])
{
  cs_lnum_t *n2o;
  CS_MALLOC(n2o, n_new, cs_lnum_t);
  for (cs_lnum_t i = 0; i < n_new; i++)
    n2o[i] = -1;

  for (cs_lnum_t i = 0; i < n_old; i++) {
    cs_lnum_t j = o2n[i];
    if (n2o[j] < 0)
      n2o[j] = i;
  }

  return n2o;
}

/*----------------------------------------------------------------------------*/
/*
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
_add_f2e_entry(cs_lnum_t             shift,
               cs_lnum_t             v1_id,
               cs_lnum_t             v2_id,
               const cs_adjacency_t *v2v,
               cs_adjacency_t       *f2e)
{
  /* Convention:  sgn = -1 => v2 < v1 otherwise sgn = 1
     Edge id corresponds to the position in v2v->idx */

  cs_lnum_t vidx, vref;
  if (v1_id < v2_id)
    f2e->sgn[shift] = 1, vidx = v1_id, vref = v2_id;
  else
    f2e->sgn[shift] = -1, vidx = v2_id, vref = v1_id;

#if defined(DEBUG) && !defined(NDEBUG)
  f2e->ids[shift] = -1;
#endif

  for (cs_lnum_t i = v2v->idx[vidx]; i < v2v->idx[vidx + 1]; i++) {
    if (v2v->ids[i] == vref) {
      f2e->ids[shift] = i;
      break;
    }
  }

#if defined(DEBUG) && !defined(NDEBUG)
  if (f2e->ids[shift] == -1)
    bft_error(__FILE__,
              __LINE__,
              0,
              " %s: edge not found (v1: %ld, v2: %ld)\n",
              __func__,
              (long)v1_id,
              (long)v2_id);
#endif
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update a global numbering array in case of entity renumbering
 *
 * parameters:
 *   n_new      <-- new number of elements
 *   n2o        <-- new to old array (same as old element ids list)
 *   global_num <-> global numbering (allocated if initially nullptr)
 *
 * returns:
 *   new global number of elements
 *----------------------------------------------------------------------------*/

cs_gnum_t
cs_mesh_algorithm_n2o_update_global_num(cs_lnum_t          n_new,
                                        const cs_lnum_t    n2o[],
                                        cs_gnum_t        **global_num)
{
  cs_gnum_t n_g_new = n_new;

  if (cs_glob_n_ranks == 1 && *global_num == nullptr)
    return n_g_new;

  fvm_io_num_t *n_io_num
    = fvm_io_num_create_from_select(n2o, *global_num, n_new);

  CS_FREE(*global_num);

  *global_num = fvm_io_num_transfer_global_num(n_io_num);

  n_g_new = fvm_io_num_get_global_count(n_io_num);

  n_io_num = fvm_io_num_destroy(n_io_num);

  return n_g_new;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Merge cells based on renumbering array.
 *
 * Interior faces separating merged cells are removed.
 *
 * \param[in, out]  m              mesh
 * \param[in]       n_new          new number of cells
 * \param[in]       c_o2n          cell old to new renumbering
 * \param[out]      i_f_n2o_pre    new-to-old map of interior faces induced by
 *                                 cell merge
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_algorithm_merge_cells(cs_mesh_t       *m,
                              cs_lnum_t        n_new,
                              const cs_lnum_t  c_o2n[],
                              cs_lnum_t       *i_f_n2o[])
{
  const cs_lnum_t n_old = m->n_cells;

  cs_lnum_t *c_n2o = _build_n2o(n_old, n_new, c_o2n);

  int  *cell_family;
  CS_MALLOC(cell_family, n_new, int);
  for (cs_lnum_t i = 0; i < n_new; i++) {
    cs_lnum_t j = c_n2o[i];
    cell_family[i] = m->cell_family[j];
  }
  CS_FREE(m->cell_family);
  m->cell_family = cell_family;
  cell_family = nullptr;

  int *indic = nullptr;
  CS_MALLOC(indic, n_new, int);
  for (cs_lnum_t c_id = 0; c_id < n_new; c_id++) {
    cs_lnum_t old_id = c_n2o[c_id];
    indic[c_id] = cs_glob_amr_info->indic_cells[old_id];
  }
  CS_FREE(cs_glob_amr_info->indic_cells);
  cs_glob_amr_info->indic_cells = indic;

  /* Update global numbering */

  m->n_g_cells
    = cs_mesh_algorithm_n2o_update_global_num(n_new,
                                              c_n2o,
                                              &(m->global_cell_num));

  CS_FREE(c_n2o);

  /* Transfer (cell-based) halo information to (face-based) mesh builder
     in case of periodicity, before operation modifying cell numbering  */

  cs_mesh_builder_t *mb = nullptr;

  if (m->halo != nullptr) {
    if (m->n_init_perio > 0) {
      const cs_gnum_t n_g_faces = m->n_g_i_faces + m->n_g_b_faces;
      int rank_id = cs::max(cs_glob_rank_id, 0);
      mb = cs_mesh_builder_create();
      cs_mesh_builder_define_block_dist(mb,
                                        rank_id,
                                        cs_glob_n_ranks,
                                        1,
                                        0,
                                        m->n_g_cells,
                                        n_g_faces,
                                        m->n_g_vertices);
      cs_mesh_to_builder_perio_faces(m, mb);
    }
    cs_halo_destroy(&(m->halo));
  }

  /* Update face references */

  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

# pragma omp for schedule(dynamic, CS_CL_SIZE)
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t i0 = m->i_face_cells[f_id][0];
    cs_lnum_t i1 = m->i_face_cells[f_id][1];
    if (i0 >= n_old)
      m->i_face_cells[f_id][0] = -1;
    else if (i0 > -1)
      m->i_face_cells[f_id][0] = c_o2n[i0];
    if (i1 >= n_old)
      m->i_face_cells[f_id][1] = -1;
    else if (i1 > -1)
      m->i_face_cells[f_id][1] = c_o2n[i1];
  }

# pragma omp for schedule(dynamic, CS_CL_SIZE)
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t i = m->b_face_cells[f_id];
    if (i > -1)
      m->b_face_cells[f_id] = c_o2n[i];
  }

  m->n_cells = n_new;
  m->n_cells_with_ghosts = n_new;

  /* We can now rebuild halos (and in case of periodicity, do so before
     faces are removed or merged, to convert face-based to cell-based
     information) */

  if (   m->n_domains > 1 || m->n_init_perio > 0
      || m->halo_type == CS_HALO_EXTENDED) {

    cs_mesh_init_halo(m, mb, m->halo_type, -1, false);

    if (mb != nullptr)
      cs_mesh_builder_destroy(&mb);
  }

  /* Remove excess interior faces */

  cs_lnum_t n_i_faces_new = 0;

  {
    CS_MALLOC(*i_f_n2o, m->n_i_faces, cs_lnum_t);
    cs_lnum_t *_i_f_n2o = *i_f_n2o;

    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
      cs_lnum_t i0 = m->i_face_cells[f_id][0];
      cs_lnum_t i1 = m->i_face_cells[f_id][1];
      if (i0 != i1) {
        _i_f_n2o[n_i_faces_new] = f_id;
        n_i_faces_new++;
      }
    }

    _update_i_face_arrays(m, n_i_faces_new, _i_f_n2o);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 *
 * \brief Update a global numbering array in case of entity renumbering
 *
 * parameters:
 *   n_old      <-- old number of elements
 *   n_g_old    <-- old global number of elements
 *   o2n_idx    <-- old to new index
 *   global_num <-> global numbering (allocated if initially nullptr)
 *
 * returns:
 *   new global number of elements
 */
/*----------------------------------------------------------------------------*/

cs_gnum_t
cs_mesh_algorithm_o2n_idx_update_global_num(cs_lnum_t          n_old,
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
  for (cs_lnum_t i = 0; i < n_old; i++)
    n_sub[i] = o2n_idx[i+1] - o2n_idx[i];

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

void
cs_mesh_algorithm_build_add_vertices_gnum(cs_mesh_t       *m,
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

cs_gnum_t
cs_mesh_algorithm_sync_edges_flag(const cs_mesh_t        *m,
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
 * \brief Define the face -> edges connectivity which is stored in a
 *        cs_adjacency_t structure
 *
 * \param[in] m    pointer to a cs_mesh_t structure
 * \param[in] v2v  pointer to the cs_adjacency_t structure
 *
 * \return a pointer to a new allocated cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

cs_adjacency_t *
cs_mesh_algorithm_build_f2e_connect(const cs_mesh_t      *m,
                                    const cs_adjacency_t *v2v)
{
  assert(v2v != nullptr);

  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_faces   = n_i_faces + n_b_faces;

  cs_adjacency_t *f2e = cs_adjacency_create(CS_ADJACENCY_SIGNED, -1, n_faces);

  /* Build index */

  for (cs_lnum_t i = 0; i < n_i_faces; i++)
    f2e->idx[i + 1] += m->i_face_vtx_idx[i + 1] - m->i_face_vtx_idx[i];
  for (cs_lnum_t i = 0; i < n_b_faces; i++)
    f2e->idx[n_i_faces + i + 1]
      += m->b_face_vtx_idx[i + 1] - m->b_face_vtx_idx[i];
  for (cs_lnum_t i = 0; i < n_faces; i++)
    f2e->idx[i + 1] += f2e->idx[i];

  assert(f2e->idx[n_faces]
         == m->i_face_vtx_idx[n_i_faces] + m->b_face_vtx_idx[n_b_faces]);

  /* Build matrix */

  CS_MALLOC(f2e->ids, f2e->idx[n_faces], cs_lnum_t);
  CS_MALLOC(f2e->sgn, f2e->idx[n_faces], short int);

  /* Interior faces */

# pragma omp parallel for if (n_i_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_i_faces; i++) {

    const cs_lnum_t  s       = m->i_face_vtx_idx[i];
    const int        n_vf    = m->i_face_vtx_idx[i + 1] - s;
    const cs_lnum_t *f2v_lst = m->i_face_vtx_lst + s;

    cs_lnum_t shift = f2e->idx[i];
    for (int j = 0; j < n_vf - 1; j++) {
      _add_f2e_entry(shift, f2v_lst[j], f2v_lst[j + 1], v2v, f2e);
      shift++;
    }
    _add_f2e_entry(shift, f2v_lst[n_vf - 1], f2v_lst[0], v2v, f2e);
  }

  /* Boundary faces */

# pragma omp parallel for if (n_b_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_b_faces; i++) {

    const cs_lnum_t  s       = m->b_face_vtx_idx[i];
    const int        n_vf    = m->b_face_vtx_idx[i + 1] - s;
    const cs_lnum_t *f2v_lst = m->b_face_vtx_lst + s;

    cs_lnum_t shift = f2e->idx[i + n_i_faces];
    for (int j = 0; j < n_vf - 1; j++) {
      _add_f2e_entry(shift, f2v_lst[j], f2v_lst[j + 1], v2v, f2e);
      shift++;
    }
    _add_f2e_entry(shift, f2v_lst[n_vf - 1], f2v_lst[0], v2v, f2e);

  } /* End of loop on border faces */

  return f2e;
}

/*----------------------------------------------------------------------------*/
