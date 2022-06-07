/*============================================================================
 * Mesh coarsening.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
#include "cs_mesh_adjacencies.h"
#include "cs_mesh_builder.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_to_builder.h"
#include "cs_order.h"
#include "cs_sort.h"
#include "cs_parall.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_coarsen.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_mesh_coarsen.c
        Mesh coarsening.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * State of current face merging (working structure)
 *----------------------------------------------------------------------------*/

typedef struct {

  cs_lnum_t    *face_vertices;     /* current triangle vertices list */
  cs_lnum_2_t  *e2v;               /* edge to vertices */
  cs_lnum_t    *e2f;               /* edge to face (first face) */

  cs_lnum_t     n_edges;           /* number of edges */
  cs_lnum_t     n_edges_max;       /* Maximum number of edges */
  cs_lnum_t     n_vertices;        /* number of edges */
  cs_lnum_t     n_vertices_max;    /* Maximum number of vertices */

} cs_mesh_face_merge_state_t;

/*============================================================================
 * Private function definitions
 *============================================================================*/

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
 * \brief Build cells equivalence id array.
 *
 * Cells can only be merged when all cells they should be merged with
 * are also flagged for merging (isotropic merging).
 *
 * The caller is responsible for freeing the returned array.
 *
 * \param[in, out]  m          mesh
 * \param[in]       cell_flag  coarsening flag for each cell (0: no 1: yes)
 * \param[out]      c_o2n      cell old to new renumbering
 *
 * \return  number of new cells
 */
/*----------------------------------------------------------------------------*/

static cs_lnum_t
_cell_equiv(cs_mesh_t  *m,
            const int   cell_flag[],
            cs_lnum_t  *c_o2n[])
{
  cs_lnum_t  *c_equiv;
  char       *c_r_level;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;

  BFT_MALLOC(c_equiv, n_cells_ext, cs_lnum_t);
  BFT_MALLOC(c_r_level, n_cells_ext, char);

  for (cs_lnum_t i = 0; i < n_cells_ext; i++) {
    c_equiv[i] = i;
    c_r_level[i] = 0;
  }

  /* Assign refinement level based on highest generation face */

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    if (m->i_face_r_gen[f_id] < 1)
      continue;
    for (cs_lnum_t i = 0; i < 2; i++) {
      cs_lnum_t c_id = m->i_face_cells[f_id][i];
      if (c_id < n_cells) {
        if (m->i_face_r_gen[f_id] > c_r_level[c_id])
          c_r_level[c_id] = m->i_face_r_gen[f_id];
      }
    }
  }

  if (m->halo != NULL)
    cs_halo_sync_untyped(m->halo,
                         CS_HALO_STANDARD,
                         1,
                         c_r_level);

  /* Now determine cells built from the same parent */
  /* TODO handle parallelism, with parent split across
     multiple ranks */

  int reloop = 0;

  do {

    reloop = 0;

    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
      if (m->i_face_r_gen[f_id] < 1)
        continue;
      cs_lnum_t c_id0 = m->i_face_cells[f_id][0];
      cs_lnum_t c_id1 = m->i_face_cells[f_id][1];
      if (   m->i_face_r_gen[f_id] == c_r_level[c_id0]
          && m->i_face_r_gen[f_id] == c_r_level[c_id1]) {
        cs_lnum_t min_equiv = CS_MIN(c_equiv[c_id0], c_equiv[c_id1]);
        if (c_equiv[c_id0] > min_equiv) {
          c_equiv[c_id0] = min_equiv;
          reloop = 1;
        }
        if (c_equiv[c_id1] > min_equiv) {
          c_equiv[c_id1] = min_equiv;
          reloop = 1;
        }
      }
    }

  } while (reloop);

  /* Now determine whether all subcells of a given parent are flagged
     for merging; otherwise do not merge */

  for (cs_lnum_t i = 0; i < n_cells; i++) {
    cs_lnum_t j = c_equiv[i];
    if (cell_flag[i] == 0)
      c_r_level[j] = 0;
  }
  for (cs_lnum_t i = 0; i < n_cells; i++) {
    cs_lnum_t j = c_equiv[i];
    if (j != i) {
      if (c_r_level[j] == 0) {
        c_r_level[i] = 0;
        c_equiv[i] = i;
      }
    }
  }

  BFT_FREE(c_r_level);

  /* Now compact renumbering array */

  cs_lnum_t n_c_new = 0;
  for (cs_lnum_t i = 0; i < n_cells; i++) {
    if (c_equiv[i] == i) {
      c_equiv[i] = n_c_new;
      n_c_new += 1;
    }
    else {
      assert(c_equiv[i] < i);
      c_equiv[i] = c_equiv[c_equiv[i]];
    }
  }

  *c_o2n = c_equiv;

  return n_c_new;
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
  BFT_MALLOC(n2o, n_new, cs_lnum_t);
  for (cs_lnum_t i = 0; i < n_new; i++)
    n2o[i] = -1;

  for (cs_lnum_t i = 0; i < n_old; i++) {
    cs_lnum_t j = o2n[i];
    if (n2o[j] < 0)
      n2o[j] = i;
  }

  return n2o;
}

/*----------------------------------------------------------------------------
 * Build indexed mapping from new to old array.
 *
 * The caller is responsible for freeing the returned array.
 *
 * parameters:
 *   n_old      <-- old number of elements
 *   n_new      <-- new number of elements
 *   o2n        <-- old to new array
 *   n2o_idx    --> new to old index
 *   n2o        --> new to old values
 *----------------------------------------------------------------------------*/

static void
_build_n2o_indexed(cs_lnum_t          n_old,
                   cs_lnum_t          n_new,
                   const cs_lnum_t    o2n[],
                   cs_lnum_t         *n2o_idx[],
                   cs_lnum_t         *n2o[])
{
  cs_lnum_t *_n2o_idx, *_n2o;

  BFT_MALLOC(_n2o_idx, n_new+1, cs_lnum_t);

  for (cs_lnum_t i = 0; i < n_new+1; i++)
    _n2o_idx[i] = 0;

  /* Count occurences */

  for (cs_lnum_t i = 0; i < n_old; i++)
    _n2o_idx[o2n[i] + 1] += 1;

  /* Transform count to index */

  for (cs_lnum_t i = 0; i < n_new; i++)
    _n2o_idx[i+1] += _n2o_idx[i];

  /* Build array */

  BFT_MALLOC(_n2o, _n2o_idx[n_new], cs_lnum_t);

  cs_lnum_t *shift;
  BFT_MALLOC(shift, n_new, cs_lnum_t);
  for (cs_lnum_t i = 0; i < n_new; i++)
    shift[i] = 0;

  for (cs_lnum_t i = 0; i < n_old; i++) {
    cs_lnum_t j = o2n[i];
    _n2o[_n2o_idx[j] + shift[j]] = i;
    shift[j] += 1;
  }

  BFT_FREE(shift);

  /* Set return values */

  *n2o_idx = _n2o_idx;
  *n2o = _n2o;
}

/*----------------------------------------------------------------------------
 * Update a global numbering array in case of entity renumbering
 *
 * parameters:
 *   n_new      <-- new number of elements
 *   n2o        <-- new to old array (same as old element ids list)
 *   global_num <-> global numbering (allocated if initially NULL)
 *
 * returns:
 *   new global number of elements
 *----------------------------------------------------------------------------*/

static cs_gnum_t
_n2o_update_global_num(cs_lnum_t          n_new,
                       const cs_lnum_t    n2o[],
                       cs_gnum_t        **global_num)
{
  cs_gnum_t n_g_new = n_new;

  if (cs_glob_n_ranks == 1 && *global_num == NULL)
    return n_g_new;

  fvm_io_num_t *n_io_num
    = fvm_io_num_create_from_select(n2o, *global_num, n_new, 0);

  BFT_FREE(*global_num);

  *global_num = fvm_io_num_transfer_global_num(n_io_num);

  n_g_new = fvm_io_num_get_global_count(n_io_num);

  n_io_num = fvm_io_num_destroy(n_io_num);

  return n_g_new;
}

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

  BFT_MALLOC(i_face_cells, n_new, cs_lnum_2_t);
  BFT_MALLOC(i_face_family, n_new, int);
  BFT_MALLOC(i_face_r_gen, n_new, char);

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

  BFT_FREE(m->i_face_r_gen);
  BFT_FREE(m->i_face_family);
  BFT_FREE(m->i_face_cells);
  m->i_face_r_gen = i_face_r_gen;
  m->i_face_cells = i_face_cells;
  m->i_face_family = i_face_family;

  /* Update global numbering */

  m->n_g_i_faces
    = _n2o_update_global_num(n_new, f_n2o, &(m->global_i_face_num));

  m->n_i_faces = n_new;

  /* Update connectivity */

  cs_lnum_t *i_face_vtx_idx, *i_face_vtx;
  BFT_MALLOC(i_face_vtx_idx, n_new+1, cs_lnum_t);
  BFT_MALLOC(i_face_vtx, m->i_face_vtx_connect_size, cs_lnum_t);

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

  BFT_FREE(m->i_face_vtx_idx);
  BFT_FREE(m->i_face_vtx_lst);

  m->i_face_vtx_idx = i_face_vtx_idx;
  m->i_face_vtx_lst = i_face_vtx;

  i_face_vtx_idx = NULL;
  i_face_vtx = NULL;

  m->i_face_vtx_connect_size = m->i_face_vtx_idx[n_new];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Merge cells based on renumbering array.
 *
 * Interior faces separating merged cells are removed.
 *
 * \param[in, out]  m      mesh
 * \param[in]       n_new  new number of cells
 * \param[in]       c_o2n  cell old to new renumbering
 */
/*----------------------------------------------------------------------------*/

static void
_merge_cells(cs_mesh_t       *m,
             cs_lnum_t        n_new,
             const cs_lnum_t  c_o2n[])
{
  const cs_lnum_t n_old = m->n_cells;

  cs_lnum_t *c_n2o = _build_n2o(n_old, n_new, c_o2n);

  int  *cell_family;
  BFT_MALLOC(cell_family, n_new, int);
  for (cs_lnum_t i = 0; i < n_new; i++) {
    cs_lnum_t j = c_n2o[i];
    cell_family[i] = m->cell_family[j];
  }
  BFT_FREE(m->cell_family);
  m->cell_family = cell_family;
  cell_family = NULL;

  /* Update global numbering */

  m->n_g_cells
    = _n2o_update_global_num(n_new, c_n2o, &(m->global_cell_num));

  BFT_FREE(c_n2o);

  /* Transfer (cell-based) halo information to (face-based) mesh builder
     in case of periodicity, before operation modifying cell numbering  */

  cs_mesh_builder_t *mb = NULL;

  if (m->halo != NULL) {
    if (m->n_init_perio > 0) {
      const cs_gnum_t n_g_faces = m->n_g_i_faces + m->n_g_b_faces;
      int rank_id = CS_MAX(cs_glob_rank_id, 0);
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

    if (mb != NULL)
      cs_mesh_builder_destroy(&mb);
  }

  /* Remove excess interior faces */

  cs_lnum_t n_i_faces_new = 0;

  {
    cs_lnum_t  *i_f_n2o;
    BFT_MALLOC(i_f_n2o, m->n_i_faces, cs_lnum_t);

    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
      cs_lnum_t i0 = m->i_face_cells[f_id][0];
      cs_lnum_t i1 = m->i_face_cells[f_id][1];
      if (i0 != i1) {
        i_f_n2o[n_i_faces_new] = f_id;
        n_i_faces_new++;
      }
    }

    _update_i_face_arrays(m, n_i_faces_new, i_f_n2o);

    BFT_FREE(i_f_n2o);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create face merge state structure.
 *
 * \param[in]  need_e2f  true if edge to faces array is needed.
 *
 * \return  face merge helper state structure
 */
/*----------------------------------------------------------------------------*/

static  cs_mesh_face_merge_state_t *
_face_merge_state_create(bool need_e2f)
{
  cs_mesh_face_merge_state_t *s;
  BFT_MALLOC(s, 1, cs_mesh_face_merge_state_t);

  s->n_edges = 0;
  s->n_edges_max = 3;
  s->n_vertices = 0;
  s->n_vertices_max = 3;

  BFT_MALLOC(s->face_vertices, s->n_vertices_max, cs_lnum_t);
  BFT_MALLOC(s->e2v, s->n_edges_max, cs_lnum_2_t);
  if (need_e2f)
    BFT_MALLOC(s->e2f, s->n_edges_max, cs_lnum_t);
  else
    s->e2f = NULL;

  return s;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy face merge state structure.
 *
 * \param[in, out]  s  pointer to face merge state structure
 */
/*----------------------------------------------------------------------------*/

static void
_face_merge_state_destroy(cs_mesh_face_merge_state_t  **s)
{
  if (*s != NULL) {
    cs_mesh_face_merge_state_t  *_s = *s;

    BFT_FREE(_s->face_vertices);
    BFT_FREE(_s->e2v);
    BFT_FREE(_s->e2f);

    BFT_FREE(*s);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Remove edge from face merge state structure.
 *
 * \param[in, out]  s        face merge state structure
 * \param[in, out]  edge_id  if of edge to remove
 */
/*----------------------------------------------------------------------------*/

static inline void
_face_merge_state_remove_edge(cs_mesh_face_merge_state_t  *s,
                              cs_lnum_t                    edge_id)
{
  assert(edge_id < s->n_edges);

  /* Swap with last */

  cs_lnum_t i = s->n_edges - 1;

  s->e2v[edge_id][0] = s->e2v[i][0];
  s->e2v[edge_id][1] = s->e2v[i][1];

  s->n_edges -= 1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Remove edge from face merge state structure with edge to face info.
 *
 * \param[in, out]  s        face merge state structure
 * \param[in, out]  edge_id  if of edge to remove
 */
/*----------------------------------------------------------------------------*/

static inline void
_face_merge_state_remove_edge_e2f(cs_mesh_face_merge_state_t  *s,
                                  cs_lnum_t                    edge_id)
{
  assert(edge_id < s->n_edges);

  /* Swap with last */

  cs_lnum_t i = s->n_edges - 1;

  s->e2v[edge_id][0] = s->e2v[i][0];
  s->e2v[edge_id][1] = s->e2v[i][1];
  s->e2f[edge_id] = s->e2f[i];

  s->n_edges -= 1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Merge existing faces into a new face
 *
 * This function also transforms some counts to indexes
 * (grouping some conversions so as to reduce the number of required loops).
 *
 * \param[in]       n_faces      number of faces to merge
 * \param[in]       face_ids     ids of faces to merge
 * \param[in]       face_orient  optional face orientation (1/-1), or NULL
 * \param[in]       f2v_idx      face->vertices index
 * \param[in]       f2v_lst      face->vertices connectivity
 * \param[in, out]  s            face merge state helper structure
 *
 * \return  new number of vertices, or -1 in case of error
 */
/*----------------------------------------------------------------------------*/

static cs_lnum_t
_build_merged_face(cs_lnum_t                    n_faces,
                   cs_lnum_t                    face_ids[],
                   short int                   *face_orient,
                   cs_lnum_t                    f2v_idx[],
                   cs_lnum_t                    f2v_lst[],
                   cs_mesh_face_merge_state_t  *s)
{
  s->n_edges = 0;
  s->n_vertices = 0;

  /* Loop over faces */

  for (cs_lnum_t f_i = 0; f_i < n_faces; f_i++) {

    cs_lnum_t f_id = face_ids[f_i];

    int orient = 1;
    if (face_orient != NULL)
      orient = face_orient[f_i];

    cs_lnum_t s_id = f2v_idx[f_id];
    cs_lnum_t e_id = f2v_idx[f_id+1];
    cs_lnum_t n_f_vtx = e_id - s_id;

    cs_lnum_t e_vtx[2];

    /* Loop over face edges */

    for (cs_lnum_t i = 0; i < n_f_vtx; i++) {

      if (orient > 0) {
        e_vtx[0] = f2v_lst[s_id + i];
        e_vtx[1] = f2v_lst[s_id + (i+1)%n_f_vtx];
      }
      else {
        e_vtx[0] = f2v_lst[s_id + (n_f_vtx-i)%n_f_vtx];
        e_vtx[1] = f2v_lst[s_id + (n_f_vtx-i-1)];
      }

      /* Compare to edges list: add if not present, cancel if present */

      bool insert = true;

      for (cs_lnum_t j = 0; j < s->n_edges; j++) {
        /* If present, must be in opposite direction */
        if (e_vtx[1] == s->e2v[j][0] && e_vtx[0] == s->e2v[j][1]) {
          _face_merge_state_remove_edge(s, j);
          insert = false;
          break;
        }
      }

      if (insert) {
        if (s->n_edges >= s->n_edges_max) {
          s->n_edges_max *= 2;
          BFT_REALLOC(s->e2v, s->n_edges_max, cs_lnum_2_t);
        }
        cs_lnum_t j = s->n_edges;
        s->e2v[j][0] = e_vtx[0];
        s->e2v[j][1] = e_vtx[1];
        s->n_edges += 1;
      }

    }

  }

  /* We should now be able to rebuild the face */

  if (s->n_edges > 1) {

    s->face_vertices[0] = s->e2v[0][0];
    s->face_vertices[1] = s->e2v[0][1];

    s->n_vertices = 2;

    _face_merge_state_remove_edge(s, 0);

    while (s->n_edges > 0) {
      bool matched = false;
      for (cs_lnum_t i = 0; i < s->n_edges; i++) {
        if (s->e2v[i][0] == s->face_vertices[s->n_vertices-1]) {
          s->n_vertices += 1;
          if (s->n_vertices > s->n_vertices_max) {
            s->n_vertices_max *= 2;
            BFT_REALLOC(s->face_vertices, s->n_vertices_max, cs_lnum_t);
          }
          s->face_vertices[s->n_vertices - 1] = s->e2v[i][1];
          _face_merge_state_remove_edge(s, i);
          matched = true;
          break;
        }
      }
      if (matched == false)
        bft_error(__FILE__, __LINE__, 0,
                  _("%s: %d edges do not constitute a closed contour."),
                  __func__, (int)n_faces);
    }

  }

  /* Verification */

  if (s->face_vertices[s->n_vertices - 1] == s->face_vertices[0])
    s->n_vertices -= 1;

  else
    bft_error(__FILE__, __LINE__, 0,
              _("%s: %d edges do not constitute a closed contour."),
              __func__, (int)n_faces);

  return s->n_vertices;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build interior faces equivalence array.
 *
 * \param[in, out]  m        mesh
 * \param[out]      i_f_o2n  face old to new renumbering
 *
 * return  number of new interior faces
 */
/*----------------------------------------------------------------------------*/

static cs_lnum_t
_i_faces_equiv(cs_mesh_t  *m,
               cs_lnum_t  *i_f_o2n[])
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_i_faces = m->n_i_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;

  /* Initialize arrays */

  cs_lnum_t  *_i_f_o2n;
  BFT_MALLOC(_i_f_o2n, n_i_faces, cs_lnum_t);

# pragma omp for
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++)
    _i_f_o2n[f_id] = f_id;

  /* Build cell->faces connectivity to group associated faces */

  cs_adjacency_t *c2f = cs_mesh_adjacency_c2f_lower(m, -1);

  cs_lnum_t *c2c;
  BFT_MALLOC(c2c, c2f->idx[n_cells], cs_lnum_t);

# pragma omp for schedule(dynamic, CS_CL_SIZE)
  for (cs_lnum_t c_id = 0; c_id < m->n_cells; c_id++) {

    cs_lnum_t s_id = c2f->idx[c_id];
    cs_lnum_t e_id = c2f->idx[c_id+1];

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t f_id = c2f->ids[j];

      cs_lnum_t c_id_0 = i_face_cells[f_id][0];
      cs_lnum_t c_id_1 = i_face_cells[f_id][1];

      /* Handle face when referenced by the adjacent cell with lowest id */
      if (c_id_0 == c_id)
        c2c[j] = c_id_1;
      else {
        assert(c_id_1 == c_id);
        c2c[j] = c_id_0;
      }

    }

    cs_sort_coupled_shell(s_id, e_id, c2c, c2f->ids);

    /* Now identify series of faces adjacent to the same cell */

    cs_lnum_t c_id_prev = -1;
    cs_lnum_t f_id_eq = -1;
    cs_lnum_t s_id_c = s_id;

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t c_id_0 = c2c[j];

      if (c_id_0 != c_id_prev) {
        f_id_eq = c2f->ids[s_id_c];
        for (cs_lnum_t k = s_id_c; k < j; k++)
          f_id_eq = CS_MIN(f_id_eq, c2f->ids[k]);
        for (cs_lnum_t k = s_id_c; k < j; k++) {
          cs_lnum_t f_id = c2f->ids[k];
          _i_f_o2n[f_id] = f_id_eq;
        }
        s_id_c = j;
        c_id_prev = c_id_0;
      }

    }

    if (s_id_c < e_id) {
      f_id_eq = c2f->ids[s_id_c];

      for (cs_lnum_t k = s_id_c; k < e_id; k++)
        f_id_eq = CS_MIN(f_id_eq, c2f->ids[k]);
      for (cs_lnum_t k = s_id_c; k < e_id; k++) {
        cs_lnum_t f_id = c2f->ids[k];
        _i_f_o2n[f_id] = f_id_eq;
      }
    }

  }

  BFT_FREE(c2c);
  cs_adjacency_destroy(&c2f);

  /* Now compact renumbering array */

  cs_lnum_t n_i_f_new = 0;
  for (cs_lnum_t i = 0; i < n_i_faces; i++) {
    if (_i_f_o2n[i] == i) {
      _i_f_o2n[i] = n_i_f_new;
      n_i_f_new += 1;
    }
    else {
      assert(_i_f_o2n[i] < i);
      _i_f_o2n[i] = _i_f_o2n[_i_f_o2n[i]];
    }
  }

  *i_f_o2n = _i_f_o2n;

  return n_i_f_new;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Merge interior faces based on renumbering array.
 *
 * \param[in, out]  m      mesh
 * \param[in]       n_new  new number of faces
 * \param[in]       f_o2n  face old to new renumbering
 */
/*----------------------------------------------------------------------------*/

static void
_merge_i_faces(cs_mesh_t       *m,
               cs_lnum_t        n_new,
               const cs_lnum_t  f_o2n[])
{
  const cs_lnum_t n_old = m->n_i_faces;

  cs_lnum_t *n2o_idx = NULL, *n2o = NULL;
  _build_n2o_indexed(n_old, n_new, f_o2n, &n2o_idx, &n2o);

  cs_lnum_t *i_face_vtx_idx, *i_face_vtx;
  BFT_MALLOC(i_face_vtx_idx, n_new+1, cs_lnum_t);
  BFT_MALLOC(i_face_vtx, m->i_face_vtx_idx[n_old], cs_lnum_t);

  i_face_vtx_idx[0] = 0;

  cs_mesh_face_merge_state_t *s = _face_merge_state_create(false);

  cs_lnum_t n_s_f_max = 0;
  short int *f_orient = NULL;

  for (cs_lnum_t i_n = 0; i_n < n_new; i_n++) {

    cs_lnum_t s_id = n2o_idx[i_n], e_id = n2o_idx[i_n+1];
    cs_lnum_t n_s_faces = e_id - s_id;

    /* If face does not need to be merged, simply copy it */

    if (n_s_faces == 1) {

      cs_lnum_t f_id = n2o[s_id];

      cs_lnum_t v_s_id_src = m->i_face_vtx_idx[f_id];
      cs_lnum_t n_f_vtx = m->i_face_vtx_idx[f_id+1] - v_s_id_src;

      cs_lnum_t v_s_id_dst = i_face_vtx_idx[i_n];

      i_face_vtx_idx[i_n + 1] = v_s_id_dst + n_f_vtx;
      for (cs_lnum_t i = 0; i < n_f_vtx; i++)
        i_face_vtx[v_s_id_dst + i] = m->i_face_vtx_lst[v_s_id_src + i];

    }

    else {

      /* Build orientation array */

      if (n_s_faces > n_s_f_max) {
        n_s_f_max = n_s_faces;
        BFT_REALLOC(f_orient, n_s_f_max, short int);
      }
      cs_lnum_t c_id0_cur = m->i_face_cells[s_id][0];
      f_orient[0] = 1;

      for (cs_lnum_t i = 0; i < n_s_faces; i++) {
        cs_lnum_t j = n2o[s_id + i];
        if (m->i_face_cells[j][0] == c_id0_cur)
          f_orient[i] = 1;
        else {
          assert(m->i_face_cells[j][1] == c_id0_cur);
          f_orient[i] = -1;
        }

      }

      _build_merged_face(n_s_faces,
                         n2o + s_id,
                         f_orient,
                         m->i_face_vtx_idx,
                         m->i_face_vtx_lst,
                         s);

      cs_lnum_t v_s_id_dst = i_face_vtx_idx[i_n];
      i_face_vtx_idx[i_n + 1] = v_s_id_dst + s->n_vertices;

      for (cs_lnum_t i = 0; i < s->n_vertices; i++)
        i_face_vtx[v_s_id_dst + i] = s->face_vertices[i];

    }

  }

  BFT_FREE(f_orient);

  m->i_face_vtx_connect_size = i_face_vtx_idx[n_new];

  _face_merge_state_destroy(&s);

  BFT_REALLOC(i_face_vtx, i_face_vtx_idx[n_new], cs_lnum_t);

  BFT_FREE(m->i_face_vtx_idx);
  BFT_FREE(m->i_face_vtx_lst);

  m->i_face_vtx_idx = i_face_vtx_idx;
  m->i_face_vtx_lst = i_face_vtx;

  i_face_vtx_idx = NULL;
  i_face_vtx = NULL;

  /* Transform indexed new->old array to simple array */

  for (cs_lnum_t i = 0; i < n_new; i++) {
    cs_lnum_t j = n2o_idx[i];
    assert(j >= i);
    n2o[i] = n2o[j];
  }

  BFT_FREE(n2o_idx);

  /* Now update related arrays */

  cs_lnum_2_t  *i_face_cells;
  int  *i_face_family;
  char *i_face_r_gen;
  BFT_MALLOC(i_face_family, n_new, int);
  BFT_MALLOC(i_face_cells, n_new, cs_lnum_2_t);
  BFT_MALLOC(i_face_r_gen, n_new, char);

  for (cs_lnum_t i = 0; i < n_new; i++) {
    cs_lnum_t j = n2o[i];
    i_face_cells[i][0] = m->i_face_cells[j][0];
    i_face_cells[i][1] = m->i_face_cells[j][1];
    i_face_family[i] = m->i_face_family[j];
    i_face_r_gen[i] = m->i_face_r_gen[j];
  }

  BFT_FREE(m->i_face_cells);
  m->i_face_cells = i_face_cells;
  i_face_cells = NULL;

  BFT_FREE(m->i_face_family);
  m->i_face_family = i_face_family;
  i_face_family = NULL;

  BFT_FREE(m->i_face_r_gen);
  m->i_face_r_gen = i_face_r_gen;
  i_face_r_gen = NULL;

  /* Update global numbering */

  m->n_g_i_faces
    = _n2o_update_global_num(n_new, n2o, &(m->global_i_face_num));

  m->n_i_faces = n_new;

  BFT_FREE(n2o);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Filter equivalence info for boundary faces.
 *
 * Boundary faces of a same cell can only be merged if joined by an edge
 * of maximum vertex level for that face. Thi topological criterion allows
 * avoiding merging sub-faces from different original faces.
 *
 * \param[in]       n_faces      number of faces to filter
 * \param[in]       face_ids     ids of faces to filter
 * \param[in]       f2v_idx      face->vertices index
 * \param[in]       f2v_lst      face->vertices connectivity
 * \param[in]       vtx_lv       vertex level
 * \param[in, out]  face_equiv   face equivalence info (face ids in,
 *                               lowest merged face id out)
 * \param[in, out]  s            face merge state helper structure
 */
/*----------------------------------------------------------------------------*/

static void
_filter_b_face_equiv(cs_lnum_t                    n_faces,
                     cs_lnum_t                    face_ids[],
                     cs_lnum_t                    f2v_idx[],
                     cs_lnum_t                    f2v_lst[],
                     char                         vtx_lv[],
                     cs_lnum_t                    face_equiv[],
                     cs_mesh_face_merge_state_t  *s)
{
  s->n_edges = 0;
  s->n_vertices = 0;

  /* Determine reference level from adjacent vertex with
     highest generation */

  char lv_min = 127;
  char lv_max = 0;

  for (cs_lnum_t f_i = 0; f_i < n_faces; f_i++) {
    cs_lnum_t f_id = face_ids[f_i];
    cs_lnum_t s_id = f2v_idx[f_id];
    cs_lnum_t e_id = f2v_idx[f_id+1];

    for (cs_lnum_t i = s_id; i < e_id; i++) {
      char lv = vtx_lv[f2v_lst[i]];
      if (lv > lv_max)
        lv_max = lv;
      else if (lv < lv_min)
        lv_min = lv;
    }
  }

  if (lv_min == lv_max)
    return;

  /* Loop over cell's faces */

  for (cs_lnum_t f_i = 0; f_i < n_faces; f_i++) {

    cs_lnum_t f_id = face_ids[f_i];

    cs_lnum_t s_id = f2v_idx[f_id];
    cs_lnum_t e_id = f2v_idx[f_id+1];
    cs_lnum_t n_f_vtx = e_id - s_id;

    cs_lnum_t e_vtx[2];

    /* Check if one of the face's edges can be merged */

    for (cs_lnum_t i = 0; i < n_f_vtx; i++) {

      e_vtx[0] = f2v_lst[s_id + i];
      e_vtx[1] = f2v_lst[s_id + (i+1)%n_f_vtx];

      /* If edge is not of maximum level, it cannot join
         2 sub-faces */

      if (vtx_lv[e_vtx[0]] != lv_max || vtx_lv[e_vtx[1]] != lv_max)
        continue;

      /* Compare to edges list: add if not present, cancel if present */

      bool insert = true;

      for (cs_lnum_t j = 0; j < s->n_edges; j++) {
        /* If present, must be in opposite direction */
        if (e_vtx[1] == s->e2v[j][0] && e_vtx[0] == s->e2v[j][1]) {
          cs_lnum_t f_id_eq = s->e2f[j];
          cs_lnum_t eq_0 = face_equiv[f_id_eq];
          cs_lnum_t eq_1 = face_equiv[f_id];
          cs_lnum_t eq_lo = CS_MIN(eq_0, eq_1);
          for (cs_lnum_t k = 0; k < n_faces; k++) {
            cs_lnum_t l = face_ids[k];
            if (face_equiv[l] == eq_0 || face_equiv[l] == eq_1)
              face_equiv[l] = eq_lo;
          }
          //_propagate_equiv(f_id, f_id_eq, face_equiv);
          /* Remove from queue */
          insert = false;
          _face_merge_state_remove_edge_e2f(s, j);
          break;
        }
      }

      if (insert) {
        if (s->n_edges >= s->n_edges_max) {
          s->n_edges_max *= 2;
          BFT_REALLOC(s->e2v, s->n_edges_max, cs_lnum_2_t);
          BFT_REALLOC(s->e2f, s->n_edges_max, cs_lnum_t);
        }
        cs_lnum_t j = s->n_edges;
        s->e2v[j][0] = e_vtx[0];
        s->e2v[j][1] = e_vtx[1];
        s->e2f[j] = f_id;
        s->n_edges += 1;
      }

    }

  }

  s->n_edges = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build boundary faces equivalence array.
 *
 * The caller is reqponsible for freeing the b_f_o2n array.
 *
 * \param[in, out]  m        mesh
 * \param[in]       c_flag   0 for unchanged cells, > 0 for coarsened ones
 * \param[out]      b_f_o2n  face old to new renumbering
 *
 * return  number of new interior faces
 */
/*----------------------------------------------------------------------------*/

static cs_lnum_t
_b_faces_equiv(cs_mesh_t  *m,
               const int   c_flag[],
               cs_lnum_t  *b_f_o2n[])
{
  const cs_lnum_t n_b_faces = m->n_b_faces;

  /* Initialize arrays */

  cs_lnum_t  *_b_f_o2n;
  BFT_MALLOC(_b_f_o2n, n_b_faces, cs_lnum_t);

# pragma omp for
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++)
    _b_f_o2n[f_id] = f_id;

  /* Build cell->faces connectivity to group associated faces */

  cs_adjacency_t *c2f = cs_mesh_adjacency_c2f_boundary(m);

  cs_mesh_face_merge_state_t *s = _face_merge_state_create(true);

  /* Now filter build faces equivalence */

  for (cs_lnum_t c_id = 0; c_id < m->n_cells; c_id++) {

    cs_lnum_t s_id = c2f->idx[c_id];
    cs_lnum_t e_id = c2f->idx[c_id+1];

    if ((c_flag[c_id] > 0) && (e_id - s_id > 1))
      _filter_b_face_equiv(e_id - s_id,
                           c2f->ids + s_id,
                           m->b_face_vtx_idx,
                           m->b_face_vtx_lst,
                           m->vtx_r_gen,
                           _b_f_o2n,
                           s);

  }

  _face_merge_state_destroy(&s);

  cs_adjacency_destroy(&c2f);

  /* Now compact renumbering array */

  cs_lnum_t n_b_f_new = 0;
  for (cs_lnum_t i = 0; i < n_b_faces; i++) {
    if (_b_f_o2n[i] == i) {
      _b_f_o2n[i] = n_b_f_new;
      n_b_f_new += 1;
    }
    else {
      assert(_b_f_o2n[i] < i);
      _b_f_o2n[i] = _b_f_o2n[_b_f_o2n[i]];
    }
  }

  *b_f_o2n = _b_f_o2n;

  return n_b_f_new;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Merge boundary faces based on renumbering array.
 *
 * \param[in, out]  m      mesh
 * \param[in]       n_new  new number of faces
 * \param[in]       f_o2n  face old to new renumbering
 */
/*----------------------------------------------------------------------------*/

static void
_merge_b_faces(cs_mesh_t       *m,
               cs_lnum_t        n_new,
               const cs_lnum_t  f_o2n[])
{
  const cs_lnum_t n_old = m->n_b_faces;

  cs_lnum_t *n2o_idx = NULL, *n2o = NULL;
  _build_n2o_indexed(n_old, n_new, f_o2n, &n2o_idx, &n2o);

  cs_lnum_t *b_face_vtx_idx, *b_face_vtx;
  BFT_MALLOC(b_face_vtx_idx, n_new+1, cs_lnum_t);
  BFT_MALLOC(b_face_vtx, m->b_face_vtx_idx[n_old], cs_lnum_t);

  b_face_vtx_idx[0] = 0;

  cs_mesh_face_merge_state_t *s = _face_merge_state_create(false);

  cs_lnum_t n_s_f_max = 0;
  short int *f_orient = NULL;

  for (cs_lnum_t i_n = 0; i_n < n_new; i_n++) {

    cs_lnum_t s_id = n2o_idx[i_n], e_id = n2o_idx[i_n+1];
    cs_lnum_t n_s_faces = e_id - s_id;

    /* If face does not need to be merged, simply copy it */

    if (n_s_faces == 1) {

      cs_lnum_t f_id = n2o[s_id];

      cs_lnum_t v_s_id_src = m->b_face_vtx_idx[f_id];
      cs_lnum_t n_f_vtx = m->b_face_vtx_idx[f_id+1] - v_s_id_src;

      cs_lnum_t v_s_id_dst = b_face_vtx_idx[i_n];

      b_face_vtx_idx[i_n + 1] = v_s_id_dst + n_f_vtx;
      for (cs_lnum_t i = 0; i < n_f_vtx; i++)
        b_face_vtx[v_s_id_dst + i] = m->b_face_vtx_lst[v_s_id_src + i];

    }

    else {

      /* Build orientation array (boundary faces always outwards oriented) */

      if (n_s_faces > n_s_f_max) {
        n_s_f_max = n_s_faces;
        BFT_REALLOC(f_orient, n_s_f_max, short int);
        for (cs_lnum_t i = 0; i < n_s_f_max; i++)
          f_orient[i] = 1;
      }

      _build_merged_face(n_s_faces,
                         n2o + s_id,
                         f_orient,
                         m->b_face_vtx_idx,
                         m->b_face_vtx_lst,
                         s);

      cs_lnum_t v_s_id_dst = b_face_vtx_idx[i_n];
      b_face_vtx_idx[i_n + 1] = v_s_id_dst + s->n_vertices;

      for (cs_lnum_t i = 0; i < s->n_vertices; i++)
        b_face_vtx[v_s_id_dst + i] = s->face_vertices[i];

    }

  }

  BFT_FREE(f_orient);

  m->b_face_vtx_connect_size = b_face_vtx_idx[n_new];

  _face_merge_state_destroy(&s);

  BFT_REALLOC(b_face_vtx, b_face_vtx_idx[n_new], cs_lnum_t);

  BFT_FREE(m->b_face_vtx_idx);
  BFT_FREE(m->b_face_vtx_lst);

  m->b_face_vtx_idx = b_face_vtx_idx;
  m->b_face_vtx_lst = b_face_vtx;

  b_face_vtx_idx = NULL;
  b_face_vtx = NULL;

  /* Transform indexed new->old array to simple array */

  for (cs_lnum_t i = 0; i < n_new; i++) {
    cs_lnum_t j = n2o_idx[i];
    assert(j >= i);
    n2o[i] = n2o[j];
  }

  BFT_FREE(n2o_idx);

  /* Now update related arrays */

  cs_lnum_t  *b_face_cells;
  int  *b_face_family;
  BFT_MALLOC(b_face_family, n_new, int);
  BFT_MALLOC(b_face_cells, n_new, cs_lnum_t);

  for (cs_lnum_t i = 0; i < n_new; i++) {
    cs_lnum_t j = n2o[i];
    b_face_cells[i] = m->b_face_cells[j];
    b_face_family[i] = m->b_face_family[j];
  }

  BFT_FREE(m->b_face_cells);
  m->b_face_cells = b_face_cells;
  b_face_cells = NULL;

  BFT_FREE(m->b_face_family);
  m->b_face_family = b_face_family;
  b_face_family = NULL;

  /* Update global numbering */

  m->n_g_b_faces
    = _n2o_update_global_num(n_new, n2o, &(m->global_b_face_num));

  m->n_b_faces = n_new;

  BFT_FREE(n2o);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Flag vertices belonging to coarsed cells.
 *
 * The caller is reqponsible for freeing the allocated array.
 *
 * \param[in, out]  m        mesh
 * \param[in]       c_flag   cell coarsening flag
 *
 * return  flag for vertices belonging to coarsened cells.
 */
/*----------------------------------------------------------------------------*/

static cs_lnum_t *
_flag_vertices(cs_mesh_t  *m,
               const int   c_flag[])
{
  const cs_lnum_t n_vtx = m->n_vertices;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_cells = m->n_cells;

  cs_lnum_t  *v_flag;
  BFT_MALLOC(v_flag, m->n_vertices, cs_lnum_t);

  for (cs_lnum_t v_id = 0; v_id < n_vtx; v_id++)
    v_flag[v_id] = 0;

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id_0 = m->i_face_cells[f_id][0];
    cs_lnum_t c_id_1 = m->i_face_cells[f_id][1];
    bool c_face = false;
    if (c_id_0 < n_cells) {
      if (c_flag[c_id_0] > 0)
        c_face = true;
    }
    if (c_id_1 < n_cells) {
      if (c_flag[c_id_1] > 0)
        c_face = true;
    }
    if (c_face) {
      cs_lnum_t s_id = m->i_face_vtx_idx[f_id];
      cs_lnum_t e_id = m->i_face_vtx_idx[f_id + 1];
      for (cs_lnum_t i = s_id; i < e_id; i++)
        v_flag[m->i_face_vtx_lst[i]] = 1;
    }
  }

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t c_id = m->b_face_cells[f_id];
    if (c_flag[c_id] > 0) {
      cs_lnum_t s_id = m->b_face_vtx_idx[f_id];
      cs_lnum_t e_id = m->b_face_vtx_idx[f_id + 1];
      for (cs_lnum_t i = s_id; i < e_id; i++)
        v_flag[m->b_face_vtx_lst[i]] = 1;
    }
  }

  cs_interface_set_t  *ifs = m->vtx_interfaces;
  if (ifs != NULL)
    cs_interface_set_max(ifs, m->n_vertices, 1, true, CS_LNUM_TYPE, v_flag);

  return v_flag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update vertex arrays.
 *
 * Adjencent elemtn connectivity is already handled at this stage, and most
 * auxiliary data not present here.
 *
 * \param[in, out]  m          mesh
 * \param[in]       n_vtx_new  number of new vertices
 * \param[in]       v_o2n      vertex old to new mapping
 */
/*----------------------------------------------------------------------------*/

static void
_update_vertices(cs_mesh_t  *m,
                 cs_lnum_t   n_vtx_new,
                 cs_lnum_t   v_o2n[])
{
  if (n_vtx_new < m->n_vertices) {

    /* Update vertex connectivity and global numbering */

    for (cs_lnum_t i = 0; i < m->n_vertices; i++) {
      const cs_lnum_t j = v_o2n[i];
      if (j != -1) {
        for (cs_lnum_t l = 0; l < 3; l++)
          m->vtx_coord[j*3 + l] = m->vtx_coord[i*3 + l];
        if (m->global_vtx_num != NULL)
          m->global_vtx_num[j] = m->global_vtx_num[i];
      }
    }

    /* Update mesh structure */

    m->n_vertices = n_vtx_new;

    BFT_REALLOC(m->vtx_coord, n_vtx_new*3, cs_real_t);

    if (m->global_vtx_num != NULL)
      BFT_REALLOC(m->global_vtx_num, n_vtx_new, cs_gnum_t);
  }

  if (m->vtx_interfaces != NULL)
    cs_interface_set_renumber(m->vtx_interfaces, v_o2n);

  /* Build an I/O numbering to compact the global numbering */

  if (cs_glob_n_ranks > 1) {

    fvm_io_num_t *tmp_num = fvm_io_num_create(NULL,
                                              m->global_vtx_num,
                                              m->n_vertices,
                                              0);

    if (m->n_vertices > 0)
      memcpy(m->global_vtx_num,
             fvm_io_num_get_global_num(tmp_num),
             m->n_vertices*sizeof(cs_gnum_t));

    m->n_g_vertices = fvm_io_num_get_global_count(tmp_num);

    assert(fvm_io_num_get_local_count(tmp_num) == (cs_lnum_t)m->n_vertices);

    tmp_num = fvm_io_num_destroy(tmp_num);

  }

  else { /* if (cs_glob_ranks == 1) */

    m->n_g_vertices = m->n_vertices;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Filter vertices belonging to coarsed cells.
 *
 * \param[in, out]  m      mesh
 * \param[in, out]  v_o2n  vertex coarsening flag in, old to new out
 */
/*----------------------------------------------------------------------------*/

static void
_filter_vertices(cs_mesh_t  *m,
                 cs_lnum_t   v_o2n[])
{
  /* Build vertices to vertices (edges) graph */

  cs_adjacency_t *v2v = cs_mesh_adjacency_v2v(m);
  const cs_lnum_t n_vertices = v2v->n_elts;

  /* Graph is symmetric (upper trianglular form), so we compute
      neighbor counts for all vertices. */

  int *v_n_count;
  BFT_MALLOC(v_n_count, n_vertices, int);

  for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++)
    v_n_count[v_id] = 0;

  /* For a single domain, simply count connections */

  if (m->n_domains == 1) {
    for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++) {
      cs_lnum_t s_id = v2v->idx[v_id];
      cs_lnum_t e_id = v2v->idx[v_id+1];
      v_n_count[v_id] += e_id - s_id;
      for (cs_lnum_t i = s_id; i < e_id; i++) {
        assert(v2v->ids[i] > v_id);
        v_n_count[v2v->ids[i]] += 1;
      }
    }
  }

  /* In parallel, this is more tricky, as some edges may connect to other
     domains. If a vertex has at most 2 or more neighbors, no rank should
     have a connection to it from a vertex whose id is neither the minimum
     or maximum adjacent neighbor, though. So we can exploit this property
     to check if the number of neighbors is more than 2 or not (the actual
     count is unimportant, only the fact that a vertex with only 2 neighbors
     is in the middle of an edge, which can be simplified, and vertices
     with less neighbors are isolated (1 neighbor is not expected). */

  else {
    cs_interface_set_t  *ifs = m->vtx_interfaces;

    cs_gnum_t *v_adj_min, *v_adj_max;
    BFT_MALLOC(v_adj_min, n_vertices, cs_gnum_t);
    BFT_MALLOC(v_adj_max, n_vertices, cs_gnum_t);
    for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++) {
      v_adj_min[v_id] = m->n_g_vertices + 2;
      v_adj_max[v_id] = 0;
    }

    for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++) {
      cs_lnum_t s_id = v2v->idx[v_id];
      cs_lnum_t e_id = v2v->idx[v_id+1];
      cs_gnum_t g_v_id = m->global_vtx_num[v_id];
      v_n_count[v_id] += e_id - s_id;
      for (cs_lnum_t i = s_id; i < e_id; i++) {
        cs_lnum_t v_id_1 = v2v->ids[i];
        cs_gnum_t g_v_id_1 = m->global_vtx_num[v_id_1];
        v_adj_min[v_id] = CS_MIN(v_adj_min[v_id], g_v_id_1);
        v_adj_min[v_id_1] = CS_MIN(v_adj_min[v_id_1], g_v_id);
        v_adj_max[v_id] = CS_MAX(v_adj_max[v_id], g_v_id_1);
        v_adj_max[v_id_1] = CS_MAX(v_adj_max[v_id_1], g_v_id);
        v_n_count[v_id_1] += 1;
      }
    }

    cs_interface_set_max(ifs, n_vertices, 1, true, CS_INT_TYPE, v_n_count);
    cs_interface_set_max(ifs, n_vertices, 1, true, CS_GNUM_TYPE, v_adj_max);
    cs_interface_set_min(ifs, n_vertices, 1, true, CS_GNUM_TYPE, v_adj_min);

    for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++) {
      cs_lnum_t s_id = v2v->idx[v_id];
      cs_lnum_t e_id = v2v->idx[v_id+1];
      cs_gnum_t g_v_id = m->global_vtx_num[v_id];
      for (cs_lnum_t i = s_id; i < e_id; i++) {
        cs_lnum_t v_id_1 = v2v->ids[i];
        cs_gnum_t g_v_id_1 = m->global_vtx_num[v_id_1];
        if (g_v_id_1 != v_adj_min[v_id] && g_v_id_1 != v_adj_max[v_id])
          v_n_count[v_id] = 1000;  /* arbitrary, articifial value, > 2 */
        if (g_v_id != v_adj_min[v_id_1] && g_v_id != v_adj_max[v_id_1])
          v_n_count[v_id] = 1000;
      }
    }

    cs_interface_set_max(ifs, n_vertices, 1, true, CS_INT_TYPE, v_n_count);

    BFT_FREE(v_adj_min);
    BFT_FREE(v_adj_max);
  }

  cs_adjacency_destroy(&v2v);

  /* Vertices already flagged as potentially removed and having
     2 neighbors at most (edge centers) will be removed */

  cs_lnum_t n_vtx_new = 0;

  for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++) {
    if (v_o2n[v_id] > 0 && v_n_count[v_id] < 3) {
      v_o2n[v_id] = -1;
    }
    else {
      v_o2n[v_id] = n_vtx_new;
      n_vtx_new++;
    }
  }

  BFT_FREE(v_n_count);

  /* Now remove flagged vertices from face->vertex adjacencies */

  {
    const cs_lnum_t n_i_faces = m->n_i_faces;
    cs_lnum_t n = 0;
    cs_lnum_t s_id = 0;
    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
      cs_lnum_t e_id = m->i_face_vtx_idx[f_id + 1];
      for (cs_lnum_t i = s_id; i < e_id; i++) {
        cs_lnum_t vtx_id_n = v_o2n[m->i_face_vtx_lst[i]];
        if (vtx_id_n > -1) {
          m->i_face_vtx_lst[n] = vtx_id_n;
          //m->i_face_vtx_lst[n] = m->i_face_vtx_lst[i];
          n++;
        }
      }
      s_id = e_id;
      m->i_face_vtx_idx[f_id + 1] = n;
    }

    m->i_face_vtx_connect_size = n;
    BFT_REALLOC(m->i_face_vtx_lst, n, cs_lnum_t);
  }

  {
    const cs_lnum_t n_b_faces = m->n_b_faces;
    cs_lnum_t n = 0;
    cs_lnum_t s_id = 0;
    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      cs_lnum_t e_id = m->b_face_vtx_idx[f_id + 1];
      for (cs_lnum_t i = s_id; i < e_id; i++) {
        cs_lnum_t vtx_id_n = v_o2n[m->b_face_vtx_lst[i]];
        if (vtx_id_n > -1) {
          m->b_face_vtx_lst[n] = vtx_id_n;
          //m->b_face_vtx_lst[n] = m->b_face_vtx_lst[i];
          n++;
        }
      }
      s_id = e_id;
      m->b_face_vtx_idx[f_id + 1] = n;
    }

    m->b_face_vtx_connect_size = n;
    BFT_REALLOC(m->b_face_vtx_lst, n, cs_lnum_t);
  }

  /* Update other vertex info */

  _update_vertices(m, n_vtx_new, v_o2n);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Coarsen flagged mesh cells.
 *
 * \param[in, out]  m           mesh
 * \param[in]       cell_flag   coarsening flag for each cell
 *                              (0: do not coarsen; 1: coarsen)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_coarsen_simple(cs_mesh_t  *m,
                       const int   cell_flag[])
{
  /* Timers:
     0: total
  */

  cs_timer_counter_t  timers[2];
  for (int i = 0; i < 2; i++)
    CS_TIMER_COUNTER_INIT(timers[i]);

  /* Build ghosts in case they are not present */

  int mv_save = m->verbosity;
  m->verbosity = -1;

  if ((m->n_domains > 1 || m->n_init_perio > 0) && m->halo == NULL) {
    cs_halo_type_t halo_type = m->halo_type;
    cs_mesh_builder_t *mb = (m == cs_glob_mesh) ? cs_glob_mesh_builder : NULL;
    cs_mesh_init_halo(m, mb, halo_type, -1, true);
  }

  /* Free data that can be rebuilt */

  cs_mesh_free_rebuildable(m, false);

  if (m->vtx_range_set != NULL)
    cs_range_set_destroy(&(m->vtx_range_set));

  m->verbosity = mv_save;

  /* Logging */

  if (m->verbosity > 0) {
    cs_log_printf(CS_LOG_DEFAULT, "\n");
    cs_log_separator(CS_LOG_DEFAULT);
    _print_mesh_counts(cs_glob_mesh, _("Mesh before coarsening"));
  }

  cs_lnum_t n_c_ini = m->n_cells;

  // cs_lnum_t n_v_ini = m->n_vertices;
  // cs_lnum_t n_f_ini = m->n_b_faces + m->n_i_faces;
  // cs_lnum_t n_b_f_ini = m->n_b_faces;

  cs_timer_t t0 = cs_timer_time();

  cs_timer_t t2 = cs_timer_time();

  cs_timer_counter_add_diff(&(timers[0]), &t0, &t2);

  /* Determine cells that should be merged and flag associated vertices */

  cs_lnum_t  *c_o2n = NULL;
  cs_lnum_t  n_c_new = _cell_equiv(m, cell_flag, &c_o2n);
  cs_lnum_t *vtx_flag = _flag_vertices(m, cell_flag);

  _merge_cells(m, n_c_new, c_o2n);

  /* Flag merged cells (> 0 for merged cells) */

  int *c_flag_n = NULL;
  BFT_MALLOC(c_flag_n, n_c_new, int);

  for (cs_lnum_t i = 0; i < n_c_new; i++)
    c_flag_n[i] = -1;
  for (cs_lnum_t i = 0; i < n_c_ini; i++)
    c_flag_n[c_o2n[i]] += 1;

  /* Now merge interior faces */

  cs_lnum_t  *i_f_o2n = NULL;
  cs_lnum_t  n_i_f_new = _i_faces_equiv(m, &i_f_o2n);

  _merge_i_faces(m, n_i_f_new, i_f_o2n);

  BFT_FREE(i_f_o2n);
  BFT_FREE(c_o2n);

  /* Then merge boundary faces */

  cs_lnum_t  *b_f_o2n = NULL;
  cs_lnum_t  n_b_f_new = _b_faces_equiv(m, c_flag_n, &b_f_o2n);

  _merge_b_faces(m, n_b_f_new, b_f_o2n);

  BFT_FREE(b_f_o2n);

  BFT_FREE(c_flag_n);

  /* Finally remove excess vertices */

  cs_lnum_t *vtx_o2n = vtx_flag;   /* Rename for clarity as role changes */
  vtx_flag = NULL;

  _filter_vertices(m, vtx_o2n);

  BFT_FREE(vtx_o2n);

  m->modified |= (CS_MESH_MODIFIED | CS_MESH_MODIFIED_BALANCE);

  bft_printf("\nWarning mesh coarsening algorithm not complete yet\n");

  /* Rebuild auxiliary information  */

  mv_save = m->verbosity;
  m->verbosity = -1;

  cs_mesh_update_auxiliary(cs_glob_mesh);

  m->verbosity = mv_save;

  /* Update structures */

  if (m->verbosity > 0) {

    _print_mesh_counts(cs_glob_mesh, _("Mesh after coarsening"));
    cs_log_printf(CS_LOG_DEFAULT, "\n");
    cs_log_separator(CS_LOG_DEFAULT);

    cs_log_printf
      (CS_LOG_PERFORMANCE,
       _("\nMesh coarsening:\n\n"
         "  Total:                                        %.3g\n"),
       (double)(timers[0].nsec*1.e-9));
    cs_log_printf(CS_LOG_PERFORMANCE, "\n");
    cs_log_separator(CS_LOG_PERFORMANCE);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Coarsen selected mesh cells.
 *
 * \param[in, out]  m           mesh
 * \param[in]       n_cells     number of selected cells
 * \param[in]       cells       list of selected cells (0 to n-1)
 *                              or NULL if no indirection is needed
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_coarsen_simple_selected(cs_mesh_t        *m,
                                cs_lnum_t         n_cells,
                                const cs_lnum_t   cells[])
{
  cs_lnum_t n_c_ini = m->n_cells;

  int *cell_flag;
  BFT_MALLOC(cell_flag, n_c_ini, int);
  for (cs_lnum_t i = 0; i < n_c_ini; i++)
    cell_flag[i] = 0;

  if (cells != NULL) {
    for (cs_lnum_t i = 0; i < n_cells; i++)
      cell_flag[cells[i]] = 1;
  }
  else {
    for (cs_lnum_t i = 0; i < n_cells; i++)
      cell_flag[i] = 1;
  }

  cs_mesh_coarsen_simple(m, cell_flag);

  BFT_FREE(cell_flag);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
