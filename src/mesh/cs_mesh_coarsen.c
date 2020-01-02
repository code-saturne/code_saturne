/*============================================================================
 * Mesh coarsening.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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
#include "fvm_triangulate.h"

#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_adjacencies.h"
#include "cs_mesh_quantities.h"
#include "cs_order.h"
#include "cs_parall.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#if 1
#include "fvm_nodal.h"
#include "fvm_nodal_order.h"
#include "fvm_nodal_from_desc.h"
#include "fvm_writer.h"
#include "cs_mesh_connect.h"
#endif

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

  cs_lnum_t     n_edges;           /* number of edges */
  int           n_edges_max;       /* Maximum number of edges */
  int           n_vertices_max;    /* Maximum number of vertices */

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
 * \param[in, out]  m           mesh
 * \param[in]       cell_flag   coarsening type for each cell
 *                              (0: none; 1: isotropic)
 * \param[out]      c_o2n       cell old to new renumbering
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
 *   n_old      <-- new number of elements
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
 * \return  face merge helper state structur
 */
/*----------------------------------------------------------------------------*/

static  cs_mesh_face_merge_state_t *
_face_merge_state_create(void)
{
  cs_mesh_face_merge_state_t *s;
  BFT_MALLOC(s, 1, cs_mesh_face_merge_state_t);

  s->n_edges_max = 3;
  s->n_vertices_max = 3;

  BFT_MALLOC(s->face_vertices, s->n_vertices_max, cs_lnum_t);
  BFT_MALLOC(s->e2v, s->n_edges_max, cs_lnum_2_t);

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
  cs_lnum_t n_vertices = 0;

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

    n_vertices = 2;

    _face_merge_state_remove_edge(s, 0);

    while (s->n_edges > 0) {
      for (cs_lnum_t i = 0; i < s->n_edges; i++) {
        if (s->e2v[i][0] == s->face_vertices[n_vertices-1]) {
          n_vertices += 1;
          if (n_vertices > s->n_vertices_max) {
            s->n_vertices_max *= 2;
            BFT_REALLOC(s->face_vertices, s->n_vertices_max, cs_lnum_t);
          }
          s->face_vertices[n_vertices - 1] = s->e2v[0][1];
          _face_merge_state_remove_edge(s, i);
          break;
        }
      }
      /* We should not arrive here */
    }

  }

  /* Verification */

  return n_vertices;
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
 * \param[in]       cell_flag   subdivision type for each cell
 *                              (0: none; 1: isotropic)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_coarsen_simple(cs_mesh_t  *m,
                       const int   cell_flag[])
{
  /* Timers:
     0: total
  */

  cs_timer_counter_t  timers[1];
  for (int i = 0; i < 2; i++)
    CS_TIMER_COUNTER_INIT(timers[i]);

  /* Build ghosts in case they are not present */

  int mv_save = m->verbosity;
  m->verbosity = -1;

  if ((m->n_domains > 1 || m->n_init_perio > 0) && m->halo == NULL) {
    cs_halo_type_t halo_type = m->halo_type;
    cs_mesh_builder_t *mb = (m == cs_glob_mesh) ? cs_glob_mesh_builder : NULL;
    cs_mesh_init_halo(m, mb, halo_type);
    cs_mesh_update_auxiliary(cs_glob_mesh);
  }

  /* Free data that will be rebuilt */

  cs_mesh_free_rebuildable(m, true);

  m->verbosity = mv_save;

  if (m->verbosity > 0) {
    cs_log_printf(CS_LOG_DEFAULT, "\n");
    cs_log_separator(CS_LOG_DEFAULT);
    _print_mesh_counts(cs_glob_mesh, _("Mesh before refinement"));
  }

  cs_lnum_t n_v_ini = m->n_vertices;
  cs_lnum_t n_f_ini = m->n_b_faces + m->n_i_faces;

  cs_lnum_t n_b_f_ini = m->n_b_faces;

  cs_timer_t t0 = cs_timer_time();

  cs_timer_t t2 = cs_timer_time();

  cs_timer_counter_add_diff(&(timers[0]), &t0, &t2);

  /* Determine cells that should be merged */

  cs_lnum_t  *c_o2n = NULL;
  cs_lnum_t  n_c_new = _cell_equiv(m, cell_flag, &c_o2n);

  _merge_cells(m, n_c_new, c_o2n);

  BFT_FREE(c_o2n);

  m->modified = CS_MAX(m->modified, 1);

  bft_printf("\nWarning mesh coarsening algorithm not complete yet\n");

  /* Rebuild ghosts */

  mv_save = m->verbosity;
  m->verbosity = -1;

  if (   m->n_domains > 1 || m->n_init_perio > 0
      || m->halo_type == CS_HALO_EXTENDED) {
    cs_halo_type_t halo_type = m->halo_type;
    assert(m == cs_glob_mesh);
    cs_mesh_builder_t *mb = (m == cs_glob_mesh) ? cs_glob_mesh_builder : NULL;
    cs_mesh_init_halo(m, mb, halo_type);
  }

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
       (double)(timers[0].wall_nsec*1.e-9));
    cs_log_printf(CS_LOG_PERFORMANCE, "\n");
    cs_log_separator(CS_LOG_PERFORMANCE);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Refine selected mesh cells.
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
