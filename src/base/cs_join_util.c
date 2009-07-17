/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 2008-2009 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *===========================================================================*/

/*============================================================================
 * Manipulation of low-level structures for the joining operations
 *===========================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *---------------------------------------------------------------------------*/

#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *---------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *---------------------------------------------------------------------------*/

#include <fvm_parall.h>
#include <fvm_io_num.h>
#include <fvm_order.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *---------------------------------------------------------------------------*/

#include "cs_join_util.h"
#include "cs_mesh.h"
#include "cs_search.h"
#include "cs_sort.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *---------------------------------------------------------------------------*/

#include "cs_join_set.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro and type definitions
 *===========================================================================*/

/*============================================================================
 * Private function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Reduce numbering for the selected border faces.
 * After this function, we have a compact global face numbering for the
 * selected faces.
 *
 * parameters:
 *   n_select_faces      <-- number of selected faces
 *   p_reduce_gnum       <-> pointer to the reduced numbering
 *   p_reduce_gnum_index <-> pointer to an index on ranks for the reduce
 *                           numbering
 *---------------------------------------------------------------------------*/

static void
_compact_face_gnum_selection(cs_int_t     n_select_faces,
                             fvm_gnum_t  *reduce_gnum[],
                             fvm_gnum_t  *reduce_gnum_index[])
{
  cs_int_t  i;

  fvm_gnum_t  shift = 0;
  fvm_gnum_t  *_reduce_gnum = *reduce_gnum;
  fvm_gnum_t  *_reduce_gnum_index = *reduce_gnum_index;

  const int  n_ranks = cs_glob_n_ranks;
  const int  local_rank = CS_MAX(cs_glob_rank_id, 0);

  assert(_reduce_gnum_index == NULL);

  BFT_MALLOC(_reduce_gnum_index, n_ranks + 1, fvm_gnum_t);

  for (i = 0; i < n_ranks; i++)
    _reduce_gnum_index[i] = 0;

  if (n_ranks > 1) {
#if defined(HAVE_MPI)
    MPI_Comm  mpi_comm = cs_glob_mpi_comm;
    fvm_gnum_t  _n_faces = n_select_faces;

    MPI_Allgather(&_n_faces, 1, FVM_MPI_GNUM,
                  &(_reduce_gnum_index[1]), 1, FVM_MPI_GNUM, mpi_comm);
#endif

    for (i = 0; i < n_ranks; i++)
      _reduce_gnum_index[i+1] += _reduce_gnum_index[i];

    shift = _reduce_gnum_index[local_rank];

  }
  else {

    assert(n_ranks == 1);
    _reduce_gnum_index[n_ranks] = (fvm_gnum_t)n_select_faces;

  }

  BFT_MALLOC(_reduce_gnum, n_select_faces, fvm_gnum_t);

  for (i = 0; i < n_select_faces; i++)
    _reduce_gnum[i] = shift + i + 1;

  /* Returns pointer */

  *reduce_gnum = _reduce_gnum;
  *reduce_gnum_index = _reduce_gnum_index;
}

/*----------------------------------------------------------------------------
 * Extract vertices involved in the current joining operation
 *
 * parameters:
 *   n_select_faces    <-- number of selected faces
 *   select_faces      <-- list of faces selected
 *   f2v_idx           <-- "face -> vertex" connect. index
 *   f2v_lst           <-- "face -> vertex" connect. list
 *   n_vertices        <-- number of vertices
 *   n_select_vertices <-> pointer to the number of selected vertices
 *   select_vertices   <-> pointer to the list of selected vertices
 *---------------------------------------------------------------------------*/

static void
_extract_vertices(cs_int_t         n_select_faces,
                  const cs_int_t  *select_faces,
                  const cs_int_t  *f2v_idx,
                  const cs_int_t  *f2v_lst,
                  cs_int_t         n_vertices,
                  cs_int_t        *n_select_vertices,
                  cs_int_t        *select_vertices[])
{
  int  i, j, face_id;

  cs_int_t  _n_select_vertices = 0;
  cs_int_t  *counter = NULL, *_select_vertices = NULL;

  if (n_select_faces > 0) {

    BFT_MALLOC(counter, n_vertices, cs_int_t);

    for (i = 0; i < n_vertices; i++)
      counter[i] = 0;

    for (i = 0; i < n_select_faces; i++) {

      face_id = select_faces[i] - 1;

      for (j = f2v_idx[face_id] - 1; j < f2v_idx[face_id+1] - 1; j++)
        counter[f2v_lst[j]-1] = 1;

    }

    for (i = 0; i < n_vertices; i++)
      _n_select_vertices += counter[i];

    BFT_MALLOC(_select_vertices, _n_select_vertices, cs_int_t);

    _n_select_vertices = 0;
    for (i = 0; i < n_vertices; i++)
      if (counter[i] == 1)
        _select_vertices[_n_select_vertices++] = i + 1;

    assert(_n_select_vertices > 0);

    BFT_FREE(counter);

  } /* End if n_select_faces > 0 */

  /* Return pointers */

  *n_select_vertices = _n_select_vertices;
  *select_vertices = _select_vertices;
}

/*----------------------------------------------------------------------------
 * Extract faces implied in the current joining operation.
 * These are faces which share at least one vertex which is in the
 * select_vertices array.
 *
 * parameters:
 *   n_vertices        <-- number of vertices in the whole mesh
 *   n_select_vertices <-- number of selected vertices
 *   select_vertices   <-- list of selected vertices
 *   n_faces           <-- number of faces in the whole mesh
 *   f2v_idx           <-- "face -> vertex" connect. index
 *   f2v_lst           <-- "face -> vertex" connect. list
 *   n_contig_faces    <-> pointer to the number of contiguous faces
 *   contig_faces      <-> pointer to the list of contiguous faces
 *---------------------------------------------------------------------------*/

static void
_extract_contig_faces(cs_int_t         n_vertices,
                      cs_int_t         n_select_vertices,
                      const cs_int_t   select_vertices[],
                      cs_int_t         n_faces,
                      const cs_int_t   f2v_idx[],
                      const cs_int_t   f2v_lst[],
                      cs_int_t        *n_contig_faces,
                      cs_int_t        *contig_faces[])
{
  cs_int_t  i, j,  vtx_id, shift;

  cs_int_t   _n_contig_faces = 0;
  cs_int_t  *_contig_faces = NULL, *counter = NULL;
  cs_int_t  *v2f_idx = NULL, *v2f_lst = NULL;

  if (n_select_vertices == 0)
    return;

  /* Reverse face -> vertex connectivity */

  BFT_MALLOC(counter, n_vertices, cs_int_t);

  for (i = 0; i < n_vertices; i++)
    counter[i] = 0;

  for (i = 0; i < n_faces; i++) {
    for (j = f2v_idx[i] - 1; j < f2v_idx[i+1] - 1; j++) {
      vtx_id = f2v_lst[j] - 1;
      counter[vtx_id] += 1;
    }
  } /* End of loop on faces */

  /* Define v2f_idx */

  BFT_MALLOC(v2f_idx, n_vertices + 1, cs_int_t);

  v2f_idx[0] = 0;
  for (i = 0; i < n_vertices; i++)
    v2f_idx[i+1] = v2f_idx[i] + counter[i];

  for (i = 0; i < n_vertices; i++)
    counter[i] = 0;

  /* Define v2f_lst */

  BFT_MALLOC(v2f_lst, v2f_idx[n_vertices], cs_int_t);

  for (i = 0; i < n_faces; i++) {

    for (j = f2v_idx[i] - 1; j < f2v_idx[i+1] - 1; j++) {

      vtx_id = f2v_lst[j] - 1;
      shift = v2f_idx[vtx_id] + counter[vtx_id];
      v2f_lst[shift] = i+1;
      counter[vtx_id] += 1;

    }

  } /* End of loop on faces */

  BFT_REALLOC(counter, n_faces, cs_int_t);

  for (i = 0; i < n_faces; i++)
    counter[i] = 0;

  /* Count the number of contiguous faces */

  for (i = 0; i < n_select_vertices; i++) {

    vtx_id = select_vertices[i] - 1;

    for (j = v2f_idx[vtx_id]; j < v2f_idx[vtx_id+1]; j++)
      counter[v2f_lst[j]-1] = 1;

  }

  for (i = 0; i < n_faces; i++)
    _n_contig_faces += counter[i];

  /* Define contig_faces */

  BFT_MALLOC(_contig_faces, _n_contig_faces, cs_int_t);

  _n_contig_faces = 0;
  for (i = 0; i < n_faces; i++) {
    if (counter[i] == 1) {
      _contig_faces[_n_contig_faces] = i+1;
      _n_contig_faces += 1;
    }
  }

  /* Free memory */

  BFT_FREE(v2f_idx);
  BFT_FREE(v2f_lst);
  BFT_FREE(counter);

  /* Return pointers */

  *n_contig_faces = _n_contig_faces;
  *contig_faces = _contig_faces;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Get the related edge id from a couple of vertex ids.
 * Done only if the run is parallel.
 *
 * parameters:
 *   v1_id   <-- first vertex id
 *   v2_id   <-- second vertex id
 *   v2v_idx <-- vertex -> vertex connectivity index
 *   v2v_lst <-- vertex -> vertex connectivity
 *
 * return:
 *   related edge_id in cs_join_edges_t structure
 *---------------------------------------------------------------------------*/

inline static cs_int_t
_get_edge_id(cs_int_t        v1_id,
             cs_int_t        v2_id,
             const cs_int_t  v2v_idx[],
             const cs_int_t  v2v_lst[])
{
  int  i, va, vb;
  cs_int_t  edge_id = -1;

  if (v1_id < v2_id)
    va = v1_id, vb = v2_id;
  else
    vb = v1_id, va = v2_id;

  for (i = v2v_idx[va]; i < v2v_idx[va+1]; i++)
    if (v2v_lst[i] == vb + 1)
      break;

  if (i != v2v_idx[va+1])
    edge_id = i;

  return  edge_id;
}

/*----------------------------------------------------------------------------
 * Get the full selection of faces related to "single" vertices.
 * Done only if the run is parallel.
 *
 * parameters:
 *   selection_tag    <-- tag to know if a vertices is in selection
 *   v1_id            <-- first vertex id
 *   v2_id            <-- second vertex id
 *   ref_v2v_idx      <-- vertex -> vertex connect. index
 *   ref_v2v_lst      <-- vertex -> vertex connect. list
 *   max_size         <-> pointer to the max allocated size of new_s_vertices
 *   n_new_s_vertices <-> number of new "single" vertices
 *   new_s_vertices   <-> list of single vertices to add
 *---------------------------------------------------------------------------*/

static void
_add_s_vertex(const fvm_gnum_t   selection_tag[],
              cs_int_t           vid1,
              cs_int_t           vid2,
              const cs_int_t     ref_v2v_idx[],
              const cs_int_t     ref_v2v_lst[],
              cs_int_t          *max_size,
              cs_int_t          *n_new_s_vertices,
              cs_int_t          *new_s_vertices[])
{
  cs_int_t  edge_id;

  _Bool  in_selection = false;
  cs_int_t  _max_size = *max_size;
  cs_int_t  _n_new_s_vertices = *n_new_s_vertices;
  cs_int_t  *_new_s_vertices = *new_s_vertices;

  if (selection_tag[vid1] > 0)
    if (selection_tag[vid2] > 0)
      in_selection = true;

  if (in_selection == true) { /* Check if this edge is in ref_edges */

    edge_id = _get_edge_id(vid1, vid2, ref_v2v_idx, ref_v2v_lst);

    if (edge_id == -1) {

      _new_s_vertices[_n_new_s_vertices++] = vid1 + 1;

      if (_n_new_s_vertices >= _max_size) {
        _max_size *= 2;
        BFT_REALLOC(_new_s_vertices, _max_size, cs_int_t);
      }

    } /* this edge is not selected */

  } /* the two vertices are in the selection */

  /* Return pointers */

  *max_size = _max_size;
  *n_new_s_vertices = _n_new_s_vertices;
  *new_s_vertices = _new_s_vertices;
}

/*----------------------------------------------------------------------------
 * Get the full selection of faces related to "single" vertices.
 *
 * Done only if the run is parallel.
 *
 * parameters:
 *   n_vertices    <-- number of vertices in the parent mesh
 *   selection_tag <-> tag to know if a vertices is in selection
 *   selection     <-> pointer to a fvm_join_selection_t structure
 *   b_f2v_idx     <-- border "face -> vertex" connectivity index
 *   b_f2v_lst     <-- border "face -> vertex" connectivity
 *   i_f2v_idx     <-- interior "face -> vertex" connectivity
 *   i_f2v_lst     <-- interior "face -> vertex" connectivity
 *---------------------------------------------------------------------------*/

static void
_add_s_vertices_from_adj_faces(cs_int_t           n_vertices,
                               fvm_gnum_t         selection_tag[],
                               cs_join_select_t  *selection,
                               const cs_int_t    *b_f2v_idx,
                               const cs_int_t    *b_f2v_lst,
                               const cs_int_t    *i_f2v_idx,
                               const cs_int_t    *i_f2v_lst)
{
  cs_int_t  i, j, fid, save, shift, s, e, n_ref_edges;

  cs_int_t  n_new_s_vertices = 0, max_size = 10;
  cs_int_t  *count = NULL, *new_s_vertices = NULL;
  cs_int_t  *ref_v2v_idx = NULL, *ref_v2v_lst = NULL;

  /* Build reference edges */

  BFT_MALLOC(ref_v2v_idx, n_vertices + 1, cs_int_t);

  for (i = 0; i < n_vertices + 1; i++)
    ref_v2v_idx[i] = 0;

  cs_join_build_edges_idx(selection->n_faces,
                          selection->faces,
                          b_f2v_idx,
                          b_f2v_lst,
                          ref_v2v_idx);

  BFT_MALLOC(count, n_vertices, cs_int_t);

  for (i = 0; i < n_vertices; i++) {
    ref_v2v_idx[i+1] += ref_v2v_idx[i];
    count[i] = 0;
  }

  BFT_MALLOC(ref_v2v_lst, ref_v2v_idx[n_vertices], cs_int_t);

  cs_join_build_edges_lst(selection->n_faces,
                          selection->faces,
                          b_f2v_idx,
                          b_f2v_lst,
                          count,
                          ref_v2v_idx,
                          ref_v2v_lst);

  /* Ordering in order to clean the list */

  for (i = 0; i < n_vertices; i++)
    cs_sort_shell(ref_v2v_idx[i],
                  ref_v2v_idx[i+1],
                  ref_v2v_lst);

  /* Delete redundancies. Clean structure. */

  save = ref_v2v_idx[0];
  shift = 0;

  for (i = 0; i < n_vertices; i++) {

    s = save;
    e = ref_v2v_idx[i+1];

    if (e - s > 0) {
      ref_v2v_lst[shift++] = ref_v2v_lst[s];
      for (j = s + 1; j < e; j++)
        if (ref_v2v_lst[j-1] != ref_v2v_lst[j])
          ref_v2v_lst[shift++] = ref_v2v_lst[j];
    }

    save = e;
    ref_v2v_idx[i+1] = shift;

  }

  n_ref_edges = ref_v2v_idx[n_vertices];
  BFT_REALLOC(ref_v2v_lst, n_ref_edges, cs_int_t);

#if 0 && defined(DEBUG) && !defined(NDEBUG) /* Dump the structure */
  for (i = 0; i < n_vertices; i++) {
    bft_printf("%9d: (%d, %d) v-v:", i+1, ref_v2v_idx[i],ref_v2v_idx[i+1]);
    for (j = ref_v2v_idx[i]; j < ref_v2v_idx[i+1]; j++)
      bft_printf(" %d ", ref_v2v_lst[j]);
    bft_printf("\n");
    bft_printf_flush();
  }
#endif

  /* Scan adjacent faces to add new "single" vertices */

  BFT_MALLOC(new_s_vertices, max_size, cs_int_t);

  /* Scan adjacent faces */

  for (i = 0; i < selection->n_b_adj_faces; i++) {

    fid = selection->b_adj_faces[i] - 1;
    s = b_f2v_idx[fid] - 1;
    e = b_f2v_idx[fid+1] - 1;

    for (j = s; j < e - 1; j++)
      _add_s_vertex(selection_tag,
                    b_f2v_lst[j]-1,
                    b_f2v_lst[j+1]-1,
                    ref_v2v_idx,
                    ref_v2v_lst,
                    &max_size,
                    &n_new_s_vertices,
                    &new_s_vertices);

    _add_s_vertex(selection_tag,
                  b_f2v_lst[e-1]-1,
                  b_f2v_lst[s]-1,
                  ref_v2v_idx,
                  ref_v2v_lst,
                  &max_size,
                  &n_new_s_vertices,
                  &new_s_vertices);

  }

  for (i = 0; i < selection->n_i_adj_faces; i++) {

    fid = selection->i_adj_faces[i] - 1;
    s = i_f2v_idx[fid] - 1;
    e = i_f2v_idx[fid+1] - 1;

    for (j = s; j < e - 1; j++)
      _add_s_vertex(selection_tag,
                    i_f2v_lst[j]-1,
                    i_f2v_lst[j+1]-1,
                    ref_v2v_idx,
                    ref_v2v_lst,
                    &max_size,
                    &n_new_s_vertices,
                    &new_s_vertices);

    _add_s_vertex(selection_tag,
                  i_f2v_lst[e-1]-1,
                  i_f2v_lst[s]-1,
                  ref_v2v_idx,
                  ref_v2v_lst,
                  &max_size,
                  &n_new_s_vertices,
                  &new_s_vertices);

  }

  BFT_REALLOC(new_s_vertices, n_new_s_vertices, cs_int_t);

  /* Add these vertices to the "single" list */

  BFT_REALLOC(selection->s_vertices,
              selection->n_s_vertices + n_new_s_vertices, cs_int_t);

  for (i = 0; i < n_new_s_vertices; i++)
    selection->s_vertices[selection->n_s_vertices + i] = new_s_vertices[i];

  selection->n_s_vertices += n_new_s_vertices;

  if (selection->n_s_vertices > 1) { /* Delete redundant elements */

    cs_sort_shell(0, selection->n_s_vertices, selection->s_vertices);

    n_new_s_vertices = 1;
    for (i = 1; i < selection->n_s_vertices; i++)
      if (selection->s_vertices[i] != selection->s_vertices[i-1])
        selection->s_vertices[n_new_s_vertices++] = selection->s_vertices[i];

    selection->n_s_vertices = n_new_s_vertices;
    BFT_REALLOC(selection->s_vertices, n_new_s_vertices, cs_int_t);

  }

  /* Free memory */

  BFT_FREE(new_s_vertices);
  BFT_FREE(ref_v2v_idx);
  BFT_FREE(ref_v2v_lst);
  BFT_FREE(count);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Compute the cell center (= barycenter of the cell vertices) for the
 * selected cells.
 *
 * If selection = NULL, the whole mesh is selected.
 *
 * The caller is responsible for freeing the returned array.
 *
 * parameters:
 *   mesh           <-- pointer to a cs_mesh_t structure
 *   n_select_cells <-- number of cells in the selection
 *   selection      <-- list of selected cell number (size: n_cells)
 *                      cell "i" selected if selection[i] > 0
 *                      cell "i" not selected if selection[i] < 0
 *                      0 < selection[i] <= n_select_cells
 *
 * returns:
 *   array of cell centers (size: n_select_cells*3, interlaced)
 *---------------------------------------------------------------------------*/

static cs_real_t *
_cell_cen_vtx(const cs_mesh_t  *mesh,
              cs_int_t          n_select_cells,
              const cs_int_t    selection[])
{
  cs_int_t  i, j, k, coord, cid, cid1, cid2, fid, vid, shift, vtx_counter;
  cs_real_t  denum;

  cs_int_t  *v_tag = NULL, *counter = NULL;
  cs_int_t  *f2v_idx = NULL, *f2v_lst = NULL, *c2f_idx = NULL, *c2f_lst = NULL;

  cs_real_t  *cell_cen = NULL;

  if (selection == NULL)
    n_select_cells = mesh->n_cells;

  if (n_select_cells == 0)
    return NULL;

  BFT_MALLOC(cell_cen, 3*n_select_cells, cs_real_t);

  /* Extract "cell -> faces" connectivity for the selected faces */
  /* Build c2f_idx */

  BFT_MALLOC(c2f_idx, n_select_cells + 1, cs_int_t);

  for (i = 0; i < n_select_cells + 1; i++)
    c2f_idx[i] = 0;

  for (i = 0; i < mesh->n_b_faces; i++) {

    cid = mesh->b_face_cells[i] - 1;

    if (selection != NULL)
      cid = selection[cid];

    if (cid > -1) /* selected cell */
      c2f_idx[cid+1] += 1;

  } /* End of loop on border faces */

  for (i = 0; i < mesh->n_i_faces; i++) {

    cid1 = mesh->i_face_cells[2*i  ] - 1;
    cid2 = mesh->i_face_cells[2*i+1] - 1;

    if (selection != NULL) {

      if (cid1 > -1 && cid1 < mesh->n_cells)
        cid1 = selection[cid1];
      else
        cid1 = -1;

      if (cid2 > -1 && cid2 < mesh->n_cells)
        cid2 = selection[cid2];
      else
        cid2 = -1;

    }

    if (cid1 > -1 && cid1 < mesh->n_cells) /* selected cell */
      c2f_idx[cid1+1] += 1;
    if (cid2 > -1 && cid2 < mesh->n_cells) /* selected cell */
      c2f_idx[cid2+1] += 1;

  } /* End of loop on interior faces */

  c2f_idx[0] = 1;
  for (i = 0; i < n_select_cells; i++)
    c2f_idx[i+1] += c2f_idx[i];

  /* Build c2f_lst */

  BFT_MALLOC(c2f_lst, c2f_idx[n_select_cells]-1, cs_int_t);
  BFT_MALLOC(counter, n_select_cells, cs_int_t);

  for (i = 0; i < n_select_cells; i++)
    counter[i] = 0;

  for (i = 0; i < mesh->n_b_faces; i++) {

    cid = mesh->b_face_cells[i] - 1;

    if (selection != NULL)
      cid = selection[cid];

    if (cid > -1) { /* selected cell */
      shift = c2f_idx[cid] - 1 + counter[cid];
      c2f_lst[shift] = i+1;
      counter[cid] += 1;
    }

  } /* End of loop on border faces */

  for (i = 0; i < mesh->n_i_faces; i++) {

    cid1 = mesh->i_face_cells[2*i  ] - 1;
    cid2 = mesh->i_face_cells[2*i+1] - 1;

    if (selection != NULL) {

      if (cid1 > -1 && cid1 < mesh->n_cells)
        cid1 = selection[cid1];
      else
        cid1 = -1;

      if (cid2 > -1 && cid2 < mesh->n_cells)
        cid2 = selection[cid2];
      else
        cid2 = -1;

    }

    if (cid1 > -1 && cid1 < mesh->n_cells) { /* selected cell */
      shift = c2f_idx[cid1] - 1 + counter[cid1];
      c2f_lst[shift] = i+1 + mesh->n_b_faces;
      counter[cid1] += 1;
    }

    if (cid2 > -1 && cid2 < mesh->n_cells) { /* selected cell */
      shift = c2f_idx[cid2] - 1 + counter[cid2];
      c2f_lst[shift] = i+1 + mesh->n_b_faces;
      counter[cid2] += 1;
    }

  } /* End of loop on interior faces */

  /* Compute cell centers for the selected vertices */

  for (i = 0; i < 3*n_select_cells; i++)
    cell_cen[i] = 0.0;

  BFT_MALLOC(v_tag, mesh->n_vertices, cs_int_t);

  for (vid = 0 ; vid < mesh->n_vertices ; vid++)
    v_tag[vid] = -1;

  /* Loop on selected cells */

  for (i = 0; i < mesh->n_cells; i++) {

    vtx_counter = 0;

    if (selection != NULL)
      cid = selection[i];
    else
      cid = i;

    if (cid > -1) { /* Loop on faces of the cell */

      for (j = c2f_idx[cid] - 1; j < c2f_idx[cid + 1] - 1; j++) {

        fvm_lnum_t  face_num = c2f_lst[j];

        /* Interior or border face */

        if (face_num > mesh->n_b_faces) {
          fid = face_num - mesh->n_b_faces - 1;
          f2v_idx = mesh->i_face_vtx_idx;
          f2v_lst = mesh->i_face_vtx_lst;
        }
        else {
          fid = face_num - 1;
          f2v_idx = mesh->b_face_vtx_idx;
          f2v_lst = mesh->b_face_vtx_lst;
        }

        /* Loop on vertices of the face */

        for (k = f2v_idx[fid] - 1; k < f2v_idx[fid + 1] - 1; k++) {

          vid = f2v_lst[k] - 1;

          if (v_tag[vid] < cid) {
            for (coord = 0 ; coord < 3 ; coord++)
              cell_cen[3*cid + coord] += mesh->vtx_coord[3*vid + coord];
            vtx_counter += 1;
            v_tag[vid] = cid;
          }

        }

      } /* End of loop on faces of the cell */

      denum = 1/(double)vtx_counter;
      for (coord = 0; coord < 3; coord++)
        cell_cen[3*cid + coord] *= denum;

    } /* End if cid >= 0 */

  } /* End of loop on cells */

  /* Free memory */

  BFT_FREE(v_tag);
  BFT_FREE(counter);
  BFT_FREE(c2f_idx);
  BFT_FREE(c2f_lst);

  /* Return result */

  return cell_cen;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Get the full selection of faces related to "single" vertices.
 * Done only if the run is parallel.
 *
 * parameters:
 *   n_b_faces  <-- number of border faces in the parent mesh
 *   b_f2v_idx  <-- border "face -> vertex" connect. index
 *   b_f2v_lst  <-- border "face -> vertex" connect. list
 *   n_i_faces  <-- number of interior faces in the parent mesh
 *   i_f2v_idx  <-- interior "face -> vertex" connect. index
 *   i_f2v_lst  <-- interior "face -> vertex" connect. list
 *   n_vertices <-- number of vertices in the parent mesh
 *   selection  <-> pointer to a fvm_join_selection_t structure
 *---------------------------------------------------------------------------*/

static void
_single_faces_extraction(cs_int_t           n_b_faces,
                         const cs_int_t    *b_f2v_idx,
                         const cs_int_t    *b_f2v_lst,
                         cs_int_t           n_i_faces,
                         const cs_int_t    *i_f2v_idx,
                         const cs_int_t    *i_f2v_lst,
                         cs_int_t           n_vertices,
                         cs_join_select_t  *selection)
{
  int i;

  /* Extract related interior and border single faces */

  _extract_contig_faces(n_vertices,
                        selection->n_s_vertices,
                        selection->s_vertices,
                        n_b_faces,
                        b_f2v_idx,
                        b_f2v_lst,
                        &(selection->n_b_s_faces),
                        &(selection->b_s_faces));

  if (selection->n_b_s_faces > 0) {
    bft_printf(_("\n  Single border faces for the joining operation:\n"));
    for (i = 0; i < selection->n_b_s_faces; i++)
      bft_printf(" %9d | %9d\n", i, selection->b_s_faces[i]);
    bft_printf("\n");
    bft_printf_flush();
  }

  _extract_contig_faces(n_vertices,
                        selection->n_s_vertices,
                        selection->s_vertices,
                        n_i_faces,
                        i_f2v_idx,
                        i_f2v_lst,
                        &(selection->n_i_s_faces),
                        &(selection->i_s_faces));

  if (selection->n_i_s_faces > 0) {
    bft_printf(_("\n  Single interior faces for the joining operation:\n"));
    for (i = 0; i < selection->n_i_s_faces; i++)
      bft_printf(" %9d | %9d\n", i, selection->i_s_faces[i]);
    bft_printf("\n");
    bft_printf_flush();
  }

}

/*----------------------------------------------------------------------------
 * Get the full selection of vertices to extract.
 *
 * Only reaquired in parallel, as some vertices may have been selected on
 * another rank and not on the local rank, but we need to take them into
 * account to have a good update of the mesh.
 *
 * parameters:
 *   b_f2v_idx    <-- border "face -> vertex" connect. index
 *   b_f2v_lst    <-- border "face -> vertex" connect. list
 *   i_f2v_idx    <-- interior "face -> vertex" connect. index
 *   i_f2v_lst    <-- interior "face -> vertex" connect. list
 *   n_vertices   <-- number of vertices in the parent mesh
 *   n_g_vertices <-- global number of vertices in the parent mesh
 *   v_gnum       <-- global vertex numbering (NULL if n_ranks = 1)
 *   join_select  <-> pointer to a cs_join_select_t structure
 *---------------------------------------------------------------------------*/

static void
_single_vertex_extraction(const cs_int_t     b_f2v_idx[],
                          const cs_int_t     b_f2v_lst[],
                          const cs_int_t     i_f2v_idx[],
                          const cs_int_t     i_f2v_lst[],
                          cs_int_t           n_vertices,
                          fvm_gnum_t         n_g_vertices,
                          const fvm_gnum_t   v_gnum[],
                          cs_join_select_t  *selection)
{
  cs_int_t  i, id, rank, vid, shift;
  fvm_gnum_t  n_g, gnum;

  int  _have_single_elts = 0, have_single_elts = 0;
  cs_int_t  n_select_vertices = 0;
  cs_int_t  *block_tag = NULL;
  int  *send_shift = NULL, *recv_shift = NULL;
  int  *send_count = NULL, *recv_count = NULL;
  fvm_gnum_t  *send_glist = NULL, *recv_glist = NULL;

  MPI_Comm  mpi_comm = cs_glob_mpi_comm;

  const int  loc_rank = CS_MAX(cs_glob_rank_id, 0);
  const int  n_ranks = cs_glob_n_ranks;
  const cs_join_block_info_t  block = cs_join_get_block_info(n_g_vertices,
                                                             n_ranks,
                                                             loc_rank);

  assert(v_gnum != NULL);
  assert(n_ranks > 1);

  /* Define a global tag on vertices by block */

  BFT_MALLOC(block_tag, block.local_size, cs_int_t);

  for (i = 0; i < (cs_int_t)block.local_size; i++)
    block_tag[i] = 0; /* Not selected by default */

  BFT_MALLOC(send_count, n_ranks, int);
  BFT_MALLOC(recv_count, n_ranks, int);

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (i = 0; i < selection->n_vertices; i++) {
    vid = selection->vertices[i] - 1;
    rank = (v_gnum[vid] - 1)/block.size;
    send_count[rank] += 1;
  }

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);

  BFT_MALLOC(send_shift, n_ranks + 1, int);
  BFT_MALLOC(recv_shift, n_ranks + 1, int);

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (rank = 0; rank < n_ranks; rank++) {
    send_shift[rank + 1] = send_shift[rank] + send_count[rank];
    recv_shift[rank + 1] = recv_shift[rank] + recv_count[rank];
  }

  assert(send_shift[n_ranks] == selection->n_vertices);

  /* Build send_list */

  BFT_MALLOC(send_glist, send_shift[n_ranks], fvm_gnum_t);
  BFT_MALLOC(recv_glist, recv_shift[n_ranks], fvm_gnum_t);

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (i = 0; i < selection->n_vertices; i++) {

    vid = selection->vertices[i] - 1;
    rank = (v_gnum[vid] - 1)/block.size;
    shift = send_shift[rank] + send_count[rank];
    send_glist[shift] = v_gnum[vid]; /* Global number of the selected vertex */
    send_count[rank] += 1;

  }

  MPI_Alltoallv(send_glist, send_count, send_shift, FVM_MPI_GNUM,
                recv_glist, recv_count, recv_shift, FVM_MPI_GNUM,
                mpi_comm);

  /* Get a synchronized tag on each block */

  for (rank = 0; rank < n_ranks; rank++) {
    for (i = recv_shift[rank]; i < recv_shift[rank+1]; i++) {

      gnum = recv_glist[i];
      id = gnum - block.first_gnum;
      block_tag[id] = 1;

    }
  } /* End of loop on ranks */

  /* Update local vertex selection from the global sync. tag array */

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (i = 0; i < n_vertices; i++) {
    rank = (v_gnum[i] - 1)/block.size;
    send_count[rank] += 1;
  }

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (rank = 0; rank < n_ranks; rank++) {
    send_shift[rank + 1] = send_shift[rank] + send_count[rank];
    recv_shift[rank + 1] = recv_shift[rank] + recv_count[rank];
  }

  assert(send_shift[n_ranks] == n_vertices);

  /* Build send_list */

  BFT_REALLOC(send_glist, send_shift[n_ranks], fvm_gnum_t);
  BFT_REALLOC(recv_glist, recv_shift[n_ranks], fvm_gnum_t);

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (i = 0; i < n_vertices; i++) {
    rank = (v_gnum[i] - 1)/block.size;
    shift = send_shift[rank] + send_count[rank];
    send_glist[shift] = v_gnum[i]; /* Global number of the vertex */
    send_count[rank] += 1;
  }

  MPI_Alltoallv(send_glist, send_count, send_shift, FVM_MPI_GNUM,
                recv_glist, recv_count, recv_shift, FVM_MPI_GNUM,
                mpi_comm);

  /* Use sync. tag by block to send back to the rank if their vertices
     are selected */

  for (i = 0; i < recv_shift[n_ranks]; i++) {
    gnum = recv_glist[i];
    id = gnum - block.first_gnum;
    recv_glist[i] = (fvm_gnum_t)block_tag[id];
  }

  MPI_Alltoallv(recv_glist, recv_count, recv_shift, FVM_MPI_GNUM,
                send_glist, send_count, send_shift, FVM_MPI_GNUM,
                mpi_comm);

  /* Free memory */

  BFT_FREE(block_tag);
  BFT_FREE(send_count);
  BFT_FREE(send_shift);
  BFT_FREE(recv_glist);

  for (i = 0; i < n_vertices; i++)
    if (send_glist[i] > 0)
      n_select_vertices += 1;

  if (n_select_vertices != selection->n_vertices)
    _have_single_elts = 1;

  MPI_Allreduce(&_have_single_elts, &have_single_elts, 1, MPI_INT, MPI_MAX, mpi_comm);

  if (have_single_elts == 1) {

    cs_int_t  *loc_tag = NULL;

    bft_printf("\n  Synchronization of the vertex selection necessary.\n");

    if (n_select_vertices != selection->n_vertices) {

      assert(n_select_vertices > selection->n_vertices);

      BFT_MALLOC(loc_tag, n_vertices, cs_int_t);

      for (i = 0; i < n_vertices; i++)
        loc_tag[i] = 0;

      for (i = 0; i < selection->n_vertices; i++)
        loc_tag[selection->vertices[i] - 1] = 1;

      selection->n_s_vertices = n_select_vertices - selection->n_vertices;
      BFT_MALLOC(selection->s_vertices, selection->n_s_vertices, cs_int_t);

      shift = 0;
      for (i = 0; i < n_vertices; i++)
        if (send_glist[i] > 0 && loc_tag[i] == 0)
          selection->s_vertices[shift++] = i+1;

      assert(shift == selection->n_s_vertices);

      BFT_FREE(loc_tag);

    } /* End if n_select_vertices != selection->n_vertices */

    _add_s_vertices_from_adj_faces(n_vertices,
                                   send_glist,
                                   selection,
                                   b_f2v_idx,
                                   b_f2v_lst,
                                   i_f2v_idx,
                                   i_f2v_lst);

    if (selection->n_s_vertices > 0) {

      bft_printf(_("\n  Single vertices for the joining operation:\n"));
      for (i = 0; i < selection->n_s_vertices; i++)
        bft_printf(" %9d | %9d\n", i, selection->s_vertices[i]);
      bft_printf("\n");
      bft_printf_flush();

    } /* End if selection->n_s_vertices > 0 */

    /* Check globally if there are single vertices */

    MPI_Allreduce(&(selection->n_s_vertices), &n_g, 1, FVM_MPI_GNUM,
                  MPI_SUM, mpi_comm);

    if (n_g > 0) {

      cs_int_t  j, k, ifs_size;

      cs_int_t  *count = NULL;
      fvm_gnum_t  *loc_s_v_gnum = NULL, *glob_s_v_gnum = NULL;
      fvm_interface_set_t  *ifs = NULL;

      selection->do_single_sync = true;

      /* NB: n_s_vertices and n_c_vertices values are assumed to be small */

      /* Define couple vertices (vertex related on a distant to a single
         vertex) */

      MPI_Allgather(&(selection->n_s_vertices), 1, MPI_INT,
                    recv_count,                 1, MPI_INT, mpi_comm);

      recv_shift[0] = 0;
      for (i = 0; i < n_ranks; i++)
        recv_shift[i+1] = recv_shift[i] + recv_count[i];

      assert(recv_shift[n_ranks] == (cs_int_t)n_g);

      BFT_MALLOC(glob_s_v_gnum, recv_shift[n_ranks], fvm_gnum_t);
      BFT_MALLOC(loc_s_v_gnum, selection->n_s_vertices, fvm_gnum_t);

      for (i = 0; i < selection->n_s_vertices; i++)
        loc_s_v_gnum[i] = v_gnum[selection->s_vertices[i]-1];

      MPI_Allgatherv(loc_s_v_gnum, selection->n_s_vertices, FVM_MPI_GNUM,
                     glob_s_v_gnum, recv_count, recv_shift, FVM_MPI_GNUM,
                     mpi_comm);

      /* Count the number of coupled vertices */

      for (rank = 0; rank < n_ranks; rank++) {
        if (rank != loc_rank) {

          for (i = recv_shift[rank]; i < recv_shift[rank+1]; i++) {

            id = cs_search_g_binary(n_vertices ,
                                    glob_s_v_gnum[i],
                                    v_gnum);

            if (id != -1) {

              _Bool  in_s_list = false;

              for (j = 0; j < selection->n_s_vertices; j++) {
                if (selection->s_vertices[j] == id + 1) {
                  in_s_list = true;
                  break;
                }
              }

              if (in_s_list == false)
                selection->n_c_vertices += 1;

            } /* id != -1 */

          }

        } /* rank != loc_rank */
      } /* End of loop on ranks */

      BFT_MALLOC(selection->c_vertices, selection->n_c_vertices, cs_int_t);
      BFT_MALLOC(selection->c_vtx_rank_lst, selection->n_c_vertices, cs_int_t);

      /* Defined coupled vertices list and its related rank */

      for (shift = 0, rank = 0; rank < n_ranks; rank++) {
        if (rank != loc_rank) {

          for (i = recv_shift[rank]; i < recv_shift[rank+1]; i++) {

            id = cs_search_g_binary(n_vertices,
                                    glob_s_v_gnum[i],
                                    v_gnum);

            if (id != -1) {

              _Bool  in_s_list = false;

              for (j = 0; j < selection->n_s_vertices; j++) {
                if (selection->s_vertices[j] == id + 1) {
                  in_s_list = true;
                  break;
                }
              }

              if (in_s_list == false) {
                selection->c_vertices[shift] = id + 1;
                selection->c_vtx_rank_lst[shift] = rank;
                shift++;
              }

            } /* id != -1 */

          }

        }  /* rank != loc_rank */

      } /* End of loop on ranks */

      assert(shift == selection->n_c_vertices);

      if (selection->n_c_vertices > 0) {

        bft_printf(_("\n  Coupled vertices for the joining operations:\n"));
        for (i = 0; i < selection->n_c_vertices; i++)
          bft_printf(" %9d | %9d | %6d\n",
                     i, selection->c_vertices[i], selection->c_vtx_rank_lst[i]);
        bft_printf("\n");
        bft_printf_flush();

      }

      /* Build index and list of distant vertices related to single
         vertices */

      ifs = fvm_interface_set_create(n_vertices,
                                     NULL,
                                     v_gnum,
                                     NULL,
                                     0,
                                     NULL,
                                     NULL,
                                     NULL);

      /* Define index for the single vertices */

      BFT_MALLOC(selection->s_vtx_idx, selection->n_s_vertices + 1, cs_int_t);

      for (i = 0; i < selection->n_s_vertices + 1; i++)
        selection->s_vtx_idx[i] = 0;

      ifs_size = fvm_interface_set_size(ifs);

      for (i = 0; i < ifs_size; i++) {

        const fvm_interface_t  *itf = fvm_interface_set_get(ifs, i);
        const fvm_lnum_t  *loc_num = fvm_interface_get_local_num(itf);
        int  itf_size = fvm_interface_size(itf);
        int  dist_rank = fvm_interface_rank(itf);

        for (j = 0; j < selection->n_s_vertices; j++) {

          gnum = loc_s_v_gnum[j];
          id = cs_search_binary(itf_size,
                                selection->s_vertices[j],
                                loc_num);

          if (id != -1) {

            _Bool  in_s_list = false;

            for (k = recv_shift[dist_rank]; k < recv_shift[dist_rank+1]; k++) {
              if (gnum == glob_s_v_gnum[k]) {
                in_s_list = true;
                break;
              }
            }

            if (in_s_list == false)
              selection->s_vtx_idx[j+1] += 1;

          } /* Id != -1 */

        } /* End of loop on single vertices */

      } /* End of loop on interfaces */

      BFT_MALLOC(count, selection->n_s_vertices, int);

      for (i = 0; i < selection->n_s_vertices; i++) {
        selection->s_vtx_idx[i+1] += selection->s_vtx_idx[i];
        count[i] = 0;
      }

      BFT_MALLOC(selection->s_vtx_rank_lst,
                 selection->s_vtx_idx[selection->n_s_vertices], cs_int_t);

      for (i = 0; i < ifs_size; i++) {

        const fvm_interface_t  *itf = fvm_interface_set_get(ifs, i);
        const fvm_lnum_t  *loc_num = fvm_interface_get_local_num(itf);
        int  itf_size = fvm_interface_size(itf);
        int  dist_rank = fvm_interface_rank(itf);

        for (j = 0; j < selection->n_s_vertices; j++) {

          gnum = loc_s_v_gnum[j];
          id = cs_search_binary(itf_size,
                                selection->s_vertices[j],
                                loc_num);

          if (id != -1) {

            _Bool  in_s_list = false;

            for (k = recv_shift[dist_rank]; k < recv_shift[dist_rank+1]; k++) {
              if (gnum == glob_s_v_gnum[k]) {
                in_s_list = true;
                break;
              }
            }

            if (in_s_list == false) {
              shift = selection->s_vtx_idx[j] + count[j];
              selection->s_vtx_rank_lst[shift] = dist_rank;
              count[j] += 1;
            }

          } /* Id != -1 */

        } /* End of loop on single vertices */

      } /* End of loop on interfaces */

      /* Free memory */

      BFT_FREE(loc_s_v_gnum);
      BFT_FREE(glob_s_v_gnum);
      BFT_FREE(count);

      ifs = fvm_interface_set_destroy(ifs);

    } /* End if selection->n_g_s_vertices > 0 */

  } /* End if have_single_elements */

  /* Free memory */

  BFT_FREE(send_glist);
  BFT_FREE(recv_count);
  BFT_FREE(recv_shift);
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Define the cell center for a selection of cells (those one which are
 * adjacent to the selected faces).
 *
 * parameters:
 *   n_select_faces  <-- number of border faces implied in the joining op.
 *   select_faces    <-- list of faces implied in the joining op.
 *   mesh            <-- pointer to cs_mesh_t structure
 *   cell_filter     --> selection array (size: mesh->n_cells)
 *   select_cell_cen --> cell center for the selected cells
 *---------------------------------------------------------------------------*/

static void
_get_select_cell_cen(cs_int_t          n_select_faces,
                     const cs_int_t    select_faces[],
                     const cs_mesh_t  *mesh,
                     cs_int_t         *cell_filter[],
                     cs_real_t        *select_cell_cen[])
{
  cs_int_t  i, cid, fid, tag;

  cs_int_t  n_select_cells = 0;
  cs_int_t  *_cell_filter = NULL;
  cs_real_t  *_select_cell_cen = NULL;

  BFT_MALLOC(_cell_filter, mesh->n_cells, cs_int_t);

  for (i = 0; i < mesh->n_cells; i++)
    _cell_filter[i] = -1;

  for (i = 0; i < n_select_faces; i++) {
    fid = select_faces[i] - 1;
    cid = mesh->b_face_cells[fid] - 1;
    _cell_filter[cid] = i + 1;
  }

  /* Order cell selection */

  for (i = 0, tag = 0; i < mesh->n_cells; i++)
    if (_cell_filter[i] > 0)
      _cell_filter[i] = tag++;

  n_select_cells = n_select_faces; /* Border faces: 1 Face -> 1 Cell */

  _select_cell_cen = _cell_cen_vtx(mesh, n_select_cells, _cell_filter);

  /* Returns pointers */

  *cell_filter = _cell_filter;
  *select_cell_cen = _select_cell_cen;
}

/*============================================================================
 * Public function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Define a set of parameters to control a contiguous distribution by block.
 *
 * parameters:
 *   n_g_elts   <-- global number of elements to treat
 *   n_ranks    <-- number of ranks in the MPI communicator related to the
 *                  cs_join_block_info_t structure to create
 *   local_rank <-- rank in the MPI communicator related to the
 *                  cs_join_block_info_t structure to create
 *
 * returns:
 *   a new defined cs_join_block_info_t structure
 *---------------------------------------------------------------------------*/

cs_join_block_info_t
cs_join_get_block_info(fvm_gnum_t  n_g_elts,
                       int         n_ranks,
                       int         local_rank)
{
  size_t  block_size, n_treat_elts;
  fvm_gnum_t  first_glob_num, last_glob_num;

  cs_join_block_info_t  block_info;

  /* Compute block_size */

  block_size = n_g_elts/n_ranks;
  if (n_g_elts % n_ranks > 0)
    block_size++;

  first_glob_num = local_rank * block_size + 1;
  last_glob_num = first_glob_num + block_size;

  if (last_glob_num > n_g_elts) {

    if (first_glob_num > n_g_elts) /* n_ranks > n_g_elts */
      n_treat_elts = 0;
    else
      n_treat_elts = n_g_elts - first_glob_num + 1;

  }
  else
    n_treat_elts = block_size;

  block_info.n_g_elts = n_g_elts;
  block_info.first_gnum = first_glob_num;

  block_info.size = block_size;
  block_info.local_size = n_treat_elts;

  block_info.n_ranks = n_ranks;
  block_info.local_rank = local_rank;

  return  block_info;
}

/*----------------------------------------------------------------------------
 * Initialize a cs_join_param_t structure.
 *
 * parameters:
 *   join_id       <-- id of the current joining operation
 *   fraction      <-- value of the fraction parameter
 *   plane         <-- value of the plane parameter
 *   rtf           <-- value of the "reduction tolerance factor" parameter
 *   ftf           <-- value of the "merge tolerance factor" parameter
 *   etf           <-- value of the "edge equiv. tolerance factor" parameter
 *   max_sub_faces <-- maximum number of sub-faces allowed during splitting
 *   tml           <-- value of the "tree max level" parameter
 *   tmb           <-- value of the "tree max boxes" parameter
 *   tmr           <-- value of the "tree max ratio" parameter
 *   verbosity     <-- level of verbosity required
 *
 * returns:
 *   a pointer to a cs_join_param_t structure
 *---------------------------------------------------------------------------*/

cs_join_param_t
cs_join_param_define(int     join_id,
                     double  fraction,
                     double  plane,
                     double  rtf,
                     double  ftf,
                     double  etf,
                     int     max_sub_faces,
                     int     tml,
                     int     tmb,
                     double  tmr,
                     int     verbosity)
{
  cs_join_param_t  param;

  param.num = join_id + 1;

  /* geometric parameters */

  /* parameter used to compute the tolerance associated to each vertex.
     Also used for finding equivalent vertices during edge intersections */

  param.fraction = fraction;

  /* parameter used to judge if two faces are in the same plane (during
     the face splitting) */

  param.plane = plane;

  /* Coef. used to reduce the tolerance during merge step:
     new tol. = tol * coef.
     Values between [0.0, 1.0[ */

  param.reduce_tol_factor = rtf;

  /* Coef. used to modify the tolerance associated to each vertex BEFORE the
     merge operation.
     If coef = 0.0 => no vertex merge
     If coef < 1.0 => reduce vertex merge
     If coef = 1.0 => no change
     If coef > 1.0 => increase vertex merge */

  param.merge_tol_coef = ftf;

  /* Coef. used to modify locally the tolerance associated to each vertex
     BEFORE adding equivalences between vertices after edge intersections.
     If coef = 0.0 => add no equivalence
     If coef < 1.0 => reduce the number of equivalences between vertices
                      sharing the same edge
     If coef = 1.0 => no change
     If coef > 1.0 => increase the number of equivalences between vertices
                      sharing the same edge. Not advised. */

   param.edge_equiv_tol_coef = etf;

  /* Parameter to switch on/off the influence of adjacent faces in the
     computation of tolerance */

   param.include_adj_faces = true;  /* Default value: true */

  /* Parameter used to define if we get vertex equivalence trough the
     comparison of vertex tolerance or through the difference of curvilinear
     abscissa of vertices on edges.
     If include_adj_faces = false => this parameter should not have any
     effect. (Not a user-defined parameter) */

   param.edge_equiv_by_tolerance = true; /* Default value: true */

   /* Maximum number of sub-faces */

   param.max_sub_faces = max_sub_faces;

   /* Deepest level reachable during tree building */

   param.tree_max_level = tml;

   /* Max. number of boxes which can be related to a leaf of the tree
      if level != tree_max_level */

   param.tree_n_max_boxes = tmb;

   /* Stop tree building if: n_linked_boxes > tree_max_box_ratio*n_init_boxes */

   param.tree_max_box_ratio = tmr;

   /* Level of display */

   param.verbosity = verbosity;

   return param;
}

/*----------------------------------------------------------------------------
 * Create and initialize a cs_join_select_t structure.
 *
 * parameters:
 *   selection_criteria <-- pointer to a cs_mesh_select_t structure
 *   verbosity          <-- level of verbosity required
 *
 * returns:
 *   pointer to a newly created cs_join_select_t structure
 *---------------------------------------------------------------------------*/

cs_join_select_t *
cs_join_select_create(const char  *selection_criteria,
                      int          verbosity)
{
  cs_int_t  i, fid, cid, shift, n_faces;

  cs_join_select_t  *selection = NULL;
  cs_mesh_t  *mesh = cs_glob_mesh;

  const int  n_ranks = cs_glob_n_ranks;

  assert(mesh != NULL);

  /* Initialize fvm_join_selection_t struct. */

  BFT_MALLOC(selection, 1, cs_join_select_t);

  selection->n_faces = 0;
  selection->n_g_faces = 0;

  selection->faces = NULL;
  selection->compact_face_gnum = NULL;
  selection->compact_rank_index = NULL;

  selection->cell_filter = NULL;
  selection->cell_cen = NULL;
  selection->cell_gnum = NULL;

  selection->n_vertices = 0;
  selection->n_g_vertices = 0;
  selection->vertices = NULL;

  selection->n_b_adj_faces = 0;
  selection->b_adj_faces = NULL;
  selection->n_i_adj_faces = 0;
  selection->i_adj_faces = NULL;

  /*
     Single elements (Only possible in parallel. It appears
     when the domain splitting has a poor quality and elements
     on the joining interface are prisms or tetraedrals)
     s = single / c = coupled
  */

  selection->do_single_sync = false;
  selection->n_s_vertices = 0;
  selection->s_vertices = NULL;
  selection->s_vtx_idx = NULL;
  selection->s_vtx_rank_lst = NULL;

  selection->n_c_vertices = 0;
  selection->c_vertices = NULL;
  selection->c_vtx_rank_lst = NULL;

  selection->n_b_s_faces = 0;
  selection->b_s_faces = NULL;
  selection->n_i_s_faces = 0;
  selection->i_s_faces = NULL;

  /* Extract selected border faces */

  BFT_MALLOC(selection->faces, mesh->n_b_faces, fvm_lnum_t);

  cs_selector_get_b_face_list(selection_criteria,
                              &(selection->n_faces),
                              selection->faces);

  BFT_REALLOC(selection->faces, selection->n_faces, fvm_lnum_t);

  /* Define cell_gnum: global numbers of the related cells */

  BFT_MALLOC(selection->cell_gnum, selection->n_faces, fvm_gnum_t);

  if (n_ranks == 1) { /* Serial treatment */

    selection->n_g_faces = selection->n_faces;

    for (i = 0; i < selection->n_faces; i++) {
      fid = selection->faces[i] - 1;
      selection->cell_gnum[i] = mesh->b_face_cells[fid];
    }

  }

#if defined(HAVE_MPI)

  if (n_ranks > 1) { /* Parallel treatment */

    MPI_Allreduce(&(selection->n_faces), &(selection->n_g_faces),
                  1, FVM_MPI_GNUM, MPI_SUM, cs_glob_mpi_comm);

    for (i = 0; i < selection->n_faces; i++) {
      fid = selection->faces[i] - 1;
      cid = mesh->b_face_cells[fid] - 1;
      selection->cell_gnum[i] = mesh->global_cell_num[cid];
    }

  }

#endif

  /* Define the cell center for a selection of cells (those one which are
     adjacent to the selected faces).
     This will be used to check the orientation of faces at the end of the
     joining operation */

  _get_select_cell_cen(selection->n_faces,
                       selection->faces,
                       mesh,
                       &(selection->cell_filter),
                       &(selection->cell_cen));

  if (verbosity > 0)
    bft_printf(_("  Number of boundary faces selected for joining: %10u\n"),
               selection->n_g_faces);

  /* Define a compact global numbering on selected border faces and
     build an index on ranks on this compact numbering */

  _compact_face_gnum_selection(selection->n_faces,
                               &(selection->compact_face_gnum),
                               &(selection->compact_rank_index));

  assert(selection->n_g_faces == selection->compact_rank_index[n_ranks]);

  /* Extract selected vertices from the selected border faces */

  _extract_vertices(selection->n_faces,
                    selection->faces,
                    mesh->b_face_vtx_idx,
                    mesh->b_face_vtx_lst,
                    mesh->n_vertices,
                    &(selection->n_vertices),
                    &(selection->vertices));

  /* Extract list of border faces contiguous to the selected vertices  */

  _extract_contig_faces(mesh->n_vertices,
                        selection->n_vertices,
                        selection->vertices,
                        mesh->n_b_faces,
                        mesh->b_face_vtx_idx,
                        mesh->b_face_vtx_lst,
                        &(selection->n_b_adj_faces),
                        &(selection->b_adj_faces));

  /* Remove border faces already defined in selection->faces */

  for (i = 0, shift = 0, n_faces = 0;
       i < selection->n_b_adj_faces && shift < selection->n_faces; i++) {
    if (selection->b_adj_faces[i] != selection->faces[shift])
      selection->b_adj_faces[n_faces++] = selection->b_adj_faces[i];
    else
      shift++;
  }

  for (;i < selection->n_b_adj_faces; i++, n_faces++)
    selection->b_adj_faces[n_faces] = selection->b_adj_faces[i];

  selection->n_b_adj_faces = n_faces;
  BFT_REALLOC(selection->b_adj_faces, n_faces, cs_int_t);

  /* Extract list of interior faces contiguous to the selected vertices */

  _extract_contig_faces(mesh->n_vertices,
                        selection->n_vertices,
                        selection->vertices,
                        mesh->n_i_faces,
                        mesh->i_face_vtx_idx,
                        mesh->i_face_vtx_lst,
                        &(selection->n_i_adj_faces),
                        &(selection->i_adj_faces));

   /* Check if there is no forgotten vertex in the selection.
      Otherwise define structures to enable future synchronization.
      Only possible in parallel. */

#if defined(HAVE_MPI)

  if (n_ranks > 1) {

    _single_vertex_extraction(mesh->b_face_vtx_idx,
                              mesh->b_face_vtx_lst,
                              mesh->i_face_vtx_idx,
                              mesh->i_face_vtx_lst,
                              mesh->n_vertices,
                              mesh->n_g_vertices,
                              mesh->global_vtx_num,
                              selection);

    if (selection->n_s_vertices > 0)
      _single_faces_extraction(mesh->n_b_faces,
                               mesh->b_face_vtx_idx,
                               mesh->b_face_vtx_lst,
                               mesh->n_i_faces,
                               mesh->i_face_vtx_idx,
                               mesh->i_face_vtx_lst,
                               mesh->n_vertices,
                               selection);

  }

#endif

  /* Display information according to the level of verbosity */

  if (verbosity > 1) {

    bft_printf("\n  Compact index on ranks for the selected faces:\n");
    for (i = 0; i < n_ranks + 1; i++)
      bft_printf(" %5d | %11u\n", i, selection->compact_rank_index[i]);
    bft_printf("\n");

    if (verbosity > 2) {

      bft_printf(_("\n  Selected faces for the joining operation:\n"));
      for (i = 0; i < selection->n_faces; i++)
        bft_printf(" %9d | %9d | %10u | %10u\n",
                   i, selection->faces[i], selection->compact_face_gnum[i],
                   selection->cell_gnum[i]);
      bft_printf("\n");

      bft_printf(_("\n  Select vertices for the joining operation:\n"));
      for (i = 0; i < selection->n_vertices; i++)
        bft_printf(" %9d | %9d\n", i, selection->vertices[i]);
      bft_printf("\n");

      if (verbosity > 3) {

        bft_printf
          (_("\n  Contiguous border faces for the joining operation:\n"));
        for (i = 0; i < selection->n_b_adj_faces; i++)
          bft_printf(" %9d | %9d\n", i, selection->b_adj_faces[i]);
        bft_printf("\n");

        bft_printf
          (_("\n  Contiguous interior faces for the joining operation:\n"));
        for (i = 0; i < selection->n_i_adj_faces; i++)
          bft_printf(" %9d | %9d\n", i, selection->i_adj_faces[i]);
        bft_printf("\n");

      } /* End if verbosity > 3 */

    } /* End if verbosity > 2 */

    bft_printf_flush();

  } /* End if verbosity > 1 */

  return  selection;
}

/*----------------------------------------------------------------------------
 * Destroy a cs_join_select_t structure.
 *
 * parameters:
 *   join_select <-- pointer to pointer to structure to destroy
 *---------------------------------------------------------------------------*/

void
cs_join_select_destroy(cs_join_select_t  **join_select)
{
  if (*join_select != NULL) {

    cs_join_select_t *_js = *join_select;

    BFT_FREE(_js->faces);
    BFT_FREE(_js->compact_face_gnum);
    BFT_FREE(_js->cell_filter);
    BFT_FREE(_js->cell_cen);
    BFT_FREE(_js->cell_gnum);
    BFT_FREE(_js->compact_rank_index);
    BFT_FREE(_js->vertices);
    BFT_FREE(_js->b_adj_faces);
    BFT_FREE(_js->i_adj_faces);

    if (_js->do_single_sync > 0) {

      BFT_FREE(_js->s_vertices);
      BFT_FREE(_js->s_vtx_idx);
      BFT_FREE(_js->s_vtx_rank_lst);
      BFT_FREE(_js->c_vertices);
      BFT_FREE(_js->c_vtx_rank_lst);
      BFT_FREE(_js->b_s_faces);
      BFT_FREE(_js->i_s_faces);

    }

    BFT_FREE(*join_select);
  }
}

/*----------------------------------------------------------------------------
 * Build vertex -> vertex index for a selection of faces.
 *
 * "v2v_idx" is already allocated to the number of vertices in the mesh.
 * At this stage, it is just a counter.
 *
 * parameters:
 *   n_faces <-- number of selected faces
 *   faces   <-- list of selected faces
 *   f2v_idx <-- face -> vertex connectivity index
 *   f2v_lst <-- face -> vertex connectivity list
 *   v2v_idx <-> index to build (already allocated and may be used again)
 *---------------------------------------------------------------------------*/

void
cs_join_build_edges_idx(cs_int_t        n_faces,
                        const cs_int_t  faces[],
                        const cs_int_t  f2v_idx[],
                        const cs_int_t  f2v_lst[],
                        cs_int_t        v2v_idx[])
{
  cs_int_t  i, j, v1, v2, fid, s, e;

  /* Loop on all selected faces. No need to loop on other faces because
     the selected vertices are all found with this only step. */

  for (i = 0; i < n_faces; i++) {

    fid = faces[i] - 1;
    s = f2v_idx[fid] - 1;
    e = f2v_idx[fid+1] - 1;

    for (j = s; j < e - 1; j++) { /* scan edges */

      v1 = f2v_lst[j];
      v2 = f2v_lst[j+1];

      if (v1 < v2)
        v2v_idx[v1] += 1;
      else if (v2 < v1)
        v2v_idx[v2] += 1;
      else
        bft_error(__FILE__, __LINE__, 0,
                  _("  Inconsistent mesh definition. Cannot build edges.\n"
                    "  Face %d has the same vertex %d twice.\n"), fid+1, v1);

    }

    /* Last edge */

    v1 = f2v_lst[e-1];
    v2 = f2v_lst[s];

    if (v1 < v2)
      v2v_idx[v1] += 1;
    else if (v2 < v1)
      v2v_idx[v2] += 1;
    else
      bft_error(__FILE__, __LINE__, 0,
                _("  Inconsistent mesh definition. Cannot build edges.\n"
                  "  Face %d has the same vertex %d twice.\n"), fid+1, v1);

  } /* End of loop on selected faces */
}

/*----------------------------------------------------------------------------
 * Build vertex -> vertex list for a selection of faces.
 * "count" and "v2v_lst" are already allocated to the number of vertices in
 * the mesh.
 *
 * parameters:
 *   n_faces <-- number of selected faces
 *   faces   <-- list of selected faces
 *   f2v_idx <-- face -> vertex connectivity index
 *   f2v_lst <-- face -> vertex connectivity list
 *   count   <-> array used to count the number of values already added
 *   v2v_idx <-- vertex -> vertex connect. index
 *   v2v_lst <-> vertex -> vertex connect. list to build (can be used again)
 *---------------------------------------------------------------------------*/

void
cs_join_build_edges_lst(cs_int_t        n_faces,
                        const cs_int_t  faces[],
                        const cs_int_t  f2v_idx[],
                        const cs_int_t  f2v_lst[],
                        cs_int_t        count[],
                        const cs_int_t  v2v_idx[],
                        cs_int_t        v2v_lst[])
{
  cs_int_t  i, j, v1_id, v2_id, fid, s, e, shift;

  for (i = 0; i < n_faces; i++) {

    fid = faces[i] - 1;
    s = f2v_idx[fid] - 1;
    e = f2v_idx[fid+1] - 1;

    for (j = s; j < e - 1; j++) { /* Scan edges */

      v1_id = f2v_lst[j] - 1;
      v2_id = f2v_lst[j+1] - 1;

      if (v1_id < v2_id) {
        shift = v2v_idx[v1_id] + count[v1_id];
        v2v_lst[shift] = v2_id + 1;
        count[v1_id] += 1;
      }
      else if (v2_id < v1_id) {
        shift = v2v_idx[v2_id] + count[v2_id];
        v2v_lst[shift] = v1_id + 1;
        count[v2_id] += 1;
      }

    }

    /* Last edge */

    v1_id = f2v_lst[e-1] - 1;
    v2_id = f2v_lst[s] - 1;

    if (v1_id < v2_id) {
      shift = v2v_idx[v1_id] + count[v1_id];
      v2v_lst[shift] = v2_id + 1;
      count[v1_id] += 1;
    }
    else if (v2_id < v1_id) {
      shift = v2v_idx[v2_id] + count[v2_id];
      v2v_lst[shift] = v1_id + 1;
      count[v2_id] += 1;
    }

  } /* End of loop on selected faces */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
