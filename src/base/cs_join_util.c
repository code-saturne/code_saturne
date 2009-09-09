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
 * Eliminate redundancies found between two lists of elements.
 * Delete elements in elts[] and keep elements in the reference list.
 *
 * parameters:
 *  n_elts      <->  number of elements in the list to clean
 *  elts        <->  list of elements in the list to clean
 *  n_ref_elts  -->  number of elements in the reference list
 *  ref_elts    -->  list of reference elements
 *---------------------------------------------------------------------------*/

static void
_clean_selection(cs_int_t   *n_elts,
                 cs_int_t   *elts[],
                 cs_int_t    n_ref_elts,
                 cs_int_t    ref_elts[])
{
  cs_int_t  i = 0, j = 0;
  cs_int_t  _n_elts = 0;
  cs_int_t  *_elts = *elts;

  while (i < *n_elts && j < n_ref_elts) {

    if (_elts[i] < ref_elts[j])
      _elts[_n_elts++] = _elts[i++];
    else if (_elts[i] > ref_elts[j])
      j++;
    else
      i++, j++;

  }

  for (;i < *n_elts; i++, _n_elts++)
    _elts[_n_elts] = _elts[i];

  BFT_REALLOC(_elts, _n_elts, cs_int_t);

  *n_elts = _n_elts;
  *elts = _elts;
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

/*----------------------------------------------------------------------------
 * Initialize a structure for the synchronization of single
 * elements
 *
 * returns
 *   a pointer to a new structure used for synchronizing single elements
 *----------------------------------------------------------------------------*/

static cs_join_sync_t *
_create_join_sync(void)
{
  cs_join_sync_t  *sync = NULL;

  BFT_MALLOC(sync, 1, cs_join_sync_t);

  sync->n_elts = 0;
  sync->n_ranks = 0;
  sync->ranks = NULL;
  sync->index = NULL;
  sync->array = NULL;

  return sync;
}

/*----------------------------------------------------------------------------
 * Destroy a structure for the synchronization of single elements.
 *
 * parameters:
 *   sync   <->  pointer to a structure used for synchronizing single elements
 *----------------------------------------------------------------------------*/

static void
_destroy_join_sync(cs_join_sync_t   **sync)
{
  cs_join_sync_t  *_sync = *sync;

  if (_sync->n_elts > 0)
    BFT_FREE(_sync->array);
  if (_sync->n_ranks > 0)
    BFT_FREE(_sync->ranks);
  BFT_FREE(_sync->index);

  BFT_FREE(_sync);

  *sync = _sync;
}

#if defined(HAVE_MPI)
/*----------------------------------------------------------------------------
 * Define a structure used for synchronizing "single" vertices.
 * Use a fvm_interface_t structure to help the build.
 *
 * parameters:
 *   interfaces     --> pointer to a fvm_interface_set_t structure
 *   var_size       --> number of elements in var buffer
 *   count          <-> counter buffer (0: not selected, 1 otherwise)
 *   related_ranks  <-> rank associated to each single vertex (size: var_size)
 *   single         <-> data about the distribution of single vertices
 *----------------------------------------------------------------------------*/

static void
_add_single_vertices(fvm_interface_set_t    *interfaces,
                     fvm_lnum_t              var_size,
                     fvm_lnum_t             *count,
                     fvm_lnum_t             *related_ranks,
                     cs_join_sync_t         *single)
{
  int  request_count, distant_rank, n_interfaces, total_size;
  int  id, ii, last_found_rank;

  int  count_size = 0;
  fvm_lnum_t  n_entities = 0;
  int  *buf = NULL, *send_buf = NULL, *recv_buf = NULL;

  MPI_Request  *request = NULL;
  MPI_Status  *status  = NULL;
  MPI_Comm  mpi_comm = cs_glob_mpi_comm;

  const int  local_rank = CS_MAX(cs_glob_rank_id, 0);
  const fvm_lnum_t  *local_num = NULL;
  const fvm_interface_t  *interface = NULL;

  assert(count != NULL);

  /* Initialize and allocate */

  n_interfaces = fvm_interface_set_size(interfaces);

  for (id = 0; id < n_interfaces; id++) {
    count_size
      += fvm_interface_size(fvm_interface_set_get(interfaces, id));
  }

  total_size = count_size;

  BFT_MALLOC(buf, total_size * 2, int);

  BFT_MALLOC(request, n_interfaces * 2, MPI_Request);
  BFT_MALLOC(status,  n_interfaces * 2, MPI_Status);

  /* Send and Receive data from distant ranks with
     non-blocking communications */

  request_count = 0;
  count_size  = 0;

  /* Receive */

  for (id = 0; id < n_interfaces; id++) {

    interface = fvm_interface_set_get(interfaces, id);
    distant_rank = fvm_interface_rank(interface);
    n_entities = fvm_interface_size(interface);

    recv_buf = buf + count_size;

    MPI_Irecv(recv_buf,
              n_entities,
              MPI_INT,
              distant_rank,
              distant_rank,
              mpi_comm,
              &(request[request_count++]));

    count_size += n_entities;

  }

  assert(count_size == total_size);

  /* Send */

  for (id = 0; id < n_interfaces; id++) {

    /* Preparation of data to send */

    interface = fvm_interface_set_get(interfaces, id);
    distant_rank   = fvm_interface_rank(interface);
    n_entities = fvm_interface_size(interface);
    local_num = fvm_interface_get_local_num(interface);

    send_buf = buf + count_size;

    for (ii = 0; ii < n_entities; ii++)
      send_buf[ii] = count[local_num[ii]-1];

    MPI_Isend(send_buf,
              n_entities,
              MPI_INT,
              distant_rank,
              local_rank,
              mpi_comm,
              &(request[request_count++]));

    count_size += n_entities;

  }

  assert(count_size == 2*total_size);

  /* Sync after each rank had received all the messages */

  MPI_Waitall(request_count, request, status);

  BFT_FREE(request);
  BFT_FREE(status);

  /* Now we had each part to count */

  count_size = 0;
  last_found_rank = -1;

  for (id = 0 ; id < n_interfaces ; id++) {

    /* Scan data */

    interface = fvm_interface_set_get(interfaces, id);
    distant_rank = fvm_interface_rank(interface);
    n_entities = fvm_interface_size(interface);
    local_num = fvm_interface_get_local_num(interface);

    recv_buf = buf + count_size;

    for (ii = 0; ii < n_entities; ii++) {

      int vtx_id = local_num[ii] - 1;

      assert(vtx_id < var_size);

      if (count[vtx_id] == 0 && recv_buf[ii] > 0) {

        if (last_found_rank != distant_rank) {
          last_found_rank = distant_rank;
          single->n_ranks++;
        }
        single->n_elts++;

      }

    }

    count_size += n_entities;

  }

  if (single->n_elts > 0) {

    int  rank_shift = 0, vtx_shift = 0;

    BFT_MALLOC(single->ranks, single->n_ranks, int);
    BFT_MALLOC(single->index, single->n_ranks + 1, int);
    BFT_MALLOC(single->array, single->n_elts, int);

    count_size = 0;
    last_found_rank = -1;
    single->index[0] = 0;

    for (id = 0 ; id < n_interfaces ; id++) {

      /* Scan data */

      interface = fvm_interface_set_get(interfaces, id);
      distant_rank = fvm_interface_rank(interface);
      n_entities = fvm_interface_size(interface);
      local_num = fvm_interface_get_local_num(interface);

      recv_buf = buf + count_size;

      for (ii = 0; ii < n_entities; ii++) {

        int  vtx_id = local_num[ii] - 1;

        if (count[vtx_id] == 0 && recv_buf[ii] > 0) {

          if (last_found_rank != distant_rank) {
            last_found_rank = distant_rank;
            single->ranks[rank_shift++] = distant_rank;
          }

          single->array[vtx_shift++] = local_num[ii];
          single->index[rank_shift] = vtx_shift;

          related_ranks[vtx_id] = distant_rank;
          count[vtx_id] = recv_buf[ii];

        }

      }

      count_size += n_entities;

    } /* End of loop on interfaces */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    bft_printf("\n  Single vertices for the joining operations:\n");
    for (i = 0, shift = 0; i < single->n_ranks; i++) {
      for (j = single->index[i]; j < single->index[i+1]; j++) {
        bft_printf(" %9d | %6d | %9d\n",
                   shift, single->ranks[i], single->array[j]);
        shift++;
      }
    }
    bft_printf("\n");
    bft_printf_flush();
#endif

  } /* End if single->n_vertices > 0 */

  BFT_FREE(buf);

}

/*-----------------------------------------------------------------------------
 * Define a structure used for synchronizing "single" vertices.
 * Use a fvm_interface_t structure to help the build.
 *
 * parameters:
 *  interfaces     --> pointer to a fvm_interface_set_t structure
 *  var_size       --> number of elements in var buffer
 *  related_ranks  <-> rank buffer for synchronization
 *  coupled        <-> pointer to a structure to build on coupled vertices
 *----------------------------------------------------------------------------*/

static void
_add_coupled_vertices(fvm_interface_set_t    *interfaces,
                      fvm_lnum_t              var_size,
                      int                    *related_ranks,
                      cs_join_sync_t         *coupled)
{
  int  id, ii;
  int  request_count, distant_rank, n_interfaces, total_size, last_found_rank;

  int  count_size = 0;
  fvm_lnum_t  n_entities = 0;
  int  *buf = NULL, *send_buf = NULL, *recv_buf = NULL;

  MPI_Request  *request = NULL;
  MPI_Status  *status  = NULL;
  MPI_Comm  mpi_comm = cs_glob_mpi_comm;

  const int  local_rank = CS_MAX(cs_glob_rank_id, 0);
  const fvm_lnum_t  *local_num = NULL;
  const fvm_interface_t  *interface = NULL;

  assert(related_ranks != NULL);

  /* Initialize and allocate */

  n_interfaces = fvm_interface_set_size(interfaces);

  for (id = 0; id < n_interfaces; id++) {
    count_size
      += fvm_interface_size(fvm_interface_set_get(interfaces, id));
  }

  total_size = count_size;

  BFT_MALLOC(buf, total_size * 2, int);

  BFT_MALLOC(request, n_interfaces * 2, MPI_Request);
  BFT_MALLOC(status,  n_interfaces * 2, MPI_Status);

  /* Send and Receive data from distant ranks with
     non-blocking communications */

  request_count = 0;
  count_size  = 0;

  /* Receive */

  for (id = 0; id < n_interfaces; id++) {

    interface = fvm_interface_set_get(interfaces, id);
    distant_rank = fvm_interface_rank(interface);
    n_entities = fvm_interface_size(interface);

    recv_buf = buf + count_size;

    MPI_Irecv(recv_buf,
              n_entities,
              MPI_INT,
              distant_rank,
              distant_rank,
              mpi_comm,
              &(request[request_count++]));

    count_size += n_entities;

  }

  assert(count_size == total_size);

  /* Send */

  for (id = 0; id < n_interfaces; id++) {

    /* Preparation of data to send */

    interface = fvm_interface_set_get(interfaces, id);
    distant_rank   = fvm_interface_rank(interface);
    n_entities = fvm_interface_size(interface);
    local_num = fvm_interface_get_local_num(interface);

    send_buf = buf + count_size;

    for (ii = 0; ii < n_entities; ii++)
      send_buf[ii] = related_ranks[local_num[ii]-1];

    MPI_Isend(send_buf,
              n_entities,
              MPI_INT,
              distant_rank,
              local_rank,
              mpi_comm,
              &(request[request_count++]));

    count_size += n_entities;

  }

  assert(count_size == 2*total_size);

  /* Sync after each rank had received all the messages */

  MPI_Waitall(request_count, request, status);

  BFT_FREE(request);
  BFT_FREE(status);

  /* Now we scan each part to build coupled */

  count_size = 0;
  coupled->n_elts = 0;
  coupled->n_ranks = 0;
  last_found_rank = -1;

  for (id = 0; id < n_interfaces; id++) {

    /* Scan data */

    interface = fvm_interface_set_get(interfaces, id);
    distant_rank   = fvm_interface_rank(interface);
    n_entities = fvm_interface_size(interface);

    recv_buf = buf + count_size;

    for (ii = 0 ; ii < n_entities ; ii++) {

      if (recv_buf[ii] == local_rank) {
        coupled->n_elts++;
        if (last_found_rank != distant_rank) {
          last_found_rank = distant_rank;
          coupled->n_ranks++;
        }
      }

    }
    count_size += n_entities;

  }

  if (coupled->n_elts > 0) {

    int  rank_shift = 0, vtx_shift = 0;

    BFT_MALLOC(coupled->array, coupled->n_elts, int);
    BFT_MALLOC(coupled->index, coupled->n_ranks + 1, int);
    BFT_MALLOC(coupled->ranks, coupled->n_ranks, int);

    coupled->index[0] = 0;

    count_size = 0;
    last_found_rank = -1;

    for (id = 0; id < n_interfaces; id++) {

      /* Retrieve data */

      interface = fvm_interface_set_get(interfaces, id);
      distant_rank   = fvm_interface_rank(interface);
      n_entities = fvm_interface_size(interface);
      local_num = fvm_interface_get_local_num(interface);

      recv_buf = buf + count_size;

      for (ii = 0; ii < n_entities; ii++) {

        if (recv_buf[ii] == local_rank) {

          if (last_found_rank != distant_rank) {
            last_found_rank = distant_rank;
            coupled->ranks[rank_shift++] = distant_rank;
          }

          assert(local_num[ii] <= var_size);
          coupled->array[vtx_shift++] = local_num[ii];
          coupled->index[rank_shift] = vtx_shift;
        }

      }
      count_size += n_entities;

    }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    bft_printf("\n  Coupled vertices for the joining operations:\n");
    for (i = 0, shift = 0; i < coupled->n_ranks; i++) {
      for (j = coupled->index[i]; j < coupled->index[i+1]; j++) {
        bft_printf(" %9d | %6d | %9d\n",
                   shift, coupled->ranks[i], coupled->array[j]);
        shift++;
      }
    }
    bft_printf("\n");
    bft_printf_flush();
#endif

  } /* End if coupled->n_elts > 0 */

  BFT_FREE(buf);

}

/*----------------------------------------------------------------------------
 * Get the related edge id from a couple of vertex ids.
 *
 * parameters:
 *  v1_id        -->  first vertex id
 *  v2_id        -->  second vertex id
 *  v2v_idx      -->  vertex -> vertex connect. index
 *  v2v_lst      -->  vertex -> vertex connect. list
 *
 * return:
 *  related edge_id in cs_join_edges_t structure
 *---------------------------------------------------------------------------*/

inline static cs_int_t
_get_edge_id(cs_int_t      v1_id,
             cs_int_t      v2_id,
             cs_int_t      v2v_idx[],
             cs_int_t      v2v_lst[])
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
 *  vertex_tag    -->  tag to know if a vertices is in selection
 *  v1_id         -->  first vertex id
 *  v2_id         -->  second vertex id
 *  ref_v2v_idx   -->  vertex -> vertex connect. index
 *  ref_v2v_lst   -->  vertex -> vertex connect. list
 *  p_tmp_size    <->  pointer to the current number of single edges
 *  p_max_size    <->  pointer to the max allocated size of new_s_vertices
 *  p_tmp_size    <->  pointer to the single edges definition
 *---------------------------------------------------------------------------*/

static void
_add_s_edge(cs_int_t      vertex_tag[],
            cs_int_t      v1_id,
            cs_int_t      v2_id,
            cs_int_t      sel_v2v_idx[],
            cs_int_t      sel_v2v_lst[],
            int          *p_tmp_size,
            int          *p_max_size,
            int          *p_tmp_edges[])
{
  int  i, a, b, edge_id;

  _Bool  is_selected = false, is_found = false;

  int  tmp_size = *p_tmp_size;
  int  max_size = *p_max_size;
  int  *tmp_edges = *p_tmp_edges;

  if (vertex_tag[v1_id] > 0)
    if (vertex_tag[v2_id] > 0)
      is_selected = true;

  if (is_selected) { /* Check if this edge is in sel_edges */

    edge_id = _get_edge_id(v1_id, v2_id, sel_v2v_idx, sel_v2v_lst);

    if (edge_id == -1) { /* Edge not found among the selected edges */

      if (v1_id < v2_id)
        a = v1_id + 1, b = v2_id + 1;
      else
        a = v2_id + 1, b = v1_id + 1;

      /* n_edges is assumed to be small */
      for (i = 0; i < tmp_size; i++) {
        if (tmp_edges[2*i] == a) {
          if (tmp_edges[2*i+1] == b) {
            is_found = true;
            break;
          }
        }
      }

      if (!is_found) {

        tmp_edges[2*tmp_size] = a;
        tmp_edges[2*tmp_size+1] = b;
        tmp_size += 1;

        if (tmp_size >= max_size) {
          max_size *= 2;
          BFT_REALLOC(tmp_edges, 2*max_size, int);
        }

      }

    } /* this edge is not among the selection */

  } /* this edge should be among the selected edges */

  /* Return pointers */

  *p_max_size = max_size;
  *p_tmp_size = tmp_size;
  *p_tmp_edges = tmp_edges;

}

/*----------------------------------------------------------------------------
 * Get the full selection of single edges. Done only if the run is parallel.
 *
 * parameters:
 *  ifs          -->  pointer to the interface set on vertices
 *  n_vertices   -->  number of vertices in the parent mesh
 *  vertex_tag   -->  tag to know if a vertices is in selection
 *  selection    -->  pointer to a fvm_join_selection_t structure
 *  b_f2v_idx    -->  border "face -> vertex" connect. index
 *  b_f2v_lst    -->  border "face -> vertex" connect. list
 *  i_f2v_idx    -->  interior "face -> vertex" connect. index
 *  i_f2v_lst    -->  interior "face -> vertex" connect. list
 *  i_face_cells -->  interior face -> cells connect.
 *  s_edges      <->  pointer to the single edges structure to define
 *---------------------------------------------------------------------------*/

static void
_add_single_edges(fvm_interface_set_t   *ifs,
                  cs_int_t               n_vertices,
                  cs_int_t               vertex_tag[],
                  cs_join_select_t      *selection,
                  cs_int_t               b_f2v_idx[],
                  cs_int_t               b_f2v_lst[],
                  cs_int_t               i_f2v_idx[],
                  cs_int_t               i_f2v_lst[],
                  cs_int_t               i_face_cells[],
                  cs_join_sync_t        *s_edges)
{
  cs_int_t  i, j, fid, save, shift, s, e, n_sel_edges;

  int  tmp_size = 0, max_size = 10;
  int  *tmp_edges = NULL;
  cs_int_t  *count = NULL, *sel_v2v_idx = NULL, *sel_v2v_lst = NULL;

  /* Build a vertex -> vertex connectivity for the selected border faces  */

  BFT_MALLOC(sel_v2v_idx, n_vertices + 1, cs_int_t);

  for (i = 0; i < n_vertices + 1; i++)
    sel_v2v_idx[i] = 0;

  cs_join_build_edges_idx(selection->n_faces,
                          selection->faces,
                          b_f2v_idx,
                          b_f2v_lst,
                          sel_v2v_idx);

  BFT_MALLOC(count, n_vertices, cs_int_t);

  for (i = 0; i < n_vertices; i++) {
    sel_v2v_idx[i+1] += sel_v2v_idx[i];
    count[i] = 0;
  }

  BFT_MALLOC(sel_v2v_lst, sel_v2v_idx[n_vertices], cs_int_t);

  cs_join_build_edges_lst(selection->n_faces,
                          selection->faces,
                          b_f2v_idx,
                          b_f2v_lst,
                          count,
                          sel_v2v_idx,
                          sel_v2v_lst);

  BFT_FREE(count);

  /* Order sub-lists and then clean repeated elements */

  for (i = 0; i < n_vertices; i++)
    cs_sort_shell(sel_v2v_idx[i], sel_v2v_idx[i+1], sel_v2v_lst);

  save = sel_v2v_idx[0];
  shift = 0;

  for (i = 0; i < n_vertices; i++) {

    s = save;
    e = sel_v2v_idx[i+1];

    if (e - s > 0) {
      sel_v2v_lst[shift++] = sel_v2v_lst[s];
      for (j = s + 1; j < e; j++)
        if (sel_v2v_lst[j-1] != sel_v2v_lst[j])
          sel_v2v_lst[shift++] = sel_v2v_lst[j];
    }

    save = e;
    sel_v2v_idx[i+1] = shift;

  }

  n_sel_edges = sel_v2v_idx[n_vertices];
  BFT_REALLOC(sel_v2v_lst, n_sel_edges, cs_int_t);

  /* Scan adjacent faces to find new "single" edges */

  BFT_MALLOC(tmp_edges, 2*max_size, int);

  for (i = 0; i < selection->n_b_adj_faces; i++) {

    fid = selection->b_adj_faces[i] - 1;
    s = b_f2v_idx[fid] - 1;
    e = b_f2v_idx[fid+1] - 1;

    for (j = s; j < e - 1; j++)
      _add_s_edge(vertex_tag,
                  b_f2v_lst[j]-1,
                  b_f2v_lst[j+1]-1,
                  sel_v2v_idx,
                  sel_v2v_lst,
                  &tmp_size,
                  &max_size,
                  &tmp_edges);

    _add_s_edge(vertex_tag,
                b_f2v_lst[e-1]-1,
                b_f2v_lst[s]-1,
                sel_v2v_idx,
                sel_v2v_lst,
                &tmp_size,
                &max_size,
                &tmp_edges);

  }

  for (i = 0; i < selection->n_i_adj_faces; i++) {

    fid = selection->i_adj_faces[i] - 1;

    if (i_face_cells[2*fid] == 0 || i_face_cells[2*fid+1] == 0) {

      s = i_f2v_idx[fid] - 1;
      e = i_f2v_idx[fid+1] - 1;

      for (j = s; j < e - 1; j++)
        _add_s_edge(vertex_tag,
                    i_f2v_lst[j]-1,
                    i_f2v_lst[j+1]-1,
                    sel_v2v_idx,
                    sel_v2v_lst,
                    &tmp_size,
                    &max_size,
                    &tmp_edges);

      _add_s_edge(vertex_tag,
                  i_f2v_lst[e-1]-1,
                  i_f2v_lst[s]-1,
                  sel_v2v_idx,
                  sel_v2v_lst,
                  &tmp_size,
                  &max_size,
                  &tmp_edges);

    } /* Face on a parallel frontier */

  }

  /* Find the related ranks for synchronization */

  if (tmp_size > 0) {

    int   n_entities, distant_rank;

    int  *edge_tag = NULL;

    const int  n_interfaces = fvm_interface_set_size(ifs);
    const fvm_lnum_t  *local_num = NULL;
    const fvm_interface_t  *interface = NULL;

    s_edges->n_elts = tmp_size;
    BFT_REALLOC(tmp_edges, 2*s_edges->n_elts, int);

    BFT_MALLOC(edge_tag, s_edges->n_elts, int);
    BFT_MALLOC(s_edges->array, 2*s_edges->n_elts, int);
    BFT_MALLOC(s_edges->index, n_interfaces + 1, int);
    BFT_MALLOC(s_edges->ranks, n_interfaces, int);

    for (i = 0; i < s_edges->n_elts; i++)
      edge_tag[i] = 0; /* Not matched */

    for (i = 0; i < n_interfaces; i++) {
      s_edges->index[i] = 0;
      s_edges->ranks[i] = -1;
    }
    s_edges->index[n_interfaces] = 0;

    for (i = 0; i < n_interfaces; i++) {

      interface = fvm_interface_set_get(ifs, i);
      distant_rank = fvm_interface_rank(interface);
      n_entities = fvm_interface_size(interface);
      local_num = fvm_interface_get_local_num(interface);

      for (j = 0; j < s_edges->n_elts; j++) {

        if (edge_tag[j] == 0) {
          if (cs_search_binary(n_entities, tmp_edges[2*j], local_num) > -1) {
            if (cs_search_binary(n_entities, tmp_edges[2*j+1], local_num) > -1) {

              s_edges->array[2*s_edges->index[i+1]] = tmp_edges[2*j];
              s_edges->array[2*s_edges->index[i+1]+1] = tmp_edges[2*j+1];
              s_edges->index[i+1] += 1;
              edge_tag[j] = 1;

              if (s_edges->n_ranks == 0 ||
                  (   s_edges->n_ranks > 0
                   && s_edges->ranks[s_edges->n_ranks-1] != distant_rank) ) {
                s_edges->ranks[s_edges->n_ranks] = distant_rank;
                s_edges->n_ranks++;
              }

            }
          }
        } /* Not matched yet */

      } /* Loop on single edges */

    } /* Loop on interfaces */

    /* sanity check */

    for (j = 0; j < s_edges->n_elts; j++)
      if (edge_tag[j] == 0)
        bft_error(__FILE__, __LINE__, 0,
                  _(" Can't find the distant rank in the interface set"
                    " for the current edge [%d, %d] (local num.)\n"),
                  s_edges->array[2*j], s_edges->array[2*j+1]);

    /* Memory management */

    BFT_REALLOC(s_edges->ranks, s_edges->n_ranks, int);
    BFT_REALLOC(s_edges->index, s_edges->n_ranks + 1, int);
    BFT_FREE(edge_tag);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    bft_printf("\n  Single edges for the joining operations:\n");
    for (i = 0, shift = 0; i < s_edges->n_ranks; i++) {
      for (j = s_edges->index[i]; j < s_edges->index[i+1]; j++) {
        bft_printf(" %9d | %6d | (%9d, %9d)\n",
                   shift, s_edges->ranks[i],
                   s_edges->array[2*j], s_edges->array[2*j+1]);
        shift++;
      }
    }
    bft_printf("\n");
    bft_printf_flush();
#endif

  }

  /* Free memory */

  BFT_FREE(tmp_edges);
  BFT_FREE(sel_v2v_idx);
  BFT_FREE(sel_v2v_lst);

}

/*----------------------------------------------------------------------------
 * Define a structure for the coupled edges. Only done if single edges have
 * been detected and only in parallel.
 *
 * parameters:
 *  ifs          -->  pointer to the interface set on vertices
 *  s_edges      -->  single edges structure used to build coupled_edges
 *  c_edges      <->  pointer to the coupled edges structure to define
 *---------------------------------------------------------------------------*/

static void
_add_coupled_edges(fvm_interface_set_t   *ifs,
                   cs_join_sync_t        *s_edges,
                   cs_join_sync_t        *c_edges)
{
  int  i, j, id, n_entities;
  int  request_count, rank_shift, distant_rank;

  int  *buf = NULL, *recv_buf = NULL, *send_buf = NULL;
  MPI_Request  *request = NULL;
  MPI_Status  *status  = NULL;
  MPI_Comm  mpi_comm = cs_glob_mpi_comm;

  const int  local_rank = CS_MAX(cs_glob_rank_id, 0);
  const fvm_lnum_t  *local_num = NULL, *distant_num = NULL;
  const int  n_interfaces = fvm_interface_set_size(ifs);
  const fvm_interface_t  *interface = NULL;

  /* Exchange number of single edges */

  BFT_MALLOC(request, n_interfaces * 2, MPI_Request);
  BFT_MALLOC(status,  n_interfaces * 2, MPI_Status);
  BFT_MALLOC(buf, 2*n_interfaces, int);

  for (i = 0; i < 2*n_interfaces; i++)
    buf[i] = 0;

  request_count = 0;

  for (id = 0; id < n_interfaces; id++) {

    interface = fvm_interface_set_get(ifs, id);
    distant_rank = fvm_interface_rank(interface);

    MPI_Irecv(&(buf[id]),
              1,
              MPI_INT,
              distant_rank,
              distant_rank,
              mpi_comm,
              &(request[request_count++]));

  }

  /* Send */

  rank_shift = 0;

  for (id = 0; id < n_interfaces; id++) {

    /* Preparation of data to send */

    interface = fvm_interface_set_get(ifs, id);
    distant_rank = fvm_interface_rank(interface);

    if (rank_shift < s_edges->n_ranks) {
      if (s_edges->ranks[rank_shift] == distant_rank) {
        buf[n_interfaces + id] =
          s_edges->index[rank_shift+1] - s_edges->index[rank_shift];
        rank_shift++;
      }
    }

    MPI_Isend(&(buf[n_interfaces + id]),
              1,
              MPI_INT,
              distant_rank,
              local_rank,
              mpi_comm,
              &(request[request_count++]));

  } /* End of loop on interfaces */

  /* Sync after each rank had received all the messages */

  MPI_Waitall(request_count, request, status);

  /* Define c_edges */

  for (i = 0; i < n_interfaces; i++) {
    if (buf[i] > 0) {
      c_edges->n_elts += buf[i];
      c_edges->n_ranks += 1;
    }
  }

  BFT_MALLOC(c_edges->ranks, c_edges->n_ranks, int);
  BFT_MALLOC(c_edges->index, c_edges->n_ranks + 1, int);
  BFT_MALLOC(c_edges->array, 2*c_edges->n_elts, int);

  for (i = 0; i < c_edges->n_ranks + 1; i++)
    c_edges->index[i] = 0;

  rank_shift = 0;

  for (i = 0; i < n_interfaces; i++) {
    if (buf[i] > 0) {

      interface = fvm_interface_set_get(ifs, i);
      distant_rank = fvm_interface_rank(interface);
      c_edges->ranks[rank_shift++] = distant_rank;
      c_edges->index[rank_shift] = buf[i];

    }
  }

  assert(rank_shift == c_edges->n_ranks);

  for (i = 0; i < c_edges->n_ranks; i++)
    c_edges->index[i+1] += c_edges->index[i];

  assert(c_edges->index[c_edges->n_ranks] == c_edges->n_elts);

  BFT_FREE(buf);

  /* Exchange couple of vertices defining coupled edges */

  request_count = 0;
  rank_shift = 0;

  for (id = 0; id < n_interfaces; id++) {

    interface = fvm_interface_set_get(ifs, id);
    distant_rank = fvm_interface_rank(interface);

    if (rank_shift < c_edges->n_ranks) {
      if (c_edges->ranks[rank_shift] == distant_rank) {

        n_entities = c_edges->index[rank_shift+1] - c_edges->index[rank_shift];
        n_entities *= 2; /* couple of vertices */
        recv_buf = c_edges->array + 2*c_edges->index[rank_shift];
        rank_shift++;

        MPI_Irecv(recv_buf, /* receive distant num */
                  n_entities,
                  MPI_INT,
                  distant_rank,
                  distant_rank,
                  mpi_comm,
                  &(request[request_count++]));

      }
    }

  } /* End of loop on interfaces */

  /* Send */

  rank_shift = 0;

  for (id = 0; id < n_interfaces; id++) {

    interface = fvm_interface_set_get(ifs, id);
    distant_rank = fvm_interface_rank(interface);

    if (rank_shift < s_edges->n_ranks) {
      if (s_edges->ranks[rank_shift] == distant_rank) {

        n_entities = s_edges->index[rank_shift+1] - s_edges->index[rank_shift];
        n_entities *= 2; /* couple of vertices */
        send_buf = s_edges->array + 2*s_edges->index[rank_shift];
        rank_shift++;

        MPI_Isend(send_buf,
                  n_entities,
                  MPI_INT,
                  distant_rank,
                  local_rank,
                  mpi_comm,
                  &(request[request_count++]));

      }
    }

  } /* End of loop on interfaces */

  /* Sync after each rank had received all the messages */

  MPI_Waitall(request_count, request, status);

  BFT_FREE(request);
  BFT_FREE(status);

  rank_shift = 0;

  /* Switch received couple vertices from distant num. to local num. */

  for (i = 0; i < n_interfaces; i++) {

    interface = fvm_interface_set_get(ifs, i);
    distant_rank = fvm_interface_rank(interface);

    if (rank_shift < c_edges->n_ranks) {
      if (c_edges->ranks[rank_shift] == distant_rank) {

        n_entities = fvm_interface_size(interface);
        local_num = fvm_interface_get_local_num(interface);
        distant_num = fvm_interface_get_distant_num(interface);

        for (j = c_edges->index[rank_shift];
             j < c_edges->index[rank_shift+1]; j++) {

          id = cs_search_binary(n_entities, c_edges->array[2*j], distant_num);
          assert(id != -1);
          c_edges->array[2*j] = local_num[id];

          id = cs_search_binary(n_entities, c_edges->array[2*j+1], distant_num);
          assert(id != -1);
          c_edges->array[2*j+1] = local_num[id];

        } /* Loop on couple edges */

        rank_shift++;

      } /* This rank is in the list of related ranks */
    }

  } /* Loop on interfaces */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (c_edges->n_elts > 0) {
    int  shift;
    bft_printf("\n  Coupled edges for the joining operations:\n");
    for (i = 0, shift = 0; i < c_edges->n_ranks; i++) {
      for (j = c_edges->index[i]; j < c_edges->index[i+1]; j++) {
        bft_printf(" %9d | %6d | (%9d, %9d)\n",
                   shift, c_edges->ranks[i],
                   c_edges->array[2*j], c_edges->array[2*j+1]);
        shift++;
      }
    }
    bft_printf("\n");
    bft_printf_flush();
  }
#endif

}

/*----------------------------------------------------------------------------
 * Get the full selection of vertices to extract. Only in parallel. Somme
 * vertices may have been selected on another rank and not on the local rank
 * but you have to take them into account to have a good update of the mesh
 *
 * parameters:
 *  b_f2v_idx       -->  border "face -> vertex" connect. index
 *  b_f2v_lst       -->  border "face -> vertex" connect. list
 *  i_f2v_idx       -->  interior "face -> vertex" connect. index
 *  i_f2v_lst       -->  interior "face -> vertex" connect. list
 *  n_vertices      -->  number of vertices in the parent mesh
 *  v_gnum          -->  global vertex numbering (NULL if n_ranks = 1)
 *  i_face_cells    -->  interior face -> cells connect.
 *  join_select     <->  pointer to a fvm_join_selection_t structure
 *---------------------------------------------------------------------------*/

static void
_get_single_elements(cs_int_t            b_f2v_idx[],
                     cs_int_t            b_f2v_lst[],
                     cs_int_t            i_f2v_idx[],
                     cs_int_t            i_f2v_lst[],
                     cs_int_t            n_vertices,
                     fvm_gnum_t          v_gnum[],
                     cs_int_t            i_face_cells[],
                     cs_join_select_t   *selection)
{
  cs_int_t  i;
  fvm_gnum_t  n_g_elts;

  cs_int_t  *vtx_tag = NULL, *related_ranks = NULL;
  fvm_interface_set_t  *ifs = NULL;

  MPI_Comm  mpi_comm = cs_glob_mpi_comm;

  assert(v_gnum != NULL);

  /* Build an interface on vertices to detect single vertices */

  ifs = fvm_interface_set_create(n_vertices,
                                 NULL,
                                 v_gnum,
                                 NULL,
                                 0,
                                 NULL,
                                 NULL,
                                 NULL);

  assert(ifs != NULL);

  /* Define a counter on vertices. 1 if selected, 0 otherwise */

  BFT_MALLOC(vtx_tag, n_vertices, cs_int_t);
  BFT_MALLOC(related_ranks, n_vertices, cs_int_t);

  for (i = 0; i < n_vertices; i++) {
    vtx_tag[i] = 0;
    related_ranks[i] = -1;
  }

  for (i = 0; i < selection->n_vertices; i++)
    vtx_tag[selection->vertices[i] - 1] = 1;

  /* Has to be done before the search of single edges because
     we need a synchronized vtx_tag */

  _add_single_vertices(ifs,
                       n_vertices,
                       vtx_tag,
                       related_ranks,
                       selection->s_vertices);

  MPI_Allreduce(&(selection->s_vertices->n_elts), &n_g_elts, 1, FVM_MPI_GNUM,
                MPI_SUM, mpi_comm);

  if (n_g_elts > 0) {

    bft_printf("\n  Global number of single vertices found: %6u\n", n_g_elts);
    bft_printf_flush();
    selection->do_single_sync = true;

    _add_coupled_vertices(ifs,
                          n_vertices,
                          related_ranks,
                          selection->c_vertices);

  } /* End if n_g_s_vertices > 0 */

  _add_single_edges(ifs,
                    n_vertices,
                    vtx_tag,
                    selection,
                    b_f2v_idx,
                    b_f2v_lst,
                    i_f2v_idx,
                    i_f2v_lst,
                    i_face_cells,
                    selection->s_edges);

  MPI_Allreduce(&(selection->s_edges->n_elts), &n_g_elts, 1, FVM_MPI_GNUM,
                MPI_SUM, mpi_comm);

  if (n_g_elts > 0) {

    bft_printf("  Global number of single edges found: %6d\n", n_g_elts);
    bft_printf_flush();
    selection->do_single_sync = true;

    _add_coupled_edges(ifs,
                       selection->s_edges,
                       selection->c_edges);

  } /* End if n_g_s_elts > 0 */

  /* Free memory */

  BFT_FREE(vtx_tag);
  BFT_FREE(related_ranks);

  ifs = fvm_interface_set_destroy(ifs);

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

        cs_int_t  face_num = c2f_lst[j];

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
     computation of tolerance (Not a user-defined parameter) */

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
  cs_int_t  i, fid, cid;

  cs_join_select_t  *selection = NULL;
  fvm_lnum_t  *order = NULL, *ordered_faces = NULL;
  cs_mesh_t  *mesh = cs_glob_mesh;

  const int  n_ranks = cs_glob_n_ranks;

  assert(mesh != NULL);

  /* Initialize cs_join_select_t struct. */

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
  selection->s_vertices = _create_join_sync();
  selection->c_vertices = _create_join_sync();
  selection->s_edges = _create_join_sync();
  selection->c_edges = _create_join_sync();

  /* Extract selected border faces */

  BFT_MALLOC(selection->faces, mesh->n_b_faces, fvm_lnum_t);

  cs_selector_get_b_face_list(selection_criteria,
                              &(selection->n_faces),
                              selection->faces);

  BFT_MALLOC(order, selection->n_faces, fvm_lnum_t);
  BFT_MALLOC(ordered_faces, selection->n_faces, fvm_lnum_t);

  fvm_order_local_allocated(selection->faces, NULL, order, selection->n_faces);

  for (i = 0; i < selection->n_faces; i++)
    ordered_faces[i] = selection->faces[order[i]];

  BFT_FREE(order);
  BFT_FREE(selection->faces);
  selection->faces = ordered_faces;

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

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf(_("\n  Selected faces for the joining operation:\n"));
  for (i = 0; i < selection->n_faces; i++)
    bft_printf(" %9d | %9d | %10u | %10u\n",
               i, selection->faces[i], selection->compact_face_gnum[i],
               selection->cell_gnum[i]);
  bft_printf("\n");
#endif

  /* Extract selected vertices from the selected border faces */

  _extract_vertices(selection->n_faces,
                    selection->faces,
                    mesh->b_face_vtx_idx,
                    mesh->b_face_vtx_lst,
                    mesh->n_vertices,
                    &(selection->n_vertices),
                    &(selection->vertices));

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf(_("\n  Select vertices for the joining operation:\n"));
  for (i = 0; i < selection->n_vertices; i++)
    bft_printf(" %9d | %9d\n", i, selection->vertices[i]);
  bft_printf("\n");
#endif

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

  _clean_selection(&(selection->n_b_adj_faces),
                   &(selection->b_adj_faces),
                   selection->n_faces,
                   selection->faces);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf(_("\n  Contiguous border faces for the joining operation:\n"));
  for (i = 0; i < selection->n_b_adj_faces; i++)
    bft_printf(" %9d | %9d\n", i, selection->b_adj_faces[i]);
  bft_printf("\n");
#endif

  /* Extract list of interior faces contiguous to the selected vertices */

  _extract_contig_faces(mesh->n_vertices,
                        selection->n_vertices,
                        selection->vertices,
                        mesh->n_i_faces,
                        mesh->i_face_vtx_idx,
                        mesh->i_face_vtx_lst,
                        &(selection->n_i_adj_faces),
                        &(selection->i_adj_faces));

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf(_("\n  Contiguous interior faces for the joining operation:\n"));
  for (i = 0; i < selection->n_i_adj_faces; i++)
    bft_printf(" %9d | %9d\n", i, selection->i_adj_faces[i]);
  bft_printf("\n");
#endif

   /* Check if there is no forgotten vertex in the selection.
      Otherwise define structures to enable future synchronization.
      Only possible in parallel. */

#if defined(HAVE_MPI)
  if (n_ranks > 1)
    _get_single_elements(mesh->b_face_vtx_idx,
                         mesh->b_face_vtx_lst,
                         mesh->i_face_vtx_idx,
                         mesh->i_face_vtx_lst,
                         mesh->n_vertices,
                         mesh->global_vtx_num,
                         mesh->i_face_cells,
                         selection);
#endif

  /* Display information according to the level of verbosity */

  if (verbosity > 1) {

    bft_printf("\n  Compact index on ranks for the selected faces:\n");
    for (i = 0; i < n_ranks + 1; i++)
      bft_printf(" %5d | %11u\n", i, selection->compact_rank_index[i]);
    bft_printf("\n");

    if (selection->do_single_sync == true) {
      bft_printf("\n Information about single/coupled elements:\n");
      bft_printf("   Number of single vertices : %6d with %3d related ranks\n",
                 selection->s_vertices->n_elts, selection->s_vertices->n_ranks);
      bft_printf("   Number of coupled vertices: %6d with %3d related ranks\n",
                 selection->c_vertices->n_elts, selection->c_vertices->n_ranks);
      bft_printf("   Number of single edges    : %6d with %3d related ranks\n",
                 selection->s_edges->n_elts, selection->s_edges->n_ranks);
      bft_printf("   Number of coupled edges   : %6d with %3d related ranks\n",
                 selection->c_edges->n_elts, selection->c_edges->n_ranks);
    }

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

    _destroy_join_sync(&(_js->s_vertices));
    _destroy_join_sync(&(_js->c_vertices));
    _destroy_join_sync(&(_js->s_edges));
    _destroy_join_sync(&(_js->c_edges));

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
