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
 * Management of conforming and non-conforming join
 *===========================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *---------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *---------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>
#include <bft_timer.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *---------------------------------------------------------------------------*/

#include <fvm_io_num.h>
#include <fvm_order.h>
#include <fvm_parall.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *---------------------------------------------------------------------------*/

#include "cs_search.h"
#include "cs_sort.h"
#include "cs_join_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *---------------------------------------------------------------------------*/

#include "cs_join_update.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Static global variables
 *===========================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *===========================================================================*/

/* Local structure to update the mesh connectivty after the merge step */

typedef struct {

  cs_int_t   n_vertices;  /* v2v_idx size - 1 */
  cs_int_t   n_edges;     /* v2v_lst size */

  cs_int_t  *v2v_idx;
  cs_int_t  *v2v_lst;

  cs_int_t  *v2v_sub_idx;
  cs_int_t  *v2v_sub_lst;  /* if -1: edge has been deleted */

} edge_builder_t;

/*============================================================================
 * Private function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Compute the cross product of two vectors.
 *
 * parameters:
 *  v1     <--  first vector
 *  v2     <--  second vector
 *  result -->  cross product v1xv2
 *
 * returns:
 *  the resulting cross product (v1 x v2)
 *----------------------------------------------------------------------------*/

inline static void
_cross_product(double   v1[],
               double   v2[],
               double   result[])
{
  result[0] = v1[1] * v2[2] - v2[1] * v1[2];
  result[1] = v2[0] * v1[2] - v1[0] * v2[2];
  result[2] = v1[0] * v2[1] - v2[0] * v1[1];
}

/*----------------------------------------------------------------------------
 * Compute the dot product of two vectors.
 *
 * parameters:
 *  v1     <--  first vector
 *  v2     <--  second vector
 *
 * returns:
 *  the resulting dot product (v1.v2)
 *----------------------------------------------------------------------------*/

inline static double
_dot_product(double   v1[],
             double   v2[])
{
  int  i;
  double  result = 0.0;

  for (i = 0; i < 3; i++)
    result += v1[i] * v2[i];

  return result;
}

/*----------------------------------------------------------------------------
 * Normalize a vector.
 *
 * parameters:
 *  v      <->  vector to normalize
 *----------------------------------------------------------------------------*/

inline static void
_normalize(double   v[])
{
  int  i;
  double  norm, inv_norm;

  norm = sqrt(_dot_product(v, v));
  inv_norm = 1/norm;

  for (i = 0; i < 3; i++)
    v[i] *= inv_norm;
}

/*----------------------------------------------------------------------------
 * Get the related edge id in edge_builder_t structure from a couple of
 * vertex ids.
 *
 * parameters:
 *  v1_id         <--  first vertex id
 *  v2_id         <--  second vertex id
 *  edge_builder  <--  pointer to an edge_builder_t structure
 *
 * return:
 *  related edge_id in edge_builder_t structure
 *---------------------------------------------------------------------------*/

inline static cs_int_t
_get_join_edge_id(cs_int_t           v1_id,
                  cs_int_t           v2_id,
                  edge_builder_t    *edge_builder)
{
  cs_int_t  i, va, vb;
  cs_int_t  edge_id = -1;

  if (v1_id < v2_id)
    va = v1_id, vb = v2_id;
  else
    vb = v1_id, va = v2_id;

  for (i = edge_builder->v2v_idx[va]; i < edge_builder->v2v_idx[va+1]; i++)
    if (edge_builder->v2v_lst[i] == vb + 1)
      break;

  if (i != edge_builder->v2v_idx[va+1])
    edge_id = i;

  return  edge_id;
}

#if defined(HAVE_MPI)
/*----------------------------------------------------------------------------
 * Retrieve the local new global numbering for the initial vertices from
 * the new vertex global numbering defined by block.
 *
 * parameters:
 *  mesh               <--  pointer of pointer to cs_mesh_t structure
 *  p_o2n_vtx_gnum     <->  in : array on blocks on the new global vertex
 *                          out: local array on the new global vertex
 *---------------------------------------------------------------------------*/

static void
_get_local_o2n_vtx_gnum(cs_mesh_t    *mesh,
                        fvm_gnum_t   *p_o2n_vtx_gnum[])
{
  cs_int_t  i, shift, rank;
  fvm_gnum_t  new_gnum;

  cs_int_t  *send_shift = NULL, *recv_shift = NULL;
  cs_int_t  *send_count = NULL, *recv_count = NULL;
  fvm_gnum_t  *send_glist = NULL, *recv_glist = NULL;
  fvm_gnum_t  *block_gnum = *p_o2n_vtx_gnum;
  fvm_gnum_t  *local_gnum = NULL;

  MPI_Comm  mpi_comm = cs_glob_mpi_comm;

  const cs_int_t  n_ranks = cs_glob_n_ranks;
  const cs_int_t  local_rank = cs_glob_rank_id;
  const cs_join_block_info_t  block_info
    = cs_join_get_block_info(mesh->n_g_vertices,
                             n_ranks,
                             local_rank);

  BFT_MALLOC(local_gnum, mesh->n_vertices, fvm_gnum_t);

  /* Request the new vtx gnum related to the initial vtx gnum */

  BFT_MALLOC(send_count, n_ranks, cs_int_t);
  BFT_MALLOC(recv_count, n_ranks, cs_int_t);

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (i = 0; i < mesh->n_vertices; i++) {
    rank = (mesh->global_vtx_num[i] - 1)/block_info.size;
    send_count[rank] += 1;
  }

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);

  BFT_MALLOC(send_shift, n_ranks + 1, cs_int_t);
  BFT_MALLOC(recv_shift, n_ranks + 1, cs_int_t);

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (rank = 0; rank < n_ranks; rank++) {
    send_shift[rank + 1] = send_shift[rank] + send_count[rank];
    recv_shift[rank + 1] = recv_shift[rank] + recv_count[rank];
  }

  /* Build send_list */

  BFT_MALLOC(send_glist, send_shift[n_ranks], fvm_gnum_t);
  BFT_MALLOC(recv_glist, recv_shift[n_ranks], fvm_gnum_t);

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (i = 0; i < mesh->n_vertices; i++) {

    rank = (mesh->global_vtx_num[i] - 1)/block_info.size;
    shift = send_shift[rank] + send_count[rank];
    send_glist[shift] = mesh->global_vtx_num[i];  /* Old global number */
    send_count[rank] += 1;

  }

  MPI_Alltoallv(send_glist, send_count, send_shift, FVM_MPI_GNUM,
                recv_glist, recv_count, recv_shift, FVM_MPI_GNUM,
                mpi_comm);

  /* Send back to the original rank the new global vertex number */

  for (rank = 0; rank < n_ranks; rank++) {

    for (i = recv_shift[rank]; i < recv_shift[rank+1]; i++) {

      shift = recv_glist[i] - block_info.first_gnum;
      new_gnum = block_gnum[shift];
      recv_glist[i] = new_gnum;

    }

  } /* End of loop on ranks */

  MPI_Alltoallv(recv_glist, recv_count, recv_shift, FVM_MPI_GNUM,
                send_glist, send_count, send_shift, FVM_MPI_GNUM,
                mpi_comm);

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (i = 0; i < mesh->n_vertices; i++) {

    rank = (mesh->global_vtx_num[i] - 1)/block_info.size;
    shift = send_shift[rank] + send_count[rank];
    local_gnum[i] = send_glist[shift];  /* New global number */
    send_count[rank] += 1;

  }

  BFT_FREE(send_count);
  BFT_FREE(send_shift);
  BFT_FREE(send_glist);
  BFT_FREE(recv_glist);
  BFT_FREE(recv_count);
  BFT_FREE(recv_shift);
  BFT_FREE(block_gnum);

  /* Return pointer */

  *p_o2n_vtx_gnum = local_gnum;

}

/*----------------------------------------------------------------------------
 * Create arrays used to synchronize "single" elements from "coupled" elements
 *
 * parameters:
 *  sel               <--  list of all implied entities in the joining op.
 *  p_n_s_ranks       <->  number of ranks associated to single elements
 *  p_s_ranks         <->  list of "single" ranks
 *  p_s_rank2vtx_idx  <->  single ranks -> single vertices index
 *  p_s_rank2vtx_lst  <->  single ranks -> single vertices list
 *  p_n_c_ranks       <->  number of ranks associated to coupled elements
 *  p_c_ranks         <->  list of "coupled" ranks
 *  p_c_rank2vtx_idx  <->  coupled ranks -> coupled vertices index
 *  p_c_rank2vtx_lst  <->  coupled ranks -> coupled vertices list
 *---------------------------------------------------------------------------*/

static void
_prepare_single_sync(cs_join_select_t    *sel,
                     cs_int_t            *p_n_s_ranks,
                     cs_int_t            *p_s_ranks[],
                     cs_int_t            *p_s_rank2vtx_idx[],
                     cs_int_t            *p_s_rank2vtx_lst[],
                     cs_int_t            *p_n_c_ranks,
                     cs_int_t            *p_c_ranks[],
                     cs_int_t            *p_c_rank2vtx_idx[],
                     cs_int_t            *p_c_rank2vtx_lst[])
{
  cs_int_t  i, j, rank, shift;

  cs_int_t  n_s_ranks = 0, n_c_ranks = 0;
  cs_int_t  *count = NULL, *s_ranks = NULL, *c_ranks = NULL;
  cs_int_t  *s_rank2vtx_idx = NULL, *s_rank2vtx_lst = NULL;
  cs_int_t  *c_rank2vtx_idx = NULL, *c_rank2vtx_lst = NULL;

  const int  n_ranks = cs_glob_n_ranks;

  BFT_MALLOC(count, n_ranks, cs_int_t);

  for (i = 0; i < n_ranks; i++)
    count[i] = 0;

  /* Define s_ranks and n_s_ranks */

  for (i = 0; i < sel->s_vtx_idx[sel->n_s_vertices]; i++)
    count[sel->s_vtx_rank_lst[i]] += 1;

  for (i = 0; i < n_ranks; i++)
    if (count[i] > 0)
      n_s_ranks++;

  BFT_MALLOC(s_ranks, n_s_ranks, cs_int_t);

  n_s_ranks = 0;
  for (i = 0; i < n_ranks; i++)
    if (count[i] > 0)
      s_ranks[n_s_ranks++] = i;

  /* Define c_ranks and n_c_ranks */

  for (i = 0; i < n_ranks; i++)
    count[i] = 0;

  for (i = 0; i < sel->n_c_vertices; i++)
    count[sel->c_vtx_rank_lst[i]] += 1;

  for (i = 0; i < n_ranks; i++)
    if (count[i] > 0)
      n_c_ranks++;

  BFT_MALLOC(c_ranks, n_c_ranks, cs_int_t);

  n_c_ranks = 0;
  for (i = 0; i < n_ranks; i++)
    if (count[i] > 0)
      c_ranks[n_c_ranks++] = i;

  /* Invert sel->s_vtx_rank_idx */

  BFT_MALLOC(s_rank2vtx_idx, n_ranks + 1, cs_int_t);

  for (i = 0; i < n_ranks + 1; i++)
    s_rank2vtx_idx[i] = 0;

  for (i = 0; i < sel->n_s_vertices; i++)
    for (j = sel->s_vtx_idx[i]; j < sel->s_vtx_idx[i+1]; j++)
      s_rank2vtx_idx[sel->s_vtx_rank_lst[j] + 1] += 1;

  for (i = 0; i < n_ranks; i++)
    s_rank2vtx_idx[i+1] += s_rank2vtx_idx[i];

  BFT_MALLOC(s_rank2vtx_lst, s_rank2vtx_idx[n_ranks], cs_int_t);

  for (i = 0; i < n_ranks; i++)
    count[i] = 0;

  for (i = 0; i < sel->n_s_vertices; i++) {
    for (j = sel->s_vtx_idx[i]; j < sel->s_vtx_idx[i+1]; j++) {

      rank = sel->s_vtx_rank_lst[j];
      shift = s_rank2vtx_idx[rank] + count[rank];
      s_rank2vtx_lst[shift] = sel->s_vertices[i];
      count[rank] += 1;

    }
  }

  /* Invert sel->c_vtx_rank_lst */

  BFT_MALLOC(c_rank2vtx_idx, n_ranks + 1, cs_int_t);

  for (i = 0; i < n_ranks + 1; i++)
    c_rank2vtx_idx[i] = 0;

  for (i = 0; i < sel->n_c_vertices; i++)
    c_rank2vtx_idx[sel->c_vtx_rank_lst[i] + 1] += 1;

  for (i = 0; i < n_ranks; i++)
    c_rank2vtx_idx[i+1] += c_rank2vtx_idx[i];

  BFT_MALLOC(c_rank2vtx_lst, c_rank2vtx_idx[n_ranks], cs_int_t);

  for (i = 0; i < n_ranks; i++)
    count[i] = 0;

  for (i = 0; i < sel->n_c_vertices; i++) {

    rank = sel->c_vtx_rank_lst[i];
    shift = c_rank2vtx_idx[rank] + count[rank];
    c_rank2vtx_lst[shift] = sel->c_vertices[i];
    count[rank] += 1;

  }

  BFT_FREE(count);

  /* Return pointers */

  *p_n_s_ranks = n_s_ranks;
  *p_s_ranks = s_ranks;
  *p_s_rank2vtx_idx = s_rank2vtx_idx;
  *p_s_rank2vtx_lst = s_rank2vtx_lst;
  *p_n_c_ranks = n_c_ranks;
  *p_c_ranks = c_ranks;
  *p_c_rank2vtx_idx = c_rank2vtx_idx;
  *p_c_rank2vtx_lst = c_rank2vtx_lst;

}

/*----------------------------------------------------------------------------
 * Return true if the current couple of vertices is a "single" edge to sync.
 *
 * parameters:
 *  vid1           <--  first vertex id
 *  vid2           <--  second vertex id
 *  tag            <--  array on vertices (1 if vertex is "single" else 0)
 *  edge_builder   <->  pointer to an edge_builder_t structure
 *---------------------------------------------------------------------------*/

inline static cs_bool_t
_is_s_edge(cs_int_t          vid1,
           cs_int_t          vid2,
           cs_int_t          tag[],
           edge_builder_t   *edge_builder)
{
  cs_int_t  e_id;

  if (tag[vid1] == 1 || tag[vid2] == 1) {

    e_id = _get_join_edge_id(vid1, vid2, edge_builder);
    assert(e_id != -1);

    if (edge_builder->v2v_sub_idx[e_id+1]-edge_builder->v2v_sub_idx[e_id] > 0)
      return false;
    else
      return true;

  }
  else
    return false;

}

/*----------------------------------------------------------------------------
 * Update elements of a cs_mesh_t structure related to the vertices after the
 * merge step.
 *
 * parameters:
 *  selection         <--  list of all implied entities in the joining op.
 *  n_bm_vertices     <--  number of vertices in mesh before the fusion step
 *  old_vtx_gnum      <--  old global vertex numbering
 *  o2n_vtx_id        <->  relation between init. and current local num.
 *  n_j_vertices      <--  number of vertices in join_mesh
 *  join2mesh_vtx_id  <->  relation between join mesh and after fusion vertex
 *  edge_builder      <->  pointer to an edge_builder_t structure
 *  mesh              <->  pointer of pointer to cs_mesh_t structure
 *---------------------------------------------------------------------------*/

static void
_sync_single_elements(cs_join_select_t    *selection,
                      cs_int_t             n_bm_vertices,
                      fvm_gnum_t           old_vtx_gnum[],
                      cs_int_t             o2n_vtx_id[],
                      cs_int_t             n_j_vertices,
                      cs_int_t             join2mesh_vtx_id[],
                      edge_builder_t      *edge_builder,
                      cs_mesh_t           *mesh)
{
  cs_int_t  i, j, s, e, id, vid, fid, rank, shift, length, request_count;

  cs_int_t  n_s_ranks = 0, n_c_ranks = 0;
  cs_int_t  *selection_tag = NULL, *s_ranks = NULL, *c_ranks = NULL;
  cs_int_t  *s_rank2vtx_idx = NULL, *s_rank2vtx_lst = NULL;
  cs_int_t  *c_rank2vtx_idx = NULL, *c_rank2vtx_lst = NULL;
  cs_int_t  *new_v2v_sub_idx = NULL, *new_v2v_sub_lst = NULL;
  fvm_gnum_t  *new_vtx_gnum = NULL;
  cs_real_t  *new_coord = NULL;

  MPI_Request  *request = NULL;
  MPI_Status   *status = NULL;
  MPI_Comm  mpi_comm = cs_glob_mpi_comm;

  const int  loc_rank = cs_glob_rank_id;
  const int  n_ranks = cs_glob_n_ranks;

  bft_printf("\n  Synchronization of the \"single\" elements after the merge"
             " step.\n");
  bft_printf_flush();

  assert(n_ranks > 1);

  _prepare_single_sync(selection,
                       &n_s_ranks,
                       &s_ranks,
                       &s_rank2vtx_idx,
                       &s_rank2vtx_lst,
                       &n_c_ranks,
                       &c_ranks,
                       &c_rank2vtx_idx,
                       &c_rank2vtx_lst);

  /* Allocate MPI buffers used for exchanging data */

  BFT_MALLOC(request, n_c_ranks + n_s_ranks, MPI_Request);
  BFT_MALLOC(status, n_c_ranks + n_s_ranks, MPI_Status);
  BFT_MALLOC(selection_tag, mesh->n_vertices, cs_int_t);

  /* Define a selection tag to find quickly "single" vertices */

  for (i = 0; i < mesh->n_vertices; i++)
    selection_tag[i] = 0;
  for (i = 0; i < selection->n_s_vertices; i++)
    selection_tag[selection->s_vertices[i]-1] = 1;

  {  /* Synchronization of vertex coordinates */

    double  *s_buf = NULL, *c_buf = NULL;

    BFT_MALLOC(s_buf, 3*s_rank2vtx_idx[n_ranks], double);

    /* Receive data from distant ranks */

    request_count = 0;

    for (i = 0; i < n_s_ranks; i++) {

      rank = s_ranks[i];
      s = s_rank2vtx_idx[rank];
      e = s_rank2vtx_idx[rank+1];
      length = 3*(e-s);

      MPI_Irecv(&(s_buf[3*s]), length, MPI_DOUBLE,
                rank, rank, mpi_comm, &(request[request_count++]));

    }

    /* We wait for posting all receives (often recommended) */

    MPI_Barrier(mpi_comm);

    /* Build c_buf = buffer to send */

    BFT_MALLOC(c_buf, 3*c_rank2vtx_idx[n_ranks], double);

    for (shift = 0, i = 0; i < c_rank2vtx_idx[n_ranks]; i++) {

      int  new_id = o2n_vtx_id[c_rank2vtx_lst[i]-1];

      c_buf[shift++] = mesh->vtx_coord[3*new_id];
      c_buf[shift++] = mesh->vtx_coord[3*new_id+1];
      c_buf[shift++] = mesh->vtx_coord[3*new_id+2];

    }

    /* Send data to distant ranks */

    for (i = 0; i < n_c_ranks; i++) {

      rank = c_ranks[i];
      s = c_rank2vtx_idx[rank];
      e = c_rank2vtx_idx[rank+1];
      length = 3*(e-s);

      MPI_Isend(&(c_buf[3*s]), length, MPI_DOUBLE,
                rank, loc_rank, mpi_comm, &(request[request_count++]));

    }

    /* Wait for all exchanges */

    MPI_Waitall(request_count, request, status);

    /* Update vertex coordinates */

    for (i = 0; i < s_rank2vtx_idx[n_ranks]; i++) {

      int  new_id = o2n_vtx_id[s_rank2vtx_lst[i]-1];

      mesh->vtx_coord[3*new_id]   = s_buf[shift++];
      mesh->vtx_coord[3*new_id+1] = s_buf[shift++];
      mesh->vtx_coord[3*new_id+2] = s_buf[shift++];

    }

    BFT_FREE(c_buf);
    BFT_FREE(s_buf);

  } /* End of vertex coordinates synchronization */

  /* Get single edges (at least connected to a single vertices.
     Loop on single border and interior faces. */

  {
    int  k, edge_id, vid1, vid2;
    fvm_gnum_t  cur, prev;

    int  n_s_edges = 0, c_sub_size = 0, s_sub_size = 0, n_new_vertices = 0;
    cs_int_t  *sub_elt_count = NULL;
    int  *s_count = NULL, *s_shift = NULL, *s_buf = NULL;
    int  *c_count = NULL, *c_shift = NULL, *c_buf = NULL;
    fvm_gnum_t  *s_edge_def = NULL, *s_sub_edge_def = NULL, *c_gbuf = NULL;
    fvm_gnum_t  *c_sub_gbuf = NULL, *s_sub_gbuf = NULL;
    cs_real_t  *c_sub_coord = NULL, *s_sub_coord = NULL;
    cs_real_t  *s_sub_coord_def = NULL;

    BFT_MALLOC(s_count, n_s_ranks, int); /* interior + border */
    BFT_MALLOC(c_count, n_c_ranks, int);

    /* Send number of edges for border and interior faces */

    for (i = 0; i < n_s_ranks; i++)
      s_count[i] = 0;
    for (i = 0; i < n_c_ranks; i++)
      c_count[i] = 0;

    for (i = 0; i < selection->n_b_s_faces; i++) {

      fid = selection->b_s_faces[i] - 1;
      s = mesh->b_face_vtx_idx[fid] - 1;
      e = mesh->b_face_vtx_idx[fid+1] - 1;

      for (j = s; j < e - 1; j++) {

        vid1 = mesh->b_face_vtx_lst[j] - 1;
        vid2 = mesh->b_face_vtx_lst[j+1] - 1;

        if (_is_s_edge(vid1, vid2, selection_tag, edge_builder))
          n_s_edges += 1;

      }

      vid1 = mesh->b_face_vtx_lst[e-1] - 1;
      vid2 = mesh->b_face_vtx_lst[s] - 1;

      if (_is_s_edge(vid1, vid2, selection_tag, edge_builder))
        n_s_edges += 1;

    }

    for (i = 0; i < selection->n_i_s_faces; i++) {

      fid = selection->i_s_faces[i] - 1;
      s = mesh->i_face_vtx_idx[fid] - 1;
      e = mesh->i_face_vtx_idx[fid+1] - 1;

      for (j = s; j < e - 1; j++) {

        vid1 = mesh->i_face_vtx_lst[j] - 1;
        vid2 = mesh->i_face_vtx_lst[j+1] - 1;

        if (_is_s_edge(vid1, vid2, selection_tag, edge_builder))
          n_s_edges += 1;

      }

      vid1 = mesh->i_face_vtx_lst[e-1] - 1;
      vid2 = mesh->i_face_vtx_lst[s] - 1;

      if (_is_s_edge(vid1, vid2, selection_tag, edge_builder))
        n_s_edges += 1;

    }

    request_count = 0;

    for (i = 0; i < n_c_ranks; i++)
      MPI_Irecv(&(c_count[i]), 1, MPI_INT, c_ranks[i], c_ranks[i],
                mpi_comm, &(request[request_count++]));

    /* We wait for posting all receives (often recommended) */

    MPI_Barrier(mpi_comm);

    /* Send data to distant ranks */

    for (i = 0; i < n_s_ranks; i++)
      MPI_Isend(&n_s_edges, 1, MPI_INT,
                s_ranks[i], loc_rank, mpi_comm, &(request[request_count++]));

    /* Wait for all exchanges */

    MPI_Waitall(request_count, request, status);

    /* Send edge definition */

    BFT_MALLOC(s_edge_def, 2*n_s_edges, fvm_gnum_t);

    /* Define edges for border faces */

    shift = 0;
    for (i = 0; i < selection->n_b_s_faces; i++) {

      fid = selection->b_s_faces[i] - 1;
      s = mesh->b_face_vtx_idx[fid] - 1;
      e = mesh->b_face_vtx_idx[fid+1] - 1;

      for (j = s; j < e - 1; j++) {

        vid1 = mesh->b_face_vtx_lst[j] - 1;
        vid2 = mesh->b_face_vtx_lst[j+1] - 1;

        if (_is_s_edge(vid1, vid2, selection_tag, edge_builder)) {
          s_edge_def[shift++] = old_vtx_gnum[vid1];
          s_edge_def[shift++] = old_vtx_gnum[vid2];
        }

      }

      vid1 = mesh->b_face_vtx_lst[e-1] - 1;
      vid2 = mesh->b_face_vtx_lst[s] - 1;

      if (_is_s_edge(vid1, vid2, selection_tag, edge_builder)) {
        s_edge_def[shift++] = old_vtx_gnum[vid1];
        s_edge_def[shift++] = old_vtx_gnum[vid2];
      }

    }

    /* Define edges for interior faces */

    for (i = 0; i < selection->n_i_s_faces; i++) {

      fid = selection->i_s_faces[i] - 1;
      s = mesh->i_face_vtx_idx[fid] - 1;
      e = mesh->i_face_vtx_idx[fid+1] - 1;

      for (j = s; j < e - 1; j++) {

        vid1 = mesh->i_face_vtx_lst[j] - 1;
        vid2 = mesh->i_face_vtx_lst[j+1] - 1;

        if (_is_s_edge(vid1, vid2, selection_tag, edge_builder)) {
          s_edge_def[shift++] = old_vtx_gnum[vid1];
          s_edge_def[shift++] = old_vtx_gnum[vid2];
        }

      }

      vid1 = mesh->i_face_vtx_lst[e-1] - 1;
      vid2 = mesh->i_face_vtx_lst[s] - 1;

      if (_is_s_edge(vid1, vid2, selection_tag, edge_builder)) {
        s_edge_def[shift++] = old_vtx_gnum[vid1];
        s_edge_def[shift++] = old_vtx_gnum[vid2];
      }

    }

    assert(shift == 2*n_s_edges);

    BFT_MALLOC(c_shift, n_c_ranks + 1, int);

    c_shift[0] = 0;
    for (rank = 0; rank < n_c_ranks; rank++)
      c_shift[rank+1] = c_shift[rank] + c_count[rank];

    BFT_MALLOC(c_gbuf, 2*c_shift[n_c_ranks], fvm_gnum_t);

    /* Exchange edge definition */

    request_count = 0;

    for (i = 0; i < n_c_ranks; i++) {
      rank = c_ranks[i];
      length = 2 * c_count[i];
      MPI_Irecv(&(c_gbuf[2*c_shift[i]]), length, FVM_MPI_GNUM,
                rank, rank, mpi_comm, &(request[request_count++]));
    }

    /* We wait for posting all receives (often recommended) */

    MPI_Barrier(mpi_comm);

    /* Send data to distant ranks */

    for (i = 0; i < n_s_ranks; i++) {
      rank = s_ranks[i];
      MPI_Isend(s_edge_def, 2*n_s_edges, FVM_MPI_GNUM,
                rank, loc_rank, mpi_comm, &(request[request_count++]));
    }

    /* Wait for all exchanges */

    MPI_Waitall(request_count, request, status);

    /* Get a number of sub-element for each received edge */

    BFT_MALLOC(c_buf, c_shift[n_c_ranks], int);

    for (rank = 0; rank < n_c_ranks; rank++) {

      for (i = c_shift[rank]; i < c_shift[rank+1]; i++) {

        fvm_gnum_t  v1_gnum = c_gbuf[2*i];
        fvm_gnum_t  v2_gnum = c_gbuf[2*i+1];

        /* Search for v1_gnum and v2_gnum in global_vtx_num */

        vid1 = cs_search_g_binary(0, n_bm_vertices-1,
                                  v1_gnum,
                                  old_vtx_gnum);

        if (vid1 > -1) {

          vid2 = cs_search_g_binary(0, n_bm_vertices-1,
                                    v2_gnum,
                                    old_vtx_gnum);

          if (vid2 > -1) {

            edge_id = _get_join_edge_id(vid1, vid2, edge_builder);

            if (edge_id != -1)
              c_buf[i] =  edge_builder->v2v_sub_idx[edge_id+1]
                        - edge_builder->v2v_sub_idx[edge_id];
            else
              c_buf[i] = 0; /* Not in edge_builder_t structure */

          }
          else
            c_buf[i] = 0;

        }
        else
          c_buf[i] = 0;

      } /* End of loop on edges received from rank */

    } /* End of loop on "coupled" ranks */

    /* Exchange number of sub elements for each edge */

    BFT_MALLOC(sub_elt_count, n_s_edges, int);
    BFT_MALLOC(s_buf, n_s_edges * n_s_ranks, int);

    for (i = 0; i < n_s_edges; i++)
      sub_elt_count[i] = 0;

    request_count = 0;

    for (i = 0; i < n_s_ranks; i++)
      MPI_Irecv(&(s_buf[i*n_s_edges]), n_s_edges, MPI_INT,
                s_ranks[i], s_ranks[i], mpi_comm, &(request[request_count++]));

    /* We wait for posting all receives (often recommended) */

    MPI_Barrier(mpi_comm);

    /* Send data to distant ranks */

    for (i = 0; i < n_c_ranks; i++)
      MPI_Isend(&(c_buf[c_shift[i]]), c_count[i], MPI_INT,
                c_ranks[i], loc_rank, mpi_comm, &(request[request_count++]));

    /* Wait for all exchanges */

    MPI_Waitall(request_count, request, status);

    /* We assume that the sub-edge definition is identical for all ranks where
       a sub-edge definition exists */

    for (rank = 0; rank < n_s_ranks; rank++) {
      for (i = 0; i < n_s_edges; i++)
        if (s_buf[rank*n_s_edges + i] > 0)
          sub_elt_count[i] = s_buf[rank*n_s_edges + i];
    }

    /* Define buffer to send whith sub-elements from "coupled" ranks to
       "single" ranks */

    for (i = 0; i < c_shift[n_c_ranks]; i++)
      c_sub_size += c_buf[i];

    BFT_MALLOC(c_sub_gbuf, c_sub_size, fvm_gnum_t);
    BFT_MALLOC(c_sub_coord, 3*c_sub_size, cs_real_t);

    shift = 0;
    for (rank = 0; rank < n_c_ranks; rank++) {

      for (i = c_shift[rank]; i < c_shift[rank+1]; i++) {

        fvm_gnum_t  v1_gnum = c_gbuf[2*i];
        fvm_gnum_t  v2_gnum = c_gbuf[2*i+1];

        vid1 = cs_search_g_binary(0, n_bm_vertices-1,
                                  v1_gnum,
                                  old_vtx_gnum);

        if (vid1 > -1) {

          vid2 = cs_search_g_binary(0, n_bm_vertices-1,
                                    v2_gnum,
                                    old_vtx_gnum);

          if (vid2 > -1) {

            edge_id = _get_join_edge_id(vid1, vid2, edge_builder);

            if (edge_id != -1) {

              for (j = edge_builder->v2v_sub_idx[edge_id];
                   j < edge_builder->v2v_sub_idx[edge_id+1]; j++) {

                vid = edge_builder->v2v_sub_lst[j] - 1;
                c_sub_gbuf[shift] = mesh->global_vtx_num[vid];
                for (k = 0; k < 3; k++)
                  c_sub_coord[3*shift+k] = mesh->vtx_coord[3*vid+k];
                shift++;

              }

            } /* edge_id != -1 */

          } /* vid2 > -1 */

        } /* vid1 > -1 */

      } /* End of loop on edges received from rank */

    } /* End of loop on ranks */

    BFT_MALLOC(s_shift, n_s_ranks + 1, cs_int_t);

    s_shift[0] = 0;
    for (rank = 0; rank < n_s_ranks; rank++) {
      for (i = 0; i < n_s_edges; i++)
        s_count[rank] += s_buf[rank*n_s_edges + i];
      s_shift[rank+1] = s_shift[rank] + s_count[rank];
    }

    BFT_MALLOC(s_sub_gbuf, s_shift[n_s_ranks], fvm_gnum_t);
    BFT_MALLOC(s_sub_coord, 3*s_shift[n_s_ranks], cs_real_t);

    s_sub_size = 0;
    for (i = 0; i < n_s_edges; i++)
      s_sub_size += sub_elt_count[i];

    BFT_MALLOC(s_sub_edge_def, s_sub_size, fvm_gnum_t);
    BFT_MALLOC(s_sub_coord_def, 3*s_sub_size, cs_real_t);

    /* Exchange sub-edge definition: global vertex number  */

    request_count = 0;

    for (i = 0; i < n_s_ranks; i++)
      MPI_Irecv(&(s_sub_gbuf[s_shift[i]]), s_count[i], FVM_MPI_GNUM,
                s_ranks[i], s_ranks[i], mpi_comm, &(request[request_count++]));

    /* We wait for posting all receives (often recommended) */

    MPI_Barrier(mpi_comm);

    /* Send data to distant ranks */

    c_sub_size = 0;
    for (i = 0; i < n_c_ranks; i++) {

      c_count[i] = 0;
      for (j = c_shift[i]; j < c_shift[i+1]; j++)
        c_count[i] += c_buf[j];

      MPI_Isend(&(c_sub_gbuf[c_sub_size]), c_count[i], FVM_MPI_GNUM,
                c_ranks[i], loc_rank, mpi_comm, &(request[request_count++]));

      c_sub_size += c_count[i];

    } /* End of loop on ranks */

    /* Wait for all exchanges */

    MPI_Waitall(request_count, request, status);

    /* Exchange sub-edge definition: vertex coordinates */

    request_count = 0;

    for (i = 0; i < n_s_ranks; i++)
      MPI_Irecv(&(s_sub_coord[3*s_shift[i]]), 3*s_count[i], FVM_MPI_COORD,
                s_ranks[i], s_ranks[i], mpi_comm, &(request[request_count++]));

    /* We wait for posting all receives (often recommended) */

    MPI_Barrier(mpi_comm);

    /* Send data to distant ranks */

    c_sub_size = 0;
    for (i = 0; i < n_c_ranks; i++) {

      c_count[i] = 0;
      for (j = c_shift[i]; j < c_shift[i+1]; j++)
        c_count[i] += 3*c_buf[j];

      MPI_Isend(&(c_sub_coord[c_sub_size]), c_count[i], FVM_MPI_COORD,
                c_ranks[i], loc_rank, mpi_comm, &(request[request_count++]));

      c_sub_size += c_count[i];

    } /* End of loop on ranks */

    /* Wait for all exchanges */

    MPI_Waitall(request_count, request, status);

    /* Free memory */

    BFT_FREE(c_count);
    BFT_FREE(c_shift);
    BFT_FREE(c_buf);
    BFT_FREE(c_gbuf);
    BFT_FREE(c_sub_gbuf);
    BFT_FREE(c_sub_coord);
    BFT_FREE(s_count);

    if (n_s_edges > 0) {

      cs_int_t  n_sub_elts, e_shift;

      cs_int_t  *inv_order = NULL;
      fvm_lnum_t  *order = NULL;

      /* Define s_sub_edge_def from the received s_sub_gbuf */

      for (rank = 0; rank < n_s_ranks; rank++) {

        e_shift = 0;
        shift = s_shift[rank];

        for (i = 0; i < n_s_edges; i++) {

          n_sub_elts = s_buf[rank*n_s_edges + i];

          if (sub_elt_count[i] == n_sub_elts) {

            for (j = 0; j < n_sub_elts; j++) {
              s_sub_edge_def[e_shift + j] = s_sub_gbuf[shift + j];
              for (k = 0; k < 3; k++)
                s_sub_coord_def[3*(e_shift+j)+k] = s_sub_coord[3*(shift+j)+k];
            }

          }
          e_shift += sub_elt_count[i];
          shift += n_sub_elts;

        } /* End of loop on edges */

        assert(e_shift == s_sub_size);

      } /* End of loop on ranks */

      BFT_FREE(s_sub_gbuf);
      BFT_FREE(s_sub_coord);
      BFT_FREE(s_count);
      BFT_FREE(s_shift);
      BFT_FREE(s_buf);

      /* Update vertices from list of sub elements. Define new vertices */

      BFT_MALLOC(order, s_sub_size, fvm_lnum_t);

      fvm_order_local_allocated(NULL, s_sub_edge_def, order, s_sub_size);

      prev = 0;
      n_new_vertices = 0;

      for (i = 0; i < s_sub_size; i++) {

        cur = s_sub_edge_def[order[i]];
        if (cur != prev) {
          prev = cur;
          id = cs_search_g_binary(0, mesh->n_vertices-1,
                                  cur,
                                  mesh->global_vtx_num);
          if (id == -1) /* Add vertex */
            n_new_vertices++;
        }

      }

      BFT_REALLOC(mesh->global_vtx_num,
                  mesh->n_vertices + n_new_vertices, fvm_gnum_t);

      BFT_REALLOC(mesh->vtx_coord,
                  3*(mesh->n_vertices + n_new_vertices), cs_real_t);

      prev = 0;
      n_new_vertices = 0;

      for (i = 0; i < s_sub_size; i++) {

        cur = s_sub_edge_def[order[i]];
        if (cur != prev) {

          prev = cur;
          id = cs_search_g_binary(0, mesh->n_vertices-1,
                                  cur,
                                  mesh->global_vtx_num);

          if (id == -1) { /* Add vertex */

            shift = mesh->n_vertices + n_new_vertices;
            mesh->global_vtx_num[shift] = cur;
            for (k = 0; k < 3; k++)
              mesh->vtx_coord[3*shift+k] = s_sub_coord_def[3*order[i]+k];
            n_new_vertices++;

          }

        }

      }

      mesh->n_vertices += n_new_vertices;
      bft_printf(_("  Add %d new vertices from the single elements sync.\n"),
                 n_new_vertices);

      /* Reorder global_vtx_num in order to have an ordered list */

      BFT_REALLOC(order, mesh->n_vertices, fvm_lnum_t);

      fvm_order_local_allocated(NULL,
                                mesh->global_vtx_num,
                                order,
                                mesh->n_vertices);

      BFT_MALLOC(new_vtx_gnum, mesh->n_vertices, fvm_gnum_t);

      for (i = 0; i < mesh->n_vertices; i++)
        new_vtx_gnum[i] = mesh->global_vtx_num[order[i]];

      BFT_FREE(mesh->global_vtx_num);
      mesh->global_vtx_num = new_vtx_gnum;

      /* Define a new mesh->vtx_coord */

      BFT_MALLOC(new_coord, 3*mesh->n_vertices, cs_real_t);

      for (i = 0; i < mesh->n_vertices; i++) {
        vid = order[i];
        for (k = 0; k < 3; k++)
          new_coord[3*i+k] = mesh->vtx_coord[3*vid+k];
      }

      BFT_FREE(mesh->vtx_coord);
      mesh->vtx_coord = new_coord;

      /* Define a new o2n_vtx_id and o2n_vtx_gnum */

      BFT_MALLOC(inv_order, mesh->n_vertices, cs_int_t);

      for (i = 0; i < mesh->n_vertices; i++) {
        j = order[i];
        inv_order[j] = i;
      }

      for (i = 0; i < n_bm_vertices; i++) {
        vid = o2n_vtx_id[i];
        o2n_vtx_id[i] = inv_order[vid];
      }

      /* Define a new join2mesh_vtx_id */

      for (i = 0; i < n_j_vertices; i++) {
        vid = join2mesh_vtx_id[i];
        join2mesh_vtx_id[i] = inv_order[vid];
      }

      /* Update edge_builder_t structure (v2v_sub_idx and v2v_sub_lst) */

      BFT_MALLOC(new_v2v_sub_idx, edge_builder->n_edges + 1, cs_int_t);

      for (i = 0; i < edge_builder->n_edges; i++)
        new_v2v_sub_idx[i+1] =  edge_builder->v2v_sub_idx[i+1]
                              - edge_builder->v2v_sub_idx[i];

      /* Update index for single border faces */

      e_shift = 0;
      for (i = 0; i < selection->n_b_s_faces; i++) {

        fid = selection->b_s_faces[i] - 1;
        s = mesh->b_face_vtx_idx[fid] - 1;
        e = mesh->b_face_vtx_idx[fid+1] - 1;

        for (j = s; j < e - 1; j++) {

          vid1 = mesh->b_face_vtx_lst[j] - 1;
          vid2 = mesh->b_face_vtx_lst[j+1] - 1;

          if (_is_s_edge(vid1, vid2, selection_tag, edge_builder)) {
            edge_id = _get_join_edge_id(vid1, vid2, edge_builder);
            new_v2v_sub_idx[edge_id+1] = sub_elt_count[e_shift++];
          }

        }

        vid1 = mesh->b_face_vtx_lst[e-1] - 1;
        vid2 = mesh->b_face_vtx_lst[s] - 1;

        if (_is_s_edge(vid1, vid2, selection_tag, edge_builder)) {
          edge_id = _get_join_edge_id(vid1, vid2, edge_builder);
          new_v2v_sub_idx[edge_id+1] = sub_elt_count[e_shift++];
        }

      }

      /* Update index for single interior faces */

      for (i = 0; i < selection->n_i_s_faces; i++) {

        fid = selection->i_s_faces[i] - 1;
        s = mesh->i_face_vtx_idx[fid] - 1;
        e = mesh->i_face_vtx_idx[fid+1] - 1;

        for (j = s; j < e - 1; j++) {

          vid1 = mesh->i_face_vtx_lst[j] - 1;
          vid2 = mesh->i_face_vtx_lst[j+1] - 1;

          if (_is_s_edge(vid1, vid2, selection_tag, edge_builder)) {
            edge_id = _get_join_edge_id(vid1, vid2, edge_builder);
            new_v2v_sub_idx[edge_id+1] = sub_elt_count[e_shift++];
          }

        }

        vid1 = mesh->i_face_vtx_lst[e-1] - 1;
        vid2 = mesh->i_face_vtx_lst[s] - 1;

        if (_is_s_edge(vid1, vid2, selection_tag, edge_builder)) {
          edge_id = _get_join_edge_id(vid1, vid2, edge_builder);
          new_v2v_sub_idx[edge_id+1] = sub_elt_count[e_shift++];
        }

      }

      new_v2v_sub_idx[0] = 0;
      for (i = 0; i < edge_builder->n_edges; i++)
        new_v2v_sub_idx[i+1] += new_v2v_sub_idx[i];

      /* Update v2v_sub_lst */

      BFT_MALLOC(new_v2v_sub_lst,
                 new_v2v_sub_idx[edge_builder->n_edges], cs_int_t);

      for (i = 0; i < edge_builder->n_edges; i++) {

        shift = new_v2v_sub_idx[i];

        for (j = edge_builder->v2v_sub_idx[i];
             j < edge_builder->v2v_sub_idx[i+1]; j++) {
          vid = edge_builder->v2v_sub_lst[j] - 1;
          new_v2v_sub_lst[shift++] = inv_order[vid] + 1;
        }

      }

      /* Update sub list for single border faces */

      e_shift = 0;
      shift = 0;

      for (i = 0; i < selection->n_b_s_faces; i++) {

        fid = selection->b_s_faces[i] - 1;
        s = mesh->b_face_vtx_idx[fid] - 1;
        e = mesh->b_face_vtx_idx[fid+1] - 1;

        for (j = s; j < e - 1; j++) {

          vid1 = mesh->b_face_vtx_lst[j] - 1;
          vid2 = mesh->b_face_vtx_lst[j+1] - 1;

          if (_is_s_edge(vid1, vid2, selection_tag, edge_builder)) {

            edge_id = _get_join_edge_id(vid1, vid2, edge_builder);

            for (k = new_v2v_sub_idx[edge_id];
                 k < new_v2v_sub_idx[edge_id+1]; k++) {

              id = cs_search_g_binary(0, mesh->n_vertices-1,
                                      s_sub_edge_def[shift++],
                                      mesh->global_vtx_num);

              assert(id != -1);
              new_v2v_sub_lst[k] = id + 1;

            }

          }

        }  /* End of loop face connect */

        vid1 = mesh->b_face_vtx_lst[e-1] - 1;
        vid2 = mesh->b_face_vtx_lst[s] - 1;

        if (_is_s_edge(vid1, vid2, selection_tag, edge_builder)) {

          edge_id = _get_join_edge_id(vid1, vid2, edge_builder);

          for (k = new_v2v_sub_idx[edge_id];
               k < new_v2v_sub_idx[edge_id+1]; k++) {

            id = cs_search_g_binary(0, mesh->n_vertices-1,
                                    s_sub_edge_def[shift++],
                                    mesh->global_vtx_num);

            assert(id != -1);
            new_v2v_sub_lst[k] = id + 1;

          }

        }

      } /* End of loop on single border faces */

      /* Update sub list for single interior faces */

      for (i = 0; i < selection->n_i_s_faces; i++) {

        fid = selection->i_s_faces[i] - 1;
        s = mesh->i_face_vtx_idx[fid] - 1;
        e = mesh->i_face_vtx_idx[fid+1] - 1;

        for (j = s; j < e - 1; j++) {

          vid1 = mesh->i_face_vtx_lst[j] - 1;
          vid2 = mesh->i_face_vtx_lst[j+1] - 1;

          if (_is_s_edge(vid1, vid2, selection_tag, edge_builder)) {

            edge_id = _get_join_edge_id(vid1, vid2, edge_builder);

            for (k = new_v2v_sub_idx[edge_id];
                 k < new_v2v_sub_idx[edge_id+1]; k++) {

              id = cs_search_g_binary(0, mesh->n_vertices-1,
                                      s_sub_edge_def[shift++],
                                      mesh->global_vtx_num);

              assert(id != -1);
              new_v2v_sub_lst[k] = id + 1;

            }

          }

        } /* End of loop face connect */

        vid1 = mesh->i_face_vtx_lst[e-1] - 1;
        vid2 = mesh->i_face_vtx_lst[s] - 1;

        if (_is_s_edge(vid1, vid2, selection_tag, edge_builder)) {

          edge_id = _get_join_edge_id(vid1, vid2, edge_builder);

          for (k = new_v2v_sub_idx[edge_id];
               k < new_v2v_sub_idx[edge_id+1]; k++) {

            id = cs_search_g_binary(0, mesh->n_vertices-1,
                                    s_sub_edge_def[shift++],
                                    mesh->global_vtx_num);

            assert(id != -1);
            new_v2v_sub_lst[k] = id + 1;

          }

        }

      } /* End of loop on single interior faces */

      BFT_FREE(order);
      BFT_FREE(inv_order);
      BFT_FREE(edge_builder->v2v_sub_idx);
      BFT_FREE(edge_builder->v2v_sub_lst);

      edge_builder->v2v_sub_idx = new_v2v_sub_idx;
      edge_builder->v2v_sub_lst = new_v2v_sub_lst;

    } /* End if n_s_edges > 0 */

    /* Free memory */

    BFT_FREE(sub_elt_count);
    BFT_FREE(s_sub_edge_def);
    BFT_FREE(s_sub_coord_def);
    BFT_FREE(s_edge_def);

  } /* Update edge_builder structure */

  /* Free memory */

  BFT_FREE(s_rank2vtx_idx);
  BFT_FREE(s_rank2vtx_lst);
  BFT_FREE(s_ranks);
  BFT_FREE(c_rank2vtx_idx);
  BFT_FREE(c_rank2vtx_lst);
  BFT_FREE(c_ranks);
  BFT_FREE(request);
  BFT_FREE(status);
  BFT_FREE(selection_tag);

}
#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Update elements of a cs_mesh_t structure implied in the joining
 * operation. Update the vertices after the merge step.
 *
 * parameters:
 *  o2n_vtx_gnum       <--  local array on the new global vertex
 *  join_mesh          <--  pointer to a local cs_join_mesh_t struct.
 *  mesh               <->  pointer of pointer to cs_mesh_t structure
 *  p_join2mesh_vtx_id -->  relation between join mesh and after vertex merge
 *  p_o2n_vtx_id       -->  relation between init. and current local num.
 *  p_old_vtx_gnum     -->  old global vertex numbering
 *---------------------------------------------------------------------------*/

static void
_update_vertices_after_merge(fvm_gnum_t         o2n_vtx_gnum[],
                             cs_join_mesh_t    *join_mesh,
                             cs_mesh_t         *mesh,
                             cs_int_t          *p_join2mesh_vtx_id[],
                             cs_int_t          *p_o2n_vtx_id[],
                             fvm_gnum_t        *p_old_vtx_gnum[])
{
  cs_int_t  i, j, k, o_id, j_id;
  fvm_gnum_t  prev, cur;

  cs_int_t  n_am_vertices = -1; /* af: after merge */
  cs_real_t  *new_vtx_coord = NULL;
  cs_int_t   *o2n_vtx_id = NULL, *join2mesh_vtx_id = NULL;
  fvm_lnum_t  *order = NULL;
  fvm_gnum_t  *new_vtx_gnum = NULL, *tmp_vtx_gnum = NULL;

  const cs_int_t  n_bm_vertices = mesh->n_vertices; /* bm: before merge */
  const cs_int_t  j_n_vertices = join_mesh->n_vertices;
  const cs_join_vertex_t  *j_vertices = join_mesh->vertices;
  const cs_int_t  n_vertices = n_bm_vertices + j_n_vertices;
  const int  n_ranks = cs_glob_n_ranks;

  /* Update initial vertices (local and global numbering) */

  BFT_MALLOC(o2n_vtx_id, n_bm_vertices, cs_int_t);
  BFT_MALLOC(join2mesh_vtx_id, j_n_vertices, cs_int_t);
  BFT_MALLOC(tmp_vtx_gnum, n_vertices, fvm_gnum_t);
  BFT_MALLOC(new_vtx_gnum, n_vertices, fvm_gnum_t);
  BFT_MALLOC(order, n_vertices, cs_int_t);

  for (i = 0; i < n_bm_vertices; i++)
    tmp_vtx_gnum[i] = o2n_vtx_gnum[i];

  for (i = 0, j = n_bm_vertices; i < j_n_vertices; i++, j++)
    tmp_vtx_gnum[j] = j_vertices[i].gnum;

  fvm_order_local_allocated(NULL, tmp_vtx_gnum, order, n_vertices);

  /* Define o2n_vtx_id and join2mesh_vtx_id arrays */

  if (order[0] < n_bm_vertices)
    prev = o2n_vtx_gnum[order[0]] + 1;
  else
    prev = j_vertices[order[0]-n_bm_vertices].gnum + 1;

  for (i = 0; i < n_vertices; i++) {

    o_id = order[i];

    if (o_id < n_bm_vertices) { /* Belongs to the initial mesh */

      cur = o2n_vtx_gnum[o_id];
      assert(cur <= mesh->n_g_vertices);

      if (cur != prev) {

        n_am_vertices++;
        prev = cur;
        o2n_vtx_id[o_id] = n_am_vertices;
        new_vtx_gnum[n_am_vertices] = cur;

      }
      else
        o2n_vtx_id[o_id] = n_am_vertices;

    }
    else { /* Belongs to the join_mesh */

      j_id = o_id - n_bm_vertices;
      cur = j_vertices[j_id].gnum;

      if (cur != prev) {

        n_am_vertices++;
        prev = cur;
        new_vtx_gnum[n_am_vertices] = cur;
        join2mesh_vtx_id[j_id] = n_am_vertices;

      }
      else
        join2mesh_vtx_id[j_id] = n_am_vertices;

    }

  } /* End of loop on vertices */

  /* n_am_vertices was up to now an id. Move to a number */
  n_am_vertices++;

  /* Partial free memory */

  BFT_FREE(tmp_vtx_gnum);
  BFT_FREE(order);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n\n Dump Old2New array (local mesh): "
             "old_n_vertices = %d - new_n_vertices = %d\n",
             n_bm_vertices, n_am_vertices);
  for (i = 0; i < n_bm_vertices; i++)
    bft_printf("Old num : %7d (%9u) => New num : %7d (%9u)\n",
               i+1, (n_ranks >1 ? mesh->global_vtx_num[i] : (fvm_gnum_t)i+1),
               o2n_vtx_id[i]+1,  o2n_vtx_gnum[i]);
  bft_printf_flush();
#endif

  /* Update global vertex information */

  mesh->n_vertices = n_am_vertices;
  mesh->n_g_vertices = mesh->n_vertices;

  BFT_REALLOC(new_vtx_gnum, n_am_vertices, fvm_gnum_t);
  *p_old_vtx_gnum = mesh->global_vtx_num;
  mesh->global_vtx_num = new_vtx_gnum;

#if defined(HAVE_MPI)
  if (n_ranks > 1) { /* Get temporary a no contiguous numbering */

    fvm_gnum_t  glob_gmax;
    fvm_gnum_t loc_gmax = new_vtx_gnum[n_am_vertices - 1];
    MPI_Comm  mpi_comm = cs_glob_mpi_comm;

    /* Find the max. global number */

    MPI_Allreduce(&loc_gmax, &glob_gmax, 1, FVM_MPI_GNUM, MPI_MAX, mpi_comm);

    mesh->n_g_vertices = glob_gmax;

  }
#endif

  /* Update vtx_coord for initial vertices */

  BFT_MALLOC(new_vtx_coord, 3*n_am_vertices, cs_real_t);

  for (i = 0; i < n_bm_vertices; i++)  /* Initialize new vertex coord. */
    for (k = 0; k < 3; k++)
      new_vtx_coord[3*o2n_vtx_id[i]+k] = mesh->vtx_coord[3*i+k];

  /* Update vtx_coord for new vertices */

  for (i = 0; i < join_mesh->n_vertices; i++)
    for (k = 0; k < 3; k++)
      new_vtx_coord[3*join2mesh_vtx_id[i]+k] = j_vertices[i].coord[k];

  BFT_FREE(mesh->vtx_coord);
  mesh->vtx_coord = new_vtx_coord;

  /* Return pointer */

  *p_join2mesh_vtx_id = join2mesh_vtx_id;
  *p_o2n_vtx_id = o2n_vtx_id;

}

/*----------------------------------------------------------------------------
 * Update elements of a cs_mesh_t structure implied in the joining
 * operation. Update the vertices after the face split step.
 *
 * parameters:
 *  join_mesh          <--  pointer to a local cs_join_mesh_t struct.
 *  mesh               <->  pointer of pointer to cs_mesh_t structure
 *  p_join2mesh_vtx_id -->  relation between join mesh and after merge vertex
 *---------------------------------------------------------------------------*/

static void
_update_vertices_after_split(cs_join_mesh_t   *join_mesh,
                             cs_mesh_t        *mesh,
                             cs_int_t         *p_join2mesh_vtx_id[])
{
  cs_int_t  i, j, k, o_id, j_id, v_id;
  fvm_gnum_t  prev, cur;

  cs_int_t  n_as_vertices = -1; /* ac: after splitting */
  fvm_lnum_t  *order = NULL;
  cs_real_t  *new_vtx_coord = NULL;
  cs_int_t   *o2n_vtx_id = NULL, *join2mesh_vtx_id = NULL;
  fvm_gnum_t  *new_vtx_gnum = NULL, *tmp_vtx_gnum = NULL;

  const cs_int_t  n_bs_vertices = mesh->n_vertices; /* bs: before splitting */
  const cs_int_t  j_n_vertices = join_mesh->n_vertices;
  const cs_join_vertex_t  *j_vertices = join_mesh->vertices;
  const cs_int_t  n_vertices = n_bs_vertices + j_n_vertices;

  /* Update initial vertices (local and global numbering) */

  BFT_MALLOC(o2n_vtx_id, n_bs_vertices, cs_int_t);
  BFT_MALLOC(join2mesh_vtx_id, j_n_vertices, cs_int_t);
  BFT_MALLOC(tmp_vtx_gnum, n_vertices, fvm_gnum_t);
  BFT_MALLOC(new_vtx_gnum, n_vertices, fvm_gnum_t);
  BFT_MALLOC(order, n_vertices, fvm_lnum_t);

  for (i = 0; i < n_bs_vertices; i++)
    tmp_vtx_gnum[i] = mesh->global_vtx_num[i];

  for (i = 0, j = n_bs_vertices; i < j_n_vertices; i++, j++)
    tmp_vtx_gnum[j] = j_vertices[i].gnum;

  fvm_order_local_allocated(NULL, tmp_vtx_gnum, order, n_vertices);

  /* Define o2n_vtx_id and join2mesh_vtx_id arrays */

  if (order[0] < n_bs_vertices)
    prev = mesh->global_vtx_num[order[0]] + 1;
  else
    prev = j_vertices[order[0]-n_bs_vertices].gnum + 1;

  for (i = 0; i < n_vertices; i++) {

    o_id = order[i];

    if (o_id < n_bs_vertices) { /* Belongs to the initial mesh */

      cur = mesh->global_vtx_num[o_id];

      if (cur != prev) {

        n_as_vertices++;
        prev = cur;
        o2n_vtx_id[o_id] = n_as_vertices;
        new_vtx_gnum[n_as_vertices] = cur;

      }
      else
        o2n_vtx_id[o_id] = n_as_vertices;

    }
    else { /* Belongs to the join_mesh */

      j_id = o_id - n_bs_vertices;
      cur = j_vertices[j_id].gnum;

      if (cur != prev) {

        n_as_vertices++;
        prev = cur;
        new_vtx_gnum[n_as_vertices] = cur;
        join2mesh_vtx_id[j_id] = n_as_vertices;

      }
      else
        join2mesh_vtx_id[j_id] = n_as_vertices;

    }

  } /* End of loop on vertices */

  /* n_as_vertices was up to now an id. Move to a number */
  n_as_vertices++;

  /* Memory management */

  BFT_FREE(tmp_vtx_gnum);
  BFT_FREE(order);
  BFT_REALLOC(new_vtx_gnum, n_as_vertices, fvm_gnum_t);
  BFT_MALLOC(new_vtx_coord, 3*n_as_vertices, cs_real_t);

  mesh->n_vertices = n_as_vertices;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n\n Dump Old2New array (local mesh): "
             "old_n_vertices = %d - new_n_vertices = %d\n",
             n_bs_vertices, n_as_vertices);
  for (i = 0; i < n_bs_vertices; i++)
    bft_printf("Old num : %7d (%9u) => New num : %7d (%9u)\n",
               i+1, (cs_glob_n_ranks >1 ? mesh->global_vtx_num[i] : (fvm_gnum_t)i+1),
               o2n_vtx_id[i]+1,  new_vtx_gnum[o2n_vtx_id[i]]);
  bft_printf_flush();
#endif

  /* Update vtx_coord for initial vertices */

  for (i = 0; i < n_bs_vertices; i++)  /* Initialize new vertex coord. */
    for (k = 0; k < 3; k++)
      new_vtx_coord[3*o2n_vtx_id[i]+k] = mesh->vtx_coord[3*i+k];

  /* Update vtx_coord for new vertices */

  for (i = 0; i < join_mesh->n_vertices; i++)
    for (k = 0; k < 3; k++)
      new_vtx_coord[3*join2mesh_vtx_id[i]+k] = j_vertices[i].coord[k];

  /* Update interior face connect. */

  for (i = 0; i < mesh->n_i_faces; i++) {
    for (j = mesh->i_face_vtx_idx[i]-1; j < mesh->i_face_vtx_idx[i+1]-1; j++) {
      v_id = mesh->i_face_vtx_lst[j] - 1;
      mesh->i_face_vtx_lst[j] = o2n_vtx_id[v_id] + 1;
    }
  }

  /* Update border face connect. */

  for (i = 0; i < mesh->n_b_faces; i++) {
    for (j = mesh->b_face_vtx_idx[i]-1; j < mesh->b_face_vtx_idx[i+1]-1; j++) {
      v_id = mesh->b_face_vtx_lst[j] - 1;
      mesh->b_face_vtx_lst[j] = o2n_vtx_id[v_id] + 1;
    }
  }

  /* Free memory */

  BFT_FREE(o2n_vtx_id);
  BFT_FREE(mesh->vtx_coord);
  BFT_FREE(mesh->global_vtx_num);

  mesh->vtx_coord = new_vtx_coord;
  mesh->global_vtx_num = new_vtx_gnum;

  /* Return pointer */

  *p_join2mesh_vtx_id = join2mesh_vtx_id;

}

/*----------------------------------------------------------------------------
 * Allocate and initialize the definition of a edge_builder_t structure.
 * Equivalent to build a local edge-based connectivity for the selected faces.
 *
 * parameters:
 *  join_select   <--  list of all implied entities in the joining operation
 *  mesh          <--  pointer to a cs_mesh_t structure
 *
 * returns:
 *  a new allocated pointer to an edge_builder_t structure.
 *---------------------------------------------------------------------------*/

static edge_builder_t *
_init_edge_builder(cs_join_select_t    *join_select,
                   cs_mesh_t           *mesh)
{
  cs_int_t  i, j, save, shift;

  cs_int_t  *count = NULL;
  edge_builder_t  *edge_builder = NULL;

  assert(sizeof(cs_int_t) == sizeof(fvm_lnum_t));

  /* Allocate and initialize edge_builder_t structure */

  BFT_MALLOC(edge_builder, 1, edge_builder_t);

  edge_builder->n_vertices = mesh->n_vertices;
  edge_builder->n_edges = 0;
  edge_builder->v2v_lst = NULL;
  edge_builder->v2v_sub_idx = NULL;
  edge_builder->v2v_sub_lst = NULL;

  BFT_MALLOC(edge_builder->v2v_idx, edge_builder->n_vertices + 1, cs_int_t);

  for (i = 0; i < edge_builder->n_vertices + 1; i++)
    edge_builder->v2v_idx[i] = 0;

  /* Build vertex -> vertex index */

  cs_join_build_edges_idx(join_select->n_faces,
                          join_select->faces,
                          mesh->b_face_vtx_idx,
                          mesh->b_face_vtx_lst,
                          edge_builder->v2v_idx);

  if(join_select->n_b_s_faces > 0)
    cs_join_build_edges_idx(join_select->n_b_s_faces,
                            join_select->b_s_faces,
                            mesh->b_face_vtx_idx,
                            mesh->b_face_vtx_lst,
                            edge_builder->v2v_idx);

  if(join_select->n_i_s_faces > 0)
    cs_join_build_edges_idx(join_select->n_i_s_faces,
                            join_select->i_s_faces,
                            mesh->i_face_vtx_idx,
                            mesh->i_face_vtx_lst,
                            edge_builder->v2v_idx);

  BFT_MALLOC(count, edge_builder->n_vertices, cs_int_t);

  for (i = 0; i < edge_builder->n_vertices; i++) {
    edge_builder->v2v_idx[i+1] += edge_builder->v2v_idx[i];
    count[i] = 0;
  }

  edge_builder->n_edges = edge_builder->v2v_idx[edge_builder->n_vertices];

  /* Build vertex -> vertex list */

  BFT_MALLOC(edge_builder->v2v_lst, edge_builder->n_edges, cs_int_t);

  /* Fill v2v_lst */

  cs_join_build_edges_lst(join_select->n_faces,
                          join_select->faces,
                          mesh->b_face_vtx_idx,
                          mesh->b_face_vtx_lst,
                          count,
                          edge_builder->v2v_idx,
                          edge_builder->v2v_lst);

  if(join_select->n_b_s_faces > 0)
    cs_join_build_edges_lst(join_select->n_b_s_faces,
                            join_select->b_s_faces,
                            mesh->b_face_vtx_idx,
                            mesh->b_face_vtx_lst,
                            count,
                            edge_builder->v2v_idx,
                            edge_builder->v2v_lst);

  if(join_select->n_i_s_faces > 0)
    cs_join_build_edges_lst(join_select->n_i_s_faces,
                            join_select->i_s_faces,
                            mesh->i_face_vtx_idx,
                            mesh->i_face_vtx_lst,
                            count,
                            edge_builder->v2v_idx,
                            edge_builder->v2v_lst);

  /* Free memory */

  BFT_FREE(count);

  /* Ordering in order to clean the list */

  for (i = 0; i < edge_builder->n_vertices; i++)
    cs_sort_shell(edge_builder->v2v_idx[i],
                  edge_builder->v2v_idx[i+1],
                  edge_builder->v2v_lst);

  /* Delete redundancies. Clean structure. */

  save = edge_builder->v2v_idx[0];
  shift = 0;

  for (i = 0; i < edge_builder->n_vertices; i++) {

    cs_int_t  start = save;
    cs_int_t  end = edge_builder->v2v_idx[i+1];

    if (end - start > 0) {

      edge_builder->v2v_lst[shift++] = edge_builder->v2v_lst[start];

      for (j = start + 1; j < end; j++)
        if (edge_builder->v2v_lst[j-1] != edge_builder->v2v_lst[j])
          edge_builder->v2v_lst[shift++] = edge_builder->v2v_lst[j];

    }

    save = end;
    edge_builder->v2v_idx[i+1] = shift;

  }

  edge_builder->n_edges = edge_builder->v2v_idx[edge_builder->n_vertices];
  BFT_REALLOC(edge_builder->v2v_lst, edge_builder->n_edges, cs_int_t);

  return edge_builder;
}

/*----------------------------------------------------------------------------
 * Get the local face connectivity for the current selected face before/after
 * the merge step.
 *
 * parameters:
 *  select_id      <--  id of the selected face in selection array
 *  o2n_vtx_gnum   <--  local array on the new global vertex
 *  join_select    <--  list of all implied entities in the joining operation
 *  join_mesh      <--  pointer to the local cs_join_mesh_t structure
 *  mesh           <--  pointer to cs_mesh_t structure
 *  bm_tmp         <--  connectivity before the merge step
 *  am_tmp         <--  connectivity after the merge step
 *---------------------------------------------------------------------------*/

static void
_get_local_faces_connect(cs_int_t            select_id,
                         fvm_gnum_t          o2n_vtx_gnum[],
                         cs_join_select_t   *join_select,
                         cs_join_mesh_t     *join_mesh,
                         cs_mesh_t          *mesh,
                         cs_int_t            bm_tmp[],
                         cs_int_t            am_tmp[])
{
  cs_int_t  i, j, k, v_id, bm_shift;
  fvm_gnum_t  new_gnum, v_gnum;

  cs_int_t  fid = join_select->faces[select_id] - 1;
  fvm_gnum_t  fgnum = join_select->compact_face_gnum[select_id];
  cs_int_t  am_s = join_mesh->face_vtx_idx[select_id] - 1;
  cs_int_t  am_e = join_mesh->face_vtx_idx[select_id+1] - 1;
  cs_int_t  n_am_face_vertices = am_e - am_s;
  cs_int_t  bm_s = mesh->b_face_vtx_idx[fid] - 1;
  cs_int_t  bm_e = mesh->b_face_vtx_idx[fid+1] - 1;
  cs_int_t  n_bm_face_vertices = bm_e - bm_s;
  cs_int_t  fst_match_id = -1;

  const cs_join_vertex_t  *vertices = join_mesh->vertices;
  const int  n_ranks = cs_glob_n_ranks;

  assert(join_mesh->face_gnum[select_id] == fgnum);

  /* Store the face connectivity before the fusion step */

  for (i = bm_s, j = 0; i < bm_e; i++, j++)
    bm_tmp[j] = mesh->b_face_vtx_lst[i] - 1;
  bm_tmp[n_bm_face_vertices] = mesh->b_face_vtx_lst[bm_s] - 1;

  /* Store the face connectivity after the fusion step */

  for (i = am_s, j = 0; i < am_e; i++, j++)
    am_tmp[j] = join_mesh->face_vtx_lst[i] - 1;
  am_tmp[n_am_face_vertices] = join_mesh->face_vtx_lst[am_s] - 1;

  /* Find position of initial vertices in the face connectivity after the
     fusion step. If not found -1 (initialized value) */

  bm_shift = 0;
  while (fst_match_id == -1 && bm_shift < n_bm_face_vertices) {

    v_id = bm_tmp[bm_shift];
    v_gnum = (n_ranks > 1 ? mesh->global_vtx_num[v_id] : (fvm_gnum_t)v_id+1);
    new_gnum = o2n_vtx_gnum[v_id];

    for (k = 0; k < n_am_face_vertices; k++) {
      if (vertices[am_tmp[k]].gnum == new_gnum) {
        fst_match_id = k;
        break;
      }
    }

    if (fst_match_id == -1)
      bm_shift++;

  }

  if (fst_match_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              _("  Cannot find the first corresponding vertex between the face"
                " connectivity before/after the fusion step.\n"
                "  Current global face number: %u\n"), fgnum);

  /* Store the face connectivity before the fusion step */

  for (i = 0; i < n_bm_face_vertices; i++) {
    j = bm_s + (bm_shift + i) % n_bm_face_vertices;
    bm_tmp[i] = mesh->b_face_vtx_lst[j] - 1;
  }
  bm_tmp[n_bm_face_vertices] = mesh->b_face_vtx_lst[bm_s + bm_shift] - 1;

  /* Store the face connectivity after the fusion step */

  for (i = 0; i < n_am_face_vertices; i++) {
    j = am_s + (fst_match_id + i) % n_am_face_vertices;
    am_tmp[i] = join_mesh->face_vtx_lst[j] - 1;
  }
  am_tmp[n_am_face_vertices] = join_mesh->face_vtx_lst[am_s + fst_match_id] -1;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n\n  Face: %d (%u)\n", fid+1, fgnum);
  bft_printf("bm_tmp:%d, n: %2d, v:", i, n_bm_face_vertices+1);
  for (j = 0; j < n_bm_face_vertices + 1; j++)
    bft_printf(" %8d (%u)", bm_tmp[j]+1, o2n_vtx_gnum[bm_tmp[j]]);
  bft_printf("\nam_tmp:%d, n: %2d, v:", i, n_am_face_vertices+1);
  for (j = 0; j < n_am_face_vertices + 1; j++)
    bft_printf(" %8d (%u)", am_tmp[j]+1, vertices[am_tmp[j]].gnum);
  bft_printf("\n");
  bft_printf_flush();
#endif

}

/*----------------------------------------------------------------------------
 * Define the new edge connectivity for the initial edges.
 *
 * parameters:
 *  join_select       <--  keep all implied entities in the joining operation
 *  join_mesh         <--  pointer to the local cs_join_mesh_t structure
 *  mesh              <--  pointer of pointer to cs_mesh_t structure
 *  o2n_vtx_gnum      <--  local array on the new global vertex
 *  join2mesh_vtx_id  <--  relation between join mesh and after fusion vertex
 *  edge_builder      <->  pointer to an edge_builder_t structure
 *---------------------------------------------------------------------------*/

static void
_complete_edge_builder(cs_join_select_t    *join_select,
                       cs_join_mesh_t      *join_mesh,
                       cs_mesh_t           *mesh,
                       fvm_gnum_t           o2n_vtx_gnum[],
                       cs_int_t             join2mesh_vtx_id[],
                       edge_builder_t      *edge_builder)
{
  cs_int_t  i, j, j1, j2, k, shift;
  cs_int_t  v1_id, v2_id, edge_id, n_subs;
  fvm_gnum_t  v1_gnum, v2_gnum;
  cs_bool_t  direct_scan, degenerated_edge;

  cs_int_t  am_max = 0, bm_max = 0;
  cs_int_t  *am_tmp = NULL, *bm_tmp = NULL;

  const cs_join_vertex_t  *vertices = join_mesh->vertices;

  /* Sanity checks */

  assert(edge_builder != NULL);
  assert(join_mesh->n_faces == join_select->n_faces);

  /* Define a list of new vertices for each initial selected edge */

  BFT_MALLOC(edge_builder->v2v_sub_idx, edge_builder->n_edges + 1, cs_int_t);

  for (i = 0; i < edge_builder->n_edges + 1; i++)
    edge_builder->v2v_sub_idx[i] = 0;

  for (i = 0; i < join_mesh->n_faces; i++)
    am_max = CS_MAX(am_max,
                    join_mesh->face_vtx_idx[i+1]-join_mesh->face_vtx_idx[i]);

  for (i = 0; i < join_select->n_faces; i++) {
    j = join_select->faces[i] - 1;
    bm_max = CS_MAX(bm_max, mesh->b_face_vtx_idx[j+1]-mesh->b_face_vtx_idx[j]);
  }

  BFT_MALLOC(am_tmp, am_max + 1, cs_int_t);
  BFT_MALLOC(bm_tmp, bm_max + 1, cs_int_t);

  /* Count the number of sub-elements to add to each initial edge */

  for (i = 0; i < join_select->n_faces; i++) {

    cs_int_t  fid = join_select->faces[i] - 1;
    fvm_gnum_t  fgnum = join_select->compact_face_gnum[i];
    cs_int_t  bm_s = mesh->b_face_vtx_idx[fid] - 1;
    cs_int_t  bm_e = mesh->b_face_vtx_idx[fid+1] - 1;
    cs_int_t  n_bm_face_vertices = bm_e - bm_s;
    cs_int_t  am_s = join_mesh->face_vtx_idx[i] - 1;
    cs_int_t  am_e = join_mesh->face_vtx_idx[i+1] - 1;
    cs_int_t  n_am_face_vertices = am_e - am_s;

    _get_local_faces_connect(i,           /* id of the face */
                             o2n_vtx_gnum,
                             join_select,
                             join_mesh,
                             mesh,
                             bm_tmp,
                             am_tmp);

    for (j = 0, j1 = 0, j2 = 0; j < n_bm_face_vertices; j++) {

      degenerated_edge = false;
      v1_id = bm_tmp[j];
      v2_id = bm_tmp[j+1];
      v1_gnum = o2n_vtx_gnum[v1_id];
      v2_gnum = o2n_vtx_gnum[v2_id];

      if (v1_gnum != vertices[am_tmp[j1]].gnum)
        j2 = j1 + 1;

      else {

        edge_id = _get_join_edge_id(v1_id, v2_id, edge_builder);

        assert(edge_id != -1);

        if (v1_gnum == v2_gnum) { /* Initial edge has been deleted */
          n_subs = 1;
          degenerated_edge = true;
        }
        else { /* Look for the next initial vertex */

          for (j2 = j1 + 1; j2 < n_am_face_vertices + 1; j2++)
            if (v2_gnum == vertices[am_tmp[j2]].gnum)
              break;

          if (j2 == n_am_face_vertices + 1) { /* Init. edge has been deleted */

            j2 = j1 + 1;
            degenerated_edge = true;
            n_subs = 1;

          }
          else
            n_subs = j2 - j1 - 1;

        }

        /* Add n_subs elements to the initial edge definition */

        if (edge_builder->v2v_sub_idx[edge_id+1] == 0)
          edge_builder->v2v_sub_idx[edge_id+1] = n_subs;

        else { /* This edge has already been scanned */

          if (edge_builder->v2v_sub_idx[edge_id+1] != n_subs) {

            /* Different result for the same edge is only possible for
               a degenerated edge */

            if (degenerated_edge == false && n_subs == 1)
              bft_error(__FILE__, __LINE__, 0,
                        _("  Face %d (%u): Edge with two different"
                          " descriptions: (%d, %d) [%u, %u]\n"
                          "  n_subs: %d - previous n_subs: %d\n"
                          "  Impossible to continue the mesh update after the"
                          " fusion operation.\n"),
                        fid+1, fgnum, v1_id+1, v2_id+1, v1_gnum, v2_gnum,
                        n_subs, edge_builder->v2v_sub_idx[edge_id+1]);
            else
              edge_builder->v2v_sub_idx[edge_id+1] =
                CS_MAX(n_subs, edge_builder->v2v_sub_idx[edge_id+1]);

          }

        }

        j1 = j2;

      } /* v1_gnum == vertices[am_tmp[j1]].gnum */

    } /* End of loop on initial face connectivity */

  } /* End of loop on selected faces */

  /* Build sub index */

  for (i = 0; i < edge_builder->n_edges; i++)
    edge_builder->v2v_sub_idx[i+1] += edge_builder->v2v_sub_idx[i];

  BFT_MALLOC(edge_builder->v2v_sub_lst,
             edge_builder->v2v_sub_idx[edge_builder->n_edges], cs_int_t);

  for (i = 0; i < edge_builder->v2v_sub_idx[edge_builder->n_edges]; i++)
    edge_builder->v2v_sub_lst[i] = -1; /* value = degenerated edge */

  /* Add sub-elements to each initial edge */

  for (i = 0; i < join_select->n_faces; i++) {

    cs_int_t  fid = join_select->faces[i] - 1;
    cs_int_t  bm_s = mesh->b_face_vtx_idx[fid] - 1;
    cs_int_t  bm_e = mesh->b_face_vtx_idx[fid+1] - 1;
    cs_int_t  n_bm_face_vertices = bm_e - bm_s;
    cs_int_t  am_s = join_mesh->face_vtx_idx[i] - 1;
    cs_int_t  am_e = join_mesh->face_vtx_idx[i+1] - 1;
    cs_int_t  n_am_face_vertices = am_e - am_s;

    _get_local_faces_connect(i,           /* id of the face */
                             o2n_vtx_gnum,
                             join_select,
                             join_mesh,
                             mesh,
                             bm_tmp,
                             am_tmp);

    for (j = 0, j1 = 0, j2 = 0; j < n_bm_face_vertices; j++) {

      v1_id = bm_tmp[j];
      v2_id = bm_tmp[j+1];
      v1_gnum = o2n_vtx_gnum[v1_id];
      v2_gnum = o2n_vtx_gnum[v2_id];

      if (v1_gnum != vertices[am_tmp[j1]].gnum)
        j2 = j1 + 1;

      else {

        if (v1_id < v2_id)
          direct_scan = true;
        else
          direct_scan = false;

        edge_id = _get_join_edge_id(v1_id, v2_id, edge_builder);

        if (v1_gnum != v2_gnum) { /* Initial edge is not degenerated */

          for (j2 = j1 + 1; j2 < n_am_face_vertices + 1; j2++)
            if (v2_gnum == vertices[am_tmp[j2]].gnum)
              break;

          if (j2 == n_am_face_vertices + 1) /* Initial edge is degenerated */
            j2 = j1 + 1;

          else {  /* Add n_subs elements to the initial edge definition */

            n_subs = j2 - j1 - 1;
            shift = edge_builder->v2v_sub_idx[edge_id];

            if (direct_scan == true) {
              for (k = 0; k < n_subs; k++)
                edge_builder->v2v_sub_lst[shift + k] =
                  join2mesh_vtx_id[am_tmp[j1 + 1 + k]] + 1;
            }
            else {
              for (k = 0; k < n_subs; k++)
                edge_builder->v2v_sub_lst[shift + k] =
                  join2mesh_vtx_id[am_tmp[j1 + n_subs - k]] + 1;
              /* j1 + 1 + n_subs - 1 - k */
            }

          }

        } /* Initial edge is not degenerated */

        j1 = j2;

      } /* v1_gnum == vertices[am_tmp[j]].gnum */

    } /* End of loop on initial face connectivity */

  } /* End of loop on selected faces */

  /* Free memory */

  BFT_FREE(am_tmp);
  BFT_FREE(bm_tmp);

}

/*----------------------------------------------------------------------------
 * Update selected face connectivity after the merge step.
 *
 * parameters:
 *  n_adj_faces     <--  number of adjacent faces
 *  adj_faces       <--  list of adjacent face numbers
 *  n_s_faces       <--  number of single faces
 *  s_faces         <--  list of "single" face numbers
 *  p_n_sel_faces   <--  pointer to the number of faces in the selection list
 *  p_sel_faces     <--  pointer to the selection list
 *---------------------------------------------------------------------------*/

static void
_get_select_face_lst(cs_int_t     n_adj_faces,
                     cs_int_t     adj_faces[],
                     cs_int_t     n_s_faces,
                     cs_int_t     s_faces[],
                     cs_int_t    *p_n_sel_faces,
                     cs_int_t    *p_sel_faces[])
{
  cs_int_t  i, shift, prev, cur, n;

  cs_int_t  *order = NULL, *tmp = NULL;

  n = n_adj_faces + n_s_faces;

  BFT_MALLOC(tmp, n, cs_int_t);

  shift = 0;
  for (i = 0; i < n_adj_faces; i++)
    tmp[shift++] = adj_faces[i];
  for (i = 0; i < n_s_faces; i++)
    tmp[shift++] = s_faces[i];

  assert(shift == n);

  /* Order selection and eliminate redundancies */

  BFT_MALLOC(order, n, cs_int_t);

  fvm_order_local_allocated(tmp, NULL, order, n);

  prev = 0;
  shift = 0;

  for (i = 0; i < n; i++) {

    cur = tmp[order[i]];
    if (cur != prev) {
      prev = cur;
      order[shift++] = cur;
    }

  }

  BFT_FREE(tmp);
  BFT_REALLOC(order, shift, cs_int_t);

  /* Return pointers */

  *p_n_sel_faces = shift;
  *p_sel_faces = order;

}

/*----------------------------------------------------------------------------
 * Update selected face connectivity after the fusion step.
 *
 * parameters:
 *  join_select       <--  cs_join_select_t struct.
 *  join_mesh         <--  cs_join_mesh_t structure
 *  join2mesh_vtx_id  <--  relation between vtx_id of join_mesh & cs_mesh_t
 *  n_faces           <--  local number of faces in the mesh
 *  p_f2v_idx         <--  face -> vertex connectivity index
 *  p_f2v_lst         <--  face -> vertex connectivity list
 *---------------------------------------------------------------------------*/

static void
_update_selected_face_connect(cs_join_select_t   *join_select,
                              cs_join_mesh_t     *join_mesh,
                              cs_int_t            join2mesh_vtx_id[],
                              cs_int_t            n_faces,
                              cs_int_t           *p_f2v_idx[],
                              cs_int_t           *p_f2v_lst[])
{
  cs_int_t  i, j, shift, v_id, select_id, join_fid;
  fvm_gnum_t  fgnum;

  cs_int_t  *new_f2v_lst = NULL, *new_f2v_idx = NULL;
  cs_int_t  *f2v_idx = *p_f2v_idx;
  cs_int_t  *f2v_lst = *p_f2v_lst;

  BFT_MALLOC(new_f2v_idx, n_faces + 1, cs_int_t);

  for (i = 0; i < n_faces + 1; i++)
    new_f2v_idx[i] = 0;

  for (i = 0, select_id = 0; i < n_faces; i++) {

    cs_bool_t  in_selection = false;

    if (select_id < join_select->n_faces)
      if (i+1 == join_select->faces[select_id])
        in_selection = true;

    if (in_selection) { /* Among selected faces */

      fgnum = join_select->compact_face_gnum[select_id];

      if (join_mesh->face_gnum[select_id] == fgnum)
        join_fid = select_id;
      else
        join_fid = cs_search_g_binary(0,
                                      join_mesh->n_faces-1,
                                      fgnum,
                                      join_mesh->face_gnum);

      assert(join_fid < join_mesh->n_faces);
      assert(join_fid != -1);

      new_f2v_idx[i+1] =  join_mesh->face_vtx_idx[join_fid+1]
                        - join_mesh->face_vtx_idx[join_fid];

      select_id++;

    }
    else /* Not a selected face. Do not update. */
      new_f2v_idx[i+1] = f2v_idx[i+1] - f2v_idx[i];

  } /* End of loop on border faces */

  new_f2v_idx[0] = 1;
  for (i = 0; i < n_faces; i++)
    new_f2v_idx[i+1] += new_f2v_idx[i];

  BFT_MALLOC(new_f2v_lst, new_f2v_idx[n_faces]-1, cs_int_t);

  for (i = 0, select_id = 0; i < n_faces; i++) {

    cs_bool_t  in_selection = false;

    if (select_id < join_select->n_faces)
      if (i+1 == join_select->faces[select_id])
        in_selection = true;

    shift = new_f2v_idx[i]-1;

    if (in_selection) { /* Among selected faces */

      fgnum = join_select->compact_face_gnum[select_id];

      if (join_mesh->face_gnum[select_id] == fgnum)
        join_fid = select_id;
      else
        join_fid = cs_search_g_binary(0,
                                      join_mesh->n_faces-1,
                                      fgnum,
                                      join_mesh->face_gnum);

      for (j = join_mesh->face_vtx_idx[join_fid]-1;
           j < join_mesh->face_vtx_idx[join_fid+1]-1; j++, shift++) {
        v_id = join_mesh->face_vtx_lst[j] - 1;
        new_f2v_lst[shift] = join2mesh_vtx_id[v_id] + 1;
      }

      select_id++;

    }
    else { /* Not a selected face. Do not update. */

      for (j = f2v_idx[i]-1; j < f2v_idx[i+1]-1; j++, shift++)
        new_f2v_lst[shift] = f2v_lst[j];

    }

  } /* End of loop on border faces */

  BFT_FREE(f2v_lst);
  BFT_FREE(f2v_idx);

  /* Return pointers */

  *p_f2v_idx = new_f2v_idx;
  *p_f2v_lst = new_f2v_lst;

}

/*----------------------------------------------------------------------------
 * Update adjacent face connectivity.
 *
 * parameters:
 *  n_adj_faces    <--  number of adjacent faces
 *  adj_faces      <--  list of adjacent face numbers
 *  edge_builder   <--  edge_builder_t structure
 *  o2n_vtx_id     <--  relation between old/new vertix id
 *  n_faces        <--  local number of faces in the mesh
 *  p_f2v_idx      <--  face -> vertex connectivity index
 *  p_f2v_lst      <--  face -> vertex connectivity list
 *---------------------------------------------------------------------------*/

static void
_update_adj_face_connect(cs_int_t          n_adj_faces,
                         cs_int_t          adj_faces[],
                         edge_builder_t   *edge_builder,
                         cs_int_t          o2n_vtx_id[],
                         cs_int_t          n_faces,
                         cs_int_t         *p_f2v_idx[],
                         cs_int_t         *p_f2v_lst[])
{
  cs_int_t  i, j, k, l, v1_id, v2_id, n_face_vertices, shift, select_id;
  cs_int_t  s, e, v_s, v_e, v_sub_s, v_sub_e, edge_id;

  cs_int_t  max = 0;
  cs_int_t  *new_f2v_lst = NULL, *new_f2v_idx = NULL, *tmp = NULL;
  cs_int_t  *f2v_idx = *p_f2v_idx;
  cs_int_t  *f2v_lst = *p_f2v_lst;

  BFT_MALLOC(new_f2v_idx, n_faces+1, cs_int_t);

  for (i = 0; i < n_faces + 1; i++)
    new_f2v_idx[i] = 0;

  for (i = 0; i < n_faces; i++)
    max = CS_MAX(max, f2v_idx[i+1] - f2v_idx[i]);

  BFT_MALLOC(tmp, max + 1, cs_int_t);

  /* first: update index (counting phase) */

  for (i = 0, select_id = 0; i < n_faces; i++) {

    cs_bool_t  in_selection = false;

    if (select_id < n_adj_faces)
      if (i+1 == adj_faces[select_id])
        in_selection = true;

    if (in_selection) { /* Among selected faces */

      cs_int_t  count = 0;

      s = f2v_idx[i] - 1;
      e = f2v_idx[i+1] - 1;
      n_face_vertices = e - s;

      for (k = 0; k < n_face_vertices; k++)
        tmp[k] = f2v_lst[s + k] - 1;
      tmp[n_face_vertices] = f2v_lst[s] - 1;

      for (j = 0; j < n_face_vertices; j++) { /* Scan edges */

        if (tmp[j] < tmp[j+1])
          v1_id = tmp[j], v2_id = tmp[j+1];
        else if (tmp[j+1] < tmp[j])
          v1_id = tmp[j+1], v2_id = tmp[j];
        else /* delete the current edge (count += 0) */
          v1_id = -1;

        if (v1_id > -1) {

          v_s = edge_builder->v2v_idx[v1_id];
          v_e = edge_builder->v2v_idx[v1_id+1];
          count += 1; /* for v1_id */

          for (edge_id = v_s; edge_id < v_e; edge_id++)
            if (edge_builder->v2v_lst[edge_id] == v2_id + 1)
              break;

          if (edge_id < v_e) {

            v_sub_s = edge_builder->v2v_sub_idx[edge_id];
            v_sub_e = edge_builder->v2v_sub_idx[edge_id + 1];

            if (v_sub_e - v_sub_s > 0)
              if (edge_builder->v2v_sub_lst[v_sub_s] != -1)
                count +=  v_sub_e - v_sub_s;

          } /* End if exist edge_id */

        } /* End if v1_id > -1 */

      } /* End of loop on face connectivity */

      assert(count > 2); /* At least a triangle */
      assert(select_id < n_adj_faces);

      new_f2v_idx[i+1] = count;
      select_id++;

    }
    else  /* Not an adjacent face. Do not update. */
      new_f2v_idx[i+1] = f2v_idx[i+1] - f2v_idx[i];

  } /* End of loop on faces */

  new_f2v_idx[0] = 1;
  for (i = 0; i < n_faces; i++)
    new_f2v_idx[i+1] += new_f2v_idx[i];

  BFT_MALLOC(new_f2v_lst, new_f2v_idx[n_faces]-1, cs_int_t);

  /* second: update list (filling phase) */

  for (i = 0, select_id = 0; i < n_faces; i++) {

    cs_bool_t  direct_scan;
    cs_bool_t  in_selection = false;

    if (select_id < n_adj_faces)
      if (i+1 == adj_faces[select_id])
        in_selection = true;

    shift = new_f2v_idx[i] - 1;

    if (in_selection) { /* Among selected faces */

      cs_int_t  count = 0;

      s = f2v_idx[i] - 1;
      e = f2v_idx[i+1] - 1;
      n_face_vertices = e - s;

      for (k = 0; k < n_face_vertices; k++)
        tmp[k] = f2v_lst[s + k] - 1;
      tmp[n_face_vertices] = f2v_lst[s] - 1;

      for (j = 0; j < n_face_vertices; j++) { /* Scan edges */

        if (tmp[j] < tmp[j+1]) {
          v1_id = tmp[j];
          v2_id = tmp[j+1];
          direct_scan = true;
        }
        else if (tmp[j+1] < tmp[j]) {
          v1_id = tmp[j+1];
          v2_id = tmp[j];
          direct_scan = false;
        }
        else /*  delete the current edge (count += 0) */
          v1_id = -1;

        if (v1_id > -1) {

          v_s = edge_builder->v2v_idx[v1_id];
          v_e = edge_builder->v2v_idx[v1_id+1];

          new_f2v_lst[shift + count] = o2n_vtx_id[tmp[j]] + 1;
          count += 1; /* for v1_id */

          for (edge_id = v_s; edge_id < v_e; edge_id++)
            if (edge_builder->v2v_lst[edge_id] == v2_id + 1)
              break;

          if (edge_id < v_e) {

            v_sub_s = edge_builder->v2v_sub_idx[edge_id];
            v_sub_e = edge_builder->v2v_sub_idx[edge_id + 1];

            if (v_sub_e - v_sub_s > 0) {

              if (edge_builder->v2v_sub_lst[v_sub_s] != -1) {

                /* Edge is in edge_builder and is not to delete */

                if (direct_scan == true)
                  for (l = v_sub_s; l < v_sub_e; l++, count++)
                    new_f2v_lst[shift + count] = edge_builder->v2v_sub_lst[l];
                else
                  for (l = v_sub_e - 1; l >= v_sub_s; l--, count++)
                    new_f2v_lst[shift + count] = edge_builder->v2v_sub_lst[l];

              } /* edge is not degenerated */

            } /* There are sub-elements to add */

          } /* End if exist edge_id */

        } /* End if v1_id > -1 */

      } /* End of loop on face connectivity */

      select_id++;

    }
    else  /* Not in selection. Do not update */
      for (j = f2v_idx[i] - 1; j < f2v_idx[i+1] - 1; j++)
        new_f2v_lst[shift++] = f2v_lst[j];

  } /* End of loop on faces */

  BFT_FREE(f2v_lst);
  BFT_FREE(f2v_idx);
  BFT_FREE(tmp);

  /* Return pointers */

  *p_f2v_idx = new_f2v_idx;
  *p_f2v_lst = new_f2v_lst;

}

/*----------------------------------------------------------------------------
 * Compute barycenter and face normal for the current face.
 *
 * parameters:
 *  n_face_vertices    <--  number of vertices defining the face
 *  face_vtx_coord     <--  coordinates of each vertex of the face
 *                         (size: n_face_vertices+1)
 *  face_barycenter    -->  barycentre of the face
 *  face_normal        -->  normal of the face (norm = face area)
 *----------------------------------------------------------------------------*/

static void
_get_face_quantity(cs_int_t     n_face_vertices,
                   cs_real_t    face_vtx_coord[],
                   cs_real_t    face_barycenter[],
                   cs_real_t    face_normal[])
{
  cs_int_t  i, coord;
  cs_real_t  v1[3], v2[3], tri_normal[3];

  cs_real_t  inv_n_face_vertices = 1/(double)n_face_vertices;

  /* Initialization */

  for (coord = 0; coord < 3; coord++) {
    face_normal[coord] = 0.0;
    face_barycenter[coord] = 0.0;
  }

  /* Compute face barycenter */

  for (i = 0; i < n_face_vertices; i++)
    for (coord = 0; coord < 3; coord++)
      face_barycenter[coord] += face_vtx_coord[3*i+coord];

  for (coord = 0; coord < 3; coord++)
    face_barycenter[coord] *= inv_n_face_vertices;

  /* Compute triangle normal and update face normal */

  for (i = 0; i < n_face_vertices; i++) {

    for (coord = 0; coord < 3; coord++) {
      v1[coord] = face_vtx_coord[3*i     + coord] - face_barycenter[coord];
      v2[coord] = face_vtx_coord[3*(i+1) + coord] - face_barycenter[coord];
    }

    _cross_product(v1, v2, tri_normal);

    for (coord = 0; coord < 3; coord++) {
      tri_normal[coord] *= 0.5;
      face_normal[coord] += tri_normal[coord];
    }

  }

  _normalize(face_normal);
}

/*----------------------------------------------------------------------------
 * Define send_rank_index and send_faces to prepare the exchange of new faces
 * between mesh structures.
 *
 * parameters:
 *  n2o_hist          <--   new -> old global face numbering
 *  gnum_rank_index   <--  index on ranks for the old global face numbering
 *  p_send_rank_index -->  index on ranks for sending face
 *  p_send_faces      -->  list of face ids to send
 *---------------------------------------------------------------------------*/

static void
_get_faces_to_send(cs_join_gset_t   *n2o_hist,
                   fvm_gnum_t        gnum_rank_index[],
                   cs_int_t         *p_send_rank_index[],
                   fvm_gnum_t       *p_send_faces[])
{
  cs_int_t  i, j, rank, shift;

  cs_int_t  reduce_size = 0;
  cs_int_t  *reduce_ids = NULL, *send_rank_index = NULL;
  fvm_gnum_t  *reduce_index = NULL, *send_faces = NULL;
  cs_join_gset_t  *distrib = NULL;

  const int  n_ranks = cs_glob_n_ranks;

  /* Sanity checks */

  assert(gnum_rank_index != NULL);
  assert(n2o_hist != NULL);
  assert(n_ranks > 1);

  distrib = cs_join_gset_create(n_ranks);

  for (i = 0; i < n_ranks; i++)
    distrib->g_elts[i] = 0; /* used to store local count */

  /* Compact init. global face distribution. Remove ranks without face
     at the begining */

  for (i = 0; i < n_ranks; i++)
    if (gnum_rank_index[i] < gnum_rank_index[i+1])
      reduce_size++;

  BFT_MALLOC(reduce_index, reduce_size+1, fvm_gnum_t);
  BFT_MALLOC(reduce_ids, reduce_size, cs_int_t);

  reduce_size = 0;
  reduce_index[0] = gnum_rank_index[0] + 1;

  for (i = 0; i < n_ranks; i++) {

    /* Add +1 to gnum_rank_index because it's an id and we work on numbers */

    if (gnum_rank_index[i] < gnum_rank_index[i+1]) {
      reduce_index[reduce_size+1] = gnum_rank_index[i+1] + 1;
      reduce_ids[reduce_size++] = i;
    }

  }

  /* Count number of ranks associated to each new face */

  for (i = 0; i < n2o_hist->n_elts; i++) {

    for (j = n2o_hist->index[i]; j < n2o_hist->index[i+1]; j++) {

      int  reduce_rank = cs_search_gindex_binary(0,
                                                 reduce_size,
                                                 n2o_hist->g_list[j],
                                                 reduce_index);

      assert(reduce_rank != -1);
      assert(reduce_rank < reduce_size);

      rank = reduce_ids[reduce_rank];
      distrib->index[rank+1] += 1;

    } /* End of loop on old faces */

  } /* End of loop on new faces */

  for (i = 0; i < n_ranks; i++)
    distrib->index[i+1] += distrib->index[i];

  BFT_MALLOC(distrib->g_list, distrib->index[n_ranks], fvm_gnum_t);

  /* Fill the list of ranks */

  for (i = 0; i < n2o_hist->n_elts; i++) {

    for (j = n2o_hist->index[i]; j < n2o_hist->index[i+1]; j++) {

      int  reduce_rank = cs_search_gindex_binary(0,
                                                 reduce_size,
                                                 n2o_hist->g_list[j],
                                                 reduce_index);

      assert(reduce_rank != -1);
      assert(reduce_rank < reduce_size);

      rank = reduce_ids[reduce_rank];
      shift = distrib->index[rank] + distrib->g_elts[rank];
      distrib->g_list[shift] = n2o_hist->g_list[j];
      distrib->g_elts[rank] += 1;

    } /* End of loop on old faces */

  } /* End of loop on new faces */

  /* Free memory */

  BFT_FREE(reduce_ids);
  BFT_FREE(reduce_index);

  /* Define arrays to return */

  BFT_MALLOC(send_rank_index, n_ranks + 1, cs_int_t);

  for (i = 0; i < n_ranks + 1; i++)
    send_rank_index[i] = distrib->index[i];

  BFT_MALLOC(send_faces, send_rank_index[n_ranks], fvm_gnum_t);

  for (i = 0; i < send_rank_index[n_ranks]; i++)
    send_faces[i] = distrib->g_list[i];

  cs_join_gset_destroy(&distrib);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n Exchange to update mesh after the face split operation:\n");
  for (i = 0; i < n_ranks; i++) {
    start = send_rank_index[i];
    end = send_rank_index[i+1];
    bft_printf(" Send to rank %5d (n = %10d):", i, end - start);
    for (j = start; j < end; j++)
      bft_printf(" %d ", send_faces[j]);
    bft_printf("\n");
  }
#endif

  /* Returns pointers */

  *p_send_rank_index = send_rank_index;
  *p_send_faces = send_faces;

}

#if defined(HAVE_MPI)
/*----------------------------------------------------------------------------
 * Get the related global cell number for each old global face number
 *
 * parameters:
 *   n_ranks        <--  number of ranks in the mpi_comm
 *   send_shift     <--  index on ranks for the face distribution
 *   send_faces     <--  list of face ids to send
 *   join_select    <--  list all local entities implied in the joining op.
 *   mpi_comm       <--  mpi communicator on which take places comm.
 *---------------------------------------------------------------------------*/

static void
_exchange_faces(cs_int_t            n_ranks,
                cs_int_t            send_shift[],
                fvm_gnum_t          send_gbuf[],
                cs_join_select_t   *join_select,
                MPI_Comm            mpi_comm)
{
  int  i, rank, fid;
  fvm_gnum_t  compact_fgnum;

  cs_int_t  *send_count = NULL, *recv_count = NULL, *recv_shift = NULL;
  fvm_gnum_t  *recv_gbuf = NULL;

  const int  loc_rank = cs_glob_rank_id;
  const fvm_gnum_t  loc_rank_s = join_select->compact_rank_index[loc_rank];
  const fvm_gnum_t  loc_rank_e = join_select->compact_rank_index[loc_rank+1];

  /* Sanity checks */

#if defined(DEBUG) && !defined(NDEBUG)
  int  n_verif_ranks;

  MPI_Comm_size(mpi_comm, &n_verif_ranks);

  assert(n_ranks == n_verif_ranks);
#endif

  assert(send_shift != NULL);
  assert(n_ranks > 1);

  /* Count the number of faces to recv */

  BFT_MALLOC(send_count, n_ranks, cs_int_t);
  BFT_MALLOC(recv_count, n_ranks, cs_int_t);

  for (i = 0; i < n_ranks; i++)
    send_count[i] = send_shift[i+1] - send_shift[i];

  /* Exchange number of elements to send */

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);

  BFT_MALLOC(recv_shift, n_ranks + 1, cs_int_t);

  /* Build index arrays */

  recv_shift[0] = 0;
  for (rank = 0; rank < n_ranks; rank++)
    recv_shift[rank+1] = recv_shift[rank] + recv_count[rank];

  BFT_MALLOC(recv_gbuf, recv_shift[n_ranks], fvm_gnum_t);

  MPI_Alltoallv(send_gbuf, send_count, send_shift, FVM_MPI_GNUM,
                recv_gbuf, recv_count, recv_shift, FVM_MPI_GNUM, mpi_comm);

  /* Get the related global cell number for each received face */

  for (rank = 0; rank < n_ranks; rank++) {

    for (i = recv_shift[rank]; i < recv_shift[rank+1]; i++) {

      compact_fgnum = recv_gbuf[i];
      fid = compact_fgnum - 1 - loc_rank_s;

      assert(loc_rank_s < compact_fgnum);
      assert(compact_fgnum <= loc_rank_e);

      recv_gbuf[i] = join_select->cell_gnum[fid];

    }

  } /* End of loop on ranks */

  /* Return values to send ranks */

  MPI_Alltoallv(recv_gbuf, recv_count, recv_shift, FVM_MPI_GNUM,
                send_gbuf, send_count, send_shift, FVM_MPI_GNUM, mpi_comm);

  /* Free memory */

  BFT_FREE(send_count);
  BFT_FREE(recv_count);
  BFT_FREE(recv_shift);
  BFT_FREE(recv_gbuf);

}
#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Get the related global cell numbers connected to the old face numbers.
 *
 * parameters:
 *  join_select      <--  list of all implied entities in the joining op.
 *  n2o_face_hist    <--  face history structure (new -> old)
 *  cell_gnum        -->  pointer to the created array
 *---------------------------------------------------------------------------*/

static void
_get_linked_cell_gnum(cs_join_select_t    *join_select,
                      cs_join_gset_t      *n2o_face_hist,
                      fvm_gnum_t          *p_cell_gnum[])
{
  cs_int_t  i, j;
  fvm_gnum_t  compact_fgnum;

  fvm_gnum_t  *cell_gnum = NULL;

  const int  n_ranks = cs_glob_n_ranks;

  BFT_MALLOC(cell_gnum,
             n2o_face_hist->index[n2o_face_hist->n_elts], fvm_gnum_t);

  if (n_ranks == 1) {

    for (i = 0; i < n2o_face_hist->n_elts; i++) {

      for (j = n2o_face_hist->index[i]; j < n2o_face_hist->index[i+1]; j++) {
        compact_fgnum = n2o_face_hist->g_list[j];
        cell_gnum[j] = join_select->cell_gnum[compact_fgnum - 1];
      }

    } /* End of loop on n2o_face_hist elements */

  }
  else {

    cs_int_t  *send_shift = NULL;
    fvm_gnum_t  *send_gbuf = NULL;

    _get_faces_to_send(n2o_face_hist,
                       join_select->compact_rank_index,
                       &send_shift,
                       &send_gbuf);

    assert(send_shift[n_ranks] == n2o_face_hist->index[n2o_face_hist->n_elts]);

#if defined(HAVE_MPI)
    _exchange_faces(n_ranks,
                    send_shift,
                    send_gbuf,
                    join_select,
                    cs_glob_mpi_comm);
#endif

    for (i = 0; i < send_shift[n_ranks]; i++)
      cell_gnum[i] = send_gbuf[i];

    /* Free memory */

    BFT_FREE(send_gbuf);
    BFT_FREE(send_shift);;

  } /* End if n_ranks > 1 */

  /* Return pointer */

  *p_cell_gnum = cell_gnum;

}

/*----------------------------------------------------------------------------
 * Update mesh structure by adding new border faces after the face cutting
 * and deleting old one.
 *
 * parameters:
 *  join_select       <--  list of all implied entities in the joining op.
 *  join_mesh         <--  pointer to the local cs_join_mesh_t structure
 *  join2mesh_vtx_id  <--  relation between vertices in join_mesh/mesh
 *  n_new_b_faces     <--  local number of border faces after the joining
 *  new_face_type     <--  type (border/interior) of new faces
 *  n2o_face_hist     <--  face history structure (new -> old)
 *  mesh              <->  pointer of pointer to cs_mesh_t structure
 *---------------------------------------------------------------------------*/

static void
_add_new_border_faces(cs_join_select_t     *join_select,
                      cs_join_mesh_t       *join_mesh,
                      cs_int_t              join2mesh_vtx_id[],
                      cs_int_t              n_new_b_faces,
                      cs_join_face_type_t   new_face_type[],
                      cs_join_gset_t       *n2o_face_hist,
                      cs_mesh_t            *mesh)
{
  cs_int_t  i, j, select_id, vid, fid, shift, n_face_vertices;
  fvm_gnum_t  compact_old_fgnum;

  cs_int_t  n_ib_faces = mesh->n_b_faces, n_fb_faces = 0;
  fvm_gnum_t  n_g_ib_faces = mesh->n_g_b_faces;
  cs_int_t  *new_f2v_idx = NULL, *new_f2v_lst = NULL;
  cs_int_t  *new_face_family = NULL, *new_face_cells = NULL;
  fvm_gnum_t  *new_fgnum = NULL;

  const int  n_ranks = cs_glob_n_ranks;
  const int  rank = cs_glob_rank_id;
  const fvm_gnum_t  rank_start = join_select->compact_rank_index[rank] + 1;
  const fvm_gnum_t  rank_end = join_select->compact_rank_index[rank+1] + 1;

  n_fb_faces = n_ib_faces + n_new_b_faces - join_select->n_faces;
  mesh->n_b_faces = n_fb_faces;
  mesh->n_g_b_faces = n_fb_faces;

  BFT_MALLOC(new_f2v_idx, n_fb_faces + 1, cs_int_t);
  BFT_MALLOC(new_face_cells, n_fb_faces, cs_int_t);
  BFT_MALLOC(new_face_family, n_fb_faces, cs_int_t);

  if (n_ranks > 1)
    BFT_MALLOC(new_fgnum, n_fb_faces, fvm_gnum_t);

  /* Delete faces included in join_selection. Add other initial faces.
      - face -> vertex index
      - face -> cells connectivity
      - face family
      - face global num. (first pass)
  */

  n_fb_faces = 0;
  for (i = 0, select_id = 0; i < n_ib_faces; i++) {

    cs_bool_t  in_selection = false;

    if (select_id < join_select->n_faces) {
      if (i+1 == join_select->faces[select_id]) {
        in_selection = true;
        select_id++;
      }
    }

    if (in_selection == false) {

      if (n_ranks > 1)
        new_fgnum[n_fb_faces] = mesh->global_b_face_num[i];

      new_face_cells[n_fb_faces] = mesh->b_face_cells[i];
      new_face_family[n_fb_faces] = mesh->b_face_family[i];

      n_fb_faces++;
      n_face_vertices = mesh->b_face_vtx_idx[i+1] - mesh->b_face_vtx_idx[i];
      new_f2v_idx[n_fb_faces] = n_face_vertices;

    }

  }

  assert(n_fb_faces == n_ib_faces - join_select->n_faces);

  /* Add faces resulting from the joining operation */

  if (n_new_b_faces > 0) {
    for (i = 0; i < join_mesh->n_faces; i++) {
      if (new_face_type[i] == CS_JOIN_FACE_BORDER) {

        shift = n2o_face_hist->index[i];
        compact_old_fgnum = n2o_face_hist->g_list[shift];

        /* Initial selected border face must be in the selection */
        assert(rank_start <= compact_old_fgnum);
        assert(compact_old_fgnum < rank_end);

        fid = join_select->faces[compact_old_fgnum - rank_start] - 1;

        if (n_ranks > 1)
          new_fgnum[n_fb_faces] = join_mesh->face_gnum[i] + n_g_ib_faces;

        new_face_cells[n_fb_faces] = mesh->b_face_cells[fid];
        new_face_family[n_fb_faces] = mesh->b_face_family[fid];

        n_fb_faces++;
        n_face_vertices =
          join_mesh->face_vtx_idx[i+1] - join_mesh->face_vtx_idx[i];
        new_f2v_idx[n_fb_faces] = n_face_vertices;

      } /* If new border face */
    }
  } /* If n_new_b_faces > 0 */

  assert(mesh->n_b_faces == n_fb_faces);

  /* Build index */

  new_f2v_idx[0] = 1;
  for (i = 0; i < n_fb_faces; i++)
    new_f2v_idx[i+1] += new_f2v_idx[i];

  BFT_MALLOC(new_f2v_lst, new_f2v_idx[n_fb_faces]-1, cs_int_t);

  /* Define the face -> vertex connectivity */

  n_fb_faces = 0;
  for (i = 0, select_id = 0; i < n_ib_faces; i++) {

    cs_bool_t  in_selection = false;

    if (select_id < join_select->n_faces) {
      if (i+1 == join_select->faces[select_id]) {
        in_selection = true;
        select_id++;
      }
    }

    if (in_selection == false) {

      shift = new_f2v_idx[n_fb_faces]-1;

      for (j = mesh->b_face_vtx_idx[i]-1; j < mesh->b_face_vtx_idx[i+1]-1; j++)
        new_f2v_lst[shift++] = mesh->b_face_vtx_lst[j];

      n_fb_faces++;

    }

  }

  if (n_new_b_faces > 0) {
    for (i = 0; i < join_mesh->n_faces; i++) {
      if (new_face_type[i] == CS_JOIN_FACE_BORDER) {

        shift = new_f2v_idx[n_fb_faces] - 1;

        for (j = join_mesh->face_vtx_idx[i]-1;
             j < join_mesh->face_vtx_idx[i+1]-1; j++) {
          vid = join_mesh->face_vtx_lst[j] - 1;
          new_f2v_lst[shift++] = join2mesh_vtx_id[vid] + 1;
        }

        n_fb_faces++;

      }  /* If new border face */
    }
  } /* If n_new_b_faces > 0 */

  if (n_ranks > 1) { /* Get a compact global face numbering */

    fvm_io_num_t  *new_io_num = fvm_io_num_create(NULL,
                                                  new_fgnum,
                                                  n_fb_faces,
                                                  0); /* Not shared */

    const fvm_gnum_t  *new_io_gnum = fvm_io_num_get_global_num(new_io_num);

    mesh->n_g_b_faces = fvm_io_num_get_global_count(new_io_num);

    for (i = 0; i < n_fb_faces; i++)
      new_fgnum[i] = new_io_gnum[i];

    fvm_io_num_destroy(new_io_num);

    BFT_FREE(mesh->global_b_face_num);
    mesh->global_b_face_num = new_fgnum;

  }

  /* Free memory */

  BFT_FREE(mesh->b_face_vtx_idx);
  BFT_FREE(mesh->b_face_vtx_lst);
  BFT_FREE(mesh->b_face_cells);
  BFT_FREE(mesh->b_face_family);

  /* Update structure */

  mesh->b_face_vtx_idx = new_f2v_idx;
  mesh->b_face_vtx_lst = new_f2v_lst;
  mesh->b_face_cells = new_face_cells;
  mesh->b_face_family = new_face_family;
  mesh->b_face_vtx_connect_size = new_f2v_idx[n_fb_faces]-1;

}

/*----------------------------------------------------------------------------
 * Update mesh structure by adding new interior faces after the face split.
 *
 * parameters:
 *  join_select       <--  list of all implied entities in the joining op.
 *  join_mesh         <--  pointer to the local cs_join_mesh_t structure
 *  join2mesh_vtx_id  <--  relation between vertices in join_mesh/mesh
 *  cell_gnum         <--  global cell num. related to each initial face
 *  n_new_i_faces     <--  local number of interior faces after the joining
 *  default_family    <--  default family num to assign to each int. faces
 *  new_face_type     <--  type (border/interior) of new faces
 *  n2o_face_hist     <--  face history structure (new -> old)
 *  mesh              <->  pointer of pointer to cs_mesh_t structure
 *---------------------------------------------------------------------------*/

static void
_add_new_interior_faces(cs_join_select_t      *join_select,
                        cs_join_mesh_t        *join_mesh,
                        cs_int_t               join2mesh_vtx_id[],
                        fvm_gnum_t             cell_gnum[],
                        cs_int_t               n_new_i_faces,
                        cs_int_t               default_family,
                        cs_join_face_type_t    new_face_type[],
                        cs_join_gset_t        *n2o_face_hist,
                        cs_mesh_t             *mesh)
{
  cs_int_t  i, j, k, vid, shift, n_face_vertices, fid[2];
  fvm_gnum_t  compact_fgnum, cgnum[2];

  cs_int_t  n_fi_faces = 0, n_ii_faces = mesh->n_i_faces;
  cs_int_t  *new_f2v_idx = mesh->i_face_vtx_idx;
  cs_int_t  *new_f2v_lst = mesh->i_face_vtx_lst;
  cs_int_t  *new_face_family = mesh->i_face_family;
  cs_int_t  *new_face_cells = mesh->i_face_cells;
  fvm_gnum_t  n_g_ii_faces = mesh->n_g_i_faces;
  fvm_gnum_t  *new_fgnum = mesh->global_i_face_num;

  const int  n_ranks = cs_glob_n_ranks;
  const int  rank = cs_glob_rank_id;
  const fvm_gnum_t  loc_rank_s = join_select->compact_rank_index[rank];
  const fvm_gnum_t  loc_rank_e = join_select->compact_rank_index[rank+1];

  n_fi_faces = n_ii_faces + n_new_i_faces;
  mesh->n_i_faces = n_fi_faces;
  mesh->n_g_i_faces = n_fi_faces;

  BFT_REALLOC(new_f2v_idx, n_fi_faces + 1, cs_int_t);
  BFT_REALLOC(new_face_cells, 2*n_fi_faces, cs_int_t);
  BFT_REALLOC(new_face_family, n_fi_faces, cs_int_t);

  /* Add faces resulting from the joining operation
     - face -> vertex index
     - face -> cells connectivity
     - face family
  */

  n_fi_faces = n_ii_faces;
  for (i = 0; i < join_mesh->n_faces; i++) {

    if (new_face_type[i] == CS_JOIN_FACE_INTERIOR) {

      for (j = n2o_face_hist->index[i], k = 0;
           j < n2o_face_hist->index[i+1]; j++, k++) {

        cgnum[k] = cell_gnum[j];
        compact_fgnum = n2o_face_hist->g_list[j];

        if (loc_rank_s < compact_fgnum && compact_fgnum <= loc_rank_e)
          fid[k] = join_select->faces[compact_fgnum - loc_rank_s - 1] - 1;
        else
          fid[k] = -1;

      }

      if (cgnum[0] < cgnum[1]) { /* Keep the same order */

        if (fid[0] > -1)
          new_face_cells[2*n_fi_faces] = mesh->b_face_cells[fid[0]];
        else
          new_face_cells[2*n_fi_faces] = 0; /* Cell is on a distant rank */

        if (fid[1] > -1)
          new_face_cells[2*n_fi_faces+1] = mesh->b_face_cells[fid[1]];
        else
          new_face_cells[2*n_fi_faces+1] = 0; /* Cell is on a distant rank */

      }
      else {

        if (fid[0] > -1)
          new_face_cells[2*n_fi_faces+1] = mesh->b_face_cells[fid[0]];
        else
          new_face_cells[2*n_fi_faces+1] = 0; /* Cell is on a distant rank */

        if (fid[1] > -1)
          new_face_cells[2*n_fi_faces] = mesh->b_face_cells[fid[1]];
        else
          new_face_cells[2*n_fi_faces] = 0; /* Cell is on a distant rank */

      }

      new_face_family[n_fi_faces] = default_family; /* Default value.
                                                       TODO: Define real family
                                                    */

      n_fi_faces++;
      n_face_vertices = join_mesh->face_vtx_idx[i+1]-join_mesh->face_vtx_idx[i];
      new_f2v_idx[n_fi_faces] = n_face_vertices;

    }

  } /* End of loop on join_mesh faces */

  assert(mesh->n_i_faces == n_fi_faces);

  /* Build index */

  for (i = n_ii_faces; i < n_fi_faces; i++)
    new_f2v_idx[i+1] += new_f2v_idx[i];

  BFT_REALLOC(new_f2v_lst, new_f2v_idx[n_fi_faces]-1, cs_int_t);

  /* Define the face -> vertex connectivity list */

  n_fi_faces = n_ii_faces;
  for (i = 0; i < join_mesh->n_faces; i++) {
    if (new_face_type[i] == CS_JOIN_FACE_INTERIOR) {

      shift = new_f2v_idx[n_fi_faces]-1;

      for (j = join_mesh->face_vtx_idx[i]-1;
           j < join_mesh->face_vtx_idx[i+1]-1; j++) {
        vid = join_mesh->face_vtx_lst[j]-1;
        new_f2v_lst[shift++] = join2mesh_vtx_id[vid] + 1;
      }
      n_fi_faces++;

    }
  }

  if (n_ranks > 1) { /* Get a compact global face numbering */

    fvm_io_num_t  *new_io_num = NULL;
    const fvm_gnum_t  *new_io_gnum = NULL;

    BFT_REALLOC(new_fgnum, mesh->n_i_faces, fvm_gnum_t);

    n_fi_faces = n_ii_faces;
    for (i = 0; i < join_mesh->n_faces; i++)
      if (new_face_type[i] == CS_JOIN_FACE_INTERIOR)
        new_fgnum[n_fi_faces++] = join_mesh->face_gnum[i] + n_g_ii_faces;

    new_io_num = fvm_io_num_create(NULL, new_fgnum, n_fi_faces, 0);
    new_io_gnum = fvm_io_num_get_global_num(new_io_num);
    mesh->n_g_i_faces = fvm_io_num_get_global_count(new_io_num);

    for (i = 0; i < n_fi_faces; i++)
      new_fgnum[i] = new_io_gnum[i];

    fvm_io_num_destroy(new_io_num);

  }

  /* Update structure */

  mesh->i_face_vtx_idx = new_f2v_idx;
  mesh->i_face_vtx_lst = new_f2v_lst;
  mesh->global_i_face_num = new_fgnum;
  mesh->i_face_cells = new_face_cells;
  mesh->i_face_family = new_face_family;
  mesh->i_face_vtx_connect_size = new_f2v_idx[n_fi_faces]-1;

}

/*----------------------------------------------------------------------------
 * Get the family number to assign by default to the interior faces.
 *
 * parameters:
 *  mesh        <--  pointer to cs_mesh_t structure
 *
 * returns:
 *  the value of the default family number.
 *---------------------------------------------------------------------------*/

static cs_int_t
_get_default_family(cs_mesh_t   *mesh)
{
  int  i, j, n_grp, grp_num, grp_idx, n_colors, new_size;

  char **groups = NULL;
  int  *colors = NULL, *items = NULL;
  int  default_family = -1;

  BFT_MALLOC(groups, mesh->n_max_family_items, char*);
  BFT_MALLOC(colors, mesh->n_max_family_items, int);

  for (i = 0; i < mesh->n_families; i++) {

    n_colors = 0;
    n_grp  = 0;

    for (j = 0; j < mesh->n_max_family_items; j++) {

      if (mesh->family_item[j * mesh->n_families + i] > 0)
        colors[n_colors++] = mesh->family_item[j * mesh->n_families + i];

      else if (mesh->family_item[j * cs_glob_mesh->n_families + i] < 0) {
        grp_num = -mesh->family_item[j * mesh->n_families + i] - 1;
        grp_idx = mesh->group_idx[grp_num];
        groups[n_grp++] = mesh->group_lst + grp_idx -1;
      }

    }

    if (n_grp == 0 && n_colors == 0)
      default_family = i+1;

  } /* End of loop on families */

  BFT_FREE(groups);
  BFT_FREE(colors);

  if (default_family == -1) { /* Add a new family */

    mesh->n_families += 1;
    default_family = mesh->n_families;

    new_size = mesh->n_max_family_items *mesh->n_families;
    BFT_MALLOC(items, new_size, cs_int_t);

    for (i = 0; i < mesh->n_families - 1; i++)
      for (j = 0; j < mesh->n_max_family_items; j++)
        items[j * mesh->n_families + i] =
          mesh->family_item[j * mesh->n_families + i];

    for (j = 0; j < mesh->n_max_family_items; j++)
      items[j * (mesh->n_families-1) + (mesh->n_families-1)] = 0;

    BFT_FREE(mesh->family_item);
    mesh->family_item = items;

  }

  return default_family;

}

/*----------------------------------------------------------------------------
 * Re-orientation of joined faces. Check and correct if needed.
 *
 * parameters:
 *  mesh             <--  pointer of pointer to cs_mesh_t structure
 *  n_old_i_faces    <--  initial number of interior faces
 *  n_old_b_faces    <--  initial number of border faces
 *  n_select_cells   <--  number of cells in the selection
 *  cell_selection   -->  selection array (size: mesh->n_cells)
 *  select_cell_cen  -->  cell center for the selected cells
 *---------------------------------------------------------------------------*/

static void
_reorient_faces(cs_mesh_t    *mesh,
                cs_int_t      n_old_i_faces,
                cs_int_t      n_old_b_faces,
                cs_int_t      n_select_cells,
                cs_int_t      cell_selection[],
                cs_real_t     select_cell_cen[])
{
  cs_int_t  i, j, k, shift, s, e, vid, cid1, cid2, n_face_vertices;
  cs_real_t  dot_prod;

  cs_int_t  max_connect = 0;
  cs_int_t  *face_connect = NULL;
  cs_real_t  face_barycenter[3], face_normal[3], vect[3];
  cs_real_t  *face_vtx_coord = NULL;

  /* Remark:
     n_select_cells == join_select->n_faces (see _get_select_cell_cen() ) */

  /* Compute max face -> vertex connect. */

  for (i = n_old_b_faces - n_select_cells; i < mesh->n_b_faces; i++)
    max_connect = CS_MAX(max_connect,
                         mesh->b_face_vtx_idx[i+1] - mesh->b_face_vtx_idx[i]);

  for (i = n_old_i_faces; i < mesh->n_i_faces; i++)
    max_connect = CS_MAX(max_connect,
                         mesh->i_face_vtx_idx[i+1] - mesh->i_face_vtx_idx[i]);

  BFT_MALLOC(face_vtx_coord, 3*(max_connect+1), cs_real_t);
  BFT_MALLOC(face_connect, max_connect, cs_int_t);

  /* Border faces */

  for (i = n_old_b_faces - n_select_cells; i < mesh->n_b_faces; i++) {

    /* Define face_vtx_coord */

    s = mesh->b_face_vtx_idx[i] - 1;
    e = mesh->b_face_vtx_idx[i+1] - 1;
    n_face_vertices = e - s;
    shift = 0;

    for (j = s; j < e; j++) {

      vid = mesh->b_face_vtx_lst[j] - 1;
      face_connect[shift] = vid + 1;
      for (k = 0; k < 3; k++)
        face_vtx_coord[3*shift+k] = mesh->vtx_coord[3*vid+k];
      shift++;

    }

    vid = mesh->b_face_vtx_lst[s] - 1;
    for (k = 0; k < 3; k++)
      face_vtx_coord[3*shift+k] = mesh->vtx_coord[3*vid+k];

    /* Compute face barycenter and face normal */

    _get_face_quantity(n_face_vertices,
                       face_vtx_coord,
                       face_barycenter,
                       face_normal);  // unitary

    /*  Cb: cell barycenter
        Fb: face barycenter
        Nf: face normal
                            --->   ->
        Good orientation if CbFb . Nf > 0
        Else if reorientation is required.
    */

    cid1 = mesh->b_face_cells[i] - 1;
    assert(cid1 > -1);
    cid1 = cell_selection[cid1];
    assert(cid1 > -1);

    for (k = 0; k < 3; k++)
      vect[k] = face_barycenter[k] - select_cell_cen[3*cid1+k];

    _normalize(vect);
    dot_prod = _dot_product(vect, face_normal);

    if (dot_prod < 0.0)
      for (j = s, k = n_face_vertices-1; j < e; j++, k--)
        mesh->b_face_vtx_lst[j] = face_connect[k];

  } /* End of loop on border faces */

  /* Interior faces */

  for (i = n_old_i_faces; i < mesh->n_i_faces; i++) {

    /* Define face_vtx_coord */

    s = mesh->i_face_vtx_idx[i] - 1;
    e = mesh->i_face_vtx_idx[i+1] - 1;
    n_face_vertices = e - s;
    shift = 0;

    for (j = s; j < e; j++) {

      vid = mesh->i_face_vtx_lst[j] - 1;
      face_connect[shift] = vid + 1;
      for (k = 0; k < 3; k++)
        face_vtx_coord[3*shift+k] = mesh->vtx_coord[3*vid+k];
      shift++;

    }

    vid = mesh->i_face_vtx_lst[s] - 1;
    for (k = 0; k < 3; k++)
      face_vtx_coord[3*shift+k] = mesh->vtx_coord[3*vid+k];

    /* Compute face barycenter and face normal */

    _get_face_quantity(n_face_vertices,
                       face_vtx_coord,
                       face_barycenter,
                       face_normal);

    /*  Cb1: cell barycenter for cell 1
        Cb2: cell barycenter for cell 2
        Fb : face barycenter
        Nf : face normal

        if only Cb1 is available:
                              ---->   ->
          Good orientation if Cb1Fb . Nf > 0

        else:
                              ----->   ->
          Good orientation if Cb1Cb2 . Nf > 0
    */

    /* If we are on a parallel frontier: i_face_cell[] = 0 */

    cid1 = mesh->i_face_cells[2*i] - 1;
    cid2 = mesh->i_face_cells[2*i+1] - 1;

    if (cid1 > -1)
      cid1 = cell_selection[cid1];
    if (cid2 > -1)
      cid2 = cell_selection[cid2];

    if (cid2 < 0) { /* Parallel frontier. cid2 not available */

      assert(cs_glob_n_ranks > 1);

      for (k = 0; k < 3; k++)
        vect[k] = face_barycenter[k] - select_cell_cen[3*cid1+k];

      _normalize(vect);
      dot_prod = _dot_product(vect, face_normal);

    }
    else if (cid1 < 0) {

      assert(cs_glob_n_ranks > 1);

      for (k = 0; k < 3; k++)
        vect[k] = select_cell_cen[3*cid2+k] - face_barycenter[k];

      _normalize(vect);
      dot_prod = _dot_product(vect, face_normal);

    }
    else {

      for (k = 0; k < 3; k++)
        vect[k] = select_cell_cen[3*cid2+k] - select_cell_cen[3*cid1+k];

      _normalize(vect);
      dot_prod = _dot_product(vect, face_normal);


    }

    if (dot_prod < 0.0) /* Re-orient */
      for (j = s, k = n_face_vertices-1; j < e; j++, k--)
        mesh->i_face_vtx_lst[j] = face_connect[k];

  } /* End of loop on interior faces */

  /* Free memory */

  BFT_FREE(face_vtx_coord);
  BFT_FREE(face_connect);

}

/*----------------------------------------------------------------------------
 * Define a new face connectivty without redundant edge definition.
 *
 * parameters:
 *  s        <--  starting index in f2v_lst
 *  e        <--  ending index in f2v_lst
 *  f2v_lst  <--  face -> vertex connectivity list
 *  connect  <--  buffer to store locally the new face connectivity
 *  kill     <--  buffer to store vertex to delete from the connectivity
 *
 * returns:
 *  new number of vertices in the face connectivity.
 *---------------------------------------------------------------------------*/

static cs_int_t
_delete_edges(cs_int_t    s,
              cs_int_t    e,
              cs_int_t    f2v_lst[],
              cs_int_t    connect[],
              cs_int_t    kill[])
{
  cs_int_t  j, k, shift, count;

  cs_int_t  connect_size = e - s;

  /* Define local connectivity */

  for (j = s, k = 0; j < e; j++, k++)
    kill[k] = 0, connect[k] = f2v_lst[j];

  connect[k] = f2v_lst[s], kill[k++] = 0;
  connect[k] = f2v_lst[s+1], kill[k++] = 0;

  /* Find degenerated edges */

  count = 1;
  while (count > 0) {

    count = 0;
    for (j = 0; j < connect_size; j++) {
      if (connect[j] == connect[j+2]) {
        count++;
        kill[j] = 1;
        kill[(j+1) % connect_size] = 1;
      }
    }

    shift = 0;
    for (j = 0; j < connect_size; j++)
      if (kill[j] == 0)
        connect[shift++] = connect[j];

    connect_size = shift;
    connect[shift++] = connect[0];
    connect[shift++] = connect[1];

    for (j = 0; j < connect_size + 2; j++)
      kill[j] = 0;

  } /* End of while */

  /* Find empty edges */

  count = 1;
  while (count > 0) {

    count = 0;
    for (j = 0; j < connect_size; j++) {
      if (connect[j] == connect[j+1]) {
        count++;
        kill[(j+1) % connect_size] = 1;
      }
    }

    shift = 0;
    for (j = 0; j < connect_size; j++)
      if (kill[j] == 0)
        connect[shift++] = connect[j];

    connect_size = shift;
    connect[shift++] = connect[0];

    for (j = 0; j < connect_size + 1; j++)
      kill[j] = 0;

  } /* End of while */

  return connect_size;
}

/*============================================================================
 * Public function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Update mesh structure (vertices + faces) after the fusion step.
 *
 * parameters:
 *  join_param         <--  set of parameters for the joining operation
 *  join_select        <--  list of all implied entities in the joining op.
 *  o2n_vtx_gnum       <--  in : array on blocks on the new global vertex
 *                          out: local array on the new global vertex
 *  join_mesh          <--  pointer to the local cs_join_mesh_t structure
 *  mesh               <->  pointer of pointer to cs_mesh_t structure
 *---------------------------------------------------------------------------*/

void
cs_join_update_mesh_after_merge(cs_join_param_t      join_param,
                                cs_join_select_t    *join_select,
                                fvm_gnum_t            o2n_vtx_gnum[],
                                cs_join_mesh_t      *join_mesh,
                                cs_mesh_t           *mesh)
{
  cs_int_t  i, j, shift, select_id, adj_id, old_id, s_id, new_num, n_sel_faces;

  cs_int_t  *o2n_vtx_id = NULL, *join2mesh_vtx_id = NULL, *sel_faces = NULL;
  fvm_gnum_t  *old_vtx_gnum = NULL;
  edge_builder_t  *edge_builder = NULL;

  const cs_int_t  n_bm_vertices = mesh->n_vertices; /* bf: before fusion */
  const cs_int_t  n_ranks = cs_glob_n_ranks;

  edge_builder = _init_edge_builder(join_select, mesh);

  /* Build an array keeping relation between old/new global vertex num. */

  if (n_ranks == 1) {

    fvm_gnum_t  *loc_vtx_gnum = NULL;

    BFT_MALLOC(loc_vtx_gnum, n_bm_vertices, fvm_gnum_t);

    /* Initialize array */

    for (i = 0; i < n_bm_vertices; i++)
      loc_vtx_gnum[i] = i+1;

    /* Update value for selected vertices */

    for (i = 0, shift = 0;
         i < n_bm_vertices && shift < join_select->n_vertices; i++) {
      if (i + 1 == join_select->vertices[shift])
        loc_vtx_gnum[i] = o2n_vtx_gnum[shift++];
    }

    BFT_FREE(o2n_vtx_gnum);
    o2n_vtx_gnum = loc_vtx_gnum;

  } /* End if serial mode */

#if defined(HAVE_MPI)
  if (n_ranks > 1)
    _get_local_o2n_vtx_gnum(mesh, &o2n_vtx_gnum);
#endif

  /* Update mesh structure. Define new vertices */

  _update_vertices_after_merge(o2n_vtx_gnum,
                               join_mesh,
                               mesh,
                               &join2mesh_vtx_id, /* size: join_mesh->n_vertices */
                               &o2n_vtx_id,       /* size: n_bm_vertices */
                               &old_vtx_gnum);    /* size: n_bm_vertices */

  /* Define the evolution of each initial edge in edge_builder_t struct. */

  _complete_edge_builder(join_select,
                         join_mesh,
                         mesh,
                         o2n_vtx_gnum,
                         join2mesh_vtx_id,
                         edge_builder);

#if defined(HAVE_MPI)
  if (join_select->do_single_sync == true)
    _sync_single_elements(join_select,
                          n_bm_vertices,
                          old_vtx_gnum,
                          o2n_vtx_id,
                          join_mesh->n_vertices,
                          join2mesh_vtx_id,
                          edge_builder,
                          mesh);
#endif

  BFT_FREE(old_vtx_gnum);

#if 0 && defined(DEBUG) && !defined(NDEBUG) /* Dump the structure */
  if (join_param.verbosity > 2) {
    int k;

    bft_printf("\n  Dump edge_builder_t structure (%p)\n", edge_builder);

    if (edge_builder != NULL) {
      bft_printf("  n_vertices: %10d\n"
                 "  n_edges   : %10d\n",
                 edge_builder->n_vertices, edge_builder->n_edges);

      for (i = 0; i < edge_builder->n_vertices; i++) {

        bft_printf("%9d - [%10.4f %10.4f %10.4f]: (%d, %d) v-v:",
                   i+1, mesh->vtx_coord[3*i], mesh->vtx_coord[3*i+1],
                   mesh->vtx_coord[3*i+2], edge_builder->v2v_idx[i],
                   edge_builder->v2v_idx[i+1]);

        for (j = edge_builder->v2v_idx[i]; j < edge_builder->v2v_idx[i+1]; j++) {
          bft_printf(" %d (", edge_builder->v2v_lst[j]);
          for (k = edge_builder->v2v_sub_idx[j];
               k < edge_builder->v2v_sub_idx[j+1]; k++)
            bft_printf("%d ", edge_builder->v2v_sub_lst[k]);
          bft_printf(") ");
        }
        bft_printf("\n");

      }
      bft_printf_flush();
    }

  }
#endif

  /* Update connectivity for the selected faces */

  _update_selected_face_connect(join_select,
                                join_mesh,
                                join2mesh_vtx_id,
                                mesh->n_b_faces,
                                &(mesh->b_face_vtx_idx),
                                &(mesh->b_face_vtx_lst));

  /* Partial free memory */

  BFT_FREE(o2n_vtx_gnum);

  /* Update adjancent border face connectivity */

  if (join_select->n_b_s_faces > 0)
    _get_select_face_lst(join_select->n_b_adj_faces,
                         join_select->b_adj_faces,
                         join_select->n_b_s_faces,
                         join_select->b_s_faces,
                         &n_sel_faces,
                         &sel_faces);

  else {
    n_sel_faces = join_select->n_b_adj_faces;
    sel_faces = join_select->b_adj_faces;
  }

  _update_adj_face_connect(n_sel_faces,
                           sel_faces,
                           edge_builder,
                           o2n_vtx_id,
                           mesh->n_b_faces,
                           &(mesh->b_face_vtx_idx),
                           &(mesh->b_face_vtx_lst));

  if (join_select->n_b_s_faces > 0)
    BFT_FREE(sel_faces);

  /* Update adjancent interior face connectivity */

  if (join_select->n_i_s_faces > 0)
    _get_select_face_lst(join_select->n_i_adj_faces,
                         join_select->i_adj_faces,
                         join_select->n_i_s_faces,
                         join_select->i_s_faces,
                         &n_sel_faces,
                         &sel_faces);

  else {
    n_sel_faces = join_select->n_i_adj_faces;
    sel_faces = join_select->i_adj_faces;
  }

  _update_adj_face_connect(n_sel_faces,
                           sel_faces,
                           edge_builder,
                           o2n_vtx_id,
                           mesh->n_i_faces,
                           &(mesh->i_face_vtx_idx),
                           &(mesh->i_face_vtx_lst));

  if (join_select->n_i_s_faces > 0)
    BFT_FREE(sel_faces);

  /* Free memory */

  BFT_FREE(edge_builder->v2v_idx);
  BFT_FREE(edge_builder->v2v_lst);
  BFT_FREE(edge_builder->v2v_sub_idx);
  BFT_FREE(edge_builder->v2v_sub_lst);
  BFT_FREE(edge_builder);

  /* Update initial face connectivity for the remaining faces */

  for (i = 0, select_id = 0, adj_id = 0, s_id = 0; i < mesh->n_b_faces; i++) {

    cs_bool_t  do_update = true;

    if (select_id < join_select->n_faces) {
      if (i+1 == join_select->faces[select_id]) {
        do_update = false; /* Already done */
        select_id++;
      }
    }

    if (adj_id < join_select->n_b_adj_faces) {
      if (i+1 == join_select->b_adj_faces[adj_id]) {
        do_update = false; /* Already done */
        adj_id++;
      }
    }

    if (s_id < join_select->n_b_s_faces) {
      if (i+1 == join_select->b_s_faces[s_id]) {
        do_update = false; /* Already done */
        s_id++;
      }
    }

    if (do_update == true) {

      for (j = mesh->b_face_vtx_idx[i]-1;
           j < mesh->b_face_vtx_idx[i+1]-1; j++) {

        old_id = mesh->b_face_vtx_lst[j] - 1;
        new_num = o2n_vtx_id[old_id] + 1;
        mesh->b_face_vtx_lst[j] = new_num;

      }

    }

  } /* End of loop on border faces */

  for (i = 0, adj_id = 0, s_id = 0; i < mesh->n_i_faces; i++) {

    cs_bool_t  do_update = true;

    if (adj_id < join_select->n_i_adj_faces) {
      if (i+1 == join_select->i_adj_faces[adj_id]) {
        do_update = false; /* Already done */
        adj_id++;
      }
    }

    if (s_id < join_select->n_i_s_faces) {
      if (i+1 == join_select->i_s_faces[s_id]) {
        do_update = false; /* Already done */
        s_id++;
      }
    }

    if (do_update == true) {

      for (j = mesh->i_face_vtx_idx[i]-1;
           j < mesh->i_face_vtx_idx[i+1]-1; j++) {

        old_id = mesh->i_face_vtx_lst[j] - 1;
        new_num = o2n_vtx_id[old_id] + 1;
        mesh->i_face_vtx_lst[j] = new_num;

      }

    }

  } /* End of loop on interior faces */

  /* Update the cs_join_select_t structure */

  for (i = 0; i < join_select->n_vertices; i++) {

    old_id = join_select->vertices[i] - 1;
    new_num = o2n_vtx_id[old_id] + 1;
    join_select->vertices[i] = new_num;

  }

  /* Post if required */

  cs_join_post_after_merge(join_param, join_select);

  /* Free memory */

  BFT_FREE(join2mesh_vtx_id);
  BFT_FREE(o2n_vtx_id);

}

/*----------------------------------------------------------------------------
 * Update mesh structure (vertices + faces) after the face split step.
 *
 * parameters:
 *  join_param        <--  set of parameters for the joining operation
 *  join_select       <--  list of all implied entities in the joining op.
 *  o2n_face_hist     <--  relation between faces before/after the joining
 *  n_select_cells    <--  number of cells in the selection
 *  cell_selection    -->  selection array (size: mesh->n_cells)
 *  select_cell_cen   -->  cell center for the selected cells
 *  join_mesh         <--  pointer to the local cs_join_mesh_t structure
 *  mesh              <->  pointer of pointer to cs_mesh_t structure
 *---------------------------------------------------------------------------*/

void
cs_join_update_mesh_after_split(cs_join_param_t      join_param,
                                cs_join_select_t    *join_select,
                                cs_join_gset_t      *o2n_face_hist,
                                cs_join_mesh_t       *join_mesh,
                                cs_mesh_t            *mesh)
{
  int  i, n_matches, default_family;

  fvm_gnum_t  n_g_new_b_faces = 0;
  cs_int_t  n_new_i_faces = 0, n_new_b_faces = 0;
  cs_int_t  n_old_i_faces = mesh->n_i_faces;
  cs_int_t  n_old_b_faces = mesh->n_b_faces;
  cs_int_t  *join2mesh_vtx_id = NULL;
  fvm_gnum_t  *cell_gnum = NULL;
  cs_join_face_type_t  *new_face_type = NULL;
  cs_join_gset_t  *n2o_face_hist = NULL;
  cs_int_t  n_select_cells = join_select->n_faces; /* because border faces */
  cs_int_t  *cell_selection = join_select->cell_filter;
  cs_real_t  *select_cell_cen = join_select->cell_cen;

  const int  n_ranks = cs_glob_n_ranks;

  assert(mesh != NULL);

  /* Invert face historic */

  n2o_face_hist = cs_join_gset_invert(o2n_face_hist);

#if defined(HAVE_MPI)
  if (n_ranks > 1) {

    cs_join_gset_t  *n2o_sync_block = NULL;
    MPI_Comm  mpi_comm = cs_glob_mpi_comm;

    n2o_sync_block = cs_join_gset_sync_by_block(join_mesh->n_g_faces,
                                                n2o_face_hist,
                                                mpi_comm);

    cs_join_gset_update_from_block(join_mesh->n_g_faces,
                                   n2o_sync_block,
                                   n2o_face_hist,
                                   mpi_comm);

    cs_join_gset_destroy(&n2o_sync_block);

  }
#endif

#if 1 && defined(DEBUG) && !defined(NDEBUG)
  if (join_param.verbosity > 3) { /* Full dump of structures */

    int  len;
    FILE  *dbg_file = NULL;
    char  *filename = NULL;

    len = strlen("JoinDBG_n2oFaceHist.dat")+1+2+4;
    BFT_MALLOC(filename, len, char);
    sprintf(filename, "Join%02dDBG_n2oFaceHist%04d.dat",
            join_param.num, cs_glob_rank_id);
    dbg_file = fopen(filename, "w");

    cs_join_gset_dump(dbg_file, n2o_face_hist);

    fflush(dbg_file);
    BFT_FREE(filename);
    fclose(dbg_file);

  }
#endif

  /* Get new subfaces evolution */

  assert(n2o_face_hist->n_elts == join_mesh->n_faces);
  BFT_MALLOC(new_face_type, join_mesh->n_faces, cs_join_face_type_t);

  for (i = 0; i < join_mesh->n_faces; i++) {

    assert(join_mesh->face_gnum[i] == n2o_face_hist->g_elts[i]);

    n_matches = n2o_face_hist->index[i+1] - n2o_face_hist->index[i];

    if (n_matches == 1) {
      n_new_b_faces += 1;
      new_face_type[i] = CS_JOIN_FACE_BORDER;
    }
    else if (n_matches == 2) {
      n_new_i_faces += 1;
      new_face_type[i] = CS_JOIN_FACE_INTERIOR;
    }
    else
      bft_error(__FILE__, __LINE__, 0,
                _("  Incompatible new face type found.\n"
                  "  n_matches is not to 1 or 2: n_matches = %d\n"),
                n_matches);

  }

  if (join_param.verbosity > 0)
    bft_printf(_("\n  Local configuration after the joining operation:\n"
                 "    Number of interior faces to add: %9d\n"
                 "    Number of border faces to add  : %9d\n"),
               n_new_i_faces, n_new_b_faces);

  if (n_ranks == 1)
    n_g_new_b_faces = n_new_b_faces;

#if defined(HAVE_MPI)
  if (n_ranks > 1) {

    MPI_Comm  mpi_comm = cs_glob_mpi_comm;

    MPI_Allreduce(&n_new_b_faces, &n_g_new_b_faces, 1, FVM_MPI_GNUM,
                  MPI_SUM, mpi_comm);

    if (cs_glob_rank_id == 0)
      bft_printf(_("\n  Global configuration after the joining operation:\n"
                   "     Global number of border faces to add : %10u\n"),
                 n_g_new_b_faces);

  }
#endif

  /* Get the family number to assign by default to the interior faces */

  default_family = _get_default_family(mesh);

  if (join_param.verbosity > 1)
    bft_printf("\n  Default family for interior joined faces: %d\n",
               default_family);

  /* Define join2mesh_vtx_id and add new vertices (only in parallel mode)
     These vertices already exist but only on other ranks. During the face
     cutting op., these vertices are used to define the new face connect. */

  _update_vertices_after_split(join_mesh, mesh, &join2mesh_vtx_id);

  /* Get associated global cell number */

  _get_linked_cell_gnum(join_select, n2o_face_hist, &cell_gnum);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n  List of linked global cell number\n");
  for (i = 0; i < n2o_face_hist->index[n2o_face_hist->n_elts]; i++)
    bft_printf(" %9d - %10u - %10u\n",
               i, n2o_face_hist->g_list[i], cell_gnum[i]);
  bft_printf_flush();
#endif

  /*  Update mesh structure:
        - Update first the interior faces because we need the global
          numbering of the initial border faces
        - Then update border faces
  */

  _add_new_interior_faces(join_select,
                          join_mesh,
                          join2mesh_vtx_id,
                          cell_gnum,
                          n_new_i_faces,
                          default_family,
                          new_face_type,
                          n2o_face_hist,
                          mesh);

  _add_new_border_faces(join_select,
                        join_mesh,
                        join2mesh_vtx_id,
                        n_new_b_faces,
                        new_face_type,
                        n2o_face_hist,
                        mesh);

  /* Re-orientaton of joined faces if needed */

  _reorient_faces(mesh,
                  n_old_i_faces,
                  n_old_b_faces,
                  n_select_cells,
                  cell_selection,
                  select_cell_cen);

  /* Update global vertex information */

  if (n_ranks == 1) {
    mesh->n_g_vertices = mesh->n_vertices;
    BFT_FREE(mesh->global_vtx_num);
  }

#if defined(HAVE_MPI)
  if (n_ranks > 1) { /* Define a new compact global vertex numbering */

    fvm_io_num_t  *vtx_io_num = fvm_io_num_create(NULL,
                                                  mesh->global_vtx_num,
                                                  mesh->n_vertices,
                                                  0); /* Not shared */
    const fvm_gnum_t  n_g_io_num = fvm_io_num_get_global_count(vtx_io_num);
    const fvm_gnum_t  *io_gnum = fvm_io_num_get_global_num(vtx_io_num);

    mesh->n_g_vertices = n_g_io_num;

    for (i = 0; i < mesh->n_vertices; i++)
      mesh->global_vtx_num[i] = io_gnum[i];

    fvm_io_num_destroy(vtx_io_num);

  }
#endif

  /* Free memory */

  BFT_FREE(new_face_type);
  BFT_FREE(cell_gnum);
  BFT_FREE(join2mesh_vtx_id);

  cs_join_gset_destroy(&n2o_face_hist);

  /* Post if required */

  cs_join_post_after_split(n_old_i_faces, n_new_i_faces,
                           n_old_b_faces, n_new_b_faces,
                           n_g_new_b_faces,
                           join_select->n_faces,
                           mesh,
                           join_param);
}

/*----------------------------------------------------------------------------
 * Clean a cs_mesh_t struct.  (delete redundant and empty edge definition)
 *
 * parameters:
 *  param     <--  set of parameters for the joining operation
 *  mesh      <--  pointer to a cs_mesh_t struct.
 *---------------------------------------------------------------------------*/

void
cs_join_update_mesh_clean(cs_join_param_t   param,
                          cs_mesh_t        *mesh)
{
  cs_int_t  i, j, s, e, n_vertices, n_init_vertices, connect_size;

  cs_int_t  connect_shift = 0;
  cs_int_t  max_connect = 0, b_size = 10, i_size = 10;
  cs_int_t  n_b_clean_faces = 0, n_i_clean_faces = 0;
  cs_int_t  *b_clean_faces = NULL, *i_clean_faces = NULL;
  cs_int_t  *kill = NULL, *connect = NULL;

  for (i = 0; i < mesh->n_b_faces; i++)
    max_connect = CS_MAX(max_connect,
                         mesh->b_face_vtx_idx[i+1] - mesh->b_face_vtx_idx[i]);

  for (i = 0; i < mesh->n_i_faces; i++)
    max_connect = CS_MAX(max_connect,
                         mesh->i_face_vtx_idx[i+1] - mesh->i_face_vtx_idx[i]);

  BFT_MALLOC(kill, max_connect + 2, cs_int_t);
  BFT_MALLOC(connect, max_connect + 2, cs_int_t);

  if (param.verbosity > 1) {
    BFT_MALLOC(b_clean_faces, b_size, cs_int_t);
    BFT_MALLOC(i_clean_faces, i_size, cs_int_t);
  }

  /* Border face treatment */

  for (i = 0; i < mesh->n_b_faces; i++) {

    s = mesh->b_face_vtx_idx[i] - 1;
    e = mesh->b_face_vtx_idx[i+1] - 1;
    n_init_vertices = e - s;
    connect_size = -1;

    n_vertices = n_init_vertices;
    while (connect_size != n_vertices) {

      connect_size = _delete_edges(s, e, mesh->b_face_vtx_lst, connect, kill);
      assert(connect_size <= n_vertices);

      if (connect_size != n_vertices) {
        n_vertices = connect_size;
        connect_size += 1;
      }
      else
        n_vertices = connect_size;

    }

    if (n_init_vertices != n_vertices) {

      if (param.verbosity > 1) {

        bft_printf(_("  Clean border face %d. New number of vertices: %d\n"),
                   i+1, n_vertices);

        if (n_b_clean_faces + 1 > b_size) {
          b_size *= 2;
          BFT_REALLOC(b_clean_faces, b_size, cs_int_t);
        }
        b_clean_faces[n_b_clean_faces] = i+1;

      }
      n_b_clean_faces++;

    }

    for (j = 0; j < n_vertices; j++)
      mesh->b_face_vtx_lst[connect_shift++] = connect[j];
    mesh->b_face_vtx_idx[i] = connect_shift;

  } /* End of loop on border faces */

  if (param.verbosity > 0)
    bft_printf(_("\n  Degenerated connectivity for %d final border faces.\n"),
               n_b_clean_faces);

  for (i = mesh->n_b_faces; i > 0; i--)
    mesh->b_face_vtx_idx[i] = mesh->b_face_vtx_idx[i-1] + 1;
  mesh->b_face_vtx_idx[0] = 1;

  BFT_REALLOC(mesh->b_face_vtx_lst, mesh->b_face_vtx_idx[mesh->n_b_faces],
              cs_int_t);

  /* Interior face treatment */

  connect_shift = 0;
  for (i = 0; i < mesh->n_i_faces; i++) {

    s = mesh->i_face_vtx_idx[i] - 1;
    e = mesh->i_face_vtx_idx[i+1] - 1;
    n_init_vertices = e - s;
    connect_size = -1;

    n_vertices = n_init_vertices;
    while (connect_size != n_vertices) {

      connect_size = _delete_edges(s, e, mesh->i_face_vtx_lst, connect, kill);
      assert(connect_size <= n_vertices);

      if (connect_size != n_vertices) {
        n_vertices = connect_size;
        connect_size += 1;
      }
      else
        n_vertices = connect_size;

    }

    if (n_init_vertices != n_vertices) {

      if (param.verbosity > 1) {

        bft_printf(_("  Clean interior face %d. New number of vertices: %d\n"),
                   i+1, n_vertices);

        if (n_i_clean_faces + 1 > i_size) {
          i_size *= 2;
          BFT_REALLOC(i_clean_faces, i_size, cs_int_t);
        }
        i_clean_faces[n_i_clean_faces] = i+1;

      }
      n_i_clean_faces++;

    }

    for (j = 0; j < n_vertices; j++)
      mesh->i_face_vtx_lst[connect_shift++] = connect[j];
    mesh->i_face_vtx_idx[i] = connect_shift;

  } /* End of loop on interior faces */

  if (param.verbosity > 0)
    bft_printf(_("  Degenerated connectivity for %d final interior faces.\n"
                 "  Mesh cleaning done.\n"), n_i_clean_faces);

  for (i = mesh->n_i_faces; i > 0; i--)
    mesh->i_face_vtx_idx[i] = mesh->i_face_vtx_idx[i-1] + 1;
  mesh->i_face_vtx_idx[0] = 1;

  BFT_REALLOC(mesh->i_face_vtx_lst, mesh->i_face_vtx_idx[mesh->n_i_faces],
              cs_int_t);

  if (param.verbosity > 1) { /* Post-treat clean faces */

    if (n_i_clean_faces > 0 || n_b_clean_faces > 0) {

      BFT_REALLOC(i_clean_faces, n_i_clean_faces, cs_int_t);
      BFT_REALLOC(b_clean_faces, n_b_clean_faces, cs_int_t);

      cs_join_post_cleaned_faces(n_i_clean_faces,
                                 i_clean_faces,
                                 n_b_clean_faces,
                                 b_clean_faces,
                                 param);

    }

    BFT_FREE(b_clean_faces);
    BFT_FREE(i_clean_faces);

  }

  /* Free memory */

  BFT_FREE(kill);
  BFT_FREE(connect);

}

/*---------------------------------------------------------------------------*/

END_C_DECLS
