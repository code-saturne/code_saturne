/*============================================================================
 * Management of conforming and non-conforming join
 *===========================================================================*/

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
 *---------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *---------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_io_num.h"
#include "fvm_periodicity.h"

#include "cs_order.h"
#include "cs_search.h"
#include "cs_sort.h"
#include "cs_join_perio.h"
#include "cs_join_post.h"
#include "cs_join_util.h"
#include "cs_mesh_group.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *---------------------------------------------------------------------------*/

#include "cs_join_update.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *===========================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *===========================================================================*/

/* Local structure to update the mesh connectivty after the merge step */

typedef struct {

  cs_lnum_t   n_vertices;  /* v2v_idx size - 1 */
  cs_lnum_t   n_edges;     /* v2v_lst size */

  cs_lnum_t  *v2v_idx;
  cs_lnum_t  *v2v_lst;

  cs_lnum_t  *v2v_sub_idx;
  cs_lnum_t  *v2v_sub_lst;  /* if -1: edge has been deleted */

} edge_builder_t;

/*============================================================================
 * Private function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Compute the cross product of two vectors.
 *
 * parameters:
 *   v1     <--  first vector
 *   v2     <--  second vector
 *   result -->  cross product v1xv2
 *
 * returns:
 *   the resulting cross product (v1 x v2)
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
 *   v1 <-- first vector
 *   v2 <-- second vector
 *
 * returns:
 *   the resulting dot product (v1.v2)
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
 *   v <-> vector to normalize
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
 * Compute face normal for the current face.
 *
 * parameters:
 *   n_face_vertices <-- number of vertices defining the face
 *   face_vtx_coord  <-- coordinates of each vertex of the face
 *                        (size: n_face_vertices+1)
 *   normal          <-> normal of the face (|normal| = 1)
 *----------------------------------------------------------------------------*/

static void
_get_face_normal(cs_lnum_t         n_face_vertices,
                 const cs_real_t   face_vtx_coord[],
                 cs_real_t         normal[])
{
  cs_lnum_t  i, coord;
  cs_real_t  v1[3], v2[3], tri_normal[3], barycenter[3];

  cs_real_t  inv_n_face_vertices = 1/(double)n_face_vertices;

  /* Initialization */

  for (coord = 0; coord < 3; coord++) {
    normal[coord] = 0.0;
    barycenter[coord] = 0.0;
  }

  /* Compute face barycenter */

  for (i = 0; i < n_face_vertices; i++)
    for (coord = 0; coord < 3; coord++)
      barycenter[coord] += face_vtx_coord[3*i+coord];

  for (coord = 0; coord < 3; coord++)
    barycenter[coord] *= inv_n_face_vertices;

  /* Compute triangle normal and update face normal */

  for (i = 0; i < n_face_vertices; i++) {

    for (coord = 0; coord < 3; coord++) {
      v1[coord] = face_vtx_coord[3*i     + coord] - barycenter[coord];
      v2[coord] = face_vtx_coord[3*(i+1) + coord] - barycenter[coord];
    }

    _cross_product(v1, v2, tri_normal);

    for (coord = 0; coord < 3; coord++) {
      tri_normal[coord] *= 0.5;
      normal[coord] += tri_normal[coord];
    }

  }

  _normalize(normal);
}

/*----------------------------------------------------------------------------
 * Get the related edge id in edge_builder_t structure from a couple of
 * vertex ids.
 *
 * parameters:
 *   v1_id        <-- first vertex id
 *   v2_id        <-- second vertex id
 *   edge_builder <-- pointer to an edge_builder_t structure
 *
 * returns:
 *   related edge_id in edge_builder_t structure
 *---------------------------------------------------------------------------*/

inline static cs_lnum_t
_get_join_edge_id(cs_lnum_t              v1_id,
                  cs_lnum_t              v2_id,
                  const edge_builder_t  *edge_builder)
{
  cs_lnum_t  i, va, vb;
  cs_lnum_t  edge_id = -1;

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
 * Update elements of a cs_mesh_t structure related to the vertices after the
 * merge step.
 *
 * parameters:
 *   selection  --> list of all implied entities in the joining op.
 *   o2n_vtx_id <-> relation between init. and current local num.
 *   mesh       <-> pointer of pointer to cs_mesh_t structure
 *---------------------------------------------------------------------------*/

static void
_sync_single_vertices(const cs_join_select_t  *selection,
                      cs_lnum_t                o2n_vtx_id[],
                      cs_mesh_t               *mesh)
{
  cs_lnum_t  i, s, e, rank, shift, length, request_count;

  double  *s_buf = NULL, *c_buf = NULL;
  cs_join_sync_t  *s_vertices = selection->s_vertices;
  cs_join_sync_t  *c_vertices = selection->c_vertices;
  MPI_Request  *request = NULL;
  MPI_Status   *status = NULL;
  MPI_Comm  mpi_comm = cs_glob_mpi_comm;

  const int  loc_rank = CS_MAX(cs_glob_rank_id, 0);

  bft_printf("\n  Synchronization of the \"single\" elements after the merge"
             " step.\n");
  bft_printf_flush();

  assert(cs_glob_n_ranks > 1);

  /* Allocate MPI buffers used for exchanging data */

  BFT_MALLOC(request, c_vertices->n_ranks + s_vertices->n_ranks, MPI_Request);
  BFT_MALLOC(status, c_vertices->n_ranks + s_vertices->n_ranks, MPI_Status);

  /* Synchronization of vertex coordinates */

  BFT_MALLOC(s_buf, 3*s_vertices->n_elts, double);

  /* Receive data from distant ranks */

  request_count = 0;

  for (i = 0; i < s_vertices->n_ranks; i++) {

    rank = s_vertices->ranks[i];
    s = s_vertices->index[i];
    e = s_vertices->index[i+1];
    length = 3*(e-s);

    MPI_Irecv(&(s_buf[3*s]),
              length,
              MPI_DOUBLE,
              rank,
              rank,
              mpi_comm,
              &(request[request_count++]));

  }

  /* We wait for posting all receives (often recommended) */

  MPI_Barrier(mpi_comm);

  /* Build c_buf = buffer to send */

  BFT_MALLOC(c_buf, 3*c_vertices->n_elts, double);

  for (shift = 0, i = 0; i < c_vertices->n_elts; i++) {

    cs_lnum_t  new_id = o2n_vtx_id[c_vertices->array[i]-1];

    c_buf[shift++] = mesh->vtx_coord[3*new_id];
    c_buf[shift++] = mesh->vtx_coord[3*new_id+1];
    c_buf[shift++] = mesh->vtx_coord[3*new_id+2];

  }

  /* Send data to distant ranks */

  for (i = 0; i < c_vertices->n_ranks; i++) {

    rank = c_vertices->ranks[i];
    s = c_vertices->index[i];
    e = c_vertices->index[i+1];
    length = 3*(e-s);

    MPI_Isend(&(c_buf[3*s]),
              length,
              MPI_DOUBLE,
              rank,
              loc_rank,
              mpi_comm,
              &(request[request_count++]));

  }

  /* Wait for all exchanges */

  MPI_Waitall(request_count, request, status);

  /* Update vertex coordinates */

  for (shift = 0, i = 0; i < s_vertices->n_elts; i++) {

    cs_lnum_t  new_id = o2n_vtx_id[s_vertices->array[i]-1];

    mesh->vtx_coord[3*new_id]   = s_buf[shift++];
    mesh->vtx_coord[3*new_id+1] = s_buf[shift++];
    mesh->vtx_coord[3*new_id+2] = s_buf[shift++];

  }

  /* Free memory */

  BFT_FREE(c_buf);
  BFT_FREE(s_buf);
  BFT_FREE(request);
  BFT_FREE(status);

}

/*----------------------------------------------------------------------------
 * Update cs_join_edges_t through the synchronization of "single" edges.
 *
 * parameters:
 *  selection        --> list of all implied entities in the joining op.
 *  n_bm_vertices    --> number of vertices in mesh before the merge step
 *  o2n_vtx_id       <-> relation between init. and current local num.
 *  n_j_vertices     --> number of vertices in join_mesh
 *  join2mesh_vtx_id <-> relation between join mesh and after merge vertex
 *  join_edges       <-> pointer to a cs_join_edges_t structure
 *  mesh             <-> pointer of pointer to cs_mesh_t structure
 *---------------------------------------------------------------------------*/

static void
_sync_single_edges(const cs_join_select_t   *selection,
                   cs_lnum_t                 n_bm_vertices,
                   cs_lnum_t                 o2n_vtx_id[],
                   cs_lnum_t                 n_j_vertices,
                   cs_lnum_t                 join2mesh_vtx_id[],
                   edge_builder_t           *edge_builder,
                   cs_mesh_t                *mesh)
{
  int  i, j, k, id, vid, rank, shift, edge_id, vid1, vid2;
  int  length, request_count, distant_rank;
  cs_gnum_t cur, prev;

  int  c_sub_size = 0, s_sub_size = 0, n_new_vertices = 0;
  cs_lnum_t  *new_v2v_sub_idx = NULL, *new_v2v_sub_lst = NULL;
  int  *s_count = NULL, *s_sub_index = NULL;
  int  *c_count = NULL, *c_sub_index = NULL;
  cs_gnum_t *c_sub_gbuf = NULL, *s_sub_gbuf = NULL;
  cs_gnum_t *new_vtx_gnum = NULL;
  cs_real_t *c_sub_coord = NULL, *s_sub_coord = NULL, *new_coord = NULL;

  MPI_Request  *request = NULL;
  MPI_Status   *status = NULL;
  MPI_Comm  mpi_comm = cs_glob_mpi_comm;

  const int  loc_rank = CS_MAX(cs_glob_rank_id, 0);
  const cs_join_sync_t  *s_edges = selection->s_edges;
  const cs_join_sync_t  *c_edges = selection->c_edges;

  assert(cs_glob_n_ranks > 1);

  /* Allocate MPI buffers used for exchanging data */

  BFT_MALLOC(request, c_edges->n_ranks + s_edges->n_ranks, MPI_Request);
  BFT_MALLOC(status, c_edges->n_ranks + s_edges->n_ranks, MPI_Status);

  /* Get a number of sub-element for each received edge */

  BFT_MALLOC(s_count, s_edges->n_elts, int);
  BFT_MALLOC(c_count, c_edges->n_elts, int);

  for (i = 0; i < s_edges->n_elts; i++)
    s_count[i] = 0;
  for (i = 0; i < c_edges->n_elts; i++)
    c_count[i] = 0;

  for (i = 0; i < c_edges->n_elts; i++) {

    vid1 = c_edges->array[2*i] - 1;
    vid2 = c_edges->array[2*i+1] - 1;
    edge_id = _get_join_edge_id(vid1, vid2, edge_builder);

    assert(edge_id != -1);

    c_count[i] =  edge_builder->v2v_sub_idx[edge_id+1]
                - edge_builder->v2v_sub_idx[edge_id];

  } /* End of loop on c_edges */

  /* Exchange number of sub elements for each edge */

  request_count = 0;

  for (i = 0; i < s_edges->n_ranks; i++) {

    int  *recv_buf;

    distant_rank = s_edges->ranks[i];
    length = s_edges->index[i+1] - s_edges->index[i];
    recv_buf = s_count + s_edges->index[i];

    MPI_Irecv(recv_buf,
              length,
              MPI_INT,
              distant_rank,
              distant_rank,
              mpi_comm,
              &(request[request_count++]));

  }

  /* We wait for posting all receives (often recommended) */

  MPI_Barrier(mpi_comm);

  /* Send data to distant ranks */

  for (i = 0; i < c_edges->n_ranks; i++) {

    int  *send_buf;

    distant_rank = c_edges->ranks[i];
    length = c_edges->index[i+1] - c_edges->index[i];
    send_buf = c_count + c_edges->index[i];

    MPI_Isend(send_buf,
              length,
              MPI_INT,
              distant_rank,
              loc_rank,
              mpi_comm,
              &(request[request_count++]));

  }

  /* Wait for all exchanges */

  MPI_Waitall(request_count, request, status);

  /* Define buffer to send whith sub-elements from "coupled" edges to
     "single" edges */

  BFT_MALLOC(c_sub_index, c_edges->n_ranks + 1, int);

  c_sub_index[0] = 0;
  for (i = 0; i < c_edges->n_ranks; i++) {
    c_sub_index[i+1] = 0;
    for (j = c_edges->index[i]; j < c_edges->index[i+1]; j++)
      c_sub_index[i+1] += c_count[j];
  }

  for (i = 0; i < c_edges->n_ranks; i++)
    c_sub_index[i+1] += c_sub_index[i];

  c_sub_size = c_sub_index[c_edges->n_ranks];

  BFT_MALLOC(c_sub_gbuf, c_sub_size, cs_gnum_t);
  BFT_MALLOC(c_sub_coord, 3*c_sub_size, cs_real_t);

  shift = 0;
  for (rank = 0; rank < c_edges->n_ranks; rank++) {

    for (i = c_edges->index[rank]; i < c_edges->index[rank+1]; i++) {

      vid1 = c_edges->array[2*i] - 1;
      vid2 = c_edges->array[2*i+1] - 1;
      edge_id = _get_join_edge_id(vid1, vid2, edge_builder);

      for (j = edge_builder->v2v_sub_idx[edge_id];
           j < edge_builder->v2v_sub_idx[edge_id+1]; j++) {

        vid = edge_builder->v2v_sub_lst[j] - 1;

        c_sub_gbuf[shift] = mesh->global_vtx_num[vid];
        for (k = 0; k < 3; k++)
          c_sub_coord[3*shift+k] = mesh->vtx_coord[3*vid+k];

        shift++;

      }

    } /* End of loop on edges received from rank */

  } /* End of loop on ranks */

  BFT_MALLOC(s_sub_index, s_edges->n_ranks + 1, int);

  s_sub_index[0] = 0;
  for (i = 0; i < s_edges->n_ranks; i++) {
    s_sub_index[i+1] = 0;
    for (j = s_edges->index[i]; j < s_edges->index[i+1]; j++)
      s_sub_index[i+1] += s_count[j];
  }

  for (i = 0; i < s_edges->n_ranks; i++)
    s_sub_index[i+1] += s_sub_index[i];

  s_sub_size = s_sub_index[s_edges->n_ranks];

  BFT_MALLOC(s_sub_gbuf, s_sub_size, cs_gnum_t);
  BFT_MALLOC(s_sub_coord, 3*s_sub_size, cs_real_t);

  /* Exchange sub-edge definition: global vertex number  */

  request_count = 0;

  for (i = 0; i < s_edges->n_ranks; i++) {

    cs_gnum_t *recv_gbuf;

    distant_rank = s_edges->ranks[i];
    length = s_sub_index[i+1] - s_sub_index[i];
    recv_gbuf = s_sub_gbuf + s_sub_index[i];

    MPI_Irecv(recv_gbuf,
              length,
              CS_MPI_GNUM,
              distant_rank,
              distant_rank,
              mpi_comm,
              &(request[request_count++]));

  }

  /* We wait for posting all receives (often recommended) */

  MPI_Barrier(mpi_comm);

  /* Send data to distant ranks */

  for (i = 0; i < c_edges->n_ranks; i++) {

    cs_gnum_t *send_gbuf;

    distant_rank = c_edges->ranks[i];
    length = c_sub_index[i+1] - c_sub_index[i];
    send_gbuf = c_sub_gbuf + c_sub_index[i];

    MPI_Isend(send_gbuf,
              length,
              CS_MPI_GNUM,
              distant_rank,
              loc_rank,
              mpi_comm,
              &(request[request_count++]));

  }

  /* Wait for all exchanges */

  MPI_Waitall(request_count, request, status);

  /* Exchange sub-edge definition: vertex coordinates */

  request_count = 0;

  for (i = 0; i < s_edges->n_ranks; i++) {

    cs_real_t  *recv_dbuf;

    distant_rank = s_edges->ranks[i];
    length = 3*(s_sub_index[i+1] - s_sub_index[i]);
    recv_dbuf = s_sub_coord + 3*s_sub_index[i];

    MPI_Irecv(recv_dbuf,
              length,
              CS_MPI_REAL,
              distant_rank,
              distant_rank,
              mpi_comm,
              &(request[request_count++]));

  }

  /* We wait for posting all receives (often recommended) */

  MPI_Barrier(mpi_comm);

  /* Send data to distant ranks */

  for (i = 0; i < c_edges->n_ranks; i++) {

    cs_real_t  *send_dbuf;

    distant_rank = c_edges->ranks[i];
    length = 3*(c_sub_index[i+1] - c_sub_index[i]);
    send_dbuf = c_sub_coord + 3*c_sub_index[i];

    MPI_Isend(send_dbuf,
              length,
              CS_MPI_REAL,
              distant_rank,
              loc_rank,
              mpi_comm,
              &(request[request_count++]));

  }

  /* Wait for all exchanges */

  MPI_Waitall(request_count, request, status);

  /* Free memory */

  BFT_FREE(c_count);
  BFT_FREE(c_sub_index);
  BFT_FREE(c_sub_gbuf);
  BFT_FREE(c_sub_coord);
  BFT_FREE(s_sub_index);
  BFT_FREE(request);
  BFT_FREE(status);

  if (s_edges->n_elts > 0) {

    cs_lnum_t *order = NULL, *inv_order = NULL;

    /* Update vertices from list of sub elements. Define new vertices */

    BFT_MALLOC(order, s_sub_size, cs_lnum_t);

    cs_order_gnum_allocated(NULL, s_sub_gbuf, order, s_sub_size);

    prev = 0;
    n_new_vertices = 0;

    for (i = 0; i < s_sub_size; i++) {

      cur = s_sub_gbuf[order[i]];
      if (cur != prev) {
        prev = cur;
        id = cs_search_g_binary(mesh->n_vertices,
                                cur,
                                mesh->global_vtx_num);
        if (id == -1) /* Add vertex */
          n_new_vertices++;
      }

    }  /* End of loop on received sub-elements */

    BFT_REALLOC(mesh->global_vtx_num,
                mesh->n_vertices + n_new_vertices, cs_gnum_t);

    BFT_REALLOC(mesh->vtx_coord,
                3*(mesh->n_vertices + n_new_vertices), cs_real_t);

    prev = 0;
    n_new_vertices = 0;

    for (i = 0; i < s_sub_size; i++) {

      cur = s_sub_gbuf[order[i]];
      if (cur != prev) {
        prev = cur;
        id = cs_search_g_binary(mesh->n_vertices,
                                cur,
                                mesh->global_vtx_num);

        if (id == -1) { /* Add vertex */

          shift = mesh->n_vertices + n_new_vertices;
          mesh->global_vtx_num[shift] = cur;
          for (k = 0; k < 3; k++)
            mesh->vtx_coord[3*shift+k] = s_sub_coord[3*order[i]+k];
          n_new_vertices++;

        }

      }

    } /* End of loop on received sub-elements */

    mesh->n_vertices += n_new_vertices;
    if (cs_glob_join_log != NULL)
      fprintf(cs_glob_join_log,
              "  Add %d new vertices from the single elements sync.\n",
              n_new_vertices);

    /* Reorder global_vtx_num in order to have an ordered list */

    BFT_REALLOC(order, mesh->n_vertices, cs_lnum_t);

    cs_order_gnum_allocated(NULL,
                            mesh->global_vtx_num,
                            order,
                            mesh->n_vertices);

    BFT_MALLOC(new_vtx_gnum, mesh->n_vertices, cs_gnum_t);

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

    BFT_MALLOC(inv_order, mesh->n_vertices, cs_lnum_t);

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

    BFT_MALLOC(new_v2v_sub_idx, edge_builder->n_edges + 1, cs_lnum_t);

    for (i = 0; i < edge_builder->n_edges; i++)
      new_v2v_sub_idx[i+1] =  edge_builder->v2v_sub_idx[i+1]
                            - edge_builder->v2v_sub_idx[i];

    for (i = 0; i < s_edges->n_elts; i++) {

      vid1 = s_edges->array[2*i] - 1;
      vid2 = s_edges->array[2*i+1] - 1;
      edge_id = _get_join_edge_id(vid1, vid2, edge_builder);
      assert(edge_id != -1);
      new_v2v_sub_idx[edge_id+1] = s_count[i];

    }

    new_v2v_sub_idx[0] = 0;
    for (i = 0; i < edge_builder->n_edges; i++)
      new_v2v_sub_idx[i+1] += new_v2v_sub_idx[i];

    /* Update v2v_sub_lst */

    BFT_MALLOC(new_v2v_sub_lst,
               new_v2v_sub_idx[edge_builder->n_edges], cs_lnum_t);

    for (i = 0; i < edge_builder->n_edges; i++) {

      shift = new_v2v_sub_idx[i];

      for (j = edge_builder->v2v_sub_idx[i];
           j < edge_builder->v2v_sub_idx[i+1]; j++) {
        vid = edge_builder->v2v_sub_lst[j] - 1;
        new_v2v_sub_lst[shift++] = inv_order[vid] + 1;
      }

    }

    /* Update sub list for single edges */

    shift = 0;

    for (i = 0; i < s_edges->n_elts; i++) {

      vid1 = s_edges->array[2*i] - 1;
      vid2 = s_edges->array[2*i+1] - 1;
      edge_id = _get_join_edge_id(vid1, vid2, edge_builder);
      assert(edge_id != -1);

      for (j = new_v2v_sub_idx[edge_id];
           j < new_v2v_sub_idx[edge_id+1]; j++) {

        id = cs_search_g_binary(mesh->n_vertices,
                                s_sub_gbuf[shift++],
                                mesh->global_vtx_num);

        assert(id != -1);
        new_v2v_sub_lst[j] = id + 1;

      }

    }

    BFT_FREE(order);
    BFT_FREE(inv_order);
    BFT_FREE(edge_builder->v2v_sub_idx);
    BFT_FREE(edge_builder->v2v_sub_lst);

    edge_builder->v2v_sub_idx = new_v2v_sub_idx;
    edge_builder->v2v_sub_lst = new_v2v_sub_lst;

  } /* End if s_edges->n_elts > 0 */

  /* Free memory */

  BFT_FREE(s_count);
  BFT_FREE(s_sub_gbuf);
  BFT_FREE(s_sub_coord);

}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Update elements of a cs_mesh_t structure implied in the joining
 * operation. Update the vertices after the merge step.
 *
 * parameters:
 *   o2n_vtx_gnum       <-- local array on the new global vertex
 *   join_mesh          <-- pointer to a local cs_join_mesh_t struct.
 *   mesh               <-> pointer to cs_mesh_t structure
 *   p_join2mesh_vtx_id --> relation between join mesh and after vertex merge
 *   p_o2n_vtx_id       --> relation between init. and current local num.
 *---------------------------------------------------------------------------*/

static void
_update_vertices_after_merge(const cs_gnum_t        o2n_vtx_gnum[],
                             const cs_join_mesh_t  *join_mesh,
                             cs_mesh_t             *mesh,
                             cs_lnum_t             *p_join2mesh_vtx_id[],
                             cs_lnum_t             *p_o2n_vtx_id[])
{
  cs_lnum_t  i, j, k, o_id, j_id;
  cs_gnum_t  prev, cur;

  cs_lnum_t   n_am_vertices = -1; /* am: after merge */
  cs_real_t  *new_vtx_coord = NULL;
  cs_lnum_t  *o2n_vtx_id = NULL, *join2mesh_vtx_id = NULL;
  cs_lnum_t  *order = NULL;
  cs_gnum_t  *new_vtx_gnum = NULL, *tmp_vtx_gnum = NULL;

  const cs_lnum_t  n_bm_vertices = mesh->n_vertices; /* bm: before merge */
  const cs_lnum_t  n_j_vertices = join_mesh->n_vertices;
  const cs_join_vertex_t  *j_vertices = join_mesh->vertices;
  const cs_lnum_t  n_vertices = n_bm_vertices + n_j_vertices;
  const int  n_ranks = cs_glob_n_ranks;

  /* Update initial vertices (local and global numbering) */

  BFT_MALLOC(o2n_vtx_id, n_bm_vertices, cs_lnum_t);
  BFT_MALLOC(join2mesh_vtx_id, n_j_vertices, cs_lnum_t);
  BFT_MALLOC(tmp_vtx_gnum, n_vertices, cs_gnum_t);
  BFT_MALLOC(new_vtx_gnum, n_vertices, cs_gnum_t);
  BFT_MALLOC(order, n_vertices, cs_lnum_t);

  for (i = 0; i < n_bm_vertices; i++)
    tmp_vtx_gnum[i] = o2n_vtx_gnum[i];

  for (i = 0, j = n_bm_vertices; i < n_j_vertices; i++, j++)
    tmp_vtx_gnum[j] = j_vertices[i].gnum;

  cs_order_gnum_allocated(NULL, tmp_vtx_gnum, order, n_vertices);

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
  if (cs_glob_join_log != NULL) {
    fprintf(cs_glob_join_log,
            "\n\n Dump Old2New array (local mesh): "
            "old_n_vertices = %d - new_n_vertices = %d\n",
            n_bm_vertices, n_am_vertices);
    for (i = 0; i < n_bm_vertices; i++)
      fprintf(cs_glob_join_log,
              "Old num : %7d (%9llu) => New num : %7d (%9llu) (%9llu)\n",
              i+1,
              (unsigned long long)(n_ranks >1 ?   mesh->global_vtx_num[i]
                                                : (cs_gnum_t)i+1),
              o2n_vtx_id[i]+1, (unsigned long long)o2n_vtx_gnum[i],
              (unsigned long long)new_vtx_gnum[o2n_vtx_id[i]]);
    fflush(cs_glob_join_log);
  }
#endif

  /* Update global vertex information */

  mesh->n_vertices = n_am_vertices;
  mesh->n_g_vertices =  new_vtx_gnum[n_am_vertices - 1];

  BFT_REALLOC(new_vtx_gnum, n_am_vertices, cs_gnum_t);
  BFT_FREE(mesh->global_vtx_num);
  mesh->global_vtx_num = new_vtx_gnum;

#if defined(HAVE_MPI)
  if (n_ranks > 1) { /* Get temporary a no contiguous numbering */

    cs_gnum_t glob_gmax;
    cs_gnum_t loc_gmax = new_vtx_gnum[n_am_vertices - 1];
    MPI_Comm  mpi_comm = cs_glob_mpi_comm;

    /* Find the max. global number */

    MPI_Allreduce(&loc_gmax, &glob_gmax, 1, CS_MPI_GNUM, MPI_MAX, mpi_comm);

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
 *   join_mesh          <-- pointer to a local cs_join_mesh_t struct.
 *   mesh               <-> pointer of pointer to cs_mesh_t structure
 *   p_join2mesh_vtx_id --> relation between join mesh and after merge vertex
 *---------------------------------------------------------------------------*/

static void
_update_vertices_after_split(const cs_join_mesh_t  *join_mesh,
                             cs_mesh_t             *mesh,
                             cs_lnum_t             *p_join2mesh_vtx_id[])
{
  cs_lnum_t  i, j, k, o_id, j_id, v_id;
  cs_gnum_t  prev, cur;

  bool tmp_global_vtx_num = false;

  cs_lnum_t  n_as_vertices = -1; /* ac: after splitting */
  cs_lnum_t  *order = NULL;
  cs_real_t  *new_vtx_coord = NULL;
  cs_lnum_t  *o2n_vtx_id = NULL, *join2mesh_vtx_id = NULL;
  cs_gnum_t  *new_vtx_gnum = NULL, *tmp_vtx_gnum = NULL;

  const cs_lnum_t  n_bs_vertices = mesh->n_vertices; /* bs: before splitting */
  const cs_lnum_t  n_j_vertices = join_mesh->n_vertices;
  const cs_join_vertex_t  *j_vertices = join_mesh->vertices;
  const cs_lnum_t  n_vertices = n_bs_vertices + n_j_vertices;

  if (mesh->global_vtx_num == NULL) {
    tmp_global_vtx_num = true;
    BFT_MALLOC(mesh->global_vtx_num, n_bs_vertices, cs_gnum_t);
    for (i = 0; i < n_bs_vertices; i++)
      mesh->global_vtx_num[i] = i+1;
  }

  /* Update initial vertices (local and global numbering) */

  BFT_MALLOC(o2n_vtx_id, n_bs_vertices, cs_lnum_t);
  BFT_MALLOC(join2mesh_vtx_id, n_j_vertices, cs_lnum_t);
  BFT_MALLOC(tmp_vtx_gnum, n_vertices, cs_gnum_t);
  BFT_MALLOC(new_vtx_gnum, n_vertices, cs_gnum_t);
  BFT_MALLOC(order, n_vertices, cs_lnum_t);

  for (i = 0; i < n_bs_vertices; i++)
    tmp_vtx_gnum[i] = mesh->global_vtx_num[i];

  for (i = 0, j = n_bs_vertices; i < n_j_vertices; i++, j++)
    tmp_vtx_gnum[j] = j_vertices[i].gnum;

  cs_order_gnum_allocated(NULL, tmp_vtx_gnum, order, n_vertices);

  /* Define o2n_vtx_id and join2mesh_vtx_id arrays */

  assert(n_vertices > 0);

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
  BFT_REALLOC(new_vtx_gnum, n_as_vertices, cs_gnum_t);
  BFT_MALLOC(new_vtx_coord, 3*n_as_vertices, cs_real_t);

  mesh->n_vertices = n_as_vertices;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (cs_glob_join_log != NULL) {
    fprintf(cs_glob_join_log,
            "\n\n Dump Old2New array (local mesh): "
            "old_n_vertices = %d - new_n_vertices = %d\n",
            n_bs_vertices, n_as_vertices);
    for (i = 0; i < n_bs_vertices; i++)
      fprintf(cs_glob_join_log,
              "Old num : %7d (%9llu) => New num : %7d (%9llu)\n",
              i+1,
              (unsigned long long)(cs_glob_n_ranks >1 ?  mesh->global_vtx_num[i]
                                                        : (cs_gnum_t)i+1),
              o2n_vtx_id[i]+1, (unsigned long long)new_vtx_gnum[o2n_vtx_id[i]]);
    fflush(cs_glob_join_log);
  }
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
    for (j = mesh->i_face_vtx_idx[i]; j < mesh->i_face_vtx_idx[i+1]; j++) {
      v_id = mesh->i_face_vtx_lst[j];
      mesh->i_face_vtx_lst[j] = o2n_vtx_id[v_id];
    }
  }

  /* Update border face connect. */

  for (i = 0; i < mesh->n_b_faces; i++) {
    for (j = mesh->b_face_vtx_idx[i]; j < mesh->b_face_vtx_idx[i+1]; j++) {
      v_id = mesh->b_face_vtx_lst[j];
      mesh->b_face_vtx_lst[j] = o2n_vtx_id[v_id];
    }
  }

  /* Free memory */

  if (tmp_global_vtx_num)
    BFT_FREE(mesh->global_vtx_num);

  BFT_FREE(o2n_vtx_id);
  BFT_FREE(mesh->vtx_coord);
  BFT_FREE(mesh->global_vtx_num);

  mesh->vtx_coord = new_vtx_coord;
  mesh->global_vtx_num = new_vtx_gnum;

  /* Return pointer */

  *p_join2mesh_vtx_id = join2mesh_vtx_id;

}

/*----------------------------------------------------------------------------
 * Update border and interior face state
 *
 * parameters:
 *  selection    <->  list of entities participating in the join operation
 *  join_mesh    <--  pointer to a cs_join_mesh_t structure
 *  mesh         <--  pointer of pointer to cs_mesh_t structure
 *  j2m_vid      <--  relation between join mesh and after merge vertex
 *---------------------------------------------------------------------------*/

static void
_update_face_state(cs_join_select_t        *selection,
                   const cs_join_mesh_t    *join_mesh,
                   const cs_mesh_t         *mesh,
                   const cs_lnum_t          j2m_vid[])
{
  int  i, j, s, e;
  cs_join_state_t  v_state, f_state;
  bool  have_new;

  assert(selection != NULL);

  cs_join_state_t  *states = NULL;

  BFT_MALLOC(states, mesh->n_vertices, cs_join_state_t);

  /* Initialize */

  for (i = 0; i < mesh->n_vertices; i++)
    states[i] = CS_JOIN_STATE_UNDEF;

  /* Update vertex state thanks to join_mesh */

  for (i = 0; i < join_mesh->n_vertices; i++)
    states[j2m_vid[i]] = join_mesh->vertices[i].state;

  /* Border faces treatment */

  for (i = 0; i < mesh->n_b_faces; i++) {

    s = mesh->b_face_vtx_idx[i];
    e = mesh->b_face_vtx_idx[i+1];
    f_state = CS_JOIN_STATE_UNDEF;
    have_new = false;

    for (j = s; j < e; j++) {

      v_state = states[mesh->b_face_vtx_lst[j]];
      f_state = CS_MAX(f_state, v_state);
      if (v_state == CS_JOIN_STATE_NEW)
        have_new = true;

    }

    if (f_state == CS_JOIN_STATE_ORIGIN && have_new == true)
      f_state = CS_JOIN_STATE_MERGE;
    else if (f_state == CS_JOIN_STATE_PERIO && have_new == true)
      f_state = CS_JOIN_STATE_PERIO_MERGE;

    selection->b_face_state[i] = f_state;

  } /* End of loop on border faces */

  /* Interior faces treatment */

  for (i = 0; i < mesh->n_i_faces; i++) {

    s = mesh->i_face_vtx_idx[i];
    e = mesh->i_face_vtx_idx[i+1];
    f_state = CS_JOIN_STATE_UNDEF;
    have_new = false;

    for (j = s; j < e; j++) {

      v_state = states[mesh->i_face_vtx_lst[j]];
      f_state = CS_MAX(f_state, v_state);
      if (v_state == CS_JOIN_STATE_NEW)
        have_new = true;

    }

    if (f_state == CS_JOIN_STATE_ORIGIN && have_new == true)
      f_state = CS_JOIN_STATE_MERGE;
    else if (f_state == CS_JOIN_STATE_PERIO && have_new == true)
      f_state = CS_JOIN_STATE_PERIO_MERGE;

    selection->i_face_state[i] = f_state;

  } /* End of loop on interior faces */

  /* Free memory */

  BFT_FREE(states);

}

/*----------------------------------------------------------------------------
 * Allocate and initialize the definition of a edge_builder_t structure.
 * Equivalent to build a local edge-based connectivity for the selected faces.
 *
 * parameters:
 *   join_select <-- list of entities participating in the join operation
 *   mesh        <-- pointer to a cs_mesh_t structure
 *
 * returns:
 *   a new allocated pointer to an edge_builder_t structure.
 *---------------------------------------------------------------------------*/

static edge_builder_t *
_init_edge_builder(const cs_join_select_t  *join_select,
                   const cs_mesh_t         *mesh)
{
  cs_lnum_t  i, j, save, shift;

  cs_lnum_t  *count = NULL;
  edge_builder_t  *edge_builder = NULL;

  assert(sizeof(cs_lnum_t) == sizeof(cs_lnum_t));

  /* Allocate and initialize edge_builder_t structure */

  BFT_MALLOC(edge_builder, 1, edge_builder_t);

  edge_builder->n_vertices = mesh->n_vertices;
  edge_builder->n_edges = 0;
  edge_builder->v2v_lst = NULL;
  edge_builder->v2v_sub_idx = NULL;
  edge_builder->v2v_sub_lst = NULL;

  BFT_MALLOC(edge_builder->v2v_idx, edge_builder->n_vertices + 1, cs_lnum_t);

  for (i = 0; i < edge_builder->n_vertices + 1; i++)
    edge_builder->v2v_idx[i] = 0;

  /* Build vertex -> vertex index */

  cs_join_build_edges_idx(join_select->n_faces,
                          join_select->faces,
                          mesh->b_face_vtx_idx,
                          mesh->b_face_vtx_lst,
                          edge_builder->v2v_idx);

  if (cs_glob_n_ranks > 1) {   /* Add edges from "single" edges */

    cs_join_sync_t  *s_edges = join_select->s_edges;

    assert(s_edges != NULL);

    if (s_edges->n_elts > 0)
      for (i = 0; i < s_edges->n_elts; i++)
        edge_builder->v2v_idx[s_edges->array[2*i]] += 1;

  }

  BFT_MALLOC(count, edge_builder->n_vertices, cs_lnum_t);

  for (i = 0; i < edge_builder->n_vertices; i++) {
    edge_builder->v2v_idx[i+1] += edge_builder->v2v_idx[i];
    count[i] = 0;
  }

  edge_builder->n_edges = edge_builder->v2v_idx[edge_builder->n_vertices];

  /* Build vertex -> vertex list */

  BFT_MALLOC(edge_builder->v2v_lst, edge_builder->n_edges, cs_lnum_t);

  /* Fill v2v_lst */

  cs_join_build_edges_lst(join_select->n_faces,
                          join_select->faces,
                          mesh->b_face_vtx_idx,
                          mesh->b_face_vtx_lst,
                          count,
                          edge_builder->v2v_idx,
                          edge_builder->v2v_lst);

  /* Add edges from "single" edges */

  if (cs_glob_n_ranks > 1) {   /* Add edges from "single" edges */

    cs_join_sync_t  *s_edges = join_select->s_edges;

    for (i = 0; i < s_edges->n_elts; i++) {

      int  vid1 = s_edges->array[2*i] - 1;
      int  vid2 = s_edges->array[2*i+1] - 1;

      assert(vid1 < vid2);
      shift = edge_builder->v2v_idx[vid1] + count[vid1];
      edge_builder->v2v_lst[shift] = vid2 + 1;
      count[vid1] += 1;

    }

  }

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

    cs_lnum_t  start = save;
    cs_lnum_t  end = edge_builder->v2v_idx[i+1];

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
  BFT_REALLOC(edge_builder->v2v_lst, edge_builder->n_edges, cs_lnum_t);

  return edge_builder;
}

/*----------------------------------------------------------------------------
 * Get the local face connectivity for the current selected face before/after
 * the merge step.
 *
 * parameters:
 *   select_id    <-- id of the selected face in selection array
 *   o2n_vtx_gnum <-- local array on the new global vertex
 *   join_select  <-- list of all implied entities in the joining operation
 *   join_mesh    <-- pointer to the local cs_join_mesh_t structure
 *   mesh         <-- pointer to cs_mesh_t structure
 *   bm_tmp       <-> connectivity before the merge step
 *   am_tmp       <-> connectivity after the merge step
 *---------------------------------------------------------------------------*/

static void
_get_local_faces_connect(cs_lnum_t                select_id,
                         const cs_gnum_t          o2n_vtx_gnum[],
                         const cs_join_select_t  *join_select,
                         const cs_join_mesh_t    *join_mesh,
                         const cs_mesh_t         *mesh,
                         cs_lnum_t                bm_tmp[],
                         cs_lnum_t                am_tmp[])
{
  cs_lnum_t  i, j, k, v_id, bm_shift;
  cs_gnum_t  new_gnum;

  cs_lnum_t  fid = join_select->faces[select_id] - 1;
  cs_gnum_t  fgnum = join_select->compact_face_gnum[select_id];
  cs_lnum_t  am_s = join_mesh->face_vtx_idx[select_id];
  cs_lnum_t  am_e = join_mesh->face_vtx_idx[select_id+1];
  cs_lnum_t  n_am_face_vertices = am_e - am_s;
  cs_lnum_t  bm_s = mesh->b_face_vtx_idx[fid];
  cs_lnum_t  bm_e = mesh->b_face_vtx_idx[fid+1];
  cs_lnum_t  n_bm_face_vertices = bm_e - bm_s;
  cs_lnum_t  fst_match_id = -1;

  const cs_join_vertex_t  *vertices = join_mesh->vertices;

  assert(join_mesh->face_gnum[select_id] == fgnum);

  /* Store the face connectivity before the merge step */

  for (i = bm_s, j = 0; i < bm_e; i++, j++)
    bm_tmp[j] = mesh->b_face_vtx_lst[i];
  bm_tmp[n_bm_face_vertices] = mesh->b_face_vtx_lst[bm_s];

  /* Store the face connectivity after the merge step */

  for (i = am_s, j = 0; i < am_e; i++, j++)
    am_tmp[j] = join_mesh->face_vtx_lst[i];
  am_tmp[n_am_face_vertices] = join_mesh->face_vtx_lst[am_s];

  /* Find position of initial vertices in the face connectivity after the
     merge step. If not found -1 (initialized value) */

  bm_shift = 0;
  while (fst_match_id == -1 && bm_shift < n_bm_face_vertices) {

    v_id = bm_tmp[bm_shift];
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
                " connectivity before/after the merge step.\n"
                "  Current global face number: %llu\n"),
              (unsigned long long)fgnum);

  /* Store the face connectivity before the merge step */

  for (i = 0; i < n_bm_face_vertices; i++) {
    j = bm_s + (bm_shift + i) % n_bm_face_vertices;
    bm_tmp[i] = mesh->b_face_vtx_lst[j];
  }
  bm_tmp[n_bm_face_vertices] = mesh->b_face_vtx_lst[bm_s + bm_shift];

  /* Store the face connectivity after the merge step */

  for (i = 0; i < n_am_face_vertices; i++) {
    j = am_s + (fst_match_id + i) % n_am_face_vertices;
    am_tmp[i] = join_mesh->face_vtx_lst[j];
  }
  am_tmp[n_am_face_vertices] = join_mesh->face_vtx_lst[am_s + fst_match_id];

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (cs_glob_join_log != NULL) {
    fprintf(cs_glob_join_log, "\n\n  Face: %d (%llu)\n", fid+1,
            (unsigned long long)fgnum);
    fprintf(cs_glob_join_log, "bm_tmp:%d, n: %2d, v:", i, n_bm_face_vertices+1);
    for (j = 0; j < n_bm_face_vertices + 1; j++)
      fprintf(cs_glob_join_log, " %8d (%llu)",
             bm_tmp[j]+1, (unsigned long long)o2n_vtx_gnum[bm_tmp[j]]);
    fprintf(cs_glob_join_log, "\nam_tmp:%d, n: %2d, v:",
           i, n_am_face_vertices+1);
    for (j = 0; j < n_am_face_vertices + 1; j++)
      fprintf(cs_glob_join_log, " %8d (%llu)",
             am_tmp[j]+1, (unsigned long long)vertices[am_tmp[j]].gnum);
    fprintf(cs_glob_join_log, "\n");
    fflush(cs_glob_join_log);
  }
#endif

}

/*----------------------------------------------------------------------------
 * Define the new edge connectivity for the initial edges.
 *
 * parameters:
 *   join_select      <-- keep all implied entities in the joining operation
 *   join_mesh        <-- pointer to the local cs_join_mesh_t structure
 *   mesh             <-- pointer of pointer to cs_mesh_t structure
 *   o2n_vtx_gnum     <-- local array on the new global vertex
 *   join2mesh_vtx_id <-- relation between join mesh and vertex after merge
 *   edge_builder     <-> pointer to an edge_builder_t structure
 *---------------------------------------------------------------------------*/

static void
_complete_edge_builder(const cs_join_select_t  *join_select,
                       const cs_join_mesh_t    *join_mesh,
                       const cs_mesh_t         *mesh,
                       const cs_gnum_t          o2n_vtx_gnum[],
                       const cs_lnum_t          join2mesh_vtx_id[],
                       edge_builder_t          *edge_builder)
{
  cs_lnum_t  i, j, j1, j2, k, shift;
  cs_lnum_t  v1_id, v2_id, edge_id, n_subs;
  cs_gnum_t  v1_gnum, v2_gnum;
  bool  direct_scan, degenerate_edge;

  cs_lnum_t  am_max = 0, bm_max = 0;
  cs_lnum_t  *am_tmp = NULL, *bm_tmp = NULL;

  const cs_join_vertex_t  *vertices = join_mesh->vertices;

  /* Sanity checks */

  assert(edge_builder != NULL);
  assert(join_mesh->n_faces == join_select->n_faces);

  /* Define a list of new vertices for each initial selected edge */

  BFT_MALLOC(edge_builder->v2v_sub_idx, edge_builder->n_edges + 1, cs_lnum_t);

  for (i = 0; i < edge_builder->n_edges + 1; i++)
    edge_builder->v2v_sub_idx[i] = 0;

  for (i = 0; i < join_mesh->n_faces; i++)
    am_max = CS_MAX(am_max,
                    join_mesh->face_vtx_idx[i+1]-join_mesh->face_vtx_idx[i]);

  for (i = 0; i < join_select->n_faces; i++) {
    j = join_select->faces[i] - 1;
    bm_max = CS_MAX(bm_max, mesh->b_face_vtx_idx[j+1]-mesh->b_face_vtx_idx[j]);
  }

  BFT_MALLOC(am_tmp, am_max + 1, cs_lnum_t);
  BFT_MALLOC(bm_tmp, bm_max + 1, cs_lnum_t);

  /* Count the number of sub-elements to add to each initial edge */

  for (i = 0; i < join_select->n_faces; i++) {

    cs_lnum_t  fid = join_select->faces[i] - 1;
    cs_gnum_t  fgnum = join_select->compact_face_gnum[i];
    cs_lnum_t  bm_s = mesh->b_face_vtx_idx[fid];
    cs_lnum_t  bm_e = mesh->b_face_vtx_idx[fid+1];
    cs_lnum_t  n_bm_face_vertices = bm_e - bm_s;
    cs_lnum_t  am_s = join_mesh->face_vtx_idx[i];
    cs_lnum_t  am_e = join_mesh->face_vtx_idx[i+1];
    cs_lnum_t  n_am_face_vertices = am_e - am_s;

    _get_local_faces_connect(i,           /* id of the face */
                             o2n_vtx_gnum,
                             join_select,
                             join_mesh,
                             mesh,
                             bm_tmp,
                             am_tmp);

    for (j = 0, j1 = 0, j2 = 0; j < n_bm_face_vertices; j++) {

      degenerate_edge = false;
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
          degenerate_edge = true;
        }
        else { /* Look for the next initial vertex */

          for (j2 = j1 + 1; j2 < n_am_face_vertices + 1; j2++)
            if (v2_gnum == vertices[am_tmp[j2]].gnum)
              break;

          if (j2 == n_am_face_vertices + 1) { /* Init. edge has been deleted */

            j2 = j1 + 1;
            degenerate_edge = true;
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
               a degenerate edge */

            if (degenerate_edge == false && n_subs == 1)
              bft_error(__FILE__, __LINE__, 0,
                        _("  Face %d (%llu): Edge with two different"
                          " descriptions: (%d, %d) [%llu, %llu]\n"
                          "  n_subs: %d - previous n_subs: %d\n"
                          "  Impossible to continue the mesh update after the"
                          " merge operation.\n"),
                        fid+1, (unsigned long long)fgnum, v1_id+1, v2_id+1,
                        (unsigned long long)v1_gnum,
                        (unsigned long long)v2_gnum,
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
             edge_builder->v2v_sub_idx[edge_builder->n_edges], cs_lnum_t);

  for (i = 0; i < edge_builder->v2v_sub_idx[edge_builder->n_edges]; i++)
    edge_builder->v2v_sub_lst[i] = -1; /* value = degenerate edge */

  /* Add sub-elements to each initial edge */

  for (i = 0; i < join_select->n_faces; i++) {

    cs_lnum_t  fid = join_select->faces[i] - 1;
    cs_lnum_t  bm_s = mesh->b_face_vtx_idx[fid];
    cs_lnum_t  bm_e = mesh->b_face_vtx_idx[fid+1];
    cs_lnum_t  n_bm_face_vertices = bm_e - bm_s;
    cs_lnum_t  am_s = join_mesh->face_vtx_idx[i];
    cs_lnum_t  am_e = join_mesh->face_vtx_idx[i+1];
    cs_lnum_t  n_am_face_vertices = am_e - am_s;

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

        if (v1_gnum != v2_gnum) { /* Initial edge is not degenerate */

          for (j2 = j1 + 1; j2 < n_am_face_vertices + 1; j2++)
            if (v2_gnum == vertices[am_tmp[j2]].gnum)
              break;

          if (j2 == n_am_face_vertices + 1) /* Initial edge is degenerate */
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

        } /* Initial edge is not degenerate */

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
 *   join_select      <-- cs_join_select_t struct.
 *   join_mesh        <-- cs_join_mesh_t structure
 *   join2mesh_vtx_id <-- relation between vtx_id of join_mesh & cs_mesh_t
 *   n_faces          <-- local number of faces in the mesh
 *   p_f2v_idx        <-> face -> vertex connectivity index
 *   p_f2v_lst        <-> face -> vertex connectivity list
 *---------------------------------------------------------------------------*/

static void
_update_selected_face_connect(const cs_join_select_t  *join_select,
                              const cs_join_mesh_t    *join_mesh,
                              const cs_lnum_t          join2mesh_vtx_id[],
                              cs_lnum_t                n_faces,
                              cs_lnum_t               *p_f2v_idx[],
                              cs_lnum_t               *p_f2v_lst[])
{
  cs_lnum_t  i, j, shift, v_id, select_id, join_fid;
  cs_gnum_t  fgnum;

  cs_lnum_t  *new_f2v_lst = NULL, *new_f2v_idx = NULL;
  cs_lnum_t  *f2v_idx = *p_f2v_idx;
  cs_lnum_t  *f2v_lst = *p_f2v_lst;

  BFT_MALLOC(new_f2v_idx, n_faces + 1, cs_lnum_t);

  for (i = 0; i < n_faces + 1; i++)
    new_f2v_idx[i] = 0;

  for (i = 0, select_id = 0; i < n_faces; i++) {

    bool  in_selection = false;

    if (select_id < join_select->n_faces)
      if (i+1 == join_select->faces[select_id])
        in_selection = true;

    if (in_selection) { /* Among selected faces */

      fgnum = join_select->compact_face_gnum[select_id];

      if (join_mesh->face_gnum[select_id] == fgnum)
        join_fid = select_id;
      else
        join_fid = cs_search_g_binary(join_mesh->n_faces,
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

  new_f2v_idx[0] = 0;
  for (i = 0; i < n_faces; i++)
    new_f2v_idx[i+1] += new_f2v_idx[i];

  BFT_MALLOC(new_f2v_lst, new_f2v_idx[n_faces], cs_lnum_t);

  for (i = 0, select_id = 0; i < n_faces; i++) {

    bool  in_selection = false;

    if (select_id < join_select->n_faces)
      if (i+1 == join_select->faces[select_id])
        in_selection = true;

    shift = new_f2v_idx[i];

    if (in_selection) { /* Among selected faces */

      fgnum = join_select->compact_face_gnum[select_id];

      if (join_mesh->face_gnum[select_id] == fgnum)
        join_fid = select_id;
      else
        join_fid = cs_search_g_binary(join_mesh->n_faces,
                                      fgnum,
                                      join_mesh->face_gnum);

      for (j = join_mesh->face_vtx_idx[join_fid];
           j < join_mesh->face_vtx_idx[join_fid+1]; j++, shift++) {
        v_id = join_mesh->face_vtx_lst[j];
        new_f2v_lst[shift] = join2mesh_vtx_id[v_id];
      }

      select_id++;

    }
    else { /* Not a selected face. Do not update. */

      for (j = f2v_idx[i]; j < f2v_idx[i+1]; j++, shift++)
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
 *   n_adj_faces  <-- number of adjacent faces
 *   adj_faces    <-- list of adjacent face numbers
 *   edge_builder <-- edge_builder_t structure
 *   o2n_vtx_id   <-- relation between old/new vertix id
 *   n_faces      <-- local number of faces in the mesh
 *   p_f2v_idx    <-- face -> vertex connectivity index
 *   p_f2v_lst    <-- face -> vertex connectivity list
 *---------------------------------------------------------------------------*/

static void
_update_adj_face_connect(cs_lnum_t              n_adj_faces,
                         const cs_lnum_t        adj_faces[],
                         const edge_builder_t  *edge_builder,
                         const cs_lnum_t        o2n_vtx_id[],
                         cs_lnum_t              n_faces,
                         cs_lnum_t             *p_f2v_idx[],
                         cs_lnum_t             *p_f2v_lst[])
{
  cs_lnum_t  i, j, k, l, v1_id, v2_id, n_face_vertices, shift, select_id;
  cs_lnum_t  s, e, v_s, v_e, v_sub_s, v_sub_e, edge_id;

  cs_lnum_t  max = 0;
  cs_lnum_t  *new_f2v_lst = NULL, *new_f2v_idx = NULL, *tmp = NULL;
  cs_lnum_t  *f2v_idx = *p_f2v_idx;
  cs_lnum_t  *f2v_lst = *p_f2v_lst;

  BFT_MALLOC(new_f2v_idx, n_faces+1, cs_lnum_t);

  for (i = 0; i < n_faces + 1; i++)
    new_f2v_idx[i] = 0;

  for (i = 0; i < n_faces; i++)
    max = CS_MAX(max, f2v_idx[i+1] - f2v_idx[i]);

  BFT_MALLOC(tmp, max + 1, cs_lnum_t);

  /* first: update index (counting phase) */

  for (i = 0, select_id = 0; i < n_faces; i++) {

    bool  in_selection = false;

    if (select_id < n_adj_faces)
      if (i+1 == adj_faces[select_id])
        in_selection = true;

    if (in_selection) { /* Among selected faces */

      cs_lnum_t  count = 0;

      s = f2v_idx[i];
      e = f2v_idx[i+1];
      n_face_vertices = e - s;

      for (k = 0; k < n_face_vertices; k++)
        tmp[k] = f2v_lst[s + k];
      tmp[n_face_vertices] = f2v_lst[s];

      for (j = 0; j < n_face_vertices; j++) { /* Scan edges */

        if (tmp[j] < tmp[j+1]) {
          v1_id = tmp[j];
          v2_id = tmp[j+1];
        }
        else if (tmp[j+1] < tmp[j]) {
          v1_id = tmp[j+1];
          v2_id = tmp[j];
        }
        else { /* delete the current edge (count += 0) */
          v1_id = -1;
          v2_id = -1;
        }

        /* edge builder->n_vertices is still equal to the initial
           number of vertices; this should not be a problem,
           as an edge always starts with a previously existing
           vertex */

        assert(v1_id < edge_builder->n_vertices);

        if (v1_id > -1 && v1_id < edge_builder->n_vertices) {

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

  new_f2v_idx[0] = 0;
  for (i = 0; i < n_faces; i++)
    new_f2v_idx[i+1] += new_f2v_idx[i];

  BFT_MALLOC(new_f2v_lst, new_f2v_idx[n_faces], cs_lnum_t);

  /* second: update list (filling phase) */

  for (i = 0, select_id = 0; i < n_faces; i++) {

    bool  direct_scan;
    bool  in_selection = false;

    if (select_id < n_adj_faces)
      if (i+1 == adj_faces[select_id])
        in_selection = true;

    shift = new_f2v_idx[i];

    if (in_selection) { /* Among selected faces */

      cs_lnum_t  count = 0;

      s = f2v_idx[i];
      e = f2v_idx[i+1];
      n_face_vertices = e - s;

      for (k = 0; k < n_face_vertices; k++)
        tmp[k] = f2v_lst[s + k];
      tmp[n_face_vertices] = f2v_lst[s];

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

        if (v1_id > -1 && v1_id < edge_builder->n_vertices) {

          v_s = edge_builder->v2v_idx[v1_id];
          v_e = edge_builder->v2v_idx[v1_id+1];

          new_f2v_lst[shift + count] = o2n_vtx_id[tmp[j]];
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
                    new_f2v_lst[shift + count] = edge_builder->v2v_sub_lst[l] - 1;
                else
                  for (l = v_sub_e - 1; l >= v_sub_s; l--, count++)
                    new_f2v_lst[shift + count] = edge_builder->v2v_sub_lst[l] - 1;

              } /* edge is not degenerate */

            } /* There are sub-elements to add */

          } /* End if exist edge_id */

        } /* End if v1_id > -1 */

      } /* End of loop on face connectivity */

      select_id++;

    }
    else  /* Not in selection. Do not update */
      for (j = f2v_idx[i]; j < f2v_idx[i+1]; j++)
        new_f2v_lst[shift++] = f2v_lst[j];

  } /* End of loop on faces */

  BFT_FREE(f2v_lst);
  BFT_FREE(f2v_idx);
  BFT_FREE(tmp);

  /* Return pointers */

  *p_f2v_idx = new_f2v_idx;
  *p_f2v_lst = new_f2v_lst;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Get the related global cell number for each old global face number
 *
 * parameters:
 *   n2o_hist     <--  new -> old global face numbering
 *   join_select  <--  list all local entities implied in the joining op.
 *   mesh         <--  pointer to the associated cs_mesh_t structure
 *   join_param   <--  set of user-defined parameter
 *   cell_gnum    <->  global cell number related to each old face
 *   face_family  <-> family number related to each old face
 *---------------------------------------------------------------------------*/

static void
_exchange_cell_gnum_and_family(const cs_join_gset_t     *n2o_hist,
                               const cs_join_select_t   *join_select,
                               const cs_mesh_t          *mesh,
                               cs_join_param_t           join_param,
                               cs_gnum_t                 cell_gnum[],
                               cs_lnum_t                 face_family[])
{
  int  rank;
  cs_lnum_t  i, j, fid, shift;
  cs_gnum_t  compact_fgnum;

  int  reduce_size = 0;
  int  *reduce_ids = NULL, *parent = NULL;
  int  *send_count = NULL, *recv_count = NULL;
  int  *send_shift = NULL, *recv_shift = NULL;
  cs_gnum_t  *recv_gbuf = NULL, *send_gbuf = NULL, *reduce_index = NULL;

  MPI_Comm  mpi_comm = cs_glob_mpi_comm;

  const int  n_ranks = cs_glob_n_ranks;
  const int  loc_rank = CS_MAX(cs_glob_rank_id, 0);
  const cs_gnum_t  *gnum_rank_index = join_select->compact_rank_index;
  const cs_gnum_t  loc_rank_s = join_select->compact_rank_index[loc_rank];
  const cs_gnum_t  loc_rank_e = join_select->compact_rank_index[loc_rank+1];

  /* Sanity checks */

  assert(gnum_rank_index != NULL);
  assert(n2o_hist != NULL);
  assert(n_ranks > 1);

  BFT_MALLOC(send_count, n_ranks, int);
  BFT_MALLOC(send_shift, n_ranks + 1, int);

  send_shift[0] = 0;
  for (i = 0; i < n_ranks; i++) {
    send_count[i] = 0;
    send_shift[i+1] = 0;
  }

  /* Compact init. global face distribution. Remove ranks without face
     at the begining */

  for (i = 0; i < n_ranks; i++)
    if (gnum_rank_index[i] < gnum_rank_index[i+1])
      reduce_size++;

  BFT_MALLOC(reduce_index, reduce_size+1, cs_gnum_t);
  BFT_MALLOC(reduce_ids, reduce_size, int);

  reduce_size = 0;

  /* Add +1 to gnum_rank_index because it's an id and we work on numbers */

  reduce_index[0] = gnum_rank_index[0] + 1;
  for (i = 0; i < n_ranks; i++) {
    if (gnum_rank_index[i] < gnum_rank_index[i+1]) {
      reduce_index[reduce_size+1] = gnum_rank_index[i+1] + 1;
      reduce_ids[reduce_size++] = i;
    }
  }

  /* Count number of ranks associated to each new face */

  for (i = 0; i < n2o_hist->n_elts; i++) {

    for (j = n2o_hist->index[i]; j < n2o_hist->index[i+1]; j++) {

      int  reduce_rank = cs_search_gindex_binary(reduce_size,
                                                 n2o_hist->g_list[j],
                                                 reduce_index);

      assert(reduce_rank < reduce_size);
      assert(reduce_rank != -1);

      rank = reduce_ids[reduce_rank];
      send_shift[rank+1] += 1;

    } /* End of loop on old faces */

  } /* End of loop on new faces */

  for (i = 0; i < n_ranks; i++)
    send_shift[i+1] += send_shift[i];

  assert(send_shift[n_ranks] == n2o_hist->index[n2o_hist->n_elts]);

  BFT_MALLOC(send_gbuf, send_shift[n_ranks]*2, cs_gnum_t);
  BFT_MALLOC(parent, send_shift[n_ranks], int);

  /* Fill the list of ranks */

  for (i = 0; i < n2o_hist->n_elts; i++) {

    for (j = n2o_hist->index[i]; j < n2o_hist->index[i+1]; j++) {

      int  reduce_rank = cs_search_gindex_binary(reduce_size,
                                                 n2o_hist->g_list[j],
                                                 reduce_index);

      assert(reduce_rank < reduce_size);
      assert(reduce_rank != -1);

      rank = reduce_ids[reduce_rank];
      shift = send_shift[rank] + send_count[rank];
      parent[shift] = j;
      send_gbuf[shift] = n2o_hist->g_list[j];
      send_count[rank] += 1;

    } /* End of loop on old faces */

  } /* End of loop on new faces */

  /* Free memory */

  BFT_FREE(reduce_ids);
  BFT_FREE(reduce_index);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (cs_glob_join_log != NULL) {
    fprintf(cs_glob_join_log,
           "\n Exchange to update mesh after the face split operation:\n");
    for (i = 0; i < n_ranks; i++) {
      int start = send_shift[i], end = send_shift[i+1];
      fprintf(cs_glob_join_log, " Send to rank %5d (n = %10d):", i, end - start);
      for (j = start; j < end; j++)
        fprintf(cs_glob_join_log, " %llu ", (unsigned long long)send_gbuf[j]);
      fprintf(cs_glob_join_log, "\n");
    }
  }
#endif

  /* Count the number of faces to recv */

  BFT_MALLOC(recv_count, n_ranks, int);
  BFT_MALLOC(recv_shift, n_ranks + 1, int);

  /* Exchange number of elements to send */

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);

  /* Build index arrays */

  recv_shift[0] = 0;
  for (rank = 0; rank < n_ranks; rank++)
    recv_shift[rank+1] = recv_shift[rank] + recv_count[rank];

  BFT_MALLOC(recv_gbuf, recv_shift[n_ranks]*2, cs_gnum_t);

  MPI_Alltoallv(send_gbuf, send_count, send_shift, CS_MPI_GNUM,
                recv_gbuf, recv_count, recv_shift, CS_MPI_GNUM, mpi_comm);

  /* Now switch from 1 to 2 entries in send and receive buffers */

  for (i = recv_shift[n_ranks] - 1; i > -1; i--)
    recv_gbuf[i*2] = recv_gbuf[i];

  /* Get the related global cell number for each received face */

  if (join_param.perio_type != FVM_PERIODICITY_NULL) {

    for (rank = 0; rank < n_ranks; rank++) {

      for (i = recv_shift[rank]; i < recv_shift[rank+1]; i++) {

        compact_fgnum = recv_gbuf[i*2];

        assert(loc_rank_s < compact_fgnum);
        assert(compact_fgnum <= loc_rank_e);

        /* Get the id in the original face selection */

        if (compact_fgnum % 2 == 0) { /* Periodic face */
          fid = join_select->faces[(compact_fgnum - loc_rank_s)/2 - 1] - 1;
          recv_gbuf[i*2] = mesh->global_cell_num[mesh->b_face_cells[fid]];
          recv_gbuf[i*2+1] = 0;
        }
        else { /* Original face */
          fid = join_select->faces[(compact_fgnum - loc_rank_s)/2] - 1;
          recv_gbuf[i*2] = mesh->global_cell_num[mesh->b_face_cells[fid]];
          recv_gbuf[i*2+1] = mesh->b_face_family[fid];
        }

      }

    } /* End of loop on ranks */

  }
  else { /* Not a periodic case */

    for (rank = 0; rank < n_ranks; rank++) {

      for (i = recv_shift[rank]; i < recv_shift[rank+1]; i++) {

        compact_fgnum = recv_gbuf[i*2];

        assert(loc_rank_s < compact_fgnum);
        assert(compact_fgnum <= loc_rank_e);

        fid = join_select->faces[compact_fgnum - 1 - loc_rank_s] - 1;
        if (mesh->b_face_cells[fid] > -1)
          recv_gbuf[i*2] = mesh->global_cell_num[mesh->b_face_cells[fid]];
        else
          recv_gbuf[i*2] = 0;
        recv_gbuf[i*2+1] = mesh->b_face_family[fid];

      }

    } /* End of loop on ranks */

  } /* End if not a periodic case */

  /* Update shifts */

  for (rank = 0; rank < n_ranks; rank++) {
    send_count[rank] *= 2;
    recv_count[rank] *= 2;
    send_shift[rank] *= 2;
    recv_shift[rank] *= 2;
  }
  send_shift[n_ranks] *= 2;
  recv_shift[n_ranks] *= 2;

  /* Return values to send ranks */

  MPI_Alltoallv(recv_gbuf, recv_count, recv_shift, CS_MPI_GNUM,
                send_gbuf, send_count, send_shift, CS_MPI_GNUM, mpi_comm);

  /* Define cell_gnum */

  for (i = 0; i < send_shift[n_ranks] / 2; i++) {
    cell_gnum[parent[i]] = send_gbuf[i*2];
    face_family[parent[i]] = send_gbuf[i*2+1];
  }

  /* Free memory */

  BFT_FREE(send_count);
  BFT_FREE(send_shift);
  BFT_FREE(send_gbuf);
  BFT_FREE(recv_count);
  BFT_FREE(recv_shift);
  BFT_FREE(recv_gbuf);
  BFT_FREE(parent);

}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Get the related global cell numbers and face families associated
 * to the old face numbers.
 *
 * parameters:
 *   join_select    <-- list of all implied entities in the joining op.
 *   join_param     <-- set of user-defined parameter
 *   n2o_face_hist  <-- face history structure (new -> old)
 *   mesh           <--  pointer to the associated cs_mesh_t structure
 *   cell_gnum      --> pointer to the created global cell number array
 *   face_family    --> pointer to the created face family array
 *---------------------------------------------------------------------------*/

static void
_get_linked_cell_gnum_and_family(const cs_join_select_t  *join_select,
                                 cs_join_param_t          join_param,
                                 const cs_join_gset_t    *n2o_face_hist,
                                 const cs_mesh_t         *mesh,
                                 cs_gnum_t               *p_cell_gnum[],
                                 cs_lnum_t               *p_face_family[])
{
  cs_lnum_t  i, j;
  cs_gnum_t  compact_fgnum;

  cs_gnum_t  *cell_gnum = NULL;
  cs_lnum_t   *face_family = NULL;

  const int  n_ranks = cs_glob_n_ranks;

  BFT_MALLOC(cell_gnum,
             n2o_face_hist->index[n2o_face_hist->n_elts],
             cs_gnum_t);

  BFT_MALLOC(face_family,
             n2o_face_hist->index[n2o_face_hist->n_elts],
             cs_lnum_t);

  if (n_ranks == 1) {

    if (join_param.perio_type != FVM_PERIODICITY_NULL) { /* Periodic case */

      for (i = 0; i < n2o_face_hist->n_elts; i++) {

        for (j = n2o_face_hist->index[i]; j < n2o_face_hist->index[i+1]; j++) {

          compact_fgnum = n2o_face_hist->g_list[j];

          if (compact_fgnum % 2 == 0) { /* Periodic face */
            cs_lnum_t fid = join_select->faces[compact_fgnum/2 - 1] - 1;
            if (mesh->global_cell_num != NULL)
              cell_gnum[j]
                = mesh->global_cell_num[mesh->b_face_cells[fid]];
            else
              cell_gnum[j] = mesh->b_face_cells[fid] + 1;
            face_family[j] = 0;
          }
          else { /* Original face */
            cs_lnum_t fid = join_select->faces[compact_fgnum/2] - 1;
            if (mesh->global_cell_num != NULL)
              cell_gnum[j]
                = mesh->global_cell_num[mesh->b_face_cells[fid]];
            else
              cell_gnum[j] = mesh->b_face_cells[fid] + 1;
            face_family[j] = mesh->b_face_family[fid];
          }

        }

      } /* End of loop on n2o_face_hist elements */

    }
    else { /* Non-periodic case */

      if (mesh->global_cell_num != NULL) {

        for (i = 0; i < n2o_face_hist->n_elts; i++) {
          for (j = n2o_face_hist->index[i]; j < n2o_face_hist->index[i+1]; j++) {
            compact_fgnum = n2o_face_hist->g_list[j];
            cs_lnum_t fid = join_select->faces[compact_fgnum - 1] - 1;
            cs_lnum_t cid = mesh->b_face_cells[fid];
            if (cid > -1)
              cell_gnum[j] = mesh->global_cell_num[cid];
            else
              cell_gnum[j] = 0;
            face_family[j] = mesh->b_face_family[fid];
          }
        } /* End of loop on n2o_face_hist elements */

      }
      else {

        for (i = 0; i < n2o_face_hist->n_elts; i++) {
          for (j = n2o_face_hist->index[i]; j < n2o_face_hist->index[i+1]; j++) {
            compact_fgnum = n2o_face_hist->g_list[j];
            cs_lnum_t fid = join_select->faces[compact_fgnum - 1] - 1;
            cell_gnum[j] = mesh->b_face_cells[fid] + 1;
            face_family[j] = mesh->b_face_family[fid];
          }
        } /* End of loop on n2o_face_hist elements */
      }

    } /* Not a periodic case */

  }

#if defined(HAVE_MPI)
  if (n_ranks > 1)
    _exchange_cell_gnum_and_family(n2o_face_hist,
                                   join_select,
                                   mesh,
                                   join_param,
                                   cell_gnum,
                                   face_family);
#endif

  /* Return pointer */

  *p_cell_gnum = cell_gnum;
  *p_face_family = face_family;
}

/*----------------------------------------------------------------------------
 * Print information before calling bft_error()
 *
 * parameters:
 *   jfnum    <--  new sub join face to reorient
 *   cgnum    <--  cgnum1 and cgnum2 related to this interior face
 *   fnum     <--  related selected face num in original face
 *   jmesh    <--  pointer to a cs_join_mesh_t structure
 *---------------------------------------------------------------------------*/

static void
_print_error_info(cs_lnum_t               jfnum,
                  const cs_gnum_t         cgnum[],
                  const cs_lnum_t         fnum[],
                  const cs_join_mesh_t   *jmesh)
{
  int  i, vid;
  int  jms = jmesh->face_vtx_idx[jfnum-1];
  int  jme = jmesh->face_vtx_idx[jfnum];

  if (cs_glob_join_log != NULL) {
    fprintf(cs_glob_join_log,
            "\n   cgnum (%llu, %llu) - fnum: (%d, %d)\n",
            (unsigned long long)cgnum[0], (unsigned long long)cgnum[1],
            fnum[0], fnum[1]);

    fprintf(cs_glob_join_log,
            "  Join Face connectivity %d (%llu): ",
            jfnum, (unsigned long long)jmesh->face_gnum[jfnum-1]);
    for (i = jms; i < jme; i++) {
      vid = jmesh->face_vtx_lst[i];
      fprintf(cs_glob_join_log, "%llu ",
              (unsigned long long)jmesh->vertices[vid].gnum);
    }
    fprintf(cs_glob_join_log, "\n");
    fflush(cs_glob_join_log);
  }

  /* Define a specific name for the output */

#if 0 && defined(DEBUG) && !defined(NDEBUG) /* Dump mesh structure */
  {
    int  len;
    FILE  *dbg_file = NULL;
    char  *fullname = NULL;

    const  int  rank_id = CS_MAX(cs_glob_rank_id, 0);

    len = strlen("JoinDBG_ErrorOrient.dat") + 4 + 1;
    BFT_MALLOC(fullname, len, char);
    sprintf(fullname, "JoinDBG_ErrorOrient%04d.dat", rank_id);

    dbg_file = fopen(fullname, "w");
    cs_join_mesh_dump(dbg_file, jmesh);
    fflush(dbg_file);
    fclose(dbg_file);
  }
#endif

  bft_error(__FILE__, __LINE__, 0,
            _("  Cannot achieve to reorient the current joined face.\n"));

}

/*----------------------------------------------------------------------------
 * Check if the orientation is the same between two faces thanks to its
 * normal
 *
 * parameters:
 *   omfnum  <--  old face number in mesh structure
 *   jmfnum  <--  new face number in join mesh structure
 *   mesh    <->  pointer to the original cs_mesh_t structure after the
 *                merge step
 *   jmesh   <--  pointer to a cs_join_mesh_t structure
 *   dtmp    <->  work buffer of double
 *   gtmp    <->  work buffer of gnum
 *
 * returns:
 *  1 if in same orientation, else -1. Return 0 if there is a problem.
 *---------------------------------------------------------------------------*/

static int
_get_geom_orient(cs_lnum_t               omfnum,
                 cs_lnum_t               jmfnum,
                 const cs_mesh_t        *mesh,
                 const cs_join_mesh_t   *jmesh,
                 double                  dtmp[])
{
  int  i, j, k, jvid, mvid;

  int  ret = 0;
  int  jmfid = jmfnum - 1, omfid = omfnum - 1;
  int  jms = jmesh->face_vtx_idx[jmfid];
  int  jme = jmesh->face_vtx_idx[jmfid+1];
  int  ms = mesh->b_face_vtx_idx[omfid];
  int  me = mesh->b_face_vtx_idx[omfid+1];
  int  jsize = jme - jms, size = me - ms;
  double  *jcoord = &(dtmp[0]);
  double  *coord = &(dtmp[3*(jsize+1)]);
  double  jnormal[3], normal[3];

  /* Fill work buffers */

  for (i = jms, j = 0; i < jme; i++, j++) {
    jvid = jmesh->face_vtx_lst[i];
    for (k = 0; k < 3; k++)
      jcoord[3*j+k] = jmesh->vertices[jvid].coord[k];
  }
  jvid = jmesh->face_vtx_lst[jms];
  for (k = 0; k < 3; k++)
    jcoord[3*j+k] = jmesh->vertices[jvid].coord[k];

  for (i = ms, j = 0; i < me; i++, j++) {
    mvid = mesh->b_face_vtx_lst[i];
    for (k = 0; k < 3; k++)
      coord[3*j+k] = mesh->vtx_coord[3*mvid+k];
  }
  mvid = mesh->b_face_vtx_lst[ms];
  for (k = 0; k < 3; k++)
    coord[3*j+k] = mesh->vtx_coord[3*mvid+k];

  _get_face_normal(jsize, jcoord, jnormal);
  _get_face_normal(size, coord, normal);

  if (_dot_product(jnormal, normal) < 0)
    ret = -1;
  else if (_dot_product(jnormal, normal) > 0)
    ret = 1;

  return ret;
}

/*----------------------------------------------------------------------------
 * Check if the orientation is the same between two faces thanks to its
 * topology
 *
 * parameters:
 *   omfnum  <--  old face number in mesh structure
 *   jmfnum  <--  new face number in join mesh structure
 *   mesh    <->  pointer to the original cs_mesh_t structure after the
 *                merge step
 *   jmesh   <--  pointer to a cs_join_mesh_t structure
 *   tmp     <->  work buffer
 *
 * returns:
 *  1 if in same orientation, else -1. If no common edge was found, return 0
 *---------------------------------------------------------------------------*/

static int
_get_topo_orient(cs_lnum_t               omfnum,
                 cs_lnum_t               jmfnum,
                 const cs_mesh_t        *mesh,
                 const cs_join_mesh_t   *jmesh,
                 cs_gnum_t               gtmp[])
{
  int  i, j, k, jvid, mvid;
  cs_gnum_t  ref1, ref2;

  int  ret = 0;
  int  jmfid = jmfnum - 1, omfid = omfnum - 1;
  int  jms = jmesh->face_vtx_idx[jmfid];
  int  jme = jmesh->face_vtx_idx[jmfid+1];
  int  ms = mesh->b_face_vtx_idx[omfid];
  int  me = mesh->b_face_vtx_idx[omfid+1];
  int  jsize = jme - jms, size = me - ms;
  cs_gnum_t  *jconnect = &(gtmp[0]);
  cs_gnum_t  *connect = &(gtmp[jsize+1]);

  /* Fill work buffers
     mesh->global_vtx_num is always allocated even if in serial run */

  for (i = jms, k = 0; i < jme; i++, k++) {
    jvid = jmesh->face_vtx_lst[i];
    jconnect[k] = jmesh->vertices[jvid].gnum;
  }
  jvid = jmesh->face_vtx_lst[jms];
  jconnect[k] = jmesh->vertices[jvid].gnum;

  for (i = ms, k = 0; i < me; i++, k++) {
    mvid = mesh->b_face_vtx_lst[i];
    connect[k] = mesh->global_vtx_num[mvid];
  }
  mvid = mesh->b_face_vtx_lst[ms];
  connect[k] = mesh->global_vtx_num[mvid];

  /* Find a common edge between the two face connectivities */

  for (i = 0; i < jsize && ret == 0; i++) {

    ref1 = jconnect[i];
    ref2 = jconnect[i+1];

    for (j = 0; j < size; j++) {
      if (ref2 == connect[j])
        if (ref1 == connect[j+1])
          ret = -1; /* Opposite scan order */
      if (ref1 == connect[j])
        if (ref2 == connect[j+1])
          ret = 1; /* Same scan order */
    }

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n JoinFace num: %d - connect: [", jmfnum);
  for (k = 0; k < jsize; k++)
    bft_printf(" %llu", (unsigned long long)(jconnect[k]));
  bft_printf(" %u]\n", jconnect[jsize]);
  bft_printf(" OldMeshFace num: %d - connect: [", omfnum);
  for (k = 0; k < size; k++)
    bft_printf(" %llu", (unsigned long long)(connect[k]));
  bft_printf(" %llu]\n", (unsigned long long)(connect[size]));
  bft_printf("ret: %d\n", ret);
  bft_printf_flush();
#endif

  return ret;
}

/*----------------------------------------------------------------------------
 * Check if the orientation is the same between two faces. Topological algo.
 *
 * parameters:
 *   jfnum    <--  new sub join face to reorient
 *   cgnum    <--  cgnum1 and cgnum2 related to this interior face
 *   fnum     <--  related selected face num in original face
 *   mesh     <->  pointer to the original cs_mesh_t structure after the
 *                 merge step
 *   jmesh    <--  pointer to a cs_join_mesh_t structure
 *   ltmp     <->  work buffer to store local connectivity
 *   dtmp     <->  work buffer to compute face normal (only in geometric algo)
 *   gtmp     <->  work buffer to store global numbering
 *---------------------------------------------------------------------------*/

static void
_reorient(cs_lnum_t               jfnum,
          cs_gnum_t               cgnum[],
          cs_lnum_t               fnum[],
          const cs_mesh_t        *mesh,
          const cs_join_mesh_t   *jmesh,
          cs_lnum_t               ltmp[],
          cs_gnum_t               gtmp[],
          double                  dtmp[])
{
  int  i, k, orient_tag;

  int  jms = jmesh->face_vtx_idx[jfnum-1];
  int  jme = jmesh->face_vtx_idx[jfnum];
  int  n_face_vertices = jme - jms;

  if (cs_glob_n_ranks == 1)
    if (fnum[0] == 0 || fnum[1] == 0)
      _print_error_info(jfnum, cgnum, fnum, jmesh);

  /* Get the right orientation. We assume that all original border faces
     have a right orientation i.e. outward normal orientation */

  if (cgnum[0] < cgnum[1]) { /* fnum[0] forces the orientation */

    if (fnum[0] > 0) { /* fnum[0] belongs to local_rank */

      orient_tag = _get_topo_orient(fnum[0], jfnum, mesh, jmesh, gtmp);

      if (orient_tag == -1) { /* fnum[0] and jfnum don't share the same
                                 orientation => Re-orient */

        for (i = jms, k = 0; i < jme; i++, k++)
          ltmp[k] = jmesh->face_vtx_lst[i];
        for (i = jms + 1, k = 1; i < jme; i++, k++)
          jmesh->face_vtx_lst[i] = ltmp[n_face_vertices - k];

      }
      else if (orient_tag == 0) { /* Edge not found in fnum[0] */

        if (fnum[1] < 1) { /* Switch to a geometrical test
                              fnum[1] should own a common edge and orientation
                              should be opposite to fnum[0].
                              So we test face orientation and check this
                              assumption. */

          orient_tag = _get_geom_orient(fnum[0], jfnum, mesh, jmesh, dtmp);

          if (orient_tag < 0) { /* Need to reorient */

            for (i = jms, k = 0; i < jme; i++, k++)
              ltmp[k] = jmesh->face_vtx_lst[i];
            for (i = jms + 1, k = 1; i < jme; i++, k++)
              jmesh->face_vtx_lst[i] = ltmp[n_face_vertices - k];

          }
          else if (orient_tag == 0)
            _print_error_info(jfnum, cgnum, fnum, jmesh);

        }
        else { /* fnum[1] > 0 */

          orient_tag = _get_topo_orient(fnum[1], jfnum, mesh, jmesh, gtmp);

          if (orient_tag == 1) { /* fnum[1] and jfnum share the same
                                    orientation => Re-orient */

            for (i = jms, k = 0; i < jme; i++, k++)
              ltmp[k] = jmesh->face_vtx_lst[i];
            for (i = jms + 1, k = 1; i < jme; i++, k++)
              jmesh->face_vtx_lst[i] = ltmp[n_face_vertices - k];

          }
          else if (orient_tag == 0) { /* Edge not found in fnum[1] */

            /* Switch to a geometrical test */

            orient_tag = _get_geom_orient(fnum[1], jfnum, mesh, jmesh, dtmp);

            if (orient_tag == 1) { /* fnum[1] and jfnum share the same
                                      orientation => Re-orient */

              for (i = jms, k = 0; i < jme; i++, k++)
                ltmp[k] = jmesh->face_vtx_lst[i];
              for (i = jms + 1, k = 1; i < jme; i++, k++)
                jmesh->face_vtx_lst[i] = ltmp[n_face_vertices - k];

            }
            else if (orient_tag == 0)
              _print_error_info(jfnum, cgnum, fnum, jmesh);

          }

        } /* fnum[1] */

      } /* Edge not found in fnum[0] */

    } /* fnum[0] > 0 */

    else { /* fnum[0] < 1 */

      if (fnum[1] > 0) {

        orient_tag = _get_topo_orient(fnum[1], jfnum, mesh, jmesh, gtmp);

        if (orient_tag == 1) { /* fnum[1] and jfnum share the same orientation
                                  => Re-orient */

          for (i = jms, k = 0; i < jme; i++, k++)
            ltmp[k] = jmesh->face_vtx_lst[i];
          for (i = jms + 1, k = 1; i < jme; i++, k++)
            jmesh->face_vtx_lst[i] = ltmp[n_face_vertices - k];

        }
        else if (orient_tag == 0) {

          /* Switch to a geometrical test */

          orient_tag = _get_geom_orient(fnum[1], jfnum, mesh, jmesh, dtmp);

          if (orient_tag == 1) { /* fnum[1] and jfnum share the same
                                    orientation => Re-orient */

            for (i = jms, k = 0; i < jme; i++, k++)
              ltmp[k] = jmesh->face_vtx_lst[i];
            for (i = jms + 1, k = 1; i < jme; i++, k++)
              jmesh->face_vtx_lst[i] = ltmp[n_face_vertices - k];

          }
          else if (orient_tag == 0) /* Geom. and topo. tests fail */
            _print_error_info(jfnum, cgnum, fnum, jmesh);

        }

      } /* fnum[1] > 0 */

      else  /* fnum[1] < 1 && fnum[0] < 1 */
        _print_error_info(jfnum, cgnum, fnum, jmesh);

    } /* fnum[0] < 1 */

  }
  else { /* cgnum[0] > cgnum[1] => fnum[1] forces the orientation  */

    if (fnum[1] > 0) { /* fnum[1] belongs to local_rank */

      orient_tag = _get_topo_orient(fnum[1], jfnum, mesh, jmesh, gtmp);

      if (orient_tag == -1) { /* fnum[1] and jfnum don't share the same
                                 orientation => re-orient */

        for (i = jms, k = 0; i < jme; i++, k++)
          ltmp[k] = jmesh->face_vtx_lst[i];
        for (i = jms + 1, k = 1; i < jme; i++, k++)
          jmesh->face_vtx_lst[i] = ltmp[n_face_vertices - k];

      }
      else if (orient_tag == 0) { /* Edge not found in fnum[1] */

        if (fnum[0] < 1) { /* Switch to a geomtrical test
                              fnum[1] should own a common edge and orientation
                              should be opposite to fnum[0].
                              So we test face orientation and check this
                              assumption. */

          orient_tag = _get_geom_orient(fnum[1], jfnum, mesh, jmesh, dtmp);

          if (orient_tag < 0) { /* Need to reorient */

            for (i = jms, k = 0; i < jme; i++, k++)
              ltmp[k] = jmesh->face_vtx_lst[i];
            for (i = jms + 1, k = 1; i < jme; i++, k++)
              jmesh->face_vtx_lst[i] = ltmp[n_face_vertices - k];

          }
          else if (orient_tag == 0)  /* Geom. and topo. tests fail */
            _print_error_info(jfnum, cgnum, fnum, jmesh);

        }
        else { /* fnum[0] > 0 */

          orient_tag = _get_topo_orient(fnum[0], jfnum, mesh, jmesh, gtmp);

          if (orient_tag == 1) { /* fnum[0] and jfnum share the same
                                    orientation => re-orient */

            for (i = jms, k = 0; i < jme; i++, k++)
              ltmp[k] = jmesh->face_vtx_lst[i];
            for (i = jms + 1, k = 1; i < jme; i++, k++)
              jmesh->face_vtx_lst[i] = ltmp[n_face_vertices - k];

          }
          else if (orient_tag == 0) { /* Edge not found in fnum[0] */

            /* Switch to a geometrical test */

            orient_tag = _get_geom_orient(fnum[0], jfnum, mesh, jmesh, dtmp);

            if (orient_tag == 1) { /* Need to reorient */

              for (i = jms, k = 0; i < jme; i++, k++)
                ltmp[k] = jmesh->face_vtx_lst[i];
              for (i = jms + 1, k = 1; i < jme; i++, k++)
                jmesh->face_vtx_lst[i] = ltmp[n_face_vertices - k];

            }
            else if (orient_tag == 0)  /* Geom. and topo. tests fail */
              _print_error_info(jfnum, cgnum, fnum, jmesh);

          }

        } /* fnum[0] > 0 */

      } /* Edge not found in fnum[1] */

    }
    else { /* fnum[1] < 1 */

      if (fnum[0] > 0) { /* fnum[0] belongs to local_rank */

        orient_tag = _get_topo_orient(fnum[0], jfnum, mesh, jmesh, gtmp);

        if (orient_tag == 1) { /* fnum[0] and jfnum share the same
                                  orientation => re-orient */

          for (i = jms, k = 0; i < jme; i++, k++)
            ltmp[k] = jmesh->face_vtx_lst[i];
          for (i = jms + 1, k = 1; i < jme; i++, k++)
            jmesh->face_vtx_lst[i] = ltmp[n_face_vertices - k];

        }
        else if (orient_tag == 0) { /* Edge not found in fnum[0] */

          /* Switch to a geometrical test */

          orient_tag = _get_geom_orient(fnum[0], jfnum, mesh, jmesh, dtmp);

          if (orient_tag == 1) { /* Need to reorient */

            for (i = jms, k = 0; i < jme; i++, k++)
              ltmp[k] = jmesh->face_vtx_lst[i];
            for (i = jms + 1, k = 1; i < jme; i++, k++)
              jmesh->face_vtx_lst[i] = ltmp[n_face_vertices - k];

          }
          else if (orient_tag == 0)  /* Geom. and topo. tests fail */
            _print_error_info(jfnum, cgnum, fnum, jmesh);

        }

      } /* fnum[0] > 0 */

    } /* fnum[1] < 1 */

  } /* cgnum[1] < cgnum[0] */

}

/*----------------------------------------------------------------------------
 * Update mesh structure by adding new border faces after the face cutting
 * and deleting old one.
 *
 * parameters:
 *   join_select      <-- list of all implied entities in the joining op.
 *   join_param       <-- set of user-defined parameter
 *   jmesh            <-- pointer to the local cs_join_mesh_t structure
 *   join2mesh_vtx_id <-- relation between vertices in join_mesh/mesh
 *   n_new_b_faces    <-- local number of border faces after the joining
 *   new_face_type    <-- type (border/interior) of new faces
 *   new_face_family  <-- family of new faces
 *   n2o_face_hist    <-- face history structure (new -> old)
 *   mesh             <-> pointer of pointer to cs_mesh_t structure
 *---------------------------------------------------------------------------*/

static void
_add_new_border_faces(const cs_join_select_t     *join_select,
                      cs_join_param_t             join_param,
                      const cs_join_mesh_t       *jmesh,
                      const cs_lnum_t             join2mesh_vtx_id[],
                      cs_lnum_t                   n_new_b_faces,
                      const cs_join_face_type_t   new_face_type[],
                      const int                   new_face_family[],
                      const cs_join_gset_t       *n2o_face_hist,
                      cs_mesh_t                  *mesh)
{
  cs_lnum_t  i, j, k, select_id, vid, fid, shift;
  cs_lnum_t  n_face_vertices, max_size, orient_tag;
  cs_gnum_t  compact_old_fgnum;

  cs_lnum_t   n_ib_faces = mesh->n_b_faces, n_fb_faces = 0;
  cs_gnum_t  n_g_ib_faces = mesh->n_g_b_faces;
  cs_lnum_t  *new_f2v_idx = NULL, *new_f2v_lst = NULL, *ltmp = NULL;
  cs_lnum_t  *_new_face_family = NULL, *new_face_cells = NULL;
  cs_gnum_t  *new_fgnum = NULL, *gtmp = NULL;

  const int  n_ranks = cs_glob_n_ranks;
  const int  rank = CS_MAX(cs_glob_rank_id, 0);
  const cs_gnum_t  rank_start = join_select->compact_rank_index[rank] + 1;
  const cs_gnum_t  rank_end = join_select->compact_rank_index[rank+1] + 1;

  n_fb_faces = n_ib_faces + n_new_b_faces - join_select->n_faces;
  mesh->n_b_faces = n_fb_faces;
  mesh->n_g_b_faces = n_fb_faces;

  BFT_MALLOC(new_f2v_idx, n_fb_faces + 1, cs_lnum_t);
  BFT_MALLOC(new_face_cells, n_fb_faces, cs_lnum_t);
  BFT_MALLOC(_new_face_family, n_fb_faces, cs_lnum_t);

  if (n_ranks > 1)
    BFT_MALLOC(new_fgnum, n_fb_faces, cs_gnum_t);

  max_size = 0;
  for (i = 0; i < n_ib_faces; i++)
    max_size = CS_MAX(max_size,
                      mesh->b_face_vtx_idx[i+1]-mesh->b_face_vtx_idx[i]);
  for (i = 0; i < jmesh->n_faces; i++)
    max_size = CS_MAX(max_size,
                      jmesh->face_vtx_idx[i+1]-jmesh->face_vtx_idx[i]);

  BFT_MALLOC(gtmp, 2*(max_size+1), cs_gnum_t);
  BFT_MALLOC(ltmp, max_size, cs_lnum_t);

  /* Delete faces included in join_selection. Add other initial faces.
      - face -> vertex index
      - face -> cells connectivity
      - face family
      - face global num. (first pass)
  */

  n_fb_faces = 0;
  for (i = 0, select_id = 0; i < n_ib_faces; i++) {

    bool  in_selection = false;

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
      _new_face_family[n_fb_faces] = mesh->b_face_family[i];

      n_fb_faces++;
      n_face_vertices = mesh->b_face_vtx_idx[i+1] - mesh->b_face_vtx_idx[i];
      new_f2v_idx[n_fb_faces] = n_face_vertices;

    }

  }

  assert(n_fb_faces == n_ib_faces - join_select->n_faces);

  /* Add faces resulting from the joining operation */

  if (n_new_b_faces > 0) {
    for (i = 0; i < jmesh->n_faces; i++) {

      if (new_face_type[i] == CS_JOIN_FACE_BORDER) {

        int  jms = jmesh->face_vtx_idx[i];
        int  jme = jmesh->face_vtx_idx[i+1];

        n_face_vertices = jme - jms;
        j = n2o_face_hist->index[i];
        compact_old_fgnum = n2o_face_hist->g_list[j];

        /* Initial selected border face must be in the selection */
        assert(rank_start <= compact_old_fgnum);
        assert(compact_old_fgnum < rank_end);

        if (join_param.perio_type != FVM_PERIODICITY_NULL) {
          fid = (compact_old_fgnum - rank_start)/2;
          fid = join_select->faces[fid] - 1;
        }
        else
          fid = join_select->faces[compact_old_fgnum - rank_start] - 1;

        new_face_cells[n_fb_faces] = mesh->b_face_cells[fid];
        _new_face_family[n_fb_faces] = new_face_family[i];

        if (n_ranks > 1)
          new_fgnum[n_fb_faces] = jmesh->face_gnum[i] + n_g_ib_faces;

        /* Check orientation: fid forces the orientation */

        orient_tag = _get_topo_orient(fid+1, i+1, mesh, jmesh, gtmp);

        if (orient_tag == -1) { /* Different orientation => re-orient*/
          for (j = jms, k = 0; j < jme; j++, k++)
            ltmp[k] = jmesh->face_vtx_lst[j];
          for (j = jms + 1, k = 1; j < jme; j++, k++)
          jmesh->face_vtx_lst[j] = ltmp[n_face_vertices - k];
        }
        else if (orient_tag == 0) { /* No common edge found */

          for (j = jms; j < jme; j++) {
            vid = jmesh->face_vtx_lst[j];
            bft_printf(" %d (%llu)", vid+1,
                       (unsigned long long)(jmesh->vertices[vid].gnum));
          }
          bft_printf("\n");
          bft_printf_flush();

          bft_error(__FILE__, __LINE__, 0,
                    _("  Cannot achieve to reorient the current joined"
                      " face with face %d (selected face).\n"), fid+1);

        }

        n_fb_faces++;
        new_f2v_idx[n_fb_faces] = n_face_vertices;

      } /* If new border face */

      else if (new_face_type[i] == CS_JOIN_FACE_MULTIPLE_BORDER) {

        int  jms = jmesh->face_vtx_idx[i];
        int  jme = jmesh->face_vtx_idx[i+1];

        n_face_vertices = jme - jms;

        new_face_cells[n_fb_faces] = -1;

        /* When a face with multiple parents is present on */
        for (j = n2o_face_hist->index[i]; j < n2o_face_hist->index[i+1]; j++) {

          compact_old_fgnum = n2o_face_hist->g_list[j];

          if (    rank_start <= compact_old_fgnum
              &&  compact_old_fgnum < rank_end) {

            /* Initial selected border face must be in the selection */

            if (join_param.perio_type != FVM_PERIODICITY_NULL) {
              fid = (compact_old_fgnum - rank_start)/2;
              fid = join_select->faces[fid] - 1;
            }
            else
              fid = join_select->faces[compact_old_fgnum - rank_start] - 1;

            if (mesh->b_face_cells[fid] > -1) {

              new_face_cells[n_fb_faces] = mesh->b_face_cells[fid];

              /* Check orientation: fid forces the orientation */
              orient_tag = _get_topo_orient(fid+1, i+1, mesh, jmesh, gtmp);

              if (orient_tag < 0) { /* Different orientation => re-orient*/
                for (j = jms, k = 0; j < jme; j++, k++)
                  ltmp[k] = jmesh->face_vtx_lst[j];
                for (j = jms + 1, k = 1; j < jme; j++, k++)
                  jmesh->face_vtx_lst[j] = ltmp[n_face_vertices - k];
              }
              else if (orient_tag == 0) { /* No common edge found */

                for (j = jms; j < jme; j++)
                  vid = jmesh->face_vtx_lst[j];

                bft_error(__FILE__, __LINE__, 0,
                          _("  Cannot achieve to reorient the current joined"
                            " face with face %d (selected face).\n"), fid+1);

              }

            }

            _new_face_family[n_fb_faces] = new_face_family[i];
          }

        }

        if (n_ranks > 1)
          new_fgnum[n_fb_faces] = jmesh->face_gnum[i] + n_g_ib_faces;

        n_fb_faces++;
        new_f2v_idx[n_fb_faces] = n_face_vertices;

      }

    } /* Loop on faces */

  } /* If n_new_b_faces > 0 */

  BFT_FREE(gtmp);
  BFT_FREE(ltmp);

  assert(mesh->n_b_faces == n_fb_faces);

  /* Build index */

  new_f2v_idx[0] = 0;
  for (i = 0; i < n_fb_faces; i++)
    new_f2v_idx[i+1] += new_f2v_idx[i];

  BFT_MALLOC(new_f2v_lst, new_f2v_idx[n_fb_faces], cs_lnum_t);

  /* Define the face -> vertex connectivity */

  n_fb_faces = 0;
  for (i = 0, select_id = 0; i < n_ib_faces; i++) {

    bool  in_selection = false;

    if (select_id < join_select->n_faces) {
      if (i+1 == join_select->faces[select_id]) {
        in_selection = true;
        select_id++;
      }
    }

    if (in_selection == false) {

      shift = new_f2v_idx[n_fb_faces];

      for (j = mesh->b_face_vtx_idx[i]; j < mesh->b_face_vtx_idx[i+1]; j++)
        new_f2v_lst[shift++] = mesh->b_face_vtx_lst[j];

      n_fb_faces++;

    }

  }

  if (n_new_b_faces > 0) {
    for (i = 0; i < jmesh->n_faces; i++) {
      if (   new_face_type[i] == CS_JOIN_FACE_BORDER
          || new_face_type[i] == CS_JOIN_FACE_MULTIPLE_BORDER) {

        shift = new_f2v_idx[n_fb_faces];

        for (j = jmesh->face_vtx_idx[i];
             j < jmesh->face_vtx_idx[i+1]; j++) {
          vid = jmesh->face_vtx_lst[j];
          new_f2v_lst[shift++] = join2mesh_vtx_id[vid];
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

    const cs_gnum_t  *new_io_gnum = fvm_io_num_get_global_num(new_io_num);

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
  mesh->b_face_family = _new_face_family;
  mesh->b_face_vtx_connect_size = new_f2v_idx[n_fb_faces];
}

/*----------------------------------------------------------------------------
 * Update mesh families based on combined faces.
 *
 * parameters:
 *   n2o_face_hist   <-- relation between faces before/after the joining
 *   join_mesh       <-> pointer to the local cs_join_mesh_t structure
 *   mesh            <-> pointer of pointer to cs_mesh_t structure
 *
 * returns:
 *   array of family numbers fo new faces
 *---------------------------------------------------------------------------*/

static int *
_update_families(const cs_join_gset_t    *n2o_face_hist,
                 cs_lnum_t                old_face_family[],
                 cs_join_mesh_t          *join_mesh,
                 cs_mesh_t               *mesh)
{
  cs_lnum_t  i, j, k;

  int  null_family = 0;
  cs_lnum_t  *face_family_idx = NULL;
  int  *face_family = NULL;
  int  *new_face_family = NULL;

  assert(mesh != NULL);
  assert(join_mesh != NULL);

  /* Get new subfaces evolution */

  assert(n2o_face_hist->n_elts == join_mesh->n_faces);

  BFT_MALLOC(face_family_idx, join_mesh->n_faces + 1, cs_lnum_t);
  BFT_MALLOC(face_family, n2o_face_hist->index[join_mesh->n_faces], int);

  /* Compact numbering (remove zeroes) */

  if (mesh->n_families > 0) {
    if (mesh->family_item[0] == 0)
      null_family = 1;
  }

  face_family_idx[0] = 0;
  k = 0;

  for (i = 0; i < join_mesh->n_faces; i++) {

    int prev_fam = 0;
    const cs_lnum_t start_id = n2o_face_hist->index[i];
    const cs_lnum_t end_id = n2o_face_hist->index[i+1];

    cs_sort_shell(start_id, end_id, old_face_family);

    for (j = start_id; j < end_id; j++) {
      if (old_face_family[j] != prev_fam) {
        face_family[k] = old_face_family[j];
        prev_fam = old_face_family[j];
        k += 1;
      }
    }

    /* Add null family if none other was added */

    if (k == face_family_idx[i])
      face_family[k++] = null_family;

    face_family_idx[i+1] = k;
  }

  BFT_REALLOC(face_family, face_family_idx[join_mesh->n_faces], int);

  /* Build new combined families if necessary and flatten element families */

  if (mesh->n_families > 0) {
    BFT_MALLOC(new_face_family, join_mesh->n_faces, int);
    cs_mesh_group_combine_classes(mesh,
                                  join_mesh->n_faces,
                                  face_family_idx,
                                  face_family,
                                  new_face_family);
  }

  BFT_FREE(face_family);
  BFT_FREE(face_family_idx);

  return new_face_family;
}

/*----------------------------------------------------------------------------
 * Update mesh structure by adding new interior faces after the face split.
 * Orient new interior faces.
 *
 * parameters:
 *   join_select      <-- list of all implied entities in the joining op.
 *   join_param       <-- set of user-defined parameter
 *   jmesh            <-- pointer to the local cs_join_mesh_t structure
 *   join2mesh_vtx_id <-- relation between vertices in join_mesh/mesh
 *   cell_gnum        <-- global cell num. related to each initial face
 *   n_new_i_faces    <-- local number of interior faces after the joining
 *   new_face_type    <-- type (border/interior) of new faces
 *   new_face_family  <-- family of new faces
 *   n2o_face_hist    <-- face history structure (new -> old)
 *   mesh             <-> pointer of pointer to cs_mesh_t structure
 *---------------------------------------------------------------------------*/

static void
_add_new_interior_faces(const cs_join_select_t     *join_select,
                        cs_join_param_t             join_param,
                        cs_join_mesh_t             *jmesh,
                        const cs_lnum_t             join2mesh_vtx_id[],
                        const cs_gnum_t             cell_gnum[],
                        cs_lnum_t                   n_new_i_faces,
                        cs_join_face_type_t         new_face_type[],
                        const int                   new_face_family[],
                        const cs_join_gset_t       *n2o_face_hist,
                        cs_mesh_t                  *mesh)
{
  cs_lnum_t  i, j, k, vid, id, shift, fnum[2], max_size;
  cs_gnum_t  compact_fgnum, cgnum[2];

  cs_gnum_t  *gtmp = NULL;
  double  *dtmp = NULL;
  cs_lnum_t  *ltmp = NULL;
  cs_lnum_t  n_fi_faces = 0, n_ii_faces = mesh->n_i_faces;
  cs_lnum_t  *new_f2v_idx = mesh->i_face_vtx_idx;
  cs_lnum_t  *new_f2v_lst = mesh->i_face_vtx_lst;
  cs_lnum_t   *_new_face_family = mesh->i_face_family;
  cs_lnum_2_t  *new_face_cells = mesh->i_face_cells;
  cs_gnum_t  n_g_ii_faces = mesh->n_g_i_faces;
  cs_gnum_t  *new_fgnum = mesh->global_i_face_num;

  const int  n_ranks = cs_glob_n_ranks;
  const int  rank = CS_MAX(cs_glob_rank_id, 0);
  const cs_gnum_t  rank_start = join_select->compact_rank_index[rank] + 1;
  const cs_gnum_t  rank_end = join_select->compact_rank_index[rank+1] + 1;

  assert(mesh->global_vtx_num != NULL); /* Even if in serial run */

  n_fi_faces = n_ii_faces + n_new_i_faces;
  mesh->n_i_faces = n_fi_faces;
  mesh->n_g_i_faces = n_fi_faces;

  BFT_REALLOC(new_f2v_idx, n_fi_faces + 1, cs_lnum_t);
  BFT_REALLOC(new_face_cells, n_fi_faces, cs_lnum_2_t);
  BFT_REALLOC(_new_face_family, n_fi_faces, cs_lnum_t);

  max_size = 0;
  for (i = 0; i < jmesh->n_faces; i++)
    if (new_face_type[i] == CS_JOIN_FACE_INTERIOR)
      max_size = CS_MAX(max_size,
                        jmesh->face_vtx_idx[i+1]-jmesh->face_vtx_idx[i]);
  for (i = 0; i < join_select->n_faces; i++) {
    id = join_select->faces[i] - 1;
    max_size = CS_MAX(max_size,
                      mesh->b_face_vtx_idx[id+1]-mesh->b_face_vtx_idx[id]);
  }

  BFT_MALLOC(dtmp, 6*(max_size+1), double);
  BFT_MALLOC(gtmp, 2*(max_size+1), cs_gnum_t);
  BFT_MALLOC(ltmp, max_size, cs_lnum_t);

  /* Add faces resulting from the joining operation
     - face -> vertex index
     - face -> cells connectivity
     - face family
  */

  n_fi_faces = n_ii_faces;
  for (i = 0; i < jmesh->n_faces; i++) {

    if (new_face_type[i] == CS_JOIN_FACE_INTERIOR) {

      cs_lnum_t  jms = jmesh->face_vtx_idx[i];
      cs_lnum_t  jme = jmesh->face_vtx_idx[i+1];
      cs_lnum_t  n_face_vertices = jme - jms;

      for (j = n2o_face_hist->index[i], k = 0;
           j < n2o_face_hist->index[i+1];
           j++, k++) {

        cgnum[k] = cell_gnum[j];
        compact_fgnum = n2o_face_hist->g_list[j];

        if (rank_start <= compact_fgnum && compact_fgnum < rank_end) {

          if (join_param.perio_type != FVM_PERIODICITY_NULL) {

            id = (compact_fgnum - rank_start)/2;
            if (compact_fgnum % 2 == 0) /* Periodic face */
              fnum[k] = -join_select->faces[id];
            else /* Original face */
              fnum[k] = join_select->faces[id];

          }
          else /* Not a periodic operation */
            fnum[k] = join_select->faces[compact_fgnum - rank_start];

        }
        else
          fnum[k] = 0;

      }

      if (cgnum[0] < cgnum[1]) { /* Keep the same order and fnum[0]
                                    forces the orientation */

        if (fnum[0] > 0)
          new_face_cells[n_fi_faces][0] = mesh->b_face_cells[fnum[0]-1];
        else
          new_face_cells[n_fi_faces][0] = -1; /* Cell is on a distant rank
                                                 or periodic */

        if (fnum[1] > 0)
          new_face_cells[n_fi_faces][1] = mesh->b_face_cells[fnum[1]-1];
        else
          new_face_cells[n_fi_faces][1] = -1; /* Cell is on a distant rank
                                                 or periodic */

      }
      else { /* cgnum[1] < cgnum[0] and fnum[1] forces the orientation  */

        if (fnum[0] > 0)
          new_face_cells[n_fi_faces][1] = mesh->b_face_cells[fnum[0]-1];
        else
          new_face_cells[n_fi_faces][1] = -1; /* Cell is on a distant rank
                                                 or periodic */

        if (fnum[1] > 0)
          new_face_cells[n_fi_faces][0] = mesh->b_face_cells[fnum[1]-1];
        else
          new_face_cells[n_fi_faces][0] = -1; /* Cell is on a distant rank
                                                 or periodic */

      }

      /* Re-orient face in order to be consistent with the new i_face_cells */

      if (fnum[0] > 0 || fnum[1] > 0)
        _reorient(i+1, cgnum, fnum, mesh, jmesh, ltmp, gtmp, dtmp);
      else
        if (   join_param.perio_type == FVM_PERIODICITY_NULL
            || cs_glob_n_ranks == 1)
          bft_error(__FILE__, __LINE__, 0,
                    _(" Incoherency found before interior face orientation"
                      " checking.\n"
                      " Join face %d (%llu) and related faces [%d, %d]\n"),
                    i+1, (unsigned long long)(jmesh->face_gnum[i]),
                    fnum[0], fnum[1]);

      _new_face_family[n_fi_faces] = new_face_family[i];

      n_fi_faces++;
      new_f2v_idx[n_fi_faces] = n_face_vertices;

    } /* New interior face */

  } /* End of loop on jmesh faces */

  BFT_FREE(dtmp);
  BFT_FREE(gtmp);
  BFT_FREE(ltmp);

  assert(mesh->n_i_faces == n_fi_faces);

  /* Build index */

  for (i = n_ii_faces; i < n_fi_faces; i++)
    new_f2v_idx[i+1] += new_f2v_idx[i];

  BFT_REALLOC(new_f2v_lst, new_f2v_idx[n_fi_faces], cs_lnum_t);

  /* Define the face -> vertex connectivity list */

  n_fi_faces = n_ii_faces;
  for (i = 0; i < jmesh->n_faces; i++) {
    if (new_face_type[i] == CS_JOIN_FACE_INTERIOR) {

      shift = new_f2v_idx[n_fi_faces];

      for (j = jmesh->face_vtx_idx[i];
           j < jmesh->face_vtx_idx[i+1]; j++) {
        vid = jmesh->face_vtx_lst[j];
        new_f2v_lst[shift++] = join2mesh_vtx_id[vid];
      }
      n_fi_faces++;

    }
  }

  if (n_ranks > 1) { /* Get a compact global face numbering */

    fvm_io_num_t *new_io_num = NULL;
    const cs_gnum_t  *new_io_gnum = NULL;

    BFT_REALLOC(new_fgnum, mesh->n_i_faces, cs_gnum_t);

    n_fi_faces = n_ii_faces;
    for (i = 0; i < jmesh->n_faces; i++) {
      if (new_face_type[i] == CS_JOIN_FACE_INTERIOR)
        new_fgnum[n_fi_faces++] = jmesh->face_gnum[i] + n_g_ii_faces;
    }

    new_io_num = fvm_io_num_create(NULL, new_fgnum, n_fi_faces, 0);
    new_io_gnum = fvm_io_num_get_global_num(new_io_num);
    mesh->n_g_i_faces = fvm_io_num_get_global_count(new_io_num);

    for (i = 0; i < n_fi_faces; i++)
      new_fgnum[i] = new_io_gnum[i];

    fvm_io_num_destroy(new_io_num);

  }

  /* Update structure */

  mesh->n_i_faces = n_fi_faces;

  mesh->i_face_vtx_idx = new_f2v_idx;
  mesh->i_face_vtx_lst = new_f2v_lst;
  mesh->global_i_face_num = new_fgnum;
  mesh->i_face_cells = new_face_cells;
  mesh->i_face_family = _new_face_family;
  mesh->i_face_vtx_connect_size = new_f2v_idx[n_fi_faces];
}

/*----------------------------------------------------------------------------
 * Delete unused vertices
 *
 * parameters:
 *   param   <-- set of parameters for the joining operation
 *   mesh    <-> cs_mesh_t structure to clean
 *---------------------------------------------------------------------------*/

static void
_clean_vertices(cs_join_param_t   param,
                cs_mesh_t        *mesh)
{
  int  i, j, k, vid;

  cs_lnum_t  n_i_vertices = mesh->n_vertices, n_f_vertices = 0;
  cs_lnum_t  *tag = NULL;

  const int  n_ranks = cs_glob_n_ranks;

  assert(mesh != NULL);

  /* Tag vertices really used in the mesh definition */

  BFT_MALLOC(tag, n_i_vertices, cs_lnum_t);

  for (i = 0; i < n_i_vertices; i++)
    tag[i] = 0;

  for (i = 0; i < mesh->n_i_faces; i++)
    for (j = mesh->i_face_vtx_idx[i]; j < mesh->i_face_vtx_idx[i+1]; j++)
      tag[mesh->i_face_vtx_lst[j]] = 1;

  for (i = 0; i < mesh->n_b_faces; i++)
    for (j = mesh->b_face_vtx_idx[i]; j < mesh->b_face_vtx_idx[i+1]; j++)
      tag[mesh->b_face_vtx_lst[j]] = 1;

  for (i = 0; i < n_i_vertices; i++)
    if (tag[i] > 0)
      n_f_vertices++;

  if (param.verbosity > 3)
    fprintf(cs_glob_join_log,
            "  Delete %d local vertices not used in mesh definition.\n",
            n_i_vertices - n_f_vertices);

  mesh->n_vertices = n_f_vertices;
  mesh->n_g_vertices = n_f_vertices;

  /* Update global vertex information */

  if (n_ranks == 1) /* not more useful for sequential computations */
    BFT_FREE(mesh->global_vtx_num);

  if (n_ranks > 1) {

    fvm_io_num_t *vtx_io_num = NULL;
    const cs_gnum_t  *io_gnum = NULL;

    /* Define a new compact global vertex numbering */

    n_f_vertices = 0;
    for (i = 0; i < n_i_vertices; i++)
      if (tag[i] > 0)
        mesh->global_vtx_num[n_f_vertices++] = mesh->global_vtx_num[i];

    BFT_REALLOC(mesh->global_vtx_num, n_f_vertices, cs_gnum_t);

    vtx_io_num = fvm_io_num_create(NULL,
                                   mesh->global_vtx_num,
                                   mesh->n_vertices,
                                   0); /* Not shared */

    io_gnum = fvm_io_num_get_global_num(vtx_io_num);
    mesh->n_g_vertices = fvm_io_num_get_global_count(vtx_io_num);

    for (i = 0; i < mesh->n_vertices; i++)
      mesh->global_vtx_num[i] = io_gnum[i];

    fvm_io_num_destroy(vtx_io_num);

  }

  /* Compact vertex coordinates */

  n_f_vertices = 0;
  for (i = 0; i < n_i_vertices; i++) {
    if (tag[i] > 0) {

      tag[i] = n_f_vertices + 1;
      for (k = 0; k < 3; k++)
        mesh->vtx_coord[3*n_f_vertices+k] = mesh->vtx_coord[3*i+k];
      n_f_vertices++;

    }
  }

  BFT_REALLOC(mesh->vtx_coord, 3*mesh->n_vertices, cs_coord_t);

  /* Update interior face connectivity */

  for (i = 0; i < mesh->n_i_faces; i++) {
    for (j = mesh->i_face_vtx_idx[i];
         j < mesh->i_face_vtx_idx[i+1]; j++) {
      vid = mesh->i_face_vtx_lst[j];
      mesh->i_face_vtx_lst[j] = tag[vid] - 1;
    }
  }

  /* Update border face connectivity */

  for (i = 0; i < mesh->n_b_faces; i++) {
    for (j = mesh->b_face_vtx_idx[i];
         j < mesh->b_face_vtx_idx[i+1]; j++) {
      vid = mesh->b_face_vtx_lst[j];
      mesh->b_face_vtx_lst[j] = tag[vid] - 1;
    }
  }

  BFT_FREE(tag);

}

/*----------------------------------------------------------------------------
 * Define a new face connectivty without redundant edge definition.
 *
 * parameters:
 *   s       <-- starting index in f2v_lst
 *   e       <-- ending index in f2v_lst
 *   f2v_lst <-- face -> vertex connectivity list
 *   connect --> buffer to store locally the new face connectivity
 *   kill    --> buffer to store vertex to delete from the connectivity
 *
 * returns:
 *   new number of vertices in the face connectivity.
 *---------------------------------------------------------------------------*/

static cs_lnum_t
_delete_edges(cs_lnum_t        s,
              cs_lnum_t        e,
              const cs_lnum_t  f2v_lst[],
              cs_lnum_t        connect[],
              cs_lnum_t        kill[])
{
  cs_lnum_t  j, k, shift, count;

  cs_lnum_t  connect_size = e - s;

  /* Define local connectivity */

  for (j = s, k = 0; j < e; j++, k++)
    kill[k] = 0, connect[k] = f2v_lst[j] + 1;

  connect[k] = f2v_lst[s] + 1, kill[k++] = 0;
  connect[k] = f2v_lst[s+1] + 1, kill[k++] = 0;

  /* Find degenerate edges */

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

/*----------------------------------------------------------------------------
 * For faces with multiple parents, when some instances are adjacent
 * to a local cell and others are not, mark isolated instances
 * (i.e. copies on ranks not owning the local cell) as discarded.
 *
 * parameters:
 *   join_param      <-- set of parameters for the joining operation
 *   join_select     <-- list of all implied entities in the joining op.
 *   n2o_face_hist   <-- relation between faces before/after the joining
 *   join_mesh       <-> pointer to the local cs_join_mesh_t structure
 *   mesh            <-> pointer of pointer to cs_mesh_t structure
 *   cell_gnum    <->  global cell number related to each old face
 *   new_face_type    <-- type (border/interior) of new faces
 *
 * returns:
 *   local number of discarded faces
 *---------------------------------------------------------------------------*/

static cs_lnum_t
_discard_extra_multiple_b_faces(cs_join_param_t          join_param,
                                const cs_join_select_t  *join_select,
                                const cs_join_gset_t    *n2o_face_hist,
                                cs_join_mesh_t          *join_mesh,
                                cs_mesh_t               *mesh,
                                cs_gnum_t                cell_gnum[],
                                cs_join_face_type_t      new_face_type[])
{
  cs_lnum_t n_discard = 0;

  /* In periodic cases, no isolated faces have been selected,
     so this stage should not be necessary */

  if (join_param.perio_type != FVM_PERIODICITY_NULL)
    return n_discard;

  const int  rank = CS_MAX(cs_glob_rank_id, 0);
  const cs_gnum_t  rank_start = join_select->compact_rank_index[rank] + 1;
  const cs_gnum_t  rank_end = join_select->compact_rank_index[rank+1] + 1;

  for (cs_lnum_t i = 0; i < join_mesh->n_faces; i++) {

    if (new_face_type[i] == CS_JOIN_FACE_MULTIPLE_BORDER) {

      cs_lnum_t n_adj_cells = 0, c_id = -1;
      cs_gnum_t c_gnum = 0;

      for (cs_lnum_t j = n2o_face_hist->index[i];
           j < n2o_face_hist->index[i+1];
           j++) {

        cs_gnum_t c_g = cell_gnum[j];
        if (c_g > 0 && c_g != c_gnum) {
          c_gnum = c_g;
          n_adj_cells += 1;
        }

        cs_gnum_t compact_old_fgnum = n2o_face_hist->g_list[j];

        /* Initial selected border face must be in the selection */
        if (    rank_start <= compact_old_fgnum
             &&  compact_old_fgnum < rank_end) {
          cs_lnum_t fid = join_select->faces[compact_old_fgnum - rank_start] - 1;

          if (mesh->b_face_cells[fid] > -1)
            c_id = mesh->b_face_cells[fid];
        }

      }

      if (n_adj_cells == 1 && c_id < 0) {
        new_face_type[i] = CS_JOIN_FACE_DISCARD;
        n_discard += 1;
      }

    }

  }

  return n_discard;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Update mesh structure (vertices + faces) after the merge step.
 *
 * parameters:
 *   join_param   <-- set of parameters for the joining operation
 *   join_select  <-- list of all implied entities in the joining op.
 *   o2n_vtx_gnum <-> in : array on blocks on the new global vertex
 *                    out: local array on the new global vertex
 *   join_mesh    <-- pointer to the local cs_join_mesh_t structure
 *   mesh         <-> pointer of pointer to cs_mesh_t structure
 *---------------------------------------------------------------------------*/

void
cs_join_update_mesh_after_merge(cs_join_param_t        join_param,
                                cs_join_select_t      *join_select,
                                cs_gnum_t              o2n_vtx_gnum[],
                                cs_join_mesh_t        *join_mesh,
                                cs_mesh_t             *mesh)
{
  cs_lnum_t  i, j, select_id, adj_id, old_id, new_num;

  cs_lnum_t  *o2n_vtx_id = NULL, *join2mesh_vtx_id = NULL;
  edge_builder_t  *edge_builder = NULL;

  assert(join_select != NULL);
  assert(join_mesh != NULL);

  const cs_lnum_t  n_bm_vertices = mesh->n_vertices; /* bm: before merge */

  edge_builder = _init_edge_builder(join_select, mesh);

  /* Update mesh structure. Define new vertices */

  _update_vertices_after_merge
    (o2n_vtx_gnum,
     join_mesh,
     mesh,
     &join2mesh_vtx_id, /* size: join_mesh->n_vertices */
     &o2n_vtx_id);      /* size: n_bm_vertices */

  /* Define the evolution of each initial edge in edge_builder_t struct. */

  _complete_edge_builder(join_select,
                         join_mesh,
                         mesh,
                         o2n_vtx_gnum,
                         join2mesh_vtx_id,
                         edge_builder);

  BFT_FREE(o2n_vtx_gnum); /* Not useful after this point */

#if defined(HAVE_MPI)
  if (join_select->do_single_sync == true) {

    _sync_single_vertices(join_select,
                          o2n_vtx_id,
                          mesh);

    _sync_single_edges(join_select,
                       n_bm_vertices,
                       o2n_vtx_id,
                       join_mesh->n_vertices,
                       join2mesh_vtx_id,
                       edge_builder,
                       mesh);

  }
#endif

#if 0 && defined(DEBUG) && !defined(NDEBUG) /* Dump the structure */
  {
    int k;

    bft_printf("\n  Dump edge_builder_t structure (%p)\n", edge_builder);

    if (edge_builder != NULL) {
      bft_printf("  n_vertices: %9d\n"
                 "  n_edges   : %9d\n"
                 "  idx_size  : %9d\n",
                 edge_builder->n_vertices, edge_builder->n_edges,
                 edge_builder->v2v_idx[edge_builder->n_vertices]);

      for (i = 0; i < edge_builder->n_vertices; i++) {

        bft_printf("%7d - %9u- [%10.4f %10.4f %10.4f]: (%d, %d) v-v:",
                   i+1, mesh->global_vtx_num[i],
                   mesh->vtx_coord[3*i], mesh->vtx_coord[3*i+1],
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

  /* Update adjacent border face connectivity */

  _update_adj_face_connect(join_select->n_b_adj_faces,
                           join_select->b_adj_faces,
                           edge_builder,
                           o2n_vtx_id,
                           mesh->n_b_faces,
                           &(mesh->b_face_vtx_idx),
                           &(mesh->b_face_vtx_lst));

  /* Update adjacent interior face connectivity */

  _update_adj_face_connect(join_select->n_i_adj_faces,
                           join_select->i_adj_faces,
                           edge_builder,
                           o2n_vtx_id,
                           mesh->n_i_faces,
                           &(mesh->i_face_vtx_idx),
                           &(mesh->i_face_vtx_lst));

  /* Free memory */

  BFT_FREE(edge_builder->v2v_idx);
  BFT_FREE(edge_builder->v2v_lst);
  BFT_FREE(edge_builder->v2v_sub_idx);
  BFT_FREE(edge_builder->v2v_sub_lst);
  BFT_FREE(edge_builder);

  /* Update initial face connectivity for the remaining faces */

  for (i = 0, select_id = 0, adj_id = 0; i < mesh->n_b_faces; i++) {

    bool  do_update = true;

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

    if (do_update == true) {

      for (j = mesh->b_face_vtx_idx[i];
           j < mesh->b_face_vtx_idx[i+1]; j++) {

        old_id = mesh->b_face_vtx_lst[j];
        mesh->b_face_vtx_lst[j] = o2n_vtx_id[old_id];

      }

    }

  } /* End of loop on border faces */

  for (i = 0, adj_id = 0; i < mesh->n_i_faces; i++) {

    bool  do_update = true;

    if (adj_id < join_select->n_i_adj_faces) {
      if (i+1 == join_select->i_adj_faces[adj_id]) {
        do_update = false; /* Already done */
        adj_id++;
      }
    }

    if (do_update == true) {

      for (j = mesh->i_face_vtx_idx[i];
           j < mesh->i_face_vtx_idx[i+1]; j++) {

        old_id = mesh->i_face_vtx_lst[j];
        mesh->i_face_vtx_lst[j] = o2n_vtx_id[old_id];

      }

    }

  } /* End of loop on interior faces */

  /* Update the cs_join_select_t structure */

  for (i = 0; i < join_select->n_vertices; i++) {

    old_id = join_select->vertices[i] - 1;
    new_num = o2n_vtx_id[old_id] + 1;
    join_select->vertices[i] = new_num;

  }

  /* Update face state */

  _update_face_state(join_select, join_mesh, mesh, join2mesh_vtx_id);

  /* Post if required */

  if (join_param.visualization > 2)
    cs_join_post_after_merge(join_param, join_select);

  /* Free memory */

  BFT_FREE(join2mesh_vtx_id);
  BFT_FREE(o2n_vtx_id);
}

/*----------------------------------------------------------------------------
 * Update mesh structure (vertices + faces) after the face split step.
 *
 * parameters:
 *   join_param      <-- set of parameters for the joining operation
 *   join_select     <-- list of all implied entities in the joining op.
 *   n2o_face_hist   <-- relation between faces before/after the joining
 *   join_mesh       <-> pointer to the local cs_join_mesh_t structure
 *   mesh            <-> pointer of pointer to cs_mesh_t structure
 *   mesh_builder    <-> pointer of pointer to cs_mesh__builder_t structure
 *---------------------------------------------------------------------------*/

void
cs_join_update_mesh_after_split(cs_join_param_t          join_param,
                                const cs_join_select_t  *join_select,
                                const cs_join_gset_t    *n2o_face_hist,
                                cs_join_mesh_t          *join_mesh,
                                cs_mesh_t               *mesh,
                                cs_mesh_builder_t       *mesh_builder)
{
  int  i, j;

  cs_gnum_t  n_g_new_b_faces = 0, n_g_multiple_bfaces = 0;
  cs_lnum_t  n_new_i_faces = 0, n_new_b_faces = 0, n_undef_faces = 0;
  cs_lnum_t  n_multiple_bfaces = 0;
  cs_lnum_t  n_old_i_faces = mesh->n_i_faces, n_old_b_faces = mesh->n_b_faces;
  cs_lnum_t  *join2mesh_vtx_id = NULL;
  cs_gnum_t  *cell_gnum = NULL;
  int  *old_face_family = NULL;
  int  *new_face_family = NULL;
  cs_join_face_type_t  *new_face_type = NULL;
  FILE *logfile = cs_glob_join_log;

  const int  n_ranks = cs_glob_n_ranks;

  assert(mesh != NULL);
  assert(join_mesh != NULL);
  assert(mesh_builder != NULL);

  /* Get associated global cell number */

  _get_linked_cell_gnum_and_family(join_select, join_param, n2o_face_hist, mesh,
                                   &cell_gnum, &old_face_family);

  /* Get new subfaces evolution */

  assert(n2o_face_hist->n_elts == join_mesh->n_faces);

  BFT_MALLOC(new_face_type, join_mesh->n_faces, cs_join_face_type_t);

  for (i = 0; i < join_mesh->n_faces; i++) {

    assert(join_mesh->face_gnum[i] == n2o_face_hist->g_elts[i]);

    int n_parents = n2o_face_hist->index[i+1] - n2o_face_hist->index[i];

    if (n_parents == 1) {
      n_new_b_faces += 1;
      new_face_type[i] = CS_JOIN_FACE_BORDER;
    }
    else if (n_parents == 2) {
      j = n2o_face_hist->index[i];
      if (join_param.perio_type != FVM_PERIODICITY_NULL) {
        n_new_i_faces += 1;
        new_face_type[i] = CS_JOIN_FACE_INTERIOR;
      }
      else {
        cs_gnum_t g0 = cell_gnum[j], g1 = cell_gnum[j+1];
        if (g0 > 0 && g1 > 0 && g0 != g1)  {
          n_new_i_faces += 1;
          new_face_type[i] = CS_JOIN_FACE_INTERIOR;
        }
        else {
          n_new_b_faces += 1;
          new_face_type[i] = CS_JOIN_FACE_MULTIPLE_BORDER;
          n_multiple_bfaces += 1;
        }
      }
    }
    else if (n_parents > 2) {

      if (join_param.verbosity > 2) {
        fprintf(logfile,
                "  Warning: Face %d (%llu) has more than two ancestors.\n"
                "  Old faces implied:",
                i+1, (unsigned long long)join_mesh->face_gnum[i]);
        for (j = n2o_face_hist->index[i]; j < n2o_face_hist->index[i+1]; j++)
          fprintf(logfile, " %llu",
                  (unsigned long long)n2o_face_hist->g_list[j]);
        fprintf(logfile, "\n");
      }

      /* Border face by default */
      n_new_b_faces += 1;
      new_face_type[i] = CS_JOIN_FACE_MULTIPLE_BORDER;
      n_multiple_bfaces += 1;

    }
    else {

      n_undef_faces += 1;
      new_face_type[i] = CS_JOIN_FACE_UNDEFINED;
      n_multiple_bfaces += 1;

      if (logfile != NULL) {
        fprintf(logfile,
                "  Warning: Face %d (%llu) has no ancestor.\n",
                i+1, (unsigned long long)join_mesh->face_gnum[i]);
        for (j = n2o_face_hist->index[i]; j < n2o_face_hist->index[i+1]; j++)
          fprintf(logfile,
                  " %llu",
                  (unsigned long long)n2o_face_hist->g_list[j]);
        fprintf(logfile, "\n");
      }
    }

  } /* End of loop on faces */

  /* For faces with multiple parents, when some instances are adjacent
     to a local cell and others are not, discard isolated instances
     (i.e. copies on ranks not owning the local cell */

  if (n_multiple_bfaces > 1) {
    cs_lnum_t n_discard = _discard_extra_multiple_b_faces(join_param,
                                                          join_select,
                                                          n2o_face_hist,
                                                         join_mesh,
                                                          mesh,
                                                          cell_gnum,
                                                          new_face_type);
    n_multiple_bfaces -= n_discard;
    n_new_b_faces -= n_discard;
  }

  if (join_param.verbosity > 2)
    fprintf(logfile,
            "\n  Local configuration after the joining operation:\n"
            "    Number of interior faces to add: %9d\n"
            "    Number of border faces to add  : %9d\n",
            n_new_i_faces, n_new_b_faces);

  if (n_ranks == 1) {
    n_g_new_b_faces = n_new_b_faces;
    n_g_multiple_bfaces = n_multiple_bfaces;
  }

  cs_gnum_t n_g_undef = n_undef_faces;

#if defined(HAVE_MPI)
  if (n_ranks > 1) {
    cs_gnum_t _loc[3], _glob[3];
    MPI_Comm  mpi_comm = cs_glob_mpi_comm;

    _loc[0] = n_new_b_faces;
    _loc[1] = n_undef_faces;
    _loc[2] = n_multiple_bfaces;

    MPI_Allreduce(_loc, _glob, 3, CS_MPI_GNUM, MPI_SUM, mpi_comm);

    n_g_new_b_faces = _glob[0];
    n_g_undef = _glob[1];
    n_g_multiple_bfaces = _glob[2];
  }
#endif

  if (n_g_undef > 0)
    bft_error(__FILE__, __LINE__, 0,
              _(" There are %llu undefined faces with no ancestor.\n"
                " Cannot continue the joining algorithm.\n"),
              (unsigned long long)n_g_undef);

  /* Print information about new mesh configuration */

  if (cs_glob_mesh->verbosity)
    bft_printf(_("\n  Global configuration after the joining operation:\n"
                 "     Global number of border faces to add: %10llu\n"),
             (unsigned long long)n_g_new_b_faces);

  if (n_g_multiple_bfaces > 0) {

    cs_lnum_t *multiple_bfaces = NULL;

    if (cs_glob_mesh->verbosity)
      bft_printf(_("     Global number of multiple border faces: %10llu\n"),
                 (unsigned long long)n_g_multiple_bfaces);

    /* Post-processing of the implied faces */

    if (n_multiple_bfaces > 0) {

      BFT_MALLOC(multiple_bfaces, n_multiple_bfaces, cs_lnum_t);

      n_multiple_bfaces = 0;
      for (i = 0; i < join_mesh->n_faces; i++)
        if (new_face_type[i] == CS_JOIN_FACE_MULTIPLE_BORDER)
          multiple_bfaces[n_multiple_bfaces++] = i+1;

    }

    cs_join_post_faces_subset("MultipleBorderFaces",
                              join_mesh,
                              n_multiple_bfaces,
                              multiple_bfaces);

    BFT_FREE(multiple_bfaces);

  }

  /* Define join2mesh_vtx_id and add new vertices (only in parallel mode)
     These vertices already exist but only on other ranks. During the face
     splitting op., these vertices are used to define the new face connect. */

  _update_vertices_after_split(join_mesh, mesh, &join2mesh_vtx_id);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n  List of linked global cell number\n");
  for (i = 0; i < n2o_face_hist->index[n2o_face_hist->n_elts]; i++)
    bft_printf(" %9d - Old face: %10u - Cell gnum: %10u\n",
               i, n2o_face_hist->g_list[i], cell_gnum[i]);
  bft_printf_flush();
#endif

  new_face_family = _update_families(n2o_face_hist,
                                     old_face_family,
                                     join_mesh,
                                     mesh);

  BFT_FREE(old_face_family);

  /*  Update mesh structure:
        - Update first the interior faces because we need the global
          numbering of the initial border faces and its connectivity
        - Then update border faces  */

  _add_new_interior_faces(join_select,
                          join_param,
                          join_mesh,
                          join2mesh_vtx_id,
                          cell_gnum,
                          n_new_i_faces,
                          new_face_type,
                          new_face_family,
                          n2o_face_hist,
                          mesh);

  _add_new_border_faces(join_select,
                        join_param,
                        join_mesh,
                        join2mesh_vtx_id,
                        n_new_b_faces,
                        new_face_type,
                        new_face_family,
                        n2o_face_hist,
                        mesh);

  BFT_FREE(new_face_family);

  if (join_param.perio_type != FVM_PERIODICITY_NULL)
    cs_join_perio_split_update(join_param,
                               n_old_i_faces,
                               new_face_type,
                               join_mesh,
                               mesh,
                               mesh_builder);

  /* Delete unused vertices and define a compact global vertex numbering */

  _clean_vertices(join_param, mesh);

  /* Free memory */

  BFT_FREE(new_face_type);
  BFT_FREE(cell_gnum);
  BFT_FREE(join2mesh_vtx_id);

  /* Post if required */

  cs_join_post_after_split(n_old_i_faces,
                           n_old_b_faces,
                           n_g_new_b_faces,
                           join_select->n_faces,
                           mesh,
                           join_param);

  /* Check face -> vertex connectivity */

  for (i = 0; i < mesh->n_i_faces; i++) {

    for (j = mesh->i_face_vtx_idx[i];
         j < mesh->i_face_vtx_idx[i+1]; j++) {

       if (mesh->i_face_vtx_lst[j] < 0 ||
           mesh->i_face_vtx_lst[j] >= mesh->n_vertices)
         bft_error(__FILE__, __LINE__, 0,
                   "  Incoherency found in face -> vertex connect.\n"
                   "  for interior face %d (%llu)\n"
                   "  vtx: %d and n_vertices = %d\n",
                   i+1, (unsigned long long)(mesh->global_i_face_num[i]),
                   mesh->i_face_vtx_lst[j], mesh->n_vertices);

    }

  }

  for (i = 0; i < mesh->n_b_faces; i++) {

    for (j = mesh->b_face_vtx_idx[i];
         j < mesh->b_face_vtx_idx[i+1]; j++) {

       if (mesh->b_face_vtx_lst[j] < 0 ||
           mesh->b_face_vtx_lst[j] >= mesh->n_vertices)
         bft_error(__FILE__, __LINE__, 0,
                   "  Incoherency found in face -> vertex connect.\n"
                   "  for border face %d (%llu)\n"
                   "  vtx: %d and n_vertices = %d\n",
                   i+1, (unsigned long long)(mesh->global_b_face_num[i]),
                   mesh->b_face_vtx_lst[j], mesh->n_vertices);

    }

  }

}

/*----------------------------------------------------------------------------
 * Clean a cs_mesh_t struct.
 * Delete redundant and empty edge definitions.
 *
 * parameters:
 *   param   <-- set of parameters for the joining operation
 *   mesh    <-> pointer to a cs_mesh_t structure
 *---------------------------------------------------------------------------*/

void
cs_join_update_mesh_clean(cs_join_param_t   param,
                          cs_mesh_t        *mesh)
{
  assert(mesh != NULL);

  cs_lnum_t  i, j, s, e, n_vertices, n_init_vertices, connect_size;

  cs_lnum_t  connect_shift = 0;
  cs_lnum_t  max_connect = 0, b_size = 10, i_size = 10;
  cs_lnum_t  n_b_clean_faces = 0, n_i_clean_faces = 0;
  cs_gnum_t  n_g_clean_faces[2] = {0, 0};
  cs_lnum_t  *b_clean_faces = NULL, *i_clean_faces = NULL;
  cs_lnum_t  *kill = NULL, *connect = NULL;
  FILE *logfile = cs_glob_join_log;

  for (i = 0; i < mesh->n_b_faces; i++)
    max_connect = CS_MAX(max_connect,
                         mesh->b_face_vtx_idx[i+1] - mesh->b_face_vtx_idx[i]);

  for (i = 0; i < mesh->n_i_faces; i++)
    max_connect = CS_MAX(max_connect,
                         mesh->i_face_vtx_idx[i+1] - mesh->i_face_vtx_idx[i]);

  BFT_MALLOC(kill, max_connect + 2, cs_lnum_t);
  BFT_MALLOC(connect, max_connect + 2, cs_lnum_t);

  if (param.visualization > 1) {
    BFT_MALLOC(b_clean_faces, b_size, cs_lnum_t);
    BFT_MALLOC(i_clean_faces, i_size, cs_lnum_t);
  }

  /* Border face treatment */

  for (i = 0; i < mesh->n_b_faces; i++) {

    s = mesh->b_face_vtx_idx[i];
    e = mesh->b_face_vtx_idx[i+1];
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

      if (param.verbosity > 2)
        fprintf(logfile,
                "  Clean boundary face %d. New number of vertices: %d\n",
                i+1, n_vertices);

      if (param.visualization > 1) {
        if (n_b_clean_faces + 1 > b_size) {
          b_size *= 2;
          BFT_REALLOC(b_clean_faces, b_size, cs_lnum_t);
        }
        b_clean_faces[n_b_clean_faces] = i+1;
      }

      n_b_clean_faces++;

    }

    for (j = 0; j < n_vertices; j++)
      mesh->b_face_vtx_lst[connect_shift++] = connect[j] - 1;
    mesh->b_face_vtx_idx[i] = connect_shift;

  } /* End of loop on border faces */

  if (param.verbosity > 2)
    fprintf(logfile,
            "\n  Degenerate connectivity for %d final local boundary faces.\n",
            n_b_clean_faces);

  for (i = mesh->n_b_faces; i > 0; i--)
    mesh->b_face_vtx_idx[i] = mesh->b_face_vtx_idx[i-1];
  mesh->b_face_vtx_idx[0] = 0;

  BFT_REALLOC(mesh->b_face_vtx_lst, mesh->b_face_vtx_idx[mesh->n_b_faces],
              cs_lnum_t);

  /* Interior face treatment */

  connect_shift = 0;
  for (i = 0; i < mesh->n_i_faces; i++) {

    s = mesh->i_face_vtx_idx[i];
    e = mesh->i_face_vtx_idx[i+1];
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

      if (param.verbosity > 2)
        fprintf(logfile,
                "  Clean interior face %d. New number of vertices: %d\n",
                i+1, n_vertices);

      if (param.visualization > 1) {
        if (n_i_clean_faces + 1 > i_size) {
          i_size *= 2;
          BFT_REALLOC(i_clean_faces, i_size, cs_lnum_t);
        }
        i_clean_faces[n_i_clean_faces] = i+1;
      }

      n_i_clean_faces++;

    }

    for (j = 0; j < n_vertices; j++)
      mesh->i_face_vtx_lst[connect_shift++] = connect[j] - 1;
    mesh->i_face_vtx_idx[i] = connect_shift;

  } /* End of loop on interior faces */

  if (param.verbosity > 2)
    fprintf(logfile,
            "  Degenerate connectivity for %d final local interior faces.\n",
            n_i_clean_faces);

  for (i = mesh->n_i_faces; i > 0; i--)
    mesh->i_face_vtx_idx[i] = mesh->i_face_vtx_idx[i-1];
  mesh->i_face_vtx_idx[0] = 0;

  BFT_REALLOC(mesh->i_face_vtx_lst, mesh->i_face_vtx_idx[mesh->n_i_faces],
              cs_lnum_t);

  n_g_clean_faces[0] = n_i_clean_faces;
  n_g_clean_faces[1] = n_b_clean_faces;

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    cs_gnum_t  buf[2];
    MPI_Allreduce(n_g_clean_faces, buf, 2, CS_MPI_GNUM, MPI_SUM,
                  cs_glob_mpi_comm);
    n_g_clean_faces[0] = buf[0];
    n_g_clean_faces[1] = buf[1];
  }
#endif

  if (param.visualization > 1) { /* Post-process clean faces */

    if (n_g_clean_faces[0] > 0 || n_g_clean_faces[1] > 0) {

      BFT_REALLOC(i_clean_faces, n_i_clean_faces, cs_lnum_t);
      BFT_REALLOC(b_clean_faces, n_b_clean_faces, cs_lnum_t);

      cs_join_post_cleaned_faces(n_i_clean_faces,
                                 i_clean_faces,
                                 n_b_clean_faces,
                                 b_clean_faces,
                                 param);

    }

    BFT_FREE(b_clean_faces);
    BFT_FREE(i_clean_faces);

  } /* visualization > 1 */

  if (param.verbosity > 0) {
    bft_printf(_("\n  Mesh cleaning done for degenerate faces.\n"
                 "    Global number of cleaned interior faces: %8llu\n"
                 "    Global number of cleaned border faces:   %8llu\n"),
               (unsigned long long)n_g_clean_faces[0],
               (unsigned long long) n_g_clean_faces[1]);
    bft_printf_flush();
  }

  if (n_g_clean_faces[0] + n_g_clean_faces[1] > 0)
    mesh->modified = 1;

  /* Free memory */

  BFT_FREE(kill);
  BFT_FREE(connect);

}

/*---------------------------------------------------------------------------*/

END_C_DECLS
