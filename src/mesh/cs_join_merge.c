/*============================================================================
 * Set of subroutines for:
 *  - merging equivalent vertices,
 *  - managing tolerance reduction
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *---------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_io_num.h"

#include "cs_block_dist.h"
#include "cs_log.h"
#include "cs_order.h"
#include "cs_search.h"
#include "cs_join_post.h"
#include "cs_parall.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *---------------------------------------------------------------------------*/

#include "cs_join_merge.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local macro definitions
 *===========================================================================*/

/* Turn on (1) or off (0) the tolerance reduc. */
#define  CS_JOIN_MERGE_TOL_REDUC  1
#define  CS_JOIN_MERGE_INV_TOL  1

/*============================================================================
 * Local structure and type definitions
 *===========================================================================*/

/*============================================================================
 * Global variable definitions
 *===========================================================================*/

/* Parameters to control the vertex merge */

enum {

  CS_JOIN_MERGE_MAX_GLOB_ITERS = 50,  /* Max. number of glob. iter. for finding
                                         equivalent vertices */
  CS_JOIN_MERGE_MAX_LOC_ITERS = 100   /* Max. number of loc. iter. for finding
                                         equivalent vertices */
};

/* Coefficient to deal with rounding approximations */

static const double  cs_join_tol_eps_coef2 = 1.0001*1.001;

/* Counter on the number of loops useful to converge for the merge operation */

static int  _glob_merge_counter = 0, _loc_merge_counter = 0;

/*============================================================================
 * Private function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Initialize counter for the merge operation
 *---------------------------------------------------------------------------*/

static void
_initialize_merge_counter(void)
{
  _glob_merge_counter = 0;
  _loc_merge_counter = 0;
}

#if 0 && defined(DEBUG) && !defined(NDEBUG)

/*----------------------------------------------------------------------------
 * Dump an cs_join_eset_t structure on vertices.
 *
 * parameters:
 *   e_set   <-- cs_join_eset_t structure to dump
 *   mesh    <-- cs_join_mesh_t structure associated
 *   logfile <-- handle to log file
 *---------------------------------------------------------------------------*/

static void
_dump_vtx_eset(const cs_join_eset_t  *e_set,
               const cs_join_mesh_t  *mesh,
               FILE                  *logfile)
{
  int  i;

  fprintf(logfile, "\n  Dump an cs_join_eset_t structure (%p)\n",
          (const void *)e_set);
  fprintf(logfile, "  n_max_equiv: %10d\n", e_set->n_max_equiv);
  fprintf(logfile, "  n_equiv    : %10d\n\n", e_set->n_equiv);

  for (i = 0; i < e_set->n_equiv; i++) {

    cs_lnum_t  v1_num = e_set->equiv_couple[2*i];
    cs_lnum_t  v2_num = e_set->equiv_couple[2*i+1];

    fprintf(logfile,
            " %10d - local: (%9d, %9d) - global: (%10llu, %10llu)\n",
            i, v1_num, v2_num,
            (unsigned long long)(mesh->vertices[v1_num-1]).gnum,
            (unsigned long long)(mesh->vertices[v2_num-1]).gnum);

  }
  fflush(logfile);
}

#endif /* Only in debug mode */

/*----------------------------------------------------------------------------
 * Compute the length of a segment between two vertices.
 *
 * parameters:
 *   v1 <-- cs_join_vertex_t structure for the first vertex of the segment
 *   v2 <-- cs_join_vertex_t structure for the second vertex of the segment
 *
 * returns:
 *    length of the segment
 *---------------------------------------------------------------------------*/

inline static cs_real_t
_compute_length(cs_join_vertex_t  v1,
                cs_join_vertex_t  v2)
{
  cs_lnum_t  k;
  cs_real_t  len = 0.0, d2 = 0.0;

  for (k = 0; k < 3; k++) {
    cs_real_t  d = v1.coord[k] - v2.coord[k];
    d2 += d * d;
  }
  len = sqrt(d2);

  return len;
}

/*----------------------------------------------------------------------------
 * Compute a new cs_join_vertex_t structure.
 *
 * parameters:
 *   curv_abs   <-- curvilinear abscissa of the intersection
 *   gnum       <-- global number associated to the new
 *                  cs_join_vertex_t structure
 *   vtx_couple <-- couple of vertex numbers defining the current edge
 *   work       <-- local cs_join_mesh_t structure
 *
 * returns:
 *   a new cs_join_vertex_t structure
 *---------------------------------------------------------------------------*/

static cs_join_vertex_t
_get_new_vertex(cs_coord_t             curv_abs,
                cs_gnum_t              gnum,
                const cs_lnum_t        vtx_couple[],
                const cs_join_mesh_t  *work)
{
  cs_lnum_t  k;
  cs_join_vertex_t  new_vtx_data;

#if defined(DEBUG) && !defined(NDEBUG)
  /* Avoid Valgrind warnings in byte copies due to padding */
  memset(&new_vtx_data, 0, sizeof(cs_join_vertex_t));
#endif

  cs_join_vertex_t  v1 = work->vertices[vtx_couple[0]-1];
  cs_join_vertex_t  v2 = work->vertices[vtx_couple[1]-1];

  assert(curv_abs >= 0.0);
  assert(curv_abs <= 1.0);

  /* New vertex features */

  new_vtx_data.state = CS_JOIN_STATE_NEW;
  new_vtx_data.gnum = gnum;
  new_vtx_data.tolerance = (1-curv_abs)*v1.tolerance + curv_abs*v2.tolerance;

  for (k = 0; k < 3; k++)
    new_vtx_data.coord[k] = (1-curv_abs)*v1.coord[k] + curv_abs*v2.coord[k];

  return new_vtx_data;
}

/*----------------------------------------------------------------------------
 * Define a tag (3 values) to globally order intersections.
 *
 * parameters:
 *   tag           <-> tag to fill
 *   e1_gnum       <-- global number for the first edge
 *   e2_gnum       <-- global number for the second edge
 *   link_vtx_gnum <-- global number of the vertex associated to the current
 *                      intersection
 *---------------------------------------------------------------------------*/

static void
_define_inter_tag(cs_gnum_t  tag[],
                  cs_gnum_t  e1_gnum,
                  cs_gnum_t  e2_gnum,
                  cs_gnum_t  link_vtx_gnum)
{
  if (e1_gnum < e2_gnum) {
    tag[0] = e1_gnum;
    tag[1] = e2_gnum;
  }
  else {
    tag[0] = e2_gnum;
    tag[1] = e1_gnum;
  }

  tag[2] = link_vtx_gnum;
}

/*----------------------------------------------------------------------------
 * Creation of new vertices.
 *
 * Update list of equivalent vertices.
 *
 * parameters:
 *   work               <-- pointer to a cs_join_mesh_t structure
 *   edges              <-- list of edges
 *   inter_set          <-- structure including data on edge intersections
 *   init_max_vtx_gnum  <-- initial max. global numbering for vertices
 *   n_iwm_vertices     <-- initial local number of vertices (work struct)
 *   n_new_vertices     <-- local number of new vertices to define
 *   p_n_g_new_vertices <-> pointer to the global number of new vertices
 *   p_new_vtx_gnum     <-> pointer to the global numbering array for the
 *                          new vertices
 *---------------------------------------------------------------------------*/

static void
_compute_new_vertex_gnum(const cs_join_mesh_t       *work,
                         const cs_join_edges_t      *edges,
                         const cs_join_inter_set_t  *inter_set,
                         cs_gnum_t                   init_max_vtx_gnum,
                         cs_lnum_t                   n_iwm_vertices,
                         cs_lnum_t                   n_new_vertices,
                         cs_gnum_t                  *p_n_g_new_vertices,
                         cs_gnum_t                  *p_new_vtx_gnum[])
{
  cs_lnum_t  i;

  cs_gnum_t  n_g_new_vertices = 0;
  cs_lnum_t  n_new_vertices_save = n_new_vertices;
  cs_lnum_t  *order = NULL;
  cs_gnum_t  *inter_tag = NULL, *adjacency = NULL, *new_vtx_gnum = NULL;
  fvm_io_num_t  *new_vtx_io_num = NULL;

  /* Define a fvm_io_num_t structure to get the global numbering
     for the new vertices.
     First, build a tag associated to each intersection */

  BFT_MALLOC(new_vtx_gnum, n_new_vertices, cs_gnum_t);
  BFT_MALLOC(inter_tag, 3*n_new_vertices, cs_gnum_t);

  n_new_vertices = 0;

  for (i = 0; i < inter_set->n_inter; i++) {

    cs_join_inter_t  inter1 = inter_set->inter_lst[2*i];
    cs_join_inter_t  inter2 = inter_set->inter_lst[2*i+1];
    cs_gnum_t  e1_gnum = edges->gnum[inter1.edge_id];
    cs_gnum_t  e2_gnum = edges->gnum[inter2.edge_id];

    if (inter1.vtx_id + 1 > n_iwm_vertices) {

      if (inter2.vtx_id + 1 > n_iwm_vertices)
        _define_inter_tag(&(inter_tag[3*n_new_vertices]),
                          e1_gnum, e2_gnum,
                          0);
      else
        _define_inter_tag(&(inter_tag[3*n_new_vertices]),
                          e1_gnum, e2_gnum,
                          (work->vertices[inter2.vtx_id]).gnum);

      n_new_vertices++;

    } /* New vertices for this intersection */

    if (inter2.vtx_id + 1 > n_iwm_vertices) {

      if (inter1.vtx_id + 1 > n_iwm_vertices)
        _define_inter_tag(&(inter_tag[3*n_new_vertices]),
                          e1_gnum, e2_gnum,
                          init_max_vtx_gnum + 1);
      else
        _define_inter_tag(&(inter_tag[3*n_new_vertices]),
                          e1_gnum, e2_gnum,
                          (work->vertices[inter1.vtx_id]).gnum);

      n_new_vertices++;

    } /* New vertices for this intersection */

  } /* End of loop on intersections */

  if (n_new_vertices != n_new_vertices_save)
    bft_error(__FILE__, __LINE__, 0,
              _("  The number of new vertices to create is not consistent.\n"
                "     Previous number: %10d\n"
                "     Current number:  %10d\n\n"),
              n_new_vertices_save, n_new_vertices);

  /* Create a new fvm_io_num_t structure */

  BFT_MALLOC(order, n_new_vertices, cs_lnum_t);

  cs_order_gnum_allocated_s(NULL, inter_tag, 3, order, n_new_vertices);

  BFT_MALLOC(adjacency, 3*n_new_vertices, cs_gnum_t);

  for (i = 0; i < n_new_vertices; i++) {

    cs_lnum_t  o_id = order[i];

    adjacency[3*i] = inter_tag[3*o_id];
    adjacency[3*i+1] = inter_tag[3*o_id+1];
    adjacency[3*i+2] = inter_tag[3*o_id+2];

  }

  BFT_FREE(inter_tag);

  if (cs_glob_n_ranks > 1) {

    const cs_gnum_t  *global_num = NULL;

    new_vtx_io_num =
      fvm_io_num_create_from_adj_s(NULL, adjacency, n_new_vertices, 3);

    n_g_new_vertices = fvm_io_num_get_global_count(new_vtx_io_num);
    global_num = fvm_io_num_get_global_num(new_vtx_io_num);

    for (i = 0; i < n_new_vertices; i++)
      new_vtx_gnum[order[i]] = global_num[i] + init_max_vtx_gnum;

    fvm_io_num_destroy(new_vtx_io_num);

  } /* End of parallel treatment */

  else {

    if (n_new_vertices > 0) {

      cs_gnum_t  new_gnum = init_max_vtx_gnum + 1;

      new_vtx_gnum[order[0]] = new_gnum;

      for (i = 1; i < n_new_vertices; i++) {

        if (adjacency[3*i] != adjacency[3*(i-1)])
          new_gnum += 1;
        else {
          if (adjacency[3*i+1] != adjacency[3*(i-1)+1])
            new_gnum += 1;
          else
            if (adjacency[3*i+2] != adjacency[3*(i-1)+2])
              new_gnum += 1;
        }

        new_vtx_gnum[order[i]] = new_gnum;

      }

    } /* End if n_new_vertices > 0 */

    n_g_new_vertices = n_new_vertices;

  } /* End of serial treatment */

  /* Free memory */

  BFT_FREE(order);
  BFT_FREE(adjacency);

  /* Return pointer */

  *p_n_g_new_vertices = n_g_new_vertices;
  *p_new_vtx_gnum = new_vtx_gnum;

}

/*----------------------------------------------------------------------------
 * Get vertex id associated to the current intersection.
 *
 * Create a new vertex id if needed. Update n_new_vertices in this case.
 *
 * parameters:
 *   inter           <-- a inter_t structure
 *   vtx_couple      <-- couple of vertex numbers defining the current edge
 *   n_init_vertices <-- initial number of vertices
 *   n_new_vertices  <-- number of new vertices created
 *
 * returns:
 *   vertex id associated to the current intersection.
 *---------------------------------------------------------------------------*/

static cs_lnum_t
_get_vtx_id(cs_join_inter_t  inter,
            const cs_lnum_t  vtx_couple[],
            cs_lnum_t        n_init_vertices,
            cs_lnum_t       *p_n_new_vertices)
{
  cs_lnum_t  vtx_id = -1;
  cs_lnum_t  n_new_vertices = *p_n_new_vertices;

  assert(inter.curv_abs >= 0.0);
  assert(inter.curv_abs <= 1.0);

  if (inter.curv_abs <= 0.0)
    vtx_id = vtx_couple[0] - 1;

  else if (inter.curv_abs >= 1.0)
    vtx_id = vtx_couple[1] - 1;

  else {

    assert(inter.curv_abs > 0 && inter.curv_abs < 1.0);
    vtx_id = n_init_vertices + n_new_vertices;
    n_new_vertices++;

  }

  assert(vtx_id != -1);

  *p_n_new_vertices = n_new_vertices;

  return vtx_id;
}


/*----------------------------------------------------------------------------
 * Test if we have to continue to spread the tag associate to each vertex
 *
 * parameters:
 *   n_vertices   <-- local number of vertices
 *   prev_vtx_tag <-- previous tag for each vertex
 *   vtx_tag      <-- tag for each vertex
 *
 * returns:
 *   true or false
 *---------------------------------------------------------------------------*/

static bool
_is_spread_not_converged(cs_lnum_t        n_vertices,
                         const cs_gnum_t  prev_vtx_tag[],
                         const cs_gnum_t  vtx_tag[])
{
  cs_lnum_t  i;

  bool  have_to_continue = true;

  for (i = 0; i < n_vertices; i++)
    if (vtx_tag[i] != prev_vtx_tag[i])
      break;

  if (i == n_vertices)
    have_to_continue = false;

  return have_to_continue;
}

/*----------------------------------------------------------------------------
 * Spread the tag associated to each vertex according the rule:
 *  Between two equivalent vertices, the tag associated to each considered
 *  vertex is equal to the minimal global number.
 *
 * parameters:
 *  n_vertices <-- local number of vertices
 *  vtx_eset   <-- structure dealing with vertices equivalences
 *  vtx_tag    <-> tag for each vertex
 *---------------------------------------------------------------------------*/

static void
_spread_tag(cs_lnum_t              n_vertices,
            const cs_join_eset_t  *vtx_eset,
            cs_gnum_t              vtx_tag[])
{
  cs_lnum_t  i, v1_id, v2_id;
  cs_gnum_t  v1_gnum, v2_gnum;
  cs_lnum_t  *equiv_lst = vtx_eset->equiv_couple;

  for (i = 0; i < vtx_eset->n_equiv; i++) {

    v1_id = equiv_lst[2*i] - 1, v2_id = equiv_lst[2*i+1] - 1;
    assert(v1_id < n_vertices);
    assert(v1_id < n_vertices);
    v1_gnum = vtx_tag[v1_id], v2_gnum = vtx_tag[v2_id];

    if (v1_gnum != v2_gnum) {

      cs_gnum_t  min_gnum = CS_MIN(v1_gnum, v2_gnum);

      vtx_tag[v1_id] = min_gnum;
      vtx_tag[v2_id] = min_gnum;
    }

  } /* End of loop on vertex equivalences */
}

/*----------------------------------------------------------------------------
 * Define an array wich keeps the new vertex id of each vertex.
 *
 * If two vertices have the same vertex id, they should merge.
 *
 * parameters:
 *   vtx_eset     <-- structure dealing with vertex equivalences
 *   n_vertices   <-- local number of vertices
 *   prev_vtx_tag <-> previous tag for each vertex
 *   vtx_tag      <-> tag for each vertex
 *---------------------------------------------------------------------------*/

static void
_local_spread(const cs_join_eset_t  *vtx_eset,
              cs_lnum_t              n_vertices,
              cs_gnum_t              prev_vtx_tag[],
              cs_gnum_t              vtx_tag[])
{
  cs_lnum_t  i;

  _loc_merge_counter++;

  _spread_tag(n_vertices, vtx_eset, vtx_tag);

  while (_is_spread_not_converged(n_vertices, prev_vtx_tag, vtx_tag)) {

    _loc_merge_counter++;

    if (_loc_merge_counter > CS_JOIN_MERGE_MAX_LOC_ITERS)
      bft_error(__FILE__, __LINE__, 0,
                _("\n  The authorized maximum number of iterations "
                  " for the merge of vertices has been reached.\n"
                  "  Local counter on iteration : %d (MAX =%d)\n"
                  "  Check the fraction parameter.\n"),
                _loc_merge_counter, CS_JOIN_MERGE_MAX_LOC_ITERS);

    for (i = 0; i < n_vertices; i++)
      prev_vtx_tag[i] = vtx_tag[i];

    _spread_tag(n_vertices, vtx_eset, vtx_tag);
  }
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Exchange local vtx_tag buffer over the ranks and update global vtx_tag
 * buffers. Apply modifications observed on the global vtx_tag to the local
 * vtx_tag.
 *
 * parameters:
 *   block_size        <-- size of block for the current rank
 *   work              <-- local cs_join_mesh_t structure which has initial
 *                         vertex data
 *   vtx_tag           <-> local vtx_tag for the local vertices
 *   glob_vtx_tag      <-> global vtx_tag affected to the local rank
 *                         (size: block_size)
 *   prev_glob_vtx_tag <-> same but for the previous iteration
 *   recv2glob         <-> buffer used to place correctly receive elements
 *   send_count        <-> buffer used to count the number of elts to send
 *   send_shift        <-> index on ranks of the elements to send
 *   send_glob_buffer  <-> buffer used to save elements to send
 *   recv_count        <-> buffer used to count the number of elts to receive
 *   recv_shift        <-> index on ranks of the elements to receive
 *   recv_glob_buffer  <-> buffer used to save elements to receive
 *
 * returns:
 *   true if we have to continue the spread, false otherwise.
 *---------------------------------------------------------------------------*/

static bool
_global_spread(cs_lnum_t              block_size,
               const cs_join_mesh_t  *work,
               cs_gnum_t              vtx_tag[],
               cs_gnum_t              glob_vtx_tag[],
               cs_gnum_t              prev_glob_vtx_tag[],
               cs_gnum_t              recv2glob[],
               int                    send_count[],
               int                    send_shift[],
               cs_gnum_t              send_glob_buffer[],
               int                    recv_count[],
               int                    recv_shift[],
               cs_gnum_t              recv_glob_buffer[])
{
  bool  ret_value;
  cs_lnum_t  i, local_value, global_value;

  cs_lnum_t  n_vertices = work->n_vertices;
  int  n_ranks = cs_glob_n_ranks;
  cs_gnum_t  _n_ranks = cs_glob_n_ranks;
  MPI_Comm  mpi_comm = cs_glob_mpi_comm;

  _glob_merge_counter++;

  /* Push modifications in local vtx_tag to the global vtx_tag */

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (i = 0; i < n_vertices; i++) {

    int  rank = (work->vertices[i].gnum - 1) % _n_ranks;
    cs_lnum_t  shift = send_shift[rank] + send_count[rank];

    send_glob_buffer[shift] = vtx_tag[i];
    send_count[rank] += 1;

  }

  MPI_Alltoallv(send_glob_buffer, send_count, send_shift, CS_MPI_GNUM,
                recv_glob_buffer, recv_count, recv_shift, CS_MPI_GNUM,
                mpi_comm);

  /* Apply update to glob_vtx_tag */

  for (i = 0; i < recv_shift[n_ranks]; i++) {
    cs_lnum_t  cur_id = recv2glob[i];
    glob_vtx_tag[cur_id] = CS_MIN(glob_vtx_tag[cur_id], recv_glob_buffer[i]);
  }

  ret_value = _is_spread_not_converged(block_size,
                                       prev_glob_vtx_tag,
                                       glob_vtx_tag);

  if (ret_value == false)
    local_value = 0;
  else
    local_value = 1;

  MPI_Allreduce(&local_value, &global_value, 1, MPI_INT, MPI_SUM, mpi_comm);

  if (global_value > 0) { /* Store the current state as the previous one
                             Update local vtx_tag */

    if (_glob_merge_counter > CS_JOIN_MERGE_MAX_GLOB_ITERS)
      bft_error(__FILE__, __LINE__, 0,
                _("\n  The authorized maximum number of iterations "
                  " for the merge of vertices has been reached.\n"
                  "  Global counter on iteration : %d (MAX =%d)\n"
                  "  Check the fraction parameter.\n"),
                _glob_merge_counter, CS_JOIN_MERGE_MAX_GLOB_ITERS);

    for (i = 0; i < block_size; i++)
      prev_glob_vtx_tag[i] = glob_vtx_tag[i];

    for (i = 0; i < recv_shift[n_ranks]; i++)
      recv_glob_buffer[i] = glob_vtx_tag[recv2glob[i]];

    MPI_Alltoallv(recv_glob_buffer, recv_count, recv_shift, CS_MPI_GNUM,
                  send_glob_buffer, send_count, send_shift, CS_MPI_GNUM,
                  mpi_comm);

    /* Update vtx_tag */

    for (i = 0; i < n_ranks; i++)
      send_count[i] = 0;

    assert(send_shift[n_ranks] == n_vertices);

    for (i = 0; i < n_vertices; i++) {

      int  rank = (work->vertices[i].gnum - 1) % _n_ranks;
      cs_lnum_t  shift = send_shift[rank] + send_count[rank];

      vtx_tag[i] = CS_MIN(send_glob_buffer[shift], vtx_tag[i]);
      send_count[rank] += 1;

    }

    return true;

  } /* End if prev_glob_vtx_tag != glob_vtx_tag */

  else
    return false; /* No need to continue */
}

/*----------------------------------------------------------------------------
 * Initialize and allocate buffers for the tag operation in parallel mode.
 *
 * parameters:
 *   n_g_vertices_to_treat <-- global number of vertices to consider for the
 *                             merge operation (existing + created vertices)
 *   work                  <-- local cs_join_mesh_t structure which has
 *                             initial vertex data
 *   p_block_size          <-> size of block for the current rank
 *   p_send_count          <-> buf. for counting the number of elts to send
 *   p_send_shift          <-> index on ranks of the elements to send
 *   p_send_glob_buf       <-> buf. for saving elements to send
 *   p_recv_count          <-> buf. for counting the number of elts to receive
 *   p_recv_shift          <-> index on ranks of the elements to receive
 *   p_recv_glob_buf       <-> buf. for storing elements to receive
 *   p_recv2glob           <-> buf. for putting correctly received elements
 *   p_glob_vtx_tag        <-> vtx_tag locally treated (size = block_size)
 *   p_prev_glob_vtx_tag   <-> idem but for the previous iteration
 *---------------------------------------------------------------------------*/

static void
_parall_tag_init(cs_gnum_t              n_g_vertices_to_treat,
                 const cs_join_mesh_t  *work,
                 cs_lnum_t             *p_block_size,
                 int                   *p_send_count[],
                 int                   *p_send_shift[],
                 cs_gnum_t             *p_send_glob_buf[],
                 int                   *p_recv_count[],
                 int                   *p_recv_shift[],
                 cs_gnum_t             *p_recv_glob_buf[],
                 cs_gnum_t             *p_recv2glob[],
                 cs_gnum_t             *p_glob_vtx_tag[],
                 cs_gnum_t             *p_prev_glob_vtx_tag[])
{
  cs_lnum_t  i;

  cs_lnum_t  n_vertices = work->n_vertices;
  cs_lnum_t  block_size = 0, left_over = 0;
  int        *send_count = NULL, *recv_count = NULL;
  int        *send_shift = NULL, *recv_shift = NULL;
  cs_gnum_t  *recv2glob = NULL;
  cs_gnum_t  *recv_glob_buffer = NULL, *send_glob_buffer = NULL;
  cs_gnum_t  *glob_vtx_tag = NULL, *prev_glob_vtx_tag = NULL;
  MPI_Comm  mpi_comm = cs_glob_mpi_comm;

  const int  n_ranks = cs_glob_n_ranks;
  const int  local_rank = CS_MAX(cs_glob_rank_id, 0);
  const cs_gnum_t  _n_ranks = n_ranks, _local_rank = local_rank;

  /* Allocate and intialize vtx_tag associated to the local rank */

  block_size = n_g_vertices_to_treat / _n_ranks;
  left_over = n_g_vertices_to_treat % _n_ranks;
  if (local_rank < left_over)
    block_size += 1;

  BFT_MALLOC(glob_vtx_tag, block_size, cs_gnum_t);
  BFT_MALLOC(prev_glob_vtx_tag, block_size, cs_gnum_t);

  for (i = 0; i < block_size; i++) {
    cs_gnum_t gi = i;
    prev_glob_vtx_tag[i] = gi*_n_ranks + _local_rank + 1;
    glob_vtx_tag[i] = gi*_n_ranks + _local_rank + 1;
  }

  /* Allocate and define send/recv_count/shift */

  BFT_MALLOC(send_count, n_ranks, int);
  BFT_MALLOC(recv_count, n_ranks, int);
  BFT_MALLOC(send_shift, n_ranks + 1, int);
  BFT_MALLOC(recv_shift, n_ranks + 1, int);

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (i = 0; i < n_vertices; i++) {
    int  rank = (work->vertices[i].gnum - 1) % _n_ranks;
    send_count[rank] += 1;
  }

  MPI_Alltoall(send_count, 1, MPI_INT,
               recv_count, 1, MPI_INT, mpi_comm);

  /* Build index */

  for (i = 0; i < n_ranks; i++) {
    send_shift[i+1] = send_shift[i] + send_count[i];
    recv_shift[i+1] = recv_shift[i] + recv_count[i];
  }

  assert(send_shift[n_ranks] == n_vertices);

  /* Allocate and define recv2glob */

  BFT_MALLOC(send_glob_buffer, send_shift[n_ranks], cs_gnum_t);
  BFT_MALLOC(recv2glob, recv_shift[n_ranks], cs_gnum_t);
  BFT_MALLOC(recv_glob_buffer, recv_shift[n_ranks], cs_gnum_t);

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (i = 0; i < n_vertices; i++) {

    int  rank = (work->vertices[i].gnum - 1) % _n_ranks;
    cs_lnum_t  shift = send_shift[rank] + send_count[rank];

    send_glob_buffer[shift] = (work->vertices[i].gnum - 1) / _n_ranks;
    send_count[rank] += 1;

  }

  MPI_Alltoallv(send_glob_buffer, send_count, send_shift, CS_MPI_GNUM,
                recv2glob, recv_count, recv_shift, CS_MPI_GNUM,
                mpi_comm);

  /* Return pointers */

  *p_block_size = block_size;
  *p_send_count = send_count;
  *p_recv_count = recv_count;
  *p_send_shift = send_shift;
  *p_recv_shift = recv_shift;
  *p_send_glob_buf = send_glob_buffer;
  *p_recv_glob_buf = recv_glob_buffer;
  *p_recv2glob = recv2glob;
  *p_glob_vtx_tag = glob_vtx_tag;
  *p_prev_glob_vtx_tag = prev_glob_vtx_tag;

}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Tag with the same number all the vertices which might be merged together
 *
 * parameters:
 *   n_g_vertices_tot <-- global number of vertices to consider for the
 *                        merge operation (existing + created vertices)
 *   vtx_eset         <-- structure dealing with vertex equivalences
 *   work             <-- local cs_join_mesh_t structure which has initial
 *                        vertex data
 *   verbosity        <-- level of detail in information display
 *   p_vtx_tag        --> pointer to the vtx_tag for the local vertices
 *---------------------------------------------------------------------------*/

static void
_tag_equiv_vertices(cs_gnum_t              n_g_vertices_tot,
                    const cs_join_eset_t  *vtx_eset,
                    const cs_join_mesh_t  *work,
                    int                    verbosity,
                    cs_gnum_t             *p_vtx_tag[])
{
  cs_lnum_t  i;

  cs_gnum_t  *vtx_tag = NULL;
  cs_gnum_t  *prev_vtx_tag = NULL;
  FILE  *logfile = cs_glob_join_log;

  const cs_lnum_t  n_vertices = work->n_vertices;
  const int  n_ranks = cs_glob_n_ranks;

  /* Local initialization : we tag each vertex by its global number */

  BFT_MALLOC(prev_vtx_tag, n_vertices, cs_gnum_t);
  BFT_MALLOC(vtx_tag, n_vertices, cs_gnum_t);

  for (i = 0; i < work->n_vertices; i++) {

    cs_gnum_t  v_gnum = work->vertices[i].gnum;

    vtx_tag[i] = v_gnum;
    prev_vtx_tag[i] = v_gnum;

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  for (i = 0; i < n_vertices; i++)
    fprintf(logfile, " Initial vtx_tag[%6d] = %9llu\n",
            i, (unsigned long long)vtx_tag[i]);
  fflush(logfile);
#endif

  /* Compute vtx_tag */

  _local_spread(vtx_eset, n_vertices, prev_vtx_tag, vtx_tag);

  if (n_ranks > 1) { /* Parallel treatment */

    bool  go_on;

    cs_lnum_t  block_size = 0;
    int        *send_count = NULL, *recv_count = NULL;
    int        *send_shift = NULL, *recv_shift = NULL;
    cs_gnum_t  *recv2glob = NULL;
    cs_gnum_t  *recv_glob_buffer = NULL, *send_glob_buffer = NULL;
    cs_gnum_t  *glob_vtx_tag = NULL, *prev_glob_vtx_tag = NULL;

#if defined(HAVE_MPI)
    _parall_tag_init(n_g_vertices_tot,
                     work,
                     &block_size,
                     &send_count, &send_shift, &send_glob_buffer,
                     &recv_count, &recv_shift, &recv_glob_buffer,
                     &recv2glob,
                     &glob_vtx_tag,
                     &prev_glob_vtx_tag);

    go_on = _global_spread(block_size,
                           work,
                           vtx_tag,
                           glob_vtx_tag,
                           prev_glob_vtx_tag,
                           recv2glob,
                           send_count, send_shift, send_glob_buffer,
                           recv_count, recv_shift, recv_glob_buffer);

    while (go_on == true) {

      /* Local convergence of vtx_tag */

      _local_spread(vtx_eset, n_vertices, prev_vtx_tag, vtx_tag);

      /* Global update and test to continue */

      go_on = _global_spread(block_size,
                             work,
                             vtx_tag,
                             glob_vtx_tag,
                             prev_glob_vtx_tag,
                             recv2glob,
                             send_count, send_shift, send_glob_buffer,
                             recv_count, recv_shift, recv_glob_buffer);

    }

    /* Partial free */

    BFT_FREE(glob_vtx_tag);
    BFT_FREE(prev_glob_vtx_tag);
    BFT_FREE(send_count);
    BFT_FREE(send_shift);
    BFT_FREE(send_glob_buffer);
    BFT_FREE(recv_count);
    BFT_FREE(recv_shift);
    BFT_FREE(recv2glob);
    BFT_FREE(recv_glob_buffer);

#endif
  } /* End of parallel treatment */

  BFT_FREE(prev_vtx_tag);

  if (verbosity > 3) {
    fprintf(logfile,
            "\n  Number of local iterations to converge on vertex"
            " equivalences: %3d\n", _loc_merge_counter);
    if (n_ranks > 1)
      fprintf(logfile,
              "  Number of global iterations to converge on vertex"
              " equivalences: %3d\n\n", _glob_merge_counter);
    fflush(logfile);
  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (logfile != NULL) {
    for (i = 0; i < n_vertices; i++)
      fprintf(logfile, " Final vtx_tag[%6d] = %9llu\n",
              i, (unsigned long long)vtx_tag[i]);
    fflush(logfile);
  }
#endif

  /* Returns pointer */

  *p_vtx_tag = vtx_tag;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Build in parallel a cs_join_gset_t structure to store all the potential
 * merges between vertices and its associated cs_join_vertex_t structure.
 *
 * parameters:
 *   work             <-- local cs_join_mesh_t structure which
 *                        has initial vertex data
 *   vtx_tag          <-- local vtx_tag for the local vertices
 *   send_count       <-> buffer used to count the number of elts to send
 *   send_shift       <-> index on ranks of the elements to send
 *   recv_count       <-> buffer used to count the number of elts to receive
 *   recv_shift       <-> index on ranks of the elements to receive
 *   p_vtx_merge_data <-> a pointer to a cs_join_vertex_t structure which
 *                        stores data about merged vertices
 *   p_merge_set      <-> pointer to a cs_join_gset_t struct. storing the
 *                        evolution of each global vtx number
 *---------------------------------------------------------------------------*/

static void
_build_parall_merge_structures(const cs_join_mesh_t    *work,
                               const cs_gnum_t          vtx_tag[],
                               int                      send_count[],
                               int                      send_shift[],
                               int                      recv_count[],
                               int                      recv_shift[],
                               cs_join_vertex_t        *p_vtx_merge_data[],
                               cs_join_gset_t         **p_merge_set)
{
  cs_lnum_t  i;

  cs_lnum_t  n_vertices = work->n_vertices;
  cs_gnum_t   *recv_gbuf = NULL, *send_gbuf = NULL;
  cs_join_vertex_t  *send_vtx_data = NULL, *recv_vtx_data = NULL;
  cs_join_gset_t  *merge_set = NULL;

  MPI_Datatype  CS_MPI_JOIN_VERTEX = cs_join_mesh_create_vtx_datatype();
  MPI_Comm  mpi_comm = cs_glob_mpi_comm;

  const int  n_ranks = cs_glob_n_ranks;
  const cs_gnum_t  _n_ranks = n_ranks;

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (i = 0; i < n_vertices; i++) {
    int  rank = (vtx_tag[i] - 1) % _n_ranks;
    send_count[rank] += 1;
  }

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);

  /* Build index */

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (i = 0; i < n_ranks; i++) {
    send_shift[i+1] = send_shift[i] + send_count[i];
    recv_shift[i+1] = recv_shift[i] + recv_count[i];
  }

  assert(send_shift[n_ranks] == n_vertices);

  /* Allocate and define recv_gbuf and send_gbuf */

  BFT_MALLOC(send_gbuf, send_shift[n_ranks], cs_gnum_t);
  BFT_MALLOC(recv_gbuf, recv_shift[n_ranks], cs_gnum_t);

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (i = 0; i < n_vertices; i++) {

    int  rank = (vtx_tag[i] - 1) % _n_ranks;
    cs_lnum_t  shift = send_shift[rank] + send_count[rank];

    send_gbuf[shift] = vtx_tag[i];
    send_count[rank] += 1;

  }

  MPI_Alltoallv(send_gbuf, send_count, send_shift, CS_MPI_GNUM,
                recv_gbuf, recv_count, recv_shift, CS_MPI_GNUM,
                mpi_comm);

  /* Allocate and build send_vtx_data, receive recv_vtx_data. */

  BFT_MALLOC(recv_vtx_data, recv_shift[n_ranks], cs_join_vertex_t);
  BFT_MALLOC(send_vtx_data, send_shift[n_ranks], cs_join_vertex_t);

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (i = 0; i < n_vertices; i++) {

    int  rank = (vtx_tag[i] - 1) % _n_ranks;
    cs_lnum_t  shift = send_shift[rank] + send_count[rank];

    send_vtx_data[shift] = work->vertices[i];
    send_count[rank] += 1;

  }

  MPI_Alltoallv(send_vtx_data, send_count, send_shift, CS_MPI_JOIN_VERTEX,
                recv_vtx_data, recv_count, recv_shift, CS_MPI_JOIN_VERTEX,
                mpi_comm);

  /* Partial free memory */

  BFT_FREE(send_vtx_data);
  BFT_FREE(send_gbuf);
  MPI_Type_free(&CS_MPI_JOIN_VERTEX);

  /* Build merge set */

  merge_set = cs_join_gset_create_from_tag(recv_shift[n_ranks], recv_gbuf);

  cs_join_gset_sort_sublist(merge_set);

  /* Free memory */

  BFT_FREE(recv_gbuf);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (cs_glob_join_log != NULL) {
    FILE *logfile = cs_glob_join_log;
    fprintf(logfile,
            "\n  Number of vertices to treat for the merge step: %d\n",
            recv_shift[n_ranks]);
    fprintf(logfile,
            "  List of vertices to treat:\n");
    for (i = 0; i < recv_shift[n_ranks]; i++) {
      fprintf(logfile, " %9d - ", i);
      cs_join_mesh_dump_vertex(logfile, recv_vtx_data[i]);
    }
    fflush(logfile);
  }
#endif

  /* Set return pointers */

  *p_merge_set = merge_set;
  *p_vtx_merge_data = recv_vtx_data;
}

/*----------------------------------------------------------------------------
 * Exchange the updated cs_join_vertex_t array over the ranks.
 *
 * parameters:
 *   work           <-- local cs_join_mesh_t structure which
 *                      has initial vertex data
 *   vtx_tag        <-- local vtx_tag for the local vertices
 *   send_count     <-- buffer used to count the number of elts to send
 *   send_shift     <-- index on ranks of the elements to send
 *   recv_count     <-- buffer used to count the number of elts to receive
 *   recv_shift     <-- index on ranks of the elements to receive
 *   vtx_merge_data <-- pointer to a cs_join_vertex_t structure
 *---------------------------------------------------------------------------*/

static void
_exchange_merged_vertices(const cs_join_mesh_t  *work,
                          const cs_gnum_t        vtx_tag[],
                          int                    send_count[],
                          int                    send_shift[],
                          int                    recv_count[],
                          int                    recv_shift[],
                          cs_join_vertex_t       vtx_merge_data[])
{
  cs_lnum_t  i;

  cs_join_vertex_t  *updated_vtx_data = NULL;
  int  n_ranks = cs_glob_n_ranks;
  cs_gnum_t  _n_ranks = n_ranks;
  MPI_Datatype  cs_mpi_join_vertex = cs_join_mesh_create_vtx_datatype();
  MPI_Comm  mpi_comm = cs_glob_mpi_comm;

  /* Allocate send_vtx_data and exchange vtx_merge_data */

  BFT_MALLOC(updated_vtx_data, send_shift[n_ranks], cs_join_vertex_t);

  MPI_Alltoallv(vtx_merge_data, recv_count, recv_shift, cs_mpi_join_vertex,
                updated_vtx_data, send_count, send_shift, cs_mpi_join_vertex,
                mpi_comm);

  /* Replace work->vertices by the updated structure after merge
     of vertices */

  assert(send_shift[n_ranks] == work->n_vertices);

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (i = 0; i < work->n_vertices; i++) {

    int  rank = (vtx_tag[i] - 1) % _n_ranks;
    cs_lnum_t  shift = send_shift[rank] + send_count[rank];

    work->vertices[i] = updated_vtx_data[shift];
    send_count[rank] += 1;

  }

  /* Free memory */

  MPI_Type_free(&cs_mpi_join_vertex);
  BFT_FREE(updated_vtx_data);
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Get the resulting cs_join_vertex_t structure after the merge of a set
 * of vertices.
 *
 * parameters:
 *   n_elts    <-- size of the set
 *   set       <-- set of vertices
 *
 * returns:
 *   a cs_join_vertex_t structure for the resulting vertex
 *---------------------------------------------------------------------------*/

static cs_join_vertex_t
_compute_merged_vertex(cs_lnum_t               n_elts,
                       const cs_join_vertex_t  set[])
{
  cs_lnum_t  i, k;
  cs_real_t  w;
  cs_join_vertex_t  mvtx;

  cs_real_t  denum = 0.0;

#if defined(DEBUG) && !defined(NDEBUG)
  /* Avoid Valgrind warnings in byte copies due to padding */
  memset(&mvtx, 0, sizeof(cs_join_vertex_t));
#endif

  assert(n_elts > 0);

  /* Initialize cs_join_vertex_t structure */

  mvtx.state = CS_JOIN_STATE_UNDEF;
  mvtx.gnum = set[0].gnum;
  mvtx.tolerance = set[0].tolerance;

  for (k = 0; k < 3; k++)
    mvtx.coord[k] = 0.0;

  /* Compute the resulting vertex */

  for (i = 0; i < n_elts; i++) {

    mvtx.tolerance = CS_MIN(set[i].tolerance, mvtx.tolerance);
    mvtx.gnum = CS_MIN(set[i].gnum, mvtx.gnum);
    mvtx.state = CS_MAX(set[i].state, mvtx.state);

    /* Compute the resulting coordinates of the merged vertices */

#if CS_JOIN_MERGE_INV_TOL
    w = 1.0/set[i].tolerance;
#else
    w = 1.0;
#endif
    denum += w;

    for (k = 0; k < 3; k++)
      mvtx.coord[k] += w * set[i].coord[k];

  }

  for (k = 0; k < 3; k++)
    mvtx.coord[k] /= denum;

  if (mvtx.state == CS_JOIN_STATE_ORIGIN)
    mvtx.state = CS_JOIN_STATE_MERGE;
  else if (mvtx.state == CS_JOIN_STATE_PERIO)
    mvtx.state = CS_JOIN_STATE_PERIO_MERGE;

  return mvtx;
}

/*----------------------------------------------------------------------------
 * Merge between identical vertices.
 *
 * Only the vertex numbering and the related tolerance may be different.
 * Store new data associated to the merged vertices in vertices array.
 *
 * parameters:
 *   param      <-- set of user-defined parameters
 *   merge_set  <-> a pointer to a cs_join_vertex_t structure which
 *                  stores data about merged vertices
 *   n_vertices <-- number of vertices in vertices array
 *   vertices   <-> array of cs_join_vertex_t structures
 *   equiv_gnum --> equivalence between id in vertices (same global number
 *                  initially or identical vertices: same coordinates)
 *---------------------------------------------------------------------------*/

static void
_pre_merge(cs_join_param_t     param,
           cs_join_gset_t     *merge_set,
           cs_join_vertex_t    vertices[],
           cs_join_gset_t    **p_equiv_gnum)
{
  cs_lnum_t  i, j, j1, j2, k, k1, k2, n_sub_elts;
  cs_real_t  deltad, deltat, limit, min_tol;
  cs_join_vertex_t  mvtx, coupled_vertices[2];

  cs_lnum_t  max_n_sub_elts = 0, n_local_pre_merge = 0;
  cs_lnum_t  *merge_index = merge_set->index;
  cs_gnum_t  *merge_list = merge_set->g_list;
  cs_gnum_t  *sub_list = NULL, *init_list = NULL;
  cs_join_gset_t  *equiv_gnum = NULL;

  const cs_real_t  pmf = param.pre_merge_factor;

  cs_join_gset_sort_sublist(merge_set);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  {
    int  len;
    FILE  *dbg_file = NULL;
    char  *filename = NULL;

    len = strlen("JoinDBG_InitMergeSet.dat")+1+2+4;
    BFT_MALLOC(filename, len, char);
    sprintf(filename, "Join%02dDBG_InitMergeSet%04d.dat",
            param.num, CS_MAX(cs_glob_rank_id, 0));
    dbg_file = fopen(filename, "w");

    cs_join_gset_dump(dbg_file, merge_set);

    fflush(dbg_file);
    BFT_FREE(filename);
    fclose(dbg_file);
  }
#endif

  /* Compute the max. size of a sub list */

  for (i = 0; i < merge_set->n_elts; i++)
    max_n_sub_elts = CS_MAX(max_n_sub_elts,
                            merge_index[i+1] - merge_index[i]);

  BFT_MALLOC(sub_list, max_n_sub_elts, cs_gnum_t);

  /* Store initial merge list */

  BFT_MALLOC(init_list, merge_index[merge_set->n_elts], cs_gnum_t);

  for (i = 0; i < merge_index[merge_set->n_elts]; i++)
    init_list[i] = merge_list[i];

  /* Apply merge */

  for (i = 0; i < merge_set->n_elts; i++) {

    cs_lnum_t  f_s = merge_index[i];
    cs_lnum_t  f_e = merge_index[i+1];

    n_sub_elts = f_e - f_s;

    for (j = f_s, k = 0; j < f_e; j++, k++)
      sub_list[k] = merge_list[j];

    for (j1 = 0; j1 < n_sub_elts - 1; j1++) {

      cs_lnum_t  v1_id = sub_list[j1];
      cs_join_vertex_t  v1 = vertices[v1_id];

      for (j2 = j1 + 1; j2 < n_sub_elts; j2++) {

        cs_lnum_t  v2_id = sub_list[j2];
        cs_join_vertex_t  v2 = vertices[v2_id];

        if (v1.gnum == v2.gnum) { /* Possible if n_ranks > 1 */

          if (sub_list[j1] < sub_list[j2])
            k1 = j1, k2 = j2;
          else
            k1 = j2, k2 = j1;

          for (k = 0; k < n_sub_elts; k++)
            if (sub_list[k] == sub_list[k2])
              sub_list[k] = sub_list[k1];

        }
        else {

          min_tol = CS_MIN(v1.tolerance, v2.tolerance);
          limit = min_tol * pmf;
          deltat = CS_ABS(v1.tolerance - v2.tolerance);

          if (deltat < limit) {

            deltad = _compute_length(v1, v2);

            if (deltad < limit) { /* Do a pre-merge */

              n_local_pre_merge++;

              if (v1.gnum < v2.gnum)
                k1 = j1, k2 = j2;
              else
                k1 = j2, k2 = j1;

              for (k = 0; k < n_sub_elts; k++)
                if (sub_list[k] == sub_list[k2])
                  sub_list[k] = sub_list[k1];

              coupled_vertices[0] = v1, coupled_vertices[1] = v2;
              mvtx = _compute_merged_vertex(2, coupled_vertices);
              vertices[v1_id] = mvtx;
              vertices[v2_id] = mvtx;

            } /* deltad < limit */

          } /* deltat < limit */

        } /* v1.gnum != v2.gnum */

      } /* End of loop on j2 */
    } /* End of loop on j1 */

    /* Update vertices */

    for (j = f_s, k = 0; j < f_e; j++, k++)
      vertices[merge_list[j]] = vertices[sub_list[k]];

    /* Update merge list */

    for (j = f_s, k = 0; j < f_e; j++, k++)
      merge_list[j] = sub_list[k];

  } /* End of loop on merge_set elements */

  /* Keep equivalences between identical vertices in equiv_gnum */

  equiv_gnum = cs_join_gset_create_by_equiv(merge_set, init_list);

  /* Clean merge set */

  cs_join_gset_clean(merge_set);

  /* Display information about the joining */

  if (param.verbosity > 0) {

    cs_gnum_t n_g_counter = n_local_pre_merge;
    cs_parall_counter(&n_g_counter, 1);

    bft_printf(_("\n  Pre-merge for %llu global element couples.\n"),
               (unsigned long long)n_g_counter);

    if (param.verbosity > 2) {
      fprintf(cs_glob_join_log, "\n  Local number of pre-merges: %d\n",
              n_local_pre_merge);
    }
  }

  /* Free memory */

  BFT_FREE(sub_list);
  BFT_FREE(init_list);

  /* Return pointer */

  *p_equiv_gnum = equiv_gnum;
}

/*----------------------------------------------------------------------------
 * Check if all vertices in the set include the ref_vertex in their tolerance.
 *
 * parameters:
 *   set_size   <-- size of set of vertices
 *   vertices   <-- set of vertices to check
 *   ref_vertex <-- ref. vertex
 *
 * returns:
 *   true if all vertices have ref_vertex in their tolerance, false otherwise
 *---------------------------------------------------------------------------*/

static bool
_is_in_tolerance(cs_lnum_t               set_size,
                 const cs_join_vertex_t  set[],
                 cs_join_vertex_t        ref_vertex)
{
  cs_lnum_t  i;

  for (i = 0; i < set_size; i++) {

    cs_real_t  d2ref = _compute_length(set[i], ref_vertex);
    cs_real_t  tolerance =  set[i].tolerance * cs_join_tol_eps_coef2;

    if (d2ref > tolerance)
      return false;

  }

  return true;
}

/*----------------------------------------------------------------------------
 * Test if we have to continue to the subset building
 *
 * parameters:
 *   set_size  <-- size of set
 *   prev_num  <-> array used to store previous subset_num
 *   new_num   <-> number associated to each vertices of the set
 *
 * returns:
 *   true or false
 *---------------------------------------------------------------------------*/

static bool
_continue_subset_building(int              set_size,
                          const cs_lnum_t  prev_num[],
                          const cs_lnum_t  new_num[])
{
  int  i;

  for (i = 0; i < set_size; i++)
    if (new_num[i] != prev_num[i])
      return true;

  return false;
}

/*----------------------------------------------------------------------------
 * Define subsets of vertices.
 *
 * parameters:
 *   set_size    <-- size of set
 *   state       <-- array keeping the state of the link
 *   subset_num  <-> number associated to each vertices of the set
 *---------------------------------------------------------------------------*/

static void
_iter_subset_building(cs_lnum_t               set_size,
                      const cs_lnum_t         state[],
                      cs_lnum_t               subset_num[])
{
  cs_lnum_t  i1, i2, k;

  for (k = 0, i1 = 0; i1 < set_size-1; i1++) {
    for (i2 = i1 + 1; i2 < set_size; i2++, k++) {

      if (state[k] == 1) { /* v1 - v2 are in tolerance each other */

        int _min = CS_MIN(subset_num[i1], subset_num[i2]);

        subset_num[i1] = _min;
        subset_num[i2] = _min;

      }

    }
  }

}

/*----------------------------------------------------------------------------
 * Define subsets of vertices.
 *
 * parameters:
 *   set_size    <-- size of set
 *   state       <-- array keeping the state of the link
 *   prev_num    <-> array used to store previous subset_num
 *   subset_num  <-> number associated to each vertices of the set
 *---------------------------------------------------------------------------*/

static void
_build_subsets(cs_lnum_t         set_size,
               const cs_lnum_t   state[],
               cs_lnum_t         prev_num[],
               cs_lnum_t         subset_num[])
{
  int  i;
  cs_lnum_t  n_loops = 0;

  /* Initialize */

  for (i = 0; i < set_size; i++) {
    subset_num[i] = i+1;
    prev_num[i] = subset_num[i];
  }

  _iter_subset_building(set_size, state, subset_num);

  while (   _continue_subset_building(set_size, prev_num, subset_num)
         && n_loops < CS_JOIN_MERGE_MAX_LOC_ITERS ) {

    n_loops++;

    for (i = 0; i < set_size; i++)
      prev_num[i] = subset_num[i];

    _iter_subset_building(set_size, state, subset_num);

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (cs_glob_join_log != NULL && n_loops >= CS_JOIN_MERGE_MAX_LOC_ITERS)
    fprintf(cs_glob_join_log,
            "WARNING max sub_loops to build subset reached\n");
#endif

}

/*----------------------------------------------------------------------------
 * Check if each subset is consistent with tolerance of vertices
 * If a transitivity is found, modify the state of the link
 * state = 1 (each other in their tolerance)
 *       = 0 (not each other in their tolerance)
 *
 * parameters:
 *   set_size    <-- size of set
 *   set         <-- pointer to a set of vertices
 *   state       <-> array keeping the state of the link
 *   subset_num  <-> number associated to each vertices of the set
 *   issues      <-> numbering of inconsistent subsets
 *   verbosity   <-- level of verbosity
 *
 * returns:
 *  number of subsets not consistent
 *---------------------------------------------------------------------------*/

static cs_lnum_t
_check_tol_consistency(cs_lnum_t               set_size,
                       const cs_join_vertex_t  set[],
                       cs_lnum_t               state[],
                       cs_lnum_t               subset_num[],
                       cs_lnum_t               issues[],
                       cs_lnum_t               verbosity)
{
  cs_lnum_t  i1, i2, j, k;

  cs_lnum_t  n_issues = 0;
  FILE  *logfile = cs_glob_join_log;

  for (k = 0, i1 = 0; i1 < set_size-1; i1++) {
    for (i2 = i1 + 1; i2 < set_size; i2++, k++) {

      if (state[k] == 0) {

        if (subset_num[i1] == subset_num[i2]) {

          if (verbosity > 4)
            fprintf(logfile,
                    " Transitivity detected between (%llu, %llu)\n",
                    (unsigned long long)set[i1].gnum,
                    (unsigned long long)set[i2].gnum);

          for (j = 0; j < n_issues; j++)
            if (issues[j] == subset_num[i1])
              break;
          if (j == n_issues)
            issues[n_issues++] = subset_num[i1];

        }
      }

    } /* End of loop on i2 */
  } /* ENd of loop on i1 */

  return  n_issues; /* Not a subset number */
}

/*----------------------------------------------------------------------------
 * Check if the merged vertex related to a subset is consistent with tolerance
 * of each vertex of the subset.
 *
 * parameters:
 *   set_size    <-- size of set
 *   subset_num  <-- number associated to each vertices of the set
 *   set         <-- pointer to a set of vertices
 *   merge_set   <-> merged vertex related to each subset
 *   work_set    <-> work array of vertices
 *
 * returns:
 *  true if all subsets are consistent otherwise false
 *---------------------------------------------------------------------------*/

static bool
_check_subset_consistency(cs_lnum_t               set_size,
                          const cs_lnum_t         subset_num[],
                          const cs_join_vertex_t  set[],
                          cs_join_vertex_t        merge_set[],
                          cs_join_vertex_t        work_set[])
{
  cs_lnum_t  i, set_id, subset_size;

  bool  is_consistent = true;

  /* Apply merged to each subset */

  for (set_id = 0; set_id < set_size; set_id++) {

    subset_size = 0;
    for (i = 0; i < set_size; i++)
      if (subset_num[i] == set_id+1)
        work_set[subset_size++] = set[i];

    if (subset_size > 0) {

      merge_set[set_id] = _compute_merged_vertex(subset_size, work_set);

      if (!_is_in_tolerance(subset_size, work_set, merge_set[set_id]))
        is_consistent = false;

    }

  } /* End of loop on subsets */

  return is_consistent;
}

/*----------------------------------------------------------------------------
 * Get position of the link between vertices i1 and i2.
 *
 * parameters:
 *   i1     <-- id in set for vertex 1
 *   i2     <-- id in set for vertex 2
 *   idx    <-- array of positions
 *
 * returns:
 *   position in an array like distances or state
 *---------------------------------------------------------------------------*/

inline static cs_lnum_t
_get_pos(cs_lnum_t       i1,
         cs_lnum_t       i2,
         const cs_lnum_t  idx[])
{
  cs_lnum_t  pos = -1;

  if (i1 < i2)
    pos = idx[i1] + i2-i1-1;
  else {
    assert(i1 != i2);
    pos = idx[i2] + i1-i2-1;
  }

  return pos;
}

/*----------------------------------------------------------------------------
 * Break equivalences for vertices implied in transitivity issue
 *
 * parameters:
 *   param       <-- parameters used to manage the joining algorithm
 *   set_size    <-- size of set
 *   set         <-- pointer to a set of vertices
 *   state       <-> array keeping the state of the link
 *   n_issues    <-- number of detected transitivity issues
 *   issues      <-- subset numbers of subset with a transitivity issue
 *   idx         <-- position of vertices couple in array like distances
 *   subset_num  <-- array of subset numbers
 *   distances   <-- array keeping the distances between vertices
 *---------------------------------------------------------------------------*/

static void
_break_equivalence(cs_join_param_t         param,
                   cs_lnum_t               set_size,
                   const cs_join_vertex_t  set[],
                   int                     state[],
                   cs_lnum_t               n_issues,
                   const int               issues[],
                   const int               idx[],
                   const int               subset_num[],
                   const double            distances[])
{
  cs_lnum_t  i, i1, i2, k;

  for (i = 0; i < n_issues; i++) {

    /* Find the weakest equivalence and break it.
       Purpose: Have the minimal number of equivalences to break
       for each subset where an inconsistency has been detected */

    int  i_save = 0;
    double rtf = -1.0, dist_save = 0.0;

    for (k = 0, i1 = 0; i1 < set_size-1; i1++) {
      for (i2 = i1 + 1; i2 < set_size; i2++, k++) {

        if (state[k] == 1) { /* v1, v2 are equivalent */

          if (   subset_num[i1] == issues[i]
              && subset_num[i2] == issues[i]) {

            /* Vertices belongs to a subset where an inconsistency
               has been found */

            double  rtf12 = distances[k]/set[i1].tolerance;
            double  rtf21 = distances[k]/set[i2].tolerance;

            assert(rtf12 < 1.0); /* Because they are equivalent */
            assert(rtf21 < 1.0);

            if (rtf12 >= rtf21) {
              if (rtf12 > rtf) {
                rtf = rtf12;
                i_save = i1;
                dist_save = distances[k];
              }
            }
            else {
              if (rtf21 > rtf) {
                rtf = rtf21;
                i_save = i2;
                dist_save = distances[k];
              }
            }

          }
        }

      } /* End of loop on i1 */
    } /* End of loop on i2 */

    if (rtf > 0.0) {

      /* Break equivalence between i_save and all vertices linked to
         i_save with a distance to i_save >= dist_save */

      for (i2 = 0; i2 < set_size; i2++) {

        if (i2 != i_save) {

          k = _get_pos(i_save, i2, idx);
          if (distances[k] >= dist_save && state[k] == 1) {

            state[k] = 0; /* Break equivalence */

            if (param.verbosity > 3)
              fprintf(cs_glob_join_log,
                      " %2d - Break equivalence between [%llu, %llu]"
                      " (dist_ref: %6.4e)\n",
                      issues[i],
                      (unsigned long long)set[i_save].gnum,
                      (unsigned long long)set[i2].gnum, dist_save);

          }
        }

      } /* End of loop on vertices */

    } /* rtf > 0.0 */

  } /* End of loop on issues */

}

/*----------------------------------------------------------------------------
 * Break equivalences between vertices until each vertex of the list has
 * the resulting vertex of the merge under its tolerance.
 *
 * parameters:
 *   param         <-- set of user-defined parameters
 *   set_size      <-- size of the set of vertices
 *   set           <-> set of vertices
 *   vbuf          <-> tmp buffer
 *   rbuf          <-> tmp buffer
 *   ibuf          <-> tmp buffer
 *
 * returns:
 *   number of loops necessary to build consistent subsets
 *---------------------------------------------------------------------------*/

static cs_lnum_t
_solve_transitivity(cs_join_param_t    param,
                    cs_lnum_t          set_size,
                    cs_join_vertex_t   set[],
                    cs_join_vertex_t   vbuf[],
                    cs_real_t          rbuf[],
                    cs_lnum_t          ibuf[])
{
  cs_lnum_t  i1, i2, k, n_issues;

  cs_lnum_t  n_loops = 0;
  bool  is_end = false;
  cs_lnum_t  *subset_num = NULL, *state = NULL, *prev_num = NULL;
  cs_lnum_t  *subset_issues = NULL, *idx = NULL;
  cs_real_t  *distances = NULL;
  cs_join_vertex_t  *merge_set = NULL, *work_set = NULL;

  /* Sanity checks */

  assert(set_size > 0);

  /* Define temporary buffers */

  subset_num = &(ibuf[0]);
  prev_num = &(ibuf[set_size]);
  subset_issues = &(ibuf[2*set_size]);
  idx = &(ibuf[3*set_size]);
  state = &(ibuf[4*set_size]);
  distances = &(rbuf[0]);
  merge_set = &(vbuf[0]);
  work_set = &(vbuf[set_size]);

  /* Compute distances between each couple of vertices among the set */

  for (k = 0, i1 = 0; i1 < set_size-1; i1++)
    for (i2 = i1 + 1; i2 < set_size; i2++, k++)
      distances[k] = _compute_length(set[i1], set[i2]);

  /* Compute initial state of equivalences between vertices */

  for (k = 0, i1 = 0; i1 < set_size-1; i1++) {
    for (i2 = i1 + 1; i2 < set_size; i2++, k++) {
      if (   set[i1].tolerance < distances[k]
          || set[i2].tolerance < distances[k])
        state[k] = 0;
      else
        state[k] = 1;
    }
  }

  idx[0] = 0;
  for (k = 1; k < set_size - 1; k++)
    idx[k] = set_size - k + idx[k-1];

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (cs_glob_join_log != NULL) {
    cs_join_dump_array(cs_glob_join_log, "double", "\nDist",
                       set_size*(set_size-1)/2, distances);
    cs_join_dump_array(cs_glob_join_log, "int", "\nInit. State",
                       set_size*(set_size-1)/2, state);
  }
#endif

  _build_subsets(set_size, state, prev_num, subset_num);

  while (is_end == false && n_loops < param.n_max_equiv_breaks) {

    n_loops++;

    n_issues = _check_tol_consistency(set_size,
                                      set,
                                      state,
                                      subset_num,
                                      subset_issues,
                                      param.verbosity);

    if (n_issues > 0)
      _break_equivalence(param,
                         set_size,
                         set,
                         state,
                         n_issues,
                         subset_issues,
                         idx,
                         subset_num,
                         distances);

    _build_subsets(set_size, state, prev_num, subset_num);

    is_end = _check_subset_consistency(set_size,
                                       subset_num,
                                       set,
                                       merge_set,
                                       work_set);

  } /* End of while */

  if (param.verbosity > 3) {

    fprintf(cs_glob_join_log,
            " Number of tolerance reductions:  %4d\n", n_loops);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    cs_join_dump_array(cs_glob_join_log, "int", "\nFinal Subset",
                       set_size, subset_num);
#endif

  }

  /* Apply merged to each subset */

  for (k = 0; k < set_size; k++)
    set[k] = merge_set[subset_num[k]-1];

  return n_loops;
}

/*----------------------------------------------------------------------------
 * Merge between vertices. Store new data associated to the merged vertices
 * in vertices.
 *
 * parameters:
 *   param      <-- set of user-defined parameters
 *   merge_set  <-> a pointer to a cs_join_vertex_t structure which
 *                  stores data about merged vertices
 *   n_vertices <-- number of vertices in vertices array
 *   vertices   <-> array of cs_join_vertex_t structures
 *---------------------------------------------------------------------------*/

static void
_merge_vertices(cs_join_param_t    param,
                cs_join_gset_t    *merge_set,
                cs_lnum_t          n_vertices,
                cs_join_vertex_t   vertices[])
{
  cs_lnum_t  i, j, k, list_size;
  cs_join_vertex_t  merged_vertex;
  bool  ok;

  cs_lnum_t  max_list_size = 0, vv_max_list_size = 0;
  cs_lnum_t  n_loops = 0, n_max_loops = 0, n_transitivity = 0;

  cs_join_gset_t  *equiv_gnum = NULL;
  cs_real_t  *rbuf = NULL;
  cs_lnum_t  *merge_index = NULL, *ibuf = NULL;
  cs_gnum_t  *merge_list = NULL, *merge_ref_elts = NULL;
  cs_gnum_t  *list = NULL;
  cs_join_vertex_t  *set = NULL, *vbuf = NULL;
  FILE  *logfile = cs_glob_join_log;

  const int  verbosity = param.verbosity;

  /* Sanity check */

  assert(param.merge_tol_coef >= 0.0);

  /* Pre-merge of identical vertices */

  _pre_merge(param, merge_set, vertices, &equiv_gnum);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  {
    int  len;
    FILE  *dbg_file = NULL;
    char  *filename = NULL;

    len = strlen("JoinDBG_MergeSet.dat")+1+2+4;
    BFT_MALLOC(filename, len, char);
    sprintf(filename, "Join%02dDBG_MergeSet%04d.dat",
            param.num, CS_MAX(cs_glob_rank_id, 0));
    dbg_file = fopen(filename, "w");

    cs_join_gset_dump(dbg_file, merge_set);

    fflush(dbg_file);
    BFT_FREE(filename);
    fclose(dbg_file);
  }
#endif /* defined(DEBUG) && !defined(NDEBUG) */

  /* Modify the tolerance for the merge operation if needed */

  if (fabs(param.merge_tol_coef - 1.0) > 1e-30) {
    for (i = 0; i < n_vertices; i++)
      vertices[i].tolerance *= param.merge_tol_coef;
  }

  /* Compute the max. size of a sub list */

  merge_index = merge_set->index;
  merge_list = merge_set->g_list;
  merge_ref_elts = merge_set->g_elts;

  for (i = 0; i < merge_set->n_elts; i++) {
    list_size = merge_index[i+1] - merge_index[i];
    max_list_size = CS_MAX(max_list_size, list_size);
  }
  vv_max_list_size = ((max_list_size-1)*max_list_size)/2;

  if (verbosity > 0) {   /* Display information */

    cs_lnum_t g_max_list_size = max_list_size;
    cs_parall_counter_max(&g_max_list_size, 1);

    if (g_max_list_size < 2) {
      cs_join_gset_destroy(&equiv_gnum);
      bft_printf(_("\n  No need to merge vertices.\n"));
      return;
    }
    else
      bft_printf(_("\n  Max size of a merge set of vertices: %llu\n"),
                 (unsigned long long)g_max_list_size);
  }

  /* Temporary buffers allocation */

  BFT_MALLOC(ibuf, 4*max_list_size + vv_max_list_size, cs_lnum_t);
  BFT_MALLOC(rbuf, vv_max_list_size, cs_real_t);
  BFT_MALLOC(vbuf, 2*max_list_size, cs_join_vertex_t);
  BFT_MALLOC(list, max_list_size, cs_gnum_t);
  BFT_MALLOC(set, max_list_size, cs_join_vertex_t);

  /* Merge set of vertices */

  for (i = 0; i < merge_set->n_elts; i++) {

    list_size = merge_index[i+1] - merge_index[i];

    if (list_size > 1) {

      for (j = 0, k = merge_index[i]; k < merge_index[i+1]; k++, j++) {
        list[j] = merge_list[k];
        set[j] = vertices[list[j]];
      }

      /* Define the resulting cs_join_vertex_t structure of the merge */

      merged_vertex = _compute_merged_vertex(list_size, set);

      /* Check if the vertex resulting of the merge is in the tolerance
         for each vertex of the list */

      ok = _is_in_tolerance(list_size, set, merged_vertex);

#if CS_JOIN_MERGE_TOL_REDUC
      if (ok == false) { /*
                            The merged vertex is not in the tolerance of
                            each vertex. This is a transitivity problem.
                            We have to split the initial set into several
                            subsets.
                         */

        n_transitivity++;

        /* Display information on vertices to merge */
        if (verbosity > 3) {
          fprintf(logfile,
                  "\n Begin merge for ref. elt: %llu - list_size: %d\n",
                  (unsigned long long)merge_ref_elts[i],
                  merge_index[i+1] - merge_index[i]);
          for (j = 0; j < list_size; j++) {
            fprintf(logfile, "%9llu -", (unsigned long long)list[j]);
            cs_join_mesh_dump_vertex(logfile, set[j]);
          }
          fprintf(logfile, "\nMerged vertex rejected:\n");
          cs_join_mesh_dump_vertex(logfile, merged_vertex);
        }

        n_loops = _solve_transitivity(param,
                                      list_size,
                                      set,
                                      vbuf,
                                      rbuf,
                                      ibuf);

        for (j = 0; j < list_size; j++)
          vertices[list[j]] = set[j];

        n_max_loops = CS_MAX(n_max_loops, n_loops);

        if (verbosity > 3) { /* Display information */
          fprintf(logfile, "\n  %3d loop(s) to get consistent subsets\n",
                  n_loops);
          fprintf(logfile, "\n End merge for ref. elt: %llu - list_size: %d\n",
                  (unsigned long long)merge_ref_elts[i],
                  merge_index[i+1] - merge_index[i]);
          for (j = 0; j < list_size; j++) {
            fprintf(logfile, "%7llu -", (unsigned long long)list[j]);
            cs_join_mesh_dump_vertex(logfile, vertices[list[j]]);
          }
          fprintf(logfile, "\n");
        }

      }
      else /* New vertex data for the sub-elements */

#endif /* CS_JOIN_MERGE_TOL_REDUC */

        for (j = 0; j < list_size; j++)
          vertices[list[j]] = merged_vertex;

    } /* list_size > 1 */

  } /* End of loop on potential merges */

  /* Apply merge to vertex initially identical */

  if (equiv_gnum != NULL) {

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    {
      int  len;
      FILE  *dbg_file = NULL;
      char  *filename = NULL;

      len = strlen("JoinDBG_EquivMerge.dat")+1+2+4;
      BFT_MALLOC(filename, len, char);
      sprintf(filename, "Join%02dDBG_EquivMerge%04d.dat",
              param.num, CS_MAX(cs_glob_rank_id, 0));
      dbg_file = fopen(filename, "w");

      cs_join_gset_dump(dbg_file, equiv_gnum);

      fflush(dbg_file);
      BFT_FREE(filename);
      fclose(dbg_file);
    }
#endif /* defined(DEBUG) && !defined(NDEBUG) */

    for (i = 0; i < equiv_gnum->n_elts; i++) {

      cs_lnum_t  start = equiv_gnum->index[i];
      cs_lnum_t  end = equiv_gnum->index[i+1];
      cs_lnum_t  ref_id = equiv_gnum->g_elts[i];

      for (j = start; j < end; j++)
        vertices[equiv_gnum->g_list[j]] = vertices[ref_id];

    }
  }

  if (verbosity > 0) {

    cs_gnum_t n_g_counter = n_transitivity;
    cs_parall_counter(&n_g_counter, 1);

    bft_printf(_("\n  Excessive transitivity for %llu set(s) of vertices.\n"),
               (unsigned long long)n_g_counter);

    if (verbosity > 1) {
      cs_lnum_t g_n_max_loops = n_max_loops;
      cs_parall_counter_max(&g_n_max_loops, 1);
      bft_printf(_("\n  Max. number of iterations to solve transitivity"
                   " excess: %llu\n"), (unsigned long long)g_n_max_loops);
    }
  }

  /* Free memory */

  BFT_FREE(ibuf);
  BFT_FREE(vbuf);
  BFT_FREE(rbuf);
  BFT_FREE(set);
  BFT_FREE(list);

  cs_join_gset_destroy(&equiv_gnum);
}

/*----------------------------------------------------------------------------
 * Keep an history of the evolution of each vertex id before/after the merge
 * operation.
 *
 * parameters:
 *   n_iwm_vertices    <-- number of vertices before intersection for the
 *                          work cs_join_mesh_t structure
 *   iwm_vtx_gnum      <-- initial global vertex num. (work mesh struct.)
 *   init_max_vtx_gnum <-- initial max. global numbering for vertices
 *   n_vertices        <-- number of vertices before merge/after intersection
 *   vertices          <-- array of cs_join_vertex_t structures
 *   p_o2n_vtx_gnum    --> distributed array by block on the new global vertex
 *                         numbering for the initial vertices (before inter.)
 *---------------------------------------------------------------------------*/

static void
_keep_global_vtx_evolution(cs_lnum_t               n_iwm_vertices,
                           const cs_gnum_t         iwm_vtx_gnum[],
                           cs_gnum_t               init_max_vtx_gnum,
                           cs_lnum_t               n_vertices,
                           const cs_join_vertex_t  vertices[],
                           cs_gnum_t              *p_o2n_vtx_gnum[])
{
  cs_lnum_t  i;

  cs_gnum_t  *o2n_vtx_gnum = NULL;

  const int  n_ranks = cs_glob_n_ranks;
  const int  local_rank = CS_MAX(cs_glob_rank_id, 0);

  assert(n_iwm_vertices <= n_vertices); /* after inter. >= init */

  if (n_ranks == 1) {

    BFT_MALLOC(o2n_vtx_gnum, n_iwm_vertices, cs_gnum_t);

    for (i = 0; i < n_iwm_vertices; i++)
      o2n_vtx_gnum[i] = vertices[i].gnum;

  }

#if defined(HAVE_MPI) /* Parallel treatment */
  if (n_ranks > 1) {

    cs_lnum_t  shift, rank;
    cs_lnum_t  block_size = 0;

    int        *send_shift = NULL, *recv_shift = NULL;
    int        *send_count = NULL, *recv_count = NULL;
    cs_gnum_t  *send_glist = NULL, *recv_glist = NULL;

    cs_block_dist_info_t  bi = cs_block_dist_compute_sizes(local_rank,
                                                           n_ranks,
                                                           1,
                                                           0,
                                                           init_max_vtx_gnum);

    MPI_Comm  mpi_comm = cs_glob_mpi_comm;

    if (bi.gnum_range[1] > bi.gnum_range[0])
      block_size = bi.gnum_range[1] - bi.gnum_range[0];

    /* Initialize o2n_vtx_gnum */

    BFT_MALLOC(o2n_vtx_gnum, block_size, cs_gnum_t);

    for (i = 0; i < block_size; i++) {
      cs_gnum_t g_id = i;
      o2n_vtx_gnum[i] = bi.gnum_range[0] + g_id;
    }

    /* Send new vtx global number to the related rank = the good block */

    BFT_MALLOC(send_count, n_ranks, int);
    BFT_MALLOC(recv_count, n_ranks, int);

    for (i = 0; i < n_ranks; i++)
      send_count[i] = 0;

    for (i = 0; i < n_iwm_vertices; i++) {
      rank = (iwm_vtx_gnum[i] - 1)/(cs_gnum_t)(bi.block_size);
      send_count[rank] += 2;
    }

    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);

    BFT_MALLOC(send_shift, n_ranks + 1, cs_lnum_t);
    BFT_MALLOC(recv_shift, n_ranks + 1, cs_lnum_t);

    send_shift[0] = 0;
    recv_shift[0] = 0;

    for (rank = 0; rank < n_ranks; rank++) {
      send_shift[rank + 1] = send_shift[rank] + send_count[rank];
      recv_shift[rank + 1] = recv_shift[rank] + recv_count[rank];
    }

    assert(send_shift[n_ranks] == 2*n_iwm_vertices);

    /* Build send_list */

    BFT_MALLOC(send_glist, send_shift[n_ranks], cs_gnum_t);
    BFT_MALLOC(recv_glist, recv_shift[n_ranks], cs_gnum_t);

    for (i = 0; i < n_ranks; i++)
      send_count[i] = 0;

    for (i = 0; i < n_iwm_vertices; i++) {

      rank = (iwm_vtx_gnum[i] - 1)/(cs_gnum_t)(bi.block_size);
      shift = send_shift[rank] + send_count[rank];

      send_glist[shift] = iwm_vtx_gnum[i];  /* Old global number */
      send_glist[shift+1] = vertices[i].gnum;   /* New global number */
      send_count[rank] += 2;

    }

    MPI_Alltoallv(send_glist, send_count, send_shift, CS_MPI_GNUM,
                  recv_glist, recv_count, recv_shift, CS_MPI_GNUM,
                  mpi_comm);

    BFT_FREE(send_count);
    BFT_FREE(send_shift);
    BFT_FREE(send_glist);
    BFT_FREE(recv_count);

    /* Update o2n_vtx_gnum */

    for (rank = 0; rank < n_ranks; rank++) {

      for (i = recv_shift[rank]; i < recv_shift[rank+1]; i+=2) {

        cs_gnum_t  o_gnum = recv_glist[i];
        cs_gnum_t  n_gnum = recv_glist[i+1];
        cs_lnum_t  id = o_gnum - bi.gnum_range[0];

#if 0 && defined(DEBUG) && !defined(NDEBUG)
        if (o2n_vtx_gnum[id] != bi.gnum_range[0] + id)
          assert(o2n_vtx_gnum[id] == n_gnum);
#endif

        o2n_vtx_gnum[id] = n_gnum;

      }

    } /* End of loop on ranks */

    BFT_FREE(recv_shift);
    BFT_FREE(recv_glist);

  }
#endif /* HAVE_MPI */

  /* Set return pointer */

  *p_o2n_vtx_gnum = o2n_vtx_gnum;
}

/*----------------------------------------------------------------------------
 * Keep a history of the evolution of each vertex id before/after the merge
 * operation for the current mesh (local point of view).
 *
 * parameters:
 *   n_vertices      <-- number of vertices before merge/after intersection
 *   vertices        <-- array of cs_join_vertex_t structures
 *   p_n_am_vertices --> number of vertices after the merge step
 *   p_o2n_vtx_id    --> array keeping the evolution of the vertex ids
 *---------------------------------------------------------------------------*/

static void
_keep_local_vtx_evolution(cs_lnum_t                n_vertices,
                          const cs_join_vertex_t   vertices[],
                          cs_lnum_t               *p_n_am_vertices,
                          cs_lnum_t               *p_o2n_vtx_id[])
{
  cs_lnum_t  i;
  cs_gnum_t  prev;

  cs_lnum_t  n_am_vertices = 0;
  cs_lnum_t  *o2n_vtx_id = NULL;
  cs_lnum_t  *order = NULL;
  cs_gnum_t  *vtx_gnum = NULL;

  if (n_vertices == 0)
    return;

  BFT_MALLOC(vtx_gnum, n_vertices, cs_gnum_t);

  for (i = 0; i < n_vertices; i++)
    vtx_gnum[i] = vertices[i].gnum;

  /* Order vertices according to their global numbering */

  BFT_MALLOC(order, n_vertices, cs_lnum_t);

  cs_order_gnum_allocated(NULL, vtx_gnum, order, n_vertices);

  /* Delete vertices sharing the same global number. Keep only one */

  BFT_MALLOC(o2n_vtx_id, n_vertices, cs_lnum_t);

  prev = vtx_gnum[order[0]];
  o2n_vtx_id[order[0]] = n_am_vertices;

  for (i = 1; i < n_vertices; i++) {

    cs_lnum_t  o_id = order[i];
    cs_gnum_t  cur = vtx_gnum[o_id];

    if (cur != prev) {
      prev = cur;
      n_am_vertices++;
      o2n_vtx_id[o_id] = n_am_vertices;
    }
    else
      o2n_vtx_id[o_id] = n_am_vertices;

  } /* End of loop on vertices */

  /* n_am_vertices is an id */
  n_am_vertices += 1;

  assert(n_am_vertices <= n_vertices); /* after merge <= after inter. */

  /* Free memory */

  BFT_FREE(order);
  BFT_FREE(vtx_gnum);

  /* Set return pointers */

  *p_n_am_vertices = n_am_vertices;
  *p_o2n_vtx_id = o2n_vtx_id;
}

/*----------------------------------------------------------------------------
 * Search for new elements to add to the definition of the current edge
 * These new sub-elements come from initial edges which are now (after the
 * merge step) sub-edge of the current edge
 * Count step
 *
 * parameters:
 *   edge_id        <-- id of the edge to deal with
 *   inter_edges    <-- structure keeping data on edge intersections
 *   edges          <-- edges definition
 *   n_iwm_vertices <-- initial number of vertices in work_mesh
 *   n_new_sub_elts --> number of new elements to add in the edge definition
 *
 * returns:
 *  number of new sub-elements related to this edge
 *---------------------------------------------------------------------------*/

static cs_lnum_t
_count_new_sub_edge_elts(cs_lnum_t                     edge_id,
                         const cs_join_inter_edges_t  *inter_edges,
                         const cs_join_edges_t        *edges,
                         cs_lnum_t                     n_iwm_vertices)
{
  cs_lnum_t  j, k, j1, j2, sub_edge_id;
  cs_lnum_t  start, end, _start, _end, v1_num, v2_num;
  bool  found;

  cs_lnum_t  n_new_sub_elts = 0;

  start = inter_edges->index[edge_id];
  end = inter_edges->index[edge_id+1];

  for (j1 = start; j1 < end-1; j1++) {

    v1_num = inter_edges->vtx_lst[j1];

    if (v1_num <= n_iwm_vertices) { /* v1 is an initial vertex */
      for (j2 = j1+1; j2 < end; j2++) {

        v2_num = inter_edges->vtx_lst[j2];

        if (v2_num <= n_iwm_vertices) { /* (v1,v2) is an initial edge */

          sub_edge_id = CS_ABS(cs_join_mesh_get_edge(v1_num,
                                                     v2_num,
                                                     edges))-1;
          assert(sub_edge_id != -1);
          _start = inter_edges->index[sub_edge_id];
          _end = inter_edges->index[sub_edge_id+1];

          for (j = _start; j < _end; j++) {

            found = false;
            for (k = j1+1; k < j2; k++)
              if (inter_edges->vtx_glst[k] == inter_edges->vtx_glst[j])
                found = true;

            if (found == false)
              n_new_sub_elts += 1;

          } /* End of loop on sub_edge_id definition */

        }

      } /* End of loop on j2 */
    }

  } /* End of loop on j1 */

  return n_new_sub_elts;
}

/*----------------------------------------------------------------------------
 * Update a cs_join_inter_edges_t structure after the merge operation.
 * cs_join_inter_edges_t structure should be not NULL.
 *
 * parameters:
 *   param          <-- user-defined parameters for the joining algorithm
 *   n_iwm_vertices <-- initial number of vertices in work_mesh
 *   o2n_vtx_id     <-- array keeping the evolution of the vertex ids
 *   edges          <-- edges definition
 *   p_inter_edges  <-> pointer to the structure keeping data on
 *                      edge intersections
 *---------------------------------------------------------------------------*/

static void
_update_inter_edges_after_merge(cs_join_param_t          param,
                                cs_lnum_t                n_iwm_vertices,
                                const cs_lnum_t          o2n_vtx_id[],
                                const cs_join_edges_t   *edges,
                                const cs_join_mesh_t    *mesh,
                                cs_join_inter_edges_t  **p_inter_edges)
{
  cs_lnum_t  i, j,k, j1, j2,  start_shift, idx_shift;
  cs_lnum_t  save, _start, _end, start, end;
  cs_lnum_t  v1_num, v2_num, v1_id, v2_id, sub_edge_id;
  cs_gnum_t  v1_gnum, v2_gnum, new_gnum, prev_gnum;
  bool  found;

  cs_lnum_t  n_adds = 0;

  cs_join_inter_edges_t  *inter_edges = *p_inter_edges;
  cs_join_inter_edges_t  *new_inter_edges = NULL;
  cs_lnum_t  n_edges = inter_edges->n_edges;
  cs_lnum_t  init_list_size = inter_edges->index[n_edges];
  FILE  *logfile = cs_glob_join_log;

  assert(n_edges == edges->n_edges);

  /* Define vtx_glst to compare global vertex numbering */

  if (inter_edges->vtx_glst == NULL)
    BFT_MALLOC(inter_edges->vtx_glst, inter_edges->index[n_edges], cs_gnum_t);

  for (i = 0; i < inter_edges->index[n_edges]; i++) {
    v1_id = inter_edges->vtx_lst[i] - 1;
    inter_edges->vtx_glst[i] = mesh->vertices[v1_id].gnum;
  }

  /* Delete redundancies of vertices sharing the same global numbering
     after the merge step and define a new index */

  idx_shift = 0;
  save = inter_edges->index[0];

  for (i = 0; i < n_edges; i++) {

    start = save;
    end = inter_edges->index[i+1];

    if (end - start > 0) {

      start_shift = start;
      v1_id = edges->def[2*i] - 1;
      v2_id = edges->def[2*i+1] - 1;
      v1_gnum = mesh->vertices[v1_id].gnum;
      v2_gnum = mesh->vertices[v2_id].gnum;
      prev_gnum = inter_edges->vtx_glst[start_shift];

      /* Don't take into account vertices with the same number as the
         first edge element */

      while (prev_gnum == v1_gnum && start_shift + 1 < end)
        prev_gnum = inter_edges->vtx_glst[++start_shift];

      if (prev_gnum != v1_gnum && start_shift < end) {

        inter_edges->vtx_lst[idx_shift] = inter_edges->vtx_lst[start_shift];
        inter_edges->abs_lst[idx_shift] = inter_edges->abs_lst[start_shift];
        inter_edges->vtx_glst[idx_shift] = inter_edges->vtx_glst[start_shift];
        idx_shift += 1;

        for (j = start_shift + 1; j < end; j++) {

          new_gnum = inter_edges->vtx_glst[j];

          /* Don't take into account redundancies and vertices with the same
             number as the second edge element */

          if (prev_gnum != new_gnum && new_gnum != v2_gnum) {
            prev_gnum = new_gnum;
            inter_edges->vtx_lst[idx_shift] = inter_edges->vtx_lst[j];
            inter_edges->abs_lst[idx_shift] = inter_edges->abs_lst[j];
            inter_edges->vtx_glst[idx_shift] = inter_edges->vtx_glst[j];
            idx_shift += 1;
          }

        }

      } /* If start_shift < end */

    } /* end - start > 0 */

    save = inter_edges->index[i+1];
    inter_edges->index[i+1] = idx_shift;

  } /* End of loop on edge intersections */

  inter_edges->max_sub_size = 0;

  for (i = 0; i < n_edges; i++)
    inter_edges->max_sub_size =
            CS_MAX(inter_edges->max_sub_size,
                       inter_edges->index[i+1] - inter_edges->index[i]);

  assert(inter_edges->index[n_edges] <= init_list_size);

  BFT_REALLOC(inter_edges->vtx_lst, inter_edges->index[n_edges], cs_lnum_t);
  BFT_REALLOC(inter_edges->abs_lst, inter_edges->index[n_edges], cs_coord_t);

#if 0 && defined(DEBUG) && !defined(NDEBUG) /* Dump local structures */
  fprintf(logfile, " AFTER REDUNDANCIES CLEAN\n");
  cs_join_inter_edges_dump(logfile, inter_edges, edges, mesh);
#endif

  /* Add new vertices from initial edges which are now sub-edges */

  for (i = 0; i < n_edges; i++)
    n_adds += _count_new_sub_edge_elts(i, inter_edges, edges, n_iwm_vertices);

  if (param.verbosity > 2)
    fprintf(logfile,
            "  Number of sub-elements to add to edge definition: %8d\n",
            n_adds);

  if (n_adds > 0) { /* Define a new inter_edges structure */

    new_inter_edges = cs_join_inter_edges_create(n_edges);

    BFT_MALLOC(new_inter_edges->vtx_lst,
               inter_edges->index[n_edges] + n_adds, cs_lnum_t);
    BFT_MALLOC(new_inter_edges->abs_lst,
               inter_edges->index[n_edges] + n_adds, cs_coord_t);

    n_adds = 0;
    idx_shift = 0;
    new_inter_edges->index[0] = 0;

    for (i = 0; i < n_edges; i++) {

      new_inter_edges->edge_gnum[i] = inter_edges->edge_gnum[i];
      start = inter_edges->index[i];
      end = inter_edges->index[i+1];

      if (start - end > 0) {

        for (j1 = start; j1 < end-1; j1++) {

          v1_num = inter_edges->vtx_lst[j1];
          new_inter_edges->vtx_lst[idx_shift] = v1_num;
          new_inter_edges->abs_lst[idx_shift] = inter_edges->abs_lst[j1];
          idx_shift++;

          if (v1_num <= n_iwm_vertices) { /* v1 is an initial vertex */
            for (j2 = j1+1; j2 < end; j2++) {

              v2_num = inter_edges->vtx_lst[j2];

              if (v2_num <= n_iwm_vertices) { /* (v1,v2) is an initial edge */

                sub_edge_id =
                  CS_ABS(cs_join_mesh_get_edge(v1_num, v2_num, edges))-1;
                assert(sub_edge_id != -1);

                _start = inter_edges->index[sub_edge_id];
                _end = inter_edges->index[sub_edge_id+1];

                for (j = _start; j < _end; j++) {

                  found = false;
                  for (k = j1+1; k < j2; k++)
                    if (inter_edges->vtx_glst[k] == inter_edges->vtx_glst[j])
                      found = true;

                  if (found == false) {

                    new_inter_edges->vtx_lst[idx_shift] =
                      inter_edges->vtx_lst[j];
                    new_inter_edges->abs_lst[idx_shift] =
                      inter_edges->abs_lst[j];
                    idx_shift++;

                  }

                } /* End of loop on sub_edge_id definition */

              }

            } /* End of loop on j2 */
          }

        } /* End of loop on j1 */

        /* Add last vertex in the previous edge definition */

        new_inter_edges->vtx_lst[idx_shift] = inter_edges->vtx_lst[end-1];
        new_inter_edges->abs_lst[idx_shift] = inter_edges->abs_lst[end-1];
        idx_shift++;

      } /* If end - start > 0 */

      new_inter_edges->index[i+1] = idx_shift;

    } /* End of loop on edges */

    cs_join_inter_edges_destroy(&inter_edges);
    inter_edges = new_inter_edges;

    inter_edges->max_sub_size = 0;

    for (i = 0; i < n_edges; i++)
      inter_edges->max_sub_size = CS_MAX(inter_edges->max_sub_size,
                                         inter_edges->index[i+1]);

  } /* End if n_adds > 0 */

#if 0 && defined(DEBUG) && !defined(NDEBUG) /* Dump local structures */
  if (logfile != NULL) {
    fprintf(logfile, " AFTER SUB ELTS ADD\n");
    cs_join_inter_edges_dump(logfile, inter_edges, edges, mesh);
  }
#endif

  /* Update cs_join_inter_edges_t structure */

  for (i = 0; i < n_edges; i++) {

    start = inter_edges->index[i];
    end = inter_edges->index[i+1];

    for (j = start; j < end; j++) {

      cs_lnum_t old_id = inter_edges->vtx_lst[j] - 1;

      inter_edges->vtx_lst[j] = o2n_vtx_id[old_id] + 1;

    }

  }

  /* Return pointer */

  *p_inter_edges = inter_edges;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Define send_rank_index and send_faces to prepare the exchange of new faces
 * between mesh structures.
 *
 * parameters:
 *   n_send          <-- number of faces to send
 *   n_g_faces       <-- global number of faces to be joined
 *   face_gnum       <-- global face number
 *   gnum_rank_index <-- index on ranks for the init. global face numbering
 *   p_n_send        --> number of face/rank couples to send
 *   p_send_rank     --> rank ids to which to send
 *   p_send_faces    --> list of face ids to send
 *---------------------------------------------------------------------------*/

static void
_get_faces_to_send(cs_lnum_t         n_faces,
                   cs_gnum_t         n_g_faces,
                   const cs_gnum_t   face_gnum[],
                   const cs_gnum_t   gnum_rank_index[],
                   cs_lnum_t        *p_n_send,
                   int              *p_send_rank[],
                   cs_lnum_t        *p_send_faces[])
{
  cs_lnum_t  i, rank, reduce_rank;
  cs_block_dist_info_t  bi;

  cs_lnum_t   reduce_size = 0, n_send = 0;
  int        *send_rank = NULL;
  cs_lnum_t  *send_faces = NULL;
  cs_lnum_t  *reduce_ids = NULL, *count = NULL;
  cs_gnum_t  *reduce_index = NULL;

  const int  local_rank = CS_MAX(cs_glob_rank_id, 0);
  const int  n_ranks = cs_glob_n_ranks;

  /* Sanity checks */

  assert(gnum_rank_index != NULL);
  assert(n_ranks > 1);

  /* Compute block_size */

  bi = cs_block_dist_compute_sizes(local_rank,
                                   n_ranks,
                                   1,
                                   0,
                                   n_g_faces);

  /* Compact init. global face distribution. Remove ranks without face
     at the begining */

  for (i = 0; i < n_ranks; i++)
    if (gnum_rank_index[i] < gnum_rank_index[i+1])
      reduce_size++;

  BFT_MALLOC(reduce_index, reduce_size+1, cs_gnum_t);
  BFT_MALLOC(reduce_ids, reduce_size, cs_lnum_t);

  reduce_size = 0;
  reduce_index[0] = gnum_rank_index[0] + 1;

  for (i = 0; i < n_ranks; i++) {

    /* Add +1 to gnum_rank_index because it's an id and we work on numbers */

    if (gnum_rank_index[i] < gnum_rank_index[i+1]) {
      reduce_index[reduce_size+1] = gnum_rank_index[i+1] + 1;
      reduce_ids[reduce_size++] = i;
    }

  }

  BFT_MALLOC(send_rank, n_faces, int);
  BFT_MALLOC(send_faces, n_faces, cs_lnum_t);

  /* Fill the list of ranks */

  n_send = 0;

  for (i = 0; i < n_faces; i++) {

    if (face_gnum[i] >= bi.gnum_range[0] && face_gnum[i] < bi.gnum_range[1]) {

      /* The current face is a "main" face for the local rank */

      reduce_rank = cs_search_gindex_binary(reduce_size,
                                            face_gnum[i],
                                            reduce_index);

      assert(reduce_rank > -1);
      assert(reduce_rank < reduce_size);

      rank = reduce_ids[reduce_rank];
      send_rank[n_send] = rank;
      send_faces[n_send] = i;
      n_send += 1;

    } /* End of loop on initial faces */

  }

  BFT_REALLOC(send_rank, n_send, int);
  BFT_REALLOC(send_faces, n_send, cs_lnum_t);

  /* Free memory */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (cs_glob_join_log != NULL) {
    FILE *logfile = cs_glob_join_log;
    for (i = 0; i < n_send; i++)
      fprintf(logfile, " %d (%llu) to rank %d\n",
              send_faces[i], (unsigned long long)face_gnum[send_faces[i]],
              send_rank[i]);
    fflush(logfile);
  }
#endif

  BFT_FREE(count);
  BFT_FREE(reduce_ids);
  BFT_FREE(reduce_index);

  /* Set return pointers */

  *p_n_send = n_send;
  *p_send_rank = send_rank;
  *p_send_faces = send_faces;
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Update local_mesh by redistributing mesh.
 * Send back to the original rank the new face and vertex description.
 *
 * parameters:
 *   gnum_rank_index  <--  index on ranks for the old global face numbering
 *   send_mesh        <--  distributed mesh on faces to join
 *   p_recv_mesh      <->  mesh on local selected faces to be joined
 *---------------------------------------------------------------------------*/

static void
_redistribute_mesh(const cs_gnum_t         gnum_rank_index[],
                   const cs_join_mesh_t   *send_mesh,
                   cs_join_mesh_t        **p_recv_mesh)
{
  cs_join_mesh_t  *recv_mesh = *p_recv_mesh;

  const int  n_ranks = cs_glob_n_ranks;

  /* sanity checks */

  assert(send_mesh != NULL);
  assert(recv_mesh != NULL);

  if (n_ranks == 1)
    cs_join_mesh_copy(&recv_mesh, send_mesh);

#if defined(HAVE_MPI)
  if (n_ranks > 1) { /* Parallel mode */

    cs_lnum_t   n_send = 0;
    int        *send_rank = NULL;
    cs_lnum_t  *send_faces = NULL;

    MPI_Comm  mpi_comm = cs_glob_mpi_comm;

    /* Free some structures of the mesh */

    cs_join_mesh_reset(recv_mesh);

    _get_faces_to_send(send_mesh->n_faces,
                       send_mesh->n_g_faces,
                       send_mesh->face_gnum,
                       gnum_rank_index,
                       &n_send,
                       &send_rank,
                       &send_faces);

    assert(n_send <= send_mesh->n_faces);

    /* Get the new face connectivity from the distributed send_mesh */

    cs_join_mesh_exchange(n_send,
                          send_rank,
                          send_faces,
                          send_mesh,
                          recv_mesh,
                          mpi_comm);

    BFT_FREE(send_faces);
    BFT_FREE(send_rank);

  }
#endif

  /* Return pointers */

  *p_recv_mesh = recv_mesh;

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Creation of new vertices.
 *
 * Update list of equivalent vertices, and assign a vertex (existing or
 * newly created) to each intersection.
 *
 * parameters:
 *   verbosity          <-- verbosity level
 *   edges              <-- list of edges
 *   work               <-> joining mesh maintaining initial vertex data
 *   inter_set          <-> cs_join_inter_set_t structure including
 *                          data on edge-edge  intersections
 *   init_max_vtx_gnum  <-- initial max. global numbering for vertices
 *   p_n_g_new_vertices <-> pointer to the global number of new vertices
 *   p_vtx_eset         <-> pointer to a structure dealing with vertex
 *                          equivalences
 *---------------------------------------------------------------------------*/

void
cs_join_create_new_vertices(int                     verbosity,
                            const cs_join_edges_t  *edges,
                            cs_join_mesh_t         *work,
                            cs_join_inter_set_t    *inter_set,
                            cs_gnum_t               init_max_vtx_gnum,
                            cs_gnum_t              *p_n_g_new_vertices,
                            cs_join_eset_t        **p_vtx_eset)
{
  cs_lnum_t  i, shift;
  double  tol_min;
  cs_join_vertex_t  v1, v2;

  cs_lnum_t  n_new_vertices = 0;
  cs_gnum_t  n_g_new_vertices = 0;
  cs_gnum_t  *new_vtx_gnum = NULL;
  cs_lnum_t  n_iwm_vertices = work->n_vertices;
  cs_join_eset_t  *vtx_equiv = *p_vtx_eset;

  /* Count the number of new vertices. Update cs_join_inter_set_t struct. */

  for (i = 0; i < inter_set->n_inter; i++) {

    cs_join_inter_t  inter1 = inter_set->inter_lst[2*i];
    cs_join_inter_t  inter2 = inter_set->inter_lst[2*i+1];

    inter1.vtx_id = _get_vtx_id(inter1,
                                &(edges->def[2*inter1.edge_id]),
                                n_iwm_vertices,
                                &n_new_vertices);

    inter2.vtx_id = _get_vtx_id(inter2,
                                &(edges->def[2*inter2.edge_id]),
                                n_iwm_vertices,
                                &n_new_vertices);

    inter_set->inter_lst[2*i] = inter1;
    inter_set->inter_lst[2*i+1] = inter2;

  } /* End of loop on intersections */

  /* Compute the global numbering for the new vertices (Take into account
     potential redundancies) */

  _compute_new_vertex_gnum(work,
                           edges,
                           inter_set,
                           init_max_vtx_gnum,
                           n_iwm_vertices,
                           n_new_vertices,
                           &n_g_new_vertices,
                           &new_vtx_gnum);

  if (verbosity > 0)
    bft_printf(_("\n  Global number of new vertices to create: %10llu\n"),
               (unsigned long long)n_g_new_vertices);

  /* Define new vertices */

  work->n_vertices += n_new_vertices;
  work->n_g_vertices += n_g_new_vertices;

  BFT_REALLOC(work->vertices, work->n_vertices, cs_join_vertex_t);

#if defined(DEBUG) && !defined(NDEBUG) /* Prepare sanity checks */
  {
    cs_join_vertex_t  incoherency;

    /* Initialize to incoherent values new vertices structures */

    incoherency.gnum = 0;
    incoherency.coord[0] = -9999.9999;
    incoherency.coord[1] = -9999.9999;
    incoherency.coord[2] = -9999.9999;
    incoherency.tolerance = -1.0;
    incoherency.state = CS_JOIN_STATE_UNDEF;

    for (i = 0; i < n_new_vertices; i++)
      work->vertices[n_iwm_vertices + i] = incoherency;

  }
#endif

  /* Fill vertices structure with new vertex definitions */

  for (i = 0; i < inter_set->n_inter; i++) {

    cs_join_inter_t  inter1 = inter_set->inter_lst[2*i];
    cs_join_inter_t  inter2 = inter_set->inter_lst[2*i+1];
    cs_lnum_t  v1_num = inter1.vtx_id + 1;
    cs_lnum_t  v2_num = inter2.vtx_id + 1;
    cs_lnum_t  equiv_id = vtx_equiv->n_equiv;

    assert(inter1.vtx_id < work->n_vertices);
    assert(inter2.vtx_id < work->n_vertices);

    /* Create new vertices if needed */

    if (v1_num > n_iwm_vertices) {

      shift = inter1.vtx_id - n_iwm_vertices;
      v1 = _get_new_vertex(inter1.curv_abs,
                           new_vtx_gnum[shift],
                           &(edges->def[2*inter1.edge_id]),
                           work);
      tol_min = v1.tolerance;

    }
    else
      tol_min = work->vertices[v1_num-1].tolerance;

    if (v2_num > n_iwm_vertices) {

      shift = inter2.vtx_id - n_iwm_vertices;
      v2 = _get_new_vertex(inter2.curv_abs,
                           new_vtx_gnum[shift],
                           &(edges->def[2*inter2.edge_id]),
                           work);
      tol_min = CS_MIN(tol_min, v2.tolerance);

    }
    else
      tol_min = CS_MIN(tol_min, work->vertices[v2_num-1].tolerance);

    /* A new vertex has a tolerance equal to the minimal tolerance
       between the two vertices implied in the intersection */

    if (v1_num > n_iwm_vertices) {
      v1.tolerance = tol_min;
      work->vertices[inter1.vtx_id] = v1;
    }
    if (v2_num > n_iwm_vertices) {
      v2.tolerance = tol_min;
      work->vertices[inter2.vtx_id] = v2;
    }

    /* Add equivalence between the two current vertices */

    cs_join_eset_check_size(equiv_id, &vtx_equiv);

    if (v1_num < v2_num) {
      vtx_equiv->equiv_couple[2*equiv_id] = v1_num;
      vtx_equiv->equiv_couple[2*equiv_id+1] = v2_num;
    }
    else {
      vtx_equiv->equiv_couple[2*equiv_id] = v2_num;
      vtx_equiv->equiv_couple[2*equiv_id+1] = v1_num;
    }

    vtx_equiv->n_equiv += 1;

  } /* End of loop on intersections */

  /* Free memory */

  BFT_FREE(new_vtx_gnum);

#if defined(DEBUG) && !defined(NDEBUG) /* Sanity checks */
  for (i = 0; i < work->n_vertices; i++) {

    cs_join_vertex_t  vtx = work->vertices[i];

    if (vtx.gnum == 0 || vtx.tolerance < -0.99)
      bft_error(__FILE__, __LINE__, 0,
                _("  Inconsistent value found in cs_join_vertex_t struct.:\n"
                  "    Vertex %d is defined by:\n"
                  "      %llu - [%7.4le, %7.4le, %7.4le] - %lg\n"),
                i, (unsigned long long)vtx.gnum,
                vtx.coord[0], vtx.coord[1], vtx.coord[2],
                vtx.tolerance);

  } /* End of loop on vertices */

#if 0
  _dump_vtx_eset(vtx_equiv, work, cs_glob_join_log);
#endif

#endif

  /* Set return pointers */

  *p_n_g_new_vertices = n_g_new_vertices;
  *p_vtx_eset = vtx_equiv;
}

/*----------------------------------------------------------------------------
 * Merge of equivalent vertices (and tolerance reduction if necessary)
 *
 * Define a new cs_join_vertex_t structure (stored in "work" structure).
 * Returns an updated cs_join_mesh_t and cs_join_edges_t structures.
 *
 * parameters:
 *   param            <-- set of user-defined parameters for the joining
 *   n_g_vertices_tot <-- global number of vertices (initial parent mesh)
 *   work             <-> pointer to a cs_join_mesh_t structure
 *   vtx_eset         <-- structure storing equivalences between vertices
 *                        (two vertices are equivalent if they are within
 *                        each other's tolerance)
 *---------------------------------------------------------------------------*/

void
cs_join_merge_vertices(cs_join_param_t        param,
                       cs_gnum_t              n_g_vertices_tot,
                       cs_join_mesh_t        *work,
                       const cs_join_eset_t  *vtx_eset)
{
  cs_gnum_t  *vtx_tags = NULL;
  cs_join_gset_t  *merge_set = NULL;

  const int  n_ranks = cs_glob_n_ranks;

  /* Initialize counters for the merge operation */

  _initialize_merge_counter();

#if 0 && defined(DEBUG) && !defined(NDEBUG) /* Dump local structures */
  _dump_vtx_eset(vtx_eset, work, cs_glob_join_log);
#endif

  if (param.verbosity > 2) {
    cs_gnum_t g_n_equiv = vtx_eset->n_equiv;
    cs_parall_counter(&g_n_equiv, 1);
    fprintf(cs_glob_join_log,
            "\n"
            "  Final number of equiv. between vertices; local: %9d\n"
            "                                          global: %9llu\n",
            vtx_eset->n_equiv, (unsigned long long)g_n_equiv);
  }

  /* Operate merge between equivalent vertices.
     Manage reduction of tolerance if necessary */

  /* Tag with the same number all the vertices which might be merged together */

  _tag_equiv_vertices(n_g_vertices_tot,
                      vtx_eset,
                      work,
                      param.verbosity,
                      &vtx_tags);

  if (n_ranks == 1) { /* Serial mode */

    /* Build a merge list */

    merge_set = cs_join_gset_create_from_tag(work->n_vertices, vtx_tags);

    /* Merge of equivalent vertices */

    _merge_vertices(param,
                    merge_set,
                    work->n_vertices,
                    work->vertices);

  }

#if defined(HAVE_MPI)
  if (n_ranks > 1) { /* Parallel mode: we work by block */

    int  *send_count = NULL, *recv_count = NULL;
    int  *send_shift = NULL, *recv_shift = NULL;
    cs_join_vertex_t  *vtx_merge_data = NULL;

    BFT_MALLOC(send_count, n_ranks, int);
    BFT_MALLOC(recv_count, n_ranks, int);
    BFT_MALLOC(send_shift, n_ranks+1, int);
    BFT_MALLOC(recv_shift, n_ranks+1, int);

    /* Build a merge list in parallel */

    _build_parall_merge_structures(work,
                                   vtx_tags,
                                   send_count, send_shift,
                                   recv_count, recv_shift,
                                   &vtx_merge_data,
                                   &merge_set);

    /* Merge of equivalent vertices for the current block */

    _merge_vertices(param,
                    merge_set,
                    recv_shift[n_ranks],
                    vtx_merge_data);

    _exchange_merged_vertices(work,
                              vtx_tags,
                              send_count, send_shift,
                              recv_count, recv_shift,
                              vtx_merge_data);

    BFT_FREE(send_count);
    BFT_FREE(send_shift);
    BFT_FREE(recv_count);
    BFT_FREE(recv_shift);
    BFT_FREE(vtx_merge_data);

  }
#endif /* HAVE_MPI */

  /* Free memory */

  BFT_FREE(vtx_tags);

  cs_join_gset_destroy(&merge_set);

  if (param.verbosity > 1)
    bft_printf(_("\n"
                 "  Merging of equivalent vertices done.\n"));
}

/*----------------------------------------------------------------------------
 * Merge of equivalent vertices (and reduction of tolerance if necessary)
 *
 * Define a new cs_join_vertex_t structure (stored in "work" structure)
 * Returns an updated cs_join_mesh_t and cs_join_edges_t structures.
 *
 * parameters:
 *   param                <-- set of user-defined parameters for the joining
 *   n_iwm_vertices       <-- initial number of vertices (work mesh struct.)
 *   iwm_vtx_gnum         <-- initial global vertex num. (work mesh struct)
 *   init_max_vtx_gnum    <-- initial max. global numbering for vertices
 *   rank_face_gnum_index <-- index on face global numbering to determine
 *                            the related rank
 *   p_mesh               <-> pointer to cs_join_mesh_t structure
 *   p_edges              <-> pointer to cs_join_edges_t structure
 *   p_inter_edges        <-> pointer to a cs_join_inter_edges_t struct.
 *   p_local_mesh         <-> pointer to a cs_join_mesh_t structure
 *   p_o2n_vtx_gnum       --> array on blocks on the new global vertex
 *                            numbering for the init. vertices (before inter.)
 *---------------------------------------------------------------------------*/

void
cs_join_merge_update_struct(cs_join_param_t          param,
                            cs_lnum_t                n_iwm_vertices,
                            const cs_gnum_t          iwm_vtx_gnum[],
                            cs_gnum_t                init_max_vtx_gnum,
                            const cs_gnum_t          rank_face_gnum_index[],
                            cs_join_mesh_t         **p_mesh,
                            cs_join_edges_t        **p_edges,
                            cs_join_inter_edges_t  **p_inter_edges,
                            cs_join_mesh_t         **p_local_mesh,
                            cs_gnum_t               *p_o2n_vtx_gnum[])
{
  cs_lnum_t  n_am_vertices = 0; /* new number of vertices after merge */
  cs_lnum_t  *o2n_vtx_id = NULL;
  cs_gnum_t  *o2n_vtx_gnum = NULL;
  cs_join_mesh_t  *mesh = *p_mesh;
  cs_join_mesh_t  *local_mesh = *p_local_mesh;
  cs_join_edges_t  *edges = *p_edges;
  cs_join_inter_edges_t  *inter_edges = *p_inter_edges;

  /* Keep an history of the evolution of each vertex */

  _keep_global_vtx_evolution(n_iwm_vertices,   /* n_vertices before inter */
                             iwm_vtx_gnum,
                             init_max_vtx_gnum,
                             mesh->n_vertices, /* n_vertices after inter */
                             mesh->vertices,
                             &o2n_vtx_gnum);   /* defined by block in // */

  _keep_local_vtx_evolution(mesh->n_vertices, /* n_vertices after inter */
                            mesh->vertices,
                            &n_am_vertices,   /* n_vertices after merge */
                            &o2n_vtx_id);

  /* Update all structures which keeps data about vertices */

  if (inter_edges != NULL) { /* The join type is not conform */

    /* Update inter_edges structure */

    _update_inter_edges_after_merge(param,
                                    n_iwm_vertices,
                                    o2n_vtx_id,  /* size of mesh->n_vertices */
                                    edges,
                                    mesh,
                                    &inter_edges);

    assert(edges->n_edges == inter_edges->n_edges);  /* Else: problem for
                                                        future synchro. */

    /* Update cs_join_mesh_t structure after the merge of vertices
       numbering of the old vertices + add new vertices */

    cs_join_mesh_update(mesh,
                        edges,
                        inter_edges->index,
                        inter_edges->vtx_lst,
                        n_am_vertices,
                        o2n_vtx_id);

  } /* End if inter_edges != NULL */

  else
    /* Update cs_join_mesh_t structure after the merge of vertices
       numbering of the old vertices + add new vertices */

    cs_join_mesh_update(mesh,
                        edges,
                        NULL,
                        NULL,
                        n_am_vertices,
                        o2n_vtx_id);

  BFT_FREE(o2n_vtx_id);

  /* Update local_mesh by redistributing mesh */

  _redistribute_mesh(rank_face_gnum_index,
                     mesh,
                     &local_mesh);

  /* Set return pointers */

  *p_mesh = mesh;
  *p_edges = edges;
  *p_inter_edges = inter_edges;
  *p_o2n_vtx_gnum = o2n_vtx_gnum;
  *p_local_mesh = local_mesh;
}

/*---------------------------------------------------------------------------*/

END_C_DECLS
