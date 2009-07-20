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
 * Set of subroutines for:
 *  - merging equivalent vertices,
 *  - managing tolerance reduction
 *===========================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

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
 * BFT library headers
 *---------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_timer.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM headers
 *---------------------------------------------------------------------------*/

#include <fvm_io_num.h>
#include <fvm_order.h>
#include <fvm_parall.h>

/*----------------------------------------------------------------------------
 * Local headers
 *---------------------------------------------------------------------------*/

#include "cs_search.h"
#include "cs_join_post.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *---------------------------------------------------------------------------*/

#include "cs_join_merge.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Structure and type definitions
 *===========================================================================*/

/*============================================================================
 * Macro definitions
 *===========================================================================*/

/* Turn on (1) or off (0) the tolerance reduc. */
#define  CS_JOIN_MERGE_TOL_REDUC  1
#define  CS_JOIN_MERGE_INV_TOL  0

/*============================================================================
 * Global variable definitions
 *===========================================================================*/

/* Parameters to control the vertex merge */

enum {

  CS_JOIN_MERGE_MAX_GLOB_ITERS = 5,  /* Max. number of glob. iter. for finding
                                        equivalent vertices */
  CS_JOIN_MERGE_MAX_LOC_ITERS = 15,  /* Max. number of loc. iter. for finding
                                        equivalent vertices */
  CS_JOIN_MERGE_MAX_REDUCTIONS = 100 /* Max. number of loc. iter. for
                                        reducing the tolerance */

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
 *   e_set <-- cs_join_eset_t structure to dump
 *   mesh  <-- cs_join_mesh_t structure associated
 *---------------------------------------------------------------------------*/

static void
_dump_vtx_eset(const cs_join_eset_t    *e_set,
               const cs_join_mesh_t    *mesh)
{
  int  i;

  bft_printf("\n  Dump an cs_join_eset_t structure (%p)\n", e_set);
  bft_printf("  n_max_equiv: %10d\n", e_set->n_max_equiv);
  bft_printf("  n_equiv    : %10d\n\n", e_set->n_equiv);

  for (i = 0; i < e_set->n_equiv; i++) {

    cs_int_t  v1_num = e_set->equiv_couple[2*i];
    cs_int_t  v2_num = e_set->equiv_couple[2*i+1];

    bft_printf(" %10d - local: (%9d, %9d) - global: (%10u, %10u)\n",
               i, v1_num, v2_num, (mesh->vertices[v1_num-1]).gnum,
               (mesh->vertices[v2_num-1]).gnum);

  }
  bft_printf_flush();
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
  cs_int_t  k;
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
_get_new_vertex(float                  curv_abs,
                fvm_gnum_t             gnum,
                const cs_int_t         vtx_couple[],
                const cs_join_mesh_t  *work)
{
  cs_int_t  k;
  cs_join_vertex_t  new_vtx_data;

  cs_join_vertex_t  v1 = work->vertices[vtx_couple[0]-1];
  cs_join_vertex_t  v2 = work->vertices[vtx_couple[1]-1];

  assert(curv_abs >= 0.0);
  assert(curv_abs <= 1.0);

  /* New global number */

  new_vtx_data.gnum = gnum;

  /* New tolerance */

  new_vtx_data.tolerance = (1-curv_abs)*v1.tolerance + curv_abs*v2.tolerance;

  /* New coordinates */

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
_define_inter_tag(fvm_gnum_t  tag[],
                  fvm_gnum_t  e1_gnum,
                  fvm_gnum_t  e2_gnum,
                  fvm_gnum_t  link_vtx_gnum)
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
 *   n_g_vertices       <-- global number of vertices (initial full mesh)
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
                         fvm_gnum_t                  n_g_vertices,
                         cs_int_t                    n_iwm_vertices,
                         cs_int_t                    n_new_vertices,
                         fvm_gnum_t                 *p_n_g_new_vertices,
                         fvm_gnum_t                 *p_new_vtx_gnum[])
{
  cs_int_t  i;

  fvm_gnum_t  n_g_new_vertices = 0;
  cs_int_t  n_new_vertices_save = n_new_vertices;
  fvm_lnum_t  *order = NULL;
  fvm_gnum_t  *inter_tag = NULL, *adjacency = NULL, *new_vtx_gnum = NULL;
  fvm_io_num_t  *new_vtx_io_num = NULL;

  /* Define a fvm_io_num_t structure to get the global numbering
     for the new vertices.
     First, build a tag associated to each intersection */

  BFT_MALLOC(new_vtx_gnum, n_new_vertices, fvm_gnum_t);
  BFT_MALLOC(inter_tag, 3*n_new_vertices, fvm_gnum_t);

  n_new_vertices = 0;

  for (i = 0; i < inter_set->n_inter; i++) {

    cs_join_inter_t  inter1 = inter_set->inter_lst[2*i];
    cs_join_inter_t  inter2 = inter_set->inter_lst[2*i+1];
    fvm_gnum_t  e1_gnum = edges->gnum[inter1.edge_id];
    fvm_gnum_t  e2_gnum = edges->gnum[inter2.edge_id];

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
                          n_g_vertices + 1);
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

  BFT_MALLOC(order, n_new_vertices, fvm_lnum_t);

  fvm_order_local_allocated_s(NULL, inter_tag, 3, order, n_new_vertices);

  BFT_MALLOC(adjacency, 3*n_new_vertices, fvm_gnum_t);

  for (i = 0; i < n_new_vertices; i++) {

    cs_int_t  o_id = order[i];

    adjacency[3*i] = inter_tag[3*o_id];
    adjacency[3*i+1] = inter_tag[3*o_id+1];
    adjacency[3*i+2] = inter_tag[3*o_id+2];

  }

  BFT_FREE(inter_tag);

  if (cs_glob_n_ranks > 1) {

    const fvm_gnum_t  *global_num = NULL;

    new_vtx_io_num =
      fvm_io_num_create_from_adj_s(NULL, adjacency, n_new_vertices, 3);

    n_g_new_vertices = fvm_io_num_get_global_count(new_vtx_io_num);
    global_num = fvm_io_num_get_global_num(new_vtx_io_num);

    for (i = 0; i < n_new_vertices; i++)
      new_vtx_gnum[order[i]] = global_num[i] + n_g_vertices;

    fvm_io_num_destroy(new_vtx_io_num);

  } /* End of parallel treatment */

  else {

    if (n_new_vertices > 0) {

      fvm_gnum_t  new_gnum = n_g_vertices + 1;

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

static cs_int_t
_get_vtx_id(cs_join_inter_t  inter,
            const cs_int_t   vtx_couple[],
            cs_int_t         n_init_vertices,
            cs_int_t        *p_n_new_vertices)
{
  cs_int_t  vtx_id = -1;
  cs_int_t  n_new_vertices = *p_n_new_vertices;

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

static cs_bool_t
_is_spread_not_converged(cs_int_t          n_vertices,
                         const fvm_gnum_t  prev_vtx_tag[],
                         const fvm_gnum_t  vtx_tag[])
{
  cs_int_t  i;

  cs_bool_t  have_to_continue = true;

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
 *  vtx_eset <-- structure dealing with vertices equivalences
 *  vtx_tag  <-> tag for each vertex
 *---------------------------------------------------------------------------*/

static void
_spread_tag(const cs_join_eset_t  *vtx_eset,
            fvm_gnum_t             vtx_tag[])
{
  cs_int_t  i;

  cs_int_t  *equiv_lst = vtx_eset->equiv_couple;

  for (i = 0; i < vtx_eset->n_equiv; i++) {

    cs_int_t  v1_id = equiv_lst[2*i] - 1;
    cs_int_t  v2_id = equiv_lst[2*i+1] - 1;
    fvm_gnum_t  v1_gnum = vtx_tag[v1_id];
    fvm_gnum_t  v2_gnum = vtx_tag[v2_id];

    if (v1_gnum != v2_gnum) {

      fvm_gnum_t  min_gnum = CS_MIN(v1_gnum, v2_gnum);

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
              cs_int_t               n_vertices,
              fvm_gnum_t             prev_vtx_tag[],
              fvm_gnum_t             vtx_tag[])
{
  cs_int_t  i;

  _loc_merge_counter++;

  _spread_tag(vtx_eset, vtx_tag);

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

    _spread_tag(vtx_eset, vtx_tag);
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

static cs_bool_t
_global_spread(cs_int_t               block_size,
               const cs_join_mesh_t  *work,
               fvm_gnum_t             vtx_tag[],
               fvm_gnum_t             glob_vtx_tag[],
               fvm_gnum_t             prev_glob_vtx_tag[],
               fvm_gnum_t             recv2glob[],
               cs_int_t               send_count[],
               cs_int_t               send_shift[],
               fvm_gnum_t             send_glob_buffer[],
               cs_int_t               recv_count[],
               cs_int_t               recv_shift[],
               fvm_gnum_t             recv_glob_buffer[])
{
  cs_bool_t  ret_value;
  cs_int_t  i, local_value, global_value;

  cs_int_t  n_vertices = work->n_vertices;
  int  n_ranks = cs_glob_n_ranks;
  MPI_Comm  mpi_comm = cs_glob_mpi_comm;

  _glob_merge_counter++;

  /* Push modifications in local vtx_tag to the global vtx_tag */

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (i = 0; i < n_vertices; i++) {

    int  rank = (work->vertices[i].gnum - 1) % n_ranks;
    cs_int_t  shift = send_shift[rank] + send_count[rank];

    send_glob_buffer[shift] = vtx_tag[i];
    send_count[rank] += 1;

  }

  MPI_Alltoallv(send_glob_buffer, send_count, send_shift, FVM_MPI_GNUM,
                recv_glob_buffer, recv_count, recv_shift, FVM_MPI_GNUM,
                mpi_comm);

  /* Apply update to glob_vtx_tag */

  for (i = 0; i < recv_shift[n_ranks]; i++) {
    cs_int_t  cur_id = recv2glob[i];
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

    MPI_Alltoallv(recv_glob_buffer, recv_count, recv_shift, FVM_MPI_GNUM,
                  send_glob_buffer, send_count, send_shift, FVM_MPI_GNUM,
                  mpi_comm);

    /* Update vtx_tag */

    for (i = 0; i < n_ranks; i++)
      send_count[i] = 0;

    assert(send_shift[n_ranks] == n_vertices);

    for (i = 0; i < n_vertices; i++) {

      int  rank = (work->vertices[i].gnum - 1) % n_ranks;
      cs_int_t  shift = send_shift[rank] + send_count[rank];

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
_parall_tag_init(fvm_gnum_t             n_g_vertices_to_treat,
                 const cs_join_mesh_t  *work,
                 cs_int_t              *p_block_size,
                 cs_int_t              *p_send_count[],
                 cs_int_t              *p_send_shift[],
                 fvm_gnum_t            *p_send_glob_buf[],
                 cs_int_t              *p_recv_count[],
                 cs_int_t              *p_recv_shift[],
                 fvm_gnum_t            *p_recv_glob_buf[],
                 fvm_gnum_t            *p_recv2glob[],
                 fvm_gnum_t            *p_glob_vtx_tag[],
                 fvm_gnum_t            *p_prev_glob_vtx_tag[])
{
  cs_int_t  i;

  cs_int_t  n_vertices = work->n_vertices;
  cs_int_t  block_size = 0, left_over = 0;
  cs_int_t  *send_count = NULL, *recv_count = NULL;
  cs_int_t  *send_shift = NULL, *recv_shift = NULL;
  fvm_gnum_t  *recv2glob = NULL;
  fvm_gnum_t  *recv_glob_buffer = NULL, *send_glob_buffer = NULL;
  fvm_gnum_t  *glob_vtx_tag = NULL, *prev_glob_vtx_tag = NULL;
  MPI_Comm  mpi_comm = cs_glob_mpi_comm;

  const int  n_ranks = cs_glob_n_ranks;
  const int  local_rank = CS_MAX(cs_glob_rank_id, 0);

  /* Allocate and intialize vtx_tag associated to the local rank */

  block_size = n_g_vertices_to_treat / n_ranks;
  left_over = n_g_vertices_to_treat % n_ranks;
  if (local_rank < left_over)
    block_size += 1;

  BFT_MALLOC(glob_vtx_tag, block_size, fvm_gnum_t);
  BFT_MALLOC(prev_glob_vtx_tag, block_size, fvm_gnum_t);

  for (i = 0; i < block_size; i++) {
    prev_glob_vtx_tag[i] = i*n_ranks + local_rank + 1;
    glob_vtx_tag[i] = i*n_ranks + local_rank + 1;
  }

  /* Allocate and define send/recv_count/shift */

  BFT_MALLOC(send_count, n_ranks, cs_int_t);
  BFT_MALLOC(recv_count, n_ranks, cs_int_t);
  BFT_MALLOC(send_shift, n_ranks + 1, cs_int_t);
  BFT_MALLOC(recv_shift, n_ranks + 1, cs_int_t);

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (i = 0; i < n_vertices; i++) {
    int  rank = (work->vertices[i].gnum - 1) % n_ranks;
    send_count[rank] += 1;
  }

  MPI_Alltoall(&(send_count[0]), 1, FVM_MPI_LNUM,
               &(recv_count[0]), 1, FVM_MPI_LNUM, mpi_comm);

  /* Build index */

  for (i = 0; i < n_ranks; i++) {
    send_shift[i+1] = send_shift[i] + send_count[i];
    recv_shift[i+1] = recv_shift[i] + recv_count[i];
  }

  assert(send_shift[n_ranks] == n_vertices);

  /* Allocate and define recv2glob */

  BFT_MALLOC(send_glob_buffer, send_shift[n_ranks], fvm_gnum_t);
  BFT_MALLOC(recv2glob, recv_shift[n_ranks], fvm_gnum_t);
  BFT_MALLOC(recv_glob_buffer, recv_shift[n_ranks], fvm_gnum_t);

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (i = 0; i < n_vertices; i++) {

    int  rank = (work->vertices[i].gnum - 1) % n_ranks;
    cs_int_t  shift = send_shift[rank] + send_count[rank];

    send_glob_buffer[shift] = (work->vertices[i].gnum - 1) / n_ranks;
    send_count[rank] += 1;

  }

  MPI_Alltoallv(send_glob_buffer, send_count, send_shift, FVM_MPI_GNUM,
                recv2glob, recv_count, recv_shift, FVM_MPI_GNUM,
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
 *   verbosity        <-- level of accuracy in information display
 *   p_vtx_tag        --> pointer to the vtx_tag for the local vertices
 *---------------------------------------------------------------------------*/

static void
_tag_equiv_vertices(fvm_gnum_t             n_g_vertices_tot,
                    const cs_join_eset_t  *vtx_eset,
                    const cs_join_mesh_t  *work,
                    int                    verbosity,
                    fvm_gnum_t            *p_vtx_tag[])
{
  cs_int_t  i;

  fvm_gnum_t  *vtx_tag = NULL;
  fvm_gnum_t  *prev_vtx_tag = NULL;

  const cs_int_t  n_vertices = work->n_vertices;
  const int  n_ranks = cs_glob_n_ranks;

  /* Local initialization : we tag each vertex by its global number */

  BFT_MALLOC(prev_vtx_tag, n_vertices, fvm_gnum_t);
  BFT_MALLOC(vtx_tag, n_vertices, fvm_gnum_t);

  for (i = 0; i < work->n_vertices; i++) {

    fvm_gnum_t  v_gnum = work->vertices[i].gnum;

    vtx_tag[i] = v_gnum;
    prev_vtx_tag[i] = v_gnum;

  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  for (i = 0; i < n_vertices; i++)
    bft_printf(" Initial vtx_tag[%6d] = %9u\n", i, vtx_tag[i]);
  bft_printf_flush();
#endif

  /* Compute vtx_tag */

  _local_spread(vtx_eset, n_vertices, prev_vtx_tag, vtx_tag);

  if (n_ranks > 1) { /* Parallel treatment */

    cs_bool_t  go_on;

    cs_int_t  block_size = 0;
    cs_int_t  *send_count = NULL, *recv_count = NULL;
    cs_int_t  *send_shift = NULL, *recv_shift = NULL;
    fvm_gnum_t  *recv2glob = NULL;
    fvm_gnum_t  *recv_glob_buffer = NULL, *send_glob_buffer = NULL;
    fvm_gnum_t  *glob_vtx_tag = NULL, *prev_glob_vtx_tag = NULL;

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

  if (verbosity > 2) {
    bft_printf(_("\n  Number of local iterations to converge on vertex"
                 " equivalences: %3d\n"), _loc_merge_counter);
    if (n_ranks > 1)
      bft_printf(_("  Number of global iterations to converge on vertex"
                   " equivalences: %3d\n\n"), _glob_merge_counter);
    bft_printf_flush();
  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  for (i = 0; i < n_vertices; i++)
    bft_printf(" Final vtx_tag[%6d] = %9u\n", i, vtx_tag[i]);
  bft_printf_flush();
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
                               const fvm_gnum_t         vtx_tag[],
                               cs_int_t                 send_count[],
                               cs_int_t                 send_shift[],
                               cs_int_t                 recv_count[],
                               cs_int_t                 recv_shift[],
                               cs_join_vertex_t        *p_vtx_merge_data[],
                               cs_join_gset_t         **p_merge_set)
{
  cs_int_t  i;

  cs_int_t  n_vertices = work->n_vertices;
  fvm_gnum_t  *recv_gbuf = NULL, *send_gbuf = NULL;
  cs_join_vertex_t  *send_vtx_data = NULL, *recv_vtx_data = NULL;
  cs_join_gset_t  *merge_set = NULL;

  MPI_Datatype  CS_MPI_JOIN_VERTEX = cs_join_mesh_create_vtx_datatype();
  MPI_Comm  mpi_comm = cs_glob_mpi_comm;

  const int  n_ranks = cs_glob_n_ranks;

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (i = 0; i < n_vertices; i++) {
    int  rank = (vtx_tag[i] - 1) % n_ranks;
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

  BFT_MALLOC(send_gbuf, send_shift[n_ranks], fvm_gnum_t);
  BFT_MALLOC(recv_gbuf, recv_shift[n_ranks], fvm_gnum_t);

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (i = 0; i < n_vertices; i++) {

    int  rank = (vtx_tag[i] - 1) % n_ranks;
    cs_int_t  shift = send_shift[rank] + send_count[rank];

    send_gbuf[shift] = vtx_tag[i];
    send_count[rank] += 1;

  }

  MPI_Alltoallv(send_gbuf, send_count, send_shift, FVM_MPI_GNUM,
                recv_gbuf, recv_count, recv_shift, FVM_MPI_GNUM,
                mpi_comm);

  /* Allocate and build send_vtx_data, receive recv_vtx_data. */

  BFT_MALLOC(recv_vtx_data, recv_shift[n_ranks], cs_join_vertex_t);
  BFT_MALLOC(send_vtx_data, send_shift[n_ranks], cs_join_vertex_t);

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (i = 0; i < n_vertices; i++) {

    int  rank = (vtx_tag[i] - 1) % n_ranks;
    cs_int_t  shift = send_shift[rank] + send_count[rank];

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
  bft_printf("\n  Number of vertices to treat for the merge step: %d\n",
             recv_shift[n_ranks]);
  bft_printf("  List of vertices to treat:\n");
  for (i = 0; i < recv_shift[n_ranks]; i++) {
    bft_printf(" %9d - ", i);
    cs_join_mesh_dump_vertex(recv_vtx_data[i]);
  }
  bft_printf_flush();
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
                          const fvm_gnum_t       vtx_tag[],
                          cs_int_t               send_count[],
                          cs_int_t               send_shift[],
                          cs_int_t               recv_count[],
                          cs_int_t               recv_shift[],
                          cs_join_vertex_t       vtx_merge_data[])
{
  cs_int_t  i;

  cs_join_vertex_t  *updated_vtx_data = NULL;
  int  n_ranks = cs_glob_n_ranks;
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

    int  rank = (vtx_tag[i] - 1) % n_ranks;
    cs_int_t  shift = send_shift[rank] + send_count[rank];

    work->vertices[i] = updated_vtx_data[shift];
    send_count[rank] += 1;

  }

  /* Free memory */

  MPI_Type_free(&cs_mpi_join_vertex);
  BFT_FREE(updated_vtx_data);
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Check if all vertices in the list include the target_vertex in their
 * tolerance.
 *
 * If check is ok, no need to apply a reduction of the tolerance.
 *
 * parameters:
 *   start      <-- index for the last vertex in id_lst
 *   end        <-- index for the last vertex in id_lst
 *   list       <-- list of id in vertices
 *   vertices   <-- pointer to an array of cs_join_vertex_t structures
 *   ref_vertex <-- vertex resulting of the merge
 *
 * returns:
 *   true if all vertices have ref_vertex in their tolerance, false otherwise
 *---------------------------------------------------------------------------*/

static cs_bool_t
_is_in_tolerance(cs_int_t                start,
                 cs_int_t                end,
                 const fvm_gnum_t        list[],
                 const cs_join_vertex_t  vertices[],
                 cs_join_vertex_t        ref_vertex)
{
  cs_int_t  i;

  for (i = start; i < end; i++) {

    cs_join_vertex_t  cur_vertex = vertices[list[i]];
    cs_real_t  d_cur_ref = _compute_length(cur_vertex, ref_vertex);
    cs_real_t  tolerance =  cur_vertex.tolerance * cs_join_tol_eps_coef2;

    if (d_cur_ref > tolerance)
      return false;

  } /* End of loop on each vertex of the list */

  return true;
}

/*----------------------------------------------------------------------------
 * Get the resulting cs_join_vertex_t structure after the merge of a list
 * of vertices.
 *
 * parameters:
 *   ref_gnum <-- global number associated to the current list of vertices
 *   start    <-- index for the first vertex in id_list
 *   end      <-- index for the last vertex in id_list
 *   list     <-- list of id in vertices array
 *   vertices <-- array of cs_join_vertex_t structures
 *
 * returns:
 *   a cs_join_vertex_t structure for the resulting vertex
 *---------------------------------------------------------------------------*/

static cs_join_vertex_t
_get_merged_vertex(fvm_gnum_t              ref_gnum,
                   cs_int_t                start,
                   cs_int_t                end,
                   const fvm_gnum_t        list[],
                   const cs_join_vertex_t  vertices[])
{
  cs_int_t  i, k;
  cs_join_vertex_t  merged_vertex;

  cs_int_t  n_elts = end - start;

  /* Initialize cs_join_vertex_t structure */

  merged_vertex.gnum = ref_gnum;
  merged_vertex.tolerance = vertices[list[start]].tolerance;

  /* Compute the resulting vertex data of the merge */

  for (i = start; i < end; i++)
    merged_vertex.tolerance = CS_MIN(vertices[list[i]].tolerance,
                                     merged_vertex.tolerance);

  assert(n_elts > 0);

  /* Compute the resulting coordinates of the merged vertices */

  for (k = 0; k < 3; k++)
    merged_vertex.coord[k] = 0.0;

#if CS_JOIN_MERGE_INV_TOL
  {
    cs_real_t  denum = 0.0;

    for (i = start; i < end; i++)
      denum += 1.0/vertices[list[i]].tolerance;

    for (k = 0; k < 3; k++) {

      for (i = start; i < end; i++)
        merged_vertex.coord[k] +=
          1.0/vertices[list[i]].tolerance * vertices[list[i]].coord[k];

      merged_vertex.coord[k] /= denum;

    } /* End of loop on coordinates */

  }
#else

  for (k = 0; k < 3; k++) {

    for (i = start; i < end; i++)
      merged_vertex.coord[k] += vertices[list[i]].coord[k];
    merged_vertex.coord[k] /= n_elts;

  }
#endif

  return merged_vertex;
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
_trivial_merge(cs_join_param_t     param,
               cs_join_gset_t     *merge_set,
               cs_join_vertex_t    vertices[],
               cs_join_gset_t    **p_equiv_gnum)
{
  cs_int_t  i, j, j1, j2, k, k1, k2, n_sub_elts;
  cs_real_t  delta;

  cs_int_t  max_n_sub_elts = 0;
  cs_int_t  *merge_index = merge_set->index;
  fvm_gnum_t  *merge_list = merge_set->g_list;
  fvm_gnum_t  *sub_list = NULL, *init_list = NULL;
  cs_join_gset_t  *equiv_gnum = NULL;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (param.verbosity > 2) {

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

  BFT_MALLOC(sub_list, max_n_sub_elts, fvm_gnum_t);

  /* Store initial merge list */

  BFT_MALLOC(init_list, merge_index[merge_set->n_elts], fvm_gnum_t);

  for (i = 0; i < merge_index[merge_set->n_elts]; i++)
    init_list[i] = merge_list[i];

  /* Apply merge */

  for (i = 0; i < merge_set->n_elts; i++) {

    cs_int_t  f_s = merge_index[i];
    cs_int_t  f_e = merge_index[i+1];

    n_sub_elts = f_e - f_s;

    for (j = f_s, k = 0; j < f_e; j++, k++)
      sub_list[k] = merge_list[j];

    for (j1 = 0; j1 < n_sub_elts - 1; j1++) {

      cs_int_t  v1_id = sub_list[j1];
      cs_join_vertex_t  v1 = vertices[v1_id];

      for (j2 = j1 + 1; j2 < n_sub_elts; j2++) {

        cs_int_t  v2_id = sub_list[j2];
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

          if (fabs(v1.tolerance - v2.tolerance) < 1e-30) {

            delta = 0.0;
            for (k = 0; k < 3; k++)
              delta += v1.coord[k] - v2.coord[k];

            if (fabs(delta) < 1e-30) { /* Identical vertices */

              if (v1.gnum < v2.gnum)
                k1 = j1, k2 = j2;
              else
                k1 = j2, k2 = j1;

              for (k = 0; k < n_sub_elts; k++)
                if (sub_list[k] == sub_list[k2])
                  sub_list[k] = sub_list[k1];

            } /* End if delta == 0.0 */

          } /* End if v1.tol != v2.tol */

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

  /* Free memory */

  BFT_FREE(sub_list);
  BFT_FREE(init_list);

  /* Return pointer */

  *p_equiv_gnum = equiv_gnum;
}

/*----------------------------------------------------------------------------
 * Apply the tolerance reduction until each vertex of the list has
 * a resulting vertex under its tolerance.
 *
 * parameters:
 *   param        <-- set of user-defined parameters
 *   start        <-- index for the last vertex in sub_list
 *   end          <-- index for the last vertex in sub_list
 *   list         <-- list of id in vertices
 *   vertices     <-> pointer to a cs_join_vertex_t structure
 *   distances    <-- list of distance between couples of vertices
 *   ref_tags     <-> list of reference tags for each vertex of sub_list
 *   work_tags    <-> list of working tags for each vertex of sub_list
 *   bool_lst     <-> list of booleans
 *   n_reductions <-> number of tolerance reduction done
 *---------------------------------------------------------------------------*/

static void
_reduce_tolerance(cs_join_param_t     param,
                  cs_int_t            start,
                  cs_int_t            end,
                  const fvm_gnum_t    list[],
                  cs_join_vertex_t    vertices[],
                  const cs_real_t     distances[],
                  fvm_gnum_t          ref_tags[],
                  fvm_gnum_t          work_tags[],
                  cs_bool_t           bool_lst[],
                  int                *n_reductions)
{
  cs_int_t  i;

  cs_bool_t reduc_tol = true;

  while (reduc_tol == true) {

    cs_int_t  i1, i2, shift;

    *n_reductions += 1;

    if (*n_reductions > CS_JOIN_MERGE_MAX_REDUCTIONS)
      bft_error(__FILE__, __LINE__, 0,
                _("  Max number of tolerance reductions has been reached.\n"
                  "  Check your joining parameters.\n"));

    /* Reduce tolerance by the constant "reduce_tol_factor" */

    for (i = start; i < end; i++)
      vertices[list[i]].tolerance *= param.reduce_tol_factor;

    /* Define a boolean list on couples of vertices which are under
       tolerance each other */

    for (shift = 0, i1 = start; i1 < end - 1; i1++) {

      cs_join_vertex_t  v1 = vertices[list[i1]];

      for (i2 = i1 + 1; i2 < end; i2++) {

        cs_join_vertex_t  v2 = vertices[list[i2]];

        if (v2.tolerance < distances[shift] || v1.tolerance < distances[shift])
          bool_lst[shift++] = false;
        else
          bool_lst[shift++] = true;

      }

    }

    /* Update tag list according to the values of bool_lst */

    for (i = start; i < end; i++)
      work_tags[i] = ref_tags[i];

    for (shift = 0, i1 = start; i1 < end - 1; i1++) {
      for (i2 = i1 + 1; i2 < end; i2++) {

        if (bool_lst[shift] == true) {

          fvm_gnum_t  _min = CS_MIN(work_tags[i1], work_tags[i2]);

          work_tags[i1] = _min;
          work_tags[i2] = _min;

        }
        shift++;

      } /* End of loop on i2 */
    } /* End of loop on i1 */

    /* Check if the reduction of the tolerance is sufficient */

    reduc_tol = false;

    for (shift = 0, i1 = start; i1 < end - 1; i1++) {
      for (i2 = i1 + 1; i2 < end; i2++) {

        if (bool_lst[shift] == true) {
          if (work_tags[i1] != work_tags[i2])
            reduc_tol = true;
        }
        else  {/* bool_lst = false */
          if (work_tags[i1] == work_tags[i2])
            reduc_tol = true;
        }
        shift++;

      } /* End of loop on i2 */
    } /* End of loop on i1 */

  } /* End of while reduc_tol = true */

  /* Store new equivalences between vertices in ref_tags */

  for (i = start; i < end; i++) {
    ref_tags[i] = work_tags[i];
    work_tags[i] = i;
  }

  { /* Order ref_tags and keep the original position in work_tags */

    int  h, j;
    cs_int_t  n_sub_elts = end - start;

    /* Compute stride */
    for (h = 1; h <= n_sub_elts/9; h = 3*h+1) ;

    /* Sort array */
    for ( ; h > 0; h /= 3) {

      for (i = start + h; i < end; i++) {

        fvm_gnum_t vr = ref_tags[i];
        fvm_gnum_t vw = work_tags[i];

        j = i;
        while ( (j >= start+h) && (vr < ref_tags[j-h]) ) {
          ref_tags[j] = ref_tags[j-h];
          work_tags[j] = work_tags[j-h];
          j -= h;
        }
        ref_tags[j] = vr;
        work_tags[j] = vw;

      } /* Loop on array elements */

    } /* End of loop on stride */

  } /* End of sort */

  for (i = start; i < end; i++)
    work_tags[i] = list[work_tags[i]];

}

/*----------------------------------------------------------------------------
 * Apply a reduction of the tolerance until each vertex of the list has
 * the resulting vertex of the merge under its tolerance.
 *
 * parameters:
 *   param         <-- set of user-defined parameters
 *   start         <-- index for the last vertex in sub_list
 *   end           <-- index for the last vertex in sub_list
 *   sub_list      <-> list of id in vertices for each vertex of the sub list
 *   vertices      <-> pointer to a cs_join_vertex_t structure
 *   distances     <-> list of distance between couples of vertices
 *   ref_tag_list  <-> list of reference tags for each vertex of sub_list
 *   work_tag_list <-> list of working tags for each vertex of sub_list
 *   bool_lst      <-> list of booleans
 *   n_reductions  <-> number of tolerance reduction done
 *---------------------------------------------------------------------------*/

static void
_merge_with_tol_reduction(cs_join_param_t    param,
                          cs_int_t           start,
                          cs_int_t           end,
                          fvm_gnum_t         sub_list[],
                          cs_join_vertex_t   vertices[],
                          cs_real_t          distances[],
                          fvm_gnum_t         ref_tags[],
                          fvm_gnum_t         work_tags[],
                          cs_bool_t          bool_list[],
                          int               *n_reductions)
{
  cs_int_t  i1, i2, i, j;

  cs_int_t  shift = 0;
  cs_int_t  n_sub_equiv = 1;
  cs_int_t  n_sub_elts = end - start;
  cs_int_t  *sub_index = NULL;

  /* Sanity checks */

  assert(param.reduce_tol_factor >= 0.0);
  assert(param.reduce_tol_factor < 1.0);

  if (param.verbosity > 2) { /* Display information */

    bft_printf(_("\n"
                 "  Merge: Local tolerance reduction for the"
                 " following vertices (global numbering):\n"
                 "        "));
    for (j = start; j < end; j++)
      bft_printf("%u  ", vertices[sub_list[j]].gnum);
    bft_printf("\n");
  }

  /* Initialize temporary buffers */

  for (i = start; i < end; i++)
    ref_tags[i] = vertices[sub_list[i]].gnum;

  for (i1 = start; i1 < end - 1; i1++)
    for (i2 = i1 + 1; i2 < end; i2++)
      distances[shift++] = _compute_length(vertices[sub_list[i1]],
                                           vertices[sub_list[i2]]);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\t\t\t  (BEGIN) Reduce tolerance\n");
  cs_join_dump_array("gnum", "sub_list", n_sub_elts, &(sub_list[start]));
  cs_join_dump_array("gnum", "ref_tags", n_sub_elts, &(ref_tags[start]));
  shift = 0;
  bft_printf("\n  distances:\n");
  for (i1 = start; i1 < end - 1; i1++) {
    bft_printf(" %10u : ", ref_tags[i1]);
    for (i2 = i1+1; i2 < end; i2++)
      bft_printf(" (%10u, %le)", ref_tags[i2], distances[shift++]);
    bft_printf("\n");
  }
  bft_printf("\n");
  bft_printf_flush();
#endif

  /* Reduce tolerance until coherent sub-equivalences appear.
     Local operation. */

  _reduce_tolerance(param,
                    start,
                    end,
                    sub_list,
                    vertices,
                    distances,
                    ref_tags,
                    work_tags,
                    bool_list,
                    n_reductions);

  /* Define an index on new sub equivalences */

  BFT_MALLOC(sub_index, n_sub_elts+1, cs_int_t);

  sub_index[0] = start;
  sub_index[n_sub_equiv] = start + 1;

  for (i = start + 1; i < end; i++) {

    if (ref_tags[i] != ref_tags[i-1]) {
      sub_index[n_sub_equiv + 1] = sub_index[n_sub_equiv] + 1;
      n_sub_equiv++;
    }
    else
      sub_index[n_sub_equiv] += 1;

  }

  assert(n_sub_equiv > 1);

  /* Loop on sub-equivalences */

  for (i = 0; i < n_sub_equiv; i++) {

    cs_bool_t  ok;
    cs_join_vertex_t  merged_vertex;

    cs_int_t  sub_start = sub_index[i];
    cs_int_t  sub_end = sub_index[i+1];

    /* Define a sub_list */

    for (j = sub_start; j < sub_end; j++)
      sub_list[j] = work_tags[j];

    if (sub_end - sub_start > 1) {

      merged_vertex = _get_merged_vertex(ref_tags[sub_start],
                                         sub_start,
                                         sub_end,
                                         sub_list,
                                         vertices);

      /* Check if the vertex resulting of the merge is in the tolerance
         for each vertex of the list. Previous op. should assume this. */

      ok = _is_in_tolerance(sub_start,
                            sub_end,
                            sub_list,
                            vertices,
                            merged_vertex);

      assert(ok == true); /* One call to _reduce_tolerance should be enough */

      for (j = sub_start; j < sub_end; j++)
        vertices[sub_list[j]] = merged_vertex;

    } /* n_sub_sub_list > 1 */

  } /* End of loop on sub-equivalences */

  BFT_FREE(sub_index);

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
                cs_int_t           n_vertices,
                cs_join_vertex_t   vertices[])
{
  cs_int_t  i, j, k, tmp_size, n_sub_elts;
  cs_join_vertex_t  merged_vertex;
  cs_bool_t  ok;

  int  n_reduction_counter = 0, n_max_tol_reductions = 0;
  cs_int_t  max_n_sub_elts = 0;

  cs_join_gset_t  *equiv_gnum = NULL;
  cs_bool_t  *bool_list = NULL;
  cs_int_t  *merge_index = NULL;
  fvm_gnum_t  *merge_list = NULL, *merge_ref_elts = NULL;
  fvm_gnum_t  *sub_list = NULL, *tmp_buffer = NULL;
  fvm_gnum_t  *ref_tags = NULL, *work_tags = NULL;
  cs_real_t  *distances = NULL;

  const int  verbosity = param.verbosity;

  /* Sanity check */

  assert(param.merge_tol_coef >= 0.0);

  /* Pre-merge of identical vertices */

  _trivial_merge(param, merge_set, vertices, &equiv_gnum);

  merge_index = merge_set->index;
  merge_list = merge_set->g_list;
  merge_ref_elts = merge_set->g_elts;

#if 0 && defined(DEBUG) && !defined(NDEBUG)

  if (verbosity > 2) {

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

  /* Compute the max. size of a sub list */

  for (i = 0; i < merge_set->n_elts; i++) {
    n_sub_elts = merge_index[i+1] - merge_index[i];
    max_n_sub_elts = CS_MAX(max_n_sub_elts, n_sub_elts);
  }

  /* Modify the tolerance for the merge operation if needed */

  if (fabs(param.merge_tol_coef - 1.0) < 1e-30) {
    for (i = 0; i < n_vertices; i++)
      vertices[i].tolerance *= param.merge_tol_coef;
  }

  /* Compute the max. size of a sub list */

  for (i = 0; i < merge_set->n_elts; i++) {
    n_sub_elts = merge_index[i+1] - merge_index[i];
    max_n_sub_elts = CS_MAX(max_n_sub_elts, n_sub_elts);
  }

  if (verbosity > 0) {   /* Display information */

    fvm_lnum_t g_max_n_sub_elts = max_n_sub_elts;
    fvm_parall_counter_max(&g_max_n_sub_elts, 1);

    if (g_max_n_sub_elts < 2) {
      bft_printf(_("\n  No need to merge vertices.\n"));
      return;
    }
    else
      bft_printf(_("\n  Max size of a merge list: %lu\n"),
                 (unsigned long)g_max_n_sub_elts);
  }

  /* Temporary buffers */

  BFT_MALLOC(tmp_buffer, 3*max_n_sub_elts, fvm_gnum_t);

  sub_list = tmp_buffer;
  ref_tags = &(tmp_buffer[max_n_sub_elts]);
  work_tags = &(tmp_buffer[2*max_n_sub_elts]);

  tmp_size = ((max_n_sub_elts-1)*max_n_sub_elts)/2;
  BFT_MALLOC(distances, tmp_size, cs_real_t);
  BFT_MALLOC(bool_list, tmp_size, cs_bool_t);

  /* Apply merge */

  for (i = 0; i < merge_set->n_elts; i++) {

    n_sub_elts = merge_index[i+1] - merge_index[i];

    if (n_sub_elts > 1) {

      for (k = merge_index[i], j = 0; k < merge_index[i+1]; k++, j++)
        sub_list[j] = merge_list[k];

      if (verbosity > 2) { /* Display information */

        bft_printf("\n Begin merge for ref. elt: %u - n_sub_elts: %d\n",
                   merge_ref_elts[i], merge_index[i+1] - merge_index[i]);
        for (j = 0; j < n_sub_elts; j++) {
          bft_printf("%9u -", sub_list[j]);
          cs_join_mesh_dump_vertex(vertices[sub_list[j]]);
        }
        bft_printf("\n");

      }

      /* Define the resulting cs_join_vertex_t structure of the merge */

      merged_vertex = _get_merged_vertex(merge_ref_elts[i],
                                         0, n_sub_elts,
                                         sub_list,
                                         vertices);

      /* Check if the vertex resulting of the merge is in the tolerance
         for each vertex of the list */

      ok = _is_in_tolerance(0, n_sub_elts, sub_list, vertices, merged_vertex);

#if CS_JOIN_MERGE_TOL_REDUC

      if (ok == false) { /* Reduction of the tolerance is necessary */

        int  n_tol_reductions = 0;

        n_reduction_counter += 1;

        _merge_with_tol_reduction(param,
                                  0,
                                  n_sub_elts,
                                  sub_list,
                                  vertices,
                                  distances,
                                  ref_tags,
                                  work_tags,
                                  bool_list,
                                  &n_tol_reductions);

        if (verbosity > 1) /* Display information */
          bft_printf(_("\n  Number of tolerance reductions: %4d\n"),
                     n_tol_reductions);

        n_max_tol_reductions = CS_MAX(n_max_tol_reductions,
                                      n_tol_reductions);

      }

      else /* New vertex data for the sub-elements */

#endif /* CS_JOIN_MERGE_TOL_REDUC */

        for (j = 0; j < n_sub_elts; j++)
          vertices[sub_list[j]] = merged_vertex;

      if (verbosity > 2) { /* Display information */

        bft_printf("\n End merge for ref. elt: %u - n_sub_elts: %d\n",
                   merge_ref_elts[i], merge_index[i+1] - merge_index[i]);
        for (j = 0; j < n_sub_elts; j++) {
          bft_printf("%7u -", sub_list[j]);
          cs_join_mesh_dump_vertex(vertices[sub_list[j]]);
        }
        bft_printf("\n");

      }

    } /* sub_list_size > 1 */

  } /* End of loop on potential merges */

  /* Apply merge to vertex initially identical */

  if (equiv_gnum != NULL) {

#if 0 && defined(DEBUG) && !defined(NDEBUG)

    if (verbosity > 2) {

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

      cs_int_t  start = equiv_gnum->index[i];
      cs_int_t  end = equiv_gnum->index[i+1];
      cs_int_t  ref_id = equiv_gnum->g_elts[i];

      for (j = start; j < end; j++)
        vertices[equiv_gnum->g_list[j]] = vertices[ref_id];

    }
  }

  if (verbosity > 0) {

    fvm_gnum_t n_g_reduction_counter = n_reduction_counter;
    fvm_parall_counter(&n_g_reduction_counter, 1);

    bft_printf(_("\n  Tolerance reduction for %lu elements.\n"),
               (unsigned long)n_g_reduction_counter);

    if (verbosity > 1) {
      fvm_lnum_t g_n_max_tol_reductions = n_max_tol_reductions;
      fvm_parall_counter_max(&g_n_max_tol_reductions, 1);
      bft_printf(_("\n  Max. number of tolerance reductions: %lu\n"),
                 (unsigned long)g_n_max_tol_reductions);
    }
  }

  /* Free memory */

  BFT_FREE(tmp_buffer);
  BFT_FREE(distances);
  BFT_FREE(bool_list);

  cs_join_gset_destroy(&equiv_gnum);
}

/*----------------------------------------------------------------------------
 * Keep an history of the evolution of each vertex id before/after the merge
 * operation.
 *
 * parameters:
 *   n_iwm_vertices   <-- number of vertices before intersection for the work
 *                        cs_join_mesh_t structure
 *   n_g_ifm_vertices <-- global number of vertices on the initial full mesh
 *   iwm_vtx_gnum     <-- initial global vertex num. (work mesh struct.)
 *   n_vertices       <-- number of vertices before merge/after intersection
 *   vertices         <-- array of cs_join_vertex_t structures
 *   p_o2n_vtx_gnum   --> distributed array by block on the new global vertex
 *                        numbering for the initial vertices (before inter.)
 *---------------------------------------------------------------------------*/

static void
_keep_global_vtx_evolution(cs_int_t                n_iwm_vertices,
                           fvm_gnum_t              n_g_ifm_vertices,
                           const fvm_gnum_t        iwm_vtx_gnum[],
                           cs_int_t                n_vertices,
                           const cs_join_vertex_t  vertices[],
                           fvm_gnum_t             *p_o2n_vtx_gnum[])
{
  cs_int_t  i;
  cs_join_block_info_t  block_info;

  int  n_ranks = cs_glob_n_ranks;
  fvm_gnum_t  *o2n_vtx_gnum = NULL;

  const int  local_rank = CS_MAX(cs_glob_rank_id, 0);

  assert(n_iwm_vertices <= n_vertices); /* after inter. >= init */

  block_info = cs_join_get_block_info(n_g_ifm_vertices,
                                      n_ranks,
                                      local_rank);

  BFT_MALLOC(o2n_vtx_gnum, block_info.local_size, fvm_gnum_t);

  if (n_ranks == 1) {

    for (i = 0; i < n_iwm_vertices; i++)
      o2n_vtx_gnum[i] = vertices[i].gnum;

    /* Return pointer */

    *p_o2n_vtx_gnum = o2n_vtx_gnum;

    return;
  }

#if defined(HAVE_MPI) /* Parallel treatment */
  {
    fvm_gnum_t  ii;
    cs_int_t  shift, rank, n_recv_elts;

    cs_int_t  *send_shift = NULL, *recv_shift = NULL;
    cs_int_t  *send_count = NULL, *recv_count = NULL;
    fvm_gnum_t  *send_glist = NULL, *recv_glist = NULL;

    MPI_Comm  mpi_comm = cs_glob_mpi_comm;

    /* Initialize o2n_vtx_gnum */

    for (ii = 0; ii < block_info.local_size; ii++)
      o2n_vtx_gnum[ii] = block_info.first_gnum + ii;

    /* Send new vtx global number to the related rank = the good block */

    BFT_MALLOC(send_count, n_ranks, cs_int_t);
    BFT_MALLOC(recv_count, n_ranks, cs_int_t);

    for (i = 0; i < n_ranks; i++)
      send_count[i] = 0;

    for (i = 0; i < n_iwm_vertices; i++) {
      rank = (iwm_vtx_gnum[i] - 1)/block_info.size;
      send_count[rank] += 2;
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

    assert(send_shift[n_ranks] == 2*n_iwm_vertices);

    /* Build send_list */

    BFT_MALLOC(send_glist, send_shift[n_ranks], fvm_gnum_t);
    BFT_MALLOC(recv_glist, recv_shift[n_ranks], fvm_gnum_t);

    for (i = 0; i < n_ranks; i++)
      send_count[i] = 0;

    for (i = 0; i < n_iwm_vertices; i++) {

      rank = (iwm_vtx_gnum[i] - 1)/block_info.size;
      shift = send_shift[rank] + send_count[rank];

      send_glist[shift] = iwm_vtx_gnum[i];  /* Old global number */
      send_glist[shift+1] = vertices[i].gnum;   /* New global number */
      send_count[rank] += 2;

    }

    MPI_Alltoallv(send_glist, send_count, send_shift, FVM_MPI_GNUM,
                  recv_glist, recv_count, recv_shift, FVM_MPI_GNUM,
                  mpi_comm);

    n_recv_elts = recv_shift[n_ranks]/2;

    BFT_FREE(send_count);
    BFT_FREE(send_shift);
    BFT_FREE(send_glist);
    BFT_FREE(recv_count);

    /* Update o2n_vtx_gnum */

    for (rank = 0; rank < n_ranks; rank++) {

      for (i = recv_shift[rank]; i < recv_shift[rank+1]; i+=2) {

        fvm_gnum_t  o_gnum = recv_glist[i];
        fvm_gnum_t  n_gnum = recv_glist[i+1];
        cs_int_t  id = o_gnum - block_info.first_gnum;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
        if (o2n_vtx_gnum[id] != block_info.first_gnum + id)
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
 *   p_n_af_vertices --> number of vertices after the merge step
 *   p_o2n_vtx_id    --> array keeping the evolution of the vertex ids
 *---------------------------------------------------------------------------*/

static void
_keep_local_vtx_evolution(cs_int_t                 n_vertices,
                          const cs_join_vertex_t   vertices[],
                          cs_int_t                *p_n_af_vertices,
                          cs_int_t                *p_o2n_vtx_id[])
{
  cs_int_t  i;
  fvm_gnum_t  prev;

  cs_int_t  n_af_vertices = 0;
  cs_int_t  *o2n_vtx_id = NULL;
  fvm_lnum_t  *order = NULL;
  fvm_gnum_t  *vtx_gnum = NULL;

  if (n_vertices == 0)
    return;

  BFT_MALLOC(vtx_gnum, n_vertices, fvm_gnum_t);

  for (i = 0; i < n_vertices; i++)
    vtx_gnum[i] = vertices[i].gnum;

  /* Order vertices according to their global numbering */

  BFT_MALLOC(order, n_vertices, fvm_lnum_t);

  fvm_order_local_allocated(NULL, vtx_gnum, order, n_vertices);

  /* Delete vertices sharing the same global number. Keep only one */

  BFT_MALLOC(o2n_vtx_id, n_vertices, cs_int_t);

  prev = vtx_gnum[order[0]];
  o2n_vtx_id[order[0]] = n_af_vertices;

  for (i = 1; i < n_vertices; i++) {

    cs_int_t  o_id = order[i];
    fvm_gnum_t  cur = vtx_gnum[o_id];

    if (cur != prev) {
      prev = cur;
      n_af_vertices++;
      o2n_vtx_id[o_id] = n_af_vertices;
    }
    else
      o2n_vtx_id[o_id] = n_af_vertices;

  } /* End of loop on vertices */

  /* n_af_vertices is an id */
  n_af_vertices += 1;

  assert(n_af_vertices <= n_vertices); /* after merge <= after inter. */

  /* Free memory */

  BFT_FREE(order);
  BFT_FREE(vtx_gnum);

  /* Set return pointers */

  *p_n_af_vertices = n_af_vertices;
  *p_o2n_vtx_id = o2n_vtx_id;
}

/*----------------------------------------------------------------------------
 * Update a cs_join_inter_edges_t structure after the merge operation.
 * cs_join_inter_edges_t structure should be not NULL.
 *
 * parameters:
 *   o2n_vtx_id    <-- array keeping the evolution of the vertex ids
 *   p_inter_edges <-> pointer to the structure keeping data on
 *                     edge intersections
 *---------------------------------------------------------------------------*/

static void
_update_inter_edges_after_merge(const cs_int_t           o2n_vtx_id[],
                                cs_join_inter_edges_t  **p_inter_edges)
{
  cs_int_t  i, j, prev_num, new_num;

  cs_int_t  *new_index = NULL;
  cs_join_inter_edges_t  *inter_edges = *p_inter_edges;
  cs_int_t  n_edges = inter_edges->n_edges;

  /* No need to update vtx_glst because it's no more used */

  if (inter_edges->vtx_glst != NULL)
    BFT_FREE(inter_edges->vtx_glst);

  /* Update cs_join_inter_edges_t structure */

  for (i = 0; i < n_edges; i++) {

    cs_int_t  start = inter_edges->index[i];
    cs_int_t  end = inter_edges->index[i+1];

    for (j = start; j < end; j++) {

      cs_int_t old_id = inter_edges->vtx_lst[j] - 1;

      inter_edges->vtx_lst[j] = o2n_vtx_id[old_id] + 1;

    }

  }

  /* Delete redundancies and define a new index */

  BFT_MALLOC(new_index, n_edges + 1, cs_int_t);

  for (i = 0; i < n_edges + 1; i++)
    new_index[i] = 0;

  for (i = 0; i < n_edges; i++) {

    cs_int_t  start = inter_edges->index[i];
    cs_int_t  end = inter_edges->index[i+1];

    if (end - start > 0) {

      prev_num = inter_edges->vtx_lst[start];
      new_index[i+1] += 1;

      for (j = start + 1; j < end; j++) {

        new_num = inter_edges->vtx_lst[j];

        if (prev_num != new_num) {
          prev_num = new_num;
          new_index[i+1] += 1;
        }

      }

    } /* end - start > 0 */

  } /* End of loop on edge intersections */

  inter_edges->max_sub_size = 0;

  for (i = 0; i < n_edges; i++) {

    inter_edges->max_sub_size = CS_MAX(inter_edges->max_sub_size,
                                       new_index[i+1]);
    new_index[i+1] += new_index[i];

  }

  for (i = 0; i < n_edges; i++) {

    cs_int_t  start = inter_edges->index[i];
    cs_int_t  end = inter_edges->index[i+1];

    if (end - start > 0) {

      cs_int_t  shift = new_index[i];

      inter_edges->vtx_lst[shift] = inter_edges->vtx_lst[start];
      inter_edges->abs_lst[shift++] = inter_edges->abs_lst[start];

      prev_num = inter_edges->vtx_lst[start];

      for (j = start + 1; j < end; j++) {

        new_num = inter_edges->vtx_lst[j];

        if (prev_num != new_num) {
          inter_edges->vtx_lst[shift] = new_num;
          inter_edges->abs_lst[shift++] = inter_edges->abs_lst[j];
          prev_num = new_num;
        }

      }

    } /* end - start  > 0 */

  } /* End of loop on edge intersections */

  /* Apply updates to the structure */

  BFT_FREE(inter_edges->index);
  inter_edges->index = new_index;

  BFT_REALLOC(inter_edges->vtx_lst, inter_edges->index[n_edges], cs_int_t);
  BFT_REALLOC(inter_edges->abs_lst, inter_edges->index[n_edges], float);

  /* Return pointer */

  *p_inter_edges = inter_edges;
}

/*----------------------------------------------------------------------------
 * Define send_rank_index and send_faces to prepare the exchange of new faces
 * between mesh structures.
 *
 * parameters:
 *   n_faces           <-- number of faces to send
 *   n_g_faces         <-- global number of faces to be joined
 *   face_gnum         <-- global face number
 *   gnum_rank_index   <-- index on ranks for the init. global face numbering
 *   p_send_rank_index --> index on ranks for sending face
 *   p_send_faces      --> list of face ids to send
 *---------------------------------------------------------------------------*/

static void
_get_faces_to_send(cs_int_t           n_faces,
                   fvm_gnum_t         n_g_faces,
                   const fvm_gnum_t   face_gnum[],
                   const fvm_gnum_t   gnum_rank_index[],
                   cs_int_t          *p_send_rank_index[],
                   cs_int_t          *p_send_faces[])
{
  cs_int_t  i, rank, shift;
  fvm_gnum_t  start_gnum, end_gnum;
  cs_join_block_info_t  block_info;

  cs_int_t  reduce_size = 0;
  cs_int_t  *send_rank_index = NULL, *send_faces = NULL;
  cs_int_t  *reduce_ids = NULL, *count = NULL;
  fvm_gnum_t  *reduce_index = NULL;

  const int  local_rank = CS_MAX(cs_glob_rank_id, 0);
  const int  n_ranks = cs_glob_n_ranks;

  /* Sanity checks */

  assert(gnum_rank_index != NULL);
  assert(n_ranks > 1);

  /* Compute block_size */

  block_info = cs_join_get_block_info(n_g_faces, n_ranks, local_rank);
  start_gnum = block_info.first_gnum;
  end_gnum = block_info.first_gnum + block_info.local_size;

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

  BFT_MALLOC(send_rank_index, n_ranks + 1, cs_int_t);

  for (i = 0; i < n_ranks + 1; i++)
    send_rank_index[i] = 0;

  /* Count number of ranks associated to each new face */

  for (i = 0; i < n_faces; i++) {

    if (face_gnum[i] >= start_gnum && face_gnum[i] < end_gnum) {

      /* The current face is a "main" face for the local rank */

      int  reduce_rank = cs_search_gindex_binary(reduce_size,
                                                 face_gnum[i],
                                                 reduce_index);

      assert(reduce_rank != -1);
      assert(reduce_rank < reduce_size);

      rank = reduce_ids[reduce_rank];
      send_rank_index[rank+1] += 1;

    }

  }

  for (i = 0; i < n_ranks; i++)
    send_rank_index[i+1] += send_rank_index[i];

  BFT_MALLOC(send_faces, send_rank_index[n_ranks], cs_int_t);
  BFT_MALLOC(count, n_faces, cs_int_t);

  for (i = 0; i < n_faces; i++)
    count[i] = 0;

  /* Fill the list of ranks */

  for (i = 0; i < n_faces; i++) {

    if (face_gnum[i] >= start_gnum && face_gnum[i] < end_gnum) {

      /* The current face is a "main" face for the local rank */

      int  reduce_rank = cs_search_gindex_binary(reduce_size,
                                                 face_gnum[i],
                                                 reduce_index);

      rank = reduce_ids[reduce_rank];
      shift = send_rank_index[rank] + count[rank];
      send_faces[shift] = i;
      count[rank] += 1;

    } /* End of loop on initial faces */

  }

  /* Free memory */

  BFT_FREE(count);
  BFT_FREE(reduce_ids);
  BFT_FREE(reduce_index);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  for (rank = 0; rank < n_ranks; rank++) {
    bft_printf(" rank %d: ", rank);
    for (i = send_rank_index[rank]; i < send_rank_index[rank+1]; i++)
      bft_printf(" %d (%u)", send_faces[i], face_gnum[send_faces[i]]);
    bft_printf("\n");
    bft_printf_flush();
  }
#endif

  /* Set return pointers */

  *p_send_rank_index = send_rank_index;
  *p_send_faces = send_faces;
}

/*----------------------------------------------------------------------------
 * Update local_mesh by redistributing mesh.
 * Send back to the original rank the new face and vertex description.
 *
 * parameters:
 *   gnum_rank_index <-- index on ranks for the old global face numbering
 *   send_mesh       <-- distributed mesh on faces to join
 *   p_recv_mesh     <-> mesh on local selected faces to be joined
 *---------------------------------------------------------------------------*/

static void
_redistribute_mesh(const fvm_gnum_t        gnum_rank_index[],
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

    cs_int_t  *send_rank_index = NULL, *send_faces = NULL;

    MPI_Comm  mpi_comm = cs_glob_mpi_comm;

    /* Free some structures of the mesh */

    cs_join_mesh_reset(recv_mesh);

    _get_faces_to_send(send_mesh->n_faces,
                       send_mesh->n_g_faces,
                       send_mesh->face_gnum,
                       gnum_rank_index,
                       &send_rank_index,
                       &send_faces);

    assert(send_rank_index[n_ranks] <= send_mesh->n_faces);

    /* Get the new face connectivity from the distributed send_mesh */

    cs_join_mesh_exchange(n_ranks,
                          send_rank_index,
                          send_faces,
                          send_mesh,
                          recv_mesh,
                          mpi_comm);

    BFT_FREE(send_faces);
    BFT_FREE(send_rank_index);

  }
#endif

  /* Return pointers */

  *p_recv_mesh = recv_mesh;

}

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
 *   n_g_vertices       <-- global number of vertices (initial parent mesh)
 *   p_n_g_new_vertices <-> pointer to the global number of new vertices
 *   p_vtx_eset         <-> pointer to a structure dealing with vertex
 *                          equivalences
 *---------------------------------------------------------------------------*/

void
cs_join_create_new_vertices(int                     verbosity,
                            const cs_join_edges_t  *edges,
                            cs_join_mesh_t         *work,
                            cs_join_inter_set_t    *inter_set,
                            fvm_gnum_t              n_g_vertices,
                            fvm_gnum_t             *p_n_g_new_vertices,
                            cs_join_eset_t        **p_vtx_eset)
{
  cs_int_t  i;

  cs_int_t  n_new_vertices = 0;
  fvm_gnum_t  n_g_new_vertices = 0;
  fvm_gnum_t  *new_vtx_gnum = NULL;
  cs_int_t  n_iwm_vertices = work->n_vertices;
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
                           n_g_vertices,
                           n_iwm_vertices,
                           n_new_vertices,
                           &n_g_new_vertices,
                           &new_vtx_gnum);

  if (verbosity > 0)
    bft_printf(_("  Number of new vertices to create: %10lu\n"),
               (unsigned long)n_g_new_vertices);

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

    for (i = 0; i < n_new_vertices; i++)
      work->vertices[n_iwm_vertices + i] = incoherency;

  }
#endif

  /* Fill vertices structure with new vertex definitions */

  for (i = 0; i < inter_set->n_inter; i++) {

    cs_join_inter_t  inter1 = inter_set->inter_lst[2*i];
    cs_join_inter_t  inter2 = inter_set->inter_lst[2*i+1];
    cs_int_t  v1_num = inter1.vtx_id + 1;
    cs_int_t  v2_num = inter2.vtx_id + 1;
    cs_int_t  equiv_id = vtx_equiv->n_equiv;

    assert(inter1.vtx_id < work->n_vertices);
    assert(inter2.vtx_id < work->n_vertices);

    /* Create new vertices if needed */

    if (v1_num > n_iwm_vertices) { /* Add a new vertex */

      cs_int_t shift = inter1.vtx_id - n_iwm_vertices;
      cs_join_vertex_t  new = _get_new_vertex(inter1.curv_abs,
                                               new_vtx_gnum[shift],
                                               &(edges->def[2*inter1.edge_id]),
                                               work);

      work->vertices[inter1.vtx_id] = new;

    }

    if (v2_num > n_iwm_vertices) { /* Add a new vertex */

      cs_int_t shift = inter2.vtx_id - n_iwm_vertices;
      cs_join_vertex_t  new = _get_new_vertex(inter2.curv_abs,
                                               new_vtx_gnum[shift],
                                               &(edges->def[2*inter2.edge_id]),
                                               work);

      work->vertices[inter2.vtx_id] = new;

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
                  "      %u - [%7.4le, %7.4le, %7.4le] - %lg\n"),
                i, vtx.gnum, vtx.coord[0], vtx.coord[1], vtx.coord[2],
                vtx.tolerance);

  } /* End of loop on vertices */

#if 0
  if (verbosity > 3)  /* Dump local structures */
    _dump_vtx_eset(vtx_equiv, work);
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
                       fvm_gnum_t             n_g_vertices_tot,
                       cs_join_mesh_t        *work,
                       const cs_join_eset_t  *vtx_eset)
{
  double  clock_start, clock_end, cpu_start, cpu_end;

  fvm_gnum_t  *vtx_tags = NULL;
  cs_join_gset_t  *merge_set = NULL;

  const int  n_ranks = cs_glob_n_ranks;

  /* Initialize counters for the merge operation */

  _initialize_merge_counter();

#if 0 && defined(DEBUG) && !defined(NDEBUG) /* Dump local structures */
  if (param.verbosity > 3)
    _dump_vtx_eset(vtx_eset, work);
#endif

  if (param.verbosity > 1) {
    fvm_gnum_t g_n_equiv = vtx_eset->n_equiv;
    fvm_parall_counter(&g_n_equiv, 1);
    bft_printf(_("\n"
                 "  Final number of equiv. between vertices; local: %9d\n"
                 "                                          global: %9lu\n"),
               vtx_eset->n_equiv, (unsigned long)g_n_equiv);
  }

  /* Operate merge between equivalent vertices.
     Manage reduction of tolerance if necessary */

  clock_start = bft_timer_wtime();
  cpu_start = bft_timer_cpu_time();

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

    cs_int_t  *send_count = NULL, *recv_count = NULL;
    cs_int_t  *send_shift = NULL, *recv_shift = NULL;
    cs_join_vertex_t  *vtx_merge_data = NULL;

    BFT_MALLOC(send_count, n_ranks, cs_int_t);
    BFT_MALLOC(recv_count, n_ranks, cs_int_t);
    BFT_MALLOC(send_shift, n_ranks+1, cs_int_t);
    BFT_MALLOC(recv_shift, n_ranks+1, cs_int_t);

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

  clock_end = bft_timer_wtime();
  cpu_end = bft_timer_cpu_time();

  cs_join_gset_destroy(&merge_set);

  if (param.verbosity > 1)
    bft_printf(_("\n"
                 "          Vertex merge (only)\n"
                 "              wall clock time:       %10.3g\n"
                 "              cpu time:              %10.3g\n"),
               clock_end - clock_start, cpu_end - cpu_start);
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
 *   n_g_ifm_vertices     <-- initial global number of vertices (full mesh)
 *   iwm_vtx_gnum         <-- initial global vertex num. (work mesh struct)
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
                            cs_int_t                 n_iwm_vertices,
                            fvm_gnum_t               n_g_ifm_vertices,
                            const fvm_gnum_t         iwm_vtx_gnum[],
                            const fvm_gnum_t         rank_face_gnum_index[],
                            cs_join_mesh_t         **p_mesh,
                            cs_join_edges_t        **p_edges,
                            cs_join_inter_edges_t  **p_inter_edges,
                            cs_join_mesh_t         **p_local_mesh,
                            fvm_gnum_t              *p_o2n_vtx_gnum[])
{
  cs_int_t  n_af_vertices = 0; /* new number of vertices after merge */
  cs_int_t  *o2n_vtx_id = NULL;
  fvm_gnum_t  *o2n_vtx_gnum = NULL;
  cs_join_mesh_t  *mesh = *p_mesh;
  cs_join_mesh_t  *local_mesh = *p_local_mesh;
  cs_join_edges_t  *edges = *p_edges;
  cs_join_inter_edges_t  *inter_edges = *p_inter_edges;

  /* Keep an history of the evolution of each vertex */

  _keep_global_vtx_evolution(n_iwm_vertices,   /* n_vertices before inter */
                             n_g_ifm_vertices,
                             iwm_vtx_gnum,
                             mesh->n_vertices, /* n_vertices after inter */
                             mesh->vertices,
                             &o2n_vtx_gnum);   /* defined by block in // */

  _keep_local_vtx_evolution(mesh->n_vertices, /* n_vertices after inter */
                            mesh->vertices,
                            &n_af_vertices,   /* n_vertices after merge */
                            &o2n_vtx_id);

  /* Update all structures which keeps data about vertices */

  if (inter_edges != NULL) { /* The join type is not conform */

    /* Update inter_edges structure */

    _update_inter_edges_after_merge(o2n_vtx_id, &inter_edges);

    assert(edges->n_edges == inter_edges->n_edges);  /* Else: problem for
                                                        future synchro. */

#if 0 && defined(DEBUG) && !defined(NDEBUG) /* Dump local structures */
    cs_join_inter_edges_dump(inter_edges, edges, mesh);
#endif

    /* Update cs_join_mesh_t structure after the merge of vertices
       numbering of the old vertices + add new vertices */

    cs_join_mesh_update(mesh,
                        edges,
                        inter_edges->index,
                        inter_edges->vtx_lst,
                        n_af_vertices,
                        o2n_vtx_id);

  } /* End if inter_edges != NULL */

  else
    /* Update cs_join_mesh_t structure after the merge of vertices
       numbering of the old vertices + add new vertices */

    cs_join_mesh_update(mesh,
                        edges,
                        NULL,
                        NULL,
                        n_af_vertices,
                        o2n_vtx_id);

  BFT_FREE(o2n_vtx_id);

  /* Update local_mesh by redistributing mesh */

  _redistribute_mesh(rank_face_gnum_index,
                     mesh,
                     &local_mesh);

  /* Clean mesh: remove degenerate and empty edges */

  cs_join_mesh_clean(mesh, param.verbosity);

  /* Define a new cs_join_edges_t structure */

  cs_join_mesh_destroy_edges(&edges);
  edges = cs_join_mesh_define_edges(mesh);

  /* Set return pointers */

  *p_mesh = mesh;
  *p_edges = edges;
  *p_inter_edges = inter_edges;
  *p_o2n_vtx_gnum = o2n_vtx_gnum;
  *p_local_mesh = local_mesh;
}

/*---------------------------------------------------------------------------*/

END_C_DECLS
