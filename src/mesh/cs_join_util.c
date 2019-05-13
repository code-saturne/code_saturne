/*============================================================================
 * Manipulation of low-level structures for the joining operations
 *===========================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *---------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_defs.h"
#include "fvm_io_num.h"

#include "cs_join_util.h"
#include "cs_file.h"
#include "cs_mesh.h"
#include "cs_order.h"
#include "cs_search.h"
#include "cs_sort.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *---------------------------------------------------------------------------*/

#include "cs_join_set.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *===========================================================================*/

int  cs_glob_join_count = 0;
int  cs_glob_n_joinings = 0;
cs_join_t  **cs_glob_join_array = NULL;

FILE  *cs_glob_join_log = NULL;

/*============================================================================
 * Macro and type definitions
 *===========================================================================*/

/* Directory name separator
   (historically, '/' for Unix/Linux, '\' for Windows, ':' for Mac
   but '/' should work for all on modern systems) */

#define DIR_SEPARATOR '/'

/*============================================================================
 * Private function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Initialize a cs_join_param_t structure.
 *
 * parameters:
 *   join_num      <-- number of the current joining operation
 *   fraction      <-- value of the fraction parameter
 *   plane         <-- value of the plane parameter
 *   perio_type    <-- periodicity type (FVM_PERIODICITY_NULL if none)
 *   perio_matrix  <-- periodicity transformation matrix
 *   verbosity     <-- level of verbosity required
 *   visualization <-- level of visualization required
 *   preprocessing <-- is joining part of the preprocessing stage ?
 *
 * returns:
 *   a pointer to a cs_join_param_t structure
 *---------------------------------------------------------------------------*/

static cs_join_param_t
_join_param_define(int                      join_num,
                   float                    fraction,
                   float                    plane,
                   fvm_periodicity_type_t   perio_type,
                   double                   perio_matrix[3][4],
                   int                      verbosity,
                   int                      visualization,
                   bool                     preprocessing)
{
  double  cplane;
  cs_join_param_t  param;

  param.num = join_num;

  /* Possible periodicity information */

  param.perio_type = perio_type;

  if (param.perio_type == FVM_PERIODICITY_NULL) {
    int i, j;
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 4; j++)
        param.perio_matrix[i][j] = 0.;
    }
  }
  else
    memcpy(param.perio_matrix, perio_matrix, sizeof(double)*12);

  /* geometric parameters */

  /* parameter used to compute the tolerance associated to each vertex.
     Also used for finding equivalent vertices during edge intersections */

  param.fraction = fraction;

  /* parameter used to judge if two faces are in the same plane (during
     the face splitting) */

  param.plane = plane;
  cplane = cos(param.plane *acos(-1.0)/180.);
  param.plane_criteria = cplane * cplane;

  /* Coef. used to modify the tolerance associated to each vertex BEFORE the
     merge operation.
     If coef = 0.0 => no vertex merge
     If coef < 1.0 => reduce vertex merge
     If coef = 1.0 => no change
     If coef > 1.0 => increase vertex merge */

  param.merge_tol_coef = 1.0;

  /* Coef. used to compute a limit under which two vertices are merged
     before the merge step.
     Default value: 1e-3; Should be small. [1e-4, 1e-2] */

   param.pre_merge_factor = 0.05;

  /* Maximum number of equivalence breaks */

   param.n_max_equiv_breaks = 500;

   /* Tolerance computation mode: tcm
      1: (default) tol = min. edge length related to a vertex * fraction
      2: tolerance is computed like in mode 1 with in addition, the
         multiplication by a coef. which is equal to the max sin(e1, e2)
         where e1 and e2 are two edges sharing the same vertex V for which
         we want to compute the tolerance
     11: like 1 but only in taking into account only the selected faces
     12: like 2 but only in taking into account only the selected faces
   */

   param.tcm = 1;

   /* Intersection computation mode: icm
      1: (default) Original algorithm. Try to clip intersection on extremity
      2: New intersection algorithm. Avoid to clip intersection on extremity
   */

   param.icm = 1;

   /* Maximum number of sub-faces for an initial selected face
      Default value: 200 */

   param.max_sub_faces = 200;

   /* Deepest level reachable during tree building
      Default value: 30  */

   param.tree_max_level = 30;

   /* Max. number of boxes which can be related to a leaf of the tree
      if level != tree_max_level
      Default value: 25  */

   param.tree_n_max_boxes = 25;

   /* Stop tree building if: n_linked_boxes > tree_max_box_ratio*n_init_boxes
      Default value: 5.0 */

   param.tree_max_box_ratio = 5.0;
   param.tree_max_box_ratio_distrib = 2.0;

   /* Level of display */

   param.verbosity = verbosity;
   param.visualization = visualization;

   /* Is joining part of preprocessing ? */

   param.preprocessing = preprocessing;

   return param;
}

/*----------------------------------------------------------------------------
 * Initialize a cs_join_stats_t structure.
 *
 * parameters:
 *   join_num      <-- number of the current joining operation
 *   fraction      <-- value of the fraction parameter
 *   plane         <-- value of the plane parameter
 *   perio_type    <-- periodicity type (FVM_PERIODICITY_NULL if none)
 *   perio_matrix  <-- periodicity transformation matrix
 *   verbosity     <-- level of verbosity required
 *   visualization <-- level of visualization required
 *   preprocessing <-- is joining part of the preprocessing stage ?
 *
 * returns:
 *   a pointer to a cs_join_param_t structure
 *---------------------------------------------------------------------------*/

static cs_join_stats_t
_join_stats_init(void)
{
  cs_join_stats_t  stats;

  memset(&stats, 0, sizeof(cs_join_stats_t));

  CS_TIMER_COUNTER_INIT(stats.t_box_build);
  CS_TIMER_COUNTER_INIT(stats.t_box_query);
  CS_TIMER_COUNTER_INIT(stats.t_inter_sort);

  CS_TIMER_COUNTER_INIT(stats.t_l_join_mesh);
  CS_TIMER_COUNTER_INIT(stats.t_edge_inter);
  CS_TIMER_COUNTER_INIT(stats.t_new_vtx);
  CS_TIMER_COUNTER_INIT(stats.t_merge_vtx);
  CS_TIMER_COUNTER_INIT(stats.t_u_merge_vtx);
  CS_TIMER_COUNTER_INIT(stats.t_split_faces);

  CS_TIMER_COUNTER_INIT(stats.t_total);

  return stats;
}

/*----------------------------------------------------------------------------
 * Initialize a structure for the synchronization of single
 * elements
 *
 * returns:
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

  if (_sync == NULL)
    return;

  if (_sync->array != NULL)
    BFT_FREE(_sync->array);
  if (_sync->ranks != NULL)
    BFT_FREE(_sync->ranks);
  BFT_FREE(_sync->index);

  BFT_FREE(_sync);

  *sync = _sync;
}

/*----------------------------------------------------------------------------
 * Reduce numbering for the selected boundary faces.
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
_compact_face_gnum_selection(cs_lnum_t   n_select_faces,
                             cs_gnum_t  *reduce_gnum[],
                             cs_gnum_t  *reduce_gnum_index[])
{
  cs_lnum_t  i;

  cs_gnum_t  shift = 0;
  cs_gnum_t  *_reduce_gnum = *reduce_gnum;
  cs_gnum_t  *_reduce_gnum_index = *reduce_gnum_index;

  const int  n_ranks = cs_glob_n_ranks;
  const int  local_rank = CS_MAX(cs_glob_rank_id, 0);

  assert(_reduce_gnum_index == NULL);

  BFT_MALLOC(_reduce_gnum_index, n_ranks + 1, cs_gnum_t);

  for (i = 0; i < n_ranks; i++)
    _reduce_gnum_index[i] = 0;

  if (n_ranks > 1) {
#if defined(HAVE_MPI)
    MPI_Comm  mpi_comm = cs_glob_mpi_comm;
    cs_gnum_t  _n_faces = n_select_faces;

    MPI_Allgather(&_n_faces, 1, CS_MPI_GNUM,
                  &(_reduce_gnum_index[1]), 1, CS_MPI_GNUM, mpi_comm);
#endif

    for (i = 0; i < n_ranks; i++)
      _reduce_gnum_index[i+1] += _reduce_gnum_index[i];

    shift = _reduce_gnum_index[local_rank];

  }
  else {

    assert(n_ranks == 1);
    _reduce_gnum_index[n_ranks] = (cs_gnum_t)n_select_faces;

  }

  BFT_MALLOC(_reduce_gnum, n_select_faces, cs_gnum_t);

  for (i = 0; i < n_select_faces; i++)
    _reduce_gnum[i] = shift + i + 1;

  /* Returns pointer */

  *reduce_gnum = _reduce_gnum;
  *reduce_gnum_index = _reduce_gnum_index;
}

/*----------------------------------------------------------------------------
 * Extract faces implied in the current joining operation.
 * These are faces which share at least one vertex which is in the
 * select_vertices array.
 *
 * parameters:
 *   n_vertices        <-- number of vertices in the whole mesh
 *   selection         <-- pointer to a selection structure
 *   n_faces           <-- number of faces in the whole mesh
 *   f2v_idx           <-- "face -> vertex" connect. index
 *   f2v_lst           <-- "face -> vertex" connect. list
 *   n_contig_faces    <-> pointer to the number of contiguous faces
 *   contig_faces      <-> pointer to the list of contiguous faces
 *---------------------------------------------------------------------------*/

static void
_extract_contig_faces(cs_lnum_t          n_vertices,
                      cs_join_select_t  *selection,
                      cs_lnum_t          n_faces,
                      const cs_lnum_t    f2v_idx[],
                      const cs_lnum_t    f2v_lst[],
                      cs_lnum_t         *n_contig_faces,
                      cs_lnum_t         *contig_faces[])
{
  cs_lnum_t  i, j,  vtx_id, shift;

  cs_lnum_t   _n_contig_faces = 0;
  cs_lnum_t  *_contig_faces = NULL, *counter = NULL;
  cs_lnum_t  *v2f_idx = NULL, *v2f_lst = NULL;

  const cs_lnum_t  n_select_vertices = selection->n_vertices;
  const cs_lnum_t  n_single_vertices = selection->s_vertices->n_elts;
  const cs_lnum_t  *select_vertices = selection->vertices;
  const cs_lnum_t  *single_vertices = selection->s_vertices->array;

  if (n_select_vertices + n_single_vertices == 0)
    return;

  /* Reverse face -> vertex connectivity */

  BFT_MALLOC(counter, n_vertices, cs_lnum_t);

  for (i = 0; i < n_vertices; i++)
    counter[i] = 0;

  for (i = 0; i < n_faces; i++) {
    for (j = f2v_idx[i]; j < f2v_idx[i+1]; j++) {
      vtx_id = f2v_lst[j];
      counter[vtx_id] += 1;
    }
  } /* End of loop on faces */

  /* Define v2f_idx */

  BFT_MALLOC(v2f_idx, n_vertices + 1, cs_lnum_t);

  v2f_idx[0] = 0;
  for (i = 0; i < n_vertices; i++)
    v2f_idx[i+1] = v2f_idx[i] + counter[i];

  for (i = 0; i < n_vertices; i++)
    counter[i] = 0;

  /* Define v2f_lst */

  BFT_MALLOC(v2f_lst, v2f_idx[n_vertices], cs_lnum_t);

  for (i = 0; i < n_faces; i++) {

    for (j = f2v_idx[i]; j < f2v_idx[i+1]; j++) {

      vtx_id = f2v_lst[j];
      shift = v2f_idx[vtx_id] + counter[vtx_id];
      v2f_lst[shift] = i+1;
      counter[vtx_id] += 1;

    }

  } /* End of loop on faces */

  BFT_REALLOC(counter, n_faces, cs_lnum_t);

  for (i = 0; i < n_faces; i++)
    counter[i] = 0;

  /* Count the number of contiguous faces */

  for (i = 0; i < n_select_vertices; i++) {

    vtx_id = select_vertices[i] - 1;

    for (j = v2f_idx[vtx_id]; j < v2f_idx[vtx_id+1]; j++)
      counter[v2f_lst[j]-1] = 1;

  }

  for (i = 0; i < n_single_vertices; i++) {

    vtx_id = single_vertices[i] - 1;

    for (j = v2f_idx[vtx_id]; j < v2f_idx[vtx_id+1]; j++)
      counter[v2f_lst[j]-1] = 1;

  }

  for (i = 0; i < n_faces; i++)
    _n_contig_faces += counter[i];

  /* Define contig_faces */

  BFT_MALLOC(_contig_faces, _n_contig_faces, cs_lnum_t);

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
 * Define a structure used for synchronizing "single" vertices.
 * Use a cs_interface_t structure to help the build.
 *
 * parameters:
 *   interfaces     <-- pointer to a cs_interface_set_t structure
 *   var_size       <-- number of elements in var buffer
 *   count          <-> counter buffer (0: not selected, 1 otherwise)
 *   related_ranks  <-> rank associated to each single vertex (size: var_size)
 *   single         <-> data about the distribution of single vertices
 *   verbosity      <-- verbosity level
 *----------------------------------------------------------------------------*/

static void
_add_single_vertices(cs_interface_set_t  *interfaces,
                     cs_lnum_t            var_size,
                     cs_lnum_t           *count,
                     cs_lnum_t           *related_ranks,
                     cs_join_sync_t      *single,
                     int                  verbosity)
{
  int  distant_rank, n_interfaces, last_found_rank;
  cs_lnum_t  id, ii;

  int  n_max_ranks = 0;
  cs_lnum_t  count_size = 0, n_entities = 0, total_size = 0, n_max_elts = 0;
  cs_lnum_t  *buf = NULL, *recv_buf = NULL;

  const cs_lnum_t  *local_id = NULL;
  const cs_interface_t  *interface = NULL;

  assert(count != NULL);
  assert(single != NULL);

  /* Initialize and allocate */

  n_interfaces = cs_interface_set_size(interfaces);

  total_size = cs_interface_set_n_elts(interfaces);

  BFT_MALLOC(buf, total_size, cs_lnum_t);

  /* Exchange with distant ranks */

  cs_interface_set_copy_array(interfaces,
                              CS_LNUM_TYPE,
                              1,
                              true,
                              count,
                              buf);

  /* Now we estimate the max. number of ranks and elements involved */

  count_size = 0;
  last_found_rank = -1;

  for (id = 0; id < n_interfaces; id++) {

    /* Scan data */

    interface = cs_interface_set_get(interfaces, id);
    distant_rank = cs_interface_rank(interface);
    n_entities = cs_interface_size(interface);
    local_id = cs_interface_get_elt_ids(interface);

    recv_buf = buf + count_size;

    for (ii = 0; ii < n_entities; ii++) {

      int  vtx_id = local_id[ii];

      assert(vtx_id < var_size);

      if (count[vtx_id] == 0 && recv_buf[ii] > 0) { /* Find a single vertex */

        if (last_found_rank != distant_rank) {
          last_found_rank = distant_rank;
          n_max_ranks++;
        }
        n_max_elts++;

      }

    }

    count_size += n_entities;

  }

  if (n_max_elts > 0) { /* We have found single vertices */

    BFT_MALLOC(single->ranks, n_max_ranks, int);
    BFT_MALLOC(single->index, n_max_ranks + 1, cs_lnum_t);
    BFT_MALLOC(single->array, n_max_elts, cs_lnum_t);

    count_size = 0;
    last_found_rank = -1;
    single->index[0] = 0;
    single->n_elts = 0;
    single->n_ranks = 0;

    for (id = 0; id < n_interfaces; id++) {

      /* Scan data */

      interface = cs_interface_set_get(interfaces, id);
      distant_rank = cs_interface_rank(interface);
      n_entities = cs_interface_size(interface);
      local_id = cs_interface_get_elt_ids(interface);

      recv_buf = buf + count_size;

      for (ii = 0; ii < n_entities; ii++) {

        int  vtx_id = local_id[ii];

        if (count[vtx_id] == 0 && recv_buf[ii] > 0) {

          if (last_found_rank != distant_rank) {
            last_found_rank = distant_rank;
            single->ranks[single->n_ranks++] = distant_rank;
          }

          single->array[single->n_elts++] = local_id[ii] + 1;
          single->index[single->n_ranks] = single->n_elts;

          related_ranks[vtx_id] = distant_rank;
          count[vtx_id] = recv_buf[ii];

        }

      }

      count_size += n_entities;

    } /* End of loop on interfaces */

    BFT_REALLOC(single->ranks, single->n_ranks, int);
    BFT_REALLOC(single->index, single->n_ranks + 1, cs_lnum_t);
    BFT_REALLOC(single->array, single->n_elts, cs_lnum_t);

  } /* End if n_max_elts > 0 */

  BFT_FREE(buf);

  if (cs_glob_join_log != NULL && verbosity > 3) {
    int  i, j;
    FILE  *logfile = cs_glob_join_log;

    fprintf(logfile,
            "\n  Single vertices for the joining operation: (%p)\n",
            (void *)single);
    fprintf(logfile,
            "  Single vertices: n_elts : %8d\n"
            "  Single vertices: n_ranks: %8d\n",
            single->n_elts, single->n_ranks);

    if (single->n_elts > 0) {
      for (i = 0; i < single->n_ranks; i++)
        for (j = single->index[i]; j < single->index[i+1]; j++)
          fprintf(logfile, " %9d | %6d | %9d\n",
                     j, single->ranks[i], single->array[j]);
      fprintf(logfile, "\n");
    }
    fflush(logfile);
  }

}

/*-----------------------------------------------------------------------------
 * Define a structure used for synchronizing "single" vertices.
 * Use a cs_interface_t structure to help the build.
 *
 * parameters:
 *   interfaces     --> pointer to a cs_interface_set_t structure
 *   var_size       --> number of elements in var buffer
 *   related_ranks  <-> rank buffer for synchronization
 *   coupled        <-> pointer to a structure to build on coupled vertices
 *   verbosity      <-- verbosity level
 *----------------------------------------------------------------------------*/

static void
_add_coupled_vertices(cs_interface_set_t  *interfaces,
                      cs_lnum_t            var_size,
                      int                 *related_ranks,
                      cs_join_sync_t      *coupled,
                      int                  verbosity)
{
  int  distant_rank, n_interfaces, last_found_rank;
  cs_lnum_t id, ii;

  cs_lnum_t  n_entities = 0, count_size = 0, total_size = 0;
  int  *buf = NULL, *recv_buf = NULL;

  cs_datatype_t int_type = (sizeof(int) == 8) ? CS_INT64 : CS_INT32;

  const int  local_rank = CS_MAX(cs_glob_rank_id, 0);
  const cs_lnum_t  *local_id = NULL;
  const cs_interface_t  *interface = NULL;

  assert(related_ranks != NULL);
  assert(coupled != NULL);

  /* Initialize and allocate */

  n_interfaces = cs_interface_set_size(interfaces);

  n_interfaces = cs_interface_set_size(interfaces);

  total_size = cs_interface_set_n_elts(interfaces);

  BFT_MALLOC(buf, total_size, int);

  /* Exchange with distant ranks */

  cs_interface_set_copy_array(interfaces,
                              int_type,
                              1,
                              true,
                              related_ranks,
                              buf);

  /* Now we scan each part to build coupled */

  count_size = 0;
  coupled->n_elts = 0;
  coupled->n_ranks = 0;
  last_found_rank = -1;

  for (id = 0; id < n_interfaces; id++) {

    /* Scan data */

    interface = cs_interface_set_get(interfaces, id);
    distant_rank   = cs_interface_rank(interface);
    n_entities = cs_interface_size(interface);

    recv_buf = buf + count_size;

    for (ii = 0; ii < n_entities; ii++) {

      if (recv_buf[ii] == local_rank) { /* Coupled vertex found */
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

    BFT_MALLOC(coupled->array, coupled->n_elts, cs_lnum_t);
    BFT_MALLOC(coupled->index, coupled->n_ranks + 1, cs_lnum_t);
    BFT_MALLOC(coupled->ranks, coupled->n_ranks, int);

    coupled->index[0] = 0;

    count_size = 0;
    last_found_rank = -1;

    for (id = 0; id < n_interfaces; id++) {

      /* Retrieve data */

      interface = cs_interface_set_get(interfaces, id);
      distant_rank   = cs_interface_rank(interface);
      n_entities = cs_interface_size(interface);
      local_id = cs_interface_get_elt_ids(interface);

      recv_buf = buf + count_size;

      for (ii = 0; ii < n_entities; ii++) {

        if (recv_buf[ii] == local_rank) {

          if (last_found_rank != distant_rank) {
            last_found_rank = distant_rank;
            coupled->ranks[rank_shift++] = distant_rank;
          }

          assert(local_id[ii] < var_size);
          coupled->array[vtx_shift++] = local_id[ii] + 1;
          coupled->index[rank_shift] = vtx_shift;
        }

      }
      count_size += n_entities;

    }

  } /* End if coupled->n_elts > 0 */

  BFT_FREE(buf);

  if (cs_glob_join_log != NULL && verbosity > 3) {

    int  i, j;
    FILE  *logfile = cs_glob_join_log;

    fprintf(logfile, "\n  Coupled vertices for the joining operation: (%p)\n",
            (void *)coupled);
    fprintf(logfile,
            "  Coupled vertices: n_elts : %8d\n"
            "  Coupled vertices: n_ranks: %8d\n",
            coupled->n_elts, coupled->n_ranks);

    if (coupled->n_elts > 0) {
      for (i = 0; i < coupled->n_ranks; i++)
        for (j = coupled->index[i]; j < coupled->index[i+1]; j++)
          fprintf(logfile, " %9d | %6d | %9d\n",
                  j, coupled->ranks[i], coupled->array[j]);
      fprintf(logfile, "\n");
   }
    fflush(logfile);

  }

}

/*----------------------------------------------------------------------------
 * Get the full selection of vertices to extract. Only in parallel. Somme
 * vertices may have been selected on another rank and not on the local rank
 * but you have to take them into account to have a good update of the mesh.
 *
 * parameters:
 *   n_vertices   <--  number of vertices in the parent mesh
 *   ifs          <--  pointer to a cs_interface_set_t struct.
 *   p_vtx_tag    <->  pointer to an array on vertices. tag=1 if selected
 *   join_select  <->  pointer to a fvm_join_selection_t structure
 *   verbosity    <--  verbosity level
 *---------------------------------------------------------------------------*/

static void
_get_missing_vertices(cs_lnum_t            n_vertices,
                      cs_interface_set_t  *ifs,
                      cs_lnum_t           *p_vtx_tag[],
                      cs_join_select_t    *selection,
                      int                  verbosity)
{
  cs_lnum_t  i;
  cs_gnum_t  n_l_elts, n_g_elts;

  cs_lnum_t  *vtx_tag = NULL, *related_ranks = NULL;

  /* Define a counter on vertices. 1 if selected, 0 otherwise */

  BFT_MALLOC(vtx_tag, n_vertices, cs_lnum_t);
  BFT_MALLOC(related_ranks, n_vertices, cs_lnum_t);

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
                       selection->s_vertices,
                       verbosity);

  n_l_elts = selection->s_vertices->n_elts;

  MPI_Allreduce(&n_l_elts, &n_g_elts, 1, CS_MPI_GNUM,
                MPI_SUM, cs_glob_mpi_comm);

  if (n_g_elts > 0) {

    bft_printf(_("\n  Global number of single vertices found: %6llu\n"),
               (unsigned long long)n_g_elts);
    bft_printf_flush();
    selection->do_single_sync = true;

    _add_coupled_vertices(ifs,
                          n_vertices,
                          related_ranks,
                          selection->c_vertices,
                          verbosity);

  } /* End if n_g_elts > 0 */

  /* Return pointers */

  *p_vtx_tag = vtx_tag;

  /* Free memory */

  BFT_FREE(related_ranks);
}

/*----------------------------------------------------------------------------
 * Define a vertex -> vertex connectivity for vertices belonging to the
 * selected boundary faces.
 *
 * parameters:
 *   n_vertices  <--  number of vertices in the parent mesh
 *   selection   <--  pointer to a fvm_join_selection_t structure
 *   b_f2v_idx   <--  boundary "face -> vertex" connect. index
 *   b_f2v_lst   <--  boundary "face -> vertex" connect. list
 *   p_v2v_idx   <->  vertex -> vertex connect. index
 *   p_v2v_lst   <->  vertex -> vertex connect. list
 *---------------------------------------------------------------------------*/

static void
_get_select_v2v_connect(cs_lnum_t              n_vertices,
                        cs_join_select_t      *selection,
                        cs_lnum_t              b_f2v_idx[],
                        cs_lnum_t              b_f2v_lst[],
                        cs_lnum_t             *p_v2v_idx[],
                        cs_lnum_t             *p_v2v_lst[])
{
  cs_lnum_t  i, j, save, s, e, n_sel_edges, shift;

  cs_lnum_t  *count = NULL, *sel_v2v_idx = NULL, *sel_v2v_lst = NULL;

  /* Build a vertex -> vertex connectivity for the selected boundary faces  */

  BFT_MALLOC(sel_v2v_idx, n_vertices + 1, cs_lnum_t);

  for (i = 0; i < n_vertices + 1; i++)
    sel_v2v_idx[i] = 0;

  cs_join_build_edges_idx(selection->n_faces,
                          selection->faces,
                          b_f2v_idx,
                          b_f2v_lst,
                          sel_v2v_idx);

  BFT_MALLOC(count, n_vertices, cs_lnum_t);

  for (i = 0; i < n_vertices; i++) {
    sel_v2v_idx[i+1] += sel_v2v_idx[i];
    count[i] = 0;
  }

  BFT_MALLOC(sel_v2v_lst, sel_v2v_idx[n_vertices], cs_lnum_t);

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
  BFT_REALLOC(sel_v2v_lst, n_sel_edges, cs_lnum_t);

  /* Return pointers */

  *p_v2v_idx = sel_v2v_idx;
  *p_v2v_lst = sel_v2v_lst;

}

/*----------------------------------------------------------------------------
 * Get the related edge id from a couple of vertex ids.
 *
 * parameters:
 *   v1_id     <-- first vertex id
 *   v2_id     <-- second vertex id
 *   v2v_idx   <-- vertex -> vertex connect. index
 *   v2v_lst   <-- vertex -> vertex connect. list
 *
 * returns:
 *   related edge_id in cs_join_edges_t structure
 *---------------------------------------------------------------------------*/

inline static cs_lnum_t
_get_edge_id(cs_lnum_t        v1_id,
             cs_lnum_t        v2_id,
             const cs_lnum_t  v2v_idx[],
             const cs_lnum_t  v2v_lst[])
{
  int  i, va, vb;

  cs_lnum_t  edge_id = -1;

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
 *   vertex_tag   <--  tag to know if a vertices is in selection
 *   v1_id        <--  first vertex id
 *   v2_id        <--  second vertex id
 *   ref_v2v_idx  <--  vertex -> vertex connect. index
 *   ref_v2v_lst  <--  vertex -> vertex connect. list
 *   p_tmp_size   <->  pointer to the current number of single edges
 *   p_max_size   <->  pointer to the max allocated size of new_s_vertices
 *   p_tmp_size   <->  pointer to the single edges definition
 *---------------------------------------------------------------------------*/

static void
_add_s_edge(cs_lnum_t        vertex_tag[],
            cs_lnum_t        v1_id,
            cs_lnum_t        v2_id,
            const cs_lnum_t  sel_v2v_idx[],
            const cs_lnum_t  sel_v2v_lst[],
            cs_lnum_t       *p_tmp_size,
            cs_lnum_t       *p_max_size,
            cs_lnum_t       *p_tmp_edges[])
{
  cs_lnum_t  i, a, b, edge_id;

  if (vertex_tag[v1_id] > 0 && vertex_tag[v2_id] > 0) {

    _Bool  is_found = false;

    cs_lnum_t  tmp_size = *p_tmp_size;
    cs_lnum_t  max_size = *p_max_size;
    cs_lnum_t *tmp_edges = *p_tmp_edges;

    edge_id = _get_edge_id(v1_id, v2_id, sel_v2v_idx, sel_v2v_lst);

    if (edge_id == -1) { /* Edge not found among the selected edges */

      assert(v1_id != v2_id);
      if (v1_id < v2_id)
        a = v1_id, b = v2_id;
      else
        a = v2_id, b = v1_id;

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
          BFT_REALLOC(tmp_edges, 2*max_size, cs_lnum_t);
        }

      }

    } /* this edge should be among the selected edges */

    /* Return pointers */

    *p_max_size = max_size;
    *p_tmp_size = tmp_size;
    *p_tmp_edges = tmp_edges;

  }
}

/*----------------------------------------------------------------------------
 * Copy selected edge definitions at parallel interfaces (only in parallel).
 *
 * The returned arrays (sel_v2v_recv_idx and sel_v2v_recv_lst) both use
 * 0 to n-1 numbering.

 * The caller is responsible for freeing those returned arrays.
 *
 * parameters:
 *   ifs              <-- pointer to the interface set on vertices
 *   sel_v2v_idx      <-- vertex -> vertex connect. index
 *   sel_v2v_lst      <-- vertex -> vertex connect. list
 *   sel_v2v_recv_idx --> pointer to received edges index
 *   sel_v2v_recv_lst --> pointer to received edges connectivity
 *---------------------------------------------------------------------------*/

static void
_copy_interface_edges(cs_interface_set_t   *ifs,
                      cs_lnum_t             sel_v2v_idx[],
                      cs_lnum_t             sel_v2v_lst[],
                      cs_lnum_t           **sel_v2v_recv_idx,
                      cs_lnum_t           **sel_v2v_recv_lst)
{
  int i;
  cs_lnum_t  j, if_shift;

  /* Exchange v2v info;
     only info relative to edges where both vertices are available
     on the distant side is sent. */

  int n_interfaces = cs_interface_set_size(ifs);

  cs_lnum_t ifs_tot_size = cs_interface_set_n_elts(ifs);

  cs_lnum_t *_sel_v2v_send_idx = NULL, *_sel_v2v_recv_idx = NULL;
  cs_lnum_t *_sel_v2v_send_lst = NULL, *_sel_v2v_recv_lst = NULL;

  BFT_MALLOC(_sel_v2v_send_idx, ifs_tot_size + 1, cs_lnum_t);
  BFT_MALLOC(_sel_v2v_recv_idx, ifs_tot_size + 1, cs_lnum_t);

  /* Counting pass for v2v_info exchange
     (send/receive indexes are prepared as counts) */

  if_shift = 0;
  _sel_v2v_send_idx[0] = 0;
  _sel_v2v_recv_idx[0] = 0;

  for (i = 0; i < n_interfaces; i++) {

    const cs_interface_t *interface = cs_interface_set_get(ifs, i);
    const cs_lnum_t *local_id = cs_interface_get_elt_ids(interface);
    const cs_lnum_t if_size = cs_interface_size(interface);

    for (j = 0; j < if_size; j++) {

      cs_lnum_t k = local_id[j];
      cs_lnum_t s_id = sel_v2v_idx[k];
      cs_lnum_t e_id = sel_v2v_idx[k+1];

      _sel_v2v_send_idx[if_shift + j + 1] = 0;

      for (cs_lnum_t l = s_id; l < e_id; l++) {
        if (cs_search_binary(if_size, sel_v2v_lst[l] - 1, local_id) > -1)
          _sel_v2v_send_idx[if_shift + j + 1] += 1;
      }

    }

    if_shift += if_size;

  }

  cs_interface_set_copy_array(ifs,
                              CS_LNUM_TYPE,
                              1,
                              false,
                              _sel_v2v_send_idx + 1,
                              _sel_v2v_recv_idx + 1);

  for (j = 0; j < ifs_tot_size; j++) {
    _sel_v2v_send_idx[j+1] += _sel_v2v_send_idx[j];
    _sel_v2v_recv_idx[j+1] += _sel_v2v_recv_idx[j];
  }

  /* Data pass for v2v_info exchange */

  cs_interface_set_add_match_ids(ifs);

  BFT_MALLOC(_sel_v2v_send_lst, _sel_v2v_send_idx[ifs_tot_size], cs_lnum_t);
  BFT_MALLOC(_sel_v2v_recv_lst, _sel_v2v_recv_idx[ifs_tot_size], cs_lnum_t);

  if_shift = 0;

  for (i = 0; i < n_interfaces; i++) {

    const cs_interface_t *interface = cs_interface_set_get(ifs, i);
    const cs_lnum_t if_size = cs_interface_size(interface);

    const cs_lnum_t *local_id = cs_interface_get_elt_ids(interface);
    const cs_lnum_t *distant_id = cs_interface_get_match_ids(interface);

    for (j = 0; j < if_size; j++) {

      cs_lnum_t k = local_id[j];
      cs_lnum_t s_id = sel_v2v_idx[k];
      cs_lnum_t e_id = sel_v2v_idx[k+1];
      cs_lnum_t v_send_size = 0;

      for (cs_lnum_t l = s_id; l < e_id; l++) {
        cs_lnum_t m = cs_search_binary(if_size, sel_v2v_lst[l] - 1, local_id);
        if (m > -1) {
          _sel_v2v_send_lst[_sel_v2v_send_idx[if_shift + j] + v_send_size]
            = distant_id[m];
          v_send_size += 1;
        }
      }

    }

    if_shift += if_size;

  }

  cs_interface_set_copy_indexed(ifs,
                                CS_LNUM_TYPE,
                                false,
                                _sel_v2v_send_idx,
                                _sel_v2v_recv_idx,
                                _sel_v2v_send_lst,
                                _sel_v2v_recv_lst);

  /* Free temporary data and set return pointers */

  BFT_FREE(_sel_v2v_send_idx);
  BFT_FREE(_sel_v2v_send_lst);

  cs_interface_set_free_match_ids(ifs);

  *sel_v2v_recv_idx = _sel_v2v_recv_idx;
  *sel_v2v_recv_lst = _sel_v2v_recv_lst;
}

/*----------------------------------------------------------------------------
 * Get the full selection of single edges. Done only if the run is parallel.
 *
 * parameters:
 *   ifs          <-- pointer to the interface set on vertices
 *   vertex_tag   <-- tag to know if a vertices is in selection
 *   selection    <-- pointer to a fvm_join_selection_t structure
 *   sel_v2v_idx  <-- vertex -> vertex connect. index
 *   sel_v2v_lst  <-- vertex -> vertex connect. list
 *   b_f2v_idx    <-- boundary "face -> vertex" connect. index
 *   b_f2v_lst    <-- boundary "face -> vertex" connect. list
 *   i_f2v_idx    <-- interior "face -> vertex" connect. index
 *   i_f2v_lst    <-- interior "face -> vertex" connect. list
 *   i_face_cells <-- interior face -> cells connect.
 *   s_edges      <-> pointer to the single edges structure to define
 *   verbosity    <-- verbosity level
 *---------------------------------------------------------------------------*/

static void
_add_single_edges(cs_interface_set_t   *ifs,
                  cs_lnum_t             vertex_tag[],
                  cs_join_select_t     *selection,
                  cs_lnum_t             sel_v2v_idx[],
                  cs_lnum_t             sel_v2v_lst[],
                  cs_lnum_t             b_f2v_idx[],
                  cs_lnum_t             b_f2v_lst[],
                  cs_lnum_t             i_f2v_idx[],
                  cs_lnum_t             i_f2v_lst[],
                  cs_lnum_2_t           i_face_cells[],
                  cs_join_sync_t       *s_edges,
                  int                   verbosity)
{
  int i;
  cs_lnum_t  j, fid, s, e;

  int  tmp_size = 0, max_size = 10;
  cs_lnum_t *sel_v2v_recv_idx = NULL, *sel_v2v_recv_lst = NULL;
  int  *tmp_edges = NULL;

  assert(s_edges != NULL);

  /* Exchange v2v info;
     only info relative to edges where both vertices are available
     on the distant side is sent. */

  _copy_interface_edges(ifs,
                        sel_v2v_idx,
                        sel_v2v_lst,
                        &sel_v2v_recv_idx,
                        &sel_v2v_recv_lst);

  /* Scan adjacent faces to find new "single" edges */

  BFT_MALLOC(tmp_edges, 2*max_size, int);

  for (i = 0; i < selection->n_b_adj_faces; i++) {

    fid = selection->b_adj_faces[i] - 1;
    s = b_f2v_idx[fid];
    e = b_f2v_idx[fid+1];

    for (j = s; j < e - 1; j++)
      _add_s_edge(vertex_tag,
                  b_f2v_lst[j],
                  b_f2v_lst[j+1],
                  sel_v2v_idx,
                  sel_v2v_lst,
                  &tmp_size,
                  &max_size,
                  &tmp_edges);

    _add_s_edge(vertex_tag,
                b_f2v_lst[e-1],
                b_f2v_lst[s],
                sel_v2v_idx,
                sel_v2v_lst,
                &tmp_size,
                &max_size,
                &tmp_edges);

  }

  for (i = 0; i < selection->n_i_adj_faces; i++) {

    fid = selection->i_adj_faces[i] - 1;

    if (i_face_cells[fid][0] == -1 || i_face_cells[fid][1] == -1) {

      s = i_f2v_idx[fid];
      e = i_f2v_idx[fid+1];

      for (j = s; j < e - 1; j++)
        _add_s_edge(vertex_tag,
                    i_f2v_lst[j],
                    i_f2v_lst[j+1],
                    sel_v2v_idx,
                    sel_v2v_lst,
                    &tmp_size,
                    &max_size,
                    &tmp_edges);

      _add_s_edge(vertex_tag,
                  i_f2v_lst[e-1],
                  i_f2v_lst[s],
                  sel_v2v_idx,
                  sel_v2v_lst,
                  &tmp_size,
                  &max_size,
                  &tmp_edges);

    } /* Face on a parallel boundary */

  }

  /* Find the related ranks for synchronization */

  if (tmp_size > 0) {

    int   distant_rank;

    cs_lnum_t  if_shift;
    int  last_found_rank = -1;
    cs_lnum_t n_s_edges = 0;

    bool *edge_tag = NULL;

    const int  n_interfaces = cs_interface_set_size(ifs);

    s_edges->n_elts = tmp_size;
    BFT_REALLOC(tmp_edges, 2*tmp_size, int);

    BFT_MALLOC(edge_tag, tmp_size, bool);
    BFT_MALLOC(s_edges->array, 2*tmp_size, cs_lnum_t);
    BFT_MALLOC(s_edges->index, n_interfaces + 1, cs_lnum_t);
    BFT_MALLOC(s_edges->ranks, n_interfaces, int);

    for (i = 0; i < tmp_size; i++)
      edge_tag[i] = false; /* Not matched */

    /* Loop on interfaces (naturally ordered by rank) */

    if_shift = 0;

    for (i = 0; i < n_interfaces; i++) {

      const cs_interface_t *interface = cs_interface_set_get(ifs, i);
      const cs_lnum_t if_size = cs_interface_size(interface);

      const cs_lnum_t *local_id = cs_interface_get_elt_ids(interface);

      distant_rank = cs_interface_rank(interface);

      /* Loop on single edges */

      for (j = 0; j < tmp_size; j++) {

        if (edge_tag[j] == false) { /* Not already treated */

          cs_lnum_t i0 = cs_search_binary(if_size,
                                          tmp_edges[2*j],
                                          local_id);

          if (i0 > -1) {

            cs_lnum_t s_id = sel_v2v_recv_idx[if_shift + i0];
            cs_lnum_t e_id = sel_v2v_recv_idx[if_shift + i0 + 1];

            /* Search on short, unsorted array */

            for (cs_lnum_t l = s_id; l < e_id; l++) {

              if (sel_v2v_recv_lst[l] == tmp_edges[2*j+1]) {

                if (last_found_rank != distant_rank) {
                  s_edges->ranks[s_edges->n_ranks] = distant_rank;
                  s_edges->index[s_edges->n_ranks] = n_s_edges;
                  s_edges->n_ranks++;
                  last_found_rank = distant_rank;
                }

                edge_tag[j] = true; /* Tag as done */
                s_edges->array[2*n_s_edges] = tmp_edges[2*j] + 1;
                s_edges->array[2*n_s_edges+1] = tmp_edges[2*j+1] + 1;
                n_s_edges++;

                break;
              }
            }

          }

        } /* Not matched yet */

      } /* End of loop on single edges */

      if_shift += if_size;

    } /* End of loop on interfaces */

    s_edges->index[s_edges->n_ranks] = n_s_edges;
    s_edges->n_elts = n_s_edges;

    /* Memory management */

    if (s_edges->n_elts != tmp_size)
      BFT_REALLOC(s_edges->array, 2*s_edges->n_elts, cs_lnum_t);

    BFT_REALLOC(s_edges->ranks, s_edges->n_ranks, int);
    BFT_REALLOC(s_edges->index, s_edges->n_ranks + 1, cs_lnum_t);
    BFT_FREE(edge_tag);

  } /* End if tmp_size > 0 */

  /* Free memory */

  BFT_FREE(tmp_edges);
  BFT_FREE(sel_v2v_recv_idx);
  BFT_FREE(sel_v2v_recv_lst);

  /* Logging */

  if (cs_glob_join_log != NULL && verbosity > 3) {

    FILE  *logfile = cs_glob_join_log;

    fprintf(logfile, "\n  Single edges for the joining operation: (%p)\n",
            (void *)s_edges);
    fprintf(logfile,
            "  Single edges: n_elts : %8d\n"
            "  Single edges: n_ranks: %8d\n",
            s_edges->n_elts, s_edges->n_ranks);

    if (s_edges->n_elts > 0) {
      for (i = 0; i < s_edges->n_ranks; i++)
        for (j = s_edges->index[i]; j < s_edges->index[i+1]; j++)
          fprintf(logfile, " %9d | %6d | (%9d, %9d)\n",
                  j, s_edges->ranks[i], s_edges->array[2*j],
                  s_edges->array[2*j+1]);
      fprintf(logfile, "\n");
    }
    fflush(logfile);
  }
}

/*----------------------------------------------------------------------------
 * Define a structure for the coupled edges. Only done if single edges have
 * been detected and only in parallel.
 *
 * parameters:
 *  ifs          <--  pointer to the interface set on vertices
 *  s_edges      <--  single edges structure used to build coupled_edges
 *  c_edges      <->  pointer to the coupled edges structure to define
 *  verbosity    <--  verbosity level
 *---------------------------------------------------------------------------*/

static void
_add_coupled_edges(cs_interface_set_t   *ifs,
                   cs_join_sync_t       *s_edges,
                   cs_join_sync_t       *c_edges,
                   int                   verbosity)
{
  cs_lnum_t  i, j, id, n_entities;
  int  request_count, rank_shift, distant_rank;

  cs_lnum_t  *buf = NULL, *recv_buf = NULL, *send_buf = NULL;
  MPI_Request  *request = NULL;
  MPI_Status  *status  = NULL;
  MPI_Comm  mpi_comm = cs_glob_mpi_comm;

  const int  local_rank = CS_MAX(cs_glob_rank_id, 0);
  const cs_lnum_t  *local_id = NULL, *distant_id = NULL;
  const int  n_interfaces = cs_interface_set_size(ifs);
  const cs_interface_t  *interface = NULL;

  assert(s_edges != NULL);
  assert(c_edges != NULL);

  /* Exchange number of single edges */

  BFT_MALLOC(request, n_interfaces * 2, MPI_Request);
  BFT_MALLOC(status,  n_interfaces * 2, MPI_Status);
  BFT_MALLOC(buf, 2*n_interfaces, cs_lnum_t);

  for (i = 0; i < 2*n_interfaces; i++)
    buf[i] = 0;

  request_count = 0;

  for (id = 0; id < n_interfaces; id++) {

    interface = cs_interface_set_get(ifs, id);
    distant_rank = cs_interface_rank(interface);

    MPI_Irecv(&(buf[id]),
              1,
              CS_MPI_LNUM,
              distant_rank,
              distant_rank,
              mpi_comm,
              &(request[request_count++]));

  }

  /* Send */

  rank_shift = 0;

  for (id = 0; id < n_interfaces; id++) {

    /* Preparation of data to send */

    interface = cs_interface_set_get(ifs, id);
    distant_rank = cs_interface_rank(interface);

    if (rank_shift < s_edges->n_ranks) {
      if (s_edges->ranks[rank_shift] == distant_rank) {
        buf[n_interfaces + id] =
          s_edges->index[rank_shift+1] - s_edges->index[rank_shift];
        rank_shift++;
      }
    }

    MPI_Isend(&(buf[n_interfaces + id]),
              1,
              CS_MPI_LNUM,
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
  BFT_MALLOC(c_edges->index, c_edges->n_ranks + 1, cs_lnum_t);
  BFT_MALLOC(c_edges->array, 2*c_edges->n_elts, cs_lnum_t);

  for (i = 0; i < c_edges->n_ranks + 1; i++)
    c_edges->index[i] = 0;

  rank_shift = 0;

  for (i = 0; i < n_interfaces; i++) {
    if (buf[i] > 0) {

      interface = cs_interface_set_get(ifs, i);
      distant_rank = cs_interface_rank(interface);
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

    interface = cs_interface_set_get(ifs, id);
    distant_rank = cs_interface_rank(interface);

    if (rank_shift < c_edges->n_ranks) {
      if (c_edges->ranks[rank_shift] == distant_rank) {

        n_entities = c_edges->index[rank_shift+1] - c_edges->index[rank_shift];
        n_entities *= 2; /* couple of vertices */
        recv_buf = c_edges->array + 2*c_edges->index[rank_shift];
        rank_shift++;

        MPI_Irecv(recv_buf, /* receive distant num */
                  n_entities,
                  CS_MPI_LNUM,
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

    interface = cs_interface_set_get(ifs, id);
    distant_rank = cs_interface_rank(interface);

    if (rank_shift < s_edges->n_ranks) {
      if (s_edges->ranks[rank_shift] == distant_rank) {

        n_entities = s_edges->index[rank_shift+1] - s_edges->index[rank_shift];
        n_entities *= 2; /* couple of vertices */
        send_buf = s_edges->array + 2*s_edges->index[rank_shift];
        rank_shift++;

        MPI_Isend(send_buf,
                  n_entities,
                  CS_MPI_LNUM,
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

  /* Switch received couple vertices from distant id. to local id. */

  cs_interface_set_add_match_ids(ifs);

  for (i = 0; i < n_interfaces; i++) {

    interface = cs_interface_set_get(ifs, i);
    distant_rank = cs_interface_rank(interface);
    distant_id = cs_interface_get_match_ids(interface);

    n_entities = cs_interface_size(interface);

    if (rank_shift < c_edges->n_ranks) {
      if (c_edges->ranks[rank_shift] == distant_rank) {

        local_id = cs_interface_get_elt_ids(interface);

        for (j = c_edges->index[rank_shift];
             j < c_edges->index[rank_shift+1];
             j++) {

          id = cs_search_binary(n_entities, c_edges->array[2*j] - 1, distant_id);
          assert(id != -1);
          c_edges->array[2*j] = local_id[id] + 1;

          id = cs_search_binary(n_entities,
                                c_edges->array[2*j+1] - 1,
                                distant_id);
          assert(id != -1);
          c_edges->array[2*j+1] = local_id[id] + 1;

        } /* Loop on couple edges */

        rank_shift++;

      } /* This rank is in the list of related ranks */
    }

  } /* Loop on interfaces */

  cs_interface_set_free_match_ids(ifs);

  if (cs_glob_join_log != NULL && verbosity > 3) {

    FILE  *logfile = cs_glob_join_log;

    fprintf(logfile, "\n  Coupled edges for the joining operation: (%p)\n",
            (void *)c_edges);
    fprintf(logfile,
            "  Coupled edges: n_elts : %8d\n"
            "  Coupled edges: n_ranks: %8d\n",
            c_edges->n_elts, c_edges->n_ranks);

    if (c_edges->n_elts > 0) {
      int  shift;
      for (i = 0, shift = 0; i < c_edges->n_ranks; i++) {
        for (j = c_edges->index[i]; j < c_edges->index[i+1]; j++) {
          fprintf(logfile, " %9d | %6d | (%9d, %9d)\n",
                  shift, c_edges->ranks[i],
                  c_edges->array[2*j], c_edges->array[2*j+1]);
          shift++;
        }
      }
      fprintf(logfile, "\n");
      fflush(logfile);
    }
  }

}

/*----------------------------------------------------------------------------
 * Get the full selection of single edges. Done only if the run is parallel.
 *
 * parameters:
 *  selection    <-> pointer to a fvm_join_selection_t structure
 *  sel_v2v_idx  <-- vertex -> vertex connect. index
 *  sel_v2v_lst  <-- vertex -> vertex connect. list
 *---------------------------------------------------------------------------*/

static void
_filter_edge_element(cs_join_select_t   *selection,
                     const cs_lnum_t     sel_v2v_idx[],
                     const cs_lnum_t     sel_v2v_lst[])
{
  cs_lnum_t  i, j, vid1, vid2, edge_id, request_count, shift, save;

  int  *c_edge_tag = NULL, *s_edge_tag = NULL;
  cs_join_sync_t  *s_edges = selection->s_edges;
  cs_join_sync_t  *c_edges = selection->c_edges;

  MPI_Request  *request = NULL;
  MPI_Status   *status = NULL;
  MPI_Comm  mpi_comm = cs_glob_mpi_comm;

  const int  loc_rank = CS_MAX(cs_glob_rank_id, 0);

  assert(cs_glob_n_ranks > 1);
  assert(c_edges != NULL);
  assert(s_edges != NULL);

  /* Allocate MPI buffers used for exchanging data */

  BFT_MALLOC(request, c_edges->n_ranks + s_edges->n_ranks, MPI_Request);
  BFT_MALLOC(status, c_edges->n_ranks + s_edges->n_ranks, MPI_Status);

  BFT_MALLOC(c_edge_tag, c_edges->n_elts, int);
  BFT_MALLOC(s_edge_tag, s_edges->n_elts, int);

  for (i = 0; i < c_edges->n_elts; i++)
    c_edge_tag[i] = 1; /* define as selected */

  for (i = 0; i < s_edges->n_elts; i++)
    s_edge_tag[i] = 1; /* define as selected */

  for (i = 0; i < c_edges->n_elts; i++) {

    vid1 = c_edges->array[2*i] - 1;
    vid2 = c_edges->array[2*i+1] - 1;
    edge_id = _get_edge_id(vid1, vid2, sel_v2v_idx, sel_v2v_lst);

    if (edge_id == -1)
      c_edge_tag[i] = 0; /* unselect this coupled edge */

  } /* End of loop on c_edges */

  /* Exchange between ranks the status of the coupled edges.
     If one receive a tag equal to 0 => delete the related single edge */

  request_count = 0;

  for (i = 0; i < s_edges->n_ranks; i++) {

    int  distant_rank = s_edges->ranks[i];
    int  length = s_edges->index[i+1] - s_edges->index[i];
    int  *recv_buf = s_edge_tag + s_edges->index[i];

    MPI_Irecv(recv_buf,
              length,
              MPI_INT,
              distant_rank,
              distant_rank,
              mpi_comm,
              &(request[request_count++]));

  }

  /* Send data to distant ranks */

  for (i = 0; i < c_edges->n_ranks; i++) {

    int  distant_rank = c_edges->ranks[i];
    int  length = c_edges->index[i+1] - c_edges->index[i];
    int  *send_buf = c_edge_tag + c_edges->index[i];

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

  /* Delete unselected elements */

  shift = 0;
  if (c_edges-> n_elts > 0) {

    save = c_edges->index[0];

    for (i = 0; i < c_edges->n_ranks; i++) {

      for (j = save; j < c_edges->index[i+1]; j++) {
        if (c_edge_tag[j] == 1) {
          c_edges->array[2*shift] = c_edges->array[2*j];
          c_edges->array[2*shift+1] = c_edges->array[2*j+1];
          shift++;
        }
      }
      save = c_edges->index[i+1];
      c_edges->index[i+1] = shift;

    }
    c_edges->n_elts = shift;

  } /* c_edges->n_elts > 0 */

  shift = 0;
  if (s_edges->n_elts > 0) {

    save = s_edges->index[0];

    for (i = 0; i < s_edges->n_ranks; i++) {

      for (j = save; j < s_edges->index[i+1]; j++) {
        if (s_edge_tag[j] == 1) {
          s_edges->array[2*shift] = s_edges->array[2*j];
          s_edges->array[2*shift+1] = s_edges->array[2*j+1];
          shift++;
        }
      }
      save = s_edges->index[i+1];
      s_edges->index[i+1] = shift;

    }
    s_edges->n_elts = shift;

  } /* s_edges->n_elts > 0 */

  /* Free memory */

  BFT_FREE(c_edge_tag);
  BFT_FREE(s_edge_tag);
  BFT_FREE(request);
  BFT_FREE(status);
}

/*----------------------------------------------------------------------------
 * Get the full selection of vertices to extract. Only in parallel. Somme
 * vertices may have been selected on another rank and not on the local rank
 * but you have to take them into account to have a good update of the mesh
 *
 * parameters:
 *  b_f2v_idx     <-- boundary "face -> vertex" connect. index
 *  b_f2v_lst     <-- boundary "face -> vertex" connect. list
 *  i_f2v_idx     <-- interior "face -> vertex" connect. index
 *  i_f2v_lst     <-- interior "face -> vertex" connect. list
 *  n_vertices    <-- number of vertices in the parent mesh
 *  vtx_tag       <-- tag on vertices. 1 if selected, 0 otherwise
 *  ifs           <-- pointer to a cs_interface_set_t struct.
 *  i_face_cells  <-- interior face -> cells connect.
 *  join_select   <-> pointer to a fvm_join_selection_t structure
 *  verbosity     <-- verbosity level
 *---------------------------------------------------------------------------*/

static void
_get_missing_edges(cs_lnum_t            b_f2v_idx[],
                   cs_lnum_t            b_f2v_lst[],
                   cs_lnum_t            i_f2v_idx[],
                   cs_lnum_t            i_f2v_lst[],
                   cs_lnum_t            n_vertices,
                   cs_lnum_t            vtx_tag[],
                   cs_interface_set_t  *ifs,
                   cs_lnum_2_t          i_face_cells[],
                   cs_join_select_t    *selection,
                   int                  verbosity)
{
  cs_gnum_t  n_l_elts, n_g_elts;

  cs_lnum_t  *sel_v2v_idx = NULL, *sel_v2v_lst = NULL;

  assert(ifs != NULL);
  assert(selection != NULL);

  /* Define single edge element */

  _get_select_v2v_connect(n_vertices,
                          selection,
                          b_f2v_idx,
                          b_f2v_lst,
                          &sel_v2v_idx,
                          &sel_v2v_lst);

  _add_single_edges(ifs,
                    vtx_tag,
                    selection,
                    sel_v2v_idx,
                    sel_v2v_lst,
                    b_f2v_idx,
                    b_f2v_lst,
                    i_f2v_idx,
                    i_f2v_lst,
                    i_face_cells,
                    selection->s_edges,
                    verbosity);

  n_l_elts = selection->s_edges->n_elts;

  MPI_Allreduce(&n_l_elts, &n_g_elts, 1, CS_MPI_GNUM,
                MPI_SUM, cs_glob_mpi_comm);

  if (n_g_elts > 0) {

    bft_printf(_("  Global number of single edges found:    %6llu\n"),
               (unsigned long long)n_g_elts);
    bft_printf_flush();
    selection->do_single_sync = true;

    _add_coupled_edges(ifs,
                       selection->s_edges,
                       selection->c_edges,
                       verbosity);

    _filter_edge_element(selection,
                         sel_v2v_idx,
                         sel_v2v_lst);

  } /* End if n_g_s_elts > 0 */

  if (cs_glob_join_log != NULL && verbosity > 3) {
    int  i, j;
    cs_join_sync_t  *s_edges = selection->s_edges;
    cs_join_sync_t  *c_edges = selection->c_edges;
    FILE  *logfile = cs_glob_join_log;

    fprintf(logfile, "\n  After possible filtering\n");

    fprintf(logfile, "\n  Coupled edges for the joining operation: (%p)\n",
            (void *)c_edges);
    fprintf(logfile,
            "  Coupled edges: n_elts : %8d\n"
            "  Coupled edges: n_ranks: %8d\n",
            c_edges->n_elts, c_edges->n_ranks);

    if (c_edges->n_elts > 0) {
      int  shift;
      for (i = 0, shift = 0; i < c_edges->n_ranks; i++) {
        for (j = c_edges->index[i]; j < c_edges->index[i+1]; j++) {
          fprintf(logfile, " %9d | %6d | (%9d, %9d)\n",
                  shift, c_edges->ranks[i],
                  c_edges->array[2*j], c_edges->array[2*j+1]);
          shift++;
        }
      }
      fprintf(logfile, "\n");
    }

    fprintf(logfile, "\n  Single edges for the joining operation: (%p)\n",
            (void *)s_edges);
    fprintf(logfile,
            "  Single edges: n_elts : %8d\n"
            "  Single edges: n_ranks: %8d\n",
            s_edges->n_elts, s_edges->n_ranks);

    if (s_edges->n_elts > 0) {
      for (i = 0; i < s_edges->n_ranks; i++)
        for (j = s_edges->index[i]; j < s_edges->index[i+1]; j++)
          fprintf(logfile, " %9d | %6d | (%9d, %9d)\n",
                  j, s_edges->ranks[i], s_edges->array[2*j],
                  s_edges->array[2*j+1]);
      fprintf(logfile, "\n");
    }

    fflush(logfile);
    bft_printf_flush();
   }

  /* Free memory */

  BFT_FREE(sel_v2v_idx);
  BFT_FREE(sel_v2v_lst);

}
#endif /* defined(HAVE_MPI) */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Create and initialize a cs_join_t structure.
 *
 * parameters:
 *   join_number   <-- number related to the joining operation
 *   sel_criteria  <-- boundary face selection criteria
 *   fraction      <-- value of the fraction parameter
 *   plane         <-- value of the plane parameter
 *   perio_type    <-- periodicity type (FVM_PERIODICITY_NULL if none)
 *   perio_matrix  <-- periodicity transformation matrix
 *   verbosity     <-- level of verbosity required
 *   visualization <-- level of visualization required
 *   preprocessing <-- is joining part of the preprocessing stage ?
 *
 * returns:
 *   a pointer to a new allocated cs_join_t structure
 *---------------------------------------------------------------------------*/

cs_join_t *
cs_join_create(int                      join_number,
               const char              *sel_criteria,
               float                    fraction,
               float                    plane,
               fvm_periodicity_type_t   perio_type,
               double                   perio_matrix[3][4],
               int                      verbosity,
               int                      visualization,
               bool                     preprocessing)
{
  size_t  l;

  cs_join_t  *join = NULL;

  /* Check main parameter values */

  if (fraction < 0.0 || fraction >= 1.0)
    bft_error(__FILE__, __LINE__, 0,
              _("Mesh joining:"
                "  Forbidden value for the fraction parameter.\n"
                "  It must be between [0.0, 1.0[ and is here: %f\n"),
              fraction);

  if (plane < 0.0 || plane >= 90.0)
    bft_error(__FILE__, __LINE__, 0,
              _("Mesh joining:"
                "  Forbidden value for the plane parameter.\n"
                "  It must be between [0, 90] and is here: %f\n"),
              plane);

  /* Initialize structure */

  BFT_MALLOC(join, 1, cs_join_t);

  join->selection = NULL;

  join->param = _join_param_define(join_number,
                                   fraction,
                                   plane,
                                   perio_type,
                                   perio_matrix,
                                   verbosity,
                                   visualization,
                                   preprocessing);

  join->stats = _join_stats_init();

  join->log_name = NULL;

  /* Copy the selection criteria for future use */

  l = strlen(sel_criteria);
  BFT_MALLOC(join->criteria, l + 1, char);
  strcpy(join->criteria, sel_criteria);

  /* Initialize log file if necessary */

  if (verbosity > 2) {
    char logname[80];
    char dir[] = "log";
    char rank_add[16] = "";
    char perio_add[16] = "";
    if (cs_file_isdir(dir) == 0) {
      if (cs_glob_rank_id < 1)
        if (cs_file_mkdir_default(dir) != 0)
          bft_error(__FILE__, __LINE__, 0,
                    _("The log directory cannot be created"));
#if defined(HAVE_MPI)
      if (cs_glob_n_ranks > 1)
        MPI_Barrier(cs_glob_mpi_comm); /* to avoid race conditions */
#endif
    }
    if (perio_type != FVM_PERIODICITY_NULL)
      strcpy(perio_add, "_perio");
    if (cs_glob_n_ranks > 1)
      sprintf(rank_add, "_r%04d", cs_glob_rank_id);
    sprintf(logname, "log%cjoin_%02d%s%s.log", DIR_SEPARATOR,
            join_number, perio_add, rank_add);
    BFT_MALLOC(join->log_name, strlen(logname) + 1, char);
    strcpy(join->log_name, logname);
  }

  return join;
}

/*----------------------------------------------------------------------------
 * Destroy a cs_join_t structure.
 *
 * parameters:
 *  join           <-> pointer to the cs_join_t structure to destroy
 *---------------------------------------------------------------------------*/

void
cs_join_destroy(cs_join_t  **join)
{
  if (*join != NULL) {

    cs_join_t  *_join = *join;

    BFT_FREE(_join->log_name);
    BFT_FREE(_join->criteria);

    BFT_FREE(_join);
    *join = NULL;

  }
}

/*----------------------------------------------------------------------------
 * Create and initialize a cs_join_select_t structure.
 *
 * parameters:
 *   selection_criteria <-- pointer to a cs_mesh_select_t structure
 *   perio_type         <-- periodicity type (FVM_PERIODICITY_NULL if none)
 *   verbosity          <-- level of verbosity required
 *
 * returns:
 *   pointer to a newly created cs_join_select_t structure
 *---------------------------------------------------------------------------*/

cs_join_select_t *
cs_join_select_create(const char              *selection_criteria,
                      fvm_periodicity_type_t   perio_type,
                      int                      verbosity)
{
  cs_lnum_t  i;

  cs_lnum_t  *vtx_tag = NULL;
  cs_join_select_t  *selection = NULL;
  cs_lnum_t  *order = NULL, *ordered_faces = NULL;
  cs_interface_set_t  *ifs = NULL;
  cs_mesh_t  *mesh = cs_glob_mesh;
  FILE  *logfile = cs_glob_join_log;

  const int  n_ranks = cs_glob_n_ranks;

  assert(mesh != NULL);

  /* Initialize cs_join_select_t struct. */

  BFT_MALLOC(selection, 1, cs_join_select_t);

  selection->n_init_b_faces = mesh->n_b_faces;
  selection->n_init_i_faces = mesh->n_i_faces;
  selection->n_init_vertices = mesh->n_vertices;

  selection->n_faces = 0;
  selection->n_g_faces = 0;

  selection->faces = NULL;
  selection->compact_face_gnum = NULL;
  selection->compact_rank_index = NULL;

  selection->n_vertices = 0;
  selection->n_g_vertices = 0;
  selection->vertices = NULL;

  selection->n_b_adj_faces = 0;
  selection->b_adj_faces = NULL;
  selection->n_i_adj_faces = 0;
  selection->i_adj_faces = NULL;

  selection->b_face_state = NULL;
  selection->i_face_state = NULL;

  selection->n_couples = 0;
  selection->per_v_couples = NULL;

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

  /* Extract selected boundary faces */

  BFT_MALLOC(selection->faces, mesh->n_b_faces, cs_lnum_t);

  cs_selector_get_b_face_num_list(selection_criteria,
                                  &(selection->n_faces),
                                  selection->faces);

  /* In case of periodicity, ensure no isolated faces are
     selected */

  if (perio_type != FVM_PERIODICITY_NULL) {
    cs_lnum_t j = 0;
    for (i = 0; i < selection->n_faces; i++) {
      cs_lnum_t f_id = selection->faces[i]-1;
      if (mesh->b_face_cells[f_id] > -1)
        selection->faces[j++] = f_id+1;
    }
    selection->n_faces = j;
  }

  BFT_MALLOC(order, selection->n_faces, cs_lnum_t);
  BFT_MALLOC(ordered_faces, selection->n_faces, cs_lnum_t);

  cs_order_gnum_allocated(selection->faces, NULL, order, selection->n_faces);

  for (i = 0; i < selection->n_faces; i++)
    ordered_faces[i] = selection->faces[order[i]];

  BFT_FREE(order);
  BFT_FREE(selection->faces);
  selection->faces = ordered_faces;

  if (n_ranks == 1)
    selection->n_g_faces = selection->n_faces;

#if defined(HAVE_MPI)
  if (n_ranks > 1) { /* Parallel treatment */

    cs_gnum_t n_l_faces = selection->n_faces;

    MPI_Allreduce(&n_l_faces, &(selection->n_g_faces),
                  1, CS_MPI_GNUM, MPI_SUM, cs_glob_mpi_comm);

  }
#endif

  if (verbosity > 0)
    bft_printf
      (_("  Global number of boundary faces selected for joining: %10llu\n"),
       (unsigned long long)selection->n_g_faces);

  /* Define a compact global numbering on selected boundary faces and
     build an index on ranks on this compact numbering */

  _compact_face_gnum_selection(selection->n_faces,
                               &(selection->compact_face_gnum),
                               &(selection->compact_rank_index));

  assert(selection->n_g_faces == selection->compact_rank_index[n_ranks]);

  /* Extract selected vertices from the selected boundary faces */

  cs_join_extract_vertices(selection->n_faces,
                           selection->faces,
                           mesh->b_face_vtx_idx,
                           mesh->b_face_vtx_lst,
                           mesh->n_vertices,
                           &(selection->n_vertices),
                           &(selection->vertices));

#if defined(HAVE_MPI)
  if (n_ranks > 1) { /* Search for missing vertices */

    assert(mesh->global_vtx_num != NULL);

    ifs = cs_interface_set_create(mesh->n_vertices,
                                  NULL,
                                  mesh->global_vtx_num,
                                  NULL,
                                  0,
                                  NULL,
                                  NULL,
                                  NULL);

    assert(ifs != NULL);

    _get_missing_vertices(mesh->n_vertices,
                          ifs,
                          &vtx_tag,
                          selection,
                          verbosity);

  }
#endif

  /* Extract list of boundary faces contiguous to the selected vertices  */

  _extract_contig_faces(mesh->n_vertices,
                        selection,
                        mesh->n_b_faces,
                        mesh->b_face_vtx_idx,
                        mesh->b_face_vtx_lst,
                        &(selection->n_b_adj_faces),
                        &(selection->b_adj_faces));

  /* Remove boundary faces already defined in selection->faces */

  cs_join_clean_selection(&(selection->n_b_adj_faces),
                          &(selection->b_adj_faces),
                          selection->n_faces,
                          selection->faces);

  /* Extract list of interior faces contiguous to the selected vertices */

  _extract_contig_faces(mesh->n_vertices,
                        selection,
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

    _get_missing_edges(mesh->b_face_vtx_idx,
                       mesh->b_face_vtx_lst,
                       mesh->i_face_vtx_idx,
                       mesh->i_face_vtx_lst,
                       mesh->n_vertices,
                       vtx_tag,
                       ifs,
                       mesh->i_face_cells,
                       selection,
                       verbosity);

    BFT_FREE(vtx_tag);
    cs_interface_set_destroy(&ifs);

  }
#endif

  /* Face state setting */

  BFT_MALLOC(selection->b_face_state, mesh->n_b_faces, cs_join_state_t);
  BFT_MALLOC(selection->i_face_state, mesh->n_i_faces, cs_join_state_t);

  for (i = 0; i < mesh->n_b_faces; i++)
    selection->b_face_state[i] = CS_JOIN_STATE_UNDEF;

  for (i = 0; i < mesh->n_i_faces; i++)
    selection->i_face_state[i] = CS_JOIN_STATE_UNDEF;

  for (i = 0; i < selection->n_faces; i++)
    selection->b_face_state[selection->faces[i]-1] = CS_JOIN_STATE_ORIGIN;

  for (i = 0; i < selection->n_b_adj_faces; i++)
    selection->b_face_state[selection->b_adj_faces[i]-1] = CS_JOIN_STATE_ORIGIN;

  for (i = 0; i < selection->n_i_adj_faces; i++)
    selection->i_face_state[selection->i_adj_faces[i]-1] = CS_JOIN_STATE_ORIGIN;

  /* Display information according to the level of verbosity */

  if (verbosity > 2) {

    assert(logfile != NULL);

    fprintf(logfile,
            "\n  Local information about selection structure:\n");
    fprintf(logfile,
            "    number of faces:               %8d\n",
            selection->n_faces);
    fprintf(logfile,
            "    number of vertices:            %8d\n",
            selection->n_vertices);
    fprintf(logfile,
            "    number of adj. boundary faces: %8d\n",
            selection->n_b_adj_faces);
    fprintf(logfile,
            "    number of adj. interior faces: %8d\n",
            selection->n_i_adj_faces);

    if (selection->do_single_sync == true) {
      fprintf(logfile,
              "\n Information on single/coupled elements:\n");
      fprintf(logfile,
              "   Number of single vertices : %6d with %3d related ranks\n",
              selection->s_vertices->n_elts, selection->s_vertices->n_ranks);
      fprintf(logfile,
              "   Number of coupled vertices: %6d with %3d related ranks\n",
              selection->c_vertices->n_elts, selection->c_vertices->n_ranks);
      fprintf(logfile,
              "   Number of single edges    : %6d with %3d related ranks\n",
              selection->s_edges->n_elts, selection->s_edges->n_ranks);
      fprintf(logfile,
              "   Number of coupled edges   : %6d with %3d related ranks\n",
              selection->c_edges->n_elts, selection->c_edges->n_ranks);
    }

    if (verbosity > 3) {
      fprintf(logfile,
              "\n  Compact index on ranks for the selected faces:\n");
      for (i = 0; i < n_ranks + 1; i++)
        fprintf(logfile,
                " %5d | %11llu\n", i,
                (unsigned long long)selection->compact_rank_index[i]);
      fprintf(logfile, "\n");
    }

    /* Debug level logging */

    if (verbosity > 4) {

      fprintf(logfile,
              "\n  Selected faces for the joining operation:\n");
      for (i = 0; i < selection->n_faces; i++)
        fprintf(logfile,
                " %9d | %9d | %10llu\n",
                i, selection->faces[i],
                (unsigned long long)selection->compact_face_gnum[i]);
      fprintf(logfile, "\n");
    }

    if (verbosity > 4) {
      fprintf(logfile,
              "\n  Selected vertices for the joining operation:\n");
      for (i = 0; i < selection->n_vertices; i++)
        fprintf(logfile,
                " %9d | %9d\n", i, selection->vertices[i]);
      fprintf(logfile, "\n");
    }

    if (verbosity > 4) {
      fprintf(logfile,
              "\n  Contiguous boundary faces for the joining operation:\n");
      for (i = 0; i < selection->n_b_adj_faces; i++)
        fprintf(logfile,
                " %9d | %9d\n", i, selection->b_adj_faces[i]);
      fprintf(logfile, "\n");
    }

    if (verbosity > 4) {
      fprintf(logfile,
              "\n  Contiguous interior faces for the joining operation:\n");
      for (i = 0; i < selection->n_i_adj_faces; i++)
        fprintf(logfile, " %9d | %9d\n", i, selection->i_adj_faces[i]);
      fprintf(logfile, "\n");
    }

    fflush(logfile);

    bft_printf_flush();

  } /* End if verbosity > 2 */

  return  selection;
}

/*----------------------------------------------------------------------------
 * Destroy a cs_join_select_t structure.
 *
 * parameters:
 *   param       <-- user-defined joining parameters
 *   join_select <-- pointer to pointer to structure to destroy
 *---------------------------------------------------------------------------*/

void
cs_join_select_destroy(cs_join_param_t     param,
                       cs_join_select_t  **join_select)
{
  if (*join_select != NULL) {

    cs_join_select_t *_js = *join_select;

    BFT_FREE(_js->faces);
    BFT_FREE(_js->compact_face_gnum);
    BFT_FREE(_js->compact_rank_index);
    BFT_FREE(_js->vertices);
    BFT_FREE(_js->b_adj_faces);
    BFT_FREE(_js->i_adj_faces);

    BFT_FREE(_js->b_face_state);
    BFT_FREE(_js->i_face_state);

    if (param.perio_type != FVM_PERIODICITY_NULL)
      BFT_FREE(_js->per_v_couples);

    _destroy_join_sync(&(_js->s_vertices));
    _destroy_join_sync(&(_js->c_vertices));
    _destroy_join_sync(&(_js->s_edges));
    _destroy_join_sync(&(_js->c_edges));

    BFT_FREE(*join_select);
    *join_select = NULL;

  }
}

/*----------------------------------------------------------------------------
 * Extract vertices from a selection of faces.
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

void
cs_join_extract_vertices(cs_lnum_t         n_select_faces,
                         const cs_lnum_t  *select_faces,
                         const cs_lnum_t  *f2v_idx,
                         const cs_lnum_t  *f2v_lst,
                         cs_lnum_t         n_vertices,
                         cs_lnum_t        *n_select_vertices,
                         cs_lnum_t        *select_vertices[])
{
  int  i, j, face_id;

  cs_lnum_t  _n_select_vertices = 0;
  cs_lnum_t  *counter = NULL, *_select_vertices = NULL;

  if (n_select_faces > 0) {

    BFT_MALLOC(counter, n_vertices, cs_lnum_t);

    for (i = 0; i < n_vertices; i++)
      counter[i] = 0;

    for (i = 0; i < n_select_faces; i++) {

      face_id = select_faces[i] - 1;

      for (j = f2v_idx[face_id]; j < f2v_idx[face_id+1]; j++)
        counter[f2v_lst[j]] = 1;

    }

    for (i = 0; i < n_vertices; i++)
      _n_select_vertices += counter[i];

    BFT_MALLOC(_select_vertices, _n_select_vertices, cs_lnum_t);

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
 * Eliminate redundancies found between two lists of elements.
 * Delete elements in elts[] and keep elements in the reference list.
 *
 * parameters:
 *  n_elts      <->  number of elements in the list to clean
 *  elts        <->  list of elements in the list to clean
 *  n_ref_elts  <--  number of elements in the reference list
 *  ref_elts    <--  list of reference elements
 *---------------------------------------------------------------------------*/

void
cs_join_clean_selection(cs_lnum_t  *n_elts,
                        cs_lnum_t  *elts[],
                        cs_lnum_t   n_ref_elts,
                        cs_lnum_t   ref_elts[])
{
  cs_lnum_t  i = 0, j = 0;
  cs_lnum_t  _n_elts = 0;
  cs_lnum_t  *_elts = *elts;

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

  BFT_REALLOC(_elts, _n_elts, cs_lnum_t);

  *n_elts = _n_elts;
  *elts = _elts;
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
cs_join_build_edges_idx(cs_lnum_t        n_faces,
                        const cs_lnum_t  faces[],
                        const cs_lnum_t  f2v_idx[],
                        const cs_lnum_t  f2v_lst[],
                        cs_lnum_t        v2v_idx[])
{
  cs_lnum_t  i, j, v1, v2, fid, s, e;

  /* Loop on all selected faces. No need to loop on other faces because
     the selected vertices are all found with this only step. */

  for (i = 0; i < n_faces; i++) {

    fid = faces[i] - 1;
    s = f2v_idx[fid];
    e = f2v_idx[fid+1];

    for (j = s; j < e - 1; j++) { /* scan edges */

      v1 = f2v_lst[j] + 1;
      v2 = f2v_lst[j+1] + 1;

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

    v1 = f2v_lst[e-1] + 1;
    v2 = f2v_lst[s] + 1;

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
cs_join_build_edges_lst(cs_lnum_t        n_faces,
                        const cs_lnum_t  faces[],
                        const cs_lnum_t  f2v_idx[],
                        const cs_lnum_t  f2v_lst[],
                        cs_lnum_t        count[],
                        const cs_lnum_t  v2v_idx[],
                        cs_lnum_t        v2v_lst[])
{
  cs_lnum_t  i, j, v1_id, v2_id, fid, s, e, shift;

  for (i = 0; i < n_faces; i++) {

    fid = faces[i] - 1;
    s = f2v_idx[fid];
    e = f2v_idx[fid+1];

    for (j = s; j < e - 1; j++) { /* Scan edges */

      v1_id = f2v_lst[j];
      v2_id = f2v_lst[j+1];

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

    v1_id = f2v_lst[e-1];
    v2_id = f2v_lst[s];

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
