/*============================================================================
 * Split faces during the joining operation
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
#include "fvm_periodicity.h"
#include "cs_block_dist.h"
#include "cs_mesh.h"
#include "cs_order.h"
#include "cs_search.h"
#include "cs_join_post.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *---------------------------------------------------------------------------*/

#include "cs_join_split.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Macro definitions
 *===========================================================================*/

/*============================================================================
 * Structure and type definitions
 *===========================================================================*/

/* enum structure to define a category for each kind of problem
   which can be encountered during the face splitting operation */

typedef enum {

  NO_SPLIT_ERROR,
  OPEN_CYCLE_ERROR,
  EDGE_TRAVERSED_TWICE_ERROR,
  LOOP_LIMIT_ERROR,
  MAX_SPLIT_ERROR

} cs_join_split_error_t;

/* ----------------------------------------------------------------- */
/* Definition of a structure used to build the new face connectivity */
/* ----------------------------------------------------------------- */

typedef struct {

  cs_lnum_t        n_faces;
  cs_lnum_t       *face_index;       /* Face -> Subface index */
  cs_join_rset_t  *subface_index;    /* Subface -> vertex connect. index */
  cs_join_rset_t  *subface_connect;  /* Subface -> vertex connect. list */

  /* The two following arrays are built after the face splitting. So, they
     don't need to be resizable */

  cs_gnum_t   n_g_subfaces;     /* Global number of subfaces after splitting */
  cs_gnum_t  *subface_gconnect; /* Subface -> glob. vertex list */
  cs_gnum_t  *subface_gnum;     /* Subface global numbering */

} face_builder_t;

/*============================================================================
 * Global variable definitions
 *===========================================================================*/

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
 *---------------------------------------------------------------------------*/

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
 *---------------------------------------------------------------------------*/

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
 * Get the norm of a vector.
 *
 * parameters:
 *  v      <->  vector to work with.
 *---------------------------------------------------------------------------*/

inline static double
_norm(double   v[])
{
  return  sqrt(_dot_product(v, v));
}

/*----------------------------------------------------------------------------
 * Compute the cosine of two vectors.
 *
 * parameters:
 *  v1     <--  first vector
 *  v2     <--  second vector
 *
 * returns:                   / \
 *  the resulting cosine of (v1.v2)
 *---------------------------------------------------------------------------*/

inline static double
_cosine(double   v1[],
        double   v2[])
{
  double  dprod = _dot_product(v1, v2);
  double  n1 = _norm(v1);
  double  n2 = _norm(v2);

  assert(n1 > 0.0);
  assert(n2 > 0.0);

  return dprod / (n1 * n2);
}

/*----------------------------------------------------------------------------
 * Compute the bounding box related to a face in which reconstruction
 * is allowed.
 *
 * parameters:
 *  m      <--  pointer to a cs_join_mesh_t structure
 *  fid    <--  face id
 *  f_min  -->  min coordinates for bounding box
 *  f_max  -->  max coordinates for bounding box
 *---------------------------------------------------------------------------*/

inline static void
_face_bbox(const cs_join_mesh_t  *m,
           int                    fid,
           double                 f_min[],
           double                 f_max[])
{
  int  i, k;

  const double  eps = 1e-5;
  const cs_join_vertex_t  *vertices = m->vertices;

  for (i = m->face_vtx_idx[fid]; i < m->face_vtx_idx[fid+1]; i++) {

    int  vid = m->face_vtx_lst[i];
    double  tol = (1 + eps) * vertices[vid].tolerance;

    for (k = 0; k < 3; k++) {

      double  v_min = vertices[vid].coord[k] - tol;
      double  v_max = vertices[vid].coord[k] + tol;

      f_min[k] = CS_MIN(f_min[k], v_min);
      f_max[k] = CS_MAX(f_max[k], v_max);

    }

  } /* End of loop on face vertices */

}

/*----------------------------------------------------------------------------
 * Locally renumber an indexed list of global elements according to an
 * ordering array.
 *
 * parameters:
 *   n_elts    <-- number of elements in index
 *   order     <-- ordering array
 *   index     <-- unordered index on global elements (index[0] = 0)
 *   glist     <-- unordered indexed list of global elements
 *   new_index --> ordered index on global elements (index[0] = 0)
 *   new_glist --> ordered indexed list of global elements
 *---------------------------------------------------------------------------*/

static void
_renumber_local_ordered_i(cs_lnum_t         n_elts,
                          const cs_lnum_t   order[],
                          const cs_lnum_t   index[],
                          cs_gnum_t         glist[],
                          cs_lnum_t        *new_index[],
                          cs_gnum_t        *new_glist[])
{
  cs_lnum_t  i, j, k, o_id;

  cs_lnum_t  *_new_index = NULL;
  cs_gnum_t   *_new_glist = NULL;

  assert(index[0] == 0); /* case index[0] = 1 coulb be coded in the future */

  /* Build a new index */

  BFT_MALLOC(_new_index, n_elts + 1, cs_lnum_t);

  for (i = 0; i < n_elts; i++) {
    o_id = order[i];
    _new_index[i+1] = index[o_id+1] - index[o_id];
  }

  _new_index[0] = 0;
  for (i = 0; i < n_elts; i++)
    _new_index[i+1] += _new_index[i];

  /* Build a new list */

  BFT_MALLOC(_new_glist, _new_index[n_elts], cs_gnum_t);

  for (i = 0; i < n_elts; i++) {

    o_id = order[i];

    for (j = index[o_id], k = _new_index[i]; j < index[o_id+1]; j++, k++)
      _new_glist[k] = glist[j];
  }

  /* Return pointers */

  *new_index = _new_index;
  *new_glist = _new_glist;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Define send_rank_index and send_faces to prepare the exchange of new faces
 * between mesh structures.
 *
 * parameters:
 *   o2n_hist        <-- old -> new global face numbering
 *   gnum_rank_index <-- index on ranks for the old global face numbering
 *   send_rank_index --> index on ranks for sending face
 *   send_faces      --> list of face ids to send
 *---------------------------------------------------------------------------*/

static void
_get_faces_to_send(const cs_join_gset_t  *o2n_hist,
                   const cs_gnum_t        gnum_rank_index[],
                   cs_lnum_t             *send_rank_index[],
                   cs_lnum_t             *send_faces[])
{
  cs_lnum_t  i, j, rank, start, end;

  cs_lnum_t  reduce_size = 0;
  cs_lnum_t  *_send_rank_index = NULL, *_send_faces = NULL, *reduce_ids = NULL;
  cs_gnum_t  *reduce_index = NULL;
  cs_join_gset_t  *new_face_rank = NULL;

  const int  n_ranks = cs_glob_n_ranks;

  /* Sanity checks */

  assert(gnum_rank_index != NULL);
  assert(o2n_hist != NULL);
  assert(n_ranks > 1);

  new_face_rank = cs_join_gset_create(n_ranks);

  for (i = 0; i < n_ranks; i++)
    new_face_rank->g_elts[i] = 0; /* used to store local count */

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

  /* Count number of ranks associated to each new face */

  for (i = 0; i < o2n_hist->n_elts; i++) {

    int  reduce_rank = cs_search_gindex_binary(reduce_size,
                                               o2n_hist->g_elts[i],
                                               reduce_index);

    assert(reduce_rank != -1);
    assert(reduce_rank < reduce_size);

    rank = reduce_ids[reduce_rank];
    new_face_rank->index[rank+1] += o2n_hist->index[i+1] - o2n_hist->index[i];

  }

  for (i = 0; i < n_ranks; i++)
    new_face_rank->index[i+1] += new_face_rank->index[i];

  BFT_MALLOC(new_face_rank->g_list, new_face_rank->index[n_ranks], cs_gnum_t);

  /* Fill the list of ranks */

  for (i = 0; i < o2n_hist->n_elts; i++) {

    int  reduce_rank = cs_search_gindex_binary(reduce_size,
                                               o2n_hist->g_elts[i],
                                               reduce_index);

    assert(reduce_rank != -1);
    assert(reduce_rank < reduce_size);

    rank = reduce_ids[reduce_rank];
    start = o2n_hist->index[i];
    end = o2n_hist->index[i+1];

    for (j = start; j < end; j++) {

      cs_lnum_t  shift =  new_face_rank->index[rank]
                        + new_face_rank->g_elts[rank];
      cs_lnum_t  new_fid = o2n_hist->g_list[j] - 1;

      new_face_rank->g_list[shift] = new_fid;
      new_face_rank->g_elts[rank] += 1;

    } /* End of loop on new faces */

  } /* End of loop on initial faces */

  /* Free memory */

  BFT_FREE(reduce_ids);
  BFT_FREE(reduce_index);

  cs_join_gset_clean(new_face_rank);

  /* Define arrays to return */

  BFT_MALLOC(_send_rank_index, n_ranks + 1, cs_lnum_t);

  for (i = 0; i < n_ranks + 1; i++)
    _send_rank_index[i] = new_face_rank->index[i];

  BFT_MALLOC(_send_faces, _send_rank_index[n_ranks], cs_lnum_t);

  for (i = 0; i < _send_rank_index[n_ranks]; i++)
    _send_faces[i] = new_face_rank->g_list[i];

  cs_join_gset_destroy(&new_face_rank);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (cs_glob_join_log != NULL) {
    fprintf(cs_glob_join_log,
            "\n Exchange to do after the splitting operation:\n");
    for (i = 0; i < n_ranks; i++) {
      start = _send_rank_index[i];
      end = _send_rank_index[i+1];
      fprintf(cs_glob_join_log,
              " Send to rank %5d (n = %10d):", i, end - start);
      for (j = start; j < end; j++)
        fprintf(cs_glob_join_log, " %d ", _send_faces[j]);
      fprintf(cs_glob_join_log, "\n");
    }
  }
#endif

  /* Return pointers */

  *send_rank_index = _send_rank_index;
  *send_faces = _send_faces;
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Define the head_edges and ext_edges lists.
 *
 * Use a permutation number "perm" to build lists with a shift.
 *
 * Note: perm < n_face_vertices
 *
 * parameters:
 *   face_id    <-- face_id of the current in the cs_join_mesh_t struct.
 *   mesh       <-- cs_join_mesh_t structure
 *   edges      <-- cs_join_edges_t structure
 *   head_edges <-> head edegs list struct.
 *   ext_edges  <-> exterior edge list struct.
 *   perm       <-- permutation shift
 *---------------------------------------------------------------------------*/

static void
_define_head_and_ext_edges(cs_lnum_t               face_id,
                           const cs_join_mesh_t   *mesh,
                           const cs_join_edges_t  *edges,
                           cs_join_rset_t         *head_edges,
                           cs_join_rset_t         *ext_edges,
                           cs_lnum_t               perm)
{
  cs_lnum_t  i, j, k, shift;
  cs_lnum_t  couple[2];

  cs_lnum_t  start_id = mesh->face_vtx_idx[face_id];
  cs_lnum_t  end_id = mesh->face_vtx_idx[face_id+1];
  cs_lnum_t  n_face_vertices = end_id - start_id;

  assert(perm < n_face_vertices);

  for (k = 0, i = start_id; i < end_id-1; i++) {

    couple[0] = mesh->face_vtx_lst[i];
    couple[1] = mesh->face_vtx_lst[i+1];

    for (j = edges->vtx_idx[couple[0]];
         j < edges->vtx_idx[couple[0]+1]; j++)
      if (edges->adj_vtx_lst[j] == couple[1])
        break;

    shift = (k + perm) % n_face_vertices;
    head_edges->array[shift] = edges->edge_lst[j];
    ext_edges->array[shift] = edges->edge_lst[j];
    k++;

  } /* End of loop on face vertices */

  /* Last couple : last vertex and first vertex */

  couple[0] = mesh->face_vtx_lst[end_id-1];
  couple[1] = mesh->face_vtx_lst[start_id];

  for (j = edges->vtx_idx[couple[0]];
       j < edges->vtx_idx[couple[0]+1]; j++)
    if (edges->adj_vtx_lst[j] == couple[1])
      break;

  shift = (k + perm) % n_face_vertices;
  head_edges->array[shift] = edges->edge_lst[j];
  ext_edges->array[shift] = edges->edge_lst[j];
  k++;

  assert(k == n_face_vertices);

  head_edges->n_elts = n_face_vertices;
  ext_edges->n_elts = n_face_vertices;

}

/*----------------------------------------------------------------------------
 * Create a face_builder_t structure.
 *
 * parameters:
 *   n_faces <--number of faces associated to the structure
 *
 * returns:
 *   a new allocated face_builder_t structure.
 *---------------------------------------------------------------------------*/

static face_builder_t *
_create_face_builder(cs_lnum_t  n_faces)
{
  cs_lnum_t  i;

  face_builder_t  *builder = NULL;

  BFT_MALLOC(builder, 1, face_builder_t);

  builder->n_faces = n_faces;

  BFT_MALLOC(builder->face_index, n_faces + 1, cs_lnum_t);

  for (i = 0; i < n_faces + 1; i++)
    builder->face_index[i] = 0;

  builder->subface_index = cs_join_rset_create(n_faces + 1);
  builder->subface_index->array[0] = 0;

  /* Initialize the size of the sub-face connectivity as if all faces drives
     to one triangular sub-face which is the minimal connectivity size
     possible */

  builder->subface_connect = cs_join_rset_create(15);
  builder->subface_gnum = NULL;
  builder->subface_gconnect = NULL;

  return builder;
}

/*----------------------------------------------------------------------------
 * Destroy a face_builder_t  structure.
 *
 * parameters:
 *   builder <-- pointer to the structure to delete
 *
 * returns:
 *   a NULL pointer.
 *---------------------------------------------------------------------------*/

static face_builder_t *
_destroy_face_builder(face_builder_t  *builder)
{
  if (builder == NULL)
    return NULL;

  BFT_FREE(builder->face_index);

  cs_join_rset_destroy(&(builder->subface_index));
  cs_join_rset_destroy(&(builder->subface_connect));

  if (builder->subface_gnum != NULL)
    BFT_FREE(builder->subface_gnum);

  if (builder->subface_gconnect != NULL)
    BFT_FREE(builder->subface_gconnect);

  BFT_FREE(builder);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Dump a face and its subface connect. from a face_builder_t structure.
 *
 * parameters:
 *   face_id <-- face_id to dump
 *   builder <-- pointer to the structure to dump
 *   logfile <-- handle to log file
 *---------------------------------------------------------------------------*/

static void
_dump_face_builder(cs_lnum_t              face_id,
                   const face_builder_t  *builder,
                   FILE                  *logfile)
{
  cs_lnum_t  i, j;
  cs_lnum_t  subface_id = 0;
  cs_lnum_t  face_s = builder->face_index[face_id];
  cs_lnum_t  face_e = builder->face_index[face_id+1];
  cs_lnum_t  n_subfaces = face_e - face_s;

  fprintf(logfile, " Face %9d (n_subfaces: %d):\n", face_id+1, n_subfaces);

  for (i = face_s; i <face_e; i++, subface_id++) {

    cs_lnum_t  subface_s = builder->subface_index->array[i];
    cs_lnum_t  subface_e = builder->subface_index->array[i+1];

    if (builder->subface_gnum == NULL)
      fprintf(logfile, "   subface %4d: (%d, %d) -",
              subface_id, subface_s, subface_e);
    else
      fprintf(logfile, "   subface %4d (%10llu) - (%d, %d) -",
              subface_id,
              (unsigned long long)builder->subface_gnum[face_s + subface_id],
              subface_s, subface_e);

    if (builder->subface_gconnect == NULL) {
      for (j = subface_s; j < subface_e; j++)
        fprintf(logfile, " %d ", builder->subface_connect->array[j]);
    }
    else {
      for (j = subface_s; j < subface_e; j++)
        fprintf(logfile, " %llu ",
                (unsigned long long)builder->subface_gconnect[j]);
    }
    fprintf(logfile, "\n");

  }
  fflush(logfile);
}

/*----------------------------------------------------------------------------
 * Step for building subfaces.
 * Find the best face (in sense of coplanarity) among faces sharing the
 * current edge to test.
 *
 * parameters:
 *   param       <-- set of parameters for the joining operation
 *   eid         <-- id of the current edge to test
 *   fid         <-- fid of the current in the cs_join_mesh_t struct.
 *   fnorm       <-- normal of the current face to rebuild
 *   adj_fnorm   <-- normal of the best adj face (if found)
 *   face_normal <-- normal for each face of the mesh
 *   e2f_idx     <-- "edge -> face" connect. index
 *   e2f_lst     <-- "edge -> face" connect. list
 *
 * returns:
 *   true if an adjacent face was found else false
 *---------------------------------------------------------------------------*/

static bool
_find_best_adj_face(cs_join_param_t         param,
                    cs_lnum_t               eid,
                    cs_lnum_t               fid,
                    double                  fnorm[3],
                    double                  adj_fnorm[3],
                    const double            face_normal[],
                    const cs_lnum_t        *e2f_idx,
                    const cs_lnum_t        *e2f_lst)
{
  cs_lnum_t  j, k, adj_fid;
  double  dprod, dprod2, test_fnorm[3];

  int  best_fid = -1;
  double  max_plane = -DBL_MAX;
  bool  respect_plane = false;

  const double  plane2 = param.plane_criteria;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
#define _DBGTST 1
  bool  tst_dbg = ( fid==492 || fid==1039 || fid==744 || fid==262 ||
                         fid==546 || fid==1057 || fid==564 ? true : false);
#else
#define _DBGTST 0
#endif

  /*
    if (f2f_connect == NULL)  What we already do...
    else  To be implemented ....
  */

  /* Loop on faces sharing this edge */

  for (j = e2f_idx[eid]; j < e2f_idx[eid+1]; j++) {

    adj_fid = e2f_lst[j] - 1;

    if (adj_fid != fid) {

      for (k = 0; k < 3; k++)
        test_fnorm[k] = face_normal[3*adj_fid+k];

      dprod = _dot_product(test_fnorm, fnorm);
      dprod2 = dprod * dprod;

#if _DBGTST && defined(DEBUG) && !defined(NDEBUG)
      if (tst_dbg && cs_glob_join_log != NULL)
        fprintf(cs_glob_join_log,
                "\tAdjFace: %d dp: %10.8e dp2: %10.8e Vs plane2: %10.8e;"
                " AdjFNorm: [%g %g %g]; FNorm: [%g %g %g]\n",
                adj_fid+1, dprod, dprod2, plane2, test_fnorm[0],
                test_fnorm[1], test_fnorm[2], fnorm[0], fnorm[1], fnorm[2]);
#endif

      if (dprod2 > plane2) { /* A better face has been found */

        respect_plane = true;
        if (dprod2 > max_plane) {
          best_fid = adj_fid;
          max_plane = dprod2;
          if (dprod > 0.0) {
            for (k = 0; k < 3; k++)
              adj_fnorm[k] = face_normal[3*best_fid+k];
          }
          else {
            for (k = 0; k < 3; k++)
              adj_fnorm[k] = -face_normal[3*best_fid+k];
          }
        }

#if _DBGTST && defined(DEBUG) && !defined(NDEBUG)
        if (tst_dbg && cs_glob_join_log != NULL)
          fprintf(cs_glob_join_log,
                  "\t-> Adj face OK (best: %d, max: %10.8e)\n",
                  best_fid+1, max_plane);
#endif

      } /* dprod2 > plane2 */

    }
    else {  /* If adj_fid == fid => the best choice */

      assert(adj_fid == fid);
      respect_plane = true;

      for (k = 0; k < 3; k++)
        adj_fnorm[k] = fnorm[k];

#if _DBGTST && defined(DEBUG) && !defined(NDEBUG)
      if (tst_dbg && cs_glob_join_log != NULL)
        fprintf(cs_glob_join_log, "\t Adj face is current face\n");
#endif

      /* Anayway, it's the best choice  */
      return respect_plane;

    }

  } /* End of loop on faces sharing this edge */

#if _DBGTST && defined(DEBUG) && !defined(NDEBUG)
  if (tst_dbg && cs_glob_join_log != NULL && respect_plane)
    fprintf(cs_glob_join_log, "\tAdjFace: %d AdjFNorm: [%5.3e %5.3e %5.3e]\n",
            best_fid+1, adj_fnorm[0], adj_fnorm[1], adj_fnorm[2]);
#endif

  return respect_plane;
}

/*----------------------------------------------------------------------------
 * Step for building subfaces.
 * Find the next edge and vertex numbers from a given edge.
 *
 * parameters:
 *   param       <-- set of parameters for the joining operation
 *   fid         <-- fid of the current in the cs_join_mesh_t struct.
 *   vid1        <-- first vertex id of the current edge
 *   vid2        <-- second vertex id of the current edge
 *   face_normal <-- normal for each face of the mesh
 *   work        <-- cs_join_mesh_t structure
 *   e2f_idx     <-- "edge -> face" connect. index
 *   e2f_lst     <-- "edge -> face" connect. list
 *   next_edge   --> pointer to the next edge number found
 *   next_vertex --> pointer to the next vertex number found
 *
 * returns:
 *   an error code (enum: cs_join_split_error_t)
 *---------------------------------------------------------------------------*/

static cs_join_split_error_t
_find_next(cs_join_param_t         param,
           cs_lnum_t               fid,
           cs_lnum_t               vid1,
           cs_lnum_t               vid2,
           cs_real_t               max_coord[3],
           cs_real_t               min_coord[3],
           const double            face_normal[],
           const cs_join_mesh_t   *work,
           const cs_join_edges_t  *edges,
           const cs_lnum_t        *e2f_idx,
           const cs_lnum_t        *e2f_lst,
           cs_lnum_t              *next_edge,
           cs_lnum_t              *next_vertex)
{
  cs_lnum_t  i, j, k;
  double  norm, dprod, adj_fnorm[3], fnorm[3], v1v2[3], v2v3[3];

  /* Look for the connected vertices and its associated edge */

  cs_lnum_t  v2v_s = edges->vtx_idx[vid2];
  cs_lnum_t  v2v_e = edges->vtx_idx[vid2+1];
  cs_lnum_t  n_connect_vertices = v2v_e - v2v_s;

  cs_lnum_t  *f2f_connect = NULL;   /* To be implemented ... */

  const cs_join_vertex_t  *vertices = work->vertices;
  const double  min_limit_cos = -1.1, max_limit_cos = 1.1;
  const double  eps_dot_prod = 1e-8;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
#define _DBGTST 1
  bool  tst_dbg = (fid==492 || fid==1039 || fid==744 || fid==262 ||
                   fid==546 || fid==1057 || fid==564 ? true : false);
#else
#define _DBGTST 0
#endif

  if (f2f_connect != NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("  face splitting with face -> face"
                " connectivity is not yet implemented\n"));

  for (j = 0; j < 3; j++)
    fnorm[j] = face_normal[3*fid+j];

  *next_edge = 0;
  *next_vertex = 0; /* To enable a check at the end */

  if (n_connect_vertices > 2) { /* Look for the edge which is
                                   the most on the left */

    cs_lnum_t  left_next_edge = -1, left_next_vertex = -1;
    cs_lnum_t  right_next_edge = -1, right_next_vertex = -1;
    cs_real_t  left_min_cos = max_limit_cos;
    cs_real_t  right_max_cos = min_limit_cos;

    for (k = 0; k < 3; k++)
      v1v2[k] = vertices[vid2].coord[k]- vertices[vid1].coord[k];
    norm = _norm(v1v2);
    for (k = 0; k < 3; k++)
      v1v2[k] /= norm;

    /* Loop on connected vertices */

    for (i = v2v_s; i < v2v_e; i++) {

      cs_lnum_t  vid3 = edges->adj_vtx_lst[i];

#if _DBGTST && defined(DEBUG) && !defined(NDEBUG)
      if (tst_dbg && cs_glob_join_log != NULL)
        fprintf(cs_glob_join_log,
                "  Next vtx to test: << %9llu >> (v2v_idx: %d)\n",
                (unsigned long long)vertices[vid3].gnum, i);
#endif

      if (vid3 != vid1) {

        bool  is_in_bbox = true;
        cs_lnum_t  connect_eid = CS_ABS(edges->edge_lst[i]) - 1;

        /* Test if the connected vertex is inside the face */

        for (k = 0; k < 3; k++)
          if (   vertices[vid3].coord[k] < min_coord[k]
              || vertices[vid3].coord[k] > max_coord[k])
            is_in_bbox = false;

        if (is_in_bbox == true) {

          bool  respect_plane = _find_best_adj_face(param,
                                                    connect_eid, fid,
                                                    fnorm, adj_fnorm,
                                                    face_normal,
                                                    e2f_idx, e2f_lst);

          if (respect_plane == true) {

            /* Continue to build the new sub-face connectivity.
               adj_edge_id is in the plane as edge_id

               We look for the edge which is the most on the
               left among all the adjacent edges.

                                        V3   |
                                         \   |   /
                       adj. edge : E2 =>  \  |  /
                                           \ | /      Left part
                                            \|/
               current edge : E1 => V1------V2------  ............
                                            /|\
                                           / | \      Right part
                                          /  |  \
                                         /   |   \

               E2 is in the left part if :
                 - angle A between E1 and E2 is in [0, Pi[
                 - i.e. sin(A) >= 0
                    ->   ->  ->            ->
                 - (E1 X E2).n  >= 0 where n is the face normal vector

                 Edge which is the most on the left is the one which has
                 the minimal cos(A)

                 If there is no edge in the left part, we choose the edge
                 in the right part with the biggest cos(A) */

            double  cprod[3], mean_normal[3], cosine;

            for (k = 0; k < 3; k++) {
              v2v3[k] = vertices[vid3].coord[k] - vertices[vid2].coord[k];
              mean_normal[k] = 0.5 * (fnorm[k] + adj_fnorm[k]);
            }

            norm = _norm(v2v3);
            for (k = 0; k < 3; k++)
              v2v3[k] /= norm;

            _cross_product(v1v2, v2v3, cprod);
            dprod = _dot_product(cprod, mean_normal);
            cosine = _cosine(v1v2, v2v3);

            if (dprod >= -eps_dot_prod * _norm(cprod)) {

#if _DBGTST && defined(DEBUG) && !defined(NDEBUG)
              if (tst_dbg && cs_glob_join_log != NULL)
                fprintf(cs_glob_join_log,
                        "\tLeft choice >> cos: %10.8e / min_cos: %10.8e\n"
                        "\t         and dprod: %10.8e / crit: %10.8e\n",
                        cosine, left_min_cos, dprod,
                        -eps_dot_prod*_norm(cprod));
#endif
              /* Left part. We choose the edge with the smallest cosine */
              if (cosine < left_min_cos) {
                left_min_cos = cosine;
                left_next_vertex = vid3 + 1;
                left_next_edge = edges->edge_lst[i];
              }

            }
            else { /* In the right part. We choose the edge with
                      the biggest cosine. */

#if _DBGTST && defined(DEBUG) && !defined(NDEBUG)
              if (tst_dbg && cs_glob_join_log != NULL)
                fprintf(cs_glob_join_log,
                        "\tRight choice >> cos: %10.8e / max_cos: %10.8e\n"
                        "\t          and dprod: %10.8e / crit: %10.8e\n",
                        cosine, right_max_cos, dprod,
                        -eps_dot_prod*_norm(cprod));
#endif

              if (cosine > right_max_cos) {
                right_max_cos = cosine;
                right_next_vertex = vid3 + 1;
                right_next_edge = edges->edge_lst[i];
              }

            } /* End if dot_prod < 0 */

          } /* End if respect_plane = true */

        } /* The connected vertex is inside the face bounding box */

      } /* vid3 != vid1 */

    } /* End of loop on connected vertices */

    if (left_min_cos < max_limit_cos) {
      *next_edge = left_next_edge;
      *next_vertex = left_next_vertex;
    }
    else if (right_max_cos > min_limit_cos) {
      *next_edge = right_next_edge;
      *next_vertex = right_next_vertex;
    }
    else if (   left_min_cos  >= max_limit_cos
             && right_max_cos <= min_limit_cos)
      return OPEN_CYCLE_ERROR;

#if _DBGTST && defined(DEBUG) && !defined(NDEBUG)
    if (tst_dbg && cs_glob_join_log != NULL)
      fprintf(cs_glob_join_log, " [Result] >> next_vtx: %llu; next_edge: %llu\n",
              (unsigned long long)vertices[*next_vertex-1].gnum,
              (unsigned long long)edges->gnum[CS_ABS(*next_edge)-1]);
#endif

  } /* End if n_connect_vertices > 2 */

  else if (n_connect_vertices == 2) {

    /* Loop on connected vertices */

    for (i = v2v_s; i < v2v_e; i++) {

      cs_lnum_t  vid3 = edges->adj_vtx_lst[i];

      if (vid3 != vid1) {
        *next_edge = edges->edge_lst[i];
        *next_vertex = vid3 + 1;
      }

    } /* End of loop on connected vertices */

#if _DBGTST && defined(DEBUG) && !defined(NDEBUG)
    if (tst_dbg && cs_glob_join_log != NULL)
      fprintf(cs_glob_join_log,
              " [Result] >> next_vtx: %llu; next_edge: %llu (ONLY 2)\n",
              (unsigned long long)vertices[*next_vertex-1].gnum,
              (unsigned long long)edges->gnum[CS_ABS(*next_edge)-1]);
#endif

  }
  else { /* No connection */

    assert(n_connect_vertices < 2);

    bft_error(__FILE__, __LINE__, 0,
              _(" Joining operation : split face %d\n"
                " Problem in the connectivity. Could not find a "
                "connection with the vertex %d\n"), fid, vid1+1);


  } /* End of test on the number of vertices connected to vid2 */

  assert(*next_edge != 0);
  assert(*next_vertex != 0);

  return NO_SPLIT_ERROR;

}

/*----------------------------------------------------------------------------
 * Split the current face into sub-faces under the "plane" tolerance
 * (check if two faces are coplanear).
 *
 * parameters:
 *   fid           <-- fid of the current in the cs_join_mesh_t struct.
 *   block_id      <-- current id in face_builder_t structure
 *   param         <-- set of user-defined parameters
 *   face_normal   <-- normal for each face of the mesh
 *   work          <-- cs_join_mesh_t structure
 *   e2f_idx       <-- "edge -> face" connect. index
 *   e2f_lst       <-- "edge -> face" connect. list
 *   builder       <-> face_builder structure
 *   head_edges    <-> pointer to a resizable set structure
 *   subface_edges <-> pointer to a resizable set structure
 *   ext_edges     <-> pointer to a resizable set structure
 *   int_edges     <-> pointer to a resizable set structure
 *
 * returns:
 *   an error code (enum: cs_join_split_error_t)
 *---------------------------------------------------------------------------*/

static cs_join_split_error_t
_split_face(cs_lnum_t               fid,
            cs_lnum_t               block_id,
            cs_join_param_t         param,
            const cs_real_t         face_normal[],
            const cs_join_mesh_t   *work,
            const cs_join_edges_t  *edges,
            const cs_lnum_t        *e2f_idx,
            const cs_lnum_t        *e2f_lst,
            face_builder_t         *builder,
            cs_join_rset_t        **head_edges,
            cs_join_rset_t        **subface_edges,
            cs_join_rset_t        **ext_edges,
            cs_join_rset_t        **int_edges)
{
  cs_lnum_t  j, k, i1, i2, i_int, i_ext;
  cs_lnum_t  first_vid, vid1, vid2;
  cs_lnum_t  subface_shift, connect_shift, connect_start;
  cs_lnum_t  next_vertex, next_edge;
  cs_join_split_error_t  status;

  cs_lnum_t  n_subfaces = 0, head_edge_shift = 0;
  cs_real_t  max_coord[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
  cs_real_t  min_coord[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
  cs_join_rset_t  *_head_edges = *head_edges;
  cs_join_rset_t  *_subface_edges = *subface_edges;
  cs_join_rset_t  *_ext_edges = *ext_edges;
  cs_join_rset_t  *_int_edges = *int_edges;

  const cs_gnum_t  *fgnum = work->face_gnum;
  const int  max_subfaces = param.max_sub_faces;
  const int  verbosity = param.verbosity;
  FILE  *logfile = cs_glob_join_log;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
#define _DBGTST 1
  bool  tst_dbg = (fid==492 || fid==1039 || fid==744 || fid==262 ||
                   fid==546 || fid==1057 || fid==564 ? true : false);
  const cs_join_vertex_t  *vertices = work->vertices;
#else
#define _DBGTST 0
#endif

  /* Bounding box:  min./max. coordinates of vertices allowed for
     reconstruction */

  _face_bbox(work, fid, min_coord, max_coord);

  /* Loop on head edges */

  subface_shift = builder->face_index[block_id];
  connect_shift = builder->subface_index->array[subface_shift];
  connect_start = connect_shift;

  assert(connect_start == builder->subface_connect->n_elts);

  while (_head_edges->n_elts > 0) { /* Loop until there is one head edge */

    _int_edges->n_elts = 0;
    head_edge_shift = 0;
    _subface_edges->n_elts = 0;

    assert(_head_edges->n_elts == _ext_edges->n_elts);

    while (head_edge_shift < _head_edges->n_elts) { /* Build a new sub-face */

      cs_lnum_t  head_edge_num = _head_edges->array[head_edge_shift];
      cs_lnum_t  edge_num = head_edge_num;
      cs_lnum_t  edge_id = CS_ABS(edge_num) - 1;

#if _DBGTST && defined(DEBUG) && !defined(NDEBUG)
      if (tst_dbg && logfile != NULL)
        fprintf(logfile,
                " fnum: %d, fgnum: %llu, head_shift: %d, edge_num: %d\n",
                fid+1, (unsigned long long)fgnum[fid],
                head_edge_shift, edge_num);
#endif

      if (edge_num > 0) {
        vid1 = edges->def[2*edge_id] - 1;
        vid2 = edges->def[2*edge_id+1] - 1;
      }
      else {
        vid1 = edges->def[2*edge_id+1] - 1;
        vid2 = edges->def[2*edge_id] - 1;
      }

      cs_join_rset_resize(&(builder->subface_connect), connect_shift + 1);
      cs_join_rset_resize(&_subface_edges, _subface_edges->n_elts);

      builder->subface_connect->array[connect_shift++] = vid1+1;
      builder->subface_connect->array[connect_shift++] = vid2+1;

      _subface_edges->array[_subface_edges->n_elts++] = edge_num;

      first_vid = vid1;

#if _DBGTST && defined(DEBUG) && !defined(NDEBUG)
      if (tst_dbg && cs_glob_join_log != NULL)
        fprintf(cs_glob_join_log, " Current edge (v1,v2): [%llu,%llu]\n",
                (unsigned long long)vertices[vid1].gnum,
                (unsigned long long)vertices[vid2].gnum);
#endif

      while (vid2 != first_vid) {

        status = _find_next(param,
                            fid,
                            vid1, vid2,
                            max_coord, min_coord,
                            face_normal,
                            work,
                            edges,
                            e2f_idx, e2f_lst,
                            &next_edge, &next_vertex);

        if (status == OPEN_CYCLE_ERROR) { /* Set return pointers */

          *head_edges = _head_edges;
          *ext_edges = _ext_edges;
          *int_edges = _int_edges;
          *subface_edges = _subface_edges;

          if (verbosity > 3)
            fprintf(logfile, " Warning: open cycle for global face %llu\n",
                    (unsigned long long)fgnum[fid]);

          return OPEN_CYCLE_ERROR; /* open cycle */

        }

        assert(status == NO_SPLIT_ERROR);

        /* Add the next edge in the sub-face definition */

        cs_join_rset_resize(&_subface_edges, _subface_edges->n_elts);
        _subface_edges->array[_subface_edges->n_elts++] = next_edge;

        if (next_vertex != first_vid + 1) {

          /* Add the next vertex in the sub-face definition */
          cs_join_rset_resize(&(builder->subface_connect), connect_shift);
          builder->subface_connect->array[connect_shift++] = next_vertex;

        }

        /* To avoid an infinite loop. We check that the next_edge is not
           already in the sub-face definition */

        for (i1 = 0; i1 < _subface_edges->n_elts - 1; i1++) {

          cs_lnum_t e1 = CS_ABS(_subface_edges->array[i1]);

          for (i2 = i1 + 1; i2 < _subface_edges->n_elts; i2++) {

            cs_lnum_t e2 = CS_ABS(_subface_edges->array[i2]);

            if (e1 == e2) { /* Returns pointers */

              *head_edges = _head_edges;
              *ext_edges = _ext_edges;
              *int_edges = _int_edges;
              *subface_edges = _subface_edges;

              if (verbosity > 3)
                fprintf(logfile, " Warning: global face %llu scanned twice\n",
                        (unsigned long long)fgnum[fid]);

              return EDGE_TRAVERSED_TWICE_ERROR;  /* Face building problem */

            }

          }
        } /* End of loop on sub-face edges to check no-redundancy */

        /* Clean _head_edges if next_edge belongs to this list */

        for (j = 0; j < _head_edges->n_elts; j++)
          if (CS_ABS(next_edge) == CS_ABS(_head_edges->array[j]))
            break;

        if (j != _head_edges->n_elts) {
          _head_edges->n_elts -= 1;
          for (k = j; k < _head_edges->n_elts; k++)
            _head_edges->array[k] = _head_edges->array[k+1];
        }

        /* Update the couple (vid1, vid2) to continue the sub-face building */

        vid1 = vid2;
        vid2 = next_vertex - 1;

        /* Test if next_edge is in _ext_edges */

        for (i_ext = 0; i_ext < _ext_edges->n_elts; i_ext++)
          if (CS_ABS(next_edge) == CS_ABS(_ext_edges->array[i_ext]))
            break;

        if (i_ext != 0 && i_ext == _ext_edges->n_elts) {

          /* next_edge is not in _ext_edges
             Test if next_edge is in the _int_edges list.
               If not : store next_edge in the _int_edges list.
               If yes : delete it. */

          for (i_int = 0; i_int < _int_edges->n_elts; i_int++)
            if (CS_ABS(next_edge) == CS_ABS(_int_edges->array[i_int]))
              break;

          if (i_int == _int_edges->n_elts) { /* Add next_edge to the list */
            cs_join_rset_resize(&_int_edges, _int_edges->n_elts);
            _int_edges->array[_int_edges->n_elts++] = next_edge;
          }
          else { /* Delete next_edge of the list */

            _int_edges->n_elts -= 1;
            for (k = i_int; k < _int_edges->n_elts; k++)
              _int_edges->array[k] = _int_edges->array[k+1];

          }

        } /* Next_edge is not an ext_edges */

      } /* End of while vid2 != first_vid */

      /* End of building the sub-face */

      builder->subface_connect->n_elts = connect_shift;
      connect_start = connect_shift;
      _subface_edges->n_elts = 0;

      subface_shift++;
      builder->subface_index->n_elts = subface_shift;
      cs_join_rset_resize(&(builder->subface_index), subface_shift);
      builder->subface_index->array[subface_shift] = connect_start;

      head_edge_shift++;
      n_subfaces++;

#if _DBGTST && defined(DEBUG) && !defined(NDEBUG)
      if (tst_dbg && cs_glob_join_log != NULL)
        fprintf(cs_glob_join_log,
                " END OF BUILDING subface %d\n\n", n_subfaces);
#endif

      if (n_subfaces > max_subfaces) { /* Set return pointers */

        *head_edges = _head_edges;
        *ext_edges = _ext_edges;
        *int_edges = _int_edges;
        *subface_edges = _subface_edges;

        if (verbosity > 3)
          fprintf(logfile,
                  " Warning: loop limit reached for global face %llu\n",
                  (unsigned long long)fgnum[fid]);

        return LOOP_LIMIT_ERROR;
      }

    } /* End of while head_shift < head_edge->n_elts */

    /* Build a _head_edges set from _int_edges. We have to invert
       edge direction to be consistent in the sub-face definition */

    cs_join_rset_resize(&_head_edges, _int_edges->n_elts);
    cs_join_rset_resize(&_ext_edges, _int_edges->n_elts);

    _head_edges->n_elts = 0;
    _ext_edges->n_elts = 0;

    for (i_int = 0; i_int < _int_edges->n_elts; i_int++) {
      _head_edges->array[_head_edges->n_elts++] = -_int_edges->array[i_int];
      _ext_edges->array[_ext_edges->n_elts++] = -_int_edges->array[i_int];
    }

  } /* End of loop on head edges */

  builder->face_index[block_id+1] = builder->face_index[block_id] + n_subfaces;

  assert(builder->face_index[block_id+1] = subface_shift);

  /* Returns pointers */

  *head_edges = _head_edges;
  *ext_edges = _ext_edges;
  *int_edges = _int_edges;
  *subface_edges = _subface_edges;

  return NO_SPLIT_ERROR;
}

/*----------------------------------------------------------------------------
 * Compare two elements in an indexed list and returns true if element in
 * position i1 is strictly greater than element in position i2.
 *
 * parameters:
 *   i1     <-- position in index for the first element
 *   i2     <-- position in index for the second element
 *   index  <-- number of values to compare for each entity
 *   number <-- pointer to numbers of entities that should be ordered.
 *              (if NULL, a default 1 to n numbering is considered)
 *
 * returns:
 *   true or false
 *---------------------------------------------------------------------------*/

inline static bool
_indexed_is_greater(size_t           i1,
                    size_t           i2,
                    const cs_lnum_t  index[],
                    const cs_gnum_t  number[])
{
  int  i;

  cs_lnum_t  i1_s = index[i1], i1_e = index[i1+1], s1 = i1_e - i1_s;
  cs_lnum_t  i2_s = index[i2], i2_e = index[i2+1], s2 = i2_e - i2_s;

  if (s1 > s2) {

    for (i = 0; i < s2; i++) {
      if (number[i1_s + i] > number[i2_s + i])
        return true;
      else if (number[i1_s + i] < number[i2_s + i])
        return false;
    }

    return true;
  }
  else { /* s1 <= s2 */

    for (i = 0; i < s1; i++) {
      if (number[i1_s + i] > number[i2_s + i])
        return true;
      else if (number[i1_s + i] < number[i2_s + i])
        return false;
    }

    return false;
  }
}

/*----------------------------------------------------------------------------
 * Define a global numbering for the subfaces after splitting faces.
 * Synchronize face connectivity.
 *
 * parameters:
 *   builder  <-> pointer to a face builder structure
 *   work     <-- cs_join_mesh_t structure
 *---------------------------------------------------------------------------*/

static void
_get_subface_gnum(face_builder_t         *builder,
                  const cs_join_mesh_t   *work)
{
  cs_lnum_t  i, j, k, shift;
  cs_gnum_t  min_val;

  cs_lnum_t  max_size = 0;
  cs_lnum_t  n_subfaces = builder->face_index[builder->n_faces];
  cs_lnum_t  *index = builder->subface_index->array;
  cs_gnum_t  *gconnect = builder->subface_gconnect;
  cs_gnum_t  *glob_list = NULL, *tmp = NULL, *vgnum = NULL;
  fvm_io_num_t *subface_io_num = NULL;

  const cs_gnum_t  *global_num = NULL;
  const cs_join_vertex_t  *vertices = work->vertices;

  assert(index != NULL);

  /* Allocate the buffer we want to define */

  BFT_MALLOC(builder->subface_gnum, n_subfaces, cs_gnum_t);

  /* Re-arrange gconnect in order to have for each subface:
      - vertex with the minimal glob. num. in first place,
      - vertex with the minimal glob. num. between the two possible one
        in second place
  */

  for (i = 0; i < n_subfaces; i++)
    max_size = CS_MAX(max_size, index[i+1] - index[i]);

  BFT_MALLOC(tmp, max_size, cs_gnum_t);
  BFT_MALLOC(glob_list, index[n_subfaces], cs_gnum_t);

  /* Build glob_list */

  for (i = 0; i < n_subfaces; i++) {

    cs_lnum_t  start = index[i], end = index[i+1];
    cs_lnum_t  n_elts = end - start;

    assert(n_elts > 1);

    for (j = start; j < end; j++)
      glob_list[j] = gconnect[j];

    min_val = glob_list[start];
    k = 0;

    for (j = start + 1; j < end; j++) {
      if (glob_list[j] < min_val) {
        min_val = glob_list[j];
        k = j - start;
      }
    }

    k = (n_elts - k) % n_elts; /* Define the permutation factor */

    /* Permutation in order to the kth element in first position */

    for (j = 0; j < n_elts; j++) {
      shift = (j + k) % n_elts;
      tmp[shift] = glob_list[start + j];
    }

    assert(tmp[0] == min_val);
    glob_list[start] = tmp[0];

    if (tmp[1] > tmp[n_elts-1]) { /* Inverse order */
      for (j = 1; j < n_elts; j++)
        glob_list[start + j] = tmp[n_elts - j];
    }
    else { /* Keep the current order */
      for (j = 1; j < n_elts; j++)
        glob_list[start + j] = tmp[j];
    }

  } /* End of loop on subfaces */

  /* Copy glob_list as the new subface global connectivity and
     use it to define a synchronized sub-face connectivity */

  BFT_MALLOC(vgnum, work->n_vertices, cs_gnum_t);

  for (i = 0; i < work->n_vertices; i++)
    vgnum[i] = vertices[i].gnum;

  for (i = 0; i < n_subfaces; i++) {

    for (j = index[i]; j < index[i+1]; j++) {

      gconnect[j] = glob_list[j];

      k = cs_search_g_binary(work->n_vertices,
                             glob_list[j],
                             vgnum);

      assert(k != -1);
      builder->subface_connect->array[j] = k+1;

    }

  } /* End of loop on subfaces */

  BFT_FREE(vgnum);

  if (cs_glob_n_ranks > 1) { /* Parallel treatment */

    /* Local ordering */

    cs_lnum_t  *order_index = NULL;
    cs_gnum_t  *order_glob_list = NULL;
    cs_lnum_t  *order = cs_order_gnum_i(NULL,
                                        glob_list,
                                        index,
                                        n_subfaces);

    _renumber_local_ordered_i(n_subfaces,
                              order,
                              index,
                              glob_list,
                              &order_index,
                              &order_glob_list);

    subface_io_num = fvm_io_num_create_from_adj_i(NULL,
                                                  order_index,
                                                  order_glob_list,
                                                  n_subfaces);

    assert(subface_io_num != NULL);
    assert(n_subfaces == fvm_io_num_get_local_count(subface_io_num));

    global_num = fvm_io_num_get_global_num(subface_io_num);

    for (i = 0; i < n_subfaces; i++)
      builder->subface_gnum[order[i]] = global_num[i];

    builder->n_g_subfaces = fvm_io_num_get_global_count(subface_io_num);

    BFT_FREE(order);
    BFT_FREE(order_index);
    BFT_FREE(order_glob_list);

    subface_io_num = fvm_io_num_destroy(subface_io_num);

  }
  else { /* Serial treatment */

    cs_gnum_t  gnum = 1;
    cs_lnum_t  *order = cs_order_gnum_i(NULL,
                                        glob_list,
                                        index,
                                        n_subfaces);

    builder->subface_gnum[order[0]] = gnum;

    for (i = 1; i < n_subfaces; i++) {
      if (_indexed_is_greater(order[i], order[i-1], index, glob_list)) gnum++;
      builder->subface_gnum[order[i]] = gnum;
    }
    builder->n_g_subfaces = gnum;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    if (cs_glob_join_log != NULL) {
      for (i = 0; i < n_subfaces; i++) {
        cs_lnum_t  start = index[i], end = index[i+1];
        cs_lnum_t  n_elts = end - start;
        fprintf(cs_glob_join_log,
                " subface %5d - gnum: %10llu - connect_size: %d - ",
                i+1, (unsigned long long)builder->subface_gnum[i], n_elts);
        for (j = start; j < end; j++)
          fprintf(cs_glob_join_log, " %8llu ",
                  (unsigned long long)glob_list[j]);
        fprintf(cs_glob_join_log, "\n");
      }
      fflush(cs_glob_join_log);
    }
#endif

    BFT_FREE(order);

  }

  if (cs_glob_mesh->verbosity > 0) {
    bft_printf(_("\n  Global number of faces after splitting: %10llu\n"),
               (unsigned long long)builder->n_g_subfaces);
    bft_printf_flush();
  }

  /* Free memory */

  BFT_FREE(tmp);
  BFT_FREE(glob_list);

}

/*----------------------------------------------------------------------------
 * Update a cs_join_mesh_t structure thanks to a face_builder_t structure.
 *
 * parameters:
 *   block_info  <-- set of paramaters defining a contiguous distribution
 *   builder     <-- pointer to the distributed face builder structure
 *   mesh        <-> pointer to the local cs_join_mesh_t structure
 *   p_o2n_hist  <-> pointer to old global face -> new local face numbering
 *---------------------------------------------------------------------------*/

static void
_update_mesh_after_split(cs_block_dist_info_t    bi,
                         face_builder_t         *builder,
                         cs_join_mesh_t        **mesh,
                         cs_join_gset_t        **p_o2n_hist)
{
  cs_lnum_t  i, j, k, id, shift, n_subfaces, o_id;
  cs_gnum_t  prev, cur;

  cs_lnum_t  n_new_faces = 0, block_size = 0;
  char  *new_mesh_name = NULL;
  cs_lnum_t  *subfaces = NULL;
  cs_lnum_t  *order = NULL;
  cs_join_gset_t  *o2n_hist = NULL;
  cs_join_mesh_t  *init_mesh = *mesh, *new_mesh = NULL;

  if (bi.gnum_range[1] > bi.gnum_range[0])
    block_size = bi.gnum_range[1] - bi.gnum_range[0];

  /* Sanity checks */

  assert(init_mesh != NULL);
  assert(builder != NULL);

  /* Create a new cs_join_mesh_t structure */

  BFT_MALLOC(new_mesh_name, strlen("AfterSplitting_n") + 5 + 1, char);
  sprintf(new_mesh_name,"%s%05d", "AfterSplitting_n",
          CS_MAX(cs_glob_rank_id, 0));

  new_mesh = cs_join_mesh_create(new_mesh_name);

  BFT_FREE(new_mesh_name);

  if (block_size != builder->n_faces)
    bft_error(__FILE__, __LINE__, 0,
              " Inconsistent values between:\n"
              "    block_size             %8ld\n"
              "    builder->n_faces:      %8ld\n"
              " These values should be equal.\n",
              (long)block_size, (long)builder->n_faces);

  /* Compute the number of new faces */

  n_subfaces = builder->face_index[builder->n_faces];

  BFT_MALLOC(order, n_subfaces, cs_lnum_t);

  cs_order_gnum_allocated(NULL, builder->subface_gnum, order, n_subfaces);

  prev = 0;
  n_new_faces = 0;

  for (i = 0; i < n_subfaces; i++) {

    cur = builder->subface_gnum[order[i]];
    if (prev != cur) {
      prev = cur;
      n_new_faces++;
    }

  } /* End of loop on subfaces */

  /* Build new cell_gnum array */

  BFT_MALLOC(subfaces, n_new_faces, cs_lnum_t);

  prev = 0;
  n_new_faces = 0;

  for (i = 0; i < n_subfaces; i++) {

    o_id = order[i];
    cur = builder->subface_gnum[o_id];

    if (prev != cur) {
      prev = cur;
      subfaces[n_new_faces] = o_id;
      n_new_faces++;
    }

  }

  /* Define the new face connectivities */

  new_mesh->n_faces = n_new_faces;
  new_mesh->n_g_faces = builder->n_g_subfaces;

  BFT_MALLOC(new_mesh->face_gnum, n_new_faces, cs_gnum_t);
  BFT_MALLOC(new_mesh->face_vtx_idx, n_new_faces + 1, cs_lnum_t);

  for (i = 0; i < n_new_faces; i++) {

      id = subfaces[i];
      new_mesh->face_gnum[i] = builder->subface_gnum[id];
      new_mesh->face_vtx_idx[i+1] =  builder->subface_index->array[id+1]
                                   - builder->subface_index->array[id];

  } /* End of loop on new faces */

  new_mesh->face_vtx_idx[0] = 0;
  for (i = 0; i < n_new_faces; i++)
    new_mesh->face_vtx_idx[i+1] += new_mesh->face_vtx_idx[i];

  BFT_MALLOC(new_mesh->face_vtx_lst,
             new_mesh->face_vtx_idx[n_new_faces],
             cs_lnum_t);

  for (i = 0; i < n_new_faces; i++) {

    id = subfaces[i];
    shift = new_mesh->face_vtx_idx[i];

    for (j = builder->subface_index->array[id], k = 0;
         j < builder->subface_index->array[id+1]; j++, k++)
      new_mesh->face_vtx_lst[shift + k] = builder->subface_connect->array[j] - 1;

  } /* End of loop on new faces */

  /* Copy vertex data */

  new_mesh->n_g_vertices = init_mesh->n_g_vertices;
  new_mesh->n_vertices = init_mesh->n_vertices;

  BFT_MALLOC(new_mesh->vertices, new_mesh->n_vertices, cs_join_vertex_t);

  for (i = 0; i < init_mesh->n_vertices; i++)
    new_mesh->vertices[i] = init_mesh->vertices[i];

  cs_join_mesh_vertex_clean(new_mesh);

  /* Free memory */

  BFT_FREE(subfaces);
  BFT_FREE(order);

  /* Create a structure in which we keep a history of global
     face numbering for each new face */

  o2n_hist = cs_join_gset_create(block_size);

  if (block_size > 0) {

    assert(builder != NULL);
    assert(builder->n_faces == block_size);

    /* Historic is a part of the data held in builder structure */

    /* store old glob. face num. */
    for (i = 0; i < block_size; i++) {
      cs_gnum_t g_id = i;
      o2n_hist->g_elts[i] = bi.gnum_range[0] + g_id;
    }

    for (i = 0; i < builder->n_faces + 1; i++)
      o2n_hist->index[i] = builder->face_index[i];

    BFT_MALLOC(o2n_hist->g_list,
               o2n_hist->index[o2n_hist->n_elts], cs_gnum_t);

    for (i = 0; i < builder->n_faces; i++) {

      for (j = builder->face_index[i]; j < builder->face_index[i+1]; j++) {

        id = cs_search_g_binary(new_mesh->n_faces,
                                builder->subface_gnum[j],
                                new_mesh->face_gnum);

        assert(id != -1);
        o2n_hist->g_list[j] = id + 1;  /* store local new face num. */

      }

    } /* End of loop on initial faces */

  } /* block.local_size > 0 */

  cs_join_mesh_destroy(&init_mesh);

  /* Set return pointer */

  *mesh = new_mesh;
  *p_o2n_hist = o2n_hist;

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Build new faces after the vertex merge operation. Split initial faces into
 * subfaces and keep the historic between initial/final faces.
 *
 * parameters:
 *   param           <-- set of user-defined parameters
 *   face_normal     <-- array of normal vector on each face
 *   edges           <-- list of edges
 *   work            <-> pointer to a cs_join_mesh_t structure
 *   old2new_history <-> relation between new faces and old one:
 *                       old global face -> new local face
 *---------------------------------------------------------------------------*/

void
cs_join_split_faces(cs_join_param_t          param,
                    const cs_real_t          face_normal[],
                    const cs_join_edges_t   *edges,
                    cs_join_mesh_t         **work,
                    cs_join_gset_t         **old2new_history)
{
  cs_lnum_t  fid, j, face_s, subface_s, block_id, vid;
  cs_gnum_t  vgnum;
  cs_join_split_error_t  code;
  cs_block_dist_info_t  bi;

  cs_lnum_t  _n_problems = 0, n_face_problems = 0, n_max_face_vertices = 6;
  cs_lnum_t  block_size = 0;
  cs_lnum_t  *e2f_idx = NULL, *e2f_lst = NULL;
  cs_join_gset_t  *_old2new_history = NULL;
  cs_join_rset_t  *open_cycle = NULL, *edge_traversed_twice = NULL;
  cs_join_rset_t  *loop_limit = NULL, *head_edges = NULL;
  cs_join_rset_t  *subface_edges = NULL, *ext_edges = NULL, *int_edges = NULL;
  face_builder_t  *builder = NULL;
  cs_join_mesh_t  *w = *work;

  const cs_lnum_t  n_init_faces = w->n_faces;
  const int  n_ranks = cs_glob_n_ranks;
  const int  local_rank = CS_MAX(cs_glob_rank_id, 0);

  assert(w != NULL);
  assert(edges != NULL);

  /* Use the cs_join_edges_t structure to build
     the "edge -> face" connectivity */

  cs_join_mesh_get_edge_face_adj(w, edges, &e2f_idx, &e2f_lst);

  /* Define buffers to manage errors */

  open_cycle = cs_join_rset_create(3);
  edge_traversed_twice = cs_join_rset_create(3);
  loop_limit = cs_join_rset_create(3);

  /* Define buffers and structures to build the new faces */

  head_edges = cs_join_rset_create(n_max_face_vertices);
  subface_edges = cs_join_rset_create(n_max_face_vertices);
  ext_edges = cs_join_rset_create(n_max_face_vertices);
  int_edges = cs_join_rset_create(n_max_face_vertices);

  /* Compute block_size */

  bi = cs_block_dist_compute_sizes(local_rank,
                                   n_ranks,
                                   1,
                                   0,
                                   w->n_g_faces);

  if (bi.gnum_range[1] > bi.gnum_range[0])
    block_size = bi.gnum_range[1] - bi.gnum_range[0];

  builder = _create_face_builder(block_size);

  /*
     We only have to treat faces for the current rank's block because the
     initial face distribution assumes that the current rank can correctly
     split only faces in its block.
     Main loop on faces.
  */

  for (fid = 0, block_id = 0; fid < n_init_faces; fid++) {

    int  block_rank = (w->face_gnum[fid] - 1)/(cs_gnum_t)(bi.block_size);

    if (block_rank == local_rank) { /* This face is a "main" face for the
                                       local rank */

      int  n_face_vertices = w->face_vtx_idx[fid+1] - w->face_vtx_idx[fid];

      if (n_face_vertices > n_max_face_vertices) { /* Manage list size */

        n_max_face_vertices = n_face_vertices;
        cs_join_rset_resize(&head_edges, n_face_vertices);
        cs_join_rset_resize(&subface_edges, n_face_vertices);
        cs_join_rset_resize(&ext_edges, n_face_vertices);
        cs_join_rset_resize(&int_edges, n_face_vertices);

      }

      /* Fill head_edges and ext_edges */

      _define_head_and_ext_edges(fid,
                                 w, edges,
                                 head_edges, ext_edges,
                                 0); /* No permutation */

      /* Store initial builder state in case of code > 0 */

      face_s = builder->face_index[block_id];
      subface_s = builder->subface_index->array[face_s];

      /* Split the current face into subfaces */

      code = _split_face(fid, block_id,
                         param,
                         face_normal,
                         w, edges,
                         e2f_idx, e2f_lst,
                         builder,
                         &head_edges, &subface_edges,
                         &ext_edges, &int_edges);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
      if (param.verbosity > 2 && code != NO_SPLIT_ERROR) {
        fprintf(cs_glob_join_log,
                "  Split face %d -> returned code: %d\n", fid+1, code);
        _dump_face_builder(block_id, builder, cs_glob_join_log);
      }
#endif

      if (code > NO_SPLIT_ERROR) { /* Manage error */

        _n_problems++;

        /* We change the starting edge for traversing edges to build
           the new face. This may be enough to solve the problem when
           the face is warped. */

        while (_n_problems < n_face_vertices && code > NO_SPLIT_ERROR) {

          /* Fill head_edges and ext_edges */

          _define_head_and_ext_edges(fid,
                                     w, edges,
                                     head_edges, ext_edges,
                                     _n_problems); /* permutation */

          /* Retrieve initial builder state */

          builder->face_index[block_id] = face_s;
          builder->subface_index->array[face_s] = subface_s;
          builder->subface_connect->n_elts = subface_s;

          /* Split the current face into subfaces */

          code = _split_face(fid, block_id,
                             param,
                             face_normal,
                             w, edges,
                             e2f_idx, e2f_lst,
                             builder,
                             &head_edges, &subface_edges,
                             &ext_edges, &int_edges);

          _n_problems++;

        } /* End of while */

        if (_n_problems >= n_face_vertices && code != NO_SPLIT_ERROR) {

          n_face_problems++;

          switch (code) {

          case OPEN_CYCLE_ERROR:
            cs_join_rset_resize(&open_cycle, open_cycle->n_elts);
            open_cycle->array[open_cycle->n_elts] = fid + 1;
            open_cycle->n_elts += 1;
            break;

          case EDGE_TRAVERSED_TWICE_ERROR:
            cs_join_rset_resize(&edge_traversed_twice,
                                edge_traversed_twice->n_elts);
            edge_traversed_twice->array[edge_traversed_twice->n_elts] = fid + 1;
            edge_traversed_twice->n_elts += 1;
            break;

          case LOOP_LIMIT_ERROR:
            cs_join_rset_resize(&loop_limit, loop_limit->n_elts);
            loop_limit->array[loop_limit->n_elts] = fid + 1;
            loop_limit->n_elts += 1;
            break;

          case NO_SPLIT_ERROR: /* To avoid a warning */
            assert(0);
            break;

          case MAX_SPLIT_ERROR: /* To avoid a warning */
            assert(0);
            break;

          }

          /* Keep the initial face connectivity */

          builder->face_index[block_id] = face_s;
          builder->face_index[block_id+1] = face_s + 1;

          /* face -> subface connectivity index update */

          cs_join_rset_resize(&(builder->subface_index), face_s+1);

          builder->subface_index->n_elts = face_s + 1;
          builder->subface_index->array[face_s] = subface_s;
          builder->subface_index->array[face_s+1] = subface_s + n_face_vertices;

          /* face -> subface connectivity list update */

          cs_join_rset_resize(&(builder->subface_connect),
                              subface_s + n_face_vertices);

          builder->subface_connect->n_elts = subface_s + n_face_vertices;

          for (j = 0; j < n_face_vertices; j++)
            builder->subface_connect->array[subface_s + j]
              = w->face_vtx_lst[j + w->face_vtx_idx[fid]] + 1;

          if (param.verbosity > 2) {
            fprintf(cs_glob_join_log,
                    "\n Keep initial connectivity for face %d [%llu]:\n",
                    fid+1, (unsigned long long)w->face_gnum[fid]);
            _dump_face_builder(block_id, builder, cs_glob_join_log);
          }

        } /* End if n_current_face_problems >= n_face_vertices */

        _n_problems = 0;

      } /* End of error managing */

      /* Ready to treat the next face in the block of the current rank */

      block_id++;

    } /* End if current face is in rank block */

  } /* End of loop on faces */

  /* Delete lists */

  cs_join_rset_destroy(&head_edges);
  cs_join_rset_destroy(&subface_edges);
  cs_join_rset_destroy(&ext_edges);
  cs_join_rset_destroy(&int_edges);

  /* Display information in case of problem during the face splitting */

  {
    cs_gnum_t  g_in_buf[4], g_out_buf[4];

    cs_gnum_t  n_g_face_problems = n_face_problems;
    cs_gnum_t  n_g_open_cycles = open_cycle->n_elts;
    cs_gnum_t  n_g_edges_twice = edge_traversed_twice->n_elts;
    cs_gnum_t  n_g_loop_limit = loop_limit->n_elts;

#if defined(HAVE_MPI)
    if (n_ranks > 1) {

      MPI_Comm  mpi_comm = cs_glob_mpi_comm;

      g_in_buf[0] = n_g_face_problems;
      g_in_buf[1] = n_g_open_cycles;
      g_in_buf[2] = n_g_edges_twice;
      g_in_buf[3] = n_g_loop_limit;

      MPI_Allreduce(g_in_buf, g_out_buf, 4, CS_MPI_GNUM, MPI_SUM, mpi_comm);

      n_g_face_problems = g_out_buf[0];
      n_g_open_cycles = g_out_buf[1];
      n_g_edges_twice = g_out_buf[2];
      n_g_loop_limit = g_out_buf[3];

    }
#endif

    if (n_g_face_problems > 0) {

      bft_printf
        (_("\n  *** WARNING ***\n"
           "  Globally, %llu problem(s) found during the face splitting\n"
           "     %12llu  open cycles,\n"
           "     %12llu  edges traversed twice,\n"
           "     %12llu  faces split into more than max_subfaces (= %d)\n\n"
           "    => Eventually modify joining parameters\n\n"),
         (unsigned long long)n_g_face_problems,
         (unsigned long long)n_g_open_cycles,
         (unsigned long long)n_g_edges_twice,
         (unsigned long long)n_g_loop_limit, param.max_sub_faces);
      bft_printf_flush();

      assert(   n_g_face_problems
             == n_g_open_cycles + n_g_edges_twice + n_g_loop_limit);

      /* post-processing of encountered problems */

      cs_join_post_init(); /* Init. post-processing anyway */

      if (n_g_open_cycles > 0)
        cs_join_post_faces_subset("OpenCycleErr",
                                  w,
                                  open_cycle->n_elts,
                                  open_cycle->array);

      if (n_g_edges_twice > 0)
        cs_join_post_faces_subset("EdgeScannedTwiceErr",
                                  w,
                                  edge_traversed_twice->n_elts,
                                  edge_traversed_twice->array);

      if (n_g_loop_limit > 0) {
        bft_printf(_("At least one original face has been cut into more than"
                     "%d subfaces\n. You can increase this parameter in"
                     "cs_user_join() or cs_user_periodicity()\n"
                     "by setting the advanced parameter max_sub_face.\n"
                     "Be cautious, as this may produce a mesh with"
                     "a poor quality.\n"), param.max_sub_faces);
        cs_join_post_faces_subset("LoopLimitErr",
                                  w,
                                  loop_limit->n_elts,
                                  loop_limit->array);
      }

    } /* End of information display */

  } /* Manage potential errors during face splitting */

  /* Free vtx_struct structure */

  BFT_FREE(e2f_idx);
  BFT_FREE(e2f_lst);

  /* Delete error management lists */

  cs_join_rset_destroy(&open_cycle);
  cs_join_rset_destroy(&edge_traversed_twice);
  cs_join_rset_destroy(&loop_limit);

  { /* Define a global number for each new sub-faces */

    int  n_subfaces = builder->face_index[builder->n_faces];
    int  sub_connect_size = builder->subface_index->array[n_subfaces];

    /* Define subface_gconnect */

    BFT_MALLOC(builder->subface_gconnect, sub_connect_size, cs_gnum_t);

    for (j = 0; j < sub_connect_size; j++) {
      vid = builder->subface_connect->array[j] - 1;
      vgnum = w->vertices[vid].gnum;
      builder->subface_gconnect[j] = vgnum;
    }

    _get_subface_gnum(builder, w);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    if (cs_glob_join_log != NULL) {
      fprintf(cs_glob_join_log, "\nFINAL BUILDER STATE\n");
      for (fid = 0; fid < builder->n_faces; fid++)
        _dump_face_builder(fid, builder, cs_glob_join_log);
    }
#endif

    BFT_FREE(builder->subface_gconnect);

  }

  /* Update cs_join_mesh_t structure and keep a relation between
     new and old faces.

     Delete redundancies in the sub-face connectivity.
     Some new interior faces can be defined twice. See the example below.

                 .-------------------------.
         First   |                         |
         initial |                         |
         face    |     .-------------------+--------.
                 |     |                   |        |
                 |     | interior sub-face |        |
                 |     | defined twice     |        |
                 |     |                   |        |
                 `-----+-------------------'        | Second
                       |                            | initial
                       |                            | face
                       `----------------------------'

     Reduce the definition of the working mesh to the set of new
     global face numbers built from the local block.
     For each new sub-face we maintain a relation between new and old
     face global number. Update also face state */

  _update_mesh_after_split(bi,
                           builder,
                           &w,
                           &_old2new_history);

  /* Free face_builder_t structure */

  builder = _destroy_face_builder(builder);

  /* Set return pointers */

  *work = w;
  *old2new_history = _old2new_history;
}

/*----------------------------------------------------------------------------
 * Update after face splitting of the local join mesh structure.
 * Send back to the original rank the new face description.
 *
 * parameters:
 *   param           <-- set of user-defined parameter
 *   work_mesh       <-- distributed mesh on faces to join
 *   gnum_rank_index <-- index on ranks for the old global face numbering
 *   o2n_hist        <-> old global face -> new local face numbering
 *   local_mesh      <-> mesh on local selected faces to be joined
 *---------------------------------------------------------------------------*/

void
cs_join_split_update_struct(const cs_join_param_t   param,
                            const cs_join_mesh_t   *work_mesh,
                            const cs_gnum_t         gnum_rank_index[],
                            cs_join_gset_t        **o2n_hist,
                            cs_join_mesh_t        **local_mesh)
{
  cs_lnum_t  i;

  cs_join_gset_t  *_o2n_hist = *o2n_hist;
  cs_join_mesh_t  *_local_mesh = *local_mesh;

  const int  n_ranks = cs_glob_n_ranks;

  /* sanity checks */

  assert(work_mesh != NULL);
  assert(_local_mesh != NULL);
  assert(_o2n_hist != NULL || (_o2n_hist == NULL && work_mesh->n_faces == 0));

  if (n_ranks == 1)
    cs_join_mesh_copy(&_local_mesh, work_mesh);

#if defined(HAVE_MPI)
  if (n_ranks > 1) { /* Parallel mode */

    cs_lnum_t  j, subface_id;
    cs_lnum_t  *send_rank_index = NULL, *send_faces = NULL;
    cs_gnum_t  *init_face_gnum = NULL;
    cs_join_gset_t  *distrib_sync_hist = NULL;
    cs_lnum_t  n_init_faces = _local_mesh->n_faces;
    cs_gnum_t  n_g_init_faces = _local_mesh->n_g_faces;

    MPI_Comm  mpi_comm = cs_glob_mpi_comm;

    if (param.perio_type != FVM_PERIODICITY_NULL) {

      n_g_init_faces *= 2;

      /* Save the initial global face numbering */

      BFT_MALLOC(init_face_gnum, 2*n_init_faces, cs_gnum_t);

      for (i = 0; i < n_init_faces; i++) {
        init_face_gnum[2*i] = _local_mesh->face_gnum[i];
        init_face_gnum[2*i+1] = _local_mesh->face_gnum[i]+1;
      }

      n_init_faces *= 2;

    }
    else {

      /* Save the initial global face numbering */

      BFT_MALLOC(init_face_gnum, n_init_faces, cs_gnum_t);

      for (i = 0; i < n_init_faces; i++)
        init_face_gnum[i] = _local_mesh->face_gnum[i];

    }

    /* Free some structures of the mesh */

    cs_join_mesh_reset(_local_mesh);

    _get_faces_to_send(_o2n_hist, /* Local subface num is stored up to now */
                       gnum_rank_index,
                       &send_rank_index,
                       &send_faces);

    /* Get the new face connectivity from the distributed work_mesh */

    cs_join_mesh_exchange(n_ranks,
                          send_rank_index,
                          send_faces,
                          work_mesh,
                          _local_mesh,
                          mpi_comm);

    BFT_FREE(send_faces);
    BFT_FREE(send_rank_index);

    /* Order face by increasing global number */

    cs_join_mesh_face_order(_local_mesh);

    /* Move to a global numbering of the subfaces */

    for (i = 0; i < _o2n_hist->n_elts; i++) {
      for (j = _o2n_hist->index[i]; j < _o2n_hist->index[i+1]; j++) {
        subface_id = _o2n_hist->g_list[j] - 1;
        _o2n_hist->g_list[j] = work_mesh->face_gnum[subface_id];
      }
    }

    /* Synchronize _o2n_hist */

    distrib_sync_hist = cs_join_gset_block_sync(n_g_init_faces,
                                                _o2n_hist,
                                                mpi_comm);

    cs_join_gset_destroy(&_o2n_hist);
    _o2n_hist = cs_join_gset_create(n_init_faces);

    for (i = 0; i < n_init_faces; i++)
      _o2n_hist->g_elts[i] = init_face_gnum[i];

    BFT_FREE(init_face_gnum);

    cs_join_gset_block_update(n_g_init_faces,
                              distrib_sync_hist,
                              _o2n_hist,
                              mpi_comm);

    cs_join_gset_destroy(&distrib_sync_hist);

  }
#endif /* HAVE_MPI */

  /* Set return pointers */

  *o2n_hist = _o2n_hist;
  *local_mesh = _local_mesh;
}

/*---------------------------------------------------------------------------*/

END_C_DECLS

