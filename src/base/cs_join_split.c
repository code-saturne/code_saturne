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
 * Split faces during the joining operation
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
 * FVM library headers
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

#include "cs_join_split.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

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

  cs_int_t         n_faces;
  cs_int_t        *face_index;       /* Face -> Subface index */
  cs_join_rset_t  *subface_index;    /* Subface -> vertex connect. index */
  cs_join_rset_t  *subface_connect;  /* Subface -> vertex connect. list */

  /* The two following arrays are built after the face splitting. So, they
     don't need to be resizable */

  fvm_gnum_t    n_g_subfaces;     /* Global number of subfaces after splitting */
  fvm_gnum_t   *subface_gconnect; /* Subface -> glob. vertex list */
  fvm_gnum_t   *subface_gnum;    /* Subface global numbering */

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
  double  result = dprod / (n1 * n2);

  return result;
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
_renumber_local_ordered_i(cs_int_t          n_elts,
                          const fvm_lnum_t  order[],
                          const cs_int_t    index[],
                          fvm_gnum_t        glist[],
                          cs_int_t         *new_index[],
                          fvm_gnum_t       *new_glist[])
{
  cs_int_t  i, j, k, o_id;

  cs_int_t  *_new_index = NULL;
  fvm_gnum_t  *_new_glist = NULL;

  if (n_elts < 1)
    return;

  assert(index[0] == 0); /* case index[0] = 1 coulb be coded in the future */

  /* Build a new index */

  BFT_MALLOC(_new_index, n_elts + 1, cs_int_t);

  for (i = 0; i < n_elts; i++) {
    o_id = order[i];
    _new_index[i+1] = index[o_id+1] - index[o_id];
  }

  _new_index[0] = 0;
  for (i = 0; i < n_elts; i++)
    _new_index[i+1] += _new_index[i];

  /* Build a new list */

  BFT_MALLOC(_new_glist, _new_index[n_elts], fvm_gnum_t);

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
                   const fvm_gnum_t       gnum_rank_index[],
                   cs_int_t              *send_rank_index[],
                   cs_int_t              *send_faces[])
{
  cs_int_t  i, j, rank, start, end;

  cs_int_t  reduce_size = 0;
  cs_int_t  *_send_rank_index = NULL, *_send_faces = NULL, *reduce_ids = NULL;
  fvm_gnum_t  *reduce_index = NULL;
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

  BFT_MALLOC(new_face_rank->g_list, new_face_rank->index[n_ranks], fvm_gnum_t);

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

      cs_int_t  shift =  new_face_rank->index[rank]
                       + new_face_rank->g_elts[rank];
      cs_int_t  new_fid = o2n_hist->g_list[j] - 1;

      new_face_rank->g_list[shift] = new_fid;
      new_face_rank->g_elts[rank] += 1;

    } /* End of loop on new faces */

  } /* End of loop on initial faces */

  /* Free memory */

  BFT_FREE(reduce_ids);
  BFT_FREE(reduce_index);

  cs_join_gset_clean(new_face_rank);

  /* Define arrays to return */

  BFT_MALLOC(_send_rank_index, n_ranks + 1, cs_int_t);

  for (i = 0; i < n_ranks + 1; i++)
    _send_rank_index[i] = new_face_rank->index[i];

  BFT_MALLOC(_send_faces, _send_rank_index[n_ranks], cs_int_t);

  for (i = 0; i < _send_rank_index[n_ranks]; i++)
    _send_faces[i] = new_face_rank->g_list[i];

  cs_join_gset_destroy(&new_face_rank);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n Exchange to do after the splitting operation:\n");
  for (i = 0; i < n_ranks; i++) {
    start = _send_rank_index[i];
    end = _send_rank_index[i+1];
    bft_printf(" Send to rank %5d (n = %10d):", i, end - start);
    for (j = start; j < end; j++)
      bft_printf(" %d ", _send_faces[j]);
    bft_printf("\n");
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
_define_head_and_ext_edges(cs_int_t                face_id,
                           const cs_join_mesh_t   *mesh,
                           const cs_join_edges_t  *edges,
                           cs_join_rset_t         *head_edges,
                           cs_join_rset_t         *ext_edges,
                           cs_int_t                perm)
{
  cs_int_t  i, j, k, shift;
  cs_int_t  couple[2];

  cs_int_t  start_id = mesh->face_vtx_idx[face_id]-1;
  cs_int_t  end_id = mesh->face_vtx_idx[face_id+1]-1;
  cs_int_t  n_face_vertices = end_id - start_id;

  assert(perm < n_face_vertices);

  for (k = 0, i = start_id; i < end_id-1; i++) {

    couple[0] = mesh->face_vtx_lst[i];
    couple[1] = mesh->face_vtx_lst[i+1];

    for (j = edges->vtx_idx[couple[0]-1];
         j < edges->vtx_idx[couple[0]]; j++)
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

  for (j = edges->vtx_idx[couple[0]-1];
       j < edges->vtx_idx[couple[0]]; j++)
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
_create_face_builder(cs_int_t  n_faces)
{
  cs_int_t  i;

  face_builder_t  *builder = NULL;

  if (n_faces == 0)
    return NULL;

  BFT_MALLOC(builder, 1, face_builder_t);

  builder->n_faces = n_faces;

  BFT_MALLOC(builder->face_index, n_faces + 1, cs_int_t);

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
 *---------------------------------------------------------------------------*/

static void
_dump_face_builder(cs_int_t               face_id,
                   const face_builder_t  *builder)
{
  cs_int_t  i, j;
  cs_int_t  subface_id = 0;
  cs_int_t  face_s = builder->face_index[face_id];
  cs_int_t  face_e = builder->face_index[face_id+1];
  cs_int_t  n_subfaces = face_e - face_s;

  bft_printf(" Face %9d (n_subfaces: %d):\n", face_id+1, n_subfaces);

  for (i = face_s; i <face_e; i++, subface_id++) {

    cs_int_t  subface_s = builder->subface_index->array[i];
    cs_int_t  subface_e = builder->subface_index->array[i+1];

    if (builder->subface_gnum == NULL)
      bft_printf("   subface %4d: (%d, %d) -",
                 subface_id, subface_s, subface_e);
    else
      bft_printf("   subface %4d (%10u) - (%d, %d) -",
                 subface_id, builder->subface_gnum[face_s + subface_id],
                 subface_s, subface_e);

    if (builder->subface_gconnect == NULL) {
      for (j = subface_s; j < subface_e; j++)
        bft_printf(" %d ", builder->subface_connect->array[j]);
    }
    else {
      for (j = subface_s; j < subface_e; j++)
        bft_printf(" %u ", builder->subface_gconnect[j]);
    }
    bft_printf("\n");

  }
  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Split the current face into sub-faces under the "plane" tolerance
 * (check if two faces are coplanear).
 *
 * parameters:
 *   face_id       <-- face_id of the current in the cs_join_mesh_t struct.
 *   block_id      <-- current id in face_builder_t structure
 *   plane         <-- tolerance parameter to check coplanearity
 *   max_subfaces  <-- maximum number of sub faces
 *   verbosity     <-- level of accuracy in information display
 *   face_normal   <-- normal for each face of the mesh
 *   work          <-- cs_join_mesh_t structure
 *   edge_face_idx <-- "edge -> face" connect. index
 *   edge_face_lst <-- "edge -> face" connect. list
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
_split_face(cs_int_t                face_id,
            cs_int_t                block_id,
            double                  plane,
            int                     max_subfaces,
            int                     verbosity,
            const cs_real_t         face_normal[],
            const cs_join_mesh_t   *work,
            const cs_join_edges_t  *edges,
            const cs_int_t         *edge_face_idx,
            const cs_int_t         *edge_face_lst,
            face_builder_t         *builder,
            cs_join_rset_t        **head_edges,
            cs_join_rset_t        **subface_edges,
            cs_join_rset_t        **ext_edges,
            cs_join_rset_t        **int_edges)
{
  cs_int_t  i, j, k, i1, i2, i_int, i_ext;
  cs_int_t  first_vid, vid1, vid2;
  cs_int_t  subface_shift, connect_shift, connect_start;
  cs_int_t   next_vertex, next_edge;
  cs_real_t  face_norm[3], v1v2_vect[3], connect_vect[3];

  cs_real_t  min_max_d = 0.0;
  cs_int_t  n_subfaces = 0, head_edge_shift = 0;
  cs_int_t  start_id = work->face_vtx_idx[face_id]-1;
  cs_int_t  end_id = work->face_vtx_idx[face_id+1]-1;

  cs_real_t  max_coord[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
  cs_real_t  min_coord[3] = {DBL_MAX, DBL_MAX, DBL_MAX};

  cs_join_rset_t  *_head_edges = *head_edges;
  cs_join_rset_t  *_subface_edges = *subface_edges;
  cs_join_rset_t  *_ext_edges = *ext_edges;
  cs_join_rset_t  *_int_edges = *int_edges;

  const cs_join_vertex_t  *vertices = work->vertices;
  const cs_real_t  min_limit_cos = -1.1, max_limit_cos = 1.1;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  cs_bool_t  tst_dbg = (verbosity > 3 ? true : false);
#endif

  /* To be implemented ... */
  cs_int_t  *face_face_connect = NULL;

  /* Min./Max. coordinates of the current face */

  for (i = start_id; i < end_id; i++) {

    cs_int_t  vid = work->face_vtx_lst[i] - 1;

    for (j = 0; j < 3; j++) {
      min_coord[j] = CS_MIN(min_coord[j], vertices[vid].coord[j]);
      max_coord[j] = CS_MAX(max_coord[j], vertices[vid].coord[j]);
    }

  }

  for (j = 0; j < 3; j++)
    min_max_d = CS_MAX(min_max_d, max_coord[j] - min_coord[j]);

  for (j = 0; j < 3; j++) {
    min_coord[j] -= 0.5 * min_max_d;
    max_coord[j] += 0.5 * min_max_d;
  }

  /* Fill the current face normal */

  for (j = 0; j < 3; j++)
    face_norm[j] = face_normal[3*face_id+j];

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

      cs_int_t  head_edge_num = _head_edges->array[head_edge_shift];
      cs_int_t  edge_num = head_edge_num;
      cs_int_t  edge_id = FVM_ABS(edge_num) - 1;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
      if (tst_dbg)
        bft_printf(" face_gnum: %u, head_shift: %d, edge_num: %d\n",
                   work->face_gnum[face_id], head_edge_shift, edge_num);
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

#if 0 && defined(DEBUG) && !defined(NDEBUG)
      if (tst_dbg)
        bft_printf("  v1: %u, v2: %u, first: %u\n", vertices[vid1].gnum,
                   vertices[vid2].gnum, vertices[first_vid].gnum);
#endif

      while (vid2 != first_vid) {

        /* Look for the connected vertices and its associated edge */

        cs_int_t  v2v_start = edges->vtx_idx[vid2];
        cs_int_t  v2v_end = edges->vtx_idx[vid2+1];
        cs_int_t  n_connect_vertices = v2v_end - v2v_start;

        next_edge = 0;
        next_vertex = 0; /* To enable a check at the end */

        if (n_connect_vertices > 2) { /* Look for the edge which is
                                         the most on the left */

          cs_int_t  left_next_edge = -1, left_next_vertex = -1;
          cs_int_t  right_next_edge = -1, right_next_vertex = -1;
          cs_real_t  left_min_cos = max_limit_cos;
          cs_real_t  right_max_cos = min_limit_cos;

          for (j = 0; j < 3; j++)
            v1v2_vect[j] = vertices[vid2].coord[j]- vertices[vid1].coord[j];

          /* Loop on connected vertices */

          for (i = v2v_start; i < v2v_end; i++) {

            cs_int_t  connect_vid = edges->adj_vtx_lst[i]-1;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
            if (tst_dbg)
              bft_printf("  (%d) vertex connected to v2: %u\n",
                         i, vertices[connect_vid].gnum);
#endif

            if (connect_vid != vid1) {

              cs_bool_t  is_inside = true;
              cs_int_t  connect_edge_id = CS_ABS(edges->edge_lst[i]) - 1;

              /* Test if the connected vertex is inside the face */

              for (k = 0; k < 3; k++)
                if (  vertices[connect_vid].coord[k] < min_coord[k]
                   || vertices[connect_vid].coord[k] > max_coord[k])
                  is_inside = false;

              if (is_inside == true) {

                cs_bool_t  found = false;

                for (k = 0; k < 3; k++)
                  connect_vect[k] =  vertices[connect_vid].coord[k]
                                   - vertices[vid2].coord[k];

                /* Loop on faces sharing this edge */

                for (j = edge_face_idx[connect_edge_id];
                     (j < edge_face_idx[connect_edge_id+1] && found == false);
                     j++) {

                  cs_int_t  adj_face_id = edge_face_lst[j] - 1;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
                  if (tst_dbg)
                    bft_printf("\t Adj face (%u) through edge (%u)\n",
                               work->face_gnum[adj_face_id],
                               edges->gnum[connect_edge_id]);
#endif

                  if (adj_face_id != face_id) {

                    cs_real_t  dprod, dprod2;
                    cs_real_t  adj_face_norm[3];

                    for (k = 0; k < 3; k++)
                      adj_face_norm[k] = face_normal[3*adj_face_id+k];

                    dprod = _dot_product(adj_face_norm, face_norm);
                    dprod2 = dprod * dprod;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
                    if (tst_dbg)
                      bft_printf("\t dp: %g, dp2: %g Vs plane: %g -"
                                 " AdjNorm: [%g %g %g] - Norm: [%g %g %g]\n",
                                 dprod, dprod2, plane, adj_face_norm[0],
                                 adj_face_norm[1], adj_face_norm[2],
                                 face_norm[0], face_norm[1], face_norm[2]);
#endif

                    if (dprod2 > plane) {

                      if (face_face_connect == NULL) {
                        found = true;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
                        if (tst_dbg)
                          bft_printf("\t -> Adj face is OK\n");
#endif
                      }
                      else { /* To be implemented ... */

                        assert(face_face_connect != NULL);
                        bft_error(__FILE__, __LINE__, 0,
                                  _("  face splitting with face -> face"
                                    " connectivity is not yet implemented\n"));

                        /* Check if adj_face_id+1 is in face-face connect. */
                      }

                    }

                  } /* End if adj_face_id != face_id */
                  else {

                    assert(adj_face_id == face_id);
                    found = true;

                  }

                } /* End of loop on faces sharing this edge */

                if (found == true) {

                  /* Continue to build the new sub-face connectivity.
                     adj_edge_id is in the plane as edge_id

                     We look for the edge which is the most on the
                     left among all the adjacent edges.

                                                 |
                                             \   |   /
                    adj. edge : E2 ==>        \  |  /
                                               \ | /      Left part
                                                \|/
                    current edge : E1 ==> -->----o------  ............
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
                     in the right part with the biggest cos(A)

                  */

                  cs_real_t  cprod[3], dprod, cosine;

                  _cross_product(v1v2_vect, connect_vect, cprod);
                  dprod = _dot_product(cprod, face_norm);
                  cosine = _cosine(v1v2_vect, connect_vect);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
                  if (tst_dbg)
                    bft_printf("\t  dot_prod: %g, cosine: %g\n",
                               dprod, cosine);
#endif

                  if (dprod > 0) { /* Left part. We choose the edge
                                      with the smallest cosine */

                    if (cosine < left_min_cos) {
                      left_min_cos = cosine;
                      left_next_vertex = connect_vid + 1;
                      left_next_edge = edges->edge_lst[i];
                    }

                  } /* dprod > 0 */

                  else { /* In the right part. We choose the edge with
                            the biggest cosine. */

                    if (cosine > right_max_cos) {
                      right_max_cos = cosine;
                      right_next_vertex = connect_vid + 1;
                      right_next_edge = edges->edge_lst[i];
                    }

                  } /* End if dot_prod < 0 */

                } /* End if found = true */

              } /* The connected vertex is inside the face bounding box */

            } /* connect_vid != vid1 */

          } /* End of loop on connected vertices */

          if (left_min_cos < max_limit_cos) {
            next_edge = left_next_edge;
            next_vertex = left_next_vertex;
          }
          else if (right_max_cos > min_limit_cos) {
            next_edge = right_next_edge;
            next_vertex = right_next_vertex;
          }
          else if (   left_min_cos  >= max_limit_cos
                   && right_max_cos <= min_limit_cos) {

            /* Set return pointers */

            *head_edges = _head_edges;
            *ext_edges = _ext_edges;
            *int_edges = _int_edges;
            *subface_edges = _subface_edges;

            if (verbosity > 1)
              bft_printf("\nWarning: open cycle for global face %u\n",
                         work->face_gnum[face_id]);

            return OPEN_CYCLE_ERROR; /* open cycle */

          }

          assert(next_edge != 0);
          assert(next_vertex != 0);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
          if (tst_dbg)
            bft_printf(" Result: next_vertex: %u - next_edge: %u\n",
                       vertices[next_vertex-1].gnum,
                       edges->gnum[CS_ABS(next_edge)-1]);
#endif

        } /* End if n_connect_vertices > 2 */

        else if (n_connect_vertices == 2) {

          /* Loop on connected vertices */

          for (i = v2v_start; i < v2v_end; i++) {

            cs_int_t  connect_vid = edges->adj_vtx_lst[i]-1;

            if (connect_vid != vid1) {
              next_edge = edges->edge_lst[i];
              next_vertex = connect_vid + 1;
            }

          } /* End of loop on connected vertices */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
          if (tst_dbg)
            bft_printf(" Result (connect only by 2):"
                       " next_vertex: %u - next_edge: %u\n",
                       vertices[next_vertex-1].gnum,
                       edges->gnum[CS_ABS(next_edge)-1]);
#endif

        }
        else {

          assert(n_connect_vertices < 2);

          bft_error(__FILE__, __LINE__, 0,
                    _(" Joining operation : split face %d\n"
                      " Problem in the connectivity. Could not find a "
                      "connection with the vertex %d\n"),
                    face_id, vid1+1);


        } /* End of test on the number of vertices connected to vid2 */

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

          cs_int_t e1 = CS_ABS(_subface_edges->array[i1]);

          for (i2 = i1 + 1; i2 < _subface_edges->n_elts; i2++) {

            cs_int_t e2 = CS_ABS(_subface_edges->array[i2]);

            if (e1 == e2) {

              /* Returns pointers */

              *head_edges = _head_edges;
              *ext_edges = _ext_edges;
              *int_edges = _int_edges;
              *subface_edges = _subface_edges;

              if (verbosity > 1)
                bft_printf("\nWarning: global face %u traversed twice.\n",
                           work->face_gnum[face_id]);

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

        /* Test if next_edges is in the _int_edges list.
             If not : store next_edge in the _int_edges list.
             If yes : delete it.
        */

        if (i_ext != 0 && i_ext == _ext_edges->n_elts) {

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

        } /* Next_edge is an interior edge */

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

#if 0 && defined(DEBUG) && !defined(NDEBUG)
      if (tst_dbg)
        bft_printf(" END OF BUILDING subface %d\n\n", n_subfaces);
#endif

      if (n_subfaces > max_subfaces) { /* Too many sub-faces */

        /* Set return pointers */

        *head_edges = _head_edges;
        *ext_edges = _ext_edges;
        *int_edges = _int_edges;
        *subface_edges = _subface_edges;

        if (verbosity > 1)
          bft_printf("\nWarning: loop limit exceeded for global face %u\n.",
                     work->face_gnum[face_id]);

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

inline static cs_bool_t
_indexed_is_greater(size_t           i1,
                    size_t           i2,
                    const cs_int_t   index[],
                    const fvm_gnum_t number[])
{
  int  i;

  cs_int_t  i1_s = index[i1], i1_e = index[i1+1], s1 = i1_e - i1_s;
  cs_int_t  i2_s = index[i2], i2_e = index[i2+1], s2 = i2_e - i2_s;

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
 *
 * parameters:
 *   builder <-> pointer to a face builder structure
 *---------------------------------------------------------------------------*/

static void
_get_subface_gnum(face_builder_t  *builder)
{
  cs_int_t  i, j, k, shift;
  fvm_gnum_t  min_val;

  cs_int_t  max_size = 0;
  cs_int_t  n_subfaces = builder->face_index[builder->n_faces];
  cs_int_t  *index = builder->subface_index->array;
  fvm_gnum_t  *gconnect = builder->subface_gconnect;
  fvm_gnum_t  *glob_list = NULL, *tmp = NULL;
  fvm_io_num_t  *subface_io_num = NULL;

  const fvm_gnum_t  *global_num = NULL;

  assert(index != NULL);
  assert(gconnect != NULL);

  /* Allocate the buffer we want to define */

  BFT_MALLOC(builder->subface_gnum, n_subfaces, fvm_gnum_t);

  /* Re-arrange gconnect in order to have for each subface:
      - vertex with the minimal glob. num. in first place,
      - vertex with the minimal glob. num. between the two possible one
        in second place
  */

  for (i = 0; i < n_subfaces; i++)
    max_size = FVM_MAX(max_size, index[i+1] - index[i]);

  BFT_MALLOC(tmp, max_size, fvm_gnum_t);
  BFT_MALLOC(glob_list, index[n_subfaces], fvm_gnum_t);

  /* Build glob_list */

  for (i = 0; i < n_subfaces; i++) {

    cs_int_t  start = index[i], end = index[i+1];
    cs_int_t  n_elts = end - start;

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

  /* Copy glob_list as the new subface global connectivity */

  for (i = 0; i < index[n_subfaces]; i++)
    gconnect[i] = glob_list[i];

  if (cs_glob_n_ranks > 1) { /* Parallel treatment */

    /* Local ordering */

    cs_int_t  *order_index = NULL;
    fvm_gnum_t  *order_glob_list = NULL;
    fvm_lnum_t  *order = fvm_order_local_i(NULL,
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

    fvm_gnum_t  gnum = 1;
    fvm_lnum_t  *order = fvm_order_local_i(NULL,
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
    for (i = 0; i < n_subfaces; i++) {

      cs_int_t  start = index[i], end = index[i+1];
      cs_int_t  n_elts = end - start;

      bft_printf(" subface %5d - gnum: %u - connect_size: %d - ",
                 i+1, builder->subface_gnum[i], n_elts);
      for (j = start; j < end; j++)
        bft_printf(" %u ", glob_list[j]);
      bft_printf("\n");

    }
    bft_printf_flush();
#endif

    BFT_FREE(order);

  }

  bft_printf(_("\n  Global number of faces after splitting: %10u\n"),
             builder->n_g_subfaces);
  bft_printf_flush();

  /* Free memory */

  BFT_FREE(tmp);
  BFT_FREE(glob_list);

}

/*----------------------------------------------------------------------------
 * Update a cs_join_mesh_t structure thanks to a face_builder_t structure.
 *
 * parameters:
 *   block_info <-- set of paramaters defining a contiguous distribution
 *   builder    <-- pointer to the distributed face builder structure
 *   mesh       <-> pointer to the local cs_join_mesh_t structure
 *---------------------------------------------------------------------------*/

static void
_update_mesh_after_split(cs_join_block_info_t    block_info,
                         face_builder_t         *builder,
                         cs_join_mesh_t        **mesh)
{
  cs_int_t  i, j, k, id, shift, n_subfaces, o_id;
  fvm_gnum_t  prev, cur;

  cs_int_t  n_new_faces = 0;
  char  *new_mesh_name = NULL;
  cs_int_t  *subfaces = NULL;
  fvm_lnum_t  *order = NULL;
  cs_join_mesh_t  *init_mesh = *mesh;
  cs_join_mesh_t  *new_mesh = NULL;

  /* Sanity checks */

  assert(init_mesh != NULL);

  /* Create a new cs_join_mesh_t structure */

  BFT_MALLOC(new_mesh_name, strlen("AfterSplitting_n") + 5 + 1, char);
  sprintf(new_mesh_name,"%s%05d", "AfterSplitting_n",
          CS_MAX(cs_glob_rank_id, 0));

  new_mesh = cs_join_mesh_create(new_mesh_name);

  BFT_FREE(new_mesh_name);

  if (builder == NULL) { /* No face to treat for the current rank */

    assert(block_info.local_size == 0);

    cs_join_mesh_destroy(&init_mesh);

    /* Return pointers */

    *mesh = new_mesh;

    return;
  }

  assert((int)block_info.local_size == builder->n_faces);

  /* Compute the number of new faces */

  n_subfaces = builder->face_index[builder->n_faces];

  BFT_MALLOC(order, n_subfaces, fvm_lnum_t);

  fvm_order_local_allocated(NULL, builder->subface_gnum, order, n_subfaces);

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

  BFT_MALLOC(subfaces, n_new_faces, cs_int_t);

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

  BFT_MALLOC(new_mesh->face_gnum, n_new_faces, fvm_gnum_t);
  BFT_MALLOC(new_mesh->face_vtx_idx, n_new_faces + 1, cs_int_t);

  for (i = 0; i < n_new_faces; i++) {

    id = subfaces[i];
    new_mesh->face_gnum[i] = builder->subface_gnum[id];
    new_mesh->face_vtx_idx[i+1] =  builder->subface_index->array[id+1]
                                 - builder->subface_index->array[id];

  } /* End of loop on new faces */

  new_mesh->face_vtx_idx[0] = 1;
  for (i = 0; i < n_new_faces; i++)
    new_mesh->face_vtx_idx[i+1] += new_mesh->face_vtx_idx[i];

  BFT_MALLOC(new_mesh->face_vtx_lst,
             new_mesh->face_vtx_idx[n_new_faces]-1,
             cs_int_t);

  for (i = 0; i < n_new_faces; i++) {

    id = subfaces[i];
    shift = new_mesh->face_vtx_idx[i] - 1;

    for (j = builder->subface_index->array[id], k = 0;
         j < builder->subface_index->array[id+1]; j++, k++)
      new_mesh->face_vtx_lst[shift + k] = builder->subface_connect->array[j];

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

  cs_join_mesh_destroy(&init_mesh);

  /* Set return pointer */

  *mesh = new_mesh;
}

/*----------------------------------------------------------------------------
 * For each new sub-face we keep a relation between new and old
 * face global number.
 *
 * parameters:
 *   mesh       <-- mesh on sub-faces after splitting
 *   block_info <-- set of paramaters defining a contiguous distribution
 *   builder    <-- pointer to the distributed face builder structure
 *
 * returns:
 *   a pointer to a cs_join_gset_t structure saving relation between new
 *   and old faces
 *---------------------------------------------------------------------------*/

static cs_join_gset_t *
_keep_history(cs_join_mesh_t        *mesh,
              cs_join_block_info_t   block_info,
              face_builder_t        *builder)
{
  cs_int_t  i, j;

  cs_join_gset_t  *o2n_hist = NULL;

  assert(mesh != NULL);

  /* Create structure in which we keep the history of each new face */

  o2n_hist = cs_join_gset_create(block_info.local_size);

  if (block_info.local_size > 0) {

    assert(builder != NULL);
    assert(builder->n_faces == (cs_int_t)block_info.local_size);

    /* Historic is a part of the data held in builder structure */

    for (i = 0; i < (cs_int_t)block_info.local_size; i++)
      /* store old glob. face num. */
      o2n_hist->g_elts[i] = block_info.first_gnum + i;

    for (i = 0; i < builder->n_faces + 1; i++)
      o2n_hist->index[i] = builder->face_index[i];

    BFT_MALLOC(o2n_hist->g_list,
               o2n_hist->index[o2n_hist->n_elts],
               fvm_gnum_t);

    for (i = 0; i < builder->n_faces; i++) {
      for (j = builder->face_index[i]; j < builder->face_index[i+1]; j++) {

        cs_int_t  id = cs_search_g_binary(mesh->n_faces,
                                          builder->subface_gnum[j],
                                          mesh->face_gnum);

        assert(id != -1);
        o2n_hist->g_list[j] = id + 1;  /* store local face num. */

      }
    }

  } /* block.local_size > 0 */

  /* Free memory */

  return  o2n_hist;
}

/*============================================================================
 * Public function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Build new faces after the vertex fusion operation. Split initial faces into
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
  cs_int_t  i, j, face_s, subface_s, block_id;
  cs_join_split_error_t  code;
  cs_join_block_info_t  block_info;

  cs_int_t  _n_problems = 0, n_face_problems = 0, n_max_face_vertices = 0;
  cs_int_t  *edge_face_idx = NULL, *edge_face_lst = NULL;
  cs_join_gset_t  *_old2new_history = NULL;
  cs_join_rset_t  *open_cycle = NULL, *edge_traversed_twice = NULL;
  cs_join_rset_t  *loop_limit = NULL, *head_edges = NULL;
  cs_join_rset_t  *subface_edges = NULL, *ext_edges = NULL, *int_edges = NULL;
  face_builder_t  *loc_builder = NULL;
  cs_join_mesh_t  *_work = *work;

  const double plane = cos(param.plane *acos(-1.0)/180.);

  const cs_int_t  n_init_faces = _work->n_faces;
  const int  n_ranks = cs_glob_n_ranks;
  const int  local_rank = CS_MAX(cs_glob_rank_id, 0);

  assert(_work != NULL);
  assert(edges != NULL);

  /* Use the cs_join_edges_t structure to build
     the "edge -> face" connectivity */

  cs_join_mesh_get_edge_face_adj(_work, edges, &edge_face_idx, &edge_face_lst);

  /* Define buffers to manage errors */

  open_cycle = cs_join_rset_create(2);
  edge_traversed_twice = cs_join_rset_create(2);
  loop_limit = cs_join_rset_create(2);

  /* Define buffers and structures to build the new faces */

  head_edges = cs_join_rset_create(5);
  subface_edges = cs_join_rset_create(5);
  ext_edges = cs_join_rset_create(5);
  int_edges = cs_join_rset_create(5);

  /* Compute block_size */

  block_info = cs_join_get_block_info(_work->n_g_faces, n_ranks, local_rank);

  loc_builder = _create_face_builder(block_info.local_size);

  /*
     We only have to treat faces for the current rank's block because the
     initial face distribution assumes that the current rank can correctly
     split only faces in its block.
     Main loop on faces.
  */

  for (i = 0, block_id = 0; i < n_init_faces; i++) {

    int  block_rank = (_work->face_gnum[i] - 1)/block_info.size;

    if (block_rank == local_rank) { /* This face is a "main" face for the
                                       local rank */

      cs_int_t  n_face_vertices =
        _work->face_vtx_idx[i+1] - _work->face_vtx_idx[i];

      /* Manage list size */

      if (n_face_vertices > n_max_face_vertices) {

        n_max_face_vertices = n_face_vertices;
        cs_join_rset_resize(&head_edges, n_face_vertices);
        cs_join_rset_resize(&subface_edges, n_face_vertices);
        cs_join_rset_resize(&ext_edges, n_face_vertices);
        cs_join_rset_resize(&int_edges, n_face_vertices);

      }

      /* Fill head_edges and ext_edges */

      _define_head_and_ext_edges(i,  /* face_id */
                                 _work,
                                 edges,
                                 head_edges,
                                 ext_edges,
                                 0); /* No permutation */

      /* Store initial loc_builder state in case of code > 0 */

      face_s = loc_builder->face_index[block_id];
      subface_s = loc_builder->subface_index->array[face_s];

      /* Split the current face into subfaces */

      code = _split_face(i,       /* face_id */
                         block_id,
                         plane,
                         param.max_sub_faces,
                         param.verbosity,
                         face_normal,
                         _work,
                         edges,
                         edge_face_idx,
                         edge_face_lst,
                         loc_builder,
                         &head_edges,
                         &subface_edges,
                         &ext_edges,
                         &int_edges);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
      if (param.verbosity > 1 && code != NO_SPLIT_ERROR) {
        bft_printf(_("  Split face: %d with returned code: %d\n"), i+1, code);
        _dump_face_builder(block_id, loc_builder);
      }
#endif

      if (code > NO_SPLIT_ERROR) { /* Manage error */

        _n_problems++;

        /* We change the starting edge for traversing edges to build
           the new face. This may be enough to solve the problem when
           the face is warped. */

        while (_n_problems < n_face_vertices && code > NO_SPLIT_ERROR) {

          /* Fill head_edges and ext_edges */

          _define_head_and_ext_edges(i, /* face_id */
                                     _work,
                                     edges,
                                     head_edges,
                                     ext_edges,
                                     _n_problems); /* permutation */

          /* Retrieve initial loc_builder state */

          loc_builder->face_index[block_id] = face_s;
          loc_builder->subface_index->array[face_s] = subface_s;
          loc_builder->subface_connect->n_elts = subface_s;

          /* Split the current face into subfaces */

          code = _split_face(i,   /* face_id */
                             block_id,
                             plane,
                             param.max_sub_faces,
                             param.verbosity,
                             face_normal,
                             _work,
                             edges,
                             edge_face_idx,
                             edge_face_lst,
                             loc_builder,
                             &head_edges,
                             &subface_edges,
                             &ext_edges,
                             &int_edges);

          _n_problems++;

        } /* End of while */

        if (_n_problems >= n_face_vertices && code != NO_SPLIT_ERROR) {

          n_face_problems++;

          switch (code) {

          case OPEN_CYCLE_ERROR:
            cs_join_rset_resize(&open_cycle, open_cycle->n_elts);
            open_cycle->array[open_cycle->n_elts] = i + 1;
            open_cycle->n_elts += 1;
            break;

          case EDGE_TRAVERSED_TWICE_ERROR:
            cs_join_rset_resize(&edge_traversed_twice,
                                edge_traversed_twice->n_elts);
            edge_traversed_twice->array[edge_traversed_twice->n_elts] = i + 1;
            edge_traversed_twice->n_elts += 1;
            break;

          case LOOP_LIMIT_ERROR:
            cs_join_rset_resize(&loop_limit, loop_limit->n_elts);
            loop_limit->array[loop_limit->n_elts] = i + 1;
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

          loc_builder->face_index[block_id] = face_s;
          loc_builder->face_index[block_id+1] = face_s + 1;

          /* face -> subface connectivity index update */

          cs_join_rset_resize(&(loc_builder->subface_index), face_s+1);

          loc_builder->subface_index->n_elts = face_s + 1;
          loc_builder->subface_index->array[face_s] = subface_s;
          loc_builder->subface_index->array[face_s+1]
            = subface_s + n_face_vertices;

          /* face -> subface connectivity list update */

          cs_join_rset_resize(&(loc_builder->subface_connect),
                              subface_s + n_face_vertices);

          loc_builder->subface_connect->n_elts = subface_s + n_face_vertices;

          for (j = 0; j < n_face_vertices; j++)
            loc_builder->subface_connect->array[subface_s + j]
              = _work->face_vtx_lst[j + _work->face_vtx_idx[i] - 1];

          if (param.verbosity > 1) {
            bft_printf("\n Keep initial connectivity for face %d (%u):\n",
                       i+1, _work->face_gnum[i]);
            _dump_face_builder(block_id, loc_builder);
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
    fvm_gnum_t  g_in_buf[4], g_out_buf[4];

    fvm_gnum_t  n_g_face_problems = n_face_problems;
    fvm_gnum_t  n_g_open_cycles = open_cycle->n_elts;
    fvm_gnum_t  n_g_edges_twice = edge_traversed_twice->n_elts;
    fvm_gnum_t  n_g_loop_limit = loop_limit->n_elts;

#if defined(HAVE_MPI)
    if (n_ranks > 1) {

      MPI_Comm  mpi_comm = cs_glob_mpi_comm;

      g_in_buf[0] = n_g_face_problems;
      g_in_buf[1] = n_g_open_cycles;
      g_in_buf[2] = n_g_edges_twice;
      g_in_buf[3] = n_g_loop_limit;

      MPI_Allreduce(g_in_buf, g_out_buf, 4, FVM_MPI_GNUM, MPI_SUM, mpi_comm);

      n_g_face_problems = g_out_buf[0];
      n_g_open_cycles = g_out_buf[1];
      n_g_edges_twice = g_out_buf[2];
      n_g_loop_limit = g_out_buf[3];

    }
#endif

    if (n_g_face_problems > 0) {

      bft_printf
        (_("  Warning: (%lu) problem(s) found during the face splitting\n"
           "     %12lu  open cycles,\n"
           "     %12lu  edges traversed twice,\n"
           "     %12lu  faces split into more than "
           "max_subfaces (= %d)\n\n"
           "  => Eventually modify joining parameters\n\n"),
         (unsigned long)n_g_face_problems,
         (unsigned long)n_g_open_cycles,
         (unsigned long)n_g_edges_twice,
         (unsigned long)n_g_loop_limit, param.max_sub_faces);
      bft_printf_flush();

      assert(   n_g_face_problems
             == n_g_open_cycles + n_g_edges_twice + n_g_loop_limit);

      /* post-processing of encountered problems */

      cs_join_post_init(); /* Init. post-processing anyway */

      if (n_g_open_cycles > 0)
        cs_join_post_faces_subset("OpenCycleErr",
                                  _work,
                                  open_cycle->n_elts,
                                  open_cycle->array);

      if (n_g_edges_twice > 0)
        cs_join_post_faces_subset("EdgeScannedTwiceErr",
                                  _work,
                                  edge_traversed_twice->n_elts,
                                  edge_traversed_twice->array);

      if (n_g_loop_limit > 0)
        cs_join_post_faces_subset("LoopLimitErr",
                                  _work,
                                  loop_limit->n_elts,
                                  loop_limit->array);

    } /* End of information display */

  } /* Manage potential errors during face splitting */

  /* Free vtx_struct structure */

  BFT_FREE(edge_face_idx);
  BFT_FREE(edge_face_lst);

  /* Define a global number for each new sub-faces */

  if (loc_builder != NULL) {

    cs_int_t  n_subfaces = loc_builder->face_index[loc_builder->n_faces];
    cs_int_t  sub_connect_size = loc_builder->subface_index->array[n_subfaces];

    /* Define subface_gconnect */

    BFT_MALLOC(loc_builder->subface_gconnect, sub_connect_size, fvm_gnum_t);

    for (i = 0; i < sub_connect_size; i++) {

      cs_int_t  vid = loc_builder->subface_connect->array[i] - 1;
      fvm_gnum_t  vgnum = _work->vertices[vid].gnum;

      loc_builder->subface_gconnect[i] = vgnum;

    }

    if (sub_connect_size > 0)
      _get_subface_gnum(loc_builder);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    bft_printf("\nFINAL BUILDER STATE\n");
    for (i = 0; i < loc_builder->n_faces; i++)
      _dump_face_builder(i, loc_builder);
#endif

    BFT_FREE(loc_builder->subface_gconnect);

  }

  /* Delete error management lists */

  cs_join_rset_destroy(&open_cycle);
  cs_join_rset_destroy(&edge_traversed_twice);
  cs_join_rset_destroy(&loop_limit);

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
  */

  /* Reduce the definition of the working mesh to the set of new
     global face numbers built from the local block. */

  _update_mesh_after_split(block_info, loc_builder, &_work);

  /* For each new sub-face we maintain a relation between new and old
     face global number */

  _old2new_history = _keep_history(_work,
                                   block_info,
                                   loc_builder);

  /* Free face_builder_t structure */

  loc_builder = _destroy_face_builder(loc_builder);

  /* Set return pointers */

  *work = _work;
  *old2new_history = _old2new_history;
}

/*----------------------------------------------------------------------------
 * Update after face splitting of the local join mesh structure.
 * Send back to the original rank the new face description.
 *
 * parameters:
 *   work_mesh       <-- distributed mesh on faces to join
 *   gnum_rank_index <-- index on ranks for the old global face numbering
 *   o2n_hist        <-> old global face -> new local face numbering
 *   local_mesh      <-> mesh on local selected faces to be joined
 *---------------------------------------------------------------------------*/

void
cs_join_split_update_struct(const cs_join_mesh_t   *work_mesh,
                            const fvm_gnum_t        gnum_rank_index[],
                            cs_join_gset_t        **o2n_hist,
                            cs_join_mesh_t        **local_mesh)
{
  cs_int_t  i;

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

    cs_int_t  j, subface_id;
    cs_int_t  *send_rank_index = NULL, *send_faces = NULL;
    fvm_gnum_t  *init_face_gnum = NULL;
    cs_join_gset_t  *distrib_sync_hist = NULL;
    cs_int_t  n_init_faces = _local_mesh->n_faces;
    fvm_gnum_t  n_g_init_faces = _local_mesh->n_g_faces;

    MPI_Comm  mpi_comm = fvm_parall_get_mpi_comm();

    /* Save the initial global face numbering */

    BFT_MALLOC(init_face_gnum, n_init_faces, fvm_gnum_t);

    for (i = 0; i < n_init_faces; i++)
      init_face_gnum[i] = _local_mesh->face_gnum[i];

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

    distrib_sync_hist = cs_join_gset_sync_by_block(n_g_init_faces,
                                                   _o2n_hist,
                                                   mpi_comm);

    cs_join_gset_destroy(&_o2n_hist);
    _o2n_hist = cs_join_gset_create(n_init_faces);

    for (i = 0; i < n_init_faces; i++)
      _o2n_hist->g_elts[i] = init_face_gnum[i];

    BFT_FREE(init_face_gnum);

    cs_join_gset_update_from_block(n_g_init_faces,
                                   distrib_sync_hist,
                                   _o2n_hist,
                                   mpi_comm);

    cs_join_gset_destroy(&distrib_sync_hist);


  }
#endif /* HAVE_MPI */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  { /* Full dump of structures */

    int  len;
    FILE  *dbg_file = NULL;
    char  *filename = NULL;

    len = strlen("JoinDBG_o2nFaceHist.dat")+1+4;
    BFT_MALLOC(filename, len, char);
    sprintf(filename, "JoinDBG_o2nFaceHist%04d.dat",
            CS_MAX(cs_glob_rank_id, 0));
    dbg_file = fopen(filename, "w");

    cs_join_gset_dump(dbg_file, _o2n_hist);

    fflush(dbg_file);
    BFT_FREE(filename);
    fclose(dbg_file);

  }
#endif

  /* Set return pointers */

  *o2n_hist = _o2n_hist;
  *local_mesh = _local_mesh;
}

/*---------------------------------------------------------------------------*/

END_C_DECLS

