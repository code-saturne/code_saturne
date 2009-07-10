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
 *============================================================================*/

/*============================================================================
 * Manipulation of a cs_join_mesh_t structure
 *===========================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *---------------------------------------------------------------------------*/

#include <assert.h>
#include <string.h>
#include <math.h>

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
#include <fvm_nodal.h>
#include <fvm_nodal_from_desc.h>
#include <fvm_nodal_order.h>
#include <fvm_order.h>
#include <fvm_parall.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *---------------------------------------------------------------------------*/

#include "cs_search.h"
#include "cs_join_post.h"
#include "cs_join_set.h"
#include "cs_join_util.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *---------------------------------------------------------------------------*/

#include "cs_join_mesh.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro and type definitions
 *===========================================================================*/

/*============================================================================
 * Private function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Compute the cross product of two vectors.
 *
 * parameters:
 *   v1 <-- first vector
 *   v2 <-- second vector
 *
 * returns:
 *   the resulting cross product (v1 x v2)
 *----------------------------------------------------------------------------*/

inline static void
_cross_product(const double  v1[],
               const double  v2[],
               double        result[])
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
_dot_product(const double  v1[],
             const double  v2[])
{
  int  i;
  double  result = 0.0;

  for (i = 0; i < 3; i++)
    result += v1[i] * v2[i];

  return result;
}

/*----------------------------------------------------------------------------
 * Find for each face of the list its related rank
 *
 * parameters:
 *   n_elts     <-- number of elements in the glob. list
 *   glob_list  <-- global numbering list (must be ordered)
 *   rank_index <-- index defining a range of global numbers related
 *                  to a rank
 *
 * returns:
 *   an array of size n_elts
 *---------------------------------------------------------------------------*/

static cs_int_t *
_get_rank_from_index(cs_int_t          n_elts,
                     const fvm_gnum_t  glob_list[],
                     const fvm_gnum_t  rank_index[])
{
  cs_int_t  i, rank;

  cs_int_t  *rank_list = NULL;

  if (n_elts == 0)
    return NULL;

  BFT_MALLOC(rank_list, n_elts, cs_int_t);

  for (i = 0, rank = 0; i < n_elts; i++) {

    for (;rank_index[rank+1] < glob_list[i]; rank++);

    assert(rank < fvm_parall_get_size());
    rank_list[i] = rank;

  } /* End of loop on elements */

  return rank_list;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Get the index on ranks and th list of faces to send from a list of global
 * faces to receive.
 *
 * parameters:
 *   n_ranks           <-- number of ranks
 *   gnum_rank_index   <-- index on ranks for the global elements
 *   n_elts            <-- number of elements to get
 *   glob_list         <-- global number of faces to get (ordered)
 *   p_send_rank_index --> index on ranks for the faces to send
 *   p_send_faces      --> list of face ids to send
 *---------------------------------------------------------------------------*/

static void
_get_send_faces(int                n_ranks,
                const fvm_gnum_t   gnum_rank_index[],
                cs_int_t           n_elts,
                const fvm_gnum_t   glob_list[],
                cs_int_t          *send_rank_index[],
                cs_int_t          *send_faces[])
{
  int  i, rank, shift;
  fvm_gnum_t  first_gface_id;

  cs_int_t  *gface_ranks = NULL, *_send_faces = NULL, *_send_rank_index = NULL;
  cs_int_t  *send_count = NULL, *recv_count = NULL, *send_shift = NULL;
  fvm_gnum_t  *gfaces_to_send = NULL, *gfaces_to_recv = NULL;

  MPI_Comm  comm = cs_glob_mpi_comm;

  const int  local_rank = cs_glob_rank_id;

  /* Sanity checks */

  assert(n_ranks > 1);
  assert(gnum_rank_index != NULL);

  /* Find for each element of the list, the rank which owns the element */

  gface_ranks = _get_rank_from_index(n_elts, glob_list, gnum_rank_index);

  first_gface_id = gnum_rank_index[local_rank];

  /* Count the number of faces for which we want the connectivity */

  BFT_MALLOC(send_count, n_ranks, cs_int_t);
  BFT_MALLOC(recv_count, n_ranks, cs_int_t);

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (i = 0; i < n_elts; i++)
    send_count[gface_ranks[i]] += 1;

  /* Exchange number of elements for which we want a connectivity */

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

  /* Build index arrays */

  BFT_MALLOC(send_shift, n_ranks + 1, cs_int_t);
  BFT_MALLOC(_send_rank_index, n_ranks + 1, cs_int_t);

  send_shift[0] = 0;
  _send_rank_index[0] = 0;

  for (rank = 0; rank < n_ranks; rank++) {

    _send_rank_index[rank+1] = _send_rank_index[rank] + recv_count[rank];
    send_shift[rank+1] = send_shift[rank] + send_count[rank];

  }

  /* Build gfaces_to_recv = glob_list but potentially in a different order.
     List of face (global numbering) for which we want the connectivity */

  BFT_MALLOC(gfaces_to_recv, send_shift[n_ranks], fvm_gnum_t);
  BFT_MALLOC(gfaces_to_send, _send_rank_index[n_ranks], fvm_gnum_t);

  assert(send_shift[n_ranks] == n_elts);

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (i = 0; i < n_elts; i++) {

    rank = gface_ranks[i];
    shift = send_count[rank] + send_shift[rank];
    gfaces_to_recv[shift] = glob_list[i];
    send_count[rank] += 1;

  }

  /* Exchange list of global num. to exchange */

  MPI_Alltoallv(gfaces_to_recv, send_count, send_shift, FVM_MPI_GNUM,
                gfaces_to_send, recv_count, _send_rank_index, FVM_MPI_GNUM,
                comm);

  BFT_MALLOC(_send_faces, _send_rank_index[n_ranks], cs_int_t);

  /* Define face ids to send */

  for (rank = 0; rank < n_ranks; rank++) {

    for (i = _send_rank_index[rank]; i < _send_rank_index[rank + 1]; i++)
      _send_faces[i] = gfaces_to_send[i] - 1 - first_gface_id;

  } /* End of loop on ranks */

  /* Free memory */

  BFT_FREE(gface_ranks);
  BFT_FREE(gfaces_to_recv);
  BFT_FREE(send_shift);
  BFT_FREE(send_count);
  BFT_FREE(recv_count);
  BFT_FREE(gfaces_to_send);

  /* Set return pointers */

  *send_rank_index = _send_rank_index;
  *send_faces = _send_faces;
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Clean the given cs_join_mesh_t structure: remove empty edges.
 *
 * parameters:
 *   mesh      <-> pointer to the cs_join_mesh_t structure to clean
 *   verbosity <-- level of display
 *---------------------------------------------------------------------------*/

static void
_remove_empty_edges(cs_join_mesh_t  *mesh,
                    int              verbosity)
{
  cs_int_t  i, j, n_face_vertices;

  cs_int_t  shift = 0, n_simplified_faces = 0;
  cs_int_t  *new_face_vtx_idx = NULL;

  BFT_MALLOC(new_face_vtx_idx, mesh->n_faces + 1, cs_int_t);

  new_face_vtx_idx[0] = 1;

  for (i = 0; i < mesh->n_faces; i++) {

    cs_int_t  s = mesh->face_vtx_idx[i] - 1;
    cs_int_t  e = mesh->face_vtx_idx[i+1] - 1;

    mesh->face_vtx_lst[shift++] = mesh->face_vtx_lst[s];

    /* Loop on face vertices */

    for (j = s + 1; j < e - 1; j++)
      if (mesh->face_vtx_lst[j] != mesh->face_vtx_lst[j+1])
        mesh->face_vtx_lst[shift++] = mesh->face_vtx_lst[j];

    if (mesh->face_vtx_lst[e-1] != mesh->face_vtx_lst[s])
      mesh->face_vtx_lst[shift++] = mesh->face_vtx_lst[e-1];

    new_face_vtx_idx[i+1] = shift + 1;

    n_face_vertices = new_face_vtx_idx[i+1] - new_face_vtx_idx[i];

    if (n_face_vertices < e - s) {

      n_simplified_faces++;
      if (verbosity > 2)
        bft_printf("  Simplified face %d (%u)\n", i+1, mesh->face_gnum[i]);

      if (n_face_vertices < 3)
        bft_error(__FILE__, __LINE__, 0,
                  _("  The simplified face has less than 3 vertices.\n"
                    "  Check your joining parameters.\n"
                    "  Face %d (%u)\n"), i+1, mesh->face_gnum[i]);
    }

  } /* End of loop on faces */

  BFT_FREE(mesh->face_vtx_idx);
  mesh->face_vtx_idx = new_face_vtx_idx;

  BFT_REALLOC(mesh->face_vtx_lst,
              new_face_vtx_idx[mesh->n_faces]-1, cs_int_t);

  if (verbosity > 1) {
    fvm_gnum_t n_g_simplified_faces = n_simplified_faces;
    fvm_parall_counter(&n_g_simplified_faces, 1);
    bft_printf(_("\n  Number of simplified faces: %lu\n"),
               (unsigned long)n_simplified_faces);
  }
}

/*----------------------------------------------------------------------------
 * Clean the given cs_join_mesh_t structure: remove degenerate edges.
 *
 * parameters:
 *   mesh      <-> pointer to the cs_join_mesh_t structure to clean
 *   verbosity <-- level of display
 *---------------------------------------------------------------------------*/

static void
_remove_degenerate_edges(cs_join_mesh_t  *mesh,
                         int              verbosity)
{
  /*
    - In the definition of faces based on new edges, a same edge may be
      traversed twice, in the opposite direction; this is due to merging
      of edges, as shown below.

       x                                      x
       |\                                     |
       | \                                    |
     a2|  \a3                               A2|
       |   \                                  |
       |    \      a4       Merge of          |
    ---s1----s2------       vertices          x
       |      \             s1 and s2        / \
       |       \                            /   \
     a1|        \a4                      A1/     \A3
       |         \                        /       \
       |          \                      /         \
       x-----------x                    x-----------x
            a5                                A4


    Face: a1 a2 a3 a4 a5            Face: A1 A2 -A2 A3 A4


   Caution:    the final configuration may be
               A2 A1 A3 A4 -A2
               where the references of edges to delete may be at the
               beginning or end of the face definition

   Remark:     several edge pairs may possibly be referenced twice,
               in the form
               ... A1 A2 -A2 -A1 ... (where the removal of A2 will
               make ... A1 -A1 ... appear); We thus run as many passes
               as necessary on a given face.
  */

  cs_int_t  i, j, k, count, n_face_vertices;

  cs_int_t  shift = 0;
  cs_int_t  n_faces = mesh->n_faces;
  cs_int_t  n_modified_faces = 0;
  fvm_gnum_t  n_g_modified_faces = 0;
  cs_join_rset_t  *tmp = NULL;
  cs_join_rset_t  *kill = NULL;

  tmp = cs_join_rset_create(8);
  kill = cs_join_rset_create(8);

  for (i = 0; i < n_faces; i++) {

    cs_int_t  start_id = mesh->face_vtx_idx[i] - 1;
    cs_int_t  end_id = mesh->face_vtx_idx[i+1] - 1;
    cs_int_t  n_init_vertices = end_id - start_id;
    cs_int_t  n_elts = n_init_vertices + 2;

    assert(n_init_vertices > 2);

    /* Build a temporary list based on the face connectivity */

    cs_join_rset_resize(&tmp, n_elts);
    cs_join_rset_resize(&kill, n_elts);

    for (j = start_id, k = 0; j < end_id; j++, k++) {
      tmp->array[k] = mesh->face_vtx_lst[j];
      kill->array[k] = 0;
    }

    tmp->array[k] = mesh->face_vtx_lst[start_id];
    kill->array[k++] = 0;
    tmp->array[k] = mesh->face_vtx_lst[start_id+1];
    kill->array[k++] = 0;

    assert(n_elts == k);
    tmp->n_elts = n_elts;
    kill->n_elts = n_elts;

    /* Find degenerate edges */

    count = 1;
    n_face_vertices = n_init_vertices;

    while (count > 0) {

      count = 0;
      for (j = 0; j < n_face_vertices; j++) {
        if (tmp->array[j] == tmp->array[j+2]) {
          count++;
          kill->array[j] = 1;
          kill->array[(j+1)%n_face_vertices] = 1;
        }
      }

      tmp->n_elts = 0;
      for (j = 0; j < n_face_vertices; j++) {
        if (kill->array[j] == 0)
          tmp->array[tmp->n_elts++] = tmp->array[j];
      }

      n_face_vertices = tmp->n_elts;
      tmp->array[tmp->n_elts++] = tmp->array[0];
      tmp->array[tmp->n_elts++] = tmp->array[1];

      kill->n_elts = tmp->n_elts;
      for (j = 0; j < kill->n_elts; j++)
        kill->array[j] = 0;

    } /* End of while */

    if (n_face_vertices != n_init_vertices) {

      n_modified_faces += 1;

      if (verbosity > 2) { /* Display the degenerate face */

        bft_printf("\n  Degenerate connectivity for face: %d [%u]:",
                   i+1, mesh->face_gnum[i]);

        bft_printf("\n    Initial def: ");
        for (j = start_id; j < end_id; j++) {
          cs_int_t  v_id = mesh->face_vtx_lst[j] - 1;
          bft_printf(" %d (%u) ", v_id+1, mesh->vertices[v_id].gnum);
        }

        bft_printf("\n    Final def:   ");
        for (j = 0; j < n_face_vertices; j++) {
          cs_int_t  v_id = tmp->array[j] - 1;
          bft_printf(" %d (%u) ", v_id+1, mesh->vertices[v_id].gnum);
        }

        bft_printf("\n");
        bft_printf_flush();
      }

    } /* End if n_face_vertices != n_init_vertices */

    for (j = 0; j < n_face_vertices; j++)
      mesh->face_vtx_lst[shift++] = tmp->array[j];
    mesh->face_vtx_idx[i] = shift;

  } /* End of loop on faces */

  fvm_parall_counter(&n_g_modified_faces, 1);

  bft_printf("  Degenerate connectivity for %lu faces.\n"
             "  Mesh cleaning done.\n",
             (unsigned long)n_g_modified_faces);

  for (i = n_faces; i > 0; i--)
    mesh->face_vtx_idx[i] = mesh->face_vtx_idx[i-1] + 1;
  mesh->face_vtx_idx[0] = 1;

  BFT_REALLOC(mesh->face_vtx_lst, mesh->face_vtx_idx[n_faces], cs_int_t);

  /* Free memory */

  cs_join_rset_destroy(&tmp);
  cs_join_rset_destroy(&kill);
}

/*----------------------------------------------------------------------------
 * Count the number of new vertices to add in the new face definition
 *
 * parameters:
 *   v1_num           <-- first vertex number
 *   v2_num           <-- second vertex number
 *   old2new          <-- indirection array between old and new numbering
 *   edges            <-- cs_join_edges_t structure
 *   edge_index       <-- edge -> new added vertex connectivity index
 *   edge_new_vtx_lst <-- edge -> new added vertex connectivity list
 *
 * returns:
 *   a number of vertices to add
 *---------------------------------------------------------------------------*/

static int
_count_new_added_vtx_to_edge(cs_int_t                v1_num,
                             cs_int_t                v2_num,
                             const cs_int_t          old2new[],
                             const cs_join_edges_t  *edges,
                             const cs_int_t          edge_index[],
                             const cs_int_t          edge_new_vtx_lst[])
{
  cs_int_t  i, edge_id, edge_num;

  cs_int_t  new_v1_num = old2new[v1_num-1] + 1;
  cs_int_t  new_v2_num = old2new[v2_num-1] + 1;
  cs_int_t  n_adds = 0;

  assert(v1_num > 0);
  assert(v2_num > 0);
  assert(new_v1_num > 0);
  assert(new_v2_num > 0);
  assert(edge_index != NULL);

  /* Find the related edge */

  edge_num = cs_join_mesh_get_edge(v1_num, v2_num, edges);
  edge_id = CS_ABS(edge_num) - 1;

  if (v1_num == v2_num)
    bft_error(__FILE__, __LINE__, 0,
              _("\n Problem in mesh connectivity.\n"
                " Detected when updating connectivity.\n"
                " Edge number: %d (%u) - (%d, %d) in old numbering.\n"),
              edge_num, edges->gnum[edge_id], v1_num, v2_num);

  /* Add the first vertex (new_v1_num) */

  n_adds = 1;

  /* Add another vertices if needed */

  for (i = edge_index[edge_id]; i < edge_index[edge_id+1]; i++) {

    cs_int_t  new_vtx_num = edge_new_vtx_lst[i];

    if (new_vtx_num != new_v1_num && new_vtx_num != new_v2_num)
      n_adds++;

  }

  return n_adds;
}

/*----------------------------------------------------------------------------
 * Add new vertex to the face -> vertex connectivity
 *
 * parameters:
 *   v1_num           <-- first vertex number
 *   v2_num           <-- second vertex number
 *   old2new          <-- indirection array between old and new numbering
 *   edges            <-- cs_join_edges_t structure
 *   edge_index       <-- edge -> new added vertex connectivity index
 *   edge_new_vtx_lst <-- edge -> new added vertex connectivity list
 *   new_face_vtx_lst <-> new face -> vertex connectivity list
 *   p_shift          <-> pointer to the shift in the connectivity list
 *---------------------------------------------------------------------------*/

static void
_add_new_vtx_to_edge(cs_int_t                v1_num,
                     cs_int_t                v2_num,
                     const cs_int_t          old2new[],
                     const cs_join_edges_t  *edges,
                     const cs_int_t          edge_index[],
                     const cs_int_t          edge_new_vtx_lst[],
                     cs_int_t                new_face_vtx_lst[],
                     cs_int_t               *p_shift)
{
  cs_int_t  new_v1_num = old2new[v1_num-1] + 1;
  cs_int_t  shift = *p_shift;

  /* Add first vertex num to the connectivity list */

  new_face_vtx_lst[shift++] = new_v1_num;

  if (edge_new_vtx_lst != NULL) {

    cs_int_t  i, edge_id, edge_num, e_start, e_end;

    cs_int_t  new_v2_num = old2new[v2_num-1] + 1;

    /* Find the related edge */

    edge_num = cs_join_mesh_get_edge(v1_num, v2_num, edges);
    edge_id = CS_ABS(edge_num) - 1;
    e_start = edge_index[edge_id];
    e_end = edge_index[edge_id+1];

    /* Add a vertex if needed */

    if (edge_num > 0) {

      for (i = e_start; i < e_end; i++) {

        cs_int_t  new_vtx_num = edge_new_vtx_lst[i];

        if (new_vtx_num != new_v1_num && new_vtx_num != new_v2_num)
          new_face_vtx_lst[shift++] = new_vtx_num;

      }
    }
    else { /* edge_num < 0 */

      for (i = e_end - 1; i > e_start - 1; i--) {

        cs_int_t  new_vtx_num = edge_new_vtx_lst[i];

        if (new_vtx_num != new_v1_num && new_vtx_num != new_v2_num)
          new_face_vtx_lst[shift++] = new_vtx_num;

      }

    } /* End if edge_num < 0 */

  } /* End if edge_new_vtx_lst != NULL */

  /* Return pointer */

  *p_shift = shift;

}

/*============================================================================
 * Public function definitions
 *===========================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Create a MPI_Datatype for the cs_join_vertex_t structure.
 *
 * returns:
 *   a MPI_Datatype associated to the cs_join_vertex_t structure
 *---------------------------------------------------------------------------*/

MPI_Datatype
cs_join_mesh_create_vtx_datatype(void)
{
  int  j;
  cs_join_vertex_t  v_data;
  MPI_Datatype  new_type;

  int  blocklengths[3] = {1, 1, 3};
  MPI_Aint  displacements[3] = {0 , 0, 0};
  MPI_Datatype  types[3] = {MPI_DOUBLE, FVM_MPI_GNUM, FVM_MPI_COORD};

  /* Initialize bbox */

  v_data.tolerance = 0.0;
  v_data.gnum = 1;
  for (j = 0; j < 3; j++)
    v_data.coord[j] = 0.;

  /* Define array of displacements */

  MPI_Get_address(&v_data, displacements);
  MPI_Get_address(&v_data.gnum, displacements + 1);
  MPI_Get_address(&v_data.coord, displacements + 2);

  displacements[2] -= displacements[0];
  displacements[1] -= displacements[0];
  displacements[0] = 0;

  /* Create new datatype */

  MPI_Type_create_struct(3, blocklengths, displacements, types, &new_type);

  MPI_Type_commit(&new_type);

  return new_type;
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Allocate and initialize a new cs_join_mesh_t structure.
 *
 * parameters:
 *   name <-- name of the mesh
 *
 * returns:
 *   a pointer to a cs_join_mesh_t structure.
 *---------------------------------------------------------------------------*/

cs_join_mesh_t *
cs_join_mesh_create(const char  *name)
{
  cs_join_mesh_t  *new_mesh = NULL;

  BFT_MALLOC(new_mesh, 1, cs_join_mesh_t);

  if (name != NULL) {

    int  len = strlen(name);

    BFT_MALLOC(new_mesh->name, len + 1, char);
    strncpy(new_mesh->name, name, len);
    new_mesh->name[len] = '\0';

  }
  else
    new_mesh->name = NULL;

  new_mesh->n_faces = 0;
  new_mesh->n_g_faces = 0;
  new_mesh->face_gnum = NULL;
  new_mesh->face_vtx_idx = NULL;
  new_mesh->face_vtx_lst = NULL;
  new_mesh->n_vertices = 0;
  new_mesh->n_g_vertices = 0;
  new_mesh->vertices = NULL;

  return new_mesh;
}

/*----------------------------------------------------------------------------
 * Get a cs_join_mesh_t structure with the given list of global faces inside.
 *
 * Exchange between ranks to get the connectivity associated to each
 * face of the global numbering list.
 *
 * parameters:
 *   mesh_name       <-- name of the created mesh
 *   n_elts          <-- number of elements in the global list
 *   glob_sel        <-- list of global elements (ordered)
 *   gnum_rank_index <-- index on ranks for the global elements
 *   local_mesh      <-- pointer to the local part of the distributed
 *                       cs_join_mesh_t structure on selected elements
 *
 * returns:
 *   a pointer to a new allocated cs_join_mesh_t structure
 *---------------------------------------------------------------------------*/

cs_join_mesh_t *
cs_join_mesh_create_from_glob_sel(const char            *mesh_name,
                                  cs_int_t               n_elts,
                                  const fvm_gnum_t       glob_sel[],
                                  const fvm_gnum_t       gnum_rank_index[],
                                  const cs_join_mesh_t  *local_mesh)
{
  cs_join_mesh_t  *new_mesh = NULL;

  const int  n_ranks = cs_glob_n_ranks;

  if (n_ranks == 1) {

    cs_int_t  i;
    cs_int_t *loc_sel = NULL;

    BFT_MALLOC(loc_sel, n_elts, cs_int_t);

    for (i = 0; i < n_elts; i++)
      loc_sel[i] = glob_sel[i];

    new_mesh =  cs_join_mesh_create_from_subset(mesh_name,
                                                n_elts,
                                                loc_sel,
                                                local_mesh);

    BFT_FREE(loc_sel);
  }

#if defined(HAVE_MPI)

  else { /* Parallel mode */

    cs_int_t  *send_rank_index = NULL, *send_faces = NULL;

    new_mesh = cs_join_mesh_create(mesh_name);

    /* Define a send list (face ids to send) from the global list of faces
       to receive. */

    _get_send_faces(n_ranks,
                    gnum_rank_index,
                    n_elts,
                    glob_sel,
                    &send_rank_index,
                    &send_faces);

    /* Get useful connectivity on ranks for computing local intersections */

    cs_join_mesh_exchange(n_ranks,
                          send_rank_index,
                          send_faces,
                          local_mesh,
                          new_mesh,
                          cs_glob_mpi_comm);

    BFT_FREE(send_faces);
    BFT_FREE(send_rank_index);
  }

#endif

  return new_mesh;
}

/*----------------------------------------------------------------------------
 * Allocate and define a cs_join_mesh_t structure relative to an extraction
 * of selected faces.
 *
 * The selection must be ordered.
 *
 * parameters:
 *   mesh_name   <-- name of the name to create
 *   subset_size <-- number of selected faces in the subset
 *   selection   <-- list of selected faces. Numbering in parent mesh
 *   parent_mesh <-- parent cs_join_mesh_t structure
 *
 * returns:
 *   a pointer to a cs_join_mesh_t structure
 *---------------------------------------------------------------------------*/

cs_join_mesh_t *
cs_join_mesh_create_from_subset(const char            *mesh_name,
                                cs_int_t               subset_size,
                                const cs_int_t         selection[],
                                const cs_join_mesh_t  *parent_mesh)
{
  cs_int_t  i, j, shift, id, parent_id, vtx_num, start, end;

  cs_int_t  n_select_vertices = 0;
  cs_int_t  *select_vertices = NULL, *counter = NULL;

  cs_join_mesh_t  *mesh = NULL;

  /* Get the selected vertices relative to the subset selection */

  BFT_MALLOC(counter, parent_mesh->n_vertices, cs_int_t);

  for (i = 0; i < parent_mesh->n_vertices; i++)
    counter[i] = 0;

  for (i = 0; i < subset_size; i++) {

    parent_id = selection[i] - 1;

    for (j = parent_mesh->face_vtx_idx[parent_id] - 1;
         j < parent_mesh->face_vtx_idx[parent_id+1] - 1; j++)
      counter[parent_mesh->face_vtx_lst[j]-1] += 1;

  }

  for (i = 0; i < parent_mesh->n_vertices; i++)
    if (counter[i] > 0)
      n_select_vertices++;

  BFT_MALLOC(select_vertices, n_select_vertices, cs_int_t);

  n_select_vertices = 0;
  for (i = 0; i < parent_mesh->n_vertices; i++)
    if (counter[i] > 0)
      select_vertices[n_select_vertices++] = i+1;

  /* Create a new cs_join_mesh_t structure */

  mesh = cs_join_mesh_create(mesh_name);

  mesh->n_faces = subset_size;

  /* Build face_vtx_idx */

  BFT_MALLOC(mesh->face_vtx_idx, mesh->n_faces + 1, cs_int_t);

  for (i = 0; i < mesh->n_faces; i++) {

    parent_id = selection[i] - 1;
    mesh->face_vtx_idx[i+1] =  parent_mesh->face_vtx_idx[parent_id+1]
                             - parent_mesh->face_vtx_idx[parent_id];

  }

  mesh->face_vtx_idx[0] = 1;
  for (i = 0; i < mesh->n_faces; i++)
    mesh->face_vtx_idx[i+1] += mesh->face_vtx_idx[i];

  BFT_MALLOC(mesh->face_vtx_lst,
             mesh->face_vtx_idx[mesh->n_faces] - 1, cs_int_t);

  /* Build face_vtx_lst */

  for (i = 0; i < mesh->n_faces; i++) {

    parent_id = selection[i] - 1;
    start = parent_mesh->face_vtx_idx[parent_id] - 1;
    end = parent_mesh->face_vtx_idx[parent_id+1] - 1;
    shift = mesh->face_vtx_idx[i] - 1;

    for (j = 0; j < end - start; j++) {

      vtx_num = parent_mesh->face_vtx_lst[start + j];

      id = cs_search_binary(0,
                            n_select_vertices,
                            vtx_num,
                            select_vertices);

      mesh->face_vtx_lst[shift + j] = id + 1;

    }

  } /* End of loop on selected faces */

 /* Define vertices */

  mesh->n_vertices = n_select_vertices;

  BFT_MALLOC(mesh->vertices, n_select_vertices, cs_join_vertex_t);

  for (i = 0; i < n_select_vertices; i++)
    mesh->vertices[i] = parent_mesh->vertices[select_vertices[i] - 1];

  /* Define global face numbering and linked global cell numbering */

  BFT_MALLOC(mesh->face_gnum, mesh->n_faces, fvm_gnum_t);

  for (i = 0; i < mesh->n_faces; i++) {

    parent_id = selection[i] - 1;
    mesh->face_gnum[i] = parent_mesh->face_gnum[parent_id];

  }

  if (cs_glob_n_ranks == 1) {

    mesh->n_g_faces = mesh->n_faces;
    mesh->n_g_vertices = mesh->n_vertices;

  }
  else {

    fvm_io_num_t  *io_num = NULL;
    fvm_gnum_t  *vtx_gnum = NULL;

    const fvm_gnum_t  *io_gnum = NULL;

    /* Get the global number of faces in the subset */

    io_num = fvm_io_num_create(selection, mesh->face_gnum, subset_size, 0);

    mesh->n_g_faces = fvm_io_num_get_global_count(io_num);

    io_num = fvm_io_num_destroy(io_num);

    /* Get the global number of vertices in the subset */

    BFT_MALLOC(vtx_gnum, mesh->n_vertices, fvm_gnum_t);

    for (i = 0; i < mesh->n_vertices; i++)
      vtx_gnum[i] = mesh->vertices[i].gnum;

    io_num = fvm_io_num_create(NULL, vtx_gnum, mesh->n_vertices, 0);

    mesh->n_g_vertices = fvm_io_num_get_global_count(io_num);

    io_gnum = fvm_io_num_get_global_num(io_num);

    for (i = 0; i < mesh->n_vertices; i++)
      mesh->vertices[i].gnum = io_gnum[i];

    io_num = fvm_io_num_destroy(io_num);

    BFT_FREE(vtx_gnum);

  }

  /* Free memory */

  BFT_FREE(counter);
  BFT_FREE(select_vertices);

  return  mesh;
}

/*----------------------------------------------------------------------------
 * Allocate and define a cs_join_mesh_t structure relative to an extraction
 * of selected faces.
 *
 * parameters:
 *   mesh_name           <-- name of the mesh to create
 *   n_faces             <-- number of selected faces
 *   n_g_faces           <-- global number of selected faces
 *   selected_faces      <-- list of selected faces (based on parent mesh)
 *   face_gnum           <-- global face numbers
 *   face_vtx_idx        <-- "face -> vertex" connectivity index
 *   face_vtx_lst        <-- "face -> vertex" connectivity
 *   n_select_vertices   <-- number of vertices used in joining operation
 *   n_g_select_vertices <-- global number of vertices used in joining
 *   select_vertices     <-- list of selected vertices (based on parent mesh)
 *   vtx_data            <-- array on data associated to selected vertices
 *
 * returns:
 *   a pointer to a cs_join_mesh_t structure
 *---------------------------------------------------------------------------*/

cs_join_mesh_t *
cs_join_mesh_create_from_extract(const char              *mesh_name,
                                 cs_int_t                 n_faces,
                                 fvm_gnum_t               n_g_faces,
                                 const cs_int_t           selected_faces[],
                                 const fvm_gnum_t         face_gnum[],
                                 const cs_int_t           face_vtx_idx[],
                                 const cs_int_t           face_vtx_lst[],
                                 cs_int_t                 n_select_vertices,
                                 fvm_gnum_t               n_g_select_vertices,
                                 const cs_int_t           selected_vertices[],
                                 const cs_join_vertex_t  *vtx_data)
{
  cs_int_t  i, j, shift, id, face_id, vtx_num, start, end;

  cs_join_mesh_t  *mesh = NULL;

  mesh = cs_join_mesh_create(mesh_name);

  /* Define face connectivity */

  mesh->n_faces = n_faces;
  mesh->n_g_faces = n_g_faces;

  /* Define face_vtx_idx */

  BFT_MALLOC(mesh->face_vtx_idx, n_faces + 1, cs_int_t);

  for (i = 0; i < n_faces; i++) {

    face_id = selected_faces[i] - 1;
    mesh->face_vtx_idx[i+1] = face_vtx_idx[face_id+1] - face_vtx_idx[face_id];

  }

  mesh->face_vtx_idx[0] = 1;
  for (i = 0; i < n_faces; i++)
    mesh->face_vtx_idx[i+1] += mesh->face_vtx_idx[i];

  BFT_MALLOC(mesh->face_vtx_lst, mesh->face_vtx_idx[n_faces] - 1, cs_int_t);

  /* Define face_vtx_lst */

  for (i = 0; i < n_faces; i++) {

    face_id = selected_faces[i] - 1;
    start = face_vtx_idx[face_id] - 1;
    end = face_vtx_idx[face_id+1] - 1;
    shift = mesh->face_vtx_idx[i] - 1;

    for (j = 0; j < end - start; j++) {

      vtx_num = face_vtx_lst[start + j];

      id = cs_search_binary(0,
                            n_select_vertices,
                            vtx_num,
                            selected_vertices);

      mesh->face_vtx_lst[shift + j] = id + 1;

    }

  } /* End of loop on selected faces */

  /* Define global face numbering */

  BFT_MALLOC(mesh->face_gnum, mesh->n_faces, fvm_gnum_t);

  if (face_gnum == NULL)
    for (i = 0; i < n_faces; i++)
      mesh->face_gnum[i] = selected_faces[i];
  else
    for (i = 0; i < n_faces; i++)
      mesh->face_gnum[i] = face_gnum[i];

  /* Define vertices */

  mesh->n_vertices = n_select_vertices;
  mesh->n_g_vertices = n_g_select_vertices;

  BFT_MALLOC(mesh->vertices, n_select_vertices, cs_join_vertex_t);

  for (i = 0; i < n_select_vertices; i++)
    mesh->vertices[i] = vtx_data[i];

  return  mesh;
}

/*----------------------------------------------------------------------------
 * Destroy a cs_join_mesh_t structure.
 *
 * parameters:
 *  mesh <->  pointer to pointer to cs_join_mesh_t structure to destroy
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_destroy(cs_join_mesh_t  **mesh)
{
  if (*mesh != NULL) {
    cs_join_mesh_t *m = *mesh;
    BFT_FREE(m->name);
    BFT_FREE(m->face_vtx_idx);
    BFT_FREE(m->face_vtx_lst);
    BFT_FREE(m->face_gnum);
    BFT_FREE(m->vertices);
    BFT_FREE(*mesh);
  }
}

/*----------------------------------------------------------------------------
 * Re-initialize an existing cs_join_mesh_t structure.
 *
 * parameters:
 *   mesh <-> pointer to a cs_join_mesh_t structure
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_reset(cs_join_mesh_t  *mesh)
{
  if (mesh == NULL)
    return;

  mesh->n_faces = 0;
  mesh->n_g_faces = 0;

  BFT_FREE(mesh->face_gnum);
  BFT_FREE(mesh->face_vtx_lst);
  BFT_FREE(mesh->face_vtx_idx);

  mesh->n_vertices = 0;
  mesh->n_g_vertices = 0;

  BFT_FREE(mesh->vertices);
}

/*----------------------------------------------------------------------------
 * Copy a cs_join_mesh_t structure into another.
 *
 * parameters:
 *   mesh     <-> pointer to a cs_join_mesh_t structure to fill
 *   ref_mesh <-- pointer to the reference
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_copy(cs_join_mesh_t        **mesh,
                  const cs_join_mesh_t   *ref_mesh)
{
  cs_int_t  i;
  cs_join_mesh_t  *_mesh = *mesh;

  if (ref_mesh == NULL) {
    cs_join_mesh_destroy(mesh);
    return;
  }

  if (_mesh == NULL)
    _mesh = cs_join_mesh_create(ref_mesh->name);

  _mesh->n_faces = ref_mesh->n_faces;
  _mesh->n_g_faces = ref_mesh->n_g_faces;

  BFT_REALLOC(_mesh->face_gnum, _mesh->n_faces, fvm_gnum_t);
  BFT_REALLOC(_mesh->face_vtx_idx, _mesh->n_faces + 1, cs_int_t);

  _mesh->face_vtx_idx[0] = 1;

  for (i = 0; i < _mesh->n_faces; i++) {
    _mesh->face_gnum[i] = ref_mesh->face_gnum[i];
    _mesh->face_vtx_idx[i+1] = ref_mesh->face_vtx_idx[i+1];
  }

  BFT_REALLOC(_mesh->face_vtx_lst,
              _mesh->face_vtx_idx[_mesh->n_faces]-1,
              cs_int_t);

  for (i = 0; i < _mesh->face_vtx_idx[_mesh->n_faces]-1; i++)
    _mesh->face_vtx_lst[i] = ref_mesh->face_vtx_lst[i];

  _mesh->n_vertices = ref_mesh->n_vertices;
  _mesh->n_g_vertices = ref_mesh->n_g_vertices;

  BFT_REALLOC(_mesh->vertices, _mesh->n_vertices, cs_join_vertex_t);

  memcpy(_mesh->vertices,
         ref_mesh->vertices,
         _mesh->n_vertices*sizeof(cs_join_vertex_t));

  /* Set return pointer */

  *mesh = _mesh;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Get the connectivity of a list of global elements distributed over the
 * ranks.
 *
 * parameters:
 *   n_ranks         <-- number of ranks in the MPI communicator
 *   send_rank_index <-- index on ranks for the face distribution
 *   send_faces      <-- list of face ids to send
 *   send_mesh       <-- pointer to the sending cs_join_mesh_t structure
 *   recv_mesh       <-> pointer to the receiving cs_join_mesh_t structure
 *   comm            <-- mpi communicator on which take places comm.
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_exchange(int                    n_ranks,
                      const cs_int_t         send_rank_index[],
                      const cs_int_t         send_faces[],
                      const cs_join_mesh_t  *send_mesh,
                      cs_join_mesh_t        *recv_mesh,
                      MPI_Comm               comm)
{
  int  i, j, rank, shift, start, end, face_id, vtx_id, vtx_count;
  int  local_rank;

  cs_int_t  n_face_to_recv = 0, n_vertices = 0, vtx_tag_size = 0;
  cs_int_t  *vtx_shift = NULL, *vtx_tag = NULL;
  cs_int_t  *send_count = NULL, *recv_count = NULL;
  cs_int_t  *send_shift = NULL, *recv_shift = NULL;
  fvm_gnum_t  *send_gbuf = NULL, *recv_gbuf = NULL;
  cs_join_vertex_t  *send_vtx_buf = NULL, *recv_vtx_buf = NULL;

  MPI_Datatype  MPI_JOIN_VERTEX = cs_join_mesh_create_vtx_datatype();

  /* Sanity checks */

#if defined(DEBUG) && !defined(NDEBUG)
  int  n_verif_ranks;

  MPI_Comm_size(comm, &n_verif_ranks);

  assert(n_ranks == n_verif_ranks);
#endif

  assert(send_mesh != NULL);
  assert(recv_mesh != NULL);
  assert(send_rank_index != NULL);
  assert(n_ranks > 1);

  MPI_Comm_rank(comm, &local_rank);

  /* Count the number of faces to recv */

  BFT_MALLOC(send_count, n_ranks, cs_int_t);
  BFT_MALLOC(recv_count, n_ranks, cs_int_t);

  for (i = 0; i < n_ranks; i++)
    send_count[i] = send_rank_index[i+1] - send_rank_index[i];

  /* Exchange number of elements to send */

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

  for (i = 0; i < n_ranks; i++)
    n_face_to_recv += recv_count[i];

  /* Update cs_join_mesh_t structure */

  recv_mesh->n_faces = n_face_to_recv;

  BFT_MALLOC(recv_mesh->face_gnum, n_face_to_recv, fvm_gnum_t);
  BFT_MALLOC(recv_mesh->face_vtx_idx, n_face_to_recv + 1, cs_int_t);

  /*
     The mesh doesn't change from a global point of view.
     It's only a redistribution of the elements according to the send_faces
     list.
  */

  recv_mesh->n_g_faces = send_mesh->n_g_faces;
  recv_mesh->n_g_vertices = send_mesh->n_g_vertices;

  /* Exchange face connect. count */

  BFT_MALLOC(vtx_tag, send_mesh->n_vertices, cs_int_t);
  BFT_MALLOC(vtx_shift, n_ranks+1, cs_int_t);

  vtx_shift[0] = 0;

  for (i = 0; i < n_ranks; i++) {
    vtx_shift[i+1] = 0;
    send_count[i] = 0;
  }

  for (rank = 0; rank < n_ranks; rank++) {

    for (i = 0; i < send_mesh->n_vertices; i++)
      vtx_tag[i] = -1;

    for (i = send_rank_index[rank]; i < send_rank_index[rank+1]; i++) {

      face_id = send_faces[i];
      start = send_mesh->face_vtx_idx[face_id] - 1;
      end = send_mesh->face_vtx_idx[face_id+1] - 1;
      n_vertices = end - start;

      for (j = start; j < end; j++) {

        vtx_id = send_mesh->face_vtx_lst[j] - 1;

        if (vtx_tag[vtx_id] < 0) {
          vtx_tag[vtx_id] = 1;
          vtx_shift[rank+1] += 1;
        }

      }

      send_count[rank] +=  1            /* face_gnum     */
                         + 1            /* n_vertices    */
                         + n_vertices;  /* face connect. */

    }

  } /* End of loop on ranks */

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

  BFT_MALLOC(send_shift, n_ranks + 1, cs_int_t);
  BFT_MALLOC(recv_shift, n_ranks + 1, cs_int_t);

  /* Build index arrays */

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (rank = 0; rank < n_ranks; rank++) {

    recv_shift[rank+1] = recv_shift[rank] + recv_count[rank];
    send_shift[rank+1] = send_shift[rank] + send_count[rank];
    vtx_shift[rank+1] += vtx_shift[rank];

  }

  /* Build send_gbuf to exchange face connectivity */

  BFT_MALLOC(send_gbuf, send_shift[n_ranks], fvm_gnum_t);
  BFT_MALLOC(recv_gbuf, recv_shift[n_ranks], fvm_gnum_t);

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (rank = 0; rank < n_ranks; rank++) {

    vtx_count = 0;

    for (i = 0; i < send_mesh->n_vertices; i++)
      vtx_tag[i] = -1;

    for (i = send_rank_index[rank]; i < send_rank_index[rank + 1]; i++) {

      shift = send_shift[rank] + send_count[rank];
      face_id = send_faces[i];

      start = send_mesh->face_vtx_idx[face_id] - 1;
      end = send_mesh->face_vtx_idx[face_id+1] - 1;
      n_vertices = end - start;

      send_gbuf[shift++] = send_mesh->face_gnum[face_id];
      send_gbuf[shift++] = n_vertices;

      for (j = start; j < end; j++) {

        vtx_id = send_mesh->face_vtx_lst[j] - 1;

        if (vtx_tag[vtx_id] < 0)
          vtx_tag[vtx_id] = vtx_count++;

        send_gbuf[shift++] = vtx_tag[vtx_id];

      }

      send_count[rank] +=  1            /* face_gnum     */
                         + 1            /* n_vertices    */
                         + n_vertices;  /* face connect. */

    }

  } /* End of loop on ranks */

  MPI_Alltoallv(send_gbuf, send_count, send_shift, FVM_MPI_GNUM,
                recv_gbuf, recv_count, recv_shift, FVM_MPI_GNUM, comm);

  BFT_FREE(send_gbuf);

  /* Scan recv_gbuf to build face->vertex connect. index */

  shift = 0;
  face_id = 0;

  while (shift < recv_shift[n_ranks]) {

    recv_mesh->face_gnum[face_id] = recv_gbuf[shift++];
    n_vertices = recv_gbuf[shift++];
    recv_mesh->face_vtx_idx[face_id+1] = n_vertices;

    face_id++;
    shift += n_vertices;

  }

  assert(n_face_to_recv == face_id);
  assert(shift == recv_shift[n_ranks]);

  recv_mesh->face_vtx_idx[0] = 1;
  for (i = 0; i < recv_mesh->n_faces; i++)
    recv_mesh->face_vtx_idx[i+1] += recv_mesh->face_vtx_idx[i];

  /* Scan recv_gbuf to build face->vertex connectivity list */

  BFT_MALLOC(recv_mesh->face_vtx_lst,
             recv_mesh->face_vtx_idx[n_face_to_recv], cs_int_t);

  vtx_tag_size = send_mesh->n_vertices;

  if (recv_mesh->face_vtx_idx[n_face_to_recv] > send_mesh->n_vertices) {

    vtx_tag_size = recv_mesh->face_vtx_idx[n_face_to_recv];
    BFT_REALLOC(vtx_tag, recv_mesh->face_vtx_idx[n_face_to_recv], cs_int_t);

  }

  /* Store vtx_shift data into send_count in order to re-use it */

  for (rank = 0; rank < n_ranks; rank++)
    send_count[rank] = vtx_shift[rank+1] - vtx_shift[rank];

  for (rank = 0; rank < n_ranks + 1; rank++)
    vtx_shift[rank] = 0;

  shift = 0;
  face_id = 0;

  for (rank = 0; rank < n_ranks; rank++) {

    vtx_count = 0;

    for (i = 0; i < vtx_tag_size; i++)
      vtx_tag[i] = -1;

    while (shift < recv_shift[rank + 1]) {

      shift += 1; /* skip face_gnum */
      n_vertices = recv_gbuf[shift++];
      start = recv_mesh->face_vtx_idx[face_id] - 1;

      for (j = 0; j < n_vertices; j++) {

        vtx_id = recv_gbuf[shift++];

        if (vtx_tag[vtx_id] < 0) {
          vtx_count++;
          vtx_tag[vtx_id] = vtx_shift[rank] + vtx_count;
        }

        recv_mesh->face_vtx_lst[start + j] = vtx_tag[vtx_id];

      }

      face_id++;

    } /* Scan recv_gbuf for the current rank */

    vtx_shift[rank+1] = vtx_shift[rank] + vtx_count;

  } /* End of loop on ranks */

  /* Exchange number of vertices to communicate */

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (rank = 0; rank < n_ranks; rank++) {
    send_shift[rank+1] = send_shift[rank] + send_count[rank];
    recv_shift[rank+1] = recv_shift[rank] + recv_count[rank];
  }

  /* Partial memory management */

  BFT_FREE(vtx_shift);
  BFT_FREE(recv_gbuf);

  BFT_MALLOC(send_vtx_buf, send_shift[n_ranks], cs_join_vertex_t);
  BFT_MALLOC(recv_vtx_buf, recv_shift[n_ranks], cs_join_vertex_t);

  /* Exchange vertex buffers */

  for (rank = 0; rank < n_ranks; rank++) {

    vtx_count = 0;

    for (i = 0; i < vtx_tag_size; i++)
      vtx_tag[i] = -1;

    for (i = send_rank_index[rank]; i < send_rank_index[rank + 1]; i++) {

      face_id = send_faces[i];
      start = send_mesh->face_vtx_idx[face_id] - 1;
      end = send_mesh->face_vtx_idx[face_id+1] - 1;
      n_vertices = end - start;

      for (j = start; j < end; j++) {

        vtx_id = send_mesh->face_vtx_lst[j] - 1;

        if (vtx_tag[vtx_id] < 0) { /* add the vertex to send_vtx_buf */

          shift = send_shift[rank] + vtx_count;
          vtx_tag[vtx_id] = vtx_count++;
          send_vtx_buf[shift] = send_mesh->vertices[vtx_id];

        }

      } /* End of loop on the face connectivity */

    } /* End of loop on faces to send to the current rank */

  } /* End of loop on ranks */

  MPI_Alltoallv(send_vtx_buf, send_count, send_shift, MPI_JOIN_VERTEX,
                recv_vtx_buf, recv_count, recv_shift, MPI_JOIN_VERTEX,
                comm);

  /* Update cs_join_mesh_t structure */

  recv_mesh->n_vertices = recv_shift[n_ranks];

  BFT_MALLOC(recv_mesh->vertices, recv_shift[n_ranks], cs_join_vertex_t);

  memcpy(recv_mesh->vertices,
         recv_vtx_buf,
         recv_shift[n_ranks]*sizeof(cs_join_vertex_t));

  /* Delete vertices which appear several times */

  cs_join_mesh_vertex_clean(recv_mesh);

  /* Free memory */

  MPI_Type_free(&MPI_JOIN_VERTEX);
  BFT_FREE(vtx_tag);
  BFT_FREE(send_count);
  BFT_FREE(recv_count);
  BFT_FREE(send_shift);
  BFT_FREE(recv_shift);
  BFT_FREE(send_vtx_buf);
  BFT_FREE(recv_vtx_buf);

}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Destroy a cs_join_edges_t structure.
 *
 * parameters:
 *   edges <->  pointer to pointer to cs_join_edges_t structure to destroy
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_destroy_edges(cs_join_edges_t  **edges)
{
  if (*edges != NULL) {

    cs_join_edges_t  *e = *edges;

    if (e->n_edges > 0) {

      BFT_FREE(e->gnum);
      BFT_FREE(e->def);

      BFT_FREE(e->vtx_idx);
      BFT_FREE(e->adj_vtx_lst);
      BFT_FREE(e->edge_lst);
    }

    BFT_FREE(*edges);
  }
}

/*----------------------------------------------------------------------------
 * Order a cs_join_mesh_t structure according to the global face numbering
 *
 * Delete redundancies.
 *
 * parameters:
 *   mesh <-> pointer to a cs_join_mesh_t structure to order
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_face_order(cs_join_mesh_t  *mesh)
{
  int  i, j, o_id;
  cs_int_t  shift, start, end, n_new_faces;
  fvm_gnum_t  prev, cur;

  cs_int_t  n_faces = mesh->n_faces;
  cs_int_t  *num_buf = NULL,  *selection = NULL;
  fvm_lnum_t  *order = NULL;
  fvm_gnum_t  *gnum_buf = NULL;

  assert(mesh != NULL);

  if (n_faces == 0)
    return;

  /* Order faces according to their global numbering */

  BFT_MALLOC(order, n_faces, fvm_lnum_t);

  fvm_order_local_allocated(NULL, mesh->face_gnum, order, n_faces);

  /* Order global face numbering */

  BFT_MALLOC(gnum_buf, n_faces, fvm_gnum_t);
  BFT_MALLOC(selection, n_faces, cs_int_t);

  for (i = 0; i < n_faces; i++)
    gnum_buf[i] = mesh->face_gnum[i];

  prev = 0;
  n_new_faces = 0;

  for (i = 0; i < n_faces; i++) {

    o_id = order[i];
    cur = gnum_buf[o_id];

    if (prev != cur) {
      prev = cur;
      selection[n_new_faces] = o_id;
      mesh->face_gnum[n_new_faces] = cur;
      n_new_faces++;
    }

  }

  mesh->n_faces = n_new_faces;

  BFT_FREE(gnum_buf);
  BFT_FREE(order);

  BFT_REALLOC(mesh->face_gnum, n_new_faces, fvm_gnum_t);
  BFT_REALLOC(selection, n_new_faces, cs_int_t);

  /* Order face -> vertex connectivity list */

  BFT_MALLOC(num_buf, mesh->face_vtx_idx[n_faces], cs_int_t);

  for (i = 0; i < mesh->face_vtx_idx[n_faces] - 1; i++)
    num_buf[i] = mesh->face_vtx_lst[i];

  shift = 0;

  for (i = 0; i < n_new_faces; i++) {

    o_id = selection[i];
    start = mesh->face_vtx_idx[o_id] - 1;
    end = mesh->face_vtx_idx[o_id+1] - 1;

    for (j = start; j < end; j++)
      mesh->face_vtx_lst[shift++] = num_buf[j];

  } /* End of loop on faces */

  BFT_REALLOC(num_buf, n_faces, cs_int_t);

  for (i = 0; i < n_faces; i++)
    num_buf[i] = mesh->face_vtx_idx[i+1] - mesh->face_vtx_idx[i];

  for (i = 0; i < n_new_faces; i++) {
    o_id = selection[i];
    mesh->face_vtx_idx[i+1] = mesh->face_vtx_idx[i] + num_buf[o_id];
  }

  /* Memory management */

  BFT_FREE(selection);
  BFT_FREE(num_buf);
  BFT_REALLOC(mesh->face_vtx_idx, n_new_faces+1, cs_int_t);
  BFT_REALLOC(mesh->face_vtx_lst, mesh->face_vtx_idx[n_new_faces], cs_int_t);
}

/*----------------------------------------------------------------------------
 * Delete vertices which appear several times (same global number) and
 * vertices which are not used in face definition.
 *
 * parameters:
 *   mesh <-> pointer to cs_join_mesh_t structure to clean
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_vertex_clean(cs_join_mesh_t  *mesh)
{
  cs_int_t  i, j, shift, n_init_vertices, n_final_vertices;
  fvm_gnum_t  prev, cur;

  fvm_lnum_t  *order = NULL;
  cs_int_t  *init2final = NULL, *tag = NULL;
  fvm_gnum_t  *gnum_buf = NULL;
  cs_join_vertex_t  *final_vertices = NULL;

  assert(mesh != NULL);

  n_init_vertices = mesh->n_vertices;

  if (n_init_vertices < 2)
    return;

  /* Count the final number of vertices */

  BFT_MALLOC(order, n_init_vertices, fvm_lnum_t);
  BFT_MALLOC(tag, n_init_vertices, cs_int_t);
  BFT_MALLOC(gnum_buf, n_init_vertices, fvm_gnum_t);

  for (i = 0; i < n_init_vertices; i++) {
    gnum_buf[i] = mesh->vertices[i].gnum;
    tag[i] = 0;
  }

  /* Tag vertices really used in the mesh definition */

  for (i = 0; i < mesh->n_faces; i++)
    for (j = mesh->face_vtx_idx[i] - 1; j < mesh->face_vtx_idx[i+1] - 1; j++)
      tag[mesh->face_vtx_lst[j] - 1] = 1;

  /* Order vertices by increasing global number */

  fvm_order_local_allocated(NULL, gnum_buf, order, n_init_vertices);

  n_final_vertices = 0;
  prev = 0;

  for (i = 0; i < n_init_vertices; i++) {

    shift = order[i];
    cur = gnum_buf[shift];

    if (prev != cur && tag[i] > 0) {
      n_final_vertices++;
      prev = cur;
    }

  }

  /* Define the final vertices structure and indirection between
     initial numbering and final numbering */

  BFT_MALLOC(final_vertices, n_final_vertices, cs_join_vertex_t);
  BFT_MALLOC(init2final, n_init_vertices, cs_int_t);

  n_final_vertices = 0;
  prev = 0;

  for (i = 0; i < n_init_vertices; i++) {

    shift = order[i];
    cur = gnum_buf[shift];

    if (prev != cur  && tag[i] > 0) {

      final_vertices[n_final_vertices++] = mesh->vertices[shift];
      prev = cur;

    }

    init2final[shift] = n_final_vertices;

  }

  BFT_FREE(mesh->vertices);

  mesh->vertices = final_vertices;
  mesh->n_vertices = n_final_vertices;

  /* Update face->vertex connectivity list */

  for (i = 0; i < mesh->n_faces; i++) {

    for (j = mesh->face_vtx_idx[i] - 1; j < mesh->face_vtx_idx[i+1] - 1; j++)
      mesh->face_vtx_lst[j] = init2final[mesh->face_vtx_lst[j]-1];

  } /* end of loop on faces */

  BFT_FREE(init2final);
  BFT_FREE(gnum_buf);
  BFT_FREE(tag);
  BFT_FREE(order);
}

/*----------------------------------------------------------------------------
 * Clean the given cs_join_mesh_t structure, removing degenerate edges.
 *
 * parameters:
 *   mesh      <-> pointer to the cs_join_mesh_t structure to clean
 *   verbosity <-- level of display
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_clean(cs_join_mesh_t  *mesh,
                   int              verbosity)
{

  /* Delete empty edge:
       These edges are generated during the merge step. If two vertices
       sharing the same edge are fused, we have to delete the resulting
       edge.
  */

  _remove_empty_edges(mesh, verbosity);

  /* Delete degenerate edge:
       These edges are generated during the merge step.


       x                                      x
       |\                                     |
       | \                                    |
     a2|  \a3                               A2|
       |   \                Merge            |
       |    \      a4         of              |
    ---s1----s2------       vertices          x
       |      \            s1 and s2         / \
       |       \                            /   \
     a1|        \a4                      A1/     \A3
       |         \                        /       \
       |          \                      /         \
       x-----------x                    x-----------x
            a5                                A4

   */

  _remove_degenerate_edges(mesh, verbosity);

}

/*----------------------------------------------------------------------------
 * Define a list of edges associated to a cs_join_mesh_t structure.
 *
 * parameters:
 *   mesh <-- pointer to a cs_join_mesh_t structure
 *
 * returns:
 *   a pointer to the new defined cs_join_edges_t structure.
 *---------------------------------------------------------------------------*/

cs_join_edges_t *
cs_join_mesh_define_edges(const cs_join_mesh_t  *mesh)
{
  int  i, j;
  cs_int_t  v1_num, v2_num, o_id1, o_id2;
  cs_int_t  edge_shift, shift, n_init_edges;
  fvm_gnum_t  v1_gnum, v2_gnum;

  fvm_lnum_t  *order = NULL;
  cs_int_t  *vtx_counter = NULL, *vtx_lst = NULL;
  fvm_gnum_t  *adjacency = NULL;
  cs_join_edges_t  *edges = NULL;

  if (mesh == NULL)
    return edges;

  /* Initialization and structure allocation */

  BFT_MALLOC(edges, 1, cs_join_edges_t);

  edges->n_edges = 0;
  edges->def = NULL;
  edges->gnum = NULL;
  edges->n_vertices = mesh->n_vertices;
  edges->vtx_idx = NULL;
  edges->adj_vtx_lst = NULL;
  edges->edge_lst = NULL;

  /* Define edges */

  n_init_edges = mesh->face_vtx_idx[mesh->n_faces] - 1;

  if (n_init_edges == 0)
    return edges;

  BFT_MALLOC(edges->def, 2*n_init_edges, cs_int_t);
  BFT_MALLOC(edges->vtx_idx, mesh->n_vertices + 1, cs_int_t);

  for (i = 0; i < mesh->n_vertices + 1; i++)
    edges->vtx_idx[i] = 0;

  /* Loop on faces to initialize edge list */

  BFT_MALLOC(vtx_lst, 2*n_init_edges, cs_int_t);
  BFT_MALLOC(adjacency, 2*n_init_edges, fvm_gnum_t);

  for (shift = 0, i = 0; i < mesh->n_faces; i++) {

    cs_int_t  start = mesh->face_vtx_idx[i] - 1;
    cs_int_t  end =  mesh->face_vtx_idx[i+1] - 1;

    for (j = start; j < end - 1; j++) {

      v1_num = mesh->face_vtx_lst[j];
      v1_gnum = (mesh->vertices[v1_num-1]).gnum;
      v2_num = mesh->face_vtx_lst[j+1];
      v2_gnum = (mesh->vertices[v2_num-1]).gnum;

      if (v1_gnum > v2_gnum) {

        vtx_lst[2*shift] = v2_num;
        adjacency[2*shift] = v2_gnum;
        vtx_lst[2*shift+1] = v1_num;
        adjacency[2*shift+1] = v1_gnum;

      }
      else {

        vtx_lst[2*shift] = v1_num;
        adjacency[2*shift] = v1_gnum;
        vtx_lst[2*shift+1] = v2_num;
        adjacency[2*shift+1] = v2_gnum;

      }

      shift++;

    } /* End of loop on n-1 first vertices */

    v1_num = mesh->face_vtx_lst[end-1];
    v1_gnum = (mesh->vertices[v1_num-1]).gnum;
    v2_num = mesh->face_vtx_lst[start];
    v2_gnum = (mesh->vertices[v2_num-1]).gnum;

    if (v1_gnum > v2_gnum) {

      vtx_lst[2*shift] = v2_num;
      adjacency[2*shift] = v2_gnum;
      vtx_lst[2*shift+1] = v1_num;
      adjacency[2*shift+1] = v1_gnum;

    }
    else {

      vtx_lst[2*shift] = v1_num;
      adjacency[2*shift] = v1_gnum;
      vtx_lst[2*shift+1] = v2_num;
      adjacency[2*shift+1] = v2_gnum;

    }

    shift++;

  } /* End of loop on faces */

  assert(shift == n_init_edges);

  BFT_MALLOC(order, n_init_edges, fvm_lnum_t);

  fvm_order_local_allocated_s(NULL, adjacency, 2, order, n_init_edges);

  /* Fill cs_join_edges_t structure */

  o_id1 = order[0];
  edges->def[0] = vtx_lst[2*o_id1];
  edges->def[1] = vtx_lst[2*o_id1+1];
  edges->vtx_idx[vtx_lst[2*o_id1]] += 1;
  edges->vtx_idx[vtx_lst[2*o_id1+1]] += 1;
  edge_shift = 1;

  for (i = 1; i < n_init_edges; i++) {

    o_id1 = order[i-1];
    o_id2 = order[i];

    if (   vtx_lst[2*o_id1]   != vtx_lst[2*o_id2]
        || vtx_lst[2*o_id1+1] != vtx_lst[2*o_id2+1]) {

      edges->vtx_idx[vtx_lst[2*o_id2]] += 1;
      edges->vtx_idx[vtx_lst[2*o_id2+1]] += 1;
      edges->def[2*edge_shift] = vtx_lst[2*o_id2];
      edges->def[2*edge_shift+1] = vtx_lst[2*o_id2+1];
      edge_shift++;

    }

  } /* End of loop on edges */

  edges->n_edges = edge_shift;
  BFT_REALLOC(edges->def, 2*edges->n_edges, cs_int_t);

  /* Build adj_vtx_lst and edge_lst */

  BFT_MALLOC(vtx_counter, mesh->n_vertices, cs_int_t);

  for (i = 0; i < mesh->n_vertices; i++) {
    edges->vtx_idx[i+1] += edges->vtx_idx[i];
    vtx_counter[i] = 0;
  }

  BFT_MALLOC(edges->adj_vtx_lst, edges->vtx_idx[mesh->n_vertices], cs_int_t);
  BFT_MALLOC(edges->edge_lst, edges->vtx_idx[mesh->n_vertices], cs_int_t);

  {
    cs_int_t  vtx_id_a, vtx_id_b, shift_a, shift_b;
    fvm_gnum_t  vtx_gnum_a, vtx_gnum_b;

    cs_int_t  cur_edge_num = 1;

    /* Initiate edge_lst and adj_vtx_lst building */

    o_id1 = order[0];

    vtx_id_a = vtx_lst[2*o_id1]-1;
    vtx_id_b = vtx_lst[2*o_id1+1]-1;

    vtx_gnum_a = (mesh->vertices[vtx_id_a]).gnum;
    vtx_gnum_b = (mesh->vertices[vtx_id_a]).gnum;

    shift_a = edges->vtx_idx[vtx_id_a];
    shift_b = edges->vtx_idx[vtx_id_b];

    vtx_counter[vtx_id_a] += 1;
    vtx_counter[vtx_id_b] += 1;

    edges->adj_vtx_lst[shift_a] = vtx_id_b + 1;
    edges->adj_vtx_lst[shift_b] = vtx_id_a + 1;

    if (vtx_gnum_a > vtx_gnum_b) {
      edges->edge_lst[shift_a] = -cur_edge_num;
      edges->edge_lst[shift_b] = cur_edge_num;
    }
    else {
      edges->edge_lst[shift_a] = cur_edge_num;
      edges->edge_lst[shift_b] = -cur_edge_num;
    }

    cur_edge_num++;

    for (i = 1; i < n_init_edges; i++) {

      o_id1 = order[i-1];
      o_id2 = order[i];

      if (   vtx_lst[2*o_id1]   != vtx_lst[2*o_id2]
          || vtx_lst[2*o_id1+1] != vtx_lst[2*o_id2+1]) {

        vtx_id_a = vtx_lst[2*o_id2]-1;
        vtx_id_b = vtx_lst[2*o_id2+1]-1;

        vtx_gnum_a = (mesh->vertices[vtx_id_a]).gnum;
        vtx_gnum_b = (mesh->vertices[vtx_id_a]).gnum;

        shift_a = edges->vtx_idx[vtx_id_a] + vtx_counter[vtx_id_a];
        shift_b = edges->vtx_idx[vtx_id_b] + vtx_counter[vtx_id_b];

        vtx_counter[vtx_id_a] += 1;
        vtx_counter[vtx_id_b] += 1;

        edges->adj_vtx_lst[shift_a] = vtx_id_b + 1;
        edges->adj_vtx_lst[shift_b] = vtx_id_a + 1;

        if (vtx_gnum_a > vtx_gnum_b) {
          edges->edge_lst[shift_a] = -cur_edge_num;
          edges->edge_lst[shift_b] = cur_edge_num;
        }
        else {
          edges->edge_lst[shift_a] = cur_edge_num;
          edges->edge_lst[shift_b] = -cur_edge_num;
        }

        cur_edge_num++;

      }

    } /* End of loop on edges */

    assert(cur_edge_num - 1 == edges->n_edges);

  } /* End of adj_vtx_lst and edge_lst building */

  /* Partial clean-up */

  BFT_FREE(vtx_lst);
  BFT_FREE(vtx_counter);

  /* Define a global numbering on edges */

  BFT_MALLOC(edges->gnum, edges->n_edges, fvm_gnum_t);
  BFT_REALLOC(adjacency, 2*edges->n_edges, fvm_gnum_t);

  for (i = 0; i < edges->n_edges; i++) {

    cs_int_t  v1_id = edges->def[2*i] - 1;
    cs_int_t  v2_id = edges->def[2*i+1] - 1;

    v1_gnum = (mesh->vertices[v1_id]).gnum;
    v2_gnum = (mesh->vertices[v2_id]).gnum;

    if (v1_gnum > v2_gnum) {
      adjacency[2*i] = v2_gnum;
      adjacency[2*i+1] = v1_gnum;
    }
    else {
      adjacency[2*i] = v1_gnum;
      adjacency[2*i+1] = v2_gnum;
    }

  } /* End of loop on edges */

  /* Order vtx_lst and build an order list into adjacency */

  fvm_order_local_allocated_s(NULL, adjacency, 2, order, edges->n_edges);

  if (cs_glob_n_ranks == 1) { /* Serial treatment */

    edges->n_g_edges = edges->n_edges;

    for (i = 0; i < edges->n_edges; i++) {

      cs_int_t  o_id = order[i];

      edges->gnum[i] = o_id+1;

    }

  }
  else { /* Parallel treatment */

    fvm_gnum_t  *order_couples = NULL;
    fvm_io_num_t  *edge_io_num = NULL;
    const fvm_gnum_t  *edges_gnum = NULL;

    assert(cs_glob_n_ranks > 1);

    BFT_MALLOC(order_couples, 2*edges->n_edges, fvm_gnum_t);

    for (i = 0; i < edges->n_edges; i++) {

      cs_int_t  o_id = order[i];

      order_couples[2*i] = adjacency[2*o_id];
      order_couples[2*i+1] = adjacency[2*o_id+1];

    } /* End of loop on edges */

    edge_io_num = fvm_io_num_create_from_adj_s(NULL,
                                               order_couples,
                                               edges->n_edges,
                                               2);

    edges->n_g_edges = fvm_io_num_get_global_count(edge_io_num);
    edges_gnum = fvm_io_num_get_global_num(edge_io_num);

    for (i = 0; i < edges->n_edges; i++)
      edges->gnum[i] = edges_gnum[i];

    /* Partial Clean-up */

    BFT_FREE(order_couples);
    fvm_io_num_destroy(edge_io_num);

  }

  /* Memory management */

  BFT_FREE(adjacency);
  BFT_FREE(order);

  /* Return pointers */

  return edges;
}

/*----------------------------------------------------------------------------
 * Get the edge number relative to a couple of vertex numbers.
 *
 * edge_num > 0 if couple is in the same order as the edge->def
 * edge_num < 0 otherwise
 *
 * parameters:
 *   v1_num <-- vertex number for the first vertex
 *   v2_num <-- vertex number for the second vertex
 *   edges  <-- pointer to a cs_join_edges_t structure
 *
 * returns:
 *   an edge number relative to the couple of vertices
 *---------------------------------------------------------------------------*/

cs_int_t
cs_join_mesh_get_edge(cs_int_t                v1_num,
                      cs_int_t                v2_num,
                      const cs_join_edges_t  *edges)
{
  cs_int_t  i;
  cs_int_t  edge_num = 0;

  assert(edges != NULL);
  assert(v1_num > 0);
  assert(v2_num > 0);

  if (edges->vtx_idx[v1_num] - edges->vtx_idx[v1_num-1] == 0)
    bft_error(__FILE__, __LINE__, 0,
              _(" The given vertex number: %d is not defined"
                " in the edge structure (edges->vtx_idx).\n"), v1_num);

  for (i = edges->vtx_idx[v1_num-1]; i < edges->vtx_idx[v1_num]; i++) {
    if (edges->adj_vtx_lst[i] == v2_num) {
      edge_num = edges->edge_lst[i];
      break;
    }
  }

  if (edge_num == 0)
    bft_error(__FILE__, __LINE__, 0,
              _(" The given couple of vertex numbers :\n"
                "   vertex 1 : %d\n"
                "   vertex 2 : %d\n"
                " is not defined in the edge structure.\n"), v1_num, v2_num);

  assert(edge_num != 0);

  return edge_num;
}

/*----------------------------------------------------------------------------
 * Re-organize the cs_join_mesh_t structure after a renumbering of
 * the vertices following the merge operation + a new description of each
 * face.
 *
 * parameters:
 *   mesh             <-> pointer to the cs_join_mesh_t structure to update
 *   edges            <-- pointer to a cs_join_edges_t structure
 *   edge_index       <-- index on edges for the new vertices
 *   edge_new_vtx_lst <-- list of new vertices for each edge
 *   n_new_vertices   <-- new local number of vertices after merge
 *   old2new          <-- array storing the relation between old/new vertex id
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_update(cs_join_mesh_t         *mesh,
                    const cs_join_edges_t  *edges,
                    const cs_int_t          edge_index[],
                    const cs_int_t          edge_new_vtx_lst[],
                    cs_int_t                n_new_vertices,
                    const cs_int_t          old2new[])
{
  cs_int_t  i, j, n_adds;

  cs_join_vertex_t  *new_vertices = NULL;
  cs_int_t  *new_face_vtx_idx = NULL, *new_face_vtx_lst = NULL;

  /* Sanity checks */

  assert(mesh != NULL);
  assert(edges != NULL);

  /* Update description and numbering for face -> vertex connectivity */

  if (edge_new_vtx_lst != NULL) {

    BFT_MALLOC(new_face_vtx_idx, mesh->n_faces + 1, cs_int_t);

    for (i = 0; i < mesh->n_faces + 1; i++)
      new_face_vtx_idx[i] = 0;

    /* Update face -> vertex connectivity.
       add new vertices between the existing one */

    for (i = 0; i < mesh->n_faces; i++) {

      cs_int_t  start_id = mesh->face_vtx_idx[i] - 1;
      cs_int_t  end_id = mesh->face_vtx_idx[i+1] - 1;

      for (j = start_id; j < end_id - 1; j++) {

        n_adds = _count_new_added_vtx_to_edge(mesh->face_vtx_lst[j],
                                              mesh->face_vtx_lst[j+1],
                                              old2new,
                                              edges,
                                              edge_index,
                                              edge_new_vtx_lst);

        new_face_vtx_idx[i+1] += n_adds;

      }

      /* Case end - start */

      n_adds = _count_new_added_vtx_to_edge(mesh->face_vtx_lst[end_id-1],
                                            mesh->face_vtx_lst[start_id],
                                            old2new,
                                            edges,
                                            edge_index,
                                            edge_new_vtx_lst);

      new_face_vtx_idx[i+1] += n_adds;

    } /* End of loop on faces */

    /* Build new face_vtx_idx */

    new_face_vtx_idx[0] = 1;
    for (i = 0; i < mesh->n_faces; i++) {

      new_face_vtx_idx[i+1] += new_face_vtx_idx[i];

      if (new_face_vtx_idx[i+1] < 3)
        bft_error(__FILE__, __LINE__, 0,
                  _(" Problem in mesh connectivity."
                    " Face: %lu\n"
                    " Problem detected during connectivity update:\n"
                    " The face is defined by less than 3 points"
                    " (excessive merging has occured).\n\n"
                    " Modify joining parameters to reduce merging"
                    " (fraction & merge).\n"),
                  (unsigned long)(mesh->face_gnum[i]));

    }

    /* Build new_face_vtx_lst */

    BFT_MALLOC(new_face_vtx_lst, new_face_vtx_idx[mesh->n_faces], cs_int_t);

  } /* End if edge_new_vtx_lst != NULL */

  else { /* edge_new_vtx_lst == NULL */

    new_face_vtx_idx = mesh->face_vtx_idx;
    new_face_vtx_lst = mesh->face_vtx_lst;

  }

  for (i = 0; i < mesh->n_faces; i++) {

    cs_int_t  start_id = mesh->face_vtx_idx[i] - 1;
    cs_int_t  end_id = mesh->face_vtx_idx[i+1] - 1;
    cs_int_t  shift = new_face_vtx_idx[i] - 1;

    for (j = start_id; j < end_id-1; j++)
      _add_new_vtx_to_edge(mesh->face_vtx_lst[j],
                           mesh->face_vtx_lst[j+1],
                           old2new,
                           edges,
                           edge_index,
                           edge_new_vtx_lst,
                           new_face_vtx_lst,
                           &shift);

    /* Case end - start */

    _add_new_vtx_to_edge(mesh->face_vtx_lst[end_id-1],
                         mesh->face_vtx_lst[start_id],
                         old2new,
                         edges,
                         edge_index,
                         edge_new_vtx_lst,
                         new_face_vtx_lst,
                         &shift);

  } /* End of loop on faces */

  if (edge_new_vtx_lst != NULL) {

    BFT_FREE(mesh->face_vtx_idx);
    BFT_FREE(mesh->face_vtx_lst);

    mesh->face_vtx_idx = new_face_vtx_idx;
    mesh->face_vtx_lst = new_face_vtx_lst;
  }

  /* Define the new_vertices structure */

  BFT_MALLOC(new_vertices, n_new_vertices, cs_join_vertex_t);

  for (i = 0; i < mesh->n_vertices; i++)
    new_vertices[old2new[i]] = mesh->vertices[i];

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n\n Dump Old2New array : "
             "n_old_vertices = %d - n_new_vertices = %d\n",
             mesh->n_vertices, n_new_vertices);
  for (i = 0; i < mesh->n_vertices; i++)
    bft_printf("Old num : %7d (%9u) => New num : %7d (%9u)\n",
               i+1         , mesh->vertices[i].gnum,
               old2new[i]+1, new_vertices[old2new[i]].gnum);
  bft_printf_flush();
#endif

  /* Update mesh structure */

  BFT_FREE(mesh->vertices);

  mesh->n_vertices = n_new_vertices;
  mesh->n_g_vertices = n_new_vertices;
  mesh->vertices = new_vertices;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    fvm_gnum_t  *vtx_gnum = NULL;
    fvm_io_num_t  *io_num = NULL;

    /* Global number of selected vertices and associated
       fvm_io_num_t structure */

    BFT_MALLOC(vtx_gnum, n_new_vertices, fvm_gnum_t);

    for (i = 0; i < n_new_vertices; i++)
      vtx_gnum[i] = (mesh->vertices[i]).gnum;

    io_num = fvm_io_num_create(NULL, vtx_gnum, n_new_vertices, 0);

    mesh->n_g_vertices = fvm_io_num_get_global_count(io_num);

    fvm_io_num_destroy(io_num);

    BFT_FREE(vtx_gnum);

  }
#endif
}

/*----------------------------------------------------------------------------
 * Compute for each face of the cs_join_mesh_t structure the face normal.
 * || face_normal || = 1 (divided by the area of the face)
 *
 * The caller is responsible for freeing the returned array.
 *
 * parameters:
 *   mesh <-- pointer to a cs_join_mesh_t structure
 *
 *                          Pi+1
 *              *---------*                   B  : barycenter of the polygon
 *             / .       . \
 *            /   .     .   \                 Pi : vertices of the polygon
 *           /     .   .     \
 *          /       . .  Ti   \               Ti : triangle
 *         *.........B.........* Pi
 *     Pn-1 \       . .       /
 *           \     .   .     /
 *            \   .     .   /
 *             \ .   T0  . /
 *              *---------*
 *            P0
 *
 *
 * returns:
 *   an array with the face normal for each face of the mesh
 *---------------------------------------------------------------------------*/

cs_real_t *
cs_join_mesh_get_face_normal(const cs_join_mesh_t  *mesh)
{
  cs_int_t  i, j, k, vid;
  double  inv_norm;

  cs_int_t  n_max_vertices = 0;
  cs_real_t  *face_vtx_coord = NULL;
  cs_real_t  *face_normal = NULL;

  if (mesh == NULL)
    return face_normal;

  if (mesh->n_faces == 0)
    return face_normal;

  BFT_MALLOC(face_normal, 3*mesh->n_faces, cs_real_t);

  for (i = 0; i < 3*mesh->n_faces; i++)
    face_normal[i] = 0.0;

  /* Compute n_max_vertices */

  for (i = 0; i < mesh->n_faces; i++)
    n_max_vertices = CS_MAX(n_max_vertices,
                            mesh->face_vtx_idx[i+1] - mesh->face_vtx_idx[i]);

  BFT_MALLOC(face_vtx_coord, 3*(n_max_vertices+1), cs_real_t);

  for (i = 0; i < mesh->n_faces; i++) {

    cs_real_t  v1[3], v2[3], tri_normal[3];

    cs_int_t  shift = 0;
    cs_int_t  s = mesh->face_vtx_idx[i] - 1;
    cs_int_t  e = mesh->face_vtx_idx[i+1] - 1;
    cs_int_t  n_face_vertices = e - s;
    double  inv_n_face_vertices = 1/(double)n_face_vertices;

    cs_real_t  bary[3] = { 0.0, 0.0, 0.0};
    cs_real_t  fnorm[3] = { 0.0, 0.0, 0.0};

    /* Fill face_vtx_coord */

    for (j = s; j < e; j++) {
      vid = mesh->face_vtx_lst[j] - 1;
      for (k = 0; k < 3; k++)
        face_vtx_coord[shift++] = mesh->vertices[vid].coord[k];
    }

    vid = mesh->face_vtx_lst[s] - 1;
    for (k = 0; k < 3; k++)
      face_vtx_coord[shift++] = mesh->vertices[vid].coord[k];

    /* Compute the barycenter of the face */

    for (j = 0; j < n_face_vertices; j++)
      for (k = 0; k < 3; k++)
        bary[k] += face_vtx_coord[3*j+k];

    for (k = 0; k < 3; k++)
      bary[k] *= inv_n_face_vertices;

    /* Loop on the triangles of the face defined by an edge of the face
       and the barycenter */

    for (j = 0; j < n_face_vertices; j++) {

      /*    Computation of the normal of each triangle Ti:
             ->            -->   -->
             N(Ti) = 1/2 ( BPi X BPi+1 )
      */

      for (k = 0; k < 3; k++) {
        v1[k] = face_vtx_coord[3*j    + k] - bary[k];
        v2[k] = face_vtx_coord[3*(j+1)+ k] - bary[k];
      }

      _cross_product(v1, v2, tri_normal);

      /*   Computation of the normal of the polygon
           => vectorial sum of normals of each triangle

           ->      n-1   ->
           N(P) =  Sum ( N(Ti) )
                   i=0
      */

      for (k = 0; k < 3; k++)
        fnorm[k] += 0.5 * tri_normal[k];

    } /* End of loop on vertices of the face */

    inv_norm = 1/sqrt(_dot_product(fnorm, fnorm));

    for (k = 0; k < 3; k++)
      face_normal[3*i+k] = inv_norm * fnorm[k];

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    bft_printf("  Face_num: %5d (%u)- face_normal [%8.4e, %8.4e, %8.4e]\n",
               i+1, mesh->face_gnum[i],
               face_normal[3*i], face_normal[3*i+1], face_normal[3*i+2]);
    bft_printf_flush();
#endif

  } /* End of loop on faces */

  /* Free memory */

  BFT_FREE(face_vtx_coord);

  return face_normal;
}

/*----------------------------------------------------------------------------
 * Allocate and define an "edge -> face" connectivity
 *
 * parameters:
 *   mesh          <-- pointer to a cs_join_mesh_t structure
 *   edges         <-- pointer to a cs_join_edges_t structure
 *   edge_face_idx --> pointer to the edge -> face connect. index
 *   edge_face_lst --> pointer to the edge -> face connect. list
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_get_edge_face_adj(const cs_join_mesh_t   *mesh,
                               const cs_join_edges_t  *edges,
                               cs_int_t               *edge_face_idx[],
                               cs_int_t               *edge_face_lst[])
{
  cs_int_t  i, j, k, edge_id, shift;
  cs_int_t  n_edges, n_faces;

  cs_int_t  n_max_vertices = 0;
  cs_int_t  *counter = NULL;
  cs_int_t  *face_connect = NULL;
  cs_int_t  *_edge_face_idx = NULL;
  cs_int_t  *_edge_face_lst = NULL;

  if (mesh == NULL || edges == NULL)
    return;

  n_edges = edges->n_edges;
  n_faces = mesh->n_faces;

  /* Compute n_max_vertices */

  for (i = 0; i < n_faces; i++)
    n_max_vertices = CS_MAX(n_max_vertices,
                            mesh->face_vtx_idx[i+1]-mesh->face_vtx_idx[i]);

  BFT_MALLOC(face_connect, n_max_vertices + 1, cs_int_t);
  BFT_MALLOC(counter, n_edges, cs_int_t);

  /* Build an edge -> face connectivity */

  BFT_MALLOC(_edge_face_idx, n_edges+1, cs_int_t);

  for (i = 0; i < n_edges+1; i++)
    _edge_face_idx[i] = 0;

  for (i = 0; i < n_edges; i++)
    counter[i] = 0;

  /* Build index */

  for (i = 0; i < n_faces; i++) {

    cs_int_t  start_id = mesh->face_vtx_idx[i] - 1;
    cs_int_t  end_id = mesh->face_vtx_idx[i+1] - 1;
    cs_int_t  n_face_vertices = end_id - start_id;

    for (j = start_id, k = 0; j < end_id; j++, k++)
      face_connect[k] = mesh->face_vtx_lst[j];
    face_connect[n_face_vertices] = mesh->face_vtx_lst[start_id];

    assert(n_face_vertices == k);

    for (j = 0; j < n_face_vertices; j++) {

      cs_int_t  vtx_id1 = face_connect[j] - 1;

      for (k = edges->vtx_idx[vtx_id1]; k < edges->vtx_idx[vtx_id1+1]; k++)
        if (edges->adj_vtx_lst[k] == face_connect[j+1])
          break;

      assert(k != edges->vtx_idx[vtx_id1+1]);

      _edge_face_idx[CS_ABS(edges->edge_lst[k])] += 1;

    } /* End of loop on vertices of the face */

  } /* End of loop on faces */

  for (i = 0; i < n_edges; i++)
    _edge_face_idx[i+1] += _edge_face_idx[i];

  BFT_MALLOC(_edge_face_lst, _edge_face_idx[n_edges], cs_int_t);

  /* Fill "edge -> face" connectivity list */

  for (i = 0; i < n_faces; i++) {

    cs_int_t  start_id = mesh->face_vtx_idx[i] - 1;
    cs_int_t  end_id = mesh->face_vtx_idx[i+1] - 1;
    cs_int_t  n_face_vertices = end_id - start_id;

    for (j = start_id, k = 0; j < end_id; j++, k++)
      face_connect[k] = mesh->face_vtx_lst[j];
    face_connect[n_face_vertices] = mesh->face_vtx_lst[start_id];

    for (j = 0; j < n_face_vertices; j++) {

      cs_int_t  vtx_id1 = face_connect[j] - 1;

      for (k = edges->vtx_idx[vtx_id1];
           k < edges->vtx_idx[vtx_id1+1]; k++)
        if (edges->adj_vtx_lst[k] == face_connect[j+1])
          break;

      edge_id = CS_ABS(edges->edge_lst[k]) - 1;
      shift = _edge_face_idx[edge_id] + counter[edge_id];
      _edge_face_lst[shift] = i+1;
      counter[edge_id] += 1;

    } /* End of loop on vertices of the face */

  } /* End of loop on faces */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf("\n DUMP EDGE -> FACE CONNECTIVITY:\n\n");

  for (i = 0; i < n_edges; i++) {

    cs_int_t  start = _edge_face_idx[i];
    cs_int_t  end = _edge_face_idx[i+1];
    cs_int_t  v1_id = edges->def[2*i] - 1;
    cs_int_t  v2_id = edges->def[2*i+1] - 1;

    bft_printf(" edge_num: %6d (%9u) [%9u - %9u]: size: %4d, faces: ",
               i+1, edges->gnum[i],
               mesh->vertices[v1_id].gnum, mesh->vertices[v2_id].gnum,
               end-start);
    for (j = start; j < end; j++)
      bft_printf(" %u ", mesh->face_gnum[_edge_face_lst[j]-1]);
    bft_printf("\n");
    bft_printf_flush();
  }
  bft_printf("\n");
#endif

  /* Set return pointers */

  *edge_face_idx = _edge_face_idx;
  *edge_face_lst = _edge_face_lst;

  /* Free memory*/

  BFT_FREE(counter);
  BFT_FREE(face_connect);
}

/*----------------------------------------------------------------------------
 * Dump a cs_join_vertex_t structure.
 *
 * parameters:
 *   vertex <-- cs_join_vertex_t structure to dump
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_dump_vertex(const cs_join_vertex_t   vertex)
{
  assert(vertex.gnum > 0);
  assert(vertex.tolerance >= 0.0);

  bft_printf(" %10u | %11.6f |  [ %12.6g %12.6g  %12.6g]\n",
             vertex.gnum, vertex.tolerance,
             vertex.coord[0], vertex.coord[1], vertex.coord[2]);
}

/*----------------------------------------------------------------------------
 * Dump a cs_join_vertex_t structure into a file.
 *
 * parameters:
 *   file   <-- pointer to a FILE structure
 *   vertex <-- cs_join_vertex_t structure to dump
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_dump_vertex_file(FILE                   *file,
                              const cs_join_vertex_t  vertex)
{
  assert(vertex.gnum > 0);
  assert(vertex.tolerance >= 0.0);

  fprintf(file," %10u | %11.6f |  [ %12.6g  %12.6g  %12.6g]\n",
          vertex.gnum, vertex.tolerance,
          vertex.coord[0], vertex.coord[1], vertex.coord[2]);
}

/*----------------------------------------------------------------------------
 * Dump a cs_join_mesh_t structure.
 *
 * parameters:
 *   mesh <-- pointer to a cs_join_mesh_t structure to dump
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_dump(const cs_join_mesh_t  *mesh)
{
  int  i, j;

  if (mesh == NULL) {
    bft_printf(_("\n\n  -- Dump a cs_join_mesh_t structure: (%p) --\n"),
               mesh);
    return;
  }

  bft_printf(_("\n\n  -- Dump a cs_join_mesh_t structure: %s (%p) --\n"),
             mesh->name, mesh);
  bft_printf(_("\n mesh->n_faces:       %9d\n"), mesh->n_faces);
  bft_printf(_(" mesh->n_g_faces:   %11u\n\n"), mesh->n_g_faces);

  if (mesh->face_vtx_idx != NULL) {

    for (i = 0; i < mesh->n_faces; i++) {

      cs_int_t  start = mesh->face_vtx_idx[i] - 1;
      cs_int_t  end = mesh->face_vtx_idx[i+1] - 1;

      bft_printf(_("\n face_id: %5d gnum: %9u - n_vertices: %4d\n"),
                 i, mesh->face_gnum[i], end-start);

      for (j = start; j < end; j++) {

        cs_int_t  vtx_id = mesh->face_vtx_lst[j]-1;
        cs_join_vertex_t  v_data = mesh->vertices[vtx_id];

        bft_printf(" %6d - %8u - [ %8.4f %8.4f %8.4f]\n",
                   vtx_id+1, v_data.gnum,
                   v_data.coord[0], v_data.coord[1], v_data.coord[2]);

      }
      bft_printf("\n");

      /* Check if there is no incoherency in the mesh definition */

      for (j = start; j < end - 1; j++) {

        cs_int_t  vtx_id1 = mesh->face_vtx_lst[j]-1;
        cs_int_t  vtx_id2 = mesh->face_vtx_lst[j+1]-1;

        if (vtx_id1 == vtx_id2) {
          bft_printf(_("  Incoherency found in the current mesh definition\n"
                       "  Face number: %d (global: %u)\n"
                       "  Vertices: local (%d, %d), global (%u, %u)"
                       " are defined twice\n"),
                     i+1, mesh->face_gnum[i], vtx_id1+1, vtx_id2+1,
                     (mesh->vertices[vtx_id1]).gnum,
                     (mesh->vertices[vtx_id2]).gnum);
          bft_printf_flush();
          assert(0);
        }

      }

      {
        cs_int_t  vtx_id1 = mesh->face_vtx_lst[end-1]-1;
        cs_int_t  vtx_id2 = mesh->face_vtx_lst[start]-1;

        if (vtx_id1 == vtx_id2) {
          bft_printf(_("  Incoherency found in the current mesh definition\n"
                       "  Face number: %d (global: %u)\n"
                       "  Vertices: local (%d, %d), global (%u, %u)"
                       " are defined twice\n"),
                     i+1, mesh->face_gnum[i], vtx_id1+1, vtx_id2+1,
                     (mesh->vertices[vtx_id1]).gnum,
                     (mesh->vertices[vtx_id2]).gnum);
          bft_printf_flush();
          assert(0);
        }
      }

    } /* End of loop on faces */

  } /* End if face_vtx_idx != NULL */

  bft_printf(_("\n Dump vertex data\n"
               "   mesh->vertices     :  %p\n"
               "   mesh->n_vertices   : %11d\n"
               "   mesh->n_g_vertices : %11u\n\n"),
             mesh->vertices, mesh->n_vertices, mesh->n_g_vertices);

  if (mesh->n_vertices > 0) {

    bft_printf(_(" Local Num | Global Num |  Tolerance  |"
                 "        Coordinates\n\n"));

    for (i = 0; i < mesh->n_vertices; i++) {

      bft_printf(" %9d |", i+1);
      cs_join_mesh_dump_vertex(mesh->vertices[i]);

    }

  }
  bft_printf("\n");
  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Dump a cs_join_mesh_t structure into a file.
 *
 * parameters:
 *   file <-- pointer to a FILE structure
 *   mesh <-- pointer to cs_join_mesh_t structure to dump
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_dump_file(FILE                  *file,
                       const cs_join_mesh_t  *mesh)
{
  int  i, j;

  if (mesh == NULL) {
    fprintf(file,
            _("\n\n  -- Dump a cs_join_mesh_t structure: (%p) --\n"),
            (const void *)mesh);
    return;
  }

  fprintf(file,_("\n\n  -- Dump a cs_join_mesh_t structure: %s (%p) --\n"),
          mesh->name, (const void *)mesh);
  fprintf(file,_("\n mesh->n_faces:     %11d\n"), mesh->n_faces);
  fprintf(file,_(" mesh->n_g_faces:   %11u\n\n"), mesh->n_g_faces);

  if (mesh->face_vtx_idx != NULL) {

    for (i = 0; i < mesh->n_faces; i++) {

      cs_int_t  start = mesh->face_vtx_idx[i] - 1;
      cs_int_t  end = mesh->face_vtx_idx[i+1] - 1;

      fprintf(file,_("\n face_id: %9d gnum: %10u - n_vertices : %4d\n"),
              i, mesh->face_gnum[i],end-start);

      for (j = start; j < end; j++) {

        cs_int_t  vtx_id = mesh->face_vtx_lst[j]-1;
        cs_join_vertex_t  v_data = mesh->vertices[vtx_id];

        fprintf(file," %9d - %10u - [ %12.6g %12.6g %12.6g]\n",
                vtx_id+1, v_data.gnum,
                v_data.coord[0], v_data.coord[1], v_data.coord[2]);

      }
      fprintf(file,"\n");

      /* Check if there is no incoherency in the mesh definition */

      for (j = start; j < end - 1; j++) {

        cs_int_t  vtx_id1 = mesh->face_vtx_lst[j]-1;
        cs_int_t  vtx_id2 = mesh->face_vtx_lst[j+1]-1;

        if (vtx_id1 == vtx_id2) {
          fprintf(file,_("  Incoherency found in the current mesh definition\n"
                         "  Face number: %d (global: %u)\n"
                         "  Vertices: local (%d, %d), global (%u, %u)"
                         " are defined twice\n"),
                  i+1, mesh->face_gnum[i], vtx_id1+1, vtx_id2+1,
                  (mesh->vertices[vtx_id1]).gnum,
                  (mesh->vertices[vtx_id2]).gnum);
          fflush(file);
          assert(0);
        }

      }

      {
        cs_int_t  vtx_id1 = mesh->face_vtx_lst[end-1]-1;
        cs_int_t  vtx_id2 = mesh->face_vtx_lst[start]-1;

        if (vtx_id1 == vtx_id2) {
          fprintf(file,
                  _("  Incoherency found in the current mesh definition\n"
                    "  Face number: %d (global: %u)\n"
                    "  Vertices: local (%d, %d), global (%u, %u)"
                    " are defined twice\n"),
                  i+1, mesh->face_gnum[i], vtx_id1+1, vtx_id2+1,
                  (mesh->vertices[vtx_id1]).gnum,
                  (mesh->vertices[vtx_id2]).gnum);
          fflush(file);
          assert(0);
        }
      }

    } /* End of loop on faces */

  } /* End if face_vtx_idx != NULL */

  fprintf(file,_("\n Dump vertex data\n"
                 "   mesh->vertices     :  %p\n"
                 "   mesh->n_vertices   : %11d\n"
                 "   mesh->n_g_vertices : %11u\n\n"),
          (const void *)mesh->vertices, mesh->n_vertices, mesh->n_g_vertices);

  if (mesh->n_vertices > 0) {

    fprintf(file,_(" Local Num | Global Num |  Tolerance  |"
                   "        Coordinates\n\n"));

    for (i = 0; i < mesh->n_vertices; i++) {
      fprintf(file," %9d |", i+1);
      cs_join_mesh_dump_vertex_file(file, mesh->vertices[i]);
    }

  }
  fprintf(file,"\n");
  fflush(file);
}

/*----------------------------------------------------------------------------
 * Dump a list of cs_join_edge_t structures.
 *
 * parameters:
 *   edges <-- cs_join_edges_t structure to dump
 *   mesh  <-- associated cs_join_mesh_t structure
 *---------------------------------------------------------------------------*/

void
cs_join_mesh_dump_edges(const cs_join_edges_t  *edges,
                        const cs_join_mesh_t   *mesh)
{
  cs_int_t  i, j;

  if (edges == NULL)
    return;

  bft_printf(_("\n  Edge connectivity used in the joining operation:\n"));
  bft_printf(_("  Number of edges    :  %8d\n"), edges->n_edges);
  bft_printf(_("  Number of vertices :  %8d\n"), edges->n_vertices);

  for (i = 0; i < edges->n_edges; i++) { /* Dump edge connectivity */

    cs_int_t  v1_id = edges->def[2*i] - 1;
    cs_int_t  v2_id = edges->def[2*i+1] - 1;
    fvm_gnum_t  v1_gnum = (mesh->vertices[v1_id]).gnum;
    fvm_gnum_t  v2_gnum = (mesh->vertices[v2_id]).gnum;

    bft_printf(_("  Edge %6d  (%7u) <Vertex> [%7u %7u]\n"),
               i+1, edges->gnum[i], v1_gnum, v2_gnum);

    /* Check coherency */

    if (v1_id == v2_id) {
      bft_printf("  Incoherency found in the current edge definition\n"
                 "  Edge number: %d\n"
                 "  Vertices: local (%d, %d), global (%u, %u)"
                 " are defined twice\n",
                 i+1, v1_id+1, v2_id+1, v1_gnum, v2_gnum);
      bft_printf_flush();
      assert(0);
    }

    if (v1_gnum == v2_gnum) {
      bft_printf("  Incoherency found in the current edge definition\n"
                 "  Edge number: %d\n"
                 "  Vertices: local (%d, %d), global (%u, %u)"
                 " are defined twice\n",
                 i+1, v1_id+1, v2_id+1, v1_gnum, v2_gnum);
      bft_printf_flush();
      assert(0);
    }

  } /* End of loop on edges */

  bft_printf(_("\n  Vertex -> Vertex connectivity :\n\n"));

  for (i = 0; i < mesh->n_vertices; i++) {

    cs_int_t  start = edges->vtx_idx[i];
    cs_int_t  end = edges->vtx_idx[i+1];

    bft_printf(_("  Vertex %6d (%7u) - %3d - "),
               i+1, (mesh->vertices[i]).gnum, end - start);

    for (j = start; j < end; j++)
      bft_printf(" [ v: %7u, e: %7u] ",
                 (mesh->vertices[edges->adj_vtx_lst[j]-1]).gnum,
                 edges->gnum[edges->edge_lst[j]]);
    bft_printf("\n");

  }

  bft_printf_flush();
}

/*---------------------------------------------------------------------------*/

END_C_DECLS
