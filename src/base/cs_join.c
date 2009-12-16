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
 * Management of conforming and non-conforming joining
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

#include <fvm_parall.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *---------------------------------------------------------------------------*/

#include "cs_mesh_quantities.h"
#include "cs_post.h"
#include "cs_join_intersect.h"
#include "cs_join_merge.h"
#include "cs_join_mesh.h"
#include "cs_join_post.h"
#include "cs_join_set.h"
#include "cs_join_split.h"
#include "cs_join_update.h"
#include "cs_join_util.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *---------------------------------------------------------------------------*/

#include "cs_join.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Structure Definitions
 *===========================================================================*/

/*============================================================================
 * Static global variables
 *===========================================================================*/

int  cs_glob_n_joinings = 0;

cs_join_t  **cs_glob_join_array = NULL;

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
 * Compute the distance between two vertices.
 *
 * parameters:
 *   id         <-- local id of the vertex to deal with
 *   quantities <-- array keeping edge vector and its length
 *
 * returns:
 *   sinus in radian between edges sharing vertex id
 *---------------------------------------------------------------------------*/

inline static double
_compute_sinus(int           id,
               const double  quantities[])
{
  int  i;

  double  sinus;
  double  norm_a, norm_b, a[3], b[3], c[3];

  for (i = 0; i < 3; i++) {
    a[i] = -quantities[4*id+i];
    b[i] = quantities[4*(id+1)+i];
  }

  norm_a = quantities[4*id+3];
  norm_b = quantities[4*(id+1)+3];

  _cross_product(a, b, c);

  sinus = _norm(c) / (norm_a * norm_b);

  return sinus;
}

/*----------------------------------------------------------------------------
 * Compute the distance between two vertices.
 *
 * parameters:
 *   a <-- coordinates of the first vertex.
 *   b <-- coordinates of the second vertex.
 *
 * returns:
 *   distance between a and b.
 *---------------------------------------------------------------------------*/

inline static double
_compute_length(const double  a[3],
                const double  b[3])
{
  double  length;

  length = sqrt(  (b[0] - a[0])*(b[0] - a[0])
                + (b[1] - a[1])*(b[1] - a[1])
                + (b[2] - a[2])*(b[2] - a[2]));

  return length;
}

/*----------------------------------------------------------------------------
 * Compute tolerance (mode 2)
 * tolerance = min[ edge length * sinus(v1v2) * fraction]
 *
 * parameters:
 *   vertex_coords    <--  coordinates of vertices.
 *   vertex_tolerance <->  local tolerance affected to each vertex and
 *                         to be updated
 *   n_faces          <--  number of selected faces
 *   face_lst         <--  list of faces selected to compute the tolerance
 *   f2v_idx          <--  "face -> vertex" connect. index
 *   f2v_lst          <--  "face -> vertex" connect. list
 *   fraction         <--  parameter used to compute the tolerance
 *---------------------------------------------------------------------------*/

static void
_compute_tolerance2(const cs_real_t   vtx_coords[],
                    double            vtx_tolerance[],
                    const cs_int_t    n_faces,
                    const cs_int_t    face_lst[],
                    const cs_int_t    f2v_idx[],
                    const cs_int_t    f2v_lst[],
                    double            fraction)
{
  int  i, j, k, coord;
  double  tolerance, sinus;
  double  a[3], b[3];

  int   n_max_face_vertices = 0;
  int  *face_connect = NULL;
  double  *edge_quantities = NULL;

  for (i = 0; i < n_faces; i++) {
    int  fid = face_lst[i] - 1;
    n_max_face_vertices = CS_MAX(n_max_face_vertices,
                                 f2v_idx[fid+1] - f2v_idx[fid]);
  }

  BFT_MALLOC(face_connect, n_max_face_vertices + 1, int);
  BFT_MALLOC(edge_quantities, 4 * (n_max_face_vertices + 1), double);

  for (i = 0; i < n_faces; i++) {

    int  face_id = face_lst[i] - 1;
    int  start = f2v_idx[face_id] - 1;
    int  end = f2v_idx[face_id + 1] - 1;
    int  n_face_vertices = end - start;

    /* Keep face connect */

    for (k = 0, j = start; j < end; j++, k++)
      face_connect[k] = f2v_lst[j] - 1;
    face_connect[k] = f2v_lst[start] - 1;

    /* Keep edge lengths and edge vectors:
        - edge_quantities[4*k+0..2] = edge vector
        - edge_quantities[4*k+3] = edge length */

    for (k = 0; k < n_face_vertices; k++) {

      for (coord = 0; coord < 3; coord++) {
        a[coord] = vtx_coords[3*face_connect[k] + coord];
        b[coord] = vtx_coords[3*face_connect[k+1] + coord];
        edge_quantities[4*(k+1)+coord] = b[coord] - a[coord];
      }
      edge_quantities[4*(k+1)+3] = _compute_length(a, b);

    }

    for (coord = 0; coord < 4; coord++)
      edge_quantities[coord] = edge_quantities[4*k+coord];

    /* Loop on the vertices of the face to update tolerance on
       each vertex */

    for (k = 0; k < n_face_vertices; k++) {

      int  vid = face_connect[k];

      tolerance = fraction * CS_MIN(edge_quantities[4*k+3],
                                    edge_quantities[4*(k+1)+3]);
      sinus = _compute_sinus(k, edge_quantities);

      vtx_tolerance[vid] = FVM_MIN(vtx_tolerance[vid], sinus*tolerance);

    }

  } /* End of loop on faces */

  BFT_FREE(face_connect);
  BFT_FREE(edge_quantities);

}

/*----------------------------------------------------------------------------
 * Compute tolerance (mode 1)
 * tolerance = shortest edge length * fraction
 *
 * parameters:
 *   vertex_coords    <--  coordinates of vertices.
 *   vertex_tolerance <->  local tolerance affected to each vertex and
 *                         to be updated
 *   n_faces          <--  number of selected faces
 *   face_lst         <--  list of faces selected to compute the tolerance
 *   face_vtx_idx     <--  "face -> vertex" connect. index
 *   face_vtx_lst     <--  "face -> vertex" connect. list
 *   fraction         <--  parameter used to compute the tolerance
 *---------------------------------------------------------------------------*/

static void
_compute_tolerance1(const cs_real_t   vtx_coords[],
                    double            vtx_tolerance[],
                    const cs_int_t    n_faces,
                    const cs_int_t    face_lst[],
                    const cs_int_t    face_vtx_idx[],
                    const cs_int_t    face_vtx_lst[],
                    double            fraction)
{
  cs_int_t  i, j, k, start, end, face_id, vtx_id1, vtx_id2;
  double  length, tolerance;
  double  a[3], b[3];

  for (i = 0; i < n_faces; i++) {

    face_id = face_lst[i] - 1;
    start = face_vtx_idx[face_id] - 1;
    end = face_vtx_idx[face_id + 1] - 1;

    /* Loop on the vertices of the face */

    for (j = start; j < end - 1; j++) {

      vtx_id1 = face_vtx_lst[j] - 1;
      vtx_id2 = face_vtx_lst[j+1] - 1;

      for (k = 0; k < 3; k++) {
        a[k] = vtx_coords[3*vtx_id1 + k];
        b[k] = vtx_coords[3*vtx_id2 + k];
      }

      length = _compute_length(a, b);
      tolerance = length * fraction;
      vtx_tolerance[vtx_id1] = FVM_MIN(vtx_tolerance[vtx_id1], tolerance);
      vtx_tolerance[vtx_id2] = FVM_MIN(vtx_tolerance[vtx_id2], tolerance);

    }

    /* Case end - start */

    vtx_id1 = face_vtx_lst[end-1] - 1;
    vtx_id2 = face_vtx_lst[start] - 1;

    for (k = 0; k < 3; k++) {
      a[k] = vtx_coords[3*vtx_id1 + k];
      b[k] = vtx_coords[3*vtx_id2 + k];
    }

    length = _compute_length(a, b);
    tolerance = length * fraction;
    vtx_tolerance[vtx_id1] = FVM_MIN(vtx_tolerance[vtx_id1], tolerance);
    vtx_tolerance[vtx_id2] = FVM_MIN(vtx_tolerance[vtx_id2], tolerance);

  } /* End of loop on faces */

}

/*----------------------------------------------------------------------------
 * Define for each vertex a tolerance which is the radius of the
 * sphere in which the vertex can be fused with another vertex.
 * This tolerance is computed from the given list of faces (interior or border)
 *
 * parameters:
 *   param            <--  set of user-defined parameters for the joining
 *   vertex_coords    <--  coordinates of vertices.
 *   vertex_tolerance <->  local tolerance affected to each vertex and
 *                         to be updated
 *   n_faces          <--  number of selected faces
 *   face_lst         <--  list of faces selected to compute the tolerance
 *   face_vtx_idx     <--  "face -> vertex" connect. index
 *   face_vtx_lst     <--  "face -> vertex" connect. list
 *---------------------------------------------------------------------------*/

static void
_get_local_tolerance(cs_join_param_t  param,
                     const cs_real_t  vtx_coords[],
                     double           vtx_tolerance[],
                     const cs_int_t   n_faces,
                     const cs_int_t   face_lst[],
                     const cs_int_t   face_vtx_idx[],
                     const cs_int_t   face_vtx_lst[])
{

  if (param.tcm % 10 == 1) {

    /* tol = min(edge length * fraction) */

    _compute_tolerance1(vtx_coords,
                        vtx_tolerance,
                        n_faces,
                        face_lst,
                        face_vtx_idx,
                        face_vtx_lst,
                        param.fraction);

  }
  else if (param.tcm % 10 == 2) {

    /* tol = min(edge length * sin(v1v2) * fraction) */

    _compute_tolerance2(vtx_coords,
                        vtx_tolerance,
                        n_faces,
                        face_lst,
                        face_vtx_idx,
                        face_vtx_lst,
                        param.fraction);

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              _("  Tolerance computation mode (%d) is not defined\n"));

}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Exchange local vertex tolerances to get a global vertex tolerance.
 *
 * parameters:
 *   n_vertices        <-- number of local selected vertices
 *   select_vtx_io_num <-- fvm_io_num_t structure for the selected vertices
 *   vertex_data       <-> data associated to each selected vertex
 *---------------------------------------------------------------------------*/

static void
_get_global_tolerance(cs_int_t             n_vertices,
                      const fvm_io_num_t  *select_vtx_io_num,
                      cs_join_vertex_t     vtx_data[])
{
  cs_int_t  i, rank, vtx_id, block_size, shift;
  fvm_gnum_t  first_vtx_gnum;

  double  *g_vtx_tolerance = NULL, *send_list = NULL, *recv_list = NULL;
  cs_int_t  *send_count = NULL, *recv_count = NULL;
  cs_int_t  *send_shift = NULL, *recv_shift = NULL;
  fvm_gnum_t  *send_glist = NULL, *recv_glist = NULL;
  fvm_gnum_t  n_g_vertices = fvm_io_num_get_global_count(select_vtx_io_num);
  const fvm_gnum_t  *io_gnum = fvm_io_num_get_global_num(select_vtx_io_num);

  MPI_Comm  mpi_comm = cs_glob_mpi_comm;
  const int  local_rank = CS_MAX(cs_glob_rank_id, 0);
  const int  n_ranks = cs_glob_n_ranks;

  /* Define a fvm_io_num_t structure on vertices */

  block_size = n_g_vertices / n_ranks;
  if (n_g_vertices % n_ranks > 0)
    block_size += 1;

  /* Count the number of vertices to send to each rank */
  /* ------------------------------------------------- */

  BFT_MALLOC(send_count, n_ranks, int);
  BFT_MALLOC(recv_count, n_ranks, int);
  BFT_MALLOC(send_shift, n_ranks + 1, int);
  BFT_MALLOC(recv_shift, n_ranks + 1, int);

  send_shift[0] = 0;
  recv_shift[0] = 0;

  for (rank = 0; rank < n_ranks; rank++)
    send_count[rank] = 0;

  for (i = 0; i < n_vertices; i++) {
    rank = (io_gnum[i] - 1)/block_size;
    send_count[rank] += 1;
  }

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, mpi_comm);

  for (rank = 0; rank < n_ranks; rank++) {
    send_shift[rank + 1] = send_shift[rank] + send_count[rank];
    recv_shift[rank + 1] = recv_shift[rank] + recv_count[rank];
  }

  assert(send_shift[n_ranks] == n_vertices);

  /* Send the global numbering for each vertex */
  /* ----------------------------------------- */

  BFT_MALLOC(send_glist, n_vertices, fvm_gnum_t);
  BFT_MALLOC(recv_glist, recv_shift[n_ranks], fvm_gnum_t);

  for (rank = 0; rank < n_ranks; rank++)
    send_count[rank] = 0;

  for (i = 0; i < n_vertices; i++) {
    rank = (io_gnum[i] - 1)/block_size;
    shift = send_shift[rank] + send_count[rank];
    send_count[rank] += 1;
    send_glist[shift] = io_gnum[i];
  }

  MPI_Alltoallv(send_glist, send_count, send_shift, FVM_MPI_GNUM,
                recv_glist, recv_count, recv_shift, FVM_MPI_GNUM, mpi_comm);

  /* Send the vertex tolerance for each vertex */
  /* ----------------------------------------- */

  BFT_MALLOC(send_list, n_vertices, double);
  BFT_MALLOC(recv_list, recv_shift[n_ranks], double);

  for (rank = 0; rank < n_ranks; rank++)
    send_count[rank] = 0;

  for (i = 0; i < n_vertices; i++) {
    rank = (io_gnum[i] - 1)/block_size;
    shift = send_shift[rank] + send_count[rank];
    send_count[rank] += 1;
    send_list[shift] = vtx_data[i].tolerance;
  }

  MPI_Alltoallv(send_list, send_count, send_shift, MPI_DOUBLE,
                recv_list, recv_count, recv_shift, MPI_DOUBLE, mpi_comm);

  /* Define the global tolerance array */

  BFT_MALLOC(g_vtx_tolerance, block_size, double);

  for (i = 0; i < block_size; i++)
    g_vtx_tolerance[i] = DBL_MAX;

  first_vtx_gnum = block_size * local_rank + 1;

  for (i = 0; i < recv_shift[n_ranks]; i++) {
    vtx_id = recv_glist[i] - first_vtx_gnum;
    g_vtx_tolerance[vtx_id] = FVM_MIN(g_vtx_tolerance[vtx_id], recv_list[i]);
  }

  /* Replace local vertex tolerance by the new computed global tolerance */

  for (i = 0; i < recv_shift[n_ranks]; i++) {
    vtx_id = recv_glist[i] - first_vtx_gnum;
    recv_list[i] = g_vtx_tolerance[vtx_id];
  }

  MPI_Alltoallv(recv_list, recv_count, recv_shift, MPI_DOUBLE,
                send_list, send_count, send_shift, MPI_DOUBLE, mpi_comm);

  for (rank = 0; rank < n_ranks; rank++)
    send_count[rank] = 0;

  for (i = 0; i < n_vertices; i++) {
    rank = (io_gnum[i] - 1)/block_size;
    shift = send_shift[rank] + send_count[rank];
    send_count[rank] += 1;
    vtx_data[i].tolerance = send_list[shift];
  }

  /* Free memory */

  BFT_FREE(recv_glist);
  BFT_FREE(send_glist);
  BFT_FREE(send_list);
  BFT_FREE(recv_list);
  BFT_FREE(recv_count);
  BFT_FREE(send_count);
  BFT_FREE(recv_shift);
  BFT_FREE(send_shift);
  BFT_FREE(g_vtx_tolerance);
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Define a cs_join_mesh_t structure from a selection of faces and its
 * related vertices.
 *
 * parameters:
 *   name       <-- mesh name of the resulting cs_join_mesh_t structure
 *   param      <-- set of user-defined parameters for the joining
 *   selection  <-> selected entities
 *   b_f2v_idx  <-- border "face -> vertex" connectivity index
 *   b_f2v_lst  <-- border "face -> vertex" connectivity
 *   i_f2v_idx  <-- interior "face -> vertex" connectivity index
 *   i_f2v_lst  <-- interior "face -> vertex" connectivity
 *   n_vertices <-- number of vertices in the parent mesh
 *   vtx_coord  <-- coordinates of vertices in parent mesh
 *   vtx_gnum   <-- global numbering of vertices
 *---------------------------------------------------------------------------*/

static cs_join_mesh_t *
_extract_mesh(const char              *name,
              const cs_join_param_t    param,
              cs_join_select_t        *selection,
              const cs_int_t           b_f2v_idx[],
              const cs_int_t           b_f2v_lst[],
              const cs_int_t           i_f2v_idx[],
              const cs_int_t           i_f2v_lst[],
              const cs_int_t           n_vertices,
              const cs_real_t          vtx_coord[],
              const fvm_gnum_t         vtx_gnum[])
{
  int  i;
  double  clock_start, clock_end, cpu_start, cpu_end;

  double  *vtx_tolerance = NULL;
  cs_join_vertex_t  *vtx_data = NULL;
  cs_join_mesh_t  *join_mesh = NULL;

  const int  n_ranks = cs_glob_n_ranks;

  clock_start = bft_timer_wtime();
  cpu_start = bft_timer_cpu_time();

  /*
     Define a tolerance around each vertex in the selection list.
     Tolerance is the radius of the sphere in which the vertex can be merged
     with another vertex. Radius is the min(fraction * edge_length) on all
     edges connected to a vertex.
     Store all data about a vertex in a cs_join_vertex_t structure.
  */

  if (param.fraction >= 1.0 || param.fraction < 0.0)
    bft_error(__FILE__, __LINE__, 0,
              _("Incompatible value for the \"fraction\" parameter.\n"
                "Value must be lower than 1.0 or greater than 0.0\n"
                "The current value is : %f\n"), param.fraction);

  if (selection->n_vertices > 0) {

    BFT_MALLOC(vtx_tolerance, n_vertices, double);

    /* Compute the tolerance for each vertex of the mesh */

    if (param.fraction > 0.0) {

      for (i = 0; i < n_vertices; i++)
        vtx_tolerance[i] = DBL_MAX;

      /* Define local tolerance */

      _get_local_tolerance(param,
                           vtx_coord,
                           vtx_tolerance,
                           selection->n_faces,
                           selection->faces,
                           b_f2v_idx,
                           b_f2v_lst);

      if (param.tcm / 10 == 0) {

        /* Update local tolerance with adjacent border faces */

        _get_local_tolerance(param,
                             vtx_coord,
                             vtx_tolerance,
                             selection->n_b_adj_faces,
                             selection->b_adj_faces,
                             b_f2v_idx,
                             b_f2v_lst);

        /* Update local tolerance with adjacent interior faces */

        _get_local_tolerance(param,
                             vtx_coord,
                             vtx_tolerance,
                             selection->n_i_adj_faces,
                             selection->i_adj_faces,
                             i_f2v_idx,
                             i_f2v_lst);

      } /* Include adjacent faces in the computation of the vertex tolerance */

    } /* End if tolerance > 0.0 */

    else
      for (i = 0; i < n_vertices; i++)
        vtx_tolerance[i] = 0.0;


    /* Initialize vtx_data array */

    BFT_MALLOC(vtx_data, selection->n_vertices, cs_join_vertex_t);

    for (i = 0; i < selection->n_vertices; i++) {

      cs_int_t  vtx_id = selection->vertices[i]-1;

      if (n_ranks > 1)
        vtx_data[i].gnum = vtx_gnum[vtx_id];
      else
        vtx_data[i].gnum = vtx_id + 1;

      vtx_data[i].coord[0] = vtx_coord[3*vtx_id];
      vtx_data[i].coord[1] = vtx_coord[3*vtx_id+1];
      vtx_data[i].coord[2] = vtx_coord[3*vtx_id+2];

      vtx_data[i].tolerance = vtx_tolerance[vtx_id];

    }

    BFT_FREE(vtx_tolerance);

  } /* End if selection->n_vertices > 0 */

#if 1 && defined(DEBUG) && !defined(NDEBUG)   /* Sanity check */
  for (i = 0; i < selection->n_vertices; i++)
    if (vtx_data[i].tolerance > (DBL_MAX - 1.))
      bft_error(__FILE__, __LINE__, 0,
                _("Incompatible value for the \"vertex tolerance\" parameter\n"
                  "Value must be lower than DBL_MAX and current value is : %f"
                  " (global numbering : %u)\n"),
                vtx_data[i].tolerance, vtx_data[i].gnum);
#endif

  /* Parallel treatment : synchro over the ranks */

#if defined(HAVE_MPI)

  if (n_ranks > 1) {

    /* Global number of selected vertices and associated
       fvm_io_num_t structure */

    fvm_io_num_t  *select_vtx_io_num = fvm_io_num_create(selection->vertices,
                                                         vtx_gnum,
                                                         selection->n_vertices,
                                                         0);

    selection->n_g_vertices = fvm_io_num_get_global_count(select_vtx_io_num);

    _get_global_tolerance(selection->n_vertices,
                          select_vtx_io_num,
                          vtx_data);

    if (param.verbosity > 1)
      bft_printf(_("  Global number of selected vertices: %11lu\n\n"),
                 (unsigned long)(selection->n_g_vertices));

    fvm_io_num_destroy(select_vtx_io_num);

  }

#endif /* defined(HAVE_MPI) */

  /* Define the join mesh structure from the selected faces and the
     related vtx_data on selected vertices */

  join_mesh = cs_join_mesh_create_from_extract(name,
                                               selection->n_faces,
                                               selection->n_g_faces,
                                               selection->faces,
                                               selection->compact_face_gnum,
                                               b_f2v_idx,
                                               b_f2v_lst,
                                               selection->n_vertices,
                                               selection->n_g_vertices,
                                               selection->vertices,
                                               vtx_data);

  if (param.verbosity > 0)
    cs_join_mesh_minmax_tol(param, join_mesh);

  /* Free memory */

  BFT_FREE(vtx_data);

  clock_end = bft_timer_wtime();
  cpu_end = bft_timer_cpu_time();

  if (param.verbosity > 2)
    bft_printf(_("\n    Definition of local joining mesh:\n"
                 "        wall clock time:            %10.3g\n"
                 "        CPU time:                   %10.3g\n"),
               clock_end - clock_start, cpu_end - cpu_start);

  if (param.verbosity > 2)
    cs_join_post_dump_mesh("LocalMesh", join_mesh, param);

  return  join_mesh;
}

/*----------------------------------------------------------------------------
 * Define a cs_join_mesh_t structure on only faces which will be
 * potentially modified by the joining operation.
 *
 * In serial mode, this is a subset of the local join mesh.
 * In parallel mode, this is a distributed subset of the global join mesh.
 *
 * Subset is a restriction on faces which could intersect each other.
 *
 * Distribution is made so that ther is a well-balanced number of faces
 * on each rank and so that faces in the mesh are spatially coherent
 * to insure no problem for differents joining operations.
 *
 * Get the associated edges and the list of possible intersections
 * between these edges.
 *
 * parameters:
 *   param                <-- set of user-defined parameter
 *   rank_face_gnum_index <-- index on face global numering to determine
 *                            the related rank
 *   local_mesh           <-- pointer to a cs_join_mesh_t structure
 *   p_work_mesh          --> pointer to the work cs_join_mesh_t structure
 *   p_work_edges         --> pointer to the cs_join_edges_t structure
 *   p_work_face_normal   --> pointer to the normal of faces defined in
 *                            the work mesh struture.
 *   p_edge_edge_vis      --> pointer to a cs_join_gset_t structure storing
 *                            the visibility between edges
 *---------------------------------------------------------------------------*/

static void
_get_work_struct(cs_join_param_t         param,
                 const fvm_gnum_t        rank_face_gnum_index[],
                 const cs_join_mesh_t   *local_mesh,
                 cs_join_mesh_t        **p_work_mesh,
                 cs_join_edges_t       **p_work_edges,
                 cs_real_t              *p_work_face_normal[],
                 cs_join_gset_t        **p_edge_edge_vis)
{
  double  clock_start, clock_end, cpu_start, cpu_end;

  cs_int_t  n_inter_faces = 0;
  char  *mesh_name = NULL;
  cs_real_t  *face_normal = NULL;
  fvm_gnum_t  *intersect_face_gnum = NULL;
  cs_join_gset_t  *face_face_vis = NULL, *edge_edge_vis = NULL;
  cs_join_mesh_t  *work_mesh = NULL;
  cs_join_edges_t  *work_edges = NULL;

  const int  n_ranks = cs_glob_n_ranks;
  const int  local_rank = CS_MAX(cs_glob_rank_id, 0);

  /*
    Build a bounding box for each selected face.
    Find intersections between bounding boxes for the whole selected mesh
    and retrieve a list (cs_join_gset_t structure) which has been
    distributed over the ranks containing this information.
  */

  face_face_vis = cs_join_intersect_faces(param, local_mesh);

  /* Define an ordered list of all implied faces without redundancy */

  /* TODO: check if this is necessary after cleanup done in face_face_vis */

  clock_start = bft_timer_wtime();
  cpu_start = bft_timer_cpu_time();

  cs_join_gset_single_order(face_face_vis,
                            &n_inter_faces,
                            &intersect_face_gnum);

  clock_end = bft_timer_wtime();
  cpu_end = bft_timer_cpu_time();

  if (param.verbosity > 1)
    bft_printf(_("\n  Sorting possible intersections between faces:\n"
                 "      wall clock time:            %10.3g\n"
                 "      CPU time:                   %10.3g\n"),
               clock_end - clock_start, cpu_end - cpu_start);

  /* Define a distributed cs_join_mesh_t structure to store the connectivity
     of the intersecting faces associated to their bounding boxes in
     face_inter list */

  if (n_ranks > 1) {

    BFT_MALLOC(mesh_name, strlen("WorkMesh_j_n") + 2 + 5 + 1, char);
    sprintf(mesh_name,"%s%02d%s%05d",
            "WorkMesh_j", param.num, "_n", local_rank);

  }
  else {

    BFT_MALLOC(mesh_name, strlen("WorkMesh_j") + 2 + 1, char);
    sprintf(mesh_name,"%s%02d", "WorkMesh_j", param.num);

  }

  work_mesh = cs_join_mesh_create_from_glob_sel(mesh_name,
                                                n_inter_faces,
                                                intersect_face_gnum,
                                                rank_face_gnum_index,
                                                local_mesh);

  /* Define a cs_join_edges_t structure associated to a cs_join_mesh_t
     structure on which we work */

  work_edges = cs_join_mesh_define_edges(work_mesh);

  /* Transform face_inter into edge_inter */

  edge_edge_vis = cs_join_intersect_face_to_edge(work_mesh,
                                                 work_edges,
                                                 face_face_vis);

  /* Define the normal vector for each selected face before any modification */

  face_normal = cs_join_mesh_get_face_normal(work_mesh);

  /* Free memory */

  BFT_FREE(mesh_name);
  BFT_FREE(intersect_face_gnum);

  cs_join_gset_destroy(&face_face_vis);

  /* Return pointers */

  *p_work_mesh = work_mesh;
  *p_work_edges = work_edges;
  *p_edge_edge_vis = edge_edge_vis;
  *p_work_face_normal = face_normal;
}

/*----------------------------------------------------------------------------
 * Build several structures useful to join faces.
 *
 * parameters:
 *   join_param         <-- set of parameters for the joining operation
 *   join_selection     <-> list of implied entities in the joining operation
 *   mesh               <-- pointer of pointer to cs_mesh_t structure
 *   p_loc_join_mesh    --> local cs_join_mesh_t structure based on local face
 *                          selection
 *   p_work_join_mesh   --> distributed and balanced cs_join_mesh_t structure
 *                          based on the global face selection
 *   p_work_join_edges  --> edges definition related to work_join_mesh
 *   p_work_face_normal --> unitary normal for the faces of work_join_mesh
 *   p_edge_edge_vis    --> list of all potential intersections between edges
 *---------------------------------------------------------------------------*/

static void
_build_join_structures(cs_join_param_t           join_param,
                       cs_join_select_t         *join_selection,
                       const cs_mesh_t          *mesh,
                       cs_join_mesh_t          **p_loc_join_mesh,
                       cs_join_mesh_t          **p_work_join_mesh,
                       cs_join_edges_t         **p_work_join_edges,
                       cs_real_t                *p_work_face_normal[],
                       cs_join_gset_t          **p_edge_edge_vis)
{
  double  clock_start, clock_end, cpu_start, cpu_end;

  char  *mesh_name = NULL;
  cs_real_t  *work_face_normal = NULL;
  cs_join_gset_t  *edge_edge_vis = NULL;
  cs_join_mesh_t  *loc_mesh = NULL, *work_mesh = NULL;
  cs_join_edges_t  *work_edges = NULL;

  clock_start = bft_timer_wtime();
  cpu_start = bft_timer_cpu_time();

  /* Define a cs_join_mesh_structure from the selected connectivity */

  if (cs_glob_n_ranks > 1) {
    BFT_MALLOC(mesh_name, strlen("LocalMesh_n") + 5 + 1, char);
    sprintf(mesh_name,"%s%05d", "LocalMesh_n", CS_MAX(cs_glob_rank_id, 0));
  }
  else {
    BFT_MALLOC(mesh_name, strlen("LocalMesh") + 1, char);
    sprintf(mesh_name,"%s", "LocalMesh");
  }

  loc_mesh = _extract_mesh(mesh_name,
                           join_param,
                           join_selection,
                           mesh->b_face_vtx_idx,
                           mesh->b_face_vtx_lst,
                           mesh->i_face_vtx_idx,
                           mesh->i_face_vtx_lst,
                           mesh->n_vertices,
                           mesh->vtx_coord,
                           mesh->global_vtx_num);

  /* Partial free memory */

  BFT_FREE(mesh_name);

  /*
    Define a cs_join_mesh_t structure on only faces which will be
    potentially modified by the joining operation.

    In serial mode, this is a subset of the local join mesh.
    In parallel mode, this is a distributed subset of the global join mesh.

    Subset is a restriction on faces which could be intersected each other.

    Distribution is made so that there is a well-balanced number of faces
    on each rank and so that faces in the mesh are spatially coherent
    to insure no problem during the different joining operations.

    Get the associated edges and the list of potential intersections
    between these edges through an edge-edge visibility.
  */

  _get_work_struct(join_param,
                   join_selection->compact_rank_index,
                   loc_mesh,
                   &work_mesh,
                   &work_edges,
                   &work_face_normal,
                   &edge_edge_vis);

  clock_end = bft_timer_wtime();
  cpu_end = bft_timer_cpu_time();

  if (join_param.verbosity > 1)
    bft_printf(_("\n  Definition of structures for the joining algorithm:\n"
                 "      wall clock time:            %10.3g\n"
                 "      CPU time:                   %10.3g\n"),
               clock_end - clock_start, cpu_end - cpu_start);
  bft_printf_flush();

  /* Return pointers */

  *p_loc_join_mesh = loc_mesh;
  *p_work_join_mesh = work_mesh;
  *p_work_join_edges = work_edges;
  *p_edge_edge_vis = edge_edge_vis;
  *p_work_face_normal = work_face_normal;

}


/*----------------------------------------------------------------------------
 * From real intersection between edges, define new vertices and/or
 * update old vertices.
 * Keep the relation between two intersecting edges through an equivalence
 * between the vertex of each edge.
 * Store also the new description of initial edges through a
 * cs_join_inter_edges_t structure and synchronize it to get all the
 * possible equivalences between vertices.
 *
 * parameters:
 *   param                <--  set of user-defined parameter for the joining
 *   work_join_mesh       <->  pointer to a cs_join_mesh_t structure
 *   work_join_edges      <--  pointer to a cs_join_edges_t structure
 *   p_edge_edge_vis      <->  pointer to a cs_join_glist_t structure
 *                             (freed here)
 *   n_g_ifm_vertices     <--  global number of vertices on the full mesh before
 *                             joining. Use to create the new glob. vertex num.
 *   p_n_g_new_vertices   -->  global number of vertices created during the
 *                             intersection of edges
 *   p_vtx_eset           -->  structure storing equivalences between vertices
 *                             Two vertices are equivalent if they are each
 *                             other in their tolerance
 *   p_inter_edges        -->  structure storing the definition of new vertices
 *                             on initial edges
 *---------------------------------------------------------------------------*/

static void
_intersect_edges(cs_join_param_t          param,
                 cs_join_mesh_t          *work_join_mesh,
                 const cs_join_edges_t   *work_join_edges,
                 cs_join_gset_t         **p_edge_edge_vis,
                 fvm_gnum_t               n_g_ifm_vertices,
                 fvm_gnum_t              *p_n_g_new_vertices,
                 cs_join_eset_t         **p_vtx_eset,
                 cs_join_inter_edges_t  **p_inter_edges)
{
  double  clock_start, clock_end, cpu_start, cpu_end;
  cs_join_type_t  join_type;

  fvm_gnum_t  n_g_new_vertices = 0;
  cs_join_inter_edges_t  *inter_edges = NULL;
  cs_join_eset_t  *vtx_eset = NULL;
  cs_join_inter_set_t  *inter_set = NULL;

  const int  n_ranks = cs_glob_n_ranks;

  clock_start = bft_timer_wtime();
  cpu_start = bft_timer_cpu_time();

  /*
     Compute the intersections between edges.
     Store the output in two data structures:
      - a cs_join_eset_t struct. to store equiv. between vertices
        issued from the same intersection
      - a cs_join_inter_set_t struct. to store detected intersections
     Return the type of the joining operation: conform or not.
  */

  join_type = cs_join_intersect_edges(param,
                                      *p_edge_edge_vis,
                                      work_join_edges,
                                      work_join_mesh,
                                      &vtx_eset,
                                      &inter_set);

  cs_join_gset_destroy(p_edge_edge_vis);

#if 0 && defined(DEBUG) && !defined(NDEBUG) /* Dump structures after inter. */
  cs_join_inter_set_dump(inter_set, work_join_edges, work_join_mesh);
#endif

  /* Synchronize join_type */

#if defined(HAVE_MPI)
  if (n_ranks > 1) {

    int  tag = (join_type == CS_JOIN_TYPE_CONFORM ? 0 : 1);
    int  sync_tag = tag;

    MPI_Allreduce(&tag, &sync_tag, 1, MPI_INT, MPI_MAX, cs_glob_mpi_comm);

    if (sync_tag == 0)
      join_type = CS_JOIN_TYPE_CONFORM;
    else
      join_type = CS_JOIN_TYPE_NO_CONFORM;

  }
#endif

  if (join_type == CS_JOIN_TYPE_CONFORM) {

    bft_printf(_("\n  Joining operation is conforming.\n"));
    bft_printf_flush();

    inter_set = cs_join_inter_set_destroy(inter_set);

  }
  else {

    assert(join_type == CS_JOIN_TYPE_NO_CONFORM);

    bft_printf(_("\n  Joining operation is non-conforming.\n"));
    bft_printf_flush();

    /* Creation of new vertices. Update list of equivalent vertices.
       Associate to each intersection a vertex (old or created) */

    cs_join_create_new_vertices(param.verbosity,
                                work_join_edges,
                                work_join_mesh,
                                inter_set,
                                n_g_ifm_vertices,
                                &n_g_new_vertices,
                                &vtx_eset);

    inter_edges = cs_join_inter_edges_define(work_join_edges, inter_set);
    inter_set = cs_join_inter_set_destroy(inter_set);

#if 0 && defined(DEBUG) && !defined(NDEBUG) /* Dump structures after inter. */
    cs_join_inter_edges_dump(inter_edges, work_join_edges, work_join_mesh);
#endif

    /* Synchronize inter_edges structure definition */

#if defined(HAVE_MPI)
    if (n_ranks > 1 ) {

      cs_join_inter_edges_t  *sync_block = NULL;

      sync_block = cs_join_inter_edges_part_to_block(work_join_mesh,
                                                     work_join_edges,
                                                     inter_edges);

      cs_join_inter_edges_block_to_part(work_join_edges->n_g_edges,
                                        sync_block,
                                        inter_edges);

      /* Add new vertices to edge description if necessary */

      cs_join_intersect_update_struct(work_join_edges,
                                      work_join_mesh,
                                      &inter_edges);

      sync_block = cs_join_inter_edges_destroy(sync_block);

      cs_join_mesh_sync_vertices(work_join_mesh);

    }
#endif

    /* Find if there are new equivalences between vertices on a same edge */

    cs_join_add_equiv_from_edges(param,
                                 work_join_mesh,
                                 work_join_edges,
                                 inter_edges,
                                 vtx_eset);

  } /* no conform joining operation */

  /* Order and delete redundant equivalences */

  cs_join_eset_clean(&vtx_eset);

  /* Memory management: final state for vtx_eset (no more equiv. to get) */

  vtx_eset->n_max_equiv = vtx_eset->n_equiv;
  BFT_REALLOC(vtx_eset->equiv_couple, 2*vtx_eset->n_equiv, cs_int_t);

  clock_end = bft_timer_wtime();
  cpu_end = bft_timer_cpu_time();

  if (param.verbosity > 0)
    bft_printf(_("\n"
                 "  Edge intersections and vertex creation:\n"
                 "    wall clock time:            %10.3g\n"
                 "    CPU time:                   %10.3g\n"),
               clock_end - clock_start, cpu_end - cpu_start);

  /* Returns pointers */

  *p_vtx_eset = vtx_eset;
  *p_inter_edges = inter_edges;
  *p_n_g_new_vertices = n_g_new_vertices;

#if 0 && defined(DEBUG) && !defined(NDEBUG) /* Dump structures after inter. */
  cs_join_inter_edges_dump(inter_edges, work_join_edges, work_join_mesh);
#endif

}

/*----------------------------------------------------------------------------
 * Merge vertices from equivalences found between vertices.
 * Update local and work structures after the merge step.
 *
 * parameters:
 *  param                <--  set of user-defined parameter for the joining
 *  join_select          <--  list of all implied entities in the
 *                            joining operation
 *  n_iwm_vertices       <--  initial number of vertices in work struct.
 *  n_g_ifm_vertices     <--  initial global number of vertices for the full
 *                            mesh
 *  n_g_new_vertices     <--  global number of vertices created with the
 *                            intersection of edges
 *  rank_face_gnum_index <--  index on face global numering to determine the
 *                            related rank
 *  vtx_equiv_set        <--  structure storing equivalences between vertices
 *                            Two vertices are equivalent if they are each
 *                            other in their tolerance
 *  inter_edges          <--  structure storing the definition of new vertices
 *                            on initial edges.
 *  p_work               <->  pointer to a cs_join_mesh_t structure
 *  p_edges              <->  pointer to a cs_join_edges_t structure
 *  p_local_join_mesh    <->  pointer to a cs_join_mesh_t structure
 *  mesh                 <->  pointer to a cs_mesh_t struct. to update
 *---------------------------------------------------------------------------*/

static void
_merge_vertices(cs_join_param_t           param,
                cs_join_select_t         *join_select,
                cs_int_t                  n_iwm_vertices,
                fvm_gnum_t                n_g_ifm_vertices,
                fvm_gnum_t                n_g_new_vertices,
                cs_join_eset_t           *vtx_eset,
                cs_join_inter_edges_t    *inter_edges,
                cs_join_mesh_t          **p_work_join_mesh,
                cs_join_edges_t         **p_work_join_edges,
                cs_join_mesh_t          **p_local_join_mesh,
                cs_mesh_t                *mesh)
{
  int  i;
  double  clock_start, clock_end, cpu_start, cpu_end;

  fvm_gnum_t  n_g_ai_vertices = n_g_ifm_vertices + n_g_new_vertices;
  fvm_gnum_t  *rank_face_gnum_index = join_select->compact_rank_index;
  fvm_gnum_t  *iwm_vtx_gnum = NULL;
  fvm_gnum_t  *o2n_vtx_gnum = NULL;
  cs_join_mesh_t  *local_join_mesh = *p_local_join_mesh;
  cs_join_mesh_t  *work_join_mesh = *p_work_join_mesh;
  cs_join_edges_t  *work_join_edges = *p_work_join_edges;

  assert(local_join_mesh != NULL);
  assert(work_join_mesh != NULL);
  assert(work_join_edges != NULL);

  clock_start = bft_timer_wtime();
  cpu_start = bft_timer_cpu_time();

  /*
    Store the initial global vertex numbering
    Initial vertices are between [0, n_init_vertices[
    Added vertices from inter. are between [n_init_vertices, n_vertices]
  */

  BFT_MALLOC(iwm_vtx_gnum, n_iwm_vertices, fvm_gnum_t);

  for (i = 0; i < n_iwm_vertices; i++)
    iwm_vtx_gnum[i] = (work_join_mesh->vertices[i]).gnum;

  /* Merge vertices */

  cs_join_merge_vertices(param,
                         n_g_ai_vertices,  /* ai: after intersection */
                         work_join_mesh,
                         vtx_eset);

  cs_join_eset_destroy(&vtx_eset);

  /*  Keep the evolution of vertex global numbering.
      Update work and local structures after vertex merge */

  cs_join_merge_update_struct(param,
                              n_iwm_vertices,
                              n_g_ifm_vertices,
                              iwm_vtx_gnum,
                              rank_face_gnum_index,
                              &work_join_mesh,
                              &work_join_edges,
                              &inter_edges,
                              &local_join_mesh,
                              &o2n_vtx_gnum);

  /* Free memory */

  BFT_FREE(iwm_vtx_gnum);

  inter_edges = cs_join_inter_edges_destroy(inter_edges);

  /* Post if required and level of verbosity is reached */

  if (param.verbosity > 2)
    cs_join_post_dump_mesh("MergeWorkMesh", work_join_mesh, param);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    cs_join_mesh_dump_edges(work_join_edges, work_join_mesh);
#endif

  /* Update cs_mesh_t structure after the vertex merge */

  cs_join_update_mesh_after_merge(param,
                                  join_select,
                                  o2n_vtx_gnum,       /* free inside */
                                  local_join_mesh,
                                  mesh);

  clock_end = bft_timer_wtime();
  cpu_end = bft_timer_cpu_time();

  if (param.verbosity > 0)
    bft_printf(_("\n"
                 "  Merge vertices:\n"
                 "    wall clock time:            %10.3g\n"
                 "    CPU time:                   %10.3g\n"),
               clock_end - clock_start, cpu_end - cpu_start);
  bft_printf_flush();

  /* Set return pointers */

  *p_local_join_mesh = local_join_mesh;
  *p_work_join_mesh = work_join_mesh;
  *p_work_join_edges = work_join_edges;
}

/*----------------------------------------------------------------------------
 * Split faces and update cs_mesh_t structure.
 *
 * parameters:
 *  param                <--  set of user-defined parameter
 *  join_select          <--  list of all implied entities in the
 *                            joining operation
 *  work_join_edges      <--  pointer to a cs_join_edges_t structure
 *  work_face_normal     <--  normal based on the original face definition
 *  rank_face_gnum_index <--  index on face global numering to determine the
 *                            related rank
 *  p_work_join_mesh     <->  pointer to a cs_join_mesh_t structure
 *  local_join_mesh      <--  pointer to a cs_join_mesh_t structure
 *  p_mesh               <->  pointer to cs_mesh_t struct.
 *---------------------------------------------------------------------------*/

static void
_split_faces(cs_join_param_t      param,
             cs_join_select_t    *join_select,
             cs_join_edges_t     *work_join_edges,
             fvm_coord_t         *work_face_normal,
             cs_join_mesh_t     **p_work_join_mesh,
             cs_join_mesh_t      *local_join_mesh,
             cs_mesh_t          **p_mesh)
{
  double  clock_start, clock_end, cpu_start, cpu_end;

  fvm_gnum_t  *rank_face_gnum_index = join_select->compact_rank_index;
  cs_join_gset_t  *old2new_hist = NULL;
  cs_mesh_t  *mesh = *p_mesh;

  clock_start = bft_timer_wtime();
  cpu_start = bft_timer_cpu_time();

  cs_join_split_faces(param,
                      work_face_normal,
                      work_join_edges,
                      p_work_join_mesh,
                      &old2new_hist);

  /* Send back to the original rank the new face description */

  cs_join_split_update_struct(*p_work_join_mesh,
                              rank_face_gnum_index,
                              &old2new_hist,
                              &local_join_mesh);

  /* Update cs_mesh_t structure after the face splitting */

  cs_join_update_mesh_after_split(param,
                                  join_select,
                                  old2new_hist,
                                  local_join_mesh,
                                  mesh);

  clock_end = bft_timer_wtime();
  cpu_end = bft_timer_cpu_time();

  /* Partial free memory */

  if (param.verbosity > 0)
    bft_printf(_("\n"
                 "  Split old faces and reconstruct new faces\n"
                 "    wall clock time:            %10.3g\n"
                 "    CPU time:                   %10.3g\n"),
               clock_end - clock_start, cpu_end - cpu_start);
  bft_printf_flush();

  /* Post if required and level of verbosity is reached */

  if (param.verbosity > 2)
    cs_join_post_dump_mesh("SplitWorkMesh", *p_work_join_mesh, param);

  /* Free memory */

  cs_join_gset_destroy(&old2new_hist);
}

/*----------------------------------------------------------------------------
 * Delete all cs_join_t structures.
 *---------------------------------------------------------------------------*/

static void
_destroy_all_joinings(void)
{
  cs_int_t  i;

  for (i = 0; i < cs_glob_n_joinings; i++) {

    cs_join_t  *join = cs_glob_join_array[i];

    BFT_FREE(join->criteria);
    BFT_FREE(join);

    cs_glob_join_array[i] = NULL;

  } /* End of loop on cs_join_t structures */

  cs_glob_n_joinings = 0;
  BFT_FREE(cs_glob_join_array);

}

/*============================================================================
 *  Public function prototypes for Fortran API
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Get the number of joining operations already defined
 *
 * Fortran Interface:
 *
 * SUBROUTINE NUMJOI
 * *****************
 *
 * INTEGER        numjoi       : --> : number of joining op. already defined
 *----------------------------------------------------------------------------*/

void CS_PROCF(numjoi, NUMJOI)
(
 cs_int_t    *numjoi
)
{
  *numjoi = cs_glob_n_joinings;
}

/*----------------------------------------------------------------------------
 * Define new boundary faces joining.
 *
 * Fortran Interface:
 *
 * SUBROUTINE DEFJO1
 * *****************
 *
 * INTEGER        numjoi           : <-- : number related to the joining op.
 * CHARACTER*     joining_criteria : <-- : boundary face selection criteria,
 * REAL           fraction         : <-- : parameter for merging vertices
 * REAL           plane            : <-- : parameter for splitting faces
 * INTEGER        verbosity        : <-- : verbosity level
 * INTEGER        joining_c_len    : <-- : length of joining_criteria
 *---------------------------------------------------------------------------*/

void CS_PROCF(defjo1, DEFJO1)
(
 cs_int_t    *numjoi,
 const char  *joining_criteria,
 cs_real_t   *fraction,
 cs_real_t   *plane,
 cs_int_t    *verbosity,
 cs_int_t    *joining_c_len
 CS_ARGF_SUPP_CHAINE
)
{
  char *_joining_criteria = NULL;

  if (joining_criteria != NULL && *joining_c_len > 0)
    _joining_criteria = cs_base_string_f_to_c_create(joining_criteria,
                                                     *joining_c_len);
  if (_joining_criteria != NULL && strlen(_joining_criteria) == 0)
    cs_base_string_f_to_c_free(&_joining_criteria);

  cs_join_add(*numjoi,
              _joining_criteria,
              *fraction,
              *plane,
              *verbosity);

  if (_joining_criteria != NULL)
    cs_base_string_f_to_c_free(&_joining_criteria);
}

/*----------------------------------------------------------------------------
 * Set advanced parameters for the joining algorithm.
 *
 * Fortran Interface:
 *
 * SUBROUTINE SETAJP
 * *****************
 *
 * INTEGER      join_num          : <-- : join number
 * REAL         mtf               : <-- : merge tolerance coefficient
 * REAL         pmf               : <-- : pre-merge factor
 * INTEGER      tcm               : <-- : tolerance computation mode
 * INTEGER      icm               : <-- : intersection computation mode
 * INTEGER      maxbrk            : <-- : max number of tolerance reduction
 * INTEGER      max_sub_faces     : <-- : max. possible number of sub-faces
 *                                        by splitting a selected face
 * INTEGER      tml               : <-- : tree max level
 * INTEGER      tmb               : <-- : tree max boxes
 * REAL         tmr               : <-- : tree max ratio
 *---------------------------------------------------------------------------*/

void CS_PROCF(setajp, SETAJP)
(
 cs_int_t    *join_num,
 cs_real_t   *mtf,
 cs_real_t   *pmf,
 cs_int_t    *tcm,
 cs_int_t    *icm,
 cs_int_t    *maxbrk,
 cs_int_t    *max_sub_faces,
 cs_int_t    *tml,
 cs_int_t    *tmb,
 cs_real_t   *tmr
 CS_ARGF_SUPP_CHAINE
)
{
  cs_int_t  i, join_id = -1;
  cs_join_t  *join = NULL;

  /* Look for the joining structure related to "join_num" */

  for (i = 0; i < cs_glob_n_joinings; i++) {

    join = cs_glob_join_array[i];
    if (*join_num == join->param.num) {
      join_id = i;
      break;
    }

  }

  if (join_id < 0)
    bft_error(__FILE__, __LINE__, 0,
              _("  Join number %d is not defined.\n"), *join_num);

  assert(join != NULL);

  cs_join_set_advanced_param(join,
                             *mtf,
                             *pmf,
                             *tcm,
                             *icm,
                             *maxbrk,
                             *max_sub_faces,
                             *tml,
                             *tmb,
                             *tmr);

}

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
 *   verbosity     <-- level of verbosity required
 *---------------------------------------------------------------------------*/

void
cs_join_add(int     join_number,
            char   *sel_criteria,
            float   fraction,
            float   plane,
            int     verbosity)
{
  size_t  l;

  cs_join_t  *join = NULL;

  /* Check parameters value */

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

   /* Allocate and initialize a cs_join_t structure */

  BFT_REALLOC(cs_glob_join_array, cs_glob_n_joinings + 1, cs_join_t *);
  BFT_MALLOC(join, 1, cs_join_t);

  join->param = cs_join_param_define(join_number,
                                     fraction,
                                     plane,
                                     verbosity);

  /* Copy the selection criteria for future use */

  l = strlen(sel_criteria);
  BFT_MALLOC(join->criteria, l + 1, char);
  strcpy(join->criteria, sel_criteria);

  /* Update global array */

  cs_glob_join_array[cs_glob_n_joinings] = join;
  cs_glob_n_joinings++;
}

/*----------------------------------------------------------------------------
 * Set advanced parameters to user-defined values.
 *
 * parameters:
 *   join           <-> pointer a to cs_join_t struct. to update
 *   mtf            <-- merge tolerance coefficient
 *   pmf            <-- pre-merge factor
 *   tcm            <-- tolerance computation mode
 *   icm            <-- intersection computation mode
 *   maxbrk         <-- max number of equivalences to break (merge step)
 *   max_sub_faces  <-- max. possible number of sub-faces by splitting a face
 *   tml            <-- tree max level
 *   tmb            <-- tree max boxes
 *   tmr            <-- tree max ratio
 *---------------------------------------------------------------------------*/

void
cs_join_set_advanced_param(cs_join_t   *join,
                           cs_real_t    mtf,
                           cs_real_t    pmf,
                           cs_int_t     tcm,
                           cs_int_t     icm,
                           cs_int_t     maxbrk,
                           cs_int_t     max_sub_faces,
                           cs_int_t     tml,
                           cs_int_t     tmb,
                           cs_real_t    tmr)
{
  /* Deepest level reachable during tree building */

  if (tml < 1)
    bft_error(__FILE__, __LINE__, 0,
              _("Mesh joining:"
                "  Forbidden value for the tml parameter.\n"
                "  It must be between > 0 and is here: %d\n"), tml);

  join->param.tree_max_level = tml;

  /* Max. number of boxes which can be related to a leaf of the tree
     if level != tree_max_level */

  if (tmb < 1)
    bft_error(__FILE__, __LINE__, 0,
              _("Mesh joining:"
                "  Forbidden value for the tmb parameter.\n"
                "  It must be between > 0 and is here: %d\n"), tmb);

  join->param.tree_n_max_boxes = tmb;

  /* Stop tree building if:
     n_linked_boxes > tree_max_box_ratio*n_init_boxes */

  if (tmr <= 0.0 )
    bft_error(__FILE__, __LINE__, 0,
              _("Mesh joining:"
                "  Forbidden value for the tmr parameter.\n"
                "  It must be between > 0.0 and is here: %f\n"), tmr);

  join->param.tree_max_box_ratio = tmr;

  /* Coef. used to modify the tolerance associated to each vertex BEFORE the
     merge operation.
     If coef = 0.0 => no vertex merge
     If coef < 1.0 => reduce vertex merge
     If coef = 1.0 => no change
     If coef > 1.0 => increase vertex merge */

  if (mtf < 0.0)
    bft_error(__FILE__, __LINE__, 0,
              _("Mesh joining:"
                "  Forbidden value for the merge tolerance factor.\n"
                "  It must be positive or nul and not: %f\n"), mtf);

  join->param.merge_tol_coef = mtf;

   /* Maximum number of equivalence breaks */

  if (maxbrk < 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Mesh joining:"
                "  Forbidden value for the max. number of tolerance breaks.\n"
                "  It must be between >= 0 and not: %d\n"), maxbrk);

  join->param.n_max_equiv_breaks = maxbrk;

  /* Pre-merge factor. This parameter is used to define a limit
     under which two vertices are merged before the merge step.
     Tolerance limit for the pre-merge = pmf * fraction
     Default value: 0.10 */

  join->param.pre_merge_factor = pmf;

  /* Tolerance computation mode */

  if ( (tcm)%10 < 1 || (tcm)%10 > 2)
    bft_error(__FILE__, __LINE__, 0,
              _("Mesh joining:"
                "  Forbidden value for the tcm parameter.\n"
                "  It must be between 1, 2 or 11, 12 and here is: %d\n"), tcm);

  join->param.tcm = tcm;

  /* Intersection computation mode */

  if (icm != 1 && icm != 2)
    bft_error(__FILE__, __LINE__, 0,
              _("Mesh joining:"
                "  Forbidden value for icm parameter.\n"
                "  It must be 1 or 2 and here is: %d\n"), icm);

  join->param.icm = icm;

  /* Maximum number of sub-faces */

  if (max_sub_faces < 1)
    bft_error(__FILE__, __LINE__, 0,
              _("Mesh joining:"
                "  Forbidden value for the maxsf parameter.\n"
                "  It must be between > 0 and here is: %d\n"), max_sub_faces);

  join->param.max_sub_faces = max_sub_faces;

}

/*----------------------------------------------------------------------------
 * Apply all the defined joining operations.
 *---------------------------------------------------------------------------*/

void
cs_join_all(void)
{
  cs_int_t  join_id;
  double  clock_start, clock_end, cpu_start, cpu_end;
  double  full_clock_start, full_clock_end, full_cpu_start, full_cpu_end;

  cs_real_t  *b_face_cog = NULL, *b_face_normal = NULL;
  cs_mesh_t  *mesh = cs_glob_mesh;

  if (cs_glob_n_joinings < 1)
    return;

  /* Sanity checks */

  assert(sizeof(cs_int_t) == sizeof(fvm_lnum_t));
  assert(sizeof(double) == sizeof(cs_real_t));

  full_clock_start = bft_timer_wtime();
  full_cpu_start = bft_timer_cpu_time();

  cs_join_post_init();

  /* Loop on each defined joining to deal with */

  for (join_id = 0; join_id < cs_glob_n_joinings; join_id++) {

    cs_join_t  *join_info = cs_glob_join_array[join_id];
    cs_join_param_t  join_param = join_info->param;
    cs_join_select_t  *join_select = NULL;

    clock_start = bft_timer_wtime();  /* Start timer */
    cpu_start = bft_timer_cpu_time();

    /* Print information into listing file */

    bft_printf(_("\n -------------------------------------------------------\n"
                 "  Joining number %d:\n\n"), join_id + 1);
    bft_printf(_("  Selection criteria: \"%s\"\n"), join_info->criteria);

    if (join_param.verbosity > 0) {
      bft_printf(_("\n"
                   "  Parameters for the joining operation:\n"
                   "    Shortest incident edge fraction:          %8.5f\n"
                   "    Maximum angle between joined face planes: %8.5f\n\n"),
                 join_param.fraction, join_param.plane);

      if (join_param.verbosity > 1)
         bft_printf(_("  Advanced join parameters:\n"
                      "    Deepest level reachable in tree building: %8d\n"
                      "    Max boxes by leaf:                        %8d\n"
                      "    Max ratio of linked boxes / init. boxes:  %8.5f\n"
                      "    Merge step tolerance multiplier:          %8.5f\n"
                      "    Pre-merge factor:                         %8.5f\n"
                      "    Tolerance computation mode:               %8d\n"
                      "    Intersection computation mode:            %8d\n"
                      "    Max. number of equiv. breaks:             %8d\n"
                      "    Max. number of subfaces by face:          %8d\n\n"),
                    join_param.tree_max_level,
                    join_param.tree_n_max_boxes,
                    join_param.tree_max_box_ratio,
                    join_param.merge_tol_coef,
                    join_param.pre_merge_factor,
                    join_param.tcm, join_param.icm,
                    join_param.n_max_equiv_breaks,
                    join_param.max_sub_faces);

      cs_mesh_print_info(mesh, _(" Before joining"));
      bft_printf("\n");
    }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    {
      int  len;
      FILE  *dbg_file = NULL;
      char  *filename = NULL;

      len = strlen("JoinDBG_InitMesh_.dat")+1+4+2;
      BFT_MALLOC(filename, len, char);
      sprintf(filename, "Join%02dDBG_InitMesh_%04d.dat",
              join_id+1, fvm_parall_get_rank());
      dbg_file = fopen(filename, "w");

      cs_mesh_dump_file(dbg_file, mesh);

      fflush(dbg_file);
      BFT_FREE(filename);
      fclose(dbg_file);
    }
#endif

    /* Build arrays and structures required for selection;
       will be destoyed after joining and rebuilt once all
       join operations are finished */

    cs_mesh_init_group_classes(mesh);

    cs_mesh_quantities_b_faces(mesh, &b_face_cog, &b_face_normal);

    cs_glob_mesh->select_b_faces
      = fvm_selector_create(mesh->dim,
                            mesh->n_b_faces,
                            mesh->class_defs,
                            mesh->b_face_family,
                            1,
                            b_face_cog,
                            b_face_normal);

    /* Get selected faces for this joining and define the related
       cs_join_face_select_t structure.
       - Compute the global number of selected faces
       - Get the adjacent faces, ...  */

    join_select = cs_join_select_create(join_info->criteria,
                                        join_param.verbosity);

    bft_printf(_("\n  Element selection successfully done.\n"));
    bft_printf_flush();

    /* Free arrays and structures needed for selection */

    BFT_FREE(b_face_cog);
    BFT_FREE(b_face_normal);

    mesh->class_defs = fvm_group_class_set_destroy(mesh->class_defs);

    if (mesh->select_b_faces != NULL)
      mesh->select_b_faces = fvm_selector_destroy(mesh->select_b_faces);
    if (mesh->class_defs != NULL)
      mesh->class_defs = fvm_group_class_set_destroy(mesh->class_defs);

    /* Now execute the joining operation */

    if (join_select->n_g_faces > 0) {

      cs_int_t  n_iwm_vertices;      /* iwm: initial work mesh */
      fvm_gnum_t  n_g_ifm_vertices;  /* ifm: initial full mesh */
      fvm_gnum_t  n_g_new_vertices;

      cs_real_t  *work_face_normal = NULL;
      cs_join_gset_t  *edge_edge_visibility = NULL;
      cs_join_mesh_t  *work_join_mesh = NULL, *local_join_mesh = NULL;
      cs_join_edges_t  *work_join_edges = NULL;
      cs_join_eset_t  *vtx_eset = NULL;
      cs_join_inter_edges_t  *inter_edges = NULL;

      _build_join_structures(join_param,
                             join_select,
                             mesh,
                             &local_join_mesh,
                             &work_join_mesh,
                             &work_join_edges,
                             &work_face_normal,
                             &edge_edge_visibility);

      n_iwm_vertices = work_join_mesh->n_vertices;
      n_g_ifm_vertices = mesh->n_g_vertices;

      if (join_param.verbosity > 2)
        bft_printf(_("\n  Number of faces to treat locally: %10d\n"),
                   work_join_mesh->n_faces);

      /*

        Define new vertices and/or update old vertices from the real
        intersection found between edges,
        Keep the relation between two intersecting edges through an
        equivalence between the vertex of each edge.
        Store also the new description of the initial edges through a
        cs_join_inter_edges_t structure and synchronize it to get all the
        possible equivalences between vertices.
        Work mesh structure is not yet fully updated by the new vertices
        because the synchronization step has to be done.

      */

      _intersect_edges(join_param,
                       work_join_mesh,
                       work_join_edges,
                       &edge_edge_visibility, /* free during this step */
                       mesh->n_g_vertices,
                       &n_g_new_vertices,
                       &vtx_eset,
                       &inter_edges);

      /*
         Merge vertices from equivalences found between vertices.
         Update work structures after the merge step.
         Keep the evolution of the global numbering of initial vertices and
         get the sync_block cs_join_inter_edges_t structure to enable the
         local structure update after the merge step.
      */

      _merge_vertices(join_param,
                      join_select,
                      n_iwm_vertices,
                      n_g_ifm_vertices,
                      n_g_new_vertices,
                      vtx_eset,         /* free during this step */
                      inter_edges,      /* free during this step */
                      &work_join_mesh,
                      &work_join_edges,
                      &local_join_mesh,
                      mesh);

      /*
         Split faces in work_join_mesh. Apply modification to the
         local_join_mesh. Keep a history between old --> new faces.
         Update cs_mesh_t structure.
      */

      _split_faces(join_param,
                   join_select,
                   work_join_edges,
                   work_face_normal,
                   &work_join_mesh,
                   local_join_mesh,
                   &mesh);

      /* Free memory */

      cs_join_mesh_destroy(&local_join_mesh);
      cs_join_mesh_destroy(&work_join_mesh);
      cs_join_mesh_destroy_edges(&work_join_edges);

      BFT_FREE(work_face_normal);

      /* Clean mesh (delete redundant edge definition) */

      cs_join_update_mesh_clean(join_param, mesh);

    }
    else
      bft_printf(_("\nStop joining algorithm: no face selected...\n"));

    /* Free memory */

    cs_join_select_destroy(&join_select);

    clock_end = bft_timer_wtime();
    cpu_end = bft_timer_cpu_time();

    bft_printf(_("\n"
                 "  Complete joining treatment for joining %2d\n"
                 "    wall clock time:            %10.3g\n"
                 "    CPU time:                   %10.3g\n"),
               join_id+1, clock_end - clock_start, cpu_end - cpu_start);
    bft_printf_flush();

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    {
      int  len;
      FILE  *dbg_file = NULL;
      char  *filename = NULL;

      len = strlen("JoinDBG_FinalMesh_.dat")+1+4+2;
      BFT_MALLOC(filename, len, char);
      sprintf(filename, "Join%02dDBG_FinalMesh_%04d.dat",
              join_id+1, fvm_parall_get_rank());
      dbg_file = fopen(filename, "w");

      cs_mesh_dump_file(dbg_file, mesh);

      fflush(dbg_file);
      BFT_FREE(filename);
      fclose(dbg_file);
    }
#endif

    if (join_param.verbosity > 0) {
      bft_printf("\n");
      cs_mesh_print_info(mesh, _(" After joining"));
      bft_printf("\n");
    }

#if defined(HAVE_MPI)   /* Synchronization */
    if (cs_glob_n_ranks > 1)
      MPI_Barrier(cs_glob_mpi_comm);
#endif

  } /* End of loop on joinings */

  /* Destroy all remaining structures relative to joining operation */

  _destroy_all_joinings();

  full_clock_end = bft_timer_wtime();
  full_cpu_end = bft_timer_cpu_time();

  bft_printf(_("\n"
               "  All joining operations successfully finished:\n"
               "\n"
               "  Time summary:\n"
               "    wall clock time:            %10.3g\n"
               "    CPU time:                   %10.3g\n\n"),
             full_clock_end - full_clock_start,
             full_cpu_end - full_cpu_start);
  bft_printf_flush();

}

/*---------------------------------------------------------------------------*/

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  {
    int  len;
    FILE  *dbg_file = NULL;
    char  *filename = NULL;

    len = strlen("JoinDBG_.dat")+ strlen(local_join_mesh->name) + 4 + 1 + 2;
    BFT_MALLOC(filename, len, char);
    sprintf(filename, "Join%02dDBG_%s%04d.dat",
            join_param.num, local_join_mesh->name, fvm_parall_get_rank());
    dbg_file = fopen(filename, "w");

    cs_join_mesh_dump_file(dbg_file, local_join_mesh);

    fflush(dbg_file);
    BFT_FREE(filename);
    fclose(dbg_file);
  }
#endif

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    cs_debug_glob_mesh_dump("FinalGlobalVertices", mesh);
#endif

END_C_DECLS
