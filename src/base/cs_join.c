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
#include "cs_join_perio.h"
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

/*============================================================================
 * Private function definitions
 *===========================================================================*/

#if 0 && defined(DEBUG) && !defined(NDEBUG)
/*----------------------------------------------------------------------------
 * Dump a cs_mesh_t structure into a file for debugging purpose
 *
 * parameters:
 *   join_num    <-- join number
 *   basename    <-- string
 *   mesh        <-- pointer to cs_mesh_t structure
 *---------------------------------------------------------------------------*/

static void
_dump_mesh(const  int          join_num,
           const  char         basename[],
           const  cs_mesh_t   *mesh)
{
  int  base_len, len;
  FILE  *dbg_file = NULL;
  char  *filename = NULL;

  base_len = strlen(basename);
  len = strlen("JoinDBG__.dat")+1+4+2+base_len;
  BFT_MALLOC(filename, len, char);
  sprintf(filename, "Join%02dDBG_%s_%04d.dat",
          join_num, basename, fvm_parall_get_rank());
  dbg_file = fopen(filename, "w");

  cs_mesh_dump_file(dbg_file, mesh);

  fflush(dbg_file);
  BFT_FREE(filename);
  fclose(dbg_file);
}

/*----------------------------------------------------------------------------
 * Dump a cs_join_gset_t structure into a file for debugging purpose
 *
 * parameters:
 *   join_num    <-- join number
 *   basename    <-- string
 *   set         <-- pointer to cs_join_gset_t structure
 *---------------------------------------------------------------------------*/

static void
_dump_gset(const  int               join_num,
           const  char              basename[],
           const  cs_join_gset_t   *set)
{
  int  base_len, len;
  FILE  *dbg_file = NULL;
  char  *filename = NULL;

  base_len = strlen(basename);
  len = strlen("JoinDBG__.dat")+1+4+2+base_len;
  BFT_MALLOC(filename, len, char);
  sprintf(filename, "Join%02dDBG_%s_%04d.dat",
          join_num, basename, fvm_parall_get_rank());
  dbg_file = fopen(filename, "w");

  cs_join_gset_dump(dbg_file, set);

  fflush(dbg_file);
  BFT_FREE(filename);
  fclose(dbg_file);
}
#endif

/*----------------------------------------------------------------------------
 * Build a structure keeping data about entities selection and modify mesh
 * in case of periodicity.
 *
 * parameters:
 *   this_join  <-- pointer to a cs_join_t structure
 *   mesh       <-> pointer to cs_mesh_t structure
 *---------------------------------------------------------------------------*/

static void
_select_entities(cs_join_t   *this_join,
                 cs_mesh_t   *mesh)
{
  cs_real_t  *b_face_cog = NULL, *b_face_normal = NULL;
  cs_join_param_t   param = this_join->param;

  const char   *selection_criteria = this_join->criteria;

  cs_mesh_init_group_classes(mesh);

  cs_mesh_quantities_b_faces(mesh, &b_face_cog, &b_face_normal);

  cs_glob_mesh->select_b_faces = fvm_selector_create(mesh->dim,
                                                     mesh->n_b_faces,
                                                     mesh->class_defs,
                                                     mesh->b_face_family,
                                                     1,
                                                     b_face_cog,
                                                     b_face_normal);

  /* Get selected faces for this joining and define the related
     cs_join_face_select_t structure.
     - Compute the global number of selected faces
     - Get the adjacent faces, ... */

  this_join->selection = cs_join_select_create(selection_criteria,
                                               param.verbosity);

  /* Free arrays and structures needed for selection */

  BFT_FREE(b_face_cog);
  BFT_FREE(b_face_normal);

  mesh->class_defs = fvm_group_class_set_destroy(mesh->class_defs);

  if (mesh->select_b_faces != NULL)
    mesh->select_b_faces = fvm_selector_destroy(mesh->select_b_faces);
  if (mesh->class_defs != NULL)
    mesh->class_defs = fvm_group_class_set_destroy(mesh->class_defs);

  /* Return selection struct. */

  bft_printf(_("\n  Element selection successfully done.\n"));
  bft_printf_flush();
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
 *   this_join          <-- pointer to a cs_join_t structure
 *   mesh               <-- pointer of pointer to cs_mesh_t structure
 *   p_loc_jmesh        --> local cs_join_mesh_t structure based on local face
 *                          selection
 *   p_work_jmesh       --> distributed and balanced cs_join_mesh_t structure
 *                          based on the global face selection
 *   p_work_join_edges  --> edges definition related to work_jmesh
 *   p_work_face_normal --> unitary normal for the faces of work_jmesh
 *   p_edge_edge_vis    --> list of all potential intersections between edges
 *---------------------------------------------------------------------------*/

static void
_build_join_structures(cs_join_t            *this_join,
                       const cs_mesh_t      *mesh,
                       cs_join_mesh_t     **p_loc_jmesh,
                       cs_join_mesh_t     **p_work_jmesh,
                       cs_join_edges_t    **p_work_join_edges,
                       cs_real_t           *p_work_face_normal[],
                       cs_join_gset_t     **p_edge_edge_vis)
{
  double  clock_start, clock_end, cpu_start, cpu_end;

  char  *mesh_name = NULL;
  cs_real_t  *work_face_normal = NULL;
  cs_join_gset_t  *edge_edge_vis = NULL;
  cs_join_mesh_t  *loc_jmesh = NULL, *work_jmesh = NULL;
  cs_join_edges_t  *work_edges = NULL;
  cs_join_param_t  param = this_join->param;
  cs_join_select_t  *selection = this_join->selection;

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

  loc_jmesh = cs_join_mesh_create_from_select(mesh_name,
                                              param,
                                              selection,
                                              mesh->b_face_vtx_idx,
                                              mesh->b_face_vtx_lst,
                                              mesh->i_face_vtx_idx,
                                              mesh->i_face_vtx_lst,
                                              mesh->n_vertices,
                                              mesh->vtx_coord,
                                              mesh->global_vtx_num);

  if (param.perio_num > 0)
    cs_join_perio_apply(this_join, loc_jmesh, mesh);

  if (param.verbosity > 0)
    cs_join_mesh_minmax_tol(param, loc_jmesh);

  clock_end = bft_timer_wtime();
  cpu_end = bft_timer_cpu_time();

  if (param.verbosity > 2) {

    bft_printf(_("\n    Definition of local joining mesh:\n"
                 "        wall clock time:            %10.3g\n"
                 "        CPU time:                   %10.3g\n"),
               clock_end - clock_start, cpu_end - cpu_start);

    cs_join_post_dump_mesh("LocalMesh", loc_jmesh, param);

  }

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

  _get_work_struct(param,
                   selection->compact_rank_index,
                   loc_jmesh,
                   &work_jmesh,
                   &work_edges,
                   &work_face_normal,
                   &edge_edge_vis);

  clock_end = bft_timer_wtime();
  cpu_end = bft_timer_cpu_time();

  if (param.verbosity > 1)
    bft_printf(_("\n  Definition of structures for the joining algorithm:\n"
                 "      wall clock time:            %10.3g\n"
                 "      CPU time:                   %10.3g\n"),
               clock_end - clock_start, cpu_end - cpu_start);
  bft_printf_flush();

  /* Return pointers */

  *p_loc_jmesh = loc_jmesh;
  *p_work_jmesh = work_jmesh;
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
 *   this_join            <--  pointer to a cs_join_t structure
 *   work_jmesh           <->  pointer to a cs_join_mesh_t structure
 *   work_join_edges      <--  pointer to a cs_join_edges_t structure
 *   p_edge_edge_vis      <->  pointer to a cs_join_glist_t structure
 *                             (freed here)
 *   init_max_vtx_gnum    <--  initial max. global numbering for vertices
 *   p_n_g_new_vertices   -->  global number of vertices created during the
 *                             intersection of edges
 *   p_vtx_eset           -->  structure storing equivalences between vertices
 *                             Two vertices are equivalent if they are each
 *                             other in their tolerance
 *   p_inter_edges        -->  structure storing the definition of new vertices
 *                             on initial edges
 *---------------------------------------------------------------------------*/

static void
_intersect_edges(cs_join_t               *this_join,
                 cs_join_mesh_t          *work_jmesh,
                 const cs_join_edges_t   *work_join_edges,
                 cs_join_gset_t         **p_edge_edge_vis,
                 fvm_gnum_t               init_max_vtx_gnum,
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
  cs_join_param_t  param = this_join->param;

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
                                      work_jmesh,
                                      &vtx_eset,
                                      &inter_set);

  cs_join_gset_destroy(p_edge_edge_vis);

#if 0 && defined(DEBUG) && !defined(NDEBUG) /* Dump structures after inter. */
  cs_join_inter_set_dump(inter_set, work_join_edges, work_jmesh);
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
                                work_jmesh,
                                inter_set,
                                init_max_vtx_gnum,
                                &n_g_new_vertices,
                                &vtx_eset);

    inter_edges = cs_join_inter_edges_define(work_join_edges, inter_set);
    inter_set = cs_join_inter_set_destroy(inter_set);

#if 0 && defined(DEBUG) && !defined(NDEBUG) /* Dump structures after inter. */
    cs_join_inter_edges_dump(inter_edges, work_join_edges, work_jmesh);
#endif

    /* Synchronize inter_edges structure definition */

#if defined(HAVE_MPI)
    if (n_ranks > 1 ) {

      cs_join_inter_edges_t  *sync_block = NULL;

      sync_block = cs_join_inter_edges_part_to_block(work_jmesh,
                                                     work_join_edges,
                                                     inter_edges);

      cs_join_inter_edges_block_to_part(work_join_edges->n_g_edges,
                                        sync_block,
                                        inter_edges);

      /* Add new vertices to edge and mesh description if necessary */

      cs_join_intersect_update_struct(param.verbosity,
                                      work_join_edges,
                                      work_jmesh,
                                      &inter_edges);

      sync_block = cs_join_inter_edges_destroy(sync_block);

      cs_join_mesh_sync_vertices(work_jmesh);

    }
#endif

    /* Find if there are new equivalences between vertices on a same edge */

    cs_join_add_equiv_from_edges(param,
                                 work_jmesh,
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
  cs_join_inter_edges_dump(inter_edges, work_join_edges, work_jmesh);
#endif

}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Retrieve the local new global numbering for the initial vertices from
 * the new vertex global numbering defined by block.
 *
 * parameters:
 *   param              <-- set of user-defined parameter
 *   select             <-- pointer to a cs_join_select_t structure
 *   mesh               <-- pointer of pointer to cs_mesh_t structure
 *   init_max_vtx_gnum  <--  initial max. global numbering for vertices
 *   p_o2n_vtx_gnum     <-> in:  array on blocks on the new global vertex
 *                          out: local array on the new global vertex
 *---------------------------------------------------------------------------*/

static void
_get_local_o2n_vtx_gnum(cs_join_param_t    param,
                        cs_join_select_t  *select,
                        cs_mesh_t         *mesh,
                        fvm_gnum_t         init_max_vtx_gnum,
                        fvm_gnum_t        *p_o2n_vtx_gnum[])
{
  cs_int_t  i, shift, rank;

  cs_int_t  *send_shift = NULL, *recv_shift = NULL;
  cs_int_t  *send_count = NULL, *recv_count = NULL;
  fvm_gnum_t  *send_glist = NULL, *recv_glist = NULL;
  fvm_gnum_t  *new_gnum_by_block = *p_o2n_vtx_gnum;
  fvm_gnum_t  *new_local_gnum = NULL;

  const int  n_ranks = cs_glob_n_ranks;
  const int  local_rank = CS_MAX(cs_glob_rank_id, 0);

  cs_join_block_info_t  block_info = cs_join_get_block_info(init_max_vtx_gnum,
                                                            n_ranks,
                                                            local_rank);

  MPI_Comm  mpi_comm = cs_glob_mpi_comm;

  if (param.perio_num > 0)
    BFT_MALLOC(new_local_gnum, mesh->n_vertices + select->n_vertices,
               fvm_gnum_t);
  else
    BFT_MALLOC(new_local_gnum, mesh->n_vertices, fvm_gnum_t);

  /* Request the new vtx gnum related to the initial vtx gnum */

  BFT_MALLOC(send_count, n_ranks, cs_int_t);
  BFT_MALLOC(recv_count, n_ranks, cs_int_t);

  for (i = 0; i < n_ranks; i++)
    send_count[i] = 0;

  for (i = 0; i < mesh->n_vertices; i++) {
    rank = (mesh->global_vtx_num[i] - 1)/block_info.size;
    send_count[rank] += 1;
  }

  if (param.perio_num > 0) {

    for (i = 0; i < select->n_vertices; i++) {
      rank = (select->per_v_couples[2*i+1] - 1)/block_info.size;
      send_count[rank] += 1;
    }

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

  if (param.perio_num > 0) {

    for (i = 0; i < select->n_vertices; i++) {
      rank = (select->per_v_couples[2*i+1] - 1)/block_info.size;
      shift = send_shift[rank] + send_count[rank];
      send_glist[shift] = select->per_v_couples[2*i+1]; /* Old global num. */
      send_count[rank] += 1;
    }

  }

  MPI_Alltoallv(send_glist, send_count, send_shift, FVM_MPI_GNUM,
                recv_glist, recv_count, recv_shift, FVM_MPI_GNUM,
                mpi_comm);

  /* Send back to the original rank the new global vertex number */

  for (rank = 0; rank < n_ranks; rank++) {

    for (i = recv_shift[rank]; i < recv_shift[rank+1]; i++) {
      shift = recv_glist[i] - block_info.first_gnum;
      recv_glist[i] = new_gnum_by_block[shift];
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
    new_local_gnum[i] = send_glist[shift];  /* New global number */
    send_count[rank] += 1;

  }

  if (param.perio_num > 0) {

    for (i = 0; i < select->n_vertices; i++) {
      rank = (select->per_v_couples[2*i+1] - 1)/block_info.size;
      shift = send_shift[rank] + send_count[rank];
      new_local_gnum[mesh->n_vertices + i] = send_glist[shift]; /* New glob num. */
      send_count[rank] += 1;
    }

  }

  BFT_FREE(send_count);
  BFT_FREE(send_shift);
  BFT_FREE(send_glist);
  BFT_FREE(recv_glist);
  BFT_FREE(recv_count);
  BFT_FREE(recv_shift);
  BFT_FREE(new_gnum_by_block);

  /* Return pointer */

  *p_o2n_vtx_gnum = new_local_gnum;
}

#endif

/*----------------------------------------------------------------------------
 * Define o2n_vtx_gnum for the current rank and in case of periodic apply
 * changes from periodic faces to original faces. Update also per_v_couples
 *
 * parameters:
 *  this_join          <-- pointer to a cs_join_t structure
 *  mesh               <-- pointer to a cs_mesh_t structure
 *  local_jmesh        <-> pointer to a local cs_join_mesh_t structure
 *  p_work_jmesh       <-> distributed join mesh struct. on which operations
 *                         take place
 *   p_work_edges      <-> join edges struct. related to work_jmesh
 *  n_g_new_vertices   <-- global number of vertices created with the
 *                         intersection of edges
 *  init_max_vtx_gnum  <-- initial max. global numbering for vertices
 *  o2n_vtx_gnum       <-> in:  array on blocks on the new global vertex
 *                         out: local array on the new global vertex
 *---------------------------------------------------------------------------*/
static void
_prepare_update_after_merge(cs_join_t          *this_join,
                            cs_mesh_t          *mesh,
                            cs_join_mesh_t     *local_jmesh,
                            cs_join_mesh_t    **p_work_jmesh,
                            cs_join_edges_t   **p_work_edges,
                            fvm_gnum_t          n_g_new_vertices,
                            fvm_gnum_t          init_max_vtx_gnum,
                            fvm_gnum_t         *p_o2n_vtx_gnum[])
{
  int  i, select_id, shift;

  fvm_gnum_t  *o2n_vtx_gnum = *p_o2n_vtx_gnum;
  cs_join_mesh_t  *work_jmesh = *p_work_jmesh;
  cs_join_edges_t  *work_edges = *p_work_edges;
  cs_join_select_t  *selection = this_join->selection;
  cs_join_param_t  param = this_join->param;

  const cs_int_t  n_ranks = cs_glob_n_ranks;

  /* Build an array keeping relation between old/new global vertex num. */

  if (n_ranks == 1) {

    fvm_gnum_t  *loc_vtx_gnum = NULL;

    BFT_MALLOC(loc_vtx_gnum, mesh->n_vertices, fvm_gnum_t);

    /* Initialize array */

    for (i = 0; i < mesh->n_vertices; i++)
      loc_vtx_gnum[i] = i+1;

    /* Update value for selected vertices */

    for (i = 0, select_id = 0;
         i < mesh->n_vertices && select_id < selection->n_vertices; i++) {

      if (i + 1 == selection->vertices[select_id]) /* Is a selected vertex */
        loc_vtx_gnum[i] = o2n_vtx_gnum[select_id++];

    }

    if (param.perio_num > 0) { /* Update per_v_couples */

      assert(selection->n_vertices == selection->n_couples);

      for (i = 0; i < selection->n_vertices; i++) {
        shift = selection->n_vertices + i;
        selection->per_v_couples[2*i] = o2n_vtx_gnum[i];
        selection->per_v_couples[2*i+1] = o2n_vtx_gnum[shift];
      }

    }

    BFT_FREE(o2n_vtx_gnum);
    o2n_vtx_gnum = loc_vtx_gnum; /* Without periodic vertices and for
                                    all vertices in cs_mesh_t structure */

  } /* End if serial mode */

#if defined(HAVE_MPI)
  if (n_ranks > 1) {

    _get_local_o2n_vtx_gnum(param,
                            selection,
                            mesh,
                            init_max_vtx_gnum,
                            &o2n_vtx_gnum);

    if (param.perio_num > 0) {

      for (i = 0; i < selection->n_vertices; i++) {
        select_id = selection->vertices[i] - 1;
        shift = i + mesh->n_vertices;
        selection->per_v_couples[2*i] = o2n_vtx_gnum[select_id];
        selection->per_v_couples[2*i+1] = o2n_vtx_gnum[shift];
      }

      BFT_REALLOC(o2n_vtx_gnum, mesh->n_vertices, fvm_gnum_t);

    }

  }
#endif

  if (param.perio_num > 0)
    cs_join_perio_merge_back(this_join,
                             local_jmesh,
                             &work_jmesh,
                             &work_edges,
                             init_max_vtx_gnum,
                             n_g_new_vertices);

  /* Return pointer */

  *p_o2n_vtx_gnum = o2n_vtx_gnum;
  *p_work_jmesh = work_jmesh;
  *p_work_edges = work_edges;

}

/*----------------------------------------------------------------------------
 * Merge vertices from equivalences found between vertices.
 * Update local and work structures after the merge step.
 *
 * parameters:
 *  this_join            <--  pointer to a cs_join_t structure
 *  n_iwm_vertices       <--  initial number of vertices in work struct.
 *  init_max_vtx_gnum    <--  initial max. global numbering for vertices
 *  n_g_new_vertices     <--  global number of vertices created with the
 *                            intersection of edges
 *  rank_face_gnum_index <--  index on face global numering to determine the
 *                            related rank
 *  vtx_equiv_set        <--  structure storing equivalences between vertices
 *                            Two vertices are equivalent if they are each
 *                            other in their tolerance
 *  inter_edges          <--  structure storing the definition of new vertices
 *                            on initial edges.
 *  p_work_jmesh         <->  pointer to a cs_join_mesh_t structure
 *  p_work_join_edges    <->  pointer to a cs_join_edges_t structure
 *  p_local_jmesh        <->  pointer to a cs_join_mesh_t structure
 *  mesh                 <->  pointer to a cs_mesh_t struct. to update
 *---------------------------------------------------------------------------*/

static void
_merge_vertices(cs_join_t                *this_join,
                cs_int_t                  n_iwm_vertices,
                fvm_gnum_t                init_max_vtx_gnum,
                fvm_gnum_t                n_g_new_vertices,
                cs_join_eset_t           *vtx_eset,
                cs_join_inter_edges_t    *inter_edges,
                cs_join_mesh_t          **p_work_jmesh,
                cs_join_edges_t         **p_work_join_edges,
                cs_join_mesh_t          **p_local_jmesh,
                cs_mesh_t                *mesh)
{
  int  i;
  double  clock_start, clock_end, cpu_start, cpu_end;

  fvm_gnum_t  new_max_vtx_gnum = init_max_vtx_gnum + n_g_new_vertices;
  fvm_gnum_t  *iwm_vtx_gnum = NULL;
  fvm_gnum_t  *o2n_vtx_gnum = NULL;
  cs_join_mesh_t  *local_jmesh = *p_local_jmesh;
  cs_join_mesh_t  *work_jmesh = *p_work_jmesh;
  cs_join_edges_t  *work_join_edges = *p_work_join_edges;
  cs_join_param_t  param = this_join->param;
  cs_join_select_t  *selection = this_join->selection;
  fvm_gnum_t  *rank_face_gnum_index = selection->compact_rank_index;

  assert(local_jmesh != NULL);
  assert(work_jmesh != NULL);
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
    iwm_vtx_gnum[i] = (work_jmesh->vertices[i]).gnum;

  /* Merge vertices */

  cs_join_merge_vertices(param,
                         new_max_vtx_gnum,
                         work_jmesh,
                         vtx_eset);

  cs_join_eset_destroy(&vtx_eset);

  /*  Keep the evolution of vertex global numbering.
      Update work and local structures after vertex merge */

  cs_join_merge_update_struct(param,
                              n_iwm_vertices,
                              iwm_vtx_gnum,
                              init_max_vtx_gnum,
                              rank_face_gnum_index,
                              &work_jmesh,
                              &work_join_edges,
                              &inter_edges,
                              &local_jmesh,
                              &o2n_vtx_gnum); /* Defined by slice in // */

  /* Free memory */

  BFT_FREE(iwm_vtx_gnum);

  inter_edges = cs_join_inter_edges_destroy(inter_edges);

  /* Post if required and level of verbosity is reached */

  if (param.verbosity > 3)
    cs_join_post_dump_mesh("MergeBeforePerioWorkMesh", work_jmesh, param);

  /* Define o2n_vtx_gnum for the current rank and
     apply back periodic transformation if needed */

  _prepare_update_after_merge(this_join,
                              mesh,
                              local_jmesh,
                              &work_jmesh,
                              &work_join_edges,
                              n_g_new_vertices,
                              init_max_vtx_gnum,
                              &o2n_vtx_gnum); /* Defined for the local rank */

  /* Update cs_mesh_t structure after the vertex merge */

  cs_join_update_mesh_after_merge(param,
                                  selection,
                                  o2n_vtx_gnum,       /* free inside */
                                  local_jmesh,
                                  mesh);

  /* Clean meshes after update (empty edges, degenerated edges,  ...) */

  cs_join_mesh_clean(work_jmesh, param.verbosity);
  cs_join_mesh_clean(local_jmesh, param.verbosity);

  /* Define a new cs_join_edges_t structure */

  cs_join_mesh_destroy_edges(&work_join_edges);
  work_join_edges = cs_join_mesh_define_edges(work_jmesh);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  cs_join_mesh_dump_edges(work_join_edges, work_jmesh);
#endif

  clock_end = bft_timer_wtime();
  cpu_end = bft_timer_cpu_time();

  if (param.verbosity > 0)
    bft_printf(_("\n"
                 "  Merge vertices:\n"
                 "    wall clock time:            %10.3g\n"
                 "    CPU time:                   %10.3g\n"),
               clock_end - clock_start, cpu_end - cpu_start);
  bft_printf_flush();

  /* Post if required and level of verbosity is reached */

  if (param.verbosity > 2)
    cs_join_post_dump_mesh("MergeWorkMesh", work_jmesh, param);

  /* Set return pointers */

  *p_local_jmesh = local_jmesh;
  *p_work_jmesh = work_jmesh;
  *p_work_join_edges = work_join_edges;
}

/*----------------------------------------------------------------------------
 * Prepare mesh update after split operation. Invert face history and update
 * local_jmesh in case of periodicity by applying back split operation.
 *
 * parameters:
 *  this_join          <-- pointer to a cs_join_t structure
 *  local_jmesh        <-> pointer to a local cs_join_mesh_t structure
 *  mesh               <-- pointer to a cs_mesh_t structure
 *  p_history          <-> pointer to the history of face splitting
 *                         in: old -> new face links
 *                         out: new -> old face links
 *---------------------------------------------------------------------------*/
static void
_prepare_update_after_split(cs_join_t          *this_join,
                            cs_join_mesh_t     *local_jmesh,
                            cs_mesh_t          *mesh,
                            cs_join_gset_t    **p_history)
{
  cs_join_param_t  param = this_join->param;
  cs_join_gset_t  *o2n_hist = *p_history, *n2o_hist = NULL;

  const cs_int_t  n_ranks = cs_glob_n_ranks;

  /* Invert face historic */

  n2o_hist = cs_join_gset_invert(o2n_hist);

#if defined(HAVE_MPI)
    if (n_ranks > 1) {

      cs_join_gset_t  *n2o_sync_block = NULL;

      MPI_Comm  mpi_comm = cs_glob_mpi_comm;

      n2o_sync_block = cs_join_gset_block_sync(local_jmesh->n_g_faces,
                                               n2o_hist,
                                               mpi_comm);

      cs_join_gset_block_update(local_jmesh->n_g_faces,
                                n2o_sync_block,
                                n2o_hist,
                                mpi_comm);

      cs_join_gset_destroy(&n2o_sync_block);

    }
#endif

  /* Build an array keeping relation between old/new global vertex num. */

  if (param.perio_num > 0)
    cs_join_perio_split_back(this_join,
                             local_jmesh,
                             mesh,
                             o2n_hist,
                             &n2o_hist);

  cs_join_gset_destroy(&o2n_hist);

  /* Return pointer */

  *p_history = n2o_hist;
}

/*----------------------------------------------------------------------------
 * Split faces and update cs_mesh_t structure.
 *
 * parameters:
 *  this_join            <-- pointer to a cs_join_t structure
 *  work_join_edges      <--  pointer to a cs_join_edges_t structure
 *  work_face_normal     <--  normal based on the original face definition
 *  rank_face_gnum_index <--  index on face global numering to determine the
 *                            related rank
 *  p_work_jmesh         <->  pointer to a cs_join_mesh_t structure
 *  local_jmesh          <--  pointer to a cs_join_mesh_t structure
 *  p_mesh               <->  pointer to cs_mesh_t struct.
 *---------------------------------------------------------------------------*/

static void
_split_faces(cs_join_t           *this_join,
             cs_join_edges_t     *work_join_edges,
             fvm_coord_t         *work_face_normal,
             cs_join_mesh_t     **p_work_jmesh,
             cs_join_mesh_t      *local_jmesh,
             cs_mesh_t          **p_mesh)
{
  double  clock_start, clock_end, cpu_start, cpu_end;

  cs_mesh_t  *mesh = *p_mesh;
  cs_join_gset_t  *history = NULL;
  cs_join_param_t  param = this_join->param;
  cs_join_select_t  *selection = this_join->selection;
  fvm_gnum_t  *rank_face_gnum_index = selection->compact_rank_index;

  clock_start = bft_timer_wtime();
  cpu_start = bft_timer_cpu_time();

  cs_join_split_faces(param,
                      work_face_normal,
                      work_join_edges,
                      p_work_jmesh,
                      &history); /* old -> new */

  /* Send back to the original rank the new face description */

  cs_join_split_update_struct(param,
                              *p_work_jmesh,
                              rank_face_gnum_index,
                              &history, /* old -> new */
                              &local_jmesh);

  _prepare_update_after_split(this_join,
                              local_jmesh,
                              mesh,
                              &history); /* in: old->new, out: new -> old */

  /* Update cs_mesh_t structure after the face splitting */

  cs_join_update_mesh_after_split(param,
                                  selection,
                                  history, /* new -> old */
                                  local_jmesh,
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
    cs_join_post_dump_mesh("SplitWorkMesh", *p_work_jmesh, param);

  /* Free memory */

  cs_join_gset_destroy(&history);

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

  return;
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
 * Add a cs_join_t structure to the list of pending joinings.
 *
 * parameters:
 *   join_number  <-- number related to the joining operation
 *   sel_criteria <-- boundary face selection criteria
 *   fraction     <-- value of the fraction parameter
 *   plane        <-- value of the plane parameter
 *   verbosity    <-- level of verbosity required
 *---------------------------------------------------------------------------*/

void
cs_join_add(int     join_number,
            char   *sel_criteria,
            float   fraction,
            float   plane,
            int     verbosity)
{
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

  cs_glob_join_array[cs_glob_n_joinings] = cs_join_create(join_number,
                                                          sel_criteria,
                                                          fraction,
                                                          plane,
                                                          0, /* No perio. */
                                                          verbosity);

  cs_glob_n_joinings++;
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

    cs_join_t  *this_join = cs_glob_join_array[join_id];
    cs_join_param_t  join_param = this_join->param;

    clock_start = bft_timer_wtime();  /* Start timer */
    cpu_start = bft_timer_cpu_time();

    /* Print information into listing file */

    bft_printf(_("\n -------------------------------------------------------\n"
                 "  Joining number %d:\n\n"), join_param.num);

    if (join_param.perio_num > 0)
      bft_printf(_("\n  Joining for periodicity: %d\n"), join_param.perio_num);

    bft_printf(_("  Selection criteria: \"%s\"\n"), this_join->criteria);

    if (join_param.verbosity > 0) {
      bft_printf(_("\n"
                   "  Parameters for the joining operation:\n"
                   "    Shortest incident edge fraction:          %8.5f\n"
                   "    Maximum angle between joined face planes: %8.5f\n\n"),
                 join_param.fraction, join_param.plane);

      bft_printf(_("  Advanced joining parameters:\n"
                   "    Verbosity level:                          %8d\n"
                   "    Deepest level reachable in tree building: %8d\n"
                   "    Max boxes by leaf:                        %8d\n"
                   "    Max ratio of linked boxes / init. boxes:  %8.5f\n"
                   "    Merge step tolerance multiplier:          %8.5f\n"
                   "    Pre-merge factor:                         %8.5f\n"
                   "    Tolerance computation mode:               %8d\n"
                   "    Intersection computation mode:            %8d\n"
                   "    Max. number of equiv. breaks:             %8d\n"
                   "    Max. number of subfaces by face:          %8d\n\n"),
                 join_param.verbosity,
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
    _dump_mesh(join_param.num, "InitMesh", mesh);
#endif

    /* Build arrays and structures required for selection;
       will be destoyed after joining and rebuilt for each new join
       operation in order to take into account mesh modification  */

    _select_entities(this_join, mesh);

    /* Now execute the joining operation */

    if (this_join->selection->n_g_faces > 0) {

      cs_int_t  n_iwm_vertices;      /* iwm: initial work mesh */
      fvm_gnum_t  init_max_vtx_gnum, n_g_new_vertices;

      cs_real_t  *work_face_normal = NULL;
      cs_join_gset_t  *edge_edge_visibility = NULL;
      cs_join_mesh_t  *work_jmesh = NULL, *local_jmesh = NULL;
      cs_join_edges_t  *work_join_edges = NULL;
      cs_join_eset_t  *vtx_eset = NULL;
      cs_join_inter_edges_t  *inter_edges = NULL;

      _build_join_structures(this_join,
                             mesh,
                             &local_jmesh,
                             &work_jmesh,
                             &work_join_edges,
                             &work_face_normal,
                             &edge_edge_visibility);

      if (join_param.verbosity > 2)
        bft_printf(_("\n  Number of faces to treat locally: %10d\n"),
                   work_jmesh->n_faces);

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

      n_iwm_vertices = work_jmesh->n_vertices;

      init_max_vtx_gnum = mesh->n_g_vertices;
      if (join_param.perio_num > 0)
        init_max_vtx_gnum += this_join->selection->n_g_vertices;

      _intersect_edges(this_join,
                       work_jmesh,
                       work_join_edges,
                       &edge_edge_visibility, /* free during this step */
                       init_max_vtx_gnum,
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

      _merge_vertices(this_join,
                      n_iwm_vertices,
                      init_max_vtx_gnum,
                      n_g_new_vertices,
                      vtx_eset,         /* free during this step */
                      inter_edges,      /* free during this step */
                      &work_jmesh,
                      &work_join_edges,
                      &local_jmesh,
                      mesh);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    _dump_mesh(join_param.num, "MeshAfterMerge", mesh);
#endif

      /*
         Split faces in work_jmesh. Apply modification to the
         local_jmesh. Keep a history between old --> new faces.
         Update cs_mesh_t structure.
      */

      _split_faces(this_join,
                   work_join_edges,
                   work_face_normal,
                   &work_jmesh,
                   local_jmesh,
                   &mesh);

      /* Free memory */

      cs_join_mesh_destroy(&local_jmesh);
      cs_join_mesh_destroy(&work_jmesh);
      cs_join_mesh_destroy_edges(&work_join_edges);

      BFT_FREE(work_face_normal);

      /* Clean mesh (delete redundant edge definition) */

      cs_join_update_mesh_clean(join_param, mesh);

    }
    else
      bft_printf(_("\nStop joining algorithm: no face selected...\n"));

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    _dump_mesh(join_param.num, "FinalMesh", mesh);
#endif

    if (join_param.verbosity > 0) {
      bft_printf("\n");
      cs_mesh_print_info(mesh, _(" After joining"));
      bft_printf("\n");
    }

    clock_end = bft_timer_wtime();
    cpu_end = bft_timer_cpu_time();

    bft_printf(_("\n"
                 "  Complete joining treatment for joining %2d\n"
                 "    wall clock time:            %10.3g\n"
                 "    CPU time:                   %10.3g\n"),
               join_param.num, clock_end - clock_start, cpu_end - cpu_start);
    bft_printf_flush();

    /* Free memory */

    cs_join_destroy(&this_join);

#if defined(HAVE_MPI)   /* Synchronization */
    if (cs_glob_n_ranks > 1)
      MPI_Barrier(cs_glob_mpi_comm);
#endif

  } /* End of loop on joinings */

  /* Destroy all remaining structures relative to joining operation */

  cs_glob_join_count += cs_glob_n_joinings;
  cs_glob_n_joinings = 0;

  BFT_FREE(cs_glob_join_array);

  if (cs_glob_n_join_perio > 0)
    cs_join_perio_transfer_builder();

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
    cs_debug_glob_mesh_dump("FinalGlobalVertices", mesh);
#endif

END_C_DECLS
