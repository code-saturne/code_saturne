/*============================================================================
 * Management of conforming and non-conforming joining
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

#include <assert.h>
#include <errno.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *---------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_periodicity.h"

#include "cs_all_to_all.h"
#include "cs_block_dist.h"
#include "cs_join_intersect.h"
#include "cs_join_merge.h"
#include "cs_join_mesh.h"
#include "cs_join_post.h"
#include "cs_join_perio.h"
#include "cs_join_set.h"
#include "cs_join_split.h"
#include "cs_join_update.h"
#include "cs_join_util.h"
#include "cs_log.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *---------------------------------------------------------------------------*/

#include "cs_join.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

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
  len = strlen("log/JoinDBG__.dat")+1+4+2+base_len;
  BFT_MALLOC(filename, len, char);
  sprintf(filename, "log%cJoin%02dDBG_%s_%04d.dat", CS_DIR_SEPARATOR,
          join_num, basename, cs_glob_rank_id);
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
  len = strlen("log/JoinDBG__.dat")+1+4+2+base_len;
  BFT_MALLOC(filename, len, char);
  sprintf(filename, "log%cJoin%02dDBG_%s_%04d.dat", CS_DIR_SEPARATOR,
          join_num, basename, cs_glob_rank_id);
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

  if (mesh->verbosity > 0) {
    bft_printf(_("\n  Element selection successfully done.\n"));
    bft_printf_flush();
  }
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
 *   stats                <-> joining statistics
 *---------------------------------------------------------------------------*/

static void
_get_work_struct(cs_join_param_t         param,
                 const cs_gnum_t         rank_face_gnum_index[],
                 const cs_join_mesh_t   *local_mesh,
                 cs_join_mesh_t        **p_work_mesh,
                 cs_join_edges_t       **p_work_edges,
                 cs_real_t              *p_work_face_normal[],
                 cs_join_gset_t        **p_edge_edge_vis,
                 cs_join_stats_t        *stats)
{
  cs_lnum_t  n_inter_faces = 0;
  char  *mesh_name = NULL;
  cs_real_t  *face_normal = NULL;
  cs_gnum_t  *intersect_face_gnum = NULL;
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

  face_face_vis = cs_join_intersect_faces(param, local_mesh, stats);

  /* Define an ordered list of all implied faces without redundancy */

  /* TODO: check if this is necessary after cleanup done in face_face_vis */

  cs_timer_t t0 = cs_timer_time();

  cs_join_gset_single_order(face_face_vis,
                            &n_inter_faces,
                            &intersect_face_gnum);

  cs_timer_t t1 = cs_timer_time();

  if (param.verbosity > 1)
    bft_printf(_("\n  Sorted possible intersections between faces.\n"));

  cs_timer_counter_add_diff(&(stats->t_inter_sort), &t0, &t1);

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
 *   mesh               <-> pointer of pointer to cs_mesh_t structure
 *   p_loc_jmesh        --> local cs_join_mesh_t structure based on local face
 *                          selection
 *   p_work_jmesh       --> distributed and balanced cs_join_mesh_t structure
 *                          based on the global face selection
 *   p_work_join_edges  --> edges definition related to work_jmesh
 *   p_work_face_normal --> unitary normal for the faces of work_jmesh
 *   p_edge_edge_vis    --> list of all potential intersections between edges
 *---------------------------------------------------------------------------*/

static void
_build_join_structures(cs_join_t           *this_join,
                       cs_mesh_t           *mesh,
                       cs_join_mesh_t     **p_loc_jmesh,
                       cs_join_mesh_t     **p_work_jmesh,
                       cs_join_edges_t    **p_work_join_edges,
                       cs_real_t           *p_work_face_normal[],
                       cs_join_gset_t     **p_edge_edge_vis)
{
  char  *mesh_name = NULL;
  cs_real_t  *work_face_normal = NULL;
  cs_join_gset_t  *edge_edge_vis = NULL;
  cs_join_mesh_t  *loc_jmesh = NULL, *work_jmesh = NULL;
  cs_join_edges_t  *work_edges = NULL;
  cs_join_param_t  param = this_join->param;
  cs_join_select_t  *selection = this_join->selection;

  cs_timer_t t0 = cs_timer_time();

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

  BFT_FREE(mesh_name);

  if (param.perio_type != FVM_PERIODICITY_NULL)
    cs_join_perio_apply(this_join, loc_jmesh, mesh);

  if (param.verbosity > 0)
    cs_join_mesh_minmax_tol(param, loc_jmesh);

  cs_timer_t t1 = cs_timer_time();

  if (param.visualization > 2)
    cs_join_post_dump_mesh("LocalMesh", loc_jmesh, param);

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
                   &edge_edge_vis,
                   &(this_join->stats));

  /* log performance info of previous step here only to simplify
     "pretty printing". */

  cs_timer_counter_add_diff(&(this_join->stats.t_l_join_mesh), &t0, &t1);

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
 *
 * returns:
 *   type of joining detected
 *---------------------------------------------------------------------------*/

static cs_join_type_t
_intersect_edges(cs_join_t               *this_join,
                 cs_join_mesh_t          *work_jmesh,
                 const cs_join_edges_t   *work_join_edges,
                 cs_join_gset_t         **p_edge_edge_vis,
                 cs_gnum_t                init_max_vtx_gnum,
                 cs_gnum_t               *p_n_g_new_vertices,
                 cs_join_eset_t         **p_vtx_eset,
                 cs_join_inter_edges_t  **p_inter_edges)
{
  cs_join_type_t  join_type = CS_JOIN_TYPE_NULL;

  cs_gnum_t  n_g_new_vertices = 0;
  cs_join_inter_edges_t  *inter_edges = NULL;
  cs_join_eset_t  *vtx_eset = NULL;
  cs_join_inter_set_t  *inter_set = NULL;
  cs_join_param_t  param = this_join->param;

  const int  n_ranks = cs_glob_n_ranks;

  cs_timer_t t0 = cs_timer_time();

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

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(this_join->stats.t_edge_inter), &t0, &t1);
  t0 = t1;

  cs_join_gset_destroy(p_edge_edge_vis);

#if 0 && defined(DEBUG) && !defined(NDEBUG) /* Dump structures after inter. */
  cs_join_inter_set_dump(inter_set, work_join_edges, work_jmesh);
#endif

  if (join_type == CS_JOIN_TYPE_NULL) {

    if (cs_glob_mesh->verbosity > 0) {
      bft_printf(_("\n  Joining operation did modify any vertices.\n"));
      bft_printf_flush();
    }

    cs_join_inter_set_destroy(&inter_set);

  }
  else if (join_type == CS_JOIN_TYPE_CONFORMING) {

    if (cs_glob_mesh->verbosity > 0) {
      bft_printf(_("\n  Joining operation is conforming.\n"));
      bft_printf_flush();
    }

    cs_join_inter_set_destroy(&inter_set);

  }
  else {

    assert(join_type == CS_JOIN_TYPE_NON_CONFORMING);

    if (cs_glob_mesh->verbosity > 0) {
      bft_printf(_("\n  Joining operation is non-conforming.\n"));
      bft_printf_flush();
    }

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
    cs_join_inter_set_destroy(&inter_set);

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

      cs_join_inter_edges_destroy(&sync_block);

      cs_join_mesh_sync_vertices(work_jmesh);

    }
#endif

    /* Find if there are new equivalences between vertices on a same edge */

    cs_join_add_equiv_from_edges(param,
                                 work_jmesh,
                                 work_join_edges,
                                 inter_edges,
                                 vtx_eset);

  } /* non conforming joining operation */

  /* Order and delete redundant equivalences */

  cs_join_eset_clean(&vtx_eset);

  /* Memory management: final state for vtx_eset (no more equiv. to get) */

  vtx_eset->n_max_equiv = vtx_eset->n_equiv;
  BFT_REALLOC(vtx_eset->equiv_couple, 2*vtx_eset->n_equiv, cs_lnum_t);

  t1 = cs_timer_time();

  cs_timer_counter_add_diff(&(this_join->stats.t_new_vtx), &t0, &t1);

  if (param.verbosity > 0)
    bft_printf(_("\n"
                 "  Edge intersections and vertex creation done.\n"));

  /* Returns pointers */

  *p_vtx_eset = vtx_eset;
  *p_inter_edges = inter_edges;
  *p_n_g_new_vertices = n_g_new_vertices;

#if 0 && defined(DEBUG) && !defined(NDEBUG) /* Dump structures after inter. */
  cs_join_inter_edges_dump(inter_edges, work_join_edges, work_jmesh);
#endif

  return join_type;
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
                        cs_gnum_t          init_max_vtx_gnum,
                        cs_gnum_t         *p_o2n_vtx_gnum[])
{
  cs_lnum_t n_tot_vertices = mesh->n_vertices;

  int *dest_rank = NULL;
  cs_gnum_t *new_gnum_by_block = *p_o2n_vtx_gnum;
  cs_gnum_t *new_local_gnum = NULL;
  cs_all_to_all_t *d = NULL;

  const int  n_ranks = cs_glob_n_ranks;
  const int  local_rank = cs_glob_rank_id;

  cs_block_dist_info_t  bi = cs_block_dist_compute_sizes(local_rank,
                                                         n_ranks,
                                                         1,
                                                         0,
                                                         init_max_vtx_gnum);

  MPI_Comm  comm = cs_glob_mpi_comm;

  if (param.perio_type != FVM_PERIODICITY_NULL)
    n_tot_vertices = mesh->n_vertices + select->n_vertices;

  /* Initialize new vertex gnum with old one */

  BFT_MALLOC(new_local_gnum, n_tot_vertices, cs_gnum_t);
  BFT_MALLOC(dest_rank, n_tot_vertices, int);

  for (cs_lnum_t i = 0; i < mesh->n_vertices; i++) {
    dest_rank[i] = (mesh->global_vtx_num[i] - 1)/(cs_gnum_t)(bi.block_size);
    assert(dest_rank[i] >= 0 && dest_rank[i] < n_ranks);
    new_local_gnum[i] = mesh->global_vtx_num[i];
  }

  if (param.perio_type != FVM_PERIODICITY_NULL) {
    for (cs_lnum_t i = 0; i < select->n_vertices; i++) {
      const cs_lnum_t j = mesh->n_vertices + i;
      dest_rank[j] = (  (select->per_v_couples[2*i+1] - 1)
                      / (cs_gnum_t)(bi.block_size)) * bi.rank_step;
      assert(dest_rank[j] >= 0 && dest_rank[j] < n_ranks);
      new_local_gnum[j] = select->per_v_couples[2*i+1];
    }
  }

  /* Send old vtx gnum to matching block */

  d = cs_all_to_all_create(n_tot_vertices,
                           0,     /* flags */
                           NULL,  /* dest_id */
                           dest_rank,
                           comm);

  cs_all_to_all_transfer_dest_rank(d, &dest_rank);

  cs_gnum_t *b_data = cs_all_to_all_copy_array(d,
                                               CS_GNUM_TYPE,
                                               1,
                                               false,  /* reverse */
                                               new_local_gnum,
                                               NULL);

  /* Request the new vtx gnum related to the initial vtx gnum */

  cs_lnum_t n_b_vtx = cs_all_to_all_n_elts_dest(d);

  for (cs_lnum_t i = 0; i < n_b_vtx; i++) {

    /* Transform old to new vertex number */
    cs_gnum_t o_shift = b_data[i] - bi.gnum_range[0];
    b_data[i] = new_gnum_by_block[o_shift];

  }

  /* Send data back */

  cs_all_to_all_copy_array(d,
                           CS_GNUM_TYPE,
                           1,
                           true,  /* reverse */
                           b_data,
                           new_local_gnum);

  BFT_FREE(b_data);

  cs_all_to_all_destroy(&d);

  /* Return pointer */

  BFT_FREE(new_gnum_by_block);
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
                            cs_gnum_t           n_g_new_vertices,
                            cs_gnum_t           init_max_vtx_gnum,
                            cs_gnum_t          *p_o2n_vtx_gnum[])
{
  int  i, select_id, shift;

  cs_gnum_t  *o2n_vtx_gnum = *p_o2n_vtx_gnum;
  cs_join_mesh_t  *work_jmesh = *p_work_jmesh;
  cs_join_edges_t  *work_edges = *p_work_edges;
  cs_join_select_t  *selection = this_join->selection;
  cs_join_param_t  param = this_join->param;

  const cs_lnum_t  n_ranks = cs_glob_n_ranks;

  /* Build an array keeping relation between old/new global vertex num. */

  if (n_ranks == 1) {

    cs_gnum_t  *loc_vtx_gnum = NULL;

    BFT_MALLOC(loc_vtx_gnum, mesh->n_vertices, cs_gnum_t);

    /* Initialize array */

    for (i = 0; i < mesh->n_vertices; i++)
      loc_vtx_gnum[i] = i+1;

    /* Update value for selected vertices */

    for (i = 0, select_id = 0;
         i < mesh->n_vertices && select_id < selection->n_vertices; i++) {

      if (i + 1 == selection->vertices[select_id]) /* Is a selected vertex */
        loc_vtx_gnum[i] = o2n_vtx_gnum[select_id++];

    }

    if (param.perio_type != FVM_PERIODICITY_NULL) { /* Update per_v_couples */

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

    if (param.perio_type != FVM_PERIODICITY_NULL) {

      for (i = 0; i < selection->n_vertices; i++) {
        select_id = selection->vertices[i] - 1;
        shift = i + mesh->n_vertices;
        selection->per_v_couples[2*i] = o2n_vtx_gnum[select_id];
        selection->per_v_couples[2*i+1] = o2n_vtx_gnum[shift];
      }

      BFT_REALLOC(o2n_vtx_gnum, mesh->n_vertices, cs_gnum_t);

    }

  }
#endif

  if (param.perio_type != FVM_PERIODICITY_NULL)
    cs_join_perio_merge_back(this_join,
                             local_jmesh,
                             mesh,
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
 *   this_join            <--  pointer to a cs_join_t structure
 *   n_iwm_vertices       <--  initial number of vertices in work struct.
 *   init_max_vtx_gnum    <--  initial max. global numbering for vertices
 *   n_g_new_vertices     <--  global number of vertices created with the
 *                             intersection of edges
 *   vtx_eset             <->  structure storing equivalences between vertices
 *                             Two vertices are equivalent if they are each
 *                             other in their tolerance; freed here after use
 *   inter_edges          <->  structure storing the definition of new vertices
 *                             on initial edges; freed here after use
 *   p_work_jmesh         <->  pointer to a cs_join_mesh_t structure
 *   p_work_join_edges    <->  pointer to a cs_join_edges_t structure
 *   p_local_jmesh        <->  pointer to a cs_join_mesh_t structure
 *   mesh                 <->  pointer to a cs_mesh_t struct. to update
 *---------------------------------------------------------------------------*/

static void
_merge_vertices(cs_join_t                *this_join,
                cs_lnum_t                 n_iwm_vertices,
                cs_gnum_t                 init_max_vtx_gnum,
                cs_gnum_t                 n_g_new_vertices,
                cs_join_eset_t          **vtx_eset,
                cs_join_inter_edges_t   **inter_edges,
                cs_join_mesh_t          **p_work_jmesh,
                cs_join_edges_t         **p_work_join_edges,
                cs_join_mesh_t          **p_local_jmesh,
                cs_mesh_t                *mesh)
{
  int  i;

  cs_gnum_t  new_max_vtx_gnum = init_max_vtx_gnum + n_g_new_vertices;
  cs_gnum_t  *iwm_vtx_gnum = NULL;
  cs_gnum_t  *o2n_vtx_gnum = NULL;
  cs_join_mesh_t  *local_jmesh = *p_local_jmesh;
  cs_join_mesh_t  *work_jmesh = *p_work_jmesh;
  cs_join_edges_t  *work_join_edges = *p_work_join_edges;
  cs_join_param_t  param = this_join->param;
  cs_join_select_t  *selection = this_join->selection;
  cs_gnum_t  *rank_face_gnum_index = selection->compact_rank_index;

  assert(local_jmesh != NULL);
  assert(work_jmesh != NULL);
  assert(work_join_edges != NULL);

  cs_timer_t t0 = cs_timer_time();

  /*
    Store the initial global vertex numbering
    Initial vertices are between [0, n_init_vertices[
    Added vertices from inter. are between [n_init_vertices, n_vertices]
  */

  BFT_MALLOC(iwm_vtx_gnum, n_iwm_vertices, cs_gnum_t);

  for (i = 0; i < n_iwm_vertices; i++)
    iwm_vtx_gnum[i] = (work_jmesh->vertices[i]).gnum;

  /* Merge vertices */

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(this_join->stats.t_u_merge_vtx), &t0, &t1);
  t0 = t1;

  cs_join_merge_vertices(param,
                         new_max_vtx_gnum,
                         work_jmesh,
                         *vtx_eset);

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(this_join->stats.t_merge_vtx), &t0, &t1);
  t0 = t1;

  cs_join_eset_destroy(vtx_eset);

  /*  Keep the evolution of vertex global numbering.
      Update work and local structures after vertex merge */

  cs_join_merge_update_struct(param,
                              n_iwm_vertices,
                              iwm_vtx_gnum,
                              init_max_vtx_gnum,
                              rank_face_gnum_index,
                              &work_jmesh,
                              &work_join_edges,
                              inter_edges,
                              &local_jmesh,
                              &o2n_vtx_gnum); /* Defined by slice in // */

  /* Free memory */

  BFT_FREE(iwm_vtx_gnum);

  cs_join_inter_edges_destroy(inter_edges);

  /* Post if required and level of verbosity is reached */

  if (param.visualization > 3)
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

  /* Clean meshes after update (empty edges, degenerate edges,  ...) */

  cs_join_mesh_clean(work_jmesh, param.verbosity);
  cs_join_mesh_clean(local_jmesh, param.verbosity);

  /* Define a new cs_join_edges_t structure */

  cs_join_mesh_destroy_edges(&work_join_edges);
  work_join_edges = cs_join_mesh_define_edges(work_jmesh);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  cs_join_mesh_dump_edges(work_join_edges, work_jmesh);
#endif

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(this_join->stats.t_u_merge_vtx), &t0, &t1);

  if (param.verbosity > 0) {
    bft_printf(_("\n"
                 "  Merge vertices and mesh update done.\n"));
    bft_printf_flush();
  }

  /* Post if required and level of verbosity is reached */

  if (param.visualization > 2)
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
 *  mesh_builder       <-- pointer to a cs_mesh_builder_t structure
 *  p_history          <-> pointer to the history of face splitting
 *                         in: old -> new face links
 *                         out: new -> old face links
 *---------------------------------------------------------------------------*/
static void
_prepare_update_after_split(cs_join_t          *this_join,
                            cs_join_mesh_t     *local_jmesh,
                            cs_mesh_t          *mesh,
                            cs_mesh_builder_t  *mesh_builder,
                            cs_join_gset_t    **p_history)
{
  cs_join_param_t  param = this_join->param;
  cs_join_gset_t  *o2n_hist = *p_history, *n2o_hist = NULL;

  const cs_lnum_t  n_ranks = cs_glob_n_ranks;

  /* Invert face historic */

  n2o_hist = cs_join_gset_invert(o2n_hist);

#if defined(HAVE_MPI)
    if (n_ranks > 1) {

      cs_join_gset_t  *n2o_sync_block = NULL;

      MPI_Comm  comm = cs_glob_mpi_comm;

      n2o_sync_block = cs_join_gset_block_sync(local_jmesh->n_g_faces,
                                               n2o_hist,
                                               comm);

      cs_join_gset_block_update(local_jmesh->n_g_faces,
                                n2o_sync_block,
                                n2o_hist,
                                comm);

      cs_join_gset_destroy(&n2o_sync_block);

    }
#endif

  /* Build an array keeping relation between old/new global vertex num. */

  if (param.perio_type != FVM_PERIODICITY_NULL)
    cs_join_perio_split_back(this_join,
                             local_jmesh,
                             mesh,
                             mesh_builder,
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
 *  join_type            <-> join type as detected so far
 *  work_join_edges      <-- pointer to a cs_join_edges_t structure
 *  work_face_normal     <-- normal based on the original face definition
 *  rank_face_gnum_index <-- index on face global numering to determine the
 *                           related rank
 *  p_work_jmesh         <-> pointer to a cs_join_mesh_t structure
 *  local_jmesh          <-- pointer to a cs_join_mesh_t structure
 *  mesh                 <-> pointer to cs_mesh_t structure
 *  mesh_builder         <-> pointer to cs_mesh_builder structure
 *---------------------------------------------------------------------------*/

static void
_split_faces(cs_join_t           *this_join,
             cs_join_type_t      *join_type,
             cs_join_edges_t     *work_join_edges,
             cs_coord_t          *work_face_normal,
             cs_join_mesh_t     **p_work_jmesh,
             cs_join_mesh_t      *local_jmesh,
             cs_mesh_t           *mesh,
             cs_mesh_builder_t   *mesh_builder)
{
  cs_join_gset_t  *history = NULL;
  cs_join_param_t  param = this_join->param;
  cs_join_select_t  *selection = this_join->selection;
  cs_gnum_t  *rank_face_gnum_index = selection->compact_rank_index;

  cs_timer_t t0 = cs_timer_time();

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
                              mesh_builder,
                              &history); /* in: old->new, out: new -> old */

  /* Update cs_mesh_t structure after the face splitting */

  cs_join_update_mesh_after_split(param,
                                  selection,
                                  history, /* new -> old */
                                  local_jmesh,
                                  mesh,
                                  mesh_builder);

  cs_timer_t t1 = cs_timer_time();

  cs_timer_counter_add_diff(&(this_join->stats.t_split_faces), &t0, &t1);

  bft_printf_flush();

  /* Post if required and level of verbosity is reached */

  if (param.visualization > 2)
    cs_join_post_dump_mesh("SplitWorkMesh", *p_work_jmesh, param);

  /* Check join type if required */

  if (cs_glob_mesh->verbosity > 0 && *join_type == CS_JOIN_TYPE_NULL) {

    int _join_type = CS_JOIN_TYPE_NULL;
    for (cs_lnum_t i = 0; i < history->n_elts; i++) {
      if (history->index[i+1] != i+1) {
        _join_type = (int)CS_JOIN_TYPE_NON_CONFORMING;
        break;
      }
      else if (local_jmesh->face_gnum[i] != history->g_elts[i]) {
        _join_type = (int)CS_JOIN_TYPE_CONFORMING;
        break;
      }
    }
    cs_parall_max(1, CS_INT_TYPE, &_join_type);

    *join_type = _join_type;

    if (*join_type == CS_JOIN_TYPE_NULL)
      bft_printf(_("\n  Joining operation is null.\n"));
    else if (*join_type == CS_JOIN_TYPE_CONFORMING)
      bft_printf(_("\n  Joining operation is conforming.\n"));
    else if (*join_type == CS_JOIN_TYPE_NON_CONFORMING)
      bft_printf(_("\n  Joining operation is non-conforming.\n"));
    bft_printf_flush();

  }

  /* Free memory */

  cs_join_gset_destroy(&history);
}

/*----------------------------------------------------------------------------
 * Print information on a mesh structure.
 *
 * parameters:
 *   mesh  <--  pointer to mesh structure.
 *   name  <--  associated name.
 *----------------------------------------------------------------------------*/

static void
_print_mesh_info(const cs_mesh_t  *mesh,
                 const char       *name)
{
  bft_printf(_(" %s\n"
               "     Number of cells:          %llu\n"
               "     Number of interior faces: %llu\n"
               "     Number of boundary faces: %llu\n"
               "     Number of vertices:       %llu\n"),
             name,
             (unsigned long long)(mesh->n_g_cells),
             (unsigned long long)(mesh->n_g_i_faces),
             (unsigned long long)(mesh->n_g_b_faces),
             (unsigned long long)(mesh->n_g_vertices));
}

/*----------------------------------------------------------------------------
 * Print initial joining info into log file
 *---------------------------------------------------------------------------*/

static void
_print_join_info(cs_mesh_t  *mesh,
                 cs_join_t  *this_join,
                 cs_join_param_t join_param)
{
  bft_printf(_("\n -------------------------------------------------------\n"
               "  Joining number %d:\n\n"), join_param.num);

  if (join_param.perio_type != FVM_PERIODICITY_NULL) {
    double m[3][4];
    memcpy(m, join_param.perio_matrix, sizeof(double)*12);
    bft_printf(_("  Periodicity type: %s\n"),
               _(fvm_periodicity_type_name[join_param.perio_type]));
    bft_printf(_("  Transformation matrix:  %12.5g %12.5g %12.5g %12.5g\n"
                 "                          %12.5g %12.5g %12.5g %12.5g\n"
                 "                          %12.5g %12.5g %12.5g %12.5g\n\n"),
               m[0][0], m[0][1], m[0][2], m[0][3],
               m[1][0], m[1][1], m[1][2], m[1][3],
               m[2][0], m[2][1], m[2][2], m[2][3]);
  }

  bft_printf(_("  Selection criteria: \"%s\"\n"), this_join->criteria);

  if (join_param.verbosity > 0) {
    bft_printf(_("\n"
                 "  Parameters for the joining operation:\n"
                 "    Shortest incident edge fraction:          %8.5f\n"
                 "    Maximum angle between joined face planes: %8.5f\n\n"),
               join_param.fraction, join_param.plane);

    bft_printf(_("  Advanced joining parameters:\n"
                 "    Verbosity level:                          %8d\n"
                 "    Visualization level:                      %8d\n"
                 "    Deepest level reachable in tree building: %8d\n"
                 "    Max boxes by leaf:                        %8d\n"
                 "    Max ratio of linked boxes / init. boxes:  %8.5f\n"
                 "    Max ratio of boxes for distribution:      %8.5f\n"
                 "    Merge step tolerance multiplier:          %8.5f\n"
                 "    Pre-merge factor:                         %8.5f\n"
                 "    Tolerance computation mode:               %8d\n"
                 "    Intersection computation mode:            %8d\n"
                 "    Max. number of equiv. breaks:             %8d\n"
                 "    Max. number of subfaces by face:          %8d\n\n"),
               join_param.verbosity,
               join_param.visualization,
               join_param.tree_max_level,
               join_param.tree_n_max_boxes,
               join_param.tree_max_box_ratio,
               join_param.tree_max_box_ratio_distrib,
               join_param.merge_tol_coef,
               join_param.pre_merge_factor,
               join_param.tcm, join_param.icm,
               join_param.n_max_equiv_breaks,
               join_param.max_sub_faces);

    _print_mesh_info(mesh, _(" Before joining"));
    bft_printf("\n");
  }
}

/*----------------------------------------------------------------------------
 * Set advanced parameters to user-defined values.
 *
 * Out-of range values are silently set to minimum acceptable values
 * where those are possible.
 *
 * parameters:
 *   join           <-> pointer to a cs_join_t structure to update
 *   mtf            <-- merge tolerance coefficient
 *   pmf            <-- pre-merge factor
 *   tcm            <-- tolerance computation mode
 *   icm            <-- intersection computation mode
 *   maxbrk         <-- max number of equivalences to break (merge step)
 *   max_sub_faces  <-- max. possible number of sub-faces by splitting a face
 *   tml            <-- tree max level
 *   tmb            <-- tree max boxes
 *   tmr            <-- tree max ratio
 *   tmr_distrib    <-- tree max ratio for distribution
 *---------------------------------------------------------------------------*/

static void
_set_advanced_param(cs_join_t   *join,
                    double       mtf,
                    double       pmf,
                    int          tcm,
                    int          icm,
                    int          maxbrk,
                    int          max_sub_faces,
                    int          tml,
                    int          tmb,
                    double       tmr,
                    double       tmr_distrib)
{
  /* Deepest level reachable during tree building */

  if (tml < 1)
    tml = 1;

  join->param.tree_max_level = tml;

  /* Max. number of boxes which can be related to a leaf of the tree
     if level != tree_max_level */

  if (tmb < 1)
    tmb = 1;

  join->param.tree_n_max_boxes = tmb;

  /* Stop tree building if:
     n_linked_boxes > tree_max_box_ratio*n_init_boxes */

  if (tmr < 1.0 )
    tmr = 1.0;

  join->param.tree_max_box_ratio = tmr;

  if (tmr_distrib < 1.0 )
    tmr_distrib = 1.0;

  join->param.tree_max_box_ratio_distrib = tmr_distrib;

  /* Coef. used to modify the tolerance associated to each vertex BEFORE the
     merge operation.
     If coef = 0.0 => no vertex merge
     If coef < 1.0 => reduce vertex merge
     If coef = 1.0 => no change
     If coef > 1.0 => increase vertex merge */

  if (mtf < 0.0)
    mtf = 0.0;

  join->param.merge_tol_coef = mtf;

   /* Maximum number of equivalence breaks */

  if (maxbrk < 0)
    maxbrk = 0;

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
                "  It must be 1, 2, 11, or 12 and not: %d\n"), tcm);

  join->param.tcm = tcm;

  /* Intersection computation mode */

  if (icm != 1 && icm != 2)
    bft_error(__FILE__, __LINE__, 0,
              _("Mesh joining:"
                "  Forbidden value for icm parameter.\n"
                "  It must be 1 or 2 and not: %d\n"), icm);

  join->param.icm = icm;

  /* Maximum number of sub-faces */

  if (max_sub_faces < 1)
    bft_error(__FILE__, __LINE__, 0,
              _("Mesh joining:"
                "  Forbidden value for the maxsf parameter.\n"
                "  It must be > 0 and not: %d\n"), max_sub_faces);

  join->param.max_sub_faces = max_sub_faces;

}

/*----------------------------------------------------------------------------
 * Log statistics and timings for a given joining.
 *
 * parameters:
 *   join   <-- pointer to a cs_join_t structure
 *---------------------------------------------------------------------------*/

static void
_join_performance_log(const cs_join_t  *this_join)
{
  char buf[80];

  const cs_join_stats_t  *stats = &(this_join->stats);

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\nJoining number %d:\n\n"), this_join->param.num);

  if (this_join->stats.n_calls > 1)
    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("\n  Number of calls (statistics are cumulative): %d:\n\n"),
                  this_join->stats.n_calls);

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  Determination of possible face intersections:\n\n"
                  "    bounding-box tree layout: %dD\n"), stats->bbox_layout);

  if (cs_glob_n_ranks > 1 || stats->n_calls > 1) {

    cs_gnum_t n = CS_MAX(stats->n_calls, 1);

    if (stats->n_calls > 1)
      strncpy(buf, _("                                   rank mean"), 79);
    else if (cs_glob_n_ranks <= 1)
      strncpy(buf, _("                                   call mean"), 79);
    else
      strncpy(buf, _("                              rank/call mean"), 79);

    cs_log_printf
      (CS_LOG_PERFORMANCE,
       _("%s      minimum      maximum\n"
         "    depth:                        %10llu | %10llu | %10llu\n"
         "    number of leaves:             %10llu | %10llu | %10llu\n"
         "    number of boxes:              %10llu | %10llu | %10llu\n"
         "    leaves over threshold:        %10llu | %10llu | %10llu\n"
         "    boxes per leaf:               %10llu | %10llu | %10llu\n"
         "    Memory footprint (kb):\n"
         "      final search structure:     %10llu | %10llu | %10llu\n"
         "      temporary search structure: %10llu | %10llu | %10llu\n\n"),
       buf,
       (unsigned long long)(stats->bbox_depth[0] / n),
       (unsigned long long)stats->bbox_depth[1],
       (unsigned long long)stats->bbox_depth[2],
       (unsigned long long)(stats->n_leaves[0] / n),
       (unsigned long long)stats->n_leaves[1],
       (unsigned long long)stats->n_leaves[2],
       (unsigned long long)(stats->n_boxes[0] / n),
       (unsigned long long)stats->n_boxes[1],
       (unsigned long long)stats->n_boxes[2],
       (unsigned long long)(stats->n_th_leaves[0] / n),
       (unsigned long long)stats->n_th_leaves[1],
       (unsigned long long)stats->n_th_leaves[2],
       (unsigned long long)(stats->n_leaf_boxes[0] / n),
       (unsigned long long)stats->n_leaf_boxes[1],
       (unsigned long long)stats->n_leaf_boxes[2],
       (unsigned long long)(stats->box_mem_final[0] / n),
       (unsigned long long)stats->box_mem_final[1],
       (unsigned long long)stats->box_mem_final[2],
       (unsigned long long)(stats->box_mem_required[0] / n),
       (unsigned long long)stats->box_mem_required[1],
       (unsigned long long)stats->box_mem_required[2]);

  }
  else
    cs_log_printf
      (CS_LOG_PERFORMANCE,
       _("    depth:                        %10llu\n"
         "    number of leaves:             %10llu\n"
         "    number of boxes:              %10llu\n"
         "    leaves over threshold:        %10llu\n"
         "    boxes per leaf:               %10llu mean [%llu min, %llu max]\n"
         "    Memory footprint (kb):\n"
         "      final search structure:     %10llu\n"
         "      temporary search structure: %10llu\n\n"),
       (unsigned long long)stats->bbox_depth[0],
       (unsigned long long)stats->n_leaves[0],
       (unsigned long long)stats->n_boxes[0],
       (unsigned long long)stats->n_th_leaves[0],
       (unsigned long long)stats->n_leaf_boxes[0],
       (unsigned long long)stats->n_leaf_boxes[1],
       (unsigned long long)stats->n_leaf_boxes[2],
       (unsigned long long)stats->box_mem_final[0],
       (unsigned long long)stats->box_mem_required[0]);

  cs_log_printf
    (CS_LOG_PERFORMANCE,
     _("  Associated times:\n"
       "    Face bounding boxes tree construction:          %10.3g\n"
       "    Face bounding boxes neighborhood query:         %10.3g\n"),
     stats->t_box_build.wall_nsec*1.e-9,
     stats->t_box_query.wall_nsec*1.e-9);

  cs_log_printf
    (CS_LOG_PERFORMANCE,
     _("    Sorting possible intersections between faces:   %10.3g\n"),
     stats->t_inter_sort.wall_nsec*1e-9);

  cs_log_printf
    (CS_LOG_PERFORMANCE,
     _("    Definition of local joining mesh:               %10.3g\n"),
     stats->t_l_join_mesh.wall_nsec*1e-9);

  cs_log_printf
    (CS_LOG_PERFORMANCE,
     _("    Edge intersections:                             %10.3g\n"),
     stats->t_edge_inter.wall_nsec*1e-9);

  cs_log_printf
    (CS_LOG_PERFORMANCE,
     _("    Creation of new vertices:                       %10.3g\n"),
     stats->t_new_vtx.wall_nsec*1e-9);

  cs_log_printf
    (CS_LOG_PERFORMANCE,
     _("    Merging vertices:                               %10.3g\n"),
     stats->t_merge_vtx.wall_nsec*1e-9);

  cs_log_printf
    (CS_LOG_PERFORMANCE,
     _("    Updating structures with vertex merging:        %10.3g\n"),
     stats->t_u_merge_vtx.wall_nsec*1e-9);

  cs_log_printf
    (CS_LOG_PERFORMANCE,
     _("    Split old faces and reconstruct new faces:      %10.3g\n"),
     stats->t_split_faces.wall_nsec*1e-9);

  cs_log_printf
    (CS_LOG_PERFORMANCE,
     _("\n"
       "  Complete treatment for joining %2d:\n"
       "    wall clock time:                                %10.3g\n"),
     this_join->param.num,
     stats->t_total.wall_nsec*1e-9);

  cs_log_printf_flush(CS_LOG_PERFORMANCE);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Add a cs_join_t structure to the list of pending joinings.
 *
 * parameters:
 *   sel_criteria  <-- boundary face selection criteria
 *   fraction      <-- value of the fraction parameter
 *   plane         <-- value of the plane parameter
 *   verbosity     <-- level of verbosity required
 *   visualization <-- level of visualization required
 *
 * returns:
 *   number (1 to n) associated with new joining
 *---------------------------------------------------------------------------*/

int
cs_join_add(const char  *sel_criteria,
            float        fraction,
            float        plane,
            int          verbosity,
            int          visualization)
{
  /* Allocate and initialize a cs_join_t structure */

  BFT_REALLOC(cs_glob_join_array, cs_glob_n_joinings + 1, cs_join_t *);

  cs_glob_join_array[cs_glob_n_joinings]
    = cs_join_create(cs_glob_n_joinings + 1,
                     sel_criteria,
                     fraction,
                     plane,
                     FVM_PERIODICITY_NULL,
                     NULL,
                     verbosity,
                     visualization,
                     true);

  cs_glob_join_count++; /* Store number of joining (without periodic ones) */
  cs_glob_n_joinings++;

  return cs_glob_n_joinings;
}

/*----------------------------------------------------------------------------
 * Set advanced parameters for the joining algorithm.
 *
 * parameters:
 *   join_num       <-> joining operation number
 *   mtf            <-- merge tolerance coefficient
 *   pmf            <-- pre-merge factor
 *   tcm            <-- tolerance computation mode
 *   icm            <-- intersection computation mode
 *   max_break      <-- max number of equivalences to break (merge step)
 *   max_sub_faces  <-- max. possible number of sub-faces by splitting a face
 *   tml            <-- tree max level
 *   tmb            <-- tree max boxes
 *   tmr            <-- tree max ratio
 *   tmr_distrib    <-- tree max ratio for distribution
 *---------------------------------------------------------------------------*/

void
cs_join_set_advanced_param(int      join_num,
                           double   mtf,
                           double   pmf,
                           int      tcm,
                           int      icm,
                           int      max_break,
                           int      max_sub_faces,
                           int      tml,
                           int      tmb,
                           double   tmr,
                           double   tmr_distrib)
{
  int  i, join_id = -1;
  cs_join_t  *join = NULL;

  /* Search for the joining structure related to join_num */

  for (i = 0; i < cs_glob_n_joinings; i++) {

    join = cs_glob_join_array[i];
    if (join_num == join->param.num) {
      join_id = i;
      break;
    }

  }

  if (join_id < 0)
    bft_error(__FILE__, __LINE__, 0,
              _("  Joining number %d is not defined.\n"), join_num);

  assert(join != NULL);

  _set_advanced_param(join,
                      mtf,
                      pmf,
                      tcm,
                      icm,
                      max_break,
                      max_sub_faces,
                      tml,
                      tmb,
                      tmr,
                      tmr_distrib);
}

/*----------------------------------------------------------------------------
 * Apply all the defined joining operations.
 *
 * parameters:
 *   preprocess <-- true if we are in the preprocessing stage
 *---------------------------------------------------------------------------*/

void
cs_join_all(bool  preprocess)
{
  int  join_id;
  double  full_clock_start, full_clock_end;

  cs_join_type_t  join_type = CS_JOIN_TYPE_NULL;
  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_mesh_builder_t  *mesh_builder = cs_glob_mesh_builder;

  if (cs_glob_n_joinings < 1)
    return;

  /* Sanity checks */

  assert(sizeof(cs_lnum_t) == sizeof(cs_lnum_t));
  assert(sizeof(double) == sizeof(cs_real_t));

  full_clock_start = cs_timer_wtime();

  /* Disable all writers until explicitely enabled for this stage */

  cs_post_disable_writer(0);

  cs_join_post_init();

  /* Loop on each defined joining to deal with */

  for (join_id = 0; join_id < cs_glob_n_joinings; join_id++) {

    if (cs_glob_join_array[join_id] == NULL)
      continue;

    cs_join_t  *this_join = cs_glob_join_array[join_id];

    cs_join_param_t  join_param = this_join->param;

    if (join_param.preprocessing != preprocess)
      continue;

    cs_timer_t t0 = cs_timer_time();  /* Start timer */

    /* Open log file if required */

    if (this_join->log_name != NULL) {
      cs_glob_join_log = fopen(this_join->log_name, "w");
      if (cs_glob_join_log == NULL)
        bft_error(__FILE__, __LINE__, errno,
                  _("Unable to open file: \"%s\" for logging."),
                  this_join->log_name);
    }

    /* Print informations into log file */

    if (mesh->verbosity > 0)
      _print_join_info(mesh, this_join, join_param);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
    _dump_mesh(join_param.num, "InitMesh", mesh);
#endif

    /* Build arrays and structures required for selection;
       will be destroyed after joining and rebuilt for each new join
       operation in order to take into account mesh modification  */

    _select_entities(this_join, mesh);

    /* Now execute the joining operation */

    if (this_join->selection->n_g_faces > 0) {

      cs_lnum_t  n_iwm_vertices;      /* iwm: initial work mesh */
      cs_gnum_t  init_max_vtx_gnum, n_g_new_vertices;

      cs_real_t  *work_face_normal = NULL;
      cs_join_gset_t  *edge_edge_visibility = NULL;
      cs_join_mesh_t  *work_jmesh = NULL, *local_jmesh = NULL;
      cs_join_edges_t  *work_join_edges = NULL;
      cs_join_eset_t  *vtx_eset = NULL;
      cs_join_inter_edges_t  *inter_edges = NULL;

      if (join_param.perio_type != FVM_PERIODICITY_NULL) {

        cs_join_perio_init(this_join, mesh, &mesh_builder);

        if (cs_glob_mesh_builder == NULL)
          cs_glob_mesh_builder = mesh_builder;
      }

      _build_join_structures(this_join,
                             mesh,
                             &local_jmesh,
                             &work_jmesh,
                             &work_join_edges,
                             &work_face_normal,
                             &edge_edge_visibility);

      if (mesh->verbosity > 0 && join_param.verbosity > 2)
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
      if (join_param.perio_type != FVM_PERIODICITY_NULL)
        init_max_vtx_gnum += this_join->selection->n_g_vertices;

      join_type
        = _intersect_edges(this_join,
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

      if (join_type != CS_JOIN_TYPE_NULL)
        _merge_vertices(this_join,
                        n_iwm_vertices,
                        init_max_vtx_gnum,
                        n_g_new_vertices,
                        &vtx_eset,        /* free during this step */
                        &inter_edges,     /* free during this step */
                        &work_jmesh,
                        &work_join_edges,
                        &local_jmesh,
                        mesh);

      else {
        cs_join_eset_destroy(&vtx_eset);
        cs_join_inter_edges_destroy(&inter_edges);
      }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
      _dump_mesh(join_param.num, "MeshAfterMerge", mesh);
#endif

      /*
         Split faces in work_jmesh. Apply modification to the
         local_jmesh. Keep a history between old --> new faces.
         Update cs_mesh_t structure.
      */

      _split_faces(this_join,
                   &join_type,
                   work_join_edges,
                   work_face_normal,
                   &work_jmesh,
                   local_jmesh,
                   mesh,
                   mesh_builder);

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

    if (mesh->verbosity > 0 && join_param.verbosity > 0) {
      bft_printf("\n");
      _print_mesh_info(mesh, _(" After joining"));
      bft_printf("\n");
    }

    /* Free temporary structures */

    cs_join_select_destroy(this_join->param, &(this_join->selection));

   /* Optional synchronization (to be safe) */

#if defined(HAVE_MPI)
    if (cs_glob_n_ranks > 1)
      MPI_Barrier(cs_glob_mpi_comm);
#endif

    /* Close log file if present */

    if (cs_glob_join_log != NULL) {
      if (fclose(cs_glob_join_log) != 0)
        bft_error(__FILE__, __LINE__, errno,
                  _("Error closing log file for joining: %d."),
                  this_join->param.num);
      cs_glob_join_log = NULL;
    }

    /* Timing */

    cs_timer_t t1 = cs_timer_time();
    cs_timer_counter_add_diff(&(this_join->stats.t_total), &t0, &t1);

    if (mesh->verbosity > 0) {
      bft_printf(_("\n"
                   "  Joining %2d completed (%.3g s).\n"),
                 join_param.num,
                 this_join->stats.t_total.wall_nsec*1e-9);
      bft_printf_flush();
    }

    /* Free memory */

    this_join->stats.n_calls += 1;

    if (join_param.preprocessing) {
      _join_performance_log(this_join);
      cs_join_destroy(&this_join);
      cs_glob_join_array[join_id] = NULL;
    }

    /* Force joining visualization level to 0 when the joining is not
       a part of preprocessing to avoid the increase of postprocessing
       meshes during the computation */

    if (!join_param.preprocessing) {
      this_join->param.visualization = 0;
    }

    /* Set mesh modification flag */

    if (join_type != CS_JOIN_TYPE_NULL)
      mesh->modified = 1;

  } /* End of loop on joinings */

  /* Destroy all remaining structures relative to joining operation
     if not needed anymore */

  for (join_id = 0; join_id < cs_glob_n_joinings; join_id++) {
    if (cs_glob_join_array[join_id] != NULL)
      break;
  }
  if (join_id >= cs_glob_n_joinings) {
    BFT_FREE(cs_glob_join_array);
    cs_glob_n_joinings = 0;
  }

  /* Re-enable writers disabled when entering this stage */

  cs_post_enable_writer(0);

  full_clock_end = cs_timer_wtime();

  if (mesh->verbosity > 0) {

    bft_printf(_("\n"
                 "  All joining operations successfully finished:\n"
                 "\n"
                 "    Wall clock time:            %10.3g\n\n"),
               full_clock_end - full_clock_start);
    bft_printf_flush();

    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("\n"
                    "Joining operations time summary:\n"
                    "  wall clock time:            %10.3g\n\n"),
                  full_clock_end - full_clock_start);

    cs_log_separator(CS_LOG_PERFORMANCE);

  }
}

/*----------------------------------------------------------------------------
 * Clear remaining memory for defined joining operations.
 *---------------------------------------------------------------------------*/

void
cs_join_finalize()
{
  bool have_log = false;

  for (int join_id = 0; join_id < cs_glob_n_joinings; join_id++) {
    if (cs_glob_join_array[join_id] != NULL) {
      have_log = true;
      _join_performance_log(cs_glob_join_array[join_id]);
      cs_join_destroy(&(cs_glob_join_array[join_id]));
    }
  }

  BFT_FREE(cs_glob_join_array);
  cs_glob_n_joinings = 0;

  if (have_log) {
    cs_log_printf(CS_LOG_PERFORMANCE, "\n");
    cs_log_separator(CS_LOG_PERFORMANCE);
  }
}

/*----------------------------------------------------------------------------
 * Flag boundary faces that will be selected for joining.
 *
 * parameters:
 *   mesh          <-- pointer to mesh structure
 *   preprocess    <-- true if we are in the preprocessing stage
 *   b_select_flag <-> false for boundary faces not selected for joining,
 *                     true for those selected
 *---------------------------------------------------------------------------*/

void
cs_join_mark_selected_faces(const cs_mesh_t  *mesh,
                            bool              preprocess,
                            bool              b_select_flag[])
{
  for (cs_lnum_t face_id = 0; face_id < mesh->n_b_faces; face_id++)
    b_select_flag[face_id] = false;

  /* Count active joinings */

  int n_joinings = 0;

  for (int join_id = 0; join_id < cs_glob_n_joinings; join_id++) {
    if (cs_glob_join_array[join_id] == NULL)
      continue;
    cs_join_t  *this_join = cs_glob_join_array[join_id];
    cs_join_param_t  join_param = this_join->param;
    if (join_param.preprocessing == preprocess)
      n_joinings++;
  }

  if (n_joinings < 1)
    return;

  /* Prepare selection structures */

  cs_lnum_t *b_face_list;
  BFT_MALLOC(b_face_list, mesh->n_b_faces, cs_lnum_t);

  cs_real_t  *b_face_cog = NULL, *b_face_normal = NULL;

  cs_mesh_quantities_b_faces(mesh, &b_face_cog, &b_face_normal);

  /* Build temporary selection structures */

  const fvm_group_class_set_t *class_defs = mesh->class_defs;
  fvm_group_class_set_t *_class_defs = NULL;

  if (class_defs == NULL) {
    _class_defs = fvm_group_class_set_create();
    class_defs = _class_defs;
  }

  fvm_selector_t *select_b_faces = fvm_selector_create(mesh->dim,
                                                       mesh->n_b_faces,
                                                       class_defs,
                                                       mesh->b_face_family,
                                                       1,
                                                       b_face_cog,
                                                       b_face_normal);

  /* Loop on each defined joining to deal with */

  for (int join_id = 0; join_id < cs_glob_n_joinings; join_id++) {

    if (cs_glob_join_array[join_id] == NULL)
      continue;

    cs_join_t  *this_join = cs_glob_join_array[join_id];

    cs_join_param_t  join_param = this_join->param;

    if (join_param.preprocessing != preprocess)
      continue;

    /* Extract and mark selected boundary faces */

    cs_lnum_t n_b_faces = 0;

    fvm_selector_get_list(select_b_faces,
                          this_join->criteria,
                          1,
                          &n_b_faces,
                          b_face_list);

    for (cs_lnum_t i = 0; i < n_b_faces; i++) {
      cs_lnum_t face_id = b_face_list[i] - 1;
      b_select_flag[face_id] = true;
    }

  } /* End of loop on joinings */

  /* Free arrays and structures needed for selection */

  BFT_FREE(b_face_cog);
  BFT_FREE(b_face_normal);

  select_b_faces = fvm_selector_destroy(select_b_faces);

  if (_class_defs != NULL)
    _class_defs = fvm_group_class_set_destroy(_class_defs);

  BFT_FREE(b_face_list);
}

/*---------------------------------------------------------------------------*/

#if 0 && defined(DEBUG) && !defined(NDEBUG)
cs_debug_glob_mesh_dump("FinalGlobalVertices", mesh);
#endif

END_C_DECLS
