/*============================================================================
 * Management of post-processing for joining operation
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
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *---------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"

#include "fvm_nodal.h"
#include "fvm_nodal_order.h"
#include "fvm_nodal_from_desc.h"
#include "fvm_writer.h"

#include "cs_file.h"
#include "cs_mesh_connect.h"
#include "cs_post.h"
#include "cs_timer_stats.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *---------------------------------------------------------------------------*/

#include "cs_join_post.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Macro and type definitions
 *===========================================================================*/

typedef struct {

  int            writer_num;    /* identifier for the related writer */
  fvm_writer_t  *writer;        /* writer used for post-processing */

} cs_join_post_t;

/* Directory name separator
   (historically, '/' for Unix/Linux, '\' for Windows, ':' for Mac
   but '/' should work for all on modern systems) */

#define DIR_SEPARATOR '/'

/*============================================================================
 * Static global variables
 *===========================================================================*/

static  cs_join_post_t  _cs_join_post_param;

static  bool            _cs_join_post_initialized = false;
static  int             _post_stage_stat_id = -1;

/*============================================================================
 * Private function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Initialize post-processing writer with same format and associated
 * options as default writer, but no time dependency, intended to
 * post elements implied in the joining operations.
 *
 * returns:
 *   id of associated writer (< 0, or 0 in case of failure)
 *----------------------------------------------------------------------------*/

static int
_init_join_writer(void)
{
  int  writer_id = cs_post_get_free_writer_id();

  /* Special case for Catalyst: if matching co-processing script is
     not available, revert to EnSight Gold format */

  int default_format_id
    = fvm_writer_get_format_id(cs_post_get_default_format());

  if (default_format_id == fvm_writer_get_format_id("Catalyst")) {
    if (! cs_file_isreg("error.py"))
      return 0;
  }

  cs_post_define_writer(writer_id,
                        "joining",
                        "postprocessing",
                        fvm_writer_format_name(default_format_id),
                        cs_post_get_default_format_options(),
                        FVM_WRITER_FIXED_MESH,
                        false,
                        false,
                        -1,
                        -1.0);

  return  writer_id;
}

/*----------------------------------------------------------------------------
 * Write a field of "double" on the vertices of the selected mesh.
 * Variable is interlaced.
 *
 * parameters:
 *  mesh      <--  mesh on which we want to write the current field.
 *  varname   <--  name of the field.
 *  dim       <--  dimension of the field to export.
 *  field     <--  variable to write.
 *---------------------------------------------------------------------------*/

static void
_post_vtx_dfield(fvm_nodal_t   *mesh,
                 const char    *varname,
                 int            dim,
                 const double  *field)
{
  fvm_writer_t  *writer = _cs_join_post_param.writer;

  cs_lnum_t  parent_num_shift[2]  = {0, 0};

  const double  *var_ptr[9] = {NULL, NULL, NULL,
                               NULL, NULL, NULL,
                               NULL, NULL, NULL};

  assert(writer != NULL);
  assert(sizeof(double) == 8);

  var_ptr[0] = field;

  fvm_writer_export_field(writer,
                          mesh,
                          varname,
                          FVM_WRITER_PER_NODE,
                          dim,
                          CS_INTERLACE,
                          0,
                          parent_num_shift,
                          CS_DOUBLE,
                          -1,
                          0.,
                          (const void **)var_ptr);
}

/*----------------------------------------------------------------------------
 * Write an integer field on the elements of the selected mesh.
 * Variable is interlaced.
 *
 * parameters:
 *  mesh      <--  mesh on which we want to write the current field.
 *  varname   <--  name of the field.
 *  dim       <--  dimension of the field to export.
 *  field     <--  variable to write.
 *---------------------------------------------------------------------------*/

static void
_post_elt_ifield(fvm_nodal_t  *mesh,
                 const char   *varname,
                 int           dim,
                 const int    *field)
{
  fvm_writer_t  *writer = _cs_join_post_param.writer;

  cs_lnum_t  parent_num_shift[2]  = {0, 0};
  cs_datatype_t  datatype = CS_DATATYPE_NULL;

  const int  *var_ptr[9] = {NULL, NULL, NULL,
                            NULL, NULL, NULL,
                            NULL, NULL, NULL};

  assert(writer != NULL);
  assert(sizeof(cs_lnum_t) == sizeof(int));

  if (sizeof(int) == 4)
    datatype = CS_INT32;
  else if (sizeof(int) == 8)
    datatype = CS_INT64;
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Size of \"int\" is not 4 or 8 bytes.\n"
                " Check the datatype of the field to export.\n"));

  var_ptr[0] = field;

  fvm_writer_export_field(writer,
                          mesh,
                          varname,
                          FVM_WRITER_PER_ELEMENT,
                          dim,
                          CS_INTERLACE,
                          0,
                          parent_num_shift,
                          datatype,
                          -1,
                          0.,
                          (const void **)var_ptr);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Create a writer to output post-processing files for a joining operation.
 *---------------------------------------------------------------------------*/

void
cs_join_post_init(void)
{
  if (_cs_join_post_initialized == true)
    return;

  _post_stage_stat_id = cs_timer_stats_id_by_name("postprocessing_stage");

  int  writer_num = _init_join_writer();

  if (writer_num != 0) {

    _cs_join_post_initialized = true;

    cs_post_activate_writer(writer_num, 1);

    _cs_join_post_param.writer = cs_post_get_writer(writer_num);
    _cs_join_post_param.writer_num = writer_num;

  }
}

/*----------------------------------------------------------------------------
 * Post-treatment of a cs_join_mesh_t structure.
 *
 * parameters:
 *   mesh_name <-- name of the mesh for the post-processing
 *   mesh      <-- pointer to a cs_join_mesh_t structure to post-process
 *---------------------------------------------------------------------------*/

void
cs_join_post_mesh(const char            *mesh_name,
                  const cs_join_mesh_t  *join_mesh)
{
  if (_cs_join_post_initialized == true)
    return;

  int t_top_id = cs_timer_stats_switch(_post_stage_stat_id);

  int  i, j;
  cs_lnum_t  n_vertices;

  const char *name = NULL;
  int  *ifield = NULL;
  double  *dfield = NULL;
  cs_gnum_t  *vertex_gnum = NULL;
  cs_real_t  *vertex_coord = NULL;
  cs_lnum_t  *parent_vtx_num = NULL;
  fvm_nodal_t  *post_mesh = NULL;
  fvm_writer_t  *writer = _cs_join_post_param.writer;

  const int  local_rank = CS_MAX(cs_glob_rank_id, 0);
  const cs_lnum_t  face_list_shift[2] = {0, join_mesh->n_faces};
  const cs_lnum_t  *face_vertex_idx[1] = {join_mesh->face_vtx_idx};
  const cs_lnum_t  *face_vertex_lst[1] = {join_mesh->face_vtx_lst};

  /* Define an fvm_nodal_mesh_t structure from a cs_join_mesh_t structure */

  /* Create an empty fvm_nodal_t structure. */

  if (mesh_name == NULL)
    name = join_mesh->name;
  else
    name = mesh_name;

  post_mesh = fvm_nodal_create(name, 3);

  /* Define fvm_nodal_t structure */

  fvm_nodal_from_desc_add_faces(post_mesh,
                                join_mesh->n_faces,
                                NULL,
                                1,
                                face_list_shift,
                                face_vertex_idx,
                                face_vertex_lst,
                                NULL,
                                NULL);

  /* Define vertex_coord for fvm_nodal_set_shared_vertices() */

  BFT_MALLOC(vertex_coord, 3*join_mesh->n_vertices, cs_real_t);

  for (i = 0; i < join_mesh->n_vertices; i++)
    for (j = 0; j < 3; j++)
      vertex_coord[3*i+j] = (join_mesh->vertices[i]).coord[j];

  fvm_nodal_set_shared_vertices(post_mesh, vertex_coord);

  /* Order faces by increasing global number */

  fvm_nodal_order_faces(post_mesh, join_mesh->face_gnum);
  fvm_nodal_init_io_num(post_mesh, join_mesh->face_gnum, 2);

  /* Order vertices by increasing global number */

  BFT_MALLOC(vertex_gnum, join_mesh->n_vertices, cs_gnum_t);

  for (i = 0; i < join_mesh->n_vertices; i++)
    vertex_gnum[i] = (join_mesh->vertices[i]).gnum;

  fvm_nodal_order_vertices(post_mesh, vertex_gnum);
  fvm_nodal_init_io_num(post_mesh, vertex_gnum, 0);

  /* Write current mesh */

  fvm_writer_export_nodal(writer, post_mesh);

  BFT_FREE(vertex_gnum);
  BFT_FREE(vertex_coord);

  /* Write rank associated to each face */

  BFT_MALLOC(ifield, join_mesh->n_faces, int);

  for (i = 0; i < join_mesh->n_faces; i++)
    ifield[i] = local_rank;

  _post_elt_ifield(post_mesh, _("Rank"), 1, ifield);

  BFT_FREE(ifield);

  /* Write vertex tolerance */

  n_vertices = fvm_nodal_get_n_entities(post_mesh, 0);

  BFT_MALLOC(parent_vtx_num, n_vertices, cs_lnum_t);
  BFT_MALLOC(dfield, n_vertices, double);

  fvm_nodal_get_parent_num(post_mesh, 0, parent_vtx_num);

  for (i = 0; i < n_vertices; i++) {

    cs_join_vertex_t  data = join_mesh->vertices[parent_vtx_num[i]-1];

    dfield[i] = data.tolerance;
  }

  _post_vtx_dfield(post_mesh, _("VtxTolerance"), 1, dfield);

  BFT_FREE(parent_vtx_num);
  BFT_FREE(dfield);

  post_mesh = fvm_nodal_destroy(post_mesh);

  cs_timer_stats_switch(t_top_id);
}

/*----------------------------------------------------------------------------
 * Post-process a subset of faces of a cs_join_mesh_t structure.
 *
 * parameters:
 *   mesh_name        <-- name of the sub-set mesh
 *   mesh             <-- pointer to the parent cs_join_mesh_t structure
 *   n_selected_faces <-- number of selected faces (size of the sub-set)
 *   selected_faces   <-- list of local number in parent mesh
 *---------------------------------------------------------------------------*/

void
cs_join_post_faces_subset(const char            *mesh_name,
                          const cs_join_mesh_t  *parent_mesh,
                          cs_lnum_t              n_select_faces,
                          const cs_lnum_t        selected_faces[])
{
  if (_cs_join_post_initialized == true)
    return;

  int t_top_id = cs_timer_stats_switch(_post_stage_stat_id);

  cs_join_mesh_t  *subset_mesh = NULL;

  assert(parent_mesh != NULL);

  subset_mesh = cs_join_mesh_create_from_subset(mesh_name,
                                                n_select_faces,
                                                selected_faces,
                                                parent_mesh);

  cs_join_post_mesh(subset_mesh->name, subset_mesh);

  cs_join_mesh_destroy(&subset_mesh);

  cs_timer_stats_switch(t_top_id);
}

/*----------------------------------------------------------------------------
 * Post-process mesh after the update following the merge operation.
 *
 * parameters:
 *   join_param  <-- set of parameters for the joining operation
 *   join_select <-- list of participating entities in the joining operation
 *---------------------------------------------------------------------------*/

void
cs_join_post_after_merge(cs_join_param_t          join_param,
                         const cs_join_select_t  *join_select)
{
  if (_cs_join_post_initialized == true)
    return;

  int t_top_id = cs_timer_stats_switch(_post_stage_stat_id);

  int  adj_mesh_id, sel_mesh_id;

  int  writer_ids[] = {_cs_join_post_param.writer_num};
  char  *mesh_name = NULL;
  fvm_nodal_t *adj_mesh = NULL, *sel_mesh = NULL;

  adj_mesh_id = cs_post_get_free_mesh_id();

  BFT_MALLOC(mesh_name, strlen("AdjacentJoinFaces_j") + 2 + 1, char);
  sprintf(mesh_name,"%s%02d", "AdjacentJoinFaces_j", join_param.num);

  adj_mesh = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                            mesh_name,
                                            false, /* include families */
                                            join_select->n_i_adj_faces,
                                            join_select->n_b_adj_faces,
                                            join_select->i_adj_faces,
                                            join_select->b_adj_faces);

  cs_post_define_existing_mesh(adj_mesh_id,
                               adj_mesh,
                               0,    /* dim_shift */
                               true, /* transfer ownership */
                               false,
                               1,
                               writer_ids);

  sel_mesh_id = cs_post_get_free_mesh_id();

  BFT_REALLOC(mesh_name, strlen("JoinFacesAfterMerge_j") + 2 + 1, char);
  sprintf(mesh_name,"%s%02d", "JoinFacesAfterMerge_j", join_param.num);

  sel_mesh = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                            mesh_name,
                                            false, /* include families */
                                            0,
                                            join_select->n_faces,
                                            NULL,
                                            join_select->faces);

  cs_post_define_existing_mesh(sel_mesh_id,
                               sel_mesh,
                               0,    /* dim_shift */
                               true, /* transfer ownership */
                               false,
                               1,
                               writer_ids);

  /* Post */

  cs_post_activate_writer(_cs_join_post_param.writer_num, 1);
  cs_post_write_meshes(NULL);

  cs_post_free_mesh(sel_mesh_id);
  cs_post_free_mesh(adj_mesh_id);

  BFT_FREE(mesh_name);

  cs_timer_stats_switch(t_top_id);
}

/*----------------------------------------------------------------------------
 * Post-process mesh after the update following the split operation.
 *
 * parameters:
 *   n_old_i_faces   <-- initial number of interior faces
 *   n_old_b_faces   <-- initial number of border faces
 *   n_g_new_b_faces <-- global number of new border faces
 *   n_select_faces  <-- number of selected faces
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   join_param      <-- set of parameters for the joining operation
 *---------------------------------------------------------------------------*/

void
cs_join_post_after_split(cs_lnum_t         n_old_i_faces,
                         cs_lnum_t         n_old_b_faces,
                         cs_gnum_t         n_g_new_b_faces,
                         cs_lnum_t         n_select_faces,
                         const cs_mesh_t  *mesh,
                         cs_join_param_t   join_param)
{
  if (join_param.visualization < 1 || _cs_join_post_initialized == false)
    return;

  int t_top_id = cs_timer_stats_switch(_post_stage_stat_id);

  cs_lnum_t  i, j;

  int  writer_ids[] = {_cs_join_post_param.writer_num};
  char  *mesh_name = NULL;
  cs_lnum_t  *post_i_faces = NULL, *post_b_faces = NULL;
  fvm_nodal_t  *post_i_mesh = NULL;
  int  post_i_mesh_id = cs_post_get_free_mesh_id();
  int  post_b_mesh_id = 0;

  const int  n_new_i_faces = mesh->n_i_faces - n_old_i_faces;
  const int  n_new_b_faces = mesh->n_b_faces - n_old_b_faces + n_select_faces;

  /* Define list of faces to post-treat */

  BFT_MALLOC(post_i_faces, n_new_i_faces, cs_lnum_t);
  BFT_MALLOC(post_b_faces, n_new_b_faces, cs_lnum_t);

  for (i = n_old_i_faces, j = 0; i < mesh->n_i_faces; i++, j++)
    post_i_faces[j] = i + 1;

  for (i = n_old_b_faces-n_select_faces, j = 0; i < mesh->n_b_faces; i++, j++)
    post_b_faces[j] = i + 1;

  BFT_MALLOC(mesh_name, strlen("InteriorJoinedFaces_j") + 2 + 1, char);
  sprintf(mesh_name,"%s%02d", "InteriorJoinedFaces_j", join_param.num);

  post_i_mesh = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                               mesh_name,
                                               false, /* include families */
                                               n_new_i_faces,
                                               0,
                                               post_i_faces,
                                               NULL);

  cs_post_define_existing_mesh(post_i_mesh_id,
                               post_i_mesh,
                               0,    /* dim_shift */
                               true, /* transfer ownership */
                               false,
                               1,
                               writer_ids);

  if (join_param.visualization > 1 && n_g_new_b_faces > 0) {

    fvm_nodal_t  *post_b_mesh = NULL;
    post_b_mesh_id = cs_post_get_free_mesh_id();

    BFT_REALLOC(mesh_name, strlen("BoundaryJoinedFaces_j") + 2 + 1, char);
    sprintf(mesh_name,"%s%02d", "BoundaryJoinedFaces_j", join_param.num);

    post_b_mesh = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                                 mesh_name,
                                                 false, /* include families */
                                                 0,
                                                 n_new_b_faces,
                                                 NULL,
                                                 post_b_faces);

    cs_post_define_existing_mesh(post_b_mesh_id,
                                 post_b_mesh,
                                 0,    /* dim_shift */
                                 true, /* transfer ownership */
                                 false,
                                 1,
                                 writer_ids);

  }

  /* Post */

  cs_post_activate_writer(_cs_join_post_param.writer_num, 1);
  cs_post_write_meshes(NULL);

  if (post_b_mesh_id != 0)
    cs_post_free_mesh(post_b_mesh_id);
  cs_post_free_mesh(post_i_mesh_id);

  BFT_FREE(post_i_faces);
  BFT_FREE(post_b_faces);
  BFT_FREE(mesh_name);

  cs_timer_stats_switch(t_top_id);
}

/*----------------------------------------------------------------------------
 * Post-process mesh after the update following the split operation.
 *
 * parameters:
 *   n_i_clean_faces <-- number of interior faces cleaned
 *   i_clean_faces   <-> list of interior face numbers (ordered on exit)
 *   n_b_clean_faces <-- number of border faces cleaned
 *   b_clean_faces   <-> list of border face numbers (ordered on exit)
 *   param           <-- set of parameters for the joining operation
 *---------------------------------------------------------------------------*/

void
cs_join_post_cleaned_faces(cs_lnum_t        n_i_clean_faces,
                           cs_lnum_t        i_clean_faces[],
                           cs_lnum_t        n_b_clean_faces,
                           cs_lnum_t        b_clean_faces[],
                           cs_join_param_t  param)
{
  if (_cs_join_post_initialized == false)
    return;

  int t_top_id = cs_timer_stats_switch(_post_stage_stat_id);

  int  writer_ids[] = {_cs_join_post_param.writer_num};
  int  post_mesh_id = cs_post_get_free_mesh_id();
  char  *name = NULL;
  fvm_nodal_t *export_mesh = NULL;

  BFT_MALLOC(name, strlen("CleanFaces_j") + 2 + 1, char);
  sprintf(name,"%s%02d", "CleanFaces_j", param.num);

  export_mesh = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                               name,
                                               false, /* include families */
                                               n_i_clean_faces,
                                               n_b_clean_faces,
                                               i_clean_faces,
                                               b_clean_faces);

  cs_post_define_existing_mesh(post_mesh_id,
                               export_mesh,
                               0,    /* dim_shift */
                               true, /* transfer ownership */
                               false,
                               1,
                               writer_ids);

  /* Output post-processing data */

  cs_post_activate_writer(_cs_join_post_param.writer_num, 1);
  cs_post_write_meshes(NULL);

  cs_post_free_mesh(post_mesh_id);

  BFT_FREE(name);

  cs_timer_stats_switch(t_top_id);
}

/*----------------------------------------------------------------------------
 * Output processor-specific post-processing data for a cs_join_mesh_t
 * structure according to the visualization level.
 *
 * parameters:
 *   basename <-- generic name for the mesh to post
 *   mesh     <-- fvm_join_mesh_t structure to post-process
 *   param    <-- fvm_join_param_t structure
 *---------------------------------------------------------------------------*/

void
cs_join_post_dump_mesh(const char            *basename,
                       const cs_join_mesh_t  *mesh,
                       cs_join_param_t        param)
{
  int  rank, len;

  cs_join_mesh_t  *tmp = NULL;
  char  *fullname = NULL;

  const  int  n_ranks = cs_glob_n_ranks;
  const  int  rank_id = CS_MAX(cs_glob_rank_id, 0);

  /* Define a specific name for the output */

  len = strlen("log/JoinDBG_.dat") + strlen(basename) + 4 + 2 + 1;
  BFT_MALLOC(fullname, len, char);
  sprintf(fullname, "log%cJoin%02dDBG_%s%04d.dat", DIR_SEPARATOR,
          param.num, basename, rank_id);

#if 0 && defined(DEBUG) && !defined(NDEBUG) /* Dump mesh structure */
  if (param.verbosity > 3) {
    FILE  *dbg_file = NULL;
    dbg_file = fopen(fullname, "w");
    cs_join_mesh_dump_file(dbg_file, mesh);
    fflush(dbg_file);
    fclose(dbg_file);
  }
#endif

  if (_cs_join_post_initialized == true && param.visualization > 3) {

    if (n_ranks == 1)
      cs_join_post_mesh(fullname, mesh);

    else { /* Parallel */

      for (rank = 0; rank < n_ranks; rank++) {

        char *mesh_name = NULL;

        BFT_MALLOC(mesh_name, strlen(basename) + 2 + 2 + 5 + 1, char);
        sprintf(mesh_name,"%s%02d%s%05d", basename, param.num, "_n", rank);

        if (rank_id == rank)
          cs_join_post_mesh(mesh_name, mesh);

        else { /* Pieces empty on other ranks */
          tmp = cs_join_mesh_create(mesh_name);
          cs_join_post_mesh(mesh_name, tmp);
          cs_join_mesh_destroy(&tmp);
        }

        BFT_FREE(mesh_name);

      } /* End of loop on ranks */
    } /* End of parallel treatment */
  }

  BFT_FREE(fullname);

#if defined(HAVE_MPI)
  if (n_ranks > 1)
    MPI_Barrier(cs_glob_mpi_comm);
#endif
}

/*---------------------------------------------------------------------------*/

END_C_DECLS
