/*============================================================================
 * Management of conforming and non-conforming joining in case of periodicity
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
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *---------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_periodicity.h"

#include "cs_mesh.h"
#include "cs_search.h"

#include "cs_join_mesh.h"
#include "cs_join_post.h"
#include "cs_join_set.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *---------------------------------------------------------------------------*/

#include "cs_join_perio.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Structure Definitions
 *===========================================================================*/

/*============================================================================
 * Global variables
 *===========================================================================*/

/*============================================================================
 * Static global variables
 *===========================================================================*/

/*============================================================================
 * Private function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Add a new periodic cs_join_t structure.
 *
 * parameters:
 *   perio_type    <-- periodicity number
 *   perio_matrix  <-- periodicity transformation matrix
 *   sel_criteria  <-- boundary face selection criteria
 *   fraction      <-- value of the fraction parameter
 *   plane         <-- value of the plane parameter
 *   verbosity     <-- level of verbosity required
 *   visualization <-- level of visualization required
 *
 * returns:
 *   the global joining number associated with the periodic joining
 *---------------------------------------------------------------------------*/

static int
_add_perio_join(fvm_periodicity_type_t  perio_type,
                double                  matrix[3][4],
                const char             *criteria,
                float                   fraction,
                float                   plane,
                int                     verbosity,
                int                     visualization)
{
   /* Allocate and initialize a cs_join_t structure */

  BFT_REALLOC(cs_glob_join_array, cs_glob_n_joinings + 1, cs_join_t *);

  cs_glob_join_array[cs_glob_n_joinings]
    = cs_join_create(cs_glob_n_joinings + 1,
                     criteria,
                     fraction,
                     plane,
                     perio_type,
                     matrix,
                     verbosity,
                     visualization,
                     true);

  cs_glob_n_joinings++;

  return cs_glob_n_joinings;
}

/*----------------------------------------------------------------------------
 * Update work_mesh by redistributing local join mesh.
 *
 * parameters:
 *   param             <--  set of user-defined parameter
 *   gnum_rank_index   <--  index on ranks for the old global face numbering
 *   local_mesh        <--  mesh on local selected faces to be joined
 *   p_work_mesh       <->  distributed mesh on faces to join
 *---------------------------------------------------------------------------*/

static void
_redistribute_mesh(cs_join_param_t         param,
                   const cs_gnum_t         gnum_rank_index[],
                   const cs_join_mesh_t   *local_mesh,
                   cs_join_mesh_t        **p_work_mesh)
{
  cs_join_mesh_t  *work_mesh = *p_work_mesh;

  const int  n_ranks = cs_glob_n_ranks;
  const int  local_rank = CS_MAX(cs_glob_rank_id, 0);

  /* sanity checks */

  assert(local_mesh != NULL);
  assert(work_mesh != NULL);

  if (n_ranks == 1)
    cs_join_mesh_copy(&work_mesh, local_mesh);

#if defined(HAVE_MPI)
  if (n_ranks > 1) { /* Parallel mode */

    int  i, n_work_faces;

    char  *mesh_name = NULL;
    cs_gnum_t  *work_faces = NULL;

    n_work_faces = work_mesh->n_faces;
    BFT_MALLOC(work_faces, n_work_faces, cs_gnum_t);

    for (i = 0; i < n_work_faces; i++)
      work_faces[i] = work_mesh->face_gnum[i];

    /* Replace current work mesh */

    cs_join_mesh_destroy(&work_mesh);

    BFT_MALLOC(mesh_name, strlen("WorkMesh_j_n") + 2 + 5 + 1, char);
    sprintf(mesh_name,"%s%02d%s%05d",
            "WorkMesh_j", param.num, "_n", local_rank);

    work_mesh = cs_join_mesh_create_from_glob_sel(mesh_name,
                                                  n_work_faces,
                                                  work_faces,
                                                  gnum_rank_index,
                                                  local_mesh);

    BFT_FREE(mesh_name);
    BFT_FREE(work_faces);

  }
#endif

  /* Return pointers */

  *p_work_mesh = work_mesh;

}

/*----------------------------------------------------------------------------
 * Delete interior faces temporary added for periodicity operation.
 * Only done if n_ranks > 1 because in serial, these periodic faces are
 * merged with initially present faces
 *
 * parameters:
 *   param <-- set of parameters for the joining operation
 *   mesh  <-> cs_mesh_t structure to update
 *---------------------------------------------------------------------------*/

static void
_perio_face_clean(cs_join_param_t      param,
                  cs_mesh_t           *mesh)
{
  int  i, j, k, shift;

  int  n_ii_faces = mesh->n_i_faces;
  int  n_fi_faces = 0;
  cs_lnum_t  *new_f2v_idx = NULL;
  int  *tag = NULL;

  assert(cs_glob_n_ranks > 1);

  BFT_MALLOC(tag, n_ii_faces, int);

  for (i = 0; i < n_ii_faces; i++) {

    if (mesh->i_face_cells[i][0] == -1 && mesh->i_face_cells[i][1] == -1)
      tag[i] = -1;
    else {
      mesh->i_face_cells[n_fi_faces][0] = mesh->i_face_cells[i][0];
      mesh->i_face_cells[n_fi_faces][1] = mesh->i_face_cells[i][1];
      n_fi_faces++;
      tag[i] = n_fi_faces;
    }

  }

  if (param.verbosity > 3)
    fprintf(cs_glob_join_log,
            "\n  Delete %d interior periodic faces locally\n",
            n_ii_faces - n_fi_faces);

  mesh->n_i_faces = n_fi_faces;
  BFT_REALLOC(mesh->i_face_cells, mesh->n_i_faces, cs_lnum_2_t);
  BFT_MALLOC(new_f2v_idx, n_fi_faces + 1, cs_lnum_t);

  n_fi_faces = 0;
  for (i = 0; i < n_ii_faces; i++) {
    if (tag[i] > 0) {
      mesh->global_i_face_num[n_fi_faces] = mesh->global_i_face_num[i];
      mesh->i_face_family[n_fi_faces] = mesh->i_face_family[i];
      new_f2v_idx[n_fi_faces + 1] =  mesh->i_face_vtx_idx[i+1]
                                   - mesh->i_face_vtx_idx[i];
      n_fi_faces++;
    }
  }

  BFT_REALLOC(mesh->global_i_face_num, mesh->n_i_faces, cs_gnum_t);
  BFT_REALLOC(mesh->i_face_family, mesh->n_i_faces, cs_lnum_t);

  /* Update interior face connectivity */

  new_f2v_idx[0] = 0;
  for (i = 0; i < n_fi_faces; i++)
    new_f2v_idx[i+1] += new_f2v_idx[i];

  n_fi_faces = 0;
  for (i = 0; i < n_ii_faces; i++) {
    if (tag[i] > 0) {
      shift = new_f2v_idx[n_fi_faces];
      for (k = 0, j = mesh->i_face_vtx_idx[i];
           j < mesh->i_face_vtx_idx[i+1]; j++, k++)
        mesh->i_face_vtx_lst[shift+k] = mesh->i_face_vtx_lst[j];
      n_fi_faces++;
    }
  }

  BFT_REALLOC(mesh->i_face_vtx_lst, new_f2v_idx[n_fi_faces], cs_lnum_t);
  BFT_FREE(mesh->i_face_vtx_idx);

  mesh->i_face_vtx_idx = new_f2v_idx;
  mesh->i_face_vtx_connect_size = new_f2v_idx[n_fi_faces];

  /* There is no need to define a new glbal interior face numbering
     because the excluded faces are always defined on an another rank */

  BFT_FREE(tag);
}

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Check if periodic joining operations are queued
 *
 * Fortran Interface:
 *
 * SUBROUTINE TSTJPE
 * *****************
 *
 * INTEGER        iperio    : <-> : do we have periodicity ?
 * INTEGER        iperot    : <-> : do we have periodicity of rotation ?
 *----------------------------------------------------------------------------*/

void CS_PROCF(tstjpe, tstjpe)
(
 cs_int_t    *iperio,
 cs_int_t    *iperot
)
{
  int i;

  for (i = 0; i < cs_glob_n_joinings; i++) {
    cs_join_param_t param = (cs_glob_join_array[i])->param;
    if (param.perio_type > FVM_PERIODICITY_NULL)
      *iperio = 1;
    if (param.perio_type > FVM_PERIODICITY_TRANSLATION)
      *iperot = 1;
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Define a translational periodicity
 *
 * parameters:
 *   sel_criteria  <-- boundary face selection criteria
 *   fraction      <-- value of the fraction parameter
 *   plane         <-- value of the plane parameter
 *   verbosity     <-- level of verbosity required
 *   visualization <-- level of visualization required
 *   trans         <-- translation vector
 *
 * returns:
 *   joining number (1 to n) associated with new periodicity
 *----------------------------------------------------------------------------*/

int
cs_join_perio_add_translation(const char    *sel_criteria,
                              double         fraction,
                              double         plane,
                              int            verbosity,
                              int            visualization,
                              const double   trans[3])
{
  double  matrix[3][4];
  fvm_periodicity_t *tmp_perio = fvm_periodicity_create(0.001);

  int join_num = 0;

  assert((trans[0]*trans[0] + trans[1]*trans[1] + trans[2]*trans[2]) > 0.0);

  /* Use temporary periodicity structure to convert parameters to matrix */

  fvm_periodicity_add_translation(tmp_perio, 1, trans);

  fvm_periodicity_get_matrix(tmp_perio, 0, matrix);

  join_num = _add_perio_join(FVM_PERIODICITY_TRANSLATION,
                             matrix,
                             sel_criteria,
                             fraction,
                             plane,
                             verbosity,
                             visualization);

  tmp_perio = fvm_periodicity_destroy(tmp_perio);

  return join_num;
}

/*----------------------------------------------------------------------------
 * Define a rotational periodicity
 *
 * parameters:
 *   perio_num     <-- number related to the periodicity
 *   sel_criteria  <-- boundary face selection criteria
 *   fraction      <-- value of the fraction parameter
 *   plane         <-- value of the plane parameter
 *   verbosity     <-- level of verbosity required
 *   visualization <-- level of visualization required
 *   theta         <-- rotation angle (in degrees)
 *   axis          <-- axis vector
 *   invariant     <-- invariant point coordinates
 *
 * returns:
 *   joining number (1 to n) associated with new periodicity
 *----------------------------------------------------------------------------*/

int
cs_join_perio_add_rotation(const char    *sel_criteria,
                           double         fraction,
                           double         plane,
                           int            verbosity,
                           int            visualization,
                           double         theta,
                           const double   axis[3],
                           const double   invariant[3])
{
  double  matrix[3][4];
  fvm_periodicity_t *tmp_perio = fvm_periodicity_create(0.001);

  int join_num = 0;

  assert(theta*theta > 0.0);
  assert((axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]) > 0.0);

  /* Use temporary periodicity structure to convert parameters to matrix */

  fvm_periodicity_add_rotation(tmp_perio,
                               1,
                               theta,
                               axis,
                               invariant);

  fvm_periodicity_get_matrix(tmp_perio, 0, matrix);

  join_num = _add_perio_join(FVM_PERIODICITY_ROTATION,
                             matrix,
                             sel_criteria,
                             fraction,
                             plane,
                             verbosity,
                             visualization);

  tmp_perio = fvm_periodicity_destroy(tmp_perio);

  /* Add a tag to indicate the use of rotation */
  cs_glob_mesh->have_rotation_perio = 1;

  return join_num;
}

/*----------------------------------------------------------------------------
 * Define a periodicity using a matrix
 *
 * parameters:
 *   perio_num     <-- number related to the periodicity
 *   sel_criteria  <-- boundary face selection criteria
 *   fraction      <-- value of the fraction parameter
 *   plane         <-- value of the plane parameter
 *   verbosity     <-- level of verbosity required
 *   visualization <-- level of visualization required
 *   matrix        <-- transformation matrix
 *
 * returns:
 *   joining number (1 to n) associated with new periodicity
 *----------------------------------------------------------------------------*/

int
cs_join_perio_add_mixed(const char    *sel_criteria,
                        double         fraction,
                        double         plane,
                        int            verbosity,
                        int            visualization,
                        double         matrix[3][4])
{
  int join_num = 0;

  join_num = _add_perio_join(FVM_PERIODICITY_MIXED,
                             matrix,
                             sel_criteria,
                             fraction,
                             plane,
                             verbosity,
                             visualization);

  /* Add a tag to indicate the use of rotation */
  cs_glob_mesh->have_rotation_perio = 1;

  return join_num;
}

/*----------------------------------------------------------------------------
 * Add periodicity information to mesh and create or update mesh builder
 * for a new periodic joining.
 *
 * parameters:
 *   this_join <-- high level join structure
 *   mesh      <-> pointer to a cs_mesh_t structure
 *   builder   <-> pointer to a cs_mesh_builder_t structure pointer
 *---------------------------------------------------------------------------*/

void
cs_join_perio_init(cs_join_t           *this_join,
                   cs_mesh_t           *mesh,
                   cs_mesh_builder_t  **builder)
{
  cs_join_param_t  param = this_join->param;
  fvm_periodicity_t  *periodicity = NULL;
  cs_mesh_builder_t  *_builder;

  int  perio_id = - 1;

  assert(mesh != NULL);

  /* Update mesh periodicity information */

  if (mesh->periodicity == NULL)
    mesh->periodicity = fvm_periodicity_create(0.001);

  periodicity = mesh->periodicity;

  mesh->n_init_perio += 1;
  if (param.perio_type > FVM_PERIODICITY_TRANSLATION)
    mesh->have_rotation_perio = 1;

  perio_id = fvm_periodicity_get_n_transforms(periodicity)/2;

  fvm_periodicity_add_by_matrix(periodicity,
                                perio_id + 1,
                                param.perio_type,
                                param.perio_matrix);

  assert(mesh->n_init_perio == perio_id + 1);

  /* Update builder */

  if (*builder != NULL)
    _builder = *builder;

  else {
    _builder = cs_mesh_builder_create();
    *builder = _builder;
  }

  _builder->n_perio += 1;

  assert(mesh->n_init_perio == _builder->n_perio);

  BFT_REALLOC(_builder->n_per_face_couples, mesh->n_init_perio, cs_lnum_t);
  BFT_REALLOC(_builder->per_face_couples, mesh->n_init_perio, cs_gnum_t *);

  _builder->n_per_face_couples[mesh->n_init_perio - 1] = 0;
  _builder->per_face_couples[mesh->n_init_perio - 1] = NULL;
}

/*----------------------------------------------------------------------------
 * Duplicate and apply transformation to the selected faces and also to
 * their related vertices. Modify compact_face_gnum to take into account
 * new periodic faces and create a periodic vertex couple list.
 *
 * parameters:
 *   this_join <-- high level join structure
 *   jmesh     <-> local join mesh struct. to duplicate and transform
 *   mesh      <-- pointer to a cs_mesh_t structure
 *---------------------------------------------------------------------------*/

void
cs_join_perio_apply(cs_join_t          *this_join,
                    cs_join_mesh_t     *jmesh,
                    const cs_mesh_t    *mesh)
{
  cs_lnum_t  i, j, k, shift;
  cs_real_t  matrix[3][4], xyz[4];

  cs_join_param_t  param = this_join->param;
  cs_join_select_t  *select = this_join->selection;
  fvm_periodicity_t  *periodicity = mesh->periodicity;

  int  perio_id = - 1;
  const int  n_ranks = cs_glob_n_ranks;
  const int  n_init_vertices = jmesh->n_vertices;
  const int  n_init_faces = jmesh->n_faces;

  /* Get mesh periodicity information */

  perio_id = fvm_periodicity_get_n_transforms(mesh->periodicity)/2 - 1;
  assert(perio_id > -1);

  fvm_periodicity_get_matrix(mesh->periodicity, 2*perio_id+1, matrix);

  /* Retrieve related transformation */

  fvm_periodicity_get_matrix(periodicity, 2*perio_id, matrix);

  /* Duplicate and transform vertices */

  jmesh->n_vertices *= 2;
  jmesh->n_g_vertices *= 2;

  BFT_REALLOC(jmesh->vertices, jmesh->n_vertices, cs_join_vertex_t);

  shift = n_init_vertices;
  for (i = 0; i < n_init_vertices; i++) {

    /* Copy tolerance, coord, state and gnum to initialize new_vtx */

    cs_join_vertex_t  new_vtx = jmesh->vertices[i];

    for (j = 0; j < 3; j++) {
      xyz[j] = new_vtx.coord[j];
      new_vtx.coord[j] = 0.0;
    }
    xyz[3] = 1;

    for (j = 0; j < 3; j++)
      for (k = 0; k < 4; k++)
        new_vtx.coord[j] += matrix[j][k]*xyz[k];

    new_vtx.state = CS_JOIN_STATE_PERIO;
    jmesh->vertices[shift++] = new_vtx;

  }

  /* Add a periodic vertex couple list */

  select->n_couples = n_init_vertices;
  BFT_MALLOC(select->per_v_couples, 2*n_init_vertices, cs_gnum_t);

  if (n_ranks > 1) { /* Global numbering update */

    cs_gnum_t  *gnum = NULL;
    fvm_io_num_t  *io_num = NULL;
    const cs_gnum_t  *io_gnum = NULL;

    BFT_MALLOC(gnum, n_init_vertices, cs_gnum_t);

    for (i = 0, shift = n_init_vertices; i < n_init_vertices; i++, shift++)
      gnum[i] = jmesh->vertices[shift].gnum;

    io_num = fvm_io_num_create(NULL, gnum, n_init_vertices, 0);
    io_gnum = fvm_io_num_get_global_num(io_num);

    for (i = 0, shift = n_init_vertices; i < n_init_vertices; i++, shift++) {
      jmesh->vertices[shift].gnum = io_gnum[i] + mesh->n_g_vertices;
      select->per_v_couples[2*i] = jmesh->vertices[i].gnum;
      select->per_v_couples[2*i+1] = jmesh->vertices[shift].gnum;
    }

    fvm_io_num_destroy(io_num);
    BFT_FREE(gnum);

  }
  else { /* Serial run */

    for (i = 0, shift = n_init_vertices; i < n_init_vertices; i++, shift++) {
      jmesh->vertices[shift].gnum = i + 1 + mesh->n_g_vertices;
      select->per_v_couples[2*i] = jmesh->vertices[i].gnum;
      select->per_v_couples[2*i+1] = jmesh->vertices[shift].gnum;
    }

  }

  /* Duplicate and transform faces */

  jmesh->n_faces *= 2;
  jmesh->n_g_faces *= 2;

  BFT_REALLOC(jmesh->face_vtx_idx, jmesh->n_faces + 1, cs_lnum_t);
  BFT_REALLOC(jmesh->face_gnum, jmesh->n_faces, cs_gnum_t);
  BFT_REALLOC(jmesh->face_vtx_lst,
              2*(jmesh->face_vtx_idx[n_init_faces]), cs_lnum_t);

  for (i = 0; i < n_init_faces; i++) {

    int  pfid = n_init_faces + i;
    int  s = jmesh->face_vtx_idx[i];
    int  e = jmesh->face_vtx_idx[i+1];
    int  ps = jmesh->face_vtx_idx[pfid];
    int  pe = jmesh->face_vtx_idx[pfid] + e - s;
    cs_gnum_t  new_gnum = 2*jmesh->face_gnum[i];

    jmesh->face_gnum[i] = new_gnum - 1;
    jmesh->face_gnum[pfid] = new_gnum;

    for (j = s, shift = ps; j < e; j++, shift++)
      jmesh->face_vtx_lst[shift] = jmesh->face_vtx_lst[j] + n_init_vertices;
    jmesh->face_vtx_idx[pfid+1] =  pe;

  }

  /* Modify cs_join_select_t structure */

  for (i = 0; i < n_ranks + 1; i++)
    select->compact_rank_index[i] *= 2;

  for (i = 0; i < select->n_faces; i++)
    select->compact_face_gnum[i] = 2*select->compact_face_gnum[i] - 1;

  /* Order faces in mesh struct. by increasing global face number.
     We have to be order in this way in order to keep an exact
     redistribution (work_mesh) */

  cs_join_mesh_face_order(jmesh);

  if (param.verbosity > 2)
    fprintf(cs_glob_join_log,
            "  Apply periodicity to the local join mesh structure\n"
            "  New number of faces to treat locally: %8d\n",
            jmesh->n_faces);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (cs_glob_join_log != NULL) {
    fprintf(cs_glob_join_log, "\n Periodic vertex couples:\n");
    for (i = 0; i < n_init_vertices; i++)
      fprintf(cs_glob_join_log, " %6d  (%9u, %9u)\n", i+1,
              select->per_v_couples[2*i],  select->per_v_couples[2*i+1]);
    fflush(cs_glob_join_log);
  }
#endif

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (cs_glob_join_log != NULL) {
    fprintf(cs_glob_join_log, "\n  Selected faces for the joining operation:\n");
    for (i = 0; i < select->n_faces; i++)
      fprintf(cs_glob_join_log, " %9d | %9d | %10llu\n",
              i, select->faces[i],
              (unsigned long long)select->compact_face_gnum[i]);
    fprintf(cs_glob_join_log, "\n");
  }
#endif

}

/*----------------------------------------------------------------------------
 * Duplicate and apply transformation to the selected faces and also to
 * their related vertices. Modify compact_face_gnum to take into account
 * new periodic faces and create a periodic vertex couple list.
 *
 * parameters:
 *   this_join          <-- pointer to a high level join structure
 *   jmesh              <-> local join mesh struct. to duplicate and transform
 *   mesh               <-- pointer to a cs_mesh_t structure
 *   p_work_jmesh       <-> distributed join mesh struct. on which operations
 *                          take place
 *   p_work_edges       <-> join edges struct. related to work_jmesh
 *   init_max_vtx_gnum  <-- initial max. global numbering for vertices
 *   n_g_new_vertices   <-- global number of vertices created during the
 *                          intersection of edges
 *---------------------------------------------------------------------------*/

void
cs_join_perio_merge_back(cs_join_t          *this_join,
                         cs_join_mesh_t     *jmesh,
                         const cs_mesh_t    *mesh,
                         cs_join_mesh_t    **p_work_jmesh,
                         cs_join_edges_t   **p_work_edges,
                         cs_gnum_t           init_max_vtx_gnum,
                         cs_gnum_t           n_g_new_vertices)
{
  cs_lnum_t  i, j, k, shift, vid, start, end, perio_start, perio_end;
  cs_lnum_t  n_new_vertices, n_init_faces;
  cs_real_t  matrix[3][4], xyz[4];
  bool  is_modified;
  cs_gnum_t  new_gnum;
  cs_join_state_t  state;

  cs_lnum_t  *new_f2v_idx = NULL, *new_f2v_lst = NULL, *vtag = NULL;
  cs_lnum_t  *linked_id = NULL;
  cs_gnum_t  *gnum = NULL;
  bool  *f_state = NULL;
  cs_join_mesh_t  *work_jmesh = *p_work_jmesh;
  cs_join_edges_t  *work_edges = *p_work_edges;
  cs_join_param_t  param = this_join->param;
  cs_join_select_t  *select = this_join->selection;

  int  perio_id = - 1;
  const int  n_ranks = cs_glob_n_ranks;

  /* Retrieve related back transformation */

  perio_id = fvm_periodicity_get_n_transforms(mesh->periodicity)/2 - 1;
  assert(perio_id > -1);

  fvm_periodicity_get_matrix(mesh->periodicity, 2*perio_id+1, matrix);

  BFT_MALLOC(linked_id, jmesh->n_vertices, cs_lnum_t);
  BFT_MALLOC(gnum, jmesh->n_vertices, cs_gnum_t);

  for (i = 0; i < jmesh->n_vertices; i++) {
    linked_id[i] = -1; /* Default: no link */
    gnum[i] = jmesh->vertices[i].gnum; /* ordered list */
  }

  /* Linked vertex id by inverse transformation */

  for (i = 0; i < select->n_couples; i++) {

    j = cs_search_g_binary(jmesh->n_vertices,
                           select->per_v_couples[2*i],
                           gnum);

    k = cs_search_g_binary(jmesh->n_vertices,
                           select->per_v_couples[2*i+1],
                           gnum);

    linked_id[k] = j;

  }

  BFT_FREE(gnum);

  /* Scan faces to detect new vertices in order to define and build a
     new face->vertex index */

  n_init_faces = jmesh->n_faces/2;

  BFT_MALLOC(f_state, jmesh->n_faces, bool);
  BFT_MALLOC(new_f2v_idx, jmesh->n_faces + 1, cs_lnum_t);
  BFT_MALLOC(vtag, jmesh->n_vertices, cs_lnum_t);

  for (i = 0; i < jmesh->n_vertices; i++)
    vtag[i] = 0;

  for (i = 0; i < n_init_faces; i++) {

    is_modified = false;
    start = jmesh->face_vtx_idx[2*i];
    end = jmesh->face_vtx_idx[2*i+1];
    perio_start = end;
    perio_end = jmesh->face_vtx_idx[2*i+2];

    for (j = perio_start; j < perio_end; j++) {

      vid = jmesh->face_vtx_lst[j];
      state = jmesh->vertices[vid].state;

      if (state == CS_JOIN_STATE_PERIO_MERGE) {

        is_modified = true;
        vtag[vid] = -1;
        assert(linked_id[vid] > -1);

      }
      else if (state == CS_JOIN_STATE_NEW) {

        is_modified = true;
        vtag[vid] = 1;
        assert(linked_id[vid] == -1);

      }
      else if (state == CS_JOIN_STATE_MERGE) {

        is_modified = true;
        if (linked_id[vid] > -1) /* Update is enough */
          vtag[vid] = -2;

        else /* New vertex for the periodic face but not for
                  the original intersected face.
                  Add a new vertex for the related original face */

          vtag[vid] = 2;

      }

    }

    if (is_modified == true)
      new_f2v_idx[2*i+1] = perio_end - perio_start;
    else
      new_f2v_idx[2*i+1] = end - start;
    new_f2v_idx[2*i+2] = perio_end - perio_start;

    f_state[2*i] = is_modified;
    f_state[2*i+1] = false;

  }

  n_new_vertices = 0;
  for (i = 0; i < jmesh->n_vertices; i++)
    if (vtag[i] > 0)
      n_new_vertices++;

  BFT_REALLOC(jmesh->vertices,
              jmesh->n_vertices + n_new_vertices,
              cs_join_vertex_t);

  /* Transform back new periodic vertices */

  n_new_vertices = 0;

  for (i = 0; i < jmesh->n_vertices; i++) {

    if (vtag[i] > 0) { /* New vertex to transform back */

      cs_join_vertex_t  new_vtx = jmesh->vertices[i];

      for (j = 0; j < 3; j++) {
        xyz[j] = new_vtx.coord[j];
        new_vtx.coord[j] = 0.0;
      }
      xyz[3] = 1;

      for (j = 0; j < 3; j++)
        for (k = 0; k < 4; k++)
          new_vtx.coord[j] += matrix[j][k]*xyz[k];

      jmesh->vertices[jmesh->n_vertices + n_new_vertices] = new_vtx;
      n_new_vertices++;
      vtag[i] = jmesh->n_vertices + n_new_vertices;

    }
    else if (vtag[i] < 0) { /* Existing vertex to update */

      cs_join_vertex_t  new_vtx = jmesh->vertices[linked_id[i]];

      assert(   new_vtx.state != CS_JOIN_STATE_MERGE
             || new_vtx.state != CS_JOIN_STATE_PERIO_MERGE);

      for (j = 0; j < 3; j++) {
        xyz[j] = jmesh->vertices[i].coord[j];
        new_vtx.coord[j] = 0.0;
      }
      xyz[3] = 1;

      for (j = 0; j < 3; j++)
        for (k = 0; k < 4; k++)
          new_vtx.coord[j] += matrix[j][k]*xyz[k];

      new_vtx.state = CS_JOIN_STATE_MERGE;
      jmesh->vertices[linked_id[i]] = new_vtx;

    }

  }

  if (n_ranks > 1) { /* Global numbering update */

    fvm_io_num_t  *io_num = NULL;
    const cs_gnum_t  *io_gnum = NULL;

    BFT_MALLOC(gnum, n_new_vertices, cs_gnum_t);

    for (i = 0, shift = jmesh->n_vertices; i < n_new_vertices; i++, shift++)
      gnum[i] = jmesh->vertices[shift].gnum;

    io_num = fvm_io_num_create(NULL, gnum, n_new_vertices, 0);
    io_gnum = fvm_io_num_get_global_num(io_num);

    for (i = 0, shift = jmesh->n_vertices; i < n_new_vertices; i++, shift++) {
      new_gnum = io_gnum[i] + init_max_vtx_gnum + n_g_new_vertices;
      jmesh->vertices[shift].gnum = new_gnum;
    }

    jmesh->n_g_vertices += fvm_io_num_get_global_count(io_num);

    fvm_io_num_destroy(io_num);
    BFT_FREE(gnum);

  }
  else { /* Serial run */

    for (i = 0, shift = jmesh->n_vertices; i < n_new_vertices; i++, shift++) {
      new_gnum = i + 1 + init_max_vtx_gnum + n_g_new_vertices;
      jmesh->vertices[shift].gnum = new_gnum;
    }

    jmesh->n_g_vertices += n_new_vertices;

  }

  /* Update face->vertex connectivity for original faces if needed
     Copy connectivity for periodic faces */

  new_f2v_idx[0] = 0;
  for (i = 0; i < jmesh->n_faces; i++)
    new_f2v_idx[i+1] += new_f2v_idx[i];

  BFT_MALLOC(new_f2v_lst, new_f2v_idx[jmesh->n_faces], cs_lnum_t);

  for (i = 0; i < n_init_faces; i++) {

    start = jmesh->face_vtx_idx[2*i];
    end = jmesh->face_vtx_idx[2*i+1];
    perio_start = end;
    perio_end = jmesh->face_vtx_idx[2*i+2];

    if (f_state[2*i] == false) { /* No modification to apply */

      for (j = start, shift = new_f2v_idx[2*i];
           j < perio_start; j++, shift++)
        new_f2v_lst[shift] = jmesh->face_vtx_lst[j];

    }
    else { /* Modification to apply from the periodic face */

      for (j = perio_start, shift = new_f2v_idx[2*i];
           j < perio_end; j++, shift++) {

        vid = jmesh->face_vtx_lst[j];
        state = jmesh->vertices[vid].state;

        if (   state == CS_JOIN_STATE_PERIO_MERGE
            || state == CS_JOIN_STATE_PERIO)
          new_f2v_lst[shift] = linked_id[vid];
        else if (state == CS_JOIN_STATE_MERGE) {
          if (linked_id[vid] > -1)
            new_f2v_lst[shift] = linked_id[vid];
          else
            new_f2v_lst[shift] = vtag[vid] - 1;
        }
        else if (state == CS_JOIN_STATE_NEW)
          new_f2v_lst[shift] = vtag[vid] - 1;
        else
          bft_error(__FILE__, __LINE__, 0,
                    _("  Vertex state (%d) is not consistent.\n"
                      "  Can not apply changes from periodic faces to"
                      " original face.\n"
                      "  Check your periodicity parameters.\n"), state);
      }

    } /* End of test on face state */

    /* Copy periodic face connectivity */

    for (j = perio_start, shift = new_f2v_idx[2*i+1];
         j < perio_end; j++, shift++)
      new_f2v_lst[shift] = jmesh->face_vtx_lst[j];

  } /* End of loop on initial faces */

  /* Add new vertices to the periodic vertex list */

  if (param.verbosity > 3)
    fprintf(cs_glob_join_log,
            "  Add locally %d new vertices for periodicity\n",
            n_new_vertices);

  shift = select->n_couples;
  select->n_couples += n_new_vertices;
  BFT_REALLOC(select->per_v_couples, 2*select->n_couples, cs_gnum_t);

  for (i = 0; i < jmesh->n_vertices; i++) {
    if (vtag[i] > 0) {
      vid = vtag[i] - 1;
      select->per_v_couples[2*shift] = jmesh->vertices[vid].gnum;
      select->per_v_couples[2*shift+1] = jmesh->vertices[i].gnum;
      shift++;
    }
  }

  BFT_FREE(vtag);
  BFT_FREE(linked_id);
  BFT_FREE(f_state);

  /* Reshape join_mesh structure */

  BFT_FREE(jmesh->face_vtx_idx);
  BFT_FREE(jmesh->face_vtx_lst);

  jmesh->face_vtx_idx = new_f2v_idx;
  jmesh->face_vtx_lst = new_f2v_lst;
  jmesh->n_vertices += n_new_vertices;

  /* Update now work_jmesh by exchanging jmesh over the ranks */

  _redistribute_mesh(param,
                     select->compact_rank_index,
                     jmesh,
                     &work_jmesh);

  /* Define a new cs_join_edges_t structure related to work_jmesh */

  cs_join_mesh_destroy_edges(&work_edges);
  work_edges = cs_join_mesh_define_edges(work_jmesh);

  /* Reshape join_mesh structure by deleting periodic faces */

  shift = 0;

  for (i = 0; i < n_init_faces; i++) {

    jmesh->face_gnum[i] = jmesh->face_gnum[2*i];

    for (j = jmesh->face_vtx_idx[2*i]; j < jmesh->face_vtx_idx[2*i+1]; j++)
      jmesh->face_vtx_lst[shift++] = jmesh->face_vtx_lst[j];

    jmesh->face_vtx_idx[i+1] = shift;

  }

  BFT_REALLOC(jmesh->face_gnum, n_init_faces, cs_gnum_t);
  BFT_REALLOC(jmesh->face_vtx_idx, n_init_faces + 1, cs_lnum_t);
  BFT_REALLOC(jmesh->face_vtx_lst, shift, cs_lnum_t);

  jmesh->n_faces = n_init_faces;
  jmesh->n_g_faces /= 2;

  /* Return pointer */

  *p_work_jmesh = work_jmesh;
  *p_work_edges = work_edges;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (cs_glob_join_log != NULL) {
    fprintf(cs_glob_join_log, "\n Periodic vertex couples:\n");
    for (i = 0; i < select->n_couples; i++)
      fprintf(cs_glob_join_log, " %6d  (%9u, %9u)\n", i+1,
              select->per_v_couples[2*i], select->per_v_couples[2*i+1]);
    fflush(cs_glob_join_log);
  }
#endif

}

/*----------------------------------------------------------------------------
 * Duplicate and apply transformation to the selected faces and also to
 * their related vertices. Update jmesh structure.
 * Define a new n2o_hist.
 *
 * parameters:
 *   this_join  <-- pointer to a high level join structure
 *   jmesh      <-> local join mesh struct. to duplicate and transform
 *   mesh       <-- pointer to a cs_mesh_t structure
 *   builder    <-- pointer to a cs_mesh_builder_t structure
 *   o2n_hist   <-- old global face -> new local face numbering
 *   p_n2o_hist <-- new global face -> old local face numbering
 *---------------------------------------------------------------------------*/

void
cs_join_perio_split_back(cs_join_t          *this_join,
                         cs_join_mesh_t     *jmesh,
                         cs_mesh_t          *mesh,
                         cs_mesh_builder_t  *builder,
                         cs_join_gset_t     *o2n_hist,
                         cs_join_gset_t    **p_n2o_hist)
{
  cs_lnum_t  i, j, k, shift, vid, fid, start, end, perio_start, perio_end;
  cs_lnum_t  n_final_faces, n1_faces, n2_faces;
  cs_lnum_t  shift1, shift2, shift3, shift4;
  cs_lnum_t  n_sub_ori, n_sub_per, n_contrib, n_couples;

  cs_lnum_t  n_vertices_to_add = 0, n_g_vertices_to_add = 0;
  cs_lnum_t  *new_f2v_idx = NULL, *new_f2v_lst = NULL;
  cs_lnum_t  *linked_id = NULL, *f_tag = NULL;
  cs_gnum_t  *gnum = NULL, *f2_gnum = NULL, *new_fgnum = NULL;
  cs_join_gset_t  *new_history = NULL, *n2o_hist = *p_n2o_hist;

  cs_join_select_t  *select = this_join->selection;

  int  perio_id = - 1;
  const int  n_ranks = cs_glob_n_ranks;

  assert(builder != NULL);

  perio_id = fvm_periodicity_get_n_transforms(mesh->periodicity)/2 - 1;
  assert(perio_id > -1);

  /* Detect periodic face to delete and associate a tag for each new face */

  BFT_MALLOC(f_tag, jmesh->n_faces, cs_lnum_t);

  assert(n2o_hist->n_elts == jmesh->n_faces);
  n_couples = 0;

  for (i = 0; i < n2o_hist->n_elts; i++) {

    start = n2o_hist->index[i];
    end = n2o_hist->index[i+1];
    n_contrib = end - start;

    if (n_contrib == 1) { /* New face remains a border face */

      if (n2o_hist->g_list[start]%2 == 0)
        f_tag[i] = 0; /* To delete without transferring back face connect. */
      else
        f_tag[i] = 1; /* To keep. Original face */

    }
    else if (n_contrib == 2) { /* New face becomes an interior face */

      if (n2o_hist->g_list[start]%2 == 0) { /* Old periodic face */

        assert(n2o_hist->g_list[start+1]%2 == 1); /* Associated face should
                                                     be an old original face */

        f_tag[i] = -1; /* To keep with transferring back face connect. */
        n_couples++;

      }
      else { /* Old original face */

        assert(n2o_hist->g_list[start]%2 == 1);
        assert(n2o_hist->g_list[start+1]%2 == 0); /* Associated face should
                                                     be an old periodic face */

        f_tag[i] = -1; /* To keep with transferring back face connect. */
        n_couples++;

      }

    }
    else { /* n_contrib > 2 => remains a border face and keep the new face
              if not all the old faces are periodic */

      bool  have_perio = true;

      f_tag[i] = 0; /* Initialize as if we want to delete the new face */

      for (j = start; j < end; j++)
        if (n2o_hist->g_list[j]%2 == 1) /* Old original face are taking part
                                           in the composition of the new
                                           face */
          have_perio = false;

      if (have_perio == false)
        f_tag[i] = 1; /* To keep without transferring back face connect. */
    }

  } /* End of loop on new faces */

  /* Loop over old faces and change the current tag for the new faces
     according to the number of subdivisions applied to the original
     faces */

  for (i = 0; i < select->n_faces; i++) {

    start = o2n_hist->index[2*i];
    end = o2n_hist->index[2*i+1];
    perio_start = end;
    perio_end = o2n_hist->index[2*i+2];

    n_sub_ori = end - start;
    n_sub_per = perio_end - perio_start;

    assert(n_sub_per > 0 && n_sub_ori > 0);

    if (n_sub_ori == 1 && n_sub_per > 1) {

      fid = cs_search_g_binary(jmesh->n_faces,
                               o2n_hist->g_list[start],
                               jmesh->face_gnum);

      f_tag[fid] = 2; /* Original face to delete and to replace by
                         a set of new sub-faces */

      for (j = perio_start; j < perio_end; j++) {

        fid = cs_search_g_binary(jmesh->n_faces,
                                 o2n_hist->g_list[j],
                                 jmesh->face_gnum);

        assert(fid > -1);

        if (f_tag[fid] == 0)
          f_tag[fid] = -2; /* Switch tag:
                              Will be deleted after periodicity application */

      }

    } /* n_sub_ori == 1 && n_sub_per > 1 */

    else if (n_sub_ori == 1 && n_sub_per == 1) {

      fid = cs_search_g_binary(jmesh->n_faces,
                               o2n_hist->g_list[perio_start],
                               jmesh->face_gnum);

      if (n2o_hist->index[fid+1] - n2o_hist->index[fid] == 2) {

        assert(f_tag[fid] == -1); /* New face to keep with transferring
                                     back face connectivity */

        fid = cs_search_g_binary(jmesh->n_faces,
                                 o2n_hist->g_list[start],
                                 jmesh->face_gnum);


        assert(fid > -1);
        f_tag[fid] = 2; /* Original face to delete and to replace by
                           a set the new sub-face */

      }

    }

  } /* End of loop on old faces */

  /* Count the final number of new subfaces */

  n_final_faces = 0;
  for (i = 0; i < jmesh->n_faces; i++) {

    if (f_tag[i] == 1 || f_tag[i] == -2)
      n_final_faces++;
    else if (f_tag[i] == -1)
      n_final_faces += 2;

  }

  /* Define a global number for each face to transfer back by periodicity */

  n2_faces = 0;
  for (i = 0; i < jmesh->n_faces; i++)
    if (f_tag[i] < 0)
      n2_faces++;

  BFT_MALLOC(f2_gnum, n2_faces, cs_gnum_t);

  if (n_ranks > 1) { /* Parallel run */

    fvm_io_num_t  *io_num = NULL;
    const cs_gnum_t  *io_gnum = NULL;

    n2_faces = 0;
    for (i = 0; i < jmesh->n_faces; i++)
      if (f_tag[i] < 0)
        f2_gnum[n2_faces++] = jmesh->face_gnum[i];

    io_num = fvm_io_num_create(NULL, f2_gnum, n2_faces, 0);
    io_gnum = fvm_io_num_get_global_num(io_num);

    for (i = 0; i < n2_faces; i++)
      f2_gnum[i] = jmesh->n_g_faces + io_gnum[i];

    fvm_io_num_destroy(io_num);

  }
  else { /* Serial run */

    n2_faces = 0;
    for (i = 0; i < jmesh->n_faces; i++) {
      if (f_tag[i] < 0) {
        f2_gnum[n2_faces] = jmesh->n_faces + n2_faces + 1;
        n2_faces++;
      }
    }

  }

  /* Define the new face -> vertex index and the new global face numbering */

  BFT_MALLOC(new_f2v_idx, n_final_faces + 1, cs_lnum_t);
  BFT_MALLOC(new_fgnum, n_final_faces, cs_gnum_t);

  new_history = cs_join_gset_create(n_final_faces);
  new_history->n_g_elts = new_history->n_elts;

  n1_faces = 0;

  for (i = 0; i < jmesh->n_faces; i++) {
    if (f_tag[i] == 1 || f_tag[i] == -1) {

      new_f2v_idx[n1_faces+1]
        = jmesh->face_vtx_idx[i+1] - jmesh->face_vtx_idx[i];
      new_fgnum[n1_faces] = jmesh->face_gnum[i];

      new_history->index[n1_faces+1] =
        n2o_hist->index[i+1] - n2o_hist->index[i];

      n1_faces++;

    }
  }

  shift = n1_faces;
  n2_faces = 0;
  for (i = 0; i < jmesh->n_faces; i++) {
    if (f_tag[i] < 0) {

      assert(f_tag[i] == -1 || f_tag[i] == -2);

      new_f2v_idx[shift + 1]
        = jmesh->face_vtx_idx[i+1] - jmesh->face_vtx_idx[i];
      new_fgnum[shift] = f2_gnum[n2_faces];

      new_history->index[shift+1] =
        n2o_hist->index[i+1] - n2o_hist->index[i];

      n2_faces++;
      shift++;

    }
  }

  assert(n1_faces + n2_faces == n_final_faces);

  /* Detect if there are new vertices to add to the jmesh definition */

  BFT_MALLOC(gnum, jmesh->n_vertices, cs_gnum_t);
  BFT_MALLOC(linked_id, jmesh->n_vertices, cs_lnum_t);

  for (i = 0; i < jmesh->n_vertices; i++) {
    linked_id[i] = -1; /* Default: no link */
    gnum[i] = jmesh->vertices[i].gnum; /* ordered list */
  }

  for (i = 0; i < select->n_couples; i++) {

    j = cs_search_g_binary(jmesh->n_vertices,
                           select->per_v_couples[2*i],
                           gnum);

    k = cs_search_g_binary(jmesh->n_vertices,
                           select->per_v_couples[2*i+1],
                           gnum);

    linked_id[k] = j;

  }

  BFT_FREE(f2_gnum);
  BFT_FREE(gnum);

  for (i = 0; i < jmesh->n_faces; i++) {
    if (f_tag[i] < 0) { /* Periodicity to transfer back */

      for (j = jmesh->face_vtx_idx[i]; j < jmesh->face_vtx_idx[i+1]; j++) {
        vid = jmesh->face_vtx_lst[j];

        if (linked_id[vid] == -1)
          linked_id[vid] = -2; /* Have to get the related id. Add a new
                                  vertex */

      }

    }
  } /* End of loop on new faces */

  n_vertices_to_add = 0;
  for (i = 0; i < jmesh->n_vertices; i++)
    if (linked_id[i] == -2)
      n_vertices_to_add += 1;

  if (n_vertices_to_add > 0) {

    cs_lnum_t  i1, i2;
    cs_real_t  matrix[3][4], xyz[4];

    BFT_REALLOC(jmesh->vertices, jmesh->n_vertices + n_vertices_to_add,
                cs_join_vertex_t);
    BFT_REALLOC(linked_id,  jmesh->n_vertices + n_vertices_to_add, cs_lnum_t);

    /* Retrieve related back transformation */

    fvm_periodicity_get_matrix(mesh->periodicity, 2*perio_id+1, matrix);

    n_vertices_to_add = 0;
    for (i1 = 0; i1 < jmesh->n_vertices; i1++) {
      if (linked_id[i1] == -2) {

        cs_join_vertex_t  new_vtx = jmesh->vertices[i1];

        i2 = jmesh->n_vertices + n_vertices_to_add;
        linked_id[i1] = i2;

        for (k = 0; k < 3; k++) {
          xyz[k] = new_vtx.coord[k];
          new_vtx.coord[k] = 0.0;
        }
        xyz[3] = 1;

        for (j = 0; j < 3; j++)
          for (k = 0; k < 4; k++)
            new_vtx.coord[j] += matrix[j][k]*xyz[k];

        jmesh->vertices[i2] = new_vtx;
        n_vertices_to_add += 1;

      }
    }

  } /* n_vertices_to_add > 0 */

  /* Define a global vertex num for the new added vertices */

#if defined(HAVE_MPI)
  if (n_ranks > 1) {

    fvm_io_num_t  *io_num = NULL;
    const cs_gnum_t  *io_gnum = NULL;

    MPI_Allreduce(&n_vertices_to_add, &n_g_vertices_to_add, 1, MPI_INT,
                  MPI_SUM, cs_glob_mpi_comm);

    /* Define global numbering for vertices to add
       and keep relation between periodic couples */

    if (n_g_vertices_to_add > 0) {

      BFT_MALLOC(gnum, n_vertices_to_add, cs_gnum_t);

      for (i = 0, shift = jmesh->n_vertices; i < n_vertices_to_add;
           i++, shift++)
        gnum[i] = jmesh->vertices[shift].gnum;

      io_num = fvm_io_num_create(NULL, gnum, n_vertices_to_add, 0);
      io_gnum = fvm_io_num_get_global_num(io_num);

      for (i = 0, shift = jmesh->n_vertices; i < n_vertices_to_add;
           i++, shift++)
        jmesh->vertices[shift].gnum = io_gnum[i] + mesh->n_g_vertices;

      n_g_vertices_to_add = fvm_io_num_get_global_count(io_num);
      jmesh->n_g_vertices += n_g_vertices_to_add;

      fvm_io_num_destroy(io_num);
      BFT_FREE(gnum);

    }

  }
#endif

  if (n_ranks == 1 && n_vertices_to_add > 0) {

    for (i = 0, shift = jmesh->n_vertices; i < n_vertices_to_add; i++, shift++)
      jmesh->vertices[shift].gnum = i + 1 + mesh->n_g_vertices;

    jmesh->n_g_vertices += n_vertices_to_add;

  }

  jmesh->n_vertices += n_vertices_to_add;

  new_f2v_idx[0] = 0;
  new_history->index[0] = 0;

  for (i = 0; i < n_final_faces; i++) {
    new_f2v_idx[i+1] += new_f2v_idx[i];
    new_history->index[i+1] += new_history->index[i];
  }

  BFT_MALLOC(new_f2v_lst, new_f2v_idx[n_final_faces], cs_lnum_t);
  BFT_MALLOC(new_history->g_list, new_history->index[new_history->n_elts],
             cs_gnum_t);

  /* Define a new face connectivity and a new face history */

  shift1 = 0; /*kept faces in face connect. */
  shift2 = new_f2v_idx[n1_faces]; /* transfered faces in face connect. */

  shift3 = 0; /* kept faces in face history */
  shift4 = new_history->index[n1_faces]; /* transfered faces in face history */

  /* Store periodic couples. First in local join numbering.
     Will move next to a global numbering (after interior face add) */

  builder->n_per_face_couples[perio_id] = n_couples;
  BFT_MALLOC(builder->per_face_couples[perio_id], 2*n_couples, cs_gnum_t);

  n2_faces = n1_faces;
  n1_faces = 0;
  n_couples = 0;

  for (i = 0; i < jmesh->n_faces; i++) {

    start = jmesh->face_vtx_idx[i];
    end = jmesh->face_vtx_idx[i+1];

    if (f_tag[i] == 1) { /* Original face to keep */

      for (j = start; j < end; j++)
        new_f2v_lst[shift1++] = jmesh->face_vtx_lst[j];

      for (j = n2o_hist->index[i]; j < n2o_hist->index[i+1]; j++)
        new_history->g_list[shift3++] = n2o_hist->g_list[j];

      n1_faces++;

    }
    else if (f_tag[i] == -1) { /* Periodic face to keep and to transfer back */

      /* Store periodic face couple */

      n1_faces++;
      n2_faces++;
      builder->per_face_couples[perio_id][2*n_couples] = n2_faces;
      builder->per_face_couples[perio_id][2*n_couples+1] = n1_faces;
      n_couples++;

      /* Define face connectivity */

      for (j = start; j < end; j++) {
        vid = jmesh->face_vtx_lst[j];
        assert(linked_id[vid] > -1);
        new_f2v_lst[shift1++] = vid;
        new_f2v_lst[shift2++] = linked_id[vid];
      }

      /* Define new old->new face history */

      for (j = n2o_hist->index[i]; j < n2o_hist->index[i+1]; j++) {

        if (n2o_hist->g_list[j]%2 == 0) { /* Periodic face */
          new_history->g_list[shift3++] = n2o_hist->g_list[j];
          /* Original face related to this periodic one */
          new_history->g_list[shift4++] = n2o_hist->g_list[j] - 1;
        }
        else { /* Original face */
          new_history->g_list[shift3++] = n2o_hist->g_list[j];
          /* Periodic face related to this original one */
          new_history->g_list[shift4++] = n2o_hist->g_list[j] + 1;
        }

      }

    }
    else if (f_tag[i] == -2) { /* Periodic face to transfer back */

      n2_faces++;

      for (j = start; j < end; j++) {
        vid = jmesh->face_vtx_lst[j];
        assert(linked_id[vid] > -1);
        new_f2v_lst[shift2++] = linked_id[vid];
      }

      for (j = n2o_hist->index[i]; j < n2o_hist->index[i+1]; j++) {

        assert(n2o_hist->g_list[j]%2 == 0); /* Periodic face */

        /* Original face related to this periodic one */
        new_history->g_list[shift4++] = n2o_hist->g_list[j] - 1;

      }

    }

  } /* End of loop on faces */

  BFT_FREE(linked_id);
  BFT_FREE(f_tag);

  /* Reshape join_mesh structure */

  jmesh->n_faces = n_final_faces;
  jmesh->n_g_faces = jmesh->n_faces;

  BFT_FREE(jmesh->face_gnum);
  BFT_FREE(jmesh->face_vtx_idx);
  BFT_FREE(jmesh->face_vtx_lst);

  jmesh->face_gnum = new_fgnum;
  jmesh->face_vtx_idx = new_f2v_idx;
  jmesh->face_vtx_lst = new_f2v_lst;

  if (n_ranks > 1) { /* Update global face number */

    fvm_io_num_t  *io_num = NULL;
    const cs_gnum_t  *io_gnum = NULL;

    io_num = fvm_io_num_create(NULL, jmesh->face_gnum, jmesh->n_faces, 0);
    io_gnum = fvm_io_num_get_global_num(io_num);

    for (i = 0; i < jmesh->n_faces; i++) {
      jmesh->face_gnum[i] = io_gnum[i];
      new_history->g_elts[i] = io_gnum[i];
    }

    jmesh->n_g_faces = fvm_io_num_get_global_count(io_num);
    new_history->n_g_elts = jmesh->n_g_faces;

    fvm_io_num_destroy(io_num);
    BFT_FREE(gnum);

  }
  else { /* Serial run */

    for (i = 0; i < jmesh->n_faces; i++) {
      jmesh->face_gnum[i] = i+1;
      new_history->g_elts[i] = i+1;
    }

  }

  /* Return pointer */

  cs_join_gset_destroy(&n2o_hist);

  *p_n2o_hist = new_history;

}

/*----------------------------------------------------------------------------
 * Define a list of coupled faces by periodicty in global numbering.
 *
 * For parallel runs:
 *  - remove isolated periodic faces in the mesh definition
 *  - define a consistent face connectivity in order to prepare the building
 *    of periodic vertex couples
 *
 * parameters:
 *   param        <-- set of parameters for the joining operation
 *   n_ii_faces   <-- initial local number of interior faces
 *   face_type    <-- type of faces in join mesh (interior or border ...)
 *   jmesh        <-- pointer to a cs_join_mesh_t structure
 *   mesh         <-> pointer to a cs_mesh_t structure
 *   mesh_builder <-> pointer to a cs_mesh_t structure
 *---------------------------------------------------------------------------*/

void
cs_join_perio_split_update(cs_join_param_t             param,
                           cs_lnum_t                   n_ii_faces,
                           const cs_join_face_type_t   face_type[],
                           const cs_join_mesh_t       *jmesh,
                           cs_mesh_t                  *mesh,
                           cs_mesh_builder_t          *mesh_builder)
{
  int  i, shift;

  cs_gnum_t  *o2n_num = NULL;

  int  perio_id = - 1;
  const int  n_j_faces = jmesh->n_faces;
  const int  n_ranks = cs_glob_n_ranks;

  assert(mesh_builder != NULL);

  perio_id = fvm_periodicity_get_n_transforms(mesh->periodicity)/2 - 1;
  assert(perio_id > -1);

  /* Initialize o2n_num */

  BFT_MALLOC(o2n_num, n_j_faces, cs_gnum_t);

  for (i = 0; i < n_j_faces; i++)
    o2n_num[i] = 0;

  if (n_ranks == 1) {

    shift = n_ii_faces + 1;
    for (i = 0; i < n_j_faces; i++)
      if (face_type[i] == CS_JOIN_FACE_INTERIOR)
        o2n_num[i] = shift++;

  }
  else { /* n_ranks > 1 */

    shift = n_ii_faces;
    for (i = 0; i < n_j_faces; i++)
      if (face_type[i] == CS_JOIN_FACE_INTERIOR)
        o2n_num[i] = mesh->global_i_face_num[shift++];

  }

  /* Apply new numbering */

  for (i = 0; i < mesh_builder->n_per_face_couples[perio_id]; i++) {

    cs_gnum_t  old1 = mesh_builder->per_face_couples[perio_id][2*i] - 1;
    cs_gnum_t  old2 = mesh_builder->per_face_couples[perio_id][2*i+1] - 1;

    assert(o2n_num[old1] > 0);
    assert(o2n_num[old2] > 0);

    mesh_builder->per_face_couples[perio_id][2*i] = o2n_num[old1];
    mesh_builder->per_face_couples[perio_id][2*i+1] = o2n_num[old2];

  }

  BFT_FREE(o2n_num);

  if (n_ranks > 1) /* Remove isolated periodic face for the current mesh */
    _perio_face_clean(param, mesh);

}

/*---------------------------------------------------------------------------*/

END_C_DECLS

