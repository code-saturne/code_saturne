/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
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
 * Management of the post-processing
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_config.h>
#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_nodal.h>
#include <fvm_writer.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_mesh_connect.h"
#include "cs_prototypes.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_post.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local types and structures
 *============================================================================*/

/* Mesh location type */
/*--------------------*/

typedef enum {

  CS_POST_LOCATION_CELL,         /* Values located at cells */
  CS_POST_LOCATION_I_FACE,       /* Values located at interior faces */
  CS_POST_LOCATION_B_FACE,       /* Values located at boundary faces */
  CS_POST_LOCATION_VERTEX        /* Values located at vertices */

} cs_post_location_t;

/* Writer structure */
/*------------------*/

/* This object is based on a choice of a case, directory, and format,
   as well as a flag for associated mesh's time dependency, and the default
   output frequency for associated variables. */

typedef struct {

  int            id;           /* Identifier (< 0 for "reservable" writer,
                                * > 0 for user writer */
  int            frequency;    /* Default output frequency */
  cs_bool_t      write_displ;  /* Write displacement field if true */

  int            active;       /* 0 if no output at current time step,
                                  1 in case of output */

  fvm_writer_t  *writer;       /* Associated FVM writer */

} cs_post_writer_t;

/* Post-processing mesh structure */
/*--------------------------------*/

/* This object manages the link between an exportable mesh and
   associated writers. */

typedef struct {

  int                     id;            /* Identifier (< 0 for "reservable"
                                            mesh, > 0 for user mesh */

  int                     ent_flag[3];   /* Presence of cells (ent_flag[0],
                                            interior faces (ent_flag[1]),
                                            or boundary faces (ent_flag[2])
                                            on one processor at least */

  int                     cat_id;        /* Optional category id for
                                            as regards variables output
                                            (-1 as base volume mesh, -2 as
                                            base boundary mesh, identical
                                            to id by default) */

  int                     alias;         /* If > -1, index in array of
                                            post-processing meshes of the
                                            first mesh sharing the same
                                            exportable mesh */

  int                     n_writers;     /* Number of associated writers */
  int                    *writer_id;     /* Array of associated writer ids */
  int                     nt_last;       /* Time step number for the last
                                            output (-1 before first output) */

  cs_int_t                n_i_faces;     /* N. associated interior faces */
  cs_int_t                n_b_faces;     /* N. associated boundary faces */

  const fvm_nodal_t      *exp_mesh;      /* Associated exportable mesh */
  fvm_nodal_t            *_exp_mesh;     /* Associated exportble mesh,
                                            if owner */

  fvm_writer_time_dep_t   mod_flag_min;  /* Minimum mesh time dependency */
  fvm_writer_time_dep_t   mod_flag_max;  /* Maximum mesh time dependency */

} cs_post_mesh_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Backup of initial vertex coordinates */

static cs_bool_t   _cs_post_deformable = false;
static cs_real_t  *_cs_post_ini_vtx_coo = NULL;

/* Flag to indicate output of domain number in parallel mode */

static cs_bool_t        _cs_post_domain = true;

/* Array of exportable meshes associated with post-processing */
/* (meshes -1 and -2 reserved, so free ids start at -2)*/

static int              _cs_post_min_mesh_id = -2;
static int              _cs_post_n_meshes = 0;
static int              _cs_post_n_meshes_max = 0;
static cs_post_mesh_t  *_cs_post_meshes = NULL;

/* Array of writers for post-processing; */
/* writers -1 (default) and -2 (show errors) are reserved */

static int                _cs_post_min_writer_id = -2;
static int                _cs_post_n_writers = 0;
static int                _cs_post_n_writers_max = 0;
static cs_post_writer_t  *_cs_post_writers = NULL;

/* Array of registered variable output functions and instances */

static int                _cs_post_nbr_var_tp = 0;
static int                _cs_post_nbr_var_tp_max = 0;

static cs_post_time_dep_var_t **_cs_post_f_var_tp = NULL;
static int                     *_cs_post_i_var_tp = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Convert cs_post_type_t datatype to fvm_datatype_t.
 *
 * parameters:
 *   type_cs <-- Code_Saturne data type
 *
 * returns
 *   corresponding FVM datatype
 *----------------------------------------------------------------------------*/

static fvm_datatype_t
_cs_post_cnv_datatype(cs_post_type_t  type_cs)
{
  fvm_datatype_t type_fvm = FVM_DATATYPE_NULL;

  switch(type_cs) {

  case CS_POST_TYPE_cs_int_t:
    if (sizeof(cs_int_t) == 4)
      type_fvm = FVM_INT32;
    else if (sizeof(cs_int_t) == 8)
      type_fvm = FVM_INT64;
    break;

  case CS_POST_TYPE_cs_real_t:
    if (sizeof(cs_real_t) == sizeof(double))
      type_fvm = FVM_DOUBLE;
    else if (sizeof(cs_real_t) == sizeof(float))
      type_fvm = FVM_FLOAT;
    break;

  case CS_POST_TYPE_int:
    if (sizeof(int) == 4)
      type_fvm = FVM_INT32;
    else if (sizeof(int) == 8)
      type_fvm = FVM_INT64;
    break;

  case CS_POST_TYPE_float:
    type_fvm = FVM_FLOAT;
    break;

  case CS_POST_TYPE_double:
    type_fvm = FVM_DOUBLE;
    break;

  default:
    assert(0);
  }

  return type_fvm;
}

/*----------------------------------------------------------------------------
 * Search for position in the array of writers of a writer with a given id.
 *
 * parameters:
 *   writer_id <-- id of writer
 *
 * returns:
 *   position in the array of writers
 *----------------------------------------------------------------------------*/

static int
_cs_post_writer_id(const int  writer_id)
{
  cs_int_t  id;

  cs_post_writer_t  *writer = NULL;

  /* Search for requested writer */

  for (id = 0; id < _cs_post_n_writers; id++) {
    writer = _cs_post_writers + id;
    if (writer->id == writer_id)
      break;
  }
  if (id >= _cs_post_n_writers)
    bft_error(__FILE__, __LINE__, 0,
              _("The requested post-processing writer number\n"
                "%d is not defined.\n"), (int)(writer_id));

  return id;
}

/*----------------------------------------------------------------------------
 * Search for position in the array of meshes of a mesh with a given id.
 *
 * parameters:
 *   mesh_id <-- id of mesh
 *
 * returns:
 *   position in the array of meshes
 *----------------------------------------------------------------------------*/

static int
_cs_post_mesh_id(int  mesh_id)
{
  int id;
  cs_post_mesh_t  *post_mesh = NULL;

  /* Search for requested mesh */

  for (id = 0; id < _cs_post_n_meshes; id++) {
    post_mesh = _cs_post_meshes + id;
    if (post_mesh->id == mesh_id)
      break;
  }
  if (id >= _cs_post_n_meshes)
    bft_error(__FILE__, __LINE__, 0,
              _("The requested post-processing mesh number\n"
                "%d is not defined.\n"), (int)mesh_id);

  return id;
}

/*----------------------------------------------------------------------------
 * Add a post-processing mesh, do basic initialization, and return a pointer
 * to the associated structure
 *
 * parameters:
 *   mesh_id <-- requested mesh id
 *
 * returns:
 *   pointer to associated structure
 *----------------------------------------------------------------------------*/

static cs_post_mesh_t *
_cs_post_add_mesh(int  mesh_id)
{
  /* local variables */

  int    i, j;

  cs_post_mesh_t  *post_mesh = NULL;

  /* Check that the requested mesh is available */

  if (mesh_id == 0)
      bft_error(__FILE__, __LINE__, 0,
                _("The requested post-processing mesh number\n"
                  "must be < 0 (reserved) or > 0 (user).\n"));

  for (i = 0; i < _cs_post_n_meshes; i++) {
    if ((_cs_post_meshes + i)->id == mesh_id)
      bft_error(__FILE__, __LINE__, 0,
                _("The requested post-processing mesh number\n"
                  "(%d) has already been assigned.\n"), (int)mesh_id);
  }

  /* Resize global array of exportable meshes */

  if (_cs_post_n_meshes == _cs_post_n_meshes_max) {

    if (_cs_post_n_meshes_max == 0)
      _cs_post_n_meshes_max = 8;
    else
      _cs_post_n_meshes_max *= 2;

    BFT_REALLOC(_cs_post_meshes,
                _cs_post_n_meshes_max,
                cs_post_mesh_t);

  }

  _cs_post_n_meshes += 1;

  if (mesh_id < _cs_post_min_mesh_id)
    _cs_post_min_mesh_id = mesh_id;

  /* Assign newly created mesh to the structure */

  post_mesh = _cs_post_meshes + _cs_post_n_meshes - 1;

  post_mesh->id = mesh_id;
  post_mesh->cat_id = mesh_id;
  post_mesh->alias = -1;

  post_mesh->n_writers = 0;
  post_mesh->writer_id = NULL;

  post_mesh->nt_last = -1;

  for (j = 0; j < 3; j++)
    post_mesh->ent_flag[j] = 0;

  post_mesh->n_i_faces = 0;
  post_mesh->n_b_faces = 0;

  post_mesh->exp_mesh = NULL;
  post_mesh->_exp_mesh = NULL;

  /* Minimum and maximum time dependency flags initially inverted,
     will be recalculated after mesh - writer associations */

  post_mesh->mod_flag_min = FVM_WRITER_TRANSIENT_CONNECT;
  post_mesh->mod_flag_max = FVM_WRITER_FIXED_MESH;

  return post_mesh;
}

/*----------------------------------------------------------------------------
 * Create a post-processing mesh; lists of cells or faces to extract are
 * sorted upon exit, whether they were sorted upon calling or not.
 *
 * The list of associated cells is only necessary if the number of cells
 * to extract is strictly greater than 0 and less than the number of cells
 * of the computational mesh.
 *
 * Lists of faces are ignored if the number of extracted cells is nonzero;
 * otherwise, if the number of boundary faces to extract is equal to the
 * number of boundary faces in the computational mesh, and the number of
 * interior faces to extract is zero, than we extrac by default the boundary
 * mesh, and the list of associated boundary faces is thus not necessary.
 *
 * parameters:
 *   post_mesh   <-> pointer to partially initialized post-processing mesh
 *   mesh_name   <-- associated mesh name
 *   n_cells     <-- number of associated cells
 *   n_i_faces   <-- number of associated interior faces
 *   n_b_faces   <-- number of associated boundary faces
 *   cell_list   <-> list of associated cells
 *   i_face_list <-> list of associated interior faces
 *   b_face_list <-> list of associated boundary faces
 *----------------------------------------------------------------------------*/

static void
_cs_post_define_mesh(cs_post_mesh_t  *post_mesh,
                     const char      *mesh_name,
                     cs_int_t         n_cells,
                     cs_int_t         n_i_faces,
                     cs_int_t         n_b_faces,
                     cs_int_t         cell_list[],
                     cs_int_t         i_face_list[],
                     cs_int_t         b_face_list[])
{
  /* local variables */

  int    i;
  int    glob_flag[5];

  int    loc_flag[5] = {1, 1, 1, 0, 0};  /* Flags 0 to 2 "inverted" compared
                                            to others so as to use a single
                                            call to
                                            MPI_Allreduce(..., MPI_MIN, ...) */

  fvm_nodal_t  *exp_mesh = NULL;
  cs_bool_t     maj_ent_flag = false;

  /* Flags:
     0: 0 if cells present, 1 if none,
     1: 0 if interior faces present, 1 if none,
     2: 0 if boundary faces present, 1 if none,
     3: 1 if all cells were selected,
     4: 1 if all boundary faces and no interior faces selected */

  if (n_cells > 0)
    loc_flag[0] = 0;
  else {
    if (n_i_faces > 0)
      loc_flag[1] = 0;
    if (n_b_faces > 0)
      loc_flag[2] = 0;
  }

  if (n_cells >= cs_glob_mesh->n_cells)
    loc_flag[3] = 1;
  else
    loc_flag[3] = 0;

  if (   n_b_faces >= cs_glob_mesh->n_b_faces
      && n_i_faces == 0)
    loc_flag[4] = 1;
  else
    loc_flag[4] = 0;

  for (i = 0; i < 5; i++)
    glob_flag[i] = loc_flag[i];

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)
    MPI_Allreduce (loc_flag, glob_flag, 5, MPI_INT, MPI_MIN,
                   cs_glob_mpi_comm);
#endif

  /* Create associated structure */

  if (glob_flag[0] == 0) {

    if (glob_flag[3] == 1)
      exp_mesh = cs_mesh_connect_cells_to_nodal(cs_glob_mesh,
                                                mesh_name,
                                                cs_glob_mesh->n_cells,
                                                NULL);
    else
      exp_mesh = cs_mesh_connect_cells_to_nodal(cs_glob_mesh,
                                                mesh_name,
                                                n_cells,
                                                cell_list);

  }
  else {

    if (glob_flag[4] == 1)
      exp_mesh = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                                mesh_name,
                                                0,
                                                cs_glob_mesh->n_b_faces,
                                                NULL,
                                                NULL);
    else
      exp_mesh = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                                mesh_name,
                                                n_i_faces,
                                                n_b_faces,
                                                i_face_list,
                                                b_face_list);

  }

  /* Global indicators of mesh entity type presence;
     updated only if the mesh is not totally empty (for time-depending
     meshes, empty at certain times, we want to know the last type
     of entity used in USMPST) */

  for (i = 0; i < 3; i++) {
    if (glob_flag[i] == 0)
      maj_ent_flag = true;
  }

  if (maj_ent_flag == true) {
    for (i = 0; i < 3; i++) {
      if (glob_flag[i] == 0)         /* Inverted glob_flag 0 to 2 logic */
        post_mesh->ent_flag[i] = 1;  /* (c.f. remark above) */
      else
        post_mesh->ent_flag[i] = 0;
    }
  }

  /* Local dimensions */

  post_mesh->n_i_faces = n_i_faces;
  post_mesh->n_b_faces = n_b_faces;

  /* Link to newly created mesh */

  post_mesh->exp_mesh = exp_mesh;
  post_mesh->_exp_mesh = exp_mesh;
}

/*----------------------------------------------------------------------------
 * Update mesh time dependency flags in case of an alias based on the
 * associated writer properties:
 *
 * A mesh's definition may not be modified if the minimum time dependency
 * flag is too low (i.e. if one of the associated writers does not allow
 * changing a mesh's topology).
 *
 * Vertex coordinates and connectivity can be freed from memory if the
 * maximum time dependency flag is low enough (i.e. if none of the associated
 * writers allows modification of the mesh, and thus its future output).
 *
 * parameters:
 *   mesh_id <-- associated mesh (alias) id
 *----------------------------------------------------------------------------*/

static void
_cs_post_mod_flag_alias(int  mesh_id)
{
  int  i;

  cs_post_mesh_t  *post_mesh = NULL;
  cs_post_mesh_t  *ref_mesh = NULL;

  /* Update reference */

  post_mesh = _cs_post_meshes + mesh_id;

  if (post_mesh->alias > -1) {

    ref_mesh = _cs_post_meshes + post_mesh->alias;

    if (post_mesh->mod_flag_min < ref_mesh->mod_flag_min)
      ref_mesh->mod_flag_min = post_mesh->mod_flag_min;

    if (post_mesh->mod_flag_max < ref_mesh->mod_flag_max)
      ref_mesh->mod_flag_max = post_mesh->mod_flag_max;

  }

  /* Update alias */

  for (i = 0; i < _cs_post_n_meshes; i++) {

    post_mesh = _cs_post_meshes + i;

    if (post_mesh->alias > -1) {

      ref_mesh = _cs_post_meshes + post_mesh->alias;

      if (post_mesh->mod_flag_min > ref_mesh->mod_flag_min)
        post_mesh->mod_flag_min = ref_mesh->mod_flag_min;

      if (post_mesh->mod_flag_max > ref_mesh->mod_flag_max)
        post_mesh->mod_flag_max = ref_mesh->mod_flag_max;
    }

  }

}

/*----------------------------------------------------------------------------
 * Divide polygons or polyhedra in simpler elements if necessary.
 *
 * parameters:
 *   post_mesh <-> pointer to post-processing mesh
 *   writer    <-- pointer to associated writer
 *----------------------------------------------------------------------------*/

static void
_cs_post_divide_poly(cs_post_mesh_t          *post_mesh,
                     const cs_post_writer_t  *writer)
{
  /* Divide polygons or polyhedra into simple elements */

  if (fvm_writer_needs_tesselation(writer->writer,
                                   post_mesh->exp_mesh,
                                   FVM_CELL_POLY) > 0)
    fvm_nodal_tesselate(post_mesh->_exp_mesh,
                        FVM_CELL_POLY,
                        NULL);

  if (fvm_writer_needs_tesselation(writer->writer,
                                   post_mesh->exp_mesh,
                                   FVM_FACE_POLY) > 0)
    fvm_nodal_tesselate(post_mesh->_exp_mesh,
                        FVM_FACE_POLY,
                        NULL);
}

/*----------------------------------------------------------------------------
 * Assemble variable values defined on a mix of interior and boundary
 * faces (with no indirection) into an array defined on a single faces set.
 *
 * The resulting variable is not interlaced.
 *
 * parameters:
 *   exp_mesh    <-- exportable mesh
 *   n_i_faces   <-- number of interior faces
 *   n_b_faces   <-- number of boundary faces
 *   var_dim     <-- varible dimension
 *   interlace   <-- for vector, interlace if 1, no interlace if 0
 *   i_face_vals <-- values at interior faces
 *   b_face_vals <-- values at boundary faces
 *   var_tmp[]   --> assembled values
 *----------------------------------------------------------------------------*/

static void
_cs_post_assmb_var_faces(const fvm_nodal_t  *exp_mesh,
                         cs_int_t            n_i_faces,
                         cs_int_t            n_b_faces,
                         int                 var_dim,
                         fvm_interlace_t     interlace,
                         const cs_real_t     i_face_vals[],
                         const cs_real_t     b_face_vals[],
                         cs_real_t           var_tmp[])
{
  cs_int_t  i, j, stride_1, stride_2;

  cs_int_t  n_elts = n_i_faces + n_b_faces;

  assert(exp_mesh != NULL);

  /* The variable is defined on interior and boundary faces of the
     post-processing mesh, and has been built using values
     at the corresponding interior and boundary faces */

  /* Boundary faces contribution */

  if (interlace == FVM_INTERLACE) {
    stride_1 = var_dim;
    stride_2 = 1;
  }
  else {
    stride_1 = 1;
    stride_2 = n_b_faces;
  }

  for (i = 0; i < n_b_faces; i++) {
    for (j = 0; j < var_dim; j++)
      var_tmp[i + j*n_elts] = b_face_vals[i*stride_1 + j*stride_2];
  }

  /* Interior faces contribution */

  if (interlace == FVM_INTERLACE) {
    stride_1 = var_dim;
    stride_2 = 1;
  }
  else {
    stride_1 = 1;
    stride_2 = n_i_faces;
  }

  for (i = 0; i < n_i_faces; i++) {
    for (j = 0; j < var_dim; j++)
      var_tmp[i + n_b_faces + j*n_elts] = i_face_vals[i*stride_1 + j*stride_2];
  }

}

/*----------------------------------------------------------------------------
 * Write parallel domain (rank) number to post-processing mesh
 *
 * parameters:
 *   writer     <-- FVM writer
 *   exp_mesh   <-- exportable mesh
 *   nt_cur_abs <-- current time step number
 *   t_cur_abs  <-- current physical time
 *---------------------------------------------------------------------------*/

static void
_cs_post_write_domain(fvm_writer_t       *writer,
                      const fvm_nodal_t  *exp_mesh,
                      int                 nt_cur_abs,
                      double              t_cur_abs)
{
  int  dim_ent;
  fvm_lnum_t  i, n_elts;
  fvm_datatype_t  datatype;

  fvm_lnum_t   dec_num_parent[1]  = {0};
  cs_int_t *domain = NULL;

  int _nt_cur_abs = -1;
  double _t_cur_abs = 0.;

  const cs_int_t   *var_ptr[1] = {NULL};

  if (cs_glob_n_ranks < 2 || _cs_post_domain == false)
    return;

  dim_ent = fvm_nodal_get_max_entity_dim(exp_mesh);
  n_elts = fvm_nodal_get_n_entities(exp_mesh, dim_ent);

  /* Prepare domain number */

  BFT_MALLOC(domain, n_elts, cs_int_t);

  for (i = 0; i < n_elts; i++)
    domain[i] = cs_glob_mesh->domain_num;

  /* Prepare post-processing output */

  if (sizeof(cs_int_t) == 4)
    datatype = FVM_INT32;
  else if (sizeof(cs_int_t) == 8)
    datatype = FVM_INT64;

  var_ptr[0] = domain;

  if (fvm_writer_get_time_dep(writer) != FVM_WRITER_FIXED_MESH) {
    _nt_cur_abs = nt_cur_abs;
    _t_cur_abs = t_cur_abs;
  }

  fvm_writer_export_field(writer,
                          exp_mesh,
                          _("parallel domain"),
                          FVM_WRITER_PER_ELEMENT,
                          1,
                          FVM_INTERLACE,
                          1,
                          dec_num_parent,
                          datatype,
                          _nt_cur_abs,
                          _t_cur_abs,
                          (const void * *)var_ptr);

  /* Free memory */

  BFT_FREE(domain);
}

/*----------------------------------------------------------------------------
 * Output a post-processing mesh using associated writers.
 *
 * parameters:
 *   post_mesh  <-> pointer to post-processing mesh
 *   nt_cur_abs <-- current time step number
 *   t_cur_abs  <-- current physical time
 *----------------------------------------------------------------------------*/

static void
_cs_post_write_mesh(cs_post_mesh_t  *post_mesh,
                    int              nt_cur_abs,
                    double           t_cur_abs)
{
  int  j;
  cs_bool_t  write_mesh;
  fvm_writer_time_dep_t  time_dep;

  cs_post_writer_t *writer = NULL;

  /* Loop on writers */

  for (j = 0; j < post_mesh->n_writers; j++) {

    writer = _cs_post_writers + post_mesh->writer_id[j];

    time_dep = fvm_writer_get_time_dep(writer->writer);

    write_mesh = false;

    if (time_dep == FVM_WRITER_FIXED_MESH) {
      if (post_mesh->nt_last < 0)
        write_mesh = true;
    }
    else {
      if (post_mesh->nt_last < nt_cur_abs && writer->active == 1)
        write_mesh = true;
    }

    /* Mesh has already been output when associated with writers
       allowing only fixed meshes; for other writers, output it */

    if (write_mesh == true && time_dep != FVM_WRITER_FIXED_MESH) {
      fvm_writer_set_mesh_time(writer->writer, nt_cur_abs, t_cur_abs);
      fvm_writer_export_nodal(writer->writer, post_mesh->exp_mesh);
    }

    if (write_mesh == true && post_mesh->id == -1)
      _cs_post_write_domain(writer->writer,
                            post_mesh->exp_mesh,
                            nt_cur_abs,
                            t_cur_abs);

  }

  if (write_mesh == true)
    post_mesh->nt_last = nt_cur_abs;

  if (   post_mesh->mod_flag_max == FVM_WRITER_FIXED_MESH
      && post_mesh->_exp_mesh != NULL)
    fvm_nodal_reduce(post_mesh->_exp_mesh, 0);
}

/*----------------------------------------------------------------------------
 * Transform an array of flags (markers) to a list
 *
 * parameters:
 *   list_size <-> size of array, then list
 *   list      <-> array of flags, then list
 *
 * returns:
 *   size of list
 *----------------------------------------------------------------------------*/

static cs_int_t
_cs_post_marker_to_list(cs_int_t  list_size,
                        cs_int_t  list[])
{
  cs_int_t  cpt, ind;

  for (cpt = 0, ind = 0; ind < list_size; ind++) {
    if (list[ind] != 0) {
      list[ind] = 0;
      list[cpt++] = ind + 1;
    }
  }

  return cpt;
}

/*----------------------------------------------------------------------------
 * Loop on post-processing meshes to output variables
 *
 * parameters:
 *   nt_cur_abs <-- current time step number
 *   t_cur_abs  <-- current physical time
 *----------------------------------------------------------------------------*/

static void
_cs_post_write_displacements(int     nt_cur_abs,
                             double  t_cur_abs)
{
  int i, j;
  cs_int_t  k, nbr_val;
  fvm_datatype_t datatype;

  fvm_lnum_t   dec_num_parent[1]  = {0};
  cs_post_mesh_t  *post_mesh = NULL;
  cs_post_writer_t  *writer = NULL;
  cs_real_t   *deplacements = NULL;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_real_t   *var_ptr[1] = {NULL};

  /* Loop on writers to check if something must be done */
  /*----------------------------------------------------*/

  if (_cs_post_deformable == false)
    return;

  for (j = 0; j < _cs_post_n_writers; j++) {
    writer = _cs_post_writers + j;
    if (writer->active == 1 && writer->write_displ == true)
      break;
  }
  if (j == _cs_post_n_writers)
    return;

  /* Compute main deformation field */
  /*--------------------------------*/

  nbr_val = mesh->n_vertices * 3;

  BFT_MALLOC(deplacements, nbr_val, cs_real_t);

  assert(mesh->n_vertices == 0 || _cs_post_ini_vtx_coo != NULL);

  for (k = 0; k < nbr_val; k++)
    deplacements[k] = mesh->vtx_coord[k] - _cs_post_ini_vtx_coo[k];

  /* Prepare post-processing */
  /*-------------------------*/

  if (sizeof(cs_real_t) == sizeof(double))
    datatype = FVM_DOUBLE;
  else if (sizeof(cs_real_t) == sizeof(float))
    datatype = FVM_FLOAT;

  var_ptr[0] = deplacements;

  /* Loop on meshes to output displacements */
  /*----------------------------------------*/

  for (i = 0; i < _cs_post_n_meshes; i++) {

    post_mesh = _cs_post_meshes + i;

    for (j = 0; j < post_mesh->n_writers; j++) {

      writer = _cs_post_writers + post_mesh->writer_id[j];

      if (writer->active == 1 && writer->write_displ == true) {

        fvm_writer_export_field(writer->writer,
                                post_mesh->exp_mesh,
                                _("displacement"),
                                FVM_WRITER_PER_NODE,
                                3,
                                FVM_INTERLACE,
                                1,
                                dec_num_parent,
                                datatype,
                                (int)nt_cur_abs,
                                (double)t_cur_abs,
                                (const void * *)var_ptr);

      }

    }

  }

  /* Free memory */

  BFT_FREE(deplacements);
}

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a writer based on Fortran data; this object is based on a choice
 * of a case, directory, and format, as well as indicator for associated
 * mesh's time dependency, and the default output frequency for associated
 * variables.
 *
 * Fortran Interface: use PSTCWR (see cs_post_util.F)
 *
 * SUBROUTINE PSTCW1 (NUMGEP, NOMCAS, NOMREP, NOMFMT, OPTFMT,
 * *****************
 *                    LNMCAS, LNMFMT, LNMREP, LOPFMT,
 *                    INDMOD, NTCHR)
 *
 * INTEGER          NUMWRI      : --> : Number of writer to create (< 0 for
 *                              :     : standard writer, > 0 for user writer)
 * CHARACTER        NOMCAS      : --> : Name of associated case
 * CHARACTER        NOMREP      : --> : Name of associated directory
 * INTEGER          NOMFMT      : --> : Name of associated format
 * INTEGER          OPTFMT      : --> : Additional format options
 * INTEGER          LNMCAS      : --> : Case name length
 * INTEGER          LNMREP      : --> : Directory name length
 * INTEGER          LNMFMT      : --> : Format name length
 * INTEGER          LOPFMT      : --> : Format options string length
 * INTEGER          INDMOD      : --> : 0 if fixed, 1 if deformable,
 *                              :     : 2 if topology changes
 * INTEGER          NTCHR       : --> : Default output frequency
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstcw1, PSTCW1)
(
 const cs_int_t  *numwri,
 const char      *nomcas,
 const char      *nomrep,
 const char      *nomfmt,
 const char      *optfmt,
 const cs_int_t  *lnmcas,
 const cs_int_t  *lnmrep,
 const cs_int_t  *lnmfmt,
 const cs_int_t  *lopfmt,
 const cs_int_t  *indmod,
 const cs_int_t  *ntchr
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
)
{
  /* local variables */

  char  *case_name;
  char  *dir_name;
  char  *nom_format;
  char  *opt_format;

  /* Copy Fortran strings to C strings */

  case_name    = cs_base_string_f_to_c_create(nomcas, *lnmcas);
  dir_name    = cs_base_string_f_to_c_create(nomrep, *lnmrep);
  nom_format = cs_base_string_f_to_c_create(nomfmt, *lnmfmt);
  opt_format = cs_base_string_f_to_c_create(optfmt, *lopfmt);

  /* Main processing */

  cs_post_add_writer(*numwri,
                     case_name,
                     dir_name,
                     nom_format,
                     opt_format,
                     *indmod,
                     *ntchr);

  /* Free temporary C strings */

  cs_base_string_f_to_c_free(&case_name);
  cs_base_string_f_to_c_free(&dir_name);
  cs_base_string_f_to_c_free(&nom_format);
  cs_base_string_f_to_c_free(&opt_format);
}


/*----------------------------------------------------------------------------
 * Create a post-processing mesh; lists of cells or faces to extract are
 * sorted upon exit, whether they were sorted upon calling or not.
 *
 * The list of associated cells is only necessary if the number of cells
 * to extract is strictly greater than 0 and less than the number of cells
 * of the computational mesh.
 *
 * Lists of faces are ignored if the number of extracted cells is nonzero;
 * otherwise, if the number of boundary faces to extract is equal to the
 * number of boundary faces in the computational mesh, and the number of
 * interior faces to extract is zero, than we extrac by default the boundary
 * mesh, and the list of associated boundary faces is thus not necessary.
 *
 * Fortran interface: use PSTCMA (see cs_post_util.F)
 *
 * SUBROUTINE PSTCM1 (NUMMAI, NOMMAI, LNMMAI,
 * *****************
 *                    NBRCEL, NBRFAC, NBRFBR, LSTCEL, LSTFAC, LSTFBR)
 *
 * INTEGER          NUMMAI      : <-- : Number of output mesh to create
 *                              :     : (< 0 for standard mesh,
 *                              :     : > 0 for user mesh)
 * CHARACTER        NOMMAI      : <-- : Name of associated output mesh
 * INTEGER          LNMMAI      : <-- : Mesh name length
 * INTEGER          NBRCEL      : <-- : Number of associated cells
 * INTEGER          NBRFAC      : <-- : Number of associated interior faces
 * INTEGER          NBRFBR      : <-- : Nulber of associated boundary faces
 * INTEGER          LSTCEL      : <-- : List of associated cells
 * INTEGER          LSTFAC      : <-- : List of associated interior faces
 * INTEGER          LSTFBR      : <-- : List of associated boundary faces
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstcm1, PSTCM1)
(
 const cs_int_t  *nummai,
 const char      *nommai,
 const cs_int_t  *lnmmai,
 const cs_int_t  *nbrcel,
 const cs_int_t  *nbrfac,
 const cs_int_t  *nbrfbr,
       cs_int_t   lstcel[],
       cs_int_t   lstfac[],
       cs_int_t   lstfbr[]
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
)
{
  /* local variables */

  char  *mesh_name = NULL;

  /* Copy Fortran strings to C strings */

  mesh_name = cs_base_string_f_to_c_create(nommai, *lnmmai);

  /* Main processing */

  cs_post_add_mesh(*nummai,
                   mesh_name,
                   *nbrcel,
                   *nbrfac,
                   *nbrfbr,
                   lstcel,
                   lstfac,
                   lstfbr);

  /* Free temporary C strings */

  cs_base_string_f_to_c_free(&mesh_name);
}

/*----------------------------------------------------------------------------
 * Create a mesh based upon the extraction of edges from an existing mesh.
 *
 * The newly created edges have no link to their parent elements, so
 * no variable referencing parent elements may be output to this mesh,
 * whose main use is to visualize "true" face edges when polygonal faces
 * are subdivided by the writer. In this way, even highly non-convex
 * faces may be visualized correctly if their edges are overlaid on
 * the surface mesh with subdivided polygons.
 *
 * Fortran interface:
 *
 * SUBROUTINE PSTEDG (NUMMAI, NUMREF)
 * *****************
 *
 * INTEGER          NUMMAI      : <-- : Number of the edges mesh to create
 * INTEGER          NUMREF      : <-- : Number of the existing mesh
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstedg, PSTEDG)
(
 const cs_int_t  *nummai,
 const cs_int_t  *numref
)
{
  cs_post_add_mesh_edges(*nummai, *numref);
}

/*----------------------------------------------------------------------------
 * Assign a category to a post-processing mesh.
 *
 * By default, each mesh is assigned a category id identical to its id.
 * The automatic variables output associated with the main volume and
 * boundary meshes will also be applied to meshes of the same categories
 * (i.e. -1 and -2 respectively, whether meshes -1 and -2 are actually
 * defined or not), so setting a user mesh's category to one of these
 * values will automatically provide the same automatic variable output to
 * the user mesh.
 *
 * Fortran interface:
 *
 * SUBROUTINE PSTCAT (NUMMAI, NUMCAT)
 * *****************
 *
 * INTEGER          NUMMAI      : <-- : Number of the alias to create
 * INTEGER          NUMCAT      : <-- : Number of the assigned category
 *                                      (-1: as volume, -2: as boundary)
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstcat, PSTCAT)
(
 const cs_int_t  *nummai,
 const cs_int_t  *numcat
)
{
  cs_post_set_mesh_category(*nummai, *numcat);
}

/*----------------------------------------------------------------------------
 * Create an alias to a post-processing mesh.
 *
 * Fortran interface:
 *
 * SUBROUTINE PSTALM (NUMMAI, NUMREF)
 * *****************
 *
 * INTEGER          NUMMAI      : <-- : Number of the alias to create
 * INTEGER          NUMREF      : <-- : Number of the associated output mesh
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstalm, PSTALM)
(
 const cs_int_t  *nummai,
 const cs_int_t  *numref
)
{
  cs_post_alias_mesh(*nummai, *numref);
}

/*----------------------------------------------------------------------------
 * Associate a writer to a post-processing mesh.
 *
 * Fortran interface:
 *
 * SUBROUTINE PSTASS (NUMMAI, NUMWRI)
 * *****************
 *
 * INTEGER          NUMMAI      : <-- : Number of the associated output mesh
 * INTEGER          NUMWRI      : <-- : Number of the writer to associate
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstass, PSTASS)
(
 const cs_int_t  *nummai,
 const cs_int_t  *numwri
)
{
  cs_post_associate(*nummai, *numwri);
}

/*----------------------------------------------------------------------------
 * Update the "active" or "inactive" flag for writers based on the current
 * time step and their default output frequency.
 *
 * Fortran interface:
 *
 * SUBROUTINE PSTNTC (NTCABS)
 * *****************
 *
 * INTEGER          NTCABS      : <-- : Current time step number
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstntc, PSTNTC)
(
 const cs_int_t  *ntcabs
)
{
  cs_post_activate_if_default(*ntcabs);
}

/*----------------------------------------------------------------------------
 * Force the "active" or "inactive" flag for a specific writer or for all
 * writers for the current time step.
 *
 * Fortran interface:
 *
 * SUBROUTINE PSTNTC (NUMWRI, INDACT)
 * *****************
 *
 * INTEGER          NUMWRI      : <-- : Writer number, or 0 for all writers
 * INTEGER          INDACT      : <-- : 0 to deactivate, 1 to activate
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstact, PSTACT)
(
 const cs_int_t  *numwri,
 const cs_int_t  *indact
)
{
  cs_post_activate_writer(*numwri, *indact);
}

/*----------------------------------------------------------------------------
 * Output post-processing meshes using associated writers.
 *
 * Fortran interface:
 *
 * SUBROUTINE PSTEMA (NTCABS, TTCABS)
 * *****************
 *
 * INTEGER          NTCABS      : <-- : Current time step number
 * DOUBLE PRECISION TTCABS      : <-- : Current physical time
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstema, PSTEMA)
(
 const cs_int_t   *ntcabs,
 const cs_real_t  *ttcabs
)
{
  cs_post_write_meshes(*ntcabs, *ttcabs);
}

/*----------------------------------------------------------------------------
 * Loop on post-processing meshes to output variables
 *
 * Fortran interface:
 *
 * SUBROUTINE PSTVAR (IDBIA0, IDBRA0,
 * *****************
 *                    NDIM,   NTCABS, NCELET, NCEL,   NFAC,   NFABOR,
 *                    NFML,   NPRFML, NNOD,   LNDFAC, LNDFBR, NCELBR,
 *                    NVAR,   NSCAL,  NPHAS,  NVLSTA, NVISBR,
 *                    NIDEVE, NRDEVE, NITUSE, NRTUSE,
 *                    IFACEL, IFABOR, IFMFBR, IFMCEL, IPRFML,
 *                    IPNFAC, NODFAC, IPNFBR, NODFBR,
 *                    IDEVEL, ITUSER, IA,
 *                    TTCABS, XYZCEN, SURFAC, SURFBO, CDGFAC, CDGFBO,
 *                    XYZNOD, VOLUME,
 *                    DT,     RTPA,   RTP,    PROPCE, PROPFA, PROPFB,
 *                    COEFA,  COEFB,
 *                    STATCE, STATIV, STATFB,
 *                    RDEVEL, RTUSER, RA)
 *
 * INTEGER          IDBIA0      : <-- : Number of first free position in IA
 * INTEGER          IDBRA0      : <-- : Number of first free position in RA
 * INTEGER          NDIM        : <-- : Spatial dimension
 * INTEGER          NTCABS      : --> : Current time step number
 * INTEGER          NCELET      : <-- : Number of extended (real + ghost) cells
 * INTEGER          NFAC        : <-- : Number of interior faces
 * INTEGER          NFABOR      : <-- : Number of boundary faces
 * INTEGER          NFML        : <-- : Number of families (group classes)
 * INTEGER          NPRFML      : <-- : Number of family properties
 * INTEGER          NNOD        : <-- : Number of vertices
 * INTEGER          LNDFAC      : <-- : Size of nodfac
 * INTEGER          LNDFBR      : <-- : Size of nodfbr
 * INTEGER          NCELBR      : <-- : Number of cells on boundary
 * INTEGER          NVAR        : <-- : Number of variables
 * INTEGER          NSCAL       : <-- : Number of scalars
 * INTEGER          NPHAS       : <-- : Number of phases
 * INTEGER          NVLSTA      : <-- : Number of statistical variables (lagr)
 * INTEGER          NVISBR      : <-- : Number of boundary stat. variables (lagr)
 * INTEGER          NIDEVE      : <-- : Size of IDEVEL integer array
 * INTEGER          NRDEVE      : <-- : Size of RDEVEL floating-point array
 * INTEGER          NITUSE      : <-- : Size of ITUSER integer array
 * INTEGER          NRTUSE      : <-- : Size of RTUSER floating-point array
 * INTEGER          IFACEL      : <-- : Interior faces -> cells connectivity
 * INTEGER          IFABOR      : <-- : Boundary faces -> cell connectivity
 * INTEGER          IFMFBR      : <-- : Boundary face families
 * INTEGER          IFMCEL      : <-- : Cell families
 * INTEGER          IPRFML      : <-- : List of family properties
 * INTEGER          IPNFAC      : <-- : Interior faces -> vertices connect. idx.
 * INTEGER          NODFAC      : <-- : Interior faces -> vertices connectivity
 * INTEGER          IPNFBR      : <-- : Boundary faces -> vertices connect. idx.
 * INTEGER          NODFBR      : <-- : Boundary faces -> vertices connectivity
 * INTEGER          IDEVEL      : <-- : IDEVEL integer array
 * INTEGER          ITUSER      : <-- : ITUSER integer array
 * INTEGER          IA          : <-- : IA integer array
 * DOUBLE PRECISION TTCABS      : <-- : Current physical time
 * DOUBLE PRECISION XYZCEN      : <-- : Points associated with cell centers
 * DOUBLE PRECISION SURFAC      : <-- : Interior face surface vectors
 * DOUBLE PRECISION SURFBO      : <-- : Boundary face surface vectors
 * DOUBLE PRECISION CDGFAC      : <-- : Interior face centers
 * DOUBLE PRECISION CDGFBO      : <-- : Boundary face vectors
 * DOUBLE PRECISION XYZNOD      : <-- : Vertex coordinates (optional)
 * DOUBLE PRECISION VOLUME      : <-- : Cell volumes
 * DOUBLE PRECISION DT          : <-- : Local time step
 * DOUBLE PRECISION RTPA        : <-- : Cell variables at previous time step
 * DOUBLE PRECISION RTP         : <-- : Cell variables
 * DOUBLE PRECISION PROPCE      : <-- : Cell physical properties
 * DOUBLE PRECISION PROPFA      : <-- : Interior face physical properties
 * DOUBLE PRECISION PROPFB      : <-- : Boundary face physical properties
 * DOUBLE PRECISION COEFA       : <-- : Boundary conditions array
 * DOUBLE PRECISION COEFB       : <-- : Boundary conditions array
 * DOUBLE PRECISION STATCE      : <-- : Cell statistics (Lagrangian)
 * DOUBLE PRECISION STATIV      : <-- : Cell variance statistics (Lagrangian)
 * DOUBLE PRECISION STATFB      : <-- : Boundary face statistics (Lagrangian)
 * DOUBLE PRECISION RDEVEL      : <-- : RDEVEL floating-point array
 * DOUBLE PRECISION RTUSER      : <-- : RTUSER floating-point array
 * DOUBLE PRECISION RA          : <-- : RA floating-point array
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstvar, PSTVAR)
(
 const cs_int_t   *idbia0,
 const cs_int_t   *idbra0,
 const cs_int_t   *ndim,
 const cs_int_t   *ntcabs,
 const cs_int_t   *ncelet,
 const cs_int_t   *ncel,
 const cs_int_t   *nfac,
 const cs_int_t   *nfabor,
 const cs_int_t   *nfml,
 const cs_int_t   *nprfml,
 const cs_int_t   *nnod,
 const cs_int_t   *lndfac,
 const cs_int_t   *lndfbr,
 const cs_int_t   *ncelbr,
 const cs_int_t   *nvar,
 const cs_int_t   *nscal,
 const cs_int_t   *nphas,
 const cs_int_t   *nvlsta,
 const cs_int_t   *nvisbr,
 const cs_int_t   *nideve,
 const cs_int_t   *nrdeve,
 const cs_int_t   *nituse,
 const cs_int_t   *nrtuse,
 const cs_int_t    ifacel[],
 const cs_int_t    ifabor[],
 const cs_int_t    ifmfbr[],
 const cs_int_t    ifmcel[],
 const cs_int_t    iprfml[],
 const cs_int_t    ipnfac[],
 const cs_int_t    nodfac[],
 const cs_int_t    ipnfbr[],
 const cs_int_t    nodfbr[],
 const cs_int_t    idevel[],
       cs_int_t    ituser[],
       cs_int_t    ia[],
 const cs_real_t  *ttcabs,
 const cs_real_t   xyzcen[],
 const cs_real_t   surfac[],
 const cs_real_t   surfbo[],
 const cs_real_t   cdgfac[],
 const cs_real_t   cdgfbo[],
 const cs_real_t   xyznod[],
 const cs_real_t   volume[],
 const cs_real_t   dt[],
 const cs_real_t   rtpa[],
 const cs_real_t   rtp[],
 const cs_real_t   propce[],
 const cs_real_t   propfa[],
 const cs_real_t   propfb[],
 const cs_real_t   coefa[],
 const cs_real_t   coefb[],
 const cs_real_t   statce[],
 const cs_real_t   stativ[],
 const cs_real_t   statfb[],
 const cs_real_t   rdevel[],
       cs_real_t   rtuser[],
       cs_real_t   ra[]
)
{
  /* local variables */

  int i, j, k;
  int dim_ent;
  cs_int_t   itypps[3];
  cs_int_t   ind_cel, ind_fac, dec_num_fbr;
  cs_int_t   n_elts, n_elts_max;
  cs_int_t   nummai, numtyp, imodif;

  cs_bool_t  active;

  cs_post_mesh_t  *post_mesh;
  cs_post_writer_t  *writer;

  cs_int_t    n_cells, n_i_faces, n_b_faces;
  cs_int_t    *cell_list, *i_face_list, *b_face_list;

  cs_int_t    *num_ent_parent = NULL;
  cs_real_t   *var_trav = NULL;
  cs_real_t   *cel_vals = NULL;
  cs_real_t   *i_face_vals = NULL;
  cs_real_t   *b_face_vals = NULL;

  /* Loop on writers to check if something must be done */
  /*----------------------------------------------------*/

  for (j = 0; j < _cs_post_n_writers; j++) {
    writer = _cs_post_writers + j;
    if (writer->active == 1)
      break;
  }
  if (j == _cs_post_n_writers)
    return;

  /* Possible modification of post-processing meshes */
  /*-------------------------------------------------*/

  n_elts_max = 0;

  for (i = 0; i < _cs_post_n_meshes; i++) {

    post_mesh = _cs_post_meshes + i;

    active = false;

    for (j = 0; j < post_mesh->n_writers; j++) {
      writer = _cs_post_writers + post_mesh->writer_id[j];
      if (writer->active == 1)
        active = true;
    }

    /* Modifiable user mesh, not an alias, active at this time step */

    if (   active == true
        && post_mesh->alias < 0
        && post_mesh->id > 0
        && post_mesh->mod_flag_min == FVM_WRITER_TRANSIENT_CONNECT) {

      const fvm_nodal_t * exp_mesh = post_mesh->exp_mesh;

      dim_ent = fvm_nodal_get_max_entity_dim(exp_mesh);
      n_elts = fvm_nodal_get_n_entities(exp_mesh, dim_ent);

      if (n_elts > n_elts_max) {
        n_elts_max = n_elts;
        BFT_REALLOC(num_ent_parent, n_elts_max, cs_int_t);
      }

      nummai = post_mesh->id;
      numtyp = post_mesh->cat_id;

      /* Get corresponding entity lists */

      fvm_nodal_get_parent_num(exp_mesh, dim_ent, num_ent_parent);

      for (k = 0; k < 3; k++)
        itypps[k] = post_mesh->ent_flag[k];

      /* Oversize lists, as the user may fill an arbitrary part of those */

      BFT_MALLOC(cell_list, cs_glob_mesh->n_cells, cs_int_t);
      BFT_MALLOC(i_face_list, cs_glob_mesh->n_i_faces, cs_int_t);
      BFT_MALLOC(b_face_list, cs_glob_mesh->n_b_faces, cs_int_t);

      n_cells = 0;
      n_i_faces = 0;
      n_b_faces = 0;

      /* Set lists to zero */

      if (dim_ent == 3)
        for (ind_cel = 0; ind_cel < cs_glob_mesh->n_cells; ind_cel++)
          cell_list[ind_cel] = 0;
      else if (dim_ent == 2) {
        for (ind_fac = 0; ind_fac < cs_glob_mesh->n_b_faces; ind_fac++)
          b_face_list[ind_fac] = 0;
        for (ind_fac = 0; ind_fac < cs_glob_mesh->n_i_faces; ind_fac++)
          i_face_list[ind_fac] = 0;
      }

      /* If the elements of the FVM mesh are divided, a same parent number
         may appear several times; we thus use a marker logic. */

      if (dim_ent == 3) {
        for (ind_cel = 0; ind_cel < n_elts; ind_cel++)
          cell_list[num_ent_parent[ind_cel] - 1] = 1;
      }

      /* For faces, the number of interior "parent" faces known by FVM
         are shifted by the total number of boundary faces
         (c.f. construction in cs_mesh_connect...()) */

      else if (dim_ent == 2) {
        dec_num_fbr = cs_glob_mesh->n_b_faces;
        for (ind_fac = 0; ind_fac < n_elts; ind_fac++) {
          if (num_ent_parent[ind_fac] > dec_num_fbr)
            i_face_list[num_ent_parent[ind_fac] - dec_num_fbr - 1] = 1;
          else
            b_face_list[num_ent_parent[ind_fac] - 1] = 1;
        }
      }

      /* Transform markers to lists */

      if (dim_ent == 3) {
        n_cells = _cs_post_marker_to_list(cs_glob_mesh->n_cells,
                                          cell_list);
      }
      else if (dim_ent == 2) {
        n_i_faces = _cs_post_marker_to_list(cs_glob_mesh->n_i_faces,
                                            i_face_list);
        n_b_faces = _cs_post_marker_to_list(cs_glob_mesh->n_b_faces,
                                            b_face_list);
      }

      /* User modification of the mesh definition */

      imodif = 0;

      CS_PROCF(usmpst, USMPST) (idbia0, idbra0, &nummai,
                                ndim, ncelet, ncel, nfac, nfabor, nfml, nprfml,
                                nnod, lndfac, lndfbr, ncelbr,
                                nvar, nscal, nphas, nvlsta,
                                &n_cells, &n_i_faces, &n_b_faces,
                                nideve, nrdeve, nituse, nrtuse, &imodif,
                                itypps, ifacel, ifabor, ifmfbr, ifmcel, iprfml,
                                ipnfac, nodfac, ipnfbr, nodfbr,
                                cell_list, i_face_list, b_face_list,
                                idevel, ituser, ia,
                                xyzcen, surfac, surfbo, cdgfac, cdgfbo, xyznod,
                                volume, dt, rtpa, rtp, propce, propfa, propfb,
                                coefa, coefb, statce,
                                cel_vals, i_face_vals, b_face_vals,
                                rdevel, rtuser, ra);

      if (imodif > 0)
        cs_post_modify_mesh(post_mesh->id,
                            n_cells,
                            n_i_faces,
                            n_b_faces,
                            cell_list,
                            i_face_list,
                            b_face_list);

      BFT_FREE(cell_list);
      BFT_FREE(i_face_list);
      BFT_FREE(b_face_list);

    }

  }

  /* We now make sure aliases are synchronized */

  for (i = 0; i < _cs_post_n_meshes; i++) {

    post_mesh = _cs_post_meshes + i;

    if (post_mesh->alias > -1) {

      const cs_post_mesh_t  *ref_mesh;

      ref_mesh = _cs_post_meshes + post_mesh->alias;

      for (j = 0; j < 3; j++)
        post_mesh->ent_flag[j] = ref_mesh->ent_flag[j];

      post_mesh->n_i_faces = ref_mesh->n_i_faces;
      post_mesh->n_b_faces = ref_mesh->n_b_faces;

    }

  }

  /* Output of meshes or vertex displacement field if necessary */
  /*------------------------------------------------------------*/

  cs_post_write_meshes(*ntcabs, *ttcabs);

  if (_cs_post_deformable == true)
    _cs_post_write_displacements(*ntcabs, *ttcabs);

  /* Output of variables by registered function instances */
  /*------------------------------------------------------*/

  for (i = 0; i < _cs_post_nbr_var_tp; i++) {
    _cs_post_f_var_tp[i](_cs_post_i_var_tp[i],
                         *ntcabs,
                         *ttcabs);
  }

  /* Output of variables associated with post-processing meshes */
  /*------------------------------------------------------------*/

  /* n_elts_max already initialized before and during the
     eventual modification of post-processing mesh definitions,
     and num_ent_parent allocated if n_elts_max > 0 */

  BFT_MALLOC(var_trav, n_elts_max * 3, cs_real_t);

  /* Main loop on post-processing meshes */

  for (i = 0; i < _cs_post_n_meshes; i++) {

    post_mesh = _cs_post_meshes + i;

    active = false;

    for (j = 0; j < post_mesh->n_writers; j++) {
      writer = _cs_post_writers + post_mesh->writer_id[j];
      if (writer->active == 1)
        active = true;
    }

    /* If the mesh is active at this time step */
    /*-----------------------------------------*/

    if (active == true) {

      const fvm_nodal_t * exp_mesh = post_mesh->exp_mesh;

      dim_ent = fvm_nodal_get_max_entity_dim(exp_mesh);
      n_elts = fvm_nodal_get_n_entities(exp_mesh, dim_ent);

      if (n_elts > n_elts_max) {
        n_elts_max = n_elts;
        BFT_REALLOC(var_trav, n_elts_max * 3, cs_real_t);
        BFT_REALLOC(num_ent_parent, n_elts_max, cs_int_t);
      }

      nummai = post_mesh->id;
      numtyp = post_mesh->cat_id;

      /* Get corresponding element lists */

      fvm_nodal_get_parent_num(exp_mesh, dim_ent, num_ent_parent);

      for (k = 0; k < 3; k++)
        itypps[k] = post_mesh->ent_flag[k];

      /* We can output variables for this time step */
      /*--------------------------------------------*/

      n_cells = 0;
      n_i_faces = 0;
      n_b_faces = 0;
      cell_list = NULL;
      i_face_list = NULL;
      b_face_list = NULL;

      /* Here list sizes are adjusted, and we point to the array filled
         by fvm_nodal_get_parent_num() if possible. */

      if (dim_ent == 3) {
        n_cells = n_elts;
        cell_list = num_ent_parent;
      }

      /* The numbers of "parent" interior faces known by FVM
         are shifted by the total number of boundary faces */

      else if (dim_ent == 2 && n_elts > 0) {

        dec_num_fbr = cs_glob_mesh->n_b_faces;

        for (ind_fac = 0; ind_fac < n_elts; ind_fac++) {
          if (num_ent_parent[ind_fac] > dec_num_fbr)
            n_i_faces++;
          else
            n_b_faces++;
        }

        /* boundary faces only: parent FVM face numbers unchanged */
        if (n_i_faces == 0) {
          b_face_list = num_ent_parent;
        }

        /* interior faces only: parents FVM face numbers shifted */
        else if (n_b_faces == 0) {
          for (ind_fac = 0; ind_fac < n_elts; ind_fac++)
            num_ent_parent[ind_fac] -= dec_num_fbr;
          i_face_list = num_ent_parent;
        }

        /* interior and boundary faces: numbers must be separated */

        else {

          BFT_MALLOC(i_face_list, n_i_faces, cs_int_t);
          BFT_MALLOC(b_face_list, n_b_faces, cs_int_t);

          n_i_faces = 0, n_b_faces = 0;

          for (ind_fac = 0; ind_fac < n_elts; ind_fac++) {
            if (num_ent_parent[ind_fac] > dec_num_fbr)
              i_face_list[n_i_faces++] = num_ent_parent[ind_fac] - dec_num_fbr;
            else
              b_face_list[n_b_faces++] = num_ent_parent[ind_fac];
          }

        }

        /* In all cases, update the number of interior and boundary faces
           (useful in case of splitting of FVM mesh elements) for functions
           called by this one */

        post_mesh->n_i_faces = n_i_faces;
        post_mesh->n_b_faces = n_b_faces;

      }

      /* Pointers to variable assembly arrays, set to NULL if unused
         (so as to provoke an immediate error in case of incorrect use) */

      cel_vals = var_trav;
      i_face_vals = cel_vals + (n_cells * 3);
      b_face_vals = i_face_vals + (n_i_faces * 3);

      if (n_cells == 0)
        cel_vals = NULL;
      if (n_i_faces == 0)
        i_face_vals = NULL;
      if (n_b_faces == 0)
        b_face_vals = NULL;

      /* Standard post-processing */

      if (numtyp < 0)
        CS_PROCF(dvvpst, DVVPST) (idbia0, idbra0, &nummai, &numtyp,
                                  ndim, ncelet, ncel, nfac, nfabor, nfml, nprfml,
                                  nnod, lndfac, lndfbr, ncelbr,
                                  nvar, nscal, nphas, nvlsta, nvisbr,
                                  &n_cells, &n_i_faces, &n_b_faces,
                                  nideve, nrdeve, nituse, nrtuse,
                                  itypps, ifacel, ifabor, ifmfbr, ifmcel, iprfml,
                                  ipnfac, nodfac, ipnfbr, nodfbr,
                                  cell_list, i_face_list, b_face_list,
                                  idevel, ituser, ia,
                                  xyzcen, surfac, surfbo, cdgfac, cdgfbo, xyznod,
                                  volume, dt, rtpa, rtp, propce, propfa, propfb,
                                  coefa, coefb, statce, stativ , statfb ,
                                  cel_vals, i_face_vals, b_face_vals,
                                  rdevel, rtuser, ra);

      /* Call to user subroutine for additional post-processing */

      CS_PROCF(usvpst, USVPST) (idbia0, idbra0, &nummai,
                                ndim, ncelet, ncel, nfac, nfabor, nfml, nprfml,
                                nnod, lndfac, lndfbr, ncelbr,
                                nvar, nscal, nphas, nvlsta,
                                &n_cells, &n_i_faces, &n_b_faces,
                                nideve, nrdeve, nituse, nrtuse,
                                itypps, ifacel, ifabor, ifmfbr, ifmcel, iprfml,
                                ipnfac, nodfac, ipnfbr, nodfbr,
                                cell_list, i_face_list, b_face_list,
                                idevel, ituser, ia,
                                xyzcen, surfac, surfbo, cdgfac, cdgfbo, xyznod,
                                volume, dt, rtpa, rtp, propce, propfa, propfb,
                                coefa, coefb, statce,
                                cel_vals, i_face_vals, b_face_vals,
                                rdevel, rtuser, ra);

      /* In case of mixed interior and boundary faces, free
         additional arrays */

      if (i_face_list != NULL && b_face_list != NULL) {
        BFT_FREE(i_face_list);
        BFT_FREE(b_face_list);
      }

    }

  }

  /* Free memory */

  BFT_FREE(num_ent_parent);
  BFT_FREE(var_trav);
}

/*----------------------------------------------------------------------------
 * Post-processing output of a variable defined on cells or faces of a mesh
 * using associated writers.
 *
 * Fortran interface; use PSTEVA (see cs_post_util.F)
 *
 * SUBROUTINE PSTEV1 (NUMMAI, NOMVAR, LNMVAR, IDIMT,  IENTLA, IVARPR,
 * *****************
 *                    NTCABS, TTCABS, VARCEL, VARFAC, VARFBR)
 *
 * INTEGER          NUMMAI      : <-- : Number of associated output mesh
 * CHARACTER        NOMVAR      : <-- : Name of associated variable
 * INTEGER          LNMVAR      : <-- : Variable name length
 * INTEGER          IDIMT       : <-- : 1 for scalar, 3 for vector
 * INTEGER          IENTLA      : <-- : If a vector, 1 for interlaced values
 *                              :     : (x1, y1, z1, x2, y2, ..., yn, zn),
 *                              :     : 0 otherwise (x1, x2, ...xn, y1, y2, ...)
 * INTEGER          IVARPR      : <-- : 1 if variable is defined on "parent"
 *                              :     : mesh, 2 if defined on output mesh
 * INTEGER          NTCABS      : <-- : Current time step number
 * DOUBLE PRECISION TTCABS      : <-- : Current physical time
 * DOUBLE PRECISION VARCEL(*)   : <-- : Cell values
 * DOUBLE PRECISION VARFAC(*)   : <-- : Interior face values
 * DOUBLE PRECISION VARFBO(*)   : <-- : Boundary face values
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstev1, PSTEV1)
(
 const cs_int_t   *nummai,
 const char       *nomvar,
 const cs_int_t   *lnmvar,
 const cs_int_t   *idimt,
 const cs_int_t   *ientla,
 const cs_int_t   *ivarpr,
 const cs_int_t   *ntcabs,
 const cs_real_t  *ttcabs,
 const cs_real_t   varcel[],
 const cs_real_t   varfac[],
 const cs_real_t   varfbr[]
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
)
{
  cs_bool_t  use_parent = false;
  cs_bool_t  interlace = false;

  char  *var_name = NULL;

  if (*ivarpr == 1)
    use_parent = true;
  else if (*ivarpr == 0)
    use_parent = false;
  else
    bft_error(__FILE__, __LINE__, 0,
              _("The PSTEVA sub-routine argument IVARPR must be\n"
                "equal to 0 or 1, and not %d.\n"), (int)(*ivarpr));

  if (*ientla == 0)
    interlace = false;
  else if (*ientla == 1)
    interlace = true;
  else
    bft_error(__FILE__, __LINE__, 0,
              _("The PSTEVA sub-routine argument IENTLA must be\n"
                "equal to 0 or 1, and not %d.\n"), (int)(*ientla));


  /* Copy Fortran strings to C strings */

  var_name = cs_base_string_f_to_c_create(nomvar, *lnmvar);

  /* Main processing */

  cs_post_write_var(*nummai,
                    var_name,
                    *idimt,
                    interlace,
                    use_parent,
                    CS_POST_TYPE_cs_real_t,
                    *ntcabs,
                    *ttcabs,
                    varcel,
                    varfac,
                    varfbr);

  /* Free temporary C strings */

  cs_base_string_f_to_c_free(&var_name);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a writer; this objects manages a case's name, directory, and format,
 * as well as associated mesh's time dependency, and the default output
 * frequency for associated variables.
 *
 * parameters:
 *   writer_id <-- number of writer to create (< 0 reserved, > 0 for user)
 *   case_name <-- associated case name
 *   dir_name  <-- associated directory name
 *   fmt_name  <-- associated format name
 *   fmt_opts  <-- associated format options
 *   mod_flag  <-- 0 if fixed, 1 if deformable, 2 if topolygy changes,
 *                 +10 add a displacement field
 *   frequency <-- default output frequency
 *----------------------------------------------------------------------------*/

void
cs_post_add_writer(int          writer_id,
                   const char  *case_name,
                   const char  *dir_name,
                   const char  *fmt_name,
                   const char  *fmt_opts,
                   cs_int_t     mod_flag,
                   cs_int_t     frequency)
{
  /* local variables */

  int    i;

  cs_post_writer_t  *writer = NULL;
  fvm_writer_time_dep_t  time_dep = FVM_WRITER_FIXED_MESH;

  /* Check that the required mesh is available */

  if (writer_id == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("The requested post-processing writer number\n"
                "must be < 0 (reserved) or > 0 (user).\n"));

  for (i = 0; i < _cs_post_n_writers; i++) {
    if ((_cs_post_writers + i)->id == writer_id)
      bft_error(__FILE__, __LINE__, 0,
                _("The requested post-processing writer number\n"
                  "(%d) has already been assigned.\n"), (int)writer_id);
  }

  /* Resize global writers array */

  if (_cs_post_n_writers == _cs_post_n_writers_max) {

    if (_cs_post_n_writers_max == 0)
      _cs_post_n_writers_max = 4;
    else
      _cs_post_n_writers_max *= 2;

    BFT_REALLOC(_cs_post_writers,
                _cs_post_n_writers_max,
                cs_post_writer_t);

  }

  if (writer_id < _cs_post_min_writer_id)
    _cs_post_min_writer_id = writer_id;

  _cs_post_n_writers += 1;

  /* Assign newly created writer to the structure */

  writer = _cs_post_writers + _cs_post_n_writers - 1;

  writer->id = writer_id;
  writer->frequency = frequency;
  writer->write_displ = false;
  writer->active = 0;

  if (mod_flag >= 10) {
    writer->write_displ = true;
    mod_flag -= 10;
  }

  if (mod_flag == 1)
    time_dep = FVM_WRITER_TRANSIENT_COORDS;
  else if (mod_flag >= 2)
    time_dep = FVM_WRITER_TRANSIENT_CONNECT;

  writer->writer = fvm_writer_init(case_name,
                                   dir_name,
                                   fmt_name,
                                   fmt_opts,
                                   time_dep);
}

/*----------------------------------------------------------------------------
 * Create a post-processing mesh; lists of cells or faces to extract are
 * sorted upon exit, whether they were sorted upon calling or not.
 *
 * The list of associated cells is only necessary if the number of cells
 * to extract is strictly greater than 0 and less than the number of cells
 * of the computational mesh.
 *
 * Lists of faces are ignored if the number of extracted cells is nonzero;
 * otherwise, if the number of boundary faces to extract is equal to the
 * number of boundary faces in the computational mesh, and the number of
 * interior faces to extract is zero, than we extrac by default the boundary
 * mesh, and the list of associated boundary faces is thus not necessary.
 *
 * parameters:
 *   mesh_id     <-- number of mesh to create (< 0 reserved, > 0 for user)
 *   mesh_name   <-- associated mesh name
 *   n_cells     <-- number of associated cells
 *   n_i_faces   <-- number of associated interior faces
 *   n_b_faces   <-- number of associated boundary faces
 *   cell_list   <-- list of associated cells
 *   i_face_list <-- list of associated interior faces
 *   b_face_list <-- list of associated boundary faces
 *----------------------------------------------------------------------------*/

void
cs_post_add_mesh(int          mesh_id,
                 const char  *mesh_name,
                 cs_int_t     n_cells,
                 cs_int_t     n_i_faces,
                 cs_int_t     n_b_faces,
                 cs_int_t     cell_list[],
                 cs_int_t     i_face_list[],
                 cs_int_t     b_face_list[])
{
  /* local variables */

  cs_post_mesh_t  *post_mesh = NULL;

  /* Add and initialize base structure */

  post_mesh = _cs_post_add_mesh(mesh_id);

  /* Create mesh and assign to structure */

  _cs_post_define_mesh(post_mesh,
                       mesh_name,
                       n_cells,
                       n_i_faces,
                       n_b_faces,
                       cell_list,
                       i_face_list,
                       b_face_list);
}

/*----------------------------------------------------------------------------
 * Create a post-processing mesh associated with an existing exportable mesh
 * representation.
 *
 * If the exportable mesh is not intended to be used elsewhere, one can choose
 * to transfer its property to the post-processing mesh, which will then
 * manage its lifecycle based on its own requirements.
 *
 * If the exportable mesh must still be shared, one must be careful to
 * maintain consistency between this mesh and the post-processing output.
 *
 * The mesh in exportable dimension may be of a lower dimension than
 * its parent mesh, if it has been projected. In this case, a
 * dim_shift value of 1 indicates that parent cells are mapped to
 * exportable faces, and faces to edges, while a dim_shift value of 2
 * would indicate that parent cells are mapped to edges.
 * This is important when variables values are exported.
 *
 * parameters:
 *   mesh_id   <-- number of mesh to create (< 0 reserved, > 0 for user)
 *   exp_mesh  <-- mesh in exportable representation (i.e. fvm_nodal_t)
 *   dim_shift <-- nonzero if exp_mesh has been projected
 *   transfer  <-- if true, ownership of exp_mesh is transferred to the
 *                 post-processing mesh
 *----------------------------------------------------------------------------*/

void
cs_post_add_existing_mesh(int           mesh_id,
                          fvm_nodal_t  *exp_mesh,
                          int           dim_shift,
                          cs_bool_t     transfer)
{
  /* local variables */

  int       i;
  int       glob_flag[3];
  cs_int_t  dec_num_fbr, ind_fac;

  int    loc_flag[3] = {1, 1, 1};  /* Flags 0 to 2 "inverted" compared
                                      to others so as to use a single
                                      call to
                                      MPI_Allreduce(..., MPI_MIN, ...) */

  int         dim_ent = 0;
  cs_bool_t   maj_ent_flag = false;
  fvm_lnum_t  n_elts = 0;

  fvm_lnum_t      *num_ent_parent = NULL;
  cs_post_mesh_t  *post_mesh = NULL;

  /* Initialization of base structure */

  post_mesh = _cs_post_add_mesh(mesh_id);

  /* Assign mesh to structure */

  post_mesh->exp_mesh = exp_mesh;

  if (transfer == true)
    post_mesh->_exp_mesh = exp_mesh;

  /* Compute number of cells and/or faces */

  dim_ent = fvm_nodal_get_max_entity_dim(exp_mesh) + dim_shift;
  n_elts = fvm_nodal_get_n_entities(exp_mesh, dim_ent);

  if (dim_ent == 3 && n_elts > 0)
    loc_flag[0] = 0;

  else if (dim_ent == 2 && n_elts > 0) {

    BFT_MALLOC(num_ent_parent, n_elts, cs_int_t);

    fvm_nodal_get_parent_num(exp_mesh, dim_ent, num_ent_parent);

    dec_num_fbr = cs_glob_mesh->n_b_faces;
    for (ind_fac = 0; ind_fac < n_elts; ind_fac++) {
      if (num_ent_parent[ind_fac] > dec_num_fbr)
        post_mesh->n_i_faces += 1;
      else
        post_mesh->n_b_faces += 1;
    }

    BFT_FREE(num_ent_parent);

    if (post_mesh->n_i_faces > 0)
      loc_flag[1] = 0;
    else if (post_mesh->n_b_faces > 0)
      loc_flag[2] = 0;

  }

  for (i = 0; i < 3; i++)
    glob_flag[i] = loc_flag[i];

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)
    MPI_Allreduce (loc_flag, glob_flag, 3, MPI_INT, MPI_MIN,
                   cs_glob_mpi_comm);
#endif

  /* Global indicators of mesh entity type presence;
     updated only if the mesh is not totally empty (for time-depending
     meshes, empty at certain times, we want to know the last type
     of entity used in USMPST) */

  for (i = 0; i < 3; i++) {
    if (glob_flag[i] == 0)
      maj_ent_flag = true;
  }

  if (maj_ent_flag == true) {
    for (i = 0; i < 3; i++) {
      if (glob_flag[i] == 0)         /* Inverted glob_flag 0 to 2 logic */
        post_mesh->ent_flag[i] = 1;  /* (c.f. remark above) */
      else
        post_mesh->ent_flag[i] = 0;
    }
  }
}

/*----------------------------------------------------------------------------
 * Create a mesh based upon the extraction of edges from an existing mesh.
 *
 * The newly created edges have no link to their parent elements, so
 * no variable referencing parent elements may be output to this mesh,
 * whose main use is to visualize "true" face edges when polygonal faces
 * are subdivided by the writer. In this way, even highly non-convex
 * faces may be visualized correctly if their edges are overlaid on
 * the surface mesh with subdivided polygons.
 *
 * parameters:
 *   edges_id <-- id of edges mesh to create (< 0 reserved, > 0 for user)
 *   base_id  <-- id of existing mesh (< 0 reserved, > 0 for user)
 *----------------------------------------------------------------------------*/

void
cs_post_add_mesh_edges(int  edges_id,
                       int  base_id)
{
  /* local variables */

  char *edges_name = NULL;
  cs_post_mesh_t *post_edges = NULL;
  fvm_nodal_t *exp_edges = NULL;

  const cs_post_mesh_t *post_base = _cs_post_meshes +_cs_post_mesh_id(base_id);
  const fvm_nodal_t *exp_mesh = post_base->exp_mesh;
  const char *exp_name = fvm_nodal_get_name(exp_mesh);

  /* Add and initialize base structure */

  post_edges = _cs_post_add_mesh(edges_id);

  /* Copy mesh edges to new mesh structure */

  BFT_MALLOC(edges_name, strlen(exp_name) + strlen(_(" edges")) + 1, char);

  strcpy(edges_name, exp_name);
  strcat(edges_name, _(" edges"));

  exp_edges = fvm_nodal_copy_edges(edges_name, exp_mesh);

  BFT_FREE(edges_name);

  /* Create mesh and assign to structure */

  post_edges->exp_mesh = exp_edges;
  post_edges->_exp_mesh = exp_edges;
}

/*----------------------------------------------------------------------------
 * Assign a category to a post-processing mesh.
 *
 * By default, each mesh is assigned a category id identical to its id.
 * The automatic variables output associated with the main volume and
 * boundary meshes will also be applied to meshes of the same categories
 * (i.e. -1 and -2 respectively, whether meshes -1 and -2 are actually
 * defined or not), so setting a user mesh's category to one of these
 * values will automatically provide the same automatic variable output to
 * the user mesh.
 *
 * parameters:
 *   mesh_id     <-- id of associated mesh
 *   category_id <-- id of mesh category (-1: as volume, -2: as boundary)
 *----------------------------------------------------------------------------*/

void
cs_post_set_mesh_category(int  mesh_id,
                          int  category_id)
{
  /* local variables */

  int _mesh_id;
  cs_post_mesh_t  *post_mesh = NULL;

  /* Get base structure (return if we do not own the mesh) */

  _mesh_id = _cs_post_mesh_id(mesh_id);
  post_mesh = _cs_post_meshes + _mesh_id;

  /* Set category */

  post_mesh->cat_id = category_id;
}

/*----------------------------------------------------------------------------
 * Create an alias to a post-processing mesh.
 *
 * An alias allows association of an extra identifier (number) to an
 * existing post-processing mesh, and thus to associate different writers
 * than those associated with the existing mesh. For example, this allows
 * outputting a set of main variables every n1 time steps with one writer,
 * and outputting a specific set of variables every n2 time time steps to
 * another post-processing set using another writer, without the overhead
 * that would be incurred by duplication of the post-processing mesh.
 *
 * An alias is thus treated in all points like its associated mesh;
 * if the definition of either one is modified, that of the other is
 * modified also.
 *
 * It is forbidden to associate an alias to another alias (as there is no
 * identified use for this, and it would make consistency checking more
 * difficult), but multiple aliases may be associated with a given mesh.
 *
 * parameters:
 *   alias_id <-- number of alias to create (< 0 reserved, > 0 for user)
 *   mesh_id  <-- id of associated mesh
 *----------------------------------------------------------------------------*/

void
cs_post_alias_mesh(int  alias_id,
                   int  mesh_id)
{
  /* local variables */

  int    indref, j;

  cs_post_mesh_t  *post_mesh = NULL;
  cs_post_mesh_t  *ref_mesh = NULL;

  /* Initial checks */

  indref = _cs_post_mesh_id(mesh_id);
  ref_mesh = _cs_post_meshes + indref;

  if (ref_mesh->alias > -1)
    bft_error(__FILE__, __LINE__, 0,
              _("The mesh %d cannot be an alias of mesh %d,\n"
                "which is itself an alias of mesh %d.\n"),
              (int)alias_id, (int)mesh_id,
              (int)((_cs_post_meshes + ref_mesh->alias)->id));

  /* Initialize base structure */

  post_mesh = _cs_post_add_mesh(alias_id);

  /* Update ref_mesh, as the adress of _cs_post_meshes may
     have been changed by reallocation */

  ref_mesh = _cs_post_meshes + indref;

  /* Links to the reference mesh */

  post_mesh->alias = indref;

  post_mesh->exp_mesh = ref_mesh->exp_mesh;

  post_mesh->mod_flag_min = ref_mesh->mod_flag_min;
  post_mesh->mod_flag_max = ref_mesh->mod_flag_max;

  for (j = 0; j < 3; j++)
    post_mesh->ent_flag[j] = ref_mesh->ent_flag[j];

  post_mesh->n_i_faces = ref_mesh->n_i_faces;
  post_mesh->n_b_faces = ref_mesh->n_b_faces;
}

/*----------------------------------------------------------------------------
 * Check for the existence of a writer of the given id.
 *
 * parameters:
 *   writer_id <-- writer id to check
 *
 * returns:
 *   true if writer with this id exists, false otherwise
 *----------------------------------------------------------------------------*/

cs_bool_t
cs_post_writer_exists(int  writer_id)
{
  /* local variables */

  int id;
  cs_post_writer_t  *writer = NULL;

  /* Search for requested mesh */

  for (id = 0; id < _cs_post_n_writers; id++) {
    writer = _cs_post_writers + id;
    if (writer->id == writer_id)
      return true;
  }

  return false;
}

/*----------------------------------------------------------------------------
 * Check for the existence of a post-processing mesh of the given id.
 *
 * parameters:
 *   mesh_id <-- mesh id to check
 *
 * returns:
 *   true if mesh with this id exists, false otherwise
 *----------------------------------------------------------------------------*/

cs_bool_t
cs_post_mesh_exists(int  mesh_id)
{
  int id;
  cs_post_mesh_t  *post_mesh = NULL;

  /* Search for requested mesh */

  for (id = 0; id < _cs_post_n_meshes; id++) {
    post_mesh = _cs_post_meshes + id;
    if (post_mesh->id == mesh_id)
      return true;
  }

  return false;
}

/*----------------------------------------------------------------------------
 * Modify an existing post-processing mesh.
 *
 * The lists of cells or faces are redefined, for example to update an
 * extracted mesh based in "interesting" zones.
 *
 * It is not necessary to use this function if a mesh is simply deformed.
 *
 * parameters:
 *   mesh_id     <-- id of mesh to modify (< 0 reserved, > 0 for user)
 *   n_cells     <-- number of associated cells
 *   n_i_faces   <-- number of associated interior faces
 *   n_b_faces   <-- number of associated boundary faces
 *   cell_list   <-- list of associated cells
 *   i_face_list <-- list of associated interior faces
 *   b_face_list <-- list of associated boundary faces
 *
 *----------------------------------------------------------------------------*/

void
cs_post_modify_mesh(int       mesh_id,
                    cs_int_t  n_cells,
                    cs_int_t  n_i_faces,
                    cs_int_t  n_b_faces,
                    cs_int_t  cell_list[],
                    cs_int_t  i_face_list[],
                    cs_int_t  b_face_list[])
{
  /* local variables */

  int i, _mesh_id;
  char  *mesh_name = NULL;
  cs_post_mesh_t  *post_mesh = NULL;
  cs_post_writer_t    *writer = NULL;

  /* Get base structure (return if we do not own the mesh) */

  _mesh_id = _cs_post_mesh_id(mesh_id);
  post_mesh = _cs_post_meshes + _mesh_id;

  if (post_mesh->_exp_mesh == NULL)
    return;

  /* Remplace base structure */

  BFT_MALLOC(mesh_name,
             strlen(fvm_nodal_get_name(post_mesh->exp_mesh)) + 1,
             char);
  strcpy(mesh_name, fvm_nodal_get_name(post_mesh->exp_mesh));

  fvm_nodal_destroy(post_mesh->_exp_mesh);
  post_mesh->exp_mesh = NULL;

  _cs_post_define_mesh(post_mesh,
                       mesh_name,
                       n_cells,
                       n_i_faces,
                       n_b_faces,
                       cell_list,
                       i_face_list,
                       b_face_list);

  BFT_FREE(mesh_name);

  /* Update possible aliases */

  for (i = 0; i < _cs_post_n_meshes; i++) {
    if ((_cs_post_meshes + i)->alias == _mesh_id)
      (_cs_post_meshes + i)->exp_mesh
        = post_mesh->exp_mesh;
  }

  /* Divide polygons or polyhedra into simple elements */

  for (i = 0; i < post_mesh->n_writers; i++) {

    writer = _cs_post_writers + post_mesh->writer_id[i];
    _cs_post_divide_poly(post_mesh, writer);

  }
}

/*----------------------------------------------------------------------------
 * Return the next "reservable" (i.e. non-user) writer id available.
 *
 * Returns:
 *   the smallest negative integer present, -1
 *----------------------------------------------------------------------------*/

int
cs_post_get_free_writer_id(void)
{
  return (_cs_post_min_writer_id - 1);
}

/*----------------------------------------------------------------------------
 * Return the next "reservable" (i.e. non-user) mesh id available.
 *
 * Returns:
 *   the smallest negative integer present, -1
 *----------------------------------------------------------------------------*/

int
cs_post_get_free_mesh_id(void)
{
  return (_cs_post_min_mesh_id - 1);
}

/*----------------------------------------------------------------------------
 * Associate a writer with a post-processing mesh.
 *
 * If the writer only allows fixed (i.e. time-independent) meshes, the
 * mesh is exported immediately. Otherwise, output is delayed until
 * cs_post_write_meshes() is called for an active time step.
 *
 * parameters:
 *   mesh_id   <-- id of associated mesh
 *   writer_id <-- id of associated writer
 *----------------------------------------------------------------------------*/

void
cs_post_associate(int  mesh_id,
                  int  writer_id)
{
  int  i;
  int  _mesh_id, _writer_id;
  fvm_writer_time_dep_t mod_flag;

  cs_post_mesh_t  *post_mesh = NULL;
  cs_post_writer_t  *writer = NULL;

  /* Search for requested mesh and writer */

  _mesh_id = _cs_post_mesh_id(mesh_id);
  _writer_id = _cs_post_writer_id(writer_id);

  post_mesh = _cs_post_meshes + _mesh_id;

  /* Check that the writer is not already associated */

  for (i = 0; i < post_mesh->n_writers; i++) {
    if (post_mesh->writer_id[i] == _writer_id)
      break;
  }

  /* If the writer is not already associated, associate it */

  if (i >= post_mesh->n_writers) {

    post_mesh->n_writers += 1;
    BFT_REALLOC(post_mesh->writer_id,
                post_mesh->n_writers,
                cs_int_t);

    post_mesh->writer_id[post_mesh->n_writers - 1] = _writer_id;
    post_mesh->nt_last = - 1;

    /* Update structure */

    writer = _cs_post_writers + _writer_id;
    mod_flag = fvm_writer_get_time_dep(writer->writer);

    if (mod_flag < post_mesh->mod_flag_min)
      post_mesh->mod_flag_min = mod_flag;
    if (mod_flag > post_mesh->mod_flag_max)
      post_mesh->mod_flag_max = mod_flag;

    _cs_post_mod_flag_alias(_mesh_id);

    /* If we must compute the vertices displacement field, we need
       to save the initial vertex coordinates */

    if (   _cs_post_deformable == false
        && _cs_post_ini_vtx_coo == NULL
        && writer->write_displ == true) {

      cs_mesh_t *maillage = cs_glob_mesh;

      if (maillage->n_vertices > 0) {
        BFT_MALLOC(_cs_post_ini_vtx_coo,
                   maillage->n_vertices * 3,
                   cs_real_t);
        memcpy(_cs_post_ini_vtx_coo,
               maillage->vtx_coord,
               maillage->n_vertices * 3 * sizeof(cs_real_t));
      }

      _cs_post_deformable = true;

    }

    /* Divide polygons or polyhedra into simple elements */

    _cs_post_divide_poly(post_mesh, writer);

    /* If the writer only allows fixed (i.e. time-independent) meshes,
       output mesh immediately */

    if (mod_flag == FVM_WRITER_FIXED_MESH) {
      fvm_writer_set_mesh_time(writer->writer, 0, 0.0);
      fvm_writer_export_nodal(writer->writer, post_mesh->exp_mesh);
    }
  }

}

/*----------------------------------------------------------------------------
 * Update "active" or "inactive" flag of writers whose output frequency
 * is a divisor of the current time step number.
 *
 * parameters:
 *   nt_cur_abs <-- current time step number
 *----------------------------------------------------------------------------*/

void
cs_post_activate_if_default(int  nt_cur_abs)
{
  int  i;
  cs_post_writer_t  *writer;

  for (i = 0; i < _cs_post_n_writers; i++) {

    writer = _cs_post_writers + i;

    if (writer->frequency > 0) {
      if (nt_cur_abs % (writer->frequency) == 0)
        writer->active = 1;
      else
        writer->active = 0;
    }
    else
      writer->active = 0;

  }
}

/*----------------------------------------------------------------------------
 * Force the "active" or "inactive" flag for a specific writer or for all
 * writers for the current time step.
 *
 * parameters:
 *   writer_id <-- writer id, or 0 for all writers
 *   activate  <-- 0 to deactivate, 1 to activate
 *----------------------------------------------------------------------------*/

void
cs_post_activate_writer(int  writer_id,
                        int  activate)
{
  int i;
  cs_post_writer_t  *writer;

  if (writer_id != 0) {
    i = _cs_post_writer_id(writer_id);
    writer = _cs_post_writers + i;
    writer->active = (activate > 0) ? 1 : 0;
  }
  else {
    for (i = 0; i < _cs_post_n_writers; i++) {
      writer = _cs_post_writers + i;
      writer->active = (activate > 0) ? 1 : 0;
    }
  }
}

/*----------------------------------------------------------------------------
 * Return a pointer to the FVM library writer associated to a writer_id.
 *
 * parameters:
 *   writer_id <-- associated writer id
 *
 * Returns:
 *  a pointer to a fvm_writer_t structure
 *----------------------------------------------------------------------------*/

fvm_writer_t *
cs_post_get_writer(cs_int_t  writer_id)
{
  int  id;
  const cs_post_writer_t  *writer = NULL;

  id = _cs_post_writer_id(writer_id);
  writer = _cs_post_writers + id;

  return writer->writer;
}

/*----------------------------------------------------------------------------
 * Output post-processing meshes using associated writers.
 *
 * parameters:
 *   nt_cur_abs <-- current time step number
 *   t_cur_abs  <-- current physical time
 *----------------------------------------------------------------------------*/

void
cs_post_write_meshes(int     nt_cur_abs,
                     double  t_cur_abs)
{
  int  i;
  cs_post_mesh_t  *post_mesh;

  /* Loops on meshes and writers */

  for (i = 0; i < _cs_post_n_meshes; i++) {

    post_mesh = _cs_post_meshes + i;

    _cs_post_write_mesh(post_mesh,
                        nt_cur_abs,
                        t_cur_abs);

  }

}

/*----------------------------------------------------------------------------
 * Output a variable defined at cells or faces of a post-processing mesh
 * using associated writers.
 *
 * parameters:
 *   mesh_id     <-- id of associated mesh
 *   var_name    <-- name of variable to output
 *   var_dim     <-- 1 for scalar, 3 for vector
 *   interlace   <-- if a vector, true for interlaced values, false otherwise
 *   use_parent  <-- true if values are defined on "parent" mesh,
 *                   false if values are defined on post-processing mesh
 *   var_type    <-- variable's data type
 *   nt_cur_abs  <-- current time step number
 *   t_cur_abs   <-- current physical time
 *   cel_vals    <-- cell values
 *   i_face_vals <-- interior face values
 *   b_face_vals <-- boundary face values
 *----------------------------------------------------------------------------*/

void
cs_post_write_var(int              mesh_id,
                  const char      *var_name,
                  cs_int_t         var_dim,
                  cs_bool_t        interlace,
                  cs_bool_t        use_parent,
                  cs_post_type_t   var_type,
                  cs_int_t         nt_cur_abs,
                  cs_real_t        t_cur_abs,
                  const void      *cel_vals,
                  const void      *i_face_vals,
                  const void      *b_face_vals)
{
  cs_int_t  i;
  int       _mesh_id;


  fvm_interlace_t      _interlace;
  fvm_datatype_t       datatype;

  size_t       dec_ptr = 0;
  int          nbr_listes_parents = 0;
  fvm_lnum_t   dec_num_parent[2]  = {0, 0};
  cs_real_t   *var_tmp = NULL;
  cs_post_mesh_t  *post_mesh = NULL;
  cs_post_writer_t    *writer = NULL;

  const void  *var_ptr[2*9] = {NULL, NULL, NULL,
                               NULL, NULL, NULL,
                               NULL, NULL, NULL,
                               NULL, NULL, NULL,
                               NULL, NULL, NULL,
                               NULL, NULL, NULL};

  /* Initializations */

  _mesh_id = _cs_post_mesh_id(mesh_id);
  post_mesh = _cs_post_meshes + _mesh_id;

  if (interlace == true)
    _interlace = FVM_INTERLACE;
  else
    _interlace = FVM_NO_INTERLACE;

  datatype =  _cs_post_cnv_datatype(var_type);

  /* Assign appropriate array to FVM for output */

  /* Case of cells */
  /*---------------*/

  if (post_mesh->ent_flag[CS_POST_LOCATION_CELL] == 1) {

    if (use_parent == true) {
      nbr_listes_parents = 1;
      dec_num_parent[0] = 0;
    }
    else
      nbr_listes_parents = 0;

    var_ptr[0] = cel_vals;
    if (interlace == false) {
      if (use_parent == true)
        dec_ptr = cs_glob_mesh->n_cells_with_ghosts;
      else
        dec_ptr = fvm_nodal_get_n_entities(post_mesh->exp_mesh, 3);
      dec_ptr *= fvm_datatype_size[datatype];
      for (i = 1; i < var_dim; i++)
        var_ptr[i] = ((const char *)cel_vals) + i*dec_ptr;
    }
  }

  /* Case of faces */
  /*---------------*/

  else if (   post_mesh->ent_flag[CS_POST_LOCATION_I_FACE] == 1
           || post_mesh->ent_flag[CS_POST_LOCATION_B_FACE] == 1) {

    /* In case of indirection, all that is necessary is to set pointers */

    if (use_parent == true) {

      nbr_listes_parents = 2;
      dec_num_parent[0] = 0;
      dec_num_parent[1] = cs_glob_mesh->n_b_faces;

      if (post_mesh->ent_flag[CS_POST_LOCATION_B_FACE] == 1) {
        if (interlace == false) {
          dec_ptr = cs_glob_mesh->n_b_faces * fvm_datatype_size[datatype];
          for (i = 0; i < var_dim; i++)
            var_ptr[i] = ((const char *)b_face_vals) + i*dec_ptr;
        }
        else
          var_ptr[0] = b_face_vals;
      }

      if (post_mesh->ent_flag[CS_POST_LOCATION_I_FACE] == 1) {
        if (interlace == false) {
          dec_ptr = cs_glob_mesh->n_i_faces * fvm_datatype_size[datatype];
          for (i = 0; i < var_dim; i++)
            var_ptr[var_dim + i] = ((const char *)i_face_vals) + i*dec_ptr;
        }
        else
          var_ptr[1] = i_face_vals;
      }

    }

    /* With no indirection, we must switch to a variable defined on two
       lists of faces to a variable defined on one list */

    else {

      nbr_listes_parents = 0;

      if (post_mesh->ent_flag[CS_POST_LOCATION_B_FACE] == 1) {

        /* Case where a variable is defined both on boundary and
           interior faces: we must switch to a single list, as
           indirection is not used */

        if (post_mesh->ent_flag[CS_POST_LOCATION_I_FACE] == 1) {

          BFT_MALLOC(var_tmp,
                     (   post_mesh->n_i_faces
                      +  post_mesh->n_b_faces) * var_dim,
                     cs_real_t);

          _cs_post_assmb_var_faces(post_mesh->exp_mesh,
                                   post_mesh->n_i_faces,
                                   post_mesh->n_b_faces,
                                   var_dim,
                                   _interlace,
                                   i_face_vals,
                                   b_face_vals,
                                   var_tmp);

          _interlace = FVM_NO_INTERLACE;

          dec_ptr = fvm_datatype_size[datatype] * (  post_mesh->n_i_faces
                                                   + post_mesh->n_b_faces);

          for (i = 0; i < var_dim; i++)
            var_ptr[i] = ((char *)var_tmp) + i*dec_ptr;

        }

        /* Case where we only have boundary faces */

        else {

          if (interlace == false) {
            dec_ptr = fvm_datatype_size[datatype] * post_mesh->n_b_faces;
            for (i = 0; i < var_dim; i++)
              var_ptr[i] = ((const char *)b_face_vals) + i*dec_ptr;
          }
          else
            var_ptr[0] = b_face_vals;
        }

      }

      /* Case where we only have interior faces */

      else if (post_mesh->ent_flag[CS_POST_LOCATION_I_FACE] == 1) {

        if (interlace == false) {
          dec_ptr = fvm_datatype_size[datatype] * post_mesh->n_i_faces;
          for (i = 0; i < var_dim; i++)
            var_ptr[i] = ((const char *)i_face_vals) + i*dec_ptr;
        }
        else
          var_ptr[0] = i_face_vals;
      }

    }

  }

  /* Effective output: loop on writers */
  /*-----------------------------------*/

  for (i = 0; i < post_mesh->n_writers; i++) {

    writer = _cs_post_writers + post_mesh->writer_id[i];

    if (writer->active == 1)
      fvm_writer_export_field(writer->writer,
                              post_mesh->exp_mesh,
                              var_name,
                              FVM_WRITER_PER_ELEMENT,
                              var_dim,
                              _interlace,
                              nbr_listes_parents,
                              dec_num_parent,
                              datatype,
                              (int)nt_cur_abs,
                              (double)t_cur_abs,
                              (const void * *)var_ptr);

  }

  /* Free memory (if both interior and boundary faces present) */

  if (var_tmp != NULL)
    BFT_FREE(var_tmp);
}

/*----------------------------------------------------------------------------
 * Output a variable defined at vertices of a post-processing mesh using
 * associated writers.
 *
 * parameters:
 *   mesh_id     <-- id of associated mesh
 *   var_name    <-- name of variable to output
 *   var_dim     <-- 1 for scalar, 3 for vector
 *   interlace   <-- if a vector, true for interlaced values, false otherwise
 *   use_parent  <-- true if values are defined on "parent" mesh,
 *                   false if values are defined on post-processing mesh
 *   var_type    <-- variable's data type
 *   nt_cur_abs  <-- current time step number
 *   t_cur_abs   <-- current physical time
 *   vtx_vals    <-- vertex values
 *----------------------------------------------------------------------------*/

void
cs_post_write_vertex_var(int              mesh_id,
                         const char      *var_name,
                         cs_int_t         var_dim,
                         cs_bool_t        interlace,
                         cs_bool_t        use_parent,
                         cs_post_type_t   var_type,
                         cs_int_t         nt_cur_abs,
                         cs_real_t        t_cur_abs,
                         const void      *vtx_vals)
{
  cs_int_t  i;
  int       _mesh_id;


  cs_post_mesh_t  *post_mesh;
  cs_post_writer_t    *writer;
  fvm_interlace_t      _interlace;
  fvm_datatype_t       datatype;

  size_t       dec_ptr = 0;
  int          nbr_listes_parents = 0;
  fvm_lnum_t   dec_num_parent[1]  = {0};

  const void  *var_ptr[9] = {NULL, NULL, NULL,
                             NULL, NULL, NULL,
                             NULL, NULL, NULL};

  /* Initializations */

  _mesh_id = _cs_post_mesh_id(mesh_id);
  post_mesh = _cs_post_meshes + _mesh_id;

  if (interlace == true)
    _interlace = FVM_INTERLACE;
  else
    _interlace = FVM_NO_INTERLACE;

  assert(   sizeof(cs_real_t) == sizeof(double)
         || sizeof(cs_real_t) == sizeof(float));

  datatype =  _cs_post_cnv_datatype(var_type);

  /* Assign appropriate array to FVM for output */

  if (use_parent == true)
    nbr_listes_parents = 1;
  else
    nbr_listes_parents = 0;

  var_ptr[0] = vtx_vals;
  if (interlace == false) {
    if (use_parent == true)
      dec_ptr = cs_glob_mesh->n_vertices;
    else
      dec_ptr =   fvm_nodal_get_n_entities(post_mesh->exp_mesh, 0)
                * fvm_datatype_size[datatype];
    for (i = 1; i < var_dim; i++)
      var_ptr[i] = ((const char *)vtx_vals) + i*dec_ptr;
  }

  /* Effective output: loop on writers */
  /*-----------------------------------*/

  for (i = 0; i < post_mesh->n_writers; i++) {

    writer = _cs_post_writers + post_mesh->writer_id[i];

    if (writer->active == 1)
      fvm_writer_export_field(writer->writer,
                              post_mesh->exp_mesh,
                              var_name,
                              FVM_WRITER_PER_NODE,
                              var_dim,
                              _interlace,
                              nbr_listes_parents,
                              dec_num_parent,
                              datatype,
                              (int)nt_cur_abs,
                              (double)t_cur_abs,
                              (const void * *)var_ptr);

  }

}

/*----------------------------------------------------------------------------
 * Update references to parent mesh of post-processing meshes in case of
 * computational mesh cell renumbering.
 *
 * This function may be called only once, after possible renumbering of cells,
 * to update existing post-processing meshes. Post-processing meshes defined
 * after renumbering will automatically be based upon the new numbering,
 * so this function will not need to be called again.
 *
 * parameters:
 *   init_cell_num <-- initial cell numbering (1 to n, new -> old)
 *----------------------------------------------------------------------------*/

void
cs_post_renum_cells(const cs_int_t  init_cell_num[])
{
  int       i;
  cs_int_t  icel;
  cs_int_t  n_elts;

  cs_int_t  *renum_ent_parent = NULL;

  cs_bool_t  a_traiter = false;

  cs_post_mesh_t   *post_mesh;
  const cs_mesh_t  *maillage = cs_glob_mesh;

  if (init_cell_num == NULL)
    return;

  /* Loop on meshes */

  for (i = 0; i < _cs_post_n_meshes; i++) {

    post_mesh = _cs_post_meshes + i;

    if (post_mesh->ent_flag[CS_POST_LOCATION_CELL] > 0)
      a_traiter = true;
  }

  if (a_traiter == true) {

    /* Prepare renumbering */

    n_elts = maillage->n_cells;

    BFT_MALLOC(renum_ent_parent, n_elts, cs_int_t);

    for (icel = 0; icel < maillage->n_cells; icel++)
      renum_ent_parent[init_cell_num[icel] - 1] = icel + 1;

    /* Effective modification */

    for (i = 0; i < _cs_post_n_meshes; i++) {

      post_mesh = _cs_post_meshes + i;

      if (   post_mesh->_exp_mesh != NULL
          && post_mesh->ent_flag[CS_POST_LOCATION_CELL] > 0) {

        fvm_nodal_change_parent_num(post_mesh->_exp_mesh,
                                    renum_ent_parent,
                                    3);

      }

    }

    BFT_FREE(renum_ent_parent);

  }

}

/*----------------------------------------------------------------------------
 * Update references to parent mesh of post-processing meshes in case of
 * computational mesh interior and/or boundary faces renumbering.
 *
 * This function may be called only once, after possible renumbering of faces,
 * to update existing post-processing meshes. Post-processing meshes defined
 * after renumbering will automatically be based upon the new numbering,
 * so this function will not need to be called again.
 *
 * parameters:
 *   init_i_face_num <-- initial interior numbering (1 to n, new -> old)
 *   init_b_face_num <-- initial boundary numbering (1 to n, new -> old)
 *----------------------------------------------------------------------------*/

void
cs_post_renum_faces(const cs_int_t  init_i_face_num[],
                    const cs_int_t  init_b_face_num[])
{
  int       i;
  cs_int_t  ifac;
  cs_int_t  n_elts;

  cs_int_t  *renum_ent_parent = NULL;

  cs_bool_t  a_traiter = false;

  cs_post_mesh_t   *post_mesh;
  const cs_mesh_t  *maillage = cs_glob_mesh;

  /* Loop on meshes */

  for (i = 0; i < _cs_post_n_meshes; i++) {

    post_mesh = _cs_post_meshes + i;

    if (   post_mesh->ent_flag[CS_POST_LOCATION_I_FACE] > 0
        || post_mesh->ent_flag[CS_POST_LOCATION_B_FACE] > 0) {
      a_traiter = true;
    }

  }

  if (a_traiter == true) {

    /* Prepare renumbering */

    n_elts = maillage->n_i_faces + maillage->n_b_faces;

    BFT_MALLOC(renum_ent_parent, n_elts, cs_int_t);

    if (init_b_face_num == NULL) {
      for (ifac = 0; ifac < maillage->n_b_faces; ifac++)
        renum_ent_parent[ifac] = ifac + 1;
    }
    else {
      for (ifac = 0; ifac < maillage->n_b_faces; ifac++)
        renum_ent_parent[init_b_face_num[ifac] - 1] = ifac + 1;
    }

    if (init_i_face_num == NULL) {
      for (ifac = 0, i = maillage->n_b_faces;
           ifac < maillage->n_i_faces;
           ifac++, i++)
        renum_ent_parent[maillage->n_b_faces + ifac]
          = maillage->n_b_faces + ifac + 1;
    }
    else {
      for (ifac = 0, i = maillage->n_b_faces;
           ifac < maillage->n_i_faces;
           ifac++, i++)
        renum_ent_parent[maillage->n_b_faces + init_i_face_num[ifac] - 1]
          = maillage->n_b_faces + ifac + 1;
    }

    /* Effective modification */

    for (i = 0; i < _cs_post_n_meshes; i++) {

      post_mesh = _cs_post_meshes + i;

      if (   post_mesh->_exp_mesh != NULL
          && (   post_mesh->ent_flag[CS_POST_LOCATION_I_FACE] > 0
              || post_mesh->ent_flag[CS_POST_LOCATION_B_FACE] > 0)) {

        fvm_nodal_change_parent_num(post_mesh->_exp_mesh,
                                    renum_ent_parent,
                                    2);

      }

    }

    BFT_FREE(renum_ent_parent);
  }

}

/*----------------------------------------------------------------------------
 * Destroy all structures associated with post-processing
 *----------------------------------------------------------------------------*/

void
cs_post_finalize(void)
{
  int i;
  cs_post_mesh_t  *post_mesh = NULL;

  /* Timings */

  for (i = 0; i < _cs_post_n_writers; i++) {
    double m_wtime = 0.0, m_cpu_time = 0.0, c_wtime = 0.0, c_cpu_time = 0.;
    fvm_writer_get_times((_cs_post_writers + i)->writer,
                         &m_wtime, &m_cpu_time, &c_wtime, &c_cpu_time);
    bft_printf(_("\n"
                 "Writing of \"%s\" (%s) summary:\n"
                 "\n"
                 "  CPU time for meshes:              %12.3f\n"
                 "  CPU time for variables:           %12.3f\n"
                 "\n"
                 "  Elapsed time for meshes:          %12.3f\n"
                 "  Elapsed time for variables:       %12.3f\n"),
               fvm_writer_get_name((_cs_post_writers + i)->writer),
               fvm_writer_get_format((_cs_post_writers + i)->writer),
               m_cpu_time, c_cpu_time, m_wtime, c_wtime);
  }

  /* Initial coordinates (if mesh is deformable) */

  if (_cs_post_ini_vtx_coo != NULL)
    BFT_FREE(_cs_post_ini_vtx_coo);

  /* Exportable meshes */

  for (i = 0; i < _cs_post_n_meshes; i++) {
    post_mesh = _cs_post_meshes + i;
    if (post_mesh->_exp_mesh != NULL)
      fvm_nodal_destroy(post_mesh->_exp_mesh);
    BFT_FREE(post_mesh->writer_id);
  }

  BFT_FREE(_cs_post_meshes);

  _cs_post_min_mesh_id = -2;
  _cs_post_n_meshes = 0;
  _cs_post_n_meshes_max = 0;

  /* Writers */

  for (i = 0; i < _cs_post_n_writers; i++)
    fvm_writer_finalize((_cs_post_writers + i)->writer);

  BFT_FREE(_cs_post_writers);

  _cs_post_n_writers = 0;
  _cs_post_n_writers_max = 0;

  /* Registered processings if necessary */

  if (_cs_post_nbr_var_tp_max > 0) {
    BFT_FREE(_cs_post_f_var_tp);
    BFT_FREE(_cs_post_i_var_tp);
  }
}

/*----------------------------------------------------------------------------
 * Initialize main post-processing writer
 *----------------------------------------------------------------------------*/

void
cs_post_init_main_writer(void)
{
  /* Default values */

  cs_int_t  indic_vol = -1, indic_brd = -1, indic_syr = -1, indic_ze = -1;
  cs_int_t  indic_mod = -1;
  char  fmtchr[32 + 1] = "";
  char  optchr[96 + 1] = "";
  cs_int_t  ntchr = -1;

  const char  nomcas[] = "chr";
  const char  nomrep_ens[] = "chr.ensight";
  const char  nomrep_def[] = ".";
  const char *nomrep = NULL;

  const cs_int_t  writer_id = -1; /* Default (main) writer id */

  /* Get parameters from Fortran COMMON blocks */

  CS_PROCF(inipst, INIPST)(&indic_vol,
                           &indic_brd,
                           &indic_syr,
                           &indic_ze,
                           &indic_mod,
                           &ntchr,
                           fmtchr,
                           optchr);

  fmtchr[32] = '\0';
  optchr[96] = '\0';

  if (indic_vol == 0 && indic_brd == 0 && indic_syr == 0)
    return;

  /* Create default writer */

  if (fmtchr[0] == 'e' || fmtchr[0] == 'E')
    nomrep = nomrep_ens;
  else
    nomrep = nomrep_def;

  cs_post_add_writer(writer_id,
                     nomcas,
                     nomrep,
                     fmtchr,
                     optchr,
                     indic_mod,
                     ntchr);

}

/*----------------------------------------------------------------------------
 * Initialize main post-processing meshes
 *----------------------------------------------------------------------------*/

void
cs_post_init_main_meshes(void)
{
  /* Default values */

  cs_int_t  indic_vol = -1, indic_brd = -1, indic_syr = -1, indic_ze = -1;
  cs_int_t  indic_mod = -1;
  char  fmtchr[32 + 1] = "";
  char  optchr[96 + 1] = "";
  cs_int_t  ntchr = -1;

  cs_int_t  mesh_id = -1;

  const cs_int_t  writer_id = -1; /* Default (main) writer id */

  /* Get parameters from Fortran COMMON blocks */

  CS_PROCF(inipst, INIPST)(&indic_vol,
                           &indic_brd,
                           &indic_syr,
                           &indic_ze,
                           &indic_mod,
                           &ntchr,
                           fmtchr,
                           optchr);

  fmtchr[32] = '\0';
  optchr[96] = '\0';

  /* Definition of post-processing meshes */

  if (cs_glob_mesh->n_i_faces > 0 || cs_glob_mesh->n_b_faces > 0) {

    /*
      If the faces -> vertices connectivity is available, we
      may rebuild the nodal connectivity for post-processing
      (usual mechanism).
    */

    if (indic_vol > 0) { /* Volume mesh */

      mesh_id = -1; /* Reserved mesh id */

      cs_post_add_mesh(mesh_id,
                       _("Fluid volume"),
                       cs_glob_mesh->n_cells,
                       0,
                       0,
                       NULL,
                       NULL,
                       NULL);

      cs_post_associate(mesh_id, writer_id);

    }

    if (indic_brd > 0) { /* Boundary mesh */

      mesh_id = -2;  /* Reserved mesh id */

      cs_post_add_mesh(mesh_id,
                       _("Boundary"),
                       0,
                       0,
                       cs_glob_mesh->n_b_faces,
                       NULL,
                       NULL,
                       NULL);

      cs_post_associate(mesh_id, writer_id);

    }

  } /* End if cs_glob_mesh->n_i_faces > 0 || cs_glob_mesh->n_b_faces > 0 */

  /*
    If we do not have the faces -> vertices connectivity, we may not
    rebuild the nodal connectivity, so we must obtain it through
    another means.

    This only happens when we have directly read a mesh in the solcom
    format, in which  the nodal connectivity has already been read and
    assigned to a post-processing mesh
    (see LETGEO and cs_maillage_solcom_lit).
  */

  else if (indic_vol > 0) {

    mesh_id = -1;

    if (cs_post_mesh_exists(mesh_id))
      cs_post_associate(mesh_id, writer_id);

  }
}

/*----------------------------------------------------------------------------
 * Initialize post-processing writer with same format and associated
 * options as default writer, but no time dependency, intended to
 * troubleshoot errors.
 *----------------------------------------------------------------------------*/

void
cs_post_init_error_writer(void)
{
  /* Default values */

  cs_int_t  indic_vol = -1, indic_brd = -1, indic_syr = -1, indic_ze = -1;
  cs_int_t  indic_mod = -1;
  char  fmtchr[32 + 1] = "";
  char  optchr[96 + 1] = "";
  cs_int_t  ntchr = -1;

  const char  nomcas[] = "error";
  const char  nomrep_ens[] = "error.ensight";
  const char  nomrep_def[] = ".";
  const char *nomrep = NULL;

  const int writer_id = -2;

  if (cs_post_writer_exists(writer_id))
    return;

  /* Get parameters from Fortran COMMON blocks */

  CS_PROCF(inipst, INIPST)(&indic_vol,
                           &indic_brd,
                           &indic_syr,
                           &indic_ze,
                           &indic_mod,
                           &ntchr,
                           fmtchr,
                           optchr);

  fmtchr[32] = '\0';
  optchr[96] = '\0';

  /* Create default writer */

  if (fmtchr[0] == 'e' || fmtchr[0] == 'E')
    nomrep = nomrep_ens;
  else
    nomrep = nomrep_def;

  cs_post_add_writer(writer_id,
                     nomcas,
                     nomrep,
                     fmtchr,
                     optchr,
                     -1, /* No time dependency here */
                     ntchr);
}

/*----------------------------------------------------------------------------
 * Initialize post-processing writer with same format and associated
 * options as default writer, but no time dependency, and associate
 * and output global volume mesh.
 *
 * This is intended to help troubleshoot errors using fields based
 * on cells.
 *
 * returns:
 *   id of error output mesh (< 0), or 0 if all writers are deactivated
 *----------------------------------------------------------------------------*/

int
cs_post_init_error_writer_cells(void)
{
  int mesh_id = 0;

  const cs_mesh_t *mesh = cs_glob_mesh;

  /* If post-processing is active, output info */
  /*-------------------------------------------*/

  if (mesh->i_face_vtx_idx != NULL || mesh->b_face_vtx_idx != NULL) {

    const int writer_id = -2;
    const cs_int_t n_cells = mesh->n_cells;

    cs_post_init_error_writer();

    mesh_id = cs_post_get_free_mesh_id();

    cs_post_add_mesh(mesh_id,
                     _("Calculation domain"),
                     n_cells,
                     0,
                     0,
                     NULL,
                     NULL,
                     NULL);

    cs_post_associate(mesh_id, writer_id);

    cs_post_activate_writer(writer_id, 1);

    cs_post_write_meshes(-1, 0.0);
  }

  return mesh_id;
}

/*----------------------------------------------------------------------------
 * Register a processing of a time-dependent variable to the call to PSTVAR.
 *
 * The instance identifier associated with the function allows registering
 * the same function several times, with a diferent identifier allowing the
 * function to select a specific operation or data.
 *
 * parameters:
 *   function    <-- function to register
 *   instance_id <-- instance id associated with this registration
 *----------------------------------------------------------------------------*/

void
cs_post_add_time_dep_var(cs_post_time_dep_var_t  *function,
                         int                      instance_id)
{
  /* Resize array of registered post-processings if necessary */

  if (_cs_post_nbr_var_tp <= _cs_post_nbr_var_tp_max) {
    if (_cs_post_nbr_var_tp_max == 0)
      _cs_post_nbr_var_tp_max = 8;
    else
      _cs_post_nbr_var_tp_max *= 2;
    BFT_REALLOC(_cs_post_f_var_tp,
                _cs_post_nbr_var_tp_max,
                cs_post_time_dep_var_t *);
    BFT_REALLOC(_cs_post_i_var_tp, _cs_post_nbr_var_tp_max, cs_int_t);
  }

  /* Add a post-processing */

  _cs_post_f_var_tp[_cs_post_nbr_var_tp] = function;
  _cs_post_i_var_tp[_cs_post_nbr_var_tp] = instance_id;

  _cs_post_nbr_var_tp += 1;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
