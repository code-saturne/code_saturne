/*============================================================================
 * Functions associated to ALE formulation
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_selector.h"

#include "cs_interface.h"

#include "cs_base.h"
#include "cs_benchmark.h"
#include "cs_gradient.h"
#include "cs_gui.h"
#include "cs_gui_mesh.h"
#include "cs_gui_output.h"
#include "cs_gradient.h"
#include "cs_gradient_perio.h"
#include "cs_join.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_matrix_default.h"
#include "cs_mesh.h"
#include "cs_mesh_coherency.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_to_builder.h"
#include "cs_multigrid.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_preprocess.h"
#include "cs_prototypes.h"
#include "cs_renumber.h"
#include "cs_time_step.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_turbomachinery.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/* Turbomachinery structure */

typedef struct {

  cs_turbomachinery_model_t  model;             /* Turbomachinery model type */

  char                      *rotor_cells_c;     /* Rotor cells selection
                                                   criteria */

  cs_mesh_t                 *reference_mesh;    /* Reference mesh (before
                                                   rotation and joining) */

  cs_lnum_t                  n_b_faces_ref;     /* Reference number of
                                                   boundary faces */

  cs_lnum_t                  n_rotor_vtx;       /* Size of rotor_vtx array */
  cs_lnum_t                 *rotor_vtx;         /* List of vertices related
                                                   to the rotor in the
                                                   reference mesh */

  int                       *cell_rotor_num;    /* cell rotation axis number */

  double  omega;
  double  angle;
  double  rotation_axis[3];
  double  rotation_invariant[3];

  bool active;

} cs_turbomachinery_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

cs_turbomachinery_t  *cs_glob_turbomachinery = NULL;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void cs_f_map_turbomachinery_module(cs_int_t    *iturbo,
                                    cs_real_t    rotax[],
                                    int        **irotce);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Function for selection of faces with joining changes.
 *
 * parameters:
 *   input    <-- pointer to input (count[0]: n_previous, count[1]: n_current)
 *   n_faces  --> number of selected faces
 *   face_ids --> array of selected face ids (0 to n-1 numbering)
 *----------------------------------------------------------------------------*/

static void
_post_error_faces_select(void         *input,
                         cs_lnum_t    *n_faces,
                         cs_lnum_t   **face_ids)
{
  cs_lnum_t face_id;

  cs_lnum_t _n_faces = 0;
  cs_lnum_t *_face_ids = NULL;

  const cs_lnum_t  *count = input;

  BFT_MALLOC(_face_ids, count[1], cs_lnum_t);

  for (face_id = count[0]; face_id < count[1]; face_id++)
    _face_ids[_n_faces++] = face_id;

  *n_faces = _n_faces;
  *face_ids = _face_ids;
}

/*----------------------------------------------------------------------------
 * Create an empty turbomachinery structure
 *----------------------------------------------------------------------------*/

static cs_turbomachinery_t *
_turbomachinery_create(void)
{
  cs_turbomachinery_t  *tbm = NULL;

  BFT_MALLOC(tbm, 1, cs_turbomachinery_t);

  tbm->rotor_cells_c = NULL;
  tbm->reference_mesh = cs_mesh_create();
  tbm->n_b_faces_ref = -1;
  tbm->cell_rotor_num = NULL;
  tbm->rotor_vtx = NULL;
  tbm->model = CS_TURBOMACHINERY_NONE;

  for (int i = 0; i < 3; i++) {
    tbm->angle = 0;
    tbm->omega = 0;
    tbm->rotation_axis[i] = 0.;
    tbm->rotation_invariant[i] = 0.;
  }

  return tbm;
}

/*----------------------------------------------------------------------------
 * Compute rotation matrix
 *
 * parameters:
 *   theta           <-- rotation angle, in radians
 *   axis            <-- components of rotation axis direction vector (3)
 *   invariant_point <-- components of invariant point (3)
 *   matrix          <-- resulting rotation matrix
 *---------------------------------------------------------------------------*/

static void
_rotation_matrix(double        theta,
                 const double  axis[3],
                 const double  invariant_point[3],
                 double        matrix[3][4])
{
  int  i, j;
  double norm;
  double direction[3];
  double rot[3][3];

  const double cost = cos(theta);
  const double sint = sin(theta);
  const double onemcost = (1.0 - cost);

  /* Compute the rotation matrix, using formula:
   *  R = (1-cos(theta))axis.transp(axis) + cos(theta)I + sin(theta)V
   *
   *           [ 0            -direction(3)  direction(2)]
   *  with V = [ direction(3)       0       -direction(1)]
   *           [-direction(2)  direction(1)       0      ]
   */

  norm = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);

  direction[0] = axis[0] / norm;
  direction[1] = axis[1] / norm;
  direction[2] = axis[2] / norm;

  /* first row of rotation maxtrix */
  rot[0][0] = onemcost*direction[0]*direction[0] + cost;
  rot[0][1] = onemcost*direction[0]*direction[1] - sint*direction[2];
  rot[0][2] = onemcost*direction[0]*direction[2] + sint*direction[1];

  /* second row of rotation maxtrix */
  rot[1][0] = onemcost*direction[1]*direction[0] + sint*direction[2];
  rot[1][1] = onemcost*direction[1]*direction[1] + cost;
  rot[1][2] = onemcost*direction[1]*direction[2] - sint*direction[0];

  /* third row of rotation maxtrix */
  rot[2][0] = onemcost*direction[2]*direction[0] - sint*direction[1];
  rot[2][1] = onemcost*direction[2]*direction[1] + sint*direction[0];
  rot[2][2] = onemcost*direction[2]*direction[2] + cost;

  /* Now compute full rotation matrix in homogeneous coordinates,
   * accounting for invariant point of coordiantes t[], with the formula:
   *
   *     [1 0 0 t[0]] [r[0][0] r[0][1] r[0][3] 0] [1 0 0 -t[0]]
   * M = [0 1 0 t[1]].[r[1][0] r[1][1] r[1][3] 0].[0 1 0 -t[1]]
   *     [0 0 1 t[2]] [r[2][0] r[2][1] r[2][3] 0] [0 0 1 -t[2]]
   *     [0 0 0 1   ] [0       0       0       1] [0 0 0  1]
   */

  for (i = 0; i < 3; i++) {       /* rotation part of matrix */
    for (j = 0; j < 3; j++) {
      matrix[i][j] = rot[i][j];
    }
  }

  for (i = 0; i < 3; i++) {
    matrix[i][3] = invariant_point[i];
    for (j = 0; j < 3; j++)
      matrix[i][3] -= rot[i][j]*invariant_point[j];
  }
}

/*----------------------------------------------------------------------------
 * Compute a matrix/vector product to apply a transformation to a vector.
 *
 * parameters:
 *   m[3][4] <-- matrix of the transformation in homogeneous coord.
 *               last line = [0; 0; 0; 1] (Not used here)
 *   c[3]    <-> coordinates
 *----------------------------------------------------------------------------*/

static inline void
_apply_vector_transfo(double     matrix[3][4],
                      cs_real_t  c[3])
{
  int  i, j;

  double  c_a[4] = {c[0], c[1], c[2], 1.}; /* homogeneous coords */
  double  c_b[3] = {0, 0, 0};

  for (i = 0; i < 3; i++)
    for (j = 0; j < 4; j++)
      c_b[i] += matrix[i][j]*c_a[j];

  for (i = 0; i < 3; i++)
    c[i] = c_b[i];
}

/*----------------------------------------------------------------------------
 * Compute a matrix/vector product to apply a rotation to a vector.
 *
 * parameters:
 *   m[3][4] <-- matrix of the transformation in homogeneous coord.
 *               last line = [0; 0; 0; 1] (Not used here)
 *   v[3]    <-> vector
 *----------------------------------------------------------------------------*/

static inline void
_apply_vector_rotation(double     m[3][4],
                       cs_real_t  v[3])
{
  double  t[3] = {v[0], v[1], v[0]};

  for (int i = 0; i < 3; i++)
    v[i] = m[i][0]*t[0] + m[i][1]*t[1] + m[i][2]*t[2];
}

/*----------------------------------------------------------------------------
 * Compute a matrix * tensor * Tmatrix product to apply a rotation to a
 * given symmetric interleaved tensor
 *
 * parameters:
 *   matrix[3][4]        --> transformation matrix in homogeneous coords.
 *                           last line = [0; 0; 0; 1] (Not used here)
 *   tensor[6]           <-> incoming (6) symmetric tensor
 *----------------------------------------------------------------------------*/

static inline void
_apply_sym_tensor_rotation(cs_real_t  matrix[3][4],
                           cs_real_t  t[6])
{
  int i, j, k, l;

  double  _t[3][3], _t0[3][3];

  _t0[0][0] = t[0];
  _t0[1][1] = t[1];
  _t0[2][2] = t[2];
  _t0[0][1] = t[3];
  _t0[1][0] = t[3];
  _t0[1][2] = t[4];
  _t0[2][1] = t[4];
  _t0[0][2] = t[5];
  _t0[2][0] = t[5];

  for (k = 0; k < 3; k++) {
    for (j = 0; j < 3; j++) {
      _t[k][j] = 0.;
      for (l = 0; l < 3; l++)
        _t[k][j] += matrix[j][l] * _t0[k][l];
    }
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      _t0[i][j] = 0.;
      for (k = 0; k < 3; k++)
        _t0[i][j] += matrix[i][k] * _t[k][j];
    }
  }

  t[0] = _t0[0][0];
  t[1] = _t0[1][1];
  t[2] = _t0[2][2];
  t[3] = _t0[0][1];
  t[3] = _t0[1][0];
  t[4] = _t0[2][1];
  t[5] = _t0[2][0];
}

/*----------------------------------------------------------------------------
 * Compute velocity relative to fixed coordinates at a given point
 *
 * parameters:
 *   omega           <-- rotation velocity
 *   axis            <-- components of rotation axis direction vector (3)
 *   invariant_point <-- components of invariant point (3)
 *   coords          <-- coordinates at point
 *   velocity        --> resulting relative velocity
 *---------------------------------------------------------------------------*/

static void
_relative_velocity(double        omega,
                   const double  axis[3],
                   const double  invariant_point[3],
                   const double  coords[3],
                   cs_real_t     velocity[3])
{
  velocity[0] = (- axis[2] * (coords[1] - invariant_point[1])
                 + axis[1] * (coords[2] - invariant_point[2])) * omega;
  velocity[1] = (  axis[2] * (coords[0] - invariant_point[0])
                 - axis[0] * (coords[2] - invariant_point[2])) * omega;
  velocity[2] = (- axis[1] * (coords[0] - invariant_point[0])
                 + axis[0] * (coords[1] - invariant_point[1])) * omega;
}

/*----------------------------------------------------------------------------
 * Duplicate a cs_mesh_t structure.
 *
 * Note that some fields which are recomputable are not copied.
 *
 * parameters:
 *   mesh      <-- reference mesh
 *   mesh_copy <-> mesh copy
 *----------------------------------------------------------------------------*/

static void
_copy_mesh(const cs_mesh_t  *mesh,
           cs_mesh_t        *mesh_copy)
{
  cs_lnum_t n_elts;

  /* General features */

  mesh_copy->dim        = mesh->dim;
  mesh_copy->domain_num = mesh->domain_num;
  mesh_copy->n_domains  = mesh->n_domains;

  /* Local dimensions */

  mesh_copy->n_cells    = mesh->n_cells;
  mesh_copy->n_i_faces  = mesh->n_i_faces;
  mesh_copy->n_b_faces  = mesh->n_b_faces;
  mesh_copy->n_vertices = mesh->n_vertices;

  mesh_copy->i_face_vtx_connect_size = mesh->i_face_vtx_connect_size;
  mesh_copy->b_face_vtx_connect_size = mesh->b_face_vtx_connect_size;

  /* Local structures */

  BFT_MALLOC(mesh_copy->vtx_coord, (3*mesh->n_vertices), cs_real_t);
  memcpy(mesh_copy->vtx_coord,
         mesh->vtx_coord,
         3*mesh->n_vertices*sizeof(cs_real_t));

  BFT_MALLOC(mesh_copy->i_face_cells, mesh->n_i_faces, cs_lnum_2_t);
  memcpy(mesh_copy->i_face_cells,
         mesh->i_face_cells,
         mesh->n_i_faces*sizeof(cs_lnum_2_t));

  if (mesh->n_b_faces > 0) {
    BFT_MALLOC(mesh_copy->b_face_cells, mesh->n_b_faces, cs_lnum_t);
    memcpy(mesh_copy->b_face_cells,
           mesh->b_face_cells,
           mesh->n_b_faces*sizeof(cs_lnum_t));
  }

  BFT_MALLOC(mesh_copy->i_face_vtx_idx, mesh->n_i_faces + 1, cs_lnum_t);
  memcpy(mesh_copy->i_face_vtx_idx,
         mesh->i_face_vtx_idx,
         (mesh->n_i_faces + 1)*sizeof(cs_lnum_t));

  BFT_MALLOC(mesh_copy->i_face_vtx_lst,
             mesh->i_face_vtx_connect_size,
             cs_lnum_t);
  memcpy(mesh_copy->i_face_vtx_lst, mesh->i_face_vtx_lst,
         mesh->i_face_vtx_connect_size*sizeof(cs_lnum_t));

  BFT_MALLOC(mesh_copy->b_face_vtx_idx, mesh->n_b_faces + 1, cs_lnum_t);
  memcpy(mesh_copy->b_face_vtx_idx,
         mesh->b_face_vtx_idx,
         (mesh->n_b_faces + 1)*sizeof(cs_lnum_t));

  if (mesh->b_face_vtx_connect_size > 0) {
    BFT_MALLOC(mesh_copy->b_face_vtx_lst,
               mesh->b_face_vtx_connect_size,
               cs_lnum_t);
    memcpy(mesh_copy->b_face_vtx_lst,
           mesh->b_face_vtx_lst,
           mesh->b_face_vtx_connect_size*sizeof(cs_lnum_t));
  }

  /* Global dimension */

  mesh_copy->n_g_cells    = mesh->n_g_cells;
  mesh_copy->n_g_i_faces  = mesh->n_g_i_faces;
  mesh_copy->n_g_b_faces  = mesh->n_g_b_faces;
  mesh_copy->n_g_vertices = mesh->n_g_vertices;

  /* Global numbering */

  if (mesh->n_g_cells != (cs_gnum_t)mesh->n_cells) {
    BFT_MALLOC(mesh_copy->global_cell_num, mesh->n_cells, cs_gnum_t);
    memcpy(mesh_copy->global_cell_num,
           mesh->global_cell_num,
           mesh->n_cells*sizeof(cs_gnum_t));
  }

  if (mesh->n_g_i_faces != (cs_gnum_t)mesh->n_i_faces) {
    BFT_MALLOC(mesh_copy->global_i_face_num, mesh->n_i_faces, cs_gnum_t);
    memcpy(mesh_copy->global_i_face_num,
           mesh->global_i_face_num,
           mesh->n_i_faces*sizeof(cs_gnum_t));
  }

  if (mesh->n_g_b_faces != (cs_gnum_t)mesh->n_b_faces) {
    BFT_MALLOC(mesh_copy->global_b_face_num, mesh->n_b_faces, cs_gnum_t);
    memcpy(mesh_copy->global_b_face_num,
           mesh->global_b_face_num,
           mesh->n_b_faces*sizeof(cs_gnum_t));
  }

  if (mesh->n_g_vertices != (cs_gnum_t)mesh->n_vertices) {
    BFT_MALLOC(mesh_copy->global_vtx_num, mesh->n_vertices, cs_gnum_t);
    memcpy(mesh_copy->global_vtx_num,
           mesh->global_vtx_num,
           mesh->n_vertices*sizeof(cs_gnum_t));
  }

  /* Parallelism and/or periodic features */
  /* and extended neighborhood features */

  mesh_copy->n_init_perio = mesh->n_init_perio;
  mesh_copy->n_transforms = mesh->n_transforms;
  mesh_copy->have_rotation_perio = mesh->have_rotation_perio;

  mesh_copy->halo_type = mesh->halo_type;

  /* modified in cs_mesh_init_halo, if necessary */
  /* ------------------------------------------- */
  /* cs_lnum_t  n_cells_with_ghosts; */
  /* cs_lnum_t  n_ghost_cells; */

  mesh_copy->n_cells_with_ghosts = mesh->n_cells_with_ghosts;
  mesh_copy->n_ghost_cells = mesh->n_ghost_cells;

  /* Re-computable connectivity features */
  /*-------------------------------------*/

  mesh_copy->n_b_cells = mesh->n_b_cells;

  BFT_MALLOC(mesh_copy->b_cells, mesh->n_b_cells, cs_lnum_t);
  memcpy(mesh_copy->b_cells, mesh->b_cells, mesh->n_b_cells*sizeof(cs_lnum_t));

  /* Group and family features */
  /*---------------------------*/

  mesh_copy->n_groups = mesh->n_groups;

  if (mesh->n_groups > 0) {
    BFT_MALLOC(mesh_copy->group_idx, mesh->n_groups + 1, cs_lnum_t);
    memcpy(mesh_copy->group_idx, mesh->group_idx,
           (mesh->n_groups + 1)*sizeof(cs_lnum_t));
    BFT_MALLOC(mesh_copy->group_lst, mesh->group_idx[mesh->n_groups] - 1, char);
    memcpy(mesh_copy->group_lst, mesh->group_lst,
           (mesh->group_idx[mesh->n_groups] - 1)*sizeof(char));
  }

  mesh_copy->n_families = mesh->n_families;
  mesh_copy->n_max_family_items = mesh->n_max_family_items;

  n_elts = mesh->n_families*mesh->n_max_family_items;
  if (n_elts > 0) {
    BFT_MALLOC(mesh_copy->family_item, n_elts , cs_lnum_t);
    memcpy(mesh_copy->family_item, mesh->family_item, n_elts*sizeof(cs_lnum_t));
  }

  BFT_MALLOC(mesh_copy->cell_family, mesh->n_cells, cs_lnum_t);
  memcpy(mesh_copy->cell_family, mesh->cell_family,
         mesh->n_cells*sizeof(cs_lnum_t));

  BFT_MALLOC(mesh_copy->i_face_family, mesh->n_i_faces, cs_lnum_t);
  memcpy(mesh_copy->i_face_family, mesh->i_face_family,
         mesh->n_i_faces*sizeof(cs_lnum_t));

  if (mesh->n_b_faces > 0) {
    BFT_MALLOC(mesh_copy->b_face_family, mesh->n_b_faces, cs_lnum_t);
    memcpy(mesh_copy->b_face_family, mesh->b_face_family,
           mesh->n_b_faces*sizeof(cs_lnum_t));
  }
}

/*----------------------------------------------------------------------------
 * Update the list of rotor vertices
 *
 * parameters:
 *   tbm <-> turbomachinery structure
 *----------------------------------------------------------------------------*/

static void
_update_rotor_vertices(cs_turbomachinery_t *tbm)
{
  cs_lnum_t  f_id, v_id;

  cs_lnum_t  n_rotor_vtx = 0;
  cs_lnum_t  *rotor_vtx = NULL;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const int  *cell_flag = tbm->cell_rotor_num;

  tbm->n_rotor_vtx = 0;
  BFT_FREE(tbm->rotor_vtx);

  BFT_MALLOC(rotor_vtx, mesh->n_vertices, cs_lnum_t);

  for (v_id = 0; v_id < mesh->n_vertices; v_id++)
    rotor_vtx[v_id] = -1;

  /* Mark from interior faces */

  for (f_id = 0; f_id < mesh->n_i_faces; f_id++) {
    cs_lnum_t c_id_0 = mesh->i_face_cells[f_id][0];
    cs_lnum_t c_id_1 = mesh->i_face_cells[f_id][1];
    if (cell_flag[c_id_0] || cell_flag[c_id_1]) {
      for (cs_lnum_t i = mesh->i_face_vtx_idx[f_id];
           i < mesh->i_face_vtx_idx[f_id+1];
           i++)
        rotor_vtx[mesh->i_face_vtx_lst[i]] = 1;
    }
  }

  /* Mark from boundary faces */

  for (f_id = 0; f_id < mesh->n_b_faces; f_id++) {
    cs_lnum_t c_id = mesh->b_face_cells[f_id];
    if (cell_flag[c_id]) {
      for (cs_lnum_t i = mesh->b_face_vtx_idx[f_id];
           i < mesh->b_face_vtx_idx[f_id+1];
           i++)
        rotor_vtx[mesh->b_face_vtx_lst[i]] = 1;
    }
  }

  /* Now transform vertex flag to list */

  n_rotor_vtx = 0;
  for (v_id = 0; v_id < mesh->n_vertices; v_id++) {
    if (rotor_vtx[v_id] > -1)
      rotor_vtx[n_rotor_vtx++] = v_id;
  }
  BFT_REALLOC(rotor_vtx, n_rotor_vtx, cs_lnum_t);

  /* Update strcuture */

  tbm->n_rotor_vtx = n_rotor_vtx;
  tbm->rotor_vtx = rotor_vtx;
}

/*----------------------------------------------------------------------------
 * Update mesh vertex positions
 *
 * parameters:
 *   mesh <-> mesh to update
 *   t    <-- associated time
 *----------------------------------------------------------------------------*/

static void
_update_geometry(cs_mesh_t  *mesh,
                 cs_real_t   t)
{
  double   matrix[3][4];

  cs_turbomachinery_t *tbm = cs_glob_turbomachinery;

  tbm->angle = tbm->omega * t;

  _rotation_matrix(tbm->angle,
                   tbm->rotation_axis,
                   tbm->rotation_invariant,
                   matrix);

  for (cs_lnum_t i = 0; i < tbm->n_rotor_vtx; i++) {
    cs_lnum_t vtx_id = tbm->rotor_vtx[i];
    _apply_vector_transfo(matrix, &(mesh->vtx_coord[3*vtx_id]));
  }
}

/*----------------------------------------------------------------------------
 * Select rotor cells
 *
 * parameters:
 *   tbm <-> turbomachinery options structure
 *----------------------------------------------------------------------------*/

static void
_select_rotor_cells(cs_turbomachinery_t  *tbm)
{
  cs_lnum_t _n_cells = 0;
  cs_lnum_t *_cell_list = NULL;

  cs_mesh_t *m = cs_glob_mesh;

  assert(tbm->rotor_cells_c != NULL);

  BFT_REALLOC(tbm->cell_rotor_num, m->n_cells_with_ghosts, int);

  for (cs_lnum_t i = 0; i < m->n_cells_with_ghosts; i++)
    tbm->cell_rotor_num[i] = 0;

  BFT_MALLOC(_cell_list, m->n_cells_with_ghosts, cs_lnum_t);

  cs_selector_get_cell_list(tbm->rotor_cells_c,
                            &_n_cells,
                            _cell_list);

  for (cs_lnum_t i = 0; i < _n_cells; i++)
    tbm->cell_rotor_num[_cell_list[i] - 1] = 1;

  BFT_FREE(_cell_list);
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------
 * Assist map turbomachinery to fortran module
 *
 * parameters:
 *   iturbo <-- turbomachinery type flag
 *   rotax  <-- axis of rotation
 *   irotce <-- pointer to cell flag for rotation
 *----------------------------------------------------------------------------*/

void
cs_f_map_turbomachinery_module(cs_int_t    *iturbo,
                               cs_real_t    rotax[],
                               cs_int_t   **irotce)
{
  if (cs_glob_turbomachinery != NULL) {

    *iturbo = cs_glob_turbomachinery->model;

    for (int i = 0; i < 3; i++)
      rotax[i] =   cs_glob_turbomachinery->rotation_axis[i]
                 * cs_glob_turbomachinery->omega;

    /* Assign rotor cells flag array to module */

    *irotce = cs_glob_turbomachinery->cell_rotor_num;

  }
  else {

    *iturbo = CS_TURBOMACHINERY_NONE;

    for (int i = 0; i < 3; i++)
      rotax[i] = 0.;

    *irotce = NULL;

  }
}

/*! \endcond (end ignore by Doxygen) */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define rotor/stator model.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_set_model(cs_turbomachinery_model_t  model)
{
  if (   model == CS_TURBOMACHINERY_NONE
      && cs_glob_turbomachinery != NULL) {
    cs_turbomachinery_finalize();
    return;
  }

  else if (cs_glob_turbomachinery == NULL)
    cs_glob_turbomachinery = _turbomachinery_create();

  cs_glob_turbomachinery->model = model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return rotor/stator model.
 */
/*----------------------------------------------------------------------------*/

cs_turbomachinery_model_t
cs_turbomachinery_get_model(void)
{
  if (cs_glob_turbomachinery == NULL)
   return CS_TURBOMACHINERY_NONE;
  else
    return cs_glob_turbomachinery->model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a rotor by its axis and cell selection criteria.
 *
 * \param[in]  cell_criteria       cell selection criteria string
 * \param[in]  rotation_velocity   rotation velocity, in radians/second
 * \param[in]  rotation_axis       rotation axis vector
 * \param[in]  rotation_invariant  rotation invariant point
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_add_rotor(const char    *cell_criteria,
                            double         rotation_velocity,
                            const double   rotation_axis[3],
                            const double   rotation_invariant[3])
{
  static int n_calls = 0;

  if (n_calls > 0)
    bft_error(__FILE__, __LINE__, 0,
              "%s: only one rotor may be defined in the current version.",
              __func__);

  n_calls += 1;

  cs_turbomachinery_t *tbm = cs_glob_turbomachinery;
  if (tbm == NULL)
    return;

  const double *v = rotation_axis;
  double len = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

  tbm->omega = rotation_velocity;
  for (int i = 0; i < 3; i++) {
    tbm->rotation_axis[i] = v[i]/len;
    tbm->rotation_invariant[i] = rotation_invariant[i];
  }

  BFT_REALLOC(tbm->rotor_cells_c, strlen(cell_criteria) + 1, char);
  strcpy(tbm->rotor_cells_c, cell_criteria);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a cs_join_t structure to the list of rotor/stator joinings.
 *
 * \param[in]  sel_criteria   boundary face selection criteria
 * \param[in]  fraction       value of the fraction parameter
 * \param[in]  plane          value of the plane parameter
 * \param[in]  verbosity      level of verbosity required
 * \param[in]  visualization  level of visualization required
 *
 * \return  number (1 to n) associated with new joining
 */
/*----------------------------------------------------------------------------*/

int
cs_turbomachinery_join_add(const char  *sel_criteria,
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
                     false);

  cs_glob_join_count++; /* Store number of joining (without periodic ones) */
  cs_glob_n_joinings++;

  return cs_glob_n_joinings;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update mesh for unsteady rotor/stator computation.
 *
 * \param[in]   t_cur_mob  current rotor time
 * \param[out]  t_elapsed  elapsed computation time
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_update_mesh(double   t_cur_mob,
                              double  *t_elapsed)
{
  double  t_start, t_end;

  cs_halo_type_t halo_type = cs_glob_mesh->halo_type;

  t_start = cs_timer_wtime();

  /* Indicates we are in the framework of turbomachinery */

  cs_glob_turbomachinery->active = true;

  /* Destroy previous global mesh and related entities */

  cs_mesh_location_finalize();
  cs_mesh_quantities_destroy(cs_glob_mesh_quantities);
  cs_mesh_destroy(cs_glob_mesh);

  /* Create new global mesh and related entities */

  cs_mesh_location_initialize();
  cs_glob_mesh = cs_mesh_create();
  cs_glob_mesh->verbosity = 0;
  _copy_mesh(cs_glob_turbomachinery->reference_mesh, cs_glob_mesh);
  cs_glob_mesh_builder = cs_mesh_builder_create();
  cs_glob_mesh_quantities = cs_mesh_quantities_create();

  /* Update geometry, if necessary */

  if (cs_glob_turbomachinery->rotor_vtx != NULL)
    _update_geometry(cs_glob_mesh, t_cur_mob);

  /* Reset the interior faces -> cells connectivity */
  /* (in order to properly build the halo of the joined mesh) */

  cs_mesh_to_builder_perio_faces(cs_glob_mesh, cs_glob_mesh_builder);

  {
    int i;
    cs_lnum_t f_id;
    cs_lnum_2_t *i_face_cells = (cs_lnum_2_t *)cs_glob_mesh->i_face_cells;
    const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
    for (f_id = 0; f_id < cs_glob_mesh->n_i_faces; f_id++) {
      for (i = 0; i < 2; i++) {
        if (i_face_cells[f_id][i] > n_cells)
          i_face_cells[f_id][i] = 0;
      }
    }
  }

  /* Join meshes and build periodicity links */

  cs_join_all(false);

  cs_lnum_t boundary_changed = 0;
  if (cs_glob_turbomachinery->n_b_faces_ref > -1) {
    if (cs_glob_mesh->n_b_faces != cs_glob_turbomachinery->n_b_faces_ref)
      boundary_changed = 1;
  }
  cs_parall_counter_max(&boundary_changed, 1);

  /* Check that joining has not added or removed boundary faces.
     Postprocess new faces appearing on boundary or inside of mesh:
     this assumes that joining appends new faces at the end of the mesh */

  if (boundary_changed) {
    const int writer_id = -2;
    const int writer_ids[] = {writer_id};
    const int mesh_id = cs_post_get_free_mesh_id();
    cs_lnum_t b_face_count[] = {cs_glob_turbomachinery->n_b_faces_ref,
                                cs_glob_mesh->n_b_faces};
    cs_gnum_t n_g_b_faces_ref = cs_glob_turbomachinery->n_b_faces_ref;
    cs_parall_counter(&n_g_b_faces_ref, 1);
    cs_post_init_error_writer();
    cs_post_define_surface_mesh_by_func(mesh_id,
                                        _("Added boundary faces"),
                                        NULL,
                                        _post_error_faces_select,
                                        NULL,
                                        b_face_count,
                                        false, /* time varying */
                                        true,  /* add groups if present */
                                        false, /* auto variables */
                                        1,
                                        writer_ids);
    cs_post_activate_writer(writer_id, 1);
    cs_post_write_meshes(NULL);
    bft_error(__FILE__, __LINE__, 0,
              _("Error in turbomachinery mesh update:\n"
                "Number of boundary faces has changed from %llu to %llu.\n"
                "There are probably unjoined faces, "
                "due to an insufficiently regular mesh;\n"
                "adjusting mesh joining parameters might help."),
              (unsigned long long)n_g_b_faces_ref,
              (unsigned long long)cs_glob_mesh->n_g_b_faces);
  }

  cs_glob_turbomachinery->n_b_faces_ref = cs_glob_mesh->n_b_faces;

  /* Initialize extended connectivity, ghost cells and other remaining
     parallelism-related structures */

  cs_mesh_init_halo(cs_glob_mesh, cs_glob_mesh_builder, halo_type);
  cs_mesh_update_auxiliary(cs_glob_mesh);

  /* Destroy the temporary structure used to build the main mesh */

  cs_mesh_builder_destroy(&cs_glob_mesh_builder);

  /* Update interior faces renumbering based on code options */

  cs_renumber_i_faces(cs_glob_mesh);

  /* Build group classes */

  cs_mesh_init_group_classes(cs_glob_mesh);

  /* Print info on mesh */

  if (cs_glob_mesh->verbosity > 0)
    cs_mesh_print_info(cs_glob_mesh, _("Mesh"));

  /* Compute geometric quantities related to the mesh */

  cs_mesh_quantities_compute(cs_glob_mesh, cs_glob_mesh_quantities);

  /* Initialize selectors and locations for the mesh */

  cs_mesh_init_selectors();
  cs_mesh_location_build(cs_glob_mesh, -1);

  /* Check coherency if debugging */

#if 0
  cs_mesh_coherency_check();
#endif

  /* Update Fortran mesh sizes and quantities */

  cs_preprocess_mesh_update_fortran();

  /* Update rotor cells flag array in case of parallelism and/or periodicity */

  if (cs_glob_mesh->halo != NULL) {

    const cs_mesh_t *m = cs_glob_mesh;

    BFT_REALLOC(cs_glob_turbomachinery->cell_rotor_num,
                m->n_cells_with_ghosts,
                int);

    cs_halo_sync_untyped(m->halo,
                         CS_HALO_EXTENDED,
                         sizeof(int),
                         cs_glob_turbomachinery->cell_rotor_num);

  }

  /* Update linear algebra APIs relative to mesh */

  cs_gradient_perio_update_mesh();
  cs_matrix_update_mesh();

  t_end = cs_timer_wtime();

  *t_elapsed = t_end - t_start;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initializations for turbomachinery computation.
 *
 * \note This function should be called before once the mesh is built,
 *       but before cs_post_init_meshes() so that postprocessing meshes are
 *       updated correctly in the transient case.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_initialize(void)
{
  /* Define model; could be moved anywhere before time loop. */

  cs_gui_turbomachinery();
  cs_user_turbomachinery();

  if (cs_glob_turbomachinery == NULL)
    return;

  cs_turbomachinery_t *tbm = cs_glob_turbomachinery;

  if (tbm->model == CS_TURBOMACHINERY_NONE)
    return;

  /* Define rotor(s) */

  cs_gui_turbomachinery_rotor();
  cs_user_turbomachinery_rotor();
  _select_rotor_cells(tbm);

  /* Build the reference mesh that duplicates the global mesh before joining */

  _copy_mesh(cs_glob_mesh, tbm->reference_mesh);

  /* Build the rotor vertices list, if required */

  if (tbm->model == CS_TURBOMACHINERY_TRANSIENT) {

    if (cs_glob_mesh->halo != NULL)
      cs_halo_sync_num(cs_glob_mesh->halo,
                       CS_HALO_EXTENDED,
                       tbm->cell_rotor_num);

    _update_rotor_vertices(tbm);

  }

  /* Complete the mesh with rotor-stator joining */

  if (cs_glob_n_joinings > 0) {
    cs_real_t t_elapsed;
    cs_turbomachinery_update_mesh(0.,&t_elapsed);
  }

  /* Adapt postprocessing options if required;
     must be called before cs_post_init_meshes(). */

  if (tbm->model == CS_TURBOMACHINERY_TRANSIENT)
    cs_post_set_changing_connectivity();

  /* Destroy the reference mesh, if required */

  if (tbm->model == CS_TURBOMACHINERY_FROZEN) {
    cs_mesh_destroy(tbm->reference_mesh);
    tbm->reference_mesh = NULL;
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free turbomachinery structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_finalize(void)
{
  if (cs_glob_turbomachinery != NULL) {

    cs_turbomachinery_t *tbm = cs_glob_turbomachinery;

    BFT_FREE(tbm->rotor_cells_c);

    BFT_FREE(tbm->cell_rotor_num);

    if (tbm->rotor_vtx != NULL)
      BFT_FREE(tbm->rotor_vtx);

    if (tbm->reference_mesh != NULL)
      cs_mesh_destroy(tbm->reference_mesh);

  }

  BFT_FREE(cs_glob_turbomachinery);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reinitialize interior face-based fields.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_reinit_i_face_fields(void)
{
  const int n_fields = cs_field_n_fields();

  for (int f_id = 0; f_id < n_fields; f_id++) {

    cs_field_t *f = cs_field_by_id(f_id);

    if (   cs_mesh_location_get_type(f->location_id)
        == CS_MESH_LOCATION_INTERIOR_FACES)
      cs_field_allocate_values(f);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Resize cell-based fields.
 *
 * This function only handles fields owning their values.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_resize_cell_fields(void)
{
  const int n_fields = cs_field_n_fields();

  const cs_halo_t *halo = cs_glob_mesh->halo;
  const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(CS_MESH_LOCATION_CELLS);
  cs_lnum_t _n_cells = n_elts[2];

  for (int f_id = 0; f_id < n_fields; f_id++) {

    cs_field_t *f = cs_field_by_id(f_id);

    if (   f->location_id == CS_MESH_LOCATION_CELLS
        && f->is_owner) {

      if (f->dim > 1 && f->interleaved == false)
        bft_error(__FILE__, __LINE__, 0,
                  "%s: fields owning their values (i.e. not mapped)\n"
                  "should be interleaved, but \"%s\" is not.",
                  __func__, f->name);

      /* Ghost cell sizes may change, but the number of main cells
         is unchanged, so a simple reallocation will do */

      for (int kk = 0; kk < f->n_time_vals; kk++) {

        BFT_REALLOC(f->vals[kk], _n_cells*f->dim, cs_real_t);

        if (halo != NULL) {

          cs_halo_sync_untyped(halo,
                               CS_HALO_EXTENDED,
                               f->dim*sizeof(cs_real_t),
                               f->vals[kk]);
          if (f->dim == 3)
            cs_halo_perio_sync_var_vect(halo,
                                        CS_HALO_EXTENDED,
                                        f->vals[kk],
                                        f->dim);
        }
      }

      f->val = f->vals[0];
      if (f->n_time_vals > 1) f->val_pre = f->vals[1];

    }

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute rotation matrix
 *
 * \param[in]   rotor_num rotor number (1 to n numbering)
 * \param[in]   theta  rotation angle, in radians
 * \param[out]  matrix resulting rotation matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_rotation_matrix(int        rotor_num,
                                  double     theta,
                                  cs_real_t  matrix[3][4])
{
  cs_turbomachinery_t *tbm = cs_glob_turbomachinery;

  if (rotor_num != 1)
    bft_error(__FILE__, __LINE__, 0,
              "%s: only one rotor may be used in the current version.",
              __func__);

  _rotation_matrix(theta,
                   tbm->rotation_axis,
                   tbm->rotation_invariant,
                   matrix);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return cell rotor number.
 *
 * Each cell may be associated with a given rotor, or rotation, with 0
 * indicating that that cell does not rotate.
 *
 * \return  array defining rotor number associated with each cell
 *          (0 for none, 1 to n otherwise)
 */
/*----------------------------------------------------------------------------*/

const int *
cs_turbomachinery_get_cell_rotor_num(void)
{
  return cs_glob_turbomachinery->cell_rotor_num;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return rotation velocity
 *
 * \param[in]   rotor_num rotor number (1 to n numbering)
 */
/*----------------------------------------------------------------------------*/

double
cs_turbomachinery_get_rotation_velocity(int rotor_num)
{
  if (rotor_num != 1)
    bft_error(__FILE__, __LINE__, 0,
              "%s: only one rotor may be used in the current version.",
              __func__);

  return cs_glob_turbomachinery->omega;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Rotation of vector and tensor fields.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_rotate_fields(const cs_real_t dt[])
{
  cs_lnum_t i;
  double matrix[3][4];

  cs_turbomachinery_t *tbm = cs_glob_turbomachinery;
  cs_real_t time_step = dt[0];

  int n_fields = cs_field_n_fields();

  _rotation_matrix(tbm->omega*time_step,
                   tbm->rotation_axis,
                   tbm->rotation_invariant,
                   matrix);

  for (int f_id = 0; f_id < n_fields; f_id++) {

    cs_field_t *f = cs_field_by_id(f_id);

    if (! (f->dim > 1 && (f->type & CS_FIELD_VARIABLE)))
      continue;

    const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(f->location_id);
    cs_lnum_t _n_elts = n_elts[2];

    if (f->dim == 3) {
      if (f->interleaved == false) {
        double v[3];
        for (i = 0; i < _n_elts; i++) {
          v[0] = f->val[i];
          v[1] = f->val[i + _n_elts];
          v[2] = f->val[i + _n_elts*2];
          _apply_vector_rotation(matrix, v);
          f->val[i]             = v[0];
          f->val[i + _n_elts]   = v[1];
          f->val[i + _n_elts*2] = v[2];
        }
      }
      else {
        for (i = 0; i < _n_elts; i++)
          _apply_vector_rotation(matrix, f->val + i*3);
      }
    }

    else if (f->dim == 6) {
      assert(f->interleaved);
        for (i = 0; i < _n_elts; i++)
          _apply_sym_tensor_rotation(matrix, f->val + i*6);
    }

    assert(f->dim == 3 || f->dim == 6);

  } /* End of loop on fields */

  /* Specific handling of Reynolds stresses */

  cs_field_t  *fr11 = cs_field_by_name("r11");
  if (fr11 != NULL) {

    cs_field_t  *fr22 = cs_field_by_name("r22");
    cs_field_t  *fr33 = cs_field_by_name("r33");
    cs_field_t  *fr12 = cs_field_by_name("r12");
    cs_field_t  *fr13 = cs_field_by_name("r13");
    cs_field_t  *fr23 = cs_field_by_name("r23");

    const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(fr11->location_id);
    cs_lnum_t _n_elts = n_elts[2];

    for (i = 0; i < _n_elts; i++) {
      double t[6];
      t[0] = fr11->val[i];
      t[1] = fr22->val[i];
      t[2] = fr33->val[i];
      t[3] = fr12->val[i];
      t[4] = fr13->val[i];
      t[5] = fr23->val[i];
      _apply_sym_tensor_rotation(matrix, t);
      fr11->val[i] = t[0];
      fr22->val[i] = t[1];
      fr33->val[i] = t[2];
      fr12->val[i] = t[3];
      fr13->val[i] = t[4];
      fr23->val[i] = t[5];
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute velocity relative to fixed coordinates at a given point
 *
 * \param[in]   rotor_num  rotor number (1 to n numbering)
 * \param[in]   coords     point coordinates
 * \param[out]  velocity   velocity relative to fixed coordinates
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_relative_velocity(int              rotor_num,
                                    const cs_real_t  coords[3],
                                    cs_real_t        velocity[3])
{
  if (rotor_num != 1)
    bft_error(__FILE__, __LINE__, 0,
              "%s: only one rotor may be used in the current version.",
              __func__);

  const cs_turbomachinery_t *tbm = cs_glob_turbomachinery;

  _relative_velocity(tbm->omega,
                     tbm->rotation_axis,
                     tbm->rotation_invariant,
                     coords,
                     velocity);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
