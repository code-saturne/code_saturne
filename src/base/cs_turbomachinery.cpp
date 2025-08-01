/*============================================================================
 * Turbomachinery modeling features.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_printf.h"

#include "fvm/fvm_selector.h"

#include "base/cs_interface.h"

#include "base/cs_base.h"
#include "base/cs_boundary_zone.h"
#include "base/cs_coupling.h"
#include "alge/cs_cell_to_vertex.h"
#include "base/cs_ext_neighborhood.h"
#include "base/cs_field.h"
#include "base/cs_function.h"
#include "alge/cs_gradient.h"
#include "gui/cs_gui.h"
#include "gui/cs_gui_mesh.h"
#include "gui/cs_gui_output.h"
#include "alge/cs_gradient.h"
#include "mesh/cs_join.h"
#include "base/cs_halo.h"
#include "base/cs_halo_perio.h"
#include "alge/cs_matrix_default.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_adjacencies.h"
#include "mesh/cs_mesh_coherency.h"
#include "mesh/cs_mesh_location.h"
#include "mesh/cs_mesh_quantities.h"
#include "mesh/cs_mesh_to_builder.h"
#include "alge/cs_multigrid.h"
#include "base/cs_parall.h"
#include "base/cs_post.h"
#include "base/cs_preprocess.h"
#include "base/cs_prototypes.h"
#include "base/cs_renumber.h"
#include "base/cs_rotation.h"
#include "base/cs_time_step.h"
#include "base/cs_timer.h"
#include "base/cs_timer_stats.h"
#include "base/cs_restart.h"
#include "base/cs_sat_coupling.h"
#include "base/cs_preprocessor_data.h"
#include "base/cs_volume_zone.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_turbomachinery.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/* Turbomachinery structure */

typedef struct {

  cs_turbomachinery_model_t  model;             /* turbomachinery model type */

  int                        n_rotors;          /* number of rotors */

  int                        n_couplings;       /* number of couplings */

  cs_rotation_t             *rotation;          /* rotation structures */
  char                     **rotor_cells_c;     /* rotor cells selection
                                                   criteria (for each rotor) */

  int                        nt_cur;            /* current time step */
  int                        n_max_join_tries;  /* maximum number of tries
                                                   for joining differences */
  double                     dt_retry;          /* time shift multiplier for
                                                   retry position */
  double                     t_cur;             /* current time for update */
  double                     t_cur_kc;          /* Kahan compensation for
                                                   current time for update */

  cs_mesh_t                 *reference_mesh;    /* reference mesh (before
                                                   rotation and joining) */

  cs_lnum_t                  n_b_faces_ref;     /* reference number of
                                                   boundary faces */

  int                       *cell_rotor_num;    /* cell rotation axis number */

  cs_real_t                 *coftur;            /* wall BC coefficient to update
                                                   wall velocity after geometry
                                                   update */

  cs_real_t                 *hfltur;            /* wall exchange coefficient
                                                   for update after that
                                                   of geometry */

  bool active;

} cs_turbomachinery_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

cs_turbomachinery_t  *_turbomachinery = nullptr;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void cs_f_map_turbomachinery_model(int  *iturbo);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update Kahan compensation.
 *
 * Called as a function to help avoid compiler optimizating Kahan
 * summation out (though aggressive or LTO optimizations might).
 *
 * new compensation term is (a - b) - c
 *
 * \param[in]  a  a
 * \param[in]  b  b
 * \param[in]  c  c
 */
/*----------------------------------------------------------------------------*/

static void
_update_kahan_compensation(double  *a,
                           double  *b,
                           double  *c)
{
  _turbomachinery->t_cur_kc = (*a - *b) - *c;
}

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
  cs_lnum_t *_face_ids = nullptr;

  const cs_lnum_t  *count = reinterpret_cast<cs_lnum_t *>(input);

  CS_MALLOC(_face_ids, count[1], cs_lnum_t);

  for (face_id = count[0]; face_id < count[1]; face_id++)
    _face_ids[_n_faces++] = face_id;

  *n_faces = _n_faces;
  *face_ids = _face_ids;
}

/*----------------------------------------------------------------------------
 * Function for tagging of coupling mesh.
 *
 * We use the rotor number as a tag
 *
 * Note: if the context pointer is non-null, it must point to valid data
 * when the selection function is called, so that value or structure
 * should not be temporary (i.e. local);
 *
 * parameters:
 *   context         <-> pointer to optional (untyped) value or structure.
 *   mesh            <-> nodal mesh which should be tagged
 *   n_points        <-- number of points to tag
 *   point_list_base <-- base numbering for point_list
 *   point_list      <-- optional indirection for points
 *   point_tag       --> point tag values (size: n_tags)
 *----------------------------------------------------------------------------*/

static void
_turbomachinery_coupling_tag(void            *context,
                             fvm_nodal_t     *mesh,
                             cs_lnum_t        n_points,
                             cs_lnum_t        point_list_base,
                             const cs_lnum_t  point_list[],
                             int             *point_tag)
{
  const cs_turbomachinery_t *tbm
    = reinterpret_cast<cs_turbomachinery_t *>(context);
  const cs_mesh_t *m = cs_glob_mesh;

  /* Tag elements (boundary faces) */

  if (mesh != nullptr) {

    cs_lnum_t n_elts = fvm_nodal_get_n_entities(mesh, 3);

    const int ent_dim = 3;
    int *elt_tag;
    cs_lnum_t *parent_num;

    CS_MALLOC(elt_tag, n_elts, int);
    CS_MALLOC(parent_num, n_elts, cs_lnum_t);

    fvm_nodal_get_parent_num(mesh, ent_dim, parent_num);
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t c_id = parent_num[i] - 1;
      elt_tag[i] = tbm->cell_rotor_num[c_id];
    }

    CS_FREE(parent_num);

    fvm_nodal_set_tag(mesh, elt_tag, ent_dim);

    CS_FREE(elt_tag);

  }

  /* Tag vertices */

  if (point_list != nullptr) {
    for (cs_lnum_t i = 0; i < n_points; i++) {
      cs_lnum_t f_id = point_list[i] - point_list_base;
      cs_lnum_t c_id = m->b_face_cells[f_id];
      point_tag[i] = tbm->cell_rotor_num[c_id];
    }
  }
  else {
    for (cs_lnum_t i = 0; i < n_points; i++) {
      cs_lnum_t c_id = m->b_face_cells[i];
      point_tag[i] = tbm->cell_rotor_num[c_id];
    }
  }
}

/*----------------------------------------------------------------------------
 * Create an empty turbomachinery structure
 *----------------------------------------------------------------------------*/

static cs_turbomachinery_t *
_turbomachinery_create(void)
{
  cs_turbomachinery_t  *tbm = nullptr;

  CS_MALLOC(tbm, 1, cs_turbomachinery_t);

  tbm->n_rotors = 0;
  tbm->rotor_cells_c = nullptr;

  CS_MALLOC(tbm->rotation, 1, cs_rotation_t); /* Null rotation at id 0 */
  cs_rotation_t *r = tbm->rotation;
  r->omega = 0;
  r->angle = 0;
  for (int i = 0; i < 3; i++) {
    r->axis[i] = 0;
    r->invariant[i] = 0;
  }

  tbm->nt_cur = 0;
  tbm->n_max_join_tries = 5;
  tbm->t_cur = 0;
  tbm->t_cur_kc = 0;
  tbm->dt_retry = 1e-2;

  tbm->reference_mesh = cs_mesh_create();
  tbm->n_b_faces_ref = -1;
  tbm->coftur = nullptr;
  tbm->hfltur = nullptr;
  tbm->cell_rotor_num = nullptr;
  tbm->model = CS_TURBOMACHINERY_NONE;
  tbm->n_couplings = 0;

  return tbm;
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
  double  t[3] = {v[0], v[1], v[2]};

  for (int i = 0; i < 3; i++)
    v[i] = m[i][0]*t[0] + m[i][1]*t[1] + m[i][2]*t[2];
}

/*----------------------------------------------------------------------------
 * Compute a matrix * tensor * Tmatrix product to apply a rotation to a
 * given symmetric tensor
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
 * In case of coupling only (i.e. no joining), only a partial copy of the
 * mesh is done.
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
  mesh_copy->time_dep   = mesh->time_dep;

  /* Local dimensions */

  mesh_copy->n_cells    = mesh->n_cells;
  mesh_copy->n_i_faces  = mesh->n_i_faces;
  mesh_copy->n_b_faces  = mesh->n_b_faces;
  mesh_copy->n_vertices = mesh->n_vertices;

  mesh_copy->i_face_vtx_connect_size = mesh->i_face_vtx_connect_size;
  mesh_copy->b_face_vtx_connect_size = mesh->b_face_vtx_connect_size;

  /* Local structures */

  CS_REALLOC(mesh_copy->vtx_coord, (3*mesh->n_vertices), cs_real_t);
  memcpy(mesh_copy->vtx_coord,
         mesh->vtx_coord,
         3*mesh->n_vertices*sizeof(cs_real_t));

  if (cs_glob_n_joinings < 1)
    return;

  CS_MALLOC(mesh_copy->i_face_cells, mesh->n_i_faces, cs_lnum_2_t);
  memcpy(mesh_copy->i_face_cells,
         mesh->i_face_cells,
         mesh->n_i_faces*sizeof(cs_lnum_2_t));

  if (mesh->n_b_faces > 0) {
    CS_MALLOC(mesh_copy->b_face_cells, mesh->n_b_faces, cs_lnum_t);
    memcpy(mesh_copy->b_face_cells,
           mesh->b_face_cells,
           mesh->n_b_faces*sizeof(cs_lnum_t));
  }

  CS_MALLOC(mesh_copy->i_face_vtx_idx, mesh->n_i_faces + 1, cs_lnum_t);
  memcpy(mesh_copy->i_face_vtx_idx,
         mesh->i_face_vtx_idx,
         (mesh->n_i_faces + 1)*sizeof(cs_lnum_t));

  CS_MALLOC(mesh_copy->i_face_vtx_lst,
            mesh->i_face_vtx_connect_size,
            cs_lnum_t);
  memcpy(mesh_copy->i_face_vtx_lst, mesh->i_face_vtx_lst,
         mesh->i_face_vtx_connect_size*sizeof(cs_lnum_t));

  CS_MALLOC(mesh_copy->b_face_vtx_idx, mesh->n_b_faces + 1, cs_lnum_t);
  memcpy(mesh_copy->b_face_vtx_idx,
         mesh->b_face_vtx_idx,
         (mesh->n_b_faces + 1)*sizeof(cs_lnum_t));

  if (mesh->b_face_vtx_connect_size > 0) {
    CS_MALLOC(mesh_copy->b_face_vtx_lst,
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

  if (mesh->global_cell_num != nullptr) {
    CS_MALLOC(mesh_copy->global_cell_num, mesh->n_cells, cs_gnum_t);
    memcpy(mesh_copy->global_cell_num,
           mesh->global_cell_num,
           mesh->n_cells*sizeof(cs_gnum_t));
  }

  if (mesh->global_i_face_num != nullptr) {
    CS_MALLOC(mesh_copy->global_i_face_num, mesh->n_i_faces, cs_gnum_t);
    memcpy(mesh_copy->global_i_face_num,
           mesh->global_i_face_num,
           mesh->n_i_faces*sizeof(cs_gnum_t));
  }

  if (mesh->global_b_face_num != nullptr) {
    CS_MALLOC(mesh_copy->global_b_face_num, mesh->n_b_faces, cs_gnum_t);
    memcpy(mesh_copy->global_b_face_num,
           mesh->global_b_face_num,
           mesh->n_b_faces*sizeof(cs_gnum_t));
  }

  if (mesh->global_vtx_num != nullptr) {
    CS_MALLOC(mesh_copy->global_vtx_num, mesh->n_vertices, cs_gnum_t);
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

  CS_MALLOC(mesh_copy->b_cells, mesh->n_b_cells, cs_lnum_t);
  memcpy(mesh_copy->b_cells, mesh->b_cells, mesh->n_b_cells*sizeof(cs_lnum_t));

  /* Group and family features */
  /*---------------------------*/

  mesh_copy->n_groups = mesh->n_groups;

  if (mesh->n_groups > 0) {
    CS_MALLOC(mesh_copy->group_idx, mesh->n_groups + 1, int);
    memcpy(mesh_copy->group_idx, mesh->group_idx,
           (mesh->n_groups + 1)*sizeof(cs_lnum_t));
    CS_MALLOC(mesh_copy->group, mesh->group_idx[mesh->n_groups], char);
    memcpy(mesh_copy->group, mesh->group,
           mesh->group_idx[mesh->n_groups]*sizeof(char));
  }

  mesh_copy->n_families = mesh->n_families;
  mesh_copy->n_max_family_items = mesh->n_max_family_items;

  n_elts = mesh->n_families*mesh->n_max_family_items;
  if (n_elts > 0) {
    CS_MALLOC(mesh_copy->family_item, n_elts , int);
    memcpy(mesh_copy->family_item, mesh->family_item, n_elts*sizeof(int));
  }

  CS_MALLOC(mesh_copy->cell_family, mesh->n_cells_with_ghosts, int);
  memcpy(mesh_copy->cell_family, mesh->cell_family,
         mesh->n_cells_with_ghosts*sizeof(int));

  CS_MALLOC(mesh_copy->i_face_family, mesh->n_i_faces, int);
  memcpy(mesh_copy->i_face_family, mesh->i_face_family,
         mesh->n_i_faces*sizeof(int));

  if (mesh->n_b_faces > 0) {
    CS_MALLOC(mesh_copy->b_face_family, mesh->n_b_faces, int);
    memcpy(mesh_copy->b_face_family, mesh->b_face_family,
           mesh->n_b_faces*sizeof(int));
  }

  if (mesh->i_face_r_gen != nullptr) {
    CS_MALLOC(mesh_copy->i_face_r_gen, mesh->n_i_faces, char);
    memcpy(mesh_copy->i_face_r_gen, mesh->i_face_r_gen,
           mesh->n_i_faces);
  }
  if (mesh->vtx_r_gen != nullptr) {
    CS_MALLOC(mesh_copy->vtx_r_gen, mesh->n_vertices, char);
    memcpy(mesh_copy->vtx_r_gen, mesh->vtx_r_gen, mesh->n_vertices);
  }
}

/*----------------------------------------------------------------------------
 * Update mesh vertex positions
 *
 * parameters:
 *   t    <-- associated time
 *----------------------------------------------------------------------------*/

static void
_update_angle(const cs_time_step_t  *ts)
{
  cs_turbomachinery_t *tbm = _turbomachinery;

  /* Now update coordinates */

  if (tbm->nt_cur < ts->nt_cur) {

    tbm->nt_cur = ts->nt_cur;

    if (ts->nt_cur > ts->nt_prev) {

      double dt = ts->dt[0];

      for (int j = 0; j < tbm->n_rotors+1; j++) {
        cs_rotation_t *r = tbm->rotation + j;
        r->angle += r->omega * dt;
        /* To avoid degrading numerical precision after many rotations,
           keep angle value in [0, 2pi] range since it is periodic. */
        r->angle = fmod(r->angle, atan(1.)*8.);
      }

      double z = dt - tbm->t_cur_kc;
      double t = tbm->t_cur + z;

      _update_kahan_compensation(&t, &(tbm->t_cur), &z);

      tbm->t_cur = t;
    }

  }
}

/*----------------------------------------------------------------------------
 * Update mesh vertex positions
 *
 * parameters:
 *   mesh <-> mesh to update
 *   dt   <-- associated time delta (0 for current, unmodified time)
 *----------------------------------------------------------------------------*/

static void
_update_geometry(cs_mesh_t  *mesh,
                 cs_real_t   dt)
{
  cs_turbomachinery_t *tbm = _turbomachinery;

  cs_lnum_t  f_id, v_id;

  cs_lnum_t  *vtx_rotor_num = nullptr;

  const int  *cell_flag = tbm->cell_rotor_num;

  CS_MALLOC(vtx_rotor_num, mesh->n_vertices, cs_lnum_t);

  for (v_id = 0; v_id < mesh->n_vertices; v_id++)
    vtx_rotor_num[v_id] = 0;

  /* Mark from interior faces */

  for (f_id = 0; f_id < mesh->n_i_faces; f_id++) {
    cs_lnum_t c_id_0 = mesh->i_face_cells[f_id][0];
    cs_lnum_t c_id_1 = mesh->i_face_cells[f_id][1];
    assert(c_id_0 > -1);
    assert(c_id_1 > -1);
    if (c_id_0 < mesh->n_cells && cell_flag[c_id_0] != 0) {
      for (cs_lnum_t i = mesh->i_face_vtx_idx[f_id];
           i < mesh->i_face_vtx_idx[f_id+1];
           i++)
        vtx_rotor_num[mesh->i_face_vtx_lst[i]] = cell_flag[c_id_0];
    }
    else if (c_id_1 < mesh->n_cells && cell_flag[c_id_1] != 0) {
      for (cs_lnum_t i = mesh->i_face_vtx_idx[f_id];
           i < mesh->i_face_vtx_idx[f_id+1];
           i++)
        vtx_rotor_num[mesh->i_face_vtx_lst[i]] = cell_flag[c_id_1];
    }
  }

  /* Mark from boundary faces */

  for (f_id = 0; f_id < mesh->n_b_faces; f_id++) {
    cs_lnum_t c_id = mesh->b_face_cells[f_id];
    if (cell_flag[c_id] != 0) {
      for (cs_lnum_t i = mesh->b_face_vtx_idx[f_id];
           i < mesh->b_face_vtx_idx[f_id+1];
           i++)
        vtx_rotor_num[mesh->b_face_vtx_lst[i]] = cell_flag[c_id];
    }
  }

  /* Now update coordinates */

  cs_real_34_t  *m;

  CS_MALLOC(m, tbm->n_rotors+1, cs_real_34_t);

  for (int j = 0; j < tbm->n_rotors+1; j++) {
    cs_rotation_t *r = tbm->rotation + j;

    cs_rotation_matrix(r->angle + r->omega*dt,
                       r->axis,
                       r->invariant,
                       m[j]);
  }

  for (v_id = 0; v_id < mesh->n_vertices; v_id++) {
    if (vtx_rotor_num[v_id] > 0)
      _apply_vector_transfo(m[vtx_rotor_num[v_id]],
                            &(mesh->vtx_coord[3*v_id]));
  }

  CS_FREE(m);
  CS_FREE(vtx_rotor_num);
}

/*----------------------------------------------------------------------------
 * Check that rotors and stators are originally disjoint
 *
 * parameters:
 *   mesh <-- mesh to check
 *----------------------------------------------------------------------------*/

static void
_check_geometry(cs_mesh_t  *mesh)
{
  cs_turbomachinery_t *tbm = _turbomachinery;

  cs_gnum_t n_errors = 0;

  const int  *cell_flag = tbm->cell_rotor_num;

  /* Mark from interior faces */

  for (cs_lnum_t f_id = 0; f_id < mesh->n_i_faces; f_id++) {
    cs_lnum_t c_id_0 = mesh->i_face_cells[f_id][0];
    cs_lnum_t c_id_1 = mesh->i_face_cells[f_id][1];
    if (cell_flag[c_id_1] - cell_flag[c_id_0] != 0)
      n_errors ++;
  }

  cs_parall_counter(&n_errors, 1);
  if (n_errors > 0)
    bft_error
      (__FILE__, __LINE__, 0,
       _("%s: some faces of the initial mesh belong to different\n"
         "rotor/stator sections.\n"
         "These sections must be initially disjoint to rotate freely."),
       __func__);
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
  cs_lnum_t *_cell_list = nullptr;

  cs_mesh_t *m = cs_glob_mesh;

  assert(tbm->rotor_cells_c != nullptr);

  CS_REALLOC(tbm->cell_rotor_num, m->n_cells_with_ghosts, int);

  for (cs_lnum_t i = 0; i < m->n_cells_with_ghosts; i++)
    tbm->cell_rotor_num[i] = 0;

  CS_MALLOC(_cell_list, m->n_cells_with_ghosts, cs_lnum_t);

  for (int r_id = 0; r_id < tbm->n_rotors; r_id++) {
    cs_selector_get_cell_list(tbm->rotor_cells_c[r_id],
                              &_n_cells,
                              _cell_list);
    cs_gnum_t _n_g_cells = _n_cells;
    cs_parall_counter(&_n_g_cells, 1);
    if (_n_g_cells == 0)
      bft_error
        (__FILE__, __LINE__, 0,
         _("%sd: The rotor %d with cell selection criteria\n"
           "  \"%s\"\n"
           "does not contain any cell.\n"
           "This rotor should be removed or its selection criteria modified."),
         __func__, r_id + 1, tbm->rotor_cells_c[r_id]);
    for (cs_lnum_t i = 0; i < _n_cells; i++)
      tbm->cell_rotor_num[_cell_list[i]] = r_id + 1;
  }

  CS_FREE(_cell_list);

  if (m->halo != nullptr)
    cs_halo_sync_untyped(m->halo,
                         CS_HALO_EXTENDED,
                         sizeof(int),
                         tbm->cell_rotor_num);

  if (tbm->model == CS_TURBOMACHINERY_TRANSIENT)
    _check_geometry(m);
}

/*----------------------------------------------------------------------------
 * Update mesh for unsteady rotor/stator computation when no joining is used.
 *
 * parameters:
 *   restart_mode  true for restart, false otherwise
 *   t_elapsed     elapsed computation time
 */
/*----------------------------------------------------------------------------*/

static void
_update_mesh_coupling(double  *t_elapsed)
{
  double  t_start, t_end;

  cs_turbomachinery_t *tbm = _turbomachinery;

  int t_stat_id = cs_timer_stats_id_by_name("mesh_processing");
  int t_top_id = cs_timer_stats_switch(t_stat_id);

  t_start = cs_timer_wtime();

  /* Indicates we are in the framework of turbomachinery */

  tbm->active = true;

  /* Cell and boundary face numberings can be moved from old mesh
     to new one, as the corresponding parts of the mesh should not change */

  _copy_mesh(tbm->reference_mesh, cs_glob_mesh);

  /* Update geometry, if necessary */

  _update_angle(cs_glob_time_step);

  if (tbm->n_rotors > 0)
    _update_geometry(cs_glob_mesh, 0);

  /* Recompute geometric quantities related to the mesh */

  cs_mesh_quantities_compute(cs_glob_mesh, cs_glob_mesh_quantities);

  /* Update linear algebra APIs relative to mesh */

  t_end = cs_timer_wtime();

  *t_elapsed = t_end - t_start;

  cs_timer_stats_switch(t_top_id);
}

/*----------------------------------------------------------------------------
 * Update mesh for unsteady rotor/stator computation.
 *
 * parameters:
 *   restart_mode  true for restart, false otherwise
 *   t_elapsed     elapsed computation time
 */
/*----------------------------------------------------------------------------*/

static void
_update_mesh(bool     restart_mode,
             double  *t_elapsed)
{
  double  t_start, t_end;

  cs_halo_type_t halo_type = cs_glob_mesh->halo_type;
  cs_turbomachinery_t *tbm = _turbomachinery;

  int t_stat_id = cs_timer_stats_id_by_name("mesh_processing");
  int t_top_id = cs_timer_stats_switch(t_stat_id);

  t_start = cs_timer_wtime();

  /* Indicates we are in the framework of turbomachinery */

  tbm->active = true;

  /* In case of simple coupling, use simpler update */

  if (cs_glob_n_joinings < 1) {
    _update_mesh_coupling(t_elapsed);
    return;
  }

  /* Cell and boundary face numberings can be moved from old mesh
     to new one, as the corresponding parts of the mesh should not change */

  cs_numbering_t *cell_numbering = nullptr;

  if (restart_mode == false) {
    cell_numbering = cs_glob_mesh->cell_numbering;
    cs_glob_mesh->cell_numbering = nullptr;
  }

  /* Destroy previous global mesh and related entities */

  cs_mesh_quantities_free_all(cs_glob_mesh_quantities);

  cs_mesh_reinit(cs_glob_mesh);

  /* Create new global mesh and related entities */

  cs_glob_mesh->verbosity = 0;
  cs_glob_mesh_builder = cs_mesh_builder_create();

  _update_angle(cs_glob_time_step);

  if (restart_mode == false) {

    int n_retry = cs::max(tbm->n_max_join_tries, 1);
    cs_lnum_t boundary_changed = 0;
    double eps_dt = 0.;

    const cs_time_step_t *ts = cs_glob_time_step;

    do {

      n_retry -= 1;

      _copy_mesh(tbm->reference_mesh, cs_glob_mesh);

      /* Update geometry, if necessary */

      if (tbm->n_rotors > 0)
        _update_geometry(cs_glob_mesh, eps_dt);

      /* Reset the interior faces -> cells connectivity */
      /* (in order to properly build the halo of the joined mesh) */

      cs_mesh_to_builder_perio_faces(cs_glob_mesh, cs_glob_mesh_builder);

      {
        int i;
        cs_lnum_t f_id;
        cs_lnum_2_t *i_face_cells = cs_glob_mesh->i_face_cells;
        const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
        for (f_id = 0; f_id < cs_glob_mesh->n_i_faces; f_id++) {
          for (i = 0; i < 2; i++) {
            if (i_face_cells[f_id][i] >= n_cells)
              i_face_cells[f_id][i] = -1;
          }
        }
      }

      /* Join meshes and build periodicity links */

      cs_join_all(false);

      boundary_changed = 0;
      if (tbm->n_b_faces_ref > -1) {
        if (cs_glob_mesh->n_b_faces != tbm->n_b_faces_ref)
          boundary_changed = 1;
      }
      cs_parall_counter_max(&boundary_changed, 1);

      /* Check that joining has not added or removed boundary faces.
         Postprocess new faces appearing on boundary or inside of mesh
         (which assumes that joining appends new faces at the end of the mesh)
         or try again with slightly different rotation angle
       */

      if (boundary_changed) {

        const char join_err_fmt[]
          = N_("Error in turbomachinery mesh update:\n"
               "Number of boundary faces has changed from %llu to %llu.\n"
               "There are probably unjoined faces, "
               "due to an insufficiently regular mesh;\n"
               "adjusting mesh joining parameters might help.");

        cs_gnum_t n_g_b_faces_ref = tbm->n_b_faces_ref;
        cs_parall_counter(&n_g_b_faces_ref, 1);

        if (n_retry < 1)  {
          const int writer_id = -2;
          const int writer_ids[] = {writer_id};
          const int mesh_id = cs_post_get_free_mesh_id();
          cs_lnum_t b_face_count[] = {tbm->n_b_faces_ref,
                                      cs_glob_mesh->n_b_faces};
          cs_post_init_error_writer();
          cs_post_define_surface_mesh_by_func(mesh_id,
                                              _("Added boundary faces"),
                                              nullptr,
                                              _post_error_faces_select,
                                              nullptr,
                                              b_face_count,
                                              false, /* time varying */
                                              true,  /* add groups if present */
                                              false, /* auto variables */
                                              1,
                                              writer_ids);
          cs_post_activate_writer(writer_id, true);
          cs_post_write_meshes(nullptr);
          bft_error(__FILE__, __LINE__, 0,
                    _(join_err_fmt),
                    (unsigned long long)n_g_b_faces_ref,
                    (unsigned long long)cs_glob_mesh->n_g_b_faces);
        }
        else {
          eps_dt += tbm->dt_retry * ts->dt[0];
          bft_printf(_(join_err_fmt),
                     (unsigned long long)n_g_b_faces_ref,
                     (unsigned long long)cs_glob_mesh->n_g_b_faces);
          bft_printf("\nTrying again with eps_dt = %lg\n", eps_dt);

          /* Reinitialize global mesh and related entities */

          cs_mesh_reinit(cs_glob_mesh);
          cs_glob_mesh->verbosity = 0;

          cs_glob_mesh_builder = cs_mesh_builder_create();

          cs_mesh_quantities_free_all(cs_glob_mesh_quantities);
        }
      }
    } while (boundary_changed && n_retry >= 0);
  }
  else {

    cs_mesh_to_builder_partition(tbm->reference_mesh,
                                 cs_glob_mesh_builder);

    if (cs_file_isreg("restart/mesh.csm"))
      cs_preprocessor_data_add_file("restart/mesh.csm", 0, nullptr, nullptr);
    else
      cs_preprocessor_data_add_file("restart/mesh", 0, nullptr, nullptr);

    cs_preprocessor_data_read_headers(cs_glob_mesh,
                                      cs_glob_mesh_builder,
                                      false);

    if (tbm->reference_mesh->n_g_cells != cs_glob_mesh->n_g_cells)
      bft_error
        (__FILE__, __LINE__, 0,
         _("Error in turbomachinery mesh restart:\n"
           "  number of cells expected/read: %llu/%llu\n"
           "Check your restart directory."),
         (unsigned long long)tbm->reference_mesh->n_g_cells,
         (unsigned long long)cs_glob_mesh->n_g_cells);

    cs_preprocessor_data_read_mesh(cs_glob_mesh,
                                   cs_glob_mesh_builder,
                                   false);
  }

  cs_glob_mesh->n_b_faces_all = cs_glob_mesh->n_b_faces;
  cs_glob_mesh->n_g_b_faces_all = cs_glob_mesh->n_g_b_faces;

  tbm->n_b_faces_ref = cs_glob_mesh->n_b_faces;

  /* Initialize extended connectivity, ghost cells and other remaining
     parallelism-related structures */

  cs_mesh_init_halo(cs_glob_mesh, cs_glob_mesh_builder, halo_type,
                    cs_glob_mesh->verbosity, true);
  cs_mesh_update_auxiliary(cs_glob_mesh);

  /* Destroy the temporary structure used to build the main mesh */

  cs_mesh_builder_destroy(&cs_glob_mesh_builder);

  /* Update numberings (cells saved from previous, faces
     faces updated; it is important that boundary faces
     renumbering produce the same result at each iteration) */

  if (restart_mode)
    cs_renumber_cells(cs_glob_mesh);
  else
    cs_glob_mesh->cell_numbering = cell_numbering;

  cs_renumber_i_faces(cs_glob_mesh);
  cs_renumber_b_faces(cs_glob_mesh);

  /* Build group classes */

  cs_mesh_init_group_classes(cs_glob_mesh);

  /* Print info on mesh */

  if (cs_glob_mesh->verbosity > 0)
    cs_mesh_print_info(cs_glob_mesh, _("Mesh"));

  /* Compute geometric quantities related to the mesh */

  cs_mesh_quantities_compute(cs_glob_mesh, cs_glob_mesh_quantities);
  cs_mesh_bad_cells_detect(cs_glob_mesh, cs_glob_mesh_quantities);
  cs_user_mesh_bad_cells_tag(cs_glob_mesh, cs_glob_mesh_quantities);

  cs_ext_neighborhood_reduce(cs_glob_mesh,
                             cs_glob_mesh_quantities);

  /* Initialize selectors and locations for the mesh */

  cs_mesh_init_selectors();
  cs_mesh_location_build(cs_glob_mesh, -1);
  cs_volume_zone_build_all(true);
  cs_boundary_zone_build_all(true);

  /* Check coherency if debugging */

#if 0
  cs_mesh_coherency_check();
#endif

  /* Update Fortran mesh sizes and quantities */

  cs_preprocess_mesh_update_fortran();

  /* Update mapping for accelerated devices */

#if defined(HAVE_ACCEL)
  cs_preprocess_mesh_update_device();
#endif

  /* Update rotor cells flag array in case of parallelism and/or periodicity */

  if (cs_glob_mesh->halo != nullptr) {

    const cs_mesh_t *m = cs_glob_mesh;
    CS_REALLOC(tbm->cell_rotor_num,
               m->n_cells_with_ghosts,
               int);

    cs_halo_sync_untyped(m->halo,
                         CS_HALO_EXTENDED,
                         sizeof(int),
                         tbm->cell_rotor_num);

  }

  cs_gradient_free_quantities();
  cs_cell_to_vertex_free();
  cs_mesh_adjacencies_update_mesh();

  /* Update linear algebra APIs relative to mesh */

  cs_matrix_update_mesh();

  t_end = cs_timer_wtime();

  *t_elapsed = t_end - t_start;

  cs_timer_stats_switch(t_top_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate the relative pressure associated to the given cells.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or nullptr if no
 *                               filtering is required
 * \param[in, out]  input        pointer to associated mesh structure
 *                               (to be cast as cs_mesh_t *) for interior
 *                               faces or vertices, unused otherwise
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

static void
_relative_pressure_f(int               location_id,
                     cs_lnum_t         n_elts,
                     const cs_lnum_t  *elt_ids,
                     void             *input,
                     void             *vals)
{
  CS_UNUSED(input);
  assert(location_id == CS_MESH_LOCATION_CELLS);

  const cs_turbomachinery_t *tbm = _turbomachinery;

  cs_real_t *p_rel = reinterpret_cast<cs_real_t *>(vals);

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_real_3_t *cell_cen = mq->cell_cen;
  const cs_real_t *cvar_pr = cs_field_by_name("pressure")->val;
  const cs_real_t *cpro_rho = cs_field_by_name("density")->val;

  if (elt_ids != nullptr) {
    for (cs_lnum_t idx = 0; idx <  n_elts; idx++) {
      cs_lnum_t i = elt_ids[idx];
      int r_num = tbm->cell_rotor_num[i];
      cs_real_t vr[3];
      if (r_num > 0) {
        cs_rotation_velocity(tbm->rotation + r_num, cell_cen[i], vr);
      }
      else {
        vr[0] = 0; vr[1] = 0; vr[2] = 0;
      }
      p_rel[idx] = cvar_pr[i] - cpro_rho[i] * 0.5 * cs_math_3_square_norm(vr);
    }
  }

  else {
    for (cs_lnum_t i = 0; i <  n_elts; i++) {
      int r_num = tbm->cell_rotor_num[i];
      cs_real_t vr[3];
      if (r_num > 0) {
        cs_rotation_velocity(tbm->rotation + r_num, cell_cen[i], vr);
      }
      else {
        vr[0] = 0; vr[1] = 0; vr[2] = 0;
      }
      p_rel[i] = cvar_pr[i] - cpro_rho[i] * 0.5 * cs_math_3_square_norm(vr);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate the relative velocity associated to the given cells.
 *
 * \param[in]       location_id  base associated mesh location id
 * \param[in]       n_elts       number of associated elements
 * \param[in]       elt_ids      ids of associated elements, or nullptr if no
 *                               filtering is required
 * \param[in, out]  input        pointer to associated mesh structure
 *                               (to be cast as cs_mesh_t *) for interior
 *                               faces or vertices, unused otherwise
 * \param[in, out]  vals         pointer to output values
 *                               (size: n_elts*dimension)
 */
/*----------------------------------------------------------------------------*/

static void
_relative_velocity_f(int               location_id,
                     cs_lnum_t         n_elts,
                     const cs_lnum_t  *elt_ids,
                     void             *input,
                     void             *vals)
{
  CS_UNUSED(input);
  assert(location_id == CS_MESH_LOCATION_CELLS);

  const cs_turbomachinery_t *tbm = _turbomachinery;

  cs_real_3_t *v_rel = reinterpret_cast<cs_real_3_t *>(vals);

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_real_3_t *cell_cen = mq->cell_cen;
  const cs_real_3_t *cvar_vel
    = (const cs_real_3_t *)cs_field_by_name("velocity")->val;

  if (elt_ids != nullptr) {
    for (cs_lnum_t idx = 0; idx <  n_elts; idx++) {
      cs_lnum_t i = elt_ids[idx];
      int r_num = tbm->cell_rotor_num[i];
      cs_real_t vr[3];
      if (r_num > 0) {
        cs_rotation_velocity(tbm->rotation + r_num, cell_cen[i], vr);
      }
      else {
        vr[0] = 0; vr[1] = 0; vr[2] = 0;
      }
      for (cs_lnum_t j = 0; j < 3; j++)
        v_rel[idx][j] = cvar_vel[i][j] - vr[j];
    }
  }

  else {
    for (cs_lnum_t i = 0; i <  n_elts; i++) {
      int r_num = tbm->cell_rotor_num[i];
      cs_real_t vr[3];
      if (r_num > 0) {
        cs_rotation_velocity(tbm->rotation + r_num, cell_cen[i], vr);
      }
      else {
        vr[0] = 0; vr[1] = 0; vr[2] = 0;
      }
      for (cs_lnum_t j = 0; j < 3; j++)
        v_rel[i][j] = cvar_vel[i][j] - vr[j];
    }
  }
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Assist map turbomachinery to fortran module
 *
 * parameters:
 *   iturbo <-- turbomachinery type flag
 *----------------------------------------------------------------------------*/

void
cs_f_map_turbomachinery_model(int  *iturbo)
{
  if (_turbomachinery != nullptr)
    *iturbo = _turbomachinery->model;
  else
    *iturbo = CS_TURBOMACHINERY_NONE;
}

/*----------------------------------------------------------------------------
 * Map turbomachinery arrays associated to wall BC update
 *----------------------------------------------------------------------------*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

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
      && _turbomachinery != nullptr) {
    cs_turbomachinery_finalize();
    return;
  }

  else if (_turbomachinery == nullptr)
    _turbomachinery = _turbomachinery_create();

  _turbomachinery->model = model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return rotor/stator model.
 */
/*----------------------------------------------------------------------------*/

cs_turbomachinery_model_t
cs_turbomachinery_get_model(void)
{
  if (_turbomachinery == nullptr)
   return CS_TURBOMACHINERY_NONE;
  else
    return _turbomachinery->model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of boundary couplings used for rotor/stator model.
 *
 * Joining-based definitions are not counted here.
 *
 * \return  number of boundary couplings used for rotor/stator model.
 */
/*----------------------------------------------------------------------------*/

int
cs_turbomachinery_get_n_couplings(void)
{
  int retval = 0;

  if (_turbomachinery != nullptr)
    retval = _turbomachinery->n_couplings;

  return retval;
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
  cs_turbomachinery_t *tbm = _turbomachinery;
  if (tbm == nullptr)
    return;

  const double *v = rotation_axis;
  double len = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

  int r_id = tbm->n_rotors;
  tbm->n_rotors += 1;

  CS_REALLOC(tbm->rotation, tbm->n_rotors + 1, cs_rotation_t);
  cs_rotation_t *r = tbm->rotation + r_id + 1;
  r->omega = rotation_velocity;
  r->angle = 0;
  for (int i = 0; i < 3; i++) {
    r->axis[i] = v[i]/len;
    r->invariant[i] = rotation_invariant[i];
  }

  CS_REALLOC(tbm->rotor_cells_c, tbm->n_rotors, char *);
  CS_MALLOC(tbm->rotor_cells_c[r_id], strlen(cell_criteria) + 1, char);
  strcpy(tbm->rotor_cells_c[r_id], cell_criteria);
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

  CS_REALLOC(cs_glob_join_array,  cs_glob_n_joinings + 1, cs_join_t *);

  cs_glob_join_array[cs_glob_n_joinings]
    = cs_join_create(cs_glob_n_joinings + 1,
                     sel_criteria,
                     fraction,
                     plane,
                     FVM_PERIODICITY_NULL,
                     nullptr,
                     verbosity,
                     visualization,
                     false);

  cs_glob_join_count++; /* Store number of joining (without periodic ones) */
  cs_glob_n_joinings++;

  /* Set mesh modification type */

  if (_turbomachinery != nullptr) {
    if (_turbomachinery->model == CS_TURBOMACHINERY_TRANSIENT) {
      cs_glob_mesh->time_dep = CS_MESH_TRANSIENT_CONNECT;
      _turbomachinery->reference_mesh->time_dep = CS_MESH_TRANSIENT_CONNECT;
    }
  }

  return cs_glob_n_joinings;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a boundary coupling to the list of rotor/stator couplings.
 *
 * \param[in]  sel_criteria   boundary face selection criteria
 * \param[in]  tolerance      value of the search tolerance
 * \param[in]  verbosity      level of verbosity required
 *
 * \return  number (1 to n) associated with new coupling
 */
/*----------------------------------------------------------------------------*/

int
cs_turbomachinery_coupling_add(const char  *sel_criteria,
                               float        tolerance,
                               int          verbosity)
{
  cs_sat_coupling_add_internal(_turbomachinery_coupling_tag,
                               _turbomachinery,
                               sel_criteria,
                               nullptr,
                               nullptr,
                               "all[]",
                               tolerance,
                               0, /* reverse */
                               verbosity);

  _turbomachinery->n_couplings += 1;

  if (_turbomachinery != nullptr) {
    if (   _turbomachinery->model == CS_TURBOMACHINERY_TRANSIENT
        && cs_glob_mesh->time_dep < CS_MESH_TRANSIENT_COORDS)
      cs_glob_mesh->time_dep = CS_MESH_TRANSIENT_COORDS;
  }

  return cs_sat_coupling_n_couplings();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update mesh for unsteady rotor/stator computation.
 *
 * \param[out]  t_elapsed     elapsed computation time
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_update_mesh(double  *t_elapsed)
{
  _update_mesh(false, t_elapsed);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update mesh for unsteady rotor/stator computation in case of restart.
 *
 * Reads mesh from checkpoint when available.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_restart_mesh(void)
{
  if (cs_turbomachinery_get_model() != CS_TURBOMACHINERY_TRANSIENT)
    return;

  if (cs_glob_time_step->nt_prev > 0) {
    double t_elapsed;
    if (cs_file_isreg("restart/mesh") || cs_file_isreg("restart/mesh.csm"))
      _update_mesh(true, &t_elapsed);
    else
      _update_mesh(false, &t_elapsed);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Definitions for turbomachinery computation.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_define(void)
{
  /* Define model; could be moved anywhere before time loop. */

  cs_gui_turbomachinery();
  cs_user_turbomachinery();

  if (_turbomachinery == nullptr)
    return;

  cs_turbomachinery_t *tbm = _turbomachinery;

  if (tbm->model == CS_TURBOMACHINERY_NONE)
    return;

  /* Define rotor(s) */

  cs_gui_turbomachinery_rotor();
  cs_user_turbomachinery_rotor();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initializations for turbomachinery computation.
 *
 * \note This function should be called after the mesh is built,
 *       but before cs_post_init_meshes() so that postprocessing meshes are
 *       updated correctly in the transient case.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_initialize(void)
{
  if (_turbomachinery == nullptr)
    return;

  cs_turbomachinery_t *tbm = _turbomachinery;

  if (tbm->model == CS_TURBOMACHINERY_NONE)
    return;

  /* Select rotor cells */

  _select_rotor_cells(tbm);

  /* Build the reference mesh that duplicates the global mesh before joining;
     first remove the boundary face numbering, as it will need to be
     rebuilt after the first joining */

  if (cs_glob_mesh->b_face_numbering != nullptr && cs_glob_n_joinings > 0)
    cs_numbering_destroy(&(cs_glob_mesh->b_face_numbering));

  _copy_mesh(cs_glob_mesh, tbm->reference_mesh);

  /* Reorder reference mesh by global number to avoid some issues with
     joining, especially in serial mode where global numbers are not
     expected to be present at joining stages */

  cs_renumber_i_faces_by_gnum(tbm->reference_mesh);
  cs_renumber_b_faces_by_gnum(tbm->reference_mesh);

  /* Complete the mesh with rotor-stator joining */

  if (cs_glob_n_joinings > 0) {
    cs_real_t t_elapsed;
    cs_turbomachinery_update_mesh(&t_elapsed);
  }

  /* Adapt postprocessing options if required;
     must be called before cs_post_init_meshes(). */

  if (tbm->model == CS_TURBOMACHINERY_TRANSIENT) {
    cs_post_set_changing_connectivity();

    CS_MALLOC(tbm->coftur, cs_glob_mesh->n_b_faces, cs_real_t);
    CS_MALLOC(tbm->hfltur, cs_glob_mesh->n_b_faces, cs_real_t);
  }

  /* Destroy the reference mesh, if required */

  if (tbm->model == CS_TURBOMACHINERY_FROZEN) {
    cs_mesh_destroy(tbm->reference_mesh);
    tbm->reference_mesh = nullptr;
  }

  /* Set global rotations pointer */

  cs_glob_rotation = tbm->rotation;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free turbomachinery structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_finalize(void)
{
  if (_turbomachinery != nullptr) {

    cs_turbomachinery_t *tbm = _turbomachinery;

    for (int i = tbm->n_rotors -1; i >= 0; i--)
      CS_FREE(tbm->rotor_cells_c[i]);
    CS_FREE(tbm->rotor_cells_c);

    CS_FREE(tbm->rotation);

    CS_FREE(tbm->cell_rotor_num);

    if (tbm->reference_mesh != nullptr)
      cs_mesh_destroy(tbm->reference_mesh);

    CS_FREE(tbm->coftur);
    CS_FREE(tbm->hfltur);

    /* Unset global rotations pointer for safety */
    cs_glob_rotation = nullptr;
  }

  CS_FREE(_turbomachinery);
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

      /* Ghost cell sizes may change, but the number of main cells
         is unchanged, so a simple reallocation will do */

      for (int kk = 0; kk < f->n_time_vals; kk++) {

        CS_REALLOC(f->vals[kk], _n_cells*f->dim, cs_real_t);

        if (halo != nullptr) {

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
      if (f->n_time_vals > 1)
        f->val_pre = f->vals[1];

      if (f->grad != nullptr) {

        CS_REALLOC(f->grad, _n_cells*f->dim*3, cs_real_t);

        if (halo != nullptr) {
          cs_halo_sync_var_strided(halo,
                                   CS_HALO_EXTENDED,
                                   f->grad,
                                   3*f->dim);

          if (f->dim == 1)
            cs_halo_perio_sync_var_vect(halo,
                                        CS_HALO_EXTENDED,
                                        f->grad,
                                        3);
          else if (f->dim == 3)
            cs_halo_perio_sync_var_tens(halo,
                                        CS_HALO_EXTENDED,
                                        f->grad);
          else if (f->dim == 6)
            cs_halo_perio_sync_var_sym_tens_grad(halo,
                                                 CS_HALO_EXTENDED,
                                                 f->grad);
        }
      }
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
  cs_turbomachinery_t *tbm = _turbomachinery;
  const cs_rotation_t *r = tbm->rotation + rotor_num;

  cs_rotation_matrix(theta, r->axis, r->invariant, matrix);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of rotors.
 *
 * Note that the number of associated rotations is n_rotors + 1, as the
 * first rotation id is reserved for the fixed portion of the domain.
 *
 * \return  number of rotors
 */
/*----------------------------------------------------------------------------*/

int
cs_turbomachinery_n_rotors(void)
{
  int retval = 1;

  if (_turbomachinery != nullptr)
    retval = _turbomachinery->n_rotors;

  return retval;
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
  return _turbomachinery->cell_rotor_num;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get arrays associated to wall BC update.
 *
 * \param[out] coftur  values of "cofimp" term before geometry update
 * \param[out] hfltur  local exchange coefficient before geometry update
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_get_wall_bc_coeffs(cs_real_t  **coftur,
                                     cs_real_t  **hfltur)
{
  *coftur = _turbomachinery->coftur;
  *hfltur = _turbomachinery->hfltur;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return rotation velocity
 *
 * \param[in]   rotor_num rotor number (1 to n numbering)
 */
/*----------------------------------------------------------------------------*/

double
cs_turbomachinery_get_rotation_velocity(int  rotor_num)
{
  cs_turbomachinery_t *tbm = _turbomachinery;

  return (tbm->rotation + rotor_num)->omega;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set rotation velocity
 *
 * param[in]  rotor_num  rotor number (1 to n numbering)
 * param[in]  omega      rotation velocity
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_set_rotation_velocity(int     rotor_num,
                                        double  omega)
{
  cs_turbomachinery_t *tbm = _turbomachinery;

  (tbm->rotation + rotor_num)->omega = omega;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set turbomachinery joining retry parameters.
 *
 * When a joing leads to a different number of boundary faces from the
 * previous position, the rotor positions may be perturbed by a small
 * quantity to try to obtain a better joining.
 *
 * param[in]  n_max_join_retries   maximum number of retries before considering
 *                                 the joining has failed
 * param[in]  dt_retry_multiplier  time step multiplier for new position retry
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_set_rotation_retry(int     n_max_join_retries,
                                     double  dt_retry_multiplier)
{
  cs_turbomachinery_t *tbm = _turbomachinery;

  tbm->n_max_join_tries = n_max_join_retries;
  tbm->dt_retry = dt_retry_multiplier;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build rotation matrices for a given time interval.
 *
 * The caller is responsible for freeing the array when not needed.
 *
 * \param[in]  dt  associated time delta (0 for current, unmodified time)
 *
 * \return  array of rotation matrices.
 */
/*----------------------------------------------------------------------------*/

cs_real_34_t *
cs_turbomachinery_get_rotation_matrices(double dt)
{
  const cs_turbomachinery_t *tbm = _turbomachinery;
  cs_real_34_t  *m;

  CS_MALLOC(m, tbm->n_rotors+1, cs_real_34_t);

  for (int j = 0; j < tbm->n_rotors+1; j++) {
    cs_rotation_t *r = tbm->rotation + j;

    cs_rotation_matrix(r->omega*dt,
                       r->axis,
                       r->invariant,
                       m[j]);
  }

  return m;
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
  cs_real_34_t  *m;

  cs_turbomachinery_t *tbm = _turbomachinery;
  cs_real_t time_step = dt[0];

  CS_MALLOC(m, tbm->n_rotors+1, cs_real_34_t);

  for (int j = 0; j < tbm->n_rotors+1; j++) {
    cs_rotation_t *r = tbm->rotation + j;
    cs_rotation_matrix(r->omega * time_step,
                       r->axis,
                       r->invariant,
                       m[j]);
  }

  int n_fields = cs_field_n_fields();

  for (int f_id = 0; f_id < n_fields; f_id++) {

    cs_field_t *f = cs_field_by_id(f_id);

    if (! (f->dim > 1 && (f->type & CS_FIELD_VARIABLE)))
      continue;

    const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(f->location_id);
    cs_lnum_t _n_elts = n_elts[2];

    if (f->dim == 3) {
      for (i = 0; i < _n_elts; i++)
        _apply_vector_rotation(m[tbm->cell_rotor_num[i]], f->val + i*3);
    }

    else if (f->dim == 6) {
      for (i = 0; i < _n_elts; i++)
        _apply_sym_tensor_rotation(m[tbm->cell_rotor_num[i]], f->val + i*6);
    }

    assert(f->dim == 3 || f->dim == 6);

  } /* End of loop on fields */

  /* Specific handling of Reynolds stresses */

  cs_field_t  *frij = cs_field_by_name("rij");
  if (frij != nullptr) {
    const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(frij->location_id);
    cs_lnum_t _n_elts = n_elts[2];
    cs_real_6_t *cvar_r = (cs_real_6_t *)(frij->val);
    for (i = 0; i < _n_elts; i++) {
      _apply_sym_tensor_rotation(m[tbm->cell_rotor_num[i]], cvar_r[i]);
    }

  }

  CS_FREE(m);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute velocity relative to fixed coordinates at a given point
 *
 * \deprecated
 * Use \ref cs_rotation_velocity for more consistent naming of this reference
 * frame velocity.
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
  cs_turbomachinery_t *tbm = _turbomachinery;
  const cs_rotation_t *r = tbm->rotation + rotor_num;

  _relative_velocity(r->omega,
                     r->axis,
                     r->invariant,
                     coords,
                     velocity);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read turbomachinery metadata from restart file.
 *
 * The mesh is handled separately.
 *
 * \param[in, out]  r  associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_restart_read(cs_restart_t  *r)
{
  cs_turbomachinery_t *tbm = _turbomachinery;

  if (tbm == nullptr)
    return;

  cs_real_t *t_angle;
  CS_MALLOC(t_angle, tbm->n_rotors+2, cs_real_t);

  t_angle[0] = tbm->t_cur;
  for (int i = 0; i < tbm->n_rotors+1; i++) {
    cs_rotation_t *rot = tbm->rotation + i;
    t_angle[i+1] = rot->angle;
  }

  int retcode = cs_restart_read_section(r,
                                        "turbomachinery:rotor_time_and_angle",
                                        CS_MESH_LOCATION_NONE,
                                        tbm->n_rotors+2,
                                        CS_TYPE_cs_real_t,
                                        t_angle);

  if (retcode == CS_RESTART_SUCCESS) {
    tbm->nt_cur = cs_glob_time_step->nt_prev;
    tbm->t_cur = t_angle[0];
    for (int i = 0; i < tbm->n_rotors+1; i++) {
      cs_rotation_t *rot = tbm->rotation + i;
      rot->angle = t_angle[i+1];
    }
  }

  CS_FREE(t_angle);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write turbomachinery metadata to checkpoint file.
 *
 * The mesh is handled separately.
 *
 * \param[in, out]  r  associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_restart_write(cs_restart_t  *r)
{
  const cs_turbomachinery_t *tbm = _turbomachinery;

  if (tbm == nullptr)
    return;

  if (tbm->model == CS_TURBOMACHINERY_NONE)
    return;

  cs_real_t *t_angle;
  CS_MALLOC(t_angle, tbm->n_rotors+2, cs_real_t);

  t_angle[0] = tbm->t_cur;
  for (int i = 0; i < tbm->n_rotors+1; i++) {
    cs_rotation_t *rot = tbm->rotation + i;
    t_angle[i+1] = rot->angle;
  }

  cs_restart_write_section(r,
                           "turbomachinery:rotor_time_and_angle",
                           CS_MESH_LOCATION_NONE,
                           tbm->n_rotors+2,
                           CS_TYPE_cs_real_t,
                           t_angle);

  CS_FREE(t_angle);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create or access function objects specific to
 *        turbomachinery models (relative_pressure, relative_velocity).
 */
/*----------------------------------------------------------------------------*/

void
cs_turbomachinery_define_functions(void)
{
  if (_turbomachinery == nullptr)
    return;

  /* Relative pressure */

  {
    cs_function_t *f
      = cs_function_define_by_func("relative_pressure",
                                   CS_MESH_LOCATION_CELLS,
                                   1,
                                   true,
                                   CS_REAL_TYPE,
                                   _relative_pressure_f,
                                   nullptr);

    const char label[] = "Rel Pressure";
    CS_MALLOC(f->label, strlen(label) + 1, char);
    strcpy(f->label, label);

    f->type = CS_FUNCTION_INTENSIVE;
    f->post_vis = CS_POST_ON_LOCATION;
  }

  /* Relative velocity */

  {
    cs_function_t *f
      = cs_function_define_by_func("relative_velocity",
                                   CS_MESH_LOCATION_CELLS,
                                   3,
                                   true,
                                   CS_REAL_TYPE,
                                   _relative_velocity_f,
                                   nullptr);

    const char label[] = "Rel Velocity";
    CS_MALLOC(f->label, strlen(label) + 1, char);
    strcpy(f->label, label);

    f->type = CS_FUNCTION_INTENSIVE;
    f->post_vis = CS_POST_ON_LOCATION;
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
