/*============================================================================
 * Handling of new particles.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

/*============================================================================
 * Functions dealing with particle tracking
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <limits.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_printf.h"
#include "bft_error.h"
#include "bft_mem.h"

#include "cs_base.h"
#include "cs_coal.h"
#include "cs_halo.h"
#include "cs_ht_convert.h"
#include "cs_log.h"
#include "cs_interface.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_adjacencies.h"
#include "cs_mesh_quantities.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_physical_model.h"
#include "cs_physical_constants.h"
#include "cs_random.h"
#include "cs_search.h"
#include "cs_thermal_model.h"
#include "cs_timer_stats.h"

#include "cs_field.h"
#include "cs_field_pointer.h"

#include "cs_lagr_clogging.h"
#include "cs_lagr_deposition_model.h"
#include "cs_lagr_roughness.h"
#include "cs_lagr_dlvo.h"
#include "cs_lagr_stat.h"
#include "cs_lagr.h"
#include "cs_lagr_tracking.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_new.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Determine accumulated surfaces of an implicit decomposition
 *        of a polygonal face.
 *
 * If the face is degenerate or has one or more concave angles,
 * the perimeter is used instead, and the returned value is negative.
 *
 * \param[in]   n_vertices     number of face vertices
 * \param[in]   vertex_ids     ids of face vertices
 * \param[in]   vertex_coords  vertex coordinates
 * \param[in]   face_center    coordinates of face center
 * \param[out]  acc_surf_r     accumulated surface ratio associated to
 *                             each edge (or negative edge lengths in
 *                             degenerate cases)
 *
 * \return
 *   total face surface (or - face perimeter if locally inverted)
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_face_sub_surfaces(cs_lnum_t        n_vertices,
                   const cs_lnum_t  vertex_ids[],
                   const cs_real_t  vertex_coords[][3],
                   const cs_real_t  face_center[3],
                   cs_real_t        acc_surf_r[])
{
  const double epsilon = 1.e-24;

  cs_real_t s = 0;

  cs_real_t v0[3], v1[3], vn[2][3] = {{0, 0, 0}, {0, 0, 0}};

  for (cs_lnum_t tri_id = 0; tri_id < n_vertices; tri_id++) {

    cs_lnum_t v_id_0 = vertex_ids[tri_id];
    cs_lnum_t v_id_1 = vertex_ids[(tri_id+1) % n_vertices];

    for (cs_lnum_t i = 0; i < 3; i++) {
      v0[i] = vertex_coords[v_id_0][i] - face_center[i];
      v1[i] = vertex_coords[v_id_1][i] - face_center[i];
    }

    cs_math_3_cross_product(v0, v1, vn[tri_id%2]);

    if (cs_math_3_dot_product(vn[0], vn[1]) < 0) {
      s = -1;
      break;
    }

    s += cs_math_3_norm(vn[tri_id%2]);

    acc_surf_r[tri_id] = s;

  }

  /* fallback in case of exterior face center */

  if (s <= 0) {

    for (cs_lnum_t tri_id = 0; tri_id < n_vertices; tri_id++) {

      cs_lnum_t v_id_0 = vertex_ids[tri_id];
      cs_lnum_t v_id_1 = vertex_ids[(tri_id+1) % n_vertices];

      for (cs_lnum_t i = 0; i < 3; i++)
        v0[i] = vertex_coords[v_id_1][i] - vertex_coords[v_id_0][i];

      s -= cs_math_3_norm(v0);

      acc_surf_r[tri_id] = s;

    }

  }

  /* Normalize */

  cs_real_t sd = CS_ABS(s);

  if (sd >= epsilon) {
    for (cs_lnum_t tri_id = 0; tri_id < n_vertices; tri_id++)
      acc_surf_r[tri_id] /= sd;
  }
  else {
    for (cs_lnum_t tri_id = 0; tri_id < n_vertices; tri_id++)
      acc_surf_r[tri_id] = 1;
  }

  /* avoid rounding error issues */

  if (s < 0)
    acc_surf_r[n_vertices-1] = -1;
  else
    acc_surf_r[n_vertices-1] = 1;

  return s;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Position a point randomly on a given face
 *
 * If the face has one or more concave angles, the point will be assigned
 * to a randomly determined edge.
 *
 * \param[in]   n_vertices     number of face vertices
 * \param[in]   vertex_ids     ids of face vertices
 * \param[in]   vertex_coords  vertex coordinates
 * \param[in]   face_center    coordinates of face center
 * \param[in]   acc_surf_r     accumulated surface ratio associated to
 *                             each edge (or negative edge lengths in
 *                             degenerate cases)
 * \param[out]  coords         coordinates of point in face
 */
/*----------------------------------------------------------------------------*/

static void
_random_point_in_face(cs_lnum_t        n_vertices,
                      const cs_lnum_t  vertex_ids[],
                      const cs_real_t  vertex_coords[][3],
                      const cs_real_t  face_center[3],
                      const cs_real_t  acc_surf_r[],
                      cs_real_t        coords[3])
{
  cs_lnum_t tri_id = 0;
  cs_real_t r[3];

  /* determine triangle to choose */

  cs_random_uniform(3, r);

  if (r[2] > 1) /* account for possible ? rounding errors */
    r[2] = 1;

  while (tri_id < n_vertices && r[2] > acc_surf_r[tri_id])
    tri_id++;

  /* randomize inside triangle */

  if (tri_id < n_vertices) {

    /* account for possible ? rounding errors */
    for (int i = 0; i < 2; i++) {
      if (r[i] < 0)
        r[i] = 0;
    }

    /* Distribution based on \cite Osada:2002 */

    r[0] = sqrt(r[0]);
    cs_real_t a[3] = {1. - r[0], r[0]*(1.-r[1]), r[0]*r[1]};

    cs_lnum_t v_id_0 = vertex_ids[tri_id];
    cs_lnum_t v_id_1 = vertex_ids[(tri_id+1) % n_vertices];

    for (cs_lnum_t j = 0; j < 3; j++)
      coords[j] =   a[0] * face_center[j]
                  + a[1] * vertex_coords[v_id_0][j]
                  + a[2] * vertex_coords[v_id_1][j];

  }

  /* fallback on edges */

  else {

    tri_id = 0;
    while (tri_id < n_vertices-1 && r[2] > -acc_surf_r[tri_id])
      tri_id++;

    cs_real_t a[2] = {r[0], 1. - r[0]};

    cs_lnum_t v_id_0 = vertex_ids[tri_id];
    cs_lnum_t v_id_1 = vertex_ids[(tri_id+1) % n_vertices];

    for (cs_lnum_t j = 0; j < 3; j++)
      coords[j] =   a[0] * vertex_coords[v_id_0][j]
                  + a[1] * vertex_coords[v_id_1][j];

  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Inject a series of particles at random positions on given faces.
 *
 * The fluid velocity and other variables and attributes are computed here.
 *
 * \param[in,out]  particles          pointer to particle set
 * \param[in]      n_faces            number of faces in zone
 * \param[in]      face_ids           ids of faces in zone
 * \param[in]      face_particle_idx  starting index of added particles
 *                                    for each face in zone
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_new(cs_lagr_particle_set_t  *particles,
            cs_lnum_t                n_faces,
            const cs_lnum_t          face_ids[],
            const cs_lnum_t          face_particle_idx[])
{
  const double d_eps = 1e-3;

  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_mesh_quantities_t *fvq  = cs_glob_mesh_quantities;

  cs_real_t  *acc_surf_r = NULL;
  cs_lnum_t   n_vertices_max = 0;

  /* Loop on faces */

  for (cs_lnum_t li = 0; li < n_faces; li++) {

    cs_lnum_t n_f_p = face_particle_idx[li+1] - face_particle_idx[li];

    if (n_f_p < 1)
      continue;

    cs_lnum_t p_s_id = particles->n_particles + face_particle_idx[li];

    const cs_lnum_t face_id = (face_ids != NULL) ? face_ids[li] : li;

    cs_lnum_t n_vertices =   mesh->b_face_vtx_idx[face_id+1]
                           - mesh->b_face_vtx_idx[face_id];

    const cs_lnum_t *vertex_ids =   mesh->b_face_vtx_lst
                                  + mesh->b_face_vtx_idx[face_id];

    if (n_vertices > n_vertices_max) {
      n_vertices_max = n_vertices*2;
      BFT_REALLOC(acc_surf_r, n_vertices_max, cs_real_t);
    }

    _face_sub_surfaces(n_vertices,
                       vertex_ids,
                       (const cs_real_3_t *)mesh->vtx_coord,
                       fvq->b_face_cog + 3*face_id,
                       acc_surf_r);

    /* distribute new particles */

    cs_lnum_t c_id = mesh->b_face_cells[face_id];
    const cs_real_t *c_cen = fvq->cell_cen + c_id*3;

    for (cs_lnum_t i = 0; i < n_f_p; i++) {

      cs_lnum_t p_id = p_s_id + i;

      cs_lagr_particles_set_lnum(particles, p_id, CS_LAGR_CELL_ID, c_id);

      cs_real_t *part_coord
        = cs_lagr_particles_attr_get_ptr<cs_real_t>(particles, p_id, CS_LAGR_COORDS);

      _random_point_in_face(n_vertices,
                            vertex_ids,
                            (const cs_real_3_t *)mesh->vtx_coord,
                            fvq->b_face_cog + 3*face_id,
                            acc_surf_r,
                            part_coord);

      /* For safety, move particle slightly inside cell */

      for (cs_lnum_t j = 0; j < 3; j++)
        part_coord[j] += (c_cen[j] - part_coord[j])*d_eps;

    }

  }

  BFT_FREE(acc_surf_r);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Inject a series of particles at random positions on given cells.
 *
 * \warning Currently works only for tri and quadrangular faces.
 *
 * The fluid velocity and other variables and attributes are computed here.
 *
 * \param[in,out]  particles          pointer to particle set
 * \param[in]      n_cells            number of cells in zone
 * \param[in]      cell_ids           ids of cells in zone
 * \param[in]      cell_particle_idx  starting index of added particles
 *                                    for each cell in zone
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_new_v(cs_lagr_particle_set_t  *particles,
              cs_lnum_t                n_cells,
              const cs_lnum_t          cell_ids[],
              const cs_lnum_t          cell_particle_idx[])
{
  const double w_eps = 1e-24;
  const double d_eps = 1e-3;

  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq  = cs_glob_mesh_quantities;

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  if (ma->cell_i_faces == NULL)
    cs_mesh_adjacencies_update_cell_i_faces();

  cs_lnum_t  *cell_subface_index = NULL;
  cs_real_t  *acc_vol_r = NULL;
  cs_real_t  *acc_surf_r = NULL;
  cs_lnum_t  n_divisions_max = 0, n_faces_max = 0;

  /* Loop on cells */

  for (cs_lnum_t li = 0; li < n_cells; li++) {

    cs_lnum_t n_c_p = cell_particle_idx[li+1] - cell_particle_idx[li];

    if (n_c_p < 1) /* ignore cells with no injected particles */
      continue;

    cs_lnum_t p_s_id = particles->n_particles +  cell_particle_idx[li];

    const cs_lnum_t cell_id = (cell_ids != NULL) ? cell_ids[li] : li;

    const cs_real_t *cell_cen = fvq->cell_cen + cell_id*3;


    const cs_lnum_t n_cell_i_faces =   ma->cell_cells_idx[cell_id+1]
                                     - ma->cell_cells_idx[cell_id];
    const cs_lnum_t n_cell_b_faces =   ma->cell_b_faces_idx[cell_id+1]
                                     - ma->cell_b_faces_idx[cell_id];

    cs_lnum_t n_cell_faces = n_cell_i_faces + n_cell_b_faces;

    if (ma->cell_hb_faces_idx != NULL) {
      n_cell_faces +=   ma->cell_hb_faces_idx[cell_id+1]
                      - ma->cell_hb_faces_idx[cell_id];
    }

    if (n_cell_faces > n_faces_max) {
      n_faces_max = n_cell_faces*2;
      BFT_REALLOC(cell_subface_index, n_faces_max+1, cs_lnum_t);
      BFT_REALLOC(acc_vol_r, n_faces_max, cs_real_t);
    }

    cell_subface_index[0] = 0;

    /* Loop on cell faces to determine volumes */

    bool fallback = false;
    cs_real_t t_vol = 0;

    for (cs_lnum_t i = 0; i < n_cell_faces; i++) {

      cs_lnum_t face_id, n_vertices;
      const cs_lnum_t *vertex_ids;
      const cs_real_t *face_cog, *face_normal;

      /* Outward normal: always well oriented for external faces,
         depend on the connectivity for internal faces */

      cs_real_t v_mult = 1;

      if (i < n_cell_i_faces) { /* Interior face */

        face_id = ma->cell_i_faces[ma->cell_cells_idx[cell_id] + i];

        if (cell_id == mesh->i_face_cells[face_id][1])
          v_mult = -1;
        cs_lnum_t vtx_s = mesh->i_face_vtx_idx[face_id];
        n_vertices = mesh->i_face_vtx_idx[face_id+1] - vtx_s;
        vertex_ids = mesh->i_face_vtx_lst + vtx_s;
        face_cog = fvq->i_face_cog + (3*face_id);
        face_normal = fvq->i_face_normal + (3*face_id);

      }
      else { /* Boundary faces */

        cs_lnum_t j = i - n_cell_i_faces;
        if (j < n_cell_b_faces)
          face_id = ma->cell_b_faces[ma->cell_b_faces_idx[cell_id] + j];

        else {
          assert(ma->cell_hb_faces_idx != NULL);
          j -= n_cell_b_faces;
          face_id = ma->cell_hb_faces[ma->cell_hb_faces_idx[cell_id] + j];
        }

        cs_lnum_t vtx_s = mesh->b_face_vtx_idx[face_id];
        n_vertices = mesh->b_face_vtx_idx[face_id+1] - vtx_s;
        vertex_ids = mesh->b_face_vtx_lst + vtx_s;
        face_cog = fvq->b_face_cog + (3*face_id);
        face_normal = fvq->b_face_normal + (3*face_id);

      }

      cell_subface_index[i+1] = cell_subface_index[i] + n_vertices;

      if (cell_subface_index[i+1] > n_divisions_max) {
        n_divisions_max = cell_subface_index[i+1]*2;
        BFT_REALLOC(acc_surf_r, n_divisions_max, cs_real_t);
      }

      cs_real_t f_surf
        = _face_sub_surfaces(n_vertices,
                             vertex_ids,
                             (const cs_real_3_t *)mesh->vtx_coord,
                             face_cog,
                             acc_surf_r + cell_subface_index[i]);

      cs_real_t fh = 0;
      if (f_surf > 0) {
        /* face normal should have length f_surf, so no need to divide here */
        for (cs_lnum_t j = 0; j < 3; j++)
          fh += (face_cog[j] - cell_cen[j]) * face_normal[j];
      }
      fh *= v_mult;

      t_vol += CS_ABS(fh);
      acc_vol_r[i] = t_vol;

      if (fh <= 0 || f_surf <= 0)
        fallback = true;

    }

    if (t_vol >= w_eps) {
      for (cs_lnum_t i = 0; i < n_cell_faces; i++)
        acc_vol_r[i] /= t_vol;
    }
    else {
      for (cs_lnum_t i = 0; i < n_cell_faces; i++)
        acc_vol_r[i] = 1;
    }
    acc_vol_r[n_cell_faces - 1] = 1;

    /* If needed, apply fallback to all faces, as in non-convex cases,
       some cones may be partially masked by inverted cones;
       weight is not based strictly on edge length in this case,
       but bias cannot be avoid in this mode anyways, so do not bother
       with extra steps. */

    if (fallback) {
      for (cs_lnum_t i = 0; i < cell_subface_index[n_cell_faces]; i++) {
        if (acc_surf_r[i] > 0)
          acc_surf_r[i] *= -1;
      }
    }

    /* distribute new particles */

    for (cs_lnum_t c_i = 0; c_i < n_c_p; c_i++) {

      cs_lnum_t p_id = p_s_id + c_i;

      cs_lagr_particles_set_lnum(particles, p_id, CS_LAGR_CELL_ID, cell_id);

      cs_real_t *part_coord
        = cs_lagr_particles_attr_get_ptr<cs_real_t>(particles, p_id, CS_LAGR_COORDS);

      /* search for matching center-to-face cone */

      cs_real_t r[2];
      cs_random_uniform(2, r);

      cs_lnum_t i = 0;
      while (i < n_cell_faces && r[0] > acc_vol_r[i])
        i++;

      cs_lnum_t face_id, n_vertices;
      const cs_lnum_t *vertex_ids;
      const cs_real_t *face_cog;

      if (i < n_cell_i_faces) { /* Interior face */

        face_id = ma->cell_i_faces[ma->cell_cells_idx[cell_id] + i];

        cs_lnum_t vtx_s = mesh->i_face_vtx_idx[face_id];
        n_vertices = mesh->i_face_vtx_idx[face_id+1] - vtx_s;
        vertex_ids = mesh->i_face_vtx_lst + vtx_s;
        face_cog = fvq->i_face_cog + (3*face_id);

      }
      else { /* Boundary faces */

        cs_lnum_t j = i - n_cell_i_faces;
        if (j < n_cell_b_faces)
          face_id = ma->cell_b_faces[ma->cell_b_faces_idx[cell_id] + j];

        else {
          assert(ma->cell_hb_faces_idx != NULL);
          j -= n_cell_b_faces;
          face_id = ma->cell_hb_faces[ma->cell_hb_faces_idx[cell_id] + j];
        }

        cs_lnum_t vtx_s = mesh->b_face_vtx_idx[face_id];
        n_vertices = mesh->b_face_vtx_idx[face_id+1] - vtx_s;
        vertex_ids = mesh->b_face_vtx_lst + vtx_s;
        face_cog = fvq->b_face_cog + (3*face_id);

      }

      _random_point_in_face(n_vertices,
                            vertex_ids,
                            (const cs_real_3_t *)mesh->vtx_coord,
                            face_cog,
                            acc_surf_r + cell_subface_index[i],
                            part_coord);

      /* In regular case, place point on segment joining cell center and
         point in cell; volume of truncated cone proportional to
         cube of distance along segment, so distribution compensates for this */

      if (fallback == false) {

        cs_real_t t = pow(r[1], 1./3.) * (1.0 - d_eps);
        for (cs_lnum_t j = 0; j < 3; j++)
          part_coord[j] += (cell_cen[j] - part_coord[j]) * (1. - t);
      }

      /* Move particle slightly towards cell center cell
         (assuming cell is star-shaped) */

      else {
        if (fvq->cell_vol[cell_id] > 0) {
          for (cs_lnum_t j = 0; j < 3; j++)
            part_coord[j] += (cell_cen[j] - part_coord[j])*d_eps;
        }
      }

    } /* end of loop on new particles */

  } /* end of loop on cells */

  BFT_FREE(acc_surf_r);
  BFT_FREE(acc_vol_r);
  BFT_FREE(cell_subface_index);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialization for new particles.
 *
 * The fluid velocity seen is computed here.
 *
 * \param[in]  particle_range  start and past-the-end ids of new particles
 *                             for this zone and class
 * \param[in]  time_id         associated time id (0: current, 1: previous)
 * \param[in]  visc_length     viscous layer thickness
 *                             (size: number of mesh boundary faces)
 * \param[in]      zis         injection data for this zone and set
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_new_particle_init(const cs_lnum_t                 particle_range[2],
                          int                             time_id,
                          const cs_real_t                 visc_length[],
                          const cs_lagr_injection_set_t  *zis)
{
  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;

  const cs_mesh_adjacencies_t  *ma = cs_glob_mesh_adjacencies;

  cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  const cs_lagr_zone_data_t  *bcs = cs_glob_lagr_boundary_conditions;

  cs_lagr_extra_module_t  *extra = cs_get_lagr_extra_module();
  const cs_coal_model_t  *coal_model = cs_glob_coal_model;

  /* Non-Lagrangian fields */

  const cs_real_t  *xashch = coal_model->xashch;
  const cs_real_t  *cp2ch  = coal_model->cp2ch;
  const cs_real_t  *xwatch = coal_model->xwatch;
  const cs_real_t  *rho0ch = coal_model->rho0ch;

  const cs_real_t *vela = extra->vel->vals[time_id];//FIXME
  const cs_real_t *cval_h = NULL, *cval_t = NULL, *cval_t_l = NULL;
  cs_real_t *_cval_t = NULL;

  cs_real_t tscl_shift = 0;

  const cs_real_3_t  *vel = NULL;
  const cs_real_t    *cvar_k = NULL;
  const cs_real_6_t  *cvar_rij = NULL;

  vel = (const cs_real_3_t *)extra->vel->vals[time_id];

  /* Initialize pointers (used to simplify future tests) */

  if (   (   cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_HEAT
          && cs_glob_lagr_specific_physics->itpvar == 1)
      || cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_COAL
      || cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_CTWR) {

    const cs_field_t *f = cs_field_by_name_try("temperature");
    if (f != NULL)
      cval_t = f->val;
    else if (   cs_glob_thermal_model->thermal_variable
             == CS_THERMAL_MODEL_ENTHALPY)
      cval_h = cs_field_by_name("enthalpy")->val;

    if (cs_glob_thermal_model->temperature_scale == CS_TEMPERATURE_SCALE_KELVIN)
      tscl_shift = - cs_physical_constants_celsius_to_kelvin;
  }

  if (cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_CTWR) {
    cval_t_l = cs_field_by_name("temp_l_r")->val;
  }

  const cs_real_t pis6 = cs_math_pi / 6.0;
  const int shape = cs_glob_lagr_model->shape;

  /* Prepare  enthalpy to temperature conversion if needed */

  if (   cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_HEAT
      && cs_glob_lagr_specific_physics->itpvar == 1
      && cval_t == NULL
      && cval_h != NULL) {

    BFT_MALLOC(_cval_t, cs_glob_mesh->n_cells_with_ghosts, cs_real_t);
    cs_ht_convert_h_to_t_cells(cval_h, _cval_t);
    cval_t = _cval_t;

  }


  /* Initialization */

  cs_real_t d2s3 = 2.0 / 3.0;

  /* Simulate instantaneous turbulent fluid velocities
     as seen by solid particles along their trajectory.
     -------------------------------------------------- */

  if (cs_glob_lagr_model->idistu == 1) {

    if (extra->cvar_rij != NULL)
      cvar_rij = (const cs_real_6_t *) extra->cvar_rij->vals[time_id];

    else if (extra->cvar_k != NULL) {
      cvar_k = (const cs_real_t *)extra->cvar_k->vals[time_id];
      if (extra->cvar_k != NULL)
        cvar_k = (const cs_real_t *)extra->cvar_k->val;
    }

    else {
      bft_error
        (__FILE__, __LINE__, 0,
         _("The Lagrangian module is incompatible with the selected\n"
           " turbulence model.\n\n"
           "Turbulent dispersion is used with:\n"
           "  cs_glob_lagr_model->idistu = %d\n"
           "And the turbulence model is iturb = %d\n\n"
           "The only turbulence models compatible with the Lagrangian model's\n"
           "turbulent dispersion are k-epsilon, Rij-epsilon, v2f, and k-omega."),
         cs_glob_lagr_model->idistu,
         extra->iturb);
    }

  }

  /* Random draws and computation of particle characteristic times */

  cs_lnum_t  n = particle_range[1] - particle_range[0];
  cs_real_3_t  *vagaus;

  BFT_MALLOC(vagaus, n, cs_real_3_t);

  if (cs_glob_lagr_model->idistu == 1 && n > 0) {
    cs_random_normal(n*3, (cs_real_t *)vagaus);
  }

  else {

    for (cs_lnum_t i = 0; i < n; i++) {
      vagaus[i][0] = 0.0;
      vagaus[i][1] = 0.0;
      vagaus[i][2] = 0.0;
    }

  }

  cs_real_33_t *eig_vec;
  cs_real_3_t *eig_val;

  BFT_MALLOC(eig_vec, n_cells, cs_real_33_t);
  BFT_MALLOC(eig_val, n_cells, cs_real_3_t);

  /* First stage: compute cell values
   * Initialisation from the mean Eulerian fluid
     ------------------------------------------- */

  if (cs_glob_lagr_model->idistu == 1) {

    cs_real_33_t *sym_rij;
    BFT_MALLOC(sym_rij, n_cells, cs_real_33_t);
    cs_real_t tol_err = 1.0e-12;

    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      if (cvar_rij != NULL) {
        sym_rij[cell_id][0][0] = cvar_rij[cell_id][0];
        sym_rij[cell_id][1][1] = cvar_rij[cell_id][1];
        sym_rij[cell_id][2][2] = cvar_rij[cell_id][2];
        sym_rij[cell_id][0][1] = cvar_rij[cell_id][3];
        sym_rij[cell_id][1][0] = cvar_rij[cell_id][3];
        sym_rij[cell_id][1][2] = cvar_rij[cell_id][4];
        sym_rij[cell_id][2][1] = cvar_rij[cell_id][4];
        sym_rij[cell_id][0][2] = cvar_rij[cell_id][5];
        sym_rij[cell_id][2][0] = cvar_rij[cell_id][5];
      }
      /* TODO do it better for EVM models */
      else {
        cs_real_t w = 0.;

        if (cvar_k != NULL)
          w = d2s3 * cvar_k[cell_id];

        sym_rij[cell_id][0][0] = w;
        sym_rij[cell_id][1][1] = w;
        sym_rij[cell_id][2][2] = w;
        sym_rij[cell_id][0][1] = 0.;
        sym_rij[cell_id][1][0] = 0.;
        sym_rij[cell_id][1][2] = 0.;
        sym_rij[cell_id][2][1] = 0.;
        sym_rij[cell_id][0][2] = 0.;
        sym_rij[cell_id][2][0] = 0.;
      }

      eig_vec[cell_id][0][0] = 1;
      eig_vec[cell_id][0][1] = 0;
      eig_vec[cell_id][0][2] = 0;
      eig_vec[cell_id][1][0] = 0;
      eig_vec[cell_id][1][1] = 1;
      eig_vec[cell_id][1][2] = 0;
      eig_vec[cell_id][2][0] = 0;
      eig_vec[cell_id][2][1] = 0;
      eig_vec[cell_id][2][2] = 1;

      cs_math_33_eig_val_vec(sym_rij[cell_id], tol_err, eig_val[cell_id], eig_vec[cell_id]);
    }

    BFT_FREE(sym_rij);
  } else {
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

      eig_vec[cell_id][0][0] = 1.;
      eig_vec[cell_id][1][1] = 1.;
      eig_vec[cell_id][2][2] = 1.;
      eig_vec[cell_id][0][1] = 0.;
      eig_vec[cell_id][0][2] = 0.;
      eig_vec[cell_id][1][0] = 0.;
      eig_vec[cell_id][1][2] = 0.;
      eig_vec[cell_id][2][0] = 0.;
      eig_vec[cell_id][2][1] = 0.;
      eig_val[cell_id][0] = 0.;
      eig_val[cell_id][1] = 0.;
      eig_val[cell_id][2] = 0.;
    }
  }

  /* Second stage: initialize particle attributes
     ------------------------------------------- */

  for (cs_lnum_t p_id = particle_range[0]; p_id < particle_range[1]; p_id++) {

    unsigned char *particle = p_set->p_buffer + p_am->extents * p_id;

    cs_lnum_t c_id  = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_CELL_ID);
    cs_lnum_t l_id = p_id - particle_range[0];

    /* Particle velocity components */

    cs_real_t *part_vel =
      cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                               CS_LAGR_VELOCITY);

    if (zis->velocity_profile == CS_LAGR_IN_IMPOSED_FLUID_VALUE) {
      for (cs_lnum_t i = 0; i < 3; i++)
        part_vel[i] = vel[c_id][i];
    }

    /* velocity as seen from fluid */

    cs_real_t  *vel_seen
      = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                 CS_LAGR_VELOCITY_SEEN);

    for (cs_lnum_t i = 0; i < 3; i++) {
      vel_seen[i] = vel[c_id][i]
                  + vagaus[l_id][0] * sqrt(eig_val[c_id][0]) * eig_vec[c_id][0][i]
                  + vagaus[l_id][1] * sqrt(eig_val[c_id][1]) * eig_vec[c_id][1][i]
                  + vagaus[l_id][2] * sqrt(eig_val[c_id][2]) * eig_vec[c_id][2][i];
    }

    /* Diameter (always set base) */

    cs_lagr_particle_set_real(particle, p_am, CS_LAGR_DIAMETER,
                              zis->diameter);

    if (zis->diameter_variance > 0.0) {

      /* Randomize diameter, ensuring we obtain a
         positive diameter in the 99,7% range */

      cs_real_t d3   = 3.0 * zis->diameter_variance;

      int i_r = 0; /* avoid infinite loop in case of very improbable
                      random series... */

      for (i_r = 0; i_r < 20; i_r++) {
        double    random;
        cs_random_normal(1, &random);

        cs_real_t diam =   zis->diameter
                         + random * zis->diameter_variance;

        if (diam > 0 && (   diam >= zis->diameter - d3
                         && diam <= zis->diameter + d3)) {
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_DIAMETER, diam);
          break;
        }
      }

    }

    /* Shape for spheroids without inertia */
    if (   shape == CS_LAGR_SHAPE_SPHEROID_STOC_MODEL
        || shape == CS_LAGR_SHAPE_SPHEROID_JEFFERY_MODEL) {

      /* Spherical radii a b c */
      auto *radii = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                             CS_LAGR_RADII);

      for (cs_lnum_t i = 0; i < 3; i++) {
        radii[i] = zis->radii[i];
      }

      /* Shape parameters */
      auto *shape_param =
        cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                 CS_LAGR_SHAPE_PARAM);

      /* Compute shape parameters from radii */
      /* FIXME valid for all spheroids only (a = b, c > a,b ) */
      cs_real_t lamb = radii[2] / radii[1];  /* FIXME do not divide by 0... */
      cs_real_t lamb_m1 = (radii[2] - radii[1]) / radii[1];
      cs_real_t _a2 = radii[0] * radii[0];
      /* TODO MF shape_param check development in series */
      cs_real_t aux1 = lamb * lamb;
      cs_real_t aux2 = aux1 -1;
      if (lamb_m1 > 1e-10) {
        cs_real_t aux3 = sqrt(aux2 - 1);
        cs_real_t kappa = -log(lamb + aux3);
        shape_param[0] = aux1/aux2 + lamb*kappa/(aux2*aux3);
        shape_param[1] = shape_param[0];
        shape_param[2] = -2./aux2 - 2.*lamb*kappa/(aux2*aux3);
        shape_param[3] = -2. * _a2 *lamb*kappa/aux3;
      }
      else if (lamb_m1 < -1e-10) {
        cs_real_t aux3 = sqrt(1. - aux2);
        cs_real_t kappa = acos(lamb);
        shape_param[0] = aux1/aux2+lamb*kappa/(-aux2*aux3);
        shape_param[1] = shape_param[0];
        shape_param[2] = -2./aux2 - 2.*lamb*kappa/(-aux2*aux3);
        shape_param[3] = 2. * _a2 * lamb*kappa/aux3;
      }
      else {
        shape_param[0] = 2.0/3.0;
        shape_param[1] = 2.0/3.0;
        shape_param[2] = 2.0/3.0;
        shape_param[3] = 2. * _a2;
      }

      if (shape == CS_LAGR_SHAPE_SPHEROID_STOC_MODEL) {
        /* Orientation */
        auto *orientation =
          cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                   CS_LAGR_ORIENTATION);
        auto *quaternion =
          cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                   CS_LAGR_QUATERNION);
        for (cs_lnum_t i = 0; i < 3; i++) {
          orientation[i] = zis->orientation[i];
        }

        /* Compute orientation from uniform orientation on a unit-sphere */
        cs_real_t theta0;
        cs_real_t phi0;
        cs_random_uniform(1, &theta0);
        cs_random_uniform(1, &phi0);
        theta0   = acos(2.0*theta0-1.0);
        phi0     = phi0*2.0*cs_math_pi;
        orientation[0] = sin(theta0)*cos(phi0);
        orientation[1] = sin(theta0)*sin(phi0);
        orientation[2] = cos(theta0);
        /* Initialize quaternions */
        quaternion[0] = 1.;
        quaternion[1] = 0.;
        quaternion[2] = 0.;
        quaternion[3] = 0.;
        /* TODO initialize other things */
      }

      if (shape == CS_LAGR_SHAPE_SPHEROID_JEFFERY_MODEL) {

        /* Euler parameters */
        auto *euler = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                               CS_LAGR_EULER);

        for (cs_lnum_t i = 0; i < 4; i++)
          euler[i] = zis->euler[i];

        /* Compute Euler angles
           (random orientation with a uniform distribution in [-1;1]) */
        cs_real_t trans_m[3][3];
        /* Generate the first two vectors */
        for (cs_lnum_t id = 0; id < 3; id++) {
          cs_random_uniform(1, &trans_m[id][0]); /* (?,0) */
          cs_random_uniform(1, &trans_m[id][1]); /* (?,1) */
          cs_random_uniform(1, &trans_m[id][2]); /* (?,2) */
          cs_real_3_t loc_vector =  {-1.+2*trans_m[id][0],
            -1.+2*trans_m[id][1],
            -1.+2*trans_m[id][2]};
          cs_real_t norm_trans_m = cs_math_3_norm( loc_vector );
          while ( norm_trans_m > 1 )
          {
            cs_random_uniform(1, &trans_m[id][0]); /* (?,0) */
            cs_random_uniform(1, &trans_m[id][1]); /* (?,1) */
            cs_random_uniform(1, &trans_m[id][2]); /* (?,2) */
            loc_vector[0] = -1.+2*trans_m[id][0];
            loc_vector[1] = -1.+2*trans_m[id][1];
            loc_vector[2] = -1.+2*trans_m[id][2];
            norm_trans_m = cs_math_3_norm( loc_vector );
          }
          for (cs_lnum_t id1 = 0; id1 < 3; id1++)
            trans_m[id][id1] = (-1.+2*trans_m[id][id1]) / norm_trans_m;
        }
        /* Correct 2nd vector (for perpendicularity to the 1st) */
        cs_real_3_t loc_vector0 =  {trans_m[0][0],
          trans_m[0][1],
          trans_m[0][2]};
        cs_real_3_t loc_vector1 =  {trans_m[1][0],
          trans_m[1][1],
          trans_m[1][2]};
        cs_real_t scal_prod = cs_math_3_dot_product(loc_vector0, loc_vector1);
        for (cs_lnum_t id = 0; id < 3; id++)
          trans_m[1][id] -= scal_prod * trans_m[0][id];
        /* Re-normalize */
        loc_vector1[0] = trans_m[1][0];
        loc_vector1[1] = trans_m[1][1];
        loc_vector1[2] = trans_m[1][2];
        cs_real_t norm_trans_m = cs_math_3_norm( loc_vector1 );
        for (cs_lnum_t id = 0; id < 3; id++)
          trans_m[1][id] /= norm_trans_m;

        /* Compute last vector (cross product of the two others) */
        loc_vector1[0] = trans_m[1][0];
        loc_vector1[1] = trans_m[1][1];
        loc_vector1[2] = trans_m[1][2];
        cs_real_3_t loc_vector2 =  {trans_m[2][0],
          trans_m[2][1],
          trans_m[2][2]};
        cs_math_3_cross_product( loc_vector0, loc_vector1, loc_vector2);
        for (cs_lnum_t id = 0; id < 3; id++)
          trans_m[2][id] = loc_vector2[id];

        /* Write Euler angles */
        cs_real_t random;
        cs_random_uniform(1, &random);
        if (random >= 0.5)
          euler[0] = pow(0.25*(trans_m[0][0]+trans_m[1][1]+trans_m[2][2]+1.),
                         0.5);
        else
          euler[0] = -pow(0.25*(trans_m[0][0]+trans_m[1][1]+trans_m[2][2]+1.),
                          0.5);
        euler[1] = 0.25 * (trans_m[2][1] - trans_m[1][2]) / euler[0];
        euler[2] = 0.25 * (trans_m[0][2] - trans_m[2][0]) / euler[0];
        euler[3] = 0.25 * (trans_m[1][0] - trans_m[0][1]) / euler[0];

        /* Compute initial angular velocity */

        /* Local reference frame */
        cs_real_t grad_vf_r[3][3];
        cs_math_33_transform_a_to_r(extra->grad_vel[c_id],
                                    trans_m,
                                    grad_vf_r);

        auto *ang_vel =
          cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                   CS_LAGR_ANGULAR_VEL);

        ang_vel[0] = 0.5*(grad_vf_r[2][1] - grad_vf_r[1][2]);
        ang_vel[1] = 0.5*(grad_vf_r[0][2] - grad_vf_r[2][0]);
        ang_vel[2] = 0.5*(grad_vf_r[0][1] - grad_vf_r[1][0]);
      }

    }

    /* Other parameters */
    cs_real_t diam = cs_lagr_particle_get_real(particle, p_am,
                                               CS_LAGR_DIAMETER);
    cs_real_t mporos = cs_glob_lagr_clogging_model->mporos;
    if (cs_glob_lagr_model->clogging == 1) {
      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_DIAMETER,
                                diam/(1.-mporos));
      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_HEIGHT, diam);
    }

    /* Other variables (mass, ...) depending on physical model  */
    cs_real_t d3 = pow(diam, 3.0);

    if (cs_glob_lagr_model->n_stat_classes > 0)
      cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_STAT_CLASS,
                                zis->cluster);

    if (   cs_glob_lagr_model->agglomeration == 1
        || cs_glob_lagr_model->fragmentation == 1) {
      cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_AGGLO_CLASS_ID,
                                zis->aggregat_class_id);
      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_AGGLO_FRACTAL_DIM,
                                zis->aggregat_fractal_dim);
    }

    /* used for 2nd order only */
    if (p_am->displ[0][CS_LAGR_TAUP_AUX] > 0)
      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_TAUP_AUX, 0.0);

    if (   cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_OFF
        || cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_HEAT) {

      if (cs_glob_lagr_model->clogging == 0)
        cs_lagr_particle_set_real(particle, p_am, CS_LAGR_MASS,
                                  zis->density * pis6 * d3);
      else
        cs_lagr_particle_set_real(particle, p_am, CS_LAGR_MASS,
                                  zis->density * pis6 * d3
                                  * pow(1.0-mporos, 3));

      if (   cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_HEAT
          && cs_glob_lagr_specific_physics->itpvar == 1) {

        if (zis->temperature_profile < 1) {
          if (cval_t != NULL)
            cs_lagr_particle_set_real(particle, p_am,
                                      CS_LAGR_FLUID_TEMPERATURE,
                                      cval_t[c_id] + tscl_shift);
        }

        /* constant temperature set, may be modified later by user function */
        else if (zis->temperature_profile == 1)
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_TEMPERATURE,
                                    zis->temperature);

        cs_lagr_particle_set_real(particle, p_am, CS_LAGR_CP,
                                  zis->cp);
        if (extra->radiative_model > 0)
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_EMISSIVITY,
                                    zis->emissivity);

      }

    }

    else if (cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_COAL) {

      int coal_id = zis->coal_number - 1;

      cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_COAL_ID, coal_id);
      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_FLUID_TEMPERATURE,
                                cval_t[c_id] + tscl_shift);

      auto *particle_temp
        = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am, CS_LAGR_TEMPERATURE);
      for (int ilayer = 0;
           ilayer < cs_glob_lagr_model->n_temperature_layers;
           ilayer++)
        particle_temp[ilayer] = zis->temperature;

      /* composition from DP_FCP */

      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_CP, cp2ch[coal_id]);

      cs_real_t mass = rho0ch[coal_id] * pis6 * d3;

      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_MASS, mass);
      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_WATER_MASS,
                                xwatch[coal_id] * mass);

      cs_real_t *particle_coal_mass
          = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                     CS_LAGR_COAL_MASS);
      cs_real_t *particle_coke_mass
        = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                   CS_LAGR_COKE_MASS);
      for (int ilayer = 0;
           ilayer < cs_glob_lagr_model->n_temperature_layers;
           ilayer++) {

        particle_coal_mass[ilayer]
          =    (1.0 - xwatch[coal_id]
                    - xashch[coal_id])
            * cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS)
            / cs_glob_lagr_model->n_temperature_layers;
        particle_coke_mass[ilayer] = 0.0;

      }

      cs_lagr_particle_set_real
        (particle, p_am,
         CS_LAGR_SHRINKING_DIAMETER,
         cs_lagr_particle_get_real(particle, p_am, CS_LAGR_DIAMETER));
      cs_lagr_particle_set_real
        (particle, p_am,
         CS_LAGR_INITIAL_DIAMETER,
         cs_lagr_particle_get_real(particle, p_am, CS_LAGR_DIAMETER));

      cs_real_t *particle_coal_density
        = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                   CS_LAGR_COAL_DENSITY);
      for (int ilayer = 0;
           ilayer < cs_glob_lagr_model->n_temperature_layers;
           ilayer++)
        particle_coal_density[ilayer] = rho0ch[coal_id];

    }

    /* Cooling tower model*/
    if (cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_CTWR) {

      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_MASS,
                                zis->density * pis6 * d3
                                );

      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_CP,
                                  zis->cp);

      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_TEMPERATURE,
                                cval_t_l[c_id]);

      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_FLUID_TEMPERATURE,
                                cval_t[c_id]+tscl_shift);
    }

    /* statistical weight */
    cs_lagr_particle_set_real(particle, p_am, CS_LAGR_STAT_WEIGHT,
                              zis->stat_weight);

    /* Fouling index */
    cs_lagr_particle_set_real(particle, p_am, CS_LAGR_FOULING_INDEX,
                              zis->fouling_index);

    /* Initialization of deposition model
     * And compute velocity fluctuation if deposition model is active */

    if (cs_glob_lagr_model->deposition == 1) {

      cs_real_t random;
      cs_random_uniform(1, &random);
      cs_lagr_particle_set_real(particle, p_am,
                                CS_LAGR_INTERF, 5.0 + 15.0 * random);
      cs_lagr_particle_set_real(particle, p_am,
                                CS_LAGR_YPLUS, 1000.0);
      cs_lagr_particle_set_lnum(particle, p_am,
                                CS_LAGR_MARKO_VALUE, -1);
      cs_lagr_particle_set_lnum(particle, p_am,
                                CS_LAGR_NEIGHBOR_FACE_ID, -1);

      /* Compute normalized wall-normal particle distance (y+) */

      cs_real_t yplus = 1000.0;
      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_YPLUS, yplus);

      for (cs_lnum_t il = ma->cell_b_faces_idx[c_id];
           il < ma->cell_b_faces_idx[c_id+1];
           il++) {

        cs_lnum_t face_id = ma->cell_b_faces[il];

        char b_type = bcs->elt_type[face_id];

        /* Test if the particle is located in a boundary cell */

        if (   b_type == CS_LAGR_DEPO1
            || b_type == CS_LAGR_DEPO2
            || b_type == CS_LAGR_DEPO_DLVO
            || b_type == CS_LAGR_REBOUND) {

          /* Calculation of the wall units  */

          cs_lnum_t  *neighbor_face_id;
          cs_real_t  *particle_yplus;

          neighbor_face_id
            = cs_lagr_particle_attr_get_ptr<cs_lnum_t>(particle, p_am,
                                                       CS_LAGR_NEIGHBOR_FACE_ID);
          particle_yplus
            = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                       CS_LAGR_YPLUS);

          cs_lagr_test_wall_cell(particle, p_am, visc_length,
                                 particle_yplus, neighbor_face_id);

        }
        else {
          cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_NEIGHBOR_FACE_ID, -1);
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_YPLUS, 0.);
        }

      }

      if (yplus < cs_lagr_particle_get_real(particle, p_am, CS_LAGR_INTERF)) {

        cs_lagr_particle_set_lnum(particle,
                                  p_am,
                                  CS_LAGR_MARKO_VALUE,
                                  CS_LAGR_COHERENCE_STRUCT_DEGEN_INNER_ZONE_DIFF);

      }

      else if (yplus > 100.0) {

        cs_lagr_particle_set_lnum(particle,
                                  p_am,
                                  CS_LAGR_MARKO_VALUE,
                                  CS_LAGR_COHERENCE_STRUCT_BULK);

      }

      else {

        cs_random_uniform(1, &random);

        if (random < 0.25)
          cs_lagr_particle_set_lnum(particle,
                                    p_am,
                                    CS_LAGR_MARKO_VALUE,
                                    CS_LAGR_COHERENCE_STRUCT_DEGEN_DIFFUSION);

        else if (random > 0.625)
          cs_lagr_particle_set_lnum(particle,
                                    p_am,
                                    CS_LAGR_MARKO_VALUE,
                                    CS_LAGR_COHERENCE_STRUCT_SWEEP);

        else /* if ((random > 0.25) && (random < 0.625)) */
          cs_lagr_particle_set_lnum(particle,
                                    p_am,
                                    CS_LAGR_MARKO_VALUE,
                                    CS_LAGR_COHERENCE_STRUCT_EJECTION);

      }

      if (yplus <= cs_lagr_particle_get_real(particle, p_am, CS_LAGR_INTERF)) {

        for (cs_lnum_t i = 0; i < 3; i++)
          vel_seen[i] = vel[c_id][i];

      }

      /* Initialization of additional "pointers"
       * for the resuspension model              */

      if (cs_glob_lagr_model->resuspension > 0) {

        cs_lagr_particle_set_real(particle, p_am, CS_LAGR_ADHESION_FORCE, 0.0);
        cs_lagr_particle_set_real(particle, p_am, CS_LAGR_ADHESION_TORQUE, 0.0);
        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_N_LARGE_ASPERITIES, 0);
        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_N_SMALL_ASPERITIES, 0);
        cs_lagr_particle_set_real(particle, p_am, CS_LAGR_DISPLACEMENT_NORM, 0.0);

      }

    }

    /* Initialization of clogging model */

    if (cs_glob_lagr_model->clogging == 1) {

      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_DEPO_TIME, 0.0);
      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_CONSOL_HEIGHT, 0.0);
      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_CLUSTER_NB_PART, 1.0);

    }

    /* Initialize the additional user variables */

    if (cs_glob_lagr_model->n_user_variables > 0) {
      cs_real_t  *user_attr
        = cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                   CS_LAGR_USER);
      for (int i = 0;
           i < cs_glob_lagr_model->n_user_variables;
           i++)
        user_attr[i] = 0.0;
    }

    cs_lagr_particles_current_to_previous(p_set, p_id);

  } /* End loop over new particles */


  /* Update weights to have the correct flow rate
     -------------------------------------------- */

  if (zis->flow_rate > 0.0 && zis->n_inject > 0) {

    cs_real_t dmass = 0.0;

    cs_lnum_t p_s_id = particle_range[0];
    cs_lnum_t p_e_id = particle_range[1];

    for (cs_lnum_t p_id = p_s_id; p_id < p_e_id; p_id++)
      dmass += cs_lagr_particles_get_real(p_set, p_id, CS_LAGR_MASS);

    cs_parall_sum(1, CS_REAL_TYPE, &dmass);

    /* Compute weights */

    if (dmass > 0.0) {
      cs_real_t s_weight =   zis->flow_rate * cs_glob_lagr_time_step->dtp
                           / dmass;
      for (cs_lnum_t p_id = p_s_id; p_id < p_e_id; p_id++)
        cs_lagr_particles_set_real(p_set, p_id, CS_LAGR_STAT_WEIGHT, s_weight);
    }

    else {

      char z_type_name[32] = "unknown";
      if (zis->location_id == CS_MESH_LOCATION_BOUNDARY_FACES)
        strncpy(z_type_name, _("boundary"), 31);
      else if (zis->location_id == CS_MESH_LOCATION_CELLS)
        strncpy(z_type_name, _("volume"), 31);
      z_type_name[31] = '\0';

      bft_error(__FILE__, __LINE__, 0,
                _("Lagrangian %s zone %d, set %d:\n"
                  " imposed flow rate is %g\n"
                  " while mass of injected particles is 0."),
                z_type_name, zis->zone_id, zis->set_id,
                (double)zis->flow_rate);

    }

  }

  BFT_FREE(_cval_t);
  BFT_FREE(vagaus);
  BFT_FREE(eig_vec);
  BFT_FREE(eig_val);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
