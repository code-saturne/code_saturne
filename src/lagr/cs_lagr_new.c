/*============================================================================
 * Handling of new particles.
 *============================================================================*/

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
#include "cs_halo.h"
#include "cs_log.h"
#include "cs_interface.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_adjacencies.h"
#include "cs_mesh_quantities.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_physical_model.h"
#include "cs_random.h"
#include "cs_search.h"
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
 *----------------------------------------------------------------------------*/

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

      cs_lagr_particles_set_lnum(particles, p_id, CS_LAGR_CELL_NUM, c_id+1);

      cs_real_t *part_coord
        = cs_lagr_particles_attr(particles, p_id, CS_LAGR_COORDS);

      _random_point_in_face(n_vertices,
                            vertex_ids,
                            (const cs_real_3_t *)mesh->vtx_coord,
                            fvq->b_face_cog + 3*face_id,
                            acc_surf_r,
                            part_coord);

      /* For safety, move particle slighty inside cell */

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
 *----------------------------------------------------------------------------*/

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

  cs_lnum_t *cell_face_idx = NULL, *cell_face_lst = NULL;

  cs_lagr_get_cell_face_connectivity(&cell_face_idx,
                                     &cell_face_lst);

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
    const cs_lnum_t n_cell_faces = cell_face_idx[cell_id+1] - cell_face_idx[cell_id];

    const cs_real_t *cell_cen = fvq->cell_cen + cell_id*3;

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

      const cs_lnum_t face_num = cell_face_lst[cell_face_idx[cell_id] + i];

      if (face_num > 0) { /* Interior face */

        face_id = face_num - 1;

        if (cell_id == mesh->i_face_cells[face_id][1])
          v_mult = -1;
        cs_lnum_t vtx_s = mesh->i_face_vtx_idx[face_id];
        n_vertices = mesh->i_face_vtx_idx[face_id+1] - vtx_s;
        vertex_ids = mesh->i_face_vtx_lst + vtx_s;
        face_cog = fvq->i_face_cog + (3*face_id);
        face_normal = fvq->i_face_normal + (3*face_id);

      }
      else { /* Boundary faces */

        assert(face_num < 0);

        face_id = -face_num - 1;

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

    for (cs_lnum_t i = 0; i < n_c_p; i++) {

      cs_lnum_t p_id = p_s_id + i;

      cs_lagr_particles_set_lnum(particles, p_id, CS_LAGR_CELL_NUM, cell_id+1);

      cs_real_t *part_coord
        = cs_lagr_particles_attr(particles, p_id, CS_LAGR_COORDS);

      /* search for matching center-to-face cone */

      cs_real_t r[2];
      cs_random_uniform(2, r);

      cs_lnum_t c_id = 0;
      while (c_id < n_cell_faces && r[0] > acc_vol_r[c_id])
        c_id++;

      cs_lnum_t face_id, n_vertices;
      const cs_lnum_t *vertex_ids;
      const cs_real_t *face_cog;

      const cs_lnum_t face_num = cell_face_lst[cell_face_idx[cell_id] + c_id];

      if (face_num > 0) { /* Interior face */

        face_id = face_num - 1;

        cs_lnum_t vtx_s = mesh->i_face_vtx_idx[face_id];
        n_vertices = mesh->i_face_vtx_idx[face_id+1] - vtx_s;
        vertex_ids = mesh->i_face_vtx_lst + vtx_s;
        face_cog = fvq->i_face_cog + (3*face_id);

      }
      else { /* Boundary faces */

        assert(face_num < 0);

        face_id = -face_num - 1;

        cs_lnum_t vtx_s = mesh->b_face_vtx_idx[face_id];
        n_vertices = mesh->b_face_vtx_idx[face_id+1] - vtx_s;
        vertex_ids = mesh->b_face_vtx_lst + vtx_s;
        face_cog = fvq->b_face_cog + (3*face_id);

      }

      _random_point_in_face(n_vertices,
                            vertex_ids,
                            (const cs_real_3_t *)mesh->vtx_coord,
                            face_cog,
                            acc_surf_r + cell_subface_index[c_id],
                            part_coord);

      /* In regular case, place point on segment joining cell center and
         point in cell; volume of truncated cone proportional to
         cube of distance along segment, so distribution compensates for this */

      if (fallback == false) {

        cs_real_t t = pow(r[1], 1./3.) * (1.0 - d_eps);
        for (cs_lnum_t j = 0; j < 3; j++)
          part_coord[j] += (cell_cen[j] - part_coord[j]) * (1. - t);
      }

      /* Move particle slighty towards cell center cell
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
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_new_particle_init(const cs_lnum_t  particle_range[2],
                          int              time_id,
                          const cs_real_t  visc_length[])
{
  cs_lagr_particle_set_t  *pset = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = pset->p_am;

  const cs_lagr_zone_data_t  *bcs = cs_glob_lagr_boundary_conditions;

  cs_lagr_extra_module_t  *extra = cs_get_lagr_extra_module();

  /* Map field arrays */

  const cs_real_3_t  *vel = NULL;
  const cs_real_t    *cvar_k = NULL;
  const cs_real_t    *cvar_r11 = NULL, *cvar_r22 = NULL, *cvar_r33 = NULL;
  const cs_real_6_t  *cvar_rij = NULL;

  vel = (const cs_real_3_t *)extra->vel->vals[time_id];

  /* Initialization */

  cs_real_t d2s3 = 2.0 / 3.0;

  /* Simulate instantaneous turbulent fluid velocities
     as seen by solid particles along their trajectory.
     -------------------------------------------------- */

  if (cs_glob_lagr_time_scheme->idistu == 1) {

    if (extra->cvar_k != NULL)
      cvar_k = (const cs_real_t *)extra->cvar_k->vals[time_id];

    else if (extra->cvar_rij != NULL)
      cvar_rij = (const cs_real_6_t *) extra->cvar_rij->vals[time_id];

    /* Deprecated irijco = 0 */
    else if (extra->cvar_r11 != NULL) {
      cvar_r11 = (const cs_real_t *)extra->cvar_r11->vals[time_id];
      cvar_r22 = (const cs_real_t *)extra->cvar_r22->vals[time_id];
      cvar_r33 = (const cs_real_t *)extra->cvar_r33->vals[time_id];
    }

    else {
      bft_error
        (__FILE__, __LINE__, 0,
         _("The Lagrangian module is incompatible with the selected\n"
           " turbulence model.\n\n"
           "Turbulent dispersion is used with:\n"
           "  cs_glob_lagr_time_scheme->idistu = %d\n"
           "And the turbulence model is iturb = %d\n\n"
           "The only turbulence models compatible with the Lagrangian model's\n"
           "turbulent dispersion are k-epsilon, Rij-epsilon, v2f, and k-omega."),
         cs_glob_lagr_time_scheme->idistu,
         extra->iturb);
    }

  }

  /* Random draws and computation of particle characteristic times */

  cs_lnum_t  n = particle_range[1] - particle_range[0];
  cs_real_3_t  *vagaus;

  BFT_MALLOC(vagaus, n, cs_real_3_t);

  if (cs_glob_lagr_time_scheme->idistu == 1 && n > 0) {
    cs_random_normal(n*3, (cs_real_t *)vagaus);
  }

  else {

    for (cs_lnum_t i = 0; i < n; i++) {
      vagaus[i][0] = 0.0;
      vagaus[i][1] = 0.0;
      vagaus[i][2] = 0.0;
    }

  }

  for (cs_lnum_t p_id = particle_range[0]; p_id < particle_range[1]; p_id++) {

    unsigned char *particle = pset->p_buffer + p_am->extents * p_id;

    cs_lnum_t iel  = cs_lagr_particle_get_cell_id(particle, p_am);
    cs_lnum_t l_id = p_id - particle_range[0];

    cs_real_t  *vel_seen
      = cs_lagr_particle_attr(particle, p_am, CS_LAGR_VELOCITY_SEEN);

    cs_real_t w = 0.;

    if (cs_glob_lagr_time_scheme->idistu == 1) {
      if (cvar_k != NULL)
        w = cvar_k[iel];
      else if (cvar_rij != NULL)
        w = 0.5*(cvar_rij[iel][0] + cvar_rij[iel][1] + cvar_rij[iel][2]);
      /* Deprecated irijco = 0 */
      else if (cvar_r11 != NULL)
        w = 0.5 * (cvar_r11[iel] + cvar_r22[iel] + cvar_r33[iel]);
    }

    cs_real_t tu = sqrt(d2s3 * w);

    for (cs_lnum_t i = 0; i < 3; i++)
      vel_seen[i] = vel[iel][i] + vagaus[l_id][i] * tu;

    cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_REBOUND_ID, -1);

    cs_lagr_particle_set_real(particle, p_am, CS_LAGR_TR_TRUNCATE, 0);

  }

  BFT_FREE(vagaus);

  /* Compute velcocity fluctuation if deposition model is active */

  if (cs_glob_lagr_model->deposition == 1) {

    const cs_mesh_adjacencies_t  *ma = cs_glob_mesh_adjacencies;

    for (cs_lnum_t p_id = particle_range[0]; p_id < particle_range[1]; p_id++) {

      unsigned char *particle = pset->p_buffer + p_am->extents * p_id;

      cs_lnum_t iel = cs_lagr_particle_get_cell_id(particle, p_am);

      /* Compute normalized wall-normal particle distance (y+) */

      cs_real_t yplus = 1000.0;
      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_YPLUS, yplus);

      for (cs_lnum_t il = ma->cell_b_faces_idx[iel];
           il < ma->cell_b_faces_idx[iel+1];
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
            = cs_lagr_particle_attr(particle, p_am, CS_LAGR_NEIGHBOR_FACE_ID);
          particle_yplus
            = cs_lagr_particle_attr(particle, p_am, CS_LAGR_YPLUS);

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

        cs_real_t random;
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

        cs_real_t *vel_seen
          = cs_lagr_particle_attr(particle, p_am, CS_LAGR_VELOCITY_SEEN);

        for (cs_lnum_t i = 0; i < 3; i++)
          vel_seen[i] = vel[iel][i];

      }

      /* No deposited particles at the injection */
      cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG,
                                CS_LAGR_PART_IN_FLOW);

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

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
