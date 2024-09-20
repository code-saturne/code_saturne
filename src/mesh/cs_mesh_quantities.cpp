/*============================================================================
 * Management of mesh quantities
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_array.h"
#include "cs_base.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_adjacencies.h"
#include "cs_porosity_from_scan.h"
#include "cs_post.h"

#include "fvm_nodal_from_desc.h"
#include "fvm_nodal_order.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*! \file cs_mesh_quantities.c
 *
 * \brief Management of mesh quantities.
 *
 * Please refer to the
 * <a href="../../theory.pdf#meshquantities"><b>geometric quantities</b></a>
 * section of the theory guide for more informations.
 */
/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to cs_mesh_quantities_t structure for the main mesh */

cs_mesh_quantities_t  *cs_glob_mesh_quantities = nullptr;

/* Choice of the algorithm for computing gravity centers of the cells */

static int _cell_cen_algorithm = 0;
static int _ajust_face_cog_compat_v11_v52 = 0;

/* Flag (mask) to activate bad cells correction
 * CS_BAD_CELLS_WARPED_CORRECTION
 * CS_FACE_DISTANCE_CLIP
 * are set as default options
 * */
unsigned cs_glob_mesh_quantities_flag =
  CS_BAD_CELLS_WARPED_CORRECTION + CS_FACE_DISTANCE_CLIP;

/* Number of computation updates */

static int _n_computations = 0;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Project solid vertices to a plane
 *
 * parameters:
 *   v0   <--  fluid vertex
 *   v1   <--  solid vertex
 *   nw   <--  plane normal
 *   pw   <--  point on the plane
 *   vout <->  projected vertex on the plane
 *----------------------------------------------------------------------------*/

static void
_proj_solid_vtx_to_plane(const cs_real_t  v0[3],
                         const cs_real_t  v1[3],
                         const cs_real_t  nw[3],
                         const cs_real_t  pw[3],
                         cs_real_t        vout[3])
{
  cs_real_t v01[3];

  for (cs_lnum_t i = 0; i < 3; i++)
    v01[i] = v1[i] - v0[i];

  cs_math_3_normalize(v01, v01);
  cs_real_t edgde_dot_nw = cs_math_3_dot_product(v01, nw);
  cs_real_t vone_dot_nw = CS_MAX(cs_math_3_dot_product(v1, nw), 0.);
  cs_real_t proj_factor; // = (  cs_math_fmax(cs_math_3_dot_product(v01, nw), 0.)
                         //     > cs_math_epzero) ?
  if (edgde_dot_nw <= 0)
    proj_factor = 0.;
  if (vone_dot_nw  > edgde_dot_nw)
    proj_factor = 1.;
  else
    proj_factor = vone_dot_nw / edgde_dot_nw;
    // cs_math_3_dot_product(v1, nw) /cs_math_3_dot_product(v01, nw) : 0.;

  for (cs_lnum_t i = 0; i < 3; i++)
    vout[i] = (v1[i] + pw[i]) - proj_factor * v01[i];
}

/*----------------------------------------------------------------------------
 * Build the geometrical matrix linear gradient correction
 *
 * parameters:
 *   m    <--  mesh
 *   fvq  <->  mesh quantities
 *----------------------------------------------------------------------------*/

static void
_compute_corr_grad_lin(const cs_mesh_t       *m,
                       cs_mesh_quantities_t  *fvq)
{
  /* Local variables */

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_with_ghosts = m->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = CS_MAX(m->n_b_faces, m->n_b_faces_all);

  const cs_alloc_mode_t amode = cs_alloc_mode_read_mostly;

  const cs_lnum_t  *b_face_cells = m->b_face_cells;
  const cs_lnum_2_t *restrict i_face_cells =
    (const cs_lnum_2_t *)m->i_face_cells;

  const int *restrict c_disable_flag = fvq->c_disable_flag;
  cs_lnum_t has_dc = fvq->has_disable_flag; /* Has cells disabled? */

  const cs_real_t *restrict cell_vol = fvq->cell_f_vol;
  const cs_real_3_t *restrict i_face_normal =
    (const cs_real_3_t *)fvq->i_f_face_normal;
  const cs_real_3_t *restrict b_face_normal =
    (const cs_real_3_t *)fvq->b_f_face_normal;
  const cs_real_3_t *restrict b_face_cog =
    (const cs_real_3_t *)fvq->b_f_face_cog;
  const cs_real_3_t *restrict i_face_cog =
    (const cs_real_3_t *)fvq->i_f_face_cog;

  if (fvq->corr_grad_lin_det == nullptr)
    BFT_MALLOC(fvq->corr_grad_lin_det, n_cells_with_ghosts, cs_real_t);
  if (fvq->corr_grad_lin == nullptr) {
    CS_MALLOC_HD(fvq->corr_grad_lin, n_cells_with_ghosts, cs_real_33_t, amode);
    cs_mem_advise_set_read_mostly(fvq->corr_grad_lin);
  }

  cs_real_t    *restrict corr_grad_lin_det = fvq->corr_grad_lin_det;
  cs_real_33_t *restrict corr_grad_lin     = fvq->corr_grad_lin;

  /* Initialization */
  for (cs_lnum_t cell_id = 0; cell_id < n_cells_with_ghosts; cell_id++) {
    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        corr_grad_lin[cell_id][i][j] = 0.;
    }
  }

  /* Internal faces contribution */
  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
    cs_lnum_t cell_id1 = i_face_cells[face_id][0];
    cs_lnum_t cell_id2 = i_face_cells[face_id][1];

    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        cs_real_t flux = i_face_cog[face_id][i] * i_face_normal[face_id][j];
        corr_grad_lin[cell_id1][i][j] += flux;
        corr_grad_lin[cell_id2][i][j] -= flux;
      }
    }
  }

  /* Boundary faces contribution */
  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    cs_lnum_t cell_id = b_face_cells[face_id];
    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        cs_real_t flux = b_face_cog[face_id][i] * b_face_normal[face_id][j];
        corr_grad_lin[cell_id][i][j] += flux;
      }
    }
  }

  /* Immersed boundaries contribution */
  if (fvq->c_w_face_cog != nullptr && fvq->c_w_face_normal != nullptr) {
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      for (cs_lnum_t i = 0; i < 3; i++) {
        for (cs_lnum_t j = 0; j < 3; j++) {
          cs_real_t flux = fvq->c_w_face_cog[cell_id*3+i]
                           * fvq->c_w_face_normal[cell_id*3+j];
          corr_grad_lin[cell_id][i][j] += flux;
        }
      }
    }
  }

  /* Matrix inversion */
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    cs_real_t dvol;
    /* Is the cell disabled (for solid or porous)? Not the case if coupled */
    if (has_dc * c_disable_flag[has_dc * cell_id] == 0)
      dvol = 1. / cell_vol[cell_id];
    else
      dvol = 0.;
    double cocg11 = corr_grad_lin[cell_id][0][0] * dvol;
    double cocg12 = corr_grad_lin[cell_id][1][0] * dvol;
    double cocg13 = corr_grad_lin[cell_id][2][0] * dvol;
    double cocg21 = corr_grad_lin[cell_id][0][1] * dvol;
    double cocg22 = corr_grad_lin[cell_id][1][1] * dvol;
    double cocg23 = corr_grad_lin[cell_id][2][1] * dvol;
    double cocg31 = corr_grad_lin[cell_id][0][2] * dvol;
    double cocg32 = corr_grad_lin[cell_id][1][2] * dvol;
    double cocg33 = corr_grad_lin[cell_id][2][2] * dvol;

    double a11 = cocg22 * cocg33 - cocg32 * cocg23;
    double a12 = cocg32 * cocg13 - cocg12 * cocg33;
    double a13 = cocg12 * cocg23 - cocg22 * cocg13;
    double a21 = cocg31 * cocg23 - cocg21 * cocg33;
    double a22 = cocg11 * cocg33 - cocg31 * cocg13;
    double a23 = cocg21 * cocg13 - cocg11 * cocg23;
    double a31 = cocg21 * cocg32 - cocg31 * cocg22;
    double a32 = cocg31 * cocg12 - cocg11 * cocg32;
    double a33 = cocg11 * cocg22 - cocg21 * cocg12;

    double det_inv = cocg11 * a11 + cocg21 * a12 + cocg31 * a13;

    if (fabs(det_inv) >= 1.e-15) {
      det_inv = 1. / det_inv;

      corr_grad_lin[cell_id][0][0] = a11 * det_inv;
      corr_grad_lin[cell_id][0][1] = a12 * det_inv;
      corr_grad_lin[cell_id][0][2] = a13 * det_inv;
      corr_grad_lin[cell_id][1][0] = a21 * det_inv;
      corr_grad_lin[cell_id][1][1] = a22 * det_inv;
      corr_grad_lin[cell_id][1][2] = a23 * det_inv;
      corr_grad_lin[cell_id][2][0] = a31 * det_inv;
      corr_grad_lin[cell_id][2][1] = a32 * det_inv;
      corr_grad_lin[cell_id][2][2] = a33 * det_inv;

      double a1 = corr_grad_lin[cell_id][0][0];
      double a2 = corr_grad_lin[cell_id][0][1];
      double a3 = corr_grad_lin[cell_id][0][2];
      double a4 = corr_grad_lin[cell_id][1][0];
      double a5 = corr_grad_lin[cell_id][1][1];
      double a6 = corr_grad_lin[cell_id][1][2];
      double a7 = corr_grad_lin[cell_id][2][0];
      double a8 = corr_grad_lin[cell_id][2][1];
      double a9 = corr_grad_lin[cell_id][2][2];

      double determinant =  a1 * (a5*a9 - a8*a6)
                          - a2 * (a4*a9 - a7*a6)
                          + a3 * (a4*a8 - a7*a5);

      corr_grad_lin_det[cell_id] = determinant;
    }
    else {
      corr_grad_lin[cell_id][0][0] = 0.;
      corr_grad_lin[cell_id][0][1] = 0.;
      corr_grad_lin[cell_id][0][2] = 0.;
      corr_grad_lin[cell_id][1][0] = 0.;
      corr_grad_lin[cell_id][1][1] = 0.;
      corr_grad_lin[cell_id][1][2] = 0.;
      corr_grad_lin[cell_id][2][0] = 0.;
      corr_grad_lin[cell_id][2][1] = 0.;
      corr_grad_lin[cell_id][2][2] = 0.;

      corr_grad_lin_det[cell_id] = 1.;
    }
  }

  if (m->halo != nullptr) {
    cs_halo_sync_var(m->halo, CS_HALO_STANDARD, corr_grad_lin_det);
    cs_halo_sync_var_strided(m->halo, CS_HALO_STANDARD,
                             (cs_real_t *)corr_grad_lin, 9);
    /* TODO handle rotational periodicity */
  }
}

/*----------------------------------------------------------------------------
 * Compute quantities associated to faces (border or internal)
 *
 * parameters:
 *   n_faces         <--  number of faces
 *   vtx_coord       <--  vertex coordinates
 *   face_vtx_idx    <--  "face -> vertices" connectivity index
 *   face_vtx        <--  "face -> vertices" connectivity
 *   face_cog        -->  coordinates of the center of gravity of the faces
 *   face_normal     -->  face surface normals
 *
 *                          Pi+1
 *              *---------*                   B  : barycenter of the polygon
 *             / .       . \
 *            /   .     .   \                 Pi : vertices of the polygon
 *           /     .   .     \
 *          /       . .  Ti   \               Ti : triangle
 *         *.........B.........* Pi
 *     Pn-1 \       . .       /
 *           \     .   .     /
 *            \   .     .   /
 *             \ .   T0  . /
 *              *---------*
 *            P0
 *----------------------------------------------------------------------------*/

static void
_compute_face_quantities(cs_lnum_t        n_faces,
                         const cs_real_t  vtx_coord[][3],
                         const cs_lnum_t  face_vtx_idx[],
                         const cs_lnum_t  face_vtx[],
                         cs_real_t        face_cog[][3],
                         cs_real_t        face_normal[][3])
{
  /* Checking */

  assert(face_cog != nullptr || n_faces == 0);
  assert(face_normal != nullptr || n_faces == 0);

  const cs_real_t one_third = 1./3.;
  const cs_real_t s_epsilon = 1.e-32; /* TODO define better "zero" threshold */

  /* Loop on faces */

# pragma omp parallel for  if (n_faces > CS_THR_MIN)
  for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++) {

    /* Define the polygon (P) according to the vertices (Pi) of the face */

    cs_lnum_t s_id = face_vtx_idx[f_id];
    cs_lnum_t e_id = face_vtx_idx[f_id + 1];

    cs_lnum_t n_face_vertices = e_id - s_id;

    if (n_face_vertices == 3) {
      const cs_lnum_t v0 = face_vtx[s_id];
      const cs_lnum_t v1 = face_vtx[s_id+1];
      const cs_lnum_t v2 = face_vtx[s_id+2];
      cs_real_t v01[3], v02[3], vn[3];
      for (cs_lnum_t i = 0; i < 3; i++)
        face_cog[f_id][i] = one_third * (  vtx_coord[v0][i]
                                         + vtx_coord[v1][i]
                                         + vtx_coord[v2][i]);
      for (cs_lnum_t i = 0; i < 3; i++)
        v01[i] = vtx_coord[v1][i] - vtx_coord[v0][i];
      for (cs_lnum_t i = 0; i < 3; i++)
        v02[i] = vtx_coord[v2][i] - vtx_coord[v0][i];
      cs_math_3_cross_product(v01, v02, vn);
      for (cs_lnum_t i = 0; i < 3; i++)
        face_normal[f_id][i] = 0.5*vn[i];
    }

    else { /* For non-triangle faces, assume a division into triangles
              joining edges and an approximate face center */

      /* Compute approximate face center coordinates for the polygon */

      cs_real_t a_center[3] = {0, 0, 0};
      cs_real_t f_center[3] = {0, 0, 0}, f_norm[3] = {0, 0, 0};

      for (cs_lnum_t j = s_id; j < e_id; j++) {
        const cs_lnum_t v0 = face_vtx[j];
        for (cs_lnum_t i = 0; i < 3; i++)
          a_center[i] += vtx_coord[v0][i];
      }

      for (cs_lnum_t i = 0; i < 3; i++)
        a_center[i] /= n_face_vertices;

      cs_real_t d_surf_d3 = 0.;

      /* In most cases, the following 2 loops can be merged into a single loop,
         but for very bad quality faces, some sub-triangles could be oriented
         differently, so we use 2 passes for safety. */

      if (n_face_vertices < 8) { /* version with local caching for most cases */

        cs_real_t vc0[3], vc1[3], vn[8][3], vtc[8][3];

        /* First pass (face normal) */

        for (cs_lnum_t tri_id = 0; tri_id < n_face_vertices; tri_id++) {

          const cs_lnum_t v0 = face_vtx[s_id + tri_id];
          const cs_lnum_t v1 = face_vtx[s_id + (tri_id+1)%n_face_vertices];

          for (cs_lnum_t i = 0; i < 3; i++) {
            vc0[i] = vtx_coord[v0][i] - a_center[i];
            vc1[i] = vtx_coord[v1][i] - a_center[i];
            /* center in the relative reference frame (shifted by a_center) */
            vtc[tri_id][i] = vc0[i] + vc1[i];
          }

          cs_math_3_cross_product(vc0, vc1, vn[tri_id]);

          for (cs_lnum_t i = 0; i < 3; i++)
            f_norm[i] += 0.5 * vn[tri_id][i];

        }

        /* Second pass (face center) */

        cs_real_t sum_w = cs_math_3_norm(f_norm);
        cs_real_t inv_norm = 1.;
        if (sum_w > s_epsilon) {
          inv_norm = 1. / sum_w;
          d_surf_d3 = one_third * inv_norm;
        }
        cs_real_t n[3];
        for (int i = 0; i < 3; i++)
          n[i] = inv_norm * f_norm[i];

        for (cs_lnum_t tri_id = 0; tri_id < n_face_vertices; tri_id++) {

          /* Projected surface of the triangle in the normal direction */
          cs_real_t w = 0.5 * cs_math_3_dot_product(vn[tri_id], n);

          for (cs_lnum_t i = 0; i < 3; i++)
            f_center[i] += w*vtc[tri_id][i];

        }

      }
      else  { /* generic version */

        cs_real_t vc0[3], vc1[3], vn[3], vtc[3];

        /* First pass (face normal) */

        for (cs_lnum_t tri_id = 0; tri_id < n_face_vertices; tri_id++) {

          const cs_lnum_t v0 = face_vtx[s_id + tri_id];
          const cs_lnum_t v1 = face_vtx[s_id + (tri_id+1)%n_face_vertices];

          for (cs_lnum_t i = 0; i < 3; i++) {
            vc0[i] = vtx_coord[v0][i] - a_center[i];
            vc1[i] = vtx_coord[v1][i] - a_center[i];
          }

          cs_math_3_cross_product(vc0, vc1, vn);

          for (cs_lnum_t i = 0; i < 3; i++)
            f_norm[i] += 0.5 * vn[i];

        }

        /* Second pass (face center) */

        cs_real_t sum_w = cs_math_3_norm(f_norm);
        cs_real_t inv_norm = 1.;
        if (sum_w > s_epsilon) {
          inv_norm = 1. / sum_w;
          d_surf_d3 = one_third * inv_norm;
        }
        cs_real_t n[3];
        for (int i = 0; i < 3; i++)
          n[i] = inv_norm * f_norm[i];

        for (cs_lnum_t tri_id = 0; tri_id < n_face_vertices; tri_id++) {

          const cs_lnum_t v0 = face_vtx[s_id + tri_id];
          const cs_lnum_t v1 = face_vtx[s_id + (tri_id+1)%n_face_vertices];

          for (cs_lnum_t i = 0; i < 3; i++) {
            vc0[i] = vtx_coord[v0][i] - a_center[i];
            vc1[i] = vtx_coord[v1][i] - a_center[i];
            /* center in the relative reference frame (shifted by a_center) */
            vtc[i] = vc0[i] + vc1[i];
          }

          cs_math_3_cross_product(vc0, vc1, vn);

          /* Projected surface of the triangle in the normal direction */
          cs_real_t w = 0.5 * cs_math_3_dot_product(vn, n);

          for (cs_lnum_t i = 0; i < 3; i++)
            f_center[i] += w*vtc[i];

        }

      }

      for (cs_lnum_t i = 0; i < 3; i++)
        face_normal[f_id][i] = f_norm[i];

      for (cs_lnum_t i = 0; i < 3; i++)
        face_cog[f_id][i] = a_center[i] + d_surf_d3 * f_center[i];

    } /* end of test on triangle */

  } /* end of loop on faces */
}

/*----------------------------------------------------------------------------
 * Compute quantities associated to faces (border or internal)
 *
 * parameters:
 *   n_faces         <--  number of faces
 *   vtx_coord       <--  vertex coordinates
 *   face_vtx_idx    <--  "face -> vertices" connectivity index
 *   face_vtx_lst    <--  "face -> vertices" connectivity list
 *   face_normal     -->  surface normal of the face
 *
 *
 *                          Pi+1
 *              *---------*                   B  : barycenter of the polygon
 *             / .       . \
 *            /   .     .   \                 Pi : vertices of the polygon
 *           /     .   .     \
 *          /       . .  Ti   \               Ti : triangle
 *         *.........B.........* Pi
 *     Pn-1 \       . .       /
 *           \     .   .     /
 *            \   .     .   /
 *             \ .   T0  . /
 *              *---------*
 *            P0
 *----------------------------------------------------------------------------*/

static void
_compute_face_normal(cs_lnum_t         n_faces,
                     const cs_real_t   vtx_coord[][3],
                     const cs_lnum_t   face_vtx_idx[],
                     const cs_lnum_t   face_vtx[],
                     cs_real_t         face_normal[][3])
{
  /* Checking */

  assert(face_normal != nullptr || n_faces == 0);

  /* Loop on faces */

# pragma omp parallel for  if (n_faces > CS_THR_MIN)
  for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++) {

    /* Define the polygon (P) according to the vertices (Pi) of the face */

    cs_lnum_t s_id = face_vtx_idx[f_id];
    cs_lnum_t e_id = face_vtx_idx[f_id + 1];

    cs_lnum_t n_face_vertices = e_id - s_id;

    if (n_face_vertices == 3) {
      const cs_lnum_t v0 = face_vtx[s_id];
      const cs_lnum_t v1 = face_vtx[s_id+1];
      const cs_lnum_t v2 = face_vtx[s_id+2];
      cs_real_t v01[3], v02[3], vn[3];
      for (cs_lnum_t i = 0; i < 3; i++)
        v01[i] = vtx_coord[v1][i] - vtx_coord[v0][i];
      for (cs_lnum_t i = 0; i < 3; i++)
        v02[i] = vtx_coord[v2][i] - vtx_coord[v0][i];
      cs_math_3_cross_product(v01, v02, vn);
      for (cs_lnum_t i = 0; i < 3; i++)
        face_normal[f_id][i] = 0.5*vn[i];
    }

    else {

      /* Compute approximate face center coordinates for the polygon */

      cs_real_t a_center[3] = {0, 0, 0};
      cs_real_t f_norm[3] = {0, 0, 0};

      for (cs_lnum_t j = s_id; j < e_id; j++) {
        const cs_lnum_t v0 = face_vtx[j];
        for (cs_lnum_t i = 0; i < 3; i++)
          a_center[i] += vtx_coord[v0][i];
      }

      for (cs_lnum_t i = 0; i < 3; i++)
        a_center[i] /= n_face_vertices;

      /* loop on edges, with implied subdivision into triangles
         defined by edge and cell center */

      cs_real_t vc0[3], vc1[3], vn[3];

      for (cs_lnum_t tri_id = 0; tri_id < n_face_vertices; tri_id++) {

        const cs_lnum_t v0 = face_vtx[s_id + tri_id];
        const cs_lnum_t v1 = face_vtx[s_id + (tri_id+1)%n_face_vertices];

        for (cs_lnum_t i = 0; i < 3; i++) {
          vc0[i] = vtx_coord[v0][i] - a_center[i];
          vc1[i] = vtx_coord[v1][i] - a_center[i];
        }

        cs_math_3_cross_product(vc0, vc1, vn);

        for (cs_lnum_t i = 0; i < 3; i++)
          f_norm[i] += vn[i];

      }

      for (cs_lnum_t i = 0; i < 3; i++)
        face_normal[f_id][i] = 0.5*f_norm[i];

    } /* end of test on triangle */

  } /* end of loop on faces */
}

/*----------------------------------------------------------------------------
 * Compute face surfaces based on face norms.
 *
 * parameters:
 *   n_faces         <--  number of faces
 *   face_norm       <--  face surface normals
 *   face_surf       -->  face surfaces
 *----------------------------------------------------------------------------*/

static void
_compute_face_surface(cs_lnum_t        n_faces,
                      const cs_real_t  face_norm[],
                      cs_real_t        face_surf[])
{
# pragma omp parallel for  if (n_faces > CS_THR_MIN)
  for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++)
    face_surf[f_id] = cs_math_3_norm(face_norm + f_id*3);
}

/*----------------------------------------------------------------------------
 * Adjust the position of face centers to obtain a contribution to cell
 * volumes consistent with that obtained using a splitting of the face
 * into triangles using face edges and an approxmate center.
 *
 * This adjustment was done by default from code_saturne 1.1 to 5.2.
 * To handle warped faces, a correction of the cell center is used so
 * that the contribution to the cell volume using Green's theorem is
 * consistent with the volume computed for sub-faces.
 *
 * It is based on using the coordinates origin instead of the cell centers
 * to compute the correction, whereas the volume computation uses the cell
 * center for better rounding behavior, and uses the cell vertices center
 * instead of the computed face center, so the consistency gain with
 * warped faces (where this is supposed to be useful) would need to be
 * checked.
 *
 * parameters:
 *   n_faces         <--  number of faces
 *   vtx_coord       <--  vertex coordinates
 *   face_vtx_idx    <--  "face -> vertices" connectivity index
 *   face_vtx        <--  "face -> vertices" connectivity list
 *   face_cog        <->  coordinates of the center of gravity of the faces
 *   face_norm       <--  face surface normals
 *
 *                          Pi+1
 *              *---------*                   B  : barycenter of the polygon
 *             / .       . \
 *            /   .     .   \                 Pi : vertices of the polygon
 *           /     .   .     \
 *          /       . .  Ti   \               Ti : triangle
 *         *.........B.........* Pi
 *     Pn-1 \       . .       /
 *           \     .   .     /
 *            \   .     .   /
 *             \ .   T0  . /
 *              *---------*
 *            P0
 *----------------------------------------------------------------------------*/

static void
_adjust_face_cog_v11_v52(cs_lnum_t         n_faces,
                         const cs_real_t   vtx_coord[][3],
                         const cs_lnum_t   face_vtx_idx[],
                         const cs_lnum_t   face_vtx[],
                         cs_real_t         face_cog[][3],
                         const cs_real_t   face_norm[][3])
{
  const cs_real_t one_third = 1./3.;

  /* Loop on faces
   --------------- */

  for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++) {

    cs_real_t tri_vol_part = 0.;

    /* Define the polygon (P) according to the vertices (Pi) of the face */

    cs_lnum_t s_id = face_vtx_idx[f_id];
    cs_lnum_t e_id = face_vtx_idx[f_id + 1];

    cs_lnum_t n_face_vertices = e_id - s_id;

    /* No warping - related correction required for triangles*/
    if (n_face_vertices < 4)
      continue;

    /* Compute approximate face center coordinates for the polygon */

    cs_real_t a_center[3] = {0, 0, 0}, f_center[3] = {0, 0, 0};

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      const cs_lnum_t v0 = face_vtx[j];
      for (cs_lnum_t i = 0; i < 3; i++)
        a_center[i] += vtx_coord[v0][i];
    }

    for (cs_lnum_t i = 0; i < 3; i++)
      a_center[i] /= n_face_vertices;

    /* loop on edges, with implied subdivision into triangles
       defined by edge and cell center */

    cs_real_t vc0[3], vc1[3], vn[3], vtc[3];
    cs_real_t sum_w = 0;

    for (cs_lnum_t tri_id = 0; tri_id < n_face_vertices; tri_id++) {

      const cs_lnum_t v0 = face_vtx[s_id + tri_id];
      const cs_lnum_t v1 = face_vtx[s_id + (tri_id+1)%n_face_vertices];

      for (cs_lnum_t i = 0; i < 3; i++) {
        vc0[i] = vtx_coord[v0][i] - a_center[i];
        vc1[i] = vtx_coord[v1][i] - a_center[i];
        vtc[i] = vtx_coord[v1][i] + vtx_coord[v0][i] + a_center[i];
      }

      cs_math_3_cross_product(vc0, vc1, vn);

      tri_vol_part += one_third * 0.5 * cs_math_3_dot_product(vtc, vn);

      cs_real_t w = 0.5 * cs_math_3_norm(vn);

      if (cs_math_3_dot_product(vn, face_norm[f_id]) < 0.0)
        w *= -1.0;

      sum_w += w;

      for (cs_lnum_t i = 0; i < 3; i++)
        f_center[i] += w * vtc[i];

    } /* End of loop on triangles of the face */

    /*
     * Compute the center of gravity G(P) of the polygon P, and
     * compute the part of volume of the polygon (before rectification)
     *
     *  -->    ->
     *  OG(P).N(P)
     *------------ */

    cs_real_t face_vol_part = 0.0;

    for (cs_lnum_t i = 0; i < 3; i++) {
      f_center[i] *= one_third / sum_w;
      face_vol_part += (f_center[i] * face_norm[f_id][i]);
    }

    cs_real_t rectif_cog =    (tri_vol_part - face_vol_part)
                            / (sum_w * sum_w);

    for (cs_lnum_t i = 0; i < 3; i++)
      face_cog[f_id][i] = f_center[i] + rectif_cog * face_norm[f_id][i];

  } /* End of loop on faces */
}

/*----------------------------------------------------------------------------
 * Refine face center computation based on initial value.
 *
 * This may be useful for warped faces; for plance faces, the face center
 * should remain identical.
 *
 * parameters:
 *   n_faces         <--  number of faces
 *   vtx_coord       <--  vertex coordinates
 *   face_vtx_idx    <--  "face -> vertices" connectivity index
 *   face_vtx        <--  "face -> vertices" connectivity
 *   face_cog        <->  coordinates of the center of gravity of the faces
 *   face_norm       <--  face surface normals
 *
 *                          Pi+1
 *              *---------*                   B  : barycenter of the polygon
 *             / .       . \
 *            /   .     .   \                 Pi : vertices of the polygon
 *           /     .   .     \
 *          /       . .  Ti   \               Ti : triangle
 *         *.........B.........* Pi
 *     Pn-1 \       . .       /
 *           \     .   .     /
 *            \   .     .   /
 *             \ .   T0  . /
 *              *---------*
 *            P0
 *----------------------------------------------------------------------------*/

static void
_refine_warped_face_centers(cs_lnum_t          n_faces,
                            const cs_real_3_t  vtx_coord[],
                            const cs_lnum_t    face_vtx_idx[],
                            const cs_lnum_t    face_vtx[],
                            cs_real_t         face_cog[][3],
                            const cs_real_t   face_norm[][3])
{
  /* Checking */

  assert(face_cog != nullptr || n_faces == 0);

  const cs_real_t one_third = 1./3.;
  const cs_real_t s_epsilon = 1.e-32; /* TODO define better "zero" threshold */

  /* Loop on faces */

  for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++) {

    /* Define the polygon (P) according to the vertices (Pi) of the face */

    cs_lnum_t s_id = face_vtx_idx[f_id];
    cs_lnum_t e_id = face_vtx_idx[f_id + 1];

    cs_lnum_t n_face_vertices = e_id - s_id;

    if (n_face_vertices > 3) {

      cs_real_t ref_size = sqrt(cs_math_3_norm(face_norm[f_id]));

      for (int ite = 0; ite < 5; ite++) {

        /* Use previous face center coordinates for the polygon */

        cs_real_t a_center[3];
        cs_real_t f_center[3] = {0, 0, 0};

        for (cs_lnum_t i = 0; i < 3; i++)
          a_center[i] = face_cog[f_id][i];

        /* loop on edges, with implied subdivision into triangles
           defined by edge and cell center */

        cs_real_t vc0[3], vc1[3], vn[3], vtc[3];
        cs_real_t sum_w = 0;

        for (cs_lnum_t tri_id = 0; tri_id < n_face_vertices; tri_id++) {

          const cs_lnum_t v0 = face_vtx[s_id + tri_id];
          const cs_lnum_t v1 = face_vtx[s_id + (tri_id+1)%n_face_vertices];

          for (cs_lnum_t i = 0; i < 3; i++) {
            vc0[i] = vtx_coord[v0][i] - a_center[i];
            vc1[i] = vtx_coord[v1][i] - a_center[i];
            vtc[i] = vtx_coord[v1][i] + vtx_coord[v0][i] + a_center[i];
          }

          cs_math_3_cross_product(vc0, vc1, vn);

          cs_real_t w = cs_math_3_norm(vn);

          if (cs_math_3_dot_product(vn, face_norm[f_id]) < 0.0)
            w *= -1.0;

          sum_w += w;

          for (cs_lnum_t i = 0; i < 3; i++)
            f_center[i] += w*vtc[i];

        }

        if (sum_w > s_epsilon) {
          for (cs_lnum_t i = 0; i < 3; i++)
            face_cog[f_id][i] = one_third * f_center[i]/sum_w;
          if (cs_math_3_distance(face_cog[f_id], a_center) / ref_size < 1e-8)
            break;
        }
        else
          break;

      } /* end of face iterations */

    } /* end of test on triangle */

  } /* end of loop on faces */
}

/*----------------------------------------------------------------------------
 * Recompute quantities associated to faces (border or internal) when
 * the quality of the mesh is not good enough
 *----------------------------------------------------------------------------*/

static void
_correct_cell_face_center(const cs_mesh_t    *mesh,
                          cs_lnum_t           n_cells_with_ghosts,
                          cs_lnum_t           n_i_faces,
                          cs_lnum_t           n_b_faces,
                          const cs_lnum_2_t   i_face_cells[],
                          const cs_lnum_t     b_face_cells[],
                          cs_real_3_t         cell_cen[],
                          cs_real_3_t         i_face_cog[],
                          cs_real_3_t         b_face_cog[],
                          cs_real_3_t         i_face_normal[],
                          cs_real_3_t         b_face_normal[])
{
  int nitmax = 500;
  cs_real_3_t *i_face_cog0, *b_face_cog0;
  cs_real_3_t *i_face_cen, *b_face_cen;

  cs_real_t *relaxf;
  cs_real_t *relaxb;

  cs_real_33_t *dxidxj;
  cs_real_t *determinant;

  BFT_MALLOC(i_face_cog0, n_i_faces, cs_real_3_t);
  BFT_MALLOC(b_face_cog0, n_b_faces, cs_real_3_t);
  BFT_MALLOC(i_face_cen, n_i_faces, cs_real_3_t);
  BFT_MALLOC(b_face_cen, n_b_faces, cs_real_3_t);
  BFT_MALLOC(relaxf, n_i_faces, cs_real_t);
  BFT_MALLOC(relaxb, n_b_faces, cs_real_t);
  BFT_MALLOC(dxidxj, n_cells_with_ghosts, cs_real_33_t);
  BFT_MALLOC(determinant, n_cells_with_ghosts, cs_real_t);

  /* Iterative process */
  for (int sweep = 0; sweep < 1; sweep++) {

    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++)
      for (cs_lnum_t i = 0; i < 3; i++)
        i_face_cog0[face_id][i] = i_face_cog[face_id][i];

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++)
      for (cs_lnum_t i = 0; i < 3; i++)
        b_face_cog0[face_id][i] = b_face_cog[face_id][i];

    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
      cs_lnum_t cell_id1 = i_face_cells[face_id][0];
      cs_lnum_t cell_id2 = i_face_cells[face_id][1];

      /* IF0 . S */
      double ps1 = cs_math_3_distance_dot_product(cell_cen[cell_id1],
                                                  i_face_cog0[face_id],
                                                  i_face_normal[face_id]);

      /* IJ . S */
      double ps2 = cs_math_3_distance_dot_product(cell_cen[cell_id1],
                                                  cell_cen[cell_id2],
                                                  i_face_normal[face_id]);

      double lambda = 0.5;
      if (CS_ABS(ps2) > 1.e-20)
        lambda = ps1 / ps2;

      lambda = CS_MAX(lambda, 1./3.);
      lambda = CS_MIN(lambda, 2./3.);

      /* F = I + lambda * IJ */
      for (cs_lnum_t i = 0; i < 3; i++)
        i_face_cen[face_id][i] = cell_cen[cell_id1][i]
                               + lambda * (  cell_cen[cell_id2][i]
                                           - cell_cen[cell_id1][i]);
    }

    /* Compute the projection of I on the boundary face: b_face_cen */
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      cs_lnum_t cell_id = b_face_cells[face_id];

      cs_real_t normal[3];
      /* Normal is vector 0 if the b_face_normal norm is too small */
      cs_math_3_normalize(b_face_normal[face_id], normal);

      if (cs_math_3_norm(normal) > 0.) {
        double lambda = cs_math_3_distance_dot_product(cell_cen[cell_id],
                                                       b_face_cog[face_id],
                                                       normal);

        for (cs_lnum_t i = 0; i < 3; i++)
          b_face_cen[face_id][i] = cell_cen[cell_id][i] + lambda * normal[i];
      }
      else {
        for (cs_lnum_t i = 0; i < 3; i++)
          b_face_cen[face_id][i] =  b_face_cog[face_id][i];
      }
    }

    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++)
      relaxf[face_id] = 1.;
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++)
      relaxb[face_id] = 1.;

    int iiter = 0;
    cs_gnum_t irelax = 0;

    do
    {
      iiter +=1;

      for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
        for (cs_lnum_t i = 0; i < 3; i++)
          i_face_cog[face_id][i]
            = (1. - relaxf[face_id]) * i_face_cog0[face_id][i]
                  + relaxf[face_id]  * i_face_cen[face_id][i];
      }

      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
        for (cs_lnum_t i = 0; i < 3; i++)
          b_face_cog[face_id][i]
            = (1. - relaxb[face_id]) * b_face_cog0[face_id][i]
                  + relaxb[face_id]  * b_face_cen[face_id][i];
      }

      for (cs_lnum_t cell_id = 0; cell_id < n_cells_with_ghosts; cell_id++) {
        for (cs_lnum_t i = 0; i < 3; i++) {
          for (cs_lnum_t j = 0; j < 3; j++)
            dxidxj[cell_id][i][j] = 0.;
        }
      }

      for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
        cs_lnum_t cell_id1 = i_face_cells[face_id][0];
        cs_lnum_t cell_id2 = i_face_cells[face_id][1];

        for (cs_lnum_t i = 0; i < 3; i++)
          for (cs_lnum_t j = 0; j < 3; j++) {
            double fluxij = i_face_cog[face_id][i] //TODO minus celli
              * i_face_normal[face_id][j];
            dxidxj[cell_id1][i][j] += fluxij;
            dxidxj[cell_id2][i][j] -= fluxij;
          }
      }

      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
        cs_lnum_t cell_id = b_face_cells[face_id];

        for (cs_lnum_t i = 0; i < 3; i++)
          for (cs_lnum_t j = 0; j < 3; j++) {
            double fluxij =  b_face_cog[face_id][i]
                           * b_face_normal[face_id][j];
            dxidxj[cell_id][i][j] += fluxij;
          }
      }

      for (cs_lnum_t cell_id = 0; cell_id <  mesh->n_cells; cell_id++) {
        double vol = (  dxidxj[cell_id][0][0]
                      + dxidxj[cell_id][1][1]
                      + dxidxj[cell_id][2][2]) / 3.;

        //FIXME
        if (vol >= 0)
          vol = CS_MAX(vol,1.e-20);

        double a1 = dxidxj[cell_id][0][0] / vol;
        double a2 = dxidxj[cell_id][0][1] / vol;
        double a3 = dxidxj[cell_id][0][2] / vol;
        double a4 = dxidxj[cell_id][1][0] / vol;
        double a5 = dxidxj[cell_id][1][1] / vol;
        double a6 = dxidxj[cell_id][1][2] / vol;
        double a7 = dxidxj[cell_id][2][0] / vol;
        double a8 = dxidxj[cell_id][2][1] / vol;
        double a9 = dxidxj[cell_id][2][2] / vol;

        determinant[cell_id] = fabs(  a1 * (a5*a9 - a8*a6)
                                    - a2 * (a4*a9 - a7*a6)
                                    + a3 * (a4*a8 - a7*a5));

        //FIXME
        determinant[cell_id] = CS_MAX(determinant[cell_id],1.e-20);
        determinant[cell_id] = CS_MIN(determinant[cell_id],1./determinant[cell_id]);

        /* M-matrix structure control */
        double mmatrice1a = a1 - CS_ABS(a2) - CS_ABS(a3);
        double mmatrice1b = a1 - CS_ABS(a4) - CS_ABS(a7);
        double mmatrice2a = a5 - CS_ABS(a4) - CS_ABS(a6);
        double mmatrice2b = a5 - CS_ABS(a2) - CS_ABS(a8);
        double mmatrice3a = a9 - CS_ABS(a7) - CS_ABS(a8);
        double mmatrice3b = a9 - CS_ABS(a6) - CS_ABS(a3);

        if (   mmatrice1a <= 0. || mmatrice1b <= 0.
            || mmatrice2a <= 0. || mmatrice2b <= 0.
            || mmatrice3a <= 0. || mmatrice3b <= 0.)
          determinant[cell_id] = 0.;
      }

      if (mesh->halo != nullptr)
        cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, determinant);

      irelax = 0;

      //FIXME test was 0.001
      cs_real_t threshold = 0.1;
      for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
        cs_lnum_t cell_id1 = i_face_cells[face_id][0];
        cs_lnum_t cell_id2 = i_face_cells[face_id][1];

        if (   determinant[cell_id1] < threshold
            || determinant[cell_id2] < threshold) {
          irelax +=1;
          relaxf[face_id] *= 0.95;
        }
      }

      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
        cs_lnum_t cell_id = b_face_cells[face_id];

        if (determinant[cell_id] < threshold) {
          irelax += 1;
          relaxb[face_id] *= 0.95;
        }
      }

      cs_parall_counter(&irelax, 1);

    } while (iiter < nitmax && irelax > 0);

  }
  BFT_FREE(i_face_cog0);
  BFT_FREE(b_face_cog0);
  BFT_FREE(i_face_cen);
  BFT_FREE(b_face_cen);
  BFT_FREE(relaxf);
  BFT_FREE(relaxb);
  BFT_FREE(dxidxj);
  BFT_FREE(determinant);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell centers and volumes.
 *
 * \param[in]   mesh         pointer to mesh structure
 * \param[in]   i_face_norm  surface normal of internal faces
 * \param[in]   i_face_cog   center of gravity of internal faces
 * \param[in]   b_face_norm  surface normal of border faces
 * \param[in]   b_face_cog   center of gravity of border faces
 * \param[out]  cell_cen     cell centers
 * \param[out]  cell_vol     cell volumes
 */
/*----------------------------------------------------------------------------*/

static void
_compute_cell_quantities(const cs_mesh_t      *mesh,
                         const cs_real_3_t     i_face_norm[],
                         const cs_real_3_t     i_face_cog[],
                         const cs_real_3_t     b_face_norm[],
                         const cs_real_3_t     b_face_cog[],
                         cs_real_3_t *restrict cell_cen,
                         cs_real_t   *restrict cell_vol)
{
  /* Mesh connectivity */

  const  cs_lnum_t  n_i_faces = mesh->n_i_faces;
  const  cs_lnum_t  n_b_faces = CS_MAX(mesh->n_b_faces, mesh->n_b_faces_all);
  const  cs_lnum_t  n_cells = mesh->n_cells;
  const  cs_lnum_t  n_cells_ext = mesh->n_cells_with_ghosts;
  const  cs_lnum_2_t  *i_face_cells
    = (const cs_lnum_2_t *)(mesh->i_face_cells);
  const  cs_lnum_t  *b_face_cells = mesh->b_face_cells;

  /* Checking */

  assert(cell_cen != nullptr);
  assert(cell_vol != nullptr);

  /* Compute approximate cell center using face centers */

  cs_real_3_t *a_cell_cen;
  BFT_MALLOC(a_cell_cen, n_cells_ext, cs_real_3_t);

  cs_mesh_quantities_cell_faces_cog(mesh,
                                    (const cs_real_t *)i_face_norm,
                                    (const cs_real_t *)i_face_cog,
                                    (const cs_real_t *)b_face_norm,
                                    (const cs_real_t *)b_face_cog,
                                    (cs_real_t *)a_cell_cen);

  /* Initialization */

  for (cs_lnum_t j = 0; j < n_cells_ext; j++)
    cell_vol[j] = 0.;

  for (cs_lnum_t j = 0; j < n_cells_ext; j++) {
    for (cs_lnum_t i = 0; i < 3; i++)
      cell_cen[j][i] = 0.;
  }

  /* Loop on interior faces
     ---------------------- */

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    /* For each cell sharing the internal face, we update
     * cell_cen and cell_area */

    cs_lnum_t c_id1 = i_face_cells[f_id][0];
    cs_lnum_t c_id2 = i_face_cells[f_id][1];

    /* Implicit subdivision of cell into face vertices-cell-center pyramids */

    if (c_id1 > -1) {
      cs_real_t pyra_vol_3 = cs_math_3_distance_dot_product(a_cell_cen[c_id1],
                                                            i_face_cog[f_id],
                                                            i_face_norm[f_id]);
      for (cs_lnum_t i = 0; i < 3; i++)
        cell_cen[c_id1][i] += pyra_vol_3 *(  0.75*i_face_cog[f_id][i]
                                           + 0.25*a_cell_cen[c_id1][i]);
      cell_vol[c_id1] += pyra_vol_3;
    }
    if (c_id2 > -1) {
      cs_real_t pyra_vol_3 = cs_math_3_distance_dot_product(i_face_cog[f_id],
                                                            a_cell_cen[c_id2],
                                                            i_face_norm[f_id]);
      for (cs_lnum_t i = 0; i < 3; i++)
        cell_cen[c_id2][i] += pyra_vol_3 *(  0.75*i_face_cog[f_id][i]
                                           + 0.25*a_cell_cen[c_id2][i]);
      cell_vol[c_id2] += pyra_vol_3;
    }

  } /* End of loop on interior faces */

  /* Loop on boundary faces
     --------------------- */

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    /* For each cell sharing a border face, we update the numerator
     * of cell_cen and cell_area */

    const cs_lnum_t c_id1 = b_face_cells[f_id];

    /* Computation of the area of the face
       (note that c_id1 == -1 may happen for isolated faces,
       which are cleaned afterwards) */

    if (c_id1 > -1) {
      cs_real_t pyra_vol_3 = cs_math_3_distance_dot_product(a_cell_cen[c_id1],
                                                            b_face_cog[f_id],
                                                            b_face_norm[f_id]);
      for (cs_lnum_t i = 0; i < 3; i++)
        cell_cen[c_id1][i] += pyra_vol_3 *(  0.75*b_face_cog[f_id][i]
                                           + 0.25*a_cell_cen[c_id1][i]);
      cell_vol[c_id1] += pyra_vol_3;
    }

  } /* End of loop on boundary faces */

  BFT_FREE(a_cell_cen);

  /* Loop on cells to finalize the computation
     ----------------------------------------- */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    for (cs_lnum_t i = 0; i < 3; i++)
      cell_cen[c_id][i] /= cell_vol[c_id];

    cell_vol[c_id] /= 3.0;

  }
}

/*----------------------------------------------------------------------------*
 * Compute new cell centers by minimizing the distance to faces
 *
 * parameters:
 *   mesh           <--  pointer to mesh structure
 *   i_face_normal  <--  surface normal of internal faces
 *   i_face_cog     <--  center of gravity of internal faces
 *   b_face_normal  <--  surface normal of border faces
 *   b_face_cog     <--  center of gravity of border faces
 *   cell_cen       -->  center of gravity of cells
 *----------------------------------------------------------------------------*/

static void
_recompute_cell_cen_face(const cs_mesh_t     *mesh,
                         const cs_real_3_t   i_face_normal[],
                         const cs_real_3_t   i_face_cog[],
                         const cs_real_3_t   b_face_normal[],
                         const cs_real_3_t   b_face_cog[],
                         cs_real_3_t         cell_cen[])
{
  const  cs_lnum_t  n_i_faces = mesh->n_i_faces;
  const  cs_lnum_t  n_b_faces = CS_MAX(mesh->n_b_faces, mesh->n_b_faces_all);

  const  cs_lnum_t  n_cells_with_ghosts = mesh->n_cells_with_ghosts;

  const  cs_lnum_2_t  *i_face_cells
    = (const cs_lnum_2_t *)(mesh->i_face_cells);
  const  cs_lnum_t  *b_face_cells = mesh->b_face_cells;

  /* First pass of verification */
  int *pb1;
  BFT_MALLOC(pb1, n_cells_with_ghosts, int);

  for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells_with_ghosts; cell_id++)
    pb1[cell_id] = 0;

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

    cs_lnum_t cell_id1 = i_face_cells[face_id][0];
    cs_lnum_t cell_id2 = i_face_cells[face_id][1];

    /* IF . S */
    double psi1 = cs_math_3_distance_dot_product(cell_cen[cell_id1],
                                                 i_face_cog[face_id],
                                                 i_face_normal[face_id]);
    /* JF . S */
    double psj1 = cs_math_3_distance_dot_product(cell_cen[cell_id2],
                                                 i_face_cog[face_id],
                                                 i_face_normal[face_id]);
    if (psi1 < 0.)
      pb1[cell_id1]++;
    if (psj1 > 0.)
      pb1[cell_id2]++;
  }

  cs_gnum_t cpt1 = 0;
  for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells; cell_id++)
    if (pb1[cell_id] > 0) cpt1++;
  cs_parall_counter(&cpt1, 1);

  if (cpt1 > 0) {
    bft_printf("Total number of cell centers on the other side of a face "
               "(before correction) = %llu / %ld\n",
               (unsigned long long)cpt1,
               (long)mesh->n_cells);

    /* Second pass */
    cs_real_33_t *a;
    cs_real_3_t  *b;
    cs_real_3_t  *cdgbis;

    BFT_MALLOC(a, n_cells_with_ghosts, cs_real_33_t);
    BFT_MALLOC(b, n_cells_with_ghosts, cs_real_3_t);
    BFT_MALLOC(cdgbis, n_cells_with_ghosts, cs_real_3_t);

    /* init matrice et second membre */
    for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells_with_ghosts; cell_id++) {
      for (cs_lnum_t i = 0; i < 3; i++) {
        b[cell_id][i] = 0.;
        for (cs_lnum_t j = 0; j < 3; j++)
          a[cell_id][i][j] = 0.;
      }
    }

    /* Contribution from interior faces */
    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

      const cs_lnum_t cell_id1 = i_face_cells[face_id][0];
      const cs_lnum_t cell_id2 = i_face_cells[face_id][1];
      const double surfn = cs_math_3_norm(i_face_normal[face_id]);

      /* If the surface of the face is < 1.e-20 computation is not done,
       * to avoid dividing by 0 and unnecessary operations.
       */
      if (!(cs_glob_mesh_quantities_flag & CS_FACE_NULL_SURFACE) ||
          surfn  > 1.e-20) {
        for (cs_lnum_t i = 0; i < 3; i++)
          for (cs_lnum_t j = 0; j < 3; j++) {
            a[cell_id1][i][j] +=   i_face_normal[face_id][i]
                                 * i_face_normal[face_id][j] / surfn;
            a[cell_id2][i][j] +=   i_face_normal[face_id][i]
                                 * i_face_normal[face_id][j] / surfn;
          }

        double ps = cs_math_3_dot_product(i_face_normal[face_id],
                                          i_face_cog[face_id]);

        for (cs_lnum_t i = 0; i < 3; i++) {
          b[cell_id1][i] += ps * i_face_normal[face_id][i] / surfn;
          b[cell_id2][i] += ps * i_face_normal[face_id][i] / surfn;
        }
      }

    }

    /* Contribution from boundary faces */
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

      cs_lnum_t cell_id = b_face_cells[face_id];
      double surfn = cs_math_3_norm(b_face_normal[face_id]);

      /* If the surface of the face is < 1.e-20 computation is not done,
       * to avoid dividing by 0 and unnecessary operations.
       */
      if (!(cs_glob_mesh_quantities_flag & CS_FACE_NULL_SURFACE) ||
          surfn > 1.e-20) {
        for (cs_lnum_t i = 0; i < 3; i++)
          for (cs_lnum_t j = 0; j < 3; j++) {
            a[cell_id][i][j] +=   b_face_normal[face_id][i]
                                * b_face_normal[face_id][j] / surfn;
          }

        double ps = cs_math_3_dot_product(b_face_normal[face_id],
                                          b_face_cog[face_id]);

        for (cs_lnum_t i = 0; i < 3; i++) {
          b[cell_id][i] += ps * b_face_normal[face_id][i] / surfn;
        }
      }

    }

    /* invert system */
    double aainv[3][3];
    double bb[3];
    for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells; cell_id++) {

      cdgbis[cell_id][0] = cell_cen[cell_id][0];
      cdgbis[cell_id][1] = cell_cen[cell_id][1];
      cdgbis[cell_id][2] = cell_cen[cell_id][2];

      if (pb1[cell_id] > 0) {

        double adim = a[cell_id][0][0] + a[cell_id][1][1] + a[cell_id][2][2];

        if (adim > 0.) {
          bb[0] = b[cell_id][0] / adim;
          bb[1] = b[cell_id][1] / adim;
          bb[2] = b[cell_id][2] / adim;

          /* Matrix inversion */
          double cocg11 = a[cell_id][0][0] / adim;
          double cocg12 = a[cell_id][0][1] / adim;
          double cocg13 = a[cell_id][0][2] / adim;
          double cocg21 = a[cell_id][1][0] / adim;
          double cocg22 = a[cell_id][1][1] / adim;
          double cocg23 = a[cell_id][1][2] / adim;
          double cocg31 = a[cell_id][2][0] / adim;
          double cocg32 = a[cell_id][2][1] / adim;
          double cocg33 = a[cell_id][2][2] / adim;

          double a11 = cocg22 * cocg33 - cocg32 * cocg23;
          double a12 = cocg32 * cocg13 - cocg12 * cocg33;
          double a13 = cocg12 * cocg23 - cocg22 * cocg13;
          double a21 = cocg31 * cocg23 - cocg21 * cocg33;
          double a22 = cocg11 * cocg33 - cocg31 * cocg13;
          double a23 = cocg21 * cocg13 - cocg11 * cocg23;
          double a31 = cocg21 * cocg32 - cocg31 * cocg22;
          double a32 = cocg31 * cocg12 - cocg11 * cocg32;
          double a33 = cocg11 * cocg22 - cocg21 * cocg12;

          double det_inv = cocg11 * a11 + cocg21 * a12 + cocg31 * a13;

          if (CS_ABS(det_inv) >= 1.e-15) {
            det_inv = 1. / det_inv;

            aainv[0][0] = a11 * det_inv;
            aainv[0][1] = a12 * det_inv;
            aainv[0][2] = a13 * det_inv;
            aainv[1][0] = a21 * det_inv;
            aainv[1][1] = a22 * det_inv;
            aainv[1][2] = a23 * det_inv;
            aainv[2][0] = a31 * det_inv;
            aainv[2][1] = a32 * det_inv;
            aainv[2][2] = a33 * det_inv;

            for (cs_lnum_t i = 0; i < 3; i++)
              cdgbis[cell_id][i] =   aainv[i][0] * bb[0]
                                   + aainv[i][1] * bb[1]
                                   + aainv[i][2] * bb[2];
          }
        }
      }
    }

    if (mesh->halo != nullptr) {
      cs_halo_sync_var_strided(mesh->halo, CS_HALO_EXTENDED,
                               (cs_real_t *)cdgbis, 3);
      if (mesh->n_init_perio > 0)
        cs_halo_perio_sync_coords(mesh->halo, CS_HALO_EXTENDED,
                                  (cs_real_t *)cdgbis);
    }

    /* Second verification */

    int *pb2;
    BFT_MALLOC(pb2, n_cells_with_ghosts, int);

    for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells_with_ghosts; cell_id++)
      pb2[cell_id] = 0;

    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

      const cs_lnum_t cell_id1 = i_face_cells[face_id][0];
      const cs_lnum_t cell_id2 = i_face_cells[face_id][1];

      /* IF . S */
      double psi1 = cs_math_3_distance_dot_product(cdgbis[cell_id1],
                                                   i_face_cog[face_id],
                                                   i_face_normal[face_id]);
      /* JF . S */
      double psj1 = cs_math_3_distance_dot_product(cdgbis[cell_id2],
                                                   i_face_cog[face_id],
                                                   i_face_normal[face_id]);
      if (psi1 < 0.)
        pb2[cell_id1]++;
      if (psj1 > 0.)
        pb2[cell_id2]++;
    }

    cs_gnum_t cpt2 = 0;
    for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells; cell_id++)
      if (pb2[cell_id] > 0)
        cpt2++;
    cs_parall_counter(&cpt2, 1);

    bft_printf("Total number of cell centers on the other side of a face "
               "(after correction) = %llu / %ld\n",
               (unsigned long long)cpt2,
               (long)mesh->n_cells);

    for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells; cell_id++) {
      if (pb1[cell_id] > 0 && pb2[cell_id] == 0) {
        cell_cen[cell_id][0] = cdgbis[cell_id][0];
        cell_cen[cell_id][1] = cdgbis[cell_id][1];
        cell_cen[cell_id][2] = cdgbis[cell_id][2];
      }
    }

    if (mesh->halo != nullptr) {
      cs_halo_sync_var_strided(mesh->halo, CS_HALO_EXTENDED,
                               (cs_real_t *)cell_cen, 3);
      if (mesh->n_init_perio > 0)
        cs_halo_perio_sync_coords(mesh->halo, CS_HALO_EXTENDED,
                                  (cs_real_t *)cell_cen);
    }

    BFT_FREE(a);
    BFT_FREE(b);
    BFT_FREE(cdgbis);
    BFT_FREE(pb2);
  }

  /* Free memory */
  BFT_FREE(pb1);
}

/*----------------------------------------------------------------------------
 * Compute the volume of cells C from their n faces F(i) and their center of
 * gravity G(Fi) where i=0, n-1
 *
 *         1    n-1
 *  G(C) = - .  Sum  Surf(Fi) G(Fi)
 *         3    i=0
 *
 * parameters:
 *   mesh           <--  pointer to mesh structure
 *   i_face_norm    <--  surface normal of internal faces
 *   i_face_cog     <--  center of gravity of internal faces
 *   b_face_norm    <--  surface normal of border faces
 *   b_face_cog     <--  center of gravity of border faces
 *   cell_cen       <--  center of gravity of cells
 *   cell_vol       -->  cells volume
 *----------------------------------------------------------------------------*/

static void
_compute_cell_volume(const cs_mesh_t   *mesh,
                     const cs_real_3_t  i_face_norm[],
                     const cs_real_3_t  i_face_cog[],
                     const cs_real_3_t  b_face_norm[],
                     const cs_real_3_t  b_face_cog[],
                     const cs_real_3_t  cell_cen[],
                     cs_real_t          cell_vol[])
{
  const cs_lnum_t  n_b_faces = CS_MAX(mesh->n_b_faces, mesh->n_b_faces_all);
  const cs_real_t  a_third = 1.0/3.0;

  /* Initialization */

  for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells_with_ghosts; cell_id++)
    cell_vol[cell_id] = 0;

  /* Loop on internal faces */

  for (cs_lnum_t fac_id = 0; fac_id < mesh->n_i_faces; fac_id++) {

    const cs_lnum_t cell_id1 = mesh->i_face_cells[fac_id][0];
    const cs_lnum_t cell_id2 = mesh->i_face_cells[fac_id][1];

    cell_vol[cell_id1] += cs_math_3_distance_dot_product(cell_cen[cell_id1],
                                                         i_face_cog[fac_id],
                                                         i_face_norm[fac_id]);
    cell_vol[cell_id2] -= cs_math_3_distance_dot_product(cell_cen[cell_id2],
                                                         i_face_cog[fac_id],
                                                         i_face_norm[fac_id]);

  }

  /* Loop on border faces */

  for (cs_lnum_t fac_id = 0; fac_id < n_b_faces; fac_id++) {

    const cs_lnum_t cell_id1 = mesh->b_face_cells[fac_id];

    cell_vol[cell_id1] += cs_math_3_distance_dot_product(cell_cen[cell_id1],
                                                         b_face_cog[fac_id],
                                                         b_face_norm[fac_id]);
  }

  /* First computation of the volume */

  for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells; cell_id++)
    cell_vol[cell_id] *= a_third;
}

/*----------------------------------------------------------------------------
 * Correction of small or negative volumes.
 *
 * Based on smoothing with neighbors; does not conserve the total volume.
 *----------------------------------------------------------------------------*/

static void
_cell_bad_volume_correction(const cs_mesh_t   *mesh,
                            cs_real_t          cell_vol[])
{
  if (mesh->halo != nullptr)
    cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, cell_vol);

  /* Iterations in order to get vol_I / max(vol_J) > critmin */

  double *vol_neib_max;
  BFT_MALLOC(vol_neib_max, mesh->n_cells_with_ghosts, double);

  for (int iter = 0; iter < 10; iter++) {

    for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells_with_ghosts; cell_id++)
      vol_neib_max[cell_id] = 0.;

    for (cs_lnum_t fac_id = 0; fac_id < mesh->n_i_faces; fac_id++) {

      const cs_lnum_t cell_id1 = mesh->i_face_cells[fac_id][0];
      const cs_lnum_t cell_id2 = mesh->i_face_cells[fac_id][1];
      const cs_real_t vol1 = cell_vol[cell_id1];
      const cs_real_t vol2 = cell_vol[cell_id2];

      if (vol2 > 0.)
        vol_neib_max[cell_id1] = CS_MAX(vol_neib_max[cell_id1], vol2);

      if (vol1 > 0.)
        vol_neib_max[cell_id2] = CS_MAX(vol_neib_max[cell_id2], vol1);
    }

    /* Previous value of 0.2 sometimes leads to computation divergence */
    /* 0.01 seems better and safer for the moment */
    const cs_real_t critmin = 0.01;

    for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells; cell_id++)
      cell_vol[cell_id] = CS_MAX(cell_vol[cell_id],
                                 critmin * vol_neib_max[cell_id]);

    if (mesh->halo != nullptr)
      cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, cell_vol);
  }

  BFT_FREE(vol_neib_max);
}

/*----------------------------------------------------------------------------
 * Compute the total, min, and max volumes of cells.
 *
 * parameters:
 *   mesh           <--  pointer to mesh structure
 *   cell_vol       <--  cells volume
 *   min_vol        -->  minimum control volume
 *   max_vol        -->  maximum control volume
 *   tot_vol        -->  total   control volume
 *----------------------------------------------------------------------------*/

static void
_cell_volume_reductions(const cs_mesh_t  *mesh,
                        cs_real_t         cell_vol[],
                        cs_real_t        *min_vol,
                        cs_real_t        *max_vol,
                        cs_real_t        *tot_vol)
{
  *min_vol =  cs_math_infinite_r;
  *max_vol = -cs_math_infinite_r;
  *tot_vol = 0.;

  for (cs_lnum_t cell_id = 0; cell_id < mesh->n_cells; cell_id++) {
    *min_vol = CS_MIN(*min_vol, cell_vol[cell_id]);
    *max_vol = CS_MAX(*max_vol, cell_vol[cell_id]);
    *tot_vol = *tot_vol + cell_vol[cell_id];
  }
}

/*----------------------------------------------------------------------------
 * Compute some distances relative to faces and associated weighting.
 *
 * parameters:
 *   n_i_faces      <--  number of interior faces
 *   n_b_faces      <--  number of border  faces
 *   i_face_cells   <--  interior "faces -> cells" connectivity
 *   b_face_cells   <--  border "faces -> cells" connectivity
 *   i_face_u_norm  <--  unit normal of interior faces
 *   i_face_norm    <--  surface normal of interior faces
 *   b_face_u_norm  <--  unit normal of border faces
 *   b_face_norm    <--  surface normal of border faces
 *   i_face_cog     <--  center of gravity of interior faces
 *   b_face_cog     <--  center of gravity of border faces
 *   cell_cen       <--  cell center
 *   cell_vol       <--  cell volume
 *   i_dist         -->  distance IJ.Nij for interior faces
 *   b_dist         -->  likewise for border faces
 *   weight         -->  weighting factor (Aij=pond Ai+(1-pond)Aj)
 *----------------------------------------------------------------------------*/

static void
_compute_face_distances(cs_lnum_t        n_i_faces,
                        cs_lnum_t        n_b_faces,
                        const cs_lnum_t  i_face_cells[][2],
                        const cs_lnum_t  b_face_cells[],
                        const cs_real_t  i_face_u_normal[][3],
                        const cs_real_t  i_face_normal[][3],
                        const cs_real_t  b_face_u_normal[][3],
                        const cs_real_t  b_face_normal[][3],
                        const cs_real_t  i_face_cog[][3],
                        const cs_real_t  b_face_cog[][3],
                        const cs_real_t  cell_cen[][3],
                        const cs_real_t  cell_vol[],
                        cs_real_t        i_dist[],
                        cs_real_t        b_dist[],
                        cs_real_t        weight[])
{
  cs_gnum_t w_count[2] = {0};

  /* Interior faces */

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

    const cs_real_t *u_normal = i_face_u_normal[face_id];

    const cs_lnum_t cell_id1 = i_face_cells[face_id][0];
    const cs_lnum_t cell_id2 = i_face_cells[face_id][1];

    /* Distance between the neighbor cell centers
     * and dot-product with the normal */
    i_dist[face_id] = cs_math_3_distance_dot_product(cell_cen[cell_id1],
                                                     cell_cen[cell_id2],
                                                     u_normal);

    if (CS_ABS(i_dist[face_id]) > 1e-20) {
      /* Distance between the face center of gravity
         and the neighbor cell center
         and dot-product with the normal */
      cs_real_t dist2f = cs_math_3_distance_dot_product(i_face_cog[face_id],
                                                        cell_cen[cell_id2],
                                                        u_normal);
      weight[face_id] = dist2f / i_dist[face_id];
    }
    else {
      weight[face_id] = 0.5;
    }

    /* Clipping of cell cell distances */
    if (cs_glob_mesh_quantities_flag & CS_FACE_DISTANCE_CLIP) {

      /* Min value between IJ and
       * (Omega_i+Omega_j)/S_ij which is exactly the distance for tetras */

      cs_real_t face_normal_norm = cs_math_3_norm(i_face_normal[face_id]);
      cs_real_t distmax = -1;

      /* If CS_FACE_NULL_SURFACE is used, only update distmax value
       * if the face surface is not 0.
       */
      if (   !(cs_glob_mesh_quantities_flag & CS_FACE_NULL_SURFACE)
          && face_normal_norm > 1.e-20)
        distmax = cs_math_fmin(cs_math_3_distance(cell_cen[cell_id1],
                                                  cell_cen[cell_id2]),
                               (  (cell_vol[cell_id1] + cell_vol[cell_id2])
                                / face_normal_norm));

      /* Previous value of 0.2 sometimes leads to computation divergence */
      /* 0.01 seems better and safer for the moment */
      const cs_real_t critmin = 0.01;
      if (i_dist[face_id] < critmin * distmax) {
        w_count[0] += 1;
        i_dist[face_id] = cs_math_fmax(i_dist[face_id], critmin * distmax);
      }

      /* Clippings due to null surface */
      if (cs_glob_mesh_quantities_flag & CS_FACE_NULL_SURFACE) {
        if (face_normal_norm <= 1.e-20)
          i_dist[face_id] = cs_math_3_distance(cell_cen[cell_id1],
                                               cell_cen[cell_id2]);

        /* i_dist is arbitrary clipped to 1. in this case.
         * FIXME: This value could be improved in future releases
         */
        if (i_dist[face_id] < 1.e-20)
          i_dist[face_id] = 1.;
      }

      /* Clipping of weighting */
      weight[face_id] = cs_math_fmax(weight[face_id], 0.001);
      weight[face_id] = cs_math_fmin(weight[face_id], 0.999);
    }
  }

  /* Boundary faces */

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

    const cs_real_t *normal = b_face_u_normal[face_id];

    const cs_lnum_t cell_id = b_face_cells[face_id];

    /* Distance between the face center of gravity
       and the neighbor cell center */
    b_dist[face_id] = cs_math_3_distance_dot_product(cell_cen[cell_id],
                                                     b_face_cog[face_id],
                                                     normal);
    /* Clipping of cell boundary distances */
    if (cs_glob_mesh_quantities_flag & CS_FACE_DISTANCE_CLIP) {

      /* Min value between IF and
       * (Omega_i)/S which is exactly the distance for tetrahedra */
      const cs_real_t face_normal_norm = cs_math_3_norm(b_face_normal[face_id]);
      cs_real_t distmax = -1;

      /* If CS_FACE_NULL_SURFACE is used, only update distmax value
       * if the face surface is not 0.
       */
      if (   !(cs_glob_mesh_quantities_flag & CS_FACE_NULL_SURFACE)
          && face_normal_norm > 1.e-20) {
        distmax = fmin(cs_math_3_distance(cell_cen[cell_id],
                                          b_face_cog[face_id]),
                       cell_vol[cell_id]/face_normal_norm);
      }

      double critmin = 0.01;
      if (b_dist[face_id] < critmin * distmax) {
        w_count[1] += 1;
        b_dist[face_id] = fmax(b_dist[face_id], critmin * distmax);
      }

      /* Clippings due to null surface */
      if (cs_glob_mesh_quantities_flag & CS_FACE_NULL_SURFACE) {
        if (face_normal_norm <= 1.e-20)
          b_dist[face_id] = cs_math_3_distance(cell_cen[cell_id],
                                               b_face_cog[face_id]);

        /* b_dist is arbirary clipped to 1. in this case.
         * FIXME: This value could be improved in future releases
         */
        if (b_dist[face_id] < 1.e-20)
          b_dist[face_id] = 1.;
      }

    }
  }

  cs_parall_counter(w_count, 2);

  if (w_count[0] > 0)
    bft_printf(_("\n"
                 "%llu faces have a too small distance between centers.\n"
                 "For these faces, the weight may be clipped.\n"),
               (unsigned long long)w_count[0]);


  if (w_count[1] > 0)
    bft_printf(_("\n"
                 "%llu boundary faces have a too small distance between\n"
                 "cell center and face center.\n"),
               (unsigned long long)w_count[1]);
}

/*----------------------------------------------------------------------------
 * Compute some vectors to handle non-orthogonalities.
 *
 * Let a face and I, J the centers of neighboring cells
 *   (only I is defined for a border face)
 *
 * The face is oriented from I to J, with Nij its normal.
 *   (border faces are oriented towards the exterior)
 * The norm of Nij is 1.
 * The face surface is Sij.
 *
 * I' and J' are defined as the orthogonal projection of I and J on the line
 * orthogonal to the face passing through the center of gravity F of the face.
 *   (only I' is defined for a border face)
 *
 * We compute here the vector I'J' for interior faces (dijpf)
 *                 the vector II'  for border faces   (diipb)
 *                 the vector OF   for interior faces (dofij)
 *
 * We also have the following formulae
 *   II' = IG - (IG.Nij)Nij
 *   JJ' = JG - (JG.Nij)Nij
 *
 * parameters:
 *   dim            <--  dimension
 *   n_i_faces      <--  number of interior faces
 *   n_b_faces      <--  number of border  faces
 *   i_face_cells   <--  interior "faces -> cells" connectivity
 *   b_face_cells   <--  border "faces -> cells" connectivity
 *   i_face_u_norm  <--  unit normal of interior faces
 *   b_face_u_norm  <--  unit normal of border faces
 *   i_face_cog     <--  center of gravity of interior faces
 *   b_face_cog     <--  center of gravity of border faces
 *   cell_cen       <--  cell center
 *   weight         <--  weighting factor (Aij=pond Ai+(1-pond)Aj)
 *   dijpf          -->  vector i'j' for interior faces
 *   diipb          -->  vector ii'  for border faces
 *   dofij          -->  vector OF   for interior faces
 *----------------------------------------------------------------------------*/

static void
_compute_face_vectors(int              dim,
                      cs_lnum_t        n_i_faces,
                      cs_lnum_t        n_b_faces,
                      const cs_lnum_t  i_face_cells[][2],
                      const cs_lnum_t  b_face_cells[],
                      const cs_real_t  i_face_u_normal[][3],
                      const cs_real_t  b_face_u_normal[][3],
                      const cs_real_t  i_face_cog[],
                      const cs_real_t  b_face_cog[],
                      const cs_real_t  cell_cen[],
                      const cs_real_t  weight[],
                      const cs_real_t  b_dist[],
                      cs_real_t        dijpf[],
                      cs_real_t        diipb[],
                      cs_real_t        dofij[])
{
  /* Interior faces */

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

    const cs_lnum_t cell_id1 = i_face_cells[face_id][0];
    const cs_lnum_t cell_id2 = i_face_cells[face_id][1];

    const cs_real_t *surfn = i_face_u_normal[face_id];

    /* ---> IJ */
    const cs_real_t vecijx
      = cell_cen[cell_id2*dim]     - cell_cen[cell_id1*dim];
    const cs_real_t vecijy
      = cell_cen[cell_id2*dim + 1] - cell_cen[cell_id1*dim + 1];
    const cs_real_t vecijz
      = cell_cen[cell_id2*dim + 2] - cell_cen[cell_id1*dim + 2];

    /* ---> DIJPP = IJ.NIJ */
    const cs_real_t dipjp = vecijx*surfn[0] + vecijy*surfn[1] + vecijz*surfn[2];

    /* ---> DIJPF = (IJ.NIJ).NIJ */
    dijpf[face_id*dim]     = dipjp*surfn[0];
    dijpf[face_id*dim + 1] = dipjp*surfn[1];
    dijpf[face_id*dim + 2] = dipjp*surfn[2];

    const cs_real_t pond = weight[face_id];

    /* ---> DOFIJ = OF */
    dofij[face_id*dim]     = i_face_cog[face_id*dim]
      - (        pond *cell_cen[cell_id1*dim]
         + (1. - pond)*cell_cen[cell_id2*dim]);

    dofij[face_id*dim + 1] = i_face_cog[face_id*dim + 1]
      - (        pond *cell_cen[cell_id1*dim + 1]
         + (1. - pond)*cell_cen[cell_id2*dim + 1]);

    dofij[face_id*dim + 2] = i_face_cog[face_id*dim + 2]
      - (        pond *cell_cen[cell_id1*dim + 2]
         + (1. - pond)*cell_cen[cell_id2*dim + 2]);
  }

  /* Boundary faces */
  cs_gnum_t w_count = 0;

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

    cs_lnum_t cell_id = b_face_cells[face_id];

    const cs_real_t *normal = b_face_u_normal[face_id];

    /* ---> IF */
    cs_real_t vec_if[3] = {
      b_face_cog[face_id*dim]     - cell_cen[cell_id*dim],
      b_face_cog[face_id*dim + 1] - cell_cen[cell_id*dim + 1],
      b_face_cog[face_id*dim + 2] - cell_cen[cell_id*dim + 2]};

    /* ---> diipb = IF - (IF.NIJ)NIJ */
    cs_math_3_orthogonal_projection(normal, vec_if, &diipb[face_id*dim]);

    /* Limiter on boundary face reconstruction */
    if (cs_glob_mesh_quantities_flag & CS_FACE_RECONSTRUCTION_CLIP) {
      cs_real_t iip = cs_math_3_norm(&diipb[face_id*dim]);

      bool is_clipped = false;
      cs_real_t corri = 1.;

      if (iip > 0.5 * b_dist[face_id]) {
        is_clipped = true;
        corri = 0.5 * b_dist[face_id] / iip;
      }

      diipb[face_id*dim]    *= corri;
      diipb[face_id*dim +1] *= corri;
      diipb[face_id*dim +2] *= corri;

      if (is_clipped)
        w_count++;
    }
  }

  cs_parall_counter(&w_count, 1);

  if (w_count > 0)
    bft_printf(_("\n"
                 "%llu boundary faces have a too large reconstruction distance.\n"
                 "For these faces, reconstruction are limited.\n"),
               (unsigned long long)w_count);
}

/*----------------------------------------------------------------------------
 * Compute some vectors to handle non-orthogonalities.
 *
 * Let a face and I, J the centers of neighboring cells
 *   (only I is defined for a border face)
 *
 * The face is oriented from I to J, with Nij its normal.
 *   (border faces are oriented towards the exterior)
 * The norm of Nij is 1.
 * The face surface is Sij.
 *
 * I' and J' are defined as the orthogonal projection of I and J on the line
 * orthogonal to the face passing through the center of gravity F of the face.
 *   (only I' is defined for a border face)
 *
 * We compute here the vector II' for interior faces (diipf)
 *                 the vector JJ' for interior faces (djjpf)
 *
 * We also have the following formulae
 *   II' = IF - (IF.Nij)Nij
 *   JJ' = JF - (JF.Nij)Nij
 *
 * parameters:
 *   n_cells        <--  number of cells
 *   n_i_faces      <--  number of interior faces
 *   i_face_cells   <--  interior "faces -> cells" connectivity
 *   i_face_u_norm  <--  unit normal of interior faces
 *   i_face_norm    <--  surface normal of interior faces
 *   i_face_cog     <--  center of gravity of interior faces
 *   cell_cen       <--  cell center
 *   cell_vol       <--  cell volume
 *   dist           <--  interior distance
 *   diipf          -->  vector ii' for interior faces
 *   djjpf          -->  vector jj' for interior faces
 *----------------------------------------------------------------------------*/

static void
_compute_face_sup_vectors(cs_lnum_t          n_cells,
                          cs_lnum_t          n_i_faces,
                          const cs_lnum_2_t  i_face_cells[],
                          const cs_real_t    i_face_u_normal[][3],
                          const cs_real_t    i_face_normal[][3],
                          const cs_real_t    i_face_cog[][3],
                          const cs_real_t    cell_cen[][3],
                          const cs_real_t    cell_vol[],
                          const cs_real_t    dist[],
                          cs_real_t          diipf[][3],
                          cs_real_t          djjpf[][3])
{
  cs_gnum_t w_count = 0;

  /* Interior faces */

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

    const cs_lnum_t cell_id1 = i_face_cells[face_id][0];
    const cs_lnum_t cell_id2 = i_face_cells[face_id][1];

    /* Normalized normal */
    const cs_real_t *u_normal = i_face_u_normal[face_id];

    /* ---> IF and JF */
    cs_real_t vec_if[3] = {
      i_face_cog[face_id][0] - cell_cen[cell_id1][0],
      i_face_cog[face_id][1] - cell_cen[cell_id1][1],
      i_face_cog[face_id][2] - cell_cen[cell_id1][2]};

    cs_real_t vec_jf[3] = {
      i_face_cog[face_id][0] - cell_cen[cell_id2][0],
      i_face_cog[face_id][1] - cell_cen[cell_id2][1],
      i_face_cog[face_id][2] - cell_cen[cell_id2][2]};

    /* ---> diipf = IF - (IF.Nij)Nij */
    cs_math_3_orthogonal_projection(u_normal, vec_if, diipf[face_id]);

    /* ---> djjpf = JF - (JF.Nij)Nij */
    cs_math_3_orthogonal_projection(u_normal, vec_jf, djjpf[face_id]);

    /* Limiter on interior face reconstruction */
    if (cs_glob_mesh_quantities_flag & CS_FACE_RECONSTRUCTION_CLIP) {

      cs_real_t surfn = cs_math_3_norm(i_face_normal[face_id]);
      cs_real_t iip   = cs_math_3_norm(diipf[face_id]);

      cs_real_t corri = 1.;
      bool is_clipped = false;

      if (0.5 * dist[face_id] < iip) {
        is_clipped = true;
        corri = 0.5 * dist[face_id] / iip;
      }

      diipf[face_id][0] *= corri;
      diipf[face_id][1] *= corri;
      diipf[face_id][2] *= corri;

      iip = cs_math_3_norm(diipf[face_id]);

      corri = 1.;

      if (0.9 * cell_vol[cell_id1] < surfn * iip) {
        is_clipped = true;
        corri = 0.9 * cell_vol[cell_id1] / (surfn * iip);
      }

      diipf[face_id][0] *= corri;
      diipf[face_id][1] *= corri;
      diipf[face_id][2] *= corri;

      cs_real_t jjp = cs_math_3_norm(djjpf[face_id]);

      cs_real_t corrj = 1.;

      if (0.5 * dist[face_id] < jjp) {
        is_clipped = true;
        corrj = 0.5 * dist[face_id] / jjp;
      }

      djjpf[face_id][0] *= corrj;
      djjpf[face_id][1] *= corrj;
      djjpf[face_id][2] *= corrj;

      jjp = cs_math_3_norm(djjpf[face_id]);

      corrj = 1.;
      if (0.9 * cell_vol[cell_id2] < surfn * jjp) {
        is_clipped = true;
        corrj = 0.9 * cell_vol[cell_id2] / (surfn * jjp);
      }

      if (is_clipped && cell_id1 < n_cells)
        w_count++;

      djjpf[face_id][0] *= corrj;
      djjpf[face_id][1] *= corrj;
      djjpf[face_id][2] *= corrj;

    }
  }

  cs_parall_counter(&w_count, 1);

  if (w_count > 0)
    bft_printf
      (_("\n"
         "%llu internal faces have a too large reconstruction distance.\n"
         "For these faces, reconstruction are limited.\n"
         "\n"),
       (unsigned long long)w_count);
}

/*----------------------------------------------------------------------------
 * Evaluate boundary thickness.
 *
 * parameters:
 *   m             <-- pointer to mesh structure
 *   m_quantities  <-- pointer to mesh quantities structures.
 *   b_thickness   --> boundary thickness
 *----------------------------------------------------------------------------*/

static void
_b_thickness(const cs_mesh_t             *m,
             const cs_mesh_quantities_t  *mq,
             cs_real_t                    b_thickness[])
{
  const cs_real_3_t  *cell_cen
    = (const cs_real_3_t  *)(mq->cell_cen);
  const cs_real_3_t  *b_face_cog
    = (const cs_real_3_t  *)(mq->b_face_cog);
  const cs_real_3_t  *b_face_normal
    = (const cs_real_3_t  *)(mq->b_face_normal);
  const cs_real_t  *b_face_surf
    = (const cs_real_t *)(mq->b_face_surf);

  for (cs_lnum_t f_id = 0; f_id < m->n_b_faces; f_id++) {
    cs_lnum_t c_id = m->b_face_cells[f_id];
    b_thickness[f_id]
      = (  (b_face_cog[f_id][0] - cell_cen[c_id][0])*b_face_normal[f_id][0]
         + (b_face_cog[f_id][1] - cell_cen[c_id][1])*b_face_normal[f_id][1]
         + (b_face_cog[f_id][2] - cell_cen[c_id][2])*b_face_normal[f_id][2])
        * 2.0 / b_face_surf[f_id];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute approximate cells centers as the mean of the given face
 *         centers weighted by the associated surfaces
 *         (considering solid face for cs_glob_porous_model = 3).
 *
 *           n-1
 *           Sum  Surf(Fi) G(Fi)
 *           i=0
 *  G(C) = -----------------------
 *           n-1
 *           Sum  Surf(Fi)
 *           i=0
 *
 * \param[in]   m                     pointer to mesh structure
 * \param[in]   i_f_face_cell_normal  fluid face normal (cell to face)
 * \param[in]   i_f_face_cog          fluid face centers of gravity
 *                                    (two values for the two sides)
 * \param[in]   b_f_face_normal       fluid surface normal of border faces
 * \param[in]   b_f_face_cog          fluid center of gravity of border faces
 * \param[in]   c_w_face_normal       solid surface normal immersed in the cell
 * \param[in]   c_w_face_cog          solid center of gravity immersed in the
 *                                    cell
 * \param[out]  a_cell_cen            approximate cell centers
 */
/*----------------------------------------------------------------------------*/

static void
_mesh_quantities_cell_faces_cog_solid
  (const cs_mesh_t   *m,
   const cs_real_t    i_f_face_cell_normal[][2][3],
   const cs_real_t    i_f_face_cog[][2][3],
   const cs_real_3_t  b_f_face_normal[],
   const cs_real_3_t  b_f_face_cog[],
   const cs_real_3_t  c_w_face_normal[],
   const cs_real_3_t  c_w_face_cog[],
   cs_real_t          a_cell_cen[][3])
{
  /* Mesh connectivity */

  const cs_lnum_t n_cells = m->n_cells;
  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
  const short int *cell_i_faces_sgn = ma->cell_i_faces_sgn;
  const cs_lnum_t *cell_i_faces = ma->cell_i_faces;
  const cs_lnum_t *cell_b_faces_idx = ma->cell_b_faces_idx;
  const cs_lnum_t *cell_b_faces = ma->cell_b_faces;

  cs_real_t *cell_area;
  BFT_MALLOC(cell_area, n_cells, cs_real_t);

  /* Loop on cells
     ------------- */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    /* Initialization */
    cell_area[c_id] = 0.;
    for (cs_lnum_t i = 0; i < 3; i++)
      a_cell_cen[c_id][i] = 0.;

    /* Interior faces */
    const cs_lnum_t s_id_i = c2c_idx[c_id];
    const cs_lnum_t e_id_i = c2c_idx[c_id+1];

    /* Loop on interior faces of cell c_id */
    for (cs_lnum_t cidx = s_id_i; cidx < e_id_i; cidx++) {
      const cs_lnum_t f_id = cell_i_faces[cidx];
      const short int sign = cell_i_faces_sgn[cidx];
      const cs_lnum_t side = (sign+2)%3;

      /* For each cell sharing the internal face, we update
       * cell_cen and cell_area */

      /* Computation of face area  */
      const cs_real_t area = cs_math_3_norm(i_f_face_cell_normal[f_id][side]);
      cell_area[c_id] += area;

      for (cs_lnum_t i = 0; i < 3; i++)
        a_cell_cen[c_id][i] += i_f_face_cog[f_id][side][i]*area;

    } /* End of loop on interior faces */

    const cs_lnum_t s_id_b = cell_b_faces_idx[c_id];
    const cs_lnum_t e_id_b = cell_b_faces_idx[c_id+1];

    for (cs_lnum_t cidx = s_id_b; cidx < e_id_b; cidx++) {
      const cs_lnum_t f_id = cell_b_faces[cidx];

      cs_real_t area = cs_math_3_norm(b_f_face_normal[f_id]);

      cell_area[c_id] += area;

      /* Computation of the numerator */

      for (cs_lnum_t i = 0; i < 3; i++)
        a_cell_cen[c_id][i] += b_f_face_cog[f_id][i]*area;
    }

    /* Loop on cells: optional immersed boundary contribution */
    if (c_w_face_normal != nullptr) {
      const cs_real_t area = cs_math_3_norm(c_w_face_normal[c_id]);

      cell_area[c_id] += area;

      for (cs_lnum_t i = 0; i < 3; i++)
        a_cell_cen[c_id][i] += c_w_face_cog[c_id][i] * area;
    }

    if (cell_area[c_id] > DBL_MIN) {
      for (cs_lnum_t i = 0; i < 3; i++) {
        a_cell_cen[c_id][i] /= cell_area[c_id];
      }
    }

  } /* End of loop on cell */

  BFT_FREE(cell_area);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute fluid cell centers and volumes, when immersed boundary
 * is enabled (considering solid face for cs_glob_porous_model = 3).
 *
 * Volume V and center of gravity G of the fluid cells C is computed from the
 * approximate center of the fluid cell Ga(C), the center of gravity G(Fi) and
 * the normal N(Fi) of the faces, including the solid face when the cell is
 * cut by a solid wall.
 *
 *         n-1
 *         Sum  (G(Fi) - Ga(C)) . N(Fi)
 *         i=0
 *  V(C) = --------------------------------
 *                        3
 *
 *         n-1
 *         Sum  0.75 * G(Fi) + 0.25 G(C)
 *         i=0
 *  G(C) = -----------------------------
 *                 3 * V(C_fluid)
 *
 * \param[in]   m                     pointer to mesh structure
 * \param[in]   i_f_face_cell_normal  interior fluid face normal vector
 *                                    (two values for the two sides of the face)
 * \param[in]   i_f_face_cog          interior fluid face center of gravity
 *                                    (two values for the two sides of the face)
 * \param[in]   b_f_face_normal       fluid surface normal of border faces
 * \param[in]   b_f_face_cog          fluid center of gravity of border faces
 * \param[in]   c_w_face_normal       solid surface normal immersed in the cell
 * \param[in]   c_w_face_cog          solid center of gravity immersed in the
 *                                    cell
 * \param[out]  cell_f_cen            fluid cell centers
 * \param[out]  cell_f_vol            fluid cell volumes
 */
/*----------------------------------------------------------------------------*/

static void
_compute_fluid_solid_cell_quantities
(
  const cs_mesh_t       *m,
  const cs_real_t        i_f_face_cell_normal[][2][3],
  const cs_real_t        i_f_face_cog[][2][3],
  const cs_real_t        b_f_face_normal[][3],
  const cs_real_t        b_f_face_cog[][3],
  const cs_real_t        c_w_face_normal[][3],
  const cs_real_t        c_w_face_cog[][3],
  cs_real_3_t  *restrict cell_f_cen,
  cs_real_t    *restrict cell_f_vol
)
{
  /* Mesh connectivity */

  //const cs_lnum_t  n_b_faces = CS_MAX(m->n_b_faces, m->n_b_faces_all);
  const cs_lnum_t  n_cells = m->n_cells;
  const cs_lnum_t  n_cells_ext = m->n_cells_with_ghosts;
  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
  const short int *cell_i_faces_sgn = ma->cell_i_faces_sgn;
  const cs_lnum_t *cell_i_faces = ma->cell_i_faces;
  const cs_lnum_t *cell_b_faces_idx = ma->cell_b_faces_idx;
  const cs_lnum_t *cell_b_faces = ma->cell_b_faces;

  const cs_ibm_cog_location_t cog_location
    = cs_glob_porosity_from_scan_opt->cog_location;

  /* Checking */

  assert(cell_f_cen != nullptr);
  assert(cell_f_vol != nullptr);

  /* Compute approximate cell center using face centers */

  cs_real_3_t *a_cell_cen;
  BFT_MALLOC(a_cell_cen, n_cells_ext, cs_real_3_t);

  // TODO merge with cs_mesh_quantities_cell_faces_cog
  _mesh_quantities_cell_faces_cog_solid(m,
                                        i_f_face_cell_normal,
                                        i_f_face_cog,
                                        b_f_face_normal,
                                        b_f_face_cog,
                                        c_w_face_normal,
                                        c_w_face_cog,
                                        a_cell_cen);

  /* Compute COG from pyramids sub-volume */

  cs_real_t *_cell_f_vol;
  BFT_MALLOC(_cell_f_vol, n_cells_ext, cs_real_t);

  cs_real_3_t *_cell_f_cen;
  BFT_MALLOC(_cell_f_cen, n_cells_ext, cs_real_3_t);

  /* Loop on cells
     ------------- */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    _cell_f_vol[c_id] = 0.;
    for (cs_lnum_t i = 0; i < 3; i++)
      _cell_f_cen[c_id][i] = 0.;

    /* Interior faces */
    const cs_lnum_t s_id_i = c2c_idx[c_id];
    const cs_lnum_t e_id_i = c2c_idx[c_id+1];

    /* For each cell sharing the internal face, we update
     * _cell_f_cen and cell_area */

    for (cs_lnum_t cidx = s_id_i; cidx < e_id_i; cidx++) {
      const cs_lnum_t f_id = cell_i_faces[cidx];
      const short int sign = cell_i_faces_sgn[cidx];
      const cs_lnum_t side = (sign+2)%3;

      /* Implicit subdivision of cell into face vertices-cell-center pyramids */
      const cs_real_t pyra_vol_3
        = sign*cs_math_3_distance_dot_product(a_cell_cen[c_id],
                                              i_f_face_cog[f_id][side],
                                              i_f_face_cell_normal[f_id][side]);

      for (cs_lnum_t i = 0; i < 3; i++)
        _cell_f_cen[c_id][i] += pyra_vol_3 *(  0.75*i_f_face_cog[f_id][side][i]
                                             + 0.25*a_cell_cen[c_id][i]);
      _cell_f_vol[c_id] += pyra_vol_3;

    } /* End of loop on interior faces */

    /* Loop on boundary faces */

    const cs_lnum_t s_id_b = cell_b_faces_idx[c_id];
    const cs_lnum_t e_id_b = cell_b_faces_idx[c_id+1];

    for (cs_lnum_t cidx = s_id_b; cidx < e_id_b; cidx++) {
      const cs_lnum_t f_id = cell_b_faces[cidx];

      if (c_id > -1) {
        cs_real_t pyra_vol_3
          = cs_math_3_distance_dot_product(a_cell_cen[c_id],
                                           b_f_face_cog[f_id],
                                           b_f_face_normal[f_id]);

        for (cs_lnum_t i = 0; i < 3; i++)
          _cell_f_cen[c_id][i] += pyra_vol_3 *(  0.75*b_f_face_cog[f_id][i]
                                               + 0.25*a_cell_cen[c_id][i]);
        _cell_f_vol[c_id] += pyra_vol_3;

      }
    }

    /* Add pyramid formed with the solid face */

    if (cs_math_3_norm(c_w_face_normal[c_id]) > 0.) {

      cs_real_t pyra_vol_3
        = cs_math_3_distance_dot_product(a_cell_cen[c_id],
                                         c_w_face_cog[c_id],
                                         c_w_face_normal[c_id]);
      for (cs_lnum_t i = 0; i < 3; i++)
        _cell_f_cen[c_id][i] += pyra_vol_3 *(  0.75*c_w_face_cog[c_id][i]
                                             + 0.25*a_cell_cen[c_id][i]);
      _cell_f_vol[c_id] += pyra_vol_3;

    }

    /* Finalize the computation */

    if (_cell_f_vol[c_id] > 0.) {

      for (cs_lnum_t i = 0; i < 3; i++)
        _cell_f_cen[c_id][i] /= _cell_f_vol[c_id];

      _cell_f_vol[c_id] /= 3.0;
      cell_f_vol[c_id] = _cell_f_vol[c_id];

      if (cog_location == CS_COG_FROM_FLUID_FACES) {
        for (cs_lnum_t i = 0; i < 3; i++)
          cell_f_cen[c_id][i] = a_cell_cen[c_id][i];
      }
      else if (cog_location == CS_COG_FROM_PYRAMID) {
        for (cs_lnum_t i = 0; i < 3; i++)
          cell_f_cen[c_id][i] = _cell_f_cen[c_id][i];
      }
      else if (cog_location == CS_COG_WITHOUT_RECONSTRUCTION_FOR_IBM_PLANE) {

        /* Parallel projection to the IBM wall from pyramids COG
           to be orthogonal to the wall (to have II' = 0) */

        if (cs_math_3_norm(c_w_face_normal[c_id]) > 0) {
          cs_real_t unit_n[3];
          cs_math_3_normalize(c_w_face_normal[c_id], unit_n);

          cs_real_t dot = cs_math_3_distance_dot_product(c_w_face_cog[c_id],
                                                         _cell_f_cen[c_id],
                                                         unit_n);

          for (cs_lnum_t i = 0; i < 3; i++)
            cell_f_cen[c_id][i] = c_w_face_cog[c_id][i] + dot * unit_n[i];
        }
      }
    }

  } /* End of loop on cells */

  /* Free memory */
  BFT_FREE(a_cell_cen);
  BFT_FREE(_cell_f_vol);
  BFT_FREE(_cell_f_cen);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute mesh quantities needed for preprocessing.
 *
 * \param[in]       m   pointer to mesh structure
 * \param[in, out]  mq  pointer to mesh quantities structures.
 */
/*----------------------------------------------------------------------------*/

static void
_compute_unit_normals(const cs_mesh_t       *m,
                      cs_mesh_quantities_t  *mq)
{
  cs_lnum_t  n_i_faces = m->n_i_faces;
  cs_lnum_t  n_b_faces = CS_MAX(m->n_b_faces, m->n_b_faces_all);

  const cs_real_3_t *i_face_normal = (const cs_real_3_t *)mq->i_face_normal;
  const cs_real_3_t *b_face_normal = (const cs_real_3_t *)mq->b_face_normal;

  const cs_alloc_mode_t amode = cs_alloc_mode_read_mostly;

  /* If this is not an update, allocate members of the structure */

  if (mq->i_face_u_normal == nullptr) {
    CS_MALLOC_HD(mq->i_face_u_normal, n_i_faces, cs_nreal_3_t, amode);
    cs_mem_advise_set_read_mostly(mq->i_face_u_normal);
  }

  if (mq->b_face_u_normal == nullptr) {
    CS_MALLOC_HD(mq->b_face_u_normal, n_b_faces, cs_nreal_3_t, amode);
    cs_mem_advise_set_read_mostly(mq->b_face_u_normal);
  }

# pragma omp parallel for  if (n_i_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_i_faces; i++) {
    cs_math_3_normalize(i_face_normal[i], mq->i_face_u_normal[i]);
  }

# pragma omp parallel for  if (n_b_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_b_faces; i++) {
    cs_math_3_normalize(b_face_normal[i], mq->b_face_u_normal[i]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Post-processes the immersed boundary (ib) planes for display
 *         on paraview.
 *
 * \param[in]       n_ib_cells      ib cell number
 * \param[in]       n_glob_vtx      total vertex number
 * \param[in]       ibcell_cells    connectivity ib_cell->cells
 * \param[in]       vtx_ids         vertex ids on both sides of a IB vertex
 *                                  (v0<v1)
 * \param[in]       w_vtx_idx       ib vertex indexes
 * \param[in]       face_vertex_idx vertex indexes of the ib faces
 * \param[in]       w_vtx           ib vertex coordinates
 */
/*----------------------------------------------------------------------------*/

static void
_post_plane_ib(const cs_lnum_t      n_ib_cells,
               const cs_lnum_t      n_glob_vtx,
               const cs_lnum_t      ibcell_cells[],
               const cs_lnum_t      vtx_ids[][2],
               const cs_lnum_t      w_vtx_idx[],
               const cs_lnum_t      face_vertex_idx[],
               const cs_real_t      w_vtx[][3])
{
  cs_lnum_t *face_vertex_ids;
  BFT_MALLOC(face_vertex_ids, n_glob_vtx, cs_lnum_t);
  for (cs_lnum_t i = 0; i < n_glob_vtx; i++)
    face_vertex_ids[i] = i;

  /* Reordering the vertices to obtain a closed contour
     -------------------------------------------------- */

  cs_lnum_t n_w_vtx = 0;
  cs_coord_3_t *edge;
  BFT_MALLOC(edge, n_glob_vtx, cs_coord_3_t);

  for (cs_lnum_t c_id_ib = 0; c_id_ib < n_ib_cells; c_id_ib++) {

    const cs_lnum_t c_id = ibcell_cells[c_id_ib];
    const int n_vtx = w_vtx_idx[c_id+1]-w_vtx_idx[c_id];

    /* Edge reference is the first 2 ib vtx by default, then connect other
       ib edges by finding the corresponding vtx ids associated to the
       second ib vertex */

    cs_lnum_t k = 0;
    const cs_lnum_t s_id = w_vtx_idx[c_id];
    const cs_lnum_t e_id = w_vtx_idx[c_id+1];

    const int n_edge = 0.5*n_vtx;
    cs_real_23_t edge_l[n_edge];
    cs_lnum_t vtx_ids_l[n_edge][2][2];

    for (cs_lnum_t j = s_id; j < e_id; j=j+2) {

      /* Store edges and associated vtx_ids of the ib cell */
      for (cs_lnum_t i = 0; i < 3; i++) {
        edge_l[k][0][i] = w_vtx[j][i];
        edge_l[k][1][i] = w_vtx[j+1][i];
      }
      vtx_ids_l[k][0][0] = vtx_ids[j][0];
      vtx_ids_l[k][0][1] = vtx_ids[j][1];

      vtx_ids_l[k][1][0] = vtx_ids[j+1][0];
      vtx_ids_l[k][1][1] = vtx_ids[j+1][1];
      k++;
    }

    /* The reference edge is the first k=0 */
    for (cs_lnum_t i = 0; i < 3; i++) {
      edge[n_w_vtx][i]   = edge_l[0][0][i];
      edge[n_w_vtx+1][i] = edge_l[0][1][i];
    }

    n_w_vtx += 2;

    /* Reference edge id = 0 */
    int edge_id = 0;
    int n_iter = 0;
    int cpt = n_edge - 1;

    /* Fist edge is x0->x1. The vtx to connect is x1. */
    cs_lnum_t vtx_to_connect[2] = {vtx_ids_l[0][1][0], vtx_ids_l[0][1][1]};

    do {

      for (int j = 0; j < n_edge; j++) {

        if (j != edge_id) {

          int new_vtx_id = -1;

          /* x0 of edge j is connected to the vertex */
          if (   vtx_ids_l[j][0][0] == vtx_to_connect[0]
              && vtx_ids_l[j][0][1] == vtx_to_connect[1]) {

            new_vtx_id = 1;
          }
          /* x1 of edge j is connected to the vertex */
          else if (   vtx_ids_l[j][1][0] == vtx_to_connect[0]
                   && vtx_ids_l[j][1][1] == vtx_to_connect[1]) {

            new_vtx_id = 0;
          }
          else {
            continue;
          }

          vtx_to_connect[0] = vtx_ids_l[j][new_vtx_id][0];
          vtx_to_connect[1] = vtx_ids_l[j][new_vtx_id][1];

          for (cs_lnum_t i = 0; i < 3; i++)
            edge[n_w_vtx][i] = edge_l[j][new_vtx_id][i];

          n_w_vtx++;
          edge_id = j;
          cpt -= 1;
          break;
        }
      }

      n_iter++;
      if (n_iter > 100)
        bft_error(__FILE__, __LINE__, 0,
                  _("Problem connecting vertices in cell_id = %d"), c_id);

    } while (cpt != 1);

  } /* End loop on ib cells */

  /* Create IBM mesh and give IBM faces and vertices */

  fvm_nodal_t *ib_mesh = fvm_nodal_create("ib_mesh", 3);
  const cs_lnum_t face_list_shift[2] = {0, n_ib_cells};

  const cs_lnum_t *face_vtx_idx[1] = {face_vertex_idx};
  const cs_lnum_t *face_vtx_ids[1] = {face_vertex_ids};

  fvm_nodal_from_desc_add_faces(ib_mesh,
                                -1,
                                n_ib_cells,
                                nullptr,
                                1,
                                face_list_shift,
                                face_vtx_idx,
                                face_vtx_ids,
                                nullptr,
                                nullptr);

  fvm_nodal_set_shared_vertices(ib_mesh,
                                (const cs_coord_t *)edge);

  /* Parallel numbering */

  fvm_io_num_t *face_io_num = nullptr;
  fvm_io_num_t *vtx_io_num = nullptr;

  if (cs_glob_n_ranks > 1) {

    /* Order faces by increasing global number */

    face_io_num = fvm_io_num_create_from_scan(n_ib_cells);
    const cs_gnum_t *face_gnum = fvm_io_num_get_global_num(face_io_num);
    fvm_nodal_order_faces(ib_mesh, face_gnum);
    fvm_nodal_init_io_num(ib_mesh, face_gnum, 2);

    /* Order vertices by increasing global number */

    vtx_io_num = fvm_io_num_create_from_scan(n_glob_vtx);
    const cs_gnum_t *vertex_gnum = fvm_io_num_get_global_num(vtx_io_num);
    fvm_nodal_order_vertices(ib_mesh, vertex_gnum);
    fvm_nodal_init_io_num(ib_mesh, vertex_gnum, 0);

  }

  /* Print ib_mesh info in listing */

  //fvm_nodal_dump(ib_mesh);

  /* Create default writer */

  fvm_writer_t *writer = nullptr;
  writer = fvm_writer_init("ib_mesh",
                           "postprocessing",
                           cs_post_get_default_format(),
                           cs_post_get_default_format_options(),
                           FVM_WRITER_FIXED_MESH);

  /* Write current mesh */

  fvm_writer_export_nodal(writer, ib_mesh);

  /* Free memory */

  fvm_writer_finalize(writer);
  ib_mesh = fvm_nodal_destroy(ib_mesh);

  if (cs_glob_n_ranks > 1) {
    fvm_io_num_destroy(vtx_io_num);
    fvm_io_num_destroy(face_io_num);
  }

  BFT_FREE(face_vertex_ids);
  BFT_FREE(edge);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Query or modification of the option for computing cell centers.
 *
 * \param[in]  algo_choice  < 0 : query
 *                            0 : computation based on face centers (default)
 *                            1 : computation by cell sub-volumes
 *
 * \return  0 or 1 according to the selected algorithm
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_quantities_cell_cen_choice(int  algo_choice)
{
  if (algo_choice > -1 && algo_choice < 2)
    _cell_cen_algorithm = algo_choice;

  return _cell_cen_algorithm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Query or modification of the option for computing face centers.
 *
 * \param[in]  algo_choice  < 0 : query
 *                            0 : standard computation
 *                            1 : use adjustment for volume
 *                                from versions 1.1 to 5.3
 *
 * \return  0 or 1 according to the selected algorithm
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_quantities_face_cog_choice(int  algo_choice)
{
  if (algo_choice > -1 && algo_choice < 2)
    _ajust_face_cog_compat_v11_v52 = algo_choice;

  return _ajust_face_cog_compat_v11_v52;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a mesh quantities structure.
 *
 * \return  pointer to created cs_mesh_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

cs_mesh_quantities_t *
cs_mesh_quantities_create(void)
{
  cs_mesh_quantities_t  *mesh_quantities = nullptr;

  BFT_MALLOC(mesh_quantities, 1, cs_mesh_quantities_t);

  mesh_quantities->cell_cen = nullptr;
  mesh_quantities->cell_f_cen = nullptr;
  mesh_quantities->cell_s_cen = nullptr;
  mesh_quantities->cell_vol = nullptr;
  mesh_quantities->cell_f_vol = nullptr;
  mesh_quantities->i_face_normal = nullptr;
  mesh_quantities->b_face_normal = nullptr;
  mesh_quantities->i_f_face_normal = nullptr;
  mesh_quantities->b_f_face_normal = nullptr;
  mesh_quantities->c_w_face_normal = nullptr;
  mesh_quantities->i_face_cog = nullptr;
  mesh_quantities->b_face_cog = nullptr;
  mesh_quantities->i_f_face_cog = nullptr;
  mesh_quantities->b_f_face_cog = nullptr;
  mesh_quantities->c_w_face_cog = nullptr;
  mesh_quantities->i_face_surf = nullptr;
  mesh_quantities->b_face_surf = nullptr;
  mesh_quantities->i_f_face_surf = nullptr;
  mesh_quantities->b_f_face_surf = nullptr;
  mesh_quantities->c_w_face_surf = nullptr;
  mesh_quantities->i_face_u_normal = nullptr;
  mesh_quantities->b_face_u_normal = nullptr;
  mesh_quantities->i_f_face_factor = nullptr;
  mesh_quantities->b_f_face_factor = nullptr;
  mesh_quantities->i_dist = nullptr;
  mesh_quantities->b_dist = nullptr;
  mesh_quantities->c_w_dist_inv = nullptr;
  mesh_quantities->weight = nullptr;
  mesh_quantities->i_f_weight = nullptr;
  mesh_quantities->dijpf = nullptr;
  mesh_quantities->diipb = nullptr;
  mesh_quantities->dofij = nullptr;
  mesh_quantities->diipf = nullptr;
  mesh_quantities->djjpf = nullptr;
  mesh_quantities->corr_grad_lin_det = nullptr;
  mesh_quantities->corr_grad_lin = nullptr;
  mesh_quantities->b_sym_flag = nullptr;
  mesh_quantities->has_disable_flag = 0;
  mesh_quantities->c_disable_flag = nullptr;
  mesh_quantities->bad_cell_flag = nullptr;

  return (mesh_quantities);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a mesh quantities structure.
 *
 * \param[in]  mq  pointer to mesh quantities structures
 *
 * \return  nullptr
 */
/*----------------------------------------------------------------------------*/

cs_mesh_quantities_t *
cs_mesh_quantities_destroy(cs_mesh_quantities_t  *mq)
{
  cs_mesh_quantities_free_all(mq);

  BFT_FREE(mq);

  return (mq);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reset a mesh quantities structure to its empty initial state.
 *
 * \param[in]  mq  pointer to mesh quantities structures
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_quantities_free_all(cs_mesh_quantities_t  *mq)
{
  CS_FREE_HD(mq->cell_cen);
  mq->cell_f_cen = nullptr;
  mq->cell_s_cen = nullptr;
  BFT_FREE(mq->cell_vol);
  mq->cell_f_vol = nullptr;

  BFT_FREE(mq->i_face_normal);
  BFT_FREE(mq->b_face_normal);
  mq->i_f_face_normal = nullptr;
  mq->b_f_face_normal = nullptr;
  mq->c_w_face_normal = nullptr;

  CS_FREE_HD(mq->i_face_cog);
  CS_FREE_HD(mq->b_face_cog);
  mq->i_f_face_cog = nullptr;
  mq->b_f_face_cog = nullptr;
  mq->c_w_face_cog = nullptr;
  CS_FREE_HD(mq->i_face_surf);
  CS_FREE_HD(mq->b_face_surf);
  mq->c_w_face_surf = nullptr;

  CS_FREE_HD(mq->i_face_u_normal);
  CS_FREE_HD(mq->b_face_u_normal);

  mq->i_f_face_factor = nullptr;
  mq->b_f_face_factor = nullptr;

  CS_FREE_HD(mq->i_dist);
  CS_FREE_HD(mq->b_dist);
  mq->c_w_dist_inv = nullptr;

  CS_FREE_HD(mq->weight);
  CS_FREE_HD(mq->i_f_weight);

  CS_FREE_HD(mq->dijpf);
  CS_FREE_HD(mq->diipb);
  CS_FREE_HD(mq->dofij);
  CS_FREE_HD(mq->diipf);
  CS_FREE_HD(mq->djjpf);

  BFT_FREE(mq->corr_grad_lin_det);
  CS_FREE_HD(mq->corr_grad_lin);
  CS_FREE_HD(mq->b_sym_flag);
  CS_FREE_HD(mq->c_disable_flag);
  BFT_FREE(mq->bad_cell_flag);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute mesh quantities needed for preprocessing.
 *
 * \param[in]       m   pointer to mesh structure
 * \param[in, out]  mq  pointer to mesh quantities structures.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_quantities_compute_preprocess(const cs_mesh_t       *m,
                                      cs_mesh_quantities_t  *mq)
{
  cs_lnum_t  n_i_faces = m->n_i_faces;
  cs_lnum_t  n_b_faces = CS_MAX(m->n_b_faces, m->n_b_faces_all);
  cs_lnum_t  n_cells_with_ghosts = m->n_cells_with_ghosts;

  const cs_alloc_mode_t amode = cs_alloc_mode_read_mostly;

  /* If this is not an update, allocate members of the structure */

  if (mq->cell_cen == nullptr) {
    CS_MALLOC_HD(mq->cell_cen, n_cells_with_ghosts*3, cs_real_t, amode);
    cs_mem_advise_set_read_mostly(mq->cell_cen);
  }

  if (mq->cell_vol == nullptr) {
    CS_MALLOC_HD(mq->cell_vol, n_cells_with_ghosts, cs_real_t, amode);
    cs_mem_advise_set_read_mostly(mq->cell_vol);
  }

  if (mq->i_face_normal == nullptr) {
    CS_MALLOC_HD(mq->i_face_normal, n_i_faces*3, cs_real_t, amode);
    cs_mem_advise_set_read_mostly(mq->i_face_normal);
  }

  if (mq->b_face_normal == nullptr) {
    CS_MALLOC_HD(mq->b_face_normal, n_b_faces*3, cs_real_t, amode);
    cs_mem_advise_set_read_mostly(mq->b_face_normal);
  }

  if (mq->i_face_cog == nullptr) {
    CS_MALLOC_HD(mq->i_face_cog, n_i_faces*3, cs_real_t, amode);
    cs_mem_advise_set_read_mostly(mq->i_face_cog);
  }

  if (mq->b_face_cog == nullptr) {
    CS_MALLOC_HD(mq->b_face_cog, n_b_faces*3, cs_real_t, amode);
    cs_mem_advise_set_read_mostly(mq->b_face_cog);
  }

  if (mq->i_face_surf == nullptr) {
    CS_MALLOC_HD(mq->i_face_surf, n_i_faces, cs_real_t, amode);
    cs_mem_advise_set_read_mostly(mq->i_face_surf);
  }

  if (mq->b_face_surf == nullptr) {
    CS_MALLOC_HD(mq->b_face_surf, n_b_faces, cs_real_t, amode);
    cs_mem_advise_set_read_mostly(mq->b_face_surf);
  }

  /* Compute face centers of gravity, normals, and surfaces */

  _compute_face_quantities(n_i_faces,
                           (const cs_real_3_t *)m->vtx_coord,
                           m->i_face_vtx_idx,
                           m->i_face_vtx_lst,
                           (cs_real_3_t *)mq->i_face_cog,
                           (cs_real_3_t *)mq->i_face_normal);

  _compute_face_surface(n_i_faces,
                        mq->i_face_normal,
                        mq->i_face_surf);

  _compute_face_quantities(n_b_faces,
                           (const cs_real_3_t *)m->vtx_coord,
                           m->b_face_vtx_idx,
                           m->b_face_vtx_lst,
                           (cs_real_3_t *)mq->b_face_cog,
                           (cs_real_3_t *)mq->b_face_normal);

  _compute_face_surface(n_b_faces,
                        mq->b_face_normal,
                        mq->b_face_surf);

  if (cs_glob_mesh_quantities_flag & CS_FACE_CENTER_REFINE) {
    _refine_warped_face_centers
      (n_i_faces,
       (const cs_real_3_t *)m->vtx_coord,
       m->i_face_vtx_idx,
       m->i_face_vtx_lst,
       (cs_real_3_t *)mq->i_face_cog,
       (const cs_real_3_t *)mq->i_face_normal);

    _refine_warped_face_centers
      (n_b_faces,
       (const cs_real_3_t *)m->vtx_coord,
       m->b_face_vtx_idx,
       m->b_face_vtx_lst,
       (cs_real_3_t *)mq->b_face_cog,
       (const cs_real_3_t *)mq->b_face_normal);
  }

  if (_ajust_face_cog_compat_v11_v52) {
    _adjust_face_cog_v11_v52
      (n_i_faces,
       (const cs_real_3_t *)m->vtx_coord,
       m->i_face_vtx_idx,
       m->i_face_vtx_lst,
       (cs_real_3_t *)mq->i_face_cog,
       (const cs_real_3_t *)mq->i_face_normal);

    _adjust_face_cog_v11_v52
      (n_b_faces,
       (const cs_real_3_t *)m->vtx_coord,
       m->b_face_vtx_idx,
       m->b_face_vtx_lst,
       (cs_real_3_t *)mq->b_face_cog,
       (const cs_real_3_t *)mq->b_face_normal);
  }

  /* Compute cell centers from face barycenters or vertices */

  bool volume_computed = false;

  switch (_cell_cen_algorithm) {

  case 0:
    cs_mesh_quantities_cell_faces_cog(m,
                                      mq->i_face_normal,
                                      mq->i_face_cog,
                                      mq->b_face_normal,
                                      mq->b_face_cog,
                                      mq->cell_cen);

    break;
  case 1:
    _compute_cell_quantities(m,
                             (const cs_real_3_t *)mq->i_face_normal,
                             (const cs_real_3_t *)mq->i_face_cog,
                             (const cs_real_3_t *)mq->b_face_normal,
                             (const cs_real_3_t *)mq->b_face_cog,
                             (cs_real_3_t *)mq->cell_cen,
                             mq->cell_vol);
    volume_computed = true;
    break;

  default:
    assert(0);

  }

  if (cs_glob_mesh_quantities_flag & CS_CELL_CENTER_CORRECTION) {
    _recompute_cell_cen_face(m,
                             (const cs_real_3_t *)(mq->i_face_normal),
                             (const cs_real_3_t *)(mq->i_face_cog),
                             (const cs_real_3_t *)(mq->b_face_normal),
                             (const cs_real_3_t *)(mq->b_face_cog),
                             (cs_real_3_t *)(mq->cell_cen));
    volume_computed = false; /* should not be different with plane faces,
                                not sure with warped faces */
  }

  /* Recompute face centers as the middle of two cell centers if possible */
  if (cs_glob_mesh_quantities_flag & CS_CELL_FACE_CENTER_CORRECTION) {

     if (m->halo != nullptr) {
       cs_halo_sync_var_strided(m->halo, CS_HALO_EXTENDED,
                                mq->cell_cen, 3);
       if (m->n_init_perio > 0)
         cs_halo_perio_sync_coords(m->halo, CS_HALO_EXTENDED,
                                   mq->cell_cen);
     }

    _correct_cell_face_center(m,
                              n_cells_with_ghosts,
                              n_i_faces,
                              n_b_faces,
                              (const cs_lnum_2_t *)(m->i_face_cells),
                              m->b_face_cells,
                              (cs_real_3_t *)(mq->cell_cen),
                              (cs_real_3_t *)(mq->i_face_cog),
                              (cs_real_3_t *)(mq->b_face_cog),
                              (cs_real_3_t *)(mq->i_face_normal),
                              (cs_real_3_t *)(mq->b_face_normal));

    volume_computed = false;

  }

  /* Compute the volume of cells */

  if (volume_computed == false)
    _compute_cell_volume(m,
                         (const cs_real_3_t *)(mq->i_face_normal),
                         (const cs_real_3_t *)(mq->i_face_cog),
                         (const cs_real_3_t *)(mq->b_face_normal),
                         (const cs_real_3_t *)(mq->b_face_cog),
                         (const cs_real_3_t *)(mq->cell_cen),
                         mq->cell_vol);

  /* Correction of small or negative volumes
     (doesn't conserve the total volume) */

  if (cs_glob_mesh_quantities_flag & CS_CELL_VOLUME_RATIO_CORRECTION)
    _cell_bad_volume_correction(m,
                                mq->cell_vol);

  /* Compute unit normals */

  _compute_unit_normals(m, mq);

  /* Synchronize geometric quantities */

  if (m->halo != nullptr) {

    cs_halo_sync_var_strided(m->halo, CS_HALO_EXTENDED,
                             mq->cell_cen, 3);
    if (m->n_init_perio > 0)
      cs_halo_perio_sync_coords(m->halo, CS_HALO_EXTENDED,
                                mq->cell_cen);

    cs_halo_sync_var(m->halo, CS_HALO_EXTENDED, mq->cell_vol);

  }

  _cell_volume_reductions(m,
                          mq->cell_vol,
                          &(mq->min_vol),
                          &(mq->max_vol),
                          &(mq->tot_vol));

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {

    cs_real_t  _min_vol, _max_vol, _tot_vol;

    MPI_Allreduce(&(mq->min_vol), &_min_vol, 1, CS_MPI_REAL,
                  MPI_MIN, cs_glob_mpi_comm);

    MPI_Allreduce(&(mq->max_vol), &_max_vol, 1, CS_MPI_REAL,
                  MPI_MAX, cs_glob_mpi_comm);

    MPI_Allreduce(&(mq->tot_vol), &_tot_vol, 1, CS_MPI_REAL,
                  MPI_SUM, cs_glob_mpi_comm);

    mq->min_vol = _min_vol;
    mq->max_vol = _max_vol;
    mq->tot_vol = _tot_vol;

  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell and faces quantities needed at the immersed boundaries.
 *
 * \param[in]       m               pointer to mesh structure
 * \param[in]       cen_points      point belonging to the immersed solid plane
 *                                  for each cell (or nullptr)
 * \param[in, out]  mq              pointer to mesh quantities structure
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_quantities_solid_compute(const cs_mesh_t       *m,
                                 const cs_real_3_t     *cen_points,
                                 cs_mesh_quantities_t  *mq)
{
  const cs_lnum_t * i_face_vtx_idx
    = (const cs_lnum_t *)m->i_face_vtx_idx;
  const cs_lnum_t * b_face_vtx_idx
    = (const cs_lnum_t *)m->b_face_vtx_idx;
  const cs_lnum_t * b_face_vtx_lst
    = (const cs_lnum_t *)m->b_face_vtx_lst;
  const cs_real_3_t *restrict vtx_coord
    = (const cs_real_3_t *)m->vtx_coord;
  const cs_real_3_t *restrict i_face_cog = (const cs_real_3_t *)mq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog = (const cs_real_3_t *)mq->b_face_cog;

  cs_real_3_t *restrict i_f_face_cog    = (cs_real_3_t *)mq->i_f_face_cog;
  cs_real_3_t *restrict b_f_face_cog    = (cs_real_3_t *)mq->b_f_face_cog;
  cs_real_3_t *restrict i_f_face_normal = (cs_real_3_t *)mq->i_f_face_normal;
  cs_real_3_t *restrict b_f_face_normal = (cs_real_3_t *)mq->b_f_face_normal;
  cs_real_3_t *restrict c_w_face_normal = (cs_real_3_t *)mq->c_w_face_normal;
  cs_real_t *restrict c_w_face_surf     = (cs_real_t *)mq->c_w_face_surf;
  cs_real_3_t *restrict c_w_face_cog    = (cs_real_3_t *)mq->c_w_face_cog;
  cs_field_t *f_poro = cs_field_by_name("porosity");

  const cs_lnum_t n_cells = m->n_cells;
  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  cs_mesh_adjacencies_update_cell_i_faces();

  const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
  const cs_lnum_t *cell_i_faces = ma->cell_i_faces;
  const cs_lnum_t *cell_b_faces_idx = ma->cell_b_faces_idx;
  const cs_lnum_t *cell_b_faces = ma->cell_b_faces;
  const short int *cell_i_faces_sgn = ma->cell_i_faces_sgn;

  /* Initialization */
  _compute_cell_quantities(m,
                          (const cs_real_3_t *)mq->i_face_normal,
                          (const cs_real_3_t *)mq->i_face_cog,
                          (const cs_real_3_t *)mq->b_face_normal,
                          (const cs_real_3_t *)mq->b_face_cog,
                          (cs_real_3_t *)mq->cell_f_cen,
                           mq->cell_vol);

  /* If no points belonging to the plane are given, stop here */
  if (cen_points == nullptr) {

    if (m->halo != nullptr) {

      cs_halo_sync_var_strided(m->halo, CS_HALO_EXTENDED,
                               mq->cell_f_cen, 3);
      if (m->n_init_perio > 0)
        cs_halo_perio_sync_coords(m->halo, CS_HALO_EXTENDED,
                                  mq->cell_f_cen);

      cs_halo_sync_var(m->halo, CS_HALO_EXTENDED, mq->cell_vol);

    }

    cs_array_real_copy(3*m->n_b_faces, mq->b_face_cog, mq->b_f_face_cog);
    return;

  }

  cs_real_23_t *i_f_face_cog_dual;
  BFT_MALLOC(i_f_face_cog_dual, m->n_i_faces, cs_real_23_t);
  cs_array_real_set_scalar(3*2*m->n_i_faces,
                           -HUGE_VAL,
                           (cs_real_t *)i_f_face_cog_dual);

  cs_real_23_t *i_f_face_cell_normal;
  BFT_MALLOC(i_f_face_cell_normal, m->n_i_faces, cs_real_23_t);

  cs_array_real_set_scalar(3*2*m->n_i_faces,
                           -HUGE_VAL,
                           (cs_real_t *)i_f_face_cell_normal);

  cs_array_real_set_scalar(3*m->n_cells_with_ghosts,
                           0.,
                           (cs_real_t *)c_w_face_cog);

  cs_array_real_set_scalar(m->n_cells_with_ghosts,
                           0.,
                           (cs_real_t *)c_w_face_surf);

  /* ib face vertex number */
  cs_lnum_t *w_vtx_idx;
  BFT_MALLOC(w_vtx_idx, m->n_cells+1, cs_lnum_t);
  w_vtx_idx[0] = 0;

  /* ib face vertex coordinates */
  cs_real_3_t *w_vtx;
  cs_lnum_t w_vtx_max_size = m->n_cells*16;
  BFT_MALLOC(w_vtx, w_vtx_max_size, cs_real_3_t);

  /* two vertices id for an ib edge */
  cs_lnum_2_t *vtx_ids;
  BFT_MALLOC(vtx_ids, w_vtx_max_size, cs_lnum_2_t);

  /* compute and store the first vertex on the solid immersed face
   * to compute its surface and COG */

  cs_real_3_t *v_w_ref;
  BFT_MALLOC(v_w_ref, m->n_cells_with_ghosts, cs_real_3_t);
  cs_array_real_set_scalar(3*m->n_cells_with_ghosts,
                           0.,
                           (cs_real_t *)v_w_ref);

  cs_lnum_t *flag_id;
  BFT_MALLOC(flag_id, m->n_cells_with_ghosts, cs_lnum_t);
  cs_array_lnum_fill_zero(m->n_cells_with_ghosts, flag_id);

  cs_real_t *sum_surf;
  BFT_MALLOC(sum_surf, m->n_cells_with_ghosts, cs_real_t);
  cs_array_real_fill_zero(m->n_cells_with_ghosts, sum_surf);

  cs_real_t *_i_f_surf;
  BFT_MALLOC(_i_f_surf, m->n_i_faces, cs_real_t);
  cs_array_real_set_scalar(m->n_i_faces, DBL_MAX, _i_f_surf);

  cs_lnum_t w_vtx_s_id = 0;

  /* ib cell number */
  cs_lnum_t n_ib_cells = 0;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    cs_lnum_t n_w_vtx = 0;

    /* First guess of fluid volume */
    mq->cell_f_vol[c_id] = (1 - mq->c_disable_flag[c_id]) * mq->cell_vol[c_id];

    /* Interior faces */
    const cs_lnum_t s_id_i = c2c_idx[c_id];
    const cs_lnum_t e_id_i = c2c_idx[c_id+1];

    /* Loop on interior faces of cell c_id */
    for (cs_lnum_t cidx = s_id_i; cidx < e_id_i; cidx++) {
      const cs_lnum_t face_id = cell_i_faces[cidx];
      const cs_lnum_t side = (cell_i_faces_sgn[cidx]+2)%3;

      /* Initialization */

      for (cs_lnum_t i = 0; i < 3; i++) {
        i_f_face_cell_normal[face_id][side][i] = 0.0;
        i_f_face_cog_dual[face_id][side][i] = 0.0;
      }

      const cs_lnum_t s_id = i_face_vtx_idx[face_id];
      const cs_lnum_t e_id = i_face_vtx_idx[face_id + 1];
      const cs_lnum_t n_face_vertices = e_id - s_id;

      const cs_lnum_t *vertex_ids = m->i_face_vtx_lst + s_id;
      /* Number of solid vertices seen from each side of the inner face */
      cs_lnum_t n_s_face_vertices = 0;
      cs_lnum_t f_vtx_id = 0;

      cs_real_t *_vc[10][3];
      cs_real_3_t *vc = (cs_real_3_t *)_vc;
      if (n_face_vertices > 10)
        BFT_MALLOC(vc, n_face_vertices, cs_real_3_t);

      cs_array_real_set_scalar(3*n_face_vertices,
                               -1.,
                               (cs_real_t *)vc);

      const cs_lnum_t max_n_f_face_vertices = 2*n_face_vertices;

      /* Build fluid vertices coordinates array */
      cs_real_t *_f_vtx_coord[20][3];
      cs_real_3_t *f_vtx_coord = (cs_real_3_t *)_f_vtx_coord;
      if (max_n_f_face_vertices > 20)
        BFT_MALLOC(f_vtx_coord, max_n_f_face_vertices, cs_real_3_t);

      cs_array_real_set_scalar(3*max_n_f_face_vertices,
                               -1.,
                               (cs_real_t *)f_vtx_coord);

      /* Loop over vertices of the face */
      for (cs_lnum_t j = 0; j < n_face_vertices; j++) {

        for (cs_lnum_t i = 0; i < 3; i++)
          vc[j][i] = vtx_coord[vertex_ids[j]][i] - cen_points[c_id][i];

        /* Evaluating if the face vertex is solid or fluid */
        /* Vertex seen from cell c_id is solid */
        cs_real_t vn[3], nw[3];
        cs_math_3_normalize(vc[j], vn);
        cs_math_3_normalize(c_w_face_normal[c_id], nw);

        if (cs_math_3_dot_product(vn, nw) > 0.)
          n_s_face_vertices += 1;

      }
      /* If the cell is completely immersed in the solid,
       * all the faces are solid.
       * The c_w_face_normal has a zero norm because it is not at the
       * immersed interface.
       */
      if (   cs_math_3_norm(c_w_face_normal[c_id]) < DBL_MIN
          && mq->cell_f_vol[c_id] < DBL_MIN)
        n_s_face_vertices = n_face_vertices;

      if (n_s_face_vertices < n_face_vertices) {
        for (cs_lnum_t i = 0; i < 3; i++)
          i_f_face_cog_dual[face_id][side][i] = i_face_cog[face_id][i];
      }

      /* Fluid face from one side */
      if (n_s_face_vertices == 0) {
        for (cs_lnum_t i = 0; i < 3; i++)
          i_f_face_cell_normal[face_id][side][i]
            = mq->i_face_normal[face_id*3+i];
      }

      /* We deal with a cell at the immersed interface */
      if (n_s_face_vertices > 0 && n_s_face_vertices < n_face_vertices) {

        /* Counter for the number of vertices at the wall immersed interface */
        cs_lnum_t s_vtx_id = -1;

        /* Store the (two) vertices at the interface */
        cs_real_t *_v_w_inteface[2][3];
        cs_real_3_t *v_w_inteface = (cs_real_3_t *)_v_w_inteface;
        cs_array_real_fill_zero(3 * 2, (cs_real_t *)v_w_inteface);

        /* 3 consecutive vertices */
        cs_lnum_3_t loc_id = {n_face_vertices - 1, 0, 1};

        /* Loop over face vertices */
        for (cs_lnum_t j = 0; j < n_face_vertices; loc_id[0] = loc_id[1], j++) {

          loc_id[1] = j;
          loc_id[2] = (j+1)%n_face_vertices;

          /* Store if the three consecutive vertices are solid or not */
          cs_real_t vc_dot_nw[3];

          for (cs_lnum_t i = 0; i < 3; i++) {
            cs_real_3_t vn, nw;
            cs_math_3_normalize(vc[loc_id[i]], vn);
            cs_math_3_normalize(c_w_face_normal[c_id], nw);

            /* Strictly positive if solid */
            vc_dot_nw[i] = cs_math_3_dot_product(vn, nw);
          }

          /* Vertex 1 is solid
           * Projecting solid vertex to the wall-plane along an edge */
          if (vc_dot_nw[1] > 0.) {

            cs_real_t vfluid[3];

            /* 1st neighbor fluid vertex */
            if (!(vc_dot_nw[0] > 0.)) {
              s_vtx_id++;
              if (s_vtx_id >= 2)
                BFT_REALLOC(v_w_inteface, s_vtx_id + 1, cs_real_3_t);

              _proj_solid_vtx_to_plane(vc[loc_id[0]], /* direction of proj */
                                       vc[loc_id[1]], /* projected vtx */
                                       c_w_face_normal[c_id],
                                       cen_points[c_id],
                                       vfluid); /* relative to origin */

              for (cs_lnum_t i = 0; i < 3; i++) {
                f_vtx_coord[f_vtx_id][i] = vfluid[i];
                /* Storing vertex for solid face computations */
                v_w_inteface[s_vtx_id][i] = vfluid[i];
              }
              if (flag_id[c_id] == 0) {
                /* Storing flag for solid face computations */
                for (cs_lnum_t i = 0; i < 3; i++)
                  v_w_ref[c_id][i] = vfluid[i];
              }
              f_vtx_id++;
              flag_id[c_id] = 1;

              if (w_vtx_s_id + n_w_vtx > w_vtx_max_size) {
                w_vtx_max_size *= 2;
                BFT_REALLOC(w_vtx, w_vtx_max_size, cs_real_3_t);
                BFT_REALLOC(vtx_ids, w_vtx_max_size, cs_lnum_2_t);
              }

              /* Coordinates of the ib vertex */
              for (cs_lnum_t i = 0; i < 3; i++)
                w_vtx[w_vtx_s_id + n_w_vtx][i] = vfluid[i];

              /* Vertex ids of the ib vertex */
              if (vertex_ids[loc_id[0]] < vertex_ids[loc_id[1]]) {
                vtx_ids[w_vtx_s_id + n_w_vtx][0] = vertex_ids[loc_id[0]];
                vtx_ids[w_vtx_s_id + n_w_vtx][1] = vertex_ids[loc_id[1]];
              }
              else {
                vtx_ids[w_vtx_s_id + n_w_vtx][0] = vertex_ids[loc_id[1]];
                vtx_ids[w_vtx_s_id + n_w_vtx][1] = vertex_ids[loc_id[0]];
              }

              n_w_vtx++;

            }
            /* 2nd neighbor fluid vertex */
            if (!(vc_dot_nw[2] > 0.)) {

              s_vtx_id++;
              if (s_vtx_id >= 2)
                BFT_REALLOC(v_w_inteface, s_vtx_id + 1, cs_real_3_t);

              _proj_solid_vtx_to_plane(vc[loc_id[2]],
                                       vc[loc_id[1]],
                                       c_w_face_normal[c_id],
                                       cen_points[c_id],
                                       vfluid);

              for (cs_lnum_t i = 0; i < 3; i++) {
                f_vtx_coord[f_vtx_id][i] = vfluid[i];
                /* Storing vertex for solid face computations */
                v_w_inteface[s_vtx_id][i] = vfluid[i];

              }
              if (flag_id[c_id] == 0) {
                /* Storing flag for solid face computations */
                for (cs_lnum_t i = 0; i < 3; i++)
                  v_w_ref[c_id][i] = vfluid[i];
              }
              f_vtx_id++;
              flag_id[c_id] = 1;

              if (w_vtx_s_id + n_w_vtx > w_vtx_max_size) {
                w_vtx_max_size *= 2;
                BFT_REALLOC(w_vtx, w_vtx_max_size, cs_real_3_t);
                BFT_REALLOC(vtx_ids, w_vtx_max_size, cs_lnum_2_t);
              }

              /* Coordinates of the ib vertex */
              for (cs_lnum_t i = 0; i < 3; i++)
                w_vtx[w_vtx_s_id + n_w_vtx][i] = vfluid[i];

              /* Vertex ids of the ib vertex */
              if (vertex_ids[loc_id[1]] < vertex_ids[loc_id[2]]) {
                vtx_ids[w_vtx_s_id + n_w_vtx][0] = vertex_ids[loc_id[1]];
                vtx_ids[w_vtx_s_id + n_w_vtx][1] = vertex_ids[loc_id[2]];
              }
              else {
                vtx_ids[w_vtx_s_id + n_w_vtx][0] = vertex_ids[loc_id[2]];
                vtx_ids[w_vtx_s_id + n_w_vtx][1] = vertex_ids[loc_id[1]];
              }

              n_w_vtx++;

            }
          }
          /* Vertex 1 is either fully fluid or exactly on the wall */
          else {
            for (cs_lnum_t i = 0; i < 3; i++)
              f_vtx_coord[f_vtx_id][i] = vc[loc_id[1]][i] + cen_points[c_id][i];
            f_vtx_id++;

          }

        } /* Vertices loop */

        /* COG and surface of the solid face */
        cs_real_t v01[3], v02[3], vn[3];
        for (cs_lnum_t i = 0; i < 3; i++) {
          v01[i] = v_w_inteface[0][i] - v_w_ref[c_id][i];
          v02[i] = v_w_inteface[1][i] - v_w_ref[c_id][i];
        }
        cs_math_3_cross_product(v01, v02, vn);
        sum_surf[c_id] += 0.5*cs_math_3_norm(vn);

        for (cs_lnum_t i = 0; i < 3; i++)
          c_w_face_cog[c_id][i] += (  v_w_ref[c_id][i]
                                    + v_w_inteface[0][i]
                                    + v_w_inteface[1][i])
                                    * cs_math_3_norm(vn);

        if (v_w_inteface != (cs_real_3_t *)_v_w_inteface)
          BFT_FREE(v_w_inteface);
      }

      cs_lnum_t n_f_face_vertices = f_vtx_id;

      /* Computing quantities of face intersected by different planes */
      if (n_f_face_vertices > 2) {

        cs_lnum_t *_f_face_pos[8];
        cs_lnum_t *f_face_pos = (cs_lnum_t *)_f_face_pos;
        if (n_f_face_vertices > 8)
          BFT_MALLOC(f_face_pos, n_f_face_vertices, cs_lnum_t);

        for (cs_lnum_t i = 0; i < n_f_face_vertices; i++)
          f_face_pos[i] = i;

        cs_lnum_t f_vtx_idx[2] = {0, n_f_face_vertices};
        cs_real_t _i_f_face_cog[1][3]  = { {0., 0., 0.} };
        cs_real_t _i_f_face_normal[1][3] = { {0., 0., 0.} };

        cs_real_t *_f_vtx_coord_l[10][3];
        cs_real_3_t *f_vtx_coord_l = (cs_real_3_t *)_f_vtx_coord_l;
        if (n_f_face_vertices > 10)
          BFT_MALLOC(f_vtx_coord_l, n_f_face_vertices, cs_real_3_t);

        for (cs_lnum_t ffv = 0; ffv < n_f_face_vertices; ffv++) {
          for (cs_lnum_t i = 0; i < 3; i++)
            f_vtx_coord_l[ffv][i] = f_vtx_coord[ffv][i];
        }

        _compute_face_quantities(1,
                                 f_vtx_coord_l,
                                 f_vtx_idx,
                                 f_face_pos,
                                 _i_f_face_cog,
                                 _i_f_face_normal);

        if (f_vtx_coord_l != (cs_real_3_t *)_f_vtx_coord_l)
          BFT_FREE(f_vtx_coord_l);

        /* Storing fluid face COG associated to each cell */
        for (cs_lnum_t i = 0; i < 3; i++) {
          i_f_face_cog_dual[face_id][side][i] = _i_f_face_cog[0][i];
          i_f_face_cell_normal[face_id][side][i] = _i_f_face_normal[0][i];
        }

        if (f_face_pos != (cs_lnum_t *)_f_face_pos)
          BFT_FREE(f_face_pos);
      }

      if (vc != (cs_real_3_t *)_vc)
        BFT_FREE(vc);

      if (f_vtx_coord != (cs_real_3_t *)_f_vtx_coord)
        BFT_FREE(f_vtx_coord);

    } /* End loop on adjacents cells */

    const cs_lnum_t s_id_b = cell_b_faces_idx[c_id];
    const cs_lnum_t e_id_b = cell_b_faces_idx[c_id+1];

    /* Loop on boundary faces of cell c_id
       ----------------------------------- */

    for (cs_lnum_t cidx = s_id_b; cidx < e_id_b; cidx++) {
      const cs_lnum_t face_id = cell_b_faces[cidx];

      /* Initialization*/
      cs_lnum_t s_id =  b_face_vtx_idx[face_id];
      cs_lnum_t e_id =  b_face_vtx_idx[face_id + 1];
      cs_lnum_t n_face_vertices = e_id - s_id;
      cs_lnum_t n_s_face_vertices = 0;
      cs_lnum_t f_vtx_count = 0;

      cs_real_t *_vc[10][3];
      cs_real_3_t *vc = (cs_real_3_t *)_vc;
      if (n_face_vertices > 10)
        BFT_MALLOC(vc, n_face_vertices, cs_real_3_t);

      cs_array_real_set_scalar(3*n_face_vertices,
                               -1.,
                               (cs_real_t *)vc);

      const cs_lnum_t max_n_f_face_vertices = 2*n_face_vertices;

      /* Build fluid vertices coordinates array */
      cs_real_t *_f_vtx_coord[20][3];
      cs_real_3_t *f_vtx_coord = (cs_real_3_t *)_f_vtx_coord;
      if (max_n_f_face_vertices > 20)
        BFT_MALLOC(f_vtx_coord, max_n_f_face_vertices, cs_real_3_t);

      cs_array_real_set_scalar(3*max_n_f_face_vertices,
                               -1.,
                               (cs_real_t *)f_vtx_coord);

      /* Loop over vertices of the boundary face */
      for (cs_lnum_t i_v = s_id; i_v < e_id; i_v++)  {
        const cs_lnum_t j = i_v - s_id;

        /* Create vectors from the points COG to the face vertices */
        for (cs_lnum_t i = 0; i < 3; i++) {
          vc[j][i] = vtx_coord[b_face_vtx_lst[i_v]][i] - cen_points[c_id][i];
        }

        /* Evaluate if the vertex is solid or fluid */
        cs_real_t vc_dot_nw = cs_math_3_dot_product(vc[j], c_w_face_normal[c_id]);

        if (   vc_dot_nw > 0.
            && cs_math_3_norm(c_w_face_normal[c_id]) > 0.)
          n_s_face_vertices++;
      }

      /* If the cell is completely immersed in the solid,
       * all the faces are solid.
       * The c_w_face_normal has a zero norm because it is not at the
       * immersed interface. */

      if (   cs_math_3_norm(c_w_face_normal[c_id]) < DBL_MIN
          && mq->cell_f_vol[c_id] < DBL_MIN)
        n_s_face_vertices = n_face_vertices;

      for (cs_lnum_t i = 0; i < 3; i++)
        b_f_face_cog[face_id][i] = b_face_cog[face_id][i];

      /* Fully fluid face */
      if (n_s_face_vertices == 0) {
        for (cs_lnum_t i = 0; i < 3; i++) {
          b_f_face_normal[face_id][i] = mq->b_face_normal[face_id*3+i];
        }
        mq->b_f_face_surf[face_id] = mq->b_face_surf[face_id];
      }

      cs_real_t vfluid[3];

      /* Initialization */

      if (   cs_math_3_norm(c_w_face_normal[c_id]) > 0.
          && n_s_face_vertices > 0
          && n_s_face_vertices < n_face_vertices) {

        /* Counter for the number of vertices at the wall immersed interface */
        cs_lnum_t s_vtx_count = -1;

        /* Store the (two) vertices at the interface */
        cs_real_t *_v_w_inteface[2][3];
        cs_real_3_t *v_w_inteface = (cs_real_3_t *)_v_w_inteface;
        cs_array_real_fill_zero(3 * 2, (cs_real_t *)v_w_inteface);

        /* 3 consecutive vertices */
        cs_lnum_t loc_id[3] = {n_face_vertices - 1, 0, 1};

        /* Loop over face vertices */
        for (cs_lnum_t j = 0; j < n_face_vertices; loc_id[0] = loc_id[1], j++) {

          loc_id[1] = j;
          loc_id[2] = (j+1)%n_face_vertices;

          cs_real_t vc_dot_nw[3] = {cs_math_3_dot_product(vc[loc_id[0]],
                                                          c_w_face_normal[c_id]),
                                    cs_math_3_dot_product(vc[loc_id[1]],
                                                          c_w_face_normal[c_id]),
                                    cs_math_3_dot_product(vc[loc_id[2]],
                                                          c_w_face_normal[c_id])};

          /* Project solid vertex to the wall-plane */
          if (vc_dot_nw[1] > 0.) {
            /* 1st neighbor fluid vertex */
            if (!(vc_dot_nw[0] > 0.)) {
              s_vtx_count++;
              if (s_vtx_count >= 2)
                BFT_REALLOC(v_w_inteface, s_vtx_count + 1, cs_real_3_t);

              _proj_solid_vtx_to_plane(vc[loc_id[0]],
                                       vc[loc_id[1]],
                                       c_w_face_normal[c_id],
                                       cen_points[c_id],
                                       vfluid);

              for (cs_lnum_t i = 0; i < 3; i++) {
                f_vtx_coord[f_vtx_count][i] = vfluid[i];
                /* Storing vertex for solid face computations */
                v_w_inteface[s_vtx_count][i] = vfluid[i];
              }
              if (flag_id[c_id] == 0) {
                /* Storing flag for solid face computations */
                for (cs_lnum_t i = 0; i < 3; i++)
                  v_w_ref[c_id][i] = vfluid[i];
                flag_id[c_id] = 1;
              }
              f_vtx_count++;

              if (w_vtx_s_id + n_w_vtx > w_vtx_max_size) {
                w_vtx_max_size *= 2;
                BFT_REALLOC(w_vtx, w_vtx_max_size, cs_real_3_t);
                BFT_REALLOC(vtx_ids, w_vtx_max_size, cs_lnum_2_t);
              }

              /* Coordinates of the ib vertex */
              for (cs_lnum_t i = 0; i < 3; i++)
                w_vtx[w_vtx_s_id + n_w_vtx][i] = vfluid[i];

              /* Vertex ids of the ib vertex */
              if (  b_face_vtx_lst[s_id+loc_id[0]]
                  < b_face_vtx_lst[s_id+loc_id[1]]) {
                vtx_ids[w_vtx_s_id + n_w_vtx][0] = b_face_vtx_lst[s_id+loc_id[0]];
                vtx_ids[w_vtx_s_id + n_w_vtx][1] = b_face_vtx_lst[s_id+loc_id[1]];
              }
              else {
                vtx_ids[w_vtx_s_id + n_w_vtx][0] = b_face_vtx_lst[s_id+loc_id[1]];
                vtx_ids[w_vtx_s_id + n_w_vtx][1] = b_face_vtx_lst[s_id+loc_id[0]];
              }

              n_w_vtx++;

            }
            /* 2nd neighbor fluid vertex */
            if (!(vc_dot_nw[2] > 0.)) {
              s_vtx_count++;
              if (s_vtx_count >= 2)
                BFT_REALLOC(v_w_inteface, s_vtx_count + 1, cs_real_3_t);

              _proj_solid_vtx_to_plane(vc[loc_id[2]],
                                       vc[loc_id[1]],
                                       c_w_face_normal[c_id],
                                       cen_points[c_id],
                                       vfluid);

              for (cs_lnum_t i = 0; i < 3; i++) {
                f_vtx_coord[f_vtx_count][i] = vfluid[i];
                /* Store vertex for solid face computations */
                v_w_inteface[s_vtx_count][i] = vfluid[i];
              }
              if (flag_id[c_id] == 0) {
                /* Store flag for solid face computations */
                for (cs_lnum_t i = 0; i < 3; i++)
                  v_w_ref[c_id][i] = vfluid[i];
              }
              f_vtx_count++;
              flag_id[c_id] = 1;

              if (w_vtx_s_id + n_w_vtx > w_vtx_max_size) {
                w_vtx_max_size *= 2;
                BFT_REALLOC(w_vtx, w_vtx_max_size, cs_real_3_t);
                BFT_REALLOC(vtx_ids, w_vtx_max_size, cs_lnum_2_t);
              }

              /* Coordinates of the ib vertex */
              for (cs_lnum_t i = 0; i < 3; i++)
                w_vtx[w_vtx_s_id + n_w_vtx][i] = vfluid[i];

              /* Vertex ids of the ib vertex */
              if (   b_face_vtx_lst[s_id+loc_id[1]]
                  < b_face_vtx_lst[s_id+loc_id[2]]) {
                vtx_ids[w_vtx_s_id + n_w_vtx][0] = b_face_vtx_lst[s_id+loc_id[1]];
                vtx_ids[w_vtx_s_id + n_w_vtx][1] = b_face_vtx_lst[s_id+loc_id[2]];
              }
              else {
                vtx_ids[w_vtx_s_id + n_w_vtx][0] = b_face_vtx_lst[s_id+loc_id[2]];
                vtx_ids[w_vtx_s_id + n_w_vtx][1] = b_face_vtx_lst[s_id+loc_id[1]];
              }

              n_w_vtx++;

            }
          }
          else { /* Fluid or wall vertex */

            for (cs_lnum_t i = 0; i < 3; i++)
              f_vtx_coord[f_vtx_count][i] = vc[loc_id[1]][i] + cen_points[c_id][i];

            f_vtx_count++;

          }
        } /* Vertices loop */

        /* COG and normal/surf of the solid face */
        cs_real_t v01[3], v02[3], vn[3];
        for (cs_lnum_t i = 0; i < 3; i++) {
          v01[i] = v_w_inteface[0][i] - v_w_ref[c_id][i];
          v02[i] = v_w_inteface[1][i] - v_w_ref[c_id][i];
        }
        cs_math_3_cross_product(v01, v02, vn);
        sum_surf[c_id] += 0.5*cs_math_3_norm(vn);

        for (cs_lnum_t i = 0; i < 3; i++)
          c_w_face_cog[c_id][i] += (  v_w_ref[c_id][i]
                                    + v_w_inteface[0][i]
                                    + v_w_inteface[1][i])
                                 * cs_math_3_norm(vn);

        if (v_w_inteface != (cs_real_3_t *)_v_w_inteface)
          BFT_FREE(v_w_inteface);
      }

      cs_lnum_t n_f_face_vertices = f_vtx_count;

      /* Computing quantities of fluid part of interface faces */
      if (n_f_face_vertices > 2) {

        cs_lnum_t *_f_face_pos[8];
        cs_lnum_t *f_face_pos = (cs_lnum_t *)_f_face_pos;
        if (n_f_face_vertices > 8)
          BFT_MALLOC(f_face_pos, n_f_face_vertices, cs_lnum_t);

        for (cs_lnum_t i = 0; i < n_f_face_vertices; i++) {
          f_face_pos[i] = i;
        }

        cs_lnum_2_t f_vtx_idx = {0, n_f_face_vertices};
        cs_real_t face_cog[1][3]  = {{0., 0., 0.}};
        cs_real_t face_normal[1][3] = {{0., 0., 0.}};

        cs_real_t *_f_vtx_coord_l[8][3];
        cs_real_3_t *f_vtx_coord_l = (cs_real_3_t *)_f_vtx_coord_l;
        if (n_f_face_vertices > 8)
          BFT_MALLOC(f_vtx_coord_l, n_f_face_vertices, cs_real_3_t);

        for (cs_lnum_t ffv = 0; ffv < n_f_face_vertices; ffv++) {
          for (cs_lnum_t i = 0; i < 3; i++)
            f_vtx_coord_l[ffv][i] = f_vtx_coord[ffv][i];
        }

        _compute_face_quantities(1,
                                 f_vtx_coord_l,
                                 f_vtx_idx,
                                 f_face_pos,
                                 face_cog,
                                 face_normal);

        if (f_vtx_coord_l != (cs_real_3_t *)_f_vtx_coord_l)
          BFT_FREE(f_vtx_coord_l);

        /* Adjust porosity and COG */

        for (cs_lnum_t i = 0; i < 3; i++) {
          b_f_face_normal[face_id][i] = face_normal[0][i];
          b_f_face_cog[face_id][i] = face_cog[0][i];
        }
        mq->b_f_face_surf[face_id] = cs_math_3_norm(b_f_face_normal[face_id]);

        if (f_face_pos != (cs_lnum_t *)_f_face_pos)
          BFT_FREE(f_face_pos);
      }

      if (vc != (cs_real_3_t *)_vc)
        BFT_FREE(vc);

      if (f_vtx_coord != (cs_real_3_t *)_f_vtx_coord)
        BFT_FREE(f_vtx_coord);

    } /* Loop on boundary faces */

    /* Set center of gravity and normal/surface of solid faces */

    if (   cs_math_3_norm(c_w_face_normal[c_id]) > 0.
        && sum_surf[c_id] > 0.) {
      c_w_face_surf[c_id] = sum_surf[c_id];

      for (cs_lnum_t i = 0; i < 3; i++) {
        c_w_face_cog[c_id][i]    *= 1./(6.*c_w_face_surf[c_id]);
        /* From now the normal is not unitary anymore */
        c_w_face_normal[c_id][i] *= c_w_face_surf[c_id];
      }

    }
    /* Plan outside the cell (it is caused due to tolerance of ple routine) */
    else {
      c_w_face_surf[c_id] = 0.;
      for (cs_lnum_t i = 0; i < 3; i++) {
        c_w_face_cog[c_id][i] = 0.;
        c_w_face_normal[c_id][i] = 0.;
      }
    }

    if (n_w_vtx > 0)
      n_ib_cells += 1;

    w_vtx_s_id += n_w_vtx;
    w_vtx_idx[c_id+1] = w_vtx_s_id;

  } /* End loop on cells */

  BFT_FREE(sum_surf);
  BFT_FREE(v_w_ref);
  BFT_FREE(flag_id);
  BFT_REALLOC(w_vtx, w_vtx_idx[m->n_cells], cs_real_3_t);

  /* Synchronize geometrics quantities with face duality */

  if (cs_glob_n_ranks > 1 || m->n_init_perio > 0) {

    /* Build global interior faces interface */

    int n_perio = m->n_init_perio;
    int *perio_num = nullptr;
    cs_lnum_t *n_perio_face_couples = nullptr;
    cs_gnum_t **perio_face_couples = nullptr;

    if (n_perio > 0) {
      BFT_MALLOC(perio_num, n_perio, int);
      for (int i = 0; i < n_perio; i++)
        perio_num[i] = i+1;

      cs_mesh_get_perio_faces(m,
                              &n_perio_face_couples,
                              &perio_face_couples);
    }

    cs_interface_set_t *f_if
      = cs_interface_set_create(m->n_i_faces,
                                nullptr,
                                m->global_i_face_num,
                                m->periodicity,
                                n_perio,
                                perio_num,
                                n_perio_face_couples,
                                (const cs_gnum_t *const *)perio_face_couples);

    if (n_perio > 0) {
      for (int i = 0; i < n_perio; i++)
        BFT_FREE(perio_face_couples[i]);

      BFT_FREE(perio_face_couples);
      BFT_FREE(n_perio_face_couples);
      BFT_FREE(perio_num);
    }

    /* Synchronization */
    cs_interface_set_max(f_if,
                         m->n_i_faces,
                         2*3,
                         true,
                         CS_REAL_TYPE,
                         (cs_real_23_t *)i_f_face_cell_normal);

    cs_interface_set_max(f_if,
                         m->n_i_faces,
                         2*3,
                         true,
                         CS_REAL_TYPE,
                         (cs_real_23_t *)i_f_face_cog_dual);

    cs_interface_set_destroy(&f_if);
  }

  /* Adjusting face porosity and COG
   * we take the smaller surface side and its COG */

  for (cs_lnum_t f_id = 0; f_id < m->n_i_faces; f_id++) {

    cs_real_t area0 = cs_math_3_norm(i_f_face_cell_normal[f_id][0]);
    cs_real_t area1 = cs_math_3_norm(i_f_face_cell_normal[f_id][1]);

    if (area0 < area1) {
      for (cs_lnum_t i = 0; i < 3; i++) {
        i_f_face_cog[f_id][i] = i_f_face_cog_dual[f_id][0][i];
        i_f_face_normal[f_id][i] = i_f_face_cell_normal[f_id][0][i];
      }
    }
    else {
      for (cs_lnum_t i = 0; i < 3; i++) {
        i_f_face_cog[f_id][i] = i_f_face_cog_dual[f_id][1][i];
        i_f_face_normal[f_id][i] = i_f_face_cell_normal[f_id][1][i];
      }
    }

    mq->i_f_face_surf[f_id] = cs_math_3_norm(i_f_face_normal[f_id]);

  }

  //TODO merge with regular function _compute_cell_quantities

  /* Compute fluid cell centers from face barycenters
   * and fluid volumes */
  _compute_fluid_solid_cell_quantities(m,
                                       (const cs_real_23_t *)i_f_face_cell_normal,
                                       (const cs_real_23_t *)i_f_face_cog_dual,
                                       (const cs_real_3_t *)mq->b_f_face_normal,
                                       (const cs_real_3_t *)mq->b_f_face_cog,
                                       (const cs_real_3_t *)mq->c_w_face_normal,
                                       (const cs_real_3_t *)mq->c_w_face_cog,
                                       (cs_real_3_t *)mq->cell_f_cen,
                                       mq->cell_f_vol);

  /* Correction of small or negative volumes
     (doesn't conserve the total volume) */

  if (cs_glob_mesh_quantities_flag & CS_CELL_VOLUME_RATIO_CORRECTION)
    _cell_bad_volume_correction(m, mq->cell_f_vol);

  /* Synchronize geometric quantities */

  if (m->halo != nullptr) {

    cs_halo_sync_var_strided(m->halo, CS_HALO_EXTENDED,
                             mq->cell_f_cen, 3);
    if (m->n_init_perio > 0)
      cs_halo_perio_sync_coords(m->halo, CS_HALO_EXTENDED,
                                mq->cell_f_cen);

    cs_halo_sync_var(m->halo, CS_HALO_EXTENDED, mq->cell_f_vol);

  }

  /* Update geometrical entities: distances to the cell center of gravity.
   * Note: only direction is used for "i_face_normal" and "b_face_normal"
   * */

  _compute_face_distances(m->n_i_faces,
                          m->n_b_faces,
                          (const cs_lnum_2_t *)(m->i_face_cells),
                          (const cs_lnum_t   *)(m->b_face_cells),
                          (const cs_real_3_t *)(mq->i_face_u_normal),
                          (const cs_real_3_t *)(mq->i_f_face_normal),
                          (const cs_real_3_t *)(mq->b_face_u_normal),
                          (const cs_real_3_t *)(mq->b_f_face_normal),
                          (const cs_real_3_t *)(mq->i_f_face_cog),
                          (const cs_real_3_t *)(mq->b_f_face_cog),
                          (const cs_real_3_t *)(mq->cell_f_cen),
                          (const cs_real_t *)(mq->cell_f_vol),
                          mq->i_dist,
                          mq->b_dist,
                          mq->weight);

  _compute_face_vectors(m->dim,
                        m->n_i_faces,
                        m->n_b_faces,
                        (const cs_lnum_2_t *)(m->i_face_cells),
                        m->b_face_cells,
                        (const cs_real_3_t *)mq->i_face_u_normal,
                        (const cs_real_3_t *)mq->b_face_u_normal,
                        mq->i_f_face_cog,
                        mq->b_f_face_cog,
                        mq->cell_f_cen,
                        mq->weight,
                        mq->b_dist,
                        mq->dijpf,
                        mq->diipb,
                        mq->dofij);

  _compute_face_sup_vectors(m->n_cells,
                            m->n_i_faces,
                            (const cs_lnum_2_t *)(m->i_face_cells),
                            (const cs_real_3_t *)(mq->i_face_u_normal),
                            (const cs_real_3_t *)(mq->i_face_normal),
                            (const cs_real_3_t *)(mq->i_f_face_cog),
                            (const cs_real_3_t *)(mq->cell_f_cen),
                            mq->cell_f_vol,
                            mq->i_dist,
                            (cs_real_3_t *)(mq->diipf),
                            (cs_real_3_t *)(mq->djjpf));

  /* Reinitializing the wall normal before correcting */
  cs_array_real_set_scalar(3*m->n_cells_with_ghosts,
                           0.,
                           (cs_real_t *)c_w_face_normal);

  cs_porosity_from_scan_opt_t *poro_from_scan = cs_glob_porosity_from_scan_opt;
  /* Update the cell porosity field value and penalize small fluid cells */
  cs_real_t poro_threshold = poro_from_scan->porosity_threshold;

  for (cs_lnum_t c_id = 0; c_id < m->n_cells; c_id++) {
    cs_real_t porosity = mq->cell_f_vol[c_id]/mq->cell_vol[c_id];
    if (porosity > 1.)
      porosity = 1.;

    f_poro->val[c_id] = porosity;
    if (porosity > poro_threshold)
      mq->c_disable_flag[c_id] = 0;

    /* Penalize ibm cells with small porosity */
    else if (porosity < poro_threshold && porosity > cs_math_epzero) {

      mq->c_disable_flag[c_id] = 1;

      f_poro->val[c_id] = 0.;
      mq->cell_f_vol[c_id] = 0.0;
      mq->c_w_face_surf[c_id] = 0.0;

      for (cs_lnum_t i = 0; i < 3; i++)
        c_w_face_normal[c_id][i] = 0.0;

      /* Interior faces */
      const cs_lnum_t s_id_i = c2c_idx[c_id];
      const cs_lnum_t e_id_i = c2c_idx[c_id+1];

      /* Loop on interior faces of cell c_id */
      for (cs_lnum_t cidx = s_id_i; cidx < e_id_i; cidx++) {
        const cs_lnum_t face_id = cell_i_faces[cidx];

        i_f_face_normal[face_id][0] = 0.;
        i_f_face_normal[face_id][1] = 0.;
        i_f_face_normal[face_id][2] = 0.;
        mq->i_f_face_surf[face_id] = 0.;
      }

      /* Boundary faces */
      const cs_lnum_t s_id_b = cell_b_faces_idx[c_id];
      const cs_lnum_t e_id_b = cell_b_faces_idx[c_id+1];

      for (cs_lnum_t cidx = s_id_b; cidx < e_id_b; cidx++) {
        const cs_lnum_t face_id = cell_b_faces[cidx];

        b_f_face_normal[face_id][0] = 0.;
        b_f_face_normal[face_id][1] = 0.;
        b_f_face_normal[face_id][2] = 0.;
      }
    }
  }

  /* Compute solid normal from fluid normals */
  cs_real_33_t *xpsn;
  BFT_MALLOC(xpsn, m->n_cells_with_ghosts, cs_real_33_t);
  cs_array_real_set_scalar(9*m->n_cells_with_ghosts,
                           0.,
                           (cs_real_t *)xpsn);

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    /* Interior faces */
    const cs_lnum_t s_id_i = c2c_idx[c_id];
    const cs_lnum_t e_id_i = c2c_idx[c_id+1];

    bool is_active_cell = false;

    /* Loop on interior faces of cell c_id */
    for (cs_lnum_t cidx = s_id_i; cidx < e_id_i; cidx++) {
      const cs_lnum_t f_id = cell_i_faces[cidx];
      const short int sign = cell_i_faces_sgn[cidx];

      const cs_real_t face_norm = cs_math_3_norm(i_f_face_normal[f_id]);
      if (face_norm > 0.0)
        is_active_cell = true;

      for (cs_lnum_t i = 0; i < 3; i++) {

        c_w_face_normal[c_id][i] -= sign*i_f_face_normal[f_id][i];

        const cs_real_t xfmxc
          = (mq->i_f_face_cog[3*f_id+i] - mq->cell_f_cen[3*c_id+i]);

        for (cs_lnum_t j = 0; j < 3; j++)
          xpsn[c_id][i][j] += sign*xfmxc*i_f_face_normal[f_id][j];
      }
    } /* End loop on adjacent cells */

    /* Penalizes IBM cell when its neighbors have zero fluid surfaces */
    if (!(is_active_cell)) {

      mq->c_disable_flag[c_id] = 1;

      /* Forced porosity to 0 */
      mq->cell_f_vol[c_id] = 0.0;
      f_poro->val[c_id] = 0.0;
      mq->c_w_face_surf[c_id] = 0.0;

      for (cs_lnum_t i = 0; i < 3; i++)
        c_w_face_normal[c_id][i] = 0.0;

      /* If no internal contribution, the cell is isolated and we
         skip the boundary contribution */
      continue;
    }

    /* Boundary faces */
    const cs_lnum_t s_id_b = cell_b_faces_idx[c_id];
    const cs_lnum_t e_id_b = cell_b_faces_idx[c_id+1];

    /* Loop on boundary faces of cell c_id */
    for (cs_lnum_t cidx = s_id_b; cidx < e_id_b; cidx++) {
      const cs_lnum_t f_id = cell_b_faces[cidx];

      for (cs_lnum_t i = 0; i < 3; i++) {
        c_w_face_normal[c_id][i] -= b_f_face_normal[f_id][i];

        const cs_real_t xfmxc
          = (mq->b_f_face_cog[3*f_id+i] - mq->cell_f_cen[3*c_id+i]);

        for (cs_lnum_t j = 0; j < 3; j++)
          xpsn[c_id][i][j] += xfmxc*b_f_face_normal[f_id][j];
      }
    } /* End loop on boundary cells */

  } /* End loop on cells */

  /* Correction of solid face center and distance to the immersed wall */

  cs_real_t *c_w_dist_inv = mq->c_w_dist_inv;

  const cs_real_t m_identity[3][3] = {{1., 0., 0.,}, {0., 1., 0.}, {0., 0., 1.}};

  for (cs_lnum_t c_id = 0; c_id < m->n_cells; c_id++) {

    cs_real_t xc[3] = {mq->cell_f_cen[3*c_id],
                       mq->cell_f_cen[3*c_id + 1],
                       mq->cell_f_cen[3*c_id + 2]};
    cs_real_t pyr_vol = cs_math_3_distance_dot_product(xc,
                                                       c_w_face_cog[c_id],
                                                       c_w_face_normal[c_id]);
    cs_real_t vol_min = cs_math_epzero*mq->cell_f_vol[c_id];

    c_w_face_surf[c_id] = cs_math_3_norm(c_w_face_normal[c_id]);

    /* Compare the volume of the pyramid (xw - xc).Sw to the fluid
     *  volume */
    if (pyr_vol > vol_min) {

      cs_real_t vc_w_f_cen[3];
      cs_real_t c_w_normal_unit[3];
      cs_math_3_normalize(c_w_face_normal[c_id], c_w_normal_unit);

      /* Correction of solid plane center xw-xc */

      cs_real_t mat[3][3];
      for (cs_lnum_t i = 0; i < 3; i++) {
        for (cs_lnum_t j = 0; j < 3; j++) {
          mat[i][j] = mq->cell_f_vol[c_id]*m_identity[i][j]
                    - xpsn[c_id][i][j];
        }
      }
      cs_math_33_3_product(mat, c_w_face_normal[c_id], vc_w_f_cen);
      cs_real_t d_w = cs_math_3_square_norm(c_w_face_normal[c_id]);
      if (d_w > DBL_MIN)
        d_w = 1./d_w;

      for (cs_lnum_t i = 0; i < 3; i++)
        vc_w_f_cen[i] *= d_w;

      for (cs_lnum_t i = 0; i < 3; i++)
        c_w_face_cog[c_id][i] = vc_w_f_cen[i] + mq->cell_f_cen[c_id*3+i];

      /* Distance to the immersed wall */

      cs_real_t d_dist = cs_math_3_dot_product(vc_w_f_cen, c_w_normal_unit);
      if (d_dist > DBL_MIN)
        d_dist = 1./d_dist;

      c_w_dist_inv[c_id] = d_dist;

      /* else c_w_dist_inv = 0 by default as it is a field */
    }
  }

  /* Synchronization */
  if (m->halo != nullptr) {
    cs_halo_sync_var_strided(m->halo, CS_HALO_EXTENDED,
                             (cs_real_t *)c_w_face_normal, 3);

    cs_halo_sync_var(m->halo, CS_HALO_EXTENDED, mq->c_w_dist_inv);
    cs_halo_sync_var(m->halo, CS_HALO_EXTENDED, mq->c_w_face_surf);

    cs_halo_sync_var_strided(m->halo, CS_HALO_EXTENDED,
                             (cs_real_t *)mq->c_w_face_cog, 3);

  }

  /* Connectivity ib_cell to cells */
  cs_lnum_t *ibcell_cells; // TODO: use only b_face_cells
  BFT_MALLOC(ibcell_cells, n_ib_cells, cs_lnum_t);

  cs_lnum_t *face_vertex_idx;
  BFT_MALLOC(face_vertex_idx, n_ib_cells + 1, cs_lnum_t);
  face_vertex_idx[0] = 0;

  cs_lnum_t n_ib_cells0 = 0;
  cs_lnum_t n_glob_vtx = 0;

  /* face_vertex_idx and ibcell_cells computation */
  for (cs_lnum_t c_id = 0; c_id < m->n_cells; c_id++) {
    const int n_vtx = w_vtx_idx[c_id+1]-w_vtx_idx[c_id];

    if (n_vtx > 0) {

      if (n_vtx < 6) {
        bft_error(__FILE__, __LINE__, 0,
                  _("Number of vertex to construct the IBM approximate"
                    " plane \n in cell id = %d equals %d."
                    " It should be greater than 6."),
                  c_id, n_vtx);
      }

      if (n_vtx%2 != 0) {
        bft_error(__FILE__, __LINE__, 0,
                  _("Vertex number must be even. Cell id = %d, n_vtx = %d"),
                  c_id, n_vtx);
      }

      ibcell_cells[n_ib_cells0] = c_id;
      n_ib_cells0 += 1;
      n_glob_vtx += 0.5*n_vtx;
      face_vertex_idx[n_ib_cells0] = 0.5*w_vtx_idx[c_id+1];
    }
  }

  _post_plane_ib(n_ib_cells,
                 n_glob_vtx,
                 ibcell_cells,
                 vtx_ids,
                 w_vtx_idx,
                 face_vertex_idx,
                 w_vtx);

  /* Free memory */
  BFT_FREE(face_vertex_idx);
  BFT_FREE(ibcell_cells);
  BFT_FREE(_i_f_surf);

  BFT_FREE(i_f_face_cell_normal);
  BFT_FREE(xpsn);

  BFT_FREE(w_vtx_idx);
  BFT_FREE(vtx_ids);
  BFT_FREE(w_vtx);

  BFT_FREE(i_f_face_cog_dual);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute mesh quantities.
 *
 * \param[in]       m   pointer to mesh structure
 * \param[in, out]  mq  pointer to mesh quantities structures.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_quantities_compute(const cs_mesh_t       *m,
                           cs_mesh_quantities_t  *mq)
{
  cs_lnum_t  dim = m->dim;
  cs_lnum_t  n_i_faces = m->n_i_faces;
  cs_lnum_t  n_b_faces = m->n_b_faces;
  cs_lnum_t  n_cells_with_ghosts = m->n_cells_with_ghosts;

  const cs_alloc_mode_t amode = cs_alloc_mode_read_mostly;

  /* Update the number of passes */

  _n_computations++;

  cs_mesh_quantities_compute_preprocess(m, mq);

  /* Fluid surfaces and volumes: point to standard quantities and
   * may be modified afterwards */
  mq->i_f_face_normal = mq->i_face_normal;
  mq->b_f_face_normal = mq->b_face_normal;
  mq->i_f_face_surf = mq->i_face_surf;
  mq->b_f_face_surf = mq->b_face_surf;
  mq->i_f_face_cog = mq->i_face_cog;
  mq->b_f_face_cog = mq->b_face_cog;

  mq->cell_f_vol = mq->cell_vol;
  mq->cell_f_cen = mq->cell_cen;
  mq->cell_s_cen = nullptr;

  /* Porous models */
  if (mq->c_disable_flag == nullptr) {
    if (mq->has_disable_flag == 1) {
      cs_lnum_t n_cells_ext = n_cells_with_ghosts;
      CS_MALLOC_HD(mq->c_disable_flag, n_cells_ext, int, amode);
      for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++)
        mq->c_disable_flag[cell_id] = 0;
    }
    else {
      CS_MALLOC_HD(mq->c_disable_flag, 1, int, amode);
      mq->c_disable_flag[0] = 0;
    }
    cs_mem_advise_set_read_mostly(mq->c_disable_flag);
  }

  mq->min_f_vol = mq->min_vol;
  mq->max_f_vol = mq->max_vol;
  mq->tot_f_vol = mq->tot_vol;

  if (mq->i_dist == nullptr) {
    CS_MALLOC_HD(mq->i_dist, n_i_faces, cs_real_t, amode);
    cs_mem_advise_set_read_mostly(mq->i_dist);
  }

  if (mq->b_dist == nullptr) {
    CS_MALLOC_HD(mq->b_dist, n_b_faces, cs_real_t, amode);
    cs_mem_advise_set_read_mostly(mq->b_dist);
  }

  if (mq->weight == nullptr) {
    CS_MALLOC_HD(mq->weight, n_i_faces, cs_real_t, amode);
    cs_mem_advise_set_read_mostly(mq->weight);
  }

  if (mq->i_f_weight == nullptr) {
    CS_MALLOC_HD(mq->i_f_weight, n_i_faces, cs_real_t, amode);
    cs_mem_advise_set_read_mostly(mq->i_f_weight);
  }

  if (mq->dijpf == nullptr) {
    CS_MALLOC_HD(mq->dijpf, n_i_faces*dim, cs_real_t, amode);
    cs_mem_advise_set_read_mostly(mq->dijpf);
  }

  if (mq->diipb == nullptr) {
    CS_MALLOC_HD(mq->diipb, n_b_faces*dim, cs_real_t, amode);
    cs_mem_advise_set_read_mostly(mq->diipb);
  }

  if (mq->dofij == nullptr) {
    CS_MALLOC_HD(mq->dofij, n_i_faces*dim, cs_real_t, amode);
    cs_mem_advise_set_read_mostly(mq->dofij);
  }

  if (mq->diipf == nullptr) {
    CS_MALLOC_HD(mq->diipf, n_i_faces*dim, cs_real_t, amode);
    cs_mem_advise_set_read_mostly(mq->diipf);
  }

  if (mq->djjpf == nullptr) {
    CS_MALLOC_HD(mq->djjpf, n_i_faces*dim, cs_real_t, amode);
    cs_mem_advise_set_read_mostly(mq->djjpf);
  }

  if (mq->b_sym_flag == nullptr) {
    CS_MALLOC_HD(mq->b_sym_flag, n_b_faces, int, amode);
    for (cs_lnum_t i = 0; i < n_b_faces; i++)
      mq->b_sym_flag[i] = 1;
    cs_mem_advise_set_read_mostly(mq->b_sym_flag);
  }

  /* Compute some distances relative to faces and associated weighting */

  _compute_face_distances(m->n_i_faces,
                          m->n_b_faces,
                          (const cs_lnum_2_t *)(m->i_face_cells),
                          m->b_face_cells,
                          (const cs_real_3_t *)(mq->i_face_u_normal),
                          (const cs_real_3_t *)(mq->i_face_normal),
                          (const cs_real_3_t *)(mq->b_face_u_normal),
                          (const cs_real_3_t *)(mq->b_face_normal),
                          (const cs_real_3_t *)(mq->i_face_cog),
                          (const cs_real_3_t *)(mq->b_face_cog),
                          (const cs_real_3_t *)(mq->cell_cen),
                          (const cs_real_t *)(mq->cell_vol),
                          mq->i_dist,
                          mq->b_dist,
                          mq->weight);

  /* Compute some vectors relative to faces to handle non-orthogonalities */

  _compute_face_vectors(dim,
                        m->n_i_faces,
                        m->n_b_faces,
                        (const cs_lnum_2_t *)(m->i_face_cells),
                        m->b_face_cells,
                        (const cs_real_3_t *)mq->i_face_u_normal,
                        (const cs_real_3_t *)mq->b_face_u_normal,
                        mq->i_face_cog,
                        mq->b_face_cog,
                        mq->cell_cen,
                        mq->weight,
                        mq->b_dist,
                        mq->dijpf,
                        mq->diipb,
                        mq->dofij);

  /* Compute additional vectors relative to faces
     to handle non-orthogonalities */

  _compute_face_sup_vectors
    (m->n_cells,
     m->n_i_faces,
     (const cs_lnum_2_t *)(m->i_face_cells),
     (const cs_real_3_t *)(mq->i_face_u_normal),
     (const cs_real_3_t *)(mq->i_face_normal),
     (const cs_real_3_t *)(mq->i_face_cog),
     (const cs_real_3_t *)(mq->cell_cen),
     mq->cell_vol,
     mq->i_dist,
     (cs_real_3_t *)(mq->diipf),
     (cs_real_3_t *)(mq->djjpf));

  /* Build the geometrical matrix linear gradient correction */
  if (cs_glob_mesh_quantities_flag & CS_BAD_CELLS_WARPED_CORRECTION)
    _compute_corr_grad_lin(m, mq);

  /* Print some information on the control volumes, and check min volume */

  if (_n_computations == 1)
    bft_printf(_(" --- Information on the volumes\n"
                 "       Minimum control volume      = %14.7e\n"
                 "       Maximum control volume      = %14.7e\n"
                 "       Total volume for the domain = %14.7e\n"),
               mq->min_vol, mq->max_vol,
               mq->tot_vol);
  else {
    if (mq->min_vol <= 0.) {
      bft_printf(_(" --- Information on the volumes\n"
                   "       Minimum control volume      = %14.7e\n"
                   "       Maximum control volume      = %14.7e\n"
                   "       Total volume for the domain = %14.7e\n"),
                 mq->min_vol, mq->max_vol,
                 mq->tot_vol);
      bft_printf(_("\nAbort due to the detection of a negative control "
                   "volume.\n"));
    }
  }
}

/*----------------------------------------------------------------------------
 * Compute min, max, and total
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   mesh_quantities <-> pointer to a cs_mesh_quantities_t structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_fluid_compute(const cs_mesh_t       *mesh,
                                 cs_mesh_quantities_t  *mesh_quantities)
{
  CS_UNUSED(mesh);
  CS_UNUSED(mesh_quantities);
}

/*----------------------------------------------------------------------------
 * Compute the total, min, and max fluid volumes of cells
 *
 * parameters:
 *   mesh            <-- pointer to mesh structure
 *   mesh_quantities <-> pointer to a mesh quantities structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_fluid_vol_reductions(const cs_mesh_t       *mesh,
                                        cs_mesh_quantities_t  *mesh_quantities)
{
  _cell_volume_reductions(mesh,
                          mesh_quantities->cell_f_vol,
                          &(mesh_quantities->min_f_vol),
                          &(mesh_quantities->max_f_vol),
                          &(mesh_quantities->tot_f_vol));

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {

    cs_real_t  _min_f_vol, _max_f_vol, _tot_f_vol;

    MPI_Allreduce(&(mesh_quantities->min_f_vol), &_min_f_vol, 1, CS_MPI_REAL,
                  MPI_MIN, cs_glob_mpi_comm);

    MPI_Allreduce(&(mesh_quantities->max_f_vol), &_max_f_vol, 1, CS_MPI_REAL,
                  MPI_MAX, cs_glob_mpi_comm);

    MPI_Allreduce(&(mesh_quantities->tot_f_vol), &_tot_f_vol, 1, CS_MPI_REAL,
                  MPI_SUM, cs_glob_mpi_comm);

    mesh_quantities->min_f_vol = _min_f_vol;
    mesh_quantities->max_f_vol = _max_f_vol;
    mesh_quantities->tot_f_vol = _tot_f_vol;

  }
#endif
}

/*----------------------------------------------------------------------------
 * Compute fluid section mesh quantities at the initial step
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   mesh_quantities <-> pointer to a cs_mesh_quantities_t structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_init_fluid_sections(const cs_mesh_t       *mesh,
                            cs_mesh_quantities_t  *mesh_quantities)
{
  cs_lnum_t  n_i_faces = mesh->n_i_faces;
  cs_lnum_t  n_b_faces = mesh->n_b_faces;

  cs_real_3_t *restrict i_face_normal =
    (cs_real_3_t *)mesh_quantities->i_face_normal;
  cs_real_3_t *restrict b_face_normal =
    (cs_real_3_t *)mesh_quantities->b_face_normal;
  cs_real_3_t *restrict i_f_face_normal =
    (cs_real_3_t *)mesh_quantities->i_f_face_normal;
  cs_real_3_t *restrict b_f_face_normal =
    (cs_real_3_t *)mesh_quantities->b_f_face_normal;

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
    mesh_quantities->i_f_face_surf[face_id]
      = mesh_quantities->i_face_surf[face_id];

    for (cs_lnum_t i = 0; i < 3; i++)
      i_f_face_normal[face_id][i] = i_face_normal[face_id][i];

    mesh_quantities->i_f_face_factor[face_id][0] = 1.;
    mesh_quantities->i_f_face_factor[face_id][1] = 1.;
  }

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    mesh_quantities->b_f_face_surf[face_id]
      = mesh_quantities->b_face_surf[face_id];

    for (cs_lnum_t i = 0; i < 3; i++)
      b_f_face_normal[face_id][i] = b_face_normal[face_id][i];

    mesh_quantities->b_f_face_factor[face_id] = 1.;
  }
}

/*----------------------------------------------------------------------------
 * Compute mesh quantities -> vectors II' and JJ'
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   mesh_quantities <-> pointer to a cs_mesh_quantities_t structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_sup_vectors(const cs_mesh_t       *mesh,
                               cs_mesh_quantities_t  *mesh_quantities)
{
  cs_lnum_t  dim = mesh->dim;
  cs_lnum_t  n_i_faces = mesh->n_i_faces;

  if (mesh_quantities->diipf == nullptr)
    BFT_MALLOC(mesh_quantities->diipf, n_i_faces*dim, cs_real_t);

  if (mesh_quantities->djjpf == nullptr)
    BFT_MALLOC(mesh_quantities->djjpf, n_i_faces*dim, cs_real_t);

  _compute_face_sup_vectors
    (mesh->n_cells,
     mesh->n_i_faces,
     (const cs_lnum_2_t *)(mesh->i_face_cells),
     (const cs_real_3_t *)(mesh_quantities->i_face_u_normal),
     (const cs_real_3_t *)(mesh_quantities->i_face_normal),
     (const cs_real_3_t *)(mesh_quantities->i_face_cog),
     (const cs_real_3_t *)(mesh_quantities->cell_cen),
     mesh_quantities->cell_vol,
     mesh_quantities->i_dist,
     (cs_real_3_t *)(mesh_quantities->diipf),
     (cs_real_3_t *)(mesh_quantities->djjpf));
}

/*----------------------------------------------------------------------------
 * Compute internal and border face normal.
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   p_i_face_normal <-> pointer to the internal face normal array
 *   p_b_face_normal <-> pointer to the border face normal array
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_face_normal(const cs_mesh_t   *mesh,
                               cs_real_t         *p_i_face_normal[],
                               cs_real_t         *p_b_face_normal[])
{
  cs_real_t  *i_face_normal = nullptr, *b_face_normal = nullptr;

  const cs_lnum_t  n_b_faces = mesh->n_b_faces;
  const cs_lnum_t  n_i_faces = mesh->n_i_faces;

  /* Internal face treatment */

  BFT_MALLOC(i_face_normal, n_i_faces * 3, cs_real_t);

  _compute_face_normal(mesh->n_i_faces,
                       (const cs_real_3_t *)mesh->vtx_coord,
                       mesh->i_face_vtx_idx,
                       mesh->i_face_vtx_lst,
                       (cs_real_3_t *)i_face_normal);

  *p_i_face_normal = i_face_normal;

  /* Boundary face treatment */

  BFT_MALLOC(b_face_normal, n_b_faces * 3, cs_real_t);

  _compute_face_normal(mesh->n_b_faces,
                       (const cs_real_3_t *)mesh->vtx_coord,
                       mesh->b_face_vtx_idx,
                       mesh->b_face_vtx_lst,
                       (cs_real_3_t *)b_face_normal);

  *p_b_face_normal = b_face_normal;
}

/*----------------------------------------------------------------------------
 * Compute interior face centers and normals.
 *
 * The corresponding arrays are allocated by this function, and it is the
 * caller's responsibility to free them when they are no longer needed.
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   p_i_face_cog    <-> pointer to the interior face center array
 *   p_i_face_normal <-> pointer to the interior face normal array
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_i_faces(const cs_mesh_t   *mesh,
                           cs_real_t         *p_i_face_cog[],
                           cs_real_t         *p_i_face_normal[])
{
  cs_real_t  *i_face_cog = nullptr, *i_face_normal = nullptr;

  BFT_MALLOC(i_face_cog, mesh->n_i_faces * mesh->dim, cs_real_t);
  BFT_MALLOC(i_face_normal, mesh->n_i_faces * mesh->dim, cs_real_t);

  _compute_face_quantities(mesh->n_i_faces,
                          (const cs_real_3_t *)mesh->vtx_coord,
                          mesh->i_face_vtx_idx,
                          mesh->i_face_vtx_lst,
                          (cs_real_3_t *)i_face_cog,
                          (cs_real_3_t *)i_face_normal);

  *p_i_face_cog = i_face_cog;
  *p_i_face_normal = i_face_normal;
}

/*----------------------------------------------------------------------------
 * Compute border face centers and normals.
 *
 * The corresponding arrays are allocated by this function, and it is the
 * caller's responsibility to free them when they are no longer needed.
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   p_b_face_cog    <-> pointer to the border face center array
 *   p_b_face_normal <-> pointer to the border face normal array
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_b_faces(const cs_mesh_t   *mesh,
                           cs_real_t         *p_b_face_cog[],
                           cs_real_t         *p_b_face_normal[])
{
  cs_real_t  *b_face_cog = nullptr, *b_face_normal = nullptr;

  BFT_MALLOC(b_face_cog, mesh->n_b_faces * mesh->dim, cs_real_t);
  BFT_MALLOC(b_face_normal, mesh->n_b_faces * mesh->dim, cs_real_t);

  _compute_face_quantities(mesh->n_b_faces,
                           (const cs_real_3_t *)mesh->vtx_coord,
                           mesh->b_face_vtx_idx,
                           mesh->b_face_vtx_lst,
                           (cs_real_3_t *)b_face_cog,
                           (cs_real_3_t *)b_face_normal);

  *p_b_face_cog = b_face_cog;
  *p_b_face_normal = b_face_normal;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute approximate cells centers as the mean of the given face
 *         centers weighted by the associated surfaces.
 *
 *           n-1
 *           Sum  Surf(Fi) G(Fi)
 *           i=0
 *  G(C) = -----------------------
 *           n-1
 *           Sum  Surf(Fi)
 *           i=0
 *
 * \param[in]   mesh         pointer to mesh structure
 * \param[in]   i_face_norm  surface normal of internal faces
 * \param[in]   i_face_cog   center of gravity of internal faces
 * \param[in]   b_face_norm  surface normal of border faces
 * \param[in]   b_face_cog   center of gravity of border faces
 * \param[out]  cell_cen     cell centers
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_quantities_cell_faces_cog(const cs_mesh_t  *mesh,
                                  const cs_real_t   i_face_norm[],
                                  const cs_real_t   i_face_cog[],
                                  const cs_real_t   b_face_norm[],
                                  const cs_real_t   b_face_cog[],
                                  cs_real_t         cell_cen[])
{
  cs_real_t  *cell_area = nullptr;

  /* Mesh connectivity */

  const  cs_lnum_t  n_i_faces = mesh->n_i_faces;
  const  cs_lnum_t  n_b_faces = mesh->n_b_faces;
  const  cs_lnum_t  n_cells = mesh->n_cells;
  const  cs_lnum_t  n_cells_with_ghosts = mesh->n_cells_with_ghosts;
  const  cs_lnum_2_t  *i_face_cells
    = (const cs_lnum_2_t *)(mesh->i_face_cells);
  const  cs_lnum_t  *b_face_cells = mesh->b_face_cells;

  /* Return if ther is not enough data (Solcom case except rediative module
     or Pre-processor 1.2.d without option "-n") */

  if (mesh->i_face_vtx_lst == nullptr && mesh->b_face_vtx_lst == nullptr)
    return;

  /* Checking */

  assert(cell_cen != nullptr);

  /* Initialization */

  BFT_MALLOC(cell_area, n_cells_with_ghosts, cs_real_t);

  for (cs_lnum_t j = 0; j < n_cells_with_ghosts; j++) {

    cell_area[j] = 0.;

    for (cs_lnum_t i = 0; i < 3; i++)
      cell_cen[3*j + i] = 0.;

  }

  /* Loop on interior faces
     ---------------------- */

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

    /* For each cell sharing the internal face, we update
     * cell_cen and cell_area */

    cs_lnum_t c_id1 = i_face_cells[f_id][0];
    cs_lnum_t c_id2 = i_face_cells[f_id][1];

    /* Computation of the area of the face */

    cs_real_t area = cs_math_3_norm(i_face_norm + 3*f_id);

    if (   !(cs_glob_mesh_quantities_flag & CS_FACE_NULL_SURFACE)
        || area > 1.e-20) {
      if (c_id1 > -1) {
        cell_area[c_id1] += area;
        for (cs_lnum_t i = 0; i < 3; i++)
          cell_cen[3*c_id1 + i] += i_face_cog[3*f_id + i]*area;
      }
      if (c_id2 > -1) {
        cell_area[c_id2] += area;
        for (cs_lnum_t i = 0; i < 3; i++)
          cell_cen[3*c_id2 + i] += i_face_cog[3*f_id + i]*area;
      }
    }

  } /* End of loop on interior faces */

  /* Loop on boundary faces
     --------------------- */

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {

    /* For each cell sharing a border face, we update the numerator
     * of cell_cen and cell_area */

    cs_lnum_t c_id1 = b_face_cells[f_id];

    /* Computation of the area of the face
       (note that c_id1 == -1 may happen for isolated faces,
       which are cleaned afterwards) */

    if (c_id1 > -1) {

      cs_real_t area = cs_math_3_norm(b_face_norm + 3*f_id);

      if (!(cs_glob_mesh_quantities_flag & CS_FACE_NULL_SURFACE) ||
          area > 1.e-20) {
        cell_area[c_id1] += area;

        /* Computation of the numerator */

        for (cs_lnum_t i = 0; i < 3; i++)
          cell_cen[3*c_id1 + i] += b_face_cog[3*f_id + i]*area;
      }

    }

  } /* End of loop on boundary faces */

  /* Loop on cells to finalize the computation of center of gravity
     -------------------------------------------------------------- */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    for (cs_lnum_t i = 0; i < 3; i++)
      cell_cen[c_id*3 + i] /= cell_area[c_id];

  }

  /* Free memory */

  BFT_FREE(cell_area);
}

/*----------------------------------------------------------------------------
 * Compute cell volumes.
 *
 * The corresponding array is allocated by this function, and it is the
 * caller's responsability to free it when they are no longer needed.
 *
 * parameters:
 *   mesh     <-- pointer to a cs_mesh_t structure
 *
 * return:
 *   pointer to newly allocated cell volumes array
 *----------------------------------------------------------------------------*/

cs_real_t *
cs_mesh_quantities_cell_volume(const cs_mesh_t  *mesh)
{
  cs_real_t *cell_vol;
  BFT_MALLOC(cell_vol, mesh->n_cells_with_ghosts, cs_real_t);

  cs_real_3_t *cell_cen;
  BFT_MALLOC(cell_cen, mesh->n_cells_with_ghosts, cs_real_3_t);

  cs_real_t  *i_face_cog = nullptr, *i_face_normal = nullptr;
  cs_real_t  *b_face_cog = nullptr, *b_face_normal = nullptr;

  cs_mesh_quantities_i_faces(mesh, &i_face_cog, &i_face_normal);
  cs_mesh_quantities_b_faces(mesh, &b_face_cog, &b_face_normal);

  _compute_cell_quantities(mesh,
                           (const cs_real_3_t *)i_face_normal,
                           (const cs_real_3_t *)i_face_cog,
                           (const cs_real_3_t *)b_face_normal,
                           (const cs_real_3_t *)b_face_cog,
                           (cs_real_3_t *)cell_cen,
                           cell_vol);

  BFT_FREE(cell_cen);
  BFT_FREE(b_face_normal);
  BFT_FREE(b_face_cog);
  BFT_FREE(i_face_normal);
  BFT_FREE(i_face_cog);

  return cell_vol;
}

/*----------------------------------------------------------------------------
 * Check that no negative volumes are present, and exit on error otherwise.
 *
 * parameters:
 *   mesh            <-- pointer to mesh structure
 *   mesh_quantities <-- pointer to mesh quantities structure
 *   allow_error     <-- 1 if errors are allowed, 0 otherwise
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_check_vol(const cs_mesh_t             *mesh,
                             const cs_mesh_quantities_t  *mesh_quantities,
                             int                          allow_error)
{
  cs_lnum_t  cell_id;

  cs_gnum_t  error_count = 0;

  for (cell_id = 0; cell_id < mesh->n_cells; cell_id++) {
    if (mesh_quantities->cell_vol[cell_id] < 0.0)
      error_count += 1;
  }

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    cs_gnum_t tot_error_count = 0;
    MPI_Allreduce(&error_count, &tot_error_count, 1, CS_MPI_GNUM, MPI_SUM,
                  cs_glob_mpi_comm);
    error_count = tot_error_count;
  }
#endif

  /* Exit with error */

  if (error_count > 0) {
    const char fmt[]
      = N_("  %llu cells have a Negative volume.\n"
           " Run mesh quality check for post-processing output.\n"
           " In case of mesh joining, this may be due to overly "
           " agressive joining parameters.");

    if (allow_error) {
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(_(fmt), (unsigned long long)error_count);
      bft_printf("\n\n");
    }
    else
      bft_error(__FILE__, __LINE__, 0,
                _(fmt), (unsigned long long)error_count);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the bounding box for cells.
 *
 * The corresponding array is allocated by this function, and it is the
 * caller's responsability to free it when they are no longer needed.
 *
 * \param[in]   m          pointer to mesh structure
 * \param[in]   tolerance  addition to local extents of each element:
 *                         extent = base_extent * (1 + tolerance)
 *
 * \return  pointer to newly allocated cell volumes array
 */
/*----------------------------------------------------------------------------*/

cs_real_6_t *
cs_mesh_quantities_cell_extents(const cs_mesh_t  *m,
                                cs_real_t         tolerance)
{
  cs_real_6_t *bbox;

  BFT_MALLOC(bbox, m->n_cells_with_ghosts, cs_real_6_t);

  for (cs_lnum_t i = 0; i < m->n_cells_with_ghosts; i++) {
    bbox[i][0] = HUGE_VAL;
    bbox[i][1] = HUGE_VAL;
    bbox[i][2] = HUGE_VAL;
    bbox[i][3] = -HUGE_VAL;
    bbox[i][4] = -HUGE_VAL;
    bbox[i][5] = -HUGE_VAL;
  }

  const cs_lnum_t n_i_faces = m->n_i_faces;

  for (cs_lnum_t i = 0; i < n_i_faces; i++) {
    cs_lnum_t c_id_0 = m->i_face_cells[i][0];
    cs_lnum_t c_id_1 = m->i_face_cells[i][1];
    cs_lnum_t s_id = m->i_face_vtx_idx[i];
    cs_lnum_t e_id = m->i_face_vtx_idx[i+1];
    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t vtx_id = m->i_face_vtx_lst[j];
      const cs_real_t *coo = m->vtx_coord + vtx_id*3;
      if (c_id_0 > -1) {
        for (cs_lnum_t k = 0; k < 3; k++) {
          bbox[c_id_0][k] = CS_MIN(bbox[c_id_0][k], coo[k]);
          bbox[c_id_0][k+3] = CS_MAX(bbox[c_id_0][k+3], coo[k]);
        }
      }
      if (c_id_1 > -1) {
        for (cs_lnum_t k = 0; k < 3; k++) {
          bbox[c_id_1][k] = CS_MIN(bbox[c_id_1][k], coo[k]);
          bbox[c_id_1][k+3] = CS_MAX(bbox[c_id_1][k+3], coo[k]);
        }
      }
    }
  }

  const cs_lnum_t n_b_faces = m->n_b_faces;

  for (cs_lnum_t i = 0; i < n_b_faces; i++) {
    cs_lnum_t c_id = m->b_face_cells[i];
    cs_lnum_t s_id = m->b_face_vtx_idx[i];
    cs_lnum_t e_id = m->b_face_vtx_idx[i+1];
    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t vtx_id = m->b_face_vtx_lst[j];
      const cs_real_t *coo = m->vtx_coord + vtx_id*3;
      if (c_id > -1) {
        for (cs_lnum_t k = 0; k < 3; k++) {
          bbox[c_id][k] = CS_MIN(bbox[c_id][k], coo[k]);
          bbox[c_id][k+3] = CS_MAX(bbox[c_id][k+3], coo[k]);
        }
      }
    }
  }

  {
    const cs_lnum_t n_cells = m->n_cells;

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      cs_real_t delta[3];

      for (cs_lnum_t i = 0; i < 3; i++)
        delta[i] = (bbox[c_id][3+i] - bbox[c_id][i]) * tolerance;

      for (cs_lnum_t i = 0; i < 3; i++) {
        bbox[c_id][i]   = bbox[c_id][i]   - delta[i];
        bbox[c_id][3+i] = bbox[c_id][3+i] + delta[i];
      }

    }
  }

  return bbox;
}

/*----------------------------------------------------------------------------
 * Return the number of times mesh quantities have been computed.
 *
 * returns:
 *   number of times mesh quantities have been computed
 *----------------------------------------------------------------------------*/

int
cs_mesh_quantities_compute_count(void)
{
  return _n_computations;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Determine local boundary thickness around each vertex.
 *
 * \param[in]   m            pointer to mesh structure
 * \param[in]   mq           pointer to mesh quantities structures.
 * \param[in]   n_passes     number of smoothing passes
 * \param[out]  b_thickness  thickness for each mesh vertex
 *                           (0 at non-boundary vertices)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_quantities_b_thickness_v(const cs_mesh_t             *m,
                                 const cs_mesh_quantities_t  *mq,
                                 int                          n_passes,
                                 cs_real_t                    b_thickness[])
{
  cs_real_t *v_sum = nullptr;
  cs_real_t *f_b_thickness = nullptr;

  BFT_MALLOC(v_sum, m->n_vertices*2, cs_real_t);

  BFT_MALLOC(f_b_thickness, m->n_b_faces*2, cs_real_t);
  _b_thickness(m, mq, f_b_thickness);

  if (n_passes < 1)
    n_passes = 1;

  for (int i = 0; i < n_passes; i++) {

    for (cs_lnum_t j = 0; j < m->n_vertices*2; j++)
      v_sum[j] = 0.;

    for (cs_lnum_t f_id = 0; f_id < m->n_b_faces; f_id++) {
      cs_lnum_t s_id = m->b_face_vtx_idx[f_id];
      cs_lnum_t e_id = m->b_face_vtx_idx[f_id+1];
      const cs_real_t f_s = mq->b_face_surf[f_id];
      for (cs_lnum_t k = s_id; k < e_id; k++) {
        cs_lnum_t v_id = m->b_face_vtx_lst[k];
        v_sum[v_id*2]   += f_s * f_b_thickness[f_id];
        v_sum[v_id*2+1] += f_s;
      }
    }

    if (m->vtx_interfaces != nullptr)
      cs_interface_set_sum(m->vtx_interfaces,
                           m->n_vertices,
                           2,
                           true,
                           CS_REAL_TYPE,
                           v_sum);

    /* Prepare for next smoothing */

    if (i < n_passes-1) {

      for (cs_lnum_t j = 0; j < m->n_b_faces*2; j++)
        f_b_thickness[j] = 0.;

      for (cs_lnum_t f_id = 0; f_id < m->n_b_faces; f_id++) {
        cs_lnum_t s_id = m->b_face_vtx_idx[f_id];
        cs_lnum_t e_id = m->b_face_vtx_idx[f_id+1];
        for (cs_lnum_t k = s_id; k < e_id; k++) {
          cs_lnum_t v_id = m->b_face_vtx_lst[k];
          f_b_thickness[f_id] += v_sum[v_id*2];
          f_b_thickness[f_id + m->n_b_faces] += v_sum[v_id*2 + 1];
        }
      }

      for (cs_lnum_t f_id = 0; f_id < m->n_b_faces; f_id++) {
        if (f_b_thickness[f_id + m->n_b_faces] > 0)
          f_b_thickness[f_id] /= f_b_thickness[f_id + m->n_b_faces];
      }

    }

  }

  BFT_FREE(f_b_thickness);

  for (cs_lnum_t j = 0; j < m->n_vertices; j++) {
    if (v_sum[j*2+1] > 0)
      b_thickness[j] = v_sum[j*2] / v_sum[j*2+1];
    else
      b_thickness[j] = 0;
  }

  BFT_FREE(v_sum);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Determine local boundary thickness around each boundary face.
 *
 * \param[in]   m            pointer to mesh structure
 * \param[in]   mq           pointer to mesh quantities structures.
 * \param[in]   n_passes     number of optional smoothing passes
 * \param[out]  b_thickness  thickness for each mesh boundary face
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_quantities_b_thickness_f(const cs_mesh_t             *m,
                                 const cs_mesh_quantities_t  *mq,
                                 int                          n_passes,
                                 cs_real_t                    b_thickness[])
{
  if (n_passes < 1)
    _b_thickness(m, mq, b_thickness);

  else {

    cs_real_t *v_b_thickness = nullptr;

    BFT_MALLOC(v_b_thickness, m->n_vertices, cs_real_t);

    cs_mesh_quantities_b_thickness_v(m,
                                     mq,
                                     n_passes,
                                     v_b_thickness);

    for (cs_lnum_t f_id = 0; f_id < m->n_b_faces; f_id++) {
      b_thickness[f_id] = 0;
      cs_lnum_t s_id = m->b_face_vtx_idx[f_id];
      cs_lnum_t e_id = m->b_face_vtx_idx[f_id+1];
      for (cs_lnum_t k = s_id; k < e_id; k++) {
        cs_lnum_t v_id = m->b_face_vtx_lst[k];
        b_thickness[f_id] += v_b_thickness[v_id];
      }
      b_thickness[f_id] /= (e_id - s_id);
    }

    BFT_FREE(v_b_thickness);

  }
}

/*----------------------------------------------------------------------------
 * Compute quantities associated to a list of faces (border or internal)
 *
 * parameters:
 *   n_faces         <--  number of faces
 *   vtx_coord       <--  vertex coordinates
 *   face_vtx_idx    <--  "face -> vertices" connectivity index
 *   face_vtx        <--  "face -> vertices" connectivity
 *   face_cog        -->  coordinates of the center of gravity of the faces
 *   face_normal     -->  face surface normals
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_compute_face_quantities(cs_lnum_t        n_faces,
                                           const cs_real_t  vtx_coord[][3],
                                           const cs_lnum_t  face_vtx_idx[],
                                           const cs_lnum_t  face_vtx[],
                                           cs_real_t        face_cog[][3],
                                           cs_real_t        face_normal[][3])
{
  _compute_face_quantities(n_faces,
                           vtx_coord,
                           face_vtx_idx,
                           face_vtx,
                           face_cog,
                           face_normal);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log mesh quantities options to setup file.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_quantities_log_setup(void)
{
  if (cs_glob_mesh_quantities_flag != 0 || _cell_cen_algorithm != 1)
    cs_log_printf(CS_LOG_SETUP, _("\n"
                                  "Mesh quantity computation options\n"
                                  "---------------------------------\n\n"));

  const char *cen_type_name[] = {N_("weighted center of face centers"),
                                 N_("center of mass")};
  cs_log_printf(CS_LOG_SETUP,
                _("  Cell centers: %s\n"),
                _(cen_type_name[_cell_cen_algorithm]));

  if (cs_glob_mesh_quantities_flag != 0) {

    const char *correction_name[] = {"CS_BAD_CELLS_WARPED_CORRECTION",
                                     "CS_BAD_CELLS_REGULARISATION",
                                     "CS_CELL_FACE_CENTER_CORRECTION",
                                     "CS_CELL_CENTER_CORRECTION",
                                     "CS_FACE_DISTANCE_CLIP",
                                     "CS_FACE_RECONSTRUCTION_CLIP",
                                     "CS_CELL_VOLUME_RATIO_CORRECTION",
                                     "CS_FACE_CENTER_REFINE"};

    cs_log_printf(CS_LOG_SETUP,
       ("\n"
        "   Mesh quantity corrections:\n"));

    for (int i = 0; i < 8; i++) {
      if (cs_glob_mesh_quantities_flag & (1 << i))
        cs_log_printf(CS_LOG_SETUP, "      %s\n", correction_name[i]);
    }

  }
}

/*----------------------------------------------------------------------------
 * Dump a cs_mesh_quantities_t structure
 *
 * parameters:
 *   mesh            <-- pointer to a cs_mesh_t structure
 *   mesh_quantities <-- pointer to a cs_mesh_quantities_t structure
 *----------------------------------------------------------------------------*/

void
cs_mesh_quantities_dump(const cs_mesh_t             *mesh,
                        const cs_mesh_quantities_t  *mesh_quantities)
{
  cs_lnum_t  i;

  const cs_lnum_t  n_cells = mesh->n_cells_with_ghosts;
  const cs_lnum_t  n_i_faces = mesh->n_i_faces;
  const cs_lnum_t  n_b_faces = mesh->n_b_faces;

  const cs_real_t  *cell_cen = mesh_quantities->cell_cen;
  const cs_real_t  *cell_vol = mesh_quantities->cell_vol;
  const cs_real_t  *i_fac_norm = mesh_quantities->i_face_normal;
  const cs_real_t  *b_fac_norm = mesh_quantities->b_face_normal;
  const cs_real_t  *i_fac_cog = mesh_quantities->i_face_cog;
  const cs_real_t  *b_fac_cog = mesh_quantities->b_face_cog;
  const cs_real_t  *i_fac_surf = mesh_quantities->i_face_surf;
  const cs_real_t  *b_fac_surf = mesh_quantities->b_face_surf;

  bft_printf("\n\nDUMP OF A MESH QUANTITIES STRUCTURE: %p\n\n",
             (const void *)mesh_quantities);

  if (mesh_quantities == nullptr)
    return;

  /* Cell data */

  bft_printf("\n\n"
             "    ---------------"
             "    Cell quantities"
             "    ---------------\n\n");

  bft_printf("Cell center coordinates:\n");
  for (i = 0; i < n_cells; i++)
    bft_printf("    < %ld >    %.3f    %.3f    %.3f\n", (long)i+1,
               cell_cen[3*i], cell_cen[3*i+1], cell_cen[3*i+2]);

  bft_printf("\nCell volume:\n");
  for (i = 0; i < n_cells; i++)
    bft_printf("    < %ld >    %.3f\n", (long)i+1, cell_vol[i]);

  /* Internal faces data */

  bft_printf("\n\n"
             "    ------------------------"
             "    Interior face quantities"
             "    ------------------------\n\n");

  bft_printf("\nInterior face normals\n");
  for (i = 0; i < n_i_faces; i++)
    bft_printf("    < %ld >    %.3f    %.3f    %.3f\n", (long)i+1,
               i_fac_norm[3*i], i_fac_norm[3*i+1], i_fac_norm[3*i+2]);

  bft_printf("\nInterior face centers\n");
  for (i = 0; i < n_i_faces; i++)
    bft_printf("    < %ld >    %.3f    %.3f    %.3f\n", (long)i+1,
               i_fac_cog[3*i], i_fac_cog[3*i+1], i_fac_cog[3*i+2]);

  bft_printf("\nInterior face surfaces\n");
  for (i = 0; i < n_i_faces; i++)
    bft_printf("    < %ld >    %.3f\n", (long)i+1, i_fac_surf[i]);

  /* Border faces data */

  bft_printf("\n\n"
             "    ------------------------"
             "    Boundary face quantities"
             "    ------------------------\n\n");

  bft_printf("\nBoundary face normals\n");
  for (i = 0; i < n_b_faces; i++)
    bft_printf("    < %ld >    %.3f    %.3f    %.3f\n", (long)i+1,
               b_fac_norm[3*i], b_fac_norm[3*i+1], b_fac_norm[3*i+2]);

  bft_printf("\nBoundary faces centers\n");
  for (i = 0; i < n_b_faces; i++)
    bft_printf("    < %ld >    %.3f    %.3f    %.3f\n", (long)i+1,
               b_fac_cog[3*i], b_fac_cog[3*i+1], b_fac_cog[3*i+2]);

  bft_printf("\nBoundary face surfaces\n");
  for (i = 0; i < n_b_faces; i++)
    bft_printf("    < %ld >    %.3f\n", (long)i+1, b_fac_surf[i]);

  bft_printf("\n\nEND OF DUMP OF MESH QUANTITIES STRUCTURE\n\n");
  bft_printf_flush();
}

/*----------------------------------------------------------------------------*/

#if 0 /* Test if face orientation is OK */

  cs_lnum_t  i, fac_id, cell_id;
  cs_real_t  cogfac[3];
  cs_real_t  cogcel[3];
  cs_real_t  normal[3];
  cs_real_t  pscal;

  for (fac_id = 0; fac_id < mesh->n_b_faces; fac_id++) {

    cell_id = mesh->b_face_cells[fac_id];
    pscal = 0;

    for (i = 0; i < 3; i++) {
      cogcel[i]  = cs_glob_mesh_quantities->cell_cen[cell_id*3 + i];
      cogfac[i]  = cs_glob_mesh_quantities->b_face_cog[fac_id*3 + i];
      normal[i] = cs_glob_mesh_quantities->b_face_normal[fac_id*3 + i];
      pscal += normal[i] * (cogfac[i] - cogcel[i]);
    }

    if (pscal < 0.0)
      printf("num_fac_brd = %d, num_cel = %d, pscal = %f\n",
             fac_id + 1, cell_id + 1, pscal);
  }

#endif

/*----------------------------------------------------------------------------*/

END_C_DECLS
