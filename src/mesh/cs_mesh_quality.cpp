/*============================================================================
 * Compute several mesh quality criteria.
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_dispatch.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"

#include "alge/cs_blas.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_adjacencies.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_post.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "mesh/cs_mesh_quality.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define CS_MESH_QUALITY_N_SUBS  10

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute the minimum and the maximum of a vector (locally).
 *
 * parameters:
 *   n_vals    <-- local number of elements
 *   var       <-- pointer to vector
 *   min       --> minimum
 *   max       --> maximum
 *----------------------------------------------------------------------------*/

static void
_compute_local_minmax(cs_lnum_t        n_vals,
                      const cs_real_t  var[],
                      cs_real_t       *min,
                      cs_real_t       *max)
{
  cs_lnum_t  i;
  cs_real_t  _min = DBL_MAX, _max = -DBL_MAX;

  for (i = 0; i < n_vals; i++) {
    _min = cs::min(_min, var[i]);
    _max = cs::max(_max, var[i]);
  }

  *min = _min;
  *max = _max;
}

/*----------------------------------------------------------------------------
 * Display the distribution of values of a real vector.
 *
 * parameters:
 *   n_steps <-- number of histogram steps
 *   var_min <-- minimum variable value
 *   var_max <-- maximum variable value
 *   count   <-> count for each histogram slice (size: n_steps)
 *               local values in, global values out
 *----------------------------------------------------------------------------*/

static void
_display_histograms(int        n_steps,
                    cs_real_t  var_min,
                    cs_real_t  var_max,
                    cs_gnum_t  count[])
{
#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    cs_gnum_t _g_count[CS_MESH_QUALITY_N_SUBS];
    cs_gnum_t *g_count = _g_count;

    if (n_steps > CS_MESH_QUALITY_N_SUBS)
      CS_MALLOC(g_count, n_steps, cs_gnum_t);

    MPI_Allreduce(count, g_count, n_steps, CS_MPI_GNUM, MPI_SUM,
                  cs_glob_mpi_comm);

    for (int i = 0; i < n_steps; i++)
      count[i] = g_count[i];

    if (n_steps > CS_MESH_QUALITY_N_SUBS)
      CS_FREE(g_count);

  }

#endif

  /* Print base min, max, and increment */

  bft_printf(_("    minimum value =         %10.5e\n"), (double)var_min);
  bft_printf(_("    maximum value =         %10.5e\n\n"), (double)var_max);

  double var_step = cs::abs(var_max - var_min) / n_steps;

  if (cs::abs(var_max - var_min) > 0.) {

    /* Number of elements in each subdivision */

    int i, j;
    for (i = 0, j = 1; i < n_steps - 1; i++, j++)
      bft_printf("    %3d : [ %10.5e ; %10.5e [ = %10llu\n",
                 i+1, var_min + i*var_step, var_min + j*var_step,
                 (unsigned long long)(count[i]));

    bft_printf("    %3d : [ %10.5e ; %10.5e ] = %10llu\n",
               n_steps,
               var_min + (n_steps - 1)*var_step,
               var_max,
               (unsigned long long)(count[n_steps - 1]));

  }
}

/*----------------------------------------------------------------------------
 * Display the distribution of values of a real vector on cells or
 * boundary faces.
 *
 * parameters:
 *   n_vals     <-- number of values
 *   var        <-- pointer to vector (size: n_vals)
 *----------------------------------------------------------------------------*/

static void
_histogram(cs_lnum_t        n_vals,
           const cs_real_t  var[])
{
  cs_gnum_t count[CS_MESH_QUALITY_N_SUBS];
  const int  n_steps = CS_MESH_QUALITY_N_SUBS;

  assert (sizeof(double) == sizeof(cs_real_t));

  /* Compute global min and max */

  cs_real_t  max, min;
  _compute_local_minmax(n_vals, var, &min, &max);

  /* Default initialization */

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {
    cs_real_t mm[2] = {min, -max}, mmg[2];
    MPI_Allreduce(mm, mmg, 2, CS_MPI_REAL, MPI_MIN,
                  cs_glob_mpi_comm);
    min = mmg[0]; max = -mmg[1];
  }

#endif

  /* Define axis subdivisions */

  for (int j = 0; j < n_steps; j++)
    count[j] = 0;

  if (cs::abs(max - min) > 0.) {

    cs_real_t step = cs::abs(max - min) / n_steps;

    /* Loop on values */

    for (int i = 0; i < n_vals; i++) {

      /* Associated subdivision */

      int j, k;
      for (j = 0, k = 1; k < n_steps; j++, k++) {
        if (var[i] < min + k*step)
          break;
      }
      count[j] += 1;

    }

  }

  _display_histograms(n_steps, min, max, count);
}

/*----------------------------------------------------------------------------
 * Display the distribution of values of a real vector on interior faces.
 *
 * parameters:
 *   mesh <-- pointer to mesh structure
 *   var  <-- pointer to vector (size: mesh->n_i_faces)
 *----------------------------------------------------------------------------*/

static void
_int_face_histogram(const cs_mesh_t  *mesh,
                    const cs_real_t   var[])
{
  cs_gnum_t count[CS_MESH_QUALITY_N_SUBS];
  const int  n_steps = CS_MESH_QUALITY_N_SUBS;

  assert (sizeof (double) == sizeof (cs_real_t));

  /* Compute global min and max */

  cs_real_t  max, min;
  _compute_local_minmax(mesh->n_i_faces, var, &min, &max);

  /* Default initialization */

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {
    cs_real_t mm[2] = {min, -max}, mmg[2];
    MPI_Allreduce(mm, mmg, 2, CS_MPI_REAL, MPI_MIN,
                  cs_glob_mpi_comm);
    min = mmg[0]; max = -mmg[1];
  }

#endif

  /* Define axis subdivisions */

  for (int j = 0; j < n_steps; j++)
    count[j] = 0;

  if (cs::abs(max - min) > 0.) {

    cs_lnum_t n_i_faces = mesh->n_i_faces;
    cs_real_t step = cs::abs(max - min) / n_steps;
    const cs_lnum_2_t *i_face_cells = mesh->i_face_cells;

    /* Loop on faces */

    for (cs_lnum_t i = 0; i < n_i_faces; i++) {

      if (i_face_cells[i][0] < mesh->n_cells) {

        /* Associated subdivision */

        int j, k;
        for (j = 0, k = 1; k < n_steps; j++, k++) {
          if (var[i] < min + k*step)
            break;
        }
        count[j] += 1;

      }

    }

  }

  _display_histograms(n_steps, min, max, count);
}

/*----------------------------------------------------------------------------
 * Compute weighting coefficient for internal faces.
 *
 * parameters:
 *   mesh             <-- pointer to mesh structure.
 *   mesh_quantities  <-- pointer to mesh quantities structures.
 *   ctx              <-- Reference to dispatch context
 *   weighting        <-> array for weigthing coefficient.
 *----------------------------------------------------------------------------*/

static void
_compute_weighting(const cs_mesh_t             *mesh,
                   const cs_mesh_quantities_t  *mesh_quantities,
                   cs_dispatch_context         &ctx,
                   cs_real_t                    weighting[])
{
  const cs_lnum_2_t *i_face_cells = mesh->i_face_cells;

  const cs_real_3_t *cell_cen = mesh_quantities->cell_cen;
  const cs_real_3_t *i_face_cog = mesh_quantities->i_face_cog;
  const cs_nreal_3_t *i_face_u_normal = mesh_quantities->i_face_u_normal;

  /* Loop on internal faces */

  ctx.parallel_for(mesh->n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {

    /* Get local number of the cells in contact with the face */

    cs_lnum_t cell1 = i_face_cells[face_id][0];
    cs_lnum_t cell2 = i_face_cells[face_id][1];

    /* Get information on mesh quantities */

    cs_real_t  v0[3], v1[3], v2[3];

    /* Compute weighting coefficient with two approaches. Keep the max value. */

    for (cs_lnum_t i = 0; i < 3; i++) {
      v0[i] = cell_cen[cell2][i] - cell_cen[cell1][i];
      v1[i] = i_face_cog[face_id][i] - cell_cen[cell1][i];
      v2[i] = cell_cen[cell2][i] - i_face_cog[face_id][i];
    }

    cs_real_t coef0 = cs_math_3_dot_product(v0, i_face_u_normal[face_id]);
    cs_real_t coef1 = cs_math_3_dot_product(v1, i_face_u_normal[face_id])/coef0;
    cs_real_t coef2 = cs_math_3_dot_product(v2, i_face_u_normal[face_id])/coef0;

    weighting[face_id] = cs::max(coef1, coef2);

  }); /* End of loop on faces */
}

/*----------------------------------------------------------------------------
 * Compute center offsetting coefficient for internal faces.
 *
 * parameters:
 *   mesh             <-- pointer to mesh structure.
 *   ma               <-- pointer to mesh adjacencies structure
 *   mesh_quantities  <-- pointer to mesh quantities structures.
 *   ctx              <-- Reference to dispatch context
 *   offsetting       <-> array for offsetting coefficient
 *----------------------------------------------------------------------------*/

static void
_compute_offsetting(const cs_mesh_t               *mesh,
                    const cs_mesh_adjacencies_t   *ma,
                    const cs_mesh_quantities_t    *mesh_quantities,
                    cs_dispatch_context          &ctx,
                    cs_real_t                     offsetting[])
{
  const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
  const cs_lnum_t *c2f = ma->cell_i_faces;
  if (c2f == nullptr) {
    cs_mesh_adjacencies_update_cell_i_faces();
    c2f = ma->cell_i_faces;
  }

  const cs_real_t *cell_vol = mesh_quantities->cell_vol;
  const cs_real_t *i_face_surf = mesh_quantities->i_face_surf;
  const cs_real_3_t *dofij = mesh_quantities->dofij;

  ctx.parallel_for(mesh->n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {

    cs_lnum_t s_id = c2c_idx[ii];
    cs_lnum_t e_id = c2c_idx[ii+1];

    /* Loop on internal faces */

    cs_real_t offset = 1.;
    cs_real_t vol = cell_vol[ii];

    if (vol > 1e-30) {
      for (cs_lnum_t i = s_id; i < e_id; i++) {
        const cs_lnum_t face_id = c2f[i];

        double of_s =   cs_math_3_norm(dofij[face_id])
                      * i_face_surf[face_id];

        offset = cs::min(offset, 1. - pow(of_s/vol, 1./3.));
      }
    }

    offsetting[ii] = offset;

  });
}

/*----------------------------------------------------------------------------
 * Compute angle between face normal and segment based on centers of the
 * adjacent cells. Evaluates a level of non-orthogonality.
 *
 * parameters:
 *   mesh             <-- pointer to mesh structure.
 *   mesh_quantities  <-- pointer to mesh quantities structures.
 *   ctx              <-- Reference to dispatch context
 *   i_face_ortho     <-> array for internal faces.
 *   b_face_ortho     <-> array for border faces.
 *----------------------------------------------------------------------------*/

static void
_compute_orthogonality(const cs_mesh_t             *mesh,
                       const cs_mesh_quantities_t  *mesh_quantities,
                       cs_dispatch_context         &ctx,
                       cs_real_t                    i_face_ortho[],
                       cs_real_t                    b_face_ortho[])
{
  const cs_lnum_2_t *i_face_cells = mesh->i_face_cells;
  const cs_lnum_t *b_face_cells = mesh->b_face_cells;

  const cs_real_3_t *cell_cen = mesh_quantities->cell_cen;
  const cs_real_3_t *b_face_cog = mesh_quantities->b_face_cog;
  const cs_nreal_3_t *i_face_u_normal = mesh_quantities->i_face_u_normal;
  const cs_nreal_3_t *b_face_u_normal = mesh_quantities->b_face_u_normal;

  const double  rad_to_deg = 180. / acos(-1.);

  /* Loop on internal faces */

  ctx.parallel_for(mesh->n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {

    cs_lnum_t cell1 = i_face_cells[face_id][0];
    cs_lnum_t cell2 = i_face_cells[face_id][1];

    /* Compute angle which evaluates the non-orthogonality. */

    cs_real_t dc[3];
    for (cs_lnum_t i = 0; i < 3; i++)
      dc[i] = cell_cen[cell2][i] - cell_cen[cell1][i];

    cs_real_t cos_alpha =   cs_math_3_dot_product(dc, i_face_u_normal[face_id])
                          / cs_math_3_norm(dc);

    cos_alpha = cs::abs(cos_alpha);
    cos_alpha = cs::min(cos_alpha, 1);

    if (cos_alpha < 1.)
      i_face_ortho[face_id] = acos(cos_alpha) * rad_to_deg;
    else
      i_face_ortho[face_id] = 0.;

  }); /* End of loop on internal faces */

  /* Loop on boundary faces */

  ctx.parallel_for(mesh->n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {

    cs_lnum_t cell1 = b_face_cells[face_id];

    /* Compute alpha: angle wich evaluate the difference with orthogonality. */

    cs_real_t dc[3];
    for (cs_lnum_t i = 0; i < 3; i++)
      dc[i] = b_face_cog[face_id][i] - cell_cen[cell1][i];

    cs_real_t cos_alpha =   cs_math_3_dot_product(dc, b_face_u_normal[face_id])
                          / cs_math_3_norm(dc);

    cos_alpha = cs::abs(cos_alpha);
    cos_alpha = cs::min(cos_alpha, 1);

    if (cos_alpha < 1.)
      b_face_ortho[face_id] = acos(cos_alpha) * rad_to_deg;
    else
      b_face_ortho[face_id] = 0.;

  }); /* End of loop on boundary faces */
}

/*----------------------------------------------------------------------------
 * Evaluate face warping angle.
 *
 * parameters:
 *   idx_start       <-- first vertex index
 *   idx_end         <-- last vertex index
 *   face_vertex_id  <-- face -> vertices connectivity
 *   face_u_normal   <-- face unit normal
 *   vertex_coords   <-- vertices coordinates
 *   face_warping    --> face warping angle
 *----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE static void
_get_face_warping(cs_lnum_t          idx_start,
                  cs_lnum_t          idx_end,
                  const cs_nreal_t   face_u_normal[],
                  const cs_lnum_t    face_vertex_id[],
                  const cs_real_t    vertex_coords[],
                  double            &face_warping)
{
  double  cos_alpha = 0.;

  const double  rad_to_deg = 180. / cs_math_pi;

  cs_lnum_t n_vtx = idx_end - idx_start;

  /* Loop on edges */

  for (cs_lnum_t idx = 0; idx < n_vtx; idx++) {

    cs_lnum_t vertex_id1 = face_vertex_id[idx_start + idx];
    cs_lnum_t vertex_id2 = face_vertex_id[idx_start + (idx+1)%n_vtx];

    /* Get vertex coordinates */

    cs_real_t  vect[3];
    for (cs_lnum_t i = 0; i < 3; i++)
      vect[i] =  vertex_coords[vertex_id2*3 + i]
               - vertex_coords[vertex_id1*3 + i];

    cs_real_t edge_cos_alpha
      =   cs_math_3_dot_product(vect, face_u_normal)
        / cs_math_3_norm(vect);

    edge_cos_alpha = cs::abs(edge_cos_alpha);
    cos_alpha = cs::max(cos_alpha, edge_cos_alpha);

  }

  cos_alpha = cs::min(cos_alpha, 1.);

  face_warping = 90. - acos(cos_alpha) * rad_to_deg;
}

/*----------------------------------------------------------------------------
 * Compute cellwise the warping error
 *  Id = 1/|c| \sum_(f \in F_c) x_f \otimes \vect{f}
 * Froebinus norm is used to get a scalar-valued quantity
 *
 * parameters:
 *   mesh             <-- pointer to mesh structure.
 *   ma               <-- pointer to mesh adjacencies structure
 *   mesh_quantities  <-- pointer to mesh quantities structures.
 *   ctx              <-- Reference to dispatch context
 *   warp_error       <-> array of values to compute
 *----------------------------------------------------------------------------*/

static void
_compute_warp_error(const cs_mesh_t              *mesh,
                    const cs_mesh_adjacencies_t  *ma,
                    const cs_mesh_quantities_t   *mesh_quantities,
                    cs_dispatch_context          &ctx,
                    cs_real_t                     warp_error[])
{
  const cs_lnum_t *c2c_idx = ma->cell_cells_idx;
  const cs_lnum_t *c2f = ma->cell_i_faces;
  if (c2f == nullptr) {
    cs_mesh_adjacencies_update_cell_i_faces();
    c2f = ma->cell_i_faces;
  }
  const short int *c2f_sgn = ma->cell_i_faces_sgn;

  const cs_lnum_t *restrict c2b_idx = ma->cell_b_faces_idx;
  const cs_lnum_t *restrict c2b = ma->cell_b_faces;

  const cs_real_3_t *cell_cen = mesh_quantities->cell_cen;
  const cs_real_t  *cell_vol = mesh_quantities->cell_vol;
  const cs_real_3_t *i_face_cog = mesh_quantities->i_face_cog;
  const cs_real_3_t *b_face_cog = mesh_quantities->b_face_cog;
  const cs_real_t *i_face_surf = mesh_quantities->i_face_surf;
  const cs_real_t *b_face_surf = mesh_quantities->b_face_surf;
  const cs_nreal_3_t *i_face_u_normal = mesh_quantities->i_face_u_normal;
  const cs_nreal_3_t *b_face_u_normal = mesh_quantities->b_face_u_normal;

  ctx.parallel_for(mesh->n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    cs_real_t vol = cell_vol[c_id];
    if (vol < 1e-30) {
      warp_error[c_id] = 1e-30;
      return;
    }

    const cs_real_t  invvol_c = 1/vol;
    const cs_real_t  *xc = cell_cen[c_id];

    cs_real_t tens[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

    /* Loop on interior faces */

    cs_lnum_t s_id = c2c_idx[c_id];
    cs_lnum_t e_id = c2c_idx[c_id+1];

    for (cs_lnum_t i = s_id; i < e_id; i++) {

      cs_lnum_t  f_id = c2f[i];
      int  sgn = c2f_sgn[i];

      const cs_real_t *xf = i_face_cog[f_id];
      cs_real_t sgn_surf[3];
      for (cs_lnum_t j = 0; j < 3; j++)
        sgn_surf[j] = sgn * i_face_u_normal[f_id][j]*i_face_surf[f_id];

      for (int ki = 0; ki < 3; ki++)
        for (int kj = 0; kj < 3; kj++)
          tens[ki][kj] += (xf[ki] - xc[ki]) * sgn_surf[kj];

    } // Loop on interior faces

    /* Loop on boundary faces */

    s_id = c2b_idx[c_id];
    e_id = c2b_idx[c_id+1];

    for (cs_lnum_t i = s_id; i < e_id; i++) {

      cs_lnum_t  f_id = c2b[i];

      const cs_real_t *xf = b_face_cog[f_id];
      cs_real_t surf[3];
      for (cs_lnum_t j = 0; j < 3; j++)
        surf[j] = b_face_u_normal[f_id][j]*b_face_surf[f_id];

      for (int ki = 0; ki < 3; ki++)
        for (int kj = 0; kj < 3; kj++)
          tens[ki][kj] += (xf[ki] - xc[ki]) * surf[kj];

    } // Loop on boundary faces

    for (int ki = 0; ki < 3; ki++)
      for (int kj = 0; kj < 3; kj++)
        tens[ki][kj] *= invvol_c;

    warp_error[c_id] = cs::abs(cs_math_33_determinant(tens) - 1.);

  }); // Loop on cells
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Evaluate face warping angle for internal and border faces..
 *
 * parameters:
 *   mesh             <-- pointer to a cs_mesh_t structure
 *   i_face_u_normal  <-- internal face unit normal
 *   b_face_u_normal  <-- boundary face unit normal
 *   i_face_warping   --> face warping angle for internal faces
 *   b_face_warping   --> face warping angle for border faces
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
cs_mesh_quality_compute_warping(const cs_mesh_t      *mesh,
                                const cs_nreal_3_t    i_face_u_normal[],
                                const cs_nreal_3_t    b_face_u_normal[],
                                cs_real_t             i_face_warping[],
                                cs_real_t             b_face_warping[])
{
  const cs_lnum_t  *i_face_vtx_idx = mesh->i_face_vtx_idx;
  const cs_lnum_t  *i_face_vtx = mesh->i_face_vtx_lst;

  const cs_real_t *vtx_coord = mesh->vtx_coord;

  assert(mesh->dim == 3);

  cs_dispatch_context ctx;

  if (cs_check_device_ptr(mesh->i_face_vtx_idx) == CS_ALLOC_HOST)
    ctx.set_use_gpu(false);

  /* Compute warping for internal faces */
  /*------------------------------------*/

  ctx.parallel_for(mesh->n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {

    /* Evaluate warping for each edge of face. Keep the max. */

    cs_lnum_t idx_start = i_face_vtx_idx[face_id];
    cs_lnum_t idx_end = i_face_vtx_idx[face_id + 1];

    _get_face_warping(idx_start,
                      idx_end,
                      i_face_u_normal[face_id],
                      i_face_vtx,
                      vtx_coord,
                      i_face_warping[face_id]);

  }); /* End of loop on internal faces */

  /* Compute warping for border faces */
  /*----------------------------------*/

  cs_mesh_quality_compute_b_face_warping(mesh,
                                         b_face_u_normal,
                                         b_face_warping);

  ctx.wait();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate face warping angle for boundary faces..
 *
 * \param[in]   mesh              pointer to a cs_mesh_t structure
 * \param[in]   b_face_u_normal   boundary face unit normal
 * \param[out]  b_face_warping    face warping angle for boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_quality_compute_b_face_warping(const cs_mesh_t     *mesh,
                                       const cs_nreal_3_t   b_face_u_normal[],
                                       cs_real_t            b_face_warping[])
{
  const cs_lnum_t  *b_face_vtx_idx = mesh->b_face_vtx_idx;
  const cs_lnum_t  *b_face_vtx = mesh->b_face_vtx_lst;
  const cs_real_t *vtx_coord = mesh->vtx_coord;

  assert(mesh->dim == 3);

  cs_dispatch_context ctx;

  if (cs_check_device_ptr(mesh->i_face_vtx_idx) == CS_ALLOC_HOST)
    ctx.set_use_gpu(false);

  ctx.parallel_for(mesh->n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {

    /* Evaluate warping for each edge */

    cs_lnum_t idx_start = b_face_vtx_idx[face_id];
    cs_lnum_t idx_end = b_face_vtx_idx[face_id + 1];

    _get_face_warping(idx_start,
                      idx_end,
                      b_face_u_normal[face_id],
                      b_face_vtx,
                      vtx_coord,
                      b_face_warping[face_id]);

  });

  ctx.wait();
}

/*----------------------------------------------------------------------------
 * Compute mesh quality indicators
 *
 * parameters:
 *   mesh             <-- pointer to mesh structure.
 *   mesh_quantities  <-- pointer to mesh quantities structure.
 *----------------------------------------------------------------------------*/

void
cs_mesh_quality(const cs_mesh_t             *mesh,
                const cs_mesh_quantities_t  *mesh_quantities)
{
  cs_lnum_t  i;

  bool  compute_volume = true;
  bool  compute_weighting = true;
  bool  compute_orthogonality = true;
  bool  compute_warping = true;
  bool  compute_thickness = true;
  bool  compute_warp_error = true;

  double  *working_array = nullptr;

  const cs_lnum_t  n_i_faces = mesh->n_i_faces;
  const cs_lnum_t  n_b_faces = mesh->n_b_faces;
  const cs_lnum_t  n_cells = mesh->n_cells;
  const cs_lnum_t  n_cells_wghosts = mesh->n_cells_with_ghosts;

  const cs_time_step_t *ts = cs_glob_time_step;

  cs_dispatch_context ctx;

  if (cs_check_device_ptr(mesh->i_face_vtx_idx) == CS_ALLOC_HOST)
    ctx.set_use_gpu(false);

  /* Check input data */

  assert(mesh_quantities->i_face_normal != nullptr || mesh->n_i_faces == 0);
  assert(mesh_quantities->i_face_cog != nullptr || mesh->n_i_faces == 0);
  assert(mesh_quantities->cell_cen != nullptr);
  assert(mesh_quantities->cell_vol != nullptr);

  /* Assume resulting fields will be exported to postprocessing meshes */

  /* TODO
     For the moment, we export the mesh at this stage; this should be moved
     once PSTEMA has been moved from CALTRI to an earlier step. */

  cs_post_activate_writer(0, true);

  bool post_fields = false;

  if (cs_post_mesh_find_next_with_cat_id(CS_POST_MESH_SURFACE, 0) != 0)
    post_fields = true;

  if (cs_post_mesh_find_next_with_cat_id(CS_POST_MESH_BOUNDARY, 0) != 0)
    post_fields = true;

  if (post_fields)
    cs_post_write_meshes(ts);

  /* Evaluate mesh quality criteria */
  /*--------------------------------*/

  /*--------------*/
  /* Face warping */
  /*--------------*/

  if (compute_warping == true) {

    cs_real_t *i_face_warping = nullptr, *b_face_warping = nullptr;

    CS_MALLOC_HD(working_array, n_i_faces + n_b_faces, cs_real_t, cs_alloc_mode);

    i_face_warping = working_array;
    b_face_warping = working_array + n_i_faces;

    cs_mesh_quality_compute_warping(mesh,
                                    mesh_quantities->i_face_u_normal,
                                    mesh_quantities->b_face_u_normal,
                                    i_face_warping,
                                    b_face_warping);

    /* Display histograms */

    if (mesh->n_g_i_faces > 0) {
      bft_printf(_("\n  Histogram of the interior faces warping:\n\n"));
      _int_face_histogram(mesh, i_face_warping);
    }

    if (mesh->n_g_b_faces > 0) {
      bft_printf(_("\n  Histogram of the boundary faces warping:\n\n"));
      _histogram(n_b_faces, b_face_warping);
    }

    /* Post processing */

    int mesh_id = cs_post_mesh_find_next_with_cat_id(CS_POST_MESH_SURFACE, 0);

    while (mesh_id != 0) {
      cs_post_write_var(mesh_id,
                        CS_POST_WRITER_ALL_ASSOCIATED,
                        "Face_Warp",
                        1,
                        false,
                        true,
                        CS_POST_TYPE_cs_real_t,
                        nullptr,
                        i_face_warping,
                        b_face_warping,
                        ts);

      mesh_id = cs_post_mesh_find_next_with_cat_id(CS_POST_MESH_SURFACE,
                                                   mesh_id);
    }

    mesh_id = cs_post_mesh_find_next_with_cat_id(CS_POST_MESH_BOUNDARY, 0);

    while (mesh_id != 0) {
      cs_post_write_var(mesh_id,
                        CS_POST_WRITER_ALL_ASSOCIATED,
                        "Face_Warp",
                        1,
                        false,
                        true,
                        CS_POST_TYPE_cs_real_t,
                        nullptr,
                        nullptr,
                        b_face_warping,
                        ts);

      mesh_id = cs_post_mesh_find_next_with_cat_id(CS_POST_MESH_BOUNDARY,
                                                   mesh_id);
    }

    CS_FREE(working_array);

  } /* End of face warping treatment */

  /*----------------------------------------------*/
  /* Weighting and center offsetting coefficients */
  /*----------------------------------------------*/

  if (compute_weighting == true) {

    /* Only defined on internal faces */

    CS_MALLOC_HD(working_array, n_i_faces + n_cells_wghosts, cs_real_t, cs_alloc_mode);

    cs_real_t *weighting = working_array;
    cs_real_t *offsetting = working_array + n_i_faces;

    _compute_weighting(mesh,
                       mesh_quantities,
                       ctx,
                       weighting);

    _compute_offsetting(mesh,
                        cs_glob_mesh_adjacencies,
                        mesh_quantities,
                        ctx,
                        offsetting);

    ctx.wait();

    /* Display histograms */

    if (mesh->n_g_i_faces > 0) {
      bft_printf(_("\n  Histogram of the interior faces "
                   "weighting coefficient:\n\n"));
      _int_face_histogram(mesh, weighting);
    }

    bft_printf(_("\n  Histogram of the cells "
                 "off-centering coefficient:\n\n"));
    _histogram(n_cells, offsetting);

    /* Post processing */

    int mesh_id = cs_post_mesh_find_next_with_cat_id(CS_POST_MESH_SURFACE, 0);

    while (mesh_id != 0) {
      cs_post_write_var(mesh_id,
                        CS_POST_WRITER_ALL_ASSOCIATED,
                        "Weighting coefficient",
                        1,
                        false,
                        true,
                        CS_POST_TYPE_cs_real_t,
                        nullptr,
                        weighting,
                        nullptr,
                        ts);

      mesh_id = cs_post_mesh_find_next_with_cat_id(CS_POST_MESH_SURFACE, mesh_id);
    }

    mesh_id = cs_post_mesh_find_next_with_cat_id(CS_POST_MESH_VOLUME, 0);

    while (mesh_id != 0) {
      cs_post_write_var(mesh_id,
                        CS_POST_WRITER_ALL_ASSOCIATED,
                        "Offset",
                        1,
                        false,
                        true,
                        CS_POST_TYPE_cs_real_t,
                        offsetting,
                        nullptr,
                        nullptr,
                        ts);

      mesh_id = cs_post_mesh_find_next_with_cat_id(CS_POST_MESH_VOLUME, mesh_id);
    }

    CS_FREE(working_array);

  } /* End of off-setting and weighting treatment */

  /*---------------------*/
  /* Angle orthogonality */
  /*---------------------*/

  if (compute_orthogonality == true) {

    CS_MALLOC_HD(working_array, n_i_faces + n_b_faces, cs_real_t, cs_alloc_mode);

    for (i = 0; i < n_i_faces + n_b_faces; i++)
      working_array[i] = 0.;

    cs_real_t *i_face_ortho = working_array;
    cs_real_t *b_face_ortho = working_array + n_i_faces;

    _compute_orthogonality(mesh,
                           mesh_quantities,
                           ctx,
                           i_face_ortho,
                           b_face_ortho);

    ctx.wait();

    /* Display histograms */

    if (mesh->n_g_i_faces > 0) {
      bft_printf(_("\n  Histogram of the interior faces "
                   "non-orthogonality coefficient (in degrees):\n\n"));
      _int_face_histogram(mesh, i_face_ortho);
    }

    if (mesh->n_g_b_faces > 0) {
      bft_printf(_("\n  Histogram of the boundary faces "
                   "non-orthogonality coefficient (in degrees):\n\n"));
      _histogram(n_b_faces, b_face_ortho);
    }

    /* Post processing */

    int mesh_id = cs_post_mesh_find_next_with_cat_id(CS_POST_MESH_SURFACE, 0);

    while (mesh_id != 0) {
      cs_post_write_var(mesh_id,
                        CS_POST_WRITER_ALL_ASSOCIATED,
                        "Non_Ortho",
                        1,
                        false,
                        true,
                        CS_POST_TYPE_cs_real_t,
                        nullptr,
                        i_face_ortho,
                        b_face_ortho,
                        ts);

      mesh_id = cs_post_mesh_find_next_with_cat_id(CS_POST_MESH_SURFACE,
                                                   mesh_id);
    }

    mesh_id = cs_post_mesh_find_next_with_cat_id(CS_POST_MESH_BOUNDARY, 0);

    while (mesh_id != 0) {
      cs_post_write_var(mesh_id,
                        CS_POST_WRITER_ALL_ASSOCIATED,
                        "Non_Ortho",
                        1,
                        false,
                        true,
                        CS_POST_TYPE_cs_real_t,
                        nullptr,
                        nullptr,
                        b_face_ortho,
                        ts);

      mesh_id = cs_post_mesh_find_next_with_cat_id(CS_POST_MESH_BOUNDARY,
                                                   mesh_id);
    }

    CS_FREE(working_array);

  } /* End of non-orthogonality treatment */

  /*-------------*/
  /* Cell volume */
  /*-------------*/

  if (compute_volume == true) {

    /* Display histograms */

    bft_printf(_("\n  Histogram of cell volumes:\n\n"));
    _histogram(n_cells, mesh_quantities->cell_vol);

    /* Post processing */

    int mesh_id = cs_post_mesh_find_next_with_cat_id(CS_POST_MESH_VOLUME, 0);

    while (mesh_id != 0) {
      cs_post_write_var(mesh_id,
                        CS_POST_WRITER_ALL_ASSOCIATED,
                        "Cell_Volume",
                        1,
                        false,
                        true,
                        CS_POST_TYPE_cs_real_t,
                        mesh_quantities->cell_vol,
                        nullptr,
                        nullptr,
                        ts);

      mesh_id = cs_post_mesh_find_next_with_cat_id(CS_POST_MESH_VOLUME, mesh_id);
    }

  } /* End of cell volume treatment */

  /*-------------------------*/
  /* boundary cell thickness */
  /*-------------------------*/

  if (compute_thickness == true && mesh->n_g_b_faces > 0) {

    cs_real_t *b_thickness;
    CS_MALLOC_HD(b_thickness, mesh->n_b_faces, cs_real_t, cs_alloc_mode);

    cs_mesh_quantities_b_thickness_f(mesh, mesh_quantities, 0, b_thickness);

    /* Display histograms */

    bft_printf(_("\n  Histogram of boundary cell thickness:\n\n"));
    _histogram(mesh->n_b_faces, b_thickness);

    /* Post processing */

    int mesh_id = cs_post_mesh_find_next_with_cat_id(CS_POST_MESH_BOUNDARY, 0);

    while (mesh_id != 0) {
      cs_post_write_var(CS_POST_MESH_BOUNDARY,
                        CS_POST_WRITER_ALL_ASSOCIATED,
                        "Boundary_Cell_Thickness",
                        1,
                        false,
                        true,
                        CS_POST_TYPE_cs_real_t,
                        nullptr,
                        nullptr,
                        b_thickness,
                        ts);

      mesh_id = cs_post_mesh_find_next_with_cat_id(CS_POST_MESH_BOUNDARY,
                                                   mesh_id);
    }

    CS_FREE(b_thickness);

  } /* End of boundary cell thickness treatment */

  /*---------------*/
  /* warping error */
  /*---------------*/

  if (compute_warp_error == true) {

    cs_real_t  *warp_error = nullptr;
    CS_MALLOC_HD(warp_error, mesh->n_cells ,cs_real_t, cs_alloc_mode);

    _compute_warp_error(mesh,
                        cs_glob_mesh_adjacencies,
                        mesh_quantities,
                        ctx,
                        warp_error);

    ctx.wait();

    /* Display histograms */

    bft_printf(_("\n  Histogram of the cellwise warping error :\n\n"));
    _histogram(n_cells, warp_error);

    /* Post processing */

    int mesh_id = cs_post_mesh_find_next_with_cat_id(CS_POST_MESH_VOLUME, 0);

    while (mesh_id != 0) {
      cs_post_write_var(mesh_id,
                        CS_POST_WRITER_ALL_ASSOCIATED,
                        "Warp_Error",
                        1,
                        false,
                        true,
                        CS_POST_TYPE_cs_real_t,
                        warp_error,
                        nullptr,
                        nullptr,
                        ts);

      mesh_id = cs_post_mesh_find_next_with_cat_id(CS_POST_MESH_VOLUME, mesh_id);
    }

    double  l2_error = cs_gres(mesh->n_cells,
                               mesh_quantities->cell_vol,
                               warp_error,
                               warp_error);

    if (l2_error > 0) /* May be < 0 in case of negative volumes */
      bft_printf(" L2-error norm induced by warping : %5.3e\n", sqrt(l2_error));

    CS_FREE(warp_error);

  } /* End of cell volume treatment */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
