/*============================================================================
 * Gradient reconstruction at boundaries.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_ext_neighborhood.h"
#include "cs_field.h"
#include "cs_halo.h"
#include "cs_internal_coupling.h"
#include "cs_gradient_priv.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_adjacencies.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gradient_boundary.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*!
 * \file cs_gradient_boundary.c
 * \brief Gradient reconstruction at boundaries and associated functions.
 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local macros
 *============================================================================*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Add compute 3x3 cocg for least squares algorithm contribution from hidden
 * faces (as pure homogeneous Neumann BC's) to a single cell
 *
 * parameters:
 *   c_id               <-- cell id
 *   cell_hb_faces_idx  <-- cells -> hidden boundary faces index
 *   cell_hb_faces      <-- cells -> hidden boundary faces adjacency
 *   b_face_surf        <-- boundary faces surfaces
 *   b_face_normal      <-- boundary faces normals
 *   cocg               <-> cocg covariance matrix for given cell
 *----------------------------------------------------------------------------*/

static inline void
_add_hb_faces_cocg_lsq_cell(cs_lnum_t        c_id,
                            const cs_lnum_t  cell_hb_faces_idx[],
                            const cs_lnum_t  cell_hb_faces[],
                            const cs_real_t  b_face_surf[],
                            const cs_real_t  b_face_normal[][3],
                            cs_cocg_t        cocg[6])

{
  cs_lnum_t s_id = cell_hb_faces_idx[c_id];
  cs_lnum_t e_id = cell_hb_faces_idx[c_id+1];

  for (cs_lnum_t i = s_id; i < e_id; i++) {

    cs_lnum_t f_id = cell_hb_faces[i];

    cs_real_t udbfs = 1. / b_face_surf[f_id];

    cs_real_t dddij[3];
    for (cs_lnum_t ll = 0; ll < 3; ll++)
      dddij[ll] = udbfs * b_face_normal[f_id][ll];

    cocg[0] += dddij[0]*dddij[0];
    cocg[1] += dddij[1]*dddij[1];
    cocg[2] += dddij[2]*dddij[2];
    cocg[3] += dddij[0]*dddij[1];
    cocg[4] += dddij[1]*dddij[2];
    cocg[5] += dddij[0]*dddij[2];

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the values of a scalar at boundary face I' positions
 *         using least-squares interpolation.
 *
 * This assumes ghost cell values which might be used are already
 * synchronized.
 *
 * A simple limiter is applied to ensure the maximum principle is preserved
 * (using non-reconstructed values in case of non-homogeneous Neumann
 * conditions).
 *
 * This function uses a local iterative approach to compute the cell gradient,
 * as handling of the boundary condition terms b in higher dimensions
 * would otherwise require solving higher-dimensional systems, often at
 * a higher cost.
 *
 * \remark
 *
 * To compute the values at I', we only need the gradient along II', so
 * in most cases, we could simply assume a Neuman BC.
 *
 * The same logic is applied as for \ref cs_gradient_boundary_iprime_lsq_s.
 *
 * \param[in]   m               pointer to associated mesh structure
 * \param[in]   fvq             pointer to associated finite volume quantities
 * \param[in]   cpl             structure associated with internal coupling,
 *                              or NULL
 * \param[in]   n_faces         number of faces at which to compute values
 * \param[in]   face_ids        ids of boundary faces at which to compute
 *                              values, or NULL for all
 * \param[in]   halo_type       halo (cell neighborhood) type
 * \param[in]   var_dim         variable dimension
 * \param[in]   clip_coeff      clipping (limiter) coefficient
 *                              (no limiter if < 0)
 * \param[in]   bc_coeff_a      boundary condition term a, or NULL
 * \param[in]   bc_coeff_b      boundary condition term b, or NULL
 * \param[in]   c_weight        cell variable weight, or NULL
 * \param[in]   var             variable values et cell centers
 * \param[out]  var_iprime      variable values et face iprime locations
 */
/*----------------------------------------------------------------------------*/

static void
_gradient_b_faces_iprime_strided_lsq(const cs_mesh_t               *m,
                                     const cs_mesh_quantities_t    *fvq,
                                     const cs_internal_coupling_t  *cpl,
                                     cs_lnum_t                      n_faces,
                                     const cs_lnum_t               *face_ids,
                                     cs_halo_type_t                 halo_type,
                                     cs_lnum_t                      var_dim,
                                     double                         clip_coeff,
                                     const cs_real_t               *bc_coeff_a,
                                     const cs_real_t               *bc_coeff_b,
                                     const cs_real_t                c_weight[],
                                     const cs_real_t                var[],
                                     cs_real_t            *restrict var_iprime)
{
  /* Initialization */

  const int n_it_max = 50;
  const cs_real_t  eps_cvg = 1e-5;
  const cs_real_t  eps_cvg2 = cs_math_pow2(eps_cvg);
  const cs_real_t  eps_dvg = 1e-2;
  const cs_real_t  eps_dvg2 = cs_math_pow2(eps_dvg);

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
  const cs_lnum_t *restrict cell_cells_idx
    = (const cs_lnum_t *restrict) ma->cell_cells_idx;
  const cs_lnum_t *restrict cell_cells_e_idx
    = (const cs_lnum_t *restrict) ma->cell_cells_e_idx;
  const cs_lnum_t *restrict cell_b_faces_idx
    = (const cs_lnum_t *restrict) ma->cell_b_faces_idx;
  const cs_lnum_t *restrict cell_cells
    = (const cs_lnum_t *restrict) ma->cell_cells;
  const cs_lnum_t *restrict cell_cells_e
    = (const cs_lnum_t *restrict) ma->cell_cells_e;
  const cs_lnum_t *restrict cell_b_faces
    = (const cs_lnum_t *restrict) ma->cell_b_faces;

  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;

  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;
  const cs_real_t *restrict b_face_surf
    = (const cs_real_t *restrict)fvq->b_face_surf;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)fvq->b_face_cog;
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;

  cs_lnum_t   hn_stride = 0;
  const bool _hn_faces[1] = {false};
  const bool  *hn_faces = _hn_faces;
  if (cpl != NULL) {
    hn_stride = 1;
    hn_faces = (const bool *)cpl->coupled_faces;
  }

  assert(var_dim <= 9);  /* Local arrays with hard-coded dimensions follow. */

  /* Loop on selected boundary faces */

  #pragma omp parallel for schedule(dynamic, CS_THR_MIN)
  for (cs_lnum_t f_idx = 0; f_idx < n_faces; f_idx++) {

    cs_lnum_t f_id = (face_ids != NULL) ? face_ids[f_idx] : f_idx;
    cs_lnum_t c_id = m->b_face_cells[f_id];

    /* Reconstruct gradients using least squares for non-orthogonal meshes */

    cs_real_t cocg[6] = {0., 0., 0., 0., 0., 0.};
    cs_real_t rhs[9][3];
    cs_real_t var_i[9], var_min[9], var_max[9];

    for (cs_lnum_t ll = 0; ll < var_dim; ll++) {
      rhs[ll][0] = 0;
      rhs[ll][1] = 0;
      rhs[ll][2] = 0;
    }

    for (cs_lnum_t ll = 0; ll < var_dim; ll++) {
      var_i[ll] = var[c_id*var_dim + ll];
      var_min[ll] = var_i[ll];
      var_max[ll] = var_i[ll];
    }

    /* Contribution from adjacent cells */

    int n_adj = (halo_type == CS_HALO_EXTENDED) ? 2 : 1;

    for (int adj_id = 0; adj_id < n_adj; adj_id++) {

      const cs_lnum_t *restrict cell_cells_p;
      cs_lnum_t s_id, e_id;

      if (adj_id == 0){
        s_id = cell_cells_idx[c_id];
        e_id = cell_cells_idx[c_id+1];
        cell_cells_p = (const cs_lnum_t *restrict)(cell_cells);
      }
      else if (cell_cells_e_idx != NULL){
        s_id = cell_cells_e_idx[c_id];
        e_id = cell_cells_e_idx[c_id+1];
        cell_cells_p = (const cs_lnum_t *restrict)(cell_cells_e);
      }
      else
        break;

      if (c_weight == NULL) {

        for (cs_lnum_t i = s_id; i < e_id; i++) {

          cs_real_t dc[3];
          cs_lnum_t c_id1 = cell_cells_p[i];
          for (cs_lnum_t ii = 0; ii < 3; ii++)
            dc[ii] = cell_cen[c_id1][ii] - cell_cen[c_id][ii];

          cs_real_t ddc = 1. / cs_math_3_square_norm(dc);

          cocg[0] += dc[0]*dc[0]*ddc;
          cocg[1] += dc[1]*dc[1]*ddc;
          cocg[2] += dc[2]*dc[2]*ddc;
          cocg[3] += dc[0]*dc[1]*ddc;
          cocg[4] += dc[1]*dc[2]*ddc;
          cocg[5] += dc[0]*dc[2]*ddc;

          const cs_real_t *var_j = var + (c_id1*var_dim);

          for (cs_lnum_t kk = 0; kk < var_dim; kk++) {
            cs_real_t pfac = (var_j[kk] - var_i[kk]) * ddc;
            for (cs_lnum_t ll = 0; ll < 3; ll++)
              rhs[kk][ll] += dc[ll] * pfac;
          }
          for (cs_lnum_t ll = 0; ll < var_dim; ll++) {
            var_min[ll] = CS_MIN(var_min[ll], var_j[ll]);
            var_max[ll] = CS_MAX(var_max[ll], var_j[ll]);
          }

        }

      }
      else {

        for (cs_lnum_t i = s_id; i < e_id; i++) {

          cs_real_t dc[3];
          cs_lnum_t c_id1 = cell_cells_p[i];
          for (cs_lnum_t ii = 0; ii < 3; ii++)
            dc[ii] = cell_cen[c_id1][ii] - cell_cen[c_id][ii];

          cs_real_t ddc = 1. / cs_math_3_square_norm(dc);

          const cs_real_t *var_j = var + (c_id1*var_dim);

          cocg[0] += dc[0]*dc[0]*ddc;
          cocg[1] += dc[1]*dc[1]*ddc;
          cocg[2] += dc[2]*dc[2]*ddc;
          cocg[3] += dc[0]*dc[1]*ddc;
          cocg[4] += dc[1]*dc[2]*ddc;
          cocg[5] += dc[0]*dc[2]*ddc;

          cs_real_t _weight =   2. * c_weight[c_id1]
                              / (c_weight[c_id] + c_weight[c_id1]);

          for (cs_lnum_t kk = 0; kk < var_dim; kk++) {
            cs_real_t pfac = (var_j[kk] - var_i[kk]) * ddc;
            for (cs_lnum_t ll = 0; ll < 3; ll++)
              rhs[kk][ll] += dc[ll] * pfac * _weight;
          }

          for (cs_lnum_t ll = 0; ll < var_dim; ll++) {
            var_min[ll] = CS_MIN(var_min[ll], var_j[ll]);
            var_max[ll] = CS_MAX(var_max[ll], var_j[ll]);
          }
        }
      }

    } /* End of contribution from interior and extended cells */

    /* Contribution from hidden boundary faces */

    if (ma->cell_hb_faces_idx != NULL)
      _add_hb_faces_cocg_lsq_cell(c_id,
                                  ma->cell_hb_faces_idx,
                                  ma->cell_hb_faces,
                                  fvq->b_face_surf,
                                  (const cs_real_3_t *)fvq->b_face_normal,
                                  cocg);

    /* Contribution from boundary faces. */

    /* Use local ids for buffering: 12 should be more than enough for all but
       the most complex (and probably bad quality) cells,
       but if we go beyond that, recompute some values in loop */

    const cs_lnum_t n_coeff_b_contrib_buf_max = 12;

    cs_lnum_t n_coeff_b_contrib = 0;
    bool coeff_b_contrib_f[n_coeff_b_contrib_buf_max];

    cs_lnum_t s_id = cell_b_faces_idx[c_id];
    cs_lnum_t e_id = cell_b_faces_idx[c_id+1];

    for (cs_lnum_t i = s_id; i < e_id; i++) {

      bool coeff_b_contrib = false;

      cs_lnum_t c_f_id = cell_b_faces[i];

      cs_real_t dif[3];
      cs_real_t ddif;

      /* For coupled faces, use pure Neumann condition */

      if (hn_faces[c_f_id * hn_stride]) {

        cs_real_t udbfs = 1. / b_face_surf[c_f_id];
        for (cs_lnum_t ll = 0; ll < 3; ll++)
          dif[ll] = udbfs * b_face_normal[c_f_id][ll];

        ddif = 1. / cs_math_3_square_norm(dif);

      }
      else {

        for (cs_lnum_t ll = 0; ll < 3; ll++)
          dif[ll] = b_face_cog[c_f_id][ll] - cell_cen[c_id][ll];

        ddif = 1. / cs_math_3_square_norm(dif);

        cs_real_t var_f[9];

        const cs_real_t *a = bc_coeff_a + (c_f_id*var_dim);
        const cs_real_t *b = bc_coeff_b + (c_f_id*var_dim*var_dim);

        for (cs_lnum_t kk = 0; kk < var_dim; kk++) {
          var_f[kk] = a[kk];
          for (cs_lnum_t ll = 0; ll < var_dim; ll++) {
            var_f[kk] += b[kk*var_dim + ll] * var_i[ll];
          }
        }

        for (cs_lnum_t ll = 0; ll < var_dim; ll++) {
          var_min[ll] = CS_MIN(var_min[ll], var_f[ll]);
          var_max[ll] = CS_MAX(var_max[ll], var_f[ll]);
          if (fabs(var_f[ll] - a[ll]) > 1e-24)
            coeff_b_contrib = true;
        }

        for (cs_lnum_t kk = 0; kk < var_dim; kk++) {
          cs_real_t pfac = (var_f[kk] - var_i[kk]) * ddif;
          for (cs_lnum_t ll = 0; ll < 3; ll++)
            rhs[kk][ll] += dif[ll] * pfac;
        }

        cs_lnum_t i_rel = i - s_id;
        if (i_rel < n_coeff_b_contrib_buf_max)
          coeff_b_contrib_f[i_rel] = coeff_b_contrib;

        if (coeff_b_contrib)
          n_coeff_b_contrib += 1;

      }

      cocg[0] += dif[0]*dif[0]*ddif;
      cocg[1] += dif[1]*dif[1]*ddif;
      cocg[2] += dif[2]*dif[2]*ddif;
      cocg[3] += dif[0]*dif[1]*ddif;
      cocg[4] += dif[1]*dif[2]*ddif;
      cocg[5] += dif[0]*dif[2]*ddif;

    } /* End of contribution from boundary faces */

    /* Invert local covariance matrix */

    cs_real_t a00 = cocg[1]*cocg[2] - cocg[4]*cocg[4];
    cs_real_t a01 = cocg[4]*cocg[5] - cocg[3]*cocg[2];
    cs_real_t a02 = cocg[3]*cocg[4] - cocg[1]*cocg[5];
    cs_real_t a11 = cocg[0]*cocg[2] - cocg[5]*cocg[5];
    cs_real_t a12 = cocg[3]*cocg[5] - cocg[0]*cocg[4];
    cs_real_t a22 = cocg[0]*cocg[1] - cocg[3]*cocg[3];

    cs_real_t det_inv = 1. / (cocg[0]*a00 + cocg[3]*a01 + cocg[5]*a02);

    cocg[0] = a00 * det_inv;
    cocg[1] = a11 * det_inv;
    cocg[2] = a22 * det_inv;
    cocg[3] = a01 * det_inv;
    cocg[4] = a12 * det_inv;
    cocg[5] = a02 * det_inv;

    /* First estimate of cell gradient.

       Caution: in this function, the local gradient storage is transposed
       relative to that used in usual global (by array) computations, so
       as to allow using a usually oversized (fixed at compilation) leading
       dimension while (9, where vectors require 3 and symmetric tensors 6). */

    cs_real_t grad[9][3];

    for (cs_lnum_t ll = 0; ll < var_dim; ll++) {

      grad[ll][0] = (  cocg[0] * rhs[ll][0]
                     + cocg[3] * rhs[ll][1]
                     + cocg[5] * rhs[ll][2]);
      grad[ll][1] = (  cocg[3] * rhs[ll][0]
                     + cocg[1] * rhs[ll][1]
                     + cocg[4] * rhs[ll][2]);
      grad[ll][2] = (  cocg[5] * rhs[ll][0]
                     + cocg[4] * rhs[ll][1]
                     + cocg[2] * rhs[ll][2]);

    }

    /* First estimation of reconstructed value at I' */

    cs_real_t  var_ip[9];

    for (cs_lnum_t ll = 0; ll < var_dim; ll++) {
      var_ip[ll] = var_i[ll] + cs_math_3_dot_product(grad[ll], diipb[f_id]);
    }

    /* Refine boundary value estimation iteratively to account for bc_coeffs */

    if (n_coeff_b_contrib > 0) {

      /* Compute norms for convergence testing. Note that we prefer to
         test convergence based on the variable at I' rather than of the
         gradient itself, as only the gradient component tangential to
         II' is used here, so convergence of all components may not be
         necessary.

         We do not know whether the value at I or F (not reconstructed)
         is the largest, so we take thev average if the 2 matching norms. */

      cs_real_t ref_norm2 = 0;

      if (var_dim == 3) {
        ref_norm2 =   cs_math_3_square_norm(var_i)
                    + cs_math_3_square_norm(var_ip);
      }
      else {
        for (cs_lnum_t ii = 0; ii < var_dim; ii++)
          ref_norm2 += cs_math_sq(var_i[ii]) + cs_math_sq(var_ip[ii]);
      }

      cs_real_t var_ip_0[9], var_ip_f[9];
      cs_real_t grad_0[9][3];
      cs_real_t var_ip_d_norm2 = 0;

      memcpy(var_ip_0, var_ip, sizeof(cs_real_t)*var_dim);
      memcpy(grad_0, grad, sizeof(cs_real_t)*var_dim*3);

      /* Iterate over boundary condition contributions. */

      for (int n_c_it = 0; n_c_it < n_it_max; n_c_it++) {

        cs_real_t rhs_c[9][3];

        for (cs_lnum_t ll = 0; ll < var_dim; ll++) {
          rhs_c[ll][0] = 0;
          rhs_c[ll][1] = 0;
          rhs_c[ll][2] = 0;
        }

        /* Loop on boundary faces */

        for (cs_lnum_t i = s_id; i < e_id; i++) {

          cs_lnum_t i_rel = i - s_id;
          if (i_rel < n_coeff_b_contrib_buf_max) {
            if (coeff_b_contrib_f[i_rel] == false)
              continue;
          }

          cs_lnum_t c_f_id = cell_b_faces[i];

          /* Remark: we could avoid recomputing dif and ddif in most cases by
             saving at least a few values from the previous step, in
             fixed-size buffers indexed by i_rel. */

          cs_real_t dif[3];
          for (cs_lnum_t ii = 0; ii < 3; ii++)
            dif[ii] = b_face_cog[c_f_id][ii] - cell_cen[c_id][ii];

          cs_real_t ddif = 1. / cs_math_3_square_norm(dif);

          /* Note that the contribution to the right-hand side from
             bc_coeff_a[c_f_id] + (bc_coeff_b -1).var[c_id] has already been
             counted, so it does not appear in the following terms.
             We should perhaps check whether the numerical sensitivity
             is lower when proceeding thus or when computing the full
             face value at each step. */

          for (cs_lnum_t ll = 0; ll < var_dim; ll++) {
            var_ip_f[ll] = cs_math_3_dot_product(grad[ll], diipb[c_f_id]);
          }

          const cs_real_t *b = bc_coeff_b + (c_f_id*var_dim*var_dim);

          for (cs_lnum_t kk = 0; kk < var_dim; kk++) {
            cs_real_t pfac = 0;
            for (cs_lnum_t ll = 0; ll < var_dim; ll++) {
              pfac += b[kk*var_dim + ll] * var_ip_f[ll] * ddif;
            }

            for (cs_lnum_t ll = 0; ll < 3; ll++)
              rhs_c[kk][ll] += dif[ll] * pfac;
          }

        } /* End of loop on boundary faces */

        /* Compute gradient correction */

        cs_real_t grad_c[9][3];

        for (cs_lnum_t ii = 0; ii < var_dim; ii++) {

          grad_c[ii][0] = (  cocg[0] * rhs_c[ii][0]
                           + cocg[3] * rhs_c[ii][1]
                           + cocg[5] * rhs_c[ii][2]);
          grad_c[ii][1] = (  cocg[3] * rhs_c[ii][0]
                           + cocg[1] * rhs_c[ii][1]
                           + cocg[4] * rhs_c[ii][2]);
          grad_c[ii][2] = (  cocg[5] * rhs_c[ii][0]
                           + cocg[4] * rhs_c[ii][1]
                           + cocg[2] * rhs_c[ii][2]);

        }

        /* Update gradient and var_ip values */

        for (cs_lnum_t ii = 0; ii < var_dim; ii++) {
          for (cs_lnum_t jj = 0; jj < 3; jj++)
            grad[ii][jj] = grad_0[ii][jj] + grad_c[ii][jj];
        }

        cs_real_t var_ip_prv[9], var_ip_d[9];

        for (cs_lnum_t ll = 0; ll < var_dim; ll++) {
          var_ip_prv[ll] = var_ip[ll];
          var_ip[ll] = var_i[ll] + cs_math_3_dot_product(grad[ll], diipb[f_id]);
          var_ip_d[ll] = var_ip_prv[ll] - var_ip[ll];
        }

        if (var_dim == 3)
          var_ip_d_norm2 = cs_math_3_square_norm(var_ip_d);
        else {
          var_ip_d_norm2 = 0;
          for (cs_lnum_t ii = 0; ii < var_dim; ii++)
            var_ip_d_norm2 += cs_math_sq(var_ip_d[ii]);
        }

#if 0
        printf("grads %d: it %d, ref_norm %g, it_norm %g\n",
               c_id, n_c_it, sqrt(ref_norm2), sqrt(var_ip_d_norm2));
#endif

        /* Compare square of distances to avoid computing square root,
           so all thresholds are squared also (i.e. 1e-10 instead of 1e-5). */
        if (var_ip_d_norm2 < 1e-38)
          break;

        else if (var_ip_d_norm2 < eps_cvg2 * ref_norm2)
          break;

      } /* End of loop on iterations */

      /* If the last correction was too large, we suspect
         the the algorithm did not converge at all/diverged,
         so we simply restore the non-reconstructed value
         (additional precaution, not encountered in testing). */

      if (var_ip_d_norm2 > eps_dvg2 * ref_norm2) {
        memcpy(var_ip, var_i, var_dim*sizeof(cs_real_t));
#if 0
        printf("%s: non-convergence for face %ld\n"
               "  use non-recontruced value\n", __func__, (long)f_id);
#endif
      }

    } /* End of iterative correction for BC coeffs */

    /* Apply simple limiter */

    if (clip_coeff >= 0) {
      for (cs_lnum_t ii = 0; ii < var_dim; ii++) {
        cs_real_t d = var_max[ii] - var_min[ii];
        var_max[ii] += d*clip_coeff;
        var_min[ii] -= d*clip_coeff;

        if (var_ip[ii] < var_min[ii])
          var_ip[ii] = var_min[ii];
        if (var_ip[ii] > var_max[ii])
          var_ip[ii] = var_max[ii];
      }
    }

    /* Save final value */

    for (cs_lnum_t ii = 0; ii < var_dim; ii++)
      var_iprime[f_idx*var_dim + ii] = var_ip[ii];

  } /* End of loop on selected faces */
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the values of a scalar at boundary face I' positions
 *         using least-squares interpolation.
 *
 * This assumes ghost cell values for the variable (var) are up-to-date.
 *
 * A simple limiter is applied to ensure the maximum principle is preserved
 * (using non-reconstructed values in case of non-homogeneous Neumann
 * conditions).
 *
 * \remark
 *
 * To compute the values at I', we only need the gradient along II', so
 * in most cases, we could simply assume a Neumann BC for a given face.
 *
 * We still use the provided BC's when possible, for the following cases:
 * - Given a non-uniform Dirichlet condition and a non-orthogonal mesh,
 *   the Dirichlet values at face centers (shifted by II' relative to I')
 *   can convey a portion of the information of the gradient along II'.
 * - For cells with multiple boundary faces, information from faces whose
 *   normals are not orthogonal to II' can also provide a significant
 *   contribution to the normal.
 *
 * For coupled faces (with internal coupling), we assume a Neumann BC,
 * as the values that could be computed at the intersections of the
 * segments joining coupled cell centers and the matching boundary face
 * do not currently account for the non-linearity associated with wall laws.
 *
 * A more precise computation for cells with multiple coupled boundary faces
 * would be possible by using multiple passes, first estimating boundary face
 * values without recontruction, then using that information for faces which
 * are not tangential to II' in a second pass. Note though that as the mesh
 * is refined, the portion of the boundary adjacent to cells with
 * multiple coupled faces (i.e. at corners and ridges) diminishes.
 *
 * \param[in]   m               pointer to associated mesh structure
 * \param[in]   fvq             pointer to associated finite volume quantities
 * \param[in]   cpl             structure associated with internal coupling,
 *                              or NULL
 * \param[in]   n_faces         number of faces at which to compute values
 * \param[in]   face_ids        ids of boundary faces at which to compute
 *                              values, or NULL for all
 * \param[in]   halo_type       halo (cell neighborhood) type
 * \param[in]   clip_coeff      clipping (limiter) coefficient
 *                              (no limiter if < 0)
 * \param[in]   bc_coeff_a      boundary condition term a, or NULL
 * \param[in]   bc_coeff_b      boundary condition term b, or NULL
 * \param[in]   c_weight        cell variable weight, or NULL
 * \param[in]   var             variable values et cell centers
 * \param[out]  var_iprime      variable values et face iprime locations
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_boundary_iprime_lsq_s(const cs_mesh_t               *m,
                                  const cs_mesh_quantities_t    *fvq,
                                  const cs_internal_coupling_t  *cpl,
                                  cs_lnum_t                      n_faces,
                                  const cs_lnum_t               *face_ids,
                                  cs_halo_type_t                 halo_type,
                                  double                         clip_coeff,
                                  const cs_real_t               *bc_coeff_a,
                                  const cs_real_t               *bc_coeff_b,
                                  const cs_real_t                c_weight[],
                                  const cs_real_t                var[],
                                  cs_real_t           *restrict  var_iprime)
{
  /* Initialization */

  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;

  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

  const cs_lnum_t *restrict cell_cells_idx
    = (const cs_lnum_t *restrict) ma->cell_cells_idx;
  const cs_lnum_t *restrict cell_cells_e_idx
    = (const cs_lnum_t *restrict) ma->cell_cells_e_idx;
  const cs_lnum_t *restrict cell_b_faces_idx
    = (const cs_lnum_t *restrict) ma->cell_b_faces_idx;
  const cs_lnum_t *restrict cell_cells
    = (const cs_lnum_t *restrict) ma->cell_cells;
  const cs_lnum_t *restrict cell_cells_e
    = (const cs_lnum_t *restrict) ma->cell_cells_e;
  const cs_lnum_t *restrict cell_b_faces
    = (const cs_lnum_t *restrict) ma->cell_b_faces;

  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;

  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;
  const cs_real_t *restrict b_face_surf
    = (const cs_real_t *restrict)fvq->b_face_surf;
  const cs_real_t *restrict b_dist
    = (const cs_real_t *restrict)fvq->b_dist;
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)fvq->diipb;

  cs_lnum_t   hn_stride = 0;
  bool _hn_faces[1] = {false};
  const bool  *hn_faces = _hn_faces;
  if (bc_coeff_a == NULL && bc_coeff_b == NULL) {
    hn_stride = 0;
    _hn_faces[0] = true;
  }
  else if (cpl != NULL) {
    hn_stride = 1;
    hn_faces = (const bool *)cpl->coupled_faces;
  }

  /* Loop on selected boundary faces */

  #pragma omp parallel for schedule(dynamic, CS_THR_MIN)
  for (cs_lnum_t f_idx = 0; f_idx < n_faces; f_idx++) {

    cs_lnum_t f_id = (face_ids != NULL) ? face_ids[f_idx] : f_idx;
    cs_lnum_t c_id = b_face_cells[f_id];

    /* Reconstruct gradients using least squares for non-orthogonal meshes */

    cs_real_t cocg[6] = {0., 0., 0., 0., 0., 0.};
    cs_real_t rhs[3] = {0, 0, 0};

    cs_real_t var_i = var[c_id];
    cs_real_t var_min = var_i, var_max = var_i;

    /* Contribution from adjacent cells */

    int n_adj = (halo_type == CS_HALO_EXTENDED) ? 2 : 1;

    for (int adj_id = 0; adj_id < n_adj; adj_id++) {

      const cs_lnum_t *restrict cell_cells_p;
      cs_lnum_t s_id, e_id;

      if (adj_id == 0){
        s_id = cell_cells_idx[c_id];
        e_id = cell_cells_idx[c_id+1];
        cell_cells_p = (const cs_lnum_t *restrict)(cell_cells);
      }
      else if (cell_cells_e_idx != NULL){
        s_id = cell_cells_e_idx[c_id];
        e_id = cell_cells_e_idx[c_id+1];
        cell_cells_p = (const cs_lnum_t *restrict)(cell_cells_e);
      }
      else
        break;

      if (c_weight == NULL) {

        for (cs_lnum_t i = s_id; i < e_id; i++) {

          cs_real_t dc[3];
          cs_lnum_t c_id1 = cell_cells_p[i];
          for (cs_lnum_t ii = 0; ii < 3; ii++)
            dc[ii] = cell_cen[c_id1][ii] - cell_cen[c_id][ii];

          cs_real_t ddc = 1. / cs_math_3_square_norm(dc);

          cocg[0] += dc[0]*dc[0]*ddc;
          cocg[1] += dc[1]*dc[1]*ddc;
          cocg[2] += dc[2]*dc[2]*ddc;
          cocg[3] += dc[0]*dc[1]*ddc;
          cocg[4] += dc[1]*dc[2]*ddc;
          cocg[5] += dc[0]*dc[2]*ddc;

          cs_real_t var_j = var[c_id1];
          var_min = CS_MIN(var_min, var_j);
          var_max = CS_MAX(var_max, var_j);

          cs_real_t pfac = (var_j - var_i) * ddc;
          for (cs_lnum_t ll = 0; ll < 3; ll++)
            rhs[ll] += dc[ll] * pfac;

        }

      }
      else {

        for (cs_lnum_t i = s_id; i < e_id; i++) {

          cs_real_t dc[3];
          cs_lnum_t c_id1 = cell_cells_p[i];
          for (cs_lnum_t ii = 0; ii < 3; ii++)
            dc[ii] = cell_cen[c_id1][ii] - cell_cen[c_id][ii];

          cs_real_t ddc = 1. / cs_math_3_square_norm(dc);

          cocg[0] += dc[0]*dc[0]*ddc;
          cocg[1] += dc[1]*dc[1]*ddc;
          cocg[2] += dc[2]*dc[2]*ddc;
          cocg[3] += dc[0]*dc[1]*ddc;
          cocg[4] += dc[1]*dc[2]*ddc;
          cocg[5] += dc[0]*dc[2]*ddc;

          cs_real_t _weight =   2. * c_weight[c_id1]
                              / (c_weight[c_id] + c_weight[c_id1]);

          cs_real_t var_j = var[c_id1];
          var_min = CS_MIN(var_min, var_j);
          var_max = CS_MAX(var_max, var_j);

          cs_real_t pfac = (var_j - var_i) * ddc;
          for (cs_lnum_t ll = 0; ll < 3; ll++)
            rhs[ll] += dc[ll] * pfac * _weight;

        }
      }

    } /* End of contribution from interior and extended cells */

    /* Contribution from hidden boundary faces */

    if (ma->cell_hb_faces_idx != NULL)
      _add_hb_faces_cocg_lsq_cell(c_id,
                                  ma->cell_hb_faces_idx,
                                  ma->cell_hb_faces,
                                  fvq->b_face_surf,
                                  (const cs_real_3_t *)fvq->b_face_normal,
                                  cocg);

    /* Contribution from boundary faces. */

    cs_lnum_t s_id = cell_b_faces_idx[c_id];
    cs_lnum_t e_id = cell_b_faces_idx[c_id+1];

    for (cs_lnum_t i = s_id; i < e_id; i++) {

      cs_lnum_t c_f_id = cell_b_faces[i];

      cs_real_t dddij[3];

      /* For coupled faces, use pure Neumann condition,
         for other faces, use regular BC's */

      cs_real_t a, b;
      if (hn_faces[c_f_id * hn_stride]) {
        a = var_i;
        b = 0;
      }
      else {
        a = bc_coeff_a[c_f_id];
        b = bc_coeff_b[c_f_id];

        /* Use unreconstructed value for limiter */
        cs_real_t var_f = a + b*var_i;
        var_min = CS_MIN(var_min, var_f);
        var_max = CS_MAX(var_max, var_f);
      }

      cs_real_t unddij = 1. / b_dist[c_f_id];
      cs_real_t udbfs = 1. / b_face_surf[c_f_id];
      cs_real_t umcbdd = (1. - b) * unddij;

      for (cs_lnum_t ll = 0; ll < 3; ll++) {
        dddij[ll] =   udbfs * b_face_normal[c_f_id][ll]
                    + umcbdd * diipb[c_f_id][ll];
      }

      cocg[0] += dddij[0]*dddij[0];
      cocg[1] += dddij[1]*dddij[1];
      cocg[2] += dddij[2]*dddij[2];
      cocg[3] += dddij[0]*dddij[1];
      cocg[4] += dddij[1]*dddij[2];
      cocg[5] += dddij[0]*dddij[2];

      cs_real_t pfac = (a + (b-1.)*var_i) * unddij;

      for (cs_lnum_t ll = 0; ll < 3; ll++)
        rhs[ll] += dddij[ll] * pfac;

    } /* End of contribution from boundary faces */

    /* Invert local covariance matrix */

    cs_real_t a00 = cocg[1]*cocg[2] - cocg[4]*cocg[4];
    cs_real_t a01 = cocg[4]*cocg[5] - cocg[3]*cocg[2];
    cs_real_t a02 = cocg[3]*cocg[4] - cocg[1]*cocg[5];
    cs_real_t a11 = cocg[0]*cocg[2] - cocg[5]*cocg[5];
    cs_real_t a12 = cocg[3]*cocg[5] - cocg[0]*cocg[4];
    cs_real_t a22 = cocg[0]*cocg[1] - cocg[3]*cocg[3];

    cs_real_t det_inv = 1. / (cocg[0]*a00 + cocg[3]*a01 + cocg[5]*a02);

    cs_real_t grad[3];

    grad[0] = (  a00 * rhs[0]
               + a01 * rhs[1]
               + a02 * rhs[2]) * det_inv;
    grad[1] = (  a01 * rhs[0]
               + a11 * rhs[1]
               + a12 * rhs[2]) * det_inv;
    grad[2] = (  a02 * rhs[0]
               + a12 * rhs[1]
               + a22 * rhs[2]) * det_inv;

    /* Finally, reconstruct value at I' */

    cs_real_t  var_ip = var_i + (  grad[0]*diipb[f_id][0]
                                 + grad[1]*diipb[f_id][1]
                                 + grad[2]*diipb[f_id][2]);

    /* Apply simple limiter */

    if (clip_coeff >= 0) {
      cs_real_t d = var_max - var_min;
      var_max += d*clip_coeff;
      var_min -= d*clip_coeff;

      if (var_ip < var_min)
        var_ip = var_min;
      if (var_ip > var_max)
        var_ip = var_max;
    }

    var_iprime[f_idx] = var_ip;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the values of a vector at boundary face I' positions
 *         using least-squares interpolation.
 *
 * This assumes ghost cell values which might be used are already
 * synchronized.
 *
 * A simple limiter is applied to ensure the maximum principle is preserved
 * (using non-reconstructed values in case of non-homogeneous Neumann
 * conditions).
 *
 * This function uses a local iterative approach to compute the cell gradient,
 * as handling of the boundary condition terms b in higher dimensions
 * would otherwise require solving higher-dimensional systems, often at
 * a higher cost.
 *
 * \remark
 *
 * To compute the values at I', we only need the gradient along II', so
 * in most cases, we could simply assume a Neuman BC.
 *
 * The same logic is applied as for \ref cs_gradient_boundary_iprime_lsq_s.
 *
 * \param[in]   m               pointer to associated mesh structure
 * \param[in]   fvq             pointer to associated finite volume quantities
 * \param[in]   cpl             structure associated with internal coupling,
 *                              or NULL
 * \param[in]   n_faces         number of faces at which to compute values
 * \param[in]   face_ids        ids of boundary faces at which to compute
 *                              values, or NULL for all
 * \param[in]   halo_type       halo (cell neighborhood) type
 * \param[in]   clip_coeff      clipping (limiter) coefficient
 *                              (no limiter if < 0)
 * \param[in]   bc_coeff_a      boundary condition term a, or NULL
 * \param[in]   bc_coeff_b      boundary condition term b, or NULL
 * \param[in]   c_weight        cell variable weight, or NULL
 * \param[in]   var             variable values et cell centers
 * \param[out]  var_iprime      variable values et face iprime locations
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_boundary_iprime_lsq_v(const cs_mesh_t               *m,
                                  const cs_mesh_quantities_t    *fvq,
                                  const cs_internal_coupling_t  *cpl,
                                  cs_lnum_t                      n_faces,
                                  const cs_lnum_t               *face_ids,
                                  cs_halo_type_t                 halo_type,
                                  double                         clip_coeff,
                                  const cs_real_t               *bc_coeff_a[3],
                                  const cs_real_t               *bc_coeff_b[3][3],
                                  const cs_real_t                c_weight[],
                                  const cs_real_t                var[][6],
                                  cs_real_t            *restrict var_iprime[3])
{
  _gradient_b_faces_iprime_strided_lsq(m,
                                       fvq,
                                       cpl,
                                       n_faces,
                                       face_ids,
                                       halo_type,
                                       3, /* var_dim */
                                       clip_coeff,
                                       (const cs_real_t *)bc_coeff_a,
                                       (const cs_real_t *)bc_coeff_b,
                                       c_weight,
                                       (const cs_real_t *)var,
                                       (cs_real_t *) var_iprime);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the values of a symmetric tensor at boundary face
 *         I' positions using least-squares interpolation.
 *
 * This assumes ghost cell values which might be used are already
 * synchronized.
 *
 * A simple limiter is applied to ensure the maximum principle is preserved
 * (using non-reconstructed values in case of non-homogeneous Neumann
 * conditions).
 *
 * This function uses a local iterative approach to compute the cell gradient,
 * as handling of the boundary condition terms b in higher dimensions
 * would otherwise require solving higher-dimensional systems, often at
 * a higher cost.
 *
 * \remark
 *
 * To compute the values at I', we only need the gradient along II', so
 * in most cases, we could simply assume a Neuman BC.
 *
 * The same logic is applied as for \ref cs_gradient_boundary_iprime_lsq_s.
 *
 * \param[in]   m               pointer to associated mesh structure
 * \param[in]   fvq             pointer to associated finite volume quantities
 * \param[in]   cpl             structure associated with internal coupling,
 *                              or NULL
 * \param[in]   n_faces         number of faces at which to compute values
 * \param[in]   face_ids        ids of boundary faces at which to compute
 *                              values, or NULL for all
 * \param[in]   halo_type       halo (cell neighborhood) type
 * \param[in]   clip_coeff      clipping (limiter) coefficient
 *                              (no limiter if < 0)
 * \param[in]   bc_coeff_a      boundary condition term a, or NULL
 * \param[in]   bc_coeff_b      boundary condition term b, or NULL
 * \param[in]   c_weight        cell variable weight, or NULL
 * \param[in]   var             variable values et cell centers
 * \param[out]  var_iprime      variable values et face iprime locations
 */
/*----------------------------------------------------------------------------*/

void
cs_gradient_boundary_iprime_lsq_t(const cs_mesh_t               *m,
                                  const cs_mesh_quantities_t    *fvq,
                                  const cs_internal_coupling_t  *cpl,
                                  cs_lnum_t                      n_faces,
                                  const cs_lnum_t               *face_ids,
                                  cs_halo_type_t                 halo_type,
                                  double                         clip_coeff,
                                  const cs_real_t               *bc_coeff_a[6],
                                  const cs_real_t               *bc_coeff_b[6][6],
                                  const cs_real_t                c_weight[],
                                  const cs_real_t                var[][6],
                                  cs_real_t            *restrict var_iprime[6])
{
  _gradient_b_faces_iprime_strided_lsq(m,
                                       fvq,
                                       cpl,
                                       n_faces,
                                       face_ids,
                                       halo_type,
                                       6, /* var_dim */
                                       clip_coeff,
                                       (const cs_real_t *)bc_coeff_a,
                                       (const cs_real_t *)bc_coeff_b,
                                       c_weight,
                                       (const cs_real_t *)var,
                                       (cs_real_t *restrict) var_iprime);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
