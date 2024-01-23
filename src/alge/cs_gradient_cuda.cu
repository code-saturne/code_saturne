/*============================================================================
 * Gradient reconstruction, CUDA implementations.
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

#include <assert.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#include <cuda_runtime_api.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"

#include "cs_base_accel.h"
#include "cs_base_cuda.h"
#include "cs_blas.h"
#include "cs_cell_to_vertex.h"
#include "cs_ext_neighborhood.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_math_cuda.cuh"
#include "cs_mesh.h"
#include "cs_mesh_adjacencies.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gradient.h"
#include "cs_gradient_priv.h"

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------*/

/*=============================================================================
 * Local macros
 *============================================================================*/

/*=============================================================================
 * Local definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Inverse a 3x3 symmetric matrix (with symmetric storage)
 *         in place, using Cramer's rule
 *
 * \param[in, out]  a   matrix to inverse
 */
/*----------------------------------------------------------------------------*/

__device__ static void
_math_6_inv_cramer_sym_in_place_cuda(cs_cocg_t in[6])
{
  auto in00 = in[1]*in[2] - in[4]*in[4];
  auto in01 = in[4]*in[5] - in[3]*in[2];
  auto in02 = in[3]*in[4] - in[1]*in[5];
  auto in11 = in[0]*in[2] - in[5]*in[5];
  auto in12 = in[3]*in[5] - in[0]*in[4];
  auto in22 = in[0]*in[1] - in[3]*in[3];

  double det_inv = 1. / (in[0]*in00 + in[3]*in01 + in[5]*in02);

  in[0] = in00 * det_inv;
  in[1] = in11 * det_inv;
  in[2] = in22 * det_inv;
  in[3] = in01 * det_inv;
  in[4] = in12 * det_inv;
  in[5] = in02 * det_inv;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign 1 to coefb values
 *
 * \param[out]  a   matrix to inverse
 */
/*----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
__global__ static void
_set_one_to_coeff_b(const cs_lnum_t  n_b_faces,
                    cs_real_t       *(bc_coeff_b)[stride][stride])
{
  cs_lnum_t c_idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (c_idx >= n_b_faces) {
    return;
  }

  cs_lnum_t f_id = c_idx / stride;
  size_t i = c_idx % stride;

  bc_coeff_b[f_id][i][i] = 1;
}

/*----------------------------------------------------------------------------
 * Recompute cocg at boundaries, using saved cocgb
 *----------------------------------------------------------------------------*/

#define INSTANTIATE_LSQ(name, stride) template void name <stride>  \
  (const cs_mesh_t               *m,                               \
   const cs_mesh_adjacencies_t   *madj,                            \
   const cs_mesh_quantities_t    *fvq,                             \
   const cs_halo_type_t           halo_type,                       \
   int                            inc,                             \
   int                            n_c_iter_max,                    \
   cs_real_t                      c_eps,                           \
   const cs_real_t (*restrict coefav)[stride],                     \
   const cs_real_t (*restrict coefbv)[stride][stride],             \
   const cs_real_t (*restrict pvar)[stride],                       \
   const cs_real_t      *restrict c_weight,                        \
   const cs_cocg_6_t    *restrict cocgb,                           \
   cs_cocg_6_t          *restrict cocg,                            \
   cs_real_t (*restrict gradv)[stride][3])

template <typename T>
__global__ static void
_compute_cocg_from_cocgb(cs_lnum_t         n_b_cells,
                         const cs_lnum_t  *b_cells,
                         T                *cocg,
                         T                *cocgb)
{
  size_t ii = blockIdx.x*blockDim.x + threadIdx.x;

  if (ii < n_b_cells) {
    cs_lnum_t c_id = b_cells[ii];
    auto a = cocg[c_id];
    auto b = cocgb[ii];
    a[0] = b[0];
    a[1] = b[1];
    a[2] = b[2];
    a[3] = b[3];
    a[4] = b[4];
    a[5] = b[5];
  }
}

/*----------------------------------------------------------------------------
 * Save COCGB at boundaries, using COCG.
 *----------------------------------------------------------------------------*/

template <typename T>
__global__ static void
_save_cocgb(cs_lnum_t         size,
            const cs_lnum_t  *b_cells,
            T                *cocg,
            T                *cocgb)
{
  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
  if (ii < size) {
    cs_lnum_t c_id = b_cells[ii];
    auto      a    = cocgb[ii];
    auto      b    = cocg[c_id];
    a[0]           = b[0];
    a[1]           = b[1];
    a[2]           = b[2];
    a[3]           = b[3];
    a[4]           = b[4];
    a[5]           = b[5];
  }
}

/*----------------------------------------------------------------------------
 * Recompute the inverse of cocg
 *----------------------------------------------------------------------------*/

template <typename T>
__global__ static void
_compute_cocg_inv(cs_lnum_t         n_cells,
                  const cs_lnum_t  *cell_ids,
                  T                *cocg)
{
  cs_lnum_t c_idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (c_idx >= n_cells)
    return;

  cs_lnum_t c_id = (cell_ids != NULL) ? cell_ids[c_idx] : c_idx;

  auto a = cocg[c_id];

  auto a00 = a[1]*a[2] - a[4]*a[4];
  auto a01 = a[4]*a[5] - a[3]*a[2];
  auto a02 = a[3]*a[4] - a[1]*a[5];
  auto a11 = a[0]*a[2] - a[5]*a[5];
  auto a12 = a[3]*a[5] - a[0]*a[4];
  auto a22 = a[0]*a[1] - a[3]*a[3];

  double det_inv = 1. / (a[0]*a00 + a[3]*a01 + a[5]*a02);

  a[0] = a00 * det_inv;
  a[1] = a11 * det_inv;
  a[2] = a22 * det_inv;
  a[3] = a01 * det_inv;
  a[4] = a12 * det_inv;
  a[5] = a02 * det_inv;
}

/*----------------------------------------------------------------------------
 * Compute the gradient using the least-squares method.
 *----------------------------------------------------------------------------*/

template <typename T>
__global__ static void
_compute_gradient_lsq_s(cs_lnum_t     n_cells,
                        cs_real_3_t  *grad,
                        T            *cocg,
                        cs_real_4_t  *rhsv)
{
  size_t c_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (c_id < n_cells) {

    auto _cocg = cocg[c_id];
    auto _rhsv = rhsv[c_id];

    grad[c_id][0] =   _cocg[0] *_rhsv[0]
                    + _cocg[3] *_rhsv[1]
                    + _cocg[5] *_rhsv[2];
    grad[c_id][1] =   _cocg[3] *_rhsv[0]
                    + _cocg[1] *_rhsv[1]
                    + _cocg[4] *_rhsv[2];
    grad[c_id][2] =   _cocg[5] *_rhsv[0]
                    + _cocg[4] *_rhsv[1]
                    + _cocg[2] *_rhsv[2];

  }
}

/*----------------------------------------------------------------------------
 * Recompute the cocg and rhsv contributions from interior and extended cells
 *----------------------------------------------------------------------------*/

template <typename T>
__global__ static void
_compute_cocg_rhsv_lsq_s_i_face(cs_lnum_t           size,
                                T                  *cocg,
                                const cs_lnum_t    *cell_cells_idx,
                                const cs_lnum_t    *cell_cells,
                                const cs_real_3_t  *cell_cen,
                                cs_real_4_t        *rhsv,
                                const cs_real_t    *c_weight)
{
  cs_lnum_t c_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (c_id < size) {

    /* Initialize COCG (RHS initialize before) */

    cocg[c_id][0] = 0; cocg[c_id][1] = 0; cocg[c_id][2] = 0;
    cocg[c_id][3] = 0; cocg[c_id][4] = 0; cocg[c_id][5] = 0;

    cs_lnum_t s_id = cell_cells_idx[c_id];
    cs_lnum_t e_id = cell_cells_idx[c_id + 1];
    cs_real_t dc[3], ddc, _weight;
    cs_lnum_t c_id1;

    /* Add contributions from neighbor cells/interior faces */

    for (cs_lnum_t i = s_id; i < e_id; i++) {
      c_id1 = cell_cells[i];

      dc[0] = cell_cen[c_id1][0] - cell_cen[c_id][0];
      dc[1] = cell_cen[c_id1][1] - cell_cen[c_id][1];
      dc[2] = cell_cen[c_id1][2] - cell_cen[c_id][2];

      ddc = 1. / (dc[0] * dc[0] + dc[1] * dc[1] + dc[2] * dc[2]);
      if (c_weight == NULL)
        _weight = 1;
      else
        _weight = 2. * c_weight[c_id1] / (c_weight[c_id] + c_weight[c_id1]);

      _weight *= (rhsv[c_id1][3] - rhsv[c_id][3]) * ddc;

      rhsv[c_id][0] += dc[0] * _weight;
      rhsv[c_id][1] += dc[1] * _weight;
      rhsv[c_id][2] += dc[2] * _weight;

      cocg[c_id][0] += dc[0] * dc[0] * ddc;
      cocg[c_id][1] += dc[1] * dc[1] * ddc;
      cocg[c_id][2] += dc[2] * dc[2] * ddc;
      cocg[c_id][3] += dc[0] * dc[1] * ddc;
      cocg[c_id][4] += dc[1] * dc[2] * ddc;
      cocg[c_id][5] += dc[0] * dc[2] * ddc;
    }
  }
}

/*----------------------------------------------------------------------------
 * Compute cocg and rhsv contributions at boundaries
 *----------------------------------------------------------------------------*/

template <typename T>
__global__ static void
_compute_cocg_rhsv_lsq_s_b_face(cs_lnum_t         n_b_cells,
                                cs_real_t         inc,
                                const cs_lnum_t   b_cells[],
                                const cs_lnum_t   cell_b_faces_idx[],
                                const cs_lnum_t   cell_b_faces[],
                                const cs_real_t   b_face_normal[][3],
                                const cs_real_t   b_face_surf[],
                                const cs_real_t   b_dist[],
                                const cs_real_t   diipb[][3],
                                const cs_real_t   coefap[],
                                const cs_real_t   coefbp[],
                                T                *restrict cocg,
                                cs_real_4_t      *restrict rhsv)
{
  cs_lnum_t b_c_idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (b_c_idx >= n_b_cells)
    return;

  cs_lnum_t c_id = b_cells[b_c_idx];

  cs_lnum_t s_id = cell_b_faces_idx[c_id];
  cs_lnum_t e_id = cell_b_faces_idx[c_id + 1];

  for (cs_lnum_t i = s_id; i < e_id; i++) {
    cs_lnum_t f_id = cell_b_faces[i];

    cs_real_t unddij = 1. / b_dist[f_id];
    cs_real_t umcbdd = (1. - coefbp[f_id]) * unddij;
    cs_real_t udbfs  = 1. / b_face_surf[f_id];
    cs_real_t dddij[3];
    dddij[0] = udbfs * b_face_normal[f_id][0] + umcbdd * diipb[f_id][0];
    dddij[1] = udbfs * b_face_normal[f_id][1] + umcbdd * diipb[f_id][1];
    dddij[2] = udbfs * b_face_normal[f_id][2] + umcbdd * diipb[f_id][2];

    unddij *= (coefap[f_id] * inc + (coefbp[f_id] - 1.) * rhsv[c_id][3]);

    rhsv[c_id][0] += dddij[0] * unddij;
    rhsv[c_id][1] += dddij[1] * unddij;
    rhsv[c_id][2] += dddij[2] * unddij;

    cocg[c_id][0] += dddij[0] * dddij[0];
    cocg[c_id][1] += dddij[1] * dddij[1];
    cocg[c_id][2] += dddij[2] * dddij[2];
    cocg[c_id][3] += dddij[0] * dddij[1];
    cocg[c_id][4] += dddij[1] * dddij[2];
    cocg[c_id][5] += dddij[0] * dddij[2];
  }
}

/*----------------------------------------------------------------------------
 * Recompute the rhsv contribution from interior and extended cells
 *----------------------------------------------------------------------------*/

__global__ static void
_compute_rhsv_lsq_s_i_face(cs_lnum_t          size,
                           const cs_lnum_t   *cell_cells_idx,
                           const cs_lnum_t   *cell_cells,
                           const cs_real_3_t *cell_cen,
                           const cs_real_t   *c_weight,
                           cs_real_4_t       *rhsv)
{
  cs_lnum_t c_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (c_id < size) {
    cs_lnum_t s_id = cell_cells_idx[c_id];
    cs_lnum_t e_id = cell_cells_idx[c_id + 1];
    cs_real_t dc[3], ddc, _weight;
    cs_lnum_t c_id1;
    for (cs_lnum_t i = s_id; i < e_id; i++) {
      c_id1 = cell_cells[i];
      if (c_weight == NULL)
        _weight = 1;
      else
        _weight = 2. * c_weight[c_id1] / (c_weight[c_id] + c_weight[c_id1]);

      dc[0] = cell_cen[c_id1][0] - cell_cen[c_id][0];
      dc[1] = cell_cen[c_id1][1] - cell_cen[c_id][1];
      dc[2] = cell_cen[c_id1][2] - cell_cen[c_id][2];

      ddc = 1. / (dc[0] * dc[0] + dc[1] * dc[1] + dc[2] * dc[2]);
      _weight *= (rhsv[c_id1][3] - rhsv[c_id][3]) * ddc;

      rhsv[c_id][0] += dc[0] * _weight;
      rhsv[c_id][1] += dc[1] * _weight;
      rhsv[c_id][2] += dc[2] * _weight;
    }
  }
}

/*----------------------------------------------------------------------------
 * Compute rhsv from contributions at boundaries
 *----------------------------------------------------------------------------*/

__global__ static void
_compute_rhsv_lsq_s_b_face(cs_lnum_t         n_b_cells,
                           cs_real_t         inc,
                           const cs_lnum_t   b_cells[],
                           const cs_lnum_t   cell_b_faces_idx[],
                           const cs_lnum_t   cell_b_faces[],
                           const cs_real_t   b_face_normal[][3],
                           const cs_real_t   b_face_surf[],
                           const cs_real_t   b_dist[],
                           const cs_real_t   diipb[][3],
                           const cs_real_t   coefap[],
                           const cs_real_t   coefbp[],
                           cs_real_4_t      *restrict rhsv)
{
  cs_lnum_t b_c_idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (b_c_idx >= n_b_cells)
    return;

  cs_lnum_t c_id = b_cells[b_c_idx];

  cs_lnum_t s_id = cell_b_faces_idx[c_id];
  cs_lnum_t e_id = cell_b_faces_idx[c_id + 1];

  for (cs_lnum_t i = s_id; i < e_id; i++) {
    cs_lnum_t f_id = cell_b_faces[i];

    cs_real_t unddij = 1. / b_dist[f_id];
    cs_real_t umcbdd = (1. - coefbp[f_id]) * unddij;
    cs_real_t udbfs  = 1. / b_face_surf[f_id];
    cs_real_t dddij[3];
    dddij[0] = udbfs * b_face_normal[f_id][0] + umcbdd * diipb[f_id][0];
    dddij[1] = udbfs * b_face_normal[f_id][1] + umcbdd * diipb[f_id][1];
    dddij[2] = udbfs * b_face_normal[f_id][2] + umcbdd * diipb[f_id][2];

    unddij *= (coefap[f_id] * inc + (coefbp[f_id] - 1.) * rhsv[c_id][3]);

    rhsv[c_id][0] += dddij[0] * unddij;
    rhsv[c_id][1] += dddij[1] * unddij;
    rhsv[c_id][2] += dddij[2] * unddij;
  }
}

/*----------------------------------------------------------------------------
 * Initialize RHSV with pvar
 *----------------------------------------------------------------------------*/

__global__ static void
_init_rhsv(cs_lnum_t         size,
           cs_real_4_t      *restrict rhsv,
           const cs_real_t  *pvar)
{
  cs_lnum_t c_id = blockIdx.x * blockDim.x + threadIdx.x;

  if (c_id < size) {
    rhsv[c_id][0] = 0.0;
    rhsv[c_id][1] = 0.0;
    rhsv[c_id][2] = 0.0;
    rhsv[c_id][3] = pvar[c_id];
  }
}

/*----------------------------------------------------------------------------
 * Compute rhsv contributions from interior faces for strided case
 *----------------------------------------------------------------------------*/

template <unsigned int blocksize, cs_lnum_t stride>
__global__ static void
_compute_rhs_lsq_strided_i_face(cs_lnum_t             n_cells,
                                const cs_lnum_t      *restrict cell_cells_idx,
                                const cs_lnum_t      *restrict cell_cells,
                                const cs_lnum_t      *restrict cell_cells_e_idx,
                                const cs_lnum_t      *restrict cell_cells_e,
                                const cs_lnum_t      *restrict cell_i_faces,
                                const short int      *restrict cell_i_faces_sgn,
                                const cs_real_3_t    *restrict cell_f_cen,
                                const cs_real_t     (*restrict pvar)[stride],
                                const cs_real_t      *restrict weight,
                                const cs_real_t      *restrict c_weight,
                                cs_real_t           (*restrict rhs)[stride][3])
{
  cs_lnum_t c_id1 = blockIdx.x * blockDim.x + threadIdx.x;
  cs_lnum_t lidx = threadIdx.x;

  if (c_id1 >= n_cells) {
    return;
  }

  // size_t c_id1 = c_id / (3*3);
  // size_t i = (c_id / 3) % 3;
  // size_t j = c_id % 3;

  __shared__ cs_real_t _rhs[blocksize][stride][3];

  for (cs_lnum_t i = 0; i < stride; i++) {
    for (cs_lnum_t j = 0; j < 3; j++) {
      _rhs[lidx][i][j] = 0.0;
    }
  }

  auto pvar1 = pvar[c_id1];

  auto cell_f_cen1 = cell_f_cen[c_id1];

  /* Contribution from standard neighborhod */

  cs_lnum_t s_id = cell_cells_idx[c_id1];
  cs_lnum_t e_id = cell_cells_idx[c_id1 + 1];

  cs_real_t dc[stride], fctb[stride], lweight, ddc;

  for (cs_lnum_t idx = s_id; idx < e_id; idx++) {
    cs_lnum_t c_id2 = cell_cells[idx];

    auto cell_f_cen2 = cell_f_cen[c_id2];

    dc[0] = cell_f_cen2[0] - cell_f_cen1[0];
    dc[1] = cell_f_cen2[1] - cell_f_cen1[1];
    dc[2] = cell_f_cen2[2] - cell_f_cen1[2];

    ddc = 1./(dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

    if (c_weight == NULL) {
      lweight = 1.;
    }
    else {
      cs_lnum_t f_id = cell_i_faces[idx];
      cs_real_t pond, denom;
      pond = (cell_i_faces_sgn[idx] > 0) ? weight[f_id] : 1. - weight[f_id];
      denom = 1. / (        pond *c_weight[c_id1]
                    + (1. - pond)*c_weight[c_id2]);
      lweight = c_weight[c_id2] * denom;
    }

    auto pvar2 = pvar[c_id2];

    for (cs_lnum_t i = 0; i < stride; i++) {
      cs_real_t pfac = (pvar2[i] - pvar1[i]) * ddc;
      for (cs_lnum_t j = 0; j < 3; j++) {
        fctb[j] = dc[j] * pfac;
        _rhs[lidx][i][j] += lweight * fctb[j];
      }
    }

  }

  /* Contribution from standard neighborhod */

  if (cell_cells_e_idx != NULL) {
    s_id = cell_cells_e_idx[c_id1];
    e_id = cell_cells_e_idx[c_id1 + 1];

    for (cs_lnum_t idx = s_id; idx < e_id; idx++) {
      cs_lnum_t c_id2 = cell_cells[idx];

      auto cell_f_cen2 = cell_f_cen[c_id2];

      dc[0] = cell_f_cen2[0] - cell_f_cen1[0];
      dc[1] = cell_f_cen2[1] - cell_f_cen1[1];
      dc[2] = cell_f_cen2[2] - cell_f_cen1[2];

      ddc = 1./(dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]);

      auto pvar2 = pvar[c_id2];

      for (cs_lnum_t i = 0; i < stride; i++) {
        cs_real_t pfac = (pvar2[i] - pvar1[i]) * ddc;
        for (cs_lnum_t j = 0; j < 3; j++) {
          _rhs[lidx][i][j] += dc[j] * pfac;
        }
      }
    }
  }

  /* Copy from shared memory */

  for (cs_lnum_t i = 0; i < stride; i++) {
    for (cs_lnum_t j = 0; j < 3; j++) {
      rhs[c_id1][i][j] = _rhs[lidx][i][j];
    }
  }
}

/*----------------------------------------------------------------------------
 * Compute base rhsv contributions from boundary faces for strided case
 *----------------------------------------------------------------------------*/

template <unsigned int blocksize, cs_lnum_t stride>
__global__ static void
_compute_rhs_lsq_strided_b_face(cs_lnum_t             n_b_cells,
                                const int             inc,
                                const cs_lnum_t      *restrict cell_b_faces_idx,
                                const cs_lnum_t      *restrict cell_b_faces,
                                const cs_lnum_t      *restrict b_cells,
                                const cs_real_3_t    *restrict cell_cen,
                                const cs_real_3_t    *restrict b_face_cog,
                                const cs_real_t      *restrict b_dist,
                                const cs_real_t (*restrict coefav)[stride],
                                const cs_real_t (*restrict coefbv)[stride][stride],
                                const cs_real_t     (*restrict pvar)[stride],
                                const cs_cocg_6_t    *restrict cocgb,
                                cs_cocg_6_t          *restrict cocg,
                                cs_real_t           (*restrict rhs)[stride][3])
{
  cs_lnum_t c_idx = blockIdx.x * blockDim.x + threadIdx.x;
  cs_lnum_t lidx = threadIdx.x;

  if (c_idx >= n_b_cells) {
    return;
  }

  cs_lnum_t c_id = b_cells[c_idx];

  for (cs_lnum_t ll = 0; ll < 6; ll++)
    cocg[c_id][ll] = cocgb[c_idx][ll];

  cs_lnum_t s_id = cell_b_faces_idx[c_id];
  cs_lnum_t e_id = cell_b_faces_idx[c_id + 1];

  __shared__ cs_real_t _rhs[blocksize][3][3];

  for (cs_lnum_t i = 0; i < stride; i++){
    for (cs_lnum_t j = 0; j < 3; j++){
      _rhs[lidx][i][j] = rhs[c_id][i][j];
    }
  }

  auto pvar_c = pvar[c_id];

  for (cs_lnum_t idx = s_id; idx < e_id; idx++) { /* loop on boundary faces */

    cs_lnum_t f_id = cell_b_faces[idx];

    auto _coefav = coefav[f_id];
    auto _coefbv = coefbv[f_id];

    cs_real_t dif[3];
    for (cs_lnum_t ll = 0; ll < 3; ll++)
      dif[ll] = b_face_cog[f_id][ll] - cell_cen[c_id][ll];

    cs_real_t ddif = 1. / cs_math_3_square_norm_cuda(dif);

    cocg[c_id][0] += dif[0]*dif[0]*ddif;
    cocg[c_id][1] += dif[1]*dif[1]*ddif;
    cocg[c_id][2] += dif[2]*dif[2]*ddif;
    cocg[c_id][3] += dif[0]*dif[1]*ddif;
    cocg[c_id][4] += dif[1]*dif[2]*ddif;
    cocg[c_id][5] += dif[0]*dif[2]*ddif;

    cs_real_t var_f[stride];

    for (cs_lnum_t kk = 0; kk < stride; kk++) {
      var_f[kk] = _coefav[kk]*inc;
      for (cs_lnum_t ll = 0; ll < stride; ll++) {
        var_f[kk] += _coefbv[ll][kk] * pvar_c[ll];
      }

      cs_real_t pfac = (var_f[kk] - pvar_c[kk]) * ddif;

      for (cs_lnum_t ll = 0; ll < 3; ll++)
        _rhs[lidx][kk][ll] += dif[ll] * pfac;
    }

  } /* loop on boundary faces */

  _math_6_inv_cramer_sym_in_place_cuda(cocg[c_id]);

  for (cs_lnum_t i = 0; i < stride; i++){
    for (cs_lnum_t j = 0; j < 3; j++){
      rhs[c_id][i][j] = _rhs[lidx][i][j];
    }
  }
}

/*----------------------------------------------------------------------------
 * Compute the gradient for strided types from the LSQ covariance,
 * using on thread per value for coalescing
 *----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
__global__ static void
_compute_gradient_lsq_strided(cs_lnum_t          n_cells,
                              cs_real_t        (*restrict grad)[stride][3],
                              cs_cocg_6_t       *restrict cocg,
                              const cs_real_t  (*restrict rhs)[stride][3])
{
  size_t t_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (t_id >= n_cells*stride*3)
    return;

  size_t c_id = t_id / (stride*3);
  size_t i = (t_id / stride) % stride;
  size_t j = t_id % 3;

  auto cocg_temp = cocg[c_id];
  cs_real_t cocg_l[3];

  cocg_l[0] = cocg_temp[5];
  cocg_l[1] = cocg_temp[4];
  cocg_l[2] = cocg_temp[2];

  if (j == 0) {
    cocg_l[0] = cocg_temp[0];
    cocg_l[1] = cocg_temp[3];
    cocg_l[2] = cocg_temp[5];
  }

  if (j == 1) {
    cocg_l[0] = cocg_temp[3];
    cocg_l[1] = cocg_temp[1];
    cocg_l[2] = cocg_temp[4];
  }

  grad[c_id][i][j] =   rhs[c_id][i][0] * cocg_l[0]
                     + rhs[c_id][i][1] * cocg_l[1]
                     + rhs[c_id][i][2] * cocg_l[2];

}

/*----------------------------------------------------------------------------
 * Correct gradient with Neumann BC's at non-orthogonal boundary cells
 * for strided cases, using fixed-point algorimth
 *----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
__global__ static void
_correct_gradient_b_strided(const cs_lnum_t              n_b_cells,
                            const int                    n_c_iter_max,
                            const cs_real_t              c_eps,
                            const cs_real_t              epzero,
                            const cs_lnum_t    *restrict b_cells,
                            const cs_lnum_t    *restrict cell_b_faces_idx,
                            const cs_lnum_t    *restrict cell_b_faces,
                            const cs_real_3_t  *restrict b_face_cog,
                            const cs_real_3_t  *restrict cell_cen,
                            const cs_real_3_t  *restrict diipb,
                            const cs_real_t   (*restrict coefbv)[stride][stride],
                            cs_cocg_6_t        *restrict cocg,
                            cs_real_t         (*restrict grad)[stride][3])
{
  size_t t_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (t_id >= n_b_cells)
    return;

  cs_lnum_t c_id = b_cells[t_id];

  cs_lnum_t s_id = cell_b_faces_idx[c_id];
  cs_lnum_t e_id = cell_b_faces_idx[c_id+1];

  auto c_grad = grad[c_id];
  auto _cocg = cocg[c_id];
  auto _cell_cen = cell_cen[c_id];

  const cs_real_t eps_dvg = 1e-2;

  cs_real_t grad_0[stride][3], grad_i[stride][3];

  for (cs_lnum_t i = 0; i < stride; i++) {
    for (cs_lnum_t j = 0; j < 3; j++) {
      grad_0[i][j] = c_grad[i][j];
      grad_i[i][j] = c_grad[i][j];
    }
  }

  cs_real_t ref_norm = 0;
  for (cs_lnum_t kk = 0; kk < stride; kk++) {
    for (cs_lnum_t ll = 0; ll < 3; ll++)
      ref_norm += cs_math_abs_cuda(c_grad[kk][ll]);
  }

  cs_real_t c_norm = 0;
  int n_c_it;

  for (n_c_it = 0; n_c_it < n_c_iter_max; n_c_it++) {

    cs_real_t dif[3], grad_c[stride][3], var_ip_f[stride];

    cs_real_t rhs_c[stride][3];
    for (cs_lnum_t ll = 0; ll < stride; ll++) {
      rhs_c[ll][0] = 0;
      rhs_c[ll][1] = 0;
      rhs_c[ll][2] = 0;
    }

    for (cs_lnum_t fidx = s_id; fidx < e_id; fidx++) {
      cs_lnum_t f_id = cell_b_faces[fidx];

      for (cs_lnum_t ii = 0; ii < 3; ii++)
        dif[ii] = b_face_cog[f_id][ii] - _cell_cen[ii];

      cs_real_t ddif = 1. / cs_math_3_square_norm_cuda(dif);

      for (cs_lnum_t ll = 0; ll < stride; ll++) {
        var_ip_f[ll] = cs_math_3_dot_product_cuda(c_grad[ll], diipb[f_id]);
      }

      auto b = coefbv[f_id];

      for (cs_lnum_t kk = 0; kk < stride; kk++) {
        cs_real_t pfac = 0;
        for (cs_lnum_t ll = 0; ll < stride; ll++) {
          pfac += b[kk][ll] * var_ip_f[ll] * ddif;
        }

        for (cs_lnum_t ll = 0; ll < 3; ll++)
          rhs_c[kk][ll] += dif[ll] * pfac;
      }

    }

    for (cs_lnum_t i = 0; i < stride; i++){
      grad_c[i][0] =   rhs_c[i][0] * _cocg[0]
                     + rhs_c[i][1] * _cocg[3]
                     + rhs_c[i][2] * _cocg[5];

      grad_c[i][1] =   rhs_c[i][0] * _cocg[3]
                     + rhs_c[i][1] * _cocg[1]
                     + rhs_c[i][2] * _cocg[4];

      grad_c[i][2] =   rhs_c[i][0] * _cocg[5]
                     + rhs_c[i][1] * _cocg[4]
                     + rhs_c[i][2] * _cocg[2];
    }

    c_norm = 0.0;
    for (cs_lnum_t ii = 0; ii < stride; ii++) {
      for (cs_lnum_t jj = 0; jj < 3; jj++) {
        c_grad[ii][jj] = grad_0[ii][jj] + grad_c[ii][jj];
        c_norm += cs_math_abs_cuda(c_grad[ii][jj] - grad_i[ii][jj]);
        grad_i[ii][jj] = c_grad[ii][jj];
      }
    }

    if (c_norm < ref_norm * c_eps || c_norm < epzero)
        break;
  }

  for (cs_lnum_t ii = 0; ii < stride; ii++) {
    for (cs_lnum_t jj = 0; jj < 3; jj++) {
      grad[c_id][ii][jj] = c_grad[ii][jj];
    }
  }

  if (c_norm > eps_dvg * ref_norm) {
    for (cs_lnum_t ii = 0; ii < stride; ii++) {
      for (cs_lnum_t jj = 0; jj < 3; jj++) {
        grad[c_id][ii][jj] = grad_0[ii][jj];
      }
    }

    n_c_it *= -1;
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Semi-private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute cell gradient using least-squares reconstruction for non-orthogonal
 * meshes (nswrgp > 1).
 *
 * Optionally, a volume force generating a hydrostatic pressure component
 * may be accounted for.
 *
 * cocg is computed to account for variable B.C.'s (flux).
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   halo_type      <-- halo type (extended or not)
 *   recompute_cocg <-- flag to recompute cocg
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   coefap         <-- B.C. coefficients for boundary face normals
 *   coefbp         <-- B.C. coefficients for boundary face normals
 *   pvar           <-- variable
 *   c_weight       <-- weighted gradient coefficient variable,
 *                      or NULL
 *   cocg           <-> associated cell covariance array (on device)
 *   cocgb          <-> saved boundary cell covariance array (on device)
 *   grad           <-> gradient of pvar (halo prepared for periodicity
 *                      of rotation)
 *----------------------------------------------------------------------------*/

void
cs_gradient_scalar_lsq_cuda(const cs_mesh_t              *m,
                            const cs_mesh_quantities_t   *fvq,
                            cs_halo_type_t                halo_type,
                            bool                          recompute_cocg,
                            cs_real_t                     inc,
                            const cs_real_t              *coefap,
                            const cs_real_t              *coefbp,
                            const cs_real_t              *pvar,
                            const cs_real_t     *restrict c_weight,
                            cs_cocg_6_t         *restrict cocg,
                            cs_cocg_6_t         *restrict cocgb,
                            cs_real_3_t         *restrict grad)
{
  const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;

  const cs_lnum_t n_cells     = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces   = m->n_b_faces;

  static bool init_cocg = true;  /* A copy from CPU would suffice, or better,
                                    a separate device computation. */

  int device_id;
  cudaGetDevice(&device_id);

  cudaStream_t stream, stream1;
  cudaStreamCreate(&stream1);
  cudaStreamCreate(&stream);

  cs_real_4_t *rhsv;
  CS_CUDA_CHECK(cudaMalloc(&rhsv, n_cells_ext * sizeof(cs_real_4_t)));

  void *_pvar_d = NULL, *_coefa_d = NULL, *_coefb_d = NULL;
  const cs_real_t *pvar_d = NULL, *coefa_d = NULL, *coefb_d = NULL;

  cs_sync_or_copy_h2d(pvar, n_cells_ext, device_id, stream1,
                      &pvar_d, &_pvar_d);

  cs_sync_or_copy_h2d(coefap, n_b_faces, device_id, stream1,
                      &coefa_d, &_coefa_d);
  cs_sync_or_copy_h2d(coefbp, n_b_faces, device_id, stream1,
                      &coefb_d, &_coefb_d);

  unsigned int blocksize = 256;
  unsigned int gridsize_b
    = (unsigned int)ceil((double)m->n_b_cells / blocksize);
  unsigned int gridsize_bf
    = (unsigned int)ceil((double)m->n_b_faces / blocksize);
  unsigned int gridsize = (unsigned int)ceil((double)m->n_cells / blocksize);
  unsigned int gridsize_ext
    = (unsigned int)ceil((double)n_cells_ext / blocksize);

  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(m->b_face_cells);
  const cs_lnum_t *restrict cell_cells_idx
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(ma->cell_cells_idx);
  const cs_lnum_t *restrict cell_cells
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(ma->cell_cells);
  const cs_lnum_t *restrict cell_cells_e_idx
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(ma->cell_cells_e_idx);
  const cs_lnum_t *restrict cell_cells_e
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(ma->cell_cells_e);

  const cs_lnum_t *restrict b_cells
    = (cs_lnum_t *restrict)cs_get_device_ptr_const_pf(m->b_cells);
  const cs_lnum_t *restrict cell_b_faces_idx
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(ma->cell_b_faces_idx);
  const cs_lnum_t *restrict cell_b_faces
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(ma->cell_b_faces);

  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)cs_get_device_ptr_const_pf(fvq->cell_cen);
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)cs_get_device_ptr_const_pf(fvq->b_face_normal);
  const cs_real_t *restrict b_face_surf
    = (const cs_real_t *restrict)cs_get_device_ptr_const_pf(fvq->b_face_surf);
  const cs_real_t *restrict b_dist
    = (const cs_real_t *restrict)cs_get_device_ptr_const_pf(fvq->b_dist);
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)cs_get_device_ptr_const_pf(fvq->diipb);
  const cs_real_t *restrict weight
    = (const cs_real_t *restrict)cs_get_device_ptr_const_pf(fvq->weight);

  cudaStreamDestroy(stream1);
  cudaStreamSynchronize(0);

  _init_rhsv<<<gridsize_ext, blocksize, 0, stream>>>
    (n_cells_ext, rhsv, pvar_d);

  /* Reconstruct gradients using least squares for non-orthogonal meshes */
  /*---------------------------------------------------------------------*/

  /* Recompute cocg and at rhsv from interior cells */

  if (init_cocg) {

    _compute_cocg_rhsv_lsq_s_i_face<<<gridsize, blocksize, 0, stream>>>
      (n_cells, cocg, cell_cells_idx, cell_cells, cell_cen, rhsv, c_weight);

    /* Contribution from extended neighborhood */
    if (halo_type == CS_HALO_EXTENDED && cell_cells_e_idx != NULL)
      _compute_cocg_rhsv_lsq_s_i_face<<<gridsize, blocksize, 0, stream>>>
        (n_cells,
         cocg,
         cell_cells_e_idx,
         cell_cells_e,
         cell_cen,
         rhsv,
         c_weight);

    _save_cocgb<<<gridsize_b, blocksize, 0, stream>>>
      (m->n_b_cells, b_cells, cocg, cocgb);

    _compute_cocg_rhsv_lsq_s_b_face<<<gridsize_bf, blocksize, 0, stream>>>
      (m->n_b_cells,
       inc,
       b_cells,
       cell_b_faces_idx,
       cell_b_faces,
       b_face_normal,
       b_face_surf,
       b_dist,
       diipb,
       coefa_d,
       coefb_d,
       cocg,
       rhsv);

    /* Invert COCG at all cells */
    _compute_cocg_inv<<<gridsize, blocksize, 0, stream>>>
      (m->n_cells, NULL, cocg);

    init_cocg = false;

  }
  else {

    if (recompute_cocg) {

      _compute_cocg_from_cocgb<<<gridsize_b, blocksize, 0, stream>>>
        (m->n_b_cells, b_cells, cocg, cocgb);

      /* Recompute cocg and rhsv at boundaries*/
      _compute_cocg_rhsv_lsq_s_b_face<<<gridsize_bf, blocksize, 0, stream>>>
        (m->n_b_cells,
         inc,
         b_cells,
         cell_b_faces_idx,
         cell_b_faces,
         b_face_normal,
         b_face_surf,
         b_dist,
         diipb,
         coefa_d,
         coefb_d,
         cocg,
         rhsv);

      /* Invert COCG at boundary cells */
      _compute_cocg_inv<<<gridsize_b, blocksize, 0, stream>>>
        (m->n_b_cells, b_cells, cocg);

    }
    else
      _compute_rhsv_lsq_s_b_face<<<gridsize_bf, blocksize, 0, stream>>>
        (m->n_b_cells,
         inc,
         b_cells,
         cell_b_faces_idx,
         cell_b_faces,
         b_face_normal,
         b_face_surf,
         b_dist,
         diipb,
         coefa_d,
         coefb_d,
         rhsv);

    _compute_rhsv_lsq_s_i_face<<<gridsize, blocksize, 0, stream>>>
      (n_cells, cell_cells_idx, cell_cells, cell_cen, c_weight, rhsv);

    /* Contribution from extended neighborhood */
    if (halo_type == CS_HALO_EXTENDED && cell_cells_e_idx != NULL)
      _compute_rhsv_lsq_s_i_face<<<gridsize, blocksize, 0, stream>>>
        (n_cells,
         cell_cells_e_idx,
         cell_cells_e,
         cell_cen,
         c_weight,
         rhsv);
  }

  /* Compute gradient */
  /*------------------*/

  void *_grad_d = NULL;
  cs_real_3_t *grad_d = NULL;

  if (cs_check_device_ptr(grad) == CS_ALLOC_HOST) {
    size_t size = n_cells_ext * sizeof(cs_real_t) * 3;
    CS_CUDA_CHECK(cudaMalloc(&_grad_d, size));
    grad_d = (cs_real_3_t *)_grad_d;
  }
  else {
    grad_d = (cs_real_3_t *)cs_get_device_ptr((void *)grad);
  }

  _compute_gradient_lsq_s<<<gridsize, blocksize, 0, stream>>>
    (n_cells, grad_d, cocg, rhsv);

  cudaStreamSynchronize(stream);
  cudaStreamDestroy(stream);

  /* Sync to host */
  if (_grad_d != NULL) {
    size_t size = n_cells_ext * sizeof(cs_real_t) * 3;
    cs_cuda_copy_d2h(grad, grad_d, size);
  }
  else
    cs_sync_d2h(grad);

  if (_pvar_d != NULL)
    CS_CUDA_CHECK(cudaFree(_pvar_d));
  if (_coefa_d != NULL)
    CS_CUDA_CHECK(cudaFree(_coefa_d));
  if (_coefb_d != NULL)
    CS_CUDA_CHECK(cudaFree(_coefb_d));

  CS_CUDA_CHECK(cudaFree(rhsv));

  /* Synchronize halos (TODO move this to device) */

  if (m->halo != NULL) {
    cs_halo_sync_var_strided(m->halo, CS_HALO_STANDARD, (cs_real_t *)grad, 3);
    if (m->have_rotation_perio)
      cs_halo_perio_sync_var_vect(m->halo, CS_HALO_STANDARD, (cs_real_t *)grad,
                                  3);
  }
}

/*----------------------------------------------------------------------------*/
/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Semi-private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute cell gradient of a vector or tensor using least-squares
 * reconstruction for non-orthogonal meshes.
 *
 * template parameters:
 *   e2n           type of assembly algorithm used
 *   stride        3 for vectors, 6 for symmetric tensors
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   madj           <-- pointer to mesh adjacencies structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   halo_type      <-- halo type (extended or not)
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   n_c_iter_max   <-- maximum number of iterations for boundary correction
 *   c_eps          <-- relative tolerance for boundary correction
 *   coefav         <-- B.C. coefficients for boundary face normals
 *   coefbv         <-- B.C. coefficients for boundary face normals
 *   pvar           <-- variable
 *   c_weight       <-- weighted gradient coefficient variable, or NULL
 *   cocgb          <-- saved boundary cell covariance array (on device)
 *   cocg           <-> cocg covariance matrix for given cell
 *   grad           --> gradient of pvar (du_i/dx_j : grad[][i][j])
 *----------------------------------------------------------------------------*/

template <cs_lnum_t stride>
void
cs_gradient_strided_lsq_cuda(const cs_mesh_t               *m,
                             const cs_mesh_adjacencies_t   *madj,
                             const cs_mesh_quantities_t    *fvq,
                             const cs_halo_type_t           halo_type,
                             int                            inc,
                             int                            n_c_iter_max,
                             cs_real_t                      c_eps,
                             const cs_real_t (*restrict coefav)[stride],
                             const cs_real_t (*restrict coefbv)[stride][stride],
                             const cs_real_t (*restrict pvar)[stride],
                             const cs_real_t      *restrict c_weight,
                             const cs_cocg_6_t    *restrict cocgb,
                             cs_cocg_6_t          *restrict cocg,
                             cs_real_t (*restrict grad)[stride][3])
{
  using grad_t = cs_real_t[stride][3];

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  int device_id;
  cudaGetDevice(&device_id);

  cudaStream_t stream;
  cudaStreamCreate(&stream);

  cudaEvent_t e_start, e_init, e_i_faces, e_b_faces, e_gradient;
  cudaEvent_t e_b_correction, e_halo, e_stop;
  float msec = 0.0f;

  if (cs_glob_timer_kernels_flag > 0) {
    CS_CUDA_CHECK(cudaEventCreate(&e_start));
    CS_CUDA_CHECK(cudaEventCreate(&e_init));
    CS_CUDA_CHECK(cudaEventCreate(&e_i_faces));
    CS_CUDA_CHECK(cudaEventCreate(&e_b_faces));
    CS_CUDA_CHECK(cudaEventCreate(&e_gradient));
    CS_CUDA_CHECK(cudaEventCreate(&e_b_correction));
    CS_CUDA_CHECK(cudaEventCreate(&e_halo));
    CS_CUDA_CHECK(cudaEventCreate(&e_stop));

    // Record the start event
    CS_CUDA_CHECK(cudaEventRecord(e_start, stream));
  }

  grad_t *grad_d = NULL, *_grad_d = NULL;
  if (cs_check_device_ptr(grad) == CS_ALLOC_HOST) {
    size_t size = n_cells_ext * sizeof(cs_real_t) * stride * 3;
    CS_CUDA_CHECK(cudaMalloc(&_grad_d, size));
    grad_d = _grad_d;
  }
  else {
    grad_d = (grad_t *)cs_get_device_ptr((void *)grad);
  }

  void *_pvar_d = NULL, *_coefa_d = NULL, *_coefb_d = NULL;
  decltype(pvar) pvar_d = NULL, coefa_d = NULL;
  decltype(coefbv) coefb_d = NULL;

  using real_3_t = cs_real_t[3];

  cs_sync_or_copy_h2d(pvar, n_cells_ext*stride, device_id, stream,
                      &pvar_d, &_pvar_d);

  cs_sync_or_copy_h2d(coefav, n_b_faces*stride, device_id, stream,
                      &coefa_d, &_coefa_d);
  cs_sync_or_copy_h2d(coefbv, n_b_faces*stride*stride, device_id, stream,
                      &coefb_d, &_coefb_d);

  const cs_lnum_t *restrict b_face_cells
    = cs_get_device_ptr_const_pf(m->b_face_cells);
  const cs_lnum_t *restrict b_cells
    = cs_get_device_ptr_const_pf(m->b_cells);
  const cs_lnum_t *restrict cell_cells_idx
    = cs_get_device_ptr_const_pf(madj->cell_cells_idx);
  const cs_lnum_t *restrict cell_cells
    = cs_get_device_ptr_const_pf(madj->cell_cells);

  const cs_lnum_t *restrict cell_cells_e_idx = NULL;
  const cs_lnum_t *restrict cell_cells_e = NULL;
  if (halo_type == CS_HALO_EXTENDED) {
    cell_cells_e_idx = cs_get_device_ptr_const_pf(madj->cell_cells_e_idx);
    cell_cells_e = cs_get_device_ptr_const_pf(madj->cell_cells_e);
  }

  const cs_lnum_t *restrict cell_b_faces_idx
    = cs_get_device_ptr_const_pf(madj->cell_b_faces_idx);
  const cs_lnum_t *restrict cell_b_faces
    = cs_get_device_ptr_const_pf(madj->cell_b_faces);
  const cs_lnum_t *restrict cell_i_faces
    = cs_get_device_ptr_const_pf(madj->cell_i_faces);
  const short int *restrict cell_i_faces_sgn
    = cs_get_device_ptr_const_pf(madj->cell_i_faces_sgn);

  const cs_real_3_t *restrict cell_f_cen
    = cs_get_device_ptr_const_pf((cs_real_3_t *)fvq->cell_f_cen);
  const cs_real_t *restrict weight
    = cs_get_device_ptr_const_pf(fvq->weight);
  const cs_real_t *restrict b_dist
    = cs_get_device_ptr_const_pf(fvq->b_dist);
  const cs_real_3_t *restrict b_f_face_cog
    = cs_get_device_ptr_const_pf((cs_real_3_t *)fvq->b_f_face_cog);
  const cs_real_3_t *restrict diipb
    = cs_get_device_ptr_const_pf((cs_real_3_t *)fvq->diipb);

  decltype(grad) rhs_d;
  CS_CUDA_CHECK(cudaMalloc(&rhs_d, n_cells_ext * sizeof(cs_real_t)*stride*3));
  cudaMemset(rhs_d, 0, n_cells_ext*sizeof(grad));

  // cs_cuda_copy_h2d(grad_d, grad, sizeof(cs_real_t) * n_cells * stride * 3);

  if (cs_glob_timer_kernels_flag > 0)
    CS_CUDA_CHECK(cudaEventRecord(e_init, stream));

  const unsigned int blocksize = 256;
  int gridsize;

  gridsize = cs_cuda_grid_size(n_cells, blocksize);
  _compute_rhs_lsq_strided_i_face<blocksize><<<gridsize, blocksize, 0, stream>>>
    (n_cells,
     cell_cells_idx,
     cell_cells,
     cell_cells_e_idx,
     cell_cells_e,
     cell_i_faces,
     cell_i_faces_sgn,
     cell_f_cen,
     pvar_d,
     weight,
     c_weight,
     rhs_d);

  if (cs_glob_timer_kernels_flag > 0)
    CS_CUDA_CHECK(cudaEventRecord(e_i_faces, stream));

  gridsize = cs_cuda_grid_size(m->n_b_cells, blocksize);
  _compute_rhs_lsq_strided_b_face<blocksize><<<gridsize, blocksize, 0, stream>>>
    (m->n_b_cells,
     inc,
     cell_b_faces_idx,
     cell_b_faces,
     b_cells,
     cell_f_cen,
     b_f_face_cog,
     b_dist,
     coefa_d,
     coefb_d,
     pvar_d,
     cocgb,
     cocg,
     rhs_d);

  if (cs_glob_timer_kernels_flag > 0)
    CS_CUDA_CHECK(cudaEventRecord(e_b_faces, stream));

  gridsize = cs_cuda_grid_size(n_cells*stride*3, blocksize);
  _compute_gradient_lsq_strided<<<gridsize, blocksize, 0, stream>>>
    (n_cells, grad_d, cocg, rhs_d);

  if (cs_glob_timer_kernels_flag > 0)
    CS_CUDA_CHECK(cudaEventRecord(e_gradient, stream));

  gridsize = cs_cuda_grid_size(m->n_b_cells, blocksize);
  _correct_gradient_b_strided<stride><<<gridsize, blocksize, 0, stream>>>
    (m->n_b_cells,
     n_c_iter_max,
     c_eps,
     cs_math_epzero,
     b_cells,
     cell_b_faces_idx,
     cell_b_faces,
     b_f_face_cog,
     cell_f_cen,
     diipb,
     coefb_d,
     cocg,
     grad_d);

  if (cs_glob_timer_kernels_flag > 0)
    CS_CUDA_CHECK(cudaEventRecord(e_b_correction, stream));

  cudaStreamSynchronize(stream);

  cs_sync_strided_gradient_halo_d(m, halo_type, grad_d);

  if (cs_glob_timer_kernels_flag > 0)
    CS_CUDA_CHECK(cudaEventRecord(e_halo, stream));

  /* Sync to host */
  if (grad_d != NULL) {
    size_t size = n_cells * sizeof(cs_real_t) * stride * 3;
    cs_cuda_copy_d2h(grad, grad_d, size);
  }
  else
    cs_sync_d2h(grad);

  if (cs_glob_timer_kernels_flag > 0) {
    CS_CUDA_CHECK(cudaEventRecord(e_stop, stream));
    CS_CUDA_CHECK(cudaEventSynchronize(e_stop));
  }

  cudaStreamSynchronize(stream);
  cudaStreamDestroy(stream);

  if (cs_glob_timer_kernels_flag > 0) {
    printf("%d: %s<%d>", cs_glob_rank_id, __func__, stride);

    msec = 0.0f;
    CS_CUDA_CHECK(cudaEventElapsedTime(&msec, e_start, e_init));
    printf(", init = %.0f", msec*1000.f);

    msec = 0.0f;
    CS_CUDA_CHECK(cudaEventElapsedTime(&msec, e_init, e_i_faces));
    printf(", i_faces = %.0f", msec*1000.f);

    msec = 0.0f;
    CS_CUDA_CHECK(cudaEventElapsedTime(&msec, e_i_faces, e_b_faces));
    printf(", b_faces = %.0f", msec*1000.f);

    msec = 0.0f;
    CS_CUDA_CHECK(cudaEventElapsedTime(&msec, e_b_faces, e_gradient));
    printf(", gradient = %.0f", msec*1000.f);

    msec = 0.0f;
    CS_CUDA_CHECK(cudaEventElapsedTime(&msec, e_gradient, e_b_correction));
    printf(", b_correction = %.0f", msec*1000.f);

    msec = 0.0f;
    CS_CUDA_CHECK(cudaEventElapsedTime(&msec, e_b_correction, e_halo));
    printf(", halo = %.0f", msec*1000.f);

    msec = 0.0f;
    CS_CUDA_CHECK(cudaEventElapsedTime(&msec, e_start, e_stop));
    printf(", total = %.0f\n", msec*1000.f);

    CS_CUDA_CHECK(cudaEventDestroy(e_start));
    CS_CUDA_CHECK(cudaEventDestroy(e_init));
    CS_CUDA_CHECK(cudaEventDestroy(e_i_faces));
    CS_CUDA_CHECK(cudaEventDestroy(e_b_faces));
    CS_CUDA_CHECK(cudaEventDestroy(e_gradient));
    CS_CUDA_CHECK(cudaEventDestroy(e_b_correction));
    CS_CUDA_CHECK(cudaEventDestroy(e_halo));
    CS_CUDA_CHECK(cudaEventDestroy(e_stop));
  }

  CS_CUDA_CHECK(cudaFree(rhs_d));

  if (_pvar_d != NULL)
    CS_CUDA_CHECK(cudaFree(_pvar_d));
  if (_coefa_d != NULL)
    CS_CUDA_CHECK(cudaFree(_coefa_d));
  if (_coefb_d != NULL)
    CS_CUDA_CHECK(cudaFree(_coefb_d));

  if (_grad_d != NULL)
    CS_CUDA_CHECK(cudaFree(_grad_d));
}

INSTANTIATE_LSQ(cs_gradient_strided_lsq_cuda, 3);
INSTANTIATE_LSQ(cs_gradient_strided_lsq_cuda, 6);

/*----------------------------------------------------------------------------*/
