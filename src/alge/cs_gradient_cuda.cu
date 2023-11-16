/*============================================================================
 * Gradient reconstruction, CUDA implementations.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
#include "cs_gradient_cuda.cuh"

#include "cs_gradient.h"
#include "cs_gradient_lsq_vector.cuh"
#include "cs_gradient_lsq_vector_gather.cuh"
#include "cs_gradient_lsq_vector_gather_v2.cuh"
#include "cs_gradient_lsq_vector_gather_v3.cuh"
#include "cs_gradient_lsq_vector_v2.cuh"
#include "cs_gradient_lsq_vector_v3.cuh"
#include "cs_gradient_priv.h"
#include "cs_reconstruct_vector_gradient_gather.cuh"
#include "cs_reconstruct_vector_gradient_gather_v2.cuh"
#include "cs_reconstruct_vector_gradient_scatter_v2.cuh"

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

/*----------------------------------------------------------------------------
 * Recompute cocg at boundaries, using saved cocgb
 *----------------------------------------------------------------------------*/

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

__global__ static void
_init_rhs_v3(cs_lnum_t         size,
           double3      *restrict rhs)
{
  cs_lnum_t c_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (c_id >= size)
    return;

  rhs[c_id] = make_double3(0.0, 0.0, 0.0);
}

__global__ static void
_compute_gradient_lsq_v_v3(cs_lnum_t           size,
                        cs_real_33_t        *restrict gradv,
                        cs_real_33_t        *restrict rhs,
                        cs_cocg_6_t         *restrict cocg)
{
  size_t c_id = blockIdx.x * blockDim.x + threadIdx.x;
  if (c_id >= size)
    return;
  auto& gradc = gradv[c_id];
  auto& rhsc = rhs[c_id];
  auto cocgc = cocg[c_id];
  for(cs_lnum_t i = 0; i < 3; i++){
    auto& gradci = gradc[i];
    auto rhsci = rhsc[i];
    gradci[0] =   rhsci[0] * cocgc[0]
                          + rhsci[1] * cocgc[3]
                          + rhsci[2] * cocgc[5];

    gradci[1] =   rhsci[0] * cocgc[3]
                        + rhsci[1] * cocgc[1]
                        + rhsci[2] * cocgc[4];

    gradci[2] =   rhsci[0] * cocgc[5]
                        + rhsci[1] * cocgc[4]
                        + rhsci[2] * cocgc[2];
  }
}

__global__ static void
_compute_gradient_lsq_b_v(cs_lnum_t         size,
                          cs_lnum_t         n_b_cells,
                          cs_lnum_t         *restrict b_cells,
                          cs_real_33_t        *restrict gradv,
                          cs_real_33_t        *restrict rhs,
                          cs_cocg_6_t         *restrict cocg,
                          cs_real_3_t *restrict b_face_normal,
                          cs_lnum_t *restrict cell_b_faces,
                          cs_lnum_t *restrict cell_b_faces_idx)
{
  size_t c_id = blockIdx.x * blockDim.x + threadIdx.x;

  cs_lnum_t _33_9_idx[9][2];
  int nn = 0;
  for (int ll = 0; ll < 3; ll++) {
    for (int mm = 0; mm < 3; mm++) {
      _33_9_idx[nn][0] = ll;
      _33_9_idx[nn][1] = mm;
      nn++;
    }
  }

   /* Loop on boundary cells */
  cs_lnum_t c_id1 = b_cells[c_id];
  cs_real_t cocgb[3][3], cocgb_v[45], rhsb_v[9], x[9];

  cocgb[0][0] = cocg[c_id][0];
  cocgb[0][1] = cocg[c_id][3];
  cocgb[0][2] = cocg[c_id][5];
  cocgb[1][0] = cocg[c_id][3];
  cocgb[1][1] = cocg[c_id][1];
  cocgb[1][2] = cocg[c_id][4];
  cocgb[2][0] = cocg[c_id][5];
  cocgb[2][1] = cocg[c_id][4];
  cocgb[2][2] = cocg[c_id][2];

  cs_lnum_t s_id = cell_b_faces_idx[c_id];
  cs_lnum_t e_id = cell_b_faces_idx[c_id+1];
  cs_lnum_t f_id;
  cs_real_3_t normal;
  cs_real_t norm, inverse_norm;

  for (cs_lnum_t index = s_id; index < e_id; index++) {

    f_id = cell_b_faces[index];

    /* Normal is vector 0 if the b_face_normal norm is too small */
    norm = sqrt(b_face_normal[index][0]*b_face_normal[index][0]
            + b_face_normal[index][1]*b_face_normal[index][1]
            + b_face_normal[index][2]*b_face_normal[index][2]);

    inverse_norm = 1. / norm;

    normal[0] = inverse_norm * b_face_normal[index][0];
    normal[1] = inverse_norm * b_face_normal[index][1];
    normal[2] = inverse_norm * b_face_normal[index][2];

    for (cs_lnum_t ii = 0; ii < 3; ii++) {
      for (cs_lnum_t jj = 0; jj < 3; jj++)
        cocgb[ii][jj] += normal[ii] * normal[jj];
    }

  }

}

/*----------------------------------------------------------------------------
 * Synchronize of copy a T type array from the host to a device.
 *
 * parameters:
 *   val_h          <-- pointer to host data
 *   n_vals         <-- number of data values
 *   device_id      <-- associated device id
 *   stream         <-- associated stream (for async prefetch only)
 *   val_d          --> matching pointer on device
 *   buf_d          --> matching allocation pointer on device (should be freed
 *                      after use if non-NULL)
 *----------------------------------------------------------------------------*/

template <typename T>
static void
_sync_or_copy_real_h2d(const  T   *val_h,
                       cs_lnum_t           n_vals,
                       int                 device_id,
                       cudaStream_t        stream,
                       const T   **val_d,
                       void              **buf_d)
{
  const T  *_val_d = NULL;
  void             *_buf_d = NULL;

  cs_alloc_mode_t alloc_mode = cs_check_device_ptr(val_h);
  size_t size = n_vals * sizeof(T);

  if (alloc_mode == CS_ALLOC_HOST) {
    CS_CUDA_CHECK(cudaMalloc(&_buf_d, size));
    cs_cuda_copy_h2d(_buf_d, val_h, size);
    _val_d = (const T *)_buf_d;
  }
  else {
    _val_d = (const T *)cs_get_device_ptr((void *)val_h);

    if (alloc_mode == CS_ALLOC_HOST_DEVICE_SHARED)
      cudaMemPrefetchAsync(val_h, size, device_id, stream);
    else
      cs_sync_h2d(val_h);
  }

  *val_d = _val_d;
  *buf_d = _buf_d;
}

/* Compute gridsize*/

unsigned int 
get_gridsize(unsigned int size, unsigned int blocksize){
  unsigned int gridsize = (unsigned int)ceil((double)size / blocksize);

  return gridsize;
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

extern "C" void
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

  _sync_or_copy_real_h2d(pvar, n_cells_ext, device_id, stream1,
                         &pvar_d, &_pvar_d);

  _sync_or_copy_real_h2d(coefap, n_b_faces, device_id, stream1,
                         &coefa_d, &_coefa_d);
  _sync_or_copy_real_h2d(coefbp, n_b_faces, device_id, stream1,
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

  /* Synchronize halos */

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
 *   madj           <-- pointer to mesh adjacencies structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   halo_type      <-- halo type (extended or not)
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   coefav         <-- B.C. coefficients for boundary face normals
 *   coefbv         <-- B.C. coefficients for boundary face normals
 *   pvar           <-- variable
 *   gradv          --> gradient of pvar (du_i/dx_j : gradv[][i][j])
 *----------------------------------------------------------------------------*/
extern "C" void
cs_lsq_vector_gradient_cuda(const cs_mesh_t               *m,
                     const cs_mesh_adjacencies_t   *madj,
                     const cs_mesh_quantities_t    *fvq,
                     const cs_halo_type_t           halo_type,
                     const int                      inc,
                     const cs_real_3_t    *restrict coefav,
                     const cs_real_33_t   *restrict coefbv,
                     const cs_real_3_t    *restrict pvar,
                     const cs_real_t      *restrict c_weight,
                     cs_cocg_6_t          *restrict cocg,
                     cs_cocg_6_t          *restrict cocgb,
                     cs_real_33_t         *restrict gradv,
                     cs_real_33_t         *restrict rhs)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces   = m->n_b_faces;
  const cs_lnum_t n_i_faces   = m->n_i_faces;


  int device_id;
  cudaGetDevice(&device_id);

  cudaStream_t stream;
  cudaStreamCreate(&stream);

  cudaEvent_t start, mem_h2d, init, i_faces, halo, b_faces, gradient, stop;
  float msec = 0.0f, msecTotal = 0.0f;
  CS_CUDA_CHECK(cudaEventCreate(&start));
  CS_CUDA_CHECK(cudaEventCreate(&mem_h2d));
  CS_CUDA_CHECK(cudaEventCreate(&init));
  CS_CUDA_CHECK(cudaEventCreate(&i_faces));
  CS_CUDA_CHECK(cudaEventCreate(&halo));
  CS_CUDA_CHECK(cudaEventCreate(&b_faces));
  CS_CUDA_CHECK(cudaEventCreate(&gradient));
  CS_CUDA_CHECK(cudaEventCreate(&stop));

  // Record the start event
  CS_CUDA_CHECK(cudaEventRecord(start, stream));

  cs_real_33_t *rhs_d;
  CS_CUDA_CHECK(cudaMalloc(&rhs_d, n_cells_ext * sizeof(cs_real_33_t)));


  cs_real_33_t *grad_d = NULL;
  CS_CUDA_CHECK(cudaMalloc(&grad_d, n_cells * sizeof(cs_real_33_t)));

  void *_pvar_d = NULL, *_coefa_d = NULL, *_coefb_d = NULL,
  *_cell_cells_idx_d = NULL;
  const cs_real_3_t *pvar_d = NULL, *coefa_d = NULL;
  const cs_real_33_t *coefb_d = NULL;
  const cs_lnum_t *cell_cells_idx_d = NULL;

  // cs_cuda_copy_h2d(rhs_d, rhs, n_cells * sizeof(cs_real_33_t));

  unsigned int blocksize = 256;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)cs_get_device_ptr_const_pf(m->i_face_cells);
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(m->b_face_cells);
  const cs_lnum_t *restrict b_cells
    = (cs_lnum_t *restrict)cs_get_device_ptr_const_pf(m->b_cells);
  const cs_lnum_t *restrict cell_cells_idx
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(madj->cell_cells_idx);
  const cs_lnum_t *restrict cell_cells_lst
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(m->cell_cells_lst);
  const cs_lnum_t *restrict cell_b_faces_idx
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(madj->cell_b_faces_idx);
  const cs_lnum_t *restrict cell_b_faces
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(madj->cell_b_faces);
  const cs_lnum_t *restrict cell_i_faces
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(madj->cell_i_faces);
  const short int *restrict cell_i_faces_sgn
    = (const short int *restrict)cs_get_device_ptr_const_pf(madj->cell_i_faces_sgn);
  const int n_i_groups = m->i_face_numbering->n_groups;
  const int n_i_threads = m->i_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
  const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

  const cs_lnum_t *restrict cell_cells
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(madj->cell_cells);
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)cs_get_device_ptr_const_pf(fvq->cell_cen);
  const cs_real_3_t *restrict cell_f_cen
    = (const cs_real_3_t *restrict)cs_get_device_ptr_const_pf(fvq->cell_f_cen);
  const cs_real_t *restrict weight
    = (const cs_real_t *restrict)cs_get_device_ptr_const_pf(fvq->weight);
  const cs_real_t *restrict b_dist
    = (const cs_real_t *restrict)cs_get_device_ptr_const_pf(fvq->b_dist);
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)cs_get_device_ptr_const_pf(fvq->b_face_normal);
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)cs_get_device_ptr_const_pf(fvq->b_f_face_cog);

  const cs_real_t *restrict cell_f_cen_1d
    = (const cs_real_t *restrict)cs_get_device_ptr_const_pf(fvq->cell_f_cen);
  const cs_lnum_t *restrict i_face_cells_1d
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(m->i_face_cells);

  cs_lnum_t stride = 3;

  // printf("n_i_thread:%d\tn_i_groups:%d\tn_cells%d\n", n_i_threads, n_i_groups, n_cells);

  _sync_or_copy_real_h2d(pvar, n_cells_ext, device_id, stream,
                         &pvar_d, &_pvar_d);

  _sync_or_copy_real_h2d(coefav, n_b_faces, device_id, stream,
                         &coefa_d, &_coefa_d);
  _sync_or_copy_real_h2d(coefbv, n_b_faces, device_id, stream,
                         &coefb_d, &_coefb_d);

  CS_CUDA_CHECK(cudaEventRecord(mem_h2d, stream));

  // _init_rhs<<<get_gridsize(n_cells_ext, blocksize), blocksize, 0, stream>>>
  //   (n_cells_ext,
  //    rhs_d);
  cudaMemset(rhs_d, 0, n_cells_ext*sizeof(cs_real_33_t));

  // _init_rhs_v2<<<get_gridsize(n_cells_ext*3*3, blocksize), blocksize, 0, stream>>>
  //   (n_cells_ext*3*3,
  //    rhs_d);

  // _init_rhs_v3<<<get_gridsize(n_cells_ext*3, blocksize), blocksize, 0, stream>>>
  //   (n_cells_ext*3,
  //    rhs_d);

  CS_CUDA_CHECK(cudaEventRecord(init, stream));
	
  
  // _compute_rhs_lsq_v_i_face_v0<<<get_gridsize(n_i_faces, blocksize), blocksize, 0, stream>>>
  //     (n_i_faces,
  //      i_face_cells, 
  //      cell_f_cen, 
  //      rhs_d, 
  //      pvar_d, 
  //      weight, 
  //      c_weight);

  // _compute_rhs_lsq_v_i_face_cf<<<get_gridsize(n_i_faces, blocksize), blocksize, 0, stream>>>
  //     (n_i_faces,
  //      i_face_cells,
  //      cell_f_cen,
  //      rhs_d,
  //      pvar_d,
  //      weight,
  //      c_weight);
  // _compute_rhs_lsq_v_i_face<<<get_gridsize(n_i_faces, blocksize), blocksize, 0, stream>>>
  //     (n_i_faces,
  //      i_face_cells,
  //      cell_f_cen,
  //      rhs_d,
  //      pvar_d,
  //      weight,
  //      c_weight);

  _compute_rhs_lsq_v_i_face_v2cf<<<get_gridsize(n_i_faces, blocksize), blocksize, 0, stream>>>
     (n_i_faces,
      i_face_cells, 
      cell_f_cen, 
      rhs_d, 
      pvar_d, 
      weight, 
      c_weight);

  // _compute_rhs_lsq_v_i_face_v3<<<get_gridsize(n_i_faces*3*3, blocksize), blocksize, 0, stream>>>
  //     (n_i_faces*3*3,
  //      i_face_cells, 
  //      cell_f_cen, 
  //      rhs_d, 
  //      pvar_d, 
  //      weight, 
  //      c_weight);
  // assert(cell_cells_idx);
  // assert(cell_cells);
  // assert(cell_f_cen);
  // assert(rhs_d);
  // assert(pvar_d);
  // assert(weight);
  // _compute_rhs_lsq_v_i_face_gather<<<get_gridsize(n_cells, blocksize), blocksize, 0, stream>>>
  //     (n_cells,
  //      cell_cells_idx,
  //      cell_cells,
  //      cell_i_faces,
  //      cell_i_faces_sgn,
  //      cell_f_cen, 
  //      rhs_d, 
  //      pvar_d, 
  //      weight, 
  //      c_weight);

  // _compute_rhs_lsq_v_i_face_gather_v2<<<get_gridsize(n_cells, blocksize), blocksize, 0, stream>>>
  //    (n_cells,
  //     cell_cells_idx,
  //     cell_cells,
  //     cell_i_faces,
  //     cell_i_faces_sgn,
  //     cell_f_cen, 
  //     rhs_d, 
  //     pvar_d, 
  //     weight, 
  //     c_weight);
  
  // _compute_rhs_lsq_v_i_face_gather_v4<<<get_gridsize(n_cells, blocksize), blocksize, 0, stream>>>
  //     (n_cells,
  //      cell_cells_idx,
  //      cell_cells,
  //      cell_i_faces,
  //      cell_i_faces_sgn,
  //      cell_f_cen, 
  //      rhs_d, 
  //      pvar_d, 
  //      weight, 
  //      c_weight);

  CS_CUDA_CHECK(cudaEventRecord(i_faces, stream));

  if(halo_type == CS_HALO_EXTENDED && cell_cells_idx != NULL){

    _compute_rhs_lsq_v_b_neighbor<<<get_gridsize(n_cells, blocksize), blocksize, 0, stream>>>
      (n_cells, 
       cell_cells_idx, 
       cell_cells, 
       cell_f_cen, 
       rhs_d, 
       pvar_d);
  }
  CS_CUDA_CHECK(cudaEventRecord(halo, stream));

  // _compute_rhs_lsq_v_b_face<<<get_gridsize(m->n_b_cells, blocksize), blocksize, 0, stream>>>
  //     (m->n_b_faces,
  //      b_face_cells,
  //      cell_f_cen,
  //      b_face_normal,
  //      rhs_d,
  //      pvar_d,
  //      b_dist,
  //      coefb_d,
  //      coefa_d,
  //      inc);

  // _compute_rhs_lsq_v_b_face_gather_stride<3, cs_real_3_t, cs_real_33_t><<<get_gridsize(m->n_b_cells, blocksize), blocksize, 0, stream>>>
  //     (m->n_b_cells,
  //      cell_b_faces_idx,
  //      cell_b_faces,
  //      b_cells,
  //      b_face_cog,
  //      cell_cen, 
  //      rhs_d, 
  //      pvar_d, 
  //      coefb_d, 
  //      coefa_d,
  //      cocg,
  //      cocgb, 
  //      inc);
    
  _compute_rhs_lsq_v_b_face_gather_v3<<<get_gridsize(m->n_b_cells, blocksize), blocksize, 0, stream>>>
     (m->n_b_cells,
      cell_b_faces_idx,
      cell_b_faces,
      b_cells,
      b_face_normal, 
      rhs_d, 
      pvar_d, 
      b_dist, 
      coefb_d, 
      coefa_d, 
      inc);

  // _compute_rhs_lsq_v_b_face_v2<<<get_gridsize(m->n_b_cells, blocksize), blocksize, 0, stream>>>
  //     (m->n_b_faces,
  //      b_face_cells, 
  //      cell_f_cen, 
  //      b_face_normal, 
  //      rhs_d, 
  //      pvar_d, 
  //      b_dist, 
  //      coefb_d, 
  //      coefa_d, 
  //      inc);

  CS_CUDA_CHECK(cudaEventRecord(b_faces, stream));


  // if (rhs_d != NULL) {
  //   size_t size = n_cells_ext * sizeof(cs_real_t) * 3 * 3;
  //   cs_cuda_copy_d2h(rhs, rhs_d, size);
  // }
  // else
  //   cs_sync_d2h(rhs);

  // /* Compute gradient */
  // /*------------------*/

  // _compute_gradient_lsq_v<<<get_gridsize(n_cells, blocksize), blocksize, 0, stream>>>
  //   (n_cells,
  //    grad_d, 
  //    rhs_d, 
  //    cocg);

  // _compute_gradient_lsq_v_v4<<<get_gridsize(n_cells, blocksize), blocksize, 0, stream>>>
  //   (n_cells,
  //    grad_d, 
  //    rhs_d, 
  //    cocg);


  // _compute_gradient_lsq_v_v5<<<get_gridsize(n_cells*3*3, blocksize), blocksize, 0, stream>>>
  //   (n_cells*3*3,
  //    gradv_d, 
  //    rhs_d, 
  //    cocg);

  _compute_gradient_lsq_v_v6<<<get_gridsize(n_cells*3*3, blocksize), blocksize, 0, stream>>>
    (n_cells*3*3,
     grad_d, 
     rhs_d, 
     cocg);

  CS_CUDA_CHECK(cudaEventRecord(gradient, stream));

  // /* Sync to host */
  if (grad_d != NULL) {
    size_t size = n_cells * sizeof(cs_real_t) * 3 * 3;
    cs_cuda_copy_d2h(gradv, grad_d, size);
  }
  else
    cs_sync_d2h(gradv);

  CS_CUDA_CHECK(cudaEventRecord(stop, stream));
  CS_CUDA_CHECK(cudaEventSynchronize(stop));

  cudaStreamSynchronize(stream);
  cudaStreamDestroy(stream);

  printf("lsq Kernels :");
  msec = 0.0f;
	CS_CUDA_CHECK(cudaEventElapsedTime(&msec, mem_h2d, init));
  printf("Kernels execution time in us: \t");
  printf("Init = %f\t", msec*1000.f);

  msec = 0.0f;
	CS_CUDA_CHECK(cudaEventElapsedTime(&msec, init, i_faces));
  printf("I_faces = %f\t", msec*1000.f);

  msec = 0.0f;
	CS_CUDA_CHECK(cudaEventElapsedTime(&msec, i_faces, halo));
  printf("Halo = %f\t", msec*1000.f);

  msec = 0.0f;
  CS_CUDA_CHECK(cudaEventElapsedTime(&msec, halo, b_faces));
  printf("B_faces = %f\t", msec*1000.f);

  msec = 0.0f;
  CS_CUDA_CHECK(cudaEventElapsedTime(&msec, b_faces, gradient));
  printf("Gradient = %f\t", msec*1000.f);

  msec = 0.0f;
  CS_CUDA_CHECK(cudaEventElapsedTime(&msec, mem_h2d, gradient));
  printf("Total kernel = %f\t", msec*1000.f);

  msec = 0.0f;
  CS_CUDA_CHECK(cudaEventElapsedTime(&msec, start, stop));
  printf("Total = %f\t", msec*1000.f);

  printf("\n");


  if (_pvar_d != NULL)
    CS_CUDA_CHECK(cudaFree(_pvar_d));
  if (_coefa_d != NULL)
    CS_CUDA_CHECK(cudaFree(_coefa_d));
  if (_coefb_d != NULL)
    CS_CUDA_CHECK(cudaFree(_coefb_d));

  CS_CUDA_CHECK(cudaFree(rhs_d));
  CS_CUDA_CHECK(cudaFree(grad_d));
  
}


__global__ static void
_compute_reconstruct_v_i_face(cs_lnum_t            size,
                          const cs_lnum_t      *i_group_index,
                          const cs_lnum_2_t      *i_face_cells,
                          const cs_real_3_t    *pvar,
                          const cs_real_t         *weight,
                          const cs_real_t      *c_weight,
                          const cs_real_33_t        *restrict r_grad,
                          cs_real_33_t        *restrict grad,
                          const cs_real_3_t *restrict dofij,
                          const cs_real_3_t *restrict i_f_face_normal)
{
  cs_lnum_t f_id = blockIdx.x * blockDim.x + threadIdx.x;

  if(f_id >= size){
    return;
  }
  cs_lnum_t c_id1, c_id2;
  cs_real_t pond, ktpond, pfaci, pfacj, rfac;

  c_id1 = i_face_cells[f_id][0];
  c_id2 = i_face_cells[f_id][1];

  pond = weight[f_id];
  ktpond = (c_weight == NULL) ?
        pond :                    // no cell weighting
        pond * c_weight[c_id1] // cell weighting active
          / (      pond * c_weight[c_id1]
            + (1.0-pond)* c_weight[c_id2]);


  for (cs_lnum_t i = 0; i < 3; i++) {
    pfaci = (1.0-ktpond) * (pvar[c_id2][i] - pvar[c_id1][i]);
    pfacj = - ktpond * (pvar[c_id2][i] - pvar[c_id1][i]);

    /* Reconstruction part */
    rfac = 0.5 * (  dofij[f_id][0]*(  r_grad[c_id1][i][0]
                                              + r_grad[c_id2][i][0])
                            + dofij[f_id][1]*(  r_grad[c_id1][i][1]
                                              + r_grad[c_id2][i][1])
                            + dofij[f_id][2]*(  r_grad[c_id1][i][2]
                                              + r_grad[c_id2][i][2]));

    for (cs_lnum_t j = 0; j < 3; j++) {
      atomicAdd(&grad[c_id1][i][j],(pfaci + rfac) * i_f_face_normal[f_id][j]);
      atomicAdd(&grad[c_id2][i][j], - ((pfacj + rfac) * i_f_face_normal[f_id][j]));

    }
  }

}


__global__ static void
_compute_reconstruct_v_b_face(cs_lnum_t            size,
                              const bool                *coupled_faces,
                              cs_lnum_t                 cpl_stride,
                              const cs_real_33_t  *restrict coefbv,
                              const cs_real_3_t   *restrict coefav,
                              const cs_real_3_t   *restrict pvar,
                              int                           inc,
                              const cs_real_3_t *restrict diipb,
                              const cs_real_33_t        *restrict r_grad,
                              cs_real_33_t        *restrict grad,
                              const cs_real_3_t *restrict b_f_face_normal,
                              const cs_lnum_t *restrict b_face_cells)
{
  cs_lnum_t f_id = blockIdx.x * blockDim.x + threadIdx.x;

  if(f_id >= size){
    return;
  }
  cs_lnum_t c_id;
  cs_real_t pfac, rfac, vecfac;

  // if (coupled_faces[f_id * cpl_stride])
  //   return;

  c_id = b_face_cells[f_id];

  for (cs_lnum_t i = 0; i < 3; i++) {

    pfac = inc*coefav[f_id][i];

    for (cs_lnum_t k = 0; k < 3; k++){
      pfac += coefbv[f_id][i][k] * pvar[c_id][k];
    }

    pfac -= pvar[c_id][i];

  //   /* Reconstruction part */
    rfac = 0.;
    for (cs_lnum_t k = 0; k < 3; k++) {
      vecfac =   r_grad[c_id][k][0] * diipb[f_id][0]
                          + r_grad[c_id][k][1] * diipb[f_id][1]
                          + r_grad[c_id][k][2] * diipb[f_id][2];
      rfac += coefbv[f_id][i][k] * vecfac;
    }

    for (cs_lnum_t j = 0; j < 3; j++)
    atomicAdd(&grad[c_id][i][j], (pfac + rfac) * b_f_face_normal[f_id][j]);

  }
}



__global__ static void
_compute_reconstruct_correction(cs_lnum_t            size,
                               cs_lnum_t            has_dc,
                               const int *restrict c_disable_flag,
                               const cs_real_t *restrict cell_f_vol,
                               cs_real_33_t        *restrict grad,
                               const cs_real_33_t *restrict corr_grad_lin,
                               bool                         test_bool
                              )
{
  cs_lnum_t c_id = blockIdx.x * blockDim.x + threadIdx.x;

  if(c_id >= size){
    return;
  }
  cs_real_t dvol;
  /* Is the cell disabled (for solid or porous)? Not the case if coupled */
  if (has_dc * c_disable_flag[has_dc * c_id] == 0)
    dvol = 1. / cell_f_vol[c_id];
  else
    dvol = 0.;


  for (cs_lnum_t i = 0; i < 3; i++) {
    for (cs_lnum_t j = 0; j < 3; j++)
      grad[c_id][i][j] *= dvol;
  }


  if (test_bool) {
    cs_real_t gradpa[3];
    // printf("dvol = %.17lg\n", dvol);
    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        gradpa[j] = grad[c_id][i][j];
        grad[c_id][i][j] = 0.;
      }

      for (cs_lnum_t j = 0; j < 3; j++){
        for (cs_lnum_t k = 0; k < 3; k++){
          grad[c_id][i][j] += corr_grad_lin[c_id][j][k] * gradpa[k];
        }
      }
    }
  }

}

/*----------------------------------------------------------------------------
 * Reconstruct the gradient of a vector using a given gradient of
 * this vector (typically lsq).
 *
 * parameters:
 *   m              <-- pointer to associated mesh structure
 *   fvq            <-- pointer to associated finite volume quantities
 *   cpl            <-- structure associated with internal coupling, or NULL
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   coefav         <-- B.C. coefficients for boundary face normals
 *   coefbv         <-- B.C. coefficients for boundary face normals
 *   pvar           <-- variable
 *   c_weight       <-- weighted gradient coefficient variable
 *   r_grad         --> gradient used for reconstruction
 *   grad           --> gradient of pvar (du_i/dx_j : grad[][i][j])
 *----------------------------------------------------------------------------*/
extern "C" void
cs_reconstruct_vector_gradient_cuda(const cs_mesh_t              *m,
                                    const cs_mesh_adjacencies_t  *madj,
                                    const cs_mesh_quantities_t   *fvq,
                                    const cs_internal_coupling_t *cpl,
                                    cs_halo_type_t                halo_type,
                                    int                           inc,
                                    const cs_real_3_t   *restrict coefav,
                                    const cs_real_33_t  *restrict coefbv,
                                    const cs_real_3_t   *restrict pvar,
                                    const cs_real_t     *restrict c_weight,
                                    const cs_real_33_t  *restrict r_grad,
                                    cs_real_33_t        *restrict grad,
                                    const bool                   *coupled_faces,
                                    cs_lnum_t                     cpl_stride,
                                    bool                          test_bool,
                                    bool                          PERF
                                    )
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces   = m->n_b_faces;
  const cs_lnum_t n_i_faces   = m->n_i_faces;

  int device_id;
  cudaGetDevice(&device_id);


  cudaStream_t stream;
  cudaStreamCreate(&stream);
  
  cudaEvent_t start, mem_h2d, init, i_faces, b_faces_1, b_faces_2, b_faces_3, stop;
  float msec = 0.0f, msec_tot;
  CS_CUDA_CHECK(cudaEventCreate(&start));
  CS_CUDA_CHECK(cudaEventCreate(&mem_h2d));
  CS_CUDA_CHECK(cudaEventCreate(&init));
  CS_CUDA_CHECK(cudaEventCreate(&i_faces));
  CS_CUDA_CHECK(cudaEventCreate(&b_faces_1));
  CS_CUDA_CHECK(cudaEventCreate(&b_faces_2));
  CS_CUDA_CHECK(cudaEventCreate(&b_faces_3));
  CS_CUDA_CHECK(cudaEventCreate(&stop));


  // Record the start event
  CS_CUDA_CHECK(cudaEventRecord(start, stream));

  cs_real_33_t *grad_d;
  CS_CUDA_CHECK(cudaMalloc(&grad_d, n_cells_ext * sizeof(cs_real_33_t)));

  void *_pvar_d = NULL, *_coefa_d = NULL, *_coefb_d = NULL,
  *_cell_cells_idx_d = NULL, *_r_grad_d = NULL;
  const cs_real_3_t *pvar_d = NULL, *coefa_d = NULL;
  const cs_real_33_t *coefb_d = NULL, *r_grad_d = NULL;
  const cs_lnum_t *cell_cells_idx_d = NULL;
  bool  *coupled_faces_d;
  CS_CUDA_CHECK(cudaMalloc(&coupled_faces_d, sizeof(bool) * 2));
  cs_cuda_copy_h2d(coupled_faces_d, coupled_faces, sizeof(bool) * 2);


  unsigned int blocksize = 256;
  unsigned int gridsize_b
    =  (unsigned int)ceil((double)m->n_b_cells / blocksize);
  unsigned int gridsize_if
    = (unsigned int)ceil((double)m->n_i_faces / blocksize);
  unsigned int gridsize_bf
    = (unsigned int)ceil((double)m->n_b_faces / blocksize);
  unsigned int gridsize
      = (unsigned int)ceil((double)m->n_cells / blocksize);
  unsigned int gridsize_init
      = (unsigned int)ceil((double)m->n_cells*3*3 / blocksize);
  unsigned int gridsize_ext
    = (unsigned int)ceil((double)n_cells_ext / blocksize);

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)cs_get_device_ptr_const_pf(m->i_face_cells);
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(m->b_face_cells);
  const cs_lnum_t *restrict cell_b_faces_idx
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(madj->cell_b_faces_idx);
  const cs_lnum_t *restrict cell_cells_lst;
    // = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(m->cell_cells_lst);
  const int n_i_groups
      = m->i_face_numbering->n_groups;
  const int n_i_threads
      = m->i_face_numbering->n_threads;
  cs_lnum_t *restrict i_group_index;
  // printf("m->i_face_numbering->group_index = ", m->i_face_numbering->group_index);
  CS_CUDA_CHECK(cudaMalloc(&i_group_index, sizeof(int)*n_i_groups * n_i_threads * 2));
  cs_cuda_copy_h2d(i_group_index, (void *)m->i_face_numbering->group_index, sizeof(int)*n_i_groups * n_i_threads * 2);
  // = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(m->i_face_numbering->group_index);
  
  const int n_b_groups
      = m->b_face_numbering->n_groups;
  const int n_b_threads
      = m->b_face_numbering->n_threads;

  cs_lnum_t *restrict b_group_index;
  CS_CUDA_CHECK(cudaMalloc(&b_group_index, sizeof(int)*n_i_groups * n_i_threads * 2));
  cs_cuda_copy_h2d(b_group_index, (void *)m->b_face_numbering->group_index, sizeof(int)*n_b_groups * n_b_threads * 2);
  // printf("Avant allocation\n");
  const cs_lnum_t *restrict cell_cells_idx
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(madj->cell_cells_idx);
  const cs_lnum_t *restrict cell_cells
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(madj->cell_cells);
  
  
  
  // if (madj->cell_i_faces == NULL) {
  cs_mesh_adjacencies_update_cell_i_faces();
  // }
  assert(madj->cell_i_faces);
  const cs_lnum_t n_cells_i_face = (madj->cell_cells_idx[n_cells]);
  cs_lnum_t *restrict cell_i_faces;
    // = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(madj->cell_i_faces);
  CS_CUDA_CHECK(cudaMalloc(&cell_i_faces, sizeof(cs_lnum_t)*n_cells_i_face));
  cs_cuda_copy_h2d(cell_i_faces, madj->cell_i_faces, sizeof(cs_lnum_t)*n_cells_i_face);
  assert(cell_i_faces);

  short int *restrict cell_i_faces_sgn;
  CS_CUDA_CHECK(cudaMalloc(&cell_i_faces_sgn, sizeof(short int)*n_cells_i_face));
  cs_cuda_copy_h2d(cell_i_faces_sgn, madj->cell_i_faces_sgn, sizeof(short int)*n_cells_i_face);

  const cs_lnum_t *restrict b_cells
    = (cs_lnum_t *restrict)cs_get_device_ptr_const_pf(m->b_cells);
  const cs_lnum_t *restrict cell_b_faces
    = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(madj->cell_b_faces);

  assert(m->b_cells);
  assert(madj->cell_b_faces);
  assert(madj->cell_b_faces_idx);
  assert(b_cells);
  assert(cell_b_faces);
  assert(cell_b_faces_idx);

  const cs_real_3_t *restrict cell_cen;
    // = (const cs_real_3_t *restrict)cs_get_device_ptr_const_pf(fvq->cell_cen);
  const cs_lnum_t *restrict cell_vol;
    // = (const cs_lnum_t *restrict)cs_get_device_ptr_const_pf(fvq->cell_vol);
  cs_real_t *restrict cell_f_vol;
  CS_CUDA_CHECK(cudaMalloc(&cell_f_vol, n_cells * sizeof(cs_real_t)));
  cs_cuda_copy_h2d(cell_f_vol, (void *)fvq->cell_f_vol, sizeof(cs_real_t)*n_cells);
    // = (const cs_real_t *restrict)cs_get_device_ptr_const_pf(fvq->cell_f_vol);
  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2)
    cell_f_vol = fvq->cell_vol;
  const cs_real_3_t *restrict cell_f_cen
    = (const cs_real_3_t *restrict)cs_get_device_ptr_const_pf(fvq->cell_f_cen);
  const cs_real_t *restrict weight
    = (const cs_real_t *restrict)cs_get_device_ptr_const_pf(fvq->weight);
  const cs_real_t *restrict b_dist;
    // = (const cs_real_t *restrict)cs_get_device_ptr_const_pf(fvq->b_dist);
  const cs_real_3_t *restrict b_face_normal;
    // = (const cs_real_3_t *restrict)cs_get_device_ptr_const_pf(fvq->b_face_normal);
  cs_real_3_t *restrict i_f_face_normal;
  // printf("fvq->i_f_face_normal = ", fvq->i_f_face_normal);
  CS_CUDA_CHECK(cudaMalloc(&i_f_face_normal, sizeof(cs_real_3_t)*n_i_faces));
  cs_cuda_copy_h2d(i_f_face_normal, (void *)fvq->i_f_face_normal, sizeof(cs_real_3_t)*n_i_faces);
  // = (const cs_real_3_t *restrict)cs_get_device_ptr_const_pf(fvq->i_f_face_normal);

  const cs_real_3_t *restrict b_f_face_normal
  = (const cs_real_3_t *restrict)cs_get_device_ptr_const_pf(fvq->b_f_face_normal);
  cs_real_3_t *restrict dofij;
  // printf("fvq->dofij = ", fvq->dofij);
  CS_CUDA_CHECK(cudaMalloc(&dofij, sizeof(cs_real_3_t)*n_i_faces));
  cs_cuda_copy_h2d(dofij, (void *)fvq->dofij, sizeof(cs_real_3_t)*n_i_faces);
  // = (const cs_real_3_t *restrict)cs_get_device_ptr_const_pf(fvq->dofij);
  const cs_real_3_t *restrict diipb
    = (const cs_real_3_t *restrict)cs_get_device_ptr_const_pf(fvq->diipb);
  cs_real_33_t *restrict corr_grad_lin;
  CS_CUDA_CHECK(cudaMalloc(&corr_grad_lin, n_cells * sizeof(cs_real_33_t)));
  cs_cuda_copy_h2d(corr_grad_lin, (void *)fvq->corr_grad_lin, sizeof(cs_real_33_t)*n_cells);
    // = (const cs_real_33_t *restrict)cs_get_device_ptr_const_pf(fvq->corr_grad_lin);
  const cs_lnum_t has_dc
      = fvq->has_disable_flag;
  int *restrict c_disable_flag;
  CS_CUDA_CHECK(cudaMalloc(&c_disable_flag, n_cells * sizeof(int)));
  cs_cuda_copy_h2d(c_disable_flag, (void *)fvq->c_disable_flag, sizeof(int)*n_cells);
    // = (const int *restrict)cs_get_device_ptr_const_pf(fvq->c_disable_flag);


  _sync_or_copy_real_h2d(pvar, n_cells_ext, device_id, stream,
    &pvar_d, &_pvar_d);

  _sync_or_copy_real_h2d(r_grad, n_cells_ext, device_id, stream,
    &r_grad_d, &_r_grad_d);

  _sync_or_copy_real_h2d(coefav, n_b_faces, device_id, stream,
        &coefa_d, &_coefa_d);
  _sync_or_copy_real_h2d(coefbv, n_b_faces, device_id, stream,
        &coefb_d, &_coefb_d);
      

  // ----------------------------Begin of Kernels part 1-------------------------------------------
  
  CS_CUDA_CHECK(cudaEventRecord(mem_h2d, stream));

  /* Initialization */

  cudaMemset(grad_d, 0, n_cells * sizeof(cs_real_33_t));

  CS_CUDA_CHECK(cudaEventRecord(init, stream));
  
  
  /* Interior faces contribution */
  // _compute_reconstruct_v_i_face<<<gridsize_if, blocksize, 0, stream>>>
  //                               (n_i_faces,
  //                               i_group_index,
  //                               i_face_cells,
  //                               pvar_d,
  //                               weight,
  //                               c_weight,
  //                               r_grad_d,
  //                               grad_d,
  //                               dofij,
  //                               i_f_face_normal);

  // _compute_reconstruct_v_i_face_v2<<<gridsize_if * 3, blocksize, 0, stream>>>
  //                               (n_i_faces * 3,
  //                               i_group_index,
  //                               i_face_cells,
  //                               pvar_d,
  //                               weight,
  //                               c_weight,
  //                               r_grad_d,
  //                               grad_d,
  //                               dofij,
  //                               i_f_face_normal);

  _compute_reconstruct_v_i_face_v2cf<<<gridsize_if * 3, blocksize, 0, stream>>>
                                (n_i_faces * 3,
                                i_group_index,
                                i_face_cells,
                                pvar_d,
                                weight,
                                c_weight,
                                r_grad_d,
                                grad_d,
                                dofij,
                                i_f_face_normal);
  
  // printf("Avant les assert dans gradient_cuda.cu\n");
  // assert(cell_cells_idx);
  // assert(cell_cells);
  // assert(weight);
  // assert(cell_i_faces);
  // assert(cell_i_faces_sgn);
  // printf("n_i_faces = %d\n", n_i_faces);
  // printf("n_cells = %d\n", n_cells);
  // for(int i = 0; i< n_i_faces; i++){
  //   printf("i = %d && weight = %f \n", i, fvq->weight[i]);
  //   printf("i = %d && c_id2 = %d \n", i, madj->cell_cells[i]);
  //   printf("i = %d && s_id = %d \n", i, madj->cell_cells_idx[i]);
  //   printf("i = %d && f_id = %d \n", i, madj->cell_i_faces_sgn[i]);
  // }
  // printf("Après les assert dans gradient_cuda.cu\n");
  // _compute_reconstruct_v_i_face_gather<<<gridsize, blocksize, 0, stream>>>
  //                                     ( n_cells,
  //                                       pvar_d,
  //                                       weight,
  //                                       c_weight,
  //                                       r_grad_d,
  //                                       grad_d,
  //                                       dofij,
  //                                       i_f_face_normal,
  //                                       cell_cells_idx,
  //                                       cell_cells,
  //                                       cell_i_faces,
  //                                       cell_i_faces_sgn);


  // _compute_reconstruct_v_i_face_gather_v2<<<gridsize * 3 * 3, blocksize, 0, stream>>>
  //                                     ( n_cells * 3 * 3,
  //                                       i_face_cells,
  //                                       pvar_d,
  //                                       weight,
  //                                       c_weight,
  //                                       r_grad_d,
  //                                       grad_d,
  //                                       dofij,
  //                                       i_f_face_normal,
  //                                       cell_cells_idx,
  //                                       cell_cells,
  //                                       cell_i_faces,
  //                                       cell_i_faces_sgn,
  //                                       n_i_faces);

  CS_CUDA_CHECK(cudaEventRecord(i_faces, stream));
  
  // ----------------------------End of Kernels part 1-------------------------------------------

  // if (grad_d != NULL) {
  //   size_t size = n_cells_ext * sizeof(cs_real_t) * 3 * 3;
  //   cs_cuda_copy_d2h(grad, grad_d, size);
  // }
  // else
  //   cs_sync_d2h(grad);
    
  // size_t size = n_cells_ext * sizeof(cs_real_t) * 3 * 3;
  // cs_cuda_copy_d2h(r_grad, r_grad_d, size);
    

  /* Contribution from coupled faces */
  // if (cpl != NULL) {
  //   cs_internal_coupling_initialize_vector_gradient(cpl, c_weight, pvar, grad);
  //   cs_internal_coupling_reconstruct_vector_gradient(cpl, r_grad, grad);
  // }
  
  // cs_cuda_copy_h2d(grad_d, grad, n_cells_ext * sizeof(cs_real_33_t));

  // _sync_or_copy_real_h2d(r_grad, n_cells_ext, device_id, stream,
  //   &r_grad_d, &_r_grad_d);

  CS_CUDA_CHECK(cudaEventRecord(b_faces_1, stream));
  
  // ----------------------------Begin of Kernels part 2-------------------------------------------
  // _compute_reconstruct_v_b_face<<<gridsize_bf, blocksize, 0, stream>>>
  //                             ( n_b_faces,
  //                               coupled_faces_d,
  //                               cpl_stride,
  //                               coefb_d,
  //                               coefa_d,
  //                               pvar_d,
  //                               inc,
  //                               diipb,
  //                               r_grad_d,
  //                               grad_d,
  //                               b_f_face_normal,
  //                               b_face_cells);


  // _compute_reconstruct_v_b_face_v2<<<gridsize_bf * 3, blocksize, 0, stream>>>
  //                             ( n_b_faces * 3,
  //                               coupled_faces_d,
  //                               cpl_stride,
  //                               coefb_d,
  //                               coefa_d,
  //                               pvar_d,
  //                               inc,
  //                               diipb,
  //                               r_grad_d,
  //                               grad_d,
  //                               b_f_face_normal,
  //                               b_face_cells);
  

  _compute_reconstruct_v_b_face_v2cf<<<gridsize_bf * 3, blocksize, 0, stream>>>
                              ( n_b_faces * 3,
                                coupled_faces_d,
                                cpl_stride,
                                coefb_d,
                                coefa_d,
                                pvar_d,
                                inc,
                                diipb,
                                r_grad_d,
                                grad_d,
                                b_f_face_normal,
                                b_face_cells);

  // _compute_reconstruct_v_b_face_gather<<<gridsize_b, blocksize, 0, stream>>>
  //                             ( m->n_b_cells,
  //                               coupled_faces_d,
  //                               cpl_stride,
  //                               coefb_d,
  //                               coefa_d,
  //                               pvar_d,
  //                               inc,
  //                               diipb,
  //                               r_grad_d,
  //                               grad_d,
  //                               b_f_face_normal,
  //                               b_cells,
  //                               cell_b_faces,
  //                               cell_b_faces_idx);


  // _compute_reconstruct_v_b_face_gather_v2<<<gridsize_b * 3, blocksize, 0, stream>>>
  //                             ( m->n_b_cells * 3,
  //                               coupled_faces_d,
  //                               cpl_stride,
  //                               coefb_d,
  //                               coefa_d,
  //                               pvar_d,
  //                               inc,
  //                               diipb,
  //                               r_grad_d,
  //                               grad_d,
  //                               b_f_face_normal,
  //                               b_cells,
  //                               cell_b_faces,
  //                               cell_b_faces_idx);

  CS_CUDA_CHECK(cudaEventRecord(b_faces_2, stream));
  
  // _compute_reconstruct_correction<<<gridsize, blocksize, 0, stream>>>
  //                             ( n_cells,
  //                               has_dc,
  //                               c_disable_flag,
  //                               cell_f_vol,
  //                               grad_d,
  //                               corr_grad_lin,
  //                               test_bool
  //                             );

  _compute_reconstruct_correction_v2<<<gridsize * 3, blocksize, 0, stream>>>
                              ( n_cells * 3,
                                has_dc,
                                c_disable_flag,
                                cell_f_vol,
                                grad_d,
                                corr_grad_lin,
                                test_bool
                              );
  CS_CUDA_CHECK(cudaEventRecord(b_faces_3, stream));

  // ----------------------------End of Kernels part 2-------------------------------------------

  /* Sync to host */
  if (grad_d != NULL) {
    size_t size = n_cells_ext * sizeof(cs_real_t) * 3 * 3;
    cs_cuda_copy_d2h(grad, grad_d, size);
  }
  else
    cs_sync_d2h(grad);
    

  CS_CUDA_CHECK(cudaEventRecord(stop, stream));
  CS_CUDA_CHECK(cudaEventSynchronize(stop));

  cudaStreamSynchronize(stream);
  cudaStreamDestroy(stream);

  if(PERF){
    printf("recconstruct Kernels times:\t");

    msec = 0.0f;
    CS_CUDA_CHECK(cudaEventElapsedTime(&msec, mem_h2d, init));
    printf("Kernels execution time in us: \t");
    printf("Init = %f\t", msec*1000.f);

    msec = 0.0f;
    CS_CUDA_CHECK(cudaEventElapsedTime(&msec, init, i_faces));
    printf("I_faces = %f\t", msec*1000.f);

    // msec = 0.0f;
    // CS_CUDA_CHECK(cudaEventElapsedTime(&msec, i_faces, b_faces_1));
    // printf("CPU part = %f\t", msec*1000.f);

    msec = 0.0f;
    CS_CUDA_CHECK(cudaEventElapsedTime(&msec, b_faces_1, b_faces_2));
    printf("B_faces = %f\t", msec*1000.f);

    msec = 0.0f;
    CS_CUDA_CHECK(cudaEventElapsedTime(&msec, b_faces_2, b_faces_3));
    printf("Correction = %f\t", msec*1000.f);

    printf("\n");

    msec_tot = 0.0f;
    msec = 0.0f;
    CS_CUDA_CHECK(cudaEventElapsedTime(&msec, mem_h2d, i_faces));
    printf("reconstruct Total kernel part 1= %f\t", msec*1000.f);
    msec_tot = msec;

    msec = 0.0f;
    CS_CUDA_CHECK(cudaEventElapsedTime(&msec, b_faces_1, b_faces_3));
    printf("Total kernel part 2= %f\t", msec*1000.f);
    msec_tot += msec;

    printf("Total kernel 1 and 2= %f\t", msec_tot*1000.f);

    msec = 0.0f;
    CS_CUDA_CHECK(cudaEventElapsedTime(&msec, start, stop));
    printf("Total = %f\t", msec*1000.f);

    printf("\n");
  }

  if (_pvar_d != NULL)
    CS_CUDA_CHECK(cudaFree(_pvar_d));
  if (_coefa_d != NULL)
    CS_CUDA_CHECK(cudaFree(_coefa_d));
  if (_coefb_d != NULL)
    CS_CUDA_CHECK(cudaFree(_coefb_d));
  if (_r_grad_d != NULL)
    CS_CUDA_CHECK(cudaFree(_r_grad_d));

  CS_CUDA_CHECK(cudaFree(cell_i_faces));
  CS_CUDA_CHECK(cudaFree(cell_i_faces_sgn));

  CS_CUDA_CHECK(cudaFree(coupled_faces_d));
  CS_CUDA_CHECK(cudaFree(i_group_index));
  CS_CUDA_CHECK(cudaFree(b_group_index));
  CS_CUDA_CHECK(cudaFree(cell_f_vol));
  CS_CUDA_CHECK(cudaFree(i_f_face_normal));
  CS_CUDA_CHECK(cudaFree(dofij));
  CS_CUDA_CHECK(cudaFree(corr_grad_lin));
  CS_CUDA_CHECK(cudaFree(c_disable_flag));
  CS_CUDA_CHECK(cudaFree(grad_d));
}
