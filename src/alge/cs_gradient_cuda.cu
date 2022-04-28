/*============================================================================
 * Gradient reconstruction, CUDA implementations.
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
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_adjacencies.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_porous_model.h"
#include "cs_prototypes.h"
#include "cs_timer.h"
#include "cs_timer_stats.h"

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

/*----------------------------------------------------------------------------
 * Synchronize of copy a cs_real_t type array from the host to a device.
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

static void
_sync_or_copy_real_h2d(const  cs_real_t   *val_h,
                       cs_lnum_t           n_vals,
                       int                 device_id,
                       cudaStream_t        stream,
                       const cs_real_t   **val_d,
                       void              **buf_d)
{
  const cs_real_t  *_val_d = NULL;
  void             *_buf_d = NULL;

  cs_alloc_mode_t alloc_mode = cs_check_device_ptr(val_h);
  size_t size = n_vals * sizeof(cs_real_t);

  if (alloc_mode == CS_ALLOC_HOST) {
    CS_CUDA_CHECK(cudaMalloc(&_buf_d, size));
    cs_cuda_copy_h2d(_buf_d, val_h, size);
    _val_d = (const cs_real_t *)_buf_d;
  }
  else {
    _val_d = (const cs_real_t *)cs_get_device_ptr((void *)val_h);

    if (alloc_mode == CS_ALLOC_HOST_DEVICE_SHARED)
      cudaMemPrefetchAsync(val_h, size, device_id, stream);
    else
      cs_sync_h2d(val_h);
  }

  *val_d = _val_d;
  *buf_d = _buf_d;
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
 *   hyd_p_flag     <-- flag for hydrostatic pressure
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   fext           <-- exterior force generating pressure
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
                            int                           hyd_p_flag,
                            cs_real_t                     inc,
                            const cs_real_3_t            *f_ext,
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

  /* External forces not handled yet */

  CS_NO_WARN_IF_UNUSED(f_ext);

  /* Reconstruct gradients using least squares for non-orthogonal meshes */
  /*---------------------------------------------------------------------*/

  /* Recompute cocg and at rhsv from interior cells */

  if (hyd_p_flag == 0 && init_cocg) {

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
  else if (hyd_p_flag == 0) {

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

  /* Case with hydrostatic pressure */
  /*--------------------------------*/

  else { /* if hyd_p_flag == 1 */

    bft_error(__FILE__, __LINE__, 0,
              "%s does not support hyd_p_flag == 1 yet", __func__);

  } /* End of test on hydrostatic pressure */

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
