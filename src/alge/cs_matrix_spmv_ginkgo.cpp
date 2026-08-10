/*============================================================================
 * Sparse Matrix SpMV operations using Ginkgo.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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
 * Standard headers
 *----------------------------------------------------------------------------*/

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <new>

/*----------------------------------------------------------------------------
 * Ginkgo headers
 *----------------------------------------------------------------------------*/

#include <ginkgo/ginkgo.hpp>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"
#include "base/cs_base_accel.h"
#include "base/cs_halo.h"
#include "base/cs_mem.h"
#include "alge/cs_matrix.h"
#include "alge/cs_matrix_priv.h"
#include "alge/cs_matrix_spmv.h"
#include "alge/cs_matrix_spmv_ginkgo.h"

/*----------------------------------------------------------------------------*/
/*! \file cs_matrix_spmv_ginkgo.cpp
 *
 * \brief Sparse Matrix SpMV operations using Ginkgo library.
 */
/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

struct cs_matrix_ginkgo_map_t {

  int  use_device;                          /* 0 for host, 1 for device */

  std::shared_ptr<gko::Executor> executor;  /* Ginkgo executor */

  /* Ginkgo LinOp CSR matrix (view) */
  std::unique_ptr<const gko::matrix::Csr<cs_real_t, cs_lnum_t>> mat_a;

  /* Ginkgo LinOp diagonal (for diagonal exclusion) */
  std::unique_ptr<const gko::matrix::Diagonal<cs_real_t>> mat_d;

  /* Mapping to multiplier and SpMV output */
  std::unique_ptr<gko::matrix::Dense<cs_real_t>> vec_x;
  std::unique_ptr<gko::matrix::Dense<cs_real_t>> vec_y;

  /* LinOp Scaling factors for y = y - D * x */
  std::unique_ptr<const gko::matrix::Dense<cs_real_t>> alpha_neg_one;
  std::unique_ptr<const gko::matrix::Dense<cs_real_t>> beta_one;

  /* Raw pointers to check for minimal updates */
  const cs_real_t  *ptr_val;                /* Raw pointer to matrix values */
  const cs_real_t  *ptr_ad;                 /* Raw pointer to diag values */
  const void       *ptr_x;                  /* Raw pointer to x values */
  const void       *ptr_y;                  /* Raw pointer to y values */

  bool  mapped;                             /* Mapping initialized */

};

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return or create the Ginkgo executor matching device/host selection.
 *
 * parameters:
 *   use_device <-- 0 for host, 1 for device (GPU)
 *
 * returns:
 *   shared pointer to Ginkgo executor
 *----------------------------------------------------------------------------*/

static std::shared_ptr<gko::Executor>
_get_executor(int  use_device)
{
  if (use_device == 1) {

#if defined(HAVE_CUDA)
    int dev_id = cs_get_device_id();
    if (dev_id >= 0) {
      cudaStream_t stream = cs_cuda_get_stream(0);
      auto allocator = std::make_shared<gko::CudaAllocator>();

      return gko::CudaExecutor::create(dev_id,
                                       gko::OmpExecutor::create(),
                                       allocator,
                                       stream);
    }
#elif defined(HAVE_HIP)
    int dev_id = cs_get_device_id();
    if (dev_id >= 0) {
      hipStream_t stream = cs_hip_get_stream(0);
      auto allocator = std::make_shared<gko::HipAllocator>(stream);

      return gko::HipExecutor::create(dev_id,
                                      gko::OmpExecutor::create(),
                                      false,
                                      allocator,
                                      stream);
    }
#endif
    bft_error(__FILE__, __LINE__, 0,
              _("Ginkgo device execution requested, but no\n"
                "device acceleration is available in this build."));
    return nullptr;

  } // use_device

  else {
#if defined(_OPENMP)
    return gko::OmpExecutor::create();
#else
    return gko::ReferenceExecutor::create();
#endif
  }
}

/*----------------------------------------------------------------------------
 * Unset matrix Ginkgo mapping.
 *
 * parameters:
 *   matrix <-> pointer to matrix structure
 *----------------------------------------------------------------------------*/

static void
_unset_ginkgo_map(cs_matrix_t  *matrix)
{
  auto *csm = static_cast<cs_matrix_ginkgo_map_t *>(matrix->ext_lib_map);

  if (csm == nullptr)
    return;

  if (csm->mapped) {
    csm->mat_a.reset();
    if (csm->mat_d != nullptr) {
      csm->mat_d.reset();
      csm->alpha_neg_one.reset();
      csm->beta_one.reset();
    }
    csm->executor.reset();
    if (csm->vec_x != nullptr)
      csm->vec_x.reset();
    if (csm->vec_y != nullptr)
      csm->vec_y.reset();
  }

  // Explicit destruction as we use placement new
  csm->~cs_matrix_ginkgo_map_t();

  CS_FREE(matrix->ext_lib_map);
  matrix->destroy_adaptor = nullptr;
}

/*----------------------------------------------------------------------------
 * Set matrix Ginkgo mapping.
 *
 * parameters:
 *   matrix     <-> pointer to matrix structure
 *   use_device <-- 0 for host, 1 for device
 *
 * returns:
 *   pointer to Ginkgo mapping structure
 *----------------------------------------------------------------------------*/

static cs_matrix_ginkgo_map_t *
_set_ginkgo_map(cs_matrix_t  *matrix,
                int           use_device)
{
  auto *csm = static_cast<cs_matrix_ginkgo_map_t *>(matrix->ext_lib_map);

  if (csm != nullptr) {
    _unset_ginkgo_map(matrix);
    csm = nullptr;
  }

  CS_MALLOC(csm, 1, cs_matrix_ginkgo_map_t);
  new(csm) cs_matrix_ginkgo_map_t();
  csm->use_device = use_device;
  csm->executor = _get_executor(use_device);
  csm->mapped = false;
  csm->ptr_val = nullptr;
  csm->ptr_ad = nullptr;
  csm->ptr_x = nullptr;
  csm->ptr_y = nullptr;
  matrix->ext_lib_map = static_cast<void *>(csm);
  matrix->destroy_adaptor = _unset_ginkgo_map;

  if (matrix->type != CS_MATRIX_CSR)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: Ginkgo SpMV only supports CSR matrices."), __func__);

  const auto *ms
    = static_cast<const cs_matrix_struct_csr_t *>(matrix->structure);
  const auto *mc
    = static_cast<const cs_matrix_coeff_t *>(matrix->coeffs);

  cs_lnum_t n_rows = ms->n_rows;
  cs_lnum_t n_cols_ext = ms->n_cols_ext;
  cs_lnum_t nnz = ms->row_index[n_rows];

  const cs_lnum_t *row_index = ms->row_index;
  const cs_lnum_t *col_id = ms->col_id;
  const cs_real_t *val = mc->val;

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  if (use_device == 1) {
    row_index = static_cast<const cs_lnum_t *>
      (cs_get_device_ptr_const(const_cast<cs_lnum_t *>(ms->row_index)));
    col_id = static_cast<const cs_lnum_t *>
      (cs_get_device_ptr_const(const_cast<cs_lnum_t *>(ms->col_id)));
    val = static_cast<const cs_real_t *>
      (cs_get_device_ptr_const(const_cast<cs_real_t *>(mc->val)));
  }
#endif

  csm->mat_a = gko::matrix::Csr<cs_real_t, cs_lnum_t>::create_const
    (csm->executor,
     gko::dim<2>{static_cast<gko::size_type>(n_rows),
                 static_cast<gko::size_type>(n_cols_ext)},
     gko::make_const_array_view(csm->executor, nnz, val),
     gko::make_const_array_view(csm->executor, nnz, col_id),
     gko::make_const_array_view(csm->executor, n_rows + 1, row_index));

  csm->ptr_val = mc->val;
  csm->mapped = true;

  return csm;
}

/*----------------------------------------------------------------------------
 * Update matrix Ginkgo mapping.
 *
 * parameters:
 *   csm       <-> Ginkgo matrix mapping
 *   matrix    <-> pointer to matrix structure
 *   x         <-- pointer to input vector (on execution location)
 *   y         <-- pointer to output vector (on execution location)
 *----------------------------------------------------------------------------*/

static void
_update_ginkgo_map(cs_matrix_ginkgo_map_t    *csm,
                   const cs_matrix_t         *matrix,
                   cs_real_t                 *x,
                   cs_real_t                 *y)
{
  assert(csm != nullptr);

  if (csm->ptr_x != static_cast<void *>(x)) {
    auto n_cols = static_cast<gko::size_type>(matrix->n_cols_ext);
    csm->vec_x = gko::matrix::Dense<cs_real_t>::create
      (csm->executor,
       gko::dim<2>{n_cols, 1},
       gko::make_array_view(csm->executor, n_cols, x),
       1);
    csm->ptr_x = static_cast<void *>(x);
  }

  if (csm->ptr_y != static_cast<void *>(y)) {
    auto n_rows = static_cast<gko::size_type>(matrix->n_rows);
    csm->vec_y = gko::matrix::Dense<cs_real_t>::create
      (csm->executor,
       gko::dim<2>{n_rows, 1},
       gko::make_array_view(csm->executor, n_rows, y),
       1);
    csm->ptr_y = static_cast<void *>(y);
  }
}

/*----------------------------------------------------------------------------
 * Initiate ghost cell synchronization for multiplying vector.
 *----------------------------------------------------------------------------*/

static cs_halo_state_t *
_pre_vector_multiply_sync_x_start(const cs_matrix_t  *matrix,
                                  cs_real_t          *restrict x)
{
  cs_halo_state_t *hs = nullptr;

  if (matrix->halo != nullptr) {
    hs = cs_halo_state_get_default();

    cs_halo_sync_pack(matrix->halo,
                      CS_HALO_STANDARD,
                      CS_REAL_TYPE,
                      matrix->db_size,
                      x,
                      nullptr,
                      hs);

    cs_halo_sync_start(matrix->halo, x, hs);
  }

  return hs;
}

#if defined(HAVE_ACCEL)

/*----------------------------------------------------------------------------
 * Initiate ghost cell synchronization on device for multiplying vector.
 *----------------------------------------------------------------------------*/

static cs_halo_state_t *
_pre_vector_multiply_sync_x_start_d(const cs_matrix_t  *matrix,
                                    cs_real_t          *restrict x)
{
  cs_halo_state_t *hs = nullptr;

  if (matrix->halo != nullptr) {
    hs = cs_halo_state_get_default();

    cs_halo_sync_pack_d(matrix->halo,
                        CS_HALO_STANDARD,
                        CS_REAL_TYPE,
                        matrix->db_size,
                        x,
                        nullptr,
                        hs);

    cs_halo_sync_start(matrix->halo, x, hs);
  }

  return hs;
}

#endif /* defined(HAVE_ACCEL) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Main part of matrix-vector product y = A.x with CSR matrix,
 *       scalar Ginkgo version.
 *
 * \param[in]   matrix        pointer to matrix structure
 * \param[in]   exclude_diag  exclude diagonal if true
 * \param[in]   use_device    run on accelerated device ?
 * \param[in]   hs            halo state if ghost cells are present
 * \param[in]   x             multiplying vector values
 * \param[out]  y             resulting vector
 */
/*----------------------------------------------------------------------------*/

static void
_matrix_spmv_ginkgo_csr_body(cs_matrix_t      *matrix,
                             bool              exclude_diag,
                             bool              use_device,
                             cs_halo_state_t  *hs,
                             cs_real_t         x[],
                             cs_real_t         y[])
{
  assert(matrix != nullptr);
  assert(matrix->type == CS_MATRIX_CSR);

  /* Map matrix if not yet done */

  auto *csm = static_cast<cs_matrix_ginkgo_map_t *>(matrix->ext_lib_map);
  const auto *mc = static_cast<const cs_matrix_coeff_t *>(matrix->coeffs);

  if (csm == nullptr || csm->use_device != use_device || csm->ptr_val != mc->val)
    csm = _set_ginkgo_map(matrix, 0);

  /* Finalize ghost cell communication */

  if (hs != nullptr)
    cs_halo_sync_wait(matrix->halo, x, hs);

  _update_ginkgo_map(csm, matrix, x, y);

  csm->mat_a->apply(csm->vec_x.get(), csm->vec_y.get());

  if (exclude_diag) {
    const cs_real_t *ad = cs_matrix_get_diagonal(matrix);
    auto n_rows = static_cast<gko::size_type>(matrix->n_rows);

    if (csm->mat_d == nullptr || csm->ptr_ad != ad) {
      csm->mat_d = gko::matrix::Diagonal<cs_real_t>::create_const
        (csm->executor,
         n_rows,
         gko::make_const_array_view(csm->executor, n_rows, ad));
      csm->ptr_ad = ad;
    }

    if (csm->alpha_neg_one == nullptr) {
      csm->alpha_neg_one = gko::initialize<gko::matrix::Dense<cs_real_t>>
        ({-1.0}, csm->executor);
      csm->beta_one = gko::initialize<gko::matrix::Dense<cs_real_t>>
        ({1.0}, csm->executor);
    }

    auto vec_x_loc = gko::matrix::Dense<cs_real_t>::create
      (csm->executor,
       gko::dim<2>{n_rows, 1},
       gko::make_array_view(csm->executor, n_rows, x),
       1);

    csm->mat_d->apply(csm->alpha_neg_one.get(),
                      vec_x_loc.get(),
                      csm->beta_one.get(),
                      csm->vec_y.get());
  }

  csm->executor->synchronize();
}

/*! \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print Ginkgo library information.
 *
 * \param[in]  log_type  destination log type
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_ginkgo_library_info(cs_log_t  log_type)
{
  auto version = gko::version_info::get();

  cs_log_printf(log_type,
                _("    Ginkgo %lu.%lu.%lud (%s)\n"),
                version.core_version.major,
                version.core_version.minor,
                version.core_version.patch,
                version.core_version.tag);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix-vector product y = A.x with CSR matrix, scalar Ginkgo version.
 *
 * \param[in]   matrix        pointer to matrix structure
 * \param[in]   exclude_diag  exclude diagonal if true
 * \param[in]   sync          synchronize ghost cells if true
 * \param[in]   x             multiplying vector values
 * \param[out]  y             resulting vector
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_ginkgo_csr(cs_matrix_t  *matrix,
                          bool          exclude_diag,
                          bool          sync,
                          cs_real_t     x[],
                          cs_real_t     y[])
{
  assert(matrix != nullptr);
  assert(matrix->type == CS_MATRIX_CSR);

  /* Ghost cell communication */

  cs_halo_state_t *hs
    = (sync) ? _pre_vector_multiply_sync_x_start(matrix, x) : nullptr;

  /* Main function */

  _matrix_spmv_ginkgo_csr_body(matrix, exclude_diag, false, hs, x, y);
}

#if defined(HAVE_ACCEL)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Matrix-vector product y = A.x with CSR matrix, scalar Ginkgo device.
 *
 * \param[in]   matrix        pointer to matrix structure
 * \param[in]   exclude_diag  exclude diagonal if true
 * \param[in]   sync          synchronize ghost cells if true
 * \param[in]   x             multiplying vector values (at execution location)
 * \param[out]  y             resulting vector (at execution location)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_spmv_ginkgo_device_csr(cs_matrix_t  *matrix,
                                 bool          exclude_diag,
                                 bool          sync,
                                 cs_real_t     x[],
                                 cs_real_t     y[])
{
  assert(matrix != nullptr);
  assert(matrix->type == CS_MATRIX_CSR);

  /* Ghost cell communication */

  cs_halo_state_t *hs = nullptr;
#if defined(HAVE_ACCEL)
  if (sync)
    hs = _pre_vector_multiply_sync_x_start_d(matrix, x);
#endif

  /* Main function */

  _matrix_spmv_ginkgo_csr_body(matrix, exclude_diag, true, hs, x, y);
}

#endif /* defined(HAVE_ACCEL) */

/*----------------------------------------------------------------------------*/
