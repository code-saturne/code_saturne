/*============================================================================
 * Various base algorithms.
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
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <chrono>

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>

#if defined(HAVE_CUDA)
#include <cub/cub.cuh>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"

#include "base/cs_base.h"
#include "base/cs_dispatch.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "base/cs_parall.h"
#include "base/cs_reducers.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_algorithm.h"

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*!
  \file cs_algorithm.cpp
        Various base algorithms.
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

#if defined(HAVE_CUDA)

/*--------------------------------------------------------------------------*/
/*!
 * \brief Functor type for selecting values greater than some cutoff
 */
/*--------------------------------------------------------------------------*/

struct cs_algorithm_greater_than
{
  cs_lnum_t cutoff;

  __host__ __device__ __forceinline__
  cs_algorithm_greater_than(cs_lnum_t  c) : cutoff(c) {}

  __host__ __device__ __forceinline__
  bool operator()(const cs_lnum_t  &a) const {
    return (a > cutoff);
  }
};

#endif // defined(HAVE_CUDA)

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

#if defined(HAVE_OPENMP)

/*--------------------------------------------------------------------------*/
/*!
 * \brief Transform a count to an index in place, OpenMP threaded version.
 *
 * For n input elements, the array size should be size n+1, to account
 * for the past-the-end count.
 *
 * \param[in]       n_threads  number of threads
 * \param[in]       n          number of elements
 * \param[in, out]  a <->      count in, index out (size: n+1)
 */
/*--------------------------------------------------------------------------*/

static void
_count_to_index_inplace_omp(int        n_threads,
                            cs_lnum_t  n,
                            cs_lnum_t  a[])
{
  assert(a != nullptr && n > 0);

  constexpr int cl_size_max = 128;  // Max expected cache line size
  constexpr int n_t_noalloc = 32;

  // Space partial counts by stride to try to avoid false sharing.
  constexpr cs_lnum_t stride = (cl_size_max / sizeof(cs_lnum_t))
                               + ((cl_size_max % sizeof(cs_lnum_t) == 0) ? 0 : 1);
  cs_lnum_t partial_sum_[cl_size_max * n_t_noalloc];
  cs_lnum_t *partial_sum = partial_sum_;
  if (n_threads > n_t_noalloc)
    CS_MALLOC(partial_sum, n_threads*stride, cs_lnum_t);

  // Compute partial sums (only reading from a, not writing to it)

  #pragma omp parallel shared(partial_sum) num_threads(n_threads)
  {
    int t_id = omp_get_thread_num();
    cs_lnum_t t_s_id, t_e_id;
    cs_parall_thread_range(n, sizeof(cs_lnum_t), t_id, n_threads,
                           &t_s_id, &t_e_id);

    cs_lnum_t l_count = 0;
    for (cs_lnum_t i = t_s_id; i < t_e_id; i++)
      l_count += a[i];

    partial_sum[t_id*stride] = l_count;
  }

  // Serial inclusive scan of partial sums.

  {
    cs_lnum_t s = partial_sum[0];
    for (int j = 1; j < n_threads; j++) {
      s += partial_sum[j*stride];
      partial_sum[j*stride] = s;
    }
  }

  // Now finalize computation of index.

  #pragma omp parallel shared(partial_sum) num_threads(n_threads)
  {
    int t_id = omp_get_thread_num();
    cs_lnum_t t_s_id, t_e_id;
    cs_parall_thread_range(n, sizeof(cs_lnum_t), t_id, n_threads,
                           &t_s_id, &t_e_id);

    cs_lnum_t s = (t_id == 0) ? 0 : partial_sum[(t_id-1)*stride];

    for (cs_lnum_t i = t_s_id; i < t_e_id; i++) {
      cs_lnum_t c = a[i];
      a[i] = s;
      s += c;
    }

    if (t_id == n_threads-1)  // Last thread handles past-the end value
      a[n] = s;
  }

  if (partial_sum != partial_sum_)
    CS_FREE(partial_sum);
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief Select and compact elements whose values are greater than a
 *        given number, OpenMP threaded version.
 *
 * \param[in]       n_threads    number of threads
 * \param[in]       n            number of elements
 * \param[in]       c            cutoff value to compare to
 * \param[in, out]  a            elements in, selected elements out
 *
 * \return  number of selected elements
 */
/*--------------------------------------------------------------------------*/

static cs_lnum_t
_select_if_gt_omp(int        n_threads,
                  cs_lnum_t  n,
                  cs_lnum_t  c,
                  cs_lnum_t  a[])
{
  cs_lnum_t n_select = 0;

  assert(a != nullptr && n > 0);

  constexpr int cl_size_max = 128;  // Max expected cache line size
  constexpr int n_t_noalloc = 32;

  // Space partial counts by stride to try to avoid false sharing.
  constexpr cs_lnum_t stride = (cl_size_max / sizeof(cs_lnum_t))
                               + ((cl_size_max % sizeof(cs_lnum_t) == 0) ? 0 : 1);
  cs_lnum_t partial_count_[cl_size_max * n_t_noalloc];
  cs_lnum_t *partial_count = partial_count_;
  if (n_threads > n_t_noalloc)
    CS_MALLOC(partial_count, n_threads*(stride+1), cs_lnum_t);

  // Compute partial sums (only reading from a, not writing to it)

  #pragma omp parallel shared(partial_count) num_threads(n_threads)
  {
    int t_id = omp_get_thread_num();
    cs_lnum_t t_s_id, t_e_id;
    cs_parall_thread_range(n, sizeof(cs_lnum_t), t_id, n_threads,
                           &t_s_id, &t_e_id);

    cs_lnum_t l_count = 0;
    for (cs_lnum_t i = t_s_id; i < t_e_id; i++) {
      cs_lnum_t ai = a[i];
      if (ai > c) {
        a[t_s_id + l_count] = ai;
        l_count += 1;
      }
    }

    partial_count[t_id*stride] = l_count;
  }

  // Serial inclusive scan of partial counts.

  {
    n_select = partial_count[0];
    for (int j = 1; j < n_threads; j++) {
      n_select += partial_count[j*stride];
      partial_count[j*stride] = n_select;
    }
  }

  // Copy portions of array which may be overlapped with memmove semantics.

  {
    for (int t_id = 1; t_id < n_threads; t_id++) {
      cs_lnum_t s_id = partial_count[(t_id-1)*stride];
      cs_lnum_t e_id = partial_count[t_id*stride];
      cs_lnum_t n_t_select = (e_id - s_id);
      cs_lnum_t t_s_id, t_e_id;
      cs_parall_thread_range(n, sizeof(cs_lnum_t), t_id, n_threads,
                             &t_s_id, &t_e_id);
      t_e_id = t_s_id + n_t_select;
      if (t_s_id > n_select)  // No overlap, do this in next parallel pass.
        break;
      if (t_e_id > n_select) {
        e_id -= (t_e_id - n_select);
        assert(e_id >= s_id);
      }
      if (e_id > s_id > 0) {
        size_t nb = (e_id-s_id)*sizeof(cs_lnum_t);
        cs_lnum_t *dst = &a[s_id];
        cs_lnum_t *src = &a[t_s_id];
        memmove(dst, src, nb);
      }
    }
  }

  // Now copy portions which may not be overlapped with memcpy semantics.

  #pragma omp parallel shared(partial_count) num_threads(n_threads)
  {
    int t_id = omp_get_thread_num();
    if (t_id > 0) {
      cs_lnum_t s_id = partial_count[(t_id-1)*stride];
      cs_lnum_t e_id = partial_count[t_id*stride];
      cs_lnum_t n_t_select = (e_id - s_id);
      cs_lnum_t t_s_id, t_e_id;
      cs_parall_thread_range(n, sizeof(cs_lnum_t), t_id, n_threads,
                             &t_s_id, &t_e_id);
      t_e_id = t_s_id + n_t_select;
      if (t_e_id > n_select) {
        if (t_s_id < n_select) {
          // A part has already been handled by memmove.
          cs_lnum_t n_rem = (t_e_id - n_select);
          cs_lnum_t n_shift = (e_id - n_rem) - s_id;
          s_id = e_id - n_rem;
          assert(s_id <= n_select);
          t_s_id += n_shift;
        }
        if (e_id > s_id > 0) {
          size_t nb = (e_id-s_id)*sizeof(cs_lnum_t);
          cs_lnum_t *dst = &a[s_id];
          cs_lnum_t *src = &a[t_s_id];
          memcpy(dst, src, nb);
        }
      }
    }
  }

  if (partial_count != partial_count_)
    CS_FREE(partial_count);

  return n_select;
}

#endif // defined(HAVE_OPENMP)

#if defined(HAVE_CUDA)

/*--------------------------------------------------------------------------*/
/*!
 * \brief Transform a count to an index in place, CUDA (cub) version.
 *
 * For n input elements, the array size should be size n+1, to account
 * for the past-the-end count.
 *
 * If temporary storage is provided by the caller, it will be used
 * it its size is sufficient, avoiding the overhead of local
 * memory allocation.
 *
 * \param[in]       stream              associated stream
 * \param[in]       n                   number of elements
 * \param[in, out]  a                   count in, index out (size: n+1)
 * \param           tmp_size_caller     size of storage provided by caller
 * \param           tmp_storage_caller  storage provided by caller
 */
/*--------------------------------------------------------------------------*/

static void
_count_to_index_inplace_cuda(cudaStream_t  stream,
                             cs_lnum_t     n,
                             cs_lnum_t     a[],
                             size_t        tmp_size_caller,
                             void         *tmp_storage_caller)
{
  // Ensure past-the-end value is initialized.
  CS_CUDA_CHECK(cudaMemsetAsync(&(a[n]), 0, sizeof(cs_lnum_t), stream));
  CS_CUDA_CHECK(cudaGetLastError());

  assert(a != nullptr);

  unsigned char *tmp_storage
    = reinterpret_cast<unsigned char *>(tmp_storage_caller);
  size_t tmp_storage_size = 0;

  CS_CUDA_CHECK(cub::DeviceScan::ExclusiveSum(nullptr, tmp_storage_size,
                                              a, a, n+1, stream));
  CS_CUDA_CHECK(cudaGetLastError());

  if (tmp_storage_size > tmp_size_caller)
    CS_MALLOC_HD(tmp_storage, tmp_storage_size, unsigned char, CS_ALLOC_DEVICE);

  CS_CUDA_CHECK(cub::DeviceScan::ExclusiveSum(tmp_storage, tmp_storage_size,
                                              a, a, n+1, stream));
  CS_CUDA_CHECK(cudaGetLastError());

  CS_CUDA_CHECK(cub::SyncStream(stream));
  CS_CUDA_CHECK(cudaGetLastError());

  if (tmp_storage != tmp_storage_caller)
    CS_FREE(tmp_storage);
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief Select and compact elements whose values are greater than a
 *        given number, CUDA (cub) version.
 *
 * If temporary storage is provided by the caller, it will be used
 * it its size is sufficient, avoiding the overhead of local
 * memory allocation.
 *
 * \param[in]       stream              associated stream
 * \param[in]       n                   number of elements
 * \param[in]       c                   cutoff value to compare to
 * \param[in, out]  a                   elements in, selected elements out
 * \param           tmp_size_caller     size of storage provided by caller
 * \param           tmp_storage_caller  storage provided by caller
 *
 * \return  number of selected elements
 */
/*--------------------------------------------------------------------------*/

static cs_lnum_t
_select_if_gt_cuda(cudaStream_t  stream,
                   cs_lnum_t     n,
                   cs_lnum_t     c,
                   cs_lnum_t     a[],
                   size_t        tmp_size_caller,
                   void         *tmp_storage_caller)
{
  assert(a != nullptr);

  unsigned char *tmp_storage
    = reinterpret_cast<unsigned char *>(tmp_storage_caller);
  size_t tmp_storage_size = 0;

  // Use existing device to host buffers to avoid managment overhead.
  cs_lnum_t *r_device, *r_host;
  {
    int64_t *r_grid;
    int stream_id = cs_cuda_get_stream_id(stream);
    if (stream_id < 0)
      stream_id = 0;
    cs_cuda_get_2_stage_reduce_buffers
      (stream_id, 1, sizeof(cs_lnum_t), 1,
       (void *&)r_grid, (void *&)r_device, (void *&)r_host);
  }

  cs_algorithm_greater_than select_op(c);

  CS_CUDA_CHECK(cub::DeviceSelect::If(nullptr, tmp_storage_size,
                                      a, r_device, n, select_op, stream));
  CS_CUDA_CHECK(cudaGetLastError());

  if (tmp_storage_size > tmp_size_caller)
    CS_MALLOC_HD(tmp_storage, tmp_storage_size, unsigned char, CS_ALLOC_DEVICE);

  CS_CUDA_CHECK(cub::DeviceSelect::If(tmp_storage, tmp_storage_size,
                                      a, r_device, n, select_op, stream));
  CS_CUDA_CHECK(cudaGetLastError());

  CS_CUDA_CHECK(cudaMemcpyAsync(r_host, r_device, sizeof(cs_lnum_t),
                                cudaMemcpyDeviceToHost, stream));

  CS_CUDA_CHECK(cudaStreamSynchronize(stream));

  if (tmp_storage != tmp_storage_caller)
    CS_FREE(tmp_storage);

  return *r_host;
}

#endif

#if defined(HAVE_HIP)

/*--------------------------------------------------------------------------*/
/*!
 * \brief Transform a count to an index in place, ROCm (rocPRIM) version.
 *
 * For n input elements, the array size should be size n+1, to account
 * for the past-the-end count.
 *
 * \param[in]       n     number of elements
 * \param[in, out]  a <-> count in, index out (size: n+1)
 */
/*--------------------------------------------------------------------------*/

static void
_count_to_index_inplace_rocm(hipStream_t  stream,
                             cs_lnum_t    n,
                             cs_lnum_t    a[])
{
  // Ensure past-the-end value is initialized.
  hipMemsetAsync(&(a[n]), 0, sizeof(cs_lnum_t), stream);

  assert(a != nullptr && n > 0);

  unsigned char *tmp_storage = nullptr;
  size_t tmp_storage_size = 0;

  rocprim::exclusive_scan(tmp_storage, tmp_storage_size,
                          a, a, 0, n+1, rocprim::plus<cs_lnum_t>(),
                          stream);

  CS_MALLOC_HD(tmp_storage, tmp_storage_size, unsigned char, CS_ALLOC_DEVICE);

  rocprim::exclusive_scan(tmp_storage, tmp_storage_size,
                          a, a, 0, n+1, rocprim::plus<cs_lnum_t>(),
                          stream);

  cub::SyncStream(stream);

  CS_FREE(tmp_storage);
}

#endif

/*--------------------------------------------------------------------------*/
/*!
 * \brief Transform a count to an index in place, serial version.
 *
 * For n input elements, the array size should be size n+1, to account
 * for the past-the-end count.
 *
 * \param[in]       n     number of elements
 * \param[in, out]  a <-> count in, index out (size: n+1)
 */
/*--------------------------------------------------------------------------*/

static void
_count_to_index_inplace_serial(cs_lnum_t  n,
                               cs_lnum_t  a[])
{
  assert(a != nullptr);

  cs_lnum_t s = 0;

  for (cs_lnum_t i = 0; i < n; i++) {
    cs_lnum_t c = a[i];
    a[i] = s;
    s += c;
  }

  a[n] = s;
}

/*--------------------------------------------------------------------------*/
/*!
 * \brief Select and compact elements whose values are greater than a
 *        given number, single-threaded version.
 *
 * \param[in]       n            number of elements
 * \param[in]       c            cutoff value to compare to
 * \param[in, out]  a            elements in, selecetd elements out
 *
 * \return  number of selected elements
 */
/*--------------------------------------------------------------------------*/

static cs_lnum_t
_select_if_gt_serial(cs_lnum_t  n,
                     cs_lnum_t  c,
                     cs_lnum_t  a[])
{
  cs_lnum_t n_select = 0;

  assert(a != nullptr && n > 0);

  // Compute partial sums (only reading from a, not writing to it)

  for (cs_lnum_t i = 0; i < n; i++) {
    cs_lnum_t ai = a[i];
    if (ai > c) {
      a[n_select] = ai;
      n_select += 1;
    }
  }

  return n_select;
}

namespace cs {
namespace algorithm {

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*--------------------------------------------------------------------------*/
/*
 * \brief Transform a count to an index in-place.
 *
 * For n input elements, the array size should be size n+1, to account
 * for the past-the-end count.
 *
 * If temporary storage is provided by the caller, it will be used
 * it its size is sufficient, avoiding the overhead of local
 * memory allocation. This is useful only for device code.
 *
 * \param[in]       ctx          associated dispatch context
 * \param[in]       n            number of elements
 * \param[in, out]  a            count in, index out (size: n+1)
 * \param           tmp_size     optional temporary memory size
 * \param           tmp_storage  optional temporary memory
 */
/*--------------------------------------------------------------------------*/

void
count_to_index(cs_dispatch_context       &ctx,
               cs_lnum_t                  n,
               cs_lnum_t                  a[],
               [[maybe_unused]]  size_t   tmp_size,
               [[maybe_unused]]  void    *tmp_storage)
{
#if defined(HAVE_CUDA)
  if (ctx.use_gpu()) {
    _count_to_index_inplace_cuda(ctx.stream(), n, a,
                                 tmp_size, tmp_storage);
    return;
  }
#endif

#if defined(HAVE_HIP)
  if (ctx.use_gpu()) {
    _count_to_index_inplace_rocm(ctx.stream(), n, a,
                                 tmp_size, tmp_storage);
    return;
  }
#endif

#if defined(HAVE_OPENMP)

  int n_threads = ctx.n_cpu_threads();
  if (n_threads == -1)
    n_threads = cs_parall_n_threads(n, CS_THR_MIN);

  if (n_threads > 1) {
    _count_to_index_inplace_omp(n_threads, n, a);
    return;
  }

#endif

  _count_to_index_inplace_serial(n, a);
}

/*--------------------------------------------------------------------------*/
/*
 * Select and compact elements whose values are greater than a given number.
 *
 * If temporary storage is provided by the caller, it will be used
 * it its size is sufficient, avoiding the overhead of local
 * memory allocation. This is useful only for device code.
 *
 * \param[in]       ctx          associated dispatch context
 * \param[in]       n            number of elements
 * \param[in]       c            cutoff value to compare to
 * \param[in, out]  a            elements in, selected elements out
 * \param           tmp_size     optional temporary memory size
 * \param           tmp_storage  optional temporary memory
 *
 * \return  number of selected elements
 */
/*--------------------------------------------------------------------------*/

cs_lnum_t
select_if_gt(cs_dispatch_context      &ctx,
             cs_lnum_t                 n,
             cs_lnum_t                 c,
             cs_lnum_t                 a[],
             [[maybe_unused]] size_t   tmp_size,
             [[maybe_unused]] void    *tmp_storage)
{
#if defined(HAVE_CUDA)
  if (ctx.use_gpu()) {
    return _select_if_gt_cuda(ctx.stream(), n, c, a,
                              tmp_size, tmp_storage);
  }
#endif

#if defined(HAVE_HIP)
  if (ctx.use_gpu()) {
    // TODO: call GPU code for this and return;
    // Currently defaults to OpenMP (requiring UVM).
  }
#endif

#if defined(HAVE_OPENMP)

  int n_threads = ctx.n_cpu_threads();
  if (n_threads == -1)
    n_threads = cs_parall_n_threads(n, CS_THR_MIN);

  if (n_threads > 1) {
    return _select_if_gt_omp(n_threads, n, c, a);
  }

#endif

  return _select_if_gt_serial(n, c, a);
}

/*--------------------------------------------------------------------------*/
/*
 * \brief Sum a local counter.
 *
 * \param[in]  ctx      associated dispatch context
 * \param[in]  n        number of elements
 * \param[in]  counter  local counter
 *
 * \return sum of local counter values
 */
/*--------------------------------------------------------------------------*/

cs_gnum_t
count_reduce_sum(cs_dispatch_context  &ctx,
                 cs_lnum_t             n,
                 short                 counter[])
{
  cs_gnum_t count = 0;

  ctx.parallel_for_reduce_sum
    (n, count, [=] CS_F_HOST_DEVICE
     (cs_lnum_t i, CS_DISPATCH_REDUCER_TYPE(cs_gnum_t) &sum) {
      sum += (cs_gnum_t)counter[i];
    });
  ctx.wait();

  return count;
}

/*----------------------------------------------------------------------------*/

} // namespace algorithm
} // namespace cs
