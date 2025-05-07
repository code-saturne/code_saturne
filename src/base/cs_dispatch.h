#ifndef CS_DISPATCH_H
#define CS_DISPATCH_H

/*============================================================================
 * Class to dispatch computation using various runtimes (OpenMP, CUDA, ...)
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

// Valid only for C++

#ifdef __cplusplus

/*----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C++ library headers
 *----------------------------------------------------------------------------*/

#include <utility>
#include <cmath>

#if defined(SYCL_LANGUAGE_VERSION)
#include <sycl/sycl.hpp>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_assert.h"
#include "base/cs_mem.h"

#ifdef __CUDACC__
#include "base/cs_base_cuda.h"
#include "base/cs_cuda_reduce.h"
#include "cs_math_cuda.cuh"
#endif

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_dispatch.h
        Dispatch computation using various runtimes (OpenMP, CUDA, ...)
*/

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#if defined(SYCL_LANGUAGE_VERSION)

#define CS_DISPATCH_REDUCER_TYPE(type) auto

#else

#define CS_DISPATCH_REDUCER_TYPE(type) type

#endif

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Algorithm classes for indirect sums from a graph's edges to its nodes,
  such as face-to cell sums. */

typedef enum {

  CS_DISPATCH_SUM_SIMPLE,  /*!< Simple sum (assumes data race-avoiding
                                numbering/coloring) */
  CS_DISPATCH_SUM_ATOMIC   /*!< Atomic sum */

} cs_dispatch_sum_type_t;

/*!
 * Provide default implementations of a cs_context based on parallel_for
 * function. This class is a mixin that use CRTP (Curiously Recurring
 * Template Pattern) to provide such functions.
 */

template <class Derived>
class cs_dispatch_context_mixin {
public:

  // Loop over n elements
  // Must be redefined by the child class
  template <class F, class... Args>
  decltype(auto)
  parallel_for(cs_lnum_t n, F&& f, Args&&... args) = delete;

  // Assembly loop over all internal faces
  template <class M, class F, class... Args>
  decltype(auto)
  parallel_for_i_faces(const M*          m,
                       F&&               f,
                       Args&&...         args);

  // Assembly loop over all boundary faces
  template <class M, class F, class... Args>
  decltype(auto)
  parallel_for_b_faces(const M*          m,
                       F&&               f,
                       Args&&...         args);

  // Parallel reduction with simple sum.
  // Must be redefined by the child class
  template <class T, class F, class... Args>
  decltype(auto)
  parallel_for_reduce_sum
    (cs_lnum_t n, T& sum, F&& f, Args&&... args) = delete;

  // Parallel reduction with reducer template.
  // Must be redefined by the child class
  template <class T, class R, class F, class... Args>
  decltype(auto)
  parallel_for_reduce
    (cs_lnum_t n, T& r, R& reducer, F&& f, Args&&... args) = delete;

  // Wait upon completion
  // Must be redefined by the child class
  template <class... Args>
  decltype(auto)
  wait(void) = delete;

  // Query sum type for assembly loop over all interior faces
  // Must be redefined by the child class
  template <class M>
  bool
  try_get_parallel_for_i_faces_sum_type(const M*                  m,
                                        cs_dispatch_sum_type_t&  st);

  // Query sum type for assembly loop over all boundary faces
  // Must be redefined by the child class
  template <class M>
  bool
  try_get_parallel_for_b_faces_sum_type(const M*                 m,
                                        cs_dispatch_sum_type_t&  st);

};

// Default implementation of parallel_for_i_faces based on parallel_for
template <class Derived>
template <class M, class F, class... Args>
decltype(auto) cs_dispatch_context_mixin<Derived>::parallel_for_i_faces
  (const M* m, F&& f, Args&&... args) {
  return static_cast<Derived*>(this)->parallel_for
                                        (m->n_i_faces,
                                         static_cast<F&&>(f),
                                         static_cast<Args&&>(args)...);
}

// Default implementation of parallel_for_b_faces based on parallel_for_sum
template <class Derived>
template <class M, class F, class... Args>
decltype(auto) cs_dispatch_context_mixin<Derived>::parallel_for_b_faces
  (const M* m, F&& f, Args&&... args) {
  return static_cast<Derived*>(this)->parallel_for
                                        (m->n_b_faces,
                                         static_cast<F&&>(f),
                                         static_cast<Args&&>(args)...);
}

// Default implementation of get interior faces sum type
template <class Derived>
template <class M>
bool cs_dispatch_context_mixin<Derived>::try_get_parallel_for_i_faces_sum_type
  ([[maybe_unused]]const M*          m,
   cs_dispatch_sum_type_t&           st) {
  st = CS_DISPATCH_SUM_SIMPLE;
  return true;
}

// Default implementation of get boundary faces sum type
template <class Derived>
template <class M>
bool cs_dispatch_context_mixin<Derived>::try_get_parallel_for_b_faces_sum_type
  ([[maybe_unused]]const M*          m,
   cs_dispatch_sum_type_t&           st) {
  st = CS_DISPATCH_SUM_SIMPLE;
  return true;
}

/*
 * cs_context to execute loops with OpenMP on the CPU
 */

class cs_host_context : public cs_dispatch_context_mixin<cs_host_context> {

private:

  cs_lnum_t  n_min_per_thread;  /*!< Minimum number of elements per thread */
  int        n_threads_;        /*!< If defined, force number of threads */

public:

  cs_host_context()
    : n_min_per_thread(CS_THR_MIN), n_threads_(-1)
  {}

  //!  Compute array index bounds for a local thread.
  //
  // \param[in]       n          size of array
  // \param[in]       type_size  element type size (or multiple)
  // \param[in, out]  s_id       start index for the current thread
  // \param[in, out]  e_id       past-the-end index for the current thread

private:

#if defined(HAVE_OPENMP)

  // Determine number of threads that should actually be used.
  int
  n_threads(cs_lnum_t   n)
  {
    int n_t = n_threads_;
    if (n_t < 0) {
      n_t = cs_glob_n_threads;
      int n_t_l = n / n_min_per_thread;
      if (n_t_l < n_t)
        n_t = n_t_l;
      if (n_t < 1)
        n_t = 1;
    }
    return n_t;
  }

#endif

  // Determine element range for current thread.
  void
  thread_range(cs_lnum_t                n,
               [[maybe_unused]] size_t  type_size,
               cs_lnum_t                &s_id,
               cs_lnum_t                &e_id)
  {
#if defined(HAVE_OPENMP)
    const int t_id = omp_get_thread_num();
    const int n_t = omp_get_num_threads();
    const cs_lnum_t t_n = (n + n_t - 1) / n_t;
    const cs_lnum_t cl_m = CS_CL_SIZE / type_size;  /* Cache line multiple */

    s_id =  t_id    * t_n;
    e_id = (t_id+1) * t_n;
    s_id = cs_align(s_id, cl_m);
    e_id = cs_align(e_id, cl_m);
    if (e_id > n) e_id = n;
#else
    s_id = 0;
    e_id = n;
#endif
  }

public:

  //! Set minimum number of elements threshold for CPU multithread execution.
  void
  set_n_min_per_cpu_thread(cs_lnum_t  n) {
    this->n_min_per_thread = n;
  }

  //! Get minimum number of elements threshold for CPU multithread execution.
  cs_lnum_t
  n_min_per_cpu_thread(void) {
    return this->n_min_per_thread;
  }

  //! Set number of threads for CPU multithread execution.
  void
  set_n_cpu_threads(int  n) {
    this->n_threads_ = n;
  }

  //! Get number of threads for CPU multithread execution (-1 if automatic)
  int
  n_cpu_threads(void) {
    return this->n_threads_;
  }

  //! Iterate using a plain omp parallel for
  template <class F, class... Args>
  bool
  parallel_for(cs_lnum_t n, F&& f, Args&&... args) {
#   pragma omp parallel for  num_threads(n_threads(n))
    for (cs_lnum_t i = 0; i < n; ++i) {
      f(i, args...);
    }
    return true;
  }

  //! Loop over the interior faces of a mesh using a specific numbering
  //! that avoids conflicts between threads.
  template <class M, class F, class... Args>
  bool
  parallel_for_i_faces(const M* m, F&& f, Args&&... args) {
    const int n_i_groups  = m->i_face_numbering->n_groups;
    const int n_i_threads = m->i_face_numbering->n_threads;
    const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;
    for (int g_id = 0; g_id < n_i_groups; g_id++) {
      #pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t f_id = i_group_index[(t_id * n_i_groups + g_id) * 2];
             f_id < i_group_index[(t_id * n_i_groups + g_id) * 2 + 1];
             f_id++) {
          f(f_id, args...);
        }
      }
    }
    return true;
  }

  //! Loop over the boundary faces of a mesh using a specific numbering
  //! that avoids conflicts between threads.
  template <class M, class F, class... Args>
  bool
  parallel_for_b_faces(const M* m, F&& f, Args&&... args) {
    const int n_b_threads = m->b_face_numbering->n_threads;
    const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

    #pragma omp parallel for
    for (int t_id = 0; t_id < n_b_threads; t_id++) {
      for (cs_lnum_t f_id = b_group_index[t_id*2];
           f_id < b_group_index[t_id*2 + 1];
           f_id++) {
        f(f_id, args...);
      }
    }
    return true;
  }

  //! Plain OpenMP parallel reduction with simple sum.
  template <class T, class F, class... Args>
  bool
  parallel_for_reduce_sum(cs_lnum_t n,
                          T&        sum,
                          F&&       f,
                          Args&&... args) {
    sum = 0;

#if 0
#   pragma omp parallel for reduction(+:sum) num_threads(n_threads(n))
    for (cs_lnum_t i = 0; i < n; ++i) {
      f(i, sum, args...);
    }
#else
#   pragma omp parallel num_threads(n_threads(n))
    {
      cs_lnum_t s_id, e_id;
      thread_range(n, 4, s_id, e_id);

      const cs_lnum_t _n = e_id - s_id;
      cs_lnum_t n_sblocks, blocks_in_sblocks;
      const cs_lnum_t block_size = 60;
      { // superblock counts
        cs_lnum_t n_blocks = (_n + block_size - 1) / block_size;
        n_sblocks = (n_blocks > 1) ? std::sqrt(n_blocks) : 1;
        cs_lnum_t n_b = block_size * n_sblocks;
        blocks_in_sblocks = (_n + n_b - 1) / n_b;
      }

      for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {
        T sum_sblock = 0;

        for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
          cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid) + s_id;
          cs_lnum_t end_id = start_id + block_size;
          if (end_id > e_id)
            end_id = e_id;
          T sum_block = 0;
          for (cs_lnum_t i = start_id; i < end_id; i++) {
            f(i, sum_block, args...);
          }
          sum_sblock += sum_block;
        }

        #pragma omp atomic
        sum += sum_sblock;
      }
    }

#endif
    return true;
  }

  //! OpenMP parallel reduction with general reducer.
  // In case the reduction involves floating-point sums,
  // we use a Superblock / block loop to reduce numerical error.
  template <class T, class R, class F, class... Args>
  bool
  parallel_for_reduce(cs_lnum_t n,
                      T&        result,
                      R&        reducer,
                      F&&       f,
                      Args&&... args) {
    reducer.identity(result);

#   pragma omp parallel num_threads(n_threads(n))
    {
      cs_lnum_t s_id, e_id;
      thread_range(n, 4, s_id, e_id);

      const cs_lnum_t _n = e_id - s_id;
      cs_lnum_t n_sblocks, blocks_in_sblocks;
      const cs_lnum_t block_size = 60;
      { // superblock counts
        cs_lnum_t n_blocks = (_n + block_size - 1) / block_size;
        n_sblocks = (n_blocks > 1) ? std::sqrt(n_blocks) : 1;
        cs_lnum_t n_b = block_size * n_sblocks;
        blocks_in_sblocks = (_n + n_b - 1) / n_b;
      }

      for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {
        T result_sblock;
        reducer.identity(result_sblock);

        for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
          cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid) + s_id;
          cs_lnum_t end_id = start_id + block_size;
          if (end_id > e_id)
            end_id = e_id;
          T result_block;
          reducer.identity(result_block);
          for (cs_lnum_t i = start_id; i < end_id; i++) {
            f(i, result_block, args...);
            reducer.combine(result_sblock, result_block);
          }
        }

        #pragma omp critical
        {
          reducer.combine(result, result_sblock);
        }
      }
    }

    return true;
  }

  //! Wait upon completion
  // No-op here as Open-MP based methods used here have implicit barriers.
  template <class... Args>
  bool
  wait(void) {
    return true;
  }

  // Get interior faces sum type associated with this context
  template <class M>
  bool
  try_get_parallel_for_i_faces_sum_type([[maybe_unused]]const M*   m,
                                        cs_dispatch_sum_type_t&    st) {
    st = CS_DISPATCH_SUM_SIMPLE;
    return true;
  }

  // Get boundary faces sum type associated with this context
  template <class M>
  bool
  try_get_parallel_for_b_faces_sum_type([[maybe_unused]]const M*   m,
                                        cs_dispatch_sum_type_t&    st) {
    st = CS_DISPATCH_SUM_SIMPLE;
    return true;
  }

};

#if defined(__CUDACC__)

/* Default kernel that loops over an integer range and calls a device functor.
   This kernel uses a grid_size-stride loop and thus guarantees that all
   integers are processed, even if the grid is smaller.
   All arguments *must* be passed by value to avoid passing CPU references
   to the GPU. */

template <class F, class... Args>
__global__ void cs_cuda_kernel_parallel_for(cs_lnum_t n, F f, Args... args) {
  // grid_size-stride loop
  for (cs_lnum_t id = blockIdx.x * blockDim.x + threadIdx.x; id < n;
       id += blockDim.x * gridDim.x) {
    f(id, args...);
  }
}

/* Default kernel that loops over an integer range and calls a device functor,
   also reducing a sum over all elements.
   This kernel uses a grid_size-stride loop and thus guarantees that all
   integers are processed, even if the grid is smaller.
   All arguments *must* be passed by value to avoid passing CPU references
   to the GPU. */

template <class T, class F, class... Args>
__global__ void
cs_cuda_kernel_parallel_for_reduce_sum(cs_lnum_t   n,
                                       T          *b_res,
                                       F           f,
                                       Args...     args) {
  // grid_size-stride loop
  extern __shared__  int p_stmp[];
  T *stmp = reinterpret_cast<T *>(p_stmp);
  const cs_lnum_t tid = threadIdx.x;

  stmp[tid] = 0;

  for (cs_lnum_t id = blockIdx.x * blockDim.x + threadIdx.x; id < n;
       id += blockDim.x * gridDim.x) {
     f(id, stmp[tid], args...);
  }

  switch (blockDim.x) {
  case 1024:
    cs_cuda_reduce_block_reduce_sum<1024, 1>(stmp, tid, b_res);
    break;
  case 512:
    cs_cuda_reduce_block_reduce_sum<512, 1>(stmp, tid, b_res);
    break;
  case 256:
    cs_cuda_reduce_block_reduce_sum<256, 1>(stmp, tid, b_res);
    break;
  case 128:
    cs_cuda_reduce_block_reduce_sum<128, 1>(stmp, tid, b_res);
    break;
  default:
    assert(0);
  }
}

/* Default kernel that loops over an integer range and calls a device functor.
   also computing a reduction over all elements.
   This kernel uses a grid_size-stride loop and thus guarantees that all
   integers are processed, even if the grid is smaller.
   All arguments *must* be passed by value to avoid passing CPU references
   to the GPU. */

template <class T, class R, class F, class... Args>
__global__ void
cs_cuda_kernel_parallel_for_reduce(cs_lnum_t   n,
                                   T          *b_res,
                                   R          &reducer,
                                   F           f,
                                   Args...     args) {
  // grid_size-stride loop
  extern __shared__  int p_stmp[];
  T *stmp = reinterpret_cast<T *>(p_stmp);
  const cs_lnum_t tid = threadIdx.x;

  reducer.identity(stmp[tid]);

  for (cs_lnum_t id = blockIdx.x * blockDim.x + threadIdx.x; id < n;
       id += blockDim.x * gridDim.x) {
    T rd;
    /* It would be safer to call reducer.identyity() here in case all
       values of rd are not set for each thread, but this might incurr
       a small performance penalty, and is redundant in most cases,
       so we consider all values of rd must be set by the caller. */
    // reducer.identity(rd);
    f(id, rd, args...);
    stmp[tid] = rd;
  }

  switch (blockDim.x) {
  case 1024:
    cs_cuda_reduce_block_reduce<1024, R>(stmp, tid, b_res);
    break;
  case 512:
    cs_cuda_reduce_block_reduce<512, R>(stmp, tid, b_res);
    break;
  case 256:
    cs_cuda_reduce_block_reduce<256, R>(stmp, tid, b_res);
    break;
  case 128:
    cs_cuda_reduce_block_reduce<128, R>(stmp, tid, b_res);
    break;
  default:
    assert(0);
  }
}

/*!
 * Context to execute loops with CUDA on the device
 */

class cs_device_context : public cs_dispatch_context_mixin<cs_device_context> {

private:

  long          grid_size_;   /*!< Associated grid size; if <= 0, each kernel
                                launch will use a grid size based on
                                the number of elements. */
  long          block_size_;  /*!< Associated block size */
  cudaStream_t  stream_;      /*!< Associated CUDA stream */
  int           device_;      /*!< Associated CUDA device id */

  bool          use_gpu_;     /*!< Run on GPU if available */

public:

  //! Constructor

  cs_device_context(void)
    : grid_size_(0), block_size_(256), stream_(cs_cuda_get_stream(0)),
      device_(0), use_gpu_(true)
  {
    device_ = cs_glob_cuda_device_id;
  }

  cs_device_context(long          grid_size,
                    long          block_size,
                    cudaStream_t  stream,
                    int           device)
    : grid_size_(grid_size), block_size_(block_size), stream_(stream),
      device_(device), use_gpu_(true)
  {}

  cs_device_context(long          grid_size,
                    long          block_size,
                    cudaStream_t  stream)
    : grid_size_(grid_size), block_size_(block_size), stream_(stream),
      device_(0), use_gpu_(true)
  {
    device_ = cs_base_cuda_get_device();
  }

  cs_device_context(long  grid_size,
                    long  block_size)
    : grid_size_(grid_size), block_size_(block_size),
       stream_(cs_cuda_get_stream(0)), device_(0), use_gpu_(true)
  {
    device_ = cs_base_cuda_get_device();
  }

  cs_device_context(cudaStream_t  stream)
    : grid_size_(0), block_size_(256), stream_(stream), device_(0),
      use_gpu_(true)
  {
    device_ = cs_base_cuda_get_device();
  }

  //! Change grid_size configuration, but keep the stream and device
  //
  // \param[in]  grid_size   CUDA grid size, or -1 for automatic choice
  // \param[in]  block_size  CUDA block size (power of 2 if reduction is used)

  void
  set_cuda_grid(long  grid_size,
                long  block_size) {
    this->grid_size_ = (grid_size > 0) ? grid_size : -1;
    this->block_size_ = block_size;
  }

  //! Change stream, but keep the grid and device configuration

  void
  set_cuda_stream(cudaStream_t stream) {
    this->stream_ = stream;
  }

  //! Change stream, but keep the grid and device configuration

  void
  set_cuda_stream(int  stream_id) {
    this->stream_ = cs_cuda_get_stream(stream_id);
  }

  //! Get associated stream

  cudaStream_t
  cuda_stream(void) {
    return this->stream_;
  }

  //! Change CUDA device

  void
  set_cuda_device(int  device) {
    this->device_ = device;
  }

  //! Set or unset execution on GPU

  void
  set_use_gpu(bool  use_gpu) {
    this->use_gpu_ = use_gpu;
  }

  //! Check whether we are trying to run on GPU

  bool
  use_gpu(void) {
    return (device_ >= 0 && use_gpu_);
  }

  //! Check preferred allocation mode depending on execution policy

  cs_alloc_mode_t
  alloc_mode(void) {
    cs_alloc_mode_t amode
      = (device_ >= 0 && use_gpu_) ? CS_ALLOC_DEVICE : CS_ALLOC_HOST;
    return (amode);
  }

  cs_alloc_mode_t
  alloc_mode(bool readable_on_cpu) {
    cs_alloc_mode_t amode = CS_ALLOC_HOST;
    if (device_ >= 0 && use_gpu_) {
      if (readable_on_cpu)
        amode = CS_ALLOC_HOST_DEVICE_SHARED;
      else
        amode = CS_ALLOC_DEVICE;
    }
    return (amode);
  }

public:

  //! Try to launch on the GPU and return false if not available
  template <class F, class... Args>
  bool
  parallel_for(cs_lnum_t n, F&& f, Args&&... args) {
    if (device_ < 0 || use_gpu_ == false) {
      return false;
    }

    long l_grid_size = grid_size_;
    if (l_grid_size < 1) {
      l_grid_size = (n % block_size_) ? n/block_size_ + 1 : n/block_size_;
    }

    if (n > 0)
      cs_cuda_kernel_parallel_for<<<l_grid_size, block_size_, 0, stream_>>>
        (n, static_cast<F&&>(f), static_cast<Args&&>(args)...);

    return true;
  }

  //! Try to launch on the GPU and return false if not available
  template <class M, class F, class... Args>
  bool
  parallel_for_i_faces(const M* m, F&& f, Args&&... args) {
    const cs_lnum_t n = m->n_i_faces;
    if (device_ < 0 || use_gpu_ == false) {
      return false;
    }

    long l_grid_size = grid_size_;
    if (l_grid_size < 1) {
      l_grid_size = (n % block_size_) ? n/block_size_ + 1 : n/block_size_;
    }

    if (n > 0)
      cs_cuda_kernel_parallel_for<<<l_grid_size, block_size_, 0, stream_>>>
        (n, static_cast<F&&>(f), static_cast<Args&&>(args)...);

    return true;
  }

  //! Launch kernel on the GPU with simple sum reduction
  //! The reduction involves an implicit wait().
  template <class T, class F, class... Args>
  bool
  parallel_for_reduce_sum(cs_lnum_t n,
                          T&        sum,
                          F&&       f,
                          Args&&... args) {
    if (device_ < 0 || use_gpu_ == false) {
      return false;
    }

    sum = 0;

    long l_grid_size = grid_size_;
    if (l_grid_size < 1) {
      l_grid_size = (n % block_size_) ? n/block_size_ + 1 : n/block_size_;
    }
    if (n == 0) {
      return true;
    }

    int stream_id = cs_cuda_get_stream_id(stream_);
    if (stream_id < 0)
      stream_id = 0;

    T *r_grid_, *r_reduce_;
    cs_cuda_get_2_stage_reduce_buffers
      (stream_id, n, sizeof(sum), l_grid_size,
       (void *&)r_grid_, (void *&)r_reduce_);

    int smem_size = block_size_ * sizeof(T);
    cs_cuda_kernel_parallel_for_reduce_sum
      <<<l_grid_size, block_size_, smem_size, stream_>>>
      (n, r_grid_, static_cast<F&&>(f), static_cast<Args&&>(args)...);

#if defined(DEBUG) || !defined(NDEBUG)
    cudaError_t retcode = cudaGetLastError();
    if (retcode != cudaSuccess)
      bft_error(__FILE__, __LINE__, 0,
                "[CUDA error] %d: %s\n"
                "with grid size %ld, block size %ld, shared memory size %d.",
                retcode, ::cudaGetErrorString(retcode),
                l_grid_size, block_size_, smem_size);
#endif

    switch (block_size_) {
    case 1024:
      cs_cuda_reduce_sum_single_block<1024, 1>
        <<<1, block_size_, 0, stream_>>>
        (l_grid_size, r_grid_, r_reduce_);
      break;
    case 512:
      cs_cuda_reduce_sum_single_block<512, 1>
        <<<1, block_size_, 0, stream_>>>
        (l_grid_size, r_grid_, r_reduce_);
      break;
    case 256:
      cs_cuda_reduce_sum_single_block<256, 1>
        <<<1, block_size_, 0, stream_>>>
        (l_grid_size, r_grid_, r_reduce_);
      break;
    case 128:
      cs_cuda_reduce_sum_single_block<128, 1>
        <<<1, block_size_, 0, stream_>>>
        (l_grid_size, r_grid_, r_reduce_);
      break;
    default:
      cs_assert(0);
    }

#if defined(DEBUG) || !defined(NDEBUG)
    retcode = cudaGetLastError();
    if (retcode != cudaSuccess)
      bft_error(__FILE__, __LINE__, 0,
                "[CUDA error] %d: %s\n"
                "with grid size %ld, block size %ld, shared memory size %d.",
                retcode, ::cudaGetErrorString(retcode),
                l_grid_size, block_size_, (int)smem_size);
#endif

    CS_CUDA_CHECK(cudaStreamSynchronize(stream_));
    CS_CUDA_CHECK(cudaGetLastError());
    sum = r_reduce_[0];

    return true;
  }

  //! Parallel reduction with general reducer.
  template <class T, class R, class F, class... Args>
  bool
  parallel_for_reduce(cs_lnum_t n,
                      T&        result,
                      R&        reducer,
                      F&&       f,
                      Args&&... args) {
    if (device_ < 0 || use_gpu_ == false) {
      return false;
    }

    reducer.identity(result);

    long l_grid_size = grid_size_;
    if (l_grid_size < 1) {
      l_grid_size = (n % block_size_) ? n/block_size_ + 1 : n/block_size_;
    }
    if (n == 0) {
      return true;
    }

    int stream_id = cs_cuda_get_stream_id(stream_);
    if (stream_id < 0)
      stream_id = 0;

    T *r_grid_, *r_reduce_;
    cs_cuda_get_2_stage_reduce_buffers
      (stream_id, n, sizeof(result), l_grid_size,
       (void *&)r_grid_, (void *&)r_reduce_);

    int l_block_size = block_size_;
    int smem_size = l_block_size * sizeof(T);
    while (smem_size > cs_glob_cuda_shared_mem_per_block) {
      // We should have a runtime failure if even blocks of size 64
      // are too large relative to the available shared memory.
      if (l_block_size < 2)
        bft_error(__FILE__, __LINE__, 0,
                  "Type of size %d exceeds capacity of "
                  "CUDA shared memory (%d).",
                  (int)sizeof(T), cs_glob_cuda_shared_mem_per_block);
      l_block_size /= 2;
      smem_size = l_block_size * sizeof(T);
    }

    cudaError_t retcode = cudaSuccess;

    cs_cuda_kernel_parallel_for_reduce<T, R>
      <<<l_grid_size, l_block_size, smem_size, stream_>>>
      (n, r_grid_, reducer, static_cast<F&&>(f),
       static_cast<Args&&>(args)...);

#if defined(DEBUG) || !defined(NDEBUG)
    retcode = cudaGetLastError();
    if (retcode != cudaSuccess)
      bft_error(__FILE__, __LINE__, 0,
                "[CUDA error] %d: %s\n"
                "with grid size %ld, block size %d, shared memory size %d.",
                retcode, ::cudaGetErrorString(retcode),
                l_grid_size, l_block_size, smem_size);
#endif

    switch (l_block_size) {
    case 1024:
      cs_cuda_reduce_single_block<1024, R>
        <<<1, l_block_size, smem_size, stream_>>>
        (l_grid_size, r_grid_, r_reduce_);
      break;
    case 512:
      cs_cuda_reduce_single_block<512, R>
        <<<1, l_block_size, smem_size, stream_>>>
        (l_grid_size, r_grid_, r_reduce_);
      break;
    case 256:
      cs_cuda_reduce_single_block<256, R>
        <<<1, l_block_size, smem_size, stream_>>>
        (l_grid_size, r_grid_, r_reduce_);
      break;
    case 128:
      cs_cuda_reduce_single_block<128, R>
        <<<1, l_block_size, smem_size, stream_>>>
        (l_grid_size, r_grid_, r_reduce_);
      break;
    default:
      cs_assert(0);
    }

#if defined(DEBUG) || !defined(NDEBUG)
    retcode = cudaGetLastError();
    if (retcode != cudaSuccess)
      bft_error(__FILE__, __LINE__, 0,
                "[CUDA error] %d: %s\n"
                "with grid size %ld, block size %d, shared memory size %d.",
                retcode, ::cudaGetErrorString(retcode),
                l_grid_size, l_block_size, (int)smem_size);
#endif

    CS_CUDA_CHECK(cudaStreamSynchronize(stream_));
    CS_CUDA_CHECK(cudaGetLastError());
    result = r_reduce_[0];

    return true;
  }

  //! Synchronize associated stream
  template <class... Args>
  bool
  wait(void) {
    if (device_ > -1 && use_gpu_) {
      CS_CUDA_CHECK(cudaStreamSynchronize(stream_));
      CS_CUDA_CHECK(cudaGetLastError());
      return true;
    }
    return false;
  }

  // Get interior faces sum type associated with this context
  template <class M>
  bool
  try_get_parallel_for_i_faces_sum_type(const M                 *m,
                                        cs_dispatch_sum_type_t  &st) {
    if (device_ < 0 || use_gpu_ == false) {
      return false;
    }

    st = CS_DISPATCH_SUM_ATOMIC;
    return true;
  }

  // Get boundary faces sum type associated with this context
  template <class M>
  bool
  try_get_parallel_for_b_faces_sum_type(const M                 *m,
                                        cs_dispatch_sum_type_t  &st) {
    if (device_ < 0 || use_gpu_ == false) {
      return false;
    }

    st = CS_DISPATCH_SUM_ATOMIC;
    return true;
  }

};

#elif defined(SYCL_LANGUAGE_VERSION)

/*! Default queue for SYCL */
#if !defined(CS_GLOB_SYCL_QUEUE_IS_DEFINED)
extern sycl::queue  cs_glob_sycl_queue;
#define CS_GLOB_SYCL_QUEUE_IS_DEFINED 1
#endif

/*!
 * Context to execute loops with SYCL on the device
 */

class cs_device_context : public cs_dispatch_context_mixin<cs_device_context> {

private:

  sycl::queue      &queue_;      /*!< Associated SYCL queue */
  bool              is_gpu;      /*!< Is the associated device à GPU ? */

  bool              use_gpu_;    /*!< Run on GPU ? */

public:

  //! Constructor

  cs_device_context(void)
    : queue_(cs_glob_sycl_queue), is_gpu(false), use_gpu_(true)
  {
    is_gpu = queue_.get_device().is_gpu();
  }

  //! Set or unset execution on GPU

  void
  set_use_gpu(bool  use_gpu) {
    this->use_gpu_ = use_gpu;
  }

  //! Check whether we are trying to run on GPU

  bool
  use_gpu(void) {
    return (is_gpu && use_gpu_);
  }

  //! Check preferred allocation mode depending on execution policy

  cs_alloc_mode_t
  alloc_mode(void) {
    cs_alloc_mode_t amode
      = (is_gpu && use_gpu_) ? CS_ALLOC_HOST_DEVICE_SHARED : CS_ALLOC_HOST;
    return (amode);
  }

  cs_alloc_mode_t
  alloc_mode([[maybe_unused]] bool readable_on_cpu) {
    cs_alloc_mode_t amode
      = (is_gpu && use_gpu_) ? CS_ALLOC_HOST_DEVICE_SHARED : CS_ALLOC_HOST;
    return (amode);
  }

public:

  //! Try to launch on the device and return false if not available
  template <class F, class... Args>
  bool
  parallel_for(cs_lnum_t n, F&& f, Args&&... args) {
    if (is_gpu == false || use_gpu_ == false) {
      return false;
    }

    queue_.parallel_for(n, static_cast<F&&>(f), static_cast<Args&&>(args)...);

    return true;
  }

  //! Try to launch on the GPU and return false if not available
  template <class M, class F, class... Args>
  bool
  parallel_for_i_faces(const M* m, F&& f, Args&&... args) {
    const cs_lnum_t n = m->n_i_faces;
    if (is_gpu == false || use_gpu_ == false) {
      return false;
    }

    queue_.parallel_for(n, static_cast<F&&>(f), static_cast<Args&&>(args)...);

    return true;
  }

  //! Launch kernel with simple sum reduction.
  template <class T, class F, class... Args>
  bool
  parallel_for_reduce_sum(cs_lnum_t n,
                          T&        sum_,
                          F&&       f,
                          Args&&... args) {
    if (is_gpu == false || use_gpu_ == false) {
      return false;
    }

    sum_ = 0;

    // TODO: use persistent allocation as we do in CUDA BLAS to avoid
    //       excess allocation/deallocation.
    T *sum_ptr = (T *)sycl::malloc_shared(sizeof(T), queue_);

    queue_.parallel_for(n,
                        sycl::reduction(sum_ptr, (T)0, sycl::plus<T>()),
                        static_cast<F&&>(f),
                        static_cast<Args&&>(args)...).wait();

    sum_ = sum_ptr[0];

    sycl::free((void *)sum_ptr, queue_);

    return true;
  }

  //! Parallel reduction with general reducer
  template <class T, class R, class F, class... Args>
  bool
  parallel_for_reduce(cs_lnum_t n,
                      T&        result,
                      R&        reducer,
                      F&&       f,
                      Args&&... args) {

    // TODO implement this
    return false;
  }

  //! Synchronize associated stream
  template <class... Args>
  bool
  wait(void) {
    if (is_gpu && use_gpu_) {
      queue_.wait();
      return true;
    }
    return false;
  }

  // Get interior faces sum type associated with this context
  template <class M>
  bool
  try_get_parallel_for_i_faces_sum_type(const M                 *m,
                                        cs_dispatch_sum_type_t  &st) {
    if (is_gpu == false || use_gpu_ == false) {
      return false;
    }

    st = CS_DISPATCH_SUM_ATOMIC;
    return true;
  }

  // Get interior faces sum type associated with this context
  template <class M>
  bool
  try_get_parallel_for_b_faces_sum_type(const M                 *m,
                                        cs_dispatch_sum_type_t  &st) {
    if (is_gpu == false || use_gpu_ == false) {
      return false;
    }

    st = CS_DISPATCH_SUM_ATOMIC;
    return true;
  }

};

#elif defined(HAVE_OPENMP_TARGET)

/*!
 * Context to execute loops with OpenMP target on the device
 */

class cs_device_context : public cs_dispatch_context_mixin<cs_device_context> {

private:

  bool              is_gpu;      /*!< Is the associated device à GPU ? */

  bool              use_gpu_;    /*!< Run on GPU ? */

public:

  //! Constructor

  cs_device_context(void)
    : is_gpu(false), use_gpu_(true)
  {
    // This should be improved for any actual use of this approach
    // beyond basic testing
    is_gpu = (omp_get_num_devices() > 1) ? true : false;
  }

  //! Set or unset execution on GPU

  void
  set_use_gpu(bool  use_gpu) {
    this->use_gpu_ = use_gpu;
  }

  //! Check whether we are trying to run on GPU

  bool
  use_gpu(void) {
    return (is_gpu && use_gpu_);
  }

  //! Check preferred allocation mode depending on execution policy

  cs_alloc_mode_t
  alloc_mode(void) {
    cs_alloc_mode_t amode
      = (is_gpu && use_gpu_) ? CS_ALLOC_HOST_DEVICE_SHARED : CS_ALLOC_HOST;
    return (amode);
  }

  cs_alloc_mode_t
  alloc_mode([[maybe_unused]] bool readable_on_cpu) {
    cs_alloc_mode_t amode
      = (is_gpu && use_gpu_) ? CS_ALLOC_HOST_DEVICE_SHARED : CS_ALLOC_HOST;
    return (amode);
  }

public:

  //! Try to launch on the device and return false if not available
  template <class F, class... Args>
  bool
  parallel_for(cs_lnum_t n, F&& f, Args&&... args) {
    if (is_gpu == false || use_gpu_ == false) {
      return false;
    }

    //! Distribute to device
#   pragma omp target teams distribute parallel for
    for (cs_lnum_t i = 0; i < n; ++i) {
      f(i, args...);
    }

    return true;
  }

  //! Launch kernel with simple sum reduction.
  template <class T, class F, class... Args>
  bool
  parallel_for_reduce_sum(cs_lnum_t n,
                          T&        sum,
                          F&&       f,
                          Args&&... args) {
    if (is_gpu == false || use_gpu_ == false) {
      return false;
    }

    sum = 0;
    //! Distribute to device
#   pragma omp target teams distribute parallel for reduction(+:sum)
    for (cs_lnum_t i = 0; i < n; ++i) {
      f(i, sum, args...);
    }

    return true;
  }

  //! Parallel reduction with general reducer.
  template <class T, class R, class F, class... Args>
  bool
  parallel_for_reduce(cs_lnum_t n,
                      T&        result,
                      R&        reducer,
                      F&&       f,
                      Args&&... args) {

    // TODO implement this
    return false;
  }

  //! Synchronize associated stream
  template <class... Args>
  bool
  wait(void) {
    return true;
  }

  // Get interior faces sum type associated with this context
  template <class M>
  bool
  try_get_parallel_for_i_faces_sum_type(const M                 *m,
                                        cs_dispatch_sum_type_t  &st) {
    if (is_gpu == false || use_gpu_ == false) {
      return false;
    }

    st = CS_DISPATCH_SUM_ATOMIC;
    return true;
  }

  // Get interior faces sum type associated with this context
  template <class M>
  bool
  try_get_parallel_for_b_faces_sum_type(const M                 *m,
                                        cs_dispatch_sum_type_t  &st) {
    if (is_gpu == false || use_gpu_ == false) {
      return false;
    }

    st = CS_DISPATCH_SUM_ATOMIC;
    return true;
  }

};

#endif  // __CUDACC__ or SYCL or defined(HAVE_OPENMP_TARGET)

/*!
 * Context to group unused options and catch missing execution paths.
 */

class cs_void_context : public cs_dispatch_context_mixin<cs_void_context> {

public:

  //! Constructor

  cs_void_context(void)
  {}

#if !defined(__CUDACC__)

  /* Fill-in for CUDA methods, so as to allow using these methods
     in final cs_dispatch_context even when CUDA is not available,
     and without requiring a static cast of the form

     static_cast<cs_device_context&>(ctx).set_use_gpu(true);
  */

  void
  set_cuda_grid([[maybe_unused]] long  grid_size,
                [[maybe_unused]] long  block_size) {
  }

  void
  set_cuda_stream([[maybe_unused]] int  stream_id) {
  }

  void
  set_cuda_device([[maybe_unused]] int  device_id) {
  }

#endif  // !defined(__CUDACC__)

#if    !defined(__CUDACC__) \
    && !defined(SYCL_LANGUAGE_VERSION) \
    && !defined(HAVE_OPENMP_TARGET)

  /* Fill-in for device methods */

  void
  set_use_gpu([[maybe_unused]] bool  use_gpu) {
  }

  //! Check whether we are trying to run on GPU

  bool
  use_gpu(void) {
    return false;
  }

  //! Check preferred allocation mode depending on execution policy

  cs_alloc_mode_t
  alloc_mode(void) {
    return CS_ALLOC_HOST;
  }

  cs_alloc_mode_t
  alloc_mode([[maybe_unused]] bool readable_on_cpu) {
    return CS_ALLOC_HOST;
  }

#endif  // ! __CUDACC__ && ! SYCL_LANGUAGE_VERSION && ! defined(HAVE_OPENMP_TARGET)

public:

  // Abort execution if no execution method is available.
  template <class F, class... Args>
  bool parallel_for([[maybe_unused]] cs_lnum_t  n,
                    [[maybe_unused]] F&&        f,
                    [[maybe_unused]] Args&&...  args) {
    cs_assert(0);
    return false;
  }

  // Abort execution if no execution method is available.
  template <class T, class F, class... Args>
  bool parallel_for_reduce_sum([[maybe_unused]] cs_lnum_t  n,
                               [[maybe_unused]] T&         sum,
                               [[maybe_unused]] F&&        f,
                               [[maybe_unused]] Args&&...  args) {
    cs_assert(0);
    return false;
  }

  // Abort execution if no execution method is available.
  template <class T, class R, class F, class... Args>
  bool parallel_for_reduce([[maybe_unused]] cs_lnum_t  n,
                           [[maybe_unused]] T&         result,
                           [[maybe_unused]] R&         reducer,
                           [[maybe_unused]] F&&        f,
                           [[maybe_unused]] Args&&...  args) {
    cs_assert(0);
    return false;
  }

  // Abort execution if no synchronization method is available.
  template <class... Args>
  bool
  wait(void) {
    cs_assert(0);
    return false;
  }

};

/*!
 * cs_context that is a combination of multiple contexts.
 * This context will try every context in order, until one actually runs.
 */

template <class... Contexts>
class cs_combined_context
  : public cs_dispatch_context_mixin<cs_combined_context<Contexts...>>,
    public Contexts... {

private:
  using mixin_t = cs_dispatch_context_mixin<cs_combined_context<Contexts...>>;

public:
  cs_combined_context() = default;
  cs_combined_context(Contexts... contexts)
    : Contexts(std::move(contexts))...
  {}

public:

  /*--------------------------------------------------------------------------*/
  /* \brief Parallel computation over interior faces.
   *
   * This method is intended for use when assembling cell values
   * with face-based computations, using the appropriate \ref cs_dispatch_sum
   * functions.
   *
   * On CPU, loops are scheduled based on the current face numbering, so
   * as to avoid thread races when summing values. On GPU, atomic sums are used.
   *
   * \tparam M  mesh type structure (templated mostly to avoid
   *            dependency to mesh definitions in lower level code)
   * \tparam F  lambda function or functor
   *
   * \param[in]  m     pointer to  mesh
   * \param[in]  f     lambda function or functor to execute
   */
  /*--------------------------------------------------------------------------*/

  template <class M, class F, class... Args>
  auto parallel_for_i_faces(const M* m, F&& f, Args&&... args) {
    bool launched = false;
    [[maybe_unused]] decltype(nullptr) try_execute[] = {
      (   launched = launched
       || Contexts::parallel_for_i_faces(m, f, args...), nullptr)...
    };
  }

  /*--------------------------------------------------------------------------*/
  /* \brief Parallel computation over boundary faces.
   *
   * This method is intended for use when assembling cell values
   * with face-based computations, using the appropriate \ref cs_dispatch_sum
   * functions.
   *
   * On CPU, loops are scheduled based on the current face numbering, so
   * as to avoid thread races when summing values. On GPU, atomic sums are used.
   *
   * \tparam M  mesh type structure (templated mostly to avoid
   *            dependency to mesh definitions in lower level code)
   * \tparam F  lambda function or functor
   *
   * \param[in]  m     pointer to  mesh
   * \param[in]  f     lambda function or functor to execute
   */
  /*--------------------------------------------------------------------------*/

  template <class M, class F, class... Args>
  auto parallel_for_b_faces(const M* m, F&& f, Args&&... args) {
    bool launched = false;
    [[maybe_unused]] decltype(nullptr) try_execute[] = {
      (   launched = launched
       || Contexts::parallel_for_b_faces(m, f, args...), nullptr)...
    };
  }

  /*--------------------------------------------------------------------------*/
  /* \brief General parallel computation over elements.
   *
   * \tparam F  lambda function or functor
   *
   * \param[in]  n  number of elements to compute
   * \param[in]  f  lambda function or functor to execute
   */
  /*--------------------------------------------------------------------------*/

  template <class F, class... Args>
  auto parallel_for(cs_lnum_t n, F&& f, Args&&... args) {
    bool launched = false;
    [[maybe_unused]] decltype(nullptr) try_execute[] = {
      (   launched = launched
       || Contexts::parallel_for(n, f, args...), nullptr)...
    };
  }

  /*--------------------------------------------------------------------------*/
  /* \brief General parallel computation over elements, with a floating-point
   *        sum reduction.
   *
   * \tparam T  reduced element type
   * \tparam F  lambda function or functor
   *
   * \param[in]  n    number of elements to compute
   * \param[in]  sum  resulting sum
   * \param[in]  f    lambda function or functor to execute
   */
  /*--------------------------------------------------------------------------*/

  template <class T, class F, class... Args>
  auto parallel_for_reduce_sum
    (cs_lnum_t n, T& sum, F&& f, Args&&... args) {
    bool launched = false;
    [[maybe_unused]] decltype(nullptr) try_execute[] = {
      (   launched = launched
       || Contexts::parallel_for_reduce_sum(n, sum, f, args...),
          nullptr)...
    };
  }

  /*--------------------------------------------------------------------------*/
  /* \brief General parallel computation over elements, with a
   *        user-defined reduction.
   *
   * \tparam T  reduced element type
   * \tparam R  reducer class
   * \tparam F  lambda function or functor
   *
   * \param[in]   n        number of elements to compute
   * \param[out]  result   resulting sum
   * \param[in]   reducer  reducer object
   * \param[in]   f        lambda function or functor to execute
   */
  /*--------------------------------------------------------------------------*/

  template <class T, class R, class F, class... Args>
  auto parallel_for_reduce
    (cs_lnum_t n, T& result, R& reducer, F&& f, Args&&... args) {
    bool launched = false;
    [[maybe_unused]] decltype(nullptr) try_execute[] = {
      (   launched = launched
       || Contexts::parallel_for_reduce(n, result, reducer, f, args...),
          nullptr)...
    };
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Wait (synchronize) until launched computations have finished.
   */
  /*--------------------------------------------------------------------------*/

  void
    wait(void) {
    bool done = false;
    [[maybe_unused]] decltype(nullptr) try_execute[] = {
      (   done = done
       || Contexts::wait(), nullptr)...
    };
  }

  /*--------------------------------------------------------------------------*/
  /*! \brief Return sum type to be used with \ref parallel_for_i_faces.
   *
   * \tparam M  mesh type structure (templated mostly to avoid
   *            dependency to mesh definitions in lower level code)
   *
   * \param[in]  m     pointer to  mesh
   *
   * \return  assembly sum type that should be used in parallel_for_i_faces
   */
  /*--------------------------------------------------------------------------*/

  template <class M>
  cs_dispatch_sum_type_t
  get_parallel_for_i_faces_sum_type(const M* m) {
    cs_dispatch_sum_type_t sum_type = CS_DISPATCH_SUM_ATOMIC;
    bool known = false;
    [[maybe_unused]] decltype(nullptr) try_query[] = {
      (   known = known
       || Contexts::try_get_parallel_for_i_faces_sum_type(m, sum_type),
          nullptr)...
    };
    return sum_type;
  }

  /*--------------------------------------------------------------------------*/
  /*! \brief Return sum type to be used with \ref parallel_for_b_faces.
   *
   * \tparam M  mesh type structure (templated mostly to avoid
   *            dependency to mesh definitions in lower level code)
   *
   * \param[in]  m     pointer to  mesh
   *
   * \return  assembly sum type that should be used in parallel_for_b_faces
   */
  /*--------------------------------------------------------------------------*/

  template <class M>
  cs_dispatch_sum_type_t
  get_parallel_for_b_faces_sum_type(const M* m) {
    cs_dispatch_sum_type_t sum_type = CS_DISPATCH_SUM_ATOMIC;
    bool known = false;
    [[maybe_unused]] decltype(nullptr) try_query[] = {
      (   known = known
       || Contexts::try_get_parallel_for_b_faces_sum_type(m, sum_type),
          nullptr)...
    };
    return sum_type;
  }

};

/*----------------------------------------------------------------------------*/
/*!
 * Default cs_dispatch_context that is a combination of (CPU) and GPU (CUDA)
 * context if available.
 */
/*----------------------------------------------------------------------------*/

class cs_dispatch_context : public cs_combined_context<
#if   defined(__CUDACC__) \
  || defined(SYCL_LANGUAGE_VERSION) \
  || defined(HAVE_OPENMP_TARGET)
  cs_device_context,
#endif
  cs_host_context,
  cs_void_context
>
{

private:
  using base_t = cs_combined_context<
#if   defined(__CUDACC__) \
   || defined(SYCL_LANGUAGE_VERSION) \
   || defined(HAVE_OPENMP_TARGET)
  cs_device_context,
#endif
  cs_host_context,
  cs_void_context
>;

public:
  using base_t::base_t;
  using base_t::operator=;

};

/*
  Remarks:

  Instantiation can simply be done using:

  `cs_dispatch_context ctx;`

  Instanciation can also be done with specific construction options,
  for example:

  `cs_dispatch_context ctx(cs_device_context(stream), {});`

  or:

  `cs_dispatch_context ctx(cs_device_context(), {});`

*/

/*=============================================================================
 * Global variable definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  sum values using a chosen dispatch sum type.
 *
 * This allows calling a unique function for assembly in lambda functions
 * called by a dispatch context for interior or boundary faces.
 *
 * \tparam  T    array type
 *
 * \param[in, out]  dest      destination
 * \param[in]       src       source value
 * \param[in]       sum_type  sum type
 */
/*----------------------------------------------------------------------------*/

#ifdef __CUDA_ARCH__   // Test whether we are on GPU or CPU...

template <typename T>
__device__ static void __forceinline__
cs_dispatch_sum(T                       *dest,
                const T                  src,
                cs_dispatch_sum_type_t   sum_type)
{
  if (sum_type == CS_DISPATCH_SUM_ATOMIC) {
#if 1
    using sum_v = assembled_value<T>;
    sum_v v;

    v.get() = src;
    sum_v::ref(*dest).conflict_free_add(-1u, v);
#else
    atomicAdd(dest, src);
#endif
  }
  else if (sum_type == CS_DISPATCH_SUM_SIMPLE) {
    *dest += src;
  }
}

#elif defined(SYCL_LANGUAGE_VERSION)

template <typename T>
inline void
cs_dispatch_sum(T                       *dest,
                const T                  src,
                cs_dispatch_sum_type_t   sum_type)
{
  if (sum_type == CS_DISPATCH_SUM_SIMPLE) {
    *dest += src;
  }
  else if (sum_type == CS_DISPATCH_SUM_ATOMIC) {
    sycl::atomic_ref<T,
                     sycl::memory_order::relaxed,
                     sycl::memory_scope::device> aref(*dest);
    aref.fetch_add(src);
  }
}

#else  // ! CUDA or SYCL

template <typename T>
inline void
cs_dispatch_sum(T                       *dest,
                const T                  src,
                cs_dispatch_sum_type_t   sum_type)
{
  if (sum_type == CS_DISPATCH_SUM_SIMPLE) {
    *dest += src;
  }
  else if (sum_type == CS_DISPATCH_SUM_ATOMIC) {
    #pragma omp atomic
    *dest += src;
  }
}

#endif // __CUDA_ARCH__

/*----------------------------------------------------------------------------*/
/*!
 * \brief  sum values using a chosen dispatch sum type.
 *
 * This allows calling a unique function for assembly in lambda functions
 * called by a dispatch context for interior or boundary faces.
 *
 * \tparam  dim  array stride
 * \tparam  T    array type
 *
 * \param[in, out]  dest      destination values
 * \param[in]       src       source values
 * \param[in]       sum_type  sum type
 */
/*----------------------------------------------------------------------------*/

#ifdef __CUDA_ARCH__   // Test whether we are on GPU or CPU...

template <size_t dim, typename T>
__device__ static void __forceinline__
cs_dispatch_sum(T                       *dest,
                const T                 *src,
                cs_dispatch_sum_type_t   sum_type)
{
  if (sum_type == CS_DISPATCH_SUM_SIMPLE) {
    for (cs_lnum_t i = 0; i < dim; i++) {
      dest[i] += src[i];
    }
  }
  else if (sum_type == CS_DISPATCH_SUM_ATOMIC) {
#if __CUDA_ARCH__ >= 700
    using sum_v = assembled_value<T, dim>;
    sum_v v;

    for (size_t i = 0; i < dim; i++) {
      v[i].get() = src[i];
    }

    sum_v &vs = reinterpret_cast<sum_v &>(*dest);
    vs.conflict_free_add(-1u, v);

    //sum_v::ref(dest).conflict_free_add(-1u, v);
#else
    for (size_t i = 0; i < dim; i++) {
      atomicAdd(&dest[i], src[i]);
    }
#endif
  }
}

#elif defined(SYCL_LANGUAGE_VERSION)

template <size_t dim, typename T>
inline void
cs_dispatch_sum(T                       *dest,
                const T                 *src,
                cs_dispatch_sum_type_t   sum_type)
{
  if (sum_type == CS_DISPATCH_SUM_SIMPLE) {
    for (size_t i = 0; i < dim; i++) {
      dest[i] += src[i];
    }
  }
  else if (sum_type == CS_DISPATCH_SUM_ATOMIC) {
    for (size_t i = 0; i < dim; i++) {
      sycl::atomic_ref<T,
                       sycl::memory_order::relaxed,
                       sycl::memory_scope::device> aref(dest[i]);
      aref.fetch_add(src[i]);
    }
  }
}

#else  // ! CUDA or SYCL

template <size_t dim, typename T>
inline void
cs_dispatch_sum(T                       *dest,
                const T                 *src,
                cs_dispatch_sum_type_t   sum_type)
{
  if (sum_type == CS_DISPATCH_SUM_SIMPLE) {
    for (size_t i = 0; i < dim; i++) {
      dest[i] += src[i];
    }
  }
  else if (sum_type == CS_DISPATCH_SUM_ATOMIC) {
    for (size_t i = 0; i < dim; i++) {
      #pragma omp atomic
      dest[i] += src[i];
    }
  }
}

#endif // __CUDA_ARCH__

#endif /* __cplusplus */

/*----------------------------------------------------------------------------*/

#endif /* CS_DISPATCH_H */
