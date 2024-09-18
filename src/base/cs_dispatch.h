#ifndef __CS_DISPATCH_H__
#define __CS_DISPATCH_H__

/*============================================================================
 * Class to dispatch computation using various runtimes (OpenMP, CUDA, ...)
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

// Valid only for C++

#ifdef __cplusplus

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C++ library headers
 *----------------------------------------------------------------------------*/

#include <utility>

#if defined(SYCL_LANGUAGE_VERSION)
#include <sycl/sycl.hpp>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_assert.h"
#include "cs_mesh.h"

#ifdef __NVCC__
#include "cs_base_cuda.h"
#include "cs_blas_cuda.h"
#include "cs_alge_cuda.cuh"
#endif

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#if defined(SYCL_LANGUAGE_VERSION)

#define CS_DISPATCH_SUM_DOUBLE auto

#else

#define CS_DISPATCH_SUM_DOUBLE double

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
  template <class F, class... Args>
  decltype(auto)
  parallel_for_i_faces(const cs_mesh_t*  m,
                       F&&               f,
                       Args&&...         args);

  // Assembly loop over all boundary faces
  template <class F, class... Args>
  decltype(auto)
  parallel_for_b_faces(const cs_mesh_t*  m,
                       F&&               f,
                       Args&&...         args);

  // Parallel reduction with simple sum.
  // Must be redefined by the child class
  template <class F, class... Args>
  decltype(auto)
  parallel_for_reduce_sum
    (cs_lnum_t n, double& sum, F&& f, Args&&... args) = delete;

  // Query sum type for assembly loop over all interior faces
  // Must be redefined by the child class
  bool
  try_get_parallel_for_i_faces_sum_type(const cs_mesh_t*         m,
                                        cs_dispatch_sum_type_t&  st);

  // Query sum type for assembly loop over all boundary faces
  // Must be redefined by the child class
  bool
  try_get_parallel_for_b_faces_sum_type(const cs_mesh_t*         m,
                                        cs_dispatch_sum_type_t&  st);

};

// Default implementation of parallel_for_i_faces based on parallel_for
template <class Derived>
template <class F, class... Args>
decltype(auto) cs_dispatch_context_mixin<Derived>::parallel_for_i_faces
  (const cs_mesh_t* m, F&& f, Args&&... args) {
  return static_cast<Derived*>(this)->parallel_for
                                        (m->n_i_faces,
                                         static_cast<F&&>(f),
                                         static_cast<Args&&>(args)...);
}

// Default implementation of parallel_for_b_faces based on parallel_for_sum
template <class Derived>
template <class F, class... Args>
decltype(auto) cs_dispatch_context_mixin<Derived>::parallel_for_b_faces
  (const cs_mesh_t* m, F&& f, Args&&... args) {
  return static_cast<Derived*>(this)->parallel_for
                                        (m->n_b_faces,
                                         static_cast<F&&>(f),
                                         static_cast<Args&&>(args)...);
}

// Default implementation of get interior faces sum type
template <class Derived>
bool cs_dispatch_context_mixin<Derived>::try_get_parallel_for_i_faces_sum_type
  ([[maybe_unused]]const cs_mesh_t*  m,
   cs_dispatch_sum_type_t&           st) {
  st = CS_DISPATCH_SUM_SIMPLE;
  return true;
}

// Default implementation of get boundary faces sum type
template <class Derived>
bool cs_dispatch_context_mixin<Derived>::try_get_parallel_for_b_faces_sum_type
  ([[maybe_unused]]const cs_mesh_t*  m,
   cs_dispatch_sum_type_t&           st) {
  st = CS_DISPATCH_SUM_SIMPLE;
  return true;
}

/*
 * cs_context to execute loops with OpenMP on the CPU
 */

class cs_host_context : public cs_dispatch_context_mixin<cs_host_context> {

private:

  cs_lnum_t  n_min_for_threads;  /*!< Run on single thread
                                   under this threshold */

public:

  cs_host_context()
    : n_min_for_threads(CS_THR_MIN)
  {}

public:

  //! Set minimum number of elements threshold for CPU multithread execution.
  void
  set_n_min_for_cpu_threads(cs_lnum_t  n) {
    this->n_min_for_threads = n;
  }

  //! Get minimum number of elements threshold for CPU multithread execution.
  cs_lnum_t
  n_min_for_cpu_threads(void) {
    return this->n_min_for_threads;
  }

  //! Iterate using a plain omp parallel for
  template <class F, class... Args>
  bool
  parallel_for(cs_lnum_t n, F&& f, Args&&... args) {
#   pragma omp parallel for  if (n >= n_min_for_threads)
    for (cs_lnum_t i = 0; i < n; ++i) {
      f(i, args...);
    }
    return true;
  }

  //! Loop over the interior faces of a mesh using a specific numbering
  //! that avoids conflicts between threads.
  template <class F, class... Args>
  bool
  parallel_for_i_faces(const cs_mesh_t* m, F&& f, Args&&... args) {
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
  template <class F, class... Args>
  bool
  parallel_for_b_faces(const cs_mesh_t* m, F&& f, Args&&... args) {
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
  template <class F, class... Args>
  bool
  parallel_for_reduce_sum(cs_lnum_t n,
                          double&   sum,
                          F&&       f,
                          Args&&... args) {
    sum = 0;
#   pragma omp parallel for reduction(+:sum) if (n >= n_min_for_threads)
    for (cs_lnum_t i = 0; i < n; ++i) {
      f(i, sum, args...);
    }
    return true;
  }

  // Get interior faces sum type associated with this context
  bool
  try_get_parallel_for_i_faces_sum_type([[maybe_unused]]const cs_mesh_t*  m,
                                        cs_dispatch_sum_type_t&           st) {
    st = CS_DISPATCH_SUM_SIMPLE;
    return true;
  }

  // Get boundary faces sum type associated with this context
  bool
  try_get_parallel_for_b_faces_sum_type([[maybe_unused]]const cs_mesh_t*  m,
                                        cs_dispatch_sum_type_t&           st) {
    st = CS_DISPATCH_SUM_SIMPLE;
    return true;
  }

};

#if defined(__NVCC__)

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

/* Default kernel that loops over an integer range and calls a device functor.
   This kernel uses a grid_size-stride loop and thus guarantees that all
   integers are processed, even if the grid is smaller.
   All arguments *must* be passed by value to avoid passing CPU references
   to the GPU. */

template <class F, class... Args>
__global__ void
cs_cuda_kernel_parallel_for_reduce_sum(cs_lnum_t   n,
                                       double     *b_res,
                                       F           f,
                                       Args...     args) {
  // grid_size-stride loop
  extern double __shared__ stmp[];
  const cs_lnum_t tid = threadIdx.x;

  stmp[tid] = 0;

  for (cs_lnum_t id = blockIdx.x * blockDim.x + threadIdx.x; id < n;
       id += blockDim.x * gridDim.x) {
     f(id, stmp[tid], args...);
  }

  switch (blockDim.x) {
  case 1024:
    cs_blas_cuda_block_reduce_sum<1024, 1>(stmp, tid, b_res);
    break;
  case 512:
    cs_blas_cuda_block_reduce_sum<512, 1>(stmp, tid, b_res);
    break;
  case 256:
    cs_blas_cuda_block_reduce_sum<256, 1>(stmp, tid, b_res);
    break;
  case 128:
    cs_blas_cuda_block_reduce_sum<128, 1>(stmp, tid, b_res);
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
    device_ = cs_base_cuda_get_device();
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

  void
  set_cuda_grid(long  grid_size,
                long  block_size) {
    this->grid_size_ = grid_size;
    this->block_size_ = block_size;
  }

  //! Change stream, but keeps the grid and device configuration

  void
  set_cuda_stream(cudaStream_t stream) {
    this->stream_ = stream;
  }

  //! Change stream, but keeps the grid and device configuration

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
  template <class F, class... Args>
  bool
  parallel_for_i_faces(const cs_mesh_t* m, F&& f, Args&&... args) {
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
  template <class F, class... Args>
  bool
  parallel_for_reduce_sum(cs_lnum_t n,
                          double&   sum,
                          F&&       f,
                          Args&&... args) {
    sum = 0;
    if (device_ < 0 || use_gpu_ == false) {
      return false;
    }

    long l_grid_size = grid_size_;
    if (l_grid_size < 1) {
      l_grid_size = (n % block_size_) ? n/block_size_ + 1 : n/block_size_;
    }
    if (n == 0) {
      sum = 0;
      return true;
    }

    double *r_grid_, *r_reduce_;
    cs_blas_cuda_get_2_stage_reduce_buffers
      (n, 1, l_grid_size, r_grid_, r_reduce_);

    int smem_size = block_size_ * sizeof(double);
    cs_cuda_kernel_parallel_for_reduce_sum
      <<<l_grid_size, block_size_, smem_size, stream_>>>
      (n, r_grid_, static_cast<F&&>(f), static_cast<Args&&>(args)...);

    switch (block_size_) {
    case 1024:
      cs_blas_cuda_reduce_single_block<1024, 1>
        <<<1, block_size_, 0, stream_>>>
        (l_grid_size, r_grid_, r_reduce_);
      break;
    case 512:
      cs_blas_cuda_reduce_single_block<512, 1>
        <<<1, block_size_, 0, stream_>>>
        (l_grid_size, r_grid_, r_reduce_);
      break;
    case 256:
      cs_blas_cuda_reduce_single_block<256, 1>
        <<<1, block_size_, 0, stream_>>>
        (l_grid_size, r_grid_, r_reduce_);
      break;
    case 128:
      cs_blas_cuda_reduce_single_block<128, 1>
        <<<1, block_size_, 0, stream_>>>
        (l_grid_size, r_grid_, r_reduce_);
      break;
    default:
      cs_assert(0);
    }

    cudaStreamSynchronize(stream_);
    sum = r_reduce_[0];

    return true;
  }

  //! Synchronize associated stream
  void
  wait(void) {
    if (device_ > -1 && use_gpu_)
      cudaStreamSynchronize(stream_);
  }

  // Get interior faces sum type associated with this context
  bool
  try_get_parallel_for_i_faces_sum_type(const cs_mesh_t         *m,
                                        cs_dispatch_sum_type_t  &st) {
    if (device_ < 0 || use_gpu_ == false) {
      return false;
    }

    st = CS_DISPATCH_SUM_ATOMIC;
    return true;
  }

  // Get boundary faces sum type associated with this context
  bool
  try_get_parallel_for_b_faces_sum_type(const cs_mesh_t         *m,
                                        cs_dispatch_sum_type_t  &st) {
    if (device_ < 0 || use_gpu_ == false) {
      return false;
    }

    st = CS_DISPATCH_SUM_ATOMIC;
    return true;
  }

};

#elif defined(SYCL_LANGUAGE_VERSION)

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
  template <class F, class... Args>
  bool
  parallel_for_i_faces(const cs_mesh_t* m, F&& f, Args&&... args) {
    const cs_lnum_t n = m->n_i_faces;
    if (is_gpu == false || use_gpu_ == false) {
      return false;
    }

    queue_.parallel_for(n, static_cast<F&&>(f), static_cast<Args&&>(args)...);

    return true;
  }

  //! Launch kernel with simple sum reduction.
  template <class F, class... Args>
  bool
  parallel_for_reduce_sum(cs_lnum_t n,
                          double&   sum_,
                          F&&       f,
                          Args&&... args) {
    sum_ = 0;
    if (is_gpu == false || use_gpu_ == false) {
      return false;
    }

    // TODO: use persistent allocation as we do in CUDA BLAS to avoid
    //       excess allocation/deallocation.
    double *sum_ptr = (double *)sycl::malloc_shared(sizeof(double), queue_);

    queue_.parallel_for(n,
                        sycl::reduction(sum_ptr, 0., sycl::plus<double>()),
                        static_cast<F&&>(f),
                        static_cast<Args&&>(args)...).wait();

    sum_ = sum_ptr[0];

    sycl::free((void *)sum_ptr, queue_);

    return true;
  }

  //! Synchronize associated stream
  void
  wait(void) {
    if (is_gpu && use_gpu_)
      queue_.wait();
  }

  // Get interior faces sum type associated with this context
  bool
  try_get_parallel_for_i_faces_sum_type(const cs_mesh_t         *m,
                                        cs_dispatch_sum_type_t  &st) {
    if (is_gpu == false || use_gpu_ == false) {
      return false;
    }

    st = CS_DISPATCH_SUM_ATOMIC;
    return true;
  }

  // Get interior faces sum type associated with this context
  bool
  try_get_parallel_for_b_faces_sum_type(const cs_mesh_t         *m,
                                        cs_dispatch_sum_type_t  &st) {
    if (is_gpu == false || use_gpu_ == false) {
      return false;
    }

    st = CS_DISPATCH_SUM_ATOMIC;
    return true;
  }

};

#endif  // __NVCC__ or SYCL

/*!
 * Context to group unused options and catch missing execution paths.
 */

class cs_void_context : public cs_dispatch_context_mixin<cs_void_context> {

public:

  //! Constructor

  cs_void_context(void)
  {}

#if !defined(__NVCC__)

  /* Fill-in for CUDA methods, so as to allow using these methods
     in final cs_dispatch_context even when CUDA is not available,
     and without requireing a static cast of the form

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

#endif  // __NVCC__

#if !defined(__NVCC__) && !defined(SYCL_LANGUAGE_VERSION)

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

  void
  wait(void) {
  }

#endif  // ! __NVCC__ && ! SYCL_LANGUAGE_VERSION

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
  template <class F, class... Args>
  bool parallel_for_reduce_sum([[maybe_unused]] cs_lnum_t  n,
                               [[maybe_unused]] double&    sum,
                               [[maybe_unused]] F&&        f,
                               [[maybe_unused]] Args&&...  args) {
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

  template <class F, class... Args>
  auto parallel_for_i_faces(const cs_mesh_t* m, F&& f, Args&&... args) {
    bool launched = false;
    [[maybe_unused]] decltype(nullptr) try_execute[] = {
      (   launched = launched
       || Contexts::parallel_for_i_faces(m, f, args...), nullptr)...
    };
  }

  template <class F, class... Args>
  auto parallel_for_b_faces(const cs_mesh_t* m, F&& f, Args&&... args) {
    bool launched = false;
    [[maybe_unused]] decltype(nullptr) try_execute[] = {
      (   launched = launched
       || Contexts::parallel_for_b_faces(m, f, args...), nullptr)...
    };
  }

  template <class F, class... Args>
  auto parallel_for(cs_lnum_t n, F&& f, Args&&... args) {
    bool launched = false;
    [[maybe_unused]] decltype(nullptr) try_execute[] = {
      (   launched = launched
       || Contexts::parallel_for(n, f, args...), nullptr)...
    };
  }

  template <class F, class... Args>
  auto parallel_for_reduce_sum
    (cs_lnum_t n, double& sum, F&& f, Args&&... args) {
    bool launched = false;
    [[maybe_unused]] decltype(nullptr) try_execute[] = {
      (   launched = launched
       || Contexts::parallel_for_reduce_sum(n, sum, f, args...),
          nullptr)...
    };
  }

  cs_dispatch_sum_type_t
  get_parallel_for_i_faces_sum_type(const cs_mesh_t* m) {
    cs_dispatch_sum_type_t sum_type = CS_DISPATCH_SUM_ATOMIC;
    bool known = false;
    [[maybe_unused]] decltype(nullptr) try_query[] = {
      (   known = known
       || Contexts::try_get_parallel_for_i_faces_sum_type(m, sum_type),
          nullptr)...
    };
    return sum_type;
  }

  cs_dispatch_sum_type_t
  get_parallel_for_b_faces_sum_type(const cs_mesh_t* m) {
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
#if defined(__NVCC__) || defined(SYCL_LANGUAGE_VERSION)
  cs_device_context,
#endif
  cs_host_context,
  cs_void_context
>
{

private:
  using base_t = cs_combined_context<
#if defined(__NVCC__) || defined(SYCL_LANGUAGE_VERSION)
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

  Instanciation can also be done with specific contruction options,
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

#endif /* __CS_DISPATCH_H__ */
