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
#include "cs_assert.h"
#include "cs_mesh.h"

#ifdef __NVCC__
#include "cs_base_cuda.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C++ library headers
 *----------------------------------------------------------------------------*/

#include <utility>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __NVCC__

#define CS_F_HOST __host__
#define CS_F_DEVICE __device__
#define CS_F_HOST_DEVICE __host__ __device__

#else

#define CS_F_HOST
#define CS_F_DEVICE
#define CS_F_HOST_DEVICE

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

} cs_dispatch_sum_t;

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
  decltype(auto) parallel_for(cs_lnum_t n, F&& f, Args&&... args) = delete;

  // Loop over n elements with sum type info
  // Must be redefined by the child class
  template <class F, class... Args>
  decltype(auto) parallel_for_sum(cs_lnum_t n, F&& f, Args&&... args) = delete;

  // Assembly loop over all internal faces
  template <class F, class... Args>
  decltype(auto) parallel_for_i_faces(const cs_mesh_t*  m,
                                      F&&               f,
                                      Args&&...         args);

  // Assembly loop over all boundary faces
  template <class F, class... Args>
  decltype(auto) parallel_for_b_faces(const cs_mesh_t*  m,
                                      F&&               f,
                                      Args&&...         args);

};

// Default implementation of parallel_for_i_faces based on parallel_for_sum
template <class Derived>
template <class F, class... Args>
decltype(auto) cs_dispatch_context_mixin<Derived>::parallel_for_i_faces
  (const cs_mesh_t* m, F&& f, Args&&... args) {
  return static_cast<Derived*>(this)->parallel_for_sum
                                        (m->n_i_faces,
                                         static_cast<F&&>(f),
                                         static_cast<Args&&>(args)...);
}

// Default implementation of parallel_for_b_faces based on parallel_for_sum
template <class Derived>
template <class F, class... Args>
decltype(auto) cs_dispatch_context_mixin<Derived>::parallel_for_b_faces
  (const cs_mesh_t* m, F&& f, Args&&... args) {
  return static_cast<Derived*>(this)->parallel_for_sum
                                        (m->n_b_faces,
                                         static_cast<F&&>(f),
                                         static_cast<Args&&>(args)...);
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

  /*! Set minimum number of elements threshold for CPU multithread execution */
  void set_n_min_for_cpu_threads(cs_lnum_t  n) {
    this->n_min_for_threads = n;
  }

  // Iterate using a plain omp parallel for
  template <class F, class... Args>
  bool
  parallel_for(cs_lnum_t n, F&& f, Args&&... args) {
#   pragma omp parallel for  if (n >= n_min_for_threads)
    for (cs_lnum_t i = 0; i < n; ++i) {
      f(i, args...);
    }
    return true;
  }

  // Iterate using a plain omp parallel for
  template <class F, class... Args>
  bool
  parallel_for_sum(cs_lnum_t n, F&& f, Args&&... args) {
#   pragma omp parallel for  if (n >= n_min_for_threads)
    for (cs_lnum_t i = 0; i < n; ++i) {
      f(i, CS_DISPATCH_SUM_ATOMIC, args...);
    }
    return true;
  }

  // Loop over the interior faces of a mesh using a specific numbering
  // that avoids conflicts between threads.
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
          f(f_id, CS_DISPATCH_SUM_SIMPLE, args...);
        }
      }
    }
    return true;
  }

  // Loop over the boundary faces of a mesh using a specific numbering
  // that avoids conflicts between threads.
  template <class F, class... Args>
  bool parallel_for_b_faces(const cs_mesh_t* m, F&& f, Args&&... args) {
    const int n_b_threads = m->b_face_numbering->n_threads;
    const cs_lnum_t *restrict b_group_index = m->b_face_numbering->group_index;

    #pragma omp parallel for
    for (int t_id = 0; t_id < n_b_threads; t_id++) {
      for (cs_lnum_t f_id = b_group_index[t_id*2];
           f_id < b_group_index[t_id*2 + 1];
           f_id++) {
        f(f_id, CS_DISPATCH_SUM_SIMPLE, args...);
      }
    }
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

/* Kernel for assembly that loops over an integer range and calls a device
   functor. This is similar to the above, but passes an extra argument
   indicating the sum type to use, allowing a templated test.
   This kernel uses a grid_size-stride loop and thus guarantees that all
   integers are processed, even if the grid is smaller.
   All arguments *must* be passed by value to avoid passing CPU references
   to the GPU. */

template <class F, class... Args>
__global__ void cs_cuda_kernel_parallel_for_sum(cs_lnum_t n, F f, Args... args) {
  // grid_size-stride loop
  for (cs_lnum_t id = blockIdx.x * blockDim.x + threadIdx.x; id < n;
       id += blockDim.x * gridDim.x) {
    f(id, CS_DISPATCH_SUM_ATOMIC, args...);
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

  cs_lnum_t  n_min_for_device_;  /*!< Run on CPU under this threshold */

public:

  //! Constructor

  cs_device_context(void)
    : grid_size_(0), block_size_(256), stream_(cs_cuda_get_stream(0)),
      device_(0), n_min_for_device_(0)
  {
    device_ = cs_base_cuda_get_device();
  }

  cs_device_context(long          grid_size,
                    long          block_size,
                    cudaStream_t  stream,
                    int           device)
    : grid_size_(grid_size), block_size_(block_size), stream_(stream),
      device_(device), n_min_for_device_(0)
  {}

  cs_device_context(long          grid_size,
                    long          block_size,
                    cudaStream_t  stream)
    : grid_size_(grid_size), block_size_(block_size), stream_(stream),
      device_(0), n_min_for_device_(0)
  {
    device_ = cs_base_cuda_get_device();
  }

  cs_device_context(long  grid_size,
                    long  block_size)
    : grid_size_(grid_size), block_size_(block_size),
       stream_(cs_cuda_get_stream(0)), device_(0), n_min_for_device_(0)
  {
    device_ = cs_base_cuda_get_device();
  }

  cs_device_context(cudaStream_t  stream)
    : grid_size_(0), block_size_(256), stream_(stream), device_(0),
      n_min_for_device_(0)
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

  //! Change CUDA device

  void
  set_cuda_device(int  device) {
    this->device_ = device;
  }

  //! Set minimum number of elements threshold for GPU execution

  void
  set_n_min_for_gpu(cs_lnum_t  n) {
    this->n_min_for_device_ = n;
  }

public:

  //! Try to launch on the GPU and return false if not available
  template <class F, class... Args>
  bool
  parallel_for(cs_lnum_t n, F&& f, Args&&... args) {
    if (device_ < 0 || n < n_min_for_device_) {
      return false;
    }

    long l_grid_size = grid_size_;
    if (l_grid_size < 1) {
      l_grid_size = (n % block_size_) ? n/block_size_ + 1 : n/block_size_;
    }

    cs_cuda_kernel_parallel_for<<<l_grid_size, block_size_, 0, stream_>>>
      (n, static_cast<F&&>(f), static_cast<Args&&>(args)...);

    return true;
  }

  //! Try to launch on the GPU and return false if not available
  template <class F, class... Args>
  bool
  parallel_for_i_faces(const cs_mesh_t* m, F&& f, Args&&... args) {
    const cs_lnum_t n = m->n_i_faces;
    if (device_ < 0 || n < n_min_for_device_) {
      return false;
    }

    long l_grid_size = grid_size_;
    if (l_grid_size < 1) {
      l_grid_size = (n % block_size_) ? n/block_size_ + 1 : n/block_size_;
    }

    cs_cuda_kernel_parallel_for_sum<<<l_grid_size, block_size_, 0, stream_>>>
      (n, static_cast<F&&>(f), static_cast<Args&&>(args)...);

    return true;
  }

  //! Synchronize associated stream
  void
  wait(void) {
    cudaStreamSynchronize(stream_);
  }

};

#elif defined(SYCL_LANGUAGE_VERSION)

/*!
 * Context to execute loops with SYCL on the device
 */


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

     static_cast<cs_device_context&>(ctx).set_n_min_for_gpu(200);
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
  set_n_min_for_gpu([[maybe_unused]] int  n) {
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
  bool parallel_for_sum([[maybe_unused]] cs_lnum_t  n,
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
    void parallel_for_sum(const cs_mesh_t* m, F&& f, Args&&... args) {
    bool launched = false;
    [[maybe_unused]] decltype(nullptr) try_execute[] = {
      (   launched = launched
       || Contexts::parallel_for_sum(m, f, args...), nullptr)...
    };
  }

  template <class F, class... Args>
    void parallel_for_i_faces(const cs_mesh_t* m, F&& f, Args&&... args) {
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
};

/*----------------------------------------------------------------------------*/
/*!
 * Default cs_dispatch_context that is a combination of (CPU) and GPU (CUDA)
 * context if available.
 */
/*----------------------------------------------------------------------------*/

class cs_dispatch_context : public cs_combined_context<
#ifdef __NVCC__
  cs_device_context,
#endif
  cs_host_context,
  cs_void_context
>
{

private:
  using base_t = cs_combined_context<
#ifdef __NVCC__
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

#endif /* __cplusplus */

#endif /* __CS_DISPATCH_H__ */
