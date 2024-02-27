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

/*----------------------------------------------------------------------------
 * Macro to remove parentheses
 *
 * CS_REMOVE_PARENTHESES(x)   -> x
 * CS_REMOVE_PARENTHESES((x)) -> x
 *----------------------------------------------------------------------------*/

#define CS_REMOVE_PARENTHESES(x) CS_REMOVE_PARENTHESES_ESCAPE(CS_PARENTHESIS x)
#define CS_REMOVE_PARENTHESES_ESCAPE_(...) CS_VANISH_ ## __VA_ARGS__
#define CS_REMOVE_PARENTHESES_ESCAPE(...) CS_REMOVE_PARENTHESES_ESCAPE_(__VA_ARGS__)
#define CS_PARENTHESIS(...) CS_PARENTHESIS __VA_ARGS__
#define CS_VANISH_CS_PARENTHESIS

#ifdef __NVCC__
#define CS_CUDA_HOST __host__
#define CS_CUDA_DEVICE __device__
#define CS_CUDA_HOST_DEVICE __host__ __device__
#else
#define CS_CUDA_HOST
#define CS_CUDA_DEVICE
#define CS_CUDA_HOST_DEVICE
#endif

/*============================================================================
 * Type definitions
 *============================================================================*/

/*!
 * Provide default implementations of a csContext based on parallel_for
 * function. This class is a mixin that use CRTP (Curiously Recurring
 * Template Pattern) to provide such functions.
 */

template <class Derived>
class csDispatchContextMixin {
public:

  // Loop over all internal faces
  template <class F, class... Args>
  decltype(auto) parallel_for_i_faces(const cs_mesh_t*  m,
                                      F&&               f,
                                      Args&&...         args);

  // Loop over all boundary faces
  template <class F, class... Args>
  decltype(auto) parallel_for_b_faces(const cs_mesh_t*  m,
                                      F&&               f,
                                      Args&&...         args);

  // Loop over all cells
  template <class F, class... Args>
  decltype(auto) parallel_for_cells(const cs_mesh_t* m, F&& f, Args&&... args);

  // Loop over n elements
  // Must be redefined by the child class
  template <class F, class... Args>
  decltype(auto) parallel_for(cs_lnum_t n, F&& f, Args&&... args) = delete;
};

// Default implementation of parallel_for_i_faces based on parallel_for
template <class Derived>
template <class F, class... Args>
decltype(auto) csDispatchContextMixin<Derived>::parallel_for_i_faces
  (const cs_mesh_t* m, F&& f, Args&&... args) {
  return static_cast<Derived*>(this)->parallel_for
                                        (m->n_i_faces,
                                         static_cast<F&&>(f),
                                         static_cast<Args&&>(args)...);
}

// Default implementation of parallel_for_b_faces based on parallel_for
template <class Derived>
template <class F, class... Args>
decltype(auto) csDispatchContextMixin<Derived>::parallel_for_b_faces
  (const cs_mesh_t* m, F&& f, Args&&... args) {
  return static_cast<Derived*>(this)->parallel_for
                                        (m->n_b_faces,
                                         static_cast<F&&>(f),
                                         static_cast<Args&&>(args)...);
}

// Default implementation of parallel_for_cells based on parallel_for
template <class Derived>
template <class F, class... Args>
decltype(auto) csDispatchContextMixin<Derived>::parallel_for_cells
  (const cs_mesh_t* m, F&& f, Args&&... args) {
  return static_cast<Derived*>(this)->parallel_for
                                        (m->n_cells,
                                         static_cast<F&&>(f),
                                         static_cast<Args&&>(args)...);
}

/*!
 * csContext to execute loops with OpenMP on the CPU
 */
class cs_host_context : public csDispatchContextMixin<cs_host_context> {

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
  bool parallel_for(cs_lnum_t n, F&& f, Args&&... args) {
#   pragma omp parallel for  if (n >= n_min_for_threads)
    for (cs_lnum_t i = 0; i < n; ++i) {
      f(i, args...);
    }
    return true;
  }

  // Loop over the interior faces of a mesh using a specific numbering
  // that avoids conflicts between threads.
  template <class F, class... Args>
  bool parallel_for_i_faces(const cs_mesh_t* m, F&& f, Args&&... args) {
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
        f(f_id, args...);
      }
    }
    return true;
  }
};

#ifdef __NVCC__

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

/*!
 * csContext to execute loops with CUDA on the device
 */
class cs_device_context : public csDispatchContextMixin<cs_device_context> {

private:

  long  grid_size;        /*!< Associated grid size; if <= 0, each kernel
                            launch will use a grid size based on
                            the number of elements. */
  long  block_size;       /*!< Associated block size */
  cudaStream_t  stream;   /*!< Associated CUDA stream */
  int   device;           /*!< Associated CUDA device id */

  cs_lnum_t  n_min_for_device;  /*!< Run on CPU under this threshold */

public:

  cs_device_context(void)
    : grid_size(0), block_size(0), stream(nullptr), device(0),
      n_min_for_device(0)
  {
    block_size = 256;
    device = cs_base_cuda_get_device();
  }

  cs_device_context(long  grid_size,
                    long  block_size,
                    cudaStream_t  stream,
                    int  device)
    : grid_size(grid_size), block_size(block_size), stream(stream),
      device(device), n_min_for_device(0)
  {}

  cs_device_context(long  grid_size,
                    long  block_size,
                    cudaStream_t  stream)
    : grid_size(grid_size), block_size(block_size), stream(stream),
      device(0), n_min_for_device(0)
  {
    block_size = block_size;
    grid_size = grid_size;
    device = cs_base_cuda_get_device();
  }

  cs_device_context(long  grid_size,
                    long  block_size)
    : grid_size(grid_size), block_size(block_size), stream(cudaStreamLegacy),
      device(0), n_min_for_device(0)
  {
    device = cs_base_cuda_get_device();
  }

  cs_device_context(cudaStream_t  stream)
    : grid_size(0), block_size(0), stream(stream), device(0),
      n_min_for_device(0)
  {
    block_size = 256;
    device = cs_base_cuda_get_device();
    printf("init device %d\n", device);
  }

  // Change grid_size configuration, but keeps the stream and device
  void set_cuda_config(long grid_size,
                       long block_size) {
    this->grid_size = grid_size;
    this->block_size = block_size;
  }

  // Change stream, but keeps the grid and device configuration
  void set_cuda_config(cudaStream_t stream) {
    this->stream = stream;
  }

  // Change grid_size configuration and stream, but keeps the device
  void set_cuda_config(long grid_size,
                       long block_size,
                       cudaStream_t stream) {
    this->grid_size = grid_size;
    this->block_size = block_size;
    this->stream = stream;
  }

  // Change grid_size configuration, stream and device
  void set_cuda_config(long grid_size,
                       long block_size,
                       cudaStream_t stream,
                       int device) {
    this->grid_size = grid_size;
    this->block_size = block_size;
    this->stream = stream;
    this->device = device;
  }

  /*! Set minimum number of elements threshold for GPU execution */
  void set_n_min_for_gpu(cs_lnum_t  n) {
    this->n_min_for_device = n;
  }

public:

  // Try to iterate on the GPU and return false if the GPU is not available
  template <class F, class... Args>
  bool parallel_for(cs_lnum_t n, F&& f, Args&&... args) {
    if (device < 0 || n < n_min_for_device) {
      return false;
    }

    long l_grid_size = grid_size;
    if (l_grid_size < 1) {
      l_grid_size = (n % block_size) ?  n/block_size + 1 : n/block_size;
    }

    cs_cuda_kernel_parallel_for<<<l_grid_size, block_size, 0, stream>>>
      (n, static_cast<F&&>(f), static_cast<Args&&>(args)...);

    return true;
  }
};

#else

/*!
 * placeholder to device context
 */
class cs_device_context : public csDispatchContextMixin<cs_device_context> {

public:

  cs_device_context(void)
  {}

  // Change grid_size configuration, but keeps the stream and device
  void set_cuda_config([[maybe_unused]] long grid_size,
                       [[maybe_unused]] long block_size) {
  }

  /*! Set minimum number of elements threshold for GPU execution */
  void set_n_min_for_gpu([[maybe_unused]] cs_lnum_t  n_min_for_device) {
  }

public:

  // Return false
  template <class F, class... Args>
  bool parallel_for(cs_lnum_t n, F&& f, Args&&... args) {
    return false;
  }
};

#endif

/*!
 * csContext that is a combination of multiple contexts.
 * This context will try every context in order, until one actually runs.
 */

template <class... Contexts>
class csCombinedContext : public csDispatchContextMixin<csCombinedContext<Contexts...>>,
                          public Contexts... {

private:
  using mixin_t = csDispatchContextMixin<csCombinedContext<Contexts...>>;

public:
  csCombinedContext() = default;
  csCombinedContext(Contexts... contexts)
    : Contexts(std::move(contexts))...
  {}

public:

  template <class F, class... Args>
    void parallel_for_i_faces(const cs_mesh_t* m, F&& f, Args&&... args) {
    bool launched = false;
    decltype(nullptr) try_execute[] = {
      (   launched = launched
       || Contexts::parallel_for_i_faces(m, f, args...), nullptr)...
    };
    // TODO: raise an error if all contexts return false
  }

  template <class F, class... Args>
    auto parallel_for_b_faces(const cs_mesh_t* m, F&& f, Args&&... args) {
    bool launched = false;
    decltype(nullptr) try_execute[] = {
      (   launched = launched
       || Contexts::parallel_for_b_faces(m, f, args...), nullptr)...
    };
    // TODO: raise an error if all contexts return false
  }

  template <class F, class... Args>
    auto parallel_for_cells_(const cs_mesh_t* m, F&& f, Args&&... args) {
    bool launched = false;
    decltype(nullptr) try_execute[] = {
      (  launched = launched
       || Contexts::parallel_for_cells(m, f, args...), nullptr)...
    };
    // TODO: raise an error if all contexts return false
  }

  template <class F, class... Args>
    auto parallel_for(cs_lnum_t n, F&& f, Args&&... args) {
    bool launched = false;
    decltype(nullptr) try_execute[] = {
      (   launched = launched
       || Contexts::parallel_for(n, f, args...), nullptr)...
    };
    // TODO: raise an error if all contexts return false
  }
};

/*!
 * Default cs_dispatch_context that is a combination of (CPU) and GPU (CUDA)
 * context if available.
 */

class cs_dispatch_context : public csCombinedContext<
  cs_device_context,
  cs_host_context
>
{

private:
  using base_t = csCombinedContext<
  cs_device_context,
  cs_host_context
>;

public:
  using base_t::base_t;
  using base_t::operator=;

};

/*=============================================================================
 * Global variable definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

#endif /* __cplusplus */

#endif /* __CS_DISPATCH_H__ */
