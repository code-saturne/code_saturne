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
#define CS_CUDA_HOST __host__
#define CS_CUDA_DEVICE __device__
#define CS_CUDA_HOST_DEVICE __host__ __device__
#else
#define CS_CUDA_HOST
#define CS_CUDA_DEVICE
#define CS_CUDA_HOST_DEVICE
#endif

// Define a functor that is implemented both on device and host.
#define CS_HOST_DEVICE_FUNCTOR(capture, args, ...) \
  ([CS_REMOVE_PARENTHESES(capture)] CS_CUDA_HOST_DEVICE args __VA_ARGS__)

/*============================================================================
 * Type definitions
 *============================================================================*/

/*!
 * Provide default implementations of a csContext based on iter function.
 * This class is a mixin that use CRTP (Curiously Recurring Template Pattern)
 * to provide such functions.
 */

template <class Derived>
class csDispatchContextMixin {
  public:
    // Iterate the mesh over the internal faces
    template <class F, class... Args>
    decltype(auto) iter_i_faces(const cs_mesh_t* m, F&& f, Args&&... args);

    // Iterate the mesh over the boundary faces
    template <class F, class... Args>
    decltype(auto) iter_b_faces(const cs_mesh_t* m, F&& f, Args&&... args);

    // Iterate the mesh over the cells
    template <class F, class... Args>
    decltype(auto) iter_cells(const cs_mesh_t* m, F&& f, Args&&... args);

    // Iterate over the integers between 0 and n
    // Must be redefined by the child class
    template <class F, class... Args>
    decltype(auto) iter(cs_lnum_t n, F&& f, Args&&... args) = delete;
};

// Default implementation of iter_i_faces based on iter
template <class Derived>
template <class F, class... Args>
decltype(auto) csDispatchContextMixin<Derived>::iter_i_faces
  (const cs_mesh_t* m, F&& f, Args&&... args) {
  return static_cast<Derived*>(this)->iter(m->n_i_faces,
                                           static_cast<F&&>(f),
                                           static_cast<Args&&>(args)...);
}

// Default implementation of iter_b_faces based on iter
template <class Derived>
template <class F, class... Args>
decltype(auto) csDispatchContextMixin<Derived>::iter_b_faces
  (const cs_mesh_t* m, F&& f, Args&&... args) {
  return static_cast<Derived*>(this)->iter(m->n_b_faces,
                                           static_cast<F&&>(f),
                                           static_cast<Args&&>(args)...);
}

// Default implementation of iter_cells based on iter
template <class Derived>
template <class F, class... Args>
decltype(auto) csDispatchContextMixin<Derived>::iter_cells
  (const cs_mesh_t* m, F&& f, Args&&... args) {
  return static_cast<Derived*>(this)->iter(m->n_cells,
                                           static_cast<F&&>(f),
                                           static_cast<Args&&>(args)...);
}

/*!
 * csContext to execute loops with OpenMP on the CPU
 */
class csOpenMpContext : public csDispatchContextMixin<csOpenMpContext> {
  public:
    // iter using a plain omp parallel for
    template <class F, class... Args>
    bool iter(cs_lnum_t n, F&& f, Args&&... args) {
      #pragma omp parallel for
      for (cs_lnum_t i = 0; i < n; ++i) {
        f(i, args...);
      }
      return true;
    }

  // loop over the internal faces of a mesh using a specific numbering
  //that avoids conflicts between threads
  template <class F, class... Args>
  bool iter_i_faces(const cs_mesh_t* m, F&& f, Args&&... args) {
    const int n_i_groups                    = m->i_face_numbering->n_groups;
    const int n_i_threads                   = m->i_face_numbering->n_threads;
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
};

#ifdef __NVCC__

/* Default kernel that loops over an integer range and calls a device functor.
   This kernel uses a grid-stride loop and thus guarantees that all integers
   are processed, even if the grid is smaller.
   All arguments *must* be passed by value to avoid passing CPU references
   to the GPU. */

template <class F, class... Args>
__global__ void cuda_kernel_iter_n(cs_lnum_t n, F f, Args... args) {
  // grid-stride loop
  for (cs_lnum_t id = blockIdx.x * blockDim.x + threadIdx.x; id < n;
       id += blockDim.x * gridDim.x) {
    f(id, args...);
  }
}

/*!
 * csContext to execute loops with CUDA on the GPU
 */
class csCudaContext : public csDispatchContextMixin<csCudaContext> {
  private:
    long grid;
    long block;
    cudaStream_t stream;
    int device;
  public:
    csCudaContext() : grid(0), block(0), stream(nullptr), device(-1) {}
    csCudaContext(long grid, long block, cudaStream_t stream, int device) : grid(grid), block(block), stream(stream), device(device) {}
    csCudaContext(long grid, long block, cudaStream_t stream) : grid(grid), block(block), stream(stream), device(0) {
      cudaGetDevice(&device);
    }
    csCudaContext(long grid, long block) : grid(grid), block(block), stream(cudaStreamLegacy), device(0) {
      cudaGetDevice(&device);
    }

    // Change grid configuration, but keeps the stream and device
    void set_cuda_config(long grid, long block) {
      this->grid = grid;
      this->block = block;
    }

    // Change grid configuration and stream, but keeps the device
    void set_cuda_config(long grid, long block, cudaStream_t stream) {
      this->grid = grid;
      this->block = block;
      this->stream = stream;
    }

    // Change grid configuration, stream and device
    void set_cuda_config(long grid, long block, cudaStream_t stream, int device) {
      this->grid = grid;
      this->block = block;
      this->stream = stream;
      this->device = device;
    }
  public:
    // Try to iterate on the GPU and return false if the GPU is not available
    template <class F, class... Args>
    bool iter(cs_lnum_t n, F&& f, Args&&... args) {
      if (device < 0) {
        return false;
      }
      cuda_kernel_iter_n<<<grid, block, 0, stream>>>(n, static_cast<F&&>(f), static_cast<Args&&>(args)...);
    }
};
#endif

/**
 * csContext that is a combination of multiple contexts.
 * This context will try every context in order, until one actually succeed to execute.
 */
template <class... Contexts>
class csCombinedContext : public csDispatchContextMixin<csCombinedContext<Contexts...>>, public Contexts... {
  private:
    using mixin_t = csDispatchContextMixin<csCombinedContext<Contexts...>>;
  public:
    csCombinedContext() = default;
    csCombinedContext(Contexts... contexts) : Contexts(std::move(contexts))... {}
  public:
    template <class F, class... Args>
    void iter_i_faces(const cs_mesh_t* m, F&& f, Args&&... args) {
      bool executed = false;
      decltype(nullptr) try_execute[] = {
        (executed = executed || Contexts::iter_i_faces(m, f, args...), nullptr)...
      };
      // TODO: raise an error if all contexts return false
    }
    template <class F, class... Args>
    auto iter_b_faces(const cs_mesh_t* m, F&& f, Args&&... args) {
      bool executed = false;
      decltype(nullptr) try_execute[] = {
        (executed = executed || Contexts::iter_b_faces(m, f, args...), nullptr)...
      };
      // TODO: raise an error if all contexts return false
    }
    template <class F, class... Args>
    auto iter_cells_(const cs_mesh_t* m, F&& f, Args&&... args) {
      bool executed = false;
      decltype(nullptr) try_execute[] = {
        (executed = executed || Contexts::iter_cells(m, f, args...), nullptr)...
      };
      // TODO: raise an error if all contexts return false
    }
    template <class F, class... Args>
    auto iter(cs_lnum_t n, F&& f, Args&&... args) {
      bool executed = false;
      decltype(nullptr) try_execute[] = {
        (executed = executed || Contexts::iter(n, f, args...), nullptr)...
      };
      // TODO: raise an error if all contexts return false
    }
};

/**
 * Default csContext that is a combination of OpenMP (CPU) and CUDA if available
 */
class csContext : public csCombinedContext<
#ifdef __NVCC__
  csCudaContext,
#endif
  csOpenMpContext
> {
  private:
    using base_t = csCombinedContext<
#ifdef __NVCC__
  csCudaContext,
#endif
  csOpenMpContext
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

#if 0
/**
 * Other examples of dispatch
 */
template <class F, class... Args>
void
cpu_iter_face(const cs_mesh_t *m, F f, Args... args)
{
  const int n_i_groups                    = m->i_face_numbering->n_groups;
  const int n_i_threads                   = m->i_face_numbering->n_threads;
  const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;

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
}

template <class F, class... Args>
void
cpu_iter_cell(const cs_mesh_t *m, F f, Args... args)
{
  const cs_lnum_t n_cells = m->n_cells;
#pragma omp parallel for
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    f(c_id, args...);
  }
}

template <class F, class... Args>
__global__ void
gpu_iter_kernel(cs_lnum_t n, F f, Args... args)
{
  for (cs_lnum_t id = blockIdx.x * blockDim.x + threadIdx.x; id < n;
       id += blockDim.x * gridDim.x) {
    f(id, args...);
  }
}
template <class F, class... Args>
void
gpu_iter_face(const cs_mesh_t *m,
              int              block,
              int              grid,
              int              stream,
              F                f,
              Args... args)
{
  gpu_iter_kernel<<<block, grid> > >(m->n_i_faces, f, args...);
}
template <class F, class... Args>
void
gpu_iter_cell(const cs_mesh_t *m,
              int              block,
              int              grid,
              int              stream,
              F                f,
              Args... args)
{
  gpu_iter_kernel<<<block, grid, 0, stream> > >(m->n_cells, f, args...);
}

// // Example:
//
// cpu_iter_cell(m, [=] (cs_lnum_t c_id) {
//   // Do some stuff with c_id
// });
//
// gpu_iter_cell(m, bloc, grid, 0, [=] __device__ (cs_lnum_t c_id) {
//   // Do some stuff with c_id
// });

template <class... Args> class Kernel {
public:
  virtual void operator()(const cs_mesh_t *m, Args... args) = 0;
  virtual ~Kernel() {}
};

template <class Derived, class... Args> class CpuFaceKernel : Kernel<Args...> {
public:
  void
  operator()(const cs_mesh_t *m, Args... args) override
  {
    const int n_i_groups                    = m->i_face_numbering->n_groups;
    const int n_i_threads                   = m->i_face_numbering->n_threads;
    const cs_lnum_t *restrict i_group_index = m->i_face_numbering->group_index;

    const cs_lnum_2_t *restrict i_face_cells
      = (const cs_lnum_2_t *restrict)m->i_face_cells;
    const cs_lnum_t *restrict b_face_cells
      = (const cs_lnum_t *restrict)m->b_face_cells;

    for (int g_id = 0; g_id < n_i_groups; g_id++) {

#pragma omp parallel for
      for (int t_id = 0; t_id < n_i_threads; t_id++) {
        for (cs_lnum_t f_id = i_group_index[(t_id * n_i_groups + g_id) * 2];
             f_id < i_group_index[(t_id * n_i_groups + g_id) * 2 + 1];
             f_id++) {
          static_cast<Derived *>(this)->call(f_id, static_cast<Args>(args)...);
        }
      }
    }
  }
};

template <class Derived, class... Args> class CpuCellKernel : Kernel<Args...> {
public:
  void
  operator()(const cs_mesh_t *m, Args... args) override
  {
    const cs_lnum_t n_cells = m->n_cells;
#pragma omp parallel for
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      static_cast<Derived *>(this)->call(c_id, static_cast<Args>(args)...);
    }
  }
};

template <class... Args> class GpuKernel : Kernel<Args...> {
protected:
  int block;
  int grid;
  int stream;

protected:
  template <class Derived>
  __global__ static void
  kernel(cs_lnum_t n, Args... args)
  {
    for (cs_lnum_t id = blockIdx.x * blockDim.x + threadIdx.x; id < n;
         id += blockDim.x * gridDim.x) {
      Derived::call(id, static_cast<Args>(args)...);
    }
  }

public:
  GpuKernel(int block, int grid, int stream)
    : block(block), grid(grid), stream(stream)
  {
  }
};

template <class Derived, class... Args>
class GpuFaceKernel : GpuKernel<Args...> {
public:
  using GpuKernel<Args...>::GpuKernel;

  void
  operator()(const cs_mesh_t *m, Args... args) override
  {
    kernel<Derived>
      <<<block, grid, 0, stream> > >(m.n_i_faces, static_cast<Args>(args)...);
  }
};

template <class Derived, class... Args>
class GpuCellKernel : GpuKernel<Args...> {
public:
  using GpuKernel<Args...>::GpuKernel;

  void
  operator()(const cs_mesh_t *m, Args... args) override
  {
    kernel<Derived>
      <<<block, grid, 0, stream> > >(m.n_cells, static_cast<Args>(args)...);
  }
};

// // Example:
//
// class MyCpuKernel : CpuFaceKernel<MyCpuKernel, ...> {
// public:
//   static void call(cs_lnum_t f_id, ...) {
//     // Do something with f_id
//   }
// }
//
// class MyGpuKernel : GpuFaceKernel<MyGpuKernel, ...> {
// public:
//   __device__ static void call(cs_lnum_t f_id, ...) {
//     // Do something with f_id
//   }
// }
//
// class MyKernel : CpuFaceKernel<MyKernel, ...>, GpuFaceKernel<MyKernel,
// ...> {
// public:
//   __host__ __device__ static void call(cs_lnum_t f_id, ...) {
//     // Do something with f_id
//   }
// }

#endif // commented out

#endif /* __cplusplus */

#endif /* __CS_DISPATCH_H__ */
