/*============================================================================
 * Array handling utilities.
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C and C++ library headers
 *----------------------------------------------------------------------------*/

#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"
#include "base/cs_base_cuda.h"

/*----------------------------------------------------------------------------*/

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public inline function prototypes
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/* Default kernel that loops over an integer range and calls a device functor.
   This kernel uses a grid_size-stride loop and thus guarantees that all
   integers are processed, even if the grid is smaller.
   All arguments *must* be passed by value to avoid passing CPU references
   to the GPU. */

template <typename T, size_t stride>
__global__ void
cuda_kernel_set_value(cs_lnum_t    n,
                      const T      ref_val,
                      const int    size_arrs,
                      T          **array_ptrs)
{
  // grid_size-stride loop
  for (cs_lnum_t id = blockIdx.x * blockDim.x + threadIdx.x; id < n;
       id += blockDim.x * gridDim.x) {
    for (int j = 0; j < size_arrs; j++) {
      for (size_t k = 0; k < stride; k++) {
        array_ptrs[j][id*stride + k] = ref_val;
      } // end loop stride
    }
  } // end loop arrays
}

template <typename T, size_t stride>
__global__ void
cuda_kernel_set_value(cs_lnum_t    n,
                      const T     *ref_val,
                      const int    size_arrs,
                      T          **array_ptrs)
{
  // grid_size-stride loop
  for (cs_lnum_t id = blockIdx.x * blockDim.x + threadIdx.x; id < n;
       id += blockDim.x * gridDim.x) {
    for (int j = 0; j < size_arrs; j++) {
      for (size_t k = 0; k < stride; k++) {
        array_ptrs[j][id*stride + k] = ref_val[k];
      } // end loop stride
    }
  } // end loop arrays
}

template <typename T, size_t stride>
__global__ void
cuda_kernel_set_value(cs_lnum_t   n,
                      const T     ref_val,
                      T          *array)
{
  // grid_size-stride loop
  for (cs_lnum_t id = blockIdx.x * blockDim.x + threadIdx.x; id < n;
       id += blockDim.x * gridDim.x) {
    for (size_t k = 0; k < stride; k++) {
      array[id*stride + k] = ref_val;
    } // end loop stride
  } // end loop arrays
}

template <typename T, size_t stride>
__global__ void
cuda_kernel_set_value(cs_lnum_t   n,
                      const T    *ref_val,
                      T          *array)
{
  // grid_size-stride loop
  for (cs_lnum_t id = blockIdx.x * blockDim.x + threadIdx.x; id < n;
       id += blockDim.x * gridDim.x) {
    for (size_t k = 0; k < stride; k++) {
      array[id*stride + k] = *ref_val;
    } // end loop stride
  } // end loop arrays
}

template <typename T, size_t stride>
__global__ void
cs_cuda_kernel_set_zero(cs_lnum_t   n,
                        T          *array)
{
  T c = static_cast<T>(0);

  // grid_size-stride loop
  for (cs_lnum_t id = blockIdx.x * blockDim.x + threadIdx.x; id < n;
       id += blockDim.x * gridDim.x) {
    for (size_t k = 0; k < stride; k++) {
      array[id*stride + k] = c;
    } // end loop stride
  } // end loop arrays
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign values to all elements of multiple arrays. ref_val is input
 *        as a pointer or an array
 *
 * Template parmeters.
 *                 T       type name
 *                 stride  1 for scalars, 3 for vectors, 6 for symetric tensors
 *                 Arrays  varadiac parameters pack
 *
 * Function parameters:
 * \param[in]      stream  cuda stream used for the operation
 * \param[in]      n_elts  total number of elements to set
 * \param[in]      ref_val value to assign
 * \param[out]     arrays  arrays to set
 */
/*----------------------------------------------------------------------------*/

template <typename T, size_t stride, typename... Arrays>
void
cs_arrays_set_value(cudaStream_t     stream,
                    bool             async,
                    cs_lnum_t        n_elts,
                    T               *ref_val,
                    Arrays&&...      arrays)
{
  const long block_size_ = 256;
  const long grid_size_ = (n_elts % block_size_) ?
                           n_elts/block_size_ + 1 : n_elts/block_size_;

  /* Explicit expansion of arrays */
  const int size_arrs = sizeof...(arrays);
  T* _array_ptrs[] = {arrays ...};

  if (size_arrs > 1) {
    cudaGraph_t graph;
    cudaGraphExec_t graph_exec = NULL;

    cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);
    for (int j = 0; j < size_arrs; j++) {
      cuda_kernel_set_value<T, stride><<<grid_size_, block_size_, 0, stream>>>
        (n_elts, ref_val, _array_ptrs[j]);
    }
    cudaStreamEndCapture(stream, &graph);

    cudaError_t status = cudaGraphInstantiate(&graph_exec, graph,
                                              nullptr, nullptr, 0);
    cudaGraphLaunch(graph_exec, stream);

    /* For multiple arrays, sync even if async is required,
       to ensure the graph is not destroyed before all kernels
       are scheduled */

    CS_CUDA_CHECK(cudaStreamSynchronize(stream));
    CS_CUDA_CHECK(cudaGetLastError());

    cudaGraphDestroy(graph);
    cudaGraphExecDestroy(graph_exec);
  }

  else {
    cuda_kernel_set_value<T, stride><<<grid_size_, block_size_, 0, stream>>>
      (n_elts, ref_val, _array_ptrs[0]);

    if (async == false) {
      CS_CUDA_CHECK(cudaStreamSynchronize(stream));
      CS_CUDA_CHECK(cudaGetLastError());
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign values to all elements of multiple arrays. ref_val is input
 *        as a scalar value.
 *
 * Template parmeters.
 *                 T       type name
 *                 stride  1 for scalars, 3 for vectors, 6 for symetric tensors
 *                 Arrays  varadiac parameters pack
 *
 * Function parameters:
 * \param[in]      stream  cuda stream used for the operation
 * \param[in]      n_elts  total number of elements to set
 * \param[in]      ref_val value to assign
 * \param[out]     arrays  arrays to set
 */
/*----------------------------------------------------------------------------*/

template <typename T, size_t stride, typename... Arrays>
void
cs_arrays_set_value(cudaStream_t     stream,
                    bool             async,
                    cs_lnum_t        n_elts,
                    T                ref_val,
                    Arrays&&...      arrays)
{
  const long block_size_ = 256;
  const long grid_size_ = (n_elts % block_size_) ?
                           n_elts/block_size_ + 1 : n_elts/block_size_;

  /* Explicit expansion of arrays */
  const int size_arrs = sizeof...(arrays);
  T* _array_ptrs[] = {arrays ...};

  if (size_arrs > 1) {
    cudaGraph_t graph;
    cudaGraphExec_t graph_exec = NULL;

    cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);

    for (int j = 0; j < size_arrs; j++) {
      cuda_kernel_set_value<T, stride><<<grid_size_, block_size_, 0, stream>>>
        (n_elts, ref_val, _array_ptrs[j]);
    }
    cudaStreamEndCapture(stream, &graph);

    cudaError_t status = cudaGraphInstantiate(&graph_exec, graph,
                                              nullptr, nullptr, 0);
    cudaGraphLaunch(graph_exec, stream);

    /* For multiple arrays, sync even if async is required,
       to ensure the graph is not destroyed before all kernels
       are scheduled */

    CS_CUDA_CHECK(cudaStreamSynchronize(stream));
    CS_CUDA_CHECK(cudaGetLastError());

    cudaGraphDestroy(graph);
    cudaGraphExecDestroy(graph_exec);
  }

  else {
    cuda_kernel_set_value<T, stride><<<grid_size_, block_size_, 0, stream>>>
      (n_elts, ref_val, _array_ptrs[0]);

    if (async == false) {
      CS_CUDA_CHECK(cudaStreamSynchronize(stream));
      CS_CUDA_CHECK(cudaGetLastError());
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign values to all elements of multiple arrays. ref_val is input
 *        as a pointer or an array
 *
 * Template parmeters.
 *                 T       type name
 *                 stride  1 for scalars, 3 for vectors, 6 for symetric tensors
 *                 Arrays  varadiac parameters pack
 *
 * Function parameters:
 * \param[in]      stream  cuda stream used for the operation
 * \param[in]      n_elts  total number of elements to set
 * \param[in]      ref_val value to assign
 * \param[out]     arrays  arrays to set
 */
/*----------------------------------------------------------------------------*/

template <typename T, size_t stride, typename... Arrays>
void
cs_arrays_set_zero(cudaStream_t     stream,
                   bool             async,
                   cs_lnum_t        n_elts,
                   Arrays&&...      arrays)
{
  const long block_size_ = 256;
  const long grid_size_ = (n_elts % block_size_) ?
                           n_elts/block_size_ + 1 : n_elts/block_size_;

  /* Explicit expansion of arrays */
  const int size_arrs = sizeof...(arrays);
  T* _array_ptrs[] = {arrays ...};

  if (size_arrs > 1) {
    cudaGraph_t graph;
    cudaGraphExec_t graph_exec = NULL;

    cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);
    for (int j = 0; j < size_arrs; j++) {
      cs_cuda_kernel_set_zero<T, stride><<<grid_size_, block_size_, 0, stream>>>
        (n_elts, _array_ptrs[j]);
    }
    cudaStreamEndCapture(stream, &graph);

    cudaError_t status = cudaGraphInstantiate(&graph_exec, graph,
                                              nullptr, nullptr, 0);
    cudaGraphLaunch(graph_exec, stream);

    /* For multiple arrays, sync even if async is required,
       to ensure the graph is not destroyed before all kernels
       are scheduled */

    CS_CUDA_CHECK(cudaStreamSynchronize(stream));
    CS_CUDA_CHECK(cudaGetLastError());

    cudaGraphDestroy(graph);
    cudaGraphExecDestroy(graph_exec);
  }

  else {
    cs_cuda_kernel_set_zero<T, stride><<<grid_size_, block_size_, 0, stream>>>
      (n_elts, _array_ptrs[0]);

    if (async == false) {
      CS_CUDA_CHECK(cudaStreamSynchronize(stream));
      CS_CUDA_CHECK(cudaGetLastError());
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy values from an array to another of the same dimensions.
 *
 * Template parmeters.
 *                 T       type name
 *
 * \param[in]   stream  cuda stream used for the operation
 * \param[in]   size    number of elements * dimension
 * \param[in]   src     source array values
 * \param[out]  dest    destination array values
 */
/*----------------------------------------------------------------------------*/

template <typename T>
void
cs_array_copy(cudaStream_t     stream,
              const cs_lnum_t  size,
              const T*         src,
              T*               dest)
{
  cudaMemcpyAsync(dest, src, size*sizeof(T),
                  cudaMemcpyDeviceToDevice, stream);
}

/*----------------------------------------------------------------------------*/
