#ifndef CS_REDUCE_OP_H
#define CS_REDUCE_OP_H

/*============================================================================
 * Reduce operatiosn (sum, min, max) using various runtimes (OpenMP, CUDA, ...)
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
 * Standard C++ library headers
 *----------------------------------------------------------------------------*/

#include <utility>
#include <cmath>

#if defined(SYCL_LANGUAGE_VERSION)
#include <sycl/sycl.hpp>
#endif

#if defined(__HIP_DEVICE_COMPILE__)
#include "hip/hip_runtime.h"
#endif
/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_assert.h"
#include "base/cs_mem.h"
#include "base/cs_atomic.h"

#ifdef __CUDACC__
#include "base/cs_base_cuda.h"
#include "base/cs_cuda_reduce.h"
#include "cs_math_cuda.cuh"
#endif

#ifdef __HIPCC__
#include "base/cs_base_hip.h"
#include "base/cs_hip_reduce.h"
#include "cs_math_hip.h"
#endif

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_reduce_op.h
        Reduce operatiosn (sum, min, max) using various runtimes (OpenMP, CUDA, ...)
*/

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Algorithm classes for indirect operations from a graph's edges to its nodes,
  such as face-to cell sums. */

typedef enum {
  /*! Simple operation (assumes data race-avoiding numbering/coloring) */
  CS_OP_SIMPLE,
  /*! Atomic operation */
  CS_OP_ATOMIC,
  // Aliases to be removed once deprecated code is deleted
  CS_DISPATCH_SUM_SIMPLE = CS_OP_SIMPLE,
  CS_DISPATCH_SUM_ATOMIC = CS_OP_ATOMIC
} cs_reduce_op_type_t;

// Aliases to be removed once deprecated code is deleted
using cs_dispatch_sum_type_t = cs_reduce_op_type_t;

/*=============================================================================
 * Global variable definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

namespace cs {
namespace op {
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
sum
(
  T                    *dest,
  const T               src,
  cs_reduce_op_type_t   sum_type = CS_OP_ATOMIC
)
{
  if (sum_type == CS_OP_ATOMIC) {
#if 1
    using sum_v = assembled_value<T>;
    sum_v v;

    v.get() = src;
    sum_v::ref(*dest).conflict_free_add(-1u, v);
#else
    atomicAdd(dest, src);
#endif
  }
  else if (sum_type == CS_OP_SIMPLE) {
    *dest += src;
  }
}

#elif defined(__HIP_DEVICE_COMPILE__)

template <typename T>
__device__ static void __forceinline__
sum
(
  T                   *dest,
  const T              src,
  cs_reduce_op_type_t  sum_type = CS_OP_ATOMIC
)
{
  if (sum_type == CS_OP_ATOMIC) {
    atomicAdd(dest, src);
  }
  else if (sum_type == CS_OP_SIMPLE) {
    *dest += src;
  }
}

#elif defined(SYCL_LANGUAGE_VERSION)

template <typename T>
inline void
sum
(
  T                   *dest,
  const T              src,
  cs_reduce_op_type_t  sum_type = CS_OP_ATOMIC
)
{
  if (sum_type == CS_OP_SIMPLE) {
    *dest += src;
  }
  else if (sum_type == CS_OP_ATOMIC) {
    sycl::atomic_ref<T,
                     sycl::memory_order::relaxed,
                     sycl::memory_scope::device> aref(*dest);
    aref.fetch_add(src);
  }
}

#else  // ! CUDA or HIP or SYCL

template <typename T>
inline void
sum
(
  T                   *dest,
  const T              src,
  cs_reduce_op_type_t  sum_type = CS_OP_ATOMIC
)
{
  if (sum_type == CS_OP_SIMPLE) {
    *dest += src;
  }
  else if (sum_type == CS_OP_ATOMIC) {
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
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
sum
(
  T                   *dest,
  const T             *src,
  cs_reduce_op_type_t  sum_type = CS_OP_ATOMIC
)
{
  if (sum_type == CS_OP_SIMPLE) {
    for (cs_lnum_t i = 0; i < dim; i++) {
      dest[i] += src[i];
    }
  }
  else if (sum_type == CS_OP_ATOMIC) {
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

#elif defined(__HIP_DEVICE_COMPILE__)

template <size_t dim, typename T>
__device__ static void __forceinline__
sum
(
  T                   *dest,
  const T             *src,
  cs_reduce_op_type_t  sum_type = CS_OP_ATOMIC)
{
  if (sum_type == CS_OP_SIMPLE) {
    for (cs_lnum_t i = 0; i < dim; i++) {
      dest[i] += src[i];
    }
  }
  else if (sum_type == CS_OP_ATOMIC) {
    for (size_t i = 0; i < dim; i++) {
      atomicAdd(&dest[i], src[i]);
    }
  }
}

#elif defined(SYCL_LANGUAGE_VERSION)

template <size_t dim, typename T>
inline void
sum
(
  T                   *dest,
  const T             *src,
  cs_reduce_op_type_t  sum_type = CS_OP_ATOMIC
)
{
  if (sum_type == CS_OP_SIMPLE) {
    for (size_t i = 0; i < dim; i++) {
      dest[i] += src[i];
    }
  }
  else if (sum_type == CS_OP_ATOMIC) {
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
sum
(
  T                   *dest,
  const T             *src,
  cs_reduce_op_type_t  sum_type = CS_OP_ATOMIC
)
{
  if (sum_type == CS_OP_SIMPLE) {
    for (size_t i = 0; i < dim; i++) {
      dest[i] += src[i];
    }
  }
  else if (sum_type == CS_OP_ATOMIC) {
    for (size_t i = 0; i < dim; i++) {
      #ifdef _OPENMP
      #pragma omp atomic
      #endif
      dest[i] += src[i];
    }
  }
}

#endif // __CUDA_ARCH__

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Minimum of values using a chosen dispatch min type.
 *
 * This allows calling a unique function for assembly in lambda functions
 * called by a dispatch context for interior or boundary faces.
 *
 * \tparam  T    value type
 *
 * \param[in, out]  dest      destination
 * \param[in]       src       source value
 * \param[in]       min_type  min type
 */
/*----------------------------------------------------------------------------*/

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)

template <typename T>
__device__ static void __forceinline__
min
(
  T                   *dest,
  const T              src,
  cs_reduce_op_type_t  min_type = CS_OP_ATOMIC
)
{
  if (min_type == CS_OP_ATOMIC) {
    cs::atomic::fetch_min(dest, src);
  }
  else if (min_type == CS_OP_SIMPLE) {
    if (src < *dest)
      *dest = src;
  }
}

#elif defined(SYCL_LANGUAGE_VERSION)

template <typename T>
inline void
min
(
  T                       *dest,
  const T                  src,
  cs_reduce_op_type_t   min_type = CS_OP_ATOMIC
)
{
  if (min_type == CS_OP_SIMPLE) {
    if (src < *dest)
      *dest = src;
  }
  else if (min_type == CS_OP_ATOMIC) {
    cs::atomic::fetch_min(dest, src);
  }
}

#else  // ! CUDA or HIP or SYCL

template <typename T>
inline void
min
(
  T                    *dest,
  const T               src,
  cs_reduce_op_type_t   min_type = CS_OP_ATOMIC
)
{
  if (min_type == CS_OP_SIMPLE) {
    if (src < *dest)
      *dest = src;
  }
  else if (min_type == CS_OP_ATOMIC) {
    cs::atomic::fetch_min(dest, src);
  }
}

#endif // __CUDA_ARCH__

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Minimum of strided array values using a chosen dispatch min type.
 *
 * \tparam  dim       array stride
 * \tparam  T         array type
 *
 * \param[in, out]  dest      destination values
 * \param[in]       src       source values
 * \param[in]       min_type  min type
 */
/*----------------------------------------------------------------------------*/

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)

template <size_t dim, typename T>
__device__ static void __forceinline__
min
(
  T                   *dest,
  const T             *src,
  cs_reduce_op_type_t  min_type = CS_OP_ATOMIC
)
{
  if (min_type == CS_OP_SIMPLE) {
    for (size_t i = 0; i < dim; i++) {
      if (src[i] < dest[i])
        dest[i] = src[i];
    }
  }
  else if (min_type == CS_OP_ATOMIC) {
    for (size_t i = 0; i < dim; i++) {
      cs::atomic::fetch_min(&dest[i], src[i]);
    }
  }
}

#elif defined(SYCL_LANGUAGE_VERSION)

template <size_t dim, typename T>
inline void
min
(
  T                   *dest,
  const T             *src,
  cs_reduce_op_type_t  min_type = CS_OP_ATOMIC
)
{
  if (min_type == CS_OP_SIMPLE) {
    for (size_t i = 0; i < dim; i++) {
      if (src[i] < dest[i])
        dest[i] = src[i];
    }
  }
  else if (min_type == CS_OP_ATOMIC) {
    for (size_t i = 0; i < dim; i++) {
      cs::atomic::fetch_min(&dest[i], src[i]);
    }
  }
}

#else  // ! CUDA or HIP or SYCL

template <size_t dim, typename T>
inline void
min
(
  T                   *dest,
  const T             *src,
  cs_reduce_op_type_t  min_type = CS_OP_ATOMIC
)
{
  if (min_type == CS_OP_SIMPLE) {
    for (size_t i = 0; i < dim; i++) {
      if (src[i] < dest[i])
        dest[i] = src[i];
    }
  }
  else if (min_type == CS_OP_ATOMIC) {
    for (size_t i = 0; i < dim; i++) {
      cs::atomic::fetch_min(&dest[i], src[i]);
    }
  }
}

#endif // __CUDA_ARCH__

} /* namespace op */
} /* namespace cs */
/*----------------------------------------------------------------------------*/

#endif /* CS_DISPATCH_H */
