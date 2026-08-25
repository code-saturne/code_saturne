#pragma once

/*============================================================================
 * Wrappers for atomic operations portability
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

#if __has_include(<version>)
#include <version>
#endif
#if defined(__cpp_lib_atomic_ref)
  #include <atomic>
#else
  #define __cpp_lib_atomic_ref 0  // To simplify further tests
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

namespace cs {
namespace atomic {

/*=============================================================================
 * Public inline functions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Atomic compare and swap
 *
 * \param[in, out]  address   address to update
 * \param[in]       compare   reference value
 * \param[in]       value     value to set
 *
 * \return  value of sum before adddition
 */
/*----------------------------------------------------------------------------*/

template <typename T>
CS_F_HOST_DEVICE static inline bool
compare_exchange(T  *address,
                 T   compare,
                 T   value)
{
  // GPU version
  #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
  return = (atomicCAS(&(*address)), compare, value) ? : true : false;

  // CPU version, C++20
  #elif __cpp_lib_atomic_ref >= 201806L  // or __cplusplus >= 202002L
  std::atomic_ref<T> address_{*address};
  return address_.compare_exchange_weak(compare,
                                        value,
                                        std::memory_order_relaxed,
                                        std::memory_order_relaxed);

  // CPU version, gcc intrinsics
  #elif (defined(__GNUC__) || defined(__clang__)) && defined(__ATOMIC_RELAXED)
  return __atomic_compare_exchange_n(address,
                                     &compare,
                                     value,
                                     false, // not weak
                                     __ATOMIC_RELAXED,
                                     __ATOMIC_RELAXED);

  // OpenMP version
  #elif defined(_OPENMP)
  bool r;
  T &address_ = *address;
  #if _OPENMP >= 202111
 #pragma omp atomic compare capture
  {
    r = (address_ == compare); if (r) address_ = value;
  }
  return r;
  #else
  #pragma omp critical
  {
    r = (address_ == compare); if (r) address_ = value;
  }
  #endif
  return r;

  // Unsupported
  #else
  #error "No known implementation for atomic fetch_add."
  return false;

  #endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Atomic fetch and add.
 *
 * \param[in, out]  sum    sum to update
 * \param[in]       value  value to add
 *
 * \return  value of sum before adddition
 */
/*----------------------------------------------------------------------------*/

template <typename T>
CS_F_HOST_DEVICE static inline T
fetch_add(T  *sum,
          T   value)
{
  // GPU version
  #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
  return atomicAdd(&(*sum), value);

  // CPU version, C++20
  #elif __cpp_lib_atomic_ref >= 201806L  // or __cplusplus >= 202002L
  std::atomic_ref<T> sum_{*sum};
  return sum_.fetch_add(value);

  // CPU version, gcc intrinsics
  #elif (defined(__GNUC__) || defined(__clang__)) && defined(__ATOMIC_RELAXED)
  return __atomic_fetch_add(sum, value, __ATOMIC_RELAXED);

  // OpenMP version
  #elif defined(_OPENMP)
  T prev;
  T &sum_ = *sum;
  #pragma omp atomic capture
  {
    prev = sum_;
    sum_ += value;
  }
  return prev;

  // Unsupported
  #else
  #error "No known implementation for atomic fetch_add."

  #endif
}

/*----------------------------------------------------------------------------*/

} // namespace atomic
} // namespace cs
