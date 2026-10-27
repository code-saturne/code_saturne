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

#include <type_traits>

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
  T old = atomicCAS(&(*address), compare, value);
  return  (old == compare) ? true : false;

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

  // SYCL GPU version
  #elif defined(SYCL_LANGUAGE_VERSION) && defined(__SYCL_DEVICE_ONLY__)

  sycl::atomic_ref<T,
                   sycl::memory_order::relaxed,
                   sycl::memory_scope::device> aref(*sum);
  aref.fetch_add(value);

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
/*!
 * \brief Atomic fetch and minimum.
 *
 * \param[in, out]  address  address of value to update
 * \param[in]       value    candidate value
 *
 * \return  value of address before min operation
 */
/*----------------------------------------------------------------------------*/

template <typename T>
CS_F_HOST_DEVICE static inline T
fetch_min(T  *address,
          T   value)
{
  // GPU version
  #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
  if constexpr (std::is_integral<T>::value) {
    if constexpr (sizeof(T) == sizeof(int)) {
      if constexpr (std::is_signed<T>::value)
        return static_cast<T>(atomicMin(reinterpret_cast<int *>(address),
                                        static_cast<int>(value) ) );
      else
        return static_cast<T>(atomicMin(reinterpret_cast<unsigned int *>(address),
                                        static_cast<unsigned int>(value) ) );
    }
    else if constexpr (sizeof(T) == sizeof(unsigned long long int)) {
      if constexpr (std::is_signed<T>::value)
        return static_cast<T>(atomicMin(reinterpret_cast<long long int *>(address),
                                        static_cast<long long int>(value) ) );
      else
        return static_cast<T>(atomicMin(reinterpret_cast<unsigned long long int *>(address),
                                        static_cast<unsigned long long int>(value) ) );
    }
    else {
      T old = *address;
      while (value < old) {
        T assumed = old;
        old = atomicCAS(address, assumed, value);
        if (assumed == old)
          break;
      }
      return old;
    }
  }
  else if constexpr (std::is_floating_point<T>::value && sizeof(T) == sizeof(float)) {
    int *address_as_i = reinterpret_cast<int *>(address);
    int old = *address_as_i, assumed;
    do {
      assumed = old;
      T cur_val = __int_as_float(assumed);
      if (cur_val <= value || !(value == value))
        return cur_val;
      old = atomicCAS(address_as_i, assumed, __float_as_int(value));
    } while (assumed != old);
    return __int_as_float(old);
  }
  else if constexpr (std::is_floating_point<T>::value && sizeof(T) == sizeof(double)) {
    unsigned long long int *address_as_ull =
      reinterpret_cast<unsigned long long int *>(address);
    unsigned long long int old = *address_as_ull, assumed;
    do {
      assumed = old;
      T cur_val = __longlong_as_double(assumed);
      if (cur_val <= value || !(value == value))
        return cur_val;
      old = atomicCAS(address_as_ull, assumed, __double_as_longlong(value));
    } while (assumed != old);
    return __longlong_as_double(old);
  }

  // SYCL GPU version
  #elif defined(SYCL_LANGUAGE_VERSION) && defined(__SYCL_DEVICE_ONLY__)

  sycl::atomic_ref<T,
                   sycl::memory_order::relaxed,
                   sycl::memory_scope::device> aref(*address);
  aref.fetch_min(value);

  // CPU version, C++20
  #elif __cpp_lib_atomic_ref >= 201806L  // or __cplusplus >= 202002L
  std::atomic_ref<T> address_{*address};
  T prev = address_.load(std::memory_order_relaxed);
  while (value < prev && !address_.compare_exchange_weak(prev,
                                                         value,
                                                         std::memory_order_relaxed,
                                                         std::memory_order_relaxed))
    ;
  return prev;

  // CPU version, gcc intrinsics
  #elif (defined(__GNUC__) || defined(__clang__)) && defined(__ATOMIC_RELAXED)
  T prev;
  __atomic_load(address, &prev, __ATOMIC_RELAXED);
  while (value < prev && !__atomic_compare_exchange(address,
                                                    &prev,
                                                    &value,
                                                    true, /* weak */
                                                    __ATOMIC_RELAXED,
                                                    __ATOMIC_RELAXED))
    ;
  return prev;

  // OpenMP version
  #elif defined(_OPENMP)
  T prev;
  T &address_ = *address;
  #if _OPENMP >= 202111
  #pragma omp atomic compare capture
  {
    prev = address_;
    if (address_ > value) address_ = value;
  }
  #else
  #pragma omp critical
  {
    prev = address_;
    if (address_ > value) address_ = value;
  }
  #endif
  return prev;

  // Unsupported
  #else
  #error "No known implementation for atomic fetch_min."
  return value;

  #endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Atomic fetch and maximum.
 *
 * \param[in, out]  address  address of value to update
 * \param[in]       value    candidate value
 *
 * \return  value of address before max operation
 */
/*----------------------------------------------------------------------------*/

template <typename T>
CS_F_HOST_DEVICE static inline T
fetch_max(T  *address,
          T   value)
{
  // GPU version
  #if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
  if constexpr (std::is_integral<T>::value) {
    if constexpr (sizeof(T) == sizeof(int)) {
      if constexpr (std::is_signed<T>::value)
        return static_cast<T>(atomicMax(reinterpret_cast<int *>(address),
                                        static_cast<int>(value) ) );
      else
        return static_cast<T>(atomicMax(reinterpret_cast<unsigned int *>(address),
                                        static_cast<unsigned int>(value) ) );
    }
    else if constexpr (sizeof(T) == sizeof(unsigned long long int)) {
      if constexpr (std::is_signed<T>::value)
        return static_cast<T>(atomicMax(reinterpret_cast<long long int *>(address),
                                        static_cast<long long int>(value) ) );
      else
        return static_cast<T>(atomicMax(reinterpret_cast<unsigned long long int *>(address),
                                        static_cast<unsigned long long int>(value) ) );
    }
    else {
      T old = *address;
      while (old < value) {
        T assumed = old;
        old = atomicCAS(address, assumed, value);
        if (assumed == old)
          break;
      }
      return old;
    }
  }
  else if constexpr (std::is_floating_point<T>::value && sizeof(T) == sizeof(float)) {
    int *address_as_i = reinterpret_cast<int *>(address);
    int old = *address_as_i, assumed;
    do {
      assumed = old;
      T cur_val = __int_as_float(assumed);
      if (cur_val >= value || !(value == value))
        return cur_val;
      old = atomicCAS(address_as_i, assumed, __float_as_int(value));
    } while (assumed != old);
    return __int_as_float(old);
  }
  else if constexpr (std::is_floating_point<T>::value && sizeof(T) == sizeof(double)) {
    unsigned long long int *address_as_ull =
      reinterpret_cast<unsigned long long int *>(address);
    unsigned long long int old = *address_as_ull, assumed;
    do {
      assumed = old;
      T cur_val = __longlong_as_double(assumed);
      if (cur_val >= value || !(value == value))
        return cur_val;
      old = atomicCAS(address_as_ull, assumed, __double_as_longlong(value));
    } while (assumed != old);
    return __longlong_as_double(old);
  }

  // SYCL GPU version
  #elif defined(SYCL_LANGUAGE_VERSION) && defined(__SYCL_DEVICE_ONLY__)

  sycl::atomic_ref<T,
                   sycl::memory_order::relaxed,
                   sycl::memory_scope::device> aref(*address);
  aref.fetch_max(value);

  // CPU version, C++20
  #elif __cpp_lib_atomic_ref >= 201806L  // or __cplusplus >= 202002L
  std::atomic_ref<T> address_{*address};
  T prev = address_.load(std::memory_order_relaxed);
  while (prev < value && !address_.compare_exchange_weak(prev,
                                                         value,
                                                         std::memory_order_relaxed,
                                                         std::memory_order_relaxed))
    ;
  return prev;

  // CPU version, gcc intrinsics
  #elif (defined(__GNUC__) || defined(__clang__)) && defined(__ATOMIC_RELAXED)
  T prev;
  __atomic_load(address, &prev, __ATOMIC_RELAXED);
  while (prev < value && !__atomic_compare_exchange(address,
                                                    &prev,
                                                    &value,
                                                    true, /* weak */
                                                    __ATOMIC_RELAXED,
                                                    __ATOMIC_RELAXED))
    ;
  return prev;

  // OpenMP version
  #elif defined(_OPENMP)
  T prev;
  T &address_ = *address;
  #if _OPENMP >= 202111
  #pragma omp atomic compare capture
  {
    prev = address_;
    if (address_ < value) address_ = value;
  }
  #else
  #pragma omp critical
  {
    prev = address_;
    if (address_ < value) address_ = value;
  }
  #endif
  return prev;

  // Unsupported
  #else
  #error "No known implementation for atomic fetch_max."
  return value;

  #endif
}

/*----------------------------------------------------------------------------*/

} // namespace atomic
} // namespace cs
