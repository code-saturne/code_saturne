/*============================================================================
 * Math functions for HIP.
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

#pragma once

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C/C++/HIP library headers
 *----------------------------------------------------------------------------*/

#include <float.h>
#include "hip/hip_runtime.h"
#include <hip/hip_runtime_api.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_base_hip.h"

/*=============================================================================
 * Semi-private inline functions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the bit mask of the threads that have the same value v.
 *
 * This partitionins the incoming threads into groups whose members have the
 * same v value.
 *
 * This reduced to calling __match_any_sync on recent HIP (>=6.2) architectures.
 *
 * \param[in]  mask
 * \param[in]  v
 *
 * \return
 */
/*----------------------------------------------------------------------------*/

template <class V>
__device__ cs_hip_mask_t
_conflict_mask(cs_hip_mask_t mask,
               V             v) noexcept
{
  // It seems that although the dunder warp functions are available for
  // ROCm 6.2, they may be disabled...
  // Documentation indicates that it is activated by default on ROCm 7.0
#if (HIP_VERSION >= 60200000 && defined(HIP_ENABLE_WARP_SYNC_BUILTINS)) \
  || (HIP_VERSION >= 70000000 && !defined(HIP_DISABLE_WARP_SYNC_BUILTINS))
  return __match_any_sync(mask, v);
#else
  cs_hip_mask_t lanemask_eq = 1u << (threadIdx.x % CS_HIP_WARP_SIZE);

  if (!(mask & lanemask_eq))
    return 0;

  cs_hip_mask_t ref, ballot;
  int leader;

  goto entry;

loop:
  mask &= ~ballot;

entry:
  // Handle different warp sizes
#if CS_HIP_WARP_SIZE == 64
// ffsll seems to have problems distinguishing between long long int and
// long long unsigned...
//  leader = __ffsll(mask) - 1;
  leader = ( mask == 0 ? -1 : __builtin_ctzll(mask) );
#else
  leader = __ffs(mask) - 1;
#endif

  ref = __shfl(v, leader);
  ballot = __ballot(v == ref);

  if (!(ballot & lanemask_eq))
    goto loop;

  return ballot;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reduction step for conflict-free addition.
 *
 * \param[in]       mask  mask for threads in warpp
 * \param[in, out]  v     value which should be reduced between threads
 *
 * \return true for rank 0, false for others.
 */
/*----------------------------------------------------------------------------*/

template <class T>
__device__ bool
_reduce_add(cs_hip_mask_t  mask,
            cs_hip_mask_t  peers,
            T             &v) noexcept
{
  int laneid = __lane_id();   // Ensure portability depending on warpSize for AMD
  cs_hip_mask_t lanemask_lt = (1u << laneid) - 1;
  cs_hip_mask_t lanemask_gt = -2u << laneid;

#if CS_HIP_WARP_SIZE == 64
  int rank = __popcll(peers & lanemask_lt);
#else
  int rank = __popc(peers & lanemask_lt);
#endif

  bool is_leader = (rank == 0);

  peers &= lanemask_gt;

  while (__ballot(peers != 0)) {
#if CS_HIP_WARP_SIZE == 64
    int next = ( peers == 0 ? -1 : __builtin_ctzll(peers) ) + 1;
// ffsll seems to have problems distinguishing between long long int and
// long long unsigned...
//    int next = __ffsll(peers);   // 1-based index
#else
    int next = __ffs(peers);   // 1-based index
#endif

    auto tmp = v.shuffle(mask, next - 1);
    if (next) {
      v.add(tmp);
    }

    peers &= __ballot(!(rank & 1));
    rank >>= 1;
  }

  return is_leader;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Class for value assembly using warp-based partial summation to
 *        reduce conflict occurences in downstream atomic addition.
 */
/*----------------------------------------------------------------------------*/

template <class T, size_t...>
class assembled_value {

private:
  T value = {};

public:
  using inner_type = T;

public:

  // __device__
  assembled_value() noexcept = default;

  __device__
  assembled_value(T val) noexcept : value(val)
  {
  }

  __device__ void
  add(const assembled_value &restrict other) restrict noexcept
  {
    value += other.value;
  }

  __device__ void
  atomic_add(const assembled_value &restrict other) restrict noexcept
  {
    atomicAdd(&value, other.value);
  }

  __device__ assembled_value
  exchange(const assembled_value &restrict other) restrict noexcept
  {
    assembled_value previous = *this;
    *this                    = other;
    return previous;
  }

  __device__ assembled_value
  atomic_exchange(const assembled_value &restrict other) restrict noexcept
  {
    return assembled_value(atomicExch(&value, other.value));
  }

  __device__ assembled_value
  shuffle(cs_hip_mask_t mask, unsigned laneid) const noexcept
  {
#if HIP_VERSION >= 60200000 && defined(HIP_ENABLE_WARP_SYNC_BUILTINS)
    return assembled_value(__shfl_sync(mask, value, laneid));
#else
    CS_UNUSED(mask);
    return assembled_value(__shfl(value, laneid));
#endif
  }

  __device__ cs_hip_mask_t
  conflict_mask(cs_hip_mask_t mask) const noexcept
  {
    return _conflict_mask(mask, (uintptr_t)this);
  }

  __device__ bool
  reduce_add(cs_hip_mask_t mask, cs_hip_mask_t peers) noexcept
  {
    return _reduce_add(mask, peers, *this);
  }

  __device__ void
  conflict_free_add(cs_hip_mask_t mask, assembled_value other) noexcept
  {
    cs_hip_mask_t peers = conflict_mask(mask);
    if (other.reduce_add(mask, peers)) {
      atomic_add(other);
    }
  }

  __device__ inner_type &
  operator*() noexcept
  {
    return value;
  }

  __device__ inner_type const &
  operator*() const noexcept
  {
    return value;
  }

  __device__ inner_type *
  operator->() noexcept
  {
    return &value;
  }

  __device__ inner_type const *
  operator->() const noexcept
  {
    return &value;
  }

  __device__ inner_type &
  get() noexcept
  {
    return value;
  }

  __device__ inner_type const &
  get() const noexcept
  {
    return value;
  }

  static __device__ assembled_value &
  ref(inner_type &r) noexcept
  {
    return reinterpret_cast<assembled_value &>(r);
  }

  static __device__ assembled_value const &
  ref(inner_type const &r) noexcept
  {
    return reinterpret_cast<assembled_value const &>(r);
  }
};

template <class T, size_t Head, size_t... Tail>
class assembled_value<T, Head, Tail...> {

private:
  assembled_value<T, Tail...> data[Head];

public:
  using inner_type = typename assembled_value<T, Tail...>::inner_type[Head];

public:
  // __device__ assembled_value() noexcept = default;

  __device__ void
  add(const assembled_value &restrict other) restrict noexcept
  {
    for (size_t i = 0; i < Head; ++i) {
      data[i].add(other.data[i]);
    }
  }

  __device__ void
  atomic_add(const assembled_value &restrict other) restrict noexcept
  {
    for (size_t i = 0; i < Head; ++i) {
      data[i].atomic_add(other.data[i]);
    }
  }

  __device__ assembled_value
  exchange(const assembled_value &restrict other) restrict noexcept
  {
    assembled_value previous;
    for (size_t i = 0; i < Head; ++i) {
      previous.data[i] = data[i].exchange(other.data[i]);
    }
    return previous;
  }

  __device__ assembled_value
  atomic_exchange(const assembled_value &restrict other) restrict noexcept
  {
    assembled_value previous;
    for (size_t i = 0; i < Head; ++i) {
      previous.data[i] = data[i].atomic_exchange(other.data[i]);
    }
    return previous;
  }

  __device__ assembled_value
  shuffle(cs_hip_mask_t mask, unsigned laneid) const noexcept
  {
    assembled_value shuffled;
    for (size_t i = 0; i < Head; ++i) {
      shuffled.data[i] = data[i].shuffle(mask, laneid);
    }
    return shuffled;
  }

  __device__ cs_hip_mask_t
  conflict_mask(cs_hip_mask_t mask) const noexcept
  {
    return _conflict_mask(mask, (uintptr_t)this);
  }

  __device__ bool
  reduce_add(cs_hip_mask_t mask, cs_hip_mask_t peers) noexcept
  {
    return _reduce_add(mask, peers, *this);
  }

  __device__ void
  conflict_free_add(cs_hip_mask_t mask, assembled_value other) noexcept
  {
    cs_hip_mask_t peers = conflict_mask(mask);
    if (other.reduce_add(mask, peers)) {
      atomic_add(other);
    }
  }

  __device__ assembled_value<T, Tail...> &
  operator[](size_t i) noexcept
  {
    return data[i];
  }

  __device__ assembled_value<T, Tail...> const &
  operator[](size_t i) const noexcept
  {
    return data[i];
  }

  __device__ inner_type &
  get() noexcept
  {
    return reinterpret_cast<inner_type &>(*this);
  }

  __device__ inner_type const &
  get() const noexcept
  {
    return reinterpret_cast<inner_type const &>(*this);
  }

  static __device__ assembled_value &
  ref(inner_type &r) noexcept
  {
    return reinterpret_cast<assembled_value &>(r);
  }

  static __device__ assembled_value const &
  ref(inner_type const &r) noexcept
  {
    return reinterpret_cast<assembled_value const &>(r);
  }
};

/*=============================================================================
 * Device function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Normalise a vector of 3 real values.
 *
 * To normalize in-place, \p vin and \p vout may point to the same array.
 *
 * \param[in]     vin           vector
 * \param[out]    vout          normalized vector
 */
/*----------------------------------------------------------------------------*/

template <typename T>
__device__ static void
cs_math_3_normalize_hip(const T  in[3],
                        T        out[3])
{
  T norm = sqrt(  in[0]*in[0]
                + in[1]*in[1]
                + in[2]*in[2]);

  T inverse_norm =  1. / norm;

  out[0] = inverse_norm * in[0];
  out[1] = inverse_norm * in[1];
  out[2] = inverse_norm * in[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the square norm of a vector of 3 real values.
 *
 * \param[in]     v             vector of 3 real values
 *
 * \return square norm of v.
 */
/*----------------------------------------------------------------------------*/

template <typename T>
__device__ static T
cs_math_3_square_norm_hip(const T  in[3])
{
  T norm = in[0]*in[0] + in[1]*in[1] + in[2]*in[2];
  return norm;
}

/*----------------------------------------------------------------------------*/
