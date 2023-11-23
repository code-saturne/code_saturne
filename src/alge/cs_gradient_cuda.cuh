/*============================================================================
 * Gradient reconstruction, CUDA implementations.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <chrono>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#include <cuda_runtime_api.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"

#include "cs_base_accel.h"
#include "cs_base_cuda.h"
#include "cs_blas.h"
#include "cs_cell_to_vertex.h"
#include "cs_ext_neighborhood.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_adjacencies.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_porous_model.h"
#include "cs_prototypes.h"
#include "cs_timer.h"
#include "cs_timer_stats.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gradient.h"
#include "cs_gradient_priv.h"

__device__ cs_real_t
cs_math_fabs_cuda(cs_real_t  x)
{
  cs_real_t ret = (x <  0) ? -x : x;

  return ret;
}

__device__ cs_real_t
cs_math_3_dot_product_cuda(const cs_real_t  u[3],
                      const cs_real_t  v[3])
{
  cs_real_t prod = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];

  return prod;
}


__device__ void cs_math_3_normalise_cuda(const cs_real_t in[3],
                                         cs_real_t out[3])
{
  cs_real_t norm = sqrt(in[0]*in[0] 
          + in[1]*in[1]
          + in[2]*in[2]);

  cs_real_t inverse_norm =  1. / norm;

  out[0] = inverse_norm * in[0];
  out[1] = inverse_norm * in[1];
  out[2] = inverse_norm * in[2];
}

__device__ cs_real_t cs_math_3_square_norm_cuda(const cs_real_t in[3]){
  cs_real_t norm = in[0]*in[0] + in[1]*in[1] + in[2]*in[2];
  return norm;
}

__device__ void _math_6_inv_cramer_sym_in_place_cuda(cs_cocg_t in[6]){
  cs_real_t in00 = in[1]*in[2] - in[4]*in[4];
  cs_real_t in01 = in[4]*in[5] - in[3]*in[2];
  cs_real_t in02 = in[3]*in[4] - in[1]*in[5];
  cs_real_t in11 = in[0]*in[2] - in[5]*in[5];
  cs_real_t in12 = in[3]*in[5] - in[0]*in[4];
  cs_real_t in22 = in[0]*in[1] - in[3]*in[3];

  cs_real_t det_inv = 1. / (in[0]*in00 + in[3]*in01 + in[5]*in02);

  in[0] = in00 * det_inv;
  in[1] = in11 * det_inv;
  in[2] = in22 * det_inv;
  in[3] = in01 * det_inv;
  in[4] = in12 * det_inv;
  in[5] = in02 * det_inv;
}


template <class V>
__device__ uint32_t _conflict_mask(uint32_t mask, V v) noexcept {
#if __CUDA_ARCH__ >= 700
  return __match_any_sync(mask, v);
#else
  uint32_t lanemask_eq = 1u << (threadIdx.x % 32);
  if (!(mask & lanemask_eq))
    return 0;
  uint32_t ref, ballot;
  int leader;
  goto entry;
loop:
  mask &= ~ballot;
entry:
  leader = __ffs(mask) - 1;
  ref = __shfl_sync(mask, v, leader);
  ballot = __ballot_sync(mask, v == ref);
  if (!(ballot & lanemask_eq))
    goto loop;
  return ballot;
#endif
}

template <class T>
__device__ bool _reduce_add(uint32_t mask, uint32_t peers, T& v) noexcept {
  int laneid = threadIdx.x % 32;
  uint32_t lanemask_lt = (1u << laneid) - 1;
  uint32_t lanemask_gt = -2u << laneid;
  int rank = __popc(peers & lanemask_lt);
  bool is_leader = rank == 0;

  peers &= lanemask_gt;
  while (__any_sync(mask, peers)) {
    int next = __ffs(peers);

    auto tmp = v.shuffle(mask, next - 1);
    if (next) {
      v.add(tmp);
    }

    peers &= __ballot_sync(mask, !(rank & 1));

    rank >>= 1;
  }

  return is_leader;
}


template <class T, size_t...>
class AtomicCell {
  private:
    T value = {};
  public:
    using inner_type = T;
  public:
    __device__ AtomicCell() noexcept = default;
    __device__ AtomicCell(T value) noexcept : value(value) {}
    __device__ void add(const AtomicCell&restrict other) restrict noexcept {
      value += other.value;
    }
    __device__ void atomic_add(const AtomicCell&restrict other) restrict noexcept {
      atomicAdd(&value, other.value);
    }
    __device__ AtomicCell exchange(const AtomicCell&restrict other) restrict noexcept {
      AtomicCell previous = *this;
      *this = other;
      return previous;
    }
    __device__ AtomicCell atomic_exchange(const AtomicCell&restrict other) restrict noexcept {
      return AtomicCell(atomicExch(&value, other.value));
    }
    __device__ AtomicCell shuffle(uint32_t mask, unsigned laneid) const noexcept {
      return AtomicCell(__shfl_sync(mask, value, laneid));
    }
    __device__ uint32_t conflict_mask(uint32_t mask) const noexcept {
      return _conflict_mask(mask, (uintptr_t)this);
    }
    __device__ bool reduce_add(uint32_t mask, uint32_t peers) noexcept {
      return _reduce_add(mask, peers, *this);
    }
    __device__ void conflict_free_add(uint32_t mask, AtomicCell other) noexcept {
      uint32_t peers = conflict_mask(mask);
      if (other.reduce_add(mask, peers)) {
        atomic_add(other);
      }
    }
    __device__ inner_type& operator*() noexcept {
      return value;
    }
    __device__ inner_type const& operator*() const noexcept {
      return value;
    }
    __device__ inner_type* operator->() noexcept {
      return &value;
    }
    __device__ inner_type const* operator->() const noexcept {
      return &value;
    }
    __device__ inner_type& get() noexcept {
      return value;
    }
    __device__ inner_type const& get() const noexcept {
      return value;
    }
    static __device__ AtomicCell& ref(inner_type& r) noexcept {
      return reinterpret_cast<AtomicCell&>(r);
    }
    static __device__ AtomicCell const& ref(inner_type const& r) noexcept {
      return reinterpret_cast<AtomicCell const&>(r);
    }
};

template <class T, size_t Head, size_t... Tail>
class AtomicCell<T, Head, Tail...> {
  private:
    AtomicCell<T, Tail...> data[Head];
  public:
    using inner_type = typename AtomicCell<T, Tail...>::inner_type[Head];
  public:
    __device__ AtomicCell() noexcept = default;
    __device__ void add(const AtomicCell&restrict other) restrict noexcept {
      for (size_t i = 0; i < Head; ++i) {
        data[i].add(other.data[i]);
      }
    }
    __device__ void atomic_add(const AtomicCell&restrict other) restrict noexcept {
      for (size_t i = 0; i < Head; ++i) {
        data[i].atomic_add(other.data[i]);
      }
    }
    __device__ AtomicCell exchange(const AtomicCell&restrict other) restrict noexcept {
      AtomicCell previous;
      for (size_t i = 0; i < Head; ++i) {
        previous.data[i] = data[i].exchange(other.data[i]);
      }
      return previous;
    }
    __device__ AtomicCell atomic_exchange(const AtomicCell&restrict other) restrict noexcept {
      AtomicCell previous;
      for (size_t i = 0; i < Head; ++i) {
        previous.data[i] = data[i].atomic_exchange(other.data[i]);
      }
      return previous;
    }
    __device__ AtomicCell shuffle(uint32_t mask, unsigned laneid) const noexcept {
      AtomicCell shuffled;
      for (size_t i = 0; i < Head; ++i) {
        shuffled.data[i] = data[i].shuffle(mask, laneid);
      }
      return shuffled;
    }
    __device__ uint32_t conflict_mask(uint32_t mask) const noexcept {
      return _conflict_mask(mask, (uintptr_t)this);
    }
    __device__ bool reduce_add(uint32_t mask, uint32_t peers) noexcept {
      return _reduce_add(mask, peers, *this);
    }
    __device__ void conflict_free_add(uint32_t mask, AtomicCell other) noexcept {
      uint32_t peers = conflict_mask(mask);
      if (other.reduce_add(mask, peers)) {
        atomic_add(other);
      }
    }
    __device__ AtomicCell<T, Tail...>& operator[](size_t i) noexcept {
      return data[i];
    }
    __device__ AtomicCell<T, Tail...> const& operator[](size_t i) const noexcept {
      return data[i];
    }
    __device__ inner_type& get() noexcept {
      return reinterpret_cast<inner_type&>(*this);
    }
    __device__ inner_type const& get() const noexcept {
      return reinterpret_cast<inner_type const&>(*this);
    }
    static __device__ AtomicCell& ref(inner_type& r) noexcept {
      return reinterpret_cast<AtomicCell&>(r);
    }
    static __device__ AtomicCell const& ref(inner_type const& r) noexcept {
      return reinterpret_cast<AtomicCell const&>(r);
    }
};
