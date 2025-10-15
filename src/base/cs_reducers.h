#ifndef CS_REDUCERS_H
#define CS_REDUCERS_H

/*============================================================================
 * Structures for reduction
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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * External library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_base.h"
#include "base/cs_math.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Structures for reduction
   ------------------------ */

// real

template<size_t stride>
struct cs_float_n {
  float r[stride];
};

template<size_t stride>
struct cs_double_n {
  double r[stride];
};

// integer

template<size_t stride>
struct cs_int_n {
  cs_lnum_t i[stride];
};

// various types

struct cs_data_1int_1float {
  cs_lnum_t i[1];
  float     r[1];
};

struct cs_data_1int_2float {
  cs_lnum_t i[1];
  float     r[2];
};

struct cs_data_2int_2float {
  cs_lnum_t i[2];
  float     r[2];
};

struct cs_data_1double_2float {
  double d[1];
  float  r[2];
};

struct cs_data_3float_3float {
  float r1[3];
  float r2[3];
};

template<size_t stride>
struct cs_data_double_n {
  double r[stride];
};

// Combined

template<size_t stride>
struct cs_data_double_int_n {
  double     r[stride];
  cs_lnum_t  i[stride];
};

/* Reduction
   --------- */

struct cs_reduce_min1r_max1r_sum2r {
  using T = cs_data_double_n<4>;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    a.r[0] =  cs_math_infinite_r;
    a.r[1] = -cs_math_infinite_r;
    a.r[2] = 0.;
    a.r[3] = 0.;
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    a.r[0] = cs::min(a.r[0], b.r[0]);
    a.r[1] = cs::max(a.r[1], b.r[1]);
    a.r[2] += b.r[2];
    a.r[3] += b.r[3];
  }
};

struct cs_reduce_sum1i_min1float {
  using T = cs_data_1int_1float;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    a.i[0] = 0;
    a.r[0] = cs_math_infinite_r;
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    a.i[0] += b.i[0];
    a.r[0] = cs::min(a.r[0], b.r[0]);
  }
};

// Max (1 real)

struct cs_reduce_max1r {

  CS_F_HOST_DEVICE void
  identity(cs_real_t &a) const {
    a = -HUGE_VAL;
  }

  CS_F_HOST_DEVICE void
  combine(volatile cs_real_t &a, volatile const cs_real_t &b) const {
    a = cs::max(a, b);
  }
};

// Min (3 reals) and Max (3 reals)

struct cs_reduce_min3float_max3float {
  using T = cs_data_3float_3float;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    a.r1[0] = cs_math_infinite_r;
    a.r1[1] = cs_math_infinite_r;
    a.r1[2] = cs_math_infinite_r;


    a.r2[0] = -cs_math_infinite_r;
    a.r2[1] = -cs_math_infinite_r;
    a.r2[2] = -cs_math_infinite_r;
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    a.r1[0] = cs::min(a.r1[0], b.r1[0]);
    a.r1[1] = cs::min(a.r1[1], b.r1[1]);
    a.r1[2] = cs::min(a.r1[2], b.r1[2]);

    a.r2[0] = cs::max(a.r2[0], b.r2[0]);
    a.r2[1] = cs::max(a.r2[1], b.r2[1]);
    a.r2[2] = cs::max(a.r2[2], b.r2[2]);
  }
};

// n_sum for double

template<size_t stride>
struct cs_reduce_sum_nr {
  using T = cs_double_n<stride>;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    for (size_t i = 0; i < stride; i++)
      a.r[i] = 0.;
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    for (size_t i = 0; i < stride; i++)
      a.r[i] += b.r[i];
  }
};

// n_sum for integer

template<size_t stride>
struct cs_reduce_sum_ni {
  using T = cs_int_n<stride>;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    for (size_t i = 0; i < stride; i++)
      a.i[i] = 0.;
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    for (size_t i = 0; i < stride; i++)
      a.i[i] += b.i[i];
  }
};

// n_min_max

template<size_t stride>
struct cs_reduce_min_max_nr {
  using T = cs_double_n<2*stride>;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    for (size_t i = 0; i < stride; i++) {
      a.r[i] = HUGE_VAL;
      a.r[stride + i] = -HUGE_VAL;
    }
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    for (size_t i = 0; i < stride; i++) {
      a.r[i] = cs::min(a.r[i], b.r[i]);
      a.r[stride + i] = cs::max(a.r[stride + i], b.r[stride + i]);
    }
  }
};

// n min_max_sum

template<size_t stride>
struct cs_reduce_min_max_sum_nr {
  using T = cs_double_n<3*stride>;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    for (size_t i = 0; i < stride; i++) {
      a.r[i] = HUGE_VAL;
      a.r[stride + i] = -HUGE_VAL;
      a.r[2*stride + i] = 0.;
    }
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    for (size_t i = 0; i < stride; i++) {
      a.r[i] = cs::min(a.r[i], b.r[i]);
      a.r[stride + i] = cs::max(a.r[stride + i], b.r[stride + i]);
      a.r[2*stride + i] += b.r[2*stride + i];
    }
  }
};

// Min (1 real), max (1 real) and sum (1 int)

struct cs_reduce_min1float_max1float_sum1int {
  using T = cs_data_1int_2float;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
      a.r[0] = HUGE_VAL;
      a.r[1] = -HUGE_VAL;
      a.i[0] = 0;
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
      a.r[0] = cs::min(a.r[0], b.r[0]);
      a.r[1] = cs::max(a.r[1], b.r[1]);
      a.i[0] += b.i[0];
  }
};

// Min (1 real), max (1 real) and sum (2 int)

struct cs_reduce_min1float_max1float_sum2int {
  using T = cs_data_2int_2float;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    a.r[0] =  cs_math_infinite_r;
    a.r[1] = -cs_math_infinite_r;
    a.i[0] = 0;
    a.i[1] = 0;
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    a.r[0] = cs::min(a.r[0], b.r[0]);
    a.r[1] = cs::max(a.r[1], b.r[1]);
    a.i[0] += b.i[0];
    a.i[1] += b.i[1];
  }
};

// Min, max (float) and sum (double)

struct cs_reduce_sum1double_min1float_max1float {
  using T = cs_data_1double_2float;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    a.d[0] = 0.0;

    a.r[0] =  cs_math_infinite_r;
    a.r[1] = -cs_math_infinite_r;
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    a.d[0] += b.d[0];

    a.r[0] = cs::min(a.r[0], b.r[0]);
    a.r[1] = cs::max(a.r[1], b.r[1]);
  }
};

template<size_t stride>
struct cs_reduce_min_max_sum_nr_with_norm {
  using T = cs_double_n<3*(stride + 1)>;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    const size_t _stride = stride + 1;
    for (size_t i = 0; i < _stride; i++) {
      a.r[i] = HUGE_VAL;
      a.r[_stride + i] = -HUGE_VAL;
      a.r[2*_stride + i] = 0.;
    }
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {

    const size_t _stride = stride + 1;

    for (size_t i = 0; i < stride; i++) {
      a.r[i] = cs::min(a.r[i], b.r[i]);
      a.r[_stride + i] = cs::max(a.r[_stride + i], b.r[_stride + i]);
      a.r[2*_stride + i] += b.r[2*_stride + i];
    }
    a.r[_stride - 1] = cs::min(a.r[_stride - 1], b.r[_stride - 1]);
    a.r[2*_stride - 1] = cs::max(a.r[2*_stride - 1], b.r[2*_stride - 1]);
    a.r[3*_stride - 1] += b.r[3*_stride - 1];
  }
};

template<size_t stride>
struct cs_reduce_min_max_weighted_sum_nr {
  using T = cs_double_n<4*stride>;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    for (size_t i = 0; i < stride; i++) {
      a.r[i] = HUGE_VAL;
      a.r[stride + i] = -HUGE_VAL;
      a.r[2*stride + i] = 0.;
      a.r[3*stride + i] = 0.;
    }
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    for (size_t i = 0; i < stride; i++) {
      a.r[i] = cs::min(a.r[i], b.r[i]);
      a.r[stride + i] = cs::max(a.r[stride + i], b.r[stride + i]);
      a.r[2*stride + i] += b.r[2*stride + i];
      a.r[3*stride + i] += b.r[3*stride + i];

    }
  }
};

template<size_t stride>
struct cs_reduce_min_max_weighted_sum_nr_with_norm {
  using T = cs_double_n<4*(stride + 1)>;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    const size_t _stride = stride + 1;
    for (size_t i = 0; i < _stride; i++) {
      a.r[i] = HUGE_VAL;
      a.r[_stride + i] = -HUGE_VAL;
      a.r[2*_stride + i] = 0.;
      a.r[3*_stride + i] = 0.;
    }
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {

    const size_t _stride = stride + 1;
    for (size_t i = 0; i < stride; i++) {
      a.r[i] = cs::min(a.r[i], b.r[i]);
      a.r[_stride + i] = cs::max(a.r[_stride + i], b.r[_stride + i]);
      a.r[2*_stride + i] += b.r[2*_stride + i];
      a.r[3*_stride + i] += b.r[3*_stride + i];
    }
    a.r[_stride - 1] = cs::min(a.r[_stride - 1], b.r[_stride - 1]);
    a.r[2*_stride - 1] = cs::max(a.r[2*_stride - 1], b.r[2*_stride - 1]);
    a.r[3*_stride - 1] += b.r[3*_stride - 1];
    a.r[4*_stride - 1] += b.r[4*_stride - 1];
  }
};

struct cs_reduce_max1float_bcast3float {
  using T = cs_float_n<4>;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    a.r[3] = -cs_math_infinite_r;

    a.r[0] = 0.;
    a.r[1] = 0.;
    a.r[2] = 0.;
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    if (a.r[3] < b.r[3]) {
      a.r[3] = b.r[3]; // a=max(a,b)

      a.r[0] = b.r[0];
      a.r[1] = b.r[1];
      a.r[2] = b.r[2];
    }
  }
};

// Strided min/max;
// The first half of the structure is used for min,
// the second for max.
// Only the first half needs to be set by the user, as combine
// handles the rest.

template<size_t stride>
struct cs_reduce_minmax_n {
  using T = cs_data_double_n<stride*2>;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    for (size_t i = 0; i < stride; i++) {
      a.r[i] = HUGE_VAL;
      a.r[stride + i] = -HUGE_VAL;
    }
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    for (size_t i = 0; i < stride; i++) {
      // Do not use stride for b, as only the first half needs to be set.
      a.r[i] = cs::min(a.r[i], b.r[i]);
      a.r[stride+i] = cs::max(a.r[stride+i], b.r[stride+i]);
    }
  }
};

// Strided min:

template<size_t stride>
struct cs_reduce_min_nr {
  using T = cs_data_double_n<stride>;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    for (size_t i = 0; i < stride; i++) {
      a.r[i] = HUGE_VAL;
    }
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    for (size_t i = 0; i < stride; i++) {
      a.r[i] = cs::min(a.r[i], b.r[i]);
    }
  }
};

// Strided min/max/loc;

template<size_t stride>
struct cs_reduce_minmaxloc_n {
  using T = cs_data_double_int_n<stride*2>;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    for (size_t i = 0; i < stride; i++) {
      a.r[i] = HUGE_VAL;
      a.r[stride + i] = -HUGE_VAL;
      a.i[i] = -1.;
      a.i[stride + i] = -1.;
    }
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    for (size_t i = 0; i < stride; i++) {
      if (a.r[i] > b.r[i]) {
        a.r[i] = b.r[i];
        a.i[i] = b.i[i];
      }
      if (a.r[stride + i] < b.r[stride + i]) {
        a.r[stride + i] = b.r[stride + i];
        a.i[stride + i] = b.i[stride + i];
      }
    }
  }
};

/*============================================================================
 * Templated function definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

#endif // __cplusplus

/*----------------------------------------------------------------------------*/

#endif /* CS_REDUCERS_H */
