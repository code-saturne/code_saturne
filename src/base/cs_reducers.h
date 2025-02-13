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

// 2 reals

struct cs_data_2r {
  double r[2];
};

struct cs_data_4r {
  double r[4];
};

struct cs_data_14r {
  cs_real_t r[14];
};

struct cs_data_3r_3r {
  cs_real_t r1[3];
  cs_real_t r2[3];
};

struct cs_data_2d {
  double d[2];
};

struct cs_data_4d {
  double d[4];
};

// 2 cs_lnum_t

struct cs_data_2i {
  cs_lnum_t i[2];
};

struct cs_data_7i {
  cs_lnum_t i[7];
};

template<size_t stride>
struct cs_double_n {
  double r[stride];

struct cs_data_3g {
  cs_gnum_t g[3];
};

struct cs_data_1d2r {
  double d[1];
  cs_real_t r[2];
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

// Max (real case)

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

// Min and max

struct cs_reduce_min1r_max1r {
  using T = cs_data_2r;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    a.r[0] =  HUGE_VAL;
    a.r[1] = -HUGE_VAL;
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    a.r[0] = cs::min(a.r[0], b.r[0]);
    a.r[1] = cs::max(a.r[1], b.r[1]);
  }
};

struct cs_reduce_min2r_max2r {
  using T = cs_data_4r;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    a.r[0] =  HUGE_VAL;
    a.r[1] = -HUGE_VAL;
    a.r[2] =  HUGE_VAL;
    a.r[3] = -HUGE_VAL;
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    a.r[0] = cs::min(a.r[0], b.r[0]);
    a.r[1] = cs::max(a.r[1], b.r[1]);
    a.r[2] = cs::min(a.r[2], b.r[2]);
    a.r[3] = cs::max(a.r[3], b.r[3]);
  }
};

struct cs_reduce_min7r_max7r {
  using T = cs_data_14r;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    for (int i = 0; i < 7; i++) {
      a.r[i] =  HUGE_VAL;
      a.r[i + 7] = -HUGE_VAL;
    }
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    for (int i = 0; i < 7; i++) {
      a.r[i] = cs::min(a.r[i], b.r[i]);
      a.r[i + 7] = cs::max(a.r[i + 7], b.r[i + 7]);
    }
  }
};

// Min (3 reals) and Max (3 reals)

struct cs_reduce_min3r_max3r {
  using T = cs_data_3r_3r;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    a.r1[0] =  cs_math_infinite_r;
    a.r1[1] =  cs_math_infinite_r;
    a.r1[2] =  cs_math_infinite_r;


    a.r2[0] = -cs_math_infinite_r;
    a.r2[1] = -cs_math_infinite_r;
    a.r2[2] = -cs_math_infinite_r;
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    a.r1[0] = cs::min(a.r1[0], b.r2[0]);
    a.r1[1] = cs::min(a.r1[1], b.r2[1]);
    a.r1[2] = cs::min(a.r1[2], b.r2[2]);

    a.r2[0] = cs::max(a.r2[0], b.r2[0]);
    a.r2[1] = cs::max(a.r2[1], b.r2[1]);
    a.r2[2] = cs::max(a.r2[2], b.r2[2]);
  }
};

// Min, max (cs_real_t), sum (double) and sum (cs_gnum_t)

struct cs_reduce_sum1d_min1r_max1r {
  using T = cs_data_1d2r;

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

// Sum and sum

struct cs_reduce_sum2i {
  using T = cs_data_2i;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    a.i[0] = 0;
    a.i[1] = 0;
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    a.i[0] += b.i[0];
    a.i[1] += b.i[1];
  }
};

struct cs_reduce_sum7i {
  using T = cs_data_7i;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    for (int j = 0; j < 7; j++)
      a.i[j] = 0;
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    for (int j = 0; j < 7; j++)
      a.i[j] += b.i[j];
  }
};

struct cs_reduce_sum3g {
  using T = cs_data_3g;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    a.g[0] = 0;
    a.g[1] = 0;
    a.g[2] = 0;
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    a.g[0] += b.g[0];
    a.g[1] += b.g[1];
    a.g[2] += b.g[2];
});

struct cs_reduce_sum2r {
  using T = cs_data_2r;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    a.r[0] = 0.;
    a.r[1] = 0.;
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    a.r[0] += b.r[0];
    a.r[1] += b.r[1];
  }
};

struct cs_reduce_sum2d {
  using T = cs_data_2d;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    a.d[0] = 0.;
    a.d[1] = 0.;
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    a.d[0] += b.d[0];
    a.d[1] += b.d[1];
  }
});

template<size_t stride>
struct cs_reduce_sum_n {
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

struct cs_reduce_sum4d {
  using T = cs_data_4d;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    a.d[0] = 0.0;
    a.d[1] = 0.0;
    a.d[2] = 0.0;
    a.d[3] = 0.0;
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    a.d[0] += b.d[0];
    a.d[1] += b.d[1];
    a.d[2] += b.d[2];
    a.d[3] += b.d[3];
  }
});

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

struct cs_reduce_max1r_bcast3r {
  using T = cs_data_4r;

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
});

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
