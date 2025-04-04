#ifndef __CS_REDUCERS_H__
#define __CS_REDUCERS_H__

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

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

END_C_DECLS

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

struct cs_reduce_sum2r {
  using T = cs_data_2r;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    a.r[0] = 0;
    a.r[1] = 0;
  }

  CS_F_HOST_DEVICE void
  combine(volatile T &a, volatile const T &b) const {
    a.r[0] += b.r[0];
    a.r[1] += b.r[1];
  }
};

template<size_t stride>
struct cs_reduce_sum_n {
  using T = cs_double_n<stride>;

  CS_F_HOST_DEVICE void
  identity(T &a) const {
    for (size_t i = 0; i < stride; i++)
      a.r[i] = 0.;;
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

/*============================================================================
 * Templated function definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

#endif /* __CS_REDUCERS_H__ */
