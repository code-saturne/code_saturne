/*============================================================================
 * Common array reduction operations.
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_mem.h"
#include "bft/bft_error.h"
#include "bft/bft_printf.h"
#include "base/cs_parall.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_array_reduce.h"
#include "base/cs_reducers.h"
#include "base/cs_dispatch.h"

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

/* Block size for superblock algorithm */

#define CS_SBLOCK_BLOCK_SIZE 60

/* Cache line multiple, in cs_real_t units */

#define CS_CL  (CS_CL_SIZE/8)

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute blocks sizes for superblock algorithm.
 *
 * parameters:
 *   n                 <-- size of array
 *   block_size        <-- block size
 *   n_sblocks         --> number of superblocks
 *   blocks_in_sblocks --> number of blocks per superblock
 *----------------------------------------------------------------------------*/

static inline void
_sbloc_sizes(cs_lnum_t   n,
             cs_lnum_t   block_size,
             cs_lnum_t  *n_sblocks,
             cs_lnum_t  *blocks_in_sblocks)
{
  cs_lnum_t n_blocks = (n + block_size - 1) / block_size;
  *n_sblocks = (n_blocks > 1) ? sqrt(n_blocks) : 1;

  cs_lnum_t n_b = block_size * *n_sblocks;
  *blocks_in_sblocks = (n + n_b - 1) / n_b;
}

/*----------------------------------------------------------------------------
 * Compute weighted sum of a 1-dimensional array
 *
 * The algorithm here is similar to that used for blas.
 *
 * parameters:
 *   n        <-- local number of elements
 *   v        <-- pointer to values (size: n)
 *   w        <-- pointer to weights (size: n)
 *
 * returns:
 *   resulting weighted sum
 *----------------------------------------------------------------------------*/

static double
_cs_real_wsum_1d(cs_lnum_t        n,
                 const cs_real_t  v[],
                 const cs_real_t  w[])
{
  double v_w_sum = 0.;

#pragma omp parallel reduction(+:v_w_sum) if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n, sizeof(cs_real_t), &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_real_t *_v = v + s_id;
    const cs_real_t *_w = w + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double s = 0.;

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
        cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);

        if (end_id > _n)
          end_id = _n;

        double c = 0.;
        for (cs_lnum_t i = start_id; i < end_id; i++)
          c += _v[i] * _w[i];

        s += c;
      }

      v_w_sum += s;
    }
  }

  return v_w_sum;
}

/*----------------------------------------------------------------------------
 * Compute weighted sum of a 1-dimensional array using an element list
 * relative to weights
 *
 * The algorithm here is similar to that used for blas.
 *
 * parameters:
 *   n        <-- local number of elements
 *   wl       <-- pointer to elements weights list
 *   v        <-- pointer to values (size: n)
 *   w        <-- pointer to weights (size: n)
 *
 * returns:
 *   resulting weighted sum
 *----------------------------------------------------------------------------*/

static double
_cs_real_wsum_1d_iw(cs_lnum_t        n,
                    const cs_lnum_t  wl[],
                    const cs_real_t  v[],
                    const cs_real_t  w[])
{
  double v_w_sum = 0.;

#pragma omp parallel reduction(+:v_w_sum) if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n, sizeof(cs_real_t), &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_lnum_t *_wl = wl + s_id;
    const cs_real_t *_v = v + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double s = 0.;

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
        cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);

        if (end_id > _n)
          end_id = _n;

        double c = 0.;
        for (cs_lnum_t i = start_id; i < end_id; i++)
          c += _v[i] * w[_wl[i]];

        s += c;
      }

      v_w_sum += s;
    }
  }

  return v_w_sum;
}
/*----------------------------------------------------------------------------
 * Compute weighted sum of a 1-dimensional array using an element list
 *
 * The algorithm here is similar to that used for blas.
 *
 * parameters:
 *   n        <-- local number of elements
 *   vl       <-- pointer to elements list
 *   v        <-- pointer to values (size: n)
 *   w        <-- pointer to weights (size: n)
 *
 * returns:
 *   resulting weighted sum
 *----------------------------------------------------------------------------*/

static double
_cs_real_wsum_1d_iv(cs_lnum_t        n,
                    const cs_lnum_t  vl[],
                    const cs_real_t  v[],
                    const cs_real_t  w[])
{
  double v_w_sum = 0.;

#pragma omp parallel reduction(+:v_w_sum) if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n, sizeof(cs_real_t), &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_lnum_t *_vl = vl + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double s = 0.;

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
        cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);

        if (end_id > _n)
          end_id = _n;

        double c = 0.;
        for (cs_lnum_t li = start_id; li < end_id; li++) {
          cs_lnum_t i = _vl[li];
          c += v[i] * w[i];
        }

        s += c;
      }

      v_w_sum += s;
    }
  }

  return v_w_sum;
}

/*----------------------------------------------------------------------------
 * Compute weighted sum and total weight (sum) of a 1-dimensional array
 *
 * The algorithm here is similar to that used for blas.
 *
 * parameters:
 *   n        <-- local number of elements
 *   v        <-- pointer to values (size: n)
 *   w        <-- pointer to weights (size: n)
 *   wsum     --> resulting weighted sum
 *   wtot     --> resulting sum of weights
 *
 *----------------------------------------------------------------------------*/

static void
_cs_real_wsum_components_1d(cs_lnum_t        n,
                            const cs_real_t  v[],
                            const cs_real_t  w[],
                            double          *wsum,
                            double          *wtot)
{
  *wsum = 0.;
  *wtot = 0.;

#pragma omp parallel if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n, sizeof(cs_real_t), &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_real_t *_v = v + s_id;
    const cs_real_t *_w = w + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    double lsum[2] = {0., 0.};

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double s[2] = {0., 0.};

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
        cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);

        if (end_id > _n)
          end_id = _n;

        double c[2] = {0., 0.};
        for (cs_lnum_t i = start_id; i < end_id; i++) {
          c[0] += _v[i] * _w[i];
          c[1] += _w[i];
        }

        s[0] += c[0];
        s[1] += c[1];
      }

      lsum[0] += s[0];
      lsum[1] += s[1];
    }

#   pragma omp critical
    {
      *wsum += lsum[0];
      *wtot += lsum[1];
    }

  }
}

/*----------------------------------------------------------------------------
 * Compute weighted sum and total weight (sum) of a 1-dimensional array
 * using an element list
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n        <-- local number of elements
 *   wl       <-- pointer to element weights list
 *   v        <-- pointer to values (size: n)
 *   w        <-- pointer to weights (size: n)
 *   wsum     --> resulting weighted sum
 *   wtot     --> resulting sum of weights
 *
 *----------------------------------------------------------------------------*/

static void
_cs_real_wsum_components_1d_iw(cs_lnum_t        n,
                               const cs_lnum_t  wl[],
                               const cs_real_t  v[],
                               const cs_real_t  w[],
                               double          *wsum,
                               double          *wtot)
{
  *wsum = 0.;
  *wtot = 0.;

#pragma omp parallel if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n, sizeof(cs_real_t), &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_lnum_t *_wl = wl + s_id;
    const cs_real_t *_v = v + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    double lsum[2] = {0., 0.};

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double s[2] = {0., 0.};

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
        cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);

        if (end_id > _n)
          end_id = _n;

        double c[2] = {0., 0.};
        for (cs_lnum_t li = start_id; li < end_id; li++) {
          cs_lnum_t i = _wl[li];
          c[0] += _v[li] * w[i];
          c[1] += w[i];
        }

        s[0] += c[0];
        s[1] += c[1];
      }

      lsum[0] += s[0];
      lsum[1] += s[1];
    }

#   pragma omp critical
    {
      *wsum += lsum[0];
      *wtot += lsum[1];
    }

  }
}

/*----------------------------------------------------------------------------
 * Compute weighted sum and total weight (sum) of a 1-dimensional array
 * using an element list
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n        <-- local number of elements
 *   vl       <-- pointer to element list
 *   v        <-- pointer to values (size: n)
 *   w        <-- pointer to weights (size: n)
 *   wsum     --> resulting weighted sum
 *   wtot     --> resulting sum of weights
 *
 *----------------------------------------------------------------------------*/

static void
_cs_real_wsum_components_1d_iv(cs_lnum_t        n,
                               const cs_lnum_t  vl[],
                               const cs_real_t  v[],
                               const cs_real_t  w[],
                               double          *wsum,
                               double          *wtot)
{
  *wsum = 0.;
  *wtot = 0.;

#pragma omp parallel if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n, sizeof(cs_real_t), &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_lnum_t *_vl = vl + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    double lsum[2] = {0., 0.};

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double s[2] = {0., 0.};

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
        cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);

        if (end_id > _n)
          end_id = _n;

        double c[2] = {0., 0.};
        for (cs_lnum_t li = start_id; li < end_id; li++) {
          cs_lnum_t i = _vl[li];
          c[0] += v[i] * w[i];
          c[1] += w[i];
        }

        s[0] += c[0];
        s[1] += c[1];
      }

      lsum[0] += s[0];
      lsum[1] += s[1];
    }

#   pragma omp critical
    {
      *wsum += lsum[0];
      *wtot += lsum[1];
    }

  }
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum) of a strided array.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n        <-- local number of elements
 *   v        <-- pointer to values (size: n)
 *   vmin     --> resulting min
 *   vmax     --> resulting max
 *   vsum     --> resulting sum
 *----------------------------------------------------------------------------*/

template <size_t stride>
static void
_cs_real_sstats(cs_dispatch_context  ctx,
                cs_lnum_t            n,
                const cs_real_t      v[],
                double              *vmin,
                double              *vmax,
                double              *vsum)
{
  struct cs_double_n<3*stride> rd;
  struct cs_reduce_min_max_sum_nr<stride> reducer;

  ctx.parallel_for_reduce(n, rd, reducer,
    [=] CS_F_HOST_DEVICE (cs_lnum_t i, cs_double_n<3*stride> &res) {
    for (size_t k = 0; k < stride; k++) {
      res.r[k] = v[stride*i + k];
      res.r[stride + k] = v[stride*i + k];
      res.r[2*stride + k] = v[stride*i + k];
    }
  });

  ctx.wait();

  for (size_t k = 0; k < stride; k++) {
    vmin[k] = rd.r[k];
    vmax[k] = rd.r[stride + k];
    vsum[k] = rd.r[2*stride + k];
  }
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum) of a strided array.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n        <-- local number of elements
 *   v        <-- pointer to values (size: n)
 *   vmin     --> resulting min
 *   vmax     --> resulting max
 *   vsum     --> resulting sum
 *----------------------------------------------------------------------------*/

template <size_t stride>
static void
_cs_real_sstats_iv(cs_dispatch_context  ctx,
                   cs_lnum_t            n,
                   const cs_lnum_t     *vl,
                   const cs_real_t      v[],
                   double              *vmin,
                   double              *vmax,
                   double              *vsum)
{
  struct cs_double_n<3*stride> rd;
  struct cs_reduce_min_max_sum_nr<stride> reducer;

  ctx.parallel_for_reduce(n, rd, reducer,
    [=] CS_F_HOST_DEVICE (cs_lnum_t i, cs_double_n<3*stride> &res) {
    for (size_t k = 0; k < stride; k++) {
      res.r[k] = v[stride*vl[i] + k];
      res.r[stride + k] = v[stride*vl[i] + k];
      res.r[2*stride + k] = v[stride*vl[i] + k];
    }
  });

  ctx.wait();

  for (size_t k = 0; k < stride; k++) {
    vmin[k] = rd.r[k];
    vmax[k] = rd.r[stride + k];
    vsum[k] = rd.r[2*stride + k];
  }
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, several sums and weighted sums)
 * of a 1-dimensional array.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n        <-- local number of elements
 *   v        <-- pointer to values (size: n)
 *   w        <-- pointer to weights (size: n)
 *   vmin     --> resulting min
 *   vmax     --> resulting max
 *   vsum     --> resulting sum
 *   wsum     --> resulting weighted sum
 *   asum     --> resulting sum of absolute values
 *   ssum     --> resulting weighted sum array
 *   wssum    --> resulting weighted sum of squared values
 *----------------------------------------------------------------------------*/

static void
_cs_real_norms_1d(cs_lnum_t         n,
                  const cs_real_t   v[],
                  const cs_real_t   w[],
                  double           *vmin,
                  double           *vmax,
                  double           *vsum,
                  double           *wsum,
                  double           *asum,
                  double           *ssum,
                  double           *wssum)
{
  *vmin = HUGE_VAL;
  *vmax = -HUGE_VAL;
  *vsum = 0.;
  *wsum = 0.;
  *asum = 0.;
  *ssum = 0.;
  *wssum = 0.;

# pragma omp parallel if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n, sizeof(cs_real_t), &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_real_t *_v = v + s_id;
    const cs_real_t *_w = w + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    cs_real_t lmin = HUGE_VAL;
    cs_real_t lmax = -HUGE_VAL;

    double lsum[5] = {0., 0., 0., 0., 0.};

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double s[5] = {0., 0., 0., 0., 0.};

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {

        cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        if (end_id > _n)
          end_id = _n;
        double c[5] = {0., 0., 0., 0., 0.};
        for (cs_lnum_t i = start_id; i < end_id; i++) {
          const cs_real_t  val = _v[i], val2 = val*val;
          c[0] += val;
          c[1] += val*_w[i];
          c[2] += fabs(val);
          c[3] += val2;
          c[4] += val2*_w[i];
          if (val < lmin)
            lmin = val;
          if (val > lmax)
            lmax = val;
        }
        s[0] += c[0];
        s[1] += c[1];
        s[2] += c[2];
        s[3] += c[3];
        s[4] += c[4];

      } /* Loop on blocks */

      lsum[0] += s[0];
      lsum[1] += s[1];
      lsum[2] += s[2];
      lsum[3] += s[3];
      lsum[4] += s[4];

    } /* Loop on super-block */

#   pragma omp critical
    {
      if (lmin < *vmin)
        *vmin = lmin;
      if (lmax > *vmax)
        *vmax = lmax;
      *vsum += lsum[0];
      *wsum += lsum[1];
      *asum += lsum[2];
      *ssum += lsum[3];
      *wssum += lsum[4];

    }

  }
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, several sums and weighted sums)
 * of a 1-dimensional array.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n        <-- local number of elements
 *   wl       <-- pointer to element weights list
 *   v        <-- pointer to values (size: n)
 *   w        <-- pointer to weights (size: n)
 *   vmin     --> resulting min
 *   vmax     --> resulting max
 *   vsum     --> resulting sum
 *   wsum     --> resulting weighted sum
 *   asum     --> resulting sum of absolute values
 *   ssum     --> resulting weighted sum array
 *   wssum    --> resulting weighted sum of squared values
 *----------------------------------------------------------------------------*/

static void
_cs_real_norms_1d_iw(cs_lnum_t         n,
                     const cs_lnum_t   wl[],
                     const cs_real_t   v[],
                     const cs_real_t   w[],
                     double           *vmin,
                     double           *vmax,
                     double           *vsum,
                     double           *wsum,
                     double           *asum,
                     double           *ssum,
                     double           *wssum)
{
  *vmin = HUGE_VAL;
  *vmax = -HUGE_VAL;
  *vsum = 0.;
  *wsum = 0.;
  *asum = 0.;
  *ssum = 0.;
  *wssum = 0.;

# pragma omp parallel if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n, sizeof(cs_real_t), &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_lnum_t *_wl = wl + s_id;
    const cs_real_t *_v = v + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    cs_real_t lmin = HUGE_VAL;
    cs_real_t lmax = -HUGE_VAL;

    double lsum[5] = {0., 0., 0., 0., 0.};

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double s[5] = {0., 0., 0., 0., 0.};

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {

        cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        if (end_id > _n)
          end_id = _n;
        double c[5] = {0., 0., 0., 0., 0.};
        for (cs_lnum_t i = start_id; i < end_id; i++) {
          const cs_real_t  val = _v[i], val2 = val*val;
          c[0] += val;
          c[1] += val*w[_wl[i]];
          c[2] += fabs(val);
          c[3] += val2;
          c[4] += val2*w[_wl[i]];
          if (val < lmin)
            lmin = val;
          if (val > lmax)
            lmax = val;
        }
        s[0] += c[0];
        s[1] += c[1];
        s[2] += c[2];
        s[3] += c[3];
        s[4] += c[4];

      } /* Loop on blocks */

      lsum[0] += s[0];
      lsum[1] += s[1];
      lsum[2] += s[2];
      lsum[3] += s[3];
      lsum[4] += s[4];

    } /* Loop on super-block */

#   pragma omp critical
    {
      if (lmin < *vmin)
        *vmin = lmin;
      if (lmax > *vmax)
        *vmax = lmax;
      *vsum += lsum[0];
      *wsum += lsum[1];
      *asum += lsum[2];
      *ssum += lsum[3];
      *wssum += lsum[4];

    }

  }
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, several sums and weighted sums)
 * of a 1-dimensional array.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n        <-- local number of elements
 *   vl       <-- pointer to element list
 *   v        <-- pointer to values (size: n)
 *   w        <-- pointer to weights (size: n)
 *   vmin     --> resulting min
 *   vmax     --> resulting max
 *   vsum     --> resulting sum
 *   wsum     --> resulting weighted sum
 *   asum     --> resulting sum of absolute values
 *   ssum     --> resulting weighted sum array
 *   wssum    --> resulting weighted sum of squared values
 *----------------------------------------------------------------------------*/

static void
_cs_real_norms_1d_iv(cs_lnum_t         n,
                     const cs_lnum_t   vl[],
                     const cs_real_t   v[],
                     const cs_real_t   w[],
                     double           *vmin,
                     double           *vmax,
                     double           *vsum,
                     double           *wsum,
                     double           *asum,
                     double           *ssum,
                     double           *wssum)
{
  *vmin = HUGE_VAL;
  *vmax = -HUGE_VAL;
  *vsum = 0.;
  *wsum = 0.;
  *asum = 0.;
  *ssum = 0.;
  *wssum = 0.;

# pragma omp parallel if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n, sizeof(cs_real_t), &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_lnum_t *_vl = vl + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    cs_real_t lmin = HUGE_VAL;
    cs_real_t lmax = -HUGE_VAL;

    double lsum[5] = {0., 0., 0., 0., 0.};

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double s[5] = {0., 0., 0., 0., 0.};

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {

        cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        if (end_id > _n)
          end_id = _n;
        double c[5] = {0., 0., 0., 0., 0.};
        for (cs_lnum_t li = start_id; li < end_id; li++) {
          cs_lnum_t i = _vl[li];
          const cs_real_t  val = v[i], val2 = val*val;
          c[0] += val;
          c[1] += val*w[i];
          c[2] += fabs(val);
          c[3] += val2;
          c[4] += val2*w[i];
          if (val < lmin)
            lmin = val;
          if (val > lmax)
            lmax = val;
        }
        s[0] += c[0];
        s[1] += c[1];
        s[2] += c[2];
        s[3] += c[3];
        s[4] += c[4];

      } /* Loop on blocks */

      lsum[0] += s[0];
      lsum[1] += s[1];
      lsum[2] += s[2];
      lsum[3] += s[3];
      lsum[4] += s[4];

    } /* Loop on super-block */

#   pragma omp critical
    {
      if (lmin < *vmin)
        *vmin = lmin;
      if (lmax > *vmax)
        *vmax = lmax;
      *vsum += lsum[0];
      *wsum += lsum[1];
      *asum += lsum[2];
      *ssum += lsum[3];
      *wssum += lsum[4];

    }

  }
}

/*----------------------------------------------------------------------------
 * Compute simple local weighted norms (l1, l2) of a 1-dimensional array.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n_src_elts  <-- number of source elements
 *   src2v_idx   <-- index array between source elements and local elements
 *   src2v_ids   <-- list of ids of local elements
 *   v           <-- pointer to array values
 *   w           <-- pointer to weights (size = src2v_idx[n_src_elts])
 *   vsum        --> (weighted) sum array (size: dim, or 4 if dim = 3)
 *   asum        --> (weighted) sum of absolute values (same size as vsum)
 *   ssum        --> (weighted) sum of squared values (same size as vsum)
 *----------------------------------------------------------------------------*/

static void
_cs_real_scatter_norms_1d(cs_lnum_t          n_src_elts,
                          const cs_lnum_t   *src2v_idx,
                          const cs_lnum_t   *src2v_ids,
                          const cs_real_t    v[],
                          const cs_real_t    w[],
                          double             vsum[],
                          double             asum[],
                          double             ssum[])
{
  *vsum = 0.;
  *asum = 0.;
  *ssum = 0.;

# pragma omp parallel if (n_src_elts > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n_src_elts, sizeof(cs_real_t), &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    double lvsum = 0., lasum = 0., lssum = 0.;

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double s[3] = {0., 0., 0.};

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {

        const cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        if (end_id > _n)
          end_id = _n;

        double c[3] = {0., 0., 0.};

        /* Loop on source elements */
        for (cs_lnum_t id = s_id + start_id; id < s_id + end_id; id++) {
          for (cs_lnum_t j = src2v_idx[id]; j < src2v_idx[id+1]; j++) {

            const cs_real_t  val = v[src2v_ids[j]];
            const cs_real_t  weight = w[j];

            c[0] += val*weight;
            c[1] += fabs(val)*weight;
            c[2] += val*val*weight;

          }
        }

        s[0] += c[0];
        s[1] += c[1];
        s[2] += c[2];

      } /* Loop on blocks */

      lvsum += s[0];
      lasum += s[1];
      lssum += s[2];

    } /* Loop on super-block */

#   pragma omp critical
    {
      *vsum += lvsum;
      *asum += lasum;
      *ssum += lssum;
    }

  } /* Block openMP */
}

/*----------------------------------------------------------------------------
 * Compute simple local weighted norms (l1, l2) of a 1-dimensional array.
 *
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n_src_elts  <-- number of source elements
 *   src2v_idx   <-- index array between source elements and local elements
 *   src2v_ids   <-- list of ids of local elements
 *   filter_list <-- list of source elements on which values are defined
 *   v           <-- pointer to array values
 *   w           <-- pointer to weights (size = src2v_idx[n_src_elts])
 *   vsum        --> (weighted) sum array (size: dim, or 4 if dim = 3)
 *   asum        --> (weighted) sum of absolute values (same size as vsum)
 *   ssum        --> (weighted) sum of squared values (same size as vsum)
 *----------------------------------------------------------------------------*/

static void
_cs_real_scatter_norms_1d_filtered(cs_lnum_t          n_src_elts,
                                   const cs_lnum_t   *src2v_idx,
                                   const cs_lnum_t   *src2v_ids,
                                   const cs_lnum_t   *filter_list,
                                   const cs_real_t    v[],
                                   const cs_real_t    w[],
                                   double             vsum[],
                                   double             asum[],
                                   double             ssum[])
{
  *vsum = 0.;
  *asum = 0.;
  *ssum = 0.;

# pragma omp parallel if (n_src_elts > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n_src_elts, sizeof(cs_real_t), &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_lnum_t *_fl = filter_list + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    double lvsum = 0., lasum = 0., lssum = 0.;

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double s[3] = {0., 0., 0.};

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {

        const cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        if (end_id > _n)
          end_id = _n;

        double c[3] = {0., 0., 0.};

        /* Loop on source elements */
        for (cs_lnum_t lid = start_id; lid < end_id; lid++) {
          const cs_lnum_t  id = _fl[lid];
          for (cs_lnum_t j = src2v_idx[id]; j < src2v_idx[id+1]; j++) {
            const cs_real_t  val = v[src2v_ids[j]];
            const cs_real_t  weight = w[j];

            c[0] += val*weight;
            c[1] += fabs(val)*weight;
            c[2] += val*val*weight;

          }
        }

        s[0] += c[0];
        s[1] += c[1];
        s[2] += c[2];

      } /* Loop on blocks */

      lvsum += s[0];
      lasum += s[1];
      lssum += s[2];

    } /* Loop on super-block */

#   pragma omp critical
    {
      *vsum += lvsum;
      *asum += lasum;
      *ssum += lssum;
    }

  } /* Block openMP */
}

/*----------------------------------------------------------------------------
 * Compute the min./max. of a 3-dimensional array.
 *
 * The algorithm here is similar to that used for blas.
 *
 * parameters:
 *   n        <-- local number of elements
 *   v        <-- pointer to values (size: n)
 *   vmin     --> resulting min array (size: dim, or 4 if dim = 3)
 *   vmax     --> resulting max array (same size as vmin)
 *----------------------------------------------------------------------------*/

static void
_cs_real_minmax_3d(cs_lnum_t           n,
                   const cs_real_3_t   v[],
                   cs_real_t          *vmin,
                   cs_real_t          *vmax)
{
  for (int j = 0; j < 4; j++) {
    vmin[j] = HUGE_VAL;
    vmax[j] = -HUGE_VAL;
  }

# pragma omp parallel if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n, sizeof(cs_real_t), &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_real_3_t *_v = v + s_id;

    cs_real_t lmin[4], lmax[4];
    for (int j = 0; j < 4; j++) {
      lmin[j] = HUGE_VAL;
      lmax[j] = -HUGE_VAL;
    }

    for (cs_lnum_t i = 0; i < _n; i++) {

      double v_norm2 = 0.;
      for (int j = 0; j < 3; j++) {
        const double  val = _v[i][j];
        v_norm2 += val*val;
        if (val < lmin[j])
          lmin[j] = val;
        if (val > lmax[j])
          lmax[j] = val;
      }
      const double v_norm = sqrt(v_norm2);

      if (v_norm < lmin[3])
        lmin[3] = v_norm;
      if (v_norm > lmax[3])
        lmax[3] = v_norm;

    }

#   pragma omp critical
    {
      for (int j = 0; j < 4; j++) {
        if (lmin[j] < vmin[j])
          vmin[j] = lmin[j];
        if (lmax[j] > vmax[j])
          vmax[j] = lmax[j];
      }
    }

  } /* openMP block */
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum, and vector norm) of a
 * strided array's components and norm.
 *
 * The array is interleaved.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   ctx      <-- cs_dispatch_context struct from higer level
 *   n        <-- local number of elements
 *   v        <-- pointer to field values (size: n*stride)
 *   vmin     --> resulting min array (size: stride + 1)
 *   vmax     --> resulting max array (size: stride + 1)
 *   vsum     --> resulting sum array (size: stride + 1)
 *----------------------------------------------------------------------------*/

template <size_t stride>
static void
_cs_real_sstats_with_norm(cs_dispatch_context  ctx,
                          cs_lnum_t            n,
                          const cs_real_t      v[],
                          double               vmin[],
                          double               vmax[],
                          double               vsum[])
{
  const size_t _stride = stride + 1;
  struct cs_double_n<3*_stride> rd;
  struct cs_reduce_min_max_sum_nr_with_norm<stride> reducer;

  ctx.parallel_for_reduce(n, rd, reducer,
    [=] CS_F_HOST_DEVICE (cs_lnum_t i, cs_double_n<3*_stride> &res) {
    for (size_t k = 0; k < stride; k++) {
      res.r[k] = v[stride*i + k];
      res.r[_stride + k] = v[stride*i + k];
      res.r[2*_stride + k] = v[stride*i + k];
    }
  });

  ctx.wait();

  for (size_t k = 0; k < _stride; k++) {
    vmin[k] = rd.r[k];
    vmax[k] = rd.r[_stride + k];
    vsum[k] = rd.r[2*_stride + k];
  }
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum) of a subset of a strided
 * array's components and norm.
 *
 * The array is interleaved.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   ctx      <-- cs_dispatch_context struct from higer level
 *   n        <-- local number of elements
 *   vl       <-- pointer to element list
 *   v        <-- pointer to field values (size: n*stride)
 *   vmin     --> resulting min array (size: stride + 1)
 *   vmax     --> resulting max array (size: stride + 1)
 *   vsum     --> resulting sum array (size: stride + 1)
 *----------------------------------------------------------------------------*/

template <size_t stride>
static void
_cs_real_sstats_with_norm_iv(cs_dispatch_context  ctx,
                             cs_lnum_t            n,
                             const cs_lnum_t      vl[],
                             const cs_real_t      v[],
                             double              *vmin,
                             double              *vmax,
                             double              *vsum)
{
  const size_t _stride = stride + 1;
  struct cs_double_n<3*_stride> rd;
  struct cs_reduce_min_max_sum_nr_with_norm<stride> reducer;

  ctx.parallel_for_reduce(n, rd, reducer,
    [=] CS_F_HOST_DEVICE (cs_lnum_t i, cs_double_n<3*_stride> &res) {
    for (size_t k = 0; k < stride; k++) {
      res.r[k] = v[stride*vl[i] + k];
      res.r[_stride + k] = v[stride*vl[i] + k];
      res.r[2*_stride + k] = v[stride*vl[i] + k];
    }
  });

  ctx.wait();

  for (size_t k = 0; k < _stride; k++) {
    vmin[k] = rd.r[k];
    vmax[k] = rd.r[_stride + k];
    vsum[k] = rd.r[2*_stride + k];
  }
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum, weighted sum) of a strided
 * array's components.
 *
 * The array is interleaved.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   ctx      <-- cs_dispatch_context struct from higer level
 *   n        <-- local number of elements
 *   v        <-- pointer to field values (size: n*stride)
 *   w        <-- pointer to weight values (size: n)
 *   vmin     --> resulting min array (size: stride)
 *   vmax     --> resulting max array (size: stride)
 *   vsum     --> resulting sum array (size: stride)
 *   wsum     --> resulting weighted sum array (size: stride)
 *----------------------------------------------------------------------------*/

template <size_t stride>
static void
_cs_real_sstats_weighted(cs_dispatch_context  ctx,
                         cs_lnum_t            n,
                         const cs_real_t      v[],
                         const cs_real_t      w[],
                         double               vmin[],
                         double               vmax[],
                         double               vsum[],
                         double               wsum[])
{
  struct cs_double_n<4*stride> rd;
  struct cs_reduce_min_max_weighted_sum_nr<stride> reducer;

  ctx.parallel_for_reduce(n, rd, reducer,
    [=] CS_F_HOST_DEVICE (cs_lnum_t i, cs_double_n<4*stride> &res) {
    for (size_t k = 0; k < stride; k++) {
      res.r[k] = v[stride*i + k];
      res.r[stride + k] = v[stride*i + k];
      res.r[2*stride + k] = v[stride*i + k];
      res.r[3*stride + k] = w[i]*v[stride*i + k];
    }
  });

  ctx.wait();

  for (size_t k = 0; k < stride; k++) {
    vmin[k] = rd.r[k];
    vmax[k] = rd.r[stride + k];
    vsum[k] = rd.r[2*stride + k];
    wsum[k] = rd.r[3*stride + k];
  }
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum, weighted sum) of a strided
 * array's components.
 *
 * The array is interleaved.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   ctx      <-- cs_dispatch_context struct from higer level
 *   n        <-- local number of elements
 *   vl       <-- pointer to element list
 *   v        <-- pointer to field values (size: n*stride)
 *   w        <-- pointer to weight values (size: n)
 *   vmin     --> resulting min array (size: stride)
 *   vmax     --> resulting max array (size: stride)
 *   vsum     --> resulting sum array (size: stride)
 *   wsum     --> resulting weighted sum array (size: stride)
 *----------------------------------------------------------------------------*/

template <size_t stride>
static void
_cs_real_sstats_weighted_iv(cs_dispatch_context  ctx,
                            cs_lnum_t            n,
                            const cs_lnum_t      vl[],
                            const cs_real_t      v[],
                            const cs_real_t      w[],
                            double               vmin[],
                            double               vmax[],
                            double               vsum[],
                            double               wsum[])
{
  struct cs_double_n<4*stride> rd;
  struct cs_reduce_min_max_weighted_sum_nr<stride> reducer;

  ctx.parallel_for_reduce(n, rd, reducer,
    [=] CS_F_HOST_DEVICE (cs_lnum_t i, cs_double_n<4*stride> &res) {
    for (size_t k = 0; k < stride; k++) {
      res.r[k] = v[stride*vl[i] + k];
      res.r[stride + k] = v[stride*vl[i] + k];
      res.r[2*stride + k] = v[stride*vl[i] + k];
      res.r[3*stride + k] = w[vl[i]]*v[stride*vl[i] + k];
    }
  });

  ctx.wait();

  for (size_t k = 0; k < stride; k++) {
    vmin[k] = rd.r[k];
    vmax[k] = rd.r[stride + k];
    vsum[k] = rd.r[2*stride + k];
    wsum[k] = rd.r[3*stride + k];
  }
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum, weighted sum) of a strided
 * array's components.
 *
 * The array is interleaved.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   ctx      <-- cs_dispatch_context struct from higer level
 *   n        <-- local number of elements
 *   wl       <-- pointer to element list
 *   v        <-- pointer to field values (size: n*stride)
 *   w        <-- pointer to weight values (size: n)
 *   vmin     --> resulting min array (size: stride)
 *   vmax     --> resulting max array (size: stride)
 *   vsum     --> resulting sum array (size: stride)
 *   wsum     --> resulting weighted sum array (size: stride)
 *----------------------------------------------------------------------------*/

template <size_t stride>
static void
_cs_real_sstats_weighted_iw(cs_dispatch_context  ctx,
                            cs_lnum_t            n,
                            const cs_lnum_t      wl[],
                            const cs_real_t      v[],
                            const cs_real_t      w[],
                            double               vmin[],
                            double               vmax[],
                            double               vsum[],
                            double               wsum[])
{
  struct cs_double_n<4*stride> rd;
  struct cs_reduce_min_max_weighted_sum_nr<stride> reducer;

  ctx.parallel_for_reduce(n, rd, reducer,
    [=] CS_F_HOST_DEVICE (cs_lnum_t i, cs_double_n<4*stride> &res) {
    for (size_t k = 0; k < stride; k++) {
      res.r[k] = v[stride*i + k];
      res.r[stride + k] = v[stride*i + k];
      res.r[2*stride + k] = v[stride*i + k];
      res.r[3*stride + k] = w[wl[i]]*v[stride*i + k];
    }
  });

  ctx.wait();

  for (size_t k = 0; k < stride; k++) {
    vmin[k] = rd.r[k];
    vmax[k] = rd.r[stride + k];
    vsum[k] = rd.r[2*stride + k];
    wsum[k] = rd.r[3*stride + k];
  }
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum, weighted sum) of a strided
 * array's components and norm.
 *
 * The array is interleaved.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   ctx      <-- cs_dispatch_context struct from higer level
 *   n        <-- local number of elements
 *   v        <-- pointer to field values (size: n*stride)
 *   w        <-- pointer to weight values (size: n)
 *   vmin     --> resulting min array (size: stride + 1)
 *   vmax     --> resulting max array (size: stride + 1)
 *   vsum     --> resulting sum array (size: stride + 1)
 *   wsum     --> resulting weighted sum array (size: stride + 1)
 *----------------------------------------------------------------------------*/

template <size_t stride>
static void
_cs_real_sstats_weighted_with_norm(cs_dispatch_context  ctx,
                                   cs_lnum_t            n,
                                   const cs_real_t      v[],
                                   const cs_real_t      w[],
                                   double               vmin[],
                                   double               vmax[],
                                   double               vsum[],
                                   double               wsum[])
{
  /* stride + 1 due to computing the vector norm */
  const size_t _stride = stride + 1;
  struct cs_double_n<4*_stride> rd;
  struct cs_reduce_min_max_weighted_sum_nr_with_norm<stride> reducer;

  ctx.parallel_for_reduce(n, rd, reducer,
    [=] CS_F_HOST_DEVICE (cs_lnum_t i, cs_double_n<4*_stride> &res) {
    double norm = 0.;
    for (size_t k = 0; k < stride; k++) {
      res.r[k] = v[stride*i + k];
      res.r[_stride + k] = v[stride*i + k];
      res.r[2*_stride + k] = v[stride*i + k];
      res.r[3*_stride + k] = w[i]*v[stride*i + k];
      norm += v[stride*i + k]*v[stride*i + k];
    }
    res.r[_stride - 1] = sqrt(norm);
    res.r[2*_stride - 1] = sqrt(norm);
    res.r[3*_stride - 1] = sqrt(norm);
    res.r[4*_stride - 1] = w[i]*sqrt(norm);
  });

  ctx.wait();

  for (size_t k = 0; k < _stride; k++) {
    vmin[k] = rd.r[k];
    vmax[k] = rd.r[_stride + k];
    vsum[k] = rd.r[2*_stride + k];
    wsum[k] = rd.r[3*_stride + k];
  }
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum, weighted sum) of a strided
 * array's components and norm.
 *
 * The array is interleaved.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   ctx      <-- cs_dispatch_context struct from higer level
 *   n        <-- local number of elements
 *   vl       <-- pointer to element list
 *   v        <-- pointer to field values (size: n*stride)
 *   w        <-- pointer to weight values (size: n)
 *   vmin     --> resulting min array (size: stride + 1)
 *   vmax     --> resulting max array (size: stride + 1)
 *   vsum     --> resulting sum array (size: stride + 1)
 *   wsum     --> resulting weighted sum array (size: stride + 1)
 *----------------------------------------------------------------------------*/

template <size_t stride>
static void
_cs_real_sstats_weighted_with_norm_iv(cs_dispatch_context  ctx,
                                      cs_lnum_t            n,
                                      const cs_lnum_t      vl[],
                                      const cs_real_t      v[],
                                      const cs_real_t      w[],
                                      double               vmin[],
                                      double               vmax[],
                                      double               vsum[],
                                      double               wsum[])
{
  /* stride + 1 due to computing the vector norm */
  const size_t _stride = stride + 1;
  struct cs_double_n<4*_stride> rd;
  struct cs_reduce_min_max_weighted_sum_nr_with_norm<stride> reducer;

  ctx.parallel_for_reduce(n, rd, reducer,
    [=] CS_F_HOST_DEVICE (cs_lnum_t i, cs_double_n<4*_stride> &res) {
    double norm = 0.;
    for (size_t k = 0; k < stride; k++) {
      res.r[k] = v[stride*vl[i] + k];
      res.r[_stride + k] = v[stride*vl[i] + k];
      res.r[2*_stride + k] = v[stride*vl[i] + k];
      res.r[3*_stride + k] = w[vl[i]]*v[stride*vl[i] + k];
      norm += v[stride*vl[i] + k]*v[stride*vl[i] + k];
    }
    res.r[_stride - 1] = sqrt(norm);
    res.r[2*_stride - 1] = sqrt(norm);
    res.r[3*_stride - 1] = sqrt(norm);
    res.r[4*_stride - 1] = w[vl[i]]*sqrt(norm);
  });

  ctx.wait();

  for (size_t k = 0; k < _stride; k++) {
    vmin[k] = rd.r[k];
    vmax[k] = rd.r[_stride + k];
    vsum[k] = rd.r[2*_stride + k];
    wsum[k] = rd.r[3*_stride + k];
  }
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum, weighted sum) of a strided
 * array's components and norm.
 *
 * The array is interleaved.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   ctx      <-- cs_dispatch_context struct from higer level
 *   n        <-- local number of elements
 *   wl       <-- pointer to element list
 *   v        <-- pointer to field values (size: n*stride)
 *   w        <-- pointer to weight values (size: n)
 *   vmin     --> resulting min array (size: stride + 1)
 *   vmax     --> resulting max array (size: stride + 1)
 *   vsum     --> resulting sum array (size: stride + 1)
 *   wsum     --> resulting weighted sum array (size: stride + 1)
 *----------------------------------------------------------------------------*/

template <size_t stride>
static void
_cs_real_sstats_weighted_with_norm_iw(cs_dispatch_context  ctx,
                                      cs_lnum_t            n,
                                      const cs_lnum_t      wl[],
                                      const cs_real_t      v[],
                                      const cs_real_t      w[],
                                      double               vmin[],
                                      double               vmax[],
                                      double               vsum[],
                                      double               wsum[])
{
  /* stride + 1 due to computing the vector norm */
  const size_t _stride = stride + 1;
  struct cs_double_n<4*_stride> rd;
  struct cs_reduce_min_max_weighted_sum_nr_with_norm<stride> reducer;

  ctx.parallel_for_reduce(n, rd, reducer,
    [=] CS_F_HOST_DEVICE (cs_lnum_t i, cs_double_n<4*_stride> &res) {
    double norm = 0.;
    for (size_t k = 0; k < stride; k++) {
      res.r[k] = v[stride*i + k];
      res.r[_stride + k] = v[stride*i + k];
      res.r[2*_stride + k] = v[stride*i + k];
      res.r[3*_stride + k] = w[wl[i]]*v[stride*i + k];
      norm += v[stride*i + k]*v[stride*i + k];
    }
    res.r[_stride - 1] = sqrt(norm);
    res.r[2*_stride - 1] = sqrt(norm);
    res.r[3*_stride - 1] = sqrt(norm);
    res.r[4*_stride - 1] = w[wl[i]]*sqrt(norm);
  });

  ctx.wait();

  for (size_t k = 0; k < _stride; k++) {
    vmin[k] = rd.r[k];
    vmax[k] = rd.r[_stride + k];
    vsum[k] = rd.r[2*_stride + k];
    wsum[k] = rd.r[3*_stride + k];
  }
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, several sums and weighted sums)
 * of a 1-dimensional array.
 *
 * The array is interleaved.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n        <-- local number of elements
 *   v        <-- pointer to field values (size: n*3)
 *   w        <-- pointer to weight values (size: n)
 *   vmin     --> resulting min array (size: 4)
 *   vmax     --> resulting max array (size: 4)
 *   vsum     --> resulting sum array (size: 4)
 *   wsum     --> resulting weighted sum array (size: 4)
 *   asum     --> resulting sum of absolute values (size: 4)
 *   ssum     --> resulting weighted sum array (size: 4)
 *   wssum    --> resulting weighted sum of squared values (size: 4)
 *----------------------------------------------------------------------------*/

static void
_cs_real_norms_3d(cs_lnum_t           n,
                  const cs_real_3_t   v[],
                  const cs_real_t     w[],
                  double              vmin[4],
                  double              vmax[4],
                  double              vsum[4],
                  double              wsum[4],
                  double              asum[4],
                  double              ssum[4],
                  double              wssum[4])
{
  for (cs_lnum_t j = 0; j < 4; j++) {
    vmin[j] = HUGE_VAL;
    vmax[j] = -HUGE_VAL;
    vsum[j] = 0.;
    wsum[j] = 0.;
    asum[j] = 0.;
    ssum[j] = 0.;
    wssum[j] = 0.;
  }

# pragma omp parallel if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n, sizeof(cs_real_t), &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_real_3_t *_v = v + s_id;
    const cs_real_t *_w = w + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    cs_real_t lmin[4], lmax[4];
    for (int j = 0; j < 4; j++) {
      lmin[j] = HUGE_VAL;
      lmax[j] = -HUGE_VAL;
    }

    double lsum[20];
    for (int j = 0; j < 20; j++) lsum[j] = 0.;

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double s[20];
      for (cs_lnum_t j = 0; j < 20; j++) s[j] = 0.;

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {

        cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        if (end_id > _n)
          end_id = _n;

        double c[20];
        for (cs_lnum_t j = 0; j < 20; j++) c[j] = 0.;

        for (cs_lnum_t i = start_id; i < end_id; i++) {
          double v_norm2 = 0.;
          for (cs_lnum_t k = 0; k < 3; k++) {
            const cs_real_t  val = _v[i][k], val2 = val*val;
            c[k]   += val;
            c[k+4] += val*_w[i];
            c[k+8] += fabs(val);
            c[k+12] += val2;
            c[k+16] += val2*_w[i];
            v_norm2 += val2;
            if (val < lmin[k])
              lmin[k] = val;
            if (val > lmax[k])
              lmax[k] = val;
          }

          const double v_norm = sqrt(v_norm2);
          c[3] += v_norm;
          c[7] += v_norm*_w[i];
          c[11] += v_norm;
          c[15] += v_norm2;
          c[19] += v_norm2*_w[i];
          if (v_norm < lmin[3])
            lmin[3] = v_norm;
          if (v_norm > lmax[3])
            lmax[3] = v_norm;
        }

        for (cs_lnum_t j = 0; j < 20; j++)
          s[j] += c[j];

      } /* Loop on blocks */

      for (cs_lnum_t j = 0; j < 20; j++)
        lsum[j] += s[j];

    } /* Loop on super-block */

#   pragma omp critical
    {
      for (cs_lnum_t j = 0; j < 4; j++) {
        if (lmin[j] < vmin[j])
          vmin[j] = lmin[j];
        if (lmax[j] > vmax[j])
          vmax[j] = lmax[j];
        vsum[j] += lsum[j];
        wsum[j] += lsum[4+j];
        asum[j] += lsum[8+j];
        ssum[j] += lsum[12+j];
        wssum[j] += lsum[16+j];
      }
    }

  }
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, several sums and weighted sums)
 * of a 1-dimensional array.
 *
 * The array is interleaved.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n        <-- local number of elements
 *   wl       <-- pointer to weights list
 *   v        <-- pointer to field values (size: n*3)
 *   w        <-- pointer to weight values (size: n)
 *   vmin     --> resulting min array (size: 4)
 *   vmax     --> resulting max array (size: 4)
 *   vsum     --> resulting sum array (size: 4)
 *   wsum     --> resulting weighted sum array (size: 4)
 *   asum     --> resulting sum of absolute values (size: 4)
 *   ssum     --> resulting weighted sum array (size: 4)
 *   wssum    --> resulting weighted sum of squared values (size: 4)
 *----------------------------------------------------------------------------*/

static void
_cs_real_norms_3d_iw(cs_lnum_t           n,
                     const cs_lnum_t     wl[],
                     const cs_real_3_t   v[],
                     const cs_real_t     w[],
                     double              vmin[4],
                     double              vmax[4],
                     double              vsum[4],
                     double              wsum[4],
                     double              asum[4],
                     double              ssum[4],
                     double              wssum[4])
{
  for (cs_lnum_t j = 0; j < 4; j++) {
    vmin[j] = HUGE_VAL;
    vmax[j] = -HUGE_VAL;
    vsum[j] = 0.;
    wsum[j] = 0.;
    asum[j] = 0.;
    ssum[j] = 0.;
    wssum[j] = 0.;
  }

# pragma omp parallel if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n, sizeof(cs_real_t), &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_lnum_t *_wl = wl + s_id;
    const cs_real_3_t *_v = v + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    cs_real_t lmin[4], lmax[4];
    for (int j = 0; j < 4; j++) {
      lmin[j] = HUGE_VAL;
      lmax[j] = -HUGE_VAL;
    }

    double lsum[20];
    for (int j = 0; j < 20; j++) lsum[j] = 0.;

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double s[20];
      for (cs_lnum_t j = 0; j < 20; j++) s[j] = 0.;

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {

        cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        if (end_id > _n)
          end_id = _n;

        double c[20];
        for (cs_lnum_t j = 0; j < 20; j++) c[j] = 0.;

        for (cs_lnum_t i = start_id; i < end_id; i++) {
          double v_norm2 = 0.;
          const cs_real_t  wi = w[_wl[i]];
          for (cs_lnum_t k = 0; k < 3; k++) {
            const cs_real_t  val = _v[i][k], val2 = val*val;
            c[k]   += val;
            c[k+4] += val*wi;
            c[k+8] += fabs(val);
            c[k+12] += val2;
            c[k+16] += val2*wi;
            v_norm2 += val2;
            if (val < lmin[k])
              lmin[k] = val;
            if (val > lmax[k])
              lmax[k] = val;
          }

          const double v_norm = sqrt(v_norm2);
          c[3] += v_norm;
          c[7] += v_norm*wi;
          c[11] += v_norm;
          c[15] += v_norm2;
          c[19] += v_norm2*wi;
          if (v_norm < lmin[3])
            lmin[3] = v_norm;
          if (v_norm > lmax[3])
            lmax[3] = v_norm;
        }

        for (cs_lnum_t j = 0; j < 20; j++)
          s[j] += c[j];

      } /* Loop on blocks */

      for (cs_lnum_t j = 0; j < 20; j++)
        lsum[j] += s[j];

    } /* Loop on super-block */

#   pragma omp critical
    {
      for (cs_lnum_t j = 0; j < 4; j++) {
        if (lmin[j] < vmin[j])
          vmin[j] = lmin[j];
        if (lmax[j] > vmax[j])
          vmax[j] = lmax[j];
        vsum[j] += lsum[j];
        wsum[j] += lsum[4+j];
        asum[j] += lsum[8+j];
        ssum[j] += lsum[12+j];
        wssum[j] += lsum[16+j];
      }
    }

  }
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, several sums and weighted sums)
 * of a 1-dimensional array.
 *
 * The array is interleaved.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n        <-- local number of elements
 *   vl       <-- pointer to weights list
 *   v        <-- pointer to field values (size: n*3)
 *   w        <-- pointer to weight values (size: n)
 *   vmin     --> resulting min array (size: 4)
 *   vmax     --> resulting max array (size: 4)
 *   vsum     --> resulting sum array (size: 4)
 *   wsum     --> resulting weighted sum array (size: 4)
 *   asum     --> resulting sum of absolute values (size: 4)
 *   ssum     --> resulting weighted sum array (size: 4)
 *   wssum    --> resulting weighted sum of squared values (size: 4)
 *----------------------------------------------------------------------------*/

static void
_cs_real_norms_3d_iv(cs_lnum_t           n,
                     const cs_lnum_t     vl[],
                     const cs_real_3_t   v[],
                     const cs_real_t     w[],
                     double              vmin[4],
                     double              vmax[4],
                     double              vsum[4],
                     double              wsum[4],
                     double              asum[4],
                     double              ssum[4],
                     double              wssum[4])
{
  for (cs_lnum_t j = 0; j < 4; j++) {
    vmin[j] = HUGE_VAL;
    vmax[j] = -HUGE_VAL;
    vsum[j] = 0.;
    wsum[j] = 0.;
    asum[j] = 0.;
    ssum[j] = 0.;
    wssum[j] = 0.;
  }

# pragma omp parallel if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n, sizeof(cs_real_t), &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_lnum_t *_vl = vl + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    cs_real_t lmin[4], lmax[4];
    for (int j = 0; j < 4; j++) {
      lmin[j] = HUGE_VAL;
      lmax[j] = -HUGE_VAL;
    }

    double lsum[20];
    for (int j = 0; j < 20; j++) lsum[j] = 0.;

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double s[20];
      for (cs_lnum_t j = 0; j < 20; j++) s[j] = 0.;

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {

        cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        if (end_id > _n)
          end_id = _n;

        double c[20];
        for (cs_lnum_t j = 0; j < 20; j++) c[j] = 0.;

        for (cs_lnum_t li = start_id; li < end_id; li++) {
          double v_norm2 = 0.;
          const cs_lnum_t  i = _vl[li];
          const cs_real_t  wi = w[i];
          for (cs_lnum_t k = 0; k < 3; k++) {
            const cs_real_t  val = v[i][k], val2 = val*val;
            c[k]   += val;
            c[k+4] += val*wi;
            c[k+8] += fabs(val);
            c[k+12] += val2;
            c[k+16] += val2*wi;
            v_norm2 += val2;
            if (val < lmin[k])
              lmin[k] = val;
            if (val > lmax[k])
              lmax[k] = val;
          }

          const double v_norm = sqrt(v_norm2);
          c[3] += v_norm;
          c[7] += v_norm*wi;
          c[11] += v_norm;
          c[15] += v_norm2;
          c[19] += v_norm2*wi;
          if (v_norm < lmin[3])
            lmin[3] = v_norm;
          if (v_norm > lmax[3])
            lmax[3] = v_norm;
        }

        for (cs_lnum_t j = 0; j < 20; j++)
          s[j] += c[j];

      } /* Loop on blocks */

      for (cs_lnum_t j = 0; j < 20; j++)
        lsum[j] += s[j];

    } /* Loop on super-block */

#   pragma omp critical
    {
      for (cs_lnum_t j = 0; j < 4; j++) {
        if (lmin[j] < vmin[j])
          vmin[j] = lmin[j];
        if (lmax[j] > vmax[j])
          vmax[j] = lmax[j];
        vsum[j] += lsum[j];
        wsum[j] += lsum[4+j];
        asum[j] += lsum[8+j];
        ssum[j] += lsum[12+j];
        wssum[j] += lsum[16+j];
      }
    }

  }
}

/*----------------------------------------------------------------------------
 * Compute simple local weighted norms (l1, l2) of a 3-dimensional array.
 * The array is interleaved.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n_src_elts  <-- number of source elements
 *   src2v_idx   <-- index array between source elements and local elements
 *   src2v_ids   <-- list of ids of local elements
 *   v           <-- pointer to array values
 *   w           <-- pointer to weights (size = src2v_idx[n_src_elts])
 *   vsum        --> (weighted) sum array (size: dim, or 4 if dim = 3)
 *   asum        --> (weighted) sum of absolute values (same size as vsum)
 *   ssum        --> (weighted) sum of squared values (same size as vsum)
 *----------------------------------------------------------------------------*/

static void
_cs_real_scatter_norms_3d(cs_lnum_t           n_src_elts,
                          const cs_lnum_t    *src2v_idx,
                          const cs_lnum_t    *src2v_ids,
                          const cs_real_3_t   v[],
                          const cs_real_t     w[],
                          double              vsum[],
                          double              asum[],
                          double              ssum[])
{
  /* Initialize quantities to return */

  for (int i = 0; i < 4; i++) vsum[i] = asum[i] = ssum[i] = 0.;

# pragma omp parallel if (n_src_elts > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n_src_elts, sizeof(cs_real_t), &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_lnum_t *_idx = src2v_idx + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    double lsum[12];
    for (int i = 0; i < 12; i++) lsum[i] = 0;

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double s[12];
      for (int i = 0; i < 12; i++) s[i] = 0;

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {

        const cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        if (end_id > _n)
          end_id = _n;

        double c[12];
        for (int i = 0; i < 12; i++) c[i] = 0;

        /* Loop on source elements */

        for (cs_lnum_t id = start_id; id < end_id; id++) {
          for (cs_lnum_t j = _idx[id]; j < _idx[id+1]; j++) {

            const cs_real_t  weight = w[j];
            const cs_lnum_t  elt_id = src2v_ids[j];

            double  v_norm2 = 0.;
            for (int k = 0; k < 3; k++) {
              const double  val = v[elt_id][k];
              c[k]   += val*weight;
              c[k+4] += fabs(val)*weight;
              c[k+8] += val*val*weight;
              v_norm2 += val*val;
            }

            const double  v_norm = sqrt(v_norm2);
            c[3] += weight*v_norm;
            c[7] += weight*v_norm;
            c[11] += weight*v_norm2;

          }
        }

        for (int i = 0; i < 12; i++) s[i] += c[i];

      } /* Loop on blocks */

      for (int i = 0; i < 12; i++) lsum[i] += s[i];

    } /* Loop on super-block */

#   pragma omp critical
    {
      for (int i = 0; i < 4; i++) {
        vsum[i] += lsum[i];
        asum[i] += lsum[i+4];
        ssum[i] += lsum[i+8];
      }
    }

  } /* Block openMP */
}

/*----------------------------------------------------------------------------
 * Compute simple local weighted norms (l1, l2) of a 3-dimensional array.
 * The array is interleaved.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n_src_elts  <-- number of source elements
 *   src2v_idx   <-- index array between source elements and local elements
 *   src2v_ids   <-- list of ids of local elements
 *   filter_list <-- list of source elements on which values are defined
 *   v           <-- pointer to array values
 *   w           <-- pointer to weights (size = src2v_idx[n_src_elts])
 *   vsum        --> (weighted) sum array (size: dim, or 4 if dim = 3)
 *   asum        --> (weighted) sum of absolute values (same size as vsum)
 *   ssum        --> (weighted) sum of squared values (same size as vsum)
 *----------------------------------------------------------------------------*/

static void
_cs_real_scatter_norms_3d_filtered(cs_lnum_t           n_src_elts,
                                   const cs_lnum_t    *src2v_idx,
                                   const cs_lnum_t    *src2v_ids,
                                   const cs_lnum_t    *filter_list,
                                   const cs_real_3_t   v[],
                                   const cs_real_t     w[],
                                   double              vsum[],
                                   double              asum[],
                                   double              ssum[])
{
  /* Initialize quantities to return */

  for (int i = 0; i < 4; i++) vsum[i] = asum[i] = ssum[i] = 0.;

# pragma omp parallel if (n_src_elts > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n_src_elts, sizeof(cs_real_t), &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_lnum_t *_fl = filter_list + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    double lsum[12];
    for (int i = 0; i < 12; i++) lsum[i] = 0;

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double s[12];
      for (int i = 0; i < 12; i++) s[i] = 0;

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {

        const cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        if (end_id > _n)
          end_id = _n;

        double c[12];
        for (int i = 0; i < 12; i++) c[i] = 0;

        /* Loop on source elements */

        for (cs_lnum_t lid = start_id; lid < end_id; lid++) {
          const cs_lnum_t  id = _fl[lid];
          for (cs_lnum_t j = src2v_idx[id]; j < src2v_idx[id+1]; j++) {

            const cs_real_t  weight = w[j];
            const cs_lnum_t  elt_id = src2v_ids[j];

            double  v_norm2 = 0.;
            for (int k = 0; k < 3; k++) {
              const double  val = v[elt_id][k];
              c[k]    += val*weight;
              c[k+4]  += fabs(val)*weight;
              c[k+8]  += val*val*weight;
              v_norm2 += val*val;
            }

            const double v_norm = sqrt(v_norm2);
            c[3]  += weight*v_norm;
            c[7]  += weight*v_norm;
            c[11] += weight*v_norm2;

          }
        }

        for (int i = 0; i < 12; i++) s[i] += c[i];

      } /* Loop on blocks */

      for (int i = 0; i < 12; i++) lsum[i] += s[i];

    } /* Loop on super-block */

#   pragma omp critical
    {
      for (int i = 0; i < 4; i++) {
        vsum[i] += lsum[i];
        asum[i] += lsum[i+4];
        ssum[i] += lsum[i+8];
      }
    }

  } /* Block openMP */
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum) of an
 * n-dimensional field array's components.
 *
 * The maximum allowed dimension is 10.
 * The array is interleaved.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n        <-- local number of elements
 *   dim      <-- local array dimension (max: 10)
 *   vl       <-- pointer to optional element list, or nullptr
 *   v        <-- pointer to field values (size: n*dim)
 *   vmin     --> resulting min array (size: dim)
 *   vmax     --> resulting max array (size: dim)
 *   vsum     --> resulting sum array (size: dim)
 *----------------------------------------------------------------------------*/

static void
_cs_real_sstats_nd(cs_lnum_t         n,
                   int               dim,
                   const cs_lnum_t  *vl,
                   const cs_real_t   v[],
                   double            vmin[],
                   double            vmax[],
                   double            vsum[])
{
  assert(dim <= 10);

  for (cs_lnum_t j = 0; j < dim; j++) {
    vmin[j] = HUGE_VAL;
    vmax[j] = -HUGE_VAL;
    vsum[j] = 0.;
  }

# pragma omp parallel if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n, sizeof(cs_real_t), &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    cs_real_t lmin[10], lmax[10];
    double lsum[10];

    for (cs_lnum_t j = 0; j < dim; j++) {
      lmin[j] = HUGE_VAL;
      lmax[j] = -HUGE_VAL;
      lsum[j] = 0;
    }

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double s[10];
      for (cs_lnum_t j = 0; j < dim; j++)
        s[j] = 0.;

      if (vl == nullptr) {

        const cs_real_t *_v = v + s_id*dim;

        for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
          cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
          cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
          if (end_id > _n)
            end_id = _n;
          double c[10];
          for (cs_lnum_t j = 0; j < dim; j++)
            c[j] = 0.;
          for (cs_lnum_t i = start_id; i < end_id; i++) {
            for (cs_lnum_t j = 0; j < dim; j++) {
              c[j] += _v[i*dim + j];
              if (_v[i*dim + j] < lmin[j])
                lmin[j] = _v[i*dim + j];
              if (_v[i*dim+j] > lmax[j])
                lmax[j] = _v[i*dim+j];
            }
          }
          for (cs_lnum_t j = 0; j < dim; j++)
            s[j] += c[j];
        }

      }
      else {

        const cs_lnum_t *_vl = vl + s_id;

        for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
          cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
          cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
          if (end_id > _n)
            end_id = _n;
          double c[10];
          for (cs_lnum_t j = 0; j < dim; j++)
            c[j] = 0.;
          for (cs_lnum_t li = start_id; li < end_id; li++) {
            cs_lnum_t i = _vl[li];
            for (cs_lnum_t j = 0; j < dim; j++) {
              c[j] += v[i*dim + j];
              if (v[i*dim + j] < lmin[j])
                lmin[j] = v[i*dim + j];
              if (v[i*dim+j] > lmax[j])
                lmax[j] = v[i*dim+j];
            }
          }
          for (cs_lnum_t j = 0; j < dim; j++)
            s[j] += c[j];
        }

      }

      for (cs_lnum_t j = 0; j < dim; j++)
        lsum[j] += s[j];

    }

#   pragma omp critical
    {
      for (cs_lnum_t j = 0; j < dim; j++) {
        if (lmin[j] < vmin[j])
          vmin[j] = lmin[j];
        if (lmax[j] > vmax[j])
          vmax[j] = lmax[j];
        vsum[j] += lsum[j];
      }
    }

  }
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum, weighted sum) of an
 * n-dimensional field array's components.
 *
 * The maximum allowed dimension is 10.
 * The array is interleaved.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n        <-- local number of elements
 *   dim      <-- local array dimension (max: 10)
 *   vl       <-- pointer to optional element list, or nullptr
 *   wl       <-- pointer to optional element list, or nullptr
 *                (ignored if vl != nullptr)
 *   v        <-- pointer to field values (size: n*dim)
 *   w        <-- pointer to weight values (size: n)
 *   vmin     --> resulting min array (size: dim)
 *   vmax     --> resulting max array (size: dim)
 *   vsum     --> resulting sum array (size: dim)
 *   wsum     --> resulting weighted sum array (size: dim)
 *----------------------------------------------------------------------------*/

static void
_cs_real_sstats_nd_w(cs_lnum_t         n,
                     int               dim,
                     const cs_lnum_t  *vl,
                     const cs_lnum_t  *wl,
                     const cs_real_t   v[],
                     const cs_real_t   w[],
                     double            vmin[],
                     double            vmax[],
                     double            vsum[],
                     double            wsum[])
{
  assert(dim <= 10);

  for (cs_lnum_t j = 0; j < dim; j++) {
    vmin[j] = HUGE_VAL;
    vmax[j] = -HUGE_VAL;
    vsum[j] = 0.;
    wsum[j] = 0.;
  }

# pragma omp parallel if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n, sizeof(cs_real_t), &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    const int dim2 = dim*2;

    cs_real_t lmin[10], lmax[10];
    double lsum[20];

    for (cs_lnum_t j = 0; j < dim; j++) {
      lmin[j] = HUGE_VAL;
      lmax[j] = -HUGE_VAL;
      lsum[j] = 0;
      lsum[j+dim] = 0;
    }

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double s[20];
      for (cs_lnum_t j = 0; j < dim2; j++)
        s[j] = 0.;

      if (vl == nullptr && wl == nullptr) {

        const cs_real_t *_v = v + s_id*dim;
        const cs_real_t *_w = w + s_id;

        for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
          cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
          cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
          if (end_id > _n)
            end_id = _n;
          double c[20];
          for (cs_lnum_t j = 0; j < dim2; j++)
            c[j] = 0.;
          for (cs_lnum_t i = start_id; i < end_id; i++) {
            for (cs_lnum_t j = 0; j < dim; j++) {
              c[j]     += _v[i*dim + j];
              c[j+dim] += _v[i*dim + j]*_w[i];
              if (_v[i*dim + j] < lmin[j])
                lmin[j] = _v[i*dim + j];
              if (_v[i*dim+j] > lmax[j])
                lmax[j] = _v[i*dim+j];
            }
          }
          for (cs_lnum_t j = 0; j < dim2; j++)
            s[j] += c[j];
        }
      }

      else if (vl == nullptr) {

        const cs_lnum_t *_wl = wl + s_id;
        const cs_real_t *_v = v + s_id*dim;

        for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
          cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
          cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
          if (end_id > _n)
            end_id = _n;
          double c[20];
          for (cs_lnum_t j = 0; j < dim2; j++)
            c[j] = 0.;
          for (cs_lnum_t i = start_id; i < end_id; i++) {
            cs_real_t wi = w[_wl[i]];
            for (cs_lnum_t j = 0; j < dim; j++) {
              c[j]     += _v[i*dim + j];
              c[j+dim] += _v[i*dim + j]*wi;
              if (_v[i*dim + j] < lmin[j])
                lmin[j] = _v[i*dim + j];
              if (_v[i*dim+j] > lmax[j])
                lmax[j] = _v[i*dim+j];
            }
          }
          for (cs_lnum_t j = 0; j < dim2; j++)
            s[j] += c[j];
        }

      }
      else { /* vl != nullptr */

        const cs_lnum_t *_vl = vl + s_id;

        for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
          cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
          cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
          if (end_id > _n)
            end_id = _n;
          double c[20];
          for (cs_lnum_t j = 0; j < dim2; j++)
            c[j] = 0.;
          for (cs_lnum_t li = start_id; li < end_id; li++) {
            cs_lnum_t i = _vl[li];
            for (cs_lnum_t j = 0; j < dim; j++) {
              c[j]     += v[i*dim + j];
              c[j+dim] += v[i*dim + j]*w[i];
              if (v[i*dim + j] < lmin[j])
                lmin[j] = v[i*dim + j];
              if (v[i*dim+j] > lmax[j])
                lmax[j] = v[i*dim+j];
            }
          }
          for (cs_lnum_t j = 0; j < dim2; j++)
            s[j] += c[j];
        }

      }

      for (cs_lnum_t j = 0; j < dim2; j++)
        lsum[j] += s[j];

    }

#   pragma omp critical
    {
      for (cs_lnum_t j = 0; j < dim; j++) {
        if (lmin[j] < vmin[j])
          vmin[j] = lmin[j];
        if (lmax[j] > vmax[j])
          vmax[j] = lmax[j];
        vsum[j] += lsum[j];
        wsum[j] += lsum[dim+j];
      }
    }

  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the min./max. of a 1-dimensional array.
 *
 * \param[in]   n        local number of elements
 * \param[in]   v        pointer to values (size: n)
 * \param[out]  vmin     minimum value
 * \param[out]  vmax     maximum value
 */
/*----------------------------------------------------------------------------*/

void
cs_array_reduce_minmax(cs_lnum_t          n,
                       const cs_real_t    v[],
                       cs_real_t         &vmin,
                       cs_real_t         &vmax)
{
  vmin = HUGE_VAL;
  vmax = -HUGE_VAL;

# pragma omp parallel if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n, sizeof(cs_real_t), &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_real_t *_v = v + s_id;

    cs_real_t lmin = HUGE_VAL;
    cs_real_t lmax = -HUGE_VAL;

    for (cs_lnum_t i = 0; i < _n; i++) {
      if (_v[i] < lmin)
        lmin = _v[i];
      if (_v[i] > lmax)
        lmax = _v[i];
    }

#   pragma omp critical
    {
      if (lmin < vmin)
        vmin = lmin;
      if (lmax > vmax)
        vmax = lmax;
    }

  } /* openMP block */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute simple local stats (minima, maxima, sum) of an
 * n-dimensional cs_real_t array's components.
 *
 * The array is interleaved.
 *
 * For arrays of dimension 3, the statistics relative to the norm
 * are also computed, and added at the end of the statistics arrays
 * (which must be size dim+1).
 *
 * The algorithm here is similar to that used for BLAS, but computes several
 * quantities simultaneously for better cache behavior
 *
 * \param[in]   ctx         cs_dispatch_context struct from higher level
 * \param[in]   n_elts      number of local elements
 * \param[in]   v_elt_list  optional list of parent elements on which values
 *                          are defined, or nullptr
 * \param[in]   v           pointer to array values
 * \param[out]  vmin        resulting min array (size: dim, or 4 if dim = 3)
 * \param[out]  vmax        resulting max array (size: dim, or 4 if dim = 3)
 * \param[out]  vsum        resulting sum array (size: dim, or 4 if dim = 3)
 */
/*----------------------------------------------------------------------------*/

void
cs_array_reduce_simple_stats_l(cs_dispatch_context  ctx,
                               const cs_lnum_t      n_elts,
                               const int            dim,
                               const cs_lnum_t     *v_elt_list,
                               const cs_real_t      v[],
                               double               vmin[],
                               double               vmax[],
                               double               vsum[])
{
  /* If all values are defined on same list */

  if (v_elt_list == nullptr) {
    if (dim == 1)
      _cs_real_sstats<1>(ctx, n_elts, v, vmin, vmax, vsum);
    else if (dim == 3)
      _cs_real_sstats_with_norm<3>(ctx, n_elts, v, vmin, vmax, vsum);
    else if (dim == 6)
      _cs_real_sstats_with_norm<6>(ctx, n_elts, v, vmin, vmax, vsum);
    else /* dim is only known at runtime, so can not be template parameter */
      _cs_real_sstats_nd(n_elts, dim, nullptr, v, vmin, vmax, vsum);
  }

  /* If values are defined on parent list */

  else {
    if (dim == 1)
      _cs_real_sstats_iv<1>(ctx, n_elts, v_elt_list, v, vmin, vmax, vsum);
    else if (dim == 3)
      _cs_real_sstats_with_norm_iv<3>(ctx, n_elts, v_elt_list,
          v, vmin, vmax, vsum);
    else if (dim == 6)
      _cs_real_sstats_with_norm_iv<6>(ctx, n_elts, v_elt_list,
          v, vmin, vmax, vsum);
    else /* dim is only known at runtime, so can not be template parameter */
      _cs_real_sstats_nd(n_elts, dim, v_elt_list, v, vmin, vmax, vsum);
  }
}

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute weighted sums of an n-dimensional cs_real_t array's
 * components.
 *
 * The maximum allowed dimension is 9 (allowing for a rank-2 tensor).
 * The array is interleaved.
 *
 * For arrays of dimension 3, the statistics relative to the norm
 * are also computed, and added at the end of the statistics arrays
 * (which must be size dim+1).
 *
 * The algorithm here is similar to that used for BLAS.
 *
 * \param[in]   n_elts      number of local elements
 * \param[in]   dim         local array dimension (max: 9)
 * \param[in]   v_elt_list  optional list of parent elements on which values
 *                          are defined, or nullptr
 * \param[in]   w_elt_list  optional list of parent elements on which weights
 *                          are defined, or nullptr; if v_elt_list is defined
 *                          (ie. non-null),then w_elt_list = v_elt_list is
 *                          assumed, so this parameter is ignored
 * \param[in]   v           pointer to array values
 * \param[in]   w           pointer to weights
 * \param[out]  wsum        resulting weighted sum array
 *                          (size: dim, or 4 if dim = 3)
 */
/*----------------------------------------------------------------------------*/

void
cs_array_reduce_wsum_l(cs_lnum_t         n_elts,
                       int               dim,
                       const cs_lnum_t  *v_elt_list,
                       const cs_lnum_t  *w_elt_list,
                       const cs_real_t   v[],
                       const cs_real_t   w[],
                       double            wsum[])
{
  /* If all values are defined on same list */

  if (v_elt_list == nullptr && w_elt_list == nullptr) {
    if (dim == 1)
      wsum[0] = _cs_real_wsum_1d(n_elts, v, w);
    else if (dim == 3)
      bft_error(__FILE__, __LINE__, 0,
                _("_cs_real_wsum_3d not implemented yet\n"));
    else
      bft_error(__FILE__, __LINE__, 0,
                _("_cs_real_wsum_nd not implemented yet\n"));
  }

  /* If weights are defined on weights list only */

  else if (v_elt_list == nullptr) { /* w_elt_list != nullptr */
    if (dim == 1)
      wsum[0] = _cs_real_wsum_1d_iw(n_elts, w_elt_list, v, w);
    else if (dim == 3)
      bft_error(__FILE__, __LINE__, 0,
                _("_cs_real_wsum_3d_iw not implemented yet\n"));
    else
      bft_error(__FILE__, __LINE__, 0,
                _("_cs_real_wsum_nd_iw not implemented yet\n"));
  }

  /* If weights are defined on parent list */

  else { /* v_elt_list != nullptr */
    if (dim == 1)
      wsum[0] = _cs_real_wsum_1d_iv(n_elts, v_elt_list, v, w);
    else if (dim == 3)
      bft_error(__FILE__, __LINE__, 0,
                _("_cs_real_wsum_3d_iv not implemented yet\n"));
    else
      bft_error(__FILE__, __LINE__, 0,
                _("_cs_real_wsum_nd_iv not implemented yet\n"));
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute weighted sums of an n-dimensional cs_real_t array's
 * components. Output is both weighted sum and sum of weights, hence allowing
 * for the computation of local, or global using parallel operations afterwards,
 * mean of the array.
 *
 * The maximum allowed dimension is 9 (allowing for a rank-2 tensor).
 * The array is interleaved.
 *
 * For arrays of dimension 3, the statistics relative to the norm
 * are also computed, and added at the end of the statistics arrays
 * (which must be size dim+1).
 *
 * The algorithm here is similar to that used for BLAS.
 *
 * \param[in]   n_elts      number of local elements
 * \param[in]   dim         local array dimension (max: 9)
 * \param[in]   v_elt_list  optional list of parent elements on which values
 *                          are defined, or nullptr
 * \param[in]   w_elt_list  optional list of parent elements on which weights
 *                          are defined, or nullptr; if v_elt_list is defined
 *                          (ie. non-null),then w_elt_list = v_elt_list is
 *                          assumed, so this parameter is ignored
 * \param[in]   v           pointer to array values
 * \param[in]   w           pointer to weights
 * \param[out]  wsum        resulting weighted sum array
 *                          (size: dim, or 4 if dim = 3)
 * \param[out]  wtot        resulting local sum of weights array
 */
/*----------------------------------------------------------------------------*/

void
cs_array_reduce_wsum_components_l(cs_lnum_t         n_elts,
                                  int               dim,
                                  const cs_lnum_t  *v_elt_list,
                                  const cs_lnum_t  *w_elt_list,
                                  const cs_real_t   v[],
                                  const cs_real_t   w[],
                                  double            wsum[],
                                  double            wtot[])
{
  /* If all values are defined on same list */

  if (v_elt_list == nullptr && w_elt_list == nullptr) {
    if (dim == 1)
      _cs_real_wsum_components_1d(n_elts, v, w, wsum, wtot);
    else if (dim == 3)
      bft_error(__FILE__, __LINE__, 0,
                _("_cs_real_wsum_3d not implemented yet\n"));
    else
      bft_error(__FILE__, __LINE__, 0,
                _("_cs_real_wsum_nd not implemented yet\n"));
  }

  /* If weights are defined on weights list only */

  else if (v_elt_list == nullptr) { /* w_elt_list != nullptr */
    if (dim == 1)
      _cs_real_wsum_components_1d_iw(n_elts, w_elt_list, v, w, wsum, wtot);
    else if (dim == 3)
      bft_error(__FILE__, __LINE__, 0,
                _("_cs_real_wsum_3d_iw not implemented yet\n"));
    else
      bft_error(__FILE__, __LINE__, 0,
                _("_cs_real_wsum_nd_iw not implemented yet\n"));
  }

  /* If weights are defined on parent list */

  else { /* v_elt_list != nullptr */
    if (dim == 1)
      _cs_real_wsum_components_1d_iv(n_elts, v_elt_list, v, w, wsum, wtot);
    else if (dim == 3)
      bft_error(__FILE__, __LINE__, 0,
                _("_cs_real_wsum_3d_iv not implemented yet\n"));
    else
      bft_error(__FILE__, __LINE__, 0,
                _("_cs_real_wsum_nd_iv not implemented yet\n"));
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute sums of an n-dimensional cs_real_t array's components.
 *
 * The maximum allowed dimension is 3. The array is interleaved.
 *
 * For arrays of dimension 3, the statistics relative to the norm
 * are also computed, and added at the end of the statistics arrays
 * (which must be size dim+1).
 *
 * The algorithm here is similar to that used for BLAS.
 *
 * \param[in]   n_elts      number of local elements
 * \param[in]   dim         local array dimension (max: 9)
 * \param[in]   v_elt_list  optional list of parent elements on which values
 *                          are defined, or nullptr
 * \param[in]   v           pointer to array values
 * \param[out]  vmin        resulting min array (size: dim, or 4 if dim = 3)
 * \param[out]  vmax        resulting max array (size: dim, or 4 if dim = 3)
 */
/*----------------------------------------------------------------------------*/

void
cs_array_reduce_minmax_l(cs_lnum_t         n_elts,
                         int               dim,
                         const cs_lnum_t  *v_elt_list,
                         const cs_real_t   v[],
                         cs_real_t         vmin[],
                         cs_real_t         vmax[])
{
  /* If all values are defined on same list */

  if (v_elt_list == nullptr) {
    if (dim == 1)
      cs_array_reduce_minmax(n_elts, v, vmin[0], vmax[0]);
    else if (dim == 3)
      _cs_real_minmax_3d(n_elts, (const cs_real_3_t *)v, vmin, vmax);
    else
      bft_error(__FILE__, __LINE__, 0,
                _("_cs_real_minmax_nd not implemented yet\n"));
  }

  /* If values are defined on parent list */

  else {
    if (dim == 1)
      bft_error(__FILE__, __LINE__, 0,
                _("cs_array_reduce_minmax_iv not implemented yet\n"));
    else if (dim == 3)
      bft_error(__FILE__, __LINE__, 0,
                _("_cs_real_minmax_3d_iv not implemented yet\n"));
    else
      bft_error(__FILE__, __LINE__, 0,
                _("_cs_real_minmax_nd_iv not implemented yet\n"));
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute simple local stats (minima, maxima, sum, weighted sum) of
 * an n-dimensional cs_real_t array's components for a given mesh location.
 *
 * The maximum allowed dimension is 10.
 * The array is interleaved.
 *
 * For arrays of dimension 3, the statistics relative to the norm
 * are also computed, and added at the end of the statistics arrays
 * (which must be size dim+1).
 *
 * The algorithm here is similar to that used for BLAS, but computes several
 * quantities simultaneously for better cache behavior
 *
 * \param[in]   ctx         cs_dispatch_context struct from higher level
 * \param[in]   n_elts      number of local elements
 * \param[in]   dim         local array dimension (max: 10)
 * \param[in]   v_elt_list  optional list of parent elements on which values
 *                          are defined, or nullptr
 * \param[in]   w_elt_list  optional list of parent elements on which weights
 *                          are defined, or nullptr; if v_elt_list is defined
 *                          (ie. non-null),then w_elt_list = v_elt_list is
 *                          assumed, so this parameter is ignored
 * \param[in]   v     pointer to array values
 * \param[in]   w     pointer to weights
 * \param[out]  vmin  resulting min array (size: dim, or 4 if dim = 3)
 * \param[out]  vmax  resulting max array (size: dim, or 4 if dim = 3)
 * \param[out]  vsum  resulting sum array (size: dim, or 4 if dim = 3)
 * \param[out]  wsum  resulting weighted sum array (size: dim, or 4 if dim = 3)
 */
/*----------------------------------------------------------------------------*/

void
cs_array_reduce_simple_stats_l_w(cs_dispatch_context  ctx,
                                 cs_lnum_t            n_elts,
                                 int                  dim,
                                 const cs_lnum_t     *v_elt_list,
                                 const cs_lnum_t     *w_elt_list,
                                 const cs_real_t      v[],
                                 const cs_real_t      w[],
                                 double               vmin[],
                                 double               vmax[],
                                 double               vsum[],
                                 double               wsum[])
{
  /* If all values are defined on same list */

  if (v_elt_list == nullptr && w_elt_list == nullptr) {
    if (dim == 1)
      _cs_real_sstats_weighted<1>(ctx, n_elts, v, w, vmin, vmax, vsum, wsum);
    else if (dim == 3)
      _cs_real_sstats_weighted_with_norm<3>(ctx, n_elts, v, w,
          vmin, vmax, vsum, wsum);
    else if (dim == 6)
      _cs_real_sstats_weighted_with_norm<6>(ctx, n_elts, v, w,
          vmin, vmax, vsum, wsum);
    else /* dim is only known at runtime, so can not be template parameter */
      _cs_real_sstats_nd_w(n_elts, dim, nullptr, nullptr, v, w,
          vmin, vmax, vsum, wsum);
  }

  /* If weights are defined on parent list */

  else if (v_elt_list == nullptr) { /* w_elt_list != nullptr */
    if (dim == 1)
      _cs_real_sstats_weighted_iw<1>(ctx, n_elts, w_elt_list,
          v, w, vmin, vmax, vsum, wsum);
    else if (dim == 3)
      _cs_real_sstats_weighted_with_norm_iw<3>(ctx, n_elts, w_elt_list, v, w,
          vmin, vmax, vsum, wsum);
    else if (dim == 6)
      _cs_real_sstats_weighted_with_norm_iw<6>(ctx, n_elts, w_elt_list, v, w,
          vmin, vmax, vsum, wsum);
    else /* dim is only known at runtime, so can not be template parameter */
      _cs_real_sstats_nd_w(n_elts, dim, nullptr, w_elt_list, v, w,
          vmin, vmax, vsum, wsum);
  }

  /* If weights are defined on parent list */

  else { /* v_elt_list != nullptr */
    if (dim == 1)
      _cs_real_sstats_weighted_iv<1>(ctx, n_elts, v_elt_list,
          v, w, vmin, vmax, vsum, wsum);
    else if (dim == 3)
      _cs_real_sstats_weighted_with_norm_iv<3>(ctx, n_elts, v_elt_list, v, w,
          vmin, vmax, vsum, wsum);
    else if (dim == 6)
      _cs_real_sstats_weighted_with_norm_iv<6>(ctx, n_elts, v_elt_list, v, w,
          vmin, vmax, vsum, wsum);
    else /* dim is only known at runtime, so can not be template parameter */
      _cs_real_sstats_nd_w(n_elts, dim, v_elt_list, nullptr, v, w,
          vmin, vmax, vsum, wsum);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute simple local stats and norms (minima, maxima, sum, weighted
 *         sum, sum of absolute values, sum of squared values and weighted sum
 *         of squared values) of an n-dimensional cs_real_t array's components
 *
 * The maximum allowed dimension is 3.
 * The array is interleaved.
 *
 * For arrays of dimension 3, the statistics relative to the norm
 * are also computed, and added at the end of the statistics arrays
 * (which must be size dim+1).
 *
 * The algorithm here is similar to that used for BLAS, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n_elts     <-- number of local elements
 *   dim        <-- local array dimension (max: 3)
 *   v_elt_list <-- optional list of parent elements on which values
 *                  are defined, or nullptr
 *   w_elt_list <-- optional list of parent elements on which weights
 *                  are defined, or nullptr; if v_elt_list is defined
 *                  (ie. non-null),then w_elt_list = v_elt_list is assumed,
 *                  so this parameter is ignored
 *   v          <-- pointer to array values
 *   w          <-- pointer to weights
 *   vmin       --> resulting min array (size: dim, or 4 if dim = 3)
 *   vmax       --> resulting max array (same size as vmin)
 *   vsum       --> resulting sum array (same size as vmin)
 *   wsum       --> resulting weighted sum array (same size as vmin)
 *   asum       --> resulting sum of absolute values (same size as vmin)
 *   ssum       --> resulting weighted sum array (same size as vmin)
 *   wssum      --> resulting weighted sum of squared values (same size as vmin)
 */
/*----------------------------------------------------------------------------*/

void
cs_array_reduce_simple_norms_l(cs_lnum_t         n_elts,
                               int               dim,
                               const cs_lnum_t  *v_elt_list,
                               const cs_lnum_t  *w_elt_list,
                               const cs_real_t   v[],
                               const cs_real_t   w[],
                               double            vmin[],
                               double            vmax[],
                               double            vsum[],
                               double            wsum[],
                               double            asum[],
                               double            ssum[],
                               double            wssum[])
{
  /* If all values are defined on same list */

  if (v_elt_list == nullptr && w_elt_list == nullptr) {
    if (dim == 1)
      _cs_real_norms_1d(n_elts, v, w,
                        vmin, vmax, vsum, wsum, asum, ssum, wssum);
    else if (dim == 3)
      _cs_real_norms_3d(n_elts, (const cs_real_3_t *)v, w,
                        vmin, vmax, vsum, wsum, asum, ssum, wssum);
    else
      bft_error(__FILE__, __LINE__, 0,
                _("_cs_real_norms_nd not implemented yet\n"));
  }

  /* If weights are defined on parent list */

  else if (v_elt_list == nullptr) { /* w_elt_list != nullptr */
    if (dim == 1)
      _cs_real_norms_1d_iw(n_elts, w_elt_list, v, w,
                           vmin, vmax, vsum, wsum, asum, ssum, wssum);
    else if (dim == 3)
      _cs_real_norms_3d_iw(n_elts, w_elt_list, (const cs_real_3_t *)v, w,
                           vmin, vmax, vsum, wsum, asum, ssum, wssum);
    else
      bft_error(__FILE__, __LINE__, 0,
                _("_cs_real_norms_nd_iw not implemented yet\n"));
  }

  /* If weights are defined on parent list */

  else { /* v_elt_list != nullptr */
    if (dim == 1)
      _cs_real_norms_1d_iv(n_elts, v_elt_list, v, w,
                           vmin, vmax, vsum, wsum, asum, ssum, wssum);
    else if (dim == 3)
      _cs_real_norms_3d_iv(n_elts, v_elt_list, (const cs_real_3_t *)v, w,
                           vmin, vmax, vsum, wsum, asum, ssum, wssum);
    else
      bft_error(__FILE__, __LINE__, 0,
                _("_cs_real_norms_nd_iv not implemented yet\n"));
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute simple local weighted norms (l1, l2) of an n-dimensional
 *         cs_real_t array's components
 *         The weight array is mandatory
 *
 * The maximum allowed dimension is 3.
 * The array is interleaved.
 *
 * For arrays of dimension 3, the statistics relative to the norm
 * are also computed, and added at the end of the statistics arrays
 * (which must be size dim+1).
 *
 * The algorithm here is similar to that used for BLAS, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n_src_elts  <-- number of source elements
 *   src2v_idx   <-- index array between source elements and local elements
 *   src2v_ids   <-- list of ids of local elements
 *   filter_list <-- optional list of source elements on which values
 *                   are defined, or nullptr
 *   dim         <-- local array dimension (max: 3)
 *   n_v_elts    <-- number of local elements in the array of values
 *   v           <-- pointer to array values
 *   w           <-- pointer to weights (size = src2v_idx[n_src_elts])
 *   vsum        --> (weighted) sum array (size: dim, or 4 if dim = 3)
 *   asum        --> (weighted) sum of absolute values (same size as vsum)
 *   ssum        --> (weighted) sum of squared values (same size as vsum)
 */
/*----------------------------------------------------------------------------*/

void
cs_array_scatter_reduce_norms_l(cs_lnum_t          n_src_elts,
                                const cs_lnum_t   *src2v_idx,
                                const cs_lnum_t   *src2v_ids,
                                const cs_lnum_t   *filter_list,
                                int                dim,
                                cs_lnum_t          n_v_elts,
                                const cs_real_t    v[],
                                const cs_real_t    w[],
                                double             vsum[],
                                double             asum[],
                                double             ssum[])
{
  CS_UNUSED(n_v_elts); /* Useful to check coherency */

  /* If all values are defined on same list */

  if (filter_list == nullptr) {
    if (dim == 1)
      _cs_real_scatter_norms_1d(n_src_elts, src2v_idx, src2v_ids,
                                v, w,
                                vsum, asum, ssum);
    else if (dim == 3)
      _cs_real_scatter_norms_3d(n_src_elts, src2v_idx, src2v_ids,
                                (const cs_real_3_t *)v, w,
                                vsum, asum, ssum);
    else
      bft_error(__FILE__, __LINE__, 0,
                _(" _cs_real_scatter_norms_nd not implemented yet\n"));
  }
  else { /* filter_list != nullptr */

    if (dim == 1)
      _cs_real_scatter_norms_1d_filtered(n_src_elts, src2v_idx, src2v_ids,
                                         filter_list, v, w,
                                         vsum, asum, ssum);
    else if (dim == 3)
      _cs_real_scatter_norms_3d_filtered(n_src_elts, src2v_idx, src2v_ids,
                                         filter_list, (const cs_real_3_t *)v, w,
                                         vsum, asum, ssum);
    else
      bft_error(__FILE__, __LINE__, 0,
                _(" _cs_real_scatter_norms_nd_filtered not implemented yet\n"));

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
