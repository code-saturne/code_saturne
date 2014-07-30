/*============================================================================
 * Common array reduction operations.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_array_reduce.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/* Minimum size for OpenMP loops (needs benchmarking to adjust) */
#define THR_MIN 128

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
 * Compute sum of a 1-dimensional array.
 *
 * The algorithm here is similar to that used for blas.
 *
 * parameters:
 *   n        <-- local number of elements
 *   v        <-- pointer to values (size: n)
 *
 * returns:
 *   resulting sum
 *----------------------------------------------------------------------------*/

static double
_cs_real_sum_1d(cs_lnum_t        n,
                const cs_real_t  v[])
{
  const cs_lnum_t block_size = 60;

  cs_lnum_t i;
  cs_lnum_t sid, bid;
  cs_lnum_t start_id, end_id;
  double c, s;

  cs_lnum_t n_blocks = n / block_size;
  cs_lnum_t n_sblocks = sqrt(n_blocks);
  cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

  double v_sum = 0.;

# pragma omp parallel private(bid, start_id, end_id, i, c, s)    \
                              reduction(+:v_sum) if (n > THR_MIN)
  {
    # pragma omp for
    for (sid = 0; sid < n_sblocks; sid++) {

      s = 0.0;

      for (bid = 0; bid < blocks_in_sblocks; bid++) {
        start_id = block_size * (blocks_in_sblocks*sid + bid);
        end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        c = 0.0;
        for (i = start_id; i < end_id; i++)
          c += v[i];
        s += c;
      }

      v_sum += s;

    }

  } /* End of OpenMP-threaded section */

  start_id = block_size * n_sblocks*blocks_in_sblocks;
  end_id = n;
  c = 0.0;
  for (i = start_id; i < end_id; i++)
    c += v[i];

  v_sum += c;

  return v_sum;
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum) of a 1-dimensional array.
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

static void
_cs_real_sstats_1d(cs_lnum_t         n,
                   const cs_real_t   v[],
                   double           *vmin,
                   double           *vmax,
                   double           *vsum)
{
  const cs_lnum_t block_size = 60;

  cs_lnum_t i;
  cs_lnum_t sid, bid;
  cs_lnum_t start_id, end_id;
  double c, s, lmin, lmax;

  cs_lnum_t n_blocks = n / block_size;
  cs_lnum_t n_sblocks = sqrt(n_blocks);
  cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

#if defined(__xlc__)
#pragma disjoint(*v, *vmin, *vmax, *vsum)
#endif

  *vmin = HUGE_VAL;
  *vmax = -HUGE_VAL;
  *vsum = 0.;

# pragma omp parallel private(bid, start_id, end_id, i, \
                              c, s, lmin, lmax) if (n > THR_MIN)
  {
    # pragma omp for
    for (sid = 0; sid < n_sblocks; sid++) {

      lmin = HUGE_VAL;
      lmax = -HUGE_VAL;
      s = 0.0;

      for (bid = 0; bid < blocks_in_sblocks; bid++) {
        start_id = block_size * (blocks_in_sblocks*sid + bid);
        end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        c = 0.0;
        for (i = start_id; i < end_id; i++) {
          c += v[i];
          if (v[i] < lmin)
            lmin = v[i];
          if (v[i] > lmax)
            lmax = v[i];
        }
        s += c;
      }

#     pragma omp critical
      {
        if (lmin < *vmin)
          *vmin = lmin;
        if (lmax > *vmax)
          *vmax = lmax;
        *vsum += s;
      }

    }
  } /* End of OpenMP-threaded section */

  start_id = block_size * n_sblocks*blocks_in_sblocks;
  end_id = n;
  c = 0.0;
  for (i = start_id; i < end_id; i++) {
    c += v[i];
    if (v[i] < *vmin)
      *vmin = v[i];
    if (v[i] > *vmax)
      *vmax = v[i];
  }

  *vsum += c;
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum) of a
 * subset of a 1-dimensional array.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n        <-- local number of elements
 *   vl       <-- pointer to elements list
 *   v        <-- pointer to element values (size: n)
 *   vmin     --> resulting min
 *   vmax     --> resulting max
 *   vsum     --> resulting sum
 *----------------------------------------------------------------------------*/

static void
_cs_real_sstats_1d_l(cs_lnum_t         n,
                     const cs_lnum_t   vl[],
                     const cs_real_t   v[],
                     double           *vmin,
                     double           *vmax,
                     double           *vsum)
{
  const cs_lnum_t block_size = 60;

  cs_lnum_t i, li;
  cs_lnum_t sid, bid;
  cs_lnum_t start_id, end_id;
  double c, s, lmin, lmax;

  cs_lnum_t n_blocks = n / block_size;
  cs_lnum_t n_sblocks = sqrt(n_blocks);
  cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

#if defined(__xlc__)
#pragma disjoint(*vl, *v, *vmin, *vmax, *vsum)
#endif

  *vmin = HUGE_VAL;
  *vmax = -HUGE_VAL;
  *vsum = 0.;

# pragma omp parallel private(bid, start_id, end_id, li, i, \
                              c, s, lmin, lmax) if (n > THR_MIN)
  {
    # pragma omp for
    for (sid = 0; sid < n_sblocks; sid++) {

      lmin = HUGE_VAL;
      lmax = -HUGE_VAL;
      s = 0.0;

      for (bid = 0; bid < blocks_in_sblocks; bid++) {
        start_id = block_size * (blocks_in_sblocks*sid + bid);
        end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        c = 0.0;
        for (li = start_id; li < end_id; li++) {
          i = vl[li];
          c += v[i];
          if (v[i] < lmin)
            lmin = v[i];
          if (v[i] > lmax)
            lmax = v[i];
        }
        s += c;
      }

#     pragma omp critical
      {
        if (lmin < *vmin)
          *vmin = lmin;
        if (lmax > *vmax)
          *vmax = lmax;
        *vsum += s;
      }

    }
  } /* End of OpenMP-threaded section */

  start_id = block_size * n_sblocks*blocks_in_sblocks;
  end_id = n;
  c = 0.0;
  for (li = start_id; li < end_id; li++) {
    i = vl[li];
    c += v[i];
    if (v[i] < *vmin)
      *vmin = v[i];
    if (v[i] > *vmax)
      *vmax = v[i];
  }

  *vsum += c;
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum, weighted sum) of a
 * 1-dimensional array.
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
 *----------------------------------------------------------------------------*/

static void
_cs_real_sstats_1d_w(cs_lnum_t         n,
                     const cs_real_t   v[],
                     const cs_real_t   w[],
                     double           *vmin,
                     double           *vmax,
                     double           *vsum,
                     double           *wsum)
{
  const cs_lnum_t block_size = 60;

  cs_lnum_t i;
  cs_lnum_t sid, bid;
  cs_lnum_t start_id, end_id;
  double c[2], s[2], lmin, lmax;

  cs_lnum_t n_blocks = n / block_size;
  cs_lnum_t n_sblocks = sqrt(n_blocks);
  cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

#if defined(__xlc__)
#pragma disjoint(*v, *w, *vmin, *vmax, *vsum, *wsum)
#endif

  *vmin = HUGE_VAL;
  *vmax = -HUGE_VAL;
  *vsum = 0.;
  *wsum = 0.;

# pragma omp parallel private(bid, start_id, end_id, i, \
                              c, s, lmin, lmax) if (n > THR_MIN)
  {
    # pragma omp for
    for (sid = 0; sid < n_sblocks; sid++) {

      lmin = HUGE_VAL;
      lmax = -HUGE_VAL;
      s[0] = 0.0;
      s[1] = 0.0;

      for (bid = 0; bid < blocks_in_sblocks; bid++) {
        start_id = block_size * (blocks_in_sblocks*sid + bid);
        end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        c[0] = 0.0;
        c[1] = 0.0;
        for (i = start_id; i < end_id; i++) {
          c[0] += v[i];
          c[1] += v[i]*w[i];
          if (v[i] < lmin)
            lmin = v[i];
          if (v[i] > lmax)
            lmax = v[i];
        }
        s[0] += c[0];
        s[1] += c[1];
      }

#     pragma omp critical
      {
        if (lmin < *vmin)
          *vmin = lmin;
        if (lmax > *vmax)
          *vmax = lmax;
        *vsum += s[0];
        *wsum += s[1];
      }

    }
  } /* End of OpenMP-threaded section */

  start_id = block_size * n_sblocks*blocks_in_sblocks;
  end_id = n;
  c[0] = 0.0;
  c[1] = 0.0;
  for (i = start_id; i < end_id; i++) {
    c[0] += v[i];
    c[1] += v[i]*w[i];
    if (v[i] < *vmin)
      *vmin = v[i];
    if (v[i] > *vmax)
      *vmax = v[i];
  }

  *vsum += c[0];
  *wsum += c[1];
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum, weighted sum) of a
 * 1-dimensional array using an element list relative to weights.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n        <-- local number of elements
 *   wl       <-- pointer to element weights list
 *   v        <-- pointer to values (size: n)
 *   w        <-- pointer to elements weights
 *   vmin     --> resulting min
 *   vmax     --> resulting max
 *   vsum     --> resulting sum
 *   wsum     --> resulting weighted sum
 *----------------------------------------------------------------------------*/

static void
_cs_real_sstats_1d_w_l(cs_lnum_t         n,
                       const cs_lnum_t   wl[],
                       const cs_real_t   v[],
                       const cs_real_t   w[],
                       double           *vmin,
                       double           *vmax,
                       double           *vsum,
                       double           *wsum)
{
  const cs_lnum_t block_size = 60;

  cs_lnum_t i;
  cs_lnum_t sid, bid;
  cs_lnum_t start_id, end_id;
  double c[2], s[2], lmin, lmax;

  cs_lnum_t n_blocks = n / block_size;
  cs_lnum_t n_sblocks = sqrt(n_blocks);
  cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

#if defined(__xlc__)
#pragma disjoint(*wl, *v, *w, *vmin, *vmax, *vsum, *wsum)
#endif

  *vmin = HUGE_VAL;
  *vmax = -HUGE_VAL;
  *vsum = 0.;
  *wsum = 0.;

# pragma omp parallel private(bid, start_id, end_id, i, \
                              c, s, lmin, lmax) if (n > THR_MIN)
  {
    # pragma omp for
    for (sid = 0; sid < n_sblocks; sid++) {

      lmin = HUGE_VAL;
      lmax = -HUGE_VAL;
      s[0] = 0.0;
      s[1] = 0.0;

      for (bid = 0; bid < blocks_in_sblocks; bid++) {
        start_id = block_size * (blocks_in_sblocks*sid + bid);
        end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        c[0] = 0.0;
        c[1] = 0.0;
        for (i = start_id; i < end_id; i++) {
          c[0] += v[i];
          c[1] += v[i]*w[wl[i]];
          if (v[i] < lmin)
            lmin = v[i];
          if (v[i] > lmax)
            lmax = v[i];
        }
        s[0] += c[0];
        s[1] += c[1];
      }

#     pragma omp critical
      {
        if (lmin < *vmin)
          *vmin = lmin;
        if (lmax > *vmax)
          *vmax = lmax;
        *vsum += s[0];
        *wsum += s[1];
      }

    }
  } /* End of OpenMP-threaded section */

  start_id = block_size * n_sblocks*blocks_in_sblocks;
  end_id = n;
  c[0] = 0.0;
  c[1] = 0.0;
  for (i = start_id; i < end_id; i++) {
    c[0] += v[i];
    c[1] += v[i]*w[wl[i]];
    if (v[i] < *vmin)
      *vmin = v[i];
    if (v[i] > *vmax)
      *vmax = v[i];
  }

  *vsum += c[0];
  *wsum += c[1];
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum, weighted sum) of a
 * subset of a 1-dimensional array.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n        <-- local number of elements
 *   vl       <-- pointer to element list
 *   v        <-- pointer to field values (size: n)
 *   w        <-- pointer to weight values (size: n)
 *   vmin     --> resulting min
 *   vmax     --> resulting max
 *   vsum     --> resulting sum
 *   wsum     --> resulting weighted sum
 *----------------------------------------------------------------------------*/

static void
_cs_real_sstats_1d_l_w(cs_lnum_t         n,
                       const cs_lnum_t   vl[],
                       const cs_real_t   v[],
                       const cs_real_t   w[],
                       double           *vmin,
                       double           *vmax,
                       double           *vsum,
                       double           *wsum)
{
  const cs_lnum_t block_size = 60;

  cs_lnum_t i, li;
  cs_lnum_t sid, bid;
  cs_lnum_t start_id, end_id;
  double c[2], s[2], lmin, lmax;

  cs_lnum_t n_blocks = n / block_size;
  cs_lnum_t n_sblocks = sqrt(n_blocks);
  cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

#if defined(__xlc__)
#pragma disjoint(*vl, *v, *w, *vmin, *vmax, *vsum, *wsum)
#endif

  *vmin = HUGE_VAL;
  *vmax = -HUGE_VAL;
  *vsum = 0.;
  *wsum = 0.;

# pragma omp parallel private(bid, start_id, end_id, li, i, \
                              c, s, lmin, lmax) if (n > THR_MIN)
  {
    # pragma omp for
    for (sid = 0; sid < n_sblocks; sid++) {

      lmin = HUGE_VAL;
      lmax = -HUGE_VAL;
      s[0] = 0.0;
      s[1] = 0.0;

      for (bid = 0; bid < blocks_in_sblocks; bid++) {
        start_id = block_size * (blocks_in_sblocks*sid + bid);
        end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        c[0] = 0.0;
        c[1] = 0.0;
        for (li = start_id; li < end_id; li++) {
          i = vl[li];
          c[0] += v[i];
          c[1] += v[i]*w[i];
          if (v[i] < lmin)
            lmin = v[i];
          if (v[i] > lmax)
            lmax = v[i];
        }
        s[0] += c[0];
        s[1] += c[1];
      }

#     pragma omp critical
      {
        if (lmin < *vmin)
          *vmin = lmin;
        if (lmax > *vmax)
          *vmax = lmax;
        *vsum += s[0];
        *wsum += s[1];
      }

    }
  } /* End of OpenMP-threaded section */

  start_id = block_size * n_sblocks*blocks_in_sblocks;
  end_id = n;
  c[0] = 0.0;
  c[1] = 0.0;
  for (li = start_id; li < end_id; li++) {
    i = vl[li];
    c[0] += v[i];
    c[1] += v[i]*w[i];
    if (v[i] < *vmin)
      *vmin = v[i];
    if (v[i] > *vmax)
      *vmax = v[i];
  }

  *vsum += c[0];
  *wsum += c[1];
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum) of a 3d
 * array's components and norm.
 *
 * The array is interleaved.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n        <-- local number of elements
 *   v        <-- pointer to field values (size: n*3)
 *   vmin     --> resulting min array (size: 4)
 *   vmax     --> resulting max array (size: 4)
 *   vsum     --> resulting sum array (size: 4)
 *----------------------------------------------------------------------------*/

static void
_cs_real_sstats_3d(cs_lnum_t          n,
                   const cs_real_3_t  v[],
                   double             vmin[4],
                   double             vmax[4],
                   double             vsum[4])
{
  const cs_lnum_t block_size = 60;

  int j;
  cs_lnum_t i, sid, bid;
  cs_lnum_t start_id, end_id;
  double v_norm;
  double c[4], s[4], lmin[4], lmax[4];

  cs_lnum_t n_blocks = n / block_size;
  cs_lnum_t n_sblocks = sqrt(n_blocks);
  cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

#if defined(__xlc__)
#pragma disjoint(*v, *vmin, *vmax, *vsum)
#endif

  for (j = 0; j < 4; j++) {
    vmin[j] = HUGE_VAL;
    vmax[j] = -HUGE_VAL;
    vsum[j] = 0.;
  }

# pragma omp parallel private(bid, start_id, end_id, i, j, \
                              c, s, v_norm, lmin, lmax) if (n > THR_MIN)
  {
    # pragma omp for
    for (sid = 0; sid < n_sblocks; sid++) {

      for (j = 0; j < 4; j++) {
        lmin[j] = HUGE_VAL;
        lmax[j] = -HUGE_VAL;
      }
      for (j = 0; j < 4; j++)
        s[j] = 0.0;

      for (bid = 0; bid < blocks_in_sblocks; bid++) {
        start_id = block_size * (blocks_in_sblocks*sid + bid);
        end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        for (j = 0; j < 4; j++)
          c[j] = 0.0;
        for (i = start_id; i < end_id; i++) {
          for (j = 0; j < 3; j++) {
            c[j]   += v[i][j];
            if (v[i][j] < lmin[j])
              lmin[j] = v[i][j];
            if (v[i][j] > lmax[j])
              lmax[j] = v[i][j];
          }
          v_norm = sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
          c[3] += v_norm;
          if (v_norm < lmin[3])
            lmin[3] = v_norm;
          if (v_norm > lmax[3])
            lmax[3] = v_norm;
        }

        for (j = 0; j < 4; j++)
          s[j] += c[j];
      }

#     pragma omp critical
      {
        for (j = 0; j < 4; j++) {
          if (lmin[j] < vmin[j])
            vmin[j] = lmin[j];
          if (lmax[j] > vmax[j])
            vmax[j] = lmax[j];
          vsum[j] += s[j];
        }
      }
    }
  } /* End of OpenMP-threaded section */

  for (j = 0; j < 4; j++)
    c[j] = 0.0;
  start_id = block_size * n_sblocks*blocks_in_sblocks;
  end_id = n;
  for (i = start_id; i < end_id; i++) {
    for (j = 0; j < 3; j++) {
      c[j]   += v[i][j];
      if (v[i][j] < vmin[j])
        vmin[j] = v[i][j];
      if (v[i][j] > vmax[j])
        vmax[j] = v[i][j];
    }
    v_norm = sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
    c[3] += v_norm;
    if (v_norm < vmin[3])
      vmin[3] = v_norm;
    if (v_norm > vmax[3])
      vmax[3] = v_norm;
  }

  for (j = 0; j < 4; j++)
    vsum[j] += c[j];
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum) of a subset of a 3d
 * array's components and norm.
 *
 * The array is interleaved.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n        <-- local number of elements
 *   vl       <-- pointer to element list
 *   v        <-- pointer to field values (size: n*3)
 *   vmin     --> resulting min array (size: 4)
 *   vmax     --> resulting max array (size: 4)
 *   vsum     --> resulting sum array (size: 4)
 *----------------------------------------------------------------------------*/

static void
_cs_real_sstats_3d_l(cs_lnum_t          n,
                     const cs_lnum_t    vl[],
                     const cs_real_3_t  v[],
                     double             vmin[4],
                     double             vmax[4],
                     double             vsum[4])
{
  const cs_lnum_t block_size = 60;

  int j;
  cs_lnum_t i, li, sid, bid;
  cs_lnum_t start_id, end_id;
  double v_norm;
  double c[4], s[4], lmin[4], lmax[4];

  cs_lnum_t n_blocks = n / block_size;
  cs_lnum_t n_sblocks = sqrt(n_blocks);
  cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

#if defined(__xlc__)
#pragma disjoint(*v, *vmin, *vmax, *vsum)
#endif

  for (j = 0; j < 4; j++) {
    vmin[j] = HUGE_VAL;
    vmax[j] = -HUGE_VAL;
    vsum[j] = 0.;
  }

# pragma omp parallel private(bid, start_id, end_id, li, i, j, \
                              c, s, v_norm, lmin, lmax) if (n > THR_MIN)
  {
    # pragma omp for
    for (sid = 0; sid < n_sblocks; sid++) {

      for (j = 0; j < 4; j++) {
        lmin[j] = HUGE_VAL;
        lmax[j] = -HUGE_VAL;
      }
      for (j = 0; j < 4; j++)
        s[j] = 0.0;

      for (bid = 0; bid < blocks_in_sblocks; bid++) {
        start_id = block_size * (blocks_in_sblocks*sid + bid);
        end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        for (j = 0; j < 4; j++)
          c[j] = 0.0;
        for (li = start_id; li < end_id; li++) {
          i = vl[li];
          for (j = 0; j < 3; j++) {
            c[j]   += v[i][j];
            if (v[i][j] < lmin[j])
              lmin[j] = v[i][j];
            if (v[i][j] > lmax[j])
              lmax[j] = v[i][j];
          }
          v_norm = sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
          c[3] += v_norm;
          if (v_norm < lmin[3])
            lmin[3] = v_norm;
          if (v_norm > lmax[3])
            lmax[3] = v_norm;
        }

        for (j = 0; j < 4; j++)
          s[j] += c[j];
      }

#     pragma omp critical
      {
        for (j = 0; j < 4; j++) {
          if (lmin[j] < vmin[j])
            vmin[j] = lmin[j];
          if (lmax[j] > vmax[j])
            vmax[j] = lmax[j];
          vsum[j] += s[j];
        }
      }
    }
  } /* End of OpenMP-threaded section */

  for (j = 0; j < 4; j++)
    c[j] = 0.0;
  start_id = block_size * n_sblocks*blocks_in_sblocks;
  end_id = n;
  for (li = start_id; li < end_id; li++) {
    i = vl[li];
    for (j = 0; j < 3; j++) {
      c[j]   += v[i][j];
      if (v[i][j] < vmin[j])
        vmin[j] = v[i][j];
      if (v[i][j] > vmax[j])
        vmax[j] = v[i][j];
    }
    v_norm = sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
    c[3] += v_norm;
    if (v_norm < vmin[3])
      vmin[3] = v_norm;
    if (v_norm > vmax[3])
      vmax[3] = v_norm;
  }

  for (j = 0; j < 4; j++)
    vsum[j] += c[j];
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum, weighted sum) of a 3d
 * field array's components and norm.
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
 *----------------------------------------------------------------------------*/

static void
_cs_real_sstats_3d_w(cs_lnum_t          n,
                     const cs_real_3_t  v[],
                     const double       w[],
                     double             vmin[4],
                     double             vmax[4],
                     double             vsum[4],
                     double             wsum[4])
{
  const cs_lnum_t block_size = 60;

  int j;
  cs_lnum_t i, sid, bid;
  cs_lnum_t start_id, end_id;
  double v_norm;
  double c[8], s[8], lmin[4], lmax[4];

  cs_lnum_t n_blocks = n / block_size;
  cs_lnum_t n_sblocks = sqrt(n_blocks);
  cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

#if defined(__xlc__)
#pragma disjoint(*v, *w, *vmin, *vmax, *vsum, *wsum)
#endif

  for (j = 0; j < 4; j++) {
    vmin[j] = HUGE_VAL;
    vmax[j] = -HUGE_VAL;
    vsum[j] = 0.;
    wsum[j] = 0.;
  }

# pragma omp parallel private(bid, start_id, end_id, i, j, \
                              c, s, v_norm, lmin, lmax) if (n > THR_MIN)
  {
    # pragma omp for
    for (sid = 0; sid < n_sblocks; sid++) {

      for (j = 0; j < 4; j++) {
        lmin[j] = HUGE_VAL;
        lmax[j] = -HUGE_VAL;
      }
      for (j = 0; j < 8; j++)
        s[j] = 0.0;

      for (bid = 0; bid < blocks_in_sblocks; bid++) {
        start_id = block_size * (blocks_in_sblocks*sid + bid);
        end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        for (j = 0; j < 8; j++)
          c[j] = 0.0;
        for (i = start_id; i < end_id; i++) {
          for (j = 0; j < 3; j++) {
            c[j]   += v[i][j];
            c[j+4] += v[i][j]*w[i];
            if (v[i][j] < lmin[j])
              lmin[j] = v[i][j];
            if (v[i][j] > lmax[j])
              lmax[j] = v[i][j];
          }
          v_norm = sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
          c[3] += v_norm;
          c[7] += v_norm*w[i];
          if (v_norm < lmin[3])
            lmin[3] = v_norm;
          if (v_norm > lmax[3])
            lmax[3] = v_norm;
        }

        for (j = 0; j < 8; j++)
          s[j] += c[j];
      }

#     pragma omp critical
      {
        for (j = 0; j < 4; j++) {
          if (lmin[j] < vmin[j])
            vmin[j] = lmin[j];
          if (lmax[j] > vmax[j])
            vmax[j] = lmax[j];
          vsum[j] += s[j];
          wsum[j] += s[4+j];
        }
      }

    }
  } /* End of OpenMP-threaded section */

  for (j = 0; j < 8; j++)
    c[j] = 0.0;
  start_id = block_size * n_sblocks*blocks_in_sblocks;
  end_id = n;
  for (i = start_id; i < end_id; i++) {
    for (j = 0; j < 3; j++) {
      c[j]   += v[i][j];
      c[j+4] += v[i][j]*w[i];
      if (v[i][j] < vmin[j])
        vmin[j] = v[i][j];
      if (v[i][j] > vmax[j])
        vmax[j] = v[i][j];
    }
    v_norm = sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
    c[3] += v_norm;
    c[7] += v_norm*w[i];
    if (v_norm < vmin[3])
      vmin[3] = v_norm;
    if (v_norm > vmax[3])
      vmax[3] = v_norm;
  }

  for (j = 0; j < 4; j++) {
    vsum[j] += c[j];
    wsum[j] += c[4+j];
  }
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum, weighted sum) of a 3d
 * field array's components and norm  using an element list relative
 * to weights.
 *
 * The array is interleaved.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n        <-- local number of elements
 *   wl       <-- pointer to weights list
 *   v        <-- pointer to element values (size: n*3)
 *   w        <-- pointer to weight values
 *   vmin     --> resulting min array (size: 4)
 *   vmax     --> resulting max array (size: 4)
 *   vsum     --> resulting sum array (size: 4)
 *   wsum     --> resulting weighted sum array (size: 4)
 *----------------------------------------------------------------------------*/

static void
_cs_real_sstats_3d_w_l(cs_lnum_t          n,
                       const cs_lnum_t    wl[],
                       const cs_real_3_t  v[],
                       const double       w[],
                       double             vmin[4],
                       double             vmax[4],
                       double             vsum[4],
                       double             wsum[4])
{
  const cs_lnum_t block_size = 60;

  int j;
  cs_lnum_t i, sid, bid;
  cs_lnum_t start_id, end_id;
  double v_norm, wi;
  double c[8], s[8], lmin[4], lmax[4];

  cs_lnum_t n_blocks = n / block_size;
  cs_lnum_t n_sblocks = sqrt(n_blocks);
  cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

#if defined(__xlc__)
#pragma disjoint(*wl, *v, *w, *vmin, *vmax, *vsum, *wsum)
#endif

  for (j = 0; j < 4; j++) {
    vmin[j] = HUGE_VAL;
    vmax[j] = -HUGE_VAL;
    vsum[j] = 0.;
    wsum[j] = 0.;
  }

# pragma omp parallel private(bid, start_id, end_id, i, j, \
                              wi, c, s, v_norm, lmin, lmax) if (n > THR_MIN)
  {
    # pragma omp for
    for (sid = 0; sid < n_sblocks; sid++) {

      for (j = 0; j < 4; j++) {
        lmin[j] = HUGE_VAL;
        lmax[j] = -HUGE_VAL;
      }
      for (j = 0; j < 8; j++)
        s[j] = 0.0;

      for (bid = 0; bid < blocks_in_sblocks; bid++) {
        start_id = block_size * (blocks_in_sblocks*sid + bid);
        end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        for (j = 0; j < 8; j++)
          c[j] = 0.0;
        for (i = start_id; i < end_id; i++) {
          wi = w[wl[i]];
          for (j = 0; j < 3; j++) {
            c[j]   += v[i][j];
            c[j+4] += v[i][j]*wi;
            if (v[i][j] < lmin[j])
              lmin[j] = v[i][j];
            if (v[i][j] > lmax[j])
              lmax[j] = v[i][j];
          }
          v_norm = sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
          c[3] += v_norm;
          c[7] += v_norm*wi;
          if (v_norm < lmin[3])
            lmin[3] = v_norm;
          if (v_norm > lmax[3])
            lmax[3] = v_norm;
        }

        for (j = 0; j < 8; j++)
          s[j] += c[j];
      }

#     pragma omp critical
      {
        for (j = 0; j < 4; j++) {
          if (lmin[j] < vmin[j])
            vmin[j] = lmin[j];
          if (lmax[j] > vmax[j])
            vmax[j] = lmax[j];
          vsum[j] += s[j];
          wsum[j] += s[4+j];
        }
      }
    }
  } /* End of OpenMP-threaded section */

  for (j = 0; j < 8; j++)
    c[j] = 0.0;
  start_id = block_size * n_sblocks*blocks_in_sblocks;
  end_id = n;
  for (i = start_id; i < end_id; i++) {
    wi = w[wl[i]];
    for (j = 0; j < 3; j++) {
      c[j]   += v[i][j];
      c[j+4] += v[i][j]*wi;
      if (v[i][j] < vmin[j])
        vmin[j] = v[i][j];
      if (v[i][j] > vmax[j])
        vmax[j] = v[i][j];
    }
    v_norm = sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
    c[3] += v_norm;
    c[7] += v_norm*wi;
    if (v_norm < vmin[3])
      vmin[3] = v_norm;
    if (v_norm > vmax[3])
      vmax[3] = v_norm;
  }

  for (j = 0; j < 4; j++) {
    vsum[j] += c[j];
    wsum[j] += c[4+j];
  }
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum, weighted sum) of a
 * subset of a 3d field array's components and norm.
 *
 * The array is interleaved.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n        <-- local number of elements
 *   vl       <-- pointer to element list
 *   v        <-- pointer to element values (size: n*3)
 *   w        <-- pointer to weight values (size: n)
 *   vmin     --> resulting min array (size: 4)
 *   vmax     --> resulting max array (size: 4)
 *   vsum     --> resulting sum array (size: 4)
 *   wsum     --> resulting weighted sum array (size: 4)
 *----------------------------------------------------------------------------*/

static void
_cs_real_sstats_3d_l_w(cs_lnum_t          n,
                       const cs_lnum_t    vl[],
                       const cs_real_3_t  v[],
                       const double       w[],
                       double             vmin[4],
                       double             vmax[4],
                       double             vsum[4],
                       double             wsum[4])
{
  const cs_lnum_t block_size = 60;

  int j;
  cs_lnum_t li, i, sid, bid;
  cs_lnum_t start_id, end_id;
  double v_norm;
  double c[8], s[8], lmin[4], lmax[4];

  cs_lnum_t n_blocks = n / block_size;
  cs_lnum_t n_sblocks = sqrt(n_blocks);
  cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

#if defined(__xlc__)
#pragma disjoint(*vl, *v, *w, *vmin, *vmax, *vsum, *wsum)
#endif

  for (j = 0; j < 4; j++) {
    vmin[j] = HUGE_VAL;
    vmax[j] = -HUGE_VAL;
    vsum[j] = 0.;
    wsum[j] = 0.;
  }

# pragma omp parallel private(bid, start_id, end_id, li, i, j, \
                              c, s, v_norm, lmin, lmax) if (n > THR_MIN)
  {
    # pragma omp for
    for (sid = 0; sid < n_sblocks; sid++) {

      for (j = 0; j < 4; j++) {
        lmin[j] = HUGE_VAL;
        lmax[j] = -HUGE_VAL;
      }
      for (j = 0; j < 8; j++)
        s[j] = 0.0;

      for (bid = 0; bid < blocks_in_sblocks; bid++) {
        start_id = block_size * (blocks_in_sblocks*sid + bid);
        end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        for (j = 0; j < 8; j++)
          c[j] = 0.0;
        for (li = start_id; li < end_id; li++) {
          i = vl[li];
          for (j = 0; j < 3; j++) {
            c[j]   += v[i][j];
            c[j+4] += v[i][j]*w[i];
            if (v[i][j] < lmin[j])
              lmin[j] = v[i][j];
            if (v[i][j] > lmax[j])
              lmax[j] = v[i][j];
          }
          v_norm = sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
          c[3] += v_norm;
          c[7] += v_norm*w[i];
          if (v_norm < lmin[3])
            lmin[3] = v_norm;
          if (v_norm > lmax[3])
            lmax[3] = v_norm;
        }

        for (j = 0; j < 8; j++)
          s[j] += c[j];
      }

#     pragma omp critical
      {
        for (j = 0; j < 4; j++) {
          if (lmin[j] < vmin[j])
            vmin[j] = lmin[j];
          if (lmax[j] > vmax[j])
            vmax[j] = lmax[j];
          vsum[j] += s[j];
          wsum[j] += s[4+j];
        }
      }
    }
  } /* End of OpenMP-threaded section */

  for (j = 0; j < 8; j++)
    c[j] = 0.0;
  start_id = block_size * n_sblocks*blocks_in_sblocks;
  end_id = n;
  for (li = start_id; li < end_id; li++) {
    i = vl[li];
    for (j = 0; j < 3; j++) {
      c[j]   += v[i][j];
      c[j+4] += v[i][j]*w[i];
      if (v[i][j] < vmin[j])
        vmin[j] = v[i][j];
      if (v[i][j] > vmax[j])
        vmax[j] = v[i][j];
    }
    v_norm = sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
    c[3] += v_norm;
    c[7] += v_norm*w[i];
    if (v_norm < vmin[3])
      vmin[3] = v_norm;
    if (v_norm > vmax[3])
      vmax[3] = v_norm;
  }

  for (j = 0; j < 4; j++) {
    vsum[j] += c[j];
    wsum[j] += c[4+j];
  }
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum) of an
 * n-dimensional field array's components.
 *
 * The maximum allowed dimension is 9 (allowing for a rank-2 tensor).
 * The array is interleaved.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n        <-- local number of elements
 *   dim      <-- local array dimension (max: 9)
 *   vl       <-- pointer to optional element list, or NULL
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
  const cs_lnum_t block_size = 60;

  int j;
  cs_lnum_t li, i, sid, bid;
  cs_lnum_t start_id, end_id;
  double c[9], s[9], lmin[9], lmax[9];

  cs_lnum_t n_blocks = n / block_size;
  cs_lnum_t n_sblocks = sqrt(n_blocks);
  cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

  assert(dim <= 9);

#if defined(__xlc__)
#pragma disjoint(*vl, *v, *vmin, *vmax, *vsum)
#endif

  for (j = 0; j < dim; j++) {
    vmin[j] = HUGE_VAL;
    vmax[j] = -HUGE_VAL;
    vsum[j] = 0.;
  }

# pragma omp parallel private(bid, start_id, end_id, li, i, j, \
                              c, s, lmin, lmax) if (n > THR_MIN)
  {
    # pragma omp for
    for (sid = 0; sid < n_sblocks; sid++) {

      for (j = 0; j < dim; j++) {
        lmin[j] = HUGE_VAL;
        lmax[j] = -HUGE_VAL;
      }
      for (j = 0; j < dim; j++)
        s[j] = 0.0;

      if (vl == NULL) {

        for (bid = 0; bid < blocks_in_sblocks; bid++) {
          start_id = block_size * (blocks_in_sblocks*sid + bid);
          end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
          for (j = 0; j < dim; j++)
            c[j] = 0.0;
          for (i = start_id; i < end_id; i++) {
            for (j = 0; j < dim; j++) {
              c[j] += v[i*dim + j];
              if (v[i*dim + j] < lmin[j])
                lmin[j] = v[i*dim + j];
              if (v[i*dim+j] > lmax[j])
                lmax[j] = v[i*dim+j];
            }
          }
          for (j = 0; j < dim; j++)
            s[j] += c[j];
        }

      }
      else {

        for (bid = 0; bid < blocks_in_sblocks; bid++) {
          start_id = block_size * (blocks_in_sblocks*sid + bid);
          end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
          for (j = 0; j < dim; j++)
            c[j] = 0.0;
          for (li = start_id; li < end_id; li++) {
            i = vl[li];
            for (j = 0; j < dim; j++) {
              c[j] += v[i*dim + j];
              if (v[i*dim + j] < lmin[j])
                lmin[j] = v[i*dim + j];
              if (v[i*dim+j] > lmax[j])
                lmax[j] = v[i*dim+j];
            }
          }
          for (j = 0; j < dim; j++)
            s[j] += c[j];
        }

      }

#     pragma omp critical
      {
        for (j = 0; j < dim; j++) {
          if (lmin[j] < vmin[j])
            vmin[j] = lmin[j];
          if (lmax[j] > vmax[j])
            vmax[j] = lmax[j];
          vsum[j] += s[j];
        }
      }

    }
  } /* End of OpenMP-threaded section */

  for (j = 0; j < dim; j++)
    c[j] = 0.0;
  start_id = block_size * n_sblocks*blocks_in_sblocks;
  end_id = n;

  if (vl == NULL) {

    for (i = start_id; i < end_id; i++) {
      for (j = 0; j < dim; j++) {
        c[j] += v[i*dim + j];
        if (v[i*dim + j] < lmin[j])
          vmin[j] = v[i*dim + j];
        if (v[i*dim+j] > lmax[j])
          vmax[j] = v[i*dim+j];
      }
    }

  }
  else {

    for (li = start_id; li < end_id; li++) {
      i = vl[i];
      for (j = 0; j < dim; j++) {
        c[j] += v[i*dim + j];
        if (v[i*dim + j] < lmin[j])
          vmin[j] = v[i*dim + j];
        if (v[i*dim+j] > lmax[j])
          vmax[j] = v[i*dim+j];
      }
    }

  }

  for (j = 0; j < dim; j++)
    vsum[j] += c[j];
}

/*----------------------------------------------------------------------------
 * Compute simple local stats (minima, maxima, sum, weighted sum) of an
 * n-dimensional field array's components.
 *
 * The maximum allowed dimension is 9 (allowing for a rank-2 tensor).
 * The array is interleaved.
 *
 * The algorithm here is similar to that used for blas, but computes several
 * quantities simultaneously for better cache behavior
 *
 * parameters:
 *   n        <-- local number of elements
 *   dim      <-- local array dimension (max: 9)
 *   vl       <-- pointer to optional element list, or NULL
 *   wl       <-- pointer to optional element list, or NULL
 *                (ignored if vl != NULL)
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
                     const double      w[],
                     double            vmin[],
                     double            vmax[],
                     double            vsum[],
                     double            wsum[])
{
  const cs_lnum_t block_size = 60;

  int j;
  cs_lnum_t li, i, sid, bid;
  cs_lnum_t start_id, end_id;
  double wi, c[18], s[18], lmin[9], lmax[9];

  cs_lnum_t n_blocks = n / block_size;
  cs_lnum_t n_sblocks = sqrt(n_blocks);
  cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

  const int dim2 = dim*2;

  assert(dim <= 9);

#if defined(__xlc__)
#pragma disjoint(*vl, *wl, *v, *w, *vmin, *vmax, *vsum, *wsum)
#endif

  for (j = 0; j < dim; j++) {
    vmin[j] = HUGE_VAL;
    vmax[j] = -HUGE_VAL;
    vsum[j] = 0.;
    wsum[j] = 0.;
  }

# pragma omp parallel private(bid, start_id, end_id, li, i, j, \
                              wi, c, s, lmin, lmax) if (n > THR_MIN)
  {
    # pragma omp for
    for (sid = 0; sid < n_sblocks; sid++) {

      for (j = 0; j < dim; j++) {
        lmin[j] = HUGE_VAL;
        lmax[j] = -HUGE_VAL;
      }
      for (j = 0; j < dim2; j++)
        s[j] = 0.0;

      if (vl == NULL && wl == NULL) {
        for (bid = 0; bid < blocks_in_sblocks; bid++) {
          start_id = block_size * (blocks_in_sblocks*sid + bid);
          end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
          for (j = 0; j < dim2; j++)
            c[j] = 0.0;
          for (i = start_id; i < end_id; i++) {
            for (j = 0; j < dim; j++) {
              c[j]     += v[i*dim + j];
              c[j+dim] += v[i*dim + j]*w[i];
              if (v[i*dim + j] < lmin[j])
                lmin[j] = v[i*dim + j];
              if (v[i*dim+j] > lmax[j])
                lmax[j] = v[i*dim+j];
            }
          }
          for (j = 0; j < dim2; j++)
            s[j] += c[j];
        }
      }
      else if (vl == NULL) {
        for (bid = 0; bid < blocks_in_sblocks; bid++) {
          start_id = block_size * (blocks_in_sblocks*sid + bid);
          end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
          for (j = 0; j < dim2; j++)
            c[j] = 0.0;
          for (i = start_id; i < end_id; i++) {
            wi = wl[i];
            for (j = 0; j < dim; j++) {
              c[j]     += v[i*dim + j];
              c[j+dim] += v[i*dim + j]*wi;
              if (v[i*dim + j] < lmin[j])
                lmin[j] = v[i*dim + j];
              if (v[i*dim+j] > lmax[j])
                lmax[j] = v[i*dim+j];
            }
          }
          for (j = 0; j < dim2; j++)
            s[j] += c[j];
        }
      }
      else { /* vl != NULL */
        for (bid = 0; bid < blocks_in_sblocks; bid++) {
          start_id = block_size * (blocks_in_sblocks*sid + bid);
          end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
          for (j = 0; j < dim2; j++)
            c[j] = 0.0;
          for (li = start_id; li < end_id; li++) {
            i = vl[i];
            for (j = 0; j < dim; j++) {
              c[j]     += v[i*dim + j];
              c[j+dim] += v[i*dim + j]*w[i];
              if (v[i*dim + j] < lmin[j])
                lmin[j] = v[i*dim + j];
              if (v[i*dim+j] > lmax[j])
                lmax[j] = v[i*dim+j];
            }
          }
          for (j = 0; j < dim2; j++)
            s[j] += c[j];
        }
      }

#     pragma omp critical
      {
        for (j = 0; j < dim; j++) {
          if (lmin[j] < vmin[j])
            vmin[j] = lmin[j];
          if (lmax[j] > vmax[j])
            vmax[j] = lmax[j];
          vsum[j] += s[j];
          wsum[j] += s[dim+j];
        }
      }

    }
  } /* End of OpenMP-threaded section */

  for (j = 0; j < dim2; j++)
    c[j] = 0.0;
  start_id = block_size * n_sblocks*blocks_in_sblocks;
  end_id = n;

  if (vl == NULL && wl == NULL) {
    for (i = start_id; i < end_id; i++) {
      for (j = 0; j < dim; j++) {
        c[j]     += v[i*dim + j];
        c[j+dim] += v[i*dim + j]*w[i];
        if (v[i*dim + j] < lmin[j])
          vmin[j] = v[i*dim + j];
        if (v[i*dim+j] > lmax[j])
          vmax[j] = v[i*dim+j];
      }
    }
  }
  else if (vl == NULL) {
    for (i = start_id; i < end_id; i++) {
      wi = wl[i];
      for (j = 0; j < dim; j++) {
        c[j]     += v[i*dim + j];
        c[j+dim] += v[i*dim + j]*wi;
        if (v[i*dim + j] < lmin[j])
          vmin[j] = v[i*dim + j];
        if (v[i*dim+j] > lmax[j])
          vmax[j] = v[i*dim+j];
      }
    }
  }
  else { /* vl != NULL */
    for (li = start_id; li < end_id; li++) {
      i = vl[i];
      for (j = 0; j < dim; j++) {
        c[j]     += v[i*dim + j];
        c[j+dim] += v[i*dim + j]*w[i];
        if (v[i*dim + j] < lmin[j])
          vmin[j] = v[i*dim + j];
        if (v[i*dim+j] > lmax[j])
          vmax[j] = v[i*dim+j];
      }
    }
  }

  for (j = 0; j < dim; j++) {
    vsum[j] += c[j];
    wsum[j] += c[dim+j];
  }
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute sums of an n-dimensional cs_real_t array's components.
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
 *                          are defined, or NULL
 * \param[in]   v           pointer to array values
 * \param[out]  vsum        resulting sum array (size: dim, or 4 if dim = 3)
 */
/*----------------------------------------------------------------------------*/

void
cs_array_reduce_sum_l(cs_lnum_t         n_elts,
                      int               dim,
                      const cs_lnum_t  *v_elt_list,
                      const cs_real_t   v[],
                      double            vsum[])
{
  /* If all values are defined on same list */

  if (v_elt_list == NULL) {
    if (dim == 1)
      vsum[0] = _cs_real_sum_1d(n_elts, v);
    else if (dim == 3)
      bft_error(__FILE__, __LINE__, 0,
                _("_cs_real_sum_3d not implemented yet\n"));
    else
      bft_error(__FILE__, __LINE__, 0,
                _("_cs_real_sum_nd not implemented yet\n"));
  }

  /* If values are defined on parent list */

  else {
    if (dim == 1)
      bft_error(__FILE__, __LINE__, 0,
                _("_cs_real_sum_1d_l not implemented yet\n"));
    else if (dim == 3)
      bft_error(__FILE__, __LINE__, 0,
                _("_cs_real_sum_3d_l not implemented yet\n"));
    else
      bft_error(__FILE__, __LINE__, 0,
                _("_cs_real_sum_nd_l not implemented yet\n"));
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute simple local stats (minima, maxima, sum) of an
 * n-dimensional cs_real_t array's components.
 *
 * The maximum allowed dimension is 9 (allowing for a rank-2 tensor).
 * The array is interleaved.
 *
 * For arrays of dimension 3, the statistics relative to the norm
 * are also computed, and added at the end of the statistics arrays
 * (which must be size dim+1).
 *
 * The algorithm here is similar to that used for BLAS, but computes several
 * quantities simultaneously for better cache behavior
 *
 * \param[in]   n_elts      number of local elements
 * \param[in]   dim         local array dimension (max: 9)
 * \param[in]   v_elt_list  optional list of parent elements on which values
 *                          are defined, or NULL
 * \param[in]   v           pointer to array values
 * \param[out]  vmin        resulting min array (size: dim, or 4 if dim = 3)
 * \param[out]  vmax        resulting max array (size: dim, or 4 if dim = 3)
 * \param[out]  vsum        resulting sum array (size: dim, or 4 if dim = 3)
 */
/*----------------------------------------------------------------------------*/

void
cs_array_reduce_simple_stats_l(cs_lnum_t         n_elts,
                               int               dim,
                               const cs_lnum_t  *v_elt_list,
                               const cs_real_t   v[],
                               double            vmin[],
                               double            vmax[],
                               double            vsum[])
{
  /* If all values are defined on same list */

  if (v_elt_list == NULL) {
    if (dim == 1)
      _cs_real_sstats_1d(n_elts, v, vmin, vmax, vsum);
    else if (dim == 3)
      _cs_real_sstats_3d(n_elts, (const cs_real_3_t *)v,
                         vmin, vmax, vsum);
    else
      _cs_real_sstats_nd(n_elts, dim, NULL, v, vmin, vmax, vsum);
  }

  /* If values are defined on parent list */

  else {
    if (dim == 1)
      _cs_real_sstats_1d_l(n_elts, v_elt_list, v, vmin, vmax, vsum);
    else if (dim == 3)
      _cs_real_sstats_3d_l(n_elts, v_elt_list, (const cs_real_3_t *)v,
                           vmin, vmax, vsum);
    else
      _cs_real_sstats_nd(n_elts, dim, v_elt_list, v, vmin, vmax, vsum);
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute simple local stats (minima, maxima, sum, weighted sum) of
 * an n-dimensional cs_real_t array's components for a given mesh location.
 *
 * The maximum allowed dimension is 9 (allowing for a rank-2 tensor).
 * The array is interleaved.
 *
 * For arrays of dimension 3, the statistics relative to the norm
 * are also computed, and added at the end of the statistics arrays
 * (which must be size dim+1).
 *
 * The algorithm here is similar to that used for BLAS, but computes several
 * quantities simultaneously for better cache behavior
 *
 * \param[in]   dim         local array dimension (max: 9)
 * \param[in]   v_elt_list  optional list of parent elements on which values
 *                          are defined, or NULL
 * \param[in]   w_elt_list  optional list of parent elements on which weights
 *                          are defined, or NULL; if v_elt_list is defined
 *                          (ie. non-NULL),then w_elt_list = v_elt_list is
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
cs_array_reduce_simple_stats_l_w(cs_lnum_t         n_elts,
                                 int               dim,
                                 const cs_lnum_t  *v_elt_list,
                                 const cs_lnum_t  *w_elt_list,
                                 const cs_real_t   v[],
                                 const cs_real_t   w[],
                                 double            vmin[],
                                 double            vmax[],
                                 double            vsum[],
                                 double            wsum[])
{
  /* If all values are defined on same list */

  if (v_elt_list == NULL && w_elt_list == NULL) {
    if (dim == 1)
      _cs_real_sstats_1d_w(n_elts, v, w, vmin, vmax, vsum, wsum);
    else if (dim == 3)
      _cs_real_sstats_3d_w(n_elts, (const cs_real_3_t *)v, w,
                           vmin, vmax, vsum, wsum);
    else
      _cs_real_sstats_nd_w(n_elts, dim, NULL, NULL, v, w,
                           vmin, vmax, vsum, wsum);
  }

  /* If weights are defined on parent list */

  else if (v_elt_list == NULL) { /* w_elt_list != NULL */
    if (dim == 1)
      _cs_real_sstats_1d_w_l(n_elts, w_elt_list, v, w, vmin, vmax, vsum, wsum);
    else if (dim == 3)
      _cs_real_sstats_3d_w_l(n_elts, w_elt_list, (const cs_real_3_t *)v, w,
                             vmin, vmax, vsum, wsum);
    else
      _cs_real_sstats_nd_w(n_elts, dim, NULL, w_elt_list, v, w,
                           vmin, vmax, vsum, wsum);
  }

  /* If weights are defined on parent list */

  else { /* v_elt_list != NULL */
    if (dim == 1)
      _cs_real_sstats_1d_l_w(n_elts, v_elt_list, v, w, vmin, vmax, vsum, wsum);
    else if (dim == 3)
      _cs_real_sstats_3d_l_w(n_elts, v_elt_list, (const cs_real_3_t *)v, w,
                             vmin, vmax, vsum, wsum);
    else
      _cs_real_sstats_nd_w(n_elts, dim, v_elt_list, NULL, v, w,
                           vmin, vmax, vsum, wsum);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
