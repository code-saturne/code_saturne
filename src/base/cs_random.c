/*============================================================================
 * Random number generation.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_random.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_random.c
        Random number generation.

  Based on the uniform, gaussian, and poisson random number generation code
  from netlib.org: lagged (-273,-607) Fibonacci; Box-Muller;
  by W.P. Petersen, IPS, ETH Zuerich.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

static struct {
  int      fill_1[1214];
  double  *buff;
  int      ptr;
  int      e_2;
} klotz0_1 = {{0}, NULL, 0, 0};

static struct {
  double   xbuff[1024];
  int      first;
  int      xptr;
  double   e_3;
} klotz1_1 = {{0}, 0, 0, 0.};

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Global variables
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Save static variables used by uniform number generator.
 *
 * \param[out]  save_block  saved state values
 */
/*----------------------------------------------------------------------------*/

static void
_random_uniform_save(cs_real_t  save_block[608])
{
  /* Saves blocks klotz0, containing seeds and
     pointer to position in seed block. The entire contents
     of klotz0 (pointer in buff, and buff) must be saved. */

  save_block[0] = (double) klotz0_1.ptr;
# pragma omp simd
  for (int i = 0; i < 607; ++i) {
    save_block[i + 1] = klotz0_1.buff[i];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Restore static variables used by uniform number generator.
 *
 * \param[in]  save_block  saved state values
 */
/*----------------------------------------------------------------------------*/

static void
_random_uniform_restore(cs_real_t  save_block[608])
{
  /* Restores block klotz0, containing seeds and pointer
     to position in seed block. The entire contents
     of klotz0 must be restored. */

  klotz0_1.ptr = (int) save_block[0];
#  pragma omp simd
  for (int i = 0; i < 607; ++i) {
    klotz0_1.buff[i] = save_block[i + 1];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Internal function for Box_Muller normal distribution.
 */
/*----------------------------------------------------------------------------*/

static void
_normal00(void)
{
  /* Local variables */
  double twopi, r1, r2, t1, t2;

  twopi = 6.2831853071795862;
  cs_random_uniform(1024, klotz1_1.xbuff);
#pragma omp simd
  for (int i = 0; i < 1023; i += 2) {
    r1 = twopi * klotz1_1.xbuff[i];
    t1 = cos(r1);
    t2 = sin(r1);
    r2 = sqrt(-2.*(log(1. - klotz1_1.xbuff[i+1])));
    klotz1_1.xbuff[i]   = t1 * r2;
    klotz1_1.xbuff[i+1] = t2 * r2;
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize random number generator.
 *
 * Generates initial seed buffer by linear congruential method.
 * Taken from Marsaglia, FSU report FSU-SCRI-87-50.
 *
 * \param[in]  seed  variable seed, with 0 < seed < 31328
 */
/*----------------------------------------------------------------------------*/

void
cs_random_seed(int  seed)
{
  int kl = 9373;
  int ij = 1802;

  /* Set pointer */

  klotz0_1.buff = (double *)(klotz0_1.fill_1);

  /* Seed should be > 0 and < 31328 */
  if (seed > 0)
    ij = seed % 31328;

  int i = ij / 177 % 177 + 2;
  int j = ij % 177 + 2;
  int k = kl / 169 % 178 + 1;
  int l = kl % 169;

  for (int ii = 0; ii < 607; ++ii) {
    double s = 0.;
    double t = .5;
    for (int jj = 1; jj <= 24; ++jj) {
      int m = i * j % 179 * k % 179;
      i = j;
      j = k;
      k = m;
      l = (l * 53 + 1) % 169;
      if (l * m % 64 >= 32) {
        s += t;
      }
      t *= (double).5;
    }
    klotz0_1.buff[ii] = s;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Uniform distribution random number generator.
 *
 * Portable lagged Fibonacci series uniform random number generator
 * with "lags" -273 und -607:
 * W.P. Petersen, IPS, ETH Zuerich, 19 Mar. 92
 *
 * \param[in]   n  number of values to compute
 * \param[out]  a  pseudo-random numbers following uniform distribution
 */
/*----------------------------------------------------------------------------*/

void
cs_random_uniform(cs_lnum_t  n,
                  cs_real_t  a[])
{
  int buffsz = 607;

  int bptr, aptr0, i, k, q, kptr;
  double t;
  int left, vl, qq, k273, k607;

  int aptr = 0;
  int nn = n;

 L1:

  if (nn <= 0)
    return;

  /* factor nn = q*607 + r */

  q = (nn - 1) / 607;
  left = buffsz - klotz0_1.ptr;

  if (q <= 1) {

    /* only one or fewer full segments */

    if (nn < left) {
      kptr = klotz0_1.ptr;
      for (i = 0; i < nn; ++i) {
        a[i + aptr] = klotz0_1.buff[kptr + i];
      }
      klotz0_1.ptr += nn;
      return;
    }
    else {
      kptr = klotz0_1.ptr;
#     pragma omp simd
      for (i = 0; i < left; ++i) {
        a[i + aptr] = klotz0_1.buff[kptr + i];
      }
      klotz0_1.ptr = 0;
      aptr += left;
      nn -= left;
      /*  buff -> buff case */
      vl = 273;
      k273 = 334;
      k607 = 0;
      for (k = 0; k < 3; ++k) {
#       pragma omp simd
        for (i = 0; i < vl; ++i) {
          t = klotz0_1.buff[k273+i]+klotz0_1.buff[k607+i];
          klotz0_1.buff[k607+i] = t - (double) ((int) t);
        }
        k607 += vl;
        k273 += vl;
        vl = 167;
        if (k == 0) {
          k273 = 0;
        }
      }
      goto L1;
    }
  } else {

    /* more than 1 full segment */

    kptr = klotz0_1.ptr;
#   pragma omp simd
    for (i = 0; i < left; ++i) {
      a[i + aptr] = klotz0_1.buff[kptr + i];
    }
    nn -= left;
    klotz0_1.ptr = 0;
    aptr += left;

/* buff -> a(aptr0) */

    vl = 273;
    k273 = 334;
    k607 = 0;
    for (k = 0; k < 3; ++k) {
      if (k == 0) {
#       pragma omp simd
        for (i = 0; i < vl; ++i) {
          t = klotz0_1.buff[k273+i]+klotz0_1.buff[k607+i];
          a[aptr + i] = t - (double) ((int) t);
        }
        k273 = aptr;
        k607 += vl;
        aptr += vl;
        vl = 167;
      } else {
#       pragma omp simd
        for (i = 0; i < vl; ++i) {
          t = a[k273 + i] + klotz0_1.buff[k607 + i];
          a[aptr + i] = t - (double) ((int) t);
        }
        k607 += vl;
        k273 += vl;
        aptr += vl;
      }
    }
    nn += -607;

    /* a(aptr-607) -> a(aptr) for last of the q-1 segments */

    aptr0 = aptr - 607;
    vl = 607;

    for (qq = 0; qq < q-2; ++qq) {
      k273 = aptr0 + 334;
#     pragma omp simd
      for (i = 0; i < vl; ++i) {
        t = a[k273 + i] + a[aptr0 + i];
        a[aptr + i] = t - (double) ((int) t);
      }
      nn += -607;
      aptr += vl;
      aptr0 += vl;
    }

    /* a(aptr0) -> buff, last segment before residual */

    vl = 273;
    k273 = aptr0 + 334;
    k607 = aptr0;
    bptr = 0;
    for (k = 0; k < 3; ++k) {
      if (k == 0) {
#       pragma omp simd
        for (i = 0; i < vl; ++i) {
          t = a[k273 + i] + a[k607 + i];
          klotz0_1.buff[bptr + i] = t - (double) ((int) t);
        }
        k273 = 0;
        k607 += vl;
        bptr += vl;
        vl = 167;
      } else {
#       pragma omp simd
        for (i = 0; i < vl; ++i) {
          t = klotz0_1.buff[k273 + i] + a[k607 + i];
          klotz0_1.buff[bptr + i] = t - (double) ((int) t);
        }
        k607 += vl;
        k273 += vl;
        bptr += vl;
      }
    }
    goto L1;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Normal distribution random number generator.
 *
 * Box-Muller method for Gaussian random numbers.
 *
 * \param[in]   n  number of values to compute
 * \param[out]  x  pseudo-random numbers following normal distribution
 */
/*----------------------------------------------------------------------------*/

void
cs_random_normal(cs_lnum_t  n,
                 cs_real_t  x[])
{
  int buffsz = 1024;

  /* Local variables */
  int left, i, nn, ptr, kptr;

  /* Box-Muller method for Gaussian random numbers */

  nn = n;
  if (nn <= 0)
    return;
  if (klotz1_1.first == 0) {
    _normal00();
    klotz1_1.first = 1;
  }
  ptr = 0;

L1:
  left = buffsz - klotz1_1.xptr;
  if (nn < left) {
    kptr = klotz1_1.xptr;
#   pragma omp simd
    for (i = 0; i < nn; ++i) {
      x[i + ptr] = klotz1_1.xbuff[kptr + i];
    }
    klotz1_1.xptr += nn;
    return;
  } else {
    kptr = klotz1_1.xptr;
#   pragma omp simd
    for (i = 0; i < left; ++i) {
      x[i + ptr] = klotz1_1.xbuff[kptr + i];
    }
    klotz1_1.xptr = 0;
    ptr += left;
    nn -= left;
    _normal00();
    goto L1;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Poisson distribution random number generator.
 *
 * q(mu,p) = exp(-mu) mu**p/p!
 *
 * \param[in]   n   number of values to compute
 * \param[in]   mu  Poisson distribution parameter
 * \param[out]  p   pseudo-random numbers following Poisson distribution
 */
/*----------------------------------------------------------------------------*/

void
cs_random_poisson(cs_lnum_t  n,
                  cs_real_t  mu,
                  int        p[])
{
  /* Local variables */
  int left, indx[1024], i, k;
  double q[1024], u[1024];
  int nsegs, p0;
  double q0;
  int ii, jj;
  int nl0;
  double pmu;

  if (n <= 0)
    return;

  pmu = exp(-mu);
  p0 = 0;

  nsegs = (n - 1) / 1024;
  left = n - (nsegs << 10);
  ++nsegs;
  nl0 = left;

  for (k = 0; k < nsegs; ++k) {

    for (i = 0; i < left; ++i) {
      indx[i] = i;
      p[p0 + i] = 0;
      q[i] = 1.;
    }

    /* Begin iterative loop on segment of p's */

  L1:

    /* Get the needed uniforms */

    cs_random_uniform(left, u);

    jj = 0;

    for (i = 0; i < left; ++i) {
      ii = indx[i];
      q0 = q[ii] * u[i];
      q[ii] = q0;
      if (q0 > pmu) {
        indx[jj++] = ii;
        ++p[p0 + ii];
      }
    }

    /* any left in this segment? */

    left = jj;
    if (left > 0) {
      goto L1;
    }

    p0  += nl0;
    nl0  = 1024;
    left = 1024;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Save static variables used by random number generator.
 *
 * \param[out]  save_block  saved state values
 */
/*----------------------------------------------------------------------------*/

void
cs_random_save(cs_real_t  save_block[1634])
{
  int i, k;

  /* The entire contents of blocks klotz0 and klotz1 must be saved. */

  if (klotz1_1.first == 0) {
    _normal00();
    klotz1_1.first = 1;
  }

  _random_uniform_save(save_block);

  save_block[608] = (double) klotz1_1.first;
  save_block[609] = (double) klotz1_1.xptr;
  k = 610;
# pragma omp simd
  for (i = 0; i < 1024; ++i) {
    save_block[i + k] = klotz1_1.xbuff[i];
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Restore static variables used by random number generator.
 *
 * \param[out]  save_block  saved state values
 */
/*----------------------------------------------------------------------------*/

void
cs_random_restore(cs_real_t  save_block[1634])
{
  int i, k;

  /* The entire contents of klotz0 and klotz1 must be restored. */

  _random_uniform_restore(save_block);
  klotz1_1.first = (int) save_block[608];
  if (klotz1_1.first == 0)
    bft_error(__FILE__, __LINE__, 0,
              "In %s, restore of uninitialized block.", __func__);
  klotz1_1.xptr = (int) save_block[609];
  k = 610;
# pragma omp simd
  for (i = 0; i < 1024; ++i) {
    klotz1_1.xbuff[i] = save_block[i + k];
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
