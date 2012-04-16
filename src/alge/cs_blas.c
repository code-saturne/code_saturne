/*============================================================================
 * Portability and fallback layer for BLAS functions
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2012 EDF S.A.

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
 * External library headers
 *----------------------------------------------------------------------------*/

#include <math.h>

#if defined(HAVE_CBLAS)
#include <cblas.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_blas.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* Minimum size for OpenMP loops (needs benchmarking to adjust) */
#define THR_MIN 128

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 *  Public function definitions for Fortran API
 *============================================================================*/

/* Return the dot product of 2 vectors: x.y */

double CS_PROCF(csdot, CSDOT)(const cs_int_t  *n,
                              const double    *x,
                              const double    *y)
{
  return cs_dot(*n, x, y);
}

/*============================================================================
 *  Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Constant times a vector plus a vector: y <-- ax + y
 *
 * parameters:
 *   n <-- size of arrays x and y
 *   a <-- multiplier for x
 *   x <-- array of floating-point values
 *   y <-- array of floating-point values
 *----------------------------------------------------------------------------*/

#if    !defined(HAVE_ACML) && !defined(HAVE_ESSL) \
    && !defined(HAVE_ATLAS) && !defined(HAVE_CBLAS) && !defined(HAVE_MKL)

void cs_axpy(cs_lnum_t      n,
             double         a,
             const double  *x,
             double        *restrict y)
{
  cs_lnum_t  i;

  if (n < 1)
    return;

# pragma omp parallel for
  for (i = 0; i < n; i++)
    y[i] += (a * x[i]);
}

#endif /* BLAS defined */

/*----------------------------------------------------------------------------
 * Return the dot product of 2 vectors: x.y
 *
 * For better precision, a superblock algorithm is used.
 *
 * parameters:
 *   n <-- size of arrays x and y
 *   x <-- array of floating-point values
 *   y <-- array of floating-point values
 *
 * returns:
 *   dot product
 *----------------------------------------------------------------------------*/

double
cs_dot(cs_lnum_t      n,
       const double  *x,
       const double  *y)
{
  const cs_lnum_t block_size = 60;

  cs_lnum_t sid, bid, i;
  cs_lnum_t start_id, end_id;
  double sdot, cdot;

  cs_lnum_t n_blocks = n / block_size;
  cs_lnum_t n_sblocks = sqrt(n_blocks);
  cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

  double dot = 0.0;

 /*
  * The algorithm used is l3superblock60, based on the article:
  * "Reducing Floating Point Error in Dot Product Using the Superblock Family
  * of Algorithms" by Anthony M. Castaldo, R. Clint Whaley, and Anthony
  * T. Chronopoulos, SIAM J. SCI. COMPUT., Vol. 31, No. 2, pp. 1156–1174
  * 2008 Society for Industrial and Applied Mathematics
  */

# pragma omp parallel for reduction(+:dot) private(bid, start_id, end_id, i, \
                                                   cdot, sdot) if (n > THR_MIN)
  for (sid = 0; sid < n_sblocks; sid++) {

    sdot = 0.0;

    for (bid = 0; bid < blocks_in_sblocks; bid++) {
      start_id = block_size * (blocks_in_sblocks*sid + bid);
      end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
      cdot = 0.0;
      for (i = start_id; i < end_id; i++)
        cdot += x[i]*y[i];
      sdot += cdot;
    }

    dot += sdot;

  }

  cdot = 0.0;
  start_id = block_size * n_sblocks*blocks_in_sblocks;
  end_id = n;
  for (i = start_id; i < end_id; i++)
    cdot += x[i]*y[i];
  dot += cdot;

  return dot;
}

/*----------------------------------------------------------------------------
 * Return 2 dot products of 2 vectors: x.x, and x.y
 *
 * The products could be computed separately, but computing them
 * simultaneously adds more optimization opportunities and possibly better
 * cache behavior.
 *
 * For better precision, a superblock algorithm is used.
 *
 * parameters:
 *   n  <-- size of arrays x and y
 *   x  <-- array of floating-point values
 *   y  <-- array of floating-point values
 *   xx --> x.x dot product
 *   xy --> x.y dot product
 *----------------------------------------------------------------------------*/

void
cs_dot_xx_xy(cs_lnum_t               n,
             const double  *restrict x,
             const double  *restrict y,
             double                 *xx,
             double                 *xy)
{
  const cs_lnum_t block_size = 60;

  cs_lnum_t sid, bid, i;
  cs_lnum_t start_id, end_id;
  double sdot_xx, sdot_xy, cdot_xx, cdot_xy;

  cs_lnum_t n_blocks = n / block_size;
  cs_lnum_t n_sblocks = sqrt(n_blocks);
  cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

  double dot_xx = 0.0;
  double dot_xy = 0.0;

#if defined(__xlc__)
#pragma disjoint(*x, *y, *xx, *xy)
#endif

# pragma omp parallel for private(bid, start_id, end_id, i, \
                                  cdot_xx, cdot_xy, sdot_xx, sdot_xy) \
                          reduction(+:dot_xx, dot_xy) if (n > THR_MIN)
  for (sid = 0; sid < n_sblocks; sid++) {

    sdot_xx = 0.0;
    sdot_xy = 0.0;

    for (bid = 0; bid < blocks_in_sblocks; bid++) {
      start_id = block_size * (blocks_in_sblocks*sid + bid);
      end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
      cdot_xx = 0.0;
      cdot_xy = 0.0;
      for (i = start_id; i < end_id; i++) {
        cdot_xx += x[i]*x[i];
        cdot_xy += x[i]*y[i];
      }
      sdot_xx += cdot_xx;
      sdot_xy += cdot_xy;
    }

    dot_xx += sdot_xx;
    dot_xy += sdot_xy;

  }

  cdot_xx = 0.0;
  cdot_xy = 0.0;
  start_id = block_size * n_sblocks*blocks_in_sblocks;
  end_id = n;
  for (i = start_id; i < end_id; i++) {
    cdot_xx += x[i]*x[i];
    cdot_xy += x[i]*y[i];
  }
  dot_xx += cdot_xx;
  dot_xy += cdot_xy;

  *xx = dot_xx;
  *xy = dot_xy;
}

/*----------------------------------------------------------------------------
 * Return 2 dot products of 3 vectors: x.y, and y.z
 *
 * The products could be computed separately, but computing them
 * simultaneously adds more optimization opportunities and possibly better
 * cache behavior.
 *
 * For better precision, a superblock algorithm is used.
 *
 * parameters:
 *   n  <-- size of arrays x and y
 *   x  <-- array of floating-point values
 *   y  <-- array of floating-point values
 *   z  <-- array of floating-point values
 *   xy --> x.y dot product
 *   yz --> y.z dot product
 *----------------------------------------------------------------------------*/

void
cs_dot_xy_yz(cs_lnum_t               n,
             const double  *restrict x,
             const double  *restrict y,
             const double  *restrict z,
             double                 *xy,
             double                 *yz)
{
  const cs_lnum_t block_size = 60;

  cs_lnum_t sid, bid, i;
  cs_lnum_t start_id, end_id;
  double sdot_xy, sdot_yz, cdot_xy, cdot_yz;

  cs_lnum_t n_blocks = n / block_size;
  cs_lnum_t n_sblocks = sqrt(n_blocks);
  cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

  double dot_xy = 0.0;
  double dot_yz = 0.0;

#if defined(__xlc__)
#pragma disjoint(*x, *y, *xy, *yz)
#endif

# pragma omp parallel for private(bid, start_id, end_id, i, \
                                  cdot_xy, cdot_yz, sdot_xy, sdot_yz) \
                                  reduction(+:dot_xy, dot_yz) if (n > THR_MIN)
  for (sid = 0; sid < n_sblocks; sid++) {

    sdot_xy = 0.0;
    sdot_yz = 0.0;

    for (bid = 0; bid < blocks_in_sblocks; bid++) {
      start_id = block_size * (blocks_in_sblocks*sid + bid);
      end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
      cdot_xy = 0.0;
      cdot_yz = 0.0;
      for (i = start_id; i < end_id; i++) {
        cdot_xy += x[i]*y[i];
        cdot_yz += y[i]*z[i];
      }
      sdot_xy += cdot_xy;
      sdot_yz += cdot_yz;
    }

    dot_xy += sdot_xy;
    dot_yz += sdot_yz;

  }

  cdot_xy = 0.0;
  cdot_yz = 0.0;
  start_id = block_size * n_sblocks*blocks_in_sblocks;
  end_id = n;
  for (i = start_id; i < end_id; i++) {
    cdot_xy += x[i]*y[i];
    cdot_yz += y[i]*z[i];
  }
  dot_xy += cdot_xy;
  dot_yz += cdot_yz;

  *xy = dot_xy;
  *yz = dot_yz;
}

/*----------------------------------------------------------------------------
 * Return 3 dot products of 3 vectors: x.y, x.y, and y.z
 *
 * The products could be computed separately, but computing them
 * simultaneously adds more optimization opportunities and possibly better
 * cache behavior.
 *
 * For better precision, a superblock algorithm is used.
 *
 * parameters:
 *   n  <-- size of arrays x and y
 *   x  <-- array of floating-point values
 *   y  <-- array of floating-point values
 *   z  <-- array of floating-point values
 *   xx --> x.y dot product
 *   xy --> x.y dot product
 *   yz --> y.z dot product
 *----------------------------------------------------------------------------*/

void
cs_dot_xx_xy_yz(cs_lnum_t               n,
                const double  *restrict x,
                const double  *restrict y,
                const double  *restrict z,
                double                 *xx,
                double                 *xy,
                double                 *yz)
{
  const cs_lnum_t block_size = 60;

  cs_lnum_t sid, bid, i;
  cs_lnum_t start_id, end_id;
  double sdot_xx, sdot_xy, sdot_yz, cdot_xx, cdot_xy, cdot_yz;

  cs_lnum_t n_blocks = n / block_size;
  cs_lnum_t n_sblocks = sqrt(n_blocks);
  cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

  double dot_xx = 0.0;
  double dot_xy = 0.0;
  double dot_yz = 0.0;

#if defined(__xlc__)
#pragma disjoint(*x, *y, *xy, *yz)
#endif

# pragma omp parallel for private(bid, start_id, end_id, i, cdot_xx, cdot_xy, \
                                  cdot_yz, sdot_xx, sdot_xy, sdot_yz) \
                          reduction(+:dot_xx, dot_xy, dot_yz) if (n > THR_MIN)
  for (sid = 0; sid < n_sblocks; sid++) {

    sdot_xx = 0.0;
    sdot_xy = 0.0;
    sdot_yz = 0.0;

    for (bid = 0; bid < blocks_in_sblocks; bid++) {
      start_id = block_size * (blocks_in_sblocks*sid + bid);
      end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
      cdot_xx = 0.0;
      cdot_xy = 0.0;
      cdot_yz = 0.0;
      for (i = start_id; i < end_id; i++) {
        cdot_xx += x[i]*x[i];
        cdot_xy += x[i]*y[i];
        cdot_yz += y[i]*z[i];
      }
      sdot_xx += cdot_xx;
      sdot_xy += cdot_xy;
      sdot_yz += cdot_yz;
    }

    dot_xx += sdot_xx;
    dot_xy += sdot_xy;
    dot_yz += sdot_yz;

  }

  cdot_xx = 0.0;
  cdot_xy = 0.0;
  cdot_yz = 0.0;
  start_id = block_size * n_sblocks*blocks_in_sblocks;
  end_id = n;
  for (i = start_id; i < end_id; i++) {
    cdot_xx += x[i]*x[i];
    cdot_xy += x[i]*y[i];
    cdot_yz += y[i]*z[i];
  }
  dot_xx += cdot_xx;
  dot_xy += cdot_xy;
  dot_yz += cdot_yz;

  *xx = dot_xx;
  *xy = dot_xy;
  *yz = dot_yz;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
