/*============================================================================
 * Portability and fallback layer for BLAS functions
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2013 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_parall.h"

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

double CS_PROCF(csdot, CSDOT)(const cs_int_t   *n,
                              const cs_real_t  *x,
                              const cs_real_t  *y)
{
  return cs_dot(*n, x, y);
}

/* Return the global residual of 2 extensiv vectors: x.y */

double CS_PROCF(csres, CSRES)(const cs_int_t   *n,
                              const cs_real_t  *vol,
                              const cs_real_t  *x,
                              const cs_real_t  *y)
{
  return cs_gres(*n, vol, x, y);
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

void cs_axpy(cs_lnum_t      n,
             double         a,
             const cs_real_t  *x,
             cs_real_t        *restrict y)
{
  cs_lnum_t  i;

  if (n < 1)
    return;

# pragma omp parallel for
  for (i = 0; i < n; i++)
    y[i] += (a * x[i]);
}

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
cs_dot(cs_lnum_t         n,
       const cs_real_t  *x,
       const cs_real_t  *y)
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
 * Return the global resildual of 2 extensive vectors:
 *  1/sum(vol) . sum(X.Y/vol)
 *
 * For better precision, a superblock algorithm is used.
 *
 * parameters:
 *   n   <-- size of arrays x and y
 *   vol <-- array of floating-point values
 *   x   <-- array of floating-point values
 *   y   <-- array of floating-point values
 *
 * returns:
 *   dot product
 *----------------------------------------------------------------------------*/

double
cs_gres(cs_lnum_t         n,
       const cs_real_t  *vol,
       const cs_real_t  *x,
       const cs_real_t  *y)
{
  const cs_lnum_t block_size = 60;

  cs_lnum_t sid, bid, i;
  cs_lnum_t start_id, end_id;
  double sdot, cdot;
  double svtot, cvtot;

  cs_lnum_t n_blocks = n / block_size;
  cs_lnum_t n_sblocks = sqrt(n_blocks);
  cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

  double dot = 0.;
  double vtot = 0.;

 /*
  * The algorithm used is l3superblock60, based on the article:
  * "Reducing Floating Point Error in Dot Product Using the Superblock Family
  * of Algorithms" by Anthony M. Castaldo, R. Clint Whaley, and Anthony
  * T. Chronopoulos, SIAM J. SCI. COMPUT., Vol. 31, No. 2, pp. 1156–1174
  * 2008 Society for Industrial and Applied Mathematics
  */

# pragma omp parallel for reduction(+:dot, vtot) private(bid, start_id, end_id, i, \
                                                   cdot, sdot, cvtot, svtot) if (n > THR_MIN)
  for (sid = 0; sid < n_sblocks; sid++) {

    sdot = 0.;
    svtot = 0.;

    for (bid = 0; bid < blocks_in_sblocks; bid++) {
      start_id = block_size * (blocks_in_sblocks*sid + bid);
      end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
      cdot = 0.;
      cvtot = 0.;
      for (i = start_id; i < end_id; i++) {
        cdot += x[i]*y[i]/vol[i];
        cvtot += vol[i];
      }
      sdot += cdot;
      svtot += cvtot;
    }

    dot += sdot;

  }

  cdot = 0.;
  cvtot = 0.;
  start_id = block_size * n_sblocks*blocks_in_sblocks;
  end_id = n;
  for (i = start_id; i < end_id; i++) {
    cdot += x[i]*y[i]/vol[i];
    cvtot += vol[i];
  }
  dot += cdot;
  vtot += cvtot;

  cs_parall_sum(1, CS_DOUBLE, &dot);
  cs_parall_sum(1, CS_DOUBLE, &vtot);

  dot /= vtot;

  return dot;
}
/*----------------------------------------------------------------------------
 * Return dot products of a vector with itself: x.x
 *
 * For better precision, a superblock algorithm is used.
 *
 * parameters:
 *   n  <-- size of arrays x and y
 *   x  <-- array of floating-point values
 *
 * returns:
 *   dot product
 *----------------------------------------------------------------------------*/

double
cs_dot_xx(cs_lnum_t         n,
          const cs_real_t  *x)
{
  const cs_lnum_t block_size = 60;

  cs_lnum_t sid, bid, i;
  cs_lnum_t start_id, end_id;
  double sdot_xx, cdot_xx;

  cs_lnum_t n_blocks = n / block_size;
  cs_lnum_t n_sblocks = sqrt(n_blocks);
  cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

  double dot_xx = 0.0;

# pragma omp parallel for private(bid, start_id, end_id, i, \
                                  cdot_xx, sdot_xx) \
                          reduction(+:dot_xx) if (n > THR_MIN)
  for (sid = 0; sid < n_sblocks; sid++) {

    sdot_xx = 0.0;

    for (bid = 0; bid < blocks_in_sblocks; bid++) {
      start_id = block_size * (blocks_in_sblocks*sid + bid);
      end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
      cdot_xx = 0.0;
      for (i = start_id; i < end_id; i++)
        cdot_xx += x[i]*x[i];
      sdot_xx += cdot_xx;
    }

    dot_xx += sdot_xx;

  }

  cdot_xx = 0.0;
  start_id = block_size * n_sblocks*blocks_in_sblocks;
  end_id = n;
  for (i = start_id; i < end_id; i++)
    cdot_xx += x[i]*x[i];
  dot_xx += cdot_xx;

  return dot_xx;
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
cs_dot_xx_xy(cs_lnum_t                    n,
             const cs_real_t  *restrict   x,
             const cs_real_t  *restrict   y,
             double                      *xx,
             double                      *xy)
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
cs_dot_xy_yz(cs_lnum_t                    n,
             const cs_real_t  *restrict   x,
             const cs_real_t  *restrict   y,
             const cs_real_t  *restrict   z,
             double                      *xy,
             double                      *yz)
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
cs_dot_xx_xy_yz(cs_lnum_t                    n,
                const cs_real_t  *restrict   x,
                const cs_real_t  *restrict   y,
                const cs_real_t  *restrict   z,
                double                      *xx,
                double                      *xy,
                double                      *yz)
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

/*----------------------------------------------------------------------------
 * Return the global dot product of 2 vectors: x.y
 *
 * In parallel mode, the local results are summed on the default
 * global communicator.
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
cs_gdot(cs_lnum_t         n,
        const cs_real_t  *x,
        const cs_real_t  *y)
{
  double retval = cs_dot(n, x, y);

  cs_parall_sum(1, CS_DOUBLE, &retval);

  return retval;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
