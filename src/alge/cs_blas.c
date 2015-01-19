/*============================================================================
 * Portability and fallback layer for BLAS functions
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Constant times a vector plus a vector: y <-- ax + y
 *
 * \param[in]       n  size of arrays x and y
 * \param[in]       a  multiplier for x
 * \param[in]       x  array of floating-point values
 * \param[in, out]  y  array of floating-point values
 */
/*----------------------------------------------------------------------------*/

void cs_axpy(cs_lnum_t         n,
             double            a,
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the dot product of 2 vectors: x.y
 *
 * \param[in]  n  size of arrays x and y
 * \param[in]  x  array of floating-point values
 * \param[in]  y  array of floating-point values
 *
 * \return  dot product
 */
/*----------------------------------------------------------------------------*/

double
cs_dot(cs_lnum_t         n,
       const cs_real_t  *x,
       const cs_real_t  *y)
{
  const cs_lnum_t block_size = 60;

  const cs_lnum_t n_blocks = n / block_size;
  const cs_lnum_t n_sblocks = sqrt(n_blocks);
  const cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

  double dot = 0.0;

 /*
  * The algorithm used is l3superblock60, based on the article:
  * "Reducing Floating Point Error in Dot Product Using the Superblock Family
  * of Algorithms" by Anthony M. Castaldo, R. Clint Whaley, and Anthony
  * T. Chronopoulos, SIAM J. SCI. COMPUT., Vol. 31, No. 2, pp. 1156–1174
  * 2008 Society for Industrial and Applied Mathematics
  */

# pragma omp parallel for reduction(+:dot) if (n > CS_THR_MIN)
  for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

    double sdot = 0.0;

    for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
      cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
      cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
      double cdot = 0.0;
      for (cs_lnum_t i = start_id; i < end_id; i++)
        cdot += x[i]*y[i];
      sdot += cdot;
    }

    dot += sdot;

  }

  double cdot = 0.0;
  cs_lnum_t start_id = block_size * n_sblocks*blocks_in_sblocks;
  cs_lnum_t end_id = n;
  for (cs_lnum_t i = start_id; i < end_id; i++)
    cdot += x[i]*y[i];
  dot += cdot;

  return dot;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the global residual of 2 extensive vectors:
 *        1/sum(vol) . sum(X.Y/vol)
 *
 * In parallel mode, the local results are summed on the default
 * global communicator.
 *
 * \param[in]  n    size of arrays x and y
 * \param[in]  vol  array of floating-point values
 * \param[in]  x    array of floating-point values
 * \param[in]  y    array of floating-point values
 *
 * \return  global residual
 */
/*----------------------------------------------------------------------------*/

double
cs_gres(cs_lnum_t         n,
        const cs_real_t  *vol,
        const cs_real_t  *x,
        const cs_real_t  *y)
{
  const cs_lnum_t block_size = 60;

  const cs_lnum_t n_blocks = n / block_size;
  const cs_lnum_t n_sblocks = sqrt(n_blocks);
  const cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

  double dot = 0.;
  double vtot = 0.;

 /*
  * The algorithm used is l3superblock60, based on the article:
  * "Reducing Floating Point Error in Dot Product Using the Superblock Family
  * of Algorithms" by Anthony M. Castaldo, R. Clint Whaley, and Anthony
  * T. Chronopoulos, SIAM J. SCI. COMPUT., Vol. 31, No. 2, pp. 1156–1174
  * 2008 Society for Industrial and Applied Mathematics
  */

# pragma omp parallel for reduction(+:dot, vtot) if (n > CS_THR_MIN)
  for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

    double sdot = 0.;
    double svtot = 0.;

    for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
      cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
      cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
      double cdot = 0.;
      double cvtot = 0.;
      for (cs_lnum_t i = start_id; i < end_id; i++) {
        cdot += x[i]*y[i]/vol[i];
        cvtot += vol[i];
      }
      sdot += cdot;
      svtot += cvtot;
    }

    dot += sdot;
    vtot += svtot;

  }

  double cdot = 0.;
  double cvtot = 0.;
  cs_lnum_t start_id = block_size * n_sblocks*blocks_in_sblocks;
  cs_lnum_t end_id = n;
  for (cs_lnum_t i = start_id; i < end_id; i++) {
    cdot += x[i]*y[i]/vol[i];
    cvtot += vol[i];
  }
  dot += cdot;
  vtot += cvtot;

  double atot[2] = {dot, vtot};
  cs_parall_sum(2, CS_DOUBLE, atot);

  dot = atot[0] / atot[1];

  return dot;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return dot products of a vector with itself: x.x
 *
 * For better precision, a superblock algorithm is used.
 *
 * \param[in]  n  size of array x
 * \param[in]  x  array of floating-point values
 *
 * \return  dot product
 */
/*----------------------------------------------------------------------------*/

double
cs_dot_xx(cs_lnum_t         n,
          const cs_real_t  *x)
{
  const cs_lnum_t block_size = 60;

  const cs_lnum_t n_blocks = n / block_size;
  const cs_lnum_t n_sblocks = sqrt(n_blocks);
  const cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

  double dot_xx = 0.0;

# pragma omp parallel for reduction(+:dot_xx) if (n > CS_THR_MIN)
  for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

    double sdot_xx = 0.0;

    for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
      cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
      cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
      double cdot_xx = 0.0;
      for (cs_lnum_t i = start_id; i < end_id; i++)
        cdot_xx += x[i]*x[i];
      sdot_xx += cdot_xx;
    }

    dot_xx += sdot_xx;

  }

  double cdot_xx = 0.0;
  cs_lnum_t start_id = block_size * n_sblocks*blocks_in_sblocks;
  cs_lnum_t end_id = n;
  for (cs_lnum_t i = start_id; i < end_id; i++)
    cdot_xx += x[i]*x[i];
  dot_xx += cdot_xx;

  return dot_xx;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return 2 dot products of 2 vectors: x.x, and x.y
 *
 * The products could be computed separately, but computing them
 * simultaneously adds more optimization opportunities and possibly better
 * cache behavior.
 *
 * For better precision, a superblock algorithm is used.
 *
 * \param[in]   n   size of arrays x and y
 * \param[in]   x   array of floating-point values
 * \param[in]   y   array of floating-point values
 * \param[out]  xx  x.x dot product
 * \param[out]  xy  x.y dot product
 */
/*----------------------------------------------------------------------------*/

void
cs_dot_xx_xy(cs_lnum_t                    n,
             const cs_real_t  *restrict   x,
             const cs_real_t  *restrict   y,
             double                      *xx,
             double                      *xy)
{
  const cs_lnum_t block_size = 60;

  const cs_lnum_t n_blocks = n / block_size;
  const cs_lnum_t n_sblocks = sqrt(n_blocks);
  const cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

  double dot_xx = 0.0;
  double dot_xy = 0.0;

#if defined(__xlc__)
#pragma disjoint(*x, *y, *xx, *xy)
#endif

# pragma omp parallel for reduction(+:dot_xx, dot_xy) if (n > CS_THR_MIN)
  for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

    double sdot_xx = 0.0;
    double sdot_xy = 0.0;

    for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
      cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
      cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
      double cdot_xx = 0.0;
      double cdot_xy = 0.0;
      for (cs_lnum_t i = start_id; i < end_id; i++) {
        cdot_xx += x[i]*x[i];
        cdot_xy += x[i]*y[i];
      }
      sdot_xx += cdot_xx;
      sdot_xy += cdot_xy;
    }

    dot_xx += sdot_xx;
    dot_xy += sdot_xy;

  }

  double cdot_xx = 0.0;
  double cdot_xy = 0.0;
  cs_lnum_t start_id = block_size * n_sblocks*blocks_in_sblocks;
  cs_lnum_t end_id = n;
  for (cs_lnum_t i = start_id; i < end_id; i++) {
    cdot_xx += x[i]*x[i];
    cdot_xy += x[i]*y[i];
  }
  dot_xx += cdot_xx;
  dot_xy += cdot_xy;

  *xx = dot_xx;
  *xy = dot_xy;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return 2 dot products of 3 vectors: x.x, and y.z
 *
 * The products could be computed separately, but computing them
 * simultaneously adds more optimization opportunities and possibly better
 * cache behavior.
 *
 * For better precision, a superblock algorithm is used.
 *
 * \param[in]   n   size of arrays x and y, and z
 * \param[in]   x   array of floating-point values
 * \param[in]   y   array of floating-point values
 * \param[in]   z   array of floating-point values
 * \param[out]  xy  x.y dot product
 * \param[out]  yz  y.z dot product
 */
/*----------------------------------------------------------------------------*/

void
cs_dot_xy_yz(cs_lnum_t                    n,
             const cs_real_t  *restrict   x,
             const cs_real_t  *restrict   y,
             const cs_real_t  *restrict   z,
             double                      *xy,
             double                      *yz)
{
  const cs_lnum_t block_size = 60;

  const cs_lnum_t n_blocks = n / block_size;
  const cs_lnum_t n_sblocks = sqrt(n_blocks);
  const cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

  double dot_xy = 0.0;
  double dot_yz = 0.0;

#if defined(__xlc__)
#pragma disjoint(*x, *y, *xy, *yz)
#endif

# pragma omp parallel for reduction(+:dot_xy, dot_yz) if (n > CS_THR_MIN)
  for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

    double sdot_xy = 0.0;
    double sdot_yz = 0.0;

    for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
      cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
      cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
      double cdot_xy = 0.0;
      double cdot_yz = 0.0;
      for (cs_lnum_t i = start_id; i < end_id; i++) {
        cdot_xy += x[i]*y[i];
        cdot_yz += y[i]*z[i];
      }
      sdot_xy += cdot_xy;
      sdot_yz += cdot_yz;
    }

    dot_xy += sdot_xy;
    dot_yz += sdot_yz;

  }

  double cdot_xy = 0.0;
  double cdot_yz = 0.0;
  cs_lnum_t start_id = block_size * n_sblocks*blocks_in_sblocks;
  cs_lnum_t end_id = n;
  for (cs_lnum_t i = start_id; i < end_id; i++) {
    cdot_xy += x[i]*y[i];
    cdot_yz += y[i]*z[i];
  }
  dot_xy += cdot_xy;
  dot_yz += cdot_yz;

  *xy = dot_xy;
  *yz = dot_yz;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return 3 dot products of 3 vectors: x.x, x.y, and y.z
 *
 * The products could be computed separately, but computing them
 * simultaneously adds more optimization opportunities and possibly better
 * cache behavior.
 *
 * For better precision, a superblock algorithm is used.
 *
 * \param[in]   n   size of arrays x and y, and z
 * \param[in]   x   array of floating-point values
 * \param[in]   y   array of floating-point values
 * \param[in]   z   array of floating-point values
 * \param[out]  xx  x.x dot product
 * \param[out]  xy  x.y dot product
 * \param[out]  yz  y.z dot product
 */
/*----------------------------------------------------------------------------*/

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

  const cs_lnum_t n_blocks = n / block_size;
  const cs_lnum_t n_sblocks = sqrt(n_blocks);
  const cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

  double dot_xx = 0.0;
  double dot_xy = 0.0;
  double dot_yz = 0.0;

#if defined(__xlc__)
#pragma disjoint(*x, *y, *xy, *yz)
#endif

# pragma omp parallel for reduction(+:dot_xx, dot_xy, dot_yz) if (n > CS_THR_MIN)
  for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

    double sdot_xx = 0.0;
    double sdot_xy = 0.0;
    double sdot_yz = 0.0;

    for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
      cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
      cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
      double cdot_xx = 0.0;
      double cdot_xy = 0.0;
      double cdot_yz = 0.0;
      for (cs_lnum_t i = start_id; i < end_id; i++) {
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

  double cdot_xx = 0.0;
  double cdot_xy = 0.0;
  double cdot_yz = 0.0;
  cs_lnum_t start_id = block_size * n_sblocks*blocks_in_sblocks;
  cs_lnum_t end_id = n;
  for (cs_lnum_t i = start_id; i < end_id; i++) {
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
/*!
 * \brief Return 5 dot products of 3 vectors: x.x, y.y, x.y, x.z, and y.z
 *
 * The products could be computed separately, but computing them
 * simultaneously adds more optimization opportunities and possibly better
 * cache behavior.
 *
 * For better precision, a superblock algorithm is used.
 *
 * \param[in]   n   size of arrays x and y, and z
 * \param[in]   x   array of floating-point values
 * \param[in]   y   array of floating-point values
 * \param[in]   z   array of floating-point values
 * \param[out]  xx  x.x dot product
 * \param[out]  yy  y.y dot product
 * \param[out]  xy  x.y dot product
 * \param[out]  xz  x.z dot product
 * \param[out]  yz  y.z dot product
 */
/*----------------------------------------------------------------------------*/

void
cs_dot_xx_yy_xy_xz_yz(cs_lnum_t                    n,
                      const cs_real_t  *restrict   x,
                      const cs_real_t  *restrict   y,
                      const cs_real_t  *restrict   z,
                      double                      *xx,
                      double                      *yy,
                      double                      *xy,
                      double                      *xz,
                      double                      *yz)
{
  const cs_lnum_t block_size = 60;

  const cs_lnum_t n_blocks = n / block_size;
  const cs_lnum_t n_sblocks = sqrt(n_blocks);
  const cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

  double dot_xx = 0.0;
  double dot_yy = 0.0;
  double dot_xy = 0.0;
  double dot_xz = 0.0;
  double dot_yz = 0.0;

# pragma omp parallel for reduction(+:dot_xx, dot_yy, dot_xy, dot_xz, dot_yz) \
                          if (n > CS_THR_MIN)
  for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

    double sdot_xx = 0.0;
    double sdot_yy = 0.0;
    double sdot_xy = 0.0;
    double sdot_xz = 0.0;
    double sdot_yz = 0.0;

    for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
      cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
      cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
      double cdot_xx = 0.0;
      double cdot_yy = 0.0;
      double cdot_xy = 0.0;
      double cdot_xz = 0.0;
      double cdot_yz = 0.0;
      for (cs_lnum_t i = start_id; i < end_id; i++) {
        cdot_xx += x[i]*x[i];
        cdot_yy += y[i]*y[i];
        cdot_xy += x[i]*y[i];
        cdot_xz += x[i]*z[i];
        cdot_yz += y[i]*z[i];
      }
      sdot_xx += cdot_xx;
      sdot_yy += cdot_yy;
      sdot_xy += cdot_xy;
      sdot_xz += cdot_xz;
      sdot_yz += cdot_yz;
    }

    dot_xx += sdot_xx;
    dot_yy += sdot_yy;
    dot_xy += sdot_xy;
    dot_xz += sdot_xz;
    dot_yz += sdot_yz;

  }

  double cdot_xx = 0.0;
  double cdot_yy = 0.0;
  double cdot_xy = 0.0;
  double cdot_xz = 0.0;
  double cdot_yz = 0.0;
  cs_lnum_t start_id = block_size * n_sblocks*blocks_in_sblocks;
  cs_lnum_t end_id = n;
  for (cs_lnum_t i = start_id; i < end_id; i++) {
    cdot_xx += x[i]*x[i];
    cdot_yy += y[i]*y[i];
    cdot_xy += x[i]*y[i];
    cdot_xz += x[i]*z[i];
    cdot_yz += y[i]*z[i];
  }
  dot_xx += cdot_xx;
  dot_yy += cdot_yy;
  dot_xy += cdot_xy;
  dot_xz += cdot_xz;
  dot_yz += cdot_yz;

  *xx = dot_xx;
  *yy = dot_yy;
  *xy = dot_xy;
  *xz = dot_xz;
  *yz = dot_yz;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the global dot product of 2 vectors: x.y
 *
 * In parallel mode, the local results are summed on the default
 * global communicator.
 *
 * For better precision, a superblock algorithm is used.
 *
 * \param[in]  n  size of arrays x and y
 * \param[in]  x  array of floating-point values
 * \param[in]  y  array of floating-point values
 *
 * \return  dot product
 */
/*----------------------------------------------------------------------------*/

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
