/*============================================================================
 * BLAS (Basic Linear Algebra Subroutine) functions
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
 * External library headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>

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

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_blas.c
        BLAS (Basic Linear Algebra Subroutines) type functions

  \enum cs_blas_reduce_t

  \brief Reduction algorithm type

  \var CS_BLAS_REDUCE_SUPERBLOCK
       Reduction based on l3superblock60 algorithm, described in
       \cite Castaldo:2008

  \var CS_BLAS_REDUCE_KAHAN
       Reduction based on Kahan's compensated summation, described in
       \cite Kahan:1965
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* Block size for superblock algorithm */

#define CS_SBLOCK_BLOCK_SIZE 60

/* Cache line multiple, in cs_real_t units */

#define CS_CL  (CS_CL_SIZE/8)

/* Vector size, in cs_real_t units (large enough for AVX512) */

#define CS_VS  8

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Function pointer types
 *----------------------------------------------------------------------------*/

typedef double
(cs_dot_t) (cs_lnum_t         n,
            const cs_real_t  *x,
            const cs_real_t  *y);

typedef double
(cs_dot_xx_t) (cs_lnum_t         n,
               const cs_real_t  *x);

typedef void
(cs_dot_xx_xy_t) (cs_lnum_t         n,
                  const cs_real_t  *x,
                  const cs_real_t  *y,
                  double           *xx,
                  double           *xy);

typedef void
(cs_dot_xy_yz_t) (cs_lnum_t         n,
                  const cs_real_t  *x,
                  const cs_real_t  *y,
                  const cs_real_t  *z,
                  double           *xy,
                  double           *yz);

typedef void
(cs_dot_xx_xy_yz_t) (cs_lnum_t         n,
                     const cs_real_t  *x,
                     const cs_real_t  *y,
                     const cs_real_t  *z,
                     double           *xx,
                     double           *xy,
                     double           *yz);

typedef void
(cs_dot_xx_yy_xy_xz_yz_t)(cs_lnum_t                    n,
                          const cs_real_t  *restrict   x,
                          const cs_real_t  *restrict   y,
                          const cs_real_t  *restrict   z,
                          double                      *xx,
                          double                      *yy,
                          double                      *xy,
                          double                      *xz,
                          double                      *yz);

typedef double
(cs_gres_t)(cs_lnum_t         n,
            const cs_real_t  *vol,
            const cs_real_t  *x,
            const cs_real_t  *y);

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute array index bounds for a local thread.
 *
 * When called inside an OpenMP parallel section, this will return the
 * start an past-the-end indexes for the array range assigned to that thread.
 * In other cases, the start index is 1, and the past-the-end index is n;
 *
 * \param[in]   n     size of array
 * \param[out]  s_id  start index for the current thread
 * \param[out]  e_id  past-the-end index for the current thread
 */
/*----------------------------------------------------------------------------*/

static void
_thread_range(cs_lnum_t   n,
              cs_lnum_t  *s_id,
              cs_lnum_t  *e_id)
{
#if defined(HAVE_OPENMP)
  int t_id = omp_get_thread_num();
  int n_t = omp_get_num_threads();
  cs_lnum_t t_n = (n + n_t - 1) / n_t;
  *s_id =  t_id    * t_n;
  *e_id = (t_id+1) * t_n;
  *s_id = cs_align(*s_id, CS_CL);
  *e_id = cs_align(*e_id, CS_CL);
  if (*e_id > n) *e_id = n;
#else
  *s_id = 0;
  *e_id = n;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute blocks sizes for superblock algorithm.
 *
 * \param[in]   n                 size of array
 * \param[in]   block_size        block size
 * \param[out]  n_sblocks         number of superblocks
 * \param[out]  blocks_in_sblocks number of blocks per superblock
 */
/*----------------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the dot product of 2 vectors: x.y
 *        using a superblock algorithm.
 *
 * \param[in]  n  size of arrays x and y
 * \param[in]  x  array of floating-point values
 * \param[in]  y  array of floating-point values
 *
 * \return  dot product
 */
/*----------------------------------------------------------------------------*/

static double
_cs_dot_superblock(cs_lnum_t         n,
                   const cs_real_t  *x,
                   const cs_real_t  *y)
{
  double dot = 0.0;

# pragma omp parallel reduction(+:dot) if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(n, &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_real_t *_x = x + s_id;
    const cs_real_t *_y = y + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double sdot = 0.0;

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
        cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        if (end_id > _n)
          end_id = _n;
        double cdot = 0.0;
        for (cs_lnum_t i = start_id; i < end_id; i++)
        cdot += _x[i]*_y[i];
        sdot += cdot;
      }

      dot += sdot;

    }

  }

  return dot;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return dot products of a vector with itself: x.x
 *        using a superblock algorithm.
 *
 * \param[in]  n  size of array x
 * \param[in]  x  array of floating-point values
 *
 * \return  dot product
 */
/*----------------------------------------------------------------------------*/

static double
_cs_dot_xx_superblock(cs_lnum_t         n,
                      const cs_real_t  *x)
{
  double dot_xx = 0.0;

# pragma omp parallel reduction(+:dot_xx) if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(n, &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_real_t *_x = x + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double sdot_xx = 0.0;

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
        cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        if (end_id > _n)
          end_id = _n;
        double cdot_xx = 0.0;
        for (cs_lnum_t i = start_id; i < end_id; i++)
          cdot_xx += _x[i]*_x[i];
        sdot_xx += cdot_xx;
      }

      dot_xx += sdot_xx;

    }

  }

  return dot_xx;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return 2 dot products of 2 vectors: x.x, and x.y
 *        using a superblock algorithm.
 *
 * The products could be computed separately, but computing them
 * simultaneously adds more optimization opportunities and possibly better
 * cache behavior.
 *
 * \param[in]   n   size of arrays x and y
 * \param[in]   x   array of floating-point values
 * \param[in]   y   array of floating-point values
 * \param[out]  xx  x.x dot product
 * \param[out]  xy  x.y dot product
 */
/*----------------------------------------------------------------------------*/

static void
_cs_dot_xx_xy_superblock(cs_lnum_t                    n,
                         const cs_real_t  *restrict   x,
                         const cs_real_t  *restrict   y,
                         double                      *xx,
                         double                      *xy)
{
  double dot_xx = 0.0, dot_xy = 0.0;

# pragma omp parallel reduction(+:dot_xx, dot_xy) if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(n, &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_real_t *_x = x + s_id;
    const cs_real_t *_y = y + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double sdot_xx = 0.0;
      double sdot_xy = 0.0;

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
        cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        if (end_id > _n)
          end_id = _n;
        double cdot_xx = 0.0;
        double cdot_xy = 0.0;
        for (cs_lnum_t i = start_id; i < end_id; i++) {
          cdot_xx += _x[i]*_x[i];
          cdot_xy += _x[i]*_y[i];
        }
        sdot_xx += cdot_xx;
        sdot_xy += cdot_xy;
      }

      dot_xx += sdot_xx;
      dot_xy += sdot_xy;

    }

  }

  *xx = dot_xx;
  *xy = dot_xy;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return 2 dot products of 3 vectors: x.x, and y.z
 *        using a superblock algorithm.
 *
 * The products could be computed separately, but computing them
 * simultaneously adds more optimization opportunities and possibly better
 * cache behavior.
 *
 * \param[in]   n   size of arrays x and y, and z
 * \param[in]   x   array of floating-point values
 * \param[in]   y   array of floating-point values
 * \param[in]   z   array of floating-point values
 * \param[out]  xy  x.y dot product
 * \param[out]  yz  y.z dot product
 */
/*----------------------------------------------------------------------------*/

static void
_cs_dot_xy_yz_superblock(cs_lnum_t                    n,
                         const cs_real_t  *restrict   x,
                         const cs_real_t  *restrict   y,
                         const cs_real_t  *restrict   z,
                         double                      *xy,
                         double                      *yz)
{
  double dot_xy = 0.0, dot_yz = 0.0;

# pragma omp parallel reduction(+:dot_xy, dot_yz) if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(n, &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_real_t *_x = x + s_id;
    const cs_real_t *_y = y + s_id;
    const cs_real_t *_z = z + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double sdot_xy = 0.0;
      double sdot_yz = 0.0;

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
        cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        if (end_id > _n)
          end_id = _n;
        double cdot_xy = 0.0;
        double cdot_yz = 0.0;
        for (cs_lnum_t i = start_id; i < end_id; i++) {
          cdot_xy += _x[i]*_y[i];
          cdot_yz += _y[i]*_z[i];
        }
        sdot_xy += cdot_xy;
        sdot_yz += cdot_yz;
      }

      dot_xy += sdot_xy;
      dot_yz += sdot_yz;

    }
  }

  *xy = dot_xy;
  *yz = dot_yz;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return 3 dot products of 3 vectors: x.x, x.y, and y.z
 *        using a superblock algorithm.

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

static void
_cs_dot_xx_xy_yz_superblock(cs_lnum_t                    n,
                            const cs_real_t  *restrict   x,
                            const cs_real_t  *restrict   y,
                            const cs_real_t  *restrict   z,
                            double                      *xx,
                            double                      *xy,
                            double                      *yz)
{
  double dot_xx = 0.0, dot_xy = 0.0, dot_yz = 0.0;

# pragma omp parallel reduction(+:dot_xx, dot_xy, dot_yz) if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(n, &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_real_t *_x = x + s_id;
    const cs_real_t *_y = y + s_id;
    const cs_real_t *_z = z + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double sdot_xx = 0.0;
      double sdot_xy = 0.0;
      double sdot_yz = 0.0;

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
        cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        if (end_id > _n)
          end_id = _n;
        double cdot_xx = 0.0;
        double cdot_xy = 0.0;
        double cdot_yz = 0.0;
        for (cs_lnum_t i = start_id; i < end_id; i++) {
          cdot_xx += _x[i]*_x[i];
          cdot_xy += _x[i]*_y[i];
          cdot_yz += _y[i]*_z[i];
        }
        sdot_xx += cdot_xx;
        sdot_xy += cdot_xy;
        sdot_yz += cdot_yz;
      }

      dot_xx += sdot_xx;
      dot_xy += sdot_xy;
      dot_yz += sdot_yz;

    }
  }

  *xx = dot_xx;
  *xy = dot_xy;
  *yz = dot_yz;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return 5 dot products of 3 vectors: x.x, y.y, x.y, x.z, and y.z
 *        using a superblock algorithm.
 *
 * The products could be computed separately, but computing them
 * simultaneously adds more optimization opportunities and possibly better
 * cache behavior.
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

static void
_cs_dot_xx_yy_xy_xz_yz_superblock(cs_lnum_t                    n,
                                  const cs_real_t  *restrict   x,
                                  const cs_real_t  *restrict   y,
                                  const cs_real_t  *restrict   z,
                                  double                      *xx,
                                  double                      *yy,
                                  double                      *xy,
                                  double                      *xz,
                                  double                      *yz)
{
  double dot_xx = 0.0, dot_yy = 0.0;
  double dot_xy = 0.0, dot_xz = 0.0, dot_yz = 0.0;

# pragma omp parallel reduction(+:dot_xx, dot_yy, dot_xy, dot_xz, dot_yz) \
                      if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(n, &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_real_t *_x = x + s_id;
    const cs_real_t *_y = y + s_id;
    const cs_real_t *_z = z + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double sdot_xx = 0.0;
      double sdot_yy = 0.0;
      double sdot_xy = 0.0;
      double sdot_xz = 0.0;
      double sdot_yz = 0.0;

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
        cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        if (end_id > _n)
          end_id = _n;
        double cdot_xx = 0.0;
        double cdot_yy = 0.0;
        double cdot_xy = 0.0;
        double cdot_xz = 0.0;
        double cdot_yz = 0.0;
        for (cs_lnum_t i = start_id; i < end_id; i++) {
          cdot_xx += _x[i]*_x[i];
          cdot_yy += _y[i]*_y[i];
          cdot_xy += _x[i]*_y[i];
          cdot_xz += _x[i]*_z[i];
          cdot_yz += _y[i]*_z[i];
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
  }

  *xx = dot_xx;
  *yy = dot_yy;
  *xy = dot_xy;
  *xz = dot_xz;
  *yz = dot_yz;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the global residual of 2 extensive vectors:
 *        1/sum(vol) . sum(X.Y.vol)
 *        using a superblock algorithm.
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

static double
_cs_gres_superblock(cs_lnum_t         n,
                    const cs_real_t  *vol,
                    const cs_real_t  *x,
                    const cs_real_t  *y)
{
  double dot = 0.;
  double vtot = 0.;

# pragma omp parallel reduction(+:dot, vtot) if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(n, &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_real_t *_vol = vol + s_id;
    const cs_real_t *_x = x + s_id;
    const cs_real_t *_y = y + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double sdot = 0.;
      double svtot = 0.;

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
        cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        if (end_id > _n)
          end_id = _n;
        double cdot = 0.;
        double cvtot = 0.;
        for (cs_lnum_t i = start_id; i < end_id; i++) {
          cdot += _x[i] * _y[i] * _vol[i];
          cvtot += _vol[i];
        }
        sdot += cdot;
        svtot += cvtot;
      }

      dot += sdot;
      vtot += svtot;

    }

  }

  double atot[2] = {dot, vtot};
  cs_parall_sum(2, CS_DOUBLE, atot);

  dot = atot[0] / atot[1];

  return dot;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the dot product of 2 vectors: x.y
 *        using Kahan summation.
 *
 * \param[in]  n  size of arrays x and y
 * \param[in]  x  array of floating-point values
 * \param[in]  y  array of floating-point values
 *
 * \return  dot product
 */
/*----------------------------------------------------------------------------*/

static double
_cs_dot_kahan(cs_lnum_t         n,
              const cs_real_t  *x,
              const cs_real_t  *y)
{
  double dot = 0.0;

# pragma omp parallel reduction(+:dot) if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(n, &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_lnum_t _nv = (_n/CS_VS)*CS_VS;
    const cs_real_t *_x = x + s_id;
    const cs_real_t *_y = y + s_id;

    double d[CS_VS], c[CS_VS];

    for (cs_lnum_t i = 0; i < CS_VS; i ++) {
      d[i] = 0;
      c[i] = 0;
    }

    /* vector multiple */

    for (cs_lnum_t i = 0; i < _nv; i += CS_VS) {
      double z[CS_VS], t[CS_VS];
      for (cs_lnum_t j = 0; j < CS_VS; j++) {
        z[j] = (_x[i+j]* _y[i+j]) - c[j];
        t[j] = d[j] + z[j];
        c[j] = (t[j] - d[j]) - z[j];
        d[j] = t[j];
      }
    }

    /* vector remainder */

    for (cs_lnum_t i = _nv; i < _n; i++) {
      double z = (_x[i]* _y[i]) - c[0];
      double t = d[0] + z;
      c[0] = (t - d[0]) - z;
      d[0] = t;
    }

    /* Kahan sum on vector */

    {
      c[0] = 0;
      for (cs_lnum_t i = 0; i < CS_VS; i++) {
        double z = d[i] - c[0];
        double t = dot + z;
        c[0] = (t - dot) - z;
        dot = t;
      }
    }
  }

  return dot;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the dot product of 2 vectors: x.x
 *        using Kahan summation.
 *
 * \param[in]  n  size of arrays x and y
 * \param[in]  x  array of floating-point values
 *
 * \return  dot product
 */
/*----------------------------------------------------------------------------*/

static double
_cs_dot_xx_kahan(cs_lnum_t         n,
                 const cs_real_t  *x)
{
  double dot = 0.0;

# pragma omp parallel reduction(+:dot) if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(n, &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_lnum_t _nv = (_n/CS_VS)*CS_VS;
    const cs_real_t *_x = x + s_id;

    double d[CS_VS], c[CS_VS];

    for (cs_lnum_t i = 0; i < CS_VS; i ++) {
      d[i] = 0;
      c[i] = 0;
    }

    /* vector multiple */

    for (cs_lnum_t i = 0; i < _nv; i += CS_VS) {
      double z[CS_VS], t[CS_VS];
      for (cs_lnum_t j = 0; j < CS_VS; j++) {
        z[j] = (_x[i+j]* _x[i+j]) - c[j];
        t[j] = d[j] + z[j];
        c[j] = (t[j] - d[j]) - z[j];
        d[j] = t[j];
      }
    }

    /* vector remainder */

    for (cs_lnum_t i = _nv; i < _n; i++) {
      double z = (_x[i]* _x[i]) - c[0];
      double t = d[0] + z;
      c[0] = (t - d[0]) - z;
      d[0] = t;
    }

    /* Kahan sum on vector */

    {
      c[0] = 0;
      for (cs_lnum_t i = 0; i < CS_VS; i++) {
        double z = d[i] - c[0];
        double t = dot + z;
        c[0] = (t - dot) - z;
        dot = t;
      }
    }
  }

  return dot;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return 2 dot products of 2 vectors: x.x, and x.y
 *        using Kahan summation.
 *
 * The products could be computed separately, but computing them
 * simultaneously adds more optimization opportunities and possibly better
 * cache behavior.
 *
 * \param[in]   n   size of arrays x and y
 * \param[in]   x   array of floating-point values
 * \param[in]   y   array of floating-point values
 * \param[out]  xx  x.x dot product
 * \param[out]  xy  x.y dot product
 */
/*----------------------------------------------------------------------------*/

static void
_cs_dot_xx_xy_kahan(cs_lnum_t                    n,
                    const cs_real_t  *restrict   x,
                    const cs_real_t  *restrict   y,
                    double                      *xx,
                    double                      *xy)
{
  double dot_xx = 0.0, dot_xy = 0.0;

# pragma omp parallel reduction(+:dot_xx, dot_xy) if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(n, &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_real_t *_x = x + s_id;
    const cs_real_t *_y = y + s_id;

    double s[2] = {0, 0};
    double c[2] = {0, 0};

    for (cs_lnum_t i = 0; i < _n; i++) {
      double t[2];
      double z[2] = {(_x[i]* _x[i]) - c[0],
                     (_x[i]* _y[i]) - c[1]};
      for (int j = 0; j < 2; j++) {
        t[j] = s[j] + z[j];
        c[j] = (t[j] - s[j]) - z[j];
        s[j] = t[j];
      }
    }

    dot_xx = s[0];
    dot_xy = s[1];
  }

  *xx = dot_xx;
  *xy = dot_xy;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return 2 dot products of 3 vectors: x.x, and y.z
 *        using Kahan summation.
 *
 * The products could be computed separately, but computing them
 * simultaneously adds more optimization opportunities and possibly better
 * cache behavior.
 *
 * \param[in]   n   size of arrays x and y, and z
 * \param[in]   x   array of floating-point values
 * \param[in]   y   array of floating-point values
 * \param[in]   z   array of floating-point values
 * \param[out]  xy  x.y dot product
 * \param[out]  yz  y.z dot product
 */
/*----------------------------------------------------------------------------*/

static void
_cs_dot_xy_yz_kahan(cs_lnum_t                    n,
                    const cs_real_t  *restrict   x,
                    const cs_real_t  *restrict   y,
                    const cs_real_t  *restrict   z,
                    double                      *xy,
                    double                      *yz)
{
  double dot_xy = 0.0, dot_yz = 0.0;

# pragma omp parallel reduction(+:dot_xy, dot_yz) if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(n, &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_real_t *_x = x + s_id;
    const cs_real_t *_y = y + s_id;
    const cs_real_t *_z = z + s_id;

    double s[2] = {0, 0};
    double c[2] = {0, 0};

    for (cs_lnum_t i = 0; i < _n; i++) {
      double t[2];
      double u[2] = {(_x[i]*_y[i]) - c[0],
                     (_y[i]*_z[i]) - c[1]};
      for (int j = 0; j < 2; j++) {
        t[j] = s[j] + u[j];
        c[j] = (t[j] - s[j]) - u[j];
        s[j] = t[j];
      }
    }

    dot_xy = s[0];
    dot_xy = s[1];
  }

  *xy = dot_xy;
  *yz = dot_yz;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return 3 dot products of 3 vectors: x.x, x.y, and y.z
 *        using Kahan summation.

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

static void
_cs_dot_xx_xy_yz_kahan(cs_lnum_t                    n,
                       const cs_real_t  *restrict   x,
                       const cs_real_t  *restrict   y,
                       const cs_real_t  *restrict   z,
                       double                      *xx,
                       double                      *xy,
                       double                      *yz)
{
  double dot_xx = 0.0, dot_xy = 0.0, dot_yz = 0.0;

# pragma omp parallel reduction(+:dot_xx, dot_xy, dot_yz) if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(n, &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_real_t *_x = x + s_id;
    const cs_real_t *_y = y + s_id;
    const cs_real_t *_z = z + s_id;

    double s[3] = {0, 0, 0};
    double c[3] = {0, 0, 0};

    for (cs_lnum_t i = 0; i < _n; i++) {
      double t[3];
      double u[3] = {(_x[i]*_x[i]) - c[0],
                     (_x[i]*_y[i]) - c[1],
                     (_y[i]*_z[i]) - c[2]};
      for (int j = 0; j < 3; j++) {
        t[j] = s[j] + u[j];
        c[j] = (t[j] - s[j]) - u[j];
        s[j] = t[j];
      }
    }

    dot_xx = s[0];
    dot_xy = s[1];
    dot_yz = s[2];
  }

  *xx = dot_xx;
  *xy = dot_xy;
  *yz = dot_yz;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return 5 dot products of 3 vectors: x.x, y.y, x.y, x.z, and y.z
 *        using Kahan summation.
 *
 * The products could be computed separately, but computing them
 * simultaneously adds more optimization opportunities and possibly better
 * cache behavior.
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

static void
_cs_dot_xx_yy_xy_xz_yz_kahan(cs_lnum_t                    n,
                             const cs_real_t  *restrict   x,
                             const cs_real_t  *restrict   y,
                             const cs_real_t  *restrict   z,
                             double                      *xx,
                             double                      *yy,
                             double                      *xy,
                             double                      *xz,
                             double                      *yz)
{
  double dot_xx = 0.0, dot_yy = 0.0;
  double dot_xy = 0.0, dot_xz = 0.0, dot_yz = 0.0;

# pragma omp parallel reduction(+:dot_xx, dot_yy, dot_xy, dot_xz, dot_yz) \
                      if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(n, &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_real_t *_x = x + s_id;
    const cs_real_t *_y = y + s_id;
    const cs_real_t *_z = z + s_id;

    double s[5] = {0, 0, 0, 0, 0};
    double c[5] = {0, 0, 0, 0, 0};

    for (cs_lnum_t i = 0; i < _n; i++) {
      double t[5];
      double u[5] = {(_x[i]*_x[i]) - c[0],
                     (_y[i]*_y[i]) - c[1],
                     (_x[i]*_y[i]) - c[2],
                     (_x[i]*_z[i]) - c[3],
                     (_y[i]*_z[i]) - c[4]};
      for (int j = 0; j < 5; j++) {
        t[j] = s[j] + u[j];
        c[j] = (t[j] - s[j]) - u[j];
        s[j] = t[j];
      }
    }

    dot_xx = s[0];
    dot_yy = s[1];
    dot_xy = s[2];
    dot_xz = s[3];
    dot_yz = s[4];
  }

  *xx = dot_xx;
  *yy = dot_yy;
  *xy = dot_xy;
  *xz = dot_xz;
  *yz = dot_yz;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the global residual of 2 extensive vectors:
 *        1/sum(vol) . sum(X.Y.vol)
 *        using Kahan summation.
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

static double
_cs_gres_kahan(cs_lnum_t         n,
               const cs_real_t  *vol,
               const cs_real_t  *x,
               const cs_real_t  *y)
{
  double dot = 0.;
  double vtot = 0.;

# pragma omp parallel reduction(+:dot, vtot) if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(n, &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_real_t *_vol = vol + s_id;
    const cs_real_t *_x = x + s_id;
    const cs_real_t *_y = y + s_id;

    double s[2] = {0, 0};
    double c[2] = {0, 0};

    for (cs_lnum_t i = 0; i < _n; i++) {
      double t[2];
      double u[2] = {(_x[i]*_y[i]*_vol[i]) - c[0],
                     _vol[i]               - c[1]};
      for (int j = 0; j < 2; j++) {
        t[j] = s[j] + u[j];
        c[j] = (t[j] - s[j]) - u[j];
        s[j] = t[j];
      }
    }

    dot  = s[0];
    vtot = s[1];
  }

  double atot[2] = {dot, vtot};
  cs_parall_sum(2, CS_DOUBLE, atot);

  dot = atot[0] / atot[1];

  return dot;
}

/*============================================================================
 * Static global function pointers
 *============================================================================*/

static cs_dot_t       *_cs_glob_dot       = _cs_dot_superblock;
static cs_dot_xx_t    *_cs_glob_dot_xx    = _cs_dot_xx_superblock;
static cs_dot_xx_xy_t *_cs_glob_dot_xx_xy = _cs_dot_xx_xy_superblock;
static cs_dot_xy_yz_t *_cs_glob_dot_xy_yz = _cs_dot_xy_yz_superblock;
static cs_dot_xx_xy_yz_t  *_cs_glob_dot_xx_xy_yz
  = _cs_dot_xx_xy_yz_superblock;
static cs_dot_xx_yy_xy_xz_yz_t  *_cs_glob_dot_xx_yy_xy_xz_yz
  = _cs_dot_xx_yy_xy_xz_yz_superblock;

static cs_gres_t      *_cs_glob_gres      = _cs_gres_superblock;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the preferred BLAS reduction algorithm family.
 *
 * This may not be enforced for all algorithms, though it should at least
 * be enforced for the most general functions such as \ref cs_dot.
 *
 * \param[in]  mode   BLAS mode to use
 */
/*----------------------------------------------------------------------------*/

void
cs_blas_set_reduce_algorithm(cs_blas_reduce_t  mode)
{
  switch(mode) {
    case CS_BLAS_REDUCE_SUPERBLOCK:
      {
        _cs_glob_dot = _cs_dot_superblock;
        _cs_glob_dot_xx = _cs_dot_xx_superblock;
        _cs_glob_dot_xx_xy = _cs_dot_xx_xy_superblock;
        _cs_glob_dot_xy_yz = _cs_dot_xy_yz_superblock;
        _cs_glob_dot_xx_xy_yz = _cs_dot_xx_xy_yz_superblock;
        _cs_glob_dot_xx_yy_xy_xz_yz = _cs_dot_xx_yy_xy_xz_yz_superblock;
        _cs_glob_gres = _cs_gres_superblock;
      }
      break;
    case CS_BLAS_REDUCE_KAHAN:
      {
        _cs_glob_dot    = _cs_dot_kahan;
        _cs_glob_dot_xx = _cs_dot_xx_kahan;
        _cs_glob_dot_xx_xy = _cs_dot_xx_xy_kahan;
        _cs_glob_dot_xy_yz = _cs_dot_xy_yz_kahan;
        _cs_glob_dot_xx_xy_yz = _cs_dot_xx_xy_yz_kahan;
        _cs_glob_dot_xx_yy_xy_xz_yz = _cs_dot_xx_yy_xy_xz_yz_kahan;
        _cs_glob_gres = _cs_gres_kahan;
      }
      break;
  }
}

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
  if (n < 1)
    return;

# pragma omp parallel if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(n, &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_real_t *_x = x + s_id;
    cs_real_t *restrict _y = y + s_id;

    for (cs_lnum_t i = 0; i < _n; i++)
      _y[i] += (a * _x[i]);
  }
}

/*----------------------------------------------------------------------------
 * Return the sum of a vector. For better precision, a superblock algorithm
 * is used.
 *
 * \param[in]  n  size of array x
 * \param[in]  x  array of floating-point values
 *
 * \return the resulting sum
 *----------------------------------------------------------------------------*/

double
cs_sum(cs_lnum_t         n,
       const cs_real_t  *x)
{
  double sum = 0.0;

# pragma omp parallel reduction(+:sum) if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(n, &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_real_t *_x = x + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double sblk_sum = 0.0;

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
        cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        if (end_id > _n)
          end_id = _n;
        double blk_sum = 0.0;
        for (cs_lnum_t i = start_id; i < end_id; i++)
          blk_sum += _x[i];
        sblk_sum += blk_sum;
      }

      sum += sblk_sum;

    } /* Loop on super-blocks */

  } /* OpenMP block */

  return sum;
}

/*----------------------------------------------------------------------------
 * Return the weighted sum of a vector. For better precision, a superblock
 * algorithm is used.
 *
 * \param[in]  n  size of array x
 * \param[in]  w  array of floating-point weights
 * \param[in]  x  array of floating-point values
 *
 * \return the resulting weighted sum
 *----------------------------------------------------------------------------*/

double
cs_weighted_sum(cs_lnum_t         n,
                const cs_real_t  *w,
                const cs_real_t  *x)
{
  double wsum = 0.0;

# pragma omp parallel reduction(+:wsum) if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(n, &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_real_t *_x = x + s_id;
    const cs_real_t *_w = w + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double sblk_sum = 0.0;

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
        cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        if (end_id > _n)
          end_id = _n;
        double blk_sum = 0.0;
        for (cs_lnum_t i = start_id; i < end_id; i++)
          blk_sum += _w[i]*_x[i];
        sblk_sum += blk_sum;
      }

      wsum += sblk_sum;

    } /* Loop on super-blocks */

  } /* OpenMP block */

  return wsum;
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
  return _cs_glob_dot(n, x, y);
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
  return _cs_glob_dot_xx(n, x);
}

/*----------------------------------------------------------------------------
 * Return weighted dot products of a vector with itself: x.x
 *
 * For better precision, a superblock algorithm is used.
 *
 * \param[in]  n  size of array x
 * \param[in]  w  array of weights
 * \param[in]  x  array of floating-point values
 *
 * \return the resulting dot product
 *----------------------------------------------------------------------------*/

double
cs_dot_wxx(cs_lnum_t         n,
           const cs_real_t  *w,
           const cs_real_t  *x)
{
  double dot_wxx = 0.0;

# pragma omp parallel reduction(+:dot_wxx) if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(n, &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const cs_real_t *_x = x + s_id;
    const cs_real_t *_w = w + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double sdot_wxx = 0.0;

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
        cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        if (end_id > _n)
          end_id = _n;
        double cdot_wxx = 0.0;
        for (cs_lnum_t i = start_id; i < end_id; i++)
          cdot_wxx += _w[i]*_x[i]*_x[i];
        sdot_wxx += cdot_wxx;
      }

      dot_wxx += sdot_wxx;

    }

  }

  return dot_wxx;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return 2 dot products of 2 vectors: x.x, and x.y
 *
 * The products could be computed separately, but computing them
 * simultaneously adds more optimization opportunities and possibly better
 * cache behavior.
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
  _cs_glob_dot_xx_xy(n, x, y, xx, xy);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return 2 dot products of 3 vectors: x.y, and y.z
 *
 * The products could be computed separately, but computing them
 * simultaneously adds more optimization opportunities and possibly better
 * cache behavior.
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
  _cs_glob_dot_xy_yz(n, x, y, z, xy, yz);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return 3 dot products of 3 vectors: x.x, x.y, and y.z
 *
 * The products could be computed separately, but computing them
 * simultaneously adds more optimization opportunities and possibly better
 * cache behavior.
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
  _cs_glob_dot_xx_xy_yz(n, x, y, z, xx, xy, yz);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return 5 dot products of 3 vectors: x.x, y.y, x.y, x.z, and y.z
 *
 * The products could be computed separately, but computing them
 * simultaneously adds more optimization opportunities and possibly better
 * cache behavior.
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
  _cs_glob_dot_xx_yy_xy_xz_yz(n, x, y, z, xx, yy, xy, xz, yz);
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
/*!
 * \brief Return the global residual of 2 extensive vectors:
 *        1/sum(vol) . sum(X.Y.vol)
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
  return _cs_glob_gres(n, vol, x, y);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
