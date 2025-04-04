#ifndef __CS_ARRAY_REDUCE_H__
#define __CS_ARRAY_REDUCE_H__

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"
#include "base/cs_dispatch.h"
#include "base/cs_reducers.h"

/*----------------------------------------------------------------------------*/

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

#ifdef __cplusplus

/*----------------------------------------------------------------------------*/
/*
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
                       cs_real_t         &vmax);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute sums of an n-dimensional cs_real_t array's components.
 *
 * The array is interleaved.
 *
 * The algorithm here is similar to that used for BLAS.
 *
 * Template parmeters.
 *              stride  1 for scalars, 3 for vectors, 6 for symetric tensors
 *
 * \param[in]   ctx         cs_dispatch_context struct from higer level
 * \param[in]   n_elts      number of local elements
 * \param[in]   v_elt_list  optional list of parent elements on which values
 *                          are defined, or nullptr
 * \param[in]   v           pointer to array values
 * \param[out]  vsum        resulting strided sum array
 * */
/*----------------------------------------------------------------------------*/

template <size_t stride>
void
cs_array_reduce_sum_l(cs_dispatch_context ctx,
                      cs_lnum_t           n_elts,
                      const cs_lnum_t    *v_elt_list,
                      const cs_real_t     v[],
                      double              vsum[])
{
  struct cs_double_n<stride> rd;
  struct cs_reduce_sum_n<stride> reducer;

  /* If all values are defined on same list */
  if (v_elt_list == nullptr) {
    ctx.parallel_for_reduce(n_elts, rd, reducer,
      [=] CS_F_HOST_DEVICE (cs_lnum_t i, cs_double_n<stride> &res) {
      for (size_t k = 0; k < stride; k++)
        res.r[k] = v[stride*i + k];
    });
  }
  /* If values are defined on parent list */
  else {
    ctx.parallel_for_reduce(n_elts, rd, reducer,
      [=] CS_F_HOST_DEVICE (cs_lnum_t i, cs_double_n<stride> &res) {
      for (size_t k = 0; k < stride; k++)
        res.r[k] = v[stride*v_elt_list[i] + k];
    });
  }
  ctx.wait();

  for (size_t k = 0; k < stride; k++)
    vsum[k] = rd.r[k];
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
 * \param[in]   dim         local array dimension (max: 10)
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
                               double               vsum[]);
#endif // __cplusplus

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
 *                          are defined, or NULL
 * \param[in]   w_elt_list  optional list of parent elements on which weights
 *                          are defined, or NULL; if v_elt_list is defined
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
                       double            wsum[]);

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
 *                          are defined, or NULL
 * \param[in]   w_elt_list  optional list of parent elements on which weights
 *                          are defined, or NULL; if v_elt_list is defined
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
                                  double            wtot[]);

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
 *                          are defined, or NULL
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
                         cs_real_t         vmax[]);

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
 * parameters:
 *   ctx        <-- cs_dispatch_context struct from higher level
 *   n_elts     <-- number of local elements
 *   dim        <-- local array dimension (max: 9)
 *   v_elt_list <-- optional list of parent elements on which values
 *                  are defined, or NULL
 *   w_elt_list <-- optional list of parent elements on which weights
 *                  are defined, or NULL; if v_elt_list is defined
 *                  (ie. non-null),then w_elt_list = v_elt_list is assumed,
 *                  so this parameter is ignored
 *   v          <-- pointer to array values
 *   w          <-- pointer to weights
 *   vmin       --> resulting min array (size: dim, or 4 if dim = 3)
 *   vmax       --> resulting max array (size: dim, or 4 if dim = 3)
 *   vsum       --> resulting sum array (size: dim, or 4 if dim = 3)
 *   wsum       --> resulting weighted sum array (size: dim, or 4 if dim = 3)
 */
/*----------------------------------------------------------------------------*/

void
cs_array_reduce_simple_stats_l_w(cs_dispatch_context ctx,
                                 cs_lnum_t           n_elts,
                                 int                 dim,
                                 const cs_lnum_t    *v_elt_list,
                                 const cs_lnum_t    *w_elt_list,
                                 const cs_real_t     v[],
                                 const cs_real_t     w[],
                                 double              vmin[],
                                 double              vmax[],
                                 double              vsum[],
                                 double              wsum[]);

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
 *                  are defined, or NULL
 *   w_elt_list <-- optional list of parent elements on which weights
 *                  are defined, or NULL; if v_elt_list is defined
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
                               double            wssum[]);

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
 *                   are defined, or NULL
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
                                double             ssum[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ARRAY_REDUCE_H__ */
