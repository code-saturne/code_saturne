/*============================================================================
 * Functions computing square norms from CDO quantities
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_blas.h"
#include "cs_math.h"
#include "cs_parall.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_blas.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

#define _dp3  cs_math_3_dot_product

/* Block size for superblock algorithm */

#define CS_SBLOCK_BLOCK_SIZE 60

/* Cache line multiple, in cs_real_t units */

#define CS_CL  (CS_CL_SIZE/8)

/*=============================================================================
 * Local static variables
 *============================================================================*/

/* Pointer to shared structures (owned by a cs_domain_t structure) */
static const cs_cdo_quantities_t  *cs_cdo_quant;
static const cs_cdo_connect_t  *cs_cdo_connect;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Sanity checkings before computing norms
 *
 * \param[in]  func_name   name of the calling function
 * \param[in]  c2x         first pointer to check
 * \param[in]  w_c2x       second pointer to check
 */
/*----------------------------------------------------------------------------*/

static inline void
_sanity_checks(const char              func_name[],
               const cs_adjacency_t   *c2x,
               const cs_real_t        *w_c2x)
{
  assert(cs_cdo_quant != NULL && cs_cdo_connect != NULL);

  if (c2x == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s: The cs_adjacency_t structure is not allocated.\n",
              func_name);

  if (w_c2x == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s: The array storing weights is not allocated.\n",
              func_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute array index bounds for a local thread.
 *        When called inside an OpenMP parallel section, this will return the
 *        start an past-the-end indexes for the array range assigned to that
 *        thread. In other cases, the start index is 1, and the past-the-end
 *        index is n;
 *
 * \param[in]       n     size of array
 * \param[in, out]  s_id  start index for the current thread
 * \param[in, out]  e_id  past-the-end index for the current thread
 */
/*----------------------------------------------------------------------------*/

static inline void
_thread_range(cs_lnum_t    n,
              cs_lnum_t   *s_id,
              cs_lnum_t   *e_id)
{
#if defined(HAVE_OPENMP)
  const int t_id = omp_get_thread_num();
  const int n_t = omp_get_num_threads();
  const cs_lnum_t t_n = (n + n_t - 1) / n_t;

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
 * \brief  Compute the weighted square norm of an array relying on a weight
 *         which is accessed by an index. Case of a scalar-valued array.
 *
 * \param[in]  size      size of the weight array
 * \param[in]  c2x       pointer to the adjacency used to access the weights
 * \param[in]  w_c2x     weight array
 * \param[in]  array     array to analyze
 * \param[in]  do_redux  perform the parallel reduction ?
 *
 * \return the square weighted L2-norm
 */
/*----------------------------------------------------------------------------*/

static double
_c2x_scalar_sqnorm(const cs_lnum_t         size,
                   const cs_adjacency_t   *c2x,
                   const cs_real_t        *w_c2x,
                   const cs_real_t        *array,
                   bool                    do_redux)
{
 /*
  * The algorithm used is l3superblock60, based on the article:
  * "Reducing Floating Point Error in Dot Product Using the Superblock Family
  * of Algorithms" by Anthony M. Castaldo, R. Clint Whaley, and Anthony
  * T. Chronopoulos, SIAM J. SCI. COMPUT., Vol. 31, No. 2, pp. 1156--1174
  * 2008 Society for Industrial and Applied Mathematics
  */

  double  l2norm = 0;

# pragma omp parallel reduction(+:l2norm) if (size > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(size, &s_id, &e_id);

    const cs_lnum_t  n = e_id - s_id;
    const cs_lnum_t  *_ids = c2x->ids + s_id;
    const cs_real_t  *_w = w_c2x + s_id;
    const cs_lnum_t  block_size = CS_SBLOCK_BLOCK_SIZE;
    const cs_lnum_t  n_blocks = (n + block_size - 1) / block_size;
    const cs_lnum_t  n_sblocks = (n_blocks > 3) ? sqrt(n_blocks) : 1;
    const cs_lnum_t  blocks_in_sblocks =
      (n + block_size*n_sblocks - 1) / (block_size*n_sblocks);

    cs_lnum_t  shift = 0;

    for (cs_lnum_t s = 0; s < n_sblocks; s++) { /* Loop on slices */

      double  s_l2norm = 0.0;

      for (cs_lnum_t b_id = 0; b_id < blocks_in_sblocks; b_id++) {

        const cs_lnum_t  start_id = shift;
        shift += block_size;
        if (shift > n)
          shift = n, b_id = blocks_in_sblocks;
        const cs_lnum_t  end_id = shift;

        double  _l2norm = 0.0;
        for (cs_lnum_t j = start_id; j < end_id; j++) {

          const cs_real_t  v = array[_ids[j]];

          _l2norm += _w[j] * v*v;

        } /* Loop on block_size */

        s_l2norm += _l2norm;

      } /* Loop on blocks */

      l2norm += s_l2norm;

    } /* Loop on super-blocks */

  } /* OpenMP block */

  /* Parallel treatment */

  if (do_redux)
    cs_parall_sum(1, CS_DOUBLE, &l2norm);

  return l2norm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the weighted square norm of an array relying on a weight
 *         which is accessed by an index. Case of a vector-valued array.
 *
 * \param[in]  size    size of the weight array
 * \param[in]  c2x     pointer to the adjacency used to access the weights
 * \param[in]  w_c2x   weight array
 * \param[in]  array   array to analyze (vector-valued)
 *
 * \return the square weighted L2-norm
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_c2x_vector_sqnorm(const cs_lnum_t         size,
                   const cs_adjacency_t   *c2x,
                   const cs_real_t        *w_c2x,
                   const cs_real_t        *array)
{
 /*
  * The algorithm used is l3superblock60, based on the article:
  * "Reducing Floating Point Error in Dot Product Using the Superblock Family
  * of Algorithms" by Anthony M. Castaldo, R. Clint Whaley, and Anthony
  * T. Chronopoulos, SIAM J. SCI. COMPUT., Vol. 31, No. 2, pp. 1156--1174
  * 2008 Society for Industrial and Applied Mathematics
  */

  double  l2norm = 0;

# pragma omp parallel reduction(+:l2norm) if (size > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(size, &s_id, &e_id);

    const cs_lnum_t  n = e_id - s_id;
    const cs_lnum_t  *_ids = c2x->ids + s_id;
    const cs_real_t  *_w = w_c2x + s_id;
    const cs_lnum_t  block_size = CS_SBLOCK_BLOCK_SIZE;
    const cs_lnum_t  n_blocks = (n + block_size - 1) / block_size;
    const cs_lnum_t  n_sblocks = (n_blocks > 3) ? sqrt(n_blocks) : 1;
    const cs_lnum_t  blocks_in_sblocks =
      (n + block_size*n_sblocks - 1) / (block_size*n_sblocks);

    cs_lnum_t  shift = 0;

    for (cs_lnum_t s = 0; s < n_sblocks; s++) { /* Loop on slices */

      double  s_l2norm = 0.0;

      for (cs_lnum_t b_id = 0; b_id < blocks_in_sblocks; b_id++) {

        const cs_lnum_t  start_id = shift;
        shift += block_size;
        if (shift > n)
          shift = n, b_id = blocks_in_sblocks;
        const cs_lnum_t  end_id = shift;

        double  _l2norm = 0.0;
        for (cs_lnum_t j = start_id; j < end_id; j++) {

          const cs_real_t  *v = array + 3*_ids[j];

          _l2norm += _w[j] * cs_math_3_square_norm(v);

        } /* Loop on block_size */

        s_l2norm += _l2norm;

      } /* Loop on blocks */

      l2norm += s_l2norm;

    } /* Loop on super-blocks */

  } /* OpenMP block */

  /* Parallel treatment */

  cs_parall_sum(1, CS_REAL_TYPE, &l2norm);

  return (cs_real_t)l2norm;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers to main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_blas_init_sharing(const cs_cdo_quantities_t    *quant,
                         const cs_cdo_connect_t       *connect)
{
  /* Assign static const pointers */

  cs_cdo_quant = quant;
  cs_cdo_connect = connect;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the square norm of an array
 *         Case of a scalar-valued array defined as a potential at primal
 *         cells. Thus, the weigth is the cell volume. The computed quantities
 *         are synchronized in parallel.
 *
 * \param[in]  array   array to analyze
 *
 * \return the square weighted L2-norm
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdo_blas_square_norm_pcsp(const cs_real_t        *array)
{
  const cs_lnum_t  n_cells = cs_cdo_quant->n_cells;

 /*
  * The algorithm used is l3superblock60, based on the article:
  * "Reducing Floating Point Error in Dot Product Using the Superblock Family
  * of Algorithms" by Anthony M. Castaldo, R. Clint Whaley, and Anthony
  * T. Chronopoulos, SIAM J. SCI. COMPUT., Vol. 31, No. 2, pp. 1156--1174
  * 2008 Society for Industrial and Applied Mathematics
  */

  double  num = 0.;

# pragma omp parallel reduction(+:num) if (n_cells > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(n_cells, &s_id, &e_id);

    const cs_lnum_t  n = e_id - s_id;
    const cs_lnum_t  block_size = CS_SBLOCK_BLOCK_SIZE;
    const cs_lnum_t  n_blocks = (n + block_size - 1) / block_size;
    const cs_lnum_t  n_sblocks = (n_blocks > 3) ? sqrt(n_blocks) : 1;
    const cs_lnum_t  blocks_in_sblocks =
      (n + block_size*n_sblocks - 1) / (block_size*n_sblocks);

    cs_lnum_t  shift = 0;

    for (cs_lnum_t s = 0; s < n_sblocks; s++) { /* Loop on slices */

      double s_num = 0.0;

      for (cs_lnum_t b_id = 0; b_id < blocks_in_sblocks; b_id++) {

        const cs_lnum_t  start_id = shift;
        shift += block_size;
        if (shift > n)
          shift = n, b_id = blocks_in_sblocks;
        const cs_lnum_t  end_id = shift;

        double _num = 0.0;
        for (cs_lnum_t j = start_id; j < end_id; j++)
          _num += cs_cdo_quant->cell_vol[j] * array[j]*array[j];

        s_num += _num;

      } /* Loop on blocks */

      num += s_num;

    } /* Loop on super-blocks */

  } /* OpenMP block */

  /* Parallel treatment */

  if (cs_glob_n_ranks > 1) {

    cs_real_t  sum = num;
    cs_parall_sum(1, CS_REAL_TYPE, &sum);
    num = sum;

  }

  return (cs_real_t)num;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the norm ||b - a||**2
 *         Case of two scalar-valued arrays a and b defined as a potential at
 *         primal cells. Thus, the weigth is the cell volume. The computed
 *         quantities are synchronized in parallel.
 *
 * \param[in]  a   first array
 * \param[in]  b   second array
 *
 * \return the value  of ||b - a||**2
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdo_blas_square_norm_pcsp_diff(const cs_real_t        *a,
                                  const cs_real_t        *b)
{
  const cs_lnum_t  n_cells = cs_cdo_quant->n_cells;

 /*
  * The algorithm used is l3superblock60, based on the article:
  * "Reducing Floating Point Error in Dot Product Using the Superblock Family
  * of Algorithms" by Anthony M. Castaldo, R. Clint Whaley, and Anthony
  * T. Chronopoulos, SIAM J. SCI. COMPUT., Vol. 31, No. 2, pp. 1156--1174
  * 2008 Society for Industrial and Applied Mathematics
  */

  double  num = 0.;

# pragma omp parallel reduction(+:num) if (n_cells > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(n_cells, &s_id, &e_id);

    const cs_lnum_t  n = e_id - s_id;
    const cs_lnum_t  block_size = CS_SBLOCK_BLOCK_SIZE;
    const cs_lnum_t  n_blocks = (n + block_size - 1) / block_size;
    const cs_lnum_t  n_sblocks = (n_blocks > 3) ? sqrt(n_blocks) : 1;
    const cs_lnum_t  blocks_in_sblocks =
      (n + block_size*n_sblocks - 1) / (block_size*n_sblocks);

    cs_lnum_t  shift = 0;

    for (cs_lnum_t s = 0; s < n_sblocks; s++) { /* Loop on slices */

      double s_num = 0.0;

      for (cs_lnum_t b_id = 0; b_id < blocks_in_sblocks; b_id++) {

        const cs_lnum_t  start_id = shift;
        shift += block_size;
        if (shift > n)
          shift = n, b_id = blocks_in_sblocks;
        const cs_lnum_t  end_id = shift;

        double _num = 0.0;
        for (cs_lnum_t j = start_id; j < end_id; j++)
          _num += cs_cdo_quant->cell_vol[j] * (b[j] - a[j])*(b[j] - a[j]);

        s_num += _num;

      } /* Loop on blocks */

      num += s_num;

    } /* Loop on super-blocks */

  } /* OpenMP block */

  /* Parallel treatment */

  if (cs_glob_n_ranks > 1) {

    cs_real_t  sums = num;
    cs_parall_sum(1, CS_REAL_TYPE, &sums);
    num = sums;

  }

  return (cs_real_t)num;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the norm  ||a - ref||**2 / || ref||**2
 *         Case of two scalar-valued arrays a and ref defined as a potential at
 *         primal cells. Thus, the weigth is the cell volume. The computed
 *         quantities are synchronized in parallel. "ndiff" stands for
 *         "normalized difference"
 *
 * \param[in]  a     array to analyze
 * \param[in]  ref   array used for normalization and difference
 *
 * \return the normalized square weighted L2-norm of the difference between the
 *         two arrays
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdo_blas_square_norm_pcsp_ndiff(const cs_real_t        *a,
                                   const cs_real_t        *ref)
{
  const cs_lnum_t  n_cells = cs_cdo_quant->n_cells;

 /*
  * The algorithm used is l3superblock60, based on the article:
  * "Reducing Floating Point Error in Dot Product Using the Superblock Family
  * of Algorithms" by Anthony M. Castaldo, R. Clint Whaley, and Anthony
  * T. Chronopoulos, SIAM J. SCI. COMPUT., Vol. 31, No. 2, pp. 1156--1174
  * 2008 Society for Industrial and Applied Mathematics
  */

  double  num = 0., denum = 0.;

# pragma omp parallel reduction(+:num) reduction(+:denum) \
  if (n_cells > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(n_cells, &s_id, &e_id);

    const cs_lnum_t  n = e_id - s_id;
    const cs_lnum_t  block_size = CS_SBLOCK_BLOCK_SIZE;
    const cs_lnum_t  n_blocks = (n + block_size - 1) / block_size;
    const cs_lnum_t  n_sblocks = (n_blocks > 3) ? sqrt(n_blocks) : 1;
    const cs_lnum_t  blocks_in_sblocks =
      (n + block_size*n_sblocks - 1) / (block_size*n_sblocks);

    cs_lnum_t  shift = 0;

    for (cs_lnum_t s = 0; s < n_sblocks; s++) { /* Loop on slices */

      double s_num = 0.0, s_denum = 0.0;

      for (cs_lnum_t b_id = 0; b_id < blocks_in_sblocks; b_id++) {

        const cs_lnum_t  start_id = shift;
        shift += block_size;
        if (shift > n)
          shift = n, b_id = blocks_in_sblocks;
        const cs_lnum_t  end_id = shift;

        double _num = 0.0, _denum = 0.0;
        for (cs_lnum_t j = start_id; j < end_id; j++) {

          const cs_real_t  vol_c = cs_cdo_quant->cell_vol[j];

          _num += vol_c * (a[j] - ref[j])*(a[j] - ref[j]);
          _denum += vol_c * ref[j] * ref[j];

        } /* Loop on block_size */

        s_num += _num;
        s_denum += _denum;

      } /* Loop on blocks */

      num += s_num;
      denum += s_denum;

    } /* Loop on super-blocks */

  } /* OpenMP block */

  /* Parallel treatment */

  if (cs_glob_n_ranks > 1) {

    cs_real_t  sums[2] = {num, denum};
    cs_parall_sum(2, CS_REAL_TYPE, sums);
    num = sums[0], denum = sums[1];

  }

  if (fabs(denum) > cs_math_zero_threshold)
    num /= denum;

  return (cs_real_t)num;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the dot product of two arrays using the classical Euclidean
 *         dot product (without weight).
 *         Case of a scalar-valued arrays defined at primal vertices.
 *         The computed quantity is synchronized in parallel.
 *
 * \param[in]  a   first array to analyze
 * \param[in]  b   second array to analyze
 *
 * \return the value of the dot product
 */
/*----------------------------------------------------------------------------*/

double
cs_cdo_blas_dotprod_vertex(const cs_real_t        *a,
                           const cs_real_t        *b)
{
  return cs_gdot(cs_cdo_quant->n_vertices, a, b);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the square norm of an array using an Euclidean 2-norm.
 *         Case of a scalar-valued array defined at primal vertices.
 *         The computed quantities are synchronized in parallel.
 *
 * \param[in]  array   array to analyze
 *
 * \return the square weighted L2-norm
 */
/*----------------------------------------------------------------------------*/

double
cs_cdo_blas_square_norm_vertex(const cs_real_t        *array)
{
  double  retval = cs_dot_xx(cs_cdo_quant->n_vertices, array);

  cs_parall_sum(1, CS_DOUBLE, &retval);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the dot product of two arrays using a weighted Euclidean dot
 *         product relying on CDO quantities.
 *         Case of a scalar-valued arrays defined as a potential at primal
 *         vertices. Thus, the weigth is the portion of dual cell (associated
 *         to a primal vertex) inside a primal cell.  The computed quantity is
 *         synchronized in parallel.
 *
 * \param[in]  a   first array to analyze
 * \param[in]  b   second array to analyze
 *
 * \return the value of the dot product
 */
/*----------------------------------------------------------------------------*/

double
cs_cdo_blas_dotprod_pvsp(const cs_real_t        *a,
                         const cs_real_t        *b)
{
  const cs_adjacency_t  *c2v = cs_cdo_connect->c2v;
  const cs_real_t  *w_c2v = cs_cdo_quant->dcell_vol;

  _sanity_checks(__func__, c2v, w_c2v);

  const cs_lnum_t  size = c2v->idx[cs_cdo_quant->n_cells];

  /*
   * The algorithm used is l3superblock60, based on the article:
   * "Reducing Floating Point Error in Dot Product Using the Superblock Family
   * of Algorithms" by Anthony M. Castaldo, R. Clint Whaley, and Anthony
   * T. Chronopoulos, SIAM J. SCI. COMPUT., Vol. 31, No. 2, pp. 1156--1174
   * 2008 Society for Industrial and Applied Mathematics
   */

  double  dp = 0;

# pragma omp parallel reduction(+:dp) if (size > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(size, &s_id, &e_id);

    const cs_lnum_t  n = e_id - s_id;
    const cs_lnum_t  *_ids = c2v->ids + s_id;
    const cs_real_t  *_w = w_c2v + s_id;
    const cs_lnum_t  block_size = CS_SBLOCK_BLOCK_SIZE;
    const cs_lnum_t  n_blocks = (n + block_size - 1) / block_size;
    const cs_lnum_t  n_sblocks = (n_blocks > 3) ? sqrt(n_blocks) : 1;
    const cs_lnum_t  blocks_in_sblocks =
      (n + block_size*n_sblocks - 1) / (block_size*n_sblocks);

    cs_lnum_t  shift = 0;

    for (cs_lnum_t s = 0; s < n_sblocks; s++) { /* Loop on slices */

      double  s_dp = 0.0;

      for (cs_lnum_t b_id = 0; b_id < blocks_in_sblocks; b_id++) {

        const cs_lnum_t  start_id = shift;
        shift += block_size;
        if (shift > n)
          shift = n, b_id = blocks_in_sblocks;
        const cs_lnum_t  end_id = shift;

        double  _dp = 0.0;
        for (cs_lnum_t j = start_id; j < end_id; j++) {

          const cs_lnum_t  id = _ids[j];

          _dp += _w[j] * a[id]*b[id];

        } /* Loop on block_size */

        s_dp += _dp;

      } /* Loop on blocks */

      dp += s_dp;

    } /* Loop on super-blocks */

  } /* OpenMP block */

  /* Parallel treatment */

  cs_parall_sum(1, CS_REAL_TYPE, &dp);

  return (cs_real_t)dp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the square norm of an array
 *         Case of a scalar-valued array defined as a potential at primal
 *         vertices. Thus, the weigth is the portion of dual cell inside each
 *         (primal cell). The computed quantities are synchronized in parallel.
 *
 * \param[in]  array   array to analyze
 *
 * \return the square weighted L2-norm
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdo_blas_square_norm_pvsp(const cs_real_t        *array)
{
  const cs_adjacency_t  *c2v = cs_cdo_connect->c2v;
  const cs_real_t  *w_c2v = cs_cdo_quant->dcell_vol;

  _sanity_checks(__func__, c2v, w_c2v);

  return _c2x_scalar_sqnorm(c2v->idx[cs_cdo_quant->n_cells], c2v, w_c2v,
                            array, true); /* do parallel sum */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the norm ||b - a||**2
 *         Case of two scalar-valued arrays a and b defined as a potential at
 *         primal vertices. Thus, the weigth is the portion of dual cell in a
 *         primal cell. The computed quantities are synchronized in parallel.
 *
 * \param[in]  a   first array
 * \param[in]  b   second array
 *
 * \return the value  of ||b - a||**2
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdo_blas_square_norm_pvsp_diff(const cs_real_t        *a,
                                  const cs_real_t        *b)
{
  const cs_adjacency_t  *c2v = cs_cdo_connect->c2v;
  const cs_real_t  *w_c2v = cs_cdo_quant->dcell_vol;

  _sanity_checks(__func__, c2v, w_c2v);

  const cs_lnum_t  size = c2v->idx[cs_cdo_quant->n_cells];

 /*
  * The algorithm used is l3superblock60, based on the article:
  * "Reducing Floating Point Error in Dot Product Using the Superblock Family
  * of Algorithms" by Anthony M. Castaldo, R. Clint Whaley, and Anthony
  * T. Chronopoulos, SIAM J. SCI. COMPUT., Vol. 31, No. 2, pp. 1156--1174
  * 2008 Society for Industrial and Applied Mathematics
  */

  double  num = 0.;

# pragma omp parallel reduction(+:num) if (size > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(size, &s_id, &e_id);

    const cs_lnum_t  n = e_id - s_id;
    const cs_lnum_t  *_ids = c2v->ids + s_id;
    const cs_real_t  *_w = w_c2v + s_id;
    const cs_lnum_t  block_size = CS_SBLOCK_BLOCK_SIZE;
    const cs_lnum_t  n_blocks = (n + block_size - 1) / block_size;
    const cs_lnum_t  n_sblocks = (n_blocks > 3) ? sqrt(n_blocks) : 1;
    const cs_lnum_t  blocks_in_sblocks =
      (n + block_size*n_sblocks - 1) / (block_size*n_sblocks);

    cs_lnum_t  shift = 0;

    for (cs_lnum_t s = 0; s < n_sblocks; s++) { /* Loop on slices */

      double s_num = 0.0;

      for (cs_lnum_t b_id = 0; b_id < blocks_in_sblocks; b_id++) {

        const cs_lnum_t  start_id = shift;
        shift += block_size;
        if (shift > n)
          shift = n, b_id = blocks_in_sblocks;
        const cs_lnum_t  end_id = shift;

        double _num = 0.0;
        for (cs_lnum_t j = start_id; j < end_id; j++) {

          const cs_lnum_t  id = _ids[j];
          const cs_real_t  diff = b[id] - a[id];

          _num += _w[j] * diff * diff;

        }

        s_num += _num;

      } /* Loop on blocks */

      num += s_num;

    } /* Loop on super-blocks */

  } /* OpenMP block */

  /* Parallel treatment */

  cs_parall_sum(1, CS_REAL_TYPE, &num);

  return (cs_real_t)num;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the square norm of an array
 *         Case of a non-interlaced scalar-valued array of stride = 2 defined as
 *         a potential at primal vertices. Thus, the weigth is the portion of
 *         dual cell (associated to a primal vertex) inside a primal cell. The
 *         computed quantity is synchronized in parallel.
 *
 * \param[in]  array   array to analyze
 *
 * \return the square weighted L2-norm
 */
/*----------------------------------------------------------------------------*/

double
cs_cdo_blas_square_norm_2pvsp(const cs_real_t        *array)
{
  const cs_adjacency_t  *c2v = cs_cdo_connect->c2v;
  const cs_real_t  *w_c2v = cs_cdo_quant->dcell_vol;

  _sanity_checks(__func__, c2v, w_c2v);

  double res = 0;

  /* Avoid two parallel sums --> parameter is equal to "false" */

  res = _c2x_scalar_sqnorm(c2v->idx[cs_cdo_quant->n_cells], c2v, w_c2v,
                           array, false); /* do not parallel sum */

  res += _c2x_scalar_sqnorm(c2v->idx[cs_cdo_quant->n_cells], c2v, w_c2v,
                            array + cs_cdo_quant->n_vertices, false);

  cs_parall_sum(1, CS_REAL_TYPE, &res);

  return res;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the dot product of two arrays using a weighted Euclidean dot
 *         product relying on CDO quantities.
 *         Case of non-interlaced scalar-valued arrays of stride = 2 defined as
 *         a potential at primal vertices. Thus, the weigth is the portion of
 *         dual cell (associated to a primal vertex) inside a primal cell. The
 *         computed quantity is synchronized in parallel.
 *
 * \param[in]  a   first array to analyze
 * \param[in]  b   second array to analyze
 *
 * \return the value of the dot product
 */
/*----------------------------------------------------------------------------*/

double
cs_cdo_blas_dotprod_2pvsp(const cs_real_t        *a,
                          const cs_real_t        *b)
{
  const cs_adjacency_t  *c2v = cs_cdo_connect->c2v;
  const cs_real_t  *w_c2v = cs_cdo_quant->dcell_vol;

  _sanity_checks(__func__, c2v, w_c2v);

  const cs_lnum_t  n_vertices = cs_cdo_quant->n_vertices;
  const cs_lnum_t  size = c2v->idx[cs_cdo_quant->n_cells];

  const cs_real_t  *a1 = a, *a2 = a1 + n_vertices;
  const cs_real_t  *b1 = b, *b2 = b1 + n_vertices;

  /*
   * The algorithm used is l3superblock60, based on the article:
   * "Reducing Floating Point Error in Dot Product Using the Superblock Family
   * of Algorithms" by Anthony M. Castaldo, R. Clint Whaley, and Anthony
   * T. Chronopoulos, SIAM J. SCI. COMPUT., Vol. 31, No. 2, pp. 1156--1174
   * 2008 Society for Industrial and Applied Mathematics
   */

  double  dp = 0;

# pragma omp parallel reduction(+:dp) if (size > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(size, &s_id, &e_id);

    const cs_lnum_t  n = e_id - s_id;
    const cs_lnum_t  *_ids = c2v->ids + s_id;
    const cs_real_t  *_w = w_c2v + s_id;
    const cs_lnum_t  block_size = CS_SBLOCK_BLOCK_SIZE;
    const cs_lnum_t  n_blocks = (n + block_size - 1) / block_size;
    const cs_lnum_t  n_sblocks = (n_blocks > 3) ? sqrt(n_blocks) : 1;
    const cs_lnum_t  blocks_in_sblocks =
      (n + block_size*n_sblocks - 1) / (block_size*n_sblocks);

    cs_lnum_t  shift = 0;

    for (cs_lnum_t s = 0; s < n_sblocks; s++) { /* Loop on slices */

      double  s_dp = 0.0;

      for (cs_lnum_t b_id = 0; b_id < blocks_in_sblocks; b_id++) {

        const cs_lnum_t  start_id = shift;
        shift += block_size;
        if (shift > n)
          shift = n, b_id = blocks_in_sblocks;
        const cs_lnum_t  end_id = shift;

        double  _dp = 0.0;
        for (cs_lnum_t j = start_id; j < end_id; j++) {

          const cs_lnum_t  id = _ids[j];

          _dp += _w[j] * (a1[id]*b1[id] + a2[id]*b2[id]);

        } /* Loop on block_size */

        s_dp += _dp;

      } /* Loop on blocks */

      dp += s_dp;

    } /* Loop on super-blocks */

  } /* OpenMP block */

  /* Parallel treatment */

  cs_parall_sum(1, CS_DOUBLE, &dp);

  return dp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the dot product of two arrays using the classical Euclidean
 *         dot product (without weight).
 *         Case of a scalar-valued arrays defined at primal faces.
 *         The computed quantity is synchronized in parallel.
 *
 * \param[in]  a   first array to analyze
 * \param[in]  b   second array to analyze
 *
 * \return the value of the dot product
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdo_blas_dotprod_face(const cs_real_t        *a,
                         const cs_real_t        *b)
{
  return cs_gdot(cs_cdo_quant->n_faces, a, b);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the square norm of an array using an Euclidean 2-norm.
 *         Case of a scalar-valued array defined at primal faces.
 *         The computed quantities are synchronized in parallel.
 *
 * \param[in]  array   array to analyze
 *
 * \return the square weighted L2-norm
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdo_blas_square_norm_face(const cs_real_t        *array)
{
  cs_real_t  retval = cs_dot_xx(cs_cdo_quant->n_faces, array);

  cs_parall_sum(1, CS_DOUBLE, &retval);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the square norm of an array
 *         Case of a scalar-valued array defined as a potential at primal
 *         faces. Thus, the weigth is the pyramid of apex the cell center and
 *         of basis the face. The computed quantities are synchronized in
 *         parallel.
 *
 * \param[in]  array   array to analyze
 *
 * \return the square weighted L2-norm
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdo_blas_square_norm_pfsp(const cs_real_t        *array)
{
  const cs_adjacency_t  *c2f = cs_cdo_connect->c2f;
  const cs_real_t  *w_c2f = cs_cdo_quant->pvol_fc;

  _sanity_checks(__func__, c2f, w_c2f);

  return _c2x_scalar_sqnorm(c2f->idx[cs_cdo_quant->n_cells], c2f, w_c2f,
                            array, true); /* do parallel sum */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the square norm of an array
 *         Case of a vector-valued array defined as a potential at primal
 *         faces. Thus, the weigth is the pyramid of apex the cell center and
 *         of basis the face. The computed quantities are synchronized in
 *         parallel.
 *
 * \param[in]  array   array to analyze (vector-valued)
 *
 * \return the square weighted L2-norm
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdo_blas_square_norm_pfvp(const cs_real_t        *array)
{
  const cs_adjacency_t  *c2f = cs_cdo_connect->c2f;
  const cs_real_t  *w_c2f = cs_cdo_quant->pvol_fc;

  _sanity_checks(__func__, c2f, w_c2f);

  return _c2x_vector_sqnorm(c2f->idx[cs_cdo_quant->n_cells], c2f, w_c2f,
                            array);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the dot product of two arrays using a weighted Euclidean
 *         dot product relying on CDO quantities.
 *         Case of a scalar-valued arrays defined as a flux at primal
 *         faces. Thus, the weigth is the pyramid of apex the cell center and
 *         of basis the face. Each face quantity is normalized by the face
 *         surface. The computed quantity is synchronized in parallel.
 *
 * \param[in]  a   first array to analyze
 * \param[in]  b   second array to analyze
 *
 * \return the value of the dot product
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdo_blas_dotprod_pfsf(const cs_real_t        *a,
                         const cs_real_t        *b)
{
  const cs_adjacency_t  *c2f = cs_cdo_connect->c2f;
  const cs_real_t  *w_c2f = cs_cdo_quant->pvol_fc;
  const cs_real_t  *i_surf = cs_cdo_quant->i_face_surf;
  const cs_real_t  *b_surf = cs_cdo_quant->b_face_surf;
  const cs_lnum_t  size = c2f->idx[cs_cdo_quant->n_cells];
  const cs_lnum_t  n_i_faces = cs_cdo_quant->n_i_faces;

  _sanity_checks(__func__, c2f, w_c2f);

  /*
   * The algorithm used is l3superblock60, based on the article:
   * "Reducing Floating Point Error in Dot Product Using the Superblock Family
   * of Algorithms" by Anthony M. Castaldo, R. Clint Whaley, and Anthony
   * T. Chronopoulos, SIAM J. SCI. COMPUT., Vol. 31, No. 2, pp. 1156--1174
   * 2008 Society for Industrial and Applied Mathematics
   */

  double  dp = 0;

# pragma omp parallel reduction(+:dp) if (size > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(size, &s_id, &e_id);

    const cs_lnum_t  n = e_id - s_id;
    const cs_lnum_t  *_ids = c2f->ids + s_id;
    const cs_real_t  *_w = w_c2f + s_id;
    const cs_lnum_t  block_size = CS_SBLOCK_BLOCK_SIZE;
    const cs_lnum_t  n_blocks = (n + block_size - 1) / block_size;
    const cs_lnum_t  n_sblocks = (n_blocks > 3) ? sqrt(n_blocks) : 1;
    const cs_lnum_t  blocks_in_sblocks =
      (n + block_size*n_sblocks - 1) / (block_size*n_sblocks);

    cs_lnum_t  shift = 0;

    for (cs_lnum_t s = 0; s < n_sblocks; s++) { /* Loop on slices */

      double  s_dp = 0.0;

      for (cs_lnum_t b_id = 0; b_id < blocks_in_sblocks; b_id++) {

        const cs_lnum_t  start_id = shift;
        shift += block_size;
        if (shift > n)
          shift = n, b_id = blocks_in_sblocks;
        const cs_lnum_t  end_id = shift;

        double  _dp = 0.0;
        for (cs_lnum_t j = start_id; j < end_id; j++) {

          const cs_lnum_t  id = _ids[j];
          const cs_real_t  osurf =
            (id < n_i_faces) ? 1./i_surf[id] : 1./b_surf[id - n_i_faces];
          const cs_real_t  va = a[id]*osurf, vb = b[id]*osurf;

          _dp += _w[j] * va*vb;

        } /* Loop on block_size */

        s_dp += _dp;

      } /* Loop on blocks */

      dp += s_dp;

    } /* Loop on super-blocks */

  } /* OpenMP block */

  /* Parallel treatment */

  cs_parall_sum(1, CS_REAL_TYPE, &dp);

  return (cs_real_t)dp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the square norm of an array
 *         Case of a scalar-valued array defined as a flux at primal
 *         faces. Thus, the weigth is the pyramid of apex the cell center and
 *         of basis the face. Each face quantity is normalized by the face
 *         surface. The computed quantities are synchronized in parallel.
 *
 * \param[in]  array   array to analyze
 *
 * \return the square weighted L2-norm
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdo_blas_square_norm_pfsf(const cs_real_t        *array)
{
  const cs_adjacency_t  *c2f = cs_cdo_connect->c2f;
  const cs_real_t  *w_c2f = cs_cdo_quant->pvol_fc;
  const cs_real_t  *i_surf = cs_cdo_quant->i_face_surf;
  const cs_real_t  *b_surf = cs_cdo_quant->b_face_surf;
  const cs_lnum_t  size = c2f->idx[cs_cdo_quant->n_cells];
  const cs_lnum_t  n_i_faces = cs_cdo_quant->n_i_faces;

  _sanity_checks(__func__, c2f, w_c2f);

  /*
   * The algorithm used is l3superblock60, based on the article:
   * "Reducing Floating Point Error in Dot Product Using the Superblock Family
   * of Algorithms" by Anthony M. Castaldo, R. Clint Whaley, and Anthony
   * T. Chronopoulos, SIAM J. SCI. COMPUT., Vol. 31, No. 2, pp. 1156--1174
   * 2008 Society for Industrial and Applied Mathematics
   */

  double  l2norm = 0;

# pragma omp parallel reduction(+:l2norm) if (size > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(size, &s_id, &e_id);

    const cs_lnum_t  n = e_id - s_id;
    const cs_lnum_t  *_ids = c2f->ids + s_id;
    const cs_real_t  *_w = w_c2f + s_id;
    const cs_lnum_t  block_size = CS_SBLOCK_BLOCK_SIZE;
    const cs_lnum_t  n_blocks = (n + block_size - 1) / block_size;
    const cs_lnum_t  n_sblocks = (n_blocks > 3) ? sqrt(n_blocks) : 1;
    const cs_lnum_t  blocks_in_sblocks =
      (n + block_size*n_sblocks - 1) / (block_size*n_sblocks);

    cs_lnum_t  shift = 0;

    for (cs_lnum_t s = 0; s < n_sblocks; s++) { /* Loop on slices */

      double  s_l2norm = 0.0;

      for (cs_lnum_t b_id = 0; b_id < blocks_in_sblocks; b_id++) {

        const cs_lnum_t  start_id = shift;
        shift += block_size;
        if (shift > n)
          shift = n, b_id = blocks_in_sblocks;
        const cs_lnum_t  end_id = shift;

        double  _l2norm = 0.0;
        for (cs_lnum_t j = start_id; j < end_id; j++) {

          const cs_lnum_t  id = _ids[j];
          const cs_real_t  surf =
            (id < n_i_faces) ? i_surf[id] : b_surf[id - n_i_faces];
          const cs_real_t  v = array[id]/surf;

          _l2norm += _w[j] * v*v;

        } /* Loop on block_size */

        s_l2norm += _l2norm;

      } /* Loop on blocks */

      l2norm += s_l2norm;

    } /* Loop on super-blocks */

  } /* OpenMP block */

  /* Parallel treatment */

  cs_parall_sum(1, CS_REAL_TYPE, &l2norm);

  return (cs_real_t)l2norm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the norm ||b - a||**2
 *         Case of a scalar-valued array defined as a flux at primal
 *         faces. Thus, the weigth is the pyramid of apex the cell center and
 *         of basis the face. Each face quantity is normalized by the face
 *         surface. The computed quantities are synchronized in parallel.
 *
 * \param[in]  a   first array
 * \param[in]  b   second array
 *
 * \return the value of ||b - a||**2
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdo_blas_square_norm_pfsf_diff(const cs_real_t        *a,
                                  const cs_real_t        *b)
{
  const cs_adjacency_t  *c2f = cs_cdo_connect->c2f;
  const cs_real_t  *w_c2f = cs_cdo_quant->pvol_fc;
  const cs_real_t  *i_surf = cs_cdo_quant->i_face_surf;
  const cs_real_t  *b_surf = cs_cdo_quant->b_face_surf;
  const cs_lnum_t  size = c2f->idx[cs_cdo_quant->n_cells];
  const cs_lnum_t  n_i_faces = cs_cdo_quant->n_i_faces;

  _sanity_checks(__func__, c2f, w_c2f);

  /*
   * The algorithm used is l3superblock60, based on the article:
   * "Reducing Floating Point Error in Dot Product Using the Superblock Family
   * of Algorithms" by Anthony M. Castaldo, R. Clint Whaley, and Anthony
   * T. Chronopoulos, SIAM J. SCI. COMPUT., Vol. 31, No. 2, pp. 1156--1174
   * 2008 Society for Industrial and Applied Mathematics
   */

  double  l2norm = 0;

# pragma omp parallel reduction(+:l2norm) if (size > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(size, &s_id, &e_id);

    const cs_lnum_t  n = e_id - s_id;
    const cs_lnum_t  *_ids = c2f->ids + s_id;
    const cs_real_t  *_w = w_c2f + s_id;
    const cs_lnum_t  block_size = CS_SBLOCK_BLOCK_SIZE;
    const cs_lnum_t  n_blocks = (n + block_size - 1) / block_size;
    const cs_lnum_t  n_sblocks = (n_blocks > 3) ? sqrt(n_blocks) : 1;
    const cs_lnum_t  blocks_in_sblocks =
      (n + block_size*n_sblocks - 1) / (block_size*n_sblocks);

    cs_lnum_t  shift = 0;

    for (cs_lnum_t s = 0; s < n_sblocks; s++) { /* Loop on slices */

      double  s_l2norm = 0.0;

      for (cs_lnum_t b_id = 0; b_id < blocks_in_sblocks; b_id++) {

        const cs_lnum_t  start_id = shift;
        shift += block_size;
        if (shift > n)
          shift = n, b_id = blocks_in_sblocks;
        const cs_lnum_t  end_id = shift;

        double  _l2norm = 0.0;
        for (cs_lnum_t j = start_id; j < end_id; j++) {

          const cs_lnum_t  id = _ids[j];
          const cs_real_t  surf =
            (id < n_i_faces) ? i_surf[id] : b_surf[id - n_i_faces];
          const cs_real_t  v = (b[id] - a[id])/surf;

          _l2norm += _w[j] * v*v;

        } /* Loop on block_size */

        s_l2norm += _l2norm;

      } /* Loop on blocks */

      l2norm += s_l2norm;

    } /* Loop on super-blocks */

  } /* OpenMP block */

  /* Parallel treatment */

  cs_parall_sum(1, CS_REAL_TYPE, &l2norm);

  return (cs_real_t)l2norm;
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
