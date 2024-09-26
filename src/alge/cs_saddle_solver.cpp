/*============================================================================
 * Solvers for saddle-point systems arising from CDO discretizations
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_OPENMP)
#include <omp.h>
#endif

/*----------------------------------------------------------------------------
 *  BFT headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_array.h"
#include "cs_blas.h"
#include "cs_cdo_solve.h"
#include "cs_log.h"
#include "cs_parameters.h"
#include "cs_saddle_system.h"


#if 0  /* Set to 1 if systems have to be exported into a binary file */
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_saddle_solver.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_saddle_solver.c
 *
 * \brief Set of functions to solve saddle-point problems
 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

/* Cache line multiple, in cs_real_t units */

#define CS_CL  (CS_CL_SIZE/8)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Generic pointer of function definition for the application of a
 *        block preconditioner.
 *        M z = r (M is stored inside the context with several specificities)
 *
 * \param[in]      solver  pointer to a saddle-point solver structure
 * \param[in]      ctx     context pointer for block preconditioner
 * \param[in]      r       rhs of the preconditioning system
 * \param[in, out] z       array to compute
 * \param[in, out] pc_wsp  work space related to the preconditioner
 *
 * \return the number of iteration performed for this step of preconditioning
 */
/*----------------------------------------------------------------------------*/

typedef int(cs_saddle_solver_pc_apply_t)(
  cs_saddle_solver_t                   *solver,
  cs_saddle_solver_context_block_pcd_t *ctx,
  cs_real_t                            *r,
  cs_real_t                            *z,
  cs_real_t                            *pc_wsp);

/*============================================================================
 * Private variables
 *============================================================================*/

/*============================================================================
 * Local private variables
 *============================================================================*/

static int  cs_saddle_solver_n_systems = 0;
static cs_saddle_solver_t  **cs_saddle_solver_systems = nullptr;

/*============================================================================
 * Static inline private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the position id in the array used to store a upper-right
 *        triangular matrix
 *
 * \param[in] m  size of the square matrix
 * \param[in] i  row id
 * \param[in] j  column id
 *
 * \return the position id in the array
 */
/*----------------------------------------------------------------------------*/

static inline int
_get_id(int  m,
        int  i,
        int  j)
{
  return j + i*m - (i*(i + 1))/2;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the scalar multiplication of a vector split into the x1 and
 *        x2 part
 *
 * \param[in]      solver  saddle-point solver structure
 * \param[in]      scalar  scaling coefficient to apply
 * \param[in, out] x       vector to update
 */
/*----------------------------------------------------------------------------*/

static inline void
_scalar_scaling(cs_saddle_solver_t  *solver,
                const cs_real_t      scalar,
                cs_real_t           *x)
{
  assert(x != nullptr);
  cs_saddle_solver_context_block_pcd_t *ctx =
    static_cast<cs_saddle_solver_context_block_pcd_t *>(solver->context);
  cs_real_t  *x1 = x, *x2 = x + ctx->b11_max_size;

  cs_array_real_scale(solver->n1_scatter_dofs, 1, nullptr, scalar, x1);
  cs_array_real_scale(solver->n2_scatter_dofs, 1, nullptr, scalar, x2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute array index bounds for a local thread.
 *
 * When called inside an OpenMP parallel section, this will return the
 * start an past-the-end indexes for the array range assigned to that thread.
 * In other cases, the start index is 1, and the past-the-end index is n;
 *
 * \param[in]  n      size of an array
 * \param[out] s_id   start index for the current thread
 * \param[out] e_id   past-the-end index for the current thread
 */
/*----------------------------------------------------------------------------*/

static inline void
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
 * \brief Dot product between two arrays with DoFs on interfaces.
 *        One assumes that input arrays are in a "scattered" distribution One
 *        switches to a "gather" view for the dot product and then moves back
 *        to a "scatter" view.  The global reduction reduction is performed.
 *
 * \param[in]      rset   pointer to a range_set structure (synchro. op.)
 * \param[in]      size   size of arrays
 * \param[in, out] x      first array
 * \param[in, out] y      second array
 *
 * \return the value of the canonical dot product between x and y
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_scatter_global_dotprod(const cs_range_set_t  *rset,
                        cs_lnum_t              size,
                        cs_real_t              x[],
                        cs_real_t              y[])
{
  CS_NO_WARN_IF_UNUSED(size); /* Avoid a compiler warning */
  assert(size == rset->n_elts[1]);

  /* x and y are scattered arrays. One assumes that values are synchronized
     across ranks (for instance by using a cs_interface_set_sum()) */

  cs_range_set_gather(rset,
                      CS_REAL_TYPE, /* type */
                      1,            /* stride (treated as scalar up to now) */
                      x,            /* in: scatter view */
                      x);           /* out: gather view */

  cs_range_set_gather(rset,
                      CS_REAL_TYPE, /* type */
                      1,            /* stride (treated as scalar up to now) */
                      y,            /* in: scatter view */
                      y);           /* out: gather view */

  cs_real_t  result = cs_gdot(rset->n_elts[0], x, y);

  cs_range_set_scatter(rset,
                       CS_REAL_TYPE,
                       1,
                       x,            /* in: gather view */
                       x);           /* out: scatter view */
  cs_range_set_scatter(rset,
                       CS_REAL_TYPE,
                       1,
                       y,            /* in: gather view */
                       y)            /* out: scatter view */ ;

  return result;
}

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the array X such that X += mu * Y
 *        where the X and Y vector are split into two parts
 *
 * \param[in]      solver  saddle-point solver structure
 * \param[in]      mu      multiplicative coefficient to apply
 * \param[in]      y       constant array
 * \param[in, out] x       resulting array to update
 */
/*----------------------------------------------------------------------------*/

static void
_add_scaled_vector(cs_saddle_solver_t  *solver,
                   const cs_real_t      mu,
                   const cs_real_t     *y,
                   cs_real_t           *x)
{
  assert(x != nullptr && y != nullptr);
  cs_saddle_solver_context_block_pcd_t *ctx =
    static_cast<cs_saddle_solver_context_block_pcd_t *>(solver->context);

  cs_real_t  *x1 = x, *x2 = x + ctx->b11_max_size;
  const cs_real_t  *y1 = y, *y2 = y + ctx->b11_max_size;

# pragma omp parallel for if (solver->n1_scatter_dofs > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < solver->n1_scatter_dofs; i++)
    x1[i] += mu*y1[i];

# pragma omp parallel for if (solver->n2_scatter_dofs > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < solver->n2_scatter_dofs; i++)
    x2[i] += mu*y2[i];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the canonical dot product between the vectors x and y when a
 *        block preconditioner is used. The synchronization is performed during
 *        the process.
 *
 * \param[in] solver  saddle-point solver structure
 * \param[in] x       first vector in a scatter state
 * \param[in] y       second vector in a scatter state
 *
 * \return the value of the canonical dot product (x,y)
 */
/*----------------------------------------------------------------------------*/

static double
_block_pcd_dot_product(cs_saddle_solver_t  *solver,
                       cs_real_t           *x,
                       cs_real_t           *y)
{
  double  dp_value = 0.;

  if (x == nullptr || y== nullptr)
    return dp_value;

  cs_saddle_solver_context_block_pcd_t *ctx =
    static_cast<cs_saddle_solver_context_block_pcd_t *>(solver->context);
  cs_real_t  *x1 = x, *x2 = x + ctx->b11_max_size;
  cs_real_t  *y1 = y, *y2 = y + ctx->b11_max_size;

  /* First part x1 and y1 whose DoFs may be shared among processes */

  const cs_range_set_t  *rset = ctx->b11_range_set;

  if (rset != nullptr) { /* Switch to a gather view to avoid summing an element
                         twice */

    cs_range_set_gather(rset,
                        CS_REAL_TYPE,/* type */
                        1,           /* stride (viewed as scalar up to now) */
                        x1,          /* in: size = n_scatter_elts */
                        x1);         /* out: size = n_gather_elts */

    cs_range_set_gather(rset,
                        CS_REAL_TYPE,/* type */
                        1,           /* stride (viewed as scalar up to now) */
                        y1,          /* in: size = n_scatter_elts */
                        y1);         /* out: size = n_gather_elts */

    dp_value = cs_dot(rset->n_elts[0], x1, y1);

    /* Move back to a scatter view */

    cs_range_set_scatter(rset,
                         CS_REAL_TYPE,
                         1,
                         x1,
                         x1);

    cs_range_set_scatter(rset,
                         CS_REAL_TYPE,
                         1,
                         y1,
                         y1);

  }
  else {

    dp_value = cs_dot(solver->n1_scatter_dofs, x1, y1);

  }

  /* Second part: x2 and y2, DoFs which are not shared */

  dp_value += cs_dot(solver->n2_scatter_dofs, x2, y2);

  cs_parall_sum(1, CS_DOUBLE, &dp_value);

  return dp_value;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the norm of a vector in a scatter view which is split into
 *        the x1 and x2 part.
 *        The synchronization is performed during the process.
 *
 * \param[in] solver  pointer to a saddle_point solver structure
 * \param[in] x       vector used for the norm computation
 *
 * \return the value of the euclidean norm of x
 */
/*----------------------------------------------------------------------------*/

static double
_norm(cs_saddle_solver_t  *solver,
      cs_real_t           *x)
{
  double  n_square_value = 0;

  if (x == nullptr)
    return n_square_value;
  assert(solver != nullptr);

  cs_saddle_solver_context_block_pcd_t *ctx =
    static_cast<cs_saddle_solver_context_block_pcd_t *>(solver->context);
  cs_real_t  *x1 = x, *x2 = x + ctx->b11_max_size;

  const cs_range_set_t  *rset = ctx->b11_range_set;

  /* Norm for the x1 DoFs (those shared among processes) */

  double  _nx1_square = 0;

  if (rset != nullptr) { /* Switch to a gather view to avoid summing an element
                         twice */

    cs_range_set_gather(rset,
                        CS_REAL_TYPE,/* type */
                        1,           /* stride (treated as scalar up to now) */
                        x1,          /* in: size = n_sles_scatter_elts */
                        x1);         /* out: size = n_sles_gather_elts */

    /* n_elts[0] corresponds to the number of element in the gather view */

    _nx1_square = cs_dot_xx(rset->n_elts[0], x1);

    cs_range_set_scatter(rset,
                         CS_REAL_TYPE,
                         1,
                         x1,
                         x1);

  }
  else
    _nx1_square = cs_dot_xx(solver->n1_scatter_dofs, x1);

  /* Norm for the x2 DoFs (not shared so that there is no need to
     synchronize) */

  double  _nx2_square = cs_dot_xx(solver->n2_scatter_dofs, x2);

  n_square_value = _nx1_square + _nx2_square;
  cs_parall_sum(1, CS_DOUBLE, &n_square_value);
  assert(n_square_value > -DBL_MIN);

  return sqrt(n_square_value);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Join the first set of elements (velocity for instance) with the
 *        second set of elements (pressure for instance) into the same
 *        array. Do the same thing for the output rhs named b. Moreover,
 *        de-interlace the components of the first set of elements. A stride
 *        equal to 3 is assumed.
 *
 * \param[in]      n1_elts  number of elements in x1
 * \param[in]      x1       first array
 * \param[in]      n2_elts  number of elements in x2
 * \param[in]      x2       second array
 * \param[in]      rhs      system rhs (interlaced)
 * \param[in, out] u        solution array to define
 * \param[in, out] b        rhs array to define
 */
/*----------------------------------------------------------------------------*/

static void
_join_x1_vector_x2_deinterlaced(cs_lnum_t         n1_elts,
                                const cs_real_t  *x1,
                                cs_lnum_t         n2_elts,
                                const cs_real_t  *x2,
                                const cs_real_t  *rhs,
                                cs_real_t        *u,
                                cs_real_t        *b)
{
  cs_real_t  *ux = u,             *bx = b;
  cs_real_t  *uy = u + n1_elts,   *by = b + n1_elts;
  cs_real_t  *uz = u + 2*n1_elts, *bz = b + 2*n1_elts;

  /* Treatment of the first set of elements. */

# pragma omp parallel for if (CS_THR_MIN > n1_elts)
  for (cs_lnum_t i1 = 0; i1 < n1_elts; i1++) {

    const cs_real_t  *_x1 = x1 + 3*i1;
    const cs_real_t  *_rhs = rhs + 3*i1;

    ux[i1] = _x1[0], bx[i1] = _rhs[0];
    uy[i1] = _x1[1], by[i1] = _rhs[1];
    uz[i1] = _x1[2], bz[i1] = _rhs[2];

  }

  /* Treatment of the second set of elements */

  const cs_lnum_t  x1_shift = 3*n1_elts;

  cs_array_real_copy(n2_elts, x2, u + x1_shift);
  cs_array_real_copy(n2_elts, rhs + x1_shift, b + x1_shift);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Join the first set of elements (flux for instance) with the
 *        second set of elements (potential for instance) into the same
 *        array. Do the same thing for the output rhs named b.
 *
 * \param[in]      n1_elts  number of elements in x1
 * \param[in]      x1       first array
 * \param[in]      n2_elts  number of elements in x2
 * \param[in]      x2       second array
 * \param[in]      rhs      system rhs (interlaced)
 * \param[in, out] u        solution array to define
 * \param[in, out] b        rhs array to define
 */
/*----------------------------------------------------------------------------*/

static void
_join_x1_scalar_x2_deinterlaced(cs_lnum_t         n1_elts,
                                const cs_real_t  *x1,
                                cs_lnum_t         n2_elts,
                                const cs_real_t  *x2,
                                const cs_real_t  *rhs,
                                cs_real_t        *u,
                                cs_real_t        *b)
{
  /* Treatment of the first set of elements. */

  cs_array_real_copy(n1_elts, x1, u);
  cs_array_real_copy(n1_elts, rhs, b);

  /* Treatment of the second set of elements */

  cs_array_real_copy(n2_elts, x2, u + n1_elts);
  cs_array_real_copy(n2_elts, rhs + n1_elts, b + n1_elts);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform a matrix-vector multiplication for the (1,1) block when the
 *        input vector is in a scatter state.

 *        Thus, one performs a scatter --> gather (before the multiplication)
 *        and a gather --> scatter operation after the multiplication. One
 *        assumes that m11x1 is allocated to the right size. No check is done.
 *
 *        The m11 matrix is assumed to have a stride equal to 1 (db_size[3]=1)
 *
 * \param[in]      rset      pointer to a cs_range_set_t structure
 * \param[in]      m11       matrix
 * \param[in, out] vec       vector
 * \param[in, out] matvec    resulting vector for the matrix-vector product
 */
/*----------------------------------------------------------------------------*/

static void
_m11_vector_multiply_allocated(const cs_range_set_t  *rset,
                               const cs_matrix_t     *m11,
                               cs_real_t             *x1,
                               cs_real_t             *m11x1)
{
  if (m11 == nullptr || x1 == nullptr)
    return;

  /* Remark:
   * n_rows = n_gather_elts <= n_scatter_elts = n_dofs (mesh view) <= n_cols
   */

  /* scatter view to gather view for the input vector */

  cs_range_set_gather(rset,
                      CS_REAL_TYPE, 1, /* type and stride */
                      x1,             /* in:  size=n1_scatter_elts */
                      x1);            /* out: size=n1_gather_elts */

  cs_range_set_gather(rset,
                      CS_REAL_TYPE, 1, /* type and stride */
                      m11x1,          /* in:  size=n1_scatter_elts */
                      m11x1);         /* out: size=n1_gather_elts */

  cs_matrix_vector_multiply(m11, x1, m11x1);

  /* gather view to scatter view (i.e. algebraic to mesh view) */

  cs_range_set_scatter(rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       x1, x1);

  cs_range_set_scatter(rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       m11x1, m11x1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the residual of the saddle-point system
 *        res1 = rhs1 - M11.x1 - M12.x2 (residual for the first block)
 *        res2 = rhs2 - M21.x1          (residual for the second block)
 *
 *        The matrix m11 is represented with 1 block.
 *        The stride is equal to 1 or 3 for the operator m21_unassembled
 *
 * \param[in]      solver  saddle-point solver structure
 * \param[in, out] x1      array for the first part
 * \param[in, out] x2      array for the second part
 * \param[in, out] res     resulting residual vector
 */
/*----------------------------------------------------------------------------*/

static void
_compute_residual(cs_saddle_solver_t *solver,
                  cs_real_t          *x1,
                  cs_real_t          *x2,
                  cs_real_t          *res)
{
  assert(res != nullptr && x1 != nullptr && x2 != nullptr && solver != nullptr);
  assert(solver->n1_dofs_by_elt == 1 || solver->n1_dofs_by_elt == 3);
  assert(solver->n2_dofs_by_elt == 1);

  const cs_real_t *rhs1 = solver->system_helper->rhs_array[0];
  const cs_real_t *rhs2 = solver->system_helper->rhs_array[1];
  const cs_saddle_solver_context_block_pcd_t *ctx =
    static_cast<const cs_saddle_solver_context_block_pcd_t *>(solver->context);
  const cs_range_set_t                       *rset = ctx->b11_range_set;

  cs_real_t *res1 = res, *res2 = res + ctx->b11_max_size;

  cs_array_real_fill_zero(ctx->b11_max_size, res1);

  /* Two parts:
   * a) res1 = rhs1 - M11.x1 - M12.x2
   * b) res2 = rhs2 - M21.x1
   */

  cs_real_t *m12x2 = nullptr;
  BFT_MALLOC(m12x2, solver->n1_scatter_dofs, cs_real_t);
  cs_array_real_fill_zero(solver->n1_scatter_dofs, m12x2);

  const cs_adjacency_t *adj = ctx->m21_adj;

  assert(solver->n2_elts == adj->n_elts);

  if (solver->n1_dofs_by_elt == 3) {
#pragma omp parallel for if (solver->n2_elts > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < solver->n2_elts; i2++) {

      const cs_real_t _x2    = x2[i2];
      cs_real_t       _m21x1 = 0.;
      for (cs_lnum_t j2 = adj->idx[i2]; j2 < adj->idx[i2 + 1]; j2++) {

        const cs_real_t *m21_vals = ctx->m21_val + 3 * j2;
        const cs_lnum_t  shift    = 3 * adj->ids[j2];
        assert(shift < solver->n1_scatter_dofs);

        _m21x1 += cs_math_3_dot_product(m21_vals, x1 + shift);

        cs_real_t *_m12x2 = m12x2 + shift;
#pragma omp critical
        {
          _m12x2[0] += m21_vals[0] * _x2;
          _m12x2[1] += m21_vals[1] * _x2;
          _m12x2[2] += m21_vals[2] * _x2;
        }

      } /* Loop on x1 elements associated to a given x2 element */

      res2[i2] = rhs2[i2] - _m21x1;

    } /* Loop on x2 elements */
  }
  else if (solver->n1_dofs_by_elt == 1) {
#pragma omp parallel for if (solver->n2_elts > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < solver->n2_elts; i2++) {

      const cs_real_t _x2    = x2[i2];
      cs_real_t       _m21x1 = 0.;
      for (cs_lnum_t j2 = adj->idx[i2]; j2 < adj->idx[i2 + 1]; j2++) {

        assert(adj->ids[j2] < solver->n1_scatter_dofs);

        _m21x1 += ctx->m21_val[j2] * x1[adj->ids[j2]];

#pragma omp critical
        {
          m12x2[adj->ids[j2]] += ctx->m21_val[j2] * _x2;
        }

      } /* Loop on x1 elements associated to a given x2 element */

      res2[i2] = rhs2[i2] - _m21x1;

    } /* Loop on x2 elements */
  }
  else {
    bft_error(__FILE__,
              __LINE__,
              0,
              "%s: Stride %d is not supported.\n",
              __func__,
              solver->n1_dofs_by_elt);
  }

  if (rset->ifs != nullptr)
    cs_interface_set_sum(rset->ifs,
                         solver->n1_scatter_dofs,
                         1,
                         false,
                         CS_REAL_TYPE, /* stride, interlaced */
                         m12x2);

  _m11_vector_multiply_allocated(rset, ctx->m11, x1, res1);

#pragma omp parallel for if (solver->n1_scatter_dofs > CS_THR_MIN)
  for (cs_lnum_t i1 = 0; i1 < solver->n1_scatter_dofs; i1++)
    res1[i1] = rhs1[i1] - res1[i1] - m12x2[i1];

  BFT_FREE(m12x2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the matrix-vector operation divided into two parts
 *        matvec1 = M11.vec1 + M12.vec2
 *        matvec2 = M21.vec1
 *
 *        The stride is equal to 1 or 3 for the operator m21_unassembled
 *
 * \param[in]      solver  pointer to a saddl-point solver structure
 * \param[in, out] vec     array to multiply
 * \param[out]     matvec  resulting vector
 */
/*----------------------------------------------------------------------------*/

static void
_matvec_product(cs_saddle_solver_t *solver, cs_real_t *vec, cs_real_t *matvec)
{
  assert(vec != nullptr && matvec != nullptr && solver != nullptr);
  assert(solver->n1_dofs_by_elt == 1 || solver->n1_dofs_by_elt == 3);
  assert(solver->n2_dofs_by_elt == 1);

  cs_saddle_solver_context_block_pcd_t *ctx =
    static_cast<cs_saddle_solver_context_block_pcd_t *>(solver->context);

  cs_real_t *v1 = vec, *v2 = vec + ctx->b11_max_size;
  cs_real_t *mv1 = matvec, *mv2 = matvec + ctx->b11_max_size;

  const cs_range_set_t *rset = ctx->b11_range_set;
  const cs_adjacency_t *adj  = ctx->m21_adj;

  /* The resulting array is composed of two parts:
   * a) mv1 = M11.v1 + M12.v2
   * b) mv2 = M21.v1
   */

  /* 1) M12.v2 and M21.v1 */

  cs_real_t *m12v2 = nullptr;
  BFT_MALLOC(m12v2, solver->n1_scatter_dofs, cs_real_t);
  cs_array_real_fill_zero(solver->n1_scatter_dofs, m12v2);

  assert(solver->n2_scatter_dofs == adj->n_elts);

  if (solver->n1_dofs_by_elt == 3) {
#pragma omp parallel for if (adj->n_elts > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < adj->n_elts; i2++) {

      const cs_real_t _v2    = v2[i2];
      cs_real_t       _m21v1 = 0.;
      for (cs_lnum_t j2 = adj->idx[i2]; j2 < adj->idx[i2 + 1]; j2++) {

        const cs_lnum_t  shift    = 3 * adj->ids[j2];
        const cs_real_t *m21_vals = ctx->m21_val + 3 * j2;

        _m21v1 += cs_math_3_dot_product(m21_vals, v1 + shift);

        cs_real_t *_m12v2 = m12v2 + shift;
#pragma omp critical
        {
          _m12v2[0] += m21_vals[0] * _v2;
          _m12v2[1] += m21_vals[1] * _v2;
          _m12v2[2] += m21_vals[2] * _v2;
        }

      } /* Loop on x1 elements associated to a given x2 element */

      mv2[i2] = _m21v1;

    } /* Loop on x2 elements */
  }
  else if (solver->n1_dofs_by_elt == 1) {
#pragma omp parallel for if (adj->n_elts > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < adj->n_elts; i2++) {

      const cs_real_t _v2    = v2[i2];
      cs_real_t       _m21v1 = 0.;
      for (cs_lnum_t j2 = adj->idx[i2]; j2 < adj->idx[i2 + 1]; j2++) {

        _m21v1 += ctx->m21_val[j2] * v1[adj->ids[j2]];

#pragma omp critical
        {
          m12v2[adj->ids[j2]] += ctx->m21_val[j2] * _v2;
        }

      } /* Loop on x1 elements associated to a given x2 element */

      mv2[i2] = _m21v1;

    } /* Loop on x2 elements */
  }
  else {
    bft_error(__FILE__,
              __LINE__,
              0,
              "%s: Stride %d is not supported.\n",
              __func__,
              solver->n1_dofs_by_elt);
  }

  if (rset->ifs != nullptr)
    cs_interface_set_sum(rset->ifs,
                         solver->n1_scatter_dofs,
                         1,
                         false,
                         CS_REAL_TYPE, /* stride, interlaced */
                         m12v2);

  _m11_vector_multiply_allocated(rset, ctx->m11, v1, mv1);

#pragma omp parallel for if (solver->n1_scatter_dofs > CS_THR_MIN)
  for (cs_lnum_t i1 = 0; i1 < solver->n1_scatter_dofs; i1++)
    mv1[i1] += m12v2[i1];

  BFT_FREE(m12v2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve the (1,1)-block (version with one matrix)
 *
 * \param[in]      solver     pointer to a saddle-point solver structure
 * \param[in]      ctx        context pointer for block preconditioner
 * \param[in]      r          rhs of the M11 sub-system (size = b11_max_size)
 * \param[in]      scatter_r  true=perform scatter op. on r after solving
 * \param[in, out] z          array to compute (size = b11_max_size)
 * \param[in]      scatter_z  true=perform scatter op. on z after solving
 *
 * \return the number of iterations performed for this step
 */
/*----------------------------------------------------------------------------*/

static int
_solve_m11_approximation(cs_saddle_solver_t                   *solver,
                         cs_saddle_solver_context_block_pcd_t *ctx,
                         cs_real_t                            *r,
                         bool                                  scatter_r,
                         cs_real_t                            *z,
                         bool                                  scatter_z)
{
  const cs_range_set_t  *rset = ctx->b11_range_set;

  cs_matrix_t  *m11 = ctx->m11;

  /* Handle parallelism: scatter --> gather transformation
   * No gather op. for z since we initialize with zero */

  cs_range_set_gather(rset,
                      CS_REAL_TYPE, 1, /* stride as if scalar-valued */
                      r,               /* in:  size = n1_scatter_dofs */
                      r);              /* out: size = rset->n_elts[0] */

  /* Compute the norm of the rhs to normalize the linear system */

  cs_lnum_t  n_x1_elts =
    (rset != nullptr) ? rset->n_elts[0] : solver->n1_scatter_dofs;

  double  r_norm = cs_dot_xx(n_x1_elts, r);
  cs_parall_sum(1, CS_DOUBLE, &r_norm);
  r_norm = sqrt(fabs(r_norm));

  /* Solve the linear solver */

  cs_array_real_fill_zero(n_x1_elts, z);

  cs_solving_info_t m11_info = { .n_it       = 0,
                                 .rhs_norm   = r_norm,
                                 .res_norm   = DBL_MAX,
                                 .derive     = DBL_MAX,
                                 .l2residual = DBL_MAX };

  const cs_param_saddle_t  *saddlep = solver->param;
  const cs_param_sles_t  *m11_slesp = saddlep->block11_sles_param;
  assert(m11_slesp != nullptr);

  /* A gather view is used inside the following step */

  cs_sles_convergence_state_t  code = cs_sles_solve(solver->main_sles,
                                                    m11,
                                                    m11_slesp->cvg_param.rtol,
                                                    m11_info.rhs_norm,
                                                    &(m11_info.n_it),
                                                    &(m11_info.res_norm),
                                                    r,
                                                    z,
                                                    0,      /* aux. size */
                                                    nullptr);  /* aux. buffers */

  /* Output information about the convergence of the resolution */

  if (m11_slesp->verbosity > 1 && cs_log_default_is_active()) {

    cs_log_printf(CS_LOG_DEFAULT, "  <%s/%20s> "
                  "cvg_code:%-d | n_iter:%3d | residual:% -8.4e",
                  __func__, m11_slesp->name, code,
                  m11_info.n_it, m11_info.res_norm);

    if (m11_slesp->verbosity > 2)
      cs_log_printf(CS_LOG_DEFAULT, " | rhs.norm % -8.4e\n", m11_info.rhs_norm);
    cs_log_printf(CS_LOG_DEFAULT, "\n");

  }

  /* Move back: gather --> scatter view */

  if (scatter_z)
    cs_range_set_scatter(rset, CS_REAL_TYPE, 1, /* type and stride */
                         z, z);
  if (scatter_r)
    cs_range_set_scatter(rset, CS_REAL_TYPE, 1, /* type and stride */
                         r, r);

  return m11_info.n_it;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve the system related to the Schur complement approximation
 *        S.z_schur = r_schur
 *
 *        This is common to the following strategies:
 *         - Block preconditioning with a Krylov solver
 *         - Uzawa-CG algorithm
 *
 * \param[in, out] solver       pointer to a saddle-point solver structure
 * \param[in, out] schur_sles   SLES related to the Schur complement (if needed)
 * \param[in]      schur_mat    matrix for the Schur complement (if needed)
 * \param[in       inv_m22      diagonal inverseof the M22 mass matrix
 * \param[in]      scaling      scaling to apply in case of hybrid approx.
 * \param[in]      r_schur      rhs of the Schur system (size = n2_scatter_dofs)
 * \param[in, out] z_schur      array to compute (size = n2_scatter_dofs)
 *
 * \return the number of iterations performed for this step of preconditioning
 */
/*----------------------------------------------------------------------------*/

static int
_solve_schur_approximation(cs_saddle_solver_t  *solver,
                           cs_sles_t           *schur_sles,
                           cs_matrix_t         *schur_mat,
                           const cs_real_t     *inv_m22,
                           cs_real_t            scaling,
                           cs_real_t           *r_schur,
                           cs_real_t           *z_schur)
{
  assert(r_schur != nullptr && z_schur != nullptr);

  int  n_iter = 0;

  const cs_lnum_t  n2_elts = solver->n2_scatter_dofs;
  const cs_param_saddle_t  *saddlep = solver->param;
  const cs_param_saddle_schur_approx_t  schur_type = saddlep->schur_approx;

  switch (schur_type) {

  case CS_PARAM_SADDLE_SCHUR_NONE:
  case CS_PARAM_SADDLE_SCHUR_IDENTITY: /* Identity => z_schur <-- r_schur */
    cs_array_real_copy(n2_elts, r_schur, z_schur);
    break;

  case CS_PARAM_SADDLE_SCHUR_MASS_SCALED:
    assert(inv_m22 != nullptr);

#   pragma omp parallel for if (n2_elts > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < n2_elts; i2++)
      z_schur[i2] = inv_m22[i2] * r_schur[i2]; /* TODO: multiply or divide ?
                                                  Uza:divide
                                                  Block pcd:multiply */
    break;

  default:
    {
      assert(schur_sles != nullptr);
      assert(schur_mat != nullptr);
      assert(saddlep->schur_sles_param != nullptr);

      /* Norm for the x2 DoFs (not shared so that there is no need to
         synchronize) */

      double  r_norm = cs_dot_xx(n2_elts, r_schur);
      cs_parall_sum(1, CS_DOUBLE, &r_norm);
      r_norm = sqrt(fabs(r_norm));

      cs_array_real_fill_zero(n2_elts, z_schur);

      n_iter += cs_cdo_solve_scalar_cell_system(n2_elts,
                                                saddlep->schur_sles_param,
                                                schur_mat,
                                                r_norm,
                                                schur_sles,
                                                z_schur,
                                                r_schur);

      if (schur_type == CS_PARAM_SADDLE_SCHUR_MASS_SCALED_DIAG_INVERSE ||
          schur_type == CS_PARAM_SADDLE_SCHUR_MASS_SCALED_LUMPED_INVERSE) {

        assert(inv_m22 != nullptr);

#       pragma omp parallel for if (n2_elts > CS_THR_MIN)
        for (cs_lnum_t i2 = 0; i2 < n2_elts; i2++)
          z_schur[i2] = scaling*z_schur[i2] + inv_m22[i2]*r_schur[i2];

      } /* Hybrid approximations */

    }
    break;

  } /* Switch on the type of Schur complement approximation */

  return n_iter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply a block diagonal preconditioning: Compute z s.t. P_d z = r
 *        This prototype complies with the cs_saddle_solver_pc_apply_t function
 *        pointer
 *
 * \param[in]      solver  pointer to a saddle-point solver structure
 * \param[in]      ctx     context pointer for block preconditioner
 * \param[in]      r       rhs of the preconditioning system
 * \param[in, out] z       array to compute
 * \param[in, out] pc_wsp  work space related to the preconditioner
 *
 * \return the number of iterations performed for this step of preconditioning
 */
/*----------------------------------------------------------------------------*/

static int
_diag_schur_pc_apply(cs_saddle_solver_t                   *solver,
                     cs_saddle_solver_context_block_pcd_t *ctx,
                     cs_real_t                            *r,
                     cs_real_t                            *z,
                     cs_real_t                            *pc_wsp)
{
  CS_UNUSED(pc_wsp);

  if (z == nullptr)
    return 0;

  assert(solver != nullptr && ctx != nullptr && r != nullptr);

  /* 1. Solve M11.z1 = r1
     ==================== */

  int  n_iter = _solve_m11_approximation(solver, ctx,
                                         r, true,  /* = scatter r */
                                         z, true); /* = scatter z */

  /* 2. Solve S.z2 = r2 (the M22 block is a Schur approx.)
     ===================================================== */

  cs_real_t  *z2 = z + ctx->b11_max_size;
  cs_real_t  *r2 = r + ctx->b11_max_size;

  n_iter += _solve_schur_approximation(solver,
                                       ctx->schur_sles,
                                       ctx->schur_matrix,
                                       ctx->m22_mass_diag,  /* inv_m22 */
                                       ctx->schur_scaling,
                                       r2,                  /* r_schur */
                                       z2);                 /* z_schur */

  return n_iter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply a lower triangular block preconditioning with a potential Schur
 *        complement approximation. Compute z such that P_l z = r
 *
 *        This prototype complies with the cs_saddle_solver_pc_apply_t function
 *        pointer
 *
 * \param[in]      solver  pointer to a saddle-point solver structure
 * \param[in]      ctx     context pointer for block preconditioner
 * \param[in]      r       rhs of the preconditioning system
 * \param[in, out] z       array to compute
 * \param[in, out] pc_wsp  work space related to the preconditioner
 *
 * \return the number of iterations performed for this step of preconditioning
 */
/*----------------------------------------------------------------------------*/

static int
_lower_schur_pc_apply(cs_saddle_solver_t                   *solver,
                      cs_saddle_solver_context_block_pcd_t *ctx,
                      cs_real_t                            *r,
                      cs_real_t                            *z,
                      cs_real_t                            *pc_wsp)
{
  if (z == nullptr)
    return 0;

  assert(solver != nullptr && ctx != nullptr && r != nullptr);
  assert(pc_wsp != nullptr);

  /* 1. Solve m11 z1 = r1
     ==================== */

  int  n_iter = _solve_m11_approximation(solver, ctx,
                                         r, true,   /* = scatter r */
                                         z, true);  /* = scatter z */

  /* 2. Build r2_tilda = r2 - m21.z1
     ============================= */

  const cs_real_t  *r2 = r + ctx->b11_max_size;
  cs_real_t  *r2_tilda = pc_wsp;

  ctx->m21_vector_multiply(solver->n2_scatter_dofs,
                           z, ctx->m21_adj, ctx->m21_val,
                           r2_tilda);

# pragma omp parallel for if (solver->n2_scatter_dofs > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < solver->n2_scatter_dofs; i2++)
    r2_tilda[i2] = r2[i2] - r2_tilda[i2];

  /* 3. Solve S z2 = r2_tilda (S -> Schur approximation for the m22 block)
     ===================================================================== */

  cs_real_t  *z2 = z + ctx->b11_max_size;

  n_iter += _solve_schur_approximation(solver,
                                       ctx->schur_sles,
                                       ctx->schur_matrix,
                                       ctx->m22_mass_diag,  /* inv_m22 */
                                       ctx->schur_scaling,
                                       r2_tilda,
                                       z2);

  return n_iter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply an upper triangular block preconditioning.
 *        Compute z such that P_up z = r
 *        This prototype complies with the cs_saddle_solver_pc_apply_t function
 *        pointer
 *
 * \param[in]      solver  pointer to a saddle-point solver structure
 * \param[in]      ctx     context pointer for block preconditioner
 * \param[in]      r       rhs of the preconditioning system
 * \param[in, out] z       array to compute
 * \param[in, out] pc_wsp  work space related to the preconditioner
 *
 * \return the number of iterations performed for this step of preconditioning
 */
/*----------------------------------------------------------------------------*/

static int
_upper_schur_pc_apply(cs_saddle_solver_t                   *solver,
                      cs_saddle_solver_context_block_pcd_t *ctx,
                      cs_real_t                            *r,
                      cs_real_t                            *z,
                      cs_real_t                            *pc_wsp)
{
  if (z == nullptr)
    return 0;

  int  n_iter = 0;

  assert(solver != nullptr && ctx != nullptr && r != nullptr);
  assert(pc_wsp != nullptr);

  cs_real_t  *z2 = z + ctx->b11_max_size;
  cs_real_t  *r2 = r + ctx->b11_max_size;

  /* 1. Solve S z2 = r2 (S -> Schur approximation for the m22 block)
     =============================================================== */

  n_iter += _solve_schur_approximation(solver,
                                       ctx->schur_sles,
                                       ctx->schur_matrix,
                                       ctx->m22_mass_diag,  /* inv_m22 */
                                       ctx->schur_scaling,
                                       r2,                  /* r_schur */
                                       z2);                 /* z_schur */

  /* 2. Build r1_tilda = r1 - Bt.z2
     ============================== */

  const cs_range_set_t  *rset = ctx->b11_range_set;

  cs_real_t  *r1_tilda = pc_wsp;

  cs_array_real_fill_zero(ctx->b11_max_size, r1_tilda);
  ctx->m12_vector_multiply(solver->n2_scatter_dofs,
                           z2, ctx->m21_adj, ctx->m21_val,
                           r1_tilda);

  if (rset->ifs != nullptr)
    cs_interface_set_sum(rset->ifs,
                         solver->n1_scatter_dofs,
                         1, false, CS_REAL_TYPE, /* stride, interlaced */
                         r1_tilda);

# pragma omp parallel for if (solver->n1_scatter_dofs > CS_THR_MIN)
  for (cs_lnum_t i1 = 0; i1 < solver->n1_scatter_dofs; i1++)
    r1_tilda[i1] = r[i1] - r1_tilda[i1];

  /* 3. Solve m11 z1 = r1_tilda
     ========================== */

  /* No need to scatter the r1_tilda array since it is not used anymore */

  n_iter += _solve_m11_approximation(solver, ctx,
                                     r1_tilda, false, /* no scatter */
                                     z, true);        /* scatter z at exit */

  return n_iter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply symmetric Gauss-Seidel block preconditioning
 *        Compute z such that P_sgs z = r
 *        This prototype complies with the cs_saddle_solver_pc_apply_t function
 *        pointer
 *
 * \param[in]      solver  pointer to a saddle-point solver structure
 * \param[in]      ctx     context pointer for block preconditioner
 * \param[in]      sbp     block-preconditioner for the Saddle-point problem
 * \param[in]      r       rhs of the preconditioning system
 * \param[in, out] z       array to compute
 * \param[in, out] pc_wsp  work space related to the preconditioner
 *
 * \return the number of iterations performed for this step of preconditioning
 */
/*----------------------------------------------------------------------------*/

static int
_sgs_schur_pc_apply(cs_saddle_solver_t                   *solver,
                    cs_saddle_solver_context_block_pcd_t *ctx,
                    cs_real_t                            *r,
                    cs_real_t                            *z,
                    cs_real_t                            *pc_wsp)
{
  if (z == nullptr)
    return 0;

  assert(solver != nullptr && ctx != nullptr && r != nullptr);
  assert(pc_wsp != nullptr);

  const cs_lnum_t  n2_elts = solver->n2_scatter_dofs;
  const cs_range_set_t  *rset = ctx->b11_range_set;
  const cs_matrix_t  *m11 = ctx->m11;

  /* 1. Solve m11 z1_hat = r1
     ======================== */

  cs_real_t  *z1_hat = pc_wsp;  /* x1_size */

  int  n_iter = _solve_m11_approximation(solver, ctx,
                                         r, true,       /* = scatter r */
                                         z1_hat, true); /* = scatter z */

  /* 2. Build r2_hat = m21.z1_hat - r2
     ================================= */

  cs_real_t  *r2_hat = pc_wsp + ctx->b11_max_size; /* x2_size */

  ctx->m21_vector_multiply(n2_elts, z1_hat, ctx->m21_adj, ctx->m21_val,
                           r2_hat); /* z1_hat has to be "scatter" */

  const cs_real_t  *r2 = r + ctx->b11_max_size;

# pragma omp parallel for if (n2_elts > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < n2_elts; i2++)
    r2_hat[i2] = r2_hat[i2] - r2[i2];

  /* 3. Solve S z2 = r2_hat (S -> Schur approximation for the m22 block)
     =================================================================== */

  cs_real_t  *z2 = z + ctx->b11_max_size;

  n_iter += _solve_schur_approximation(solver,
                                       ctx->schur_sles,
                                       ctx->schur_matrix,
                                       ctx->m22_mass_diag,  /* inv_m22 */
                                       ctx->schur_scaling,
                                       r2_hat,              /* r_schur */
                                       z2);                 /* z_schur */

  /* 4. Build r1_tilda = m11.z1_hat + m12.z2
     ======================================= */

  cs_real_t  *r1_tilda = r2_hat; /* x1_size */

  /* scatter view to gather view */

  cs_range_set_gather(rset, CS_REAL_TYPE, 1, /* stride */
                      z1_hat,   /* in:  size=n1_scatter_elts */
                      z1_hat);  /* out: size=n1_gather_elts */

  cs_matrix_vector_multiply(m11, z1_hat, r1_tilda);

  /* gather to scatter view (i.e. algebraic to mesh view) */

  cs_range_set_scatter(rset, CS_REAL_TYPE, 1, /* type and stride */
                       r1_tilda,
                       r1_tilda);

  /* Update += of r1_tilda inside the following function */

  ctx->m12_vector_multiply(n2_elts, z2, ctx->m21_adj, ctx->m21_val,
                           r1_tilda);

  /* 5. Solve M11.z1 = r1_tilda
     ========================== */

  n_iter += _solve_m11_approximation(solver, ctx,
                                     r1_tilda, false, /* = no scatter */
                                     z, false);       /* = no scatter z */

  /* Last update z1 = -z1 + 2 z1_hat (still in gather wiew for both arrays) */

# pragma omp parallel for if (rset->n_elts[0] > CS_THR_MIN)
  for (cs_lnum_t i1 = 0; i1 < rset->n_elts[0]; i1++)
    z[i1] = 2*z1_hat[i1] - z[i1];

  /* Move back: gather --> scatter view */

  cs_range_set_scatter(rset, CS_REAL_TYPE, 1, /* type and stride */
                       z,
                       z);

  return n_iter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply one iteration of Uzawa algo. with a potential Schur approx.
 *        Compute z such that P_uza z = r
 *        This prototype complies with the cs_saddle_solver_pc_apply_t function
 *        pointer
 *
 * \param[in]      solver  pointer to a saddle-point solver structure
 * \param[in]      ctx     context pointer for block preconditioner
 * \param[in]      r       rhs of the preconditioning system
 * \param[in, out] z       array to compute
 * \param[in, out] pc_wsp  work space related to the preconditioner
 *
 * \return the number of iterations performed for this step of preconditioning
 */
/*----------------------------------------------------------------------------*/

static int
_uza_schur_pc_apply(cs_saddle_solver_t                   *solver,
                    cs_saddle_solver_context_block_pcd_t *ctx,
                    cs_real_t                            *r,
                    cs_real_t                            *z,
                    cs_real_t                            *pc_wsp)
{
  if (z == nullptr)
    return 0;

  assert(solver != nullptr && ctx != nullptr);
  assert(r != nullptr);
  assert(pc_wsp != nullptr);

  const cs_lnum_t  n1_elts = solver->n1_scatter_dofs;
  const cs_lnum_t  n2_elts = solver->n2_scatter_dofs;
  const cs_range_set_t  *rset = ctx->b11_range_set;

  /* 1. Solve m11 z1 = r1
     ==================== */

  int  n_iter = _solve_m11_approximation(solver, ctx,
                                         r, true,   /* = scatter r */
                                         z, true);  /* = scatter z */

  /* 2. Build r2_hat = r2 - m21.z1
     ============================= */

  cs_real_t  *r2_hat = pc_wsp;  /* x2_size */

  ctx->m21_vector_multiply(n2_elts, z, ctx->m21_adj, ctx->m21_val,
                           r2_hat);

  const cs_real_t  *r2 = r + ctx->b11_max_size;

# pragma omp parallel for if (n2_elts > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < n2_elts; i2++)
    r2_hat[i2] = r2[i2] - r2_hat[i2];

  /* 3. Solve S z2 = r2_hat (S -> Schur approximation for the m22 block)
     =================================================================== */

  cs_real_t  *z2 = z + ctx->b11_max_size;

  n_iter += _solve_schur_approximation(solver,
                                       ctx->schur_sles,
                                       ctx->schur_matrix,
                                       ctx->m22_mass_diag,  /* inv_m22 */
                                       ctx->schur_scaling,
                                       r2_hat,
                                       z2);

  /* 4. Build r1_tilda = m12.z2 (minus to apply after)
     ========================== */

  cs_real_t  *r1_tilda = pc_wsp + ctx->b22_max_size; /* x1_size */

  cs_array_real_fill_zero(ctx->b11_max_size, r1_tilda);
  ctx->m12_vector_multiply(n2_elts, z2, ctx->m21_adj, ctx->m21_val,
                           r1_tilda);

  if (rset->ifs != nullptr)
    cs_interface_set_sum(rset->ifs,
                         n1_elts,
                         1, false, CS_REAL_TYPE, /* stride, interlaced */
                         r1_tilda);

  /* 5. Solve m11 z1_tilda = r1_tilda
     ================================ */

  cs_real_t  *z1_tilda = pc_wsp + ctx->b11_max_size + ctx->b22_max_size;

  /* No need to scatter the r1_tilda array since it is not used anymore */

  n_iter += _solve_m11_approximation(solver, ctx,
                                     r1_tilda, false, /* no scatter */
                                     z1_tilda, true); /* scatter at exit */

  /* 6. Build r2_tilda = m21.z1_tilda
     ================================ */

  cs_real_t  *r2_tilda = z1_tilda + ctx->b11_max_size; /* x2_size */

  ctx->m21_vector_multiply(n2_elts, z1_tilda, /* z1_tilda has to be "scatter" */
                           ctx->m21_adj, ctx->m21_val,
                           r2_tilda);

  /* 7. Update solutions
     =================== */

  double denum = cs_dot(n2_elts, z2, r2_tilda);
  assert(fabs(denum)>0);
  double zeta =  cs_dot(n2_elts, z2, r2_hat)/denum;

# pragma omp parallel for if (n1_elts > CS_THR_MIN)
  for (cs_lnum_t i1 = 0; i1 < n1_elts; i1++)
    z[i1] += zeta * z1_tilda[i1];

# pragma omp parallel for if (n2_elts > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < n2_elts; i2++)
    z2[i2] *= -zeta;

  return n_iter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply a block diagonal preconditioning: Compute z s.t. P_d z = r
 *        This is equivalent to apply the identity matrix
 *
 *        This prototype complies with the cs_saddle_solver_pc_apply_t function
 *        pointer
 *
 * \param[in]      solver  pointer to a saddle-point solver structure
 * \param[in]      ctx     context pointer for block preconditioner
 * \param[in]      r       rhs of the preconditioning system
 * \param[in, out] z       array to compute
 * \param[in, out] pc_wsp  work space related to the preconditioner
 *
 * \return the number of iteration performed for this step of preconditioning
 */
/*----------------------------------------------------------------------------*/

static int
_no_pc_apply(cs_saddle_solver_t                   *solver,
             cs_saddle_solver_context_block_pcd_t *ctx,
             cs_real_t                            *r,
             cs_real_t                            *z,
             cs_real_t                            *pc_wsp)
{
  CS_NO_WARN_IF_UNUSED(solver);
  CS_NO_WARN_IF_UNUSED(pc_wsp);

  if (z == nullptr)
    return 0;
  assert(r != nullptr);

  /* z <-- r */

  cs_array_real_copy(ctx->b11_max_size + ctx->b22_max_size, r, z);

  return 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the function pointer on the preconditioner (pc) and its workspace
 *
 * \param[in]  solver    saddle-point solver structure
 * \param[in]  ctx       context structure for block preconditioner
 * \param[out] wsp_size  size of the workspace dedicated to the preconditioner
 * \param[out] p_wsp     double pointer on the workspace buffer
 *
 * \return a pointer to the preconditioner function to apply
 */
/*----------------------------------------------------------------------------*/

static cs_saddle_solver_pc_apply_t *
_set_pc_by_block(const cs_saddle_solver_t                   *solver,
                 const cs_saddle_solver_context_block_pcd_t *ctx,
                 cs_lnum_t                                  *wsp_size,
                 cs_real_t                                 **p_wsp)
{
  /* Initialization */

  *wsp_size = 0;
  *p_wsp = nullptr;

  if (ctx == nullptr) /* No preconditioning */
    return  _no_pc_apply;

  const cs_param_saddle_t  *saddlep = solver->param;

  if (saddlep->schur_approx == CS_PARAM_SADDLE_N_SCHUR_APPROX)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid Schur approximation\n", __func__);

  switch (saddlep->precond) {

  case CS_PARAM_SADDLE_PRECOND_NONE:
    return _no_pc_apply; /* No preconditioning */

  case CS_PARAM_SADDLE_PRECOND_DIAG:
    return _diag_schur_pc_apply;

  case CS_PARAM_SADDLE_PRECOND_LOWER:
    *wsp_size = ctx->b22_max_size;
    BFT_MALLOC(*p_wsp, *wsp_size, cs_real_t);

    return _lower_schur_pc_apply;

  case CS_PARAM_SADDLE_PRECOND_SGS:
    *wsp_size = 2*ctx->b11_max_size;
    BFT_MALLOC(*p_wsp, *wsp_size, cs_real_t);

    return _sgs_schur_pc_apply;

  case CS_PARAM_SADDLE_PRECOND_UPPER:
    *wsp_size = ctx->b11_max_size;
    BFT_MALLOC(*p_wsp, *wsp_size, cs_real_t);

    return _upper_schur_pc_apply;

  case CS_PARAM_SADDLE_PRECOND_UZAWA:
    *wsp_size = 2*(ctx->b11_max_size + ctx->b22_max_size);
    BFT_MALLOC(*p_wsp, *wsp_size, cs_real_t);

    return _uza_schur_pc_apply;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid block preconditioner", __func__);
    break;

  } /* Switch on type of preconditioner */

  return nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Test if one needs to do one more ALU iteration. Case of an ALU
 *        algorithm in an incremental formulation.
 *
 * \param[in, out] algo         structure to manage an iterative algorithm
 * \param[in]      l2norm_incr  value of the weighted L2 norm on the increment
 * \param[in]      n2_dofs      number of DoFs for the second block
 * \param[in]      weights      weights to apply in the computation of the norm
 * \param[in]      res2         array of residual values for the second block
 *
 * \return the convergence state
 */
/*----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_alu_incr_cvg_test(cs_iter_algo_t    *algo,
                   cs_real_t          l2norm_incr,
                   const cs_lnum_t    n2_dofs,
                   const cs_real_t   *weights,
                   const cs_real_t   *res2)
{
  /* Compute the new residual based on the norm of the divergence constraint */

  cs_real_t  res2_square = cs_dot_wxx(n2_dofs, weights, res2);
  cs_parall_sum(1, CS_DOUBLE, &res2_square);
  assert(res2_square > -DBL_MIN);
  cs_real_t  l2norm_res2 = sqrt(res2_square);

  cs_iter_algo_update_residual(algo, fmax(l2norm_incr, l2norm_res2));

  /* Update the convergence status */

  cs_sles_convergence_state_t
    cvg_status = cs_iter_algo_update_cvg_tol_auto(algo);

  cs_iter_algo_log_cvg(algo, "# ALUi");

  if (algo->verbosity > 1 && cs_log_default_is_active())
    cs_log_printf(CS_LOG_DEFAULT,
                  "### ALUi.It%03d | l2norm_res2:%10.4e; l2norm_incr:%10.4e\n",
                  cs_iter_algo_get_n_iter(algo), l2norm_res2, l2norm_incr);

  return cvg_status;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Test if one needs to do one more GKB iteration.
 *
 * \param[in]      gamma  value of the augmentation scaling
 * \param[in, out] algo   structure to manage an iterative algorithm
 * \param[in, out] ctx    context structure for the GKB algorithm
 *
 * \return the convergence state
 */
/*----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_gkb_cvg_test(double                          gamma,
              cs_iter_algo_t                 *algo,
              cs_saddle_solver_context_gkb_t *ctx)
{
  /* n = n_algo_iter + 1 since a sum on the square values of zeta is performed
     to estimate the residual in energy norm and the current number of
     iterations has not been updated yet (n_algo_iter = 0 after the first
     resolution done inside _gkb_init_solution()). The number of iterations is
     incremented at the end of the current function, inside the call to
     cs_iter_algo_update_cvg_tol_given() */

  int  n_algo_iter = cs_iter_algo_get_n_iter(algo);

  /* Update the sum of the squared values of zeta (used for renormalization) */

  cs_real_t  z2 = ctx->zeta*ctx->zeta;

  ctx->zeta_square_sum += z2;
  ctx->zeta_array[n_algo_iter % ctx->zeta_size] = z2;

  /* Compute the relative energy norm for the definition of the tolerance
     threshold. The normalization arises from an iterative estimation of the
     initial error in the energy norm */

  int  n = (n_algo_iter < ctx->zeta_size) ? n_algo_iter + 1 : ctx->zeta_size;
  cs_real_t  err2_energy = 0.;
  for (int i = 0; i < n; i++)
    err2_energy += ctx->zeta_array[i];

  double  residual_norm = sqrt(err2_energy);

  /* if n_algo_iter = 0, res0 is automatically set.  For GKB, the first
   * estimation can be rough that's why an update is made at the second
   * resolution too */

  cs_iter_algo_update_residual(algo, residual_norm);
  if (n_algo_iter == 1)
    cs_iter_algo_set_initial_residual(algo, residual_norm);

  /* Compute the tolerance */

  double  tau =
    (gamma > 0) ? sqrt(gamma)*algo->cvg_param.rtol : algo->cvg_param.rtol;
  tau *= sqrt(ctx->zeta_square_sum);

  /* Set the convergence status and increment the number of iteration */

  cs_sles_convergence_state_t  cvg_status = cs_iter_algo_get_cvg_status(algo);

  if (cvg_status == CS_SLES_CONVERGED) /* early convergence */
    cs_iter_algo_update_cvg_tol_given(algo,
                                      fmax(tau, algo->cvg_param.atol));
  else
    cvg_status =
      cs_iter_algo_update_cvg_tol_given(algo,
                                        fmax(tau, algo->cvg_param.atol));

  cs_iter_algo_log_cvg(algo, "# GKB");

  return cvg_status;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the weighted norm related to the (2,2) block
 *
 * \param[in] n2_dofs  number of DoFs (second block)
 * \param[in] w22      weights (diagonal matrix)
 * \param[in] b        base array
 *
 * \return the value of the weighted norm
 */
/*----------------------------------------------------------------------------*/

static inline double
_gkb_block22_weighted_norm(const cs_lnum_t   n2_dofs,
                           const cs_real_t  *w22,
                           const cs_real_t  *b)
{
  double beta2 = cs_dot_wxx(n2_dofs, w22, b);

  /* Parallel synchronization */

  cs_parall_sum(1, CS_DOUBLE, &beta2);
  assert(beta2 > -DBL_MIN);

  /* Keep the value of beta = ||b||_{W22} */

  return sqrt(beta2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Transform the initial saddle-point problem. The x1 unknown is
 *        modified and is stored in x1_tilda as well as the RHS related to the
 *        first block and stored in b1_tilda
 *
 * \param[in, out] solver  pointer to a saddle_point solver structure
 * \param[in, out] ctx     pointer to the GKB context structure
 * \param[in]      x1      initial guess for the (1,1)-block
 */
/*----------------------------------------------------------------------------*/

static void
_gkb_transform_system(cs_saddle_solver_t             *solver,
                      cs_saddle_solver_context_gkb_t *ctx,
                      const cs_real_t                *x1)
{
  assert(ctx != nullptr);
  cs_cdo_system_helper_t  *sh = solver->system_helper;
  cs_iter_algo_t  *algo = solver->algo;

  const cs_param_saddle_t  *saddlep = solver->param;
  const cs_param_saddle_context_gkb_t *ctxp =
    static_cast<cs_param_saddle_context_gkb_t *>(saddlep->context);
  const double  gamma = ctxp->augmentation_scaling;

  const cs_lnum_t  n1_dofs = solver->n1_scatter_dofs;
  const cs_lnum_t  n2_dofs = solver->n2_scatter_dofs;
  const cs_range_set_t  *rset = cs_cdo_system_get_range_set(sh, 0);
  const cs_matrix_t  *m11 = cs_cdo_system_get_matrix(sh, 0);

  /* Transformation of the initial saddle-point system or copy */

  cs_real_t  *rhs1 = sh->rhs_array[0];
  cs_real_t  *rhs2 = sh->rhs_array[1];

  /* Compute rhs_tilda = rhs1 + gamma.M12.M22^-1.rhs2 */

  if (gamma > 0) {

    /* Store temporary inside q = gamma * rhs2 * inv_m22 */

#   pragma omp parallel for if (n2_dofs > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++)
      ctx->q[i2] = rhs2[i2]*ctx->inv_m22[i2];

    /* Build m12q = M12.rhs_tilda */

    cs_array_real_fill_zero(n1_dofs, ctx->m12q);
    ctx->m12_vector_multiply(n2_dofs, ctx->q, ctx->m21_adj, ctx->m21_val,
                             ctx->m12q);

    /* RHS reduction is delayed */

#   pragma omp parallel for if (n1_dofs > CS_THR_MIN)
    for (cs_lnum_t i1 = 0; i1 < n1_dofs; i1++)
      ctx->rhs_tilda[i1] = rhs1[i1] + gamma*ctx->m12q[i1];

  }
  else
    cs_array_real_copy(n1_dofs, rhs1, ctx->rhs_tilda);

  /* Compute w = M11^-1.(rhs1 + gamma.M12.M22^-1.rhs2) */

  cs_real_t  normalization = ctx->square_norm_b11(ctx->rhs_tilda);
  normalization = (fabs(normalization) > FLT_MIN) ? sqrt(normalization) : 1.0;

  cs_array_real_fill_zero(n1_dofs, ctx->v);

  cs_sles_t  *init_sles =
    (ctxp->dedicated_init_sles) ? ctx->init_sles : solver->main_sles;
  assert(init_sles != nullptr);
  cs_param_sles_t  *init_slesp = cs_param_saddle_get_init_sles_param(saddlep);

  int  n_iter = cs_cdo_solve_scalar_system(n1_dofs,
                                           init_slesp,
                                           m11,
                                           rset,
                                           normalization,
                                           true, /* rhs_redux, */
                                           init_sles,
                                           ctx->v,
                                           ctx->rhs_tilda);

  cs_iter_algo_update_inner_iters(algo, n_iter);

  /* Compute x1_tilda := x1 - v */

# pragma omp parallel for if (n1_dofs > CS_THR_MIN)
  for (cs_lnum_t i1 = 0; i1 < n1_dofs; i1++)
    ctx->x1_tilda[i1] = x1[i1] - ctx->v[i1];

  /* Compute rhs_tilda := rhs2 - M21.v */

  ctx->m21_vector_multiply(n2_dofs, ctx->v, ctx->m21_adj, ctx->m21_val,
                           ctx->m21v);

# pragma omp parallel for if (n2_dofs > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++)
    ctx->rhs_tilda[i2] = rhs2[i2] - ctx->m21v[i2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute a first solution array.
 *
 * \param[in, out] solver  pointer to a saddle_point solver structure
 * \param[in, out] ctx     pointer to the GKB context structure
 * \param[in, out] x2      vector solution for the second block
 */
/*----------------------------------------------------------------------------*/

static void
_gkb_init_solution(cs_saddle_solver_t             *solver,
                   cs_saddle_solver_context_gkb_t *ctx,
                   cs_real_t                      *x2)
{
  assert(ctx != nullptr);
  cs_cdo_system_helper_t  *sh = solver->system_helper;
  cs_iter_algo_t  *algo = solver->algo;

  const cs_param_saddle_t  *saddlep = solver->param;
  const cs_lnum_t  n1_dofs = solver->n1_scatter_dofs;
  const cs_lnum_t  n2_dofs = solver->n2_scatter_dofs;
  const cs_range_set_t  *rset = cs_cdo_system_get_range_set(sh, 0);
  const cs_matrix_t  *m11 = cs_cdo_system_get_matrix(sh, 0);

  /* Compute the two following quantities:
   *  beta := ||rhs_tilta||_{M22^-1}
   *  q := {M22^-1}(rhs_tilda)/beta */

  ctx->beta = _gkb_block22_weighted_norm(n2_dofs, ctx->inv_m22, ctx->rhs_tilda);

  /* Store M11^-1.(rhs1 + gamma.M12.M22^-1.rhs2) in rhs_tilda which is not
   * useful anymore */

  if (fabs(ctx->beta) < FLT_MIN) {

    cs_array_real_copy(n1_dofs, ctx->v, ctx->rhs_tilda);
    cs_iter_algo_set_cvg_status(algo, CS_SLES_CONVERGED);
    return;

  }
  else {

    const double  scaling = 1./ctx->beta;
#   pragma omp parallel for if (n2_dofs > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++)
      ctx->q[i2] = scaling * ctx->rhs_tilda[i2] * ctx->inv_m22[i2];

    cs_array_real_copy(n1_dofs, ctx->v, ctx->rhs_tilda);

  }

  /* Solve M11.v = w = M12.q
   * w is updated in the next function */

  cs_array_real_fill_zero(n1_dofs, ctx->w);
  ctx->m12_vector_multiply(n2_dofs, ctx->q, ctx->m21_adj, ctx->m21_val,
                           ctx->w);

  if (rset->ifs != nullptr)
    cs_interface_set_sum(rset->ifs,
                         /* n_elts, stride, interlaced */
                         n1_dofs, 1, false, CS_REAL_TYPE,
                         ctx->w);

  cs_real_t  normalization = ctx->square_norm_b11(ctx->w);
  normalization = (fabs(normalization) > FLT_MIN) ? sqrt(normalization) : 1.0;

  cs_array_real_fill_zero(n1_dofs, ctx->v);

  int  n_iter = cs_cdo_solve_scalar_system(n1_dofs,
                                           saddlep->block11_sles_param,
                                           m11,
                                           rset,
                                           normalization,
                                           false, /* rhs_redux */
                                           solver->main_sles,
                                           ctx->v,
                                           ctx->w);

  cs_iter_algo_update_inner_iters(algo, n_iter);

  ctx->alpha = _scatter_global_dotprod(rset, n1_dofs, ctx->v, ctx->w);
  assert(ctx->alpha > -DBL_MIN);
  ctx->alpha = sqrt(ctx->alpha);

  const double  ov_alpha = 1./ctx->alpha;

  ctx->zeta = ctx->beta * ov_alpha;

  /* Initialize auxiliary vectors and first update of the solution vectors */

# pragma omp parallel for if (n1_dofs > CS_THR_MIN)
  for (cs_lnum_t i1 = 0; i1 < n1_dofs; i1++) {
    ctx->w[i1] *= ov_alpha; /* = M11.v */
    ctx->v[i1] *= ov_alpha;
    ctx->x1_tilda[i1] = ctx->zeta * ctx->v[i1];
  }

# pragma omp parallel for if (n2_dofs > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++) {
    ctx->d[i2] = ctx->q[i2] * ov_alpha;
    x2[i2] = -ctx->zeta * ctx->d[i2];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Test if one needs to do one more ALU iteration. Case of an ALU
 *        algorithm in an incremental formulation.
 *
 * \param[in, out] algo         structure to manage an iterative algorithm
 * \param[in]      l2norm_incr  value of the weighted L2 norm on the increment
 * \param[in]      n2_dofs      number of DoFs for the second block
 * \param[in]      weights      weights to apply in the computation of the norm
 * \param[in]      res2         array of residual values for the second block
 *
 * \return the convergence state
 */
/*----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_uzawa_cg_cvg_test(cs_iter_algo_t  *algo)
{
  /* Update the convergence status */

  cs_sles_convergence_state_t
    cvg_status = cs_iter_algo_update_cvg_tol_auto(algo);

  cs_iter_algo_log_cvg(algo, "# UZACG");

  return cvg_status;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize a cs_saddle_solver_t structure. The context and
 *        algo pointers can be defined in a second step if needed.
 *
 * \param[in] n1_elts         number of elements associated to the (1,1)-block
 * \param[in] n1_dofs_by_elt  number of DoFs by elements in the (1,1)-block
 * \param[in] n2_elts         number of elements associated to the (2,2)-block
 * \param[in] n2_dofs_by_elt  number of DoFs by elements in the (2,2)-block
 * \param[in] saddlep         set of parameters for the saddle-point solver
 * \param[in] sh              pointer to a system helper structure
 * \param[in] main_sles       pointer to the main SLES structure related to
 *                            this saddle-point problem
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

static cs_saddle_solver_t *
_create_saddle_solver(cs_lnum_t                 n1_elts,
                      int                       n1_dofs_by_elt,
                      cs_lnum_t                 n2_elts,
                      int                       n2_dofs_by_elt,
                      const cs_param_saddle_t  *saddlep,
                      cs_cdo_system_helper_t   *sh,
                      cs_sles_t                *main_sles)
{
  if (saddlep == nullptr)
    return nullptr;
  if (sh == nullptr)
    return nullptr;

  cs_saddle_solver_t  *solver = nullptr;

  BFT_MALLOC(solver, 1, cs_saddle_solver_t);

  solver->param = saddlep;        /* shared */
  solver->system_helper = sh;     /* shared */
  solver->main_sles = main_sles;

  solver->do_setup = true;

  /* Scatter view point */

  solver->n1_elts = n1_elts;
  solver->n1_dofs_by_elt = n1_dofs_by_elt;
  solver->n2_elts = n2_elts;
  solver->n2_dofs_by_elt = n2_dofs_by_elt;
  solver->n1_scatter_dofs = n1_elts * n1_dofs_by_elt;
  solver->n2_scatter_dofs = n2_elts;

  if (n1_dofs_by_elt != 1 && n1_dofs_by_elt != 3)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid parameter. This case is not handled up to now.",
              __func__);

  if (n2_dofs_by_elt != 1)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid parameter. This case is not handled up to now.",
              __func__);

  /* Structures set later */

  solver->context = nullptr;
  solver->algo = nullptr;

  /* Monitoring */

  solver->n_calls = 0;
  solver->n_iter_min = INT_MAX;
  solver->n_iter_max = 0;
  solver->n_iter_tot = 0;

  return solver;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the monitoring performance of a saddle-point solver
 *
 * \param[in] solver  pointer to the solver to log
 */
/*----------------------------------------------------------------------------*/

static void
_log_monitoring(const cs_saddle_solver_t  *solver)
{
  if (solver == nullptr)
    return;

  const cs_param_saddle_t  *saddlep = solver->param;

  int  n_calls = solver->n_calls;
  int  n_it_min = solver->n_iter_min;
  int  n_it_max = solver->n_iter_max;
  int  n_it_mean = 0;

  if (n_it_min < 0)
    n_it_min = 0;

  if (n_calls > 0)
    n_it_mean = (int)(solver->n_iter_tot/((unsigned long long)n_calls));

  cs_log_printf(CS_LOG_PERFORMANCE, "\nSummary of resolutions for \"%s\"\n",
                cs_param_saddle_get_name(saddlep));

  cs_log_printf(CS_LOG_PERFORMANCE,
                "\n  Saddle solver type:            %s\n",
                cs_param_saddle_get_type_name(saddlep->solver));

  if (n_it_mean == 0)
    cs_log_printf(CS_LOG_PERFORMANCE, "\n  No resolution\n");

  else
    cs_log_printf(CS_LOG_PERFORMANCE,
                  "  Number of calls:               %12d\n"
                  "  Minimum number of iterations:  %12d\n"
                  "  Maximum number of iterations:  %12d\n"
                  "  Mean number of iterations:     %12d\n",
                  n_calls, n_it_min, n_it_max, n_it_mean);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the number of saddle-point systems which have been added
 *
 * \return the current number of systems
 */
/*----------------------------------------------------------------------------*/

int
cs_saddle_solver_get_n_systems(void)
{
  return cs_saddle_solver_n_systems;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a pointer to a saddle-point solver from its id
 *
 * \param[in] id  id of the saddle-point system
 *
 * \return a pointer to a saddle-point solver structure
 */
/*----------------------------------------------------------------------------*/

cs_saddle_solver_t *
cs_saddle_solver_by_id(int  id)
{
  if (id < 0 || id >= cs_saddle_solver_n_systems)
    return nullptr;
  assert(cs_saddle_solver_systems != nullptr);

  return cs_saddle_solver_systems[id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a new solver for solving a saddle-point problem.
 *
 * \param[in] n1_elts         number of elements associated to the (1,1)-block
 * \param[in] n1_dofs_by_elt  number of DoFs by elements in the (1,1)-block
 * \param[in] n2_elts         number of elements associated to the (2,2)-block
 * \param[in] n2_dofs_by_elt  number of DoFs by elements in the (2,2)-block
 * \param[in] saddlep         set of parameters for the saddle-point solver
 * \param[in] sh              pointer to a system helper structure
 * \param[in] main_sles       pointer to the main SLES structure related to
 *                            this saddle-point problem
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_saddle_solver_t *
cs_saddle_solver_add(cs_lnum_t                 n1_elts,
                     int                       n1_dofs_by_elt,
                     cs_lnum_t                 n2_elts,
                     int                       n2_dofs_by_elt,
                     const cs_param_saddle_t  *saddlep,
                     cs_cdo_system_helper_t   *sh,
                     cs_sles_t                *main_sles)
{
  BFT_REALLOC(cs_saddle_solver_systems,
              cs_saddle_solver_n_systems + 1,
              cs_saddle_solver_t *);

  cs_saddle_solver_t  *solver = _create_saddle_solver(n1_elts,
                                                      n1_dofs_by_elt,
                                                      n2_elts,
                                                      n2_dofs_by_elt,
                                                      saddlep,
                                                      sh,
                                                      main_sles);
  assert(solver != nullptr);

  /* Update static variables */

  cs_saddle_solver_systems[cs_saddle_solver_n_systems] = solver;
  cs_saddle_solver_n_systems += 1;

  return solver;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free all remaining structures related to saddle-point solvers
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_finalize(void)
{
  BFT_FREE(cs_saddle_solver_systems);
  cs_saddle_solver_n_systems = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free a cs_saddle_solver_t structure
 *
 * \param[in, out] p_solver  double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_free(cs_saddle_solver_t  **p_solver)
{
  if (p_solver == nullptr)
    return;

  cs_saddle_solver_t  *solver = *p_solver;

  if (solver == nullptr)
    return;

  const cs_param_saddle_t  *saddlep = solver->param;

  switch (saddlep->solver) {

  case CS_PARAM_SADDLE_SOLVER_ALU:
    {
    cs_saddle_solver_context_alu_t *ctx =
      static_cast<cs_saddle_solver_context_alu_t *>(solver->context);

    cs_saddle_solver_context_alu_free(&ctx);

    BFT_FREE(ctx);
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_GCR:
  case CS_PARAM_SADDLE_SOLVER_MINRES:
    {
    cs_saddle_solver_context_block_pcd_t *ctx =
      static_cast<cs_saddle_solver_context_block_pcd_t *>(solver->context);

    cs_saddle_solver_context_block_pcd_free(&ctx);

    BFT_FREE(ctx);
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_GKB:
    {
    cs_saddle_solver_context_gkb_t *ctx =
      static_cast<cs_saddle_solver_context_gkb_t *>(solver->context);

    cs_saddle_solver_context_gkb_free(&ctx);

    BFT_FREE(ctx);
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_NOTAY_TRANSFORM:
    {
    cs_saddle_solver_context_notay_t *ctx =
      static_cast<cs_saddle_solver_context_notay_t *>(solver->context);

    BFT_FREE(ctx);
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_UZAWA_CG:
    {
    cs_saddle_solver_context_uzawa_cg_t *ctx =
      static_cast<cs_saddle_solver_context_uzawa_cg_t *>(solver->context);

    cs_saddle_solver_context_uzawa_cg_free(&ctx);

    BFT_FREE(ctx);
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_SIMPLE:
    {
    cs_saddle_solver_context_simple_t *ctx =
      static_cast<cs_saddle_solver_context_simple_t *>(solver->context);

    cs_saddle_solver_context_simple_free(&ctx);

    BFT_FREE(ctx);
    }
    break;


  default:
    break; /* Nothing to do */
  }

  BFT_FREE(solver);
  *p_solver = nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free/reset only a part of a cs_saddle_solver_t structure
 *
 * \param[in, out] solver  double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_clean(cs_saddle_solver_t  *solver)
{
  if (solver == nullptr)
    return;

  /* Free data related to the main SLES (a new setup stage will be mandatory) */

  cs_sles_free(solver->main_sles);

  /* Remark: Memory cleaning on the system helper (matrix coefficient/RHS) is
     done in the calling code in order to allow one to do complex algorithm
     where the system is used several times */

  /* Additional memory cleaning depending on the type of solver */

  const cs_param_saddle_t  *saddlep = solver->param;

  switch (saddlep->solver) {

  case CS_PARAM_SADDLE_SOLVER_ALU:
    {
    cs_saddle_solver_context_alu_t *ctx =
      static_cast<cs_saddle_solver_context_alu_t *>(solver->context);

    cs_saddle_solver_context_alu_clean(ctx);
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_GCR:
  case CS_PARAM_SADDLE_SOLVER_MINRES:
    {
    cs_saddle_solver_context_block_pcd_t *ctx =
      static_cast<cs_saddle_solver_context_block_pcd_t *>(solver->context);

    cs_saddle_solver_context_block_pcd_clean(ctx);
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_GKB:
    {
    cs_saddle_solver_context_gkb_t *ctx =
      static_cast<cs_saddle_solver_context_gkb_t *>(solver->context);

    cs_saddle_solver_context_gkb_clean(ctx);
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_UZAWA_CG:
    {
    cs_saddle_solver_context_uzawa_cg_t *ctx =
      static_cast<cs_saddle_solver_context_uzawa_cg_t *>(solver->context);

    cs_saddle_solver_context_uzawa_cg_clean(ctx);
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_SIMPLE:
    {
    cs_saddle_solver_context_simple_t *ctx =
      static_cast<cs_saddle_solver_context_simple_t *>(solver->context);

    cs_saddle_solver_context_simple_clean(ctx);
    }
    break;
  default:
  case CS_PARAM_SADDLE_SOLVER_FGMRES:
  case CS_PARAM_SADDLE_SOLVER_MUMPS:
  case CS_PARAM_SADDLE_SOLVER_NOTAY_TRANSFORM:
    break; /* Nothing to do */
  }

  solver->do_setup = true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the current monitoring state with n_iter
 *
 * \param[in, out] solver  pointer to a saddle solver structure
 * \param[in]      n_iter  number of iterations needed for a new call
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_update_monitoring(cs_saddle_solver_t  *solver,
                                   unsigned             n_iter)
{
  if (solver == nullptr)
    return;

  solver->n_calls += 1;
  solver->n_iter_tot += n_iter;

  if (solver->n_iter_min > n_iter)
    solver->n_iter_min = n_iter;

  if (solver->n_iter_max < n_iter)
    solver->n_iter_max = n_iter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the monitoring performance for all defined saddle-point solvers
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_log_monitoring(void)
{
  if (cs_saddle_solver_n_systems == 0)
    return;

  cs_log_printf(CS_LOG_PERFORMANCE,
                "\nMonitoring of the performance of saddle-point systems\n");

  for (int i = 0; i < cs_saddle_solver_n_systems; i++)
    _log_monitoring(cs_saddle_solver_systems[i]);

  cs_log_separator(CS_LOG_PERFORMANCE);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the lumped matrix the inverse of the diagonal of the
 *        (1,1)-block matrix. The storage of a matrix is in a gather view and
 *        the resulting array is in scatter view.
 *
 * \param[in]      solver     solver for saddle-point problems
 * \param[in]      m11        matrix related to the (1,1) block
 * \param[in]      b11_rset   range set structure for the (1,1) block
 * \param[in, out] xtra_sles  pointer to an extra SLES structure
 * \param[out]     n_iter     number of iterations for this operation
 *
 * \return a pointer to the computed array (scatter view)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_saddle_solver_m11_inv_lumped(cs_saddle_solver_t     *solver,
                                const cs_matrix_t      *m11,
                                const cs_range_set_t   *b11_rset,
                                cs_sles_t              *xtra_sles,
                                int                    *n_iter)
{
  const cs_lnum_t  b11_size = solver->n1_scatter_dofs;

  cs_real_t  *inv_lumped_m11 = nullptr;
  BFT_MALLOC(inv_lumped_m11, b11_size, cs_real_t);
  cs_array_real_fill_zero(b11_size, inv_lumped_m11);

  cs_real_t  *rhs = nullptr;
  BFT_MALLOC(rhs, b11_size, cs_real_t);
  cs_array_real_set_scalar(b11_size, 1., rhs);

  cs_param_sles_t  *slesp = cs_param_saddle_get_xtra_sles_param(solver->param);
  assert(slesp != nullptr);

  /* Normalization of rhs */

  const double  normalization = sqrt(1.0*b11_size);

  /* Solve m11.x = 1 */

  *n_iter = cs_cdo_solve_scalar_system(b11_size,
                                       slesp,
                                       m11,
                                       b11_rset,
                                       normalization,
                                       false, /* No need for rhs reducing */
                                       xtra_sles,
                                       inv_lumped_m11,
                                       rhs);

  /* Partial memory free */

  BFT_FREE(rhs);

  return inv_lumped_m11;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the resulting vector of the operation m12*x2
 *        The stride is equal to 3 for the operator m21 (unassembled)
 *        x2 corresponds to a "scatter" view
 *
 *        This is an update operation. Be careful that the resulting array has
 *        to be initialized.
 *
 * \param[in]      n2_dofs  number of DoFs for x2
 * \param[in]      x2       array for the second set
 * \param[in]      m21_adj  adjacency related to the M21 operator
 * \param[in]      m21_val  array associated to the M21 operator (unassembled)
 * \param[in, out] m12x2    resulting array (have to be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_m12_multiply_vector(cs_lnum_t              n2_dofs,
                                     const cs_real_t       *x2,
                                     const cs_adjacency_t  *m21_adj,
                                     const cs_real_t       *m21_val,
                                     cs_real_t             *m12x2)
{
  assert(n2_dofs == m21_adj->n_elts);

# pragma omp parallel for if (n2_dofs > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++) {

    const cs_real_t  _x2 = x2[i2];
    for (cs_lnum_t j = m21_adj->idx[i2]; j < m21_adj->idx[i2+1]; j++) {

      const cs_real_t  *m21_vals = m21_val + 3*j;
      cs_real_t  *_m12x2 = m12x2 + 3*m21_adj->ids[j];

#     pragma omp critical
      {
        _m12x2[0] += m21_vals[0] * _x2;
        _m12x2[1] += m21_vals[1] * _x2;
        _m12x2[2] += m21_vals[2] * _x2;
      }

    } /* Loop on x1 elements associated to a given x2 DoF */

  } /* Loop on x2 DoFs */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the resulting vector of the operation m12*x2
 *        The stride is equal to 1 for the operator m21 (unassembled)
 *        x2 corresponds to a "scatter" view.
 *
 *        This is an update operation. Be careful that the resulting array has
 *        to be initialized.
 *
 * \param[in]      n2_elts  number of elements for x2 (not DoFs)
 * \param[in]      x2       array for the second set
 * \param[in]      m21_adj  adjacency related to the M21 operator
 * \param[in]      m21_val  array associated to the M21 operator (unassembled)
 * \param[in, out] m12x2    resulting array (have to be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_m12_multiply_scalar(cs_lnum_t              n2_elts,
                                     const cs_real_t       *x2,
                                     const cs_adjacency_t  *m21_adj,
                                     const cs_real_t       *m21_val,
                                     cs_real_t             *m12x2)
{
  assert(n2_elts == m21_adj->n_elts);

# pragma omp parallel for if (n2_elts > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < n2_elts; i2++) {

    const cs_real_t  _x2 = x2[i2];
    for (cs_lnum_t j = m21_adj->idx[i2]; j < m21_adj->idx[i2+1]; j++)
#     pragma omp atomic
      m12x2[m21_adj->ids[j]] += _x2 * m21_val[j];

  } /* Loop on x2 DoFs */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the resulting vector of the operation m21*x1
 *        The stride is equal to 3 for the operator m21 operator
 *        x1 corresponds to a "scatter" view
 *
 * \param[in]      n2_dofs  number of (scatter) DoFs for (2,2)-block
 * \param[in]      x1       array for the first part
 * \param[in]      m21_adj  adjacency related to the M21 operator
 * \param[in]      m21_val  values associated to the M21 operator (unassembled)
 * \param[in, out] m21x1    resulting vector (have to be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_m21_multiply_vector(cs_lnum_t              n2_dofs,
                                     const cs_real_t       *x1,
                                     const cs_adjacency_t  *m21_adj,
                                     const cs_real_t       *m21_val,
                                     cs_real_t             *m21x1)
{
  assert(n2_dofs == m21_adj->n_elts);

# pragma omp parallel for if (n2_dofs > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++) {

    cs_real_t  _m21x1 = 0.;
    for (cs_lnum_t j = m21_adj->idx[i2]; j < m21_adj->idx[i2+1]; j++)
      _m21x1 += cs_math_3_dot_product(m21_val + 3*j, x1 + 3*m21_adj->ids[j]);

    m21x1[i2] = _m21x1;

  } /* Loop on x2 elements */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the resulting vector of the operation m21*x1
 *        The stride is equal to 1 for the operator m21 operator
 *        x1 corresponds to a "scatter" view
 *
 * \param[in]      n2_dofs  number of (scatter) DoFs for (2,2)-block
 * \param[in]      x1       array for the first part
 * \param[in]      m21_adj  adjacency related to the M21 operator
 * \param[in]      m21_val  values associated to the M21 operator (unassembled)
 * \param[in, out] m21x1    resulting vector (have to be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_m21_multiply_scalar(cs_lnum_t              n2_dofs,
                                     const cs_real_t       *x1,
                                     const cs_adjacency_t  *m21_adj,
                                     const cs_real_t       *m21_val,
                                     cs_real_t             *m21x1)
{
  assert(n2_dofs == m21_adj->n_elts);

# pragma omp parallel for if (n2_dofs > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++) {

    cs_real_t  _m21x1 = 0.;
    for (cs_lnum_t j = m21_adj->idx[i2]; j < m21_adj->idx[i2+1]; j++)
      _m21x1 += m21_val[j] * x1[m21_adj->ids[j]];

    m21x1[i2] = _m21x1;

  } /* Loop on x2 elements */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize the context structure for an algorithm related
 *        to the ALU algorithm
 *
 * \param[in, out] solver  pointer to a saddle-point solver structure
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_alu_create(cs_saddle_solver_t  *solver)
{
  if (solver == nullptr)
    return;

  cs_saddle_solver_context_alu_t *ctx = nullptr;
  BFT_MALLOC(ctx, 1, cs_saddle_solver_context_alu_t);

  ctx->inv_m22 = nullptr;
  ctx->res2 = nullptr;
  ctx->m21x1 = nullptr;

  ctx->b1_tilda = nullptr;
  ctx->rhs = nullptr;

  /* Function pointers */

  ctx->square_norm_b11 = nullptr;
  ctx->m12_vector_multiply = nullptr;
  ctx->m21_vector_multiply = nullptr;

  /* Extra SLES if needed */

  const cs_param_saddle_t  *saddlep = solver->param;
  const cs_param_saddle_context_alu_t *ctxp =
    static_cast<cs_param_saddle_context_alu_t *>(saddlep->context);

  ctx->init_sles = nullptr;
  if (ctxp->dedicated_init_sles)
    ctx->init_sles = cs_sles_find_or_add(-1, ctxp->init_sles_param->name);

  /* Rk: Setting the function for computing the b11 norm is done by the calling
   * since it depends on the type of system to solve and especially the
   * localization and stride of the (1,1)-block DoFs
   */

  /* Quick access to the (2,1)-block in an unassembled way. This is also the
   * transposed part of the (1,2)-block by definition of the saddle-point
   * system */

  const cs_cdo_system_helper_t  *sh = solver->system_helper;
  const cs_cdo_system_block_t  *block21 = sh->blocks[1];
  const cs_cdo_system_ublock_t  *b21_ublock =
    static_cast<const cs_cdo_system_ublock_t *>(block21->block_pointer);

  ctx->m21_val = b21_ublock->values;
  ctx->m21_adj = b21_ublock->adjacency;

  /* Set the context structure */

  solver->context = ctx;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free main memory consuming part of the context structure associated
 *        to an ALU algorithm
 *
 * \param[in, out] ctx  pointer to the context structure to clean
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_alu_clean(cs_saddle_solver_context_alu_t *ctx)
{
  if (ctx == nullptr)
    return;

  cs_sles_free(ctx->init_sles);

  BFT_FREE(ctx->inv_m22);
  BFT_FREE(ctx->res2);
  BFT_FREE(ctx->m21x1);
  BFT_FREE(ctx->b1_tilda);
  BFT_FREE(ctx->rhs);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the context structure associated to an ALU algorithm
 *
 * \param[in, out] p_ctx  double pointer to the context structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_alu_free
(
 cs_saddle_solver_context_alu_t  **p_ctx
)
{
  cs_saddle_solver_context_alu_t *ctx = *p_ctx;

  if (ctx == nullptr)
    return;

  cs_saddle_solver_context_alu_clean(ctx);

  BFT_FREE(ctx);
  *p_ctx = nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize the context structure for a block preconditioner
 *        used in combination with a Krylov solver such as MINRES or GCR.
 *
 * \param[in]      b22_max_size  max. size for the second part of unknows
 * \param[in, out] solver        pointer to a saddle-point solver structure
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_block_pcd_create(cs_lnum_t            b22_max_size,
                                          cs_saddle_solver_t  *solver)
{
  if (solver == nullptr)
    return;

  cs_saddle_solver_context_block_pcd_t *ctx = nullptr;
  BFT_MALLOC(ctx, 1, cs_saddle_solver_context_block_pcd_t);

  /* Sanity checks on the system helper have already been done */

  const cs_cdo_system_helper_t  *sh = solver->system_helper;
  const cs_cdo_system_block_t  *vel_block = sh->blocks[0];
  const cs_cdo_system_dblock_t  *vel_dblock =
    static_cast<const cs_cdo_system_dblock_t *>(vel_block->block_pointer);

  ctx->m12_vector_multiply = nullptr;
  ctx->m21_vector_multiply = nullptr;

  /* One assumes that up to now, the (1,1)-block is associated to only one
     matrix */

  ctx->m11 = nullptr; /* To be set later */
  ctx->b11_range_set = vel_dblock->range_set;

  ctx->b11_max_size = 0; /* To be set latter */
  ctx->b22_max_size = b22_max_size;

  /* Quick access to the (2,1)-blockin an unassembled way. This is also the
   * transposed part of the (1,2)-block by definition of the saddle-point
   * system */

  const cs_cdo_system_block_t  *block21 = sh->blocks[1];
  const cs_cdo_system_ublock_t *b21_ublock =
    static_cast<const cs_cdo_system_ublock_t *>(block21->block_pointer);

  ctx->m21_val = b21_ublock->values;
  ctx->m21_adj = b21_ublock->adjacency;

  /* The following members are defined by the higher-level functions since it
     may depends on the discretization and modelling choices */

  ctx->schur_sles = nullptr;
  ctx->xtra_sles = nullptr;
  ctx->schur_matrix = nullptr;
  ctx->schur_scaling = 1.0;
  ctx->schur_diag = nullptr;
  ctx->schur_xtra = nullptr;
  ctx->m11_inv_diag = nullptr;
  ctx->m22_mass_diag = nullptr;

  /* Schur approximation depending on the settings of the block
     preconditioner */

  const cs_param_saddle_t  *saddlep = solver->param;

  switch (saddlep->schur_approx) {

  case CS_PARAM_SADDLE_SCHUR_DIAG_INVERSE:
  case CS_PARAM_SADDLE_SCHUR_LUMPED_INVERSE:
  case CS_PARAM_SADDLE_SCHUR_MASS_SCALED_DIAG_INVERSE:
  case CS_PARAM_SADDLE_SCHUR_MASS_SCALED_LUMPED_INVERSE:
    {
      cs_param_sles_t  *schur_slesp =
        cs_param_saddle_get_schur_sles_param(saddlep);
      assert(schur_slesp != nullptr);

      ctx->schur_sles = cs_sles_find_or_add(-1, schur_slesp->name);

      cs_param_sles_t  *xtra_slesp =
        cs_param_saddle_get_xtra_sles_param(saddlep);

      if (xtra_slesp != nullptr)
        ctx->xtra_sles = cs_sles_find_or_add(-1, xtra_slesp->name);
    }
    break;

  default:
    /*  CS_PARAM_SADDLE_SCHUR_NONE,
        CS_PARAM_SADDLE_SCHUR_IDENTITY,
        CS_PARAM_SADDLE_SCHUR_MASS_SCALED */

    ctx->schur_sles = nullptr;
    break;
  }

  /* Set the context */

  solver->context = ctx;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the main memory-consuming part of the context structure for a
 *        block preconditioner used in combination with a Krylov solver such as
 *        MINRES or GCR.
 *
 * \param[in, out] ctx  pointer to the context structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_block_pcd_clean(
  cs_saddle_solver_context_block_pcd_t *ctx)
{
  if (ctx == nullptr)
    return;

  /* Remove the setup data in SLES. The pointer to the following SLES will be
     still valid */

  cs_sles_free(ctx->schur_sles);
  cs_sles_free(ctx->xtra_sles);

  BFT_FREE(ctx->schur_diag);
  BFT_FREE(ctx->schur_xtra);
  BFT_FREE(ctx->m22_mass_diag);
  BFT_FREE(ctx->m11_inv_diag);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the context structure for a block preconditioner used in
 *        combination with a Krylov solver such as MINRES or GCR.
 *
 * \param[in, out] p_ctx  double pointer to the context structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_block_pcd_free
(
 cs_saddle_solver_context_block_pcd_t  **p_ctx
 )
{
  cs_saddle_solver_context_block_pcd_t *ctx = *p_ctx;

  if (ctx == nullptr)
    return;

  cs_saddle_solver_context_block_pcd_clean(ctx);

  BFT_FREE(ctx);
  *p_ctx = nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize the context structure for the GKB algorithm
 *
 * \param[in, out] solver  pointer to a saddle-point solver structure
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_gkb_create(cs_saddle_solver_t  *solver)
{
  if (solver == nullptr)
    return;

  cs_saddle_solver_context_gkb_t *ctx = nullptr;
  BFT_MALLOC(ctx, 1, cs_saddle_solver_context_gkb_t);

  ctx->alpha = 0.0;
  ctx->beta = 0.0;
  ctx->zeta = 0.0;

  ctx->zeta_size = 0;
  ctx->zeta_array = nullptr;
  ctx->zeta_square_sum = 0;

  ctx->q = nullptr;
  ctx->d = nullptr;
  ctx->m21v = nullptr;
  ctx->inv_m22 = nullptr;
  /* m22 is set later (shared pointer) */

  ctx->w = nullptr;
  ctx->v = nullptr;
  ctx->m12q = nullptr;
  ctx->x1_tilda = nullptr;

  ctx->rhs_tilda = nullptr;

  const cs_param_saddle_t  *saddlep = solver->param;
  cs_param_saddle_context_gkb_t *ctxp =
    static_cast<cs_param_saddle_context_gkb_t *>(saddlep->context);

  /* Extra SLES if needed */

  ctx->init_sles = nullptr;
  if (ctxp->dedicated_init_sles)
    ctx->init_sles = cs_sles_find_or_add(-1, ctxp->init_sles_param->name);

  /* Rk: Setting the function for computing the b11 norm is done by the calling
   * since it depends on the type of system to solve and especially the
   * localization and stride of the (1,1)-block DoFs
   */

  ctx->square_norm_b11 = nullptr;
  ctx->m12_vector_multiply = nullptr;
  ctx->m21_vector_multiply = nullptr;

  /* Quick access to the (2,1)-block in an unassembled way. This is also the
   * transposed part of the (1,2)-block by definition of the saddle-point
   * system */

  const cs_cdo_system_helper_t  *sh = solver->system_helper;
  const cs_cdo_system_block_t  *block21 = sh->blocks[1];
  const cs_cdo_system_ublock_t  *b21_ublock =
    static_cast<const cs_cdo_system_ublock_t *>(block21->block_pointer);

  ctx->m21_val = b21_ublock->values;
  ctx->m21_adj = b21_ublock->adjacency;

  /* Set the context structure */

  solver->context = ctx;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the main memory consuming part of the context structure
 *        associated to the GKB algorithm
 *
 * \param[in, out] ctx  pointer to the context structure to clean
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_gkb_clean(cs_saddle_solver_context_gkb_t *ctx)
{
  if (ctx == nullptr)
    return;

  cs_sles_free(ctx->init_sles);

  ctx->alpha = 0.0;
  ctx->beta = 0.0;
  ctx->zeta = 0.0;

  ctx->zeta_size = 0;
  ctx->zeta_square_sum = 0;

  BFT_FREE(ctx->zeta_array);
  BFT_FREE(ctx->q);
  BFT_FREE(ctx->d);
  BFT_FREE(ctx->m21v);
  BFT_FREE(ctx->inv_m22);

  BFT_FREE(ctx->w);
  BFT_FREE(ctx->v);
  BFT_FREE(ctx->m12q);
  BFT_FREE(ctx->x1_tilda);

  BFT_FREE(ctx->rhs_tilda);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the context structure associated to the GKB algorithm
 *
 * \param[in, out] p_ctx  double pointer to the context structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_gkb_free
(
 cs_saddle_solver_context_gkb_t  **p_ctx
)
{
  cs_saddle_solver_context_gkb_t *ctx = *p_ctx;

  if (ctx == nullptr)
    return;

  cs_saddle_solver_context_gkb_clean(ctx);

  BFT_FREE(ctx);
  *p_ctx = nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize the context structure for the algorithm relying
 *        on the Notay's algebraic transformation
 *
 * \param[in, out] solver  pointer to a saddle-point solver structure
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_notay_create(cs_saddle_solver_t  *solver)
{
  if (solver == nullptr)
    return;

  cs_saddle_solver_context_notay_t *ctx = nullptr;
  BFT_MALLOC(ctx, 1, cs_saddle_solver_context_notay_t);

  /* Function pointer */

  ctx->m12_vector_multiply = nullptr;

  /* Quick access to the (2,1)-block in an unassembled way. This is also the
   * transposed part of the (1,2)-block by definition of the saddle-point
   * system */

  const cs_cdo_system_helper_t  *sh = solver->system_helper;
  const cs_cdo_system_block_t  *block21 = sh->blocks[1];
  const cs_cdo_system_ublock_t  *b21_ublock =
    static_cast<const cs_cdo_system_ublock_t *>(block21->block_pointer);

  ctx->m21_val = b21_ublock->values;
  ctx->m21_adj = b21_ublock->adjacency;

  /* Set the context structure */

  solver->context = ctx;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize the context structure for an algorithm related
 *        to the Uzawa-CG algorithm
 *
 * \param[in]      b22_max_size  max. size for the second part of unknows
 * \param[in, out] solver        pointer to a saddle-point solver structure
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_uzawa_cg_create(cs_lnum_t            b22_max_size,
                                         cs_saddle_solver_t  *solver)
{
  if (solver == nullptr)
    return;

  cs_saddle_solver_context_uzawa_cg_t *ctx = nullptr;
  BFT_MALLOC(ctx, 1, cs_saddle_solver_context_uzawa_cg_t);

  ctx->alpha = 0.0;

  ctx->res2 = nullptr;
  ctx->m21x1 = nullptr;
  ctx->gk = nullptr;

  ctx->b1_tilda = nullptr;
  ctx->dzk = nullptr;
  ctx->rhs = nullptr;

  /* Rk: Setting the function for computing the b11 norm is done by the calling
   * since it depends on the type of system to solve and especially the
   * localization and stride of the (1,1)-block DoFs
   */

  const cs_cdo_system_helper_t  *sh = solver->system_helper;
  const cs_cdo_system_block_t  *b11 = sh->blocks[0];
  const cs_cdo_system_dblock_t  *b11_dblock =
    static_cast<const cs_cdo_system_dblock_t *>(b11->block_pointer);

  /* One assumes that up to now, the (1,1)-block is associated to only one
     matrix */

  ctx->m11 = nullptr; /* To be set later */
  ctx->b11_range_set = b11_dblock->range_set;
  ctx->b11_max_size = 0; /* To be set later */
  ctx->b22_max_size = b22_max_size;

  /* Function pointers */

  ctx->square_norm_b11 = nullptr;
  ctx->m12_vector_multiply = nullptr;
  ctx->m21_vector_multiply = nullptr;

  /* Quick access to the (2,1)-block in an unassembled way. This is also the
   * transposed part of the (1,2)-block by definition of the saddle-point
   * system */

  const cs_cdo_system_block_t  *b21 = sh->blocks[1];
  const cs_cdo_system_ublock_t *b21_ublock =
    static_cast<const cs_cdo_system_ublock_t *>(b21->block_pointer);
  ;

  ctx->m21_val = b21_ublock->values;
  ctx->m21_adj = b21_ublock->adjacency;

  /* The following members are defined by the higher-level functions since it
     may depends on the discretization and modelling choices */

  ctx->schur_diag = nullptr;
  ctx->schur_xtra = nullptr;
  ctx->schur_matrix = nullptr;
  ctx->schur_sles = nullptr;
  ctx->xtra_sles = nullptr;
  ctx->init_sles = nullptr;

  ctx->m11_inv_diag = nullptr;
  ctx->inv_m22 = nullptr;

  /* Schur approximation depending on the settings of the block
     preconditioner */

  const cs_param_saddle_t  *saddlep = solver->param;
  const cs_param_saddle_context_uzacg_t *ctxp =
    static_cast<cs_param_saddle_context_uzacg_t *>(saddlep->context);

  ctx->init_sles = nullptr;
  if (ctxp->dedicated_init_sles)
    ctx->init_sles = cs_sles_find_or_add(-1, ctxp->init_sles_param->name);

  switch (saddlep->schur_approx) {

  case CS_PARAM_SADDLE_SCHUR_DIAG_INVERSE:
  case CS_PARAM_SADDLE_SCHUR_LUMPED_INVERSE:
  case CS_PARAM_SADDLE_SCHUR_MASS_SCALED_DIAG_INVERSE:
  case CS_PARAM_SADDLE_SCHUR_MASS_SCALED_LUMPED_INVERSE:
    {
      cs_param_sles_t  *schur_slesp =
        cs_param_saddle_get_schur_sles_param(saddlep);
      assert(schur_slesp != nullptr);

      ctx->schur_sles = cs_sles_find_or_add(-1, schur_slesp->name);

      cs_param_sles_t  *xtra_slesp =
        cs_param_saddle_get_xtra_sles_param(saddlep);

      if (xtra_slesp != nullptr)
        ctx->xtra_sles = cs_sles_find_or_add(-1, xtra_slesp->name);
    }
    break;

  default:
    /*  CS_PARAM_SADDLE_SCHUR_NONE,
        CS_PARAM_SADDLE_SCHUR_IDENTITY,
        CS_PARAM_SADDLE_SCHUR_MASS_SCALED */

    ctx->schur_sles = nullptr;
    break;
  }

  /* Set the context structure */

  solver->context = ctx;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free main memory consuming part of the context structure associated
 *        to a Uzawa-CG algorithm
 *
 * \param[in, out] ctx  pointer to the context structure to clean
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_uzawa_cg_clean(
  cs_saddle_solver_context_uzawa_cg_t *ctx)
{
  if (ctx == nullptr)
    return;

  BFT_FREE(ctx->res2);
  BFT_FREE(ctx->m21x1);
  BFT_FREE(ctx->gk);
  BFT_FREE(ctx->b1_tilda);
  BFT_FREE(ctx->dzk);
  BFT_FREE(ctx->rhs);

  /* Remove the setup data in SLES. The pointer to the following SLES will be
     still valid */

  cs_sles_free(ctx->schur_sles);
  cs_sles_free(ctx->xtra_sles);
  cs_sles_free(ctx->init_sles);

  BFT_FREE(ctx->schur_diag);
  BFT_FREE(ctx->schur_xtra);
  BFT_FREE(ctx->m11_inv_diag);
  BFT_FREE(ctx->inv_m22);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the context structure associated to a Uzawa-CG algorithm
 *
 * \param[in, out] p_ctx  double pointer to the context structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_uzawa_cg_free
(
 cs_saddle_solver_context_uzawa_cg_t  **p_ctx
)
{
  cs_saddle_solver_context_uzawa_cg_t *ctx = *p_ctx;

  if (ctx == nullptr)
    return;

  cs_saddle_solver_context_uzawa_cg_clean(ctx);

  BFT_FREE(ctx);
  *p_ctx = nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize the context structure for an algorithm related
 *        to the SIMPLE-like algorithm
 *
 * \param[in]      b22_max_size  max. size for the second part of unknows
 * \param[in, out] solver        pointer to a saddle-point solver structure
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_simple_create(cs_lnum_t            b22_max_size,
                                       cs_saddle_solver_t  *solver)
{
  if (solver == nullptr)
    return;

  cs_saddle_solver_context_simple_t *ctx = nullptr;
  BFT_MALLOC(ctx, 1, cs_saddle_solver_context_simple_t);

  ctx->m21x1 = nullptr;
  ctx->rhs = nullptr;
  ctx->b1_tilda = nullptr;

  const cs_cdo_system_helper_t  *sh = solver->system_helper;
  const cs_cdo_system_block_t  *b11 = sh->blocks[0];
  const cs_cdo_system_dblock_t  *b11_dblock =
    static_cast<const cs_cdo_system_dblock_t *>(b11->block_pointer);

  /* One assumes that up to now, the (1,1)-block is associated to only one
     matrix */

  ctx->m11 = nullptr; /* To be set later */
  ctx->b11_range_set = b11_dblock->range_set;
  ctx->b11_max_size = 0; /* To be set later */
  ctx->b22_max_size = b22_max_size;

  /* Function pointers */

  ctx->square_norm_b11 = nullptr;
  ctx->m12_vector_multiply = nullptr;
  ctx->m21_vector_multiply = nullptr;

  /* Quick access to the (2,1)-block in an unassembled way. This is also the
   * transposed part of the (1,2)-block by definition of the saddle-point
   * system */

  const cs_cdo_system_block_t  *b21 = sh->blocks[1];
  const cs_cdo_system_ublock_t *b21_ublock =
    static_cast<const cs_cdo_system_ublock_t *>(b21->block_pointer);
  ;

  ctx->m21_val = b21_ublock->values;
  ctx->m21_adj = b21_ublock->adjacency;

  /* The following members are defined by the higher-level functions since it
     may depends on the discretization and modelling choices */

  ctx->schur_diag = nullptr;
  ctx->schur_xtra = nullptr;
  ctx->schur_matrix = nullptr;
  ctx->schur_sles = nullptr;
  ctx->xtra_sles = nullptr;
  ctx->init_sles = nullptr;

  ctx->m11_inv_diag = nullptr;
  ctx->inv_m22 = nullptr;

  /* Schur approximation depending on the settings of the block
     preconditioner */

  const cs_param_saddle_t  *saddlep = solver->param;
  const cs_param_saddle_context_simple_t *ctxp =
    static_cast<const cs_param_saddle_context_simple_t *>(saddlep->context);

  ctx->init_sles = nullptr;
  if (ctxp->dedicated_init_sles)
    ctx->init_sles = cs_sles_find_or_add(-1, ctxp->init_sles_param->name);

  switch (saddlep->schur_approx) {

  case CS_PARAM_SADDLE_SCHUR_DIAG_INVERSE:
  case CS_PARAM_SADDLE_SCHUR_LUMPED_INVERSE:
  case CS_PARAM_SADDLE_SCHUR_MASS_SCALED_DIAG_INVERSE:
  case CS_PARAM_SADDLE_SCHUR_MASS_SCALED_LUMPED_INVERSE:
    {
      cs_param_sles_t  *schur_slesp =
        cs_param_saddle_get_schur_sles_param(saddlep);
      assert(schur_slesp != nullptr);

      ctx->schur_sles = cs_sles_find_or_add(-1, schur_slesp->name);

      cs_param_sles_t  *xtra_slesp =
        cs_param_saddle_get_xtra_sles_param(saddlep);

      /* Share the same sles of velocity block*/
      if (xtra_slesp != nullptr)
        ctx->xtra_sles = solver->main_sles;
    }
    break;

  default:
    /*  CS_PARAM_SADDLE_SCHUR_NONE,
        CS_PARAM_SADDLE_SCHUR_IDENTITY,
        CS_PARAM_SADDLE_SCHUR_MASS_SCALED */

    ctx->schur_sles = nullptr;
    break;
  }

  /* Set the context structure */

  solver->context = ctx;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free main memory consuming part of the context structure associated
 *        to a SIMPLE algorithm
 *
 * \param[in, out] ctx  pointer to the context structure to clean
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_simple_clean(cs_saddle_solver_context_simple_t *ctx)
{
  if (ctx == nullptr)
    return;

  BFT_FREE(ctx->m21x1);
  BFT_FREE(ctx->b1_tilda);
  BFT_FREE(ctx->rhs);

  /* Remove the setup data in SLES. The pointer to the following SLES will be
     still valid */

  cs_sles_free(ctx->schur_sles);
  cs_sles_free(ctx->xtra_sles);
  cs_sles_free(ctx->init_sles);

  BFT_FREE(ctx->schur_diag);
  BFT_FREE(ctx->schur_xtra);
  BFT_FREE(ctx->m11_inv_diag);
  BFT_FREE(ctx->inv_m22);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the context structure associated to a SIMPLE algorithm
 *
 * \param[in, out] p_ctx  double pointer to the context structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_simple_free
(
 cs_saddle_solver_context_simple_t  **p_ctx
)
{
  cs_saddle_solver_context_simple_t *ctx = *p_ctx;

  if (ctx == nullptr)
    return;

  cs_saddle_solver_context_simple_clean(ctx);

  BFT_FREE(ctx);
  *p_ctx = nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply the Augmented Lagrangian-Uzawa algorithm to a saddle point
 *        problem (the system is stored in a hybrid way). The stride is equal
 *        to 1 for the matrix (db_size[3] = 1) and the vector.
 *
 * \param[in, out] solver  pointer to a saddle_point solver structure
 * \param[in, out] x1      array for the first part
 * \param[in, out] x2      array for the second part
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_alu_incr(cs_saddle_solver_t  *solver,
                          cs_real_t           *x1,
                          cs_real_t           *x2)
{
  assert(solver != nullptr);

  const cs_param_saddle_t  *saddlep = solver->param;
  assert(saddlep->solver == CS_PARAM_SADDLE_SOLVER_ALU);
  const cs_param_saddle_context_alu_t *ctxp =
    static_cast<cs_param_saddle_context_alu_t *>(saddlep->context);
  const double  gamma = ctxp->augmentation_scaling;
  const cs_param_sles_t  *init_slesp = ctxp->init_sles_param;

  cs_iter_algo_t  *algo = solver->algo;

  /* Prepare the solution and rhs arrays given to the solver */

  cs_cdo_system_helper_t  *sh = solver->system_helper;
  cs_saddle_solver_context_alu_t *ctx =
    static_cast<cs_saddle_solver_context_alu_t *>(solver->context);
  assert(ctx != nullptr);

  cs_sles_t  *init_sles =
    ctxp->dedicated_init_sles ? ctx->init_sles : solver->main_sles;
  assert(init_sles != nullptr);

  const cs_lnum_t  n1_dofs = solver->n1_scatter_dofs;
  const cs_lnum_t  n2_dofs = solver->n2_scatter_dofs;
  const cs_range_set_t  *rset = cs_cdo_system_get_range_set(sh, 0);
  const cs_matrix_t  *m11 = cs_cdo_system_get_matrix(sh, 0);

  /* Workspace */

  cs_real_t  *b1 = sh->rhs;
  cs_real_t  *b2 = sh->rhs + n1_dofs;

  /* ------------------ */
  /* --- ALGO BEGIN --- */

  /* Initialization
   * ==============
   *
   * Transformation of the initial right-hand side
   * --------------------------------------------- */

  cs_real_t  *btilda_c = ctx->m21x1;

# pragma omp parallel for if (n2_dofs > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++)
    btilda_c[i2] = ctx->inv_m22[i2]*b2[i2];

  cs_saddle_system_b12_matvec(sh, btilda_c, ctx->b1_tilda,
                              true); /* reset b1_tilda */

  if (rset->ifs != nullptr) {

    cs_interface_set_sum(rset->ifs,
                         /* n_elts, stride, interlaced */
                         n1_dofs, 1, false, CS_REAL_TYPE,
                         ctx->b1_tilda);

    cs_interface_set_sum(rset->ifs,
                         /* n_elts, stride, interlaced */
                         n1_dofs, 1, false, CS_REAL_TYPE,
                         b1);

  }

  /* Build b1_tilda = b1 + gamma*m12.W^-1.b_c */

# pragma omp parallel for if (n1_dofs > CS_THR_MIN)
  for (cs_lnum_t i1 = 0; i1 < n1_dofs; i1++) {
    ctx->b1_tilda[i1] *= gamma;
    ctx->b1_tilda[i1] += b1[i1];
  }

  /* Compute the RHS for the Uzawa system: rhs = b1_tilda - b12.x2 */

  cs_saddle_system_b12_matvec(sh, x2, ctx->rhs, true);

  if (rset->ifs != nullptr)
    cs_interface_set_sum(rset->ifs,
                         /* n_elts, stride, interlaced */
                         n1_dofs, 1, false, CS_REAL_TYPE,
                         ctx->rhs);

# pragma omp parallel for if (n1_dofs > CS_THR_MIN)
  for (cs_lnum_t i1 = 0; i1 < n1_dofs; i1++) {
    ctx->rhs[i1] *= -1;
    ctx->rhs[i1] += ctx->b1_tilda[i1];
  }

  /* Solve AL.u_f = rhs
   * Modify the tolerance in order to be more accurate on this step since
   * the accuracy at this step has an influence on the global accuracy
   */

  cs_real_t  normalization = ctx->square_norm_b11(ctx->rhs);

  normalization = (fabs(normalization) > FLT_MIN) ? sqrt(normalization) : 1.0;

  int  n_iter = cs_cdo_solve_scalar_system(n1_dofs,
                                           init_slesp,
                                           m11,
                                           rset,
                                           normalization,
                                           false, /* rhs_redux */
                                           init_sles,
                                           x1,
                                           ctx->rhs);

  cs_iter_algo_update_inner_iters(algo, n_iter);
  cs_iter_algo_set_normalization(algo, normalization);

#if 0 /* Export the saddle-point system for analysis */
  cs_dbg_binary_dump_system(init_slesp->name, m11, ctx->rhs, x1);
#endif

  /* Main loop */
  /* ========= */

  cs_real_t  *x1_incr = ctx->b1_tilda;
  cs_real_t  l2norm_x1_incr = normalization;

  /* Compute the vector m21x1 = b2 - m21.x1 which corresponds to the residual
   * vector of the second block. It will be used as a stopping criteria */

  cs_saddle_system_b21_matvec(sh, x1, ctx->m21x1);

  /* Update x2 = x2 - gamma * (m21.x1 - b2).
   * Recall that m21 corresponds often to the operator -div
   * Compute the RHS for the Uzawa system:
   *   rhs = -gamma*m12.m22^-1.(m21.x1 - b2)
   */

# pragma omp parallel for if (n2_dofs > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++) {
    ctx->m21x1[i2] -= b2[i2];
    ctx->res2[i2] = ctx->inv_m22[i2] * ctx->m21x1[i2];
    x2[i2] += gamma * ctx->res2[i2];
  }

  while (_alu_incr_cvg_test(algo,
                            l2norm_x1_incr,
                            n2_dofs,
                            ctx->inv_m22, ctx->m21x1) == CS_SLES_ITERATING) {

    /* Continue building the RHS
     * m12_vector_multiply() updates the resulting array that's why one
     * initializes this array first */

    cs_saddle_system_b12_matvec(sh, ctx->res2, ctx->rhs, true);

    if (rset->ifs != nullptr)
      cs_interface_set_sum(rset->ifs,
                           /* n_elts, stride, interlaced */
                           n1_dofs, 1, false, CS_REAL_TYPE,
                           ctx->rhs);

    cs_array_real_scale(n1_dofs, 1, nullptr, -gamma, ctx->rhs);

    /* Solve AL.u_f = rhs */

    cs_array_real_fill_zero(n1_dofs, x1_incr);

    n_iter = cs_cdo_solve_scalar_system(n1_dofs,
                                        saddlep->block11_sles_param,
                                        m11,
                                        rset,
                                        l2norm_x1_incr, /* normalization */
                                        false,          /* rhs_redux */
                                        solver->main_sles,
                                        x1_incr,
                                        ctx->rhs);

#if 0 /* Export the saddle-point system for analysis */
    if (cs_iter_algo_get_n_iter(algo) == 1)
      cs_dbg_binary_dump_system(saddlep->block11_sles_param->name, m11,
                                ctx->rhs, x1_incr);
#endif

    cs_iter_algo_update_inner_iters(algo, n_iter);

    l2norm_x1_incr = ctx->square_norm_b11(x1_incr);
    l2norm_x1_incr = (fabs(l2norm_x1_incr) > FLT_MIN) ?
      sqrt(l2norm_x1_incr) : 1.0;

    /* Update the x1 array */

#   pragma omp parallel for if (n1_dofs > CS_THR_MIN)
    for (cs_lnum_t i1 = 0; i1 < n1_dofs; i1++)
      x1[i1] += x1_incr[i1];

    /* Update the m21.x1 array */

    cs_saddle_system_b21_matvec(sh, x1, ctx->m21x1);

    /* Update the x2 array
     *   x2 = x2 - gamma * (m21.x1 - b2). Recall that m21 = -div
     * Prepare the computation of the RHS for the next Uzawa system:
     *   rhs = -gamma*m12.m22^-1.(m21.x1 - b2)
     */

#   pragma omp parallel for if (n2_dofs > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++) {
      ctx->m21x1[i2] -= b2[i2];
      ctx->res2[i2] = ctx->inv_m22[i2] * ctx->m21x1[i2];
      x2[i2] += gamma * ctx->res2[i2];
    }

  } /* End of ALU iterations */

  /* --- ALGO END --- */
  /* ---------------- */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply the Notay's transformation algorithm to solve a saddle point
 *        problem (the system is stored in a monolithic way).
 *
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] x1      array for the first part
 * \param[in, out] x2      array for the second part
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_notay(cs_saddle_solver_t  *solver,
                       cs_real_t           *x1,
                       cs_real_t           *x2)
{
  assert(solver != nullptr);

  const cs_param_saddle_t  *saddlep = solver->param;
  const cs_cdo_system_helper_t  *sh = solver->system_helper;
  const cs_range_set_t  *rset = cs_cdo_system_get_range_set(sh, 0);
  const cs_matrix_t  *matrix = cs_cdo_system_get_matrix(sh, 0);
  const cs_lnum_t  n_cols = cs_matrix_get_n_columns(matrix);
  const cs_lnum_t  n1_elts = solver->n1_elts;
  const cs_lnum_t  n1_dofs = solver->n1_scatter_dofs;
  const cs_lnum_t  n2_dofs = solver->n2_scatter_dofs;
  const cs_lnum_t  n_scatter_dofs = n1_dofs + n2_dofs;

  assert(n2_dofs == solver->n2_elts);

  /* Prepare the solution and rhs arrays given to the solver */

  cs_real_t  *sol = nullptr;
  BFT_MALLOC(sol, CS_MAX(n_cols, n_scatter_dofs), cs_real_t);

  cs_real_t  *b = nullptr;
  BFT_MALLOC(b, n_scatter_dofs, cs_real_t);

  if (solver->n1_dofs_by_elt == 3 && solver->n2_dofs_by_elt == 1)
    _join_x1_vector_x2_deinterlaced(n1_elts, x1,
                                    n2_dofs, x2,
                                    sh->rhs, sol, b);
  else
    bft_error(__FILE__, __LINE__, 0,
              "%s: Case not handled yet.\n", __func__);

  /* Transformation of system proposed by Notay (The system is also modified in
     the dedicated PETSc hook function) */

# pragma omp parallel for if (CS_THR_MIN > n2_dofs)
  for (cs_lnum_t i2 = n1_dofs; i2 < n_scatter_dofs; i2++)
    b[i2] = -1.0*b[i2];

  /* Handle parallelism: switch from a scatter to a gather view */

  cs_cdo_solve_prepare_system(1,     /* stride */
                              false, /* interlace (managed here) */
                              n_scatter_dofs,
                              rset,
                              true,  /* rhs_redux */
                              sol, b);

  /* Solve the linear solver */

  const cs_param_sles_t  *slesp = saddlep->block11_sles_param;

  /* Warning: rtol associated to the saddle solver is not used in the case of
     the Notay's algebraic transformation */

  cs_real_t  rtol = slesp->cvg_param.rtol;
  cs_field_t  *fld = nullptr;

  cs_solving_info_t  sinfo;
  if (slesp->field_id > -1) {
    fld = cs_field_by_id(slesp->field_id);
    cs_field_get_key_struct(fld, cs_field_key_id("solving_info"), &sinfo);
  }

  sinfo.n_it = 0;
  sinfo.res_norm = DBL_MAX;
  sinfo.rhs_norm = 1.0; /* No renormalization by default (TODO) */

  cs_sles_convergence_state_t  code = cs_sles_solve(solver->main_sles,
                                                    matrix,
                                                    rtol,
                                                    sinfo.rhs_norm,
                                                    &(sinfo.n_it),
                                                    &(sinfo.res_norm),
                                                    b,
                                                    sol,
                                                    0,      /* aux. size */
                                                    nullptr);  /* aux. buffers */

  /* Store metadata for monitoring */

  cs_iter_algo_default_t *algo_ctx =
    static_cast<cs_iter_algo_default_t *>(solver->algo->context);

  algo_ctx->cvg_status = code;
  algo_ctx->normalization = sinfo.rhs_norm;
  algo_ctx->res = sinfo.res_norm;
  algo_ctx->n_algo_iter = sinfo.n_it;

  /* sol is computed and stored in a "gather" view. Switch to a "scatter"
     view */

  cs_range_set_scatter(rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       sol, sol);

#if defined(DEBUG) && !defined(NDEBUG) && CS_SADDLE_SOLVER_DBG > 1
  cs_range_set_scatter(rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       b, b);

  cs_dbg_fprintf_system(saddlep->name,
                        -1,
                        CS_SADDLE_SOLVER_DBG,
                        sol, b, n1_dofs);
#endif

  /* Switch from sol (not interlaced) to x1 and x2
   * Copy the part of the solution array related to x2 */

  cs_array_real_copy(n2_dofs, sol + n1_dofs, x2);

  /* Compute M12(x2) */

  cs_real_t  *m12_x2 = nullptr, *mat_diag = nullptr;

  BFT_MALLOC(m12_x2, n1_dofs, cs_real_t);

  /* Warning: Do not reset the array inside the function since there is a trick
     in the way that the system_helper is built (there are two blocks but the
     system is stored in the first block and the second block is only used for
     the matrix-vector operation in an unssembled way) */

  cs_array_real_fill_zero(n1_dofs, m12_x2);
  cs_saddle_system_b12_matvec(sh, x2, m12_x2, false);

  /* Perform the parallel synchronization. */

  if (rset->ifs != nullptr)
    cs_interface_set_sum(rset->ifs,
                         /* n_elts, stride, interlaced */
                         n1_dofs, 1, false, CS_REAL_TYPE,
                         m12_x2);

  /* Retrieve the diagonal of the matrix in a "scatter" view */

  BFT_MALLOC(mat_diag, n_scatter_dofs, cs_real_t);

  /* diag is stored in a "gather view". Switch to a "scatter view" to make
     the change of variable */

  cs_range_set_scatter(rset,
                       CS_REAL_TYPE,
                       1,         /* treated as scalar-valued up to now */
                       cs_matrix_get_diagonal(matrix), /* gathered view */
                       mat_diag);                      /* scatter view */

  const cs_param_saddle_context_notay_t *ctxp =
    static_cast<const cs_param_saddle_context_notay_t *>(saddlep->context);
  const double  alpha = ctxp->scaling_coef;
  const cs_real_t  *dx = mat_diag,             *solx = sol;
  const cs_real_t  *dy = mat_diag + n1_elts,   *soly = sol + n1_elts;
  const cs_real_t  *dz = mat_diag + 2*n1_elts, *solz = sol + 2*n1_elts;

# pragma omp parallel for if (CS_THR_MIN > n1_elts)             \
  shared(dx, dy, dz, solx, soly, solz) firstprivate(n1_elts)
  for (cs_lnum_t i1 = 0; i1 < n1_elts; i1++) {
    x1[3*i1  ] = solx[i1] - alpha * m12_x2[3*i1  ]/dx[i1];
    x1[3*i1+1] = soly[i1] - alpha * m12_x2[3*i1+1]/dy[i1];
    x1[3*i1+2] = solz[i1] - alpha * m12_x2[3*i1+2]/dz[i1];
  }

  BFT_FREE(m12_x2);
  BFT_FREE(mat_diag);
  BFT_FREE(sol);
  BFT_FREE(b);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply the MINRES algorithm to a saddle point problem (the system is
 *        stored in a hybrid way).
 *        The stride is equal to 1 for the matrix (db_size[3] = 1) and the
 *        vector.
 *
 * \param[in, out] solver  pointer to a saddle_point solver structure
 * \param[in, out] x1      array for the first part
 * \param[in, out] x2      array for the second part
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_minres(cs_saddle_solver_t  *solver,
                        cs_real_t           *x1,
                        cs_real_t           *x2)
{
  assert(solver != nullptr);

  const cs_param_saddle_t  *saddlep = solver->param;

  assert(saddlep->solver == CS_PARAM_SADDLE_SOLVER_MINRES);

  /* Prepare the solution and rhs arrays given to the solver */

  cs_cdo_system_helper_t  *sh = solver->system_helper;
  cs_saddle_solver_context_block_pcd_t *ctx =
    static_cast<cs_saddle_solver_context_block_pcd_t *>(solver->context);
  assert(ctx != nullptr);

  /* Workspace */

  const cs_lnum_t  ssys_size = ctx->b11_max_size + ctx->b22_max_size;
  cs_lnum_t  wsp_size = 7*ssys_size;
  cs_real_t  *wsp = nullptr;

  BFT_MALLOC(wsp, wsp_size, cs_real_t);

  /* Avoid calling an OpenMP initialization on the full buffer (first touch
     paragdim could slows down the memory access when used) */

  memset(wsp, 0, wsp_size*sizeof(cs_real_t));

  cs_real_t  *v     = wsp;
  cs_real_t  *vold  = wsp +   ssys_size;
  cs_real_t  *w     = wsp + 2*ssys_size;
  cs_real_t  *wold  = wsp + 3*ssys_size;
  cs_real_t  *z     = wsp + 4*ssys_size;
  cs_real_t  *zold  = wsp + 5*ssys_size;
  cs_real_t  *mz    = wsp + 6*ssys_size;

  /* Set pointer for the block preconditioning */

  cs_lnum_t  pc_wsp_size = 0;
  cs_real_t  *pc_wsp = nullptr;
  cs_saddle_solver_pc_apply_t  *pc_apply = _set_pc_by_block(solver,
                                                            ctx,
                                                            &pc_wsp_size,
                                                            &pc_wsp);

  /* ------------------ */
  /* --- ALGO BEGIN --- */

  const cs_lnum_t  n1_dofs = solver->n1_scatter_dofs;
  const cs_lnum_t  n2_dofs = solver->n2_scatter_dofs;
  const cs_range_set_t  *rset = ctx->b11_range_set;

  cs_iter_algo_t  *algo = solver->algo;
  cs_real_t  *rhs1 = sh->rhs_array[0];

  /* The RHS is not reduced by default */

  if (rset->ifs != nullptr)
    cs_interface_set_sum(rset->ifs,
                         n1_dofs,
                         1, false, CS_REAL_TYPE, /* stride, interlaced */
                         rhs1);

  /* Compute the first residual: v = b - M.x */

  cs_saddle_system_residual(sh, x1, x2,
                            v,                      /* v1 */
                            v + ctx->b11_max_size); /* v2 */

  /* Apply preconditioning: M.z = v */

  int  last_inner_iter = pc_apply(solver, ctx, v, z, pc_wsp);

  cs_iter_algo_update_inner_iters(algo, last_inner_iter);

  double  residual_norm = _norm(solver, v); /* ||v|| */

  cs_iter_algo_update_residual(algo, residual_norm);
  cs_iter_algo_set_normalization(algo, residual_norm);

  /* dp = eta = <v, z>; beta = sqrt(dp) */

  double  dp = _block_pcd_dot_product(solver, v, z);
  double  beta = sqrt(fabs(dp));
  double  eta = beta;

  /* Initialization */

  double  betaold = 1;
  double  c = 1.0, cold = 1.0, s = 0.0, sold = 0.0;
  cs_sles_convergence_state_t  cvg_status = CS_SLES_ITERATING;

  /* --- MAIN LOOP --- */

  while (cvg_status == CS_SLES_ITERATING) {

    /* z = z * ibeta; */

    assert(fabs(beta) > 0.);
    const double  ibeta = 1./beta;
    _scalar_scaling(solver, ibeta, z);

    /* Compute the matrix-vector product M.z = mz */

    cs_saddle_system_matvec(sh,
                            z, z + ctx->b11_max_size,    /* z1, z2 */
                            mz, mz + ctx->b11_max_size); /* mz1, mz2 */

    /* alpha = <z, mz> */

    const double  alpha =  _block_pcd_dot_product(solver, z, mz);
    const double  alpha_ibeta = alpha * ibeta;
    const double  beta_ibetaold = beta/betaold;

    /* v(k+1) = mz(k) - alpha*v(k) - beta v(k-1) */

#   pragma omp parallel for if (n1_dofs > CS_THR_MIN)
    for (cs_lnum_t i1 = 0; i1 < n1_dofs; i1++) {
      const cs_real_t  _v = v[i1], _vold = vold[i1];
      v[i1] = mz[i1] - alpha_ibeta*_v - beta_ibetaold*_vold;
      vold[i1] = _v;
    }

    cs_real_t  *v2 = v + ctx->b11_max_size;
    cs_real_t  *v2old = vold + ctx->b11_max_size;
    const cs_real_t  *mz2 = mz + ctx->b11_max_size;

#   pragma omp parallel for if (n2_dofs > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++) {
      const cs_real_t  _v = v2[i2], _vold = v2old[i2];
      v2[i2] = mz2[i2] - alpha_ibeta*_v - beta_ibetaold*_vold;
      v2old[i2] = _v;
    }

    /* Apply preconditionning: M.z(k+1) = v(k+1) */

    cs_array_real_copy(ssys_size, z, zold);

    last_inner_iter += pc_apply(solver, ctx, v, z, pc_wsp);

    cs_iter_algo_update_inner_iters(algo, last_inner_iter);

    /* New value for beta: beta = sqrt(<v, z>) */

    betaold = beta;
    beta = sqrt(fabs(_block_pcd_dot_product(solver, v, z)));

    /* QR factorization */

    double rho0 = c*alpha - cold*s*betaold;
    double rho1 = sqrt(rho0*rho0 + beta*beta);
    double rho2 = s*alpha + cold*c*betaold;
    double rho3 = sold*betaold;

    /* Givens rotation (update c and s)*/

    assert(fabs(rho1) > DBL_MIN);
    const double  irho1 = 1./rho1;
    cold = c, sold = s;
    c = rho0*irho1;
    s = beta*irho1;

    /* w(k+1) = irho1 * ( z(k) - rho2*w(k) - rho3 w(k-1) )*/

#   pragma omp parallel for if (n1_dofs > CS_THR_MIN)
    for (cs_lnum_t i1 = 0; i1 < n1_dofs; i1++) {
      const cs_real_t  _w = w[i1], _wold = wold[i1];
      w[i1] = irho1 * (zold[i1] - rho2*_w - rho3*_wold);
      wold[i1] = _w;
    }

    cs_real_t  *w2 = w + ctx->b11_max_size;
    cs_real_t  *w2old = wold + ctx->b11_max_size;
    const cs_real_t  *z2old = zold + ctx->b11_max_size;

#   pragma omp parallel for if (n2_dofs > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++) {
      const cs_real_t  _w = w2[i2], _wold = w2old[i2];
      w2[i2] = irho1 * (z2old[i2] - rho2*_w - rho3*_wold);
      w2old[i2] = _w;
    }

    /* Update the solution vector */
    /* x1(k+1) = x1(k) + c*eta*w(k+1) */

    const double  ceta = c*eta;

#   pragma omp parallel for if (n1_dofs > CS_THR_MIN)
    for (cs_lnum_t i1 = 0; i1 < n1_dofs; i1++)
      x1[i1] = x1[i1] + ceta*w[i1];

#   pragma omp parallel for if (n2_dofs > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++)
      x2[i2] = x2[i2] + ceta*w2[i2];

    /* Compute the current residual */

    residual_norm *= fabs(s);

    cs_iter_algo_update_residual(algo, residual_norm);

    /* Last updates */

    eta = -s*eta;

    /* Check the convergence criteria */

    cvg_status = cs_iter_algo_update_cvg_tol_auto(algo);

    cs_iter_algo_log_cvg(algo, "# Saddle.Minres");

  } /* main loop */

  /* --- ALGO END --- */
  /* ---------------- */

  if (saddlep->verbosity > 1 && cs_log_default_is_active()) {

    /* Compute the real residual norm at exit */

    cs_saddle_system_residual(sh, x1, x2,
                              v,                      /* v1 */
                              v + ctx->b11_max_size); /* v2 */

    residual_norm = _norm(solver, v); /* ||v|| */

    cs_log_printf(CS_LOG_DEFAULT,
                  " %s: Residual norm at exit= %6.4e in %d iterations\n",
                  __func__, residual_norm, cs_iter_algo_get_n_iter(algo));

  }

  /* Free temporary workspace */

  BFT_FREE(wsp);
  BFT_FREE(pc_wsp);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply the GCR algorithm to a saddle point problem (the system is
 *        stored in a hybrid way).
 *        The stride is equal to 1 for the matrix in the (1,1) block
 *        (db_size[3] = 1).
 *
 *        This algorithm is taken from 2010 Notay's paper:
 *        "An aggregation-based algebraic multigrid method" ETNA (vol. 37)
 *
 * \param[in, out] solver  pointer to a saddle_point solver structure
 * \param[in, out] x1      array for the first part
 * \param[in, out] x2      array for the second part
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_gcr(cs_saddle_solver_t  *solver,
                     cs_real_t           *x1,
                     cs_real_t           *x2)
{
  assert(solver != nullptr);

  const cs_param_saddle_t  *saddlep = solver->param;
  assert(saddlep->solver == CS_PARAM_SADDLE_SOLVER_GCR);

  const cs_param_saddle_context_block_krylov_t *ctxp =
    static_cast<cs_param_saddle_context_block_krylov_t *>(saddlep->context);
  const int  restart = ctxp->n_stored_directions;

  /* Prepare the solution and rhs arrays given to the solver */

  cs_saddle_solver_context_block_pcd_t *ctx =
    static_cast<cs_saddle_solver_context_block_pcd_t *>(solver->context);
  assert(ctx != nullptr);

  cs_cdo_system_helper_t  *sh = solver->system_helper;

  /* Workspace */

  const int  triangular_size = (restart*(restart+1))/2;
  double  *gamma = nullptr;
  BFT_MALLOC(gamma, triangular_size, double);
  memset(gamma, 0, triangular_size*sizeof(double));

  double  *alpha = nullptr, *beta = nullptr;
  BFT_MALLOC(alpha, 2*restart, double);
  memset(alpha, 0, 2*restart*sizeof(double));
  beta = alpha + restart;

  const cs_lnum_t  ssys_size = ctx->b11_max_size + ctx->b22_max_size;
  cs_lnum_t  wsp_size = (2 + 2*restart)*ssys_size;
  cs_real_t  *wsp = nullptr;
  BFT_MALLOC(wsp, wsp_size, cs_real_t);
  memset(wsp, 0, wsp_size*sizeof(cs_real_t));

  cs_real_t  *zsave = wsp;
  cs_real_t  *csave = wsp   +   restart*ssys_size;
  cs_real_t  *c_tmp = wsp   + 2*restart*ssys_size; /* size = ssys_size */
  cs_real_t  *r     = c_tmp +   ssys_size;

  /* Set pointer for the block preconditioning */

  cs_lnum_t  pc_wsp_size = 0;
  cs_real_t  *pc_wsp = nullptr;
  cs_saddle_solver_pc_apply_t  *pc_apply = _set_pc_by_block(solver,
                                                            ctx,
                                                            &pc_wsp_size,
                                                            &pc_wsp);

  /* ------------------ */
  /* --- ALGO BEGIN --- */

  const cs_lnum_t  n1_dofs = solver->n1_scatter_dofs;
  const cs_lnum_t  n2_dofs = solver->n2_scatter_dofs;
  const cs_range_set_t  *rset = ctx->b11_range_set;

  cs_iter_algo_t  *algo = solver->algo;
  cs_real_t  *rhs1 = sh->rhs_array[0];

  /* The RHS is not reduced by default */

  if (rset->ifs != nullptr)
    cs_interface_set_sum(rset->ifs,
                         n1_dofs,
                         1, false, CS_REAL_TYPE, /* stride, interlaced */
                         rhs1);

  /* Compute the first residual: r = b - M.x */

  _compute_residual(solver, x1, x2, r);

  double  residual_norm = _norm(solver, r); /* ||r|| */

  cs_iter_algo_update_residual(algo, residual_norm);
  cs_iter_algo_set_normalization(algo, residual_norm);

  /* --- MAIN LOOP --- */

  int  _restart = restart;
  cs_sles_convergence_state_t  cvg_status = CS_SLES_ITERATING;

  while (cvg_status == CS_SLES_ITERATING) {

    if (cs_iter_algo_get_n_iter(algo) > 0)
      _compute_residual(solver, x1, x2, r);

    for (int j = 0; j < _restart; j++) {

      /* Apply preconditioning: M.z = r */

      cs_real_t  *zj = zsave + j*ssys_size;
      int  inner_iter = pc_apply(solver, ctx, r, zj, pc_wsp);

      cs_iter_algo_update_inner_iters(algo, inner_iter);

      /* Compute the matrix-vector product M.z = cj
       * cj plays the role of the temporary buffer during the first part of the
       * algorithm. During the second part, one builds the final state for cj
       */

      cs_real_t  *cj = csave + j*ssys_size;
      _matvec_product(solver, zj, cj);

      for (int i = 0; i < j; i++) {

        /* Compute and store gamma_ij (i < j)/ */

        cs_real_t  *ci = csave + i*ssys_size;
        const double  gamma_ij = _block_pcd_dot_product(solver, ci, cj);

        gamma[_get_id(restart, i,j)] = gamma_ij;

        /* cj = cj - gamma_ij*ci */

        _add_scaled_vector(solver, -gamma_ij, ci, cj);

      } /* i < j */

      /* Case i == j
       * gamma_jj = sqrt(<cj, cj>) */

      const double  gamma_jj = _norm(solver, cj);
      gamma[_get_id(restart, j,j)] = gamma_jj;

      /* Compute cj = 1/gamma_jj * cj  */

      assert(fabs(gamma_jj) > 0);
      _scalar_scaling(solver, 1./gamma_jj, cj);

      /* Compute alpha, store it and use ot to update the residual
       * r = r - alpha_j*c_j */

      const double  _alpha = _block_pcd_dot_product(solver, r, cj);
      alpha[j] = _alpha;

      _add_scaled_vector(solver, -_alpha, cj, r);

      /* New residual norm */

      residual_norm = _norm(solver, r); /* ||r|| */

      cs_iter_algo_update_residual(algo, residual_norm);

      /* Check convergence */

      cvg_status = cs_iter_algo_update_cvg_tol_auto(algo);

      cs_iter_algo_log_cvg(algo, "# Saddle.GCR");

      if (cvg_status != CS_SLES_ITERATING)
        _restart = j + 1; /* Stop the loop on j */

    } /* j < _restart */

    /* Compute the solution vector as a linear combination
     *
     * Compute the beta array s.t. gamma * beta = alpha (one recalls that gamma
     * is an upper left triangular matrix of size _restart)
     */

    for (int ki = _restart - 1; ki > -1; ki--) {
      double  _beta = 0;
      for (int kj = ki + 1; kj < _restart; kj++)
        _beta += gamma[_get_id(restart, ki,kj)] * beta[kj];

      beta[ki] = 1./gamma[_get_id(restart, ki, ki)]*(alpha[ki] - _beta);
    }

    cs_real_t  *update = c_tmp;
    cs_array_real_fill_zero(ssys_size, update);

    for (int k = 0; k < _restart; k++)
      _add_scaled_vector(solver, beta[k], zsave + k*ssys_size, update);

    const cs_real_t  *u1 = update;
#   pragma omp parallel for if (n1_dofs > CS_THR_MIN)
    for (cs_lnum_t i1 = 0; i1 < n1_dofs; i1++)
      x1[i1] += u1[i1];

    const cs_real_t  *u2 = u1 + ctx->b11_max_size;
#   pragma omp parallel for if (n2_dofs > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++)
      x2[i2] += u2[i2];

  } /* Until convergence */

  /* --- ALGO END --- */
  /* ---------------- */

  if (saddlep->verbosity > 1 && cs_log_default_is_active()) {

    /* Compute the real residual norm at exit */

    _compute_residual(solver, x1, x2, r);
    residual_norm = _norm(solver, r); /* ||r|| */
    cs_log_printf(CS_LOG_DEFAULT,
                  " %s: Residual norm at exit= %6.4e in %d iterations\n",
                  __func__, residual_norm, cs_iter_algo_get_n_iter(algo));

  }

  /* Free temporary workspace */

  BFT_FREE(gamma);
  BFT_FREE(alpha);
  BFT_FREE(wsp);
  BFT_FREE(pc_wsp);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply the GKB algorithm to a saddle point problem (the system is
 *        stored in a hybrid way).
 *        The stride is equal to 1 for the matrix in the (1,1) block
 *        (db_size[3] = 1).
 *
 * \param[in, out] solver  pointer to a saddle_point solver structure
 * \param[in, out] x1      array for the first part
 * \param[in, out] x2      array for the second part
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_gkb_inhouse(cs_saddle_solver_t  *solver,
                             cs_real_t           *x1,
                             cs_real_t           *x2)
{
  assert(solver != nullptr);

  const cs_param_saddle_t  *saddlep = solver->param;
  const cs_param_saddle_context_gkb_t *ctxp =
    static_cast<cs_param_saddle_context_gkb_t *>(saddlep->context);
  const double  gamma = ctxp->augmentation_scaling;

  cs_iter_algo_t  *algo = solver->algo;
  cs_cdo_system_helper_t  *sh = solver->system_helper;
  cs_saddle_solver_context_gkb_t *ctx =
    static_cast<cs_saddle_solver_context_gkb_t *>(solver->context);

  assert(saddlep->solver == CS_PARAM_SADDLE_SOLVER_GKB);
  assert(ctx != nullptr);

  /* ------------------ */
  /* --- ALGO BEGIN --- */

  const cs_lnum_t  n1_dofs = solver->n1_scatter_dofs;
  const cs_lnum_t  n2_dofs = solver->n2_scatter_dofs;
  const cs_lnum_t  n2_elts = solver->n2_elts;
  const cs_range_set_t  *rset = cs_cdo_system_get_range_set(sh, 0);
  const cs_matrix_t  *m11 = cs_cdo_system_get_matrix(sh, 0);

  /* Transformation of the initial saddle-point system */

  _gkb_transform_system(solver, ctx, x1);

  /* Initialization of the algorithm (first update of the solution array)
   * Compute the current values for alpha, beta and zeta
   * Compute v, w, d
   */

  _gkb_init_solution(solver, ctx, x2);

  /* Main loop */
  /* ========= */

  while (CS_SLES_ITERATING == _gkb_cvg_test(gamma, algo, ctx)) {

    /* Compute g = M21.v - alpha.M22.q
     *         beta = (g, M22.g)
     */

    ctx->m21_vector_multiply(n2_dofs, ctx->v, ctx->m21_adj, ctx->m21_val,
                             ctx->m21v);

#   pragma omp parallel for if (n2_dofs > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++) {
      ctx->m21v[i2] *= ctx->inv_m22[i2];
      ctx->m21v[i2] -= ctx->alpha * ctx->q[i2];
    }

    /* Update beta */

    ctx->beta = _gkb_block22_weighted_norm(n2_dofs, ctx->m22, ctx->m21v);

    const double  scaling = 1./ctx->beta;
#   pragma omp parallel for if (n2_dofs > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++)
      ctx->q[i2] = scaling * ctx->m21v[i2];

    if (fabs(ctx->beta) < FLT_MIN) {
      cs_iter_algo_set_cvg_status(algo, CS_SLES_CONVERGED);
      break;
    }

    /* Solve M11.v = w = M12.q
     * M12.q is updated in the next function */

    cs_array_real_fill_zero(n1_dofs, ctx->m12q);
    ctx->m12_vector_multiply(n2_elts, ctx->q, ctx->m21_adj, ctx->m21_val,
                             ctx->m12q);

    if (rset->ifs != nullptr)
      cs_interface_set_sum(rset->ifs,
                           /* n_elts, stride, interlaced */
                           n1_dofs, 1, false, CS_REAL_TYPE,
                           ctx->m12q);

#   pragma omp parallel for if (n1_dofs > CS_THR_MIN)
    for (cs_lnum_t i1 = 0; i1 < n1_dofs; i1++) {
      ctx->w[i1] *= -ctx->beta;
      ctx->w[i1] += ctx->m12q[i1];
    }

    cs_real_t  normalization = ctx->alpha;

    cs_array_real_fill_zero(n1_dofs, ctx->v);

    int  n_iter = cs_cdo_solve_scalar_system(n1_dofs,
                                             saddlep->block11_sles_param,
                                             m11,
                                             rset,
                                             normalization,
                                             false, /* rhs_redux */
                                             solver->main_sles,
                                             ctx->v,
                                             ctx->w);

    cs_iter_algo_update_inner_iters(algo, n_iter);

    /* Compute alpha */

    ctx->alpha = _scatter_global_dotprod(rset, n1_dofs, ctx->v, ctx->w);
    assert(ctx->alpha > -DBL_MIN);
    ctx->alpha = sqrt(ctx->alpha);

    const double  ov_alpha = 1./ctx->alpha;

    /* zeta(k+1) = -beta/alpha * zeta(k) */

    ctx->zeta *= -ctx->beta * ov_alpha;

    /* Update vectors and solutions */

#   pragma omp parallel for if (n1_dofs > CS_THR_MIN)
    for (cs_lnum_t i1 = 0; i1 < n1_dofs; i1++) {
      ctx->v[i1] *= ov_alpha;
      ctx->x1_tilda[i1] += ctx->zeta * ctx->v[i1];

      /* Last step: w(k+1) = 1/alpha(k+1) * (M12.q - beta*w(k)) */

      ctx->w[i1] *= ov_alpha;
    }

#   pragma omp parallel for if (n2_dofs > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++) {
      ctx->d[i2] = ov_alpha * (ctx->q[i2] - ctx->beta*ctx->d[i2]);
      x2[i2] -= ctx->zeta * ctx->d[i2];
    }

  } /* End of the main loop on the GKB algorithm */

  /* Return to the initial velocity formulation
   * x1 : = x1_tilda + M11^-1.(rhs1 + gamma.M12.M22^-1.rhs2)
   * where M11^-1.(rhs1 + gamma.M12.M22^-1.rhs2) is stored in rhs_tilda */

# pragma omp parallel for if (n1_dofs > CS_THR_MIN)
  for (cs_lnum_t i1 = 0; i1 < n1_dofs; i1++)
    x1[i1] = ctx->x1_tilda[i1] + ctx->rhs_tilda[i1];

  /* --- ALGO END --- */
  /* ---------------- */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply an (external) solver to solve a saddle point problem (the
 *        system is stored in a monolithic way).
 *
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] x1      array for the first part
 * \param[in, out] x2      array for the second part
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_sles_full_system(cs_saddle_solver_t  *solver,
                                  cs_real_t           *x1,
                                  cs_real_t           *x2)
{
  assert(solver != nullptr);

  const cs_cdo_system_helper_t  *sh = solver->system_helper;
  const cs_range_set_t  *rset = cs_cdo_system_get_range_set(sh, 0);
  const cs_matrix_t  *matrix = cs_cdo_system_get_matrix(sh, 0);
  const cs_lnum_t  n_cols = cs_matrix_get_n_columns(matrix);
  const cs_lnum_t  n1_elts = solver->n1_elts;
  const cs_lnum_t  n1_dofs = solver->n1_scatter_dofs;
  const cs_lnum_t  n2_dofs = solver->n2_scatter_dofs;
  const cs_lnum_t  n_scatter_dofs = n1_dofs + n2_dofs;

  /* Prepare the solution and rhs arrays given to the solver */

  cs_real_t  *sol = nullptr;
  BFT_MALLOC(sol, CS_MAX(n_cols, n_scatter_dofs), cs_real_t);

  cs_real_t  *b = nullptr;
  BFT_MALLOC(b, n_scatter_dofs, cs_real_t);

  assert(solver->n2_elts == n2_dofs);

  if (solver->n1_dofs_by_elt == 3)
    _join_x1_vector_x2_deinterlaced(n1_elts, x1,
                                    n2_dofs, x2,
                                    sh->rhs, sol, b);
  else if (solver->n1_dofs_by_elt == 1)
    _join_x1_scalar_x2_deinterlaced(n1_elts, x1,
                                    n2_dofs, x2,
                                    sh->rhs, sol, b);
  else
    bft_error(__FILE__, __LINE__, 0,
              "%s: Case not handled yet.\n", __func__);

  /* Handle parallelism: switch from a scatter to a gather view */

  cs_cdo_solve_prepare_system(1,     /* stride */
                              false, /* interlace (managed here) */
                              n_scatter_dofs,
                              rset,
                              true,  /* rhs_redux */
                              sol, b);

  /* Solve the linear solver */

  const cs_param_saddle_t  *saddlep = solver->param;
  const cs_param_sles_t  *slesp = saddlep->block11_sles_param;

  cs_real_t  rtol = saddlep->cvg_param.rtol;

  cs_solving_info_t  sinfo;
  cs_field_t  *fld = nullptr;
  if (slesp->field_id > -1) {
    fld = cs_field_by_id(slesp->field_id);
    cs_field_get_key_struct(fld, cs_field_key_id("solving_info"), &sinfo);
  }

  sinfo.n_it = 0;
  sinfo.res_norm = DBL_MAX;
  sinfo.rhs_norm = 1.0; /* No renormalization by default (TODO) */

  cs_sles_convergence_state_t  code = cs_sles_solve(solver->main_sles,
                                                    matrix,
                                                    rtol,
                                                    sinfo.rhs_norm,
                                                    &(sinfo.n_it),
                                                    &(sinfo.res_norm),
                                                    b,
                                                    sol,
                                                    0,      /* aux. size */
                                                    nullptr);  /* aux. buffers */

#if 0 /* Export the saddle-point system for analysis */
  cs_dbg_binary_dump_system(saddlep->name, matrix, b, sol);
#endif

  /* Store metadata for monitoring */

  cs_iter_algo_default_t *algo_ctx =
    static_cast<cs_iter_algo_default_t *>(solver->algo->context);

  algo_ctx->cvg_status = code;
  algo_ctx->normalization = sinfo.rhs_norm;
  algo_ctx->res = sinfo.res_norm;
  algo_ctx->n_algo_iter = sinfo.n_it;

  /* sol is computed and stored in a "gather" view. Switch to a "scatter"
     view */

  cs_range_set_scatter(rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       sol, sol);

#if defined(DEBUG) && !defined(NDEBUG) && CS_SADDLE_SOLVER_DBG > 1
  cs_range_set_scatter(rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       b, b);

  cs_dbg_fprintf_system(saddlep->name,
                        -1,
                        CS_SADDLE_SOLVER_DBG,
                        sol, b, n1_dofs);
#endif

  /* Switch from sol (not interlaced) to x1 and x2
   * Copy the part of the solution array related to x2 */

  cs_array_real_copy(n2_dofs, sol + n1_dofs, x2);

  if (solver->n1_dofs_by_elt == 3) {

    const cs_real_t *solx = sol, *soly = sol + n1_elts, *solz = sol + 2*n1_elts;

#   pragma omp parallel for if (CS_THR_MIN > n1_elts)
    for (cs_lnum_t f = 0; f < n1_elts; f++) {
      x1[3*f  ] = solx[f];
      x1[3*f+1] = soly[f];
      x1[3*f+2] = solz[f];
    }

  }
  else if (solver->n1_dofs_by_elt == 1) /* n1_dofs == n1_elts */
    cs_array_real_copy(n1_dofs, sol, x1);

  if (slesp->field_id > -1)
    cs_field_set_key_struct(fld, cs_field_key_id("solving_info"), &sinfo);

  BFT_FREE(sol);
  BFT_FREE(b);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply the Uzawa-CG algorithm to solve a saddle point problem (the
 *        system is stored in a hybrid way). This algorithm is based on Koko's
 *        paper "Uzawa conjugate gradient method for the Stokes problem: Matlab
 *        implementation with P1-iso-P2/ P1 finite element"
 *
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] x1      array for the first part
 * \param[in, out] x2      array for the second part
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_uzawa_cg(cs_saddle_solver_t  *solver,
                          cs_real_t           *x1,
                          cs_real_t           *x2)
{
  assert(solver != nullptr);

  const cs_param_saddle_t  *saddlep = solver->param;

  cs_iter_algo_t  *algo = solver->algo;
  cs_cdo_system_helper_t  *sh = solver->system_helper;
  cs_saddle_solver_context_uzawa_cg_t *ctx =
    static_cast<cs_saddle_solver_context_uzawa_cg_t *>(solver->context);
  cs_param_saddle_context_uzacg_t *ctxp =
    static_cast<cs_param_saddle_context_uzacg_t *>(saddlep->context);

  const cs_param_sles_t  *init_slesp = ctxp->init_sles_param;

  assert(saddlep->solver == CS_PARAM_SADDLE_SOLVER_UZAWA_CG);
  assert(ctx != nullptr);
  assert(init_slesp != nullptr);

  /* ------------------ */
  /* --- ALGO BEGIN --- */

  const cs_lnum_t  n1_dofs = solver->n1_scatter_dofs;
  const cs_lnum_t  n2_dofs = solver->n2_scatter_dofs;
  const cs_lnum_t  n2_elts = solver->n2_elts;
  const cs_range_set_t  *rset = ctx->b11_range_set;
  const cs_matrix_t  *m11 = ctx->m11;

  cs_real_t  *rhs1 = sh->rhs_array[0];
  cs_real_t  *rhs2 = sh->rhs_array[1];

  /* Compute the first RHS: A.u0 = rhs = b_f - B^t.p_0 to solve */

  cs_array_real_fill_zero(n1_dofs, ctx->rhs);
  ctx->m12_vector_multiply(n2_elts, x2, ctx->m21_adj, ctx->m21_val,
                           ctx->rhs);

  if (rset->ifs != nullptr) {

    cs_interface_set_sum(rset->ifs,
                         /* n_elts, stride, interlaced */
                         n1_dofs, 1, false, CS_REAL_TYPE,
                         ctx->rhs);

    cs_interface_set_sum(rset->ifs,
                         /* n_elts, stride, interlaced */
                         n1_dofs, 1, false, CS_REAL_TYPE,
                         rhs1);

  }

# pragma omp parallel for if (n1_dofs > CS_THR_MIN)
  for (cs_lnum_t i1 = 0; i1 < n1_dofs; i1++)
    ctx->rhs[i1] = rhs1[i1] - ctx->rhs[i1];

  /* Initial normalization from the newly computed rhs */

  double  normalization = ctx->square_norm_b11(ctx->rhs);

  normalization = (fabs(normalization) > FLT_MIN) ? sqrt(normalization) : 1.0;

  /* Compute the first velocity guess
   * Modify the tolerance in order to be more accurate on this step */

  cs_sles_t  *init_sles =
    (ctxp->dedicated_init_sles) ? ctx->init_sles : solver->main_sles;
  assert(init_sles != nullptr);

  int  n_iter = cs_cdo_solve_scalar_system(n1_dofs,
                                           init_slesp,
                                           m11,
                                           rset,
                                           normalization,
                                           false, /* rhs_redux */
                                           init_sles,
                                           x1,
                                           ctx->rhs);

  cs_iter_algo_update_inner_iters(algo, n_iter);

#if 0 /* Export the saddle-point system for analysis */
  cs_dbg_binary_dump_system(saddlep->xtra_sles_param->name, m11, ctx->rhs, x1);
#endif

  /* Set pointers used in this algorithm */

  cs_real_t  *gk = ctx->gk;
  cs_real_t  *dk = ctx->res2;
  cs_real_t  *rk = ctx->m21x1;
  cs_real_t  *wk = ctx->b1_tilda;
  cs_real_t  *dwk = ctx->dzk;
  cs_real_t  *zk = ctx->rhs;

  /* Compute the first residual rk0 (in fact the velocity divergence) */

  ctx->m21_vector_multiply(n2_dofs, x1, ctx->m21_adj, ctx->m21_val, rk);

# pragma omp parallel for if (n2_dofs > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++)
    rk[i2] = rhs2[i2] - rk[i2];

  /* Solve S.zk = rk */

  n_iter = _solve_schur_approximation(solver,
                                      ctx->schur_sles,
                                      ctx->schur_matrix,
                                      ctx->inv_m22,
                                      ctx->alpha,
                                      rk,                  /* r_schur */
                                      zk);                 /* z_schur */

  cs_iter_algo_update_inner_iters(algo, n_iter);

  /* Compute g0 s.t. g0 = alpha zk + nu Mp^-1 r0 */

# pragma omp parallel for if (n2_dofs > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++)
    gk[i2] = ctx->alpha*zk[i2] + ctx->inv_m22[i2]*rk[i2];

  /* dk0 <-- gk0 */

  cs_array_real_copy(n2_dofs, gk, dk);

  double  beta_denum = cs_gdot(n2_dofs, rk, gk);
  double  res_norm =  sqrt(fabs(beta_denum));

  cs_iter_algo_set_normalization(algo, res_norm);
  cs_iter_algo_update_residual(algo, res_norm);

  /* Main loop knowing g0, r0, d0, u0, p0 */
  /* ------------------------------------ */

  while (_uzawa_cg_cvg_test(algo) == CS_SLES_ITERATING) {

    /* Sensitivity step: Compute wk as the solution of M11.wk = M12.dk */

    /* Define the rhs for this system */

    cs_array_real_fill_zero(n1_dofs, ctx->rhs);
    ctx->m12_vector_multiply(n2_elts, dk, ctx->m21_adj, ctx->m21_val, ctx->rhs);

    if (rset->ifs != nullptr)
      cs_interface_set_sum(rset->ifs,
                           /* n_elts, stride, interlaced */
                           n1_dofs, 1, false, CS_REAL_TYPE,
                           ctx->rhs);

    normalization = ctx->square_norm_b11(ctx->rhs);
    normalization = (fabs(normalization) > FLT_MIN) ? sqrt(normalization) : 1.0;

    /* Solve M11.wk = M12.dk (In our context, M12 should be -B^t (and M21 = -B
     * which is the divergence operator) This implies a sign modification
     * during the update step.
     */

    cs_array_real_fill_zero(n1_dofs, wk);

    n_iter = cs_cdo_solve_scalar_system(n1_dofs,
                                        saddlep->block11_sles_param,
                                        m11,
                                        rset,
                                        normalization,
                                        false, /* rhs_redux -->already done */
                                        solver->main_sles,
                                        wk,
                                        ctx->rhs);

    cs_iter_algo_update_inner_iters(algo, n_iter);

     /* M21 is equal to -B s.t. (-B)(-w) = dwk has the right sign */

    ctx->m21_vector_multiply(n2_dofs, wk, ctx->m21_adj, ctx->m21_val,
                             dwk);

  /* Solve the Schur complement approximation: smat.zk = dwk */

    n_iter = _solve_schur_approximation(solver,
                                        ctx->schur_sles,
                                        ctx->schur_matrix,
                                        ctx->inv_m22,        /* inv_m22 */
                                        ctx->alpha,
                                        dwk,                 /* r_schur */
                                        zk);                 /* z_schur */

    cs_iter_algo_update_inner_iters(algo, n_iter);

#   pragma omp parallel for if (n2_dofs > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++)
      zk[i2] = ctx->alpha*zk[i2] + ctx->inv_m22[i2]*dwk[i2];

    /* Updates
     *  - Compute the rho_factor = <rk,gk> / <gk, dwk>
     *  - x1(k+1) = x1(k) + rho_factor * wk  --> --wk
     *  - x2(k+1) = x2(k) - rho_factor * dk
     *  - g(k+1)  = gk    - rho_factor * zk
     */

    double  denum = cs_gdot(n2_dofs, gk, dwk);
    assert(fabs(denum) > 0);
    const double  rho_factor = cs_gdot(n2_dofs, rk, gk) / denum;

#   pragma omp parallel for if (n1_dofs > CS_THR_MIN)
    for (cs_lnum_t i1 = 0; i1 < n1_dofs; i1++)
      x1[i1] += rho_factor * wk[i1]; /* --wk */

#   pragma omp parallel for if (n2_dofs > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++) {
      x2[i2] -= rho_factor * dk[i2];
      gk[i2] -= rho_factor * zk[i2];
      rk[i2] -= rho_factor * dwk[i2];
    }

    /* Conjugate gradient direction: update d(k+1) */

    double  beta_num = cs_gdot(n2_dofs, rk, gk);
    const double  beta_factor = beta_num/beta_denum;
    beta_denum = beta_num; /* Next time it will be used at the denumerator */

    cs_iter_algo_update_residual(algo, sqrt(beta_num));

    /* dk <-- gk + beta_factor * dk */

#   pragma omp parallel for if (n2_dofs > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++)
      dk[i2] = gk[i2] + beta_factor*dk[i2];

  } /* End of main loop */

  /* --- ALGO END --- */
  /* ---------------- */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply the SIMPLE-like algorithm to solve a saddle point problem
 *
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] x1      array for the first part
 * \param[in, out] x2      array for the second part
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_simple(cs_saddle_solver_t  *solver,
                        cs_real_t           *x1,
                        cs_real_t           *x2)
{
  assert(solver != nullptr);

  const cs_param_saddle_t  *saddlep = solver->param;

  cs_iter_algo_t  *algo = solver->algo;
  cs_cdo_system_helper_t  *sh = solver->system_helper;
  cs_saddle_solver_context_simple_t *ctx =
    static_cast<cs_saddle_solver_context_simple_t *>(solver->context);
  cs_param_saddle_context_simple_t *ctxp =
    static_cast<cs_param_saddle_context_simple_t *>(saddlep->context);

  const cs_param_sles_t  *init_slesp = ctxp->init_sles_param;

  assert(saddlep->solver == CS_PARAM_SADDLE_SOLVER_SIMPLE);
  assert(ctx != nullptr);
  assert(init_slesp != nullptr);

  /* ------------------ */
  /* --- ALGO BEGIN --- */

  const cs_lnum_t  n1_dofs = solver->n1_scatter_dofs;
  const cs_lnum_t  n2_dofs = solver->n2_scatter_dofs;
  const cs_lnum_t  n2_elts = solver->n2_elts;
  const cs_range_set_t  *rset = ctx->b11_range_set;
  const cs_matrix_t  *m11 = ctx->m11;

  cs_real_t  *rhs1 = sh->rhs_array[0];
  cs_real_t  *rhs2 = sh->rhs_array[1];

  /* Compute the first RHS: A.u0 = rhs = b_f - B^t.p_0 to solve */

  cs_array_real_fill_zero(n1_dofs, ctx->rhs);
  ctx->m12_vector_multiply(n2_elts, x2, ctx->m21_adj, ctx->m21_val,
                           ctx->rhs);

  if (rset->ifs != nullptr) {

    cs_interface_set_sum(rset->ifs,
                         /* n_elts, stride, interlaced */
                         n1_dofs, 1, false, CS_REAL_TYPE,
                         ctx->rhs);

    cs_interface_set_sum(rset->ifs,
                         /* n_elts, stride, interlaced */
                         n1_dofs, 1, false, CS_REAL_TYPE,
                         rhs1);

  }

# pragma omp parallel for if (n1_dofs > CS_THR_MIN)
  for (cs_lnum_t i1 = 0; i1 < n1_dofs; i1++)
    ctx->rhs[i1] = rhs1[i1] - ctx->rhs[i1];

  /* Initial normalization from the newly computed rhs */

  double  normalization = ctx->square_norm_b11(ctx->rhs);

  normalization = (fabs(normalization) > FLT_MIN) ? sqrt(normalization) : 1.0;

  /* Compute the first velocity guess
   * Modify the tolerance in order to be more accurate on this step */

  cs_sles_t  *init_sles =
    (ctxp->dedicated_init_sles) ? ctx->init_sles : solver->main_sles;
  assert(init_sles != nullptr);

  int  n_iter = cs_cdo_solve_scalar_system(n1_dofs,
                                           init_slesp,
                                           m11,
                                           rset,
                                           normalization,
                                           false, /* rhs_redux */
                                           init_sles,
                                           x1,
                                           ctx->rhs);

  cs_iter_algo_update_inner_iters(algo, n_iter);

  /* Set pointers used in this algorithm */

  cs_real_t  *rk = ctx->m21x1;

  /* Compute the residual rk (in fact the negative velocity divergence) */

  ctx->m21_vector_multiply(n2_dofs, x1, ctx->m21_adj, ctx->m21_val, rk);

# pragma omp parallel for if (n2_dofs > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++)
    rk[i2] = rhs2[i2] + rk[i2];

  /* Solve S.dx2 = rk */

  cs_real_t *dx2 = nullptr;
  BFT_MALLOC(dx2, n2_elts, cs_real_t);

  n_iter = _solve_schur_approximation(solver,
                                      ctx->schur_sles,
                                      ctx->schur_matrix,
                                      ctx->inv_m22,
                                      1.0,
                                      rk,
                                      dx2);

  /* Update x1 and x2 */

  cs_array_real_fill_zero(n1_dofs, ctx->rhs);

  /* Calculate Grad(dx2) */
  ctx->m12_vector_multiply(solver->n2_scatter_dofs,
                           dx2, ctx->m21_adj, ctx->m21_val,
                           ctx->rhs);

  if (rset->ifs != NULL) {

    cs_interface_set_sum(rset->ifs,
                         /* n_elts, stride, interlaced */
                         n1_dofs, 1, false, CS_REAL_TYPE,
                         ctx->rhs);
  }

  cs_real_t *m11_inv = ctx->m11_inv_diag;

# pragma omp parallel for if (n1_dofs > CS_THR_MIN)
  for (cs_lnum_t i1 = 0; i1 < n1_dofs; i1++)
    x1[i1] -= m11_inv[i1]*ctx->rhs[i1];
# pragma omp parallel for if (n2_dofs > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++)
    x2[i2] += dx2[i2];

  BFT_FREE(dx2);

  /* --- ALGO END --- */
  /* ---------------- */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
