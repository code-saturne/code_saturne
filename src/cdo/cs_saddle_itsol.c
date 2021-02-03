/*============================================================================
 * In-house iterative solvers defined by blocks and associated to CDO
 * discretizations
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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

#include "cs_blas.h"
#include "cs_equation_common.h"
#include "cs_log.h"
#include "cs_parall.h"
#include "cs_parameters.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_saddle_itsol.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_saddle_itsol.c
 *
 * \brief  In-house iterative solvers defined by blocks and associated to CDO
 *         discretizations
 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic pointer of function definition for the application of a
 *         preconditioner.
 *         M z = r (M is stored inside ssys with several specificities)
 *
 * \param[in]      ssys     pointer to a cs_saddle_system_t structure
 * \param[in]      sbp      block-preconditioner for the Saddle-point problem
 * \param[in]      r        rhs of the preconditioning system
 * \param[in, out] z        array to compute
 * \param[in, out] pc_wsp   work space related to the preconditioner
 *
 * \return the number of iteration performed for this step of preconditioning
 */
/*----------------------------------------------------------------------------*/

typedef int
(cs_saddle_pc_apply_t)(cs_saddle_system_t          *ssys,
                       cs_saddle_block_precond_t   *sbp,
                       cs_real_t                   *r,
                       cs_real_t                   *z,
                       cs_real_t                   *pc_wsp);

/*============================================================================
 * Private variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Test if one needs one more iteration. The residual criterion is
 *        computed inside the main algorithm.
 *
 * \param[in, out] info     pointer to a Uzawa builder structure
 *
 * \return true (one moe iteration) otherwise false
 */
/*----------------------------------------------------------------------------*/

static bool
_cvg_test(cs_iter_algo_info_t       *info)
{
  /* Increment the number of algo. iterations */
  info->n_algo_iter += 1;

  const cs_real_t  prev_res = info->res;
  const double  epsilon = fmax(info->rtol*info->res0, info->atol);

  /* Set the convergence status */
  if (info->res < epsilon)
    info->cvg = CS_SLES_CONVERGED;

  else if (info->n_algo_iter >= info->n_max_algo_iter)
    info->cvg = CS_SLES_MAX_ITERATION;

  else if (info->res > info->dtol * prev_res)
    info->cvg = CS_SLES_DIVERGED;

  else if (info->res > info->dtol * info->res0)
    info->cvg = CS_SLES_DIVERGED;

  else
    info->cvg = CS_SLES_ITERATING;

  if (info->verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT,
                  "<Krylov.It%02d> res %5.3e | %4d %6d cvg%d | fit.eps %5.3e\n",
                  info->n_algo_iter, info->res, info->last_inner_iter,
                  info->n_inner_iter, info->cvg, epsilon);

  if (info->cvg == CS_SLES_ITERATING)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the scalar multiplication of a vector split into the x1 and
 *        x2 part
 *
 * \param[in]     ssys     pointer to a cs_saddle_system_t structure
 * \param[in]     scalar   scaling coefficient to apply
 * \param[in]     x        vector used for the norm computation
 */
/*----------------------------------------------------------------------------*/

static void
_scalar_scaling(cs_saddle_system_t   *ssys,
                cs_real_t             scalar,
                cs_real_t            *x)
{
  assert(x != NULL);
  cs_real_t  *x1 = x, *x2 = x + ssys->max_x1_size;

# pragma omp parallel for if (ssys->x1_size > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < ssys->x1_size; i++)
    x1[i] *= scalar;

# pragma omp parallel for if (ssys->x2_size > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < ssys->x2_size; i++)
    x2[i] *= scalar;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the canonical dot product between the vectors x and y
 *        The synchronization is performed during the process.
 *
 * \param[in]     ssys    pointer to a cs_saddle_system_t structure
 * \param[in]     x       first vector
 * \param[in]     y       second vector
 *
 * \return the value of the canonical dot product (x,y)
 */
/*----------------------------------------------------------------------------*/

static double
_dot_product(cs_saddle_system_t   *ssys,
             cs_real_t            *x,
             cs_real_t            *y)
{
  double  dp_value = 0.;

  if (x == NULL || y== NULL)
    return dp_value;

  const cs_range_set_t  *rset = ssys->rset;

  cs_real_t  *x1 = x, *x2 = x + ssys->max_x1_size;
  cs_real_t  *y1 = y, *y2 = y + ssys->max_x1_size;

  /* First part x1 and y1 whose DoFs are shared among processes */
  if (rset != NULL) {

    cs_range_set_gather(rset,
                        CS_REAL_TYPE,/* type */
                        1,           /* stride (treated as scalar up to now) */
                        x1,          /* in: size = n_sles_scatter_elts */
                        x1);         /* out: size = n_sles_gather_elts */
    cs_range_set_gather(rset,
                        CS_REAL_TYPE,/* type */
                        1,           /* stride (treated as scalar up to now) */
                        y1,          /* in: size = n_sles_scatter_elts */
                        y1);         /* out: size = n_sles_gather_elts */

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

  dp_value += cs_dot(ssys->x2_size, x2, y2);

  cs_parall_sum(1, CS_DOUBLE, &dp_value);

  return dp_value;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the norm of a vector split into the x1 and x2 part
 *        The synchronization is performed during the process.
 *
 * \param[in]     ssys      pointer to a cs_saddle_system_t structure
 * \param[in]     x       vector used for the norm computation
 *
 * \return the value of the euclidean norm of x
 */
/*----------------------------------------------------------------------------*/

static double
_norm(cs_saddle_system_t   *ssys,
      cs_real_t            *x)
{
  double  n_square_value = 0;

  if (x == NULL)
    return n_square_value;
  assert(ssys != NULL);

  cs_real_t  *x1 = x, *x2 = x + ssys->max_x1_size;
  const cs_range_set_t  *rset = ssys->rset;

  /* Norm for the x1 DoFs (those shared among processes) */
  double  _nx1_square = 0;

  if (rset != NULL) {
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

  /* Norm for the x2 DoFs (not shared so that there is no need to
     synchronize) */
  double  _nx2_square = cs_dot_xx(ssys->x2_size, x2);

  n_square_value = _nx1_square + _nx2_square;
  cs_parall_sum(1, CS_DOUBLE, &n_square_value);
  assert(n_square_value > -DBL_MIN);

  return sqrt(n_square_value);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the resulting vector of the operation m12*x2
 *        The stride is equal to 3 for the operator m21_unassembled
 *        x2 corresponds to a "scatter" view
 *
 * \param[in]      ssys      pointer to a cs_saddle_system_t structure
 * \param[in]      x2        array for the second part
 * \param[in, out] m12x2     resulting vector (have to be allocated)
 */
/*----------------------------------------------------------------------------*/

static void
_m12_vector_multiply(cs_saddle_system_t   *ssys,
                     const cs_real_t      *x2,
                     cs_real_t            *m12x2)
{
  const cs_adjacency_t  *adj = ssys->m21_adjacency;

  assert(ssys->x2_size == adj->n_elts);

# pragma omp parallel for if (ssys->x2_size > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < ssys->x2_size; i2++) {

    const cs_real_t  _x2 = x2[i2];
    for (cs_lnum_t j = adj->idx[i2]; j < adj->idx[i2+1]; j++) {

      const cs_real_t  *m21_vals = ssys->m21_unassembled + 3*j;
      cs_real_t  *_m12x2 = m12x2 + 3*adj->ids[j];

#     pragma omp critical
      {
        _m12x2[0] += m21_vals[0] * _x2;
        _m12x2[1] += m21_vals[1] * _x2;
        _m12x2[2] += m21_vals[2] * _x2;
      }

    } /* Loop on x1 elements associated to a given x2 element */

  } /* Loop on x2 elements */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the resulting vector of the operation m21*x1
 *        The stride is equal to 3 for the operator m21_unassembled
 *        x1 corresponds to a "scatter" view
 *
 * \param[in]      ssys      pointer to a cs_saddle_system_t structure
 * \param[in]      x1        array for the first part
 * \param[in, out] m21x1     resulting vector (have to be allocated)
 */
/*----------------------------------------------------------------------------*/

static void
_m21_vector_multiply(cs_saddle_system_t   *ssys,
                     const cs_real_t      *x1,
                     cs_real_t            *m21x1)
{
  const cs_adjacency_t  *adj = ssys->m21_adjacency;

  assert(ssys->x2_size == adj->n_elts);
# pragma omp parallel for if (ssys->x2_size > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < ssys->x2_size; i2++) {

    cs_real_t  _m21x1 = 0.;
    for (cs_lnum_t j = adj->idx[i2]; j < adj->idx[i2+1]; j++)
      _m21x1 += cs_math_3_dot_product(ssys->m21_unassembled + 3*j,
                                      x1 + 3*adj->ids[j]);

    m21x1[i2] = _m21x1;

  } /* Loop on x2 elements */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the residual divided into two parts
 *        res1 = rhs1 - M11.x1 - M12.x2
 *        res2 = rhs2 - M21.x1
 *
 *        The matrix m11 is represented with 1 block.
 *        The stride is equal to 3 for the operator m21_unassembled
 *
 * \param[in]      ssys      pointer to a cs_saddle_system_t structure
 * \param[in, out] x1        array for the first part
 * \param[in, out] x2        array for the second part
 * \param[out]     res       resulting residual vector
 */
/*----------------------------------------------------------------------------*/

static void
_compute_residual_3(cs_saddle_system_t   *ssys,
                    cs_real_t            *x1,
                    cs_real_t            *x2,
                    cs_real_t            *res)
{
  assert(res != NULL && x1 != NULL && x2 != NULL && ssys != NULL);
  assert(ssys->m21_stride == 3);
  assert(ssys->n_m11_matrices == 1);

  const cs_range_set_t  *rset = ssys->rset;

  cs_real_t  *res1 = res, *res2 = res + ssys->max_x1_size;

  /* Two parts:
   * a) rhs1 - M11.x1 -M12.x2
   * b) rhs2 - M21.x1
   */
  cs_real_t  *m12x2 = NULL;
  BFT_MALLOC(m12x2, ssys->x1_size, cs_real_t);
  memset(m12x2, 0, ssys->x1_size*sizeof(cs_real_t));

  const cs_adjacency_t  *adj = ssys->m21_adjacency;

  assert(ssys->x2_size == adj->n_elts);
# pragma omp parallel for if (ssys->x2_size > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < ssys->x2_size; i2++) {

    const cs_real_t  _x2 = x2[i2];
    cs_real_t  _m21x1 = 0.;
    for (cs_lnum_t j = adj->idx[i2]; j < adj->idx[i2+1]; j++) {

      const cs_lnum_t  shift = 3*adj->ids[j];
      assert(shift < ssys->x1_size);
      const cs_real_t  *m21_vals = ssys->m21_unassembled + 3*j;
      cs_real_t  *_m12x2 = m12x2 + shift;

      _m21x1 += cs_math_3_dot_product(m21_vals, x1 + shift);

#     pragma omp critical
      {
        _m12x2[0] += m21_vals[0] * _x2;
        _m12x2[1] += m21_vals[1] * _x2;
        _m12x2[2] += m21_vals[2] * _x2;
      }

    } /* Loop on x1 elements associated to a given x2 element */

    res2[i2] = ssys->rhs2[i2] - _m21x1;

  } /* Loop on x2 elements */

  const cs_matrix_t  *m11 = ssys->m11_matrices[0];
  cs_matrix_vector_multiply_gs_allocated(rset, m11, x1, res1);

# pragma omp parallel for if (ssys->x1_size > CS_THR_MIN)
  for (cs_lnum_t i1 = 0; i1 < ssys->x1_size; i1++)
    res1[i1] = ssys->rhs1[i1] - res1[i1] - m12x2[i1];

  /* Free memory */
  BFT_FREE(m12x2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the matrix-vector operation divided into two parts
 *        matvec1 = M11.vec1 + M12.vec2
 *        matvec2 = M21.vec1
 *
 *        The stride is equal to 3 for the operator m21_unassembled
 *
 * \param[in]      ssys      pointer to a cs_saddle_system_t structure
 * \param[in, out] vec       array to multiply
 * \param[out]     matvec    resulting vector
 */
/*----------------------------------------------------------------------------*/

static void
_matvec_product(cs_saddle_system_t   *ssys,
                cs_real_t            *vec,
                cs_real_t            *matvec)
{
  /* Sanity checks */
  assert(vec != NULL && matvec != NULL && ssys != NULL);
  assert(ssys->m21_stride == 3);
  assert(ssys->n_m11_matrices == 1);

  const cs_range_set_t  *rset = ssys->rset;

  cs_real_t  *v1 = vec, *v2 = vec + ssys->max_x1_size;
  cs_real_t  *mv1 = matvec, *mv2 = matvec + ssys->max_x1_size;

  /* Two parts:
   * a) mv1 = M11.v1 + M12.v2
   * b) mv2 = M21.v1
   */
  const cs_matrix_t  *m11 = ssys->m11_matrices[0];
  cs_matrix_vector_multiply_gs_allocated(rset, m11, v1, mv1);

  /* 1) M12.v2 and M21.v1 */
  const cs_adjacency_t  *adj = ssys->m21_adjacency;

  cs_real_t  *m12v2 = NULL;
  BFT_MALLOC(m12v2, ssys->x1_size, cs_real_t);
  memset(m12v2, 0, ssys->x1_size*sizeof(cs_real_t));

  assert(ssys->x2_size == adj->n_elts);
# pragma omp parallel for if (ssys->x2_size > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < ssys->x2_size; i2++) {

    const cs_real_t  _v2 = v2[i2];
    cs_real_t  _m21v1 = 0.;
    for (cs_lnum_t j = adj->idx[i2]; j < adj->idx[i2+1]; j++) {

      const cs_lnum_t  shift = 3*adj->ids[j];
      const cs_real_t  *m21_vals = ssys->m21_unassembled + 3*j;
      cs_real_t  *_m12v2 = m12v2 + shift;

      _m21v1 += cs_math_3_dot_product(m21_vals, v1 + shift);

#     pragma omp critical
      {
        _m12v2[0] += m21_vals[0] * _v2;
        _m12v2[1] += m21_vals[1] * _v2;
        _m12v2[2] += m21_vals[2] * _v2;
      }

    } /* Loop on x1 elements associated to a given x2 element */

    mv2[i2] = _m21v1;

  } /* Loop on x2 elements */

  if (rset->ifs != NULL)
    cs_interface_set_sum(rset->ifs,
                         ssys->x1_size,
                         1, false, CS_REAL_TYPE, /* stride, interlaced */
                         mv1);

# pragma omp parallel for if (ssys->x1_size > CS_THR_MIN)
  for (cs_lnum_t i1 = 0; i1 < ssys->x1_size; i1++)
    mv1[i1] += m12v2[i1];

  /* Free temporary memory */
  BFT_FREE(m12v2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Apply a diagonal block preconditioning within a Schur complement
 *         approximation devised in Elman'99 SIAM J. SCI. COMP.
 *
 * \param[in]      ssys    pointer to a cs_saddle_system_t structure
 * \param[in]      sbp     Block-preconditioner for the Saddle-point problem
 * \param[in]      r       rhs of the preconditioning system
 * \param[in, out] z       array to compute
 * \param[in, out] pc_wsp  work space related to the preconditioner
 *
 * \return the number of iterations performed for this step of preconditioning
 */
/*----------------------------------------------------------------------------*/

static int
_elman_schur_pc_apply(cs_saddle_system_t          *ssys,
                      cs_saddle_block_precond_t   *sbp,
                      cs_real_t                   *r,
                      cs_real_t                   *z,
                      cs_real_t                   *pc_wsp)
{
  if (z == NULL)
    return 0;

  /* Sanity checks */
  assert(ssys != NULL && sbp != NULL);
  assert(ssys->n_m11_matrices == 1);
  assert(r != NULL);

  cs_range_set_t  *rset = ssys->rset;
  cs_matrix_t  *m11 = ssys->m11_matrices[0];

  /* Prepare solving (handle parallelism) scatter --> gather transformation
   * stride = 1 for scalar-valued */
  cs_equation_prepare_system(1, ssys->x1_size, m11, rset,
                             false, /* no reduction of the rhs */
                             z,     /* unknown */
                             r);    /* rhs */

  /* Compute the norm of r standing for the rhs (gather view)
   * n_elts[0] corresponds to the number of element in the gather view
   */
  double  r_norm = cs_dot_xx(rset->n_elts[0], r);
  cs_parall_sum(1, CS_DOUBLE, &r_norm);
  r_norm = sqrt(fabs(r_norm));

  /* Solve the linear solver */
  memset(z, 0, sizeof(cs_real_t)*ssys->x1_size);

  cs_solving_info_t  m11_info =
    {.n_it = 0, .res_norm = DBL_MAX, .rhs_norm = r_norm};

  assert(sbp->m11_slesp != NULL);
  cs_param_sles_t  *m11_slesp = sbp->m11_slesp;
  cs_sles_convergence_state_t  code = cs_sles_solve(sbp->m11_sles,
                                                    m11,
                                                    CS_HALO_ROTATION_IGNORE,
                                                    m11_slesp->eps,
                                                    m11_info.rhs_norm,
                                                    &(m11_info.n_it),
                                                    &(m11_info.res_norm),
                                                    r,
                                                    z,
                                                    0,      /* aux. size */
                                                    NULL);  /* aux. buffers */

  /* Output information about the convergence of the resolution */
  int n_iter = m11_info.n_it;
  if (m11_slesp->verbosity > 1)
    cs_log_printf(CS_LOG_DEFAULT, "  <%20s/sles_cvg_code=%-d> n_iters %3d |"
                  " residual % -8.4e | normalization % -8.4e\n",
                  m11_slesp->name, code,
                  m11_info.n_it, m11_info.res_norm, m11_info.rhs_norm);

  /* Move back: gather --> scatter view */
  cs_range_set_scatter(ssys->rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       z, z);
  cs_range_set_scatter(ssys->rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       r, r);

  /* Block m22 (or the Schur approx.) */

  /* Norm for the x2 DoFs (not shared so that there is no need to
     synchronize) */
  cs_real_t  *z2 = z + ssys->max_x1_size;
  cs_real_t  *r2 = r + ssys->max_x1_size;

  r_norm = cs_dot_xx(ssys->x2_size, r2);
  cs_parall_sum(1, CS_DOUBLE, &r_norm);
  r_norm = sqrt(fabs(r_norm));

  /* S^-1 \approx (BBt)^-1 * BABt * (BBt)^-1
   * a) Solve BBt z2 = r2
   */
  memset(z2, 0, sizeof(cs_real_t)*ssys->x2_size);
  n_iter += cs_equation_solve_scalar_cell_system(ssys->x2_size,
                                                 sbp->schur_slesp,
                                                 sbp->schur_matrix,
                                                 r_norm,
                                                 sbp->schur_sles,
                                                 z2,
                                                 r2);

  /* b) m12z2 = m12.z2 */
  cs_real_t  *m12z2 = pc_wsp;

  memset(m12z2, 0, ssys->x1_size*sizeof(cs_real_t));
  /* Update += of m12z2 inside the following function */
  _m12_vector_multiply(ssys, z2, m12z2);

  if (rset->ifs != NULL)
    cs_interface_set_sum(rset->ifs,
                         ssys->x1_size,
                         1, false, CS_REAL_TYPE, /* stride, interlaced */
                         m12z2);

  /* b) m11 * m12z2 = mmz */

  /* scatter view to gather view */
  if (rset != NULL)
    cs_range_set_gather(rset,
                        CS_REAL_TYPE, /* type */
                        1,            /* stride */
                        m12z2,        /* in:  size=n_sles_scatter_elts */
                        m12z2);       /* out: size=n_sles_gather_elts */

  cs_real_t  *mmz = pc_wsp + ssys->max_x1_size;
  memset(mmz, 0, sizeof(cs_real_t)*ssys->max_x1_size);

  cs_matrix_vector_multiply(CS_HALO_ROTATION_IGNORE, m11, m12z2, mmz);

  /* gather to scatter view (i.e. algebraic to mesh view) */
  cs_range_set_scatter(rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       mmz, mmz);

  /* c) m21 * mmz = r2_tilda */
  assert(ssys->x1_size >= ssys->x2_size);
  cs_real_t  *r2_tilda = m12z2;

  _m21_vector_multiply(ssys, mmz, r2_tilda);

  /* d) Solve BBt z2 = r2_tilda */
  r_norm = cs_dot_xx(ssys->x2_size, r2_tilda);
  cs_parall_sum(1, CS_DOUBLE, &r_norm);
  r_norm = sqrt(fabs(r_norm));

  memset(z2, 0, sizeof(cs_real_t)*ssys->x2_size);

  n_iter += cs_equation_solve_scalar_cell_system(ssys->x2_size,
                                                 sbp->schur_slesp,
                                                 sbp->schur_matrix,
                                                 r_norm,
                                                 sbp->schur_sles,
                                                 z2,
                                                 r2_tilda);

  return n_iter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply experimental preconditioning: Compute z s.t. P_d z = r
 *
 * \param[in]      ssys     pointer to a cs_saddle_system_t structure
 * \param[in]      sbp      block-preconditioner for the Saddle-point problem
 * \param[in]      r        rhs of the preconditioning system
 * \param[in, out] z        array to compute
 * \param[in, out] pc_wsp   work space related to the preconditioner
 *
 * \return the number of iterations performed for this step of preconditioning
 */
/*----------------------------------------------------------------------------*/

static int
_gs_schur_pc_apply(cs_saddle_system_t          *ssys,
                   cs_saddle_block_precond_t   *sbp,
                   cs_real_t                   *r,
                   cs_real_t                   *z,
                   cs_real_t                   *pc_wsp)
{
  if (z == NULL)
    return 0;

  /* Sanity checks */
  assert(ssys != NULL && sbp != NULL);
  assert(ssys->n_m11_matrices == 1);
  assert(r != NULL);

  cs_range_set_t  *rset = ssys->rset;
  cs_matrix_t  *m11 = ssys->m11_matrices[0];

  assert(sbp->m11_slesp != NULL);
  assert(sbp->schur_slesp != NULL);

  cs_real_t  *z1_hat = pc_wsp;
  memset(z1_hat, 0, sizeof(cs_real_t)*ssys->x1_size);

  /* Prepare solving (handle parallelism) scatter --> gather transformation
   * stride = 1 for scalar-valued */
  cs_equation_prepare_system(1, ssys->x1_size, m11, rset,
                             false,    /* no reduction of the rhs */
                             z1_hat,   /* unknown */
                             r);       /* rhs */

  /* Compute the norm of r standing for the rhs (gather view)
   * n_elts[0] corresponds to the number of element in the gather view
   */
  double  r_norm = cs_dot_xx(rset->n_elts[0], r);
  cs_parall_sum(1, CS_DOUBLE, &r_norm);
  r_norm = sqrt(fabs(r_norm));

  /* Solve the linear solver */
  cs_solving_info_t  m11_info =
    {.n_it = 0, .res_norm = DBL_MAX, .rhs_norm = r_norm};

  assert(sbp->m11_slesp != NULL);
  cs_param_sles_t  *m11_slesp = sbp->m11_slesp;
  cs_sles_convergence_state_t  code = cs_sles_solve(sbp->m11_sles,
                                                    m11,
                                                    CS_HALO_ROTATION_IGNORE,
                                                    m11_slesp->eps,
                                                    m11_info.rhs_norm,
                                                    &(m11_info.n_it),
                                                    &(m11_info.res_norm),
                                                    r,
                                                    z1_hat,
                                                    0,      /* aux. size */
                                                    NULL);  /* aux. buffers */

  /* Output information about the convergence of the resolution */
  int n_iter = m11_info.n_it;
  if (m11_slesp->verbosity > 1)
    cs_log_printf(CS_LOG_DEFAULT, "  <%20s/sles_cvg_code=%-d> n_iters %3d |"
                  " residual % -8.4e | normalization % -8.4e\n",
                  m11_slesp->name, code,
                  m11_info.n_it, m11_info.res_norm, m11_info.rhs_norm);

  /* Move back: gather --> scatter view */
  cs_range_set_scatter(ssys->rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       z1_hat, z1_hat);
  cs_range_set_scatter(ssys->rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       r, r);

  /* Block m22 (or the Schur approx.) */
  /* Norm for the x2 DoFs (not shared so that there is no need to
     synchronize) */
  cs_real_t  *z2 = z + ssys->max_x1_size;

  /* c) r2_hat = m21 z1_hat
   *    r2_hat = r2_hat - r2 */
  cs_real_t  *r2_hat = pc_wsp + ssys->max_x1_size;

  _m21_vector_multiply(ssys, z1_hat, r2_hat);

  const cs_real_t  *r2 = r + ssys->max_x1_size;

# pragma omp parallel for if (ssys->x2_size > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < ssys->x2_size; i2++)
    r2_hat[i2] = r2_hat[i2] - r2[i2];

  r_norm = cs_dot_xx(ssys->x2_size, r2_hat);
  cs_parall_sum(1, CS_DOUBLE, &r_norm);
  r_norm = sqrt(fabs(r_norm));

  memset(z2, 0, sizeof(cs_real_t)*ssys->x2_size);
  n_iter += cs_equation_solve_scalar_cell_system(ssys->x2_size,
                                                 sbp->schur_slesp,
                                                 sbp->schur_matrix,
                                                 r_norm,
                                                 sbp->schur_sles,
                                                 z2,
                                                 r2_hat);

  /* Third step: New resolution with m11
   * m11 z1 = r1_tilda = m11 z1_hat + m12 z2
   * z1 = -z1 + 2*z1_hat
   */
  cs_real_t  *r1_tilda = r2_hat;

  /* scatter view to gather view */
  cs_range_set_gather(rset,
                      CS_REAL_TYPE, /* type */
                      1,            /* stride */
                      z1_hat,       /* in:  size=n_sles_scatter_elts */
                      z1_hat);      /* out: size=n_sles_gather_elts */

  cs_matrix_vector_multiply(CS_HALO_ROTATION_IGNORE, m11, z1_hat, r1_tilda);

  /* gather to scatter view (i.e. algebraic to mesh view) */
  cs_range_set_scatter(rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       r1_tilda, r1_tilda);

  /* Update += of r1_tilda inside the following function */
  _m12_vector_multiply(ssys, z2, r1_tilda);

  /* Prepare solving (handle parallelism) scatter --> gather transformation
   * stride = 1 for scalar-valued */
  cs_equation_prepare_system(1, ssys->x1_size, m11, rset,
                             false,      /* no reduction of the rhs */
                             z,          /* unknown */
                             r1_tilda);  /* rhs */

  /* Compute the norm of r standing for the rhs (gather view)
   * n_elts[0] corresponds to the number of element in the gather view
   */
  r_norm = cs_dot_xx(rset->n_elts[0], r1_tilda);
  cs_parall_sum(1, CS_DOUBLE, &r_norm);
  r_norm = sqrt(fabs(r_norm));

  /* Second solve of the linear system m11 */
  m11_info.n_it = 0, m11_info.res_norm = DBL_MAX, m11_info.rhs_norm = r_norm;

  code = cs_sles_solve(sbp->m11_sles,
                       m11,
                       CS_HALO_ROTATION_IGNORE,
                       m11_slesp->eps,
                       m11_info.rhs_norm,
                       &(m11_info.n_it),
                       &(m11_info.res_norm),
                       r1_tilda,
                       z,
                       0,      /* aux. size */
                       NULL);  /* aux. buffers */

  /* Output information about the convergence of the resolution */
  n_iter += m11_info.n_it;
  if (m11_slesp->verbosity > 1)
    cs_log_printf(CS_LOG_DEFAULT, "  <%20s/sles_cvg_code=%-d> n_iters %3d |"
                  " residual % -8.4e | normalization % -8.4e\n",
                  m11_slesp->name, code,
                  m11_info.n_it, m11_info.res_norm, m11_info.rhs_norm);

  /* Update z1 = -z1 + 2 z1_hat (still in gather wiew for both arrays) */
# pragma omp parallel for if (rset->n_elts[0] > CS_THR_MIN)
  for (cs_lnum_t i1 = 0; i1 < rset->n_elts[0]; i1++)
    z[i1] = 2*z1_hat[i1] - z[i1];

  /* Move back: gather --> scatter view */
  cs_range_set_scatter(ssys->rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       z, z);

  return n_iter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply a block diagonal preconditioning: Compute z s.t. P_d z = r
 *
 * \param[in]      ssys    pointer to a cs_saddle_system_t structure
 * \param[in]      sbp     Block-preconditioner for the Saddle-point problem
 * \param[in]      r       rhs of the preconditioning system
 * \param[in, out] z       array to compute
 * \param[in, out] pc_wsp  work space related to the preconditioner
 *
 * \return the number of iterations performed for this step of preconditioning
 */
/*----------------------------------------------------------------------------*/

static int
_diag_schur_pc_apply(cs_saddle_system_t          *ssys,
                     cs_saddle_block_precond_t   *sbp,
                     cs_real_t                   *r,
                     cs_real_t                   *z,
                     cs_real_t                   *pc_wsp)
{
  if (z == NULL)
    return 0;

  /* Sanity checks */
  assert(ssys != NULL && sbp != NULL);
  assert(ssys->n_m11_matrices == 1);
  assert(r != NULL);
  assert(sbp->schur_type == CS_PARAM_SCHUR_IDENTITY ||
         sbp->schur_type == CS_PARAM_SCHUR_DIAG_INVERSE ||
         sbp->schur_type == CS_PARAM_SCHUR_LUMPED_INVERSE);

  cs_range_set_t  *rset = ssys->rset;
  cs_matrix_t  *m11 = ssys->m11_matrices[0];

  /* Prepare solving (handle parallelism) scatter --> gather transformation
   * stride = 1 for scalar-valued */
  cs_equation_prepare_system(1, ssys->x1_size, m11, rset,
                             false, /* no reduction of the rhs */
                             z,     /* unknown */
                             r);    /* rhs */

  /* Compute the norm of r standing for the rhs (gather view)
   * n_elts[0] corresponds to the number of element in the gather view
   */
  double  r_norm = cs_dot_xx(rset->n_elts[0], r);
  cs_parall_sum(1, CS_DOUBLE, &r_norm);
  r_norm = sqrt(fabs(r_norm));

  /* Solve the linear solver */
  memset(z, 0, sizeof(cs_real_t)*ssys->x1_size);

  cs_solving_info_t  m11_info =
    {.n_it = 0, .res_norm = DBL_MAX, .rhs_norm = r_norm};

  assert(sbp->m11_slesp != NULL);
  cs_param_sles_t  *m11_slesp = sbp->m11_slesp;
  cs_sles_convergence_state_t  code = cs_sles_solve(sbp->m11_sles,
                                                    m11,
                                                    CS_HALO_ROTATION_IGNORE,
                                                    m11_slesp->eps,
                                                    m11_info.rhs_norm,
                                                    &(m11_info.n_it),
                                                    &(m11_info.res_norm),
                                                    r,
                                                    z,
                                                    0,      /* aux. size */
                                                    NULL);  /* aux. buffers */

  /* Output information about the convergence of the resolution */
  int n_iter = m11_info.n_it;
  if (m11_slesp->verbosity > 1)
    cs_log_printf(CS_LOG_DEFAULT, "  <%20s/sles_cvg_code=%-d> n_iters %3d |"
                  " residual % -8.4e | normalization % -8.4e\n",
                  m11_slesp->name, code,
                  m11_info.n_it, m11_info.res_norm, m11_info.rhs_norm);

  /* Move back: gather --> scatter view */
  cs_range_set_scatter(ssys->rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       z, z);
  cs_range_set_scatter(ssys->rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       r, r);

  /* Block m22 (or the Schur approx.) */
  if (sbp->schur_sles == NULL ||
      sbp->schur_type == CS_PARAM_SCHUR_IDENTITY) /* Precond = Identity */
    memcpy(z + ssys->max_x1_size,
           r + ssys->max_x1_size, sizeof(cs_real_t)*ssys->x2_size);

  else if (sbp->schur_type == CS_PARAM_SCHUR_MASS_SCALED) {

    assert(sbp->mass22_diag != NULL); /* inverse of the scaled mass matrix */

    const cs_real_t  *r2 = r + ssys->max_x1_size;
    cs_real_t  *z2 = z + ssys->max_x1_size;
#   pragma omp parallel for if (ssys->x2_size > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < ssys->x2_size; i2++)
      z2[i2] = sbp->mass22_diag[i2]*r2[i2];

  }
  else {

    /* Norm for the x2 DoFs (not shared so that there is no need to
       synchronize) */
    cs_real_t  *z2 = z + ssys->max_x1_size;
    cs_real_t  *r2 = r + ssys->max_x1_size;

    r_norm = cs_dot_xx(ssys->x2_size, r2);
    cs_parall_sum(1, CS_DOUBLE, &r_norm);
    r_norm = sqrt(fabs(r_norm));

    memset(z2, 0, sizeof(cs_real_t)*ssys->x2_size);
    n_iter += cs_equation_solve_scalar_cell_system(ssys->x2_size,
                                                   sbp->schur_slesp,
                                                   sbp->schur_matrix,
                                                   r_norm,
                                                   sbp->schur_sles,
                                                   z2,
                                                   r2);

    if (sbp->schur_type == CS_PARAM_SCHUR_MASS_SCALED_DIAG_INVERSE ||
        sbp->schur_type == CS_PARAM_SCHUR_MASS_SCALED_LUMPED_INVERSE) {

      assert(sbp->mass22_diag != NULL); /* inverse of the scaled mass matrix */

#     pragma omp parallel for if (ssys->x2_size > CS_THR_MIN)
      for (cs_lnum_t i2 = 0; i2 < ssys->x2_size; i2++)
        z2[i2] = sbp->schur_scaling*z2[i2] + sbp->mass22_diag[i2]*r2[i2];

    }

  }

  return n_iter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply a block diagonal preconditioning: Compute z s.t. P_d z = r
 *
 * \param[in]      ssys    pointer to a cs_saddle_system_t structure
 * \param[in]      sbp     Block-preconditioner for the Saddle-point problem
 * \param[in]      r       rhs of the preconditioning system
 * \param[in, out] z       array to compute
 * \param[in, out] pc_wsp  work space related to the preconditioner
 *
 * \return the number of iteration performed for this step of preconditioning
 */
/*----------------------------------------------------------------------------*/

static int
_no_pc_apply(cs_saddle_system_t          *ssys,
             cs_saddle_block_precond_t   *sbp,
             cs_real_t                   *r,
             cs_real_t                   *z,
             cs_real_t                   *pc_wsp)
{
  CS_UNUSED(sbp);
  CS_UNUSED(pc_wsp);

  if (z == NULL)
    return 0;

  assert(r != NULL);

  /* First part */
  memcpy(z, r, sizeof(cs_real_t)*ssys->x1_size);

  /* Second part */
  memcpy(z + ssys->max_x1_size,
         r + ssys->max_x1_size, sizeof(cs_real_t)*ssys->x2_size);

  return 0;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize a cs_saddle_block_precond_t structure
 *
 * \param[in] block_type  type of block preconditioner
 * \param[in] schur_type  type of Schur approximation
 * \param[in] m11_slesp   pointer to the settings for the M11 block
 * \param[in] m11_sles    pointer to the cs_sles_t struct. for the M11 block
 *
 * \return a pointer to the new allocated strcuture
 */
/*----------------------------------------------------------------------------*/

cs_saddle_block_precond_t *
cs_saddle_block_precond_create(cs_param_precond_block_t    block_type,
                               cs_param_schur_approx_t     schur_type,
                               cs_param_sles_t            *m11_slesp,
                               cs_sles_t                  *m11_sles)
{
  cs_saddle_block_precond_t  *sbp = NULL;

  BFT_MALLOC(sbp, 1, cs_saddle_block_precond_t);

  sbp->block_type = block_type;
  sbp->schur_type = schur_type;

  /* Block 11 settings */
  sbp->m11_slesp = m11_slesp;
  sbp->m11_sles = m11_sles;

  /* According to the settings an approximation of the Schur complement may be
     requested*/
  sbp->schur_matrix = NULL;
  sbp->schur_slesp = NULL;
  sbp->schur_sles = NULL;
  sbp->schur_scaling = 1.;

  /* Native arrays for the Schur matrix */
  sbp->schur_diag = NULL;
  sbp->schur_xtra = NULL;

  /* Other approximations */
  sbp->mass22_diag = NULL;
  sbp->m11_inv_diag = NULL;

  return sbp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free a cs_saddle_block_precond_t structure
 *
 * \param[in, out] p_sbp  double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_block_precond_free(cs_saddle_block_precond_t  **p_sbp)
{
  if (p_sbp == NULL || *p_sbp == NULL)
    return;

  cs_saddle_block_precond_t  *sbp = *p_sbp;

  BFT_FREE(sbp->schur_diag);
  BFT_FREE(sbp->schur_xtra);

  BFT_FREE(sbp->mass22_diag);
  BFT_FREE(sbp->m11_inv_diag);

  /* Free the main pointer */
  BFT_FREE(sbp);
  *p_sbp = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply the MINRES algorithm to a saddle point problem (the system is
 *        stored in a hybrid way). Please refer to cs_saddle_system_t structure
 *        definition.
 *        The stride is equal to 1 for the matrix (db_size[3] = 1) and the
 *        vector.
 *
 * \param[in]      ssys    pointer to a cs_saddle_system_t structure
 * \param[in]      sbp     Block-preconditioner for the Saddle-point problem
 * \param[in, out] x1      array for the first part
 * \param[in, out] x2      array for the second part
 * \param[in, out] info    pointer to a cs_iter_algo_info_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_minres(cs_saddle_system_t          *ssys,
                 cs_saddle_block_precond_t   *sbp,
                 cs_real_t                   *x1,
                 cs_real_t                   *x2,
                 cs_iter_algo_info_t         *info)
{
  /* Workspace */
  const cs_lnum_t  ssys_size = ssys->max_x1_size + ssys->x2_size;
  const size_t  n_bytes = sizeof(cs_real_t)*ssys_size;
  cs_lnum_t  wsp_size = 7*ssys_size;
  cs_real_t  *wsp = NULL;

  BFT_MALLOC(wsp, wsp_size, cs_real_t);
  memset(wsp, 0, wsp_size*sizeof(cs_real_t));

  cs_real_t  *v     = wsp;
  cs_real_t  *vold  = wsp + ssys_size;
  cs_real_t  *w     = wsp + 2*ssys_size;
  cs_real_t  *wold  = wsp + 3*ssys_size;
  cs_real_t  *z     = wsp + 4*ssys_size;
  cs_real_t  *zold  = wsp + 5*ssys_size;
  cs_real_t  *mz    = wsp + 6*ssys_size;

  /* Set pointer for the block preconditioning */
  cs_lnum_t  pc_wsp_size = 0;
  cs_real_t  *pc_wsp = NULL;
  cs_saddle_pc_apply_t  *pc_apply = NULL;

  if (sbp == NULL) /* No preconditioning */
    pc_apply = _no_pc_apply;

  else {

    switch (sbp->block_type) {

    case CS_PARAM_PRECOND_BLOCK_NONE:
      /* No preconditioning */
      pc_apply = _no_pc_apply;
      break;
    case CS_PARAM_PRECOND_BLOCK_DIAG:
      switch (sbp->schur_type) {

      case CS_PARAM_SCHUR_DIAG_INVERSE:
      case CS_PARAM_SCHUR_IDENTITY:
      case CS_PARAM_SCHUR_LUMPED_INVERSE:
      case CS_PARAM_SCHUR_MASS_SCALED:
      case CS_PARAM_SCHUR_MASS_SCALED_DIAG_INVERSE:
      case CS_PARAM_SCHUR_MASS_SCALED_LUMPED_INVERSE:
        pc_apply = _diag_schur_pc_apply;
        break;
      case CS_PARAM_SCHUR_ELMAN:
        pc_apply = _elman_schur_pc_apply;
        pc_wsp_size = 2*ssys->max_x1_size;
        BFT_MALLOC(pc_wsp, pc_wsp_size, cs_real_t);
        break;

      default:
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Invalid Schur approximation.", __func__);
      }
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                "%s: Invalid block preconditioner.", __func__);
    }

  }

  /* --- ALGO BEGIN --- */

  /* Compute the first residual: v = b - M.x */
  _compute_residual_3(ssys, x1, x2, v);

  /* Apply preconditioning: M.z = v */
  info->n_inner_iter += pc_apply(ssys, sbp, v, z, pc_wsp);

  info->res0 = _norm(ssys, v); /* ||v|| */
  info->res = info->res0;

  /* dp = eta = <r, z>; beta = sqrt(dp) */
  double  dp = _dot_product(ssys, v, z);
  double  beta = sqrt(fabs(dp));
  double  eta = beta;

  /* Initialization */
  double  betaold = 1;
  double  c=1.0, cold=1.0, s=0.0, sold=0.0;

  while (info->cvg == CS_SLES_ITERATING) {

    /* z = z * ibeta; */
    assert(fabs(beta) > 0.);
    const double  ibeta = 1./beta;
    _scalar_scaling(ssys, ibeta, z);

    /* Compute the matrix-vector product M.z = mz */
    _matvec_product(ssys, z, mz);

    /* alpha = <z, mz> */
    const double  alpha =  _dot_product(ssys, z, mz);
    const double  alpha_ibeta = alpha * ibeta;
    const double  beta_ibetaold = beta/betaold;

    /* v(k+1) = mz(k) - alpha*v(k) - beta v(k-1) */
#   pragma omp parallel for if (ssys->x1_size > CS_THR_MIN)
    for (cs_lnum_t i1 = 0; i1 < ssys->x1_size; i1++) {
      const cs_real_t  _v = v[i1], _vold = vold[i1];
      v[i1] = mz[i1] - alpha_ibeta*_v - beta_ibetaold*_vold;
      vold[i1] = _v;
    }

    cs_real_t  *v2 = v + ssys->max_x1_size;
    cs_real_t  *v2old = vold + ssys->max_x1_size;
    const cs_real_t  *mz2 = mz + ssys->max_x1_size;

#   pragma omp parallel for if (ssys->x2_size > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < ssys->x2_size; i2++) {
      const cs_real_t  _v = v2[i2], _vold = v2old[i2];
      v2[i2] = mz2[i2] - alpha_ibeta*_v - beta_ibetaold*_vold;
      v2old[i2] = _v;
    }

    /* Apply preconditionning: M.z(k+1) = v(k+1) */
    memcpy(zold, z, n_bytes);
    info->n_inner_iter += pc_apply(ssys, sbp, v, z, pc_wsp);

    /* New value for beta: beta = sqrt(<r, z>) */
    betaold = beta;
    beta = sqrt(fabs(_dot_product(ssys, v, z)));

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
#   pragma omp parallel for if (ssys->x1_size > CS_THR_MIN)
    for (cs_lnum_t i1 = 0; i1 < ssys->x1_size; i1++) {
      const cs_real_t  _w = w[i1], _wold = wold[i1];
      w[i1] = irho1 * (zold[i1] - rho2*_w - rho3*_wold);
      wold[i1] = _w;
    }

    cs_real_t  *w2 = w + ssys->max_x1_size;
    cs_real_t  *w2old = wold + ssys->max_x1_size;
    const cs_real_t  *z2old = zold + ssys->max_x1_size;

#   pragma omp parallel for if (ssys->x2_size > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < ssys->x2_size; i2++) {
      const cs_real_t  _w = w2[i2], _wold = w2old[i2];
      w2[i2] = irho1 * (z2old[i2] - rho2*_w - rho3*_wold);
      w2old[i2] = _w;
    }

    /* Update the solution vector */
    /* x1(k+1) = x1(k) + c*eta*w(k+1) */
    const double  ceta = c*eta;

#   pragma omp parallel for if (ssys->x1_size > CS_THR_MIN)
    for (cs_lnum_t i1 = 0; i1 < ssys->x1_size; i1++)
      x1[i1] = x1[i1] + ceta*w[i1];

#   pragma omp parallel for if (ssys->x2_size > CS_THR_MIN)
    for (cs_lnum_t i2 = 0; i2 < ssys->x2_size; i2++)
      x2[i2] = x2[i2] + ceta*w2[i2];

    /* Compute the current residual */
    info->res *= fabs(s);

    /* Last updates */
    eta = -s*eta;

    /* Check the convergence criteria */
    _cvg_test(info);

  } /* main loop */

  /* --- ALGO END --- */

#if 1
  /* Compute the real residual norm at exit */
  _compute_residual_3(ssys, x1, x2, v);
  beta = _norm(ssys, v); /* ||v|| */
  printf(" >> Residual norm at exit= %6.4e\n", beta);
  cs_log_printf(CS_LOG_DEFAULT, " >> Residual norm at exit= %6.4e\n", beta);
#endif

  /* Free temporary workspace */
  BFT_FREE(wsp);
  BFT_FREE(pc_wsp);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform a matrix-vector multiplication in case of scatter-view array
 *        as input parameter.  Thus, one performs a scatter --> gather (before
 *        the multiplication) and a gather --> scatter operation after the
 *        multiplication.  One assumes that matvec is allocated to the right
 *        size. No check is done.
 *
 *        The stride is equal to 1 for the matrix (db_size[3] = 1) and the
 *        vector
 *
 * \param[in]      rset      pointer to a cs_range_set_t structure
 * \param[in]      mat       matrix
 * \param[in, out] vec       vector
 * \param[in, out] matvec    resulting vector for the matrix-vector product
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_vector_multiply_gs_allocated(const cs_range_set_t      *rset,
                                       const cs_matrix_t         *mat,
                                       cs_real_t                 *vec,
                                       cs_real_t                 *matvec)
{
  if (mat == NULL || vec == NULL)
    return;

  /* Handle the input array
   * n_rows = n_gather_elts <= n_scatter_elts = n_dofs (mesh view) <= n_cols
   */

  /* scatter view to gather view for the input vector */
  cs_range_set_gather(rset,
                      CS_REAL_TYPE,  /* type */
                      1,             /* stride */
                      vec,           /* in:  size=n_sles_scatter_elts */
                      vec);          /* out: size=n_sles_gather_elts */

  cs_matrix_vector_multiply(CS_HALO_ROTATION_IGNORE, mat, vec, matvec);

  /* gather view to scatter view (i.e. algebraic to mesh view) */
  if (rset != NULL) {
    cs_range_set_scatter(rset,
                         CS_REAL_TYPE, 1, /* type and stride */
                         vec, vec);
    cs_range_set_scatter(rset,
                         CS_REAL_TYPE, 1, /* type and stride */
                         matvec, matvec);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform a matrix-vector multiplication in case of scatter-view array
 *        as input parameter.  Thus, one performs a scatter --> gather (before
 *        the multiplication) and a gather --> scatter operation after the
 *        multiplication.  The output parameter matvec is not allocated. A
 *        check on the size is done for the input array.
 *
 *        The stride is equal to 1 for the matrix (db_size[3] = 1) and the
 *        vector
 *
 * \param[in]      rset      pointer to a cs_range_set_t structure
 * \param[in]      mat       matrix
 * \param[in]      vec_len   size of vec
 * \param[in, out] vec       vector of real numbers
 * \param[out]     matvec    resulting vector for the matrix-vector product
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_vector_multiply_gs(const cs_range_set_t      *rset,
                             const cs_matrix_t         *mat,
                             cs_lnum_t                  vec_len,
                             cs_real_t                 *vec,
                             cs_real_t                **p_matvec)
{
  if (mat == NULL || vec == NULL)
    return;

  const cs_lnum_t  n_cols = cs_matrix_get_n_columns(mat);

  /* Handle the input array
   * n_rows = n_gather_elts <= n_scatter_elts = n_dofs (mesh view) <= n_cols
   */
  cs_real_t  *vecx = NULL;
  if (n_cols > vec_len) {
    BFT_MALLOC(vecx, n_cols, cs_real_t);
    memcpy(vecx, vec, sizeof(cs_real_t)*vec_len);
  }
  else
    vecx = vec;

  /* scatter view to gather view */
  if (rset != NULL)
    cs_range_set_gather(rset,
                        CS_REAL_TYPE,  /* type */
                        1,             /* stride */
                        vecx,          /* in:  size=n_sles_scatter_elts */
                        vecx);         /* out: size=n_sles_gather_elts */

  /* Handle the output array */
  cs_real_t  *matvec = NULL;
  BFT_MALLOC(matvec, n_cols, cs_real_t);

  cs_matrix_vector_multiply(CS_HALO_ROTATION_IGNORE, mat, vecx, matvec);

  /* gather to scatter view (i.e. algebraic to mesh view) */
  if (rset != NULL) {
    cs_range_set_scatter(rset,
                         CS_REAL_TYPE, 1, /* type and stride */
                         vecx, vec);
    cs_range_set_scatter(rset,
                         CS_REAL_TYPE, 1, /* type and stride */
                         matvec, matvec);
  }

  /* Free allocated memory if needed */
  if (vecx != vec) {
    assert(n_cols > vec_len);
    BFT_FREE(vecx);
  }

  /* return the resulting array */
  *p_matvec = matvec;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
