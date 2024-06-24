/*============================================================================
 * Set of functions to manage high-level iterative algorithms
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

#include <float.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_array.h"
#include "cs_log.h"
#include "cs_parall.h"
#include "cs_sdm.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_iter_algo.h"

/*----------------------------------------------------------------------------*/
BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_iter_algo.c
 *
 * \brief Set of functions to manage high-level iterative algorithms such as
 *        Uzawa, Golub-Kahan Bi-orthogonalization, block preconditioner or
 *        Picard and Anderson algorithms which may incorporate inner iterative
 *        solvers
 */

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

#define CS_ITER_ALGO_DBG      0

/*============================================================================
 * Private variables
 *============================================================================*/

/*============================================================================
 * Inline static private functions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reset the context of a cs_iter_algo_t structure when this context is
 *        associated to the "default" type
 *
 * \param[in, out] c       pointer to a cs_iter_algo_default_t structure
 */
/*----------------------------------------------------------------------------*/

static inline void
_reset_default(cs_iter_algo_default_t    *c)
{
  if (c == nullptr)
    return;

  c->cvg_status = CS_SLES_ITERATING;

  c->res0 = cs_math_big_r;
  c->prev_res = cs_math_big_r;
  c->res = cs_math_big_r;

  c->n_algo_iter = 0;
  c->n_inner_iter = 0;
  c->last_inner_iter = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reset the context of a cs_iter_algo_t structure when this context is
 *        associated to an Anderson acceleration algorithm
 *
 * \param[in, out] c       pointer to a cs_iter_algo_aac_t structure
 */
/*----------------------------------------------------------------------------*/

static inline void
_reset_anderson(cs_iter_algo_aac_t    *c)
{
  if (c == nullptr)
    return;

  c->cvg_status = CS_SLES_ITERATING;

  c->res0 = cs_math_big_r;
  c->prev_res = cs_math_big_r;
  c->res = cs_math_big_r;

  c->n_algo_iter = 0;
  c->n_inner_iter = 0;
  c->last_inner_iter = 0;

  c->n_dir = 0;
}

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute an estimation of the condition number of the matrix R based
 *         on the diagonal entries
 *
 * \param[in]  n_rows  number of rows to consider
 * \param[in]  R       pointer to the matrix R
 *
 * \return the estimation of the condition number
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_condition_number(int               n_rows,
                  const cs_sdm_t   *R)
{
  int  n_cols = R->n_cols;
  cs_real_t  vmax = fabs(R->val[0]);
  cs_real_t  vmin = fabs(R->val[0]);

  for (int i = 1; i < n_rows; i++){

    vmax = fmax(vmax, fabs(R->val[i*(1+n_cols)]));
    vmin = fmin(vmin, fabs(R->val[i*(1+n_cols)]));

  }

  return vmax/fmax(vmin, cs_math_epzero);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Erase the first row of the dg set of arrays and then shift each
 *         superior row below
 *
 * \param[in, out] c        pointer to the Anderson algo. structure
 */
/*----------------------------------------------------------------------------*/

static void
_anderson_shift_dg(cs_iter_algo_aac_t   *c)
{
  assert(c->dg != nullptr);
  const size_t  dir_size = sizeof(cs_real_t)*c->n_elts;

  for (int i = 0; i < c->n_dir - 1; i++) {

    cs_real_t *dg_i = c->dg + i*c->n_elts;
    cs_real_t *dg_ip1 = c->dg + (i+1)*c->n_elts;

    memcpy(dg_i, dg_ip1, dir_size); /* dg_i <= dg_(i+1) */

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the Q.R factorization
 *
 * \param[in]      m        number of columns to consider
 * \param[in]      n_elts   size of one column
 * \param[in, out] Q        pointer to the set of Q arrays
 * \param[in, out] R        pointer to the matrix R
 */
/*----------------------------------------------------------------------------*/

static void
_qrdelete(int          m,
          cs_lnum_t    n_elts,
          cs_real_t   *Q,
          cs_sdm_t    *R)
{
  assert(R->n_rows == R->n_cols); /* square matrix */

  const int  n_cols = R->n_cols;

  for (int i = 0; i < m-1; i++) {

    cs_real_t  *R_i = R->val + i*n_cols;
    cs_real_t  *R_ip1 = R->val + (i+1)*n_cols;

    double  tempR = sqrt(R_i[i+1]*R_i[i+1] + R_ip1[i+1]*R_ip1[i+1]);
    const double  c = R_i[i+1]/tempR;
    const double  s = R_ip1[i+1]/tempR;

     R_i[i+1] = tempR;
     R_ip1[i+1] = 0.0;

     if (i < m-2) {

       for (int j = i+2; j < m; j++) {
         tempR = c*R_i[j] + s*R_ip1[j];
         R_ip1[j] = -s*R_i[j] + c*R_ip1[j];
         R_i[j] = tempR;
       }

     }

     cs_real_t *Q_i = Q + i*n_elts;
     cs_real_t *Q_ip1 = Q + (i+1)*n_elts;

     for (cs_lnum_t l = 0; l < n_elts; l++) {

       double  tempQ = c*Q_i[l] + s*Q_ip1[l];

       Q_ip1[l] = -s*Q_i[l] + c*Q_ip1[l];
       Q_i[l] = tempQ;

     }

  } /* Loop on m */

  /* The Q and R matrices are not resized.
   * A block version should be used */

  for (int i = 0; i < m-1; i++) {

    cs_real_t  *R_i = R->val + i*n_cols;
    for (int j = i; j < m-1; j++)
      R_i[j] = R_i[j+1];

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the coefficient used in the linear combination for the
 *        Anderson acceleration
 *
 * \param[in, out] c         pointer to the Anderson algo. structure
 * \param[in]      dotprod   function to compute a dot product
 */
/*----------------------------------------------------------------------------*/

static void
_anderson_compute_gamma(cs_iter_algo_aac_t          *c,
                        cs_cdo_blas_dotprod_t       *dotprod)
{
  int  n_cols = c->R->n_cols;

  for (int i = c->n_dir-1; i >= 0; i--) {

    /* Solve R gamma = (fcur, Q_i) where fcur is equal to fold since one saves
       fcur into fold at the end of the step 1 of the Anderson algorithm */

    double  s = dotprod(c->fold, c->Q + i*c->n_elts);

    for (int j = i+1; j < c->n_dir; j++)
      s -= c->R->val[i*n_cols+j] * c->gamma[j];

    c->gamma[i] = s/c->R->val[i*(n_cols+1)];

  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_ITER_ALGO_DBG > 0
  cs_log_printf(CS_LOG_DEFAULT, " %s:l%d >> Gamma coefficients:",
                __func__, __LINE__);
  for (int i = 0; i < c->n_dir; i++)
    cs_log_printf(CS_LOG_DEFAULT, " (%d:=%5.3e)", i, c->gamma[i]);
  cs_log_printf(CS_LOG_DEFAULT, "\n");
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the coefficient used in the linear combination for the
 *        Anderson acceleration
 *
 * \param[in, out] c    pointer to the structure managing the Anderson algo.
 * \param[in, out] x    array of values to update
 */
/*----------------------------------------------------------------------------*/

static void
_anderson_damping(cs_iter_algo_aac_t    *c,
                  cs_real_t             *x)
{
  const double  omb = 1.0 - c->param.beta;
  const int  m_max = c->param.n_max_dir;
  const cs_real_t  *Rval = c->R->val;

  for (int j = 0; j < c->n_dir; j++) {

    double  R_gamma = 0.;
    for (int k = j; k < c->n_dir; k++)  /* R is upper diag */
      R_gamma += Rval[j*m_max+k] * c->gamma[k];

    /* x = x - (1-beta)*(fold - Q*R*gamma) */

    const cs_real_t  *Qj = c->Q + j*c->n_elts;  /* get row j */

#   pragma omp parallel if (c->n_elts > CS_THR_MIN)
    for (cs_lnum_t l = 0; l < c->n_elts; l++)
      x[l] += omb * (Qj[l] * R_gamma - c->fold[l]);

  } /* Loop on directions */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate arrays needed by the "Anderson acceleration" algorithm
 *
 * \param[in, out] c    pointer to the context structure for the Anderson algo.
 */
/*----------------------------------------------------------------------------*/

static void
_allocate_anderson_arrays(cs_iter_algo_aac_t     *c)
{
  if (c == nullptr)
    return;

  const int  n_max_dir = c->param.n_max_dir;

  BFT_MALLOC(c->fold, c->n_elts, cs_real_t);
  cs_array_real_fill_zero(c->n_elts, c->fold);

  BFT_MALLOC(c->df, c->n_elts, cs_real_t);
  cs_array_real_fill_zero(c->n_elts, c->df);

  BFT_MALLOC(c->gold, c->n_elts, cs_real_t);
  cs_array_real_fill_zero(c->n_elts, c->gold);

  BFT_MALLOC(c->dg, c->n_elts*n_max_dir, cs_real_t);
  cs_array_real_fill_zero(c->n_elts*n_max_dir, c->dg);

  BFT_MALLOC(c->gamma, n_max_dir, cs_real_t);
  memset(c->gamma, 0, sizeof(cs_real_t)*n_max_dir);

  BFT_MALLOC(c->Q, c->n_elts*n_max_dir, cs_real_t);
  cs_array_real_fill_zero(c->n_elts*n_max_dir, c->Q);

  c->R = cs_sdm_square_create(n_max_dir);
  cs_sdm_square_init(n_max_dir, c->R);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize by default a new cs_iter_algo_t structure
 *
 * \param[in] type    type of iterative algorithm
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_iter_algo_t *
cs_iter_algo_create(cs_iter_algo_type_t    type)
{
  cs_iter_algo_t *algo = nullptr;

  BFT_MALLOC(algo, 1, cs_iter_algo_t);

  algo->type = type;
  algo->verbosity = 0;
  algo->cvg_param.atol = 1e-15;
  algo->cvg_param.rtol = 1e-6;
  algo->cvg_param.dtol = 1e3;
  algo->cvg_param.n_max_iter = 2500;

  if (type & CS_ITER_ALGO_DEFAULT) {

    cs_iter_algo_default_t *c = nullptr;

    BFT_MALLOC(c, 1, cs_iter_algo_default_t);
    _reset_default(c);
    c->normalization = 1.0;

    algo->context = c;

  }
  else if (type & CS_ITER_ALGO_ANDERSON) {

    cs_iter_algo_aac_t *c = nullptr;

    BFT_MALLOC(c, 1, cs_iter_algo_aac_t);

    /* Parameters (default settings) */

    c->param.n_max_dir = 5;
    c->param.starting_iter = 3;
    c->param.max_cond = 500;
    c->param.beta = 0;
    c->param.dp_type = CS_PARAM_DOTPROD_EUCLIDEAN;

    /* Other variables */

    c->n_elts = 0;    /* Should be set latter */
    c->n_dir = 0;

    /* Work arrays and structures */

    c->fold  = nullptr;
    c->df    = nullptr;
    c->gold  = nullptr;
    c->dg    = nullptr;
    c->gamma = nullptr;
    c->Q     = nullptr;
    c->R     = nullptr;

    _reset_anderson(c);
    c->normalization = 1.0;

    algo->context = c;

  }
  else
    algo->context = nullptr;

  return algo;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a new cs_iter_algo_t structure with the given settings
 *
 * \param[in] type       type of iterative algorithm
 * \param[in] verbosity  level of information to print
 * \param[in] cvg_param  set of parameters driving the convergence of the
 *                       iterative algorithm
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_iter_algo_t *
cs_iter_algo_create_with_settings(cs_iter_algo_type_t     type,
                                  int                     verbosity,
                                  cs_param_convergence_t  cvg_param)
{
  cs_iter_algo_t  *algo = cs_iter_algo_create(type);

  algo->verbosity = verbosity;
  algo->cvg_param.atol = cvg_param.atol;
  algo->cvg_param.rtol = cvg_param.rtol;
  algo->cvg_param.dtol = cvg_param.dtol;
  algo->cvg_param.n_max_iter = cvg_param.n_max_iter;

  return algo;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free a cs_iter_algo_t structure
 *
 * \param[in, out] p_algo     double pointer on the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_free(cs_iter_algo_t   **p_algo)
{
  if (p_algo == nullptr)
    return;

  cs_iter_algo_t  *algo = *p_algo;
  if (algo == nullptr)
    return;

  if (algo->type & CS_ITER_ALGO_DEFAULT) {

    cs_iter_algo_default_t *c
      = static_cast<cs_iter_algo_default_t *>(algo->context);

    BFT_FREE(c);

  }
  else if (algo->type & CS_ITER_ALGO_ANDERSON) {

    cs_iter_algo_aac_t *c = static_cast<cs_iter_algo_aac_t *>(algo->context);

    cs_iter_algo_release_anderson_arrays(c);

    BFT_FREE(c);

  }

  BFT_FREE(algo);
  *p_algo = nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reset a cs_iter_algo_t structure
 *
 * \param[in, out] algo           pointer to a cs_iter_algo_t
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_reset(cs_iter_algo_t      *algo)
{
  if (algo == nullptr)
    return;

  if (algo->type & CS_ITER_ALGO_DEFAULT) {

    cs_iter_algo_default_t *c
      = static_cast<cs_iter_algo_default_t *>(algo->context);
    _reset_default(c);

  }
  else if (algo->type & CS_ITER_ALGO_ANDERSON) {

    cs_iter_algo_aac_t *c = static_cast<cs_iter_algo_aac_t *>(algo->context);
    _reset_anderson(c);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the members (arrays and matrix) associated to the context
 *        structure of an Anderson acceleration
 *
 * \param[in, out] c   pointer to an Anderson context structure
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_release_anderson_arrays(cs_iter_algo_aac_t      *c)
{
  if (c == nullptr)
    return;

  BFT_FREE(c->fold);
  BFT_FREE(c->df);
  BFT_FREE(c->gold);
  BFT_FREE(c->dg);
  BFT_FREE(c->gamma);
  BFT_FREE(c->Q);

  c->R = cs_sdm_free(c->R);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the verbosity of the given iterative algorithm
 *
 * \param[in, out] algo       pointer to the structure to update
 * \param[in]      verbosity  level of information to print
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_set_verbosity(cs_iter_algo_t     *algo,
                           int                 verbosity)
{
  if (algo == nullptr)
    return;

  algo->verbosity = verbosity;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the criteria related to the convergence of the given iterative
 *        algorithm
 *
 * \param[in, out] algo       pointer to the structure to update
 * \param[in]      cvg_param  set of parameters driving the convergence of the
 *                            iterative algorithm
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_set_cvg_param(cs_iter_algo_t         *algo,
                           cs_param_convergence_t  cvg_param)
{
  if (algo == nullptr)
    return;

  algo->cvg_param.atol = cvg_param.atol;
  algo->cvg_param.rtol = cvg_param.rtol;
  algo->cvg_param.dtol = cvg_param.dtol;
  algo->cvg_param.n_max_iter = cvg_param.n_max_iter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the final tolerance used to check the convergence of the algorithm
 *        This tolerance should take into account rtol and atol for instance
 *
 * \param[in, out] algo    pointer to the structure to update
 * \param[in]      tol     tolerance to apply
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_set_tolerance(cs_iter_algo_t   *algo,
                           double            tol)
{
  if (algo == nullptr)
    return;

  if (algo->type & CS_ITER_ALGO_DEFAULT) {

    cs_iter_algo_default_t *c
      = static_cast<cs_iter_algo_default_t *>(algo->context);
    c->tol = tol;

  }
  else if (algo->type & CS_ITER_ALGO_ANDERSON) {

    cs_iter_algo_aac_t *c = static_cast<cs_iter_algo_aac_t *>(algo->context);
    c->tol = tol;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the initial residual used to detect a divergence
 *
 * \param[in, out] algo     pointer to the structure to update
 * \param[in]      value    value of the initial residual
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_set_initial_residual(cs_iter_algo_t   *algo,
                                  double            value)
{
  if (algo == nullptr)
    return;

  if (algo->type & CS_ITER_ALGO_DEFAULT) {

    cs_iter_algo_default_t *c
      = static_cast<cs_iter_algo_default_t *>(algo->context);
    c->res0 = value;

  }
  else if (algo->type & CS_ITER_ALGO_ANDERSON) {

    cs_iter_algo_aac_t *c = static_cast<cs_iter_algo_aac_t *>(algo->context);
    c->res0 = value;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the normalization to apply when checking the convergence of the
 *        algorithm.
 *
 * \param[in, out] algo     pointer to the structure to update
 * \param[in]      value    normalization to apply
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_set_normalization(cs_iter_algo_t   *algo,
                               double            value)
{
  if (algo == nullptr)
    return;

  if (algo->type & CS_ITER_ALGO_DEFAULT) {

    cs_iter_algo_default_t *c
      = static_cast<cs_iter_algo_default_t *>(algo->context);
    c->normalization = value;

  }
  else if (algo->type & CS_ITER_ALGO_ANDERSON) {

    cs_iter_algo_aac_t *c = static_cast<cs_iter_algo_aac_t *>(algo->context);
    c->normalization = value;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the convergence status of the given structure
 *
 * \param[in, out] algo        pointer to the structure to update
 * \param[in]      cvg_status  status to set
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_set_cvg_status(cs_iter_algo_t               *algo,
                            cs_sles_convergence_state_t   cvg_status)
{
  if (algo == nullptr)
    return;

  if (algo->type & CS_ITER_ALGO_DEFAULT) {

    cs_iter_algo_default_t *c
      = static_cast<cs_iter_algo_default_t *>(algo->context);
    c->cvg_status = cvg_status;

  }
  else if (algo->type & CS_ITER_ALGO_ANDERSON) {

    cs_iter_algo_aac_t *c = static_cast<cs_iter_algo_aac_t *>(algo->context);
    c->cvg_status = cvg_status;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a new cs_iter_algo_aac_t structure for Anderson acceleration
 *
 * \param[in, out] algo        pointer to the structure to update
 * \param[in]      aac_param   set of parameters for the Anderson acceleration
 * \param[in]      n_elts      number of elements by direction
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_set_anderson_param(cs_iter_algo_t             *algo,
                                cs_iter_algo_param_aac_t    aac_param,
                                cs_lnum_t                   n_elts)
{
  if (algo == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: Empty structure.", __func__);
  if (cs_flag_test(algo->type, CS_ITER_ALGO_ANDERSON) == false)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of iterative algorithm.", __func__);

  cs_iter_algo_aac_t *c = static_cast<cs_iter_algo_aac_t *>(algo->context);
  assert(c != nullptr);

  c->param.n_max_dir = aac_param.n_max_dir;
  c->param.starting_iter = aac_param.starting_iter;
  c->param.max_cond = aac_param.max_cond;
  c->param.beta = aac_param.beta;
  c->param.dp_type = aac_param.dp_type;

  c->n_elts = n_elts;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the current number of iterations done
 *
 * \param[in, out] algo       pointer to the structure to examine
 *
 * \return the number of iterations done
 */
/*----------------------------------------------------------------------------*/

int
cs_iter_algo_get_n_iter(const cs_iter_algo_t        *algo)
{
  if (algo == nullptr)
    return 0;

  if (algo->type & CS_ITER_ALGO_DEFAULT) {

    cs_iter_algo_default_t *c
      = static_cast<cs_iter_algo_default_t *>(algo->context);
    return c->n_algo_iter;

  }
  else if (algo->type & CS_ITER_ALGO_ANDERSON) {

    cs_iter_algo_aac_t *c = static_cast<cs_iter_algo_aac_t *>(algo->context);
    return c->n_algo_iter;

  }
  else
    return 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the cumulated number of inner iterations done
 *
 * \param[in, out] algo       pointer to the structure to examine
 *
 * \return the number of iterations done
 */
/*----------------------------------------------------------------------------*/

int
cs_iter_algo_get_n_inner_iter(const cs_iter_algo_t        *algo)
{
  if (algo == nullptr)
    return 0;
  if (cs_flag_test(algo->type, CS_ITER_ALGO_TWO_LEVEL) == false)
    return 0;

  if (algo->type & CS_ITER_ALGO_DEFAULT) {

    cs_iter_algo_default_t *c
      = static_cast<cs_iter_algo_default_t *>(algo->context);

    return c->n_inner_iter;

  }
  else if (algo->type & CS_ITER_ALGO_ANDERSON) {

    cs_iter_algo_aac_t *c = static_cast<cs_iter_algo_aac_t *>(algo->context);

    return c->n_inner_iter;

  }
  else
    return 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the last computed residual
 *
 * \param[in, out] algo       pointer to the structure to examine
 *
 * \return the number of iterations done
 */
/*----------------------------------------------------------------------------*/

double
cs_iter_algo_get_residual(const cs_iter_algo_t        *algo)
{
  if (algo == nullptr)
    return cs_math_big_r;

  if (algo->type & CS_ITER_ALGO_DEFAULT) {

    cs_iter_algo_default_t *c
      = static_cast<cs_iter_algo_default_t *>(algo->context);
    return c->res;

  }
  else if (algo->type & CS_ITER_ALGO_ANDERSON) {

    cs_iter_algo_aac_t *c = static_cast<cs_iter_algo_aac_t *>(algo->context);
    return c->res;

  }
  else
    return cs_math_big_r;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the last convergence state
 *
 * \param[in, out] algo       pointer to the structure to examine
 *
 * \return the convergence status
 */
/*----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_iter_algo_get_cvg_status(const cs_iter_algo_t        *algo)
{
  if (algo == nullptr)
    return CS_SLES_ITERATING;

  if (algo->type & CS_ITER_ALGO_DEFAULT) {

    cs_iter_algo_default_t *c
      = static_cast<cs_iter_algo_default_t *>(algo->context);
    return c->cvg_status;

  }
  else if (algo->type & CS_ITER_ALGO_ANDERSON) {

    cs_iter_algo_aac_t *c = static_cast<cs_iter_algo_aac_t *>(algo->context);
    return c->cvg_status;

  }
  else
    return CS_SLES_ITERATING;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the normalization to apply when computing the tolerance threshold
 *
 * \param[in] algo   pointer to a cs_iter_algo_t structure
 */
/*----------------------------------------------------------------------------*/

double
cs_iter_algo_get_normalization(const cs_iter_algo_t    *algo)
{
  if (algo == nullptr)
    return 1.0;

  if (algo->type & CS_ITER_ALGO_DEFAULT) {

    cs_iter_algo_default_t *c
      = static_cast<cs_iter_algo_default_t *>(algo->context);
    return c->normalization;

  }
  else if (algo->type & CS_ITER_ALGO_ANDERSON) {

    cs_iter_algo_aac_t *c = static_cast<cs_iter_algo_aac_t *>(algo->context);
    return c->normalization;

  }
  else
    return 1.0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the set of parameters for an Anderson algorithm
 *
 * \param[in, out] algo      pointer to a cs_iter_algo_t structure
 *
 * \return a cs_iter_algo_param_aac_t structure
 */
/*----------------------------------------------------------------------------*/

cs_iter_algo_param_aac_t
cs_iter_algo_get_anderson_param(cs_iter_algo_t         *algo)
{
  if (algo == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: Empty structure.", __func__);
  if (cs_flag_test(algo->type, CS_ITER_ALGO_ANDERSON) == false)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid type of iterative algorithm.", __func__);

  cs_iter_algo_aac_t *c = static_cast<cs_iter_algo_aac_t *>(algo->context);
  assert(c != nullptr);

  return c->param;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply one more iteration of the Anderson acceleration
 *
 * \param[in, out] algo           pointer to a cs_iter_algo_t structure
 * \param[in, out] cur_iterate    current iterate
 * \param[in]      pre_iterate    previous iterate
 * \param[in]      dotprod        function to compute a dot product
 * \param[in]      sqnorm         function to compute a square norm
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_update_anderson(cs_iter_algo_t               *algo,
                             cs_real_t                    *cur_iterate,
                             const cs_real_t              *pre_iterate,
                             cs_cdo_blas_dotprod_t        *dotprod,
                             cs_cdo_blas_square_norm_t    *sqnorm)
{
  if (algo == nullptr)
    return;

  if (cs_flag_test(algo->type, CS_ITER_ALGO_ANDERSON) == false)
    return; /* Nothing to do. This is not an Anderson acceleration algorithm */

  cs_iter_algo_aac_t  *c = (cs_iter_algo_aac_t *)algo->context;

  if (c == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Anderson acceleration context not allocated.\n", __func__);

  /* Check if anderson has to begin */

  const int  shifted_iter = c->n_algo_iter - c->param.starting_iter;
  if (shifted_iter < 0)
    return;

  if (shifted_iter == 0 && c->fold == nullptr)
    _allocate_anderson_arrays(c);

  /*
   * Notation details:
   * ----------------
   * gcur = Gfunc(pre_iterate) = cur_iterate
   * fcur = gcur - pre_iterate = cur_iterate - pre_iterate
   *
   * dg = gcur - gold
   * df = fcur - fold
   */

  const int  m_max = c->param.n_max_dir;
  const cs_real_t  *gcur = cur_iterate;

  /* Step 1: Update arrays useful for the Q.R factorization (df and dg)
   * ------- */

  if (shifted_iter == 0) { /* The first iteration is simpler */

    assert(c->n_dir == 0);

    /* Set fold and gold */

#   pragma omp parallel if (c->n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < c->n_elts; i++) {
      c->fold[i] = gcur[i] - pre_iterate[i];
      c->gold[i] = gcur[i];
    }

    return;  /* Nothing else to do. Algorithm has been initialized. This first
              * step is similar to a Picard iteration */

  }
  else {

    assert(shifted_iter > 0);

    /* Shift dg (set of rows) to the right position (two cases) */

    cs_real_t  *dg = c->dg + c->n_dir*c->n_elts;

    if (c->n_dir == m_max) {

      /* Erase the first row to retrieve some space and then
       * shift each superior row below */

      _anderson_shift_dg(c);
      dg = c->dg + (c->n_dir-1)*c->n_elts;

    }

    /* Set dg, df, fold and gold */

#   pragma omp parallel if (c->n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < c->n_elts; i++) {

      cs_real_t  fcur = gcur[i] - pre_iterate[i];

      dg[i] = gcur[i] - c->gold[i];
      c->df[i] = fcur - c->fold[i];
      c->fold[i] = fcur;
      c->gold[i] = gcur[i];

    }

    c->n_dir++;

  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_ITER_ALGO_DBG > 1
  cs_log_printf(CS_LOG_DEFAULT, " %s:l%d >> Case n_dir=%d\n",
                __func__, __LINE__, c->n_dir);
#endif

  /* Step 2: Update the Q.R factorization
   * ------- */

  cs_real_t  *Rval = c->R->val;

  if (c->n_dir == 1) {

    const double  df_norm = sqrt(sqnorm(c->df));
    const cs_real_t  coef = 1.0/df_norm;

    Rval[0] = df_norm; /* R(0,0) = |df|_L2 */

#   pragma omp parallel if (c->n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < c->n_elts; i++)
      c->Q[i] = c->df[i]*coef; /* Q(0) = df/|df|_L2 */

  }
  else {

    if (c->n_dir > m_max) { /* Remove the first column and last line in the
                                Q.R factorization */

      _qrdelete(m_max, c->n_elts, c->Q, c->R);
      c->n_dir--;

    }

    for (int j = 0; j < c->n_dir-1; j++) {

      cs_real_t *Qj = c->Q + j*c->n_elts; /* get the row j */

      const double  prod = dotprod(c->df, Qj);

      Rval[j*m_max + (c->n_dir-1)] = prod;  /* R(j, n_dir) = Qj*df */

#     pragma omp parallel if (c->n_elts > CS_THR_MIN)
      for (cs_lnum_t l = 0; l < c->n_elts; l++)
        c->df[l] -= prod*Qj[l];  /* update df = df - R(j, n_dir)*Qj */

    }

    const double  df_norm = sqrt(sqnorm(c->df));
    const cs_real_t  coef = 1.0/df_norm;

    /* R(n_dir,n_dir) = |df|_L2 */

    Rval[(c->n_dir-1)*m_max+(c->n_dir-1)] = df_norm;

    /* Set the row n_dir-1 of Q */

    cs_real_t *q_n_dir = c->Q + (c->n_dir-1)*c->n_elts;

#   pragma omp parallel if (c->n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < c->n_elts; i++)
      q_n_dir[i] = c->df[i]*coef;  /* Q(n_dir, :) = df/|df|_L2 */

  } /* n_dir > 1 */

  /* Step 3: Improve the condition number of R if needed
   * ------- */

  if (c->param.max_cond > 0) {

#if defined(DEBUG) && !defined(NDEBUG) && CS_ITER_ALGO_DBG > 0
    cs_log_printf(CS_LOG_DEFAULT, " %s:l%d >> Init  cond: %5.3e; n_dir=%d\n",
                  __func__, __LINE__,
                  _condition_number(c->n_dir, c->R), c->n_dir);
#endif

    while (c->n_dir > 1 &&
           _condition_number(c->n_dir, c->R) > c->param.max_cond) {

      _qrdelete(m_max, c->n_elts, c->Q, c->R);
      c->n_dir--;

    }

#if defined(DEBUG) && !defined(NDEBUG) && CS_ITER_ALGO_DBG > 0
    cs_log_printf(CS_LOG_DEFAULT, " %s:l%d >> Final cond: %5.3e; n_dir=%d\n",
                  __func__, __LINE__,
                  _condition_number(c->n_dir, c->R), c->n_dir);
#endif

  }

  /* Step 4: Solve the least square problem (upper triangular solve)
   * ------- */

  _anderson_compute_gamma(c, dotprod);

  /* Step 5: Update cur_iterate by cur_iterate = cur_iterate[i] - gamma.dg
   * ------- */

  for (int j = 0; j < c->n_dir; j++) {

    const cs_real_t  *dg_j = c->dg + j*c->n_elts;
    const double  gamma_j = c->gamma[j];

#   pragma omp parallel if (c->n_elts > CS_THR_MIN)
    for (cs_lnum_t l = 0; l < c->n_elts; l++)
      cur_iterate[l] -= gamma_j * dg_j[l];

  }

  /* Step 6: Damping of cur_iterate
   * ------- */

  if ((c->param.beta > 0.0) && (c->param.beta < 1.0))
    _anderson_damping(c, cur_iterate);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the counting of inner iterations in two-level algorithms
 *
 * \param[in, out] algo               pointer to the structure to update
 * \param[in]      n_last_inner_iter  last number of inner loop iterations
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_update_inner_iters(cs_iter_algo_t   *algo,
                                int               n_last_inner_iter)
{
  if (algo == nullptr)
    return;
  if (cs_flag_test(algo->type, CS_ITER_ALGO_TWO_LEVEL) == false)
    return;

  if (algo->type & CS_ITER_ALGO_DEFAULT) {

    cs_iter_algo_default_t *c
      = static_cast<cs_iter_algo_default_t *>(algo->context);

    c->last_inner_iter = n_last_inner_iter;
    c->n_inner_iter += n_last_inner_iter;

  }
  else if (algo->type & CS_ITER_ALGO_ANDERSON) {

    cs_iter_algo_aac_t *c = static_cast<cs_iter_algo_aac_t *>(algo->context);

    c->last_inner_iter = n_last_inner_iter;
    c->n_inner_iter += n_last_inner_iter;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the value of the residual and its associated members.
 *
 * \param[in, out] algo      pointer to a cs_iter_algo_t structure
 * \param[in]      res       new value of the residual
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_update_residual(cs_iter_algo_t      *algo,
                             double               res)
{
  if (algo == nullptr)
    return;
  assert(res > -DBL_MIN);

  if (algo->type & CS_ITER_ALGO_DEFAULT) {

    cs_iter_algo_default_t *c
      = static_cast<cs_iter_algo_default_t *>(algo->context);

    c->prev_res = c->res;
    c->res = res;

    if (c->n_algo_iter < 1) /* Store the first residual to detect a possible
                               divergence of the algorithm */
      c->res0 = res;

  }
  else if (algo->type & CS_ITER_ALGO_ANDERSON) {

    cs_iter_algo_aac_t *c = static_cast<cs_iter_algo_aac_t *>(algo->context);

    c->prev_res = c->res;
    c->res = res;

    if (c->n_algo_iter < 1) /* Store the first residual to detect a possible
                               divergence of the algorithm */
      c->res0 = res;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the convergence status and the number of iterations. The
 *        tolerance threshold is given so that it has to be computed before
 *        calling this function.
 *
 * \param[in, out] algo      pointer to a cs_iter_algo_t structure
 * \param[in]      tol       tolerance threshold to apply
 *
 * \return the convergence state
 */
/*----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_iter_algo_update_cvg_tol_given(cs_iter_algo_t         *algo,
                                  double                  tol)
{
  cs_iter_algo_set_tolerance(algo, tol);

  if (algo->type & CS_ITER_ALGO_DEFAULT) {

    cs_iter_algo_default_t *c
      = static_cast<cs_iter_algo_default_t *>(algo->context);

    /* Increment the number of iterations */

    c->n_algo_iter += 1;

    /* Set the convergence status */

    if (c->res < c->tol)
      c->cvg_status = CS_SLES_CONVERGED;

    else if (c->n_algo_iter >= algo->cvg_param.n_max_iter)
      c->cvg_status = CS_SLES_MAX_ITERATION;

    else if (c->res > algo->cvg_param.dtol * c->prev_res ||
             c->res > algo->cvg_param.dtol * c->res0)
      c->cvg_status = CS_SLES_DIVERGED;

    else
      c->cvg_status = CS_SLES_ITERATING;

    return c->cvg_status;

  }
  else if (algo->type & CS_ITER_ALGO_ANDERSON) {

    cs_iter_algo_aac_t *c = static_cast<cs_iter_algo_aac_t *>(algo->context);

    /* Increment the number of iterations */

    c->n_algo_iter += 1;

    /* Set the convergence status */

    if (c->res < c->tol)
      c->cvg_status = CS_SLES_CONVERGED;

    else if (c->n_algo_iter >= algo->cvg_param.n_max_iter)
      c->cvg_status = CS_SLES_MAX_ITERATION;

    else if (c->res > algo->cvg_param.dtol * c->prev_res ||
             c->res > algo->cvg_param.dtol * c->res0)
      c->cvg_status = CS_SLES_DIVERGED;

    else
      c->cvg_status = CS_SLES_ITERATING;

    return c->cvg_status;

  }
  else
    return CS_SLES_ITERATING;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the convergence status and the number of iterations. The
 *        tolerance threshold is automatically computed by a default formula
 *        relying on the relative tolerance scaled by the normalization factor
 *        and the absolute tolerance.
 *
 * \param[in, out] algo      pointer to a cs_iter_algo_t structure
 *
 * \return the convergence state
 */
/*----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_iter_algo_update_cvg_tol_auto(cs_iter_algo_t         *algo)
{
  /* Set the tolerance criterion (computed at each call if the normalization is
     modified between two successive calls) */

  double  tol = fmax(algo->cvg_param.rtol*cs_iter_algo_get_normalization(algo),
                     algo->cvg_param.atol);

  return cs_iter_algo_update_cvg_tol_given(algo, tol);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the convergence of the iterative algorithm
 *
 * \param[in] algo     pointer to the iterative algorithm structure
 * \param[in] label    label to specify the log
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_log_cvg(cs_iter_algo_t      *algo,
                     const char          *label)
{
  if (algo == nullptr)
    return;

  if (algo->verbosity < 1)
    return;

  if (cs_log_default_is_active() == false)
    return;

  if (algo->type & CS_ITER_ALGO_DEFAULT) {

    cs_iter_algo_default_t *c
      = static_cast<cs_iter_algo_default_t *>(algo->context);

    if (algo->type & CS_ITER_ALGO_TWO_LEVEL) {

      if (c->n_algo_iter == 1)
        cs_log_printf(CS_LOG_DEFAULT,
                      "##%s.It    | %10s  %10s  %10s  %10s\n", label,
                      "Algo.res", "Inner.iter", "Cumul.iter", "Tolerance");

      cs_log_printf(CS_LOG_DEFAULT,
                    "##%s.It%03d | %10.4e  %10d  %10d  %10.4e\n", label,
                    c->n_algo_iter, c->res, c->last_inner_iter, c->n_inner_iter,
                    c->tol);

    }
    else {

      if (c->n_algo_iter == 1)
        cs_log_printf(CS_LOG_DEFAULT,
                      "##%s.It    | %10s  %10s\n",
                      label, "Algo.res", "Tolerance");

      cs_log_printf(CS_LOG_DEFAULT,
                    "##%s.It%03d | %10.4e  %10.4e\n",
                    label, c->n_algo_iter, c->res, c->tol);

    }

  }
  else if (algo->type & CS_ITER_ALGO_ANDERSON) {

    cs_iter_algo_aac_t *c = static_cast<cs_iter_algo_aac_t *>(algo->context);

    if (algo->type & CS_ITER_ALGO_TWO_LEVEL) {

      if (c->n_algo_iter == 1)
        cs_log_printf(CS_LOG_DEFAULT,
                      "##%s.It    | %10s  %10s  %10s  %10s\n", label,
                      "Algo.res", "Inner.iter", "Cumul.iter", "Tolerance");

      cs_log_printf(CS_LOG_DEFAULT,
                    "##%s.It%03d | %10.4e  %10d  %10d  %10.4e\n", label,
                    c->n_algo_iter, c->res, c->last_inner_iter, c->n_inner_iter,
                    c->tol);

    }
    else {

      if (c->n_algo_iter == 1)
        cs_log_printf(CS_LOG_DEFAULT,
                      "##%s.It    | %10s  %10s\n",
                      label, "Algo.res", "Tolerance");

      cs_log_printf(CS_LOG_DEFAULT,
                    "##%s.It%03d | %10.4e  %10.4e\n",
                    label, c->n_algo_iter, c->res, c->tol);
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if something wrong happens during the iterative process after
 *        one new iteration
 *
 * \param[in] func_name    name of the calling function
 * \param[in] eq_name      name of the equation being solved
 * \param[in] algo_name    name of the iterative algo. used
 * \param[in] algo         pointer to the iterative algorithm structure
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_check_warning(const char          *func_name,
                           const char          *eq_name,
                           const char          *algo_name,
                           cs_iter_algo_t      *algo)
{
  if (algo == nullptr)
    return;

  cs_sles_convergence_state_t  cvg_status = cs_iter_algo_get_cvg_status(algo);

  if (cvg_status == CS_SLES_DIVERGED) {

    int  n_iter = cs_iter_algo_get_n_iter(algo);
    double  res = cs_iter_algo_get_residual(algo);

    bft_error(__FILE__, __LINE__, 0,
              "%s: %s algorithm divergence detected.\n"
              "%s: Equation \"%s\" can not be solved correctly.\n"
              "%s: Last iteration=%d; last residual=%5.3e\n",
              func_name, algo_name, func_name, eq_name, func_name, n_iter, res);

  }
  else if (cvg_status == CS_SLES_MAX_ITERATION) {

    int  n_iter = cs_iter_algo_get_n_iter(algo);
    double  res = cs_iter_algo_get_residual(algo);

    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_WARNINGS,
                  " %s: %s algorithm reaches the max. number of iterations"
                  " when solving equation \"%s\"\n"
                  " %s: max_iter=%d; last residual=%5.3e\n",
                  func_name, algo_name, eq_name, func_name, n_iter, res);

  }
}

/*----------------------------------------------------------------------------*/
END_C_DECLS
