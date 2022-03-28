/*============================================================================
 * Set of functions to manage high-level iterative algorithms
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
 * \brief Set of functions to handle the management of high-level iterative
 *        algorithms such as Uzawa, Golub-Kahan Bi-orthogonalization, block
 *        preconditioner or Picard algorithms which incorporates inner
 *        iterative solvers
 */

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

/*! \struct _cs_iter_algo_aa_t
 *  \brief Context structure for the algorithm called Anderson acceleration
 *
 *  Set of parameters and arrays to manage the Anderson acceleration
 */

struct _cs_iter_algo_aa_t {

/*!
 * \var param
 * Set of parameters driving the behavior of the Anderson acceleration
 */

  cs_iter_algo_param_aa_t    param;


  cs_lnum_t                  n_elts;

  /*!
   * @}
   * @name Work quantities (temporary)
   * @{
   *
   * \var n_dir
   * Number of directions currently at stake
   *
   * \var fold
   * Previous values for f
   *
   * \var df
   * Difference between the current and previous values of f
   *
   * \var gold
   * Previous values for g
   *
   * \var dg
   * Difference between the current and previous values of g
   *
   * \var Q
   * Matrix Q in the Q.R factorization (seen as a bundle of vectors)
   *
   * \var R
   * Matrix R in the Q.R factorization (small dense matrix)
   *
   * \var gamma
   * Coefficients used to define the linear combination
   *
   *@}
   */

  int             n_dir;

  cs_real_t      *fold;
  cs_real_t      *df;

  cs_real_t      *gold;
  cs_real_t      *dg;

  cs_real_t      *Q;
  cs_sdm_t       *R;

  cs_real_t      *gamma;

};

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

#define CS_ITER_ALGO_DBG      0

/*============================================================================
 * Private variables
 *============================================================================*/

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
 * \param[in, out] aa        pointer to the Anderson algo. structure
 */
/*----------------------------------------------------------------------------*/

static void
_aa_shift_dg(cs_iter_algo_aa_t   *aa)
{
  assert(aa->dg != NULL);
  const size_t  dir_size = sizeof(cs_real_t)*aa->n_elts;

  for (int i = 0; i < aa->n_dir - 1; i++) {

    cs_real_t *dg_i = aa->dg + i*aa->n_elts;
    cs_real_t *dg_ip1 = aa->dg + (i+1)*aa->n_elts;

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
 * \param[in, out] aa       pointer to the structure managing the Anderson algo.
 * \param[in]      dotprod  function to compute a dot product
 */
/*----------------------------------------------------------------------------*/

static void
_aa_compute_gamma(cs_iter_algo_aa_t           *aa,
                  cs_cdo_blas_dotprod_t       *dotprod)
{
  int  n_cols = aa->R->n_cols;

  for (int i = aa->n_dir-1; i >= 0; i--) {

    /* Solve R gamma = (fcur, Q_i) where fcur is equal to fold since one saves
       fcur into fold at the end of the step 1 of the Anderson algorithm */

    double  s = dotprod(aa->fold, aa->Q + i*aa->n_elts);

    for (int j = i+1; j < aa->n_dir; j++)
      s -= aa->R->val[i*n_cols+j] * aa->gamma[j];

    aa->gamma[i] = s/aa->R->val[i*(n_cols+1)];

  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_ITER_ALGO_DBG > 0
  cs_log_printf(CS_LOG_DEFAULT, " %s:l%d >> Gamma coefficients:",
                __func__, __LINE__);
  for (int i = 0; i < aa->n_dir; i++)
    cs_log_printf(CS_LOG_DEFAULT, " (%d:=%5.3e)", i, aa->gamma[i]);
  cs_log_printf(CS_LOG_DEFAULT, "\n");
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the coefficient used in the linear combination for the
 *        Anderson acceleration
 *
 * \param[in, out] aa    pointer to the structure managing the Anderson algo.
 * \param[in, out] x     array of values to update
 */
/*----------------------------------------------------------------------------*/

static void
_aa_damping(cs_iter_algo_aa_t    *aa,
            cs_real_t            *x)
{
  const double  omb = 1.0 - aa->param.beta;
  const int  m_max = aa->param.n_max_dir;
  const cs_real_t  *Rval = aa->R->val;

  for (int j = 0; j < aa->n_dir; j++) {

    double  R_gamma = 0.;
    for (int k = j; k < aa->n_dir; k++)  /* R is upper diag */
      R_gamma += Rval[j*m_max+k] * aa->gamma[k];

    /* x = x - (1-beta)*(fold - Q*R*gamma) */

    const cs_real_t  *Qj = aa->Q + j*aa->n_elts;  /* get row j */

#   pragma omp parallel if (aa->n_elts > CS_THR_MIN)
    for (cs_lnum_t l = 0; l < aa->n_elts; l++)
      x[l] += omb * (Qj[l] * R_gamma - aa->fold[l]);

  } /* Loop on directions */
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a new cs_iter_algo_t structure
 *
 * \param[in] param     main set of parameters driving the iterative algorithm
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_iter_algo_t *
cs_iter_algo_create(cs_iter_algo_param_t    param)
{
  cs_iter_algo_t  *ia = NULL;

  BFT_MALLOC(ia, 1, cs_iter_algo_t);

  ia->param.verbosity = param.verbosity;
  ia->param.atol = param.atol;
  ia->param.rtol = param.rtol;
  ia->param.dtol = param.dtol;
  ia->param.n_max_algo_iter = param.n_max_algo_iter;

  ia->normalization = 1.0;
  ia->context = NULL;

  cs_iter_algo_reset(ia);       /* default initialization */

  return ia;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if something wrong happens during the iterative process
 *         after one new iteration
 *
 * \param[in] func_name    name of the calling function
 * \param[in] eq_name      name of the equation being solved
 * \param[in] algo_name    name of the iterative algo. used
 * \param[in] ia           pointer to the iterative algo. structure
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_post_check(const char          *func_name,
                        const char          *eq_name,
                        const char          *algo_name,
                        cs_iter_algo_t      *ia)
{
  if (ia == NULL)
    return;

  if (ia->cvg == CS_SLES_DIVERGED)
    bft_error(__FILE__, __LINE__, 0,
              "%s: %s algorithm divergence detected.\n"
              "%s: Equation \"%s\" can not be solved correctly.\n"
              "%s: Last iteration=%d; last residual=%5.3e\n",
              func_name, algo_name,
              func_name, eq_name,
              func_name, ia->n_algo_iter, ia->res);

  else if (ia->cvg == CS_SLES_MAX_ITERATION) {

    cs_base_warn(__FILE__, __LINE__);
    bft_printf(" %s: %s algorithm reaches the max. number of iterations"
               " when solving equation \"%s\"\n"
               " %s: max_iter=%d; last residual=%5.3e\n",
               func_name, algo_name, eq_name,
               func_name, ia->param.n_max_algo_iter, ia->res);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the convergence state and the number of iterations
 *
 * \param[in, out] ia      pointer to a cs_iter_algo_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_update_cvg(cs_iter_algo_t         *ia)
{
  /* Set the tolerance criterion (computed at each call if the normalization is
     modified between two successive calls) */

  ia->tol = fmax(ia->param.rtol*ia->normalization, ia->param.atol);

  /* Increment the number of Picard iterations */

  ia->n_algo_iter += 1;

  /* Set the convergence status */

  if (ia->res < ia->tol)
    ia->cvg = CS_SLES_CONVERGED;

  else if (ia->n_algo_iter >= ia->param.n_max_algo_iter)
    ia->cvg = CS_SLES_MAX_ITERATION;

  else if (ia->res > ia->param.dtol * ia->prev_res ||
           ia->res > ia->param.dtol * ia->res0)
    ia->cvg = CS_SLES_DIVERGED;

  else
    ia->cvg = CS_SLES_ITERATING;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reset a cs_iter_algo_t structure in case of a non-linear algorothm
 *
 * \param[in]       type   type of non-linear algorithm
 * \param[in, out]  algo   pointer to a cs_iter_algo_t
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_reset_nl(cs_param_nl_algo_t   nl_algo_type,
                      cs_iter_algo_t      *algo)
{
  if (algo == NULL)
    return;

  cs_iter_algo_reset(algo);  /* common to all algorithm linear or non-linear */

  if (nl_algo_type == CS_PARAM_NL_ALGO_ANDERSON) {

    cs_iter_algo_aa_t  *aa = algo->context;
    assert(aa != NULL);
    aa->n_dir = 0;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a new cs_iter_algo_aa_t structure
 *
 * \param[in] aap             set of parameters for the Anderson acceleration
 * \param[in] n_elts          number of elements by direction
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_iter_algo_aa_t *
cs_iter_algo_aa_create(cs_iter_algo_param_aa_t    aap,
                       cs_lnum_t                  n_elts)
{
  cs_iter_algo_aa_t  *aa = NULL;

  BFT_MALLOC(aa, 1, cs_iter_algo_aa_t);

  /* Parameters */

  aa->param.n_max_dir = aap.n_max_dir;
  aa->param.starting_iter = aap.starting_iter;
  aa->param.max_cond = aap.max_cond;
  aa->param.beta = aap.beta;
  aa->param.dp_type = aap.dp_type;

  /* Other variables */

  aa->n_elts = n_elts;
  aa->n_dir = 0;

  /* Work arrays and structures */

  aa->fold = NULL;
  aa->df = NULL;
  aa->gold = NULL;
  aa->dg = NULL;
  aa->gamma = NULL;
  aa->Q = NULL;
  aa->R = NULL;

  return aa;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the set of parameters for an Anderson algorithm
 *
 * \param[in, out] ia      pointer to a cs_iter_algo_t structure
 *
 * \return a cs_iter_algo_param_aa_t structure
 */
/*----------------------------------------------------------------------------*/

cs_iter_algo_param_aa_t
cs_iter_algo_get_anderson_param(cs_iter_algo_t         *ia)
{
  if (ia != NULL) {

    cs_iter_algo_aa_t  *aa = ia->context;

    if (aa != NULL)
      return aa->param;

  }

  /* Define a set of parameters by default */

  cs_iter_algo_param_aa_t  aap = {
    .n_max_dir = 4,
    .starting_iter = 2,
    .max_cond = 500,
    .beta = 0.,
    .dp_type = CS_PARAM_DOTPROD_EUCLIDEAN };

  return aap;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate arrays useful for the Anderson acceleration
 *
 * \param[in, out] aa     pointer to the structure managing the Anderson algo.
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_aa_allocate_arrays(cs_iter_algo_aa_t  *aa)
{
  if (aa == NULL)
    return;

  int  n_max_dir = aa->param.n_max_dir;
  size_t  s = sizeof(cs_real_t)*aa->n_elts;

  BFT_MALLOC(aa->fold, aa->n_elts, cs_real_t);
  memset(aa->fold, 0, s);

  BFT_MALLOC(aa->df, aa->n_elts, cs_real_t);
  memset(aa->df, 0, s);

  BFT_MALLOC(aa->gold, aa->n_elts, cs_real_t);
  memset(aa->gold, 0, s);

  BFT_MALLOC(aa->dg, aa->n_elts*n_max_dir, cs_real_t);
  memset(aa->dg, 0, s*n_max_dir);

  BFT_MALLOC(aa->gamma, n_max_dir, cs_real_t);
  memset(aa->gamma, 0, sizeof(cs_real_t)*n_max_dir);

  BFT_MALLOC(aa->Q, aa->n_elts*n_max_dir, cs_real_t);
  memset(aa->Q, 0, s*n_max_dir);

  aa->R = cs_sdm_square_create(n_max_dir);
  cs_sdm_square_init(n_max_dir, aa->R);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free arrays used during the Anderson acceleration
 *
 * \param[in, out] aa     pointer to the structure managing the Anderson algo.
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_aa_free_arrays(cs_iter_algo_aa_t  *aa)
{
  if (aa == NULL)
    return;

  BFT_FREE(aa->fold);
  BFT_FREE(aa->df);

  BFT_FREE(aa->gold);
  BFT_FREE(aa->dg);

  BFT_FREE(aa->gamma);

  BFT_FREE(aa->Q);
  aa->R = cs_sdm_free(aa->R);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_iter_algo_aa_t structure used to manage the Anderson
 *         acceleration
 *
 * \param[in, out]  info
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_aa_free(cs_iter_algo_t  *info)
{
  if (info == NULL)
    return;

  cs_iter_algo_aa_t *aa = (cs_iter_algo_aa_t *)info->context;

  if (aa == NULL)
    return;

  cs_iter_algo_aa_free_arrays(aa);

  BFT_FREE(aa);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Apply one more iteration of the Anderson acceleration
 *
 * \param[in, out] ia            pointer to a cs_iter_algo_t structure
 * \param[in, out] cur_iterate   current iterate
 * \param[in]      pre_iterate   previous iterate
 * \param[in]      dotprod       function to compute a dot product
 * \param[in]      sqnorm        function to compute a square norm
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_aa_update(cs_iter_algo_t              *ia,
                       cs_real_t                   *cur_iterate,
                       const cs_real_t             *pre_iterate,
                       cs_cdo_blas_dotprod_t       *dotprod,
                       cs_cdo_blas_square_norm_t   *sqnorm)
{
  if (ia == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Structure not allocated.\n", __func__);

  cs_iter_algo_aa_t  *aa = (cs_iter_algo_aa_t *)ia->context;
  assert(aa != NULL);

  /* Check if anderson has to begin */

  const int shifted_iter = ia->n_algo_iter - aa->param.starting_iter;
  if (shifted_iter < 0)
    return;

  if (shifted_iter == 0 && aa->fold == NULL)
    cs_iter_algo_aa_allocate_arrays(aa);

  /*
   * Notation details:
   * ----------------
   * gcur = Gfunc(pre_iterate) = cur_iterate
   * fcur = gcur - pre_iterate = cur_iterate - pre_iterate
   *
   * dg = gcur - gold
   * df = fcur - fold
   */

  const int  m_max = aa->param.n_max_dir;
  const cs_real_t  *gcur = cur_iterate;

  /* Step 1: Update arrays useful for the Q.R factorization (df and dg)
   * -------
   */

  if (shifted_iter == 0) {

    /* The first iteration is simpler than the other one */

    assert(aa->n_dir == 0);

    /* Set fold and gold */

#   pragma omp parallel if (aa->n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < aa->n_elts; i++) {
      aa->fold[i] = gcur[i] - pre_iterate[i];
      aa->gold[i] = gcur[i];
    }

    return;  /* Nothing else to do. Algorithm has been initialized. This first
              * step is similar to a Picard iteration */

  }
  else {

    assert(shifted_iter > 0);

    /* Shift dg (set of rows) to the right position (two cases) */

    cs_real_t  *dg = aa->dg + aa->n_dir*aa->n_elts;

    if (aa->n_dir == m_max) {

      /* Erase the first row to retrieve some space and then
       * shift each superior row below */

      _aa_shift_dg(aa);
      dg = aa->dg + (aa->n_dir-1)*aa->n_elts;

    }

    /* Set dg, df, fold and gold */

#   pragma omp parallel if (aa->n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < aa->n_elts; i++) {

      cs_real_t  fcur = gcur[i] - pre_iterate[i];

      dg[i] = gcur[i] - aa->gold[i];
      aa->df[i] = fcur - aa->fold[i];
      aa->fold[i] = fcur;
      aa->gold[i] = gcur[i];

    }

    aa->n_dir++;

  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_ITER_ALGO_DBG > 1
  cs_log_printf(CS_LOG_DEFAULT, " %s:l%d >> Case n_dir=%d\n",
                __func__, __LINE__, aa->n_dir);
#endif

  /* Step 2: Update the Q.R factorization
   * -------
   */

  cs_real_t  *Rval = aa->R->val;

  if (aa->n_dir == 1) {

    const double  df_norm = sqrt(sqnorm(aa->df));
    const cs_real_t  coef = 1.0/df_norm;

    Rval[0] = df_norm; /* R(0,0) = |df|_L2 */

#   pragma omp parallel if (aa->n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < aa->n_elts; i++)
      aa->Q[i] = aa->df[i]*coef; /* Q(0) = df/|df|_L2 */

  }
  else {

    if (aa->n_dir > m_max) {   /* Remove the first column and last line in the
                                   Q.R factorization */

      _qrdelete(m_max, aa->n_elts, aa->Q, aa->R);
      aa->n_dir--;

    }

    for (int j = 0; j < aa->n_dir-1; j++) {

      cs_real_t *Qj = aa->Q + j*aa->n_elts; /* get the row j */

      const double  prod = dotprod(aa->df, Qj);

      Rval[j*m_max + (aa->n_dir-1)] = prod;  /* R(j, n_dir) = Qj*df */

#     pragma omp parallel if (aa->n_elts > CS_THR_MIN)
      for (cs_lnum_t l = 0; l < aa->n_elts; l++)
        aa->df[l] -= prod*Qj[l];  /* update df = df - R(j, n_dir)*Qj */

    }

    const double  df_norm = sqrt(sqnorm(aa->df));
    const cs_real_t  coef = 1.0/df_norm;

    /* R(n_dir,n_dir) = |df|_L2 */

    Rval[(aa->n_dir-1)*m_max+(aa->n_dir-1)] = df_norm;

    /* Set the row n_dir-1 of Q */

    cs_real_t *q_n_dir = aa->Q + (aa->n_dir-1)*aa->n_elts;

#   pragma omp parallel if (aa->n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < aa->n_elts; i++)
      q_n_dir[i] = aa->df[i]*coef;  /* Q(n_dir, :) = df/|df|_L2 */

  } /* n_dir > 1 */

  /* Step 3: Improve the condition number of R if needed
   * -------
   */

  if (aa->param.max_cond > 0) {

#if defined(DEBUG) && !defined(NDEBUG) && CS_ITER_ALGO_DBG > 0
    cs_log_printf(CS_LOG_DEFAULT, " %s:l%d >> Init  cond: %5.3e; n_dir=%d\n",
                  __func__, __LINE__,
                  _condition_number(aa->n_dir, aa->R), aa->n_dir);
#endif

    while (aa->n_dir > 1 &&
           _condition_number(aa->n_dir, aa->R) > aa->param.max_cond) {

      _qrdelete(m_max, aa->n_elts, aa->Q, aa->R);
      aa->n_dir--;

    }

#if defined(DEBUG) && !defined(NDEBUG) && CS_ITER_ALGO_DBG > 0
    cs_log_printf(CS_LOG_DEFAULT, " %s:l%d >> Final cond: %5.3e; n_dir=%d\n",
                  __func__, __LINE__,
                  _condition_number(aa->n_dir, aa->R), aa->n_dir);
#endif

  }

  /* Step 4: Solve the least square problem (upper triangular solve)
   * -------
   */

  _aa_compute_gamma(aa, dotprod);

  /* Step 5: Update cur_iterate by cur_iterate = cur_iterate[i] - gamma.dg
   * -------
   */

  for (int j = 0; j < aa->n_dir; j++) {

    const cs_real_t *dg_j = aa->dg + j*aa->n_elts;
    const double  gamma_j = aa->gamma[j];

#   pragma omp parallel if (aa->n_elts > CS_THR_MIN)
    for (cs_lnum_t l = 0; l < aa->n_elts; l++)
      cur_iterate[l] -= gamma_j * dg_j[l];

  }

  /* Step 4: Damping of cur_iterate
   * -------
   */

  if ((aa->param.beta > 0.0) && (aa->param.beta < 1.0))
    _aa_damping(aa, cur_iterate);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
