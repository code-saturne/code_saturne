/*============================================================================
 * Set of functions to manage high-level iterative algorithms
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
   * \var is_initialized
   * Boolean to store if the initialization step has been done
   *
   * \var is_activated
   * Boolean to store is the Anderson is activated (in case of delay for the
   * activation of the Anderson acceleration)
   *
   * \var perform_first_iter
   * Boolean to store if one has to perform the first Anderson step which is
   * a bit different from the other iterations
   *
   * \var n_aa_iter
   * Current of Anderson acceleration performed
   *
   * \var n_dir
   * Number of directions currently at stake
   *
   * \var fval
   * Current values for f
   *
   * \var fold
   * Previous values for f
   *
   * \var df
   * Difference between the current and previous values of f
   *
   * \var gval
   * Current values for g
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

  bool            is_activated;
  bool            is_initialized;
  bool            perform_first_iter;

  int             n_aa_iter;
  int             n_dir;
  int             n_dg_cols;

  cs_real_t      *fval;
  cs_real_t      *fold;
  cs_real_t      *df;

  cs_real_t      *gval;
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

#define CS_ITER_ALGO_DBG      2

/*============================================================================
 * Private variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Remove one row during the Anderson acceleration process
 *
 * \param[in, out] aa    pointer to the structure managing the Anderson algo.
 */
/*----------------------------------------------------------------------------*/

static void
_aa_shift_rows(cs_iter_algo_aa_t    *aa)
{
  size_t  row_size = sizeof(cs_real_t)*aa->n_elts;

  aa->n_dg_cols--;

  for (int j = 0; j < aa->n_dg_cols; j++) {

    cs_real_t *dg_j = aa->dg + j*aa->n_elts;
    cs_real_t *dg_jp1 = aa->dg + (j+1)*aa->n_elts;

    memcpy(dg_j, dg_jp1, row_size); /* dg_j <= dg_(j+1) */

  }

  memcpy(aa->dg + aa->n_dg_cols*aa->n_elts, aa->dg, row_size);
}

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

    const double  temp = sqrt(R_i[i+1]*R_i[i+1] + R_ip1[i+1]*R_ip1[i+1]);
    const double  c = R_i[i+1]/temp;
    const double  s = R_ip1[i+1]/temp;

     R_i[i+1] = temp;
     R_ip1[i+1] = 0.0;

     if (i < m-2) { //diff % matlab

       for (int j = i+2; j < m; j++) {
         R_ip1[j] = -s*R_i[j] + c*R_ip1[j];
         R_i[j] = c*R_i[j] + s*R_ip1[j];
       }

     }

     cs_real_t *Q_i = Q + i*n_elts;
     cs_real_t *Q_ip1 = Q + (i+1)*n_elts;

     for (cs_lnum_t l = 0; l < n_elts; l++) {
       Q_ip1[l] = -s*Q_i[l] + c*Q_ip1[l];
       Q_i[l] = c*Q_i[l] + s*Q_ip1[l];
     }

  } /* Loop on m */

  cs_sdm_t *R_temp = cs_sdm_square_create(n_cols) ;

  cs_sdm_square_init(n_cols, R_temp);

  /* The Q and R matrices are not resized.
   * A block version should be used */

  for (int i = 0; i < m-1; i++) {

    const cs_real_t  *R_i = R->val + i*n_cols;
    cs_real_t  *Rt_i = R_temp->val + i*n_cols;

    for (int j = i; j < m-1; j++)
      Rt_i[j] = R_i[j+1];

  }

  cs_sdm_copy(R, R_temp);
  R_temp = cs_sdm_free(R_temp);

  cs_real_t *Q_imax = Q + (m-1)*n_elts;
  memset(Q_imax, 0, sizeof(cs_real_t)*n_elts);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the coefficient used in the linear combination for the
 *        Anderson acceleration
 *
 * \param[in, out] aa        pointer to the structure managing the Anderson algo.
 * \param[in]      dotprod   function to compute a dot product
 */
/*----------------------------------------------------------------------------*/

static void
_aa_compute_gamma(cs_iter_algo_aa_t           *aa,
                  cs_cdo_blas_dotprod_t       *dotprod)
{
  int  n_cols = aa->R->n_cols;

  for (int i = aa->n_dir-1; i >= 0; i--) {

    /* s = (fval, Q_i) */

    double  s = dotprod(aa->fval, aa->Q + i*aa->n_elts);

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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a new cs_iter_algo_info_t structure
 *
 * \param[in] verbosity    set the level of information printed
 * \param[in] n_max_iter   maximal number of iteration
 * \param[in] atol         absolute tolerance
 * \param[in] rtol         relative tolerance
 * \param[in] dtol         divergence tolerance
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_iter_algo_info_t *
cs_iter_algo_create(int          verbosity,
                    int          n_max_iter,
                    double       atol,
                    double       rtol,
                    double       dtol)
{
  cs_iter_algo_info_t  *iai = NULL;

  BFT_MALLOC(iai, 1, cs_iter_algo_info_t);

  iai->verbosity = verbosity;
  iai->atol = atol;
  iai->rtol = rtol;
  iai->dtol = dtol;
  iai->n_max_algo_iter = n_max_iter;
  iai->normalization = 1.0;

  iai->context = NULL;

  cs_iter_algo_reset(iai);

  return iai;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if something wrong happens during the iterative process
 *         after one new iteration
 *
 * \param[in] func_name    name of the calling function
 * \param[in] eq_name      name of the equation being solved
 * \param[in] algo_name    name of the iterative algo. used
 * \param[in] iai          pointer to the iterative algo. structure
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_post_check(const char            *func_name,
                        const char            *eq_name,
                        const char            *algo_name,
                        cs_iter_algo_info_t   *iai)
{
  if (iai == NULL)
    return;

  if (iai->cvg == CS_SLES_DIVERGED)
    bft_error(__FILE__, __LINE__, 0,
              "%s: %s algorithm divergence detected.\n"
              "%s: Equation \"%s\" can not be solved correctly.\n"
              "%s: Last iteration=%d; last residual=%5.3e\n",
              func_name, algo_name,
              func_name, eq_name,
              func_name, iai->n_algo_iter, iai->res);

  else if (iai->cvg == CS_SLES_MAX_ITERATION) {

    cs_base_warn(__FILE__, __LINE__);
    bft_printf(" %s: %s algorithm reaches the max. number of iterations"
               " when solving equation \"%s\"\n"
               " %s: max_iter=%d; last residual=%5.3e\n",
               func_name, algo_name, eq_name,
               func_name, iai->n_max_algo_iter, iai->res);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the convergence state and the number of iterations
 *
 * \param[in, out] iai      pointer to a cs_iter_algo_info_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_update_cvg(cs_iter_algo_info_t         *iai)
{
  /* Set the tolerance criterion (computed at each call if the normalization is
     modified between two successive calls) */

  iai->tol = fmax(iai->rtol*iai->normalization, iai->atol);

  /* Increment the number of Picard iterations */

  iai->n_algo_iter += 1;

  /* Set the convergence status */

  if (iai->res < iai->tol)
    iai->cvg = CS_SLES_CONVERGED;

  else if (iai->n_algo_iter >= iai->n_max_algo_iter)
    iai->cvg = CS_SLES_MAX_ITERATION;

  else if (iai->res > iai->dtol*iai->prev_res || iai->res > iai->dtol*iai->res0)
    iai->cvg = CS_SLES_DIVERGED;

  else
    iai->cvg = CS_SLES_ITERATING;
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

  aa->n_elts = n_elts;

  /* Monitoring variables */

  aa->is_initialized = false;
  aa->is_activated = false;
  aa->perform_first_iter = false;

  aa->n_aa_iter = 0;
  aa->n_dir = 0;
  aa->n_dg_cols = 0;

  /* Work arrays and structures */

  aa->fval = NULL;
  aa->fold = NULL;
  aa->df = NULL;
  aa->gval = NULL;
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
 * \param[in, out] iai      pointer to a cs_iter_algo_info_t structure
 *
 * \return a cs_iter_algo_param_aa_t structure
 */
/*----------------------------------------------------------------------------*/

cs_iter_algo_param_aa_t
cs_iter_algo_get_anderson_param(cs_iter_algo_info_t         *iai)
{
  cs_iter_algo_param_aa_t  aap = {
    .n_max_dir = 4,
    .starting_iter = 2,
    .max_cond = 500,
    .beta = 0.,
    .dp_type = CS_PARAM_DOTPROD_EUCLIDEAN };

  if (iai != NULL) {

    cs_iter_algo_aa_t  *aa = iai->context;

    if (aa != NULL)
      return aa->param;

  }

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

  BFT_MALLOC(aa->fval, aa->n_elts, cs_real_t);
  memset(aa->fval, 0, s);

  BFT_MALLOC(aa->fold, aa->n_elts, cs_real_t);
  memset(aa->fold, 0, s);

  BFT_MALLOC(aa->df, aa->n_elts, cs_real_t);
  memset(aa->df, 0, s);

  BFT_MALLOC(aa->gval, aa->n_elts, cs_real_t);
  memset(aa->gval, 0, s);

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

  BFT_FREE(aa->fval);
  BFT_FREE(aa->fold);
  BFT_FREE(aa->df);

  BFT_FREE(aa->gval);
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
cs_iter_algo_aa_free(cs_iter_algo_info_t  *info)
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
 * \param[in, out] iai           pointer to a cs_iter_algo_info_t structure
 * \param[in, out] cur_iterate   current iterate
 * \param[in]      dotprod       function to compute a dot product
 * \param[in]      sqnorm        function to compute a square norm
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_aa_update(cs_iter_algo_info_t         *iai,
                       cs_real_t                   *cur_iterate,
                       cs_cdo_blas_dotprod_t       *dotprod,
                       cs_cdo_blas_square_norm_t   *sqnorm)
{
  if (iai == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Structure not allocated.\n", __func__);

  cs_iter_algo_aa_t  *aa = (cs_iter_algo_aa_t *)iai->context;
  assert(aa != NULL);

  const size_t  dir_size = sizeof(cs_real_t)*aa->n_elts;

  /* Check if anderson has to begin */

  if (iai->n_algo_iter >= aa->param.starting_iter-1)
    aa->is_activated = true;
  else
    aa->is_activated = false;

  if (!aa->is_activated)
    return;

  if (!aa->is_initialized) {

    if (aa->gval == NULL) /* Allocate arrays */
      cs_iter_algo_aa_allocate_arrays(aa);

    memcpy(aa->gval, cur_iterate, dir_size); /* set gval = cur_iterate */
    memcpy(aa->fval, cur_iterate, dir_size); /* set fval = cur_iterate */

    aa->perform_first_iter = true;
    aa->is_initialized = true;

#if defined(DEBUG) && !defined(NDEBUG) && CS_ITER_ALGO_DBG > 1
    cs_log_printf(CS_LOG_DEFAULT, " %s:l%d >> Initialize Anderson acc.\n",
                  __func__, __LINE__);
#endif

    return;
  }

  /* Anderson acceleration has been initialized */

  int  m_max = aa->param.n_max_dir;

  memcpy(aa->gval, cur_iterate, dir_size); /* set gval = cur_iterate */

  /* set fval = gval - cur_iterate */

# pragma omp parallel if (aa->n_elts > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < aa->n_elts; i++)
    aa->fval[i] = cur_iterate[i] - aa->fval[i];

  if (!aa->perform_first_iter) {

    /* Not the first iteration of the Anderson algo */

#   pragma omp parallel if (aa->n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < aa->n_elts; i++) {

      aa->df[i] = aa->fval[i] - aa->fold[i];   /* df = fval - fold */
      aa->gold[i] = aa->gval[i] - aa->gold[i]; /* gold = gval - gold */

    }

    /* dg append gold */

    if (aa->n_dg_cols == m_max)
      _aa_shift_rows(aa);

    cs_real_t *dg_max = aa->dg + aa->n_dg_cols*aa->n_elts;
    memcpy(dg_max, aa->gold, dir_size);

    aa->n_dir++;
    aa->n_dg_cols++;

#if defined(DEBUG) && !defined(NDEBUG) && CS_ITER_ALGO_DBG > 1
    cs_log_printf(CS_LOG_DEFAULT, " %s:l%d >> Not the first iter. (n_dir: %d)\n",
                  __func__, __LINE__, aa->n_dir);
#endif

  } /* perform_first_iter == false */

  memcpy(aa->fold, aa->fval, dir_size); /* fold = fval */
  memcpy(aa->gold, aa->gval, dir_size); /* gold = gval */

  if (aa->n_dir == 1) {

    double  square_norm_df = sqnorm(aa->df);

    aa->R->val[0] = sqrt(square_norm_df); /* R(0,0) = |df|_L2 */

    cs_real_t  *Q0 = aa->Q; /* get row 0 */
    const cs_real_t  coef = 1.0/aa->R->val[0];

#   pragma omp parallel if (aa->n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < aa->n_elts; i++)
      Q0[i] = aa->df[i]*coef; /* Q(0) = df/|df|_L2 */

#if defined(DEBUG) && !defined(NDEBUG) && CS_ITER_ALGO_DBG > 1
    cs_log_printf(CS_LOG_DEFAULT, " %s:l%d >> Case n_dir=1\n",
                  __func__, __LINE__);
#endif

  }
  else if (aa->n_dir > 1) {

    if (aa->n_dir == (m_max+1)) {

      _qrdelete(m_max, aa->n_elts, aa->Q, aa->R);
      aa->n_dir--;

    }

    for (int j = 0; j < aa->n_dir-1; j++) {

      cs_real_t *Qj = aa->Q + j*aa->n_elts; /* get row j */

      const double  prod = dotprod(aa->df, Qj);

      aa->R->val[j*m_max + (aa->n_dir-1)] = prod;  /* R(j, n_dir) = Qj*df */

#     pragma omp parallel if (aa->n_elts > CS_THR_MIN)
      for (cs_lnum_t l = 0; l < aa->n_elts; l++)
        aa->df[l] -= prod*Qj[l];  /* update df = df - R(j, n_dir)*Qj */

    }

    double  square_norm_df = sqnorm(aa->df);

    /* R(n_dir,n_dir) = |df|_L2 */

    aa->R->val[(aa->n_dir-1)*m_max+(aa->n_dir-1)] = sqrt(square_norm_df);

    const cs_real_t  coef = 1.0/aa->R->val[(aa->n_dir-1)*m_max+(aa->n_dir-1)];

    /* Get row n_dir-1 */

    cs_real_t *q_n_dir = aa->Q + (aa->n_dir-1)*aa->n_elts;

#   pragma omp parallel if (aa->n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < aa->n_elts; i++)
      q_n_dir[i] = aa->df[i]*coef;  /* Q(n_dir-1) = df/|df|_L2 */

#if defined(DEBUG) && !defined(NDEBUG) && CS_ITER_ALGO_DBG > 1
    cs_log_printf(CS_LOG_DEFAULT, " %s:l%d >> Case n_dir=%d\n",
                  __func__, __LINE__, aa->n_dir);
#endif

  } /* n_dir > 1 */

  /* Improve the conditioning if necessary */

  if ((aa->param.max_cond > 0) && (aa->n_dir > 1)) {

    /* Compute the first condition number of R */

    cs_real_t  cond = _condition_number(aa->n_dir, aa->R);

#if defined(DEBUG) && !defined(NDEBUG) && CS_ITER_ALGO_DBG > 0
    cs_real_t  init_cond = cond;
    int  n_init_dir = aa->n_dir;
#endif

    while ((cond > aa->param.max_cond) && (aa->n_dir > 1)) {

      _qrdelete(aa->n_dir, aa->n_elts, aa->Q, aa->R);
      _aa_shift_rows(aa); /* shift rows of dg */
      aa->n_dir--;

      /* Evaluate the new condition number */

      cond = _condition_number(aa->n_dir, aa->R);

    } /* End of while */

#if defined(DEBUG) && !defined(NDEBUG) && CS_ITER_ALGO_DBG > 0
    if (n_init_dir > aa->n_dir)
      cs_log_printf(CS_LOG_DEFAULT, " %s:l%d >>\n"
                    "   %6d -->   %6d directions to improve conditioning.\n"
                    " %5.2e --> %5.2e for the condition number\n",
                    __func__, __LINE__, n_init_dir, aa->n_dir, init_cond, cond);
#endif

  } /* if drop tol */

  /* Solve the least square problem (upper triangular solve) */

  _aa_compute_gamma(aa, dotprod);

  /* Update cur_iterate by cur_iterate = gval[i] - gamma.dg */

  for (int j = 0; j < aa->n_dir; j++) {

    const cs_real_t *dg_j = aa->dg + j*aa->n_elts;

#   pragma omp parallel if (aa->n_elts > CS_THR_MIN)
    for (cs_lnum_t l = 0; l < aa->n_elts; l++)
      cur_iterate[l] -= aa->gamma[j] * dg_j[l];

  }

  if ((aa->param.beta > 0.0) && (aa->param.beta < 1.0)) {

    /* Relaxation */

    const double  omb = 1.0 - aa->param.beta;

    for (int j = 0; j < aa->n_dir; j++) {

      double  R_gamma = 0.;
      for (int k = j; k < aa->n_dir; k++)  /* R is upper diag */
        R_gamma += aa->R->val[j*m_max+k] * aa->gamma[k];

      /* cur_iterate = cur_iterate - (1-beta)*(fval-Q*R*gamma) */

      const cs_real_t  *Qj = aa->Q + j*aa->n_elts;  /* get row j */
#     pragma omp parallel if (aa->n_elts > CS_THR_MIN)
      for (cs_lnum_t l = 0; l < aa->n_elts; l++)
        cur_iterate[l] += omb * (Qj[l] * R_gamma - aa->fval[l]);

    }

  } /* Relaxation */

  memcpy(aa->fval, cur_iterate, dir_size); /* fval = cur_iterate */

  aa->perform_first_iter = false;
  aa->n_aa_iter += 1;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
