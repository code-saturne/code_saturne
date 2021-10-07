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

#include "cs_blas.h"
#include "cs_cdo_sqnorm.h"
#include "cs_math.h"
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
   * Booelan to store is the Anderson is activated (in case of delay for the
   * activation of the Anderson acceleration)
   *
   * \var n_aa_iter
   *
   *
   * \var mAA
   *
   *
   * \var rowNumber
   *
   *
   * \var AAFirst
   *
   *
   *
   *
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

  int             n_aa_iter;
  int             mAA;
  int             rowNumber;
  int             AAFirst;

  cs_real_t      *fval;
  cs_real_t      *fold;
  cs_real_t      *df;

  cs_real_t      *gval;
  cs_real_t      *gold;
  cs_real_t      *DG;

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
 * \brief Remove rows during the Anderson process
 *
 * \param[in, out] aa    pointer to the structure managing the Anderson algo.
 */
/*----------------------------------------------------------------------------*/

static void
_shiftRows_aa(cs_iter_algo_aa_t    *aa)
{
  size_t  row_size = sizeof(cs_real_t)*aa->n_elts;

  for (int j = 0; j < aa->rowNumber-1; j++) { /* delayed by 1 */

    cs_real_t *DGj = aa->DG + j*aa->n_elts;
    cs_real_t *DG_jp1 = aa->DG + (j+1)*aa->n_elts;

    memcpy(DGj, DG_jp1, row_size); /* DGj <= DG_(j+1) */

  }

  cs_real_t  *DG_max = aa->DG + (aa->rowNumber-1)*aa->n_elts;
  memcpy(DG_max, aa->DG, row_size);

  aa->rowNumber--;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute an estimation of the condition number of the matrix R
 *
 * \param[in]  R       pointer to the matrix R
 * \param[in]  n_rows  number of rows to consider
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

    vmax = fmax(vmax, fabs(R->val[i+i*n_cols]));
    vmin = fmin(vmin, fabs(R->val[i+i*n_cols]));

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
  const int  n_cols = R->n_cols;

  for (int i = 0; i < m-1; i++) {

    assert(R->n_rows == R->n_cols); /* square matrix */

    cs_real_t  *_Ri = R->val + i*n_cols;
    cs_real_t  *_Rip1 = R->val + (i+1)*n_cols;

    const double  temp = sqrt(_Ri[i+1]*_Ri[i+1] + _Rip1[i+1]*_Rip1[i+1]);
    const double  c = _Ri[i+1]/temp;
    const double  s = _Rip1[i+1]/temp;

     _Ri[i+1] = temp;
     _Rip1[i+1] = 0.0;

     if (i < m-2) { //diff % matlab

       for (int j = i+2; j < m; j++) {
         _Rip1[j] = -s*_Ri[j] + c*_Rip1[j];
         _Ri[j] = c*_Ri[j] + s*_Rip1[j];
       }

     }

     cs_real_t *Qi = Q + i*n_elts;
     cs_real_t *Qip1 = Q + (i+1)*n_elts;

     for (cs_lnum_t l = 0; l < n_elts; l++) {
       Qip1[l] = -s*Qi[l] + c*Qip1[l];
       Qi[l] = c*Qi[l] + s*Qip1[l];
     }

  } /* Loop on m */

  cs_sdm_t *R_temp = cs_sdm_square_create(R->n_cols) ;
  cs_sdm_square_init(R->n_cols, R_temp);

  /* Diff % matlab the Q and R matrices are not resized
   * A block version should be used */

  for (int i = 0; i < m-1; i++) {

    const cs_real_t  *_Ri = R->val + i*n_cols;
    cs_real_t  *_Rti = R_temp->val + i*n_cols;

    for (int j = i; j < m-1; j++)
      _Rti[j] = _Ri[j+1];

  }

  cs_sdm_copy(R, R_temp);
  R_temp = cs_sdm_free(R_temp);

  cs_real_t *Q_imax = Q + (m-1)*n_elts;
  memset(Q_imax, 0, sizeof(cs_real_t)*n_elts);
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
 * \param[in] n_max_dir       max. number of directions stored
 * \param[in] starting_iter   iteration at which the algorithm starts
 * \param[in] droptol         tolerance under which terms are dropped
 * \param[in] beta            relaxation coefficient
 * \param[in] n_elts          number of elements by direction
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_iter_algo_aa_t *
cs_iter_algo_aa_create(int                   n_max_dir,
                       int                   starting_iter,
                       cs_real_t             droptol,
                       cs_real_t             beta,
                       cs_lnum_t             n_elts)
{
  cs_iter_algo_aa_t  *aa = NULL;

  BFT_MALLOC(aa, 1, cs_iter_algo_aa_t);

  /* Parameters */

  aa->param.n_max_dir = n_max_dir;
  aa->param.starting_iter = starting_iter;
  aa->param.droptol = droptol;
  aa->param.beta = beta;

  aa->n_elts = n_elts;

  /* Monitoring variables */

  aa->is_initialized = false;
  aa->is_activated = false;

  aa->n_aa_iter = 0;
  aa->mAA = 0;
  aa->rowNumber = 0;
  aa->AAFirst = 0;

  /* Work arrays and structures */

  aa->fval = NULL;
  aa->fold = NULL;
  aa->df = NULL;
  aa->gval = NULL;
  aa->gold = NULL;
  aa->DG = NULL;
  aa->gamma = NULL;
  aa->Q = NULL;
  aa->R = NULL;

  return aa;
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

  BFT_MALLOC(aa->DG, aa->n_elts*n_max_dir, cs_real_t);
  memset(aa->DG, 0, s*n_max_dir);

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
  BFT_FREE(aa->DG);

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
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_aa_update(cs_iter_algo_info_t         *iai,
                       cs_real_t                   *cur_iterate)
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

    if (aa->gval == NULL)
      cs_iter_algo_aa_allocate_arrays(aa);

    memcpy(aa->gval, cur_iterate, dir_size); /* set gval = cur_iterate */
    memcpy(aa->fval, cur_iterate, dir_size); /* set fval = cur_iterate */

    aa->AAFirst = 1;
    aa->is_initialized = true;

  }
  else if (aa->is_initialized) {

    int  _mMax = aa->param.n_max_dir;

    memcpy(aa->gval, cur_iterate, dir_size); /* set gval = cur_iterate */

    /* set fval = gval - cur_iterate */

#   pragma omp parallel if (aa->n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < aa->n_elts; i++)
      aa->fval[i] = cur_iterate[i] - aa->fval[i];

    if (aa->AAFirst == 0) {

#     pragma omp parallel if (aa->n_elts > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < aa->n_elts; i++) {

        aa->df[i] = aa->fval[i] - aa->fold[i];   /* df = fval - fold */
        aa->gold[i] = aa->gval[i] - aa->gold[i]; /* gold = gval - gold */

      }

      /* DG append(gold) */

      if (aa->rowNumber == _mMax)
        _shiftRows_aa(aa);

      cs_real_t *DG_max = aa->DG + aa->rowNumber*aa->n_elts;
      memcpy(DG_max, aa->gold, dir_size);

      /* fin DG_.append(gold_) */

      aa->rowNumber++;
      aa->mAA++;  /* size of DG */

    } /* AAFirst == 0 */

    memcpy(aa->fold, aa->fval, dir_size); /* fold = fval */
    memcpy(aa->gold, aa->gval, dir_size); /* gold = gval */

    if (aa->mAA == 1) {

      double  square_norm_df = cs_dot_xx(aa->n_elts, aa->df);

      cs_parall_sum(1, CS_DOUBLE, &square_norm_df);

      aa->R->val[0] = sqrt(square_norm_df); /* R(0,0) = |df|_L2 */

      cs_real_t  *Q0 = aa->Q; /* get row 0 */
      const cs_real_t  coef = 1.0/aa->R->val[0];

#     pragma omp parallel if (aa->n_elts > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < aa->n_elts; i++)
        Q0[i] = aa->df[i]*coef; /* Q(0) = df/|df|_L2 */

    }
    else if (aa->mAA > 1) {

      if (aa->mAA == (_mMax+1)) {

        _qrdelete(_mMax, aa->n_elts, aa->Q, aa->R);
        aa->mAA = aa->mAA-1;

      }

      for (int j = 0; j < aa->mAA-1; j++) {

        cs_real_t *Qj = aa->Q + j*aa->n_elts; /* get row j */

        const double  prod = cs_gdot(aa->n_elts, aa->df, Qj);

        aa->R->val[j*_mMax + (aa->mAA-1)] = prod;  /* R(j, mAA) = Qj*df */

#       pragma omp parallel if (aa->n_elts > CS_THR_MIN)
        for (cs_lnum_t l = 0; l < aa->n_elts; l++)
          aa->df[l] -= prod*Qj[l];  /* update df = df - R(j, mAA)*Qj */

      }

      double  square_norm_df = cs_dot_xx(aa->n_elts, aa->df);

      cs_parall_sum(1, CS_DOUBLE, &square_norm_df);

      /* R(mAA,mAA) = |df|_L2 */

      aa->R->val[(aa->mAA-1)*_mMax+(aa->mAA-1)] = sqrt(square_norm_df);

      const cs_real_t  coef = 1.0/aa->R->val[(aa->mAA-1)*_mMax+(aa->mAA-1)];

      cs_real_t *Q_mAA = aa->Q + (aa->mAA-1)*aa->n_elts;  /* get row mAA-1 */

#     pragma omp parallel if (aa->n_elts > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < aa->n_elts; i++)
        Q_mAA[i] = aa->df[i]*coef;  /* Q(mAA-1) = df/|df|_L2 */

    }

    /* Drop residuals to improve conditioning if necessary */

    if ((aa->param.droptol > 0) && (aa->mAA > 0)) {

      /* Compute the first condition number of R */

      cs_real_t  cond = _condition_number(aa->rowNumber, aa->R);

      while ((cond > aa->param.droptol) && (aa->mAA>1)) {

        _qrdelete(aa->mAA, aa->n_elts, aa->Q, aa->R);
        _shiftRows_aa(aa); /* shift rows of DG */
        aa->mAA--;

        /* Evaluate the new condition number */

        cond = _condition_number(aa->rowNumber, aa->R);

      } /* End of while */

    } /* if drop tol */

    /* Solve least square problem */
    /* -------------------------- */

    /* Upper Triangular Solve */

    for (int i = aa->mAA-1; i >= 0; i--) {

      double s = cs_gdot(aa->n_elts, aa->fval, aa->Q + i*aa->n_elts);
      for (int j = i+1; j < aa->mAA; j++)
        s -= aa->R->val[i*_mMax+j] * aa->gamma[j];

      aa->gamma[i] = s/aa->R->val[i*_mMax+i];

    }

    /* Update cur_iterate by cur_iterate = gval[i] - gamma.DG */

    for (int j = 0; j < aa->mAA; j++) {

      const cs_real_t *DGj = aa->DG + j*aa->n_elts;
      for (cs_lnum_t l = 0; l < aa->n_elts; l++)
        cur_iterate[l] -= aa->gamma[j] * DGj[l];

    }

    if ((aa->param.beta > 0.0) && (aa->param.beta < 1.0)) {

      /* Relaxation */

      const double  omb = 1.0 - aa->param.beta;

      for (int j = 0; j < aa->mAA; j++) {

        double  R_gamma = 0.;
        for (int k = j; k < aa->mAA; k++)  /* R is upper diag */
          R_gamma += aa->R->val[j*_mMax+k] * aa->gamma[k];

        /* cur_iterate = cur_iterate - (1-beta)*(fval-Q*R*gamma) */

        const cs_real_t  *Qj = aa->Q + j*aa->n_elts;  /* get row j */
#       pragma omp parallel if (aa->n_elts > CS_THR_MIN)
        for (cs_lnum_t l = 0; l < aa->n_elts; l++)
          cur_iterate[l] += omb * (Qj[l] * R_gamma - aa->fval[l]);

      }

    } /* Relaxation */

    memcpy(aa->fval, cur_iterate, dir_size); /* fval = cur_iterate */

    aa->AAFirst = 0;
    aa->n_aa_iter = aa->n_aa_iter + 1;

  } /* aa->is_initialized */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
