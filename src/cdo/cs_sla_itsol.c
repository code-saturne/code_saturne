/*============================================================================
 * Iterative solvers algorithms for Sparse Linear Algebra (SLA)
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_cdo_toolbox.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_sla.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_user_simu.c

  \brief Set main parameters for the current simulation

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define  SLA_ITSOL_DEBUG 0
#define  SLA_STOPPING_LENGTH 10

typedef struct { /* Monitoring of the convergence */

  int       size;      /* System size to monitor */

  double    init_res;
  double    best_res;
  int       best_ite;

  double    res_norm;
  double    sol_norm;
  double    rhs_norm;
  double    mat_norm;

  /* Advanced parameters */
  double    stag_threshold;
  int       stag_iter;
  double    div_factor;

} _cvg_monitor_t;

/*============================================================================
 * Private constant variables
 *============================================================================*/

/* Temporary buffers for linear algebra */
static size_t  _sla_itsol_main_size = 0;
static size_t  _sla_itsol_aux_size = 0;
static size_t  _sla_itsol_pcd_size = 0;
static cs_tmpbuf_t  *_sla_itsol_main_buffer = NULL;
static cs_tmpbuf_t  *_sla_itsol_aux_buffer = NULL;
static cs_tmpbuf_t  *_sla_itsol_pcd_buffer = NULL;

/*! \endcond (end ignore by Doxygen) */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*---------------------------------------------------------------------------*
 * Return the sign of a double
 *---------------------------------------------------------------------------*/

inline static int
_sign(double a)
{
  if (a < 0)
    return -1;
  else
    return  1;
}

/*----------------------------------------------------------------------------
 * dot prod with a stride
 *----------------------------------------------------------------------------*/

inline static double
_ddot(int  n, const double  *x, int  incx, const double  *y, int  incy)
{
  int  i, j, k;
  int  inc_x = CS_ABS(incx), inc_y = CS_ABS(incy);
  double  sum = 0;

  if (n < 0) return sum;

  if (inc_x == 1 && inc_y == 1)
    for (i = 0; i < n; i++)  sum += (x[i] * y[i]);
  else
    for (i = 0, j = 0, k = 0; i < n; i++, j += inc_x, k += inc_y)
      sum += (x[j] * y[k]);

  return sum;
}

/*----------------------------------------------------------------------------*
 * y <-- a*x + y (more optimized thanks to const and restrict)
 *----------------------------------------------------------------------------*/

inline static void
_daxpy(int               n,
       double            a,
       const double     *x,
       int               incx,
       double  *restrict y,
       int               incy)
{
  int  i, j, k;
  int  inc_x = CS_ABS(incx), inc_y = CS_ABS(incy);

  if (n < 0) return;

  if (inc_x == 1 && inc_y == 1)
    for (i = 0; i < n; i++) y[i] += (a * x[i]);
  else
    for (i = 0, j = 0, k = 0; i < n; i++, j += inc_x, k += inc_y)
      y[k] += (a * x[j]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute buffer size for GMRES iterative solver
 *
 * \param[in]     nx           system size
 * \param[inout]  aux_size     pointer to auxiliary size
 * \param[inout]  krylov_size  pointer to another auxiliary size
 */
/*----------------------------------------------------------------------------*/

static void
_get_krylov_sizes(int       nx,
                  size_t   *aux_size,
                  size_t   *krylov_size)
{
  int  ksize = 3;

  const  int  krylov_maxsize = 200;  /* Parameter */

  ksize = CS_MAX(ksize, (krylov_maxsize < (int)sqrt(nx)*1.73 ?
                         krylov_maxsize : (int)sqrt(nx)*1.73));

  *aux_size = 3*ksize + (ksize-1)*(nx + ksize);
  *krylov_size = ksize;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  SSOR preconditionning  (CSR with no symmetry storage needed)
 *
 * \param[in]    size    dimension of the linear problem to solve
 * \param[in]    dinv    diagonal inverse
 * \param[in]    A       pointer to a matrix structure
 * \param[in]    rk      residual at step k
 * \param[inout] pk      preconditionned array at step k

 * \param[in]    solver  pointer to a solver information structure
 */
/*----------------------------------------------------------------------------*/

static void
_ssor_precond(int                     size,
              const double           *dinv,
              const cs_sla_matrix_t  *A,
              const double           *rk,
              double                 *pk)
{
  int  i, j, colid;
  double  sum;

  const cs_lnum_t  *col = A->col, *idx = A->idx, *didx = A->didx;
  const double  *val = A->val;

  /* Sanity check */
  assert(A->type == CS_SLA_MAT_CSR);
  assert(A->flag & CS_SLA_MATRIX_LU);

  /* Forward solve */
  for (i = 0; i < size; i++) {
    sum = 0;
    for (j = idx[i]; j < didx[i]; j++) {  /* Lower part */
      colid = col[j]-1;
      sum += val[j]*pk[colid]*dinv[colid];
    }
    pk[i] = rk[i] - sum;
  }

  /* Backward solve */
  for (i = size-1; i > -1 ; i--) {
    sum = 0;
    for (j = didx[i]+1; j < idx[i+1]; j++) /* Upper part */
      sum += val[j]*pk[col[j]-1];
    pk[i] = dinv[i]*(pk[i] - sum);
  }

}

/*---------------------------------------------------------------------------*
 * Chebyshev  preconditionning   beta <-- spectral radius
 *---------------------------------------------------------------------------*/

static void
_spectral_preconditionning(int                     size,
                           double                  beta,
                           const cs_sla_matrix_t  *a,
                           const double           *rk,
                           double                 *gk)
{
  int ii;

  const double  alpha = beta/4.0;
  const double  c0 = 1.0/sqrt(beta*alpha);
  const double  sqab = sqrt(alpha/beta);
  const double  c1 = -c0*((1-sqab)/(1+sqab));
  const double  gamma = 2*c1/(beta-alpha);
  const double  omega = c1*((beta+alpha)/(beta-alpha)) - c0/2.0;

  cs_sla_matvec(a, rk, &gk, true);

  for (ii = 0; ii < size; ii++)
    gk[ii] = gamma*gk[ii] - omega*rk[ii];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Preconditionning switch routine
 *
 * \param[in]    size     dimension of the linear problem to solve
 * \param[in]    A        pointer to a matrix structure
 * \param[in]    rk       residual at step k
 * \param[in]    precond  pointer to a solver information structure
 * \param[inout] pk       preconditionned array at step k
 */
/*----------------------------------------------------------------------------*/

static void
_preconditionning(int                        size,
                  const cs_sla_matrix_t     *A,
                  const double               rk[],
                  cs_param_precond_type_t    precond,
                  double                     pk[])
{
  int  i;
  double  beta;

  switch(precond) {

  case CS_PARAM_PRECOND_NONE:
    memcpy(pk, rk, size*sizeof(double));
    break;
  case CS_PARAM_PRECOND_DIAG:
    {
      double  *invdiag = (double *)_sla_itsol_pcd_buffer->buf;
      for (i = 0; i < size; i++) pk[i] = rk[i]*invdiag[i];
    }
    break;
  case CS_PARAM_PRECOND_SSOR:
    {
      double  *invdiag = (double *)_sla_itsol_pcd_buffer->buf;
      _ssor_precond(size, invdiag, A, rk, pk);
    }
    break;
  case CS_PARAM_PRECOND_CHEBY:
    cs_sla_get_mat_spectral_radius(A,
                                   25,     /* Number max. of iterations */
                                   1.e-3,  /* Conv. epsilon */
                                   &beta);
    _spectral_preconditionning(size, beta, A, rk, pk);
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              " Invalid preconditionner. Stop the execution.\n");

  } /* End of switch */

}

/*---------------------------------------------------------------------------*
 * Test if we have to iterate one more time
 *---------------------------------------------------------------------------*/

static _Bool
_is_cvg(const cs_param_itsol_t    param,
        const double              residual[],
        const double              solution[],
        cs_sla_sumup_t           *ret,
        _cvg_monitor_t           *cvgm)
{
  _Bool  one_more = true, print_msg = false, print_warning = false;
  double  conv_val = param.eps;

  cvgm->res_norm = cs_euclidean_norm(cvgm->size, residual);
  cvgm->sol_norm = cs_euclidean_norm(cvgm->size, solution);
  ret->residual = cvgm->res_norm;

  if (cvgm->init_res > 0.1*DBL_MAX) { /* not initialized */
    cvgm->init_res = cvgm->res_norm;
    cvgm->best_res = cvgm->res_norm;
    cvgm->best_ite = ret->iter;
  }

  if (cvgm->best_res > cvgm->res_norm) { /* Update best residual */
    cvgm->best_res = cvgm->res_norm;
    cvgm->best_ite = ret->iter;
  }

  if (param.resid_normalized)
    conv_val *= cvgm->rhs_norm;

  if (cvgm->res_norm < cs_get_eps_machine() || cvgm->res_norm < conv_val) {
    one_more = false, print_msg = true;
    ret->code = CS_SLA_CODE_CVG;
  }

  if (ret->iter >= param.n_max_iter) {
    ret->code = CS_SLA_CODE_STOP;
    one_more = false, print_warning = true, print_msg = true;
  }

  /* Is there a divergence ? */
  if (cvgm->res_norm > cvgm->div_factor * cvgm->init_res) {
    ret->code = CS_SLA_CODE_DVG;
    one_more = false, print_warning = true, print_msg = true;
  }

  /* Is there a stagnation ? */
  if (cvgm->res_norm > cvgm->stag_threshold*cvgm->best_res)  {
    if (ret->iter - cvgm->best_ite > cvgm->stag_iter) {
      ret->code = CS_SLA_CODE_STAG;
      one_more = false;
      printf("\n WARNING  >>  STAGNATION DU SOLVEUR\n"
             " best res. %8.5e (ite: %d) current res. %8.5e (ite: %d)\n",
             cvgm->best_res, cvgm->best_ite, cvgm->res_norm, ret->iter);
    }
  }

  if (ret->iter % param.output_freq == 0)
    print_msg = true;

  if (print_msg)
    printf(" CONV | %5d | %8.5e | %8.5e | %8.5e |\n",
           ret->iter, ret->residual, conv_val, cvgm->sol_norm);

  if (print_warning)
    printf("\n WARNING  >>  NON CONVERGENCE DU SOLVEUR >> CODE: %d\n",
           ret->code);

  return one_more;
}

/*---------------------------------------------------------------------------*
 * Compute system information for building a stopping criterion
 *---------------------------------------------------------------------------*/

static _cvg_monitor_t
_init_solver_info(const cs_sla_matrix_t    *matrix,
                  const double              rhs[])
{
  const int size = matrix->n_cols;

  _cvg_monitor_t  cvgm;

  /* Initialize the convegence monitor */
  cvgm.size = size;
  cvgm.res_norm = DBL_MAX;
  cvgm.init_res = DBL_MAX;
  cvgm.best_res = DBL_MAX;
  cvgm.best_ite = -1;
  cvgm.sol_norm = DBL_MAX;
  cvgm.rhs_norm = cs_euclidean_norm(size, rhs);
  cvgm.mat_norm = cs_sla_get_matrix_norm(matrix);

  /* Advanced Parameters */
  cvgm.stag_threshold = 1.20; /* between 1.0 and 2.0 */
  cvgm.stag_iter = 200;
  cvgm.div_factor = 1e3;

  printf("# |A|: %8.5e ; |b|: %8.5e\n"
         "#     |  iter |  residual   |  criterion  |   norm(X)   |\n",
         cvgm.mat_norm, cvgm.rhs_norm);

  return cvgm;
}

/*----------------------------------------------------------------------------
 * Transform using Givens rotations a system Ax=b where A is an upper
 * triangular matrix of size n*(n-1) with a lower diagonal into an equivalent
 * system A'x=b' where A' is an upper triangular matrix
 *
 * parameters:
 *  a             <-> dense matrix
 *                    (vector of size a_size*a_size elements
 *                      a(i,j)=a[i+j*size_a])
 *                     input : A / output: A'
 *  a_size        <-- matrix dimension
 *  b             <-> pre-allocated vector of a_size elements
 *                     input : b / output : b'
 *  givens_coeff  <-> pre-allocated vector of a_size elements
 *                     input : previous Givens coefficients
 *                     output : updated Givens coefficients
 *  update_rank    <-- rank of first non null coefficient on lower diagonal
 *  end_update     <-- rank of last non null coefficient on lower diagonal
 *  aux_vectors    ---  optional work array
 *----------------------------------------------------------------------------*/

static int
_givens_rot_update(double  *a,
                   int      a_size,
                   double  *b,
                   double  *givens_coeff,
                   int      update_rank,
                   int      end_update)
{
  int  i, j;
  double  norm;
  double _aux;

  for (i = 0; i < update_rank; ++i) {
    for(j = update_rank; j < end_update; ++j) {

      _aux =  givens_coeff[i]*a[j*a_size + i]
             + givens_coeff[i + a_size] * a[j*a_size + i+1];

      a[j*a_size + i+1] =  givens_coeff[i] * a[i+1 + j*a_size]
                         - givens_coeff[i + a_size] * a[j*a_size + i];

      a[j*a_size + i] = _aux;

    }
  }

  for (i = update_rank; i < end_update ; ++i) {

    norm = pow(a[i*a_size+i],2) + pow(a[i*a_size + i + 1],2);
    norm = sqrt(norm);

    givens_coeff[a_size+i] = a[i*a_size + i + 1 ]/norm;
    givens_coeff[i] = a[i*a_size + i]/norm;

    b[i+1] = -b[i]*givens_coeff[a_size+i];
    b[i] = b[i]*givens_coeff[i];

    for(j = i; j < end_update; j++) {
      _aux =  givens_coeff[i]*a[j*a_size+i]
            + givens_coeff[a_size+i]*a[j*a_size+i+1];
      if (j == i)
        a[j*a_size+i+1]=0;
      else
        a[i+1 + j*a_size] =  givens_coeff[i]*a[i+1+j*a_size]
                           - givens_coeff[a_size+i]*a[i+j*a_size];

      a[j*a_size + i] = _aux;
    }

  }

  return 0;
}

/*----------------------------------------------------------------------------
 * Compute solution of Ax=b where A is an upper triangular matrix
 *
 * parameters:
 *  a         <-- dense upper triangular matrix A
 *                (vector of size glob_size*glob_size elements
 *                  a(i,j) = a[i+j*glob_size])
 *  a_size    <-- system size
 *  glob_size <-- a_size + halo size
 *  b         <-- pre-alocated vector of a_size elements
 *  x         --> system solution, pre-alocated vector of a_size elements
 *
 *          |         |          |   |x1|   |b1|
 *          |         |          |   |x2|   |b2|
 *          |    A    |   halo   |   |x3|   |b3|
 *          |         |          | * |x4| = |b4|
 *          |_________|          |   |x5|   |b5|
 *          |                    |   |h |   |h |
 *          |     halo           |   |a |   |a |
 *          |                    |   |l |   |l |
 *          |                    |   |o |   |o |
 *----------------------------------------------------------------------------*/

static int
_solve_diag_sup_halo(double  *a,
                     int      a_size,
                     int      glob_size,
                     double  *b,
                     double  *x)
{
  int  i, j;

  for (i = a_size-1; i > -1; i--) {
    x[i] = b[i];
    for (j = i+1; j < a_size; ++j)
      x[i] = x[i] - a[j*glob_size + i] * x[j];
    x[i] /= a[i*glob_size + i];
  }

  return 0;
}

/*----------------------------------------------------------------------------
 * Solution of (ad+ax).vx = Rhs using diagonal preconditioned GMRES.
 * On entry, vx is considered initialized.
 *----------------------------------------------------------------------------*/

static cs_sla_sumup_t
_gmres(const cs_sla_matrix_t   *a,
       const double            *rhs,
       double                  *vx,
       cs_param_itsol_t         param)
{
  int  cvg;
  int  n_rows, ii, jj, check_freq, l_iter, l_old_iter, scaltest;
  size_t  ui, uk, krylov_size, aux_size;
  double  beta, scal_prod, residue, epsi;

  unsigned  n_iter = 1;
  double  *_givens_coeff = NULL, *_beta = NULL;
  double  *_krylov_vectors = NULL, *_H_matrix = NULL, *ad_inv = NULL;
  double  *dk = NULL, *gk = NULL, *bk = NULL, *fk = NULL, *krk = NULL;

  /* JF : modif for ilu */
  cs_sla_matrix_t *M = NULL;

  cs_sla_sumup_t  ret = {CS_SLA_CODE_CONTINUE, 0, DBL_MAX};
  _cvg_monitor_t  cvgm = _init_solver_info(a, rhs);

  /* Preliminary calculations */
  n_rows = a->n_rows;
  _get_krylov_sizes(n_rows, &aux_size, &krylov_size);

  check_freq = (int)(krylov_size/5) + 1;
  epsi = 1.e-15;
  scaltest = 0;
  residue = -1;

  ad_inv = (double *)_sla_itsol_pcd_buffer->buf;
  dk     = (double *)_sla_itsol_main_buffer->buf;
  gk     = (double *)_sla_itsol_main_buffer->buf + n_rows;
  bk     = (double *)_sla_itsol_main_buffer->buf + 2*n_rows;
  fk     = (double *)_sla_itsol_main_buffer->buf + 3*n_rows;

  _krylov_vectors = (double *)_sla_itsol_aux_buffer->buf;
  _H_matrix       = (double *)_sla_itsol_aux_buffer->buf
    + (krylov_size-1)*n_rows;
  _givens_coeff   = (double *)_sla_itsol_aux_buffer->buf
    + (krylov_size-1)*(n_rows+krylov_size);
  _beta           = (double *)_sla_itsol_aux_buffer->buf
    + (krylov_size-1)*(n_rows+krylov_size) + 2*krylov_size;

  /* Initialization */
  for (ui = 0; ui < krylov_size*(krylov_size-1); ui++)
    _H_matrix[ui] = 0.;

  if (param.precond == CS_PARAM_PRECOND_DIAG ||
      param.precond == CS_PARAM_PRECOND_SSOR) {
    cs_sla_matrix_get_diag(a, &ad_inv);
    for (ii = 0; ii < n_rows; ii++) ad_inv[ii] = 1./ad_inv[ii];
  }

  BFT_MALLOC(M, 1, cs_sla_matrix_t);
  *M = *a;

  cvg = 0;
  while (cvg == 0) {

    /*  Compute  rk <- a*vx (vx = x0)
                 rk <- rhs - rk (r0=b-A*x0) */

    cs_sla_matvec(a, vx, &dk, true);
    for (ii = 0; ii < n_rows; ii++)
      dk[ii] = rhs[ii] - dk[ii];

    /* beta = ||r0|| */
    beta = sqrt(cs_dp(n_rows, dk, dk));
    scal_prod = beta;
    _beta[0] = beta;
    for (ui = 1; ui < krylov_size ; ui++)
      _beta[ui] = 0.;

    /* lap */
    l_iter = 0;
    l_old_iter = 0;

    for (ui = 0; ui < krylov_size - 1; ui++) {

      ret.iter = n_iter;

      /* krk = k ieme col of _krylov_vector = vi */
      krk = _krylov_vectors + ui*n_rows;

      for (jj = 0; jj < n_rows; jj++)
        krk[jj] = dk[jj] / scal_prod;

      _preconditionning(n_rows, M, krk, param.precond, gk);

      /* Compute w=dk <- A*vj */
      cs_sla_matvec(a, gk, &dk, true);

      for (uk = 0; uk < ui+1; uk++) {

        /* Compute h(k,i)=<w,vi>=<dk,vi> */
        _H_matrix[ui*krylov_size+uk] = cs_dp(n_rows, dk,
                                             (_krylov_vectors + uk*n_rows));

        /* Compute w=dk <- w - h(i,k)*vi */
        _daxpy(n_rows, -_H_matrix[ui*krylov_size + uk],
               (_krylov_vectors + uk*n_rows),
               1, dk, 1);

      } /* End of loop on uk */

      /* Compute h(i+1,i) = sqrt<w,w> */
      scal_prod = sqrt(cs_dp(n_rows, dk, dk));
      _H_matrix[ui*krylov_size + ui + 1] = scal_prod;

      if (scal_prod < epsi) scaltest = 1;

      if ((l_iter+1) % check_freq == 0 ||
          (size_t)l_iter == krylov_size - 2    ||
          scaltest == 1) {

        /* H matrix into diag sup matrix */
        _givens_rot_update(_H_matrix, krylov_size,
                           _beta, _givens_coeff,
                           l_old_iter, l_iter+1);

        /* solve diag sup system */
        _solve_diag_sup_halo(_H_matrix, l_iter+1, krylov_size, _beta ,gk);

        l_old_iter = l_iter+1;

        for (jj = 0; jj < n_rows; jj++)
          fk[jj]=_ddot(l_iter+1, _krylov_vectors + jj, n_rows, gk, 1);

        _preconditionning(n_rows, M, fk, param.precond, gk);

        for (jj = 0; jj < n_rows; jj++)
          fk[jj] = vx[jj] + gk[jj];

        cs_sla_matvec(a, fk, &bk, true);

        /* Compute residual = | Ax-b |_1 */
        residue=0.;
        for (jj = 0; jj < n_rows; jj++)
          bk[jj] = rhs[jj] - bk[jj];

        if (!_is_cvg(param, bk, fk, &ret, &cvgm)) cvg = 1;

      } /* Endif test convergence */

      n_iter++;
      l_iter++;

      if (cvg == 1 || (size_t)l_iter == krylov_size - 1 || scaltest == 1) {
        for (jj = 0; jj < n_rows; jj++)
          vx[jj] = fk[jj];
        break;
      }

    } /* End of loop on Krylov size */

  } /* End of while (cvg) */

  BFT_FREE(M);

  printf("FIN GMRES\n"
         "cvg: %d, n_iter = %d, lap_iter = %d, kylov_size =  %lu "
         " problem size = %d\n",
         cvg, n_iter-1, l_iter, krylov_size - 1, n_rows);

  return ret;
}

/*---------------------------------------------------------------------------*
 * Bi-Conjuguate gradient
 *---------------------------------------------------------------------------*/

static cs_sla_sumup_t
_bicgstab(const cs_sla_matrix_t    *A,
          const double              b[],
          double                    x[],    /* in/out */
          cs_param_itsol_t          param)
{
  int  i;
  double  alpha, beta, omega, rho0, rho1, t;

  double  *restrict rk = NULL, *r0 = NULL;
  double  *restrict pk = NULL, *Apk = NULL;
  double  *restrict sk = NULL, *Ask = NULL;

  cs_sla_sumup_t  ret = {CS_SLA_CODE_CONTINUE, 0, DBL_MAX};
  _cvg_monitor_t  cvgm = _init_solver_info(A, b);

  const int  size = A->n_rows;

  rk  = (double *)_sla_itsol_main_buffer->buf;
  r0  = (double *)_sla_itsol_main_buffer->buf + size;
  pk  = (double *)_sla_itsol_main_buffer->buf + 2*size;
  Apk = (double *)_sla_itsol_main_buffer->buf + 3*size;
  sk  = (double *)_sla_itsol_main_buffer->buf + 4*size;
  Ask = (double *)_sla_itsol_main_buffer->buf + 5*size;

  /* Compute initial residual: rk = b - Ax */
  cs_sla_amxby(-1.0, A, x, 1.0, b, &r0);

  memcpy(rk, r0, size*sizeof(double));
  memcpy(pk, r0, size*sizeof(double));

  rho0 = cs_dp(size, rk, r0);

  /* Loops until convergence */
  while (_is_cvg(param, rk, x, &ret, &cvgm)) {

    ret.iter += 1;

    /* Compute alpha */
    cs_sla_matvec(A, pk, &Apk, true);
    t = cs_dp(size, Apk, r0);
    assert(fabs(t) > DBL_MIN);
    alpha = rho0 / t;

    /* Compute sk */
    for (i = 0; i < size; i++) sk[i] = rk[i] - alpha * Apk[i];

    /* Compute omega */
    cs_sla_matvec(A, sk, &Ask, true);
    t = cs_dp(size, Ask, Ask);
    assert(fabs(t) > DBL_MIN);
    omega = cs_dp(size, Ask, sk) / t;

    /* Update vectors */
    for (i = 0; i < size; i++) {
      x[i] += alpha * pk[i] + omega * sk[i];
      rk[i] = sk[i] - omega * Ask[i];
    }

    /* Compute beta */
    rho1 = cs_dp(size, rk, r0);
    assert(fabs(rho0) > DBL_MIN);
    beta = (alpha * rho1) /(omega * rho0);
    rho0 = rho1;

    /* Update p */
    for (i = 0; i < size; i++)
      pk[i] = rk[i] + beta * (pk[i] - omega * Apk[i]);

  } /* End of while */

  return ret;
}

/*---------------------------------------------------------------------------*
 * Conjuguate gradient (With conditionning computation). No preconditionning
 *---------------------------------------------------------------------------*/

static cs_sla_sumup_t
_cg2(const cs_sla_matrix_t  *A,
     const double           *b,
     double                  x[],    /* in/out */
     cs_param_itsol_t        param)
{
  int  i;
  double  alphak, betak, Apkpk, rkrk0, rkrk1;
  double  lambda_min, lambda_max, kappa;

  double  *rk = NULL, *restrict pk = NULL, *Apk = NULL;

  cs_sla_sumup_t  ret = {CS_SLA_CODE_CONTINUE, 0, DBL_MAX};
  _cvg_monitor_t  cvgm = _init_solver_info(A, b);

  const int  size = A->n_rows; /* symmetric => n_rows = n_cols */

  /* Assign memory */
  rk  = (double *)_sla_itsol_main_buffer->buf;
  pk  = (double *)_sla_itsol_main_buffer->buf + size;
  Apk = (double *)_sla_itsol_main_buffer->buf + 2*size;

  /* Initialize algo. and compute initial residual: rk = b - A*x */
  cs_sla_amxby(-1.0, A, x, 1.0, b, &rk);

  memcpy(pk, rk, size*sizeof(double));

  rkrk0 = cs_dp(size, rk, rk);
  lambda_min = DBL_MAX;
  lambda_max = -DBL_MAX;
  kappa = 1;

  /* Loops until convergence */
  while (_is_cvg(param, rk, x, &ret, &cvgm)) {

    ret.iter += 1;

    cs_sla_matvec(A, pk, &Apk, true);
    Apkpk = cs_dp(size, Apk, pk);
    if (Apkpk > 0)
      alphak = rkrk0 / Apkpk;
    else {
      alphak = 0.0;
      bft_error(__FILE__, __LINE__, 0,
                " The matrix is not Symmetric Positivite Definite (SPD)\n"
                " < pk, pk >_A = %g\n"
                " CG algorithm cannot converged with this kind of matrix.\n"
                " Please choose another solver.", Apkpk);
    }

    /* Update vectors */
    for (i = 0; i < size; i++) {
      x[i] = x[i] + alphak * pk[i];
      rk[i] = rk[i] - alphak * Apk[i];
    }

    rkrk1 = cs_dp(size, rk, rk);
    betak = rkrk1 / rkrk0;
    rkrk0 = rkrk1;

    for (i = 0; i < size; i++)
      pk[i] = rk[i] + betak * pk[i];

    lambda_min = CS_MIN(lambda_min, alphak);
    lambda_max = CS_MAX(lambda_max, alphak);

  } /* End of while */

#if SLA_ITSOL_DEBUG > 0
  kappa = lambda_max / lambda_min;
  printf("lambda_min: %8.5e, lambda_max: %8.5e, kappa: %8.5e\n",
         lambda_min, lambda_max, kappa);
#endif

  return ret;
}

/*---------------------------------------------------------------------------*
 * Conjuguate gradient: no preconditionning                                  *
 *---------------------------------------------------------------------------*/

static cs_sla_sumup_t
_cg(const cs_sla_matrix_t  *A,
    const double           *b,
    double                 *x,    /* in/out */
    cs_param_itsol_t        param)

{
  int  i;
  double  alphak, betak, rkrk, Apkpk;

  double  *rk = NULL, *restrict pk = NULL, *Apk = NULL;

  cs_sla_sumup_t  ret = {CS_SLA_CODE_CONTINUE, 0, DBL_MAX};
  _cvg_monitor_t  cvgm = _init_solver_info(A, b);

  const int  size = A->n_rows; /* symmetric => n_rows = n_cols */

  /* Assign memory  */
  rk  = (double *)_sla_itsol_main_buffer->buf;
  pk  = (double *)_sla_itsol_main_buffer->buf + size;
  Apk = (double *)_sla_itsol_main_buffer->buf + 2*size;

  /* Initialize algo. and compute initial residual: rk = b - A*x */
  cs_sla_amxby(-1.0, A, x, 1.0, b, &rk);

  memcpy(pk, rk, size*sizeof(double));

  /* Loops until convergence */
  while (_is_cvg(param, rk, x, &ret, &cvgm)) {

    ret.iter += 1;

    rkrk = cs_dp(size, rk, rk);
    cs_sla_matvec(A, pk, &Apk, true);
    Apkpk = cs_dp(size, Apk, pk);
    if (Apkpk > 0)
      alphak = rkrk / Apkpk;
    else
      bft_error(__FILE__, __LINE__, 0,
                " The matrix is not Symmetric Positivite Definite (SPD)\n"
                " < pk, pk >_A = %g\n"
                " CG algorithm cannot converged with this kind of matrix.\n"
                " Please choose another solver.", Apkpk);

    /* Update vectors */
    for (i = 0; i < size; i++) {
      x[i] = x[i] + alphak * pk[i];
      rk[i] = rk[i] - alphak * Apk[i];
    }

    betak = cs_dp(size, rk, rk) / rkrk;
    for (i = 0; i < size; i++)
      pk[i] = rk[i] + betak * pk[i];

  } /* End of while */

  return ret;
}

/*---------------------------------------------------------------------------*
 * Solve Ax = b
 * with a BiConjuguate gradient algorithm with a preconditionning
 *---------------------------------------------------------------------------*/

static cs_sla_sumup_t
_pcd_bicgstab(const cs_sla_matrix_t  *A,
              const double           *b,
              double                 *x,    /* in/out */
              cs_param_itsol_t        param)
{
  int  i;
  double  alpha, beta, omega, rho0, rho1, t;

  double  *rk = NULL, *res0 = NULL, *pk = NULL, *Apk = NULL;
  double  *sk = NULL, *Ask = NULL, *ppk = NULL, *spk = NULL;
  double  *pdiag  = NULL;

  cs_sla_sumup_t  ret = {CS_SLA_CODE_CONTINUE, 0, DBL_MAX};
  _cvg_monitor_t  cvgm = _init_solver_info(A, b);

  const int size = A->n_rows;

  pdiag = (double *)_sla_itsol_pcd_buffer->buf;
  res0 = (double *)_sla_itsol_main_buffer->buf;
  rk   = (double *)_sla_itsol_main_buffer->buf + size;
  pk   = (double *)_sla_itsol_main_buffer->buf + 2*size;
  Apk  = (double *)_sla_itsol_main_buffer->buf + 3*size;
  sk   = (double *)_sla_itsol_main_buffer->buf + 4*size;
  Ask  = (double *)_sla_itsol_main_buffer->buf + 5*size;
  ppk  = (double *)_sla_itsol_main_buffer->buf + 6*size;
  spk  = (double *)_sla_itsol_main_buffer->buf + 7*size;

  /* Compute initial residual: rk = b - Ax */
  cs_sla_amxby(-1.0, A, x, 1.0, b, &res0);

  /* Initialization */
  memcpy(rk, res0, size*sizeof(double));
  memcpy(pk, res0, size*sizeof(double));

  rho1 = rho0 = 1.0;
  alpha = 0.0;
  omega = 1.0;
  for (i = 0; i < size; i++)  /* In order to avoid valgrind warnings */
    Apk[i] = 0.0;

  while(_is_cvg(param, rk, x, &ret, &cvgm)) {

    ret.iter++;

    /* Compute rho and beta */
    rho1 = cs_dp(size, rk, res0);
    assert(fabs(rho0) > DBL_MIN);
    assert(fabs(omega) > DBL_MIN);
    beta = (alpha * rho1) /(omega * rho0);
    rho0 = rho1;

    /* Compute pk */
    for (i = 0; i < size; i++)
      pk[i] = rk[i] + beta * (pk[i] - omega * Apk[i]);

    /* Compute alpha */
    _preconditionning(size, A, pk, param.precond, ppk);
    cs_sla_matvec(A, ppk, &Apk, true);
    t = cs_dp(size, Apk, res0);
    assert(fabs(t) > DBL_MIN);
    alpha = rho1 / t;

    /* Compute sk */
    for (i = 0; i < size; i++)
      sk[i] = rk[i] - alpha * Apk[i];

    /* Compute omega */
    _preconditionning(size, A, sk, param.precond, spk);
    cs_sla_matvec(A, spk, &Ask, true);
    t = cs_dp(size, Ask, Ask);
    assert(fabs(t) > DBL_MIN);
    omega = cs_dp(size, Ask, sk) / t;

    /* Update vectors */
    for (i = 0; i < size; i++) {
      x[i] += alpha * ppk[i] + omega * spk[i];
      rk[i] = sk[i] - omega * Ask[i];
    }

  } /* End of while */

  return ret;
}

/*---------------------------------------------------------------------------*
 * Conjuguate gradient with preconditionning techniques
 *---------------------------------------------------------------------------*/

static cs_sla_sumup_t
_pcd_cg(const cs_sla_matrix_t   *A,
        const double            *b,
        double                  *x,   /* in/out */
        cs_param_itsol_t         param)
{
  int  i;
  double  alphak = 0., betak, rkzk, Apkpk;

  cs_sla_sumup_t  ret = {CS_SLA_CODE_CONTINUE, 0, DBL_MAX};
  _cvg_monitor_t  cvgm = _init_solver_info(A, b);

  const int  size = A->n_rows;

  double  *rk  = (double *)_sla_itsol_main_buffer->buf;
  double  *zk  = (double *)_sla_itsol_main_buffer->buf + size;
  double  *pk  = (double *)_sla_itsol_main_buffer->buf + 2*size;
  double  *Apk = (double *)_sla_itsol_main_buffer->buf + 3*size;

  /* Compute initial residual: rk = b - Ax */
  cs_sla_amxby(-1.0, A, x, 1.0, b, &rk);

  /* Initialize */
  _preconditionning(size, A, rk, param.precond, zk);
  memcpy(pk, zk, size*sizeof(double));

  /* Loops until convergence */
  while (_is_cvg(param, rk, x, &ret, &cvgm)) {

    ret.iter += 1;

    /* Compute alphak */
    rkzk = cs_dp(size, rk, zk);
    cs_sla_matvec(A, pk, &Apk, true);
    Apkpk = cs_dp(size, Apk, pk);

    if (Apkpk > 0.0)
      alphak = rkzk / Apkpk;
    else
      bft_error(__FILE__, __LINE__, 0,
                " The matrix is not Symmetric Positivite Definite (SPD)\n"
                " < pk, pk >_A = %g\n"
                " CG algorithm cannot converged with this kind of matrix.\n"
                " Please choose another solver.", Apkpk);

    /* Update vectors */
    for (i = 0; i < size; i++) {
      x[i] += alphak * pk[i];
      rk[i] -= alphak * Apk[i];
    }
    _preconditionning(size, A, rk, param.precond, zk);
    betak = cs_dp(size, rk, zk) / rkzk;
    for (i = 0; i < size; i++)
      pk[i] = zk[i] + betak * pk[i];

  } /* End of while */

  return ret;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a sparse linear system A*x = rhs
 *         Matrix A can re-order during the solving process.
 *
 * \param[in]     itsol_info   set of parameters for the iterative solver
 * \param[in]     A            matrix to invert
 * \param[in]     rhs          right hand side
 * \param[inout]  x            initial guess (in) / solution (out)
 *
 * \return  a summary of the convergence monitoring
 */
/*----------------------------------------------------------------------------*/

cs_sla_sumup_t
cs_sla_solve(const cs_param_itsol_t    itsol_info,
             const cs_sla_matrix_t    *A,
             const double              rhs[],
             double                   *x)
{
  int  i, nx;

  cs_sla_sumup_t  ret;

  /* Sanity checks */
  assert(A != NULL);
  assert(rhs != NULL);
  assert(x != NULL);
  assert(_sla_itsol_main_buffer->buf != NULL);

  nx = A->n_rows;

  switch(itsol_info.solver) {

  case CS_PARAM_ITSOL_CG:
    switch(itsol_info.precond) {

    case CS_PARAM_PRECOND_NONE:
      ret = _cg2(A, rhs, x, itsol_info);
      break;

    case CS_PARAM_PRECOND_CHEBY:
      assert(_sla_itsol_pcd_buffer->buf != NULL);
      ret = _pcd_cg(A, rhs, x, itsol_info);
      break;

    case CS_PARAM_PRECOND_DIAG:
    case CS_PARAM_PRECOND_SSOR:
      assert(_sla_itsol_pcd_buffer->buf != NULL);
      assert(A->flag & CS_SLA_MATRIX_LU);

      double  *pdiag = (double *)_sla_itsol_pcd_buffer->buf;

      /* Build diagonal inverse */
      cs_sla_matrix_get_diag(A, &pdiag);
      for (i = 0; i < nx; i++) {
        assert(fabs(pdiag[i]) > DBL_MIN);
        pdiag[i] = 1.0/pdiag[i];
      }

      ret = _pcd_cg(A, rhs, x, itsol_info);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " Invalid preconditionner. Stop Execution.\n");
    }
    break;

  case CS_PARAM_ITSOL_BICG:
    switch(itsol_info.precond) {
    case CS_PARAM_PRECOND_NONE:
      ret = _bicgstab(A, rhs, x, itsol_info);
      break;

    case CS_PARAM_PRECOND_CHEBY:
      assert(_sla_itsol_pcd_buffer->buf != NULL);
      ret = _pcd_bicgstab(A, rhs, x, itsol_info);
      break;

    case CS_PARAM_PRECOND_SSOR:
    case CS_PARAM_PRECOND_DIAG:
      assert(A->flag & CS_SLA_MATRIX_LU);
      assert(_sla_itsol_pcd_buffer->buf != NULL);
      double  *pdiag = (double *)_sla_itsol_pcd_buffer->buf;

      /* Build diagonal inverse */
      cs_sla_matrix_get_diag(A, &pdiag);
      for (i = 0; i < nx; i++) {
        assert(fabs(pdiag[i]) > DBL_MIN);
        pdiag[i] = 1.0/pdiag[i];
      }
      ret = _pcd_bicgstab(A, rhs, x, itsol_info);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " Invalid preconditionner. Stop Execution.\n");
    }
    break;

  case CS_PARAM_ITSOL_GMRES:
    switch(itsol_info.precond) {
    case CS_PARAM_PRECOND_NONE:
    case CS_PARAM_PRECOND_SSOR:
    case CS_PARAM_PRECOND_DIAG:
    case CS_PARAM_PRECOND_CHEBY:
      assert(A->flag & CS_SLA_MATRIX_LU);
      ret = _gmres(A, rhs, x, itsol_info);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " Invalid preconditionner. Stop Execution.\n");
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Invalid solver. Stop Execution.\n");

  } /* End of switch solver */

  return ret;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Evaluate the sizes of temporary buffers used in iterative solvers
 *
 * \param[in]  refsize        reference size to evaluate buffer dimension
 * \param[in]  info           information about the iterative solver to use
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_itsol_update_buffer_sizes(size_t                   refsize,
                                 const cs_param_itsol_t   info)
{
  size_t  _main_size = 0, _pcd_size = 0, _aux_size = 0;

  switch(info.solver) {

  case CS_PARAM_ITSOL_CG:
    switch(info.precond) {
    case CS_PARAM_PRECOND_NONE:
      _main_size = 3 * sizeof(double) * refsize;
      break;

    case CS_PARAM_PRECOND_CHEBY:
    case CS_PARAM_PRECOND_DIAG:
    case CS_PARAM_PRECOND_SSOR:
      _main_size = 4 * sizeof(double) * refsize;
      _pcd_size = sizeof(double) * refsize;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " Unknown preconditionning type. Stop Execution.\n");
    }
    break;

  case CS_PARAM_ITSOL_BICG:
    switch(info.precond) {
    case CS_PARAM_PRECOND_NONE:
      _main_size = 6 * sizeof(double) * refsize;
      break;

    case CS_PARAM_PRECOND_CHEBY:
    case CS_PARAM_PRECOND_DIAG:
    case CS_PARAM_PRECOND_SSOR:
      _main_size = 8 * sizeof(double) * refsize;
      _pcd_size = sizeof(double) * refsize;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " Unknown preconditionning type. Stop Execution.\n");
    }
    break;

  case CS_PARAM_ITSOL_GMRES:
    {
      size_t  auxsize, krylov_size;

      _get_krylov_sizes(refsize, &auxsize, &krylov_size);
      _main_size = 4 * sizeof(double) * refsize;
      _aux_size = sizeof(double) * auxsize;
      _pcd_size = sizeof(double) * refsize;
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Unknown solver type. Stop Execution.\n");

  } /* End of switch on solver type */

  /* Update */
  _sla_itsol_main_size = CS_MAX(_sla_itsol_main_size, _main_size);
  _sla_itsol_pcd_size = CS_MAX(_sla_itsol_pcd_size, _pcd_size);
  _sla_itsol_aux_size = CS_MAX(_sla_itsol_aux_size, _aux_size);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate temporary buffers used in iterative solvers
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_itsol_alloc_buffers(void)
{
  cs_tmpbuf_alloc(_sla_itsol_main_size, &_sla_itsol_main_buffer);
  cs_tmpbuf_alloc(_sla_itsol_aux_size,  &_sla_itsol_aux_buffer);
  cs_tmpbuf_alloc(_sla_itsol_pcd_size,  &_sla_itsol_pcd_buffer);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free temporary buffers used in iterative solvers
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_itsol_free_buffers(void)
{
  _sla_itsol_main_buffer = cs_tmpbuf_free(_sla_itsol_main_buffer);
  _sla_itsol_aux_buffer = cs_tmpbuf_free(_sla_itsol_aux_buffer);
  _sla_itsol_pcd_buffer = cs_tmpbuf_free(_sla_itsol_pcd_buffer);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of a monitoring code
 *
 * \param[in] code     type of code
 *
 * \return the associated code name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_sla_get_code_name(cs_sla_code_t  code)
{
  switch (code) {
  case CS_SLA_CODE_CVG:
  case CS_SLA_CODE_CVG_EPS:
  case CS_SLA_CODE_STOP:
    return  "Convergence";
    break;

  case CS_SLA_CODE_STAG:
    return  "Stagnation";
    break;

  case CS_SLA_CODE_DVG:
    return  "Divergence";
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid monitoring code. Stop execution."));
  }

  return "";
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
