#ifndef CS_ITER_ALGO_H
#define CS_ITER_ALGO_H

/*============================================================================
 * Set of functions to manage high-level iterative algorithms
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "alge/cs_param_sles.h"
#include "alge/cs_sles.h"
#include "base/cs_math.h"
#include "base/cs_param_types.h"
#include "cdo/cs_cdo_blas.h"
#include "cdo/cs_sdm.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* \enum cs_iter_algo_type_bit_t
 * \brief Bit values for the definition of the type of iterative algorithms
 *
 * \var CS_ITER_ALGO_DEFAULT
 * Iterative algorithms like linear solvers for saddle-point system, fixed
 * point...
 *
 * \var CS_ITER_ALGO_ANDERSON
 * Iterative algorithms relying on the Anderson acceleration
 *
 * \var CS_ITER_ALGO_TWO_LEVEL
 * Iterative algorithms using an inner/outer iterative loops
 */

enum cs_iter_algo_type_bit_t {

  CS_ITER_ALGO_DEFAULT     = 1 << 0, /* =  1 */
  CS_ITER_ALGO_ANDERSON    = 1 << 1, /* =  2 */
  CS_ITER_ALGO_TWO_LEVEL   = 1 << 2, /* =  4 */

};

typedef cs_flag_t  cs_iter_algo_type_t;

/* Structure used as context for the iterative algorithm considered by
   default */

struct cs_iter_algo_default_t {

  /*!
   * \var cvg_status
   * Converged, iterating or diverged status
   *
   * \var normalization
   * Value of the normalization for the relative tolerance.
   *
   * The stopping criterion is such that res < rtol * normalization. By default,
   * the normalization is set to 1.
   *
   * \var tol
   * Tolerance computed as tol = max(atol, normalization*rtol) where atol and
   * rtol are respectively the absolute and relative tolerance associated to
   * the algorithm
   *
   * \var prev_res
   * Value of the previous residual achieved during the iterative process
   *
   * \var res
   * Value of the residual for the iterative algorithm
   *
   * \var res0
   * Value of the first residual of the iterative process. This is used for
   * detecting the divergence of the algorithm.
   *
   * \var n_algo_iter
   * Current number of iterations for the algorithm
   *
   * \var n_inner_iter
   * Curent cumulated number of inner iterations (sum over the outer iterations)
   *
   * \var last_inner_iter
   * Last number of iterations for the inner solver
   */

  cs_sles_convergence_state_t      cvg_status;
  double                           normalization;
  double                           tol;

  double                           prev_res;
  double                           res;
  double                           res0;

  int                              n_algo_iter;
  int                              n_inner_iter;
  int                              last_inner_iter;

};

/* Structures used when an Anderson acceleration is considered */
/* ----------------------------------------------------------- */

/*! \struct cs_iter_algo_param_aac_t
 *
 *  \brief Structure storing all the parameters to drive the algorithm called
 *  Anderson acceleration
 */

struct cs_iter_algo_param_aac_t {

/*!
 * \var n_max_dir
 * Maximum number of directions
 *
 * \var starting_iter
 * Anderson acceleration starts at this iteration number
 *
 * \var max_cond
 * Tolerance under which terms are dropped in order to improve the conditioning
 * number of the QR factorization. If < 0, there is no cleaning.
 *
 * \var beta
 * Value of the relaxation coefficient. Default value is equal to 1 meaning
 * that there is no relaxation.
 */

  int      n_max_dir;
  int      starting_iter;
  double   max_cond;
  double   beta;

};

/*! \struct cs_iter_algo_aac_t
 *  \brief Context structure for the algorithm called Anderson acceleration
 *
 *  Set of parameters and arrays to manage the Anderson acceleration
 */

struct cs_iter_algo_aac_t {

  /*!
   * @name Main parameters
   * @{
   *
   * \var param
   * Set of parameters driving the behavior of the Anderson acceleration
   */

  cs_iter_algo_param_aac_t        param;

  /*!
   * \var cvg_status
   * Converged, iterating or diverged status
   *
   * \var normalization
   * Value of the normalization for the relative tolerance.
   *
   * The stopping criterion is such that res < rtol * normalization. By default,
   * the normalization is set to 1.
   *
   * \var tol
   * Tolerance computed as tol = max(atol, normalization*rtol) where atol and
   * rtol are respectively the absolute and relative tolerance associated to
   * the algorithm
   *
   * \var prev_res
   * Value of the previous residual achieved during the iterative process
   *
   * \var res
   * Value of the residual for the iterative algorithm
   *
   * \var res0
   * Value of the first residual of the iterative process. This is used for
   * detecting the divergence of the algorithm.
   *
   * \var n_algo_iter
   * Current number of iterations for the algorithm (outer iterations)
   *
   * \var n_inner_iter
   * Curent cumulated number of inner iterations (sum over the outer iterations)
   *
   * \var last_inner_iter
   * Last number of iterations for the inner solver
   */

  cs_sles_convergence_state_t    cvg_status;
  double                         normalization;
  double                         tol;

  double                         prev_res;
  double                         res;
  double                         res0;

  int                            n_algo_iter;
  int                            n_inner_iter;
  int                            last_inner_iter;

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

  cs_lnum_t       n_elts; /*<! Size of an array for one direction */
  int             n_dir;

  cs_real_t      *fold;
  cs_real_t      *df;
  cs_real_t      *gold;
  cs_real_t      *dg;

  cs_real_t      *Q;
  cs_sdm_t       *R;

  cs_real_t      *gamma;

};

/* ++++++++++++++ */
/* Main structure */
/* ++++++++++++++ */

/*! \struct cs_iter_algo_t
 *  \brief Structure to handle the convergence of an iterative algorithm
 *
 *  Metadata to manage an iterative algorithm such as Picard or Uzawa for
 *  instance. This structure can handle a two-level (i.e. outer/inner)
 *  iterative algorithm since the notion of inner and outer iterations is
 *  defined. Nevertheless, only the outer iterative algorithm is managed
 *  (information about inner iterations are only for monitoring purposes).
 */

struct cs_iter_algo_t {

/*!
 * \var type
 * Type of iterative algorithm to consider. Useful to cast on-the-fly the
 * context structure
 *
 * \var cvg_param
 * structure storing the main settings
 *
 * \var verbosity
 * Level of printed information
 *
 * \var header_done
 * Is the header line printed ? true or false
 *
 * \var context
 * pointer to structure cast on the fly
 */

  cs_iter_algo_type_t       type;
  cs_param_convergence_t    cvg_param;
  int                       verbosity;
  bool                      header_done;

  void                     *context;

};

/*============================================================================
 * Public function prototypes
 *============================================================================*/

cs_iter_algo_t *
cs_iter_algo_create(cs_iter_algo_type_t type);

cs_iter_algo_t *
cs_iter_algo_create_with_settings(cs_iter_algo_type_t    type,
                                  int                    verbosity,
                                  cs_param_convergence_t cvg_param);

void
cs_iter_algo_free(cs_iter_algo_t **p_algo);

void
cs_iter_algo_reset(cs_iter_algo_t *algo);

void
cs_iter_algo_release_anderson_arrays(cs_iter_algo_aac_t *c);

void
cs_iter_algo_set_verbosity(cs_iter_algo_t *algo,
                           int             verbosity);

void
cs_iter_algo_set_cvg_param(cs_iter_algo_t         *algo,
                           cs_param_convergence_t  cvg_param);

void
cs_iter_algo_set_tolerance(cs_iter_algo_t *algo,
                           double          tol);

void
cs_iter_algo_set_initial_residual(cs_iter_algo_t *algo,
                                  double          value);

void
cs_iter_algo_set_normalization(cs_iter_algo_t *algo,
                               double          value);

void
cs_iter_algo_set_cvg_status(cs_iter_algo_t              *algo,
                            cs_sles_convergence_state_t  cvg_status);

void
cs_iter_algo_set_anderson_param(cs_iter_algo_t           *algo,
                                cs_iter_algo_param_aac_t  aac_param,
                                cs_lnum_t                 n_elts);

int
cs_iter_algo_get_n_iter(const cs_iter_algo_t *algo);

int
cs_iter_algo_get_n_inner_iter(const cs_iter_algo_t *algo);

double
cs_iter_algo_get_residual(const cs_iter_algo_t *algo);

cs_sles_convergence_state_t
cs_iter_algo_get_cvg_status(const cs_iter_algo_t *algo);

double
cs_iter_algo_get_normalization(const cs_iter_algo_t *algo);

cs_iter_algo_param_aac_t
cs_iter_algo_get_anderson_param(cs_iter_algo_t *algo);

void
cs_iter_algo_update_anderson(cs_iter_algo_t            *algo,
                             cs_real_t                 *cur_iterate,
                             const cs_real_t           *pre_iterate,
                             cs_cdo_blas_dotprod_t     *dotprod,
                             cs_cdo_blas_square_norm_t *sqnorm);

void
cs_iter_algo_update_inner_iters(cs_iter_algo_t *algo,
                                int             n_last_inner_iter);

void
cs_iter_algo_update_residual(cs_iter_algo_t *algo,
                             double          res);

cs_sles_convergence_state_t
cs_iter_algo_update_cvg_tol_given(cs_iter_algo_t *algo,
                                  double          tol);

cs_sles_convergence_state_t
cs_iter_algo_update_cvg_tol_auto(cs_iter_algo_t *algo);

void
cs_iter_algo_log_cvg(cs_iter_algo_t *algo,
                     const char     *label);
void
cs_iter_algo_check_warning(const char     *func_name,
                           const char     *eq_name,
                           const char     *algo_name,
                           cs_iter_algo_t *algo);

/*----------------------------------------------------------------------------*/

#endif /* CS_ITER_ALGO_H */
