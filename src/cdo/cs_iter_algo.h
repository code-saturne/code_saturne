#ifndef __CS_ITER_ALGO_H__
#define __CS_ITER_ALGO_H__

/*============================================================================
 * Set of functions to manage high-level iterative algorithms
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdo_blas.h"
#include "cs_math.h"
#include "cs_param_sles.h"
#include "cs_param_types.h"
#include "cs_sles.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

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

typedef enum {

  CS_ITER_ALGO_DEFAULT             = 1 << 0, /* =  1 */
  CS_ITER_ALGO_ANDERSON            = 1 << 1, /* =  2 */
  CS_ITER_ALGO_TWO_LEVEL           = 1 << 2, /* =  4 */

} cs_iter_algo_type_bit_t;

typedef cs_flag_t  cs_iter_algo_type_t;

/* Structure used as context for the iterative algorithm considered by
   default */

typedef struct {

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

} cs_iter_algo_default_t;

/* Structures used when an Anderson acceleration is considered */
/* ----------------------------------------------------------- */

/*! \struct cs_iter_algo_param_aac_t
 *
 *  \brief Structure storing all the parameters to drive the algorithm called
 *  Anderson acceleration
 */

typedef struct {

/*!
 * \var n_max_dir
 * Maximum number of directions
 *
 * \var starting_iter
 * Anderson acceleration starts at this iteration number
 *
 * \var max_cond
 * Tolerance under which terms are dropped in order to improve the
 * conditionning number of the QR factorization
 *
 * \var beta
 * Value of the relaxation coefficient (if this is equal to zero then there is
 * non relaxation to perform)
 *
 * \var dp_type
 * Type of dot product to apply (usual Euclidean 2-norm or CDO-based one)
 */

  int                        n_max_dir;
  int                        starting_iter;
  double                     max_cond;
  double                     beta;
  cs_param_dotprod_type_t    dp_type;

} cs_iter_algo_param_aac_t;

/*! \struct cs_iter_algo_aac_t
 *  \brief Context structure for the algorithm called Anderson acceleration
 *
 *  Set of parameters and arrays to manage the Anderson acceleration
 */

typedef struct {

  /*!
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

} cs_iter_algo_aac_t;

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

typedef struct {

/*!
 * \var type
 * Type of iterative algorithm to consider. Useful to cast on-the-fly the
 * context structure
 *
 * \var cvg_param
 * structure storing the main settings
 *
 * \var context
 * pointer to structure cast on the fly
 *
 * \var verbosity
 * Level of printed information
 *
 */

  cs_iter_algo_type_t              type;
  cs_param_sles_cvg_t              cvg_param;
  int                              verbosity;

  void                            *context;

} cs_iter_algo_t;

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
cs_iter_algo_create(cs_iter_algo_type_t    type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a new cs_iter_algo_t structure with the given settings
 *
 * \param[in] type       type of iterative algorithm
 * \param[in] verbosity  level of information to print
 * \param[in] param      set of parameters driving the convergence of the
 *                       iterative algorithm
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_iter_algo_t *
cs_iter_algo_create_with_settings(cs_iter_algo_type_t    type,
                                  int                    verbosity,
                                  cs_param_sles_cvg_t    cvg_param);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free a cs_iter_algo_t structure
 *
 * \param[in, out] p_algo     double pointer on the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_free(cs_iter_algo_t   **p_algo);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reset a cs_iter_algo_t structure
 *
 * \param[in, out] algo           pointer to a cs_iter_algo_t
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_reset(cs_iter_algo_t      *algo);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reset a cs_iter_algo_t structure
 *
 * \param[in, out] algo           pointer to a cs_iter_algo_t
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_release_anderson_arrays(cs_iter_algo_aac_t      *c);

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
                           int                 verbosity);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the criteria related the convergence of the given iterative
 *        algorithm
 *
 * \param[in, out] algo       pointer to the structure to update
 * \param[in]      param      set of parameters driving the convergence of the
 *                            iterative algorithm
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_set_cvg_param(cs_iter_algo_t        *algo,
                           cs_param_sles_cvg_t    cvg_param);

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
                           double            tol);

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
                                  double            value);

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
                               double            value);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the convergence status of the given structure
 *
 * \param[in, out] algo     pointer to the structure to update
 * \param[in]      value    normalization to apply
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_set_cvg_status(cs_iter_algo_t                  *algo,
                            cs_sles_convergence_state_t      cvg_status);

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
                                cs_lnum_t                   n_elts);

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
cs_iter_algo_get_n_iter(const cs_iter_algo_t        *algo);

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
cs_iter_algo_get_n_inner_iter(const cs_iter_algo_t        *algo);

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
cs_iter_algo_get_residual(const cs_iter_algo_t        *algo);

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
cs_iter_algo_get_cvg_status(const cs_iter_algo_t        *algo);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the normalization to apply when computing the tolerance threshold
 *
 * \param[in] algo   pointer to a cs_iter_algo_t structure
 */
/*----------------------------------------------------------------------------*/

double
cs_iter_algo_get_normalization(const cs_iter_algo_t    *algo);

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
cs_iter_algo_get_anderson_param(cs_iter_algo_t         *algo);

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
                             cs_cdo_blas_square_norm_t    *sqnorm);

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
                                int               n_last_inner_iter);

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
                             double               res);

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
                                  double                  tol);

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
cs_iter_algo_update_cvg_tol_auto(cs_iter_algo_t         *algo);

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
                     const char          *label);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if something wrong happens during the iterative process after
 *        one new iteration
 *
 * \param[in] func_name    name of the calling function
 * \param[in] eq_name      name of the equation being solved
 * \param[in] algo_name    name of the iterative algo. used
 * \param[in] algo         pointer to the iterative algo. structure
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_check_warning(const char            *func_name,
                           const char            *eq_name,
                           const char            *algo_name,
                           cs_iter_algo_t        *algo);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ITER_ALGO_H__ */
