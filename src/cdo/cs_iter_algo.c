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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_cdo_sqnorm.h"

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

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_ITER_ALGO_DBG      0

/*============================================================================
 * Private variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/


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
cs_iter_algo_define(int          verbosity,
                    int          n_max_iter,
                    double       atol,
                    double       rtol,
                    double       dtol)
{
  cs_iter_algo_info_t  *info = NULL;

  BFT_MALLOC(info, 1, cs_iter_algo_info_t);

  info->verbosity = verbosity;
  info->atol = atol;
  info->rtol = rtol;
  info->dtol = dtol;
  info->n_max_algo_iter = n_max_iter;
  info->normalization = 1.0;

  cs_iter_algo_reset(info);

  return info;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if something wrong happens during the iterative process
 *
 * \param[in] func_name    name of the calling function
 * \param[in] eq_name      name of the equation being solved
 * \param[in] algo_name    name of the iterative algo. used
 * \param[in] iai          pointer to the iterative algo. structure
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_check(const char            *func_name,
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
 * \brief  Test if one has to do one more Picard iteration.
 *         Test if performed on the relative norm on the increment between
 *         two iterations but also on the divergence.
 *
 * \param[in]      pre_iterate    previous state of the mass flux iterate
 * \param[in]      cur_iterate    current state of the mass flux iterate
 * \param[in]      div_l2_norm    L2 norm of the velocity divergence
 * \param[in, out] iai            pointer to a cs_iter_algo_info_t structure
 *
 * \return the convergence state
 */
/*----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_iter_algo_navsto_fb_picard_cvg(const cs_real_t             *pre_iterate,
                                  const cs_real_t             *cur_iterate,
                                  cs_real_t                    div_l2_norm,
                                  cs_iter_algo_info_t         *iai)
{
  cs_real_t  previous_res;

  /* Set the tolerance criterion (one at each call if the normalization is
     modified between two successive calls) */

  iai->tol = fmax(iai->rtol*iai->normalization, iai->atol);

  /* Compute the norm of the difference between the two mass fluxes (the
     current one and the previous one) */

  if (iai->n_algo_iter == 0) {

    /* Store the first residual to detect a divergence */

    iai->res0 = cs_cdo_sqnorm_pfsf_diff(pre_iterate, cur_iterate);
    assert(iai->res0 > -DBL_MIN);
    iai->res0 = sqrt(iai->res0);
    previous_res = iai->res0;
    iai->res = iai->res0;

  }
  else {

    /* Compute the norm of the difference between the two mass fluxes (the
       current one and the previous one) */

    previous_res = iai->res;
    iai->res = cs_cdo_sqnorm_pfsf_diff(pre_iterate, cur_iterate);
    assert(iai->res > -DBL_MIN);
    iai->res = sqrt(iai->res);

  }

  /* Increment the number of Picard iterations */

  iai->n_algo_iter += 1;

  /* Set the convergence status */

  if (iai->res < iai->tol)
    iai->cvg = CS_SLES_CONVERGED;

  else if (iai->n_algo_iter >= iai->n_max_algo_iter)
    iai->cvg = CS_SLES_MAX_ITERATION;

  else if (iai->res > iai->dtol*previous_res || iai->res > iai->dtol*iai->res0)
    iai->cvg = CS_SLES_DIVERGED;

  else
    iai->cvg = CS_SLES_ITERATING;

  if (iai->verbosity > 0) {
    if (iai->n_algo_iter == 1)
      cs_iter_algo_navsto_print_header("## Picard");
    cs_iter_algo_navsto_print("## Picard", iai, div_l2_norm);
  }

  return iai->cvg;
}

/*----------------------------------------------------------------------------*/


END_C_DECLS
