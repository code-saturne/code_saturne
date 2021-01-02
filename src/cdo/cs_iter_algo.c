/*============================================================================
 * Set of routines to handle the management of high-level iterative algorithms
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

#include <bft_mem.h>

#include "cs_evaluate.h"

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
 * \brief Set of routines to handle the management of high-level iterative
 *        algorithms such as Uzawa, Golub-Kahan Bi-orthogonalization or
 *        Picard algorithms which incorporates inner solvers
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

  cs_iter_algo_reset(info);

  return info;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if one has to do one more Picard iteration.
 *         Test if performed on the relative norm on the increment between
 *         two iterations but also on the divergence.
 *
 * \param[in]      connect        set of additional connectivities for CDO
 * \param[in]      quant          set of additional geometrical quantities
 * \param[in]      pre_iterate    previous state of the mass flux iterate
 * \param[in]      cur_iterate    current state of the mass flux iterate
 * \param[in]      div_l2_norm    L2 norm of the velocity divergence
 * \param[in, out] a_info         pointer to a cs_iter_algo_info_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_iter_algo_navsto_fb_picard_cvg(const cs_cdo_connect_t      *connect,
                                  const cs_cdo_quantities_t   *quant,
                                  const cs_real_t             *pre_iterate,
                                  const cs_real_t             *cur_iterate,
                                  cs_real_t                    div_l2_norm,
                                  cs_iter_algo_info_t         *a_info)
{
  const cs_real_t  pre_picard_res = a_info->res;

  /* Storage of the initial residual to build a relative tolerance */
  if (a_info->n_algo_iter == 0) {

    a_info->res0 = cs_evaluate_delta_square_wc2x_norm(pre_iterate,
                                                      cur_iterate,
                                                      connect->c2f,
                                                      quant->pvol_fc);
    assert(a_info->res0 > -DBL_MIN);
    a_info->res0 = sqrt(a_info->res0);
    a_info->res = a_info->res0;
    a_info->tol = fmax(a_info->rtol*a_info->res0, a_info->atol);

  }
  else {

    a_info->res = cs_evaluate_delta_square_wc2x_norm(pre_iterate,
                                                     cur_iterate,
                                                     connect->c2f,
                                                     quant->pvol_fc);
    assert(a_info->res > -DBL_MIN);
    a_info->res = sqrt(a_info->res);

  }

  /* Increment the number of Picard iterations */
  a_info->n_algo_iter += 1;

  /* Set the convergence status */
  if (a_info->res < a_info->tol)
    a_info->cvg = CS_SLES_CONVERGED;

  else if (a_info->n_algo_iter >= a_info->n_max_algo_iter)
    a_info->cvg = CS_SLES_MAX_ITERATION;

  else if (a_info->res > a_info->dtol * pre_picard_res ||
           a_info->res > a_info->dtol * a_info->res0)
    a_info->cvg = CS_SLES_DIVERGED;

  else
    a_info->cvg = CS_SLES_ITERATING;

  if (a_info->verbosity > 0) {
    if (a_info->n_algo_iter == 1)
      cs_iter_algo_navsto_print_header("## Picard");
    cs_iter_algo_navsto_print("## Picard", a_info, div_l2_norm);
  }

}

/*----------------------------------------------------------------------------*/


END_C_DECLS
