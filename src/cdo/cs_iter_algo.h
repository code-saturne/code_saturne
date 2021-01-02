#ifndef __CS_ITER_ALGO_H__
#define __CS_ITER_ALGO_H__

/*============================================================================
 * Routines to handle iterative algorithm such augmented Lagrangian, Picard or
 * Anderson algorithms
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"
#include "cs_math.h"
#include "cs_sles.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \struct cs_iter_algo_info_t
 *  \brief Set of information related to the convergence of the iterative
 *         algorithm (Picard or Uzawa for instance)
 *
 * \var verbosity
 * Level of printed information
 *
 * \var atol
 * absolute tolerance
 *
 * \var rtol
 * relative tolerance
 *
 * \var dtol
 * tolerance to detect a divergence of the algorithm. Not used if < 0
 *
 * \var cvg
 * converged, iterating or diverged status
 *
 * \var res
 * value of the residual for the iterative algorithm
 *
 * \var res0
 * Initial value of the residual for the iterative algorithm
 *
 * \var tol
 * tolerance computed as tol = max(atol, res0*rtol) where
 * atol and rtol are respectively the absolute and relative tolerance associated
 * to the algorithm
 *
 * \var n_algo_iter
 * number of iterations for the algorithm (outer iterations)
 *
 * \var n_max_algo_iter
 * maximal number of iterations for the algorithm
 *
 * \var n_inner_iter
 * cumulated number of inner iterations (sum over the outer iterations)
 *
 * \var last_inner_iter
 * last number of iterations for the inner solver
 */

typedef struct {

  int                              verbosity;

  /* Set of tolerances to drive the convergence of the iterative algorithm */
  double                           atol;
  double                           rtol;
  double                           dtol;

  /* Variable convergence indicators */
  cs_sles_convergence_state_t      cvg;
  double                           res;
  double                           res0;
  double                           tol;

  int                              n_algo_iter;
  int                              n_max_algo_iter;
  int                              n_inner_iter;
  int                              last_inner_iter;

} cs_iter_algo_info_t;

/*============================================================================
 * Inline static public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reset a cs_iter_algo_info_t structure
 *
 * \param[in, out]  info   pointer to a cs_iter_algo_info_t
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_iter_algo_reset(cs_iter_algo_info_t    *info)
{
  if (info == NULL)
    return;

  info->cvg = CS_SLES_ITERATING;
  info->res = cs_math_big_r;
  info->n_algo_iter = 0;
  info->n_inner_iter = 0;
  info->last_inner_iter = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print header before dumping information gathered in the structure
 *         cs_iter_algo_info_t
 *
 * \param[in]  algo_name     name of the algorithm
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_iter_algo_navsto_print_header(const char   *algo_name)
{
  assert(algo_name != NULL);
  cs_log_printf(CS_LOG_DEFAULT,
                "%12s.It  -- Algo.Res   Inner  Cumul  ||div(u)||  Tolerance\n",
                algo_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print header before dumping information gathered in the structure
 *         cs_iter_algo_info_t
 *
 * \param[in]  algo_name     name of the algorithm
 * \param[in]  info          pointer to cs_iter_algo_info_t structure
 * \param[in]  div_l2        l2 norm of the divergence
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_iter_algo_navsto_print(const char                    *algo_name,
                          const cs_iter_algo_info_t     *info,
                          double                         div_l2)
{
  assert(algo_name != NULL);
  cs_log_printf(CS_LOG_DEFAULT,
                "%12s.It%02d-- %5.3e  %5d  %5d  %6.4e  %6.4e\n",
                algo_name, info->n_algo_iter, info->res,
                info->last_inner_iter, info->n_inner_iter, div_l2, info->tol);
  cs_log_printf_flush(CS_LOG_DEFAULT);
}

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
                    double       dtol);

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
                                  cs_iter_algo_info_t         *a_info);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ITER_ALGO_H__ */
