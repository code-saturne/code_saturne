/*============================================================================
 * Set of pratical functions dedicated to the groundwater flow module
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gwf_toolbox.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the values of the inventory for each tracer at the time
 *        time_eval given as parameter from an initial inventory. The evolution
 *        of the inventory follows the Bateman's solution without source term.
 *
 *        The decay coefficient is automatically retrieved from the data
 *        settings given when creating the decay chain. One assumes that the
 *        inventory values are given in mol. The time unit (s, hour, year...)
 *        has to be consistent with the unit given as parameter at the
 *        definition of the decay chain and with the value of time_eval.
 *
 * \param[in]       time_eval   time at which one evaluates the new inventory
 * \param[in]       tdc         decay chain to consider
 * \param[in]       init_inv    initial inventory for each tracer in the chain
 * \param[in, out]  eval_inv    resulting inventory at the given time
 */
/*----------------------------------------------------------------------------*/

void
cs_gwf_toolbox_bateman(double                              time_eval,
                       const cs_gwf_tracer_decay_chain_t  *tdc,
                       const double                        init_inv[],
                       double                              eval_inv[])
{
  if (tdc == nullptr)
    return;
  if (tdc->n_tracers < 1)
    return;

  /* Collect all decay coefficients in one array (either in the stack or malloc
     if there are more than 6 tracers in the decay chain */

  double *lambda = nullptr, *_lambda = nullptr;
  double  lambda6[6];

  if (tdc->n_tracers < 6)
    lambda = lambda6;
  else {
    BFT_MALLOC(_lambda, tdc->n_tracers, double);
    lambda = _lambda;
  }

  for (int it = 0; it < tdc->n_tracers; it++) {
    cs_gwf_tracer_default_context_t *tc
      = (cs_gwf_tracer_default_context_t *)tdc->tracers[it]->context;
    lambda[it] = tc->decay_coef;
  }

  for (int it = 0; it < tdc->n_tracers; it++) {

    eval_inv[it] = 0.; /* Reset the result array */

    for (int i = 0; i < it+1; i++) {

      double  a = 1.;
      for (int l = i; l < it; l++) a *= lambda[l];

      double  b = 0.;
      for (int j = i; j < it+1; j++) {

        double  c = 1.;
        for (int k = i; k < it+1; k++)
          if (k != j)
            c *= lambda[k]-lambda[j];

        b += exp(-lambda[j]*time_eval)/c;

      } /* Loop on j */

      eval_inv[it] += init_inv[i]*a*b;

    } /* Loop on i */

  } /* Loop on tracers */

  if (lambda != lambda6)
    BFT_FREE(_lambda);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
