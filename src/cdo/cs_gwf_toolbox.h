#ifndef __CS_GWF_TOOLBOX_H__
#define __CS_GWF_TOOLBOX_H__

/*============================================================================
 * Set of pratical functions dedicated to the groundwater flow module
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

#include "cs_defs.h"
#include "cs_gwf_tracer.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
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
                       double                              eval_inv[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GWF_TOOLBOX_H__ */
