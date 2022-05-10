#ifndef __CS_LAGR_AUX_MEAN_FLUID_QUANTITIES_H__
#define __CS_LAGR_AUX_MEAN_FLUID_QUANTITIES_H__

/*============================================================================
 * Functions and types for lagrangian field gradients
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute auxilary mean fluid quantities.
 *
 *  - Lagragian time
 *  - gradient of total pressure
 *  - velocity gradient
 *  - Lagragian time gradient
 *
 * \param[out]  lagr_time          Lagrangian time scale
 * \param[out]  grad_pr            pressure gradient
 * \param[out]  grad_vel           velocity gradient
 * \param[out]  grad_lagr_time     Lagrangian time gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_aux_mean_fluid_quantities(cs_field_t    *lagr_time,
                                  cs_real_3_t   *gradpr,
                                  cs_real_33_t  *grad_vel,
                                  cs_real_3_t   *grad_lagr_time);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_AUX_MEAN_FLUID_QUANTITIES_H__ */
