#ifndef __CS_GWF_PHYSICAL_PROPERTIES_H__
#define __CS_GWF_PHYSICAL_PROPERTIES_H__

/*============================================================================
 * Physical properties management for groundwater flow module.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_gwf_parameters.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Update delay of all transported species (user scalars)
 *
 * Species transport is delayed by retention in solid phase. This delay is
 * computed as follows:
 *  R = 1 + rho K_d/theta
 *  where R is the delay factor, rho the soil density,
 * K_d the contaminant distribution coefficient and theta the moisture content
 * (saturation).
 *---------------------------------------------------------------------------*/

void
cs_gwf_delay_update(void);

/*----------------------------------------------------------------------------
 * Add first-order decay to implicit part of source term array
 *
 * Corresponding EDO for decay phenomenon is:
 * dc/dt = -decay_rate * c
 *
 * parameters:
 * f_id    --> field index of scalar on which decay is set
 * ts_imp  <-> implicit part of the source term
 *----------------------------------------------------------------------------*/

void
cs_gwf_decay_rate(const int        f_id,
                  cs_real_t       *ts_imp);

/*----------------------------------------------------------------------------
 * Update sorbed concentration for scalars with kinetic sorption.
 *
 * It is estimated by the following analytical expression :
 * S^{n+1} = S^n exp(- (k^{-} + decay_rate) * dt)
 *          - C^n * k^{+}/(k^{-} + decay_rate)
 *            (exp(- (k^{-} + decay_rate) * dt) - 1)
 *
 * parameters:
 * f_id    --> field index of scalar which properties are updated
 *----------------------------------------------------------------------------*/

void
cs_gwf_sorbed_concentration_update(const int f_id);

/*----------------------------------------------------------------------------
 * Clip liquid concentration to solubility index.
 *
 * parameters:
 * f_id    --> field index of scalar which properties are updated
 *----------------------------------------------------------------------------*/

void
cs_gwf_precipitation(const int f_id);

/*----------------------------------------------------------------------------
 * Take into account kinetic chemical reaction in evolution equation
 * of total liquid concentration.
 *
 * parameters:
 * f_id    --> field index of scalar which properties are updated
 * ts_imp  <-> implicit part of the source term
 * ts_exp  <-> explicit part of the source term
 *----------------------------------------------------------------------------*/

void
cs_gwf_kinetic_reaction(const int   f_id,
                        cs_real_t  *ts_imp,
                        cs_real_t  *ts_exp);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GWF_PHYSICAL_PROPERTIES_H__ */
