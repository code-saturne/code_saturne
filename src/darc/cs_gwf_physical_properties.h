#ifndef __CS_GWF_PHYSICAL_PROPERTIES_H__
#define __CS_GWF_PHYSICAL_PROPERTIES_H__

/*============================================================================
 * Physical properties management for groundwater flow module.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
 * Update physical properties of the Ground Water Flow module.
 *
 * Species transport is delayed by retention in solid phase.
 * This delay is computed as follows:
 * R = 1 + rho K_d/theta ;
 * where R is the delay factor, rho the soil density,
 * K_d the contaminant distribution coefficient and theta
 * the moisture content (saturation).
 *
 *---------------------------------------------------------------------------*/

void
cs_gwf_physical_properties(void);

/*----------------------------------------------------------------------------
 * Update sorbed concentration for scalars with kinetic sorption.
 *
 * It is estimated by the following analytical expression :
 * S^{n+1} = S^n exp(- k^{-} dt) - C^n * k^{+}/k^{-}
 * (exp(- k^{-} dt) - 1)
 *
 * parameters:
 * c_scal -->   concentration field
 * kp     -->   kplus field
 * km     -->   kminus field
 * sorb   <--   sorbed concentration field
 *----------------------------------------------------------------------------*/

void
cs_gwf_sorbed_concentration_update(cs_real_t *c_scal,
                                   cs_real_t *kp,
                                   cs_real_t *km,
                                   cs_real_t *sorb);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_GWF_PHYSICAL_PROPERTIES_H__ */
