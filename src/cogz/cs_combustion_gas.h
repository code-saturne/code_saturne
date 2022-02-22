#ifndef __CS_COMBUSTION_GAS_H__
#define __CS_COMBUSTION_GAS_H__

/*============================================================================
 * Gas combustion model.
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_field.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \brief Compute molar and mass fractions of
 *         elementary species Ye, Xe (fuel, O2, CO2, H2O, N2) from
 *         global species Yg (fuel, oxidant, products)
 *
 * \param[in]     yg            global mass fraction
 * \param[out]    ye            elementary mass fraction
 * \param[out]    xe            elementary molar fraction
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_gas_yg2xye(const cs_real_t  yg[],
                         cs_real_t        ye[],
                         cs_real_t        xe[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_COMBUSTION_GAS_H__ */
