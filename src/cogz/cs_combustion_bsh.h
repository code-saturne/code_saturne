#ifndef __CS_COMBUSTION_BSH_HEADERS_H__
#define __CS_COMBUSTION_BSH_HEADERS_H__

/*============================================================================
 * Burke Schumann combustion model.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

#include "cogz/cs_combustion_gas.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Global variables
 *============================================================================*/

/*! Burke Schumann combustion model thermal coefficients */

extern cs_real_t coeff_therm[7][2][5];

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the fluid properties from the Burke-Schumann combustion model.
 */
/*----------------------------------------------------------------------------*/

void
cs_compute_burke_schumann_properties(cs_real_t  z_m_0,
                                     cs_real_t  zvar_0,
                                     cs_real_t  xr_m_0,
                                     cs_real_t *phi_t);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the parameters needed for the Burke-Schumann combustion
 * model.
 */
/*----------------------------------------------------------------------------*/

void
cs_burke_schumann(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculate enthalpy using the Burke-Schumann model.
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_compute_burke_schumann_enthalpy
(
  cs_real_t t,
  cs_real_t yspec[CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS]
);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_COMBUSTION_BSH_HEADERS_H__ */
