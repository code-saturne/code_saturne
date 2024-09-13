#ifndef CS_COAL_PHYSICAL_PROPERTIES_H
#define CS_COAL_PHYSICAL_PROPERTIES_H

/*============================================================================
 * Coal combustion model: physical properties
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Compute physical properties in gaseous phase.
 *
 * \param[in]     f1m     mean of tracer 1 mvl [chx1m+co]
 * \param[in]     f2m     mean of tracer 2 mvl [chx2m+co]
 * \param[in]     f3m     mean of tracer 3 (oxydant 1)
 * \param[in]     f4m     mean of tracer 4 (oxydant 2)
 * \param[in]     f5m     mean of tracer 5 (oxydant 3)
 * \param[in]     f6m     mean of tracer 6 (humidity)
 * \param[in]     f7m     mean of tracer 7 (C + O2)
 * \param[in]     f8m     mean of tracer 8 (C + CO2)
 * \param[in]     f9m     mean of tracer 9 (C + H2O)
 * \param[in]     fvp2m   f1f2 variance
 * \param[in]     enth    enthalpy in \f$ j . kg^{-1} \f$  either of
 *                        the gas or of the mixture
 * \param[in]     enthox  oxydant enthalpy
 * \param[out]    rom1    gas density
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_physprop1(const  cs_real_t  f1m[],
                  const  cs_real_t  f2m[],
                  const  cs_real_t  f3m[],
                  const  cs_real_t  f4m[],
                  const  cs_real_t  f5m[],
                  const  cs_real_t  f6m[],
                  const  cs_real_t  f7m[],
                  const  cs_real_t  f8m[],
                  const  cs_real_t  f9m[],
                  const  cs_real_t  fvp2m[],
                  const  cs_real_t  enth[],
                  const  cs_real_t  enthox[],
                  cs_real_t         rom1[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* CS_COAL_PHYSICAL_PROPERTIES */
