#ifndef CS_COMBUSTION_D3P_H
#define CS_COMBUSTION_D3P_H

/*============================================================================
 * Infinitely fast 3-point combustion model
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
 * Standard library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"
#include "base/cs_field.h"
#include "pprt/cs_physical_model.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
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
/*!
 * \brief Initialize specific fields for infinitely fast 3-point chemistry
 *        combustion model, stage 0.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_d3p_fields_init0(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize specific fields for infinitely fast 3-point chemistry
 *        combustion model, stage 1.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_d3p_fields_init1(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Integration of thermodynamic variablesfunction of mixture fraction.
 *
 * Remark: temperature integration could be weighted by Cp.
 *
 * \param[in]     indpdf        indicator for pdf integration or mean value
 * \param[in]     dirmin        Dirac's peak value at \f$ f_{min} \f$
 * \param[in]     dirmax        Dirac's peak value at \f$ f_{max} \f$
 * \param[in]     fdeb          abscissa of rectangle low boundary
 * \param[in]     ffin          abscissa of rectangle high boundary
 * \param[in]     hrec          rectangle height
 * \param[in]     w1            work array
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_d3p_integration(int        indpdf[],
                              cs_real_t  dirmin[],
                              cs_real_t  dirmax[],
                              cs_real_t  fdeb[],
                              cs_real_t  ffin[],
                              cs_real_t  hrec[],
                              cs_real_t  w1[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* CS_COMBUSTION_D3P_H */
