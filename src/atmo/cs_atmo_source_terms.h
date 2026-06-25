#ifndef CS_ATMO_SOURCE_TERMS_H
#define CS_ATMO_SOURCE_TERMS_H

/*============================================================================
 * Main for atmospheric source terms
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_base.h"
#include "base/cs_field.h"

/*----------------------------------------------------------------------------*/

/*============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Phase change source terms - Exchange terms between the injected
 *        liquid and the water vapor phase in the bulk, humid air
 *
 * \param[in]     f_id          field id
 * \param[in,out] exp_st        Explicit source term
 * \param[in,out] imp_st        Implicit source term
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_source_term(int              f_id,
                    cs_real_t        exp_st[],
                    cs_real_t        imp_st[]);

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Additional right-hand side source terms for scalar equations
 *        taking into account dry and humid atmospheric variables.
 *        If 1D atmospheric radiative module is used additional source terms for
 *        the thermal scalar equation to take into account the radiative forcing.
 *
 * \param[in]     f_id          field id
 * \param[in,out] exp_st        Explicit source term
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_scalar_source_term(int              f_id,
                           cs_real_t        exp_st[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Additional right-hand side source terms
 *         for momentum equation in case of free inlet
 *
 * \param[in,out] exp_st        Explicit source term
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_source_term_for_inlet(cs_real_3_t        exp_st[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Compute aerosol cloud droplets nucleation when using the atmospheric
 * humid model using a microphysical model.
 *
 * It is taken into account as an additional step split from advection-diffusion
 * equation, hence the droplet number is first clipped if necessary.
 *
 * \param[out]  f_ntdrp droplet number field in 1/cm**3
 * \param[in]   rom     density of air in kg/m**3
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_nuclea(cs_field_t        *f_ntdrp,
                       const cs_real_t   *rom);

/*----------------------------------------------------------------------------*/

#endif /* CS_ATMO_SOURCE_TERMS_H */
