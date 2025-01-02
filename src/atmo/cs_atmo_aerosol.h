#ifndef __CS_ATMO_AEROSOL_H__
#define __CS_ATMO_AEROSOL_H__

/*============================================================================
 * Main for atmospheric aerosols related functions
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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

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
 * \brief This function initializes the external aerosol code
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_initialize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function finalizes the external aerosol code.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function fills the given array with aerosol concentrations
 *        and numbers from the external aerosol code..
 *
 * \param[out]  array  aerosol concentrations and numbers
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_get_aero(cs_real_t  *array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function fills the given array with gas concentrations from
 *        the external aerosol code.
 *
 * \param[out]  array  gas concentrations
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_get_gas(cs_real_t  *array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function computes a time step of gaseous chemistry and aerosols
 *        dynamic using the external aerosol code.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_time_advance(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute aerosol cloud droplets nucleation when using the atmospheric
 * humid model using a microphysical model.
 *
 * It is taken into account as an additional step split from advection-diffusion
 * equation, hence the droplet number is first clipped if necessary.
 *
 * \param[out]  nc      droplet number (scalar) in 1/cm**3
 * \param[in]   rom     density of air in kg/m**3
 * \param[in]   qldia   mass fraction of liquid water in kg/kg
 * \param[in]   pphy    true pressure in pascals
 * \param[in]   refrad  radiative cooling
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_nuclea(cs_real_t         *nc,
                       const cs_real_t   *rom,
                       const cs_real_t   *qldia,
                       const cs_real_t   *pphy,
                       const cs_real_t   *refrad);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ATMO_AEROSOL_H__ */
