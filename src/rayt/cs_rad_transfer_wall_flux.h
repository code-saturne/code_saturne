#ifndef __CS_RAD_TRANSFER_WALL_FLUX_H__
#define __CS_RAD_TRANSFER_WALL_FLUX_H__

/*============================================================================
 * Radiation solver operations.
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Wall temperature computation with flux balance.
 *
 * \param[in]  isothp   list of isothermal boundaries
 * \param[in]  tmin     minimum allowed temperature (clip to this value)
 * \param[in]  tmax     maximum allowed temperature (clip to this value)
 * \param[in]  tx       temperature relaxtion parameter
 * \param[in]  qincip   radiative flux density at boundaries
 * \param[in]  textp    exterior boundary temperature in degrees C
 * \param[in]  xlamp    thermal conductivity coefficient of wall faces (w/m/k)
 * \param[in]  epap     wall thickness (m)
 * \param[in]  epsp     wall emissivity
 * \param[in]  hfconp   boundary fluid exchange coefficient
 * \param[in]  flconp   boundary convective flux density
 * \param[in]  tempkp   temperature in Kelvin
 * \param[out] twall    wall temperature in Kelvin
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_compute_wall_t(int         isothp[],
                               cs_real_t   tmin,
                               cs_real_t   tmax,
                               cs_real_t   tx,
                               cs_real_t   qincip[],
                               cs_real_t   textp[],
                               cs_real_t   xlamp[],
                               cs_real_t   epap[],
                               cs_real_t   epsp[],
                               cs_real_t   hfconp[],
                               cs_real_t   flconp[],
                               cs_real_t   tempkp[],
                               cs_real_t   twall[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RAD_TRANSFER_WALL_FLUX_H__ */
