#ifndef __CS_RAD_TRANSFER_WALL_FLUX_H__
#define __CS_RAD_TRANSFER_WALL_FLUX_H__

/*============================================================================
 * Radiation solver operations.
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
 * \param[in]  nvar     number of variable BC's
 * \param[in]  ivart    variable id of thermal variable
 * \param[in]  isothp   list of isothermal boundaries
 * \param[in]  rcodcl   boundary condition values
 *                        rcodcl[0] = Dirichlet value
 *                        rcodcl[1] = exchange coefficient value.
 *                                    (infinite if no exchange)
 *                        rcodcl[2] = flux density value (negative if gain)
 *                                    in w/m2 or rugosity height (m)
 *                                    if icodcl=6,
 *                                    - for velocities, (vistl+visct)*gradu
 *                                    - for pressure, dt*gradp,
 *                                    - for scalars,
 *                                      cp*(viscls+visct/turb_schmidt)*gradt
 * \param[out] tparop   wall temperature in Kelvin
 * \param[in]  qincip   radiative flux density at boundaries
 * \param[in]  textp    exterior boundary temperature in degrees C
 * \param[in]  tintp    interior boundary temperature in degrees C
 * \param[in]  xlamp    thermal conductivity coefficient of wall faces (w/m/k)
 * \param[in]  epap     wall thickness (m)
 * \param[in]  epsp     wall emissivity
 * \param[in]  hfconp   boundary fluid exchange coefficient
 * \param[in]  flconp   boundary convective flux density
 * \param[in]  tempkp   temperature in Kelvin
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_wall_flux(int         nvar,
                          int         ivart,
                          int         isothp[],
                          cs_real_t  *tmin,
                          cs_real_t  *tmax,
                          cs_real_t  *tx,
                          cs_real_t  *rcodcl,
                          cs_real_t   tparop[],
                          cs_real_t   qincip[],
                          cs_real_t   textp[],
                          cs_real_t   tintp[],
                          cs_real_t   xlamp[],
                          cs_real_t   epap[],
                          cs_real_t   epsp[],
                          cs_real_t   hfconp[],
                          cs_real_t   flconp[],
                          cs_real_t   tempkp[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RAD_TRANSFER_WALL_FLUX_H__ */
