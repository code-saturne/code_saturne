#ifndef __CS_RAD_TRANSFER_MODAK_H__
#define __CS_RAD_TRANSFER_MODAK_H__

/*============================================================================
 * Radiation solver operations.
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
 * \brief Compute radiative properties of a gas based on temperature,
 *        composition of products in co2, h2o and soot
 *        using linear regressions established by Modak.
 *
 * \param[out] ck     medium's absorbtion coefficient (0 if transparent)
 * \param[in]  pco2   CO2 partial pressure
 * \param[in]  ph2o   H2O partial pressure
 * \param[in]  fv     soot volume fraction
 * \param[in]  temp   temperature
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_modak(cs_real_t        ck[],
                      const cs_real_t  pco2[],
                      const cs_real_t  ph2o[],
                      const cs_real_t  fv[],
                      const cs_real_t  temp[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RAD_TRANSFER_MODAK_H__ */
