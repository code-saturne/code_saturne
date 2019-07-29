#ifndef __CS_RAD_TRANSFER_ABSORPTION_H__
#define __CS_RAD_TRANSFER_ABSORPTION_H__

/*============================================================================
 * Radiative transfer module, absorption coefficient.
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute absorption coefficient for gas mix and particles for
 *        pulverized coal or other specific physical models.
 *
 * For the P-1 model, this function also checks whether the medium's optical
 * length is at least of the order of unity.
 *
 * \param[in]   tempk  gas phase temperature at cells (in Kelvin)
 * \param[out]  cpro_caki0 Medium (gas) Absorption coefficient
 * \param[out]  kgas   radiation coefficients of the gray gases at cells
 *                      (per gas)
 * \param[out]  agas   weights of the gray gases at cells (per gas)
 * \param[out]  agasb  weights of the gray gases at boundary faces (per gas)
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_absorption(const cs_real_t  tempk[],
                           cs_real_t        cpro_cak0[],
                           cs_real_t        kgas[],
                           cs_real_t        agas[],
                           cs_real_t        agasb[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check the absorption validity fo the P-1 approximation.
 *
 * For the P-1 model, the medium's optical length should be at least of
 * the order of unity.
 *
 * \param[in]  cpro_cak  absorption coefficient values (at cells)
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_absorption_check_p1(const cs_real_t  cpro_cak[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RAD_TRANSFER_ABSORPTION_H__ */
