#ifndef __CS_RAD_TRANSFER_FSCK_H__
#define __CS_RAD_TRANSFER_FSCK_H__

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
 *  \brief Determine the radiation coefficients of the FSCK model
 *         as well as the corresponding weights.
 *
 * \param[in]     pco2        CO2 volume fraction
 * \param[in]     ph2o        H2O volume fraction
 * \param[in]     teloc       gas temperature
 * \param[out]    kloc        radiation coefficient of the i different gases
 * \param[out]    aloc        weights of the i different gases in cells
 * \param[out]    alocb       weights of the i different gases at boundaries
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_fsck(const cs_real_t  *restrict pco2,
                     const cs_real_t  *restrict ph2o,
                     const cs_real_t  *restrict teloc,
                     cs_real_t        *restrict kloc,
                     cs_real_t        *restrict aloc,
                     cs_real_t        *restrict alocb);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RAD_TRANSFER_FSCK_H__ */
