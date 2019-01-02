#ifndef __CS_LAGR_ADHESION_H__
#define __CS_LAGR_ADHESION_H__

/*============================================================================
 * Functions and types for the Lagrangian module
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
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Calculation of the adhesion force and adhesion energy
 *
 * \param[in]  ip               particle number
 * \param[in]  tempf            thermal scalar value at current time step
 * \param[out] adhesion_energ   particle adhesion energy
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_adh(cs_lnum_t   ip,
            cs_real_t   tempf,
            cs_real_t  *adhesion_energ);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Calculation of the adhesion force and adhesion energy
 *
 * \param[in]  tempf            thermal scalar value at current time step
 * \param[out] adhesion_energ   particle adhesion energy
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_adh_pp(cs_real_t   rpart,
               cs_real_t   tempf,
               cs_real_t  *adhesion_energ,
               cs_real_t  *adhesion_force);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_ADHESION_H__ */
