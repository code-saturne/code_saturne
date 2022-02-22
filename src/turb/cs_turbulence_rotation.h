#ifndef __CS_TURBULENCE_ROTATION_H__
#define __CS_TURBULENCE_ROTATION_H__

/*============================================================================
 * Computing rotation/curvature correction.
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute rotation/curvature correction for eddy-viscosity models.
 *
 * This function is called for the linear eddy viscosity RANS models,
 * when irccor = 1 is verified.
 *
 * \param[in]   dt      time step (per cell)
 * \param[out]  rotfct  rotation function of Spalart-Shur correction
 *                      at cell center
 * \param[out]  ce2rc   modified ce2 coeficient of Cazalbou correction
 *                      at cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rotation_correction(const cs_real_t   dt[],
                                  cs_real_t         rotfct[],
                                  cs_real_t         ce2rc[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TURBULENCE_ROTATION_H__ */
