#ifndef __CS_DRIFT_CONVECTIVE_FLUX_H__
#define __CS_DRIFT_CONVECTIVE_FLUX_H__

/*============================================================================
 * Compute the modified convective flux for scalars with a drift.
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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Compute the modified convective flux for scalars with a drift.
 *
 * \param[in]     f_sc          drift scalar field
 * \param[in,out] i_mass_flux   scalar mass flux at interior face centers
 * \param[in,out] b_mass_flux   scalar mass flux at boundary face centers
 * \param[in,out] divflu        divergence of drift flux
 */
/*----------------------------------------------------------------------------*/

void
cs_drift_convective_flux(cs_field_t  *f_sc,
                         cs_real_t    i_mass_flux[],
                         cs_real_t    b_mass_flux[],
                         cs_real_t    divflu[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_DRIFT_CONVECTIVE_FLUX_H__ */
