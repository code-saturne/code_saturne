#ifndef __CS_TURBULENCE_V2F_H__
#define __CS_TURBULENCE_V2F_H__

/*============================================================================
 * V2F turbulence model.
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
 * \brief  Solve the V2F phi-model equations.
 *
 * Solve the \f$\phi\f$ and diffusion for \f$ \overline{f} \f$
 * as part of the V2F phi-model
 *
 * \param[in]     ncesmp        number of cells with mass source term
 * \param[in]     icetsm        index of cells with mass source term
 * \param[in]     itypsm        mass source type for the variables
 *                              size: [nvar][ncesmp]
 * \param[in]     dt            time step (per cell)
 * \param[in]     smacel        values of the variables associated to the
 *                              mass source (for the pressure variable,
 *                              smacel is the mass flux)
 *                              size: [nvar][ncesmp]
 * \param[in]     prdv2f        v2f production term
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_v2f(cs_lnum_t         ncesmp,
                  cs_lnum_t         icetsm[],
                  int               itypsm[],
                  const cs_real_t  *dt,
                  cs_real_t         smacel[],
                  const cs_real_t   prdv2f[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of turbulent viscosity for the V2F-phi model.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_v2f_phi_mu_t(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of turbulent viscosity for the V2F-BL model.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_v2f_bl_v2k_mu_t(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TURBULENCE_V2F_H__ */
