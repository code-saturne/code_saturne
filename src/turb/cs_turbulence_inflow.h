
#ifndef __CS_TURBULENCE_INFLOW_H__
#define __CS_TURBULENCE_INFLOW_H__

/*============================================================================
 * Turbulent inflow generation
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
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
 * \brief Define mass injection for turbulent quantities based
 *        on k and epsilon values.
 *
 * \param[in]  zone_name  name of zone to which injection should be added
 * \param[in]  k          turbulent kinetic energy
 * \param[in]  eps        turbulent dissipation
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_inflow_volume_mass_injection_k_eps(const char  *zone_name,
                                                 double       k,
                                                 double       eps);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define mass injection for turbulent quantities based
 *        on a hydraulic diameter and reference velocity.
 *
 * \param[in]  zone_name  name of zone to which injection should be added
 * \param[in]  uref2      square of the reference flow velocity
 * \param[in]  dh         hydraulic diameter \f$ D_H \f$
 * \param[in]  rho        mass density \f$ \rho \f$
 * \param[in]  mu         dynamic viscosity \f$ \nu \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_inflow_volume_mass_injection_ke_hyd_diam(const char  *zone_name,
                                                       double       uref2,
                                                       double       dh,
                                                       double       rho,
                                                       double       mu);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TURBULENCE_INFLOW_H__ */
