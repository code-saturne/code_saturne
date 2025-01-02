#ifndef __CS_COMPUTE_THERMO_PRESSURE_DENSITY_H__
#define __CS_COMPUTE_THERMO_PRESSURE_DENSITY_H__

/*============================================================================
 * Thermodynamic pressure and density.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the density \f$ \rho^{n+1}\f$ with the
 *        \f$ \rho^{n-\frac{1}{2}} \f$ density with the state law
 *        and a thermodynamic pressure \f$ p_{ther}^{n+1} \f$
 *        estimated from the integral over the total
 *        fluid domain of the mass conservation equation.
 */
/*----------------------------------------------------------------------------*/

void
cs_compute_thermo_pressure_density(void);


/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_COMPUTE_THERMO_PRESSURE_DENSITY_H__ */
