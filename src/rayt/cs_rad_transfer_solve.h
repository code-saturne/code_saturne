#ifndef __CS_RAD_TRANSFER_DOM_H__
#define __CS_RAD_TRANSFER_DOM_H__

/*============================================================================
 * Radiation solver boundary conditions treatment.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve the radiative transfer equation.
 *
 *  Two types of method are available:
 *  - Discretes Ordinates Methods (DOM)
 *  - P-1 approximation (only recommended for pulverized coal)
 *
 *  \param[in, out]  bc_type       boundary face types
 *  \param[in]       nclacp        number of pulverized coal classes
 *  \param[in]       nclafu        number of fuel classes
 *  \param[in]       dt            time step (per cell)
 *  \param[in]       cp2fol        fuel oil liquid CP
 *  \param[in]       cp2ch         pulverized coal CP's
 *  \param[in]       ichcor        pulverized coal indirection
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_solve(int               bc_type[],
                      int               nclacp,
                      int               nclafu,
                      const cs_real_t   dt[],
                      cs_real_t         cp2fol,
                      const cs_real_t   cp2ch[],
                      const int         ichcor[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RAD_TRANSFER_DOM_H__ */
