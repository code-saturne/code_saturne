#ifndef __CS_TURBUBELNCE_KE_H__
#define __CS_TURBUBELNCE_KE_H__

/*============================================================================
 * k-epsilon turbulence model.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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
 * \brief Solve the k-epsilon equations.
 *
 * Solve the \f$ k - \varepsilon \f$  for incompressible flows
 * or slightly compressible flows for one time step.
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
 * \param[out]    prdv2f        v2f production term
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_ke(cs_lnum_t        ncesmp,
                 cs_lnum_t        icetsm[],
                 int              itypsm[],
                 const cs_real_t *dt,
                 cs_real_t        smacel[],
                 cs_real_t       *prdv2f);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TURBUBELNCE_KE_H__ */
