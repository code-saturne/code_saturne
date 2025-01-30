#ifndef __CS_LAGR_LAGCAR_H__
#define __CS_LAGR_LAGCAR_H__

/*============================================================================
 * Compute particle characteristics: Tp, TL et PI
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

#include "base/cs_defs.h"

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
 * \brief Compute particle characteristics: Tp, TL and PI as well as
 * covariance and variance tensors for the Stochastic model
 *
 * \param[in] iprev             time step indicator for fields
 *                                0: use fields at current time step
 *                                1: use fields at previous time step
 * \param[in]  phase_id         carrier phase id
 * \param[in]  dt               time step (per cell)
 * \param[out] taup             dynamic characteristic time
 * \param[out] tlag             fluid characteristic time
 * \param[out] piil             term in integration of up sde
 * \param[out] bx               turbulence characteristics
 * \param[out] tempct           thermal characteristic time
 * \param[out] beta             for the extended scheme
 * \param[in]  gradvf           fluid velocity gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_car(int              iprev,
            int              phase_id,
            const cs_real_t  dt[],
            cs_real_t        *taup,
            cs_real_3_t      *tlag,
            cs_real_3_t      *piil,
            cs_real_33_t     *bx,
            cs_real_t        tempct[],
            cs_real_3_t      *beta,
            cs_real_33_t     gradvf[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_LAGCAR_H__ */
