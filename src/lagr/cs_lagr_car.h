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
 * \brief Compute particle characteristics (except force_p): Tp, TL and PI
 * as well as covariance and variance tensors for the stochastic model
 *
 * \param[in] iprev                time step indicator for fields
 *                                   0: use fields at current time step
 *                                   1: use fields at previous time step
 * \param[in]  phase_id            carrier phase id
 * \param[in]  ip                  particle index in set
 * \param[in]  nor                 current step id (for 2nd order scheme)
 * \param[in]  dt_part             time step associated to the particle
 * \param[out] taup                dynamic characteristic time
 * \param[out] tlag                fluid characteristic Lagrangian time scale
 * \param[out] piil                term in integration of up sde
 * \param[out] bx                  turbulence characteristics
 * \param[out] tempct              thermal characteristic time
 * \param[out] beta                for the extended scheme
 * \param[out] vagaus              gaussian random variables
 * \param[out] br_gaus             gaussian random variables
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_car(int                iprev,
            int                phase_id,
            cs_lnum_t          ip,
            int                nor,
            const cs_real_t    dt_part,
            cs_real_t         *taup,
            cs_real_3_t        tlag,
            cs_real_3_t        piil,
            cs_real_33_t       bx,
            cs_real_2_t        tempct,
            cs_real_3_t        beta,
            cs_real_3_t       *vagaus,
            cs_real_6_t        br_gaus);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute external force impacting the particle
 *
 * \param[in]  dt_part             time step associated to the particle
 * \param[in]  ip                  particle index in set
 * \param[in]  taup                dynamic characteristic time
 * \param[in]  tlag                fluid characteristic Lagrangian time scale
 * \param[in]  piil                term in integration of up sde
 * \param[in]  bx                  turbulence characteristics
 * \param[in]  tempct              thermal characteristic time
 * \param[in]  beta                for the extended scheme
 * \param[in]  tsfext              info for return coupling source terms
 * \param[in]  vagaus              gaussian random variables
 * \param[in]  force_p             user external force field (m/s^2)$
 *
 */

void
cs_lagr_get_force_p(const cs_real_t    dt_part,
                    cs_lnum_t          ip,
                    cs_real_t         *taup,
                    cs_real_3_t       *tlag,
                    cs_real_3_t       *piil,
                    cs_real_33_t      *bx,
                    cs_real_t          tsfext,
                    cs_real_3_t       *vagaus,
                    cs_real_3_t        force_p);
/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_LAGCAR_H__ */
