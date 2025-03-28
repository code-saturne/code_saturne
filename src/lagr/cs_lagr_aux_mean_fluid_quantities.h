#ifndef __CS_LAGR_AUX_MEAN_FLUID_QUANTITIES_H__
#define __CS_LAGR_AUX_MEAN_FLUID_QUANTITIES_H__

/*============================================================================
 * Functions and types for lagrangian field gradients
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
 * \brief Compute anisotropic fluid quantities for complete model (modpl == 1).
 *
 *  - Anisotropic Lagragian time
 *  - Anisotropic Diffusion matrix
 *  - Anisotropic gradient of Lagrangian time in the relativ basis
 *
 * \param[in]   iprev                  time step indicator for fields
 *                                       0: use fields at current time step
 *                                       1: use fields at previous time step
 * \param[in]   phase_id               carrier phase id
 * \param[out]  anisotropic_lagr_time  anisotropic Lagragian time scale (modcpl)
 * \param[out]  anisotropic_bx         anisotropic diffusion term (if modcpl)
 * \param[out]  grad_lagr_time_r_et    anisotropic Lagrangian time gradient in
 *                                     relative basis
 * \param[in]   grad_lagr_time         Lagrangian time gradient
 *
 */
/*----------------------------------------------------------------------------*/

void
compute_anisotropic_prop(int            iprev,
                         int            phase_id,
                         cs_real_3_t   *anisotropic_lagr_time,
                         cs_real_3_t   *anisotropic_bx,
                         cs_real_3_t   *grad_lagr_time_r_et,
                         cs_real_3_t   *grad_lagr_time);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute gradient of particle covariance.
 *
 *  - particle velocity and particle velocity seen covariance
 *  - particle velocity seen variance
 *
 * \param[in]  phase_id        carrier phase id
 * \param[out] grad_cov_skp    gradient of particle velocity and
 *                             particle velocity seen covariance
 *
 * \param[out] grad_cov_sk     gradient of particle velocity seen covariance
 */
/*----------------------------------------------------------------------------*/

void
compute_particle_covariance_gradient(int phase_id,
                                     cs_real_3_t *grad_cov_skp[9],
                                     cs_real_3_t *grad_cov_sk[6]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute auxilary mean fluid quantities.
 *
 *  - Lagrangian time (if modcpl == 1 also anisotropic values)
 *  - gradient of total pressure
 *  - velocity gradient
 *  - temperature gradient
 *  - Lagrangian time gradient (also gradient of anisotropic values if needed)
 *  - diffusion matix
 *
 * \param[in]   iprev                  time step indicator for fields
 *                                      0: use fields at current time step
 *                                      1: use fields at previous time step
 * \param[in]   phase_id               carrier phase id
 * \param[out]  lagr_time              Lagrangian time scale
 * \param[out]  grad_pr                pressure gradient
 * \param[out]  grad_vel               velocity gradient
 * \param[out]  grad_tempf             fluid temperature gradient
 * \param[out]  grad_lagr_time         Lagrangian time gradient
 * \param[out]  anisotropic_lagr_time  anisotropic Lagragian time scale (modcpl)
 * \param[out]  anisotropic_bx         anisotropic diffusion term (if modcpl)
 * \param[out]  grad_lagr_time_r_et    anisotropic Lagrangian time gradient in
 *                                     relative basis
 * \param[out]  grad_lagr_time         Lagrangian time gradient
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_aux_mean_fluid_quantities(int            iprev,
                                  int            phase_id,
                                  cs_field_t    *lagr_time,
                                  cs_real_3_t   *grad_pr,
                                  cs_real_33_t  *grad_vel,
                                  cs_real_3_t   *grad_tempf,
                                  cs_real_3_t   *anisotropic_lagr_time,
                                  cs_real_3_t   *anisotropic_bx,
                                  cs_real_3_t   *grad_lagr_time_r_et,
                                  cs_real_3_t   *grad_lagr_time);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_AUX_MEAN_FLUID_QUANTITIES_H__ */
