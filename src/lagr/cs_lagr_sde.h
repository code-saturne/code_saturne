#ifndef __CS_LAGR_LAGESP_H__
#define __CS_LAGR_LAGESP_H__

/*============================================================================
 * Functions and types for LAGESP
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

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \brief Integration of SDEs by 1st order time scheme for one particle
 *
 * \param[in]  p_id      particle index in set
 * \param[in]  dt_part   remaining time step associated to the particle
 * \param[in]  nor       current step id (for 2nd order scheme)
 * \param[in]  taup      dynamic characteristic time
 * \param[in]  tlag      lagrangian fluid characteristic time
 * \param[in]  piil      term in integration of up sdes
 * \param[in]  bx        turbulence characteristics
 * \param[in]  vagaus    gaussian random variables
 * \param[in]  brgaus    gaussian random variables
 * \param[in]  force_p   forces per mass unit on particles (m/s^2)
 * \param[in]  beta      proportional to the gradient of T_lag
 * \param[out] terbru    Diffusion coefficient accounting for Brownian
 *                       (molecular) effect
 */
/*----------------------------------------------------------------------------*/
void
cs_sde_vels_pos_1_st_order_time_integ(cs_lnum_t                       p_id,
                                      cs_real_t                       dt_part,
                                      int                             nor,
                                      const cs_real_t                *taup,
                                      const cs_real_3_t              *tlag,
                                      const cs_real_3_t              *piil,
                                      const cs_real_33_t             *bx,
                                      const cs_real_3_t              *vagaus,
                                      const cs_real_6_t               brgaus,
                                      const cs_real_3_t               force_p,
                                      const cs_real_3_t               beta,
                                      cs_real_t                      *terbru);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Integration of particle equations of motion:
 *
 * - Standard Model : First or second order
 * - Deposition submodel (Guingo & Minier, 2007) if needed
 *
 * \param[in]  p_id      particle index in set
 * \param[in]  dt_part   remaining time step associated to the particle
 * \param[in]  nor       current step id (for 2nd order scheme)
 * \param[in]  taup      dynamic characteristic time
 * \param[in]  tlag      fluid characteristic time
 * \param[in]  piil      term in integration of U-P SDEs
 * \param[in]  bx        turbulence characteristics
 * \param[out] tsfext    info for return coupling source terms
 * \param[out] force_p   forces per mass unit on particles (m/s^2)
 * \param[out] terbru    Diffusion coefficient accounting for Brownian
 *                       (molecular) effect
 * \param[in]  vislen    nu/u* = y/y+
 * \param[in]  beta      proportional to the gradient of T_lag
 * \param[out] vagaus    gaussian random variables
 * \param[out] brgaus    gaussian random variables
 * \param[in]  nresnew
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_sde(cs_lnum_t                        p_id,
            cs_real_t                        dt_part,
            int                              nor,
            const cs_real_t                 *taup,
            const cs_real_3_t               *tlag,
            const cs_real_3_t               *piil,
            const cs_real_33_t              *bx,
            cs_real_t                       *tsfext,
            const cs_real_3_t                force_p,
            cs_real_t                       *terbru,
            const cs_real_t                  vislen[],
            const cs_real_3_t                beta,
            cs_real_3_t                     *vagaus,
            cs_real_6_t                      brgaus,
            cs_lnum_t                       *nresnew);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Integration of a stochastic differential equation (SDE) for
 *        a user particle variable (attribute).
 *
 * \f[
 *  \frac{dV}{dt} = \frac{V - PIP}{TCARAC}
 * \f]
 *
 * When there is interaction with a boundary face, the integration
 * degenerates to order 1 (even if the 2nd order scheme is active).
 *
 * \param[in]  attr    attribute/variable
 * \param[in]  p_id    particle id
 * \param[in]  nor     current step id (for 2nd order scheme)
 * \param[in]  dt_part remaining time step associated to the particle
 * \param[in]  tcarac  variable characteristic time
 * \param[in]  pip     right-hand side associated with SDE
 *----------------------------------------------------------------------------*/

void
cs_lagr_sde_attr(cs_lagr_attribute_t   attr,
                 const cs_lnum_t       p_id,
                 int                   nor,
                 const cs_real_t       dt_part,
                 cs_real_t             tcarac,
                 cs_real_t             pip);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_LAGESP_H__ */
