#ifndef __CS_LAGR_LAGESP_H__
#define __CS_LAGR_LAGESP_H__

/*============================================================================
 * Functions and types for LAGESP
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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Integration of particle equations of motion:
 *
 * - Standard Model : First or second order
 * - Deposition submodel (Guingo & Minier, 2008) if needed
 *
 * \param[in]  dt_p      lagrangian time step
 * \param[in]  taup      dynamic characteristic time
 * \param[in]  tlag      fluid characteristic time
 * \param[in]  piil      terme in P-U SDE integration
 * \param[in]  bx        turbulence characteristics
 * \param[out] tsfext    info for return coupling source terms
 * \param[in]  gradpr    pressure gradient
 * \param[in]  gradvf    fluid velocity gradient
 * \param[out] terbru    Diffusion coefficient accounting for Brownian
 *                       (molecular) effect
 * \param[in]  vislen    nu/u* = y/y+
 * \param[in]  beta      proportional to the gradient of T_lag
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_sde(cs_real_t           dt_p,
            const cs_real_t     taup[],
            const cs_real_3_t   tlag[],
            const cs_real_3_t   piil[],
            const cs_real_33_t  bx[],
            cs_real_t           tsfext[],
            const cs_real_3_t   gradpr[],
            const cs_real_33_t  gradvf[],
            cs_real_t           terbru[],
            const cs_real_t     vislen[],
            const cs_real_3_t   beta[],
            cs_lnum_t          *nresnew);

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
 * \param[in]  tcarac  variable characteristic time
 * \param[in]  pip     right-hand side associated with SDE
 *----------------------------------------------------------------------------*/

void
cs_lagr_sde_attr(cs_lagr_attribute_t   attr,
                 cs_real_t            *tcarac,
                 cs_real_t            *pip);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_LAGESP_H__ */
