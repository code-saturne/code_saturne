#ifndef CS_TURBULENCE_RIT_H
#define CS_TURBULENCE_RIT_H

/*============================================================================
 * Turbulence transport equation.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add the divergence of turbulent flux to a scalar transport equation.
 *
 * \param[in]   field_id  transported field id
 * \param[in]   xcpp      Cp
 * \param[out]  vistet    diffusivity tensor
 * \param[out]  rhs       right hand side to update
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_turbulence_rit_div(const int        field_id,
                      const cs_real_t  xcpp[],
                      cs_real_t        vistet[][6],
                      cs_real_t        rhs[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief CS_TURB_RIJ_SOURCE_TS_EXPONENTIAL: exact frozen-tau integration of the
 *        coupled Rotta-Monin/Boussinesq-buoyancy source subsystem over
 *        one cell, over a time step dt.
 *
 * Restricted to pure Rotta closure (crij2 == 0). All resonance cases
 * (C_R = C_theta, C_R + C_theta = 2, C_theta = 1, C_R = 1, and
 * combinations) are handled via their exact limiting formula.
 *
 * \param[in]   cr        Rotta constant C_R (crij1)
 * \param[in]   ctheta    thermal relaxation ratio C_theta (= 1/rvarfl)
 * \param[in]   ceps2     Ceps2 (epsilon destruction constant)
 * \param[in]   beta      thermal expansion coefficient beta_theta
 * \param[in]   grav      gravity vector
 * \param[in]   dt        time step
 * \param[in]   r0        Rij at start of step (6 components)
 * \param[in]   qtheta0   turbulent heat flux at start of step
 * \param[in]   theta2_0  temperature variance at start of step
 * \param[in]   eps0      epsilon at start of step
 * \param[out]  r1        Rij at end of step
 * \param[out]  qtheta1   turbulent heat flux at end of step
 * \param[out]  theta2_1  temperature variance at end of step
 * \param[out]  eps1      epsilon at end of step
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rit_source_step_frozen_tau(cs_real_t          cr,
                                         cs_real_t          ctheta,
                                         cs_real_t          ceps2,
                                         cs_real_t          beta,
                                         const cs_real_t    grav[3],
                                         cs_real_t          dt,
                                         const cs_real_6_t  r0,
                                         const cs_real_3_t  qtheta0,
                                         cs_real_t          theta2_0,
                                         cs_real_t          eps0,
                                         cs_real_6_t        r1,
                                         cs_real_3_t        qtheta1,
                                         cs_real_t         *theta2_1,
                                         cs_real_t         *eps1);

/*----------------------------------------------------------------------------*/
/*!
 * \brief CS_TURB_RIJ_SOURCE_TS_VAR_TAU: variable-tau trajectory integration
 *        of the coupled Rotta-Monin/Boussinesq-buoyancy source subsystem
 *        over one cell, over a time step dt. Restricted to pure Rotta
 *        closure (crij2 == 0) and to a linear epsilon-buoyancy closure
 *        (cs_turb_ce4 == 0).
 *
 * \param[in]   cr        Rotta constant C_R (crij1)
 * \param[in]   ctheta    thermal relaxation ratio C_theta (= 1/rvarfl)
 * \param[in]   ceps2     Ceps2 (epsilon destruction constant)
 * \param[in]   ceps3     Ceps3 (linear epsilon-buoyancy constant,
 *                        cs_turb_ce3)
 * \param[in]   beta      thermal expansion coefficient beta_theta
 * \param[in]   grav      gravity vector
 * \param[in]   dt        time step
 * \param[in]   r0        Rij at start of step (6 components)
 * \param[in]   qtheta0   turbulent heat flux at start of step
 * \param[in]   theta2_0  temperature variance at start of step
 * \param[in]   eps0      epsilon at start of step
 * \param[out]  r1        Rij at end of step
 * \param[out]  qtheta1   turbulent heat flux at end of step
 * \param[out]  theta2_1  temperature variance at end of step
 * \param[out]  eps1      epsilon at end of step
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rit_source_step_variable_tau(cs_real_t          cr,
                                           cs_real_t          ctheta,
                                           cs_real_t          ceps2,
                                           cs_real_t          ceps3,
                                           cs_real_t          beta,
                                           const cs_real_t    grav[3],
                                           cs_real_t          dt,
                                           const cs_real_6_t  r0,
                                           const cs_real_3_t  qtheta0,
                                           cs_real_t          theta2_0,
                                           cs_real_t          eps0,
                                           cs_real_6_t        r1,
                                           cs_real_3_t        qtheta1,
                                           cs_real_t         *theta2_1,
                                           cs_real_t         *eps1);

/*----------------------------------------------------------------------------*/

#endif /* CS_TURBULENCE_RIT_H */
