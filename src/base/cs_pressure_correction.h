#ifndef __CS_PRESSURE_CORRECTION_H__
#define __CS_PRESSURE_CORRECTION_H__

/*============================================================================
 * Pressure correction.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform the pressure correction step of the Navier-Stokes equations
 *        for incompressible or slightly compressible flows.
 *
 * This function solves the following Poisson equation on the pressure:
 * \f[
 *     D \left( \Delta t, \delta p \right) =
 * \divs \left( \rho \vect{\widetilde{u}}\right)
 *     - \Gamma^n
 *     + \dfrac{\rho^n - \rho^{n-1}}{\Delta t}
 * \f]
 * The mass flux is then updated as follows:
 * \f[
 *  \dot{m}^{n+1}_\ij = \dot{m}^{n}_\ij
 *                    - \Delta t \grad_\fij \delta p \cdot \vect{S}_\ij
 * \f]
 *
 * \Remark:
 * - an iterative process is used to solve the Poisson equation.
 * - if the arak coefficient is set to 1, the the Rhie & Chow filter is
 *   activated.
 *
 * Please refer to the
 * <a href="../../theory.pdf#resopv"><b>resopv</b></a>
 * section of the theory guide for more information.
 *
 * \param[in]       iterns    Navier-Stokes iteration number
 * \param[in]       nfbpcd    number of faces with condensation source term
 * \param[in]       ncmast    number of cells with condensation source terms
 * \param[in]       ifbpcd    index of faces with condensation source term
 * \param[in]       ltmast    index of cells with condensation source terms
 * \param[in]       isostd    indicator of standard outlet and index
 *                            of the reference outlet face
 * \param[in]       vel       velocity
 * \param[in, out]  da_uu     velocity matrix
 * \param[in]       coefav    boundary condition array for the variable
 *                            (explicit part)
 * \param[in]       coefbv    boundary condition array for the variable
 *                            (implicit part)
 * \param[in]       coefa_dp  boundary conditions for the pressure increment
 * \param[in]       coefb_dp  boundary conditions for the pressure increment
 * \param[in]       spcond    variable value associated to the condensation
 *                            source term (for ivar=ipr, spcond is the
 *                            flow rate
 *                            \f$ \Gamma_{s,cond}^n \f$)
 * \param[in]       svcond    variable value associated to the condensation
 *                            source term (for ivar=ipr, svcond is the flow rate
 *                            \f$ \Gamma_{v, cond}^n \f$)
 * \param[in]       frcxt     external forces making hydrostatic pressure
 * \param[in]       dfrcxt    variation of the external forces
 *                            composing the hydrostatic pressure
 * \param[in]       i_visc    visc*surface/dist aux faces internes
 * \param[in]       b_visc    visc*surface/dist aux faces de bord
 */
/*----------------------------------------------------------------------------*/

void
cs_pressure_correction(int        iterns,
                       cs_lnum_t  nfbpcd,
                       cs_lnum_t  ncmast,
                       cs_lnum_t  ifbpcd[nfbpcd],
                       cs_lnum_t  ltmast[],
                       int        isostd[],
                       cs_real_t  vel[restrict][3],
                       cs_real_t  da_uu[restrict][6],
                       cs_real_t  coefav[restrict][3],
                       cs_real_t  coefbv[restrict][3][3],
                       cs_real_t  coefa_dp[restrict],
                       cs_real_t  coefb_dp[restrict],
                       cs_real_t  spcond[restrict],
                       cs_real_t  svcond[restrict],
                       cs_real_t  frcxt[restrict][3],
                       cs_real_t  dfrcxt[restrict][3],
                       cs_real_t  i_visc[restrict],
                       cs_real_t  b_visc[restrict]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PRESSURE_CORRECTION_H__ */
