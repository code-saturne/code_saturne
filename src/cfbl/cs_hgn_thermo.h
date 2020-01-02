#ifndef __CS_HGN_THERMO_H__
#define __CS_HGN_THERMO_H__

/*============================================================================
 * Thermodynamic for the compressible homogeneous two-phase flow model
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of the temperature at saturation with respect to the
 *        pressure.
 *
 * Compute the temperature at saturation \f$T_{sat}\f$ with respect to the
 * pressure \f$P\f$. It corresponds to the temperature for which
 * the pressure, temperature and chemical potential are equal in both phases.
 * It reduces to solve the equality of the phasic chemical potential:
 * \f[
 *   P \rightarrow T_{sat} \/ \mu_1(T_{sat},P)=\mu_2(T_{sat},P)
 * \f]
 * This equality is solved using a Newton method.
 *
 * \param[in]     pr       pressure
 *
 * \return temperature at saturation
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_hgn_thermo_saturation_temp(cs_real_t pr);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation mixture pressure and temperature from volume, mass,
 *        energy fractions, as well as specific energy and specific volume.
 *
 * Following relations are used, that rely on phasic pressures and temperatures:
 * \f[
 *   \dfrac{1}{T} = \dfrac{\dd s}{\dd e}_{|\tau,\alpha,y,z}
 * \f]
 * \f[
 *   \dfrac{P}{T} = \dfrac{\dd s}{\dd \tau}_{|e,\alpha,y,z}
 * \f]
 *
 * \param[in]     alpha   volume fraction
 * \param[in]     y       mass fraction
 * \param[in]     z       energy fraction
 * \param[in]     e       specific energy
 * \param[in]     v       specific volume
 * \param[out]    ptp     pointer to mixture temperature
 * \param[out]    ppr     pointer to mixture pressure
 */
/*----------------------------------------------------------------------------*/

void
cs_hgn_thermo_pt(cs_real_t  alpha,
                 cs_real_t  y,
                 cs_real_t  z,
                 cs_real_t  e,
                 cs_real_t  v,
                 cs_real_t *ptp,
                 cs_real_t *ppr);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of the square of the sound speed in the mixture.
 *
 * The sound speed may be computed through the Hessian matrices of the specific
 * phasic entropies in the plane (\f$\tau\f$, \f$e\f$).
 * \f$\tau\f$ stands for specific volume, and \f$e\f$ for specific energy.
 * The sound speed is here estimated using the plane (\f$\tau\f$, \f$s\f$).
 * \f$s\f$ stands for specific entropy.
 * We use the definition
 * \f\[
 *   c^2 = -\tau^2 \der{P}{\tau}_{|s,y}
 * \f\].
 * This relation is computed by a finite difference.
 *
 * \param[in]     alpha     volume fraction
 * \param[in]     y         mass fraction
 * \param[in]     z         energy fraction
 * \param[in]     P         pressure
 * \param[in]     v         specific volume
 *
 * \return square of the sound speed.
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_hgn_thermo_c2(cs_real_t alpha,
                 cs_real_t y,
                 cs_real_t z,
                 cs_real_t P,
                 cs_real_t v);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of the specific internal energy with respect to the
 *        volume (\f$\alpha\f$), mass (\f$y\f$) and energy (\f$z\f$) fractions,
 *        as well as the pressure and the specific volume \f$\tau\f$.
 *
 * It uses a quasi-Newton method to solve:
 * \f[
 *   \mathcal{P}(\alpha, y, z, e, \tau) = P
 * \f]
 *
 * \param[in]     alpha     the volume fraction
 * \param[in]     y         the mass fraction
 * \param[in]     z         the energy fraction
 * \param[in]     pr         the pressure
 * \param[in]     v         the specific volume
 *
 * \return specific internal energy.
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_hgn_thermo_ie(cs_real_t alpha,
                 cs_real_t y,
                 cs_real_t z,
                 cs_real_t pr,
                 cs_real_t v);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of the equilibrium fractions.
 *
 * The equilibrium fractions correspond to the definition of the mixture for
 * which one gets the pressure, temperature and chemical potential equilibrium.
 *
 * They are computed by using a Dichotomy algorithm on the function
 * characterizing the equilibrium (two forms available).
 *
 * The search for the equilibrium point is done in plane (P,T). Dichotomy is
 * performed on the pressure along the saturation curve.
 *
 * \param[in]  e          specific internal energy
 * \param[in]  v          specific volume
 * \param[out] palpha_eq  pointer to equilibrium volume fraction
 * \param[out] py_eq      pointer to equilibrium mass fraction
 * \param[out] pz_eq      pointer to equilibrium energy fraction
 */
/*----------------------------------------------------------------------------*/

void
cs_hgn_thermo_eq(cs_real_t  e,
                 cs_real_t  v,
                 cs_real_t *palpha_eq,
                 cs_real_t *py_eq,
                 cs_real_t *pz_eq);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_HGN_THERMO_H__ */
