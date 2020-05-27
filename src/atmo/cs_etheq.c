/*============================================================================
 *  Compute etheta and eq variable knowing the saturation.
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

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

#include "cs_air_props.h"
#include "cs_atmo.h"
#include "cs_base.h"
#include "cs_math.h"
#include "cs_physical_constants.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_etheq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_etheq.c

  Compute etheta and eq variable knowing the saturation.
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Compute etheta and eq variable knowing the saturation.
 *
 * \param[in]       pphy        pressure [Pa]
 * \param[in]       thetal      Liquid potential temperature
 * \param[in]       qw          total water amount
 * \param[in]       qldia       mean liquid water content
 * \param[in]       xnebdia     nebulosity after the diagnostic
 * \param[in]       xnn         second order moment "n" <s'ql'>/(2*sigma_s**2)
 * \param[out]      etheta      sensible heat part of buoyancy flux
 * \param[out]      eq          latent heat part of buoyancy flux
 */
/*----------------------------------------------------------------------------*/

void
cs_etheq(cs_real_t   pphy,
         cs_real_t   thetal,
         cs_real_t   qw,
         cs_real_t   qldia,
         cs_real_t   xnebdia,
         cs_real_t   xnn,
         cs_real_t  *etheta,
         cs_real_t  *eq)
{
  const cs_fluid_properties_t *phys_pro = cs_glob_fluid_properties;

  cs_real_t rair = phys_pro->r_pg_cnst;
  cs_real_t p0 = phys_pro->p0;
  cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

  /* Most probable and simplest case
     =============================== */

  cs_real_t rvsra = phys_pro->rvsra;
  if (qldia <= 0. || cs_glob_atmo_option->subgrid_model == 0) {
    *etheta = 1.;
    *eq     = (rvsra - 1.)*thetal;
    return;
  }

  /* Initialization
     ============== */

  cs_real_t cp = phys_pro->cp0;
  cs_real_t clatev = phys_pro->clatev;

  cs_real_t lscp = clatev/cp;

  *etheta = 1.;
  *eq     = (rvsra - 1.)*thetal;

  /* Sub-grid model
     ============== */

  /* standard k-eps */

  /* Compute "liquid" temperature */
  cs_real_t tl = thetal*pow(p0/pphy, -rair/cp);

  /* Constants */
  cs_real_t beta1 = cs_math_pow2(clatev)/(cp*(rair*rvsra)*cs_math_pow2(tl));

  /* Compute saturation ql */
  cs_real_t qsl = cs_air_yw_sat(tl-tkelvi, pphy);

  cs_real_t a = 1./(1. + beta1*qsl);

  cs_real_t alpha1 = qsl*beta1*pow(pphy/p0, rair/cp)/lscp;

  /* Potential temperature (dry) */
  cs_real_t theta = thetal + (clatev/cp)*pow(p0/pphy, rair/cp)*qldia;

  /* Compute "thermodynamic" temperature */
  cs_real_t t = theta*pow(p0/pphy, -rair/cp);

  /* Compute q */
  cs_real_t q = qw - qldia;

  /* Computesaturation q */
  cs_real_t qs = cs_air_yw_sat(t-tkelvi, pphy);

  /* Calcul de de */
  cs_real_t de = (clatev/cp)*pow(p0/pphy, rair/cp) - rvsra*theta;

  /* Constants for moddis = 3 */
  cs_real_t beta2  = cs_math_pow2(clatev)/(cp*(rair*rvsra)*cs_math_pow2(t));
  cs_real_t aa     = 1./(1. + beta2*qs);
  cs_real_t alpha2 = qs*(beta2*cp/clatev)*pow(pphy/p0, rair/cp);

  /* New computation of d for moddis = 3 */
  cs_real_t dd =   (clatev/cp)*pow(p0/pphy, rair/cp)
                 * (1. + (rvsra - 1. )*q - qldia) - rvsra*theta;

  /* Nebulosity */
  cs_real_t rn = xnebdia;

  if (cs_glob_atmo_option->subgrid_model == 1) {
    /* Bechtold et al. 1995 */
    *etheta = 1. - a*alpha1*de*xnn;
    *eq     = (rvsra - 1.)*theta + a*de*xnn;
  }
  else if (cs_glob_atmo_option->subgrid_model == 2) {
    /* Bouzereau et al. 2004 */
    *etheta = 1. + (rvsra - 1.)*q - qldia - a*alpha1*dd*xnn;
    *eq     = (rvsra - 1.)*theta + a*dd*xnn;
  }
  else if (cs_glob_atmo_option->subgrid_model == 3) {
    /* Cuijpers et Duynkerke 1993, etc.
     * Linear interpolation between saturated and non-saturated cases
     * (coefficient of partial nebulosity r) */
    *etheta = 1. + (rvsra - 1.)*q - rn*(qldia + aa*alpha2*dd);
    *eq     = (rvsra - 1.)*theta + rn*aa*dd;
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
