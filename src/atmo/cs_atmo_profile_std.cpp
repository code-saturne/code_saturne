/*============================================================================
 * Compute standard atmospheric profile
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

#include "bft/bft_mem.h"

#include "atmo/cs_air_props.h"
#include "base/cs_base.h"
#include "base/cs_math.h"
#include "base/cs_physical_constants.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "atmo/cs_atmo_profile_std.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_atmo_profile_std.cpp

  Compute standard atmospheric profile.
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
/*!
 * \brief compute standard atmospheric profile (Holton p 374)
 *
 * Note:
 * - standard pressure is 101325Pa
 * - standard temperature is 288.15K
 * - standard density is 1.225 kg/m^3
 *
 * \param[in]       z_ref      reference altitude in m
 * \param[in]       p_ref      reference pressure
 * \param[in]       t_ref      reference temperature
 * \param[in]       z          absolute altitude in m
 * \param[out]      p          pressure in Pa
 * \param[out]      t          temperature in K
 * \param[out]      r          density in kg/m3
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_profile_std(cs_real_t   z_ref,
                    cs_real_t   p_ref,
                    cs_real_t   t_ref,
                    cs_real_t   z,
                    cs_real_t  *p,
                    cs_real_t  *t,
                    cs_real_t  *r)
{
  cs_real_t rair = cs_glob_fluid_properties->r_pg_cnst;
  cs_real_t zt = 11000.;
  cs_real_t g = cs_math_3_norm(cs_glob_physical_constants->gravity);
  cs_real_t a = -6.5e-3;

  if (z <= zt) {
    *t = t_ref + a*(z-z_ref);
    *p = p_ref*pow((*t)/t_ref, -g/rair/a);
    *r = (*p)/rair/(*t);
  }
  else {
    cs_real_t t11 = t_ref + a*(zt - z_ref);
    cs_real_t p11 = p_ref*pow(t11/t_ref,-g/rair/a);
    *t = t11;
    *p = p11*exp(-g/rair/t11*(z - zt));
    *r = (*p)/rair/(*t);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
