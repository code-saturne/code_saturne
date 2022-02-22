/*============================================================================
 * Compute standard atmospheric profile
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
#include "cs_base.h"
#include "cs_math.h"
#include "cs_physical_constants.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_atmo_profile_std.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_atmo_profile_std.c

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
 * \param[in]       z          absolute altitude in m
 * \param[out]      p          pressure in Pa
 * \param[out]      t          temperature in K
 * \param[out]      r          density in kg/m3
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_profile_std(cs_real_t   z,
                    cs_real_t  *p,
                    cs_real_t  *t,
                    cs_real_t  *r)
{
  cs_real_t rair = cs_glob_fluid_properties->r_pg_cnst;
  cs_real_t p0 = cs_glob_fluid_properties->p0;
  cs_real_t t0 = cs_glob_fluid_properties->t0;
  cs_real_t zt = 11000.;
  cs_real_t g = cs_math_3_norm(cs_glob_physical_constants->gravity);
  cs_real_t a = -6.5e-3;

  if (z <= zt){
    *t = t0 + a*z;
    *p = p0*pow((*t)/t0, -g/rair/a);
    *r = (*p)/rair/(*t);
  }
  else {
    cs_real_t t11 = t0 + a*zt;
    cs_real_t p11 = p0*pow(t11/t0,-g/rair/a);
    *t = t11;
    *p = p11*exp(-g/rair/t11*(z - zt));
    *r = (*p)/rair/(*t);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
