/*============================================================================
 * Temporal and z-axis interpolation for meteorological profiles
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

#include "cs_intprf.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_turb_kw_model.c

  * Temporal and z-axis interpolation for meteorological profiles
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
 * \brief Temporal and z-axis interpolation for meteorological profiles
 *
 * An optimized linear interpolation is used.
 *
 * \param[in]   nprofz      total number of measure points
 * \param[in]   nproft      total number of time values
 * \param[in]   profz       z coordinates of measure points
 * \param[in]   proft       physical times of dataset acquisition
 * \param[in]   profv       measured values
 * \param[in]   xz          interpolation elevation
 * \param[in]   t           interpolation time
 *
 * \return  interpolated value
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_intprf(int              nprofz,
          int              nproft,
          const cs_real_t  profz[],
          const cs_real_t  proft[],
          const cs_real_t  profv[],
          cs_real_t        xz,
          cs_real_t        t)
{
  /* Time interpolation
     ------------------ */

  int it, it1, it2;
  cs_real_t alphat;
  if (t <= proft[0]) {
    it1 = 0;
    it2 = 0;
    alphat = 1.;
  }
  else if (t >= proft[nproft - 1]) {
    it1 = nproft - 1;
    it2 = nproft - 1;
    alphat = 1.;
  }
  else {  /* else nproft > 1 */
    it = 0;
    while (t > proft[it + 1]) {
      it++;
    }
    it1 = it;
    it2 = it + 1;
    alphat = (proft[it2] - t)/(proft[it2] - proft[it1]);
  }

  /* Z interpolation
     --------------- */

  int iz, iz1, iz2;
  cs_real_t alphaz;
  if (xz <= profz[0]) {
    iz1 = 0;
    iz2 = 0;
    alphaz = 1.;
  }
  else if (xz >= profz[nprofz - 1]) {
    iz1 = nprofz - 1;
    iz2 = nprofz - 1;
    alphaz = 1.;
  }
  else { /* else nprofz > 1 */
    iz = 0;
    while (xz > profz[iz + 1]) {
      iz++;
    }
    iz1 = iz;
    iz2 = iz + 1;
    alphaz = (profz[iz2] - xz)/(profz[iz2] - profz[iz1]);
  }

  /* Interpolation
     ------------- */

  cs_real_t var1 =         alphaz *profv[it1*nprofz + iz1]
                   + (1. - alphaz)*profv[it1*nprofz + iz2];

  cs_real_t var2 =         alphaz *profv[it2*nprofz + iz1]
                   + (1. - alphaz)*profv[it2*nprofz + iz2];

  cs_real_t var = alphat*var1 + (1. - alphat)*var2;

  return var;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Z-axis interpolation for meteorological profiles
 *
 * An optimized linear interpolation is used.
 *
 * \param[in]   nprofz      total number of measure points
 * \param[in]   profz       z coordinates of measure points
 * \param[in]   profv       measured values
 * \param[in]   xz          interpolation elevation
 * \param[out]  z_lv        index of lower and upel level (optional, size 2),
 *                          or NULL
 * \param[out]  var         interpolated value
 */
/*----------------------------------------------------------------------------*/

void
cs_intprz(int               nprofz,
          const cs_real_t   profz[],
          const cs_real_t   profv[],
          cs_real_t         xz,
          int              *z_lv,
          cs_real_t        *var)
{
  int iz0 = 0, iz1 = 0;

  /* Z interpolation
     --------------- */

  int iz;
  cs_real_t alphaz;
  if (xz <= profz[0]) {
    alphaz = 1.;
  }
  else if (xz >= profz[nprofz - 1]) {
    iz0 = nprofz - 1;
    iz1 = nprofz - 1;
    alphaz = 1.;
  }
  else { /* else nprofz > 1 */
    iz = 0;
    while (xz > profz[iz + 1]) {
      iz++;
    }
    iz0 = iz;
    iz1 = iz + 1;
    alphaz = (profz[iz1] - xz)/(profz[iz1] - profz[iz0]);
  }

  /* Interpolation
     ------------- */

  if (z_lv != NULL) {
    z_lv[0] = iz0;
    z_lv[1] = iz1;
  }

  *var = alphaz*profv[iz0] + (1. - alphaz)*profv[iz1];
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
