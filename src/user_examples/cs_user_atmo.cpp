/*============================================================================
 * User-defined functions specific to amospheric flow models.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_atmo.cpp
 *
 * \brief User-defined functions specific to amospheric flow models.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Fill in vertical profiles of atmospheric properties prior to solving
 *        1D radiative transfers.
 *
 * \param[in, out] preray        pressure vertical profile
 * \param[in, out] temray        real temperature vertical profile
 * \param[in, out] romray        density vertical profile
 * \param[in, out] qvray         water vapor content vertical profile
 * \param[in, out] qlray         water liquid content vertical profile
 * \param[in, out] ncray         droplets density vertical profile
 * \param[in, out] aeroso        aerosol concentration vertical profile
 */
/*----------------------------------------------------------------------------*/

void
cs_user_atmo_1d_rad_prf(cs_real_t   preray[],
                        cs_real_t   temray[],
                        cs_real_t   romray[],
                        cs_real_t   qvray[],
                        cs_real_t   qlray[],
                        cs_real_t   ncray[],
                        cs_real_t   aeroso[])
{
  /*! [humid_aerosols_atmo] */

  cs_atmo_option_t *at_opt = cs_glob_atmo_option;
  const cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();

  cs_real_t *zvert = at_opt->rad_1d_z;
  cs_real_t *zray  = at_opt->rad_1d_zray;

  const int kvert = at_opt->rad_1d_nlevels;
  const int kmx   = at_opt->rad_1d_nlevels_max;

  const cs_real_t rvsra = phys_pro->rvsra;
  const cs_real_t rair = phys_pro->r_pg_cnst;
  const cs_real_t gz = cs_glob_physical_constants->gravity[2];
  const cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

  aeroso[0] = 10.0;

  for (int k = 1; k < kvert; k++) {
    zray[k] = zvert[k];

    const cs_real_t rhum = rair * (1.0 + (rvsra-1.0)*qvray[k]);
    const cs_real_t tmean = 0.50 * (temray[k-1] + temray[k]) + tkelvi;
    const cs_real_t rap = -fabs(gz) * (zray[k]-zray[k-1]) / rhum / tmean;
    preray[k] = preray[k-1] * exp(rap);

    if (zray[k] < 50.0) {
      aeroso[k] = aeroso[0];
    }
    else {
      aeroso[k] = aeroso[0]*exp(-(zray[k]-50.0) / 1.25e3);
      if (aeroso[k] < 5.0)
        aeroso[k] = 5.0;
    }
  }

  /* Filling the additional levels above meshed domain
     (at these levels, pressure, temperature, density profiles have been
     initialized with standard atmosphere profiles) */

  for (int k = kvert; k < kmx; k++) {
    zray[k] = zvert[k];

    /* read meteo data for temperature, water wapor and liquid content in
       upper layers for example to fill temray, qvray, qlray */

    const cs_real_t rhum = rair*(1.0+(rvsra-1.0)*qvray[k]);
    const cs_real_t tmean = 0.50*(temray[k-1]+temray[k]) + tkelvi;
    const cs_real_t rap = -fabs(gz)*(zray[k]-zray[k-1]) / rhum / tmean;
    preray[k] = preray[k-1]*exp(rap);
    romray[k] = preray[k] / (temray[k]+tkelvi) / rhum;

    /* nc not known above the meshed domain
       droplets radius is assumed of mean volume radius = 5 microns */
    ncray[k]
      = 1.e-6*(3.0*romray[k]*qlray[k])/(4.0*cs_math_pi*1.e3*pow(5.e-6, 3.0));

    // similarly, aerosol concentration not known
    aeroso[k] = aeroso[0]*exp(-(zray[k]-50.0) / 1.25e3);
    if (aeroso[k] < 5.0)
      aeroso[k] = 5.0;
  }

  /*! [humid_aerosols_atmo] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
