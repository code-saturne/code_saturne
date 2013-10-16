/*============================================================================
 * Compute properties for water with Freesteam
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2013 EDF S.A.

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
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_physical_properties.h"

#if defined(HAVE_FREESTEAM)
#include <freesteam/steam_ph.h>
#include <freesteam/steam_pT.h>
#include <freesteam/steam_ps.h>
#include <freesteam/steam_pu.h>
#include <freesteam/steam_pv.h>
#include <freesteam/steam_Ts.h>
#include <freesteam/steam_Tx.h>

#include <freesteam/region1.h>
#include <freesteam/region2.h>
#include <freesteam/region3.h>
#include <freesteam/region4.h>
#endif

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute properties with freeteam in a defined thermal plane.
 *----------------------------------------------------------------------------*/

void
CS_PROCF (csfspp, CSFSPP)(cs_thermo_plane_type_t *thermo_plane,
                          cs_property_type_t *property,
                          const    int *const ncel,
                          double *var1,
                          double *var2,
                          double *val)
{
#if defined(HAVE_FREESTEAM)
  if (*thermo_plane == PLANE_PH)
    for (int i = 0; i < *ncel; i++) {
      SteamState S0 = freesteam_set_ph(var1[i], var2[i]);
      switch (*property) {
      case PRESSURE:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice : you choose to work in %s plane\n"), "ph");
        break;
      case TEMPERATURE:
        val[i] = freesteam_T(S0);
        break;
      case ENTHALPY:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice : you choose to work in %s plane\n"), "ph");
        break;
      case ENTROPY:
        val[i] = freesteam_s(S0);
        break;
      case ISOBARIC_HEAT_CAPACITY:
        val[i] = freesteam_cp(S0);
        break;
      case ISOCHORIC_HEAT_CAPACITY:
        val[i] = freesteam_cv(S0);
        break;
      case SPECIFIC_VOLUME:
        val[i] = freesteam_v(S0);
        break;
      case DENSITY:
        val[i] = freesteam_rho(S0);
        break;
      case INTERNAL_ENERGY:
        val[i] = freesteam_u(S0);
        break;
      case QUALITY:
        val[i] = freesteam_x(S0);
        break;
      case THERMAL_CONDUCTIVITY:
        val[i] = freesteam_k(S0);
        break;
      case DYNAMIC_VISCOSITY:
        val[i] = freesteam_mu(S0);
        break;
      case SPEED_OF_SOUND:
        val[i] = freesteam_w(S0);
        break;
      }
    }
  else if (*thermo_plane == PLANE_PT)
    for (int i = 0; i < *ncel; i++) {
      SteamState S0 = freesteam_set_pT(var1[i], var2[i]);
      switch (*property) {
      case PRESSURE:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice : you choose to work in %s plane\n"), "pT");
        break;
      case TEMPERATURE:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice : you choose to work in %s plane\n"), "pT");
        break;
      case ENTHALPY:
        val[i] = freesteam_h(S0);
        break;
      case ENTROPY:
        val[i] = freesteam_s(S0);
        break;
      case ISOBARIC_HEAT_CAPACITY:
        val[i] = freesteam_cp(S0);
        break;
      case ISOCHORIC_HEAT_CAPACITY:
        val[i] = freesteam_cv(S0);
        break;
      case SPECIFIC_VOLUME:
        val[i] = freesteam_v(S0);
        break;
      case DENSITY:
        val[i] = freesteam_rho(S0);
        break;
      case INTERNAL_ENERGY:
        val[i] = freesteam_u(S0);
        break;
      case QUALITY:
        val[i] = freesteam_x(S0);
        break;
      case THERMAL_CONDUCTIVITY:
        val[i] = freesteam_k(S0);
        break;
      case DYNAMIC_VISCOSITY:
        val[i] = freesteam_mu(S0);
        break;
      case SPEED_OF_SOUND:
        val[i] = freesteam_w(S0);
        break;
      }
    }
  else if (*thermo_plane == PLANE_PS)
    for (int i = 0; i < *ncel; i++) {
      SteamState S0 = freesteam_set_ps(var1[i], var2[i]);
      switch (*property) {
      case PRESSURE:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice : you choose to work in %s plane\n"), "ps");
        break;
      case TEMPERATURE:
        val[i] = freesteam_T(S0);
        break;
      case ENTHALPY:
        val[i] = freesteam_h(S0);
        break;
      case ENTROPY:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice : you choose to work in %s plane\n"), "ps");
        break;
      case ISOBARIC_HEAT_CAPACITY:
        val[i] = freesteam_cp(S0);
        break;
      case ISOCHORIC_HEAT_CAPACITY:
        val[i] = freesteam_cv(S0);
        break;
      case SPECIFIC_VOLUME:
        val[i] = freesteam_v(S0);
        break;
      case DENSITY:
        val[i] = freesteam_rho(S0);
        break;
      case INTERNAL_ENERGY:
        val[i] = freesteam_u(S0);
        break;
      case QUALITY:
        val[i] = freesteam_x(S0);
        break;
      case THERMAL_CONDUCTIVITY:
        val[i] = freesteam_k(S0);
        break;
      case DYNAMIC_VISCOSITY:
        val[i] = freesteam_mu(S0);
        break;
      case SPEED_OF_SOUND:
        val[i] = freesteam_w(S0);
        break;
      }
    }
  else if (*thermo_plane == PLANE_PU)
    for (int i = 0; i < *ncel; i++) {
      SteamState S0 = freesteam_set_pu(var1[i], var2[i]);
      switch (*property) {
      case PRESSURE:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice : you choose to work in %s plane\n"), "pu");
        break;
      case TEMPERATURE:
        val[i] = freesteam_T(S0);
        break;
      case ENTHALPY:
        val[i] = freesteam_h(S0);
        break;
      case ENTROPY:
        val[i] = freesteam_s(S0);
        break;
      case ISOBARIC_HEAT_CAPACITY:
        val[i] = freesteam_cp(S0);
        break;
      case ISOCHORIC_HEAT_CAPACITY:
        val[i] = freesteam_cv(S0);
        break;
      case SPECIFIC_VOLUME:
        val[i] = freesteam_v(S0);
        break;
      case DENSITY:
        val[i] = freesteam_rho(S0);
        break;
      case INTERNAL_ENERGY:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice : you choose to work in %s plane\n"), "pu");
        break;
      case QUALITY:
        val[i] = freesteam_x(S0);
        break;
      case THERMAL_CONDUCTIVITY:
        val[i] = freesteam_k(S0);
        break;
      case DYNAMIC_VISCOSITY:
        val[i] = freesteam_mu(S0);
        break;
      case SPEED_OF_SOUND:
        val[i] = freesteam_w(S0);
        break;
      }
    }
  else if (*thermo_plane == PLANE_PV)
    for (int i = 0; i < *ncel; i++) {
      SteamState S0 = freesteam_set_pv(var1[i], var2[i]);
      switch (*property) {
      case PRESSURE:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice : you choose to work in %s plane\n"), "pv");
        break;
      case TEMPERATURE:
        val[i] = freesteam_T(S0);
        break;
      case ENTHALPY:
        val[i] = freesteam_h(S0);
        break;
      case ENTROPY:
        val[i] = freesteam_s(S0);
        break;
      case ISOBARIC_HEAT_CAPACITY:
        val[i] = freesteam_cp(S0);
        break;
      case ISOCHORIC_HEAT_CAPACITY:
        val[i] = freesteam_cv(S0);
        break;
      case SPECIFIC_VOLUME:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice : you choose to work in %s plane\n"), "pv");
        break;
      case DENSITY:
        val[i] = freesteam_rho(S0);
        break;
      case INTERNAL_ENERGY:
        val[i] = freesteam_u(S0);
        break;
      case QUALITY:
        val[i] = freesteam_x(S0);
        break;
      case THERMAL_CONDUCTIVITY:
        val[i] = freesteam_k(S0);
        break;
      case DYNAMIC_VISCOSITY:
        val[i] = freesteam_mu(S0);
        break;
      case SPEED_OF_SOUND:
        val[i] = freesteam_w(S0);
        break;
      }
    }
  else if (*thermo_plane == PLANE_TS)
    for (int i = 0; i < *ncel; i++) {
      SteamState S0 = freesteam_set_Ts(var1[i], var2[i]);
      switch (*property) {
      case PRESSURE:
        val[i] = freesteam_p(S0);
        break;
      case TEMPERATURE:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice : you choose to work in %s plane\n"), "Ts");
        break;
      case ENTHALPY:
        val[i] = freesteam_h(S0);
        break;
      case ENTROPY:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice : you choose to work in %s plane\n"), "Ts");
        break;
      case ISOBARIC_HEAT_CAPACITY:
        val[i] = freesteam_cp(S0);
        break;
      case ISOCHORIC_HEAT_CAPACITY:
        val[i] = freesteam_cv(S0);
        break;
      case SPECIFIC_VOLUME:
        val[i] = freesteam_v(S0);
        break;
      case DENSITY:
        val[i] = freesteam_rho(S0);
        break;
      case INTERNAL_ENERGY:
        val[i] = freesteam_u(S0);
        break;
      case QUALITY:
        val[i] = freesteam_x(S0);
        break;
      case THERMAL_CONDUCTIVITY:
        val[i] = freesteam_k(S0);
        break;
      case DYNAMIC_VISCOSITY:
        val[i] = freesteam_mu(S0);
        break;
      case SPEED_OF_SOUND:
        val[i] = freesteam_w(S0);
        break;
      }
    }
  else if (*thermo_plane == PLANE_TX)
    for (int i = 0; i < *ncel; i++) {
      SteamState S0 = freesteam_set_Tx(var1[i], var2[i]);
      switch (*property) {
      case PRESSURE:
        val[i] = freesteam_p(S0);
        break;
      case TEMPERATURE:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice : you choose to work in %s plane\n"), "Tx");
        break;
      case ENTHALPY:
        val[i] = freesteam_h(S0);
        break;
      case ENTROPY:
        val[i] = freesteam_s(S0);
        break;
      case ISOBARIC_HEAT_CAPACITY:
        val[i] = freesteam_cp(S0);
        break;
      case ISOCHORIC_HEAT_CAPACITY:
        val[i] = freesteam_cv(S0);
        break;
      case SPECIFIC_VOLUME:
        val[i] = freesteam_v(S0);
        break;
      case DENSITY:
        val[i] = freesteam_rho(S0);
        break;
      case INTERNAL_ENERGY:
        val[i] = freesteam_u(S0);
        break;
      case QUALITY:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice : you choose to work in %s plane\n"), "Tx");
        break;
      case THERMAL_CONDUCTIVITY:
        val[i] = freesteam_k(S0);
        break;
      case DYNAMIC_VISCOSITY:
        val[i] = freesteam_mu(S0);
        break;
      case SPEED_OF_SOUND:
        val[i] = freesteam_w(S0);
        break;
      }
    }
#else
        bft_error(__FILE__, __LINE__, 0, _("You try to use freesteam table without it\n"));
#endif
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
