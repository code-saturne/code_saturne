/*============================================================================
 * Compute properties for water with Freesteam
 *============================================================================*/

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

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute properties with Freesteam in a defined thermal plane.
 *
 * \param[in]   thermo_plane  thermodynamic plane
 * \param[in]   property      property queried
 * \param[in]   n_vals        number of values
 * \param[in]   var1          values on first plane axis
 * \param[in]   var2          values on second plane axis
 * \param[out]  val           resulting property values
 *
 * \return  floating point value associated with the key id for this field
 */
/*----------------------------------------------------------------------------*/

void
cs_phys_prop_freesteam(cs_phys_prop_thermo_plane_type_t   thermo_plane,
                       cs_phys_prop_type_t                property,
                       const cs_lnum_t                    n_vals,
                       const cs_real_t                    var1[],
                       const cs_real_t                    var2[],
                       cs_real_t                          val[])
{
#if defined(HAVE_FREESTEAM)
  if (thermo_plane == CS_PHYS_PROP_PLANE_PH) {
    for (cs_lnum_t i = 0; i < n_vals; i++) {
      SteamState S0 = freesteam_set_ph(var1[i], var2[i]);
      switch (property) {
      case CS_PHYS_PROP_PRESSURE:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice: you chose to work in the %s plane."), "ph");
        break;
      case CS_PHYS_PROP_TEMPERATURE:
        val[i] = freesteam_T(S0);
        break;
      case CS_PHYS_PROP_ENTHALPY:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice: you chose to work in the %s plane."), "ph");
        break;
      case CS_PHYS_PROP_ENTROPY:
        val[i] = freesteam_s(S0);
        break;
      case CS_PHYS_PROP_ISOBARIC_HEAT_CAPACITY:
        val[i] = freesteam_cp(S0);
        break;
      case CS_PHYS_PROP_ISOCHORIC_HEAT_CAPACITY:
        val[i] = freesteam_cv(S0);
        break;
      case CS_PHYS_PROP_SPECIFIC_VOLUME:
        val[i] = freesteam_v(S0);
        break;
      case CS_PHYS_PROP_DENSITY:
        val[i] = freesteam_rho(S0);
        break;
      case CS_PHYS_PROP_INTERNAL_ENERGY:
        val[i] = freesteam_u(S0);
        break;
      case CS_PHYS_PROP_QUALITY:
        val[i] = freesteam_x(S0);
        break;
      case CS_PHYS_PROP_THERMAL_CONDUCTIVITY:
        val[i] = freesteam_k(S0);
        break;
      case CS_PHYS_PROP_DYNAMIC_VISCOSITY:
        val[i] = freesteam_mu(S0);
        break;
      case CS_PHYS_PROP_SPEED_OF_SOUND:
        val[i] = freesteam_w(S0);
        break;
      }
    }
  }
  else if (thermo_plane == CS_PHYS_PROP_PLANE_PT) {
    for (cs_lnum_t i = 0; i < n_vals; i++) {
      SteamState S0 = freesteam_set_pT(var1[i], var2[i]);
      switch (property) {
      case CS_PHYS_PROP_PRESSURE:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice: you chose to work in the %s plane."), "pT");
        break;
      case CS_PHYS_PROP_TEMPERATURE:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice: you chose to work in the %s plane."), "pT");
        break;
      case CS_PHYS_PROP_ENTHALPY:
        val[i] = freesteam_h(S0);
        break;
      case CS_PHYS_PROP_ENTROPY:
        val[i] = freesteam_s(S0);
        break;
      case CS_PHYS_PROP_ISOBARIC_HEAT_CAPACITY:
        val[i] = freesteam_cp(S0);
        break;
      case CS_PHYS_PROP_ISOCHORIC_HEAT_CAPACITY:
        val[i] = freesteam_cv(S0);
        break;
      case CS_PHYS_PROP_SPECIFIC_VOLUME:
        val[i] = freesteam_v(S0);
        break;
      case CS_PHYS_PROP_DENSITY:
        val[i] = freesteam_rho(S0);
        break;
      case CS_PHYS_PROP_INTERNAL_ENERGY:
        val[i] = freesteam_u(S0);
        break;
      case CS_PHYS_PROP_QUALITY:
        val[i] = freesteam_x(S0);
        break;
      case CS_PHYS_PROP_THERMAL_CONDUCTIVITY:
        val[i] = freesteam_k(S0);
        break;
      case CS_PHYS_PROP_DYNAMIC_VISCOSITY:
        val[i] = freesteam_mu(S0);
        break;
      case CS_PHYS_PROP_SPEED_OF_SOUND:
        val[i] = freesteam_w(S0);
        break;
      }
    }
  }
  else if (thermo_plane == CS_PHYS_PROP_PLANE_PS) {
    for (cs_lnum_t i = 0; i < n_vals; i++) {
      SteamState S0 = freesteam_set_ps(var1[i], var2[i]);
      switch (property) {
      case CS_PHYS_PROP_PRESSURE:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice: you chose to work in the %s plane."), "ps");
        break;
      case CS_PHYS_PROP_TEMPERATURE:
        val[i] = freesteam_T(S0);
        break;
      case CS_PHYS_PROP_ENTHALPY:
        val[i] = freesteam_h(S0);
        break;
      case CS_PHYS_PROP_ENTROPY:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice: you chose to work in the %s plane."), "ps");
        break;
      case CS_PHYS_PROP_ISOBARIC_HEAT_CAPACITY:
        val[i] = freesteam_cp(S0);
        break;
      case CS_PHYS_PROP_ISOCHORIC_HEAT_CAPACITY:
        val[i] = freesteam_cv(S0);
        break;
      case CS_PHYS_PROP_SPECIFIC_VOLUME:
        val[i] = freesteam_v(S0);
        break;
      case CS_PHYS_PROP_DENSITY:
        val[i] = freesteam_rho(S0);
        break;
      case CS_PHYS_PROP_INTERNAL_ENERGY:
        val[i] = freesteam_u(S0);
        break;
      case CS_PHYS_PROP_QUALITY:
        val[i] = freesteam_x(S0);
        break;
      case CS_PHYS_PROP_THERMAL_CONDUCTIVITY:
        val[i] = freesteam_k(S0);
        break;
      case CS_PHYS_PROP_DYNAMIC_VISCOSITY:
        val[i] = freesteam_mu(S0);
        break;
      case CS_PHYS_PROP_SPEED_OF_SOUND:
        val[i] = freesteam_w(S0);
        break;
      }
    }
  }
  else if (thermo_plane == CS_PHYS_PROP_PLANE_PU) {
    for (cs_lnum_t i = 0; i < n_vals; i++) {
      SteamState S0 = freesteam_set_pu(var1[i], var2[i]);
      switch (property) {
      case CS_PHYS_PROP_PRESSURE:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice: you chose to work in the %s plane."), "pu");
        break;
      case CS_PHYS_PROP_TEMPERATURE:
        val[i] = freesteam_T(S0);
        break;
      case CS_PHYS_PROP_ENTHALPY:
        val[i] = freesteam_h(S0);
        break;
      case CS_PHYS_PROP_ENTROPY:
        val[i] = freesteam_s(S0);
        break;
      case CS_PHYS_PROP_ISOBARIC_HEAT_CAPACITY:
        val[i] = freesteam_cp(S0);
        break;
      case CS_PHYS_PROP_ISOCHORIC_HEAT_CAPACITY:
        val[i] = freesteam_cv(S0);
        break;
      case CS_PHYS_PROP_SPECIFIC_VOLUME:
        val[i] = freesteam_v(S0);
        break;
      case CS_PHYS_PROP_DENSITY:
        val[i] = freesteam_rho(S0);
        break;
      case CS_PHYS_PROP_INTERNAL_ENERGY:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice: you chose to work in the %s plane."), "pu");
        break;
      case CS_PHYS_PROP_QUALITY:
        val[i] = freesteam_x(S0);
        break;
      case CS_PHYS_PROP_THERMAL_CONDUCTIVITY:
        val[i] = freesteam_k(S0);
        break;
      case CS_PHYS_PROP_DYNAMIC_VISCOSITY:
        val[i] = freesteam_mu(S0);
        break;
      case CS_PHYS_PROP_SPEED_OF_SOUND:
        val[i] = freesteam_w(S0);
        break;
      }
    }
  }
  else if (thermo_plane == CS_PHYS_PROP_PLANE_PV) {
    for (cs_lnum_t i = 0; i < n_vals; i++) {
      SteamState S0 = freesteam_set_pv(var1[i], var2[i]);
      switch (property) {
      case CS_PHYS_PROP_PRESSURE:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice: you chose to work in the %s plane."), "pv");
        break;
      case CS_PHYS_PROP_TEMPERATURE:
        val[i] = freesteam_T(S0);
        break;
      case CS_PHYS_PROP_ENTHALPY:
        val[i] = freesteam_h(S0);
        break;
      case CS_PHYS_PROP_ENTROPY:
        val[i] = freesteam_s(S0);
        break;
      case CS_PHYS_PROP_ISOBARIC_HEAT_CAPACITY:
        val[i] = freesteam_cp(S0);
        break;
      case CS_PHYS_PROP_ISOCHORIC_HEAT_CAPACITY:
        val[i] = freesteam_cv(S0);
        break;
      case CS_PHYS_PROP_SPECIFIC_VOLUME:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice: you chose to work in the %s plane."), "pv");
        break;
      case CS_PHYS_PROP_DENSITY:
        val[i] = freesteam_rho(S0);
        break;
      case CS_PHYS_PROP_INTERNAL_ENERGY:
        val[i] = freesteam_u(S0);
        break;
      case CS_PHYS_PROP_QUALITY:
        val[i] = freesteam_x(S0);
        break;
      case CS_PHYS_PROP_THERMAL_CONDUCTIVITY:
        val[i] = freesteam_k(S0);
        break;
      case CS_PHYS_PROP_DYNAMIC_VISCOSITY:
        val[i] = freesteam_mu(S0);
        break;
      case CS_PHYS_PROP_SPEED_OF_SOUND:
        val[i] = freesteam_w(S0);
        break;
      }
    }
  }
  else if (thermo_plane == CS_PHYS_PROP_PLANE_TS) {
    for (cs_lnum_t i = 0; i < n_vals; i++) {
      SteamState S0 = freesteam_set_Ts(var1[i], var2[i]);
      switch (property) {
      case CS_PHYS_PROP_PRESSURE:
        val[i] = freesteam_p(S0);
        break;
      case CS_PHYS_PROP_TEMPERATURE:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice: you chose to work in the %s plane."), "Ts");
        break;
      case CS_PHYS_PROP_ENTHALPY:
        val[i] = freesteam_h(S0);
        break;
      case CS_PHYS_PROP_ENTROPY:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice: you chose to work in the %s plane."), "Ts");
        break;
      case CS_PHYS_PROP_ISOBARIC_HEAT_CAPACITY:
        val[i] = freesteam_cp(S0);
        break;
      case CS_PHYS_PROP_ISOCHORIC_HEAT_CAPACITY:
        val[i] = freesteam_cv(S0);
        break;
      case CS_PHYS_PROP_SPECIFIC_VOLUME:
        val[i] = freesteam_v(S0);
        break;
      case CS_PHYS_PROP_DENSITY:
        val[i] = freesteam_rho(S0);
        break;
      case CS_PHYS_PROP_INTERNAL_ENERGY:
        val[i] = freesteam_u(S0);
        break;
      case CS_PHYS_PROP_QUALITY:
        val[i] = freesteam_x(S0);
        break;
      case CS_PHYS_PROP_THERMAL_CONDUCTIVITY:
        val[i] = freesteam_k(S0);
        break;
      case CS_PHYS_PROP_DYNAMIC_VISCOSITY:
        val[i] = freesteam_mu(S0);
        break;
      case CS_PHYS_PROP_SPEED_OF_SOUND:
        val[i] = freesteam_w(S0);
        break;
      }
    }
  }
  else if (thermo_plane == CS_PHYS_PROP_PLANE_TX) {
    for (cs_lnum_t i = 0; i < n_vals; i++) {
      SteamState S0 = freesteam_set_Tx(var1[i], var2[i]);
      switch (property) {
      case CS_PHYS_PROP_PRESSURE:
        val[i] = freesteam_p(S0);
        break;
      case CS_PHYS_PROP_TEMPERATURE:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice: you chose to work in the %s plane."), "Tx");
        break;
      case CS_PHYS_PROP_ENTHALPY:
        val[i] = freesteam_h(S0);
        break;
      case CS_PHYS_PROP_ENTROPY:
        val[i] = freesteam_s(S0);
        break;
      case CS_PHYS_PROP_ISOBARIC_HEAT_CAPACITY:
        val[i] = freesteam_cp(S0);
        break;
      case CS_PHYS_PROP_ISOCHORIC_HEAT_CAPACITY:
        val[i] = freesteam_cv(S0);
        break;
      case CS_PHYS_PROP_SPECIFIC_VOLUME:
        val[i] = freesteam_v(S0);
        break;
      case CS_PHYS_PROP_DENSITY:
        val[i] = freesteam_rho(S0);
        break;
      case CS_PHYS_PROP_INTERNAL_ENERGY:
        val[i] = freesteam_u(S0);
        break;
      case CS_PHYS_PROP_QUALITY:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice: you chose to work in the %s plane."), "Tx");
        break;
      case CS_PHYS_PROP_THERMAL_CONDUCTIVITY:
        val[i] = freesteam_k(S0);
        break;
      case CS_PHYS_PROP_DYNAMIC_VISCOSITY:
        val[i] = freesteam_mu(S0);
        break;
      case CS_PHYS_PROP_SPEED_OF_SOUND:
        val[i] = freesteam_w(S0);
        break;
      }
    }
  }
#else
  bft_error(__FILE__, __LINE__, 0,
            _("Freesteam support not available in this build."));
#endif
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
