/*============================================================================
 * Steady laminar flamelet gas combustion model.
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

/*----------------------------------------------------------------------------
 * Standard library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_array.h"
#include "base/cs_array_reduce.h"
#include "base/cs_base.h"
#include "base/cs_boundary_conditions.h"
#include "base/cs_dispatch.h"
#include "cdo/cs_equation_param.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_operator.h"
#include "base/cs_field_pointer.h"
#include "base/cs_physical_constants.h"
#include "base/cs_log.h"
#include "turb/cs_les_filter.h"
#include "base/cs_math.h"
#include "base/cs_parall.h"
#include "base/cs_physical_properties.h"
#include "base/cs_restart.h"
#include "base/cs_restart_default.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "rayt/cs_rad_transfer.h"
#include "turb/cs_turbulence_model.h"

#include "pprt/cs_combustion_model.h"
#include "pprt/cs_physical_model.h"

#include "cogz/cs_combustion_gas.h"
#include "cogz/cs_combustion_boundary_conditions.h"
#include "cogz/cs_combustion_ht_convert.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cogz/cs_combustion_slfm.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_combustion_slfm.cpp
        Steady laminar flamelet gas combustion model.
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief interpolate the phiscal properties in SLFM
 *
 * \param[in]     update_rad  indicator of whether update radiation properties
 * \param[in]     zm          mean value of mixture fraction
 * \param[in]     zvar        variance of mixture fraction
 * \param[in]     kim         mean value of scalar dissipation rate
 * \param[in]     xrm         mean value of heat loss
 * \param[in,out] phim        pointer to filtered physical properties output
 * \param[in,out] rad_work    pointer to filtered radiation properties output
 */
/*----------------------------------------------------------------------------*/

static void
_filtered_physical_prop(const bool        update_rad,
                        const cs_real_t   zm,
                        const cs_real_t   zvar,
                        const cs_real_t   kim,
                        const cs_real_t   xrm,
                        cs_real_t        *phim,
                        cs_real_t        *rad_work)
{
  assert(phim != nullptr);
  assert(rad_work != nullptr);

  const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;
  const cs_rad_transfer_params_t *rt_params = cs_glob_rad_transfer_params;

  const int nlibvar = cm->nlibvar;
  const int nzm = cm->nzm;
  const int nzvar = cm->nzvar;
  const int nxr = cm->nxr;
  const int nki = cm->nki;
  const int nwsgg = rt_params->nwsgg;

  const int flamelet_zm = cm->flamelet_zm;
  const int flamelet_zvar = cm->flamelet_zvar;
  const int flamelet_ki = cm->flamelet_ki;
  const int flamelet_xr = cm->flamelet_xr;

  const cs_real_t *flamelet_library = cm->flamelet_library;
  const cs_real_t *radiation_library = cm->radiation_library;

  if (rt_params->type > CS_RAD_TRANSFER_NONE)
    assert(radiation_library != nullptr);

  cs_real_t *xdata = nullptr;
  cs_real_t *phi_2 = nullptr, *phi_3 = nullptr, *phi_4 = nullptr;
  cs_real_t *rad_work2 = nullptr, *rad_work3 = nullptr, *rad_work4 = nullptr;

  CS_MALLOC(phi_4, nzvar*nki*nxr*nlibvar, cs_real_t);
  CS_MALLOC(phi_3,       nki*nxr*nlibvar, cs_real_t);
  CS_MALLOC(phi_2,           nxr*nlibvar, cs_real_t);

  if (update_rad) {
    CS_MALLOC(rad_work4, nzvar*nki*nxr*nwsgg*2, cs_real_t);
    CS_MALLOC(rad_work3,       nki*nxr*nwsgg*2, cs_real_t);
    CS_MALLOC(rad_work2,           nxr*nwsgg*2, cs_real_t);
  }

  CS_MALLOC(xdata, nzm, cs_real_t);

  // Interpolate mixture fraction
  for (int iz = 0; iz < nzm; iz++)
    xdata[iz] = flamelet_library[iz*nzvar*nki*nxr*nlibvar + flamelet_zm];

  int n_elements = nzvar*nki*nxr*nlibvar, n_elements_rad = nzvar*nki*nxr*nwsgg*2;
  if (xdata[0] >= zm) {
    for (int i = 0; i < n_elements; i++)
      phi_4[i] = flamelet_library[i];

    if (update_rad)
      for (int i = 0; i < n_elements_rad; i++)
        rad_work4[i] = radiation_library[i];
  }
  else if (xdata[nzm - 1] <= zm) {
    for (int i = 0; i < n_elements; i++)
      phi_4[i] = flamelet_library[(nzm - 1)*n_elements + i];
    if (update_rad)
      for (int i = 0; i < n_elements_rad; i++)
        rad_work4[i] = radiation_library[(nzm - 1)*n_elements_rad + i];
  }
  else {
    int data_idx = 0;
    while (   xdata[data_idx] < zm
           && data_idx < nzm - 1) {
      data_idx++;
    }

    if (data_idx > 0) {
      cs_real_t weight_zm = (zm - xdata[data_idx - 1])/(xdata[data_idx] - xdata[data_idx - 1]);
      for (int i = 0; i < n_elements; i++)
        phi_4[i] = (1.0 - weight_zm)*flamelet_library[(data_idx - 1)*n_elements + i]
                 + weight_zm*flamelet_library[data_idx*n_elements + i];

      if (update_rad)
        for (int i = 0; i < n_elements_rad; i++)
          rad_work4[i] = (1.0 - weight_zm)*radiation_library[(data_idx - 1)*n_elements_rad + i]
                       + weight_zm*radiation_library[data_idx*n_elements_rad + i];
    }
  }

  // Interpolate mixture fraction variance

  CS_REALLOC(xdata, nzvar, cs_real_t);

  for (int izvar = 0; izvar < nzvar; izvar++)
    xdata[izvar] = phi_4[izvar*nki*nxr*nlibvar + flamelet_zvar];

  n_elements = nki*nxr*nlibvar, n_elements_rad = nki*nxr*nwsgg*2;
  if (xdata[0] >= zvar) {
    for (int i = 0; i < n_elements; i++)
      phi_3[i] = phi_4[i];
    if (update_rad)
      for (int i = 0; i < n_elements_rad; i++)
        rad_work3[i] = rad_work4[i];
  }
  else if (xdata[nzvar - 1] <= zvar) {
    for (int i = 0; i < n_elements; i++)
      phi_3[i] = phi_4[(nzvar - 1)*n_elements + i];
    if (update_rad)
      for (int i = 0; i < n_elements_rad; i++)
        rad_work3[i] = rad_work4[(nzvar - 1)*n_elements_rad + i];
  }
  else {
    int data_idx = 0;
    while (   xdata[data_idx] < zvar
           && data_idx < nzvar - 1) {
      data_idx++;
    }

    if (data_idx > 0) {
      cs_real_t weight_zvar = (zvar - xdata[data_idx - 1])/(xdata[data_idx] - xdata[data_idx - 1]);
      for (int i = 0; i < n_elements; i++)
        phi_3[i] = (1.0 - weight_zvar)*phi_4[(data_idx - 1)*n_elements + i]
                 + weight_zvar*phi_4[data_idx*n_elements + i];

      if (update_rad)
        for (int i = 0; i < n_elements_rad; i++)
          rad_work3[i] = (1.0 - weight_zvar)*rad_work4[(data_idx - 1)*n_elements_rad + i]
                       + weight_zvar*rad_work4[data_idx*n_elements_rad + i];
    }
  }

  // Interpolate scalar dissipation rate

  CS_REALLOC(xdata, nki, cs_real_t);

  for (int iki = 0; iki < nki; iki++)
    xdata[iki] = phi_3[iki*nxr*nlibvar + flamelet_ki];

  n_elements = nxr*nlibvar, n_elements_rad = nxr*nwsgg*2;
  if (xdata[0] >= kim) {
    for (int i = 0; i < n_elements; i++)
      phi_2[i] = phi_3[i];
    if (update_rad)
      for (int i = 0; i < n_elements_rad; i++)
        rad_work2[i] = rad_work3[i];
  }
  else if (xdata[nki - 1] <= kim) {
    for (int i = 0; i < n_elements; i++)
      phi_2[i] = phi_3[(nki - 1)*n_elements + i];
    if (update_rad)
      for (int i = 0; i < n_elements_rad; i++)
        rad_work2[i] = rad_work3[(nki - 1)*n_elements_rad + i];
  }
  else {
    int data_idx = 0;
    while (   xdata[data_idx] < kim
           && data_idx < nki - 1) {
      data_idx++;
    }

    if (data_idx > 0) {
      cs_real_t weight_kim = (kim - xdata[data_idx - 1])/(xdata[data_idx] - xdata[data_idx - 1]);
      for (int i = 0; i < n_elements; i++)
        phi_2[i] = (1.0 - weight_kim)*phi_3[(data_idx - 1)*n_elements + i]
                 + weight_kim*phi_3[data_idx*n_elements + i];

      if (update_rad)
        for (int i = 0; i < n_elements_rad; i++)
          rad_work2[i] = (1.0 - weight_kim)*rad_work3[(data_idx - 1)*n_elements_rad + i]
                       + weight_kim*rad_work3[data_idx*n_elements_rad + i];
    }
  }

  // Interpolate heat loss

  CS_REALLOC(xdata, nxr, cs_real_t);

  for (int ixr = 0; ixr < nxr; ixr++)
    xdata[ixr] = phi_2[ixr*nlibvar + flamelet_xr];

  n_elements = nlibvar, n_elements_rad = nwsgg*2;
  if (xdata[0] >= xrm) {
    for (int i = 0; i < n_elements; i++)
      phim[i] = phi_2[i];
    if (update_rad)
      for (int i = 0; i < n_elements_rad; i++)
        rad_work[i] = rad_work2[i];
  }
  else if (xdata[nxr - 1] <= xrm) {
    for (int i = 0; i < n_elements; i++)
      phim[i] = phi_2[(nxr - 1)*n_elements + i];
    if (update_rad)
      for (int i = 0; i < n_elements_rad; i++)
        rad_work[i] = rad_work2[(nxr - 1)*n_elements_rad + i];
  }
  else {
    int data_idx = 0;
    while (   xdata[data_idx] < xrm
           && data_idx < nxr - 1) {
      data_idx++;
    }

    if (data_idx > 0) {
      cs_real_t weight_xrm = (xrm - xdata[data_idx - 1])/(xdata[data_idx] - xdata[data_idx - 1]);
      for (int i = 0; i < n_elements; i++)
        phim[i] = (1.0 - weight_xrm)*phi_2[(data_idx - 1)*n_elements + i]
                + weight_xrm*phi_2[data_idx*n_elements + i];

      if (update_rad)
        for (int i = 0; i < n_elements_rad; i++)
          rad_work[i] = (1.0 - weight_xrm)*rad_work2[(data_idx - 1)*n_elements_rad + i]
                      + weight_xrm*rad_work2[data_idx*n_elements_rad + i];
    }
  }

  CS_FREE(xdata);
  if (update_rad) {
    CS_FREE(rad_work2);
    CS_FREE(rad_work3);
    CS_FREE(rad_work4);
  }
  CS_FREE(phi_2);
  CS_FREE(phi_3);
  CS_FREE(phi_4);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief interpolate the density
 *
 * \param[in]     zm          mean value of mixture fraction
 * \param[in]     zvar        variance of mixture fraction
 * \param[in]     kim         mean value of scalar dissipation rate
 * \param[in]     xrm         mean value of heat loss
 * \param[in,out] rhom        filtered density output
 */
/*----------------------------------------------------------------------------*/

static void
_filtered_density(const cs_real_t   zm,
                  const cs_real_t   zvar,
                  const cs_real_t   kim,
                  const cs_real_t   xrm,
                  cs_real_t        &rhom)
{

  const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  const int nzm = cm->nzm;
  const int nzvar = cm->nzvar;
  const int nxr = cm->nxr;
  const int nki = cm->nki;

  const cs_real_t *rho_library = cm->rho_library;

  cs_real_t *xdata = nullptr;
  cs_real_t *phi_2 = nullptr, *phi_3 = nullptr, *phi_4 = nullptr;

  // Only select zm, zvar, xr, rho, ki from the library
  const int n_var_local = 5;

  CS_MALLOC(phi_4, nzvar*nki*nxr*n_var_local, cs_real_t);
  CS_MALLOC(phi_3,       nki*nxr*n_var_local, cs_real_t);
  CS_MALLOC(phi_2,           nxr*n_var_local, cs_real_t);

  const int rho_zm_idx = 0;
  const int rho_zvar_idx = 1;
  const int rho_xr_idx = 2;
  const int rho_idx = 3;
  const int rho_ki_idx = 4;

  CS_MALLOC(xdata, nzm, cs_real_t);

  // Interpolate mixture fraction
  for (int iz = 0; iz < nzm; iz++)
    xdata[iz] = rho_library[iz*nzvar*nki*nxr*n_var_local + rho_zm_idx];

  int n_elements = nzvar*nki*nxr*n_var_local;
  if (xdata[0] >= zm) {
    for (int i = 0; i < n_elements; i++)
      phi_4[i] = rho_library[i];
  }
  else if (xdata[nzm - 1] <= zm) {
    for (int i = 0; i < n_elements; i++)
      phi_4[i] = rho_library[(nzm - 1)*n_elements + i];
  }
  else {
    int data_idx = 0;
    while (   xdata[data_idx] < zm
           && data_idx < nzm - 1) {
      data_idx++;
    }

    if (data_idx > 0) {
      cs_real_t weight_zm = (zm - xdata[data_idx - 1])/(xdata[data_idx] - xdata[data_idx - 1]);
      for (int i = 0; i < n_elements; i++)
        phi_4[i] = (1.0 - weight_zm)*rho_library[(data_idx - 1)*n_elements + i]
                 + weight_zm*rho_library[data_idx*n_elements + i];
    }
  }

  // Interpolate mixture fraction variance

  CS_REALLOC(xdata, nzvar, cs_real_t);

  for (int izvar = 0; izvar < nzvar; izvar++)
    xdata[izvar] = phi_4[izvar*nki*nxr*n_var_local + rho_zvar_idx];

  n_elements = nki*nxr*n_var_local;
  if (xdata[0] >= zvar) {
    for (int i = 0; i < n_elements; i++)
      phi_3[i] = phi_4[i];
  }
  else if (xdata[nzvar - 1] <= zvar) {
    for (int i = 0; i < n_elements; i++)
      phi_3[i] = phi_4[(nzvar - 1)*n_elements + i];
  }
  else {
    int data_idx = 0;
    while (   xdata[data_idx] < zvar
           && data_idx < nzvar - 1) {
      data_idx++;
    }

    if (data_idx > 0) {
      cs_real_t weight_zvar = (zvar - xdata[data_idx - 1])/(xdata[data_idx] - xdata[data_idx - 1]);
      for (int i = 0; i < n_elements; i++)
        phi_3[i] = (1.0 - weight_zvar)*phi_4[(data_idx - 1)*n_elements + i]
                 + weight_zvar*phi_4[data_idx*n_elements + i];
    }
  }

  // Interpolate scalar dissipation rate

  CS_REALLOC(xdata, nki, cs_real_t);

  for (int iki = 0; iki < nki; iki++)
    xdata[iki] = phi_3[iki*nxr*n_var_local + rho_ki_idx];

  n_elements = nxr*n_var_local;
  if (xdata[0] >= kim) {
    for (int i = 0; i < n_elements; i++)
      phi_2[i] = phi_3[i];
  }
  else if (xdata[nki - 1] <= kim) {
    for (int i = 0; i < n_elements; i++)
      phi_2[i] = phi_3[(nki - 1)*n_elements + i];
  }
  else {
    int data_idx = 0;
    while (   xdata[data_idx] < kim
           && data_idx < nki - 1) {
      data_idx++;
    }

    if (data_idx > 0) {
      cs_real_t weight_kim = (kim - xdata[data_idx - 1])/(xdata[data_idx] - xdata[data_idx - 1]);
      for (int i = 0; i < n_elements; i++)
        phi_2[i] = (1.0 - weight_kim)*phi_3[(data_idx - 1)*n_elements + i]
                 + weight_kim*phi_3[data_idx*n_elements + i];
    }
  }

  // Interpolate heat loss

  CS_REALLOC(xdata, nxr, cs_real_t);

  for (int ixr = 0; ixr < nxr; ixr++)
    xdata[ixr] = phi_2[ixr*n_var_local + rho_xr_idx];

  n_elements = n_var_local;
  if (xdata[0] >= xrm) {
    rhom = phi_2[rho_idx];
  }
  else if (xdata[nxr - 1] <= xrm) {
    rhom = phi_2[(nxr - 1)*n_elements + rho_idx];
  }
  else {
    int data_idx = 0;
    while (   xdata[data_idx] < xrm
           && data_idx < nxr - 1) {
      data_idx++;
    }

    if (data_idx > 0) {
      cs_real_t weight_xrm = (xrm - xdata[data_idx - 1])/(xdata[data_idx] - xdata[data_idx - 1]);
      rhom = (1.0 - weight_xrm)*phi_2[(data_idx - 1)*n_elements + rho_idx]
           + weight_xrm*phi_2[data_idx*n_elements + rho_idx];
    }
  }

  CS_FREE(xdata);
  CS_FREE(phi_2);
  CS_FREE(phi_3);
  CS_FREE(phi_4);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief interpolate the phiscal properties in FPV model
 *
 * \param[in]     update_rad  indicator of whether update radiation properties
 * \param[in]     zm          mean value of mixture fraction
 * \param[in]     zvar        variance of mixture fraction
 * \param[in]     progm       mean value of progress variable
 * \param[in]     xrm         mean value of heat loss
 * \param[in,out] phim        pointer to filtered physical properties output
 * \param[in,out] rad_work    pointer to filtered radiation properties output
 */
/*----------------------------------------------------------------------------*/

static void
_filtered_physical_prop_progvar(const bool        update_rad,
                                const cs_real_t   zm,
                                const cs_real_t   zvar,
                                const cs_real_t   progm,
                                const cs_real_t   xrm,
                                cs_real_t        *phim,
                                cs_real_t        *rad_work)
{
  assert(phim != nullptr);
  assert(rad_work != nullptr);

  const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;
  const cs_rad_transfer_params_t *rt_params = cs_glob_rad_transfer_params;

  const int nlibvar = cm->nlibvar;
  const int nzm = cm->nzm;
  const int nzvar = cm->nzvar;
  const int nxr = cm->nxr;
  const int nki = cm->nki;
  const int nwsgg = rt_params->nwsgg;

  const int flamelet_zm = cm->flamelet_zm;
  const int flamelet_zvar = cm->flamelet_zvar;
  const int flamelet_c = cm->flamelet_c;
  const int flamelet_xr = cm->flamelet_xr;

  const cs_real_t *flamelet_library = cm->flamelet_library;
  const cs_real_t *radiation_library = cm->radiation_library;

  if (rt_params->type > CS_RAD_TRANSFER_NONE)
    assert(radiation_library != nullptr);

  cs_real_t *xdata = nullptr;
  cs_real_t *phi_2 = nullptr, *phi_3 = nullptr, *phi_4 = nullptr;
  cs_real_t *rad_work2 = nullptr, *rad_work3 = nullptr, *rad_work4 = nullptr;

  CS_MALLOC(phi_4, nzvar*nxr*nki*nlibvar, cs_real_t);
  CS_MALLOC(phi_3,       nxr*nki*nlibvar, cs_real_t);
  CS_MALLOC(phi_2,           nki*nlibvar, cs_real_t);

  if (update_rad) {
    CS_MALLOC(rad_work4, nzvar*nxr*nki*nwsgg*2, cs_real_t);
    CS_MALLOC(rad_work3,       nxr*nki*nwsgg*2, cs_real_t);
    CS_MALLOC(rad_work2,           nki*nwsgg*2, cs_real_t);
  }

  CS_MALLOC(xdata, nzm, cs_real_t);

  // Interpolate mixture fraction
  for (int iz = 0; iz < nzm; iz++)
    xdata[iz] = flamelet_library[iz*nzvar*nxr*nki*nlibvar + flamelet_zm];

  int n_elements = nzvar*nxr*nki*nlibvar, n_elements_rad = nzvar*nxr*nki*nwsgg*2;
  if (xdata[0] >= zm) {
    for (int i = 0; i < n_elements; i++)
      phi_4[i] = flamelet_library[i];

    if (update_rad)
      for (int i = 0; i < n_elements_rad; i++)
        rad_work4[i] = radiation_library[i];
  }
  else if (xdata[nzm - 1] <= zm) {
    for (int i = 0; i < n_elements; i++)
      phi_4[i] = flamelet_library[(nzm - 1)*n_elements + i];
    if (update_rad)
      for (int i = 0; i < n_elements_rad; i++)
        rad_work4[i] = radiation_library[(nzm - 1)*n_elements_rad + i];
  }
  else {
    int data_idx = 0;
    while (   xdata[data_idx] < zm
           && data_idx < nzm - 1) {
      data_idx++;
    }

    if (data_idx > 0) {
      cs_real_t weight_zm = (zm - xdata[data_idx - 1])/(xdata[data_idx] - xdata[data_idx - 1]);
      for (int i = 0; i < n_elements; i++)
        phi_4[i] = (1.0 - weight_zm)*flamelet_library[(data_idx - 1)*n_elements + i]
                 + weight_zm*flamelet_library[data_idx*n_elements + i];

      if (update_rad)
        for (int i = 0; i < n_elements_rad; i++)
          rad_work4[i] = (1.0 - weight_zm)*radiation_library[(data_idx - 1)*n_elements_rad + i]
                       + weight_zm*radiation_library[data_idx*n_elements_rad + i];
    }
  }

  // Interpolate mixture fraction variance

  CS_REALLOC(xdata, nzvar, cs_real_t);

  for (int izvar = 0; izvar < nzvar; izvar++)
    xdata[izvar] = phi_4[izvar*nxr*nki*nlibvar + flamelet_zvar];

  n_elements = nxr*nki*nlibvar, n_elements_rad = nxr*nki*nwsgg*2;
  if (xdata[0] >= zvar) {
    for (int i = 0; i < n_elements; i++)
      phi_3[i] = phi_4[i];
    if (update_rad)
      for (int i = 0; i < n_elements_rad; i++)
        rad_work3[i] = rad_work4[i];
  }
  else if (xdata[nzvar - 1] <= zvar) {
    for (int i = 0; i < n_elements; i++)
      phi_3[i] = phi_4[(nzvar - 1)*n_elements + i];
    if (update_rad)
      for (int i = 0; i < n_elements_rad; i++)
        rad_work3[i] = rad_work4[(nzvar - 1)*n_elements_rad + i];
  }
  else {
    int data_idx = 0;
    while (   xdata[data_idx] < zvar
           && data_idx < nzvar - 1) {
      data_idx++;
    }

    if (data_idx > 0) {
      cs_real_t weight_zvar = (zvar - xdata[data_idx - 1])/(xdata[data_idx] - xdata[data_idx - 1]);
      for (int i = 0; i < n_elements; i++)
        phi_3[i] = (1.0 - weight_zvar)*phi_4[(data_idx - 1)*n_elements + i]
                 + weight_zvar*phi_4[data_idx*n_elements + i];

      if (update_rad)
        for (int i = 0; i < n_elements_rad; i++)
          rad_work3[i] = (1.0 - weight_zvar)*rad_work4[(data_idx - 1)*n_elements_rad + i]
                       + weight_zvar*rad_work4[data_idx*n_elements_rad + i];
    }
  }

  // Interpolate heat loss depending on scalar disspation rate
  // One has to interpolate them in a coupled manner when considering
  // FPV, due to strong coupling between these two coordinates

  CS_REALLOC(xdata, nxr*nki, cs_real_t);

  cs_real_t *weight_xrm = nullptr;
  CS_MALLOC(weight_xrm, nki, cs_real_t);

  for (int ixr = 0; ixr < nxr; ixr++)
    for (int iki = 0; iki < nki; iki++)
      xdata[ixr*nki + iki] = phi_3[ixr*nki*nlibvar + iki*nlibvar + flamelet_xr];

  for (int iki = 0; iki < nki; iki++) {

    if (xdata[iki] >= xrm) {
      for (int i = 0; i < nlibvar; i++)
        phi_2[iki*nlibvar + i] = phi_3[iki*nlibvar + i];
      if (update_rad)
        for (int i = 0; i < 2*nwsgg; i++)
          rad_work2[iki*2*nwsgg + i] = rad_work3[iki*2*nwsgg + i];
    }
    else if (xdata[(nxr - 1)*nki + iki] <= xrm) {
      for (int i = 0; i < nlibvar; i++)
        phi_2[iki*nlibvar + i] = phi_3[(nxr - 1)*nki*nlibvar + iki*nlibvar + i];
      if (update_rad)
        for (int i = 0; i < 2*nwsgg; i++)
          rad_work2[iki*2*nwsgg + i] =
            rad_work3[(nxr - 1)*nki*2*nwsgg + iki*2*nwsgg + i];
    }
    else {
      int data_idx = 0;
      while (   xdata[data_idx*nki + iki] < xrm
             && data_idx < nxr - 1) {
        data_idx++;
      }

      if (data_idx > 0) {
        weight_xrm[iki] = (xrm - xdata[(data_idx - 1)*nki + iki])
          / (xdata[data_idx*nki + iki] - xdata[(data_idx - 1)*nki + iki]);
        for (int i = 0; i < nlibvar; i++)
          phi_2[iki*nlibvar + i] =
            (1.0 - weight_xrm[iki])*phi_3[(data_idx - 1)*nki*nlibvar + iki*nlibvar + i]
            + weight_xrm[iki]*phi_3[data_idx*nki*nlibvar + iki*nlibvar + i];
        if (update_rad)
          for (int i = 0; i < 2*nwsgg; i++)
            rad_work2[iki*2*nwsgg + i] =
              (1.0 - weight_xrm[iki])*rad_work3[(data_idx - 1)*nki*2*nwsgg + iki*nwsgg + i]
              + weight_xrm[iki]*rad_work3[data_idx*nki*2*nwsgg + iki*2*nwsgg + i];
      }
    }
  }

  CS_FREE(weight_xrm);

  // Interpolate progress variable

  CS_REALLOC(xdata, nki, cs_real_t);

  for (int iki = 0; iki < nki; iki++)
    xdata[iki] = phi_2[iki*nlibvar + flamelet_c];

  n_elements = nlibvar, n_elements_rad = nwsgg*2;
  if (xdata[0] >= progm) {
    for (int i = 0; i < n_elements; i++)
      phim[i] = phi_2[i];
    if (update_rad)
      for (int i = 0; i < n_elements_rad; i++)
        rad_work[i] = rad_work2[i];
  }
  else if (xdata[nxr - 1] <= progm) {
    for (int i = 0; i < n_elements; i++)
      phim[i] = phi_2[(nki - 1)*n_elements + i];
    if (update_rad)
      for (int i = 0; i < n_elements_rad; i++)
        rad_work[i] = rad_work2[(nki - 1)*n_elements_rad + i];
  }
  else {
    int data_idx = 0;
    while (   xdata[data_idx] < progm
           && data_idx < nki - 1) {
      data_idx++;
    }

    if (data_idx > 0) {
      cs_real_t weight_c = (progm - xdata[data_idx - 1])/(xdata[data_idx] - xdata[data_idx - 1]);
      for (int i = 0; i < n_elements; i++)
        phim[i] = (1.0 - weight_c)*phi_2[(data_idx - 1)*n_elements + i]
                + weight_c*phi_2[data_idx*n_elements + i];

      if (update_rad)
        for (int i = 0; i < n_elements_rad; i++)
          rad_work[i] = (1.0 - weight_c)*rad_work2[(data_idx - 1)*n_elements_rad + i]
                      + weight_c*rad_work2[data_idx*n_elements_rad + i];
    }
  }

  CS_FREE(xdata);
  if (update_rad) {
    CS_FREE(rad_work2);
    CS_FREE(rad_work3);
    CS_FREE(rad_work4);
  }
  CS_FREE(phi_2);
  CS_FREE(phi_3);
  CS_FREE(phi_4);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief interpolate the density in the FPV model
 *
 * \param[in]     zm          mean value of mixture fraction
 * \param[in]     zvar        variance of mixture fraction
 * \param[in]     progm       mean value of progress variable
 * \param[in]     xrm         mean value of heat loss
 * \param[in,out] rhom        filtered density output
 */
/*----------------------------------------------------------------------------*/

static void
_filtered_density_progvar(const cs_real_t   zm,
                          const cs_real_t   zvar,
                          const cs_real_t   progm,
                          const cs_real_t   xrm,
                          cs_real_t        &rhom)
{

  const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  const int nzm = cm->nzm;
  const int nzvar = cm->nzvar;
  const int nxr = cm->nxr;
  const int nki = cm->nki;

  const cs_real_t *rho_library = cm->rho_library;

  cs_real_t *xdata = nullptr;
  cs_real_t *phi_2 = nullptr, *phi_3 = nullptr, *phi_4 = nullptr;

  // Only select zm, zvar, xr, rho, ki from the library
  const int n_var_local = 5;

  CS_MALLOC(phi_4, nzvar*nxr*nki*n_var_local, cs_real_t);
  CS_MALLOC(phi_3,       nxr*nki*n_var_local, cs_real_t);
  CS_MALLOC(phi_2,           nki*n_var_local, cs_real_t);

  const int rho_zm_idx = 0;
  const int rho_zvar_idx = 1;
  const int rho_xr_idx = 2;
  const int rho_idx = 3;
  const int rho_c_idx = 4;

  CS_MALLOC(xdata, nzm, cs_real_t);

  // Interpolate mixture fraction
  for (int iz = 0; iz < nzm; iz++)
    xdata[iz] = rho_library[iz*nzvar*nxr*nki*n_var_local + rho_zm_idx];

  int n_elements = nzvar*nxr*nki*n_var_local;
  if (xdata[0] >= zm) {
    for (int i = 0; i < n_elements; i++)
      phi_4[i] = rho_library[i];
  }
  else if (xdata[nzm - 1] <= zm) {
    for (int i = 0; i < n_elements; i++)
      phi_4[i] = rho_library[(nzm - 1)*n_elements + i];
  }
  else {
    int data_idx = 0;
    while (   xdata[data_idx] < zm
           && data_idx < nzm - 1) {
      data_idx++;
    }

    if (data_idx > 0) {
      cs_real_t weight_zm = (zm - xdata[data_idx - 1])/(xdata[data_idx] - xdata[data_idx - 1]);
      for (int i = 0; i < n_elements; i++)
        phi_4[i] = (1.0 - weight_zm)*rho_library[(data_idx - 1)*n_elements + i]
                 + weight_zm*rho_library[data_idx*n_elements + i];
    }
  }

  // Interpolate mixture fraction variance

  CS_REALLOC(xdata, nzvar, cs_real_t);

  for (int izvar = 0; izvar < nzvar; izvar++)
    xdata[izvar] = phi_4[izvar*nxr*nki*n_var_local + rho_zvar_idx];

  n_elements = nxr*nki*n_var_local;
  if (xdata[0] >= zvar) {
    for (int i = 0; i < n_elements; i++)
      phi_3[i] = phi_4[i];
  }
  else if (xdata[nzvar - 1] <= zvar) {
    for (int i = 0; i < n_elements; i++)
      phi_3[i] = phi_4[(nzvar - 1)*n_elements + i];
  }
  else {
    int data_idx = 0;
    while (   xdata[data_idx] < zvar
           && data_idx < nzvar - 1) {
      data_idx++;
    }

    if (data_idx > 0) {
      cs_real_t weight_zvar = (zvar - xdata[data_idx - 1])/(xdata[data_idx] - xdata[data_idx - 1]);
      for (int i = 0; i < n_elements; i++)
        phi_3[i] = (1.0 - weight_zvar)*phi_4[(data_idx - 1)*n_elements + i]
                 + weight_zvar*phi_4[data_idx*n_elements + i];
    }
  }

  // Interpolate heat loss depending on scalar disspation rate
  // One has to interpolate them in a coupled manner when considering
  // FPV, due to strong coupling between these two coordinates

  CS_REALLOC(xdata, nxr*nki, cs_real_t);

  cs_real_t *weight_xrm = nullptr;
  CS_MALLOC(weight_xrm, nki, cs_real_t);

  for (int ixr = 0; ixr < nxr; ixr++)
    for (int iki = 0; iki < nki; iki++)
      xdata[ixr*nki + iki] = phi_3[ixr*nki*n_var_local + iki*n_var_local + rho_xr_idx];

  for (int iki = 0; iki < nki; iki++) {

    if (xdata[iki] >= xrm) {
      for (int i = 0; i < n_var_local; i++)
        phi_2[iki*n_var_local + i] = phi_3[iki*n_var_local + i];
    }
    else if (xdata[(nxr - 1)*nki + iki] <= xrm) {
      for (int i = 0; i < n_var_local; i++)
        phi_2[iki*n_var_local + i] = phi_3[(nxr - 1)*nki*n_var_local + iki*n_var_local + i];
    }
    else {
      int data_idx = 0;
      while (   xdata[data_idx*nki + iki] < xrm
             && data_idx < nxr - 1) {
        data_idx++;
      }

      if (data_idx > 0) {
        weight_xrm[iki] = (xrm - xdata[(data_idx - 1)*nki + iki])
          / (xdata[data_idx*nki + iki] - xdata[(data_idx - 1)*nki + iki]);
        for (int i = 0; i < n_var_local; i++)
          phi_2[iki*n_var_local + i] =
            (1.0 - weight_xrm[iki])*phi_3[(data_idx - 1)*nki*n_var_local + iki*n_var_local + i]
            + weight_xrm[iki]*phi_3[data_idx*nki*n_var_local + iki*n_var_local + i];
        }
    }
  }

  CS_FREE(weight_xrm);

  // Interpolate progress variable

  CS_REALLOC(xdata, nki, cs_real_t);

  for (int iki = 0; iki < nki; iki++)
    xdata[iki] = phi_2[iki*n_var_local + rho_c_idx];

  n_elements = n_var_local;
  if (xdata[0] >= progm) {
    rhom = phi_2[rho_idx];
  }
  else if (xdata[nki - 1] <= progm) {
    rhom = phi_2[(nki - 1)*n_elements + rho_idx];
  }
  else {
    int data_idx = 0;
    while (   xdata[data_idx] < progm
           && data_idx < nki - 1) {
      data_idx++;
    }

    if (data_idx > 0) {
      cs_real_t weight_c = (progm - xdata[data_idx - 1])
                         / (xdata[data_idx] - xdata[data_idx - 1]);
      rhom = (1.0 - weight_c)*phi_2[(data_idx - 1)*n_elements + rho_idx]
           + weight_c*phi_2[data_idx*n_elements + rho_idx];
    }
  }

  CS_FREE(xdata);
  CS_FREE(phi_2);
  CS_FREE(phi_3);
  CS_FREE(phi_4);
}

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize specific fields for slfm gas combustion model.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_slfm_fields_init(void) {

  // Only when not a restart
  if (cs_restart_present())
    return;

  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;
  cs_real_t *cvar_fm = cm->fm->val;
  cs_array_real_set_scalar(n_cells_ext, 0.0, cvar_fm);

  if (cm->mode_fp2m == 0) {
    cs_real_t *cvar_fp2m = cm->fp2m->val;
    cs_array_real_set_scalar(n_cells_ext, 0.0, cvar_fp2m);
  }
  else {
    cs_real_t *cvar_fsqm = cm->fsqm->val;
    cs_real_t *cpro_recvr = cm->recvr->val;
    cs_array_real_set_scalar(n_cells_ext, 0.0, cvar_fsqm);
    cs_array_real_set_scalar(n_cells_ext, 0.0, cpro_recvr);
  }

  if (CS_F_(h) != nullptr) {
    cs_real_t *cvar_scalt = CS_F_(h)->val;
    cs_array_real_set_scalar(n_cells_ext, cm->hinoxy, cvar_scalt);
  }

  if (cm->type%100 >= 2) {
    cs_real_t *cvar_prog = cm->pvm->val;
    cs_real_t *cpro_prog =
      cm->ym[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES - 1]->val;
    cs_array_copy(n_cells_ext, cvar_prog, cpro_prog);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute physical properties for steady laminar flamelet model.
 *
 * \param[in]     iterns     Navier-Stokes sub-iterations indicator:
 *                           - if strictly negative, indicate that this
 *                                function is called outside Navier-Stokes loop
 *                           - if positive, Navier-Stokes iteration number.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_slfm_physical_properties(int   iterns)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  const int flamelet_temp = cm->flamelet_temp;
  const int flamelet_vis = cm->flamelet_vis;
  const int flamelet_dt = cm->flamelet_dt;
  const int flamelet_rho = cm->flamelet_rho;
  const int flamelet_c = cm->flamelet_c;
  const int flamelet_omg_c = cm->flamelet_omg_c;
  const int flamelet_hrr = cm->flamelet_hrr;
  const int flamelet_temp2 = cm->flamelet_temp2;
  const int *flamelet_species = cm->flamelet_species;

  const cs_real_t hinfue = cm->hinfue, hinoxy = cm->hinoxy;

  const cs_real_t *cvar_fm = CS_F_(fm)->val;

  cs_real_t *cvar_scalt = nullptr, *cpro_xr = nullptr, *fp2m = nullptr;
  if (cm->mode_fp2m == 0)
    fp2m = cm->fp2m->val;
  else {
    fp2m = cm->recvr->val;
    cs_combustion_reconstruct_variance(cvar_fm,
                                       cm->fsqm->val,
                                       fp2m);
  }

  cs_host_context ctx;

  if (cm->type%2 == 1) {

    cvar_scalt = CS_F_(h)->val;
    cpro_xr = cm->xr->val;

    ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
      cs_real_t had = cvar_fm[c_id]*hinfue + (1.0-cvar_fm[c_id])*hinoxy;
      cpro_xr[c_id] = cs::max(-(cvar_scalt[c_id] - had), 0.0);
    });
  }
  else {
    CS_MALLOC(cpro_xr, n_cells, cs_real_t);
    cs_array_real_fill_zero(n_cells, cpro_xr);
  }

  cs_real_t *cvar_progvar = nullptr,*cpro_omegac = nullptr,*cpro_totki = nullptr;

  if (cm->type%100 >= 2) {
    cvar_progvar = cm->pvm->val;
    cpro_omegac = cm->omgc->val;

    ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
      cs_real_t cmax, cmid, cmin;
      cs_combustion_slfm_max_mid_min_progvar(cvar_fm[c_id], &cmax, &cmid, &cmin);
      cvar_progvar[c_id] = cs::min(cvar_progvar[c_id], cmax);
    });
  }
  else {
    cpro_totki = cm->totki->val;
    cs_combustion_scalar_dissipation_rate(cm->fm,
                                          CS_F_(rho)->val_pre,
                                          fp2m,
                                          cpro_totki);
  }

  cs_real_t *cpro_temp = CS_F_(t)->val;
  cs_real_t *cpro_viscl = CS_F_(mu)->val;
  cs_real_t *cpro_tem2 = cm->t2m->val;
  cs_real_t *cpro_progvar = cm->ym[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES - 1]->val;
  cs_real_t *cpro_hrr = cm->hrr->val;

  cs_rad_transfer_params_t *rt_params = cs_glob_rad_transfer_params;
  const int nwsgg = rt_params->nwsgg;

  const cs_time_step_t *ts = cs_glob_time_step;

  bool update_rad = false;
  if (  rt_params->type > CS_RAD_TRANSFER_NONE
      &&ts->nt_cur % rt_params->time_control.interval_nt == 0) {
    update_rad = true;
  }

/*==============================================================================*/

  cs_real_t **cpro_species;
  cs_real_t **cpro_viscls;

  const int n_fields = cs_field_n_fields();
  const int keysca = cs_field_key_id("scalar_id");
  const int kivisl = cs_field_key_id("diffusivity_id");

  if (iterns < 0) { // outside the rho(y)-v-p coupling

    int nscapp = 0;
    for (int ii = 0; ii < n_fields; ii++) {

      cs_field_t *f_scal = cs_field_by_id(ii);

      if (!(f_scal->type & CS_FIELD_VARIABLE))
        continue;
      if (cs_field_get_key_int(f_scal, keysca) <= 0)
        continue;
      if (f_scal->type & CS_FIELD_USER)
        continue;

      nscapp += 1;
    }

    CS_MALLOC(cpro_viscls, nscapp, cs_real_t *);

    nscapp = 0;
    for (int i = 0; i < n_fields; i++) {

      cs_field_t *f_scal = cs_field_by_id(i);

      if (!(f_scal->type & CS_FIELD_VARIABLE))
        continue;
      if (cs_field_get_key_int(f_scal, keysca) <= 0)
        continue;
      if (f_scal->type & CS_FIELD_USER)
        continue;
      const int ifcvsl = cs_field_get_key_int(f_scal, kivisl);
      if (ifcvsl > -1)
        cpro_viscls[i] = cs_field_by_id(ifcvsl)->val;
    }

    CS_MALLOC(cpro_species, cm->n_gas_fl, cs_real_t *);
    for (int i = 0; i < cm->n_gas_fl; i++)
      cpro_species[i] = cm->ym[i]->val;

    cs_real_t **cpro_kg = nullptr;
    cs_real_t **cpro_emi = nullptr;
    cs_real_t *rad_work = nullptr;

    if (update_rad) {
      CS_MALLOC(cpro_kg, nwsgg, cs_real_t*);
      CS_MALLOC(cpro_emi, nwsgg, cs_real_t*);
      CS_MALLOC(rad_work, 2*nwsgg, cs_real_t);

      char f_name[64];

      for (int ig = 0; ig < nwsgg; ig++) {
        f_name[63]='\0';
        snprintf(f_name, 64, "spectral_absorption_coeff_%d",ig + 1);
        cpro_kg[ig] = cs_field_by_name(f_name)->val;
        f_name[63]='\0';
        snprintf(f_name, 64, "spectral_emission_%d",ig + 1);
        cpro_emi[ig] = cs_field_by_name(f_name)->val;
      }
    }

    cs_real_t* phim;
    CS_MALLOC(phim, cm->nlibvar, cs_real_t);

    ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
      if (cm->type%100 < 2) {
        _filtered_physical_prop(update_rad,
                                cvar_fm[c_id],
                                fp2m[c_id],
                                cpro_totki[c_id],
                                cpro_xr[c_id],
                                phim,
                                rad_work);
      }
      else {
        _filtered_physical_prop_progvar(update_rad,
                                        cvar_fm[c_id],
                                        fp2m[c_id],
                                        cvar_progvar[c_id],
                                        cpro_xr[c_id],
                                        phim,
                                        rad_work);
      }

      cpro_temp[c_id] = phim[flamelet_temp];
      cpro_viscl[c_id] = phim[flamelet_vis];

      if (cm->n_gas_fl > 0) {
        for (int i = 0; i < cm->n_gas_fl; i++) {
          cpro_species[i][c_id] = phim[flamelet_species[i]];
        }
      }

      if (nscapp > 0) {
        for (int i = 0; i < nscapp; i++)
          cpro_viscls[i][c_id] = phim[flamelet_dt]*phim[flamelet_rho];
      }

      if (flamelet_c > -1)
        cpro_progvar[c_id] = phim[flamelet_c];
      if (flamelet_hrr > -1)
        cpro_hrr[c_id] = phim[flamelet_hrr];
      if (flamelet_temp2 > -1)
        cpro_tem2[c_id] = phim[flamelet_temp2];

      if (cm->type%100 >= 2)
        cpro_omegac[c_id] = phim[flamelet_omg_c];

      if (update_rad) {
        for (int ig = 0; ig < nwsgg; ig++) {
          cpro_kg[ig][c_id] = rad_work[ig];
          cpro_emi[ig][c_id] = rad_work[ig+1];
        }
      }
    });

    CS_FREE(phim);
    CS_FREE(cpro_kg);
    CS_FREE(cpro_emi);
    CS_FREE(rad_work);
  }
  else { // Inside the rho(y)-v-p coupling

    cs_real_t *cpro_rho = CS_F_(rho)->val;

    cs_real_t *rho_eos = nullptr;
    CS_MALLOC(rho_eos, n_cells_ext, cs_real_t);
    cs_array_real_fill_zero(n_cells_ext, rho_eos);

    ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
      if (cm->type%100 < 2) {
        _filtered_density(cvar_fm[c_id],
                          fp2m[c_id],
                          cpro_totki[c_id],
                          cpro_xr[c_id],
                          rho_eos[c_id]);
      }
      else {
        _filtered_density_progvar(cvar_fm[c_id],
                                  fp2m[c_id],
                                  cvar_progvar[c_id],
                                  cpro_xr[c_id],
                                  rho_eos[c_id]);
      }

    });
    cs_les_filter(1, rho_eos, cpro_rho);

    CS_FREE(rho_eos);
    cs_combustion_boundary_conditions_density();

  }

  CS_FREE(cpro_viscls);
  CS_FREE(cpro_species);
  CS_FREE(cpro_xr);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute scalar dissipation rate for steady laminar flamelet model.
 *
 * \param[in]     f          pointer to the scalar field used to compute the
 *                           dissipation rate
 * \param[in]     cpro_rho   pointer to the density
 * \param[in]     fp2m       the variance associated to the scalar field
 * \param[in,out] cpro_totki pointer to the scalar dissipation rate
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_scalar_dissipation_rate(const cs_field_t   *f,
                                      const cs_real_t    *cpro_rho,
                                      const cs_real_t    *fp2m,
                                      cs_real_t          *cpro_totki)
{
  assert(f != nullptr);
  assert(cpro_rho != nullptr);
  assert(fp2m != nullptr);
  assert(cpro_totki != nullptr);

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_real_t *cell_vol = cs_glob_mesh_quantities->cell_vol;

  const int kivisl = cs_field_key_id("diffusivity_id");
  const int ifcvsl = cs_field_get_key_int(f, kivisl);
  cs_real_t *cpro_viscls = nullptr;

  if (ifcvsl > -1)
    cpro_viscls = cs_field_by_id(ifcvsl)->val;

  const int key_turb_diff = cs_field_key_id("turbulent_diffusivity_id");
  const int t_dif_id = cs_field_get_key_int(f, key_turb_diff);
  cs_real_t *cpro_turb_diff = nullptr;

  if (t_dif_id > -1)
   cpro_turb_diff = cs_field_by_id(t_dif_id)->val;

  const cs_real_t coef_k = 7.0e-2;
  cs_host_context ctx;
  if (cs_glob_turb_model->model == CS_TURB_LES_SMAGO_DYN) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
      const cs_real_t delta_les = cs_turb_xlesfl *
                                  pow(cs_turb_ales*cell_vol[c_id], cs_turb_bles);

      cs_real_t dnom = coef_k*delta_les*delta_les*cpro_rho[c_id];
      cpro_totki[c_id] =   (cpro_viscls[c_id] + cpro_turb_diff[c_id])
                         / dnom*fp2m[c_id];
    });
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reconstruct scalar variance in case of transporting 2nd order moment.
 *
 * \param[in]     fm         pointer to the mixture fraction
 * \param[in]     fsqm       pointer to the 2nd order moment of mixture fraction
 * \param[in,out] recvr      pointer to the reconstructed variance
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_reconstruct_variance(const cs_real_t   *fm,
                                   const cs_real_t   *fsqm,
                                   cs_real_t         *recvr)
{
  assert(fm != nullptr);
  assert(fsqm != nullptr);
  assert(recvr != nullptr);

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  cs_host_context ctx;
  ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
    recvr[c_id] = fsqm[c_id] - fm[c_id]*fm[c_id];
    recvr[c_id] = cs::max(cs::min(recvr[c_id], fm[c_id]*(1.0 - fm[c_id])), 0.0);
  });
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve maximal, middle and minimal values of progress variable
 * respectively on stable/unstable/mixing branches with a given mixture fraction
 *   - Neither heat loss nor variance is considered
 *
 * \param[in]     zm         mixture fraction value
 * \param[in]     cmax       maximal value of progress variable at given zm
 * \param[in]     cmid       middle value of progress variable at given zm
 * \param[in]     cmin       minimal value of progress variable at given zm
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_slfm_max_mid_min_progvar(const cs_real_t   zm,
                                       cs_real_t        *cmax,
                                       cs_real_t        *cmid,
                                       cs_real_t        *cmin)
{
  const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;
  const int nzm = cm->nzm;
  const int nzvar = cm->nzvar;
  const int nxr = cm->nxr;
  const int nlibvar = cm->nlibvar;
  const int nki = cm->nki;

  const int ikimid = cm->ikimid;
  const int flamelet_zm = cm->flamelet_zm;
  const int flamelet_c = cm->flamelet_c;

  cs_real_t *flamelet_library = cm->flamelet_library;

  cs_real_t *xdata = nullptr;
  CS_MALLOC(xdata, nzm, cs_real_t);

  // Interpolate mixture fraction
  for (int iz = 0; iz < nzm; iz++)
    xdata[iz] = flamelet_library[iz*nzvar*nxr*nki*nlibvar + flamelet_zm];

  if (xdata[0] >= zm) {
    *cmax = flamelet_library[(nki - 1)*nlibvar + flamelet_c];
    *cmid = flamelet_library[ikimid*nlibvar + flamelet_c];
    *cmin = flamelet_library[flamelet_c];
  }
  else if (xdata[nzm - 1] < nzm) {
    *cmax = flamelet_library[(nzm - 1)*nzvar*nxr*nki*nlibvar + (nki - 1)*nlibvar + flamelet_c];
    *cmid = flamelet_library[(nzm - 1)*nzvar*nxr*nki*nlibvar + ikimid*nlibvar + flamelet_c];
    *cmin = flamelet_library[(nzm - 1)*nzvar*nxr*nki*nlibvar + flamelet_c];
  }
  else {

    int data_idx = 0;
    while (   xdata[data_idx] < zm
           && data_idx < nzm - 1) {
      data_idx++;
    }

    int idx0, idx1;
    if (data_idx > 0) {
      cs_real_t weight_zm = (zm - xdata[data_idx - 1])
                          / (xdata[data_idx] - xdata[data_idx - 1]);

      idx0 = (data_idx - 1)*nzvar*nxr*nki*nlibvar + (nki - 1)*nlibvar + flamelet_c;
      idx1 = data_idx*nzvar*nxr*nki*nlibvar + (nki - 1)*nlibvar + flamelet_c;
      *cmax = (1.0 - weight_zm)*flamelet_library[idx0]
            + weight_zm*flamelet_library[idx1];

      idx0 = (data_idx - 1)*nzvar*nxr*nki*nlibvar + ikimid*nlibvar + flamelet_c;
      idx1 = data_idx*nzvar*nxr*nki*nlibvar + ikimid*nlibvar + flamelet_c;
      *cmid = (1.0 - weight_zm)*flamelet_library[idx0]
            + weight_zm*flamelet_library[idx1];

      idx0 = (data_idx - 1)*nzvar*nxr*nki*nlibvar + flamelet_c;
      idx1 = data_idx*nzvar*nxr*nki*nlibvar + flamelet_c;

      *cmin = (1.0 - weight_zm)*flamelet_library[idx0]
            + weight_zm*flamelet_library[idx1];
    }
  }

  CS_FREE(xdata);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Defines the source terms for the soot mass fraction and the precursor
 *        number for soot model of Moss et al for one time step.
 *
 *  The equations read: \f$ rovsdt \delta a = smbrs \f$
 *
 *  \f$ rovsdt \f$ et \f$ smbrs \f$ could already contain source term
 *  and don't have to be erased but incremented.
 *
 *  For stability sake, only positive terms should be add in \f$ rovsdt \f$.
 *  There is no constrain for \f$ smbrs \f$.
 *
 *  For a source term written \f$ S_{exp} + S_{imp} a \f$, source terms are:
 *           \f$ smbrs  = smbrs  + S_{exp} + S_{imp} a \f$
 *           \f$ rovsdt = rovsdt + \max(-S_{imp},0) \f$
 *
 *  Here are set \f$ rovsdt \f$ and \f$ smbrs \f$ containning \f$ \rho \Omega \f$
 *   - \f$ smbrs \f$ in \f$ kg_a.s^{-1} \f$ (ex: for velocity:
 *     \f$ kg.m.s^{-2} \f$, for temperature: \f$ kg.C.s^{-1} \f$,
 *     for enthalpy: \f$ J.s^{-1} \f$)
 *   - \f$ rovsdt \f$ en \f$ kg.s^{-1} \f$
 *
 * \param[in]      f_sc          pointer to scalar field
 * \param[in,out]  smbrs         explicit right hand side
 * \param[in,out]  rovsdt        implicit terms
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_slfm_source_terms(cs_field_t  *f_sc,
                                cs_real_t    smbrs[],
                                cs_real_t    rovsdt[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;

  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_real_t *cell_vol = mq_g->cell_vol;
  const cs_real_t *cell_f_vol = fvq->cell_vol;

  const int *bc_type = cs_glob_bc_type;

  const int kivisl  = cs_field_key_id("diffusivity_id");
  const int key_turb_diff = cs_field_key_id("turbulent_diffusivity_id");

  /* Initialization
     -------------- */

  /* Coef. of SGS kinetic energy used for the variance dissipation computation */
  const cs_real_t coef_k = 7.0e-2;

  cs_field_t *f_fm = CS_F_(fm);

  cs_real_t *scal_pre = f_sc->val_pre;
  cs_equation_param_t *eqp_sc = cs_field_get_equation_param(f_sc);

  const int ifcvsl = cs_field_get_key_int(f_sc, kivisl);
  const cs_real_t *viscls = NULL;
  if (ifcvsl >= 0)
    viscls = cs_field_by_id(ifcvsl)->val;

  /* Logging
     ------- */

  if (eqp_sc->verbosity >= 1)
    cs_log_printf(CS_LOG_DEFAULT,
                  _("Specific physics source terms for the variable %s\n"),
                  f_sc->name);

  cs_field_t *f_recvr = cs_field_by_name_try("reconstructed_fp2m");
  cs_field_t *f_fp2m = CS_F_(fp2m);

  cs_real_t *turb_diff = NULL;
  if (cs_glob_turb_model->model == 41) {
    /* Retrieve turbulent diffusivity value for the mixture fraction */
    const int t_dif_id = cs_field_get_key_int(f_fm, key_turb_diff);
    if (t_dif_id > -1)
      turb_diff = cs_field_by_id(t_dif_id)->val;

    /* --- Cuenot et al.:
     * STE: Prod := 0
     *      Disp := - (D + Dtur)/(C_k * Delta_les**2)*fp2m
     * VTE: Prod := 2*rho*(D + Dtur)*|grad(Z)|**2
     *      Disp := - (D + Dtur)/(C_k * Delta_les**2)*fp2m
     *
     * --- Pierce:
     * Progress variable equation:
     *      Prod := flamelet_lib(fm, fp2m, ki, progvar)
     */

    /* For the moment, this model for source computation
       is only available in LES */

    if (f_sc == f_fp2m) {

      cs_real_t *fp2m = f_fp2m->val_pre;

      /* Allocate a temporary array for the gradient reconstruction */
      cs_real_3_t *grad;
      CS_MALLOC(grad, n_cells_ext, cs_real_3_t);

      cs_real_t *coefa_p, *coefb_p;
      CS_MALLOC(coefa_p, n_b_faces, cs_real_t);
      CS_MALLOC(coefb_p, n_b_faces, cs_real_t);

      /* Homogeneous Neumann on convective inlet on the
         production term for the variance */

      cs_real_t *coefap = f_fm->bc_coeffs->a;
      cs_real_t *coefbp = f_fm->bc_coeffs->b;

      /* Overwrite diffusion at inlets */
      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
        coefa_p[f_id] = coefap[f_id];
        coefb_p[f_id] = coefbp[f_id];

        if (bc_type[f_id] == CS_CONVECTIVE_INLET) {
          coefap[f_id] = 0.;
          coefbp[f_id] = 1.;
        }

      }

      cs_field_gradient_scalar(f_fm,
                               true, /* use_previous_t */
                               1,    /* inc */
                               grad);

      /* Put back the value */
      for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
        coefap[f_id] = coefa_p[f_id];
        coefbp[f_id] = coefb_p[f_id];
      }

      CS_FREE(coefa_p);
      CS_FREE(coefb_p);

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

        const cs_real_t delta_les
          = cs_turb_xlesfl * pow(cs_turb_ales*cell_vol[c_id], cs_turb_bles);

        const cs_real_t cexp
          =   2.0 * (turb_diff[c_id] + viscls[c_id]) * cell_f_vol[c_id]
            * cs_math_3_dot_product(grad[c_id], grad[c_id])
            - ((turb_diff[c_id] + viscls[c_id])/ (coef_k
            * cs_math_pow2(delta_les)) * fp2m[c_id]) * cell_f_vol[c_id];

        const cs_real_t cimp = 0.;
        smbrs[c_id]  += cexp + cimp * scal_pre[c_id];
        rovsdt[c_id] += cs::max(-cimp, 0.);
      }

      CS_FREE(grad);
    }
    else if (f_sc == cs_field_by_name_try("mixture_fraction_2nd_moment")) {

      cs_real_t *fp2m = f_recvr->val;

      // Not necessary
      //cs_combustion_reconstruct_variance(f_fm->val_pre,
      //                                   f_sc->val_pre,
      //                                   fp2m);

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

        const cs_real_t delta_les
          = cs_turb_xlesfl *pow(cs_turb_ales*cell_vol[c_id], cs_turb_bles);

        const cs_real_t cexp = - (  (turb_diff[c_id] + viscls[c_id])
                                  / (coef_k*cs_math_pow2(delta_les))*fp2m[c_id])
                               * cell_f_vol[c_id];

        const cs_real_t cimp = 0.;
        smbrs[c_id]  += cexp + cimp * scal_pre[c_id];
        rovsdt[c_id] += cs::max(-cimp, 0.);
      }
    }

  } /* End test on model = 41 */

  if (cs_glob_physical_model_flag[CS_COMBUSTION_SLFM] >= 2) {

    if (f_sc == cs_field_by_name_try("progress_variable")) {
      cs_real_t *omega_c = cs_field_by_name("omega_c")->val;

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const cs_real_t cexp = omega_c[c_id];
        const cs_real_t cimp = 0.;

        smbrs[c_id]  += cexp + cimp * scal_pre[c_id];
        rovsdt[c_id] += cs::max(-cimp, 0.);
      }

    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

