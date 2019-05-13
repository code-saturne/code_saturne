/*============================================================================
 * Compute properties for water with Freesteam
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
#if defined(HAVE_EOS)
#include "cs_eos.hxx"
#endif

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

#if defined(HAVE_COOLPROP)
#include "cs_coolprop.hxx"
#endif

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/* Directory name separator
   (historically, '/' for Unix/Linux, '\' for Windows, ':' for Mac
   but '/' should work for all on modern systems) */

#define DIR_SEPARATOR '/'

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Thermal table structure */

typedef struct {

  char        *material;             /* material choice (water, ...) */
  char        *method;               /* method choice
                                        (cathare, thetis, freesteam, ...) */
  char        *reference;            /* reference (automatic) */
  char        *phas;                 /* phase choice (liquid/gas) */
  int          type;                 /* 0 for user
                                      * 1 for freesteam
                                      * 2 for eos
                                      * 3 for coolprop
                                      */

  cs_phys_prop_thermo_plane_type_t   thermo_plane;

  int          temp_scale;           /* temperature scale if needed
                                      * 1 for kelvin
                                      * 2 for Celsius */

} cs_thermal_table_t;

/*----------------------------------------------------------------------------
 * Function pointer types
 *----------------------------------------------------------------------------*/

typedef void
(cs_eos_create_t)(char *EOSMethod,
                  char *EOSRef);

typedef void
(cs_eos_destroy_t)(void);

typedef void
(cs_phys_prop_eos_t)(cs_phys_prop_thermo_plane_type_t   thermo_plane,
                     cs_phys_prop_type_t                property,
                     const cs_lnum_t                    n_vals,
                     double                             var1[],
                     double                             var2[],
                     cs_real_t                          val[]);

typedef void
(cs_phys_prop_coolprop_t)(char                              *CoolPropMaterial,
                          cs_phys_prop_thermo_plane_type_t   thermo_plane,
                          cs_phys_prop_type_t                property,
                          const cs_lnum_t                    n_vals,
                          const cs_real_t                    var1[],
                          const cs_real_t                    var2[],
                          cs_real_t                          val[]);

/*============================================================================
 * Static global variables
 *============================================================================*/

cs_thermal_table_t *cs_glob_thermal_table = NULL;

#if defined(HAVE_DLOPEN) && defined(HAVE_EOS)

static void                     *_cs_eos_dl_lib = NULL;
static cs_eos_create_t          *_cs_eos_create = NULL;
static cs_eos_destroy_t         *_cs_eos_destroy = NULL;
static cs_phys_prop_eos_t       *_cs_phys_prop_eos = NULL;

#endif

#if defined(HAVE_COOLPROP)

static cs_phys_prop_coolprop_t  *_cs_phys_prop_coolprop = NULL;
#if defined(HAVE_PLUGINS)
static void                     *_cs_coolprop_dl_lib = NULL;
#endif

#endif

/*----------------------------------------------------------------------------
 * Create an empty thermal_table structure
 *----------------------------------------------------------------------------*/

static cs_thermal_table_t *
_thermal_table_create(void)
{
  cs_thermal_table_t  *tt = NULL;

  BFT_MALLOC(tt, 1, cs_thermal_table_t);

  tt->material     = NULL;
  tt->method       = NULL;
  tt->reference    = NULL;
  tt->phas         = NULL;
  tt->type         = 0;
  tt->temp_scale   = 0;
  tt->thermo_plane = CS_PHYS_PROP_PLANE_PH;

  return tt;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define thermal table.
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_table_set(const char *material,
                     const char *method,
                     const char *phas,
                     const char *reference,
                     cs_phys_prop_thermo_plane_type_t thermo_plane,
                     const int   temp_scale)
{
  if (cs_glob_thermal_table == NULL)
    cs_glob_thermal_table = _thermal_table_create();

  BFT_MALLOC(cs_glob_thermal_table->material,  strlen(material) +1,  char);
  BFT_MALLOC(cs_glob_thermal_table->reference, strlen(reference) +1, char);
  BFT_MALLOC(cs_glob_thermal_table->phas,      strlen(phas) +1,      char);
  strcpy(cs_glob_thermal_table->material,  material);
  strcpy(cs_glob_thermal_table->reference, reference);
  strcpy(cs_glob_thermal_table->phas,      phas);

  if (strcmp(method, "freesteam") == 0 ||
      strcmp(material, "user_material") == 0) {
    BFT_MALLOC(cs_glob_thermal_table->method,    strlen(method) +1,    char);
    strcpy(cs_glob_thermal_table->reference, reference);
    if (strcmp(method, "freesteam") == 0)
      cs_glob_thermal_table->type = 1;
    else
      cs_glob_thermal_table->type = 0;
  }
  else if (strcmp(method, "CoolProp") == 0) {
    BFT_MALLOC(cs_glob_thermal_table->method,    strlen(method) +1,    char);
    strcpy(cs_glob_thermal_table->reference, reference);
    cs_glob_thermal_table->type = 3;
#if defined(HAVE_COOLPROP)
#if defined(HAVE_PLUGINS)
    {
      /* Open from shared library */
      _cs_coolprop_dl_lib = cs_base_dlopen_plugin("cs_coolprop");

      /* Load symbols from shared library */

      /* Function pointers need to be double-casted so as to first convert
         a (void *) type to a memory address and then convert it back to the
         original type. Otherwise, the compiler may issue a warning.
         This is a valid ISO C construction. */

      _cs_phys_prop_coolprop = (cs_phys_prop_coolprop_t *)  (intptr_t)
        cs_base_get_dl_function_pointer(_cs_coolprop_dl_lib,
                                        "cs_phys_prop_coolprop",
                                        true);
    }
#else
    _cs_phys_prop_coolprop = (cs_phys_prop_coolprop_t *)  (intptr_t)
      cs_phys_prop_coolprop;
#endif
#endif
  }
  else {
    BFT_MALLOC(cs_glob_thermal_table->method,    strlen(method) +5,    char);
    strcpy(cs_glob_thermal_table->method, "EOS_");
    strcat(cs_glob_thermal_table->method, method);
    cs_glob_thermal_table->type = 2;
#if defined(HAVE_EOS)
    {
      /* Open from shared library */
      _cs_eos_dl_lib = cs_base_dlopen_plugin("cs_eos");

      /* Function pointers need to be double-casted so as to first convert
         a (void *) type to a memory address and then convert it back to the
         original type. Otherwise, the compiler may issue a warning.
         This is a valid ISO C construction. */

      _cs_eos_create = (cs_eos_create_t *)  (intptr_t)
        cs_base_get_dl_function_pointer(_cs_eos_dl_lib, "cs_eos_create", true);
      _cs_eos_destroy = (cs_eos_destroy_t *)  (intptr_t)
        cs_base_get_dl_function_pointer(_cs_eos_dl_lib, "cs_eos_destroy", true);
      _cs_phys_prop_eos = (cs_phys_prop_eos_t *)  (intptr_t)
        cs_base_get_dl_function_pointer(_cs_eos_dl_lib, "cs_phys_prop_eos", true);

      _cs_eos_create(cs_glob_thermal_table->method,
                     cs_glob_thermal_table->reference);
    }
#endif
  }
  cs_glob_thermal_table->thermo_plane = thermo_plane;
  cs_glob_thermal_table->temp_scale = temp_scale;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief finalize thermal table.
 */
/*----------------------------------------------------------------------------*/

void
cs_thermal_table_finalize(void)
{
  if (cs_glob_thermal_table != NULL) {
#if defined(HAVE_EOS)
    if (cs_glob_thermal_table->type == 2) {
      _cs_eos_destroy();
      cs_base_dlclose("cs_eos", _cs_eos_dl_lib);
      _cs_eos_create = NULL;
      _cs_eos_destroy = NULL;
      _cs_phys_prop_eos = NULL;
    }
#endif
#if defined(HAVE_COOLPROP) && defined(HAVE_PLUGINS)
    if (cs_glob_thermal_table->type == 3) {
      cs_base_dlclose("cs_coolprop", _cs_coolprop_dl_lib);
      _cs_phys_prop_coolprop = NULL;
    }
#endif
    BFT_FREE(cs_glob_thermal_table->material);
    BFT_FREE(cs_glob_thermal_table->method);
    BFT_FREE(cs_glob_thermal_table->phas);
    BFT_FREE(cs_glob_thermal_table->reference);
    BFT_FREE(cs_glob_thermal_table);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute a physical property.
 *
 * For values var1 and var2, we can use a stride so that accesses for a given
 * element with id i will be of the form: var[i*stride]; this allows regular
 * access with stride 1, and access to constant variables stored as a
 * single-valued array with a stride of 0.
 *
 * \param[in]   property      property queried
 * \param[in]   n_vals        number of values
 * \param[in]   var1_stride   stride between successive values of var1
 * \param[in]   var2_stride   stride between successive values of var2
 * \param[in]   var1          values on first plane axis
 * \param[in]   var2          values on second plane axis
 * \param[out]  val           resulting property values
 */
/*----------------------------------------------------------------------------*/

void
cs_phys_prop_compute(cs_phys_prop_type_t          property,
                     cs_lnum_t                    n_vals,
                     cs_lnum_t                    var1_stride,
                     cs_lnum_t                    var2_stride,
                     const cs_real_t              var1[],
                     const cs_real_t              var2[],
                     cs_real_t                    val[])
{
  cs_lnum_t        _n_vals = n_vals;
  cs_real_t         _var2_c_single[1];
  cs_real_t        *_var1_c = NULL, *_var2_c = NULL;
  const cs_real_t  *var1_c = var1, *var2_c = var2;

  if (n_vals < 1)
    return;

  /* Adapt to different strides to optimize for constant arrays */

  if (var1_stride == 0 && var2_stride == 0)
    _n_vals = 1;

  if (var1_stride == 0 && n_vals > 1) {
    BFT_MALLOC(_var1_c, n_vals, cs_real_t);
    for (cs_lnum_t ii = 0; ii < n_vals; ii++)
      _var1_c[ii] = var1[0];
    var1_c = _var1_c;
  }

  if (cs_glob_thermal_table->temp_scale == 2) {
    if (_n_vals == 1) {
      _var2_c_single[0] = var2[0] + 273.15;
      var2_c = _var2_c_single;
    }
    else {
      BFT_MALLOC(_var2_c, n_vals, cs_real_t);
      for (cs_lnum_t ii = 0; ii < n_vals; ii++)
        _var2_c[ii] = var2[ii*var2_stride] + 273.15;
      var2_c = _var2_c;
    }
  }
  else {
    if (var2_stride == 0 && n_vals > 1) {
      BFT_MALLOC(_var2_c, n_vals, cs_real_t);
      for (cs_lnum_t ii = 0; ii < n_vals; ii++)
        _var2_c[ii] = var2[0];
      var2_c = _var2_c;
    }
  }

  /* Compute proper */

  if (cs_glob_thermal_table->type == 1) {
    cs_phys_prop_freesteam(cs_glob_thermal_table->thermo_plane,
                           property,
                           _n_vals,
                           var1_c,
                           var2_c,
                           val);
  }
#if defined(HAVE_EOS)
  else if (cs_glob_thermal_table->type == 2) {
    _cs_phys_prop_eos(cs_glob_thermal_table->thermo_plane,
                      property,
                      _n_vals,
                      var1_c,
                      var2_c,
                      val);
  }
#endif
#if defined(HAVE_COOLPROP)
  else if (cs_glob_thermal_table->type == 3) {
    _cs_phys_prop_coolprop(cs_glob_thermal_table->material,
                           cs_glob_thermal_table->thermo_plane,
                           property,
                           _n_vals,
                           var1_c,
                           var2_c,
                           val);
  }
#endif
  BFT_FREE(_var1_c);
  BFT_FREE(_var2_c);

  /* In case of single value, apply to all */

  if (_n_vals == 1) {
    cs_real_t val_const = val[0];
    for (cs_lnum_t ii = 0; ii < n_vals; ii++)
      val[ii] = val_const;
  }
}

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
                  _("bad choice: you choose to work in the %s plane."), "ph");
        break;
      case CS_PHYS_PROP_TEMPERATURE:
        val[i] = freesteam_T(S0);
        break;
      case CS_PHYS_PROP_ENTHALPY:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice: you choose to work in the %s plane."), "ph");
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
                  _("bad choice: you choose to work in the %s plane."), "pT");
        break;
      case CS_PHYS_PROP_TEMPERATURE:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice: you choose to work in the %s plane."), "pT");
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
                  _("bad choice: you choose to work in the %s plane."), "ps");
        break;
      case CS_PHYS_PROP_TEMPERATURE:
        val[i] = freesteam_T(S0);
        break;
      case CS_PHYS_PROP_ENTHALPY:
        val[i] = freesteam_h(S0);
        break;
      case CS_PHYS_PROP_ENTROPY:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice: you choose to work in the %s plane."), "ps");
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
                  _("bad choice: you choose to work in the %s plane."), "pu");
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
                  _("bad choice: you choose to work in the %s plane."), "pu");
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
                  _("bad choice: you choose to work in the %s plane."), "pv");
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
                  _("bad choice: you choose to work in the %s plane."), "pv");
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
                  _("bad choice: you choose to work in the %s plane."), "Ts");
        break;
      case CS_PHYS_PROP_ENTHALPY:
        val[i] = freesteam_h(S0);
        break;
      case CS_PHYS_PROP_ENTROPY:
        bft_error(__FILE__, __LINE__, 0,
                  _("bad choice: you choose to work in the %s plane."), "Ts");
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
                  _("bad choice: you choose to work in the %s plane."), "Tx");
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
                  _("bad choice: you choose to work in the %s plane."), "Tx");
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
  CS_UNUSED(thermo_plane);
  CS_UNUSED(property);
  CS_UNUSED(n_vals);
  CS_UNUSED(var1);
  CS_UNUSED(var2);
  CS_UNUSED(val);

  bft_error(__FILE__, __LINE__, 0,
            _("Freesteam support not available in this build."));
#endif
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
