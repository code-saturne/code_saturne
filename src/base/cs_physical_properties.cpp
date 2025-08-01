/*============================================================================
 * Physical properties computation and management.
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
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"
#include "base/cs_log.h"
#include "base/cs_mem.h"
#include "base/cs_timer.h"

#include "cdo/cs_property.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_physical_properties.h"
#if defined(HAVE_EOS)
#include "cs_eos.hxx"
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

/* Data elements in crystal router */

typedef enum {

  CS_PHYS_PROP_TABLE_USER,
  CS_PHYS_PROP_TABLE_EOS,
  CS_PHYS_PROP_TABLE_COOLPROP

} cs_phys_prop_table_type_t;

/* Thermal table structure */

typedef struct {

  char        *material;             /* material choice (water, ...) */
  char        *method;               /* method choice
                                        (cathare, thetis, CoolProp, ...) */
  cs_phys_prop_table_type_t  type;   /* 0 for user
                                      * 2 for EOS
                                      * 3 for CoolProp  */

  cs_phys_prop_thermo_plane_type_t   thermo_plane;

  cs_temperature_scale_t             temp_scale;  /* temperature scale  */

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
                     const double                       var1[],
                     const double                       var2[],
                     cs_real_t                          val[]);

typedef void
(cs_phys_prop_coolprop_t)(const char                        *coolprop_material,
                          const char                        *coolprop_backend,
                          cs_phys_prop_thermo_plane_type_t   thermo_plane,
                          cs_phys_prop_type_t                property,
                          const cs_lnum_t                    n_vals,
                          const cs_real_t                    var1[],
                          const cs_real_t                    var2[],
                          cs_real_t                          val[]);

typedef void
(cs_finalize_t)(void);

/*============================================================================
 * Static global variables
 *============================================================================*/

cs_thermal_table_t *cs_glob_thermal_table = nullptr;

static cs_timer_counter_t   _physprop_lib_t_tot;   /* Total time in physical
                                                      property library calls */

#if defined(HAVE_DLOPEN) && defined(HAVE_EOS)

static void                     *_cs_eos_dl_lib = nullptr;
static cs_eos_create_t          *_cs_eos_create = nullptr;
static cs_eos_destroy_t         *_cs_eos_destroy = nullptr;
static cs_phys_prop_eos_t       *_cs_phys_prop_eos = nullptr;

#endif

#if defined(HAVE_COOLPROP)

static char                     *_cs_coolprop_backend = nullptr;
static cs_phys_prop_coolprop_t  *_cs_phys_prop_coolprop = nullptr;
#if defined(HAVE_PLUGINS)
static void                     *_cs_coolprop_dl_lib = nullptr;
#endif

#endif

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create an empty thermal_table structure
 *----------------------------------------------------------------------------*/

static cs_thermal_table_t *
_thermal_table_create(void)
{
  cs_thermal_table_t  *tt = nullptr;

  CS_MALLOC(tt, 1, cs_thermal_table_t);

  tt->material     = nullptr;
  tt->method       = nullptr;
  tt->type         = CS_PHYS_PROP_TABLE_USER;
  tt->temp_scale   = CS_TEMPERATURE_SCALE_NONE;
  tt->thermo_plane = CS_PHYS_PROP_PLANE_PH;

  return tt;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get xdef of a property on a given zone.
 *
 * \param[in] pty    pointer to cs_property_t
 * \param[in] zname  name of the zone. Can be null for 'all_cells'
 *
 * \return pointer to corresponding cs_xdef_t
 *
 */
/*----------------------------------------------------------------------------*/

static cs_xdef_t *
_get_property_def_on_zone(const cs_property_t   *pty,
                          const char            *zname)
{
  cs_xdef_t *def = nullptr;

  const int z_id = cs_volume_zone_id_by_name(zname);
  for (int i = 0; i < pty->n_definitions; i++) {
    if (pty->defs[i]->z_id == z_id) {
      def = pty->defs[i];
      break;
    }
  }

  return def;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief update values of a definition
 *
 * \param[in] def       pointer to cs_xdef_t structure
 * \param[in] new_vals  array of new values
 *
 */
/*----------------------------------------------------------------------------*/

static void
_update_def_values(cs_xdef_t         *def,
                   const cs_real_t   *new_vals)
{
  cs_real_t *_context = (cs_real_t *)def->context;

  for (int i = 0; i < def->dim; i++)
    _context[i] = new_vals[i];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a new physical property using the cs_property_t structure.
 *
 * \param[in] name   name of the property
 * \param[in] dim    dimension of the property (1->scalar, 3->vector,..)
 * \param[in] refval Reference value
 *
 * \return pointer to the newly created cs_property_t structure.
 *
 */
/*----------------------------------------------------------------------------*/

static cs_property_t *
_physical_property_create(const char      *name,
                          const int        dim,
                          const cs_real_t  refval)
{
  cs_property_t *pty = cs_property_by_name(name);
  if (pty != nullptr)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: property '%s' is already defined.\n"),
                name);

  if (dim == 1)
    pty = cs_property_add(name, CS_PROPERTY_ISO);
  else if (dim == 3)
    pty = cs_property_add(name, CS_PROPERTY_ORTHO);
  else if (dim == 6)
    pty = cs_property_add(name, CS_PROPERTY_ANISO_SYM);
  else if (dim == 9)
    pty = cs_property_add(name, CS_PROPERTY_ANISO);
  else
    bft_error(__FILE__, __LINE__, 0,
              _("Error: for property '%s', dimension %d not supported.\n"),
              name, dim);

  cs_property_set_reference_value(pty, refval);

  return pty;
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
cs_thermal_table_set(const char                        *material,
                     const char                        *method,
                     const char                        *reference,
                     cs_phys_prop_thermo_plane_type_t   thermo_plane,
                     cs_temperature_scale_t             temp_scale)
{
  if (cs_glob_thermal_table == nullptr) {
    cs_glob_thermal_table = _thermal_table_create();
    CS_TIMER_COUNTER_INIT(_physprop_lib_t_tot);
  }

  CS_MALLOC(cs_glob_thermal_table->material, strlen(material) +1, char);
  strcpy(cs_glob_thermal_table->material,  material);

  if (strcmp(material, "user_material") == 0) {
    CS_MALLOC(cs_glob_thermal_table->method, strlen(method) +1, char);
    strcpy(cs_glob_thermal_table->method, method);
    cs_glob_thermal_table->type = CS_PHYS_PROP_TABLE_USER;
  }
  else if (strcmp(method, "CoolProp") == 0) {
    cs_glob_thermal_table->type = CS_PHYS_PROP_TABLE_COOLPROP;
#if defined(HAVE_COOLPROP)
    cs_physical_properties_set_coolprop_backend(_cs_coolprop_backend);
#if defined(HAVE_PLUGINS)
    {
      /* Open from shared library */
      _cs_coolprop_dl_lib = cs_base_dlopen_plugin("cs_coolprop");

      /* Load symbols from shared library */

      /* Function pointers need to be double-cast so as to first convert
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
    CS_MALLOC(cs_glob_thermal_table->method, strlen(method) +5, char);
    strcpy(cs_glob_thermal_table->method, "EOS_");
    strcat(cs_glob_thermal_table->method, method);
    cs_glob_thermal_table->type = CS_PHYS_PROP_TABLE_EOS;
#if defined(HAVE_EOS)
#if defined(HAVE_DLOPEN)
    {
      const char _reference_default[] = "";
      const char *_reference = (reference != nullptr) ?
                                 reference : _reference_default;

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
                     const_cast<char *>(_reference));
    }
#else
    bft_error(__FILE__, __LINE__, 0,
              "EOS only available if when dlloader is enabled in code_saturne.");
#endif
#else
    CS_UNUSED(reference);
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
  if (cs_glob_thermal_table != nullptr) {

    if (cs_glob_thermal_table->method != nullptr)
      cs_log_printf(CS_LOG_PERFORMANCE,
                    _("\n"
                      "Physical property computations:\n"
                      "  Total elapsed time for  %s:  %.3f s\n"),
                    cs_glob_thermal_table->method,
                    _physprop_lib_t_tot.nsec*1e-9);

#if defined(HAVE_EOS) /* always a plugin */
    if (cs_glob_thermal_table->type == CS_PHYS_PROP_TABLE_EOS) {
      _cs_eos_destroy();
      cs_base_dlclose("cs_eos", _cs_eos_dl_lib);
      _cs_eos_create = nullptr;
      _cs_eos_destroy = nullptr;
      _cs_phys_prop_eos = nullptr;
    }
#endif
#if defined(HAVE_COOLPROP) && defined(HAVE_PLUGINS)
    if (cs_glob_thermal_table->type == CS_PHYS_PROP_TABLE_COOLPROP) {
      cs_finalize_t *coolprop_finalize
        = (cs_finalize_t *)  (intptr_t)
             cs_base_get_dl_function_pointer(_cs_coolprop_dl_lib,
                                             "cs_coolprop_finalize",
                                             true);
      coolprop_finalize();
      cs_base_dlclose("cs_coolprop", _cs_coolprop_dl_lib);
      _cs_phys_prop_coolprop = nullptr;
    }
    CS_FREE(_cs_coolprop_backend);
#endif
    CS_FREE(cs_glob_thermal_table->material);
    CS_FREE(cs_glob_thermal_table->method);
    CS_FREE(cs_glob_thermal_table);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get current thermal table plane
 */
/*----------------------------------------------------------------------------*/

cs_phys_prop_thermo_plane_type_t
cs_thermal_table_get_thermo_plane(void)
{
  if (cs_glob_thermal_table != nullptr) {
    return cs_glob_thermal_table->thermo_plane;
  }
  else {
    return CS_PHYS_PROP_PLANE_PH; /* default */
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get backend set for CoolProp.
 *
 * Returns null if CoolProp is not used or backend not set yet.
 *
 * \return  pointer to CoolProp backend.
 */
/*----------------------------------------------------------------------------*/

const char *
cs_physical_properties_get_coolprop_backend(void)
{
#if defined(HAVE_COOLPROP)
  return _cs_coolprop_backend;
#else
  return nullptr;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set backend for CoolProp.
 *
 * Ignored if CoolProp is not used.
 *
 * When called from user-defined functions, this should be set from
 * cs_user_model rather than cs_user_parameters, as some reference property
 * values may be computed before calling cs_user_parameters.
 *
 * A few primary backends in CoolProp are:
 *
 * - "HEOS": The Helmholtz Equation of State backend for use with pure and
 *   pseudo-pure fluids, and mixtures, all of which are based on multi-parameter
 *   Helmholtz Energy equations of state.
 * - "REFPROP": only if REFPROP library is available
 *   (set ALTERNATIVE_REFPROP_PATH environment variable if needed)
 * - "INCOMP": Incompressible backend (for pure fluids)
 * - "TTSE&XXXX": TTSE backend, with tables generated using the XXXX backend
 *   where XXXX is one of the base backends("HEOS", "REFPROP", etc.)
 * - "BICUBIC&XXXX": Bicubic backend, with tables generated using the XXXX
 *   backend where XXXX is one of the base backends("HEOS", "REFPROP", etc.)
 *
 * \param[in]  backend  backend name; "HEOS" used if nullptr.
 */
/*----------------------------------------------------------------------------*/

void
cs_physical_properties_set_coolprop_backend(const char  *backend)
{
#if defined(HAVE_COOLPROP)

  if (backend == nullptr)
    backend = "HEOS";

  size_t l = strlen(backend);
  CS_REALLOC(_cs_coolprop_backend, l+1, char);
  strcpy(_cs_coolprop_backend, backend);

  if (cs_glob_thermal_table != nullptr) {

    CS_REALLOC(cs_glob_thermal_table->method,
               strlen("CoolProp ()") + strlen(_cs_coolprop_backend) + 3,
               char);
    strcpy(cs_glob_thermal_table->method, "CoolProp (");
    strcat(cs_glob_thermal_table->method, _cs_coolprop_backend);
    strcat(cs_glob_thermal_table->method, ")");

  }

#else
  CS_UNUSED(backend);
#endif
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
  cs_real_t        *_var1_c = nullptr, *_var2_c = nullptr;
  const cs_real_t  *var1_c = var1, *var2_c = var2;

  if (n_vals < 1)
    return;

  /* Adapt to different strides to optimize for constant arrays */

  if (var1_stride == 0 && var2_stride == 0)
    _n_vals = 1;

  if (var1_stride == 0 && n_vals > 1) {
    CS_MALLOC(_var1_c, n_vals, cs_real_t);
    for (cs_lnum_t ii = 0; ii < n_vals; ii++)
      _var1_c[ii] = var1[0];
    var1_c = _var1_c;
  }

  if (   cs_glob_thermal_table->temp_scale == CS_TEMPERATURE_SCALE_CELSIUS
      && (cs_glob_thermal_table->thermo_plane == CS_PHYS_PROP_PLANE_PT
          || cs_glob_thermal_table->thermo_plane == CS_PHYS_PROP_PLANE_TS
          || cs_glob_thermal_table->thermo_plane == CS_PHYS_PROP_PLANE_TX)) {
    if (_n_vals == 1) {
      _var2_c_single[0] = var2[0] + 273.15;
      var2_c = _var2_c_single;
    }
    else {
      CS_MALLOC(_var2_c, n_vals, cs_real_t);
      for (cs_lnum_t ii = 0; ii < n_vals; ii++)
        _var2_c[ii] = var2[ii*var2_stride] + 273.15;
      var2_c = _var2_c;
    }
  }
  else {
    if (var2_stride == 0 && n_vals > 1) {
      CS_MALLOC(_var2_c, n_vals, cs_real_t);
      for (cs_lnum_t ii = 0; ii < n_vals; ii++)
        _var2_c[ii] = var2[0];
      var2_c = _var2_c;
    }
  }

  cs_timer_t t0 = cs_timer_time();

  /* Compute proper */

#if defined(HAVE_EOS) /* always a plugin */
  if (cs_glob_thermal_table->type == CS_PHYS_PROP_TABLE_EOS) {
    _cs_phys_prop_eos(cs_glob_thermal_table->thermo_plane,
                      property,
                      _n_vals,
                      var1_c,
                      var2_c,
                      val);
  }
#endif
#if defined(HAVE_COOLPROP)
  if (cs_glob_thermal_table->type == CS_PHYS_PROP_TABLE_COOLPROP) {
    _cs_phys_prop_coolprop(cs_glob_thermal_table->material,
                           _cs_coolprop_backend,
                           cs_glob_thermal_table->thermo_plane,
                           property,
                           _n_vals,
                           var1_c,
                           var2_c,
                           val);
  }
#endif

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&_physprop_lib_t_tot, &t0, &t1);

  CS_FREE(_var1_c);
  CS_FREE(_var2_c);

  /* In case of single value, apply to all */

  if (_n_vals == 1) {
    cs_real_t val_const = val[0];
    for (cs_lnum_t ii = 0; ii < n_vals; ii++)
      val[ii] = val_const;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get reference value of a physical property
 *
 * \param[in] name  property name
 *
 * \return reference value (cs_real_t)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_physical_property_get_ref_value(const char  *name)
{

  const cs_property_t *pty = cs_property_by_name(name);
  if (pty == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: property '%s' does not exist\n"),
              name);

  return pty->ref_value;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set reference value for a physical property
 *
 * \param[in] name  property name
 * \param[in] val   new value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_physical_property_set_ref_value(const char      *name,
                                   const cs_real_t  val)
{
  cs_property_t *pty = cs_property_by_name(name);
  if (pty == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: property '%s' does not exist\n"),
              name);

  cs_property_set_reference_value(pty, val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get property reference values for a given zone
 *
 * \param[in] name    property name
 * \param[in] zname   zone name
 * \param[in] retval  array of values to return
 */
/*----------------------------------------------------------------------------*/

void
cs_physical_property_get_zone_values(const char  *name,
                                     const char  *zname,
                                     cs_real_t    retval[])
{
  const cs_property_t *pty = cs_property_by_name(name);
  if (pty == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: property '%s' does not exist\n"),
              name);

  cs_xdef_t *def = _get_property_def_on_zone(pty, zname);
  if (def == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              _("Error: property '%s' does not have a definition for "
                "zone '%s'\n"),
              name, zname);

  /* Sanity check */
  assert(def->type == CS_XDEF_BY_VALUE);

  if (pty->type & CS_PROPERTY_ISO) {
    const cs_real_t *_context = (cs_real_t *)def->context;
    retval[0] = _context[0];

  } else if (pty->type & CS_PROPERTY_ORTHO) {
    const cs_real_t *_context = (cs_real_t *)def->context;
    for (int j = 0; j < 3; j++)
      retval[j] = _context[j];

  } else if (pty->type & CS_PROPERTY_ANISO_SYM) {
    const cs_real_t *_context = (cs_real_t *)def->context;
    for (int j = 0; j < 6; j++)
      retval[j] = _context[j];

  } else if (pty->type & CS_PROPERTY_ANISO) {
    const cs_real_3_t *_context = (cs_real_3_t *)def->context;
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        retval[3*j + k] = _context[j][k];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a physical property
 *
 * \param[in] name    property name
 * \param[in] dim     property dimension
 * \param[in] refval  reference value
 */
/*----------------------------------------------------------------------------*/

void
cs_physical_property_create(const char      *name,
                            const int        dim,
                            const cs_real_t  refval)
{
  _physical_property_create(name, dim, refval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a property definition on a given zone using a single value
 *
 * \param[in] name   property name
 * \param[in] zname  zone name
 * \param[in] dim    property dimension
 * \param[in] val    reference value for the zone
 */
/*----------------------------------------------------------------------------*/

void
cs_physical_property_define_from_value(const char       *name,
                                       const char       *zname,
                                       const int         dim,
                                       const cs_real_t   val)
{
  cs_property_t *pty = cs_property_by_name(name);
  if (pty == nullptr)
    pty = _physical_property_create(name, dim, 0.);

  if (dim == 1) {
    cs_property_def_iso_by_value(pty, zname, val);
  }
  else if (dim == 3) {
    cs_real_t dvals[3] = {val, val, val};
    cs_property_def_ortho_by_value(pty, zname, dvals);
  } else if (dim == 6) {
    cs_real_t dvals[6] = {val, val, val, val, val, val};
    cs_property_def_aniso_sym_by_value(pty, zname, dvals);
  }
  else if (dim == 9) {
    cs_real_t dvals[3][3] = { {val, 0., 0.},
                              {0., val, 0.},
                              {0., 0., val} };
    cs_property_def_aniso_by_value(pty, zname, dvals);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a property multi-diemnsional definition on a given zone
 *
 * \param[in] name   property name
 * \param[in] zname  zone name
 * \param[in] dim    property dimension (>1)
 * \param[in] vals   array of values to set
 */
/*----------------------------------------------------------------------------*/

void
cs_physical_property_define_from_values(const char  *name,
                                        const char  *zname,
                                        const int    dim,
                                        cs_real_t    vals[])
{
  assert(dim > 1 && vals != nullptr);

  cs_property_t *pty = cs_property_by_name(name);

  if (pty == nullptr)
    pty = _physical_property_create(name, dim, 0.);

  if (dim == 3)
    cs_property_def_ortho_by_value(pty, zname, vals);
  else if (dim == 6)
    cs_property_def_aniso_sym_by_value(pty, zname, vals);
  else if (dim == 9) {
    cs_real_3_t *vals2use = (cs_real_3_t *)vals;
    cs_property_def_aniso_by_value(pty, zname, vals2use);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a property definition based on a cs_field_t.
 *
 * The field is created if needed
 *
 * \param[in] name          property name
 * \param[in] type_flag     field type flag
 * \param[in] location_id   location id flag
 * \param[in] dim           field dimension
 * \param[in] has_previous  does the field has val_pre
 */
/*----------------------------------------------------------------------------*/

void
cs_physical_property_define_from_field(const char  *name,
                                       int          type_flag,
                                       int          location_id,
                                       int          dim,
                                       bool         has_previous)
{
  cs_property_t *pty = cs_property_by_name(name);
  if (pty == nullptr)
    pty = _physical_property_create(name, dim, 0.);

  cs_field_t *f = cs_field_by_name_try(name);
  if (f == nullptr)
    f = cs_field_create(name, type_flag, location_id, dim, has_previous);

  cs_property_def_by_field(pty, f);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return id of field associated to property
 *
 * \param[in] name  property name
 *
 * \return field id (int)
 */
/*----------------------------------------------------------------------------*/

int
cs_physical_property_field_id_by_name(const char  *name)
{
  int retval = -1;

  cs_field_t *f = cs_field_by_name_try(name);

  if (f != nullptr)
    retval = f->id;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update reference values for a property on a given zone
 *
 * \param[in] name   property name
 * \param[in] zname  zone name
 * \param[in] vals   array of values to set
 */
/*----------------------------------------------------------------------------*/

void
cs_physical_property_update_zone_values(const char       *name,
                                        const char       *zname,
                                        const cs_real_t   vals[])
{
  cs_property_t *pty = cs_property_by_name(name);

  cs_xdef_t *def = _get_property_def_on_zone(pty, zname);

  _update_def_values(def, vals);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
