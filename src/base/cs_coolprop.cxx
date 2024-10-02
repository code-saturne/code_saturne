/*============================================================================
 * Equation of state
 *============================================================================*/

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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * CoolProp library headers
 *----------------------------------------------------------------------------*/

extern "C" {
#include "CoolPropLib.h"
}

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"
#include "cs_fp_exception.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_coolprop.hxx"

/*----------------------------------------------------------------------------*/

/* Thermal table structure */

typedef struct {

  long         abstract_state_handle;  /* handle to abstract state */

  char         backend[32];            /* backend name */
  char         fluid_names[32];        /* fluid names */

} cs_coolprop_state_t;

/*=============================================================================
 *  Global variables
 *============================================================================*/

static int  _n_states = 0;
static cs_coolprop_state_t *_states = nullptr;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computes physical properties in (P,h,Yi) for compressible flow.
 *
 * \param[in]   coolprop_material  CoolProp material
 * \param[in]   coolprop_backend   CoolProp backend ("HEOS" by default,
 *                                 "SRK" for cubic, "TTSE&HEOS" or
 *                                 "BICUBIC&HEOS" for tabulated)
 * \param[in]   thermo_plane       type of thermal plane
 * \param[in]   property           type of property to compute
 * \param(in]   n_vals             size of variable and property arrays
 * \param[in]   var1               first variable of thermodynamic plane
 *                                 (pressure)
 * \param[in]   var2               second variable of thermodynamic plane
 * \param[out]  val                computed property values
 */
/*----------------------------------------------------------------------------*/

void
cs_phys_prop_coolprop(char                              *coolprop_material,
                      const char                        *coolprop_backend,
                      cs_phys_prop_thermo_plane_type_t   thermo_plane,
                      cs_phys_prop_type_t                property,
                      const cs_lnum_t                    n_vals,
                      const cs_real_t                    var1[],
                      const cs_real_t                    var2[],
                      cs_real_t                          val[])
{
  const int mb_size = 511;
  char message_buffer[mb_size + 1];

  const char *fluid_names = coolprop_material;

  char name1[] = "P";

  char backend[32];
  if (coolprop_backend != nullptr)
    strncpy(backend, coolprop_backend, 32);
  else
    strncpy(backend, "HEOS", 32);
  backend[31] = '\0';

  long errcode = 0;
  long state_handle = -1;

  /* API choice for easier testing:
     - 0: high level interface.
     - 1: low level interfaces with array data.
     - 2: low level interface with individual calls per element.
  */
  int api_choice = 0;

  /* The high level interface is not compatible with tabulated backends */
  if (api_choice == 0) {
    if (   strncmp(backend, "TTSE", 4) == 0
        || strncmp(backend, "BICUBIC", 7) == 0)
      api_choice = 1;
  }

  if (api_choice > 0) {

    for (int i = 0; i < _n_states; i++) {
      if (   strncmp(backend, _states[i].backend, 31) == 0
             && strncmp(fluid_names, _states[i].fluid_names, 31) == 0) {
        state_handle = _states[i].abstract_state_handle;
        break;
      }
    }

    if (state_handle < 0) {
      state_handle = AbstractState_factory(backend, fluid_names,
                                           &errcode, message_buffer, mb_size);
      if (errcode != 0)
        bft_error(__FILE__, __LINE__, 0, "%s", message_buffer);

      double fractions[] = {1.0};
      AbstractState_set_fractions(state_handle, fractions, 1,
                                  &errcode, message_buffer, mb_size);
      if (errcode != 0)
        bft_error(__FILE__, __LINE__, 0, "%s", message_buffer);

      AbstractState_specify_phase(state_handle, "phase_unknown",
                                  &errcode, message_buffer, mb_size);
      if (errcode != 0)
        bft_error(__FILE__, __LINE__, 0, "%s", message_buffer);

      BFT_REALLOC(_states, _n_states+1, cs_coolprop_state_t);

      _states[_n_states].abstract_state_handle = state_handle;
      memset(_states[_n_states].backend, 0, 32);
      strncpy(_states[_n_states].backend, backend, 31);
      memset(_states[_n_states].fluid_names, 0, 32);
      strncpy(_states[_n_states].fluid_names, fluid_names, 31);
      _n_states += 1;
    }
  }

  char name2[2];
  long input_pair = -1;
  if (thermo_plane == CS_PHYS_PROP_PLANE_PH) {
    input_pair = get_input_pair_index("PH_INPUTS");
    strncpy(name2, "H", 2);
  }
  else if (thermo_plane == CS_PHYS_PROP_PLANE_PT) {
    input_pair = get_input_pair_index("PT_INPUTS");
    strncpy(name2, "T", 2);
  }
  name2[1] = '\0';

  char outputs[32];
  switch (property) {
    case CS_PHYS_PROP_PRESSURE:
      bft_error(__FILE__, __LINE__, 0,
                _("bad choice: you can't choose pressure as an output variable"));
      break;
    case CS_PHYS_PROP_TEMPERATURE:
      /* temperature is in K */
      strncpy(outputs, "T", 31);
      break;
    case CS_PHYS_PROP_ENTHALPY:
      strncpy(outputs, "H", 31);
      break;
    case CS_PHYS_PROP_ENTROPY:
      strncpy(outputs, "S", 31);
      break;
    case CS_PHYS_PROP_ISOBARIC_HEAT_CAPACITY:
      strncpy(outputs, "CPMASS", 31);
      break;
    case CS_PHYS_PROP_ISOCHORIC_HEAT_CAPACITY:
      strncpy(outputs, "CVMASS", 31);
      break;
    case CS_PHYS_PROP_SPECIFIC_VOLUME:
      bft_error(__FILE__, __LINE__, 0,
                _("bad choice: specific volume not available yet"));
      break;
    case CS_PHYS_PROP_DENSITY:
      strncpy(outputs, "D", 31);
      break;
    case CS_PHYS_PROP_INTERNAL_ENERGY:
      strncpy(outputs, "U", 31);
      break;
    case CS_PHYS_PROP_QUALITY:
      bft_error(__FILE__, __LINE__, 0,
                _("bad choice: quality not available yet"));
      break;
    case CS_PHYS_PROP_THERMAL_CONDUCTIVITY:
      strncpy(outputs, "L", 31);
      break;
    case CS_PHYS_PROP_DYNAMIC_VISCOSITY:
      strncpy(outputs, "V", 31);
      break;
    case CS_PHYS_PROP_SPEED_OF_SOUND:
      strncpy(outputs, "SPEED_OF_SOUND", 31);
      break;
  }
  outputs[31] = '\0';
  long output_key = get_param_index(outputs);

  double *prop1, *prop2, *result;
  double *_prop1 = nullptr, *_prop2 = nullptr, *_result = nullptr;

  if (sizeof(cs_real_t) != sizeof(double)) {
    BFT_MALLOC(_prop1, n_vals, double);
    BFT_MALLOC(_prop2, n_vals, double);
    BFT_MALLOC(_result, n_vals, double);
    for (cs_lnum_t i = 0; i < n_vals; i++) {
      _prop1[i] = var1[i];
      _prop2[i] = var1[i];
    }
    prop1 = _prop1;
    prop2 = _prop2;
    result = _result;
  }
  else {
    prop1 = const_cast<double *>(var1);
    prop2 = const_cast<double *>(var2);
    result = val;
  }

  cs_fp_exception_disable_trap();

  switch(api_choice) {
  case 0:
    {
      long resdim1 = n_vals, resdim2 = 1;
      double fractions[1] = {1.0};

      PropsSImulti(outputs,
                   name1, prop1, n_vals,
                   name2, prop2, n_vals,
                   backend,
                   fluid_names,
                   fractions, 1,
                   result, &resdim1, &resdim2);

      if (resdim1 != n_vals) {
        get_global_param_string("errstring", message_buffer, mb_size);
        bft_error(__FILE__, __LINE__, 0,
                  _("CoolProp was unable to compute some fluid properties with\n"
                    "input variable names: \"%s\", \"%s\" and backend \"%s\".\n\n"
                    "Fluid(s) considered: \"%s\".\n\n"
                    "%s."),
                  name1, name2, backend, fluid_names, message_buffer);
      }
      assert(resdim2 == 1);
    }
    break;

  case 1:
    {
      AbstractState_update_and_1_out(state_handle, input_pair,
                                     prop1, prop2, n_vals, output_key,
                                     result,
                                     &errcode, message_buffer, mb_size);
      if (errcode != 0)
        bft_error(__FILE__, __LINE__, 0,
                  _("CoolProp was unable to compute property \"%s\" with\n"
                    "input variable names: \"%s\", \"%s\" and backend \"%s\".\n\n"
                    "Fluid(s) considered: \"%s\".\n\n"
                    "%s."),
                  outputs, name1, name2, backend, fluid_names, message_buffer);
    }
    break;

  case 2:
    {
      for (int i = 0; i < n_vals; i++) {
        AbstractState_update(state_handle, input_pair,
                             prop1[i], prop2[i],
                             &errcode, message_buffer, mb_size);
        if (errcode == 0)
          result[i] = AbstractState_keyed_output(state_handle, output_key,
                                                 &errcode, message_buffer, mb_size);
        if (errcode != 0)
          bft_error(__FILE__, __LINE__, 0,
                    _("CoolProp was unable to compute property \"%s\" with\n"
                      "input variable names: \"%s\", \"%s\" and backend \"%s\".\n\n"
                      "Fluid(s) considered: \"%s\".\n\n"
                      "%s."),
                    outputs, name1, name2, backend, fluid_names, message_buffer);
      }
    }
    break;

  default:
    assert(0);
  }

  cs_fp_exception_restore_trap();

  if (_result != nullptr) {
    for (int i = 0; i < n_vals; i++)
      val[i] = result[i];
    BFT_FREE(_result);
  }

  BFT_FREE(_prop1);
  BFT_FREE(_prop2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free state for Coolprop plugin.
 */
/*----------------------------------------------------------------------------*/

void
cs_coolprop_finalize(void)
{
  for (int i = 0; i < _n_states; i++) {
    long errcode = 0;
    char message_buffer[512];

    AbstractState_free(_states[i].abstract_state_handle,
                       &errcode, message_buffer, 511);
  }

  BFT_FREE(_states);
  _n_states = 0;
}

END_C_DECLS

/*----------------------------------------------------------------------------*/



