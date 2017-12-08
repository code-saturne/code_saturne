/*============================================================================
 * Equation of state
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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
 * EOS library headers
 *----------------------------------------------------------------------------*/

#include <EOS/API/EOS.hxx>
#include <EOS/API/EOS_Field.hxx>
#include <EOS/API/EOS_Fields.hxx>
#include <EOS/API/EOS_Error_Field.hxx>

#if defined(EOS_PRE_V1_6)
  #include <EOS/API/EOS_enums.hxx>
#else
#include <EOS/API/EOS_enums.hxx>
  #include <EOS/API/EOS_properties.hxx>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_eos.hxx"

/*----------------------------------------------------------------------------*/

/*=============================================================================
 *  Global variables
 *============================================================================*/

/* Pointer on a EOS instances array */

NEPTUNE::EOS *eos;

/*============================================================================
 * Private function definitions
 *============================================================================*/
#ifdef __cplusplus
extern "C"
#endif
void
_eos_error_code(const char *function_name,
                NEPTUNE::EOS_Error return_code)
{
    // Generic error code returned by EOS.
    // For a field, the generic error is the worst
    // error of each computed point.
    // When the computation involves several steps (numerical differentiation,
    // function composition, etc...), the returned code is the worst code
    // ever seen in during the computation. When Newton iterations are used
    // the error obtained during the last iteration is returned.

    switch (return_code)
    {
        case NEPTUNE::error:
            bft_error(__FILE__, __LINE__, 0,
                      _("EOS error in call function %s\n"), function_name);
            break;

        case NEPTUNE::bad:
            bft_error(__FILE__, __LINE__, 0,
                      _("EOS bad return in call function %s\n"), function_name);
            break;

        case NEPTUNE::ok:
            bft_printf("EOS not in good agreement in call function %s\n",
                       function_name);
            break;

        case NEPTUNE::good:
            break;

        default:
            break;
    }
    return;
}


/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define new EOS instances
 *
 * parameters:
 *   EOSMethod  <-- table for EOS thermo properties (CATHARE, THETIS, ...)
 *   EOSRef     <-- reference table for EOS thermo properties
 *----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C"
#endif
void
cs_eos_create(char *EOSMethod,
              char *EOSRef)
{
    eos = new NEPTUNE::EOS(EOSMethod, EOSRef);
}

/*----------------------------------------------------------------------------
 * Delete EOS instances
 *----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C"
#endif
void
cs_eos_destroy(void)
{
    delete eos;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computes physical properties in (P,h,Yi) for compressible flow.
 *
 * \param[in]  thermo_plane type of thermal plane
 * \param[in]  property     type of property to compute
 * \param[in]  n_vals       size of properties arrays
 * \param[in]  var1         array of pressure
 * \param[in]  var2         array of thermal properties
 * \param[out] val          array of property
 */
/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C"
#endif
void
cs_phys_prop_eos(cs_phys_prop_thermo_plane_type_t   thermo_plane,
                 cs_phys_prop_type_t                property,
                 const cs_lnum_t                    n_vals,
                 double                             var1[],
                 double                             var2[],
                 cs_real_t                          val[])
{
  NEPTUNE::EOS_Fields input(2);

  //- Two inputs : pressure and thermal variable
  NEPTUNE::EOS_Field eos_pressure("Pressure", "p", n_vals, var1);
  input[0] = eos_pressure;

  if (thermo_plane == CS_PHYS_PROP_PLANE_PH) {
    NEPTUNE::EOS_Field eos_h("Enthalpy", "h", n_vals, var2);
    input[1] = eos_h;
  }
  else if (thermo_plane == CS_PHYS_PROP_PLANE_PT) {
    /* temperature is in K */
    NEPTUNE::EOS_Field eos_T("Temperature", "T", n_vals, var2);
    input[1] = eos_T;
  }
  else {
    bft_error(__FILE__, __LINE__, 0,
              _("bad choice: you chose to work in the %i plane with EOS."), thermo_plane);
  }

  NEPTUNE::EOS_Fields output(1);
  NEPTUNE::EOS_Field *eos_o;

  switch (property) {
    case CS_PHYS_PROP_PRESSURE:
      bft_error(__FILE__, __LINE__, 0,
                _("bad choice: you can't choose pressure as an output variable"));
      break;
    case CS_PHYS_PROP_TEMPERATURE:
      /* temperature is in K */
      eos_o = new NEPTUNE::EOS_Field("Temperature", "T", n_vals, val);
      break;
    case CS_PHYS_PROP_ENTHALPY:
      eos_o = new NEPTUNE::EOS_Field("Enthalpy", "h", n_vals, val);
      break;
    case CS_PHYS_PROP_ENTROPY:
      eos_o = new NEPTUNE::EOS_Field("entropy", "s", n_vals, val);
      break;
    case CS_PHYS_PROP_ISOBARIC_HEAT_CAPACITY:
      eos_o = new NEPTUNE::EOS_Field("SpecificHeatPressure", "cp", n_vals, val);
      break;
    case CS_PHYS_PROP_ISOCHORIC_HEAT_CAPACITY:
      bft_error(__FILE__, __LINE__, 0,
                _("bad choice: isochoric heat capacity not available yet"));
      break;
    case CS_PHYS_PROP_SPECIFIC_VOLUME:
      bft_error(__FILE__, __LINE__, 0,
                _("bad choice: specific volume not available yet"));
      break;
    case CS_PHYS_PROP_DENSITY:
      eos_o = new NEPTUNE::EOS_Field("rho", "rho", n_vals, val);
      break;
    case CS_PHYS_PROP_INTERNAL_ENERGY:
      bft_error(__FILE__, __LINE__, 0,
                _("bad choice: internal energy not available yet"));
      break;
    case CS_PHYS_PROP_QUALITY:
      bft_error(__FILE__, __LINE__, 0,
                _("bad choice: quality not available yet"));
      break;
    case CS_PHYS_PROP_THERMAL_CONDUCTIVITY:
      eos_o = new NEPTUNE::EOS_Field("lambda", "lambda", n_vals, val);
      break;
    case CS_PHYS_PROP_DYNAMIC_VISCOSITY:
      eos_o = new NEPTUNE::EOS_Field("mu", "mu", n_vals, val);
      break;
    case CS_PHYS_PROP_SPEED_OF_SOUND:
      bft_error(__FILE__, __LINE__, 0,
                _("bad choice: speed of sound not available yet"));
      break;
  }

  output[0] = *eos_o;

  int *error;
  BFT_MALLOC(error, n_vals, int);
  NEPTUNE::EOS_Error_Field eos_error_field(n_vals, error);

  NEPTUNE::EOS_Error eos_error = eos->compute(input,
                                              output,
                                              eos_error_field);

  _eos_error_code("cs_phys_prop_eos", eos_error);
  BFT_FREE(error);

}

/*----------------------------------------------------------------------------*/
