/*============================================================================
 * Base thermal model data.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_parall.h"
#include "cs_mesh_location.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_thermal_model.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_thermal_model.c
        base thermal model data.

  \struct cs_thermal_model_t

  \brief Thermal model descriptor.

  Members of this thermal model are publicly accessible, to allow for concise
  syntax, as it is expected to be used in many places.

  \var  cs_thermal_model_t::itherm
        Thermal model
           - 0: no thermal model
           - 1: temperature
           - 2: enthalpy
           - 3: total energy (only for compressible module)\n
        When a particular physics module is activated (gas combustion,
        pulverised coal, electricity or compressible), the user must not
        modify \ref itherm (the choice is made automatically: the solved
        variable is either the enthalpy or the total energy). The user is
        also reminded that, in the case of a coupling with SYRTHES, the
        solved thermal variable should be the temperature (\ref itherm = 1).
        More precisely, everything is designed in the code to allow for the
        running of a calculation coupled with SYRTHES with the enthalpy as
        thermal variable (the correspondence and conversion is then specified
        by the user in the subroutine \ref usthht). However this case has
        never been used in practice and has therefore not been tested. With the
        compressible model, it is possible to carry out calculations coupled
        with SYRTHES, although the thermal scalar represents the total energy
        and not the temperature.
  \var  cs_thermal_model_t::itpscl
        Temperature scale
        - 0: none
        - 1: Kelvin
        - 2: Celsius
        The distinction between \ref itpscl = 1 or 2 is useful only in case of
        radiation modelling. For calculations without radiation modelling,
        use \ref itpscl = 1 for the temperature.\n
        Useful if and only if \ref dimens::nscal "nscal" \f$\geqslant\f$ 1.
  \var  cs_thermal_model_t::iscalt
        Index of the thermal scalar (temperature, energy or enthalpy).\n

        The index of the corresponding variable is isca(iscalt)
        If \ref iscalt = -1, neither the temperature nor the enthalpy is
        represented by a scalar. When a specific physics module is activated
        (gas combustion, pulverised coal, electricity or compressible), the user
        must not modify \ref iscalt (the choice is made automatically). In the
        case of the compressible module, \ref iscalt does not correspond to
        the temperature nor enthalpy but to the total energy}.\n Useful if
        and only if \ref dimens::nscal "nscal" \f$\geqslant\f$ 1.

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* main thermal model structure and associated pointer */

static cs_thermal_model_t  _thermal_model = {
  .itherm = -999,
  .itpscl = 1,
  .iscalt = -1};

const cs_thermal_model_t  *cs_glob_thermal_model = &_thermal_model;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_thermal_model_get_pointers(int     **itherm,
                                int     **itpscl,
                                int     **iscalt);

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------
 * Get pointers to members of the global thermal model structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   itherm --> pointer to cs_glob_thermal_model->itherm
 *   itpscl --> pointer to cs_glob_thermal_model->itpscl
 *   iscalt --> pointer to cs_glob_thermal_model->iscalt
 *----------------------------------------------------------------------------*/

void
cs_f_thermal_model_get_pointers(int     **itherm,
                                int     **itpscl,
                                int     **iscalt)
{
  *itherm = &(_thermal_model.itherm);
  *itpscl = &(_thermal_model.itpscl);
  *iscalt = &(_thermal_model.iscalt);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *!
 * \brief Return thermal field (temperature, enthalpy, total energy according to
 *        thermal model).
 *
 * \return   pointer to thermal field
 *----------------------------------------------------------------------------*/

cs_field_t *
cs_thermal_model_field(void)
{
  cs_field_t *th_f;
  switch (_thermal_model.itherm) {
  case CS_THERMAL_MODEL_TEMPERATURE:
    th_f = CS_F_(t);
    break;
  case CS_THERMAL_MODEL_ENTHALPY:
    th_f = CS_F_(h);
    break;
  case CS_THERMAL_MODEL_TOTAL_ENERGY:
    th_f = CS_F_(energy);
    break;
  default:
    th_f = NULL;
  }

  return th_f;
}

/*----------------------------------------------------------------------------
 *!
 * \brief Provide acces to cs_glob_thermal_model
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_thermal_model_t *
cs_get_glob_thermal_model(void)
{
  return &_thermal_model;
}

/*----------------------------------------------------------------------------
 *!
 * \brief Print the thermal model structure to setup.log.
 *
 *----------------------------------------------------------------------------*/

void
cs_thermal_model_log_setup(void)
{
  cs_log_printf
    (CS_LOG_SETUP,
     _("\n"
       "Thermal model options\n"
       "---------------------\n\n"
       "  Continuous phase:\n\n"
       "    itherm:      %14d (0: no thermal model)\n"
       "                                (1: temperature)\n"
       "                                (2: enthalpy)\n"
       "                                (3: total energy)\n"
       "    itpscl:      %14d (0: none)\n"
       "                                (1: temperature in kelvin)\n"
       "                                (2: temperature in celsius)\n"
       "    iscalt:      %14d (thermal scalar number)\n"),
       cs_glob_thermal_model->itherm,
       cs_glob_thermal_model->itpscl,
       cs_glob_thermal_model->iscalt);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
