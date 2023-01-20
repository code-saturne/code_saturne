/*============================================================================
 * Base thermal model data.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

  \var  cs_thermal_model_t::thermal_variable
        Thermal variable solved for this physical model.

        When a particular physics module is activated (gas combustion,
        pulverized coal, electricity or compressible), the user must not
        modify \ref thermal_variable (the choice is made automatically:
        the solved variable is either the enthalpy or the total energy).

  \var  cs_thermal_model_t::itherm
        \deprecated alias/old name for thermal_variable

  \var  cs_thermal_model_t::temperature_scale
        Temperature scale
        The specification of the temperature scale in a consistent
        manner with the values used (initial and boundary conditions)
        is especially important in case of radiation modelling.

  \var  cs_thermal_model_t::itpscl
        \deprecated alias/old name for temperature_scale

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
  .itpscl = 1};

const cs_thermal_model_t  *cs_glob_thermal_model = &_thermal_model;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_thermal_model_get_pointers(int     **itherm,
                                int     **itpscl);

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
 *----------------------------------------------------------------------------*/

void
cs_f_thermal_model_get_pointers(int     **itherm,
                                int     **itpscl)
{
  *itherm = &(_thermal_model.itherm);
  *itpscl = &(_thermal_model.itpscl);
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
    th_f = CS_F_(e_tot);
    break;
  default:
    th_f = NULL;
  }

  return th_f;
}

/*----------------------------------------------------------------------------
 *!
 * \brief Provide access to cs_glob_thermal_model
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
  int itherm = cs_glob_thermal_model->itherm;
  int itpscl = cs_glob_thermal_model->itpscl;

  cs_log_printf(CS_LOG_SETUP,
                ("\n"
                 "Thermal model options\n"
                 "---------------------\n\n"
                 "  Continuous phase:\n\n"));

  const char *itherm_value_str[]
    = {N_("no thermal model"),
       N_("temperature)"),
       N_("enthalpy"),
       N_("total energy")};

  const char *itpscl_value_str[]
    = {N_("none"),
       N_("temperature in Kelvin"),
       N_("temperature in Celsius")};

  cs_log_printf(CS_LOG_SETUP,
                ("    Thermal model\n"));
  cs_log_printf(CS_LOG_SETUP,
                _("    itherm:    %d (%s)\n"),
                itherm, _(itherm_value_str[itherm]));

  cs_log_printf(CS_LOG_SETUP,
                ("    Temperature scale\n"));
  cs_log_printf(CS_LOG_SETUP,
                _("    itpscl:    %d (%s)\n"),
                itpscl, _(itpscl_value_str[itpscl]));

  cs_field_t *tf = cs_thermal_model_field();
  if (tf != NULL)
    cs_log_printf
      (CS_LOG_SETUP,
       _("    Thermal variable solved: %s (field id %d)\n"),
       tf->name, tf->id);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
