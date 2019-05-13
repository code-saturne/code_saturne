/*============================================================================
 * Base gas mix data.
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
#include "cs_parall.h"
#include "cs_physical_model.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gas_mix.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_gas_mix.c
        Base gas mix data.
*/
/*----------------------------------------------------------------------------*/

/*! \struct cs_gas_mix_t

  \brief Gas mix descriptor.

  Members of this structure are publicly accessible, to allow for
  concise syntax, as they are expected to be used in many places.

  \var  cs_gas_mix_t::n_species
        number of species in the gas mix
*/

/*----------------------------------------------------------------------------*/

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

/* main gas mix options and associated pointer */

static cs_gas_mix_t _gas_mix = {
  .n_species = 0,
  .sp_id_to_f_id = NULL};

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Global variables
 *============================================================================*/

const cs_gas_mix_t  *cs_glob_gas_mix = &_gas_mix;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_gas_mix_get_pointers(int **nscasp);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get pointers to members of the global gas mix structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   nscasp --> pointer to cs_glob_gas_mix->n_species
 *----------------------------------------------------------------------------*/

void
cs_f_gas_mix_get_pointers(int  **nscasp)
{
  *nscasp = &(_gas_mix.n_species);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a species field to the gas mix (set of fields).
 *
 * \param[in]   f_id   field id of an already created scalar model field
 */
/*----------------------------------------------------------------------------*/

void
cs_gas_mix_add_species(int f_id)
{
  if (cs_glob_physical_model_flag[CS_GAS_MIX] == -1)
    bft_error(__FILE__, __LINE__, 0,
              _("No gas species can be added."
                " The gas mix model is not enabled.\n"));

  cs_field_t *f = cs_field_by_id(f_id);
  if (   strcmp(f->name, "y_o2") != 0
      && strcmp(f->name, "y_n2") != 0
      && strcmp(f->name, "y_he") != 0
      && strcmp(f->name, "y_h2") != 0)
                bft_error(__FILE__, __LINE__, 0,
                          _("Only the species having the following field names "
                            "can be added to a gas mix:\n"
                            "y_o2, y_n2, y_he, y_h2\n"));

  _gas_mix.n_species++;
  BFT_REALLOC(_gas_mix.sp_id_to_f_id, _gas_mix.n_species, int);

  int sp_id = _gas_mix.n_species-1;
  _gas_mix.sp_id_to_f_id[sp_id] = f_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free array mapping gas mix species ids to field ids.
 */
/*----------------------------------------------------------------------------*/

void
cs_gas_mix_finalize(void)
{
  if (cs_glob_physical_model_flag[CS_GAS_MIX] == -1)
    bft_error(__FILE__, __LINE__, 0,
              _("The gas mix model is not enabled. Nothing to free.\n"));

  BFT_FREE(_gas_mix.sp_id_to_f_id);
  _gas_mix.n_species = 0;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
