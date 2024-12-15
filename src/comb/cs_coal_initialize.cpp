/*============================================================================
 * Coal combustion model variables initialization
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

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_printf.h"

#include "cs_array.h"
#include "cs_field.h"
#include "cs_math.h"
#include "cs_mem.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_restart.h"
#include "cs_thermal_model.h"

#include "cs_coal.h"
#include "cs_coal_ht_convert.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_coal_initialize.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_coal_initialize.cpp

  \brief Coal combustion model.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type and macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize specific field values for pulverized coal combustion.
 *
 * Here, the density and viscosity have in initialized using reference values
 * or have been read from a restart file.
 *
 * The full viscosity and Cp values are available only if they have been read
 * from a restart file.
 *
 * Physical properties such as density, viscosity, Cp are handled in
 * the dedicated cs_coal_physprop function.
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_fields_initialize(void)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  const cs_coal_model_t *cm = cs_glob_coal_model;

  /* Initialization
     -------------- */

  // Massic fraction of gas
  cs_real_t *cpro_x1 = cs_field_by_name("x_c")->val;
  cs_real_t *bpro_x1 = cs_field_by_name("b_x_c")->val;

  const cs_field_t *fld_th = cs_thermal_model_field();
  cs_real_t *cvar_scalt = fld_th->val;
  cs_real_t *cvar_xch = cs_field_by_name("x_c_h")->val;

  cs_real_t *cvar_yco2 = nullptr, *cvar_hox = nullptr;

  if (cm->ieqco2 >= 1) {
    cvar_yco2 = cs_field_by_id(cm->iyco2)->val;
  }

  if (cm->ieqnox == 1) {
    cvar_hox  = cs_field_by_id(cm->ihox)->val;
  }

  /* Variable initialization
     ----------------------- */

  const cs_time_step_t *time_step = cs_glob_time_step;
  if (time_step->nt_cur == 0) {

    // All the domain is filled with the first oxidizer at TINITK

    // Computation of H1INIT and H2INIT

    cs_real_t t1init = cs_glob_fluid_properties->t0;

    // Transported variables for the mix (solid+carrying gas)^2

    cs_real_t coefe[CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS];
    for (int ige = 0; ige < CS_COMBUSTION_COAL_MAX_ELEMENTARY_COMPONENTS; ige++)
      coefe[ige] = 0.;

    // Oxidizer are mix of O2, N2 (air), CO2 and H2O (recycled exhaust)
    // the composition of the fisrt oxidiser is taken in account

    const int ioxy = 0;
    const int ico2 = cm->ico2 - 1;
    const int ih2o = cm->ih2o - 1;
    const int in2  = cm->in2 - 1;
    const int io2  = cm->io2 - 1;

    const cs_real_t dmas =   cm->wmole[io2]  * cm->oxyo2[ioxy]
                           + cm->wmole[in2]  * cm->oxyn2[ioxy]
                           + cm->wmole[ih2o] * cm->oxyh2o[ioxy]
                           + cm->wmole[ico2] * cm->oxyco2[ioxy];

    coefe[io2]  = cm->wmole[io2]  * cm->oxyo2[ioxy ] / dmas;
    coefe[ih2o] = cm->wmole[ih2o] * cm->oxyh2o[ioxy] / dmas;
    coefe[ico2] = cm->wmole[ico2] * cm->oxyco2[ioxy] / dmas;
    coefe[in2]  = 1.0 - coefe[io2] - coefe[ih2o] - coefe[ico2];

    cs_real_t h1init = cs_coal_ht_convert_t_to_h_gas_by_yi(t1init, coefe);

    cs_arrays_set_value<cs_real_t, 1>(n_cells, h1init,
                                      cvar_scalt, cvar_xch);

    // Transported variables for the mix (passive scalars, variance)

    if (cm->ieqco2 >= 1) {
      double wmco2 = cm->wmole[ico2];
      double yco2 = cm->oxyco2[ioxy] * wmco2 / dmas;

      cs_arrays_set_value<cs_real_t, 1>(n_cells, yco2, cvar_yco2);
    }

    if (cm->ieqnox == 1) {
      cs_arrays_set_value<cs_real_t, 1>(n_cells, h1init, cvar_hox);
    }

  }

  // Initialization of the continuous mass fraction
  cs_arrays_set_value<cs_real_t, 1>(n_cells, 1., cpro_x1);
  cs_arrays_set_value<cs_real_t, 1>(n_b_faces, 1., bpro_x1);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
