/*============================================================================
 * Gas combustion model: physical properties computation
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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_boundary.h"
#include "base/cs_boundary_conditions.h"
#include "base/cs_dispatch.h"
#include "base/cs_field.h"
#include "base/cs_field_pointer.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_physical_constants.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"

#include "cogz/cs_combustion_bsh.h"
#include "cogz/cs_combustion_gas.h"
#include "cogz/cs_combustion_boundary_conditions.h"
#include "pprt/cs_combustion_model.h"
#include "pprt/cs_physical_model.h"
#include "rayt/cs_rad_transfer.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cogz/cs_combustion_physical_properties.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_combustion_physical_properties.cpp

  \brief Gas combustion model physical properties computation.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type and macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Compute physical properties for the 3-point chemistry
 *        combustion model.
 *
 * \param[in, out]   mbrom    filling indicator of romb
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_physical_properties_update_d3p(void)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  cs_rad_transfer_model_t rt_model = cs_glob_rad_transfer_params->type;

  const cs_real_t *cvar_fm = CS_F_(fm)->val;
  const cs_real_t *cvar_fp2m = CS_F_(fp2m)->val;

  cs_real_t *cpro_temp = CS_F_(t)->val;
  cs_real_t *cpro_rho = CS_F_(rho)->val;

  cs_real_t *cvar_h = nullptr;
  if (cs_glob_physical_model_flag[CS_COMBUSTION_3PT] == 1)
    cvar_h = CS_F_(h)->val;

  cs_real_t *cpro_t4m = nullptr;
  if (rt_model != CS_RAD_TRANSFER_NONE)
    cpro_t4m = cs_field_by_name("temperature_4")->val;

  cs_real_t *cpro_t2m = cs_field_by_name("temperature_2")->val;

  cs_real_t *cpro_fuel = cs_field_by_name("ym_fuel")->val;
  cs_real_t *cpro_oxyd = cs_field_by_name("ym_oxyd")->val;
  cs_real_t *cpro_prod = cs_field_by_name("ym_prod")->val;

  // Compute enthalpy defect, Kir

  cs_host_context ctx;

  const cs_real_t hinfue = cm->hinfue, hinoxy = cm->hinoxy;

  ctx.parallel_for(n_cells, [=] CS_F_HOST (cs_lnum_t c_id) {
    cs_real_t kir;
    cs_real_t phi_t[CS_BSH_NVAR_TURB];
    cs_real_t c_fm = cvar_fm[c_id];
    if (rt_model != CS_RAD_TRANSFER_NONE) {
      cs_real_t had = c_fm*hinfue + (1.-c_fm)*hinoxy;
      kir = cs_math_fmax(-(cvar_h[c_id] - had), 0.0);
    }
    else
      kir = 0.;

    cs_compute_burke_schumann_properties(c_fm, cvar_fp2m[c_id], kir, phi_t);

    cpro_rho[c_id]  = phi_t[3];
    cpro_temp[c_id] = phi_t[4];

    cpro_fuel[c_id] = phi_t[5];
    cpro_oxyd[c_id] = phi_t[6];
    cpro_prod[c_id] = phi_t[7];

    cpro_t2m[c_id] = phi_t[8];

    if (cpro_t4m != nullptr)
      cpro_t4m[c_id]  = phi_t[9];
  });

  cs_combustion_boundary_conditions_density();
}

/*----------------------------------------------------------------------------*/

END_C_DECLS