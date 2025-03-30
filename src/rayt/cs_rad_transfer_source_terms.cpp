/*============================================================================
 * Radiation solver source term integration for thermal scalar
 *============================================================================*/

/* This file is part of code_saturne, a general-purpose CFD tool.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA. */

/*----------------------------------------------------------------------------*/

#include "base/cs_defs.h"
#include "base/cs_math.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_mem.h"
#include "bft/bft_printf.h"

#include "base/cs_field_pointer.h"
#include "base/cs_parameters.h"
#include "base/cs_physical_constants.h"
#include "pprt/cs_physical_model.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "rayt/cs_rad_transfer.h"
#include "base/cs_thermal_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "rayt/cs_rad_transfer_source_terms.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file cs_rad_transfer_source_terms.cpp */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions for fortran API
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
 * \brief Implicit and explicit radiative source terms for thermal scalar.
 *
 * \param[in,out]  rhs     work array for right hand side
 * \param[in,out]  fimp    work array for unsteady term
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_source_terms(cs_real_t  rhs[],
                             cs_real_t  fimp[])
{
  if (      cs_glob_thermal_model->thermal_variable
         == CS_THERMAL_MODEL_TEMPERATURE
      ||    cs_glob_thermal_model->thermal_variable
         == CS_THERMAL_MODEL_ENTHALPY) {

   const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
   const cs_real_t *cell_vol = cs_glob_mesh_quantities->cell_vol;

    /* Implicit part   */
    cs_real_t *rad_st_impl = CS_FI_(rad_ist, 0)->val;
    /* Explicit part   */
    cs_real_t *rad_st_expl = CS_FI_(rad_est, 0)->val;

    /* Conversion Temperature -> potential Temperature (theta) */
    if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 0) {
      const cs_real_t *pottemp = CS_F_(t)->val;
      const cs_real_t *tempc = cs_field_by_name("real_temperature")->val;
      const cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        rad_st_impl[c_id] = cs::max(-rad_st_impl[c_id], 0.0);
        cs_real_t cor_factor = pottemp[c_id] / (tempc[c_id] + tkelvi);
        fimp[c_id] += cor_factor * rad_st_impl[c_id] * cell_vol[c_id];
        rhs[c_id] += cor_factor * rad_st_expl[c_id] * cell_vol[c_id];
      }
    }
    else {
      /* FIXME ? in Renuda commit fro rfsck model, implicit part
         (rad_st_imp and fimp were ignored when
         cs_glob_rad_transfer_params->imrcfsk == 1
         but this did not seem consistent */
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        rad_st_impl[c_id] = cs::max(-rad_st_impl[c_id], 0.0);
        fimp[c_id] += rad_st_impl[c_id] * cell_vol[c_id];
        rhs[c_id] += rad_st_expl[c_id] * cell_vol[c_id];
      }
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
