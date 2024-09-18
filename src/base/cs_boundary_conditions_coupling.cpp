/*============================================================================
 * Update boundary conditions for thermal field.
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_1d_wall_thermal.h"
#include "cs_1d_wall_thermal_check.h"
#include "cs_boundary_conditions.h"
#include "cs_cf_thermo.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_ht_convert.h"
#include "cs_parameters.h"
#include "cs_prototypes.h"
#include "cs_rad_transfer.h"
#include "cs_thermal_model.h"
#include "cs_wall_condensation.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_boundary_conditions_coupling.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_boundary_conditions_coupling.c
        Update boundary conditions for thermal field..
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Input data for 1D wall thermal coupling.
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_coupling_t_in(void)
{
  const cs_lnum_t nfpt1d = cs_glob_1d_wall_thermal->nfpt1d;
  const cs_lnum_t *ifpt1d = cs_glob_1d_wall_thermal->ifpt1d;
  const cs_real_t *tppt1d = cs_glob_1d_wall_thermal->tppt1d;

  /* Update boundary temperature field for radiative transfer
     or wall condensation.
 */
  if (   (nfpt1d > 0)
      && (   cs_glob_rad_transfer_params->type > 0
          || cs_glob_wall_condensation->icondb > -1)) {

    cs_real_t *b_temp = CS_F_(t_b)->val;
    for (cs_lnum_t ii = 0; ii < nfpt1d; ii++) {
      const cs_lnum_t face_id = ifpt1d[ii] - 1;

      if (   cs_glob_bc_type[face_id] != CS_SMOOTHWALL
          && cs_glob_bc_type[face_id] != CS_ROUGHWALL)
        continue;

      b_temp[face_id] = tppt1d[ii];
      if (cs_glob_thermal_model->itpscl == CS_TEMPERATURE_SCALE_CELSIUS)
        b_temp[face_id] -= cs_physical_constants_celsius_to_kelvin;
    }

  }

  cs_field_build_bc_codes_all();
  cs_field_t *th_f = cs_thermal_model_field();

  for (cs_lnum_t ii = 0; ii < nfpt1d; ii++) {
    const cs_lnum_t face_id = ifpt1d[ii] - 1;

    if (  (   th_f->bc_coeffs->icodcl[face_id] != 1
           && th_f->bc_coeffs->icodcl[face_id] != 5
           && th_f->bc_coeffs->icodcl[face_id] != 6 )
        && (   cs_glob_bc_type[face_id] == CS_SMOOTHWALL
            || cs_glob_bc_type[face_id] == CS_ROUGHWALL))
      th_f->bc_coeffs->icodcl[face_id] = 5;

    th_f->bc_coeffs->rcodcl1[face_id] = tppt1d[ii];
    th_f->bc_coeffs->rcodcl2[face_id] = cs_math_infinite_r;
    th_f->bc_coeffs->rcodcl3[face_id] = 0.0;
  }

  if (cs_glob_thermal_model->itherm != CS_THERMAL_MODEL_ENTHALPY)
    return;

  for (cs_lnum_t ii = 0; ii < nfpt1d; ii++) {
    const cs_lnum_t face_id = ifpt1d[ii] - 1;
    th_f->bc_coeffs->icodcl[face_id] = -th_f->bc_coeffs->icodcl[face_id];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief output data for 1D wall thermal coupling
 *
 * \param[in, out]  hbord  exchange coefficients for boundary
 * \param[in, out]  tbord  boundary temperature
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_coupling_t_out(cs_real_t  hbord[],
                                      cs_real_t  tbord[])
{
  const cs_mesh_t *m = cs_glob_mesh;

  const cs_lnum_t nfpt1d = cs_glob_1d_wall_thermal->nfpt1d;
  const cs_lnum_t *ifpt1d = cs_glob_1d_wall_thermal->ifpt1d;

  /* Conversion to temperature for enthalpy or energy
   *  (check for surface couplings to make sure it is needed)
   *
   * In enthalpy formulation, transform to temperatures for SYRTHES
   *  To conserve flux Phi = (lambda/d     ) Delta T
   *                or Phi = (lambda/(d Cp)) Delta H
   *  recall      hbord = lambda/d. */

  if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_ENTHALPY) {

    cs_real_t *wa = nullptr;

    // Temperature near boundary faces
    BFT_MALLOC(wa, m->n_b_faces, cs_real_t);
    cs_ht_convert_h_to_t_faces(tbord, wa);

    for (cs_lnum_t ii = 0; ii < nfpt1d; ii++) {
      const cs_lnum_t face_id = ifpt1d[ii] - 1;
      tbord[face_id] = wa[face_id];
    }
    BFT_FREE(wa);

  }
  else if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_TOTAL_ENERGY) {

    const cs_real_t cv0 = cs_glob_fluid_properties->cv0;
    cs_real_t *wa = nullptr, *cpro_cv = nullptr;
    const cs_real_t *cpro_rho = CS_F_(rho)->val;
    const cs_real_3_t *vel = (const cs_real_3_t *)CS_F_(vel)->val;

    const int icv = cs_glob_fluid_properties->icv;
    if (icv > -1)
      cpro_cv = cs_field_by_id(icv)->val;

    // Epsilon sup for perfect gas at cells
    BFT_MALLOC(wa, m->n_cells_with_ghosts, cs_real_t);
    cs_cf_thermo_eps_sup(cpro_rho, wa, m->n_cells);

    for (cs_lnum_t ii = 0; ii < nfpt1d; ii++) {
      const cs_lnum_t face_id = ifpt1d[ii] - 1;
      const cs_lnum_t c_id = m->b_face_cells[face_id];
      const cs_real_t cvt
        = tbord[face_id] - (0.5*cs_math_3_square_norm(vel[c_id]) - wa[c_id]);
      if (cpro_cv != nullptr)
        tbord[face_id] = cvt/cpro_cv[c_id];
      else
        tbord[face_id] = cvt/cv0;
    }

    BFT_FREE(wa);
  }

  // Update external boudary condition
  cs_user_1d_wall_thermal(3);

  cs_1d_wall_thermal_check(3);

  // Coupling with radiative transfer
  if (cs_glob_rad_transfer_params->type > 0) {
    for (cs_lnum_t ii = 0; ii < nfpt1d; ii++) {
      const cs_lnum_t face_id = ifpt1d[ii] - 1;

      if (   cs_glob_bc_type[face_id] != CS_SMOOTHWALL
          && cs_glob_bc_type[face_id] != CS_ROUGHWALL)
        continue;

      cs_1d_wall_thermal_solve(ii, tbord[face_id], hbord[face_id]);
    }
  }
  // Without coupling with radiative transfer
  else {
    for (cs_lnum_t ii = 0; ii < nfpt1d; ii++) {
      const cs_lnum_t face_id = ifpt1d[ii] - 1;
      cs_1d_wall_thermal_solve(ii, tbord[face_id], hbord[face_id]);
    }

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
