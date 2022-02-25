/*============================================================================
 * Mixture-length turbulence model.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

#include "cs_array.h"
#include "cs_base.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_turbulence_ml.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*! \brief Calculation of turbulent viscosity for
 *         a model of length of simple mixture
 *
 * \f[ \mu_T = \rho (\kappa L)^2 \cdot \sqrt{2 S_{ij} S_{ij}} \f]
 * \f[ S_{ij} = \dfrac{\der{u_i}{x_j} + \der{u_j}{x_i}}{2}\f]
 *
 * Edge face types are available at previous time step (except at the first time
 * step, when the itypfb and itrifb tables have not been filled).
 !*/
/*----------------------------------------------------------------------------*/

void
cs_turbulence_ml_mu_t(void)

{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  /* Initialization */

  cs_real_t *visct =  CS_F_(mu_t)->val;
  const cs_real_t *crom = CS_F_(rho)->val;

  const cs_real_t xlomlg = cs_glob_turb_rans_model->xlomlg;
  const cs_real_t coef = cs_math_pow2(cs_turb_xkappa*xlomlg)*sqrt(2.0);

  cs_real_33_t *gradv;
  BFT_MALLOC(gradv, n_cells_ext, cs_real_33_t);

  cs_field_gradient_vector(CS_F_(vel),
                           false, // no use_previous_t
                           0,    // inc
                           gradv);

  /* Compute S11^2+S22^2+S33^2+2*(S12^2+S13^2+S23^2),
     then dynamic viscosity */

# pragma omp parallel for if(n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++) {
    visct[c_id] =   cs_math_pow2(gradv[c_id][0][0])
                  + cs_math_pow2(gradv[c_id][1][1])
                  + cs_math_pow2(gradv[c_id][2][2])
                  +0.5 * (  cs_math_pow2(gradv[c_id][0][1]+gradv[c_id][1][0])
                          + cs_math_pow2(gradv[c_id][0][2]+gradv[c_id][2][0])
                          + cs_math_pow2(gradv[c_id][1][2]+gradv[c_id][2][1]));

    visct[c_id] = crom[c_id]*coef*sqrt(visct[c_id]);
  }

  /* Free temporary array */

  BFT_FREE(gradv);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
