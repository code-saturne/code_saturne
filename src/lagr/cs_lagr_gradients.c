/*============================================================================
 * Methods for lagrangian gradients
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

/*============================================================================
 * Functions dealing with lagrangian gradients
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <limits.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_mesh.h"
#include "cs_parameters.h"

#include "cs_physical_constants.h"
#include "cs_stokes_model.h"
#include "cs_field_operator.h"

#include "cs_physical_model.h"

#include "cs_lagr.h"
#include "cs_lagr_tracking.h"
#include "cs_lagr_stat.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_gradients.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond (DOXYGEN_SHOULD_SKIP_THIS) */

/*============================================================================
 * Static global variables
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute gradients.
 *
 *  - pressure gradient / rho
 *  - velocity gradient
 *
 * \param[in]   time_id   0: current time, 1: previous
 * \param[out]  gradpr    pressure gradient
 * \param[out]  gradvf    velocity gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_gradients(int            time_id,
                  cs_real_3_t   *gradpr,
                  cs_real_33_t  *gradvf)
{
  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;

  cs_real_t   ro0 = cs_glob_fluid_properties->ro0;
  cs_real_3_t grav = {cs_glob_physical_constants->gx,
                      cs_glob_physical_constants->gy,
                      cs_glob_physical_constants->gz};

  /* ====================================================================
   * 0. Parameters for gradient computation
   * ====================================================================   */

  cs_lnum_t inc = 1;
  cs_lnum_t iccocg = 1;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;
  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_var_cal_opt_t var_cal_opt;

  /* ====================================================================
   * 1. Compute pressure gradient
   * ====================================================================   */

  /* FIXME for iphydr = 1 and 2     */
  int key_cal_opt_id = cs_field_key_id("var_cal_opt");

  /* Get the calculation option from the pressure field */
  cs_field_get_key_struct(extra->pressure, key_cal_opt_id, &var_cal_opt);

  cs_gradient_type_by_imrgra(var_cal_opt.imrgra,
                             &gradient_type,
                             &halo_type);

  // FIXME not coherent with iphydr...
  cs_field_gradient_scalar(extra->pressure,
                           time_id,
                           gradient_type,
                           halo_type,
                           inc,
                           iccocg,
                           gradpr);

  /* Warning, in standard calculation, the computed pressure is
   * the hydrostatic pressure and not the real one */

  if (   cs_glob_stokes_model->iphydr == 0
      && cs_glob_physical_model_flag[CS_COMPRESSIBLE] < 0) {

    for (cs_lnum_t iel = 0; iel < cs_glob_mesh->n_cells; iel++) {

      for (int id = 0; id < 3; id++)
        gradpr[iel][id] += ro0 * grav[id];

    }

  }

  /* ====================================================================   */
  /* 2. Compute velocity gradient   */
  /* ====================================================================   */

  if (   cs_glob_lagr_time_scheme->modcpl > 0
      && cs_glob_time_step->nt_cur >= cs_glob_lagr_time_scheme->modcpl) {

    cs_field_get_key_struct(extra->vel, key_cal_opt_id, &var_cal_opt);

    cs_gradient_type_by_imrgra(var_cal_opt.imrgra,
                               &gradient_type,
                               &halo_type);

    cs_field_gradient_vector(extra->vel,
                             time_id,
                             gradient_type,
                             halo_type,
                             inc,
                             gradvf);

  }

  return;

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
