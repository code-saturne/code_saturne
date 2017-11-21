/*============================================================================
 * Methods for lagrangian gradients
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

#include "bft_mem.h"

#include "cs_mesh.h"
#include "cs_parameters.h"

#include "cs_physical_constants.h"
#include "cs_stokes_model.h"
#include "cs_turbulence_model.h"
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
 *  - gradient of total pressure
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
  cs_lnum_t n_cells_with_ghosts = cs_glob_mesh->n_cells_with_ghosts;
  cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;

  cs_real_t   ro0 = cs_glob_fluid_properties->ro0;
  cs_real_3_t grav    = {cs_glob_physical_constants->gravity[0],
                         cs_glob_physical_constants->gravity[1],
                         cs_glob_physical_constants->gravity[2]};

  cs_real_t *wpres = NULL;

  /* Hydrostatic pressure algorithm? */
  int hyd_p_flag = cs_glob_stokes_model->iphydr;

  cs_real_3_t *f_ext = NULL;
  if (hyd_p_flag == 1)
    f_ext = (cs_real_3_t *)(cs_field_by_name("volume_forces")->val);

  cs_real_t *solved_pres
    = time_id ? extra->pressure->val_pre : extra->pressure->val;

  /* retrieve 2/3 rho^{n} k^{n} from solved pressure field for EVM models */
  // FIXME if time_id = 1, we don't have k^{n-1}
  if (   cs_glob_turb_model->itytur == 2
      || cs_glob_turb_model->itytur == 5
      || cs_glob_turb_model->itytur == 6) {
    BFT_MALLOC(wpres, n_cells_with_ghosts, cs_real_t);

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      wpres[c_id] =  solved_pres[c_id]
                   -  2./3. * extra->cromf->val[c_id]
                    * extra->cvar_k->val_pre[c_id];
    }
  }
  else {
    wpres = solved_pres;
  }

  /* Parameters for gradient computation
   * =================================== */

  int tr_dim = 0;
  cs_lnum_t inc = 1;
  cs_lnum_t iccocg = 1;
  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;
  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_var_cal_opt_t var_cal_opt;

  int key_cal_opt_id = cs_field_key_id("var_cal_opt");

  /* Get the calculation option from the pressure field */

  cs_field_get_key_struct(extra->pressure, key_cal_opt_id, &var_cal_opt);

  cs_gradient_type_by_imrgra(var_cal_opt.imrgra,
                             &gradient_type,
                             &halo_type);

  cs_real_t *weight = NULL;
  cs_internal_coupling_t  *cpl = NULL;
  int w_stride = 1;

  if (var_cal_opt.iwgrec == 1) {
    /* Weighted gradient coefficients */
    int key_id = cs_field_key_id("gradient_weighting_id");
    int diff_id = cs_field_get_key_int(extra->pressure, key_id);
    if (diff_id > -1) {
      cs_field_t *weight_f = cs_field_by_id(diff_id);
      weight = weight_f->val;
      w_stride = weight_f->dim;
    }
    /* Internal coupling structure */
    key_id = cs_field_key_id_try("coupling_entity");
    if (key_id > -1) {
      int coupl_id = cs_field_get_key_int(extra->pressure, key_id);
      if (coupl_id > -1)
        cpl = cs_internal_coupling_by_id(coupl_id);
    }
  } else if (var_cal_opt.iwgrec == 0) {
    if (var_cal_opt.idiff > 0) {
      int key_id = cs_field_key_id_try("coupling_entity");
      if (key_id > -1) {
        int coupl_id = cs_field_get_key_int(extra->pressure, key_id);
        if (coupl_id > -1)
          cpl = cs_internal_coupling_by_id(coupl_id);
      }
    }
  }

  /* Compute pressure gradient
   * ========================= */

  cs_gradient_scalar("Work array",
                     gradient_type,
                     halo_type,
                     inc,
                     iccocg,
                     var_cal_opt.nswrgr,
                     tr_dim,
                     hyd_p_flag,
                     w_stride,
                     var_cal_opt.iwarni,
                     var_cal_opt.imligr,
                     var_cal_opt.epsrgr,
                     var_cal_opt.extrag,
                     var_cal_opt.climgr,
                     f_ext,
                     extra->pressure->bc_coeffs->a,
                     extra->pressure->bc_coeffs->b,
                     wpres,
                     weight,
                     cpl,
                     gradpr);

  if (wpres != solved_pres)
    BFT_FREE(wpres);

  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] < 0) {
    for (cs_lnum_t iel = 0; iel < cs_glob_mesh->n_cells; iel++) {
      for (cs_lnum_t id = 0; id < 3; id++)
        gradpr[iel][id] += ro0 * grav[id];
    }
  }

  /* Compute velocity gradient
     ========================= */

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
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
