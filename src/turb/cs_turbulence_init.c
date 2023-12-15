/*============================================================================
 * Turbulence variables initialization for various models.
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
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_array.h"
#include "cs_assert.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_gradient.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_parall.h"
#include "cs_mesh_location.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_turbulence_ke.h"
#include "cs_turbulence_model.h"
#include "cs_turbulence_rij.h"

#include "cs_turbulence_init.h"
#include "cs_log_iteration.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_turbulence_init.c
        Turbulence variables initialization for various models.
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
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief First pass at initialization of turbulence variables.
 *
 * If the reference velocity is negative (incorrect or not initialized),
 * values of k, Rij, epsilon, and omega are set to -10*cs_math_big_r.
 * We will later test if those values were modified by user initialization
 * or reading a restart file.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_init_by_ref_quantities(void)
{
  const cs_turb_model_t  *turb_model = cs_glob_turb_model;
  if (turb_model == NULL)
    return;

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  const cs_real_t uref = cs_glob_turb_ref_values->uref;
  const cs_real_t almax = cs_glob_turb_ref_values->almax;

  if (turb_model->itytur == 2 || turb_model->itytur == 5) {

    cs_real_t *cvar_k = CS_F_(k)->val;
    cs_real_t *cvar_ep = CS_F_(eps)->val;

    cs_real_t k_ini = -cs_math_big_r, ep_ini = -cs_math_big_r;

    if (uref >= 0.) {
      k_ini = 1.5 * cs_math_pow2(0.02*uref);
      ep_ini = pow(k_ini, 1.5) * cs_turb_cmu / almax;
    }

    cs_array_real_set_scalar(n_cells_ext, k_ini, cvar_k);
    cs_array_real_set_scalar(n_cells_ext, ep_ini, cvar_ep);

    if (uref >= 0.)
      cs_turbulence_ke_clip(-1, n_cells, 1);

    if (turb_model->iturb == CS_TURB_V2F_PHI) {
      cs_real_t *cvar_phi = CS_F_(phi)->val;
      cs_real_t *cvar_fb = CS_F_(f_bar)->val;
      cs_array_real_set_scalar(n_cells_ext, 2./3., cvar_phi);
      cs_array_real_set_scalar(n_cells_ext, 0., cvar_fb);
    }
    else if (turb_model->iturb == CS_TURB_V2F_BL_V2K) {
      cs_real_t *cvar_phi = CS_F_(phi)->val;
      cs_real_t *cvar_al = CS_F_(alp_bl)->val;
      cs_array_real_set_scalar(n_cells_ext, 2./3., cvar_phi);
      cs_array_real_set_scalar(n_cells_ext, 1., cvar_al);
    }

  }
  else if (turb_model->itytur == 3)
    cs_turbulence_rij_init_by_ref_quantities(uref, almax);

  else if (turb_model->iturb == CS_TURB_K_OMEGA) {

    cs_real_t *cvar_k = CS_F_(k)->val;
    cs_real_t *cvar_omg = CS_F_(omg)->val;

    cs_real_t k_ini = -cs_math_big_r, omg_ini = -cs_math_big_r;

    if (uref >= 0.) {
      /* classical eps = k^1.5/Cmu/almax and omega = eps/Cmu/k;
         values always positive, so no need to clip */
      k_ini = 1.5 * cs_math_pow2(0.02*uref);
      omg_ini = pow(k_ini, 0.5) / almax;
    }

    cs_array_real_set_scalar(n_cells_ext, k_ini, cvar_k);
    cs_array_real_set_scalar(n_cells_ext, omg_ini, cvar_omg);

  }

  else if (turb_model->iturb == CS_TURB_SPALART_ALLMARAS) {

    cs_real_t *cvar_nusa = CS_F_(nusa)->val;

    cs_real_t nusa_ini = -cs_math_big_r;

    if (uref >= 0.) {
      /* classical eps = k^1.5/Cmu/almax and nusa = Cmu*k**2/eps;
         values always positive, so no need to clip */
      nusa_ini = sqrt(1.5) * (0.02*uref) * almax;
    }

    cs_array_real_set_scalar(n_cells_ext, nusa_ini, cvar_nusa);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Clipping of turbulence variables initialization
 *
 * If the user has prescribed "correct" values (in the sense that k, eps,
 * Rii, ...  > 0), we clip these values if needed for consistency.
 *
 * If the prescribed values are manifestly incorrect (negative values),
 * we assume this is an error, print warnings, and return a positive
 * error count.
 *
 * The same logic is applied to computation restarts.
 */
/*----------------------------------------------------------------------------*/

int
cs_turbulence_init_clip_and_verify(void)
{
  int n_errors = 0;

  const cs_turb_model_t  *turb_model = cs_glob_turb_model;
  if (turb_model == NULL)
    return n_errors;

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  if (turb_model->itytur == 2 || turb_model->itytur == 5) {

    cs_real_t *cvar_k = CS_F_(k)->val;
    cs_real_t *cvar_ep = CS_F_(eps)->val;

    cs_real_t f_min[2] = {HUGE_VAL, HUGE_VAL};

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      f_min[0] = cs_math_fmin(f_min[0], cvar_k[c_id]);
      f_min[1] = cs_math_fmin(f_min[1], cvar_ep[c_id]);
    }
    cs_parall_min(2, CS_REAL_TYPE, f_min);

    cs_real_t xekmin = f_min[0], xepmin = f_min[1];

    if (xekmin >= 0. && xepmin >= 0.)
      cs_turbulence_ke_clip(-1, n_cells, 1);
    else {
      n_errors += 1;
      cs_log_warning(_("Error in variables initialization\n"
                       "----\n\n"
                       "  Negative or null turbulence\n"
                       "    Minimum value of k  = %g\n"
                       "    Minimum value of epsilon  = %g\n\n"
                       "  Check user-defined initialization, restart data,\n"
                       "  and value of reference velocity (uref)\n"),
                     xekmin, xepmin);
    }

    /* For v2-f, phi-fbar or BL-v2/k, check that phi is between 0 and 2. */

    if (turb_model->itytur == 5) {

      cs_real_t *cvar_phi = CS_F_(phi)->val;

      int n_minmax = 1;
      cs_real_t v_min[2] = {HUGE_VAL, HUGE_VAL};
      cs_real_t v_max[2] = {-HUGE_VAL, -HUGE_VAL};

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        v_min[0] = cs_math_fmin(v_min[0], cvar_phi[c_id]);
        v_max[0] = cs_math_fmax(v_max[0], cvar_phi[c_id]);
      }
      if (turb_model->iturb == CS_TURB_V2F_BL_V2K) {
        cs_real_t *cvar_al = CS_F_(alp_bl)->val;
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          v_min[1] = cs_math_fmin(v_min[1], cvar_al[c_id]);
          v_max[1] = cs_math_fmax(v_max[1], cvar_al[c_id]);
        }
        n_minmax = 2;
      }

      cs_parall_min(n_minmax, CS_REAL_TYPE, v_min);
      cs_parall_max(n_minmax, CS_REAL_TYPE, v_max);

      /* For consistency with cs_turbulence_v2f.c:_clip_v2f,
         we clip only to 0 and not to 2 */
      if (v_min[0] < 0.) {
        n_errors += 1;
        cs_log_warning
          (_("Error in variables initialization\n"
             "----\n\n"
             "  Phi variable (v2f phi_bar or BL-v2k) out of [0, 2]bounds\n"
             "    Minimum value  = %g\n"
             "    Maximum value  = %g\n\n"
             "  Check user-defined initialization, restart data,\n"
             "  and value of reference velocity (uref)\n"),
           v_min[0], v_max[0]);
      }

      /* For v2-f, BL-v2/k, also check that alpha is between 0 and 1. */

      if (turb_model->iturb == CS_TURB_V2F_BL_V2K) {
        if (v_min[1] < 0. || v_max[1] > 1.) {
          n_errors += 1;
          cs_log_warning
            (_("Error in variables initialization\n"
               "----\n\n"
               "  Alpha variable (v2f BL-v2k) out of [0, 1] bounds\n"
               "    Minimum value  = %g\n"
               "    Maximum value  = %g\n\n"
               "  Check user-defined initialization, restart data,\n"
               "  and value of reference velocity (uref)\n"),
             v_min[1], v_max[1]);
        }
      }

    } /* End of case of itytur == 5 */

  }

  else if (turb_model->itytur == 3) {

    cs_real_t v_min[5] = {HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL};
    cs_real_t al_max[1] = {-HUGE_VAL};

    cs_real_6_t *cvar_rij = (cs_real_6_t *)(CS_F_(rij)->val);
    cs_real_t *cvar_ep = CS_F_(eps)->val;

    int n_minmax = 4;

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      v_min[0] = cs_math_fmin(v_min[0], cvar_rij[c_id][0]);
      v_min[1] = cs_math_fmin(v_min[1], cvar_rij[c_id][1]);
      v_min[2] = cs_math_fmin(v_min[2], cvar_rij[c_id][2]);
      v_min[3] = cs_math_fmin(v_min[3], cvar_ep[c_id]);
    }

    if (turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM) {
      cs_real_t *cvar_al = CS_F_(alp_bl)->val;
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        v_min[4] = cs_math_fmin(v_min[4], cvar_al[c_id]);
        al_max[0] = cs_math_fmax(al_max[0], cvar_al[c_id]);
      }
      n_minmax = 5;
    }

    cs_parall_min(n_minmax, CS_REAL_TYPE, v_min);

    if (   v_min[0] >= 0 && v_min[1] >= 0
        && v_min[2] >= 0 && v_min[3] >= 0) {
      if (cs_glob_turb_rans_model->irijco == 0)
        cs_turbulence_rij_clip(-1, n_cells);
    }
    else {
      n_errors += 1;
      cs_log_warning
        (_("Error in variables initialization\n"
           "----\n\n"
           "  Negative or null turbulence\n"
           "    Minimum value of R11 = %g\n"
           "    Minimum value of R22 = %g\n"
           "    Minimum value of R33 = %g\n"
           "    Maximum value of epsilon = %g\n\n"
           "  Check user-defined initialization, restart data,\n"
           "  and value of reference velocity (uref)\n"),
         v_min[0], v_min[1], v_min[3], v_min[3]);
    }

    if (turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM) {
      cs_parall_max(1, CS_REAL_TYPE, al_max);

      if (v_min[4] < 0. || al_max[0] > 1.) {
        n_errors += 1;
        cs_log_warning
          (_("Error in variables initialization\n"
             "----\n\n"
             "  Alpha variable (EBRSM) out of [0, 1] bounds\n"
             "    Minimum value  = %g\n"
             "    Maximum value  = %g\n\n"
             "  Check user-defined initialization, restart data,\n"
             "  and value of reference velocity (uref)\n"),
           v_min[4], al_max[0]);
      }
    }

  }
  else if (turb_model->iturb == CS_TURB_K_OMEGA) {

    cs_real_t *cvar_k = CS_F_(k)->val;
    cs_real_t *cvar_omg = CS_F_(omg)->val;

    cs_real_t f_min[2] = {HUGE_VAL, HUGE_VAL};

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      f_min[0] = cs_math_fmin(f_min[0], cvar_k[c_id]);
      f_min[1] = cs_math_fmin(f_min[1], cvar_omg[c_id]);
    }
    cs_parall_min(2, CS_REAL_TYPE, f_min);

    /* For k-omega, clip only to 0 */
    if (f_min[0] < 0. || f_min[1] <= 0.) {
      n_errors += 1;
      cs_log_warning(_("Error in variables initialization\n"
                       "----\n\n"
                       "  Negative or null turbulence\n"
                       "    Minimum value of k  = %g\n"
                       "    Minimum value of omega  = %g\n\n"
                       "  Check user-defined initialization, restart data,\n"
                       "  and value of reference velocity (uref)\n"),
                     f_min[0], f_min[1]);
    }

  }
  else if (turb_model->iturb == CS_TURB_SPALART_ALLMARAS) {

    cs_real_t *cvar_nusa = CS_F_(nusa)->val;

    cs_real_t nu_min = HUGE_VAL;

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      nu_min = cs_math_fmin(nu_min, cvar_nusa[c_id]);
    }
    cs_parall_min(1, CS_REAL_TYPE, &nu_min);

    /* For Spalart-Allmaras, clip only to 0 */

    if (nu_min < 0) {
      n_errors += 1;
      cs_log_warning(_("Error in variables initialization\n"
                       "----\n\n"
                       "  Negative or null turbulence\n"
                       "    Minimum value of nu = %g\n"
                       "  Check user-defined initialization, restart data,\n"
                       "  and value of reference velocity (uref)\n"),
                     nu_min);
    }

  }

  return n_errors;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
