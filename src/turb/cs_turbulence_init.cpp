/*============================================================================
 * Turbulence variables initialization for various models.
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

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_array.h"
#include "base/cs_assert.h"
#include "base/cs_field.h"
#include "base/cs_field_pointer.h"
#include "alge/cs_gradient.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "base/cs_parall.h"
#include "mesh/cs_mesh_location.h"
#include "base/cs_reducers.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "turb/cs_turbulence_ke.h"
#include "turb/cs_turbulence_model.h"
#include "turb/cs_turbulence_rij.h"

#include "turb/cs_turbulence_init.h"
#include "base/cs_log_iteration.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_turbulence_init.cpp
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
  if (turb_model == nullptr)
    return;

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  const cs_real_t uref = cs_glob_turb_ref_values->uref;
  const cs_real_t almax = cs_glob_turb_ref_values->almax;

  cs_dispatch_context ctx;

  if (turb_model->itytur == 2 || turb_model->itytur == 5) {

    cs_span<cs_real_t> cvar_k = CS_F_(k)->get_vals_s();
    cs_span<cs_real_t> cvar_ep = CS_F_(eps)->get_vals_s();

    cs_real_t k_ini = -cs_math_big_r, ep_ini = -cs_math_big_r;

    if (uref >= 0.) {
      k_ini = 1.5 * cs_math_pow2(0.02*uref);
      ep_ini = pow(k_ini, 1.5) * cs_turb_cmu / almax;
    }

    ctx.parallel_for(n_cells_ext, CS_LAMBDA (cs_lnum_t e_id) {
      cvar_k[e_id] = k_ini;
      cvar_ep[e_id] = ep_ini;
    });

    if (uref >= 0.)
      cs_turbulence_ke_clip(-1, n_cells, 1);

    if (turb_model->model == CS_TURB_V2F_PHI) {
      cs_span<cs_real_t> cvar_phi = CS_F_(phi)->get_vals_s();
      cs_span<cs_real_t> cvar_fb = CS_F_(f_bar)->get_vals_s();
      ctx.parallel_for(n_cells_ext, CS_LAMBDA (cs_lnum_t e_id) {
        cvar_phi[e_id] = 2./3.;
        cvar_fb[e_id]  = 0.;
      });
    }
    else if (turb_model->model == CS_TURB_V2F_BL_V2K) {
      cs_span<cs_real_t> cvar_phi = CS_F_(phi)->get_vals_s();
      cs_span<cs_real_t> cvar_al = CS_F_(alp_bl)->get_vals_s();
      ctx.parallel_for(n_cells_ext, CS_LAMBDA (cs_lnum_t e_id) {
        cvar_phi[e_id] = 2./3.;
        cvar_al[e_id]  = 1.;
      });
    }

  }
  else if (turb_model->order == CS_TURB_SECOND_ORDER)
    cs_turbulence_rij_init_by_ref_quantities(uref, almax);

  else if (turb_model->model == CS_TURB_K_OMEGA) {
    cs_span<cs_real_t> cvar_k = CS_F_(k)->get_vals_s();
    cs_span<cs_real_t> cvar_omg = CS_F_(omg)->get_vals_s();

    cs_real_t k_ini = -cs_math_big_r, omg_ini = -cs_math_big_r;

    if (uref >= 0.) {
      /* classical eps = k^1.5/Cmu/almax and omega = eps/Cmu/k;
         values always positive, so no need to clip */
      k_ini = 1.5 * cs_math_pow2(0.02*uref);
      omg_ini = pow(k_ini, 0.5) / almax;
    }

    ctx.parallel_for(n_cells_ext, CS_LAMBDA (cs_lnum_t e_id) {
      cvar_k[e_id] = k_ini;
      cvar_omg[e_id] = omg_ini;
    });

  }

  else if (turb_model->model == CS_TURB_SPALART_ALLMARAS) {

    cs_span<cs_real_t> cvar_nusa = CS_F_(nusa)->get_vals_s();

    cs_real_t nusa_ini = -cs_math_big_r;

    if (uref >= 0.) {
      /* classical eps = k^1.5/Cmu/almax and nusa = Cmu*k**2/eps;
         values always positive, so no need to clip */
      nusa_ini = sqrt(1.5) * (0.02*uref) * almax;
    }

    cvar_nusa.set_to_val(ctx, nusa_ini);

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
  if (turb_model == nullptr)
    return n_errors;

  cs_dispatch_context ctx;

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  if (turb_model->itytur == 2 || turb_model->itytur == 5) {

    cs_span<cs_real_t> cvar_k = CS_F_(k)->get_vals_s();
    cs_span<cs_real_t> cvar_ep = CS_F_(eps)->get_vals_s();

    struct cs_data_double_n<2> rd_ke;
    struct cs_reduce_min_nr<2> reducer_ke;

    ctx.parallel_for_reduce(n_cells, rd_ke, reducer_ke, CS_LAMBDA
                            (cs_lnum_t c_id, cs_data_double_n<2> &res)
    {
      res.r[0] = cvar_k[c_id];
      res.r[1] = cvar_ep[c_id];
    });
    ctx.wait();

    cs_parall_min(2, CS_REAL_TYPE, &rd_ke.r[0]);

    cs_real_t xekmin = rd_ke.r[0], xepmin = rd_ke.r[1];

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

      struct cs_data_double_n<4> rd;
      struct cs_reduce_minmax_n<2> reducer;

      cs_span<cs_real_t> cvar_phi = CS_F_(phi)->get_vals_s();

      int n_minmax = 1;
      if (turb_model->model == CS_TURB_V2F_BL_V2K) {
        n_minmax = 2;

        cs_span<cs_real_t> cvar_al = CS_F_(alp_bl)->get_vals_s();

        ctx.parallel_for_reduce(n_cells, rd, reducer, CS_LAMBDA
                                (cs_lnum_t c_id, cs_data_double_n<4> &res)
        {
          res.r[0] = cvar_phi[c_id];
          res.r[1] = cvar_al[c_id];
          res.r[2] = res.r[0];
          res.r[3] = res.r[1];
        });
      }
      else {
        ctx.parallel_for_reduce(n_cells, rd, reducer, CS_LAMBDA
                                (cs_lnum_t c_id, cs_data_double_n<4> &res)
        {
          res.r[0] = cvar_phi[c_id];
          res.r[2] = res.r[0];
        });
      }

      ctx.wait();

      cs_parall_min(n_minmax, CS_REAL_TYPE, &rd.r[0]);
      cs_parall_max(n_minmax, CS_REAL_TYPE, &rd.r[2]);

      /* For consistency with cs_turbulence_v2f.c:_clip_v2f,
         we clip only to 0 and not to 2 */
      if (rd.r[0] < 0.) {
        n_errors += 1;
        cs_log_warning
          (_("Error in variables initialization\n"
             "----\n\n"
             "  Phi variable (v2f phi_bar or BL-v2k) out of [0, 2]bounds\n"
             "    Minimum value  = %g\n"
             "    Maximum value  = %g\n\n"
             "  Check user-defined initialization, restart data,\n"
             "  and value of reference velocity (uref)\n"),
           rd.r[0], rd.r[2]);
      }

      /* For v2-f, BL-v2/k, also check that alpha is between 0 and 1. */

      if (turb_model->model == CS_TURB_V2F_BL_V2K) {
        if (rd.r[1] < 0. || rd.r[3] > 1.) {
          n_errors += 1;
          cs_log_warning
            (_("Error in variables initialization\n"
               "----\n\n"
               "  Alpha variable (v2f BL-v2k) out of [0, 1] bounds\n"
               "    Minimum value  = %g\n"
               "    Maximum value  = %g\n\n"
               "  Check user-defined initialization, restart data,\n"
               "  and value of reference velocity (uref)\n"),
             rd.r[1], rd.r[3]);
        }
      }

    } /* End of case of itytur == 5 */

  }

  else if (turb_model->order == CS_TURB_SECOND_ORDER) {

    cs_span_2d<cs_real_t> cvar_rij = CS_F_(rij)->get_vals_t();
    cs_span<cs_real_t> cvar_ep = CS_F_(eps)->get_vals_s();

    struct cs_data_double_n<4> rd_rij;
    struct cs_reduce_min_nr<4> reducer_rij;

    struct cs_data_double_n<2> rd_alpha;
    struct cs_reduce_minmax_n<1> reducer_alpha;

    /* First compute min over Rij-epsilon */
    ctx.parallel_for_reduce(n_cells, rd_rij, reducer_rij, CS_LAMBDA
                            (cs_lnum_t c_id, cs_data_double_n<4> &res)
    {
      res.r[0] = cvar_rij(c_id, 0);
      res.r[1] = cvar_rij(c_id, 1);
      res.r[2] = cvar_rij(c_id, 2);
      res.r[3] = cvar_ep[c_id];
    });

    /* If EBRSM, compute min/max of alpha */
    if (turb_model->model == CS_TURB_RIJ_EPSILON_EBRSM) {
      cs_span<cs_real_t> cvar_al = CS_F_(alp_bl)->get_vals_s();
      ctx.parallel_for_reduce(n_cells, rd_alpha, reducer_alpha, CS_LAMBDA
                              (cs_lnum_t c_id, cs_data_double_n<2> &res)
      {
        res.r[0] = cvar_al[c_id];
        res.r[1] = cvar_al[c_id];
      });
    }

    ctx.wait();

    cs_parall_min(4, CS_REAL_TYPE, &rd_rij.r[0]);

    if (   rd_rij.r[0] >= 0 && rd_rij.r[1] >= 0
        && rd_rij.r[2] >= 0 && rd_rij.r[3] >= 0) {
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
         rd_rij.r[0], rd_rij.r[1], rd_rij.r[2], rd_rij.r[3]);
    }

    if (turb_model->model == CS_TURB_RIJ_EPSILON_EBRSM) {
      cs_real_t al_minmax[2] = {-rd_alpha.r[0], rd_alpha.r[1]};
      cs_parall_max(2, CS_REAL_TYPE, al_minmax);

      cs_real_t al_min = -al_minmax[0];
      cs_real_t al_max = al_minmax[1];

      if (al_min < 0. || al_max > 1.) {
        n_errors += 1;
        cs_log_warning
          (_("Error in variables initialization\n"
             "----\n\n"
             "  Alpha variable (EBRSM) out of [0, 1] bounds\n"
             "    Minimum value  = %g\n"
             "    Maximum value  = %g\n\n"
             "  Check user-defined initialization, restart data,\n"
             "  and value of reference velocity (uref)\n"),
           al_min, al_max);
      }
    }

  }
  else if (turb_model->model == CS_TURB_K_OMEGA) {

    cs_span<cs_real_t> cvar_k = CS_F_(k)->get_vals_s();
    cs_span<cs_real_t> cvar_omg = CS_F_(omg)->get_vals_s();

    struct cs_data_double_n<2> rd_komg;
    struct cs_reduce_min_nr<2> reducer_komg;

    ctx.parallel_for_reduce(n_cells, rd_komg, reducer_komg, CS_LAMBDA
                            (cs_lnum_t c_id, cs_data_double_n<2> &res)
    {
      res.r[0] = cvar_k[c_id];
      res.r[1] = cvar_omg[c_id];
    });

    ctx.wait();
    cs_parall_min(2, CS_REAL_TYPE, &rd_komg.r[0]);

    /* For k-omega, clip only to 0 */
    if (rd_komg.r[0] < 0. || rd_komg.r[1] <= 0.) {
      n_errors += 1;
      cs_log_warning(_("Error in variables initialization\n"
                       "----\n\n"
                       "  Negative or null turbulence\n"
                       "    Minimum value of k  = %g\n"
                       "    Minimum value of omega  = %g\n\n"
                       "  Check user-defined initialization, restart data,\n"
                       "  and value of reference velocity (uref)\n"),
                     rd_komg.r[0], rd_komg.r[1]);
    }

  }
  else if (turb_model->model == CS_TURB_SPALART_ALLMARAS) {

    cs_span<cs_real_t> cvar_nusa = CS_F_(nusa)->get_vals_s();

    struct cs_data_double_n<1> rd_sa;
    struct cs_reduce_min_nr<1> reducer_sa;

    ctx.parallel_for_reduce(n_cells, rd_sa, reducer_sa, CS_LAMBDA
                            (cs_lnum_t c_id, cs_data_double_n<1> &res)
    {
      res.r[0] = cvar_nusa[c_id];
    });
    ctx.wait();
    cs_parall_min(1, CS_REAL_TYPE, &rd_sa.r[0]);

    /* For Spalart-Allmaras, clip only to 0 */

    if (rd_sa.r[0] < 0) {
      n_errors += 1;
      cs_log_warning(_("Error in variables initialization\n"
                       "----\n\n"
                       "  Negative or null turbulence\n"
                       "    Minimum value of nu = %g\n"
                       "  Check user-defined initialization, restart data,\n"
                       "  and value of reference velocity (uref)\n"),
                     rd_sa.r[0]);
    }

  }

  return n_errors;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
