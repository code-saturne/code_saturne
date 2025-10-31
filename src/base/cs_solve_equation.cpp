/*============================================================================
 * Solve the convection diffusion equation.
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

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_array.h"
#include "base/cs_assert.h"
#include "atmo/cs_at_data_assim.h"
#include "atmo/cs_atmo.h"
#include "atmo/cs_atmo_aerosol.h"
#include "atmo/cs_atmo_source_terms.h"
#include "atmo/cs_atmo_profile_std.h"
#include "alge/cs_blas.h"
#include "base/cs_boundary_conditions.h"
#include "alge/cs_bw_time_diff.h"
#include "cfbl/cs_cf_model.h"
#include "comb/cs_coal.h"
#include "pprt/cs_combustion_model.h"
#include "ctwr/cs_ctwr.h"
#include "ctwr/cs_ctwr_source_terms.h"
#include "alge/cs_divergence.h"
#include "base/cs_dispatch.h"
#include "base/cs_drift_convective_flux.h"
#include "elec/cs_elec_model.h"
#include "base/cs_equation_iterative_solve.h"
#include "alge/cs_face_viscosity.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_operator.h"
#include "base/cs_field_pointer.h"
#include "alge/cs_gradient.h"
#include "gui/cs_gui.h"
#include "atmo/cs_intprf.h"
#include "lagr/cs_lagr.h"
#include "lagr/cs_lagr_precipitation_model.h"
#include "base/cs_math.h"
#include "base/cs_mass_source_terms.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_parall.h"
#include "base/cs_param_types.h"
#include "base/cs_physical_constants.h"
#include "pprt/cs_physical_model.h"
#include "base/cs_prototypes.h"
#include "rayt/cs_rad_transfer.h"
#include "rayt/cs_rad_transfer_source_terms.h"
#include "base/cs_sat_coupling.h"
#include "base/cs_scalar_clipping.h"
#include "base/cs_syr_coupling.h"
#include "base/cs_thermal_model.h"
#include "base/cs_time_step.h"
#include "turb/cs_turbulence_model.h"
#include "turb/cs_turbulence_rij.h"
#include "turb/cs_turbulence_rit.h"
#include "base/cs_velocity_pressure.h"
#include "base/cs_vof.h"
#include "base/cs_volume_mass_injection.h"
#include "base/cs_wall_condensation.h"
#include "base/cs_wall_functions.h"

#include "cogz/cs_combustion_gas.h"
#include "cogz/cs_steady_laminar_flamelet_source_terms.h"
#include "cogz/cs_soot_model.h"
#include "comb/cs_coal_source_terms.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_solve_equation.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_solve_equation.cpp

  \brief This subroutine performs the solving the convection/diffusion
         equation (with eventually source terms and/or drift) for a field
         quantity over a time step.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type and macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Prototypes for Fortran functions and variables.
 *============================================================================*/

void
cs_f_ebutss(int        f_id,
            cs_real_t  smbrs[],
            cs_real_t  rovsdt[]);

void
cs_f_lwctss(int        f_id,
            cs_real_t  smbrs[],
            cs_real_t  rovsdt[]);

/* Use legacy macro type to maintain compatibility with legacy user files */

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*============================================================================
 * Fortran function prototypes for subroutines from field.f90.
 *============================================================================*/

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize turbulent diffusivity for SGDH model.
 *----------------------------------------------------------------------------*/

static void
_init_sgdh_diff(const cs_field_t *f,
                cs_real_t        *sgdh_diff)
{
  const cs_real_t *visct = CS_F_(mu_t)->val;

  /* Retrieve turbulent Schmidt value for current scalar */
  const int ksigmas = cs_field_key_id("turbulent_schmidt");
  const cs_real_t turb_schmidt = cs_field_get_key_double(f, ksigmas);

  /* If turbulent Schmidt is variable, id of the corresponding field */
  cs_real_t *cpro_turb_schmidt = nullptr;
  const int key_turb_schmidt = cs_field_key_id("turbulent_schmidt_id");
  const int t_scd_id = cs_field_get_key_int(f, key_turb_schmidt);
  if (t_scd_id > -1)
    cpro_turb_schmidt = cs_field_by_id(t_scd_id)->val;

  /* Retrieve turbulent diffusivity value for current scalar */
  cs_real_t *cpro_turb_diff = nullptr;
  const int key_turb_diff = cs_field_key_id("turbulent_diffusivity_id");
  const int t_dif_id = cs_field_get_key_int(f, key_turb_diff);

  /* Variable turbulent diffusivity */
  if (t_dif_id > -1) {
    cpro_turb_diff = cs_field_by_id(t_dif_id)->val;
    cs_array_real_copy(cs_glob_mesh->n_cells, cpro_turb_diff, sgdh_diff);
  }
  /* Variable Schmidt number */
  else if (cpro_turb_schmidt != nullptr)
    for (cs_lnum_t c_id = 0; c_id < cs_glob_mesh->n_cells; c_id++)
      sgdh_diff[c_id] =   cs_math_fmax(visct[c_id], 0.)
                        / cpro_turb_schmidt[c_id];
  else
    for (cs_lnum_t c_id = 0; c_id < cs_glob_mesh->n_cells; c_id++)
      sgdh_diff[c_id] =   cs_math_fmax(visct[c_id], 0.)
                        / turb_schmidt;
}

/*----------------------------------------------------------------------------
 * Compute production and dissipation terms for the variance of another scalar.
 *
 * TODO this should perhaps be moved to turbulence modelling.
 *----------------------------------------------------------------------------*/

static void
_production_and_dissipation_terms(const cs_field_t  *f,
                                  const int          ifcvsl,
                                  const int          st_prv_id,
                                  const cs_real_t    rvarfl,
                                  const cs_real_t    visls_0,
                                  const cs_real_t    thetv,
                                  const cs_real_t    dt[],
                                  const cs_real_t   *xcpp,
                                  const cs_real_t   *crom,
                                  const cs_real_t   *viscl,
                                  const cs_real_t   *cvar_var,
                                  const cs_real_t   *sgdh_diff,
                                  const cs_real_t   *cpro_viscls,
                                  cs_real_t         *rhs,
                                  cs_real_t         *fimp,
                                  cs_real_t         *cpro_st,
                                  cs_real_t         *cpro_tsscal)
{
  const cs_turb_model_t *turb_model = cs_glob_turb_model;

  if (! (   turb_model->itytur == 2
         || turb_model->itytur == 3
         || turb_model->itytur == 5
         || turb_model->model == CS_TURB_K_OMEGA))
    return;

  const cs_field_t *f_fm = nullptr;  /* First moment field of
                                     which f is the variance */
  {
    const int kscavr = cs_field_key_id("first_moment_id");
    const int ivarsc = cs_field_get_key_int(f, kscavr);

    if (ivarsc > -1)
      f_fm = cs_field_by_id(ivarsc);

    else {
      bft_error(__FILE__, __LINE__, 0,
                _("%s: field %s is not a variance of another scalar\n"
                  "The calculation will not be run."),
                __func__, f->name);
      return;
    }
  }

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  const cs_real_t *cell_f_vol = cs_glob_mesh_quantities->cell_vol;

  const cs_real_t *cvara_var = f->val_pre;
  cs_real_t *cvara_var_fm = f_fm->val_pre;

  /* Get the turbulent flux model for the scalar and its variance */

  const int kturt = cs_field_key_id("turbulent_flux_model");

  int variance_turb_flux_model = 0, variance_turb_flux_model_type = 0;

  if (f_fm != nullptr) {
    variance_turb_flux_model = cs_field_get_key_int(f_fm, kturt);
    variance_turb_flux_model_type = variance_turb_flux_model / 10;
  }

  /* Homogeneous Neumann on convective inlet on the production term
     for the variance */

  const cs_real_t *coefap = f_fm->bc_coeffs->a;
  const cs_real_t *coefbp = f_fm->bc_coeffs->b;

  cs_field_bc_coeffs_t bc_coeffs_loc;
  cs_field_bc_coeffs_init(&bc_coeffs_loc);
  CS_MALLOC(bc_coeffs_loc.a, n_b_faces, cs_real_t);
  CS_MALLOC(bc_coeffs_loc.b, n_b_faces, cs_real_t);

  cs_real_t *coefa_p = bc_coeffs_loc.a;
  cs_real_t *coefb_p = bc_coeffs_loc.b;

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    coefa_p[face_id] = coefap[face_id];
    coefb_p[face_id] = coefbp[face_id];
    if (cs_glob_bc_type[face_id] == CS_CONVECTIVE_INLET) {
      coefa_p[face_id] = 0.;
      coefb_p[face_id] = 1.;
    }
  }

  cs_real_3_t *grad;
  CS_MALLOC(grad, n_cells_ext, cs_real_3_t);

  cs_field_gradient_scalar_array(f_fm->id,
                                 1,     /* inc */
                                 &bc_coeffs_loc,
                                 cvara_var_fm,
                                 grad);

  CS_FREE(coefa_p);
  CS_FREE(coefb_p);

  /* Production Term
     --------------- */

  cs_field_t *f_produc = cs_field_by_double_composite_name_try
                           ("algo:", f->name, "_production");

  /* NB: diffusivity is clipped to 0 because in LES, it might be negative.
   * Problematic
   * for the variance even if it is strange to use variance and LES...
   * and a test above prohibits use of this function in LES */

  /* Time extrapolation (2nd order) */
  if (st_prv_id >= 0) {

    if (variance_turb_flux_model_type >= 1) {
      const cs_real_3_t *xut
        = (const cs_real_3_t *)cs_field_by_composite_name(f_fm->name,
                                                          "turbulent_flux")->val;
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const cs_real_t prod = - 2 * cs_math_3_dot_product(grad[c_id],
                                                           xut[c_id]);
        if (f_produc != nullptr )
          f_produc->val[c_id] = prod;
        cpro_st[c_id] += xcpp[c_id] * cell_f_vol[c_id] * crom[c_id]
                         *  prod;

      }
    }
    /* SGDH model */
    else {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const cs_real_t prod = 2. * sgdh_diff[c_id] / crom[c_id]
                         * cs_math_3_dot_product(grad[c_id], grad[c_id]);
        if (f_produc != nullptr )
          f_produc->val[c_id] = prod;
        cpro_st[c_id] += xcpp[c_id] * cell_f_vol[c_id] * crom[c_id]
                         *  prod;
      }
    }
  }

  /* Without time extrapolation... */
  else {

    if (variance_turb_flux_model_type >= 1) {
      const cs_real_3_t *xut
        = (const cs_real_3_t *)cs_field_by_composite_name(f_fm->name,
                                                          "turbulent_flux")->val;
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const cs_real_t cprovol = xcpp[c_id] * cell_f_vol[c_id] * crom[c_id];
        /* Special time stepping to ensure positivity of the variance */
        const cs_real_t prod = -2*cs_math_3_dot_product(grad[c_id], xut[c_id]);
        if (f_produc != nullptr )
          f_produc->val[c_id] = prod;
        rhs[c_id] += cs_math_fmax(prod * cprovol, 0.);

        /* Implicit "production" term when negative, but check if the
         * variance is non-zero */
        if (   (cvar_var[c_id] > cs_math_epzero*cs_math_fabs(prod*dt[c_id]))
            && (prod*cprovol < 0.)) {
          fimp[c_id] -= prod * cprovol / cvar_var[c_id];
          rhs[c_id] += prod * cprovol;
        }
      }
    }
    /* SGDH model */
    else {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const cs_real_t prod =   2 * sgdh_diff[c_id] / crom[c_id]
                            * cs_math_3_dot_product(grad[c_id], grad[c_id]);
        if (f_produc != nullptr )
          f_produc->val[c_id] = prod;
        rhs[c_id] +=   xcpp[c_id] * cell_f_vol[c_id] * crom[c_id]
                       * prod;
      }
    }

    /* Production term for a variance  TODO compute ustdy when isso2t >0 */
    if (cs_glob_velocity_pressure_model->idilat >= 4)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cpro_tsscal[c_id] +=   2*xcpp[c_id] * sgdh_diff[c_id] * cell_f_vol[c_id]
                             * cs_math_3_dot_product(grad[c_id], grad[c_id]);
      }
  }

  CS_FREE(grad);

  /* Dissipation term
     ---------------- */

  cs_field_t *f_dissip = cs_field_by_double_composite_name_try
    ("algo:", f->name, "_dissipation");

  cs_real_6_t *cvara_rij = nullptr;
  cs_real_t *cvara_k= nullptr, *cvara_ep = nullptr;
  cs_real_t *cvara_omg= nullptr, *cvar_al = nullptr;

  const cs_real_t thetap = (st_prv_id >= 0) ? thetv : 1;

  if (turb_model->itytur == 2 || turb_model->itytur == 5) {
    cvara_k = CS_F_(k)->val_pre;
    cvara_ep = CS_F_(eps)->val_pre;
  }
  else if (turb_model->itytur == 3) {
    cvara_ep = CS_F_(eps)->val_pre;
    cvara_rij = (cs_real_6_t*)CS_F_(rij)->val_pre;
    if (   variance_turb_flux_model == 11
        || variance_turb_flux_model == 21
        || variance_turb_flux_model == 31) {
      cvar_al = cs_field_by_composite_name(f_fm->name, "alpha")->val;
    }
  }
  else if (turb_model->model == CS_TURB_K_OMEGA) {
    cvara_k = CS_F_(k)->val_pre;
    cvara_omg = CS_F_(omg)->val_pre;
  }

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    cs_real_t xe = 0, xk = 0, alpha_theta = 1;

    if (turb_model->itytur == 2 || turb_model->itytur == 5) {
      xk = cvara_k[c_id];
      xe = cvara_ep[c_id];
    }
    else if (turb_model->itytur == 3) {
      xe = cvara_ep[c_id];
      xk = 0.5 * (cvara_rij[c_id][0] + cvara_rij[c_id][1] + cvara_rij[c_id][2]);
      /* Implicit term -1/Rh * Eps/k * Variance
       *  with
       * Rh = R = 0.8 for SGDH
       * Rh = (1-alpha_T) * Pr + R * alpha_T
       * with - R = 0.5
       *      - alpha_T = 1.0 for GGDH/DFM/AFM */
    }
    else if (turb_model->model == CS_TURB_K_OMEGA) {
      xk = cvara_k[c_id];
      xe = cs_turb_cmu*xk*cvara_omg[c_id];
    }

    if (cvar_al != nullptr)
      alpha_theta = cvar_al[c_id];

    cs_real_t prdtl= viscl[c_id] * xcpp[c_id];
    if (ifcvsl > -1)
      prdtl /= cpro_viscls[c_id];
    else
      prdtl /= visls_0;

    const cs_real_t xr = (1.0 - alpha_theta)*prdtl + alpha_theta*rvarfl;
    const cs_real_t cprovol = xcpp[c_id] * crom[c_id] * cell_f_vol[c_id];
    const cs_real_t dissip_freq = xe / (xk * xr);
    const cs_real_t dissip = dissip_freq * cvara_var[c_id];
    if (f_dissip != nullptr)
      f_dissip->val[c_id] = dissip;
    /* The diagonal receives eps/Rk, (*theta possibly) */
    fimp[c_id] += dissip_freq * cprovol * thetap;
    /* The right hand side receives the dissipation */
    rhs[c_id] -= dissip * cprovol;
  }

}

/*----------------------------------------------------------------------------
 * Compute velocity of diffusion facet for scalar.
 *
 * The caller is responsible for freeing the returned weighb, weighf, and
 * viscce arrays.
 *----------------------------------------------------------------------------*/

static void
_diffusion_terms_scalar(const cs_field_t           *f,
                        const cs_equation_param_t  *eqp,
                        cs_real_t                   visls_0,
                        const cs_real_t             xcpp[],
                        const cs_real_t             sgdh_diff[],
                        const cs_real_t            *cpro_viscls,
                        cs_real_t                   viscf[],
                        cs_real_t                   viscb[],
                        cs_real_t                   rhs[],
                        cs_real_t                 **weighb,
                        cs_real_2_t               **weighf,
                        cs_real_6_t               **viscce)
{
  const cs_turb_model_t *turb_model = cs_glob_turb_model;

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  /* Get the turbulent flux model for the scalar */
  const int kturt = cs_field_key_id("turbulent_flux_model");
  const int scalar_turb_flux_model = cs_field_get_key_int(f, kturt);
  const int scalar_turb_flux_model_type = scalar_turb_flux_model / 10;

  cs_real_6_t *vistet;
  CS_MALLOC_HD(vistet, n_cells_ext, cs_real_6_t, cs_alloc_mode);

  cs_real_t *_weighb = nullptr;
  cs_real_2_t *_weighf = nullptr;
  cs_real_6_t *_viscce = nullptr;

  /* AFM model or DFM models: add div(Cp*rho*T'u') to rhs
   * Compute T'u' for GGDH */
  if (scalar_turb_flux_model_type >= 1) {
    if (   scalar_turb_flux_model == 11
        || scalar_turb_flux_model == 21
        || scalar_turb_flux_model == 31) {
      cs_field_t *f_alpha = cs_field_by_composite_name(f->name, "alpha");
      cs_turbulence_rij_solve_alpha(f_alpha->id, -1, cs_turb_xclt);
    }
    cs_turbulence_rij_transport_div_tf(f->id, xcpp, vistet, rhs);
  }

  /* Scalar diffusivity */
  cs_real_t *cpro_wgrec_s = nullptr;
  cs_real_6_t *cpro_wgrec_v = nullptr;

  /* weighting field for gradient */
  if (eqp->iwgrec == 1) {
    const int kwgrec = cs_field_key_id_try("gradient_weighting_id");
    const int iflwgr = cs_field_get_key_int(f, kwgrec);
    cs_field_t *f_g = cs_field_by_id(iflwgr);
    if (f_g->dim > 1)
      cpro_wgrec_v = (cs_real_6_t *)f_g->val;
    else
      cpro_wgrec_s = f_g->val;
  }

  if (eqp->idften & CS_ISOTROPIC_DIFFUSION) {
    cs_real_t *w1;
    CS_MALLOC_HD(w1, n_cells_ext, cs_real_t, cs_alloc_mode);

    int idifftp = eqp->idifft;
    if (scalar_turb_flux_model_type == 3)
      idifftp = 0;

    if (cpro_viscls == nullptr)
      cs_array_set_value_real(n_cells, 1, visls_0, w1);
    else
      cs_array_real_copy(n_cells, cpro_viscls, w1);

    if (idifftp > 0) {
      if (f == CS_F_(fsm)) {
        // TODO check if we can simply set visls_0 to 0 for this- field instead.
        cs_array_set_value_real(n_cells, 1, 0., w1);
      }

#     pragma omp parallel for if(n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        w1[c_id] += xcpp[c_id]*sgdh_diff[c_id];
      }
    }

    if (cpro_wgrec_s != nullptr) {
      cs_array_real_copy(n_cells, w1, cpro_wgrec_s);
      cs_halo_sync_var(m->halo, CS_HALO_STANDARD, cpro_wgrec_s);
    }

    cs_face_viscosity(m,
                      fvq,
                      eqp->imvisf,
                      w1,
                      viscf,
                      viscb);


    CS_FREE(w1);
  }

  /* Symmetric tensor diffusivity (GGDH) */
  else if (eqp->idften & CS_ANISOTROPIC_DIFFUSION) {
    CS_MALLOC_HD(_weighb, n_b_faces, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(_weighf, n_i_faces, cs_real_2_t, cs_alloc_mode);
    CS_MALLOC_HD(_viscce, n_cells_ext, cs_real_6_t, cs_alloc_mode);
    const cs_real_6_t *visten = nullptr;
    const int kctheta = cs_field_key_id("turbulent_flux_ctheta");
    const cs_real_t ctheta = cs_field_get_key_double(f, kctheta);

    if (turb_model->model != CS_TURB_RIJ_EPSILON_EBRSM) {
      cs_field_t * f_vis = cs_field_by_name("anisotropic_turbulent_viscosity");
      visten = (cs_real_6_t *)f_vis->val;
    }
    /* EBRSM and (GGDH or AFM) */
    else {
      cs_field_t * f_vis
        = cs_field_by_name("anisotropic_turbulent_viscosity_scalar");
      visten = (cs_real_6_t *)f_vis->val;
    }

    if (   (scalar_turb_flux_model == 11)
        || (scalar_turb_flux_model == 21)
        || (scalar_turb_flux_model == 20)) {
      if (cpro_viscls == nullptr) {
#       pragma omp parallel for if(n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          const cs_real_t tmp = eqp->idifft*xcpp[c_id];
          for (cs_lnum_t ii = 0; ii < 3; ii++)
            _viscce[c_id][ii] = tmp*vistet[c_id][ii] + visls_0;
          for (cs_lnum_t ii = 3; ii < 6; ii++)
            _viscce[c_id][ii] = tmp*vistet[c_id][ii];
        }
      }
      else {
#       pragma omp parallel for if(n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          const cs_real_t tmp = eqp->idifft*xcpp[c_id];
          for (cs_lnum_t ii = 0; ii < 3; ii++)
            _viscce[c_id][ii] = tmp*vistet[c_id][ii] + cpro_viscls[c_id];
          for (cs_lnum_t ii = 3; ii < 6; ii++)
            _viscce[c_id][ii] = tmp*vistet[c_id][ii];
        }
      }
    }
    else {
      if (cpro_viscls == nullptr) {
#       pragma omp parallel for if(n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          const cs_real_t temp = eqp->idifft*xcpp[c_id]*ctheta/cs_turb_csrij;
          for (cs_lnum_t ii = 0; ii < 3; ii++)
            _viscce[c_id][ii] = temp*visten[c_id][ii] + visls_0;
          for (cs_lnum_t ii = 3; ii < 6; ii++)
            _viscce[c_id][ii] = temp*visten[c_id][ii];
        }
      }
      else {
#       pragma omp parallel for if(n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          const cs_real_t tmp = eqp->idifft*xcpp[c_id]*ctheta/cs_turb_csrij;
          for (cs_lnum_t ii = 0; ii < 3; ii++)
            _viscce[c_id][ii] = tmp*visten[c_id][ii] + cpro_viscls[c_id];
          for (cs_lnum_t ii = 3; ii < 6; ii++)
            _viscce[c_id][ii] = tmp*visten[c_id][ii];
        }
      }
    }

    /* Weighting for gradient */
    if (cpro_wgrec_v != nullptr) {
      cs_array_real_copy(n_cells*6,
                         (const cs_real_t *)_viscce,
                         (cs_real_t *)cpro_wgrec_v);
      cs_halo_sync_r(m->halo, false, cpro_wgrec_v);
    }

    cs_face_anisotropic_viscosity_scalar(m,
                                         fvq,
                                         _viscce,
                                         eqp->verbosity,
                                         _weighf,
                                         _weighb,
                                         viscf,
                                         viscb);
  }

  CS_FREE_HD(vistet);

  *weighb = _weighb;
  *weighf = _weighf;
  *viscce = _viscce;
}

/*----------------------------------------------------------------------------
 * Compute velocity of diffusion facet for vector.
 *
 * The caller is responsible for freeing the returned weighb, weighf, and
 * viscce arrays.
 *----------------------------------------------------------------------------*/

static void
_diffusion_terms_vector(const cs_field_t            *f,
                        const cs_equation_param_t   *eqp,
                        cs_real_t                    viscf[],
                        cs_real_t                    viscb[],
                        cs_real_t                  **weighb,
                        cs_real_2_t                **weighf,
                        cs_real_6_t                **viscce)
{
  const cs_turb_model_t *turb_model = cs_glob_turb_model;

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  cs_real_t *_weighb = nullptr;
  cs_real_2_t *_weighf = nullptr;
  cs_real_6_t *_viscce = nullptr;

  const cs_real_t *visct = CS_F_(mu_t)->val;

  /* Viscosity and diffusivity*/
  const cs_real_t *cpro_viscls = nullptr;
  const int kivisl = cs_field_key_id("diffusivity_id");
  const int kvisl0 = cs_field_key_id("diffusivity_ref");
  const int ifcvsl = cs_field_get_key_int(f, kivisl);
  const cs_real_t visls_0 = cs_field_get_key_double(f, kvisl0);
  if (ifcvsl > -1)
    cpro_viscls = cs_field_by_id(ifcvsl)->val;

  /* Weighting field for gradient */
  cs_real_t *cpro_wgrec_s = nullptr;
  cs_real_6_t *cpro_wgrec_v = nullptr;

  if (eqp->iwgrec == 1) {
    const int kwgrec = cs_field_key_id_try("gradient_weighting_id");
    const int iflwgr = cs_field_get_key_int(f, kwgrec);
    cs_field_t *f_g = cs_field_by_id(iflwgr);
    if (f_g->dim > 1)
      cpro_wgrec_v = (cs_real_6_t *)f_g->val;
    else
      cpro_wgrec_s = f_g->val;
  }

  /* Scalar diffusivity */
  if (eqp->idften & CS_ISOTROPIC_DIFFUSION) {

    const int ksigmas = cs_field_key_id("turbulent_schmidt");
    const cs_real_t turb_schmidt = cs_field_get_key_double(f, ksigmas);

    const int idifftp = eqp->idifft;
    cs_real_t *w1;
    CS_MALLOC(w1, n_cells_ext, cs_real_t);

    if (cpro_viscls == nullptr) {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        w1[c_id] =   visls_0
                   + idifftp * fmax(visct[c_id], 0) / turb_schmidt;
      }
    }
    else {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        w1[c_id] =   cpro_viscls[c_id]
                   + idifftp * fmax(visct[c_id], 0) / turb_schmidt;
      }
    }

    /* Weighting for gradient */
    if (cpro_wgrec_s != nullptr) {
      cs_array_real_copy(n_cells, w1, cpro_wgrec_s);
      cs_halo_sync_var(m->halo, CS_HALO_STANDARD, cpro_wgrec_s);
    }

    cs_face_viscosity(m,
                      fvq,
                      eqp->imvisf,
                      w1,
                      viscf,
                      viscb);

    CS_FREE(w1);
  }

  /* Symmetric tensor diffusivity (GGDH) */
  else if (eqp->idften & CS_ANISOTROPIC_DIFFUSION) {

    CS_MALLOC_HD(_weighb, n_b_faces, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(_weighf, n_i_faces, cs_real_2_t, cs_alloc_mode);
    CS_MALLOC_HD(_viscce, n_cells_ext, cs_real_6_t, cs_alloc_mode);

    const cs_real_6_t *visten = nullptr;
    const int kctheta = cs_field_key_id("turbulent_flux_ctheta");
    const cs_real_t ctheta = cs_field_get_key_double(f, kctheta);

    if (turb_model->model != CS_TURB_RIJ_EPSILON_EBRSM) {
      cs_field_t * f_vis = cs_field_by_name("anisotropic_turbulent_viscosity");
      visten = (cs_real_6_t *)f_vis->val;
    }

    /* EBRSM and (GGDH or AFM) */
    else {
      cs_field_t * f_vis
        = cs_field_by_name("anisotropic_turbulent_viscosity_scalar");
      visten = (cs_real_6_t *)f_vis->val;
    }

    if (cpro_viscls == nullptr) {
      const cs_real_t tmp = eqp->idifft*ctheta/cs_turb_csrij;
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        for (cs_lnum_t ii = 0; ii < 3; ii++)
          _viscce[c_id][ii] = tmp*visten[c_id][ii] + visls_0;
        for (cs_lnum_t ii = 3; ii < 6; ii++)
          _viscce[c_id][ii] = tmp*visten[c_id][ii];
      }
    }
    else {
      const cs_real_t tmp = eqp->idifft*ctheta/cs_turb_csrij;
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        for (cs_lnum_t ii = 0; ii < 3; ii++)
          _viscce[c_id][ii] = tmp*visten[c_id][ii] + cpro_viscls[c_id];
        for (cs_lnum_t ii = 3; ii < 6; ii++)
          _viscce[c_id][ii] = tmp*visten[c_id][ii];
      }
    }

    /* Weighting for gradient */
    if (cpro_wgrec_v != nullptr) {
      cs_array_real_copy(6*n_cells,
                         (cs_real_t *)_viscce,
                         (cs_real_t *)cpro_wgrec_v);
      cs_halo_sync_r(m->halo, false, cpro_wgrec_v);
    }

    cs_face_anisotropic_viscosity_scalar(m,
                                         fvq,
                                         _viscce,
                                         eqp->verbosity,
                                         _weighf,
                                         _weighb,
                                         viscf,
                                         viscb);
  }

  *weighb = _weighb;
  *weighf = _weighf;
  *viscce = _viscce;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the source terms for scalars which are part of
 * specific physics models. Source terms are defined over one time step.
 *
 * Warning: source terms are treated differently from the way they are
 *          in cs_user_source_terms.
 * fimp*d(var) = smbrs is solved. fimp and rhs already hold possible user
 * source terms values and thus have to be incremented (and not overwritten).
 *
 * For stability reasons, only positive terms are added to fimp, while there
 * are no such constrains on values to be added to rhs.
 *
 * In the case of a source term of the form cexp + cimp*var, the source term
 * should be implemented as follows:
 * \f[
 *   smbrs  = smbrs  + cexp + cimp*var
 * \f]
 * \f[
 *   rovsdt = rovsdt + max(-cimp,0)
 * \f]
 *
 *rovsdt and smbrs are provided here respectively in kg/s and in kg/s*[scalar].
 * Examples:
 *   velocity \f$ kg m/s^2 \f$
 *   temperature \f$ kg K/s \f$
 *   enthalpy \f$ J/s \f$
 *
 * \param[in]       f          pointer to field structure
 * \param[in, out]  iterns     explicit source terms
 * \param[in, out]  fimp       implicit source terms
 */
/*----------------------------------------------------------------------------*/

static void
_physical_model_source_terms(cs_field_t        *f,
                             cs_real_t          rhs[],
                             cs_real_t          fimp[])
{
  const int *pm_flag = cs_glob_physical_model_flag;

  /* Gas combustion */

  if (cs_glob_combustion_gas_model != nullptr) {
    cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;
    int cm_type = cm->type / 100;
    if (cm_type == CS_COMBUSTION_SLFM)
      cs_steady_laminar_flamelet_source_terms(f, rhs, fimp);

    if (cm->isoot >= 1) {
      if (f == cm->fsm || f == cm->npm)
        cs_soot_production(f->id, rhs, fimp);
    }

    if (cm_type == CS_COMBUSTION_EBU)
      cs_f_ebutss(f->id, rhs, fimp);
    else if (cm_type == CS_COMBUSTION_LW)
      cs_f_lwctss(f->id, rhs, fimp);
  }

  /* Pulverized coal combustion */

  else if (pm_flag[CS_COMBUSTION_COAL] >= 0)
    cs_coal_source_terms_scalar(f, rhs, fimp);

  /* Atmospheric version */

  if (pm_flag[CS_ATMOSPHERIC] >= 0) {
    cs_atmo_scalar_source_term(f->id, rhs);

    cs_atmo_option_t *at_opt = cs_glob_atmo_option;
    if (at_opt->rain == true){
      cs_atmo_source_term(f->id, rhs, fimp);
    }
  }

  /*! Electric arcs, Joule effect ionic conduction */

  if (   pm_flag[CS_JOULE_EFFECT] > 0
      || pm_flag[CS_ELECTRIC_ARCS] > 0) {
    const cs_mesh_t *m = cs_glob_mesh;
    const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
    cs_elec_source_terms(m, fvq, f->id, rhs);
  }

  /*! Cooling towers */

  if (pm_flag[CS_COOLING_TOWERS] > 0)
    cs_ctwr_source_term(f->id, rhs, fimp);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve the convection/diffusion equation (with optional source
 *        terms and/or drift) for a scalar quantity over a time step.
 *
 * \param[in]     f          pointer to field structure
 * \param[in]     iterns     Navier-Stokes iteration number
 * \param[in]     itspdv     indicator to compute production/dissipation
 *                           terms for a variance:
 *                           - 0: no
 *                           - 1: yes
 * \param         viscf      visc*surface/dist at internal faces (work array)
 * \param         viscb      visc*surface/dist at boundary faces (work array)
 */
/*----------------------------------------------------------------------------*/

void
cs_solve_equation_scalar(cs_field_t        *f,
                         int                iterns,
                         int                itspdv,
                         cs_real_t          viscf[],
                         cs_real_t          viscb[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_mesh_quantities_t *mq_g = cs_glob_mesh_quantities_g;
  const cs_time_step_t *ts = cs_glob_time_step;
  const cs_fluid_properties_t *fluid_props = cs_glob_fluid_properties;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_real_t *volume = mq_g->cell_vol;
  const cs_real_t *cell_f_vol = fvq->cell_vol;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *)fvq->cell_cen;

  const cs_turb_model_t *turb_model = cs_glob_turb_model;
  const cs_wall_condensation_t *wall_condensation = cs_glob_wall_condensation;

  const cs_real_t *dt = CS_F_(dt)->val;

  /* Initialization
     -------------- */

  cs_real_t *cvar_var = f->val;
  cs_real_t *cvara_var = f->val_pre;

  const int keysca = cs_field_key_id("scalar_id");
  const int keyvar = cs_field_key_id("variable_id");

  const int ivar = cs_field_get_key_int(f, keyvar);
  const int iscal = cs_field_get_key_int(f, keysca);

  /* We might have a time step multiplier */

  cs_real_t *dtr = nullptr;
  {
    const int keycdt = cs_field_key_id("time_step_factor");
    const cs_real_t cdtvar = cs_field_get_key_double(f, keycdt);

    if (fabs(cdtvar - 1.0) > cs_math_epzero) {
      CS_MALLOC(dtr, n_cells_ext, cs_real_t);
#     pragma omp parallel for if(n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
        dtr[c_id] = dt[c_id] * cdtvar;
      dt = dtr;
    }
  }

  bool is_thermal_model_field = (f == cs_thermal_model_field());

  /* Physical quantities */

  const cs_real_t *crom = CS_F_(rho)->val;
  const cs_real_t *croma = CS_F_(rho)->val_pre;

  const int kromsl = cs_field_key_id_try("density_id");
  int icrom_scal = cs_field_get_key_int(f, kromsl);
  if (icrom_scal > 0) {
    crom = cs_field_by_id(icrom_scal)->val;
    croma = cs_field_by_id(icrom_scal)->val;
  }
  else {
    icrom_scal = CS_F_(rho)->id;
  }

  /* Key id for buoyant field (inside the Navier Stokes loop) */

  const int key_coupled_with_vel_p = cs_field_key_id_try("coupled_with_vel_p");
  const int coupled_with_vel_p_fld = cs_field_get_key_int(f, key_coupled_with_vel_p);

  if (   (coupled_with_vel_p_fld == 1 && iterns == -1)
      || (coupled_with_vel_p_fld == 0 && iterns != -1)) {
    return;
  }

  cs_equation_param_t *eqp = cs_field_get_equation_param(f);

  if (eqp->verbosity > 0)
    bft_printf(" ** SOLVING VARIABLE %s\n"
               "    ----------------\n\n", f->name);

  /* Source terms
     ============ */

  cs_real_t *fimp, *rhs;
  CS_MALLOC_HD(rhs, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(fimp, n_cells_ext, cs_real_t, cs_alloc_mode);
  //CS_MALLOC(fimp, n_cells_ext, cs_real_t);

  cs_array_real_fill_zero(n_cells, rhs);
  cs_array_real_fill_zero(n_cells, fimp);

  /* User source terms
   * ----------------- */

  /* Source terms for GUI */
  if (is_thermal_model_field == false)
    cs_gui_scalar_source_terms(f,
                               cvar_var,
                               rhs,
                               fimp);
  else
    cs_gui_thermal_source_terms(f, cvar_var, rhs, fimp);

  cs_user_source_terms(cs_glob_domain,
                       f->id,
                       rhs,
                       fimp);

  /* Coupling between multiple code_saturne instances */
  const int nbrcpl = cs_sat_coupling_n_couplings();
  if (nbrcpl > 0)
    cs_sat_coupling_exchange_at_cells(f, rhs, fimp);

  /* Store the source terms for convective limiter
     or time extrapolation for buoyant scalar */

  cs_real_t *cpro_scal_st = nullptr, *cproa_scal_st = nullptr;

  const int kst = cs_field_key_id_try("source_term_id");
  const int kstprv = cs_field_key_id_try("source_term_prev_id");
  const int st_id = cs_field_get_key_int(f, kst);
  const int st_prv_id = cs_field_get_key_int(f, kstprv);

  if (st_id > -1) {
    cpro_scal_st = cs_field_by_id(st_id)->val;
    /* Handle parallelism and periodicity */
    if (m->halo != nullptr)
      cs_halo_sync_var(m->halo, CS_HALO_STANDARD, cpro_scal_st);
  }

  if (st_prv_id > -1)
    cproa_scal_st = cs_field_by_id(st_prv_id)->val;

  if (   (eqp->ibdtso > 1) && (ts->nt_cur  > ts->nt_ini)
      && (   cs_glob_time_step_options->idtvar == CS_TIME_STEP_CONSTANT
          || cs_glob_time_step_options->idtvar ==  CS_TIME_STEP_ADAPTIVE))
    cs_backward_differentiation_in_time(f, rhs, fimp);

  /* Skip first time step after restart if previous values have not been read. */
  if (eqp->ibdtso < 0)
    eqp->ibdtso = - eqp->ibdtso;

  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 0) {

    /* Nudging towards optimal interpolation for current scalar */
    const int kopint = cs_field_key_id_try("opt_interp_id");
    const int f_oi_id =  cs_field_get_key_int(f, kopint);
    if (f_oi_id > -1)
      cs_at_data_assim_source_term(f->id, rhs, fimp);

    /* Atmospheric chemistry
     * In case of a semi-coupled resolution, computation of the explicit
     * chemical source term to be considered during dynamical resolution */

    const cs_atmo_chemistry_t *atmo_chem = cs_glob_atmo_chemistry;

    if (   atmo_chem->model > 0
        && atmo_chem->chemistry_sep_mode  == 2
        && (ts->nt_cur >= ts->nt_ini)) {
      const int *isca_chem = atmo_chem->species_to_scalar_id;
      const int nespg = atmo_chem->n_species;
      if ((isca_chem[0] <= iscal) && (iscal <= isca_chem[nespg-1]))
        cs_atmo_chem_source_terms(iscal, rhs, fimp);
    }

  }

  /* Precipitation/dissolution for Lagrangian module
   * Calculation of source terms due to precipitation and dissolution phenomena
   * FIXME: test on scalar id since original commit in 2017 is a really dumb
   * idea.  */
  if (cs_glob_lagr_model->precipitation == 1) {
    if (iscal == 1)
      cs_lagr_precipitation_mass_st(ts->dt_ref, crom, cvar_var, rhs);
    cs_assert(iscal == 1);  /* Force error in other cases to fix this */
  }

  /* If we extrapolate source terms:
   *   rhs receives -theta TS from previous time step.
   *   rhs receives the part of the source term which depends on the variable.
   *   At order 2, we assume that fimp provided by the user is < 0, so
   *     we implicit the term (so fimp*cvar_prev goes in rhs).
   *   In standard case, adapt to sign of fimp, but fimp*cvar_prev still
   *     goes in rhs (no other choice). */

   /* Map the source term pointer depending on whether it's a bouyant scalar */
  cs_real_t *cpro_st = nullptr;
  if (st_prv_id > -1) {
    if (iterns == -1)
      cpro_st = cproa_scal_st;
    else
      cpro_st = cpro_scal_st;
  }

  const cs_real_t thetv = eqp->theta;
  const int kthetss = cs_field_key_id_try("st_exp_extrapolated");
  const cs_real_t thets = cs_field_get_key_double(f, kthetss);

  if (st_prv_id > -1) {
#   pragma omp parallel for if(n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      const cs_real_t smbexp = cproa_scal_st[c_id];
      /* If the scalar is not buoyant no need of saving the current source term,
       * save directly the previous one */
      cpro_st[c_id] = rhs[c_id];
      /* Source term of the previous time step.
       * We assume -fimp > 0: we implicit the user source term (the rest) */
      rhs[c_id] = fimp[c_id]*cvara_var[c_id] - thets * smbexp;
      /* Diagonal */
      fimp[c_id] *= -thetv;
    }
  }

  /* If we do not extrapolate the ST: */
  else {
#   pragma omp parallel for if(n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      /* User source term */
      rhs[c_id] += fimp[c_id] * cvara_var[c_id];
      /* Diagonal */
      fimp[c_id] = cs_math_fmax(-fimp[c_id], 0.0);
    }
  }

  /* Add thermodynamic pressure variation for the low-Mach algorithm */
  /* NB: cs_thermal_model_field is the Enthalpy. */
  const int ipthrm = fluid_props->ipthrm;
  const int idilat = cs_glob_velocity_pressure_model->idilat;
  if ((idilat == 3 || ipthrm == 1) && is_thermal_model_field) {
    const cs_real_t pther = fluid_props->pther;
    const cs_real_t pthera = fluid_props->pthera;
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      rhs[c_id] += (pther - pthera) / dt[c_id]*cell_f_vol[c_id];
  }

  /* If solving the internal energy equation:
     Add to the RHS -pdiv(u) + tau:u
     In case of no specific EOS, skip this */

  const cs_thermal_model_t *th_model = cs_glob_thermal_model;
  const cs_cf_model_t *th_cf_model = cs_glob_cf_model;

  const cs_real_t *temp = nullptr, *tempa = nullptr, *cpro_yw = nullptr;
  cs_real_t *cpro_yv = nullptr, *xcvv = nullptr;
  if (   th_cf_model->ieos != CS_EOS_NONE
      && is_thermal_model_field
      && (   th_model->thermal_variable == CS_THERMAL_MODEL_TEMPERATURE
          || th_model->thermal_variable == CS_THERMAL_MODEL_INTERNAL_ENERGY)) {

    CS_MALLOC(xcvv, n_cells_ext, cs_real_t);

    /* Compute cv */

    cs_field_t *f_cv = cs_field_by_name_try("isobaric_heat_capacity");
    if (f_cv != nullptr) {
      cs_thermal_model_cv(f_cv->val);
      cs_array_real_copy(n_cells, f_cv->val, xcvv);
    }
    else
      cs_array_set_value_real(n_cells, 1, 1., xcvv);

    /* Note that the temperature is a postprocessing field when the thermal model
       is based on the internal energy. */

    if (th_model->thermal_variable == CS_THERMAL_MODEL_INTERNAL_ENERGY) {

      const cs_field_t *f_t = cs_field_by_name_try("temperature");

      if (f_t != nullptr) {
        temp  = f_t->val;
        tempa = f_t->val_pre;
      }

      if (th_cf_model->ieos == CS_EOS_MOIST_AIR) {
        const cs_field_t *f_yv = cs_field_by_name_try("yv");
        const cs_field_t *f_yw = cs_field_by_name_try("yw");

        if (f_yv != nullptr) {
          cpro_yv = f_yv->val;
        }

        if (f_yw != nullptr) {
          cpro_yw = f_yw->val;
        }
      }

    }
    else if (th_model->thermal_variable == CS_THERMAL_MODEL_TEMPERATURE) {
      temp  = cvar_var;
      tempa = cvara_var;
    }

    cs_real_t *vistot;
    cs_real_33_t *gradv;

    CS_MALLOC_HD(vistot, n_cells_ext, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(gradv, n_cells_ext, cs_real_33_t, cs_alloc_mode);

    /* Total viscosity */

    const cs_real_t *viscl = CS_F_(mu)->val;

    if (turb_model->itytur == 3)
      cs_array_real_copy(n_cells, viscl, vistot);
    else {
      const cs_real_t *visct = CS_F_(mu_t)->val;
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        vistot[c_id] = viscl[c_id] + visct[c_id];
    }

    cs_field_gradient_vector(CS_F_(vel), false, 1, gradv);

    /* Compute p.div(u) */

    /* Parallelism and periodicity: vel already synchronized;
       gradient (below) is always synchronized. */

    cs_halo_sync_var(m->halo, CS_HALO_STANDARD, vistot);

    /* Get pressure gradient, including the pressure increment gradient */

    cs_thermal_model_pdivu(rhs);

    cs_thermal_model_dissipation(vistot, gradv, rhs);

    CS_FREE_HD(vistot);
    CS_FREE_HD(gradv);
  }
  if (is_thermal_model_field) {

    /* Kinetic source term when needed */

    if (    iterns >= 1
        &&  cs_glob_physical_model_flag[CS_COMPRESSIBLE] < 0
        && (   th_model->thermal_variable == CS_THERMAL_MODEL_TEMPERATURE
            || th_model->thermal_variable == CS_THERMAL_MODEL_INTERNAL_ENERGY)) {

      cs_thermal_model_add_kst(rhs);
    }

    /* CFL related to the internal energy equation */

    if (th_model->thermal_variable == CS_THERMAL_MODEL_TEMPERATURE) {
        //|| th_model->thermal_variable == CS_THERMAL_MODEL_INTERNAL_ENERGY) { TODO
      cs_field_t *f_cflt = cs_field_by_name_try("cfl_t");

      if (f_cflt != nullptr) {
        cs_real_t *cflt = f_cflt->val;

        /* Only implemented for the ideal gas equation of state. */

        if (th_cf_model->ieos == CS_EOS_IDEAL_GAS) {
          const cs_real_3_t *vel = (const cs_real_3_t *)CS_F_(vel)->val;

          const int kimasf = cs_field_key_id_try("inner_mass_flux_id");
          const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
          const cs_real_t *i_massflux
            = cs_field_by_id(cs_field_get_key_int(f, kimasf) )->val;
          const cs_real_t *b_massflux
            = cs_field_by_id(cs_field_get_key_int(f, kbmasf) )->val;

          /* precaution for croma */
          cs_halo_sync_var(m->halo, CS_HALO_STANDARD, CS_F_(rho)->val_pre);

          cs_thermal_model_cflt(croma,
                                temp,
                                tempa,
                                xcvv,
                                vel,
                                i_massflux,
                                b_massflux,
                                cflt);
        }
      }

    }

  } /* if (thermal_model_field) */

  /* Volume coupling with Syrthes; order 2 not handled. */
  if (is_thermal_model_field)
    cs_syr_coupling_volume_source_terms(f->id, rhs, fimp);

  /* Specific physical models; order 2 not handled. */
  if (cs_glob_physical_model_flag[CS_PHYSICAL_MODEL_FLAG] > 0)
    _physical_model_source_terms(f, rhs, fimp);

  /*  Rayonnement
   *  2nd order not handled */
  cs_real_t *cpro_tsscal = nullptr;
  if (idilat > 3) {
    char fname[128];
    snprintf(fname, 128, "%s_dila_st", f->name); fname[127] = '\0';
    cpro_tsscal = cs_field_by_name_try(fname)->val;
  }

  if (cs_glob_rad_transfer_params->type > CS_RAD_TRANSFER_NONE) {

    if (is_thermal_model_field) {
      cs_rad_transfer_source_terms(rhs, fimp);

      /* Store the explicit radiative source term */
      if (idilat > 3) {
        const cs_real_t *cpro_tsre1 = cs_field_by_name("rad_st")->val;
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          cpro_tsscal[c_id] += cpro_tsre1[c_id]*cell_f_vol[c_id];
      }
    }

    /* Pulverized coal; order 2 not handled */
    if (cs_glob_coal_model != nullptr) {
      const cs_coal_model_t *cm = cs_glob_coal_model;
      const int nclacp = cm->nclacp;
      const int isca_ih21
        = cs_field_get_key_int(cs_field_by_id(cm->ih2[0]), keyvar);
      const int isca_ih2nl
        = cs_field_get_key_int(cs_field_by_id(cm->ih2[nclacp-1]), keyvar);
      if ((isca_ih21 <= ivar) && (ivar <= isca_ih2nl))
        cs_coal_rad_transfer_st(f, rhs, fimp);

      if (f == cs_field_by_name_try("x_c_h")) {
        const cs_real_t *cpro_tsre1 = cs_field_by_name("rad_st")->val;
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          rhs[c_id] += volume[c_id]*cpro_tsre1[c_id];

        for (int icla = 0; icla < nclacp; icla++) {
          char f_rad[64], f_xp[64];
          snprintf(f_xp, 64, "x_p_%02d", icla+1); f_xp[63] = '\0';
          snprintf(f_rad, 64, "rad_st_%02d", icla+2); f_rad[63] = '\0';
          const cs_real_t *cpro_tsre = cs_field_by_name_try(f_rad)->val;
          const cs_real_t *cpro_x2icla = cs_field_by_name_try(f_xp)->val;
          for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
            rhs[c_id] -= volume[c_id]*cpro_tsre[c_id]*cpro_x2icla[c_id];
        }
      }
    }

  }

  /* When solving the Temperature, we solve:
   *  rho*cp*Vol*dT/dt + ... */
  int imucpp = 0, iscacp;

  cs_real_t *xcpp = nullptr;
  CS_MALLOC_HD(xcpp, n_cells_ext, cs_real_t, cs_alloc_mode);

  const int kscavr = cs_field_key_id("first_moment_id");
  const int iscavr = cs_field_get_key_int(f, kscavr);
  const int kscacp = cs_field_key_id("is_temperature");
  if (iscavr > 0)
    iscacp = cs_field_get_key_int(cs_field_by_id(iscavr), kscacp);
  else
    iscacp = cs_field_get_key_int(f, kscacp);

  /* Multiplier for thermal fields (Cp/Cv) */

  if (iscacp == 1 || iscacp == 2)
    imucpp = 1;

  if (imucpp == 0) {
    cs_array_set_value_real(n_cells, 1, 1., xcpp);
  }
  else if (imucpp == 1) {
    if (CS_F_(cp) != nullptr)
      cs_array_real_copy(n_cells, CS_F_(cp)->val, xcpp);
    else
      cs_array_set_value_real(n_cells, 1, fluid_props->cp0, xcpp);
    if (iscacp == 2) {
      cs_real_t rair = fluid_props->r_pg_cnst;
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        xcpp[c_id] -= rair;
    }
  }
  cs_halo_sync_var(m->halo, CS_HALO_STANDARD, xcpp);

  /* Lagrangien (couplage retour thermique)
   * Ordre 2 non pris en compte */
  if ((cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING) &&
      (cs_glob_lagr_source_terms->ltsthe == 1)) {

    if (   is_thermal_model_field
        && (th_model->thermal_variable == CS_THERMAL_MODEL_TEMPERATURE)) {
      cs_real_t *ste = cs_field_by_name("lagr_st_temperature")->val;
      cs_real_t *sti = cs_field_by_name("lagr_st_imp_temperature")->val;
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        rhs[c_id] += ste[c_id] * cell_f_vol[c_id];
        fimp[c_id] += cs_math_fmax(sti[c_id], 0.0) * cell_f_vol[c_id];
      }
    }
    /* When solving in enthalpy, no clear way to implicit Lagrangian source
     * term */
    if (   is_thermal_model_field
        && th_model->thermal_variable == CS_THERMAL_MODEL_ENTHALPY) {
      cs_real_t *ste = cs_field_by_name("lagr_st_temperature")->val;
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        rhs[c_id] += ste[c_id] * cell_f_vol[c_id];
    }

  }

  /* Mass source term.

     TODO this section could be replaced by lower-level code from
     cs_volume_mass_injection.c for the specific field being handled
     (i.e. using the direct xdef-based evaluation).
     This would allow removing the mst_type and mst_val arrays and
     associated dimensions. */

  if (eqp->n_volume_mass_injections > 0) {
    cs_lnum_t n_elts = 0;
    const cs_lnum_t *elt_ids = nullptr;
    int *mst_type_sc = nullptr;
    cs_real_t *mst_val_sc = nullptr, *mst_val_p = nullptr;

    cs_volume_mass_injection_get_arrays(f,
                                        &n_elts,
                                        &elt_ids,
                                        &mst_type_sc,
                                        &mst_val_sc,
                                        &mst_val_p);

    cs_real_t *srcmas;
    CS_MALLOC(srcmas, n_elts, cs_real_t);

    /* When treating the Temperature, the equation is multiplied by Cp */
    for (cs_lnum_t c_idx = 0; c_idx < n_elts; c_idx++) {
      if (mst_val_p[c_idx] > 0.0) {
        const cs_lnum_t id = elt_ids[c_idx];
        srcmas[c_idx] = mst_val_p[c_idx] * xcpp[id];
      }
      else {
       srcmas[c_idx] = 0.;
      }
    }

    /* If we extrapolate STs we put Gamma Pinj in cpro_st;
       Otherwise, we put it directly in rhs */
    cs_real_t *gapinj = (st_prv_id >= 0) ? cpro_st : rhs;

    /* Increment rhs by -Gamma.cvara_var and fimp by Gamma */
    cs_mass_source_terms(1,
                         1,
                         n_elts,
                         elt_ids,
                         mst_type_sc,
                         cell_f_vol,
                         cvara_var,
                         mst_val_sc,
                         srcmas,
                         rhs,
                         fimp,
                         gapinj);

    CS_FREE(srcmas);
  }

  /* Condensation source terms for the scalars
   * associated to a surface zone (icondb=0) and
   * condensation source terms for the scalars
   * associated to a volumic zone (icondv=0)
   * taking into account the metal mass
   * structures condensation modelling */

  if (wall_condensation->icondb == 0 || wall_condensation->icondv == 0) {

    /* Wall condensation */

     cs_wall_condensation_source_terms(f,
                                       xcpp,
                                       cvara_var,
                                       rhs,
                                       fimp);
  }

  /* viscosity and diffusivity*/
  const cs_real_t *cpro_viscls = nullptr;
  const cs_real_t *viscl = CS_F_(mu)->val;
  const int kivisl = cs_field_key_id("diffusivity_id");
  const int kvisl0 = cs_field_key_id("diffusivity_ref");
  const int ifcvsl = cs_field_get_key_int(f, kivisl);
  const cs_real_t visls_0 = cs_field_get_key_double(f, kvisl0);
  if (ifcvsl > -1)
    cpro_viscls = cs_field_by_id(ifcvsl)->val;

  /* Initialize turbulent diffusivity for SGDH model */
  cs_real_t *sgdh_diff;
  CS_MALLOC_HD(sgdh_diff, n_cells_ext, cs_real_t, cs_alloc_mode);
  //CS_MALLOC(sgdh_diff, n_cells_ext, cs_real_t);

  _init_sgdh_diff(f, sgdh_diff);

  /* If the current scalar is the variance of an other scalar,
   * production and dissipation terms are added. */
  const int krvarfl = cs_field_key_id("variance_dissipation");
  const cs_real_t rvarfl = cs_field_get_key_double(f, krvarfl);
  if (itspdv == 1)
    _production_and_dissipation_terms(f,
                                      ifcvsl,
                                      st_prv_id,
                                      rvarfl,
                                      visls_0,
                                      thetv,
                                      dt,
                                      xcpp,
                                      crom,
                                      viscl,
                                      cvar_var,
                                      sgdh_diff,
                                      cpro_viscls,
                                      rhs,
                                      fimp,
                                      cpro_st,
                                      cpro_tsscal);

  if (st_prv_id > -1) {
    const cs_real_t thetp1 = 1 + thets;
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      rhs[c_id] += thetp1 * cpro_st[c_id];
  }

  /* Compressible algorithm
   * or Low Mach compressible algos with mass flux prediction */
  const cs_real_t *pcrom;
  if (  cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1
      || (   cs_glob_velocity_pressure_model->idilat > 1
          && cs_glob_velocity_pressure_param->ipredfl == 1
          && fluid_props->irovar == 1)) {
    pcrom = croma;
  }
  /* Low Mach compressible algos (conservative in time).
   * Same algo. for Volume of Fluid method */
  else if (   (   cs_glob_velocity_pressure_model->idilat > 1
               || cs_glob_vof_parameters->vof_model > 0)
           && fluid_props->irovar == 1
           && (icrom_scal == CS_F_(rho)->id)) {
    if (iterns == 1)
      pcrom = cs_field_by_id(icrom_scal)->vals[2];
    else
      pcrom = cs_field_by_id(icrom_scal)->val_pre;
  }
  /* Deprecated algo or constant density */
  else {
    pcrom = cs_field_by_id(icrom_scal)->val;
  }

  /* Face diffusion "speed"

   * We take max(mu_t, 0) as in dynamic LES mu_t can be negative
   * (clipping over (mu + mu_t)). We could have taken max(K + K_t, 0)
   * but this would allow negative K_t negative, which is considered
   * non-physical. */
  cs_real_t *weighb = nullptr;
  cs_real_2_t *weighf = nullptr;
  cs_real_6_t *viscce = nullptr;

  if (eqp->idiff >= 1) {
    _diffusion_terms_scalar(f,
                            eqp,
                            visls_0,
                            xcpp,
                            sgdh_diff,
                            cpro_viscls,
                            viscf,
                            viscb,
                            rhs,
                            &weighb,
                            &weighf,
                            &viscce);
  }
  else {
    cs_array_real_fill_zero(m->n_i_faces, viscf);
    cs_array_real_fill_zero(n_b_faces, viscb);
  }

  /* Add Rusanov fluxes */
  if (cs_glob_turb_rans_model->irijnu == 2) {
    cs_real_t *ipro_rusanov = cs_field_by_name("i_rusanov_diff")->val;
    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {
      viscf[face_id] += 0.5 * ipro_rusanov[face_id];
    }
    //TODO add boundary terms?
  }

  CS_FREE(sgdh_diff);

  if (eqp->istat == 1) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      fimp[c_id] += xcpp[c_id]*pcrom[c_id]*cell_f_vol[c_id]/dt[c_id];
  }

  /* Scalar with a Drift:
   * compute the convective flux
   *---------------------------- */

  const int keydri = cs_field_key_id_try("drift_scalar_model");
  const int iscdri = cs_field_get_key_int(f, keydri);

  /* Interior mass flux */
  const int kimasf = cs_field_key_id_try("inner_mass_flux_id");
  const int iflmas = cs_field_get_key_int(f, kimasf);
  cs_real_t *imasfl = cs_field_by_id(iflmas)->val;

  /* Boundary mass flux */
  const int kbmasf = cs_field_key_id_try("boundary_mass_flux_id");
  const int iflmab = cs_field_get_key_int(f, kbmasf);
  cs_real_t *bmasfl = cs_field_by_id(iflmab)->val;

  if (iscdri > 0)
    cs_drift_convective_flux(f,
                             imasfl,
                             bmasfl,
                             fimp,
                             rhs);

  /* Cancel RHS in disabled cells in case spurious terms were added
     by "generic" code */

  if (fvq->has_disable_flag) {
    cs_dispatch_context ctx;
    int *c_disable_flag = fvq->c_disable_flag;
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      if (c_disable_flag[c_id])
        rhs[c_id] = 0;
    });
    ctx.wait();
  }

  /* Solve
   * ===== */

  cs_real_t *cvark_var = nullptr, *wcvark_var = nullptr;
  if (iterns >= 1) {
    CS_MALLOC_HD(wcvark_var, n_cells_ext, cs_real_t, cs_alloc_mode);
    cvark_var = wcvark_var;
    cs_array_real_copy(n_cells_ext, cvar_var, cvark_var);
  }
  else {
    cvark_var = f->val;
  }

  int iescap = 0, icvflb = 0;


  /* all boundary convective flux with upwind */
  cs_real_t normp = -1.0;

  cs_real_t *dpvar;
  CS_MALLOC_HD(dpvar, n_cells_ext, cs_real_t, cs_alloc_mode);

  cs_equation_iterative_solve_scalar(cs_glob_time_step_options->idtvar,
                                     iterns,
                                     f->id,
                                     nullptr,
                                     iescap,
                                     imucpp,
                                     normp,
                                     eqp,
                                     cvara_var,
                                     cvark_var,
                                     f->bc_coeffs,
                                     imasfl,
                                     bmasfl,
                                     viscf,
                                     viscb,
                                     viscf,
                                     viscb,
                                     viscce,
                                     weighf,
                                     weighb,
                                     icvflb,
                                     nullptr,
                                     fimp,
                                     rhs,
                                     cvar_var,
                                     dpvar,
                                     xcpp,
                                     nullptr);

  CS_FREE_HD(dpvar);

  /* May be allocated in _diffusion_terms_scalar */
  CS_FREE_HD(weighb);
  CS_FREE_HD(weighf);
  CS_FREE_HD(viscce);

  /* When solving internal energy, compute the temperature */

  if (   is_thermal_model_field
      && th_model->thermal_variable == CS_THERMAL_MODEL_INTERNAL_ENERGY) {

    cs_real_t *tempk = cs_field_by_name("temperature")->val;

    /* Perfect gas, compute temperature from the internal energy */

    if (th_cf_model->ieos == CS_EOS_IDEAL_GAS) {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        tempk[c_id] = cvar_var[c_id] / xcvv[c_id];
    }

    /* Humid air module */

    else if (th_cf_model->ieos == CS_EOS_MOIST_AIR) {
      /* TODO Other Antoine law coefficients */

      /* Case of no saturation in previous iteration;
         Newton method
         First argument is the method used */

      /* precaution for cvar_pr */
      cs_halo_sync_var(m->halo, CS_HALO_STANDARD, CS_F_(p)->val);

      cs_thermal_model_newton_t(1,
                                nullptr,
                                cvar_var,
                                CS_F_(p)->val,
                                CS_F_(p)->val_pre,
                                cpro_yw,
                                cpro_yv,
                                tempk);

    }

  }

  CS_FREE(xcvv);

  /* Clipping, finalization of the computation of some terms and log
   * =============================================================== */

  cs_scalar_clipping(f);

  /* Atmospheric module:
   * finalize number of droplets due to nucleation for humid atmosphere */

  if (   cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_HUMID
      && cs_glob_atmo_option->sedimentation_model == 1
      && strcmp(f->name, "number_of_droplets") == 0
      && cs_glob_atmo_option->nucleation_model > 0) {

    const cs_real_t *cpro_liqwt = cs_field_by_name("liquid_water")->val;
    /* Test minimum liquid water to carry out nucleation */
    cs_real_t qliqmax = 0;
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      qliqmax = cs_math_fmax(cpro_liqwt[c_id], qliqmax);

    cs_parall_max(1, CS_REAL_TYPE, &qliqmax);

    if (qliqmax > 1e-8) {
      /* First : diagnose the droplet number
       * nucleation : when liquid water present calculate the
       * number of condensation nucleii (ncc) and if the droplet number (nc)
       * is smaller than ncc set it to ncc. */
      cs_real_t *pphy;
      CS_MALLOC(pphy, n_cells_ext, cs_real_t);

      if (cs_glob_atmo_option->meteo_profile == 0) {
        cs_real_t _pphy, dum;
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_atmo_profile_std(0., /* z_ref */
                              fluid_props->p0,
                              fluid_props->t0,
                              cell_cen[c_id][2], &_pphy, &dum, &dum);
          pphy[c_id] = _pphy;
        }
      }
      /* calculate pressure from meteo file */
      else if (cs_glob_atmo_option->meteo_profile == 1) {
        int nbmett = cs_glob_atmo_option->met_1d_nlevels_t;
        int nbmetm = cs_glob_atmo_option->met_1d_ntimes;
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          pphy[c_id] = cs_intprf(nbmett,
                                 nbmetm,
                                 cs_glob_atmo_option->z_temp_met,
                                 cs_glob_atmo_option->time_met,
                                 cs_glob_atmo_option->hyd_p_met,
                                 cell_cen[c_id][2],
                                 ts->t_cur);
      }
      else {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          pphy[c_id] = cs_field_by_name("meteo_pressure")->val[c_id];
      }

      cs_real_t *cvar_ntdrp = cs_field_by_name_try("number_of_droplets")->val;
      const cs_real_t *cpro_rad_cool = cs_field_by_name("radiative_cooling")->val;

      cs_atmo_aerosol_nuclea(cvar_ntdrp,
                             crom,
                             cpro_liqwt,
                             pphy,
                             cpro_rad_cool);
      CS_FREE(pphy);

    } // qliqmax.gt.1.e-8
  } // for humid atmosphere physics only

  if ((idilat > 3) && (itspdv == 1)) {

    cs_real_t xe = 0, xk = 0;

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      if (   turb_model->itytur == 2
          || turb_model->itytur == 5) {
        xk = CS_F_(k)->val_pre[c_id];
        xe = CS_F_(eps)->val_pre[c_id];
      }
      else if (turb_model->itytur == 3) {
        const cs_real_6_t *cvara_rij = (const cs_real_6_t*)CS_F_(rij)->val_pre;
        xe = CS_F_(eps)->val_pre[c_id];
        xk = 0.5*(cvara_rij[c_id][0] + cvara_rij[c_id][1] + cvara_rij[c_id][2]);
      }
      else if (turb_model->model == CS_TURB_K_OMEGA) {
        xk = CS_F_(k)->val_pre[c_id];
        xe = cs_turb_cmu*CS_F_(omg)->val_pre[c_id];
      }
      const cs_real_t rhovst =   xcpp[c_id] * crom[c_id] * xe / (xk * rvarfl)
                               * cell_f_vol[c_id];

      cpro_tsscal[c_id] -= rhovst * cvar_var[c_id];

    }
  }

  /* Store the implicit part of the radiative source term */
  if (   idilat > 3
      && cs_glob_rad_transfer_params->type > CS_RAD_TRANSFER_NONE
      && is_thermal_model_field) {
    const cs_real_t *cpro_tsri1 = cs_field_by_name("rad_st_implicit")->val;
#   pragma omp parallel for if(n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      const cs_real_t  dvar = cvar_var[c_id]-cvara_var[c_id];
      cpro_tsscal[c_id] -= cpro_tsri1[c_id]*dvar*cell_f_vol[c_id];
    }
  }

  /* Explicit balance
   * (See cs_equation_iterative_solve_scalar: remove the increment)
   * This should be valid with thr theta scheme on source terms. */
  if (eqp->verbosity > 1) {
    int ibcl = 0;
    if (eqp->nswrsm > 1)
      ibcl = 1;
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      rhs[c_id] -= eqp->istat*xcpp[c_id]*(pcrom[c_id]/dt[c_id])
        *cell_f_vol[c_id]*(cvar_var[c_id]-cvara_var[c_id])*ibcl;
    const cs_real_t sclnor = sqrt(cs_gdot(n_cells,rhs,rhs));

    bft_printf("%s: EXPLICIT BALANCE = %14.5e\n\n",f->name, sclnor);
  }

  /* Log in case of velocity/pressure inner iterations */
  if ((iterns > 0) && (eqp->verbosity > 1)) {
    cs_real_t *errork;
    CS_MALLOC(errork, n_cells_ext, cs_real_t);
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      errork[c_id] = cvar_var[c_id] - cvark_var[c_id];

    double l2errork = sqrt(cs_gres(n_cells, cell_f_vol, errork, errork));
    CS_FREE(errork);

    double l2norm = sqrt(cs_gres(n_cells, cell_f_vol, cvara_var, cvara_var));
    double dl2norm = 1.0;
    if (l2norm > 0.)
      dl2norm = 1.0/l2norm;

    bft_printf("Inner iteration scalar %10d iter=%10d L2 error = %12.4e\n"
               "L2 normalized error %e12.4 nomr %e12.4\n",f->id, iterns,
               l2errork, l2errork*dl2norm, l2norm);
  }

  CS_FREE(dtr);

  CS_FREE_HD(fimp);
  CS_FREE_HD(rhs);
  CS_FREE_HD(wcvark_var);
  CS_FREE_HD(xcpp);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve the convection/diffusion equation (with optional source
 *        terms and/or drift) for a vectorial quantity over a time step..
 *
 * \param[in]     f          pointer to field structure
 * \param[in]     iterns     Navier-Stokes iteration number
 * \param         viscf      visc*surface/dist at internal faces (work array)
 * \param         viscb      visc*surface/dist at boundary faces (work array)
 */
/*----------------------------------------------------------------------------*/

void
cs_solve_equation_vector(cs_field_t       *f,
                         int               iterns,
                         cs_real_t         viscf[],
                         cs_real_t         viscb[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_fluid_properties_t *fluid_props = cs_glob_fluid_properties;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_real_t *cell_f_vol = fvq->cell_vol;

  const cs_real_t *dt = CS_F_(dt)->val;

  /* Initialization
     -------------- */

  cs_real_3_t *cvar_var = (cs_real_3_t *)f->val;
  cs_real_3_t *cvara_var = (cs_real_3_t *)f->val_pre;

  /* If the vector is buoyant, it is inside the Navier Stokes loop, and
     so iterns >=1; otherwise it is outside of the loop and iterns = -1. */

  const int key_coupled_with_vel_p = cs_field_key_id_try("coupled_with_vel_p");
  const int coupled_with_vel_p_fld = cs_field_get_key_int(f, key_coupled_with_vel_p);

  if (   (coupled_with_vel_p_fld == 1 && iterns == -1)
      || (coupled_with_vel_p_fld == 0 && iterns != -1)) {
    return;
  }

  /* We might have a time step multiplier */

  cs_real_t *dtr = nullptr;
  {
    const int keycdt = cs_field_key_id("time_step_factor");
    const cs_real_t cdtvar = cs_field_get_key_double(f, keycdt);

    if (fabs(cdtvar - 1.0) > cs_math_epzero) {
      CS_MALLOC(dtr, n_cells_ext, cs_real_t);
#     pragma omp parallel for if(n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
        dtr[c_id] = dt[c_id] * cdtvar;
      dt = dtr;
    }
  }

  /* Physical quantities */

  const cs_real_t *croma = CS_F_(rho)->val;
  if (CS_F_(rho)->val_pre != nullptr)
    croma = CS_F_(rho)->val_pre;

  cs_equation_param_t *eqp = cs_field_get_equation_param(f);

  if (eqp->verbosity > 0)
    bft_printf(" ** SOLVING VARIABLE %s\n"
               "    ----------------\n\n", f->name);

  /* Source terms
     ============ */

  cs_real_3_t *rhs;
  cs_real_33_t *fimp;

  CS_MALLOC(fimp, n_cells_ext, cs_real_33_t);
  CS_MALLOC(rhs, n_cells_ext, cs_real_3_t);

  cs_array_real_fill_zero(9*n_cells_ext, (cs_real_t *)fimp);
  cs_array_real_fill_zero(3*n_cells_ext, (cs_real_t *)rhs);

  /* User source terms
   * ----------------- */

  cs_user_source_terms(cs_glob_domain,
                       f->id,
                       (cs_real_t *)rhs,
                       (cs_real_t *)fimp);

  /* Coupling between multiple code_saturne instances */
  const int nbrcpl = cs_sat_coupling_n_couplings();
  if (nbrcpl > 0)
    cs_sat_coupling_exchange_at_cells(f,
                                      (cs_real_t *)rhs,
                                      (cs_real_t *)fimp);

  /* Store the source terms for convective limiter */
  cs_real_3_t *cpro_vect_st = nullptr;

  const int kst = cs_field_key_id_try("source_term_id");
  const int kstprv = cs_field_key_id_try("source_term_prev_id");
  const int st_id = cs_field_get_key_int(f, kst);
  const int st_prv_id = cs_field_get_key_int(f, kstprv);

  if (st_id > -1) {
    cs_field_t *f_st = cs_field_by_id(st_id);
    if (f_st->dim != 3)
      bft_error(__FILE__, __LINE__, 0,
                _("%s: source term field %s\n"
                  "for field %s should be dimension 3."),
                __func__, f_st->name, f->name);
    cpro_vect_st = (cs_real_3_t *)f_st->val;
    cs_array_real_copy(3*n_cells,
                       (const cs_real_t *)rhs,
                       (cs_real_t *)cpro_vect_st);
    /* Handle parallelism and periodicity */
    bool on_device = cs_mem_is_device_ptr(cpro_vect_st);
    cs_halo_sync_r(m->halo, on_device, cpro_vect_st);
  }

  /* If we extrapolate source terms:
   *   rhs receives -theta TS from previous time step.
   *   rhs receives the part of the source term which depends on the variable.
   *   At order 2, we assume that fimp provided by the user is < 0, so
   *     we implicit the term (so fimp*cvar_prev goes in rhs).
   *   In standard case, adapt to sign of fimp, but fimp*cvar_prev still
   *     goes in rhs (no other choice). */

  cs_real_3_t *cproa_vect_st = nullptr;
  const cs_real_t thetv = eqp->theta;
  const int kthetss = cs_field_key_id_try("st_exp_extrapolated");
  const cs_real_t thets = cs_field_get_key_double(f, kthetss);

  if (st_prv_id > -1) {
    cproa_vect_st = (cs_real_3_t *)cs_field_by_id(st_prv_id)->val;
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      cs_real_t smbexp[3];
      for (cs_lnum_t j = 0; j < 3; ++j) {
        smbexp[j] = cproa_vect_st[c_id][j];
        /* User explicit source terms */
        cproa_vect_st[c_id][j] = rhs[c_id][j];
        rhs[c_id][j] =   fimp[c_id][j][j] * cvara_var[c_id][j]
                         - thets*smbexp[j];
        /* Diagonal */
        fimp[c_id][j][j] = -thetv*fimp[c_id][j][j];
      }
    }
  }

  /* If we do not extrapolate the ST: */
  else {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        /* User explicit source terms */
        rhs[c_id][j] += fimp[c_id][j][j] * cvara_var[c_id][j];
        /* Diagonal */
        fimp[c_id][j][j] = cs_math_fmax(-fimp[c_id][j][j], 0);
      }
    }
  }

  /* Specific physical models; order 2 not handled */

  if (cs_glob_physical_model_flag[CS_PHYSICAL_MODEL_FLAG] > 0) {
    cs_atmo_option_t *at_opt = cs_glob_atmo_option;
    if (at_opt->rain == true){
      cs_atmo_source_term(f->id, (cs_real_t *)rhs, (cs_real_t *)fimp);
    }
    if (cs_glob_physical_model_flag[CS_COOLING_TOWERS] > 0)
      cs_ctwr_source_term(f->id, (cs_real_t *)rhs, (cs_real_t *)fimp);

    if (   cs_glob_physical_model_flag[CS_JOULE_EFFECT] > 0
        || cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] > 0)
      cs_elec_source_terms_v(m, fvq, f->id, rhs);

  }

  /* Mass source term */

  if (eqp->n_volume_mass_injections > 0) {
    cs_lnum_t n_elts = 0;
    const cs_lnum_t *elt_ids = nullptr;
    int *mst_type_v = nullptr;
    cs_real_t *mst_val_v = nullptr, *mst_val_p = nullptr;

    cs_volume_mass_injection_get_arrays(f,
                                        &n_elts,
                                        &elt_ids,
                                        &mst_type_v,
                                        &mst_val_v,
                                        &mst_val_p);

    /* If we extrapolate source terms we put -Gamma Pinj in cproa_vect_st;
       if we do not extrapolate we put directly in SMBRV */
    cs_real_3_t *gavinj = (st_prv_id > -1) ? cproa_vect_st : rhs;

    /* We increment rhs by -gamma vara and fimp by gamma */
    cs_mass_source_terms(1,
                         3,
                         n_elts,
                         elt_ids,
                         mst_type_v,
                         cell_f_vol,
                         (const cs_real_t *)cvara_var,
                         mst_val_v,
                         mst_val_p,
                         (cs_real_t *)rhs,
                         (cs_real_t *)fimp,
                         (cs_real_t *)gavinj);
  }

  if (st_prv_id > -1) {
    const cs_real_t thetp1 = 1.0 + thets;
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        rhs[c_id][j] += thetp1 * cproa_vect_st[c_id][j];
    }
  }

  /* Compressible algorithm
   * or Low Mach compressible algos with mass flux prediction */

  const cs_real_t *pcrom = nullptr;

  if (   cs_glob_physical_model_flag[CS_COMPRESSIBLE] > -1
      || (   cs_glob_velocity_pressure_model->idilat > 1
          && cs_glob_velocity_pressure_param->ipredfl == 1
          && fluid_props->irovar == 1)) {
    pcrom = croma;
  }

  /* Low Mach compressible algos (conservative in time).
   * Same algo. for Volume of Fluid method */
  else if (   (   cs_glob_velocity_pressure_model->idilat > 1
               || cs_glob_vof_parameters->vof_model > 0)
           && fluid_props->irovar == 1) {
    if (iterns == 1)
      pcrom = CS_F_(rho)->vals[2];
    else
      pcrom = CS_F_(rho)->val_pre;
  }
  /* Deprecated algo or constant density */
  else {
    pcrom = CS_F_(rho)->val;
  }

  /* Face diffusion "speed"

   * We take max(mu_t, 0) as in dynamic LES mu_t can be negative
   * (clipping over (mu + mu_t)). We could have taken max(K + K_t, 0)
   * but this would allow negative K_t negatif, which is considered
   * non-physical. */

  cs_real_t *weighb = nullptr;
  cs_real_2_t *weighf = nullptr;
  cs_real_6_t *viscce = nullptr;

  if (eqp->idiff >= 1) {
    _diffusion_terms_vector(f,
                            eqp,
                            viscf,
                            viscb,
                            &weighb,
                            &weighf,
                            &viscce);
  }
  else {
    cs_array_real_fill_zero(m->n_i_faces, viscf);
    cs_array_real_fill_zero(n_b_faces, viscb);
  }

  /* Add Rusanov fluxes */
  if (cs_glob_turb_rans_model->irijnu == 2) {
    cs_real_t *ipro_rusanov = cs_field_by_name("i_rusanov_diff")->val;
    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {
      viscf[face_id] += 0.5 * ipro_rusanov[face_id];
    }

    const cs_real_3_t *restrict b_face_normal
      = (const cs_real_3_t *)fvq->b_face_normal;
    cs_real_t *bpro_rusanov = cs_field_by_name("b_rusanov_diff")->val;
    // cs_real_3_t  *coefap = (cs_real_3_t *)f->bc_coeffs->a;
    // cs_real_33_t *coefbp = (cs_real_33_t *)f->bc_coeffs->b;
    // cs_real_3_t  *cofafp = (cs_real_3_t *)f->bc_coeffs->af;
    cs_real_33_t *cofbfp = (cs_real_33_t *)f->bc_coeffs->bf;

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      cs_real_t n[3];
      cs_math_3_normalize(b_face_normal[face_id], n);

      for (cs_lnum_t i = 0; i < 3; i++) {
        for (cs_lnum_t j = 0; j < 3; j++) {
          cs_real_t bij = bpro_rusanov[face_id] * (n[i]*n[j]);
          cofbfp[face_id][i][j] +=  bij;
          //TODO? cofafp[face_id][i] -= bij * coefap[face_id][j];
        }
      }
    }
  }

  if (eqp->istat == 1) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        fimp[c_id][j][j] += pcrom[c_id] * cell_f_vol[c_id] / dt[c_id];
    }
  }

  /* Scalar with a Drift: compute the convective flux
   * ------------------------------------------------ */

  const int keydri = cs_field_key_id_try("drift_scalar_model");
  const int iscdri = cs_field_get_key_int(f, keydri);

  /* Interior mass flux */
  const int kimasf = cs_field_key_id_try("inner_mass_flux_id");
  const int iflmas = cs_field_get_key_int(f, kimasf);
  cs_real_t *imasfl = cs_field_by_id(iflmas)->val;

  /* Boundary mass flux */
  const int kbmasf = cs_field_key_id_try("boundary_mass_flux_id");
  const int iflmab = cs_field_get_key_int(f, kbmasf);
  cs_real_t *bmasfl = cs_field_by_id(iflmab)->val;

  if (iscdri > 0)
    cs_drift_convective_flux(f,
                             imasfl,
                             bmasfl,
                             (cs_real_t *) fimp,
                             (cs_real_t *) rhs);

  /* Cancel RHS in disabled cells in case spurious terms were added
     by "generic" code */

  if (fvq->has_disable_flag) {
    cs_dispatch_context ctx;
    int *c_disable_flag = fvq->c_disable_flag;
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      if (c_disable_flag[c_id]) {
        for (cs_lnum_t j = 0; j < 3; j++)
          rhs[c_id][j] = 0;
      }
    });
    ctx.wait();
  }

  /* Solve
     ===== */

  int iescap = 0, icvflb = 0, ivissv = 0;

  cs_equation_iterative_solve_vector(cs_glob_time_step_options->idtvar,
                                     iterns,
                                     f->id,
                                     nullptr,
                                     ivissv,
                                     iescap,
                                     eqp,
                                     cvara_var,
                                     cvara_var,
                                     f->bc_coeffs,
                                     imasfl,
                                     bmasfl,
                                     viscf,
                                     viscb,
                                     viscf,
                                     viscb,
                                     nullptr,
                                     nullptr,
                                     viscce,
                                     weighf,
                                     weighb,
                                     icvflb,
                                     nullptr,
                                     fimp,
                                     rhs,
                                     cvar_var,
                                     nullptr);

  /* May be allocated in _diffusion_terms_scalar */
  CS_FREE(weighb);
  CS_FREE(weighf);
  CS_FREE(viscce);

  CS_FREE(fimp);

  if (eqp->verbosity > 1) {
    cs_real_t ibcl = 0;
    if (eqp->nswrsm > 1)
      ibcl = 1;

    if (eqp->istat == 1) {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        for (cs_lnum_t j = 0; j < 3; j++) {
          rhs[c_id][j] -=   (pcrom[c_id] / dt[c_id]) * cell_f_vol[c_id]
                            * (cvar_var[c_id][j] - cvara_var[c_id][j]) * ibcl;
        }
      }
    }
    const cs_real_t sclnor = sqrt(cs_gdot(3*n_cells,
                                          (const cs_real_t *)rhs,
                                          (const cs_real_t *)rhs));

    bft_printf("%s: EXPLICIT BALANCE = %14.5e\n\n",f->name, sclnor);
  }

  CS_FREE(rhs);
  CS_FREE(dtr);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
