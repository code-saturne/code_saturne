/*============================================================================
 * Solve the convection diffusion equation.
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

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_array.h"
#include "cs_assert.h"
#include "cs_at_data_assim.h"
#include "cs_atmo.h"
#include "cs_atmo_aerosol.h"
#include "cs_atmo_profile_std.h"
#include "cs_blas.h"
#include "cs_boundary_conditions.h"
#include "cs_bw_time_diff.h"
#include "cs_cf_model.h"
#include "cs_coal.h"
#include "cs_combustion_model.h"
#include "cs_ctwr.h"
#include "cs_divergence.h"
#include "cs_drift_convective_flux.h"
#include "cs_elec_model.h"
#include "cs_equation_iterative_solve.h"
#include "cs_face_viscosity.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_gradient.h"
#include "cs_gui.h"
#include "cs_intprf.h"
#include "cs_lagr.h"
#include "cs_lagr_precipitation_model.h"
#include "cs_math.h"
#include "cs_mass_source_terms.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_param_types.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_prototypes.h"
#include "cs_rad_transfer.h"
#include "cs_rad_transfer_source_terms.h"
#include "cs_sat_coupling.h"
#include "cs_scalar_clipping.h"
#include "cs_syr_coupling.h"
#include "cs_thermal_model.h"
#include "cs_time_step.h"
#include "cs_turbulence_model.h"
#include "cs_turbulence_rij.h"
#include "cs_turbulence_rit.h"
#include "cs_velocity_pressure.h"
#include "cs_vof.h"
#include "cs_wall_condensation.h"
#include "cs_wall_functions.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_solve_equation.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_convection_diffusion.c

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

extern int cs_glob_darcy_unsteady;

void
cs_f_get_dimens(int  *nvar,
                int  *nscal);

/* Use legacy macro type to maintain compatibility with legacy user files */

void
CS_PROCF(ustssc, USTSSC)(const int        *nvar,
                         const int        *nscal,
                         const cs_lnum_t  *ncepdp,
                         const cs_lnum_t  *ncesmp,
                         const int        *iscal,
                         const cs_lnum_t  *icepdc,
                         const cs_lnum_t  *icetsm,
                         cs_lnum_t        *itypsm,
                         const cs_real_t   dt[],
                         const cs_real_t  *ckupdc,
                         const cs_real_t  *smacel,
                         cs_real_t         crvexp[],
                         cs_real_t         crvimp[]);

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

void
cs_f_solve_equation_scalar(int              f_id,
                           cs_lnum_t        ncesmp,
                           int              iterns,
                           int              itspdv,
                           const cs_lnum_t  icetsm[],
                           int              itypsm[],
                           cs_real_t        smacel[],
                           cs_real_t        viscf[],
                           cs_real_t        viscb[]);

void
cs_f_solve_equation_vector(int              f_id,
                           cs_lnum_t        ncesmp,
                           int              iterns,
                           const cs_lnum_t  icetsm[],
                           int              itypsm[],
                           cs_real_t        smacel[],
                           cs_real_t        viscf[],
                           cs_real_t        viscb[]);

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
  cs_real_t *cpro_turb_schmidt = NULL;
  const int key_turb_schmidt = cs_field_key_id("turbulent_schmidt_id");
  const int t_scd_id = cs_field_get_key_int(f, key_turb_schmidt);
  if (t_scd_id > -1)
    cpro_turb_schmidt = cs_field_by_id(t_scd_id)->val;

  /* Retrieve turbulent diffusivity value for current scalar */
  cs_real_t *cpro_turb_diff = NULL;
  const int key_turb_diff = cs_field_key_id("turbulent_diffusivity_id");
  const int t_dif_id = cs_field_get_key_int(f, key_turb_diff);

  /* Variable turbulent diffusivity */
  if (t_dif_id > -1) {
    cpro_turb_diff = cs_field_by_id(t_dif_id)->val;
    cs_array_real_copy(cs_glob_mesh->n_cells, cpro_turb_diff, sgdh_diff);
  }
  /* Variable Schmidt number */
  else if (cpro_turb_schmidt != NULL)
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
         || turb_model->iturb == CS_TURB_K_OMEGA))
    return;

  const cs_field_t *f_fm = NULL;  /* First moment field of
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

  const cs_real_t *cell_f_vol = cs_glob_mesh_quantities->cell_f_vol;

  const cs_real_t *cvara_var = f->val_pre;
  cs_real_t *cvara_var_fm = f_fm->val_pre;

  /* Get the turbulent flux model for the scalar and its variance */

  const int kturt = cs_field_key_id("turbulent_flux_model");

  int variance_turb_flux_model = 0, variance_turb_flux_model_type = 0;

  if (f_fm != NULL) {
    variance_turb_flux_model = cs_field_get_key_int(f_fm, kturt);
    variance_turb_flux_model_type = variance_turb_flux_model / 10;
  }

  /* Homogeneous Neumann on convective inlet on the production term
     for the variance */

  const cs_real_t *coefap = f_fm->bc_coeffs->a;
  const cs_real_t *coefbp = f_fm->bc_coeffs->b;

  cs_real_t *coefa_p, *coefb_p;
  BFT_MALLOC(coefa_p, n_b_faces, cs_real_t);
  BFT_MALLOC(coefb_p, n_b_faces, cs_real_t);

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    coefa_p[face_id] = coefap[face_id];
    coefb_p[face_id] = coefbp[face_id];
    if (cs_glob_bc_type[face_id] == CS_CONVECTIVE_INLET) {
      coefa_p[face_id] = 0.;
      coefb_p[face_id] = 1.;
    }
  }

  cs_real_3_t *grad;
  BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);

  cs_field_gradient_scalar_array(f_fm->id,
                                 1,     /* inc */
                                 coefa_p,
                                 coefb_p,
                                 cvara_var_fm,
                                 grad);

  BFT_FREE(coefa_p);
  BFT_FREE(coefb_p);

  /* Production Term
     --------------- */

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
        cpro_st[c_id] -=   2 * xcpp[c_id] * cell_f_vol[c_id] * crom[c_id]
                         * cs_math_3_dot_product(grad[c_id], xut[c_id]);
      }
    }
    /* SGDH model */
    else {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cpro_st[c_id] +=   2 * xcpp[c_id] * sgdh_diff[c_id] * cell_f_vol[c_id]
                         * cs_math_3_dot_product(grad[c_id], grad[c_id]);
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
        rhs[c_id] += cs_math_fmax(prod*cprovol, 0);

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
        rhs[c_id] +=   2 * xcpp[c_id] * sgdh_diff[c_id] * cell_f_vol[c_id]
                       * cs_math_3_dot_product(grad[c_id], grad[c_id]);
      }
    }

    /* Production term for a variance  TODO compute ustdy when isso2t >0 */
    if (cs_glob_velocity_pressure_model->idilat >= 4)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cpro_tsscal[c_id] +=   2*xcpp[c_id] * sgdh_diff[c_id] * cell_f_vol[c_id]
                             * cs_math_3_dot_product(grad[c_id], grad[c_id]);
      }
  }

  BFT_FREE(grad);

  /* Dissipation term
     ---------------- */

  cs_real_6_t *cvara_rij = NULL;
  cs_real_t *cvara_k= NULL, *cvara_ep = NULL, *cvara_omg= NULL, *cvar_al = NULL;

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
  else if (turb_model->iturb == CS_TURB_K_OMEGA) {
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
    else if (turb_model->iturb == CS_TURB_K_OMEGA) {
      xk = cvara_k[c_id];
      xe = cs_turb_cmu*xk*cvara_omg[c_id];
    }

    if (cvar_al != NULL)
      alpha_theta = cvar_al[c_id];

    cs_real_t prdtl= viscl[c_id] * xcpp[c_id];
    if (ifcvsl > -1)
      prdtl /= cpro_viscls[c_id];
    else
      prdtl /= visls_0;

    const cs_real_t xr = (1.0 - alpha_theta)*prdtl + alpha_theta*rvarfl;
    const cs_real_t rhovst = xcpp[c_id]*crom[c_id]*xe/(xk*xr)*cell_f_vol[c_id];
    /* The diagonal receives eps/Rk, (*theta possibly) */
    fimp[c_id] += rhovst*thetap;
    /* The right hand side receives the dissipation */
    rhs[c_id] -= rhovst*cvara_var[c_id];
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
                        cs_real_t                   w1[],
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
  BFT_MALLOC(vistet, n_cells_ext, cs_real_6_t);

  cs_real_t *_weighb = NULL;
  cs_real_2_t *_weighf = NULL;
  cs_real_6_t *_viscce = NULL;

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
  cs_real_t *cpro_wgrec_s = NULL;
  cs_real_6_t *cpro_wgrec_v = NULL;

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
    int idifftp = eqp->idifft;
    if (scalar_turb_flux_model_type == 3)
      idifftp = 0;

    if (cpro_viscls == NULL)
      cs_array_set_value_real(n_cells, 1, visls_0, w1);
    else
      cs_array_real_copy(n_cells, cpro_viscls, w1);

    if (idifftp > 0) {
#     pragma omp parallel for if(n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        w1[c_id] += xcpp[c_id]*sgdh_diff[c_id];
      }
    }

    if (cpro_wgrec_s != NULL) {
      cs_array_real_copy(n_cells, w1, cpro_wgrec_s);
      cs_halo_sync_var(m->halo, CS_HALO_STANDARD, cpro_wgrec_s);
    }

    cs_face_viscosity(m,
                      fvq,
                      eqp->imvisf,
                      w1,
                      viscf,
                      viscb);
  }

  /* Symmetric tensor diffusivity (GGDH) */
  else if (eqp->idften & CS_ANISOTROPIC_DIFFUSION) {
    BFT_MALLOC(_weighb, n_b_faces, cs_real_t);
    BFT_MALLOC(_weighf, n_i_faces, cs_real_2_t);
    BFT_MALLOC(_viscce, n_cells_ext, cs_real_6_t);
    const cs_real_6_t *visten = NULL;
    const int kctheta = cs_field_key_id("turbulent_flux_ctheta");
    const cs_real_t ctheta = cs_field_get_key_double(f, kctheta);

    if (turb_model->iturb != CS_TURB_RIJ_EPSILON_EBRSM) {
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
      if (cpro_viscls == NULL) {
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
      if (cpro_viscls == NULL) {
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
    if (cpro_wgrec_v != NULL) {
      cs_array_real_copy(n_cells*6,
                         (const cs_real_t *)_viscce,
                         (cs_real_t *)cpro_wgrec_v);
      cs_mesh_sync_var_sym_tens(cpro_wgrec_v);
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

  BFT_FREE(vistet);

  *weighb = (cs_real_t *)_weighb;
  *weighf = (cs_real_2_t *)_weighf;
  *viscce = (cs_real_6_t *)_viscce;
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

  cs_real_t *_weighb = NULL;
  cs_real_2_t *_weighf = NULL;
  cs_real_6_t *_viscce = NULL;

  const cs_real_t *visct = CS_F_(mu_t)->val;

  /* Viscosity and diffusivity*/
  const cs_real_t *cpro_viscls = NULL;
  const int kivisl = cs_field_key_id("diffusivity_id");
  const int kvisl0 = cs_field_key_id("diffusivity_ref");
  const int ifcvsl = cs_field_get_key_int(f, kivisl);
  const cs_real_t visls_0 = cs_field_get_key_double(f, kvisl0);
  if (ifcvsl > -1)
    cpro_viscls = cs_field_by_id(ifcvsl)->val;

  /* Weighting field for gradient */
  cs_real_t *cpro_wgrec_s = NULL;
  cs_real_6_t *cpro_wgrec_v = NULL;

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
    BFT_MALLOC(w1, n_cells_ext, cs_real_t);

    if (cpro_viscls == NULL) {
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
    if (cpro_wgrec_s != NULL) {
      cs_array_real_copy(n_cells, w1, cpro_wgrec_s);
      cs_halo_sync_var(m->halo, CS_HALO_STANDARD, cpro_wgrec_s);
    }

    cs_face_viscosity(m,
                      fvq,
                      eqp->imvisf,
                      w1,
                      viscf,
                      viscb);

    BFT_FREE(w1);
  }

  /* Symmetric tensor diffusivity (GGDH) */
  else if (eqp->idften & CS_ANISOTROPIC_DIFFUSION) {

    BFT_MALLOC(_weighb, n_b_faces, cs_real_t);
    BFT_MALLOC(_weighf, n_i_faces, cs_real_2_t);
    BFT_MALLOC(_viscce, n_cells_ext, cs_real_6_t);

    const cs_real_6_t *visten = NULL;
    const int kctheta = cs_field_key_id("turbulent_flux_ctheta");
    const cs_real_t ctheta = cs_field_get_key_double(f, kctheta);

    if (turb_model->iturb != CS_TURB_RIJ_EPSILON_EBRSM) {
      cs_field_t * f_vis = cs_field_by_name("anisotropic_turbulent_viscosity");
      visten = (cs_real_6_t *)f_vis->val;
    }

    /* EBRSM and (GGDH or AFM) */
    else {
      cs_field_t * f_vis
        = cs_field_by_name("anisotropic_turbulent_viscosity_scalar");
      visten = (cs_real_6_t *)f_vis->val;
    }

    if (cpro_viscls == NULL) {
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
    if (cpro_wgrec_v != NULL) {
      cs_array_real_copy(6*n_cells,
                         (cs_real_t *)_viscce,
                         (cs_real_t *)cpro_wgrec_v);
      cs_mesh_sync_var_sym_tens(cpro_wgrec_v);
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

  *weighb = (cs_real_t *)_weighb;
  *weighf = (cs_real_2_t *)_weighf;
  *viscce = (cs_real_6_t *)_viscce;
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

void
cs_f_solve_equation_scalar(int              f_id,
                           cs_lnum_t        ncesmp,
                           int              iterns,
                           int              itspdv,
                           const cs_lnum_t  icetsm[],
                           int              itypsm[],
                           cs_real_t        smacel[],
                           cs_real_t        viscf[],
                           cs_real_t        viscb[])
{
  cs_field_t *f = cs_field_by_id(f_id);

  cs_solve_equation_scalar(f,
                           ncesmp,
                           iterns,
                           itspdv,
                           icetsm,
                           itypsm,
                           smacel,
                           viscf,
                           viscb);
}

void
cs_f_solve_equation_vector(int              f_id,
                           cs_lnum_t        ncesmp,
                           int              iterns,
                           const cs_lnum_t  icetsm[],
                           int              itypsm[],
                           cs_real_t        smacel[],
                           cs_real_t        viscf[],
                           cs_real_t        viscb[])
{
  cs_field_t *f = cs_field_by_id(f_id);

  cs_solve_equation_vector(f,
                           ncesmp,
                           iterns,
                           icetsm,
                           itypsm,
                           smacel,
                           viscf,
                           viscb);
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
 * \param[in]     ncesmp     number of cells with mass source term
 * \param[in]     iterns     Navier-Stokes iteration number
 * \param[in]     itspdv     indicator to compute production/dissipation
 *                           terms for a variance:
 *                           - 0: no
 *                           - 1: yes
 * \param[in]     icetsm     index of cells with mass source term
 * \param[in]     itypsm     type of mass source term for the variables
 * \param[in]     smacel     variable value associated to the mass source
 *                           term (for ivar=ipr, smacel is the mass flux
 *                           \f$ \Gamma^n \f$)
 * \param         viscf      visc*surface/dist at internal faces (work array)
 * \param         viscb      visc*surface/dist at boundary faces (work array)
 */
/*----------------------------------------------------------------------------*/

void
cs_solve_equation_scalar(cs_field_t        *f,
                         cs_lnum_t          ncesmp,
                         int                iterns,
                         int                itspdv,
                         const cs_lnum_t    icetsm[],
                         int                itypsm[],
                         cs_real_t          smacel[],
                         cs_real_t          viscf[],
                         cs_real_t          viscb[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_time_step_t *ts = cs_glob_time_step;
  const cs_fluid_properties_t *fluid_props = cs_glob_fluid_properties;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_real_t *volume = fvq->cell_vol;
  const cs_real_t *cell_f_vol = fvq->cell_f_vol;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;

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

  cs_real_t *dtr = NULL;
  {
    const int keycdt = cs_field_key_id("time_step_factor");
    const cs_real_t cdtvar = cs_field_get_key_double(f, keycdt);

    if (fabs(cdtvar - 1.0) > cs_math_epzero) {
      BFT_MALLOC(dtr, n_cells_ext, cs_real_t);
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
    croma = cs_field_by_id(icrom_scal)->val_pre;
  }
  else {
    icrom_scal = CS_F_(rho)->id;
  }

  /* Key id for buoyant field (inside the Navier Stokes loop) */

  const int key_is_buoyant = cs_field_key_id_try("is_buoyant");
  const int is_buoyant_fld = cs_field_get_key_int(f, key_is_buoyant);

  if (   (is_buoyant_fld == 1 && iterns == -1)
      || (is_buoyant_fld == 0 && iterns != -1)) {
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
  BFT_MALLOC(fimp, n_cells_ext, cs_real_t);

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

  {
    int nvar = 0, nscal = 0;
    cs_lnum_t ncepdp = 0;
    cs_lnum_t *icepdc = NULL;
    cs_real_t *ckupdc = NULL;

    cs_f_get_dimens(&nvar, & nscal);

    CS_PROCF(ustssc, USTSSC)(&nvar, &nscal, &ncepdp, &ncesmp, &iscal,
                             icepdc, icetsm, itypsm, dt,
                             ckupdc, smacel, rhs, fimp);
  }

  cs_user_source_terms(cs_glob_domain,
                       f->id,
                       rhs,
                       fimp);

  /* Coupling between multiple code_saturne instances */
  const int nbrcpl = cs_sat_coupling_n_couplings();
  if (nbrcpl > 0)
    cs_sat_coupling_exchange_at_cells(f->id, rhs, fimp);

  /* Store the source terms for convective limiter
     or time extrapolation for buoyant scalar */

  cs_real_t *cpro_scal_st = NULL, *cproa_scal_st = NULL;

  const int kst = cs_field_key_id_try("source_term_id");
  const int kstprv = cs_field_key_id_try("source_term_prev_id");
  const int st_id = cs_field_get_key_int(f, kst);
  const int st_prv_id = cs_field_get_key_int(f, kstprv);

  if (st_id > -1) {
    cpro_scal_st = cs_field_by_id(st_id)->val;
    /* Handle parallelism and periodicity */
    if (m->halo != NULL)
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
  cs_real_t *cpro_st = NULL;
  if (st_prv_id > -1) {
    if (iterns == -1)
      cpro_st = cproa_scal_st;
    else
      cpro_st = cpro_scal_st;
  }

  const cs_real_t thetv = eqp->thetav;
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

  const cs_real_t *temp = NULL, *tempa = NULL, *cpro_yw = NULL;
  cs_real_t *cpro_yv = NULL, *xcvv = NULL;
  if (   th_cf_model->ieos != CS_EOS_NONE
      && is_thermal_model_field
      && (   th_model->thermal_variable == CS_THERMAL_MODEL_TEMPERATURE
          || th_model->thermal_variable == CS_THERMAL_MODEL_INTERNAL_ENERGY)) {

    BFT_MALLOC(xcvv, n_cells_ext, cs_real_t);

    /* Compute cv */

    cs_field_t *f_cv = cs_field_by_name_try("isobaric_heat_capacity");
    if (f_cv != NULL) {
      cs_thermal_model_cv(f_cv->val);
      cs_array_real_copy(n_cells, f_cv->val, xcvv);
    }
    else
      cs_array_set_value_real(n_cells, 1, 1., xcvv);

    /* Note that the temperature is a postprocessing field when the thermal model
       is based on the internal energy. */

    if (th_model->thermal_variable == CS_THERMAL_MODEL_INTERNAL_ENERGY) {

      const cs_field_t *f_t = cs_field_by_name_try("temperature");

      if (f_t != NULL) {
        temp  = f_t->val;
        tempa = f_t->val_pre;
      }

      if (th_cf_model->ieos == CS_EOS_MOIST_AIR) {
        const cs_field_t *f_yv = cs_field_by_name_try("yv");
        const cs_field_t *f_yw = cs_field_by_name_try("yw");

        if (f_yv != NULL) {
          cpro_yv = f_yv->val;
        }

        if (f_yw != NULL) {
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

    BFT_MALLOC(vistot, n_cells_ext, cs_real_t);
    BFT_MALLOC(gradv, n_cells_ext, cs_real_33_t);

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

    BFT_FREE(vistot);
    BFT_FREE(gradv);
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

    if (   th_model->thermal_variable == CS_THERMAL_MODEL_TEMPERATURE
        || th_model->thermal_variable == CS_THERMAL_MODEL_INTERNAL_ENERGY) {
      cs_field_t *f_cflt = cs_field_by_name_try("cfl_t");

      if (f_cflt != NULL) {
        cs_real_t *cflt = f_cflt->val;
        cs_array_set_value_real(n_cells, 1, 0.0, cflt);

        /* Only implemented for the ideal gas equation of state. */

        if (th_cf_model->ieos == CS_EOS_IDEAL_GAS) {
          const cs_real_3_t *vel = (const cs_real_3_t *)CS_F_(vel)->val;

          const int kimasf = cs_field_key_id_try("inner_mass_flux_id");
          const int iflmas = cs_field_get_key_int(f, kimasf);
          cs_real_t *imasfl = cs_field_by_id(iflmas)->val;

          /* precaution for croma */
          cs_halo_sync_var(m->halo, CS_HALO_STANDARD, CS_F_(rho)->val_pre);

          cs_thermal_model_cflt(croma,
                                temp,
                                tempa,
                                xcvv,
                                vel,
                                imasfl,
                                cflt);
        }
      }

    }

  } /* if (thermal_model_field) */

  /* Volume coupling with Syrthes; order 2 not handled. */
  if (is_thermal_model_field)
    cs_syr_coupling_volume_source_terms(f->id, rhs, fimp);

  /* Specific physical models; order 2 not handled. */
  if (cs_glob_physical_model_flag[CS_PHYSICAL_MODEL_FLAG] > 0) {
    cs_physical_model_source_terms(iscal, rhs, fimp);

    /*! Electric arcs, Joule effect ionic conduction */
    if (   cs_glob_physical_model_flag[CS_JOULE_EFFECT] > 0
        || cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] > 0)
      cs_elec_source_terms(m, fvq, f->id, rhs);

    /*! Cooling towers */
    if (cs_glob_physical_model_flag[CS_COOLING_TOWERS] > 0)
      cs_ctwr_source_term(f->id, rhs, fimp);
  }

  /*  Rayonnement
   *  Ordre 2 non pris en compte */
  cs_real_t *cpro_tsscal = NULL;
  if (idilat > 3) {
    char fname[128];
    snprintf(fname, 128, "%s_dila_st", f->name); fname[127] = '\0';
    cpro_tsscal = cs_field_by_name_try(fname)->val;
  }

  if (cs_glob_rad_transfer_params->type > CS_RAD_TRANSFER_NONE) {

    if (is_thermal_model_field) {
      cs_rad_transfer_source_terms(rhs, fimp);

      /* Conversion Temperature -> potential Temperature (theta)
       * TODO FIXME check rhs and fimp are correctly initialized... */
      if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 0) {
        const cs_real_t *pottemp = CS_F_(t)->val;
        const cs_real_t *tempc = cs_field_by_name("real_temperature")->val;
        const cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_real_t cor_factor = pottemp[c_id] / (tempc[c_id] + tkelvi);
          rhs[c_id] *= cor_factor;
          fimp[c_id] *= cor_factor;
        }
      }

      /* Store the explicit radiative source term */
      if (idilat > 3) {
        const cs_real_t *cpro_tsre1 = cs_field_by_name("rad_st")->val;
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          cpro_tsscal[c_id] += cpro_tsre1[c_id]*cell_f_vol[c_id];
      }
    }

    /* Pulverized coal; order 2 not handled */
    if (cs_glob_physical_model_flag[CS_COMBUSTION_COAL] > -1) {
      const int nclacp = cs_glob_combustion_model->coal->nclacp;
      const int isca_ih21 = cs_field_get_key_int(CS_FI_(h2, 0), keyvar);
      const int isca_ih2nl = cs_field_get_key_int(CS_FI_(h2, nclacp-1), keyvar);
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

  cs_real_t *xcpp = NULL;
  BFT_MALLOC(xcpp, n_cells_ext, cs_real_t);

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
    if (CS_F_(cp) != NULL)
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
        && (   th_model->thermal_variable == CS_THERMAL_MODEL_TEMPERATURE
            || th_model->thermal_variable == CS_THERMAL_MODEL_ENTHALPY)) {
      const int itste = cs_glob_lagr_source_terms->itste - 1;
      const int itsti = cs_glob_lagr_source_terms->itsti - 1;
      cs_real_t *ste = cs_glob_lagr_source_terms->st_val + itste*n_cells_ext;
      cs_real_t *sti = cs_glob_lagr_source_terms->st_val + itsti*n_cells_ext;
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        rhs[c_id] += ste[c_id];
        fimp[c_id] += xcpp[c_id]*cs_math_fmax(sti[c_id], 0.0);
      }
    }
  }

  /* Mass source term.

     TODO this section could be replaced by lower-level code from
     cs_volume_mass_injection.c for the specific field being handled
     (i.e. using the direct xdef-based evaluation).
     This would allow removing the itypsm and smacel arrays and
     associated dimensions. */

  cs_real_t *w1;
  BFT_MALLOC(w1, n_cells_ext, cs_real_t);

  const int ipr = cs_field_get_key_int(CS_F_(p), keyvar);
  if (ncesmp > 0) {
    cs_real_t *srcmas;
    int *_itypsm = itypsm + (ivar-1)*ncesmp;
    cs_real_t *_smacel_ipr = smacel + (ipr-1)*ncesmp;
    cs_real_t *_smacel_ivar = smacel + (ivar-1)*ncesmp;

    BFT_MALLOC(srcmas, ncesmp, cs_real_t);

    /* When treating the Temperature, the equation is multiplied by Cp */
    for (cs_lnum_t c_id = 0; c_id < ncesmp; c_id++) {
      if ((_smacel_ipr[c_id] > 0.0) && (_itypsm[c_id] == 1)) {
        const cs_lnum_t id = icetsm[c_id] - 1;
        srcmas[c_id] = _smacel_ipr[c_id]* xcpp[id];
      }
      else {
       srcmas[c_id] = 0.;
      }
    }

    /* Increment rhs by -Gamma.cvara_var and fimp by Gamma */
    cs_mass_source_terms(1,
                         1,
                         ncesmp,
                         icetsm,
                         _itypsm,
                         cell_f_vol,
                         cvara_var,
                         _smacel_ivar,
                         srcmas,
                         rhs,
                         fimp,
                         w1);
    BFT_FREE(srcmas);

    /* If we extrapolate STs we put Gamma Pinj in cpro_st */
    if (st_prv_id >= 0)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        cpro_st[c_id] += w1[c_id];
    /* Otherwise, we put it directly in rhs */
    else
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        rhs[c_id] += w1[c_id];
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
  const cs_real_t *cpro_viscls = NULL;
  const cs_real_t *viscl = CS_F_(mu)->val;
  const int kivisl = cs_field_key_id("diffusivity_id");
  const int kvisl0 = cs_field_key_id("diffusivity_ref");
  const int ifcvsl = cs_field_get_key_int(f, kivisl);
  const cs_real_t visls_0 = cs_field_get_key_double(f, kvisl0);
  if (ifcvsl > -1)
    cpro_viscls = cs_field_by_id(ifcvsl)->val;

  /* Initialize turbulent diffusivity for SGDH model */
  cs_real_t *sgdh_diff;
  BFT_MALLOC(sgdh_diff, n_cells_ext, cs_real_t);

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
           && fluid_props->irovar == 1) {
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
  cs_real_t *weighb = NULL;
  cs_real_2_t *weighf = NULL;
  cs_real_6_t *viscce = NULL;

  if (eqp->idiff >= 1) {
    _diffusion_terms_scalar(f,
                            eqp,
                            visls_0,
                            xcpp,
                            sgdh_diff,
                            cpro_viscls,
                            w1,
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

  BFT_FREE(w1);
  BFT_FREE(sgdh_diff);

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

  if (iscdri > 0) {
    cs_real_t *divflu;
    BFT_MALLOC(divflu, n_cells_ext, cs_real_t);

    cs_drift_convective_flux(f,
                             imasfl,
                             bmasfl,
                             divflu);

    const int iconvp = eqp->iconv;
    const cs_real_t thetap = eqp->thetav;

    /*  NB: if the porosity module is switched on, the porosity is already
     * taken into account in divflu */

    /* mass aggregation term */
#   pragma omp parallel for if(n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      fimp[c_id] += iconvp*thetap*divflu[c_id];
      rhs[c_id] -= iconvp*divflu[c_id]*cvara_var[c_id];
    }

    BFT_FREE(divflu);
  }

  /* Solve
   * ===== */

  cs_real_t *cvark_var = NULL, *wcvark_var = NULL;
  if (iterns >= 1) {
    BFT_MALLOC(wcvark_var, n_cells_ext, cs_real_t);
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
                                     NULL,
                                     iescap,
                                     imucpp,
                                     normp,
                                     eqp,
                                     cvara_var,
                                     cvark_var,
                                     f->bc_coeffs->a,
                                     f->bc_coeffs->b,
                                     f->bc_coeffs->af,
                                     f->bc_coeffs->bf,
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
                                     NULL,
                                     fimp,
                                     rhs,
                                     cvar_var,
                                     dpvar,
                                     xcpp,
                                     NULL);

  CS_FREE_HD(dpvar);
  if (weighb != NULL) {
    BFT_FREE(weighb);
    BFT_FREE(weighf);
    BFT_FREE(viscce);
  }

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
                                NULL,
                                cvar_var,
                                CS_F_(p)->val,
                                CS_F_(p)->val_pre,
                                cpro_yw,
                                cpro_yv,
                                tempk);

    }

  }

  BFT_FREE(xcvv);

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
      BFT_MALLOC(pphy, n_cells_ext, cs_real_t);

      if (cs_glob_atmo_option->meteo_profile == 0) {
        cs_real_t _pphy, dum;
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_atmo_profile_std(cell_cen[c_id][2], &_pphy, &dum, &dum);
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
      BFT_FREE(pphy);

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
      else if (turb_model->iturb == CS_TURB_K_OMEGA) {
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
    BFT_MALLOC(errork, n_cells_ext, cs_real_t);
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      errork[c_id] = cvar_var[c_id] - cvark_var[c_id];

    double l2errork = sqrt(cs_gres(n_cells, cell_f_vol, errork, errork));
    BFT_FREE(errork);

    double l2norm = sqrt(cs_gres(n_cells, cell_f_vol, cvara_var, cvara_var));
    double dl2norm = 1.0;
    if (l2norm > 0.)
      dl2norm = 1.0/l2norm;

    bft_printf("Inner iteration scalar %10d iter=%10d L2 error = %12.4e\n"
               "L2 normalized error %e12.4 nomr %e12.4\n",f->id, iterns,
               l2errork, l2errork*dl2norm, l2norm);
  }

  BFT_FREE(wcvark_var);
  BFT_FREE(xcpp);
  CS_FREE_HD(rhs);
  BFT_FREE(fimp);
  BFT_FREE(dtr);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve the convection/diffusion equation (with optional source
 *        terms and/or drift) for a vectorial quantity over a time step..
 *
 * \param[in]     f          pointer to field structure
 * \param[in]     ncesmp     number of cells with mass source term
 * \param[in]     iterns     Navier-Stokes iteration number
 * \param[in]     icetsm     index of cells with mass source term
 * \param[in]     itypsm     type of mass source term for the variables
 * \param[in]     smacel     variable value associated to the mass source
 *                           term (for ivar=ipr, smacel is the mass flux
 *                           \f$ \Gamma^n \f$)
 * \param         viscf      visc*surface/dist at internal faces (work array)
 * \param         viscb      visc*surface/dist at boundary faces (work array)
 */
/*----------------------------------------------------------------------------*/

void
cs_solve_equation_vector(cs_field_t       *f,
                         const cs_lnum_t   ncesmp,
                         int               iterns,
                         const cs_lnum_t   icetsm[],
                         int               itypsm[],
                         cs_real_t         smacel[],
                         cs_real_t         viscf[],
                         cs_real_t         viscb[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_fluid_properties_t *fluid_props = cs_glob_fluid_properties;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_real_t *cell_f_vol = fvq->cell_f_vol;

  const cs_real_t *dt = CS_F_(dt)->val;

  /* Initialization
     -------------- */

  cs_real_3_t *cvar_var = (cs_real_3_t *)f->val;
  cs_real_3_t *cvara_var = (cs_real_3_t *)f->val_pre;

  const int keyvar = cs_field_key_id("variable_id");
  const int ivar = cs_field_get_key_int(f, keyvar);

  /* If the vector is buoyant, it is inside the Navier Stokes loop, and
     so iterns >=1; otherwise it is outside of the loop and iterns = -1. */

  const int key_is_buoyant = cs_field_key_id_try("is_buoyant");
  const int is_buoyant_fld = cs_field_get_key_int(f, key_is_buoyant);

  if (   (is_buoyant_fld == 1 && iterns == -1)
      || (is_buoyant_fld == 0 && iterns != -1)) {
    return;
  }

  /* We might have a time step multiplier */

  cs_real_t *dtr = NULL;
  {
    const int keycdt = cs_field_key_id("time_step_factor");
    const cs_real_t cdtvar = cs_field_get_key_double(f, keycdt);

    if (fabs(cdtvar - 1.0) > cs_math_epzero) {
      BFT_MALLOC(dtr, n_cells_ext, cs_real_t);
#     pragma omp parallel for if(n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
        dtr[c_id] = dt[c_id] * cdtvar;
      dt = dtr;
    }
  }

  /* Physical quantities */

  const cs_real_t *croma = CS_F_(rho)->val;
  if (CS_F_(rho)->val_pre != NULL)
    croma = CS_F_(rho)->val_pre;

  cs_equation_param_t *eqp = cs_field_get_equation_param(f);

  if (eqp->verbosity > 0)
    bft_printf(" ** SOLVING VARIABLE %s\n"
               "    ----------------\n\n", f->name);

  /* Source terms
     ============ */

  cs_real_3_t *rhs;
  cs_real_33_t *fimp;

  BFT_MALLOC(fimp, n_cells_ext, cs_real_33_t);
  BFT_MALLOC(rhs, n_cells_ext, cs_real_3_t);

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
    cs_sat_coupling_exchange_at_cells(f->id,
                                      (cs_real_t *)rhs,
                                      (cs_real_t *)fimp);

  /* Store the source terms for convective limiter */
  cs_real_3_t *cpro_vect_st = NULL;

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
    cs_mesh_sync_var_vect((cs_real_t *)cpro_vect_st);
  }

  /* If we extrapolate source terms:
   *   rhs receives -theta TS from previous time step.
   *   rhs receives the part of the source term which depends on the variable.
   *   At order 2, we assume that fimp provided by the user is < 0, so
   *     we implicit the term (so fimp*cvar_prev goes in rhs).
   *   In standard case, adapt to sign of fimp, but fimp*cvar_prev still
   *     goes in rhs (no other choice). */

  cs_real_3_t *cproa_vect_st = NULL;
  const cs_real_t thetv = eqp->thetav;
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

    if (cs_glob_physical_model_flag[CS_COOLING_TOWERS] > 0)
      cs_ctwr_source_term(f->id, (cs_real_t *)rhs, (cs_real_t *)fimp);

    if (   cs_glob_physical_model_flag[CS_JOULE_EFFECT] > 0
        || cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] > 0)
      cs_elec_source_terms_v(m, fvq, f->id, rhs);

  }

  /* Mass source term */

  if (ncesmp > 0) {
    cs_real_3_t *gavinj;
    BFT_MALLOC(gavinj, n_cells_ext, cs_real_3_t);

    const int ipr = cs_field_get_key_int(CS_F_(p), keyvar);

    int *_itypsm = itypsm + (ivar-1)*ncesmp;
    cs_real_t *_smacel_ipr = smacel + (ipr-1)*ncesmp;
    cs_real_t *_smacel_ivar = smacel + (ivar-1)*ncesmp;

    /* We increment SMBRV by -Gamma RTPA and FIMP by Gamma */
    cs_mass_source_terms(1,
                         3,
                         ncesmp,
                         icetsm,
                         _itypsm,
                         cell_f_vol,
                         (const cs_real_t *)cvara_var,
                         _smacel_ivar,
                         _smacel_ipr,
                         (cs_real_t *)rhs,
                         (cs_real_t *)fimp,
                         (cs_real_t *)gavinj);

    /* If we extrapolate source terms we put -Gamma Pinj in cproa_vect_st */
    if (st_prv_id > -1) {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        for (cs_lnum_t j = 0; j < 3; j++)
          cproa_vect_st[c_id][j] += gavinj[c_id][j];
      }
    }
    /* If we do not extrapolate we put directly in SMBRV */
    else {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        for (cs_lnum_t j = 0; j < 3; j++)
          rhs[c_id][j] += gavinj[c_id][j];
      }
    }
    BFT_FREE(gavinj);
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

  const cs_real_t *pcrom = NULL;

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

  cs_real_t *weighb = NULL;
  cs_real_2_t *weighf = NULL;
  cs_real_6_t *viscce = NULL;

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
      = (const cs_real_3_t *restrict)fvq->b_face_normal;
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

  if (iscdri > 0) {
    cs_real_t *divflu;
    BFT_MALLOC(divflu, n_cells_ext, cs_real_t);

    cs_drift_convective_flux(f,
                             imasfl,
                             bmasfl,
                             divflu);

    const int iconvp = eqp->iconv;
    const cs_real_t thetap = eqp->thetav;

    /* NB: if the porosity module is swiched on, the the porosity is already
     * taken into account in divflu */

    /* mass aggregation term */
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        fimp[c_id][j][j] += iconvp * thetap * divflu[c_id];
        rhs[c_id][j] -= iconvp * divflu[c_id] * cvara_var[c_id][j];
      }
    }

    BFT_FREE(divflu);
  }

  /* Solve
     ===== */

  int iescap = 0, icvflb = 0, ivissv = 0;

  cs_equation_iterative_solve_vector(cs_glob_time_step_options->idtvar,
                                     iterns,
                                     f->id,
                                     NULL,
                                     ivissv,
                                     iescap,
                                     eqp,
                                     cvara_var,
                                     cvara_var,
                                     (const cs_real_3_t *)f->bc_coeffs->a,
                                     (const cs_real_33_t *)f->bc_coeffs->b,
                                     (const cs_real_3_t *)f->bc_coeffs->af,
                                     (const cs_real_33_t *)f->bc_coeffs->bf,
                                     imasfl,
                                     bmasfl,
                                     viscf,
                                     viscb,
                                     viscf,
                                     viscb,
                                     NULL,
                                     NULL,
                                     viscce,
                                     weighf,
                                     weighb,
                                     icvflb,
                                     NULL,
                                     fimp,
                                     rhs,
                                     cvar_var,
                                     NULL);

  if (weighb != NULL) {
    BFT_FREE(weighb);
    BFT_FREE(weighf);
    BFT_FREE(viscce);
  }

  BFT_FREE(fimp);

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

  BFT_FREE(rhs);
  BFT_FREE(dtr);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
