/*============================================================================
 * k-w turbulence model.
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

#include "cs_array.h"
#include "cs_balance.h"
#include "cs_blas.h"
#include "cs_halo.h"
#include "cs_base.h"
#include "cs_dispatch.h"
#include "cs_equation.h"
#include "cs_equation_iterative_solve.h"
#include "cs_face_viscosity.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_gradient.h"
#include "cs_lagr.h"
#include "cs_log.h"
#include "cs_log_iteration.h"
#include "cs_mass_source_terms.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_physical_constants.h"
#include "cs_porous_model.h"
#include "cs_prototypes.h"
#include "cs_rotation.h"
#include "cs_thermal_model.h"
#include "cs_time_step.h"
#include "cs_turbulence_model.h"
#include "cs_turbulence_rotation.h"
#include "cs_volume_mass_injection.h"
#include "cs_velocity_pressure.h"
#include "cs_wall_functions.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_turbulence_kw.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_turbulence_kw.c

  Solve the \f$ k - \omega \f$ SST for incompressible flows
  or slightly compressible flows for one time step.
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
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
 * \brief Solve the k-omega equations.
 *
 * Solve the \f$ k - \omega \f$ SST for incompressible flows
 * or slightly compressible flows for one time step.
 *
 * \param[in]     phase_id      turbulent phase id (-1 for single phase flow)
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_kw(int phase_id)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  cs_domain_t  *domain = cs_glob_domain;

  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_real_t *volume  = fvq->cell_vol;
  const cs_real_t *cell_f_vol = fvq->cell_f_vol;
  const cs_real_t *distb = fvq->b_dist;
  const cs_lnum_t *b_face_cells = m->b_face_cells;

  const cs_time_scheme_t *time_scheme = cs_glob_time_scheme;
  const cs_real_t thets  = time_scheme->thetst;

  const cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();
  cs_real_t viscl0 = phys_pro->viscl0; /* reference molecular viscosity */
  cs_real_t ro0 = phys_pro->ro0; /* reference density */
  const cs_real_t uref = cs_glob_turb_ref_values->uref;

  cs_field_t *f_k = CS_F_(k);
  cs_field_t *f_omg = CS_F_(omg);
  cs_field_t *f_vel = CS_F_(vel);
  cs_field_t *f_eps = CS_F_(eps);
  cs_field_t *f_mu = CS_F_(mu);
  cs_field_t *f_mut = CS_F_(mu_t);
  cs_field_t *f_rho = CS_F_(rho);
  cs_field_t *f_rhob = CS_F_(rho_b);

  if (phase_id >= 0) {
    f_k = CS_FI_(k, phase_id);
    f_omg = CS_FI_(omg, phase_id);
    f_vel = CS_FI_(vel, phase_id);
    f_eps = CS_FI_(eps, phase_id);
    f_mu = CS_FI_(mu, phase_id);
    f_mut = CS_FI_(mu_t, phase_id);
    f_rho = CS_FI_(rho, phase_id);
    f_rhob = CS_FI_(rho_b, phase_id);
  }

  const cs_real_t *dt = CS_F_(dt)->val;
  const cs_real_t grav[3] = {cs_glob_physical_constants->gravity[0],
                             cs_glob_physical_constants->gravity[1],
                             cs_glob_physical_constants->gravity[2]};

  cs_dispatch_context ctx;

  /* Initialization
     ============== */

  /* Turbulence model constants */

  const cs_real_t cmu = cs_turb_cmu;
  const cs_real_t ckwgm1 = cs_turb_ckwgm1;
  const cs_real_t ckwgm2 = cs_turb_ckwgm2;
  const cs_real_t ckwbt1 = cs_turb_ckwbt1;
  const cs_real_t ckwbt2 = cs_turb_ckwbt2;
  const cs_real_t ckwsw2 = cs_turb_ckwsw2;
  const cs_real_t xkappa = cs_turb_xkappa;

  const int hybrid_turb = cs_glob_turb_model->hybrid_turb;

  /* Allocate temporary arrays for the turbulence resolution */

  cs_real_t *viscf, *viscb;
  CS_MALLOC_HD(viscf, n_i_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(viscb, n_b_faces, cs_real_t, cs_alloc_mode);

  cs_real_t *smbrk, *smbrw;
  CS_MALLOC_HD(smbrk, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(smbrw, n_cells_ext, cs_real_t, cs_alloc_mode);

  cs_real_t *tinstk, *tinstw, *xf1;
  CS_MALLOC_HD(tinstk, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(tinstw, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(xf1, n_cells_ext, cs_real_t, cs_alloc_mode);

  /* Allocate work arrays */

  cs_real_t *w1, *dpvar, *gdkgdw, *prodk, *prodw;
  cs_field_t *f_tke_prod = cs_field_by_name_try("algo:k_production");
  cs_field_t *f_tke_buoy = cs_field_by_name_try("algo:k_buoyancy");
  CS_MALLOC_HD(dpvar, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(w1, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(gdkgdw, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(prodk, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(prodw, n_cells_ext, cs_real_t, cs_alloc_mode);

  cs_real_t *s2pw2 = nullptr;
  cs_real_t *maxgdsv = nullptr;
  cs_real_t *d2uidxi2 = nullptr;

  if (hybrid_turb == CS_HYBRID_DDES) {
    /* DDES hybrid model */
    CS_MALLOC_HD(s2pw2, n_cells_ext, cs_real_t, cs_alloc_mode);
  }
  else if (hybrid_turb == CS_HYBRID_SAS) {
    /* SAS hybrid model */
    CS_MALLOC_HD(maxgdsv, n_cells_ext, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(d2uidxi2, n_cells_ext, cs_real_t, cs_alloc_mode);
  }

  const cs_real_t *cvisct = (const cs_real_t *)f_mut->val;
  const cs_real_t *viscl = (const cs_real_t *)f_mu->val;

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  const cs_real_t *i_massflux
    = cs_field_by_id(cs_field_get_key_int(f_vel, kimasf))->val;
  const cs_real_t *b_massflux
    = cs_field_by_id(cs_field_get_key_int(f_vel, kbmasf))->val;

  cs_real_t *crom = (cs_real_t *)f_rho->val;
  cs_real_t *cromo = (cs_real_t *)f_rho->val;
  cs_real_t *bromo = (cs_real_t *)f_rhob->val;

  cs_real_t *cpro_pcvto = (cs_real_t *)f_mut->val;
  cs_real_t *cpro_pcvlo = (cs_real_t *)f_mu->val;

  cs_real_t *cvar_k = (cs_real_t *)f_k->val;
  cs_real_t *cvara_k = (cs_real_t *)f_k->val_pre;
  cs_real_t *cvar_omg = (cs_real_t *)f_omg->val;
  cs_real_t *cvara_omg = (cs_real_t *)f_omg->val_pre;

  cs_real_t *coefa_k = (cs_real_t *)f_k->bc_coeffs->a;
  cs_real_t *coefb_k = (cs_real_t *)f_k->bc_coeffs->b;
  cs_real_t *coefaf_k = (cs_real_t *)f_k->bc_coeffs->af;
  cs_real_t *coefbf_k = (cs_real_t *)f_k->bc_coeffs->bf;

  cs_real_t *coefa_o = (cs_real_t *)f_omg->bc_coeffs->a;
  cs_real_t *coefb_o = (cs_real_t *)f_omg->bc_coeffs->b;
  cs_real_t *coefaf_o = (cs_real_t *)f_omg->bc_coeffs->af;
  cs_real_t *coefbf_o = (cs_real_t *)f_omg->bc_coeffs->bf;

  cs_real_t *tke_prod = nullptr;
  cs_real_t *tke_buoy = nullptr;
  if (f_tke_prod != nullptr)
    tke_prod = f_tke_prod->val;
  if (f_tke_buoy != nullptr)
    tke_buoy = f_tke_buoy->val;

  const cs_equation_param_t *eqp_k
    = cs_field_get_equation_param_const(f_k);

  const cs_equation_param_t *eqp_w
    = cs_field_get_equation_param_const(f_omg);

  cs_equation_param_t eqp_k_loc = *eqp_k;
  eqp_k_loc.idften = CS_ISOTROPIC_DIFFUSION;

  cs_equation_param_t eqp_w_loc = *eqp_w;
  eqp_w_loc.idften = CS_ISOTROPIC_DIFFUSION;

  const cs_equation_param_t *eqp_u
    = cs_field_get_equation_param_const(f_vel);

  cs_field_t *f_s2kw = cs_field_by_name_try("s2");
  cs_field_t *f_divukw = cs_field_by_name_try("vel_gradient_trace");

  if (phase_id >= 0) {
    char f_name[64]; /* should be much larger than needed */

    snprintf(f_name, 63, "s2_%d", phase_id + 1);
    f_name[63] = '\0';
    f_s2kw = cs_field_by_name(f_name);

    snprintf(f_name, 63, "vel_gradient_trace_%d", phase_id + 1);
    f_divukw = cs_field_by_name(f_name);
  }

  const cs_real_t *cpro_s2kw = f_s2kw->val;
  const cs_real_t *cpro_divukw = f_divukw->val;

  int kstprv = cs_field_key_id("source_term_prev_id");
  int istprv =  cs_field_get_key_int(f_k, kstprv);
  cs_real_t *c_st_k_p = nullptr;
  cs_real_t *c_st_omg_p = nullptr;

  if (istprv >= 0) {
    c_st_k_p = cs_field_by_id(istprv)->val;
    istprv = cs_field_get_key_int(f_eps, kstprv);
    if (istprv >= 0) {
      c_st_omg_p = cs_field_by_id(istprv)->val;
    }
    if (istprv >= 0)
      istprv = 1;
  }

  /* Time extrapolation? */

  int key_t_ext_id = cs_field_key_id("time_extrapolated");

  int iroext, iviext;
  if (istprv >= 0) {
    iroext = cs_field_get_key_int(f_rho, key_t_ext_id);
    if (iroext > 0) {
      cromo = (cs_real_t *)f_rho->val_pre;
      bromo = (cs_real_t *)f_rhob->val_pre;
    }
    iviext =  cs_field_get_key_int(f_mu, key_t_ext_id);
    if (iviext > 0)
      cpro_pcvlo = (cs_real_t *)f_mu->val_pre;
    iviext = cs_field_get_key_int(f_mut, key_t_ext_id);
    if (iviext > 0) {
      cpro_pcvto = (cs_real_t *)f_mut->val_pre;
    }
  }

  const cs_real_t *w_dist = cs_field_by_name("wall_distance")->val;

  if (eqp_k->verbosity >= 1)
    cs_log_printf(CS_LOG_DEFAULT,
                  "\n"
                  "  ** Solving k-omega\n"
                  "     ---------------\n");

  /* Take user source terms into account
   * ===================================
   *
   * explicit parts stored in: smbrk, smbrw
   * implicit parts stored in: usimpk, usimpw */

  cs_real_t *usimpk;
  cs_real_t *usimpw;
  CS_MALLOC_HD(usimpk, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(usimpw, n_cells_ext, cs_real_t, cs_alloc_mode);

  cs_arrays_set_value<cs_real_t, 1>(n_cells, 0.,
                                    smbrk, smbrw,
                                    usimpk, usimpw,
                                    tinstk, tinstw);

  cs_user_source_terms(domain,
                       f_k->id,
                       smbrk,
                       usimpk);

  cs_user_source_terms(domain,
                       f_omg->id,
                       smbrw,
                       usimpw);

  if (cs_glob_porous_model == 3) {
    cs_immersed_boundary_wall_functions(f_k->id, smbrk, usimpk);
    cs_immersed_boundary_wall_functions(f_omg->id, smbrw, usimpw);
  }

  /* If source terms are extrapolated over time */

  if (istprv >= 0) {

    const cs_real_t thetak = eqp_k->theta;
    const cs_real_t thetaw = eqp_w->theta;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      /* Recover the value at time (n-1) */
      cs_real_t tuexpk = c_st_k_p[c_id];
      cs_real_t tuexpw = c_st_omg_p[c_id];
      /* Save the values for the next time-step */
      c_st_k_p[c_id] = smbrk[c_id];
      c_st_omg_p[c_id] = smbrw[c_id];
      /* Explicit Part */
      smbrk[c_id] = - thets*tuexpk;
      smbrw[c_id] = - thets*tuexpw;
      /* It is assumed that (-usimpk > 0) and though this term is implicit */
      smbrk[c_id] = usimpk[c_id]*cvara_k[c_id] + smbrk[c_id];
      smbrw[c_id] = usimpw[c_id]*cvara_omg[c_id] + smbrw[c_id];
      /* Implicit part */
      tinstk[c_id] -= usimpk[c_id]*thetak;
      tinstw[c_id] -= usimpw[c_id]*thetaw;
    });
  }

  /* If no extrapolation over time */
  else {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      /* Explicit Part */
      smbrk[c_id] += usimpk[c_id]*cvara_k[c_id];
      smbrw[c_id] += usimpw[c_id]*cvara_omg[c_id];
      /* Implicit part */
      tinstk[c_id] += cs_math_fmax(-usimpk[c_id], 0.);
      tinstw[c_id] += cs_math_fmax(-usimpw[c_id], 0.);
    });
  }

  const cs_real_t d2s3 = 2./3.;
  const cs_real_t epz2 = cs_math_pow2(cs_math_epzero);

  /* Compute dk/dxj.dw/dxj stored in gdkgdw
     ====================================== */

  /* Allocate temporary arrays for gradients calculation */

  cs_real_3_t *gradk, *grado;
  CS_MALLOC_HD(gradk, n_cells_ext, cs_real_3_t, cs_alloc_mode);
  CS_MALLOC_HD(grado, n_cells_ext, cs_real_3_t, cs_alloc_mode);

  bool use_previous_t = true;

  cs_field_gradient_scalar(f_k,
                           use_previous_t,
                           1,     /* inc */
                           gradk);

  cs_field_gradient_scalar(f_omg,
                           use_previous_t,
                           1,     /* inc */
                           grado);

  /* Initialization of work arrays in case of Hybrid turbulence modelling */

  cs_real_t *hybrid_fd_coeff = nullptr;
  cs_real_t *sas_source_term = nullptr;
  cs_real_t *psi       = nullptr;
  cs_real_t *htles_t         = nullptr;
  cs_real_t *htles_kwsst_f1  = nullptr;
  if (hybrid_turb != CS_HYBRID_NONE) {
    /* For all hybrid model, possibility to have hybrid scheme */
    hybrid_fd_coeff = cs_field_by_name("hybrid_blend")->val;
  }

  if (hybrid_turb == CS_HYBRID_DDES) {
    /* DDES hybrid model (Menter et al.)
       Computation of dU_i/dx_j * dU_i/dx_j */

    /* Computation of the velocity gradient */

    cs_real_33_t *gradv = nullptr, *_gradv = nullptr;
    {
      cs_field_t *f_vg = cs_field_by_name_try("algo:velocity_gradient");

      if (f_vel->grad != nullptr)
        gradv = (cs_real_33_t *)f_vel->grad;
      else if (f_vg != nullptr)
        gradv = (cs_real_33_t *)f_vg->val;
      else {
        CS_MALLOC_HD(_gradv, n_cells_ext, cs_real_33_t, cs_alloc_mode);
        gradv = _gradv;
      }
    }

    if (f_vel->grad == nullptr) {

      cs_field_gradient_vector(f_vel,
                               true,  /* use_previous_t */
                               1,     /* inc */
                               gradv);
    }
    else
      gradv = (cs_real_33_t *)f_vel->grad;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      s2pw2[c_id] =   cs_math_3_square_norm(gradv[c_id][0])
                    + cs_math_3_square_norm(gradv[c_id][1])
                    + cs_math_3_square_norm(gradv[c_id][2]);
    });

    ctx.wait();

    CS_FREE_HD(_gradv);
  }
  else if (hybrid_turb == CS_HYBRID_SAS) {

    sas_source_term = cs_field_by_name("hybrid_sas_source_term")->val;

    /* Scale Adaptive hybrid model (Menter et al.)
       Computation of max(|grad k|^2 /k^2 , |grad w|^2/w^2 )
       and div(grad (U) ) */

    cs_real_3_t *cvar_vel = (cs_real_3_t *)f_vel->val;
    const cs_real_3_t *cvar_vela = (const cs_real_3_t *)f_vel->val_pre;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

      cs_real_t xgdk2sk =   cs_math_3_square_norm(gradk[c_id])
                          / cs_math_pow2(cvar_k[c_id]);

      cs_real_t xgdw2sw =  cs_math_3_square_norm(grado[c_id])
                         / cs_math_pow2(cvar_omg[c_id]);

      maxgdsv[c_id] = CS_MAX(xgdk2sk, xgdw2sw);

      /* Viscosity for the following computation of the velocity Laplacian */
      w1[c_id] = 1.;

    });

    ctx.wait();

    cs_face_viscosity(m,
                      fvq,
                      eqp_u->imvisf,
                      w1,
                      viscf,
                      viscb);

    int icvflb = 0;
    int ivisep = 0;

    cs_equation_param_t eqp_u_loc = *eqp_u;
    eqp_u_loc.iconv = 0;

    cs_real_3_t *vel_laplacian;
    CS_MALLOC_HD(vel_laplacian, n_cells_ext, cs_real_3_t, cs_alloc_mode);
    cs_arrays_set_value<cs_real_t, 1>(3*n_cells_ext, 0., (cs_real_t *)vel_laplacian);

    cs_balance_vector(cs_glob_time_step_options->idtvar,
                      -1, /* f_id */
                      0,  /* imasac */
                      1,  /* inc */
                      ivisep,
                      &eqp_u_loc,
                      cvar_vel,
                      cvar_vela,
                      f_vel->bc_coeffs,
                      i_massflux,
                      b_massflux,
                      viscf,
                      viscb,
                      nullptr,
                      nullptr,
                      nullptr,
                      nullptr,
                      nullptr,
                      icvflb,
                      nullptr,
                      nullptr,
                      nullptr,
                      vel_laplacian);

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cs_real_t v_lap[3];
      for (cs_lnum_t i = 0; i < 3; i++)
        v_lap[i] = vel_laplacian[c_id][i] / volume[c_id];
      d2uidxi2[c_id] = cs_math_3_square_norm(v_lap);
    });

    ctx.wait();

    CS_FREE_HD(vel_laplacian);

  }
  else if (hybrid_turb == CS_HYBRID_HTLES) {
    /* HTLES model [Manceau; 2018] */

    psi = cs_field_by_name("htles_psi")->val;
    htles_t   = cs_field_by_name("htles_t")->val;

    htles_kwsst_f1 = cs_field_by_name("f1_kwsst")->val;
  }

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    gdkgdw[c_id] = cs_math_3_dot_product(gradk[c_id], grado[c_id]);
  });

  ctx.wait();

  /* Free memory */
  CS_FREE_HD(gradk);
  CS_FREE_HD(grado);

  /* Compute the weight f1 (stored in xf1)
     ===================================== */

  /* For the 2nd order case, we use the values at n because the term
     in (1-f1)*gdkgdw will be a property. So we still have some "constants"
     intervening in the terms in n+1/2 (ex sigma_k for diffusion) computed
     from f at n.
     -> but the effect on "constants" is small
     -> to remember if we really use second-order in time in k-omega */

  {
    const cs_real_t ckwc1 = cs_turb_ckwc1;
    const cs_real_t k_istat = eqp_k->istat;
    const cs_real_t w_istat = eqp_w->istat;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cs_real_t ro = cromo[c_id];
      cs_real_t xnu = cpro_pcvlo[c_id]/ro;
      cs_real_t xk = cvara_k[c_id];
      cs_real_t xw  = cvara_omg[c_id];
      cs_real_t cdkw = 2*ro/ckwsw2/xw*gdkgdw[c_id];
      cs_real_t xeps = cmu*xw*xk;
      cs_real_t visct = cpro_pcvto[c_id];
      cdkw = fmax(cdkw, 1.e-20);
      cs_real_t distf = fmax(w_dist[c_id], cs_math_epzero);
      cs_real_t xarg1 = fmax(sqrt(xk)/cmu/xw/distf,
                             500.*xnu/xw/cs_math_pow2(distf));
      xarg1 = fmin(xarg1,
                   4.*ro*xk/ckwsw2/cdkw/cs_math_pow2(distf));
      xf1[c_id] = tanh(cs_math_pow4(xarg1));
      if (hybrid_turb == CS_HYBRID_HTLES) {
        htles_kwsst_f1[c_id] = xf1[c_id];
      }

      /* Unsteady terms (stored in tinstk and tinstw)
         ============================================ */

      cs_real_t romvsd = crom[c_id]*cell_f_vol[c_id]/dt[c_id];
      tinstk[c_id] += k_istat*romvsd;
      tinstw[c_id] += w_istat*romvsd;

    /* Compute production terms
     * ========================
     *
     * stored in: prodk,prodw
     * At the end of this step, we keep gdkgdw, xf1, prodk, tinstW */

      /* k / (mu_T * omega) , clipped to 1 if mu_t is zero */
      if (hybrid_turb == CS_HYBRID_HTLES) {
        xeps = cmu*xw*xk*psi[c_id];
      }

      cs_real_t k_dmut_dom;
      if (visct*xw <= cs_math_epzero*xk)
        k_dmut_dom = 1.;
      else
        k_dmut_dom = xk / (visct*xw);

      prodw[c_id] = visct*cpro_s2kw[c_id] - d2s3*ro*xk*cpro_divukw[c_id];

      /* The negative part is implicit */
      cs_real_t xxf1   = xf1[c_id];
      cs_real_t xgamma = xxf1*ckwgm1 + (1.-xxf1)*ckwgm2;
      tinstw[c_id] += fmax(  d2s3*ro*cell_f_vol[c_id]
                           * (ro*xgamma*k_dmut_dom)*cpro_divukw[c_id], 0.);

      /* Take the min between prodw and the low Reynolds one */
      if (prodw[c_id] > ckwc1*ro*xeps) {
        prodk[c_id] = ckwc1*ro*xeps;
      }
      else {
        prodk[c_id] = prodw[c_id];
        tinstk[c_id] += fmax(d2s3*cell_f_vol[c_id]*ro*cpro_divukw[c_id], 0.);
      }

      /* Save production for post processing */
      if (f_tke_prod != nullptr)
        tke_prod[c_id] = prodk[c_id] / ro;

    });
  }

  ctx.wait();

  /* Take into account rotation/curvature correction, if necessary
     ============================================================ */

  /* Spalart-Shur correction: the production terms are multiplied by a
     'rotation function' */

  cs_real_t *rotfct = nullptr;
  if (cs_glob_turb_rans_model->irccor == 1) {

    /* Allocate an array for the rotation function */
    CS_MALLOC_HD(rotfct, n_cells, cs_real_t, cs_alloc_mode);

    /* Compute the rotation function (gdkgdw array not used) */
    cs_turbulence_rotation_correction(dt, rotfct, gdkgdw);

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      prodk[c_id] *= rotfct[c_id];
      prodw[c_id] *= rotfct[c_id];
    });

    /* rotfct array is used later in case of renforced coupling (ikecou = 1).
       The array is deallocated at the end of the subroutine. */

    ctx.wait();
  }

  /* Compute buoyancy terms
   * ======================
   *
   * stored in: prodk, prodw, grad_dot_g */

  cs_real_t *grad_dot_g = nullptr;
  CS_MALLOC_HD(grad_dot_g, n_cells_ext, cs_real_t, cs_alloc_mode);

  if (cs_glob_turb_rans_model->has_buoyant_term == 1) {

    /* Allocate a temporary array for the gradient calculation */

    cs_real_3_t *grad;
    CS_MALLOC_HD(grad, n_cells_ext, cs_real_3_t, cs_alloc_mode);

    /* Buoyancy production
       ------------------- */

    /* prodk=min(P,c1*eps)+G */
    /* prodw=P+(1-ce3)*G */

    cs_real_t prdtur = 1;

    cs_field_t *f_thm = cs_thermal_model_field();
    if (phase_id >= 0)
      f_thm = CS_FI_(h_tot, phase_id);

    if (f_thm != nullptr) {
      prdtur = cs_field_get_key_double(f_thm,
                                       cs_field_key_id("turbulent_schmidt"));
    }

    /* Buoyant term:     G = Beta*g*GRAD(T)/PrT/ro
       Here is computed: G =-g*GRAD(ro)/PrT/ro */

    /* Boussinesq approximation, only for the thermal scalar for the moment */

    if (cs_glob_velocity_pressure_model->idilat == 0) {

      cs_field_gradient_scalar(f_thm,
                               true,
                               1,    /* inc */
                               grad);

      /* FIXME make it dependant on the scalar and use coupled_with_vel_p field */
      cs_real_t *cpro_beta = cs_field_by_name("thermal_expansion")->val;

      /* - Beta grad(T) . g / Pr_T */
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        grad_dot_g[c_id] = - cpro_beta[c_id]
                             * (cs_math_3_dot_product(grad[c_id], grav) / prdtur);
      });

    }
    else {

      /* BCs on rho: Dirichlet ROMB
         NB: viscb is used as COEFB */

      cs_arrays_set_value<cs_real_t, 1>(n_b_faces, 0., viscb);

      cs_halo_type_t halo_type = CS_HALO_STANDARD;
      cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

      cs_gradient_type_by_imrgra(eqp_k->imrgra,
                                 &gradient_type,
                                 &halo_type);

      cs_field_bc_coeffs_t bc_coeffs_loc;
      cs_field_bc_coeffs_init(&bc_coeffs_loc);
      bc_coeffs_loc.a = bromo;
      bc_coeffs_loc.b = viscb;

      cs_gradient_scalar("cromo_grad",
                         gradient_type,
                         halo_type,
                         1,     /* inc */
                         eqp_k->nswrgr,
                         0,
                         1,             /* w_stride */
                         eqp_k->verbosity,
                         static_cast<cs_gradient_limit_t>(eqp_k->imligr),
                         eqp_k->epsrgr,
                         eqp_k->climgr,
                         nullptr,
                         &bc_coeffs_loc,
                         cromo,
                         nullptr,
                         nullptr, /* internal coupling */
                         grad);

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        cs_real_t rho = cromo[c_id];
        grad_dot_g[c_id] = cs_math_3_dot_product(grad[c_id], grav)/(rho*prdtur);
      });
    }

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cs_real_t visct = cpro_pcvto[c_id];
      cs_real_t rho = cromo[c_id];
      prodw[c_id] += visct*cs_math_fmax(-grad_dot_g[c_id], 0.);
      prodk[c_id] -= visct*grad_dot_g[c_id];

      /* Implicit Buoyant terms when negative */
      tinstk[c_id]
        += CS_MAX(cell_f_vol[c_id]*visct/cvara_k[c_id]*grad_dot_g[c_id], 0.);

      /* Save for post processing */
      if (f_tke_buoy != nullptr)
        tke_buoy[c_id] = -visct*grad_dot_g[c_id]/rho;

    });

    ctx.wait();

    /* Free memory */
    CS_FREE_HD(grad);
  }
  else {
    cs_arrays_set_value<cs_real_t, 1>(n_cells, 0., grad_dot_g);
  }

  /* Finalization of explicit and implicit source terms
   * ==================================================
   *
   * Terms are saved in smbrk, smbrw
   * At the end of this step, we keep smbrk, smbrw, gdkgdw */

  /* Standard k-w SST RANS model */

  /* Note: Multiple loops not needing halo exchanges follow here;
           kernel fusion should be possible here. */
  {
    if (hybrid_turb == CS_HYBRID_NONE) {

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

        cs_real_t visct  = cpro_pcvto[c_id];
        cs_real_t ro     = cromo[c_id];
        cs_real_t xk     = cvara_k[c_id];
        cs_real_t xw     = cvara_omg[c_id];
        cs_real_t xxf1   = xf1[c_id];
        cs_real_t xgamma = xxf1 * ckwgm1 + (1. - xxf1)*ckwgm2;
        cs_real_t xbeta  = xxf1 * ckwbt1 + (1. - xxf1)*ckwbt2;

        cs_real_t prodw_dmut;
        if (visct <= cs_math_epzero*prodw[c_id])
          prodw_dmut = 0.;
        else
          prodw_dmut = prodw[c_id]/visct;

        smbrk[c_id] += cell_f_vol[c_id]*(prodk[c_id] - cmu*ro*xw*xk);

        smbrw[c_id] += cell_f_vol[c_id]*(  ro*xgamma*prodw_dmut
                                         - xbeta*ro*cs_math_pow2(xw)
                                         + 2.*ro/xw*(1.-xxf1)/ckwsw2
                                           *gdkgdw[c_id]);

        tinstw[c_id] +=   cell_f_vol[c_id]
                        * fmax(  -2.*ro/cs_math_pow2(xw)*(1. - xxf1)
                               / ckwsw2*gdkgdw[c_id], 0.);
      });
    }

    /* DES or DDES mode for k-w SST */

    else if (   hybrid_turb == CS_HYBRID_DES
             || hybrid_turb == CS_HYBRID_DDES) {

      const cs_real_t cddes = cs_turb_cddes;

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

        cs_real_t xfd = 0;

        cs_real_t visct  = cpro_pcvto[c_id];
        cs_real_t xnu    = cpro_pcvlo[c_id];
        cs_real_t ro    = cromo[c_id];
        cs_real_t xk     = cvara_k[c_id];
        cs_real_t xw     = cvara_omg[c_id];
        cs_real_t xxf1   = xf1[c_id];
        cs_real_t xgamma = xxf1 * ckwgm1 + (1. - xxf1)*ckwgm2;
        cs_real_t xbeta  = xxf1 * ckwbt1 + (1. - xxf1)*ckwbt2;
        cs_real_t xdist  = CS_MAX(w_dist[c_id], cs_math_epzero);
        cs_real_t xlt    = sqrt(xk) / (cmu*xw);
        cs_real_t xdelta = pow(cell_f_vol[c_id], 1./3.);

        if (hybrid_turb == CS_HYBRID_DES)
          /* Detached Eddy Simulation - DES - mode */
          xfd = 1.;
        else if (hybrid_turb == CS_HYBRID_DDES) {
          /* Delayed Detached Eddy Simulation - DDES - mode */
          cs_real_t xs2pw2 = fmax(sqrt(s2pw2[c_id]), cs_math_epzero);
          cs_real_t xrd =   (visct + xnu)
                          / (ro * xs2pw2 * cs_math_pow2(xkappa)
                                         * cs_math_pow2(xdist));
          xfd    = 1. - tanh(pow(8.*xrd, 3));
        }

        cs_real_t xdiff  = CS_MAX(xlt - (cddes * xdelta), 0.);
        cs_real_t fhybr  = xlt/(xlt - xfd*xdiff);
        /* fhybr is stored so that it can be used after */
        w1[c_id]= fhybr;
        hybrid_fd_coeff[c_id] = xfd;

        /* Storage of the Fd coefficient and check if
           DDES is activated for post-processing
           NB: for RANS zones (L_RANS < L_LES) Fd is clipped to 0 */
        if (cddes*xdelta >= xlt)
          hybrid_fd_coeff[c_id] = 0.;

        cs_real_t prodw_dmut;
        if (visct <= cs_math_epzero*prodw[c_id])
          prodw_dmut = 0.;
        else
          prodw_dmut = prodw[c_id]/visct;

        smbrk[c_id] += cell_f_vol[c_id]*(prodk[c_id] - cmu*ro*xw*xk*fhybr);

        smbrw[c_id] += cell_f_vol[c_id]*(  ro*xgamma*prodw_dmut
                                         - xbeta*ro*cs_math_pow2(xw)
                                         + 2.*ro/xw*(1. - xxf1)
                                           /ckwsw2*gdkgdw[c_id]);

        tinstw[c_id] +=   cell_f_vol[c_id]
                        * fmax(-2.*ro/cs_math_pow2(xw)*(1.-xxf1)
                               /ckwsw2*gdkgdw[c_id], 0.);
      });

    }

    /* Scale Adaptive mode for k-w SST */

    else if (hybrid_turb == CS_HYBRID_SAS) {

      const cs_real_t csas = cs_turb_csas;
      const cs_real_t csas_eta2 = cs_turb_csas_eta2;

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

        cs_real_t visct  = cpro_pcvto[c_id];
        cs_real_t ro     = cromo[c_id];
        cs_real_t xk     = cvara_k[c_id];
        cs_real_t xw     = cvara_omg[c_id];
        cs_real_t xxf1   = xf1[c_id];
        cs_real_t xgamma = xxf1 * ckwgm1 + (1. - xxf1)*ckwgm2;
        cs_real_t xbeta  = xxf1 * ckwbt1 + (1. - xxf1)*ckwbt2;

        /* Computation of the Scale Adaptive model source term xqsas
           that is added to the omega equation */
        cs_real_t xsij2 = cpro_s2kw[c_id];
        cs_real_t lvkmin =   csas
                           * sqrt(  csas_eta2 * xkappa
                                  / ((xbeta/cmu)-xgamma))
                           * pow(volume[c_id], 1./3.);
        cs_real_t lvk = xkappa*sqrt(xsij2 / d2uidxi2[c_id]);
        cs_real_t lvksas = CS_MAX(lvkmin, lvk);
        cs_real_t lmod = sqrt(xk)/(pow(cmu,0.25)*xw);

        cs_real_t xqsas = fmax(  ro * csas_eta2 * xkappa * xsij2
                               * cs_math_pow2(lmod/lvksas)
                               - 6.*ro*xk*maxgdsv[c_id], 0.);

        cs_real_t prodw_dmut;
        if (visct <= cs_math_epzero*prodw[c_id])
          prodw_dmut = 0.;
        else
          prodw_dmut = prodw[c_id]/visct;

        sas_source_term[c_id] = xqsas;

        smbrk[c_id] += cell_f_vol[c_id]*(prodk[c_id] - cmu*ro*xw*xk);

        smbrw[c_id] += cell_f_vol[c_id]*(ro*xgamma*prodw_dmut
                                         - xbeta*ro*cs_math_pow2(xw)
                                         + 2.*ro/xw*(1.-xxf1)
                                                /ckwsw2*gdkgdw[c_id]
                                         + xqsas);

        tinstw[c_id] +=   cell_f_vol[c_id]
                        * fmax(-2.*ro/cs_math_pow2(xw)*(1. - xxf1)
                               /ckwsw2*gdkgdw[c_id], 0.);
      });
    }

    /* HTLES mode for k-w SST */

    else if (hybrid_turb == CS_HYBRID_HTLES) {

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

        cs_real_t visct  = cpro_pcvto[c_id];
        cs_real_t ro     = cromo[c_id];
        cs_real_t xk     = cvara_k[c_id];
        cs_real_t xw     = cvara_omg[c_id];
        cs_real_t xxf1   = xf1[c_id];
        cs_real_t xgamma = xxf1 * ckwgm1 + (1. - xxf1)*ckwgm2;
        cs_real_t xbeta  = xxf1 * ckwbt1 + (1. - xxf1)*ckwbt2;

        /* HTLES method */
        cs_real_t xpsi  = psi[c_id];
        cs_real_t xt    = htles_t[c_id];
        cs_real_t xep   = xk/xt;

        cs_real_t prodw_dmut;
        if (visct <= cs_math_epzero*prodw[c_id])
          prodw_dmut = 0.;
        else
          prodw_dmut = prodw[c_id]/visct;

        smbrk[c_id] += cell_f_vol[c_id]*(prodk[c_id] - ro*xep);

        smbrw[c_id] += cell_f_vol[c_id]*(  ro*xgamma*prodw_dmut/xpsi
                                         - xbeta*ro*cs_math_pow2(xw)
                                         + 2.*ro/xpsi/xw*(1.-xxf1)
                                           /ckwsw2*gdkgdw[c_id]);

        tinstw[c_id] +=   cell_f_vol[c_id]
                        * fmax(  -2.*ro/xpsi/cs_math_pow2(xw)*(1. - xxf1)
                               / ckwsw2*gdkgdw[c_id], 0.);
      });

    } /* End of test on hybrid models */

    /* If the solving of k-omega is uncoupled,
       negative source terms are implicited */

    if (cs_glob_turb_rans_model->ikecou == 0) {

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        cs_real_t xw    = cvara_omg[c_id];
        cs_real_t xxf1  = xf1[c_id];
        cs_real_t xbeta = xxf1 * ckwbt1 + (1. - xxf1)*ckwbt2;
        cs_real_t ro = crom[c_id];

        cs_real_t fhybr;
        if (   hybrid_turb == CS_HYBRID_DES
            || hybrid_turb == CS_HYBRID_DDES)
          fhybr = w1[c_id];
        else if (hybrid_turb == CS_HYBRID_HTLES)
          fhybr = 1./(cmu*xw*htles_t[c_id]);
        else
          fhybr = 1.;

        tinstk[c_id] += cell_f_vol[c_id]*cmu*ro*xw*fhybr;
        tinstw[c_id] += 2.*cell_f_vol[c_id]*xbeta*ro*xw;
      });
    }

  } /* End of potential kernel fusion section */

  ctx.wait();

  /* Free memory */
  CS_FREE_HD(gdkgdw);

  /* Account for Lagrangian 2-way coupling source terms
     -------------------------------------------------- */

  /* 2nd order not handled */

  if (cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING) {

    const cs_lagr_source_terms_t  *lag_st = cs_glob_lagr_source_terms;

    const cs_real_t ce4 = cs_turb_ce4;

    if (lag_st->ltsdyn == 1) {

      cs_real_t *lag_st_k = cs_field_by_name("lagr_st_k")->val;
      cs_real_t *lag_st_i = cs_field_by_name("lagr_st_imp_velocity")->val;

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

        /* Implicit and explicit source terms on k */
        cs_real_t st_k = cell_f_vol[c_id] * lag_st_k[c_id];
        cs_real_t st_i = cell_f_vol[c_id] * lag_st_i[c_id];
        smbrk[c_id] += st_k;

        /* Explicit source terms on omega: directly use CE4 constant from
           k-epsilon (without justification). To be explored if necessary */

        smbrw[c_id] +=   ce4 * st_k * cromo[c_id]
                       / cpro_pcvto[c_id];

        /* Implicit source terms on k,
         * Note: we implicit lag_st_k if negative */
        tinstk[c_id] += fmax(-st_k/cvara_k[c_id], 0.);
        tinstk[c_id] += CS_MAX(-st_i, 0.);

        /* Implicit source terms on omega */
        tinstw[c_id] += CS_MAX(-ce4 * st_k / cvara_k[c_id], 0.);

      });

    }
  }

  ctx.wait();

  /* Mass source terms (Implicit and explicit parts)
     =============================================== */

  if (   eqp_k->n_volume_mass_injections > 0
      || eqp_w->n_volume_mass_injections > 0) {

    int *mst_type = nullptr;
    cs_lnum_t n_elts = 0;
    const cs_lnum_t *elt_ids = nullptr;
    cs_real_t *mst_val = nullptr, *mst_val_p = nullptr;
    cs_real_t *gapinj_k = nullptr, *gapinj_omg = nullptr;

    /* If we extrapolate the source terms we put Gamma Pinj in c_st */
    if (istprv >= 0) {
      gapinj_k = c_st_k_p;
      gapinj_omg = c_st_omg_p;
    }
    /* Otherwise we put it directly in smbr */
    else {
      gapinj_k = smbrk;
      gapinj_omg = smbrw;
    }

    /* We increment SMBRS by -Gamma.var_prev and ROVSDT by Gamma */

    /* For k */
    cs_volume_mass_injection_get_arrays(f_k,
                                        &n_elts,
                                        &elt_ids,
                                        &mst_type,
                                        &mst_val,
                                        &mst_val_p);
    cs_mass_source_terms(1,
                         1,
                         n_elts,
                         elt_ids,
                         mst_type,
                         cell_f_vol,
                         cvara_k,
                         mst_val,
                         mst_val_p,
                         smbrk,
                         tinstk,
                         gapinj_k);

    /* For omg */
    cs_volume_mass_injection_get_arrays(f_omg,
                                        &n_elts,
                                        &elt_ids,
                                        &mst_type,
                                        &mst_val,
                                        &mst_val_p);
    cs_mass_source_terms(1,
                         1,
                         n_elts,
                         elt_ids,
                         mst_type,
                         cell_f_vol,
                         cvara_omg,
                         mst_val,
                         mst_val_p,
                         smbrw,
                         tinstw,
                         gapinj_omg);

  }

  /* Finalize source terms */
  if (istprv >= 0) {
    cs_real_t thetp1 = 1. + thets;
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      smbrk[c_id] += thetp1 * c_st_k_p[c_id];
      smbrw[c_id] += thetp1 * c_st_omg_p[c_id];
    });
  }

  /* Re-set Boundary conditions flux coefficient for k and omega
     The definition of cofaf requires hint=(mu+muT/sigma)/distb where sigma
     is not constant in the k-omega model (and not directly accessible)
     *===================================================================== */

  const cs_real_t ckwsk1 = cs_turb_ckwsk1;
  const cs_real_t ckwsk2 = cs_turb_ckwsk2;
  const cs_real_t ckwsw1 = cs_turb_ckwsw1;

  ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {

    cs_lnum_t c_id = b_face_cells[face_id];

    /* Physical Properties */
    cs_real_t visclc = viscl[c_id];
    cs_real_t visctc = cvisct[c_id];

    cs_real_t xxf1 = xf1[c_id];

    /* k */
    cs_real_t sigma = xxf1*ckwsk1 + (1. - xxf1)*ckwsk2;
    cs_real_t hint = (visclc + visctc/sigma)/distb[face_id];

    /* Translate coefa into cofaf and coefb into cofbf */
    coefaf_k[face_id] = -hint*coefa_k[face_id];
    coefbf_k[face_id] = hint*(1. - coefb_k[face_id]);

    /* Omega */
    sigma = xxf1*ckwsw1 + (1. - xxf1)*ckwsw2;
    hint = (visclc + visctc/sigma)/distb[face_id];

    /* Translate coefa into cofaf and coefb into cofbf */
    coefaf_o[face_id] = -hint*coefa_o[face_id];
    coefbf_o[face_id] = hint*(1. - coefb_o[face_id]);

  });

  ctx.wait();

  /* Account for convection/diffusion in the right-hand side for the reinforced
   * k-omega coupling (ikecou == 1)
   * ========================================================================== */

  /* Terms are saved in w5 AND w6, THEN ADDED tosmbrk, smbrw
   * At the end of this step, we keep w3-6, smbrk, smbrw, usimpk */

  cs_real_t *w5 = nullptr;
  cs_real_t *w6 = nullptr;
  cs_real_t *w7 = nullptr;

  if (cs_glob_turb_rans_model->ikecou == 1) {

    CS_MALLOC_HD(w5, n_cells_ext, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(w6, n_cells_ext, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(w7, n_cells_ext, cs_real_t, cs_alloc_mode);

    cs_arrays_set_value<cs_real_t, 1>(n_cells, 0., w5, w6);

    /* Handle k */

    if (eqp_k->idiff >=  1) {

      const cs_lnum_t k_idifft = eqp_k->idifft;

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        cs_real_t xxf1 = xf1[c_id];
        cs_real_t sigma = xxf1*ckwsk1 + (1. - xxf1)*ckwsk2;
        w7[c_id] = viscl[c_id] + k_idifft*cvisct[c_id]/sigma;
      });

      ctx.wait();

      cs_face_viscosity(m,
                        fvq,
                        eqp_k->imvisf,
                        w7,
                        viscf,
                        viscb);

    }
    else {
      cs_arrays_set_value<cs_real_t, 1>(n_i_faces, 0., viscf);
      cs_arrays_set_value<cs_real_t, 1>(n_b_faces, 0., viscb);
    }

    cs_balance_scalar(cs_glob_time_step_options->idtvar,
                      f_k->id,
                      0,     /* imucpp */
                      1,     /* imasac */
                      1,     /* inc */
                      &eqp_k_loc,
                      cvara_k,
                      cvara_k,
                      f_k->bc_coeffs,
                      i_massflux,
                      b_massflux,
                      viscf,
                      viscb,
                      nullptr,
                      nullptr,
                      nullptr,
                      nullptr,
                      0, /* boundary convective upwind flux */
                      nullptr,
                      w5);

    if (eqp_k->verbosity >= 2) {
      cs_log_printf(CS_LOG_DEFAULT,
                    " Variable %s: EXPLICIT BALANCE =  %12.5e\n",
                    cs_field_get_label(f_k),
                    sqrt(cs_gdot(n_cells, smbrk, smbrk)));
    }

    /* Handle omega */

    if (eqp_w->idiff >= 1) {

      const cs_lnum_t w_idifft = eqp_w->idifft;

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        cs_real_t xxf1 = xf1[c_id];
        cs_real_t sigma = xxf1*ckwsw1 + (1.-xxf1)*ckwsw2;
        w7[c_id] = viscl[c_id] + w_idifft*cvisct[c_id]/sigma;
      });

      ctx.wait();

      cs_face_viscosity(m,
                        fvq,
                        eqp_w->imvisf,
                        w7,
                        viscf,
                        viscb);

    }
    else {
      cs_arrays_set_value<cs_real_t, 1>(n_i_faces, 0., viscf);
      cs_arrays_set_value<cs_real_t, 1>(n_b_faces, 0., viscb);
    }

    cs_balance_scalar(cs_glob_time_step_options->idtvar,
                      f_omg->id,
                      0,     /* imucpp */
                      1,     /* imasac */
                      1,     /* inc */
                      &eqp_w_loc,
                      cvara_omg ,
                      cvara_omg ,
                      f_omg->bc_coeffs,
                      i_massflux,
                      b_massflux,
                      viscf,
                      viscb,
                      nullptr,
                      nullptr,
                      nullptr,
                      nullptr,
                      0, /* boundary convective upwind flux */
                      nullptr,
                      w6);

    if (eqp_w->verbosity >= 2) {
      cs_log_printf(CS_LOG_DEFAULT,
                    " Variable %s: EXPLICIT BALANCE =  %12.5e\n",
                    cs_field_get_label(f_omg),
                    sqrt(cs_gdot(n_cells, smbrw, smbrw)));
    }

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      smbrk[c_id] += w5[c_id];
      smbrw[c_id] += w6[c_id];
    });
  }

  /* k-omega coupling (ikecou == 1)
   * ============================== */

  /* Ordre 2 not handled */
  if (cs_glob_turb_rans_model->ikecou == 1) {

    /* Take into account, if necessary, the Spalart-Shur rotation/curvature
       correction of the production term */
    if (cs_glob_turb_rans_model->irccor == 1) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        w1[c_id] = rotfct[c_id];
      });
    } else {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        w1[c_id] = 1.;
      });
    }

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

      cs_real_t ro = crom[c_id];

      /* Coupled solution */

      cs_real_t romvsd     = 1./(ro*volume[c_id]);
      smbrk[c_id] = smbrk[c_id]*romvsd;
      smbrw[c_id] = smbrw[c_id]*romvsd;
      cs_real_t divp23 = d2s3*cs_math_fmax(cpro_divukw[c_id], 0.);
      cs_real_t produc = w1[c_id]*cpro_s2kw[c_id] + grad_dot_g[c_id];
      cs_real_t xk = cvara_k[c_id];
      cs_real_t xw = cvara_omg[c_id];
      cs_real_t xxf1 = xf1[c_id];
      cs_real_t xgamma = xxf1*ckwgm1 + (1. - xxf1)*ckwgm2;
      cs_real_t xbeta = xxf1*ckwbt1 + (1. - xxf1)*ckwbt2;

      cs_real_t a11 =   1./dt[c_id] - 1./xw*fmin(produc, 0.)
                      + divp23+cmu*xw;
      cs_real_t a12 = cmu*xk;
      cs_real_t a21 = 0.;
      cs_real_t a22 = 1./dt[c_id] + xgamma*divp23 + 2.*xbeta*xw;

      cs_real_t unsdet = 1./(a11*a22 - a12*a21);

      cs_real_t deltk = (a22*smbrk[c_id] - a12*smbrw[c_id])*unsdet;
      cs_real_t deltw = (-a21*smbrk[c_id] + a11*smbrw[c_id])*unsdet;

      /* New source term for cs_convection_diffusion_solve */

      romvsd = ro * cell_f_vol[c_id] / dt[c_id];

      smbrk[c_id] = romvsd*deltk;
      smbrw[c_id] = romvsd*deltw;

    });

    /* Substract convection/diffusion at time n from smbrk and  smbrw
       if they were computed */
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      smbrk[c_id] = smbrk[c_id] - w5[c_id];
      smbrw[c_id] = smbrw[c_id] - w6[c_id];
    });

    ctx.wait();

    /* Free memory */
    CS_FREE_HD(w5);
    CS_FREE_HD(w6);
    CS_FREE_HD(w7);
  }

  /* Solve for turbulent kinetic energy (k)
     ---------------------------------- */

  /* Face viscosity */
  if (eqp_k->idiff >= 1) {

    const cs_lnum_t k_idifft = eqp_k->idifft;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cs_real_t xxf1 = xf1[c_id];
      cs_real_t sigma = xxf1*ckwsk1 + (1.-xxf1)*ckwsk2;
      w1[c_id] = viscl[c_id] + k_idifft*cvisct[c_id]/sigma;
    });

    ctx.wait();

    cs_face_viscosity(m,
                      fvq,
                      eqp_k->imvisf,
                      w1,
                      viscf,
                      viscb);

  } else {

    cs_arrays_set_value<cs_real_t, 1>(n_i_faces, 0., viscf);
    cs_arrays_set_value<cs_real_t, 1>(n_b_faces, 0., viscb);

  }

  /* Solve k */

  cs_equation_iterative_solve_scalar(cs_glob_time_step_options->idtvar,
                                     1, /* init */
                                     f_k->id,
                                     f_k->name,
                                     0, /* iescap */
                                     0, /* imucpp */
                                     -1, /* normp */
                                     &eqp_k_loc,
                                     cvara_k,
                                     cvara_k,
                                     f_k->bc_coeffs,
                                     i_massflux,
                                     b_massflux,
                                     viscf,
                                     viscb,
                                     viscf,
                                     viscb,
                                     nullptr,
                                     nullptr,
                                     nullptr,
                                     0, /* boundary convective upwind flux */
                                     nullptr,
                                     tinstk,
                                     smbrk,
                                     cvar_k,
                                     dpvar,
                                     nullptr,
                                     nullptr);

  /* Solve for Omega */

  /* Face viscosity */
  if (eqp_w->idiff >=  1) {

    const cs_lnum_t w_idifft = eqp_w->idifft;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cs_real_t xxf1 = xf1[c_id];
      cs_real_t sigma = xxf1*ckwsw1 + (1. - xxf1)*ckwsw2;
      w1[c_id] = viscl[c_id] + w_idifft*cvisct[c_id]/sigma;
    });

    ctx.wait();

    cs_face_viscosity(m,
                      fvq,
                      eqp_w->imvisf,
                      w1,
                      viscf,
                      viscb);
  }
  else {

    cs_arrays_set_value<cs_real_t, 1>(n_i_faces, 0., viscf);
    cs_arrays_set_value<cs_real_t, 1>(n_b_faces, 0., viscb);

  }

  /* Solve omega */

  cs_equation_iterative_solve_scalar(cs_glob_time_step_options->idtvar,
                                     1,    /* init */
                                     f_omg->id,
                                     f_omg->name,
                                     0,   /* iescap */
                                     0,   /* imucpp */
                                     -1,  /* normp */
                                     &eqp_w_loc,
                                     cvara_omg,
                                     cvara_omg,
                                     f_omg->bc_coeffs,
                                     i_massflux,
                                     b_massflux,
                                     viscf,
                                     viscb,
                                     viscf,
                                     viscb,
                                     nullptr,
                                     nullptr,
                                     nullptr,
                                     0, /* boundary convective upwind flux */
                                     nullptr,
                                     tinstw,
                                     smbrw,
                                     cvar_omg,
                                     dpvar,
                                     nullptr,
                                     nullptr);

  /* Clip values
     ============ */

  /* Compute Min/Max before clipping, for logging */

  cs_real_t *cvar_var = nullptr;
  cs_real_t vrmax[2];
  cs_real_t vrmin[2];

  const double l_threshold = 1.e12;
  for (int ii = 0; ii < 2; ii++) {
    if (ii == 0)
      cvar_var = cvar_k;
    else if (ii == 1)
      cvar_var = cvar_omg;

    vrmin[ii] =  l_threshold;
    vrmax[ii] = -l_threshold;
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      cs_real_t var = cvar_var[c_id];
      vrmin[ii] = CS_MIN(vrmin[ii], var);
      vrmax[ii] = CS_MAX(vrmax[ii], var);
    }
  }

  /* Simply clip k and omega by absolute value */

  int kclipp = cs_field_key_id("clipping_id");

  int clip_k_id = cs_field_get_key_int(f_k, kclipp);
  cs_real_t *cpro_k_clipped = nullptr;
  if (clip_k_id >= 0) {
    cpro_k_clipped = cs_field_by_id(clip_k_id)->val;
  }

  int clip_w_id =  cs_field_get_key_int(f_omg, kclipp);
  cs_real_t *cpro_w_clipped = nullptr;
  if (clip_w_id >= 0) {
    cpro_w_clipped = cs_field_by_id(clip_w_id)->val;
  }

  if (clip_k_id >= 0)
    cs_arrays_set_value<cs_real_t, 1>(n_cells, 0., cpro_k_clipped);
  if (clip_w_id >= 0)
    cs_arrays_set_value<cs_real_t, 1>(n_cells, 0., cpro_w_clipped);

  cs_lnum_t iclipk[1] = {0};
  cs_lnum_t iclipw = 0;
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_real_t xk = cvar_k[c_id];
    cs_real_t xw = cvar_omg[c_id];
    if (fabs(xk) <= epz2) {
      iclipk[0] = iclipk[0] + 1;
      if (clip_k_id >= 0)
        cpro_k_clipped[c_id] = epz2 - xk;
      cvar_k[c_id] = epz2;
    } else if (xk <= 0.) {
      iclipk[0] = iclipk[0] + 1;
      if (clip_k_id >= 0)
        cpro_k_clipped[c_id] = - xk;
      cvar_k[c_id] = -xk;
    }
    if (fabs(xw) <= epz2) {
      iclipw = iclipw + 1;
      if (clip_w_id >= 0)
        cpro_w_clipped[c_id] = epz2 - xw;
      cvar_omg[c_id] = epz2;
    } else if  (xw <= 0.) {
      iclipw = iclipw + 1;
      if (clip_w_id >= 0)
        cpro_w_clipped[c_id] = - xw;
      cvar_omg[c_id] = -xw;
    }
  }

  /* Save number of clippings for log */
  cs_lnum_t iclpkmx[1] = {0};
  int id;
  cs_real_t v_min[1];
  cs_real_t v_max[1];
  cs_lnum_t iclip;
  for (int ii = 0; ii < 2; ii++ ) {
    if (ii == 0) {
      id = f_k->id;
      v_min[0] = vrmin[ii];
      v_max[0] = vrmax[ii];
      iclip = iclipk[0];
    } else if (ii == 1) {
      id = f_omg->id;
      v_min[0] = vrmin[ii];
      v_max[0] = vrmax[ii];
      iclip = iclipw;
    }

    cs_log_iteration_clipping_field(id,
                                    iclip,
                                    0,
                                    v_min,
                                    v_max,
                                    iclipk,
                                    iclpkmx);
  }

  /* Advanced reinitialization
     ========================= */

  /* Automatic reinitialization at the end of the first iteration:
     wall distance y^+ is computed with -C log(1-alpha), where C=CL*Ceta*L*kappa,
     then y so we have an idea of the wall distance in complex geometries.
     Then U is initialized with a Reichard layer,
     Epsilon by 1/(kappa y), clipped next to the wall at its value for y^+=15
     k is given by a blending between eps/(2 nu)*y^2 and utau/sqrt(Cs_Turb_Cmu)
     The blending function is chosen so that the asymptotic behavior
     and give the correct peak of k (not the same blending than for the EBRSM
     because k profile is not the same for k-omega).
     For omega, far from the wall we take eps/Cs_Turb_Cmu/k, but next to the wall,
     omega solution is enforced. */

  /*TODO FIXME: why not just before? Are the BCs uncompatible? */

  if (   cs_glob_time_step->nt_cur == 1
      && cs_glob_turb_rans_model->reinit_turb == 1) {

    cs_real_3_t *vel = (cs_real_3_t *)f_vel->val;

    cs_real_t utaurf = 0.05*uref;
    cs_real_t nu0 = viscl0 / ro0;

    const cs_lnum_t *c_disable_flag = nullptr;
    if (fvq->has_disable_flag)
      c_disable_flag = fvq->c_disable_flag;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      /* In case of coupled fluid/solid simulation
       * avoid reinitializing turbulence in the solid region */
      if (c_disable_flag != nullptr) {
        if (c_disable_flag[c_id])
          return; // return from lambda == continue loop
      }

      /* Compute the velocity magnitude */
      cs_real_t xunorm = cs_math_3_norm(vel[c_id]);

      cs_real_t ya = w_dist[c_id];
      cs_real_t ypa = ya*utaurf/nu0;
      /* Velocity magnitude is imposed (limited only),
         the direction is conserved */
      cs_real_t limiter;
      if (xunorm <= 1.e-12*uref) {
        limiter = 1.;
      }
      else {
        limiter = fmin(utaurf/xunorm * (2.5*log(1. + 0.4*ypa)
                                        + 7.8*(1. - exp(-ypa/11.)
                                               - (ypa/11.)*exp(-0.33*ypa))),
                       1.);
      }

      vel[c_id][0] = limiter*vel[c_id][0];
      vel[c_id][1] = limiter*vel[c_id][1];
      vel[c_id][2] = limiter*vel[c_id][2];

      cs_real_t ut2 = 0.05 * uref;
      cs_real_t xeps =   cs_math_pow3(utaurf)
                       * fmin(1. / (xkappa * 15. * nu0 / utaurf),
                              1. / (xkappa * ya));
      cvar_k[c_id] = xeps /2. / nu0 * cs_math_pow2(ya)
                                    * cs_math_pow2(exp(-ypa/25.))
                    + cs_math_pow2(ut2) / sqrt(cmu)
                    * cs_math_pow2(1. - exp(-ypa / 25.));

      cvar_omg[c_id] = cs_math_pow3(ut2) / (xkappa * 15. * nu0 / ut2)
                                         / (cs_math_pow2(ut2)/sqrt(cmu))
                                         / cmu;
    });

  } /* End of test on turbulence reinitialization */

  ctx.wait();

  /* Cleanup */

  CS_FREE_HD(viscf);
  CS_FREE_HD(viscb);
  CS_FREE_HD(smbrk);
  CS_FREE_HD(smbrw);
  CS_FREE_HD(tinstk);
  CS_FREE_HD(tinstw);
  CS_FREE_HD(dpvar);
  CS_FREE_HD(w1);
  CS_FREE_HD(rotfct);
  CS_FREE_HD(xf1);

  CS_FREE_HD(grad_dot_g);
  CS_FREE_HD(usimpk);
  CS_FREE_HD(usimpw);

  CS_FREE_HD(prodk);
  CS_FREE_HD(prodw);

  CS_FREE_HD(s2pw2);
  CS_FREE_HD(maxgdsv);
  CS_FREE_HD(d2uidxi2);
}

/*----------------------------------------------------------------------------*/
/*! \brief Calculation of turbulent viscosity for
 *         the \f$ k - \omega \f$ SST model.
 *
 * \f[ \mu_T = \rho A1 \dfrac{k}{\max(A1 \omega; \; S f_2)} \f]
 * with
 * \f[ S = \sqrt{  2 S_{ij} S_{ij}} \f]
 * \f[ S_{ij} = \dfrac{\der{u_i}{x_j} + \der{u_j}{x_i}}{2}\f]
 *
 * and \f$ f_2 = \tanh(arg2^2) \f$
 * \f[ arg2^2 = \max(2 \dfrac{\sqrt{k}}{C_\mu \omega y}; \;
 *                   500 \dfrac{\nu}{\omega y^2}) \f]
 * where \f$ y \f$ is the distance to the wall.
 *
 * \f$ \divs{\vect{u}} \f$ is calculated at the same time than \f$ S \f$
 * for use in cs_turbulence_kw.
 *
 * \param[in]     phase_id      turbulent phase id (-1 for single phase flow)
 *
 !*/
/*----------------------------------------------------------------------------*/

void
cs_turbulence_kw_mu_t(int phase_id)
{
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_lnum_t n_cells = mesh->n_cells;
  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t ntcabs = cs_glob_time_step->nt_cur;

  cs_dispatch_context ctx;

  /* Initialization
     ============== */

  cs_field_t *f_k = CS_F_(k);
  cs_field_t *f_omg = CS_F_(omg);
  cs_field_t *f_vel = CS_F_(vel);
  cs_field_t *f_mu = CS_F_(mu);
  cs_field_t *f_mut = CS_F_(mu_t);
  cs_field_t *f_rho = CS_F_(rho);

  if (phase_id >= 0) {
    f_k = CS_FI_(k, phase_id);
    f_omg = CS_FI_(omg, phase_id);
    f_vel = CS_FI_(vel, phase_id);
    f_mu = CS_FI_(mu, phase_id);
    f_mut = CS_FI_(mu_t, phase_id);
    f_rho = CS_FI_(rho, phase_id);
  }

  cs_real_t *visct = f_mut->val;

  const cs_real_t *viscl = f_mu->val;
  const cs_real_t *crom  = f_rho->val;
  const cs_real_t *cvar_k  = (const cs_real_t *)f_k->val;
  const cs_real_t *cvar_omg = (const cs_real_t *)f_omg->val;

  const cs_real_t *w_dist = cs_field_by_name("wall_distance")->val;

  /* Compute the scalar s2kw rate SijSij and the trace of the velocity
   * gradient
   *   (Sij^D) (Sij^D)  is stored in    s2kw (deviatoric s2kw tensor rate)
   *   tr(Grad u)       is stored in    divukw
   * ======================================================================= */

  cs_real_33_t *gradv = nullptr, *_gradv = nullptr;
  {
    cs_field_t *f_vg = cs_field_by_name_try("algo:velocity_gradient");

    if (f_vel->grad != nullptr)
      gradv = (cs_real_33_t *)f_vel->grad;
    else if (f_vg != nullptr)
      gradv = (cs_real_33_t *)f_vg->val;
    else {
      CS_MALLOC_HD(_gradv, n_cells_ext, cs_real_33_t, cs_alloc_mode);
      gradv = _gradv;
    }
  }

  /* Computation of the velocity gradient */
  if (f_vel->grad == nullptr)
    cs_field_gradient_vector(f_vel,
                             true,  /* use_previous_t */
                             1,     /* inc */
                             gradv);

  /* s2kw = Strain rate of the deviatoric part of the s2kw tensor
   *      = 2 (Sij^D).(Sij^D)
   * divukw   = trace of the velocity gradient
   *          = dudx + dvdy + dwdz
   * ============================================================== */

  const cs_real_t d1s3 = 1./3.;
  const cs_real_t d2s3 = 2./3.;

  cs_field_t *f_s2kw = cs_field_by_name_try("s2");
  cs_field_t *f_divukw = cs_field_by_name_try("vel_gradient_trace");

  if (phase_id >= 0) {
    char f_name[64]; /* should be much larger than needed */

    snprintf(f_name, 63, "s2_%d", phase_id + 1);
    f_name[63] = '\0';
    f_s2kw = cs_field_by_name(f_name);

    snprintf(f_name, 63, "vel_gradient_trace_%d", phase_id + 1);
    f_divukw = cs_field_by_name(f_name);
  }

  cs_real_t *cpro_s2kw = f_s2kw->val;
  cs_real_t *cpro_divukw = f_divukw->val;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
   cpro_s2kw[c_id] = 2.0 *(  cs_math_pow2(  d2s3*gradv[c_id][0][0]
                                          - d1s3*gradv[c_id][1][1]
                                          - d1s3*gradv[c_id][2][2])
                           + cs_math_pow2(- d1s3*gradv[c_id][0][0]
                                          + d2s3*gradv[c_id][1][1]
                                          - d1s3*gradv[c_id][2][2])
                           + cs_math_pow2(- d1s3*gradv[c_id][0][0]
                                          - d1s3*gradv[c_id][1][1]
                                          + d2s3*gradv[c_id][2][2]))
                     + cs_math_pow2(gradv[c_id][0][1] + gradv[c_id][1][0])
                     + cs_math_pow2(gradv[c_id][0][2] + gradv[c_id][2][0])
                     + cs_math_pow2(gradv[c_id][1][2] + gradv[c_id][2][1]);

   cpro_divukw[c_id] =   gradv[c_id][0][0]
                       + gradv[c_id][1][1]
                       + gradv[c_id][2][2];
  });

  /* Calculation of viscosity
   * ======================== */

  cs_real_t *psi   = nullptr;
  cs_real_t *blend = nullptr;
  const int hybrid_turb = cs_glob_turb_model->hybrid_turb;
  if (hybrid_turb == CS_HYBRID_HTLES) {
    psi = cs_field_by_name("htles_psi")->val;
    blend = cs_field_by_name("hybrid_blend")->val;
  }

  const cs_real_t cmu = cs_turb_cmu;
  const cs_real_t ckwa1 = cs_turb_ckwa1;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    const cs_real_t xk = cvar_k[c_id];
    cs_real_t xw = cvar_omg[c_id];
    const cs_real_t rom = crom[c_id];
    const cs_real_t xmu = viscl[c_id];

    // Wall distance
    const cs_real_t xdist = fmax(w_dist[c_id], cs_math_epzero);

    cs_real_t xf2;

    // FIXME should be a check on xw...
    if (xk > 0.0) {
      // Wall distance has no value at the first time step, we consider it as infinite
      if (ntcabs == 1) {
        xf2 = 0.0;
      }
      else {
        cs_real_t xarg2 = fmax(2.0 * sqrt(xk) / cmu / xw / xdist,
                               500.0 * xmu / rom / xw / cs_math_pow2(xdist));
        xf2 = tanh(cs_math_pow2(xarg2));
      }

      if (hybrid_turb == CS_HYBRID_HTLES) {
        xw *= psi[c_id];
        xf2 *= CS_MAX(cs_math_epzero, 1. - blend[c_id]);
      }
      visct[c_id] =   rom*ckwa1*xk
                    / fmax(ckwa1*xw, sqrt(cpro_s2kw[c_id])*xf2);
    }
    else {
      visct[c_id] = 1.e-30;
    }
  });

  ctx.wait();

  CS_FREE_HD(_gradv);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
