/*============================================================================
 * k-epsilon turbulence model.
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

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_array.h"
#include "alge/cs_balance.h"
#include "alge/cs_blas.h"
#include "alge/cs_convection_diffusion.h"
#include "base/cs_dispatch.h"
#include "cdo/cs_equation.h"
#include "base/cs_equation_iterative_solve.h"
#include "alge/cs_face_viscosity.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_pointer.h"
#include "base/cs_field_operator.h"
#include "base/cs_mem.h"
#include "alge/cs_gradient.h"
#include "lagr/cs_lagr.h"
#include "base/cs_log.h"
#include "base/cs_log_iteration.h"
#include "base/cs_mass_source_terms.h"
#include "base/cs_math.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_parall.h"
#include "base/cs_physical_constants.h"
#include "pprt/cs_physical_model.h"
#include "base/cs_porous_model.h"
#include "base/cs_prototypes.h"
#include "base/cs_reducers.h"
#include "base/cs_rotation.h"
#include "base/cs_solid_zone.h"
#include "base/cs_thermal_model.h"
#include "base/cs_time_step.h"
#include "turb/cs_turbulence_model.h"
#include "turb/cs_turbulence_rotation.h"
#include "base/cs_volume_mass_injection.h"
#include "base/cs_velocity_pressure.h"
#include "base/cs_wall_functions.h"

/* Atmospheric model headers
   (we should not need to call them here, we should be more modular */
#include "atmo/cs_at_data_assim.h"
#include "atmo/cs_atprke.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "turb/cs_turbulence_ke.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_turbulence_ke.cpp

  Solve the \f$ k - \varepsilon \f$  for incompressible flows
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of the E term of the \f$ \varepsilon\f$ equation
 *       (BL-V2/K model)
 *
 *  \param[in]     phase_id      turbulent phase id (-1 for single phase flow)
 *  \param[in,out] w1            E-term
 */
/*----------------------------------------------------------------------------*/

static void
_tsepls(int       phase_id,
        cs_real_t w1[])
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t    n_cells = m->n_cells;
  const cs_lnum_t    n_cells_ext = m->n_cells_with_ghosts;
  const cs_real_t   *volume  = fvq->cell_vol;
  const cs_lnum_t   *b_face_cells = m->b_face_cells;
  const cs_lnum_2_t *i_face_cells = m->i_face_cells;
  const cs_real_t   *restrict weight = fvq->weight;
  const cs_real_3_t *i_face_normal = (const cs_real_3_t *) fvq->i_face_normal;
  const cs_real_3_t *b_face_normal = (const cs_real_3_t *) fvq->b_face_normal;

  /* Initialization
   * ============== */

  cs_array_real_fill_zero(n_cells, w1);

  cs_real_33_t *w7;
  CS_MALLOC_HD(w7, n_cells_ext, cs_real_33_t, cs_alloc_mode);

  /* Calculation of the term d2Ui/dxkdxj*d2Ui/dxkdxj
   * ================================================ */

  cs_field_t *f_vel = CS_F_(vel);

  if (phase_id >= 0)
    f_vel = CS_FI_(vel, phase_id);

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

  /* Loop over u, v, w components:
     TODO interleave */

  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  for (cs_lnum_t i = 0; i < 3; i++) {

    ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t ii = 0; ii < 3; ii++) {
        for (cs_lnum_t jj = 0; jj < 3; jj++)
          w7[c_id][ii][jj] = 0.;
      }
    });

    ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
      cs_real_t duidxk[3], njsj[3];

      cs_lnum_t c_id0 = i_face_cells[face_id][0];
      cs_lnum_t c_id1 = i_face_cells[face_id][1];
      cs_real_t pnd = weight[face_id];

      duidxk[0] = pnd * gradv[c_id0][i][0] + (1. - pnd) * gradv[c_id1][i][0];
      duidxk[1] = pnd * gradv[c_id0][i][1] + (1. - pnd) * gradv[c_id1][i][1];
      duidxk[2] = pnd * gradv[c_id0][i][2] + (1. - pnd) * gradv[c_id1][i][2];
      njsj[0]   = i_face_normal[face_id][0];
      njsj[1]   = i_face_normal[face_id][1];
      njsj[2]   = i_face_normal[face_id][2];

      for (cs_lnum_t j = 0; j < 3; j++){
        cs_real_t c_w7_0[3], c_w7_1[3];
        for (cs_lnum_t k = 0; k < 3; k++){
          c_w7_0[k] =  duidxk[k]*njsj[j];
          c_w7_1[k] = -duidxk[k]*njsj[j];
        }
        if (c_id0 < n_cells)
          cs_dispatch_sum<3>(w7[c_id0][j], c_w7_0, i_sum_type);
        if (c_id1 < n_cells)
          cs_dispatch_sum<3>(w7[c_id1][j], c_w7_1, i_sum_type);
      }
    });

    ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
      cs_real_t duidxk[3], njsj[3];

      cs_lnum_t c_id0 = b_face_cells[face_id];

      duidxk[0] = gradv[c_id0][i][0];
      duidxk[1] = gradv[c_id0][i][1];
      duidxk[2] = gradv[c_id0][i][2];
      njsj[0]   = b_face_normal[face_id][0];
      njsj[1]   = b_face_normal[face_id][1];
      njsj[2]   = b_face_normal[face_id][2];
      for (cs_lnum_t j = 0; j < 3; j++){
        cs_real_t c_w7[3];
        for (cs_lnum_t k = 0; k < 3; k++) {
          c_w7[k] = duidxk[k]*njsj[j];
        }
        cs_dispatch_sum<3>(w7[c_id0][j], c_w7, b_sum_type);
      }
    });

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cs_real_t w_temp = 0.;
      for (cs_lnum_t k = 0; k < 3; k++) {
        for (cs_lnum_t j = 0; j < 3; j++) {
          w_temp += cs_math_pow2(w7[c_id][j][k] / volume[c_id]);
        }
      }
      w1[c_id] += w_temp;
    });

    ctx.wait();

  }

  /* Free memory */
  CS_FREE(_gradv);

  CS_FREE(w7);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve the k-epsilon equations.
 *
 * Solve the \f$ k - \varepsilon \f$  for incompressible flows
 * or slightly compressible flows for one time step.
 *
 * \param[in]     phase_id      turbulent phase id (-1 for single phase flow)
 * \param[out]    prdv2f        v2f production term
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_ke(int              phase_id,
                 cs_real_t       *prdv2f)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;
  cs_mesh_quantities_t  *mq_g = cs_glob_mesh_quantities_g;

  cs_domain_t  *domain = cs_glob_domain;

  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_real_t *volume  = mq_g->cell_vol;
  const cs_real_t *cell_f_vol = fvq->cell_vol;
  const cs_real_t *distb = fvq->b_dist;
  const cs_lnum_t *b_face_cells = m->b_face_cells;

  const cs_time_scheme_t *time_scheme = cs_glob_time_scheme;
  const cs_real_t thets  = time_scheme->thetst;

  const cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();
  cs_real_t viscl0 = phys_pro->viscl0; /* reference pressure */
  cs_real_t ro0 = phys_pro->ro0; /* reference density */
  const cs_real_t uref = cs_glob_turb_ref_values->uref;
  /*const cs_real_t grav[3] = {cs_glob_physical_constants->gravity[0],
                             cs_glob_physical_constants->gravity[1],
                             cs_glob_physical_constants->gravity[2]};*/

  const cs_real_t *grav = cs_glob_physical_constants->gravity;
#if defined(HAVE_ACCEL)
  cs_real_t *_grav = nullptr;
  if (cs_get_device_id() > -1) {
    CS_MALLOC_HD(_grav, 3, cs_real_t, cs_alloc_mode);
    for (int i = 0; i < 3; i++) {
      _grav[i] = cs_glob_physical_constants->gravity[i];
    }

    cs_mem_advise_set_read_mostly(_grav);

    grav = _grav;
  }
#endif

  cs_dispatch_context ctx;

  cs_field_t *f_k = CS_F_(k);
  cs_field_t *f_eps = CS_F_(eps);
  cs_field_t *f_phi = CS_F_(phi);
  cs_field_t *f_alpbl = CS_F_(alp_bl);

  const cs_real_t *dt = CS_F_(dt)->val;

  if (phase_id >= 0) {
    f_k = CS_FI_(k, phase_id);
    f_eps = CS_FI_(eps, phase_id);
    f_phi = CS_FI_(phi, phase_id);
    f_alpbl = CS_FI_(alp_bl, phase_id);
  }

  cs_real_t sigmak = cs_field_get_key_double(f_k,
                                             cs_field_key_id("turbulent_schmidt"));

  /* If turbulent Schmidt is variable, id of the corresponding field */
  int sigmak_id = cs_field_get_key_int(f_k,
                                       cs_field_key_id("turbulent_schmidt_id"));

  cs_real_t *cpro_sigmak = nullptr;
  if (sigmak_id >= 0)
    cpro_sigmak = cs_field_by_id(sigmak_id)->val;

  cs_real_t sigmae = cs_field_get_key_double(f_eps,
                                             cs_field_key_id("turbulent_schmidt"));

  /* If turbulent Schmidt is variable, id of the corresponding field */
  int sigmae_id = cs_field_get_key_int(f_eps,
                                       cs_field_key_id("turbulent_schmidt_id"));

  cs_real_t *cpro_sigmae = nullptr;
  if (sigmae_id >= 0)
    cpro_sigmae = cs_field_by_id(sigmae_id)->val;

  /* Initialization of work arrays in case of Hybrid turbulence modelling
   * ==================================================================== */

  cs_real_t *htles_psi       = nullptr;
  cs_real_t *htles_t         = nullptr;
  cs_real_t *hybrid_fd_coeff = nullptr;
  if (cs_glob_turb_model->hybrid_turb == 4) {
    htles_psi = cs_field_by_name("htles_psi")->val;
    htles_t   = cs_field_by_name("htles_t")->val;
    hybrid_fd_coeff = cs_field_by_name("hybrid_blend")->val;
  }

  /* Initialization
     ============== */

  const cs_turb_model_type_t model
    = (cs_turb_model_type_t)(cs_glob_turb_model->model);

  /* Allocate temporary arrays for the turbulence resolution */

  cs_real_t *viscf, *viscb;
  CS_MALLOC_HD(viscf, n_i_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(viscb, n_b_faces, cs_real_t, cs_alloc_mode);

  cs_real_t *smbrk, *smbre;
  CS_MALLOC_HD(smbrk, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(smbre, n_cells_ext, cs_real_t, cs_alloc_mode);

  cs_real_t *prodk, *gk;
  CS_MALLOC_HD(prodk, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(gk, n_cells_ext, cs_real_t, cs_alloc_mode);

  cs_real_t *tinstk, *tinste, *divu, *strain_sq;
  CS_MALLOC_HD(tinstk, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(tinste, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(divu, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(strain_sq, n_cells_ext, cs_real_t, cs_alloc_mode);

  /* Allocate work arrays */

  cs_real_t *w1, *w3, *w4;

  CS_MALLOC_HD(w1, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(w3, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(w4, n_cells_ext, cs_real_t, cs_alloc_mode);

  cs_real_t *usimpk, *usexpk, *usexpe, *usimpe, *dpvar;
  CS_MALLOC_HD(usimpk, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(usexpk, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(usexpe, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(usimpe, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(dpvar, n_cells_ext, cs_real_t, cs_alloc_mode);

  cs_real_t *prdtke = nullptr, *prdeps = nullptr;
  cs_real_t *strain = nullptr, *sqrt_k = nullptr;
  cs_real_3_t *grad_sqk = nullptr, *grad_s = nullptr;
  cs_real_t *e_term = nullptr;

  if (model == CS_TURB_K_EPSILON) {
    CS_MALLOC_HD(prdtke, n_cells_ext, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(prdeps, n_cells_ext, cs_real_t, cs_alloc_mode);
  }
  else if (model == CS_TURB_K_EPSILON_LS) {
    CS_MALLOC_HD(sqrt_k, n_cells_ext, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(strain, n_cells_ext, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(grad_sqk, n_cells_ext, cs_real_3_t, cs_alloc_mode);
    CS_MALLOC_HD(grad_s, n_cells_ext, cs_real_3_t, cs_alloc_mode);
  }
  else if (model == CS_TURB_V2F_BL_V2K) {
    CS_MALLOC_HD(e_term, n_cells_ext, cs_real_t, cs_alloc_mode);
  }

  /* Map field arrays */
  cs_field_t *f_mut = CS_F_(mu_t);
  cs_field_t *f_mu = CS_F_(mu);
  cs_field_t *f_vel = CS_F_(vel);

  if (phase_id >= 0) {
    f_mut = CS_FI_(mu_t, phase_id);
    f_mu = CS_FI_(mu, phase_id);
    f_vel = CS_FI_(vel, phase_id);
  }

  const cs_real_t *cvisct = (const cs_real_t *)f_mut->val;
  const cs_real_t *viscl = (const cs_real_t *)f_mu->val;

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  int iflmas =  cs_field_get_key_int(f_vel, kimasf);
  int iflmab =  cs_field_get_key_int(f_vel, kbmasf);

  const cs_real_t *imasfl =  cs_field_by_id(iflmas)->val;
  const cs_real_t *bmasfl =  cs_field_by_id(iflmab)->val;

  cs_real_t *cvar_k   =  f_k->val;
  cs_real_t *cvara_k  =  f_k->val_pre;
  cs_real_t *cvar_ep  =  f_eps->val;
  cs_real_t *cvara_ep =  f_eps->val_pre;
  cs_real_t *cvara_phi = nullptr;
  if (   model == CS_TURB_V2F_PHI
      || model == CS_TURB_V2F_BL_V2K){
    cvara_phi = f_phi->val_pre;
  }
  cs_real_t *cvara_al = nullptr;
  if (model == CS_TURB_V2F_BL_V2K) {
    cvara_al = f_alpbl->val_pre;
  }
  const cs_real_t *w_dist = nullptr;
  if (model == CS_TURB_K_EPSILON_QUAD) {
    w_dist =  cs_field_by_name("wall_distance")->val;
  }

  int kstprv = cs_field_key_id("source_term_prev_id");
  int istprv =  cs_field_get_key_int(f_k, kstprv);
  cs_real_t *c_st_k_p = nullptr;
  cs_real_t *c_st_eps_p = nullptr;
  if (istprv >= 0) {
    c_st_k_p = cs_field_by_id(istprv)->val;
    istprv = cs_field_get_key_int(f_eps, kstprv);
    if (istprv >= 0){
      c_st_eps_p = cs_field_by_id(istprv)->val;
    }
    if (istprv >= 0)
      istprv = 1;
  }

  cs_field_t *f_rho = CS_F_(rho);
  cs_field_t *f_rhob = CS_F_(rho_b);

  if (phase_id >= 0) {
    f_rho = CS_FI_(rho, phase_id);
    f_rhob = CS_FI_(rho_b, phase_id);
  }

  cs_real_t *crom  = f_rho->val;
  cs_real_t *cromo = f_rho->val;
  cs_real_t *bromo = f_rhob->val;

  cs_real_t *cpro_pcvto = f_mut->val;
  cs_real_t *cpro_pcvlo = f_mu->val;

  /* Time extrapolation? */

  int key_t_ext_id = cs_field_key_id("time_extrapolated");

  if (istprv >= 0) {
    int iroext = cs_field_get_key_int(f_rho, key_t_ext_id);
    if (iroext > 0) {
      cromo = f_rho->val_pre;
      bromo = f_rhob->val_pre;
    }
    int iviext =  cs_field_get_key_int(f_mu, key_t_ext_id);
    if (iviext > 0)
      cpro_pcvlo = f_mu->val_pre;
    iviext = cs_field_get_key_int(f_mut, key_t_ext_id);
    if (iviext > 0) {
      cpro_pcvto = f_mut->val_pre;
    }
  }

  const cs_equation_param_t *eqp_k
    = cs_field_get_equation_param_const(f_k);

  const cs_equation_param_t *eqp_eps
    = cs_field_get_equation_param_const(f_eps);

  if (eqp_k->verbosity >= 1) {
    if (model == CS_TURB_K_EPSILON) {
      cs_log_printf(CS_LOG_DEFAULT,
                    "\n"
                    "  ** Solving k-epsilon\n"
                    "     -----------------\n");
    }
    else if (model == CS_TURB_K_EPSILON_LIN_PROD) {
      cs_log_printf(CS_LOG_DEFAULT,
                    "\n"
                    "  ** Solving k-epsilon with linear production\n"
                    "     -----------------\n");
    }
    else if (model == CS_TURB_K_EPSILON_LS) {
      cs_log_printf(CS_LOG_DEFAULT,
                    "\n"
                    "  ** Solving k-epsilon, Launder-Sharma\n"
                    "     -----------------\n");
    }
    else if (model == CS_TURB_K_EPSILON_QUAD) {
      cs_log_printf(CS_LOG_DEFAULT,
                    "\n"
                    "  ** Solving k-epsilon, quadratic Baglietto\n"
                    "     -----------------\n");
    }
    else {
      cs_log_printf(CS_LOG_DEFAULT,
                    "\n"
                    "  ** Solving v2f (k and epsilon)\n"
                    "     ---------------------------\n");
    }
  }

  /* Production term due to buoyancy */

  /* dissip_buo_mdl:
   * 0: G taken in epsilon equation only when positive (default)
   * 1: G taken whatever its sign */
  const int dissip_buo_mdl = cs_glob_turb_rans_model->dissip_buo_mdl;

  /* For the model with linear production, sqrt(Cmu) is required */

  const int hybrid_turb = cs_glob_turb_model->hybrid_turb;

  const cs_real_t cmu = cs_turb_cmu;
  const cs_real_t sqrcmu = sqrt(cs_turb_cmu);
  const cs_real_t xkappa = cs_turb_xkappa;

  const cs_real_t ce1 = cs_turb_ce1;
  const cs_real_t ce2 = cs_turb_ce2;
  const cs_real_t ce3 = cs_turb_ce3;
  const cs_real_t cnl1 = cs_turb_cnl1;
  const cs_real_t cnl2 = cs_turb_cnl2;
  const cs_real_t cnl3 = cs_turb_cnl3;
  const cs_real_t cnl4 = cs_turb_cnl4;
  const cs_real_t cnl5 = cs_turb_cnl5;

  const cs_real_t d2s3 = 2./3.;
  const cs_real_t d1s3 = 1./3.;
  const cs_real_t d1s2 = 1./2.;

  /* Advanced reinitialization
     ========================= */

  /* Automatic reinitialization at the end of the first iteration:
     wall distance y^+ is computed with -C log(1-alpha), where C=CL*Ceta*L*kappa,
     then y so we have an idea of the wall distance in complex geometries.
     Then U is initialized with a Reichard layer,
     Epsilon by 1/(kappa y), clipped next to the wall at its value for y^+=15.
     k is given by a blending between eps/(2 nu)*y^2 and utau/sqrt(Cmu)
     The blending function is chosen so that the asymptotic behavior
     and give the correct peak of k (not the same blending than for the EBRSM
     because k profile is not the same for k-omega). */

  /*TODO FIXME: Are the BCs incompatible? */

  if (   cs_glob_time_step->nt_cur == 1
      && cs_glob_turb_rans_model->reinit_turb == 1
      && model == CS_TURB_V2F_BL_V2K) {

    cs_real_t *cvar_al = (cs_real_t *)f_alpbl->val;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      /* y+ is bounded by 400, because in the Reichard profile,
         it corresponds to saturation (u>uref) */
      cvar_al[c_id] = cs::max(fmin(cvar_al[c_id], 1. - exp(-400./50.)),
                              0.);
    });

    cs_field_current_to_previous(f_alpbl); // FIXME why here ?

    cvar_ep = f_eps->val;
    cs_real_3_t *vel = (cs_real_3_t *)f_vel->val;

    cs_real_t utaurf = 0.05*uref;
    cs_real_t nu0 = viscl0 / ro0;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

      /* Compute the velocity magnitude */
      cs_real_t xunorm = cs_math_3_norm(vel[c_id]);

      /* Compute YA, therefore alpha is given by 1-exp(-YA/(50 nu/utau))
           NB: y^+ = 50 give the best compromise */
      cs_real_t ya = -log(1. - cvar_al[c_id])*50.*nu0/utaurf;
      cs_real_t ypa = ya/(nu0/utaurf);

      /* Velocity magnitude is imposed (limited only), the direction is
         conserved */
      cs_real_t limiter = 1.;
      if (xunorm > 1.e-12*uref) {
        limiter = cs::min(utaurf/xunorm*   (2.5*log(1. + 0.4*ypa)
                                         + 7.8*(1. - exp(-ypa/11.)
                                                -(ypa/11.)*exp(-0.33*ypa))),
                          1.);
      }

      vel[c_id][0] = limiter*vel[c_id][0];
      vel[c_id][1] = limiter*vel[c_id][1];
      vel[c_id][2] = limiter*vel[c_id][2];

      cs_real_t ut2 = 0.05*uref;

      cvar_ep[c_id] =   cs_math_pow3(utaurf)
                      * fmin(1./(xkappa*15.*nu0/utaurf),
                             1./(xkappa*ya));
      cvar_k[c_id] =     cvar_ep[c_id]/2./nu0*cs_math_pow2(ya)
                       * cs_math_pow2(exp(-ypa/25.))
                     +   cs_math_pow2(ut2)/sqrt(cmu)
                       * cs_math_pow2(1.-exp(-ypa/25.));
    });
    ctx.wait();

    cs_field_current_to_previous(f_vel);
    cs_field_current_to_previous(f_k);
    cs_field_current_to_previous(f_eps); /*TODO phi ? */
  }

  /* Compute the scalar strain rate squared S2 =2SijSij and the trace of
     the velocity gradient
     (Sij^D) (Sij^D)  is stored in    strain_sq (deviatoric strain tensor rate)
     tr(Grad u)       is stored in    divu
     ======================================================================= */

  /* Allocate temporary arrays for gradients calculation */
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

  const cs_real_t k_istat = eqp_k->istat;
  const cs_real_t eps_istat = eqp_eps->istat;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    /* strain_sq = Stain rate of the deviatoric part of the strain tensor
       = 2 (Sij^D).(Sij^D)
       divu   = trace of the velocity gradient
       = dudx + dvdy + dwdz */

    strain_sq[c_id] = 2.
      * (  cs_math_pow2(  d2s3*gradv[c_id][0][0]
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

    divu[c_id] = cs_math_33_trace(gradv[c_id]);

    /* Unsteady terms (stored in tinstk and tinste)
       ============================================ */

    cs_real_t romvsd = crom[c_id] * cell_f_vol[c_id] / dt[c_id];
    tinstk[c_id] = k_istat*romvsd;
    tinste[c_id] = eps_istat*romvsd;

  });

  /* Compute the square root of the strain and the turbulent
     kinetic energy for Launder-Sharma k-epsilon source terms */
  if (model == CS_TURB_K_EPSILON_LS) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      strain[c_id] = sqrt(strain_sq[c_id]);
      sqrt_k[c_id] = sqrt(cs::abs(cvar_k[c_id]));
    });
  }

  /* Compute the first part of the production term: muT (S^D)**2
   * Going out of the step we keep strain_sq, divu,
   *============================================================
   *
   * For the Linear Production k-epsilon model,
   * the production term is assumed to be asymptotic in S and
   * not in mu_TxS**2 */

  cs_field_t *f_tke_prod = cs_field_by_name_try("algo:k_production");

  if (model == CS_TURB_K_EPSILON_LIN_PROD) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cs_real_t rho    = cromo[c_id];
      cs_real_t xs     = sqrt(strain_sq[c_id]);
      cs_real_t xk   = cvar_k[c_id];
      cs_real_t cmueta = cs::min(cmu*cvara_k[c_id]/cvara_ep[c_id]*xs,
                                 sqrcmu);
      prodk[c_id] = rho*cmueta*xs*cvara_k[c_id]
                  - d2s3*rho*xk*divu[c_id];
    });
  }
  else if (model == CS_TURB_K_EPSILON_QUAD) {

    /* Turbulent production for the quadratic Baglietto k-epsilon model  */

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

      cs_real_t xstrai[3][3], xrotac[3][3];

      cs_real_t visct = cpro_pcvto[c_id];
      cs_real_t rho  = crom[c_id];
      cs_real_t xeps  = cvar_ep[c_id];
      cs_real_t xk    = cvar_k[c_id];
      cs_real_t xttke = xk/xeps;

      /* Sij */
      xstrai[0][0] = gradv[c_id][0][0];
      xstrai[1][0] = d1s2*(gradv[c_id][0][1]+gradv[c_id][1][0]);
      xstrai[0][1] = xstrai[1][0];
      xstrai[2][0] = d1s2*(gradv[c_id][0][2]+gradv[c_id][2][0]);
      xstrai[0][2] = xstrai[2][0];
      xstrai[1][1] = gradv[c_id][1][1];
      xstrai[2][1] = d1s2*(gradv[c_id][1][2]+gradv[c_id][2][1]);
      xstrai[1][2] = xstrai[2][1];
      xstrai[2][2] = gradv[c_id][2][2];

      /* omegaij */
      xrotac[0][0] = 0.;
      xrotac[1][0] = d1s2*(gradv[c_id][0][1]-gradv[c_id][1][0]);
      xrotac[2][0] = d1s2*(gradv[c_id][0][2]-gradv[c_id][2][0]);
      xrotac[0][1] = -xrotac[1][0];
      xrotac[1][1] = 0.;
      xrotac[2][1] = d1s2*(gradv[c_id][1][2]-gradv[c_id][2][1]);
      xrotac[0][2] = -xrotac[2][0];
      xrotac[1][2] = -xrotac[2][1];
      xrotac[2][2] = 0.;

      cs_real_t skskjsji = 0.;
      cs_real_t wkskjsji = 0.;
      cs_real_t skiwjksji = 0.;
      cs_real_t wkwjksji = 0.;
      cs_real_t sijsij = 0.;
      cs_real_t wijwij = 0.;

      for (cs_lnum_t ii = 0; ii < 3; ii++) {
        for (cs_lnum_t jj = 0; jj < 3; jj++) {
          sijsij += xstrai[jj][ii]*xstrai[jj][ii];
          wijwij += xrotac[jj][ii]*xrotac[jj][ii];
          for (int kk = 0; kk < 3; kk++) {
            skskjsji  += xstrai[kk][ii]*xstrai[jj][kk]*xstrai[ii][jj];
            wkskjsji  += xrotac[kk][ii]*xstrai[jj][kk]*xstrai[ii][jj];
            skiwjksji += xstrai[ii][kk]*xrotac[kk][jj]*xstrai[ii][jj];
            wkwjksji  += xrotac[kk][ii]*xrotac[kk][jj]*xstrai[ii][jj];
          }
        }
      }

      cs_real_t xss  = xttke*sqrt(2.*sijsij);
      cs_real_t xcmu = d2s3/(3.9 + xss);

      /* Evaluating "constants" */
      cs_real_t xqc1 = cnl1/((cnl4 + cnl5*cs_math_pow3(xss))*xcmu);
      cs_real_t xqc2 = cnl2/((cnl4 + cnl5*cs_math_pow3(xss))*xcmu);
      cs_real_t xqc3 = cnl3/((cnl4 + cnl5*cs_math_pow3(xss))*xcmu);

      /* Evaluating the turbulent production */
      prodk[c_id] =   visct*strain_sq[c_id]
                    - d2s3*rho*xk*divu[c_id]
                    - 4.*xqc1*visct*xttke* (skskjsji - d1s3*sijsij*divu[c_id])
                    - 4.*xqc2*visct*xttke* (wkskjsji + skiwjksji)
                    - 4.*xqc3*visct*xttke* (wkwjksji - d1s3*wijwij*divu[c_id]);

    }); /* End loop on cells */

    /* End test on specific k-epsilon model
       In the general case Pk = -RijSij = mu_t*Sij Sij - rho k divu */
  }
  else {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cs_real_t rho  = crom[c_id];
      cs_real_t xk   = cvar_k[c_id];
      prodk[c_id] = cpro_pcvto[c_id] * strain_sq[c_id]
                  - d2s3*rho*xk*divu[c_id];
    });
  }

  /* Save production for post processing */
  if (f_tke_prod != nullptr) {
    cs_real_t *tke_prod = f_tke_prod->val;
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      tke_prod[c_id] = prodk[c_id] / crom[c_id];
    });
  }


  /* Take into account rotation/curvature correction, if necessary
     ============================================================ */

  /* Cazalbou correction: the Ceps2 coefficient of destruction term of epsislon
     is modified by rotation and curvature */

  /* Allocate an array for the modified Ceps2 coefficient */
  cs_real_t *ce1rc, *ce2rc;
  CS_MALLOC_HD(ce1rc, n_cells, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(ce2rc, n_cells, cs_real_t, cs_alloc_mode);

  /* Initialization of C_eps1 and C_eps2 k-epsilon models */
  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    ce2rc[c_id] = ce2;
    ce1rc[c_id] = ce1;
  });

  /* Compute the modified Ceps2 coefficient in case of rotation
   * (w1 array not used) */
  if (cs_glob_turb_rans_model->irccor == 1)
    cs_turbulence_rotation_correction(dt, w1, ce2rc);

  /* Launder-Sharma k-epsilon model
    or Baglietto quadratic k-epsilon model */
  if ( model == CS_TURB_K_EPSILON_LS
    || model == CS_TURB_K_EPSILON_QUAD) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cs_real_t rho  = crom[c_id];
      cs_real_t xeps = cvar_ep[c_id];
      cs_real_t xk   = cvar_k[c_id];
      ce2rc[c_id] *= (1. - 0.3*exp(-cs_math_pow2( rho*cs_math_pow2(xk)
                                                /viscl[c_id]
                                                /xeps)));
    });
  }

  if (model == CS_TURB_K_EPSILON_QUAD) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cs_real_t rho  = crom[c_id];
      cs_real_t xeps = cvar_ep[c_id];
      cs_real_t xk   = cvar_k[c_id];

      cs_real_t xdist = cs::max(w_dist[c_id], cs_math_epzero);
      cs_real_t xrey  = rho*xdist*sqrt(xk)/viscl[c_id];
      cs_real_t xpk   = prodk[c_id];
      cs_real_t xpkp
        =   1.33*(1. - 0.3*exp(-cs_math_pow2(rho*cs_math_pow2(xk)
                                             /viscl[c_id]/xeps)))
          * (xpk + 2.*viscl[c_id]*xk/cs_math_pow2(xdist))
          * exp(-3.75e-3*cs_math_pow2(xrey));

      ce1rc[c_id] *= (1. + xpkp/fmax(xpk, 1.e-10));
    });
  }

  /* Calculation of Ceps1* */
  if (model == CS_TURB_V2F_PHI) {
    const cs_real_t cv2fa1 = cs_turb_cv2fa1;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

      cs_real_t dsqphi = sqrt(1./fmax(cvara_phi[c_id], cs_math_epzero));
      /* Ceps1* = 1.4 * (1 - cv2fa1 sqrt(1/phi)) */
      ce1rc[c_id] *= (1. + cv2fa1*dsqphi);
    });
  }

  cs_real_t *coefap = nullptr, *coefbp = nullptr;
  cs_real_t *cofafp = nullptr, *cofbfp = nullptr;

  cs_real_t *w10 = nullptr;

  /* Calculation of Ceps2* */
  if (model == CS_TURB_V2F_BL_V2K) {

    CS_MALLOC_HD(w10, n_cells_ext, cs_real_t, cs_alloc_mode);

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      w3[c_id] = cpro_pcvto[c_id] / cromo[c_id] / sigmak;
    });
    ctx.wait();

    cs_face_viscosity(m,
                      fvq,
                      eqp_k->imvisf,
                      w3,
                      viscf,
                      viscb);

    coefap = f_k->bc_coeffs->a;
    coefbp = f_k->bc_coeffs->b;
    cofafp = f_k->bc_coeffs->af;
    cofbfp = f_k->bc_coeffs->bf;

    /* Translate coefa into cofaf and coefb into cofbf */
    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {

      cs_lnum_t c_id = b_face_cells[face_id];

      cs_real_t hint = w3[c_id]/distb[face_id];

      /* Translate coefa into cofaf and coefb into cofbf */
      cofafp[face_id] = -hint*coefap[face_id];
      cofbfp[face_id] = hint*(1.-coefbp[face_id]);

    });

    cs_diffusion_potential(f_k->id,
                           m,
                           fvq,
                           1,     /* init */
                           1,     /* inc */
                           eqp_k->imrgra,
                           eqp_k->nswrgr,
                           eqp_k->imligr,
                           0,     /* iphydp */
                           eqp_k->iwgrec,
                           eqp_k->verbosity,
                           eqp_k->epsrgr,
                           eqp_k->climgr,
                           nullptr,
                           cvara_k,
                           f_k->bc_coeffs,
                           viscf,
                           viscb,
                           w3,
                           w10);

    const cs_real_t cpale4 = cs_turb_cpale4;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      w10[c_id] = -w10[c_id]/volume[c_id]/cvara_ep[c_id];
      w10[c_id] = tanh(pow(fabs(w10[c_id]), 1.5));
      cs_real_t xcr_m1 = 1.;
      /* HTLES method */
      if (hybrid_turb == 4)
        xcr_m1 = 1. - hybrid_fd_coeff[c_id];
      ce2rc[c_id] *= (1.
                        -  (ce2-cpale4) / ce2
                          * w10[c_id]
                          * cs_math_pow3(cvara_al[c_id])
                          * xcr_m1);
    });

  }

  ctx.wait();
  CS_FREE_HD(w10);

  /* Compute the buoyancy term
   ======================================================== */

  /* Buoyant term for the Atmospheric module
     (function of the potential temperature) */
  if (   cs_glob_turb_rans_model->has_buoyant_term == 1
      && cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 1) {

    cs_atmo_buoyancy_ke_prod(gk);

    /* --- Buoyancy term     G = Beta*g.Grad(scalar)/prdtur/rho
       Here is computed  G =-g.grad(rho)/prdtur/rho */
  }
  else if (cs_glob_turb_rans_model->has_buoyant_term == 1) {

    /* Allocate a temporary for the gradient calculation */
    cs_real_3_t *grad;
    cs_real_t  *grad_dot_g;
    CS_MALLOC_HD(grad, n_cells_ext, cs_real_3_t, cs_alloc_mode);
    CS_MALLOC_HD(grad_dot_g, n_cells_ext, cs_real_t, cs_alloc_mode);

    cs_real_t prdtur = 1;

    cs_field_t *f_thm = cs_thermal_model_field();
    if (phase_id >= 0)
      f_thm = CS_FI_(h_tot, phase_id);

    if (f_thm != nullptr) {
      prdtur = cs_field_get_key_double(f_thm,
                                       cs_field_key_id("turbulent_schmidt"));
    }

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
        grad_dot_g[c_id]
          = - cpro_beta[c_id] * (  cs_math_3_dot_product(grad[c_id], grav)
                                 / prdtur);
      });

    }
    else {

      /* BCs on rho: Dirichlet ROMB
         NB: viscb is used as COEFB */

      ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
        viscb[face_id] = 0.;
      });
      ctx.wait();

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
                         1,     /* w_stride */
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
        grad_dot_g[c_id] =   cs_math_3_dot_product(grad[c_id], grav)
                           / (cromo[c_id]*prdtur);
      });
    }

    /* Production term due to buoyancy */

    cs_field_t *f_tke_buoy = cs_field_by_name_try("algo:k_buoyancy");
    cs_real_t *tke_buoy = nullptr;
    if (f_tke_buoy != nullptr)
      tke_buoy = f_tke_buoy->val;

    /* smbr* store mu_TxS**2 */
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cs_real_t rho   = cromo[c_id];
      cs_real_t visct = cpro_pcvto[c_id];

      /* Explicit Buoyant terms */
      gk[c_id] = -visct*grad_dot_g[c_id];

     /* Save for post processing */
      if (f_tke_buoy != nullptr)
        tke_buoy[c_id] = gk[c_id]/rho;
    });

    ctx.wait();

    /* Free memory */
    CS_FREE_HD(grad);
    CS_FREE_HD(grad_dot_g);
  }
  else {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      gk[c_id] = 0.;
    });
  }

  /* In v2f, we store the production in prdv2f which will be complete further
     for containing the complete production term */
  if (cs_glob_turb_model->itytur == 5) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      prdv2f[c_id] = prodk[c_id] + gk[c_id];
    });
  }

  /* Only for the bl-v2/k model, calculation of E
   * ========================================================
   *      The terms are stored in          e_term
   *      viscf, viscb
   *      Going out of the step we keep e_term */

  if (model == CS_TURB_V2F_BL_V2K) {

    /* Calculation of 2*Ceps3*(1-alpha)^3*nu*nut/eps*d2Ui/dxkdxj*d2Ui/dxkdxj:
       (i.e. E term / k)           : it is stored in e_term */

    cs_real_t *w12 = nullptr;
    /* Allocate a work array */
    CS_MALLOC_HD(w12, n_cells_ext, cs_real_t, cs_alloc_mode);

    _tsepls(phase_id, w12);

    const cs_real_t cpale3 = cs_turb_cpale3;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cs_real_t rho  = cromo[c_id];
      cs_real_t xnu  = cpro_pcvlo[c_id]/rho;
      cs_real_t xnut = cpro_pcvto[c_id]/rho;
      cs_real_t xeps = cvara_ep[c_id];

      e_term[c_id] = 2.*xnu*xnut*w12[c_id]*cpale3/xeps
                    *cs_math_pow3(1.-cvara_al[c_id]);

    });
    ctx.wait();

    /* Free memory */
    CS_FREE_HD(w12);

  }

  /* Only for the Launder-Sharma model, calculation of E and D terms
   *      The terms are stored in          grad_sqk, grad_s
   * ================================================================*/

  if (model == CS_TURB_K_EPSILON_LS) {

    /* Gradient of square root of k
     * -----------------------------*/

    coefap = (cs_real_t *)f_k->bc_coeffs->a;
    coefbp = (cs_real_t *)f_k->bc_coeffs->b;

    /* For all usual type of boundary faces (wall, inlet, sym, outlet):
       - coefa for sqrt(k) is the sqrt of the coefa for k,
       - coefb is the same as for k */

    cs_field_bc_coeffs_t bc_coeffs_sqk_loc;
    cs_field_bc_coeffs_shallow_copy(f_k->bc_coeffs, &bc_coeffs_sqk_loc);
    CS_MALLOC_HD(bc_coeffs_sqk_loc.a, n_b_faces, cs_real_t, cs_alloc_mode);

    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      bc_coeffs_sqk_loc.a[face_id] = sqrt(coefap[face_id]);
    });
    ctx.wait();

    cs_halo_type_t halo_type = CS_HALO_STANDARD;
    cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

    cs_gradient_type_by_imrgra(eqp_k->imrgra,
                               &gradient_type,
                               &halo_type);

    cs_gradient_scalar("grad_sqk",
                       gradient_type,
                       halo_type,
                       1,     /* inc */
                       eqp_k->nswrgr,
                       0,
                       1,     /* w_stride */
                       eqp_k->verbosity,
                       static_cast<cs_gradient_limit_t>(eqp_k->imligr),
                       eqp_k->epsrgr,
                       eqp_k->climgr,
                       nullptr,
                       &bc_coeffs_sqk_loc,
                       sqrt_k,
                       nullptr,
                       nullptr,  /* internal coupling */
                       grad_sqk);

    cs_field_bc_coeffs_free_copy(f_k->bc_coeffs, &bc_coeffs_sqk_loc);

    /* Gradient of the Strain (grad S)
       --------------------------------------- */

    cs_field_bc_coeffs_t bc_coeffs_sqs_loc;
    cs_field_bc_coeffs_init(&bc_coeffs_sqs_loc);

    CS_MALLOC_HD(bc_coeffs_sqs_loc.a, n_b_faces, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(bc_coeffs_sqs_loc.b, n_b_faces, cs_real_t, cs_alloc_mode);

    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      bc_coeffs_sqs_loc.a[face_id] = 0.;
      bc_coeffs_sqs_loc.b[face_id] = 1.;
    });
    ctx.wait();

    cs_gradient_type_by_imrgra(eqp_k->imrgra,
                               &gradient_type,
                               &halo_type);

    cs_gradient_scalar("grad_s",
                       gradient_type,
                       halo_type,
                       1,     /* inc */
                       eqp_k->nswrgr,
                       0,
                       1,     /* w_stride */
                       eqp_k->verbosity,
                       static_cast<cs_gradient_limit_t>(eqp_k->imligr),
                       eqp_k->epsrgr,
                       eqp_k->climgr,
                       nullptr,
                       &bc_coeffs_sqs_loc,
                       strain,
                       nullptr,
                       nullptr, /* internal coupling */
                       grad_s);

    CS_FREE_HD(bc_coeffs_sqs_loc.a);
    CS_FREE_HD(bc_coeffs_sqs_loc.b);

  }

  /* Finalization of explicit and implicit source terms
   * ===================================================
   *
   * smbre = epsilon/k (ce1 prod + ce3 G - rho ce2 epsilon)
   * smbrk =                prod +     G - rho epsilon */

  /* If we extrapolate the source terms and rho, we need here rho^n
     and visct, we need here visct^n */

  if (cs_glob_turb_model->itytur == 2) {

    /* Stores the production terms for the k-epsilon coupling option */
    if (model == CS_TURB_K_EPSILON) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        prdtke[c_id] = prodk[c_id] + gk[c_id];
        // Note: only to be consistent with previouls formulation
        // but not consistant with variant of ce3
        prdeps[c_id] = prodk[c_id] + cs::max(gk[c_id], 0.);
      });
    }
    ctx.wait();

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

      cs_real_t rho = cromo[c_id];
      cs_real_t d_tt = cvara_ep[c_id]/cvara_k[c_id];

      smbrk[c_id] =   cell_f_vol[c_id]
                    * (prodk[c_id] + gk[c_id] - rho*cvara_ep[c_id]);

      cs_real_t ce3_loc = ce3;
      if (dissip_buo_mdl == 0 && gk[c_id] < 0.)
        ce3_loc = 0.;

      smbre[c_id] =   cell_f_vol[c_id]
                    * d_tt * (  ce1rc[c_id]*prodk[c_id] + ce3_loc*gk[c_id]
                              - ce2rc[c_id]*rho*cvara_ep[c_id]);

    });

    if (model == CS_TURB_K_EPSILON_LS) {
      const int ikecou = cs_glob_turb_rans_model->ikecou;

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        cs_real_t rho   = cromo[c_id];
        cs_real_t mu   = cpro_pcvlo[c_id];
        cs_real_t xnut  = cpro_pcvto[c_id]/rho;

        cs_real_t mu_gradk2 = cell_f_vol[c_id] *  2. * mu
                        * cs_math_3_square_norm(grad_sqk[c_id]);

        if (ikecou == 0)
          tinstk[c_id] += mu_gradk2/cvara_k[c_id];

        smbrk[c_id] -= mu_gradk2;

        smbre[c_id] +=    cell_f_vol[c_id] * 2. * mu * xnut
                        * cs_math_3_square_norm(grad_s[c_id]);
      });
    }

    /* If the solving of k-epsilon is uncoupled,
     * negative source terms are implicited */

    if (cs_glob_turb_rans_model->ikecou == 0) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        cs_real_t xeps = cvara_ep[c_id];
        cs_real_t xk   = cvara_k[c_id];
        cs_real_t rho = crom[c_id];
        cs_real_t d_tt = xeps / xk;
        tinstk[c_id] += cell_f_vol[c_id] * (rho*d_tt
                                            + cs::max(d2s3*rho*divu[c_id], 0.)
                                            + cs::max(- gk[c_id] / xk, 0.));

        cs_real_t ce3_loc = ce3;
        if (dissip_buo_mdl == 0 && gk[c_id] < 0.)
          ce3_loc = 0.;

        tinste[c_id] += cell_f_vol[c_id] * (
                             ce2rc[c_id]*rho * d_tt
                             /* Note: d_tt = eps/k
                              *      -(eps/k) * k div(u)
                              * gives
                              *       div(u)
                              * as implicit contribution for eps
                              */
                           + cs::max(d2s3*ce1rc[c_id]*rho*divu[c_id], 0.)
                           + ce3_loc * d_tt * cs::max(-gk[c_id]/xeps, 0.));
      });
    }

  }
  else if (model == CS_TURB_V2F_PHI) {

    const cs_real_t cv2fct = cs_turb_cv2fct;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

      cs_real_t rho  = cromo[c_id];
      cs_real_t xnu  = cpro_pcvlo[c_id]/rho;
      cs_real_t xeps = cvara_ep[c_id];
      cs_real_t xk   = cvara_k[c_id];
      cs_real_t ttke = xk / xeps;
      cs_real_t ttmin = cv2fct*sqrt(xnu/xeps);
      cs_real_t d_tt = 1. / fmax(ttke, ttmin);

      /* Explicit part */
      smbrk[c_id] = cell_f_vol[c_id] *
                    ( prodk[c_id] + gk[c_id] - rho*cvara_ep[c_id]);

      // TODO propose other models for G_epsilon (dissip_buo_mdl)
      smbre[c_id] = cell_f_vol[c_id] * d_tt
                  *( ce1rc[c_id]*(prodk[c_id] + cs::max(gk[c_id], 0.))
                   - ce2rc[c_id]*rho*xeps);

      /* Implicit part */
      tinstk[c_id] += cell_f_vol[c_id]
        * ( rho / fmax(ttke, cs_math_epzero * ttmin)
          + d2s3 * rho * cs::max(divu[c_id], 0.)
          + cs::max(-gk[c_id] / xk, 0.));

      tinste[c_id] += cell_f_vol[c_id] * d_tt * (rho * ce2rc[c_id]
                      + rho*cs::max(d2s3*ce1rc[c_id]*ttke*divu[c_id], 0.));
    });

  }
  else if (model == CS_TURB_V2F_BL_V2K) {

    const cs_real_t cpalct = cs_turb_cpalct;

    if (cs_glob_turb_model->hybrid_turb == 0) {

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

        cs_real_t rho   = cromo[c_id];
        cs_real_t xnu   = cpro_pcvlo[c_id]/rho;
        cs_real_t xeps  = cvara_ep[c_id];
        cs_real_t xk    = cvara_k[c_id];
        cs_real_t ttke  = xk / xeps;
        cs_real_t ttmin = cpalct*sqrt(xnu/xeps);
        cs_real_t d_tt  = 1. / sqrt(cs_math_pow2(ttke) + cs_math_pow2(ttmin));

        /* Explicit part */
        smbrk[c_id] = cell_f_vol[c_id] * ( prodk[c_id] + gk[c_id]
                                          - rho*xeps
                                          - rho*e_term[c_id]*xk);

        // TODO propose other models for G_epsilon (dissip_buo_mdl)
        smbre[c_id] =   cell_f_vol[c_id]
                      * d_tt*( ce1rc[c_id]
                              * (prodk[c_id] + cs::max(gk[c_id], 0.))
                             - ce2rc[c_id]*rho*xeps);

        /* Implicit part */
        tinstk[c_id] += cell_f_vol[c_id] * rho
          * ( 1./fmax(ttke, cs_math_epzero * ttmin)
            + cs::max(d2s3 * divu[c_id], 0.)
            + e_term[c_id]); /* Note that e_term is positive */
        /* Note that ce2rc is positive */
        tinste[c_id] +=  rho * cell_f_vol[c_id] * d_tt
          * (   ce2rc[c_id]
              + cs::max(d2s3*ce1*ttke*divu[c_id], 0.));
      });
    }
    else if (cs_glob_turb_model->hybrid_turb == 4) {
      /* HTLES */

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

        cs_real_t rho   = cromo[c_id];
        cs_real_t xnu   = cpro_pcvlo[c_id]/rho;
        cs_real_t xeps  = cvara_ep[c_id];
        cs_real_t xk    = cvara_k[c_id];

        /* HTLES method */
        cs_real_t xpsi  = htles_psi[c_id];
        cs_real_t xtm   = htles_t[c_id];
        cs_real_t xepsh = xk/xtm;

        /* Modif. definition of Kolmogorov time scale */
        cs_real_t ttke  = xk / xeps;
        cs_real_t ttmin = cpalct*sqrt(xnu/(xpsi*xeps));
        cs_real_t d_tt    = 1. / sqrt(cs_math_pow2(ttke) + cs_math_pow2(ttmin));

        /* Explicit part */
        smbrk[c_id] = cell_f_vol[c_id] * ( prodk[c_id] + gk[c_id]
                                          - rho*xepsh
                                          - rho*e_term[c_id]*xk);

        // TODO propose other models for G_epsilon (dissip_buo_mdl)
        smbre[c_id] =   cell_f_vol[c_id] * d_tt
                      * ( ce1rc[c_id]*(prodk[c_id]+cs::max(gk[c_id], 0.))
                        - ce2rc[c_id]*rho*xeps);

        /* Implicit part */
        tinstk[c_id] +=   rho*cell_f_vol[c_id]
                        / cs::max(xtm, cs_math_epzero * ttmin);

        tinstk[c_id] += cs::max(d2s3*rho*cell_f_vol[c_id]*divu[c_id], 0.);
        /* Note that e_term is positive */
        tinstk[c_id] += e_term[c_id]*rho*cell_f_vol[c_id];
        /* Note that ce2rc is positive */
        tinste[c_id] += rho*cell_f_vol[c_id] * d_tt * ( ce2rc[c_id]
                          + cs::max(d2s3 * ce1rc[c_id] * ttke * divu[c_id],
                                    0.));
      });
    }

  }

  ctx.wait();

  /* Take user source terms into account
   *
   *    The square of the scalar strain rate (strain_sq)
   *    and the trace of the velocity gradient (divu) are available.
   *    The part to be explicit is stored in       usexpk, usexpe
   *    The part to be implicit is stored in       usimpk, usimpe
   *    Going out of the step we keep              strain_sq, divu,
   * ==========================================================================*/

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    usimpk[c_id] = 0.;
    usimpe[c_id] = 0.;
    usexpk[c_id] = 0.;
    usexpe[c_id] = 0.;
  });
  ctx.wait();

  cs_user_source_terms(domain,
                       f_k->id,
                       usexpk,
                       usimpk);

  cs_user_source_terms(domain,
                       f_eps->id,
                       usexpe,
                       usimpe);

  if (cs_glob_porous_model == 3) {
    cs_immersed_boundary_wall_functions(f_k->id, usexpk, usimpk);
    cs_immersed_boundary_wall_functions(f_eps->id, usexpe, usimpe);
  }

  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 0) {

    int k_interp_id = cs_field_key_id("opt_interp_id");

    /* Nudging towards optimal interpolation for k */
    if (cs_field_get_key_int(f_k, k_interp_id) >= 0)
      cs_at_data_assim_source_term(f_k->id, usexpk, usimpk);

    /* Nudging towards optimal interpolation for epsilon */
    if (cs_field_get_key_int(f_eps, k_interp_id) >= 0)
      cs_at_data_assim_source_term(f_k->id, usexpe, usimpe);

  }

  /* If source terms are extrapolated over time */
  if (istprv >= 0) {

    cs_real_t thetak = eqp_k->theta;
    cs_real_t thetae = eqp_eps->theta;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

      /* Recover the value at time (n-1) */
      cs_real_t tuexpk = c_st_k_p[c_id];
      cs_real_t tuexpe = c_st_eps_p[c_id];

      /* Save the values for the next time-step */
      c_st_k_p[c_id]   = smbrk[c_id] + usexpk[c_id];
      c_st_eps_p[c_id] = smbre[c_id] + usexpe[c_id];

      /* Explicit Part */
      smbrk[c_id] = - thets*tuexpk;
      smbre[c_id] = - thets*tuexpe;
      /* It is assumed that (-usimpk > 0) and though this term is implicit */
      smbrk[c_id] += usimpk[c_id]*cvara_k[c_id];
      smbre[c_id] += usimpe[c_id]*cvara_ep[c_id];

      /* Implicit part */
      tinstk[c_id] -= usimpk[c_id]*thetak;
      tinste[c_id] -= usimpe[c_id]*thetae;

    });

    /* If no extrapolation over time */
  }
  else {

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      /* Explicit part */
      smbrk[c_id] += usimpk[c_id]*cvara_k[c_id] + usexpk[c_id];
      smbre[c_id] += usimpe[c_id]*cvara_ep[c_id] + usexpe[c_id];

      /* Implicit part */
      tinstk[c_id] += cs::max(-usimpk[c_id], 0.);
      tinste[c_id] += cs::max(-usimpe[c_id], 0.);
    });

  }

  /* Account for Lagrangian 2-way coupling source terms
     -------------------------------------------------- */

  /* 2nd order not handled */

  if (cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING) {

    const cs_lagr_source_terms_t  *lag_st = cs_glob_lagr_source_terms;

    if (lag_st->ltsdyn == 1) {

      cs_real_t *lag_st_k = cs_field_by_name("lagr_st_k")->val;
      cs_real_t *lag_st_i = cs_field_by_name("lagr_st_imp_velocity")->val;

      const cs_real_t ce4 = cs_turb_ce4;

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

        /* Implicit and explicit source terms on k */
        cs_real_t st_k = cell_f_vol[c_id] * lag_st_k[c_id];
        cs_real_t st_i = cell_f_vol[c_id] * lag_st_i[c_id];
        smbrk[c_id] += st_k;

        /* Explicit source terms on epsilon */

        smbre[c_id] += ce4 * st_k * cvara_ep[c_id]
                                  / cvara_k[c_id];

        /* Implicit source terms on k,
         * Note: we implicit lag_st_k if negative */
        tinstk[c_id] += cs::max(-st_k/cvara_k[c_id], 0.);
        tinstk[c_id] += cs::max(-st_i, 0.);

        /* Implicit source terms on epsilon */
        tinste[c_id] += cs::max(-ce4 * st_k / cvara_k[c_id], 0.);

      });

    }

  }

  ctx.wait();

  /* Mass source terms (Implicit and explicit parts)
   * Going out of the step we keep divu,  smbrk, smbre
   * ================================================= */

  if (   eqp_k->n_volume_mass_injections > 0
      || eqp_eps->n_volume_mass_injections > 0) {

    int *mst_type = nullptr;
    cs_lnum_t n_elts = 0;
    const cs_lnum_t *elt_ids = nullptr;
    cs_real_t *mst_val = nullptr, *mst_val_p = nullptr;
    cs_real_t *gapinj_k = nullptr, *gapinj_eps = nullptr;

    /* If we extrapolate the source terms we put Gamma Pinj in c_st */
    if (istprv >= 0) {
      gapinj_k = c_st_k_p;
      gapinj_eps = c_st_eps_p;
    }
    /* Otherwise we put it directly in smbr */
    else {
      gapinj_k = smbrk;
      gapinj_eps = smbre;
    }

    /* We increment smbrs with -Gamma.var_prev and rovsdt with Gamma */

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

    /* For eps */
    cs_volume_mass_injection_get_arrays(f_eps,
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
                         cvara_ep,
                         mst_val,
                         mst_val_p,
                         smbre,
                         tinste,
                         gapinj_eps);

  }

  /* Taking into account the terms of conv/diff in the second member for the
   * strengthened coupling k-epsilon (kecou == 1)
   *      Work table                       w4
   *      The terms are stored in usexpk and usexpe, then added to smbrk, smbre
   *      Going out of the step we keep    divu,
   *      smbrk, smbre
   *      usexpk, usexpe
   * ==========================================================================*/

  if (cs_glob_turb_rans_model->ikecou == 1) {

    ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      usexpk[c_id] = 0.;
      usexpe[c_id] = 0.;
    });
    ctx.wait();

    /* Handle k */

    coefap = f_k->bc_coeffs->a;
    coefbp = f_k->bc_coeffs->b;
    cofafp = f_k->bc_coeffs->af;
    cofbfp = f_k->bc_coeffs->bf;

    if (eqp_k->idiff >=  1) {

      if (eqp_k->idifft >= 1) {
        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          w4[c_id] = viscl[c_id] + cvisct[c_id]/sigmak;
        });
        ctx.wait();

        cs_face_viscosity(m,
                          fvq,
                          eqp_k->imvisf,
                          w4,
                          viscf,
                          viscb);
      }
      else {
        cs_face_viscosity(m,
                          fvq,
                          eqp_k->imvisf,
                          const_cast<cs_real_t *>(viscl),
                          viscf,
                          viscb);
      }

    }
    else {
      ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
        viscf[face_id] = 0.;
      });
      ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
        viscb[face_id] = 0.;
      });
      ctx.wait();
    }
    cs_equation_param_t eqp_k_loc = *eqp_k;
    eqp_k_loc.idften = CS_ISOTROPIC_DIFFUSION;

    cs_balance_scalar(cs_glob_time_step_options->idtvar,
                      f_k->id,
                      0,     /* imucpp */
                      1,     /* imasac */
                      1,     /* inc */
                      &eqp_k_loc,
                      cvara_k,
                      cvara_k,
                      f_k->bc_coeffs,
                      imasfl,
                      bmasfl,
                      viscf,
                      viscb,
                      nullptr,
                      nullptr,
                      nullptr,
                      nullptr,
                      0, /* boundary convective flux with upwind */
                      nullptr,
                      usexpk);

    if (eqp_k->verbosity >= 2) {
      cs_log_printf(CS_LOG_DEFAULT,
                    " Variable %s: explicit balance = %12.5e\n",
                    cs_field_get_label(f_k),
                    sqrt(cs_gdot(n_cells, smbrk, smbrk)));
    }

    /* Handle epsilon */

    coefap = f_eps->bc_coeffs->a;
    coefbp = f_eps->bc_coeffs->b;
    cofafp = f_eps->bc_coeffs->af;
    cofbfp = f_eps->bc_coeffs->bf;

    if (eqp_eps->idiff >=  1) {

      if (eqp_eps->idifft == 1) {
        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          w4[c_id] = viscl[c_id] + cvisct[c_id]/sigmae;
        });
        ctx.wait();

        cs_face_viscosity(m, fvq, eqp_eps->imvisf,
                          w4,
                          viscf, viscb);
      }
      else {
        cs_face_viscosity(m, fvq, eqp_eps->imvisf,
                          const_cast<cs_real_t *>(viscl),
                          viscf, viscb);
      }

    }
    else {
      ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
        viscf[face_id] = 0.;
      });
      ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
        viscb[face_id] = 0.;
      });
      ctx.wait();
    }

    cs_equation_param_t eqp_eps_loc = *eqp_eps;
    eqp_eps_loc.idften = CS_ISOTROPIC_DIFFUSION;

    cs_balance_scalar(cs_glob_time_step_options->idtvar,
                      f_eps->id,
                      0,    /* imucpp */
                      1,    /* imasac */
                      1,    /* inc */
                      &eqp_eps_loc,
                      cvara_ep,
                      cvara_ep,
                      f_eps->bc_coeffs,
                      imasfl,
                      bmasfl,
                      viscf,
                      viscb,
                      nullptr,
                      nullptr,
                      nullptr,
                      nullptr,
                      0, /* boundary convective flux with upwind */
                      nullptr,
                      usexpe);

    if (eqp_eps->verbosity >= 2) {
      cs_log_printf(CS_LOG_DEFAULT,
                    " Variable %s: EXPLICIT BALANCE =  %12.5e\n",
                    cs_field_get_label(f_eps),
                    sqrt(cs_gdot(n_cells, smbre, smbre)));
    }

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      smbrk[c_id] += usexpk[c_id];
      smbre[c_id] += usexpe[c_id];
    });

  }

  /* k-Epsilon coupling (kecou == 1)
   * =============================== */

  /* Second order is not taken into account */
  if (cs_glob_turb_rans_model->ikecou == 1) {

    if (model == CS_TURB_K_EPSILON) {

      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

        cs_real_t rho   = crom[c_id];
        cs_real_t visct = cpro_pcvto[c_id];

        /* Coupled solving */
        cs_real_t romvsd      = 1./(rho*volume[c_id]);
        cs_real_t divp23      = d2s3*fmax(divu[c_id], 0.);
        cs_real_t d_tt        = cvara_ep[c_id]/cvara_k[c_id];

        smbrk[c_id] *= romvsd;
        smbre[c_id] *= romvsd;

        cs_real_t a11 = 1./dt[c_id] - 2.*cvara_k[c_id]/cvara_ep[c_id]
                                      *cmu*cs::min(prdtke[c_id]/visct, 0.)
                                    + divp23;
        cs_real_t a12 = 1.;
        cs_real_t a21 = - ce1rc[c_id]*cmu*prdeps[c_id]/visct
                        - ce2rc[c_id] * d_tt * d_tt;
        cs_real_t a22 = 1./dt[c_id] + ce1rc[c_id]*divp23
                                    + 2.*ce2rc[c_id] * d_tt;

        cs_real_t unsdet = 1./(a11*a22 -a12*a21);

        cs_real_t deltk = ( a22*smbrk[c_id] - a12*smbre[c_id])*unsdet;
        cs_real_t delte = (-a21*smbrk[c_id] + a11*smbre[c_id])*unsdet;

        /* New source term for the iterative process */
        romvsd = rho*cell_f_vol[c_id]/dt[c_id];

        smbrk[c_id] = romvsd*deltk;
        smbre[c_id] = romvsd*delte;

        /* we remove the convection/diffusion at time n from smbrk and smbre
           if they were calculated */
        smbrk[c_id] -= usexpk[c_id];
        smbre[c_id] -= usexpe[c_id];

      });

      /* In cs_parameters_check.c we block the
         model != CS_TURB_K_EPSILON / kecou = 1
         combination */
    }
    else
      bft_error(__FILE__, __LINE__, 0,
              _("ikecou=1 is not validated with this turbulent model\n"
                "---------------------------------------------------"));

  }

  /* Finalization of the Right Hand Side when activating 2nd time order
   * ================================================================== */

  if (istprv >= 0) {
    cs_real_t thetp1 = 1. + thets;
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      smbrk[c_id] += thetp1 * c_st_k_p[c_id];
      smbre[c_id] += thetp1 * c_st_eps_p[c_id];
    });
  }

  ctx.wait();

  cs_solid_zone_set_zero_on_cells(1, smbrk);

  /* Solve for turbulent kinetic energy (k)
     ---------------------------------- */

   /* We use       smbrk, tinstk
    * Work array   w1 */

  coefap = f_k->bc_coeffs->a;
  coefbp = f_k->bc_coeffs->b;
  cofafp = f_k->bc_coeffs->af;
  cofbfp = f_k->bc_coeffs->bf;

  /* Face viscosity */
  if (eqp_k->idiff >= 1) {

    const cs_real_t viscl_mult
      = (model == CS_TURB_V2F_BL_V2K) ? 0.5 : 1.0;

    if (eqp_k->idifft == 1) {
      /* Variable Schmidt number */
      if (sigmak_id >= 0) {
        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          w1[c_id] = viscl[c_id]*viscl_mult + cvisct[c_id]/cpro_sigmak[c_id];
        });
      }
      else {
        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          w1[c_id] = viscl[c_id]*viscl_mult + cvisct[c_id]/sigmak;
        });
      }
    }
    else {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        w1[c_id] = viscl[c_id]*viscl_mult;
      });
    }

    ctx.wait();

    cs_face_viscosity(m,
                      fvq,
                      eqp_k->imvisf,
                      w1,
                      viscf,
                      viscb);

  }
  else {
    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      viscf[face_id] = 0.;
    });
    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      viscb[face_id] = 0.;
    });
    ctx.wait();
  }

  /* Solve k */

  cs_equation_param_t eqp_k_loc = *eqp_k;

  cs_equation_iterative_solve_scalar(cs_glob_time_step_options->idtvar,
                                     1,  /* init */
                                     f_k->id,
                                     f_k->name,
                                     0,  /* iescap */
                                     0,  /* imucpp */
                                     -1, /* normp */
                                     &eqp_k_loc,
                                     cvara_k,
                                     cvara_k,
                                     f_k->bc_coeffs,
                                     imasfl,
                                     bmasfl,
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

  /* Solve for epsilon */

  coefap = f_eps->bc_coeffs->a;
  coefbp = f_eps->bc_coeffs->b;
  cofafp = f_eps->bc_coeffs->af;
  cofbfp = f_eps->bc_coeffs->bf;

  /* Face viscosity */
  if (eqp_eps->idiff >= 1) {

    const cs_real_t viscl_mult
      = (model == CS_TURB_V2F_BL_V2K) ? 0.5 : 1.0;

    if (eqp_eps->idifft == 1) {
      /* Variable Schmidt number */
      if (sigmae_id >= 0) {
        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          w1[c_id] = viscl[c_id]*viscl_mult + cvisct[c_id]/cpro_sigmae[c_id];
        });
      }
      else {
        ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
          w1[c_id] = viscl[c_id]*viscl_mult + cvisct[c_id]/sigmae;
        });
      }
    }
    else {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        w1[c_id] = viscl[c_id]*viscl_mult;
      });
    }

    ctx.wait();

    cs_face_viscosity(m,
                      fvq,
                      eqp_eps->imvisf,
                      w1,
                      viscf,
                      viscb);
  }
  else {
    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      viscf[face_id] = 0.;
    });
    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      viscb[face_id] = 0.;
    });
    ctx.wait();
  }

  /* Solve epsilon */

  cs_equation_param_t eqp_eps_loc = *eqp_eps;

  cs_equation_iterative_solve_scalar(cs_glob_time_step_options->idtvar,
                                     1,    /* init */
                                     f_eps->id,
                                     f_eps->name,
                                     0,   /* iescap */
                                     0,   /* imucpp */
                                     -1,  /* normp */
                                     &eqp_eps_loc,
                                     cvara_ep,
                                     cvara_ep,
                                     f_eps->bc_coeffs,
                                     imasfl,
                                     bmasfl,
                                     viscf,
                                     viscb,
                                     viscf,
                                     viscb,
                                     nullptr,
                                     nullptr,
                                     nullptr,
                                     0, /* boundary convective upwind flux */
                                     nullptr,
                                     tinste,
                                     smbre,
                                     cvar_ep,
                                     dpvar,
                                     nullptr,
                                     nullptr);

  /* Clip values
     ============ */

  cs_turbulence_ke_clip(phase_id,
                        n_cells,
                        1);

  /* Free memory */
  CS_FREE_HD(_gradv);
  CS_FREE_HD(viscf);
  CS_FREE_HD(viscb);
  CS_FREE_HD(usimpk);
  CS_FREE_HD(smbrk);
  CS_FREE_HD(smbre);
  CS_FREE_HD(prodk);
  CS_FREE_HD(gk);
  CS_FREE_HD(divu);
  CS_FREE_HD(strain_sq);
  CS_FREE_HD(w1);
  CS_FREE_HD(w3);
  CS_FREE_HD(w4);
  CS_FREE_HD(usexpk);
  CS_FREE_HD(usexpe);
  CS_FREE_HD(usimpe);
  CS_FREE_HD(dpvar);
  CS_FREE_HD(ce2rc);
  CS_FREE_HD(ce1rc);

  CS_FREE(tinstk);
  CS_FREE(tinste);

#if defined(HAVE_ACCEL)
  CS_FREE_HD(_grav);
#endif

  CS_FREE_HD(e_term);
  if (model == CS_TURB_K_EPSILON){
    CS_FREE_HD(prdtke);
    CS_FREE_HD(prdeps);
  }
  if (model == CS_TURB_K_EPSILON_LS) {
    CS_FREE_HD(sqrt_k);
    CS_FREE_HD(strain);
    CS_FREE_HD(grad_sqk);
    CS_FREE_HD(grad_s);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Clipping of the turbulent kinetic energy and turbulent dissipation.
 *
 * \param[in]     phase_id   turbulent phase id (-1 for single phase flow)
 * \param[in]     n_cells    number of cells
 * \param[in]     iclip      indicator = 0 if viscl0 is used
 *                           otherwise viscl is used.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_ke_clip(int        phase_id,
                      cs_lnum_t  n_cells,
                      int        iclip)
{
  cs_turb_ref_values_t *turb_ref_values = cs_get_glob_turb_ref_values();
  cs_real_t almax = turb_ref_values->almax;

  const cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();

  cs_real_t viscl0 = phys_pro->viscl0; /* reference pressure */
  cs_real_t ro0    = phys_pro->ro0; /* reference density */

  cs_field_t *f_k = CS_F_(k);
  cs_field_t *f_eps = CS_F_(eps);
  cs_field_t *f_rho = CS_F_(rho);
  cs_field_t *f_mu = CS_F_(mu);

  if (phase_id >= 0) {
    f_k = CS_FI_(k, phase_id);
    f_eps = CS_FI_(eps, phase_id);
    f_rho = CS_FI_(rho, phase_id);
    f_mu = CS_FI_(mu, phase_id);
  }

  cs_real_t *crom    = (cs_real_t *)f_rho->val;
  cs_real_t *cvar_k  = (cs_real_t *)f_k->val;
  cs_real_t *cvar_ep = (cs_real_t *)f_eps->val;
  cs_real_t *viscl   =  (cs_real_t *)f_mu->val;

  cs_dispatch_context ctx;

  const cs_equation_param_t *eqp
    = cs_field_get_equation_param_const(f_k);

  /* Small value to avoid exactly zero values */

  const double epz2 = cs_math_pow2(cs_math_epzero);

  /* Postprocess clippings? */

  int key_clipping_id = cs_field_key_id("clipping_id");

  int clip_k_id = cs_field_get_key_int(f_k, key_clipping_id);
  cs_real_t *cpro_k_clipped = nullptr;
  if (clip_k_id >= 0) {
    cpro_k_clipped = cs_field_by_id(clip_k_id)->val;
  }

  int clip_e_id = cs_field_get_key_int(f_eps, key_clipping_id);
  cs_real_t *cpro_e_clipped = nullptr;
  if (clip_e_id >= 0) {
    cpro_e_clipped = cs_field_by_id(clip_e_id)->val;
  }

  /* Save min and max for log
   * ======================== */

  struct cs_data_4r rd;
  struct cs_reduce_min2r_max2r reducer;

  ctx.parallel_for_reduce
    (n_cells, rd, reducer,
     [=] CS_F_HOST_DEVICE (cs_lnum_t c_id, cs_data_4r &res) {
    res.r[0] = cvar_k[c_id];
    res.r[1] = cvar_k[c_id];
    res.r[2] = cvar_ep[c_id];
    res.r[3] = cvar_ep[c_id];
  });

  ctx.wait();

  if (cpro_k_clipped != nullptr) {
    cs_array_real_fill_zero(n_cells, cpro_k_clipped);
  }
  if (cpro_e_clipped != nullptr) {
    cs_array_real_fill_zero(n_cells, cpro_e_clipped);
  }

  /* Detect values ouside "physical" bounds,
   * only for warning or when ICLKEP = 1
   * ===================================== */

  cs_gnum_t iclpke = 0;

  const cs_lnum_t ke_clip_method = cs_glob_turb_rans_model->iclkep;

  if (eqp->verbosity >= 2 || cs_glob_turb_rans_model->iclkep == 1) {

    if (iclip == 1) {

      cs_real_t xkm = 1296.*sqrt(cs_turb_cmu)/cs_math_pow2(almax);
      cs_real_t xepm = 46656.*cs_turb_cmu/cs_math_pow4(almax);

      ctx.parallel_for_reduce_sum
        (n_cells, iclpke, [=] CS_F_HOST_DEVICE
         (cs_lnum_t c_id,
          CS_DISPATCH_REDUCER_TYPE(cs_gnum_t) &sum) {

        cs_real_t xk = cvar_k[c_id];
        cs_real_t xe = cvar_ep[c_id];
        cs_real_t xkmin = xkm * cs_math_pow2(viscl[c_id] / crom[c_id]);
        cs_real_t xepmin = xepm * cs_math_pow3(viscl[c_id] / crom[c_id]);
        if (xk <= xkmin || xe <= xepmin) {
          if (ke_clip_method == 1) {
            if (clip_k_id >= 0)
              cpro_k_clipped[c_id] = xkmin - xk;
            cvar_k[c_id]  = xkmin;
            if (clip_e_id >= 0)
              cpro_e_clipped[c_id] = xepmin - xe;
            cvar_ep[c_id] = xepmin;
          }
          sum += (cs_gnum_t)1;
        }
      });

    }
    else if (iclip == 0) {

      cs_real_t xkm = 1296. * sqrt(cs_turb_cmu) / cs_math_pow2(almax)
        * cs_math_pow2(viscl0/ro0);
      cs_real_t xepm = 46656. * cs_turb_cmu/cs_math_pow4(almax)
        * cs_math_pow3(viscl0/ro0);

      ctx.parallel_for_reduce_sum
        (n_cells, iclpke, [=] CS_F_HOST_DEVICE
         (cs_lnum_t c_id,
          CS_DISPATCH_REDUCER_TYPE(cs_gnum_t) &sum) {

        cs_real_t xk = cvar_k[c_id];
        cs_real_t xe = cvar_ep[c_id];
        if (xk <= xkm || xe <= xepm) {
          if (ke_clip_method == 1) {
            cvar_k[c_id]  = xkm;
            if (clip_k_id >= 0)
              cpro_k_clipped[c_id] = xkm - xk;
            cvar_ep[c_id] = xepm;
            if (clip_e_id >= 0)
              cpro_e_clipped[c_id] = xepm - xe;
          }
          sum += (cs_gnum_t)1;
        }
      });

    }
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Call of %s with option = %d"),
                __func__, iclip);

    ctx.wait();

    /* logging */

    if (eqp->verbosity >= 2) {

      cs_parall_sum(1, CS_GNUM_TYPE, &iclpke);

      cs_log_printf(CS_LOG_DEFAULT,
                    "\n "
                    "%llu k-epsilon values beyond the scales based on almax\n",
                    (unsigned long long)iclpke);
    }

  }

  /* "standard" clipping ICLKEP = 0
   * ============================== */

  struct cs_data_2i rd_sum;
  rd_sum.i[0] = 0;
  rd_sum.i[1] = 0;
  struct cs_reduce_sum2i reducer_sum;

  if (cs_glob_turb_rans_model->iclkep == 0) {

    ctx.parallel_for_reduce
      (n_cells, rd_sum, reducer_sum, [=] CS_F_HOST_DEVICE
       (cs_lnum_t c_id, cs_data_2i &res) {

      res.i[0] = 0, res.i[1] = 0;

      cs_real_t xk = cvar_k[c_id];
      cs_real_t xe = cvar_ep[c_id];
      if (fabs(xk) <= epz2) {
        res.i[0] = 1;
        if (clip_k_id >= 0)
          cpro_k_clipped[c_id] = epz2 - cvar_k[c_id];
        cvar_k[c_id] = cs::max(cvar_k[c_id],epz2);
      }
      else if(xk <= 0.) {
        res.i[0] = 1;
        if (clip_k_id >= 0)
          cpro_k_clipped[c_id] = -xk;
        cvar_k[c_id] = -xk;
      }
      if (fabs(xe) <= epz2) {
        res.i[1] = 1;
        if (clip_e_id >= 0)
          cpro_e_clipped[c_id] = epz2 - cvar_ep[c_id];
        cvar_ep[c_id] = cs::max(cvar_ep[c_id], epz2);
      }
      else if(xe <= 0.) {
        res.i[1] = 1;
        if (clip_e_id >= 0)
          cpro_e_clipped[c_id] = - xe;
        cvar_ep[c_id] = - xe;
      }
    });

    ctx.wait();

  }

  cs_lnum_t iclpmx[1] = {0};
  int id;

  for (int ii = 0; ii < 2; ii++ ) {
    if (ii == 0)
      id = f_k->id;
    else if (ii == 1)
      id = f_eps->id;

    cs_log_iteration_clipping_field(id,
                                    rd_sum.i[ii],
                                    0,
                                    rd.r + 2*ii,
                                    rd.r + 2*ii + 1,
                                    rd_sum.i + ii,
                                    iclpmx);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of turbulent viscosity for the K-epsilon model.
 *
 * \param[in]     phase_id   turbulent phase id (-1 for single phase flow)
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_ke_mu_t(int  phase_id)
{
  cs_dispatch_context ctx;

  /* Map field arrays */

  cs_field_t *f_k = CS_F_(k);
  cs_field_t *f_eps = CS_F_(eps);
  cs_field_t *f_rho = CS_F_(rho);
  cs_field_t *f_mu = CS_F_(mu);
  cs_field_t *f_mut = CS_F_(mu_t);

  if (phase_id >= 0) {
    f_k = CS_FI_(k, phase_id);
    f_eps = CS_FI_(eps, phase_id);
    f_rho = CS_FI_(rho, phase_id);
    f_mu = CS_FI_(mu, phase_id);
    f_mut = CS_FI_(mu_t, phase_id);
  }

  const cs_real_t cmu = cs_turb_cmu;

  cs_real_t *visct = f_mut->val;
  const cs_real_t *viscl =  f_mu->val;
  const cs_real_t *crom = f_rho->val;
  const cs_real_t *cvar_k = f_k->val;
  const cs_real_t *cvar_ep = f_eps->val;

  /* Launder-Sharma */

  if (cs_glob_turb_model->model == CS_TURB_K_EPSILON_LS) {

    const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

    /* Initialization
     * ============== */

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      const cs_real_t xk = cvar_k[c_id];
      const cs_real_t xe = cvar_ep[c_id];
      const cs_real_t xmu = viscl[c_id];
      const cs_real_t xmut = crom[c_id] * cs_math_pow2(xk) / xe;
      const cs_real_t xfmu = exp(-3.4 / cs_math_pow2(1. + xmut/xmu/50.));

      visct[c_id] = cmu * xfmu * xmut;
    });

  }

  /* Baglietto model */

  else if (cs_glob_turb_model->model == CS_TURB_K_EPSILON_QUAD) {
    cs_turbulence_ke_q_mu_t(phase_id);
  }

  /* Standard and linear-production */

  else {

    const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      visct[c_id] =   crom[c_id] * cmu
                    * cs_math_pow2(cvar_k[c_id]) / cvar_ep[c_id];
    });

  }

  ctx.wait();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of turbulent viscosity for
 *        the non-linear quadratic K-epsilon from
 *        Baglietto et al. (2005)
 *
 * \param[in]     phase_id   turbulent phase id (-1 for single phase flow)
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_ke_q_mu_t(int phase_id)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  cs_dispatch_context ctx;

  /* Initialization
   * ============== */

  /* Map field arrays */

  cs_field_t *f_k = CS_F_(k);
  cs_field_t *f_eps = CS_F_(eps);
  cs_field_t *f_vel = CS_F_(vel);
  cs_field_t *f_rho = CS_F_(rho);
  cs_field_t *f_mu = CS_F_(mu);
  cs_field_t *f_mut = CS_F_(mu_t);

  if (phase_id >= 0) {
    f_k = CS_FI_(k, phase_id);
    f_eps = CS_FI_(eps, phase_id);
    f_vel = CS_FI_(vel, phase_id);
    f_rho = CS_FI_(rho, phase_id);
    f_mu = CS_FI_(mu, phase_id);
    f_mut = CS_FI_(mu_t, phase_id);
  }

  cs_real_t *visct = f_mut->val;
  const cs_real_t *viscl =  f_mu->val;
  const cs_real_t *crom = f_rho->val;
  const cs_real_t *cvar_k = f_k->val;
  const cs_real_t *cvar_ep = f_eps->val;

  const cs_real_t *w_dist = cs_field_by_name("wall_distance")->val;

  /* Calculation of velocity gradient and of
   *   S2 = S11**2+S22**2+S33**2+2*(S12**2+S13**2+S23**2)
   * ==================================================== */

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

  /* Compute turbulent viscosity
   * =========================== */

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    const cs_real_t s11 = gradv[c_id][0][0];
    const cs_real_t s22 = gradv[c_id][1][1];
    const cs_real_t s33 = gradv[c_id][2][2];
    const cs_real_t dudy = gradv[c_id][0][1];
    const cs_real_t dudz = gradv[c_id][0][2];
    const cs_real_t dvdx = gradv[c_id][1][0];
    const cs_real_t dvdz = gradv[c_id][1][2];
    const cs_real_t dwdx = gradv[c_id][2][0];
    const cs_real_t dwdy = gradv[c_id][2][1];

    const cs_real_t s2 =   cs_math_pow2(s11) + cs_math_pow2(s22) + cs_math_pow2(s33)
                         + 0.5*cs_math_pow2(dudy+dvdx)
                         + 0.5*cs_math_pow2(dudz+dwdx)
                         + 0.5*cs_math_pow2(dvdz+dwdy);

    const cs_real_t xk = cvar_k[c_id];
    const cs_real_t xe = cvar_ep[c_id];
    const cs_real_t xrom = crom[c_id];
    const cs_real_t xmu = viscl[c_id];
    const cs_real_t xdist = fmax(w_dist[c_id], 1.e-10);

    const cs_real_t xmut = xrom*cs_math_pow2(xk)/xe;
    const cs_real_t xrey = xdist*sqrt(xk)*xrom/xmu;
    const cs_real_t xttke = xk/xe;
    const cs_real_t xss = xttke*sqrt(2.0*s2);

    const cs_real_t xfmu = 1.0 - exp(- 2.9e-2*sqrt(xrey)
                                     - 1.1e-4*cs_math_pow2(xrey));
    const cs_real_t xcmu = 2.0 / 3.0 / (3.90 + xss);

    visct[c_id] = xcmu*xfmu*xmut;
  });

  ctx.wait();

  CS_FREE_HD(_gradv);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of non linear terms of the quadratic k-epsilon model
 *        (Baglietto et al.)
 *
 * \param[in]   phase_id  turbulent phase id (-1 for single phase flow)
 * \param[out]  rij       non linear terms of quadratic Boussinesq approximation
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_ke_q(int          phase_id,
                   cs_real_6_t  rij[])
{
  cs_dispatch_context ctx;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  cs_field_t *f_k = CS_F_(k);
  cs_field_t *f_eps = CS_F_(eps);
  cs_field_t *f_vel = CS_F_(vel);
  cs_field_t *f_mut = CS_F_(mu_t);

  if (phase_id >= 0) {
    f_k = CS_FI_(k, phase_id);
    f_eps = CS_FI_(eps, phase_id);
    f_vel = CS_FI_(vel, phase_id);
    f_mut = CS_FI_(mu_t, phase_id);
  }

  const cs_real_t *visct   = f_mut->val;
  const cs_real_t *cvar_k  = f_k->val;
  const cs_real_t *cvar_ep = f_eps->val;

  /* Initialization
   * ============== */

  cs_real_33_t *gradv = nullptr, *_gradv = nullptr;

  if (f_vel->grad == nullptr) {
    CS_MALLOC(_gradv, n_cells_ext, cs_real_33_t);
    gradv = _gradv;

    /* Computation of the velocity gradient */

    cs_field_gradient_vector(f_vel,
                             true,  // use_previous_t
                             1,     // inc
                             gradv);
  }
  else
    gradv = (cs_real_33_t *)f_vel->grad;

  const cs_real_t cnl1 = cs_turb_cnl1;
  const cs_real_t cnl2 = cs_turb_cnl2;
  const cs_real_t cnl3 = cs_turb_cnl3;
  const cs_real_t cnl4 = cs_turb_cnl4;
  const cs_real_t cnl5 = cs_turb_cnl5;

  const cs_real_t d1s2 = 0.5, d2s3 = 2./3;

  /*  Computation
   *===============*/

  ctx.parallel_for(n_cells_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    cs_real_t xrij[3][3];
    cs_real_t xstrai[3][3], xrotac[3][3], sikskj[3][3];
    cs_real_t wikskj[3][3], skiwjk[3][3], wikwjk[3][3];

    const cs_real_t xvisct = visct[c_id];
    const cs_real_t xeps   = cvar_ep[c_id];
    const cs_real_t xk     = cvar_k[c_id];
    const cs_real_t xttke  = xk/xeps;

    /* Sij */
    xstrai[0][0] = gradv[c_id][0][0];
    xstrai[0][1] = d1s2*(gradv[c_id][0][1]+gradv[c_id][1][0]);
    xstrai[0][2] = d1s2*(gradv[c_id][0][2]+gradv[c_id][2][0]);
    xstrai[1][0] = xstrai[0][1];
    xstrai[1][1] = gradv[c_id][1][1];
    xstrai[1][2] = d1s2*(gradv[c_id][1][2]+gradv[c_id][2][1]);
    xstrai[2][0] = xstrai[0][2];
    xstrai[2][1] = xstrai[1][2];
    xstrai[2][2] = gradv[c_id][2][2];

    /* omegaij */
    xrotac[0][0] = 0;
    xrotac[0][1] = d1s2*(-gradv[c_id][0][1]+gradv[c_id][1][0]);
    xrotac[0][2] = d1s2*(-gradv[c_id][0][2]+gradv[c_id][2][0]);
    xrotac[1][0] = -xrotac[0][1];
    xrotac[1][1] = 0;
    xrotac[1][2] = d1s2*(-gradv[c_id][1][2]+gradv[c_id][2][1]);
    xrotac[2][0] = -xrotac[0][2];
    xrotac[2][1] = -xrotac[1][2];
    xrotac[2][2] = 0;

    cs_real_t sijsij = 0;
    for (cs_lnum_t ii = 0; ii < 3; ii++) {
      for (cs_lnum_t jj = 0; jj < 3; jj++) {
        sijsij += xstrai[jj][ii]*xstrai[jj][ii];
        sikskj[jj][ii] = 0.0;
        wikskj[jj][ii] = 0.0;
        skiwjk[jj][ii] = 0.0;
        wikwjk[jj][ii] = 0.0;
        for (cs_lnum_t kk = 0; kk < 3; kk++) {
          sikskj[jj][ii] += xstrai[kk][ii]*xstrai[jj][kk];
          wikskj[jj][ii] += xrotac[kk][ii]*xstrai[jj][kk];
          skiwjk[jj][ii] += xstrai[ii][kk]*xrotac[kk][jj];
          wikwjk[jj][ii] += xrotac[kk][ii]*xrotac[kk][jj];
        }
      }
    }

    const cs_real_t xss = xttke*sqrt(2.*sijsij);
    const cs_real_t xcmu = d2s3/(3.9 + xss);

    /* Evaluating "constants". */
    const cs_real_t xqc1
      = cnl1/((cnl4+cnl5*cs_math_pow3(xss))*xcmu);
    const cs_real_t xqc2
      = cnl2/((cnl4+cnl5*cs_math_pow3(xss))*xcmu);
    const cs_real_t xqc3
      = cnl3/((cnl4+cnl5*cs_math_pow3(xss))*xcmu);

    for (cs_lnum_t ii = 0; ii < 3; ii++) {
      for (cs_lnum_t jj = 0; jj < 3; jj++) {
        xrij[jj][ii] =   4*xqc1*xvisct*xttke*sikskj[jj][ii]
                       + 4*xqc2*xvisct*xttke*(wikskj[jj][ii]+skiwjk[jj][ii])
                       + 4*xqc3*xvisct*xttke*wikwjk[jj][ii];
      }
    }

    rij[c_id][0] = xrij[0][0];
    rij[c_id][1] = xrij[1][1];
    rij[c_id][2] = xrij[2][2];
    rij[c_id][3] = xrij[1][0];
    rij[c_id][4] = xrij[2][1];
    rij[c_id][5] = xrij[2][0];
  });

  ctx.wait();

  /* Free memory */
  CS_FREE(_gradv);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
