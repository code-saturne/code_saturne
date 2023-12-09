/*============================================================================
 * k-epsilon turbulence model.
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
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_array.h"
#include "cs_balance.h"
#include "cs_blas.h"
#include "cs_convection_diffusion.h"
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
#include "cs_parall.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_porous_model.h"
#include "cs_prototypes.h"
#include "cs_rotation.h"
#include "cs_solid_zone.h"
#include "cs_thermal_model.h"
#include "cs_time_step.h"
#include "cs_turbulence_model.h"
#include "cs_turbulence_rotation.h"
#include "cs_velocity_pressure.h"
#include "cs_wall_functions.h"

/* Atmospheric model headers
   (we should not need to call them here, we should be more modular */
#include "cs_at_data_assim.h"
#include "cs_atprke.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_turbulence_ke.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_turbulence_ke.c

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
 *  \param[in,out] w1            work array to store the E-term
 */
/*----------------------------------------------------------------------------*/

static void
_tsepls(int       phase_id,
        cs_real_t w1[])
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t    n_b_faces = m->n_b_faces;
  const cs_lnum_t    n_i_faces = m->n_i_faces;
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
  BFT_MALLOC(w7, n_cells_ext, cs_real_33_t);

  /* Calculation of the term d2Ui/dxkdxj*d2Ui/dxkdxj
   * ================================================ */

  cs_field_t *f_vel = CS_F_(vel);

  if (phase_id >= 0)
    f_vel = CS_FI_(vel, phase_id);

  cs_real_33_t *gradv = NULL;

  if (f_vel->grad == NULL) {
    BFT_MALLOC(gradv, n_cells_ext, cs_real_33_t);

    /* Computation of the velocity gradient */

    cs_field_gradient_vector(f_vel,
                             true,  /* use_previous_t */
                             1,     /* inc */
                             gradv);
  } else
    gradv = (cs_real_33_t *)f_vel->grad;

  /* Loop over u, v, w components:
     TODO interleave */

  for (cs_lnum_t isou = 0; isou < 3; isou++) {

    cs_array_real_fill_zero(9*n_cells_ext, (cs_real_t *)w7);

#   pragma omp parallel for if(n_i_faces > CS_THR_MIN)
    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

      cs_real_t duidxk[3], njsj[3];

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];
      cs_real_t pnd = weight[face_id];

      duidxk[0] = pnd * gradv[ii][isou][0] + (1. - pnd) * gradv[jj][isou][0];
      duidxk[1] = pnd * gradv[ii][isou][1] + (1. - pnd) * gradv[jj][isou][1];
      duidxk[2] = pnd * gradv[ii][isou][2] + (1. - pnd) * gradv[jj][isou][2];
      njsj[0]   = i_face_normal[face_id][0];
      njsj[1]   = i_face_normal[face_id][1];
      njsj[2]   = i_face_normal[face_id][2];
      for (cs_lnum_t k = 0; k < 3; k++){
        for (cs_lnum_t j = 0; j < 3; j++){
          w7[ii][j][k] += duidxk[k]*njsj[j];
          w7[jj][j][k] -= duidxk[k]*njsj[j];
        }
      }

    }

#   pragma omp parallel for if(n_b_faces > CS_THR_MIN)
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

      cs_real_t duidxk[3], njsj[3];

      cs_lnum_t ii = b_face_cells[face_id];

      duidxk[0] = gradv[ii][isou][0];
      duidxk[1] = gradv[ii][isou][1];
      duidxk[2] = gradv[ii][isou][2];
      njsj[0]   = b_face_normal[face_id][0];
      njsj[1]   = b_face_normal[face_id][1];
      njsj[2]   = b_face_normal[face_id][2];
      for (cs_lnum_t k = 0; k < 3; k++){
        for (cs_lnum_t j = 0; j < 3; j++){
          w7[ii][j][k] += duidxk[k]*njsj[j];
        }
      }

    }

#   pragma omp parallel for if(n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      cs_real_t w_temp = 0.;
      for (cs_lnum_t k = 0; k < 3; k++){
        for (cs_lnum_t j = 0; j < 3; j++){
          w_temp += cs_math_pow2(w7[c_id][j][k] / volume[c_id]);
        }
      }
      w1[c_id] += w_temp;
    }

  }

  /* Free memory */
  if (f_vel->grad == NULL)
    BFT_FREE(gradv);

  BFT_FREE(w7);
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
 * \param[in]     ncesmp        number of cells with mass source term
 * \param[in]     icetsm        index of cells with mass source term
 * \param[in]     itypsm        mass source type for the variables
 *                              size: [nvar][ncesmp]
 * \param[in]     dt            time step (per cell)
 * \param[in]     smacel        values of the variables associated to the
 *                              mass source (for the pressure variable,
 *                              smacel is the mass flux)
 *                              size: [nvar][ncesmp]
 * \param[out]    prdv2f        v2f production term
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_ke(int              phase_id,
                 cs_lnum_t        ncesmp,
                 cs_lnum_t        icetsm[],
                 int              itypsm[],
                 const cs_real_t  dt[],
                 cs_real_t        smacel[],
                 cs_real_t       *prdv2f)
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
  cs_real_t viscl0 = phys_pro->viscl0; /* reference pressure */
  cs_real_t ro0 = phys_pro->ro0; /* reference density */
  const cs_real_t uref = cs_glob_turb_ref_values->uref;
  const cs_real_t *grav = cs_glob_physical_constants->gravity;

  const int var_key_id = cs_field_key_id("variable_id");

  cs_field_t *f_k = CS_F_(k);
  cs_field_t *f_eps = CS_F_(eps);
  cs_field_t *f_phi = CS_F_(phi);
  cs_field_t *f_alpbl = CS_F_(alp_bl);

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

  cs_real_t *cpro_sigmak = NULL;
  if (sigmak_id >= 0)
    cpro_sigmak = cs_field_by_id(sigmak_id)->val;

  cs_real_t sigmae = cs_field_get_key_double(f_eps,
                                             cs_field_key_id("turbulent_schmidt"));

  /* If turbulent Schmidt is variable, id of the corresponding field */
  int sigmae_id = cs_field_get_key_int(f_eps,
                                       cs_field_key_id("turbulent_schmidt_id"));

  cs_real_t *cpro_sigmae = NULL;
  if (sigmae_id >= 0)
    cpro_sigmae = cs_field_by_id(sigmae_id)->val;

  /* Initialization of work arrays in case of Hybrid turbulence modelling
   * ==================================================================== */

  cs_real_t *htles_psi       = NULL;
  cs_real_t *htles_t         = NULL;
  cs_real_t *hybrid_fd_coeff = NULL;
  if (cs_glob_turb_model->hybrid_turb == 4) {
    htles_psi = cs_field_by_name("htles_psi")->val;
    htles_t   = cs_field_by_name("htles_t")->val;
    hybrid_fd_coeff = cs_field_by_name("hybrid_blend")->val;
  }

  /* Initialilization
     ================ */

  /* Allocate temporary arrays for the turbulence resolution */

  cs_real_t *viscf, *viscb;
  BFT_MALLOC(viscf, n_i_faces, cs_real_t);
  BFT_MALLOC(viscb, n_b_faces, cs_real_t);

  cs_real_t *smbrk, *smbre;
  BFT_MALLOC(smbrk, n_cells_ext, cs_real_t);
  BFT_MALLOC(smbre, n_cells_ext, cs_real_t);

  cs_real_t *tinstk, *tinste, *divu, *strain_sq;
  BFT_MALLOC(tinstk, n_cells_ext, cs_real_t);
  BFT_MALLOC(tinste, n_cells_ext, cs_real_t);
  BFT_MALLOC(divu, n_cells_ext, cs_real_t);
  BFT_MALLOC(strain_sq, n_cells_ext, cs_real_t);

  /* Allocate work arrays */

  cs_real_t *w1, *w2, *w3, *w4, *w5;

  BFT_MALLOC(w1, n_cells_ext, cs_real_t);
  BFT_MALLOC(w2, n_cells_ext, cs_real_t);
  BFT_MALLOC(w3, n_cells_ext, cs_real_t);
  BFT_MALLOC(w4, n_cells_ext, cs_real_t);
  BFT_MALLOC(w5, n_cells_ext, cs_real_t);

  cs_real_t *usimpk, *w7, *w8, *usimpe, *dpvar;
  BFT_MALLOC(usimpk, n_cells_ext, cs_real_t);
  BFT_MALLOC(w7, n_cells_ext, cs_real_t);
  BFT_MALLOC(w8, n_cells_ext, cs_real_t);
  BFT_MALLOC(usimpe, n_cells_ext, cs_real_t);
  BFT_MALLOC(dpvar, n_cells_ext, cs_real_t);

  cs_real_t *prdtke = NULL, *prdeps = NULL, *sqrt_k = NULL, *strain = NULL;
  cs_real_3_t *grad_sqk = NULL, *grad_s = NULL;
  cs_real_t *coefa_sqk = NULL, *coefb_sqk = NULL;
  cs_real_t *coefa_sqs = NULL, *coefb_sqs = NULL;
  cs_real_t *w10 = NULL, *w11 = NULL;

  if (cs_glob_turb_model->iturb == CS_TURB_K_EPSILON) {
    BFT_MALLOC(prdtke, n_cells_ext, cs_real_t);
    BFT_MALLOC(prdeps, n_cells_ext, cs_real_t);
  }
  else if (cs_glob_turb_model->iturb == CS_TURB_K_EPSILON_LS) {
    BFT_MALLOC(sqrt_k, n_cells_ext, cs_real_t);
    BFT_MALLOC(strain, n_cells_ext, cs_real_t);
    BFT_MALLOC(grad_sqk, n_cells_ext, cs_real_3_t);
    BFT_MALLOC(grad_s, n_cells_ext, cs_real_3_t);
    BFT_MALLOC(coefa_sqk, n_b_faces, cs_real_t);
    BFT_MALLOC(coefb_sqk, n_b_faces, cs_real_t);
    BFT_MALLOC(coefa_sqs, n_b_faces, cs_real_t);
    BFT_MALLOC(coefb_sqs, n_b_faces, cs_real_t);
  }
  else if (cs_glob_turb_model->iturb == CS_TURB_V2F_BL_V2K) {
    BFT_MALLOC(w10, n_cells_ext, cs_real_t);
    BFT_MALLOC(w11, n_cells_ext, cs_real_t);
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
  cs_real_t *cvara_phi = NULL;
  if (   cs_glob_turb_model->iturb == CS_TURB_V2F_PHI
      || cs_glob_turb_model->iturb == CS_TURB_V2F_BL_V2K){
    cvara_phi = f_phi->val_pre;
  }
  cs_real_t *cvara_al = NULL;
  if (cs_glob_turb_model->iturb == CS_TURB_V2F_BL_V2K) {
    cvara_al = f_alpbl->val_pre;
  }
  const cs_real_t *w_dist = NULL;
  if (cs_glob_turb_model->iturb == CS_TURB_K_EPSILON_QUAD) {
    w_dist =  cs_field_by_name("wall_distance")->val;
  }

  int kstprv = cs_field_key_id("source_term_prev_id");
  int istprv =  cs_field_get_key_int(f_k, kstprv);
  cs_real_t *c_st_k_p = NULL;
  cs_real_t *c_st_eps_p = NULL;
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
    if (cs_glob_turb_model->iturb == CS_TURB_K_EPSILON) {
      cs_log_printf(CS_LOG_DEFAULT,
                    "\n"
                    "  ** Solving k-epsilon\n"
                    "     -----------------\n");
    }
    else if (cs_glob_turb_model->iturb == CS_TURB_K_EPSILON_LIN_PROD) {
      cs_log_printf(CS_LOG_DEFAULT,
                    "\n"
                    "  ** Solving k-epsilon with linear production\n"
                    "     -----------------\n");
    }
    else if (cs_glob_turb_model->iturb == CS_TURB_K_EPSILON_LS) {
      cs_log_printf(CS_LOG_DEFAULT,
                    "\n"
                    "  ** Solving k-epsilon, Launder-Sharma\n"
                    "     -----------------\n");
    }
    else {
      cs_log_printf(CS_LOG_DEFAULT,
                    "\n"
                    "  ** Solving v2f (k and epsilon)\n"
                    "     ---------------------------\n");
    }
  }

  /* For the model with linear production, sqrt(Cmu) is required */
  const cs_real_t sqrcmu = sqrt(cs_turb_cmu);

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
      && cs_glob_turb_model->iturb == CS_TURB_V2F_BL_V2K) {

    cs_real_t *cvar_al = (cs_real_t *)f_alpbl->val;

#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      /* y+ is bounded by 400, because in the Reichard profile,
         it corresponds to saturation (u>uref) */
      cvar_al[c_id] = fmax(fmin(cvar_al[c_id], 1. - exp(-400./50.)),
                           0.);
    }

    cs_field_current_to_previous(f_alpbl);

    cs_real_3_t *grad = NULL;

    if (f_alpbl->grad == NULL) {
      BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);

    /* Compute the gradient of Alpha */

    cs_field_gradient_scalar(f_alpbl,
                             false,  /* use_previous_t */
                             1,      /* inc */
                             grad);
    } else
      grad = (cs_real_3_t *)f_alpbl->grad;

    cvar_ep = f_eps->val;
    cs_real_3_t *vel = (cs_real_3_t *)f_vel->val;

    cs_real_t utaurf = 0.05*uref;
    cs_real_t nu0 = viscl0 / ro0;

#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

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
        limiter = fmin(utaurf/xunorm*(2.5*log(1. + 0.4*ypa)
                                      +7.8*(1. - exp(-ypa/11.)
                                            -(ypa/11.)*exp(-0.33*ypa))),
                       1.);
      }

      vel[c_id][0] = limiter*vel[c_id][0];
      vel[c_id][1] = limiter*vel[c_id][1];
      vel[c_id][2] = limiter*vel[c_id][2];

      cs_real_t ut2 = 0.05*uref;

      cvar_ep[c_id] =   cs_math_pow3(utaurf)
                      * fmin(1./(cs_turb_xkappa*15.*nu0/utaurf),
                             1./(cs_turb_xkappa*ya));
      cvar_k[c_id] =     cvar_ep[c_id]/2./nu0*cs_math_pow2(ya)
                       * cs_math_pow2(exp(-ypa/25.))
                     +   cs_math_pow2(ut2)/sqrt(cs_turb_cmu)
                       * cs_math_pow2(1.-exp(-ypa/25.));
    }

    cs_field_current_to_previous(f_vel);
    cs_field_current_to_previous(f_k);
    cs_field_current_to_previous(f_eps); /*TODO phi ? */

    if (f_alpbl->grad == NULL)
      BFT_FREE(grad);
  }

  /* Compute the scalar strain rate squared S2 =2SijSij and the trace of
     the velocity gradient
     (Sij^D) (Sij^D)  is stored in    strain_sq (deviatoric strain tensor rate)
     tr(Grad u)       is stored in    divu
     ======================================================================= */

  /* Allocate temporary arrays for gradients calculation */
  cs_real_33_t *gradv = NULL;

  if (f_vel->grad == NULL) {
    BFT_MALLOC(gradv, n_cells_ext, cs_real_33_t);

    /* Computation of the velocity gradient */

    cs_field_gradient_vector(f_vel,
                             true,  /* use_previous_t */
                             1,     /* inc */
                             gradv);
  } else
    gradv = (cs_real_33_t *)f_vel->grad;

  /* strain_sq = Stain rate of the deviatoric part of the strain tensor
     = 2 (Sij^D).(Sij^D)
     divu   = trace of the velocity gradient
     = dudx + dvdy + dwdz */

# pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

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

    divu[c_id] = gradv[c_id][0][0] + gradv[c_id][1][1] + gradv[c_id][2][2];
  }

  /* Compute the square root of the strain and the turbulent
     kinetic energy for Launder-Sharma k-epsilon source terms */
  if (cs_glob_turb_model->iturb == CS_TURB_K_EPSILON_LS) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      strain[c_id] = sqrt(strain_sq[c_id]);
      sqrt_k[c_id] = sqrt(fabs(cvar_k[c_id]));
    }
  }

  /* Unsteady terms (stored in tinstk and tinste)
     ============================================ */

# pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_real_t romvsd = crom[c_id] * cell_f_vol[c_id] / dt[c_id];
    tinstk[c_id] = eqp_k->istat*romvsd;
    tinste[c_id] = eqp_eps->istat*romvsd;
  }

  /* Compute the first part of the production term: muT (S^D)**2
   * Going out of the step we keep strain_sq, divu,
   *============================================================
   *
   * For the Linear Production k-epsilon model,
   * the production term is assumed to be asymptotiy in S and
   * not in mu_TxS**2 */

  cs_field_t *f_tke_prod = cs_field_by_name_try("tke_production");

  if (cs_glob_turb_model->iturb == CS_TURB_K_EPSILON_LIN_PROD) {
#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      cs_real_t rho    = cromo[c_id];
      cs_real_t xs     = sqrt(strain_sq[c_id]);
      cs_real_t cmueta = fmin(cs_turb_cmu*cvara_k[c_id]/cvara_ep[c_id]*xs,
                              sqrcmu);
      smbrk[c_id] = rho*cmueta*xs*cvara_k[c_id];
      smbre[c_id] = smbrk[c_id];
    }
    /* Save production for post processing */
    if (f_tke_prod != NULL) {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        f_tke_prod->val[c_id] = smbrk[c_id] / cromo[c_id];
    }
  }
  else if (cs_glob_turb_model->iturb == CS_TURB_K_EPSILON_QUAD) {

    /* Turbulent production for the quadratic Baglietto k-epsilon model  */

#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      cs_real_t xstrai[3][3], xrotac[3][3];

      cs_real_t visct = cpro_pcvto[c_id];
      cs_real_t xeps  = cvar_ep[c_id];
      cs_real_t xk    = cvar_k[c_id];
      cs_real_t xttke = xk/xeps;

      /* Sij */
      xstrai[0][0] = gradv[c_id][0][0];
      xstrai[1][0] = d1s2*(gradv[c_id][0][1]+gradv[c_id][1][0]);
      xstrai[2][0] = d1s2*(gradv[c_id][0][2]+gradv[c_id][2][0]);
      xstrai[0][1] = xstrai[1][0];
      xstrai[1][1] = gradv[c_id][1][1];
      xstrai[2][1] = d1s2*(gradv[c_id][1][2]+gradv[c_id][2][1]);
      xstrai[0][2] = xstrai[2][0];
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
      cs_real_t xqc1 = cs_turb_cnl1/((  cs_turb_cnl4
                                      + cs_turb_cnl5*cs_math_pow3(xss))*xcmu);
      cs_real_t xqc2 = cs_turb_cnl2/((  cs_turb_cnl4
                                      + cs_turb_cnl5*cs_math_pow3(xss))*xcmu);
      cs_real_t xqc3 = cs_turb_cnl3/((  cs_turb_cnl4
                                      + cs_turb_cnl5*cs_math_pow3(xss))*xcmu);

      /* Evaluating the turbulent production */
      smbrk[c_id] =   visct*strain_sq[c_id]
                    - 4.*xqc1*visct*xttke* (skskjsji - d1s3*sijsij*divu[c_id])
                    - 4.*xqc2*visct*xttke* (wkskjsji + skiwjksji)
                    - 4.*xqc3*visct*xttke* (wkwjksji - d1s3*wijwij*divu[c_id]);
      smbre[c_id] = smbrk[c_id];
    } /* End loop on cells */

    /* Save production for post processing */
    if (f_tke_prod != NULL) {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        f_tke_prod->val[c_id] = smbrk[c_id] / crom[c_id];
    }

    /* End test on specific k-epsilon model
       In the general case Pk = mu_t*SijSij */
  }
  else {
#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      smbrk[c_id] = cpro_pcvto[c_id] * strain_sq[c_id];
      smbre[c_id] = smbrk[c_id];
    }
    /* Save production for post processing */
    if (f_tke_prod != NULL) {
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        f_tke_prod->val[c_id] = smbrk[c_id] / crom[c_id];
    }
  }

  /* Take into account rotation/curvature correction, if necessary
     ============================================================ */

  /* Cazalbou correction: the Ceps2 coefficient of destruction term of epsislon
     is modified by rotation and curvature */

  /* Allocate an array for the modified Ceps2 coefficient */
  cs_real_t *ce1rc, *ce2rc;
  BFT_MALLOC(ce1rc, n_cells, cs_real_t);
  BFT_MALLOC(ce2rc, n_cells, cs_real_t);

  if (cs_glob_turb_rans_model->irccor == 1) {

    /* Compute the modified Ceps2 coefficient (w1 array not used) */
    /* cs_turbulence_rotation_correction(dt, w1, ce2rc); */

    cs_turbulence_rotation_correction(dt, w1, ce2rc);

  }
  else {
    if (cs_glob_turb_model->itytur == 2) {

      /* Launder-Sharma k-epsilon model */
      if (cs_glob_turb_model->iturb == CS_TURB_K_EPSILON_LS) {
#       pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_real_t rho  = crom[c_id];
          cs_real_t xeps = cvar_ep[c_id];
          cs_real_t xk   = cvar_k[c_id];
          ce2rc[c_id] = (1. - 0.3*exp(-cs_math_pow2( rho*cs_math_pow2(xk)
                                                    /viscl[c_id]
                                                    /xeps)))*cs_turb_ce2;
          ce1rc[c_id] = cs_turb_ce1;
        }
      }

      /* Baglietto quadratic k-epsilon model */
      else if (cs_glob_turb_model->iturb == CS_TURB_K_EPSILON_QUAD) {
#       pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_real_t rho  = crom[c_id];
          cs_real_t xeps = cvar_ep[c_id];
          cs_real_t xk   = cvar_k[c_id];
          ce2rc[c_id] = (1. - 0.3*exp(-cs_math_pow2( rho*cs_math_pow2(xk)
                                                    /viscl[c_id]
                                                    /xeps)))*cs_turb_ce2;

          cs_real_t xdist = fmax(w_dist[c_id], cs_math_epzero);
          cs_real_t xrey  = rho*xdist*sqrt(xk)/viscl[c_id];
          cs_real_t xpk   = smbrk[c_id] - d2s3*rho*xk*divu[c_id];
          cs_real_t xpkp
            =   1.33*(1. - 0.3*exp(-cs_math_pow2(rho*cs_math_pow2(xk)
                                                 /viscl[c_id]/xeps)))
              * (xpk + 2.*viscl[c_id]*xk/cs_math_pow2(xdist))
              * exp(-3.75e-3*cs_math_pow2(xrey));

          ce1rc[c_id] = (1. + xpkp/fmax(xpk, 1.e-10))*cs_turb_ce1;
        }
      }

      /* Other k-epsilon models */
      else {
#       pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          ce2rc[c_id] = cs_turb_ce2;
          ce1rc[c_id] = cs_turb_ce1;
        }
      }
    }

    else if (cs_glob_turb_model->iturb == CS_TURB_V2F_PHI) {
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        ce2rc[c_id] = cs_turb_cv2fe2;
        ce1rc[c_id] = cs_turb_ce1;
      }
    }
    else if (cs_glob_turb_model->iturb == CS_TURB_V2F_BL_V2K) {
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        ce2rc[c_id] = cs_turb_ccaze2;
        ce1rc[c_id] = cs_turb_ce1;
      }
    }
  }

  /* ce2rc array is used all along this function, and deallocated at the end. */

  /* Compute the buoyancy term
   * The mass sources receive production and gravity terms
   * Work arrays                      viscb
   * The mass sources are stored in   smbrk, smbre
   * Going out of the step we keep    smbrk, smbre,
   * divu,
   ======================================================== */

  /* Buoyant term for the Atmospheric module
     (function of the potential temperature) */
  if (   cs_glob_turb_rans_model->igrake == 1
      && cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 1) {

    cs_atmo_buoyancy_ke_prod(tinstk, smbrk, smbre);

    /* --- Buoyancy term     G = Beta*g.Grad(scalar)/prdtur/rho
       Here is computed  G =-g.grad(rho)/prdtur/rho */
  }
  else if (cs_glob_turb_rans_model->igrake == 1) {

    /* Allocate a temporary for the gradient calculation */
    cs_real_3_t *grad;
    cs_real_t  *grad_dot_g;
    BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);
    BFT_MALLOC(grad_dot_g, n_cells_ext, cs_real_t);

    cs_real_t prdtur = 1;

    cs_field_t *f_thm = cs_thermal_model_field();
    if (phase_id >= 0)
      f_thm = CS_FI_(h_tot, phase_id);

    if (f_thm != NULL) {
      prdtur = cs_field_get_key_double(f_thm,
                                       cs_field_key_id("turbulent_schmidt"));
    }

    /* Boussinesq approximation, only for the thermal scalar for the moment */

    if (cs_glob_velocity_pressure_model->idilat == 0) {

      cs_field_gradient_scalar(f_thm,
                               true,
                               1,    /* inc */
                               grad);

      /* FIXME make it dependant on the scalar and use is_buoyant field */
      cs_real_t *cpro_beta = cs_field_by_name("thermal_expansion")->val;

      /* - Beta grad(T) . g / Pr_T */
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        grad_dot_g[c_id]
          = - cpro_beta[c_id] * (  cs_math_3_dot_product(grad[c_id], grav)
                                 / prdtur);
      }

    }
    else {

      /* BCs on rho: Dirichlet ROMB
         NB: viscb is used as COEFB */

      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
        viscb[face_id] = 0.;
      }

      cs_halo_type_t halo_type = CS_HALO_STANDARD;
      cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

      cs_gradient_type_by_imrgra(eqp_k->imrgra,
                                 &gradient_type,
                                 &halo_type);

      cs_gradient_scalar("cromo_grad",
                         gradient_type,
                         halo_type,
                         1,     /* inc */
                         eqp_k->nswrgr,
                         0,
                         1,     /* w_stride */
                         eqp_k->verbosity,
                         eqp_k->imligr,
                         eqp_k->epsrgr,
                         eqp_k->climgr,
                         NULL,
                         bromo,
                         viscb,
                         cromo,
                         NULL,
                         NULL, /* internal coupling */
                         grad);

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        grad_dot_g[c_id] =   cs_math_3_dot_product(grad[c_id], grav)
                           / (cromo[c_id]*prdtur);
      }
    }

    /* Production term due to buoyancy
       smbrk = P+G
       smbre = P+(1-ce3)*G */

    cs_field_t *f_tke_buoy = cs_field_by_name_try("tke_buoyancy");

    /* smbr* store mu_TxS**2 */
#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      cs_real_t rho   = cromo[c_id];
      cs_real_t visct = cpro_pcvto[c_id];
      cs_real_t xeps  = cvara_ep[c_id];
      cs_real_t xk    = cvara_k[c_id];
      cs_real_t ttke  = xk / xeps;

      /* Implicit Buoyant terms when negative */
      tinstk[c_id] += fmax(rho*cell_f_vol[c_id]*cs_turb_cmu*ttke*grad_dot_g[c_id], 0.);

      /* Explicit Buoyant terms */
      smbre[c_id] = smbre[c_id] + visct*fmax(- grad_dot_g[c_id], 0.);
      smbrk[c_id] = smbrk[c_id] - visct*grad_dot_g[c_id];
      /* Save for post processing */
      if (f_tke_buoy != NULL)
        f_tke_buoy->val[c_id] = -visct*grad_dot_g[c_id]/rho;
    }

    /* Free memory */
    BFT_FREE(grad);
    BFT_FREE(grad_dot_g);

  }

  /* In v2f, we store the production in prdv2f which will be complete further
     for containing the complete production term */
  if (cs_glob_turb_model->itytur == 5) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      prdv2f[c_id] = smbrk[c_id];
    }
  }

  /* Only for the bl-v2/k model, calculation of E and Ceps2*
   * ========================================================
   *      The terms are stored in          w10, w11
   *      Work arrays                      w2, w3
   *      viscf, viscb
   *      Going out of the step we keep w10, w11 */

  cs_real_t *coefap = NULL, *coefbp = NULL, *cofafp = NULL, *cofbfp = NULL;
  cs_real_t *w12 = NULL;

  if (cs_glob_turb_model->iturb == CS_TURB_V2F_BL_V2K) {

    /* Calculation of Ceps2*: it is stored in w10 */

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      w3[c_id] = cpro_pcvto[c_id] / cromo[c_id] / sigmak;
    }

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

    cs_real_t hint;
    /* Translate coefa into cofaf and coefb into cofbf */
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

      cs_lnum_t c_id = b_face_cells[face_id];

      hint = w3[c_id]/distb[face_id];

      /* Translate coefa into cofaf and coefb into cofbf */
      cofafp[face_id] = -hint*coefap[face_id];
      cofbfp[face_id] = hint*(1.-coefbp[face_id]);

    }

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
                           NULL,
                           cvara_k,
                           coefap,
                           coefbp,
                           cofafp,
                           cofbfp,
                           viscf,
                           viscb,
                           w3,
                           w10);

#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      w10[c_id] = -w10[c_id]/volume[c_id]/cvara_ep[c_id];
      w10[c_id] = tanh(pow(fabs(w10[c_id]), 1.5));
      if (cs_glob_turb_model->hybrid_turb == 4) {
        /* HTLES method */
        cs_real_t xcr = hybrid_fd_coeff[c_id];
        w10[c_id] = cs_turb_cpale2*(1.
                                    -  (cs_turb_cpale2-cs_turb_cpale4)
                                      / cs_turb_cpale2*w10[c_id]
                                      * cs_math_pow3(cvara_al[c_id])
                                      * (1. - xcr));
      }
      else {
        w10[c_id] = cs_turb_cpale2*(1.
                                    -  (cs_turb_cpale2-cs_turb_cpale4)
                                      / cs_turb_cpale2*w10[c_id]
                                      * cs_math_pow3(cvara_al[c_id]));
      }
    }

    /* Calculation of 2*Ceps3*(1-alpha)^3*nu*nut/eps*d2Ui/dxkdxj*d2Ui/dxkdxj:
       (i.e. E term / k)           : it is stored in w11 */

    /* Allocate a work array */
    BFT_MALLOC(w12, n_cells_ext, cs_real_t);

    _tsepls(phase_id,
            w12);

#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      cs_real_t rho  = cromo[c_id];
      cs_real_t xnu  = cpro_pcvlo[c_id]/rho;
      cs_real_t xnut = cpro_pcvto[c_id]/rho;
      cs_real_t xeps = cvara_ep[c_id];

      w11[c_id] = 2.*xnu*xnut*w12[c_id]*cs_turb_cpale3/xeps
                    *cs_math_pow3(1.-cvara_al[c_id]);
    }

    /* Add the Cazalbou rotation/curvature correction if necessary */
    if (cs_glob_turb_rans_model->irccor == 1) {
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        w10[c_id] *= ce2rc[c_id]/cs_turb_ccaze2;
        w11[c_id] *= ce2rc[c_id]/cs_turb_ccaze2;
      }
    }

    /* Free memory */
    BFT_FREE(w12);

  }

  /* Only for the Launder-Sharma model, calculation of E and D terms
   *      The terms are stored in          grad_sqk, grad_s
   * ================================================================*/

  if (cs_glob_turb_model->iturb == CS_TURB_K_EPSILON_LS) {

    /* Gradient of square root of k
     * -----------------------------*/

    coefap = (cs_real_t *)f_k->bc_coeffs->a;
    coefbp = (cs_real_t *)f_k->bc_coeffs->b;

    /* For all usual type of boundary faces (wall, inlet, sym, outlet):
       - coefa for sqrt(k) is the sqrt of the coefa for k,
       - coefb is the same as for k */

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      coefa_sqk[face_id] = sqrt(coefap[face_id]);
      coefb_sqk[face_id] = coefbp[face_id];
    }

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
                       eqp_k->imligr,
                       eqp_k->epsrgr,
                       eqp_k->climgr,
                       NULL,
                       coefa_sqk,
                       coefb_sqk,
                       sqrt_k,
                       NULL,
                       NULL,  /* internal coupling */
                       grad_sqk);

    /* Gradient of the Strain (grad S)
       --------------------------------------- */

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      coefa_sqs[face_id] = 0.;
      coefb_sqs[face_id] = 1.;
    }

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
                       eqp_k->imligr,
                       eqp_k->epsrgr,
                       eqp_k->climgr,
                       NULL,
                       coefa_sqs,
                       coefb_sqs,
                       strain,
                       NULL,
                       NULL, /* internal coupling */
                       grad_s);

  }

  /* Finalization of explicit and implicit source terms
   * ===================================================
   *
   * smbre = ceps1 epsilon/k (prod + g ) - rho0 volume epsilon epsilon/k
   * smbrk =                  prod + g   - rho0 volume epsilon */

  /* If we extrapolate the source terms and rho, we need here rho^n
     and visct, we need here visct^n */

  if (cs_glob_turb_model->itytur == 2) {

    /* Stores the production terms for the k-epsilon coupling option */
    if (cs_glob_turb_model->iturb == CS_TURB_K_EPSILON) {
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        prdtke[c_id] = smbrk[c_id];
        prdeps[c_id] = smbre[c_id];
      }
    }

#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      cs_real_t rho   = cromo[c_id];

      smbrk[c_id] =   cell_f_vol[c_id]
                    * (smbrk[c_id] - rho*cvara_ep[c_id]
                                   - d2s3*rho*cvara_k[c_id]*divu[c_id]);

      smbre[c_id] =   cell_f_vol[c_id]
                    * (  cvara_ep[c_id]/cvara_k[c_id]
                       * (  ce1rc[c_id]*smbre[c_id]
                          - ce2rc[c_id]*rho*cvara_ep[c_id])
                       - d2s3*rho*ce1rc[c_id]*cvara_ep[c_id]*divu[c_id]);

    }

    if (cs_glob_turb_model->iturb == CS_TURB_K_EPSILON_LS) {
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cs_real_t rho   = cromo[c_id];
        cs_real_t mu   = cpro_pcvlo[c_id];
        cs_real_t xnut  = cpro_pcvto[c_id]/rho;

        cs_real_t mu_gradk2 = cell_f_vol[c_id] *  2. * mu
                        * cs_math_3_square_norm(grad_sqk[c_id]);

        if (cs_glob_turb_rans_model->ikecou == 0)
          tinstk[c_id] += mu_gradk2/cvara_k[c_id];
        smbrk[c_id]  -= mu_gradk2;

        smbre[c_id] +=    cell_f_vol[c_id] * 2. * mu * xnut
                        * cs_math_3_square_norm(grad_s[c_id]);
      }
    }

    /* If the solving of k-epsilon is uncoupled,
     * negative source terms are implicited */

    if (cs_glob_turb_rans_model->ikecou == 0) {
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cs_real_t xeps = cvara_ep[c_id];
        cs_real_t xk   = cvara_k[c_id];
        cs_real_t rho = crom[c_id];
        cs_real_t ttke = xk / xeps;
        tinstk[c_id] += rho*cell_f_vol[c_id]/ttke
                     + fmax(d2s3*rho*cell_f_vol[c_id]*divu[c_id], 0.);
        tinste[c_id] += ce2rc[c_id]*rho*cell_f_vol[c_id]/ttke
                     + fmax(d2s3*ce1rc[c_id]*rho*cell_f_vol[c_id]*divu[c_id], 0.);
      }
    }

  }
  else if (cs_glob_turb_model->iturb == CS_TURB_V2F_PHI) {

#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      cs_real_t rho  = cromo[c_id];
      cs_real_t xnu  = cpro_pcvlo[c_id]/rho;
      cs_real_t xeps = cvara_ep[c_id];
      cs_real_t xk   = cvara_k[c_id];
      cs_real_t xphi = fmax(cvara_phi[c_id], cs_math_epzero);
      cs_real_t ceps1= 1.4*(1. + cs_turb_cv2fa1*sqrt(1./xphi));
      cs_real_t ttke = xk / xeps;
      cs_real_t ttmin = cs_turb_cv2fct*sqrt(xnu/xeps);
      cs_real_t tt = fmax(ttke, ttmin);

      /* Explicit part */
      smbrk[c_id] = cell_f_vol[c_id]*
                    ( smbrk[c_id] - rho*cvara_ep[c_id]
                      - d2s3*rho*cvara_k[c_id]*divu[c_id] );

      smbre[c_id] = cell_f_vol[c_id]*
                    ( 1./tt*(ceps1*smbre[c_id] - ce2rc[c_id]*rho*xeps)
                      - d2s3*rho*ceps1*xk*divu[c_id] );

      /* We store the part with Pk in prdv2f which will be reused in resv2f */
      prdv2f[c_id] = prdv2f[c_id] - d2s3*rho*cvara_k[c_id]*divu[c_id];
      /*FIXME this term should be removed */

      /* Implicit part */
      tinstk[c_id] += rho*cell_f_vol[c_id]/fmax(ttke, cs_math_epzero * ttmin);

      tinstk[c_id] += fmax(d2s3*rho*cell_f_vol[c_id]*divu[c_id], 0.);
      tinste[c_id] +=   ce2rc[c_id]*rho*cell_f_vol[c_id]/tt
                      + fmax(d2s3*ceps1*ttke/tt*rho*cell_f_vol[c_id]*divu[c_id],
                             0.);
    }

  }
  else if (cs_glob_turb_model->iturb == CS_TURB_V2F_BL_V2K) {

    if (cs_glob_turb_model->hybrid_turb == 0) {

#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

        cs_real_t rho   = cromo[c_id];
        cs_real_t xnu   = cpro_pcvlo[c_id]/rho;
        cs_real_t xeps  = cvara_ep[c_id];
        cs_real_t xk    = cvara_k[c_id];
        cs_real_t ttke  = xk / xeps;
        cs_real_t ttmin = cs_turb_cpalct*sqrt(xnu/xeps);
        cs_real_t tt    = sqrt(cs_math_pow2(ttke) + cs_math_pow2(ttmin));

        /* Explicit part */
        smbrk[c_id] = cell_f_vol[c_id] * ( smbrk[c_id]
                                          - rho*xeps
                                          - rho*w11[c_id]*xk
                                          - d2s3*rho*xk*divu[c_id]);

        smbre[c_id] =   cell_f_vol[c_id]
                      * ( 1./tt*(cs_turb_cpale1*smbre[c_id] - w10[c_id]*rho*xeps)
                         - d2s3*rho*cs_turb_cpale1*xk/tt*divu[c_id]);

        /* We store the part with Pk in prdv2f which will be reused in resv2f */
        prdv2f[c_id] = prdv2f[c_id] - d2s3*rho*cvara_k[c_id]*divu[c_id];
        /*FIXME this term should be removed */

        /* Implicit part */
        tinstk[c_id] += rho*cell_f_vol[c_id]/fmax(ttke, cs_math_epzero * ttmin);

        tinstk[c_id] += fmax(d2s3*rho*cell_f_vol[c_id]*divu[c_id], 0.);
        /* Note that w11 is positive */
        tinstk[c_id] += w11[c_id]*rho*cell_f_vol[c_id];
        /* Note that w10 is positive */
        tinste[c_id] +=   w10[c_id]*rho*cell_f_vol[c_id]/tt
                        + fmax(  d2s3*cs_turb_cpale1*ttke
                               / tt*rho*cell_f_vol[c_id]*divu[c_id],
                               0.);
      }
    }
    else if (cs_glob_turb_model->hybrid_turb == 4) {
      /* HTLES */

#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

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
        cs_real_t ttmin = cs_turb_cpalct*sqrt(xnu/(xpsi*xeps));
        cs_real_t tt    = sqrt(cs_math_pow2(ttke) + cs_math_pow2(ttmin));

        /* Explicit part */
        smbrk[c_id] = cell_f_vol[c_id] * ( smbrk[c_id]
                                          - rho*xepsh
                                          - rho*w11[c_id]*xk
                                          - d2s3*rho*xk*divu[c_id]);

        smbre[c_id] =   cell_f_vol[c_id]
                      * ( 1./tt*(cs_turb_cpale1*smbre[c_id] - w10[c_id]*rho*xeps)
                         - d2s3*rho*cs_turb_cpale1*xk/tt*divu[c_id]);

        /* We store the part with Pk in prdv2f which will be reused in resv2f */
        prdv2f[c_id] = prdv2f[c_id] - d2s3*rho*cvara_k[c_id]*divu[c_id];
        /*FIXME this term should be removed */

        /* Implicit part */
        tinstk[c_id] += rho*cell_f_vol[c_id]/fmax(xtm, cs_math_epzero * ttmin);

        tinstk[c_id] += fmax(d2s3*rho*cell_f_vol[c_id]*divu[c_id], 0.);
        /* Note that w11 is positive */
        tinstk[c_id] += w11[c_id]*rho*cell_f_vol[c_id];
        /* Note that w10 is positive */
        tinste[c_id] +=   w10[c_id]*rho*cell_f_vol[c_id]/tt
                        + fmax(  d2s3*cs_turb_cpale1*ttke
                               / tt*rho*cell_f_vol[c_id]*divu[c_id],
                               0.);
      }
    }

  }

  /* Take user source terms into account
   *
   *    The square of the scalar strain rate (strain_sq)
   *    and the trace of the velocity gradient (divu) are available.
   *    The part to be explicit is stored in       w7, w8
   *    The part to be implicit is stored in       usimpk, usimpe
   *    Going out of the step we keep              strain_sq, divu,
   * ==========================================================================*/

# pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    usimpk[c_id] = 0.;
    usimpe[c_id] = 0.;
    w7[c_id]     = 0.;
    w8[c_id]     = 0.;
  }

  cs_user_source_terms(domain,
                       f_k->id,
                       w7,
                       usimpk);

  cs_user_source_terms(domain,
                       f_eps->id,
                       w8,
                       usimpe);

  if (cs_glob_porous_model == 3) {
    cs_immersed_boundary_wall_functions(f_k->id, w7, usimpk);
    cs_immersed_boundary_wall_functions(f_eps->id, w8, usimpe);
  }

  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 0) {

    int k_interp_id = cs_field_key_id("opt_interp_id");

    /* Nudging towards optimal interpolation for k */
    if (cs_field_get_key_int(f_k, k_interp_id) >= 0)
      cs_at_data_assim_source_term(f_k->id, w7, usimpk);

    /* Nudging towards optimal interpolation for epsilon */
    if (cs_field_get_key_int(f_eps, k_interp_id) >= 0)
      cs_at_data_assim_source_term(f_k->id, w8, usimpe);

  }

  /* If source terms are extrapolated over time */
  if (istprv >= 0) {

    cs_real_t thetak = eqp_k->thetav;
    cs_real_t thetae = eqp_eps->thetav;

#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      /* Recover the value at time (n-1) */
      cs_real_t tuexpk = c_st_k_p[c_id];
      cs_real_t tuexpe = c_st_eps_p[c_id];

      /* Save the values for the next time-step */
      c_st_k_p[c_id]   = smbrk[c_id] + w7[c_id];
      c_st_eps_p[c_id] = smbre[c_id] + w8[c_id];

      /* Explicit Part */
      smbrk[c_id] = - thets*tuexpk;
      smbre[c_id] = - thets*tuexpe;
      /* It is assumed that (-usimpk > 0) and though this term is implicit */
      smbrk[c_id] += usimpk[c_id]*cvara_k[c_id];
      smbre[c_id] += usimpe[c_id]*cvara_ep[c_id];

      /* Implicit part */
      tinstk[c_id] += - usimpk[c_id]*thetak;
      tinste[c_id] += - usimpe[c_id]*thetae;

    }

    /* If no extrapolation over time */
  }
  else {

#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      /* Explicit part */
      smbrk[c_id] += usimpk[c_id]*cvara_k[c_id] + w7[c_id];
      smbre[c_id] += usimpe[c_id]*cvara_ep[c_id] + w8[c_id];

      /* Implicit part */
      tinstk[c_id] += fmax(-usimpk[c_id], 0.);
      tinste[c_id] += fmax(-usimpe[c_id], 0.);
    }

  }

  /* Account for Lagrangian 2-way coupling source terms
     -------------------------------------------------- */

  /* 2nd order not handled */

  if (cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING) {

    const cs_lagr_source_terms_t  *lag_st = cs_glob_lagr_source_terms;

    if (lag_st->ltsdyn == 1) {

      cs_real_t *lag_st_k = lag_st->st_val + (lag_st->itske-1)*n_cells_ext;
      cs_real_t *lag_st_i = lag_st->st_val + (lag_st->itsli-1)*n_cells_ext;

#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

        /* Implicit and explicit source terms on k */
        smbrk[c_id] += lag_st_k[c_id];

        /* Explicit source terms on epsilon */

        smbre[c_id] += cs_turb_ce4 * lag_st_k[c_id] * cvara_ep[c_id]
                                                    / cvara_k[c_id];

        /* Implicit source terms on k */
        tinstk[c_id] += fmax(-lag_st_i[c_id], 0.);

        /* Implicit source terms on omega */
        tinste[c_id] += fmax(-cs_turb_ce4 * lag_st_k[c_id] / cvara_k[c_id], 0.);

      }

    }

  }

  /* Mass source terms (Implicit and explicit parts)
   * Going out of the step we keep divu,  smbrk, smbre
   * ================================================= */

  if (ncesmp > 0) {

#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      w2[c_id] = 0.;
      w3[c_id] = 0.;
    }

    int ivar_k = cs_field_get_key_int(f_k, var_key_id) - 1;
    int ivar_eps = cs_field_get_key_int(f_eps, var_key_id) - 1;
    int ivar_p = cs_field_get_key_int(CS_F_(p), var_key_id) - 1;

    /* We increment smbrs with -Gamma.var_prev and rovsdt with Gamma */
    /* ivar = k; */

    cs_mass_source_terms(1,
                         1,
                         ncesmp,
                         icetsm,
                         itypsm + ncesmp*ivar_k,
                         cell_f_vol,
                         cvara_k,
                         smacel + ncesmp*ivar_k,
                         smacel + ncesmp*ivar_p,
                         smbrk,
                         w2,
                         w4);

    /* ivar = eps; */
    cs_mass_source_terms(1,
                         1,
                         ncesmp,
                         icetsm,
                         itypsm + ncesmp*ivar_eps,
                         cell_f_vol,
                         cvara_ep,
                         smacel + ncesmp*ivar_eps,
                         smacel + ncesmp*ivar_p,
                         smbre,
                         w3,
                         w5);

    /* If we extrapolate the source terms we put Gamma Pinj in c_st */
    if (istprv >= 0) {
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        c_st_k_p[c_id]   += w4[c_id];
        c_st_eps_p[c_id] += w5[c_id];
      }
    }

    /* Otherwise we put it directly in smbr */
    else {
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        smbrk[c_id] += w4[c_id];
        smbre[c_id] += w5[c_id];
      }
    }

    /* Implicit part */
#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      tinstk[c_id] += w2[c_id];
      tinste[c_id] += w3[c_id];
    }

  }

  /* Taking into account the terms of conv/diff in the second member for the
   * strengthened coupling k-epsilon (kecou == 1)
   *      Work table                       w4, w5
   *      The terms are stored in          w7 et w8, then added to smbrk, smbre
   *      Going out of the step we keep    divu,
   *      smbrk, smbre
   *      w7, w8
   * ==========================================================================*/

  if (cs_glob_turb_rans_model->ikecou == 1) {

#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      w7 [c_id] = 0.;
      w8 [c_id] = 0.;
    }

    /* Handle k */

    coefap = f_k->bc_coeffs->a;
    coefbp = f_k->bc_coeffs->b;
    cofafp = f_k->bc_coeffs->af;
    cofbfp = f_k->bc_coeffs->bf;

    if (eqp_k->idiff >=  1) {

#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        w4[c_id] = viscl[c_id] + eqp_k->idifft*cvisct[c_id]/sigmak;
      }

      cs_face_viscosity(m,
                        fvq,
                        eqp_k->imvisf,
                        w4,
                        viscf,
                        viscb);

    }
    else {

      for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
        viscf[face_id] = 0.;
      }
      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
        viscb[face_id] = 0.;
      }

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
                      coefap,
                      coefbp,
                      cofafp,
                      cofbfp,
                      imasfl,
                      bmasfl,
                      viscf,
                      viscb,
                      NULL,
                      NULL,
                      NULL,
                      NULL,
                      0, /* boundary convective flux with upwind */
                      NULL,
                      w7);

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
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        w4[c_id] = viscl[c_id]
                 + eqp_eps->idifft*cvisct[c_id]/sigmae;
      }

      cs_face_viscosity(m,
                        fvq,
                        eqp_eps->imvisf,
                        w4,
                        viscf,
                        viscb);

    }
    else {

      for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
        viscf[face_id] = 0.;
      }
      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
        viscb[face_id] = 0.;
      }

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
                      coefap,
                      coefbp,
                      cofafp,
                      cofbfp,
                      imasfl,
                      bmasfl,
                      viscf,
                      viscb,
                      NULL,
                      NULL,
                      NULL,
                      NULL,
                      0, /* boundary convective flux with upwind */
                      NULL,
                      w8);

    if (eqp_eps->verbosity >= 2) {
      cs_log_printf(CS_LOG_DEFAULT,
                    " Variable %s: EXPLICIT BALANCE =  %12.5e\n",
                    cs_field_get_label(f_eps),
                    sqrt(cs_gdot(n_cells, smbre, smbre)));
    }

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      smbrk[c_id] += w7[c_id];
      smbre[c_id] += w8[c_id];
    }

  }

  /* k-Epsilon coupling (kecou == 1)
   * ================================*/

  /* Second order is not taken into account */
  if (cs_glob_turb_rans_model->ikecou == 1) {

    if (cs_glob_turb_model->iturb == CS_TURB_K_EPSILON) {

#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

        cs_real_t rho   = crom[c_id];
        cs_real_t visct = cpro_pcvto[c_id];

        /* Coupled solving */
        cs_real_t romvsd      = 1./(rho*volume[c_id]);
        cs_real_t divp23      = d2s3*fmax(divu[c_id], 0.);
        cs_real_t epssuk      = cvara_ep[c_id]/cvara_k[c_id];

        smbrk[c_id] *= romvsd;
        smbre[c_id] *= romvsd;

        cs_real_t a11 = 1./dt[c_id] - 2.*cvara_k[c_id]/cvara_ep[c_id]
                                    *cs_turb_cmu*fmin(prdtke[c_id]/visct, 0.)
                                    + divp23;
        cs_real_t a12 = 1.;
        cs_real_t a21 = - ce1rc[c_id]*cs_turb_cmu*prdeps[c_id]/visct
                        - ce2rc[c_id]*epssuk*epssuk;
        cs_real_t a22 = 1./dt[c_id] + ce1rc[c_id]*divp23
                                    + 2.*ce2rc[c_id]*epssuk;

        cs_real_t unsdet = 1./(a11*a22 -a12*a21);

        cs_real_t deltk = ( a22*smbrk[c_id] - a12*smbre[c_id])*unsdet;
        cs_real_t delte = (-a21*smbrk[c_id] + a11*smbre[c_id])*unsdet;

        /* New source term for the iterative process */
        romvsd = rho*cell_f_vol[c_id]/dt[c_id];

        smbrk[c_id] = romvsd*deltk;
        smbre[c_id] = romvsd*delte;

      }

      /* we remove the convection/diffusion at time n from smbrk and smbre
         if they were calculated */
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        smbrk[c_id] = smbrk[c_id] - w7[c_id];
        smbre[c_id] = smbre[c_id] - w8[c_id];
      }

      /* In verini we block the
         cs_glob_turb_model->iturb != CS_TURB_K_EPSILON / kecou = 1
         combination */
    }
    else
      bft_error(__FILE__, __LINE__, 0,
              _("ikecou=1 is not validated with this turbulent model\n"
                "---------------------------------------------------"));

  }

  /* Finalization of the Right Hand Side when activating 2nd time order
   * ===================================================================*/

  if (istprv >= 0) {
    cs_real_t thetp1 = 1. + thets;
#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      smbrk[c_id] += thetp1 * c_st_k_p[c_id];
      smbre[c_id] += thetp1 * c_st_eps_p[c_id];
    }
  }

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

    if (cs_glob_turb_model->iturb == CS_TURB_V2F_BL_V2K) {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        w1[c_id] = viscl[c_id]*0.5;
    }
    else {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        w1[c_id] = viscl[c_id];
    }

    /* Variable Schmidt number */
    if (sigmak_id >= 0) {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        w1[c_id] += eqp_k->idifft*cvisct[c_id]/cpro_sigmak[c_id];
    }
    else {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        w1[c_id] += eqp_k->idifft*cvisct[c_id]/sigmak;
    }

    cs_face_viscosity(m,
                      fvq,
                      eqp_k->imvisf,
                      w1,
                      viscf,
                      viscb);

  }
  else {

    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
      viscf[face_id] = 0.;
    }
    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      viscb[face_id] = 0.;
    }

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
                                     coefap,
                                     coefbp,
                                     cofafp,
                                     cofbfp,
                                     imasfl,
                                     bmasfl,
                                     viscf,
                                     viscb,
                                     viscf,
                                     viscb,
                                     NULL,
                                     NULL,
                                     NULL,
                                     0, /* boundary convective upwind flux */
                                     NULL,
                                     tinstk,
                                     smbrk,
                                     cvar_k,
                                     dpvar,
                                     NULL,
                                     NULL);

  /* Solve for epsilon */

  coefap = f_eps->bc_coeffs->a;
  coefbp = f_eps->bc_coeffs->b;
  cofafp = f_eps->bc_coeffs->af;
  cofbfp = f_eps->bc_coeffs->bf;

  /* Face viscosity */
  if (eqp_eps->idiff >= 1) {
    if (cs_glob_turb_model->iturb == CS_TURB_V2F_BL_V2K) {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        w1[c_id] = viscl[c_id]*0.5;
    } else {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        w1[c_id] = viscl[c_id];
    }

    /* Variable Schmidt number */
    if (sigmae_id >= 0) {
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        w1[c_id] += eqp_eps->idifft*cvisct[c_id]/cpro_sigmae[c_id];
    }
    else {
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        w1[c_id] += eqp_eps->idifft*cvisct[c_id]/sigmae;
    }

    cs_face_viscosity(m,
                      fvq,
                      eqp_eps->imvisf,
                      w1,
                      viscf,
                      viscb);
  }
  else {

    cs_array_real_fill_zero(n_i_faces, viscf);
    cs_array_real_fill_zero(n_b_faces, viscb);

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
                                     coefap,
                                     coefbp,
                                     cofafp,
                                     cofbfp,
                                     imasfl,
                                     bmasfl,
                                     viscf,
                                     viscb,
                                     viscf,
                                     viscb,
                                     NULL,
                                     NULL,
                                     NULL,
                                     0, /* boundary convective upwind flux */
                                     NULL,
                                     tinste,
                                     smbre,
                                     cvar_ep,
                                     dpvar,
                                     NULL,
                                     NULL);

  /* Clip values
     ============ */

  cs_turbulence_ke_clip(phase_id,
                        n_cells,
                        1);

  /* Free memory */
  if (f_vel->grad == NULL)
    BFT_FREE(gradv);
  BFT_FREE(viscf);
  BFT_FREE(viscb);
  BFT_FREE(usimpk);
  BFT_FREE(smbrk);
  BFT_FREE(smbre);
  BFT_FREE(tinstk);
  BFT_FREE(tinste);
  BFT_FREE(divu);
  BFT_FREE(strain_sq);
  BFT_FREE(w1);
  BFT_FREE(w2);
  BFT_FREE(w3);
  BFT_FREE(w4);
  BFT_FREE(w5);
  BFT_FREE(w7);
  BFT_FREE(w8);
  BFT_FREE(usimpe);
  BFT_FREE(dpvar);
  BFT_FREE(ce2rc);
  BFT_FREE(ce1rc);

  if (cs_glob_turb_model->iturb == CS_TURB_V2F_BL_V2K){
    BFT_FREE(w10);
    BFT_FREE(w11);
  }
  if (cs_glob_turb_model->iturb == CS_TURB_K_EPSILON){
    BFT_FREE(prdtke);
    BFT_FREE(prdeps);
  }
  if (cs_glob_turb_model->iturb == CS_TURB_K_EPSILON_LS) {
    BFT_FREE(sqrt_k);
    BFT_FREE(strain);
    BFT_FREE(grad_sqk);
    BFT_FREE(grad_s);
    BFT_FREE(coefa_sqk);
    BFT_FREE(coefb_sqk);
    BFT_FREE(coefa_sqs);
    BFT_FREE(coefb_sqs);
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

  const cs_equation_param_t *eqp
    = cs_field_get_equation_param_const(f_k);

  int iwarnk = eqp->verbosity;

  /* Small value to avoid exactly zero values */

  const double epz2 = cs_math_pow2(cs_math_epzero);

  /* Postprocess clippings? */

  int key_clipping_id = cs_field_key_id("clipping_id");

  int clip_k_id = cs_field_get_key_int(f_k, key_clipping_id);
  cs_real_t *cpro_k_clipped = NULL;
  if (clip_k_id >= 0) {
    cpro_k_clipped = cs_field_by_id(clip_k_id)->val;
  }

  int clip_e_id = cs_field_get_key_int(f_eps, key_clipping_id);
  cs_real_t *cpro_e_clipped = NULL;
  if (clip_e_id >= 0) {
    cpro_e_clipped = cs_field_by_id(clip_e_id)->val;
  }

  /* Save min and max for log
   * ======================== */

  const cs_real_t l_threshold = 1.e12;
  cs_real_t *cvar_var = NULL;
  cs_real_t var;
  cs_real_t vmax[2], vmin[2];

  for (int ii = 0; ii < 2; ii++ ) {
    if (ii == 0)
      cvar_var = cvar_k;
    else if (ii == 1)
      cvar_var = cvar_ep;

    vmin[ii] =  l_threshold;
    vmax[ii] = -l_threshold;

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      var = cvar_var[c_id];
      vmin[ii] = CS_MIN(vmin[ii], var);
      vmax[ii] = CS_MAX(vmax[ii], var);
    }
  }

  if (cpro_k_clipped != NULL) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      cpro_k_clipped[c_id] = 0.;
  }
  if (cpro_e_clipped != NULL) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      cpro_e_clipped[c_id] = 0.;
  }

  /* Detect values ouside "physical" bounds,
   * only for warning or when ICLKEP = 1
   * ===================================== */

  cs_gnum_t iclpke = 0;
  cs_lnum_t iclpmn[2] = {0, 0};
  cs_real_t xk, xe, xkmin, xepmin, xkm, xepm;

  if (iwarnk >= 2 || cs_glob_turb_rans_model->iclkep == 1) {

    if (iclip == 1) {

      xkm = 1296.*sqrt(cs_turb_cmu)/cs_math_pow2(almax);
      xepm = 46656.*cs_turb_cmu/cs_math_pow4(almax);
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        xk = cvar_k[c_id];
        xe = cvar_ep[c_id];
        xkmin = xkm * cs_math_pow2(viscl[c_id] / crom[c_id]);
        xepmin = xepm * cs_math_pow3(viscl[c_id] / crom[c_id]);
        if (xk <= xkmin || xe <= xepmin) {
          if (cs_glob_turb_rans_model->iclkep == 1) {
            if (clip_k_id >= 0)
              cpro_k_clipped[c_id] = xkmin - xk;
            cvar_k[c_id]  = xkmin;
            if (clip_e_id >= 0)
              cpro_e_clipped[c_id] = xepmin - xe;
            cvar_ep[c_id] = xepmin;
          }
          iclpke += 1;
        }
      }

    }
    else if (iclip == 0) {

      xkmin = 1296. * sqrt(cs_turb_cmu) / cs_math_pow2(almax)
                    * cs_math_pow2(viscl0/ro0);
      xepmin = 46656. * cs_turb_cmu/cs_math_pow4(almax)
                      * cs_math_pow3(viscl0/ro0);
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++){
        xk = cvar_k[c_id];
        xe = cvar_ep[c_id];
        if (xk <= xkmin || xe <= xepmin) {
          if (cs_glob_turb_rans_model->iclkep == 1) {
            cvar_k[c_id]  = xkmin;
            if (clip_k_id >= 0)
              cpro_k_clipped[c_id] = xkmin - xk;
            cvar_ep[c_id] = xepmin;
            if (clip_e_id >= 0)
              cpro_e_clipped[c_id] = xepmin - xe;
          }
          iclpke += 1;
        }
      }

    }
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Call of %s with option = %d"),
                __func__, iclip);

    /* save clip counts for log */

    if (cs_glob_turb_rans_model->iclkep == 1) {
      iclpmn[0] = iclpke;
      iclpmn[1] = iclpke;
    }

    /* logging */

    if (iwarnk >= 2) {
      cs_parall_sum(1, CS_GNUM_TYPE, &iclpke);

      cs_log_printf(CS_LOG_DEFAULT,
                    "\n "
                    "%llu k-epsilon values beyond the scales based on almax\n",
                    (unsigned long long)iclpke);
    }

  }

  /* "standard" clipping ICLKEP = 0
   * ============================== */

  cs_lnum_t iclpk2 = 0, iclpe2 = 0;

  if (cs_glob_turb_rans_model->iclkep == 0) {

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      xk = cvar_k[c_id];
      xe = cvar_ep[c_id];
      if (fabs(xk) <= epz2) {
        iclpk2 = iclpk2 + 1;
        if (clip_k_id >= 0)
          cpro_k_clipped[c_id] = epz2 - cvar_k[c_id];
        cvar_k[c_id] = CS_MAX(cvar_k[c_id],epz2);
      }
      else if(xk <= 0.) {
        iclpk2 = iclpk2 + 1;
        if (clip_k_id >= 0)
          cpro_k_clipped[c_id] = -xk;
        cvar_k[c_id] = -xk;
      }
      if (fabs(xe) <= epz2) {
        iclpe2 = iclpe2 + 1;
        if (clip_e_id >= 0)
          cpro_e_clipped[c_id] = epz2 - cvar_ep[c_id];
        cvar_ep[c_id] = CS_MAX(cvar_ep[c_id], epz2);
      }
      else if(xe <= 0.) {
        iclpe2 = iclpe2 + 1;
        if (clip_e_id >= 0)
          cpro_e_clipped[c_id] = - xe;
        cvar_ep[c_id] = - xe;
      }
    }

    /* save clip counts for log */

    iclpmn[0] = iclpk2;
    iclpmn[1] = iclpe2;
  }

  cs_lnum_t iclpmx[1] = {0};
  int id;

  for (int ii = 0; ii < 2; ii++ ) {
    if (ii == 0)
      id = f_k->id;
    else if (ii == 1)
      id = f_eps->id;

    cs_log_iteration_clipping_field(id,
                                    iclpmn[ii],
                                    0,
                                    vmin + ii,
                                    vmax + ii,
                                    iclpmn + ii,
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

  cs_real_t *visct = f_mut->val;
  const cs_real_t *viscl =  f_mu->val;
  const cs_real_t *crom = f_rho->val;
  const cs_real_t *cvar_k = f_k->val;
  const cs_real_t *cvar_ep = f_eps->val;

  /* Launder-Sharma */

  if (cs_glob_turb_model->iturb == CS_TURB_K_EPSILON_LS) {

    const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
    const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

    /* Initialization
     * ============== */

#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      const cs_real_t xk = cvar_k[c_id];
      const cs_real_t xe = cvar_ep[c_id];
      const cs_real_t xmu = viscl[c_id];
      const cs_real_t xmut = crom[c_id] * cs_math_pow2(xk) / xe;
      const cs_real_t xfmu = exp(-3.4 / cs_math_pow2(1. + xmut/xmu/50.));

      visct[c_id] = cs_turb_cmu * xfmu * xmut;
    }

  }

  /* Baglietto model */

  else if (cs_glob_turb_model->iturb == CS_TURB_K_EPSILON_QUAD) {
    cs_turbulence_ke_q_mu_t(phase_id);
  }

  /* Standard and linear-production */

  else {

    const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

#   pragma omp parallel for if(n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      visct[c_id] =   crom[c_id] * cs_turb_cmu
                    * cs_math_pow2(cvar_k[c_id]) / cvar_ep[c_id];
    }

  }
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

  cs_real_t *s2;
  BFT_MALLOC(s2, n_cells_ext, cs_real_t);

  /* Calculation of velocity gradient and of
   *   S2 = S11**2+S22**2+S33**2+2*(S12**2+S13**2+S23**2)
   * ==================================================== */

  cs_real_33_t *gradv = NULL;

  if (f_vel->grad == NULL) {
    BFT_MALLOC(gradv, n_cells_ext, cs_real_33_t);

    /* Computation of the velocity gradient */

    cs_field_gradient_vector(f_vel,
                             false,     // no use_previous_t
                             1,         // inc
                             gradv);
  } else
    gradv = (cs_real_33_t *)f_vel->grad;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    const cs_real_t s11 = gradv[c_id][0][0];
    const cs_real_t s22 = gradv[c_id][1][1];
    const cs_real_t s33 = gradv[c_id][2][2];
    const cs_real_t dudy = gradv[c_id][0][1];
    const cs_real_t dudz = gradv[c_id][0][2];
    const cs_real_t dvdx = gradv[c_id][1][0];
    const cs_real_t dvdz = gradv[c_id][1][2];
    const cs_real_t dwdx = gradv[c_id][2][0];
    const cs_real_t dwdy = gradv[c_id][2][1];

    s2[c_id] =   cs_math_pow2(s11) + cs_math_pow2(s22) + cs_math_pow2(s33)
               + 0.5*cs_math_pow2(dudy+dvdx)
               + 0.5*cs_math_pow2(dudz+dwdx)
               + 0.5*cs_math_pow2(dvdz+dwdy);
  }

  if (f_vel->grad == NULL)
    BFT_FREE(gradv);

  /* Compute turbulent viscosity
   * =========================== */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    const cs_real_t xk = cvar_k[c_id];
    const cs_real_t xe = cvar_ep[c_id];
    const cs_real_t xrom = crom[c_id];
    const cs_real_t xmu = viscl[c_id];
    const cs_real_t xdist = fmax(w_dist[c_id], 1.e-10);

    const cs_real_t xmut = xrom*cs_math_pow2(xk)/xe;
    const cs_real_t xrey = xdist*sqrt(xk)*xrom/xmu;
    const cs_real_t xttke = xk/xe;
    const cs_real_t xss = xttke*sqrt(2.0*s2[c_id]);

    const cs_real_t xfmu = 1.0 - exp(- 2.9e-2*sqrt(xrey)
                                     - 1.1e-4*cs_math_pow2(xrey));
    const cs_real_t xcmu = 2.0 / 3.0 / (3.90 + xss);

    visct[c_id] = xcmu*xfmu*xmut;
  }

  BFT_FREE(s2);
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

  cs_real_33_t *gradv = NULL;

  if (f_vel->grad == NULL) {
    BFT_MALLOC(gradv, n_cells_ext, cs_real_33_t);

    /* Computation of the velocity gradient */

    cs_field_gradient_vector(f_vel,
                             true, // use_previous_t
                             1,    // inc
                             gradv);
  } else
    gradv = (cs_real_33_t *)f_vel->grad;

  const cs_real_t d1s2 = 0.5, d2s3 = 2./3;

  /*  Computation
   *===============*/

# pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {

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
      = cs_turb_cnl1/((cs_turb_cnl4+cs_turb_cnl5*cs_math_pow3(xss))*xcmu);
    const cs_real_t xqc2
      = cs_turb_cnl2/((cs_turb_cnl4+cs_turb_cnl5*cs_math_pow3(xss))*xcmu);
    const cs_real_t xqc3
      = cs_turb_cnl3/((cs_turb_cnl4+cs_turb_cnl5*cs_math_pow3(xss))*xcmu);

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
  }

  /* Free memory */
  if (f_vel->grad == NULL)
    BFT_FREE(gradv);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
