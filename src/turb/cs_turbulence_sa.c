/*============================================================================
 * nu_tilda turbulence model.
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
#include "cs_domain.h"
#include "cs_equation.h"
#include "cs_equation_iterative_solve.h"
#include "cs_face_viscosity.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
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
#include "cs_time_step.h"
#include "cs_turbulence_model.h"
#include "cs_turbulence_rotation.h"
#include "cs_wall_functions.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_turbulence_sa.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_turbulence_sa.c

  Solving op the equation of \f$ \tilde{\nu} \f$, which is the scalar
  quantity defined by the Spalart-Allmaras model for 1 time-step.
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
 * \brief Compute the vorticity omega, the trace of the velocity gradient
 *        and the gradient of nusa.
 *
 *  \param[out]  vor_ome   vorticity omega
 *  \param[out]  tr_gr_u   trace of the velocity gradient
 *  \param[out]  tr_gr_nu  trace of the gradient of nusa
 */
/*----------------------------------------------------------------------------*/

static void
_vort_trace(cs_real_t   vort[],
            cs_real_t   tr_gr_u[],
            cs_real_t   tr_gr_nu[])
{
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  /* Allocate a temporary for the gradient calculation */

  cs_real_33_t *grad_vel;
  BFT_MALLOC(grad_vel, n_cells_ext, cs_real_33_t);

  cs_field_gradient_vector(CS_F_(vel),
                           true, // use_previous_t
                           1,    // inc
                           grad_vel);

  for (cs_lnum_t i = 0; i < n_cells; i++) {

    vort[i] =   cs_math_pow2(grad_vel[i][0][1] - grad_vel[i][1][0])
              + cs_math_pow2(grad_vel[i][0][2] - grad_vel[i][2][0])
              + cs_math_pow2(grad_vel[i][1][2] - grad_vel[i][2][1]);

    tr_gr_u[i] = grad_vel[i][0][0] + grad_vel[i][1][1] + grad_vel[i][2][2];

  }

  BFT_FREE(grad_vel);

  cs_real_3_t *grad_nu;
  BFT_MALLOC(grad_nu, n_cells_ext, cs_real_3_t);

  cs_field_gradient_scalar(CS_F_(nusa),
                           true,   /* use_previous_t */
                           1,      /* inc */
                           grad_nu);

  for (cs_lnum_t i = 0; i < n_cells; i++) {
    tr_gr_nu[i] = cs_math_3_square_norm(grad_nu[i]);
  }

  BFT_FREE(grad_nu);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the Source terms
 *
 *  \param[in]      dt           time step (per cell)
 *  \param[in]      b_face_type  boundary face type
 *  \param[in]      tr_gr_nu     trace of the gradient of field nusa
 *  \param[in]      vort         vorticity
 *  \param[in]      cpro_rho_o   density (at current or previous time step)
 *  \param[in]      cpro_viscl   laminar viscosity
 *  \param[in,out]  st_exp       explicit source terms
 *  \param[in,out]  st_imp       implicit user source terms
 */
/*----------------------------------------------------------------------------*/

static void
_src_terms(const cs_real_t    dt[],
           const int          b_face_type[],
           const cs_real_t    tr_gr_nu[],
           const cs_real_t    vort[],
           const cs_real_t    cpro_rho_o[],
           const cs_real_t    cpro_viscl[],
           cs_real_t          st_exp[],
           cs_real_t          st_imp[])
{
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  const cs_real_t *b_dist = cs_glob_mesh_quantities->b_dist;
  const cs_real_t *cell_f_vol = cs_glob_mesh_quantities->cell_f_vol;

  const cs_real_t *w_dist = cs_field_by_name("wall_distance")->val;
  const cs_real_t *cvara_nusa = CS_F_(nusa)->val_pre;

  const cs_real_t dsigma = 1.0 / cs_turb_csasig;

  cs_real_t dsa0 = -999.0;
  cs_real_t hssa = -999.0;

  /* Here, we only handle the case where all the walls have uniform roughness
     (we take the average as a precaution).
     To extend it, we should be able to associate every fluid cell to a boundary
     face (and then give it the appropriate roughness value). This could be done
     at the cost of using a diffusion equation. */

  const cs_field_t  *f_r = cs_field_by_name_try("boundary_roughness");
  if (f_r != NULL) {
    const cs_real_t *b_roughness = f_r->val;
    const cs_real_t *coefbp = CS_F_(nusa)->bc_coeffs->b;

    cs_real_t s[2] = {0, 0};
    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      if (b_face_type[f_id] == CS_SMOOTHWALL && b_roughness[f_id] > 0) {
        const cs_real_t cofbnu = coefbp[f_id];
        /* Roughness of the wall */
        s[0] += b_dist[f_id] * cofbnu/(1.0 - cofbnu);  /* dsa0 */
        s[1] += 1;
      }
    }

    cs_parall_sum(2, CS_REAL_TYPE, s);
    if (s[1] > 0) {
      dsa0 = s[0]/s[1];
      hssa = exp(8.50*cs_turb_xkappa)*dsa0;
    }
  }

  /* Take into account the Spalart-Shur rotation/curvature correction,
     if necessary
     => variable production term coefficient (csab1) */
  cs_real_t *csab1r;
  BFT_MALLOC(csab1r, n_cells, cs_real_t);

  if (cs_glob_turb_rans_model->irccor == 1) {
    cs_real_t *w1;
    BFT_MALLOC(w1, n_cells_ext, cs_real_t);

    /* Compute the rotation function (w1 array not used) */
    cs_turbulence_rotation_correction(dt, csab1r, w1);

    for (cs_lnum_t i = 0; i < n_cells; i++) {
      csab1r[i] *= cs_turb_csab1;
    }

    BFT_FREE(w1);
  }
  else {
    cs_array_real_set_scalar(n_cells, cs_turb_csab1, csab1r);
  }

  /* To avoid numerical problem, constant used to prevent taussa from
     being negative (see Oliver TA 2008) */
  const cs_real_t cst2 = 0.7;
  const cs_real_t cst3 = 0.9;

  /* If source terms are extrapolated, rho is rho^n
     visct is visct^n (visct not used here) */

  const cs_real_t cv13 = cs_math_pow3(cs_turb_csav1);

  for (cs_lnum_t i = 0; i < n_cells; i++) {

    const cs_real_t rho = cpro_rho_o[i];

    /* kinematic viscosity */
    const cs_real_t nu0 = cpro_viscl[i]/rho;

    /* We have to know if there is any rough wall */
    cs_real_t distbf = w_dist[i];

    /* viscosity of SA */
    const cs_real_t nusa = cvara_nusa[i];

    cs_real_t chi = nusa/nu0;
    /* If we have a rough wall */
    if (dsa0 > -998) {
      distbf += dsa0;
      chi += 0.50* hssa/distbf;
    }
    const cs_real_t chi3 = cs_math_pow3(chi);
    const cs_real_t fv1 = chi3/(chi3 + cv13);
    const cs_real_t fv2 = 1.0 - nusa /(nu0 + nusa*fv1);

    /* Numerical fix to prevent taussa to be smaller than 0 */
    const cs_real_t sbar = nusa/cs_math_pow2(cs_turb_xkappa*distbf)*fv2;
    const cs_real_t omega = sqrt(vort[i]);
    cs_real_t taussa;

    if (sbar >= -cst2*omega)
      taussa = omega + sbar;
    else
      taussa = omega * (1.0 + (cs_math_pow2(cst2)*omega + cst3*sbar)
                               / ((cst3 - 2.*cst2) * omega-sbar));

    /* Compute fw */
    cs_real_t rsa;
    if (nusa >= 10.0*taussa*cs_math_pow2(cs_turb_xkappa*distbf))
      rsa = 10.0;
    else
      rsa = nusa / (taussa*cs_math_pow2(cs_turb_xkappa*distbf));

    const cs_real_t rsa_6 = rsa*rsa*rsa * rsa*rsa*rsa;
    const cs_real_t gsa = rsa + cs_turb_csaw2*(rsa_6 - rsa);
    const cs_real_t gsa_6 = gsa*gsa*gsa * gsa*gsa*gsa;

    const cs_real_t csaw3_6 = cs_math_pow2(cs_math_pow3(cs_turb_csaw3));

    const cs_real_t fw = gsa * pow((1.0 + csaw3_6) / (gsa_6 + csaw3_6),
                                   1.0/6.0);

    st_exp[i] =   cell_f_vol[i] * rho
                * (  dsigma * cs_turb_csab2 * tr_gr_nu[i]
                   + csab1r[i] * taussa * nusa
                   - cs_turb_csaw1 * fw * cs_math_pow2(nusa/distbf));

    /* Implicitation of the negative source terms of the SA equation.
       NB: this term may be negative, and if so, then we explicit it. */

    st_imp[i] = (cs_math_fmax((  cs_turb_csaw1*fw*nusa/cs_math_pow2(distbf)
                               - csab1r[i]*taussa),
                              0.0))
                * rho * cell_f_vol[i];

  }

  BFT_FREE(csab1r);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Clipping of nusa for the Spalart-Allmaras model.
 *
 *  \param[in]  n_cells  number of cells
 */
/*----------------------------------------------------------------------------*/

static void
_clip(cs_lnum_t  n_cells)
{
  cs_real_t *cvar_nusa = CS_F_(nusa)->val;

  int key_clipping_id = cs_field_key_id("clipping_id");

  /* Postprocess clippings? */
  int clip_nusa_id = cs_field_get_key_int(CS_F_(nusa), key_clipping_id);
  cs_real_t *cpro_nusa_clipped = NULL;
  if (clip_nusa_id >= 0) {
    cpro_nusa_clipped = cs_field_by_id(clip_nusa_id)->val;
    cs_array_real_fill_zero(n_cells, cpro_nusa_clipped);
  }

  /* Save min and max for log */

  cs_lnum_t iclpmx = 0, iclpmn = 0;
  cs_real_t xnu_min = cs_math_big_r, xnu_max = -cs_math_big_r;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_real_t xnu = cvar_nusa[c_id];
    xnu_min = CS_MIN(xnu_min, xnu);
    xnu_max = CS_MAX(xnu_max, xnu);
  }

  /* "Standard" clipping  NUSA > 0 */

  cs_lnum_t iclpnu = 0;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_real_t xnu = cvar_nusa[c_id];
    if (xnu < 0) {
      if (clip_nusa_id > -1)
        cpro_nusa_clipped[c_id] = - xnu;
      iclpnu += 1;
      cvar_nusa[c_id] = 0;
    }
  }

  cs_log_iteration_clipping_field(CS_F_(nusa)->id,
                                  iclpnu,
                                  0,
                                  &xnu_min,
                                  &xnu_max,
                                  &iclpmn,
                                  &iclpmx);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve the \f$ \tilde{\nu} \f$ equations.
 *
 * Solve the equation of \f$ \tilde{\nu} \f$, which is the scalar
 * quantity defined by the Spalart-Allmaras model for one time-step.
 *
 * \param[in]     ncesmp        number of cells with mass source term
 * \param[in]     icetsm        index of cells with mass source term
 * \param[in]     itypsm        mass source type for the variables
 *                              size: [nvar][ncesmp]
 * \param[in]     dt            time step (per cell)
 * \param[in]     smacel        values of the variables associated to the
 *                              mass source (for the pressure variable,
 *                              smacel is the mass flux)
 *                              size: [nvar][ncesmp]
 * \param[in]  b_face_type      boundary face types
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_sa(cs_lnum_t        ncesmp,
                 cs_lnum_t        icetsm[],
                 int              itypsm[],
                 const cs_real_t  dt[],
                 cs_real_t        smacel[],
                 const int        b_face_type[])
{
  cs_domain_t  *domain = cs_glob_domain;

  const cs_mesh_t *m = domain->mesh;
  const cs_mesh_quantities_t  *fvq = domain->mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_i_faces = m->n_i_faces;

  const cs_real_t *cell_f_vol = fvq->cell_f_vol;

  const cs_equation_param_t *eqp_nusa
    = cs_field_get_equation_param_const(CS_F_(nusa));

  const cs_real_t *cpro_rho = CS_F_(rho)->val;
  const cs_real_t *cpro_rho_o = CS_F_(rho)->val;
  const cs_real_t *cpro_viscl = CS_F_(mu)->val;

  int key_t_ext_id = cs_field_key_id("time_extrapolated");
  int kstprv = cs_field_key_id("source_term_prev_id");

  cs_real_t *c_st_nusa_p = NULL;

  int istprv = cs_field_get_key_int(CS_F_(nusa), kstprv);
  if (istprv >= 0) {
    c_st_nusa_p = cs_field_by_id(istprv)->val;
    if (cs_field_get_key_int(CS_F_(rho), key_t_ext_id) > 0)
      cpro_rho_o = CS_F_(rho)->val_pre;
  }

  if (eqp_nusa->verbosity >= 1)
    cs_log_printf(CS_LOG_DEFAULT,
                  "\n"
                  "  ** Solving Spalart-Allmaras\n"
                  "     ------------------------\n");

  cs_real_t *cvar_nusa = CS_F_(nusa)->val;
  cs_real_t *cvara_nusa = CS_F_(nusa)->val_pre;

  cs_real_t *vort;      /* vorticity omega*/
  cs_real_t *tr_gr_u;   /* trace of the velocity gradient */
  cs_real_t *tr_gr_nu;  /* trace of the gradient of field nusa */

  cs_real_t *imp_sa, *rhs_sa;

  BFT_MALLOC(vort, n_cells_ext, cs_real_t);
  BFT_MALLOC(tr_gr_u, n_cells_ext, cs_real_t);
  BFT_MALLOC(tr_gr_nu, n_cells_ext, cs_real_t);
  BFT_MALLOC(rhs_sa, n_cells_ext, cs_real_t);
  BFT_MALLOC(imp_sa, n_cells_ext, cs_real_t);

  /* Compute the vorticity omega, the trace of the velocity gradient
     and the gradient of nusa */

  _vort_trace(vort, tr_gr_u, tr_gr_nu);

  /* Compute the buoyant term:
     gravity is not taken into account at the moment */

  /* Source terms are finalized, stored in st_exp */

  _src_terms(dt,
             b_face_type,
             tr_gr_nu,
             vort,
             cpro_rho_o,
             cpro_viscl,
             rhs_sa,
             imp_sa);

  BFT_FREE(vort);
  BFT_FREE(tr_gr_u);
  BFT_FREE(tr_gr_nu);

  /* Take user source terms into account */
  /*!
   * omega² = vort and the trace of the velocity gradient = tr_gr_u
   * are available
   * The explicit part is stored in st_exp
   * The implicit part is stored in st_imp
   !*/
  cs_real_t *st_exp, *st_imp;

  BFT_MALLOC(st_imp, n_cells_ext, cs_real_t);
  BFT_MALLOC(st_exp, n_cells_ext, cs_real_t);

  for (cs_lnum_t i = 0; i < n_cells; i++) {
    st_exp[i] = 0;
    st_imp[i] = 0;
  }

  cs_user_source_terms(domain,
                       CS_F_(nusa)->id,
                       st_exp,
                       st_imp);

  if (cs_glob_porous_model == 3)
    cs_immersed_boundary_wall_functions(CS_F_(nusa)->id,
                                        st_exp,
                                        st_imp);

  /* User source terms and d/dt(rho) and div(rho u) are taken into account
     stored in ext_term */

  /* If source terms are extrapolated */
  if (istprv >= 0) {
    const cs_time_scheme_t *time_scheme = cs_glob_time_scheme;
    const cs_real_t thetst = time_scheme->thetst;
    const cs_real_t thetv  = eqp_nusa->thetav;

    for (cs_lnum_t i = 0; i < n_cells; i++) {
      const cs_real_t tuexpn = c_st_nusa_p[i];
      c_st_nusa_p[i] = rhs_sa[i] + st_exp[i];

      /* Extrapolated explicit source terms */
      rhs_sa[i] = cvara_nusa[i]*st_imp[i] - thetst*tuexpn;

      /* Implicit user source terms*/
      /* Here it is assumed that -tsimp > 0. That is why it is implicited */
      imp_sa[i] -= st_imp[i]*thetv;
    }

  }
  else {
    for (cs_lnum_t i = 0; i < n_cells; i++) {
      rhs_sa[i] += cvara_nusa[i]*st_imp[i] + st_exp[i];
      imp_sa[i] += cs_math_fmax(st_imp[i], cs_math_epzero);
    }
  }

  BFT_FREE(st_exp);
  BFT_FREE(st_imp);

  for (cs_lnum_t i = 0; i < n_cells; i++) {
    const cs_real_t romvsd = cpro_rho[i]*cell_f_vol[i] / dt[i];
    /* imp_sa already contains the negative implicited source term */
    imp_sa[i] += eqp_nusa->istat*romvsd;
  }

  /* Explicit mass source terms */

  /* Gamma.var_prev is stored in w_1 */
  cs_real_t *w_1;
  BFT_MALLOC(w_1, n_cells_ext, cs_real_t);

  if (ncesmp > 0) {
    const int var_key_id = cs_field_key_id("variable_id");
    cs_lnum_t ivar = cs_field_get_key_int(CS_F_(nusa), var_key_id) - 1;
    cs_lnum_t ipr = cs_field_get_key_int(CS_F_(p), var_key_id) - 1;

    cs_mass_source_terms(1,
                         1,
                         ncesmp,
                         icetsm,
                         itypsm + ncesmp*ivar,
                         cell_f_vol,
                         cvara_nusa,
                         smacel + ncesmp*ivar,
                         smacel + ncesmp*ipr,
                         rhs_sa,
                         imp_sa,
                         w_1);

    /* Explicit part: Gamma Pinj
       (if we extrapolate source terms, Gamma.var_prev is stored in prev. TS) */
    if (istprv >= 0) {
      for (cs_lnum_t i = 0; i < n_cells; i++) {
        c_st_nusa_p[i] +=  w_1[i];
      }
    }
    else {
      for (cs_lnum_t i = 0; i < n_cells; i++) {
        rhs_sa[i] += w_1[i];
      }
    }

  }

  /* Finalization of the extrapolated explicit source terms */

  if (istprv >= 0) {
    const cs_time_scheme_t *time_scheme = cs_glob_time_scheme;
    const cs_real_t thetst  = time_scheme->thetst;
    const cs_real_t thetp1 = 1.0 + thetst;

    for (cs_lnum_t i = 0; i < n_cells; i++) {
      rhs_sa[i] += thetp1 * c_st_nusa_p[i];
    }
  }

  /* Solving of the transport equation on nusa */

  cs_real_t *viscf, *viscb;
  BFT_MALLOC(viscf, n_i_faces, cs_real_t);
  BFT_MALLOC(viscb, n_b_faces, cs_real_t);

  if (eqp_nusa->idiff >= 1) {
    const cs_real_t dsigma = 1.0 / cs_turb_csasig;
    const int idifft = eqp_nusa->idifft;

    /* diffusibility: 1/sigma*(mu_laminar + rho*nusa) */
    for (cs_lnum_t i = 0; i < n_cells; i++) {
      w_1[i] = dsigma *(cpro_viscl[i] + idifft*cvara_nusa[i]*cpro_rho[i]);
    }

    cs_face_viscosity(m,
                      fvq,
                      cs_glob_space_disc->imvisf,
                      w_1,
                      viscf,
                      viscb);
  }
  else {
    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++)
      viscf[f_id] = 0.0;
    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++)
      viscb[f_id] = 0.0;
  }

  BFT_FREE(w_1);

  /* Solving */

  const cs_real_t *coefap = CS_F_(nusa)->bc_coeffs->a;
  const cs_real_t *coefbp = CS_F_(nusa)->bc_coeffs->b;
  const cs_real_t *cofafp = CS_F_(nusa)->bc_coeffs->af;
  const cs_real_t *cofbfp = CS_F_(nusa)->bc_coeffs->bf;

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  const int iflmas = cs_field_get_key_int(CS_F_(nusa), kimasf);
  const int iflmab = cs_field_get_key_int(CS_F_(nusa), kbmasf);
  const cs_real_t *imasfl = cs_field_by_id(iflmas)->val;
  const cs_real_t *bmasfl = cs_field_by_id(iflmab)->val;

  cs_real_t *dpvar;
  BFT_MALLOC(dpvar, n_cells_ext, cs_real_t);

  cs_equation_param_t _eqp_nusa = *eqp_nusa;

  cs_equation_iterative_solve_scalar(cs_glob_time_step_options->idtvar,
                                     1,    /* init */
                                     CS_F_(nusa)->id,
                                     CS_F_(nusa)->name,
                                     0,   /* iescap */
                                     0,   /* imucpp */
                                     -1,  /* normp */
                                     &_eqp_nusa,
                                     cvara_nusa,
                                     cvara_nusa,
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
                                     imp_sa,
                                     rhs_sa,
                                     cvar_nusa,
                                     dpvar,
                                     NULL,
                                     NULL);

  /* Clip values */
  _clip(n_cells);

  BFT_FREE(rhs_sa);
  BFT_FREE(imp_sa);
  BFT_FREE(dpvar);
  BFT_FREE(viscf);
  BFT_FREE(viscb);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of turbulent viscosity for
 *        the Spalart-Allmaras model.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_sa_mu_t(void)
{
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_lnum_t n_cells = mesh->n_cells;

  const cs_real_t cv13 = cs_math_pow3(cs_turb_csav1);

  cs_field_t *f_nusa = CS_F_(nusa);
  cs_field_t *f_mu = CS_F_(mu);
  cs_field_t *f_mut = CS_F_(mu_t);
  cs_field_t *f_rho = CS_F_(rho);

  cs_real_t *visct = f_mut->val;

  const cs_real_t *viscl = f_mu->val;
  const cs_real_t *crom  = f_rho->val;
  const cs_real_t *cvar_nusa  = (const cs_real_t *)f_nusa->val;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id ++) {

    cs_real_t xrom = crom[c_id];
    cs_real_t nusa = cvar_nusa[c_id];
    cs_real_t xi3 = cs_math_pow3(xrom * nusa / viscl[c_id]);
    cs_real_t fv1 = xi3 / (xi3 + cv13);

    visct[c_id] = xrom * nusa * fv1;

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
