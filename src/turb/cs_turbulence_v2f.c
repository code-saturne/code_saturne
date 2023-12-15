/*============================================================================
 * V2F turbulence model.
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

#include "cs_array.h"
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
#include "cs_porous_model.h"
#include "cs_prototypes.h"
#include "cs_time_step.h"
#include "cs_turbulence_model.h"
#include "cs_wall_functions.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_turbulence_v2f.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_turbulence_v2f.c

  Solve the \f$\phi\f$ and diffusion for \f$ \overline{f} \f$
  as part of the V2F phi-model.
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
 * \brief Clipping of \f$ v^2\f$ and \f$ \phi \f$ for the Bl v2/k
 * turbulence model (no clipping on \f$ f \f$).
 *
 *  \param[in]  n_cells       number of cells
 *  \param[in]  verbosity     verbosity level
 */
/*----------------------------------------------------------------------------*/

static void
_clip_v2f(cs_lnum_t  n_cells,
          int        verbosity)
{
  cs_real_t *cvar_phi = CS_F_(phi)->val;

  int kclipp = cs_field_key_id("clipping_id");

  cs_gnum_t nclp[2] =  {0, 0};  /* Min and max clipping values respectively */

  /* Postprocess clippings ? */
  cs_real_t *cpro_phi_clipped = NULL;
  int clip_phi_id = cs_field_get_key_int(CS_F_(phi), kclipp);
  if (clip_phi_id > -1) {
    cpro_phi_clipped = cs_field_by_id(clip_phi_id)->val;
    cs_array_real_fill_zero(n_cells, cpro_phi_clipped);
  }

  cs_real_t *cvar_al = NULL, *cpro_a_clipped = NULL;
  if (cs_glob_turb_model->iturb == CS_TURB_V2F_BL_V2K) {
    cvar_al = CS_F_(alp_bl)->val;
    int  clip_a_id = cs_field_get_key_int(CS_F_(alp_bl), kclipp);
    if (clip_a_id > -1) {
      cpro_a_clipped = cs_field_by_id(clip_a_id)->val;
      cs_array_real_fill_zero(n_cells, cpro_a_clipped);
    }
  }

  /* Store Min and Max for logging */

  cs_real_t vmin[2] = {cs_math_big_r, cs_math_big_r};
  cs_real_t vmax[2] = {-cs_math_big_r, -cs_math_big_r};

  for (cs_lnum_t i = 0; i < n_cells; i++) {
    cs_real_t var = cvar_phi[i];
    vmin[0] = cs_math_fmin(vmin[0], var);
    vmax[0] = cs_math_fmax(vmax[0], var);
  }

  if (cs_glob_turb_model->iturb == CS_TURB_V2F_BL_V2K) {
    vmin[1] = cs_math_big_r;
    vmax[1] = -cs_math_big_r;
    for (cs_lnum_t i = 0; i < n_cells; i++) {
      cs_real_t var = cvar_al[i];
      vmin[1] = cs_math_fmin(vmin[1], var);
      vmax[1] = cs_math_fmax(vmax[1], var);
    }
  }

  cs_parall_min(2, CS_REAL_TYPE, vmin);
  cs_parall_max(2, CS_REAL_TYPE, vmax);

  /* For the phi-fbar and BL-v2 / k model, identification of the phi values
     greater than 2 and phi clipping in absolute value for negative values
     ---------------------------------------------------------------------- */

  /* Identification of values ​​greater than 2, for display only */

  if (verbosity > 1) {
    for (cs_lnum_t i = 0; i < n_cells; i++) {
      if (cvar_phi[i] > 2)
        nclp[1] += 1;
    }
  }

  /* Clipping in absolute value for negative values */

  for (cs_lnum_t i = 0; i < n_cells; i++) {
    if (cvar_phi[i] < 0) {
      if (cpro_phi_clipped != NULL)
        cpro_phi_clipped[i] = cvar_phi[i];
      cvar_phi[i] = -cvar_phi[i];
      nclp[0] += 1;
    }
  }

  cs_parall_counter(nclp, 2);

  if (nclp[1] > 0)
    cs_log_warning(_("Warning: variable phi, maximum physical value of 2 "
                     "exceeded for, %llu cells"),
                  (unsigned long long)nclp[1]);

  cs_lnum_t nclpmn[1] = {nclp[0]}, nclpmx[1] = {nclp[1]};
  cs_log_iteration_clipping_field(CS_F_(phi)->id,
                                  nclpmn[0], 0,
                                  vmin, vmax,
                                  nclpmn, nclpmx);

  /* For the BL-v2 / k model, clipping from alpha to 0 for negative values
     and a 1 for values ​​greater than 1.
     --------------------------------------------------------------------- */

  if (cs_glob_turb_model->iturb == CS_TURB_V2F_BL_V2K) {

    nclp[0] = 0; nclp[1] = 0;

    for (cs_lnum_t i = 0; i < n_cells; i++) {
      if (cvar_al[i] < 0) {
        if (cpro_a_clipped != NULL)
          cpro_a_clipped[i] = -cvar_al[i];
        cvar_al[i] = 0;
        nclp[0] += 1;
      }
      if (cvar_al[i] > 1) {
        if (cpro_a_clipped != NULL)
          cpro_a_clipped[i] = 1.0 - cvar_al[i];
        cvar_al[i] = 1.0;
        nclp[1] += 1;
      }
    }

    cs_parall_counter(nclp, 2);

    nclpmn[0] = nclp[0], nclpmx[0] = nclp[1];
    cs_log_iteration_clipping_field(CS_F_(alp_bl)->id,
                                    nclpmn[0], nclpmx[0],
                                    vmin+1, vmax+1,
                                    nclpmn, nclpmx);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the term grad(phi).grad(k)
 *
 *  \param[in]      n_cells       number of cells
 *  \param[in]      n_cells_ext   number of extended (local + ghost) cells
 *  \param[in,out]  grad_pk       grad(phi).grad(k) values at cells
 */
/*----------------------------------------------------------------------------*/

static void
_gradfi_dot_gradk(const  cs_lnum_t   n_cells,
                  const  cs_lnum_t   n_cells_ext,
                  cs_real_t          grad_pk[])
{
  /* Allocate temporary arrays gradients */
  cs_real_3_t *grad_phi;
  cs_real_3_t *grad_k;

  BFT_MALLOC(grad_phi, n_cells_ext, cs_real_3_t);
  BFT_MALLOC(grad_k, n_cells_ext, cs_real_3_t);

  cs_field_gradient_scalar(CS_F_(phi),
                           true,     /* use previous t */
                           1,        /* not on increment */
                           grad_phi);

  cs_field_gradient_scalar(CS_F_(k),
                           true,     /* use previous t */
                           1,        /* not on increment */
                           grad_k);

  for (cs_lnum_t i = 0; i < n_cells; i++) {
    grad_pk[i] = cs_math_3_dot_product(grad_phi[i], grad_k[i]);
  }

  /* Free memory */
  BFT_FREE(grad_phi);
  BFT_FREE(grad_k);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve the quation of f_bar / alpha
 *
 * \param[in]          istprv      integer value associated with the key of
 *                                 source terms at previous time step
 * \param[in]          crom        density at current time step
 * \param[in]          cromo       density at previous (or current) time step
 * \param[in]          viscl       lam. visc. at current time step
 * \param[in]          cpro_pcvlo  lam. visc. at previous (or current) time step
 * \param[in]          prdv2f      prod of turbulence for the v2f
 * \param[in]          grad_pk     array to store the term grad(phi).grad(k)
 * \param[in]          i_massflux  mass flux at interior faces
 * \param[in]          b_massflux  mass flux at boundary faces
 * \param[in, out]     rhs         explicit source term
 * \param[in, out]     rovsdt      implicit part of the source term
 * \param[in, out]     c_st_a_p    second-order previous source term
 */
/*----------------------------------------------------------------------------*/

static void
_solve_eq_fbr_al(const int         istprv,
                 const cs_real_t   crom[],
                 const cs_real_t   cromo[],
                 const cs_real_t   viscl[],
                 const cs_real_t   cpro_pcvlo[],
                 const cs_real_t   prdv2f[],
                 const cs_real_t   grad_pk[],
                 const cs_real_t   i_massflux[],
                 const cs_real_t   b_massflux[],
                 cs_real_t         rhs[],
                 cs_real_t         rovsdt[],
                 cs_real_t         c_st_a_p[])
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_lnum_t  n_i_faces = m->n_i_faces;
  const cs_lnum_t  n_b_faces = m->n_b_faces;
  const cs_lnum_t  n_cells = m->n_cells;
  const cs_lnum_t  n_cells_ext = m->n_cells_with_ghosts;

  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  cs_real_t  *w2, *visel, *w5, *viscf, *viscb;

  /* Allocate temporary arrays */
  BFT_MALLOC(w2, n_cells_ext, cs_real_t);
  BFT_MALLOC(visel, n_cells_ext, cs_real_t);
  BFT_MALLOC(w5, n_cells_ext, cs_real_t);
  BFT_MALLOC(viscf, n_i_faces, cs_real_t);
  BFT_MALLOC(viscb, n_b_faces, cs_real_t);

  cs_field_t *f = NULL;

  if (cs_glob_turb_model->iturb == CS_TURB_V2F_PHI) {
    f = CS_F_(f_bar);
  }
  else if (cs_glob_turb_model->iturb == CS_TURB_V2F_BL_V2K) {
    f = CS_F_(alp_bl);
  }

  cs_real_t *cvar_var = f->val;
  const cs_real_t *cvara_var = f->val_pre;

  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);

  if  (eqp->verbosity >= 1) {
    const char *label = cs_field_get_label(f);
    cs_log_printf(CS_LOG_DEFAULT,
                  _("           Solving variable %s\n\n"), label);
  }

  const cs_time_scheme_t *time_scheme = cs_glob_time_scheme;
  const cs_real_t thets = time_scheme->thetst;
  const cs_real_t thetv = eqp->thetav;

  for (cs_lnum_t i = 0; i < n_cells; i++) {
    rhs[i] = 0;
    rovsdt[i] = 0;
  }

  /* Initialization of work arrays in case of HTLES */
  cs_real_t *htles_psi = NULL;
  cs_real_t *htles_r   = NULL;
  if (cs_glob_turb_model->hybrid_turb == 4) {
    htles_psi = cs_field_by_name("htles_psi")->val;
    htles_r   = cs_field_by_name("htles_r")->val;
  }

  /* User source terms
     ----------------- */

  cs_user_source_terms(cs_glob_domain,
                       f->id,
                       rhs,
                       rovsdt);

  if (cs_glob_porous_model == 3)
    cs_immersed_boundary_wall_functions(f->id, rhs, rovsdt);

  /* If we extrapolate the source terms */
  if (istprv >= 0) {
    for (cs_lnum_t i = 0; i < n_cells; i++) {
      const cs_real_t tuexpe = c_st_a_p[i];
      /* For the future and the next step time
         We put a mark "-" because in fact we solve
         \f$-\div{\grad {\dfrac{\overline{f}}{\alpha}}} = ... \f$ */
      c_st_a_p[i] = - rhs[i];
      /* Second member of the previous step time
         we implicit the user source term (the rest)*/
      rhs[i] = - rovsdt[i]*cvara_var[i] - thets*tuexpe;
      /* Diagonal */
      rovsdt[i] *= thetv;
    }
  }
  else {
    for (cs_lnum_t i = 0; i < n_cells; i++) {
      /*  We put a mark "-" because in fact we solve
          \f$-\div\{\grad{\dfrac{\overline{f}}{\alpha}}} = ...\f$
          We solve by conjugated gradient, so we do not impose the mark
          of rovsdt */
      rhs[i] = - rovsdt[i]*cvara_var[i] - rhs[i];
    }
  }

  /*  Source term of f_bar or alpha
      ----------------------------- */

  /*! For f_bar (phi_fbar)
   *    \f[ rhs =   \dfrac{1}{L^2*(f_b + \dfrac{1}{T(C1-1)(phi-2/3)}
   *               - \dfrac{C2 Pk}{k \rho}
   *               -2 \dfrac{\nu}{k\grad{\phi}\cdot \grad{k}
   *               -\nu \div{\grad{\phi} )} \f]
   *  For alpha (BL-V2/K)
   *    \f$rhs = \dfrac{1}{L^2 (\alpha^3 - 1)} \f$
   *  In fact we put a mark "-" because the solved equation is
   *     \f[ -\div{\grad{ \dfrac{\overline{f}}{\alpha}}} = rhs \f]
   */

  for (cs_lnum_t i = 0; i < n_cells; i++) {
    visel[i] = 1.0;
  }

  int imvisf = cs_glob_space_disc->imvisf;
  cs_face_viscosity(m,
                    fvq,
                    imvisf,
                    visel,
                    viscf,
                    viscb);

  /* Translate coefa into cofaf and coefb into cofbf */
  cs_real_t *coefap = CS_F_(phi)->bc_coeffs->a;
  cs_real_t *coefbp = CS_F_(phi)->bc_coeffs->b;
  cs_real_t *cofafp = CS_F_(phi)->bc_coeffs->af;
  cs_real_t *cofbfp = CS_F_(phi)->bc_coeffs->bf;

  const cs_real_t *distb = fvq->b_dist;
  const cs_lnum_t *b_face_cells = m->b_face_cells;

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

    const cs_lnum_t c_id = b_face_cells[face_id];

    const cs_real_t hint = visel[c_id]/distb[face_id];

    /* Translate coefa into cofaf and coefb into cofbf */
    cofafp[face_id] = -hint*coefap[face_id];
    cofbfp[face_id] = hint*(1.-coefbp[face_id]);

  }

  const cs_equation_param_t *eqp_phi
    = cs_field_get_equation_param_const(CS_F_(phi));

  cs_diffusion_potential(f->id,
                         m,
                         fvq,
                         1,     /* init */
                         1,     /* inc */
                         eqp_phi->imrgra,
                         eqp_phi->nswrgr,
                         eqp_phi->imligr,
                         0,     /* iphydp */
                         eqp_phi->iwgrec,
                         eqp_phi->verbosity,
                         eqp_phi->epsrgr,
                         eqp_phi->climgr,
                         NULL,
                         CS_F_(phi)->val_pre,
                         coefap,
                         coefbp,
                         cofafp,
                         cofbfp,
                         viscf,
                         viscb,
                         visel,
                         w2);

  /* Explicit term, stores ke temporarily in w5
     w2 is already multiplied by the volume which already contains
     a mark "-" (coming from cs_diffusion_potential) */

  const cs_real_t *cvara_k = CS_F_(k)->val_pre;
  const cs_real_t *cvara_ep = CS_F_(eps)->val_pre;
  const cs_real_t *cvara_phi = CS_F_(phi)->val_pre;

  const cs_real_t *cell_f_vol = fvq->cell_f_vol;

  const cs_real_t d2s3 = 2.0/3.0;

  if (cs_glob_turb_model->iturb == CS_TURB_V2F_PHI) {

    for (cs_lnum_t i = 0; i < n_cells; i++) {
      /* Compute the time scale*/
      const cs_real_t x_k = cvara_k[i];
      const cs_real_t x_e = cvara_ep[i];
      const cs_real_t x_rho = cromo[i];
      const cs_real_t x_nu = cpro_pcvlo[i]/x_rho;
      const cs_real_t tt_ke = x_k / x_e;
      const cs_real_t tt_min = cs_turb_cv2fct*sqrt(x_nu/x_e);
      const cs_real_t time_scale = cs_math_fmax(tt_ke, tt_min);

      const cs_real_t c1 =  cs_turb_cv2fc1 - 1.;
      const cs_real_t x_p = cvara_phi[i] - d2s3;
      const cs_real_t v2f = cs_turb_cv2fc2*prdv2f[i]/x_rho/x_k;
      const cs_real_t x_es = 2.*x_nu/x_e/time_scale * grad_pk[i];
      w5[i] = -cell_f_vol[i]*(c1*x_p/time_scale - v2f - x_es) - x_nu*w2[i];
    }

  }
  else if (cs_glob_turb_model->iturb == CS_TURB_V2F_BL_V2K) {

    /* FIXME the computed terms are not used */

    for (cs_lnum_t i = 0; i < n_cells; i++) {
#if 0
      /* Compute the time scale*/
      const cs_real_t x_k = cvara_k[i];
      const cs_real_t x_e = cvara_ep[i];
      const cs_real_t x_rho = crom[i];
      const cs_real_t x_nu = cpro_pcvlo[i];
      const cs_real_t tt_ke = x_k/x_e;
      const cs_real_t tt_min = cs_turb_cpalct*sqrt(x_nu/x_e);
      const cs_real_t time_scale = sqrt(  cs_math_pow2(tt_ke)
                                        + cs_math_pow2(tt_min));
#endif
      w5[i] = cell_f_vol[i];
    }
  }

  /* If we extrapolate the source term: */
  cs_real_t thetap;

  if (istprv >= 0) {
    cs_real_t thetp1 = 1.0 + thets;
    for (cs_lnum_t i = 0; i < n_cells; i++) {
      c_st_a_p[i] += w5[i];
      rhs[i] += thetp1*c_st_a_p[i];
    }
    thetap = thetv;
  }
  /* Otherwise: rhs */
  else {
    for (cs_lnum_t i = 0; i < n_cells; i++) {
      rhs[i] += w5[i];
    }
    thetap = 1.0;
  }

  cs_real_t l2 = 1;
  /* Implicit term */
  for (cs_lnum_t i = 0; i < n_cells; i++) {

    /* Compute L^2 */
    const cs_real_t x_k = cvara_k[i];
    const cs_real_t x_e = cvara_ep[i];
    const cs_real_t x_nu = viscl[i]/crom[i];

    if (cs_glob_turb_model->iturb == CS_TURB_V2F_PHI) {
      const cs_real_t ll_ke = pow(x_k, 1.5)/x_e;
      const cs_real_t ll_min
        = cs_turb_cv2fet*pow(cs_math_pow3(x_nu)/x_e, 0.25);
      l2 = cs_math_pow2(cs_turb_cv2fcl*cs_math_fmax(ll_ke, ll_min));
    }
    else if (cs_glob_turb_model->iturb == CS_TURB_V2F_BL_V2K) {
      if (cs_glob_turb_model->hybrid_turb == 4) {
      /* HTLES method */
        const cs_real_t x_psi  = htles_psi[i];
        const cs_real_t x_r    = htles_r[i];
        const cs_real_t ll_ke  = pow(x_k, 1.5)/(x_psi*x_e);
        const cs_real_t ll_min = pow(x_r, 1.5)
          * cs_turb_cpalet*pow(cs_math_pow3(x_nu)/(x_psi*x_e), 0.25);
        l2 =   cs_math_pow2(cs_turb_cpalcl)
            * (cs_math_pow2(ll_ke) + cs_math_pow2(ll_min));
      }
      else {
        const cs_real_t ll_ke = pow(x_k, 1.5)/x_e;
        const cs_real_t ll_min
          = cs_turb_cpalet*pow(cs_math_pow3(x_nu)/x_e, 0.25);
        l2 =   cs_math_pow2(cs_turb_cpalcl)
            * (cs_math_pow2(ll_ke) + cs_math_pow2(ll_min));
      }
    }
    rhs[i] = (- cell_f_vol[i] * cvara_var[i] + rhs[i]) / l2;

    /* Matrix*/
    rovsdt[i] = (rovsdt[i] + cell_f_vol[i]*thetap) / l2;
  }

  /* Effective resolution in the equation of f_bar or alpha
     ------------------------------------------------------ */

  coefap = f->bc_coeffs->a;
  coefbp = f->bc_coeffs->b;
  cofafp = f->bc_coeffs->af;
  cofbfp = f->bc_coeffs->bf;

  /* Solve current variable */

  cs_equation_param_t eqp_loc = *eqp;

  eqp_loc.istat  = -1;
  eqp_loc.icoupl = -1;
  eqp_loc.idifft = -1;
  eqp_loc.iwgrec = 0;     /* Warning, may be overwritten if a field */
  eqp_loc.thetav = thetv;
  eqp_loc.blend_st = 0;   /* Warning, may be overwritten if a field */

  cs_equation_iterative_solve_scalar(cs_glob_time_step_options->idtvar,
                                     1,    /* init */
                                     f->id,
                                     f->name,
                                     0,   /* iescap */
                                     0,   /* imucpp */
                                     -1,  /* normp */
                                     &eqp_loc,
                                     cvara_var,
                                     cvara_var,
                                     coefap,
                                     coefbp,
                                     cofafp,
                                     cofbfp,
                                     i_massflux,
                                     b_massflux,
                                     viscf,
                                     viscb,
                                     viscf,
                                     viscb,
                                     NULL,
                                     NULL,
                                     NULL,
                                     0,  /* boundary convective upwind flux */
                                     NULL,
                                     rovsdt,
                                     rhs,
                                     cvar_var,
                                     w2, /* dpvar work array */
                                     NULL,
                                     NULL);

  /* Free memory */
  BFT_FREE(w2);
  BFT_FREE(visel);
  BFT_FREE(w5);
  BFT_FREE(viscf);
  BFT_FREE(viscb);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve the equation of phi
 *
 * \param[in]        istprv       value associated with the key of source
 *                                terms at previous time step
 * \param[in]        ncesmp       number of cells with mass source term
 * \param[in]        icetsm       number of cells with mass source
 * \param[in]        itypsm       type of masss source for the variables
 * \param[in]        dt           time step (per cell)
 * \param[in]        smacel       value of variables associated to the
 *                                mass source
 * \param[in]        crom         density at current time step
 * \param[in]        cromo        density at previous (or current) time step
 * \param[in]        viscl        lam. visc. at current time step
 * \param[in]        prdv2f       storage table of term
 * \param[in]        grad_pk      Array to store the term grad(phi).grad(k)
 *                                prod of turbulence for the v2f
 * \param[in]        i_massflux   mass flux at interior faces
 * \param[in]        b_massflux   mass flux at boundary faces
 * \param[in, out]   rhs          explicit source term
 * \param[in, out]   rovsdt       implicit part of the source term
 * \param[in, out]   c_st_phi_p
 */
/*----------------------------------------------------------------------------*/

static void
_solve_eq_phi(const int           istprv,
              cs_lnum_t           ncesmp,
              cs_lnum_t           icetsm[],
              int                 itypsm[],
              const cs_real_t    *dt,
              cs_real_t           smacel[],
              const cs_real_t     crom[],
              const cs_real_t     cromo[],
              const cs_real_t     viscl[],
              const cs_real_t     prdv2f[],
              const cs_real_t     grad_pk[],
              const cs_real_t     i_massflux[],
              const cs_real_t     b_massflux[],
              cs_real_t           rhs[],
              cs_real_t           rovsdt[],
              cs_real_t           c_st_phi_p[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t *b_face_cells = m->b_face_cells;
  const cs_lnum_t  n_i_faces = m->n_i_faces;
  const cs_lnum_t  n_b_faces = m->n_b_faces;
  const cs_lnum_t  n_cells = m->n_cells;
  const cs_lnum_t  n_cells_ext = m->n_cells_with_ghosts;

  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;
  const cs_real_t *cell_f_vol = fvq->cell_f_vol;

  cs_real_t *w2, *viscf, *viscb;

  /* Allocate temporary arrays */
  BFT_MALLOC(w2, n_cells_ext, cs_real_t);
  BFT_MALLOC(viscf, n_i_faces, cs_real_t);
  BFT_MALLOC(viscb, n_b_faces, cs_real_t);

  cs_field_t *f_phi = CS_F_(phi);

  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f_phi);

  if (eqp->verbosity >= 1) {
    const char *label = cs_field_get_label(f_phi);
    cs_log_printf(CS_LOG_DEFAULT,
                  _("           Solving variable %s\n\n"), label);
  }

  cs_real_t *cvar_phi = f_phi->val;
  cs_real_t *cvara_phi = f_phi->val_pre;
  const cs_real_t *cvara_k = CS_F_(k)->val_pre;
  const cs_real_t *cvara_ep = CS_F_(eps)->val_pre;

  cs_real_t *visct = CS_F_(mu_t)->val;
  cs_real_t *cpro_pcvto = CS_F_(mu_t)->val;
  if (istprv >= 0) {
    int key_t_ext_id = cs_field_key_id("time_extrapolated");
    int iviext = cs_field_get_key_int(CS_F_(mu_t), key_t_ext_id);
    if (iviext > 0) {
      cpro_pcvto = CS_F_(mu_t)->val_pre;
    }
  }

  const cs_time_scheme_t *time_scheme = cs_get_glob_time_scheme();
  const cs_real_t thets  = time_scheme->thetst;
  const cs_real_t thetv = eqp->thetav;

  for (cs_lnum_t i = 0; i < n_cells; i++) {
    rhs[i] = 0;
    rovsdt[i] = 0;
  }

  /* Initialization of work arrays in case of HTLES */
  cs_real_t *htles_psi = NULL;
  if (cs_glob_turb_model->hybrid_turb == 4) {
    htles_psi = cs_field_by_name("htles_psi")->val;
  }

  /* User source terms
     ----------------- */

  cs_user_source_terms(cs_glob_domain,
                       f_phi->id,
                       rhs,
                       rovsdt);

  if (cs_glob_porous_model == 3)
    cs_immersed_boundary_wall_functions(f_phi->id, rhs, rovsdt);

  /* If we extrapolate the source terms */
  if (istprv >= 0) {
    for (cs_lnum_t i = 0; i < n_cells; i++) {
      const cs_real_t tuexpe = c_st_phi_p[i];  /* save for exchange */

      /* For the following and the next time step */
      c_st_phi_p[i] = rhs[i];
      /* Right-hand side of the previous step time
         We suppose -rovsdt > 0: we implicit
         the user source term (the rest)*/
      rhs[i] = rovsdt[i]*cvara_phi[i] - thets*tuexpe;
      /* Diagonal */
      rovsdt[i] *= -thetv;
    }
  }
  else {
    for (cs_lnum_t i = 0; i < n_cells; i++) {
      rhs[i] += rovsdt[i]*cvara_phi[i];
      rovsdt[i] = cs_math_fmax(-rovsdt[i], 0.);
    }
  }

  /* Mass source term
     ---------------- */

  if (ncesmp > 0) {

    const int var_key_id = cs_field_key_id("variable_id");
    int ivar_phi = cs_field_get_key_int(f_phi, var_key_id) - 1;
    int ivar_pr = cs_field_get_key_int(CS_F_(p), var_key_id) - 1;

    /* We increment rhs by -Gamma.var_prev and rovsdt by Gamma */
    cs_mass_source_terms(1, /* iterns*/
                         1,  /* dim*/
                         ncesmp,
                         icetsm,
                         itypsm + ncesmp*ivar_phi,
                         cell_f_vol,
                         cvara_phi,
                         smacel + ncesmp*ivar_phi,
                         smacel + ncesmp*ivar_pr,
                         rhs,
                         rovsdt,
                         w2);

    /* If we extrapolate the source term we put Gamma Pinj in the prev. TS */
    if (istprv >= 0) {
      for (cs_lnum_t i = 0; i < n_cells; i++) {
        c_st_phi_p[i] += w2[i];
      }
    }
    else{
      for (cs_lnum_t i = 0; i < n_cells; i++) {
        rhs[i] += w2[i];
      }
    }
  }

  /* Mass accumulation term \f$ -\dfrad{dRO}{dt}VOLUME \f$
   * and unstable over time term
   * ----------------------------------------------------- */

  /* Adding the matrix diagonal */

  if (eqp->istat == 1) {
    for (cs_lnum_t i = 0; i < n_cells; i++)
      rovsdt[i] += (crom[i] / dt[i]) * cell_f_vol[i];
  }

  /* Source term of phi
   * \f$ \phi_fbar\f$:
   * \f[ rhs = \rho f_bar - \dfrac{\phi}{k} P_k +\dfrac{2}{k}
   *                         \dfrac{\mu_t}{\sigma_k} \grad{\phi} \cdot \grad{k} \f]
   * BL-V2/K:
   * \f[ rhs = \rho \alpha f_h + \rho (1-\alpha^p) f_w - \dfrac{\phi}{k} P_k
   *            +\dfrac{2}{k} \dfrac{\mu_t}{\sigma_k} \grad{\phi} \cdot \grad{k} \f]
   *    with \f$ f_w=-\dfrac{\epsilon}{2} \cdot \dfrac{\phi}{k} \f$ and
   *         \f$ f_h = \dfrac{1}{T} \cdot
   *                            (C1-1+C2 \dfrac{P_k}{\epsilon \rho} (2/3-\phi) \f$
   */

  /* Explicit term, store temporarily in w2 */

  const cs_real_t d2s3 = 2.0/3.0;
  cs_real_t sigmak = cs_field_get_key_double
                       (CS_F_(k), cs_field_key_id("turbulent_schmidt"));

  if (cs_glob_turb_model->iturb == CS_TURB_V2F_PHI) {
    const cs_real_t *cvar_fb = CS_F_(f_bar)->val;

    for (cs_lnum_t i = 0; i < n_cells; i++) {
      const cs_real_t x_k = cvara_k[i];
      const cs_real_t x_rho = cromo[i];

      /* The term in f_bar is taken at the current and not previous time step
       * ... a priori better.
       * Remark: if we keep this choice, we have to modify the case
       *         of the second-order (which need the previous value time step
       *         for extrapolation). */
      w2[i] = cell_f_vol[i] * (  x_rho*cvar_fb[i]
                               + 2./x_k*cpro_pcvto[i]/sigmak*grad_pk[i]);
    }
  }
  else if (cs_glob_turb_model->iturb == CS_TURB_V2F_BL_V2K) {
    const cs_real_t *cvara_al = CS_F_(alp_bl)->val_pre;

    for (cs_lnum_t i = 0; i < n_cells; i++) {
      const cs_real_t x_k = cvara_k[i];
      const cs_real_t x_e = cvara_ep[i];
      const cs_real_t x_rho = cromo[i];
      const cs_real_t x_nu = viscl[i] / crom[i];

      if (cs_glob_turb_model->hybrid_turb == 4) {
        /* HTLES method */
        const cs_real_t x_psi = htles_psi[i];
        const cs_real_t tt_ke = x_k/(x_psi*x_e);
        const cs_real_t tt_min = cs_turb_cpalct*sqrt(x_nu/(x_psi*x_e));
        const cs_real_t tt = sqrt(cs_math_pow2(tt_ke) + cs_math_pow2(tt_min));
        const cs_real_t fhomog = -1.0/tt*(  cs_turb_cpalc1-1.0
                                          +   cs_turb_cpalc2*prdv2f[i]
                                            / (x_psi*x_e)/x_rho)
                                        *(cvara_phi[i] - d2s3);
        w2[i] = cell_f_vol[i] * (cs_math_pow3(cvara_al[i])*fhomog*x_rho
                               + 2./x_k *cpro_pcvto[i]/sigmak*grad_pk[i]);
        /* FIXME implicit negative w1 and fhomog */
      }
      else {
        const cs_real_t tt_ke = x_k/x_e;
        const cs_real_t tt_min = cs_turb_cpalct*sqrt(x_nu/x_e);
        const cs_real_t tt = sqrt(cs_math_pow2(tt_ke) + cs_math_pow2(tt_min));
        const cs_real_t fhomog = -1.0/tt*(  cs_turb_cpalc1-1.0
                                          + cs_turb_cpalc2*prdv2f[i]/x_e/x_rho)
                                        *(cvara_phi[i] - d2s3);
        w2[i] = cell_f_vol[i] * (cs_math_pow3(cvara_al[i])*fhomog*x_rho
                               + 2./x_k *cpro_pcvto[i]/sigmak*grad_pk[i]);
        /* FIXME implicit negative w1 and fhomog */
      }
    }
  }

  /* If we extrapolate the source term: prev. ST */

  cs_real_t thetap = (istprv >= 0) ? thetv : 1.;

  if (istprv >= 0) {
    const cs_real_t  thetp1 = 1.0 + thets;
    for (cs_lnum_t i = 0; i < n_cells; i++) {
      c_st_phi_p[i] += w2[i];
      rhs[i] += thetp1*c_st_phi_p[i];
    }
  }
  /* Otherwise: rhs */
  else {
    for (cs_lnum_t i = 0; i < n_cells; i++) {
      rhs[i] += w2[i];
    }
  }

  /* Implict term (rhs) and matrix (rovsdt) */

  if (cs_glob_turb_model->iturb == CS_TURB_V2F_PHI) {

    for (cs_lnum_t i = 0; i < n_cells; i++) {
      const cs_real_t prdv2f_m = cs_math_fmax(prdv2f[i], 0.0);
      rhs[i] -= cell_f_vol[i] * prdv2f[i] * cvara_phi[i] / cvara_k[i];
      rovsdt[i] += cell_f_vol[i] * prdv2f_m / cvara_k[i] * thetap;
    }
  }

  else if (cs_glob_turb_model->iturb == CS_TURB_V2F_BL_V2K) {
    const cs_real_t *cvara_al = CS_F_(alp_bl)->val_pre;

    for (cs_lnum_t i = 0; i < n_cells; i++) {

      const cs_real_t x_rho = cromo[i];
      const cs_real_t prdv2f_m = cs_math_fmax(prdv2f[i], 0.0);
      const cs_real_t al_3 = cs_math_pow3(cvara_al[i]);

      if (cs_glob_turb_model->hybrid_turb == 4) {
        /* HTLES method */
        const cs_real_t x_psi = htles_psi[i];
        rhs[i] -= cell_f_vol[i] * (  prdv2f[i] + x_rho*(x_psi*cvara_ep[i])/2.
                                   * (1. - al_3))
                                * cvar_phi[i] / cvara_k[i];
        rovsdt[i] += cell_f_vol[i] * (  prdv2f_m + x_rho*(x_psi*cvara_ep[i])/2.
                                      * (1. - al_3))
                                   / cvara_k[i] * thetap;
      }
      else {
        rhs[i] -= cell_f_vol[i] * (  prdv2f[i] + x_rho*cvara_ep[i]/2.
                                   * (1. - al_3))
                                * cvar_phi[i] / cvara_k[i];
        rovsdt[i] += cell_f_vol[i] * (  prdv2f_m + x_rho*cvara_ep[i]/2.
                                      * (1. - al_3))
                                   / cvara_k[i] * thetap;
      }
    }
  }

  /* Diffusion terms
     --------------- */

  /* Viscosity
   * Normally, in the phi-model equations, only turbulent viscosity
   *  turbulente takes place in phi diffusion (the term with mu disappeared
   *  passing from \f$f\f$ to \f$ \overline{f})\f$. But as it stands,
   *  it makes the calculation unstable (because \f$\mu_t\f$ tends towards 0
   *  at the wall what decouples \f$ \phi \f$ of its boundary condition and
   *  the molecular diffusion term is integred in \f$ \overline{f} \f$,
   *  it is as if it was treated as explicit).
   *  -> we add artificially diffusion (knowing that as k=0, the phi value
   *  does not matter)
   */

  cs_real_t *coefap = f_phi->bc_coeffs->a;
  cs_real_t *coefbp = f_phi->bc_coeffs->b;
  cs_real_t *cofafp = f_phi->bc_coeffs->af;
  cs_real_t *cofbfp = f_phi->bc_coeffs->bf;

  if (eqp->idiff >= 1) {

    if (cs_glob_turb_model->iturb == CS_TURB_V2F_PHI) {
      for (cs_lnum_t i = 0; i < n_cells; i++) {
        w2[i] = viscl[i] + visct[i]/sigmak;
      }
    }
    else if (cs_glob_turb_model->iturb == CS_TURB_V2F_BL_V2K) {
      for (cs_lnum_t i = 0; i < n_cells; i++) {
        w2[i] = viscl[i]/2. + visct[i]/sigmak; /* FIXME */
      }
    }

    int imvisf = cs_glob_space_disc->imvisf;
    cs_face_viscosity(m,
                      fvq,
                      imvisf,
                      w2,
                      viscf,
                      viscb);

    const cs_real_t *distb = fvq->b_dist;

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

      const cs_lnum_t c_id = b_face_cells[face_id];

      const cs_real_t hint = w2[c_id]/distb[face_id];

      /* Translate coefa into cofaf and coefb into cofbf */
      cofafp[face_id] = -hint*coefap[face_id];
      cofbfp[face_id] = hint*(1.-coefbp[face_id]);

    }

  }
  else {

    for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
      viscf[face_id] = 0.;
    }

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      viscb[face_id] = 0.;

      /* Translate coefa into cofaf and coefb into cofbf */
      cofafp[face_id] = 0.;
      cofbfp[face_id] = 0.;
    }

  }

  /* Effective resolution of the phi equation
     ---------------------------------------- */

  if (istprv >= 0) {
    const cs_real_t thetp1 = 1.0 + thets;
    for (cs_lnum_t i = 0; i < n_cells; i++) {
      rhs[i] += thetp1 * c_st_phi_p[i];
    }
  }

  coefap = f_phi->bc_coeffs->a;
  coefbp = f_phi->bc_coeffs->b;
  cofafp = f_phi->bc_coeffs->af;
  cofbfp = f_phi->bc_coeffs->bf;

  cs_equation_param_t eqp_loc = *eqp;

  cs_equation_iterative_solve_scalar(cs_glob_time_step_options->idtvar,
                                     1,    /* init */
                                     f_phi->id,
                                     f_phi->name,
                                     0,   /* iescap */
                                     0,   /* imucpp */
                                     -1,  /* normp */
                                     &eqp_loc,
                                     cvara_phi,
                                     cvara_phi,
                                     coefap,
                                     coefbp,
                                     cofafp,
                                     cofbfp,
                                     i_massflux,
                                     b_massflux,
                                     viscf,
                                     viscb,
                                     viscf,
                                     viscb,
                                     NULL,
                                     NULL,
                                     NULL,
                                     0,   /* boundary convective upwind flux */
                                     NULL,
                                     rovsdt,
                                     rhs,
                                     cvar_phi,
                                     w2,  /* dpvar work array */
                                     NULL,
                                     NULL);

  BFT_FREE(w2);
  BFT_FREE(viscf);
  BFT_FREE(viscb);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the V2F phi-model equations.
 *
 * Solve the \f$\phi\f$ and diffusion for \f$ \overline{f} \f$
 * as part of the V2F phi-model
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
 * \param[in]     prdv2f        v2f production term
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_v2f(cs_lnum_t         ncesmp,
                  cs_lnum_t         icetsm[],
                  int               itypsm[],
                  const cs_real_t  *dt,
                  cs_real_t         smacel[],
                  const cs_real_t   prdv2f[])
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  const cs_equation_param_t *eqp_phi
    = cs_field_get_equation_param_const(CS_F_(phi));

  cs_real_t *c_st_a_p = NULL;
  cs_real_t *c_st_phi_p = NULL;

  cs_real_t *rhs;
  cs_real_t *rovsdt;

  BFT_MALLOC(rhs, n_cells_ext, cs_real_t);
  BFT_MALLOC(rovsdt, n_cells_ext, cs_real_t);

  /* Map field arrays */

  int kstprv = cs_field_key_id("source_term_prev_id");
  int istprv = cs_field_get_key_int(CS_F_(phi), kstprv);

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  int iflmas =  cs_field_get_key_int(CS_F_(phi), kimasf);
  int iflmab =  cs_field_get_key_int(CS_F_(phi), kbmasf);

  const cs_real_t *imasfl =  cs_field_by_id(iflmas)->val;
  const cs_real_t *bmasfl =  cs_field_by_id(iflmab)->val;

  const cs_real_t *crom = CS_F_(rho)->val;
  const cs_real_t *cromo = CS_F_(rho)->val;
  const cs_real_t *viscl = CS_F_(mu)->val;
  const cs_real_t *cpro_pcvlo = CS_F_(mu)->val;

  if (istprv >= 0) {
    int key_t_ext_id = cs_field_key_id("time_extrapolated");
    int iroext = cs_field_get_key_int(CS_F_(rho), key_t_ext_id);
    if (iroext > 0) {
      cromo = CS_F_(rho)->val_pre;
    }

    int iviext = cs_field_get_key_int(CS_F_(mu), key_t_ext_id);
    if (iviext > 0) {
      cpro_pcvlo = CS_F_(mu)->val_pre;
    }
  }

  /* 2nd-order previous source terms */

  if (istprv >= 0) {
    c_st_phi_p = cs_field_by_id(istprv)->val;
    if (cs_glob_turb_model->iturb == CS_TURB_V2F_PHI) {
      istprv = cs_field_get_key_int(CS_F_(f_bar), kstprv);
      if (istprv >= 0)
        c_st_a_p = cs_field_by_id(istprv)->val;
    }
  }
  else if (cs_glob_turb_model->iturb == CS_TURB_V2F_BL_V2K) {
    istprv = cs_field_get_key_int(CS_F_(alp_bl), kstprv);
    if (istprv >= 0)
      c_st_a_p = cs_field_by_id(istprv)->val;
  }

  if (eqp_phi->verbosity >= 1) {
    cs_log_printf(CS_LOG_DEFAULT,
                  "\n"
                  "   ** Solving V2F (phi and f_bar/alpha)\n"
                  "      ---------------------------------\n\n");
  }

  /* Compute grad(phi).grad(k) */

  cs_real_t *grad_pk;
  BFT_MALLOC(grad_pk, n_cells_ext, cs_real_t);

  _gradfi_dot_gradk(n_cells, n_cells_ext, grad_pk);

  /* Solve the equation of f_bar / alpha */

  _solve_eq_fbr_al(istprv,
                   crom,
                   cromo,
                   viscl,
                   cpro_pcvlo,
                   prdv2f,
                   grad_pk,
                   imasfl,
                   bmasfl,
                   rhs,
                   rovsdt,
                   c_st_a_p);

  /* Solve the equation of phi */

  _solve_eq_phi(istprv,
                ncesmp,
                icetsm,
                itypsm,
                dt,
                smacel,
                crom,
                cromo,
                viscl,
                prdv2f,
                grad_pk,
                imasfl,
                bmasfl,
                rhs,
                rovsdt,
                c_st_phi_p);

  BFT_FREE(grad_pk);
  BFT_FREE(rhs);
  BFT_FREE(rovsdt);

  /* Clipping */

  _clip_v2f(n_cells, eqp_phi->verbosity);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of turbulent viscosity for the V2F-phi model.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_v2f_phi_mu_t(void)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  cs_real_t *visct =  CS_F_(mu_t)->val;

  const cs_real_t *viscl = CS_F_(mu)->val;
  const cs_real_t *crom = CS_F_(rho)->val;
  const cs_real_t *cvar_k = CS_F_(k)->val;
  const cs_real_t *cvar_ep = CS_F_(eps)->val;
  const cs_real_t *cvar_phi = CS_F_(phi)->val;

  /* HTLES method */
  if (cs_glob_turb_model->hybrid_turb == 4) {

    // cs_real_t *psi = cs_field_by_name("htles_psi")->val;
    // cs_real_t *blend = cs_field_by_name("hybrid_blend")->val;

    //TODO VD
    bft_error(__FILE__, __LINE__, 0,
              _("%s: not implemented for hybrid_turb = 4."), __func__);
  }
  else {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      const cs_real_t xk   = cvar_k[c_id];
      const cs_real_t xe   = cvar_ep[c_id];
      const cs_real_t xrom = crom[c_id];
      const cs_real_t xnu  = viscl[c_id]/xrom;

      const cs_real_t ttke = xk / xe;
      const cs_real_t ttmin = cs_turb_cv2fct * sqrt(xnu/xe);
      const cs_real_t tt = fmax(ttke, ttmin);

      visct[c_id] = cs_turb_cmu*xrom*tt*cvar_phi[c_id]*xk;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of turbulent viscosity for the V2F-BL model.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_v2f_bl_v2k_mu_t(void)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  /* Initialization
   * ============== */

  /* Map field arrays */

  cs_real_t *visct =  CS_F_(mu_t)->val;
  const cs_real_t *viscl  =  (const cs_real_t *)CS_F_(mu)->val;
  const cs_real_t *crom  = CS_F_(rho)->val;

  const cs_real_t *cvar_k = CS_F_(k)->val;
  const cs_real_t *cvar_ep = CS_F_(eps)->val;
  const cs_real_t *cvar_phi = CS_F_(phi)->val;

  /* Calculation of velocity gradient and of
   *       S2 = S11^2+S22^2+S33^2+2*(S12^2+S13^2+S23^2)
   * ================================================== */

  cs_real_t *s2;
  cs_real_33_t *gradv;

  /* Allocate arrays */
  BFT_MALLOC(s2, n_cells_ext, cs_real_t);
  BFT_MALLOC(gradv, n_cells_ext, cs_real_33_t);

  cs_field_gradient_vector(CS_F_(vel),
                           false, // no use_previous_t
                           1,    // inc
                           gradv);

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

    s2[c_id] =  2.*(cs_math_pow2(s11)+cs_math_pow2(s22)+cs_math_pow2(s33))
              + cs_math_pow2(dudy+dvdx)
              + cs_math_pow2(dudz+dwdx)
              + cs_math_pow2(dvdz+dwdy);
    s2[c_id] = sqrt(fmax(s2[c_id], 1e-10));
  }

  /* Free memory */
  BFT_FREE(gradv);

  /* Calculation of viscosity
   * ========================= */

  const cs_real_t f1 = 0.6 / sqrt(3.) / cs_turb_cmu;

  /* HTLES method */
  if (cs_glob_turb_model->hybrid_turb == 4) {

    cs_real_t *psi = cs_field_by_name("htles_psi")->val;
    cs_real_t *blend = cs_field_by_name("hybrid_blend")->val;

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      const cs_real_t xk   = cvar_k[c_id];
      const cs_real_t xe   = cvar_ep[c_id];
      const cs_real_t xrom = crom[c_id];
      const cs_real_t xnu  = viscl[c_id]/xrom;

      const cs_real_t ttke = xk / xe;
      /* Modif. definition Kolmogorov Length scale */
      const cs_real_t ttmin = cs_turb_cv2fct * sqrt(xnu/(psi[c_id] * xe));

      /*  ft1 is not taken into account in LES mode (1/(1-xrc)->infty) */
      cs_real_t xfs2 = (s2[c_id] * CS_MAX(cs_math_epzero, 1. - blend[c_id]));
      const cs_real_t ft1 = f1 / xfs2;
      const cs_real_t ft2
        = sqrt(cs_math_pow2(ttke) + cs_math_pow2(ttmin))*cvar_phi[c_id];

      visct[c_id] = cs_turb_cmu*xrom*cvar_k[c_id]*fmin(ft1, ft2) / psi[c_id];
    }

  }
  else {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      const cs_real_t xk   = cvar_k[c_id];
      const cs_real_t xe   = cvar_ep[c_id];
      const cs_real_t xrom = crom[c_id];
      const cs_real_t xnu  = viscl[c_id]/xrom;

      const cs_real_t ttke = xk/xe;
      const cs_real_t ttmin = cs_turb_cpalct*sqrt(xnu/xe);

      /* We initially have:
       * ttlim = 0.6/cvar_phi(iel)/sqrt(3)/cmu/s2(iel)
       * tt = min(ttlim, sqrt(ttke^2 + ttmin^2))
       * visct(iel) = cmu*xrom*tt*cvar_phi(iel)*cvar_k(iel)
       *
       * When tt = ttlim, tt in
       *   visct(iel) = cmu*xrom*tt*cvar_phi(iel)*cvar_k(iel)
       * cvar_phi appears in both numerator and denominator,
       * and can be eliminated. */

      const cs_real_t ft1 = f1/s2[c_id];
      const cs_real_t ft2
        = sqrt(cs_math_pow2(ttke) + cs_math_pow2(ttmin))*cvar_phi[c_id];

      visct[c_id] = cs_turb_cmu*xrom*cvar_k[c_id]*fmin(ft1, ft2);

    }
  }

  /* Free memory */
  BFT_FREE(s2);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
