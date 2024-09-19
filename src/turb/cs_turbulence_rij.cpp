/*============================================================================
 * Rij-epsilon turbulence model.
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <chrono>

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_array.h"
#include "cs_base.h"
#include "cs_blas.h"
#include "cs_boundary_conditions.h"
#include "cs_dispatch.h"
#include "cs_domain.h"
#include "cs_equation_iterative_solve.h"
#include "cs_equation_param.h"
#include "cs_face_viscosity.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_gradient.h"
#include "cs_lagr.h"
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
#include "cs_timer.h"
#include "cs_turbomachinery.h"
#include "cs_turbulence_bc.h"
#include "cs_turbulence_model.h"
#include "cs_volume_mass_injection.h"
#include "cs_velocity_pressure.h"
#include "cs_wall_functions.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_turbulence_rij.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Tensor to vector (t2v) and vector to tensor (v2t) mask arrays */

const cs_lnum_t _t2v[3][3] = {{0, 3, 5}, {3, 1, 4}, {5, 4, 2}};

const cs_lnum_t _iv2t[6] = {0, 1, 2, 0, 1, 0};
const cs_lnum_t _jv2t[6] = {0, 1, 2, 1, 2, 2};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute all eigenvalues of a 3x3 symmetric matrix
 *        with symmetric storage.
 *
 * Based on: Oliver K. Smith "eigenvalues of a symmetric 3x3 matrix",
 *           Communication of the ACM (April 1961)
 *           (Wikipedia article entitled "Eigenvalue algorithm")
 *
 * \param[in]  m          3x3 symmetric matrix (m11, m22, m33, m12, m23, m13)
 * \param[out] eig_vals   size 3 vector
 */
/*----------------------------------------------------------------------------*/

CS_F_HOST_DEVICE void
_sym_33_eigen(const cs_real_t  m[6],
              cs_real_t        eig_vals[3])
{
  constexpr cs_real_t c_1ov3 = 1./3.;
  constexpr cs_real_t c_1ov6 = 1./6.;

  cs_real_t  e, e1, e2, e3;

  cs_real_t  p1 = cs_math_3_square_norm((const cs_real_t *)(m+3));
  cs_real_t  d2 = cs_math_3_square_norm((const cs_real_t *)m);

  cs_real_t  tr = (m[0] + m[1] + m[2]);
  cs_real_t  tr_third = c_1ov3 * tr;

  e1 = m[0] - tr_third, e2 = m[1] - tr_third, e3 = m[2] - tr_third;
  cs_real_t  p2 = e1*e1 + e2*e2 + e3*e3 + 2.*p1;

  cs_real_t  p = sqrt(p2*c_1ov6);

  if (p1 > cs_math_epzero*d2 && p > 0.) { /* m is not diagonal */

    cs_real_6_t  n;
    cs_real_t  ovp = 1./p;

    for (int  i = 0; i < 3; i++) {
      /* Diagonal */
      n[i] = ovp * (m[i] - tr_third);
      /* Extra diagonal */
      n[3 + i] = ovp * m[3 + i];
    }

    /* r should be between -1 and 1 but truncation error and bad conditioning
       can lead to slightly under/over-shoot */
    cs_real_t  r = 0.5 * cs_math_sym_33_determinant(n);

    cs_real_t  cos_theta, cos_theta_2pi3;
    if (r <= -1.) {
      cos_theta = 0.5; /* theta = pi/3; */
      cos_theta_2pi3 = -1.;
    }
    else if (r >= 1.) {
      cos_theta = 1.; /* theta = 0.; */
      cos_theta_2pi3 = -0.5;
    }
    else {
      cos_theta = cos(c_1ov3*acos(r));
      cos_theta_2pi3 = cos(c_1ov3*(acos(r) + 2.*cs_math_pi));
    }

    /* eigenvalues computed should satisfy e1 < e2 < e3 */
    e3 = tr_third + 2.*p*cos_theta;
    e1 = tr_third + 2.*p*cos_theta_2pi3;
    e2 = tr - e1 -e3; // since tr(m) = e1 + e2 + e3

  }
  else { /* m is diagonal */

    e1 = m[0], e2 = m[1], e3 = m[2];

  } /* diagonal or not */

  if (e3 < e2) e = e3, e3 = e2, e2 = e;
  if (e3 < e1) e = e3, e3 = e1, e1 = e2, e2 = e;
  else {
    if (e2 < e1) e = e2, e2 = e1, e1 = e;
  }

  /* Return values */
  eig_vals[0] = e1;
  eig_vals[1] = e2;
  eig_vals[2] = e3;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief returns the value of A with the sign of B.
 *
 *  \param[in]  a  real A
 *  \param[in]  b  real B
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_sign(cs_real_t  a,
      cs_real_t  b)

{
  cs_real_t sgn = (b < 0) ? - 1 : 1;

  return (sgn * cs_math_fabs(a));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Clipping of alpha in the framwork of the Rij-EBRSM model.
 *
 *  \param[in]  f_id          field id of alpha variable
 *  \param[in]  n_cells       number of cells
 *  \param[in]  alpha_min     minimum acceptable value for alpha
 */
/*----------------------------------------------------------------------------*/

static void
_clip_alpha(const int          f_id,
            const cs_lnum_t    n_cells,
            const cs_real_t    alpha_min[])
{
  cs_real_t *cvar_al = cs_field_by_id(f_id)->val;

  int kclipp = cs_field_key_id("clipping_id");
  cs_lnum_t nclp[2] =  {0, 0};  /* Min and max clipping values respectively */

  /* Postprocess clippings ? */
  cs_real_t *cpro_a_clipped = nullptr;
  int clip_a_id = cs_field_get_key_int(cs_field_by_id(f_id), kclipp);
  if (clip_a_id > -1) {
    cpro_a_clipped = cs_field_by_id(clip_a_id)->val;
    cs_array_real_fill_zero(n_cells, cpro_a_clipped);
  }

  /* Store local min and max for logging */
  cs_real_t vmin[1] = {cs_math_big_r};
  cs_real_t vmax[1] = {-cs_math_big_r};

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_real_t var = cvar_al[c_id];
    vmin[0] = cs_math_fmin(vmin[0], var);
    vmax[0] = cs_math_fmax(vmax[0], var);
  }

  /* Clipping (edit to avoid exactly zero values) */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    if (cvar_al[c_id] < alpha_min[c_id]) {
      if (clip_a_id > -1)
        cpro_a_clipped[c_id] = alpha_min[c_id]-cvar_al[c_id];
      nclp[0] += 1;
      cvar_al[c_id] = alpha_min[c_id];
    }
    else if (cvar_al[c_id] >= 1) {
      if (clip_a_id > -1)
        cpro_a_clipped[c_id] = cvar_al[c_id]-1.0;
      nclp[1] += 1;
      cvar_al[c_id] = 1.0;
    }
  }

  cs_lnum_t iclpmn[1] = {nclp[0]}, iclpmx[1] = {nclp[1]};
  cs_log_iteration_clipping_field(f_id,
                                  iclpmn[0],
                                  iclpmx[0],
                                  vmin, vmax,
                                  iclpmn, iclpmx);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute (rank-local) minima and maxima of Rij and epsilon.
 *
 * \param[in]   n_cells   number of cells
 * \param[in]   cvar_rij  Rij values
 * \param[in]   cvar_ep   epsilon values
 * \param[out]  vmin      local minima for Rij (0-5) and epsilon (6)
 * \param[out]  vmax      local maxima for Rij (0-5) and epsilon (6)
 */
/*----------------------------------------------------------------------------*/

static void
_rij_min_max(cs_lnum_t        n_cells,
             const cs_real_t  cvar_rij[][6],
             const cs_real_t  cvar_ep[],
             cs_real_t        vmin[7],
             cs_real_t        vmax[7])
{
  for (int ii = 0; ii < 7; ii++) {
    vmin[ii] = cs_math_big_r;
    vmax[ii] = -cs_math_big_r;
  }

# pragma omp parallel if(n_cells > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n_cells, sizeof(cs_real_t), &s_id, &e_id);

    cs_real_t t_vmin[7], t_vmax[7];
    for (int ii = 0; ii < 7; ii++) {
      t_vmin[ii] = cs_math_big_r;
      t_vmax[ii] = -cs_math_big_r;
    }

    for (cs_lnum_t c_id = s_id; c_id < e_id; c_id++) {
      for (cs_lnum_t ii = 0; ii < 6; ii++) {
        t_vmin[ii] = cs_math_fmin(t_vmin[ii], cvar_rij[c_id][ii]);
        t_vmax[ii] = cs_math_fmax(t_vmax[ii], cvar_rij[c_id][ii]);
      }
    }
    for (cs_lnum_t c_id = s_id; c_id < e_id; c_id++) {
      t_vmin[6] = cs_math_fmin(t_vmin[6], cvar_ep[c_id]);
      t_vmax[6] = cs_math_fmax(t_vmax[6], cvar_ep[c_id]);
    }

    #pragma omp critical
    {
      for (int ii = 0; ii < 7; ii++) {
        vmin[ii] = cs_math_fmin(vmin[ii], t_vmin[ii]);
        vmax[ii] = cs_math_fmax(vmax[ii], t_vmax[ii]);
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute \f$ \overline{\vect{u}' \rho'}\f$
 *
 * \param[in]   phase_id  turbulent phase id (-1 for single phase flow)
 * \param[out]  up_rhop   correlation fluctuation of velocity and density
 */
/*----------------------------------------------------------------------------*/

static void
_compute_up_rhop(int                 phase_id,
                 cs_real_t           up_rhop[][3])
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  cs_dispatch_context ctx;

  cs_field_t *f_rij = CS_F_(rij);
  cs_field_t *f_eps = CS_F_(eps);
  cs_field_t *f_rho = CS_F_(rho);
  cs_field_t *f_rho_b = CS_F_(rho_b);

  if (phase_id >= 0) {
    f_rij = CS_FI_(rij, phase_id);
    f_eps = CS_FI_(eps, phase_id);
    f_rho = CS_FI_(rho, phase_id);
    f_rho_b = CS_FI_(rho_b, phase_id);
  }

  const cs_real_t *cvara_ep = f_eps->val_pre;
  const cs_real_6_t *cvara_rij = (const cs_real_6_t *)f_rij->val_pre;
  cs_real_t *rho = f_rho->val;
  cs_real_t *rho_b = f_rho_b->val;

  const cs_equation_param_t *eqp
    = cs_field_get_equation_param_const(f_rij);

  cs_real_t cons = -1.5 * cs_turb_cmu;

  const cs_field_t *thf = cs_thermal_model_field();
  if (thf != nullptr) {
    const int ksigmas = cs_field_key_id("turbulent_schmidt");
    const cs_real_t turb_schmidt = cs_field_get_key_double(thf, ksigmas);
    cons = -1.5 * cs_turb_cmu / turb_schmidt;
  }
  const cs_velocity_pressure_model_t  *vp_model
    = cs_glob_velocity_pressure_model;
  const int idilat = vp_model->idilat;

  if (phase_id >= 0)
    thf = CS_FI_(h_tot, phase_id);

  const int kturt = cs_field_key_id("turbulent_flux_model");
  int turb_flux_model = 0;
  if (thf != nullptr)
    turb_flux_model = cs_field_get_key_int(thf, kturt);
  const int turb_flux_model_type = turb_flux_model / 10;

  cs_field_t *f_hf = nullptr;
  if (thf != nullptr)
    f_hf = cs_field_by_composite_name_try(thf->name, "turbulent_flux");

  /* We want to use the computed turbulent heat fluxes (when they exist)
     when we are not dealing with atmo cases */
  if (   (f_hf != nullptr && cs_glob_physical_model_flag[CS_ATMOSPHERIC] < 1)
      || turb_flux_model_type == 3
      || turb_flux_model_type == 2 ) {

    /* Using thermal fluxes
     * (only for thermal scalar for the moment) */
    cs_field_t *f_beta = cs_field_by_name_try("thermal_expansion");

    if (f_beta != nullptr) {
      /* Value of the corresponding turbulent flux */
      cs_real_3_t *xut = (cs_real_3_t *)f_hf->val;

      const cs_real_t *beta = f_beta->val;

      /* finalize rho'u' = -rho beta T'u' = -rho beta C k/eps R.gradT */
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

        const cs_real_t factor = - beta[c_id] * rho[c_id];

        for (cs_lnum_t i = 0; i < 3; i++)
          up_rhop[c_id][i] = factor * xut[c_id][i];

      });
      ctx.wait();
    }
    else { // beta = 0
      cs_arrays_set_value<cs_real_t, 1>(3*n_cells, 0., (cs_real_t *)up_rhop);
    }
  }

  /* Note that the buoyant term is normally expressed in term of
   * (u'T') or (u'rho') here modelled with a GGDH:
   *   (u'rho') = C * k/eps * R_ij Grad_j(rho)
   */

  else {

    /* Boussinesq or anelastic approximation:
     * buoyant terms in function of gradT */
    if (   idilat == 0
        || cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 1) {

      /* TODO take all buoyant scalars */

      if (thf != nullptr) {

        /* Using thermal fluxes
         * (only for thermal scalar for the moment) */
        cs_field_t *f_beta = cs_field_by_name_try("thermal_expansion");

        if (f_beta != nullptr) {

          /* finalize rho'u' = -rho beta T'u' = -rho beta C k/eps R.gradT */
          const cs_real_t *beta = f_beta->val;

          cs_real_3_t *gradt;
          CS_MALLOC_HD(gradt, n_cells_ext, cs_real_3_t, cs_alloc_mode);

          cs_field_gradient_scalar(thf,
                                   false,  /* use current (not previous) value */
                                   1,      /* inc */
                                   gradt);

          ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
            cs_real_t rit[3];
            cs_math_sym_33_3_product(cvara_rij[c_id], gradt[c_id], rit);

            /* factor = - rho beta C k/eps */
            const cs_real_t factor = - beta[c_id] * rho[c_id] * cons
                                       * cs_math_6_trace(cvara_rij[c_id])
                                       / (2. * cvara_ep[c_id]);

            for (cs_lnum_t i = 0; i < 3; i++)
              up_rhop[c_id][i] = factor * rit[i];
          });
          ctx.wait();

          CS_FREE_HD(gradt);
        }

        else { // beta = 0
          cs_arrays_set_value<cs_real_t, 1>(3*n_cells, 0., (cs_real_t *)up_rhop);
        }
      }

    }

    /* Not Boussinesq nor anelastic */
    else {

      cs_real_3_t *gradro;
      CS_MALLOC_HD(gradro, n_cells_ext, cs_real_3_t, cs_alloc_mode);

      /* Boundary conditions: Dirichlet romb
       * We impose in Dirichlet (coefa) the value romb */
      cs_real_t *coefb;
      CS_MALLOC_HD(coefb, n_b_faces, cs_real_t, cs_alloc_mode);

      cs_arrays_set_value<cs_real_t, 1>(n_b_faces, 0., coefb);

      cs_field_bc_coeffs_t bc_coeffs_loc;
      cs_field_bc_coeffs_init(&bc_coeffs_loc);
      bc_coeffs_loc.a = rho_b;
      bc_coeffs_loc.b = coefb;

      /* Compute gradient */
      cs_halo_type_t halo_type = CS_HALO_STANDARD;
      cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;
      cs_gradient_type_by_imrgra(eqp->imrgra,
                                 &gradient_type,
                                 &halo_type);

      cs_gradient_scalar("density",
                         gradient_type,
                         halo_type,
                         1,             /* inc */
                         eqp->nswrgr,
                         0,             /* iphydp */
                         1,             /* w_stride */
                         eqp->verbosity,
                         (cs_gradient_limit_t)(eqp->imligr),
                         eqp->epsrgr,
                         eqp->climgr,
                         nullptr,          /* f_ext */
                         &bc_coeffs_loc,
                         rho,
                         nullptr,         /* c_weight */
                         nullptr,         /* cpl */
                         gradro);

      /* finalize rho'u' */
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

        cs_real_t rit[3];
        cs_math_sym_33_3_product(cvara_rij[c_id], gradro[c_id], rit);

        const cs_real_t k_ov_eps =   cs_math_6_trace(cvara_rij[c_id])
                                   / (2.*cvara_ep[c_id]);

        for (cs_lnum_t i = 0; i < 3; i++)
          up_rhop[c_id][i] = cons * k_ov_eps * rit[i];

      });

      ctx.wait();

      CS_FREE_HD(gradro);
      CS_FREE_HD(coefb);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Terms of wall echo for \f$ R_{ij} \f$
 *        \f$var  = R_{11} \: R_{22} \: R_{33} \: R_{12} \: R_{13} \: R_{23}\f$
 *        \f$comp =  1 \:  2 \:  3 \:  4 \:  5 \:  6\f$
 *
 * \param[in]     phase_id    turbulent phase id (-1 for single phase flow)
 * \param[in]     produc      production
 * \param[in,out] rhs         work array for right-hand-side
 */
/*----------------------------------------------------------------------------*/

static void
_rij_echo(int              phase_id,
          const cs_real_t  produc[][6],
          cs_real_t        rhs[][6])
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  const cs_real_t *cell_f_vol = cs_glob_mesh_quantities->cell_f_vol;

  cs_dispatch_context ctx;

  cs_field_t *f_rij = CS_F_(rij);
  cs_field_t *f_eps = CS_F_(eps);
  cs_field_t *f_rho = CS_F_(rho);

  if (phase_id >= 0) {
    f_rij = CS_FI_(rij, phase_id);
    f_eps = CS_FI_(eps, phase_id);
    f_rho = CS_FI_(rho, phase_id);
  }

  const cs_real_t *cvara_ep = (const cs_real_t *)f_eps->val_pre;
  const cs_real_6_t *cvara_rij = (const cs_real_6_t *)f_rij->val_pre;

  cs_real_t *cromo = nullptr;
  int key_t_ext_id = cs_field_key_id("time_extrapolated");
  int iroext = cs_field_get_key_int(f_rho, key_t_ext_id);
  if ((cs_glob_time_scheme->isto2t > 0) && (iroext > 0))
    cromo = f_rho->val_pre;
  else
    cromo = f_rho->val;

  const cs_real_t d2s3 = 2./3;
  const cs_real_t cmu075 = pow(cs_turb_cmu, 0.75);

  /* Calculation in the orthogonal straight cells in corresponding walls
   * ------------------------------------------------------------------- */

  const cs_real_t *w_dist = cs_field_by_name("wall_distance")->val;

  cs_real_3_t *grad;
  CS_MALLOC_HD(grad, n_cells_ext, cs_real_3_t, cs_alloc_mode);

  /* Current gradient */
  cs_field_gradient_scalar(cs_field_by_name("wall_distance"),
                           false,  /* use_previous_t */
                           1,      /* inc */
                           grad);

  /* Tensor indexing */

  const cs_real_t crijp1 = cs_turb_crijp1;
  const cs_real_t crijp2 = cs_turb_crijp2;
  const cs_real_t crij2 = cs_turb_crij2;
  const cs_real_t xkappa = cs_turb_xkappa;

  const auto t2v = _t2v;
  const auto iv2t = _iv2t;
  const auto jv2t = _jv2t;
  const cs_real_t tdeltij[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    cs_real_t w6[6] = {0., 0., 0., 0., 0., 0.};

    cs_real_t n[3];
    cs_math_3_normalize(grad[c_id], n);
    /* Production and k */
    cs_real_t prod_k = 0.5 * cs_math_6_trace(produc[c_id]);

    const cs_real_t xk = 0.5 * cs_math_6_trace(cvara_rij[c_id]);
    cs_real_t epsk = cvara_ep[c_id] / xk;

    for (cs_lnum_t k = 0; k < 3; k++) {
      for (cs_lnum_t l = 0; l < 3; l++) {

        const cs_lnum_t kl = t2v[k][l];

        /* Terms with R km and Phi km;

           nk * nl * deltij[jj] in original expression,
           with deltij Kronecker delta so 1 for terms 0-2, 0 for terms 3-6.
           (we only have diagonal contributions) */

        const cs_real_t c_tr
          = n[k] * n[l]
            * (crijp1 * cvara_rij[c_id][kl] * epsk
               - (  crijp2 * crij2
                  * (produc[c_id][kl] - d2s3 * prod_k * tdeltij[k][l])));

        w6[0] += c_tr;
        w6[1] += c_tr;
        w6[2] += c_tr;

      } /* End of loop on l */

      for (cs_lnum_t ij = 0; ij < 6; ij++) {

        cs_lnum_t i = iv2t[ij];
        cs_lnum_t j = jv2t[ij];

        const cs_lnum_t ki = t2v[k][i];
        const cs_lnum_t kj = t2v[k][j];

        w6[ij] +=   1.5 * n[k]
                    * (- (  crijp1
                           * (  cvara_rij[c_id][ki] * n[j]
                              + cvara_rij[c_id][kj] * n[i])
                           * epsk)
                       + (  crijp2 * crij2
                          * (  ( produc[c_id][ki]
                                -d2s3 * prod_k * tdeltij[k][i]) * n[j]
                             + ( produc[c_id][kj]
                                -d2s3 * prod_k * tdeltij[k][j]) * n[i])));

      } /* End of loop on ij */

    } /* End of loop on k */

    /* Distance to the wall and amortization function:
     *   For each calculation mode: same code, test
     *   Apart from the loop */

    const cs_real_t distxn = cs_math_fmax(w_dist[c_id], cs_math_epzero); //FIXME
    const cs_real_t trrij = 0.5 * cs_math_6_trace(cvara_rij[c_id]);
    cs_real_t bb =   cmu075 * pow(trrij, 1.5)
                   / (xkappa * cvara_ep[c_id] * distxn);
    bb = cs_math_fmin(bb, 1.0);

    for (cs_lnum_t ij = 0; ij < 6; ij++)
      rhs[c_id][ij] += cromo[c_id] * cell_f_vol[c_id] * w6[ij] * bb;

  }); /* End of loop on cells */

  ctx.wait();

  CS_FREE_HD(grad);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Gravity terms for \f$R_{ij}\f$
 *
 * \param[in]   f_rij     pointer to Rij field
 * \param[in]   up_rhop   work array for \f$ \vect{u}'\rho' \f$
 * \param[in]   grav      gravity
 * \param[in]   st_prv_id id of the previous source term in case of
 *                        time extrapolation
 * \param[in]   c_st_prv  pointer to the previous source term id
 * \param[out]  fimp      implicit contribution
 * \param[out]  rhs       explicit right-hand side
 */
/*----------------------------------------------------------------------------*/

static void
_gravity_st_rij(const cs_field_t  *f_rij,
                const cs_real_t    up_rhop[][3],
                const cs_real_t    grav[],
                int                st_prv_id,
                cs_real_6_t       *c_st_prv,
                cs_real_66_t      *fimp,
                cs_real_6_t       *rhs)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_real_t *cell_f_vol = mq->cell_f_vol;

  const int coupled_components = cs_glob_turb_rans_model->irijco;

  const cs_real_t g = cs_math_3_norm(grav);
  const cs_real_t o_m_crij3 = (1. - cs_turb_crij3);
  const cs_real_t crij3 = cs_turb_crij3;

  const auto t2v = _t2v;

  cs_real_6_t *_buoyancy = nullptr, *cpro_buoyancy = nullptr;
  cs_field_t *f_buo = cs_field_by_name_try("algo:rij_buoyancy");

  if (f_buo != nullptr) {
    cpro_buoyancy = (cs_real_6_t*)f_buo->val;
  }
  else {
    CS_MALLOC_HD(_buoyancy, n_cells_ext, cs_real_6_t, cs_alloc_mode);
    cpro_buoyancy = _buoyancy;
  }
  cs_dispatch_context ctx;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

     cs_real_t gij[3][3];
     for (cs_lnum_t i = 0; i < 3; i++) {
       for (cs_lnum_t j = 0; j < 3; j++)
         gij[i][j] = up_rhop[c_id][i]*grav[j] + up_rhop[c_id][j]*grav[i];
     }

     const cs_real_t gkks3 = cs_math_33_trace(gij) / 3.;

     cpro_buoyancy[c_id][0] = gij[0][0] * o_m_crij3 + crij3*gkks3;
     cpro_buoyancy[c_id][1] = gij[1][1] * o_m_crij3 + crij3*gkks3;
     cpro_buoyancy[c_id][2] = gij[2][2] * o_m_crij3 + crij3*gkks3;
     cpro_buoyancy[c_id][3] = gij[0][1] * o_m_crij3;
     cpro_buoyancy[c_id][4] = gij[1][2] * o_m_crij3;
     cpro_buoyancy[c_id][5] = gij[0][2] * o_m_crij3;

     /* If we extrapolate the source terms: previous ST */
     if (st_prv_id > -1) {
       for (cs_lnum_t ij = 0; ij < 6; ij++)
         c_st_prv[c_id][ij] += cpro_buoyancy[c_id][ij] * cell_f_vol[c_id];
     }
     /* Otherwise RHS */
     else {
       for (cs_lnum_t ij = 0; ij < 6; ij++)
        rhs[c_id][ij] += cpro_buoyancy[c_id][ij] * cell_f_vol[c_id];
     }

  });

  /* Partial implicitation in case of coupled component solver */

  if (coupled_components != 0) {

    const cs_real_6_t *cvara_var = (const cs_real_6_t *)f_rij->val_pre;

    cs_lnum_t solid_stride = 1;
    int *c_is_solid_zone_flag = cs_solid_zone_flag(cs_glob_mesh);
    const int c_is_solid_ref[1] = {0};
    if (c_is_solid_zone_flag == nullptr) {
      if (cs_alloc_mode > CS_ALLOC_HOST) {
        CS_MALLOC_HD(c_is_solid_zone_flag, 1, int, cs_alloc_mode);
        c_is_solid_zone_flag[0] = 0;
      }
      solid_stride = 0;
    }
    const int *c_is_solid = c_is_solid_zone_flag;
    if (c_is_solid == nullptr)
      c_is_solid = c_is_solid_ref;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

      const cs_lnum_t ind = solid_stride*c_id;
      if (c_is_solid[ind])
        return;  /* return from lambda function == continue in loop */

      cs_real_t implmat2add[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

      const cs_real_t gkks3 = cs_math_3_dot_product(up_rhop[c_id], grav) / 3.;

      /* Compute inverse matrix of R^n
         (scaling by tr(R) for numerical stability) */
      cs_real_t matrn[6];
      const cs_real_t trrij = 0.5 * cs_math_6_trace(cvara_var[c_id]);
      for (int ij = 0; ij < 6; ij++)
        matrn[ij] = cvara_var[c_id][ij] / trrij;

      cs_real_t oo_matrn[6];
      cs_math_sym_33_inv_cramer(matrn, oo_matrn);
      for (int ij = 0; ij < 6; ij++)
        oo_matrn[ij] /= trrij;

      if (gkks3 <= 0) {
        /* Term "C3 tr(G) Id"
           Compute inverse matrix of R^n
           (scaling by tr(R) for numerical stability) */

        for (cs_lnum_t i = 0; i < 3; i++) {
          for (cs_lnum_t j = 0; j < 3; j++) {
            cs_lnum_t ij = t2v[i][j];
            implmat2add[i][j] = -0.5 * gkks3 * oo_matrn[ij];
          }
        }
      }

      /* Identity constant, remove the negative eigen value of Gij
       * which is:
       * g.rho'u' - g ||rho'u'||
       */
      cs_real_t impl_id_cst =   g * cs_math_3_norm(up_rhop[c_id])
                              - cs_math_3_dot_product(grav, up_rhop[c_id]);

      for (cs_lnum_t i = 0; i < 3; i++) {
        for (cs_lnum_t j = 0; j < 3; j++) {
          cs_lnum_t ij = t2v[i][j];
          implmat2add[i][j] += 0.5 * impl_id_cst * oo_matrn[ij];
        }
      }

      /* Compute the 6x6 matrix A which verifies
       * A.R = M.R + R.M^t */
      cs_real_t impl_drsm[6][6];
      cs_math_reduce_sym_prod_33_to_66(implmat2add, impl_drsm);

      for (cs_lnum_t ij = 0; ij < 6; ij++) {
        for (cs_lnum_t kl = 0; kl < 6; kl++)
          fimp[c_id][ij][kl] += cell_f_vol[c_id] * impl_drsm[ij][kl];
      }

    }); /* End of loop on cells */

    ctx.wait();

    CS_FREE_HD(c_is_solid_zone_flag);
  } /* End of test on coupled components */

  ctx.wait();

  CS_FREE_HD(_buoyancy);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Gravity terms for \f$\epsilon\f$.
 *
 *  Terms for epsilon:
 *      rom*volumr*deps/dt =
 *                     ... + CEPS1*(EPS/K)*Max(0,(Gkk/2))*volume
 *
 * \param[in]       phase_id    turbulent phase id (-1 for single phase flow)
 * \param[in]       up_rhop     work array for \f$ \vect{u}'\rho' \f$
 * \param[in]       cell_f_vol  cell fluid volume
 * \param[in]       grav        gravity
 * \param[in, out]  rhs         work array for right-hand side
 */
/*----------------------------------------------------------------------------*/

static void
_gravity_st_epsilon(int              phase_id,
                    const cs_real_t  up_rhop[][3],
                    const cs_real_t  cell_f_vol[],
                    const cs_real_t  grav[],
                    cs_real_t        rhs[])
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const int kturt = cs_field_key_id("turbulent_flux_model");
  const int krvarfl = cs_field_key_id("variance_dissipation");

  cs_dispatch_context ctx;

  cs_field_t *f_rij = CS_F_(rij);
  cs_field_t *f_eps = CS_F_(eps);
  cs_field_t *f_rho = CS_F_(rho);
  cs_field_t *f_mu = CS_F_(mu);
  if (phase_id >= 0) {
    f_rij = CS_FI_(rij, phase_id);
    f_eps = CS_FI_(eps, phase_id);
    f_rho = CS_FI_(rho, phase_id);
    f_mu = CS_FI_(mu, phase_id);
  }

  cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();
  cs_real_t cp0 = fp->cp0;

  /* Specific heat for Prandtl calculation */
  const cs_field_t *f_cp = cs_field_by_name_try("specific_heat");
  cs_real_t *cpro_cp = nullptr;
  if (f_cp != nullptr)
    cpro_cp = f_cp->val;

  const cs_real_6_t *cvara_rij = (const cs_real_6_t *)f_rij->val_pre;
  const cs_real_t *cvara_ep = f_eps->val_pre;
  const cs_real_t *crom = f_rho->val;
  const cs_real_t *viscl = f_mu->val;

  const cs_turb_model_type_t iturb
    = (cs_turb_model_type_t)cs_glob_turb_model->iturb;
  int dissip_buo_mdl = cs_glob_turb_rans_model->dissip_buo_mdl;

  const cs_real_t xct = cs_turb_xct;
  const cs_real_t ce1 = cs_turb_ce1;
  const cs_real_t ce3 = cs_turb_ce3;

  cs_field_t *f_alpha_theta = nullptr;

  cs_real_t rvarfl = 0.0;
  int turb_flux_model =  -1;

  cs_lnum_t l_viscls = 0; /* stride for uniform/local viscosity access */
  cs_real_t _visls_0 = -1;
  cs_real_t *viscls = nullptr;

  cs_field_t *f_t = cs_thermal_model_field();
  if (f_t != nullptr){
    turb_flux_model =  cs_field_get_key_int(f_t, kturt);

    cs_field_t *f_t_var = cs_field_get_variance(f_t);
    if (f_t_var != nullptr)
      rvarfl = cs_field_get_key_double(f_t_var, krvarfl);

    const int kivisl = cs_field_key_id("diffusivity_id");
    int ifcvsl = cs_field_get_key_int(f_t, kivisl);
    if (ifcvsl > -1) {
      viscls = cs_field_by_id(ifcvsl)->val;
      l_viscls = 1;
    }
    else {
      const int kvisls0 = cs_field_key_id("diffusivity_ref");
      _visls_0 = cs_field_get_key_double(f_t, kvisls0);
      viscls = &_visls_0;
      l_viscls = 0;
    }

    f_alpha_theta = cs_field_by_composite_name_try(f_t->name, "alpha");
  }

  cs_real_t *alpha_theta = nullptr;
  if (f_alpha_theta != nullptr)
    alpha_theta = f_alpha_theta->val;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    cs_real_t xcp = cp0;
    if (cpro_cp != nullptr)
      xcp = cpro_cp[c_id];

    const cs_real_t g_up_rhop = cs_math_3_dot_product(grav, up_rhop[c_id]);

    const cs_real_t tke = 0.5 * cs_math_6_trace(cvara_rij[c_id]);

    cs_real_t time_scale = tke / cvara_ep[c_id];

    cs_real_t buoyancy_constant = ce1;

    if (iturb == CS_TURB_RIJ_EPSILON_EBRSM) {

      /* Calculation of the Durbin time scale */
      const cs_real_t xttkmg
        = xct*sqrt(viscl[c_id] / crom[c_id] / cvara_ep[c_id]);

      time_scale = cs_math_fmax(time_scale, xttkmg);

      /* Calculation of the mixed time scale for EB thermal
       *  flux models (Dehoux, 2017) */
      if (  (turb_flux_model == 11)  /* EB-GGDH */
         || (turb_flux_model == 21)  /* EB-AFM */
         || (turb_flux_model == 31)) /* EB-DFM */ {

        cs_real_t prdtl = viscl[c_id]*xcp/viscls[l_viscls*c_id];
        cs_real_t xr
          = (1.0 - alpha_theta[c_id])*prdtl + alpha_theta[c_id]*rvarfl;

        buoyancy_constant = ce3;

        time_scale = time_scale * sqrt(xr)/sqrt(prdtl);
      }
    }

    if (dissip_buo_mdl == 1) {
      rhs[c_id] +=   buoyancy_constant
                   * g_up_rhop/time_scale * cell_f_vol[c_id];
    }
    else {
      rhs[c_id] +=   buoyancy_constant
                   * cs_math_fmax(0., g_up_rhop/time_scale) * cell_f_vol[c_id];
    }
  });

  ctx.wait();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prepare the resolution of the coupled Reynolds stress components
 *        \f$ R_{ij} - \varepsilon \f$ RANS (LRR) turbulence model.
 *
 * \param[in]   f_rij     pointer to Rij field
 * \param[in]   phase_id  turbulent phase id (-1 for single phase flow)
 * \param[in]   gradv     work array for the velocity grad term
 * \param[in]   produc    work array for production
 * \param[in]   up_rhop   work array for \f$ \vect{u}'\rho' \f$
 * \param[in]   grav      gravity
 * \param[out]  viscf     visc*surface/dist at internal faces
 * \param[out]  viscb     visc*surface/dist at edge faces
 * \param[out]  viscce    Daly Harlow diffusion term
 * \param[out]  rhs       right hand side
 * \param[out]  fimp      Implicit term (containing unsteady term)
 * \param[out]  weighf    working array
 * \param[out]  weighb    working array
 */
/*----------------------------------------------------------------------------*/

static void
_pre_solve_lrr(const cs_field_t  *f_rij,
               int                phase_id,
               const cs_real_t    gradv[][3][3],
               const cs_real_t    produc[][6],
               const cs_real_t    up_rhop[][3],
               const cs_real_t    grav[],
               cs_real_t          viscf[],
               cs_real_t          viscb[],
               cs_real_t          viscce[][6],
               cs_real_t          rhs[][6],
               cs_real_t          fimp[][6][6],
               cs_real_t          weighf[][2],
               cs_real_t          weighb[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_real_t *cell_f_vol = fvq->cell_f_vol;

  cs_field_t *f_eps = CS_F_(eps);
  cs_field_t *f_rho = CS_F_(rho);
  cs_field_t *f_mu = CS_F_(mu);

  if (phase_id >= 0) {
    f_eps = CS_FI_(eps, phase_id);
    f_rho = CS_FI_(rho, phase_id);
    f_mu = CS_FI_(mu, phase_id);
  }

  const cs_real_t *crom = f_rho->val;
  const cs_real_t *viscl = f_mu->val;
  const cs_real_t *cvara_ep = f_eps->val_pre;
  const cs_real_6_t*cvara_var = (const cs_real_6_t *)f_rij->val_pre;

  const cs_equation_param_t *eqp
    = cs_field_get_equation_param_const(f_rij);

  if (eqp->verbosity >= 1) {
    bft_printf(" Solving the variable %s\n",
               cs_field_get_label(f_rij));
  }

  cs_real_6_t *c_st_prv = nullptr;
  int kstprv = cs_field_key_id("source_term_prev_id");
  int st_prv_id = cs_field_get_key_int(f_rij, kstprv);
  if (st_prv_id > -1)
    c_st_prv = (cs_real_6_t *)cs_field_by_id(st_prv_id)->val;

  /* Time extrapolation ? */
  cs_real_t *cromo = nullptr;
  int key_t_ext_id = cs_field_key_id("time_extrapolated");
  int iroext = cs_field_get_key_int(f_rho, key_t_ext_id);
  if ((iroext > 0) && (st_prv_id > -1))
    cromo = f_rho->val_pre;
  else
    cromo = f_rho->val;

  /* Coefficient of the "Coriolis-type" term */
  const int icorio = cs_glob_physical_constants->icorio;
  const cs_turbomachinery_model_t tm_model = cs_turbomachinery_get_model();
  cs_real_t ccorio = 0;
  if (icorio == 1)
    ccorio = 2; /* Relative velocity formulation */
  else if (tm_model == CS_TURBOMACHINERY_TRANSIENT)
     ccorio = 1;

  const cs_real_t d1s2 = 0.5;
  const cs_real_t d1s3 = 1./3;
  const cs_real_t d2s3 = 2./3;

  const cs_real_t crij1 = cs_turb_crij1;
  const cs_real_t crij2 = cs_turb_crij2;
  const cs_real_t crijeps = cs_turb_crij_eps;

  const auto t2v = _t2v;
  const cs_real_t tdeltij[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
  const cs_real_t vdeltij[6] = {1, 1, 1, 0, 0, 0};

  cs_lnum_t solid_stride = 1;
  int *c_is_solid_zone_flag = cs_solid_zone_flag(cs_glob_mesh);
  const int c_is_solid_ref[1] = {0};
  const int *c_is_solid = c_is_solid_zone_flag;
  if (c_is_solid == nullptr) {
    c_is_solid = c_is_solid_ref;
    solid_stride = 0;
  }

  /* Production, Pressure-Strain correlation, dissipation
   * ---------------------------------------------------- */

# pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    if (c_is_solid[solid_stride*c_id])
      continue;

    cs_real_t impl_lin_cst = 0, impl_id_cst = 0;

    cs_real_t trprod = 0.5 * cs_math_6_trace(produc[c_id]);
    cs_real_t trrij  = 0.5 * cs_math_6_trace(cvara_var[c_id]);

    cs_real_t xaniso[3][3], xstrai[3][3], xrotac[3][3];

    /* aij */
    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        xaniso[i][j] =   cvara_var[c_id][t2v[i][j]] / trrij
                       - d2s3 * tdeltij[i][j];
      }
    }

    /* Sij */
    xstrai[0][0] = gradv[c_id][0][0];
    xstrai[0][1] = d1s2 * (gradv[c_id][0][1] + gradv[c_id][1][0]);
    xstrai[0][2] = d1s2 * (gradv[c_id][0][2] + gradv[c_id][2][0]);
    xstrai[1][0] = xstrai[0][1];
    xstrai[1][1] = gradv[c_id][1][1];
    xstrai[1][2] = d1s2 * (gradv[c_id][1][2] + gradv[c_id][2][1]);
    xstrai[2][0] = xstrai[0][2];
    xstrai[2][1] = xstrai[1][2];
    xstrai[2][2] = gradv[c_id][2][2];

    /* omegaij */
    xrotac[0][0] = 0;
    xrotac[0][1] = d1s2 * (gradv[c_id][0][1] - gradv[c_id][1][0]);
    xrotac[0][2] = d1s2 * (gradv[c_id][0][2] - gradv[c_id][2][0]);
    xrotac[1][0] = -xrotac[0][1];
    xrotac[1][1] = 0;
    xrotac[1][2] = d1s2 * (gradv[c_id][1][2] - gradv[c_id][2][1]);
    xrotac[2][0] = -xrotac[0][2];
    xrotac[2][1] = -xrotac[1][2];
    xrotac[2][2] = 0;

    /* Computation of implicit components */
    cs_real_t sym_strain[6] = {xstrai[0][0],
                               xstrai[1][1],
                               xstrai[2][2],
                               xstrai[1][0],
                               xstrai[2][1],
                               xstrai[2][0]};

    /* Compute inverse matrix of R^n
       (scaling by tr(R) for numerical stability) */
    cs_real_t matrn[6];
    for (cs_lnum_t ij = 0; ij < 6; ij++)
      matrn[ij] = cvara_var[c_id][ij] / trrij;

    cs_real_t oo_matrn[6];
    cs_math_sym_33_inv_cramer(matrn, oo_matrn);
    for (cs_lnum_t ij = 0; ij < 6; ij++)
      oo_matrn[ij] /= trrij;

    /* Compute the maximal eigenvalue (in terms of norm!) of S */
    cs_real_t eigen_vals[3];
    cs_math_sym_33_eigen(sym_strain, eigen_vals);
    cs_real_t eigen_max = cs_math_fabs(eigen_vals[0]);
    for (cs_lnum_t i = 1; i < 3; i++)
      eigen_max = cs_math_fmax(cs_math_fabs(eigen_max),
                               cs_math_fabs(eigen_vals[i]));

    /* Constant for the dissipation */
    cs_real_t ceps_impl = d1s3 * cvara_ep[c_id];

    /* Identity constant */
    impl_id_cst = -d1s3 * crij2 * cs_math_fmin(trprod, 0);

    /* Linear constant */
    impl_lin_cst = eigen_max * (1.0 - crij2); /* Production + Phi2 */

    cs_real_t implmat2add[3][3];
    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        cs_lnum_t ij = t2v[i][j];
        implmat2add[i][j] =   (1.0 - crij2) * xrotac[i][j]
                              + impl_lin_cst * vdeltij[ij]
                              + impl_id_cst * d1s2 * oo_matrn[ij]
                              + ceps_impl * oo_matrn[ij];
      }
    }

    /* Compute the 6x6 matrix A which verifies
     * A.R = M.R + R.M^t */
    cs_real_t impl_drsm[6][6];
    cs_math_reduce_sym_prod_33_to_66(implmat2add, impl_drsm);

    /* Rotating frame of reference => "absolute" vorticity */
    if (icorio == 1) {
      cs_real_t matrot[3][3];
      const cs_rotation_t *r = cs_glob_rotation + 1;
      cs_rotation_add_coriolis_t(r, 1., matrot);
      for (cs_lnum_t i = 0; i < 3; i++) {
        for (cs_lnum_t j = 0; j < 3; j++)
          xrotac[i][j] -= matrot[i][j];
      }
    }

    for (cs_lnum_t ij = 0; ij < 6; ij++) {
      cs_lnum_t i = _iv2t[ij];
      cs_lnum_t j = _jv2t[ij];

      /* Explicit terms */
      cs_real_t pij = (1.-crij2) * produc[c_id][t2v[j][i]];
      cs_real_t phiij1 = -cvara_ep[c_id]*crij1*xaniso[j][i];
      cs_real_t phiij2 = d2s3*crij2*trprod*vdeltij[ij];
      cs_real_t epsij = -d2s3*crijeps*cvara_ep[c_id]*vdeltij[ij];

      if (st_prv_id > -1) {
        c_st_prv[c_id][ij] +=   cromo[c_id] * cell_f_vol[c_id]
                              * (pij + phiij1 + phiij2 + epsij);
      }
      else {
        rhs[c_id][ij] +=   cromo[c_id] * cell_f_vol[c_id]
                         * (pij + phiij1 + phiij2 + epsij);

        /* Implicit terms */
        fimp[c_id][ij][ij] +=   crom[c_id] * cell_f_vol[c_id] / trrij
                                * (crij1 * cvara_ep[c_id]);

        for (cs_lnum_t kl = 0; kl < 6; kl++)
          fimp[c_id][ij][kl] +=   crom[c_id] * cell_f_vol[c_id]
                                  * impl_drsm[ij][kl];
      }
    }

  } /* end loop on cells */

  BFT_FREE(c_is_solid_zone_flag);

  /* Coriolis terms in the Phi1 and production
   * ----------------------------------------- */

  cs_real_6_t *w2;
  BFT_MALLOC(w2, n_cells_ext, cs_real_6_t);

  if ((icorio == 1) || (tm_model == 1)) {

    const int *irotce = cs_turbomachinery_get_cell_rotor_num();

#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      int rot_id = icorio;
      if (tm_model == 1) {
        rot_id = irotce[c_id];
        if (rot_id < 1)
          continue;
      }

      cs_real_t matrot[3][3];
      const cs_rotation_t *r = cs_glob_rotation + rot_id;
      cs_rotation_coriolis_t(r, 1., matrot);

      /* Compute Gij: (i,j) component of the Coriolis production */
      for (cs_lnum_t ij = 0; ij < 6; ij++) {
        cs_lnum_t i = _iv2t[ij];
        cs_lnum_t j = _jv2t[ij];

        w2[c_id][ij] = 0.;
        for (cs_lnum_t k = 0; k < 3; k++)
          w2[c_id][ij] -=   ccorio
                           * (  matrot[k][i] * cvara_var[c_id][t2v[k][j]]
                              + matrot[k][j] * cvara_var[c_id][t2v[k][i]]);
      }

      /* Coriolis contribution in the Phi1 term: (1-C2/2)Gij */
      if (icorio == 1) {
        for (cs_lnum_t ij = 0; ij < 6; ij++)
          w2[c_id][ij] =   crom[c_id] * cell_f_vol[c_id]
                         * (1.-0.5*crij2) * w2[c_id][ij];
      }

      /* If source terms are extrapolated */
      if (st_prv_id > -1) {
        for (cs_lnum_t ij = 0; ij < 6; ij++)
          c_st_prv[c_id][ij] += w2[c_id][ij];
      }
      /* Otherwise, directly in rhs */
      else {
        for (cs_lnum_t ij = 0; ij < 6; ij++)
          rhs[c_id][ij] += w2[c_id][ij];
      }

    } /* End of loop on cells */

  } /* End for Coriolis */

  /* Wall echo terms
   * --------------- */

  if (cs_glob_turb_rans_model->irijec == 1) { // todo

    cs_array_real_fill_zero(6*n_cells, (cs_real_t*)w2);

    _rij_echo(phase_id, produc, w2);

    /* If we extrapolate the source terms: c_st_prv */
    if (st_prv_id > -1)
      cs_axpy(n_cells*6, 1., (cs_real_t *)w2, (cs_real_t *)c_st_prv);

    /* Otherwise rhs */
    else
      cs_axpy(n_cells*6, 1., (cs_real_t *)w2, (cs_real_t *)rhs);

  }

  BFT_FREE(w2);

  /* Buoyancy source term
   * -------------------- */

  if (cs_glob_turb_rans_model->has_buoyant_term == 1)
    _gravity_st_rij(f_rij, up_rhop, grav, st_prv_id, c_st_prv, fimp, rhs);

  /* Diffusion term (Daly Harlow: generalized gradient hypothesis method)
   * -------------------------------------------------------------------- */

  /* Symmetric tensor diffusivity (GGDH) */
  if (eqp->idften & CS_ANISOTROPIC_RIGHT_DIFFUSION) {

    const cs_field_t *f_a_t_visc
      = cs_field_by_name("anisotropic_turbulent_viscosity");
    const cs_real_6_t *visten = (const cs_real_6_t *)f_a_t_visc->val;

#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      for (cs_lnum_t i = 0; i < 3; i++)
        viscce[c_id][i] = visten[c_id][i] + viscl[c_id];
      for (cs_lnum_t i = 3; i < 6; i++)
        viscce[c_id][i] = visten[c_id][i];
    }

    cs_face_anisotropic_viscosity_scalar(m,
                                         fvq,
                                         viscce,
                                         eqp->verbosity,
                                         weighf,
                                         weighb,
                                         viscf,
                                         viscb);

  }

  /* Scalar diffusivity */
  else {

    cs_real_t *w1;
    BFT_MALLOC(w1, n_cells_ext, cs_real_t);

    if (eqp->idifft == 1) {
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const cs_real_t trrij = .5 * cs_math_6_trace(cvara_var[c_id]);

        const cs_real_t rctse =   crom[c_id] * cs_turb_csrij
                                * cs_math_pow2(trrij) / cvara_ep[c_id];

        w1[c_id] = viscl[c_id] + rctse;
      }
    }
    else
      cs_array_real_copy(n_cells, viscl, w1);

    cs_face_viscosity(m,
                      fvq,
                      eqp->imvisf,
                      w1,
                      viscf,
                      viscb);

    BFT_FREE(w1);

  }
}

/*----------------------------------------------------------------------------*/
/*!/
 * \brief Prepare the resolution of the segregated Reynolds stress components
 *        in the  \f$ R_{ij} - \varepsilon \f$ RANS (LRR) turbulence model.
 *
 * \param[in]   f_rij     pointer to Rij field
 * \param[in]   phase_id  turbulent phase id (-1 for single phase flow)
 * \param[in]   produc    work array for production
 * \param[in]   up_rhop   work array for \f$ \vect{u}'\rho' \f$
 * \param[in]   grav      gravity
 * \param[out]  viscf     visc*surface/dist at internal faces
 * \param[out]  viscb     visc*surface/dist at edge faces
 * \param[out]  viscce    Daly Harlow diffusion term
 * \param[out]  rhs       right hand side
 * \param[out]  fimp      Implicit term (containing unsteady term)
 * \param[out]  weighf    working array
 * \param[out]  weighb    working array
 */
/*----------------------------------------------------------------------------*/

static void
_pre_solve_lrr_sg(const cs_field_t  *f_rij,
                  int                phase_id,
                  const cs_real_t    produc[][6],
                  const cs_real_t    up_rhop[][3],
                  const cs_real_t    grav[],
                  cs_real_t          viscf[],
                  cs_real_t          viscb[],
                  cs_real_t          viscce[][6],
                  cs_real_t          rhs[][6],
                  cs_real_t          fimp[][6][6],
                  cs_real_t          weighf[][2],
                  cs_real_t          weighb[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_real_t *cell_f_vol = fvq->cell_f_vol;

  cs_field_t *f_eps = CS_F_(eps);
  cs_field_t *f_rho = CS_F_(rho);
  cs_field_t *f_mu = CS_F_(mu);

  if (phase_id >= 0) {
    f_eps = CS_FI_(eps, phase_id);
    f_rho = CS_FI_(rho, phase_id);
    f_mu = CS_FI_(mu, phase_id);
  }

  const cs_real_t *crom = f_rho->val;
  const cs_real_t *viscl = f_mu->val;
  const cs_real_t *cvara_ep = f_eps->val_pre;
  const cs_real_6_t *cvara_var = (const cs_real_6_t *)f_rij->val_pre;

  const cs_equation_param_t *eqp
    = cs_field_get_equation_param_const(f_rij);

  if (eqp->verbosity >= 1) {
    bft_printf(" Solving the variable %s\n",
               cs_field_get_label(f_rij));
  }

  cs_real_6_t *c_st_prv = nullptr;
  int kstprv = cs_field_key_id("source_term_prev_id");
  int st_prv_id = cs_field_get_key_int(f_rij, kstprv);
  if (st_prv_id > -1)
    c_st_prv = (cs_real_6_t *)cs_field_by_id(st_prv_id)->val;

  /* Time extrapolation? */
  cs_real_t  *cromo = nullptr;
  int key_t_ext_id = cs_field_key_id("time_extrapolated");
  int iroext = cs_field_get_key_int(f_rho, key_t_ext_id);
  if ((iroext > 0) && (st_prv_id > -1))
    cromo = f_rho->val_pre;
  else
    cromo = f_rho->val;

  const cs_real_t thetv = eqp->theta;

  /* Coefficient of the "Coriolis-type" term */
  const int icorio = cs_glob_physical_constants->icorio;
  const cs_turbomachinery_model_t tm_model = cs_turbomachinery_get_model();
  cs_real_t ccorio = 0;
  if (icorio) {
    if (icorio == 1)
      ccorio = 2; /* Relative velocity formulation */
    else {
      if (tm_model == CS_TURBOMACHINERY_FROZEN)
        ccorio = 1; /* Mixed relative/absolute velocity formulation */
    }
  }

  const cs_real_t d1s3 = 1./3;
  const cs_real_t d2s3 = 2./3;

  const cs_real_t crij1 = cs_turb_crij1;
  const cs_real_t crij2 = cs_turb_crij2;
  const cs_real_t crijeps = cs_turb_crij_eps;

  const auto t2v = _t2v;
  const cs_real_t vdeltij[6] = {1, 1, 1, 0, 0, 0};

  cs_lnum_t solid_stride = 1;
  int *c_is_solid_zone_flag = cs_solid_zone_flag(cs_glob_mesh);
  const int c_is_solid_ref[1] = {0};
  const int *c_is_solid = c_is_solid_zone_flag;
  if (c_is_solid == nullptr) {
    c_is_solid = c_is_solid_ref;
    solid_stride = 0;
  }

  /* Production, Pressure-Strain correlation, dissipation
   * ---------------------------------------------------- */

  /* Source term:

   *  (1-crij2) Pij (for all components of Rij)

   *  deltaij*(2/3.crij2.p+2/3.crij1.epsilon)
   *                (diagonal terms for R11, R22 et R33)

   *  -deltaij*2/3*epsilon

   * If we extrapolate the source terms
   * We modify the implicit part:
   * In phi1, we will only take rho crij1 epsilon/k and not
   *                            rho crij1 epsilon/k (1-2/3 deltaij)
   * It allows to keep  k^n instead of (R11^(n+1)+R22^n+R33^n)
   * if we extrapolate the source terms. */

  if (st_prv_id > -1) {

#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      if (c_is_solid[solid_stride*c_id])
        continue;

      /* Half-traces of Prod and R */
      const cs_real_t trprod = 0.5 * cs_math_6_trace(produc[c_id]);
      const cs_real_t trrij = 0.5 * cs_math_6_trace(cvara_var[c_id]);

      for (cs_lnum_t ij = 0; ij < 6; ij++) {
        /* Calculation of Prod+Phi1+Phi2-Eps
         *  = rhoPij - C1rho eps/k(Rij-2/3k dij)
         *           - C2rho(Pij-1/3Pkk dij) -2/3rho eps dij
         * In c_st_prv:
         *  = rhoPij-C1rho eps/k(-2/3k dij)-C2rho(Pij-1/3Pkk dij)-2/3rho eps dij
         *  = rho(2/3dij[C2 Pkk/2+(C1-1)eps)]+(1-C2)Pij */
        c_st_prv[c_id][ij] +=  cromo[c_id] * cell_f_vol[c_id]
                              * (  vdeltij[ij]*d2s3
                                 * (   crij2 * trprod
                                    + (crij1-crijeps) * cvara_ep[c_id])
                                 + (1.-crij2)*produc[c_id][ij]);

        /*  In rhs = -C1rho eps/k(Rij) = rho {-C1eps/kRij} */
        rhs[c_id][ij] +=   crom[c_id]*cell_f_vol[c_id]
                         * (-crij1 * cvara_ep[c_id] / trrij
                                   * cvara_var[c_id][ij]);

        /* Calculation of the implicit part coming from Phil = C1rho eps/k(1) */
        fimp[c_id][ij][ij] +=   crom[c_id] * cell_f_vol[c_id]
                                * crij1 * cvara_ep[c_id] / trrij * thetv;
      }

      /* If we want to implicit a part of -C1rho eps/k(-2/3k dij)
       * FIXME: check if we want to use this or if it should be removed;
       *        previously "isoluc = 2", never called.
       */

      if (false) {
        for (cs_lnum_t i = 0; i < 3; i++) {
          cs_real_t t1 = d1s3 * crij1 * cvara_ep[c_id] / trrij;

          /*  We substract cromo = -C1rho eps/k(-1/3Rij dij) */
          c_st_prv[c_id][i] -=    cromo[c_id] * cell_f_vol[c_id]
                                 * t1 * cvara_var[c_id][i];

          /* We add to rhs = -C1rho eps/k(-1/3Rij) */
          rhs[c_id][i] +=   crom[c_id] * cell_f_vol[c_id]
                           * t1 * cvara_var[c_id][i];

          /* We add to fimp = C1rho eps/k(-1/3) */
          fimp[c_id][i][i] += crom[c_id] * cell_f_vol[c_id] * t1;
        }
      }
    }

  }

  /* If we do not extrapolate the source terms */
  else {

#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      if (c_is_solid[solid_stride*c_id])
        continue;

      /* Half-traces of Prod and R */
      const cs_real_t trprod = 0.5 * cs_math_6_trace(produc[c_id]);
      const cs_real_t trrij = 0.5 * cs_math_6_trace(cvara_var[c_id]);

      for (cs_lnum_t ij = 0; ij < 6; ij++) {
        /* Calculation of Prod+Phi1+Phi2-Eps
         *  = rhoPij - C1rho eps/k(Rij-2/3k dij)
         *           - C2rho(Pij-1/3Pkk dij) -2/3rho eps dij
         *  = rho{2/3dij[C2 Pkk/2+(C1-1)eps)]+(1-C2)Pij */
        rhs[c_id][ij] +=  crom[c_id] * cell_f_vol[c_id]
                        * (  vdeltij[ij] * d2s3
                           * (   crij2 * trprod
                              + (crij1-crijeps) * cvara_ep[c_id])
                           + (1-cs_turb_crij2) * produc[c_id][ij]
                           - crij1 * cvara_ep[c_id]/trrij * cvara_var[c_id][ij]);

        /* Calculation of the implicit part coming from Phi1
         *  = C1rho eps/k(1-1/3 dij) */
        fimp[c_id][ij][ij] +=   crom[c_id] * cell_f_vol[c_id]
                                * (1-d1s3 *vdeltij[ij]) * crij1
                                * cvara_ep[c_id]/trrij;
      }
    }

  }

  BFT_FREE(c_is_solid_zone_flag);

  /* Coriolis terms in the Phi1 and production
   * -----------------------------------------*/

  cs_real_6_t *w2;
  BFT_MALLOC(w2, n_cells_ext, cs_real_6_t);

  if (icorio == 1 || tm_model == CS_TURBOMACHINERY_FROZEN) {
    const int *irotce = cs_turbomachinery_get_cell_rotor_num();

#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      cs_real_t matrot[3][3];

      int rot_id = icorio;
      if (tm_model == CS_TURBOMACHINERY_FROZEN) {
        rot_id = irotce[c_id];
        if (rot_id < 1)
          continue;
      }

      const cs_rotation_t *r = cs_glob_rotation + rot_id;
      cs_rotation_coriolis_t(r, 1., matrot);

      /* Compute Gij: (i,j) component of the Coriolis production */
      for (cs_lnum_t ij = 0; ij < 6; ij++) {
        cs_lnum_t i = _iv2t[ij];
        cs_lnum_t j = _jv2t[ij];

        w2[c_id][ij] = 0.;
        for (cs_lnum_t kk = 0; kk < 3; kk++)
          w2[c_id][ij] -=   ccorio
                           * (  matrot[kk][i]*cvara_var[c_id][t2v[kk][j]]
                              + matrot[kk][j]*cvara_var[c_id][t2v[kk][i]]);
      }

      /* Coriolis contribution in the Phi1 term: (1-C2/2)Gij */
      if (icorio == 1) {
        for (cs_lnum_t ij = 0; ij < 6; ij++)
          w2[c_id][ij] *= crom[c_id] * cell_f_vol[c_id ]* (1.-0.5*cs_turb_crij2);
      }

      /* If source terms are extrapolated */
      if (st_prv_id > -1) {
        for (cs_lnum_t ij = 0; ij < 6; ij++)
          c_st_prv[c_id][ij] += w2[c_id][ij];
      }
      /* Otherwise, directly in RHS */
      else {
        for (cs_lnum_t ij = 0; ij < 6; ij++)
          rhs[c_id][ij] += w2[c_id][ij];
      }

    } /* End of loop on cells */

  }

  /* Wall echo terms
   * --------------- */

  if (cs_glob_turb_rans_model->irijec == 1) { // todo

    cs_array_real_fill_zero(6*n_cells, (cs_real_t*)w2);

    _rij_echo(phase_id, produc, w2);

    /* If we extrapolate the source terms: c_st_prv */
    if (st_prv_id > -1)
      cs_axpy(n_cells*6, 1., (cs_real_t *)w2, (cs_real_t *)c_st_prv);

    /* Otherwise rhs */
    else
      cs_axpy(n_cells*6, 1., (cs_real_t *)w2, (cs_real_t *)rhs);

  }

  BFT_FREE(w2);

  /* Buoyancy source term
   * -------------------- */

  if (cs_glob_turb_rans_model->has_buoyant_term == 1)
    _gravity_st_rij(f_rij, up_rhop, grav, st_prv_id, c_st_prv, fimp, rhs);

  /* Diffusion term (Daly Harlow: generalized gradient hypothesis method)
   * -------------------------------------------------------------------- */

  /* Symmetric tensor diffusivity (GGDH) */
  if (eqp->idften & CS_ANISOTROPIC_RIGHT_DIFFUSION) {

    const cs_field_t *f_a_t_visc
      = cs_field_by_name("anisotropic_turbulent_viscosity");
    const cs_real_6_t *visten = (const cs_real_6_t *)f_a_t_visc->val;

#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      for (cs_lnum_t i = 0; i < 3; i++)
        viscce[c_id][i] = visten[c_id][i] + viscl[c_id];
      for (cs_lnum_t i = 3; i < 6; i++)
        viscce[c_id][i] = visten[c_id][i];
    }

    cs_face_anisotropic_viscosity_scalar(m,
                                         fvq,
                                         viscce,
                                         eqp->verbosity,
                                         weighf,
                                         weighb,
                                         viscf,
                                         viscb);

  }

  /* Scalar diffusivity */
  else {

    cs_real_t *w1;
    BFT_MALLOC(w1, n_cells_ext, cs_real_t);

    if (eqp->idifft == 1) {
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        const cs_real_t trrij = .5 * cs_math_6_trace(cvara_var[c_id]);

        const cs_real_t rctse =   crom[c_id] * cs_turb_csrij
                                * cs_math_pow2(trrij) / cvara_ep[c_id];

        w1[c_id] = viscl[c_id] + rctse;
      }
    }
    else
      cs_array_real_copy(n_cells, viscl, w1);

    cs_face_viscosity(m,
                      fvq,
                      eqp->imvisf,
                      w1,
                      viscf,
                      viscb);

    BFT_FREE(w1);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prepare the resolution ofthe coupled Reynolds stress components
 *        in the \f$ R_{ij} - \varepsilon \f$ RANS (SSG) turbulence model.
 *
 * \param[in]   f_rij     pointer to Rij field
 * \param[in]   phase_id  turbulent phase id (-1 for single phase flow)
 * \param[in]   gradv     work array for the velocity grad term
 * \param[in]   produc    work array for production
 * \param[in]   up_rhop   work array for \f$ \vect{u}'\rho' \f$
 * \param[in]   grav      gravity
 * \param[out]  viscf     visc*surface/dist at internal faces
 * \param[out]  viscb     visc*surface/dist at edge faces
 * \param[out]  viscce    Daly Harlow diffusion term
 * \param[out]  rhs       right hand side
 * \param[out]  fimp      Implicit term (containing unsteady term)
 * \param[out]  weighf    working array
 * \param[out]  weighb    working array
 */
/*----------------------------------------------------------------------------*/

static void
_pre_solve_ssg(const cs_field_t  *f_rij,
               int                phase_id,
               const cs_real_t    gradv[][3][3],
               const cs_real_t    produc[][6],
               const cs_real_t    up_rhop[][3],
               const cs_real_t    grav[],
               cs_real_t          viscf[],
               cs_real_t          viscb[],
               cs_real_t          viscce[][6],
               cs_real_t          rhs[][6],
               cs_real_t          fimp[][6][6],
               cs_real_t          weighf[][2],
               cs_real_t          weighb[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_real_t *cell_f_vol = fvq->cell_f_vol;

  cs_dispatch_context ctx;

  cs_field_t *f_eps = CS_F_(eps);
  cs_field_t *f_alpbl = CS_F_(alp_bl);
  cs_field_t *f_rho = CS_F_(rho);
  cs_field_t *f_mu = CS_F_(mu);
  cs_field_t *f_mut = CS_F_(mu_t);

  if (phase_id >= 0) {
    f_eps = CS_FI_(eps, phase_id);
    f_alpbl = CS_FI_(alp_bl, phase_id);
    f_rho = CS_FI_(rho, phase_id);
    f_mu = CS_FI_(mu, phase_id);
    f_mut = CS_FI_(mu_t, phase_id);
  }

  const cs_real_t *crom = f_rho->val;
  const cs_real_t *viscl = f_mu->val;
  const cs_real_t *visct = f_mut->val;
  const cs_real_t *cvara_ep = f_eps->val_pre;
  const cs_real_6_t *cvara_var = (const cs_real_6_t *)f_rij->val_pre;

  cs_real_t *cvar_al = nullptr;
  if (cs_glob_turb_model->iturb != CS_TURB_RIJ_EPSILON_SSG)
    cvar_al = f_alpbl->val;

  const cs_equation_param_t *eqp
    = cs_field_get_equation_param_const(f_rij);

  if (eqp->verbosity >= 1) {
    bft_printf(" Solving the variable %s\n ",
               cs_field_get_label(f_rij));
  }

  const int coupled_components = cs_glob_turb_rans_model->irijco;

  cs_real_6_t *c_st_prv = nullptr;
  int kstprv = cs_field_key_id("source_term_prev_id");
  int st_prv_id = cs_field_get_key_int(f_rij, kstprv);
  if (st_prv_id > -1)
    c_st_prv = (cs_real_6_t *)cs_field_by_id(st_prv_id)->val;

  /* Time extrapolation ? */
  cs_real_t *cromo = nullptr;
  int key_t_ext_id = cs_field_key_id("time_extrapolated");
  int iroext = cs_field_get_key_int(f_rho, key_t_ext_id);
  if ((st_prv_id > -1) && (iroext > 0))
    cromo = f_rho->val_pre;
  else
    cromo = f_rho->val;

  /* Coefficient of the "Coriolis-type" term */
  const int icorio = cs_glob_physical_constants->icorio;
  cs_turbomachinery_model_t iturbo = cs_turbomachinery_get_model();
  cs_real_t ccorio = 0;
  if (icorio == 1)
    ccorio = 2; // Relative velocity formulation
  else if (iturbo == CS_TURBOMACHINERY_FROZEN)
    ccorio = 1;

  const cs_real_t d1s2 = 0.5;
  const cs_real_t d1s3 = 1./3;
  const cs_real_t d2s3 = 2./3;

  const cs_real_t cebmr1 = cs_turb_cebmr1;
  const cs_real_t cebmr2 = cs_turb_cebmr2;
  const cs_real_t cebmr3 = cs_turb_cebmr3;
  const cs_real_t cebmr4 = cs_turb_cebmr4;
  const cs_real_t cebmr5 = cs_turb_cebmr5;

  const cs_real_t cebms1 = cs_turb_cebms1;

  const cs_real_t cmu = cs_turb_cmu;

  const cs_real_t csrij  = cs_turb_csrij;
  const cs_real_t crijeps = cs_turb_crij_eps;

  const cs_real_t cssgr1 = cs_turb_cssgr1;
  const cs_real_t cssgr2 = cs_turb_cssgr2;
  const cs_real_t cssgr3 = cs_turb_cssgr3;
  const cs_real_t cssgr4 = cs_turb_cssgr4;
  const cs_real_t cssgr5 = cs_turb_cssgr5;

  const cs_real_t cssgs1 = cs_turb_cssgs1;
  const cs_real_t cssgs2 = cs_turb_cssgs2;

  const cs_turb_model_type_t iturb
    = (cs_turb_model_type_t)cs_glob_turb_model->iturb;

  const auto t2v = _t2v;
  const auto iv2t = _iv2t;
  const auto jv2t = _jv2t;
  const cs_real_t tdeltij[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
  const cs_real_t vdeltij[6] = {1, 1, 1, 0, 0, 0};

  /* Production, Pressure-Strain correlation, dissipation, Coriolis
   * -------------------------------------------------------------- */

  /* Source term
   *  -rho*epsilon*( Cs1*aij + Cs2*(aikajk -1/3*aijaij*deltaij))
   *  -Cr1*P*aij + Cr2*rho*k*sij - Cr3*rho*k*sij*sqrt(aijaij)
   *  +Cr4*rho*k(aik*sjk+ajk*sik-2/3*akl*skl*deltaij)
   *  +Cr5*rho*k*(aik*rjk + ajk*rik)
   *  -2/3*epsilon*deltaij */

  cs_real_3_t *grad_al = nullptr;
  if (cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM) {
    CS_MALLOC_HD(grad_al, n_cells_ext, cs_real_3_t, cs_alloc_mode);
    cs_field_gradient_scalar(f_alpbl, true, 1, grad_al);
  }

  const int *irotce = cs_turbomachinery_get_cell_rotor_num();

  cs_real_t *w1, *w2;
  CS_MALLOC_HD(w1, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(w2, n_cells_ext, cs_real_t, cs_alloc_mode);

  cs_lnum_t solid_stride = 1;
  int *c_is_solid_zone_flag = cs_solid_zone_flag(cs_glob_mesh);
  const int c_is_solid_ref[1] = {0};
  if (c_is_solid_zone_flag == nullptr) {
    if (cs_alloc_mode > CS_ALLOC_HOST) {
      CS_MALLOC_HD(c_is_solid_zone_flag, 1, int, cs_alloc_mode);
      c_is_solid_zone_flag[0] = 0;
    }
    solid_stride = 0;
  }
  const int *c_is_solid = c_is_solid_zone_flag;
  if (c_is_solid == nullptr)
    c_is_solid = c_is_solid_ref;

  const cs_rotation_t *rotation = cs_glob_rotation;

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    const cs_lnum_t ind = solid_stride*c_id;
    if (c_is_solid[ind])
      return; // return from lambda function == continue in loop

    cs_real_t xnal[3] = {0, 0, 0};

    cs_real_t matrot[3][3];
    cs_real_t impl_drsm[6][6];

    /* EBRSM: compute the magnitude of the Alpha gradient */

    if (iturb == CS_TURB_RIJ_EPSILON_EBRSM) {
      cs_math_3_normalize(grad_al[c_id], xnal);
    }

    cs_real_t xrij[3][3], xprod[3][3];
    cs_real_t xaniso[3][3], xstrai[3][3], xrotac[3][3];

    const cs_real_t trrij  = 0.5 * cs_math_6_trace(cvara_var[c_id]);

    /* Rij */

    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        xrij[i][j] = cvara_var[c_id][t2v[i][j]];
      }
    }

    /* Pij */
    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        xprod[i][j] = produc[c_id][t2v[i][j]];
      }
    }

    /* aij */
    for (cs_lnum_t i = 0; i < 3; i++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        xaniso[i][j] =   xrij[i][j] / trrij
                       - d2s3 * tdeltij[i][j];
      }
    }

    /* Sij */
    xstrai[0][0] = gradv[c_id][0][0];
    xstrai[0][1] = d1s2 * (gradv[c_id][0][1] + gradv[c_id][1][0]);
    xstrai[0][2] = d1s2 * (gradv[c_id][0][2] + gradv[c_id][2][0]);
    xstrai[1][0] = xstrai[0][1];
    xstrai[1][1] = gradv[c_id][1][1];
    xstrai[1][2] = d1s2 * (gradv[c_id][1][2] + gradv[c_id][2][1]);
    xstrai[2][0] = xstrai[0][2];
    xstrai[2][1] = xstrai[1][2];
    xstrai[2][2] = gradv[c_id][2][2];

    /* omegaij */
    xrotac[0][0] = 0;
    xrotac[0][1] = d1s2 * (gradv[c_id][0][1] - gradv[c_id][1][0]);
    xrotac[0][2] = d1s2 * (gradv[c_id][0][2] - gradv[c_id][2][0]);
    xrotac[1][0] = -xrotac[0][1];
    xrotac[1][1] = 0;
    xrotac[1][2] = d1s2 * (gradv[c_id][1][2] - gradv[c_id][2][1]);
    xrotac[2][0] = -xrotac[0][2];
    xrotac[2][1] = -xrotac[1][2];
    xrotac[2][2] = 0;

    int rot_id = icorio;
    if (iturbo == CS_TURBOMACHINERY_FROZEN)
      rot_id = irotce[c_id];

    /* Rotating frame of reference => "absolute" vorticity */

    if (rot_id >= 1) {
      const cs_rotation_t *r = rotation + rot_id;

      cs_rotation_coriolis_t(r, 1., matrot);
      for (cs_lnum_t i = 0; i < 3; i++) {
        for (cs_lnum_t j = i; j < 3; j++) {
          for (cs_lnum_t k = 0; k < 3; k++)
              xprod[j][i] -= ccorio * (  matrot[k][i] * xrij[k][j]
                                       + matrot[k][j] * xrij[k][i]);
        }
      }

      xprod[0][1] = xprod[1][0]; /* Ensure symmetry (j<i not computed) */
      xprod[0][2] = xprod[2][0];
      xprod[1][2] = xprod[2][1];

    }

    const cs_real_t trprod = 0.5 * (xprod[0][0] + xprod[1][1] + xprod[2][2]);

    /* aii = aijaij */
    cs_real_t aii = 0, aklskl = 0;

    for (cs_lnum_t j = 0; j < 3; j++) {
      for (cs_lnum_t i = 0; i < 3; i++) {
        aii += cs_math_pow2(xaniso[j][i]);       /* aij.aij */
        aklskl += xaniso[j][i] * xstrai[j][i];   /* aij.Sij */
      }
    }

    if (coupled_components != 0) {

      /* Computation of implicit components */
      cs_real_t sym_strain[6] = {xstrai[0][0],
                                 xstrai[1][1],
                                 xstrai[2][2],
                                 xstrai[1][0],
                                 xstrai[2][1],
                                 xstrai[2][0]};

      /* Compute inverse matrix of R^n
         (scaling by tr(R) for numerical stability) */
      cs_real_t matrn[6];
      for (cs_lnum_t ij = 0; ij < 6; ij++)
        matrn[ij] = cvara_var[c_id][ij] / trrij;

      cs_real_t oo_matrn[6];
      cs_math_sym_33_inv_cramer(matrn, oo_matrn);
      for (cs_lnum_t ij = 0; ij < 6; ij++)
        oo_matrn[ij] /= trrij;

      cs_real_t impl_lin_cst = 0;
      cs_real_t impl_id_cst = 0;

      /* Compute the maximal eigenvalue (in terms of norm!) of S */
      cs_real_t eigen_vals[3];
      _sym_33_eigen(sym_strain, eigen_vals);
      cs_real_t eigen_max = cs_math_fabs(eigen_vals[0]);
      for (cs_lnum_t i = 1; i < 3; i++)
        eigen_max = cs_math_fmax(cs_math_fabs(eigen_max),
                                 cs_math_fabs(eigen_vals[i]));

      /* Constant for the dissipation */
      const cs_real_t ceps_impl = d1s3 * cvara_ep[c_id];

      if (iturb == CS_TURB_RIJ_EPSILON_SSG) {

        /* Identity constant for phi3 */
        const cs_real_t cphi3impl = cs_math_fabs(cssgr2 - cssgr3*sqrt(aii));

        /* Identity constant */
        impl_id_cst = - d2s3 * cssgr1 * cs_math_fmin(trprod, 0)
                      - d1s3 * cssgs2 * cvara_ep[c_id] * aii
                      + cphi3impl * trrij * eigen_max
                      + 2. * d2s3 * cssgr4 * trrij * eigen_max
                      + d2s3 * trrij * cssgr4 * cs_math_fmax(aklskl, 0);

        /* Linear constant */
        impl_lin_cst = eigen_max * (crijeps + cssgr4 + cssgr5);

        cs_real_t implmat2add[3][3];
        for (cs_lnum_t i = 0; i < 3; i++) {
          for (cs_lnum_t j = 0; j < 3; j++) {
            const cs_lnum_t _ij = t2v[i][j];
            implmat2add[i][j] =   xrotac[i][j]
                                + impl_lin_cst * vdeltij[_ij]
                                + impl_id_cst * d1s2 * oo_matrn[_ij]
                                + ceps_impl * oo_matrn[_ij];
          }
        }
        /* Compute the 6x6 matrix A which verifies
         * A.R = M.R + R.M^t */
        cs_math_reduce_sym_prod_33_to_66(implmat2add, impl_drsm);

      }
      else { /* iturb == CS_TURB_RIJ_EPSILON_EBRSM */

        const cs_real_t alpha3 = cs_math_pow3(cvar_al[c_id]);

        /* Phi3 constant */
        const cs_real_t cphi3impl = cs_math_fabs(cebmr2 - cebmr3*sqrt(aii));

        /* PhiWall + epsilon_wall constants for EBRSM */
        const cs_real_t cphiw_impl = 6 * (1 - alpha3)
          * crijeps * cvara_ep[c_id] / trrij;

        /* The implicit components of Phi (pressure-velocity fluctuations)
         * are split into the linear part (A*R) and Id part (A*Id). */

        /* Identity constant */
        impl_id_cst =   alpha3*(-d2s3 * cebmr1 * cs_math_fmin(trprod, 0)
                      + cphi3impl * trrij * eigen_max
                      + 2 * d2s3 * cebmr4 * trrij * eigen_max
                      + d2s3 * trrij * cebmr4 * cs_math_fmax(aklskl, 0));

        /* Linear constant */
        impl_lin_cst
          = eigen_max * (1 + cebmr4*alpha3 + cebmr5*alpha3) + cphiw_impl;

        cs_real_t implmat2add[3][3];
        for (cs_lnum_t i = 0; i < 3; i++) {
          for (cs_lnum_t j = 0; j < 3; j++) {
            const cs_lnum_t _ij = t2v[i][j];
            implmat2add[i][j] =   xrotac[i][j]
                                + impl_lin_cst * vdeltij[_ij]
                                + impl_id_cst * d1s2 * oo_matrn[_ij]
                                + alpha3 * ceps_impl * oo_matrn[_ij];
          }
        }
        /* Compute the 6x6 matrix A which verifies
         * A.R = M.R + R.M^t */
        cs_math_reduce_sym_prod_33_to_66(implmat2add, impl_drsm);
      }

    } /* end if irijco != 0 */

    /* Rotating frame of reference => "absolute" vorticity */

    if (icorio == 1) {
      for (cs_lnum_t i = 0; i < 3; i++)
        for (cs_lnum_t j = 0; j < 3; j++) {
          xrotac[i][j] -= matrot[i][j];
      }
    }

    for (cs_lnum_t ij = 0; ij < 6; ij++) {
      cs_lnum_t i = iv2t[ij];
      cs_lnum_t j = jv2t[ij];

      cs_real_t aiksjk = 0, aikrjk = 0, aikakj = 0;

      for (cs_lnum_t k = 0; k < 3; k++) {
        // aiksjk = aik.Sjk+ajk.Sik
        aiksjk +=   xaniso[i][k] * xstrai[j][k]
                  + xaniso[j][k] * xstrai[i][k];
        // aikrjk = aik.Omega_jk + ajk.omega_ik
        aikrjk +=   xaniso[i][k] * xrotac[j][k]
                  + xaniso[j][k] * xrotac[i][k];
        // aikakj = aik*akj
        aikakj += xaniso[i][k] * xaniso[j][k];
      }

      /* If we extrapolate the source terms (rarely), we put everything
       * in the previous ST.
       * We do not implicit the term with Cs1*aij or Cr1*P*aij.
       * Otherwise, we put all in rhs and we can implicit Cs1*aij
       * and Cr1*P*aij. Here we store the right-hand-side and the
       * implicit term in w1 and w2, to avoid the test (st_prv_id >= 0)
       * in the loop on cells.
       * In the term with w1, which is set to be extrapolated, we use cromo.
       * The implicitation of the two terms can also be done in the case of
       * extrapolation, by isolating those two terms and by putting it in
       * the RHS but not in the prev. ST and by using ipcrom ....
       * to be modified if needed. */

      if (iturb == CS_TURB_RIJ_EPSILON_SSG) {

        /* Explicit terms */
        const cs_real_t pij =     xprod[j][i];
        const cs_real_t phiij1 =   -cvara_ep[c_id]
                                 * (  cssgs1 * xaniso[j][i]
                                    + cssgs2 * (aikakj - d1s3*vdeltij[ij]*aii));
        const cs_real_t phiij2 =   -cssgr1 * trprod * xaniso[j][i]
                                 +   trrij * xstrai[j][i]
                                   * (cssgr2 - cssgr3*sqrt(aii))
                                 + cssgr4*trrij * (  aiksjk
                                                   - d2s3 * vdeltij[ij] * aklskl)
                                 + cssgr5*trrij * aikrjk;
        const cs_real_t epsij = -d2s3*crijeps * cvara_ep[c_id] * vdeltij[ij];

        w1[c_id] =   cromo[c_id] * cell_f_vol[c_id]
                   * (pij + phiij1 + phiij2 + epsij);

        /* Implicit terms */
        w2[c_id] =  crom[c_id] * cell_f_vol[c_id] / trrij
                  * (  cssgs1 * cvara_ep[c_id]
                     + cssgr1 * cs_math_fmax(trprod, 0));

      }
      else { /* cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM */

        /* Compute the explicit term
         * Compute the terms near the walls and almost homogeneous
         * to phi and epsilon.
         * Compute the term near the wall \f$ \Phi_{ij}^w \f$ --> w3
         *
         * Phiw = -5.0 * (eps/k) * [R*Xn + Xn^T*R - 0.5*tr(Xn*R)*(Xn+Id)] */
        const cs_real_t xnnd = d1s2 * (xnal[i]*xnal[j] + vdeltij[ij]);

        cs_real_t phiijw = 0;
        for (cs_lnum_t kk = 0; kk < 3; kk++) {
          phiijw += xrij[kk][i] * xnal[j] * xnal[kk];
          phiijw += xrij[kk][j] * xnal[i] * xnal[kk];
          for (cs_lnum_t ll = 0; ll < 3; ll++)
            phiijw -= xrij[ll][kk] * xnal[kk] * xnal[ll] * xnnd;
        }
        phiijw = -5. * cvara_ep[c_id] / trrij * phiijw;

        /* Compute the almost homogeneous term \f$ \phi_{ij}^h \f$ */
        const cs_real_t phiij1 = -cvara_ep[c_id] * cebms1 * xaniso[j][i];
        const cs_real_t phiij2
          =   -cebmr1 * trprod * xaniso[j][i]
            +   trrij * (  xstrai[j][i] * (cebmr2 - cebmr3 * sqrt(aii))
                         + cebmr4 * (aiksjk - d2s3 * vdeltij[ij] * aklskl)
                         + cebmr5 * aikrjk);

        /* Compute \f $\e_{ij}^w \f$ (Rotta model)
         * Rij/k*epsilon */
        const cs_real_t epsijw = xrij[j][i] / trrij * cvara_ep[c_id];

        /* Compute \e_{ij}^h */
        const cs_real_t epsij = d2s3 * cvara_ep[c_id] * vdeltij[ij];

        /* Compute explicit ST of the Rij equation
         *  \f[ P_{ij} + (1-\alpha^3)\Phi_{ij}^w + \alpha^3\Phi_{ij}^h
         * - (1-\alpha^3)\e_{ij}^w   - \alpha^3\e_{ij}^h  ]\f$ --> W1 */
        const cs_real_t  alpha3 = cs_math_pow3(cvar_al[c_id]);

        w1[c_id] =    crom[c_id] * cell_f_vol[c_id]
                   * (  xprod[j][i]
                      + (1-alpha3) * phiijw + alpha3 * (phiij1+phiij2)
                      - (1-alpha3) * epsijw - alpha3 * epsij);

        /* Implicit terms */
        cs_real_t  epsijw_imp = 0; // FIXME
        if (coupled_components == 0)
          epsijw_imp = 6. * (1.-alpha3)*crijeps * cvara_ep[c_id] / trrij;

        /* The term below corresponds to the implicit part of SSG
         * in the context of elliptical weighting, it is multiplied by
         * \f$ \alpha^3 \f$*/
        w2[c_id] = crom[c_id] * cell_f_vol[c_id]
                   * (  cebms1 * cvara_ep[c_id] / trrij * alpha3
                      + cebmr1 * cs_math_fmax(trprod/trrij, 0) * alpha3
                      + epsijw_imp);

      } /* End of test on turbulence model */

      if (st_prv_id > -1) {
        c_st_prv[c_id][ij] += w1[c_id];
      }
      else {
        rhs[c_id][ij] += w1[c_id];
        fimp[c_id][ij][ij] += w2[c_id];

        if (coupled_components != 0) {
          for (cs_lnum_t jj = 0; jj < 6; jj++)
            fimp[c_id][ij][jj] +=   crom[c_id] * cell_f_vol[c_id]
                                    * impl_drsm[ij][jj];
        }
      }
    } /* End of loop on ij */

  }); /* end loop on cells */
  ctx.wait();

  CS_FREE_HD(grad_al);
  CS_FREE_HD(w2);

  /* Buoyancy source term
   * -------------------- */

  if (cs_glob_turb_rans_model->has_buoyant_term == 1)
    _gravity_st_rij(f_rij, up_rhop, grav, st_prv_id, c_st_prv, fimp, rhs);

  CS_FREE_HD(c_is_solid_zone_flag);

  /* Diffusion term (Daly Harlow: generalized gradient hypothesis method)
   * -------------------------------------------------------------------- */

  /* Symmetric tensor diffusivity (GGDH) */
  if (eqp->idften & CS_ANISOTROPIC_RIGHT_DIFFUSION) {

    const cs_field_t *f_a_t_visc
      = cs_field_by_name("anisotropic_turbulent_viscosity");
    const cs_real_6_t *visten = (const cs_real_6_t *)f_a_t_visc->val;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t i = 0; i < 3; i++)
        viscce[c_id][i] = visten[c_id][i] + viscl[c_id];
      for (cs_lnum_t i = 3; i < 6; i++)
        viscce[c_id][i] = visten[c_id][i];
    });
    ctx.wait();

    cs_face_anisotropic_viscosity_scalar(m,
                                         fvq,
                                         viscce,
                                         eqp->verbosity,
                                         weighf,
                                         weighb,
                                         viscf,
                                         viscb);

  }

  /* Scalar diffusivity */
 else {

    if (eqp->idifft == 1) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        w1[c_id] = viscl[c_id] + (csrij * visct[c_id] / cmu);
      });
      ctx.wait();
    }
    else
      cs_array_copy<cs_real_t>(n_cells, viscl, w1);

    cs_face_viscosity(m,
                      fvq,
                      eqp->imvisf,
                      w1,
                      viscf,
                      viscb);
  }

  CS_FREE_HD(w1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve epsilon for the \f$ R_{ij} - \varepsilon \f$ RANS
 *        turbulence model.
 *
 * \param[in]   phase_id  turbulent phase id (-1 for single phase flow)
 * \param[in]   gradv     work array for the term grad
 *                        of velocity only for iturb=31
 * \param[in]   produc    work array for production (without
 *                        rho volume) only for iturb=30
 * \param[in]   up_rhop   work array for \f$ \vect{u}'\rho' \f$
 *                        source terms or mass rate
 * \param[in]   grav      gravity
 * \param[in]   viscf     visc*surface/dist at internal faces
 * \param[in]   viscb     visc*surface/dist at edge faces
 * \param[out]  rhs       right hand side
 * \param[out]  fimp      Implicit term (containing unsteady term)
 !*/
/*----------------------------------------------------------------------------*/

static void
_solve_epsilon(int              phase_id,
               const cs_real_t  gradv[][3][3],
               const cs_real_t  produc[][6],
               const cs_real_t  up_rhop[][3],
               const cs_real_t  grav[],
               cs_real_t        viscf[],
               cs_real_t        viscb[],
               cs_real_t        rhs[],
               cs_real_t        fimp[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  cs_dispatch_context ctx;

  const cs_real_t *cell_f_vol = fvq->cell_f_vol;

  cs_field_t *f_rij = CS_F_(rij);
  cs_field_t *f_eps = CS_F_(eps);
  cs_field_t *f_alpbl = CS_F_(alp_bl);
  cs_field_t *f_vel = CS_F_(vel);
  cs_field_t *f_rho = CS_F_(rho);
  cs_field_t *f_mu = CS_F_(mu);
  cs_field_t *f_mut = CS_F_(mu_t);

  if (phase_id >= 0) {
    f_rij = CS_FI_(rij, phase_id);
    f_eps = CS_FI_(eps, phase_id);
    f_alpbl = CS_FI_(alp_bl, phase_id);
    f_vel = CS_FI_(vel, phase_id);
    f_rho = CS_FI_(rho, phase_id);
    f_mu = CS_FI_(mu, phase_id);
    f_mut = CS_FI_(mu_t, phase_id);
  }

  const cs_real_t *dt = CS_F_(dt)->val;
  const cs_real_t *crom = f_rho->val;
  const cs_real_t *viscl = f_mu->val;
  const cs_real_t *visct = f_mut->val;
  const cs_real_t *cvara_ep = f_eps->val_pre;
  const cs_real_6_t *cvara_rij = (const cs_real_6_t *)f_rij->val_pre;

  const cs_real_t *cvar_al = nullptr;
  if (cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM)
    cvar_al = (const cs_real_t *)(f_alpbl->val);

  cs_real_t *cvar_ep = f_eps->val;

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  int iflmas =  cs_field_get_key_int(f_vel, kimasf);
  int iflmab =  cs_field_get_key_int(f_vel, kbmasf);

  const cs_real_t *imasfl = cs_field_by_id(iflmas)->val;
  const cs_real_t *bmasfl = cs_field_by_id(iflmab)->val;

  const int ksigmas = cs_field_key_id("turbulent_schmidt");
  const cs_real_t sigmae = cs_field_get_key_double(f_eps, ksigmas);

  const cs_equation_param_t *eqp
    = cs_field_get_equation_param(f_eps);

  if (eqp->verbosity >= 1) {
    bft_printf(" Solving the variable %s\n",
               cs_field_get_label(f_eps));
  }

  int kstprv = cs_field_key_id("source_term_prev_id");
  int st_prv_id = cs_field_get_key_int(f_eps, kstprv);
  cs_real_t *c_st_prv = nullptr, *cromo = nullptr;
  if (st_prv_id > -1)
    c_st_prv = cs_field_by_id(st_prv_id)->val;

  int key_t_ext_id = cs_field_key_id("time_extrapolated");
  int iroext = cs_field_get_key_int(f_rho, key_t_ext_id);
  if ((iroext > 0) && (st_prv_id > -1))
    cromo = f_rho->val_pre;
  else
    cromo = f_rho->val;

  const cs_real_t ce1 = cs_turb_ce1;
  const cs_real_t ceps2 = cs_turb_ce2;
  const cs_real_t xct = cs_turb_xct;

  /* S as Source, V as Variable */
  const cs_time_scheme_t *time_scheme = cs_glob_time_scheme;
  const cs_real_t thets = time_scheme->thetst;
  const cs_real_t thetv = eqp->theta;

  cs_arrays_set_value<cs_real_t, 1>(n_cells, 0., rhs, fimp);

  /* Work arrays */
  cs_real_t *w1;
  CS_MALLOC_HD(w1, n_cells_ext, cs_real_t, cs_alloc_mode);

  /* User source terms
   * ----------------- */

  cs_user_source_terms(cs_glob_domain,
                       f_eps->id,
                       rhs,
                       fimp);

  if (cs_glob_porous_model == 3)
    cs_immersed_boundary_wall_functions(f_eps->id, rhs, fimp);

  /* If we extrapolate the source terms */
  if (st_prv_id > -1) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      /* Save for exchange */
      const cs_real_t tuexpe = c_st_prv[c_id];
      /* For the continuation and the next time step */
      c_st_prv[c_id] = rhs[c_id];
      /* RHS of previous time step
       *   We assume -fimp > 0: we implicit
       *   the user source term (the rest)  */
      rhs[c_id] = fimp[c_id]*cvara_ep[c_id] - thets*tuexpe;
      /* Diagonal */
      fimp[c_id] = -thetv*fimp[c_id];
    });
  }
  else {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      rhs[c_id]  += fimp[c_id]*cvara_ep[c_id];
      fimp[c_id]  = cs_math_fmax(-fimp[c_id], cs_math_zero_threshold);
    });
  }

  /* Lagrangian source terms
   * ----------------------- */

  /* Second order is not taken into account  */

  if (   cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING
      && (cs_glob_lagr_source_terms->ltsdyn == 1)) {

    const cs_real_t ce4 = cs_turb_ce4;

    const cs_real_6_t *lagr_st_rij
      = (const cs_real_6_t *)cs_field_by_name("lagr_st_rij")->val;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      /* Source terms with epsilon */
      const cs_real_t st_eps
        = -0.5 * cell_f_vol[c_id] * cs_math_6_trace(lagr_st_rij[c_id]);

      /* k */
      const cs_real_t k = 0.5 * cs_math_6_trace(cvara_rij[c_id]);

      /* equiv:       cs_turb_ce4 * st_eps * eps / (k / eps) */
      rhs[c_id]  += ce4 * st_eps / k;

      /* equiv:                    -cs_turb_ce4 * st_eps * / (k/eps) */
      fimp[c_id] += cs_math_fmax(-ce4 * st_eps / k * cvara_ep[c_id],
                                 cs_math_zero_threshold);
    });
  }
  ctx.wait();

  /* Mass source term
   * ---------------- */

  if (eqp->n_volume_mass_injections > 0) {

    int *mst_type = nullptr;
    cs_lnum_t n_elts = 0;
    const cs_lnum_t *elt_ids = nullptr;
    cs_real_t *mst_val = nullptr, *mst_val_p = nullptr;

    /* If we extrapolate the source terms, we put Gamma Pinj in c_st_prv;
        Otherwise we put it directly in rhs */
    cs_real_t *gapinj = (st_prv_id > -1) ? c_st_prv : rhs;

    /* We increment rhs with -Gamma.var_prev. and fimp with Gamma */

    cs_volume_mass_injection_get_arrays(f_eps,
                                        &n_elts,
                                        &elt_ids,
                                        &mst_type,
                                        &mst_val,
                                        &mst_val_p);
    cs_mass_source_terms(1, /* iterns*/
                         1, /* dim */
                         n_elts,
                         elt_ids,
                         mst_type,
                         cell_f_vol,
                         cvara_ep,
                         mst_val,
                         mst_val_p,
                         rhs,
                         fimp,
                         gapinj);

  }

  /* Unsteady term
   * ------------- */

  if (eqp->istat == 1) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      fimp[c_id] += (crom[c_id] / dt[c_id]) * cell_f_vol[c_id];
    });
  }

  /* Production (rho * Ce1 * epsilon / k * P)
   *    Dissipation (rho*Ce2.epsilon/k*epsilon)
   * ------------------------------------------ */

  cs_real_t thetap = (st_prv_id > -1) ? thetv : 1.;

  cs_real_t *cprod;
  CS_MALLOC_HD(cprod, n_cells_ext, cs_real_t, cs_alloc_mode);

  /* Calculation the production trace, depending we are in standard
   * Rij or in SSG (use of produc or grdvit) */

  if (cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_LRR) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cprod[c_id] = 0.5*(produc[c_id][0] + produc[c_id][1] + produc[c_id][2]);
    });
  }
  else {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cprod[c_id] = -(  cvara_rij[c_id][0] * gradv[c_id][0][0]
                      + cvara_rij[c_id][3] * gradv[c_id][0][1]
                      + cvara_rij[c_id][5] * gradv[c_id][0][2]
                      + cvara_rij[c_id][3] * gradv[c_id][1][0]
                      + cvara_rij[c_id][1] * gradv[c_id][1][1]
                      + cvara_rij[c_id][4] * gradv[c_id][1][2]
                      + cvara_rij[c_id][5] * gradv[c_id][2][0]
                      + cvara_rij[c_id][4] * gradv[c_id][2][1]
                      + cvara_rij[c_id][2] * gradv[c_id][2][2]);
    });
  }

  /* EBRSM */
  if (cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM) {

    const cs_real_t xa1 = cs_turb_xa1;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      /* Half-traces */
      const cs_real_t trprod = cprod[c_id];
      const cs_real_t trrij = 0.5 * cs_math_6_trace(cvara_rij[c_id]);

      /* Calculation of the Durbin time scale */
      const cs_real_t xttke = trrij / cvara_ep[c_id];
      const cs_real_t xttkmg
        = xct*sqrt(viscl[c_id] / crom[c_id] / cvara_ep[c_id]);
      const cs_real_t xttdrb = cs_math_fmax(xttke, xttkmg);

      const cs_real_t prdeps = trprod / cvara_ep[c_id];
      const cs_real_t alpha3 = cs_math_pow3(cvar_al[c_id]);

      /* Production (explicit) */
      /* Compute of C_eps_1' */
      w1[c_id] =   cromo[c_id] * cell_f_vol[c_id] * ce1
                 * (1.+xa1*(1-alpha3)*prdeps) * trprod / xttdrb;

      /* Dissipation (implicit) */
      const cs_real_t crom_vol = crom[c_id] * cell_f_vol[c_id];
      rhs[c_id] -= crom_vol * ceps2 * cvara_ep[c_id] / xttdrb;
      fimp[c_id] += ceps2 * crom_vol * thetap / xttdrb;
    });

  }

  /* SSG and LRR */
  else {

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      /* Half-traces */
      const cs_real_t trprod = cprod[c_id];
      const cs_real_t trrij = 0.5 * cs_math_6_trace(cvara_rij[c_id]);
      const cs_real_t xttke = trrij / cvara_ep[c_id];

      /* Production (explicit, a part might be implicit) */
      const cs_real_t cromo_vol = cromo[c_id] * cell_f_vol[c_id];
      fimp[c_id]
        += cs_math_fmax(- cromo_vol * ce1 * trprod / trrij,
                        0.0);
      w1[c_id] = cromo_vol * ce1 / xttke * trprod;

      /* Dissipation (implicit) */
      const cs_real_t crom_vol = crom[c_id] * cell_f_vol[c_id];
      rhs[c_id]  -= crom_vol * ceps2 * cs_math_pow2(cvara_ep[c_id]) / trrij;
      fimp[c_id] += ceps2 * crom_vol / xttke * thetap;
    });

  }

  /* Extrapolation of source terms (2nd order in time) */
  if (st_prv_id > -1) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      c_st_prv[c_id] += w1[c_id];
    });
  }
  else {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      rhs[c_id] += w1[c_id];
    });
  }

  ctx.wait();

  /* Buoyancy term
   * ------------- */

  /* FIXME use beta ... WARNING */
  if (cs_glob_turb_rans_model->has_buoyant_term == 1) {

    /* Extrapolation of source terms (2nd order in time) */
    if (st_prv_id > -1)
      _gravity_st_epsilon(phase_id, up_rhop, cell_f_vol, grav, c_st_prv);

    else
      _gravity_st_epsilon(phase_id, up_rhop, cell_f_vol, grav, rhs);

  }

  /* Diffusion term (Daly Harlow: generalized gradient hypothesis method)
   * -------------------------------------------------------------------- */

  cs_real_t *weighb;
  cs_real_6_t *viscce;
  cs_real_2_t *weighf;

  CS_MALLOC_HD(weighb, n_b_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(weighf, n_i_faces, cs_real_2_t, cs_alloc_mode);
  CS_MALLOC_HD(viscce, n_cells_ext, cs_real_6_t, cs_alloc_mode);

  /* Symmetric tensor diffusivity (GGDH) */
  if (eqp->idften & CS_ANISOTROPIC_DIFFUSION) {

    const cs_field_t *f_a_t_visc
      = cs_field_by_name("anisotropic_turbulent_viscosity");
    const cs_real_6_t *visten = (const cs_real_6_t *)f_a_t_visc->val;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t i = 0; i < 3; i++)
        viscce[c_id][i] = visten[c_id][i] / sigmae + viscl[c_id];
      for (cs_lnum_t i = 3; i < 6; i++)
        viscce[c_id][i] = visten[c_id][i] / sigmae;
    });
    ctx.wait();

    cs_face_anisotropic_viscosity_scalar(m,
                                         fvq,
                                         viscce,
                                         eqp->verbosity,
                                         weighf,
                                         weighb,
                                         viscf,
                                         viscb);

  }

  /* Scalar diffusivity */
  else {

    if (eqp->idifft == 1) {
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        w1[c_id] = viscl[c_id] + visct[c_id]/sigmae;
      });
      ctx.wait();
    }
    else
      cs_array_copy<cs_real_t>(n_cells, viscl, w1);

    cs_face_viscosity(m,
                      fvq,
                      eqp->imvisf,
                      w1,
                      viscf,
                      viscb);
  }

  /* Solving
   * ------- */

  if (st_prv_id > -1) {
    const cs_real_t thetp1 = 1.+thets;
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      rhs[c_id] += thetp1*c_st_prv[c_id];
    });
    ctx.wait();
  }

  cs_solid_zone_set_zero_on_cells(1, rhs);

  /* Get boundary conditions coefficients */

  cs_field_bc_coeffs_t *bc_coeffs_eps = f_eps->bc_coeffs;

  cs_equation_param_t eqp_loc = *eqp;

  eqp_loc.istat  = -1;
  eqp_loc.icoupl = -1;
  eqp_loc.idifft = -1;
  eqp_loc.iwgrec = 0;     /* Warning, may be overwritten if a field */
  eqp_loc.theta = thetv;
  eqp_loc.blend_st = 0;   /* Warning, may be overwritten if a field */

  cs_real_t *dpvar;
  CS_MALLOC_HD(dpvar, n_cells_ext, cs_real_t, cs_alloc_mode);

  cs_equation_iterative_solve_scalar(cs_glob_time_step_options->idtvar,
                                     1,   /* init */
                                     f_eps->id,
                                     nullptr,
                                     0,   /* iescap */
                                     0,   /* imucpp */
                                     -1,  /* normp  */
                                     &eqp_loc,
                                     cvara_ep,
                                     cvara_ep,
                                     bc_coeffs_eps,
                                     imasfl,
                                     bmasfl,
                                     viscf,
                                     viscb,
                                     viscf,
                                     viscb,
                                     viscce,
                                     weighf,
                                     weighb,
                                     0,  /* boundary convective upwind flux */
                                     nullptr,
                                     fimp,
                                     rhs,
                                     cvar_ep,
                                     dpvar,
                                     nullptr,
                                     nullptr);

  /* Free memory */

  CS_FREE_HD(dpvar);
  CS_FREE_HD(w1);
  CS_FREE_HD(cprod);
  CS_FREE_HD(viscce);
  CS_FREE_HD(weighb);
  CS_FREE_HD(weighf);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \brief Solve the \f$ R_{ij} - \epsilon \f$ for incompressible flows or
 *         slightly compressible flows for one time step.
 *
 * Please refer to the
 * <a href="../../theory.pdf#rijeps"><b>\f$ R_{ij} - \epsilon \f$ model</b></a>
 * section of the theory guide for more informations, as well as the
 * <a href="../../theory.pdf#turrij"><b>turrij</b></a> section.
 *
 * \param[in]     phase_id     turbulent phase id (-1 for single phase flow)
 !*/
/*-----------------------------------------------------------------------------*/

void
cs_turbulence_rij(int phase_id)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_turb_model_t *turb_model = cs_glob_turb_model;
  const cs_turb_rans_model_t *turb_rans_model = cs_glob_turb_rans_model;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  cs_dispatch_context ctx;
  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

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

  const cs_real_t *cell_f_vol = fvq->cell_f_vol;

  cs_field_t *f_rij = CS_F_(rij);
  cs_field_t *f_eps = CS_F_(eps);
  cs_field_t *f_alpbl = CS_F_(alp_bl);
  cs_field_t *f_vel = CS_F_(vel);
  cs_field_t *f_rho = CS_F_(rho);

  if (phase_id >= 0) {
    f_rij = CS_FI_(rij, phase_id);
    f_eps = CS_FI_(eps, phase_id);
    f_alpbl = CS_FI_(alp_bl, phase_id);
    f_vel = CS_FI_(vel, phase_id);
    f_rho = CS_FI_(rho, phase_id);
  }

  cs_real_t *cvar_ep = f_eps->val;
  cs_real_6_t *cvar_rij = (cs_real_6_t *)f_rij->val;

  const cs_real_t *dt = CS_F_(dt)->val;
  const cs_real_t *crom = f_rho->val;
  const cs_real_6_t *cvara_rij = (const cs_real_6_t *)f_rij->val_pre;

  const auto t2v = _t2v;
  const auto iv2t = _iv2t;
  const auto jv2t = _jv2t;
  const cs_real_t tdeltij[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
  const cs_real_t vdeltij[6] = {1, 1, 1, 0, 0, 0};

  cs_real_6_t *rhs;
  cs_real_66_t *fimp;
  CS_MALLOC_HD(rhs, n_cells_ext, cs_real_6_t, cs_alloc_mode);
  CS_MALLOC_HD(fimp,  n_cells_ext, cs_real_66_t, cs_alloc_mode);

  /* Source terms for Rij
   * -------------------- */

  cs_arrays_set_value<cs_real_t, 1>(6*n_cells, 0., (cs_real_t*)rhs);
  cs_arrays_set_value<cs_real_t, 1>(36*n_cells, 0., (cs_real_t*)fimp);

  cs_user_source_terms(cs_glob_domain,
                       f_rij->id,
                       (cs_real_t*)rhs,
                       (cs_real_t*)fimp);

  /* Time extrapolation ? */
  cs_real_6_t *c_st_prv = nullptr;
  int kstprv = cs_field_key_id("source_term_prev_id");
  int st_prv_id = cs_field_get_key_int(f_rij, kstprv);
  if (st_prv_id > -1)
    c_st_prv = (cs_real_6_t *)cs_field_by_id(st_prv_id)->val;

  const cs_equation_param_t *eqp
    = cs_field_get_equation_param_const(f_rij);

  const cs_time_step_t *time_step = cs_glob_time_step;
  const cs_time_scheme_t *time_scheme = cs_glob_time_scheme;

  const cs_real_t thets = time_scheme->thetst;
  const cs_real_t thetv = eqp->theta;

  if (c_st_prv != nullptr) {

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t ij = 0; ij < 6; ij++) {
        const cs_real_t tuexpr = c_st_prv[c_id][ij];
        /* For continuation and the next time step */
        c_st_prv[c_id][ij] = rhs[c_id][ij];

        rhs[c_id][ij] = -thets*tuexpr;
        /* Right hand side of the previous time step
         * We suppose -fimp > 0: we implicit
         * the user source term (the rest) */
        for (cs_lnum_t kl = 0; kl < 6; kl++) {
          rhs[c_id][ij] += fimp[c_id][ij][kl] * cvara_rij[c_id][kl];
          /* Diagonal */
          fimp[c_id][ij][kl] *= -thetv;
        }
      }
    });

  }
  else {

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t ij = 0; ij < 6; ij++) {
        for (cs_lnum_t kl = 0; kl < 6; kl++) {
          rhs[c_id][ij] += fimp[c_id][ij][kl] * cvara_rij[c_id][kl];
          fimp[c_id][ij][kl] = cs_math_fmax(-fimp[c_id][ij][kl], 0.);
        }
      }
    });

  }

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  const int iflmas =  cs_field_get_key_int(f_vel, kimasf);
  const int iflmab =  cs_field_get_key_int(f_vel, kbmasf);

  const cs_real_t *imasfl = cs_field_by_id(iflmas)->val;
  const cs_real_t *bmasfl = cs_field_by_id(iflmab)->val;

  const cs_real_t ro0 = cs_glob_fluid_properties->ro0;
  const cs_real_t viscl0 = cs_glob_fluid_properties->viscl0;

  if (eqp->verbosity >= 1) {
    const char *f_label = cs_field_get_label(f_rij);
    switch(turb_model->iturb) {
      case CS_TURB_RIJ_EPSILON_LRR:
        bft_printf(" ** Solving Rij-EPSILON LRR %s\n"
                   "    -----------------------\n", f_label);
        break;
    case CS_TURB_RIJ_EPSILON_SSG:
      bft_printf(" ** Solving Rij-EPSILON SSG %s\n"
                 "    -----------------------\n", f_label);
        break;
    case CS_TURB_RIJ_EPSILON_EBRSM:
      bft_printf(" ** Solving Rij-EPSILON EBRSM %s\n"
                 "    -------------------------\n", f_label);
      break;
    default:
      assert(0);
    }
  }

  /* Allocate or map temporary arrays for the turbulence resolution,
     depending on user options */

  cs_real_3_t *up_rhop = nullptr;

  cs_real_6_t *cpro_press_correl = nullptr;
  {
    cs_field_t *f_psc = cs_field_by_name_try
                          ("algo:rij_pressure_strain_correlation");
    if (f_psc != nullptr)
      cpro_press_correl = (cs_real_6_t *)f_psc->val;
  }

  cs_real_6_t *produc = nullptr, *_produc = nullptr;
  {
    cs_field_t *f_rij_p = cs_field_by_name_try
                             ("algo:rij_production");
    if (f_rij_p != nullptr)
      produc = (cs_real_6_t *)f_rij_p->val;
    else {
      CS_MALLOC_HD(_produc, n_cells_ext, cs_real_6_t, cs_alloc_mode);
      produc = _produc;
    }
  }

  cs_real_t *viscf, *viscb;
  CS_MALLOC_HD(viscf, n_i_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(viscb, n_b_faces, cs_real_t, cs_alloc_mode);

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

  /* Advanced initialization for EBRSM
   * -------------------------------- */

  /* Automatic reinitialization at the end of the first iteration:
   * wall distance y^+ is computed with -C log(1-alpha), where C=CL*Ceta*L*kappa,
   * then y so we have an idea of the wall distance in complex geometries.
   * Then U is initialized with a Reichard layer,
   * Epsilon by 1/(kappa y), clipped next to the wall at its value for y^+=15
   * k is given by a blending between eps/(2 nu)*y^2 and utau/sqrt(Cmu).
   * The blending function is chosen so that the asymptotic behavior
   * and give the correct peak of k. */

  /* TODO FIXME Are the BCs uncompatible ? */
  if (   time_step->nt_cur == 1
      && turb_rans_model->reinit_turb == 1
      && turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM) {
    cs_real_t *cvar_al = f_alpbl->val;
    cs_real_3_t *vel = (cs_real_3_t *)f_vel->val;

    const cs_real_t uref = cs_glob_turb_ref_values->uref;
    const cs_real_t utaurf = 0.05 * uref;
    const cs_real_t nu0 = viscl0 / ro0;
    const cs_real_t xkappa = cs_turb_xkappa;

    cs_real_3_t *grad;
    CS_MALLOC_HD(grad, n_cells_ext, cs_real_3_t, cs_alloc_mode);

    /* Compute the gradient of Alpha */
    cs_field_gradient_scalar(f_alpbl,
                             false,  /* use_previous_t */
                             1,      /* inc */
                             grad);

    constexpr cs_real_t c_1ov3 = 1./3.;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

      /* Velocity magnitude */
      cs_real_t xunorm = cs_math_3_norm(vel[c_id]);

      /* y+ is bounded by 400, because in the Reichard profile,
       * it corresponds to saturation (u>uref) */
      cvar_al[c_id] = cs_math_fmax(cs_math_fmin(cvar_al[c_id], (1.-exp(-8.))),
                                   0.);
      /* Magnitude and unit vector of the Alpha gradient */
      cs_real_t xnal[3];
      cs_real_t xnoral = cs_math_3_norm(grad[c_id]);
      if (xnoral <= cs_math_epzero / pow(cell_f_vol[c_id], c_1ov3)) {
        for (int i = 0; i < 3; i++)
          xnal[i] = 1.0 / sqrt(3.0);
      }
      else {
        for (int i = 0; i < 3; i++)
          xnal[i] = grad[c_id][i]/xnoral;
      }

      const cs_real_t alpha3 = cs_math_pow3(cvar_al[c_id]);

      /* Compute YA, therefore alpha is given by 1-exp(-YA/(50 nu/utau))
       * NB: y^+ = 50 give the best compromise */
      const cs_real_t ya = -log(1.-cvar_al[c_id]) * 50 * nu0 / utaurf;
      const cs_real_t ypa = ya / (nu0 / utaurf);
      /* Velocity magnitude is imposed (limited only), the direction is
       * conserved */
      cs_real_t limiter = 1.0;
      if (xunorm > 1.e-12 * uref)
        limiter = cs_math_fmin(   utaurf / xunorm
                               * (  2.5 * log(1.+0.4*ypa)
                                  + 7.80*(  1. - exp(-ypa/11.0)
                                          - (ypa/11.)*exp(-0.33*ypa))),
                                  1.0);

      for (cs_lnum_t i = 0; i < 3; i++)
        vel[c_id][i] = limiter * vel[c_id][i];

      const cs_real_t ut2 = 0.050 * uref;
      cvar_ep[c_id] =   cs_math_pow3(utaurf)
                      * cs_math_fmin(1. / (xkappa * 15.0 * nu0 / utaurf),
                                     1. / (xkappa * ya));
      const cs_real_t tke =     cvar_ep[c_id] * 0.5 / nu0*cs_math_pow2(ya)
                              * cs_math_pow2(exp(-ypa/25.))
                            +   cs_math_pow2(ut2) / 0.3
                              * cs_math_pow2((1.-exp(-ypa/25.)));

      for (cs_lnum_t ij = 0; ij < 6; ij++) {
        cs_lnum_t i = iv2t[ij];
        cs_lnum_t j = jv2t[ij];
        cvar_rij[c_id][ij] =      alpha3  * 2./3 * tke * vdeltij[ij]
                            + (1.-alpha3)*(vdeltij[ij] - xnal[i]*xnal[j])*tke;
      }

    }); /* End of loop on cells */

    ctx.wait();

    CS_FREE_HD(grad);

    cs_field_current_to_previous(f_vel);
    cs_field_current_to_previous(f_rij);
  }

  /* Compute the velocity gradient */

  if (f_vel->grad == nullptr)
    cs_field_gradient_vector(f_vel, true, 1, gradv);

  /* Compute the production term for Rij
   * ----------------------------------- */

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {

    /* Pij = - (Rik dUj/dXk + dUi/dXk Rkj)
     * Pij is stored as (P11, P22, P33, P12, P23, P13) */
    for (cs_lnum_t ij = 0; ij < 6; ij++) {
      cs_lnum_t i = iv2t[ij];
      cs_lnum_t j = jv2t[ij];

      produc[c_id][ij] = - (  cvara_rij[c_id][t2v[i][0]] * gradv[c_id][j][0]
                            + cvara_rij[c_id][t2v[i][1]] * gradv[c_id][j][1]
                            + cvara_rij[c_id][t2v[i][2]] * gradv[c_id][j][2]
                            + gradv[c_id][i][0] * cvara_rij[c_id][t2v[0][j]]
                            + gradv[c_id][i][1] * cvara_rij[c_id][t2v[1][j]]
                            + gradv[c_id][i][2] * cvara_rij[c_id][t2v[2][j]]);
    }

  });
  ctx.wait();

  /* Compute the pressure correlation  term for Rij
   * ----------------------------------------------
   * Phi,ij = Phi1,ij+Phi2,ij
   * Phi,ij = -C1 k/eps (Rij-2/3k dij) - C2 (Pij-2/3P dij)
   * Phi,ij is stored as (Phi11, Phi22, Phi33, Phi12, Phi23, Phi13)

   * TODO : coherency with the model */

  if (cpro_press_correl != nullptr)  {
    const cs_real_t d2s3 = 2./3;

    const cs_real_t crij1 = cs_turb_crij1;
    const cs_real_t crij2 = cs_turb_crij2;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      const cs_real_t k = 0.5 * cs_math_6_trace(cvara_rij[c_id]);
      const cs_real_t p = 0.5 * cs_math_6_trace(produc[c_id]);
      for (cs_lnum_t ij = 0; ij < 6; ij++)
        cpro_press_correl[c_id][ij]
          = - crij1 * cvar_ep[c_id] / k
                    * (cvara_rij[c_id][ij] - d2s3*k* vdeltij[ij])
            - crij2 * (produc[c_id][ij] - d2s3 * p * vdeltij[ij]);
    });
  }

  /* Compute the density gradient for buoyant terms
   * ---------------------------------------------- */

  if (turb_rans_model->has_buoyant_term == 1) {
    CS_MALLOC_HD(up_rhop, n_cells_ext, cs_real_3_t, cs_alloc_mode);
    _compute_up_rhop(phase_id, up_rhop);
  }

  /* Prepare to solve Rij, in a manner similar
     to that of cs_solve_equation_scalar.
   * ========================================= */

  if (cs_glob_porous_model == 3)
    cs_immersed_boundary_wall_functions(f_rij->id,
                                        (cs_real_t*)rhs,
                                        (cs_real_t*)fimp);

  /* Lagrangian source terms
   * ----------------------- */

  /* Second order is not taken into account  */

  if (   cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING
      && cs_glob_lagr_source_terms->ltsdyn == 1) {
    const cs_real_6_t *lagr_st_rij
      = (const cs_real_6_t *)cs_field_by_name_try("lagr_st_rij")->val;

    cs_real_t *lag_st_i = cs_field_by_name("lagr_st_imp_velocity")->val;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cs_real_t st_i = cell_f_vol[c_id] * lag_st_i[c_id];
      for (cs_lnum_t ij = 0; ij < 6; ij++) {
        rhs[c_id][ij] += cell_f_vol[c_id] * lagr_st_rij[c_id][ij];
        fimp[c_id][ij][ij] += cs_math_fmax(-st_i, 0.);
      }
    });
    ctx.wait();
  }

  /* Mass source terms
   *------------------ */

  if (eqp->n_volume_mass_injections > 0) {
    int *mst_type = nullptr;
    cs_lnum_t n_elts = 0;
    const cs_lnum_t *elt_ids = nullptr;
    cs_real_t *mst_val = nullptr, *mst_val_p = nullptr;

    cs_volume_mass_injection_get_arrays(f_rij,
                                        &n_elts,
                                        &elt_ids,
                                        &mst_type,
                                        &mst_val,
                                        &mst_val_p);

    /* If we extrapolate the source terms, we put Gamma Pinj in c_st_prv;
        Otherwise we put it directly in rhs */
    cs_real_6_t *gatinj = (st_prv_id > -1) ? c_st_prv : rhs;

    /* We increment rhs with -Gamma.var_prev. and rovsdr with Gamma */

    cs_mass_source_terms(1, /* iterns*/
                         6, /* dim */
                         n_elts,
                         elt_ids,
                         mst_type,
                         cell_f_vol,
                         (const cs_real_t*)cvara_rij,
                         mst_val,
                         mst_val_p,
                         (cs_real_t*)rhs,
                         (cs_real_t*)fimp,
                         (cs_real_t*)gatinj);
  }

  /* Unsteady term
   * ------------- */
  if (eqp->istat == 1) {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (int ij = 0; ij < 6; ij++)
        fimp[c_id][ij][ij] += (crom[c_id]/dt[c_id])*cell_f_vol[c_id];
    });
    ctx.wait();
  }

  /* Terms specific to Rij-epsilon model
   * ----------------------------------- */

  cs_real_t *weighb;
  cs_real_6_t *viscce;
  cs_real_2_t *weighf;

  CS_MALLOC_HD(weighb, n_b_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(weighf, n_i_faces, cs_real_2_t, cs_alloc_mode);
  CS_MALLOC_HD(viscce, n_cells_ext, cs_real_6_t, cs_alloc_mode);

  if (turb_model->iturb == CS_TURB_RIJ_EPSILON_LRR) {
    if (turb_rans_model->irijco == 1)
      _pre_solve_lrr(f_rij, phase_id, gradv,
                     produc, up_rhop, grav,
                     viscf, viscb, viscce,
                     rhs, fimp,
                     weighf, weighb);
    else
      _pre_solve_lrr_sg(f_rij, phase_id,
                        produc, up_rhop, grav,
                        viscf, viscb, viscce,
                        rhs, fimp,
                        weighf, weighb);

  }
  else { /* if (   turb_model->iturb == CS_TURB_RIJ_EPSILON_SSG
                || turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM) */
    _pre_solve_ssg(f_rij, phase_id, gradv,
                   produc, up_rhop, grav,
                   viscf, viscb, viscce,
                   rhs, fimp,
                   weighf, weighb);
  }

  cs_real_6_t  *coefap = (cs_real_6_t *)f_rij->bc_coeffs->a;
  cs_real_6_t  *cofafp = (cs_real_6_t *)f_rij->bc_coeffs->af;
  cs_real_66_t *cofbfp = (cs_real_66_t *)f_rij->bc_coeffs->bf;

  /* Add Rusanov fluxes */
  if (cs_glob_turb_rans_model->irijnu == 2) {
    cs_real_t *ipro_rusanov = cs_field_by_name("i_rusanov_diff")->val;
    ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      viscf[face_id] = cs_math_fmax(0.5 * ipro_rusanov[face_id],
                                    viscf[face_id]);
    });

    const cs_real_3_t *restrict b_face_u_normal
      = (const cs_real_3_t *)fvq->b_face_u_normal;
    cs_real_t *b_lam = cs_field_by_name("b_rusanov_diff")->val;

    ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
      const cs_real_t *n = b_face_u_normal[face_id];
      cs_real_t bf[6][6];

      for (cs_lnum_t ij = 0; ij < 6; ij++) {
        cs_lnum_t i = iv2t[ij];
        cs_lnum_t j = jv2t[ij];
        for (cs_lnum_t kl = 0; kl < 6; kl++) {
          cs_lnum_t k = iv2t[kl];
          cs_lnum_t l = jv2t[kl];
          bf[ij][kl] = b_lam[face_id] * n[l]
                       * (  n[i] * (tdeltij[j][k] - n[j]*n[k])
                          + n[j] * (tdeltij[i][k] - n[i]*n[k]));
        }
      }

      for (cs_lnum_t i = 0; i < 6; i++) {
        for (cs_lnum_t j = 0; j < 6; j++) {
          cofbfp[face_id][i][j] += bf[i][j];
          cofafp[face_id][i] -= bf[i][j] * coefap[face_id][j];
        }
      }

    });

  }

  /* Solve Rij
   * --------- */

  if (c_st_prv != nullptr) {
    const cs_real_t thetp1 = 1. + thets;
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t ij = 0; ij < 6; ij++)
        rhs[c_id][ij] += thetp1 * c_st_prv[c_id][ij];
    });
  }

  ctx.wait();

  cs_solid_zone_set_zero_on_cells(6, (cs_real_t *)rhs);

  /* All boundary convective flux with upwind */
  int icvflb = 0;

  cs_equation_param_t eqp_loc = *eqp;
  eqp_loc.istat  = -1;
  eqp_loc.iwgrec = 0;     /* Warning, may be overwritten if a field */
  eqp_loc.theta = thetv;
  eqp_loc.blend_st = 0;   /* Warning, may be overwritten if a field */

  cs_equation_iterative_solve_tensor(cs_glob_time_step_options->idtvar,
                                     f_rij->id,
                                     nullptr,
                                     &eqp_loc,
                                     cvara_rij,
                                     cvara_rij,
                                     f_rij->bc_coeffs,
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
                                     cvar_rij);

  CS_FREE_HD(viscce);
  CS_FREE_HD(fimp);
  CS_FREE_HD(rhs);

  /* Solve Epsilon
   * ------------- */

  {
    cs_real_t *_rhs, *_fimp;
    CS_MALLOC_HD(_rhs, n_cells_ext, cs_real_t, cs_alloc_mode);
    CS_MALLOC_HD(_fimp, n_cells_ext, cs_real_t, cs_alloc_mode);

    _solve_epsilon(phase_id,
                   gradv,
                   produc,
                   up_rhop,
                   grav,
                   viscf,
                   viscb,
                   _rhs,
                   _fimp);

    CS_FREE_HD(_fimp);
    CS_FREE_HD(_rhs);
  }

  /* Clipping
   * -------- */

  cs_turbulence_rij_clip(phase_id, n_cells);

  /* Free memory */

  CS_FREE_HD(_gradv);
  produc = nullptr;
  CS_FREE_HD(_produc);

  CS_FREE_HD(up_rhop);

  CS_FREE_HD(viscf);
  CS_FREE_HD(viscb);
  CS_FREE_HD(weighb);
  CS_FREE_HD(weighf);
#if defined(HAVE_ACCEL)
  CS_FREE_HD(_grav);
#endif

  if (cs_glob_timer_kernels_flag > 0) {
    ctx.wait();
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();

    std::chrono::microseconds elapsed;
    printf("%d: %s", cs_glob_rank_id, __func__);

    elapsed = std::chrono::duration_cast
      <std::chrono::microseconds>(t_stop - t_start);
    printf(", total = %ld\n", elapsed.count());
  }
}

/*----------------------------------------------------------------------------*/
/*! \brief Solve the equation on alpha in the framework of the Rij-EBRSM model.
 *
 * Also called for alpha of scalars for EB-DFM.
 *
 * \param[in]  f_id          field id of alpha variable
 * \param[in]  phase_id      turbulent phase id (-1 for single phase flow)
 * \param[in]  c_durbin_l    constant for the Durbin length
 !*/
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_solve_alpha(int        f_id,
                              int        phase_id,
                              cs_real_t  c_durbin_l)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  cs_dispatch_context ctx;
  cs_dispatch_sum_type_t i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);
  cs_dispatch_sum_type_t b_sum_type = ctx.get_parallel_for_b_faces_sum_type(m);

  const cs_real_t *distb = fvq->b_dist;
  const cs_real_t *cell_f_vol = fvq->cell_f_vol;
  const cs_lnum_t *b_face_cells = m->b_face_cells;
  const cs_lnum_2_t *i_face_cells = m->i_face_cells;

  cs_real_t *cvar_al = cs_field_by_id(f_id)->val;

  cs_field_t *f_rho = CS_F_(rho);
  cs_field_t *f_mu = CS_F_(mu);
  cs_field_t *f_rij = CS_F_(rij);
  cs_field_t *f_eps = CS_F_(eps);
  cs_field_t *f_vel = CS_F_(vel);

  if (phase_id >= 0) {
    f_rho = CS_FI_(rho, phase_id);
    f_mu = CS_FI_(mu, phase_id);
    f_rij = CS_FI_(rij, phase_id);
    f_eps = CS_FI_(eps, phase_id);
    f_vel = CS_FI_(vel, phase_id);
  }

  const cs_real_t *crom = f_rho->val;
  const cs_real_t *viscl = f_mu->val;
  const cs_real_t *cvara_ep = f_eps->val_pre;
  const cs_real_t *cvara_al = cs_field_by_id(f_id)->val_pre;
  const cs_real_6_t *cvara_rij = (const cs_real_6_t *)f_rij->val_pre;

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  int iflmas =  cs_field_get_key_int(f_vel, kimasf);
  int iflmab =  cs_field_get_key_int(f_vel, kbmasf);

  const cs_real_t *imasfl = cs_field_by_id(iflmas)->val;
  const cs_real_t *bmasfl = cs_field_by_id(iflmab)->val;

  const cs_real_t d1s2 = 0.50;
  const cs_real_t d1s4 = 0.25;
  const cs_real_t d3s2 = 1.50;
  const cs_real_t uref = cs_glob_turb_ref_values->uref;

  /* Resolving the equation of alpha
     =============================== */

  const cs_equation_param_t *eqp
    = cs_field_get_equation_param(cs_field_by_id(f_id));

  if (eqp->verbosity == 1) {
    cs_log_printf(CS_LOG_DEFAULT,
                  _(" Solving the variable %s\n"),
                  cs_field_get_label(cs_field_by_id(f_id)));
  }

  /* Allocate temporary arrays */
  cs_real_t *fimp, *rhs;
  CS_MALLOC_HD(rhs, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(fimp, n_cells_ext, cs_real_t, cs_alloc_mode);

  /* Source term of alpha
   *  \f$ rhs = \dfrac{1}{L^2 (\alpha)} - \dfrac{1}{L^2}\f$
   * In fact there is a mark "-" because the solved equation is
   *   \f$-\div{\grad {alpha}} = rhs \f$
   *================================================================*/

  /* Matrix */

  const cs_real_t thetv = eqp->theta;
  cs_real_t thetap = (cs_glob_time_scheme->isto2t > 0) ? thetv : 1.0;

  // FIXME the source term extrapolation is not well done!!!!

  /* For automatic initialization, the length scale is fixed at L^+ =50 */

  if (   cs_glob_time_step->nt_cur == 1
      && cs_glob_turb_rans_model->reinit_turb == 1) {

    const cs_real_t ro0 = cs_glob_fluid_properties->ro0;
    const cs_real_t viscl0 = cs_glob_fluid_properties->viscl0;
    const cs_real_t xlldrb = 50.0 * viscl0 / ro0 / (0.050*uref);
    const cs_real_t l2 = cs_math_pow2(xlldrb);

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      /* Explicit term */
      rhs[c_id] = cell_f_vol[c_id]*(1.0-cvara_al[c_id]) / l2;

      /* Implicit term */
      fimp[c_id] = cell_f_vol[c_id]*thetap / l2;
    });

  }
  else {

    const cs_real_t xceta = cs_turb_xceta;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      const cs_real_t xk
        = d1s2 * (cvara_rij[c_id][0] + cvara_rij[c_id][1] + cvara_rij[c_id][2]);
      const cs_real_t xnu = viscl[c_id] / crom[c_id];

      /* Integral length scale */
      const cs_real_t xllke = pow(xk, d3s2) / cvara_ep[c_id];

      /* Kolmogorov length scale */
      const cs_real_t xllkmg =   xceta
                               * pow(cs_math_pow3(xnu)/cvara_ep[c_id], d1s4);

      /* Durbin length scale */
      const cs_real_t xlldrb = c_durbin_l*cs_math_fmax(xllke, xllkmg);

      const cs_real_t l2 = cs_math_pow2(xlldrb);

      /* Explicit term */
      rhs[c_id] = cell_f_vol[c_id]*(1.0-cvara_al[c_id]) / l2;

      /* Implicit term */
      fimp[c_id] = cell_f_vol[c_id]*thetap / l2;
    });

  }

  ctx.wait();

  /* Calculation of viscf and viscb for cs_equation_iterative_solve_scalar. */
  cs_real_t *w1, *viscf, *viscb;
  CS_MALLOC_HD(w1, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(viscf, n_i_faces, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(viscb, n_b_faces, cs_real_t, cs_alloc_mode);

  cs_arrays_set_value<cs_real_t, 1>(n_cells, 1., w1);

  cs_face_viscosity(m,
                    fvq,
                    eqp->imvisf,
                    w1,
                    viscf,
                    viscb);

  CS_FREE_HD(w1);

  cs_solid_zone_set_zero_on_cells(1, rhs);

  /* Effective resolution of the equation of alpha
     ============================================= */

  cs_field_bc_coeffs_t *bc_coeffs = cs_field_by_id(f_id)->bc_coeffs;

  cs_equation_param_t eqp_loc = *eqp;

  eqp_loc.istat  = -1;
  eqp_loc.icoupl = -1;
  eqp_loc.idifft = -1;
  eqp_loc.iwgrec = 0;     /* Warning, may be overwritten if a field */
  eqp_loc.theta = thetv;
  eqp_loc.blend_st = 0;   /* Warning, may be overwritten if a field */

  cs_real_t *dpvar;
  CS_MALLOC_HD(dpvar, n_cells_ext, cs_real_t, cs_alloc_mode);

  cs_equation_iterative_solve_scalar(cs_glob_time_step_options->idtvar,
                                     1, /* init */
                                     f_id,
                                     nullptr,
                                     0, /* iescap */
                                     0, /* imucpp */
                                     -1, /* normp */
                                     &eqp_loc,
                                     cvara_al,
                                     cvara_al,
                                     bc_coeffs,
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
                                     fimp,
                                     rhs,
                                     cvar_al,
                                     dpvar,
                                     nullptr,
                                     nullptr);

  CS_FREE_HD(dpvar);
  CS_FREE_HD(rhs);

  cs_solid_zone_set_scalar_on_cells(1., cvar_al);

  /* Clipping
     ======== */

  cs_real_t *alpha_min;
  CS_MALLOC_HD(alpha_min, n_cells_ext, cs_real_t, cs_alloc_mode);

  /* Compute a first estimator of the minimal value of alpha per cell.
   * This is deduced from "alpha/L^2 - div(grad alpha) = 1/L^2" and assuming that
   * boundary cell values are 0. This value is thefore non zero but
   * much smaller than the wanted value. */

  cs_array_copy<cs_real_t>(n_cells, (const cs_real_t *)fimp, alpha_min);

  ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
    const cs_lnum_t c_id0 = i_face_cells[face_id][0];
    const cs_lnum_t c_id1 = i_face_cells[face_id][1];
    if (c_id0 < n_cells)
      cs_dispatch_sum(&alpha_min[c_id0], viscf[face_id], i_sum_type);
    if (c_id1 < n_cells)
      cs_dispatch_sum(&alpha_min[c_id1], viscf[face_id], i_sum_type);
  });

  ctx.parallel_for_b_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t  face_id) {
    const cs_lnum_t c_id = b_face_cells[face_id];
    cs_dispatch_sum(&alpha_min[c_id], viscb[face_id]/distb[face_id], b_sum_type);
  });

  ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
    alpha_min[c_id] = fimp[c_id]/alpha_min[c_id];
  });

  ctx.wait();

  _clip_alpha(f_id, n_cells, alpha_min);

  CS_FREE_HD(alpha_min);
  CS_FREE_HD(fimp);
  CS_FREE_HD(viscf);
  CS_FREE_HD(viscb);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize Rij-epsilon variables based on reference quantities.
 *
 * If uref is not provided (0 or negative), values are set at a large
 * negative value (-cs_math_big_r) to allow for later checks.
 *
 * \param[in]  uref    characteristic flow velocity
 * \param[in]  almax   characteristic macroscopic length of the domain
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_init_by_ref_quantities(cs_real_t  uref,
                                         cs_real_t  almax)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  cs_dispatch_context ctx;

  cs_real_t *cvar_ep = CS_F_(eps)->val;
  cs_real_6_t *cvar_rij = (cs_real_6_t *)CS_F_(rij)->val;

  /* With reference velocity */

  if (uref > 0) {
    const cs_real_t tr_ii = cs_math_pow2(0.02 * uref);
    const cs_real_t k = 0.5 * (3. * tr_ii);  /* trace of tensor with
                                                tr_ii diagonal) */
    const cs_real_t ep = pow(k, 1.5) * cs_turb_cmu / almax;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cvar_rij[c_id][0] = tr_ii;
      cvar_rij[c_id][1] = tr_ii;
      cvar_rij[c_id][2] = tr_ii;
      cvar_rij[c_id][3] = 0;
      cvar_rij[c_id][4] = 0;
      cvar_rij[c_id][5] = 0;
      cvar_ep[c_id] = ep;
    });
    ctx.wait();

    cs_turbulence_rij_clip(-1, n_cells);
  }

  /* Without reference velocity */
  else {
    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      for (cs_lnum_t i = 0; i < 6; i++)
        cvar_rij[c_id][i] = -cs_math_big_r;
      cvar_ep[c_id] = -cs_math_big_r;
    });
    ctx.wait();
  }

  cs_solid_zone_set_zero_on_cells(6, (cs_real_t *)cvar_rij);
  cs_solid_zone_set_scalar_on_cells(1e-12, cvar_ep);
  //cs_solid_zone_set_zero_on_cells(1, cvar_ep);

  /* For EBRSM, initialize alpha */

  if (cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM) {
    cs_real_t *cvar_al = CS_F_(alp_bl)->val;
    cs_arrays_set_value<cs_real_t, 1>(n_cells, 1., cvar_al);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Clip the turbulent Reynods stress tensor and the turbulent
 *        dissipation (coupled components version).
 *
 * \param[in]  phase_id   turbulent phase id (-1 for single phase flow)
 * \param[in]  n_cells    number of cells
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_clip(int        phase_id,
                       cs_lnum_t  n_cells)
{
  cs_field_t *f_rij = CS_F_(rij);
  cs_field_t *f_eps = CS_F_(eps);

  if (phase_id >= 0) {
    f_rij = CS_FI_(rij, phase_id);
    f_eps = CS_FI_(eps, phase_id);
  }

  cs_real_t *cvar_ep = (cs_real_t *)f_eps->val;
  cs_real_6_t *cvar_rij = (cs_real_6_t *)f_rij->val;
  const cs_real_t *cvara_ep = (const cs_real_t *)f_eps->val_pre;

  int kclipp = cs_field_key_id("clipping_id");

  /* Post-process clippings ? */

  cs_real_t *cpro_eps_clipped = nullptr;
  cs_real_6_t *cpro_rij_clipped = nullptr;
  int clip_e_id = cs_field_get_key_int(f_eps, kclipp);
  if (clip_e_id > -1) {
    cpro_eps_clipped = cs_field_by_id(clip_e_id)->val;
    cs_array_real_fill_zero(n_cells, cpro_eps_clipped);
  }
  int clip_r_id = cs_field_get_key_int(f_rij, kclipp);
  if (clip_r_id > -1) {
    cpro_rij_clipped = (cs_real_6_t *)cs_field_by_id(clip_r_id)->val;
    cs_array_real_fill_zero(n_cells*f_rij->dim,
                            (cs_real_t *)cpro_rij_clipped);
  }

  /* Compute and store Min Max values for the log. */

  cs_real_t vmin[7], vmax[7];
  _rij_min_max(n_cells, cvar_rij, cvar_ep, vmin, vmax);

  /* Clipping (modified to avoid exactly zero values) */

  const cs_real_t varrel = 1.1;
  const cs_real_t eigen_tol = 1.e-4;
  const cs_real_t epz2 = cs_math_pow2(cs_math_epzero);

  cs_lnum_t icltot = 0;
  cs_lnum_t iclrij[6] = {0, 0, 0, 0, 0, 0};
  cs_lnum_t iclep[1] = {0};

  /* Compute the maximal value of each of the diagonal components over the
   * entire domain. A reference value "rijref", used to determine if a value is
   * small is then calculated as: rijref = (r11max + r22max + r33max)/3.
   * New test is rii < epzero*rijref  */

  cs_real_t rijmax[3] = {0, 0, 0};
  for (cs_lnum_t i = 0; i < 3; i++) {
    if (vmax[i] > rijmax[i])
      rijmax[i] = vmax[i];
  }

  cs_parall_max(3, CS_REAL_TYPE, rijmax);

  const cs_real_t trref = rijmax[0] + rijmax[1] + rijmax[2];
  const cs_real_t rijref = cs_math_fmax(trref/3., cs_math_epzero);

  cs_lnum_t solid_stride = 1;
  int *c_is_solid_zone_flag = cs_solid_zone_flag(cs_glob_mesh);
  const int c_is_solid_ref[1] = {0};
  const int *c_is_solid = c_is_solid_zone_flag;
  if (c_is_solid == nullptr) {
    c_is_solid = c_is_solid_ref;
    solid_stride = 0;
  }

# pragma omp parallel if(n_cells > CS_THR_MIN)
  {
    cs_lnum_t t_icltot = 0;
    cs_lnum_t t_iclrij[6] = {0, 0, 0, 0, 0, 0};
    cs_lnum_t t_iclep[1] = {0};

    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n_cells, sizeof(cs_real_t), &s_id, &e_id);

    for (cs_lnum_t c_id = s_id; c_id < e_id; c_id++) {

      int is_clipped = 0;

      /* Special case for solid cells (which are set to 0 but should
         not count as clippings) */

      if (c_is_solid[solid_stride*c_id]) {
        for (cs_lnum_t ij = 0; ij < 6; ij++)
          cvar_rij[c_id][ij] = 0;
        cvar_ep[c_id] = 1e-12;
        continue;
      }

      /* Check if R is positive and ill-conditioned (since the former
       * will induce the latter after clipping...) */

      const cs_real_t trrij = cs_math_6_trace(cvar_rij[c_id]);

      if (trrij <= cs_math_epzero*trref) {
        for (cs_lnum_t i = 0; i < 3; i++) {
          if (cpro_rij_clipped != nullptr) {
            cpro_rij_clipped[c_id][i]
              = cvar_rij[c_id][i] - cs_math_epzero*rijref;
            cpro_rij_clipped[c_id][i+3] = cvar_rij[c_id][i+3];
          }

          cvar_rij[c_id][i] = cs_math_epzero*rijref;
          cvar_rij[c_id][i+3] = 0.0;

          t_iclrij[i]++;
          t_iclrij[i+3]++;
        }

        is_clipped = 1;
      }
      else {
        cs_real_t tensor[6];
        for (cs_lnum_t ij = 0; ij < 6; ij++)
          tensor[ij] = cvar_rij[c_id][ij]/trrij;

        cs_real_t eigen_vals[3];
        _sym_33_eigen(tensor, eigen_vals);

        cs_real_t eigen_min = eigen_vals[0];
        cs_real_t eigen_max = eigen_vals[0];
        for (cs_lnum_t i = 1; i < 3; i++) {
          eigen_min = cs_math_fmin(eigen_min, eigen_vals[i]);
          eigen_max = cs_math_fmax(eigen_max, eigen_vals[i]);
        }

        /* If negative eigenvalue, return to isotropy */

        if (   (eigen_min <= eigen_tol*eigen_max)
            || (eigen_min < cs_math_epzero)) {

          is_clipped = 1;

          eigen_min = cs_math_fmin(eigen_min, -eigen_tol);
          cs_real_t eigen_offset
            = cs_math_fmin(-eigen_min/(1.0/3.0-eigen_min)+0.1, 1.0);

          for (cs_lnum_t ij = 0; ij < 6; ij++) {
            cvar_rij[c_id][ij] = (1.0-eigen_offset)*cvar_rij[c_id][ij];

            if (ij < 3)
              cvar_rij[c_id][ij] += trrij*(eigen_offset+eigen_tol)/3.;

            if (cpro_rij_clipped != nullptr)
              cpro_rij_clipped[c_id][ij] = eigen_offset*cvar_rij[c_id][ij];

            t_iclrij[ij]++;
          }
        }
      }

      /* Epsilon */

      if (cs_math_fabs(cvar_ep[c_id]) < epz2) {
        t_iclep[0]++;
        if (cpro_eps_clipped != nullptr)
          cpro_eps_clipped[c_id] = cs_math_fabs(cvar_ep[c_id]-epz2);
        cvar_ep[c_id] = cs_math_fmax(cvar_ep[c_id],epz2);
      }
      else if (cvar_ep[c_id] <= 0) {
        t_iclep[0]++;
        if (cpro_eps_clipped != nullptr)
          cpro_eps_clipped[c_id] = 2*cs_math_fabs(cvar_ep[c_id]);
        cvar_ep[c_id] = cs_math_fmin(cs_math_fabs(cvar_ep[c_id]),
                                     varrel*cs_math_fabs(cvara_ep[c_id]));
      }

      /* Enforce Cauchy Schwartz inequality (only for x, y, z directions) */

      cs_real_t cvar_var1, cvar_var2;
      for (cs_lnum_t ij = 3; ij < 6; ij++) {
        if (ij == 3) {
          cvar_var1 = cvar_rij[c_id][0];
          cvar_var2 = cvar_rij[c_id][1];
        }
        else if (ij == 4) {
          cvar_var1 = cvar_rij[c_id][1];
          cvar_var2 = cvar_rij[c_id][2];
        }
        else if (ij == 5) {
          cvar_var1 = cvar_rij[c_id][0];
          cvar_var2 = cvar_rij[c_id][2];
        }

        const cs_real_t rijmin = sqrt(cvar_var1*cvar_var2);
        if (rijmin < cs_math_fabs(cvar_rij[c_id][ij])) {
          is_clipped = 1;
          if (cpro_rij_clipped != nullptr)
            cpro_rij_clipped[c_id][ij] = cvar_rij[c_id][ij];
          cvar_rij[c_id][ij] =   _sign(1., cvar_rij[c_id][ij])
                               * rijmin/(1.+cs_math_epzero);
          t_iclrij[ij]++;
        }
      }
      t_icltot += is_clipped;

    }  /* End of loop on cells */

    /* Sum over threads */

    for (cs_lnum_t ij = 0; ij < 6; ij++) {
      #pragma omp atomic
      iclrij[ij] += t_iclrij[ij];
    }

    #pragma omp atomic
    iclep[0] += t_iclep[0];

    #pragma omp atomic
    icltot += t_icltot;
  }

  BFT_FREE(c_is_solid_zone_flag);

  /* Store number of clippings for logging */

  cs_lnum_t iclrij_max[6] = {0, 0, 0, 0, 0, 0}, iclep_max[1] = {0};

  cs_log_iteration_clipping_field(f_rij->id, icltot, 0,
                                  vmin, vmax, iclrij, iclrij_max);

  cs_log_iteration_clipping_field(f_eps->id, iclep[0], 0,
                                  vmin+6, vmax+6, iclep, iclep_max);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the turbulent viscosity for the Reynolds Stress model.
 *
 * \param[in]  phase_id  turbulent phase id (-1 for single phase flow)
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_mu_t(int  phase_id)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  cs_dispatch_context ctx;

  /* Initialization
   * ============== */

  /* Map field arrays */

  const cs_field_t *f_rij = CS_F_(rij);
  const cs_field_t *f_eps = CS_F_(eps);
  const cs_field_t *f_alpbl = CS_F_(alp_bl);
  const cs_field_t *f_rho = CS_F_(rho);
  cs_field_t *f_mut = CS_F_(mu_t);

  if (phase_id >= 0) {
    f_rij = CS_FI_(rij, phase_id);
    f_eps = CS_FI_(eps, phase_id);
    f_alpbl = CS_FI_(alp_bl, phase_id);
    f_rho = CS_FI_(rho, phase_id);
    f_mut = CS_FI_(mu_t, phase_id);
  }

  const cs_real_6_t *cvar_rij = (const cs_real_6_t *)f_rij->val;
  const cs_real_t *cvar_ep = f_eps->val;
  const cs_real_t *crom = f_rho->val;
  cs_real_t *visct = f_mut->val;

  const cs_real_t cmu = cs_turb_cmu;

  /* EBRSM case */

  if (cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM) {

    cs_real_3_t *grad_al = nullptr;
    CS_MALLOC_HD(grad_al, n_cells_ext, cs_real_3_t, cs_alloc_mode);
    cs_field_gradient_scalar(f_alpbl, true, 1, grad_al);

    const cs_real_t *cvar_al = f_alpbl->val;

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cs_real_t xrij[3][3];
      xrij[0][0] = cvar_rij[c_id][0];
      xrij[1][1] = cvar_rij[c_id][1];
      xrij[2][2] = cvar_rij[c_id][2];
      xrij[1][0] = cvar_rij[c_id][3];
      xrij[2][1] = cvar_rij[c_id][4];
      xrij[2][0] = cvar_rij[c_id][5];
      xrij[0][1] = xrij[1][0];
      xrij[0][2] = xrij[2][0];
      xrij[1][2] = xrij[2][1];

      /* Compute the magnitude of the Alpha gradient */
      cs_real_t xnal[3];
      cs_math_3_normalize(grad_al[c_id], xnal);

      cs_real_t alpha3 = cs_math_pow3(cvar_al[c_id]);

      cs_real_t xk = 0.5 * cs_math_33_trace(xrij);
      cs_real_t xe = cvar_ep[c_id];

      /* We compute the normal Reynolds Stresses */

      cs_real_t xrnn = 0;
      for (cs_lnum_t i = 0; i < 3; i++) {
        for (cs_lnum_t j = 0; j < 3; j++)
          xrnn += xrij[i][j]*xnal[j]*xnal[i];
      }
      xrnn = (1.-alpha3)*xrnn + alpha3*xk;
      xrnn = cs_math_fmax(xrnn, 1.e-12);

      visct[c_id] = crom[c_id] * cmu * xrnn * xk / xe;
    });

    ctx.wait();
    CS_FREE_HD(grad_al);

  }

  /* SSG and LRR */

  else {

    ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
      cs_real_t xk = 0.5 * cs_math_6_trace(cvar_rij[c_id]);
      cs_real_t xe = cvar_ep[c_id];

      visct[c_id] = crom[c_id] * cmu * cs_math_pow2(xk) / xe;
    });

    ctx.wait();
  }

  /* Zero turbulent viscosity for solid cells */

  cs_solid_zone_set_zero_on_cells(1, visct);
}

/*----------------------------------------------------------------------------*/
/*! \brief Compute Rusanov equivalent diffusivity of the model.
 !*/
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_compute_rusanov(void)
{
  if (cs_glob_turb_rans_model->irijnu != 2)
    return;

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *)mq->i_face_normal;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *)mq->b_face_normal;
  const cs_lnum_2_t *i_face_cells = m->i_face_cells;
  const cs_lnum_t *b_face_cells = m->b_face_cells;

  const int *bc_type = cs_glob_bc_type;

  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_i_faces = m->n_i_faces;

  cs_real_t *ipro_rusanov = cs_field_by_name("i_rusanov_diff")->val;
  cs_real_t *bpro_rusanov = cs_field_by_name("b_rusanov_diff")->val;
  cs_real_6_t *cvar_rij = (cs_real_6_t *)(CS_F_(rij)->val);
  cs_real_t *cvar_rho = CS_F_(rho)->val;

  cs_dispatch_context ctx;

  /* TODO should depend on the model */
  ctx.parallel_for(n_i_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
    const cs_lnum_t c_id0 = i_face_cells[face_id][0];
    const cs_lnum_t c_id1 = i_face_cells[face_id][1];

    /* Note: warning the normal has the surface in it, it is done on purpose */
    cs_real_t r_nn_0 = cs_math_3_sym_33_3_dot_product(i_face_normal[face_id],
                                                      cvar_rij[c_id0],
                                                      i_face_normal[face_id]);
    r_nn_0 *= cs_math_pow2(cvar_rho[c_id0]); // to have rho in it
    cs_real_t r_nn_1 = cs_math_3_sym_33_3_dot_product(i_face_normal[face_id],
                                                      cvar_rij[c_id1],
                                                      i_face_normal[face_id]);
    r_nn_1 *= cs_math_pow2(cvar_rho[c_id1]); // to have rho in it

    cs_real_t rnn = cs_math_fmax(cs_math_fabs(r_nn_0), cs_math_fabs(r_nn_1));

    /* The part of U.n is already in the material upwind scheme */
    ipro_rusanov[face_id] = sqrt(2.0 * rnn);
  });

  ctx.parallel_for(n_b_faces, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
    const cs_lnum_t c_id0 = b_face_cells[face_id];

    /* Note: warning the normal has the surface in it, it is done on purpose */
    cs_real_t r_nn_0 = cs_math_3_sym_33_3_dot_product(b_face_normal[face_id],
                                                      cvar_rij[c_id0],
                                                      b_face_normal[face_id]);
    r_nn_0 *= cs_math_pow2(cvar_rho[c_id0]); // to have rho in it

    /* The part of U.n is already in the material upwind scheme */
    if (   bc_type[face_id] == CS_SMOOTHWALL || bc_type[face_id] == CS_ROUGHWALL
        || bc_type[face_id] == CS_SYMMETRY)
      bpro_rusanov[face_id] = sqrt(2.*cs_math_fabs(r_nn_0));
    else
      bpro_rusanov[face_id] = 0.;

  });

  ctx.wait();
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
