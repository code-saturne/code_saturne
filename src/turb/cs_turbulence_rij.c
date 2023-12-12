/*============================================================================
 * Rij-epsilon turbulence model.
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

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
#include "cs_turbomachinery.h"
#include "cs_turbulence_bc.h"
#include "cs_turbulence_model.h"
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

static const cs_lnum_t _iv2t[6] = {0, 1, 2, 0, 1, 0};
static const cs_lnum_t _jv2t[6] = {0, 1, 2, 1, 2, 2};

static const cs_lnum_t _t2v[3][3] = {{0, 3, 5},
                                     {3, 1, 4},
                                     {5, 4, 2}};

/*============================================================================
 * Private function definitions
 *============================================================================*/

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
  cs_real_t sgn = (b <  0) ? - 1 : 1;

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
  cs_gnum_t nclp[2] =  {0, 0};  /* Min and max clipping values respectively */

  /* Postprocess clippings ? */
  cs_real_t *cpro_a_clipped = NULL;
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
 * \brief Terms of wall echo for \f$ R_{ij} \f$
 *        \f$var  = R_{11} \: R_{22} \: R_{33} \: R_{12} \: R_{13} \: R_{23}\f$
 *        \f$comp =  1 \:  2 \:  3 \:  4 \:  5 \:  6\f$
 *
 * \param[in]    phase_id    turbulent phase id (-1 for single phase flow)
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

  cs_real_t *cromo = NULL;
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
  BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);

  cs_real_t *produk, *epsk;
  BFT_MALLOC(epsk, n_cells_ext, cs_real_t);
  BFT_MALLOC(produk, n_cells_ext, cs_real_t);

  cs_real_6_t *w6;
  BFT_MALLOC(w6, n_cells_ext, cs_real_6_t);

  /* Current gradient */
  cs_field_gradient_scalar(cs_field_by_name("wall_distance"),
                           false,  /* use_previous_t */
                           1,      /* inc */
                           grad);

# pragma omp parallel for if(n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    const cs_real_t norm2 = cs_math_3_square_norm(grad[c_id]);
    const cs_real_t xnorm = cs_math_fmax(sqrt(norm2), cs_math_epzero);

    /* Normalization (warning, the gradient may be sometimes equal to 0) */
    for (cs_lnum_t i = 0; i < 3; i++)
      grad[c_id][i] = -grad[c_id][i] / xnorm;

    /* Production and k */
    produk[c_id] = 0.5 * cs_math_6_trace(produc[c_id]);

    const cs_real_t xk = 0.5 * cs_math_6_trace(cvara_rij[c_id]);
    epsk[c_id] = cvara_ep[c_id] / xk;
  }

  /* Tensor indexing */

  const cs_lnum_t m_33_to_6[3][3] = {{0, 3, 5},
                                     {3, 1, 4},
                                     {5, 4, 2}};

  const cs_lnum_t i_6_to_33[6] = {0, 1, 2, 0, 0, 1};
  const cs_lnum_t j_6_to_33[6] = {0, 1, 2, 1, 2, 2};

  const cs_real_t kr_33[3][3] = {{1., 0., 0.},
                                 {0., 1., 0.},
                                 {0., 0., 1.}};

  const cs_real_t crijp1 = cs_turb_crijp1;
  const cs_real_t crijp2 = cs_turb_crijp2;
  const cs_real_t crij2 = cs_turb_crij2;

  /* Calculation of work variables
   * ----------------------------- */

  cs_array_real_fill_zero(6*n_cells, (cs_real_t *)w6);

# pragma omp parallel for if(n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    for (cs_lnum_t kk = 0; kk < 3; kk++) {

      for (cs_lnum_t mm = 0; mm < 3; mm++) {

        const cs_real_t vnk = grad[c_id][kk];
        const cs_real_t vnm = grad[c_id][mm];

        const cs_lnum_t i_km = m_33_to_6[kk][mm];
        const cs_real_t deltkm = kr_33[kk][mm];

        /* Terms with R km and Phi km;

           vnk * vnm * deltij[jj] in original expression,
           with deltij Kronecker delta so 1 for terms 0-2, 0 for terms 3-6.
           (we only have diagonal contributions) */

        const cs_real_t c_tr
          = vnk * vnm
            * (crijp1 * cvara_rij[c_id][i_km] * epsk[c_id]
               - (  crijp2 * crij2
                  * (produc[c_id][i_km] -d2s3 *produk[c_id] * deltkm)));

        w6[c_id][0] += c_tr;
        w6[c_id][1] += c_tr;
        w6[c_id][2] += c_tr;

      } /* End of loop on mm */

      for (cs_lnum_t isub = 0; isub < 6; isub++) {

        cs_lnum_t ii = i_6_to_33[isub];
        cs_lnum_t jj = j_6_to_33[isub];

        const cs_lnum_t i_ki = m_33_to_6[kk][ii];
        const cs_lnum_t i_kj = m_33_to_6[kk][jj];

        const cs_real_t deltki = kr_33[kk][ii];
        const cs_real_t deltkj = kr_33[kk][jj];

        const cs_real_t vnk = grad[c_id][kk];
        const cs_real_t vni = grad[c_id][ii];
        const cs_real_t vnj = grad[c_id][jj];

        w6[c_id][isub] +=   1.5 * vnk
                          * (- (  crijp1
                                 * (  cvara_rij[c_id][i_ki]*vnj
                                    + cvara_rij[c_id][i_kj]*vni)
                                 * epsk[c_id])
                             + (  crijp2 * crij2
                                * (  ( produc[c_id][i_ki]
                                      -d2s3 * produk[c_id] * deltki)*vnj
                                   + ( produc[c_id][i_kj]
                                      -d2s3 * produk[c_id] * deltkj)*vni)));

      } /* End of loop on isub */

    } /* End of loop on kk */

  } /* End of loop on cells */

  /* Distance to the wall and amortization function: W3
   *   For each calculation mode: same code, test
   *   Apart from the loop */

# pragma omp parallel for if(n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    const cs_real_t distxn = cs_math_fmax(w_dist[c_id], cs_math_epzero); //FIXME
    const cs_real_t trrij = 0.5 * cs_math_6_trace(cvara_rij[c_id]);
    cs_real_t bb =   cmu075 * pow(trrij, 1.5)
                   / (cs_turb_xkappa * cvara_ep[c_id] * distxn);
    bb = cs_math_fmin(bb, 1.0);

    for (cs_lnum_t isub = 0; isub < 6; isub++)
      rhs[c_id][isub] += cromo[c_id] * cell_f_vol[c_id] * w6[c_id][isub] * bb;
  }

  BFT_FREE(w6);
  BFT_FREE(grad);
  BFT_FREE(epsk);
  BFT_FREE(produk);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Gravity terms for terms
 *        For \f$R_{ij}\f$
 *
 * \param[in]   phase_id    turbulent phase id (-1 for single phase flow)
 * \param[in]   gradro      work array for \f$ \grad{\rho} \f$
 * \param[out]  buoyancy    buoyancy term
 */
/*----------------------------------------------------------------------------*/

static void
_gravity_st_rij(int              phase_id,
                const cs_real_t  gradro[][3],
                cs_real_t        buoyancy[][6])
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  cs_field_t *f_rij = CS_F_(rij);
  cs_field_t *f_eps = CS_F_(eps);

  if (phase_id >= 0) {
    f_rij = CS_FI_(rij, phase_id);
    f_eps = CS_FI_(eps, phase_id);
  }

  const cs_real_t *cvara_ep = (const cs_real_t *)f_eps->val_pre;
  const cs_real_6_t *cvara_rij = (const cs_real_6_t *)f_rij->val_pre;

  const cs_real_t uns3 = 1./3;

  cs_real_t cons = -1.5 * cs_turb_cmu;

  const cs_field_t *tf = cs_thermal_model_field();
  if (tf != NULL) {
    const int ksigmas = cs_field_key_id("turbulent_schmidt");
    const cs_real_t turb_schmidt = cs_field_get_key_double(tf, ksigmas);
    cons = -1.5 * cs_turb_cmu / turb_schmidt;
  }

  const cs_real_t *grav = cs_glob_physical_constants->gravity;
  const cs_real_t o_m_crij3 = (1. - cs_turb_crij3);

# pragma omp parallel for if(n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    cs_real_t rit[3];
    cs_math_sym_33_3_product(cvara_rij[c_id], gradro[c_id], rit);

     const cs_real_t k_ov_eps =   cs_math_6_trace(cvara_rij[c_id])
                                / (2*cvara_ep[c_id]);

     cs_real_t gij[3][3];
     for (cs_lnum_t i = 0; i < 3; i++) {
       for (cs_lnum_t j = 0; j < 3; j++)
         gij[i][j] = cons*k_ov_eps* (rit[i]*grav[j] + rit[j]*grav[i]);
     }

     const cs_real_t gkks3 = uns3*(gij[0][0] + gij[1][1] + gij[2][2]);

     buoyancy[c_id][0] = gij[0][0] * o_m_crij3 + cs_turb_crij3*gkks3;
     buoyancy[c_id][1] = gij[1][1] * o_m_crij3 + cs_turb_crij3*gkks3;
     buoyancy[c_id][2] = gij[2][2] * o_m_crij3 + cs_turb_crij3*gkks3;
     buoyancy[c_id][3] = gij[0][1] * o_m_crij3;
     buoyancy[c_id][4] = gij[1][2] * o_m_crij3;
     buoyancy[c_id][5] = gij[0][2] * o_m_crij3;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Gravity terms for \f$\epsilon\f$.
 *
 *  Terms for epsilon:
 *      rom*volumr*deps/dt =
 *                     ... + CEPS1*(EPS/K)*Max(0,(Gkk/2))*volume
 *            With Gij = -(1.5 cmu/PrT) (k/eps) (Rit Gj + Rjt Gi)
 *                 Rit = Rik drom/dxk (sum on k)
 *            We simplify (eps/k) by noting
 *                GijP = -(1.5 cmu/PrT)         (Rit Gj + Rjt Gi)
 *      rom*volume*deps/dt =
 *                     ... + CEPS1*        Max(0,(GkkP/2))*volume
 *
 * \param[in]       phase_id    turbulent phase id (-1 for single phase flow)
 * \param[in]       gradro      density gradient \f$ \grad{\rho} \f$
 * \param[in]       cell_f_vol  cell fluid volume
 * \param[in, out]  rhs         work array for right-hand side
 */
/*----------------------------------------------------------------------------*/

static void
_gravity_st_epsilon(int              phase_id,
                    const cs_real_t  gradro[][3],
                    const cs_real_t  cell_f_vol[],
                    cs_real_t        rhs[])
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  cs_field_t *f_rij = CS_F_(rij);
  if (phase_id >= 0)
    f_rij = CS_FI_(rij, phase_id);

  const cs_real_6_t *cvara_rij = (const cs_real_6_t *)f_rij->val_pre;

  cs_real_t cons = -1.5 * cs_turb_cmu;

  const cs_field_t *tf = cs_thermal_model_field();
  if (tf != NULL) {
    const int ksigmas = cs_field_key_id("turbulent_schmidt");
    const cs_real_t turb_schmidt = cs_field_get_key_double(tf, ksigmas);
    cons = -1.5 * cs_turb_cmu / turb_schmidt;
  }

  const cs_real_t *g = cs_glob_physical_constants->gravity;

# pragma omp parallel for if(n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    const cs_real_t g_rij_gradro
      = cs_math_3_sym_33_3_dot_product(g, cvara_rij[c_id], gradro[c_id]);

    // FIXME for EB-DFM and EBRSM

    const cs_real_t aa = 0.0;
    const cs_real_t bb = cons * g_rij_gradro;

    rhs[c_id] += cs_turb_ce1 * cs_math_fmax(aa, bb) * cell_f_vol[c_id];
  }
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
 * \param[in]   gradro    work array for grad rom
 *                        (without rho volume) only for iturb=30
 * \param[out]  viscf     visc*surface/dist at internal faces
 * \param[out]  viscb     visc*surface/dist at edge faces
 * \param[out]  viscce    Daly Harlow diffusion term
 * \param[out]  rhs       working array
 * \param[out]  rovsdt    working array
 * \param[out]  weighf    working array
 * \param[out]  weighb    working array
 */
/*----------------------------------------------------------------------------*/

static void
_pre_solve_lrr(const cs_field_t  *f_rij,
               int                phase_id,
               const cs_real_t    gradv[][3][3],
               const cs_real_t    produc[][6],
               const cs_real_t    gradro[][3],
               cs_real_t          viscf[],
               cs_real_t          viscb[],
               cs_real_t          viscce[][6],
               cs_real_t          rhs[][6],
               cs_real_t          rovsdt[][6][6],
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

  if (eqp->iwarni >= 1) {
    bft_printf(" Solving the variable %s\n",
               cs_field_get_label(f_rij));
  }

  cs_real_6_t *c_st_prv = NULL;
  int kstprv = cs_field_key_id("source_term_prev_id");
  int st_prv_id = cs_field_get_key_int(f_rij, kstprv);
  if (st_prv_id > -1)
    c_st_prv = (cs_real_6_t *)cs_field_by_id(st_prv_id)->val;

  /* Time extrapolation ? */
  cs_real_t *cromo = NULL;
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

  const cs_real_t deltij[6] = {1, 1, 1, 0, 0, 0};

  cs_lnum_t solid_stride = 1;
  const int c_is_solid_ref[1] = {0};
  int *c_is_solid = cs_solid_zone_flag(cs_glob_mesh);
  if (c_is_solid == NULL) {
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
    for (cs_lnum_t ii = 0; ii < 3; ii++) {
      for (cs_lnum_t jj = 0; jj < 3; jj++) {
        xaniso[ii][jj] =   cvara_var[c_id][_t2v[ii][jj]] / trrij
                         - d2s3 * cs_math_33_identity[ii][jj];
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
    for (cs_lnum_t ii = 0; ii < 6; ii++)
      matrn[ii] = cvara_var[c_id][ii] / trrij;

    cs_real_t oo_matrn[6];
    cs_math_sym_33_inv_cramer(matrn, oo_matrn);
    for (cs_lnum_t ii = 0; ii < 6; ii++)
      oo_matrn[ii] /= trrij;

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
        cs_lnum_t _ij = _t2v[i][j];
        implmat2add[i][j] =   (1.0 - crij2) * xrotac[i][j]
                              + impl_lin_cst * deltij[_ij]
                              + impl_id_cst * d1s2 * oo_matrn[_ij]
                              + ceps_impl * oo_matrn[_ij];
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
      cs_real_t pij = (1.-crij2) * produc[c_id][_t2v[j][i]];
      cs_real_t phiij1 = -cvara_ep[c_id]*crij1*xaniso[j][i];
      cs_real_t phiij2 = d2s3*crij2*trprod*deltij[ij];
      cs_real_t epsij = -d2s3*cvara_ep[c_id]*deltij[ij];

      if (st_prv_id > -1) {
        c_st_prv[c_id][ij] +=   cromo[c_id] * cell_f_vol[c_id]
                              * (pij + phiij1 + phiij2 + epsij);
      }
      else {
        rhs[c_id][ij] +=   cromo[c_id] * cell_f_vol[c_id]
                         * (pij + phiij1 + phiij2 + epsij);

        /* Implicit terms */
        rovsdt[c_id][ij][ij] +=   crom[c_id] * cell_f_vol[c_id] / trrij
                                * (crij1 * cvara_ep[c_id]);

        for (cs_lnum_t jj = 0; jj < 6; jj++)
          rovsdt[c_id][ij][jj] +=   crom[c_id] * cell_f_vol[c_id]
                                  * impl_drsm[ij][jj];
      }
    }

  } /* end loop on cells */

  if (c_is_solid != c_is_solid_ref)
    BFT_FREE(c_is_solid);

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
      for (cs_lnum_t iii = 0; iii < 6; iii++) {
        cs_lnum_t ii = _iv2t[iii];
        cs_lnum_t jj = _jv2t[iii];

        w2[c_id][iii] = 0.;
        for (cs_lnum_t kk = 0; kk < 3; kk++)
          w2[c_id][iii] -=   ccorio
                           * (  matrot[kk][ii] * cvara_var[c_id][_t2v[kk][jj]]
                              + matrot[kk][jj] * cvara_var[c_id][_t2v[kk][ii]]);
      }

      /* Coriolis contribution in the Phi1 term: (1-C2/2)Gij */
      if (icorio == 1) {
        for (cs_lnum_t ii = 0; ii < 6; ii++)
          w2[c_id][ii] =   crom[c_id] * cell_f_vol[c_id]
                         * (1.-0.5*crij2) * w2[c_id][ii];
      }

      /* If source terms are extrapolated */
      if (st_prv_id > -1) {
        for (cs_lnum_t ii = 0; ii < 6; ii++)
          c_st_prv[c_id][ii] += w2[c_id][ii];
      }
      /* Otherwise, directly in rhs */
      else {
        for (cs_lnum_t ii = 0; ii < 6; ii++)
          rhs[c_id][ii] += w2[c_id][ii];
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

  if (cs_glob_turb_rans_model->igrari == 1) {

    cs_real_6_t *_buoyancy = NULL, *cpro_buoyancy = NULL;
    cs_field_t *f_buo = cs_field_by_name_try("rij_buoyancy");

    if (f_buo != NULL) {
      cpro_buoyancy = (cs_real_6_t*)f_buo->val;
    }
    else {
      BFT_MALLOC(_buoyancy, n_cells_ext, cs_real_6_t);
      cpro_buoyancy = _buoyancy;
    }

    _gravity_st_rij(phase_id, gradro, cpro_buoyancy);

    /* If we extrapolate the source terms: previous ST */
    if (st_prv_id > -1) {
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        for (cs_lnum_t ii = 0; ii < 6; ii++)
          c_st_prv[c_id][ii] += cpro_buoyancy[c_id][ii] * cell_f_vol[c_id];
      }
    }
    /* Otherwise RHS */
    else {
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        for (cs_lnum_t ii = 0; ii < 6; ii++)
          rhs[c_id][ii] += cpro_buoyancy[c_id][ii] * cell_f_vol[c_id];
      }
    }

    BFT_FREE(_buoyancy);
  }

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
                                         eqp->iwarni,
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
 * \param[in]   gradro    work array for grad rom
 *                        (without rho volume)
 * \param[out]  viscf     visc*surface/dist at internal faces
 * \param[out]  viscb     visc*surface/dist at edge faces
 * \param[out]  viscce    Daly Harlow diffusion term
 * \param[out]  rhs       working array
 * \param[out]  rovsdt    working array
 * \param[out]  weighf    working array
 * \param[out]  weighb    working array
 */
/*----------------------------------------------------------------------------*/

static void
_pre_solve_lrr_sg(const cs_field_t  *f_rij,
                  int                phase_id,
                  const cs_real_t    produc[][6],
                  const cs_real_t    gradro[][3],
                  cs_real_t          viscf[],
                  cs_real_t          viscb[],
                  cs_real_t          viscce[][6],
                  cs_real_t          rhs[][6],
                  cs_real_t          rovsdt[][6][6],
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

  if (eqp->iwarni >= 1) {
    bft_printf(" Solving the variable %s\n",
               cs_field_get_label(f_rij));
  }

  cs_real_6_t *c_st_prv = NULL;
  int kstprv = cs_field_key_id("source_term_prev_id");
  int st_prv_id = cs_field_get_key_int(f_rij, kstprv);
  if (st_prv_id > -1)
    c_st_prv = (cs_real_6_t *)cs_field_by_id(st_prv_id)->val;

  /* Time extrapolation? */
  cs_real_t  *cromo = NULL;
  int key_t_ext_id = cs_field_key_id("time_extrapolated");
  int iroext = cs_field_get_key_int(f_rho, key_t_ext_id);
  if ((iroext > 0) && (st_prv_id > -1))
    cromo = f_rho->val_pre;
  else
    cromo = f_rho->val;

  const cs_real_t thetv = eqp->thetav;

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

  const cs_real_6_t deltij = {1, 1, 1, 0, 0, 0};

  cs_lnum_t solid_stride = 1;
  const int c_is_solid_ref[1] = {0};
  int *c_is_solid = cs_solid_zone_flag(cs_glob_mesh);
  if (c_is_solid == NULL) {
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

      for (cs_lnum_t ii = 0; ii < 6; ii++) {
        /* Calculation of Prod+Phi1+Phi2-Eps
         *  = rhoPij - C1rho eps/k(Rij-2/3k dij)
         *           - C2rho(Pij-1/3Pkk dij) -2/3rho eps dij
         * In c_st_prv:
         *  = rhoPij-C1rho eps/k(-2/3k dij)-C2rho(Pij-1/3Pkk dij)-2/3rho eps dij
         *  = rho{2/3dij[C2 Pkk/2+(C1-1)eps)]+(1-C2)Pij */
        c_st_prv[c_id][ii] +=  cromo[c_id] * cell_f_vol[c_id]
                              * (  deltij[ii]*d2s3
                                 * (   crij2 * trprod
                                    + (crij1-1.) * cvara_ep[c_id])
                                 + (1.-crij2)*produc[c_id][ii]);

        /*  In rhs = -C1rho eps/k(Rij) = rho {-C1eps/kRij} */
        rhs[c_id][ii] +=   crom[c_id]*cell_f_vol[c_id]
                         * (-crij1 * cvara_ep[c_id] / trrij
                                   * cvara_var[c_id][ii]);

        /* Calculation of the implicit part coming from Phil = C1rho eps/k(1) */
        rovsdt[c_id][ii][ii] +=   crom[c_id] * cell_f_vol[c_id]
                                * crij1 * cvara_ep[c_id] / trrij * thetv;
      }

      /* If we want to implicit a part of -C1rho eps/k(-2/3k dij)
       * FIXME: check if we want to use this or if it should be removed;
       *        previously "isoluc = 2", never called.
       */

      if (false) {
        for (cs_lnum_t ii = 0; ii < 3; ii++) {
          cs_real_t t1 = d1s3 * crij1 * cvara_ep[c_id] / trrij;

          /*  We substract cromo = -C1rho eps/k(-1/3Rij dij) */
          c_st_prv[c_id][ii] -=    cromo[c_id] * cell_f_vol[c_id]
                                 * t1 * cvara_var[c_id][ii];

          /* We add to rhs = -C1rho eps/k(-1/3Rij) */
          rhs[c_id][ii] +=   crom[c_id] * cell_f_vol[c_id]
                           * t1 * cvara_var[c_id][ii];

          /* We add to rovsdt = C1rho eps/k(-1/3) */
          rovsdt[c_id][ii][ii] += crom[c_id] * cell_f_vol[c_id] * t1;
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

      for (cs_lnum_t ii = 0; ii < 6; ii++) {
        /* Calculation of Prod+Phi1+Phi2-Eps
         *  = rhoPij - C1rho eps/k(Rij-2/3k dij)
         *           - C2rho(Pij-1/3Pkk dij) -2/3rho eps dij
         *  = rho{2/3dij[C2 Pkk/2+(C1-1)eps)]+(1-C2)Pij */
        rhs[c_id][ii] +=  crom[c_id] * cell_f_vol[c_id]
                        * (  deltij[ii] * d2s3
                           * (   crij2 * trprod
                              + (crij1-1.) * cvara_ep[c_id])
                           + (1-cs_turb_crij2) * produc[c_id][ii]
                           - crij1 * cvara_ep[c_id]/trrij * cvara_var[c_id][ii]);

        /* Calculation of the implicit part coming from Phi1
         *  = C1rho eps/k(1-1/3 dij) */
        rovsdt[c_id][ii][ii] +=   crom[c_id] * cell_f_vol[c_id]
                                * (1-d1s3 *deltij[ii]) * crij1
                                * cvara_ep[c_id]/trrij;
      }
    }

  }

  if (c_is_solid != c_is_solid_ref)
    BFT_FREE(c_is_solid);

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
      for (cs_lnum_t iii = 0; iii < 6; iii++) {
        cs_lnum_t ii = _iv2t[iii];
        cs_lnum_t jj = _jv2t[iii];

        w2[c_id][iii] = 0.;
        for (cs_lnum_t kk = 0; kk < 3; kk++)
          w2[c_id][iii] -=   ccorio
                           * (  matrot[kk][ii]*cvara_var[c_id][_t2v[kk][jj]]
                              + matrot[kk][jj]*cvara_var[c_id][_t2v[kk][ii]]);
      }

      /* Coriolis contribution in the Phi1 term: (1-C2/2)Gij */
      if (icorio == 1) {
        for (cs_lnum_t ii = 0; ii < 6; ii++)
          w2[c_id][ii] *= crom[c_id] * cell_f_vol[c_id ]* (1.-0.5*cs_turb_crij2);
      }

      /* If source terms are extrapolated */
      if (st_prv_id > -1) {
        for (cs_lnum_t ii = 0; ii < 6; ii++)
          c_st_prv[c_id][ii] += w2[c_id][ii];
      }
      /* Otherwise, directly in RHS */
      else {
        for (cs_lnum_t ii = 0; ii < 6; ii++)
          rhs[c_id][ii] += w2[c_id][ii];
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

  if (cs_glob_turb_rans_model->igrari == 1) {

    cs_real_6_t *_buoyancy = NULL, *cpro_buoyancy = NULL;

    cs_field_t *f_buo = cs_field_by_name_try("rij_buoyancy");
    if (f_buo != NULL) {
      cpro_buoyancy = (cs_real_6_t*)f_buo->val;
    }
    else {
      BFT_MALLOC(_buoyancy, n_cells_ext, cs_real_6_t);
      cpro_buoyancy = _buoyancy;
    }

    _gravity_st_rij(phase_id, gradro, cpro_buoyancy);

    /* If we extrapolate the source terms: previous ST */
    if (st_prv_id > -1) {
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        for (cs_lnum_t ii = 0; ii < 6; ii++)
          c_st_prv[c_id][ii] += cpro_buoyancy[c_id][ii] * cell_f_vol[c_id];
      }
    }
    /* Otherwise RHS */
    else {
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        for (cs_lnum_t ii = 0; ii < 6; ii++)
          rhs[c_id][ii] += cpro_buoyancy[c_id][ii] * cell_f_vol[c_id];
      }
    }

    BFT_FREE(_buoyancy);
  }

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
                                         eqp->iwarni,
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
 * \param[in]   gradro    work array for grad rom
 *                        (without rho volume) only for iturb=30
 * \param[out]  viscf     visc*surface/dist at internal faces
 * \param[out]  viscb     visc*surface/dist at edge faces
 * \param[out]  viscce    Daly Harlow diffusion term
 * \param[out]  rhs       working array
 * \param[out]  rovsdt    working array
 * \param[out]  weighf    working array
 * \param[out]  weighb    working array
 */
/*----------------------------------------------------------------------------*/

static void
_pre_solve_ssg(const cs_field_t  *f_rij,
               int                phase_id,
               const cs_real_t    gradv[][3][3],
               const cs_real_t    produc[][6],
               const cs_real_t    gradro[][3],
               cs_real_t          viscf[],
               cs_real_t          viscb[],
               cs_real_t          viscce[][6],
               cs_real_t          rhs[][6],
               cs_real_t          rovsdt[][6][6],
               cs_real_t          weighf[][2],
               cs_real_t          weighb[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_real_t *cell_f_vol = fvq->cell_f_vol;

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

  cs_real_t *cvar_al = NULL;
  if (cs_glob_turb_model->iturb != CS_TURB_RIJ_EPSILON_SSG)
    cvar_al = f_alpbl->val;

  const cs_equation_param_t *eqp
    = cs_field_get_equation_param_const(f_rij);

  if (eqp->iwarni >= 1) {
    bft_printf(" Solving the variable %s\n ",
               cs_field_get_label(f_rij));
  }

  const int coupled_components = cs_glob_turb_rans_model->irijco;

  cs_real_6_t *c_st_prv = NULL;
  int kstprv = cs_field_key_id("source_term_prev_id");
  int st_prv_id = cs_field_get_key_int(f_rij, kstprv);
  if (st_prv_id > -1)
    c_st_prv = (cs_real_6_t *)cs_field_by_id(st_prv_id)->val;

  /* Time extrapolation ? */
  cs_real_t *cromo = NULL;
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

  const cs_real_t crij3  = cs_turb_crij3;
  const cs_real_t csrij  = cs_turb_csrij;

  const cs_real_t cssgr1 = cs_turb_cssgr1;
  const cs_real_t cssgr2 = cs_turb_cssgr2;
  const cs_real_t cssgr3 = cs_turb_cssgr3;
  const cs_real_t cssgr4 = cs_turb_cssgr4;
  const cs_real_t cssgr5 = cs_turb_cssgr5;

  const cs_real_t cssgs1 = cs_turb_cssgs1;
  const cs_real_t cssgs2 = cs_turb_cssgs2;

  const cs_real_6_t deltij = {1, 1, 1, 0, 0, 0};

  /* Production, Pressure-Strain correlation, dissipation, Coriolis
   * -------------------------------------------------------------- */

  /* Source term
   *  -rho*epsilon*( Cs1*aij + Cs2*(aikajk -1/3*aijaij*deltaij))
   *  -Cr1*P*aij + Cr2*rho*k*sij - Cr3*rho*k*sij*sqrt(aijaij)
   *  +Cr4*rho*k(aik*sjk+ajk*sik-2/3*akl*skl*deltaij)
   *  +Cr5*rho*k*(aik*rjk + ajk*rik)
   *  -2/3*epsilon*deltaij */

  cs_real_3_t *grad_al = NULL;
  if (cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM) {
    BFT_MALLOC(grad_al, n_cells_ext, cs_real_3_t);
    cs_field_gradient_scalar(f_alpbl, true, 1, grad_al);
  }

  const int *irotce = cs_turbomachinery_get_cell_rotor_num();

  cs_real_t cons = -1.5 * cs_turb_cmu;

  cs_field_t *f_thm = cs_thermal_model_field();
  if (f_thm != NULL) {
    double turb_schmidt
      = cs_field_get_key_double(f_thm, cs_field_key_id("turbulent_schmidt"));
    cons = -1.5 * cs_turb_cmu / turb_schmidt;
  }

  cs_real_t *w1, *w2;
  BFT_MALLOC(w1, n_cells_ext, cs_real_t);
  BFT_MALLOC(w2, n_cells_ext, cs_real_t);

  cs_lnum_t solid_stride = 1;
  const int c_is_solid_ref[1] = {0};
  const int *c_is_solid = cs_solid_zone_flag(cs_glob_mesh);
  if (c_is_solid == NULL) {
    c_is_solid = c_is_solid_ref;
    solid_stride = 0;
  }

# pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    if (c_is_solid[solid_stride*c_id])
      continue;

    cs_real_t xnal[3] = {0, 0, 0};

    cs_real_t matrot[3][3];
    cs_real_t impl_drsm[6][6];

    /* EBRSM: compute the magnitude of the Alpha gradient */

    if (cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM) {
      cs_math_3_normalize(grad_al[c_id], xnal);
    }

    cs_real_t xrij[3][3], xprod[3][3];
    cs_real_t xaniso[3][3], xstrai[3][3], xrotac[3][3];

    const cs_real_t trrij  = 0.5 * cs_math_6_trace(cvara_var[c_id]);

    /* Rij */

    for (cs_lnum_t ii = 0; ii < 3; ii++) {
      for (cs_lnum_t jj = 0; jj < 3; jj++) {
        xrij[ii][jj] = cvara_var[c_id][_t2v[ii][jj]];
      }
    }

    /* Pij */
    for (cs_lnum_t ii = 0; ii < 3; ii++) {
      for (cs_lnum_t jj = 0; jj < 3; jj++)
        xprod[ii][jj] = produc[c_id][_t2v[ii][jj]];
    }

    /* aij */
    for (cs_lnum_t ii = 0; ii < 3; ii++) {
      for (cs_lnum_t jj = 0; jj < 3; jj++) {
        xaniso[ii][jj] =   xrij[ii][jj] / trrij
                         - d2s3 * cs_math_33_identity[ii][jj];
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
      const cs_rotation_t *r = cs_glob_rotation + rot_id;

      cs_rotation_coriolis_t(r, 1., matrot);
      for (cs_lnum_t ii = 0; ii < 3; ii++) {
        for (cs_lnum_t jj = ii; jj < 3; jj++) {
          for (cs_lnum_t kk = 0; kk < 3; kk++)
              xprod[jj][ii] -= ccorio * (  matrot[kk][ii] * xrij[kk][jj]
                                         + matrot[kk][jj] * xrij[kk][ii]);
        }
      }

      xprod[0][1] = xprod[1][0]; /* Ensure symmetry (probably not necessary) */
      xprod[0][2] = xprod[2][0];
      xprod[1][2] = xprod[2][1];

    }

    const cs_real_t trprod = 0.5 * (xprod[0][0] + xprod[1][1] + xprod[2][2]);

    /* aii = aijaij */
    cs_real_t aii = 0, aklskl = 0;

    for (cs_lnum_t jj = 0; jj < 3; jj++) {
      for (cs_lnum_t ii = 0; ii < 3; ii++) {
        aii += cs_math_pow2(xaniso[jj][ii]);         /* aij.aij */
        aklskl += xaniso[jj][ii] * xstrai[jj][ii];   /* aij.Sij */
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
      for (cs_lnum_t ii = 0; ii < 6; ii++)
        matrn[ii] = cvara_var[c_id][ii] / trrij;

      cs_real_t oo_matrn[6];
      cs_math_sym_33_inv_cramer(matrn, oo_matrn);
      for (cs_lnum_t ii = 0; ii < 6; ii++)
        oo_matrn[ii] /= trrij;

      cs_real_t impl_lin_cst = 0;
      cs_real_t impl_id_cst = 0;

      /* Compute the maximal eigenvalue (in terms of norm!) of S */
      cs_real_t eigen_vals[3];
      cs_math_sym_33_eigen(sym_strain, eigen_vals);
      cs_real_t eigen_max = cs_math_fabs(eigen_vals[0]);
      for (cs_lnum_t i = 1; i < 3; i++)
        eigen_max = cs_math_fmax(cs_math_fabs(eigen_max),
                                 cs_math_fabs(eigen_vals[i]));

      /* Constant for the dissipation */
      const cs_real_t ceps_impl = d1s3 * cvara_ep[c_id];

      if (cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_SSG) {

        /* Identity constant for phi3 */
        const cs_real_t cphi3impl = cs_math_fabs(cssgr2 - cssgr3*sqrt(aii));

        /* Identity constant */
        impl_id_cst = - d2s3 * cssgr1 * cs_math_fmin(trprod, 0)
                      - d1s3 * cssgs2 * cvara_ep[c_id] * aii
                      + cphi3impl * trrij * eigen_max
                      + 2. * d2s3 * cssgr4 * trrij * eigen_max
                      + d2s3 * trrij * cssgr4 * cs_math_fmax(aklskl, 0);

        /* Linear constant */
        impl_lin_cst = eigen_max * (1. + cssgr4 + cssgr5);

        cs_real_t implmat2add[3][3];
        for (cs_lnum_t i = 0; i < 3; i++) {
          for (cs_lnum_t j = 0; j < 3; j++) {
            const cs_lnum_t _ij = _t2v[i][j];
            implmat2add[i][j] =   xrotac[i][j]
                                + impl_lin_cst * deltij[_ij]
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
        const cs_real_t cphiw_impl = 6 * (1-alpha3) * cvara_ep[c_id] / trrij;

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
            const cs_lnum_t _ij = _t2v[i][j];
            implmat2add[i][j] =   xrotac[i][j]
                                + impl_lin_cst * deltij[_ij]
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
      cs_lnum_t i = _iv2t[ij];
      cs_lnum_t j = _jv2t[ij];

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

      if (cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_SSG) {

        /* Explicit terms */
        const cs_real_t pij =     xprod[j][i];
        const cs_real_t phiij1 =   -cvara_ep[c_id]
                                 * (  cssgs1 * xaniso[j][i]
                                    + cssgs2 * (aikakj - d1s3*deltij[ij]*aii));
        const cs_real_t phiij2 =   -cssgr1 * trprod * xaniso[j][i]
                                 +   trrij * xstrai[j][i]
                                   * (cssgr2 - cssgr3*sqrt(aii))
                                 + cssgr4*trrij * (  aiksjk
                                                   - d2s3 * deltij[ij] * aklskl)
                                 + cssgr5*trrij * aikrjk;
        const cs_real_t epsij = -d2s3 * cvara_ep[c_id] * deltij[ij];

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
        const cs_real_t xnnd = d1s2 * (xnal[i]*xnal[j] + deltij[ij]);

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
                         + cebmr4 * (aiksjk - d2s3 * deltij[ij] * aklskl)
                         + cebmr5 * aikrjk);

        /* Compute \f $\e_{ij}^w \f$ (Rotta model)
         * Rij/k*epsilon */
        const cs_real_t epsijw = xrij[j][i] / trrij * cvara_ep[c_id];

        /* Compute \e_{ij}^h */
        const cs_real_t epsij = d2s3 * cvara_ep[c_id] * deltij[ij];

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
          epsijw_imp = 6. * (1.-alpha3) * cvara_ep[c_id] / trrij;

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
        rovsdt[c_id][ij][ij] += w2[c_id];

        if (coupled_components != 0) {
          for (cs_lnum_t jj = 0; jj < 6; jj++)
            rovsdt[c_id][ij][jj] +=   crom[c_id] * cell_f_vol[c_id]
                                    * impl_drsm[ij][jj];
        }
      }
    } /* End of loop on ij */

  } /* end loop on cells */

  BFT_FREE(grad_al);

  BFT_FREE(w2);

  /* Buoyancy source term
   * -------------------- */

  if (cs_glob_turb_rans_model->igrari == 1) {

    cs_real_6_t *_buoyancy = NULL, *cpro_buoyancy = NULL;
    cs_field_t *f_buo = cs_field_by_name_try("rij_buoyancy");

    if (f_buo != NULL) {
      cpro_buoyancy = (cs_real_6_t*)f_buo->val;
    }
    else {
      BFT_MALLOC(_buoyancy, n_cells_ext, cs_real_6_t);
      cpro_buoyancy = _buoyancy;
    }

    _gravity_st_rij(phase_id, gradro, cpro_buoyancy);

    /* If we extrapolate the source terms: previous ST */
    if (st_prv_id > -1) {
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        for (cs_lnum_t ii = 0; ii < 6; ii++)
          c_st_prv[c_id][ii] += cpro_buoyancy[c_id][ii] * cell_f_vol[c_id];
    }
    /* Otherwise RHS */
    else {
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        for (cs_lnum_t ii = 0; ii < 6; ii++) {
          rhs[c_id][ii] += cpro_buoyancy[c_id][ii] * cell_f_vol[c_id];
        }
      }
    }

    BFT_FREE(_buoyancy);

    if (coupled_components != 0) {

      const cs_real_t *grav = cs_glob_physical_constants->gravity;

#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

        if (c_is_solid[solid_stride*c_id])
          continue;

        const cs_real_t trrij = 0.5 * cs_math_6_trace(cvara_var[c_id]);
        const cs_real_t k_ov_eps = trrij / cvara_ep[c_id];

        cs_real_t implmat2add[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

        cs_real_t gkks3 = 0;
        for (cs_lnum_t jj = 0; jj < 3; jj++) {
          for (cs_lnum_t ii = 0; ii < 3; ii++)
            gkks3 += grav[ii] * gradro[c_id][jj] * cvara_var[c_id][_t2v[jj][ii]];
        }
        gkks3 *= cons * k_ov_eps * crij3 * d2s3;

        if (gkks3 <= 0) {
          /* Term "C3 tr(G) Id"
             Compute inverse matrix of R^n
             (scaling by tr(R) for numerical stability) */

          cs_real_t matrn[6];
          for (int ii = 0; ii < 6; ii++)
            matrn[ii] = cvara_var[c_id][ii] / trrij;

          cs_real_t oo_matrn[6];
          cs_math_sym_33_inv_cramer(matrn, oo_matrn);
          for (int ii = 0; ii < 6; ii++)
            oo_matrn[ii] /= trrij;

          for (cs_lnum_t jj = 0; jj < 3; jj++) {
            for (cs_lnum_t ii = 0; ii < 3; ii++) {
              cs_lnum_t iii = _t2v[jj][ii];
              implmat2add[jj][ii] = -0.5 * gkks3 * oo_matrn[iii];
            }
          }
        }

        const cs_real_t gradchk = cs_math_3_dot_product(grav, gradro[c_id]);

        if (gradchk > 0) {
          /* Implicit term written as:
           *   Po . R^n+1 + R^n+1 . Po^t
           * with Po proportional to "g (x) Grad rho" */
          const cs_real_t gradro_impl = cons * (1.-crij3) * k_ov_eps;
          for (cs_lnum_t i = 0; i < 3; i++) {
            for (cs_lnum_t j = 0; j < 3; j++)
              implmat2add[i][j] -= gradro_impl * grav[i] * gradro[c_id][j];
          }
        }

        /* Compute the 6x6 matrix A which verifies
         * A.R = M.R + R.M^t */
        cs_real_t impl_drsm[6][6];
        cs_math_reduce_sym_prod_33_to_66(implmat2add, impl_drsm);

        for (cs_lnum_t ii = 0; ii < 6; ii++) {
          for (cs_lnum_t jj = 0; jj < 6; jj++)
            rovsdt[c_id][ii][jj] += cell_f_vol[c_id] * impl_drsm[ii][jj];
        }

      } /* End of loop on cells */

    } /* End of test on coupled components */

  } /* End for buoyancy source term */

  if (c_is_solid != c_is_solid_ref)
    BFT_FREE(c_is_solid);

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
                                         eqp->iwarni,
                                         weighf,
                                         weighb,
                                         viscf,
                                         viscb);

  }

  /* Scalar diffusivity */
 else {

    if (eqp->idifft == 1) {
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        w1[c_id] = viscl[c_id] + (csrij * visct[c_id] / cs_turb_cmu);
    }
    else
      cs_array_real_copy(n_cells, viscl, w1);

    cs_face_viscosity(m,
                      fvq,
                      eqp->imvisf,
                      w1,
                      viscf,
                      viscb);
  }

  BFT_FREE(w1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve epsilon for the \f$ R_{ij} - \varepsilon \f$ RANS
 *        turbulence model.
 *
 * \param[in]     phase_id    turbulent phase id (-1 for single phase flow)
 * \param[in]     ncesmp      number of cells with mass source term
 * \param[in]     icetsm      index of cells with mass source term
 * \param[in]     itypsm      type of mass source term for each variable
 * \param[in]     gradv       work array for the term grad
 *                            of velocity only for iturb=31
 * \param[in]     produc      work array for production (without
 *                            rho volume) only for iturb=30
 * \param[in]     gradro      work array for \f$ \grad{rom} \f$
 * \param[in]     smacel      value associated to each variable in the mass
 *                            source terms or mass rate
 * \param[in]     viscf       visc*surface/dist at internal faces
 * \param[in]     viscb       visc*surface/dist at edge faces
 * \param[in]     rhs         working array
 * \param[in]     rovsdt      working array
 !*/
/*----------------------------------------------------------------------------*/

static void
_solve_epsilon(int              phase_id,
               cs_lnum_t        ncesmp,
               cs_lnum_t        icetsm[],
               int              itypsm[],
               const cs_real_t  gradv[][3][3],
               const cs_real_t  produc[][6],
               const cs_real_t  gradro[][3],
               cs_real_t        smacel[],
               cs_real_t        viscf[],
               cs_real_t        viscb[],
               cs_real_t        rhs[],
               cs_real_t        rovsdt[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

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

  const cs_real_t *cvar_al = NULL;
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

  if (eqp->iwarni >= 1) {
    bft_printf(" Solving the variable %s\n",
               cs_field_get_label(f_eps));
  }

  int kstprv = cs_field_key_id("source_term_prev_id");
  int st_prv_id = cs_field_get_key_int(f_eps, kstprv);
  cs_real_t *c_st_prv = NULL, *cromo = NULL;
  if (st_prv_id > -1)
    c_st_prv = cs_field_by_id(st_prv_id)->val;

  int key_t_ext_id = cs_field_key_id("time_extrapolated");
  int iroext = cs_field_get_key_int(f_rho, key_t_ext_id);
  if ((iroext > 0) && (st_prv_id > -1))
    cromo = f_rho->val_pre;
  else
    cromo = f_rho->val;

  /* S as Source, V as Variable */
  const cs_time_scheme_t *time_scheme = cs_glob_time_scheme;
  const cs_real_t thets = time_scheme->thetst;
  const cs_real_t thetv = eqp->thetav;

  cs_array_real_fill_zero(n_cells, rhs);
  cs_array_real_fill_zero(n_cells, rovsdt);

  /* Work arrays */
  cs_real_t *w1;
  BFT_MALLOC(w1, n_cells_ext, cs_real_t);

  /* User source terms
   * ----------------- */

  cs_user_source_terms(cs_glob_domain,
                       f_eps->id,
                       rhs,
                       rovsdt);

  if (cs_glob_porous_model == 3)
    cs_immersed_boundary_wall_functions(f_eps->id, rhs, rovsdt);

  /* If we extrapolate the source terms */
  if (st_prv_id > -1) {
#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      /* Save for exchange */
      const cs_real_t tuexpe = c_st_prv[c_id];
      /* For the continuation and the next time step */
      c_st_prv[c_id] = rhs[c_id];
      /* RHS of previous time step
       *   We assume -rovsdt > 0: we implicit
       *   the user source term (the rest)  */
      rhs[c_id] = rovsdt[c_id]*cvara_ep[c_id] - thets*tuexpe;
      /* Diagonal */
      rovsdt[c_id] = -thetv*rovsdt[c_id];
    }
  }
  else {
#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      rhs[c_id]   += rovsdt[c_id]*cvara_ep[c_id];
      rovsdt[c_id]  = cs_math_fmax(-rovsdt[c_id], cs_math_zero_threshold);
    }
  }

  /* Lagrangian source terms
   * ----------------------- */

  /* Second order is not taken into account  */

  if (   cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING
      && (cs_glob_lagr_source_terms->ltsdyn == 1)) {

    const cs_real_6_t *lagr_st_rij
      = (const cs_real_6_t *)cs_field_by_name("rij_st_lagr")->val;

#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      /* Source terms with epsilon */
      const cs_real_t st_eps = -0.50 * cs_math_6_trace(lagr_st_rij[c_id]);

      /* k */
      const cs_real_t k = 0.50 * cs_math_6_trace(cvara_rij[c_id]);

      /* equiv:       cs_turb_ce4 * st_eps * eps / (k / eps) */
      rhs[c_id]   += cs_turb_ce4 * st_eps / k;

      /* equiv:                    -cs_turb_ce4 * st_eps * / (k/eps) */
      rovsdt[c_id] += cs_math_fmax(-cs_turb_ce4 * st_eps / k * cvara_ep[c_id],
                                   cs_math_zero_threshold);
   }
  }

  /* Mass source term
   * ---------------- */

  if (ncesmp > 0) {

    const int var_key_id = cs_field_key_id("variable_id");
    const int ivar_eps = cs_field_get_key_int(f_eps, var_key_id)-1;
    const int ivar_pr = cs_field_get_key_int(CS_F_(p), var_key_id)-1;

    /* We increment rhs with -Gamma.var_prev. and rovsdt with Gamma */

    cs_mass_source_terms(1, /* iterns*/
                         1, /* dim */
                         ncesmp,
                         icetsm,
                         itypsm + ncesmp*ivar_eps,
                         cell_f_vol,
                         cvara_ep,
                         smacel + ncesmp*ivar_eps,
                         smacel + ncesmp*ivar_pr,
                         rhs,
                         rovsdt,
                         w1);

    /* If we extrapolate the source terms, we put Gamma Pinj in c_st_prv */
    if (st_prv_id > -1) {
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        c_st_prv[c_id] += w1[c_id];
    }

    /* Otherwise we put it directly in rhs */
    else {
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        rhs[c_id] += w1[c_id];
    }
  }

  /* Unsteady term
   * ------------- */

  if (eqp->istat == 1) {
#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      rovsdt[c_id] += (crom[c_id] / dt[c_id]) * cell_f_vol[c_id];
    }
  }

  /* Production (rho * Ce1 * epsilon / k * P)
   *    Dissipation (rho*Ce2.epsilon/k*epsilon)
   * ------------------------------------------ */

  cs_real_t thetap = (st_prv_id > -1) ? thetv : 1.;

  cs_real_t ceps2;
  if (cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_LRR) {
    ceps2 = cs_turb_ce2;
  }
  else if (cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_SSG) {
    ceps2 = cs_turb_cssge2;
  }
  else {
    ceps2 = cs_turb_cebme2;
  }

  cs_real_t *cprod;
  BFT_MALLOC(cprod, n_cells_ext, cs_real_t);

  /* Calculation the production trace, depending we are in standard
   * Rij or in SSG (use of produc or grdvit) */

  if (cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_LRR) {
#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      cprod[c_id] = 0.5*(produc[c_id][0] + produc[c_id][1] + produc[c_id][2]);
  }
  else {
#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      cprod[c_id] = -(  cvara_rij[c_id][0] * gradv[c_id][0][0]
                      + cvara_rij[c_id][3] * gradv[c_id][0][1]
                      + cvara_rij[c_id][5] * gradv[c_id][0][2]
                      + cvara_rij[c_id][3] * gradv[c_id][1][0]
                      + cvara_rij[c_id][1] * gradv[c_id][1][1]
                      + cvara_rij[c_id][4] * gradv[c_id][1][2]
                      + cvara_rij[c_id][5] * gradv[c_id][2][0]
                      + cvara_rij[c_id][4] * gradv[c_id][2][1]
                      + cvara_rij[c_id][2] * gradv[c_id][2][2]);
  }

  /* EBRSM */
  if (cs_glob_turb_model->iturb ==  CS_TURB_RIJ_EPSILON_EBRSM) {

#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      /* Half-traces */
      const cs_real_t trprod = cprod[c_id];
      const cs_real_t trrij = 0.5 * cs_math_6_trace(cvara_rij[c_id]);

      /* Calculation of the Durbin time scale */
      const cs_real_t xttke = trrij / cvara_ep[c_id];
      const cs_real_t xttkmg
        = cs_turb_xct*sqrt(viscl[c_id] / crom[c_id] / cvara_ep[c_id]);
      const cs_real_t xttdrb = cs_math_fmax(xttke, xttkmg);

      const cs_real_t prdeps = trprod / cvara_ep[c_id];
      const cs_real_t alpha3 = cs_math_pow3(cvar_al[c_id]);

      /* Production (explicit) */
      /* Compute of C_eps_1' */
      w1[c_id] =   cromo[c_id] * cell_f_vol[c_id] * cs_turb_ce1
                 * (1+cs_turb_xa1*(1-alpha3)*prdeps) * trprod / xttdrb;

      /* Dissipation (implicit) */
      const cs_real_t crom_vol = crom[c_id] * cell_f_vol[c_id];
      rhs[c_id] -= crom_vol * ceps2 * cvara_ep[c_id] / xttdrb;
      rovsdt[c_id] += ceps2 * crom_vol * thetap / xttdrb;
    }

  }

  /* SSG and LRR */
  else {

#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      /* Half-traces */
      const cs_real_t trprod = cprod[c_id];
      const cs_real_t trrij = 0.5 * cs_math_6_trace(cvara_rij[c_id]);
      const cs_real_t xttke = trrij / cvara_ep[c_id];

      /* Production (explicit, a part might be implicit) */
      const cs_real_t cromo_vol = cromo[c_id] * cell_f_vol[c_id];
      rovsdt[c_id]
        += cs_math_fmax(- cromo_vol * cs_turb_ce1 * trprod / trrij,
                        0.0);
      w1[c_id] = cromo_vol * cs_turb_ce1 / xttke * trprod;

      /* Dissipation (implicit) */
      const cs_real_t crom_vol = crom[c_id] * cell_f_vol[c_id];
      rhs[c_id]   -= crom_vol * ceps2 * cs_math_pow2(cvara_ep[c_id]) / trrij;
      rovsdt[c_id] += ceps2 * crom_vol / xttke * thetap;
    }

  }

  /* Extrapolation of source terms (2nd order in time) */
  if (st_prv_id > -1) {
#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      c_st_prv[c_id] += w1[c_id];
  }
  else {
#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      rhs[c_id] += w1[c_id];
  }

  /* Buoyancy term
   * ------------- */

  /* FIXME use beta ... WARNING */
  if (cs_glob_turb_rans_model->igrari == 1) {

    /* Extrapolation of source terms (2nd order in time) */
    if (st_prv_id > -1)
      _gravity_st_epsilon(phase_id, gradro, cell_f_vol, c_st_prv);

    else
      _gravity_st_epsilon(phase_id, gradro, cell_f_vol, rhs);

  }

  /* Diffusion term (Daly Harlow: generalized gradient hypothesis method)
   * -------------------------------------------------------------------- */

  cs_real_t *weighb;
  cs_real_6_t *viscce;
  cs_real_2_t *weighf;

  BFT_MALLOC(weighb, n_b_faces, cs_real_t);
  BFT_MALLOC(weighf, n_i_faces, cs_real_2_t);
  BFT_MALLOC(viscce, n_cells_ext, cs_real_6_t);

  /* Symmetric tensor diffusivity (GGDH) */
  if (eqp->idften & CS_ANISOTROPIC_DIFFUSION) {

    const cs_field_t *f_a_t_visc
      = cs_field_by_name("anisotropic_turbulent_viscosity");
    const cs_real_6_t *visten = (const cs_real_6_t *)f_a_t_visc->val;

#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      for (cs_lnum_t i = 0; i < 3; i++)
        viscce[c_id][i] = visten[c_id][i] / sigmae + viscl[c_id];
      for (cs_lnum_t i = 3; i < 6; i++)
        viscce[c_id][i] = visten[c_id][i] / sigmae;
    }

    cs_face_anisotropic_viscosity_scalar(m,
                                         fvq,
                                         viscce,
                                         eqp->iwarni,
                                         weighf,
                                         weighb,
                                         viscf,
                                         viscb);

  }

  /* Scalar diffusivity */
  else {

    if (eqp->idifft == 1) {
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        w1[c_id] = viscl[c_id] + visct[c_id]/sigmae;
    }
    else
      cs_array_real_copy(n_cells, viscl, w1);

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
#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      rhs[c_id] += thetp1*c_st_prv[c_id];
  }

  cs_solid_zone_set_zero_on_cells(1, rhs);

  /* Get boundary conditions coefficients */

  cs_real_t *coefap = f_eps->bc_coeffs->a;
  cs_real_t *coefbp = f_eps->bc_coeffs->b;
  cs_real_t *cofafp = f_eps->bc_coeffs->af;
  cs_real_t *cofbfp = f_eps->bc_coeffs->bf;

  cs_equation_param_t eqp_loc = *eqp;

  eqp_loc.istat  = -1;
  eqp_loc.icoupl = -1;
  eqp_loc.idifft = -1;
  eqp_loc.iwgrec = 0;     /* Warning, may be overwritten if a field */
  eqp_loc.thetav = thetv;
  eqp_loc.blend_st = 0;   /* Warning, may be overwritten if a field */

  cs_real_t *dpvar;
  BFT_MALLOC(dpvar, n_cells_ext, cs_real_t);

  cs_equation_iterative_solve_scalar(cs_glob_time_step_options->idtvar,
                                     1,   /* init */
                                     f_eps->id,
                                     NULL,
                                     0,   /* iescap */
                                     0,   /* imucpp */
                                     -1,  /* normp  */
                                     &eqp_loc,
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
                                     viscce,
                                     weighf,
                                     weighb,
                                     0,  /* boundary convective upwind flux */
                                     NULL,
                                     rovsdt,
                                     rhs,
                                     cvar_ep,
                                     dpvar,
                                     NULL,
                                     NULL);

  /* Free memory */

  BFT_FREE(dpvar);
  BFT_FREE(w1);
  BFT_FREE(cprod);
  BFT_FREE(viscce);
  BFT_FREE(weighb);
  BFT_FREE(weighf);
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
 * \param[in]     ncesmp       number of cells with mass source term
 * \param[in]     icetsm       index of cells with mass source term
 * \param[in]     itypsm       mass source type for the variables
 * \param[in]     smacel       values of the variables associated to the
 *                             mass source
 *                             (for ivar=ipr, smacel is the mass flux)
 !*/
/*-----------------------------------------------------------------------------*/

void
cs_turbulence_rij(int          phase_id,
                  cs_lnum_t    ncesmp,
                  cs_lnum_t    icetsm[],
                  int          itypsm[],
                  cs_real_t    smacel[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_turb_model_t *turb_model = cs_glob_turb_model;
  const cs_turb_rans_model_t *turb_rans_model = cs_glob_turb_rans_model;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_real_t *cell_f_vol = fvq->cell_f_vol;

  cs_field_t *f_rij = CS_F_(rij);
  cs_field_t *f_eps = CS_F_(eps);
  cs_field_t *f_alpbl = CS_F_(alp_bl);
  cs_field_t *f_vel = CS_F_(vel);
  cs_field_t *f_rho = CS_F_(rho);
  cs_field_t *f_rhob = CS_F_(rho_b);

  if (phase_id >= 0) {
    f_rij = CS_FI_(rij, phase_id);
    f_eps = CS_FI_(eps, phase_id);
    f_alpbl = CS_FI_(alp_bl, phase_id);
    f_vel = CS_FI_(vel, phase_id);
    f_rho = CS_FI_(rho, phase_id);
    f_rhob = CS_FI_(rho_b, phase_id);
  }

  cs_real_t *cvar_ep = f_eps->val;
  cs_real_6_t *cvar_rij = (cs_real_6_t *)f_rij->val;

  const cs_real_t *dt = CS_F_(dt)->val;
  const cs_real_t *crom = f_rho->val;
  const cs_real_6_t *cvara_rij = (const cs_real_6_t *)f_rij->val_pre;

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  const int iflmas =  cs_field_get_key_int(f_vel, kimasf);
  const int iflmab =  cs_field_get_key_int(f_vel, kbmasf);

  const cs_real_t *imasfl = cs_field_by_id(iflmas)->val;
  const cs_real_t *bmasfl = cs_field_by_id(iflmab)->val;

  /* Time extrapolation ? */
  cs_real_6_t *c_st_prv = NULL;
  int kstprv = cs_field_key_id("source_term_prev_id");
  int st_prv_id = cs_field_get_key_int(f_rij, kstprv);
  if (st_prv_id > -1)
    c_st_prv = (cs_real_6_t *)cs_field_by_id(st_prv_id)->val;

  const cs_equation_param_t *eqp
    = cs_field_get_equation_param_const(f_rij);

  const cs_time_step_t *time_step = cs_glob_time_step;
  const cs_time_scheme_t *time_scheme = cs_glob_time_scheme;

  const cs_real_t thets = time_scheme->thetst;
  const cs_real_t thetv = eqp->thetav;
  const cs_real_t ro0 = cs_glob_fluid_properties->ro0;
  const cs_real_t viscl0 = cs_glob_fluid_properties->viscl0;

  if (eqp->iwarni >= 1) {
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

  cs_real_3_t *gradro = NULL;

  cs_real_6_t *cpro_press_correl = NULL;
  {
    cs_field_t *f_psc = cs_field_by_name_try("rij_pressure_strain_correlation");
    if (f_psc != NULL)
      cpro_press_correl = (cs_real_6_t *)f_psc->val;
  }

  cs_real_6_t *produc = NULL, *_produc = NULL;
  if (cs_field_by_name_try("rij_production") != NULL) {
    produc = (cs_real_6_t *)cs_field_by_name_try("rij_production")->val;
  }
  else {
    BFT_MALLOC(_produc, n_cells_ext, cs_real_6_t);
    produc = _produc;
  }

  cs_real_t *viscf, *viscb;
  BFT_MALLOC(viscf, n_i_faces, cs_real_t);
  BFT_MALLOC(viscb, n_b_faces, cs_real_t);

  cs_real_33_t *gradv = NULL, *_gradv = NULL;
  {
    cs_field_t *f_vg = cs_field_by_name_try("velocity_gradient");

    if (f_vel->grad != NULL)
      gradv = (cs_real_33_t *)f_vel->grad;
    else if (f_vg != NULL)
      gradv = (cs_real_33_t *)f_vg->val;
    else {
      BFT_MALLOC(_gradv, n_cells_ext, cs_real_33_t);
      gradv = _gradv;
    }
  }

  cs_real_6_t *smbrts;
  cs_real_66_t *rovsdtts;
  BFT_MALLOC(smbrts, n_cells_ext, cs_real_6_t);
  BFT_MALLOC(rovsdtts,  n_cells_ext, cs_real_66_t);

  /* Advanced initialiation for EBRSM
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
    BFT_MALLOC(grad, n_cells_ext, cs_real_3_t);

    /* Compute the gradient of Alpha */
    cs_field_gradient_scalar(f_alpbl,
                             false,  /* use_previous_t */
                             1,      /* inc */
                             grad);

#   pragma omp parallel for if(n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      /* Velocity magnitude */
      cs_real_t xunorm = cs_math_3_norm(vel[c_id]);

      /* y+ is bounded by 400, because in the Reichard profile,
       * it corresponds to saturation (u>uref) */
      cvar_al[c_id] = cs_math_fmax(cs_math_fmin(cvar_al[c_id], (1.-exp(-8.))),
                                   0.);
      /* Magnitude and unit vector of the Alpha gradient */
      cs_real_t xnal[3];
      cs_real_t xnoral = cs_math_3_norm(grad[c_id]);
      if (xnoral <= cs_math_epzero / pow(cell_f_vol[c_id], cs_math_1ov3)) {
        for (int ii = 0; ii < 3; ii++)
          xnal[ii] = 1.0 / sqrt(3.0);
      }
      else {
        for (int ii = 0; ii < 3; ii++)
          xnal[ii] = grad[c_id][ii]/xnoral;
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

      for (cs_lnum_t ii = 0; ii < 3; ii++)
        vel[c_id][ii] = limiter * vel[c_id][ii];

      const cs_real_t ut2 = 0.050 * uref;
      cvar_ep[c_id] =   cs_math_pow3(utaurf)
                      * cs_math_fmin(1. / (xkappa * 15.0 * nu0 / utaurf),
                                     1. / (xkappa * ya));
      const cs_real_t tke =     cvar_ep[c_id] * 0.5 / nu0*cs_math_pow2(ya)
                              * cs_math_pow2(exp(-ypa/25))
                            +   cs_math_pow2(ut2) / 0.3
                              * cs_math_pow2((1.-exp(-ypa/25.)));

      for (int ii = 0; ii < 3; ii++)
        cvar_rij[c_id][ii] =      alpha3  * 2./3 * tke
                            + (1.-alpha3) * (1. - cs_math_pow2(xnal[ii])) * tke;
      cvar_rij[c_id][3] = -(1.-alpha3)*(xnal[0]*xnal[1])*tke;
      cvar_rij[c_id][4] = -(1.-alpha3)*(xnal[1]*xnal[2])*tke;
      cvar_rij[c_id][5] = -(1.-alpha3)*(xnal[0]*xnal[2])*tke;

    } /* End of loop on cells */

    BFT_FREE(grad);

    cs_field_current_to_previous(f_vel);
    cs_field_current_to_previous(f_rij);
  }

  /* Compute the velocity gradient */

  if (f_vel->grad == NULL)
    cs_field_gradient_vector(f_vel, true, 1, gradv);

  /* Compute the production term for Rij
   * ----------------------------------- */

# pragma omp parallel for if(n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    /* Pij = - (Rik dUk/dXj + dUk/dXi Rkj)
     * Pij is stored as (P11, P22, P33, P12, P23, P13) */
    produc[c_id][0] = -2. * (  cvara_rij[c_id][0] * gradv[c_id][0][0]
                             + cvara_rij[c_id][3] * gradv[c_id][0][1]
                             + cvara_rij[c_id][5] * gradv[c_id][0][2]);

    produc[c_id][1] = -2. * (  cvara_rij[c_id][3] * gradv[c_id][1][0]
                             + cvara_rij[c_id][1] * gradv[c_id][1][1]
                             + cvara_rij[c_id][4] * gradv[c_id][1][2]);

    produc[c_id][2] = -2. * (  cvara_rij[c_id][5] * gradv[c_id][2][0]
                             + cvara_rij[c_id][4] * gradv[c_id][2][1]
                             + cvara_rij[c_id][2] * gradv[c_id][2][2]);

    produc[c_id][3] = - (  cvara_rij[c_id][3] * gradv[c_id][0][0]
                         + cvara_rij[c_id][1] * gradv[c_id][0][1]
                         + cvara_rij[c_id][4] * gradv[c_id][0][2])
                      - (  cvara_rij[c_id][0] * gradv[c_id][1][0]
                         + cvara_rij[c_id][3] * gradv[c_id][1][1]
                         + cvara_rij[c_id][5] * gradv[c_id][1][2]);

    produc[c_id][4] = - (  cvara_rij[c_id][5] * gradv[c_id][1][0]
                         + cvara_rij[c_id][4] * gradv[c_id][1][1]
                         + cvara_rij[c_id][2] * gradv[c_id][1][2])
                      - (  cvara_rij[c_id][3] * gradv[c_id][2][0]
                         + cvara_rij[c_id][1] * gradv[c_id][2][1]
                         + cvara_rij[c_id][4] * gradv[c_id][2][2]);

    produc[c_id][5] = - (  cvara_rij[c_id][5] * gradv[c_id][0][0]
                         + cvara_rij[c_id][4] * gradv[c_id][0][1]
                         + cvara_rij[c_id][2] * gradv[c_id][0][2])
                      - (  cvara_rij[c_id][0] * gradv[c_id][2][0]
                         + cvara_rij[c_id][3] * gradv[c_id][2][1]
                         + cvara_rij[c_id][5] * gradv[c_id][2][2]);

  }

  /* Compute the pressure correlation  term for Rij
   * ----------------------------------------------
   * Phi,ij = Phi1,ij+Phi2,ij
   * Phi,ij = -C1 k/eps (Rij-2/3k dij) - C2 (Pij-2/3P dij)
   * Phi,ij is stored as (Phi11, Phi22, Phi33, Phi12, Phi23, Phi13)

   * TODO : coherency with the model */

  if (cpro_press_correl != NULL)  {
    const cs_real_t d2s3 = 2./3;

    const cs_real_t crij1 = cs_turb_crij1;
    const cs_real_t crij2 = cs_turb_crij2;

#   pragma omp parallel for if(n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      const cs_real_t k = 0.5 * cs_math_6_trace(cvara_rij[c_id]);
      const cs_real_t p = 0.5 * cs_math_6_trace(produc[c_id]);
      for (cs_lnum_t ii = 0; ii < 3; ii++)
        cpro_press_correl[c_id][ii]
          = - crij1 * cvar_ep[c_id] / k * (cvara_rij[c_id][ii] - d2s3*k)
            - crij2 * (produc[c_id][ii] - d2s3 * p);
      for (cs_lnum_t ii = 3; ii < 6; ii++)
        cpro_press_correl[c_id][ii]
          = - crij1 * cvar_ep[c_id] /  k * (cvara_rij[c_id][ii])
            - crij2 * (produc[c_id][ii]);
    }
  }

  /* Compute the density gradient for buoyant terms
   * ---------------------------------------------- */

  /* Note that the buoyant term is normally expressed in temr of
   * (u'T') or (u'rho') here modelled with a GGDH:
   *   (u'rho') = C * k/eps * R_ij Grad_j(rho)

   * Buoyant term for the Atmospheric module
   * (function of the potential temperature) */

  if (turb_rans_model->igrari == 1) {

    if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 1) {

      const cs_field_t *thf = cs_thermal_model_field();
      if (thf != NULL) {

        const cs_real_t *cvara_scalt = thf->val_pre;

        /* TODO check if we really want cromo (as below)
           or simply crom) */
        const cs_real_t *cromo = f_rho->val;

        BFT_MALLOC(gradro, n_cells_ext, cs_real_3_t);

        cs_field_gradient_scalar(thf,
                                 true,   /* use_previous_t */
                                 1,      /* inc */
                                 gradro);

        /* gradro stores: - rho grad(theta)/theta
         * grad(rho) and grad(theta) have opposite signs */

#       pragma omp parallel for if(n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          for (cs_lnum_t ii = 0; ii < 3; ii++)
            gradro[c_id][ii] *= -cromo[c_id] / cvara_scalt[c_id];
        }
      }

    }
    else { /* Other models */

      const cs_velocity_pressure_model_t  *vp_model
        = cs_glob_velocity_pressure_model;
      const int idilat = vp_model->idilat;

      BFT_MALLOC(gradro, n_cells_ext, cs_real_3_t);

      if (idilat == 0) {
        const cs_field_t *thf = cs_thermal_model_field();
        if (phase_id >= 0)
          thf = CS_FI_(h_tot, phase_id);

        if (thf != NULL) {

          cs_field_gradient_scalar(thf,
                                   false,  /* use current (not previous) value */
                                   1,      /* inc */
                                   gradro);

          const cs_real_t *cpro_beta
            = cs_field_by_name("thermal_expansion")->val;

          /* gradro stores: - beta grad(T)
           * grad(rho) and grad(T) have opposite signs */

#         pragma omp parallel for if(n_cells > CS_THR_MIN)
          for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
            for (cs_lnum_t ii = 0; ii < 3; ii++)
              gradro[c_id][ii] = -ro0 * cpro_beta[c_id] * gradro[c_id][ii];
          }

        }
      }
      else {

        /* Boundary conditions: Dirichlet romb
         * We use viscb to store the relative coefficient of rom
         * We impose in Dirichlet (coefa) the value romb */

        cs_array_real_fill_zero(n_b_faces, viscb);

        int key_t_ext_id = cs_field_key_id("time_extrapolated");
        int iroext = cs_field_get_key_int(f_rho, key_t_ext_id);

        int t_i = (cs_glob_time_scheme->isto2t > 0 && iroext > 0) ? 1 : 0;
        cs_real_t *cromo = f_rho->vals[t_i];
        const cs_real_t *bromo = f_rhob->vals[t_i];

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
                           eqp->iwarni,
                           eqp->imligr,
                           eqp->epsrgr,
                           eqp->climgr,
                           NULL,          /* f_ext */
                           bromo,
                           viscb,
                           cromo,
                           NULL,         /* c_weight */
                           NULL,         /* cpl */
                           gradro);

      }

    }

  } /* End of test on igrari */

  /* Prepare to solve Rij, in a manner similar
     to that of cs_solve_equation_scalar.
   * ========================================= */

  /* Source terms for Rij
   * -------------------- */

  cs_array_real_fill_zero(6*n_cells, (cs_real_t*)smbrts);
  cs_array_real_fill_zero(36*n_cells, (cs_real_t*)rovsdtts);

  cs_user_source_terms(cs_glob_domain,
                       f_rij->id,
                       (cs_real_t*)smbrts,
                       (cs_real_t*)rovsdtts);

  if (cs_glob_porous_model == 3)
    cs_immersed_boundary_wall_functions(f_rij->id,
                                        (cs_real_t*)smbrts,
                                        (cs_real_t*)rovsdtts);

  if (c_st_prv != NULL) {

#   pragma omp parallel for if(n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      for (cs_lnum_t ii = 0; ii < 6; ii++) {
        const cs_real_t tuexpr = c_st_prv[c_id][ii];
        /* For continuation and the next time step */
        c_st_prv[c_id][ii] = smbrts[c_id][ii];

        smbrts[c_id][ii] = -thets*tuexpr;
        /* Right hand side of the previous time step
         * We suppose -rovsdt > 0: we implicit
         * the user source term (the rest) */
        for (cs_lnum_t jj = 0; jj < 6; jj++) {
          smbrts[c_id][ii] += rovsdtts[c_id][ii][jj] * cvara_rij[c_id][jj];
          /* Diagonal */
          rovsdtts[c_id][ii][jj] *= -thetv;
        }
      }
    }

  }
  else {

#   pragma omp parallel for if(n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      for (cs_lnum_t ii = 0; ii < 6; ii++) {
        for (cs_lnum_t jj = 0; jj < 6; jj++)
          smbrts[c_id][ii] += rovsdtts[c_id][ii][jj] * cvara_rij[c_id][jj];
        rovsdtts[c_id][ii][ii] = cs_math_fmax(-rovsdtts[c_id][ii][ii], 0.);
      }
    }

  }

  /* Lagrangian source terms
   * ----------------------- */

  /* Second order is not taken into account  */

  if (   cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING
      && cs_glob_lagr_source_terms->ltsdyn == 1) {
    const cs_real_6_t *lagr_st_rij
      = (const cs_real_6_t *)cs_field_by_name_try("rij_st_lagr")->val;

    const cs_lagr_source_terms_t  *lag_st = cs_glob_lagr_source_terms;
    cs_real_t *tslagi = lag_st->st_val + (lag_st->itsli-1)*n_cells_ext;

#   pragma omp parallel for if(n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      for (cs_lnum_t ii = 0; ii < 6; ii++) {
        smbrts[c_id][ii] += lagr_st_rij[c_id][ii];
        rovsdtts[c_id][ii][ii] += cs_math_fmax(-tslagi[c_id],
                                               cs_math_zero_threshold);
      }
    }
  }

  /* Mass source terms
   *------------------ */

  if (ncesmp > 0) {
    const int var_key_id = cs_field_key_id("variable_id");
    int ivar = cs_field_get_key_int(f_rij, var_key_id)-1;
    int ivar_pr = cs_field_get_key_int(CS_F_(p), var_key_id)-1;

    cs_real_6_t *gatinj;
    BFT_MALLOC(gatinj, n_cells_ext, cs_real_6_t);

    /* We increment smbrts with -Gamma.var_prev. and rovsdr with Gamma */
    cs_mass_source_terms(1, /* iterns*/
                         6, /* dim */
                         ncesmp,
                         icetsm,
                         itypsm + ivar,
                         cell_f_vol,
                         (cs_real_t*)cvara_rij,
                         smacel + ivar,
                         smacel + ivar_pr,
                         (cs_real_t*)smbrts,
                         (cs_real_t*)rovsdtts,
                         (cs_real_t*)gatinj);

    /* If we extrapolate the source terms, we put Gamma Pinj in c_st_prv */
    if (st_prv_id > -1) {
#     pragma omp parallel for if(n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        for (cs_lnum_t ii = 0; ii < 6; ii++)
          c_st_prv[c_id][ii] += gatinj[c_id][ii];
    }
    /* Otherwise we put it directly in smbrts */
    else {
#     pragma omp parallel for if(n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        for (cs_lnum_t ii = 0; ii < 6; ii++)
          smbrts[c_id][ii] += gatinj[c_id][ii];
    }

    BFT_FREE(gatinj);
  }

  /* Unsteady term
   * ------------- */
  if (eqp->istat == 1) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      for (int ii = 0; ii < 6; ii++)
        rovsdtts[c_id][ii][ii] += (crom[c_id]/dt[c_id])*cell_f_vol[c_id];
    }
  }

  /* Terms specific to Rij-epsilon model
   * ----------------------------------- */

  cs_real_t *weighb;
  cs_real_6_t *viscce;
  cs_real_2_t *weighf;

  BFT_MALLOC(weighb, n_b_faces, cs_real_t);
  BFT_MALLOC(weighf, n_i_faces, cs_real_2_t);
  BFT_MALLOC(viscce, n_cells_ext, cs_real_6_t);

   if (turb_model->iturb == CS_TURB_RIJ_EPSILON_LRR) {
     if (turb_rans_model->irijco == 1)
       _pre_solve_lrr(f_rij, phase_id, gradv,
                      produc, gradro,
                      viscf, viscb, viscce,
                      smbrts, rovsdtts,
                      weighf, weighb);
     else
       _pre_solve_lrr_sg(f_rij, phase_id,
                         produc, gradro,
                         viscf, viscb, viscce,
                         smbrts, rovsdtts,
                         weighf, weighb);

   }
   else { /* if (   turb_model->iturb == CS_TURB_RIJ_EPSILON_SSG
                 || turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM) */
     _pre_solve_ssg(f_rij, phase_id, gradv,
                    produc, gradro,
                    viscf, viscb, viscce,
                    smbrts, rovsdtts,
                    weighf, weighb);
   }

   cs_real_6_t  *coefap = (cs_real_6_t *)f_rij->bc_coeffs->a;
   cs_real_66_t *coefbp = (cs_real_66_t *)f_rij->bc_coeffs->b;
   cs_real_6_t  *cofafp = (cs_real_6_t *)f_rij->bc_coeffs->af;
   cs_real_66_t *cofbfp = (cs_real_66_t *)f_rij->bc_coeffs->bf;

   /* Add Rusanov fluxes */
   if (cs_glob_turb_rans_model->irijnu == 2) {
     cs_real_t *ipro_rusanov = cs_field_by_name("i_rusanov_diff")->val;
     for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
       viscf[face_id] = fmax(0.5 * ipro_rusanov[face_id], viscf[face_id]);
     }

     const cs_real_3_t *restrict b_face_normal
       = (const cs_real_3_t *restrict)fvq->b_face_normal;
     cs_real_t *b_lam = cs_field_by_name("b_rusanov_diff")->val;

     for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
       cs_real_t n[3];
       cs_math_3_normalize(b_face_normal[face_id], n);
       cs_real_66_t bf;
       const cs_real_t kr_33[3][3] = {{1., 0., 0.},
                                      {0., 1., 0.},
                                      {0., 0., 1.}};

       for (cs_lnum_t ij = 0; ij < 6; ij++) {
         cs_lnum_t i = _iv2t[ij];
         cs_lnum_t j = _jv2t[ij];
         for (cs_lnum_t kl = 0; kl < 6; kl++) {
           cs_lnum_t k = _iv2t[kl];
           cs_lnum_t l = _jv2t[kl];
           bf[ij][kl] = b_lam[face_id] * n[l] *(
                 n[i] * (kr_33[j][k] - n[j] * n[k])
               + n[j] * (kr_33[i][k] - n[i] * n[k])
               );
         }
       }

       for (cs_lnum_t i = 0; i < 6; i++) {
         for (cs_lnum_t j = 0; j < 6; j++) {
           cofbfp[face_id][i][j] +=  bf[i][j];
           cofafp[face_id][i] -= bf[i][j] * coefap[face_id][j];
         }
       }

     }

   }

   /* Solve Rij
    * --------- */

   if (c_st_prv != NULL) {
     const cs_real_t thetp1 = 1. + thets;
#    pragma omp parallel for if(n_cells > CS_THR_MIN)
     for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
       for (cs_lnum_t ii = 0; ii < 6; ii++)
         smbrts[c_id][ii] += thetp1 * c_st_prv[c_id][ii];
   }

   cs_solid_zone_set_zero_on_cells(6, (cs_real_t *)smbrts);

   /* All boundary convective flux with upwind */
   int icvflb = 0;

   cs_equation_param_t eqp_loc = *eqp;
   eqp_loc.istat  = -1;
   eqp_loc.iwgrec = 0;     /* Warning, may be overwritten if a field */
   eqp_loc.thetav = thetv;
   eqp_loc.blend_st = 0;   /* Warning, may be overwritten if a field */

   cs_equation_iterative_solve_tensor(cs_glob_time_step_options->idtvar,
                                      f_rij->id,
                                      NULL,
                                      &eqp_loc,
                                      cvara_rij,
                                      cvara_rij,
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
                                      viscce,
                                      weighf,
                                      weighb,
                                      icvflb,
                                      NULL,
                                      rovsdtts,
                                      smbrts,
                                      cvar_rij);

   BFT_FREE(viscce);
   BFT_FREE(rovsdtts);
   BFT_FREE(smbrts);

   /* Solve Epsilon
    * ------------- */

   {
     cs_real_t *smbr, *rovsdt;
     BFT_MALLOC(smbr, n_cells_ext, cs_real_t);
     BFT_MALLOC(rovsdt, n_cells_ext, cs_real_t);

     _solve_epsilon(phase_id,
                    ncesmp,
                    icetsm,
                    itypsm,
                    gradv,
                    produc,
                    gradro,
                    smacel,
                    viscf,
                    viscb,
                    smbr,
                    rovsdt);

     BFT_FREE(rovsdt);
     BFT_FREE(smbr);
   }

   /* Clipping
    * -------- */

   cs_turbulence_rij_clip(phase_id, n_cells);

   /* Free memory */

   gradv = NULL;
   BFT_FREE(_gradv);
   produc = NULL;
   BFT_FREE(_produc);

   BFT_FREE(gradro);

   BFT_FREE(viscf);
   BFT_FREE(viscb);
   BFT_FREE(gradv);
   BFT_FREE(weighb);
   BFT_FREE(weighf);
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

  if (eqp->iwarni == 1) {
    cs_log_printf(CS_LOG_DEFAULT,
                  _(" Solving the variable %s\n"),
                  cs_field_get_label(cs_field_by_id(f_id)));
  }

  cs_real_t *rovsdt, *rhs;

  /* Allocate temporary arrays */
  BFT_MALLOC(rhs, n_cells_ext, cs_real_t);
  BFT_MALLOC(rovsdt, n_cells_ext, cs_real_t);

  cs_array_real_fill_zero(n_cells, rhs);
  cs_array_real_fill_zero(n_cells, rovsdt);

  /* Source term of alpha
   *  \f$ rhs = \dfrac{1}{L^2 (\alpha)} - \dfrac{1}{L^2}\f$
   * In fact there is a mark "-" because the solved equation is
   *   \f$-\div{\grad {alpha}} = rhs \f$
   *================================================================*/

  /* Matrix */

  const cs_real_t thetv = eqp->thetav;
  cs_real_t thetap = (cs_glob_time_scheme->isto2t > 0) ? thetv : 1.0;

  // FIXME the source term extrapolation is not well done!!!!

  /* For automatic initialization, the length scale is fixed at L^+ =50 */

  if (   cs_glob_time_step->nt_cur == 1
      && cs_glob_turb_rans_model->reinit_turb == 1) {

    const cs_real_t ro0 = cs_glob_fluid_properties->ro0;
    const cs_real_t viscl0 = cs_glob_fluid_properties->viscl0;
    const cs_real_t xlldrb = 50.0 * viscl0 / ro0 / (0.050*uref);
    const cs_real_t l2 = cs_math_pow2(xlldrb);

#   pragma omp parallel for if(n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      /* Explicit term */
      rhs[c_id] = cell_f_vol[c_id]*(1.0-cvara_al[c_id]) / l2;

      /* Implicit term */
      rovsdt[c_id] = (rovsdt[c_id]+cell_f_vol[c_id]*thetap) / l2;
    }

  }
  else {

#   pragma omp parallel for if(n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      const cs_real_t xk
        = d1s2 * (cvara_rij[c_id][0] + cvara_rij[c_id][1] + cvara_rij[c_id][2]);
      const cs_real_t xnu = viscl[c_id] / crom[c_id];

      /* Integral length scale */
      const cs_real_t xllke = pow(xk, d3s2) / cvara_ep[c_id];

      /* Kolmogorov length scale */
      const cs_real_t xllkmg =   cs_turb_xceta
                               * pow(cs_math_pow3(xnu)/cvara_ep[c_id], d1s4);

      /* Durbin length scale */
      const cs_real_t xlldrb = c_durbin_l*cs_math_fmax(xllke, xllkmg);

      const cs_real_t l2 = cs_math_pow2(xlldrb);

      /* Explicit term */
      rhs[c_id] = cell_f_vol[c_id]*(1.0-cvara_al[c_id]) / l2;

      /* Implicit term */
      rovsdt[c_id] = (rovsdt[c_id]+cell_f_vol[c_id]*thetap) / l2;
    }

  }

  /* Calculation of viscf and viscb for cs_equation_iterative_solve_scalar. */
  cs_real_t *w1, *viscf, *viscb;
  BFT_MALLOC(w1, n_cells_ext, cs_real_t);
  BFT_MALLOC(viscf, n_i_faces, cs_real_t);
  BFT_MALLOC(viscb, n_b_faces, cs_real_t);

  cs_array_real_set_scalar(n_cells, 1., w1);

  cs_face_viscosity(m,
                    fvq,
                    eqp->imvisf,
                    w1,
                    viscf,
                    viscb);

  BFT_FREE(w1);

  cs_solid_zone_set_zero_on_cells(1, rhs);

  /* Effective resolution of the equation of alpha
     ============================================= */

  const cs_real_t *coefap = cs_field_by_id(f_id)->bc_coeffs->a;
  const cs_real_t *coefbp = cs_field_by_id(f_id)->bc_coeffs->b;
  const cs_real_t *cofafp = cs_field_by_id(f_id)->bc_coeffs->af;
  const cs_real_t *cofbfp = cs_field_by_id(f_id)->bc_coeffs->bf;

  cs_equation_param_t eqp_loc = *eqp;

  eqp_loc.istat  = -1;
  eqp_loc.icoupl = -1;
  eqp_loc.idifft = -1;
  eqp_loc.iwgrec = 0;     /* Warning, may be overwritten if a field */
  eqp_loc.thetav = thetv;
  eqp_loc.blend_st = 0;   /* Warning, may be overwritten if a field */

  cs_real_t *dpvar;
  BFT_MALLOC(dpvar, n_cells_ext, cs_real_t);

  cs_equation_iterative_solve_scalar(cs_glob_time_step_options->idtvar,
                                     1, /* init */
                                     f_id,
                                     NULL,
                                     0, /* iescap */
                                     0, /* imucpp */
                                     -1, /* normp */
                                     &eqp_loc,
                                     cvara_al,
                                     cvara_al,
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
                                     rovsdt,
                                     rhs,
                                     cvar_al,
                                     dpvar,
                                     NULL,
                                     NULL);

  BFT_FREE(dpvar);
  BFT_FREE(rhs);

  cs_solid_zone_set_scalar_on_cells(1., cvar_al);

  /* Clipping
     ======== */

  cs_real_t *alpha_min;
  BFT_MALLOC(alpha_min, n_cells_ext, cs_real_t);

  /* Compute a first estimator of the minimal value of alpha per cell.
   * This is deduced from "alpha/L^2 - div(grad alpha) = 1/L^2" and assuming that
   * boundary cell values are 0. This value is thefore non zero but
   * much smaller than the wanted value. */

  cs_array_real_fill_zero(n_cells_ext, alpha_min);
  cs_array_real_copy(n_cells, rovsdt, alpha_min);

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
    const cs_lnum_t ii = i_face_cells[face_id][0];
    const cs_lnum_t jj = i_face_cells[face_id][1];
    alpha_min[ii] += viscf[face_id];
    alpha_min[jj] += viscf[face_id];
  }

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    const cs_lnum_t ii = b_face_cells[face_id];
    alpha_min[ii] += viscb[face_id]/distb[face_id];
  }

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    alpha_min[c_id] = rovsdt[c_id]/alpha_min[c_id];

  _clip_alpha(f_id, n_cells, alpha_min);

  BFT_FREE(alpha_min);
  BFT_FREE(rovsdt);
  BFT_FREE(viscf);
  BFT_FREE(viscb);
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

  cs_real_t *cvar_ep = CS_F_(eps)->val;
  cs_real_6_t *cvar_rij = (cs_real_6_t *)CS_F_(rij)->val;

  /* With reference velocity */

  if (uref > 0) {
    const cs_real_t tr_ii = cs_math_pow2(0.02 * uref);
    const cs_real_t k = 0.5 * (3. * tr_ii);  /* trace of tensor with
                                                tr_ii diagonal) */
    const cs_real_t ep = pow(k, 1.5) * cs_turb_cmu / almax;

#   pragma omp parallel if(n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      cvar_rij[c_id][0] = tr_ii;
      cvar_rij[c_id][1] = tr_ii;
      cvar_rij[c_id][2] = tr_ii;
      cvar_rij[c_id][3] = 0;
      cvar_rij[c_id][4] = 0;
      cvar_rij[c_id][5] = 0;
      cvar_ep[c_id] = ep;
    }

    cs_turbulence_rij_clip(-1, n_cells);
  }

  /* Without reference velocity */
  else {
#   pragma omp parallel if(n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      for (cs_lnum_t i = 0; i < 6; i++)
        cvar_rij[c_id][i] = -cs_math_big_r;
      cvar_ep[c_id] =  -cs_math_big_r;
    }
  }

  cs_solid_zone_set_zero_on_cells(6, (cs_real_t *)cvar_rij);
  cs_solid_zone_set_scalar_on_cells(1e-12, cvar_ep);
  //cs_solid_zone_set_zero_on_cells(1, cvar_ep);

  /* For EBRSM, initialize alpha */

  if (cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM) {
    cs_real_t *cvar_al = CS_F_(alp_bl)->val;
    cs_array_real_set_scalar(n_cells, 1., cvar_al);
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

  cs_real_t *cpro_eps_clipped = NULL;
  cs_real_6_t *cpro_rij_clipped = NULL;
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
  for (cs_lnum_t ii = 0; ii < 3; ii++) {
    if (vmax[ii] > rijmax[ii])
      rijmax[ii] = vmax[ii];
  }

  cs_parall_max(3, CS_REAL_TYPE, rijmax);

  const cs_real_t trref = rijmax[0] + rijmax[1] + rijmax[2];
  const cs_real_t rijref = cs_math_fmax(trref/3., cs_math_epzero);

  cs_lnum_t solid_stride = 1;
  const int c_is_solid_ref[1] = {0};
  int *c_is_solid = cs_solid_zone_flag(cs_glob_mesh);
  if (c_is_solid == NULL) {
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
        for (cs_lnum_t ii = 0; ii < 6; ii++)
          cvar_rij[c_id][ii] = 0;
        cvar_ep[c_id] = 1e-12;
        continue;
      }

      /* Check if R is positive and ill-conditioned (since the former
       * will induce the latter after clipping ...*/

      const cs_real_t trrij = cs_math_6_trace(cvar_rij[c_id]);

      if (trrij <= cs_math_epzero*trref) {
        for (cs_lnum_t ii = 0; ii < 3; ii++) {
          if (cpro_rij_clipped != NULL) {
            cpro_rij_clipped[c_id][ii]
              = cvar_rij[c_id][ii] - cs_math_epzero*rijref;
            cpro_rij_clipped[c_id][ii+3] = cvar_rij[c_id][ii+3];
          }

          cvar_rij[c_id][ii] = cs_math_epzero*rijref;
          cvar_rij[c_id][ii+3] = 0.0;

          t_iclrij[ii]++;
          t_iclrij[ii+3]++;
        }

        is_clipped = 1;
      }
      else {
        cs_real_t tensor[6];
        for (cs_lnum_t ii = 0; ii < 6; ii++)
          tensor[ii] = cvar_rij[c_id][ii]/trrij;

        cs_real_t eigen_vals[3];
        cs_math_sym_33_eigen(tensor, eigen_vals);

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

          for (cs_lnum_t ii = 0; ii < 6; ii++) {
            cvar_rij[c_id][ii] = (1.0-eigen_offset)*cvar_rij[c_id][ii];

            if (ii < 3)
              cvar_rij[c_id][ii] += trrij*(eigen_offset+eigen_tol)/3.;

            if (cpro_rij_clipped != NULL)
              cpro_rij_clipped[c_id][ii] = eigen_offset*cvar_rij[c_id][ii];

            t_iclrij[ii]++;
          }
        }
      }

      /* Epsilon */

      if (cs_math_fabs(cvar_ep[c_id]) < epz2) {
        t_iclep[0]++;
        if (cpro_eps_clipped != NULL)
          cpro_eps_clipped[c_id] = cs_math_fabs(cvar_ep[c_id]-epz2);
        cvar_ep[c_id] = cs_math_fmax(cvar_ep[c_id],epz2);
      }
      else if (cvar_ep[c_id] <= 0) {
        t_iclep[0]++;
        if (cpro_eps_clipped != NULL)
          cpro_eps_clipped[c_id] = 2*cs_math_fabs(cvar_ep[c_id]);
        cvar_ep[c_id] = cs_math_fmin(cs_math_fabs(cvar_ep[c_id]),
                                     varrel*cs_math_fabs(cvara_ep[c_id]));
      }

      /* Enforce Cauchy Schwartz inequality (only for x, y, z directions) */

      cs_real_t cvar_var1, cvar_var2;
      for (cs_lnum_t ii = 3; ii < 6; ii++) {
        if (ii == 3) {
          cvar_var1 = cvar_rij[c_id][0];
          cvar_var2 = cvar_rij[c_id][1];
        }
        else if (ii == 4) {
          cvar_var1 = cvar_rij[c_id][1];
          cvar_var2 = cvar_rij[c_id][2];
        }
        else if (ii == 5) {
          cvar_var1 = cvar_rij[c_id][0];
          cvar_var2 = cvar_rij[c_id][2];
        }

        const cs_real_t rijmin = sqrt(cvar_var1*cvar_var2);
        if (rijmin < cs_math_fabs(cvar_rij[c_id][ii])) {
          is_clipped = 1;
          if (cpro_rij_clipped != NULL)
            cpro_rij_clipped[c_id][ii] = cvar_rij[c_id][ii];
          cvar_rij[c_id][ii] =   _sign(1., cvar_rij[c_id][ii])
                               * rijmin/(1.+cs_math_epzero);
          t_iclrij[ii]++;
        }
      }
      t_icltot += is_clipped;

    }  /* End of loop on cells */

    /* Sum over threads */

    for (cs_lnum_t i = 0; i < 6; i++) {
      #pragma omp atomic
      iclrij[i] += t_iclrij[i];
    }

    #pragma omp atomic
    iclep[0] += t_iclep[0];

    #pragma omp atomic
    icltot += t_icltot;
  }

  if (c_is_solid != c_is_solid_ref)
    BFT_FREE(c_is_solid);

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

  /* EBRSM case */

  if (cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM) {

    cs_real_3_t *grad_al = NULL;
    BFT_MALLOC(grad_al, n_cells_ext, cs_real_3_t);
    cs_field_gradient_scalar(f_alpbl, true, 1, grad_al);

    const cs_real_t *cvar_al = f_alpbl->val;

#   pragma omp parallel if(n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
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

      cs_real_t xk = 0.5 * (xrij[0][0] + xrij[1][1] + xrij[2][2]);
      cs_real_t xe = cvar_ep[c_id];

      /* We compute the normal Reynolds Stresses */

      cs_real_t xrnn = 0;
      for (cs_lnum_t ii = 0; ii < 3; ii++) {
        for (cs_lnum_t jj = 0; jj < 3; jj++)
          xrnn += xrij[ii][jj]*xnal[jj]*xnal[ii];
      }
      xrnn = (1.-alpha3)*xrnn + alpha3*xk;
      xrnn = cs_math_fmax(xrnn, 1.e-12);

      visct[c_id] = crom[c_id] * cs_turb_cmu * xrnn * xk / xe;
    }

    BFT_FREE(grad_al);

  }

  /* SSG and LRR */

  else {

#   pragma omp parallel if(n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      cs_real_t xk = 0.5 * cs_math_6_trace(cvar_rij[c_id]);
      cs_real_t xrnn = cs_math_fmax(xk, 1.e-12);
      cs_real_t xe = cvar_ep[c_id];

      visct[c_id] = crom[c_id] * cs_turb_cmu * xrnn * xk / xe;
    }

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
    = (const cs_real_3_t *restrict)mq->i_face_normal;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)mq->b_face_normal;
  const cs_lnum_2_t *i_face_cells = m->i_face_cells;
  const cs_lnum_t *b_face_cells = m->b_face_cells;
  const int *bc_type = cs_glob_bc_type;

  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_i_faces = m->n_i_faces;

  cs_real_t *ipro_rusanov = cs_field_by_name("i_rusanov_diff")->val;
  cs_real_t *bpro_rusanov = cs_field_by_name("b_rusanov_diff")->val;
  cs_real_6_t *cvar_rij = (cs_real_6_t *)(CS_F_(rij)->val);

  /* TODO should depend on the model */
  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
    const cs_lnum_t c_id0 = i_face_cells[face_id][0];
    const cs_lnum_t c_id1 = i_face_cells[face_id][1];

    /* Note: warning the normal has the surface in it, it is done on purpose */
    cs_real_t r_nn_0 = cs_math_3_sym_33_3_dot_product(i_face_normal[face_id],
                                                      cvar_rij[c_id0],
                                                      i_face_normal[face_id]);
    r_nn_0 *= cs_math_pow2(CS_F_(rho)->val[c_id0]); // to have rho in it
    cs_real_t r_nn_1 = cs_math_3_sym_33_3_dot_product(i_face_normal[face_id],
                                                      cvar_rij[c_id1],
                                                      i_face_normal[face_id]);
    r_nn_1 *= cs_math_pow2(CS_F_(rho)->val[c_id1]); // to have rho in it

    cs_real_t rnn = fmax(fabs(r_nn_0), fabs(r_nn_1));

    /* The part of U.n is already in the material upwind scheme */
    ipro_rusanov[face_id] = sqrt(2.0 * rnn);
  }

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    const cs_lnum_t c_id0 = b_face_cells[face_id];
    cs_real_3_t n;
    /* Warning normalized */
    cs_math_3_normalize(b_face_normal[face_id], n);

    /* Note: warning the normal has the surface in it, it is done on purpose */
    cs_real_t r_nn_0 = cs_math_3_sym_33_3_dot_product(b_face_normal[face_id],
                                                      cvar_rij[c_id0],
                                                      b_face_normal[face_id]);
    r_nn_0 *= cs_math_pow2(CS_F_(rho)->val[c_id0]); // to have rho in it

    /* The part of U.n is already in the material upwind scheme */
    if (bc_type[face_id] == CS_SMOOTHWALL ||bc_type[face_id] == CS_ROUGHWALL
        || bc_type[face_id] == CS_SYMMETRY)
      bpro_rusanov[face_id] = sqrt(2.*fabs(r_nn_0));
    else
      bpro_rusanov[face_id] = 0.;

  }

}

/*----------------------------------------------------------------------------*/

 END_C_DECLS
