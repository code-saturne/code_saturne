/*============================================================================
 * Rij-epsilon turbulence model.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
#include "cs_domain.h"
#include "cs_equation_iterative_solve.h"
#include "cs_equation_param.h"
#include "cs_face_viscosity.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_lagr.h"
#include "cs_log_iteration.h"
#include "cs_mass_source_terms.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_prototypes.h"
#include "cs_thermal_model.h"
#include "cs_time_step.h"
#include "cs_turbulence_bc.h"
#include "cs_turbulence_model.h"

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
 *  Global variables
 *============================================================================*/

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
    cs_array_set_value_real(n_cells, 1, 0, cpro_a_clipped);
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
 * \param[in]       gradro      density gradient \f$ \grad{\rho} \f$
 * \param[in]       cell_f_vol  cell fluid volume
 * \param[in, out]  rhs         work array for right-hand side
 */
/*----------------------------------------------------------------------------*/

static void
_gravity_st_epsilon(const cs_real_t  gradro[][3],
                    const cs_real_t  cell_f_vol[],
                    cs_real_t        rhs[])
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  const cs_real_6_t *cvara_rij = (const cs_real_6_t *)CS_F_(rij)->val_pre;

  cs_real_t cons = -1.5*cs_turb_cmu;

  const cs_field_t *tf = cs_thermal_model_field();
  if (tf != NULL) {
    const int ksigmas = cs_field_key_id("turbulent_schmidt");
    const cs_real_t turb_schmidt = cs_field_get_key_double(tf, ksigmas);
    cons = -1.5*cs_turb_cmu/turb_schmidt;
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

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve of epsilon for \f$ R_{ij} - \varepsilon \f$ RANS
 *        turbulence model.
 *
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
/*-------------------------------------------------------------------------------*/

void
cs_turbulence_rij_solve_eps(cs_lnum_t            ncesmp,
                            cs_lnum_t            icetsm[],
                            int                  itypsm[],
                            const cs_real_33_t   gradv[],
                            const cs_real_6_t    produc[],
                            const cs_real_3_t    gradro[],
                            cs_real_t            smacel[],
                            cs_real_t            viscf[],
                            cs_real_t            viscb[],
                            cs_real_t            rhs[],
                            cs_real_t            rovsdt[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_real_t *cell_f_vol = fvq->cell_f_vol;

  const cs_real_t *dt = CS_F_(dt)->val;
  const cs_real_t *crom = CS_F_(rho)->val;
  const cs_real_t *viscl = CS_F_(mu)->val;
  const cs_real_t *visct = CS_F_(mu_t)->val;
  const cs_real_t *cvara_ep = CS_F_(eps)->val_pre;
  const cs_real_6_t *cvara_rij = (const cs_real_6_t *)CS_F_(rij)->val_pre;

  const cs_real_t *cvar_al = NULL;
  if (cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM)
    cvar_al = (const cs_real_t *)(CS_F_(alp_bl)->val);

  cs_real_t *cvar_ep = CS_F_(eps)->val;

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  int iflmas =  cs_field_get_key_int(CS_F_(vel), kimasf);
  int iflmab =  cs_field_get_key_int(CS_F_(vel), kbmasf);

  const cs_real_t *imasfl = cs_field_by_id(iflmas)->val;
  const cs_real_t *bmasfl = cs_field_by_id(iflmab)->val;

  const int ksigmas = cs_field_key_id("turbulent_schmidt");
  const cs_real_t sigmae = cs_field_get_key_double(CS_F_(eps), ksigmas);

  const cs_equation_param_t *eqp
    = cs_field_get_equation_param(CS_F_(eps));

  if (eqp->iwarni >= 1) {
    bft_printf(" Solving the variable%s \n ",
               cs_field_get_label(CS_F_(eps)));
  }

  int kstprv = cs_field_key_id("source_term_prev_id");
  int st_prv_id = cs_field_get_key_int(CS_F_(eps), kstprv);
  cs_real_t *c_st_prv = NULL, *cromo = NULL;
  if (st_prv_id > -1)
    c_st_prv = cs_field_by_id(st_prv_id)->val;

  int key_t_ext_id = cs_field_key_id("time_extrapolated");
  int iroext = cs_field_get_key_int(CS_F_(rho), key_t_ext_id);
  if ((iroext > 0) && (st_prv_id > -1))
    cromo = CS_F_(rho)->val_pre;
  else
    cromo = CS_F_(rho)->val;

  /* S as Source, V as Variable */
  const cs_time_scheme_t *time_scheme = cs_glob_time_scheme;
  const cs_real_t thets = time_scheme->thetst;
  const cs_real_t thetv = eqp->thetav;

  cs_array_set_value_real(n_cells, 1, 0, rhs);
  cs_array_set_value_real(n_cells, 1, 0, rovsdt);

  /* Work arrays */
  cs_real_t *w1;
  BFT_MALLOC(w1, n_cells_ext, cs_real_t);

  /* User source terms
   * ================= */

  cs_user_source_terms(cs_glob_domain,
                       CS_F_(eps)->id,
                       rhs,
                       rovsdt);

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
   * ======================= */

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
   * ================ */

  if (ncesmp > 0) {

    const int var_key_id = cs_field_key_id("variable_id");
    const int ivar_eps = cs_field_get_key_int(CS_F_(eps), var_key_id)-1;
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
   * ============== */

  if (eqp->istat == 1) {
#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      rovsdt[c_id] += (crom[c_id] / dt[c_id]) * cell_f_vol[c_id];
    }
  }

  /* Production (rho * Ce1 * epsilon / k * P)
   *    Dissipation (rho*Ce2.epsilon/k*epsilon)
   * ========================================== */

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
   * ============= */

  /* FIXME use beta ... WARNING */
  if (cs_glob_turb_rans_model->igrari == 1) {

    /* Extrapolation of source terms (2nd order in time) */
    if (st_prv_id > -1)
      _gravity_st_epsilon(gradro, cell_f_vol, c_st_prv);

    else
      _gravity_st_epsilon(gradro, cell_f_vol, rhs);

  }

  /* Diffusion term (Daly Harlow: generalized gradient hypothesis method)
   *===================================================================== */

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
  /* Scalar diffusivity */
  }
  else {

    if (eqp->idifft == 1) {
#     pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        w1[c_id] = viscl[c_id] + visct[c_id]/sigmae;
    }
    else
      cs_array_copy_real(n_cells, 1, viscl, w1);

    cs_face_viscosity(m,
                      fvq,
                      eqp->imvisf,
                      w1,
                      viscf,
                      viscb);
  }

  /* Solving
   * ======= */

  if (st_prv_id > -1) {
    const cs_real_t thetp1 = 1.+thets;
#   pragma omp parallel for if(n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      rhs[c_id] += thetp1*c_st_prv[c_id];
  }

  /* Translate coefa into cofaf and coefb into cofbf */

  cs_real_t *coefap = CS_F_(eps)->bc_coeffs->a;
  cs_real_t *coefbp = CS_F_(eps)->bc_coeffs->b;
  cs_real_t *cofafp = CS_F_(eps)->bc_coeffs->af;
  cs_real_t *cofbfp = CS_F_(eps)->bc_coeffs->bf;

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
                                     CS_F_(eps)->id,
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

/*----------------------------------------------------------------------------*/
/*! \brief Solve the equation on alpha in the framework of the Rij-EBRSM model.
 *
 * Also called for alpha of scalars for EB-DFM.
 *
 * \param[in]  f_id          field id of alpha variable
 * \param[in]  c_durbin_l    constant for the Durbin length
 !*/
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_solve_alpha(int        f_id,
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

  const cs_real_t *crom = CS_F_(rho)->val;
  const cs_real_t *viscl = CS_F_(mu)->val;
  const cs_real_t *cvara_ep = CS_F_(eps)->val_pre;
  const cs_real_t *cvara_al = cs_field_by_id(f_id)->val_pre;
  const cs_real_6_t *cvara_rij = (const cs_real_6_t *)CS_F_(rij)->val_pre;

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  int iflmas =  cs_field_get_key_int(CS_F_(vel), kimasf);
  int iflmab =  cs_field_get_key_int(CS_F_(vel), kbmasf);

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

  cs_array_set_value_real(n_cells, 1, 0, rhs);
  cs_array_set_value_real(n_cells, 1, 0, rovsdt);

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

  cs_array_set_value_real(n_cells, 1, 1, w1);

  cs_face_viscosity(m,
                    fvq,
                    eqp->imvisf,
                    w1,
                    viscf,
                    viscb);

  BFT_FREE(w1);

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

  /* Clipping
     ======== */

  cs_real_t *alpha_min;
  BFT_MALLOC(alpha_min, n_cells_ext, cs_real_t);

  /* Compute a first estimator of the minimal value of alpha per cell.
   * This is deduced from "alpha/L^2 - div(grad alpha) = 1/L^2" and assuming that
   * boundary cell values are 0. This value is thefore non zero but
   * much smaller than the wanted value. */

  cs_array_set_value_real(n_cells_ext, 1, 0, alpha_min);
  cs_array_copy_real(n_cells, 1, rovsdt, alpha_min);

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
 * \brief Gravity terms for terms
 *        For \f$R_{ij}\f$
 *
 * \param[in]   gradro    work array for \f$ \grad{\rho} \f$
 * \param[out]  buoyancy  buoyancy term
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_grav_st(const cs_real_t  gradro[][3],
                          cs_real_t        buoyancy[][6])
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  const cs_real_t *cvara_ep = (const cs_real_t *)CS_F_(eps)->val_pre;
  const cs_real_6_t *cvara_rij = (const cs_real_6_t *)CS_F_(rij)->val_pre;

  cs_real_t cons = -1.5*cs_turb_cmu;
  const cs_real_t uns3 = 1./3;

  const cs_field_t *tf = cs_thermal_model_field();
  if (tf != NULL) {
    const int ksigmas = cs_field_key_id("turbulent_schmidt");
    const cs_real_t turb_schmidt = cs_field_get_key_double(tf, ksigmas);
    cons = -1.5*cs_turb_cmu/turb_schmidt;
  }

  const cs_real_t *grav = cs_glob_physical_constants->gravity;
  const cs_real_t o_m_crij3 = (1. - cs_turb_crij3);

# pragma omp parallel for if(n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    cs_real_t rit[3];
    cs_math_sym_33_3_product(cvara_rij[c_id], gradro[c_id], rit);

     const cs_real_t kseps =   cs_math_6_trace(cvara_rij[c_id])
                             / (2*cvara_ep[c_id]);

     cs_real_t gij[3][3];
     for (cs_lnum_t i = 0; i < 3; i++) {
       for (cs_lnum_t j = 0; j < 3; j++)
         gij[i][j] = cons*kseps* (rit[i]*grav[j] + rit[j]*grav[i]);
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
 * \brief Terms of wall echo for \f$ R_{ij} \f$
 *        \f$var  = R_{11} \: R_{22} \: R_{33} \: R_{12} \: R_{13} \: R_{23}\f$
 *        \f$comp =  1 \:  2 \:  3 \:  4 \:  5 \:  6\f$
 *
 * \param[in]     produc  production
 * \param[in,out] rhs     work array for right-hand-side
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_echo(const cs_real_t  produc[][6],
                       cs_real_t        rhs[][6])
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  const cs_real_t *cell_f_vol = cs_glob_mesh_quantities->cell_f_vol;

  const cs_real_t *cvara_ep = (const cs_real_t *)CS_F_(eps)->val_pre;
  const cs_real_6_t *cvara_rij = (const cs_real_6_t *)CS_F_(rij)->val_pre;

  cs_real_t *cromo = NULL;
  int key_t_ext_id = cs_field_key_id("time_extrapolated");
  int iroext = cs_field_get_key_int(CS_F_(rho), key_t_ext_id);
  if ((cs_glob_time_scheme->isto2t > 0) && (iroext > 0))
    cromo = CS_F_(rho)->val_pre;
  else
    cromo = CS_F_(rho)->val;

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
    for (int i = 0; i < 3; i++)
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
   *============================== */

  cs_array_set_value_real(n_cells, 6, 0, (cs_real_t *)w6);

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
 * \brief Clip the turbulent Reynods stress tensor and the turbulent
 *        dissipation (coupled components version).
 *
 * \param[in]  n_cells  number of cells
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_clip(cs_lnum_t  n_cells)
{
  cs_real_t *cvar_ep = (cs_real_t *)CS_F_(eps)->val;
  cs_real_6_t *cvar_rij = (cs_real_6_t *)CS_F_(rij)->val;
  const cs_real_t *cvara_ep = (const cs_real_t *)CS_F_(eps)->val_pre;

  int kclipp = cs_field_key_id("clipping_id");

  /* Post-process clippings ? */

  cs_real_t *cpro_eps_clipped = NULL;
  cs_real_6_t *cpro_rij_clipped = NULL;
  int clip_e_id = cs_field_get_key_int(CS_F_(eps), kclipp);
  if (clip_e_id > -1) {
    cpro_eps_clipped = cs_field_by_id(clip_e_id)->val;
    cs_array_set_value_real(n_cells, 1, 0, cpro_eps_clipped);
  }
  int clip_r_id = cs_field_get_key_int(CS_F_(rij), kclipp);
  if (clip_r_id > -1) {
    cpro_rij_clipped = (cs_real_6_t *)cs_field_by_id(clip_r_id)->val;
    cs_array_set_value_real(n_cells,
                            CS_F_(rij)->dim, 0,
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
  const cs_real_t rijref = cs_math_fmax(trref/3, cs_math_epzero);

# pragma omp parallel if(n_cells > CS_THR_MIN)
  {
    cs_lnum_t t_icltot = 0;
    cs_lnum_t t_iclrij[6] = {0, 0, 0, 0, 0, 0};
    cs_lnum_t t_iclep[1] = {0};

    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n_cells, sizeof(cs_real_t), &s_id, &e_id);

    for (cs_lnum_t c_id = s_id; c_id < e_id; c_id++) {

      int is_clipped = 0;

      /* Check if R is positive and ill-conditionned (since the former
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
        for (int ii = 0; ii < 6; ii++)
          tensor[ii] = cvar_rij[c_id][ii]/trrij;

        cs_real_t eigen_vals[3];
        cs_math_sym_33_eigen(tensor, eigen_vals);

        cs_real_t eigen_min = eigen_vals[0];
        cs_real_t eigen_max = eigen_vals[0];
        for (int i = 1; i < 3; i++) {
          eigen_min = cs_math_fmin(eigen_min, eigen_vals[i]);
          eigen_max = cs_math_fmax(eigen_max, eigen_vals[i]);
        }

        /* If negative eigenvalue, return to isotropy */

        if (   (eigen_min <= eigen_tol*eigen_max)
            || (eigen_min < cs_math_epzero)) {

          is_clipped = 1;

          eigen_min = cs_math_fmin(eigen_min, -eigen_tol);
          cs_real_t eigen_offset
            = cs_math_fmin(-eigen_min/(1.0/3-eigen_min)+0.1, 1.0);

          for (cs_lnum_t ii = 0; ii < 6; ii++) {
            cvar_rij[c_id][ii] = (1.0-eigen_offset)*cvar_rij[c_id][ii];

            if (ii < 3)
              cvar_rij[c_id][ii] += trrij*(eigen_offset+eigen_tol)/3;

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

      /* Enforce Cauchy Schwarz inequality (only for x, y, z directions) */

      cs_real_t cvar_var1, cvar_var2;
      for (cs_lnum_t ii = 3; ii < 6; ii++) {
        if (ii == 3) {
          cvar_var1 = cvar_rij[c_id][0];
          cvar_var2 = cvar_rij[c_id][1];
        }
        else if (ii == 4) {
          cvar_var1 = cvar_rij[c_id][0];
          cvar_var2 = cvar_rij[c_id][2];
        }
        else if (ii == 5) {
          cvar_var1 = cvar_rij[c_id][1];
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

    for (int i = 0; i < 6; i++) {
      #pragma omp atomic
      iclrij[i] += t_iclrij[i];
    }

    #pragma omp atomic
    iclep[0] += t_iclep[0];

    #pragma omp atomic
    icltot += t_icltot;
  }

  /* Store number of clippings for logging */

  cs_lnum_t iclrij_max[6] = {0, 0, 0, 0, 0, 0}, iclep_max[1] = {0};

  cs_log_iteration_clipping_field(CS_F_(rij)->id, icltot, 0,
                                  vmin, vmax, iclrij, iclrij_max);

  cs_log_iteration_clipping_field(CS_F_(eps)->id, iclep[0], 0,
                                  vmin+6, vmax+6, iclep, iclep_max);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Clip the turbulent Reynods stress tensor and the turbulent
 *        dissipation (segregated version)
 *
 * \param[in]  n_cells  number of cells
 * \param[in]  iclip    if 0, viscl0 is used; otherwise viscl is used.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_clip_sg(cs_lnum_t  n_cells,
                          int        iclip)
{
  cs_real_t *cvar_ep = (cs_real_t *)CS_F_(eps)->val;
  cs_real_6_t *cvar_rij = (cs_real_6_t *)CS_F_(rij)->val;

  const cs_real_6_t *cvara_rij = (const cs_real_6_t *)CS_F_(rij)->val_pre;
  const cs_real_t *cvara_ep = (const cs_real_t *)CS_F_(eps)->val_pre;

  /* Store Min and Max for logging. */

  cs_real_t vmin[7], vmax[7];
  _rij_min_max(n_cells, cvar_rij, cvar_ep, vmin, vmax);

  /* Clipping (modified to avoid exactly zero values). */

  const cs_real_t epz2 = cs_math_pow2(cs_math_epzero);

  cs_lnum_t iclrij[6] = {0, 0, 0, 0, 0, 0};
  cs_lnum_t iclep[1] = {0};

  if (iclip == 1) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      for (int ii = 0; ii < 3; ii++) {
        if (cvar_rij[c_id][ii] <= epz2) {
          iclrij[ii]++;
          cvar_rij[c_id][ii] = epz2;
        }
      }
    }

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      if (cs_math_fabs(cvar_ep[c_id]) <= epz2) {
        iclep[0]++;
        cvar_ep[c_id] = cs_math_fmax(cvar_ep[c_id], epz2);
      }
      else if (cvar_ep[c_id] <= 0.0) {
        iclep[0]++;
        cvar_ep[c_id] = cs_math_fabs(cvar_ep[c_id]);
      }
    }
  }
  else {
    const cs_real_t varrel = 1.1;
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      for (cs_lnum_t ii = 0; ii < 3; ii++) {
        if (cs_math_fabs(cvar_rij[c_id][ii]) <= epz2) {
          iclrij[ii]++;
          cvar_rij[c_id][ii] = cs_math_fmax(cvar_rij[c_id][ii], epz2);
        }
        else if (cvar_rij[c_id][ii] <= 0) {
          iclrij[ii]++;
          cvar_rij[c_id][ii] = cs_math_fmin(cs_math_fabs(cvar_rij[c_id][ii]),
                                            varrel*fabs(cvara_rij[c_id][ii]));
        }
      }
    }

    iclep[0] = 0;
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      if (cs_math_fabs(cvar_ep[c_id]) < epz2) {
        iclep[0]++;
        cvar_ep[c_id] = cs_math_fmax(cvar_ep[c_id], epz2);
      }
      else if(cvar_ep[c_id] <= 0) {
        iclep[0]++;
        cvar_ep[c_id] = cs_math_fmin(cs_math_fabs(cvar_ep[c_id]),
                                     varrel*cs_math_fabs(cvara_ep[c_id]));
      }
    }
  }

  /* Enforce Cauchy Schwarz inequality (only for x, y, z directions) */
  for (cs_lnum_t ii = 3; ii < 6; ii++) {
    cs_lnum_t ii1, ii2;
    if (ii == 3) {
      ii1 = 0;
      ii2 = 1;
    }
    else if (ii == 4) {
      ii1 = 1;
      ii2 = 2;
    }
    else if (ii == 5) {
      ii1 = 0;
      ii2 = 2;
    }
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      const cs_real_t rijmin = sqrt(cvar_rij[c_id][ii1]*cvar_rij[c_id][ii2]);
      if (rijmin < cs_math_fabs(cvar_rij[c_id][ii])) {
        cvar_rij[c_id][ii] = _sign(1., cvar_rij[c_id][ii])*rijmin;
        iclrij[ii]++;
      }
    }
  }

  /* Store number of clippings for logging */

  cs_lnum_t icltot = 0;
  for (int ii = 0; ii < 6; ii++) {
    icltot += iclrij[ii];
  }

  cs_lnum_t iclrij_max[6] = {0, 0, 0, 0, 0, 0}, iclep_max[1] = {0};

  cs_log_iteration_clipping_field(CS_F_(rij)->id, icltot, 0,
                                  vmin, vmax, iclrij, iclrij_max);

  cs_log_iteration_clipping_field(CS_F_(eps)->id, iclep[0], 0,
                                  vmin+6, vmax+6, iclep, iclep_max);
}

/*----------------------------------------------------------------------------*/

 END_C_DECLS
